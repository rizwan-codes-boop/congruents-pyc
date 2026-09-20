"""Self-contained preparation, ownership and native-boundary checks."""
import gc
import hashlib
import json
import importlib.util
import math
import tempfile
import unittest
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from unittest.mock import patch

import numpy as np

from congruents import Context, Galaxy, Grid, Preparation, load_catalogue, write_catalogue
from congruents.physics.ionisation import ionisation
from congruents.galaxy import properties as galaxy_properties
from congruents.preparation import Table, TableData, cache_filename, rename_cache_files
from congruents.quadrature import integrate

ROOT = Path(__file__).resolve().parents[1]
GALAXIES = load_catalogue(ROOT/"input/cat_nt.txt")

class SerialTests(unittest.TestCase):
    """Check analytic kernels and serial determinism without native model physics."""
    def test_process_module_ownership_and_dispatch(self):
        from congruents import tables, grids
        from congruents.physics import radiation, ionisation, inverse_compton, bremsstrahlung, synchrotron
        self.assertEqual(grids.log_grid.__module__, "congruents.grids")
        self.assertEqual(radiation.photon.__module__, "congruents.physics.radiation")
        self.assertEqual(ionisation.ionisation.__module__, "congruents.physics.ionisation")
        self.assertFalse((ROOT/"src/congruents/serial.py").exists())
        config = (1e-16, 1e8, .002, 1e8, 1e-16, 1e-7, 1e-12, 1e6)
        for kind, module in (("emission", inverse_compton), ("gamma", inverse_compton),
                             ("bs", bremsstrahlung), ("sy", synchrotron)):
            with self.subTest(kind=kind), patch.object(module, "generate", return_value=object()) as call:
                self.assertIs(tables.generate(kind, "3000", 0., 4, 4, config), call.return_value)
                call.assert_called_once()

    def test_quadrature_analytic(self):
        for f, low, high, expected in (
            (lambda x: x**4, 0., 1., .2),
            (math.sin, 0., math.pi, 2.),
            (lambda x: math.exp(-100*x), 0., 1., .01),
            (lambda x: 0., -1., 1., 0.),
        ):
            self.assertAlmostEqual(integrate(f, low, high), expected, places=12)
        with self.assertRaises(ValueError):
            integrate(math.sin, 1., 0.)
        with self.assertRaises(RuntimeError):
            integrate(lambda x: float("nan"), 0., 1.)

    def test_table_interpolation(self):
        with Table(TableData([1.,2.], [1.,2.], [2.,3.,3.,4.]), {}) as t:
            self.assertEqual(t.evaluate([1.,1.5,2.,3.], [1.,1.5,2.,1.]), [2.,3.,4.,0.])
            with self.assertRaises(ValueError):
                t.evaluate([1.], [1.,2.])
        with Table(TableData([1.,2.], [1.], [2.,4.]), {}) as t:
            self.assertEqual(t.evaluate([0.,1.5,3.]), [0.,3.,0.])
        for data in (TableData([2.,1.],[1.],[1.,2.]),
                     TableData([1.,2.],[1.],[-1.,2.]),
                     TableData([1.,2.],[1.],[float("nan"),2.])):
            with self.assertRaises(ValueError):
                Table(data, {})

    def test_ionisation_reference_values(self):
        expected = [-2.0513467838531972e-16, -3.7936953668288616e-16,
                    -4.955261088812638e-16]
        np.testing.assert_allclose(ionisation([.001, 1., 100.], 1.), expected,
                                   rtol=1e-14, atol=0.)

    def test_preparation_thread_setting_independence(self):
        galaxies = load_catalogue(ROOT/"input/cat_nt.txt")
        with Preparation(galaxies, Grid(8,8), threads=1) as a, \
             Preparation(galaxies, Grid(8,8), threads=4) as b:
            self.assertEqual(a.properties, b.properties)
            self.assertEqual(a.config, b.config)
            self.assertEqual(a.temperature_bounds, b.temperature_bounds)
            for i, p in enumerate(a.properties):
                self.assertTrue(all(math.isfinite(v) and v > 0 for v in p.values()))
                factor = 3. if math.log10(galaxies[i].sfr_msun_per_year /
                                           galaxies[i].stellar_mass_msun) > -10 else 1.5
                self.assertEqual(p["halo_magnetic_field_gauss"], p["magnetic_field_gauss"]/factor)


class ContextTests(unittest.TestCase):

    """Exercise configuration lifetime, unit conventions and failure recovery."""
    def test_lifecycle(self):
        for _ in range(1000):
            with Context() as context:
                self.assertEqual(context.ionisation_loss([], 1.), [])
            context.close()
            with self.assertRaises(RuntimeError):
                context.ionisation_loss([1.], 1.)

    def test_validation_survives(self):
        with Context() as context:
            for energies, density in [([0.], 1.), ([.0001], 1.),
                ([float("nan")], 1.), ([float("inf")], 1.),
                ([1.], -1.), ([1.], float("nan"))]:
                with self.assertRaises(ValueError):
                    context.ionisation_loss(energies, density)
            self.assertLess(context.ionisation_loss([1.], 1.)[0], 0)
        for threads in (0, -1, 2**40):
            with self.assertRaises(ValueError):
                Context(threads)
        with self.assertRaises(TypeError):
            Context(1.5)

    def test_missing_library(self):
        with Context(library=ROOT/"does-not-exist.so") as context:
            self.assertLess(context.ionisation_loss([1.], 1.)[0], 0)
            with self.assertRaises(FileNotFoundError):
                _ = context.openmp_enabled

    def test_numerical_error_recovery(self):
        with Context() as context:
            with self.assertRaises(RuntimeError):
                context.ionisation_loss([1.e308], 1.)
            self.assertLess(context.ionisation_loss([1.], 1.)[0], 0)

    def test_concurrent_contexts(self):
        def evaluate(threads):
            """Evaluate one independent context to exercise concurrent-call isolation."""
            with Context(threads) as context:
                return context.ionisation_loss([.01, 1., 100.] * 100, 10.)
        expected = evaluate(1)
        with ThreadPoolExecutor(max_workers=4) as pool:
            for actual in pool.map(evaluate, [1, 2, 3, 4] * 4):
                self.assertEqual(actual, expected)


    def test_version_and_density_scaling(self):
        for threads in (1, 4):
            with Context(threads) as context:
                self.assertEqual(context.version, "0.4.0")
                energy = [.001, .01, 1., 100., 1e5]
                unit = np.array(context.ionisation_loss(energy, 1.))
                for density in (0., 1e-3, 1., 1e3):
                    np.testing.assert_allclose(context.ionisation_loss(energy, density),
                                               unit*density, rtol=1e-14, atol=0.)


class PreparationTests(unittest.TestCase):
    """Check catalogue/table/cache contracts using deliberately small grids."""
    def assertClose(self, actual, expected, tolerance=1e-12):
        """Compare floating-point sequences with the test suite numerical tolerance."""
        self.assertEqual(len(actual), len(expected))
        for i, (a, b) in enumerate(zip(actual, expected)):
            self.assertTrue(math.isclose(a, b, rel_tol=tolerance, abs_tol=1e-280),
                            (i, a, b))


    def test_catalogue_roundtrip_and_validation(self):
        self.assertEqual(len(GALAXIES), 11)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/"cat.txt"
            write_catalogue(path, GALAXIES)
            self.assertEqual(load_catalogue(path), GALAXIES)
            path.write_text("n_gal\n2\nz Mstar__Msol Re__kpc SFR__Msolyrm1\n0 1 1 1\n")
            with self.assertRaises(ValueError):
                load_catalogue(path)
        for row in ((0, 1, 0, 1), (0, 1, 1, 0), (-1, 1, 1, 1), (0, float("nan"), 1, 1)):
            with self.assertRaises(ValueError):
                Galaxy(*row)
        for n in (0, 1, 4097):
            with self.assertRaises(ValueError):
                Grid(n, 8)
        with self.assertRaises(TypeError):
            Grid(2.5, 8)

    @unittest.skipUnless(importlib.util.find_spec("astropy"), "optional Astropy extra")
    def test_astropy_units(self):
        from astropy import units as u
        g = Galaxy.from_quantities(.01, 1e10*u.Msun, 2000*u.pc, 2*u.Msun/u.yr)
        self.assertEqual(g.radius_kpc, 2.)
        with self.assertRaises(u.UnitConversionError):
            Galaxy.from_quantities(.01, 1e10*u.Msun, 2*u.s, 2*u.Msun/u.yr)


    def test_cache_cold_warm_and_corruption(self):
        with tempfile.TemporaryDirectory() as tmp:
            with Preparation(GALAXIES, Grid(8, 8), cache=tmp) as p:
                p.prepare_tables()
                self.assertEqual((p.cache_hits, p.cache_misses), (0, 26))
                reference = {k: t.snapshot() for k, t in p._tables.items()}
            with Preparation(GALAXIES, Grid(8, 8), cache=tmp) as p:
                p.prepare_tables()
                self.assertEqual((p.cache_hits, p.cache_misses), (26, 0))
                self.assertEqual({k: t.snapshot() for k, t in p._tables.items()}, reference)
            # Deliberately corrupt one payload byte: checksum rejection must occur
            # before any damaged lookup values can influence a spectrum.
            path = next(Path(tmp).glob("*.cgt"))
            raw = bytearray(path.read_bytes())
            raw[-1] ^= 1
            path.write_bytes(raw)
            with Preparation(GALAXIES, Grid(8, 8), cache=tmp) as p:
                with self.assertRaises(ValueError):
                    p.prepare_tables()
            # Grid changes produce distinct keys, not reuse of wrong-shaped tables.
            with Preparation(GALAXIES, Grid(9, 8), cache=tmp) as p:
                p.table("bs")
                self.assertEqual(p.cache_misses, 1)

    def test_readable_cache_migration_and_progress(self):
        messages = []
        with tempfile.TemporaryDirectory() as tmp:
            with Preparation(GALAXIES, Grid(8,8), cache=tmp, progress=messages.append) as p:
                expected = p.table("bs").snapshot()
                meta = p.table("bs").metadata
            named = next(Path(tmp).glob("*.cgt"))
            self.assertEqual(named.name,"bremsstrahlung_8x8.cgt")
            original = named.read_bytes()
            key = hashlib.sha256(json.dumps(meta,sort_keys=True,allow_nan=False).encode()).hexdigest()
            old = named.with_name(key+".cgt")
            named.rename(old)
            messages.clear()
            with patch("congruents.tables.generate", side_effect=AssertionError("Regenerated cache")):
                with Preparation(GALAXIES, Grid(8,8), cache=tmp, progress=messages.append) as p:
                    self.assertEqual(p.table("bs").snapshot(), expected)
                    self.assertEqual((p.cache_hits,p.cache_misses),(1,0))
            self.assertFalse(old.exists())
            self.assertEqual(named.read_bytes(),original)
            self.assertTrue(any("loading cached bremsstrahlung" in m for m in messages))
            self.assertFalse(any("generating" in m for m in messages))
            self.assertTrue(any("Table 1 ready" in m for m in messages))
            named.rename(old)
            self.assertEqual(rename_cache_files(tmp),[(old,named)])
            self.assertEqual(named.read_bytes(),original)
            self.assertEqual(rename_cache_files(tmp),[])
            old.write_bytes(original)
            with self.assertRaises(FileExistsError):
                rename_cache_files(tmp)
            self.assertEqual(old.read_bytes(),original)
            self.assertEqual(named.read_bytes(),original)

    def test_cache_labels_and_generation_messages(self):
        messages=[]
        with Preparation(GALAXIES,Grid(4,4),progress=messages.append) as p:
            for kind,field,temp,prefix in (
                    ("emission","CMB",2.7,"ic_emission_CMB_T2.7K_4x4.cgt"),
                    ("gamma","FIR",25.,"ic_energy_transfer_FIR_T25K_4x4.cgt"),
                    ("sy","3000",0.,"synchrotron_4.cgt")):
                table=p.table(kind,field,temp)
                name=cache_filename(table.metadata)
                self.assertEqual(name,prefix)
                changed=dict(table.metadata,implementation_sha256="different")
                self.assertEqual(cache_filename(changed),name)
        self.assertEqual(sum("generating" in m for m in messages),3)
        self.assertEqual(sum("ready:" in m for m in messages),3)
        with self.assertRaises(TypeError):
            Preparation(GALAXIES,progress=True)

    def test_named_cache_rejects_changed_fingerprint(self):
        messages=[]
        with tempfile.TemporaryDirectory() as tmp:
            with Preparation(GALAXIES,Grid(4,4),cache=tmp) as p:
                p.table("bs")
            with Preparation(GALAXIES,Grid(4,4),cache=tmp,progress=messages.append) as p:
                p._preparation_hash="changed-physics"
                p.table("bs")
                self.assertEqual((p.cache_hits,p.cache_misses),(0,1))
            self.assertTrue(any("incompatible metadata" in message for message in messages))
            self.assertEqual(len(list(Path(tmp).glob("*.cgt"))),1)
            with Preparation(GALAXIES,Grid(4,4),cache=tmp) as p:
                p._preparation_hash="changed-physics"
                p.table("bs")
                self.assertEqual((p.cache_hits,p.cache_misses),(1,0))

    def test_lifetime_and_invalid_inputs(self):
        p = Preparation(GALAXIES, Grid(8, 8))
        with self.assertRaises(ValueError):
            p.radiation(0, [float("nan")])
        with self.assertRaises(ValueError):
            p.table("emission", "CMB", -1)
        t = p.table("bs")
        p.close()
        p.close()
        with self.assertRaises(RuntimeError):
            t.snapshot()
        with self.assertRaises(RuntimeError):
            p.radiation(0, [1e-9])
        # Legacy single-plane temperatures are rejected instead of underflowing.
        with Preparation([GALAXIES[0]], Grid(8, 8)) as single:
            with self.assertRaises(ValueError):
                single.temperature_grid("FIR")

    def test_numpy_ownership_and_views(self):
        with Preparation(GALAXIES, Grid(8, 8)) as p:
            table = p.table("bs")
            self.assertEqual(table.values.shape, (8, 8))
            self.assertEqual(table.values.dtype, np.dtype("float64"))
            self.assertTrue(table.values.flags.c_contiguous)
            for buffer in (table._x, table._y, table._values):
                self.assertTrue(buffer.flags.owndata)
                self.assertFalse(buffer.flags.writeable)
            with self.assertRaises(ValueError):
                table.values[0, 0] = 1.
            with self.assertRaises(ValueError):
                table.values.flags.writeable = True
            retained = table.values
            expected = retained.copy()
        gc.collect()
        # Closing the table handle must not invalidate retained NumPy storage.
        np.testing.assert_array_equal(retained, expected)

    def test_cache_import_numpy_lifetime(self):
        with Preparation(GALAXIES, Grid(2, 2)) as p:
            x = np.array([1., 2.])
            y = np.array([1., 2.])
            z = np.array([2., 3., 3., 4.])
            table = Table(TableData(x, y, z), {})
            x[:] = 99.
            del x, y, z
            gc.collect()
            with table:
                self.assertEqual(table.evaluate([1.5], [1.5]), [3.])
                self.assertTrue(table._values.flags.owndata)


    def test_radiation_all_galaxies(self):
        with Preparation(GALAXIES, Grid(8, 8)) as p:
            for i in range(len(GALAXIES)):
                result = p.radiation(i, [1e-12, 1e-9])
                self.assertTrue(all(v > 0 and math.isfinite(v) for v in result["urad_ub_ev_cm3"]))
                self.assertTrue(all(v >= 0 and math.isfinite(v)
                    for field in result["fields_cm3_gev"].values() for v in field))

    def test_table_families_and_combined_ic(self):
        with Preparation(GALAXIES, Grid(16, 16)) as p:
            self.assertEqual(p.prepare_tables(), 26)
            for table in p._tables.values():
                data = table.snapshot()
                self.assertEqual(len(data.values), len(data.x)*len(data.y))
                self.assertTrue(np.isfinite(data.values).all())
                self.assertTrue(np.all(np.array(data.values) >= 0))
                self.assertTrue(np.all(np.diff(data.x) > 0))
                self.assertTrue(np.all(np.diff(data.y) > 0))
            emission, gamma = (p.table(k).snapshot() for k in ("emission", "gamma"))
            self.assertClose([emission.x[0], gamma.x[0]], [p.config[0], p.config[4]])
            self.assertNotEqual(emission.values, gamma.values)
    def test_table_generation_never_calls_native_code(self):
        # Property preparation now correctly belongs to the native galaxy loop.
        # After that stage, serial table generation must not load any C code.
        with Preparation(GALAXIES,Grid(8,8)) as p:
            with patch("ctypes.CDLL",side_effect=AssertionError("Tables loaded C")):
                self.assertEqual(p.prepare_tables(),26)
                self.assertTrue(p.radiation(0,[1e-9])["fields_cm3_gev"])

    def test_property_validation_ownership_and_halo_threshold(self):
        rows = np.array([[0., 1e10, 2., 1.], [0., 1e10, 2., 1.001]])
        before = rows.copy()
        values = galaxy_properties(rows)
        np.testing.assert_array_equal(rows, before)
        self.assertEqual(values.shape, (2, 10))
        self.assertTrue(values.flags.c_contiguous)
        self.assertTrue(values.flags.owndata)
        self.assertFalse(values.flags.writeable)
        self.assertEqual(values[0, 9], values[0, 2]/1.5)
        self.assertEqual(values[1, 9], values[1, 2]/3.)
        for bad in ([], [0., 1., 1., 1.], [[21., 1., 1., 1.]],
                    [[0., 0., 1., 1.]], [[0., 1., np.nan, 1.]]):
            with self.assertRaises(ValueError):
                galaxy_properties(bad)
        with self.assertRaisesRegex(RuntimeError, r"indices \[0\]"):
            galaxy_properties([[0., 1e10, 1e300, 1.]])
        self.assertTrue(np.isfinite(galaxy_properties(rows)).all())
