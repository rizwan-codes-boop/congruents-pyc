"""Self-contained observer, exports and optional-comparison checks."""
from pathlib import Path
import subprocess
import json
import unittest
import tempfile
import os
import sys
import numpy as np
from congruents import load_catalogue
from congruents.observer import luminosity_distance_mpc, distance_factor_cm2, observer_sed
from congruents import Preparation, Grid, SolverGrid
from congruents.pipeline import run, COMPONENT_FILES
from congruents.attenuation import EBLTable
from compare_outputs import compare, FILES

ROOT = Path(__file__).resolve().parents[1]


class ObserverTests(unittest.TestCase):


    """Test observer transforms, complete export layout and thread equivalence."""
    def test_cli_refuses_existing_directory_before_work(self):
        with tempfile.TemporaryDirectory() as temp:
            proc = subprocess.run([sys.executable,"-m","congruents",str(ROOT/"input/cat_nt.txt"),temp],
                cwd=ROOT,env={**os.environ,"PYTHONPATH":str(ROOT/"src")},capture_output=True,text=True)
            self.assertNotEqual(proc.returncode,0)
            self.assertIn("Output already exists",proc.stderr)

    def test_ebl_nodes_and_invalid_file(self):
        table = EBLTable(ROOT/"input/tau_Eg_z_Franceschini.txt")
        for i,z in enumerate(table.z):
            got = table.optical_depth(table.energy/1e9*(1+z),[z])[0]
            np.testing.assert_allclose(got,table.values[i],rtol=1e-12,atol=1e-12)
        self.assertTrue(np.all(table.optical_depth([1e-16,1e8],[0.,.001,4.])>=0))
        with tempfile.NamedTemporaryFile(mode="w") as f:
            f.write("2 2 0 1");f.flush()
            with self.assertRaises(ValueError):
                EBLTable(f.name)


    def test_required_attenuation_and_validation(self):
        for z in ([0.],[-1.],[np.nan],[],.1):
            with self.assertRaises(ValueError):
                luminosity_distance_mpc(z)
        with self.assertRaises(TypeError):
            observer_sed([1,2],[[1,1]],[.1])
        for e, source, tau in (([2,1],[[1,1]],[[0,0]]),
                              ([1,2],[[1,-1]],[[0,0]]),
                              ([1,2],[[1,1]],[[0,-1]]),
                              ([1,2],[[1,1]],[0,0])):
            with self.assertRaises(ValueError):
                observer_sed(e,source,[.1],tau_internal=tau,tau_ebl=[[0,0]])

    def test_explicit_zero_attenuation_diagnostic(self):
        got = observer_sed([1,2,4],[[1,2,4]],[1.],
                           tau_internal=[[0,0,0]],tau_ebl=[[0,0,0]])
        np.testing.assert_allclose(got[0],np.array([2,16,0])*distance_factor_cm2([1.])[0])

    def test_geometry_against_fixed_quadrature(self):
        from numpy.polynomial.legendre import leggauss
        from congruents.constants import PC
        z = np.r_[[g.redshift for g in load_catalogue(ROOT/"input/cat_nt.txt")], .1, 1., 3.]
        nodes, weights = leggauss(64)
        x = (nodes[None, :]+1)*z[:, None]/2
        integral = z/2*np.sum(weights/np.sqrt(.3*(1+x)**3+.7), axis=1)
        expected = (1+z)*299792.458/70*integral
        np.testing.assert_allclose(luminosity_distance_mpc(z), expected, rtol=1e-12)
        np.testing.assert_allclose(distance_factor_cm2(z),
                                  (1+z)**2/(4*np.pi*(expected*1e6*PC)**2), rtol=1e-12)

    def test_ebl_analytic_extrapolation_and_cap(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)/"ebl.txt"
            # Synthetic tau(z, E) = z + E/1e9, including extrapolation.
            path.write_text("2 2\n0 1\n1e9 2e9\n1 2 2 3\n")
            table = EBLTable(path)
            z, e = np.array([0., .1, 4.]), np.array([1e-16, 1., 12., 1e8])
            expected = z[:, None]+np.minimum(e[None, :]/(1+z[:, None]), 1e6)
            np.testing.assert_allclose(table.optical_depth(e, z), expected, atol=1e-12)
        for name in ("Dominguez", "Gilmore"):
            table = EBLTable(ROOT/"input"/f"tau_Eg_z_{name}.txt")
            for i, z in enumerate(table.z):
                np.testing.assert_allclose(table.optical_depth(table.energy/1e9*(1+z), [z])[0],
                                           table.values[i], rtol=1e-12, atol=1e-12)

    def test_end_to_end_exports_and_threads(self):
        # Two resolutions catch axis/shape assumptions; neither asserts that a
        # coarse model is physically converged or matches digitised observations.
        for cells, photons in ((16, 8), (32, 16)):
            with self.subTest(cells=cells), tempfile.TemporaryDirectory() as temp:
                destination = Path(temp)/"output"
                with Preparation(load_catalogue(ROOT/"input/cat_nt.txt"), Grid(cells, cells)) as p:
                    result = run(p, SolverGrid(cells, photons, cells), threads=4)
                    if cells == 16:
                        single = run(p, SolverGrid(cells, photons, cells), threads=1)
                        for name in result.components:
                            np.testing.assert_array_equal(result.components[name], single.components[name])
                        for name in result.source.diagnostics:
                            np.testing.assert_array_equal(result.source.diagnostics[name], single.source.diagnostics[name])
                result.export(destination)
                exported=json.loads((destination/"metadata.json").read_text())
                for key in ("python_source_sha256", "python_source_files_sha256", "python_source_hash_scheme"):
                    self.assertEqual(exported[key],result.source.metadata[key])
                self.assertIn("physics/inverse_compton.py",exported["python_source_files_sha256"])
                self.assertEqual(len(FILES), 39)
                self.assertTrue(all((destination/name).is_file() for name in FILES))
                for name, filename in COMPONENT_FILES.items():
                    values = result.components[name]
                    self.assertEqual(values.shape, (11, photons))
                    self.assertTrue(np.isfinite(values).all())
                    self.assertTrue(np.all(values >= 0))
                    np.testing.assert_allclose(np.loadtxt(destination/filename, skiprows=1), values,
                                               rtol=6e-7, atol=1e-280)
                for name, values in result.source.diagnostics.items():
                    filename = destination/("tau_loss/"+name+".txt" if name.startswith("tau_loss") else name+".txt")
                    self.assertFalse(np.isnan(values).any())
                    if not name.startswith("tau_loss"):
                        self.assertTrue(np.isfinite(values).all())
                    np.testing.assert_allclose(np.loadtxt(filename, skiprows=0 if name.startswith("CR_specs") else 1),
                                               values, rtol=6e-7, atol=1e-280)
                for filename, values in {
                    "T_CR.txt": result.source.kinetic_energy_gev,
                    "E_gam.txt": result.source.photon_energy_gev,
                    "fcal.txt": result.source.transport["fcal"],
                    "tau_ff.txt": result.source.emission["tau_free_free"],
                    "distmod.txt": result.distance_factor,
                    "L_gamma.txt": result.gamma_luminosity,
                    "tau_gg.txt": result.tau_internal,
                    "tau_EBL.txt": result.tau_ebl,
                    "spec_total_photons.txt": result.total_photon_sed,
                }.items():
                    np.testing.assert_allclose(np.loadtxt(destination/filename, skiprows=1),
                                               values, rtol=6e-7, atol=1e-280)
                photon_components = [values for name, values in result.components.items()
                                     if name not in ("neutrino", "pion_full_calorimetry")]
                np.testing.assert_array_equal(result.total_photon_sed, sum(photon_components))
                with self.assertRaises(FileExistsError):
                    result.export(destination)


class ComparisonTests(unittest.TestCase):
    """Use synthetic files to test the comparator itself, not astrophysical agreement."""
    def setUp(self):
        """Create isolated temporary output fixtures for comparison tests."""
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.reference, self.actual = Path(tmp.name)/"reference", Path(tmp.name)/"actual"
        for directory in (self.reference, self.actual):
            for name in FILES:
                path = directory/name
                path.parent.mkdir(parents=True, exist_ok=True)
                header = "" if name.startswith("CR_specs") else "header\n"
                path.write_text(header+"0 1 inf\n2 3 4\n")

    def test_equal_and_invalid_outputs(self):
        self.assertTrue(all(r["passed"] for r in compare(self.reference, self.actual).values()))
        path = self.actual/"T_CR.txt"
        for values in ("1e-290 1 inf\n2 3 4", "0 1 -inf\n2 3 4", "0 nan inf\n2 3 4",
                       "0 2 inf\n2 3 4", "0 1", ""):
            with self.subTest(values=values):
                path.write_text("header\n"+values)
                self.assertFalse(compare(self.reference, self.actual)["T_CR.txt"]["passed"])
        path.unlink()
        self.assertFalse(compare(self.reference, self.actual)["T_CR.txt"]["passed"])
        with self.assertRaises(ValueError):
            compare(self.reference, self.reference)

    def test_optional_cli_status(self):
        import json
        command = [sys.executable, str(ROOT/"tests/compare_outputs.py")]
        for args in ([], [str(self.reference/"absent"), str(self.actual)]):
            proc = subprocess.run(command+args, capture_output=True, text=True)
            self.assertEqual(proc.returncode, 2)
            self.assertEqual(json.loads(proc.stdout)["status"], "not_run")
        proc = subprocess.run(command+[str(self.reference), str(self.actual)], capture_output=True, text=True)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(json.loads(proc.stdout)["status"], "passed")
        (self.actual/"T_CR.txt").unlink()
        proc = subprocess.run(command+[str(self.reference), str(self.actual)], capture_output=True, text=True)
        self.assertEqual(proc.returncode, 1)
        self.assertEqual(json.loads(proc.stdout)["status"], "failed")

    def test_vector_orientation_without_flattening_spectra(self):
        for name in ("T_CR.txt", "E_gam.txt", "distmod.txt", "L_gamma.txt"):
            (self.reference/name).write_text("header\n1 2 3\n")
            (self.actual/name).write_text("header\n1\n2\n3\n")
            self.assertTrue(compare(self.reference, self.actual)[name]["passed"])
        name = "spec_pi.txt"
        (self.reference/name).write_text("header\n1 2 3\n")
        (self.actual/name).write_text("header\n1\n2\n3\n")
        self.assertFalse(compare(self.reference, self.actual)[name]["passed"])
