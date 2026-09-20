"""Self-contained solver checks; independent output parity is opt-in."""
import ctypes as ct
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
from congruents import Preparation, Grid, SolverGrid, load_catalogue, solve
from congruents.solver import _load
from congruents._native import IP, pointer, check
from congruents.preparation import PROPERTY_NAMES
from congruents.provenance import python_source_provenance

ROOT = Path(__file__).resolve().parents[1]
FILES = (
    ("transport", ("fcal","Dp","Dd","Dh","Q1","Q2","protons")),
    ("electrons", ("e1d","e2d","e1h","e2h")),
    ("emission", ("ic1d","ic2d","bs1d","bs2d","sy1d","sy2d","ic1h","ic2h",
                  "sy1h","sy2h","ff","tau_ff","pi","pi_fcal1","nu")),
)


class SolverTests(unittest.TestCase):
    """All-case source checks; grid/precision variants are regression tests, not convergence."""
    def test_recursive_source_provenance(self):
        """Verify recursive hashes track physics edits, paths and file additions without absolute paths."""
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            (root/"physics").mkdir()
            (root/"model.py").write_text("value = 1\n")
            physics=root/"physics"/"model.py"
            physics.write_text("value = 2\n")
            before=python_source_provenance(root)
            self.assertEqual(set(before["python_source_files_sha256"]),{"model.py","physics/model.py"})
            self.assertEqual(before,python_source_provenance(root))
            (root/"ignored.pyc").write_bytes(b"bytecode")
            self.assertEqual(before,python_source_provenance(root))
            physics.write_text("value = 3\n")
            changed=python_source_provenance(root)
            self.assertNotEqual(before["python_source_sha256"],changed["python_source_sha256"])
            physics.rename(root/"physics"/"renamed.py")
            renamed=python_source_provenance(root)
            self.assertNotEqual(changed["python_source_sha256"],renamed["python_source_sha256"])
            physics.write_text("value = 4\n")
            self.assertNotEqual(renamed["python_source_sha256"],python_source_provenance(root)["python_source_sha256"])
            physics.unlink()
            self.assertEqual(renamed,python_source_provenance(root))

    def test_output_provenance_includes_physics(self):
        """Ensure saved source metadata carries the same complete package manifest as the live result."""
        meta=self.result.metadata
        expected=python_source_provenance()
        for key,value in expected.items():
            self.assertEqual(meta[key],value)
        self.assertEqual(meta["solver_python_sha256"],meta["python_source_sha256"])
        self.assertIn("physics/inverse_compton.py",meta["python_source_files_sha256"])
        path=self.folder/"provenance.npz"
        self.result.save(path)
        with np.load(path,allow_pickle=False) as data:
            saved=json.loads(str(data["metadata_json"]))
            self.assertEqual(saved["python_source_files_sha256"],expected["python_source_files_sha256"])
    @classmethod
    def setUpClass(cls):
        """Prepare a shared small-grid result and temporary workspace for solver regressions."""
        cls.tmp = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.tmp.cleanup)
        cls.folder = Path(cls.tmp.name)
        cls.catalogue = load_catalogue(ROOT/"input/cat_nt.txt")
        cls.p = Preparation(cls.catalogue, Grid(16,16))
        cls.addClassCleanup(cls.p.close)
        # Reuse one all-case reduced-grid solution across tests to limit runtime.
        # Independent external-output parity is handled by compare_outputs.py.
        cls.grid = SolverGrid(16,8,16)
        cls.result = solve(cls.p, cls.grid, threads=4)


    def test_thread_count_equivalence(self):
        serial = solve(self.p, self.grid, threads=1)
        for group, _ in FILES:
            for name, values in getattr(self.result, group).items():
                np.testing.assert_array_equal(values, getattr(serial, group)[name])

    def test_native_transport_and_borrowed_arrays(self):
        rows = np.ascontiguousarray([g.as_row() for g in self.p.galaxies])
        props = np.ascontiguousarray([[p[k] for k in PROPERTY_NAMES] for p in self.p.properties])
        before = [v.copy() for v in (rows,props)]
        out,cp = np.zeros((11,7,16)),np.empty(11)
        status = np.zeros(11,dtype=np.int32)
        lib=_load()
        for v in (rows,props): v.flags.writeable=False
        check(lib.cg_transport(4,11,16,pointer(rows),pointer(props),pointer(self.result.kinetic_energy_gev),
              pointer(out),pointer(cp),status.ctypes.data_as(IP)),status)
        for i,(name,v) in enumerate(self.result.transport.items()):
            if i != 5: np.testing.assert_array_equal(out[:,i],v)
        for v,copy in zip((rows,props),before): np.testing.assert_array_equal(v,copy)

    def test_parallel_boundary_and_removed_implementations(self):
        source=(ROOT/"csrc/solver.c").read_text()
        for name in ("properties_one","transport_one","spectra_one","CRe_steadystate_solve", "tau_gg_gal_BW"):
            self.assertIn(name,source)
        self.assertIn("#pragma omp parallel for",source)
        self.assertNotIn("PyObject",source)
        for name in ("steady_state","emission","hadronic","cubature","solver_inputs","diagnostics"):
            self.assertFalse((ROOT/"src/congruents"/(name+".py")).exists())
        native="\n".join(p.read_text() for p in (ROOT/"csrc").rglob("*.h"))
        self.assertNotIn("init_do_2D_IC",native)
        self.assertFalse((ROOT/"CR_spectra").exists())

    def test_save_and_readonly(self):
        path = self.folder/"result.npz"
        self.result.save(path)
        with np.load(path, allow_pickle=False) as data:
            np.testing.assert_array_equal(data["electrons__primary_disc"],
                                          self.result.electrons["primary_disc"])
            metadata = json.loads(str(data["metadata_json"]))
            self.assertTrue(metadata["legacy_table_precision"])
            self.assertEqual(metadata["solver_cells"], 16)
            self.assertEqual(len(metadata["catalogue"]), 11)
        with self.assertRaises(ValueError):
            self.result.electrons["primary_disc"][0,0] = 0.

    def test_validation_and_recovery(self):
        for arguments in ((3,8,16), (16,8,501), (16.,8,16)):
            with self.assertRaises((TypeError,ValueError)):
                SolverGrid(*arguments)
        with self.assertRaises(ValueError):
            solve(self.p, self.grid, threads=0)
        with self.assertRaises(FileNotFoundError):
            solve(self.p, self.grid, library=self.folder/"absent.dylib")
        # Invalid native buffers fail before entering any OpenMP worker.
        lib = _load()
        self.assertEqual(lib.cg_transport(1,0,0,None,None,None,None,None,None), 1)
        # The coarse-grid halo solve produces non-finite integrands on this
        # reference catalogue. Fail safely, then prove a valid run still works.
        with Preparation(self.catalogue, Grid(8,8)) as p:
            with self.assertRaisesRegex(RuntimeError, "catalogue indices"):
                solve(p, SolverGrid(8,8,8), threads=4, legacy_table_precision=False)
        recovered = solve(self.p,self.grid,threads=4)
        np.testing.assert_array_equal(recovered.electrons["primary_disc"],
                                      self.result.electrons["primary_disc"])

    def test_private_gsl_and_serial_generation_not_exported(self):
        lib = _load()
        for name in ("gsl_set_error_handler", "init_do_2D_IC", "CRe_steadystate_solve"):
            with self.assertRaises(AttributeError):
                getattr(lib,name)

    def assert_valid_result(self, result, cosmic_rays, photons):
        """Check every source component for the expected shape and physically allowed values."""
        for group, names in FILES:
            arrays = getattr(result, group)
            self.assertEqual(len(arrays), len(names))
            for name, values in arrays.items():
                with self.subTest(group=group, component=name):
                    self.assertEqual(values.shape, (11, photons if group == "emission" else cosmic_rays))
                    self.assertTrue(np.isfinite(values).all())
                    self.assertTrue(np.all(values >= 0))

    def test_all_components_and_precision_modes(self):
        self.assert_valid_result(self.result, 16, 8)
        self.assert_valid_result(solve(self.p, self.grid, threads=4, legacy_table_precision=False), 16, 8)

    def test_denser_grid_all_components(self):
        with Preparation(self.catalogue, Grid(32, 32)) as p:
            result = solve(p, SolverGrid(32, 16, 32), threads=4)
        self.assert_valid_result(result, 32, 16)

    def test_preparation_library_override(self):
        with self.assertRaises(FileNotFoundError):
            Preparation(self.catalogue,Grid(16,16),library=self.folder/"absent.so")

    def test_backend_abi_and_error_validation(self):
        lib=_load()
        self.assertEqual(lib.cg_solver_abi(),5)
        self.assertTrue(lib.cg_solver_openmp_enabled())
        self.assertEqual(lib.cg_properties(1,0,None,None,None),1)
        self.assertEqual(lib.cg_spectra(1,None,None),1)

    def test_no_table_mutation(self):
        before={key:table.snapshot() for key,table in self.p._tables.items()}
        solve(self.p,self.grid,threads=4)
        self.assertEqual(before,{key:table.snapshot() for key,table in self.p._tables.items()})
