"""Python-owned configuration and serial ionisation-loss interface."""
import threading
from .physics.ionisation import ionisation


class Context:
    """Python configuration; no native allocation or library load on creation.

    Ionisation accepts total electron energies in GeV and density in cm^-3.
    The optional library path selects the solver, loaded only on demand.
    Querying openmp_enabled explicitly loads that solver to inspect its build.
    """
    def __init__(self, threads=1, library=None):
        """Validate worker configuration without loading or allocating the native backend."""
        if isinstance(threads, bool) or not isinstance(threads, int):
            raise TypeError("threads must be an integer")
        if not 1 <= threads <= 2147483647:
            raise ValueError("threads must fit a positive C int")
        self.threads = threads
        self.library = library
        self._lock = threading.RLock()
        self._closed = False

    @property
    def version(self):
        """Return the Python interface version."""
        return "0.4.0"

    @property
    def openmp_enabled(self):
        """Report whether the selected native library was built with OpenMP.

        This query loads the library on demand; constructing Context or evaluating
        serial ionisation losses does not require a compiled backend.
        """
        from .solver import _load
        return bool(_load(self.library).cg_solver_openmp_enabled())

    def _ensure_open(self):
        """Reject operations after this object has been closed."""
        if self._closed:
            raise RuntimeError("Context is closed")

    def ionisation_loss(self, energies_gev, density_cm3):
        """Return neutral-medium electron dE/dt values in GeV/s, negative for losses.

        energies_gev is an iterable of total electron energies; density_cm3 is a
        scalar gas number density. Returns a list in input order and rejects
        invalid inputs or use after close().
        """
        with self._lock:
            self._ensure_open()
            values = tuple(float(value) for value in energies_gev)
            return ionisation(values, float(density_cm3))

    def close(self):
        """Mark this Python configuration closed; no native allocation is released."""
        with self._lock:
            self._closed = True

    def __enter__(self):
        """Verify that the resource is open before entering its managed scope."""
        self._ensure_open()
        return self

    def __exit__(self, *args):
        """Close managed resources when leaving the scope, including on exceptions."""
        self.close()
