"""Borrowed-array ABI for galaxy-parallel native physics (ABI 5).

This module contains declarations and validation, not scientific equations.
Libraries load only when a native operation is requested. All calls are
synchronous; NumPy owners stay alive until every OpenMP worker has returned.
"""
import ctypes as ct
from pathlib import Path
import sys
import threading
import numpy as np

LOCK = threading.RLock()
DP = ct.POINTER(ct.c_double)
IP = ct.POINTER(ct.c_int)


class TableInput(ct.Structure):
    """ctypes layout matching cg_table_input in csrc/solver.h; pointers are borrowed."""
    _fields_ = [("nx",ct.c_size_t),("ny",ct.c_size_t),
                ("x",DP),("y",DP),("values",DP),("temperature",ct.c_double)]


class RunInput(ct.Structure):
    """ctypes layout matching cg_run; field order and C types are ABI-sensitive."""
    _fields_ = ([(n,ct.c_size_t) for n in ("n","ne","np","ns","ncmb","nfir")]
                + [(n,DP) for n in ("rows","props","kinetic","electron","photon","cp",
                                    "transport","electrons","emission","internal")]
                + [(n,ct.POINTER(TableInput)) for n in ("ic","gamma","bs","sy")]
                + [("target_low",ct.c_double),("target_high",ct.c_double),("diagnostics",ct.c_int)]
                + [(n,DP) for n in ("loss","critical","budget","radio","escape")])


def pointer(array):
    """Return a borrowed C double pointer without copying or validating the buffer.

    The caller must supply a C-contiguous float64 NumPy array and keep its
    storage alive for the entire synchronous native call.
    """
    return array.ctypes.data_as(DP)


def descriptor(table, temperature=0.):
    """Describe a read-only table for C without copying its axes or values.

    Values have row-major shape (ny, nx). Temperature is in kelvin for CMB/FIR
    planes and zero otherwise. Keep the owning Table alive until C returns.
    """
    return TableInput(len(table.x),len(table.y),pointer(table.x),pointer(table.y),
                      pointer(table.values),temperature)


def load(path=None):
    """Load the selected native library, verify its ABI and declare ctypes signatures."""
    suffix = "dylib" if sys.platform == "darwin" else "so"
    path = Path(path) if path else Path(__file__).resolve().parents[2]/"build"/f"libcongruents_solver.{suffix}"
    if not path.is_file():
        raise FileNotFoundError(f"Build the OpenMP backend first: make solver ({path})")
    lib = ct.CDLL(str(path.resolve()))
    lib.cg_solver_abi.argtypes = []
    lib.cg_solver_abi.restype = ct.c_uint
    if lib.cg_solver_abi() != 5:
        raise RuntimeError("Unsupported solver ABI; rebuild the native library")
    lib.cg_solver_openmp_enabled.argtypes = []
    lib.cg_solver_openmp_enabled.restype = ct.c_int
    lib.cg_properties.argtypes = [ct.c_int,ct.c_size_t,DP,DP,IP]
    lib.cg_properties.restype = ct.c_int
    lib.cg_transport.argtypes = [ct.c_int,ct.c_size_t,ct.c_size_t,DP,DP,DP,DP,DP,IP]
    lib.cg_transport.restype = ct.c_int
    lib.cg_spectra.argtypes = [ct.c_int,ct.POINTER(RunInput),IP]
    lib.cg_spectra.restype = ct.c_int
    return lib


def check(code,statuses):
    """Raise on batch or per-galaxy failures so partial results cannot be accepted."""
    if code:
        raise RuntimeError(f"Native batch rejected input or resources (status {code})")
    failed = np.flatnonzero(statuses)
    if failed.size:
        raise RuntimeError(f"Native calculation failed for catalogue indices {failed.tolist()}: "
                           f"statuses {statuses[failed].tolist()}; outputs discarded")


def threads_value(threads):
    """Validate that the worker count fits a positive C integer."""
    if isinstance(threads,bool) or not isinstance(threads,int):
        raise TypeError("threads must be an integer")
    if not 1 <= threads <= 2147483647:
        raise ValueError("threads must be a positive C integer")
    return threads
