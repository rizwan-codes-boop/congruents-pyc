"""Native galaxy-property interface in the improvised legacy code's native units."""
import math
import numpy as np
from . import constants as k

PROPERTY_NAMES = ("height_pc", "density_cm3", "magnetic_field_gauss",
                  "gas_dispersion_kms", "area_pc2", "gas_surface_msun_pc2",
                  "sfr_surface_msun_yr_pc2", "stellar_surface_msun_pc2",
                  "dust_temperature_k", "halo_magnetic_field_gauss")


def properties(rows, threads=1, library=None):
    """Return a read-only (n_galaxies, 10) array of native galaxy properties.

    Input rows are (redshift, stellar mass in Msun, radius in kpc, SFR in
    Msun/yr). Output column names and units are given by PROPERTY_NAMES.
    Python owns the contiguous buffers; C evaluates the property equations
    with galaxy-only OpenMP parallelism. Any failed row raises RuntimeError.
    """
    from ._native import LOCK, IP, load, pointer, check, threads_value
    threads_value(threads)
    rows = np.ascontiguousarray(rows,dtype=np.float64)
    if (rows.ndim != 2 or rows.shape[1] != 4 or not 1 <= len(rows) <= 100000
            or not np.isfinite(rows).all() or np.any(rows[:,0]<0)
            or np.any(rows[:,0]>20) or np.any(rows[:,1:]<=0)):
        raise ValueError("Require 1..100000 valid four-column galaxy rows")
    out = np.empty((len(rows),10))
    status = np.zeros(len(rows),dtype=np.int32)
    with LOCK:
        lib=load(library)
        check(lib.cg_properties(threads,len(rows),pointer(rows),pointer(out),status.ctypes.data_as(IP)),status)
    out.flags.writeable=False
    return out
