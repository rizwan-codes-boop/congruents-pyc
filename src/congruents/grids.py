"""Numerical grid construction shared by preparation and source calculations."""
import math
import numpy as np

def log_grid(n, low, high):
    """Return n logarithmically spaced float64 samples, preserving scalar arithmetic.

    This is the public/table grid; the electron solver has its own explicitly
    fused grid construction in the native steady-state module.
    """
    if isinstance(n, bool) or not isinstance(n, int) or n < 2:
        raise ValueError("Grid size must be an integer >= 2")
    if not 0 < low < high or not math.isfinite(high):
        raise ValueError("Grid bounds must be finite, positive and increasing")
    lo = math.log(low)
    step = (math.log(high)-lo)/(n-1)
    return np.array([math.exp(lo+step*i) for i in range(n)], dtype=np.float64)
