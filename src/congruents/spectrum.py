"""Legacy log-energy quadrature of linearly interpolated differential spectra."""
import math
import numpy as np
from .quadrature import integrate


def integrate_spectrum(energy, rate, low=None, high=None):
    """Integrate E * rate(E) dE over a bounded GeV energy axis.

    For rate in GeV^-1 s^-1 the result is GeV/s. With x=log(E), E^2
    appears from the energy weighting and dE=E dx. Interpolation of rate
    remains linear in E, not log-log, and bounds cannot extrapolate.
    """
    e, y = np.asarray(energy,float), np.asarray(rate,float)
    if (e.ndim!=1 or len(e)<2 or y.shape!=e.shape or not np.isfinite(e).all() or
            not np.isfinite(y).all() or np.any(e<=0) or np.any(np.diff(e)<=0)):
        raise ValueError("Invalid spectrum")
    low, high = e[0] if low is None else low, e[-1] if high is None else high
    if not e[0]<=low<high<=e[-1]:
        raise ValueError("Integration bounds outside spectrum")
    return integrate(lambda x: math.exp(x)**2*float(np.interp(math.exp(x),e,y)),
                     math.log(low),math.log(high),rtol=1e-6)
