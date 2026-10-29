"""Exact energy moments of piecewise-linear differential spectra."""
import numpy as np


def integrate_spectrum(energy, rate, low=None, high=None):
    """Integrate E * rate(E) dE over a bounded GeV energy axis.

    For rate in GeV^-1 s^-1 the result is GeV/s. Integrate each linear
    segment analytically; interpolation is in E, not log-log, and bounds
    cannot extrapolate.
    """
    e, y = np.asarray(energy,float), np.asarray(rate,float)
    if (e.ndim!=1 or len(e)<2 or y.shape!=e.shape or not np.isfinite(e).all() or
            not np.isfinite(y).all() or np.any(e<=0) or np.any(np.diff(e)<=0)):
        raise ValueError("Invalid spectrum")
    low, high = e[0] if low is None else low, e[-1] if high is None else high
    if not e[0]<=low<high<=e[-1]:
        raise ValueError("Integration bounds outside spectrum")
    # Integrate the piecewise-linear representation exactly, splitting at knots.
    x=np.r_[low,e[(e>low)&(e<high)],high]
    value=np.interp(x,e,y)
    return float(np.sum(np.diff(x)/6*((2*x[:-1]+x[1:])*value[:-1]+
                                     (x[:-1]+2*x[1:])*value[1:]),dtype=np.longdouble))
