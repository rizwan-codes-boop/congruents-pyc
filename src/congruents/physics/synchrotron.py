"""Isotropic pitch-angle averaged synchrotron kernel R(x), dimensionless.

The finite upper integration cutoff preserves the reference convention.
"""
import math
import numpy as np
from .. import constants as k
from ..quadrature import integrate
from scipy.special import kv
from ..grids import log_grid

def kernel(x):
    """Return R(x) for x=nu/nu_c(sin(theta)=1).

    Swapping the integrals in int_0^(pi/2) sin(theta)^2 F(x/sin(theta))
    dtheta gives x int_x^infinity K_(5/3)(t) sqrt(1-(x/t)^2) dt.
    Its integral over x is 16*pi/(27*sqrt(3)), giving the isotropic
    4/3 Thomson cooling coefficient with the native emission prefactor.
    The exponential tail beyond t=150 is neglected.
    """
    if not math.isfinite(x) or x <= 0:
        raise ValueError("Require a finite positive synchrotron coordinate")
    if x >= 150.:
        return 0.
    low = math.log(x)
    return x*integrate(lambda z: math.exp(z)*float(kv(5./3., math.exp(z)))
                      *math.sqrt(max(0., -math.expm1(2*(low-z)))),
                      low, math.log(150.))

def generate(nx, config):
    """Return the dimensionless synchrotron kernel on nx logarithmic samples.

    Uses config[6:8] for the x range. The Bessel integral ends at 150;
    samples beyond this cutoff are zero.
    Returns x, a dummy singleton y axis, and values with shape (1, nx).
    """
    x = log_grid(nx,config[6],config[7])
    values = [kernel(float(a)) for a in x]
    return x,np.array([1.]),np.array(values).reshape(1,nx)
