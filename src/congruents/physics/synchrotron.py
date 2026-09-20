"""Synchrotron kernel table: x integral K_(5/3), dimensionless.

The finite upper integration cutoff preserves the reference convention.
"""
import math
import numpy as np
from .. import constants as k
from ..quadrature import integrate
from scipy.special import kv
from ..grids import log_grid

def generate(nx, config):
    """Return the dimensionless synchrotron kernel on nx logarithmic samples.

    Uses config[6:8] for the x range. The Bessel integral ends at the smaller
    of the upper x bound and 150; samples beyond this cutoff are zero.
    Returns x, a dummy singleton y axis, and values with shape (1, nx).
    """
    x = log_grid(nx,config[6],config[7])
    end = min(math.log(config[7]),math.log(150.))
    values = [a*integrate(lambda z: math.exp(z)*float(kv(5./3.,math.exp(z))),
                         math.log(a),end) if math.log(a)<end else 0. for a in x]
    return x,np.array([1.]),np.array(values).reshape(1,nx)
