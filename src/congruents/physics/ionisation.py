"""Standalone neutral-medium electron ionisation losses; GeV/s."""
import math
import numpy as np
from .. import constants as k
from ..quadrature import integrate

def ionisation(energies, density):
    """Return negative dE/dt in GeV/s for total electron energies in GeV.

    Density is hydrogen number density in cm^-3. The reference neutral-medium
    mixture and excitation energies are fixed, not inferred from Astropy.
    """
    if not math.isfinite(density) or density < 0:
        raise ValueError("density must be finite and non-negative")
    out = []
    for e in energies:
        if not math.isfinite(e) or e <= k.ME:
            raise ValueError("Expected finite TOTAL electron energy above rest mass")
        v = -9./4.*k.C*k.SIGMA_MB*k.MB_CM2*k.ME*density*1.1*(
            math.log(e/k.ME)+.91*2./3.*math.log(k.ME*1e9/15.)+
            2.*.09*2./3.*math.log(k.ME*1e9/41.5))
        if not math.isfinite(v):
            raise RuntimeError("Non-finite ionisation result")
        out.append(v)
    return out
