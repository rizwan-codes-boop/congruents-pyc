"""Inverse Compton emission and energy-transfer tables, s^-1 GeV^-1.

Target radiation is supplied by radiation.photon; only the requested table
kind selects the scattering kernel. Energy arguments are GeV.
"""
import math
import numpy as np
from .. import constants as k
from ..quadrature import integrate
from ..grids import log_grid
from .radiation import photon

def _ic_kernel(electron, target, emitted, density):
    """Dimensionless scattering kernel times target photon density.

    All energies are GeV; q bounds enforce the reference scattering domain.
    The cross-section/rate prefactor is applied by generate, not here.
    """
    if electron == emitted:
        return 0.
    gamma = 4*target/k.ME*electron/k.ME
    q = emitted/(gamma*(electron-emitted))
    if (k.ME/(2*electron))**2 < q < 1:
        return density*(2*q*math.log(q)+(1+2*q)*(1-q)+(gamma*q)**2*(1-q)/(2*(1+gamma*q)))
    return 0.

def generate(kind, field, temperature, nx, ny, config):
    """Generate an IC rate plane in s^-1 GeV^-1 by serial photon integration.

    kind is emission (x = emitted photon energy) or gamma (x = electron energy
    transfer). The y axis is initial total electron energy; all energies are
    GeV. Return x, y, values with values.shape == (ny, nx). field identifies
    the target spectrum and temperature is in kelvin for CMB/FIR.
    config uses the bound pairs documented by tables.generate, which validates
    requests before dispatch. Forbidden transitions remain zero.
    """
    x = log_grid(nx,config[4] if kind=="gamma" else config[0],
                 config[3] if kind=="gamma" else config[1])
    y = log_grid(ny,config[2],config[3])
    z = np.zeros((ny,nx))
    for i, out in enumerate(x):
        for j, electron in enumerate(y):
            out,electron = float(out),float(electron)
            if kind == "emission":
                # C fmax/fmin ignore NaN bounds when photon >= electron;
                # the kernel itself is zero there.
                if out >= electron:
                    continue
                low=max(math.log(config[4]),math.log(out*k.ME**2/(4*electron*(electron-out))))
                high=min(math.log(config[5]),math.log(out*electron/(electron-out)))
            else:
                if out >= electron:
                    continue
                low,high = math.log(config[4]),math.log(config[5])
            if low < high:
                # Evaluate the target-photon scattering contribution at one log-energy sample.
                def integrand(logtarget):
                    target=math.exp(logtarget)
                    # For a transition table, DeltaE = E_gamma - E_target. Emission tables
                    # instead hold E_gamma fixed. Do not select this branch by field index.
                    emitted=out+target if kind=="gamma" else out
                    return _ic_kernel(electron,target,emitted,photon(field,temperature,target))
                result=integrate(integrand,low,high)
                z[j,i]=.75*k.SIGMA_MB*k.MB_CM2*k.C/(electron/k.ME)**2*max(result,0)
    return x,y,z
