"""Bremsstrahlung screened differential cross-section tables, mb/GeV."""
import math
import numpy as np
from .. import constants as k
from ..grids import log_grid

def generate(nx, ny, config):
    """Return a screened differential cross-section plane in mb/GeV.

    x contains emitted photon energies and y total electron energies, both in
    GeV; values has shape (ny, nx). Entries with photon energy at or above
    electron energy are zero. config uses tables.generate bound conventions;
    validation is performed by that dispatcher.
    """
    x = log_grid(nx,config[0],config[1])
    y = log_grid(ny,config[2],config[3])
    z = np.zeros((ny,nx))
    # Reference bremsstrahlung screening functions: interpolate the tabulated
    # branch in delta, then use the analytic branch for delta > 2.
    delta_axis = [0,.01,.02,.05,.1,.2,.5,1,2,5,10]
    phi1 = [45.79,45.43,45.09,44.11,42.64,40.16,34.97,29.97,24.73,18.09,13.65]
    phi2 = [44.46,44.38,44.24,43.65,42.49,40.19,34.93,29.78,24.34,17.28,12.41]
    for i, out in enumerate(x):
        for j, electron in enumerate(y):
            out,electron = float(out),float(electron)
            if out >= electron:
                continue
            delta = out*k.ME/(4*k.ALPHA*electron*(electron-out))
            if delta > 2:
                a=b=4*(math.log(2*electron/k.ME*((electron-out)/out))-.5)
            else:
                a,b = float(np.interp(delta,delta_axis,phi1)),float(np.interp(delta,delta_axis,phi2))
            z[j,i] = 3/(8*math.pi)*k.SIGMA_MB*k.ALPHA*max(
                (1+(1-out/electron)**2)*a-2/3*(1-out/electron)*b,0)/out
    return x,y,z
