"""Serial Python attenuation, retaining improvised legacy code table and kernel rules."""
import hashlib
import math
from pathlib import Path
import numpy as np
from . import constants as k


class EBLTable:
    """Validated external optical-depth grid: redshift rows and energy columns.

    The file energy axis is eV, unlike model GeV axes. The content checksum
    is recorded in run metadata so changing EBL input is visible.
    """
    def __init__(self, path):
        """Read, validate and fingerprint the EBL optical-depth grid."""
        self.path = Path(path)
        raw = self.path.read_bytes()
        self.sha256 = hashlib.sha256(raw).hexdigest()
        tokens = [t for line in raw.decode().splitlines()
                  if not line.lstrip().startswith("#") for t in line.split()]
        if len(tokens)<2:
            raise ValueError("Missing EBL dimensions")
        nz, ne = int(tokens[0]), int(tokens[1])
        if nz < 2 or ne < 2 or len(tokens) != 2+nz+ne+nz*ne:
            raise ValueError("Invalid EBL table dimensions/token count")
        data = np.array(tokens[2:], dtype=float)
        self.z, self.energy = data[:nz], data[nz:nz+ne]
        self.values = data[nz+ne:].reshape(nz,ne)
        if (not np.isfinite(data).all() or np.any(np.diff(self.z)<=0) or
                np.any(np.diff(self.energy)<=0) or np.any(self.energy<=0) or
                np.any(self.z<0) or np.any(self.values<0)):
            raise ValueError("Invalid EBL coordinates/values")
        for a in (self.z,self.energy,self.values):
            a.flags.writeable = False

    def optical_depth(self, source_energy_gev, redshifts):
        """Bilinear extrapolation, E/(1+z) in eV capped at 1e15; clip tau>=0.

        This deliberately preserves spectra.c's source-grid indexing rather
        than resampling optical depths alongside the emitted spectrum.
        """
        e, z = np.asarray(source_energy_gev,float), np.asarray(redshifts,float)
        if (e.ndim!=1 or z.ndim!=1 or not np.isfinite(e).all() or
                not np.isfinite(z).all() or np.any(e<=0) or np.any(z<0)):
            raise ValueError("Invalid EBL query")
        x = np.minimum(e[None,:]/(1+z[:,None])*1e9,1e15)
        # Clipping the cell indices chooses an edge cell for extrapolation; u/v
        # are not clipped. This differs intentionally from zero-outside IC/BS tables.
        ix = np.clip(np.searchsorted(self.energy,x,side="right")-1,0,len(self.energy)-2)
        iy = np.clip(np.searchsorted(self.z,z,side="right")-1,0,len(self.z)-2)[:,None]
        u = (x-self.energy[ix])/(self.energy[ix+1]-self.energy[ix])
        v = (z[:,None]-self.z[iy])/(self.z[iy+1]-self.z[iy])
        a = self.values
        return np.maximum(0.,(1-u)*(1-v)*a[iy,ix]+u*(1-v)*a[iy,ix+1]+
                          (1-u)*v*a[iy+1,ix]+u*v*a[iy+1,ix+1])
