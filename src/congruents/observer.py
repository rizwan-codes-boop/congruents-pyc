"""Observer-frame conversion with explicit source/observer optical-depth axes."""
import numpy as np
from scipy.integrate import quad
from .constants import PC


def _redshifts(values):
    """Validate positive redshifts before any luminosity-distance division."""
    z = np.asarray(values, dtype=np.float64)
    if z.ndim != 1 or not z.size or not np.isfinite(z).all() or np.any(z <= 0):
        raise ValueError("Require a nonempty vector of finite positive redshifts")
    return z


def luminosity_distance_mpc(redshifts):
    """Return luminosity distances in Mpc for a vector of positive redshifts.

    Uses a flat cosmology with H0 = 70 km/s/Mpc, Omega_m = 0.3 and
    Omega_Lambda = 0.7, without radiation. Independent-distance inputs are
    not supported, so zero-redshift flux calculations are rejected.
    """
    z = _redshifts(redshifts)
    return np.array([(1+v)*(299792.458/70)*quad(
        lambda x: 1/np.sqrt(.3*(1+x)**3+.7), 0., v,
        epsabs=0., epsrel=1e-8)[0] for v in z])


def distance_factor_cm2(redshifts):
    """Legacy distmod=(1+z)^2/(4*pi*dL^2), in cm^-2 (not magnitudes)."""
    z = _redshifts(redshifts)
    return (1+z)**2/(4*np.pi*(luminosity_distance_mpc(z)*1e6*PC)**2)


def observer_sed(energy_gev, source_rate, redshifts, *, tau_internal, tau_ebl):
    """One component [galaxy, energy] -> E² dN/dE, GeV cm^-2 s^-1.

    Evaluate source rate and internal tau at E_source=(1+z)*E_observer.
    EBL tau must already be evaluated on the observer energy axis. Linear
    interpolation is used and spectra are zero outside the source domain.
    Tau arrays are required, finite, nonnegative and exactly the source shape;
    passing explicit zeros is an unattenuated diagnostic, not final output.
    No input array is modified. No C call, native loop, or new physics is added.
    """
    e = np.asarray(energy_gev, dtype=np.float64)
    z = _redshifts(redshifts)
    if (e.ndim != 1 or len(e) < 2 or not np.isfinite(e).all() or
            np.any(e <= 0) or np.any(np.diff(e) <= 0)):
        raise ValueError("Require finite positive increasing photon energies")
    arrays = []
    for name, value in (("source_rate",source_rate), ("tau_internal",tau_internal),
                        ("tau_ebl",tau_ebl)):
        a = np.asarray(value, dtype=np.float64)
        if a.shape != (len(z),len(e)) or not np.isfinite(a).all() or np.any(a < 0):
            raise ValueError(f"{name} must be a finite nonnegative [galaxy, energy] array")
        arrays.append(a)
    source, internal, ebl = arrays
    shifted = np.array([np.interp(e*(1+v),e,row,left=0.,right=0.)
                        for v,row in zip(z,source)])
    shifted_internal = np.array([np.interp(e*(1+v),e,row)
                                 for v,row in zip(z,internal)])
    result = shifted*np.exp(-shifted_internal-ebl)*distance_factor_cm2(z)[:,None]*e**2
    if not np.isfinite(result).all():
        raise ValueError("Observer spectrum overflowed; inputs are outside supported range")
    result.flags.writeable = False
    return result
