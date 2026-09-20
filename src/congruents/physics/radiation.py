"""Stellar, UV, FIR and CMB photon fields and energy-density diagnostics."""
import sys
import math
import numpy as np
from .. import constants as k
from ..quadrature import integrate

def temperature_grid(field, temperatures):
    """Return reference CMB/FIR temperature nodes in kelvin.

    Nominal steps determine the number of planes, then rounded endpoints set
    the actual spacing. A single plane cannot supply an interpolation bracket.
    """
    if field == "CMB":
        low, high, step = k.TCMB, temperatures[0], .5
        limits = (math.floor(k.TCMB*10)/10, math.ceil(high*10)/10)
    elif field == "FIR":
        low, high, step = temperatures[1], temperatures[2], 5.
        limits = (math.floor(low), math.ceil(high))
    else:
        raise ValueError("Temperature grid is only for CMB/FIR")
    n = math.ceil((high-low)/step)+1
    if n < 2 or n > 4096 or limits[0] <= 0 or limits[1] <= limits[0]:
        raise ValueError("Legacy temperature grid requires at least two planes")
    delta = (limits[1]-limits[0])/(n-1)
    return tuple([limits[0]] + [limits[0]+i*delta for i in range(1,n-1)] + [limits[1]])

def photon(field, temp, energy):
    """Undiluted differential photon density in cm^-3 GeV^-1.

    energy is in GeV and temp in kelvin; named stellar fields fix their own
    temperatures. UV is piecewise in wavelength (microns); FIR is a modified
    blackbody. Galaxy luminosity dilution is applied separately.
    """
    if field == "UV":
        wavelength = k.C*k.H/energy*1e4
        if .134 < wavelength <= .246:
            return 2.373/energy**2*k.ERG_GEV*wavelength**(-.6678)
        if .110 < wavelength <= .134:
            return 68.25/energy**2*k.ERG_GEV*wavelength
        if .0912 < wavelength <= .110:
            return 1.287e5/energy**2*k.ERG_GEV*wavelength**4.4172
        return 0.
    if field in ("3000","4000","7500"):
        temp = float(field)
    # Avoid exponential overflow in the Wien tail; its photon density is negligible.
    exponent = energy/(k.KB*temp)
    if exponent > math.log(sys.float_info.max):
        return 0.
    value = 8*math.pi*energy**2/(k.H*k.C)**3/(math.exp(exponent)-1)
    return value*energy/(2e12*k.H) if field == "FIR" else value

def dilution(galaxy, dust):
    """Dimensionless luminosity weights for 3000/4000/7500 K, UV and FIR.

    Use the reference 2*pi*Re^2*c geometry, not a spherical 4*pi convention.
    CMB needs no luminosity weight: its temperature follows galaxy redshift.
    """
    mass, radius, sfr = galaxy.stellar_mass_msun, galaxy.radius_kpc, galaxy.sfr_msun_per_year
    old = 10**(.8480565*math.log10(mass/.56)+1.521623)
    young = (10**(.7969616*math.log10(sfr/.56)+9.007323) if math.log10(sfr)>-2.6
             else 10**(.9867204*math.log10(sfr/.56)+9.476592))
    fir = 10**(1.096548*math.log10(sfr/.56)+9.710084)
    densities = (k.ARAD*3000**4,k.ARAD*4000**4,k.ARAD*7500**4,4450.1668,
                 24.8863*8*math.pi/((k.H*k.C)**3*2e12*k.H)*(k.KB*dust)**5)
    luminosities = (.574*old,.426*old,.763*young,.237*young,fir)
    return tuple(l*k.LSUN/(u*2*math.pi*(radius*1e3*k.PC)**2*k.C)
                 for l,u in zip(luminosities,densities))

def radiation(galaxy, properties, energies, config):
    """Return seven photon-density spectra and eight energy-density diagnostics.

    Spectra follow RADIATION_NAMES; energy densities are eV/cm^3 in the
    Urad_Ub output order. The two old-stellar diagnostic fields deliberately
    use unsplit luminosity, so summing those diagnostics is not the total field.
    """
    dust, cmb = properties["dust_temperature_k"], k.TCMB*(1+galaxy.redshift)
    d = dilution(galaxy,dust)
    def components(e):
        """Build the diluted photon-field components at one target energy."""
        c = photon("CMB",cmb,e)
        f = d[4]*photon("FIR",dust,e)
        a,b,s,v = (d[i]*photon(field,0,e) for i,field in enumerate(("3000","4000","7500","UV")))
        # Legacy diagnostic helpers use unsplit old luminosity.
        return (c,f,a/.574,b/.426,s,v,c+f+a+b+s+v)
    energies = tuple(float(x) for x in energies)
    if any(not math.isfinite(e) or e<=0 for e in energies):
        raise ValueError("Photon energies must be finite and positive")
    fields = tuple(zip(*(components(e) for e in energies))) if energies else ((),)*7
    # Energy density is integral E*n(E)dE; log-energy quadrature adds another E.
    # Multiplication by 1e9 converts GeV/cm^3 to the output eV/cm^3.
    integrals = [integrate(lambda x: math.exp(x)**2*components(math.exp(x))[j],
                           math.log(config[4]),math.log(config[5]),rtol=1e-6)*1e9 for j in range(7)]
    magnetic = properties["magnetic_field_gauss"]**2/(8*math.pi)/k.GEV_ERG*1e9
    return fields, (integrals[6],magnetic,*integrals[:6])
