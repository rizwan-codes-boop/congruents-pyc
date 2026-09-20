"""Shared table bounds, validation and dispatch to individual physical processes."""
import math
from . import constants as k
from .physics import inverse_compton, bremsstrahlung, synchrotron

def bounds(galaxies, properties):
    """Derive shared table bounds from all catalogue cases.

    Return config=(photon min/max, electron min/max, target min/max, SY x min/max)
    and temperature extrema. Energies are GeV; synchrotron x is dimensionless.
    Bounds include both disc and halo magnetic-field ranges.
    """
    maxb = max(p["magnetic_field_gauss"] for p in properties)
    minb = min(p["halo_magnetic_field_gauss"] for p in properties)
    maxz = max(g.redshift for g in galaxies)
    lowe, highe = 1e-3+k.ME, 1e8+k.ME
    common = 2*math.pi**2*k.ME_G**2*k.C**3
    sx0 = common/(3*k.E_ESU*maxb*k.H_ERG)*(1e-16*k.ME)/highe**2*.1
    sx1 = common/(3*k.E_ESU*minb/10*k.H_ERG)*(1e8*(1+maxz)*k.ME)/lowe**2*10
    config = (1e-16, 1e8, lowe, highe, 1.59362*k.KB*k.TCMB*1e-4, 1e-7, sx0, sx1)
    ts = (k.TCMB*(1+maxz), min(p["dust_temperature_k"] for p in properties),
          max(p["dust_temperature_k"] for p in properties))
    return config, ts

def generate(kind, field, temperature, nx, ny, config):
    """Validate a table request and return NumPy x, y and values arrays.

    config contains min/max pairs for photon, total electron and target-photon
    energies in GeV, followed by dimensionless synchrotron x bounds.
    Kinds emission and gamma produce IC rates in s^-1 GeV^-1; gamma denotes
    electron energy transfer, not a separate gamma-ray emission component.
    Kind bs produces mb/GeV and sy a dimensionless kernel. Values have shape
    (ny, nx), except sy, whose y axis is a singleton and values shape (1, nx).
    CMB/FIR temperatures are in kelvin. Generation is serial Python.
    """
    from .inputs import Grid
    Grid(nx, ny)
    if kind not in ("emission","gamma","bs","sy"):
        raise ValueError("Unknown table kind")
    if field not in ("3000", "4000", "7500", "UV", "CMB", "FIR"):
        raise ValueError("Unknown radiation field")
    if field in ("CMB","FIR") and (not math.isfinite(temperature) or not 2<=temperature<=1000):
        raise ValueError("Invalid field temperature")
    if len(config)!=8 or any(not math.isfinite(v) or v<=0 for v in config):
        raise ValueError("Invalid table bounds")
    if any(config[i]>=config[i+1] for i in (0,2,4,6)):
        raise ValueError("Table bounds must increase")
    if kind == "sy":
        return synchrotron.generate(nx, config)
    if kind == "bs":
        return bremsstrahlung.generate(nx, ny, config)
    return inverse_compton.generate(kind, field, temperature, nx, ny, config)
