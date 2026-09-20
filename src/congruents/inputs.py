"""Named catalogue inputs; native units are explicit, never guessed."""
from dataclasses import dataclass
import math
from pathlib import Path

@dataclass(frozen=True)
class Galaxy:
    """One immutable catalogue case in the explicit units named by its fields.

    Redshift zero is allowed for source calculations; observer fluxes need z > 0.
    """
    redshift: float
    stellar_mass_msun: float
    radius_kpc: float
    sfr_msun_per_year: float

    def __post_init__(self):
        """Enforce finite redshift and positive mass, radius and star-formation rate."""
        values = self.as_row()
        if not all(math.isfinite(x) for x in values):
            raise ValueError("Galaxy fields must be finite")
        if not 0 <= values[0] <= 20 or any(x <= 0 for x in values[1:]):
            raise ValueError("Require 0 <= z <= 20 and positive mass, radius and SFR")

    def as_row(self):
        """Return (redshift, mass [Msun], effective radius [kpc], SFR [Msun/yr])."""
        return (float(self.redshift), float(self.stellar_mass_msun),
                float(self.radius_kpc), float(self.sfr_msun_per_year))

    @classmethod
    def from_quantities(cls, redshift, stellar_mass, radius, sfr):
        """Accept Astropy quantities; incompatible units raise UnitConversionError."""
        from astropy import units as u
        return cls(float(redshift), stellar_mass.to_value(u.Msun),
                   radius.to_value(u.kpc), sfr.to_value(u.Msun / u.yr))

@dataclass(frozen=True)
class Grid:
    """Lookup-plane dimensions, distinct from the particle and solver grids.

    nx indexes emitted photon energy or energy transfer; ny indexes total
    electron energy. Synchrotron uses nx points and a singleton y axis.
    """
    nx: int = 1000
    ny: int = 1000

    def __post_init__(self):
        """Reject noninteger dimensions and lookup planes exceeding the configured limits."""
        for value in (self.nx, self.ny):
            if isinstance(value, bool) or not isinstance(value, int):
                raise TypeError("Grid sizes must be integers")
            if not 2 <= value <= 4096:
                raise ValueError("Grid sizes must be between 2 and 4096")
        if self.nx * self.ny > 4194304:
            raise ValueError("At most 4194304 cells per plane")

def load_catalogue(path):
    """Read the counted, four-column catalogue without sorting or filtering rows.

    Order is scientifically significant: transport retains a special index-10
    calorimetry factor. Malformed counts, headers and nonphysical inputs fail early.
    """
    lines = Path(path).read_text().splitlines()
    if len(lines) < 3 or lines[0].strip() != "n_gal":
        raise ValueError("Expected n_gal, count, column header, then four-column rows")
    try:
        count = int(lines[1])
    except ValueError as exc:
        raise ValueError("Invalid galaxy count") from exc
    if lines[2].split() != ["z", "Mstar__Msol", "Re__kpc", "SFR__Msolyrm1"]:
        raise ValueError("Unexpected catalogue columns/units")
    rows = [line.split() for line in lines[3:] if line.strip()]
    if not 1 <= count <= 100000 or len(rows) != count or any(len(r) != 4 for r in rows):
        raise ValueError("Galaxy count or four-column row shape mismatch")
    return tuple(Galaxy(*(float(x) for x in row)) for row in rows)

def write_catalogue(path, galaxies):
    """Write the counted four-column catalogue, preserving order and precision.

    galaxies must be a nonempty iterable of Galaxy objects. The destination is
    overwritten if it exists; parent directories are not created.
    """
    galaxies = tuple(galaxies)
    if not galaxies:
        raise ValueError("Catalogue must not be empty")
    body = ["n_gal", str(len(galaxies)), "z Mstar__Msol Re__kpc SFR__Msolyrm1"]
    body.extend(" ".join(format(x, ".17g") for x in g.as_row()) for g in galaxies)
    Path(path).write_text("\n".join(body) + "\n")
