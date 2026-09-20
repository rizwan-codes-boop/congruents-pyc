"""Plot completed model exports in the Figure 9 panel arrangement.

This script displays model components, not digitised paper data or a fitted
comparison. It reads/writes the working directory's output folder and opens
a Matplotlib window when executed. See README for selecting a run directory.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# Paths are relative to the working directory, not this script. Point output
# at the desired completed run before plotting; no spectra are computed here.
OUT = Path("output")


def load_vector(filename):
    """Read an exported axis after its one-line text header."""
    return np.loadtxt(OUT / filename, skiprows=1)


def load_spectra(filename):
    """Read observer-frame E^2-weighted spectra, with catalogue rows unchanged."""
    return np.loadtxt(OUT / filename, skiprows=1)


energy = load_vector("E_gam.txt")

pion = load_spectra("spec_pi.txt")
brems_primary = load_spectra("spec_BS_1_z1.txt")
brems_secondary = load_spectra("spec_BS_2_z1.txt")

ic_primary = (
    load_spectra("spec_IC_1_z1.txt")
    + load_spectra("spec_IC_1_z2.txt")
)
ic_secondary = (
    load_spectra("spec_IC_2_z1.txt")
    + load_spectra("spec_IC_2_z2.txt")
)

sync_primary = (
    load_spectra("spec_SY_1_z1.txt")
    + load_spectra("spec_SY_1_z2.txt")
)
sync_secondary = (
    load_spectra("spec_SY_2_z1.txt")
    + load_spectra("spec_SY_2_z2.txt")
)

freefree = load_spectra("spec_FF.txt")

# Sum physical photon components only. Neutrinos and the full-calorimetry
# comparison curve are excluded to avoid mixing channels or double counting.
total = (
    pion
    + brems_primary
    + brems_secondary
    + ic_primary
    + ic_secondary
    + sync_primary
    + sync_secondary
    + freefree
)

# Each tuple is (panel label, zero-based catalogue row, displayed log sSFR).
# The sSFR text is an annotation, not a rescaling of the computed spectra.
galaxies = [
    ("NGC 2146", 4, -9.28),
    ("NGC 2403", 3, -9.59),
    ("M82", 5, -9.78),
    ("NGC 253", 2, -9.73),
    ("SMC", 9, -10.1),
    ("M33", 7, -10.4),
]

# Limits read from the three panel rows of MNRAS 523, 2608, Fig. 9.
X_LIMITS = (1.0e-16, 1.0e5)
ENERGY_TICKS = [1e-15, 1e-12, 1e-9, 1e-6, 1e-3, 1e0, 1e3]
# Figure 9 shares the vertical scale within each row. Consequently, M33's
# synchrotron components lie mostly below the displayed bottom-row range.
ROW_Y_LIMITS = (
    (3.0e-13, 3.0e-9),
    (3.0e-13, 3.0e-9),
    (3.0e-12, 3.0e-8),
)
ROW_Y_TICKS = (
    (1.0e-12, 1.0e-11, 1.0e-10, 1.0e-9),
    (1.0e-12, 1.0e-11, 1.0e-10, 1.0e-9),
    (1.0e-11, 1.0e-10, 1.0e-9, 1.0e-8),
)

PLANCK_GEV_S = 4.135667696e-24


def energy_to_frequency_ghz(energy_gev):
    """Convert GeV to GHz through E=h*nu, for display only."""
    return energy_gev / (PLANCK_GEV_S * 1.0e9)


def frequency_ghz_to_energy(frequency_ghz):
    """Inverse display transform from GHz to GeV."""
    return frequency_ghz * 1.0e9 * PLANCK_GEV_S

fig, axes = plt.subplots(
    3, 2,
    figsize=(9.2, 6.8),
    sharex=True,
    sharey="row",
    gridspec_kw={"hspace": 0.0, "wspace": 0.0},
)

for panel, (ax, (name, index, log_ssfr)) in enumerate(zip(axes.flat, galaxies)):
    row = panel // 2
    ax.loglog(energy, total[index], color="C0", lw=2, label="Total")
    ax.loglog(energy, pion[index], "--", color="C0", label=r"$\pi^0$ decay")

    ax.loglog(
        energy, brems_primary[index],
        "--", color="#b5a800", label="Brems. primary"
    )
    ax.loglog(
        energy, brems_secondary[index],
        "-.", color="#b5a800", label="Brems. secondary"
    )

    ax.loglog(
        energy, ic_primary[index],
        "--", color="magenta", label="IC primary"
    )
    ax.loglog(
        energy, ic_secondary[index],
        "-.", color="magenta", label="IC secondary"
    )

    ax.loglog(
        energy, sync_primary[index],
        "--", color="#8b4935", label="Synch. primary"
    )
    ax.loglog(
        energy, sync_secondary[index],
        "-.", color="#8b4935", label="Synch. secondary"
    )

    ax.loglog(
        energy, freefree[index],
        "--", color="orangered", label="Thermal free-free"
    )

    ax.axvline(
        PLANCK_GEV_S * 1.49e9,
        color="red",
        ls="--",
        lw=1,
    )

    ax.text(
        PLANCK_GEV_S * 1.49e9,
        0.97,
        "1.49 GHz",
        transform=ax.get_xaxis_transform(),
        color="red",
        fontsize=7,
        ha="right",
        va="top",
        rotation=90,
    )

    ax.text(0.46, 0.92, name, transform=ax.transAxes,
            fontsize=9, fontweight="bold", va="top")
    ax.text(0.46, 0.79, rf"$\log(\mathrm{{sSFR}})={log_ssfr:g}$",
            transform=ax.transAxes, fontsize=8, va="top")

    ax.set_xlim(*X_LIMITS)
    ax.set_xticks(ENERGY_TICKS)
    ax.set_ylim(*ROW_Y_LIMITS[row])
    ax.set_yticks(ROW_Y_TICKS[row])
    ax.tick_params(which="both", direction="in", top=True, right=True)

    # Match the paper: numerical y labels appear only on the left column.
    if panel % 2:
        ax.tick_params(axis="y", labelleft=False)

# The paper shows frequency only above the top row. These ticks are equivalent
# to the energy scale through E = h nu.
for ax in axes[0]:
    frequency_axis = ax.secondary_xaxis(
        "top",
        functions=(energy_to_frequency_ghz, frequency_ghz_to_energy),
    )
    frequency_axis.set_xscale("log")
    frequency_axis.set_xticks([1e1, 1e4, 1e7, 1e10, 1e13, 1e16, 1e19])
    frequency_axis.set_xlabel(r"$\nu\ [\mathrm{GHz}]$", labelpad=1)
    frequency_axis.tick_params(which="both", direction="in", labelsize=8)

for ax in axes[-1]:
    ax.set_xlabel(r"$E_\gamma\ [\mathrm{GeV}]$")

for ax in axes[:, 0]:
    ax.set_ylabel(
        r"$E_\gamma^2\phi_\gamma$"
        r"$\ [\mathrm{GeV\,s^{-1}\,cm^{-2}}]$"
    )

handles, labels = axes[0, 0].get_legend_handles_labels()
fig.legend(
    handles,
    labels,
    loc="lower center",
    ncol=3,
    fontsize=7,
    frameon=True,
)

fig.subplots_adjust(left=0.10, right=0.99, top=0.91, bottom=0.18)
fig.savefig("output/mnras2023_figure9_model.pdf")
fig.savefig("output/mnras2023_figure9_model.png", dpi=300)

plt.show()
