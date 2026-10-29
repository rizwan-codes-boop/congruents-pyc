# CONGRUENTS

CONGRUENTS calculates cosmic-ray transport and multiwavelength emission from
galaxies. It produces source spectra, observer-frame spectra, optical depths
and energy-loss diagnostics. The implementation originated from the improvised
legacy code; the target scientific specification is Roth et al. (2023).

## Scientific reference policy

The immediate objective is a mathematically consistent interpretation of the
published model, including its assumptions and approximations. Retain changes
supported by a derivation and independent tests, rather than reverting them
solely to reproduce a suspected printed error. Do not change normalization merely
to fit a published curve. If published statements conflict or a normalization
convention is unresolved, document that uncertainty without claiming equivalence.
Later-paper changes must be identified separately from the Roth et al. (2023) model.

The current local scientific-audit revision is **not yet a paper-matched
release**. Its revised cooling stencil follows the linear reconstruction in B4;
its synchrotron emission uses an explicit isotropic pitch-angle average. Both
have independent analytic checks and are retained, but are not literal copies
of B6 and equation 42 respectively.
Numerical safeguards do not establish fidelity to the publication. Development
analyses and test programs are maintained outside this production folder.

Run the commands below from the repository root, with your Python environment
activated.

## Installation

Requirements:

- Python 3.9 or newer, NumPy and SciPy.
- A C11 compiler with OpenMP support, GSL and cubature.
- Matplotlib for plotting; Astropy is optional for input-unit conversion.

Install the Python package and plotting dependency:

```sh
python -m pip install -e .
python -m pip install matplotlib
# Optional Astropy support:
python -m pip install -e '.[units]'
```

Use an editable source installation; this repository does not supply a standalone
binary wheel. The Python installation does not compile the native library. Build it using
one of the configurations below, and rebuild after changing C source files.

### Building the native library

Provide GSL and cubature static libraries compiled with position-independent
code. For GSL, use `--with-pic --disable-shared` when configuring its build;
compile cubature with `-fPIC`.

Replace the paths below with your installation locations:

```sh
make solver CC=gcc \
  CPPFLAGS="-I/GSL_PREFIX/include -I/CUBATURE_PREFIX/include" \
  LDLIBS="/GSL_PREFIX/lib/libgsl.a /GSL_PREFIX/lib/libgslcblas.a /CUBATURE_PREFIX/lib/libcubature.a -lm"
```

On macOS, select an OpenMP-enabled GNU compiler rather than Apple Clang; for
example, a Homebrew GCC 14 installation can be selected with
`CC="$(brew --prefix gcc@14)/bin/gcc-14"`. The libraries must match your machine's
architecture. On an HPC system, load the required compiler/library modules and
use no more threads than your job allocation permits.

The build creates `build/libcongruents_solver.dylib` on macOS or
`build/libcongruents_solver.so` on Linux. Keep this file: the model needs it.

Local dependencies and compiled binaries are not tracked by Git. A fresh clone
on another machine needs the dependency installation described above; it does
not need another CONGRUENTS repository.

## Run the model

### Small-grid trial

This small-grid run checks that the program works. Do not use its spectra for
scientific interpretation.

```sh
PYTHONPATH=src python -u -m congruents input/cat_nt.txt output/smoke \
  --threads 4 --tables 32 --cosmic-rays 128 --photons 32 --cells 64 \
  --allow-unconverged
```

### Production run

```sh
PYTHONPATH=src python -u -m congruents input/cat_nt.txt output/production \
  --threads 4 --cache data/production-cache
```

The output directory must not already exist. Files are written after the
calculations finish; an empty output folder during computation is normal.

| Option | Default | Purpose |
|---|---|---|
| `--threads` | 1 | Number of OpenMP galaxy workers |
| `--tables` | 1000 | Points per lookup-table axis; synchrotron is one-dimensional |
| `--cosmic-rays` | 1000 | Cosmic-ray energy samples |
| `--photons` | 500 | Photon energy samples |
| `--cells` | 500 | Electron solver cells |
| `--cache` | None | Directory for saving and reusing lookup tables |
| `--ebl` | Bundled Franceschini table | Alternative EBL input file |
| `--quiet` | Off | Suppress stage and table progress messages |
| `--allow-unconverged` | Off | Diagnostic runs only: permit unmet native quadrature tolerances, recording them in metadata |

Use `python -m congruents --help` for the command-line options.

Lookup-table generation is serial and can take substantial time on a cold run.
`--threads` speeds up the galaxy calculations, not table generation. Useful
parallelism is limited by both the number of catalogue cases and available CPUs.

Progress messages identify the table being loaded or generated, its completion
time, the galaxy-calculation stage and output writing. They do not provide an
ETA or percentage within a table or galaxy solve.

### Repeat a run using the cache

Keep the cache directory and choose a new output directory:

```sh
PYTHONPATH=src python -u -m congruents input/cat_nt.txt output/production-repeat \
  --threads 4 --cache data/production-cache
```

Compatible tables are reused; galaxy calculations run again. Cache filenames
describe their contents, for example `ic_emission_CMB_T2.7K_1000x1000.cgt`.
Metadata and checksums inside each file determine whether it can be reused.
Source, dependency or configuration changes can trigger regeneration, including
some code reorganisations that do not change the physics. An incompatible table
is replaced only after the new one is complete. Corrupt matching caches raise
an error.

Use separate cache directories when retaining different physics configurations.
To force a cold run, use a new, empty cache directory. Deleting results alone
does not clear the lookup-table cache.


## Input catalogue

`input/cat_nt.txt` contains 11 model cases, including adjusted cases. Its format is:

```text
n_gal
<number of rows>
z Mstar__Msol Re__kpc SFR__Msolyrm1
<redshift> <stellar mass> <effective radius> <star-formation rate>
...
```

Mass is in solar masses, radius in kpc and star-formation rate in solar masses
per year. Values must be finite; mass, radius and SFR must be positive.
Observer-frame runs require positive redshift, with a supported maximum of 20.

**Preserve the bundled catalogue order.** The eleventh row (index 10) receives
an additional multiplier of 0.1 on its proton-calorimetry fraction. This is not selected by galaxy name;
review this prescription before using a different catalogue. Temperature-table
construction also requires a valid multi-plane CMB/FIR range, so some small or
single-galaxy catalogues are rejected.

## Plot the results

After a completed production run, run the following from the repository root:

```sh
python plot_figure9.py output/production
```

Pass the completed run directory directly; no spectra are copied or recalculated.
If omitted, the directory defaults to `output`. For headless use, set
`MPLBACKEND=Agg` and add `--no-show`.

The script displays the plot and saves these files inside the selected run:

- `mnras2023_figure9_model.png`
- `mnras2023_figure9_model.pdf`

Plotting does not rerun the model. The panel selection assumes the bundled
catalogue order and shows six cases in the Figure 9 arrangement. It does not
overlay digitised paper curves. On a headless machine, prefix the plotting
command with `MPLBACKEND=Agg` to save files without opening a window.

## Output files

| Files | Contents and units |
|---|---|
| `E_gam.txt`, `T_CR.txt`, `E_e.txt` | Photon energy, proton kinetic energy and electron total energy axes, GeV |
| `spec_*.txt` | Observer-frame E²-weighted component spectra, GeV cm⁻² s⁻¹ |
| `spec_total_photons.txt` | Sum of photon components, excluding neutrinos and the full-calorimetry comparison |
| `source.npz` | Full-precision source arrays, named by component; electron energies are total energies in GeV |
| `tau_gg.txt`, `tau_EBL.txt`, `tau_ff.txt` | Dimensionless internal, EBL and free-free optical depths |
| `CR_specs.txt`, `CR_specs_inj.txt` | Particle populations (GeV⁻¹) and injection rates (GeV⁻¹ s⁻¹) |
| `gal_data.txt` | Columns: h (pc), nH (cm⁻³), B (G), gas dispersion (km/s), area (pc²), gas surface density (M☉/pc²), SFR surface density (M☉/yr/pc²), stellar surface density (M☉/pc²), dust T (K) |
| `Urad_Ub.txt` | Total radiation, magnetic, CMB, FIR, 3000 K, 4000 K, 7500 K and UV energy densities, eV/cm³ |
| `tau_loss/` | Process-specific loss times, seconds |
| `E_loss_nucrit.txt` | Loss rates at the 1.49-GHz critical electron energy, GeV/s |
| `E_loss_leptons.txt` | Injected powers (GeV/s) followed by dimensionless emission/escape fractions |
| `L_radio.txt`, `L_gamma.txt` | Radio components at 1.49 GHz (W/Hz), and 0.1–100 GeV gamma-ray luminosity (GeV/s) |
| `distmod.txt` | Distance conversion factor, cm⁻²; not a magnitude distance modulus |
| `metadata.json` | Catalogue, grids, settings, provenance, elapsed time and cache counts |

In component filenames, `1`/`2` identify primary/secondary electrons and
`z1`/`z2` identify disc/halo. Most spectral arrays contain one row per catalogue
case. Proton rows use `T_CR.txt`; electron rows use `E_e.txt`, not `T_CR.txt`.
`CR_specs.txt` contains five rows per case (proton, primary/secondary disc,
primary/secondary halo); `CR_specs_inj.txt` contains four electron-injection rows
per case. Text output uses six-decimal scientific notation. Infinite loss times
can represent zero tabulated loss rates.

Source and observer metadata include a recursive Python source manifest, including
`physics/`, with relative file paths, per-file SHA256 hashes and an aggregate
digest. The native library has its own hash. Existing output metadata is not
rewritten when the code changes.

## Model settings to check before scientific use

- Cosmology is flat, with H0 = 70 km/s/Mpc, Ωm = 0.3 and ΩΛ = 0.7.
- The default EBL dataset is Franceschini. Alternative bundled tables can be
  selected with `--ebl`.
- Photon attenuation uses internal optical depths at emitted energy and EBL
  optical depths at observed energy. Neutrinos receive neither photon opacity.
- Halo bremsstrahlung emission is not included. Halo loss-time tables and
  the halo solver use plasma losses; the disc uses neutral ionisation losses.
- Default IC-emission, bremsstrahlung and synchrotron tables use
  six-decimal scientific-notation rounding; IC energy-transfer tables do not.
- Adaptive integrals have finite evaluation budgets. Production observer runs
  reject unmet native spectral-worker integration tolerances by default.
  The lower-level `solve()` API defaults to warning and records the audit;
  pass `require_converged=True` to reject such results. This audit does not
  certify quadrature accuracy independently. Transport normalization and Python
  table generation now reject estimates that miss their requested tolerance.
- `E_loss_leptons.txt` is a legacy diagnostic, **not a closed energy budget**:
  it includes attenuated emission, placeholder columns and disc-to-halo transfer.
  Its fractions must not be summed as a conservation test.

### Corrected physics revision (October 2026)

This revision changes predictions relative to the improvised legacy code:

- Synchrotron cooling uses the isotropic `4/3` Thomson coefficient (Roth et al.
  2023, equation 26). Emission uses the isotropic pitch-angle kernel `R(x)`
  and the corresponding critical frequency and normalization, with a beta-squared
  factor to match the integrated cooling power. The spectral shape remains the
  relativistic synchrotron approximation, not a cyclotron model.
- CMB/FIR IC tables use ordinary linear temperature interpolation, including
  exact endpoints, instead of reversed weights.
- Momentum-power-law injection uses `dp/dT = (T+m)/p`; its energy normalization
  is recalculated using the same spectrum (equation 14).
- Young-stellar luminosity follows both branches and the scaled-SFR threshold
  of equation 13. Radiation diagnostic components now sum to their reported total.
- Halo loss-time diagnostics use plasma losses, consistent with the solver.
- Observer conversion evaluates source spectra and internal opacity at
  (1+z)E_obs, while EBL opacity is evaluated at E_obs. Photon attenuation is no
  longer applied to neutrinos. The existing luminosity-distance convention is
  retained and independently tested.
- Electron injection and diffusion use kinetic energy derived from the electron
  total-energy axis, rather than the separate proton kinetic-energy grid.
- Bremsstrahlung photons cannot consume electron rest energy. Low-frequency
  Planck/free-free factors and high-energy pair cross sections use stable forms.
- The continuous-cooling lower-face slope follows the linear reconstruction in
  Roth et al. (2023), Appendix B4. This corrects the inconsistent plus sign in
  the printed B6/legacy stencil. The upper boundary has zero incoming flux.
  Independent analytic cooling/escape solutions test this substantial change;
  it should also receive scientific review.
- Nonlocal electron transition moments integrate the piecewise-linear lookup
  table exactly in transferred energy, avoiding missed narrow IC peaks.
  Energy-weighted spectral integrals likewise split exactly at interpolation knots.
- IC and bremsstrahlung cooling-time diagnostics now integrate transferred
  energy down to a final electron energy of m_e c^2, including losses below the
  tracked grid. They no longer report an artificial infinity at the grid's lower
  edge. These exact table moments replace the poorly resolved final-energy
  quadrature; a constant-rate analytic test checks both losses and density scaling.
- Solver failure handling now audits integration error estimates and LU residuals,
  records negative populations before clipping, and rejects materially negative
  solutions. Zero injection and zero high-energy tails are handled without
  log(0) extrapolation. These are numerical safeguards, not adjusted physics.
- Electron-axis export (`E_e.txt`), numerical audit metadata and the native ABI
  were updated. Independent development checks exercise analytic processes, threading and
  failure recovery; these are archived separately. Existing outputs retain
  their original provenance and are not regenerated by a source update.

For radio observations, compare like components: a measured total radio flux
includes thermal free-free emission as well as synchrotron. The Figure 9 output
is not a dust-continuum prediction; millimetre dust measurements cannot validate
its synchrotron component directly.

Rebuild the native library before running. Interface revision 8 rejects old
libraries because table conventions, electron axes and numerical diagnostics changed. Source fingerprints
invalidate old table caches; do not copy old tables over newly generated ones.
Existing output files are not updated automatically. Save corrected runs in a
separate output directory. Metadata records the physics revision and source hashes.
Exact agreement with old external reference outputs is no longer expected;
agreement must be evaluated separately using matching grids and conventions.

Injection remains normalized to the stated energy budget, using the adopted
SN rate of 0.0132168 per solar mass formed, 10^50 erg per SN in protons and
2 x 10^49 erg in electrons. The last value follows the paper's 2% statement,
not its conflicting parenthetical 10^49 erg. For unit SFR, independent kinetic-
energy integration gives coefficients 5.36423e42 (protons) and 2.37584e41
(electrons) multiplying p^-2.2 dp/dT exp(-T/Tcut), with p numerically in GeV/c.
The proton coefficient is 13.76% below the paper's quoted 6.22e42; this is not
removed by a curve-fitting multiplier. These checks validate the integral for
the adopted SN rate, not the IMF derivation of that rate.

The hadronic implementation uses Np = Qp fcal tau_loss, where
tau_loss^-1 = (1.4/1.17) nH sigma0 eta c, sigma0 = 40 mb and eta = 0.5.
Combining this population with nH c beta dsigma/dE gives a kernel coefficient
beta / [(1.4/1.17) eta sigma0], or 1.67143 beta times the literal coefficient
in equation 39, at fixed injection and differential cross section. A loss-time
description is dimensionally consistent and distinguishes energy-loss time
from collision time, but the composition/nuclear-enhancement bookkeeping has
not been independently closed. This normalization is retained, not certified
or arbitrarily rescaled. No blanket Chen & Totani normalization factor is applied.
Pair absorption retains the head-on approximation, and the Geant4 hadronic model
is extrapolated at high energies. EBL extrapolation and the existing 1-PeV cap
also require care outside the tabulated domain.

Native metadata records quadrature misses, LU backward errors, negative-cell
counts and their particle/total-energy fractions. A negative energy fraction
above 1e-10 rejects the calculation; smaller tails are clipped and reported.
This numerical safeguard is not a physical error tolerance. A componentwise
post-clipping residual can be large in negligible tail cells, so inspect it
together with the recorded energy fraction. Very coarse grids can legitimately
fail these checks and must not be used as scientific results.

The revision is **not yet scientifically certified**. Full intrinsic energy
closure, production-grid convergence in every axis, IMF/hadronic normalization
and independent observational comparisons remain required.
CANDELS selection and cosmic-background weighting belong to the later survey
implementation, not this individual-galaxy model.

## Source layout

- `src/congruents/physics/`: radiation fields and serial process kernels.
- `src/congruents/`: catalogue handling, table preparation, Python API and outputs.
- `csrc/physics/`: per-process kernels used by native galaxy calculations.
- `csrc/solver.c`: galaxy calculations and OpenMP dispatch.
- `input/`: catalogue and EBL datasets.
- `plot_figure9.py`: plotting from completed outputs.
- `build/`, `data/`, `output/`: generated library, caches and results; ignored by Git.

## License

GPL-2.0. See [LICENSE](LICENSE) and the source-file attribution notices.
