# CONGRUENTS

CONGRUENTS calculates cosmic-ray transport and multiwavelength emission from
galaxies. It produces source spectra, observer-frame spectra, optical depths
and energy-loss diagnostics. The model follows the improvised legacy code.

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

The Python installation does not compile the native library. Build it using
one of the configurations below, and rebuild after changing C source files.

### Using the existing local dependency installation

If the sibling `CONGRUENTS-c` folder contains the configured compiler and
libraries:

```sh
make DEPENDENCY_ROOT=../CONGRUENTS-c solver
```

This uses dependencies from that folder, not its model code or outputs.

### Building on another machine

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
The dependency build recipe is also available in
[the CI configuration](.github/workflows/interface.yml).

## Run the model

### Quick installation check

This small-grid run checks that the program works. Do not use its spectra for
scientific interpretation.

```sh
PYTHONPATH=src python -u -m congruents input/cat_nt.txt output/smoke \
  --threads 4 --tables 16 --cosmic-rays 16 --photons 8 --cells 16
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
a fixed proton-calorimetry factor of 0.1. This is not selected by galaxy name;
review this prescription before using a different catalogue. Temperature-table
construction also requires a valid multi-plane CMB/FIR range, so some small or
single-galaxy catalogues are rejected.

## Plot the results

After a completed production run, run the following from the repository root:

```sh
repo_dir="$PWD"
plot_dir=$(mktemp -d)
ln -s "$repo_dir/output/production" "$plot_dir/output"
(
  cd "$plot_dir"
  python "$repo_dir/plot_figure9.py"
)
```

The plotting script expects an `output` directory relative to its working
directory. The temporary link selects the run without copying or moving data.
For another run, replace `output/production` in the link command.

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
| `E_gam.txt`, `T_CR.txt` | Photon and cosmic-ray kinetic-energy axes, GeV |
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
case. `CR_specs.txt` contains five rows per case (proton, primary/secondary disc,
primary/secondary halo); `CR_specs_inj.txt` contains four electron-injection rows
per case. Text output uses six-decimal scientific notation. Infinite loss times
can represent zero tabulated loss rates.

Source and observer metadata include a recursive Python source manifest, including
`physics/`, with relative file paths, per-file SHA256 hashes and an aggregate
digest. The native library has its own hash. Existing output metadata is not
rewritten when the code changes.

## Tests and reference comparison

Tests are optional and do not run as part of a normal model calculation.
With the existing local dependency installation:

```sh
make DEPENDENCY_ROOT=../CONGRUENTS-c test PYTHON=python
```

For another dependency installation, pass the same compiler, `CPPFLAGS` and
`LDLIBS` used to build the library. The tests cover table caching, numerical
helpers, output formats, error handling and single-/multi-thread agreement.

To compare with independently generated improvised legacy code results:

```sh
python tests/compare_outputs.py /path/to/reference/output output/production
```

Match catalogue order, grids, physics, EBL choice and table precision. The
comparison checks 39 files with relative tolerance `3e-6`, absolute tolerance
`1e-280` and exact reference zeros. It requires external reference outputs;
these are not needed to run the model or the self-contained tests.
Passing tests is not proof of numerical convergence or agreement with observations.
Production-scale Linux/HPC validation is not established.

## Model settings to check before scientific use

- Cosmology is flat, with H0 = 70 km/s/Mpc, Ωm = 0.3 and ΩΛ = 0.7.
- The default EBL dataset is Franceschini. Alternative bundled tables can be
  selected with `--ebl`.
- Observer processing currently applies attenuation to all exported components,
  including neutrinos. Do not interpret this as a physical neutrino-attenuation model.
- Halo bremsstrahlung emission is not included. Halo loss-time tables use
  ionisation losses, while the halo solver and critical-energy diagnostics use
  plasma losses.
- Temperature interpolation retains reversed weighting from the improvised
  legacy code. Default IC-emission, bremsstrahlung and synchrotron tables use
  six-decimal scientific-notation rounding; IC energy-transfer tables do not.
- Adaptive integrals have finite evaluation budgets. A returned spectrum does
  not guarantee that every integral reached its requested tolerance.

## Source layout

- `src/congruents/physics/`: radiation fields and serial process kernels.
- `src/congruents/`: catalogue handling, table preparation, Python API and outputs.
- `csrc/physics/`: per-process kernels used by native galaxy calculations.
- `csrc/solver.c`: galaxy calculations and OpenMP dispatch.
- `input/`: catalogue and EBL datasets.
- `tests/`: self-contained checks and optional reference comparison.
- `plot_figure9.py`: plotting from completed outputs.
- `build/`, `data/`, `output/`: generated library, caches and results; ignored by Git.

## License

GPL-2.0. See [LICENSE](LICENSE) and the source-file attribution notices.
