"""Compare 39 outputs against external improvised legacy code results.

Usage: python tests/compare_outputs.py C_OUTPUT PYTHON_OUTPUT
Exit nonzero on missing files, shape/NaN/zero/infinity mismatches or tolerance
failures. Missing reference directories report not_run (exit 2), never passed.
Reference-only debugging exports and Python-only provenance are ignored.
"""
import argparse
import json
import hashlib
import warnings
from pathlib import Path
import numpy as np

FILES = ["T_CR.txt","E_gam.txt","fcal.txt","gal_data.txt","Urad_Ub.txt",
         "tau_ff.txt","E_loss_leptons.txt","E_loss_nucrit.txt","CR_specs.txt",
         "CR_specs_inj.txt","L_radio.txt","distmod.txt","L_gamma.txt",
         "spec_pi.txt","spec_pi_fcal1.txt","spec_nu.txt","spec_FF.txt"]
FILES += [f"spec_{p}_{population}_z{zone}.txt" for p in ("IC","BS","SY")
          for population in (1,2) for zone in (1,2) if p!="BS" or zone==1]
FILES += [f"tau_loss/tau_loss_z{zone}_{p}.txt" for zone in (1,2) for p in ("SY","BS","IC","DI","IO")]
FILES += [f"tau_loss/tau_loss_protons_{p}.txt" for p in ("PP","DI")]
VECTOR_FILES = {"T_CR.txt", "E_gam.txt", "distmod.txt", "L_gamma.txt"}


def compare(reference, actual, tolerance=3e-6):
    """Compare every required numerical file; never fit or rescale either run.

    Callers must match catalogue order, grids, physics and precision first.
    NaNs fail; equal signed infinities can represent infinite loss times.
    The relative-error summary excludes zeros/infinities, but pass/fail does not.
    """
    if not np.isfinite(tolerance) or tolerance < 0:
        raise ValueError("Tolerance must be finite and nonnegative")
    if Path(reference).resolve() == Path(actual).resolve():
        raise ValueError("Reference and actual directories must be different")
    report = {}
    for name in FILES:
        try:
            skip = 0 if name.startswith("CR_specs") else 1
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", UserWarning)
                a = np.loadtxt(Path(reference)/name,skiprows=skip,ndmin=2)
                b = np.loadtxt(Path(actual)/name,skiprows=skip,ndmin=2)
            # These outputs are vectors; reference writers use both row and
            # column layouts. Never flatten galaxy-by-energy matrices.
            if name in VECTOR_FILES and 1 in a.shape and 1 in b.shape:
                a, b = a.reshape(-1), b.reshape(-1)
            if not a.size or not b.size or a.shape!=b.shape or np.isnan(a).any() or np.isnan(b).any():
                raise ValueError("Empty data, shape mismatch or NaN")
            close = np.isclose(a,b,rtol=tolerance,atol=1e-280)
            # A reference zero must remain zero, even below the absolute floor.
            close &= (a!=0)|(b==0)
            finite = np.isfinite(a)&np.isfinite(b)&(a!=0)
            relative = np.abs((b[finite]-a[finite])/a[finite])
            report[name] = {"passed":bool(close.all()),"values":int(a.size),
                            "reference_sha256":hashlib.sha256((Path(reference)/name).read_bytes()).hexdigest(),
                            "actual_sha256":hashlib.sha256((Path(actual)/name).read_bytes()).hexdigest(),
                            "mismatches":int((~close).sum()),
                            "max_relative_error":float(relative.max()) if relative.size else 0.}
        except (OSError,ValueError) as exc:
            report[name] = {"passed":False,"error":str(exc)}
    return report


def main():
    """Parse reference/result paths and report comparison status without modifying outputs."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path, nargs="?")
    parser.add_argument("actual", type=Path, nargs="?")
    args = parser.parse_args()
    if args.reference is None or args.actual is None or not args.reference.is_dir():
        print(json.dumps({"status":"not_run", "reason":
            "Supply an existing external reference output directory and an interface output directory."}))
        return 2
    try:
        files = compare(args.reference, args.actual)
    except ValueError as exc:
        print(json.dumps({"status":"failed", "reason":str(exc)}))
        return 1
    passed = all(r["passed"] for r in files.values())
    print(json.dumps({"status":"passed" if passed else "failed",
        "reference":str(args.reference.resolve()), "actual":str(args.actual.resolve()),
        "rtol":3e-6, "atol":1e-280,
        "scope":"Numerical comparison only. Caller must ensure matching catalogue order, grids, physics and precision; reference run settings are not automatically verified.",
        "files":files}, indent=2))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
