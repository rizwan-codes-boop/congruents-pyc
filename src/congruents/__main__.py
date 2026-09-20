"""Run a catalogue without changing its scientific parameters or row order."""
import argparse
import time
from pathlib import Path
from . import Preparation, Grid, SolverGrid, load_catalogue, run


def main():
    """Validate CLI configuration, run one ordered catalogue and export a new run.

    Default grids are production-sized; --threads controls the native galaxy workers,
    not serial table generation. No test or external reference comparison is run here.
    """
    parser = argparse.ArgumentParser(description="CONGRUENTS source and observer spectra")
    parser.add_argument("catalogue",type=Path)
    parser.add_argument("output",type=Path,help="New output directory; existing paths are refused")
    parser.add_argument("--threads",type=int,default=1)
    parser.add_argument("--tables",type=int,default=1000)
    parser.add_argument("--cosmic-rays",type=int,default=1000)
    parser.add_argument("--photons",type=int,default=500)
    parser.add_argument("--cells",type=int,default=500)
    parser.add_argument("--cache",type=Path)
    parser.add_argument("--ebl",type=Path)
    parser.add_argument("--quiet",action="store_true",help="Suppress stage and table progress messages")
    args = parser.parse_args()
    # Reject overwrite before spending time generating tables or spectra.
    if args.output.exists():
        parser.error("Output already exists; choose a new directory")
    grid = SolverGrid(args.cosmic_rays,args.photons,args.cells)
    start = time.perf_counter()
    with Preparation(load_catalogue(args.catalogue),Grid(args.tables,args.tables),
                     threads=args.threads,cache=args.cache,
                     progress=None if args.quiet else lambda message: print(message,flush=True)) as preparation:
        print(f"Running {len(preparation.galaxies)} galaxies; tables={args.tables}, cells={args.cells}",flush=True)
        result = run(preparation,grid,args.threads,ebl_path=args.ebl)
        result.metadata.update(wall_seconds=time.perf_counter()-start,
                               cache_hits=preparation.cache_hits,cache_misses=preparation.cache_misses)
        preparation.report(f"Writing outputs to {args.output}")
        result.export(args.output)
        preparation.report(f"Cache summary: {preparation.cache_hits} loaded, {preparation.cache_misses} generated")
    print(f"Saved {args.output.resolve()} ({time.perf_counter()-start:.2f} seconds)")


if __name__ == "__main__":
    main()
