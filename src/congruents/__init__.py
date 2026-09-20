"""Public API for Python serial workflow and complete galaxy-parallel native physics.

Importing these names does not load the shared library. solve returns source
components; run adds diagnostics, attenuation and observer-frame outputs.
"""
from .model import Context
from .inputs import Galaxy, Grid, load_catalogue, write_catalogue
from .preparation import Preparation
from .solver import SolverGrid, SolverResult, solve
from .pipeline import RunResult, run
__all__ = ["Context", "Galaxy", "Grid", "load_catalogue", "write_catalogue", "Preparation",
           "SolverGrid", "SolverResult", "solve", "RunResult", "run"]
