"""Reproducible source-file manifests for source and observer outputs."""
import hashlib
import json
from pathlib import Path


def python_source_provenance(root=None):
    """Return a location-independent manifest and its versioned aggregate digest.

    Only Python source is included: bytecode, caches and outputs are excluded.
    The optional root supports testing against isolated source trees. Call once
    at the start of a solve and reuse that snapshot for every exported format.
    """
    root = Path(root) if root is not None else Path(__file__).parent
    files = {path.relative_to(root).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
             for path in sorted(root.rglob("*.py")) if path.is_file()}
    if not files:
        raise ValueError("No Python source files found for provenance")
    encoded = json.dumps(files, sort_keys=True, separators=(",", ":")).encode()
    return {"python_source_sha256": hashlib.sha256(encoded).hexdigest(),
            "python_source_files_sha256": files,
            "python_source_hash_scheme": "relative-path-sha256-manifest-v1"}
