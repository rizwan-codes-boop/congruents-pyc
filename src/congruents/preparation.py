"""Python-owned catalogue properties, radiation fields and lookup tables."""
from array import array
import bisect
from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import struct
import sys
import tempfile
import threading
import time
import math
from . import tables, constants
from .physics import radiation
import scipy
import numpy as np

from .galaxy import properties as galaxy_properties, PROPERTY_NAMES
from .inputs import Galaxy, Grid
from .model import Context

FIELDS = ("3000", "4000", "7500", "UV", "CMB", "FIR")
KINDS = ("emission", "gamma", "bs", "sy")
RADIATION_NAMES = ("CMB", "FIR", "3000", "4000", "7500", "UV", "total")

@dataclass(frozen=True)
class TableData:
    """Serializable axes and row-major values; constructor inputs may also be arrays."""
    x: tuple
    y: tuple
    values: tuple  # Snapshot order is flattened [iy*nx+ix]; Table also accepts [ny,nx].

class Table:
    """Read-only Python/NumPy tables; serial interpolation stays in Python."""
    def __init__(self, data, metadata):
        """Copy table axes and values into validated, immutable owned NumPy arrays."""
        self.metadata = dict(metadata)
        self._lock = threading.RLock()
        self._closed = False
        self._x = np.array(data.x, dtype=np.float64, copy=True)
        self._y = np.array(data.y, dtype=np.float64, copy=True)
        if self._x.ndim != 1 or self._y.ndim != 1:
            raise ValueError("Table axes must be one-dimensional")
        nx, ny = self._x.size, self._y.size
        if not 2 <= nx <= 4096 or not 1 <= ny <= 4096 or nx*ny > 4194304:
            raise ValueError("Invalid table shape")
        for a in (self._x,self._y):
            if not np.isfinite(a).all() or (a<=0).any() or (np.diff(a)<=0).any():
                raise ValueError("Table axes must be finite, positive and increasing")
        values = np.asarray(data.values, dtype=np.float64)
        if values.size != nx*ny or not np.isfinite(values).all() or (values<0).any():
            raise ValueError("Invalid table values")
        self._values = np.empty((ny,nx),dtype=np.float64)
        self._values[:] = values.reshape(ny,nx)
        for a in (self._x,self._y,self._values):
            a.flags.writeable = False

    def _ensure_open(self):
        """Reject operations after this object has been closed."""
        if self._closed:
            raise RuntimeError("Table is closed")

    def _view(self, a):
        """Return a read-only array view while preserving the table lifetime checks."""
        with self._lock:
            self._ensure_open()
            return a.view()

    @property
    def x(self):
        """Expose the read-only first table axis."""
        return self._view(self._x)

    @property
    def y(self):
        """Expose the read-only second table axis, including the synchrotron singleton."""
        return self._view(self._y)

    @property
    def values(self):
        """Expose the read-only row-major kernel values."""
        return self._view(self._values)

    def snapshot(self):
        """Copy table content into immutable tuples suitable for binary-cache export."""
        with self._lock:
            self._ensure_open()
            return TableData(tuple(self._x),tuple(self._y),tuple(self._values.ravel()))

    def evaluate(self, x, y=None):
        """Return interpolated values for paired coordinates as a Python list.

        Interpolation is linear in the stored coordinates, not their logarithms.
        A singleton y axis selects 1-D interpolation; otherwise x and y must be
        equal-length finite vectors, not mesh axes. Out-of-domain queries return
        zero. Coordinate/value units follow the table metadata.
        """
        with self._lock:
            self._ensure_open()
            x=np.asarray(tuple(x),dtype=np.float64)
            if x.ndim!=1 or not np.isfinite(x).all():
                raise ValueError("Coordinates must be finite one-dimensional arrays")
            if len(self._y)==1:
                return np.interp(x,self._x,self._values[0],left=0,right=0).tolist()
            if y is None:
                raise ValueError("A two-dimensional table needs electron coordinates")
            y=np.asarray(tuple(y),dtype=np.float64)
            if y.shape!=x.shape or not np.isfinite(y).all():
                raise ValueError("Coordinate lengths must match and be finite")
            valid=(x>=self._x[0])&(x<=self._x[-1])&(y>=self._y[0])&(y<=self._y[-1])
            out=np.zeros_like(x)
            xv,yv=x[valid],y[valid]
            ix=np.clip(np.searchsorted(self._x,xv,side="right")-1,0,len(self._x)-2)
            iy=np.clip(np.searchsorted(self._y,yv,side="right")-1,0,len(self._y)-2)
            u=(xv-self._x[ix])/(self._x[ix+1]-self._x[ix])
            v=(yv-self._y[iy])/(self._y[iy+1]-self._y[iy])
            z=self._values
            out[valid]=(1-v)*((1-u)*z[iy,ix]+u*z[iy,ix+1])+v*((1-u)*z[iy+1,ix]+u*z[iy+1,ix+1])
            return out.tolist()

    def close(self):
        """Invalidate this table handle without changing arrays retained by existing readers."""
        with self._lock:
            self._closed=True

    def __enter__(self):
        """Verify that the resource is open before entering its managed scope."""
        self._ensure_open()
        return self

    def __exit__(self,*args):
        """Close managed resources when leaving the scope, including on exceptions."""
        self.close()

def _save_cache(path, metadata, data):
    """Write little-endian float64 axes/payload with versioned metadata and SHA256.

    The temporary file shares the destination filesystem so os.replace is
    atomic. This format stores numerical data only; it does not unpickle code.
    """
    payload = array("d", data.x + data.y + data.values)
    if sys.byteorder != "little":
        payload.byteswap()
    raw = payload.tobytes()
    header = json.dumps({"metadata": metadata, "nx": len(data.x), "ny": len(data.y),
                         "sha256": hashlib.sha256(raw).hexdigest()}, sort_keys=True).encode()
    path.parent.mkdir(parents=True, exist_ok=True)
    # Atomic replacement; an interrupted writer cannot leave a partial cache.
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as stream:
        tmp = Path(stream.name)
        try:
            stream.write(b"CGT1" + struct.pack("<I", len(header)) + header + raw)
        except BaseException:
            tmp.unlink(missing_ok=True)
            raise
    try:
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)

class _CacheMismatch(ValueError):
    """A readable cache name exists but belongs to another configuration."""


def _read_cache(path, expected):
    """Validate schema, metadata, exact byte count and checksum before importing.

    Corrupt or mismatched existing files raise instead of being silently used.
    """
    with path.open("rb") as stream:
        prefix = stream.read(8)
        if len(prefix) != 8 or prefix[:4] != b"CGT1":
            raise ValueError("Invalid cache header")
        length = struct.unpack("<I", prefix[4:])[0]
        if length > 65536:
            raise ValueError("Oversized cache header")
        header = json.loads(stream.read(length))
        if header["metadata"] != expected:
            raise _CacheMismatch("Cache metadata mismatch")
        nx, ny = header["nx"], header["ny"]
        expected_ny = 1 if expected["kind"] == "sy" else expected["ny"]
        if nx != expected["nx"] or ny != expected_ny:
            raise ValueError("Cache shape mismatch")
        size = (nx + ny + nx * ny) * 8
        raw = stream.read(size + 1)
        if len(raw) != size or hashlib.sha256(raw).hexdigest() != header["sha256"]:
            raise ValueError("Cache payload length/checksum mismatch")
    values = array("d")
    values.frombytes(raw)
    if sys.byteorder != "little":
        values.byteswap()
    return TableData(tuple(values[:nx]), tuple(values[nx:nx+ny]), tuple(values[nx+ny:]))

def cache_filename(metadata):
    """Readable name only; internal metadata determines compatibility."""
    kind, field = metadata["kind"], metadata["field"]
    if kind not in KINDS or field not in FIELDS:
        raise ValueError("Unknown cache kind/field")
    nx, ny = metadata["nx"], metadata["ny"]
    if type(nx) is not int or type(ny) is not int or not (2 <= nx <= 4096 and 2 <= ny <= 4096):
        raise ValueError("Invalid cache dimensions")
    label = {"emission":"ic_emission", "gamma":"ic_energy_transfer",
             "bs":"bremsstrahlung", "sy":"synchrotron"}[kind]
    if kind in ("emission", "gamma"):
        label += "_" + field
        if field in ("CMB", "FIR"):
            temperature = float(metadata["temperature_k"])
            if not math.isfinite(temperature) or not 2 <= temperature <= 1000:
                raise ValueError("Invalid cache temperature")
            label += f"_T{temperature:.8g}K"
    shape = str(nx) if kind == "sy" else f"{nx}x{ny}"
    return f"{label}_{shape}.cgt"


def _rename_without_overwrite(source, target):
    """Move a cache name within one filesystem, refusing an existing destination."""
    os.link(source, target)
    source.unlink()


def rename_cache_files(directory):
    """Remove old filename hashes, retaining bytes and original fingerprints.

    Validate every candidate and destination before moving any. A historical
    fingerprint is not upgraded: incompatible tables stay incompatible.
    No native library is loaded and no table is generated by this operation.
    """
    plan = []
    targets = set()
    for path in sorted(Path(directory).glob("*.cgt")):
        old_hash = path.stem.rsplit("_", 1)[-1]
        if len(old_hash) != 64 or any(c not in "0123456789abcdef" for c in old_hash):
            continue
        if path.is_symlink():
            raise ValueError(f"Refusing cache symlink: {path}")
        with path.open("rb") as stream:
            prefix = stream.read(8)
            if len(prefix) != 8 or prefix[:4] != b"CGT1":
                raise ValueError(f"Invalid cache header: {path}")
            length = struct.unpack("<I", prefix[4:])[0]
            if length > 65536:
                raise ValueError("Oversized cache header")
            meta = json.loads(stream.read(length))["metadata"]
        target = path.with_name(cache_filename(meta))
        key = hashlib.sha256(json.dumps(meta, sort_keys=True, allow_nan=False).encode()).hexdigest()
        if old_hash != key:
            raise ValueError(f"Cache filename/metadata mismatch: {path}")
        if meta.get("schema") != 2:
            raise ValueError("Unsupported cache schema")
        _read_cache(path, meta)
        if target.exists() or target.is_symlink() or target in targets:
            raise FileExistsError(target)
        targets.add(target)
        plan.append((path, target))
    for source, target in plan:
        _rename_without_overwrite(source, target)
    return plan


class Preparation:
    """Own reusable per-field tables for a validated catalogue.

    Defaults reproduce the production grid sizes; use Grid(8,8) ONLY for smoke
    tests. Cache files are new versioned binary artifacts, never legacy text
    files. Galaxy properties use the native OpenMP stage; tables are generated
    serially in Python after those properties establish catalogue-wide bounds.
    """
    def __init__(self, galaxies, grid=None, threads=1, cache=None, library=None, progress=None):
        """Validate the catalogue, derive native galaxy properties and establish cache identity."""
        if progress is not None and not callable(progress):
            raise TypeError("progress must be callable or None")
        self._progress = progress
        self._started = time.perf_counter()
        self.galaxies = tuple(galaxies)
        if not self.galaxies or len(self.galaxies) > 100000:
            raise ValueError("Require 1..100000 galaxies")
        if not all(isinstance(g, Galaxy) for g in self.galaxies):
            raise TypeError("Expected Galaxy objects")
        self.grid = grid if grid is not None else Grid()
        if not isinstance(self.grid, Grid):
            raise TypeError("grid must be Grid")
        self._lock = threading.RLock()
        self._context = Context(threads, library)
        self._tables = {}
        self.cache = None if cache is None else Path(cache)
        self.cache_hits = self.cache_misses = 0
        # The cache identity hashes literal source bytes plus NumPy/SciPy versions.
        # Even documentation edits invalidate keys; incompatible tables are rebuilt.
        source = b"".join(Path(__file__).with_name(name).read_bytes() for name in
                          ("tables.py","grids.py","constants.py","quadrature.py","preparation.py","galaxy.py"))
        source += b"".join(p.read_bytes() for p in sorted(
            Path(__file__).with_name("physics").glob("*.py")))
        from ._native import load
        native = Path(load(library)._name).read_bytes()
        self._preparation_hash = hashlib.sha256(source+native+
                          (np.__version__+scipy.__version__).encode()).hexdigest()
        self.report(f"Calculating properties for {len(self.galaxies)} catalogue cases")
        output = galaxy_properties([g.as_row() for g in self.galaxies], threads, library)
        self.properties = tuple(dict(zip(PROPERTY_NAMES, row)) for row in output)
        self.config, self.temperature_bounds = tables.bounds(self.galaxies,self.properties)
        self.report("Galaxy properties ready")

    def report(self, message):
        """Report Python-side milestones only; never called from native workers."""
        if self._progress is not None:
            self._progress(f"[{time.perf_counter()-self._started:.1f}s] {message}")

    def _ensure_open(self):
        """Reject operations after this object has been closed."""
        if self._context._closed:
            raise RuntimeError("Preparation is closed")

    def temperature_grid(self, field):
        """Return kelvin nodes for a catalogue-wide CMB or FIR table family."""
        with self._lock:
            self._ensure_open()
            return radiation.temperature_grid(field,self.temperature_bounds)

    def radiation(self,index,energies):
        """Return radiation diagnostics for one zero-based catalogue index.

        energies is a one-dimensional sequence of photon energies in GeV.
        fields_cm3_gev maps component names to number-density spectra in
        cm^-3 GeV^-1. urad_ub_ev_cm3 contains energy densities in eV/cm^3,
        ordered total radiation, magnetic, CMB, FIR, 3000 K, 4000 K, 7500 K, UV.
        """
        with self._lock:
            self._ensure_open()
            fields,urad=radiation.radiation(self.galaxies[index],self.properties[index],energies,self.config)
            return {"fields_cm3_gev":dict(zip(RADIATION_NAMES,fields)),"urad_ub_ev_cm3":urad}

    def table(self, kind, field="3000", temperature=0.):
        """Lazily return a shared read-only field table, using the optional binary cache.

        BS and SY do not depend on radiation field, so their field/temperature
        keys are canonicalized. The Preparation object owns the returned table.
        """
        with self._lock:
            self._ensure_open()
            if kind not in KINDS or field not in FIELDS:
                raise ValueError("Unknown table kind/field")
            if kind in ("bs", "sy"):
                field, temperature = "3000", 0.
            if field not in ("CMB", "FIR"):
                temperature = 0.
            meta = {"schema": 2, "implementation_sha256": self._preparation_hash,
                    "kind": kind, "field": field, "temperature_k": float(temperature),
                    "nx": self.grid.nx, "ny": self.grid.ny, "config": list(self.config),
                    "units": "dimensionless" if kind == "sy" else "mb/GeV" if kind == "bs" else "1/(s GeV)",
                    "x_axis": "x" if kind == "sy" else "DeltaE_GeV" if kind == "gamma" else "photon_GeV",
                    "y_axis": "unused" if kind == "sy" else "electron_total_GeV"}
            encoded = json.dumps(meta, sort_keys=True, allow_nan=False).encode()
            key = hashlib.sha256(encoded).hexdigest()
            if key in self._tables:
                return self._tables[key]
            name = cache_filename(meta)
            label = Path(name).stem
            path = None if self.cache is None else self.cache / name
            candidates = [] if self.cache is None else [path, self.cache / (key + ".cgt"),
                self.cache / (Path(name).stem + "_" + key + ".cgt")]
            existing = next((candidate for candidate in candidates if candidate.exists()), None)
            started = time.perf_counter()
            number = len(self._tables) + 1
            data = None
            if existing is not None:
                self.report(f"Table {number}: loading cached {label}")
                try:
                    data = _read_cache(existing, meta)
                except _CacheMismatch:
                    self.report(f"Table {number}: incompatible metadata for {label}; regenerating")
            generated = data is None
            if not generated:
                if existing != path:
                    _rename_without_overwrite(existing, path)
                table = Table(data, meta)
                self.cache_hits += 1
            else:
                self.report(f"Table {number}: generating {label} (serial Python; no compatible cache)")
                x,y,z=tables.generate(kind,field,temperature,self.grid.nx,self.grid.ny,self.config)
                table=Table(TableData(x,y,z),meta)
                self.cache_misses += 1
            try:
                if path is not None and generated:
                    _save_cache(path, meta, table.snapshot())
                    self.report(f"Saved {path.name}")
            except BaseException:
                table.close()
                raise
            self._tables[key] = table
            self.report(f"Table {number} ready: {label} ({time.perf_counter()-started:.1f}s)")
            return table

    def prepare_tables(self):
        """Prepare all 14 logical table families; CMB/FIR contain multiple planes."""
        with self._lock:
            for kind in ("emission", "gamma"):
                for field in FIELDS[:4]:
                    self.table(kind, field)
                for field in ("CMB", "FIR"):
                    for temperature in self.temperature_grid(field):
                        self.table(kind, field, temperature)
            self.table("bs")
            self.table("sy")
            return len(self._tables)

    def close(self):
        """Invalidate owned table handles and configuration; safe to call repeatedly."""
        with self._lock:
            for table in self._tables.values():
                table.close()
            self._tables.clear()
            self._context.close()

    def __enter__(self):
        """Verify that the resource is open before entering its managed scope."""
        self._ensure_open()
        return self

    def __exit__(self, *args):
        """Close managed resources when leaving the scope, including on exceptions."""
        self.close()

    def __del__(self):
        """Release owned tables if object destruction occurs before explicit close."""
        if hasattr(self, "_tables"):
            self.close()
