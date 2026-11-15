"""Runtime oracle: Verdict, the library ``vtkMeshQuality``/meshmagick call.

``docs/mesh_quality.md`` §5 used to record the meshmagick arm of the acceptance criterion as
*not run*: meshmagick's ``print_quality()`` delegates to ``vtk.vtkMeshQuality``, no VTK wheel
exists for loongarch64, and re-deriving Verdict's algebra by hand would be the third
implementation of the same formulas -- the opposite of a cross-check.  Verdict is the code
behind ``vtkMeshQuality``, so compiling it (``tools/verdict_oracle/build.py``) gives an
oracle with a *different* producer than this tree and with the same definitions meshmagick
would have used.

This module is the Python side of that shared library:

* :func:`available` / :func:`library_path` -- is an oracle built, and where;
* :func:`provenance` -- which VTK reference and compiler built it (``provenance.json``);
* :func:`metric_names` / :func:`metrics` -- Verdict's numbers, per element.

The oracle is optional, exactly like ``trimesh`` and ``gmsh`` in
:mod:`tensorcad.meshquality.reference`: nothing in the library imports it at module load,
``available()`` is False on a host without a compiler, and the callers record that as
``not run`` rather than as a pass.

The library is found through, in order: the ``TENSORCAD_VERDICT_ORACLE`` environment
variable, ``tools/verdict_oracle/build/libverdict_oracle.so`` in this checkout, and
``~/.cache/tensorcad-verdict/``.
"""
from __future__ import annotations

import ctypes
import json
import os
from pathlib import Path
from typing import Iterable, Sequence

__all__ = [
    "KINDS",
    "LIBRARY_ENV",
    "available",
    "info",
    "library_path",
    "metric_names",
    "metrics",
    "provenance",
]

#: Cell kinds the wrapper exposes, in the order the build assigns them.
KINDS: tuple[str, ...] = ("tri", "quad", "tet", "hex")

LIBRARY_ENV = "TENSORCAD_VERDICT_ORACLE"
LIBRARY_NAME = "libverdict_oracle.so"

_REPO_ROOT = Path(__file__).resolve().parents[2]
_CANDIDATES = (
    _REPO_ROOT / "tools" / "verdict_oracle" / "build" / LIBRARY_NAME,
    Path.home() / ".cache" / "tensorcad-verdict" / LIBRARY_NAME,
)

_cache: dict[str, object] = {}


def library_path() -> Path | None:
    """Path of the compiled oracle, or ``None`` when this host has not built one."""
    override = os.environ.get(LIBRARY_ENV)
    candidates: Iterable[Path] = (Path(override),) if override else _CANDIDATES
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    return None


def _library() -> ctypes.CDLL | None:
    cached = _cache.get("library")
    if cached is not None or _cache.get("missing"):
        return cached  # type: ignore[return-value]
    path = library_path()
    if path is None:
        _cache["missing"] = True
        return None
    lib = ctypes.CDLL(str(path))
    lib.verdict_oracle_info.restype = ctypes.c_char_p
    lib.verdict_oracle_info.argtypes = []
    lib.verdict_oracle_kind_name.restype = ctypes.c_char_p
    lib.verdict_oracle_kind_name.argtypes = [ctypes.c_int]
    lib.verdict_oracle_kind_count.restype = ctypes.c_int
    lib.verdict_oracle_kind_count.argtypes = []
    lib.verdict_oracle_nodes_per_element.restype = ctypes.c_int
    lib.verdict_oracle_nodes_per_element.argtypes = [ctypes.c_int]
    lib.verdict_oracle_metric_count.restype = ctypes.c_int
    lib.verdict_oracle_metric_count.argtypes = [ctypes.c_int]
    lib.verdict_oracle_metric_name.restype = ctypes.c_char_p
    lib.verdict_oracle_metric_name.argtypes = [ctypes.c_int, ctypes.c_int]
    lib.verdict_oracle_evaluate.restype = ctypes.c_int
    lib.verdict_oracle_evaluate.argtypes = [
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
    ]
    # Reject a library that is not this wrapper (a stale .so must not be mistaken for it).
    info = lib.verdict_oracle_info().decode()
    if not info.startswith("verdict-oracle/"):
        raise RuntimeError(f"{path} is not the Verdict oracle (info={info!r})")
    _cache["library"] = lib
    _cache["path"] = path
    return lib


def available() -> bool:
    """True when a Verdict oracle is built on this host."""
    return _library() is not None


def info() -> str | None:
    """The library's self-identification string, or ``None`` when unavailable."""
    lib = _library()
    return lib.verdict_oracle_info().decode() if lib else None


def provenance() -> dict[str, object]:
    """The ``provenance.json`` written by the build, plus how the library was found."""
    path = library_path()
    record: dict[str, object] = {"library_path": str(path) if path else None}
    if path is None:
        record["built"] = False
        return record
    record["built"] = True
    sidecar = path.parent / "provenance.json"
    if sidecar.is_file():
        record.update(json.loads(sidecar.read_text()))
    return record


def _kind_index(kind: str) -> int:
    try:
        return KINDS.index(kind)
    except ValueError:
        raise ValueError(f"unknown cell kind {kind!r}; the oracle knows {', '.join(KINDS)}") from None


def metric_names(kind: str) -> tuple[str, ...]:
    """Verdict metric names the oracle exposes for ``kind``."""
    lib = _library()
    if lib is None:
        raise RuntimeError("no Verdict oracle built; run tools/verdict_oracle/build.py")
    index = _kind_index(kind)
    count = lib.verdict_oracle_metric_count(index)
    return tuple(
        lib.verdict_oracle_metric_name(index, i).decode() for i in range(count)
    )


def metrics(
    kind: str, points: Sequence[Sequence[float]], cells: Sequence[Sequence[int]]
) -> dict[str, list[float]]:
    """Verdict's per-element metrics for ``cells`` of ``kind``.

    Returns one list per metric name (see :func:`metric_names`), in element order.  Point
    and cell ordering is the caller's: the oracle is fed the nodes in the element's own
    order, exactly as :mod:`tensorcad.meshquality.metrics` defines it.
    """
    lib = _library()
    if lib is None:
        raise RuntimeError("no Verdict oracle built; run tools/verdict_oracle/build.py")
    index = _kind_index(kind)
    nodes = lib.verdict_oracle_nodes_per_element(index)
    if nodes < 0:
        raise ValueError(f"the oracle does not know the node count of {kind!r}")
    names = metric_names(kind)
    out: dict[str, list[float]] = {name: [] for name in names}
    coords = (ctypes.c_double * (nodes * 3))()
    values = (ctypes.c_double * len(names))()
    for cell in cells:
        if len(cell) != nodes:
            raise ValueError(f"{kind} element has {len(cell)} nodes, expected {nodes}")
        for slot, node in enumerate(cell):
            point = points[node]
            coords[slot * 3 + 0] = float(point[0])
            coords[slot * 3 + 1] = float(point[1])
            coords[slot * 3 + 2] = float(point[2])
        if lib.verdict_oracle_evaluate(index, coords, values) != 0:
            raise RuntimeError(f"Verdict oracle refused a {kind} element")
        for name in names:
            out[name].append(0.0)  # placeholder replaced below
        for i, name in enumerate(names):
            out[name][-1] = values[i]
    return out