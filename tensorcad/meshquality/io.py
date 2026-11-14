"""JSON exchange format for a mesh-quality check.

A quality gate in CI must be able to read a mesh without importing a mesh library, and
a solver must be able to hand one over without a VTK dependency. The spec is
deliberately tiny::

    {
      "points": [[x, y, z], ...],
      "elements": {
        "tri": [[i, j, k], ...],
        "hex": [[i0, ..., i7], ...]
      },
      "metadata": {"source": "...", "units": "m"}      # optional, echoed back
    }

Indexing is 0-based and the node order of each cell kind is the VTK/Gmsh order (the one
:mod:`tensorcad.meshquality.metrics` documents). Anything else -- a solver's native
ordering -- belongs in the caller's adapter, not in this file.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping, Sequence

__all__ = ["MeshSpec", "load_spec", "write_spec"]

from .metrics import ELEMENT_KINDS, NODES_PER_CELL


@dataclass(frozen=True)
class MeshSpec:
    """A mesh in the exchange format: points, per-kind connectivity, optional metadata."""

    points: tuple[tuple[float, float, float], ...]
    elements: dict[str, tuple[tuple[int, ...], ...]] = field(default_factory=dict)
    metadata: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {
            "points": [list(p) for p in self.points],
            "elements": {kind: [list(cell) for cell in cells] for kind, cells in self.elements.items()},
            "metadata": dict(self.metadata),
        }


def _as_cells(kind: str, raw: Sequence[Sequence[int]], count: int) -> tuple[tuple[int, ...], ...]:
    expected = NODES_PER_CELL[kind]
    cells: list[tuple[int, ...]] = []
    for position, cell in enumerate(raw):
        if len(cell) != expected:
            raise ValueError(f"{kind} cell {position} has {len(cell)} nodes, expected {expected}")
        indices = tuple(int(i) for i in cell)
        if any(i < 0 or i >= count for i in indices):
            raise ValueError(f"{kind} cell {position} indexes outside 0..{count - 1}: {indices}")
        cells.append(indices)
    return tuple(cells)


def load_spec(path: str | Path) -> MeshSpec:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    raw_points = payload.get("points")
    if not raw_points:
        raise ValueError("mesh spec needs a non-empty 'points' array")
    points = tuple((float(p[0]), float(p[1]), float(p[2])) for p in raw_points)
    elements: dict[str, tuple[tuple[int, ...], ...]] = {}
    for kind, raw in payload.get("elements", {}).items():
        if kind not in ELEMENT_KINDS:
            raise ValueError(f"unknown cell kind {kind!r}; known kinds: {', '.join(ELEMENT_KINDS)}")
        elements[kind] = _as_cells(kind, raw, len(points))
    if not elements:
        raise ValueError("mesh spec has no elements")
    return MeshSpec(points, elements, dict(payload.get("metadata", {})))


def write_spec(path: str | Path, spec: MeshSpec) -> Path:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(spec.to_dict()) + "\n", encoding="utf-8")
    return target


def spec_from_mapping(payload: Mapping[str, Any]) -> MeshSpec:
    """Build a spec from an in-memory mapping (same keys as the JSON file)."""
    import tempfile

    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as handle:
        json.dump(dict(payload), handle)
        name = handle.name
    try:
        return load_spec(name)
    finally:
        Path(name).unlink(missing_ok=True)