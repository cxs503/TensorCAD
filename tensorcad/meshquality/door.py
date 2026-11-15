"""Admission control: one verdict out of topology, element quality and face non-orthogonality.

Why this module exists
----------------------
:mod:`tensorcad.meshquality.report` gates *elements* and
:mod:`tensorcad.meshquality.topology` gates the *shell*; they were shipped as two
separate answers, and a pipeline that wanted one door had to know both call
signatures and re-implement the combination. Board audit ``docs/mesh-quality-lib-audit-m17404.md``
finding F2 recorded exactly that: the gate "exists as ``quality_gate``/``gate_report``
returning a verdict dict … so it *can* be a door", but nothing composed it.

:func:`admit` is that composition, and nothing else: it adds no measure and no
threshold of its own. Every :class:`~tensorcad.meshquality.report.Check` it returns is
the object the underlying gate produced, with its own name, limit and detail, so a
failure still says *which* threshold on *which* cell kind failed by *how much*. The
three families it asks for are the three that make a mesh unusable for a solver:

* **shell** (triangles only) -- watertight, edge/vertex-manifold, consistent winding,
  outward orientation, optional genus bound;
* **elements** -- the :data:`~tensorcad.meshquality.report.DEFAULT_GATE` limits on the
  cell kinds present (``tri``, ``quad``, ``tet``, ``hex``, ``prism``, ``pyramid``);
* **interior faces** (``tet``/``hex``) -- the OpenFOAM ``checkMesh`` non-orthogonality
  angle between the cell-centre vector and the face normal.

A gate on a cell kind the mesh does not contain stays vacuous
(:meth:`~tensorcad.meshquality.report.QualityGate.evaluate` skips it); the door keeps
that behaviour, so a surface-only mesh is not failed for having no hexahedra.

This module is deliberately a *composition*: if the verdict is wrong, the bug is in one
of the three gates and is reproducible there. ``tools/mesh_quality_door.py`` is the CI
entry point on top of it.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

from .report import DEFAULT_GATE, Check, QualityGate, measure
from .topology import (
    FaceNonOrthogonality,
    SurfaceTopology,
    gate_non_orthogonality,
    gate_surface,
    non_orthogonality,
    surface_topology,
)

__all__ = ["CELL_KINDS_WITH_INTERNAL_FACES", "DoorVerdict", "admit"]

Vec3 = tuple[float, float, float]

#: Cell kinds whose interior faces the non-orthogonality measure supports. Prisms and
#: pyramids are in the list because a boundary-layer prism column and a hex/tet transition
#: pyramid both have interior faces the solver cares about; leaving them out made ``admit``
#: skip the measure for those meshes, which reads as a PASS in the verdict.
CELL_KINDS_WITH_INTERNAL_FACES = ("tet", "hex", "prism", "pyramid")


@dataclass(frozen=True)
class DoorVerdict:
    """Every check the door applied, in the order it applied them."""

    checks: tuple[Check, ...]
    payload: Mapping[str, Any]

    @property
    def passed(self) -> bool:
        return all(check.passed for check in self.checks)

    @property
    def failures(self) -> tuple[Check, ...]:
        return tuple(check for check in self.checks if not check.passed)

    def to_dict(self) -> dict[str, Any]:
        return {
            "passed": self.passed,
            "checks": [check.to_dict() for check in self.checks],
            "failures": [check.name for check in self.failures],
            "measurements": dict(self.payload),
        }

    def write(self, path: str | Path) -> Path:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(self.to_dict(), indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return target


def admit(
    points: Sequence[Sequence[float]],
    *,
    triangles: Sequence[Sequence[int]] | None = None,
    cells: Mapping[str, Sequence[Sequence[int]]] | None = None,
    gate: QualityGate = DEFAULT_GATE,
    max_non_orthogonality: float = 70.0,
    require_watertight: bool = True,
    require_manifold: bool = True,
    require_outward: bool = True,
    max_genus: float | None = None,
) -> DoorVerdict:
    """Every check the repository's canonical library can apply to one mesh.

    ``triangles`` are the shell faces (0-based node indices); ``cells`` maps a cell kind
    (``"tet"``, ``"hex"``, ``"prism"``, ``"pyramid"``, or ``"quad"`` for a 2-D panel) to
    its connectivity. Either may be omitted; omitting both is an error rather than a
    silent pass, because a door that admits an empty mesh is not a door.

    ``max_non_orthogonality`` is in degrees and follows the ``checkMesh`` convention
    (``70`` is the usual warning level, ``75`` the usual failure level).
    """
    if not triangles and not cells:
        raise ValueError("admit() needs triangles, cells, or both; an empty mesh is not admitted")

    checks: list[Check] = []
    payload: dict[str, Any] = {}

    if triangles:
        topology: SurfaceTopology = surface_topology(points, triangles)
        checks.extend(
            gate_surface(
                topology,
                require_watertight=require_watertight,
                require_manifold=require_manifold,
                require_outward=require_outward,
                max_genus=max_genus,
            )
        )
        surface_quality = measure(points, {"tri": triangles})
        surface_report = gate.evaluate(surface_quality)
        checks.extend(surface_report.checks)
        payload["topology"] = topology.to_dict()
        payload["summary"] = surface_quality.to_dict()

    if cells:
        volume_quality = measure(points, dict(cells))
        for kind in CELL_KINDS_WITH_INTERNAL_FACES:
            if kind not in cells:
                continue
            faces: FaceNonOrthogonality = non_orthogonality(points, cells[kind], kind)
            checks.extend(gate_non_orthogonality(faces, max_non_orthogonality))
            payload[f"non_orthogonality_{kind}"] = faces.to_dict()
        checks.extend(gate.evaluate(volume_quality).checks)
        if "summary" in payload:  # a mesh carrying both a shell and cells
            payload["summary"].update(volume_quality.to_dict())
        else:
            payload["summary"] = volume_quality.to_dict()

    return DoorVerdict(checks=tuple(checks), payload=payload)