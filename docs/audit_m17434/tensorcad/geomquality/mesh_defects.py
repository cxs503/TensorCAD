"""Mesh-side naming of the facet-scale defects the BRep route names.

``tensorcad.geomquality`` audits a *BRep* -- patches, control nets, seams -- and names
fourteen defect codes in that vocabulary (``docs/geomquality-pipeline.md``). Downstream
solvers never see a BRep: they see a triangulation. The 2026-11-09 cross-route audit
(``docs/audit_crosscheck.md``) measured what that triangulation gives away: six topology
families survive as open/non-manifold/duplicate deltas, but ``SLIVER``, ``MICRO_EDGE``,
``CURVATURE_JUMP`` and ``TANGENT_JUMP`` all collapsed into the single sentence "the shell
is open", and the report recorded the gap in those words -- *"口径缺失（如实记录，未粉饰）"*.

This module closes the part of that gap that a mesh can close, and does so with limits
read off the corpus rather than picked. It gives three facet-scale codes:

=============  ==============================================  ==================
code           criterion                                       measured on the
                                                               corpus (clean →
                                                               injected), see
                                                               ``examples/mesh_route_naming``
=============  ==============================================  ==================
``MICRO_EDGE`` edge shorter than the weld tolerance (1e-6 m,   0.0-0.4 % → 9.6-9.9 %
               the same tolerance :mod:`tensorcad.geomquality`  of the edges
               uses for a micro edge) *and* more than 1 % of
               the edges
``SLIVER``     face area at or below ``(1e-4 L)**2`` -- the     0.0-1.2 % → 20.5-21.2 %
               BRep route's own sliver expression -- *and*      of the faces
               more than 2 % of the faces (the same 2 % the
               mesh-quality gate in :mod:`tensorcad.validation`
               uses)
``ISOLATED_PATCH``  more than one face component                 1 → 2 components
=============  ==============================================  ==================

Why shares and not counts. The BRep route fires on a *single* micro edge or sliver, which
is right for a patch: a CAD face is a design object, and one sub-tolerance edge in it is a
defect. A tessellation is not. All three reference hulls produce a few faces whose area is
at or below ``(1e-4 L)**2`` at the stem, where the surface closes to zero width (12 such
faces out of 988 on the Series 60 parent; 0 on the Wigley), and two of them produce a
handful of sub-micrometre edges at the same place. A count-based detector would therefore
fire on a *clean* hull, which is the failure mode this repository treats as worse than a
miss. The share limits above sit an order of magnitude above the worst clean measurement
and an order of magnitude below the weakest injected one; both ends are printed by the
example so a reader can see the margin rather than trust it.

What the module deliberately does **not** do:

* It does not re-implement the topology codes -- ``FREE_EDGE``, ``NON_MANIFOLD``,
  ``ORIENTATION``, ``DUPLICATE_FACE``, ``OPEN_TOP`` belong to
  :func:`tensorcad.geometry_audit.audit_watertightness`, and duplicating that vocabulary
  would create exactly the two-dialects problem the cross-route audit was written about.
  Open edges and non-manifold edges are reported as *context* numbers, never as codes.
* It does not claim to name ``CURVATURE_JUMP`` or ``TANGENT_JUMP``. Both produce a mesh
  whose facet statistics are indistinguishable at the corpus density (2214 interior edges,
  2210 both; 376 vs 408 edges above 5 degrees), because the information that separates them
  lives in the second derivative of the surface, which a 16x8 triangulation does not carry.
  The example reports them as detected-but-unnamed (``ISOLATED_PATCH`` is what is left
  once ``MICRO_EDGE`` and ``SLIVER`` are excluded), which is a true statement about the
  tessellation and not a claim about the patch.
* It does not repair geometry. :func:`repair_facet_defects` removes facet-scale damage and
  *refuses* when the repair would open the shell -- a repair that trades a sliver for a
  hole is not a repair.
"""
from __future__ import annotations

import math
from collections import defaultdict
from dataclasses import dataclass
from typing import Any, Mapping

from .. import mesh as meshlib
from ..geometry import Mesh

__all__ = [
    "DEFAULT_SLIVER_AREA_FACTOR",
    "DEFAULT_SLIVER_SHARE_LIMIT",
    "DEFAULT_MICRO_EDGE_SHARE_LIMIT",
    "DEFAULT_WELD_TOLERANCE",
    "MeshDefect",
    "MeshDefectReport",
    "RepairResult",
    "name_mesh_defects",
    "repair_facet_defects",
]

#: Weld tolerance in metres -- the value the BRep route calls a micro edge and the CAD
#: import path calls a weld. A caller working in millimetres must pass their own.
DEFAULT_WELD_TOLERANCE = 1e-6
#: Sliver area is ``(factor * L)**2``, the expression ``tensorcad.geomquality.brep`` uses.
DEFAULT_SLIVER_AREA_FACTOR = 1e-4
#: Share of faces at or below the sliver area above which ``SLIVER`` is reported.
DEFAULT_SLIVER_SHARE_LIMIT = 0.02
#: Share of edges at or below the weld tolerance above which ``MICRO_EDGE`` is reported.
DEFAULT_MICRO_EDGE_SHARE_LIMIT = 0.01

#: Repair hint per code, in the BRep route's register (what to do, not what is wrong).
REPAIR_HINTS: Mapping[str, str] = {
    "MICRO_EDGE": "weld the vertex pair at the weld tolerance, then drop what collapsed",
    "SLIVER": "collapse the near-zero-area facet to an edge or re-tessellate the patch",
    "ISOLATED_PATCH": "sew the island back into the shell (or delete it if it is debris)",
}


@dataclass(frozen=True)
class MeshDefect:
    """One mesh-side defect: the code, how many facets carry it, and the number behind it."""

    code: str
    count: int
    measure: float
    limit: float
    detail: str

    @property
    def repair(self) -> str:
        return REPAIR_HINTS.get(self.code, "")

    def to_dict(self) -> dict[str, Any]:
        return {
            "code": self.code,
            "count": self.count,
            "measure": self.measure,
            "limit": self.limit,
            "detail": self.detail,
            "repair": self.repair,
        }


@dataclass(frozen=True)
class MeshDefectReport:
    """Codes plus the metrics they were read off, so a verdict can be re-derived."""

    defects: tuple[MeshDefect, ...]
    metrics: Mapping[str, Any]
    declared_length: float | None
    weld_tolerance: float
    sliver_area_limit: float | None
    sliver_share_limit: float
    micro_edge_share_limit: float

    @property
    def codes(self) -> tuple[str, ...]:
        return tuple(sorted({defect.code for defect in self.defects}))

    @property
    def clean(self) -> bool:
        return not self.defects

    def named(self, code: str) -> bool:
        return any(defect.code == code for defect in self.defects)

    def defect(self, code: str) -> MeshDefect | None:
        for item in self.defects:
            if item.code == code:
                return item
        return None

    def to_dict(self) -> dict[str, Any]:
        return {
            "codes": list(self.codes),
            "defects": [defect.to_dict() for defect in self.defects],
            "metrics": dict(self.metrics),
            "declared_length": self.declared_length,
            "weld_tolerance": self.weld_tolerance,
            "sliver_area_limit": self.sliver_area_limit,
            "sliver_share_limit": self.sliver_share_limit,
            "micro_edge_share_limit": self.micro_edge_share_limit,
        }


def _edge_faces(mesh: Mesh) -> dict[tuple[int, int], list[int]]:
    table: dict[tuple[int, int], list[int]] = defaultdict(list)
    for index, face in enumerate(mesh.faces):
        a, b, c = face
        for key in ((a, b), (b, c), (c, a)):
            table[(key[0], key[1]) if key[0] < key[1] else (key[1], key[0])].append(index)
    return table


def _component_sizes(mesh: Mesh, edge_faces: Mapping[tuple[int, int], list[int]]) -> list[int]:
    """Face components under *shared edge* adjacency, largest first.

    Shared vertices deliberately do not join components: two patches that touch at a single
    point are two surfaces a solver has to handle, not one.
    """
    adjacency: dict[int, set[int]] = defaultdict(set)
    for faces in edge_faces.values():
        for a in faces:
            for b in faces:
                if a != b:
                    adjacency[a].add(b)
    seen: set[int] = set()
    sizes: list[int] = []
    for start in range(len(mesh.faces)):
        if start in seen:
            continue
        stack = [start]
        seen.add(start)
        size = 0
        while stack:
            current = stack.pop()
            size += 1
            for neighbour in adjacency[current]:
                if neighbour not in seen:
                    seen.add(neighbour)
                    stack.append(neighbour)
        sizes.append(size)
    sizes.sort(reverse=True)
    return sizes


def name_mesh_defects(
    mesh: Mesh,
    *,
    declared_length: float | None = None,
    weld_tolerance: float = DEFAULT_WELD_TOLERANCE,
    sliver_area_factor: float = DEFAULT_SLIVER_AREA_FACTOR,
    sliver_share_limit: float = DEFAULT_SLIVER_SHARE_LIMIT,
    micro_edge_share_limit: float = DEFAULT_MICRO_EDGE_SHARE_LIMIT,
) -> MeshDefectReport:
    """Name the facet-scale defects of a triangulation, in the BRep route's vocabulary.

    ``declared_length`` is the model's length in metres, the same declared quantity the
    cross-route audit found the mesh route cannot work without: it fixes the sliver area
    ``(factor * L)**2``. Without it the sliver criterion is *not* evaluated and the report
    says so (``sliver_area_limit is None``), rather than silently inventing a scale.
    """
    if not mesh.faces:
        raise ValueError("name_mesh_defects needs at least one face")
    if weld_tolerance <= 0.0:
        raise ValueError("weld_tolerance must be positive")
    if declared_length is not None and declared_length <= 0.0:
        raise ValueError("declared_length must be positive")

    edge_faces = _edge_faces(mesh)
    lengths = [
        math.dist(mesh.vertices[a], mesh.vertices[b]) for (a, b) in edge_faces
    ]
    areas = list(meshlib.triangle_areas(mesh))
    components = _component_sizes(mesh, edge_faces)

    micro_edges = [length for length in lengths if length <= weld_tolerance]
    micro_share = len(micro_edges) / len(lengths) if lengths else 0.0

    sliver_area_limit = (
        (sliver_area_factor * declared_length) ** 2 if declared_length is not None else None
    )
    if sliver_area_limit is None:
        sliver_faces: list[float] = []
        sliver_share = 0.0
    else:
        sliver_faces = [area for area in areas if area <= sliver_area_limit]
        sliver_share = len(sliver_faces) / len(areas)

    open_edges = sum(1 for faces in edge_faces.values() if len(faces) == 1)
    nonmanifold_edges = sum(1 for faces in edge_faces.values() if len(faces) > 2)

    defects: list[MeshDefect] = []
    if micro_share > micro_edge_share_limit:
        defects.append(
            MeshDefect(
                code="MICRO_EDGE",
                count=len(micro_edges),
                measure=micro_share,
                limit=micro_edge_share_limit,
                detail=(
                    f"{len(micro_edges)} edges <= {weld_tolerance:g} m "
                    f"(shortest {min(micro_edges):.3g} m), "
                    f"{100.0 * micro_share:.3f} % of {len(lengths)} edges "
                    f"against a {100.0 * micro_edge_share_limit:g} % limit -- below the "
                    f"weld tolerance, so a solver sees a zero-length edge"
                ),
            )
        )
    if sliver_area_limit is not None and sliver_share > sliver_share_limit:
        defects.append(
            MeshDefect(
                code="SLIVER",
                count=len(sliver_faces),
                measure=sliver_share,
                limit=sliver_share_limit,
                detail=(
                    f"{len(sliver_faces)} faces with area <= {sliver_area_limit:.3g} m^2 "
                    f"(smallest {min(sliver_faces):.3g} m^2), "
                    f"{100.0 * sliver_share:.3f} % of {len(areas)} faces "
                    f"against a {100.0 * sliver_share_limit:g} % limit -- "
                    f"a facet thinner than the surface it represents"
                ),
            )
        )
    if len(components) > 1:
        minor = components[1:]
        defects.append(
            MeshDefect(
                code="ISOLATED_PATCH",
                count=len(minor),
                measure=float(len(components)),
                limit=1.0,
                detail=(
                    f"{len(components)} face components (largest {components[0]} faces, "
                    f"others {minor}) -- a surface island: the shell is not one connected "
                    f"surface, and the class of damage behind the island is not recoverable "
                    f"from the facets"
                ),
            )
        )

    metrics: dict[str, Any] = {
        "faces": len(mesh.faces),
        "vertices": len(mesh.vertices),
        "edges": len(lengths),
        "interior_edges": len(lengths) - open_edges - nonmanifold_edges,
        "open_edges": open_edges,
        "nonmanifold_edges": nonmanifold_edges,
        "components": components,
        "micro_edges": len(micro_edges),
        "micro_edge_share": micro_share,
        "shortest_edge": min(lengths),
        "median_edge": sorted(lengths)[len(lengths) // 2],
        "sliver_faces": len(sliver_faces),
        "sliver_share": sliver_share,
        "smallest_face_area": min(areas),
        "median_face_area": sorted(areas)[len(areas) // 2],
        "surface_area": sum(areas),
        "signed_volume": meshlib.total_signed_volume(mesh),
    }
    return MeshDefectReport(
        defects=tuple(defects),
        metrics=metrics,
        declared_length=declared_length,
        weld_tolerance=weld_tolerance,
        sliver_area_limit=sliver_area_limit,
        sliver_share_limit=sliver_share_limit,
        micro_edge_share_limit=micro_edge_share_limit,
    )


@dataclass(frozen=True)
class RepairResult:
    """The repaired mesh, or the original one plus the reason the repair was refused."""

    mesh: Mesh
    accepted: bool
    reason: str
    codes_before: tuple[str, ...]
    codes_after: tuple[str, ...]
    merged_vertices: int
    removed_faces: int
    volume_relative_change: float
    area_relative_change: float
    before: MeshDefectReport
    after: MeshDefectReport

    def to_dict(self) -> dict[str, Any]:
        return {
            "accepted": self.accepted,
            "reason": self.reason,
            "codes_before": list(self.codes_before),
            "codes_after": list(self.codes_after),
            "merged_vertices": self.merged_vertices,
            "removed_faces": self.removed_faces,
            "volume_relative_change": self.volume_relative_change,
            "area_relative_change": self.area_relative_change,
        }


def _relative_change(before: float, after: float) -> float:
    scale = max(abs(before), 1e-300)
    return abs(after - before) / scale


def repair_facet_defects(
    mesh: Mesh,
    *,
    declared_length: float | None = None,
    weld_tolerance: float = DEFAULT_WELD_TOLERANCE,
    sliver_area_factor: float = DEFAULT_SLIVER_AREA_FACTOR,
    sliver_share_limit: float = DEFAULT_SLIVER_SHARE_LIMIT,
    micro_edge_share_limit: float = DEFAULT_MICRO_EDGE_SHARE_LIMIT,
    volume_tolerance: float = 1e-9,
    area_tolerance: float = 1e-9,
) -> RepairResult:
    """Weld micro edges, drop what collapsed, and refuse if that opens the shell.

    The contract is deliberately narrow, because "repair" is where a geometry pipeline
    silently changes the shape it was asked to measure:

    1. vertices closer than ``weld_tolerance`` are merged, and faces whose area falls to or
       below the sliver area are dropped -- the only two operations performed;
    2. the result is accepted only if it has **no new** codes, **no more** open edges than
       before, and changes the signed volume and surface area by less than
       ``volume_tolerance`` / ``area_tolerance`` relative;
    3. otherwise the *original* mesh is returned with ``accepted=False`` and the number that
       tripped the refusal, so a caller cannot mistake a refused repair for a silent one.

    Facet-scale damage only: a curvature jump, a ripple in the control net or a fold is a
    geometry problem, and welding vertices does not touch it.
    """
    before = name_mesh_defects(
        mesh,
        declared_length=declared_length,
        weld_tolerance=weld_tolerance,
        sliver_area_factor=sliver_area_factor,
        sliver_share_limit=sliver_share_limit,
        micro_edge_share_limit=micro_edge_share_limit,
    )
    welded = meshlib.weld(mesh, tolerance=weld_tolerance)
    drop_limit = before.sliver_area_limit if before.sliver_area_limit is not None else 0.0
    trimmed = meshlib.remove_degenerate_faces(welded, tolerance=drop_limit)
    if not trimmed.faces:
        # Dropping every face is not a repair: it would make every later measurement
        # vacuous, so it is refused before the metrics are even read.
        return RepairResult(
            mesh=mesh,
            accepted=False,
            reason="the sliver limit dropped every face -- original mesh returned",
            codes_before=before.codes,
            codes_after=before.codes,
            merged_vertices=0,
            removed_faces=0,
            volume_relative_change=0.0,
            area_relative_change=0.0,
            before=before,
            after=before,
        )
    after = name_mesh_defects(
        trimmed,
        declared_length=declared_length,
        weld_tolerance=weld_tolerance,
        sliver_area_factor=sliver_area_factor,
        sliver_share_limit=sliver_share_limit,
        micro_edge_share_limit=micro_edge_share_limit,
    )
    volume_change = _relative_change(before.metrics["signed_volume"], after.metrics["signed_volume"])
    area_change = _relative_change(before.metrics["surface_area"], after.metrics["surface_area"])
    new_codes = sorted(set(after.codes) - set(before.codes))
    opens = after.metrics["open_edges"] > before.metrics["open_edges"]

    reason_parts: list[str] = []
    if new_codes:
        reason_parts.append("repair introduced " + ", ".join(new_codes))
    if opens:
        reason_parts.append(
            f"open edges {before.metrics['open_edges']} -> {after.metrics['open_edges']}"
        )
    if volume_change > volume_tolerance:
        reason_parts.append(f"volume moved {volume_change:.3g} (limit {volume_tolerance:g})")
    if area_change > area_tolerance:
        reason_parts.append(f"area moved {area_change:.3g} (limit {area_tolerance:g})")
    leftovers = sorted(after.codes)
    if not reason_parts and leftovers:
        # Welding vertices and dropping collapsed facets cannot sew an island or re-fair a
        # surface; saying "accepted" while the code survives would read as a repair.
        reason_parts.append(
            "no facet-scale repair available for " + ", ".join(leftovers)
            + " -- that is geometry, not a facet"
        )
    accepted = not reason_parts
    if accepted:
        reason = (
            "welded "
            f"{len(mesh.vertices) - len(trimmed.vertices)} vertices and dropped "
            f"{len(mesh.faces) - len(trimmed.faces)} faces; codes "
            f"{list(before.codes)} -> {list(after.codes) if after.codes else ['clean']}"
        )
        result_mesh = trimmed
    else:
        reason = "; ".join(reason_parts) + " -- original mesh returned"
        result_mesh = mesh
    return RepairResult(
        mesh=result_mesh,
        accepted=accepted,
        reason=reason,
        codes_before=before.codes,
        codes_after=after.codes,
        merged_vertices=len(mesh.vertices) - len(trimmed.vertices),
        removed_faces=len(mesh.faces) - len(trimmed.faces),
        volume_relative_change=volume_change,
        area_relative_change=area_change,
        before=before,
        after=after,
    )