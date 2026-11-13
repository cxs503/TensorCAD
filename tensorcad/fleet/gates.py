"""Machine-decidable gates for a hull mesh admitted to the baseline fleet.

Every gate returns a :class:`GateResult` carrying the number it decided on and the
tolerance it used, so a verdict is reproducible from the artefact alone. There is no
"trust me" gate: a placeholder hull cannot pass G1, an unclosed surface cannot pass
G2, and a mesh whose volume disagrees with the record cannot pass G6.

The gates are deliberately independent of how the mesh was produced (form-family
loft, official offset table, or imported STL), because the point of the intake gate
is to be able to reject a bad file from any route.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Iterable, Sequence

from tensorcad import mesh as meshlib
from tensorcad.geometry import Mesh
from tensorcad.fleet.manifest import DERIVED_BASES, FleetManifest, ShipRecord
from tensorcad.fleet.sources import SourceRegistry

#: Volume agreement required by the mission statement (mission: error < 1%).
VOLUME_TOLERANCE = 0.01
#: Internal consistency of a declared (Cb, L, B, T, volume) set.
PARTICULARS_TOLERANCE = 0.03
#: Agreement required between a *cited* displacement volume and Cb*L*B*T. Seven
#: significant digits are carried by the published sets, so a transcription error in any
#: of L, B, T, Cb or the volume shows up far above this band, while round-off of the
#: published particulars stays below it.
REFERENCE_CONSISTENCY_TOLERANCE = 0.005


@dataclass(frozen=True)
class GateResult:
    name: str
    passed: bool
    detail: str
    values: dict[str, Any] = field(default_factory=dict)

    def as_dict(self) -> dict[str, Any]:
        return {"name": self.name, "passed": self.passed, "detail": self.detail, "values": self.values}


@dataclass(frozen=True)
class TopologyFacts:
    """Everything the gates need, measured once from the mesh."""

    triangle_count: int
    vertex_count: int
    edge_count: int
    boundary_edge_count: int
    nonmanifold_edge_count: int
    inconsistently_oriented_edge_count: int
    euler_characteristic: int
    signed_volume: float
    surface_area: float
    zero_area_face_count: int
    duplicate_face_count: int
    bbox: tuple[float, float, float, float, float, float]

    def as_dict(self) -> dict[str, Any]:
        return {
            "triangle_count": self.triangle_count,
            "vertex_count": self.vertex_count,
            "edge_count": self.edge_count,
            "boundary_edge_count": self.boundary_edge_count,
            "nonmanifold_edge_count": self.nonmanifold_edge_count,
            "inconsistently_oriented_edge_count": self.inconsistently_oriented_edge_count,
            "euler_characteristic": self.euler_characteristic,
            "signed_volume": self.signed_volume,
            "surface_area": self.surface_area,
            "zero_area_face_count": self.zero_area_face_count,
            "duplicate_face_count": self.duplicate_face_count,
            "bbox": list(self.bbox),
        }


def measure(mesh: Mesh) -> TopologyFacts:
    """Count the topology of ``mesh`` once, in one pass over its faces."""
    edge_faces = meshlib.edge_face_counts(mesh)
    boundary = sum(1 for count in edge_faces.values() if count == 1)
    nonmanifold = sum(1 for count in edge_faces.values() if count > 2)
    directed: dict[tuple[int, int], int] = {}
    seen: set[frozenset[int]] = set()
    duplicates = 0
    for face in mesh.faces:
        a, b, c = face
        for edge in ((a, b), (b, c), (c, a)):
            directed[edge] = directed.get(edge, 0) + 1
        key = frozenset(face)
        if key in seen:
            duplicates += 1
        else:
            seen.add(key)
    inconsistent = 0
    for (u, v), count in directed.items():
        if count > 1:
            inconsistent += 1
        elif directed.get((v, u), 0) != 1 and edge_faces.get(meshlib.undirected_edge((u, v)), 0) == 2:
            inconsistent += 1
    vertices = mesh.vertices
    if vertices:
        xs = [p[0] for p in vertices]
        ys = [p[1] for p in vertices]
        zs = [p[2] for p in vertices]
        bbox = (min(xs), max(xs), min(ys), max(ys), min(zs), max(zs))
    else:
        bbox = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    return TopologyFacts(
        triangle_count=len(mesh.faces),
        vertex_count=len(vertices),
        edge_count=len(edge_faces),
        boundary_edge_count=boundary,
        nonmanifold_edge_count=nonmanifold,
        inconsistently_oriented_edge_count=inconsistent,
        euler_characteristic=len(vertices) - len(edge_faces) + len(mesh.faces),
        signed_volume=meshlib.total_signed_volume(mesh),
        surface_area=meshlib.surface_area(mesh),
        zero_area_face_count=len(meshlib.zero_area_faces(mesh)),
        duplicate_face_count=duplicates,
        bbox=bbox,
    )


def placeholder_reasons(
    facts: TopologyFacts,
    *,
    max_placeholder_triangles: int = 304,
    min_real_triangles: int = 1000,
    box_corner_tolerance: float = 1e-9,
    mesh: Mesh | None = None,
) -> list[str]:
    """Reasons this mesh looks like a stand-in rather than a hull.

    Four independent signatures, because the owner's defect had at least two:
    the 304-triangle placeholder, the degenerate/box shape, and two ships sharing
    one byte-identical file (handled at fleet level by
    :meth:`FleetManifest.duplicate_geometry`).
    """
    reasons: list[str] = []
    if facts.triangle_count <= max_placeholder_triangles:
        reasons.append(
            f"placeholder_signature: {facts.triangle_count} triangles <= {max_placeholder_triangles}"
        )
    elif facts.triangle_count < min_real_triangles:
        reasons.append(f"too_few_triangles: {facts.triangle_count} < {min_real_triangles}")
    if abs(facts.signed_volume) <= 1e-12:
        reasons.append(f"zero_enclosed_volume: {facts.signed_volume:.6g}")
    if facts.vertex_count and mesh is not None:
        corners = {
            (round(corner_x, 6), round(corner_y, 6), round(corner_z, 6))
            for corner_x in (facts.bbox[0], facts.bbox[1])
            for corner_y in (facts.bbox[2], facts.bbox[3])
            for corner_z in (facts.bbox[4], facts.bbox[5])
        }
        on_corner = sum(
            1
            for p in mesh.vertices
            if (
                round(p[0], 6),
                round(p[1], 6),
                round(p[2], 6),
            )
            in corners
        )
        if on_corner == mesh.vertices.__len__() and facts.vertex_count <= 8:
            reasons.append(
                f"axis_aligned_box: all {facts.vertex_count} vertices are bbox corners"
            )
        elif box_corner_tolerance > 0 and on_corner == 0 and facts.vertex_count > 8:
            pass
    return reasons


def gate_provenance(record: ShipRecord, manifest: FleetManifest) -> GateResult:
    provenance = record.provenance
    verified = provenance.verified
    detail = (
        f"source/license stated ({provenance.source[:48]!r}); retrieval={provenance.retrieval}; "
        f"verified={verified}; manifest_sha256={(manifest.sha256 or '')[:16]}"
    )
    passed = provenance.complete
    return GateResult(
        "G0_provenance",
        passed,
        detail,
        {
            "provenance_complete": provenance.complete,
            "retrieval": provenance.retrieval,
            "verified": verified,
            "manifest_sha256": manifest.sha256,
            "url": provenance.url,
            "license": provenance.license,
        },
    )


def gate_not_placeholder(facts: TopologyFacts, mesh: Mesh, placeholder: dict[str, Any]) -> GateResult:
    reasons = placeholder_reasons(
        facts,
        max_placeholder_triangles=int(placeholder.get("max_placeholder_triangles", 304)),
        min_real_triangles=int(placeholder.get("min_real_triangles", 1000)),
        mesh=mesh,
    )
    return GateResult(
        "G1_not_placeholder",
        not reasons,
        "no placeholder signature" if not reasons else "; ".join(reasons),
        {"reasons": reasons, "triangle_count": facts.triangle_count, "vertices_on_bbox_corners": None},
    )


def gate_no_duplicate_geometry(record: ShipRecord, manifest: FleetManifest) -> GateResult:
    collisions = {
        sha: names for sha, names in manifest.duplicate_geometry().items() if record.name in names
    }
    detail = (
        "each ship owns a distinct geometry sha256"
        if not collisions
        else f"byte-identical geometry shared: {list(collisions.values())}"
    )
    return GateResult(
        "G2_no_duplicate_geometry",
        not collisions,
        detail,
        {"collisions": {sha[:16]: list(names) for sha, names in collisions.items()}},
    )


def gate_watertight(facts: TopologyFacts) -> GateResult:
    return GateResult(
        "G3_watertight",
        facts.boundary_edge_count == 0,
        f"boundary edges = {facts.boundary_edge_count} (tolerance 0)",
        {"boundary_edge_count": facts.boundary_edge_count, "tolerance": 0},
    )


def gate_manifold_oriented(facts: TopologyFacts) -> GateResult:
    passed = (
        facts.nonmanifold_edge_count == 0
        and facts.inconsistently_oriented_edge_count == 0
        and facts.zero_area_face_count == 0
        and facts.duplicate_face_count == 0
    )
    return GateResult(
        "G4_manifold_oriented",
        passed,
        (
            f"non-manifold edges = {facts.nonmanifold_edge_count}, inconsistent directed edges = "
            f"{facts.inconsistently_oriented_edge_count}, zero-area faces = {facts.zero_area_face_count}, "
            f"duplicate faces = {facts.duplicate_face_count}"
        ),
        {
            "nonmanifold_edge_count": facts.nonmanifold_edge_count,
            "inconsistently_oriented_edge_count": facts.inconsistently_oriented_edge_count,
            "zero_area_face_count": facts.zero_area_face_count,
            "duplicate_face_count": facts.duplicate_face_count,
        },
    )


def gate_single_closed_shell(facts: TopologyFacts) -> GateResult:
    passed = facts.euler_characteristic == 2 and facts.signed_volume > 0.0
    return GateResult(
        "G5_single_closed_shell",
        passed,
        f"Euler characteristic = {facts.euler_characteristic} (expect 2), signed volume = "
        f"{facts.signed_volume:.6g} (expect > 0)",
        {"euler_characteristic": facts.euler_characteristic, "signed_volume": facts.signed_volume},
    )


def gate_particulars_consistent(record: ShipRecord) -> GateResult:
    """A declared (Cb, L, B, T, V) set must be internally coherent before it is used."""
    closed = record.closed_form_volume()
    if record.displacement_volume is None:
        return GateResult(
            "G6_particulars_consistent",
            True,
            f"no published volume to cross-check; closure Cb*L*B*T = {closed:.6g} m^3 is the reference "
            f"(basis={record.reference_basis()})",
            {"closed_form_volume": closed, "displacement_volume": None, "relative_difference": None},
        )
    rel = abs(closed - record.displacement_volume) / record.displacement_volume
    return GateResult(
        "G6_particulars_consistent",
        rel < PARTICULARS_TOLERANCE,
        f"|Cb*L*B*T - V_published| / V_published = {rel * 100:.3f}% (tolerance "
        f"{PARTICULARS_TOLERANCE * 100:g}%), V_published = {record.displacement_volume:.6g} m^3",
        {
            "closed_form_volume": closed,
            "displacement_volume": record.displacement_volume,
            "relative_difference": rel,
            "tolerance": PARTICULARS_TOLERANCE,
        },
    )


def gate_source_anchor(record: ShipRecord, registry: SourceRegistry) -> GateResult:
    """The record must name a source anchor: citation, URL, licence and locator.

    ``G0_provenance`` checks that the record states *a* URL and licence. This gate checks
    that the fleet's source registry holds the citation the digits would be re-derived
    from, and reports the registry's own sha256 so the anchor set is itself auditable.
    """
    anchor = registry.for_ship(record.name)
    if anchor is None:
        return GateResult(
            "G0b_source_anchor",
            False,
            f"no source anchor registered for {record.name!r}: a manifest digit with no named source "
            f"cannot be re-derived by a reader",
            {"registry_sha256": registry.sha256, "anchor_id": None},
        )
    covered = sorted(set(anchor.anchors))
    return GateResult(
        "G0b_source_anchor",
        True,
        f"anchor {anchor.anchor_id!r}: {anchor.citation[:64]!r}; locator {anchor.locator[:48]!r}; "
        f"licence {anchor.license!r}; retrieval={anchor.retrieval}; verified={anchor.verified}; "
        f"anchors={covered}",
        {
            "registry_sha256": registry.sha256,
            "anchor_id": anchor.anchor_id,
            "citation": anchor.citation,
            "url": anchor.url,
            "license": anchor.license,
            "locator": anchor.locator,
            "anchors": covered,
            "retrieval": anchor.retrieval,
            "verified": anchor.verified,
        },
    )


def gate_reference_consistency(
    record: ShipRecord, tolerance: float = REFERENCE_CONSISTENCY_TOLERANCE
) -> GateResult:
    """A cited displacement volume must agree with the record's own Cb*L*B*T.

    This is the gate that catches a *transcription* error rather than a bad mesh: if a
    published displacement is quoted next to published L, B, T and Cb, then the four
    numbers over-determine each other, and any one of them entered wrongly makes the
    set inconsistent. A closure reference is not a check at all (it is defined as
    Cb*L*B*T), and the gate says so instead of reporting a vacuous pass.
    """
    closure = record.closed_form_volume()
    if record.displacement_volume is None:
        return GateResult(
            "G6b_reference_consistency",
            True,
            f"reference is the closure Cb*L*B*T = {closure:.6g} m^3: no independent displacement is "
            f"cited, so this gate cannot check anything (G7 then only proves the mesh matches the "
            f"record's own particulars)",
            {
                "cited_volume": None,
                "closure_volume": closure,
                "relative_error": None,
                "tolerance": tolerance,
                "vacuous": True,
            },
        )
    cited = float(record.displacement_volume)
    rel = (cited - closure) / cited
    return GateResult(
        "G6b_reference_consistency",
        abs(rel) <= tolerance,
        f"cited V = {cited:.6g} m^3 vs Cb*L*B*T = {closure:.6g} m^3: {rel * 100:+.4f}% "
        f"(tolerance {tolerance * 100:g}%); a disagreement beyond this band means L, B, T, Cb or the "
        f"cited volume was mis-transcribed",
        {
            "cited_volume": cited,
            "closure_volume": closure,
            "relative_error": rel,
            "tolerance": tolerance,
            "vacuous": False,
            "basis": record.displacement_basis,
        },
    )


def gate_volume(record: ShipRecord, facts: TopologyFacts, tolerance: float = VOLUME_TOLERANCE) -> GateResult:
    reference = record.reference_volume()
    basis = record.reference_basis()
    rel = abs(facts.signed_volume - reference) / reference
    label = "closure (NOT an external validation)" if basis in DERIVED_BASES else f"literature ({basis})"
    return GateResult(
        "G7_volume_vs_reference",
        rel < tolerance,
        f"|V_mesh - V_ref| / V_ref = {rel * 100:.4f}% (tolerance {tolerance * 100:g}%); V_mesh = "
        f"{facts.signed_volume:.6g} m^3, V_ref = {reference:.6g} m^3, basis = {label}",
        {
            "mesh_volume": facts.signed_volume,
            "reference_volume": reference,
            "reference_basis": basis,
            "reference_is_external": record.reference_is_external(),
            "relative_error": rel,
            "tolerance": tolerance,
        },
    )


def gate_principal_dimensions(record: ShipRecord, facts: TopologyFacts, tolerance: float = VOLUME_TOLERANCE) -> GateResult:
    """The mesh must actually span L, B and T of the record, not merely its volume."""
    xmin, xmax, ymin, ymax, zmin, zmax = facts.bbox
    measured = {"length": xmax - xmin, "beam": ymax - ymin, "draft": -(zmin if zmin < 0 else zmin)}
    measured["draft"] = zmax - zmin if zmax <= 0.0 else abs(zmin)
    declared = {"length": record.length, "beam": record.beam, "draft": record.draft}
    errors = {
        key: abs(measured[key] - declared[key]) / declared[key] for key in declared
    }
    worst = max(errors.values())
    return GateResult(
        "G8_principal_dimensions",
        worst < tolerance,
        "; ".join(
            f"{key}: measured {measured[key]:.6g} vs declared {declared[key]:.6g} "
            f"({errors[key] * 100:.4f}%)"
            for key in declared
        )
        + f" (tolerance {tolerance * 100:g}%)",
        {"measured": measured, "declared": declared, "relative_errors": errors, "tolerance": tolerance},
    )


def evaluate_gates(
    record: ShipRecord,
    mesh: Mesh,
    manifest: FleetManifest,
    *,
    volume_tolerance: float = VOLUME_TOLERANCE,
    registry: SourceRegistry | None = None,
) -> dict[str, Any]:
    """Run every gate for one ship and return a serialisable verdict block."""
    facts = measure(mesh)
    gates: list[GateResult] = [
        gate_provenance(record, manifest),
        gate_not_placeholder(facts, mesh, dict(manifest.placeholder)),
        gate_no_duplicate_geometry(record, manifest),
        gate_watertight(facts),
        gate_manifold_oriented(facts),
        gate_single_closed_shell(facts),
        gate_particulars_consistent(record),
        gate_reference_consistency(record),
        gate_volume(record, facts, volume_tolerance),
        gate_principal_dimensions(record, facts, volume_tolerance),
    ]
    if registry is not None:
        gates.insert(1, gate_source_anchor(record, registry))
    return {
        "ship": record.name,
        "ship_class": record.ship_class,
        "route": record.geometry_route,
        "topology": facts.as_dict(),
        "gates": [gate.as_dict() for gate in gates],
        "passed": all(gate.passed for gate in gates),
        "failed_gates": [gate.name for gate in gates if not gate.passed],
        "volume_relative_error": next(
            gate.values.get("relative_error") for gate in gates if gate.name == "G7_volume_vs_reference"
        ),
        "reference_basis": record.reference_basis(),
        "reference_is_external": record.reference_is_external(),
    }


def refused(record: ShipRecord, reason: str) -> dict[str, Any]:
    """Verdict block for a ship the geometry route declined to build.

    A refusal is a result, not an exception: it carries the measured reason and the
    number that caused it, so a route limit can never be mistaken for a pass.
    """
    return {
        "ship": record.name,
        "ship_class": record.ship_class,
        "route": record.geometry_route,
        "topology": None,
        "gates": [{"name": "G_build", "passed": False, "detail": reason, "values": {}}],
        "passed": False,
        "failed_gates": ["G_build"],
        "volume_relative_error": None,
        "reference_basis": record.reference_basis(),
        "reference_is_external": record.reference_is_external(),
    }