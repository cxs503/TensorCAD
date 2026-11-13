"""Hard-fail machine gate for mesh artifacts.

A surface that looks like a hull in a viewer is not necessarily a surface a
solver can use. Two failures that this module exists to catch are on record in
this repository: a **half model** (open along the symmetry plane, so it has no
inside and no displacement) and a **low-resolution grid** (a few hundred
facets, whose integrated volume is several percent off the analytic value even
though every topology check passes). Neither is visible from a thumbnail, and a
triangle-soup format such as STL cannot even represent the topology that would
make them obvious.

:func:`tensorcad.validation.inspect_geometry` *reports* these conditions;
:func:`gate_mesh` is the same idea as a gate: every structural finding is
promoted to an error, component count / Euler characteristic / optional
reference volume / resolution are added, and :func:`enforce_gate` refuses to
return a mesh that does not pass. Every reference it compares against comes
from outside the mesh (analytic displacement, requested component count, the
genus-0 Euler characteristic), so a self-consistent but wrong export cannot
satisfy the gate.

Standard library only, like the rest of the core.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

from . import mesh as meshlib
from .geometry import Mesh
from .validation import ERROR, Check, inspect_geometry

#: Default tolerance for the integrated volume against a reference value.
DEFAULT_VOLUME_TOLERANCE = 0.02


class MeshGateError(ValueError):
    """Raised by :func:`enforce_gate` when an artifact fails the machine gate."""


@dataclass(frozen=True)
class MeshGateReport:
    """Checks and measurements produced by one gate pass."""

    checks: tuple[Check, ...]
    measurements: dict[str, Any]

    @property
    def passed(self) -> bool:
        return all(check.passed for check in self.checks)

    @property
    def failures(self) -> tuple[Check, ...]:
        return tuple(check for check in self.checks if not check.passed)

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": "tensorcad.mesh.gate/1",
            "passed": self.passed,
            "failures": [check.name for check in self.failures],
            "checks": [check.to_dict() for check in self.checks],
            "measurements": self.measurements,
        }


def _components(mesh: Mesh) -> int:
    """Number of connected components over the vertices used by faces."""
    parent: dict[int, int] = {}

    def find(vertex: int) -> int:
        root = vertex
        while parent[root] != root:
            root = parent[root]
        while parent[vertex] != root:
            parent[vertex], vertex = root, parent[vertex]
        return root

    for face in mesh.faces:
        for vertex in face:
            parent.setdefault(vertex, vertex)
        for u, v in ((face[0], face[1]), (face[1], face[2])):
            root_u, root_v = find(u), find(v)
            if root_u != root_v:
                parent[root_u] = root_v
    return len({find(vertex) for vertex in parent})


def _edge_lengths(mesh: Mesh) -> dict[tuple[int, int], float]:
    lengths: dict[tuple[int, int], float] = {}
    for face in mesh.faces:
        for a, b in ((face[0], face[1]), (face[1], face[2]), (face[2], face[0])):
            key = (a, b) if a <= b else (b, a)
            if key in lengths:
                continue
            pa, pb = mesh.vertices[a], mesh.vertices[b]
            lengths[key] = math.dist(pa, pb)
    return lengths


def _facet_aspect(mesh: Mesh) -> float:
    """Worst ``longest_edge^2 / area`` over the faces (1.155 for an equilateral).

    A sliver is a triangle with a tiny area for its longest edge, and slivers are
    what make a solver's gradient estimates ill-conditioned. The value is scale
    free, so it can be compared across hull sizes.
    """
    worst = 0.0
    for face in mesh.faces:
        ax, ay, az = mesh.vertices[face[0]]
        bx, by, bz = mesh.vertices[face[1]]
        cx, cy, cz = mesh.vertices[face[2]]
        ux, uy, uz = bx - ax, by - ay, bz - az
        vx, vy, vz = cx - ax, cy - ay, cz - az
        nx, ny, nz = uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx
        area = 0.5 * math.sqrt(nx * nx + ny * ny + nz * nz)
        if area <= 0.0:
            continue
        longest = max(
            (bx - ax) ** 2 + (by - ay) ** 2 + (bz - az) ** 2,
            (cx - bx) ** 2 + (cy - by) ** 2 + (cz - bz) ** 2,
            (ax - cx) ** 2 + (ay - cy) ** 2 + (az - cz) ** 2,
        )
        worst = max(worst, longest / area)
    return worst


def measure_gated(mesh: Mesh) -> dict[str, Any]:
    """Everything the gate compares against, in one JSON-safe record."""
    counts = meshlib.edge_face_counts(mesh)
    boundary = sum(1 for count in counts.values() if count == 1)
    non_manifold = sum(1 for count in counts.values() if count > 2)
    edges = len(counts)
    lengths = _edge_lengths(mesh)
    measurements: dict[str, Any] = {
        "units": mesh.units,
        "vertices": len(mesh.vertices),
        "triangles": len(mesh.faces),
        "edges": edges,
        "euler_characteristic": len(mesh.vertices) - edges + len(mesh.faces),
        "components": _components(mesh),
        "boundary_edges": boundary,
        "non_manifold_edges": non_manifold,
        "zero_area_faces": len(meshlib.zero_area_faces(mesh)),
        "signed_volume": meshlib.total_signed_volume(mesh),
        "surface_area": meshlib.surface_area(mesh),
        "max_edge_length": max(lengths.values()) if lengths else 0.0,
        "max_facet_aspect": _facet_aspect(mesh),
    }
    if mesh.vertices:
        box = mesh.bounding_box()
        measurements["bounding_box_size"] = list(box.sizes)
        longest = max(box.sizes)
        measurements["max_relative_edge_length"] = (
            measurements["max_edge_length"] / longest if longest > 0.0 else math.inf
        )
    return measurements


def gate_mesh(
    mesh: Mesh,
    *,
    reference_volume: float | None = None,
    volume_tolerance: float = DEFAULT_VOLUME_TOLERANCE,
    expected_components: int = 1,
    expected_euler: int | None = None,
    min_faces: int | None = None,
    max_aspect: float | None = None,
    max_relative_edge_length: float | None = None,
) -> MeshGateReport:
    """Run every machine check and return the evidence, without raising.

    ``reference_volume`` must come from outside the mesh (for the parametric hull
    :func:`tensorcad.hull.analytic_volume`); a mesh that merely agrees with
    itself is not evidence. ``expected_euler`` defaults to ``2 * components``,
    the characteristic of a closed orientable genus-0 shell, which is what a hull
    with no through-holes must be.
    """
    if volume_tolerance < 0.0:
        raise ValueError("volume_tolerance must be non-negative")
    if expected_components < 1:
        raise ValueError("expected_components must be >= 1")
    if reference_volume is not None and (
        not math.isfinite(reference_volume) or reference_volume <= 0.0
    ):
        raise ValueError("reference_volume must be a positive finite number")
    if min_faces is not None and min_faces < 1:
        raise ValueError("min_faces must be >= 1")

    measurements = measure_gated(mesh)
    target_euler = 2 * expected_components if expected_euler is None else expected_euler

    checks: list[Check] = []
    for check in inspect_geometry(mesh):
        # inspect_geometry's Euler test is genus-0 / single-component specific;
        # the gate re-derives it below against the requested values.
        if check.name == "euler_characteristic_genus_zero":
            continue
        checks.append(Check(check.name, check.passed, check.detail, severity=ERROR))

    components = int(measurements["components"])
    checks.append(
        Check(
            "component_count",
            components == expected_components,
            f"{components} connected components, expected {expected_components}",
        )
    )
    euler = int(measurements["euler_characteristic"])
    checks.append(
        Check(
            "euler_characteristic",
            euler == target_euler,
            f"V-E+F = {measurements['vertices']}-{measurements['edges']}+"
            f"{measurements['triangles']} = {euler}, expected {target_euler} "
            f"({expected_components} closed genus-0 component(s))",
        )
    )

    if reference_volume is not None:
        measured = float(measurements["signed_volume"])
        relative = measured / reference_volume - 1.0
        measurements["reference_volume"] = reference_volume
        measurements["relative_volume_error"] = relative
        checks.append(
            Check(
                "reference_volume",
                abs(relative) <= volume_tolerance,
                f"integrated {measured:.6g} {mesh.units}^3 against reference "
                f"{reference_volume:.6g} ({100.0 * relative:+.3f} %, "
                f"limit {100.0 * volume_tolerance:g} %)",
            )
        )

    if min_faces is not None:
        triangles = int(measurements["triangles"])
        checks.append(
            Check(
                "facet_count",
                triangles >= min_faces,
                f"{triangles} triangles, at least {min_faces} required "
                "(a coarse export is not a small valid export)",
            )
        )

    if max_relative_edge_length is not None:
        relative_edge = float(measurements.get("max_relative_edge_length", math.inf))
        checks.append(
            Check(
                "facet_size",
                relative_edge <= max_relative_edge_length,
                f"longest facet edge is {relative_edge:.4g} of the bounding-box "
                f"diagonal-relative length, limit {max_relative_edge_length:g}",
            )
        )

    if max_aspect is not None:
        aspect = float(measurements["max_facet_aspect"])
        checks.append(
            Check(
                "facet_aspect",
                aspect <= max_aspect,
                f"worst facet has longest_edge^2/area = {aspect:.4g}, limit {max_aspect:g}",
            )
        )

    return MeshGateReport(checks=tuple(checks), measurements=measurements)


def enforce_gate(mesh: Mesh, *, label: str = "mesh", **kwargs: Any) -> MeshGateReport:
    """:func:`gate_mesh`, but raise :class:`MeshGateError` instead of returning a failure."""
    report = gate_mesh(mesh, **kwargs)
    if report.passed:
        return report
    failures = "; ".join(
        f"{check.name}: {check.detail}" if check.detail else check.name
        for check in report.failures
    )
    raise MeshGateError(f"{label} failed the machine gate: {failures}")
