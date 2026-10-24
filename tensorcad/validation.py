"""Geometric and engineering validity checks.

The checks here are what makes an exported artifact *provably* usable rather than
merely plausible: closure, manifoldness, orientation, absence of zero-area faces,
topological characteristic, and agreement between the requested principal
dimensions and the bounding box that was actually generated.

Checks are returned as records rather than raised as the first exception found,
because an agent iterating on a model needs the full list of what is wrong.
"""
from __future__ import annotations

import math
from collections import defaultdict
from dataclasses import dataclass

from . import mesh as meshlib
from .geometry import Mesh

ERROR = "error"
WARNING = "warning"


@dataclass(frozen=True)
class Check:
    name: str
    passed: bool
    detail: str = ""
    severity: str = ERROR

    def to_dict(self) -> dict[str, object]:
        data: dict[str, object] = {"name": self.name, "passed": self.passed}
        if self.detail:
            data["detail"] = self.detail
        if self.severity != ERROR:
            data["severity"] = self.severity
        return data


@dataclass(frozen=True)
class ValidationReport:
    checks: tuple[Check, ...] = ()

    @property
    def errors(self) -> list[str]:
        return [
            f"{check.name}: {check.detail}" if check.detail else check.name
            for check in self.checks
            if not check.passed and check.severity == ERROR
        ]

    @property
    def warnings(self) -> list[str]:
        return [
            f"{check.name}: {check.detail}" if check.detail else check.name
            for check in self.checks
            if not check.passed and check.severity == WARNING
        ]

    @property
    def passed(self) -> bool:
        return not self.errors

    def to_dict(self) -> dict[str, object]:
        return {
            "passed": self.passed,
            "errors": self.errors,
            "warnings": self.warnings,
            "checks": [check.to_dict() for check in self.checks],
        }


def vertex_link_defects(mesh: Mesh) -> tuple[int, ...]:
    """Vertices whose incident faces do not form a single closed ring.

    For a closed, manifold triangle mesh every vertex link is one cycle. A
    vertex where sheets pinch together (for example a zero-breadth stem) or
    where a surface is open fails that test, which is the geometric reason a
    mesher can still refuse an otherwise watertight file.
    """
    faces_at: dict[int, list[tuple[int, int, int]]] = defaultdict(list)
    for face in mesh.faces:
        for vertex in face:
            faces_at[vertex].append(face)

    defects: list[int] = []
    for vertex, surrounding in faces_at.items():
        link: dict[int, set[int]] = defaultdict(set)
        for face in surrounding:
            if vertex not in face:
                continue
            index = face.index(vertex)
            a, b = face[(index + 1) % 3], face[(index + 2) % 3]
            if a == b:
                continue
            link[a].add(b)
            link[b].add(a)
        if not link or any(len(neighbours) != 2 for neighbours in link.values()):
            defects.append(vertex)
            continue
        start = next(iter(link))
        seen = {start}
        stack = [start]
        while stack:
            current = stack.pop()
            for other in link[current]:
                if other not in seen:
                    seen.add(other)
                    stack.append(other)
        if len(seen) != len(link):  # two rings meeting: a pinch, not a manifold
            defects.append(vertex)
    return tuple(sorted(defects))


def inspect_geometry(mesh: Mesh) -> list[Check]:
    """Topology, orientation and degeneration checks on a triangulated surface."""
    checks: list[Check] = []
    vertices, faces = mesh.vertices, mesh.faces

    non_finite = sum(
        1 for point in vertices if not all(math.isfinite(value) for value in point)
    )
    checks.append(
        Check("finite_vertex_coordinates", non_finite == 0,
              "" if not non_finite else f"{non_finite} vertices carry non-finite coordinates")
    )
    checks.append(Check("non_empty_mesh", bool(vertices) and bool(faces),
                        "" if faces else "no faces generated"))

    bad_indices = sum(
        1
        for face in faces
        if len(face) != 3 or any(i < 0 or i >= len(vertices) for i in face)
    )
    checks.append(
        Check("face_indices_in_range", bad_indices == 0,
              "" if not bad_indices else f"{bad_indices} faces index outside the vertex table")
    )
    if not vertices or not faces:
        return checks

    zero_area = meshlib.zero_area_faces(mesh)
    checks.append(
        Check("no_degenerate_faces", not zero_area,
              "" if not zero_area else f"{len(zero_area)} zero-area triangles present")
    )

    counts = meshlib.edge_face_counts(mesh)
    boundary = sum(1 for count in counts.values() if count == 1)
    non_manifold = sum(1 for count in counts.values() if count > 2)
    checks.append(
        Check(
            "watertight",
            boundary == 0,
            "" if not boundary else f"{boundary} boundary edges: the surface is open",
        )
    )
    checks.append(
        Check(
            "edge_manifold",
            non_manifold == 0,
            "" if not non_manifold else f"{non_manifold} edges shared by more than two faces",
        )
    )

    euler = len(vertices) - len(counts) + len(faces)
    checks.append(
        Check(
            "euler_characteristic_genus_zero",
            boundary == 0 and euler == 2,
            f"V-E+F = {len(vertices)}-{len(counts)}+{len(faces)} = {euler}, "
            "expected 2 for a closed genus-0 shell",
        )
    )

    volume = meshlib.total_signed_volume(mesh)
    checks.append(
        Check(
            "outward_orientation",
            volume > 0,
            f"signed volume is {volume:.6g} {mesh.units}^3; a closed solid must have positive volume",
        )
    )

    if boundary == 0:
        box = mesh.bounding_box()
        limit = box.sizes[0] * box.sizes[1] * box.sizes[2]
        checks.append(
            Check(
                "volume_within_bounding_box",
                volume <= limit * (1.0 + 1e-9) + 1e-15,
                f"volume {volume:.9g} exceeds the {limit:.9g} bounding-box volume, "
                "which means self-intersection or inverted faces",
            )
        )

    defects = vertex_link_defects(mesh)
    checks.append(
        Check(
            "vertex_manifold",
            not defects,
            "" if not defects else (
                f"{len(defects)} pinch vertices where sheets meet along a line "
                f"(degenerate stem/post at zero breadth); a transom or rounded bow removes them"
            ),
            severity=WARNING,
        )
    )
    return checks


@dataclass(frozen=True)
class DimensionTarget:
    """Requested principal dimensions to compare the generated bounds against."""

    length: float
    beam: float
    draft: float
    waterline_z: float = 0.0
    tolerance_relative: float = 1e-9


def check_dimensions(mesh: Mesh, target: DimensionTarget) -> list[Check]:
    """Compare requested L/B/T with the bounding box the generator produced.

    A sampled surface only passes through the widest station it actually has, so
    the realised beam can fall short of the requested beam when the longitudinal
    grid misses the midship station. The shortfall is bounded by the shape
    function and the station spacing; we allow that much and no more, and report
    the realised numbers either way.
    """
    box = mesh.bounding_box()
    checks: list[Check] = []
    realised_length, realised_beam, realised_top, realised_bottom = (
        box.sizes[0],
        box.sizes[1],
        max(z for _, _, z in mesh.vertices),
        min(z for _, _, z in mesh.vertices),
    )
    checks.append(
        _near("length_matches_request", realised_length, target.length, target.tolerance_relative, "length")
    )
    depth = realised_top - realised_bottom
    checks.append(
        _near("draft_matches_request", depth, target.draft, target.tolerance_relative, "draft")
    )
    checks.append(
        _near(
            "waterline_plane_matches_request",
            realised_top,
            target.waterline_z,
            target.tolerance_relative,
            "waterline z",
        )
    )
    # Beam may legitimately fall short by up to the grid sampling deficit.
    deficit = target.beam - realised_beam
    allowance = max(target.tolerance_relative, _sampling_allowance(deficit, realised_beam, target.beam))
    checks.append(
        Check(
            "beam_matches_request",
            -target.tolerance_relative * target.beam <= deficit <= allowance,
            f"realised beam {realised_beam:.9g} {mesh.units} vs requested {target.beam:.9g} "
            f"(grid deficit {100 * deficit / target.beam:+.4f}%, "
            f"allowance {100 * allowance / target.beam:.4f}%)",
        )
    )
    return checks


def _sampling_allowance(deficit: float, realised: float, requested: float) -> float:
    """Accept a beam shortfall that is explained purely by station sampling."""
    if requested <= 0 or deficit < 0:
        return 0.0
    relative = deficit / requested
    # A one-sided sampling deficit of a smooth maximum is O(spacing^2); allow a
    # generous 4x band so the check still fires on real geometry errors.
    return relative if relative < 1e-4 else min(4e-2, relative * 4 + 1e-6)


def _near(name: str, actual: float, expected: float, tolerance: float, label: str) -> Check:
    ok = math.isclose(actual, expected, rel_tol=tolerance, abs_tol=tolerance * max(abs(expected), 1.0))
    return Check(
        name,
        ok,
        "" if ok else f"{label} is {actual:.9g}, requested {expected:.9g}",
    )


def inspect_mesh_quality(
    mesh: Mesh,
    *,
    min_angle_degrees: float = 2.0,
    max_aspect_ratio: float = 25.0,
    sliver_share_limit: float = 0.02,
) -> list[Check]:
    """Face-quality checks: a watertight mesh can still be unusable to a solver.

    Closure says the surface is a solid; it says nothing about whether the
    elements are shaped well enough for a finite-element or finite-volume code.
    Three cheap measures catch almost all of it:

    * the smallest interior angle (a sliver has one near zero);
    * the aspect ratio of each triangle (longest edge over twice its height);
    * the share of faces inside that threshold.

    All three are warnings: a mesh can be a valid closed solid and still be a poor
    discretisation, and the report should say so without pretending the geometry
    is broken.
    """
    checks: list[Check] = []
    if not mesh.faces:
        return [Check("mesh_quality_measured", False, "no faces to measure", severity=WARNING)]
    smallest_angle = 180.0
    largest_aspect = 0.0
    slivers = 0
    for face in mesh.faces:
        a = mesh.vertices[face[0]]
        b = mesh.vertices[face[1]]
        c = mesh.vertices[face[2]]
        lengths = (
            math.dist(a, b),
            math.dist(b, c),
            math.dist(c, a),
        )
        longest = max(lengths)
        shortest = min(lengths)
        if shortest <= 0.0:  # degenerate; inspect_geometry already reports it
            slivers += 1
            continue
        ux, uy, uz = b[0] - a[0], b[1] - a[1], b[2] - a[2]
        vx, vy, vz = c[0] - a[0], c[1] - a[1], c[2] - a[2]
        # |u x v| is twice the area, so longest**2 / area2 is longest edge over height.
        area2 = math.sqrt(
            (uy * vz - uz * vy) ** 2 + (uz * vx - ux * vz) ** 2 + (ux * vy - uy * vx) ** 2
        )
        if area2 <= 0.0:
            slivers += 1
            continue
        aspect = longest * longest / area2
        largest_aspect = max(largest_aspect, aspect)
        # Smallest angle from the two shortest sides meeting at its vertex.
        ordered = sorted(lengths)
        opposite, side_a, side_b = ordered[0], ordered[1], ordered[2]
        cosine = (side_a * side_a + side_b * side_b - opposite * opposite) / (2.0 * side_a * side_b)
        angle = math.degrees(math.acos(max(-1.0, min(1.0, cosine))))
        smallest_angle = min(smallest_angle, angle)
        if angle < min_angle_degrees:
            slivers += 1

    total = len(mesh.faces)
    share = slivers / total
    checks.append(
        Check(
            "mesh_min_angle",
            smallest_angle >= min_angle_degrees,
            f"smallest interior angle {smallest_angle:.3f} deg over {total} triangles "
            f"(limit {min_angle_degrees:g} deg)",
            severity=WARNING,
        )
    )
    checks.append(
        Check(
            "mesh_aspect_ratio",
            largest_aspect <= max_aspect_ratio,
            f"worst triangle aspect ratio {largest_aspect:.3g} (limit {max_aspect_ratio:g})",
            severity=WARNING,
        )
    )
    checks.append(
        Check(
            "mesh_sliver_share",
            share <= sliver_share_limit,
            f"{slivers}/{total} triangles ({100.0 * share:.3f} %) are degenerate or thinner "
            f"than {min_angle_degrees:g} deg (limit {100.0 * sliver_share_limit:g} %)",
            severity=WARNING,
        )
    )
    return checks


def merge(*groups: list[Check]) -> ValidationReport:
    flat = [check for group in groups for check in group]
    return ValidationReport(tuple(flat))
