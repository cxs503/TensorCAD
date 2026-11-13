"""Independent references for :mod:`tensorcad.midsurface`: an O(n^2) pairing oracle.

The module under test accelerates its wall search with a **uniform bucket grid**
(:class:`tensorcad.midsurface._FaceGrid`): for each face it collects the triangles of a
small neighbourhood of cells and ray-casts only those. Acceleration structures are the
classic place where a geometry kernel loses a hit -- a candidate pruned at a cell
boundary is invisible in every area or thickness check that follows, because the loss
shows up as an "unpaired face", which is a legal outcome of the very same function. A
test that only compares areas against a closed form therefore cannot see it.

So this module never imports ``_FaceGrid``, never imports ``_ray_triangle`` and never
uses Moller-Trumbore. It answers the same question the brute-force way:

* every face is compared with **every** other face (no spatial pruning at all), and
* the ray/triangle test is written as *plane intersection + barycentric position from
  subtriangle areas*, a different formulation from Moller-Trumbore's determinant
  solve, so a sign or branch error in one is not reproduced in the other.

Two policies are offered because the docstring of :func:`extract_midsurface` fixes one
of them and only one:

``pairs_first_hit``
    The documented policy: take the *nearest surface the ray strikes*, then accept it
    as a wall partner only if it is antiparallel within the angle tolerance and its
    distance is inside the thickness band. A closer but non-antiparallel surface
    therefore *shadows* the wall behind it and the face stays unpaired.
``pairs_nearest_antiparallel``
    The idealised reading -- the nearest surface that *would* qualify. The difference
    between the two dictionaries is exactly the set of faces hidden behind another
    surface; it is measured (not assumed) by
    ``tests/test_midsurface_oracle.py``.

The closed forms (plate / exact offset cylinder / exact offset sphere) are analytic and
share no code with the module: they exist so the extracted *area* has a reference that
is not another mesh operation.

Standard library only, like the package itself.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

from tensorcad.geometry import Face, Mesh, Vec3

Vec = tuple[float, float, float]


# --------------------------------------------------------------------------- #
# vector helpers (deliberately separate from the module's)
# --------------------------------------------------------------------------- #
def _sub(a: Vec, b: Vec) -> Vec:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _dot(a: Vec, b: Vec) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _cross(a: Vec, b: Vec) -> Vec:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _norm(a: Vec) -> float:
    return math.sqrt(_dot(a, a))


def _unit(a: Vec) -> Vec:
    length = _norm(a)
    if length == 0.0:
        raise ValueError("zero-length vector")
    return (a[0] / length, a[1] / length, a[2] / length)


def triangle_area(a: Vec, b: Vec, c: Vec) -> float:
    return 0.5 * _norm(_cross(_sub(b, a), _sub(c, a)))


def face_normal(vertices: Sequence[Vec], face: Face) -> Vec:
    """Unnormalised normal by the right-hand rule (zero for a degenerate face)."""
    return _cross(
        _sub(vertices[face[1]], vertices[face[0]]),
        _sub(vertices[face[2]], vertices[face[0]]),
    )


def face_centroid(vertices: Sequence[Vec], face: Face) -> Vec:
    return (
        (vertices[face[0]][0] + vertices[face[1]][0] + vertices[face[2]][0]) / 3.0,
        (vertices[face[0]][1] + vertices[face[1]][1] + vertices[face[2]][1]) / 3.0,
        (vertices[face[0]][2] + vertices[face[1]][2] + vertices[face[2]][2]) / 3.0,
    )


# --------------------------------------------------------------------------- #
# ray / triangle, by plane intersection + area barycentrics
# --------------------------------------------------------------------------- #
def ray_triangle_distance(
    origin: Vec, direction: Vec, a: Vec, b: Vec, c: Vec, *, area_relative_tolerance: float = 1e-12
) -> float | None:
    """Distance along ``direction`` (unit length) to the triangle, or ``None``.

    The triangle's plane is intersected first; the hit point is then located by
    comparing the areas of the three subtriangles it forms with the corners against
    the whole triangle's area. The point is inside iff the three subtriangle areas add
    up to the whole one -- a different test from Moller-Trumbore's edge equations.
    """
    normal = _cross(_sub(b, a), _sub(c, a))
    denominator = _dot(direction, normal)
    if denominator == 0.0:
        return None  # ray parallel to the triangle's plane
    distance = _dot(_sub(a, origin), normal) / denominator
    if distance <= 0.0:
        return None  # the triangle is behind the origin
    hit = (origin[0] + direction[0] * distance,
           origin[1] + direction[1] * distance,
           origin[2] + direction[2] * distance)
    total = triangle_area(a, b, c)
    if total == 0.0:
        return None
    parts = triangle_area(hit, b, c) + triangle_area(a, hit, c) + triangle_area(a, b, hit)
    if abs(parts - total) > area_relative_tolerance * total + 1e-300:
        return None  # outside the triangle
    return distance


# --------------------------------------------------------------------------- #
# the two pairing oracles
# --------------------------------------------------------------------------- #
def _face_frames(mesh: Mesh) -> list[tuple[Vec, Vec]]:
    """``(centroid, unit outward normal)`` per face; the normal is zero if degenerate."""
    frames: list[tuple[Vec, Vec]] = []
    for face in mesh.faces:
        centroid = face_centroid(mesh.vertices, face)
        raw = face_normal(mesh.vertices, face)
        frames.append((centroid, _unit(raw) if _norm(raw) > 0.0 else (0.0, 0.0, 0.0)))
    return frames


def _scan(
    mesh: Mesh,
    *,
    max_thickness: float,
    min_thickness: float,
    angle_tolerance_deg: float,
    require_antiparallel: bool,
) -> dict[int, tuple[int, float]]:
    """All-pairs scan; nearest hit (optionally nearest *qualifying* hit) per face."""
    frames = _face_frames(mesh)
    cos_limit = math.cos(math.radians(angle_tolerance_deg))
    scale = max(_norm(_sub(max(mesh.vertices), min(mesh.vertices))), 1e-12)
    result: dict[int, tuple[int, float]] = {}
    for index, (centroid, normal) in enumerate(frames):
        if _norm(normal) == 0.0:
            continue
        inward = (-normal[0], -normal[1], -normal[2])
        best: tuple[float, int] | None = None
        for candidate, candidate_face in enumerate(mesh.faces):
            if candidate == index:
                continue
            if require_antiparallel:
                other = frames[candidate][1]
                if _norm(other) == 0.0 or _dot(normal, other) > -cos_limit:
                    continue
            distance = ray_triangle_distance(
                centroid, inward,
                mesh.vertices[candidate_face[0]],
                mesh.vertices[candidate_face[1]],
                mesh.vertices[candidate_face[2]],
            )
            if distance is None or distance <= 1e-12 * scale:
                continue
            if best is None or distance < best[0]:
                best = (distance, candidate)
        if best is None:
            continue
        distance, partner = best
        if distance < min_thickness or distance > max_thickness:
            continue
        result[index] = (partner, distance)
    return result


def pairs_first_hit(
    mesh: Mesh,
    *,
    max_thickness: float,
    min_thickness: float = 0.0,
    angle_tolerance_deg: float = 20.0,
) -> dict[int, tuple[int, float]]:
    """The *documented* policy: nearest struck surface, then the band and normal filters."""
    return _scan(
        mesh,
        max_thickness=max_thickness,
        min_thickness=min_thickness,
        angle_tolerance_deg=angle_tolerance_deg,
        require_antiparallel=False,
    )


def pairs_nearest_antiparallel(
    mesh: Mesh,
    *,
    max_thickness: float,
    min_thickness: float = 0.0,
    angle_tolerance_deg: float = 20.0,
) -> dict[int, tuple[int, float]]:
    """The idealised policy: nearest surface that is antiparallel *and* inside the band."""
    return _scan(
        mesh,
        max_thickness=max_thickness,
        min_thickness=min_thickness,
        angle_tolerance_deg=angle_tolerance_deg,
        require_antiparallel=True,
    )


def shadowed_faces(
    mesh: Mesh,
    *,
    max_thickness: float,
    min_thickness: float = 0.0,
    angle_tolerance_deg: float = 20.0,
) -> tuple[int, ...]:
    """Faces that a closer, non-antiparallel surface hides from its wall partner."""
    first = pairs_first_hit(
        mesh, max_thickness=max_thickness, min_thickness=min_thickness,
        angle_tolerance_deg=angle_tolerance_deg,
    )
    ideal = pairs_nearest_antiparallel(
        mesh, max_thickness=max_thickness, min_thickness=min_thickness,
        angle_tolerance_deg=angle_tolerance_deg,
    )
    return tuple(sorted(set(ideal) - set(first)))


# --------------------------------------------------------------------------- #
# analytic references: areas the extraction is allowed to be compared against
# --------------------------------------------------------------------------- #
def plate_mid_area(length: float, breadth: float) -> float:
    """A flat wall's mid-surface is the mid-plane: ``L * B`` exactly."""
    return length * breadth


def cylinder_offset_area(radius: float, thickness: float, length: float) -> float:
    """Exact offset surface area of a tube skin at ``t/2``: ``2 pi (R - t/2) L``."""
    return 2.0 * math.pi * (radius - 0.5 * thickness) * length


def cylinder_surplus(radius: float, thickness: float) -> float:
    """Area surplus of per-facet translation of a cylinder: ``(t/2) * kappa = t/(2R)``."""
    return 0.5 * thickness * (1.0 / radius)


def sphere_offset_area(radius: float, thickness: float) -> float:
    """Exact offset surface area of a spherical skin at ``t/2``: ``4 pi (R - t/2)^2``."""
    return 4.0 * math.pi * (radius - 0.5 * thickness) ** 2


def sphere_surplus(radius: float, thickness: float) -> float:
    """Area surplus of per-facet translation of a sphere: ``(t/2) * (2/R)``."""
    return 0.5 * thickness * (2.0 / radius)


def mesh_area(mesh: Mesh) -> float:
    return sum(triangle_area(mesh.vertices[f[0]], mesh.vertices[f[1]], mesh.vertices[f[2]])
               for f in mesh.faces)


def relative_error(measured: float, reference: float) -> float:
    """``(measured - reference) / reference`` -- signed, so a shortfall is not hidden."""
    if reference == 0.0:
        raise ValueError("no relative error against a zero reference")
    return (measured - reference) / reference


# --------------------------------------------------------------------------- #
# relabelling, to expose face-order dependence
# --------------------------------------------------------------------------- #
def relabel(mesh: Mesh, permutation: Sequence[int]) -> Mesh:
    """Rebuild ``mesh`` with ``face i`` becoming ``face permutation[i]``.

    Geometry is untouched -- only the order in which an algorithm visits faces changes.
    A per-face operation must not care about that order; the mid-surface emission rule
    ("one triangle per pair") is exactly the kind of bookkeeping that can.
    """
    if sorted(permutation) != list(range(len(mesh.faces))):
        raise ValueError("permutation must be a bijection of the face indices")
    faces = [mesh.faces[permutation[i]] for i in range(len(mesh.faces))]
    return Mesh(vertices=mesh.vertices, faces=tuple(faces), units=mesh.units,
                coordinate_frame=mesh.coordinate_frame)