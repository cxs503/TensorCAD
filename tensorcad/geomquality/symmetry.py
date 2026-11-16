"""Symmetry closure: turn a half-hull patch set into a full, watertight hull.

Why this stage exists
---------------------
Publisher ship IGES files (KCS, KVLCC2, ONRT, ...) are commonly exported as the
*starboard half* of the hull: the surface ends on the centre plane ``y = const``, which is
a modelling cut, not a leak.  No local repair may close that boundary -- capping it would
build a wall through the ship -- and the watertight gate can therefore never pass, however
good the seam stage is.  Measured on ``work/assets/kcs.igs`` at a 20 mm seam budget,
58.6 % of the surviving open-boundary length lies on the centre plane; that is the number
this module acts on.

The stage detects a symmetry plane on the *open boundary* (not on a guess about the file),
mirrors the mesh across it, reverses the winding of the mirrored copy and welds the seam,
so the seam vertices are shared rather than duplicated.  Detection is required to be
one-sided: a plane that also has material on its far side is not a symmetry cut, and the
stage refuses instead of producing a doubled solid.

Everything the stage decided travels in its ``record``: which plane, what fraction of the
boundary it carried, how many vertices were snapped, and the open boundary before/after.
"""
from __future__ import annotations

from typing import Any, Sequence

from .stitch import free_edges

Vec3 = tuple[float, float, float]
Triangle = tuple[int, int, int]

AXES = {0: "x", 1: "y", 2: "z"}

#: A plane must carry at least this fraction of the open-boundary *length* to be accepted
#: as the model's symmetry cut.  The floor is deliberately below one half: on a publisher
#: hull the centre plane carries the long cut chords, but the defect gaps between patches
#: add length elsewhere, so the true cut measures 47 % on ``kcs.igs`` (over 8x the runner-up
#: candidate, see ``docs/iges_pipe_symmetry_closure.md``).  Dominance over the runner-up,
#: below, is what separates a cut from a coincidence.
DEFAULT_MIN_COVERAGE = 0.4
#: The winning plane must carry at least this multiple of the runner-up's length share.
DEFAULT_MIN_DOMINANCE = 3.0


def _distance(vertex: Vec3, plane: tuple[int, float]) -> float:
    axis, value = plane
    return abs(vertex[axis] - value)


def _edge_length(a: Vec3, b: Vec3) -> float:
    return sum((a[i] - b[i]) ** 2 for i in range(3)) ** 0.5


def detect_symmetry_plane(vertices: Sequence[Vec3], faces: Sequence[Triangle], *,
                          band: float | None = None,
                          min_coverage: float = DEFAULT_MIN_COVERAGE,
                          min_dominance: float = DEFAULT_MIN_DOMINANCE
                          ) -> dict[str, Any] | None:
    """Return the plane the open boundary lies on, or ``None``.

    ``band`` is the distance in model units within which a boundary edge counts as lying
    on the plane; it defaults to ``1e-3`` of the model extent.  A candidate plane is only
    accepted when (a) it carries at least ``min_coverage`` of the open-boundary length and
    (b) the mesh does not cross it -- every vertex is on one side, so mirroring cannot
    overlap existing material.
    """
    edges = free_edges(faces)
    if not edges or not vertices:
        return None
    lo = [min(vertex[axis] for vertex in vertices) for axis in range(3)]
    hi = [max(vertex[axis] for vertex in vertices) for axis in range(3)]
    extent = max(hi[axis] - lo[axis] for axis in range(3))
    if extent <= 0.0:
        return None
    if band is None:
        band = 1e-3 * extent

    lengths = [_edge_length(vertices[a], vertices[b]) for a, b in edges]
    total = sum(lengths)
    if total <= 0.0:
        return None

    candidates: list[dict[str, Any]] = []
    for axis in range(3):
        for side, value in (("lo", lo[axis]), ("hi", hi[axis])):
            # one-sidedness: no vertex may lie more than ``band`` beyond the plane
            if side == "lo":
                crossing = sum(1 for vertex in vertices if vertex[axis] < value - band)
            else:
                crossing = sum(1 for vertex in vertices if vertex[axis] > value + band)
            if crossing:
                continue
            carried = sum(length for (a, b), length in zip(edges, lengths)
                          if abs(0.5 * (vertices[a][axis] + vertices[b][axis]) - value) <= band)
            coverage = carried / total
            candidates.append({
                "axis": axis,
                "axis_name": AXES[axis],
                "value": value,
                "side": side,
                "coverage": coverage,
                "carried_length": carried,
                "band": band,
                "extent": extent,
            })
    if not candidates:
        return None
    candidates.sort(key=lambda item: item["coverage"], reverse=True)
    best = candidates[0]
    runner_up = candidates[1]["coverage"] if len(candidates) > 1 else 0.0
    best["runner_up_coverage"] = runner_up
    if best["coverage"] < min_coverage:
        return None
    if runner_up > 0.0 and best["coverage"] < min_dominance * runner_up:
        return None
    return best


def mirror_across_plane(vertices: Sequence[Vec3], faces: Sequence[Triangle], *,
                        plane: tuple[int, float], weld_tolerance: float
                        ) -> tuple[list[Vec3], list[Triangle], dict[str, Any]]:
    """Mirror the mesh across an axis plane, reverse the copy's winding and weld the seam.

    Vertices within ``weld_tolerance`` of the plane are snapped *onto* it first, so the
    seam closes exactly instead of leaving a sub-tolerance crack: the mirrored copy then
    reuses the snapped vertex instead of duplicating it.  Reflection reverses
    handedness, so the mirrored triangles are emitted with reversed winding to keep the
    normals outward.
    """
    from .trim_mesh import weld_vertices

    axis, value = plane
    snapped = 0
    base: list[Vec3] = []
    for vertex in vertices:
        if abs(vertex[axis] - value) <= weld_tolerance:
            point = list(vertex)
            point[axis] = value
            base.append(tuple(point))
            snapped += 1
        else:
            base.append(tuple(vertex))

    # a plane vertex is shared with its mirror image; an off-plane vertex gets a new index
    index: dict[tuple[int, ...], int] = {}
    for position, point in enumerate(base):
        index.setdefault(tuple(point), position)
    out_vertices = list(base)
    out_faces: list[Triangle] = [tuple(face) for face in faces]

    mirrored = 0
    for face in faces:
        a, b, c = (base[i] for i in face)
        if a[axis] == value and b[axis] == value and c[axis] == value:
            continue  # degenerate: a face entirely on the cut plane has no mirror
        image: list[int] = []
        for point in (a, b, c):
            reflected = list(point)
            reflected[axis] = 2.0 * value - point[axis]
            key = tuple(reflected)
            position = index.get(key)
            if position is None:
                position = len(out_vertices)
                out_vertices.append(key)
                index[key] = position
            image.append(position)
        # reflect then reverse: the mirror of (a, b, c) with outward normals is (a', c', b')
        out_faces.append((image[0], image[2], image[1]))
        mirrored += 1

    merged_vertices = 0
    if weld_tolerance > 0.0:
        out_vertices, out_faces, merged = weld_vertices(out_vertices, out_faces, weld_tolerance)
        merged_vertices = merged

    record = {
        "applied": True,
        "axis": plane[0],
        "axis_name": AXES.get(plane[0], str(plane[0])),
        "value": value,
        "weld_tolerance": weld_tolerance,
        "vertices_snapped_to_plane": snapped,
        "triangles_mirrored": mirrored,
        "vertices_merged_after_mirror": merged_vertices,
        "free_edges_before": len(free_edges(faces)),
        "free_edges_after": len(free_edges(out_faces)),
    }
    return out_vertices, out_faces, record


def symmetry_stage(vertices: Sequence[Vec3], faces: Sequence[Triangle], *,
                   enabled: bool, weld_tolerance: float,
                   band: float | None = None,
                   min_coverage: float = DEFAULT_MIN_COVERAGE,
                   min_dominance: float = DEFAULT_MIN_DOMINANCE,
                   ) -> tuple[list[Vec3], list[Triangle], dict[str, Any]]:
    """Detect-and-mirror in one step; when ``enabled`` is false, only measure.

    Measurement-only is the default because mirroring changes the model's extent: a caller
    must be able to ask "is this a half hull?" without being handed a full hull.
    """
    plane = detect_symmetry_plane(vertices, faces, band=band, min_coverage=min_coverage,
                                 min_dominance=min_dominance)
    record: dict[str, Any] = {"requested": enabled, "detected": plane}
    if not enabled:
        record["applied"] = False
        record["note"] = "measurement only; pass symmetry_closure=True to mirror the half"
        return list(vertices), [tuple(face) for face in faces], record
    if plane is None:
        record["applied"] = False
        record["note"] = ("no plane carries a dominant share of the open boundary; "
                          "mesh kept")
        return list(vertices), [tuple(face) for face in faces], record
    out_vertices, out_faces, mirror_record = mirror_across_plane(
        vertices, faces, plane=(plane["axis"], plane["value"]),
        weld_tolerance=weld_tolerance)
    record.update(mirror_record)
    record["coverage"] = plane["coverage"]
    return out_vertices, out_faces, record