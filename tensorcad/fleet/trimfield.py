"""Trim-aware sampling of IGES B-spline patches for the half-breadth field route.

Why this exists
---------------
The ``iges_field`` route (:func:`tensorcad.fleet.build.iges_field_mesh`) used to sample every
type-128 B-spline surface over its **full parametric domain** and ignore the type-142/144 trim
loops that say which part of that domain is actually ship.  For an untrimmed export (the KCS
SVA file) the basis surface *is* the hull, so nothing showed.  For the KVLCC2 archive -- an
IBM CATIA export carrying 38 basis surfaces wrapped in 38 *trimmed* patches -- it does not
hold: the basis rectangles run past the hull, so the half-breadth field picked up
``max |y|`` from sheet that is not hull surface.  Measured on the published KVLCC2 file: the
untouched route put ``|y| ~ 29 m`` (full half-beam) at every depth from the keel to the
waterline, i.e. a slab, and returned ``∇ = 132 294 m³`` against the published 312 622 m³
(``-57.7 %``) -- a trim defect, not a hull or reference defect.

What the module does
--------------------
* :func:`point_in_polygon` -- even-odd ray cast, the containment test in ``(u, v)``.
* :func:`trim_polygons` -- per trimmed patch, decide which of the two curve pointers of each
  type-142 entity lives in parameter space (``IT`` field ordering is not reliable in the
  wild) by measuring the fraction of the candidate polygon that falls inside the basis
  surface's own ``(u, v)`` domain; the winner is recorded, never defaulted.
* :func:`sample_trimmed_patches` -- sample each basis surface on a grid dense enough for
  ``spacing`` metres of chord and keep only the samples whose ``(u, v)`` is inside the trim
  loops (outer loop plus any island loops, combined with the even-odd rule).

Everything returned is traceable: each patch reports its DE pointers, the chosen field, the
inside fraction and how many samples survived, so a reviewer can recompute any number.
"""

from __future__ import annotations

import math
from typing import Any, Iterable, Mapping, Sequence

Point3 = tuple[float, float, float]
Point2 = tuple[float, float]

#: A candidate trim polygon must fall inside the basis surface's domain this often to be
#: believed; the loser pointer is a model-space curve and lands mostly outside.
MIN_INSIDE_FRACTION = 0.9


def point_in_polygon(x: float, y: float, polygon: Sequence[Point2]) -> bool:
    """Even-odd ray cast in ``(u, v)``.  Boundary points are not decided (return False)."""
    inside = False
    count = len(polygon)
    if count < 3:
        return False
    j = count - 1
    for i in range(count):
        xi, yi = polygon[i]
        xj, yj = polygon[j]
        if (yi > y) != (yj > y):
            t = (y - yi) / (yj - yi) if yj != yi else 0.0
            if x < xi + t * (xj - xi):
                inside = not inside
        j = i
    return inside


def _inside_fraction(polygon: Sequence[Point2], u0: float, u1: float, v0: float, v1: float) -> float:
    if not polygon:
        return 0.0
    hits = sum(1 for u, v in polygon if u0 - 1e-9 <= u <= u1 + 1e-9 and v0 - 1e-9 <= v <= v1 + 1e-9)
    return hits / len(polygon)


def _curve_pointer_candidates(curves: Iterable[Any]) -> list[int]:
    """Both curve pointers of each type-142 entity, de-duplicated, in file order."""
    seen: list[int] = []
    for curve in curves:
        for attribute in ("curve_uv", "curve_3d", "curve_3d_crv", "curve"):
            pointer = getattr(curve, attribute, None)
            if isinstance(pointer, int) and pointer > 0 and pointer not in seen:
                seen.append(pointer)
    return seen


def trim_polygons(iges: Any, *, samples: int = 200) -> list[dict[str, Any]]:
    """One record per trimmed patch: basis surface, chosen pointer, ``(u, v)`` loops.

    ``iges`` is an :class:`tensorcad.iges.IgesFile`.  Patches whose basis surface or trim
    pointer cannot be resolved are returned with ``verdict == "undecided"`` and no polygons,
    so the caller can decide to fall back rather than silently sampling a whole basis sheet.
    """
    records: list[dict[str, Any]] = []
    for index, (wrapper, basis, curves) in enumerate(iges.trimmed_surfaces()):
        record: dict[str, Any] = {
            "patch": index,
            "surface_de": getattr(wrapper, "surface", None),
            "boundary_de": getattr(wrapper, "boundary", None),
            "basis_bound": basis is not None,
            "candidate_pointers": _curve_pointer_candidates(curves),
            "pointer": None,
            "inside_fraction": None,
            "loops": [],
            "verdict": "undecided",
        }
        if basis is None or not record["candidate_pointers"]:
            records.append(record)
            continue
        u0, u1 = basis.domain_u
        v0, v1 = basis.domain_v
        scores: dict[int, tuple[float, list[Point2]]] = {}
        for pointer in record["candidate_pointers"]:
            try:
                polygon = iges.trim_polygon_uv(pointer, samples=samples)
            except Exception:  # a pointer that is not a curve at all
                continue
            scores[pointer] = (_inside_fraction(polygon, u0, u1, v0, v1), polygon)
        if not scores:
            records.append(record)
            continue
        best_pointer, (best_fraction, _) = max(scores.items(), key=lambda item: item[1][0])
        record["inside_fraction"] = round(best_fraction, 4)
        if best_fraction < MIN_INSIDE_FRACTION:
            record["basis"] = basis
            records.append(record)
            continue
        loops: list[list[Point2]] = []
        for pointer in record["candidate_pointers"]:
            fraction, polygon = scores[pointer]
            if fraction >= MIN_INSIDE_FRACTION and polygon:
                loops.append(polygon)
        record["pointer"] = best_pointer
        record["loops"] = loops
        record["verdict"] = "decided"
        record["basis"] = basis
        records.append(record)
    return records


def sample_trimmed_patches(
    records: Sequence[Mapping[str, Any]],
    *,
    spacing: float = 1.0,
) -> tuple[list[Point3], dict[str, Any]]:
    """Sample only the trimmed part of every basis surface.

    ``spacing`` is a chord length in model units; the number of samples per direction comes
    from the patch's control-point extent, so a 300 m hull gets a fine grid and a small
    patch does not waste points.
    """
    points: list[Point3] = []
    used = 0
    skipped = 0
    kept_samples = 0
    dropped_samples = 0
    per_patch: list[dict[str, Any]] = []
    for record in records:
        basis = record.get("basis")
        loops = record.get("loops") or []
        if basis is None or not loops:
            skipped += 1
            continue
        (x0, y0, z0), (x1, y1, z1) = basis.bounding_box()
        extent = math.dist((x0, y0, z0), (x1, y1, z1))
        count = max(4, min(240, int(math.ceil(max(extent, 1e-9) / max(spacing, 1e-9)))))
        u0, u1 = basis.domain_u
        v0, v1 = basis.domain_v
        patch_kept = 0
        for j in range(count + 1):
            v = v0 + (v1 - v0) * j / count
            for i in range(count + 1):
                u = u0 + (u1 - u0) * i / count
                hits = sum(1 for loop in loops if point_in_polygon(u, v, loop))
                if hits % 2 == 0:
                    dropped_samples += 1
                    continue
                points.append(basis.evaluate(u, v))
                patch_kept += 1
        kept_samples += patch_kept
        used += 1
        per_patch.append({"patch": record["patch"], "grid": count + 1, "kept": patch_kept})
    report = {
        "mode": "trimmed",
        "patches_used": used,
        "patches_skipped": skipped,
        "kept_samples": kept_samples,
        "dropped_samples": dropped_samples,
        "per_patch": per_patch,
    }
    return points, report