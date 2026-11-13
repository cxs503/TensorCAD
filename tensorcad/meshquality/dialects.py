"""The three aspect-ratio dialects this repository speaks, related by exact algebra.

Why this module exists
----------------------
"aspect ratio" is not one number. Three of them are live in this tree, and two of
them are produced by modules that each call themselves canonical:

============ ============================================ ===================================
dialect      definition                                   produced by
============ ============================================ ===================================
``legacy``   ``l_max**2 / (2 A)``                         ``mesh_conditioning``,
                                                          ``validation.measure_quality``
``edge_form````l_max * p / (4 sqrt(3) A)``                ``geomquality.quality``
``verdict``  ``R / (2 r)``                                ``meshquality.metrics``
============ ============================================ ===================================

For the same 45-45-90 triangle (the diagonal half of a unit square) they read
``2.0``, ``1.3938469`` and ``1.2071068`` -- a spread of 1.66x, so a threshold copied
from one dialect to another does not mean the same thing. This was found independently
twice: as finding F1 of ``docs/mesh-quality-lib-audit-m17404.md`` (legacy vs
``edge_form``) and, while writing this module, as a disagreement between the two
canonical modules themselves (``edge_form`` vs ``verdict``, ratio ``sqrt(3) R / l_max``).

The relations are exact, not empirical. With ``p`` the perimeter, ``l_max`` the longest
edge, ``A`` the area and ``R`` the circumradius:

* ``edge_form = legacy * alpha``, ``alpha = p / (2 sqrt(3) l_max)``
* ``verdict   = edge_form * beta``, ``beta = sqrt(3) R / l_max``
* ``verdict   = legacy * alpha * beta``

and the bounds follow from three facts of Euclidean geometry:

* ``l_max <= p/2`` and ``l_max >= p/3`` give ``alpha in [1/sqrt(3), sqrt(3)/2]``
  (``0.5773502691896258`` to ``0.8660254037844386``), the minimum at a degenerate
  collinear triangle and the maximum at the equilateral one;
* ``l_max <= 2 R`` gives ``beta >= sqrt(3)/2``, the minimum at a right triangle whose
  hypotenuse is the diameter (the 45-45-90 case above) and no upper bound, because a
  sliver has an unbounded circumradius;
* both minima are attainable separately, so ``verdict / legacy >= 1/2`` is the safe
  combined bound.

All three equal 1 for the equilateral triangle -- which is why "ideal value 1" is the
one statement every dialect agrees on, and why a normalised-ideal test cannot tell them
apart.

What the bounds are good for
----------------------------
A gate that is *stricter or equal* to the shipped one can be built by proving the
transfer, not by copying the number. :func:`conservative_limit` returns the limit in a
target dialect that can never pass an element the source gate would fail:

* source ``legacy`` 25 -> ``edge_form`` ``25/sqrt(3) = 14.4338`` -> ``verdict`` ``12.5``.

Anything above those (for instance the 25 that finding F1 quotes) *can* pass elements
the legacy gate rejects; :func:`dialect_audit` counts exactly how many on a given mesh
instead of arguing about it.

Provenance honesty
------------------
Which of ``edge_form`` and ``verdict`` VTK's ``vtkMeshQuality`` calls "Aspect Ratio"
(and therefore what ``meshmagick`` reports) could not be checked on this host: ``vtk``,
``gmsh`` and ``meshmagick`` are all absent and the loong64 package index carries none of
them (see ``docs/mesh_quality_dialects.md``). Both are exposed under neutral names here;
the attribution is recorded as unverified, not asserted.
"""
from __future__ import annotations

import math
from typing import Iterable, Mapping, Sequence

__all__ = [
    "ASPECT_DIALECTS",
    "ALPHA_MAX",
    "ALPHA_MIN",
    "BETA_MIN",
    "conservative_limit",
    "dialect_audit",
    "aspect_dialects",
    "tri_geometry",
    "legacy_aspect_from_metrics",
    "transfer_ratio_bounds",
]

Vec3 = tuple[float, float, float]

SQRT3 = math.sqrt(3.0)

#: ``alpha = edge_form / legacy`` range, attained at collinear / equilateral triangles.
ALPHA_MIN = 1.0 / SQRT3
ALPHA_MAX = SQRT3 / 2.0

#: ``beta = verdict / edge_form`` lower bound, attained by any right triangle whose
#: hypotenuse is the circumdiameter (e.g. 45-45-90). No upper bound: slivers.
BETA_MIN = SQRT3 / 2.0

#: ``verdict / legacy`` lower bound (product of two bounds that are not attained
#: simultaneously, so the true infimum is strictly larger; this is the safe statement).
COMBINED_MIN = 0.5

_INF = math.inf

#: Machine-readable description of each dialect: definition, producer, ideal value.
ASPECT_DIALECTS: Mapping[str, Mapping[str, object]] = {
    "legacy": {
        "definition": "l_max**2 / (2*A)",
        "reading": "longest edge over twice the height of that edge",
        "produced_by": ("tensorcad.mesh_conditioning._aspect_ratio", "tensorcad.validation.measure_quality"),
        # l_max**2/(2A) = l**2/(2*(sqrt(3)/4) l**2) = 2/sqrt(3) on the equilateral
        # triangle: this dialect is NOT normalised to 1, so its limit 25 sits on a
        # scale whose ideal element reads 1.1547.
        "ideal": 2.0 / SQRT3,
        "shipped_limit": 25.0,
    },
    "edge_form": {
        "definition": "l_max * (l_a + l_b + l_c) / (4*sqrt(3)*A)",
        "reading": "edge-length form of the Verdict 'Aspect Ratio' label (attribution unverified here)",
        "produced_by": ("tensorcad.geomquality.quality.triangle_quality",),
        "ideal": 1.0,
        "shipped_limit": 25.0,
    },
    "verdict": {
        "definition": "R / (2*r)",
        "reading": "circumradius over twice the inradius",
        "produced_by": ("tensorcad.meshquality.metrics.tri_metrics['aspect_ratio']",),
        "ideal": 1.0,
        "shipped_limit": None,
    },
}


def _sub(a: Vec3, b: Vec3) -> Vec3:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _cross(a: Vec3, b: Vec3) -> Vec3:
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def _norm(a: Vec3) -> float:
    return math.sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2])


def tri_geometry(
    a: Vec3, b: Vec3, c: Vec3
) -> dict[str, float]:
    """Area, side lengths, perimeter, circumradius and inradius of one triangle.

    Degenerate (zero-area) triangles return ``inf`` for the radii and the ratios that
    need them -- never ``nan``, so a caller can sort and compare.
    """
    sides = (_norm(_sub(b, a)), _norm(_sub(c, b)), _norm(_sub(a, c)))
    area = 0.5 * _norm(_cross(_sub(b, a), _sub(c, a)))
    perimeter = sides[0] + sides[1] + sides[2]
    longest = max(sides)
    if area <= 0.0:
        return {
            "area": 0.0,
            "l_max": longest,
            "perimeter": perimeter,
            "circumradius": _INF,
            "inradius": 0.0,
        }
    # R = a b c / (4A); r = 2A / p. Written with the product form so the relative error
    # is the same order as the area's, not worse.
    circumradius = (sides[0] * sides[1] * sides[2]) / (4.0 * area)
    return {
        "area": area,
        "l_max": longest,
        "perimeter": perimeter,
        "circumradius": circumradius,
        "inradius": 2.0 * area / perimeter if perimeter > 0.0 else 0.0,
    }


def legacy_aspect_from_metrics(area: float, longest: float) -> float:
    """``l_max**2 / (2A)`` -- the single definition behind ``mesh_conditioning``/``validation``.

    The arithmetic is spelled ``longest * longest / (2.0 * area)`` on purpose: that is
    the expression both production modules used before this module existed, so a
    delegation to it is bit-for-bit identical rather than merely close.
    """
    if area <= 0.0:
        return _INF
    return longest * longest / (2.0 * area)


def _aspects(g: Mapping[str, float]) -> dict[str, float]:
    area = g["area"]
    longest = g["l_max"]
    perimeter = g["perimeter"]
    legacy = legacy_aspect_from_metrics(area, longest)
    # alpha = p / (2 sqrt(3) l_max) is defined by edge lengths alone, so a collinear
    # triangle still has a value (its minimum, 1/sqrt(3)); beta needs the circumradius
    # and is infinite when the triangle is flat.
    alpha = perimeter / (2.0 * SQRT3 * longest) if longest > 0.0 else _INF
    if area <= 0.0:
        return {"legacy": legacy, "edge_form": _INF, "verdict": _INF, "alpha": alpha, "beta": _INF}
    edge_form = longest * perimeter / (4.0 * SQRT3 * area)
    verdict = g["circumradius"] / (2.0 * g["inradius"]) if g["inradius"] > 0.0 else _INF
    return {
        "legacy": legacy,
        "edge_form": edge_form,
        "verdict": verdict,
        "alpha": alpha,
        "beta": SQRT3 * g["circumradius"] / longest if longest > 0.0 else _INF,
    }


def aspect_dialects(
    points: Sequence[Sequence[float]], triangles: Iterable[Sequence[int]]
) -> dict[str, tuple[float, ...]]:
    """All three dialect values per triangle, plus the two transfer factors."""
    out: dict[str, list[float]] = {k: [] for k in ("legacy", "edge_form", "verdict", "alpha", "beta")}
    for face in triangles:
        a, b, c = (tuple(float(x) for x in points[i]) for i in face)  # type: ignore[misc]
        values = _aspects(tri_geometry(a, b, c))  # type: ignore[arg-type]
        for key, value in values.items():
            out[key].append(value)
    return {key: tuple(value) for key, value in out.items()}


def transfer_ratio_bounds() -> dict[str, dict[str, float]]:
    """The proved bounds, as numbers a test can assert against."""
    return {
        "alpha_edge_form_over_legacy": {"min": ALPHA_MIN, "max": ALPHA_MAX, "attained": 1.0},
        "beta_verdict_over_edge_form": {"min": BETA_MIN, "max": math.inf, "attained": 1.0},
        "verdict_over_legacy": {"min": COMBINED_MIN, "max": math.inf, "attained": None},
    }


def conservative_limit(limit: float, *, target: str, source: str = "legacy") -> float:
    """The limit in ``target`` dialect that can never pass what ``source`` at ``limit`` rejects.

    Proof, not measurement: every element satisfies ``target <= bound * source`` with the
    proved constants, so ``target <= limit / bound`` implies ``source <= limit``. The
    returned limit is therefore stricter than a naive copy of the number -- by exactly the
    factor the two dialects differ by, in the worst case.
    """
    pairs = {
        # target values: limit / (worst-case bound on target/source)
        ("legacy", "edge_form"): limit * ALPHA_MIN,        # limit / sqrt(3)
        ("edge_form", "legacy"): limit / ALPHA_MIN,        # limit * sqrt(3)
        ("legacy", "verdict"): limit * COMBINED_MIN,       # limit / 2
        ("verdict", "legacy"): limit / COMBINED_MIN,       # limit * 2
        ("edge_form", "verdict"): limit * BETA_MIN,
        # verdict/edge_form is unbounded above (sliver circumradius), so no finite
        # edge_form limit can be shown to be at least as strict; 0.0 is the honest answer.
        ("verdict", "edge_form"): 0.0,
    }
    if source == target:
        return limit
    if (source, target) in pairs:
        return pairs[(source, target)]
    raise ValueError(f"unknown dialect pair {source!r} -> {target!r} ")


def dialect_audit(
    meshes: Mapping[str, tuple[Sequence[Sequence[float]], Sequence[Sequence[int]]]],
    *,
    limit: float = 25.0,
) -> dict[str, object]:
    """Per-mesh dialect spread and the gate flips a copied limit would introduce.

    ``flips`` counts elements that pass an ``edge_form`` (or ``verdict``) gate pinned at
    ``limit`` while failing the ``legacy`` gate at the same number -- the measurable
    consequence of adopting the wrong dialect. ``conservative`` repeats the count with
    :func:`conservative_limit`, which is 0 by construction and therefore a self-check.
    """
    report: dict[str, object] = {
        "limit": limit,
        "conservative_limit": {
            "edge_form": conservative_limit(limit, target="edge_form"),
            "verdict": conservative_limit(limit, target="verdict"),
        },
        "bounds": transfer_ratio_bounds(),
        "meshes": {},
    }
    for name, (points, triangles) in meshes.items():
        values = aspect_dialects(points, triangles)
        legacy = values["legacy"]
        finite = [v for v in values["alpha"] if math.isfinite(v)]
        entry: dict[str, object] = {
            "triangles": len(legacy),
            "max_legacy": max(legacy) if legacy else None,
            "max_edge_form": max(values["edge_form"]) if legacy else None,
            "max_verdict": max(values["verdict"]) if legacy else None,
            "alpha_min": min(finite) if finite else None,
            "alpha_max": max(finite) if finite else None,
            "flips": {},
        }
        for dialect in ("edge_form", "verdict"):
            naive = sum(
                1
                for got, ref in zip(values[dialect], legacy)
                if got <= limit and ref > limit
            )
            safe = conservative_limit(limit, target=dialect)
            guarded = sum(
                1
                for got, ref in zip(values[dialect], legacy)
                if got <= safe and ref > limit
            )
            entry["flips"][dialect] = {  # type: ignore[index]
                "at_copied_limit": naive,
                "at_conservative_limit": guarded,
            }
        report["meshes"][name] = entry  # type: ignore[index]
    return report