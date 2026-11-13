"""Pins the aspect-ratio dialect algebra and its threshold-transfer guarantee.

The properties asserted here are of three kinds, in increasing strength:

* **identity** -- the library's three dialect values equal, to rounding, the numbers the
  three production modules already produce (``mesh_conditioning``, ``validation``,
  ``geomquality.quality``). That is what makes the delegation in ``validation`` a
  no-op on the numbers and closes finding F4 of
  ``docs/mesh-quality-lib-audit-m17404.md``;
* **bounds** -- ``alpha = edge_form/legacy`` stays inside ``[1/sqrt(3), sqrt(3)/2]`` and
  ``beta = verdict/edge_form`` stays ``>= sqrt(3)/2`` on random triangles, with both
  bounds *attained* by the triangles the proofs say attain them (so the bounds are
  tight, not merely safe);
* **consequence** -- a gate in a foreign dialect pinned at the copied limit 25 passes
  elements the legacy gate rejects (the auditor's counterexample is reproduced
  independently here), while :func:`conservative_limit` passes none on any mesh.

Every number below is computed, none is transcribed from a document.
"""
from __future__ import annotations

import math
import random

import pytest

from tensorcad.geometry import Mesh, MeshBuilder
from tensorcad.geomquality.quality import triangle_quality
from tensorcad.hull import HullParameters, build_mesh
from tensorcad.mesh_conditioning import _aspect_ratio
from tensorcad.meshquality import tri_metrics
from tensorcad.meshquality.dialects import (
    ALPHA_MAX,
    ALPHA_MIN,
    ASPECT_DIALECTS,
    BETA_MIN,
    COMBINED_MIN,
    aspect_dialects,
    conservative_limit,
    dialect_audit,
    legacy_aspect_from_metrics,
    transfer_ratio_bounds,
    tri_geometry,
)
from tensorcad.validation import measure_quality

SQRT3 = math.sqrt(3.0)


# --------------------------------------------------------------------------------------
# meshes: a unit box surface, a hull surface, and seeded random triangle soups
# --------------------------------------------------------------------------------------
def box_mesh() -> Mesh:
    builder = MeshBuilder()
    for point in [
        (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0), (0.0, 1.0, 1.0),
    ]:
        builder.add_vertex(point)
    for a, b, c, d in ((0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)):
        builder.add_triangle(a, b, c)
        builder.add_triangle(a, c, d)
    return builder.build()


def hull_mesh() -> Mesh:
    mesh, _ = build_mesh(HullParameters(length=10.0, beam=2.0, draft=0.8, nx=21, nz=9))
    return mesh


def random_triangles(count: int = 400, seed: int = 20261113) -> tuple[list[tuple[float, float, float]], list[tuple[int, int, int]]]:
    """Random triangles spanning well-shaped to nearly collinear (fixed seed, reproducible)."""
    rng = random.Random(seed)
    points: list[tuple[float, float, float]] = []
    triangles: list[tuple[int, int, int]] = []
    for _ in range(count):
        base = (rng.uniform(-2.0, 2.0), rng.uniform(-2.0, 2.0), rng.uniform(-2.0, 2.0))
        scale = 10.0 ** rng.uniform(-3.0, 3.0)
        p1 = tuple(base[i] + scale * rng.uniform(-1.0, 1.0) for i in range(3))
        p2 = tuple(base[i] + scale * rng.uniform(-1.0, 1.0) for i in range(3))
        start = len(points)
        points.extend((base, p1, p2))  # type: ignore[arg-type]
        triangles.append((start, start + 1, start + 2))
    return points, triangles


def equilateral(edge: float = 1.0) -> tuple[tuple[float, float, float], ...]:
    return ((0.0, 0.0, 0.0), (edge, 0.0, 0.0), (0.5 * edge, SQRT3 / 2.0 * edge, 0.0))


def right_isosceles() -> tuple[tuple[float, float, float], ...]:
    """45-45-90: the triangle on which all three dialects disagree."""
    return ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0))


# --------------------------------------------------------------------------------------
# 1. identity with the three production modules
# --------------------------------------------------------------------------------------
def test_legacy_helper_is_bit_identical_to_the_shipped_expression():
    """``mesh_conditioning`` and ``legacy_aspect_from_metrics`` must agree exactly."""
    mesh = hull_mesh()
    points = list(mesh.vertices)
    for face in mesh.faces:
        a, b, c = (points[i] for i in face)
        # mesh_conditioning's helper takes twice the area and a normal it does not use
        cross = (
            (b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1]),
            (b[2] - a[2]) * (c[0] - a[0]) - (b[0] - a[0]) * (c[2] - a[2]),
            (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]),
        )
        area2 = math.sqrt(sum(component * component for component in cross))
        shipped = _aspect_ratio([a, b, c], (0.0, 0.0, 1.0), area2)
        mine = legacy_aspect_from_metrics(0.5 * area2, max(
            math.dist(a, b), math.dist(b, c), math.dist(c, a)
        ))
        assert mine == shipped or abs(mine - shipped) <= 1e-15 * max(1.0, abs(shipped))


def test_validation_measure_quality_unchanged_by_the_delegation():
    """``validation`` reports the same worst aspect ratio through the shared definition."""
    mesh = hull_mesh()
    measured = measure_quality(mesh)
    points = list(mesh.vertices)
    worst = 0.0
    for face in mesh.faces:
        values = aspect_dialects(points, [tuple(face)])["legacy"][0]
        worst = max(worst, values)
    assert measured.max_aspect_ratio == pytest.approx(worst, rel=1e-15, abs=0.0)


def test_edge_form_matches_geomquality_triangle_quality():
    """The canonical ``geomquality`` number is the ``edge_form`` dialect, to rounding."""
    points, triangles = random_triangles(120)
    mine = aspect_dialects(points, triangles)["edge_form"]
    for (i, j, k), value in zip(triangles, mine):
        theirs = triangle_quality(points[i], points[j], points[k])
        assert value == pytest.approx(theirs.aspect_ratio, rel=1e-12, abs=1e-12)


def test_verdict_matches_meshquality_tri_metrics():
    """The newer ``meshquality`` number is the ``verdict`` dialect, to rounding."""
    points, triangles = random_triangles(120)
    mine = aspect_dialects(points, triangles)["verdict"]
    theirs = tri_metrics(points, triangles)["aspect_ratio"]
    for a, b in zip(mine, theirs):
        assert a == pytest.approx(b, rel=1e-12, abs=1e-12)


def test_the_two_canonical_modules_disagree_on_the_same_element():
    """Regression pin for the defect this module was written to fix.

    On a 45-45-90 triangle the two modules that each call themselves canonical report
    different aspect ratios; the factor is ``beta = sqrt(3) R / l_max = sqrt(3)/2`` here,
    the minimum the bound allows.
    """
    points = list(right_isosceles())
    values = aspect_dialects(points, [(0, 1, 2)])
    canonical_new = tri_metrics(points, [(0, 1, 2)])["aspect_ratio"][0]
    canonical_old = triangle_quality(*points).aspect_ratio  # type: ignore[arg-type]
    assert canonical_new == pytest.approx(values["verdict"][0], rel=1e-12)
    assert canonical_old == pytest.approx(values["edge_form"][0], rel=1e-12)
    assert canonical_new / canonical_old == pytest.approx(BETA_MIN, rel=1e-12)
    assert canonical_old != pytest.approx(canonical_new, rel=1e-6)


# --------------------------------------------------------------------------------------
# 2. bounds, and that they are tight
# --------------------------------------------------------------------------------------
def test_only_two_dialects_are_normalised_to_one_at_the_ideal_element():
    """``edge_form`` and ``verdict`` read 1 on the equilateral; ``legacy`` reads 2/sqrt(3).

    This is the reason "ideal value 1" cannot be used to tell dialects apart, and the
    reason the shipped limit 25 is not on the same scale as a normalised metric.
    """
    values = aspect_dialects(list(equilateral()), [(0, 1, 2)])
    for dialect in ("edge_form", "verdict"):
        assert values[dialect][0] == pytest.approx(1.0, abs=1e-12)
    assert values["legacy"][0] == pytest.approx(2.0 / SQRT3, rel=1e-12)
    assert ASPECT_DIALECTS["legacy"]["ideal"] == pytest.approx(2.0 / SQRT3, rel=1e-12)
    assert values["alpha"][0] == pytest.approx(ALPHA_MAX, rel=1e-12)
    assert values["beta"][0] == pytest.approx(1.0, rel=1e-12)


def test_alpha_stays_inside_the_proved_band_and_attains_both_ends():
    points, triangles = random_triangles()
    alpha = [a for a in aspect_dialects(points, triangles)["alpha"] if math.isfinite(a)]
    assert alpha and min(alpha) >= ALPHA_MIN - 1e-15 and max(alpha) <= ALPHA_MAX + 1e-15
    # attained at degeneracy: a collinear triangle has l_max = p/2, so alpha = 1/sqrt(3)
    collinear = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (3.0, 0.0, 0.0)]
    assert aspect_dialects(collinear, [(0, 1, 2)])["alpha"][0] == pytest.approx(ALPHA_MIN, rel=1e-12)
    assert aspect_dialects(list(equilateral()), [(0, 1, 2)])["alpha"][0] == pytest.approx(ALPHA_MAX, rel=1e-12)


def test_beta_is_bounded_below_and_the_minimum_is_a_right_triangle():
    points, triangles = random_triangles()
    beta = [b for b in aspect_dialects(points, triangles)["beta"] if math.isfinite(b)]
    assert beta and min(beta) >= BETA_MIN - 1e-12
    # l_max = 2R exactly when the longest edge is a diameter, i.e. a right triangle
    assert aspect_dialects(list(right_isosceles()), [(0, 1, 2)])["beta"][0] == pytest.approx(
        BETA_MIN, rel=1e-12
    )


def test_verdict_over_legacy_never_below_the_combined_bound():
    points, triangles = random_triangles()
    values = aspect_dialects(points, triangles)
    ratios = [
        v / l
        for v, l in zip(values["verdict"], values["legacy"])
        if math.isfinite(v) and math.isfinite(l) and l > 0.0
    ]
    assert ratios and min(ratios) >= COMBINED_MIN - 1e-12


def test_transfer_bounds_are_reported_as_numbers():
    bounds = transfer_ratio_bounds()
    assert bounds["alpha_edge_form_over_legacy"]["min"] == pytest.approx(1.0 / SQRT3)
    assert bounds["alpha_edge_form_over_legacy"]["max"] == pytest.approx(SQRT3 / 2.0)
    assert bounds["beta_verdict_over_edge_form"]["min"] == pytest.approx(SQRT3 / 2.0)
    assert bounds["verdict_over_legacy"]["min"] == pytest.approx(0.5)


# --------------------------------------------------------------------------------------
# 3. the consequence: copied limits relax the gate, conservative limits do not
# --------------------------------------------------------------------------------------
def test_auditor_counterexample_reproduced_independently():
    """The exact element quoted in finding F1, recomputed here from scratch."""
    points = [(0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (1.0, 0.0741, 0.0)]
    values = aspect_dialects(points, [(0, 1, 2)])
    assert values["legacy"][0] == pytest.approx(26.9906, abs=5e-4)
    assert values["edge_form"][0] == pytest.approx(15.6044, abs=5e-4)
    assert values["legacy"][0] > 25.0 > values["edge_form"][0]  # the flip


def test_copied_limit_flips_elements_and_conservative_limit_flips_none():
    legacy_limit = 25.0
    points = list(hull_mesh().vertices)
    faces = [tuple(face) for face in hull_mesh().faces]
    values = aspect_dialects(points, faces)

    def flips(dialect: str, limit: float) -> int:
        return sum(
            1
            for got, ref in zip(values[dialect], values["legacy"])
            if got <= limit and ref > legacy_limit
        )

    # at least one mesh in this set must contain an element where a copied limit is
    # unsafe, otherwise the test would pass vacuously: check the random soup too
    random_points, random_faces = random_triangles(400)
    random_values = aspect_dialects(random_points, random_faces)
    random_flips = sum(
        1
        for got, ref in zip(random_values["edge_form"], random_values["legacy"])
        if got <= legacy_limit and ref > legacy_limit
    )
    assert random_flips > 0
    assert flips("edge_form", conservative_limit(legacy_limit, target="edge_form")) == 0
    assert flips("verdict", conservative_limit(legacy_limit, target="verdict")) == 0


def test_conservative_limit_direction_and_values():
    assert conservative_limit(25.0, target="edge_form") == pytest.approx(25.0 / SQRT3)
    assert conservative_limit(25.0, target="verdict") == pytest.approx(12.5)
    assert conservative_limit(25.0, target="legacy", source="edge_form") == pytest.approx(25.0 * SQRT3)
    assert conservative_limit(25.0, target="verdict", source="edge_form") == pytest.approx(25.0 * BETA_MIN)
    # no finite edge_form limit is provably stricter than a verdict gate
    assert conservative_limit(25.0, target="edge_form", source="verdict") == 0.0
    assert conservative_limit(25.0, target="legacy", source="legacy") == 25.0
    with pytest.raises(ValueError):
        conservative_limit(25.0, target="nonsense")


@pytest.mark.parametrize("dialect", ["edge_form", "verdict"])
def test_conservative_gate_is_a_subset_of_the_legacy_gate(dialect: str):
    """The proof, checked on every element of two meshes and on random triangles."""
    meshes = {"box": (list(box_mesh().vertices), [tuple(f) for f in box_mesh().faces])}
    hull = hull_mesh()
    meshes["hull"] = (list(hull.vertices), [tuple(f) for f in hull.faces])
    meshes["random"] = random_triangles(300)
    limit = conservative_limit(25.0, target=dialect)
    for name, (points, faces) in meshes.items():
        values = aspect_dialects(points, faces)
        for got, ref in zip(values[dialect], values["legacy"]):
            if got <= limit:
                assert ref <= 25.0 + 1e-12, name


# --------------------------------------------------------------------------------------
# 4. degenerate input, scale invariance, round trip
# --------------------------------------------------------------------------------------
def test_degenerate_triangles_are_infinite_not_nan():
    points = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (3.0, 0.0, 0.0)]
    values = aspect_dialects(points, [(0, 1, 2)])
    for dialect in ("legacy", "edge_form", "verdict"):
        assert values[dialect][0] == math.inf
        assert not math.isnan(values[dialect][0])
    geometry = tri_geometry(*points)  # type: ignore[arg-type]
    assert geometry["area"] == 0.0


def test_dialects_are_scale_invariant():
    points, triangles = random_triangles(60)
    base = aspect_dialects(points, triangles)
    scaled_points = [tuple(1e3 * c for c in p) for p in points]
    scaled = aspect_dialects(scaled_points, triangles)
    for dialect in ("legacy", "edge_form", "verdict"):
        for a, b in zip(base[dialect], scaled[dialect]):
            if math.isfinite(a) and math.isfinite(b):
                assert a == pytest.approx(b, rel=1e-12)


def test_dialect_audit_is_self_consistent():
    hull = hull_mesh()
    audit = dialect_audit({"hull": (list(hull.vertices), [tuple(f) for f in hull.faces])})
    entry = audit["meshes"]["hull"]
    assert entry["triangles"] == len(hull.faces)
    for dialect in ("edge_form", "verdict"):
        assert entry["flips"][dialect]["at_conservative_limit"] == 0
    assert audit["conservative_limit"]["edge_form"] == pytest.approx(25.0 / SQRT3)
    assert set(audit["meshes"]["hull"]["flips"]) == {"edge_form", "verdict"}


def test_dialect_metadata_names_all_three_and_their_producers():
    assert set(ASPECT_DIALECTS) == {"legacy", "edge_form", "verdict"}
    for name, entry in ASPECT_DIALECTS.items():
        assert entry["ideal"] == pytest.approx(1.0 if name != "legacy" else 2.0 / SQRT3, rel=1e-12)
        assert entry["produced_by"]
        assert entry["definition"]