"""Tests for fore/aft LCB control in :mod:`tensorcad.cb_sections`.

The claim under test is the same shape as the coefficient claim in
``tests/test_cb_sections.py``, applied to a fourth quantity: a hull built with
``lcb_offset = o`` must carry its centre of buoyancy at ``o`` (fraction of L, forward
positive) while ``Cb``, ``Cm`` and ``Cwp`` stay exact. Every assertion is made against
an integration of the emitted mesh -- the divergence-theorem volume and centroid --
never against an intermediate of the construction.

The second, equally important claim is that the feature is *opt-in*: ``lcb_offset = 0``
has to reproduce the shipped symmetric family bit for bit, so this work cannot drift
the geometry every earlier round of the rebuild campaign certifies.
"""
from __future__ import annotations

import pytest

from tensorcad import mesh as meshlib
from tensorcad.geometry import Mesh
from tensorcad.cb_sections import (
    LCB_OFFSET_LIMIT,
    CalibratedHullParameters,
    LcbRequestUnreachable,
    build_calibrated_hull,
    longitudinal_offset,
    realised_coefficients,
    realised_lcb_fraction,
    section_side_signs,
    solve_lcb_tilt,
)

#: (L, B, T, Cb, Cm, Cwp, published LCB) -- the two benchmark ships whose LCB is in the
#: repository's own fleet manifest (``tensorcad/fleet/data/fleet_manifest.json``,
#: retrieval=offline_recollection) plus three interpolation points.
BENCHMARKS = [
    ("KCS", 230.0, 32.2, 10.8, 0.651, 0.971, 0.7674, -0.0148),
    ("KVLCC2", 320.0, 58.0, 20.8, 0.8098, 0.9956, 0.8700, 0.0348),
]

CASES = [
    (200.0, 30.0, 10.0, 0.55, 0.95, 0.75),
    (200.0, 30.0, 10.0, 0.65, 0.98, 0.80),
    (280.0, 40.0, 15.0, 0.72, 0.99, 0.84),
    (320.0, 58.0, 20.8, 0.8098, 0.9956, 0.87),
]

TARGETS = [-0.03, -0.0148, 0.0, 0.02, 0.0348]

#: Tolerance on the realised LCB, as a fraction of L. Round-off on the symmetric
#: family; the skewed family is a bisection, and on fine hulls the profile exponent of
#: the vanishing end sections clamps at its bound, which costs a few 1e-5 L.
LCB_TOLERANCE = 2.0e-4


def _euler(mesh) -> int:
    edges = set()
    for face in mesh.faces:
        for i in range(3):
            a, b = face[i], face[(i + 1) % 3]
            edges.add((min(a, b), max(a, b)))
    return len(mesh.vertices) - len(edges) + len(mesh.faces)


@pytest.mark.parametrize("length,beam,draft,cb,cm,cwp", CASES)
@pytest.mark.parametrize("target", TARGETS)
def test_requested_lcb_is_realised_on_the_mesh(length, beam, draft, cb, cm, cwp, target):
    hull = build_calibrated_hull(
        CalibratedHullParameters(
            length, beam, draft, cb, cm, cwp, nx=81, nz=21, lcb_offset=target
        )
    )
    measured = realised_coefficients(hull.mesh, length=length, beam=beam, draft=draft)
    assert measured["lcb_fraction"] == pytest.approx(target, abs=LCB_TOLERANCE)


@pytest.mark.parametrize("length,beam,draft,cb,cm,cwp", CASES)
@pytest.mark.parametrize("target", TARGETS)
def test_coefficients_stay_exact_under_tilt(length, beam, draft, cb, cm, cwp, target):
    """Moving the centroid must not cost Cb, Cm or Cwp -- that is the engineering point."""
    hull = build_calibrated_hull(
        CalibratedHullParameters(
            length, beam, draft, cb, cm, cwp, nx=81, nz=21, lcb_offset=target
        )
    )
    coefficients = realised_coefficients(
        hull.mesh, length=length, beam=beam, draft=draft
    )["coefficients"]
    assert abs(coefficients["Cb"] - cb) / cb < 2e-3
    assert abs(coefficients["Cm"] - cm) / cm < 1e-6
    assert abs(coefficients["Cwp"] - cwp) / cwp < 1e-6


@pytest.mark.parametrize("cwp", [0.75, 0.80, 0.84])
def test_only_the_centroid_moves_not_the_midship_section(cwp):
    """``Cm`` is a *maximum* sectional area; a tilted hull must still peak at midship."""
    hull = build_calibrated_hull(
        CalibratedHullParameters(200.0, 30.0, 10.0, 0.65, 0.98, cwp, nx=81, nz=21,
                                 lcb_offset=0.03)
    )
    assert hull.station_weights[40] == pytest.approx(max(hull.station_weights), rel=1e-12)
    assert realised_coefficients(
        hull.mesh, length=200.0, beam=30.0, draft=10.0
    )["max_section_area_m2"] == pytest.approx(0.98 * 30.0 * 10.0, rel=1e-9)


def test_zero_offset_reproduces_the_symmetric_family_bit_for_bit():
    """The default has to be the shipped round-3 surface, not a re-derivation of it."""
    base = CalibratedHullParameters(200.0, 30.0, 10.0, 0.65, 0.98, 0.80, nx=41, nz=13)
    explicit = CalibratedHullParameters(
        200.0, 30.0, 10.0, 0.65, 0.98, 0.80, nx=41, nz=13, lcb_offset=0.0
    )
    first = build_calibrated_hull(base)
    second = build_calibrated_hull(explicit)
    assert first.mesh.vertices == second.mesh.vertices
    assert first.mesh.faces == second.mesh.faces
    assert first.lcb_tilt == 0.0
    assert first.fullness_exponent == second.fullness_exponent


def test_symmetric_hull_has_zero_lcb_but_the_moment_route_is_a_real_measurement():
    """On the symmetric family the answer is known a priori -- check the instrument."""
    hull = build_calibrated_hull(
        CalibratedHullParameters(200.0, 30.0, 10.0, 0.72, 0.99, 0.84, nx=81, nz=21)
    )
    assert realised_lcb_fraction(hull.mesh) == pytest.approx(0.0, abs=1e-12)
    # A translation of the same solid must not change the answer: the offset is
    # referred to the mesh's own mid-length, not to wherever x = 0 happens to be.
    source = hull.mesh
    shifted = Mesh(
        vertices=tuple((x + 150.0, y - 20.0, z + 3.0) for x, y, z in source.vertices),
        faces=source.faces,
        units=source.units,
        coordinate_frame=source.coordinate_frame,
    )
    assert realised_lcb_fraction(shifted) == pytest.approx(0.0, abs=1e-12)
    # ... and it must not be silently zero either: an asymmetric copy has to read back.
    tilted = build_calibrated_hull(
        CalibratedHullParameters(200.0, 30.0, 10.0, 0.72, 0.99, 0.84, nx=41, nz=13,
                                 lcb_offset=0.025)
    )
    assert abs(realised_lcb_fraction(tilted.mesh)) > 1.0e-3


def test_mesh_moment_and_section_moment_agree():
    """Two independent routes to the same number: the triangles and the law that built them."""
    for target in (-0.03, -0.01, 0.01, 0.0348):
        hull = build_calibrated_hull(
            CalibratedHullParameters(280.0, 40.0, 15.0, 0.72, 0.99, 0.84, nx=81, nz=21,
                                     lcb_offset=target)
        )
        mesh_offset = realised_lcb_fraction(hull.mesh)
        section_offset = longitudinal_offset(hull.station_weights)
        assert mesh_offset == pytest.approx(section_offset, abs=1e-6)


def test_lcb_is_monotone_in_the_request():
    offsets = []
    for target in (-0.03, -0.015, 0.0, 0.015, 0.03):
        hull = build_calibrated_hull(
            CalibratedHullParameters(200.0, 30.0, 10.0, 0.65, 0.98, 0.80, nx=41, nz=13,
                                     lcb_offset=target)
        )
        offsets.append(realised_lcb_fraction(hull.mesh))
    assert offsets == sorted(offsets)
    assert all(later > earlier for earlier, later in zip(offsets, offsets[1:]))


@pytest.mark.parametrize("target", [-0.03, 0.0348])
def test_tilted_hull_is_still_a_closed_solid(target):
    hull = build_calibrated_hull(
        CalibratedHullParameters(320.0, 58.0, 20.8, 0.8098, 0.9956, 0.87, nx=81, nz=21,
                                 lcb_offset=target)
    )
    mesh = hull.mesh
    assert meshlib.is_watertight(mesh)
    assert meshlib.is_edge_manifold(mesh)
    assert _euler(mesh) == 2
    assert meshlib.total_signed_volume(mesh) > 0.0
    assert meshlib.zero_area_faces(mesh) == ()


@pytest.mark.parametrize("name,length,beam,draft,cb,cm,cwp,published", BENCHMARKS)
def test_published_benchmark_lcb_is_hit(name, length, beam, draft, cb, cm, cwp, published):
    """The point of the feature: a real ship's published LCB, not a synthetic one.

    ``name`` is carried for the failure message; the reference value is the fleet
    manifest's ``lcb_fraction`` for that hull.
    """
    hull = build_calibrated_hull(
        CalibratedHullParameters(length, beam, draft, cb, cm, cwp, nx=81, nz=21,
                                 lcb_offset=published)
    )
    measured = realised_coefficients(hull.mesh, length=length, beam=beam, draft=draft)
    assert measured["lcb_fraction"] == pytest.approx(published, abs=LCB_TOLERANCE), name
    assert abs(measured["coefficients"]["Cb"] - cb) / cb < 2e-3, name


def test_parallel_body_and_keel_rise_survive_the_tilt():
    """The two optional hull-form features must compose with LCB control."""
    hull = build_calibrated_hull(
        CalibratedHullParameters(320.0, 58.0, 20.8, 0.65, 0.98, 0.80, nx=81, nz=21,
                                 parallel_body_fraction=0.3, end_draft_ratio=0.5,
                                 lcb_offset=0.0348)
    )
    measured = realised_coefficients(hull.mesh, length=320.0, beam=58.0, draft=20.8)
    assert measured["lcb_fraction"] == pytest.approx(0.0348, abs=LCB_TOLERANCE)
    assert abs(measured["coefficients"]["Cb"] - 0.65) / 0.65 < 2e-3
    assert meshlib.is_watertight(hull.mesh)


@pytest.mark.parametrize("offset", [LCB_OFFSET_LIMIT, -LCB_OFFSET_LIMIT, 0.4, float("nan")])
def test_out_of_range_offset_is_refused_before_building(offset):
    params = CalibratedHullParameters(200.0, 30.0, 10.0, 0.65, 0.98, 0.80,
                                      lcb_offset=offset)
    problems = params.problems()
    assert any("lcb_offset" in text for text in problems)
    with pytest.raises(ValueError):
        build_calibrated_hull(params)


def test_unreachable_offset_raises_instead_of_returning_a_wrong_hull():
    """A full hull at the edge of the feasible band cannot be moved 4.9 % of L.

    Silently returning a hull whose LCB is not the requested one is exactly the
    failure mode this module exists to remove, so the refusal is part of the contract.
    """
    params = CalibratedHullParameters(200.0, 30.0, 10.0, 0.85, 0.995, 0.90, nx=41,
                                      nz=13, lcb_offset=0.049)
    with pytest.raises(LcbRequestUnreachable):
        build_calibrated_hull(params)


def test_solve_lcb_tilt_returns_zero_without_searching():
    """The short-circuit is what makes the default bit-for-bit; pin it."""
    taper = CalibratedHullParameters(200.0, 30.0, 10.0, 0.65, 0.98, 0.80).taper()
    keel = (1.0,) * len(taper)
    assert solve_lcb_tilt(0.0, taper, keel, 0.80, 0.65 / 0.98) == 0.0


def test_side_signs_are_symmetric_and_centred():
    """Ties the whole construction to one convention: +1 aft, -1 forward, 0 at midship."""
    assert section_side_signs(5) == (1.0, 1.0, 0.0, -1.0, -1.0)
    assert section_side_signs(6) == (1.0, 1.0, 1.0, -1.0, -1.0, -1.0)
    for nx in (3, 5, 21, 41, 81):
        signs = section_side_signs(nx)
        assert sum(signs) == 0.0
        assert signs == tuple(-sign for sign in reversed(signs))