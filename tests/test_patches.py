"""Tests for :mod:`tensorcad.patches` -- the multi-patch stitcher.

The acceptance anchors are analytic and independent of the stitcher:

* a **box** of six bilinear patches. A bilinear patch is exact under tessellation, so
  the closed shell's volume and area are ``a*b*c`` and ``2(ab+bc+ca)`` to round-off and
  any error would be a seam error, not a discretisation error;
* a **torus** of 16 rational biquadratic patches against ``V = 2 pi^2 R r^2`` and
  ``A = 4 pi^2 R r`` (second-order convergence expected, the surface is curved);
* a **double cone** of 8 patches with two collinear apex poles against
  ``V = 2 pi R^2 H / 3`` and ``A = 2 pi R sqrt(R^2 + H^2)``, which is the pole-collapse
  case: the apex sides are claimed by no seam and must weld to one vertex each.

Every claim below is a measured number from a run of this file.
"""
from __future__ import annotations

import math

import pytest

from tensorcad import mesh as meshlib
from tensorcad.nurbs import NURBSCurve, circular_arc
from tensorcad.patches import (
    PatchNetwork,
    PatchNetworkError,
    Seam,
    Stitched,
    bilinear_patch,
    derive_seams,
    revolve_patch,
    seam_network,
    stitch,
)

BOX = (2.0, 3.0, 4.0)


# ---------------------------------------------------------------------------
# fixtures: three networks, built from analytic geometry only
# ---------------------------------------------------------------------------
def box_network(a: float = BOX[0], b: float = BOX[1], c: float = BOX[2]) -> PatchNetwork:
    """Six bilinear patches, wound so that every normal points out of the box."""
    faces = (
        ((0.0, 0.0, 0.0), (0.0, b, 0.0), (a, 0.0, 0.0)),        # z = 0, normal -z
        ((0.0, 0.0, c), (a, 0.0, 0.0), (0.0, b, 0.0)),          # z = c, normal +z
        ((0.0, 0.0, 0.0), (a, 0.0, 0.0), (0.0, 0.0, c)),        # y = 0, normal -y
        ((0.0, b, 0.0), (0.0, 0.0, c), (a, 0.0, 0.0)),          # y = b, normal +y
        ((0.0, 0.0, 0.0), (0.0, 0.0, c), (0.0, b, 0.0)),        # x = 0, normal -x
        ((a, 0.0, 0.0), (0.0, b, 0.0), (0.0, 0.0, c)),          # x = a, normal +x
    )
    surfaces = tuple(bilinear_patch(origin, du, dv) for origin, du, dv in faces)
    return seam_network(surfaces, derive_seams(surfaces))


def meridian_arc(major: float, minor: float, start: float, sweep: float) -> NURBSCurve:
    """A ``(rho, z)`` arc of the torus tube: the tube circle centred at ``rho = major``."""
    arc = circular_arc(minor, sweep, start_angle=start)
    return NURBSCurve(
        degree=arc.degree,
        knots=arc.knots,
        control_points=tuple((major + point[0], point[1], 0.0) for point in arc.control_points),
        weights=arc.weights,
    )


def torus_network(major: float = 2.0, minor: float = 0.5) -> PatchNetwork:
    """Sixteen patches: four exact 90-degree tube arcs revolved in four 90-degree arcs."""
    profiles = [meridian_arc(major, minor, -0.5 * math.pi + k * 0.5 * math.pi, 0.5 * math.pi)
                for k in range(4)]
    surfaces = tuple(
        revolve_patch(profiles[tube], quadrant * 0.5 * math.pi, 0.5 * math.pi)
        for quadrant in range(4)
        for tube in range(4)
    )
    return seam_network(surfaces, derive_seams(surfaces))


def double_cone_network(radius: float = 1.0, half_height: float = 1.0) -> PatchNetwork:
    """Eight patches: two cones apex-to-apex, sharing their base circle at ``z = 0``."""
    up = NURBSCurve(
        degree=1,
        knots=(0.0, 0.0, 1.0, 1.0),
        control_points=((radius, 0.0, 0.0), (0.0, half_height, 0.0)),
    )
    down = NURBSCurve(
        degree=1,
        knots=(0.0, 0.0, 1.0, 1.0),
        control_points=((0.0, -half_height, 0.0), (radius, 0.0, 0.0)),
    )
    surfaces = tuple(
        revolve_patch(profile, quadrant * 0.5 * math.pi, 0.5 * math.pi)
        for profile in (up, down)
        for quadrant in range(4)
    )
    return seam_network(surfaces, derive_seams(surfaces))


def torus_volume(major: float = 2.0, minor: float = 0.5) -> float:
    return 2.0 * math.pi**2 * major * minor**2


def torus_area(major: float = 2.0, minor: float = 0.5) -> float:
    return 4.0 * math.pi**2 * major * minor


# ---------------------------------------------------------------------------
# the box: exact volume, exact area, G0 seams that measure 90 degrees
# ---------------------------------------------------------------------------
class TestBox:
    def test_the_derived_network_claims_every_side_of_every_patch(self):
        network = box_network()
        assert len(network.surfaces) == 6
        assert len(network.seams) == 12
        assert network.unclaimed_sides() == ()

    def test_the_closed_box_has_an_exact_volume_and_area(self):
        a, b, c = BOX
        stitched = stitch(box_network(), nu=4, nv=4)
        report = stitched.solid_report(
            reference_volume=a * b * c, reference_area=2 * (a * b + b * c + c * a)
        )
        assert report["watertight"] is True
        assert report["edge_manifold"] is True
        assert report["boundary_edges"] == 0
        assert report["non_manifold_edges"] == 0
        assert report["components"] == 1
        assert report["zero_area_faces"] == 0
        # A bilinear patch is exact under tessellation: the volume is exact, not merely
        # converged, so the tolerance is round-off rather than a mesh error.
        assert report["signed_volume_m3"] == pytest.approx(a * b * c, rel=1e-15)
        assert report["surface_area_m2"] == pytest.approx(
            2 * (a * b + b * c + c * a), rel=1e-15
        )
        assert report["volume_relative_error"] == pytest.approx(0.0, abs=1e-15)

    def test_every_box_seam_is_positionally_closed_to_the_last_bit(self):
        stitched = stitch(box_network(), nu=6, nv=6)
        assert len(stitched.seams) == 12
        assert all(record.gap_max == 0.0 for record in stitched.seams)
        assert all(record.reversed is False for record in stitched.seams)

    def test_a_box_seam_is_G0_only_and_the_report_says_ninety_degrees(self):
        stitched = stitch(box_network(), nu=6, nv=6)
        angles = [record.angle_max_deg for record in stitched.seams]
        assert all(angle == pytest.approx(90.0, abs=1e-9) for angle in angles)
        assert stitched.seam_report()["angle_max_deg"] == pytest.approx(90.0, abs=1e-9)

    def test_welding_is_index_level_and_counted(self):
        # 6 patches x (nu+1)(nv+1) slots, minus the welded classes.
        stitched = stitch(box_network(), nu=4, nv=4)
        slots = 6 * 5 * 5
        assert stitched.welded_nodes > 0
        assert len(stitched.mesh.vertices) == slots - stitched.welded_nodes


# ---------------------------------------------------------------------------
# the torus: a curved network, tangent seams, second-order convergence
# ---------------------------------------------------------------------------
class TestTorus:
    def test_the_network_is_closed_with_thirty_two_seams_and_no_open_side(self):
        network = torus_network()
        stitched = stitch(network, nu=4, nv=4)
        assert len(network.surfaces) == 16
        assert len(network.seams) == 32
        assert stitched.open_sides == ()
        assert stitched.pole_sides == ()
        assert stitched.mesh and meshlib.is_watertight(stitched.mesh)

    def test_seams_are_closed_and_tangent(self):
        stitched = stitch(torus_network(), nu=8, nv=8)
        assert max(record.gap_max for record in stitched.seams) < 1e-12
        # A torus is smooth, so a seam that only welds (and does not change the patch
        # data) must already be G1: the dihedral angle is round-off, not a dihedral.
        assert max(record.angle_max_deg for record in stitched.seams) < 1e-4
        assert all(record.inset_normals == 0 for record in stitched.seams)

    def test_volume_and_area_converge_at_second_order_to_the_analytic_torus(self):
        errors = []
        for samples in (8, 16, 32):
            report = stitch(torus_network(), nu=samples, nv=samples).solid_report(
                reference_volume=torus_volume(), reference_area=torus_area()
            )
            assert report["watertight"] is True
            errors.append(abs(report["volume_relative_error"]))
        assert errors[0] < 2e-2
        assert errors[-1] < 1.5e-3
        for coarse, fine in zip(errors, errors[1:]):
            # Halving the spacing must cut the error by ~4 (second order); the band is
            # wide enough for the change of regime between the two parameter
            # directions but narrow enough to reject first order (~2).
            assert 3.2 < coarse / fine < 5.0

    def test_the_surface_area_is_below_the_analytic_value_under_refinement(self):
        report = stitch(torus_network(), nu=24, nv=24).solid_report(
            reference_area=torus_area()
        )
        assert -2e-3 < report["area_relative_error"] < 0.0


# ---------------------------------------------------------------------------
# the double cone: pole collapse
# ---------------------------------------------------------------------------
class TestPoles:
    def test_both_apexes_are_unclaimed_sides_that_collapse(self):
        network = double_cone_network()
        stitched = stitch(network, nu=4, nv=4)
        assert len(network.seams) == 12          # 4 upper + 4 lower meridian seams, 4 base
        assert len(stitched.pole_sides) == 8     # the four apexes of each cone
        assert stitched.open_sides == ()
        assert meshlib.is_watertight(stitched.mesh)
        assert meshlib.is_edge_manifold(stitched.mesh)

    def test_the_closed_double_cone_matches_its_analytic_volume_and_area(self):
        radius, height = 1.0, 1.0
        errors = []
        for samples in (8, 16, 32):
            report = stitch(double_cone_network(radius, height), nu=samples, nv=samples).solid_report(
                reference_volume=2.0 * math.pi * radius**2 * height / 3.0,
                reference_area=2.0 * math.pi * radius * math.hypot(radius, height),
            )
            errors.append(abs(report["volume_relative_error"]))
            assert report["watertight"] is True
        assert errors[-1] < 1e-3
        for coarse, fine in zip(errors, errors[1:]):
            assert 3.2 < coarse / fine < 5.0

    def test_a_node_next_to_a_pole_reports_an_inset_normal_instead_of_failing(self):
        # Each meridian seam of a cone ends at the apex, where the Jacobian vanishes.
        stitched = stitch(double_cone_network(), nu=4, nv=4)
        assert sum(record.inset_normals for record in stitched.seams) == 8
        apex_seams = [record for record in stitched.seams if record.inset_normals]
        assert all(record.gap_max < 1e-12 for record in apex_seams)

    def test_claiming_a_pole_side_as_a_seam_is_refused(self):
        radius, height = 1.0, 1.0
        up = NURBSCurve(
            degree=1,
            knots=(0.0, 0.0, 1.0, 1.0),
            control_points=((radius, 0.0, 0.0), (0.0, height, 0.0)),
        )
        surfaces = tuple(
            revolve_patch(up, quadrant * 0.5 * math.pi, 0.5 * math.pi) for quadrant in range(4)
        )
        network = seam_network(surfaces, (Seam((0, "v1"), (1, "v1")),))
        with pytest.raises(PatchNetworkError, match="pole"):
            stitch(network, nu=4, nv=4)


# ---------------------------------------------------------------------------
# direction, openness, and the negative controls that must stay honest
# ---------------------------------------------------------------------------
class TestDirectionAndHonesty:
    def test_a_seam_whose_neighbour_runs_backwards_is_declared_and_welded(self):
        first = bilinear_patch((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0))
        second = bilinear_patch((1.0, 1.0, 0.0), (0.0, -1.0, 0.0), (1.0, 0.0, 0.0))
        seams = derive_seams((first, second))
        assert len(seams) == 1
        assert seams[0].reversed is True
        stitched = stitch(PatchNetwork((first, second), seams), nu=4, nv=4)
        assert stitched.seams[0].gap_max < 1e-15
        assert stitched.welded_nodes == 5        # nu + 1 nodes of the shared edge
        # Two squares sharing one edge are not a solid.
        assert meshlib.is_watertight(stitched.mesh) is False
        assert len(stitched.open_sides) == 6

    def test_the_same_pair_without_the_reversed_flag_is_reported_as_open_or_mismatched(self):
        first = bilinear_patch((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0))
        second = bilinear_patch((1.0, 1.0, 0.0), (0.0, -1.0, 0.0), (1.0, 0.0, 0.0))
        wrong = stitch(PatchNetwork((first, second), (Seam((0, "u1"), (1, "v0")),)), nu=4, nv=4)
        # The welding still happens (it is a declaration), so the honest signal is the
        # measured gap: declaring the direction wrongly takes the ends away.
        assert wrong.seams[0].reversed is False
        assert wrong.seams[0].gap_max == pytest.approx(1.0, abs=1e-12)
        assert wrong.mesh is not None

    def test_a_perturbed_patch_keeps_its_topology_and_loses_its_geometry(self):
        delta = 1e-3
        base = box_network()
        surfaces = list(base.surfaces)
        surfaces[5] = bilinear_patch(
            (BOX[0] + delta, 0.0, 0.0), (0.0, BOX[1], 0.0), (0.0, 0.0, BOX[2])
        )
        # (a) deriving the topology from the geometry now refuses the four sides of the
        # moved patch, because its corners are delta away from its neighbours'.
        derived = PatchNetwork(tuple(surfaces), derive_seams(tuple(surfaces)))
        assert len(derived.unclaimed_sides()) >= 4
        # (b) with the declaration kept from the unperturbed network, welding is an
        # index operation, so the shell stays watertight even though the two patches
        # disagree by delta.  Topology is not a proof of geometry; the seam report is.
        stitched = stitch(PatchNetwork(tuple(surfaces), base.seams), nu=4, nv=4)
        report = stitched.solid_report(reference_volume=BOX[0] * BOX[1] * BOX[2])
        assert report["watertight"] is True
        assert max(record.gap_max for record in stitched.seams) == pytest.approx(delta, rel=1e-9)
        assert abs(report["volume_relative_error"]) > 0.0

    def test_a_missing_patch_leaves_a_real_hole(self):
        network = box_network()
        open_box = PatchNetwork(network.surfaces[:5], tuple(
            seam for seam in network.seams if 5 not in (seam.first[0], seam.second[0])
        ))
        stitched = stitch(open_box, nu=4, nv=4)
        assert stitched.open_sides != ()
        assert meshlib.is_watertight(stitched.mesh) is False
        assert stitched.solid_report()["boundary_edges"] > 0

    def test_an_unclaimed_seam_that_coincides_bitwise_still_closes_by_position(self):
        # A second, undecorated welder exists: MeshBuilder keys vertices on the exact
        # coordinate tuple, so a seam that was left unclaimed but whose two sides
        # evaluate identically is still welded.  The declaration-level signal is
        # open_sides, and that is what a caller must read -- the mesh alone would
        # report a closed shell for an undeclared network.
        network = box_network()
        unclaimed = PatchNetwork(network.surfaces, network.seams[:-1])
        stitched = stitch(unclaimed, nu=4, nv=4)
        assert len(stitched.open_sides) == 2
        assert meshlib.is_watertight(stitched.mesh) is True

    def test_a_u_run_seam_against_a_v_run_seam_needs_matching_sample_counts(self):
        first = bilinear_patch((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0))
        second = bilinear_patch((1.0, 1.0, 0.0), (0.0, -1.0, 0.0), (1.0, 0.0, 0.0))
        network = PatchNetwork((first, second), (Seam((0, "u1"), (1, "v0")),))
        with pytest.raises(PatchNetworkError, match="nu == nv"):
            stitch(network, nu=4, nv=8)


# ---------------------------------------------------------------------------
# network validation and the constructors
# ---------------------------------------------------------------------------
class TestNetworkValidation:
    def test_a_side_claimed_twice_is_refused(self):
        surfaces = box_network().surfaces
        with pytest.raises(PatchNetworkError, match="claimed by seams"):
            PatchNetwork(surfaces, (Seam((0, "u0"), (1, "v0")), Seam((0, "u0"), (2, "v0"))))

    def test_an_out_of_range_patch_index_is_refused(self):
        surfaces = box_network().surfaces
        with pytest.raises(PatchNetworkError, match="the network has 6 patches"):
            PatchNetwork(surfaces, (Seam((0, "u0"), (9, "v0")),))

    def test_a_network_needs_a_patch(self):
        with pytest.raises(PatchNetworkError, match="at least one patch"):
            PatchNetwork(())

    def test_a_seam_needs_a_real_side_name(self):
        with pytest.raises(PatchNetworkError, match="side must be one of"):
            Seam((0, "w0"), (1, "v0"))

    def test_a_non_positive_sample_count_is_refused(self):
        with pytest.raises(PatchNetworkError, match="at least 1 x 1 samples"):
            stitch(box_network(), nu=0)


class TestConstructors:
    def test_a_bilinear_patch_is_a_parallelogram_with_the_declared_normal(self):
        patch = bilinear_patch((1.0, 2.0, 3.0), (2.0, 0.0, 0.0), (0.0, 3.0, 0.0))
        assert patch.evaluate(0.0, 0.0) == (1.0, 2.0, 3.0)
        assert patch.evaluate(1.0, 1.0) == (3.0, 5.0, 3.0)
        normal = patch.normal(0.5, 0.5)
        assert normal == pytest.approx((0.0, 0.0, 1.0), abs=1e-15)

    def test_a_revolved_patch_has_the_analytic_torus_patch_area(self):
        # Area of the patch phi in [0, pi/2], theta in [0, pi/2] of a torus:
        #   dA = (R + r cos(phi)) r dphi dtheta
        #   A  = dtheta * r * (R dphi + r (sin(phi1) - sin(phi0)))
        major, minor, sweep = 2.0, 0.5, 0.5 * math.pi
        patch = revolve_patch(meridian_arc(major, minor, 0.0, sweep), 0.0, sweep)
        exact = sweep * minor * (major * sweep + minor * math.sin(sweep))
        assert abs(patch.area() - exact) / exact < 1e-12
        # Not a sixteenth of the torus: the tube area is not uniform in phi.
        assert torus_area() / 16.0 < exact < 1.2 * torus_area() / 16.0

    def test_a_revolved_patch_rejects_a_sweep_beyond_a_half_turn(self):
        with pytest.raises(Exception):
            revolve_patch(meridian_arc(2.0, 0.5, 0.0, 0.5 * math.pi), 0.0, 1.5 * math.pi)


class TestReportShape:
    def test_the_seam_report_carries_a_row_per_seam_and_the_worst_numbers(self):
        stitched: Stitched = stitch(torus_network(), nu=4, nv=4)
        report = stitched.seam_report()
        assert report["patches"] == 16
        assert report["seam_count"] == 32
        assert len(report["seams"]) == 32
        assert report["gap_max_m"] == max(row["gap_max_m"] for row in report["seams"])
        assert report["angle_max_deg"] == max(row["angle_max_deg"] for row in report["seams"])
        assert report["welded_nodes"] == stitched.welded_nodes
        row = report["seams"][0]
        assert set(row) == {
            "label", "first", "second", "reversed", "samples",
            "gap_max_m", "gap_mean_m", "angle_max_deg", "angle_mean_deg", "inset_normals",
        }

    def test_the_solid_report_declines_to_invent_a_reference(self):
        report = stitch(box_network(), nu=4, nv=4).solid_report()
        assert "reference_volume_m3" not in report
        assert "volume_relative_error" not in report