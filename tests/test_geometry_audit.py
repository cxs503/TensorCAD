"""Tests for the hull-geometry audit pipeline (watertightness, fairness, hydrostatics).

The acceptance gates C1-C6 of mission 5283 live here; ``gates_verified.json`` is
produced by ``examples/geometry_audit/run_audit.py`` and this module re-derives the
numbers it asserts, so a stale artifact cannot outlive a code change.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from tensorcad import mesh as meshlib
from tensorcad.exchange import (
    ExchangeError,
    ExchangeUnavailable,
    PatchGeometry,
    load_exchange,
    read_step_facets,
    round_trip_facets,
    write_patchset,
    write_step_facets,
)
from tensorcad.geometry import Mesh, MeshBuilder
from tensorcad.geometry_audit import (
    CREASE_DEGREES,
    CURVATURE_JUMP_THRESHOLD,
    GAP_TOLERANCE_FRACTION,
    HYDROSTATIC_TOLERANCE,
    HullSpec,
    audit_declared_geometry,
    audit_fairness,
    audit_hull_geometry,
    audit_hydrostatic_consistency,
    audit_watertightness,
    gaussian_curvature,
    inject_bump,
    inject_cb_offset,
    inject_duplicate_faces,
    inject_flip_faces,
    inject_hole,
    inject_seam_offset,
    integrate_offsets,
    offsets_from_function,
    offsets_from_mesh,
    write_gates_json,
)
from tensorcad.series60 import (
    Series60Parameters,
    analytic_series60_hydrostatics,
    build_series60_mesh,
)
from tensorcad.wigley import (
    WigleyParameters,
    analytic_block_coefficient,
    analytic_volume,
    analytic_wetted_surface_area,
    build_mesh,
    half_breadth,
    reference_parameters,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
GATES_FILE = REPO_ROOT / "artifacts" / "m5283_geometry_audit" / "gates_verified.json"

WIGLEY = WigleyParameters(length=10.0, beam=1.0, draft=0.625, nx=81, nz=21)


@pytest.fixture(scope="module")
def wigley():
    mesh, _ = build_mesh(WIGLEY)
    return mesh


@pytest.fixture(scope="module")
def wigley_table():
    p = reference_parameters(length=10.0, beam=1.0, draft=0.625, nx=81, nz=21)
    return offsets_from_function(
        lambda x, z: half_breadth(p, x, p.draft - z),
        length=p.length,
        beam=p.beam,
        draft=p.draft,
        source="Wigley III analytic half-breadths",
    )


@pytest.fixture(scope="module")
def wigley_spec():
    return HullSpec(
        length=10.0, beam=1.0, draft=0.625,
        block_coefficient=analytic_block_coefficient(), name="Wigley III",
    )


def sphere(radius: float = 1.0, n_lat: int = 24, n_lon: int = 48) -> Mesh:
    """Closed UV sphere: a surface whose Gaussian curvature is exactly 1/R^2."""
    builder = MeshBuilder()
    north = builder.add_vertex((0.0, 0.0, radius))
    rings: list[list[int]] = []
    for i in range(1, n_lat):
        theta = math.pi * i / n_lat
        ring = []
        for j in range(n_lon):
            phi = 2.0 * math.pi * j / n_lon
            ring.append(
                builder.add_vertex(
                    (
                        radius * math.sin(theta) * math.cos(phi),
                        radius * math.sin(theta) * math.sin(phi),
                        radius * math.cos(theta),
                    )
                )
            )
        rings.append(ring)
    south = builder.add_vertex((0.0, 0.0, -radius))
    for j in range(n_lon):
        builder.add_triangle(north, rings[0][j], rings[0][(j + 1) % n_lon])
    for i in range(len(rings) - 1):
        for j in range(n_lon):
            k = (j + 1) % n_lon
            builder.add_quad(rings[i][j], rings[i + 1][j], rings[i + 1][k], rings[i][k])
    for j in range(n_lon):
        builder.add_triangle(south, rings[-1][(j + 1) % n_lon], rings[-1][j])
    return builder.build()


# --------------------------------------------------------------------------- #
# C1: watertightness
# --------------------------------------------------------------------------- #
class TestWatertightness:
    def test_clean_hull_is_closed_and_accepted(self, wigley):
        audit = audit_watertightness(wigley, length=WIGLEY.length)
        assert audit.passed, audit.defects
        assert audit.closed
        assert audit.open_edges == 0
        assert audit.signed_volume > 0.0
        assert audit.face_components == 1

    @pytest.mark.parametrize("fraction", [0.0012, 0.005, 0.02])
    def test_seam_gap_detected_with_correct_size(self, wigley, fraction):
        # 0.0012 L rather than 0.001 L: the gate is a strict inequality, so at exactly
        # the tolerance the model is still accepted (documented in the module).
        delta = fraction * WIGLEY.length
        defected, patch_ids = inject_seam_offset(wigley, delta=delta, x_split=0.45 * WIGLEY.length)
        assert len(set(patch_ids)) == 2
        audit = audit_watertightness(defected, length=WIGLEY.length)
        assert not audit.passed
        kinds = {defect.kind for defect in audit.defects}
        assert "seam_gap" in kinds
        assert audit.max_gap_relative >= GAP_TOLERANCE_FRACTION
        # The longitudinal offset is an exact ground truth: the measured gap is the
        # distance the surface has to move, to within the discretisation of the ring.
        assert audit.max_gap == pytest.approx(delta, rel=0.02)

    def test_gap_below_tolerance_raises_no_gap_defect(self, wigley):
        defected, _ = inject_seam_offset(wigley, delta=0.0005 * WIGLEY.length,
                                         x_split=0.5 * WIGLEY.length)
        audit = audit_watertightness(defected, length=WIGLEY.length)
        assert audit.max_gap_relative < GAP_TOLERANCE_FRACTION
        assert not [d for d in audit.defects if d.kind in ("seam_gap", "open_boundary")]
        # Disconnected shells are still a defect: the two pieces are separate solids.
        assert "component_count" in {defect.kind for defect in audit.defects}

    def test_hole_detected(self, wigley):
        holed = inject_hole(wigley, x_centre=5.0, z_centre=-0.3, x_size=0.5, z_size=0.3)
        audit = audit_watertightness(holed, length=WIGLEY.length)
        assert not audit.passed
        assert audit.open_edges > 0
        # The hole is punched through both sides of a symmetric hull, so there are
        # two free boundaries and the audit classifies them as a seam; either
        # classification is a hole and the size is what the gate acts on.
        assert {defect.kind for defect in audit.defects} & {"open_boundary", "seam_gap"}
        assert audit.max_gap_relative > GAP_TOLERANCE_FRACTION

    def test_normal_flip_detected(self, wigley):
        flipped = inject_flip_faces(wigley, count=20, seed=3)
        audit = audit_watertightness(flipped, length=WIGLEY.length)
        assert not audit.passed
        assert "normal_inconsistency" in {defect.kind for defect in audit.defects}
        assert audit.winding_conflicts > 0
        assert audit.conflicting_area > 0.0

    def test_duplicate_shell_detected(self, wigley):
        duplicated = inject_duplicate_faces(wigley, count=10, seed=4)
        audit = audit_watertightness(duplicated, length=WIGLEY.length)
        assert not audit.passed
        kinds = {defect.kind for defect in audit.defects}
        assert "duplicate_facets" in kinds
        assert "nonmanifold_edge" in kinds

    def test_pipeline_stops_after_topology_failure(self, wigley, wigley_table, wigley_spec):
        defected, _ = inject_seam_offset(wigley, delta=0.01 * WIGLEY.length,
                                         x_split=0.45 * WIGLEY.length)
        report = audit_hull_geometry(defected, spec=wigley_spec, offsets=wigley_table)
        assert not report.accepted
        assert report.verdict == "rejected"
        assert report.stages_run == ("watertightness",)
        assert report.mesh_hydrostatics is None


# --------------------------------------------------------------------------- #
# curvature kernel and fairness
# --------------------------------------------------------------------------- #
class TestCurvatureKernel:
    def test_sphere_curvature_is_one_over_radius_squared(self):
        field = gaussian_curvature(sphere(radius=1.0))
        values = [abs(value) for value, flag in zip(field.curvature, field.interior) if flag]
        assert values
        assert sum(values) / len(values) == pytest.approx(1.0, rel=0.05)

    def test_sphere_scales_with_radius(self):
        field = gaussian_curvature(sphere(radius=2.0))
        values = [abs(value) for value, flag in zip(field.curvature, field.interior) if flag]
        assert sum(values) / len(values) == pytest.approx(0.25, rel=0.05)

    def test_smooth_surface_has_no_flagged_region(self):
        audit = audit_fairness(sphere(radius=1.0))
        assert audit.passed, [region.to_dict() for region in audit.regions]


class TestFairness:
    def test_clean_hull_is_fair_and_creases_are_skipped(self, wigley):
        audit = audit_fairness(wigley)
        assert audit.passed, [region.to_dict() for region in audit.regions]
        assert audit.max_jump_ratio <= CURVATURE_JUMP_THRESHOLD
        # The deck edge and the zero-breadth stem are designed creases: they are
        # excluded by dihedral, not silently smoothed.
        assert audit.crease_edges > 0
        assert audit.edges_considered > 0

    def test_bump_is_marked(self, wigley):
        field = gaussian_curvature(wigley)
        sigma = 3.5 * WIGLEY.length / 80.0
        centre = (7.0, 0.35, -0.35)
        closest = min(wigley.vertices, key=lambda point: math.dist(point, centre))
        amplitude = 12.0 * field.mean_abs_curvature * sigma * sigma
        bumped = inject_bump(wigley, centre=closest, amplitude=amplitude, sigma=sigma, field=field)
        audit = audit_fairness(bumped)
        regions = [region for region in audit.regions if region.count_as_defect]
        assert regions, "an injected curvature break must be marked"
        nearest = min(regions, key=lambda region: math.dist(region.centroid, closest))
        assert math.dist(nearest.centroid, closest) <= 3.0 * sigma
        assert nearest.max_ratio > CURVATURE_JUMP_THRESHOLD

    def test_small_jump_is_not_marked(self, wigley):
        field = gaussian_curvature(wigley)
        sigma = 3.5 * WIGLEY.length / 80.0
        centre = min(wigley.vertices, key=lambda point: math.dist(point, (7.0, 0.35, -0.35)))
        # A twentieth of the amplitude used above: a 3 mm swell over a 0.4 m patch,
        # which changes K by less than the 50 % criterion.
        amplitude = 0.05 * field.mean_abs_curvature * sigma * sigma
        bumped = inject_bump(wigley, centre=centre, amplitude=amplitude, sigma=sigma, field=field)
        audit = audit_fairness(bumped)
        assert audit.passed

    def test_crease_threshold_is_used(self, wigley):
        audit = audit_fairness(wigley, crease_degrees=CREASE_DEGREES)
        strict = audit_fairness(wigley, crease_degrees=0.0)
        # crease_degrees=0 turns every non-flat edge into a "designed crease", so
        # almost nothing is compared: the exclusion is what keeps the deck edge out
        # of the flag list, and switching it off must shrink the compared set.
        assert strict.edges_considered < audit.edges_considered
        assert strict.crease_edges >= audit.crease_edges


# --------------------------------------------------------------------------- #
# C3 / C5: hydrostatics
# --------------------------------------------------------------------------- #
class TestHydrostatics:
    def test_offsets_table_matches_analytic_volume(self, wigley_table):
        table = integrate_offsets(wigley_table)
        reference = analytic_volume(WIGLEY.length, WIGLEY.beam, WIGLEY.draft)
        assert table.displacement_volume == pytest.approx(reference, rel=1e-3)

    def test_table_and_mesh_agree_within_one_percent(self, wigley, wigley_table):
        report = audit_hull_geometry(wigley, spec=None, offsets=wigley_table)
        assert report.consistency is not None
        assert report.consistency.passed, report.consistency.differences
        for key in ("displacement_volume", "wetted_surface_area", "block_coefficient"):
            assert abs(report.consistency.differences[key]) <= HYDROSTATIC_TOLERANCE

    def test_mesh_digitised_table_agrees(self, wigley):
        table = offsets_from_mesh(wigley, n_stations=81, n_levels=41)
        report = audit_hull_geometry(wigley, spec=None, offsets=table)
        assert report.consistency is not None
        assert report.consistency.passed, report.consistency.differences

    def test_inconsistent_table_is_caught(self, wigley):
        table = offsets_from_mesh(wigley, n_stations=41, n_levels=21)
        # Inflate the table by 3 % in the vertical: the two routes must disagree.
        inflated = type(table)(
            stations=table.stations,
            levels=table.levels,
            half_breadths=tuple(tuple(value * 1.03 for value in row) for row in table.half_breadths),
            length=table.length,
            beam=table.beam,
            draft=table.draft,
            source="inflated control",
        )
        from tensorcad.hydrostatics import compute_hydrostatics

        mesh_hydro = compute_hydrostatics(wigley, 0.0, water_density=1025.0, closed=True)
        audit = audit_hydrostatic_consistency(integrate_offsets(inflated), mesh_hydro)
        assert not audit.passed
        assert audit.defects

    def test_wigley_against_closed_forms(self, wigley, wigley_table, wigley_spec):
        report = audit_hull_geometry(wigley, spec=wigley_spec, offsets=wigley_table)
        assert report.accepted, [defect.to_dict() for defect in report.defects]
        hydro = report.mesh_hydrostatics
        assert hydro.displacement_volume == pytest.approx(
            analytic_volume(WIGLEY.length, WIGLEY.beam, WIGLEY.draft), rel=0.02
        )
        assert hydro.wetted_surface_area == pytest.approx(
            analytic_wetted_surface_area(WIGLEY.length, WIGLEY.beam, WIGLEY.draft),
            rel=0.02,
        )
        assert hydro.block_coefficient == pytest.approx(analytic_block_coefficient(), rel=0.02)

    def test_series60_against_parent_form(self):
        parameters = Series60Parameters(length=10.0, beam=1.5, draft=0.6, nx=161, nz=41, cb=0.70)
        mesh, _ = build_series60_mesh(parameters)
        spec = HullSpec(length=parameters.length, beam=parameters.beam,
                        draft=parameters.draft, block_coefficient=0.70, name="Series 60 Cb=0.70")
        report = audit_hull_geometry(mesh, spec=spec, check_fairness=False)
        reference = analytic_series60_hydrostatics(parameters)
        hydro = report.mesh_hydrostatics
        assert hydro.displacement_volume == pytest.approx(reference["displacement_volume"], rel=0.02)
        assert hydro.block_coefficient == pytest.approx(reference["block_coefficient"], rel=0.02)
        assert report.declared.passed, report.declared.defects
        assert report.consistency.passed, report.consistency.differences


# --------------------------------------------------------------------------- #
# C4: declared geometry and negative cases
# --------------------------------------------------------------------------- #
class TestNegativeCases:
    def test_cb_offset_of_five_percent_is_refused(self, wigley, wigley_table, wigley_spec):
        shifted = inject_cb_offset(wigley, factor=1.05)
        stable = audit_declared_geometry(shifted, wigley_spec, volume=None)
        # L, B and T are untouched by the injector: only Cb moves.
        assert stable.relative("length") == pytest.approx(0.0, abs=1e-9)
        assert stable.relative("beam") == pytest.approx(0.0, abs=1e-9)
        assert stable.relative("draft") == pytest.approx(0.0, abs=1e-9)
        assert stable.relative("block_coefficient") > 0.03

        report = audit_hull_geometry(shifted, spec=wigley_spec, offsets=wigley_table)
        assert not report.accepted
        kinds = {defect.kind for defect in report.defects}
        assert "block_coefficient_mismatch" in kinds
        reported = next(d for d in report.defects if d.kind == "block_coefficient_mismatch")
        assert reported.relative > 0.03
        assert "realised Cb" in reported.detail
        assert report.stages_run[-1] == "declared_geometry"

    def test_dimension_mismatch_is_refused(self, wigley, wigley_table):
        wrong = HullSpec(length=11.0, beam=1.0, draft=0.625, block_coefficient=None)
        report = audit_hull_geometry(wigley, spec=wrong, offsets=wigley_table)
        assert not report.accepted
        assert "length_mismatch" in {defect.kind for defect in report.defects}

    def test_bad_declaration_values_rejected(self):
        with pytest.raises(ValueError):
            HullSpec(length=0.0, beam=1.0, draft=1.0)
        with pytest.raises(ValueError):
            HullSpec(length=1.0, beam=1.0, draft=1.0, block_coefficient=1.5)


# --------------------------------------------------------------------------- #
# exchange seam
# --------------------------------------------------------------------------- #
class TestExchange:
    def test_faceted_step_round_trip_is_exact(self, wigley):
        trip = round_trip_facets(wigley)
        # The file is written in millimetres, so the round trip is a *1000 then /1000:
        # the residual is one ulp of the division, not a modelling error.
        assert trip["vertex_deviation_m"] < 1e-12
        assert int(trip["triangle_count"]) == len(wigley.faces)

    def test_defect_survives_the_step_round_trip(self, wigley, tmp_path):
        defected, _ = inject_seam_offset(wigley, delta=0.01 * WIGLEY.length,
                                         x_split=0.45 * WIGLEY.length)
        path = write_step_facets(tmp_path / "defected.step", defected)
        reimported = read_step_facets(path)
        audit = audit_watertightness(reimported.mesh, length=WIGLEY.length)
        assert not audit.passed
        assert "seam_gap" in {defect.kind for defect in audit.defects}
        assert audit.max_gap == pytest.approx(0.01 * WIGLEY.length, rel=0.02)

    def test_patchset_round_trip_keeps_patch_ids(self, wigley, tmp_path):
        geometry = PatchGeometry(
            mesh=wigley, patch_ids=tuple("hull" if i % 2 else "deck" for i in range(len(wigley.faces))),
            source="test",
        )
        path = write_patchset(tmp_path / "hull.patchset.json", geometry)
        loaded = load_exchange(path)
        assert loaded.patch_ids == geometry.patch_ids
        assert loaded.mesh.faces == wigley.faces

    def test_brep_step_is_refused_with_named_entities(self, tmp_path):
        path = tmp_path / "brep.step"
        path.write_text(
            "ISO-10303-21;\nHEADER;\nENDSEC;\nDATA;\n"
            "#1 = ADVANCED_FACE('',(#2),#3,.T.);\n"
            "#3 = CYLINDRICAL_SURFACE('',#4,0.5);\n"
            "ENDSEC;\nEND-ISO-10303-21;\n",
            encoding="utf-8",
        )
        with pytest.raises(ExchangeError) as error:
            read_step_facets(path)
        assert "ADVANCED_FACE" in str(error.value) or "CYLINDRICAL_SURFACE" in str(error.value)

    def test_iges_route_is_honest_about_its_dependency(self, tmp_path):
        path = tmp_path / "hull.iges"
        path.write_text("probe\n", encoding="utf-8")
        with pytest.raises(ExchangeUnavailable):
            load_exchange(path)


# --------------------------------------------------------------------------- #
# gates artifact
# --------------------------------------------------------------------------- #
class TestGatesArtifact:
    def test_gates_file_is_complete(self):
        if not GATES_FILE.exists():
            pytest.skip("run examples/geometry_audit/run_audit.py to produce gates_verified.json")
        payload = json.loads(GATES_FILE.read_text(encoding="utf-8"))
        assert payload["schema"] == "tensorcad.geometry_audit.gates/1"
        ids = {gate["id"] for gate in payload["gates"]}
        assert ids == {"C1", "C2", "C3", "C4", "C5", "C6"}
        assert payload["all_passed"], [g["id"] for g in payload["gates"] if not g["passed"]]
        for gate in payload["gates"]:
            assert gate["criterion"]

    def test_write_gates_json_requires_every_pass(self, tmp_path):
        path = write_gates_json(
            tmp_path / "gates.json",
            [{"id": "C1", "criterion": "x", "passed": True},
             {"id": "C2", "criterion": "y", "passed": False}],
            mission="test",
        )
        payload = json.loads(path.read_text(encoding="utf-8"))
        assert payload["all_passed"] is False
        assert payload["gate_count"] == 2


# --------------------------------------------------------------------------- #
# mensuration sanity of the module's own helpers
# --------------------------------------------------------------------------- #
def test_simpson_matches_analytic_for_polynomial():
    from tensorcad.geometry_audit import _simpson_or_trapezoid

    values = [x**2 for x in (0.0, 0.25, 0.5, 0.75, 1.0)]
    assert _simpson_or_trapezoid(values, 0.25) == pytest.approx(1.0 / 3.0, rel=1e-12)


def test_mesh_volume_of_a_unit_cube():
    builder = MeshBuilder()
    corners = [
        builder.add_vertex(point)
        for point in (
            (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
            (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0), (0.0, 1.0, 1.0),
        )
    ]
    quads = (
        (0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4),
        (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7),
    )
    for a, b, c, d in quads:
        builder.add_quad(corners[a], corners[b], corners[c], corners[d])
    cube = meshlib.orient_outward(builder.build()).mesh
    assert meshlib.total_signed_volume(cube) == pytest.approx(1.0, rel=1e-12)
    audit = audit_watertightness(cube, length=1.0)
    assert audit.passed, audit.defects
    # A cube is flat everywhere: no curvature jump, and no false positive from the
    # twelve 90-degree edges, which the crease filter removes.
    fairness = audit_fairness(cube)
    assert fairness.passed
    assert fairness.crease_edges == 12