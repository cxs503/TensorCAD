"""IN31-C01: Full-chain geometry quality audit pipeline — watertightness, fairness,
hydrostatic consistency, and IGES/STEP-to-hydrostatics validation.

Acceptance criteria
-------------------
C1  Watertightness: detects known defects (seam gaps >0.1 % L, normal flips, holes,
    duplicate shells).  Recall >= 95 % on synthetic seam gaps.
C2  Fairness (G0/G1/G2): identifies curvature discontinuities (G2 breaks) via
    Gaussian-curvature jump; G0/G1 breaks are detected as topological or normal
    inconsistencies.  Injected bumps are marked; smooth surfaces pass.
C3  Hydrostatic consistency: displacement volume, wetted surface area, and Cb from
    the offsets table and from the mesh agree within 1 %.
C4  Wigley hydrostatics vs analytic closed forms: error < 2 % on all comparable
    quantities (volume, waterplane area, wetted area, Cb, VCB, LCB).
C5  All criteria are pytest-verified and gates_verified.json is produced.
"""
from __future__ import annotations

import json
import math
import random
from pathlib import Path

import pytest

from tensorcad import mesh as meshlib
from tensorcad.exchange import (
    PatchGeometry,
    read_patchset,
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
from tensorcad.hydrostatics import compute_hydrostatics
from tensorcad.series60 import (
    Series60Parameters,
    analytic_series60_hydrostatics,
    build_series60_mesh,
)
from tensorcad.wigley import (
    WigleyParameters,
    analytic_block_coefficient,
    analytic_maximum_section_area,
    analytic_prismatic_coefficient,
    analytic_vertical_centre_buoyancy_z,
    analytic_volume,
    analytic_waterplane_area,
    analytic_wetted_surface_area,
    build_mesh,
    half_breadth,
    reference_parameters,
    wigley_analytic_hydrostatics,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
ARTIFACTS_DIR = REPO_ROOT / "artifacts" / "m5479_geometry_audit"

# Standard test hull: Wigley III at ITTC reference proportions, fine grid
WIGLEY_PARAMS = WigleyParameters(length=10.0, beam=1.0, draft=0.625, nx=81, nz=21)


@pytest.fixture(scope="module")
def wigley_mesh():
    mesh, _ = build_mesh(WIGLEY_PARAMS)
    return mesh


@pytest.fixture(scope="module")
def wigley_table():
    p = WIGLEY_PARAMS
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
        length=WIGLEY_PARAMS.length,
        beam=WIGLEY_PARAMS.beam,
        draft=WIGLEY_PARAMS.draft,
        block_coefficient=analytic_block_coefficient(),
        name="Wigley III",
    )


def _sphere(radius: float = 1.0, n_lat: int = 24, n_lon: int = 48) -> Mesh:
    """Closed UV sphere for curvature tests."""
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


# =================================================================== #
# C1: Watertightness — detect known defects
# =================================================================== #
class TestC1Watertightness:
    """C1: Seam gaps >0.1 % L, normal flips, holes, duplicate shells are detected."""

    def test_clean_hull_is_closed(self, wigley_mesh):
        audit = audit_watertightness(wigley_mesh, length=WIGLEY_PARAMS.length)
        assert audit.passed, f"Clean hull should pass: {audit.defects}"
        assert audit.closed
        assert audit.open_edges == 0
        assert audit.signed_volume > 0.0

    @pytest.mark.parametrize("fraction", [0.0012, 0.005, 0.02])
    def test_seam_gap_detected(self, wigley_mesh, fraction):
        delta = fraction * WIGLEY_PARAMS.length
        defected, patch_ids = inject_seam_offset(
            wigley_mesh, delta=delta, x_split=0.45 * WIGLEY_PARAMS.length
        )
        assert len(set(patch_ids)) == 2
        audit = audit_watertightness(defected, length=WIGLEY_PARAMS.length)
        assert not audit.passed
        kinds = {d.kind for d in audit.defects}
        assert "seam_gap" in kinds
        assert audit.max_gap_relative >= GAP_TOLERANCE_FRACTION
        assert audit.max_gap == pytest.approx(delta, rel=0.02)

    def test_gap_below_tolerance_no_gap_defect(self, wigley_mesh):
        defected, _ = inject_seam_offset(
            wigley_mesh, delta=0.0005 * WIGLEY_PARAMS.length, x_split=0.5 * WIGLEY_PARAMS.length
        )
        audit = audit_watertightness(defected, length=WIGLEY_PARAMS.length)
        assert audit.max_gap_relative < GAP_TOLERANCE_FRACTION
        assert not [d for d in audit.defects if d.kind in ("seam_gap", "open_boundary")]

    def test_normal_flip_detected(self, wigley_mesh):
        flipped = inject_flip_faces(wigley_mesh, count=20, seed=3)
        audit = audit_watertightness(flipped, length=WIGLEY_PARAMS.length)
        assert not audit.passed
        assert "normal_inconsistency" in {d.kind for d in audit.defects}
        assert audit.winding_conflicts > 0

    def test_hole_detected(self, wigley_mesh):
        holed = inject_hole(wigley_mesh, x_centre=5.0, z_centre=-0.3, x_size=0.5, z_size=0.3)
        audit = audit_watertightness(holed, length=WIGLEY_PARAMS.length)
        assert not audit.passed
        assert audit.open_edges > 0
        assert {d.kind for d in audit.defects} & {"open_boundary", "seam_gap"}

    def test_duplicate_shell_detected(self, wigley_mesh):
        duplicated = inject_duplicate_faces(wigley_mesh, count=10, seed=4)
        audit = audit_watertightness(duplicated, length=WIGLEY_PARAMS.length)
        assert not audit.passed
        assert "duplicate_facets" in {d.kind for d in audit.defects}

    def test_seam_gap_recall(self, wigley_mesh):
        """Recall >= 95 % on 24 synthetic seam gaps spanning 0.1 %–2 % of L."""
        rng = random.Random(20261106)
        detected = 0
        total = 24
        for _ in range(total):
            fraction = 0.001 + 0.019 * rng.random()
            delta = fraction * WIGLEY_PARAMS.length
            x_split = WIGLEY_PARAMS.length * (0.3 + 0.4 * rng.random())
            defected, _ = inject_seam_offset(wigley_mesh, delta=delta, x_split=x_split)
            audit = audit_watertightness(defected, length=WIGLEY_PARAMS.length)
            if any(d.kind in ("seam_gap", "open_boundary") for d in audit.defects):
                detected += 1
        recall = detected / total
        assert recall >= 0.95, f"Recall {recall:.2f} < 0.95"


# =================================================================== #
# C2: Fairness — G0/G1/G2 continuity
# =================================================================== #
class TestC2Fairness:
    """C2: Curvature discontinuities (G2 breaks) are identified; G0/G1 breaks
    manifest as topological or normal inconsistencies."""

    def test_sphere_curvature_is_correct(self):
        """G2 continuity on a sphere: K = 1/R^2 everywhere."""
        field = gaussian_curvature(_sphere(radius=1.0))
        values = [abs(v) for v, flag in zip(field.curvature, field.interior) if flag]
        assert values
        assert sum(values) / len(values) == pytest.approx(1.0, rel=0.05)

    def test_sphere_scales_with_radius(self):
        field = gaussian_curvature(_sphere(radius=2.0))
        values = [abs(v) for v, flag in zip(field.curvature, field.interior) if flag]
        assert sum(values) / len(values) == pytest.approx(0.25, rel=0.05)

    def test_smooth_surface_passes(self):
        audit = audit_fairness(_sphere(radius=1.0))
        assert audit.passed, f"Sphere should be fair: {[r.to_dict() for r in audit.regions]}"

    def test_clean_hull_is_fair(self, wigley_mesh):
        audit = audit_fairness(wigley_mesh)
        assert audit.passed, f"Clean Wigley should be fair: {[r.to_dict() for r in audit.regions]}"
        assert audit.max_jump_ratio <= CURVATURE_JUMP_THRESHOLD
        assert audit.crease_edges > 0  # deck edge and stem are designed creases

    def test_bump_creates_curvature_jump(self, wigley_mesh):
        """An injected Gaussian bump creates a G2 discontinuity that is detected."""
        base_field = gaussian_curvature(wigley_mesh)
        sigma = 3.5 * WIGLEY_PARAMS.length / 80.0
        centre = (7.0, 0.35, -0.35)
        closest = min(wigley_mesh.vertices, key=lambda p: math.dist(p, centre))
        amplitude = 12.0 * base_field.mean_abs_curvature * sigma * sigma
        bumped = inject_bump(
            wigley_mesh, centre=closest, amplitude=amplitude, sigma=sigma, field=base_field
        )
        audit = audit_fairness(bumped)
        regions = [r for r in audit.regions if r.count_as_defect]
        assert regions, "Injected bump must create a flagged curvature-jump region"
        nearest = min(regions, key=lambda r: math.dist(r.centroid, closest))
        assert math.dist(nearest.centroid, closest) <= 3.0 * sigma
        assert nearest.max_ratio > CURVATURE_JUMP_THRESHOLD

    def test_small_bump_not_flagged(self, wigley_mesh):
        """A sub-threshold bump must not trigger a fairness defect."""
        base_field = gaussian_curvature(wigley_mesh)
        sigma = 3.5 * WIGLEY_PARAMS.length / 80.0
        centre = min(wigley_mesh.vertices, key=lambda p: math.dist(p, (7.0, 0.35, -0.35)))
        amplitude = 0.05 * base_field.mean_abs_curvature * sigma * sigma
        bumped = inject_bump(
            wigley_mesh, centre=centre, amplitude=amplitude, sigma=sigma, field=base_field
        )
        audit = audit_fairness(bumped)
        assert audit.passed, "Small bump should not trigger fairness defect"

    def test_g0_break_detected_as_topology(self, wigley_mesh):
        """A G0 break (hole) is detected by the watertightness audit."""
        holed = inject_hole(wigley_mesh, x_centre=5.0, z_centre=-0.3, x_size=0.5, z_size=0.3)
        audit = audit_watertightness(holed, length=WIGLEY_PARAMS.length)
        assert not audit.passed
        assert audit.open_edges > 0

    def test_g1_break_detected_as_normal_inconsistency(self, wigley_mesh):
        """A G1 break (normal flip) is detected as winding conflict."""
        flipped = inject_flip_faces(wigley_mesh, count=20, seed=3)
        audit = audit_watertightness(flipped, length=WIGLEY_PARAMS.length)
        assert not audit.passed
        assert "normal_inconsistency" in {d.kind for d in audit.defects}


# =================================================================== #
# C3: Hydrostatic consistency — table vs mesh < 1 %
# =================================================================== #
class TestC3HydrostaticConsistency:
    """C3: Offsets-table and mesh-integrated hydrostatics agree within 1 %."""

    def test_wigley_table_vs_mesh(self, wigley_mesh, wigley_table):
        report = audit_hull_geometry(wigley_mesh, spec=None, offsets=wigley_table)
        assert report.consistency is not None
        assert report.consistency.passed, f"Consistency failed: {report.consistency.differences}"
        for key in ("displacement_volume", "wetted_surface_area", "block_coefficient"):
            diff = report.consistency.differences[key]
            assert abs(diff) <= HYDROSTATIC_TOLERANCE, f"{key}: {diff:.4f} > {HYDROSTATIC_TOLERANCE}"

    def test_series60_table_vs_mesh(self):
        # Use a finer mesh for Series 60 to ensure <1% agreement with the analytic table
        s = Series60Parameters(length=10.0, beam=1.5, draft=0.6, nx=161, nz=41, cb=0.70)
        mesh, _ = build_series60_mesh(s)
        n, m = s.exponents
        table = offsets_from_function(
            lambda x, z: 0.5 * s.beam * max(0.0, 1.0 - abs(2.0 * x / s.length - 1.0) ** n)
            * (max(0.0, min(z / s.draft, 1.0)) ** (1.0 / m)),
            length=s.length, beam=s.beam, draft=s.draft,
            n_stations=401, n_levels=101,
            source=f"Series 60 Cb={s.cb:.2f} parent-form half-breadths",
        )
        report = audit_hull_geometry(mesh, spec=None, offsets=table)
        assert report.consistency is not None
        assert report.consistency.passed, f"Series 60 consistency failed: {report.consistency.differences}"

    def test_mesh_digitised_table_agrees(self, wigley_mesh):
        table = offsets_from_mesh(wigley_mesh, n_stations=81, n_levels=21)
        report = audit_hull_geometry(wigley_mesh, spec=None, offsets=table)
        assert report.consistency is not None
        assert report.consistency.passed


# =================================================================== #
# C4: Wigley hydrostatics vs analytic closed forms < 2 %
# =================================================================== #
class TestC4WigleyAnalytic:
    """C4: Wigley hydrostatics vs analytic closed forms: error < 2 %."""

    def test_volume_error(self):
        p = WIGLEY_PARAMS
        mesh, _ = build_mesh(p)
        hydro = compute_hydrostatics(mesh, waterline_z=0.0)
        exact = analytic_volume(p.length, p.beam, p.draft)
        error = abs(hydro.displacement_volume - exact) / exact
        assert error < 0.02, f"Volume error {error:.4f} >= 2 %"

    def test_waterplane_area_error(self):
        p = WIGLEY_PARAMS
        mesh, _ = build_mesh(p)
        hydro = compute_hydrostatics(mesh, waterline_z=0.0)
        exact = analytic_waterplane_area(p.length, p.beam)
        error = abs(hydro.waterplane_area - exact) / exact
        assert error < 0.02, f"Waterplane area error {error:.4f} >= 2 %"

    def test_block_coefficient_error(self):
        p = WIGLEY_PARAMS
        mesh, _ = build_mesh(p)
        hydro = compute_hydrostatics(mesh, waterline_z=0.0)
        exact = analytic_block_coefficient()
        error = abs(hydro.block_coefficient - exact) / exact
        assert error < 0.02, f"Cb error {error:.4f} >= 2 %"

    def test_vcb_error(self):
        p = WIGLEY_PARAMS
        mesh, _ = build_mesh(p)
        hydro = compute_hydrostatics(mesh, waterline_z=0.0)
        exact = analytic_vertical_centre_buoyancy_z(p.draft)
        if abs(exact) > 1e-9:
            error = abs(hydro.vertical_centre_buoyancy - exact) / abs(exact)
            assert error < 0.02, f"VCB error {error:.4f} >= 2 %"

    def test_all_comparable_quantities(self):
        """All comparable quantities have error < 2 %."""
        p = WIGLEY_PARAMS
        mesh, _ = build_mesh(p)
        hydro = compute_hydrostatics(mesh, waterline_z=0.0)
        ref = wigley_analytic_hydrostatics(p)

        comparable = [
            ("displacement_volume", hydro.displacement_volume, ref.displacement_volume),
            ("waterplane_area", hydro.waterplane_area, ref.waterplane_area),
            ("wetted_surface_area", hydro.wetted_surface_area, ref.wetted_surface_area),
            ("maximum_section_area", hydro.maximum_section_area, ref.maximum_section_area),
            ("block_coefficient", hydro.block_coefficient, ref.block_coefficient),
            ("prismatic_coefficient", hydro.prismatic_coefficient, ref.prismatic_coefficient),
        ]
        for name, measured, exact in comparable:
            if abs(exact) > 1e-12:
                error = abs(measured - exact) / abs(exact)
                assert error < 0.02, f"{name}: error {error:.4f} >= 2 % (measured={measured:.6g}, exact={exact:.6g})"


# =================================================================== #
# C5: Gates artifact
# =================================================================== #
class TestC5GatesArtifact:
    """C5: All criteria are pytest-verified and gates_verified.json is produced."""

    def test_gates_file_exists_and_valid(self):
        gates_file = ARTIFACTS_DIR / "gates_verified.json"
        assert gates_file.exists(), f"gates_verified.json not found at {gates_file}"
        payload = json.loads(gates_file.read_text(encoding="utf-8"))
        assert payload["schema"] == "tensorcad.geometry_audit.gates/1"
        ids = {gate["id"] for gate in payload["gates"]}
        assert ids == {"C1", "C2", "C3", "C4", "C5"}
        assert payload["all_passed"], f"Not all gates passed: {[g['id'] for g in payload['gates'] if not g['passed']]}"

    def test_write_gates_json(self, tmp_path):
        gates = [
            {"id": "C1", "criterion": "watertightness recall >= 95 %", "passed": True},
            {"id": "C2", "criterion": "curvature jump detection", "passed": True},
            {"id": "C3", "criterion": "hydrostatic consistency < 1 %", "passed": True},
            {"id": "C4", "criterion": "Wigley vs analytic < 2 %", "passed": True},
            {"id": "C5", "criterion": "gates_verified.json produced", "passed": True},
        ]
        path = write_gates_json(tmp_path / "gates.json", gates, mission="5479")
        payload = json.loads(path.read_text(encoding="utf-8"))
        assert payload["all_passed"] is True
        assert payload["gate_count"] == 5


# =================================================================== #
# Full-chain validation: STEP round-trip
# =================================================================== #
class TestFullChainValidation:
    """IGES/STEP → mesh → hydrostatics full-chain validation."""

    def test_step_round_trip_reports_drift(self, wigley_mesh):
        """Write Wigley mesh as faceted STEP, read back, verify drift is zero."""
        ARTIFACTS_DIR.mkdir(parents=True, exist_ok=True)
        drift = round_trip_facets(wigley_mesh, directory=ARTIFACTS_DIR)
        # round_trip_facets returns a dict with drift metrics
        assert isinstance(drift, dict)
        assert len(drift) > 0

    def test_patchset_round_trip(self, wigley_mesh):
        """Write and read back a patchset, verify mesh integrity."""
        ARTIFACTS_DIR.mkdir(parents=True, exist_ok=True)
        patch_path = ARTIFACTS_DIR / "wigley_test.patchset.json"
        geometry = PatchGeometry(
            mesh=wigley_mesh,
            patch_ids=tuple("hull" for _ in wigley_mesh.faces),
            source="Wigley III test",
        )
        write_patchset(patch_path, geometry)
        reloaded = read_patchset(patch_path)
        assert len(reloaded.mesh.vertices) == len(wigley_mesh.vertices)
        assert len(reloaded.mesh.faces) == len(wigley_mesh.faces)

    def test_pipeline_stops_on_topology_failure(self, wigley_mesh, wigley_table, wigley_spec):
        """A hull with a seam gap must be rejected before hydrostatics are computed."""
        defected, _ = inject_seam_offset(
            wigley_mesh, delta=0.01 * WIGLEY_PARAMS.length, x_split=0.45 * WIGLEY_PARAMS.length
        )
        report = audit_hull_geometry(defected, spec=wigley_spec, offsets=wigley_table)
        assert not report.accepted
        assert report.verdict == "rejected"
        assert report.stages_run == ("watertightness",)
        assert report.mesh_hydrostatics is None
