"""High-order geometry gates (J1-J7) for ``tensorcad.highorder``.

The gates are stated as claims that can fail, and each one is paired with the evidence
that its claim is not vacuous:

======  ==========================================================================
J1      Gauss rule exactness order is ``2p - 1`` and curved-patch quadrature decays at
        least that fast (p = 1..4). Paired with J6, which shows a faceted patch does
        not.
J2      Normals and principal curvatures agree with the closed-form ellipsoid/sphere
        truth to ``1e-3`` of ``kappa_max``, and the exact Bezier-extraction route
        agrees with a finite-difference route to the same tolerance.
J3      Volume of the isoparametric element set (order 4, 2x2 subdivisions of every
        knot span) matches the analytic sphere/ellipsoid volume to ``1e-6``.
J4      Ellipsoid hull ``Cb`` and wetted surface ``S`` from three routes -- NURBS
        quadrature, STL tessellation, voxel grid -- agree with the AGM/Carlson
        elliptic-integral value to ``0.5%``.
J5      The ``(x, w, n, |J|)`` contract round-trips bit-identically (NPZ/HDF5/JSON),
        its content hash is stable, and the integrals recomputed from the arrays equal
        the direct ones to ``1e-14``.
J6      Negative control: the linearized patch drops to second-order convergence
        (<= 2) and its volume error exceeds ``1e-3`` -- so the J1/J3 numbers above are
        not artefacts of the harness.
J7      This module runs inside the 120 s CI budget.
======  ==========================================================================

Everything is deterministic: no random sampling, no wall-clock-dependent choices. The
parameter grids are fixed constants in this file, so a failure reproduces.
"""

from __future__ import annotations

import json
import math
import time
from pathlib import Path

import pytest

from tensorcad.highorder import contract as contract
from tensorcad.highorder import elliptic, quadrature
from tensorcad.highorder import nurbs_geom as geom
from tensorcad.nurbs import NURBSError, NURBSSurface, sphere_surface

#: Committed benchmark directory guarded by :class:`TestCommittedEvidence`.
CASE = "cad_highorder_geometry"

#: J7 needs this module's own runtime, not "time since collection". Pytest imports every conftest
#: and test module *before* running anything, so a timer started at import would charge my module
#: for tests that ran before it. The first test of the module stamps the clock instead, and the
#: last test compares against that stamp.
_TIMING: dict[str, float] = {}


@pytest.fixture(autouse=True)
def _module_clock() -> None:
    _TIMING.setdefault("first_test", time.perf_counter())

#: Ellipsoid used throughout: distinct semi-axes, so no accidental symmetry hides a bug.
ELLIPSOID = (1.0, 0.6, 0.42)

#: Tolerance for the curvature gates (fraction of ``kappa_max``).
CURVATURE_TOLERANCE = 1e-3


def bezier_patch(degree: int) -> NURBSSurface:
    """A curved non-rational Bezier patch of the requested degree (deterministic)."""
    p = degree
    net: list[tuple[tuple[float, float, float], ...]] = []
    for i in range(p + 1):
        s = i / p if p else 0.0
        row: list[tuple[float, float, float]] = []
        for j in range(p + 1):
            t = j / p if p else 0.0
            z = 0.35 * math.sin(math.pi * s) * math.sin(math.pi * t) + 0.15 * s * t
            row.append((s, 0.7 * t, z))
        net.append(tuple(row))
    knots = (0.0,) * (p + 1) + (1.0,) * (p + 1)
    return NURBSSurface(degree=(p, p), knots=(knots, knots), control_net=tuple(net))


def patch_bounds(surface: NURBSSurface, samples: int = 33) -> dict[str, float]:
    """Axis-aligned bounds of the patch, sampled on a parameter grid."""
    (u_low, u_high), (v_low, v_high) = surface.domain
    xs: list[float] = []
    ys: list[float] = []
    zs: list[float] = []
    for i in range(samples):
        u = u_low + (u_high - u_low) * i / (samples - 1)
        for j in range(samples):
            v = v_low + (v_high - v_low) * j / (samples - 1)
            point = surface.evaluate(u, v)
            xs.append(point[0])
            ys.append(point[1])
            zs.append(point[2])
    return {
        "L": max(xs) - min(xs),
        "B": max(ys) - min(ys),
        "T": max(zs) - min(zs),
    }


# ----------------------------------------------------------------------------------
# J1 -- quadrature order
# ----------------------------------------------------------------------------------
class TestJ1QuadratureOrder:
    """The rule's exactness order and the observed convergence order of the integrals."""

    @pytest.mark.parametrize("order", [1, 2, 3, 4])
    def test_exactness_order_is_two_p_minus_one(self, order: int) -> None:
        """``p`` points integrate degree ``2p - 1`` exactly and degree ``2p`` not at all."""
        probe = quadrature.gauss_exactness_probe(order)
        assert probe["degree_exact"] == 2 * order - 1
        assert probe["exact_at_limit"], probe
        assert probe["fails_beyond_limit"], probe

    @pytest.mark.parametrize("degree", [1, 2, 3, 4])
    def test_curved_patch_area_converges_at_least_at_2p_minus_1(self, degree: int) -> None:
        """Observed order of the area integral on a degree-``p`` patch is >= 2p - 1."""
        surface = bezier_patch(degree)
        field = lambda u, v: quadrature.area_element(surface, u, v)  # noqa: E731
        study = quadrature.convergence_study(
            surface, field, orders=(2, 3, 4, 6, 8), reference_order=32
        )
        assert study["fitted_order"] is not None
        assert study["fitted_order"] >= 2 * degree - 1, study
        # Above the round-off floor the error must fall monotonically: a rule that is
        # refined but not converging is a bug the fitted slope alone could hide.
        measurable = [
            row["relative_error"]
            for row in study["rows"]
            if row["relative_error"] > 10.0 * quadrature.NOISE_FLOOR
        ]
        assert len(measurable) >= 2, study["rows"]
        assert measurable == sorted(measurable, reverse=True), measurable

    def test_observations_refuse_to_report_order_at_the_noise_floor(self) -> None:
        """A pair whose finer error is round-off yields ``None``, not a fabricated order."""
        orders = quadrature.observed_orders([2, 4, 8], [1e-3, 1e-9, 0.0])
        assert orders[0] is not None and orders[0] > 1.0
        assert orders[1] is None

    def test_knot_spans_are_the_distinct_positive_intervals(self) -> None:
        knots = (0.0, 0.0, 0.0, 0.25, 0.5, 0.5, 1.0, 1.0, 1.0)
        spans = quadrature.knot_spans(knots, 2)
        assert spans == ((0.0, 0.25), (0.25, 0.5), (0.5, 1.0))


# ----------------------------------------------------------------------------------
# J2 -- normals and curvature
# ----------------------------------------------------------------------------------
class TestJ2NormalsAndCurvature:
    """Analytic differential geometry against closed-form truth and against differences."""

    def test_knot_insertion_preserves_the_surface(self) -> None:
        """The enabling operation for exact second derivatives must not move the patch."""
        surface = sphere_surface(1.0)
        inserted = geom.insert_knot(surface, 0, 0.1, 2)
        inserted = geom.insert_knot(inserted, 1, 1.4, 1)
        for u, v in ((0.02, 0.3), (0.1, 1.4), (0.55, 1.9)):
            assert inserted.evaluate(u, v) == pytest.approx(surface.evaluate(u, v), abs=1e-14)

    def test_restricted_patch_reproduces_the_parent_parameters(self) -> None:
        surface = geom.ellipsoid_surface(*ELLIPSOID)
        half = geom.half_ellipsoid_patch(*ELLIPSOID)
        for u, v in ((0.0, 1.0), (0.25, 1.3), (1.0, 2.0), (0.7, 1.75)):
            assert half.evaluate(u, v) == pytest.approx(surface.evaluate(u, v), abs=1e-14)

    def test_bezier_extraction_unlocks_second_derivatives(self) -> None:
        """The parent refuses ``order=2`` in a C0 direction; the extracted span does not."""
        surface = sphere_surface(1.0)
        with pytest.raises(NURBSError):
            surface.derivatives(0.5, 1.4, 2)
        span = geom.span_patch(surface, 0.5, 1.4)
        assert span.derivatives(0.5, 1.4, 2)[(2, 0)] != (0.0, 0.0, 0.0)

    def test_sphere_curvature_is_the_inverse_radius(self) -> None:
        radius = 1.5
        surface = sphere_surface(radius)
        kappa = 1.0 / radius
        points = geom.curvature_scan(surface, grid=5)
        worst = max(
            abs(abs(value) - kappa) / kappa
            for point in points
            for value in point.principal_curvatures
        )
        assert worst <= CURVATURE_TOLERANCE, worst
        # Sign convention: with the outward normal S_u x S_v, a sphere curves away from it,
        # so the patch reports -1/R while the geometric curvature magnitude is +1/R.
        assert all(value < 0.0 for point in points for value in point.principal_curvatures)

    def test_ellipsoid_curvature_matches_the_closed_form(self) -> None:
        a, b, c = ELLIPSOID
        surface = geom.ellipsoid_surface(a, b, c)
        points = geom.curvature_scan(surface, grid=5)
        kappa_max = max(
            max(geom.ellipsoid_principal_curvatures(a, b, c, point.point)) for point in points
        )
        worst_curvature = 0.0
        worst_normal = 0.0
        for point in points:
            closed = sorted(geom.ellipsoid_principal_curvatures(a, b, c, point.point), reverse=True)
            measured = sorted((abs(value) for value in point.principal_curvatures), reverse=True)
            for left, right in zip(measured, closed, strict=True):
                worst_curvature = max(worst_curvature, abs(left - right) / kappa_max)
            exact_normal = geom.ellipsoid_normal(a, b, c, point.point)
            deviation = math.sqrt(sum((x - y) ** 2 for x, y in zip(point.normal, exact_normal, strict=True)))
            worst_normal = max(worst_normal, deviation)
        assert worst_curvature <= CURVATURE_TOLERANCE, worst_curvature
        assert worst_normal <= CURVATURE_TOLERANCE, worst_normal

    def test_finite_difference_route_confirms_the_analytic_one(self) -> None:
        """The tangent field differentiated numerically must land on the same curvature."""
        a, b, c = ELLIPSOID
        surface = geom.ellipsoid_surface(a, b, c)
        points = geom.curvature_scan(surface, grid=4)
        kappa_max = max(
            max(geom.ellipsoid_principal_curvatures(a, b, c, point.point)) for point in points
        )
        worst = 0.0
        for point in points:
            difference = geom.finite_difference_geometry(surface, point.u, point.v)
            closed = sorted(geom.ellipsoid_principal_curvatures(a, b, c, point.point), reverse=True)
            measured = sorted((abs(value) for value in difference.principal_curvatures), reverse=True)
            for left, right in zip(measured, closed, strict=True):
                worst = max(worst, abs(left - right) / kappa_max)
        assert worst <= CURVATURE_TOLERANCE, worst
        assert worst > 0.0  # a difference scheme that returns the analytic answer is suspect

    def test_linear_direction_is_refused_rather_than_zeroed(self) -> None:
        surface = bezier_patch(1)
        with pytest.raises(NURBSError):
            geom.differential_geometry(surface, 0.5, 0.5)
        assert geom.finite_difference_geometry(surface, 0.5, 0.5).second_derivatives == (
            "finite_difference"
        )


# ----------------------------------------------------------------------------------
# J3 -- curved element volume
# ----------------------------------------------------------------------------------
class TestJ3CurvedElementVolume:
    """The exported element set must integrate the volume to 1e-6 relative."""

    def test_sphere_volume_from_order4_elements_on_subdivided_spans(self) -> None:
        radius = 1.0
        surface = sphere_surface(radius)
        elements = quadrature.curved_elements(surface, order=4, subdivisions=(2, 2))
        assert len(elements) == 24  # 6 knot spans x 2 x 2 subdivisions
        assert {element.degree for element in elements} == {(2, 2)}
        volume = math.fsum(element.volume_contribution for element in elements)
        exact = 4.0 / 3.0 * math.pi * radius**3
        assert abs(volume - exact) / exact <= 1e-6

    def test_ellipsoid_volume_from_order4_elements_on_subdivided_spans(self) -> None:
        a, b, c = ELLIPSOID
        surface = geom.ellipsoid_surface(a, b, c)
        elements = quadrature.curved_elements(surface, order=4, subdivisions=(2, 2))
        volume = math.fsum(element.volume_contribution for element in elements)
        exact = 4.0 / 3.0 * math.pi * a * b * c
        assert abs(volume - exact) / exact <= 1e-6

    def test_element_weights_sum_to_the_patch_area(self) -> None:
        surface = sphere_surface(1.0)
        elements = quadrature.curved_elements(surface, order=4, subdivisions=(2, 2))
        area = math.fsum(element.area for element in elements)
        direct = quadrature.surface_area(surface, order=12)
        assert abs(area - direct) / direct <= 1e-6

    def test_a_single_element_per_span_is_not_enough_and_says_so(self) -> None:
        """The 1e-6 gate needs a subdivided span; one order-4 element misses it."""
        surface = sphere_surface(1.0)
        elements = quadrature.curved_elements(surface, order=4, subdivisions=(1, 1))
        volume = math.fsum(element.volume_contribution for element in elements)
        exact = 4.0 / 3.0 * math.pi
        assert abs(volume - exact) / exact > 1e-6  # honest about where the gate sits


# ----------------------------------------------------------------------------------
# J4 -- ellipsoid Cb and wetted surface against elliptic integrals
# ----------------------------------------------------------------------------------
class TestJ4EllipsoidCrossValidation:
    """Three routes to Cb and S: NURBS quadrature, STL facets, voxel grid, vs elliptic."""

    def test_geometry_gate_precedes_the_calculations(self) -> None:
        """Cb/L/B/T mutual check to 1% -- a failed geometry gate must stop the Cb claim."""
        a, b, c = ELLIPSOID
        hull = geom.half_ellipsoid_patch(a, b, c)
        bounds = patch_bounds(hull)
        reference = geom.hydrostatics_reference(a, b, c)
        for key, expected in (("L", 2 * a), ("B", 2 * b), ("T", c)):
            assert abs(bounds[key] - expected) / expected < 0.01, (key, bounds[key], expected)
        assert abs(reference["Cb"] - math.pi / 6.0) <= 1e-12

    def test_wetted_surface_nurbs_vs_elliptic(self) -> None:
        a, b, c = ELLIPSOID
        hull = geom.half_ellipsoid_patch(a, b, c)
        area = quadrature.surface_area(hull, order=16)
        reference = geom.hydrostatics_reference(a, b, c)
        relative = abs(area - reference["wetted_surface"]) / reference["wetted_surface"]
        assert relative <= 5e-3, relative
        assert relative <= 1e-9, relative  # the NURBS route is exact to round-off, not 0.5%

    def test_wetted_surface_stl_vs_elliptic(self) -> None:
        a, b, c = ELLIPSOID
        hull = geom.half_ellipsoid_patch(a, b, c)
        area, _ = quadrature.faceted_area_and_volume(hull, 96)
        reference = geom.hydrostatics_reference(a, b, c)["wetted_surface"]
        assert abs(area - reference) / reference <= 5e-3

    def test_cb_voxel_vs_elliptic(self) -> None:
        from tensorcad.geometry_field import EllipsoidField
        from tensorcad.voxel import to_voxel_grid

        a, b, c = ELLIPSOID
        grid = to_voxel_grid(EllipsoidField((a, b, c), (0.0, 0.0, 0.0)), resolution=96, padding=2)
        half_volume = 0.5 * grid.volume()  # the grid is symmetric about z = 0
        cb = half_volume / (2 * a * 2 * b * c)
        reference = geom.hydrostatics_reference(a, b, c)["Cb"]
        assert abs(cb - reference) / reference <= 5e-3

    def test_cb_and_volume_from_the_nurbs_patch(self) -> None:
        a, b, c = ELLIPSOID
        hull = geom.half_ellipsoid_patch(a, b, c)
        volume = quadrature.surface_volume(hull, order=16)
        reference = geom.hydrostatics_reference(a, b, c)
        assert abs(volume - reference["volume"]) / reference["volume"] <= 1e-9
        cb = volume / (2 * a * 2 * b * c)
        assert abs(cb - reference["Cb"]) / reference["Cb"] <= 1e-9

    def test_elliptic_integrals_agree_across_their_three_routes(self) -> None:
        """AGM, Carlson and direct quadrature of the defining integral must agree."""
        for k in (0.0, 0.3, 0.7, 0.95):
            nodes, weights = quadrature.gauss_rule(40)
            # E(k) = int_0^{pi/2} sqrt(1 - k^2 sin^2) dtheta; K(k) = int dtheta / sqrt(...)
            half_pi = 0.5 * math.pi
            second = 0.0
            first = 0.0
            for node, weight in zip(nodes, weights, strict=True):
                theta = 0.5 * half_pi * (node + 1.0)
                scale = 0.5 * half_pi
                radicand = 1.0 - k * k * math.sin(theta) ** 2
                second += weight * scale * math.sqrt(radicand)
                first += weight * scale / math.sqrt(radicand)
            assert elliptic.ellipe_agm(k) == pytest.approx(second, rel=1e-12)
            assert elliptic.ellipe(k) == pytest.approx(second, rel=1e-10)
            assert elliptic.ellipk_agm(k) == pytest.approx(first, rel=1e-12)
            assert elliptic.ellipk(k) == pytest.approx(first, rel=1e-10)

    def test_triaxial_formula_degenerates_to_the_spheroid_formula(self) -> None:
        for polar in (0.4, 0.9):
            assert elliptic.triaxial_ellipsoid_area(1.0, 1.0, polar) == pytest.approx(
                elliptic.spheroid_area(1.0, polar), rel=1e-12
            )
        assert elliptic.spheroid_area(1.0, 1.0) == pytest.approx(4.0 * math.pi, rel=1e-15)


# ----------------------------------------------------------------------------------
# J5 -- the element contract
# ----------------------------------------------------------------------------------
@pytest.fixture(scope="module")
def elements() -> tuple[quadrature.CurvedElement, ...]:
    """One patch of curved elements, shared by the contract gates (they never mutate it)."""
    return quadrature.curved_elements(sphere_surface(1.0), order=4, subdivisions=(1, 1))


class TestJ5Contract:
    """Write, read, recompute: the arrays must survive three formats bit for bit."""

    def test_payload_header_states_the_contract(self, elements) -> None:
        payload = contract.build_payload(elements, provenance={"generator": "tests"})
        header = payload["header"]
        assert header["schema"] == contract.ELEMENTS_SCHEMA
        assert header["degree"] == [2, 2]
        assert header["quadrature_order"] == 4
        assert header["rows"] == 6 * 16
        assert header["content_hash"] == contract.content_hash(payload["arrays"], {
            key: value for key, value in header.items() if key != "content_hash"
        })

    @pytest.mark.parametrize("suffix", [".npz", ".h5", ".json"])
    def test_roundtrip_is_bit_identical(self, elements, tmp_path, suffix: str) -> None:
        pytest.importorskip("numpy")
        if suffix in (".h5",):
            pytest.importorskip("h5py")
        payload = contract.build_payload(elements, provenance={"generator": "tests"})
        target = tmp_path / f"elements{suffix}"
        contract.write_contract(target, elements, provenance={"generator": "tests"})
        report = contract.verify_roundtrip(target, payload)
        assert report["all_bit_identical"], report["fields"]
        assert report["hash_match"]
        assert all(entry["max_abs_difference"] == 0.0 for entry in report["fields"].values())
        assert report["recomputed_max_abs_difference"] <= 1e-14

    def test_recomputed_integrals_match_the_direct_ones(self, elements) -> None:
        payload = contract.build_payload(elements)
        recomputed = contract.recompute_integrals(payload["arrays"])
        area = math.fsum(element.area for element in elements)
        volume = math.fsum(element.volume_contribution for element in elements)
        assert abs(recomputed["area"] - area) <= 1e-14 * area
        assert abs(recomputed["volume"] - volume) <= 1e-14 * abs(volume)

    def test_npz_bytes_are_reproducible(self, elements, tmp_path) -> None:
        pytest.importorskip("numpy")
        first = tmp_path / "a.npz"
        second = tmp_path / "b.npz"
        contract.write_contract(first, elements, provenance={"generator": "tests"})
        contract.write_contract(second, elements, provenance={"generator": "tests"})
        assert first.read_bytes() == second.read_bytes()

    def test_mixed_degrees_are_rejected(self) -> None:
        low = quadrature.curved_elements(bezier_patch(1), order=2)
        high = quadrature.curved_elements(bezier_patch(2), order=2)
        with pytest.raises(ValueError):
            contract.build_payload(tuple(low) + tuple(high))


# ----------------------------------------------------------------------------------
# J6 -- the negative control
# ----------------------------------------------------------------------------------
class TestJ6NegativeControl:
    """A faceted patch must lose the p-order and miss the volume gate -- or J1/J3 are void."""

    def test_faceted_patch_converges_at_second_order_only(self) -> None:
        radius = 1.0
        surface = sphere_surface(radius)
        study = quadrature.tessellation_convergence(
            surface,
            (8, 16, 32, 64),
            exact_area=4.0 * math.pi * radius**2,
            exact_volume=4.0 / 3.0 * math.pi * radius**3,
        )
        assert study["area_fitted_order"] is not None and study["area_fitted_order"] <= 2.05
        assert study["volume_fitted_order"] is not None and study["volume_fitted_order"] <= 2.05
        # ... and therefore it cannot meet the p = 2 gate of J1 (2p - 1 = 3).
        assert study["area_fitted_order"] < 3.0

    def test_faceted_volume_error_breaks_the_j3_gate(self) -> None:
        surface = sphere_surface(1.0)
        study = quadrature.tessellation_convergence(
            surface,
            (8, 16),
            exact_area=4.0 * math.pi,
            exact_volume=4.0 / 3.0 * math.pi,
        )
        errors = [row["volume_relative_error"] for row in study["rows"]]
        assert all(error > 1e-3 for error in errors), errors

    def test_faceted_normals_are_outward_but_the_area_is_short(self) -> None:
        """The faceted shell is inscribed, so it always under-estimates -- no lucky signs."""
        surface = sphere_surface(1.0)
        area, volume = quadrature.faceted_area_and_volume(surface, 32)
        assert area < 4.0 * math.pi
        assert volume < 4.0 / 3.0 * math.pi


# ----------------------------------------------------------------------------------
# committed evidence
# ----------------------------------------------------------------------------------
BENCHMARK = Path(__file__).resolve().parents[1] / "benchmarks" / "verified" / CASE


class TestCommittedEvidence:
    """The committed artifacts must still be what the code computes today."""

    def test_gates_artifact_claims_every_gate(self) -> None:
        payload = json.loads((BENCHMARK / "gates_verified.json").read_text(encoding="utf-8"))
        assert payload["schema"] == "tensorcad.gates/1.0"
        assert payload["passed"] is True
        assert set(payload["gates"]) == {f"J{index}" for index in range(1, 8)}
        assert all(entry["passed"] for entry in payload["gates"].values())

    def test_result_artifact_matches_a_fresh_computation(self) -> None:
        result = json.loads((BENCHMARK / "result.json").read_text(encoding="utf-8"))
        reference = geom.hydrostatics_reference(*ELLIPSOID)
        recorded = result["gates"]["J4"]["reference"]
        assert recorded["wetted_surface"] == pytest.approx(reference["wetted_surface"], rel=1e-13)
        assert recorded["Cb"] == pytest.approx(reference["Cb"], rel=1e-13)
        assert result["gates"]["J4"]["passed"] is True

    def test_contract_artifact_is_readable_and_keeps_its_hash(self) -> None:
        pytest.importorskip("numpy")
        artifact = BENCHMARK / "artifacts" / "contract_sphere_r1.npz"
        payload = contract.read_contract(artifact)
        assert payload["header"]["rows"] == 24 * 16
        recomputed = contract.recompute_integrals(payload["arrays"])
        assert abs(recomputed["area"] - 4.0 * math.pi) / (4.0 * math.pi) <= 1e-6
        expected = contract.content_hash(
            payload["arrays"],
            {key: value for key, value in payload["header"].items() if key != "content_hash"},
        )
        assert expected == payload["header"]["content_hash"]


# ----------------------------------------------------------------------------------
# J7 -- CI budget
# ----------------------------------------------------------------------------------
def test_j7_ci_budget() -> None:
    """This module must finish inside the 120 s CI budget (measured, not assumed)."""
    elapsed = time.perf_counter() - _TIMING["first_test"]
    assert elapsed <= 120.0, f"high-order geometry suite took {elapsed:.1f} s"


def test_module_imports_only_the_standard_library() -> None:
    """The core must import without numpy: a solver host may have nothing installed."""
    import subprocess
    import sys

    code = (
        "import sys\n"
        "class Block:\n"
        "    def find_module(self, name, path=None):\n"
        "        if name.split('.')[0] in {'numpy', 'h5py', 'scipy'}:\n"
        "            raise ImportError('blocked: ' + name)\n"
        "        return None\n"
        "sys.meta_path.insert(0, Block())\n"
        "import tensorcad.highorder as h\n"
        "print(h.quadrature.gauss_rule(2)[0])\n"
    )
    completed = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=False
    )
    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip().startswith("(")


def test_report_numbers_are_reproducible() -> None:
    """Spot-check the numbers quoted in the mission report (guards silent drift)."""
    a, b, c = ELLIPSOID
    reference = geom.hydrostatics_reference(a, b, c)
    assert reference["Cb"] == pytest.approx(math.pi / 6.0, rel=1e-14)
    assert reference["wetted_surface"] == pytest.approx(
        elliptic.triaxial_ellipsoid_area(a, b, c) / 2.0, rel=1e-15
    )
    surface = sphere_surface(1.0)
    elements = quadrature.curved_elements(surface, order=4, subdivisions=(2, 2))
    volume = math.fsum(element.volume_contribution for element in elements)
    assert abs(volume - 4.0 / 3.0 * math.pi) / (4.0 / 3.0 * math.pi) == pytest.approx(
        2.06e-7, rel=0.05
    )
    payload = contract.build_payload(elements, provenance={"case": "sphere_r1"})
    assert payload["header"]["rows"] == 24 * 16
    assert len(json.dumps(payload["header"])) > 0  # header is JSON-serialisable as promised
