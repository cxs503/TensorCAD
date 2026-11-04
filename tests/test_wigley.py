"""Pins the Wigley hull form: closed forms, sampled convergence, exports.

Two kinds of assertion live here and they fail for different reasons:

* the **club of exact identities** (``Cb = 4/9``, ``V = 4LBT/9``, ``VCB = 5T/8``)
  would fail if the algebra or the integrator drifted;
* the **convergence claims** (second-order approach from below, offsets on the
  analytic surface) would fail if the generator stopped sampling the surface it
  documents.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

import tensorcad
from tensorcad import mesh as meshlib
from tensorcad.wigley import (
    ANALYTIC_COMPARABLE_QUANTITIES,
    REFERENCE_CASE,
    WigleyParameters,
    analytic_block_coefficient,
    analytic_maximum_section_area,
    analytic_midship_coefficient,
    analytic_prismatic_coefficient,
    analytic_relative_errors,
    analytic_volume,
    analytic_waterplane_area,
    convergence_study,
    export_wigley,
    generate_wigley,
    half_breadth,
    reference_parameters,
    wetted_area_evidence,
    wigley_analytic_hydrostatics,
)

REPO_ROOT = Path(__file__).resolve().parents[1]


# --------------------------------------------------------------- closed forms


def test_classic_wigley_coefficients_are_the_published_rationals() -> None:
    """Cb 4/9, Cp 2/3, Cm 2/3, Cwp 2/3 -- the numbers every Wigley paper quotes."""
    assert analytic_block_coefficient() == pytest.approx(4.0 / 9.0, rel=1e-15)
    assert analytic_prismatic_coefficient() == pytest.approx(2.0 / 3.0, rel=1e-15)
    assert analytic_midship_coefficient() == pytest.approx(2.0 / 3.0, rel=1e-15)
    # Cb = Cp * Cm, the identity that ties the three coefficients together.
    assert analytic_block_coefficient() / analytic_prismatic_coefficient() == pytest.approx(
        analytic_midship_coefficient(), rel=1e-15
    )
    reference = wigley_analytic_hydrostatics()
    assert reference.block_coefficient == pytest.approx(4.0 / 9.0, rel=1e-15)
    assert reference.prismatic_coefficient == pytest.approx(2.0 / 3.0, rel=1e-15)
    assert reference.midship_coefficient == pytest.approx(2.0 / 3.0, rel=1e-15)
    assert reference.waterplane_coefficient == pytest.approx(2.0 / 3.0, rel=1e-15)


def test_closed_form_volume_and_areas() -> None:
    """V = 4LBT/9, Awp = 2LB/3, Am = 2BT/3 for L=10, B=1, T=0.625."""
    length, beam, draft = 10.0, 1.0, 0.625
    assert analytic_volume(length, beam, draft) == pytest.approx(4.0 * length * beam * draft / 9.0, rel=1e-15)
    assert analytic_volume(length, beam, draft) == pytest.approx(2.7777777777777, rel=1e-12)
    assert analytic_waterplane_area(length, beam) == pytest.approx(2.0 * length * beam / 3.0, rel=1e-15)
    assert analytic_maximum_section_area(beam, draft) == pytest.approx(2.0 * beam * draft / 3.0, rel=1e-15)


def test_closed_forms_hold_for_general_exponents() -> None:
    """The identities are algebraic, so they survive n != m != 2 as well."""
    for n, m in ((2.0, 2.0), (3.0, 1.5), (1.0, 1.0), (2.5, 4.0)):
        parameters = reference_parameters(longitudinal_exponent=n, vertical_exponent=m)
        reference = wigley_analytic_hydrostatics(parameters)
        box = parameters.length * parameters.beam * parameters.draft
        assert reference.displacement_volume == pytest.approx(
            analytic_block_coefficient(n, m) * box, rel=1e-15
        )
        assert reference.waterplane_area == pytest.approx(
            parameters.length * parameters.beam * n / (n + 1.0), rel=1e-15
        )
        assert reference.maximum_section_area == pytest.approx(
            parameters.beam * parameters.draft * m / (m + 1.0), rel=1e-15
        )
        assert reference.prismatic_coefficient * reference.midship_coefficient == pytest.approx(
            reference.block_coefficient, rel=1e-12
        )


def test_vertical_centre_of_buoyancy_is_five_eighths_of_draft() -> None:
    """VCB = T(1 - (m+1)/(2(m+2))) above the keel: 5T/8 at m = 2, i.e. 3T/8 below the WL."""
    draft = 0.625
    reference = wigley_analytic_hydrostatics()
    assert reference.vertical_centre_buoyancy + draft == pytest.approx(5.0 * draft / 8.0, rel=1e-15)
    assert reference.vertical_centre_buoyancy == pytest.approx(-3.0 * draft / 8.0, rel=1e-15)
    # m = 1 puts the centroid at T/3 below the waterline (triangle in the section).
    assert wigley_analytic_hydrostatics(
        reference_parameters(vertical_exponent=1.0)
    ).vertical_centre_buoyancy == pytest.approx(-draft / 3.0, rel=1e-15)


def test_analytic_longitudinal_centres_are_midship() -> None:
    """The form is fore-aft symmetric, so LCB = LCF = L/2 identically."""
    reference = wigley_analytic_hydrostatics(reference_parameters(length=12.5))
    assert reference.longitudinal_centre_buoyancy == pytest.approx(6.25, rel=1e-15)
    assert reference.longitudinal_centre_flotation == pytest.approx(6.25, rel=1e-15)


def test_wetted_area_quadrature_is_self_converged() -> None:
    """No primitive exists, so the reference has to prove itself by refinement."""
    evidence = wetted_area_evidence()
    assert evidence["relative_change_fine_to_finer"] < 1e-12
    assert evidence["relative_change_coarse_to_fine"] < 1e-9
    # A crude sanity band: the flank is longer than 2*L*T and shorter than the
    # rectangle that would hold a box of the same draft and beam.
    parameters = reference_parameters()
    lower = 2.0 * parameters.length * parameters.draft
    assert lower < wigley_analytic_hydrostatics(parameters).wetted_surface_area < 4.0 * lower


# ------------------------------------------------------------------- geometry


def test_half_breadth_closes_at_both_ends_and_the_keel() -> None:
    parameters = reference_parameters()
    assert half_breadth(parameters, 0.0, 0.0) == 0.0
    assert half_breadth(parameters, parameters.length, 0.0) == 0.0
    assert half_breadth(parameters, 0.5 * parameters.length, parameters.draft) == 0.0
    assert half_breadth(parameters, 0.5 * parameters.length, 0.0) == pytest.approx(
        0.5 * parameters.beam, rel=1e-15
    )


def test_generated_surface_is_a_closed_oriented_solid() -> None:
    model = generate_wigley()
    assert meshlib.is_watertight(model.mesh)
    assert meshlib.is_edge_manifold(model.mesh)
    assert meshlib.zero_area_faces(model.mesh) == ()
    assert meshlib.total_signed_volume(model.mesh) > 0.0  # outward-oriented
    box = model.mesh.bounding_box()
    assert box.maximum[0] - box.minimum[0] == pytest.approx(model.parameters.length, rel=1e-12)
    assert box.maximum[1] - box.minimum[1] == pytest.approx(model.parameters.beam, rel=1e-12)
    assert box.minimum[2] == pytest.approx(-model.parameters.draft, rel=1e-12)
    assert box.maximum[2] == 0.0  # exactly the design waterline
    assert model.validation.passed
    assert model.validation.errors == []


def test_euler_characteristic_of_the_closed_surface() -> None:
    """V - E + F = 2 for the lofted and capped solid (undirected edges)."""
    mesh = generate_wigley().mesh
    edges = {
        (min(a, b), max(a, b))
        for face in mesh.faces
        for a, b in ((face[0], face[1]), (face[1], face[2]), (face[2], face[0]))
    }
    assert len(mesh.vertices) - len(edges) + len(mesh.faces) == 2


def test_generation_is_deterministic() -> None:
    first = generate_wigley(WigleyParameters(nx=21, nz=9))
    second = generate_wigley(WigleyParameters(nx=21, nz=9))
    assert first.mesh.vertices == second.mesh.vertices
    assert first.mesh.faces == second.mesh.faces


def test_waterline_offsets_sit_on_the_analytic_surface() -> None:
    """The loft interpolates only along x, so offsets are exact at sampled stations."""
    parameters = reference_parameters()
    model = generate_wigley(parameters)
    assert len(model.offsets) == parameters.nx
    assert model.offsets[0][1] == 0.0 and model.offsets[-1][1] == 0.0
    for x, width in model.offsets:
        assert width == pytest.approx(2.0 * half_breadth(parameters, x, 0.0), abs=1e-15)
    assert max(width for _, width in model.offsets) == pytest.approx(parameters.beam, rel=1e-15)


def test_even_nx_reports_the_unachievable_beam_as_a_warning() -> None:
    """No station lands on midship, so the realised beam is short -- and says so."""
    model = generate_wigley(reference_parameters(nx=40))
    assert model.validation.passed  # a warning, not an error
    assert any("requested_beam_not_achievable" in warning for warning in model.validation.warnings)
    assert max(width for _, width in model.offsets) < model.parameters.beam


# ------------------------------------------------------ sampled vs closed forms


def test_mesh_hydrostatics_track_the_closed_forms() -> None:
    parameters = reference_parameters()
    model = generate_wigley(parameters)
    reference = wigley_analytic_hydrostatics(parameters)
    hydrostatics = model.hydrostatics
    assert hydrostatics is not None
    # Inscribed chords: low, but by a small and systematic amount.
    for name in ("displacement_volume", "waterplane_area", "wetted_surface_area"):
        measured = float(getattr(hydrostatics, name))
        exact = float(getattr(reference, name))
        assert measured < exact, name
        assert abs(measured / exact - 1.0) < 0.01, name
    assert hydrostatics.block_coefficient == pytest.approx(4.0 / 9.0, rel=5e-3)
    assert hydrostatics.longitudinal_centre_buoyancy == pytest.approx(0.5 * parameters.length, abs=1e-9)
    assert hydrostatics.vertical_centre_buoyancy == pytest.approx(
        reference.vertical_centre_buoyancy, rel=1e-2
    )
    assert hydrostatics.maximum_section_area == pytest.approx(reference.maximum_section_area, rel=1e-2)
    assert hydrostatics.displacement_mass == pytest.approx(
        hydrostatics.displacement_volume * parameters.water_density, rel=1e-12
    )


def test_analytic_relative_errors_payload() -> None:
    report = analytic_relative_errors(reference_parameters(nx=81, nz=21))
    assert report["status"] == "measured"
    assert set(report["relative_error"]) == set(ANALYTIC_COMPARABLE_QUANTITIES)
    assert report["max_abs_relative_error"] < 1e-3
    # Extensive quantities are low (inscribed chords) and the symmetric-family
    # longitudinal centres are exact; VCB is the one quantity whose signed error is
    # positive, because the chords lose more volume near the broad waterline than
    # near the keel and the centroid therefore sits slightly high.
    for name in ("displacement_volume", "waterplane_area", "wetted_surface_area", "maximum_section_area"):
        assert report["relative_error"][name] < 0.0, name
    assert report["relative_error"]["longitudinal_centre_buoyancy"] == 0.0
    assert report["relative_error"]["longitudinal_centre_flotation"] == 0.0
    assert abs(report["relative_error"]["vertical_centre_buoyancy"]) < 1e-3


def test_convergence_is_second_order_and_from_below() -> None:
    study = convergence_study(((21, 9), (41, 13), (81, 21)))
    errors = [row["relative_error"]["displacement_volume"] for row in study["rows"]]
    assert all(value < 0.0 for value in errors)
    assert abs(errors[0]) > abs(errors[1]) > abs(errors[2])
    assert 1.8 <= study["observed_order"]["displacement_volume"] <= 3.0
    # Halving h four times (21 -> 81 levels) must cut the error by ~16x.
    assert abs(errors[0]) / abs(errors[2]) > 8.0


def test_the_smooth_form_beats_the_sine_family_on_the_same_grid() -> None:
    """The claim that motivates the family: finite keel angle converges faster.

    ``sin(pi x/L)**0.7 * sqrt(z/T)`` has an infinite section derivative at the
    keel; at the same 41x13 grid its volume error is ~0.88 %, the Wigley's ~0.24 %.
    """
    from tensorcad.hull import HullParameters
    from tensorcad.hull import analytic_relative_errors as hull_errors

    wigley = analytic_relative_errors(reference_parameters())["relative_error"]["displacement_volume"]
    hull = hull_errors(HullParameters())["relative_error"]["displacement_volume"]
    assert abs(wigley) < 0.5 * abs(hull)


# ---------------------------------------------------------------- parameters


@pytest.mark.parametrize(
    "overrides",
    [
        {"length": 0.0},
        {"beam": -1.0},
        {"draft": 0.0},
        {"nx": 2},
        {"nz": 2},
        {"water_density": 0.0},
        {"longitudinal_exponent": 0.5},
        {"vertical_exponent": 9.0},
        {"longitudinal_exponent": float("nan")},
    ],
)
def test_invalid_parameters_are_rejected(overrides: dict) -> None:
    parameters = WigleyParameters(**overrides)
    errors, _ = parameters.problems()
    assert errors
    with pytest.raises(ValueError):
        generate_wigley(parameters)


def test_reference_case_is_the_ittc_proportions() -> None:
    parameters = reference_parameters()
    assert parameters.length / parameters.beam == pytest.approx(float(REFERENCE_CASE["length_to_beam"]))
    assert parameters.beam / parameters.draft == pytest.approx(float(REFERENCE_CASE["beam_to_draft"]))
    assert parameters.grid_is_symmetric
    assert parameters.problems() == ([], [])


def test_generator_accepts_dimension_triples() -> None:
    from_triple = generate_wigley((12.0, 1.2, 0.75))
    from_parameters = generate_wigley(WigleyParameters(length=12.0, beam=1.2, draft=0.75))
    assert from_triple.mesh.vertices == from_parameters.mesh.vertices


# ------------------------------------------------------------------ artifacts


def test_export_writes_obj_offsets_metadata_and_stl(tmp_path: Path) -> None:
    model = generate_wigley(WigleyParameters(nx=21, nz=9))
    paths = export_wigley(model, tmp_path, stl=True)
    for key in ("mesh", "offsets", "metadata", "stl"):
        assert Path(paths[key]).exists() and Path(paths[key]).stat().st_size > 0
    metadata = json.loads(Path(paths["metadata"]).read_text(encoding="utf-8"))
    assert metadata["hull_form"] == "wigley"
    assert metadata["analytic"]["block_coefficient"] == pytest.approx(4.0 / 9.0, rel=1e-15)
    # 21x9 is the coarsest grid in the study, so the sampled Cb is ~0.6 % low.
    assert metadata["hydrostatics"]["block_coefficient"] == pytest.approx(4.0 / 9.0, rel=2e-2)
    # STL is a triangle soup with no topology: it has to re-weld before it can be
    # called watertight, and the recovered volume must match the mesh of record.
    recovered = meshlib.weld(tensorcad.read_stl(Path(paths["stl"])))
    assert meshlib.is_watertight(recovered)
    assert meshlib.total_signed_volume(recovered) == pytest.approx(
        meshlib.total_signed_volume(model.mesh), rel=1e-9
    )


def test_public_surface_exports_the_wigley_api() -> None:
    for name in (
        "WigleyParameters",
        "WigleyModel",
        "WigleyAnalytic",
        "generate_wigley",
        "wigley_analytic_hydrostatics",
        "wigley_reference_parameters",
        "wigley_convergence_study",
        "wigley_analytic_relative_errors",
        "export_wigley",
        "WIGLEY_REFERENCE_CASE",
    ):
        assert name in tensorcad.__all__, name
        assert hasattr(tensorcad, name), name


def test_committed_results_document_has_no_failed_claim() -> None:
    """``docs/wigley/results.json`` is generated by the example; keep it honest."""
    path = REPO_ROOT / "docs" / "wigley" / "results.json"
    if not path.exists():  # pragma: no cover - artifact may not be checked out
        pytest.skip("docs/wigley/results.json not present")
    data = json.loads(path.read_text(encoding="utf-8"))
    assert data["failed_claims"] == []
    assert data["topology"]["euler_characteristic"] == 2
    assert data["analytic"]["block_coefficient"] == pytest.approx(4.0 / 9.0, rel=1e-15)
    assert data["convergence"]["rows"][-1]["max_abs_relative_error"] < 1e-3
    assert math.isfinite(data["offset_on_surface_max_deviation_m"])
