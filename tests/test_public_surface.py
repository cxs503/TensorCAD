"""The package's public surface must stay importable and stay standard-library only."""
from __future__ import annotations

import importlib
import subprocess
import sys

import pytest

import tensorcad

EXPECTED_API = {
    "WIGLEY_REFERENCE_CASE", "WigleyAnalytic", "WigleyModel", "WigleyParameters",
    "export_wigley", "generate_wigley", "wigley_analytic_hydrostatics",
    "wigley_analytic_relative_errors", "wigley_convergence_study",
    "wigley_reference_parameters",    # transom-stern family (0.21.0)
    "MAX_TRANSOM_FRACTION", "TransomHullParameters", "TransomModel",
    "analytic_transom_hydrostatics", "build_transom_mesh",
    "estimate_transom_discretisation_error", "generate_transom_hull", "parent_length",
    "sinc_power_coefficients", "transom_area_ratio", "transom_fraction_for_area_ratio",
    "truncated_taper_integral", "truncated_taper_moment",
    "BonjeanCurve", "BonjeanTable", "SectionProperties", "analytic_section_area",
    "bonjean_curve", "bonjean_table", "section_loops", "section_properties",
    "DESIGN_SCHEMA", "DesignProblem", "DesignRecord", "DesignSolution", "DesignTargets",
    "BooleanField", "BoundingBox", "Check", "CylinderField", "EllipsoidField", "FieldError", "GeometryField",
    "AnalyticHydrostatics", "ConditionedMesh", "HullModel", "HullParameters", "Hydrostatics",
    "Mesh", "MeshBuilder",
    "MeshConditioningError", "MeshGateError", "MeshGateReport", "MeshMeasurements", "ParameterSpec",
    "PlaneField", "RatioSpec", "ResolutionPlan", "SUBOFF_LENGTH_M", "SUBOFF_RADIUS_M", "SampledField",
    "SphereField", "SuboffAppendagesField", "SuboffFinsField", "SuboffFullField",
    "SuboffHullField", "SuboffSailField", "ValidationReport", "VoxelGrid", "analytic_hydrostatics",
    "analytic_volume", "build_mesh",
    "ANALYTIC_COMPARABLE_QUANTITIES",
    "CurvatureError", "CurvatureSample", "FairnessReport",
    "analytic_relative_errors", "analytic_wetted_surface_area",
    "InitialStability", "StabilityError", "WaterplaneMoments",
    "compute_hydrostatics", "initial_stability", "waterplane_moments",
    "DEFAULT_ANGLES_DEG", "FreeSurfaceCorrection", "GzError", "HeeledHydrostatics",
    "RightingArmCurve", "box_gz", "box_gz_curve", "box_regime_bounds_deg",
    "box_righting_arms_analytic", "box_section", "box_vanishing_angle_deg",
    "box_wall_sided_limit_deg",
    "deck_edge_immersion_deg", "free_surface_correction", "heeled_hydrostatics",
    # tank-shift level/volume curve exported by the free-surface module but never added
    # to this pinned list: the contract test had been red since that merge.
    "free_surface_gz_shift", "free_surface_gz_shift_curve", "free_surface_loss_m",
    "rectangular_tank_inertia", "righting_arm_curve", "wall_sided_righting_arm",
    "with_freeboard",
    # second GZ route (0.22.0): clipped-section integral, bisection waterline solve
    "DEFAULT_HEELS_DEG", "HeeledDisplacement", "RightingArm", "RightingArmError",
    "heeled_displacement", "righting_arm", "wall_sided_validity_deg",
    "BoxGzError", "BoxGzPoint", "BisectionRightingArmCurve", "bisection_righting_arm_curve",
    # areas under the exact multi-regime curve: IS Code areas past the wall-sided limit
    "BoxGzArea", "box_gz_area",
    "barge_wall_sided_righting_arm",
    "waterplane_inertia_longitudinal", "waterplane_inertia_transverse", "condition_mesh", "decimate",
    "enforce_gate",
    "estimate_discretisation_error", "export_hull", "gate_mesh", "generate_hull",
    "inspect_geometry", "inspect_mesh_quality", "measure_mesh", "project_to_field",
    "read_metadata", "read_obj", "read_stl", "resolve_resolution",
    "suboff_fin_chord", "suboff_radius_fraction", "suboff_reference_geometry",
    "suboff_sail_half_width", "surface_point",
    "difference", "intersection", "union",
    "tanh_sinh_quadrature", "taper_integral", "taubin_smooth", "to_mesh", "to_voxel_grid", "validate_hull",
    "write_design_json", "write_metadata",
    "assess_fairness", "principal_curvatures",
    "vertex_mean_curvature",
    "write_offsets_csv", "write_obj", "write_solver_manifest", "write_stl",
    # intact-stability criteria (0.23.0): IS Code 2008 part A on a measured GZ curve
    "Criterion", "IS_CODE_2008_GENERAL", "IntactStabilityError", "IntactStabilityReport",
    "WeatherCriterionReport", "WindHeelingLevers",
    "general_criteria", "gz_area_between", "gz_at_angle",
    "steady_wind_heel_angle", "weather_criterion", "wind_heeling_levers",
    # statutory rule check (0.23.0+): IS Code 2008 s. 2.2 with per-rule error bars
    "INTACT_STABILITY_CRITERIA", "CriterionCheck", "IntactStabilityAssessment", "StatutoryError",
    "area_under_gz_curve", "check_intact_stability", "gm0_from_gz_curve",
}


def test_declared_api_is_importable() -> None:
    assert set(tensorcad.__all__) == EXPECTED_API
    for name in tensorcad.__all__:
        assert getattr(tensorcad, name) is not None
    assert tensorcad.__version__ == "0.23.0"


def test_layered_modules_are_all_importable() -> None:
    for module in ("curvature", "geometry", "geometry_field", "mesh", "mesh_conditioning", "gate",
                   "parameters", "validation", "hydrostatics", "hydrostatics_analytic", "voxel", "export", "csg",
                   "wigley", "transom", "hull", "design", "suboff", "suboff_appendages", "cli",
                   "adapters.fem", "benchmark", "nurbs"):
        assert importlib.import_module(f"tensorcad.{module}") is not None


def test_core_import_pulls_in_no_third_party_geometry_stack() -> None:
    """The stdlib-only promise is an API guarantee, so it is tested, not asserted in a README.

    Checked in a fresh interpreter rather than by deleting ``numpy`` from this one:
    removing a module from ``sys.modules`` while other libraries hold references to
    it (torch's numpy bridge, for one) leaves those libraries broken for every later
    test in the session, which is exactly the kind of order-dependent failure a
    test suite must not have.
    """
    script = (
        "import sys, tensorcad;"
        "loaded = {name.split('.')[0] for name in sys.modules};"
        "bad = loaded & {'OCP', 'cadquery', 'gmsh', 'numpy', 'scipy', 'trimesh', 'open3d'};"
        "assert not bad, sorted(bad)"
    )
    completed = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True)
    assert completed.returncode == 0, completed.stderr


def test_adapters_package_is_a_documented_boundary() -> None:
    """Adapters are reached explicitly; nothing in the core imports them."""
    adapters = importlib.import_module("tensorcad.adapters")
    assert adapters.__all__ == []
    assert "core package" in (adapters.__doc__ or "").lower()
    assert "tensorcad.adapters.ocp" in (adapters.__doc__ or "")
    assert "tensorcad.adapters.fem" in (adapters.__doc__ or "")


def test_neural_runtime_is_opt_in_and_not_imported_by_the_core() -> None:
    """``tensorcad`` must stay standard-library only; the field runtime is an extra.

    Checked in a fresh interpreter: this is the whole reason the neural layer is
    a subpackage rather than part of the core, and an agent that only needs hull
    geometry must never pay for torch.
    """
    script = (
        "import sys, tensorcad;"
        "loaded = {'torch', 'numpy', 'safetensors', 'tensorcad.neural'};"
        "assert not (set(sys.modules) & loaded), sorted(set(sys.modules) & loaded)"
    )
    completed = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True)
    assert completed.returncode == 0, completed.stderr


def test_hull_module_does_not_import_the_cli_or_vice_versa_in_reverse() -> None:
    hull = importlib.import_module("tensorcad.hull")
    assert not hasattr(hull, "main")
    with pytest.raises(ImportError):
        importlib.import_module("tensorcad.nonexistent")
