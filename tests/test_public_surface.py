"""The package's public surface must stay importable and stay standard-library only."""
from __future__ import annotations

import importlib
import subprocess
import sys

import pytest

import tensorcad

EXPECTED_API = {
    "DESIGN_SCHEMA", "DesignProblem", "DesignRecord", "DesignSolution", "DesignTargets",
    "BoundingBox", "Check", "CylinderField", "EllipsoidField", "FieldError", "GeometryField",
    "ConditionedMesh", "HullModel", "HullParameters", "Hydrostatics", "Mesh", "MeshBuilder",
    "MeshConditioningError", "MeshMeasurements", "ParameterSpec",
    "PlaneField", "RatioSpec", "SUBOFF_LENGTH_M", "SUBOFF_RADIUS_M", "SampledField",
    "SphereField", "SuboffAppendagesField", "SuboffFinsField", "SuboffFullField",
    "SuboffHullField", "SuboffSailField", "ValidationReport", "build_mesh",
    "compute_hydrostatics", "condition_mesh", "decimate", "estimate_discretisation_error",
    "export_hull", "generate_hull",
    "inspect_geometry", "inspect_mesh_quality", "measure_mesh", "project_to_field",
    "read_metadata", "read_obj",
    "suboff_fin_chord", "suboff_radius_fraction", "suboff_reference_geometry",
    "suboff_sail_half_width", "surface_point",
    "taubin_smooth", "to_mesh", "validate_hull", "write_design_json", "write_metadata",
    "write_offsets_csv", "write_obj", "write_solver_manifest", "write_stl",
}


def test_declared_api_is_importable() -> None:
    assert set(tensorcad.__all__) == EXPECTED_API
    for name in tensorcad.__all__:
        assert getattr(tensorcad, name) is not None
    assert tensorcad.__version__ == "0.8.0"


def test_layered_modules_are_all_importable() -> None:
    for module in ("geometry", "geometry_field", "mesh", "mesh_conditioning", "parameters",
                   "validation", "hydrostatics", "export", "hull", "design", "suboff",
                   "suboff_appendages", "cli", "adapters.fem"):
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
