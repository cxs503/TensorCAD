"""The package's public surface must stay importable and stay standard-library only."""
from __future__ import annotations

import importlib
import subprocess
import sys

import pytest

import tensorcad

EXPECTED_API = {
    "BoundingBox", "Check", "CylinderField", "EllipsoidField", "FieldError", "GeometryField",
    "HullModel", "HullParameters", "Hydrostatics", "Mesh", "MeshBuilder", "ParameterSpec",
    "PlaneField", "RatioSpec", "SampledField", "SphereField", "ValidationReport", "build_mesh",
    "compute_hydrostatics", "estimate_discretisation_error", "export_hull", "generate_hull",
    "inspect_geometry", "inspect_mesh_quality", "read_metadata", "read_obj", "surface_point",
    "to_mesh", "validate_hull", "write_metadata", "write_offsets_csv", "write_obj", "write_stl",
}


def test_declared_api_is_importable() -> None:
    assert set(tensorcad.__all__) == EXPECTED_API
    for name in tensorcad.__all__:
        assert getattr(tensorcad, name) is not None
    assert tensorcad.__version__ == "0.4.0"


def test_layered_modules_are_all_importable() -> None:
    for module in ("geometry", "geometry_field", "mesh", "parameters", "validation",
                   "hydrostatics", "export", "hull", "cli"):
        assert importlib.import_module(f"tensorcad.{module}") is not None


def test_core_import_pulls_in_no_third_party_geometry_stack() -> None:
    """The stdlib-only promise is an API guarantee, so it is tested, not asserted in a README."""
    for name in list(sys.modules):
        if name.split(".")[0] in {"OCP", "cadquery", "gmsh", "numpy", "scipy", "trimesh", "open3d"}:
            del sys.modules[name]
    importlib.reload(tensorcad)
    loaded = {name.split(".")[0] for name in sys.modules}
    assert not loaded & {"OCP", "cadquery", "gmsh", "numpy", "scipy", "trimesh", "open3d"}


def test_adapters_package_is_an_empty_boundary_placeholder() -> None:
    adapters = importlib.import_module("tensorcad.adapters")
    assert adapters.__all__ == []
    assert "core package" in (adapters.__doc__ or "").lower()


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
