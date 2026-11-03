"""TensorCAD: deterministic, agent-callable engineering geometry.

Public API at a glance::

    from tensorcad import HullParameters, generate_hull, export_hull

    model = generate_hull(HullParameters(length=12, beam=2.4, draft=0.9))
    model.hydrostatics.displacement_volume          # integrated from the mesh
    model.validation.warnings                       # what to look at next
    export_hull(model, "results/hull", convergence=True)

Layering: :mod:`tensorcad.geometry` and :mod:`tensorcad.mesh` are generic
contracts and algorithms, :mod:`tensorcad.parameters` and
:mod:`tensorcad.validation` are the checking layers,
:mod:`tensorcad.hydrostatics` integrates properties off a surface, and
:mod:`tensorcad.hull` is the one ship-specific model. Optional solver bridges
live in :mod:`tensorcad.adapters`.
"""

from .asymmetric_hull import (
    AsymmetricHullModel,
    AsymmetricHullParameters,
)
from .asymmetric_hull import (
    analytic_hydrostatics as analytic_asymmetric_hydrostatics,
)
from .asymmetric_hull import (
    generate_hull as generate_asymmetric_hull,
)
from .curvature import (
    CurvatureError,
    CurvatureSample,
    FairnessReport,
    assess_fairness,
    principal_curvatures,
    vertex_mean_curvature,
)
from .design import (
    DESIGN_SCHEMA,
    DesignProblem,
    DesignRecord,
    DesignSolution,
    DesignTargets,
    write_design_json,
    write_solver_manifest,
)
from .export import (
    export_hull,
    read_metadata,
    read_obj,
    read_stl,
    write_metadata,
    write_obj,
    write_offsets_csv,
    write_stl,
)
from .gate import MeshGateError, MeshGateReport, enforce_gate, gate_mesh
from .geometry import BoundingBox, Mesh, MeshBuilder
from .geometry_field import (
    CylinderField,
    EllipsoidField,
    FieldError,
    GeometryField,
    PlaneField,
    SampledField,
    SphereField,
    surface_point,
    to_mesh,
)
from .hull import (
    ANALYTIC_COMPARABLE_QUANTITIES,
    HullModel,
    HullParameters,
    ResolutionPlan,
    analytic_hydrostatics,
    analytic_relative_errors,
    analytic_volume,
    analytic_wetted_surface_area,
    build_mesh,
    estimate_discretisation_error,
    generate_hull,
    resolve_resolution,
    validate_hull,
)
from .hydrostatics import Hydrostatics, compute_hydrostatics
from .mesh_conditioning import (
    ConditionedMesh,
    MeshConditioningError,
    MeshMeasurements,
    condition_mesh,
    decimate,
    measure_mesh,
    project_to_field,
    taubin_smooth,
)
from .parameters import ParameterSpec, RatioSpec
from .suboff import (
    SUBOFF_LENGTH_M,
    SUBOFF_RADIUS_M,
    SuboffHullField,
    suboff_radius_fraction,
    suboff_reference_geometry,
)
from .suboff_appendages import (
    SuboffAppendagesField,
    SuboffFinsField,
    SuboffFullField,
    SuboffSailField,
    suboff_fin_chord,
    suboff_sail_half_width,
)
from .validation import (
    Check,
    ValidationReport,
    inspect_geometry,
    inspect_mesh_quality,
)
from .voxel import VoxelGrid, to_voxel_grid

__all__ = [
    "ANALYTIC_COMPARABLE_QUANTITIES",
    "DESIGN_SCHEMA",
    "SUBOFF_LENGTH_M",
    "SUBOFF_RADIUS_M",
    "AsymmetricHullModel",
    "AsymmetricHullParameters",
    "BoundingBox",
    "Check",
    "ConditionedMesh",
    "CurvatureError",
    "CurvatureSample",
    "CylinderField",
    "DesignProblem",
    "DesignRecord",
    "DesignSolution",
    "DesignTargets",
    "EllipsoidField",
    "FairnessReport",
    "FieldError",
    "GeometryField",
    "HullModel",
    "HullParameters",
    "Hydrostatics",
    "Mesh",
    "MeshBuilder",
    "MeshConditioningError",
    "MeshGateError",
    "MeshGateReport",
    "MeshMeasurements",
    "ParameterSpec",
    "PlaneField",
    "RatioSpec",
    "ResolutionPlan",
    "SampledField",
    "SphereField",
    "SuboffAppendagesField",
    "SuboffFinsField",
    "SuboffFullField",
    "SuboffHullField",
    "SuboffSailField",
    "ValidationReport",
    "VoxelGrid",
    "analytic_asymmetric_hydrostatics",
    "analytic_hydrostatics",
    "analytic_relative_errors",
    "analytic_volume",
    "analytic_wetted_surface_area",
    "assess_fairness",
    "build_mesh",
    "compute_hydrostatics",
    "condition_mesh",
    "decimate",
    "enforce_gate",
    "estimate_discretisation_error",
    "export_hull",
    "gate_mesh",
    "generate_asymmetric_hull",
    "generate_hull",
    "inspect_geometry",
    "inspect_mesh_quality",
    "measure_mesh",
    "principal_curvatures",
    "project_to_field",
    "read_metadata",
    "read_obj",
    "read_stl",
    "resolve_resolution",
    "suboff_fin_chord",
    "suboff_radius_fraction",
    "suboff_reference_geometry",
    "suboff_sail_half_width",
    "surface_point",
    "taubin_smooth",
    "to_mesh",
    "to_voxel_grid",
    "validate_hull",
    "vertex_mean_curvature",
    "write_design_json",
    "write_metadata",
    "write_obj",
    "write_offsets_csv",
    "write_solver_manifest",
    "write_stl",
]
__version__ = "0.13.0"
