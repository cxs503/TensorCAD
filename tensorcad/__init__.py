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

from .export import (
    export_hull,
    read_metadata,
    read_obj,
    write_metadata,
    write_obj,
    write_offsets_csv,
    write_stl,
)
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
    HullModel,
    HullParameters,
    build_mesh,
    estimate_discretisation_error,
    generate_hull,
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

__all__ = [
    "SUBOFF_LENGTH_M",
    "SUBOFF_RADIUS_M",
    "BoundingBox",
    "Check",
    "ConditionedMesh",
    "CylinderField",
    "EllipsoidField",
    "FieldError",
    "GeometryField",
    "HullModel",
    "HullParameters",
    "Hydrostatics",
    "Mesh",
    "MeshBuilder",
    "MeshConditioningError",
    "MeshMeasurements",
    "ParameterSpec",
    "PlaneField",
    "RatioSpec",
    "SampledField",
    "SphereField",
    "SuboffAppendagesField",
    "SuboffFinsField",
    "SuboffFullField",
    "SuboffHullField",
    "SuboffSailField",
    "ValidationReport",
    "build_mesh",
    "compute_hydrostatics",
    "condition_mesh",
    "decimate",
    "estimate_discretisation_error",
    "export_hull",
    "generate_hull",
    "inspect_geometry",
    "inspect_mesh_quality",
    "measure_mesh",
    "project_to_field",
    "read_metadata",
    "read_obj",
    "suboff_fin_chord",
    "suboff_radius_fraction",
    "suboff_reference_geometry",
    "suboff_sail_half_width",
    "surface_point",
    "taubin_smooth",
    "to_mesh",
    "validate_hull",
    "write_metadata",
    "write_obj",
    "write_offsets_csv",
    "write_stl",
]
__version__ = "0.7.0"
