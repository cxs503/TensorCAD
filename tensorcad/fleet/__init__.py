"""Baseline-hull fleet: registry, geometry routes, intake gates, run harness.

Entry point (CLI)::

    python -m tensorcad.fleet.cli --out-dir runs/m17335_fleet

Library entry point::

    from tensorcad.fleet import load_manifest, run_fleet
    report = run_fleet("runs/m17335_fleet")

The fleet is six standard hulls (KCS, KVLCC2, ONRT, JBC, DTC, Series 60). Each one is
either built from a stated source and admitted by machine-decidable gates, or refused
with the number that caused the refusal. Nothing in this package claims that a
form-family hull *is* the benchmark ship: the record's provenance and the reference
basis travel with every number into the artefacts.

Geometry routes
---------------
* ``series60_form`` — Series-60 section law (Cb ∈ [0.55, 0.85])
* ``parametric_sections`` — generalized Beta section-area law (any Cb < Cm ≤ 1)
* ``offset_table`` — lofted station × waterline half-breadth table
* ``imported_stl`` — binary/ASCII STL, welded and oriented on intake
* ``iges_field`` — IGES B-spline patches → half-breadth field → watertight solid
"""

from tensorcad.fleet.build import (
    SERIES60_CB_WINDOW,
    BuiltGeometry,
    RouteRefusal,
    build_ship_geometry,
    iges_field_mesh,
    loft_station_table,
    offset_table_mesh,
    parametric_table,
    read_stl,
    sample_offset_table,
    stl_mesh,
)
from tensorcad.fleet.controls import (
    clipped_cb_record,
    duplicate_manifest,
    duplicate_stl_pair,
    placeholder_box,
    wrong_volume_record,
)
from tensorcad.fleet.gates import (
    VOLUME_TOLERANCE,
    TopologyFacts,
    evaluate_gates,
    gate_volume,
    measure,
    placeholder_reasons,
    refused,
)
from tensorcad.fleet.manifest import (
    FleetManifest,
    ManifestError,
    Provenance,
    ShipRecord,
    fleet_sha_index,
    load_manifest,
    sha256_file,
)
from tensorcad.fleet.run import run_fleet

__all__ = [
    "SERIES60_CB_WINDOW",
    "VOLUME_TOLERANCE",
    "BuiltGeometry",
    "FleetManifest",
    "ManifestError",
    "Provenance",
    "RouteRefusal",
    "ShipRecord",
    "TopologyFacts",
    "build_ship_geometry",
    "clipped_cb_record",
    "duplicate_manifest",
    "duplicate_stl_pair",
    "evaluate_gates",
    "fleet_sha_index",
    "gate_volume",
    "iges_field_mesh",
    "load_manifest",
    "loft_station_table",
    "measure",
    "offset_table_mesh",
    "parametric_table",
    "placeholder_box",
    "placeholder_reasons",
    "read_stl",
    "refused",
    "run_fleet",
    "sample_offset_table",
    "sha256_file",
    "stl_mesh",
    "wrong_volume_record",
]
