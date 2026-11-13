"""Mesh quality as a library: metrics, statistics, and a machine gate.

Why this package exists
-----------------------
The repository already measured mesh quality in three places, each with its own
vocabulary: :func:`tensorcad.validation.inspect_mesh_quality` (minimum angle, an aspect
ratio and a sliver share, triangles only), :meth:`tensorcad.volume_mesh.VolumeMesh.quality`
(a scaled Jacobian, hexahedra and tetrahedra) and
:mod:`tensorcad.geomquality.mesh_defects` (facet-scale defect codes). A solver asking
"is this mesh usable?" had to know all three and could not compare them.

This package is the single library behind that question:

* :mod:`tensorcad.meshquality.metrics` -- per-element metrics for triangles,
  quadrilaterals, tetrahedra and hexahedra, with the Gmsh and Verdict definitions and
  an explicit normalisation for every number;
* :mod:`tensorcad.meshquality.report` -- outlier-resistant statistics and a threshold
  gate that returns data, not prose;
* :mod:`tensorcad.meshquality.io` -- a JSON exchange format so a solver can hand a mesh
  over without a mesh library;
* :mod:`tensorcad.meshquality.cli` -- the CI entry point (exit 1 on a failed gate);
* :mod:`tensorcad.meshquality.topology` -- the non-per-element half: watertight,
  edge/vertex-manifold, duplicate faces, winding, outward orientation, Euler
  characteristic, shells, and the OpenFOAM internal-face non-orthogonality;
* :mod:`tensorcad.meshquality.dialects` -- the three aspect-ratio dialects this tree
  speaks (``legacy``, ``edge_form``, ``verdict``), the exact algebra relating them, and
  the threshold transfer that keeps a gate at least as strict as the one it replaces;
* :mod:`tensorcad.meshquality.reference` -- second implementations and external oracles
  (Gmsh, trimesh, and the VTK/Verdict definitions meshmagick uses) that the benchmark
  compares against, so the library is never its own judge.

Cross-validation is the contract, not an extra: see ``docs/mesh_quality.md`` and
``examples/mesh_quality/run_benchmark.py`` for the numbers, including what could not be
run in this environment and why.
"""
from __future__ import annotations

from .dialects import ASPECT_DIALECTS, aspect_dialects, conservative_limit, dialect_audit
from .metrics import ELEMENT_KINDS, METRIC_NAMES, ElementQuality, hex_metrics, quad_metrics, tet_metrics, tri_metrics
from .topology import (
    FaceNonOrthogonality,
    SurfaceTopology,
    gate_non_orthogonality,
    gate_surface,
    non_orthogonality,
    surface_topology,
)
from .report import (
    DEFAULT_GATE,
    DEFAULT_LIMITS,
    Check,
    GateReport,
    MeshQuality,
    MetricLimit,
    QualityGate,
    SummaryStats,
    measure,
    summary_stats,
)

__all__ = [
    "DEFAULT_GATE",
    "DEFAULT_LIMITS",
    "ELEMENT_KINDS",
    "FaceNonOrthogonality",
    "METRIC_NAMES",
    "Check",
    "ASPECT_DIALECTS",
    "aspect_dialects",
    "conservative_limit",
    "dialect_audit",
    "ElementQuality",
    "GateReport",
    "MeshQuality",
    "MetricLimit",
    "QualityGate",
    "SummaryStats",
    "SurfaceTopology",
    "gate_non_orthogonality",
    "gate_surface",
    "hex_metrics",
    "measure",
    "non_orthogonality",
    "quad_metrics",
    "summary_stats",
    "surface_topology",
    "tet_metrics",
    "tri_metrics",
]