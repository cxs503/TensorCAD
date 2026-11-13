"""Geometry-quality pipeline: CAD file in, gated hydrostatics out.

The chain this package implements is the one a shipyard actually needs and almost
never has end to end:

``IGES / STEP / native B-rep`` -> ``edge-face adjacency graph`` -> ``watertightness
0/1`` -> ``declared L/B/T and Cb`` -> ``seam continuity G0/G1/G2`` ->
``hydrostatics by patch quadrature`` -> ``grid convergence`` -> ``machine gate``.

Design rules, all of them deliberate:

* **Nothing is trusted because it renders.**  Topology is decided on the adjacency
  graph (:mod:`.brep`), hydrostatics on the surface integrals (:mod:`.hydro`), and the
  comparison is always against a reference that does not come from the model itself.
* **Two exchange routes, on purpose.**  STEP carries explicit topology; IGES does not,
  so the reader has to sew patch boundaries and a seam crack survives as a free edge
  (:mod:`.io`).  Reading both is the check that the audit sees what the file says.
* **Thresholds are calibrated, then scored on held-out cases** (:mod:`.continuity`),
  and a threshold that fails to separate is reported, not tuned away.
* **Every defect code carries a repair hint** (:mod:`.defects`), because a gate that
  only says "no" is not usable by the person who has to fix the model.

The acceptance numbers live in :data:`tensorcad.geomquality.pipeline.CRITERIA` and the
verdict is computed, never written by hand.
"""
from __future__ import annotations

from .brep import (
    REPAIR_HINTS,
    TOPOLOGY_CODES,
    BRep,
    Defect,
    Edge,
    EdgeRef,
    Face,
    audit_topology,
    merge_breps,
    mirror,
    sew,
)
from .continuity import (
    DEFAULT_ANGLE_LIMIT_DEG,
    DEFAULT_CURVATURE_LIMIT,
    SeamReport,
    audit_seams,
    calibrate,
    classify,
    continuity_evidence,
    curvature_distribution,
    label_seam_dataset,
    score,
    seam_metrics,
    seam_metrics_adaptive,
)
from .defects import DEFECT_CASES, DefectCase, detect_all, run_library
from .quality import (
    DEFAULT_CRITERIA,
    IDEAL_VALUES,
    METRIC_COMPARISON,
    METRIC_DEFINITIONS,
    MeshQualityReport,
    QualityCheck,
    QualitySummary,
    QualityVerdict,
    TriangleQuality,
    mesh_quality,
    quality_gate,
    triangle_quality,
)
from .mesh_defects import (
    MeshDefect,
    MeshDefectReport,
    RepairResult,
    name_mesh_defects,
    repair_facet_defects,
)
from .hydro import (
    DEFAULT_ORDER,
    PatchHydrostatics,
    convergence_study,
    facet_hydrostatics,
    patch_hydrostatics,
    richardson_gci,
)
from .io import (
    ExchangeError,
    brep_from_dict,
    brep_to_dict,
    load_cad,
    read_brep_json,
    read_iges,
    read_step,
    round_trip,
    save_cad,
    write_brep_json,
    write_iges,
    write_step,
)
from .pipeline import (
    CRITERIA,
    CriterionResult,
    GateError,
    criterion,
    enforce,
    report_criteria,
    run_evidence,
    verify_cad,
    worst_hydrostatic_error,
)

__all__ = [
    "CRITERIA",
    "DEFAULT_ANGLE_LIMIT_DEG",
    "DEFAULT_CURVATURE_LIMIT",
    "DEFAULT_ORDER",
    "DEFECT_CASES",
    "REPAIR_HINTS",
    "TOPOLOGY_CODES",
    "BRep",
    "MeshDefect",
    "MeshDefectReport",
    "RepairResult",
    "CriterionResult",
    "Defect",
    "DefectCase",
    "Edge",
    "EdgeRef",
    "ExchangeError",
    "Face",
    "GateError",
    "PatchHydrostatics",
    "SeamReport",
    "audit_seams",
    "audit_topology",
    "brep_from_dict",
    "brep_to_dict",
    "calibrate",
    "classify",
    "continuity_evidence",
    "convergence_study",
    "criterion",
    "curvature_distribution",
    "detect_all",
    "enforce",
    "facet_hydrostatics",
    "label_seam_dataset",
    "load_cad",
    "merge_breps",
    "name_mesh_defects",
    "DEFAULT_CRITERIA",
    "IDEAL_VALUES",
    "METRIC_COMPARISON",
    "METRIC_DEFINITIONS",
    "MeshQualityReport",
    "QualityCheck",
    "QualitySummary",
    "QualityVerdict",
    "TriangleQuality",
    "mesh_quality",
    "quality_gate",
    "triangle_quality",
    "mirror",
    "patch_hydrostatics",
    "read_brep_json",
    "read_iges",
    "read_step",
    "repair_facet_defects",
    "report_criteria",
    "richardson_gci",
    "round_trip",
    "run_evidence",
    "run_library",
    "save_cad",
    "score",
    "seam_metrics",
    "seam_metrics_adaptive",
    "sew",
    "verify_cad",
    "worst_hydrostatic_error",
    "write_brep_json",
    "write_iges",
    "write_step",
]
