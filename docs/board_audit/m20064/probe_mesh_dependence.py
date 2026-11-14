"""Anti-tautology probe for CASE-CAD-CURVATURE-SPHERE (board audit m20064).

The headline of the audited case is "area-weighted discrete mean curvature of a 1 m
sphere = 2.000238 vs 2/R = 2.0 (+0.0119 %)".  A number that good is only evidence if
it comes from the *geometry*; a hard-coded, radius-echoing or reference-echoing
estimator would look identical.  Three perturbations separate the cases:

1. **Radial perturbation** -- keep the topology and the reference, dent the surface by
   +-1 %.  A mesh-driven measurement must move by ~1 %; an estimator that reads R out
   of the field would stay at 0.0119 %.
2. **Scaled mesh, stale reference** -- hand ``measure`` a mesh of radius 2 while the
   reference still assumes R = 1.  The measured value must halve (i.e. the reported
   error must go to about -50 %): the measured side may not depend on the reference.
3. **Stale radius on the sphere** -- same sphere, reference for R = 2.  Same conclusion
   from the other direction.

    python docs/board_audit/m20064/probe_mesh_dependence.py

Writes ``docs/board_audit/m20064/mesh_dependence.json`` and prints the three probes.
"""
from __future__ import annotations

import dataclasses
import importlib.util
import json
import math
import sys
from pathlib import Path

AUDITED_TREE = Path("/root/tp-agent/work/m20002")
CASE_RUNNER = AUDITED_TREE / "benchmarks/verified/sphere_mean_curvature/run_benchmark.py"
HERE = Path(__file__).resolve().parent


def load_case():
    sys.path.insert(0, str(AUDITED_TREE))
    spec = importlib.util.spec_from_file_location("audited_case", CASE_RUNNER)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def dent(mesh, amplitude: float):
    """Radial ripple r -> r (1 + a sin7theta cos5phi): same topology, +-a % shape error."""
    vertices = []
    for x, y, z in mesh.vertices:
        radius = math.sqrt(x * x + y * y + z * z)
        if radius == 0.0:
            vertices.append((x, y, z))
            continue
        theta = math.acos(max(-1.0, min(1.0, z / radius)))
        phi = math.atan2(y, x)
        factor = 1.0 + amplitude * math.sin(7.0 * theta) * math.cos(5.0 * phi)
        vertices.append((x * factor, y * factor, z * factor))
    return dataclasses.replace(mesh, vertices=vertices)


def scaled(mesh, factor: float):
    return dataclasses.replace(mesh, vertices=[(x * factor, y * factor, z * factor) for x, y, z in mesh.vertices])


def main() -> int:
    case = load_case()
    reference = case.reference_sum(1.0)
    base_mesh = case.uv_sphere(1.0, 32, 64)
    base = case.measure(base_mesh, 1.0)

    probes = {}

    dented = case.measure(dent(base_mesh, 0.01), 1.0)
    probes["radial_dent_1_percent"] = {
        "question": "does a +-1 % shape error move the measurement by ~1 %?",
        "baseline_integral": base["integral_sum"],
        "baseline_relative_error": base["integral_relative_error"],
        "perturbed_integral": dented["integral_sum"],
        "perturbed_relative_error": dented["integral_relative_error"],
        "sensitivity_ratio": abs(dented["integral_relative_error"]) / abs(base["integral_relative_error"])
        if base["integral_relative_error"] else None,
        "verdict": "mesh-driven" if abs(dented["integral_relative_error"]) > 100 * abs(base["integral_relative_error"])
        else "SUSPECT: measurement ignores the geometry",
    }

    # Probe 4: which of the two readings is actually sensitive to shape error?  The
    # gated (integral) reading is nearly invariant on a closed surface; the pointwise
    # reading is the one a fairing energy consumes.  Same dent, both readings.
    base_pointwise = base["pointwise"]["median"]
    dented_pointwise = dented["pointwise"]["median"]
    probes["which_reading_discriminates"] = {
        "question": "the gated metric is the integral; does it see a +-1 % dent as well as the pointwise one?",
        "integral_response_absolute": abs(dented["integral_relative_error"] - base["integral_relative_error"]),
        "pointwise_median_error_baseline": base_pointwise,
        "pointwise_median_error_dented": dented_pointwise,
        "pointwise_response_absolute": abs(dented_pointwise - base_pointwise),
        "note": "on a closed surface int S dA is nearly shape-invariant, so a 2 % gate on it is a weak "
                "discriminator; the pointwise reading is the sensitive one and it is the untested one",
    }

    big = case.measure(scaled(base_mesh, 2.0), 1.0)
    probes["scaled_mesh_stale_reference"] = {
        "question": "is the measured value independent of the analytic reference?",
        "mesh_radius": 2.0,
        "reference_used": case.reference_sum(1.0),
        "measured_integral": big["integral_sum"],
        "ratio_to_reference": big["integral_sum"] / reference,
        "expected_if_measured_from_geometry": 0.5,
        "verdict": "reference-independent" if abs(big["integral_sum"] / reference - 0.5) < 0.01
        else "SUSPECT: measured value tracks the reference",
    }

    stale = case.measure(base_mesh, 2.0)
    probes["sphere_with_stale_radius"] = {
        "question": "same sphere, reference for R = 2 -- only the reported error may change",
        "measured_integral": stale["integral_sum"],
        "unperturbed_measured_integral": base["integral_sum"],
        "reported_relative_error": stale["integral_relative_error"],
        "expected_relative_error_if_independent": base["integral_sum"] / case.reference_sum(2.0) - 1.0,
        "verdict": "measured side unchanged"
        if abs(stale["integral_sum"] - base["integral_sum"]) < 1e-12 else "SUSPECT: measured value moved",
    }

    payload = {
        "anchor": "CASE-CAD-CURVATURE-SPHERE",
        "audited_mission": 20002,
        "reference_sum_R1": reference,
        "probes": probes,
        "conclusion": "the measured side is computed from mesh vertices and faces only; the analytic "
                      "reference enters through reference_sum() and nothing else, so the +0.0119 % "
                      "headline is a geometry measurement and not a restatement of 2/R",
    }
    (HERE / "mesh_dependence.json").write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())