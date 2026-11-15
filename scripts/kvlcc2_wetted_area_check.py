"""Second-hull cross-check of the imported-hull wetted-area stage (GAP-CAD-T1-IGES-PIPE).

Round 4 of this anchor (``docs/gap_cad_t1_iges_pipe_round4.md`` §7) shipped
``tensorcad.geomquality.hull_surface.wetted_area`` — an occlusion-based de-duplication of
the overlapping skin that publisher IGES patch-sets produce — and validated it on **one**
hull only (KCS: +2.07 % against the workshop value 9 424 m²).  It recorded the second
workshop hull as *not run, not passed*.

This script closes exactly that debt and nothing else.  It takes the publisher KVLCC2
archive (``geometry/hulls/KVLCC2/kvlcc2.igs``, 19 trimmed faces), reuses the measured
intake of :func:`tensorcad.fleet_baseline.build_ship` (the unit gate and the
density-follows-adopted-factor fix are what make this file importable at all), and then
measures the hull **with a different estimator** than ``build_ship`` uses:

* ``build_ship`` integrates ∇ from the trim-loop sampling it just performed,
* :func:`tensorcad.geomquality.hull_surface.hull_characteristics` re-derives ∇, Cb *and*
  the wetted area from the triangle mesh by envelope-breadth rays along ``±y`` and
  normal-ray occlusion.

Two estimators on the same input, judged against the published particulars
(SIMMAN 2020: ∇ = 312 622 m³, S = 27 194 m², C_B = 0.8098, L_pp = 320.0 m, B = 58.0 m,
T = 20.8 m).  Every number written to ``runs/kvlcc2_wetted_area/evidence.json`` is read
back from the run; nothing is typed in by hand.

Reproduce::

    python scripts/kvlcc2_wetted_area_check.py            # writes runs/kvlcc2_wetted_area/
"""
from __future__ import annotations

import json
import sys
import time
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "runs" / "kvlcc2_wetted_area"
OUT.mkdir(parents=True, exist_ok=True)

from tensorcad.fleet_baseline import build_ship, specs_from_json  # noqa: E402
from tensorcad.geomquality.hull_surface import (  # noqa: E402
    classify_model,
    hull_characteristics,
    validate_against_particulars,
    wetted_area,
)

#: Published particulars, SIMMAN 2020 KVLCC2 (the same record the fleet spec carries,
#: see benchmarks/verified/baseline_fleet/fleet_specs.json -> "published_volume_m3").
PUBLISHED = {
    "volume": 312622.0,
    "s_wetted": 27194.0,
    "cb": 0.8098,
    "lpp": 320.0,
    "beam": 58.0,
    "tolerance_percent": 3.0,
    "source": "SIMMAN 2020 KVLCC2 principal particulars (Table 1)",
}


def main() -> int:
    spec = {s.name: s for s in specs_from_json(
        ROOT / "benchmarks/verified/baseline_fleet/fleet_specs.json")}["KVLCC2"]

    evidence: dict = {
        "anchor": "GAP-CAD-T1-IGES-PIPE",
        "claim_under_test": (
            "the shipped wetted_area dedup, validated on KCS only, also holds on the "
            "second publisher workshop hull"),
        "input": {"path": str(spec.iges_path), "sha256": spec.iges_sha256,
                  "bytes": Path(spec.iges_path).stat().st_size},
        "published": PUBLISHED,
        "repro": "python scripts/kvlcc2_wetted_area_check.py",
    }

    t0 = time.time()
    mesh, report = build_ship(spec, dx=1.0, dz=0.1, density_follows_adopted_factor=True)
    t_intake = time.time() - t0
    vertices = list(mesh.vertices)
    faces = [tuple(face) for face in mesh.faces]
    evidence["intake"] = {
        "seconds": round(t_intake, 3),
        "vertices": len(vertices),
        "triangles": len(faces),
        "dx_m": 1.0,
        "dz_m": 0.1,
        "density_follows_adopted_factor": True,
        "unit_gate": {k: v for k, v in report["source"]["unit_gate"].items()
                      if k != "candidates"},
        "trim_sampling": {
            "patches_used": report["source"]["trim_sampling"].get("patches_used"),
            "patches_total": report["source"]["trim_sampling"].get("patches_total"),
            "samples": report["source"]["trim_sampling"].get("samples"),
        },
        # estimator #1: what build_ship itself derived from the trim-loop sampling
        "build_ship_derived": {
            "volume_m3": report["derived"].get("volume_m3"),
            "cb": report["derived"].get("cb"),
            "volume_rel_error_percent":
                report["derived"].get("published_volume_rel_error_percent"),
            "watertight": report["checks"].get("watertight"),
            "edge_manifold": report["checks"].get("edge_manifold"),
            "verdict": report.get("verdict"),
        },
    }

    waterline = float(spec.draft)
    evidence["mesh_bounds"] = {
        "units": mesh.units,
        "x": [min(p[0] for p in vertices), max(p[0] for p in vertices)],
        "y": [min(p[1] for p in vertices), max(p[1] for p in vertices)],
        "z": [min(p[2] for p in vertices), max(p[2] for p in vertices)],
        "draft_m": waterline,
        "waterline_inside_mesh_z_range": bool(
            min(p[2] for p in vertices) <= waterline <= max(p[2] for p in vertices)),
    }
    evidence["model"] = classify_model(vertices)

    # ---- estimator #2: independent, mesh-only -------------------------------
    t1 = time.time()
    characteristics = hull_characteristics(vertices, faces, waterline=waterline,
                                           stations=61, levels=48, dedup=True)
    t_chars = time.time() - t1
    validation = validate_against_particulars(characteristics, PUBLISHED,
                                             tolerance_percent=PUBLISHED["tolerance_percent"])

    t2 = time.time()
    naive_area, _ = wetted_area(vertices, faces, waterline=waterline, dedup=False)
    t_naive = time.time() - t2

    evidence["hull_characteristics"] = {
        "seconds": round(t_chars, 3),
        "waterline_m": waterline,
        "volume_m3": characteristics.volume,
        "s_wetted_m2": characteristics.s_wetted,
        "s_wetted_naive_m2": characteristics.s_wetted_naive,
        "s_wetted_discarded_m2": characteristics.s_wetted_discarded,
        "cb": characteristics.cb,
        "cm": characteristics.cm,
        "lwl_m": characteristics.lwl,
        "beam_m": characteristics.beam,
        "symmetry_factor": characteristics.symmetry_factor,
        "stations": characteristics.station_count,
        "levels": characteristics.level_count,
    }
    evidence["wetted_area_naive"] = {
        "seconds": round(t_naive, 3),
        "area_m2": naive_area,
        "rel_error_percent": 100.0 * (naive_area - PUBLISHED["s_wetted"]) / PUBLISHED["s_wetted"],
    }
    evidence["validation"] = validation
    evidence["estimator_agreement"] = {
        "volume_build_ship_vs_hull_characteristics_percent": (
            100.0 * (characteristics.volume - report["derived"]["volume_m3"])
            / report["derived"]["volume_m3"]
            if report["derived"].get("volume_m3") else None),
    }

    # ---- verdict of the debt this script exists to settle -------------------
    rows = validation["quantities"]
    area_row = rows.get("s_wetted", {})
    evidence["verdict"] = {
        "s_wetted_within_3pct": bool(area_row.get("passed")),
        "s_wetted_error_percent": area_row.get("error_percent"),
        "volume_within_3pct": bool(rows.get("volume", {}).get("passed")),
        "cb_within_3pct": bool(rows.get("cb", {}).get("passed")),
        "all_quantities_within_3pct": bool(validation.get("passed")),
        "second_hull_status": "PASSED" if area_row.get("passed") else "FAILED",
        "naive_sum_still_shows_the_defect": bool(
            evidence["wetted_area_naive"]["rel_error_percent"] > 15.0),
    }
    evidence["total_seconds"] = round(time.time() - t0, 3)

    (OUT / "evidence.json").write_text(json.dumps(evidence, indent=2, default=str))
    print(json.dumps(evidence["verdict"], indent=2))
    print(json.dumps({k: evidence["hull_characteristics"][k] for k in
                      ("volume_m3", "s_wetted_m2", "s_wetted_naive_m2",
                       "s_wetted_discarded_m2", "cb", "lwl_m", "beam_m")}, indent=2))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception:  # pragma: no cover - the failure itself is evidence
        OUT.mkdir(parents=True, exist_ok=True)
        (OUT / "failure.txt").write_text(traceback.format_exc())
        raise