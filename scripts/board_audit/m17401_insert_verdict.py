#!/usr/bin/env python3
"""Board-audit verdict for m17349 (GAP-CAD-B1-STEP-READER, TensorCAD).

Idempotent: one row per (repo, scope, verdict_by) — re-running updates in place,
so the maturity_gate daemon re-judges the same row instead of stacking verdicts.

Every evidence value is a real pointer that the gate can resolve:
  impl_files  -> absolute files in the audited worktree (m17349)
  run_artifacts/benchmark -> copies shipped in *this* worktree under
                             artifacts/board_audit/m17401/ (relative path resolves
                             against /root/tp-agent/work/m17401)

Level claimed: 3 (system integration).  L4 is refused, with reasons recorded in
`l4_refused_because` — see docs/board_audit/m17401/AUDIT.md.

Run:
  python scripts/board_audit/m17401_insert_verdict.py [--dry-run]
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import psycopg

DSN = ("host=127.0.0.1 port=5432 dbname=tensorplatform user=foundation "
       "password=uGTN-KwuQv_7l8Tw0MMWUx_kf-nF2Rzz")
REPO = "TensorCAD"
SCOPE = "m17349 GAP-CAD-B1-STEP-READER"
MISSION = 17401              # this audit mission
SUBJECT = 17349              # audited mission
BY = "board-audit"
LEVEL = 3

HERE = Path(__file__).resolve().parents[2]          # worktree root
AUD = HERE / "artifacts/board_audit/m17401"
SUB = Path("/root/tp-agent/work/m17349")


def rel(p: Path) -> str:
    return str(p.relative_to(HERE))


def build() -> tuple[int, dict, str]:
    v = json.loads((AUD / "verify.json").read_text())
    rerun = v["C1_rerun"]
    indep = v["C2_volume_independent"]
    real = v["C4_real_third_party_file"]
    bmk = v["C7_benchmark_table"]

    benchmark = {
        "criterion": "mesh volume vs closed-form nominal within 1 % (the mission's own acceptance line)",
        "reference_source": "closed forms recomputed by the auditor from the geometry table in "
                            "docs/gap_cad_b1_step_reader.md; NOT a CAD-kernel nominal",
        "cross_check": "volume re-integrated by the auditor's own divergence-theorem script on the "
                       "re-exported OBJ (indep_volume values below)",
        "cases": [{"case": n, "ref_m3": indep[n]["closed_form_m3"], "got_m3": indep[n]["volume_m3"],
                   "error_percent": round(indep[n]["error_percent"], 6),
                   "subject_reported_error_percent": round(rerun["per_file"][n + ".step"]["error_percent"], 6),
                   "agree": abs(indep[n]["error_percent"] - rerun["per_file"][n + ".step"]["error_percent"]) < 1e-6,
                   "passes_1pct": abs(indep[n]["error_percent"]) < 1.0,
                   "watertight": rerun["per_file"][n + ".step"]["watertight"]}
                  for n in indep],
        "convergence": {"angular_steps": [r["angular_step"] for r in rerun["convergence"]],
                        "ratio_per_halving_plate": [
                            round(abs(rerun["convergence"][i]["plate_with_hole.step"]
                                      / rerun["convergence"][i + 1]["plate_with_hole.step"]), 4)
                            for i in range(len(rerun["convergence"]) - 1)],
                        "observed_order": "ratio -> 4.00 at the finest pair (3.69/3.85/3.92/4.00), "
                                          "i.e. second order asymptotically"},
        "benchmark_table_rows_for_anchor": bmk["n_rows"],
        "benchmark_table_note": bmk["note"],
    }

    evidence = {
        "impl_files": [str(SUB / "tensorcad/step_brep.py"), str(SUB / "tensorcad/exchange.py")],
        "impl_detail": {
            "tensorcad/step_brep.py": "1837 lines, pure Python ISO-10303-21 parser + analytic "
                                      "PLANE/CYLINDRICAL/CONICAL/SPHERICAL/TOROIDAL tessellator + "
                                      "divergence-theorem volume/watertight/Euler checks (no eval, "
                                      "no random, no filename special-casing)",
            "tensorcad/exchange.py:622": "read_step() dispatches FACETED_BREP -> read_step_facets, "
                                         "everything else -> step_brep.read_step_brep().to_patch_geometry()",
            "impl_origin": "delivered by m17349 commits d6ebb27 + 849e2e7 on branch agent/m17349 "
                           "(not yet in /root/tensor-ecosystem/TensorCAD main; the audit ran the "
                           "branch worktree)",
        },
        "entry": str(SUB / "tensorcad/exchange.py"),
        "entry_detail": "public seam: exchange.read_step(:622) / load_exchange(:639) -> "
                        "step_brep.read_step_brep; CLI: python -m tensorcad.step_brep "
                        "{census|mesh} <file>; audit driver: tools/step_reader_audit.py",
        "run_artifacts": [rel(AUD / "verify.json"), rel(AUD / "rerun_step_reader_audit.log"),
                          rel(AUD / "pytest_auditor_rerun.log"),
                          rel(AUD / "refs/volume_validation.json"),
                          rel(AUD / "refs/convergence.json"),
                          rel(AUD / "refs/coverage.json")],
        "run_artifacts_detail": {
            "rerun": "auditor re-ran the subject's driver in the subject worktree: exit=0, "
                     "all_pass=true, 6.2 s, and the two JSON artifacts stayed byte-identical "
                     "(sha256 volume_validation=%s)" % rerun["sha256"]["volume_validation.json"],
            "independent": "auditor re-integrated the exported OBJ meshes with its own "
                           "divergence script (scripts/board_audit/m17401_verify_claims.py C2)",
        },
        "tests": ["tests/test_step_brep_reader.py::test_volume_matches_closed_form_nominal_within_one_percent",
                  "tests/test_step_brep_reader.py::test_mesh_is_watertight_and_euler_matches_topology",
                  "tests/test_step_brep_reader.py::test_refinement_converges_second_order",
                  "tests/test_step_brep_reader.py::test_face_with_hole_is_tessellated_not_refused",
                  "tests/test_step_brep_reader.py::test_unsupported_surface_is_refused_by_name_not_guessed",
                  "tests/test_step_brep_reader.py::test_real_nx_file_parses_and_its_inch_unit_is_honoured",
                  "tests/test_step_brep_reader.py::test_real_wireframe_file_is_refused_as_a_solid",
                  "tests/test_step_brep_reader.py::test_faceted_file_still_routes_to_the_faceted_reader",
                  "tests/test_exchange_verification.py (routing table asserts the REAL level)"],
        "tests_detail": "auditor re-ran: 17 passed (test_step_brep_reader.py) + 26 passed "
                        "(test_exchange_verification.py) = 43 passed, matching the subject's claim "
                        "(log: artifacts/board_audit/m17401/pytest_auditor_rerun.log)",
        "integration_callers": ["tensorcad/exchange.py"],
        "integration_detail": {
            "call_chain": "exchange.read_step -> step_brep.read_step_brep (exchange.py:634) -> "
                          "StepBrepImport.to_patch_geometry -> PatchGeometry.mesh",
            "consumers_in_repo": "load_exchange(exchange.py:639) -> examples/geometry_audit/run_audit.py:469 "
                                 "(reimported mesh feeds audit_watertightness); "
                                 "tests/test_step_brep_reader.py:17 imports exchange.read_step",
            "not_integrated": "no solver/mesher consumer: nothing in tensorcad/{solvers,mesh,cfd}* calls "
                              "exchange.read_step, so the STEP->solve path is still open. This is why "
                              "the verdict is L3 (library seam) and not L4.",
        },
        "benchmark": benchmark,
        "repro_cmd": "cd /root/tp-agent/work/m17349 && python tools/step_reader_audit.py "
                     "# exit 0, ~6 s, writes artifacts/gap-cad-b1-step-reader/{volume_validation,"
                     "convergence}.json + 3 figures; full claim re-verification: "
                     "python scripts/board_audit/m17401_verify_claims.py --out artifacts/board_audit/"
                     "m17401/verify.json",
        "external_anchor_state": "WITHHELD: the only real third-party file (NX 7.5 / ST-Developer AP203, "
                                 "2283 instances, inch LENGTH_UNIT 39.3700787) is a wireframe "
                                 "(GEOMETRICALLY_BOUNDED_WIREFRAME_SHAPE_REPRESENTATION), 0 B-rep solids; "
                                 "the reader refuses it with 'no MANIFOLD_SOLID_BREP'. No third-party CAD "
                                 "nominal exists on linux/loong64 (no OCP wheel) -> no external anchor, "
                                 "so L5 is refused and L4 is refused (see l4_refused_because).",
        "l4_refused_because": [
            "acceptance asked for 3 REAL STEP files -> mesh with volume vs CAD nominal <1 %; delivered "
            "real-file meshes = 0/3 (the 4 meshed solids are fixtures authored by tools/step_samples.py, "
            "provenance.json says so verbatim: 'not downloaded CAD parts')",
            "the reference volumes are closed forms written in the same repo (tools/step_samples.py), "
            "not a CAD kernel's nominal, so the reference is analytic-exact but self-authored",
            "tensor_benchmark_cases has 0 rows for anchor GAP-CAD-B1-STEP-READER (only TensorCAD anchor "
            "registered is GAP-CAD-T1-BASELINE-FLEET) -> no machine-adjudicated benchmark case exists",
        ],
        "claims_verified": {
            "43 tests passed": True,
            "real file 2283 instances / 11 B_SPLINE_CURVE_WITH_KNOTS / inch scale 39.37007874015748": True,
            "error magnitudes +0.0129 / -0.2341 / -0.2341 / 0.0000 %": True,
            "watertight with 0 boundary edges, Euler 0/2/0/2": True,
            "fixtures regenerate byte-identically (4/4 sha256 equal)": all(
                x["identical"] for x in v["C3_fixture_determinism"].values()),
            "unsupported surface refused by name (mutated file probe)": v["C5_unsupported_refusal"]["refusal"]
            is not None,
        },
        "findings": [
            {"severity": "blocking-for-L4", "id": "F1",
             "what": "external/real-CAD benchmark missing: 0/3 real STEP solids meshed, nominals self-authored",
             "evidence": "C4 refusal message; provenance.json note; tensor_benchmark_cases 0 rows for the anchor"},
            {"severity": "minor", "id": "F2",
             "what": "StepBrepImport.complete=True while watertight=False and volume_m3=None on an "
                     "internally inconsistent file (8 boundary edges): a consumer keying on `complete` "
                     "alone could mistake an unclosed mesh for a finished import",
             "evidence": "C6 probe /tmp/m17401_probe/inconsistent_scale.step -> "
                         "{complete:true, watertight:false, boundary_edges:8, volume_m3:null}",
             "suggested_fix": "document that `complete` means 'every face tessellated', or make it "
                              "`complete = faces_ok and watertight`"},
            {"severity": "info", "id": "F3",
             "what": "headline convergence ratio 3.999 is the finest halving; the sequence is "
                     "3.69/3.85/3.92/4.00, so 'second order' is supported asymptotically but the "
                     "single quoted number is the best pair",
             "evidence": "convergence.json records, ratio table in evidence.benchmark.convergence"},
            {"severity": "info", "id": "F4",
             "what": "anchor has no row in tensor_benchmark_cases, so future audits of this anchor have "
                     "no DB judge to read",
             "evidence": "verify.json C7 n_rows=0; TensorCAD anchors present: GAP-CAD-T1-BASELINE-FLEET"},
        ],
        "auditor_checks_not_run": [
            "cross-check against a third-party CAD kernel (OCP/FreeCAD): absent on this machine (no wheel)",
            "B_SPLINE_SURFACE / swept / trimmed-composite surfaces: refused by design, not exercised",
            "hull-scale solids (>1e3 m parts, large coordinates): only the +1 km origin-shift probe was run",
            "full test suite of the repo (only the two relevant suites were re-run)",
        ],
    }
    reason = (
        "L3: real AP203/AP214 analytic B-rep reader (1837 lines, no third-party parser) wired into the "
        "exchange seam and validated on 4 fixtures against analyst-recomputed closed forms "
        "(|err| <= 0.234 % < 1 %) with independent OBJ re-integration, deterministic re-run and 43 "
        "passing tests. L4 refused: the acceptance criterion (3 REAL STEP files -> mesh vs CAD nominal) "
        "is met 0/3 - the only real third-party file is a wireframe and is honestly refused, the "
        "reference volumes are self-authored closed forms, and tensor_benchmark_cases has no row for "
        "this anchor. Subject's own numbers all reproduced exactly; the gap is coverage, not honesty."
    )
    return LEVEL, evidence, reason


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    level, evidence, reason = build()
    ev = {**evidence, "reason": reason}
    blob = json.dumps(ev, ensure_ascii=False, indent=1)
    if args.dry_run:
        print(blob)
        return 0
    with psycopg.connect(DSN) as conn:
        conn.execute("DELETE FROM tensor_maturity_verdicts WHERE repo=%s AND scope=%s AND verdict_by=%s",
                     (REPO, SCOPE, BY))
        row = conn.execute(
            "INSERT INTO tensor_maturity_verdicts (repo, scope, mission_id, level, evidence, verdict_by) "
            "VALUES (%s,%s,%s,%s,%s::jsonb,%s) RETURNING id, gate_state",
            (REPO, SCOPE, MISSION, level, blob, BY)).fetchone()
        conn.commit()
    print(f"inserted verdict id={row[0]} level={level} gate_state={row[1]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())