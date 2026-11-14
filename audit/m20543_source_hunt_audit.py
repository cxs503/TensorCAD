#!/usr/bin/env python3
"""Board audit of m20473 (GAP-CAD-T1-SOURCE-HUNT) — independent, self-verifying re-run.

This is a *delivery audit*, not a re-do of the mission.  Every number in the verdict and in
the accompanying report is produced here, from the artefacts the mission left behind, and
written to ``artifacts/m20543/evidence_recheck.json``.  Nothing is hand-copied.

Checks
------
1. existence of the claimed evidence files (impl / entry / run artefacts / tests / integration)
2. the mission's own repro command, re-executed into a scratch directory
3. the mission's test suite, re-executed (the mission itself reported 16 passed / 1 failed)
4. the benchmark row, re-derived from the *record on disk* (ref vs got, tolerance)
5. independence of the benchmark reference (is the SIMMAN digit actually on the SIMMAN host?)
6. staleness of the graded record against its own raw cache (timestamps + digest re-computation)
7. which of the six standard hulls actually have geometry provenance
8. whether the delivery reached the repository / the platform benchmark tables

Usage
-----
    python audit/m20543_source_hunt_audit.py \
        --tree /root/tp-agent/work/m20473 \
        --out artifacts/m20543
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
ANCHOR = "GAP-CAD-T1-SOURCE-HUNT"
MISSION = 20473
TOLERANCE = 0.03

# every artefact the mission claims to have produced, keyed by the evidence slot it fills
CLAIMED = {
    "impl_files": [
        "tensorcad/fleet/provenance.py",
        "tools/baseline_fleet/hunt_sources.py",
        "benchmarks/verified/baseline_fleet/source_hunt/run_benchmark.py",
    ],
    "entry": ["tools/baseline_fleet/hunt_sources.py"],
    "run_artifacts": [
        "artifacts/baseline_fleet/source_hunt.json",
        "artifacts/baseline_fleet/source_hunt_benchmark.json",
        "artifacts/baseline_fleet/source_hunt_benchmark.md",
    ],
    "tests": ["tests/test_source_hunt.py"],
    "integration_callers": [
        "tensorcad/fleet/official.py",
        "tools/baseline_fleet/official_fleet.py",
    ],
}

REPRO_BENCH = "benchmarks/verified/baseline_fleet/source_hunt/run_benchmark.py"
REPRO_HUNT = "tools/baseline_fleet/hunt_sources.py --out artifacts/baseline_fleet --cache artifacts/baseline_fleet/raw --fetch-geometry"

# the six standard hulls the anchor names; what "six-model" means for grading
HULLS = ["KCS", "KVLCC2", "JBC", "ONRT", "DTC", "SERIES60"]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run(cmd: list[str], cwd: Path, timeout: int = 600) -> dict[str, Any]:
    proc = subprocess.run(cmd, cwd=str(cwd), capture_output=True, text=True, timeout=timeout)
    return {
        "cmd": " ".join(cmd),
        "cwd": str(cwd),
        "returncode": proc.returncode,
        "stdout_tail": proc.stdout.strip().splitlines()[-6:],
        "stderr_tail": proc.stderr.strip().splitlines()[-6:],
    }


def present(tree: Path, rel: str) -> dict[str, Any]:
    path = tree / rel
    return {
        "path": rel,
        "present": path.is_file(),
        "bytes": path.stat().st_size if path.is_file() else 0,
        "sha256": sha256_file(path) if path.is_file() and path.stat().st_size else "",
    }


def check_evidence(tree: Path) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for slot, files in CLAIMED.items():
        out[slot] = [present(tree, rel) for rel in files]
        out[slot + "_all_present"] = all(item["present"] for item in out[slot])
    return out


def check_committed(repo: Path, files: list[str]) -> dict[str, Any]:
    """Are the deliverable files tracked by git (i.e. really delivered) or worktree-only?"""
    result: dict[str, Any] = {}
    for rel in files:
        proc = subprocess.run(
            ["git", "ls-files", "--error-unmatch", rel],
            cwd=str(repo),
            capture_output=True,
            text=True,
        )
        result[rel] = {"tracked": proc.returncode == 0}
    for rel in files:
        result[rel]["modified_uncommitted"] = False
        if result[rel]["tracked"]:
            diff = subprocess.run(
                ["git", "diff", "--quiet", "HEAD", "--", rel], cwd=str(repo), capture_output=True
            )
            result[rel]["modified_uncommitted"] = diff.returncode == 1
    head = subprocess.run(
        ["git", "log", "-1", "--format=%H %ci %s"], cwd=str(repo), capture_output=True, text=True
    )
    result["head"] = head.stdout.strip()
    # the delivery counts as landed only if every deliverable is tracked *and* carries no
    # uncommitted modification; pre-existing files that were merely edited do not count as landed
    result["delivered_to_head"] = all(
        entry["tracked"] and not entry["modified_uncommitted"]
        for key, entry in result.items()
        if isinstance(entry, dict)
    )
    result["new_files_untracked"] = [k for k, v in result.items() if isinstance(v, dict) and not v["tracked"]]
    result["files_modified_uncommitted"] = [
        k for k, v in result.items() if isinstance(v, dict) and v.get("modified_uncommitted")
    ]
    return result


def check_benchmark(tree: Path, scratch: Path) -> dict[str, Any]:
    out_dir = scratch / "bench"
    run_result = run([sys.executable, REPRO_BENCH, "--out", str(out_dir)], tree)
    payload_path = out_dir / "source_hunt_benchmark.json"
    payload = json.loads(payload_path.read_text()) if payload_path.is_file() else {}
    rows = []
    for row in payload.get("rows", []):
        ref = float(row["reference"])
        got = float(row["measured"])
        error = abs(got - ref) / abs(ref) if ref else float("inf")
        rows.append(
            {
                "kind": row["kind"],
                "ship": row["ship"],
                "quantity": row["quantity"],
                "reference": ref,
                "reference_basis": row.get("reference_basis", ""),
                "measured": got,
                "recomputed_relative_error": error,
                "claimed_relative_error": row.get("relative_error"),
                "within_3pct_recomputed": error <= TOLERANCE,
                "page_sha256": row.get("page_sha256", ""),
            }
        )
    return {
        "repro": REPRO_BENCH,
        "execution": run_result,
        "rows": rows,
        "n_rows": len(rows),
        "n_within_3pct": sum(1 for r in rows if r["within_3pct_recomputed"]),
        "all_within_3pct": bool(rows) and all(r["within_3pct_recomputed"] for r in rows),
        "unreferenced": payload.get("unreferenced", []),
        "kinds_present": sorted({r["kind"] for r in rows}),
    }


def check_reference_independence(tree: Path, record: dict, bench: dict) -> dict[str, Any]:
    """Is the benchmark's reference digit really published by a *different* host?

    The single scored row compares a model-scale denominator derived from the T2015 page
    (full/model waterline length column pair) against the figure printed on the SIMMAN page.
    Independence therefore means: the reference literal must be found in the SIMMAN bytes,
    and the measured value must be derivable from the T2015 bytes only.
    """
    raw = tree / "artifacts" / "baseline_fleet" / "raw"
    simman_files = sorted(raw.glob("simman_summary*.bin"))
    t2015_files = sorted(raw.glob("t2015_onrt_gc*.bin"))
    out: dict[str, Any] = {
        "simman_cache_files": [p.name for p in simman_files],
        "t2015_onrt_cache_files": [p.name for p in t2015_files],
    }
    row = bench["rows"][0] if bench["rows"] else None
    if row and simman_files:
        text = simman_files[0].read_text(errors="ignore")
        literal = f"{row['reference']:.4f}"
        out["reference_literal"] = literal
        out["reference_found_on_simman_host"] = literal in text
        # the measured denominator must be derivable from the T2015 page alone
        block = (record.get("scale_factors") or {}).get(row["ship"], {})
        full = block.get("full", {}).get("length_m")
        model = block.get("model", {}).get("length_m")
        if full and model:
            out["t2015_column_pair"] = {"full_m": full, "model_m": model}
            out["derived_from_t2015"] = full / model
            out["derived_matches_measured"] = abs(full / model - row["measured"]) < 1e-9
            out["t2015_page_digest_matches_record"] = (
                block.get("full", {}).get("page_sha256")
                == block.get("model", {}).get("page_sha256")
            )
        out["different_hosts"] = bool(
            "simman" in str(block.get("reference_source_key", ""))
            and "t2015" in str(block.get("model", {}).get("source_key", ""))
        )
    return out


def check_record_staleness(tree: Path, record: dict) -> dict[str, Any]:
    """The graded record vs the raw cache it is supposed to describe.

    The record was written at ``generated_utc``; the geometry archives carry their own
    ``retrieved_utc``.  Anything retrieved after the record was written cannot be reflected in
    its acceptance verdicts — that is staleness, and it is checkable without opinions.
    """
    raw = tree / "artifacts" / "baseline_fleet" / "raw"
    generated = record.get("generated_utc", "")
    metas = sorted(raw.glob("*.meta.json"))
    later: list[dict[str, Any]] = []
    for meta_path in metas:
        meta = json.loads(meta_path.read_text())
        if generated and meta.get("retrieved_utc", "") > generated:
            later.append(
                {
                    "meta": meta_path.name,
                    "url": meta.get("url"),
                    "http_status": meta.get("http_status"),
                    "retrieved_utc": meta.get("retrieved_utc"),
                }
            )
    # what does the *record* say about the geometry sources?
    said = {}
    for key, entry in (record.get("sources") or {}).items():
        if any(tag in key for tag in ("gap", "iges", "geometry", "dtc")):
            said[key] = {
                "status_class": entry.get("status_class"),
                "http_status": entry.get("http_status"),
                "sha256": entry.get("sha256", ""),
                "url": entry.get("url"),
            }
    return {
        "record_generated_utc": generated,
        "sources_retrieved_after_record": later,
        "record_says_about_geometry": said,
    }


def check_geometry_bytes(tree: Path, record: dict) -> dict[str, Any]:
    """Do we actually hold publisher geometry, and do the two hosts agree?

    The acceptance block claims ONRT/JBC geometry is WITHHELD ("archive unreachable").  The
    raw cache is the counter-evidence: compute the digests and compare them with the record.
    """
    raw = tree / "artifacts" / "baseline_fleet" / "raw"
    out: dict[str, Any] = {"archives": [], "iges_extracted": []}
    for pattern, label in (
        ("t2015_onrt_wo_gap*.bin", "ONRT (T2015 host)"),
        ("simman_onrt_wo_gap*.bin", "ONRT (SIMMAN host)"),
        ("t2015_jbc_iges*.bin", "JBC IGES zip (T2015 host)"),
    ):
        for path in sorted(raw.glob(pattern)):
            out["archives"].append(
                {"label": label, "file": path.name, "bytes": path.stat().st_size, "sha256": sha256_file(path)}
            )
    for path in sorted(raw.glob("*.igs")):
        text = path.read_text(errors="ignore")
        out["iges_extracted"].append(
            {
                "file": path.name,
                "bytes": path.stat().st_size,
                "lines": text.count("\n"),
                "sha256": sha256_file(path),
                "igs_footer": "S      1" in text.splitlines()[0] if text else False,
            }
        )
    by_label: dict[str, set[str]] = {}
    for archive in out["archives"]:
        by_label.setdefault(archive["label"].split(" (")[0], set()).add(archive["sha256"])
    out["onrt_cross_host_identical"] = (
        len(by_label.get("ONRT", set())) == 1 and len(by_label.get("ONRT", set())) == 1
        and sum(1 for a in out["archives"] if a["label"].startswith("ONRT")) == 2
    )
    out["acceptance_said"] = {
        key: (value.get("verdict") if isinstance(value, dict) else value)
        for key, value in (record.get("acceptance") or {}).items()
    }
    out["acceptance_evidence_said"] = {
        "dtc": (record.get("acceptance", {}).get("dtc_and_onrt_geometry", {}).get("evidence", {}).get("dtc", {}) or {}).get("verdict"),
        "jbc": (record.get("acceptance", {}).get("dtc_and_onrt_geometry", {}).get("evidence", {}).get("jbc", {}) or {}).get("verdict"),
        "onrt": (record.get("acceptance", {}).get("dtc_and_onrt_geometry", {}).get("evidence", {}).get("onrt", {}) or {}).get("verdict"),
    }
    return out


def check_hull_coverage(record: dict) -> dict[str, Any]:
    """Per hull: do we have particulars, a page digest, a scale factor, geometry?"""
    particulars = record.get("particulars") or {}
    scales = record.get("scale_factors") or {}
    geometry = record.get("geometry") or {}
    raw = json.dumps(record, ensure_ascii=False).lower()
    coverage = {}
    for hull in HULLS:
        key = hull.lower()
        coverage[hull] = {
            "particulars_grounded": hull in particulars,
            "page_sha256": (particulars.get(hull, {}) or {}).get("page_sha256", ""),
            "scale_factor": hull in scales,
            "scale_relative_error": (scales.get(hull, {}) or {}).get("relative_error"),
            "geometry_entry": any(key in str(k).lower() for k in geometry),
            "named_in_record": key in raw,
        }
    return {
        "hulls": coverage,
        "particulars_count": len(particulars),
        "scale_factor_count": len(scales),
        "geometry_entry_count": len(geometry),
    }


def check_integration(tree: Path, record: dict) -> dict[str, Any]:
    official = (tree / "tensorcad/fleet/official.py").read_text()
    fleet = (tree / "tools/baseline_fleet/official_fleet.py").read_text()
    out = {
        "def_line": next(
            (f"tensorcad/fleet/official.py:{i}" for i, line in enumerate(official.splitlines(), 1)
             if line.startswith("def verify_page_digests")),
            "",
        ),
        "call_line": next(
            (f"tools/baseline_fleet/official_fleet.py:{i}" for i, line in enumerate(fleet.splitlines(), 1)
             if "verify_page_digests(" in line and "import" not in line),
            "",
        ),
        "e2e_path_run_by_mission": False,  # mission's own result says extract was never run
    }
    # execute the seam for real
    proc = subprocess.run(
        [
            sys.executable,
            "-c",
            "import json;from tensorcad.fleet.official import verify_page_digests as v;"
            "r=json.load(open('artifacts/baseline_fleet/source_hunt.json'));print(json.dumps(v(r)))",
        ],
        cwd=str(tree),
        capture_output=True,
        text=True,
    )
    out["seam_returncode"] = proc.returncode
    if proc.returncode == 0:
        payload = json.loads(proc.stdout.strip().splitlines()[-1])
        out["seam_result"] = {
            "checked": payload.get("checked"),
            "problems": payload.get("problems"),
            "uncovered": payload.get("uncovered"),
        }
    else:
        out["seam_result"] = {"error": proc.stderr.strip().splitlines()[-3:]}
    return out


def check_platform_referee() -> dict[str, Any]:
    """The platform benchmark tables are the only referee: what do they say about this anchor?"""
    query = (
        "SELECT coalesce(json_agg(row_to_json(t)),'[]') FROM ("
        "SELECT case_name,metric,tolerance,status FROM tensor_benchmark_cases "
        "WHERE anchor_slug='GAP-CAD-T1-SOURCE-HUNT') t;"
    )
    cases = subprocess.run(
        ["psql", "-U", "postgres", "-d", "tensorplatform", "-tAc", query],
        capture_output=True,
        text=True,
    )
    results = subprocess.run(
        [
            "psql", "-U", "postgres", "-d", "tensorplatform", "-tAc",
            "SELECT coalesce(count(*),0) FROM tensorbench_results WHERE mission_id=20473;",
        ],
        capture_output=True,
        text=True,
    )
    reports = subprocess.run(
        [
            "psql", "-U", "postgres", "-d", "tensorplatform", "-tAc",
            "SELECT count(*) FROM tensor_reports WHERE mission_id=20473;",
        ],
        capture_output=True,
        text=True,
    )
    try:
        case_rows = json.loads(cases.stdout.strip() or "[]")
    except json.JSONDecodeError:
        case_rows = []
    return {
        "benchmark_cases_for_anchor": case_rows,
        "n_benchmark_cases": len(case_rows),
        "n_tensorbench_results_for_mission": int((results.stdout.strip() or "0") or 0),
        "n_platform_reports_for_mission": int((reports.stdout.strip() or "0") or 0),
        "psql_ok": cases.returncode == 0 and results.returncode == 0 and reports.returncode == 0,
    }


def run_hunt_verify(tree: Path, scratch: Path) -> dict[str, Any]:
    """Re-run the mission's offline digest re-check against the cache it shipped.

    ``--verify-only`` re-reads the record and recomputes every cached digest; it must not touch
    the source tree, so the record is copied into a scratch directory first.
    """
    hunt_dir = scratch / "hunt"
    hunt_dir.mkdir(parents=True, exist_ok=True)
    source = tree / "artifacts" / "baseline_fleet" / "source_hunt.json"
    (hunt_dir / "source_hunt.json").write_bytes(source.read_bytes())
    return run(
        [sys.executable, "tools/baseline_fleet/hunt_sources.py",
         "--out", str(hunt_dir), "--cache", str(tree / "artifacts/baseline_fleet/raw"),
         "--verify-only"],
        tree,
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tree", default="/root/tp-agent/work/m20473")
    parser.add_argument("--out", default="artifacts/m20543")
    args = parser.parse_args()
    tree = Path(args.tree).resolve()
    out_dir = (REPO_ROOT / args.out).resolve()
    scratch = out_dir / "scratch"
    scratch.mkdir(parents=True, exist_ok=True)

    record = json.loads((tree / "artifacts/baseline_fleet/source_hunt.json").read_text())
    evidence: dict[str, Any] = {
        "anchor": ANCHOR,
        "audited_mission": MISSION,
        "audit_mission": 20543,
        "tree": str(tree),
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "evidence_files": check_evidence(tree),
        "committed": check_committed(tree, [f for files in CLAIMED.values() for f in files]),
        "tests": run([sys.executable, "-m", "pytest", "tests/test_source_hunt.py", "-q"], tree),
        "hunt_verify_only": run_hunt_verify(tree, scratch),
        "benchmark": check_benchmark(tree, scratch),
        "record_verification_block": record.get("verification", {}),
        "repro_commands": {"benchmark": REPRO_BENCH, "hunt": REPRO_HUNT},
    }
    evidence["reference_independence"] = check_reference_independence(tree, record, evidence["benchmark"])
    evidence["record_staleness"] = check_record_staleness(tree, record)
    evidence["geometry_bytes"] = check_geometry_bytes(tree, record)
    evidence["hull_coverage"] = check_hull_coverage(record)
    evidence["integration"] = check_integration(tree, record)
    evidence["platform_referee"] = check_platform_referee()

    # ---- derived judgements (no hand-typed numbers) ---------------------------- #
    tests_out = "\n".join(evidence["tests"]["stdout_tail"])
    match = re.search(r"(\d+) failed, (\d+) passed", tests_out)
    evidence["tests_summary"] = (
        {"failed": int(match.group(1)), "passed": int(match.group(2))} if match
        else {"failed": None, "passed": None, "raw": tests_out}
    )
    evidence["verdict_basis"] = {
        "L2_impl_real": evidence["evidence_files"]["impl_files_all_present"]
        and evidence["evidence_files"]["entry_all_present"]
        and evidence["evidence_files"]["run_artifacts_all_present"],
        "L2_repro_executes": evidence["benchmark"]["execution"]["returncode"] == 0,
        "L2_benchmark_row_within_3pct": evidence["benchmark"]["all_within_3pct"],
        "L2_benchmark_row_is_not_a_computed_capability": len(evidence["benchmark"]["rows"]) == 1
        and evidence["benchmark"]["kinds_present"] == ["cross_host"],
        "L3_seam_exists": bool(evidence["integration"]["call_line"]),
        "L3_e2e_extract_run_by_mission": evidence["integration"]["e2e_path_run_by_mission"],
        "L4_platform_benchmark_registered": evidence["platform_referee"]["n_benchmark_cases"] > 0
        and evidence["platform_referee"]["n_tensorbench_results_for_mission"] > 0,
        "delivery_committed_to_git": evidence["committed"]["delivered_to_head"],
        "platform_report_published": evidence["platform_referee"]["n_platform_reports_for_mission"] > 0,
        "record_grade_stale_vs_cache": bool(
            evidence["record_staleness"]["sources_retrieved_after_record"]
        ),
        "red_test_present": evidence["tests_summary"]["failed"] not in (0, None),
    }
    (out_dir / "evidence_recheck.json").write_text(json.dumps(evidence, indent=1, ensure_ascii=False))
    print(json.dumps(evidence["verdict_basis"], indent=1))
    print("tests:", evidence["tests_summary"])
    print("benchmark rows:", evidence["benchmark"]["n_rows"], "within3%:", evidence["benchmark"]["n_within_3pct"])
    print("ref independence:", {k: v for k, v in evidence["reference_independence"].items() if not isinstance(v, dict)})
    print("stale sources:", [s["url"] for s in evidence["record_staleness"]["sources_retrieved_after_record"]])
    print("platform:", evidence["platform_referee"]["n_benchmark_cases"],
          evidence["platform_referee"]["n_tensorbench_results_for_mission"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())