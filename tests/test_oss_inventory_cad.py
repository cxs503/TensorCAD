"""Invariants for the FEAT-OSS-INVENTORY-CAD round-1 census.

These tests are the machine gate for the survey deliverables.  They assert *properties
the census must satisfy*, not that some script ran:

* the census is large enough and every entry is pinned to a commit;
* tiers are from the declared set and the tier legend is respected;
* every licence verdict is accompanied by the licence file it was read from;
* every lineage pointer is either ``path:line-range`` (citable) or explicitly
  flagged as unverified - a bare file name may never masquerade as a pointer;
* every gap anchor has an acceptance criterion, a feature point and a priority;
* the rendered report cannot drift from ``stats.json`` (numbers are interpolated);
* the artefacts are hash-pinned by ``manifest.json``.

Run: ``pytest -q tests/test_oss_inventory_cad.py`` (offline; no network).
"""
from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ART = ROOT / "artifacts" / "oss-inventory-cad-v1"
POINTER_RE = re.compile(r"^[\w./\-]+\.[A-Za-z0-9_]+:\d+(?:-\d+)?$")
LICENSE_CLASSES = {"permissive", "copyleft-weak", "copyleft-strong", "restricted",
                   "no-license", "unknown"}
TRANSFERS = {"fits", "adapt", "reference"}


def _json(name: str) -> dict:
    return json.loads((ART / name).read_text())


def test_census_size_and_pinning():
    reg = _json("registry.json")["repositories"]
    surveyed = [r for r in reg if r["status"] == "surveyed"]
    assert len(surveyed) >= 20, f"only {len(surveyed)} repositories surveyed"
    for r in surveyed:
        assert re.fullmatch(r"[0-9a-f]{40}", r["commit"] or ""), f"{r['repo']}: commit not pinned"
        assert r["url"] == f"https://github.com/{r['repo']}"
        assert r["n_source_paths"] and r["n_source_paths"] > 0


def test_failed_repos_are_recorded_not_dropped():
    reg = _json("registry.json")["repositories"]
    stats = _json("stats.json")
    failed = [r for r in reg if r["status"] != "surveyed"]
    assert len(failed) == stats["n_failed"]
    assert len(reg) == stats["n_surveyed"] + stats["n_failed"]
    # the deliberate negative control must appear as a failure, with a reason
    neg = [r for r in reg if r["repo"].startswith("tensor-ecosystem-nonexistent/")]
    assert neg and neg[0]["status"] != "surveyed" and neg[0]["clone_error"]


def test_tiers_are_declared():
    pool = {e["repo"]: e for e in
            json.loads((ROOT / "tools/oss_inventory_cad/repos.json").read_text())["repos"]}
    for r in _json("registry.json")["repositories"]:
        assert r["tier"] in (1, 2, 3), r["repo"]
        assert r["tier"] == pool[r["repo"]]["tier"], f"{r['repo']}: tier drifted from repos.json"


def test_license_verdicts_are_backed_by_a_file():
    for r in _json("registry.json")["repositories"]:
        if r["status"] != "surveyed":
            continue          # clone failures carry no licence verdict at all
        assert r["license_class"] in LICENSE_CLASSES, r["repo"]
        if r["license_class"] in ("permissive", "copyleft-weak", "copyleft-strong",
                                  "restricted"):
            assert r["license_tokens"], f"{r['repo']}: class without licence tokens"
            assert r["license_file"], f"{r['repo']}: class without licence file"
        if r["license_class"] == "no-license":
            assert not r["license_tokens"] and r["license_file"] is None
        # reuse policy must be a deterministic function of the licence class
        expect = {"permissive": "fits", "copyleft-weak": "adapt"}.get(r["license_class"],
                                                                     "reference")
        assert r["reuse_policy"] == expect, f"{r['repo']}: {r['license_class']}->{r['reuse_policy']}"


def test_lineage_pointers_are_citable_or_flagged():
    rows = _json("lineage.json")["rows"]
    assert len(rows) >= 100, f"only {len(rows)} lineage rows"
    verified = 0
    for row in rows:
        assert row["feature"].startswith("F")
        assert row["transfer"] in TRANSFERS
        if row["pointer_status"] == "blob_verified":
            assert POINTER_RE.fullmatch(row["pointer"] or ""), row["pointer"]
            verified += 1
        else:
            assert row["pointer"] is None, "unverified rows must not carry a pointer"
            assert row["evidence_kind"] == "path_signature_only"
    stats = _json("stats.json")
    assert verified == stats["n_pointer_verified"]
    assert stats["pointer_verified_pct"] >= 60.0


def test_every_feature_point_is_swept_and_attributed():
    fps = json.loads((ROOT / "tools/oss_inventory_cad/feature_points.json").read_text())
    stats = _json("stats.json")
    assert len(fps["feature_points"]) == stats["n_features"] == 16
    for fid, cov in stats["feature_coverage"].items():
        assert cov["n_repos"] > 0, f"{fid}: no repository implements this function point"
        assert cov["tensorcad_status"] in ("real", "partial", "absent")
        assert cov["tensorcad_evidence"], f"{fid}: status without code evidence"
        assert cov["n_pointer_verified"] <= cov["n_repos"]


def test_gap_anchors_are_actionable():
    anchors = _json("gap_anchors.json")["anchors"]
    assert len(anchors) >= 6, f"only {len(anchors)} GAP-CAD anchors"
    ids = [a["anchor"] for a in anchors]
    assert len(set(ids)) == len(ids)
    for a in anchors:
        assert re.fullmatch(r"GAP-CAD-\d{3}", a["anchor"])
        assert a["current_status"] in ("partial", "absent")
        assert len(a["acceptance"]) > 40, f"{a['anchor']}: acceptance criterion too thin"
        assert a["priority"] in ("P1", "P2")
        for s in a["sources"]:
            assert s["transfer"] in TRANSFERS
            assert POINTER_RE.fullmatch(s["pointer"] or "") , s["pointer"]
    # a real (REAL) feature point must never be proposed as a gap
    stats = _json("stats.json")
    for a in anchors:
        assert stats["feature_coverage"][a["feature"]]["tensorcad_status"] != "real"


def test_report_is_generated_from_stats_not_hand_copied():
    stats = _json("stats.json")
    report = (ROOT / "docs/oss_inventory_cad/REPORT-2026-11-11.md").read_text()
    for token in (f"{stats['n_surveyed']} 仓", f"{stats['n_lineage_rows']} 条谱系行",
                  f"{stats['n_gap_anchors']} 个"):
        assert token in report, f"report does not carry {token!r} from stats.json"
    hull = stats["probe"]["hull_analytic_check"]
    assert str(hull["got"]) in report and str(hull["rel_error_pct"]) in report
    # the <3% verdict must be stated explicitly, per the reporting discipline
    assert "是否 <3%" in report


def test_manifest_pins_every_artifact():
    manifest = _json("manifest.json")
    assert manifest["artifacts"], "manifest is empty"
    for rel, sha in manifest["artifacts"].items():
        p = ART / rel
        assert p.is_file(), rel
        assert hashlib.sha256(p.read_bytes()).hexdigest() == sha, f"{rel}: hash drift"
    for rel, sha in manifest["docs"].items():
        p = ROOT / rel
        assert hashlib.sha256(p.read_bytes()).hexdigest() == sha, f"{rel}: hash drift"


def test_machine_blocks_obey_the_ingestion_contract():
    text = (ART / "machine_blocks.txt").read_text()
    inv = re.findall(r"^\[OSS-INV\]$(.*?)^\[/OSS-INV\]$", text, re.M | re.S)
    sol = re.findall(r"^\[OSS-SOL\]$(.*?)^\[/OSS-SOL\]$", text, re.M | re.S)
    assert inv and sol, "machine blocks missing"
    for body in inv:
        for key in ("repo:", "lang:", "license:", "tier:", "commit:", "summary:"):
            assert key in body, f"OSS-INV missing {key}"
    for body in sol:
        assert "transfer:" in body
        val = re.search(r"^transfer:\s*(\S+)", body, re.M).group(1)
        assert val in TRANSFERS, f"transfer vocabulary drift: {val}"
        assert re.search(r"^slug:\s*GAP-CAD-\d{3}$", body, re.M), "anchor namespace drift"
        files = re.search(r"^files:\s*(.+)$", body, re.M).group(1)
        assert POINTER_RE.match(files.split(",")[0].strip()), files