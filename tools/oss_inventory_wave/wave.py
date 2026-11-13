#!/usr/bin/env python3
"""Open-source technology-inventory *wave* engine (domain-agnostic).

Why this exists
---------------
Every inventory wave in this ecosystem (LBM, FVM, CAD, SPH, BEM, DG, ...) has so far
re-implemented the same protocol from scratch, so the waves drifted: one wave verifies
upstream metadata live, another quotes remembered tags; one wave attaches file:line
pointers to feature claims, another attaches prose.  This module is the protocol as
*code*, with **zero domain semantics**: every domain-specific fact (repository pool,
function points, gap anchors, narrative) lives in ``domains/<name>/*.json`` and is
injected at run time.  ``domains/dg`` is the first instantiation.

Protocol (what a wave artifact must satisfy)
-------------------------------------------
1. ``survey``  - for every pool entry, fetch live GitHub metadata and the recursive
   file tree at the default branch; classify the licence from the API-reported SPDX id
   (measured, never remembered).  Unreachable/rate-limited entries are recorded with
   ``status: failed`` plus the error text, and a pool entry flagged ``negative_control``
   exists precisely to prove that failures are recorded instead of being dropped.
2. feature evidence - "repository X implements feature F" is derived from X's *own*
   file names: each function point carries a path signature and the survey counts
   matching paths (tests/docs/vendor excluded), keeping up to three sample paths.
3. ``build``   - assemble ``registry.json`` (tiers, evidence policy, counts, per-repo
   measured metadata, per-feature path hits) and merge the curated gap anchors after
   re-verifying their *target-repository* evidence with a fresh grep at a pinned commit.
   Markdown views are rendered from the JSON only - never hand edited.
4. ``verify``  - re-run the evidence greps and the markdown round-trip; a gap whose
   evidence no longer holds is reported as a failure, not silently kept.
5. ``plot``    - coverage heat map of pool x function point, from the registry only.

Usage
-----
    python tools/oss_inventory_wave/wave.py survey --domain dg [--budget 40] [--refresh]
    python tools/oss_inventory_wave/wave.py build  --domain dg [--target /path/TensorDG]
    python tools/oss_inventory_wave/wave.py verify --domain dg
    python tools/oss_inventory_wave/wave.py plot   --domain dg
    python tools/oss_inventory_wave/wave.py all    --domain dg --target /path/TensorDG
"""
from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import json
import os
import re
import subprocess
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TOOL = Path(__file__).resolve().parent
DOMAINS = TOOL / "domains"
API = "https://api.github.com"

#: Directories that hold *users* of a library, not the library: keeping them would let
#: ``tests/test_booleans.py`` masquerade as a boolean kernel implementation.
EXCLUDE_SEG = {
    "test", "tests", "testing", "doc", "docs", "documentation", "example", "examples",
    "sample", "samples", "demo", "demos", "benchmark", "benchmarks", "third_party",
    "thirdparty", "3rdparty", "extern", "external", "submodules", "node_modules",
    "vendor", ".git", "ci", ".github", "contrib", "cmake", "build",
}
SOURCE_EXT = {
    ".py", ".h", ".hpp", ".hh", ".hxx", ".cpp", ".cxx", ".cc", ".c", ".inl", ".ipp",
    ".txx", ".ts", ".js", ".rs", ".java", ".cs", ".f90", ".f", ".f95", ".jl", ".go",
    ".m", ".mm", ".cu", ".cl", ".glsl", ".hip",
}


# --------------------------------------------------------------------------- #
# domain configuration
# --------------------------------------------------------------------------- #
def load_domain(name: str) -> dict:
    ddir = DOMAINS / name
    if not ddir.is_dir():
        raise SystemExit(f"unknown domain {name!r}: {ddir} does not exist")
    cfg = json.loads((ddir / "wave.json").read_text())
    cfg["_dir"] = ddir
    cfg["_name"] = name
    cfg["_pool"] = json.loads((ddir / "pool.json").read_text())["repos"]
    cfg["_gaps"] = json.loads((ddir / "gaps.json").read_text())["gaps"]
    pairings = json.loads((ddir / "gaps.json").read_text()).get("pairings", {})
    cfg["_pairings"] = pairings
    cfg["_docs"] = ROOT / cfg["docs_dir"]
    cfg["_art"] = ROOT / cfg["artifacts_dir"]
    return cfg


def _now() -> str:
    return _dt.datetime.now(_dt.timezone.utc).replace(microsecond=0).isoformat()


# --------------------------------------------------------------------------- #
# GitHub API access (cached, budget aware, failure recording)
# --------------------------------------------------------------------------- #
def _cache_path(cfg: dict, key: str) -> Path:
    digest = hashlib.sha256(key.encode()).hexdigest()[:16]
    return cfg["_art"] / "raw" / "api_cache" / f"{digest}.json"


def gh_get(cfg: dict, path: str, refresh: bool = False) -> tuple[dict | None, str | None]:
    """GET api.github.com<path>; returns (payload, error).  Cached on disk."""
    cache = _cache_path(cfg, path)
    if cache.exists() and not refresh:
        blob = json.loads(cache.read_text())
        return blob.get("payload"), blob.get("error")
    req = urllib.request.Request(
        API + path,
        headers={
            "Accept": "application/vnd.github+json",
            "User-Agent": "oss-inventory-wave/1.0",
        },
    )
    payload, error = None, None
    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            payload = json.load(resp)
    except urllib.error.HTTPError as exc:
        error = f"HTTP {exc.code}"
        if exc.code == 403:
            reset = exc.headers.get("X-RateLimit-Reset")
            error += f" rate_or_permission_limited reset={reset}"
    except Exception as exc:  # noqa: BLE001 - network failure must be recorded, not raised
        error = f"{type(exc).__name__}: {exc}"
    cache.parent.mkdir(parents=True, exist_ok=True)
    cache.write_text(json.dumps({"path": path, "fetched_at": _now(),
                                 "payload": payload, "error": error}, sort_keys=True))
    return payload, error


def rate_budget(cfg: dict) -> tuple[int, str]:
    payload, error = gh_get(cfg, "/rate_limit")
    if not payload:
        return 0, error or "rate_limit_unavailable"
    core = payload.get("resources", {}).get("core", payload.get("rate", {}))
    return int(core.get("remaining", 0)), str(core.get("reset", ""))


# --------------------------------------------------------------------------- #
# tree helpers
# --------------------------------------------------------------------------- #
def _keep_path(path: str) -> bool:
    parts = path.split("/")
    if any(p.lower() in EXCLUDE_SEG for p in parts[:-1]):
        return False
    return Path(path).suffix.lower() in SOURCE_EXT


def feature_hits(paths: list[str], signature: str) -> dict:
    rx = re.compile(signature, re.I)
    matched = [p for p in paths if rx.search(p)]
    return {"n": len(matched), "samples": matched[:3]}


# --------------------------------------------------------------------------- #
# survey
# --------------------------------------------------------------------------- #
def survey(cfg: dict, budget: int = 60, refresh: bool = False) -> dict:
    raw = cfg["_art"] / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    remaining, reset = rate_budget(cfg)
    started = time.time()
    done = []
    for entry in cfg["_pool"]:
        rid, repo = entry["id"], entry["repo"]
        if entry.get("role") == "negative_control":
            # fetched on purpose: the failure must end up in the registry
            pass
        meta, err = gh_get(cfg, f"/repos/{repo}", refresh=refresh)
        if meta is None:
            if entry.get("alt_host"):
                # honest census: an upstream that is not on GitHub is recorded as such
                err = f"not_on_github (upstream hosted at {entry['alt_host']}); {err}"
            out = {"schema": "oss-inventory-wave/raw-v1", "id": rid, "repo": repo,
                   "tier": entry["tier"], "role": entry.get("role", "upstream"),
                   "fetched_at": _now(), "source": "github_api", "status": "failed",
                   "error": err, "metadata": None, "features": {},
                   "alt_host": entry.get("alt_host")}
            (raw / f"{rid}.json").write_text(json.dumps(out, indent=1, sort_keys=True))
            done.append((rid, "failed", err))
            continue
        tree_payload, tree_err = gh_get(cfg, f"/repos/{repo}/git/trees/{meta['default_branch']}?recursive=1",
                                        refresh=refresh)
        paths = []
        tree_status = "ok"
        if tree_payload and tree_payload.get("tree") is not None:
            paths = [n["path"] for n in tree_payload["tree"] if n.get("type") == "blob"]
        else:
            tree_status = "failed"
        kept = [p for p in paths if _keep_path(p)]
        feats = {f["id"]: feature_hits(kept, f["signature"]) for f in cfg["features"]}
        out = {
            "schema": "oss-inventory-wave/raw-v1", "id": rid, "repo": repo,
            "tier": entry["tier"], "role": entry.get("role", "upstream"),
            "fetched_at": _now(), "source": "github_api", "status": "ok",
            "metadata": {
                "full_name": meta["full_name"],
                "description": meta.get("description"),
                "stars": meta.get("stargazers_count"),
                "forks": meta.get("forks_count"),
                "language": meta.get("language"),
                "license_spdx": (meta.get("license") or {}).get("spdx_id"),
                "default_branch": meta.get("default_branch"),
                "pushed_at": meta.get("pushed_at"),
                "archived": meta.get("archived"),
                "html_url": meta.get("html_url"),
            },
            "tree": {"status": tree_status, "error": tree_err,
                     "n_blobs": len(paths), "n_source_paths": len(kept)},
            "features": feats,
        }
        (raw / f"{rid}.json").write_text(json.dumps(out, indent=1, sort_keys=True))
        done.append((rid, "ok", None))
        if time.time() - started > budget:
            print(f"  survey stopped by budget after {len(done)} entries "
                  f"({sum(1 for _, s, _ in done if s == 'ok')} ok)")
            break
        remaining, reset = remaining - 2, reset  # metadata + tree per entry
        if remaining <= 2:
            print(f"  survey stopped: unauthenticated quota exhausted (reset={reset})")
            break
    return {"entries": len(done), "ok": sum(1 for _, s, _ in done if s == "ok"),
            "failed": sum(1 for _, s, _ in done if s != "ok")}


# --------------------------------------------------------------------------- #
# target-repository evidence (grep at a pinned commit)
# --------------------------------------------------------------------------- #
def target_commit(cfg: dict, target: Path | None) -> str:
    if target and (target / ".git").exists():
        try:
            sha = subprocess.run(["git", "-C", str(target), "rev-parse", "HEAD"],
                                 capture_output=True, text=True, check=True).stdout.strip()
            if sha:
                return sha
        except Exception:  # noqa: BLE001
            pass
    return cfg["target_commit_pinned"]


def grep_target(target: Path, probes: list[dict], scope: str) -> list[dict]:
    files = []
    for dirpath, dirnames, filenames in os.walk(target / scope):
        dirnames[:] = [d for d in dirnames if d not in {"__pycache__", ".git"}]
        for name in sorted(filenames):
            if Path(name).suffix.lower() in {".py", ".h", ".hpp", ".cpp", ".c", ".f90", ".jl"}:
                files.append(Path(dirpath) / name)
    results = []
    for probe in probes:
        rx = re.compile(probe["pattern"], re.I)
        samples, hits = [], 0
        for path in files:
            try:
                text = path.read_text(errors="replace")
            except OSError:
                continue
            for lineno, line in enumerate(text.splitlines(), 1):
                if rx.search(line):
                    hits += 1
                    if len(samples) < 3:
                        rel = path.relative_to(target).as_posix()
                        samples.append(f"{rel}:{lineno}: {line.strip()[:140]}")
        expect = probe["expect"]
        satisfied = (hits > 0) if expect == "present" else (hits == 0)
        results.append({"pattern": probe["pattern"], "expect": expect, "hits": hits,
                        "samples": samples, "satisfied": satisfied,
                        "note": probe.get("note", "")})
    return results


def verify_gaps(cfg: dict, target_override: str | None = None) -> tuple[list[dict], str, str]:
    """Re-grep the target repo; returns (gaps with evidence, commit, source)."""
    target = None
    for cand in [target_override, os.environ.get(cfg["target_path_env"]), cfg.get("target_path_default")]:
        if cand and Path(cand, cfg["target_scope"].split("/")[0]).is_dir():
            target = Path(cand)
            break
    commit = target_commit(cfg, target)
    source = "live_grep"
    if target is None:
        source = "committed_evidence"
        committed = cfg["_docs"] / "evidence" / "target_grep.json"
        if committed.exists():
            blob = json.loads(committed.read_text())
            if blob.get("commit") != commit:
                print(f"  WARNING: committed evidence pinned at {blob.get('commit')[:12]}, "
                      f"wave.json pins {commit[:12]}")
            return blob["gaps"], blob.get("commit", commit), source
        raise SystemExit("no target checkout and no committed evidence: pass --target")
    gaps = []
    for gap in cfg["_gaps"]:
        probes = grep_target(target, gap["target_evidence"], cfg["target_scope"])
        gaps.append({"id": gap["id"], "probes": probes,
                     "satisfied": all(p["satisfied"] for p in probes)})
    return gaps, commit, source


# --------------------------------------------------------------------------- #
# build
# --------------------------------------------------------------------------- #
def licence_verdict(cfg: dict, metadata: dict | None) -> tuple[str, str]:
    """Map a measured SPDX id onto (licence_class, portability verdict).

    The verdict is the conservative reading: code from a strong-copyleft upstream may
    be read for algorithm design but never copied into this ecosystem's repositories.
    """
    policy = cfg["licence_policy"]
    if not metadata:
        return "unknown", "no_metadata"
    spdx = metadata.get("license_spdx") or "None"
    if spdx in policy["port_ok"]:
        return "permissive", "port_ok"
    if spdx in policy["adapt_allowed_file_level_copyleft"]:
        return "weak_copyleft", "adapt_at_file_level"
    if spdx in policy["read_only_copyleft"]:
        return "strong_copyleft", "read_only_reimplement"
    return "unresolved", "read_only_until_resolved"


def _raw_entries(cfg: dict) -> list[dict]:
    raw = cfg["_art"] / "raw"
    out = []
    for entry in cfg["_pool"]:
        path = raw / f"{entry['id']}.json"
        if path.exists():
            out.append(json.loads(path.read_text()))
    return out


def build(cfg: dict, target: str | None = None) -> dict:
    entries = _raw_entries(cfg)
    evidence, commit, ev_source = verify_gaps(cfg, target)
    by_id = {e["id"]: e for e in entries}
    ev_by_id = {g["id"]: g for g in evidence}

    repos = []
    for entry in cfg["_pool"]:
        e = by_id.get(entry["id"])
        rec = {
            "id": entry["id"], "repo": entry["repo"], "tier": entry["tier"],
            "role": entry.get("role", "upstream"), "focus": entry["focus"],
            "why": entry["why"], "gaps": entry.get("gaps", []),
            "status": (e or {}).get("status", "not_surveyed"),
            "error": (e or {}).get("error"),
            "metadata": (e or {}).get("metadata"),
            "tree": (e or {}).get("tree"),
            "feature_paths": (e or {}).get("features", {}),
        }
        rec["licence_class"], rec["portability"] = licence_verdict(cfg, rec["metadata"])
        repos.append(rec)

    gaps = []
    for gap in cfg["_gaps"]:
        ev = ev_by_id[gap["id"]]["probes"] if gap["id"] in ev_by_id else []
        refs = []
        for rid in gap["refs"]:
            rec = by_id.get(rid)
            if rec is None:
                continue
            refs.append({"id": rid, "repo": rec["repo"], "tier": rec["tier"],
                         "status": rec["status"],
                         "feature_paths": rec.get("features", {})})
        gaps.append({
            "id": gap["id"], "title": gap["title"], "priority": gap["priority"],
            "effort": gap["effort"], "feature_points": gap["feature_points"],
            "status": gap["status"], "current_state": gap["current_state"],
            "target_evidence": ev,
            "evidence_satisfied": all(p["satisfied"] for p in ev) if ev else False,
            "refs": refs, "acceptance": gap["acceptance"], "upstream_note": gap["upstream_note"],
            "external_anchor": gap["external_anchor"], "impact": gap["impact"],
        })

    up = [r for r in repos if r["role"] == "upstream"]
    tiers: dict[str, int] = {}
    for r in up:
        tiers[r["tier"]] = tiers.get(r["tier"], 0) + 1
    registry = {
        "schema_version": "1.0",
        "generated_by": "tools/oss_inventory_wave/wave.py build",
        "generated_at": _now(),
        "domain": cfg["domain"],
        "feature_slug": cfg["feature_slug"],
        "target_repo": cfg["target_repo"],
        "target_commit": commit,
        "target_evidence_source": ev_source,
        "target_scope": cfg["target_scope"],
        "tiers": cfg["tiers"],
        "evidence_policy": cfg["evidence_policy"],
        "counts": {
            "pool_entries": len(repos),
            "upstream_repos": len(up),
            "by_tier": tiers,
            "surveyed_ok": sum(1 for r in up if r["status"] == "ok"),
            "survey_failed": sum(1 for r in up if r["status"] == "failed"),
            "not_surveyed": sum(1 for r in up if r["status"] == "not_surveyed"),
            "negative_controls": sum(1 for r in repos if r["role"] == "negative_control"),
            "negative_controls_recorded_failed": sum(
                1 for r in repos if r["role"] == "negative_control" and r["status"] == "failed"),
            "feature_points": len(cfg["features"]),
            "by_portability": {k: sum(1 for r in up if r["portability"] == k) for k in
                               ("port_ok", "adapt_at_file_level", "read_only_reimplement",
                                "read_only_until_resolved", "no_metadata")},
            "gaps": len(gaps),
            "gaps_evidence_satisfied": sum(1 for g in gaps if g["evidence_satisfied"]),
        },
        "features": cfg["features"],
        "repos": repos,
        "gaps": gaps,
    }
    docs = cfg["_docs"]
    docs.mkdir(parents=True, exist_ok=True)
    (docs / "registry.json").write_text(json.dumps(registry, indent=1, sort_keys=False) + "\n")
    (docs / "evidence").mkdir(exist_ok=True)
    (docs / "evidence" / "target_grep.json").write_text(json.dumps(
        {"schema": "oss-inventory-wave/target-grep-v1", "domain": cfg["domain"],
         "target_repo": cfg["target_repo"], "commit": commit, "scope": cfg["target_scope"],
         "source": ev_source, "gaps": evidence}, indent=1, sort_keys=True) + "\n")
    (docs / "registry.md").write_text(render_registry_md(registry))
    (docs / "gap-proposals.md").write_text(render_gaps_md(registry))
    template = cfg["_dir"] / "readme_template.md"
    if template.exists():
        (docs / "README.md").write_text(render_readme(registry, template.read_text()))
    return registry


# --------------------------------------------------------------------------- #
# markdown rendering (deterministic, JSON only)
# --------------------------------------------------------------------------- #
TIER_TITLES = {"T1": "T1 可吸收", "T2": "T2 参考设计", "T3": "T3 背景"}


def _tier_table(registry: dict) -> str:
    rows = ["| Tier | 数量 | 定义 |", "|---|---:|---|"]
    for key in ("T1", "T2", "T3"):
        rows.append(f"| `{key}` {TIER_TITLES[key]} | {registry['counts']['by_tier'].get(key, 0)} "
                    f"| {registry['tiers'][key]} |")
    return "\n".join(rows)


def render_registry_md(registry: dict) -> str:
    c = registry["counts"]
    lines = [
        f"# {registry['domain']} 开源技术盘点 registry（{registry['feature_slug']}）",
        "",
        "> 程序化生成，勿手改。源：`docs/oss_inventory_dg/registry.json` "
        f"（`{registry['generated_by']}`，生成时间 {registry['generated_at']}）。",
        f"> 目标仓 `{registry['target_repo']}` @ `{registry['target_commit']}`"
        f"（证据来源：{registry['target_evidence_source']}）。",
        "",
        "## 1 计数",
        "",
        "| 项 | 值 |",
        "|---|---:|",
    ]
    for key in ("pool_entries", "upstream_repos", "surveyed_ok", "survey_failed",
                "not_surveyed", "negative_controls", "negative_controls_recorded_failed",
                "feature_points", "gaps", "gaps_evidence_satisfied"):
        lines.append(f"| `{key}` | {c[key]} |")
    lines += ["", "## 2 Tier 分级", "", _tier_table(registry), "",
              "## 3 证据政策", ""]
    for key, val in registry["evidence_policy"].items():
        lines.append(f"- **`{key}`**：{val}")
    lines += ["", "## 4 仓位（普查结果）", "",
              "| id | repo | tier | status | stars | licence(SPDX) | 可用性判定 | pushed_at | 源码文件数 |",
              "|---|---|---|---|---:|---|---|---|---:|"]
    for r in registry["repos"]:
        m = r["metadata"] or {}
        stars = m.get("stars", "—")
        lic = m.get("license_spdx") or "—"
        pushed = (m.get("pushed_at") or "—")[:10]
        nsrc = (r["tree"] or {}).get("n_source_paths", "—")
        lines.append(f"| `{r['id']}` | `{r['repo']}` | {r['tier']} | {r['status']} "
                     f"| {stars} | {lic} | `{r['portability']}` | {pushed} | {nsrc} |")
    lines += ["", "## 5 功能点 × 仓覆盖矩阵", "",
              "单元格=该仓源码路径中命中此功能点签名路径数（tests/docs/vendor 已剔除）；"
              "`—`=该仓未取到树。", ""]
    feats = registry["features"]
    header = "| repo | " + " | ".join(f["id"] for f in feats) + " |"
    lines += [header, "|---|" + "---|" * len(feats)]
    for r in registry["repos"]:
        if r["role"] != "upstream":
            continue
        cells = []
        for f in feats:
            hit = (r["feature_paths"] or {}).get(f["id"])
            cells.append("—" if not r["feature_paths"] else str(hit["n"]))
        lines.append(f"| `{r['id']}` | " + " | ".join(cells) + " |")
    lines += ["", "功能点定义：", ""]
    for f in feats:
        lines.append(f"- **`{f['id']}`** {f['title']} — 路径签名 `{f['signature']}`")
    failed = [r for r in registry["repos"] if r["status"] == "failed"]
    lines += ["", "## 6 取数失败（如实记录，不静默丢弃）", ""]
    if failed:
        for r in failed:
            lines.append(f"- `{r['id']}` (`{r['repo']}`) — `{r['error']}`")
    else:
        lines.append("- 无")
    return "\n".join(lines) + "\n"


def render_gaps_md(registry: dict) -> str:
    order = {"P0": 0, "P1": 1, "P2": 2}
    gaps = sorted(registry["gaps"], key=lambda g: (order.get(g["priority"], 9), g["id"]))
    lines = [f"# {registry['target_repo']} 缺口锚点提议（GAP-{registry['domain']}-*，本波 {len(gaps)} 条）",
             "",
             "> 程序化生成，勿手改。源：`registry.json` 的 `gaps` 段"
             f"（生成时间 {registry['generated_at']}）。",
             f"> 现状证据一律在本仓可 grep：`{registry['target_repo']}` @ `{registry['target_commit']}`，"
             f"scope=`{registry['target_scope']}`。",
             "> 验收判据必须是机器可裁的数值门；外部锚必须落在文献/解析解，不得用本仓实现自证。",
             "", "## 排序总表", "",
             "| # | 锚点 | 优先级 | 参考实现（Tier） | 工作量 | 现状证据 | 判据要点 |",
             "|---|---|---|---|---|---|---|"]
    for i, g in enumerate(gaps, 1):
        refs = ", ".join(f"`{r['id']}`({r['tier']})" for r in g["refs"][:4])
        short = g["acceptance"].split("；")[0]
        ok = "✅" if g["evidence_satisfied"] else "❌"
        lines.append(f"| {i} | `{g['id']}` | {g['priority']} | {refs} | {g['effort']} "
                     f"| {ok} {g['status']} | {short} |")
    lines += ["", "## 逐条详述", ""]
    for g in gaps:
        lines += [f"### `{g['id']}` — {g['title']}", "",
                  f"- **优先级**：{g['priority']}；**工作量**：{g['effort']}；"
                  f"**功能点**：{', '.join('`'+f+'`' for f in g['feature_points'])}",
                  f"- **现状（{registry['target_repo']}）**：{g['current_state']}",
                  "- **现状证据（grep @ "
                  f"{registry['target_commit'][:12]}）**："]
        for p in g["target_evidence"]:
            flag = "命中" if p["hits"] else "零命中"
            lines.append(f"  - `{p['pattern']}` 期望={p['expect']} → {flag} {p['hits']} 处 "
                         f"{'✔' if p['satisfied'] else '✘'}")
            for s in p["samples"]:
                lines.append(f"    - `{s}`")
        lines.append("- **参考实现**：")
        for r in g["refs"]:
            lanes = []
            for f in g["feature_points"]:
                hit = (r.get("feature_paths") or {}).get(f)
                if hit and hit["n"]:
                    lanes.append(f"{f}×{hit['n']} {hit['samples'][0]}")
            lane = "；".join(lanes) if lanes else "本次未取到该功能点的路径命中"
            lines.append(f"  - `{r['repo']}` ({r['tier']}, {r['status']}) — {lane}")
        lines += [f"- **验收判据（机器可裁）**：{g['acceptance']}",
                  f"- **上游要点**：{g['upstream_note']}",
                  f"- **外部锚**：{g['external_anchor']}",
                  f"- **影响**：{g['impact']}", ""]
    return "\n".join(lines) + "\n"


def render_readme(registry: dict, template: str) -> str:
    east = {"by_tier": registry["counts"]["by_tier"],
            "by_portability": registry["counts"]["by_portability"],
            "counts": registry["counts"], "domain": registry["domain"],
            "target_repo": registry["target_repo"],
            "target_commit": registry["target_commit"],
            "target_commit_short": registry["target_commit"][:12],
            "feature_slug": registry["feature_slug"],
            "generated_at": registry["generated_at"],
            "n_gaps_p0": sum(1 for g in registry["gaps"] if g["priority"] == "P0"),
            "pairs": registry["_pairings"] if "_pairings" in registry else {}}
    out = template
    for key, val in east.items():
        if isinstance(val, dict):
            for k2, v2 in val.items():
                out = out.replace("{{" + f"{key}.{k2}" + "}}", str(v2))
        else:
            out = out.replace("{{" + key + "}}", str(val))
    return out


# --------------------------------------------------------------------------- #
# plot
# --------------------------------------------------------------------------- #
def plot(cfg: dict) -> list[Path]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    registry = json.loads((cfg["_docs"] / "registry.json").read_text())
    feats = [f["id"] for f in registry["features"]]
    rows = [r for r in registry["repos"] if r["role"] == "upstream"]
    mat = np.zeros((len(rows), len(feats)))
    for i, r in enumerate(rows):
        for j, f in enumerate(feats):
            mat[i, j] = ((r["feature_paths"] or {}).get(f) or {}).get("n", 0)
    figdir = cfg["_docs"] / "figures"
    figdir.mkdir(parents=True, exist_ok=True)
    paths = []
    fig, ax = plt.subplots(figsize=(10.5, 0.34 * len(rows) + 1.9))
    shown = np.log1p(mat)
    im = ax.imshow(shown, cmap="viridis", aspect="auto")
    ax.set_xticks(range(len(feats)))
    ax.set_xticklabels(feats, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([f"{r['id']} ({r['tier']})" for r in rows], fontsize=7)
    for i in range(len(rows)):
        for j in range(len(feats)):
            if mat[i, j]:
                ax.text(j, i, int(mat[i, j]), ha="center", va="center", fontsize=6,
                        color="white" if shown[i, j] < shown.max() * 0.6 else "black")
    ax.set_title(f"{registry['domain']} inventory: function-point coverage "
                 f"(source-path hits per repo, log colour)")
    fig.colorbar(im, ax=ax, label="log1p(paths)")
    fig.tight_layout()
    p = figdir / "coverage_heatmap.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    paths.append(p)

    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    tiers = ["T1", "T2", "T3"]
    vals = [registry["counts"]["by_tier"].get(t, 0) for t in tiers]
    ok = [sum(1 for r in rows if r["tier"] == t and r["status"] == "ok") for t in tiers]
    ax.bar(tiers, vals, color="#4c72b0", label="pool")
    ax.bar(tiers, ok, color="#55a868", label="surveyed (live metadata)")
    for i, (v, o) in enumerate(zip(vals, ok)):
        ax.text(i, v + 0.1, f"{o}/{v}", ha="center", fontsize=9)
    ax.set_ylabel("repositories")
    ax.set_title(f"{registry['domain']} pool by tier — surveyed / planned")
    ax.legend(fontsize=8)
    fig.tight_layout()
    p2 = figdir / "tier_survey.png"
    fig.savefig(p2, dpi=150)
    plt.close(fig)
    paths.append(p2)
    return paths


# --------------------------------------------------------------------------- #
# verify
# --------------------------------------------------------------------------- #
def verify(cfg: dict, target: str | None = None) -> dict:
    import tempfile

    registry = json.loads((cfg["_docs"] / "registry.json").read_text())
    problems = []
    up = [r for r in registry["repos"] if r["role"] == "upstream"]
    if len(up) < 20:
        problems.append(f"census has {len(up)} upstream repos, acceptance requires >=20")
    for r in registry["repos"]:
        if r["tier"] not in registry["tiers"]:
            problems.append(f"{r['id']}: invalid tier {r['tier']}")
    for g in registry["gaps"]:
        if len(g["refs"]) < 2:
            problems.append(f"{g['id']}: fewer than 2 reference implementations")
        if not g["evidence_satisfied"]:
            problems.append(f"{g['id']}: target evidence probes unsatisfied")
        if not re.search(r"\d", g["acceptance"]):
            problems.append(f"{g['id']}: acceptance criterion carries no numeric threshold")
        if registry["target_repo"].lower() in g["external_anchor"].lower():
            problems.append(f"{g['id']}: external anchor cites the target repo itself")
    if registry["counts"]["negative_controls"] and not \
            registry["counts"]["negative_controls_recorded_failed"]:
        problems.append("negative control was not recorded as failed")
    # markdown round-trip: views must equal a fresh render of the JSON
    if (cfg["_docs"] / "registry.md").read_text() != render_registry_md(registry):
        problems.append("registry.md is not reproducible from registry.json (hand edited?)")
    if (cfg["_docs"] / "gap-proposals.md").read_text() != render_gaps_md(registry):
        problems.append("gap-proposals.md is not reproducible from registry.json (hand edited?)")
    # live evidence re-check when a checkout is available
    evidence, commit, source = verify_gaps(cfg, target)
    if not all(g["satisfied"] for g in evidence):
        problems.append("re-grep of the target repository no longer satisfies the probes")
    if source == "live_grep" and commit != registry["target_commit"]:
        problems.append(f"target commit moved: registry {registry['target_commit'][:12]} "
                        f"vs live {commit[:12]} (re-run build)")
    return {"problems": problems, "evidence_source": source, "commit": commit}


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("command", choices=["survey", "build", "plot", "verify", "all"])
    ap.add_argument("--domain", default="dg")
    ap.add_argument("--budget", type=int, default=60, help="survey wall-clock budget (s)")
    ap.add_argument("--refresh", action="store_true", help="bypass the API cache")
    ap.add_argument("--target", default=None, help="target repository checkout")
    args = ap.parse_args(argv)
    cfg = load_domain(args.domain)

    if args.command in ("survey", "all"):
        print(f"survey {cfg['feature_slug']}: {len(cfg['_pool'])} pool entries")
        print(" ", survey(cfg, budget=args.budget, refresh=args.refresh))
    if args.command in ("build", "all"):
        reg = build(cfg, args.target)
        print(f"build: {reg['counts']['upstream_repos']} upstream repos, "
              f"{reg['counts']['gaps']} gaps, docs -> {cfg['docs_dir']}")
    if args.command in ("plot", "all"):
        for p in plot(cfg):
            print(f"figure {p.relative_to(ROOT)} ({p.stat().st_size} B)")
    if args.command in ("verify", "all"):
        res = verify(cfg, args.target)
        for p in res["problems"]:
            print("PROBLEM:", p)
        print("verify:", "OK" if not res["problems"] else f"{len(res['problems'])} problem(s)",
              f"(evidence source={res['evidence_source']}, commit={res['commit'][:12]})")
        return 0 if not res["problems"] else 1
    return 0


if __name__ == "__main__":
    sys.exit(main())