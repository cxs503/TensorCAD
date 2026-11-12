"""Repository patrol: upstream diff-set + commit-theme evolution inventory.

Read-only. Every number in the generated report comes from the JSON artifacts this
module writes; nothing is typed by hand.

Two lenses are collected:

* **upstream lens** -- does a remote upstream exist, and if so how does the local
  ``HEAD`` differ from it (ahead/behind, changed files, live-vs-fetched staleness)?
* **evolution lens** -- the last ``window`` *first-parent* commits of the target
  ref, i.e. what actually landed on the integration branch. For true merges the
  content is measured against the first parent (exactly what the merge brought
  in); the packaged side-branch commits are recorded as well. Commits whose
  subject *claims* a merge but which have a single parent are reported as
  ``squash`` landings -- they carry no ancestry link back to the mission branch.

The theme/risk vocabulary lives in ``themes.json`` next to this file so that
classification is data, not code.

Usage:
    python tools/repo_patrol/patrol.py --ref origin/main --window 15 \
        --outdir artifacts/repo-patrol/upstream-diffset-v1
    python tools/repo_patrol/patrol.py --offline     # skip the live ls-remote probe
"""

# ruff: noqa: E501  (long lines are mostly Chinese report prose)
from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

SCHEMA = "repo-patrol/2"
DEFAULT_CONFIG = Path(__file__).with_name("themes.json")
VOLATILE_KEYS = ("generated_at", "repo_root")


class PatrolError(RuntimeError):
    """A git command failed or produced output we cannot interpret."""


# --------------------------------------------------------------------------- #
# git plumbing
# --------------------------------------------------------------------------- #
def run_git(repo: Path | str, *args: str, check: bool = True) -> str:
    proc = subprocess.run(
        ["git", "-C", str(repo), *args],
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    if check and proc.returncode != 0:
        raise PatrolError(f"git {' '.join(args)} failed: {proc.stderr.strip()}")
    return proc.stdout


def _count(repo: Path, rev_range: str) -> int:
    return int(run_git(repo, "rev-list", "--count", rev_range).strip() or 0)


def _resolve(repo: Path, ref: str) -> str | None:
    return run_git(repo, "rev-parse", "--verify", "--quiet", ref, check=False).strip() or None


def _parse_numstat(text: str) -> tuple[int, int, list[tuple[str, int | None, int | None]]]:
    """Return (insertions, deletions, [(path, ins, del)]); binary files count as None."""
    ins_total = dele_total = 0
    rows: list[tuple[str, int | None, int | None]] = []
    for line in text.splitlines():
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        raw_ins, raw_del, path = parts[0], parts[1], "\t".join(parts[2:])
        if raw_ins == "-" or raw_del == "-":
            rows.append((path, None, None))
            continue
        ins, dele = int(raw_ins), int(raw_del)
        ins_total += ins
        dele_total += dele
        rows.append((path, ins, dele))
    return ins_total, dele_total, rows


_SHORTSTAT = re.compile(r"(\d+) files? changed(?:, (\d+) insertions?\(\+\))?(?:, (\d+) deletions?\(-\))?")


def _parse_shortstat(text: str) -> dict:
    line = next((ln for ln in text.splitlines() if "changed" in ln), "")
    m = _SHORTSTAT.search(line)
    if not m:
        return {"files": 0, "insertions": 0, "deletions": 0}
    return {
        "files": int(m.group(1)),
        "insertions": int(m.group(2) or 0),
        "deletions": int(m.group(3) or 0),
    }


def _diff_numstat(repo: Path, base: str, head: str):
    return _parse_numstat(run_git(repo, "diff", "--numstat", "--no-renames", base, head))


def _diff_shortstat(repo: Path, base: str, head: str) -> dict:
    return _parse_shortstat(run_git(repo, "diff", "--shortstat", "--no-renames", base, head))


# --------------------------------------------------------------------------- #
# classification (data-driven)
# --------------------------------------------------------------------------- #
def load_config(path: Path | str = DEFAULT_CONFIG) -> dict:
    cfg = json.loads(Path(path).read_text(encoding="utf-8"))
    for key in ("areas", "topics", "risk_rules"):
        if key not in cfg:
            raise PatrolError(f"theme config missing '{key}'")
    return cfg


def classify_area(path: str, cfg: dict) -> str:
    for rule in cfg["areas"]:
        for prefix in rule["prefixes"]:
            if path.startswith(prefix):
                return rule["area"]
    return "other"


_TOKEN_SPLIT = re.compile(r"[^0-9a-z]+")


def _tokens(text: str) -> list[str]:
    return [t for t in _TOKEN_SPLIT.split(text.lower()) if t]


def classify_topics(subject: str, paths: list[str], cfg: dict) -> list[str]:
    """Multi-label topic assignment.

    Evidence comes from the commit subject plus the *stems* of changed files. Test
    files contribute their subject stem (``test_foo`` -> ``foo``) so that a pure
    test-infrastructure commit is not mistaken for the feature it exercises.
    """
    corpus: list[str] = list(_tokens(subject))
    for path in paths:
        if classify_area(path, cfg) == "repo-meta":
            # Bookkeeping files (CHANGELOG/README) are edited by nearly every
            # landing; their stems would tag every commit with the same topic.
            continue
        stem = path.rsplit("/", 1)[-1].rsplit(".", 1)[0]
        if stem.startswith("test_"):
            stem = stem[len("test_"):]
        corpus.extend(_tokens(stem))
    blob = " ".join(corpus)
    hits = [rule["topic"] for rule in cfg["topics"] if any(kw in blob for kw in rule["keywords"])]
    return hits or ["unclassified"]


def parse_gate_filtered(subject: str, cfg: dict) -> int | None:
    """Files declared dropped by the delivery gate, as self-reported in a subject."""
    rules = cfg["risk_rules"]
    if rules["gate_filtered_pattern"] not in subject:
        return None
    m = re.search(rules["gate_filtered_count_pattern"], subject)
    return int(m.group(1)) if m else 0


def parse_merge_claim(subject: str, cfg: dict) -> dict:
    """What the subject *claims* about provenance (independently of git ancestry)."""
    m = re.search(cfg["risk_rules"]["merge_subject_pattern"], subject)
    if not m:
        return {}
    return {"pr": int(m.group(1)), "mission": f"m{m.group(2)}", "merge_lane": m.group(3)}


def parse_mission(subject: str, cfg: dict) -> str | None:
    m = re.search(cfg["risk_rules"]["mission_pattern"], subject)
    return f"m{m.group(1)}" if m else None


def landing_mode(commit: dict) -> str:
    """How the commit reached the integration branch."""
    if commit["is_merge"]:
        return "true-merge" if commit["merge_claim"] else "unclaimed-merge"
    if commit["merge_claim"]:
        return "squash"
    return "direct"


def commit_risks(commit: dict, cfg: dict) -> list[str]:
    risks: list[str] = []
    areas = commit["areas"]
    if commit["is_merge"]:
        risks.append("merge_commit")
    if commit["landing_mode"] == "squash":
        risks.append("squash_landing_no_ancestry")
    if commit["landing_mode"] == "squash" and commit["files_changed"] == 0:
        risks.append("empty_landing")
    if areas.get("core") and not areas.get("tests"):
        risks.append("core_without_tests")
    if commit["insertions"] + commit["deletions"] > cfg["risk_rules"]["large_churn_lines"]:
        risks.append("large_churn")
    if commit.get("gate_filtered_files"):
        risks.append("gate_filtered_delivery")
    if commit["files_changed"] == 0:
        risks.append("empty_delta")
    if not (areas.get("core") or areas.get("tests") or areas.get("tools")):
        risks.append("no_code_delta")
    generic = set(cfg["risk_rules"].get("low_signal_tokens", []))
    # digits (mission ids, "6 个文件") carry no descriptive signal
    tail = [t for t in _tokens(re.sub(r"^m\d+\s*:", "", commit["subject"])) if not t.isdigit()]
    if not tail or set(tail) <= generic:
        risks.append("low_signal_subject")
    return risks


# --------------------------------------------------------------------------- #
# collectors
# --------------------------------------------------------------------------- #
def collect_upstream(repo: Path, ref: str, head: str, network: bool) -> dict:
    remotes = [r for r in run_git(repo, "remote").split() if r]
    if not remotes:
        return {"verdict": "NO_UPSTREAM", "remotes": [], "note": "no git remote configured"}
    remote = ref.split("/", 1)[0] if "/" in ref else remotes[0]
    if remote not in remotes:
        remote = remotes[0]
    url = run_git(repo, "remote", "get-url", remote, check=False).strip()
    upstream_sha = _resolve(repo, ref)
    live_sha = None
    if network:
        out = run_git(repo, "ls-remote", remote, "HEAD", check=False).strip()
        live_sha = out.split()[0] if out else None
    if upstream_sha is None:
        return {
            "verdict": "NO_UPSTREAM",
            "remotes": remotes,
            "remote": remote,
            "url": url,
            "note": f"ref {ref} not present locally",
        }
    ahead = _count(repo, f"{ref}..{head}")
    behind = _count(repo, f"{head}..{ref}")
    verdict = (
        "IN_SYNC"
        if ahead == 0 and behind == 0
        else "DIVERGED"
        if ahead and behind
        else "AHEAD_ONLY"
        if ahead
        else "BEHIND"
    )
    # "diffset" is directional: what upstream carries that HEAD does not (head -> ref),
    # and separately what HEAD carries that upstream does not (ref -> head).
    diff_files: list[str] = []
    diffstat = {"insertions": 0, "deletions": 0, "files": 0}
    if behind:
        ins, dele, rows = _diff_numstat(repo, head, ref)
        diff_files = [p for p, _, _ in rows]
        diffstat = {"insertions": ins, "deletions": dele, "files": len(rows)}
    our_files: list[str] = []
    our_stat = {"insertions": 0, "deletions": 0, "files": 0}
    if ahead:
        ins, dele, rows = _diff_numstat(repo, ref, head)
        our_files = [p for p, _, _ in rows]
        our_stat = {"insertions": ins, "deletions": dele, "files": len(rows)}
    return {
        "verdict": verdict,
        "remotes": remotes,
        "remote": remote,
        "url": url,
        "upstream_ref": ref,
        "upstream_sha": upstream_sha,
        "head_sha": head,
        "live_remote_sha": live_sha,
        "live_probe": "ls-remote" if network else "skipped",
        "fetched_ref_is_live": (live_sha == upstream_sha) if live_sha else None,
        "ahead": ahead,
        "behind": behind,
        "diffset_files": diff_files,
        "diffset_stat": diffstat,
        "our_delta_files": our_files,
        "our_delta_stat": our_stat,
    }


def collect_commit(repo: Path, sha: str, cfg: dict, max_files: int = 400) -> dict:
    parents = run_git(repo, "rev-list", "--parents", "-n", "1", sha).split()[1:]
    meta = run_git(repo, "show", "-s", "--format=%an%x1f%ad%x1f%s", "--date=short", sha).rstrip("\n")
    author, date, subject = [*meta.split("\x1f"), "", "", ""][:3]
    if parents:
        base = parents[0]
        ins, dele, rows = _diff_numstat(repo, base, sha)
        short = _diff_shortstat(repo, base, sha)
    else:  # root commit
        ins, dele, rows = _parse_numstat(run_git(repo, "show", "--numstat", "--format=", sha))
        short = _parse_shortstat(run_git(repo, "show", "--shortstat", "--format=", sha))
    paths = [p for p, _, _ in rows]
    areas: dict[str, int] = {}
    for path in paths:
        area = classify_area(path, cfg)
        areas[area] = areas.get(area, 0) + 1
    claim = parse_merge_claim(subject, cfg)
    packaged: list[str] = []
    if len(parents) > 1:
        packaged = run_git(repo, "log", "--format=%h %s", f"{parents[0]}..{parents[1]}").splitlines()
    packaged_gate = [parse_gate_filtered(ln, cfg) for ln in packaged]
    commit = {
        "sha": sha,
        "short": sha[:9],
        "author": author,
        "date": date,
        "subject": subject,
        "parents": parents,
        "is_merge": len(parents) > 1,
        "merge_claim": claim,
        "landing_mode": "",
        "files_changed": len(paths),
        "files": paths[:max_files],
        "files_truncated": len(paths) > max_files,
        "insertions": ins,
        "deletions": dele,
        "shortstat_agrees": (
            short["files"] == len(paths)
            and short["insertions"] == ins
            and short["deletions"] == dele
        ),
        "areas": areas,
        "topics": classify_topics(subject, paths, cfg),
        "mission": claim.get("mission") or parse_mission(subject, cfg),
        "pr": claim.get("pr"),
        "merge_lane": claim.get("merge_lane"),
        "gate_filtered_files": parse_gate_filtered(subject, cfg),
        "packaged_commits": packaged[:40],
        "packaged_count": len(packaged),
        "packaged_gate_filtered": sum(g for g in packaged_gate if g),
        "packaged_missions": sorted(
            {m for m in (parse_mission(ln, cfg) for ln in packaged) if m}
        ),
    }
    commit["landing_mode"] = landing_mode(commit)
    commit["risks"] = commit_risks(commit, cfg)
    return commit


def collect_window(repo: Path, ref: str, window: int, cfg: dict) -> dict:
    shas = run_git(repo, "rev-list", "--first-parent", f"--max-count={window}", ref).split()
    commits = [collect_commit(repo, sha, cfg) for sha in shas]
    totals = {
        "commits": len(commits),
        "merge_commits": sum(1 for c in commits if c["is_merge"]),
        "squash_landings": sum(1 for c in commits if c["landing_mode"] == "squash"),
        "empty_landings": sum(1 for c in commits if c["landing_mode"] == "squash" and c["files_changed"] == 0),
        "mission_commits": sum(1 for c in commits if c["mission"]),
        "unique_missions": sorted({c["mission"] for c in commits if c["mission"]}),
        "gate_filtered_files_reported": sum(
            (c["gate_filtered_files"] or 0) + c["packaged_gate_filtered"] for c in commits
        ),
        "insertions": sum(c["insertions"] for c in commits),
        "deletions": sum(c["deletions"] for c in commits),
        "files_touched_sum": sum(c["files_changed"] for c in commits),
        "files_touched_unique": len({p for c in commits for p in c["files"]}),
        "areas": {},
        "topics": {},
        "authors": {},
        "risks": {},
        "landing_modes": {},
    }
    for commit in commits:
        for area, n in commit["areas"].items():
            totals["areas"][area] = totals["areas"].get(area, 0) + n
        for topic in commit["topics"]:
            totals["topics"][topic] = totals["topics"].get(topic, 0) + 1
        totals["authors"][commit["author"]] = totals["authors"].get(commit["author"], 0) + 1
        totals["landing_modes"][commit["landing_mode"]] = (
            totals["landing_modes"].get(commit["landing_mode"], 0) + 1
        )
        for risk in commit["risks"]:
            totals["risks"][risk] = totals["risks"].get(risk, 0) + 1
    span = None
    if shas:
        oldest_parents = run_git(repo, "rev-list", "--parents", "-n", "1", shas[-1]).split()[1:]
        if oldest_parents:
            span = _diff_shortstat(repo, f"{shas[-1]}^", ref)
        else:  # window reaches the root commit: diff against the empty tree
            empty_tree = run_git(repo, "hash-object", "-t", "tree", "/dev/null").strip()
            span = _diff_shortstat(repo, empty_tree, ref)
    return {
        "commits": commits,
        "totals": totals,
        "oldest": shas[-1] if shas else None,
        "span_stat": span,
    }


def collect_crosschecks(repo: Path, ref: str, window: int, window_data: dict) -> list[dict]:
    """Recompute window quantities through an independent git path.

    Each check compares the tool's own aggregation against a git plumbing command
    that does not share its code path. Agreement verifies the inventory mechanics
    (parsing, counting, range semantics); it does not verify the theme labels.
    """
    commits = window_data["commits"]
    shas = [c["sha"] for c in commits]
    checks: list[dict] = []

    def add(name: str, a_label: str, a, b_label: str, b, tol: int = 0) -> None:
        ok = a is not None and b is not None and abs(a - b) <= tol
        checks.append(
            {
                "name": name,
                "method_a": a_label,
                "value_a": a,
                "method_b": b_label,
                "value_b": b,
                "delta": (a - b) if (a is not None and b is not None) else None,
                "match": bool(ok),
            }
        )

    add(
        "window_commit_count",
        "len(rev-list --first-parent)",
        len(commits),
        "rev-list --count --first-parent --max-count",
        int(run_git(repo, "rev-list", "--count", "--first-parent", f"--max-count={window}", ref).strip()),
    )
    if shas:
        no_walk = run_git(repo, "rev-list", "--no-walk=unsorted", "--merges", *shas).split()
        add(
            "window_merge_count",
            "sum(is_merge) from rev-list --parents",
            sum(1 for c in commits if c["is_merge"]),
            "rev-list --no-walk --merges <window shas>",
            len(no_walk),
        )
        add(
            "per_commit_numstat_vs_shortstat",
            "len(commits)",
            len(commits),
            "commits whose --numstat == --shortstat",
            sum(1 for c in commits if c["shortstat_agrees"]),
        )
        parented = [c for c in commits if c["parents"]]
        add(
            "window_insertions_shortstat",
            "sum(numstat insertions, non-root)",
            sum(c["insertions"] for c in parented),
            "sum(diff --shortstat insertions, non-root)",
            sum(_diff_shortstat(repo, f"{c['sha']}^1", c["sha"])["insertions"] for c in parented),
        )
    return checks


def collect_churn_span(repo: Path, ref: str, window_data: dict) -> dict:
    """Window churn vs single span diff.

    These two are *not* expected to be equal: lines added and later removed inside
    the window cancel out in the span diff. The gap therefore measures rewrite
    churn, and is reported as a metric rather than an agreement check.
    """
    totals = window_data["totals"]
    span = window_data["span_stat"] if window_data["span_stat"] else {"insertions": 0, "deletions": 0, "files": 0}
    per_commit = totals["insertions"] + totals["deletions"]
    span_total = span["insertions"] + span["deletions"]
    return {
        "per_commit_lines": per_commit,
        "span_lines": span_total,
        "span_insertions": span["insertions"],
        "span_deletions": span["deletions"],
        "amplification": round(per_commit / span_total, 4) if span_total else None,
        "cancelled_lines": per_commit - span_total,
        "note": "per-commit sums > span diff because intra-window add/remove pairs cancel",
    }


def _days_between(earlier: str, later: str) -> int | None:
    try:
        a = _dt.date.fromisoformat(earlier)
        b = _dt.date.fromisoformat(later)
    except ValueError:
        return None
    return (b - a).days


def collect_unmerged(repo: Path, ref: str, cfg: dict, head_date: str | None = None) -> list[dict]:
    """Remote branches carrying work that is not an ancestor of the target ref.

    Topics are taken from the *branch's own commit subjects*, not from the full
    diff against the target: a long-lived branch touches enough files that a
    path-based corpus saturates every topic and stops discriminating.
    """
    out = run_git(repo, "for-each-ref", "--format=%(refname:short)", "refs/remotes").split()
    branches = [b for b in out if b and not b.endswith("/HEAD")]
    unmerged = []
    for branch in branches:
        ahead = _count(repo, f"{ref}..{branch}")
        if ahead == 0:
            continue
        meta = run_git(repo, "show", "-s", "--format=%ad%x1f%s", "--date=short", branch).rstrip("\n")
        date, _, subject = meta.partition("\x1f")
        subjects = run_git(repo, "log", "--format=%s", "--max-count=40", f"{ref}..{branch}").splitlines()
        ins, dele, rows = _diff_numstat(repo, ref, branch)
        paths = [p for p, _, _ in rows]
        areas: dict[str, int] = {}
        for path in paths:
            area = classify_area(path, cfg)
            areas[area] = areas.get(area, 0) + 1
        unmerged.append(
            {
                "branch": branch,
                "namespace": branch.split("/")[0],
                "commits_ahead": ahead,
                "last_date": date,
                "last_subject": subject,
                "days_behind_head": _days_between(date, head_date) if head_date else None,
                "files_changed": len(paths),
                "insertions": ins,
                "deletions": dele,
                "areas": areas,
                "topics": classify_topics(" ".join(subjects), [], cfg),
                "missions": sorted({m for m in (parse_mission(s, cfg) for s in subjects) if m})[:10],
            }
        )
    return sorted(unmerged, key=lambda b: (b["last_date"], b["commits_ahead"]), reverse=True)


# --------------------------------------------------------------------------- #
# payload / artifacts
# --------------------------------------------------------------------------- #
def collect(
    repo: Path | str,
    ref: str = "origin/main",
    window: int = 15,
    network: bool = True,
    config_path: Path | str = DEFAULT_CONFIG,
) -> dict:
    repo = Path(repo)
    cfg = load_config(config_path)
    head = run_git(repo, "rev-parse", "HEAD").strip()
    branch = run_git(repo, "rev-parse", "--abbrev-ref", "HEAD").strip()
    head_date = run_git(repo, "show", "-s", "--format=%ad", "--date=short", head).strip()
    upstream = collect_upstream(repo, ref, head, network)
    target = ref if _resolve(repo, ref) else "HEAD"
    window_data = collect_window(repo, target, window, cfg)
    crosschecks = collect_crosschecks(repo, target, window, window_data)
    churn_span = collect_churn_span(repo, target, window_data)
    unmerged = collect_unmerged(repo, target, cfg, head_date) if upstream["verdict"] != "NO_UPSTREAM" else []

    totals = window_data["totals"]
    window_risks: list[str] = []
    if len(totals["authors"]) <= 1 and totals["commits"] >= cfg["risk_rules"]["bus_factor_min_commits"]:
        window_risks.append("single_author_window")
    if unmerged:
        window_risks.append("unmerged_upstream_branches")
    if upstream["verdict"] not in ("NO_UPSTREAM", "IN_SYNC"):
        window_risks.append("local_diverges_from_upstream")
    if not all(c["match"] for c in crosschecks):
        window_risks.append("crosscheck_mismatch")
    if totals["squash_landings"]:
        window_risks.append("squash_lane_present")
    if totals["empty_landings"]:
        window_risks.append("empty_landings_present")
    if totals["gate_filtered_files_reported"]:
        window_risks.append("gate_filtered_deliverables")

    return {
        "schema": SCHEMA,
        "repo_root": str(repo),
        "repo": {
            "head": head,
            "head_short": head[:9],
            "head_date": head_date,
            "branch": branch,
            "target_ref": target,
            "target_sha": _resolve(repo, target),
        },
        "params": {
            "ref": ref,
            "window": window,
            "network": network,
            "config": Path(config_path).name,
            "config_sha256": hashlib.sha256(Path(config_path).read_bytes()).hexdigest(),
            "lens": "first-parent",
        },
        "upstream": upstream,
        "window": {"commits": window_data["commits"], "totals": totals, "window_risks": window_risks},
        "churn_span": churn_span,
        "unmerged_branches": unmerged,
        "crosschecks": crosschecks,
        "evidence_grades": {
            "upstream_sync_verdict": "REAL" if network else "WITHHELD",
            "commit_metadata_churn": "REAL",
            "landing_mode_classification": "REAL",
            "theme_classification": "DERIVED",
            "gate_filtered_counts": "WITHHELD",
            "unmerged_branch_inventory": "REAL",
            "crosscheck_agreement": "REAL" if all(c["match"] for c in crosschecks) else "FAIL",
        },
    }


def canonical(payload: dict) -> str:
    stripped = {k: v for k, v in payload.items() if k not in VOLATILE_KEYS}
    return json.dumps(stripped, sort_keys=True, ensure_ascii=False, separators=(",", ":"))


def summary_of(payload: dict) -> dict:
    up = payload["upstream"]
    w = payload["window"]["totals"]
    return {
        "schema": payload["schema"],
        "head": payload["repo"]["head_short"],
        "target_ref": payload["repo"]["target_ref"],
        "upstream_verdict": up["verdict"],
        "ahead": up.get("ahead"),
        "behind": up.get("behind"),
        "live_probe": up.get("live_probe"),
        "fetched_ref_is_live": up.get("fetched_ref_is_live"),
        "window": payload["params"]["window"],
        "commits": w["commits"],
        "merge_commits": w["merge_commits"],
        "squash_landings": w["squash_landings"],
        "empty_landings": w["empty_landings"],
        "mission_commits": w["mission_commits"],
        "unique_missions": len(w["unique_missions"]),
        "insertions": w["insertions"],
        "deletions": w["deletions"],
        "files_touched_unique": w["files_touched_unique"],
        "gate_filtered_files_reported": w["gate_filtered_files_reported"],
        "themes": dict(sorted(w["topics"].items(), key=lambda kv: (-kv[1], kv[0]))),
        "areas": dict(sorted(w["areas"].items(), key=lambda kv: (-kv[1], kv[0]))),
        "landing_modes": dict(sorted(w["landing_modes"].items(), key=lambda kv: (-kv[1], kv[0]))),
        "authors": w["authors"],
        "risks": dict(sorted(w["risks"].items(), key=lambda kv: (-kv[1], kv[0]))),
        "window_risks": payload["window"]["window_risks"],
        "churn_amplification": payload["churn_span"]["amplification"],
        "unmerged_branches": len(payload["unmerged_branches"]),
        "crosschecks_ok": sum(1 for c in payload["crosschecks"] if c["match"]),
        "crosschecks_total": len(payload["crosschecks"]),
        "data_sha256": hashlib.sha256(canonical(payload).encode()).hexdigest(),
    }


def write_artifacts(payload: dict, outdir: Path | str, generated_at: str | None = None) -> dict:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    files = {}
    for name, body in (("patrol.json", payload), ("summary.json", summary_of(payload))):
        raw = json.dumps(body, indent=2, ensure_ascii=False, sort_keys=True).encode()
        (outdir / name).write_bytes(raw)
        files[name] = hashlib.sha256(raw).hexdigest()
    manifest = {
        "schema": SCHEMA,
        "files": files,
        "data_sha256": hashlib.sha256(canonical(payload).encode()).hexdigest(),
        "generated_at": generated_at or "",
        "generator": "tools/repo_patrol/patrol.py",
    }
    (outdir / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return manifest


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--repo", default=".")
    ap.add_argument("--ref", default="origin/main")
    ap.add_argument("--window", type=int, default=15)
    ap.add_argument("--outdir", default="artifacts/repo-patrol/upstream-diffset-v1")
    ap.add_argument("--config", default=str(DEFAULT_CONFIG))
    ap.add_argument("--offline", action="store_true", help="skip the live ls-remote probe")
    args = ap.parse_args(argv)
    payload = collect(args.repo, args.ref, args.window, network=not args.offline, config_path=args.config)
    manifest = write_artifacts(payload, args.outdir)
    print(json.dumps({**summary_of(payload), "manifest": manifest["files"]}, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
