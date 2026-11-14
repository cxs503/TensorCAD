"""Board cross-check audit of mission m20416 (anchor XVAL-GEOM-VOLCROSS).

Read-only with respect to the audited worktree: every number below is either
read from the audited mission's own artifact, queried from the platform
referee table, or measured by compiling / importing a *copy* of the delivered
tree under /tmp.  Nothing here repairs the audited tree.

Outputs (in this directory):
  evidence_m20416.json   machine-readable evidence bundle
  fig_m20416.png         deviation-per-route bars + loft convergence
"""

from __future__ import annotations

import json
import math
import shutil
import subprocess
import sys
from pathlib import Path

WORK = Path("/root/tp-agent/work/m20416")           # audited mission worktree
HERE = Path(__file__).resolve().parent
PROBE = Path("/tmp/m20464_probe")                    # disposable copy
ANCHOR = "XVAL-GEOM-VOLCROSS"
MISSION = 20416
BAND = 0.003          # the module's own band
OWNER_TOL = 0.15      # owner acceptance standard


def psql(sql: str) -> str:
    out = subprocess.run(
        ["psql", "-U", "postgres", "-d", "tensorplatform", "-t", "-A", "-F", "\t", "-c", sql],
        capture_output=True, text=True, check=False,
    )
    return out.stdout.strip()


def referee_rows():
    total = psql("SELECT count(*) FROM tensor_benchmark_cases")
    rows = psql(
        "SELECT case_name, metric, tolerance FROM tensor_benchmark_cases "
        f"WHERE anchor_slug='{ANCHOR}'"
    )
    return {
        "table_rows_total": int(total or 0),
        "rows_for_anchor": 0 if not rows else len(rows.splitlines()),
        "rows": [] if not rows else rows.splitlines(),
        "referee_present": bool(rows),
    }


def mission_record():
    raw = psql(f"SELECT result::text FROM tensor_dev_missions WHERE id={MISSION}")
    d = json.loads(raw)
    ev = d.get("evidence", {})
    return {
        "declared_status": ev.get("finished_status"),
        "fail_stage": d.get("fail_stage"),
        "accepted": d.get("accepted"),
        "agent_self_reported_checks": d.get("checks"),
        "head_commit": ev.get("head_commit"),
        "branch": ev.get("branch"),
        "server_quality_audit": ev.get("server_quality_audit"),
        "tool_calls": d.get("tool_calls"),
        "on_disk_changes": d.get("on_disk_changes"),
    }


def deliverable_state():
    st = {}
    for rel in ("tensorcad/volume_xval.py", "tests/test_volume_xval.py",
                "tools/volume_xval_study.py", "docs/volume_cross_validation.md",
                "artifacts/volume_xval_study.json", "docs/volume_cross_validation.png"):
        p = WORK / rel
        st[rel] = {
            "exists": p.exists(),
            "bytes": p.stat().st_size if p.exists() else 0,
        }
    # tracked vs untracked: does the repository contain the deliverable?
    git = subprocess.run(["git", "status", "--porcelain"], cwd=WORK,
                         capture_output=True, text=True, check=False).stdout
    st["git_status_porcelain"] = git.strip().splitlines()
    st["untracked_deliverables"] = [l[3:] for l in git.splitlines() if l.startswith("??")]
    # the repository HEAD of the audited worktree does not contain the module
    have = subprocess.run(["git", "cat-file", "-e", "HEAD:tensorcad/volume_xval.py"],
                          cwd=WORK, capture_output=True, text=True, check=False)
    st["in_repo_head"] = have.returncode == 0
    return st


def compile_state():
    r = subprocess.run([sys.executable, "-m", "py_compile", "tensorcad/volume_xval.py"],
                       cwd=WORK, capture_output=True, text=True, check=False)
    entry = subprocess.run([sys.executable, "-m", "py_compile", "tools/volume_xval_study.py"],
                           cwd=WORK, capture_output=True, text=True, check=False)
    test = subprocess.run([sys.executable, "-m", "py_compile", "tests/test_volume_xval.py"],
                          cwd=WORK, capture_output=True, text=True, check=False)
    # as-delivered test run (collection only: fast, proves importability)
    col = subprocess.run([sys.executable, "-m", "pytest", "tests/test_volume_xval.py",
                          "-q", "--collect-only"], cwd=WORK, capture_output=True,
                         text=True, check=False,
                         env={"PYTHONPATH": ".", "PATH": "/usr/local/bin:/usr/bin:/bin"})
    return {
        "impl_compiles": r.returncode == 0,
        "impl_error": (r.stderr or "").strip().splitlines()[-1:] or [],
        "entry_compiles": entry.returncode == 0,
        "test_file_compiles": test.returncode == 0,
        "as_delivered_pytest_collect_rc": col.returncode,
        "as_delivered_pytest_collect_tail": (col.stdout + col.stderr).strip().splitlines()[-3:],
    }


def probe_copy():
    """Copy the delivered tree to /tmp and apply the *minimal mechanical* repair
    needed to make the module importable, then measure whether the functional
    substance is intact.  The audited tree itself is not touched."""
    if PROBE.exists():
        shutil.rmtree(PROBE)
    PROBE.mkdir(parents=True)
    for sub in ("tensorcad", "tests", "tools", "docs", "artifacts"):
        shutil.copytree(WORK / sub, PROBE / sub)
    src_path = PROBE / "tensorcad/volume_xval.py"
    src = src_path.read_text()
    start = src.index("def _append_sweeps(")
    head, tail = src[:start], src[start:]
    if "    station\n    if not sweeps:" in tail:
        tail = tail.replace(
            '    station\n    if not sweeps:\n        return\n    station = sweeps["stations"]\n        level',
            '    station = sweeps["stations"]\n    level',
        )
        body_at = tail.index('    station = sweeps["stations"]')
        head2, body = tail[:body_at], tail[body_at:]
        body = "\n".join(
            (l[4:] if l.startswith("        ") and not l.startswith("    station = sweeps") else l)
            for l in body.split("\n")
        )
        tail = head2 + body
        src_path.write_text(head + tail)
        repaired = True
    else:
        repaired = False
    comp = subprocess.run([sys.executable, "-m", "py_compile", "tensorcad/volume_xval.py"],
                          cwd=PROBE, capture_output=True, text=True, check=False)
    tests = subprocess.run([sys.executable, "-m", "pytest", "tests/test_volume_xval.py", "-q"],
                           cwd=PROBE, capture_output=True, text=True, check=False,
                           env={"PYTHONPATH": ".", "PATH": "/usr/local/bin:/usr/bin:/bin"},
                           timeout=900)
    last = (tests.stdout or "").strip().splitlines()[-1:] or []
    # the study tool rewrites artifacts/volume_xval_study.json in the probe copy
    study = subprocess.run([sys.executable, "tools/volume_xval_study.py"], cwd=PROBE,
                           capture_output=True, text=True, check=False,
                           env={"PYTHONPATH": ".", "PATH": "/usr/local/bin:/usr/bin:/bin"},
                           timeout=1800)
    return {
        "repair_applied": repaired,
        "repair_kind": "remove 1 stray token + 1 duplicated guard, dedent 20-line block by 4 "
                       "(mechanical; no logic added or removed)",
        "compiles_after_repair": comp.returncode == 0,
        "pytest_after_repair": last,
        "pytest_rc": tests.returncode,
        "study_rc": study.returncode,
        "study_tail": (study.stdout or "").strip().splitlines()[-1:],
    }


def reproduce_artifact():
    """Compare the mission artifact against the one the probe copy regenerated."""
    a = json.loads((WORK / "artifacts/volume_xval_study.json").read_text())
    b_path = PROBE / "artifacts/volume_xval_study.json"
    if not b_path.exists():
        return {"compared": False}
    b = json.loads(b_path.read_text())

    def flat(d, p=""):
        if isinstance(d, dict):
            for k, v in d.items():
                yield from flat(v, f"{p}/{k}")
        elif isinstance(d, list):
            for i, v in enumerate(d):
                yield from flat(v, f"{p}[{i}]")
        else:
            yield p, d

    fa, fb = dict(flat(a)), dict(flat(b))
    numeric = string = 0
    diffs = []
    for k, x in fa.items():
        if k not in fb:
            continue
        y = fb[k]
        if isinstance(x, float) and isinstance(y, float):
            if math.isnan(x) and math.isnan(y):
                continue
            numeric += 1
            if abs(x - y) > 1e-12 * max(1.0, abs(x)):
                diffs.append({"leaf": k, "artifact": x, "reproduced": y})
        else:
            string += 1
            if x != y:
                diffs.append({"leaf": k, "artifact": str(x)[:90], "reproduced": str(y)[:90]})
    return {
        "compared": True,
        "numeric_leaves_compared": numeric,
        "string_leaves_compared": string,
        "numeric_differences": sum(1 for d in diffs if isinstance(d["artifact"], (int, float))),
        "differences": diffs[:10],
    }


def numbers_from_artifact():
    a = json.loads((WORK / "artifacts/volume_xval_study.json").read_text())
    out = {"band": a["checks"]["band"], "ships": {}, "sweeps": {}, "loft": {}}
    for c in a["checks"]["checks"]:
        ref = c["reference"]
        rows = {"reference": {"engine": ref["engine"], "volume_m3": ref["volume_m3"],
                              "evidence": ref["evidence"]}, "routes": []}
        Vref = ref["volume_m3"]
        for r in c["routes"]:
            v = r["volume_m3"]
            rows["routes"].append({
                "engine": r["engine"], "role": r["role"], "evidence": r["evidence"],
                "volume_m3": None if isinstance(v, float) and math.isnan(v) else v,
                "rel_dev_vs_ref": None if (v is None or (isinstance(v, float) and math.isnan(v))
                                           or not Vref) else (v - Vref) / Vref,
            })
        rows["engine_spread"] = c.get("engine_spread")
        rows["worst_deviation"] = c.get("worst_deviation")
        rows["within_owner_15pct"] = None if not Vref else abs(c.get("worst_deviation") or 0) < OWNER_TOL
        out["ships"][c["ship"]] = rows
    # per-quantity error verdict, computed here, never transcribed
    for name, s in out["ships"].items():
        ref = s["reference"]["volume_m3"]
        for r in s["routes"]:
            if r["volume_m3"] is not None and ref and not math.isnan(ref):
                r["abs_err_pct"] = abs(r["volume_m3"] - ref) / abs(ref) * 100.0
    out["attribution"] = {k: {"mechanism": v.get("mechanism"),
                              "engine_spread": v.get("engine_spread"),
                              "worst_deviation": v.get("worst_deviation")}
                          for k, v in a["attribution"].items()}
    out["loft_convergence"] = a.get("loft_convergence")
    out["sweeps"] = a.get("sweeps")
    return out


def integration_callers():
    import re
    hits = []
    for p in WORK.rglob("*.py"):
        rel = p.relative_to(WORK).as_posix()
        if rel in ("tensorcad/volume_xval.py", "tests/test_volume_xval.py",
                   "tools/volume_xval_study.py"):
            continue
        try:
            txt = p.read_text(errors="ignore")
        except OSError:
            continue
        if re.search(r"volume_xval|cross_validate_(hull|standard_hulls)", txt):
            hits.append(rel)
    api = subprocess.run(["grep", "-rn", "volume_xval", "--include=*.py",
                          "tensorcad/api/"], cwd=WORK, capture_output=True, text=True,
                         check=False).stdout.strip().splitlines()
    return {"repo_callers_outside_deliverables": hits,
            "rest_api_callers": api,
            "wired_into_platform": bool(hits or api)}


def figure(ev):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ships = list(ev["numbers"]["ships"])
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
    ax = axes[0]
    width = 0.2
    labels = ["self surface sum", "trimesh STL", "ray stabbing"]
    for i, ship in enumerate(ships):
        routes = [r for r in ev["numbers"]["ships"][ship]["routes"] if r.get("abs_err_pct") is not None]
        vals = [r["abs_err_pct"] for r in routes][:3]
        xs = [i + (j - 1) * width for j in range(len(vals))]
        ax.bar(xs, vals, width=width, label=ship)
        for x, v in zip(xs, vals):
            ax.text(x, v * 1.03, f"{v:.3f}", ha="center", va="bottom", fontsize=7)
    ax.axhline(BAND * 100, color="k", ls="--", lw=0.9, label="module band 0.3 %")
    ax.set_xticks(range(len(ships)))
    ax.set_xticklabels(ships)
    ax.set_ylabel("|V_route − V_reference| / V_reference  [%]")
    ax.set_title("m20416 delivered artifact: per-route deviation (log scale)")
    ax.set_yscale("log")
    ax.legend(fontsize=7)

    ax = axes[1]
    loft = ev["numbers"].get("loft_convergence") or {}
    plotted = False
    for ship, blk in loft.items():
        if not isinstance(blk, dict):
            continue
        rows = blk.get("rows") or []
        if not rows:
            continue
        try:
            h = [min(r.get("nx", 0), 1) and 1.0 / r["nx"] for r in rows]
            e = [abs(r["relative_error"]) * 100 for r in rows]
        except Exception:
            continue
        ax.loglog(h, e, "o-", label=f"{ship} (obs. order {blk.get('order_nx', float('nan')):.2f})")
        plotted = True
    if plotted:
        ax.set_xlabel("1/nx  (station spacing / L)")
        ax.set_ylabel("|error| vs reference  [%]")
        ax.set_title("loft convergence harness (needs a referee/anchor)")
        ax.legend(fontsize=7)
        ax.grid(True, which="both", alpha=0.3)
    fig.suptitle("XVAL-GEOM-VOLCROSS cross-check (m20464 audit of m20416)", fontsize=10)
    fig.tight_layout()
    fig.savefig(HERE / "fig_m20416.png", dpi=150)


def main():
    ev = {
        "mission": MISSION,
        "anchor": ANCHOR,
        "audited_worktree": str(WORK),
        "referee": referee_rows(),
        "mission_record": mission_record(),
        "deliverables": deliverable_state(),
        "compile": compile_state(),
        "numbers": numbers_from_artifact(),
        "integration": integration_callers(),
    }
    ev["probe"] = probe_copy()
    ev["reproduction"] = reproduce_artifact()
    (HERE / "evidence_m20416.json").write_text(json.dumps(ev, indent=1, default=str))
    try:
        figure(ev)
    except Exception as exc:                                    # pragma: no cover
        print("figure failed:", exc)
    print(json.dumps({k: v for k, v in ev.items() if k != "numbers"}, indent=1, default=str)[:4000])


if __name__ == "__main__":
    main()