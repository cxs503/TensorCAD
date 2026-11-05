#!/usr/bin/env python3
"""Figures for the PR #122 CI-red cure (mission m4390).

Reads ``ci_snapshot.json`` (verbatim GitHub Actions + ruff JSON, captured
programmatically -- no hand-typed numbers) and writes two PNGs:

  fig_lint_before_after.png -- ruff diagnostics per rule, before / after the cure
  fig_ci_matrix.png        -- per-Python-version CI job verdict, before / after

Usage: python3 generate_figures.py [--snapshot ci_snapshot.json]
"""
from __future__ import annotations

import argparse
import collections
import json
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = pathlib.Path(__file__).resolve().parent
COLORS = {"failure": "#c0392b", "cancelled": "#95a5a6", "success": "#27ae60",
          "queued": "#f39c12", "in_progress": "#f39c12", None: "#ecf0f1"}


def lint_figure(snap: dict, out: pathlib.Path) -> None:
    before = collections.Counter(d["code"] for d in snap["lint"]["before"])
    after = collections.Counter(d["code"] for d in snap["lint"]["after"])
    codes = sorted(before)
    files = sorted({d["file"] for d in snap["lint"]["before"]})
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))

    x = range(len(codes))
    ax1.bar([i - 0.2 for i in x], [before[c] for c in codes], 0.4, label="before (b061e9a)", color="#c0392b")
    ax1.bar([i + 0.2 for i in x], [after[c] for c in codes], 0.4, label="after (fce741e)", color="#27ae60")
    ax1.set_xticks(list(x))
    ax1.set_xticklabels(codes)
    ax1.set_ylabel("ruff diagnostics")
    ax1.set_title("ruff check . : diagnostics by rule")
    for i, c in enumerate(codes):
        ax1.text(i - 0.2, before[c] + 0.05, str(before[c]), ha="center", fontsize=9)
        ax1.text(i + 0.2, 0.05, str(after[c]), ha="center", fontsize=9)
    ax1.set_ylim(0, max(before.values()) + 0.6)
    ax1.legend()
    ax1.grid(axis="y", alpha=0.3)

    ax2.axis("off")
    rows = [(f, sum(1 for d in snap["lint"]["before"] if d["file"] == f),
             sum(1 for d in snap["lint"]["after"] if d["file"] == f)) for f in files]
    rows.append(("<total>", len(snap["lint"]["before"]), len(snap["lint"]["after"])))
    tbl = ax2.table(cellText=[[f, b, a] for f, b, a in rows],
                    colLabels=["file", "before", "after"], loc="center", cellLoc="left")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9)
    tbl.scale(1, 1.6)
    ax2.set_title("files carrying the diagnostics")
    fig.suptitle("PR #122 lint gate: every diagnostic lived in .m3753_scratch/", fontsize=12)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def ci_figure(snap: dict, out: pathlib.Path) -> None:
    runs = snap["runs"]
    versions = sorted({j["name"].replace("test (", "").rstrip(")") for j in runs["before"]["jobs"]})
    fig, ax = plt.subplots(figsize=(6.6, 3.6))
    for col, label in enumerate(("before", "after")):
        jobs = {j["name"].replace("test (", "").rstrip(")"): j for j in runs[label]["jobs"]}
        for row, ver in enumerate(versions):
            j = jobs.get(ver, {})
            state = j.get("conclusion") or j.get("status")
            ax.add_patch(plt.Rectangle((col, row), 0.9, 0.9,
                                       color=COLORS.get(state, "#ecf0f1"), ec="white"))
            step = j.get("failed_step")
            txt = f"{state}" + (f"\n@ {step}" if step else "")
            ax.text(col + 0.45, row + 0.45, txt, ha="center", va="center",
                    fontsize=8, color="white" if state in COLORS else "black")
    ax.set_xlim(0, 2)
    ax.set_ylim(0, len(versions))
    ax.set_xticks([0.45, 1.45])
    ax.set_xticklabels([f"before\n{runs['before']['head_sha'][:7]} (b061e9a)",
                        f"after\n{runs['after']['head_sha'][:7]} (fce741e)"])
    ax.set_yticks([i + 0.45 for i in range(len(versions))])
    ax.set_yticklabels([f"py {v}" for v in versions])
    ax.set_title("tests workflow matrix, PR #122")
    ax.grid(False)
    for s in ax.spines.values():
        s.set_visible(False)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--snapshot", default=str(HERE / "ci_snapshot.json"))
    args = ap.parse_args()
    snap = json.loads(pathlib.Path(args.snapshot).read_text())
    lint_figure(snap, HERE / "fig_lint_before_after.png")
    ci_figure(snap, HERE / "fig_ci_matrix.png")
    print("wrote", HERE / "fig_lint_before_after.png", HERE / "fig_ci_matrix.png")


if __name__ == "__main__":
    main()
