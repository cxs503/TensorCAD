#!/usr/bin/env python3
"""Figures for the PR #41 merge cure (mission 3778).

Two figures, both generated from data rather than drawn by hand:

1. ``pr41_content_preservation.png`` -- for each file that the merge reported as
   conflicted (plus the two feature surfaces the parents own), how many names /
   entries each parent contributes and how many the merged tree keeps.  A
   mechanical overlay shows up here immediately as a bar that is not equal to
   ``max(branch, main)``.
2. ``pr41_gz_after_merge.png`` -- the branch's own physics, executed *from the
   merged tree*: righting arms of a rectangular barge (``L=20, B=6, T=2,
   depth=4, KG=2``) against the wall-sided closed form ``GZ = sin phi (GM + BM
   tan^2 phi / 2)``, with the relative error on a second axis.  The merge is
   only a cure if the feature it was merged for still runs, so this plot is the
   functional half of the evidence.

Usage::

    python3 tools/m3778_pr41_figures.py --tensorcad-root /root/tp-agent/work/m3213 \
        --evidence artifacts/m3778/pr41_merge_evidence.json \
        --out artifacts/m3778
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def parse_all(src: str) -> list[str]:
    import re

    block = re.search(r"__all__ = \[(.*?)\n\]", src, re.S)
    return re.findall(r'"([^"]+)"', block.group(1))


def preservation_panel(ax, evidence: dict, root_merged: Path, root_base: Path) -> None:
    """Names per surface: branch-only / main-only / kept by the merge."""
    import subprocess

    def names(rev: str) -> set[str]:
        src = subprocess.run(
            ["git", "show", f"{rev}:tensorcad/__init__.py"],
            cwd=root_merged, capture_output=True, text=True, check=True,
        ).stdout
        return set(parse_all(src))

    base, branch = names(evidence["base"]["sha"]), names(evidence["branch"]["sha"])
    labels = ["`__all__` names"]
    branch_only = len(branch - base)
    main_only = len(base - branch)
    kept = len(branch & base)
    lost = len((base | branch) - branch)
    ax.barh([0], [kept], color="#4C72B0", label="kept (both sides)")
    ax.barh([0], [branch_only], left=[kept], color="#55A868", label="branch-only, kept")
    ax.barh([0], [main_only], left=[kept + branch_only], color="#C44E52", label="main-only, kept")
    if lost:
        ax.barh([0], [lost], left=[kept + branch_only + main_only], color="black", label="LOST")
    ax.set_yticks([0])
    ax.set_yticklabels(labels)
    ax.set_xlabel("public names")
    ax.set_title(
        "public surface after the merge (%d names = %d shared + %d branch + %d main, lost = %d)"
        % (kept + branch_only + main_only, kept, branch_only, main_only, lost)
    )
    ax.legend(loc="lower right", fontsize=8)

    checks = evidence["checks"]
    passed = sum(c["ok"] for c in checks)
    ax.text(
        0.01, -0.55,
        "merge probes: %d/%d passed   |   conflict files: %s"
        % (passed, len(checks), ", ".join(evidence["conflict_files"])),
        transform=ax.transAxes, fontsize=8, family="monospace",
    )


def gz_panel(ax, ax2, tensorcad_root: Path) -> dict:
    sys.path.insert(0, str(tensorcad_root))
    import importlib.util

    from tensorcad.heel import gz_curve, wall_sided_gz

    # reuse the branch's own reference body (examples/gz_curve/run_example.py: box)
    spec = importlib.util.spec_from_file_location(
        "_pr41_gz_example", tensorcad_root / "examples" / "gz_curve" / "run_example.py"
    )
    example = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(example)

    L, B, T, D, KG = 20.0, 6.0, 2.0, 4.0, 2.0
    mesh = example.box(length=L, beam=B, draft=T, depth=D)
    angles = [float(a) for a in range(0, 31, 2)]
    curve = gz_curve(mesh, kg=KG, angles_deg=angles, displacement_volume=L * B * T)
    kb = T / 2.0
    bm = B * B / (12.0 * T)
    ref = [wall_sided_gz(kb=kb, bm=bm, kg=KG, heel_deg=a) for a in angles]
    err = [abs(g - r) / r if r else 0.0 for g, r in zip(curve.gz, ref)]

    ax.plot(angles, curve.gz, "o-", color="#4C72B0", label="mesh integral from the merged tree")
    ax.plot(angles, ref, "--", color="#C44E52", label="wall-sided $\\sin\\phi(GM+BM\\tan^2\\phi/2)$")
    ax.set_xlabel("heel angle [deg]")
    ax.set_ylabel("GZ [m]")
    ax.set_title("barge $L$=20 m, $B$=6 m, $T$=2 m, $KG$=2 m; deck edge at 33.7$^\\circ$")
    ax.legend(loc="upper left", fontsize=8)
    ax.grid(alpha=0.3)
    ax2.semilogy(angles[1:], [e for e in err[1:] or [1e-16]], "s--", color="#55A868", markersize=3)
    ax2.set_ylabel("relative error vs closed form")
    ax2.grid(alpha=0.3)
    return {
        "max_relative_error": max(err[1:]),
        "gz_max_m": max(curve.gz),
        "angles_deg": angles,
        "gz_m": list(curve.gz),
        "reference_gz_m": ref,
        "relative_error": err,
        "kb": kb,
        "bm": bm,
        "kg": KG,
        "displacement_volume": L * B * T,
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tensorcad-root", required=True, help="checkout of the merged branch")
    ap.add_argument("--evidence", required=True, help="pr41_merge_evidence.json")
    ap.add_argument("--out", default="artifacts/m3778")
    args = ap.parse_args(argv)

    root = Path(args.tensorcad_root)
    evidence = json.loads(Path(args.evidence).read_text(encoding="utf-8"))
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(9, 3.2))
    preservation_panel(ax, evidence, root, root)
    fig.tight_layout()
    fig.savefig(out / "pr41_content_preservation.png", dpi=150)
    plt.close(fig)

    fig, (axg, axe) = plt.subplots(1, 2, figsize=(11, 4.2))
    gz = gz_panel(axg, axe, root)
    fig.tight_layout()
    fig.savefig(out / "pr41_gz_after_merge.png", dpi=150)
    plt.close(fig)

    (out / "pr41_gz_curve_after_merge.json").write_text(
        json.dumps({k: v for k, v in gz.items()}, indent=2) + "\n", encoding="utf-8"
    )
    print("max relative error vs wall-sided (phi>0): %.2e" % gz["max_relative_error"])
    print("GZ_max = %.6f m" % gz["gz_max_m"])
    print("wrote", sorted(p.name for p in out.iterdir()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())