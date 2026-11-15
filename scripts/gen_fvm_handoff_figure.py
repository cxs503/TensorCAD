"""Figure for docs/fvm_handoff.md: patch-test exactness, frame price, harmonic error."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tensorcad import fvm_handoff as fh  # noqa: E402

OUT = Path("artifacts/fvm_handoff")
payload = json.loads((OUT / "results.json").read_text())
ran = [c for c in payload["cases"] if c["status"] == "ran"]


def pick(name: str, field_prefix: str) -> dict:
    """Same case name carries several checks; the field name says which one."""
    match = [c for c in ran if c["case"] == name and c["metrics"].get("field", "").startswith(field_prefix)]
    if not match:
        raise SystemExit(f"no {field_prefix!r} result for {name!r}")
    return match[0]


barge_small = pick("barge_2x0.4x0.3_9x5x2", "linear")
barge_big = pick("barge_2x0.4x0.3_17x9x3", "linear")
harm_small = pick("barge_2x0.4x0.3_9x5x2", "sin(")
harm_big = pick("barge_2x0.4x0.3_17x9x3", "sin(")
hull = pick("hull_10x1x0.6_13x7x2", "linear")
frame = [c for c in ran if "yup_relative_linf" in c["metrics"]][0]

fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.4))

# (a) patch-test exactness vs the cost of the frame mix-up, log scale
labels = [
    "patch\nbarge 768",
    "patch\nbarge 4608",
    "patch\nhull 1728",
    "same table\nread y-up",
]
values = [
    barge_small["metrics"]["relative_linf"],
    barge_big["metrics"]["relative_linf"],
    hull["metrics"]["relative_linf"],
    frame["metrics"]["yup_relative_linf"],
]
colours = ["#2a7f62", "#2a7f62", "#2a7f62", "#b23b3b"]
bars = axes[0].bar(labels, values, color=colours)
axes[0].set_yscale("log")
axes[0].set_ylabel("relative $L_\\infty$ error of a linear field [-]")
axes[0].set_title("(a) linear-field patch test vs y-up frame mistake")
axes[0].axhline(1e-9, ls="--", c="k", lw=0.8, label="patch bound 1e-9")
for bar, value in zip(bars, values):
    axes[0].annotate(f"{value:.1e}", (bar.get_x() + bar.get_width() / 2, value),
                     ha="center", va="bottom", fontsize=9)
axes[0].set_ylim(1e-15, 1.0)
axes[0].legend(fontsize=8)

# (b) harmonic error under refinement
cells = [harm_small["metrics"]["cells"], harm_big["metrics"]["cells"]]
linf = [harm_small["metrics"]["relative_linf"], harm_big["metrics"]["relative_linf"]]
l2 = [harm_small["metrics"]["relative_l2"], harm_big["metrics"]["relative_l2"]]
axes[1].plot(cells, linf, "o-", label="$L_\\infty$")
axes[1].plot(cells, l2, "s-", label="$L_2$")
axes[1].axhline(0.03, ls="--", c="k", lw=0.8, label="3 % bar")
axes[1].set_xscale("log")
axes[1].set_xlabel("tetrahedra [-]")
axes[1].set_ylabel("relative error [-]")
axes[1].set_title("(b) harmonic field vs separation of variables")
for x, y in zip(cells, l2):
    axes[1].annotate(f"{y*100:.2f} %", (x, y), textcoords="offset points", xytext=(0, -14), ha="center", fontsize=8)
for x, y in zip(cells, linf):
    axes[1].annotate(f"{y*100:.2f} %", (x, y), textcoords="offset points", xytext=(0, 6), ha="center", fontsize=8)
axes[1].legend(fontsize=8)
axes[1].grid(alpha=0.3)

# (c) the solution itself: computed vs exact along the box diagonal
bundle = OUT / "barge_17x9x3"
solved = fh._solve(fh._open_bundle(bundle),
                   lambda p: fh.exact_harmonic(p, length=2.0, beam=0.4, draft=0.3))
centres = solved["cell_centers"]
t = centres[:, 0] / 2.0
order = np.argsort(t)
sample = int(len(t) / 6)
pick = order[::sample]
axes[2].plot(t[pick], solved["reference_values"][pick], "k.", ms=6, label="exact")
axes[2].plot(t[pick], solved["cell_values"][pick], "ro", ms=4, mfc="none", label="FVM (cell centre)")
axes[2].set_xlabel("normalised station $x/L$")
axes[2].set_ylabel("$\\phi$ at cell centres [-]")
axes[2].set_title("(c) harmonic field, 4608 tetrahedra")
axes[2].legend(fontsize=8)
axes[2].grid(alpha=0.3)

fig.tight_layout()
target = OUT / "figure.png"
fig.savefig(target, dpi=140)
print("wrote", target, "cells(b)", cells, "linf", linf, "l2", l2)
