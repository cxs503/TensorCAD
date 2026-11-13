#!/usr/bin/env python3
"""Figures for docs/baseline_fleet.md, drawn from the gate's own JSON reports."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

OUT = Path("artifacts/baseline_fleet/figs")
ROOT = Path("artifacts/baseline_fleet")


def main() -> int:
    audit = json.loads((ROOT / "placeholder_audit.json").read_text())
    control = json.loads((ROOT / "control_real_mesh.json").read_text())
    published = json.loads((ROOT / "published.json").read_text())
    ships = [(s["ship"], s["mesh"] or {}, s["grade"]) for s in audit["ships"]]
    ships.sort(key=lambda item: item[1].get("triangles", 0))

    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    names = [s[0] for s in ships]
    tris = [max(s[1].get("triangles", 0), 1) for s in ships]
    colours = ["#c0392b" if s[2] != "VERIFIED" and s[2] != "WITHHELD" else "#2c7fb8" for s in ships]
    ax.barh(names, tris, color=colours)
    ax.axvline(audit["min_hull_triangles"], ls="--", c="k",
               label=f"gate = {audit['min_hull_triangles']} facets")
    for i, (name, mesh, grade) in enumerate(ships):
        ax.text(tris[i] * 1.15, i, f"{mesh.get('triangles', 0)} ({grade})", va="center", fontsize=8)
    ax.set_xscale("log")
    ax.set_xlabel("triangles in *_Hull.stl (log scale)")
    ax.set_title("Fleet facet gate: the six placeholders sit 2+ decades under real tessellations")
    ax.legend(loc="lower right", fontsize=8)
    ax.grid(axis="x", alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "facet_gate.png", dpi=170)

    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    common = [s for s in ships if s[0] in published]
    labels = [s[0] for s in common]
    measured = [abs(s[1]["volume_m3"]) for s in common]
    reference = [published[s[0]] for s in common]
    x = range(len(labels))
    ax.bar([i - 0.2 for i in x], measured, width=0.4, label="placeholder |V| (m³)")
    ax.bar([i + 0.2 for i in x], reference, width=0.4, label="published displacement (m³)")
    for i, (m, r) in enumerate(zip(measured, reference)):
        ax.text(i, max(m, r) * 1.03, f"error {(abs(m - r) / r) * 100:.0f} %", ha="center", fontsize=8)
    ax.set_xticks(list(x), labels)
    ax.set_ylabel("volume (m³)")
    ax.set_title("Volume gate: placeholder boxes vs the publishers' own displacement figures")
    ax.legend(fontsize=8)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "volume_gate.png", dpi=170)
    print("wrote", OUT / "facet_gate.png", OUT / "volume_gate.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
