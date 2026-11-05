#!/usr/bin/env python3
"""Figures for the wetted-surface uncertainty study (reads summary.json only)."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"figure.dpi": 130, "font.size": 9, "axes.grid": True, "grid.alpha": 0.3})


def _load(data_dir: Path) -> dict:
    return json.loads((data_dir / "summary.json").read_text())


def fig_mesh_gci(summary: dict, out: Path) -> Path:
    mesh = summary["mesh"]
    rows = [r for r in mesh["levels"] if "p_observed" in r]
    fig, ax = plt.subplots(figsize=(5.0, 3.6))
    h = [r["h_m"] for r in rows]
    err = [abs(r["error_relative"]) for r in rows]
    gci = [r["gci_relative_95"] for r in rows]
    ax.loglog(h, err, "o-", color="#1f4e79", label="|S_mesh − S_exact| / S_exact")
    ax.errorbar(
        h,
        err,
        yerr=[[min(gci[i], 0.9 * err[i]) for i in range(len(h))], [0.0] * len(h)],
        fmt="none",
        ecolor="#c00000",
        capsize=3,
        label="Roache GCI (Fs=1.25) as 95 % band",
    )
    p = mesh["finest_p_observed"]
    ref = err[0] * (h[0] / min(h)) ** p
    ax.loglog(
        [min(h), h[0]],
        [err[0] * (h[0] / h[0]) ** p, ref],
        "k--",
        lw=0.8,
        label=f"fitted order p = {p:.2f}",
    )
    ax.set_xlabel("facet size h = L/(nx-1)  [m]")
    ax.set_ylabel("relative area error")
    ax.set_title("Mesh convergence and GCI bracket (wetted area)")
    ax.legend(fontsize=7, loc="lower right")
    fig.tight_layout()
    path = out / "ws_uncertainty_mesh_gci.png"
    fig.savefig(path)
    plt.close(fig)
    return path


def fig_budget(summary: dict, out: Path) -> Path:
    budget = summary["budget"]
    comps = sorted(budget["components"], key=lambda c: abs(c["contribution_relative"]))
    labels = [c["name"] for c in comps]
    values = [abs(c["contribution_relative"]) for c in comps]
    fig, ax = plt.subplots(figsize=(5.4, 3.2))
    ax.barh(labels, values, color="#2e75b6")
    ax.set_xscale("log")
    ax.set_xlabel("standard contribution to u(dS/S), log scale")
    bias = abs(summary["closure"]["total_relative"])
    ax.axvline(bias, color="#c00000", ls="--", label=f"modelling bias |dS/S| = {bias*100:.2f} %")
    u95 = budget["expanded_95"] * 100.0
    ax.axvline(budget["expanded_95"], color="#548235", ls=":", label=f"U95 = {u95:.3f} %")
    ax.set_title("Uncertainty budget of the closure error")
    ax.legend(fontsize=7)
    fig.tight_layout()
    path = out / "ws_uncertainty_budget.png"
    fig.savefig(path)
    plt.close(fig)
    return path


def fig_errorbars(summary: dict, out: Path) -> Path:
    rows = summary["family_sweep"]
    cb = [r["cb"] for r in rows]
    total = [100.0 * r["total_relative"] for r in rows]
    u95 = [100.0 * r["u95_relative"] for r in rows]
    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    ax.errorbar(cb, total, yerr=u95, fmt="o-", color="#1f4e79", capsize=3,
                label="closure error dS/S with 95 % bar")
    ax.axhline(3.0, color="#c00000", ls="--", lw=0.8, label="3 % acceptance limit")
    ax.axhline(-3.0, color="#c00000", ls="--", lw=0.8)
    ax.axhline(0.0, color="k", lw=0.6)
    design = summary["design_cb"]
    ax.axvline(design, color="#7f7f7f", ls=":", lw=0.8, label=f"design Cb = {design:.4f}")
    ax.set_xlabel("block coefficient Cb of the analytic family")
    ax.set_ylabel("closure error [%]")
    ax.set_title("Closure error with error bars (tensorcad.hull family)")
    ax.legend(fontsize=7)
    fig.tight_layout()
    path = out / "ws_uncertainty_closure_bars.png"
    fig.savefig(path)
    plt.close(fig)
    return path


def fig_propagation(summary: dict, out: Path) -> Path:
    rows = summary["propagation"]
    v = [r["speed_kn"] for r in rows]
    y = [100.0 * r["drt_relative"] for r in rows]
    u = [100.0 * r["u95_drt_relative"] for r in rows]
    fig, ax = plt.subplots(figsize=(5.2, 3.4))
    ax.errorbar(v, y, yerr=u, fmt="o-", color="#7030a0", capsize=3, label="dRt/Rt with 95 % bar")
    ax.axhline(3.0, color="#c00000", ls="--", lw=0.8, label="3 % limit")
    ax.axhline(-3.0, color="#c00000", ls="--", lw=0.8)
    ax.set_xlabel("speed [kn]")
    ax.set_ylabel("propagated resistance error [%]")
    ax2 = ax.twinx()
    ax2.plot(v, [r["pi"] for r in rows], "s--", color="#548235", ms=3, lw=0.8, label="pi = Rf/Rt")
    ax2.set_ylabel("pi = Rf/Rt")
    ax2.grid(False)
    lines = ax.get_legend_handles_labels()
    lines2 = ax2.get_legend_handles_labels()
    ax.legend(lines[0] + lines2[0], lines[1] + lines2[1], fontsize=7, loc="lower right")
    ax.set_title("Wet-area error propagated to total resistance")
    fig.tight_layout()
    path = out / "ws_uncertainty_rt_band.png"
    fig.savefig(path)
    plt.close(fig)
    return path


def make_figures(data_dir: Path, figure_dir: Path) -> list[Path]:
    summary = _load(data_dir)
    figure_dir.mkdir(parents=True, exist_ok=True)
    paths = [
        fig_mesh_gci(summary, figure_dir),
        fig_budget(summary, figure_dir),
        fig_errorbars(summary, figure_dir),
        fig_propagation(summary, figure_dir),
    ]
    print(json.dumps({"figures": [str(p) for p in paths]}, indent=1))
    return paths


if __name__ == "__main__":

    repo = Path(__file__).resolve().parents[1]
    make_figures(repo / "data/wetted-uncertainty-2026-11-04", repo / "figures")
