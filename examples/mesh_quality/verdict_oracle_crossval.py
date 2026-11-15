#!/usr/bin/env python3
"""Cross-validation of the mesh-quality metric module against the Verdict oracle.

Verdict is the element-quality library behind ``vtk.vtkMeshQuality``, i.e. behind
``meshmagick.print_quality()``.  ``tensorcad.meshquality.verdict_oracle`` exposes a build of
those exact sources (``tools/verdict_oracle/build.py``), so this script puts two
independently produced sets of numbers on meshes this repository generates.

How a pair is judged
--------------------
* ``identity``  same definition: the difference must be round-off,
  ``|a - b| <= TOL * max(1, |b|)``.
* ``factor``    same shape measure under two normalisations (tetrahedra: this tree
  normalises by the regular tetrahedron, Verdict by the orthogonal one), so the *ratio*
  is checked against the exact algebraic factor instead of the values;
* ``relative``  same physical quantity measured by two constructions (hexahedron volume:
  divergence theorem over six faces here, centre Jacobian in Verdict).  The relative
  difference is reported and the pair is informational, not a gate: on an affine cell the
  two constructions are identical, on a warped cell they are different questions.

Elements are excluded and counted when either side cannot score them: non-positive measure,
non-finite value, or Verdict's ``VERDICT_DBL_MAX`` sentinel (Verdict declines to score a
cell whose corner Jacobian is non-positive; this tree instead reports the negative value so
an inverted cell stays visible -- the exclusion counts are what makes that difference
auditable).

Outputs: ``artifacts/mesh_quality/verdict_crossval.{json,md,png}``.  Exit status 0 when
every gating pair passes, 1 otherwise.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import run_benchmark  # noqa: E402

from tensorcad.meshquality import verdict_oracle  # noqa: E402
from tensorcad.meshquality.metrics import (  # noqa: E402
    hex_metrics,
    quad_metrics,
    tet_metrics,
    tri_metrics,
)

ARTIFACT_DIR = REPO_ROOT / "artifacts" / "mesh_quality"

TOL = 1e-6
SENTINEL = 1e29  # Verdict's VERDICT_DBL_MAX, its "I will not score this" return

#: Meshes whose hexahedral cells have non-planar faces: the divergence-theorem volume and
#: Verdict's centre-Jacobian volume are then different questions, so their identity pair is
#: informational and the relative pair carries the number.
WARPED_HEX_MESHES = {"hull_hex"}

#: tetrahedron scaled Jacobian: this tree divides by the regular tetrahedron's value.
TET_SJ_FACTOR = 1.0 / math.sqrt(2.0)

#: (ours, theirs, mode, gating) per kind.
PAIRS: dict[str, tuple[tuple[str, str, str, bool], ...]] = {
    "tri": (
        ("area", "area", "identity", True),
        ("scaled_jacobian", "scaled_jacobian", "identity", True),
        ("min_angle", "min_angle", "identity", True),
        ("max_angle", "max_angle", "identity", True),
        ("edge_ratio", "edge_ratio", "identity", True),
    ),
    "quad": (
        ("area", "area", "identity", True),
        ("scaled_jacobian", "scaled_jacobian", "identity", True),
        ("min_angle", "min_angle", "identity", True),
        ("max_angle", "max_angle", "identity", True),
        ("edge_ratio", "edge_ratio", "identity", True),
    ),
    "tet": (
        ("volume", "volume", "identity", True),
        ("edge_ratio", "edge_ratio", "identity", True),
        ("scaled_jacobian", "scaled_jacobian", "factor", True),
    ),
    "hex": (
        ("volume", "volume", "identity", True),
        ("scaled_jacobian", "scaled_jacobian", "identity", True),
        ("edge_ratio", "edge_ratio", "identity", True),
        ("volume", "volume", "relative", False),
    ),
}

#: Metrics both libraries compute but define differently -- kept visible, never compared.
DIFFERENT_BY_DEFINITION = {
    "tet.min_angle": "this tree reports tetrahedron dihedral angles, Verdict reports face angles",
    "tet.min_sicn": "Verdict's condition number is ||J||_F||J^-1||_F/d; min_sicn is "
                    "d|det J|^(2/d)/||J||_F^2 -- they coincide only for isotropic cells",
    "hex.min_angle": "Verdict has no hex face-angle metric; this tree's hex angles are face angles",
    "skewness": "this tree uses Fluent equiangle skew on the cell's own angle set; Verdict's "
                "*_equiangle_skew applies its per-kind reference angles",
    "hex.volume(warped)": "divergence theorem over the six (non-planar) faces vs the trilinear "
                          "Jacobian at the cell centre: identical for affine cells, different "
                          "questions otherwise -- reported as the 'relative' pair",
}


def tensorcad_quality(kind: str, points, cells):
    return {
        "tri": tri_metrics,
        "quad": quad_metrics,
        "tet": tet_metrics,
        "hex": hex_metrics,
    }[kind](points, cells)


def judge(ours: list[float], theirs: list[float], measure: list[float], mode: str) -> dict:
    diffs: list[float] = []
    ratios: list[float] = []
    excluded = 0
    for a, b, m in zip(ours, theirs, measure):
        if not (math.isfinite(a) and math.isfinite(b)) or m <= 0.0 or abs(b) >= SENTINEL:
            excluded += 1
            continue
        if mode == "factor":
            if a == 0.0 or b == 0.0:
                excluded += 1
                continue
            ratios.append(a / b)
            continue
        if mode == "relative":
            if b == 0.0:
                excluded += 1
                continue
            diffs.append(abs(a - b) / abs(b))
            continue
        diffs.append(abs(a - b) / max(1.0, abs(b)))
    if mode == "factor":
        worst = max((abs(r - TET_SJ_FACTOR) / TET_SJ_FACTOR for r in ratios), default=float("nan"))
        ratio_range = (min(ratios), max(ratios)) if ratios else (float("nan"), float("nan"))
    else:
        worst = max(diffs, default=float("nan"))
        ratio_range = (float("nan"), float("nan"))
    return {
        "mode": mode,
        "worst": worst,
        "median": (sorted(diffs)[len(diffs) // 2] if diffs else float("nan")),
        "compared": len(ratios) if mode == "factor" else len(diffs),
        "excluded": excluded,
        "ratio_range": ratio_range,
    }


def main() -> int:
    if not verdict_oracle.available():
        print("no Verdict oracle built: run python tools/verdict_oracle/build.py", file=sys.stderr)
        return 1
    meshes = run_benchmark.ideal_elements() | run_benchmark.build_meshes()
    results: dict = {
        "schema": "tensorcad.meshquality.verdict_crossval/2",
        "tolerance": TOL,
        "verdict_sentinel": SENTINEL,
        "tet_scaled_jacobian_factor": TET_SJ_FACTOR,
        "oracle": verdict_oracle.provenance(),
        "oracle_info": verdict_oracle.info(),
        "different_by_definition": DIFFERENT_BY_DEFINITION,
        "meshes": {},
    }
    gate_worst = 0.0
    gate_where = ""
    rows: list[str] = [
        "| mesh.kind | pair (this tree <-> Verdict) | mode | worst | compared | excluded | verdict |",
        "|---|---|---|---:|---:|---:|---|",
    ]
    for mesh_name, (points, elements) in meshes.items():
        for kind, cells in elements.items():
            quality = tensorcad_quality(kind, points, cells)
            theirs_all = verdict_oracle.metrics(kind, points, cells)
            measure = list(quality.metrics["volume" if kind in ("tet", "hex") else "area"])
            record: dict = {"kind": kind, "cells": len(cells), "points": len(points), "pairs": {}}
            for ours_name, theirs_name, mode, gating in PAIRS[kind]:
                verdict_values = theirs_all[theirs_name]
                if ours_name == "volume" and mode == "relative":
                    # the affine-volume identity is already the gating pair above; here we
                    # only quote the residual on the warped cells.
                    pass
                verdict_row = judge(
                    list(quality.metrics[ours_name]), verdict_values, measure, mode
                )
                verdict_row["theirs"] = theirs_name
                if (kind == "hex" and ours_name == "volume" and mode == "identity"
                        and mesh_name in WARPED_HEX_MESHES):
                    gating = False
                verdict_row["gating"] = gating
                ok = gating and (verdict_row["worst"] <= TOL)
                verdict_row["pass"] = ok
                label = f"{ours_name}<->{theirs_name}"
                record["pairs"][f"{label}:{mode}"] = verdict_row
                rows.append(
                    f"| `{mesh_name}.{kind}` | `{label}` | {mode} | "
                    f"{verdict_row['worst']:.3e} | {verdict_row['compared']} | "
                    f"{verdict_row['excluded']} | {'PASS' if gating and ok else ('INFO' if not gating else 'FAIL')} |"
                )
                if gating and verdict_row["worst"] > gate_worst:
                    gate_worst = verdict_row["worst"]
                    gate_where = f"{mesh_name}.{kind}.{label}({mode})"
            results["meshes"][f"{mesh_name}.{kind}"] = record

    results["worst_gating"] = {"where": gate_where, "worst": gate_worst}
    results["pass"] = gate_worst <= TOL
    ARTIFACT_DIR.mkdir(parents=True, exist_ok=True)
    (ARTIFACT_DIR / "verdict_crossval.json").write_text(json.dumps(results, indent=2) + "\n")
    (ARTIFACT_DIR / "verdict_crossval.md").write_text("\n".join(rows) + "\n")
    try:
        _plot(results)
    except Exception as exc:  # plotting is a convenience, not the result
        print(f"plot skipped: {exc}", file=sys.stderr)
    print(f"worst gating pair: {gate_where} = {gate_worst:.3e} (tolerance {TOL:.0e})")
    print(f"artifact: {ARTIFACT_DIR / 'verdict_crossval.json'}")
    return 0 if results["pass"] else 1


def _plot(results: dict) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels: list[str] = []
    values: list[float] = []
    colors: list[str] = []
    for mesh_key, record in results["meshes"].items():
        for pair, row in record["pairs"].items():
            labels.append(f"{mesh_key}\n{pair.split(':')[0]}")
            values.append(max(row["worst"], 1e-18))
            colors.append("#35618f" if row["gating"] else "#b07d3a")
    fig, ax = plt.subplots(figsize=(max(10.0, 0.36 * len(labels)), 4.8))
    ax.bar(range(len(values)), values, color=colors)
    ax.set_yscale("log")
    ax.axhline(TOL, color="#c0392b", linestyle="--", linewidth=1.0, label=f"tolerance {TOL:.0e}")
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=90, fontsize=5.5)
    ax.set_ylabel("worst scaled difference  (identity: |a-b|/max(1,|b|))")
    ax.set_title("Mesh-quality metrics vs Verdict (the library behind vtkMeshQuality/meshmagick)")
    ax.legend(loc="upper left", fontsize=8)
    fig.tight_layout()
    fig.savefig(ARTIFACT_DIR / "verdict_crossval.png", dpi=170)


if __name__ == "__main__":
    raise SystemExit(main())