#!/usr/bin/env python
"""Generate the 6-format × 3-geometry round-trip matrix artifacts.

Produces:
- artifacts/m21651/results.json  (machine-readable)
- artifacts/m21651/matrix.md     (human-readable table)
- artifacts/m21651/figures/matrix_heatmap.png  (visual)
"""
from __future__ import annotations

import json
import tempfile
from pathlib import Path

import numpy as np

from tensorcad import mesh_formats as mf

BODIES = (mf.cube(), mf.icosphere(1.0, 2), mf.torus())
FORMATS = mf.MESH_FORMATS
ART_DIR = Path(__file__).resolve().parents[1] / "artifacts" / "m21651"
FIG_DIR = ART_DIR / "figures"


def run_matrix() -> tuple[list[dict], dict]:
    results = []
    for body in BODIES:
        for fmt in FORMATS:
            tol = 1e-5 if fmt == "stl" else 0.0
            with tempfile.TemporaryDirectory() as td:
                path = Path(td) / f"{body.name}.{fmt}"
                report = mf.roundtrip(body, fmt, path, tol=tol)
                results.append(
                    {
                        "body": body.name,
                        "format": fmt,
                        "pass": report["pass"],
                        "volume_ref": report["volume_ref"],
                        "volume_got": report["volume_got"],
                        "volume_rel_err": report["volume_rel_err"],
                        "vertices": report["vertices"],
                        "faces": report["faces"],
                        "euler": report["euler"],
                        "boundary_edges": report["boundary_edges"],
                        "bijection": report["bijection"],
                        "flipped_faces": report["flipped_faces"],
                        "bytes": report["bytes"],
                    }
                )

    summary = {
        "total": len(results),
        "pass": sum(1 for r in results if r["pass"]),
        "fail": sum(1 for r in results if not r["pass"]),
        "worst_own_reader_rel_err": max(r["volume_rel_err"] for r in results),
        "own_reader_topology_pass": sum(
            1 for r in results if r["bijection"] and r["faces"][0] == r["faces"][1]
        ),
        "own_reader_all_within_owner_threshold": all(
            r["volume_rel_err"] < 1e-5 for r in results
        ),
        "cells": sum(1 for r in results if r["pass"]),
    }
    return results, summary


def write_results_json(results: list[dict], summary: dict) -> None:
    ART_DIR.mkdir(parents=True, exist_ok=True)
    data = {
        "summary": summary,
        "results": results,
        "acceptance": {
            "implementation_threshold": 1e-5,
            "description": "6 formats × 3 geometries = 18 combinations, all pass",
        },
    }
    (ART_DIR / "results.json").write_text(json.dumps(data, indent=2), encoding="utf-8")


def write_matrix_md(results: list[dict], summary: dict) -> None:
    lines = [
        "# Native six-format matrix: round-trip fidelity",
        "",
        "## Acceptance gate",
        "",
        f"- **Total combinations**: {summary['total']}",
        f"- **Passed**: {summary['pass']}",
        f"- **Failed**: {summary['fail']}",
        f"- **Worst relative volume error**: {summary['worst_own_reader_rel_err']:.2e}",
        f"- **Topology preserved (all 18)**: {'✅' if summary['own_reader_topology_pass'] == 18 else '❌'}",
        "",
        "## Matrix table",
        "",
        "| Body | Format | Pass | Vol(ref) | Vol(got) | Rel.Err | V(in→out) | F(in→out) | Euler | Bytes |",
        "|------|--------|------|----------|----------|---------|-----------|-----------|-------|-------|",
    ]
    for r in results:
        status = "✅" if r["pass"] else "❌"
        lines.append(
            f"| {r['body']} | {r['format']} | {status} "
            f"| {r['volume_ref']:.6f} | {r['volume_got']:.6f} "
            f"| {r['volume_rel_err']:.2e} "
            f"| {r['vertices'][0]}→{r['vertices'][1]} "
            f"| {r['faces'][0]}→{r['faces'][1]} "
            f"| {r['euler'][0]}→{r['euler'][1]} "
            f"| {r['bytes']} |"
        )
    lines.append("")
    lines.append("## Notes")
    lines.append("")
    lines.append("- STL is float32 by format definition; tolerance 1e-5 applied.")
    lines.append("- All other formats are float64; exact bijection required.")
    lines.append("- CGNS written in SIDS-to-HDF5 dialect with CGNSLibraryVersion node.")
    lines.append("- SU2 written with triangles in NELEM block (volume mode).")
    lines.append("- Gmsh written in MSH 2.2 ASCII format.")
    lines.append("- VTK written as UNSTRUCTURED_GRID (meshio-readable).")
    lines.append("")
    (ART_DIR / "matrix.md").write_text("\n".join(lines), encoding="utf-8")


def generate_heatmap(results: list[dict]) -> bytes:
    """Generate a heatmap PNG showing volume relative error across the matrix."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors

    body_names = [b.name for b in BODIES]
    fmt_names = list(FORMATS)

    err_matrix = np.zeros((len(body_names), len(fmt_names)))
    for r in results:
        i = body_names.index(r["body"])
        j = fmt_names.index(r["format"])
        err_matrix[i, j] = r["volume_rel_err"]

    fig, ax = plt.subplots(figsize=(10, 5))
    norm = mcolors.LogNorm(vmin=1e-15, vmax=1e-5)
    im = ax.imshow(err_matrix, cmap="YlOrRd", norm=norm, aspect="auto")

    ax.set_xticks(range(len(fmt_names)))
    ax.set_xticklabels(fmt_names, fontsize=11)
    ax.set_yticks(range(len(body_names)))
    ax.set_yticklabels(body_names, fontsize=11)

    for i in range(len(body_names)):
        for j in range(len(fmt_names)):
            val = err_matrix[i, j]
            if val == 0:
                label = "0"
            elif val < 1e-12:
                label = "<1e-12"
            else:
                label = f"{val:.1e}"
            ax.text(j, i, label, ha="center", va="center", fontsize=9,
                    color="white" if val > 1e-8 else "black")

    ax.set_xlabel("Format", fontsize=12)
    ax.set_ylabel("Geometry", fontsize=12)
    ax.set_title("Round-trip volume relative error (6 formats × 3 geometries)", fontsize=13)
    cbar = plt.colorbar(im, ax=ax, label="Relative error (log scale)")
    plt.tight_layout()

    from io import BytesIO
    buf = BytesIO()
    plt.savefig(buf, format="png", dpi=150)
    plt.close(fig)
    return buf.getvalue()


def write_figure(results: list[dict]) -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    png_data = generate_heatmap(results)
    (FIG_DIR / "matrix_heatmap.png").write_bytes(png_data)


def main() -> None:
    results, summary = run_matrix()
    write_results_json(results, summary)
    write_matrix_md(results, summary)
    write_figure(results)
    print(f"Matrix: {summary['pass']}/{summary['total']} passed")
    print(f"Worst rel err: {summary['worst_own_reader_rel_err']:.2e}")
    print(f"Artifacts written to {ART_DIR}")


if __name__ == "__main__":
    main()
