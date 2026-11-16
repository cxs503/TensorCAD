"""Generate the foreign-dialect exchange evidence (table + JSON + figure).

Writes, under ``runs/m22862_foreign_dialects/``:

* ``matrix.json`` -- one record per (geometry, dialect, producer) with the
  measured volume, relative error, tolerance, cell count, winding flag, the
  byte-level identity facts and the sha256 of each file;
* ``matrix.md``  -- the same as the acceptance table;
* ``files/``     -- every exported mesh, so a third party can re-read the exact
  bytes without re-running this code;
* ``consumer.json`` -- the trimesh (third-party) verdicts on TensorCAD-written
  files;

and under ``artifacts/m22862_foreign_dialects/`` the figure.

Usage::

    python scripts/mesh_interop_foreign_matrix.py
"""
from __future__ import annotations

import hashlib
import json
import platform
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np  # noqa: E402

from tensorcad.mesh_foreign_dialects import (  # noqa: E402
    STL_BINARY_RELATIVE_TOLERANCE,
    foreign_matrix_markdown,
    run_foreign_matrix,
    trimesh_consumer_rows,
)
from tensorcad.mesh_interop_thirdparty import trimesh_available, trimesh_version  # noqa: E402

RUN_DIR = ROOT / "runs" / "m22862_foreign_dialects"
FIG_DIR = ROOT / "artifacts" / "m22862_foreign_dialects"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 16), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _figure(rows: list[dict], consumer: list[dict]) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    FIG_DIR.mkdir(parents=True, exist_ok=True)
    geometries = sorted({r["geometry"] for r in rows})
    dialects = ["stl_binary", "gmsh22", "vtu"]
    colours = {"stl_binary": "#c0392b", "gmsh22": "#2980b9", "vtu": "#27ae60"}
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.5, 5.2))

    # ---- left: exchange error per dialect, log scale, versus tolerance ------
    floor = 1e-18
    width = 0.26
    for i, dialect in enumerate(dialects):
        xs, ys = [], []
        for j, geom in enumerate(geometries):
            vals = [
                max(r["volume_relative_error"], floor)
                for r in rows
                if r["geometry"] == geom and r["dialect"] == dialect and r["cell_kind"] == "tri"
            ]
            xs.append(j + (i - 1) * width)
            ys.append(max(vals) if vals else np.nan)
        ax1.bar(xs, ys, width=width, color=colours[dialect], label=dialect)
    ax1.set_yscale("log")
    ax1.set_ylim(floor, 1e-3)
    ax1.axhline(1e-12, color="k", ls="--", lw=1.2, label="float64 tolerance 1e-12")
    ax1.axhline(
        STL_BINARY_RELATIVE_TOLERANCE,
        color="#c0392b",
        ls=":",
        lw=1.4,
        label=f"binary STL float32 bound {STL_BINARY_RELATIVE_TOLERANCE:.0e}",
    )
    ax1.set_xticks(range(len(geometries)))
    ax1.set_xticklabels(geometries)
    ax1.set_ylabel("|V_out - V_ref| / V_ref  (log)")
    ax1.set_title("Foreign-dialect ingest: volume error\n(worst producer per cell)")
    ax1.legend(fontsize=8, loc="upper left")
    ax1.grid(axis="y", alpha=0.3)

    # ---- right: third-party (trimesh) volume vs reference -------------------
    if consumer:
        markers = {"stl": "o", "stl_binary": "s", "obj": "^"}
        for row in consumer:
            ref = row["volume_reference"]
            ax2.scatter(
                ref,
                row["volume_thirdparty"],
                marker=markers[row["dialect"]],
                s=70,
                label=f"trimesh reads {row['dialect']}",
            )
        lo = min(r["volume_reference"] for r in consumer) * 0.9
        hi = max(r["volume_reference"] for r in consumer) * 1.1
        ax2.plot([lo, hi], [lo, hi], "k--", lw=1, label="identity")
        ax2.set_xscale("log")
        ax2.set_yscale("log")
        ax2.set_xlim(lo, hi)
        ax2.set_ylim(lo, hi)
        ax2.set_xlabel("reference volume [m^3]")
        ax2.set_ylabel("trimesh volume [m^3]")
        ax2.set_title(
            "Independent implementation (trimesh "
            f"{trimesh_version()}) reading TensorCAD files"
        )
        handles, labels = ax2.get_legend_handles_labels()
        unique = dict(zip(labels, handles))
        ax2.legend(unique.values(), unique.keys(), fontsize=8, loc="upper left")
        ax2.grid(alpha=0.3, which="both")
    else:  # pragma: no cover - environment dependent
        ax2.text(0.5, 0.5, "trimesh not available", ha="center")

    fig.suptitle(
        "GAP-CAD-B2-MESH-FORMAT-MATRIX: foreign dialects (binary STL / MSH 2.2 / VTU) "
        "and a real third-party reader",
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    png = FIG_DIR / "foreign_dialects.png"
    fig.savefig(png, dpi=140)
    plt.close(fig)
    return png


def main() -> int:
    if RUN_DIR.exists():
        shutil.rmtree(RUN_DIR)
    (RUN_DIR / "files").mkdir(parents=True, exist_ok=True)
    rows = run_foreign_matrix(workdir=RUN_DIR / "files")
    consumer = trimesh_consumer_rows(workdir=RUN_DIR / "consumer")

    for path in sorted((RUN_DIR / "files").glob("*")):
        for row in rows:
            tag = f"{row['geometry']}.{row['producer']}.{row['dialect']}"
            if path.stem.startswith(tag):
                row["sha256"] = _sha256(path)
                row["file"] = path.name
    for path in sorted((RUN_DIR / "consumer").glob("*")):
        for row in consumer:
            tag = f"{row['geometry']}.consumer.{row['dialect']}"
            if path.stem.startswith(tag):
                row["sha256"] = _sha256(path)
                row["file"] = path.name

    (RUN_DIR / "matrix.json").write_text(json.dumps(rows, indent=1))
    (RUN_DIR / "matrix.md").write_text(foreign_matrix_markdown(rows))
    (RUN_DIR / "consumer.json").write_text(json.dumps(consumer, indent=1))
    png = _figure(rows, consumer)

    meta = {
        "mission": 22862,
        "anchor": "GAP-CAD-B2-MESH-FORMAT-MATRIX",
        "rows": len(rows),
        "rows_ok": sum(1 for r in rows if r["ok"]),
        "consumer_rows": len(consumer),
        "consumer_rows_ok": sum(1 for r in consumer if r["ok"]),
        "trimesh": trimesh_version() if trimesh_available() else None,
        "python": platform.python_version(),
        "platform": platform.platform(),
        "stl_binary_tolerance": STL_BINARY_RELATIVE_TOLERANCE,
        "figure": str(png.relative_to(ROOT)),
    }
    (RUN_DIR / "meta.json").write_text(json.dumps(meta, indent=1))
    print(json.dumps(meta, indent=1))
    print((RUN_DIR / "matrix.md").read_text())
    return 0 if meta["rows_ok"] == meta["rows"] and meta["consumer_rows_ok"] == meta["consumer_rows"] else 1


if __name__ == "__main__":
    raise SystemExit(main())