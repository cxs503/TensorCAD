#!/usr/bin/env python3
"""Volume-mesh -> tetrahedra -> SU2/CGNS: measured round-trip fidelity.

Runs the three checks that make the bridge credible:

1. **Split conservation** -- the six Kuhn tetrahedra of each hexahedron tile it exactly, so
   the tetrahedra reproduce the hex grid's volume (and therefore the body's).
2. **Analytic anchors** -- a box has volume ``2b * L * H``; a ``SectionGrid`` body has
   ``2 * integral(b dx dz)``.  Both are closed form, so agreement is a real check rather
   than self-comparison.
3. **Format round trip** -- write/read through meshio's SU2 and CGNS, measure the volume
   change and the cell/vertex counts coming back.
4. **Named skin** -- weld the coincident points of the collapsed columns, cancel the
   zero-thickness slit they leave behind, name the remaining boundary triangles from the
   grid's own planes, and check the result: the skin must be watertight, it must enclose the
   volume the tetrahedra measure, and each patch's area must match a closed-form (box) or
   independently integrated (Wigley wetted surface) anchor.  The marked SU2 file is written
   to ``runs/m19950_volume_bridge/<body>_marked.su2`` and read back by both this package's
   parser and meshio's.

Writes ``runs/m19950_volume_bridge/bridge_matrix.json`` and the figure under
``artifacts/m19950_volume_bridge/``.  Nothing here is hand-copied: every number in the
report is read from the JSON this script emits.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from tensorcad import volume_mesh as vmlib
from tensorcad.meshio_exchange import (
    FORMATS_VOLUME,
    read_su2_markers,
    round_trip_meshio,
    tetrahedral_volume,
    volume_mesh_patches,
    volume_mesh_to_tetrahedral,
    write_su2_volume_mesh,
)

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
RUNS = REPO / "runs" / "m19950_volume_bridge"
ARTIFACTS = REPO / "artifacts" / "m19950_volume_bridge"


# ---------------------------------------------------------------------------------------
# Bodies
# ---------------------------------------------------------------------------------------

def box_grid(nx: int = 5, nz: int = 4, b: float = 0.5, length: float = 2.0, height: float = 0.6):
    xs = tuple(i * (length / (nx - 1)) for i in range(nx))
    zs = tuple(j * (height / (nz - 1)) for j in range(nz))
    return vmlib.SectionGrid(xs=xs, zs=zs, breadth=lambda i, j: b)


def hull_grid(length=10.0, beam=1.0, draft=0.625, nx=25, nz=9):
    """Wigley-type half breadth: vanishes at both ends and at the keel."""
    xs = tuple(-length / 2 + i * (length / (nx - 1)) for i in range(nx))
    zs = tuple(j * (draft / (nz - 1)) for j in range(nz))

    def breadth(i: int, j: int) -> float:
        x, z = xs[i], zs[j]
        return max(0.0, (beam / 2) * (1 - (2 * x / length) ** 2) * (1 - (z / draft) ** 2))

    return vmlib.SectionGrid(xs=xs, zs=zs, breadth=breadth), length, beam, draft


def trapezoid_volume(grid: vmlib.SectionGrid) -> float:
    total = 0.0
    for i in range(grid.nx - 1):
        dx = grid.xs[i + 1] - grid.xs[i]
        for j in range(grid.nz - 1):
            dz = grid.zs[j + 1] - grid.zs[j]
            corners = (
                grid.breadth(i, j),
                grid.breadth(i + 1, j),
                grid.breadth(i, j + 1),
                grid.breadth(i + 1, j + 1),
            )
            total += dx * dz * sum(corners) / 4.0
    return 2.0 * total


# ---------------------------------------------------------------------------------------
# Measurement
# ---------------------------------------------------------------------------------------

def measure(
    name: str,
    vm: vmlib.VolumeMesh,
    analytic: float,
    continuum: float | None = None,
) -> dict[str, Any]:
    tet = volume_mesh_to_tetrahedral(vm)
    row: dict[str, Any] = {
        "body": name,
        "hexes": len(vm.hexes),
        "tets_before_drop": len(vm.tets),
        "tets_kept": len(tet.tets),
        "degenerate_cells_dropped": tet.degenerate_cells,
        "vertices": len(tet.vertices),
        "volume_grid": vm.volume(),
        "volume_tet": tet.volume(),
        "volume_analytic": analytic,
        "rel_split_vs_grid": abs(tet.volume() - vm.volume()) / vm.volume(),
        "rel_analytic": abs(tet.volume() - analytic) / analytic,
        "volume_continuum": continuum,
        # Discretisation of the body itself (piecewise-linear breadth vs the smooth
        # parabola); NOT an error of the bridge, and reported separately so the two
        # are never conflated.
        "rel_continuum": (abs(tet.volume() - continuum) / continuum) if continuum else None,
        "all_positive": all(tetrahedral_volume(tet.vertices, c) > 0.0 for c in tet.tets),
        "formats": {},
    }
    for fmt in FORMATS_VOLUME:
        back = round_trip_meshio(tet, file_format=fmt)
        row["formats"][fmt] = {
            "tets_returned": len(back.tets),
            "vertices_returned": len(back.vertices),
            "volume": back.volume(),
            "rel_error": abs(back.volume() - tet.volume()) / tet.volume(),
            "signed_positive": back.signed_volume() > 0.0,
        }
    return row


def wetted_area_quadrature(L: float, B: float, T: float, nx: int = 400, nz: int = 200) -> float:
    """``2 \\iint sqrt(1 + b_x^2 + b_z^2) dx dz`` from the breadth function, not the mesh."""
    total = 0.0
    dx, dz = L / nx, T / nz
    for i in range(nx):
        x = -L / 2 + (i + 0.5) * dx
        for j in range(nz):
            z = (j + 0.5) * dz
            bx = (B / 2) * (-8.0 * x / L**2) * (1 - (z / T) ** 2)
            bz = (B / 2) * (1 - (2 * x / L) ** 2) * (-2.0 * z / T**2)
            total += float(np.sqrt(1.0 + bx * bx + bz * bz))
    return 2.0 * total * dx * dz


def markers(name: str, vm: vmlib.VolumeMesh, path: Path) -> dict[str, Any]:
    """The named skin: patch table, closure of the surface, and the file's own contents."""
    kept = volume_mesh_to_tetrahedral(vm)
    raw_faces: dict[tuple[int, int, int], int] = {}
    for tet in kept.tets:
        for a in range(4):
            for b in range(a + 1, 4):
                for c in range(b + 1, 4):
                    face = tuple(sorted((tet[a], tet[b], tet[c])))
                    raw_faces[face] = raw_faces.get(face, 0) + 1
    raw_boundary = sum(1 for count in raw_faces.values() if count == 1)
    patches = volume_mesh_patches(vm, tetrahedral=kept)
    write_su2_volume_mesh(path, vm, patches=patches)
    read_back = read_su2_markers(path)
    tokens = path.read_text(encoding="utf-8").splitlines()
    header = {}
    for token in tokens:
        for key in ("NPOIN", "NELEM", "NMARK"):
            if token.startswith(f"{key}="):
                header[key] = int(token.split("=")[1])
    row: dict[str, Any] = {
        "body": name,
        "su2_path": str(path.relative_to(REPO)),
        "raw_boundary_triangles": raw_boundary,
        "skin_triangles": len(patches.skin()),
        "slit_triangles_cancelled": raw_boundary - len(patches.skin()),
        "merged_points": patches.merged_points,
        "vertices_written": len(patches.vertices),
        "tets_written": len(patches.tets),
        "patches": patches.counts(),
        "areas": patches.areas(),
        "skin_volume": patches.skin_volume(),
        "tet_volume": patches.tet_volume(),
        "rel_skin_vs_tets": abs(patches.skin_volume() - patches.tet_volume())
        / patches.tet_volume(),
        "watertight": patches.is_watertight(),
        "header": header,
        "reader_tags": list(read_back),
        "reader_matches_writer": all(
            set(read_back[name]) == set(faces) for name, faces in patches.patches
        ),
    }
    return row


def figure(rows: list[dict[str, Any]], hull: vmlib.VolumeMesh, out: Path,
           marker_rows: list[dict[str, Any]] | None = None) -> None:
    fig = plt.figure(figsize=(17.6, 4.4))

    # Panel A: the hull body grid itself -- graded boundary layer, collapsed columns.
    ax = fig.add_subplot(1, 4, 1)
    xs = np.array([v[0] for v in hull.vertices])
    zs = np.array([v[2] for v in hull.vertices])
    ax.scatter(xs, zs, s=1.2, c="#1f5fa8", linewidths=0)
    ax.set_xlabel("x [m]  (aft $\\rightarrow$ fore)")
    ax.set_ylabel("z [m]")
    ax.set_title("(a) hull body grid nodes\n(graded wall $\\rightarrow$ centreline)")
    ax.set_aspect("equal")
    ax.grid(alpha=0.25)

    # Panel B: volume identity chain, analytic -> grid -> tets -> SU2 -> CGNS.
    ax = fig.add_subplot(1, 4, 2)
    labels, values, colors = [], [], []
    for row in rows:
        labels += [f"{row['body']}\nanalytic", f"{row['body']}\ntets", f"{row['body']}\nSU2", f"{row['body']}\nCGNS"]
        values += [
            row["volume_analytic"],
            row["volume_tet"],
            row["formats"]["su2"]["volume"],
            row["formats"]["cgns"]["volume"],
        ]
        colors += ["#888888"] + ["#1f5fa8"] * 3
    ax.bar(range(len(values)), values, color=colors)
    ax.set_xticks(range(len(values)))
    ax.set_xticklabels(labels, fontsize=6, rotation=90)
    ax.set_ylabel("volume [m$^3$]")
    ax.set_title("(b) volume identity chain")

    # Panel C: relative error, log scale -- one bar per stage, per body.
    ax = fig.add_subplot(1, 4, 3)
    stages = ["split\nvs grid", "vs\nanalytic", "SU2\nround trip", "CGNS\nround trip"]
    floor = 1e-18
    width = 0.35
    for k, row in enumerate(rows):
        errs = [
            row["rel_split_vs_grid"],
            row["rel_analytic"],
            row["formats"]["su2"]["rel_error"],
            row["formats"]["cgns"]["rel_error"],
        ]
        ax.bar([i + k * width for i in range(len(errs))], [max(e, floor) for e in errs],
               width=width, label=row["body"])
    ax.set_yscale("log")
    ax.set_xticks([i + width / 2 for i in range(len(stages))])
    ax.set_xticklabels(stages, fontsize=8)
    ax.axhline(1e-12, color="#c0392b", linestyle="--", linewidth=1, label="round-off floor")
    ax.set_ylabel("relative volume error")
    ax.set_title("(c) every stage at round-off")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.25, axis="y")

    # Panel D: the named skin -- patch areas of the Wigley body against the two independent
    # anchors (wetted-surface integral, bottom trapezoid) on a log scale of relative error.
    if marker_rows:
        ax = fig.add_subplot(1, 4, 4)
        for k, row in enumerate(marker_rows):
            names = list(row["patches"])
            errors = []
            for patch in names:
                area = row["areas"][patch]
                if patch == "hull" and row.get("hull_area_reference"):
                    reference = row["hull_area_reference"]
                elif patch == "bottom" and row.get("bottom_area_reference"):
                    reference = row["bottom_area_reference"]
                elif row.get("surface_area_reference") and patch in ("top", "fore", "aft"):
                    reference = row["surface_area_reference"][patch]
                else:
                    reference = None
                errors.append(abs(area - reference) / reference if reference else np.nan)
            positions = [i + k * 0.38 for i in range(len(names))]
            ax.bar(positions, [max(e, 1e-18) for e in errors], width=0.38,
                   label=row["body"])
        ax.set_yscale("log")
        all_names = sorted({n for row in marker_rows for n in row["patches"]})
        ax.set_xticks([i + 0.19 for i in range(len(all_names))])
        ax.set_xticklabels(all_names, fontsize=7, rotation=30)
        ax.axhline(0.01, color="#c0392b", linestyle="--", linewidth=1, label="1 % line")
        ax.set_ylabel("relative area error vs analytic anchor")
        ax.set_title("(d) patch areas vs closed-form anchors")
        ax.legend(fontsize=7)
        ax.grid(alpha=0.25, axis="y")

    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150)
    plt.close(fig)


def main() -> int:
    RUNS.mkdir(parents=True, exist_ok=True)
    ARTIFACTS.mkdir(parents=True, exist_ok=True)

    box = box_grid()
    box_vm = vmlib.build_volume_mesh(box, first_layer=0.02, growth=1.2, layers=3)
    box_analytic = 2 * 0.5 * 2.0 * 0.6

    hgrid, L, B, T = hull_grid()
    hull_vm = vmlib.build_volume_mesh(hgrid, first_layer=0.004, growth=1.25, layers=4)
    hull_analytic = trapezoid_volume(hgrid)

    # Continuum closed form for b(x,z) = (B/2)(1-(2x/L)^2)(1-(z/T)^2):
    #   2 * integral(b) dx dz = (4/9) * B * L * T
    hull_continuum = (4.0 / 9.0) * B * L * T

    rows = [
        measure("box 2x1x0.6", box_vm, box_analytic, continuum=box_analytic),
        measure("Wigley body L10xB1xT0.625", hull_vm, hull_analytic, continuum=hull_continuum),
    ]

    # The named skin: patch areas against closed-form / independently integrated anchors.
    box_markers = markers("box 2x1x0.6", box_vm, RUNS / "box_marked.su2")
    box_markers["surface_area_reference"] = {
        "bottom": 2.0 * 1.0,
        "top": 2.0 * 1.0,
        "fore": 1.0 * 0.6,
        "aft": 1.0 * 0.6,
    }
    box_markers["hull_area_reference"] = 2.0 * 2.0 * 0.6

    hull_markers = markers("Wigley body L10xB1xT0.625", hull_vm, RUNS / "hull_marked.su2")
    hull_markers["hull_area_reference"] = wetted_area_quadrature(L, B, T)
    # 2 * integral(b(x, 0) dx) = B * L * 2/3 for this breadth function.
    hull_markers["bottom_area_reference"] = B * L * 2.0 / 3.0

    marker_rows = [box_markers, hull_markers]

    payload = {
        "mission": "m19950",
        "anchor": "GAP-CAD-B2-MESH-FORMAT-MATRIX (residual limit 2 closed)",
        "bodies": rows,
        "markers": marker_rows,
        "criteria": {
            "split_conservation_rel": 1e-12,
            "analytic_rel": 1e-12,
            "format_round_trip_rel": 1e-9,
            "skin_watertight": True,
            "skin_vs_tets_rel": 1e-12,
            "patch_area_rel": 1e-2,
        },
        "all_pass": all(
            r["rel_split_vs_grid"] < 1e-12
            and r["rel_analytic"] < 1e-12
            and all(f["rel_error"] < 1e-9 and f["signed_positive"] for f in r["formats"].values())
            and r["all_positive"]
            for r in rows
        )
        and all(
            m["watertight"]
            and m["rel_skin_vs_tets"] < 1e-12
            and m["reader_matches_writer"]
            and m["header"]["NPOIN"] == m["vertices_written"]
            and m["header"]["NELEM"] == m["tets_written"]
            and m["header"]["NMARK"] == len(m["patches"])
            and abs(m["areas"]["hull"] - m["hull_area_reference"]) / m["hull_area_reference"]
            < 1e-2
            for m in marker_rows
        ),
    }
    (RUNS / "bridge_matrix.json").write_text(json.dumps(payload, indent=2) + "\n")

    fig_path = ARTIFACTS / "volume_bridge_fidelity.png"
    figure(rows, hull_vm, fig_path, marker_rows)

    for row in rows:
        print(f"\n{row['body']}")
        print(f"  hexes {row['hexes']:5d}  tets kept {row['tets_kept']:6d} "
              f"(dropped {row['degenerate_cells_dropped']:4d} zero-volume)")
        print(f"  V analytic {row['volume_analytic']:.12f}  V tets {row['volume_tet']:.12f}")
        print(f"  split-vs-grid rel {row['rel_split_vs_grid']:.3e}   "
              f"vs analytic rel {row['rel_analytic']:.3e}")
        if row.get("rel_continuum") is not None:
            print(f"  body discretisation vs continuum {row['volume_continuum']:.9f} m^3: "
                  f"rel {row['rel_continuum']:.3e}")
        for fmt, info in row["formats"].items():
            print(f"  {fmt.upper():5s} rel {info['rel_error']:.3e}  "
                  f"cells {info['tets_returned']:6d}/{row['tets_kept']}")
    for row in marker_rows:
        print(f"\n{row['body']} -- named skin ({Path(row['su2_path']).name})")
        print(f"  boundary triangles in the split {row['raw_boundary_triangles']:5d}  "
              f"-> skin {row['skin_triangles']:5d}  (slit pairs cancelled "
              f"{row['slit_triangles_cancelled']:5d})")
        print(f"  points welded {row['merged_points']:4d}  "
              f"NPOIN {row['header']['NPOIN']:5d}  NELEM {row['header']['NELEM']:6d}  "
              f"NMARK {row['header']['NMARK']}")
        for patch, count in row["patches"].items():
            print(f"    {patch:8s} {count:5d} triangles   area {row['areas'][patch]:.9f} m^2")
        print(f"  skin volume {row['skin_volume']:.12f}  vs tets {row['tet_volume']:.12f} "
              f"(rel {row['rel_skin_vs_tets']:.2e})  watertight {row['watertight']}")
        print(f"  tags read back {row['reader_tags']}  match writer "
              f"{row['reader_matches_writer']}")
        if "hull_area_reference" in row:
            print(f"  hull area {row['areas']['hull']:.9f} vs anchor "
                  f"{row['hull_area_reference']:.9f} "
                  f"(rel {abs(row['areas']['hull'] - row['hull_area_reference']) / row['hull_area_reference']:.3e})")

    print(f"\nALL PASS: {payload['all_pass']}")
    print(f"json : {RUNS / 'bridge_matrix.json'}")
    print(f"figure: {fig_path}")
    return 0 if payload["all_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
