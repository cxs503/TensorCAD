"""Reproducible evidence generator for the `tensorcad.intersect` verification.

Writes `artifacts/m6300_intersect/`:
  * `case_matrix.json`  -- analytic cases with the library verdict
  * `differential.json`  -- randomised comparison against the exact rational routes
  * `mesh_controls.json` -- mesh-level negative/positive controls and pruning checks
  * `pytest.txt`         -- the suite log (filled in by the run script)
  * `fig_intersect.png`  -- figure for the report
"""
from __future__ import annotations

import json
import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "tests"))

from intersect_oracle import (  # noqa: E402
    triangles_intersect_clip3d,
    triangles_intersect_exact,
)

from tensorcad.geometry import Mesh  # noqa: E402
from tensorcad.intersect import find_self_intersections, triangles_intersect  # noqa: E402

OUT = ROOT / "artifacts" / "m6300_intersect"
OUT.mkdir(parents=True, exist_ok=True)

DYADIC = (-1.0, -0.5, 0.0, 0.5, 1.0)
THIRDS = (-1.0, -1.0 / 3.0, 0.0, 1.0 / 3.0, 1.0)

FLAT = ((0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (0.0, 2.0, 0.0))
PIERCING = ((0.5, 0.5, -1.0), (0.5, 0.5, 1.0), (1.5, 1.5, 0.5))
LENS_A = ((1.0, -1.0, 0.0), (-0.5, 0.5, 0.5), (1.0, -1.0, -0.5))
LENS_B = ((0.5, -0.5, -0.5), (0.0, 0.0, 0.5), (1.0, -1.0, 0.5))
STAR_UP = ((0.0, 1.0, 0.0), (-1.0, -0.5, 0.0), (1.0, -0.5, 0.0))
STAR_DOWN = ((0.0, -1.0, 0.0), (1.0, 0.5, 0.0), (-1.0, 0.5, 0.0))
TOUCH_PLANE = ((1.0, 0.0, -1.0), (0.5, -0.5, -0.5), (0.0, 1.0, 0.0))
TOUCH_PLANE_B = ((1.0, -0.5, 0.0), (-0.5, 1.0, 0.5), (0.5, 1.0, -0.5))

CASES = [
    ("piercing pair", FLAT, PIERCING, True),
    ("far apart", FLAT, ((5.0, 5.0, 0.0), (6.0, 5.0, 0.0), (5.0, 6.0, 0.0)), False),
    ("1 mm gap blade", FLAT, ((0.5, 0.5, 1e-3), (0.5, 0.5, 2.0), (1.5, 1.5, 1.5)), False),
    ("shared vertex", FLAT, ((0.0, 0.0, 0.0), (0.0, 0.0, 1.0), (1.0, 1.0, 1.0)), True),
    ("shared edge", FLAT, ((0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (0.0, 1.0, 1.0)), True),
    ("coplanar lens (regression)", LENS_A, LENS_B, True),
    ("star of david (regression)", STAR_UP, STAR_DOWN, True),
    ("vertex resting on plane (regression)", TOUCH_PLANE, TOUCH_PLANE_B, True),
    ("coplanar disjoint", LENS_A, ((3.0, -3.0, 0.0), (4.0, -4.0, 0.0), (3.0, -3.0, 1.0)), False),
    ("parallel offset planes", FLAT, tuple((x, y, z + 0.5) for x, y, z in FLAT), False),
    ("collinear sliver", FLAT, ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0)), False),
]


def rand_tri(rng, grid):
    return tuple(tuple(rng.choice(grid) for _ in range(3)) for _ in range(3))


def differential(n=2000, seed=20261107, grid=DYADIC, exact_only=False):
    rng = random.Random(seed)
    rows = {"n": n, "grid": "dyadic" if grid is DYADIC else "thirds",
            "library_vs_line_clip": 0, "library_vs_halfspace": 0,
            "line_vs_halfspace": 0, "scored": 0}
    for _ in range(n):
        first = rand_tri(rng, grid)
        second = rand_tri(rng, grid)
        lib = triangles_intersect(first, second)
        line = triangles_intersect_exact(first, second)
        clip = triangles_intersect_clip3d(first, second)
        if line != clip:
            rows["line_vs_halfspace"] += 1
            if exact_only:
                continue
        rows["scored"] += 1
        rows["library_vs_line_clip"] += lib != line
        rows["library_vs_halfspace"] += lib != clip
    return rows


def octahedron() -> Mesh:
    v = ((1.0, 0, 0), (-1.0, 0, 0), (0, 1.0, 0), (0, -1.0, 0), (0, 0, 1.0), (0, 0, -1.0))
    f = ((0, 2, 4), (2, 1, 4), (1, 3, 4), (3, 0, 4), (2, 0, 5), (1, 2, 5), (3, 1, 5), (0, 3, 5))
    return Mesh(vertices=v, faces=f)


def crossing_sheets() -> Mesh:
    return Mesh(
        vertices=((0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (0.0, 2.0, 0.0),
                  (0.5, 0.5, -1.0), (0.5, 0.5, 1.0), (1.5, 1.5, 0.5)),
        faces=((0, 1, 2), (3, 4, 5)),
    )


def warped_grid(rows: int, columns: int) -> Mesh:
    vertices = []
    for r in range(rows + 1):
        for c in range(columns + 1):
            z = 0.6 if (0 < r < rows and 0 < c < columns and (r + c) % 2 == 0) else 0.0
            vertices.append((float(c), float(r), z))
    stride = columns + 1
    faces = []
    for r in range(rows):
        for c in range(columns):
            a = r * stride + c
            faces.append((a, a + 1, a + stride + 1))
            faces.append((a, a + stride + 1, a + stride))
    return Mesh(vertices=tuple(vertices), faces=tuple(faces))


def brute_force(mesh: Mesh):
    rows = []
    for i in range(len(mesh.faces)):
        for j in range(i + 1, len(mesh.faces)):
            if set(mesh.faces[i]) & set(mesh.faces[j]):
                continue
            if triangles_intersect([mesh.vertices[k] for k in mesh.faces[i]],
                                   [mesh.vertices[k] for k in mesh.faces[j]]):
                rows.append((i, j))
    return tuple(sorted(rows))


def mesh_controls():
    out = {}
    out["clean_octahedron_pairs"] = len(find_self_intersections(octahedron()))
    out["clean_octahedron_cells4_pairs"] = len(find_self_intersections(octahedron(), cells_per_axis=4))
    out["crossing_sheets_pairs"] = list(find_self_intersections(crossing_sheets()))
    mesh = warped_grid(6, 6)
    expected = brute_force(mesh)
    grid_found = find_self_intersections(mesh, cells_per_axis=8)
    out["warped_grid_faces"] = len(mesh.faces)
    out["warped_grid_brute_force"] = len(expected)
    out["warped_grid_pruned"] = len(grid_found)
    out["warped_grid_pruning_agrees"] = tuple(grid_found) == expected
    return out


def figure(case_rows, diff_dyadic, diff_thirds):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon

    fig = plt.figure(figsize=(11.0, 4.2))
    ax1 = fig.add_subplot(1, 3, 1)
    ax1.add_patch(Polygon([(p[0], p[2]) for p in LENS_A], closed=True, alpha=0.45,
                          edgecolor="tab:blue", facecolor="tab:blue"))
    ax1.add_patch(Polygon([(p[0], p[2]) for p in LENS_B], closed=True, alpha=0.45,
                          edgecolor="tab:red", facecolor="tab:red"))
    ax1.set_title("coplanar lens overlap\nno vertex inside the other", fontsize=9)
    ax1.set_xlabel("x (m)")
    ax1.set_ylabel("z (m)  [plane x+y=0]")
    ax1.grid(alpha=0.3)

    ax2 = fig.add_subplot(1, 3, 2)
    labels = [name for name, *_ in case_rows]
    verdicts = [1 if row["matches"] else 0 for row in case_rows]
    colours = ["tab:green" if v else "tab:gray" for v in verdicts]
    ax2.barh(range(len(labels)), verdicts, color=colours)
    ax2.set_yticks(range(len(labels)))
    ax2.set_yticklabels(labels, fontsize=7)
    ax2.set_xlim(0, 1.25)
    ax2.set_title("analytic case matrix\n(1 = intersects, all as constructed)", fontsize=9)
    ax2.grid(alpha=0.3, axis="x")

    ax3 = fig.add_subplot(1, 3, 3)
    names = ["dyadic\n(exact in binary)", "thirds\n(rounded)"]
    lib_vs_clip = [diff_dyadic["library_vs_halfspace"], diff_thirds["library_vs_halfspace"]]
    lib_vs_line = [diff_dyadic["library_vs_line_clip"], diff_thirds["library_vs_line_clip"]]
    xs = range(len(names))
    ax3.bar([x - 0.18 for x in xs], lib_vs_clip, width=0.34, label="vs half-space clip", color="tab:blue")
    ax3.bar([x + 0.18 for x in xs], lib_vs_line, width=0.34, label="vs plane-line clip", color="tab:orange")
    for x, v in zip(xs, lib_vs_clip):
        ax3.text(x - 0.18, v + 0.2, str(v), ha="center", fontsize=8)
    for x, v in zip(xs, lib_vs_line):
        ax3.text(x + 0.18, v + 0.2, str(v), ha="center", fontsize=8)
    ax3.set_xticks(list(xs))
    ax3.set_xticklabels(names, fontsize=8)
    ax3.set_ylabel("disagreeing pairs / 2000 resp. 600")
    ax3.set_title("differential test vs exact routes", fontsize=9)
    ax3.legend(fontsize=7)
    ax3.grid(alpha=0.3, axis="y")

    fig.tight_layout()
    path = OUT / "fig_intersect.png"
    fig.savefig(path, dpi=150)
    print("figure", path)


if __name__ == "__main__":
    case_rows = []
    for name, first, second, expected in CASES:
        got = triangles_intersect(first, second)
        case_rows.append({
            "case": name,
            "expected": expected,
            "library": got,
            "matches": got == expected,
            "oracle_line_clip": triangles_intersect_exact(first, second),
            "oracle_clip3d": triangles_intersect_clip3d(first, second),
        })
    (OUT / "case_matrix.json").write_text(json.dumps(case_rows, indent=2))

    diff_dyadic = differential()
    diff_thirds = differential(n=600, seed=7, grid=THIRDS)
    (OUT / "differential.json").write_text(json.dumps(
        {"dyadic": diff_dyadic, "thirds": diff_thirds}, indent=2))

    controls = mesh_controls()
    (OUT / "mesh_controls.json").write_text(json.dumps(controls, indent=2))

    print(json.dumps({"cases_all_match": all(r["matches"] for r in case_rows),
                      "dyadic": diff_dyadic, "thirds": diff_thirds,
                      "controls": controls}, indent=1))
    figure(case_rows, diff_dyadic, diff_thirds)