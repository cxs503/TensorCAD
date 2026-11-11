"""Figure for ``docs/nurbs_surface_degree_elevation.md`` (four panels).

Reads the artifact produced by ``examples/nurbs_surface_degree_elevation/run_example.py``
and adds the parametric deviation field, which needs a fresh evaluation.

    python tools/nurbs_surface_elevation_figure.py \
        --results examples/nurbs_surface_degree_elevation/outputs/results.json \
        --out examples/nurbs_surface_degree_elevation/outputs/nurbs_surface_degree_elevation.png
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from tensorcad.nurbs import sphere_surface

SPHERE_RADIUS = 2.2


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--results", default="examples/nurbs_surface_degree_elevation/outputs/results.json"
    )
    parser.add_argument(
        "--out",
        default="examples/nurbs_surface_degree_elevation/outputs/nurbs_surface_degree_elevation.png",
    )
    args = parser.parse_args()

    data = json.loads(Path(args.results).read_text())
    cylinder_rows = data["cylinder"]
    sphere_rows = data["sphere"]
    patch_rows = data["non_rational_patch"]["rows"]

    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.5))

    # (a) point deviation, log scale
    ax = axes[0][0]
    groups = [
        ("cylinder", cylinder_rows, "tab:blue"),
        ("sphere", sphere_rows, "tab:orange"),
        ("patch", patch_rows, "tab:green"),
    ]
    positions, values, labels, colors = [], [], [], []
    index = 0
    for _name, rows, color in groups:
        for entry in rows:
            positions.append(index)
            values.append(max(entry["max_point_deviation_m"], 1e-18))
            labels.append(f"{entry['times_u']},{entry['times_v']}")
            colors.append(color)
            index += 1
        index += 0.6
    ax.bar(positions, values, color=colors, width=0.8)
    ax.set_yscale("log")
    ax.axhline(2.2e-16, color="black", linestyle=":", linewidth=1.0, label="machine eps (2.2e-16)")
    ax.set_xticks(list(positions))
    ax.set_xticklabels(labels, rotation=90, fontsize=7)
    ax.set_ylabel(r"$\max |S_{\rm elev}(u,v) - S(u,v)|$  [m]")
    ax.set_title("(a) the elevated patch traces the same surface")
    ax.legend(fontsize=7, loc="lower right")
    for x0, group_name in (
        (0, "cylinder"),
        (len(cylinder_rows) + 0.6, "sphere"),
        (len(cylinder_rows) + len(sphere_rows) + 1.2, "bicubic patch"),
    ):
        ax.text(x0, 3e-14, group_name, fontsize=8)
    ax.set_ylim(1e-18, 1e-12)

    # (b) analytic area error
    ax = axes[0][1]
    entries = [(f"cyl {e['times_u']},{e['times_v']}", e["area_relative_error"], "tab:blue")
               for e in cylinder_rows]
    entries += [(f"sph {e['times_u']},{e['times_v']}", e["area_relative_error"], "tab:orange")
                for e in sphere_rows]
    ax.bar(range(len(entries)), [max(e[1], 1e-18) for e in entries],
           color=[e[2] for e in entries])
    ax.set_yscale("log")
    ax.axhline(2.2e-16, color="black", linestyle=":", linewidth=1.0)
    ax.set_xticks(range(len(entries)))
    ax.set_xticklabels([e[0] for e in entries], rotation=90, fontsize=7)
    ax.set_ylabel(r"$|A_{\rm elev}/A_{\rm analytic} - 1|$")
    ax.set_title(r"(b) area against the analytic $2\pi R H$ / $4\pi R^2$")

    # (c) deviation field on the parametric domain (sphere elevated to (4, 3))
    ax = axes[1][0]
    sphere = sphere_surface(SPHERE_RADIUS)
    elevated = sphere.elevate_degree(2, 1)
    us = np.linspace(*sphere.domain[0], 61)
    vs = np.linspace(*sphere.domain[1], 41)
    field = np.array(
        [
            [math.dist(elevated.evaluate(float(u), float(v)), sphere.evaluate(float(u), float(v)))
             for u in us]
            for v in vs
        ]
    )
    image = ax.pcolormesh(us, vs, field, cmap="viridis", shading="auto")
    fig.colorbar(image, ax=ax, label=r"$|S_{\rm elev} - S|$  [m]")
    ax.set_xlabel("$u$")
    ax.set_ylabel("$v$")
    ax.set_title(r"(c) deviation field, exact sphere elevated to $(4,3)$"
                 f" (max {field.max():.2e} m)")

    # (d) representation growth
    ax = axes[1][1]
    names = ["cylinder u", "cylinder v", "sphere u", "sphere v", "patch u", "patch v"]
    before = [7, 2, 7, 5, 4, 4]
    after = [16, 4, 13, 9, 6, 7]
    x = np.arange(len(names))
    ax.bar(x - 0.2, before, width=0.4, label="original", color="tab:gray")
    ax.bar(x + 0.2, after, width=0.4, label="elevated (times = 3)", color="tab:red")
    ax.set_xticks(x)
    ax.set_xticklabels(names, fontsize=8, rotation=20)
    ax.set_ylabel("control points in that direction")
    ax.set_title("(d) minimal growth: points + times x spans")
    ax.legend(fontsize=8)

    fig.tight_layout()
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
