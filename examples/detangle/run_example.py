"""Measure the self-intersection heal of ``tensorcad.detangle`` and fail loudly on a slip.

Writes ``artifacts/detangle/results.json``, ``plot.png`` and ``report.md``; exits non-zero
when any acceptance claim stops holding. Every number in ``docs/detangle.md`` is read from
that artifact.

The fixtures are analytic up to tessellation: a box (volume exactly ``l*w*h``), a UV sphere
(``4/3 pi r^3``, and the fixture's own tessellation error is reported next to it) and the
parametric hull, whose analytic volume the hull audit already records. The defect is then
*injected* -- one vertex pushed until the surface it owns pierces the opposite side -- and
the heal has to remove the penetration without moving the enclosed volume.
"""
from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

from tensorcad.detangle import excavation_faces, heal_self_intersections
from tensorcad.geometry import Mesh
from tensorcad.hull import HullParameters, analytic_volume, build_mesh
from tensorcad.intersect import find_self_intersections
from tensorcad.mesh_repair import measure_mesh, repair_mesh

HERE = Path(__file__).resolve().parent
OUT = HERE.parents[1] / "artifacts" / "detangle"

LENGTH, BEAM, DRAFT = 10.0, 2.0, 0.8
VOLUME_BUDGET = 0.01


# --------------------------------------------------------------------------- fixtures


def _outward(vertices, tri, centre):
    ax, ay, az = vertices[tri[0]]
    bx, by, bz = vertices[tri[1]]
    cx, cy, cz = vertices[tri[2]]
    ux, uy, uz = bx - ax, by - ay, bz - az
    vx, vy, vz = cx - ax, cy - ay, cz - az
    nx, ny, nz = uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx
    mx, my, mz = (ax + bx + cx) / 3 - centre[0], (ay + by + cy) / 3 - centre[1], (
        az + bz + cz
    ) / 3 - centre[2]
    return (tri[0], tri[2], tri[1]) if nx * mx + ny * my + nz * mz < 0.0 else tri


def box_mesh(n: int = 2, size=(1.0, 1.0, 1.0)) -> Mesh:
    lx, ly, lz = size
    centre = (0.5 * lx, 0.5 * ly, 0.5 * lz)
    vertices: list[tuple[float, float, float]] = []
    lookup: dict[tuple[int, int, int], int] = {}

    def vid(i, j, k):
        if (i, j, k) not in lookup:
            lookup[(i, j, k)] = len(vertices)
            vertices.append((lx * i / n, ly * j / n, lz * k / n))
        return lookup[(i, j, k)]

    faces: list[tuple[int, int, int]] = []
    for axis in range(3):
        u, v = [a for a in range(3) if a != axis]
        for side in (0, 1):
            for a in range(n):
                for b in range(n):
                    quad = []
                    for du, dv in ((0, 0), (1, 0), (1, 1), (0, 1)):
                        index = [0, 0, 0]
                        index[axis] = side * n
                        index[u] = a + du
                        index[v] = b + dv
                        quad.append(vid(*index))
                    for tri in ((quad[0], quad[1], quad[2]), (quad[0], quad[2], quad[3])):
                        faces.append(_outward(vertices, tri, centre))
    return Mesh(vertices=tuple(vertices), faces=tuple(faces))


def sphere_mesh(radius: float = 0.5, rings: int = 12, segments: int = 24) -> Mesh:
    vertices = [(0.0, 0.0, radius)]
    for ring in range(1, rings):
        theta = math.pi * ring / rings
        for segment in range(segments):
            phi = 2.0 * math.pi * segment / segments
            vertices.append(
                (
                    radius * math.sin(theta) * math.cos(phi),
                    radius * math.sin(theta) * math.sin(phi),
                    radius * math.cos(theta),
                )
            )
    vertices.append((0.0, 0.0, -radius))
    south = len(vertices) - 1
    centre = (0.0, 0.0, 0.0)

    def index(ring, segment):
        return 1 + (ring - 1) * segments + (segment % segments)

    faces: list[tuple[int, int, int]] = []
    for segment in range(segments):
        faces.append(_outward(vertices, (0, index(1, segment), index(1, segment + 1)), centre))
    for ring in range(1, rings - 1):
        for segment in range(segments):
            a, b = index(ring, segment), index(ring, segment + 1)
            c, d = index(ring + 1, segment + 1), index(ring + 1, segment)
            faces.append(_outward(vertices, (a, b, c), centre))
            faces.append(_outward(vertices, (a, c, d), centre))
    for segment in range(segments):
        faces.append(
            _outward(
                vertices,
                (south, index(rings - 1, segment + 1), index(rings - 1, segment)),
                centre,
            )
        )
    return Mesh(vertices=tuple(vertices), faces=tuple(faces))


def poke(mesh: Mesh, vertex: int, target) -> Mesh:
    """Push one vertex to ``target``: the classic fold that pierces the far side."""
    vertices = list(mesh.vertices)
    vertices[vertex] = target
    return Mesh(
        vertices=tuple(vertices),
        faces=mesh.faces,
        units=mesh.units,
        coordinate_frame=mesh.coordinate_frame,
    )


def box_apex(mesh: Mesh) -> int:
    """The vertex at (0.5, 0, 0.5): an interior vertex of the box's ``y = 0`` face."""
    return next(
        index
        for index, point in enumerate(mesh.vertices)
        if abs(point[1]) < 1e-12 and abs(point[0] - 0.5) < 1e-12 and abs(point[2] - 0.5) < 1e-12
    )


def poked_box(n: int, depth: float) -> Mesh:
    box = box_mesh(n)
    return poke(box, box_apex(box), (0.5, depth, 0.5))


def drop_faces(mesh: Mesh, indices) -> Mesh:
    doomed = set(indices)
    return Mesh(
        vertices=mesh.vertices,
        faces=tuple(face for index, face in enumerate(mesh.faces) if index not in doomed),
        units=mesh.units,
        coordinate_frame=mesh.coordinate_frame,
    )


# --------------------------------------------------------------------------- plotting


def _edges(points, faces):
    """Triangle edge segments, for a flat projection."""
    segments = []
    for face in faces:
        corners = [points[index] for index in face]
        for k in range(3):
            segments.append((corners[k], corners[(k + 1) % 3]))
    return segments


def make_plot(before, after, intersecting, healed, rows) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection

    figure, axes = plt.subplots(1, 3, figsize=(15.0, 4.6))

    def draw(axis, mesh, title, highlight=()):
        points = [(point[0], point[1]) for point in mesh.vertices]
        base = _edges(points, mesh.faces)
        axis.add_collection(LineCollection(base, colors="#bbbbbb", linewidths=0.4))
        if highlight:
            hot = _edges(points, [mesh.faces[index] for index in highlight])
            axis.add_collection(LineCollection(hot, colors="#d62728", linewidths=1.4))
        axis.autoscale()
        axis.set_aspect("equal")
        axis.set_title(title, fontsize=9)
        axis.set_xlabel("x [m]", fontsize=8)
        axis.set_ylabel("y [m]", fontsize=8)

    draw(axes[0], before, f"input: {len(intersecting)} facet(s) in a fold")
    draw(axes[1], before, "the fold (red)", highlight=sorted(intersecting))
    draw(axes[2], healed, "healed: 0 penetrating pairs")

    names = [row["case"] for row in rows if row.get("volume_relative_change") is not None]
    changes = [row["volume_relative_change"] for row in rows if row.get("volume_relative_change") is not None]
    if names:
        figure.tight_layout()
        figure.savefig(OUT / "plot.png", dpi=160)
        figure.clf()
        figure, axes2 = plt.subplots(figsize=(7.0, 4.0))
        bars = axes2.bar(range(len(names)), [abs(value) * 100.0 for value in changes], color="#1f77b4")
        axes2.axhline(VOLUME_BUDGET * 100.0, color="#d62728", linestyle="--", label="1 % budget")
        axes2.set_xticks(range(len(names)))
        axes2.set_xticklabels(names, rotation=25, ha="right", fontsize=7)
        axes2.set_ylabel("|volume change| [%]")
        axes2.set_title("enclosed volume after the heal", fontsize=10)
        axes2.legend()
        for bar, value in zip(bars, changes):
            axes2.annotate(
                f"{value * 100.0:.3f}",
                (bar.get_x() + bar.get_width() / 2, bar.get_height()),
                ha="center",
                va="bottom",
                fontsize=7,
            )
        axes2.grid(axis="y", alpha=0.3)
        figure.tight_layout()
        figure.savefig(OUT / "plot_volume.png", dpi=160)
    plt.close("all")


# --------------------------------------------------------------------------- the runs


def run() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []
    claims: list[dict] = []
    failures: list[str] = []

    def claim(name, ok, detail):
        claims.append({"case": name, "claim": detail, "ok": bool(ok)})
        if not ok:
            failures.append(f"{name}: {detail}")

    # ---------------------------------------------------------------- box, poked through
    for n in (4, 6):
        box = box_mesh(n)
        control = measure_mesh(box)
        claim(f"box n={n} control", abs(control.signed_volume - 1.0) < 1e-12, f"volume {control.signed_volume!r} m3 vs 1")
        spiked = poke(box, box_apex(box), (0.5, 1.4, 0.5))
        pairs = find_self_intersections(spiked, limit=64)
        started = time.perf_counter()
        result = heal_self_intersections(spiked, reference_volume=1.0)
        elapsed = time.perf_counter() - started
        claim(f"box n={n} detected", len(pairs) > 0, f"{len(pairs)} penetrating pair(s) injected")
        claim(f"box n={n} healed", result.healed, result.reason)
        claim(
            f"box n={n} clean after",
            result.after.self_intersections == 0 and result.after.closed and result.after.manifold,
            f"self intersections {result.after.self_intersections}, closed {result.after.closed}, "
            f"manifold {result.after.manifold}",
        )
        claim(
            f"box n={n} volume",
            result.volume_relative_change is not None
            and result.volume_relative_change <= VOLUME_BUDGET,
            f"volume {result.volume_after!r} m3 vs 1, change {result.volume_relative_change!r}",
        )
        rows.append(
            {
                "case": f"box n={n}",
                "faces_in": control.faces,
                "pairs_before": len(pairs),
                "pairs_after": result.after.self_intersections,
                "rings": result.rings_used,
                "faces_removed": result.faces_removed,
                "faces_dropped": result.dropped_faces,
                "holes_capped": result.holes.filled if result.holes else 0,
                "volume_reference": 1.0,
                "volume_after": result.volume_after,
                "volume_relative_change": result.volume_relative_change,
                "seconds": elapsed,
                "healed": result.healed,
            }
        )

    # ---------------------------------------------------------------- sphere
    sphere = sphere_mesh()
    sphere_control = measure_mesh(sphere)
    exact = 4.0 / 3.0 * math.pi * 0.5**3
    claim(
        "sphere control",
        abs(sphere_control.signed_volume - exact) / exact < 0.03,
        f"tessellated {sphere_control.signed_volume!r} m3 vs exact {exact!r} "
        f"({(sphere_control.signed_volume - exact) / exact:.3%})",
    )
    apex = next(index for index, point in enumerate(sphere.vertices) if point[2] > 0.4 and point[0] > 0.1)
    spiked = poke(sphere, apex, (0.0, 0.0, -0.6))
    pairs = find_self_intersections(spiked, limit=64)
    started = time.perf_counter()
    result = heal_self_intersections(spiked, reference_volume=sphere_control.signed_volume)
    elapsed = time.perf_counter() - started
    claim("sphere detected", len(pairs) > 0, f"{len(pairs)} penetrating pair(s) injected")
    claim("sphere healed", result.healed, result.reason)
    claim(
        "sphere clean after",
        result.after.self_intersections == 0 and result.after.closed,
        f"self intersections {result.after.self_intersections}, closed {result.after.closed}",
    )
    claim(
        "sphere volume",
        result.volume_relative_change is not None and result.volume_relative_change <= VOLUME_BUDGET,
        f"change {result.volume_relative_change!r} of the tessellated volume",
    )
    rows.append(
        {
            "case": "sphere r=0.5",
            "faces_in": sphere_control.faces,
            "pairs_before": len(pairs),
            "pairs_after": result.after.self_intersections,
            "rings": result.rings_used,
            "faces_removed": result.faces_removed,
            "faces_dropped": result.dropped_faces,
            "holes_capped": result.holes.filled if result.holes else 0,
            "volume_reference": sphere_control.signed_volume,
            "volume_exact": exact,
            "volume_after": result.volume_after,
            "volume_relative_change": result.volume_relative_change,
            "seconds": elapsed,
            "healed": result.healed,
        }
    )

    # ---------------------------------------------------------------- parametric hull
    hull = build_mesh(HullParameters(length=LENGTH, beam=BEAM, draft=DRAFT))[0]
    hull_control = measure_mesh(hull)
    hull_analytic = analytic_volume(LENGTH, BEAM, DRAFT)
    claim("hull control", hull_control.closed, f"{hull_control.faces} facets, closed {hull_control.closed}")
    keel = min(range(len(hull.vertices)), key=lambda index: hull.vertices[index][2])
    spiked = poke(hull, keel, (hull.vertices[keel][0], hull.vertices[keel][1], hull.vertices[keel][2] + 1.6))
    pairs = find_self_intersections(spiked, limit=64)
    started = time.perf_counter()
    result = heal_self_intersections(spiked, reference_volume=hull_control.signed_volume)
    elapsed = time.perf_counter() - started
    claim("hull detected", len(pairs) > 0, f"{len(pairs)} penetrating pair(s) injected at the keel")
    claim("hull healed", result.healed, result.reason)
    claim(
        "hull clean after",
        result.after.self_intersections == 0
        and result.after.closed
        and result.after.manifold
        and result.after.consistently_wound,
        f"self intersections {result.after.self_intersections}, closed {result.after.closed}, "
        f"manifold {result.after.manifold}, wound {result.after.consistently_wound}",
    )
    claim(
        "hull volume",
        result.volume_relative_change is not None and result.volume_relative_change <= VOLUME_BUDGET,
        f"change {result.volume_relative_change!r}",
    )
    rows.append(
        {
            "case": "hull 41x13",
            "faces_in": hull_control.faces,
            "pairs_before": len(pairs),
            "pairs_after": result.after.self_intersections,
            "rings": result.rings_used,
            "faces_removed": result.faces_removed,
            "faces_dropped": result.dropped_faces,
            "holes_capped": result.holes.filled if result.holes else 0,
            "volume_reference": hull_control.signed_volume,
            "volume_analytic": hull_analytic,
            "volume_after": result.volume_after,
            "volume_relative_change": result.volume_relative_change,
            "volume_error_vs_analytic": (result.volume_after - hull_analytic) / hull_analytic,
            "seconds": elapsed,
            "healed": result.healed,
        }
    )
    healed_hull = result

    # ---------------------------------------------------------------- idempotence
    again = heal_self_intersections(healed_hull.mesh, reference_volume=hull_control.signed_volume)
    claim(
        "idempotent",
        not again.healed and again.after.self_intersections == 0 and len(again.mesh.faces) == len(healed_hull.mesh.faces),
        f"second pass healed={again.healed}, faces {len(again.mesh.faces)} vs {len(healed_hull.mesh.faces)}",
    )

    # ---------------------------------------------------------------- the pipeline flag
    spiked_hull = poke(
        hull, keel, (hull.vertices[keel][0], hull.vertices[keel][1], hull.vertices[keel][2] + 0.6)
    )
    plain = repair_mesh(spiked_hull)
    with_flag = repair_mesh(
        spiked_hull,
        heal_self_intersections=True,
        reference_volume=hull_control.signed_volume,
    )
    claim(
        "repair_mesh refuses a fold by default",
        not plain.accepted and plain.after.self_intersections > 0,
        f"accepted {plain.accepted}, self intersections {plain.after.self_intersections}",
    )
    claim(
        "repair_mesh heals with the flag",
        with_flag.accepted and with_flag.after.self_intersections == 0,
        f"accepted {with_flag.accepted}, self intersections {with_flag.after.self_intersections}, "
        f"steps {[step.name for step in with_flag.steps]}",
    )
    rows.append(
        {
            "case": "repair_mesh flag",
            "pairs_before": plain.after.self_intersections,
            "pairs_after": with_flag.after.self_intersections,
            "accepted_without_flag": plain.accepted,
            "accepted_with_flag": with_flag.accepted,
        }
    )

    # ---------------------------------------------------------------- refusals
    open_mesh = drop_faces(box_mesh(2), [0, 1])
    refused = heal_self_intersections(poke(open_mesh, box_apex(box_mesh(2)), (0.5, 1.4, 0.5)))
    claim("refuses an open input", not refused.healed, refused.reason)

    tight = heal_self_intersections(poked_box(4, 1.4), max_hole_edges=4)
    claim(
        "a 4-edge cap limit falls back to the fan, then the volume guard refuses",
        not tight.healed,
        f"healed={tight.healed}: {tight.reason}",
    )

    coarse = heal_self_intersections(poked_box(2, 1.4), reference_volume=1.0)
    claim("refuses a non-local cut", not coarse.healed, coarse.reason)

    dent = poke(box_mesh(2), box_apex(box_mesh(2)), (0.5, 0.6, 0.5))
    untouched = heal_self_intersections(dent)
    claim(
        "a dent is not a fold",
        not untouched.healed and untouched.mesh is dent,
        untouched.reason,
    )

    band = excavation_faces(poked_box(4, 1.4), rings=1)
    claim("excavation preview", len(band) > 0, f"{len(band)} facet(s) in the 1-ring band")
    claim("clean input is a no-op", len(excavation_faces(hull)) == 0, "no penetration, empty band")

    rows.append(
        {
            "case": "refusals",
            "open_input": refused.reason,
            "cap_limit": tight.reason,
            "non_local": coarse.reason,
            "dent": untouched.reason,
        }
    )

    # ---------------------------------------------------------------- artifacts
    box_case = poked_box(4, 1.4)
    box_pairs = find_self_intersections(box_case, limit=64)
    box_healed = heal_self_intersections(box_case, reference_volume=1.0)
    fold = {index for pair in box_pairs for index in pair}
    make_plot(box_case, box_case, fold, box_healed.mesh, rows)

    payload = {
        "rows": rows,
        "claims": claims,
        "failures": failures,
        "volume_budget": VOLUME_BUDGET,
        "figures": ["plot.png", "plot_volume.png"],
        "python": sys.version.split()[0],
    }
    (OUT / "results.json").write_text(json.dumps(payload, indent=2, sort_keys=True))

    lines = ["# self-intersection heal: measured claims", ""]
    lines.append("| case | claim | ok |")
    lines.append("|---|---|---|")
    for row in claims:
        lines.append(f"| {row['case']} | {row['claim']} | {'yes' if row['ok'] else '**NO**'} |")
    lines.append("")
    lines.append("| case | facets | pairs before | pairs after | rings | facets cut | holes | volume change |")
    lines.append("|---|---|---|---|---|---|---|---|")
    for row in rows:
        change = row.get("volume_relative_change")
        lines.append(
            f"| {row['case']} | {row.get('faces_in', '')} | {row.get('pairs_before', '')} | "
            f"{row.get('pairs_after', '')} | {row.get('rings', '')} | {row.get('faces_removed', '')} | "
            f"{row.get('holes_capped', '')} | {'—' if change is None else f'{change:.3%}'} |"
        )
    (OUT / "report.md").write_text("\n".join(lines) + "\n")

    print("\n".join(lines))
    if failures:
        print(f"\n{len(failures)} claim(s) FAILED")
        for failure in failures:
            print(f"  - {failure}")
        return 1
    print(f"\nall {len(claims)} claims hold")
    return 0


if __name__ == "__main__":
    raise SystemExit(run())
