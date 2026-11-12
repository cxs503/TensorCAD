"""Measure the repair path of ``tensorcad.mesh_repair`` and fail loudly if a claim slips.

Writes ``artifacts/mesh_repair/results.json``, ``plot.png`` and ``report.md``; exits non-zero
when any acceptance claim stops holding. Every number in ``docs/mesh_repair.md`` is read from
that artifact.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

from tensorcad.geometry import Mesh
from tensorcad.hull import HullParameters, analytic_volume, generate_hull
from tensorcad.mesh_repair import (
    fill_holes,
    manifoldize,
    measure_mesh,
    orient_consistently,
    repair_mesh,
)

HERE = Path(__file__).resolve().parent
OUT = HERE.parents[1] / "artifacts" / "mesh_repair"

LENGTH, BEAM, DRAFT = 10.0, 2.0, 0.8


def _outward(vertices, tri, centre):
    ax, ay, az = vertices[tri[0]]
    bx, by, bz = vertices[tri[1]]
    cx, cy, cz = vertices[tri[2]]
    ux, uy, uz = bx - ax, by - ay, bz - az
    vx, vy, vz = cx - ax, cy - ay, cz - az
    nx, ny, nz = uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx
    mx = (ax + bx + cx) / 3.0 - centre[0]
    my = (ay + by + cy) / 3.0 - centre[1]
    mz = (az + bz + cz) / 3.0 - centre[2]
    return (tri[0], tri[2], tri[1]) if nx * mx + ny * my + nz * mz < 0.0 else tri


def box_mesh(n=2, size=(1.0, 1.0, 1.0)):
    lx, ly, lz = size
    centre = (0.5 * lx, 0.5 * ly, 0.5 * lz)
    vertices, lookup = [], {}

    def vid(i, j, k):
        if (i, j, k) not in lookup:
            lookup[(i, j, k)] = len(vertices)
            vertices.append((lx * i / n, ly * j / n, lz * k / n))
        return lookup[(i, j, k)]

    faces = []
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


def sphere_mesh(radius=0.5, rings=12, segments=24):
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

    faces = []
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


def drop_faces(mesh, indices):
    doomed = set(indices)
    return Mesh(
        vertices=mesh.vertices,
        faces=tuple(face for index, face in enumerate(mesh.faces) if index not in doomed),
    )


def run():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    failures = []

    def claim(name, ok, detail):
        rows.append({"case": name, "claim": detail, "ok": bool(ok)})
        if not ok:
            failures.append(f"{name}: {detail}")

    # ---------------------------------------------------------------- box controls
    box = box_mesh(2)
    box_volume = measure_mesh(box).signed_volume
    claim("box control", abs(box_volume - 1.0) < 1e-12, f"closed box volume {box_volume!r} m3 vs 1")

    for label, doomed in (
        ("box: 1 facet removed", [0]),
        ("box: quad removed (2 facets)", [0, 1]),
        ("box: 2 holes", [0, 54]),
    ):
        holed = drop_faces(box, doomed)
        result = repair_mesh(holed, reference_volume=1.0)
        claim(
            label,
            result.accepted
            and result.after.closed
            and result.after.manifold
            and abs(result.after.signed_volume - 1.0) < 1e-12,
            f"accepted={result.accepted} faces {result.before.faces}->{result.after.faces} "
            f"volume {result.after.signed_volume:.12f} m3 (rel change "
            f"{result.volume_relative_change:.3e})",
        )
        rows[-1]["volume_relative_change"] = result.volume_relative_change
        rows[-1]["filled_holes"] = result.holes.filled if result.holes else 0

    # ---------------------------------------------------------------- curved surface
    sphere = sphere_mesh()
    analytic = 4.0 / 3.0 * math.pi * 0.5**3
    tessellated = measure_mesh(sphere).signed_volume
    fill = fill_holes(drop_faces(sphere, [3, 4]))
    repaired = measure_mesh(fill.mesh)
    claim(
        "sphere: 2 facets removed",
        repaired.closed and repaired.self_intersections == 0,
        f"closed={repaired.closed} volume {repaired.signed_volume:.6f} m3 vs tessellated "
        f"{tessellated:.6f} ({(repaired.signed_volume / tessellated - 1) * 100:+.3f} %) vs "
        f"analytic {analytic:.6f} ({(repaired.signed_volume / analytic - 1) * 100:+.2f} %)",
    )
    rows[-1]["planar_deviation_m"] = fill.holes[0]["planar_deviation_m"]

    # ---------------------------------------------------------------- hull
    parameters = HullParameters(length=LENGTH, beam=BEAM, draft=DRAFT, nx=41, nz=13)
    hull = generate_hull(parameters).mesh
    reference = analytic_volume(LENGTH, BEAM, DRAFT)
    hull_before = measure_mesh(hull)
    holed = drop_faces(hull, [7, 8, 9])
    result = repair_mesh(holed, reference_volume=reference)
    claim(
        "hull 41x13: 3 facets removed",
        result.accepted
        and result.after.closed
        and result.after.self_intersections == 0
        and result.volume_relative_change < 0.01,
        f"faces {result.before.faces}->{result.after.faces}, open edges "
        f"{result.before.open_edges}->{result.after.open_edges}, volume "
        f"{result.after.signed_volume:.6f} m3 vs analytic {reference:.6f} "
        f"({(result.after.signed_volume / reference - 1) * 100:+.4f} %), volume move "
        f"{result.volume_relative_change:.3e}",
    )
    hull_repaired = result.mesh

    # ---------------------------------------------------------------- winding
    inward = Mesh(vertices=box.vertices, faces=tuple((f[0], f[2], f[1]) for f in box.faces))
    fixed, step = orient_consistently(inward)
    claim(
        "box: inside-out shell",
        measure_mesh(fixed).signed_volume > 0 and step.changed == len(box.faces),
        f"volume {measure_mesh(inward).signed_volume!r} -> "
        f"{measure_mesh(fixed).signed_volume!r} m3, {step.changed} facets flipped",
    )

    one_bad = Mesh(
        vertices=box.vertices,
        faces=tuple([(box.faces[0][0], box.faces[0][2], box.faces[0][1])] + list(box.faces[1:])),
    )
    fixed, step = orient_consistently(one_bad)
    claim(
        "box: one backwards facet",
        measure_mesh(fixed).winding_conflicts == 0 and step.changed == 1,
        f"conflicts {measure_mesh(one_bad).winding_conflicts} -> "
        f"{measure_mesh(fixed).winding_conflicts}, {step.changed} facet flipped",
    )

    # ---------------------------------------------------------------- non-manifold
    duplicated = Mesh(vertices=box.vertices, faces=box.faces + (box.faces[0],))
    result = repair_mesh(duplicated, reference_volume=1.0)
    claim(
        "box: duplicated facet",
        result.accepted and result.after.nonmanifold_edges == 0
        and abs(result.after.signed_volume - 1.0) < 1e-12,
        f"non-manifold edges {measure_mesh(duplicated).nonmanifold_edges} -> "
        f"{result.after.nonmanifold_edges}, volume {result.after.signed_volume!r} m3",
    )

    # ---------------------------------------------------------------- refusals
    vertices = list(sphere.vertices)
    vertices[0] = (0.0, 0.0, -0.9)
    punctured = Mesh(vertices=tuple(vertices), faces=sphere.faces)
    refused = repair_mesh(punctured)
    claim(
        "sphere: pole pulled through (self-intersecting)",
        (not refused.accepted)
        and "self-intersecting" in refused.reason
        and refused.mesh.faces == punctured.faces,
        f"accepted={refused.accepted}; {refused.reason[:120]}",
    )

    capped = repair_mesh(drop_faces(box_mesh(6), range(12)), reference_volume=1.0, max_hole_edges=6)
    claim(
        "box: 12-facet hole against a 6-edge cap",
        (not capped.accepted) and capped.holes is not None and len(capped.holes.refused) == 1,
        f"accepted={capped.accepted}; {capped.reason[:120]}",
    )

    payload = {
        "fixtures": {
            "box_faces": len(box.faces),
            "box_volume_m3": box_volume,
            "sphere_faces": len(sphere.faces),
            "sphere_tessellated_volume_m3": tessellated,
            "sphere_analytic_volume_m3": analytic,
            "hull_faces": len(hull.faces),
            "hull_surface_volume_m3": hull_before.signed_volume,
            "hull_analytic_volume_m3": reference,
        },
        "claims": rows,
        "failures": failures,
    }
    (OUT / "results.json").write_text(json.dumps(payload, indent=2) + "\n")

    _write_report(payload)
    _write_figure(holed, hull_repaired, rows)
    if failures:
        print("FAILED CLAIMS:")
        for item in failures:
            print(" -", item)
        return 1
    print(f"{len(rows)} claims held; artifacts in {OUT}")
    return 0


def _write_report(payload):
    lines = ["# Mesh repair: measured claims", ""]
    lines.append("| case | claim | held |")
    lines.append("|---|---|---|")
    for row in payload["claims"]:
        lines.append(f"| {row['case']} | {row['claim']} | {'yes' if row['ok'] else 'NO'} |")
    lines.append("")
    lines.append("Fixtures:")
    for key, value in payload["fixtures"].items():
        lines.append(f"* `{key}` = {value}")
    (OUT / "report.md").write_text("\n".join(lines) + "\n")


def _write_figure(holed, repaired, rows):
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as error:  # pragma: no cover - plotting is optional
        print("plot skipped:", error)
        return
    fig = plt.figure(figsize=(11, 4.6))
    for panel, (mesh, title) in enumerate(
        ((holed, "before: hull with 3 facets removed"), (repaired, "after: repair_mesh"))
    ):
        axis = fig.add_subplot(1, 3, panel + 1, projection="3d")
        x = [point[0] for point in mesh.vertices]
        y = [point[1] for point in mesh.vertices]
        z = [point[2] for point in mesh.vertices]
        axis.plot_trisurf(x, y, z, triangles=mesh.faces, linewidth=0.15, color="#9ecae1")
        axis.set_title(title, fontsize=9)
        axis.set_box_aspect((2.0, 1.0, 0.5))
    axis = fig.add_subplot(1, 3, 3)
    labels = [row["case"].replace(": ", "\n") for row in rows if "volume_relative_change" in row]
    values = [max(row["volume_relative_change"], 1e-18) for row in rows if "volume_relative_change" in row]
    axis.bar(range(len(values)), values, color="#3182bd")
    axis.set_yscale("log")
    axis.axhline(1e-2, color="black", linestyle="--", linewidth=1, label="1 % acceptance bar")
    axis.set_xticks(range(len(values)))
    axis.set_xticklabels(labels, fontsize=6)
    axis.set_ylabel("|volume change| / reference")
    axis.legend(fontsize=7)
    axis.set_title("volume moved by the fill", fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / "plot.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    sys.exit(run())