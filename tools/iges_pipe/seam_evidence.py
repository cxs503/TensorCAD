"""Evidence run for the seam stage (GAP-CAD-T1-IGES-PIPE, seam/closure sub-gap).

Reproduce from the repository root::

    python tools/iges_pipe/seam_evidence.py --out artifacts/gap_t1_iges_seam

Every number written to ``evidence.json`` is read back from the run that just happened:
the publisher files are imported through :func:`tensorcad.geomquality.trim_mesh.import_cad`,
then swept over the seam budget and measured (free edges, non-manifold edges, splits,
enclosed volume change, where the residual open boundary sits).  The figure is rendered
from that same JSON, so the plot cannot drift from the table.
"""
from __future__ import annotations

import argparse
import collections
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from tensorcad.geomquality.trim_mesh import import_cad, weld_vertices  # noqa: E402
from tensorcad.geomquality.stitch import (  # noqa: E402
    edge_use_counts,
    free_edges,
    seam_report,
    stitch_t_junctions,
)

ASSETS = Path("/srv/tensor/agent/assets/geometry/hulls")
FILES = [
    ("KCS", ASSETS / "KCS" / "KCS_hull_SVA.igs"),
    ("JBC_boss", ASSETS / "JBC" / "jbc-boss.igs"),
]
BUDGETS = (0.0, 5e-3, 1.5e-2, 5e-2)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def signed_volume(vertices, faces) -> float:
    total = 0.0
    for a, b, c in faces:
        ax, ay, az = vertices[a]
        bx, by, bz = vertices[b]
        cx, cy, cz = vertices[c]
        total += (ax * (by * cz - bz * cy) - ay * (bx * cz - bz * cx)
                  + az * (bx * cy - by * cx)) / 6.0
    return total


def residual_classification(vertices, faces) -> dict[str, int]:
    """Where the open boundary still is after stitching: centreplane / top / elsewhere."""
    open_edges = free_edges(faces)
    if not open_edges:
        return {"centreplane_y0": 0, "top_zmax": 0, "other": 0}
    zmax = max(v[2] for v in vertices)
    counts = collections.Counter()
    for index in {i for edge in open_edges for i in edge}:
        x, y, z = vertices[index]
        if abs(y) < 1e-6:
            counts["centreplane_y0"] += 1
        elif z > zmax - 1e-6:
            counts["top_zmax"] += 1
        else:
            counts["other"] += 1
    return dict(counts)


def measure(name: str, path: Path, budgets=BUDGETS) -> dict:
    base = import_cad(path, weld_tolerance=1e-3, repair=False)
    vertices = list(base.mesh.vertices)
    faces = list(base.mesh.faces)
    units = base.report["source"]["metres_per_unit"]
    length = max(max(v[0] for v in vertices) - min(v[0] for v in vertices),
                 max(v[2] for v in vertices) - min(v[2] for v in vertices))
    entry = {
        "ship": name,
        "file": str(path),
        "sha256": sha256(path),
        "bytes": path.stat().st_size,
        "metres_per_unit": units,
        "length_m": length,
        "triangles_tessellated": len(faces),
        "entity_census": base.report["import"]["entity_census"],
        "seam_deviation_m": seam_report(vertices, faces, tolerance=1e-3)[
            "distance_to_nearest_free_chord_quantiles"],
        "sweep": [],
    }
    for budget in budgets:
        current_vertices, current_faces = vertices, faces
        record = None
        if budget > 0.0:
            current_vertices, current_faces, record = stitch_t_junctions(
                vertices, faces, tolerance=budget)
        current_vertices, current_faces, _ = weld_vertices(
            current_vertices, current_faces, max(budget, 1e-3))
        uses = edge_use_counts(current_faces)
        entry["sweep"].append({
            "seam_budget_m": budget,
            "budget_over_length": budget / length if length else None,
            "free_edges": len(free_edges(current_faces)),
            "non_manifold_edges": sum(1 for count in uses.values() if count > 2),
            "splits": (record or {}).get("splits_total", 0),
            "triangles": len(current_faces),
            "volume_m3": signed_volume(current_vertices, current_faces),
            "residual_open_boundary": residual_classification(current_vertices, current_faces),
        })
    baseline = entry["sweep"][0]
    for row in entry["sweep"]:
        row["free_edges_closed_percent"] = (
            100.0 * (baseline["free_edges"] - row["free_edges"]) / baseline["free_edges"]
            if baseline["free_edges"] else 0.0)
        row["volume_change_m3"] = row["volume_m3"] - baseline["volume_m3"]
    return entry


def figure(evidence: dict, out_png: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ships = list(evidence["ships"].keys())
    figure, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for ship in ships:
        rows = evidence["ships"][ship]["sweep"]
        x = [row["budget_over_length"] * 1e5 for row in rows]  # 1e-5 L units
        axes[0].plot(x, [row["free_edges"] for row in rows], marker="o", label=ship)
    axes[0].set_xlabel("seam budget  [1e-5 · L]")
    axes[0].set_ylabel("free edges after stitching")
    axes[0].set_title("Stitching closes seams but plateaus")
    axes[0].grid(alpha=0.3)
    axes[0].legend()

    ship = ships[0]
    rows = evidence["ships"][ship]["sweep"]
    labels = [f"{row['seam_budget_m'] * 1e3:.0f} mm" for row in rows]
    keys = ["centreplane_y0", "top_zmax", "other"]
    bottom = [0.0] * len(rows)
    for key in keys:
        values = [row["residual_open_boundary"].get(key, 0) for row in rows]
        axes[1].bar(labels, values, bottom=bottom, label=key)
        bottom = [b + v for b, v in zip(bottom, values)]
    axes[1].set_title(f"{ship}: residual open boundary is the model's real openings")
    axes[1].set_ylabel("free-edge endpoints")
    axes[1].legend(fontsize=8)
    axes[1].grid(alpha=0.3, axis="y")
    figure.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(out_png, dpi=140)
    plt.close(figure)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="artifacts/gap_t1_iges_seam")
    args = parser.parse_args(argv)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    evidence: dict = {"tool": "tools/iges_pipe/seam_evidence.py", "ships": {}, "skipped": {}}
    for ship, path in FILES:
        if not path.exists():
            evidence["skipped"][ship] = f"not staged: {path}"
            continue
        evidence["ships"][ship] = measure(ship, path)

    (out / "evidence.json").write_text(json.dumps(evidence, indent=2) + "\n", encoding="utf-8")
    if evidence["ships"]:
        figure(evidence, out / "seam_stitch.png")
    for ship, entry in evidence["ships"].items():
        print(f"{ship}: {entry['triangles_tessellated']} triangles, seam deviation "
              f"q50={entry['seam_deviation_m']['q50'] * 1e3:.2f} mm")
        for row in entry["sweep"]:
            print(f"  budget {row['seam_budget_m'] * 1e3:5.1f} mm -> free {row['free_edges']:4d} "
                  f"({row['free_edges_closed_percent']:5.1f} % closed), splits {row['splits']:4d}, "
                  f"non-manifold {row['non_manifold_edges']:3d}, "
                  f"dV {row['volume_change_m3']:+.1f} m3, residual {row['residual_open_boundary']}")
    print(f"wrote {out / 'evidence.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())