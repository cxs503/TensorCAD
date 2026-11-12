"""Second probe: the signals the first sweep did not carry (signed dihedral, crease
chain magnitudes, curvature-jump rate). Focused on the four families the mesh route
could not name.

    python work/probe_mesh_defect_naming2.py
"""
from __future__ import annotations

import json
import math
import sys
from collections import defaultdict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from tensorcad.geomquality import baselines, defects  # noqa: E402

FAMILIES = ("CLEAN", "SLIVER", "MICRO_EDGE", "TANGENT_JUMP", "CURVATURE_JUMP", "ROUGHNESS")


def _sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def _norm(a):
    return math.sqrt(_dot(a, a))


def analyse(mesh, length: float) -> dict[str, object]:
    verts = list(mesh.vertices)
    faces = [tuple(f) for f in mesh.faces]
    edge_faces: dict[tuple[int, int], list[int]] = defaultdict(list)
    for index, face in enumerate(faces):
        for i in range(3):
            a, b = face[i], face[(i + 1) % 3]
            edge_faces[(a, b) if a < b else (b, a)].append(index)

    normals = []
    for face in faces:
        n = _cross(_sub(verts[face[1]], verts[face[0]]), _sub(verts[face[2]], verts[face[0]]))
        m = _norm(n)
        normals.append((n[0] / m, n[1] / m, n[2] / m) if m > 0 else (0.0, 0.0, 0.0))

    signed = []  # convex positive
    for (a, b), adjacent in edge_faces.items():
        if len(adjacent) != 2:
            continue
        f1, f2 = adjacent
        n1, n2 = normals[f1], normals[f2]
        cosine = max(-1.0, min(1.0, _dot(n1, n2)))
        magnitude = math.degrees(math.acos(cosine))
        # the vertex of f2 that is not on the shared edge decides the side
        other = next(v for v in faces[f2] if v not in (a, b))
        side = _dot(n1, _sub(verts[other], verts[a]))
        signed.append(math.copysign(magnitude, side))
    signed.sort()
    convex = [x for x in signed if x > 0.0]
    concave = [x for x in signed if x < 0.0]

    def pct(values, fraction):
        if not values:
            return 0.0
        return values[min(len(values) - 1, max(0, int(round(fraction * (len(values) - 1)))))]

    magnitudes = sorted(abs(x) for x in signed)
    # crease chains: connected components of edges above a threshold
    def chains(threshold: float) -> dict[str, float]:
        hot = {key for key, adj in edge_faces.items() if len(adj) == 2 and abs(_angle(normals, adj)) > threshold}
        seen = set()
        sizes = []
        for start in hot:
            if start in seen:
                continue
            stack = [start]
            seen.add(start)
            size = 0
            while stack:
                key = stack.pop()
                size += 1
                for vertex in key:
                    for other_key, adj in edges_by_vertex[vertex]:
                        if other_key in hot and other_key not in seen:
                            seen.add(other_key)
                            stack.append(other_key)
            sizes.append(size)
        sizes.sort(reverse=True)
        return {
            f"chains_gt{threshold:g}": len(sizes),
            f"chain_max_gt{threshold:g}": sizes[0] if sizes else 0,
            f"chain_sum_gt{threshold:g}": sum(sizes),
        }

    edges_by_vertex: dict[int, list] = defaultdict(list)
    for key in edge_faces:
        for vertex in key:
            edges_by_vertex[vertex].append((key, edge_faces[key]))

    def _angle(norm, adj):
        cosine = max(-1.0, min(1.0, _dot(norm[adj[0]], norm[adj[1]])))
        return math.degrees(math.acos(cosine))

    result = {
        "interior_edges": len(signed),
        "convex_edges": len(convex),
        "concave_edges": len(concave),
        "concave_share": len(concave) / len(signed) if signed else 0.0,
        "signed_p99": pct(signed, 0.99),
        "signed_p01": pct(signed, 0.01),
        "median_abs": pct(magnitudes, 0.5),
        "p95_abs": pct(magnitudes, 0.95),
        "gt0p5": sum(1 for x in magnitudes if x > 0.5),
        "gt2": sum(1 for x in magnitudes if x > 2.0),
        "gt5": sum(1 for x in magnitudes if x > 5.0),
        "gt15": sum(1 for x in magnitudes if x > 15.0),
    }
    for threshold in (5.0, 15.0, 30.0):
        result.update(chains(threshold))
    return result


def main() -> int:
    catalogue = baselines.specimens()
    cases = {case.code: case for case in defects.DEFECT_CASES}
    out = {}
    for name in sorted(catalogue):
        specimen = catalogue[name]()
        for family in FAMILIES:
            brep = specimen.brep if family == "CLEAN" else cases[family].inject(specimen.brep)
            key = f"{name}|{family}"
            out[key] = analyse(brep.tessellate(16, 8), specimen.length)
            print(key, json.dumps(out[key]))
    Path("/tmp/probe_mesh_defect_naming2.json").write_text(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())