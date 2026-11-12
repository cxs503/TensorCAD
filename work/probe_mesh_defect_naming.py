"""Measure mesh-only signals on the geomquality corpus (3 hulls x 14 injections).

Read-only probe behind ``tensorcad.geomquality.mesh_defects``: it exists so the
thresholds in that module were read off measurements rather than guessed.

    python work/probe_mesh_defect_naming.py            # table to stdout + JSON
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from tensorcad.geomquality import baselines, defects  # noqa: E402


def _norm(v):
    return math.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)


def _sub(verts, a, b):
    return (verts[a][0] - verts[b][0], verts[a][1] - verts[b][1], verts[a][2] - verts[b][2])


def _cross(u, v):
    return (
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    )


def signals(mesh, length: float) -> dict[str, float]:
    verts = mesh.vertices
    faces = [tuple(f) for f in mesh.faces]
    edge_faces: dict[tuple[int, int], list[int]] = {}
    for index, face in enumerate(faces):
        for i in range(3):
            a, b = face[i], face[(i + 1) % 3]
            key = (a, b) if a < b else (b, a)
            edge_faces.setdefault(key, []).append(index)

    lengths = sorted(_norm(_sub(verts, a, b)) for (a, b) in edge_faces)
    median_edge = lengths[len(lengths) // 2] if lengths else 0.0
    areas = []
    rel_altitudes = []
    min_angles = []
    for face in faces:
        e1 = _sub(verts, face[1], face[0])
        e2 = _sub(verts, face[2], face[0])
        e3 = _sub(verts, face[2], face[1])
        area = 0.5 * _norm(_cross(e1, e2))
        areas.append(area)
        sides = sorted([_norm(e1), _norm(e2), _norm(e3)])
        longest = sides[2]
        rel_altitudes.append((2.0 * area / longest) / longest if longest > 0 else 0.0)
        if sides[0] > 0.0 and sides[1] > 0.0:
            cosine = (sides[1] ** 2 + sides[2] ** 2 - sides[0] ** 2) / (2.0 * sides[1] * sides[2])
            min_angles.append(math.degrees(math.acos(max(-1.0, min(1.0, cosine)))))
        else:
            min_angles.append(0.0)
    sorted_areas = sorted(areas)
    median_area = sorted_areas[len(sorted_areas) // 2] if sorted_areas else 0.0

    normals = []
    for face in faces:
        e1 = _sub(verts, face[1], face[0])
        e2 = _sub(verts, face[2], face[0])
        cr = _cross(e1, e2)
        n = _norm(cr)
        normals.append((cr[0] / n, cr[1] / n, cr[2] / n) if n > 0 else (0.0, 0.0, 0.0))

    dihedrals = []
    for adjacent in edge_faces.values():
        if len(adjacent) == 2:
            a, b = normals[adjacent[0]], normals[adjacent[1]]
            cosine = max(-1.0, min(1.0, a[0] * b[0] + a[1] * b[1] + a[2] * b[2]))
            dihedrals.append(math.degrees(math.acos(cosine)))
    dihedrals.sort()

    def pct(values, fraction):
        if not values:
            return 0.0
        return values[min(len(values) - 1, max(0, int(round(fraction * (len(values) - 1)))))]

    return {
        "faces": len(faces),
        "open_edges": sum(1 for fs in edge_faces.values() if len(fs) != 2),
        "min_edge_rel_median": (lengths[0] / median_edge) if median_edge else 0.0,
        "edges_lt_weld_tol": sum(1 for x in lengths if x < 1e-6 * length),
        "edges_lt_1e-3_median": sum(1 for x in lengths if x < 1e-3 * median_edge),
        "median_edge_rel_length": median_edge / length,
        "min_rel_altitude": min(rel_altitudes) if rel_altitudes else 0.0,
        "faces_rel_alt_lt_1e-3": sum(1 for x in rel_altitudes if x < 1e-3),
        "faces_rel_alt_lt_1e-6": sum(1 for x in rel_altitudes if x < 1e-6),
        "faces_area_lt_1e-6_median": sum(1 for x in areas if x < 1e-6 * median_area),
        "min_angle_deg": min(min_angles) if min_angles else 0.0,
        "faces_angle_lt_1": sum(1 for x in min_angles if x < 1.0),
        "faces_angle_lt_5": sum(1 for x in min_angles if x < 5.0),
        "faces_angle_lt_10": sum(1 for x in min_angles if x < 10.0),
        "dihedral_p50": pct(dihedrals, 0.5),
        "dihedral_p95": pct(dihedrals, 0.95),
        "dihedral_p999": pct(dihedrals, 0.999),
        "dihedral_max": dihedrals[-1] if dihedrals else 0.0,
        "dihedral_gt20": sum(1 for x in dihedrals if x > 20),
        "dihedral_gt40": sum(1 for x in dihedrals if x > 40),
        "dihedral_gt80": sum(1 for x in dihedrals if x > 80),
    }


def main() -> int:
    catalogue = baselines.specimens()
    out = {}
    for name in sorted(catalogue):
        specimen = catalogue[name]()
        for case in [None, *defects.DEFECT_CASES]:
            brep = specimen.brep if case is None else case.inject(specimen.brep)
            key = f"{name}|{case.code if case else 'CLEAN'}"
            out[key] = signals(brep.tessellate(16, 8), specimen.length)
            print(key, json.dumps(out[key]))
    Path("/tmp/probe_mesh_defect_naming.json").write_text(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())