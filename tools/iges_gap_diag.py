"""Open-boundary diagnosis for the IGES/STEP → watertight mesh import pipeline.

Why this tool exists
--------------------
``import_iges`` reports *how many* free edges survive repair.  It does not report *why*
they survive, so a residual open boundary can only be argued about, not fixed.  This tool
classifies every surviving free edge by the repair operation that would close it:

``coincident``   both endpoints have a coincident twin vertex (``<= eps``) elsewhere in the
                 mesh -- a weld at ``eps`` closes it, so a weld that ran at a tighter
                 tolerance is the defect;
``t_junction``   one endpoint lies on the interior of another free edge, which the seam
                 stage can split;
``gap``          the edge's closest free edge is ``d`` away (``d > eps``) -- the file is
                 genuinely open there and no local repair may close it.

The classification is what decides whether a residual free edge is a bug in the repair
chain or a property of the input, and it is the number quoted in
``docs/iges_pipe_seam_study.md``.

Usage::

    python tools/iges_gap_diag.py work/assets/kcs.igs --seam-tolerance 0.2
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from tensorcad.geomquality.stitch import free_edges, project_on_segment  # noqa: E402
from tensorcad.geomquality.trim_mesh import import_iges  # noqa: E402


def _sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _norm(a):
    return _dot(a, a) ** 0.5


def _dist(a, b):
    return _norm(_sub(a, b))


def classify(vertices, faces, *, eps: float, seam: float) -> dict:
    """Split the surviving free edges into coincident / t_junction / gap."""
    edges = free_edges(faces)
    # bucketed nearest-vertex index at ``eps`` (a uniform hash is enough: eps is tiny
    # against the model extent, so each vertex has O(1) neighbours)
    from collections import defaultdict
    cell = max(eps, 1e-12) * 2.0
    buckets: dict[tuple[int, int, int], list[int]] = defaultdict(list)
    for index, point in enumerate(vertices):
        key = (int(point[0] // cell), int(point[1] // cell), int(point[2] // cell))
        buckets[key].append(index)

    def twin(point, exclude):
        key = (int(point[0] // cell), int(point[1] // cell), int(point[2] // cell))
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for index in buckets.get((key[0] + dx, key[1] + dy, key[2] + dz), ()):
                        if index != exclude and _dist(vertices[index], point) <= eps:
                            return index
        return None

    coincident = t_junction = 0
    gaps: list[float] = []
    for a, b in edges:
        ta, tb = twin(vertices[a], a), twin(vertices[b], b)
        if ta is not None and tb is not None:
            coincident += 1
            continue
        # is either endpoint interior to another free chord within the seam budget?
        hit = False
        for other_a, other_b in edges:
            if {other_a, other_b} & {a, b}:
                continue
            for endpoint in (a, b):
                distance, _t = project_on_segment(vertices[endpoint],
                                                  vertices[other_a], vertices[other_b])
                if distance <= seam:
                    hit = True
                    break
            if hit:
                break
        if hit:
            t_junction += 1
            continue
        best = float("inf")
        mid = tuple((vertices[a][i] + vertices[b][i]) / 2.0 for i in range(3))
        for other_a, other_b in edges:
            if {other_a, other_b} & {a, b}:
                continue
            distance, _t = project_on_segment(mid, vertices[other_a], vertices[other_b])
            best = min(best, distance)
        gaps.append(best if best < float("inf") else -1.0)

    gaps.sort()
    quantiles = {}
    if gaps:
        for q in (0.0, 0.1, 0.5, 0.9, 1.0):
            quantiles[f"q{int(q * 100):02d}"] = gaps[min(len(gaps) - 1, int(q * (len(gaps) - 1)))]
    return {
        "free_edges": len(edges),
        "coincident_closable_by_weld": coincident,
        "closable_by_t_junction_split": t_junction,
        "genuine_gaps": len(gaps),
        "gap_quantiles": quantiles,
        "eps": eps,
        "seam_tolerance": seam,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("path")
    parser.add_argument("--seam-tolerance", type=float, default=None,
                        help="seam budget handed to import_iges (None = measure only)")
    parser.add_argument("--eps", type=float, default=1e-6)
    parser.add_argument("--json", type=Path, default=None)
    args = parser.parse_args(argv)

    result = import_iges(args.path, seam_tolerance=args.seam_tolerance)
    report = result.report
    verdict = classify(result.mesh.vertices, result.mesh.faces,
                       eps=args.eps, seam=args.seam_tolerance or 0.0)
    payload = {
        "source": str(args.path),
        "gates": {name: {"passed": gate.get("passed"), "value": gate.get("value")}
                  for name, gate in report.get("gates", {}).items()},
        "after_repair": report.get("after_repair"),
        "classification": verdict,
    }
    print(json.dumps(payload, indent=1))
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(payload, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())