"""Drive the mesh-derived voxelisation ladder (10^6--10^7 triangles).

``voxelize`` in the main ladder samples an analytic field; this one voxelises an
explicit triangle mesh through ``MeshSDFField``, which is the regime an imported
CAD mesh lives in. See ``tensorcad/perf_voxel_mesh.py`` for the measurement
contract and ``docs/perf_scale_mesh_voxelize.md`` for the results.

Usage::

    # mechanism smoke (seconds)
    python scripts/perf_voxelize_mesh.py --bases 41:13 --levels 0,1 \
        --grids 8:4:2,16:8:4 --out /tmp/voxel_mesh.json

    # the real ladder (minutes; submit through the compute scheduler)
    python scripts/perf_voxelize_mesh.py --bases 41:13 --levels 0,2,3,4,5,6 \
        --grids 16:8:4,24:12:6,32:16:8 \
        --output benchmarks/perf_scale/artifacts/ladder_voxelize_mesh.json
"""
from __future__ import annotations

import argparse
import json
import shlex
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tensorcad.perf_voxel_mesh import (  # noqa: E402
    DEFAULT_BASES,
    DEFAULT_GRIDS,
    DEFAULT_LEVELS,
    check_section,
    render_section,
    run_voxelize_mesh_benchmark,
)


def _pairs(text: str) -> list[tuple[int, int]]:
    out = []
    for item in text.split(","):
        item = item.strip()
        if not item:
            continue
        nx, nz = item.split(":")
        out.append((int(nx), int(nz)))
    return out


def _ints(text: str) -> list[int]:
    return [int(piece) for piece in text.split(",") if piece.strip()]


def _grids(text: str) -> list[tuple[int, int, int]]:
    out = []
    for item in text.split(","):
        item = item.strip()
        if not item:
            continue
        parts = item.split(":")
        if len(parts) != 3:
            raise SystemExit(f"grid {item!r} is not nx:ny:nz")
        out.append(tuple(int(p) for p in parts))  # type: ignore[arg-type]
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--bases", default=",".join(f"{a}:{b}" for a, b in DEFAULT_BASES))
    parser.add_argument("--levels", default=",".join(str(v) for v in DEFAULT_LEVELS))
    parser.add_argument("--grids", default=",".join(":".join(str(v) for v in g) for g in DEFAULT_GRIDS))
    parser.add_argument("--band-factor", type=float, default=1.0)
    parser.add_argument("--max-triangles", type=int, default=30_000_000)
    parser.add_argument("--budget-seconds", type=float, default=600.0)
    parser.add_argument("--output", "--out", dest="output", default=None)
    parser.add_argument("--doc", default=None, help="write/refresh the Markdown section here")
    parser.add_argument("--check", action="store_true", help="verify --doc matches a re-run")
    args = parser.parse_args(argv)

    artifact = run_voxelize_mesh_benchmark(
        bases=_pairs(args.bases),
        levels_ladder=_ints(args.levels),
        grids=_grids(args.grids),
        band_factor=args.band_factor,
        max_triangles=args.max_triangles,
        budget_seconds=args.budget_seconds,
        command="python " + " ".join(shlex.quote(a) for a in (sys.argv[0:1] + (argv if argv is not None else sys.argv[1:]))),
    )

    if args.output:
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        Path(args.output).write_text(json.dumps(artifact, indent=2, ensure_ascii=False) + "\n")
        print(f"wrote {args.output}")

    section = render_section(artifact)
    if args.doc:
        doc = Path(args.doc)
        if args.check:
            problems = check_section(doc.read_text(), artifact)
            for problem in problems:
                print(f"DRIFT: {problem}", file=sys.stderr)
            print("section matches artifact" if not problems else "section drifted")
            return 0 if not problems else 1
        text = doc.read_text()
        begin, end = "<!-- perf-meshvoxel:section begin -->", "<!-- perf-meshvoxel:section end -->"
        if begin in text and end in text:
            head, rest = text.split(begin, 1)
            tail = rest.split(end, 1)[1]
            text = f"{head}{begin}\n{section}\n{end}{tail}"
        else:
            text = f"{text.rstrip()}\n\n{begin}\n{section}\n{end}\n"
        doc.write_text(text)
        print(f"updated section in {args.doc}")

    print(f"verdict={artifact['verdict']} gates_checked={artifact['gates_checked']} "
          f"not_run={artifact['gates_not_run']} wall={artifact['seconds']:.1f}s")
    for gate in artifact["gates"]:
        print(f"  {gate['name']:<28} {gate['status']:<9} {gate['detail']}")
    for warning in artifact["cluster_warnings"]:
        print(f"  WARN {warning}")
    return 0 if artifact["verdict"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
