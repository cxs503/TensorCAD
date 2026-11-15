"""Gmsh differential, driven by the Gmsh *source* oracle.

`tools/mesh_quality_external_diff.py` is the general door: it compares
:mod:`tensorcad.meshquality` against any oracle that speaks the ``W <kind> <metric>
<index> <value>`` protocol. This script is the Gmsh arm of that door, run against
``gmsh_source_oracle`` -- Gmsh's own ``qualityMeasures.cpp`` compiled from a Gmsh
release source tree (see ``tools/mesh_quality_oracles/build_gmsh_source_oracle.sh``).

Why a separate entry point: the ``libgmsh``-based ``gmsh_oracle`` reports
``minSICN``/``minEdge``/``maxEdge`` as well, while the source oracle deliberately
reports only the measures it can take from ``qualityMeasures.cpp`` itself
(``gamma``, ``eta``, tet ``volume``). Keeping the two mappings apart means neither
arm can silently claim a pair the other measured.

Usage::

    python tools/mesh_quality_gmsh_source_diff.py --oracle /tmp/gso_bin
    python tools/mesh_quality_gmsh_source_diff.py --oracle ... --require-oracle

Exit codes: 0 every compared pair within 1e-6; 1 a pair disagreed; 2 the oracle did
not run (and ``--require-oracle`` was given).
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "tools"))

import mesh_quality_external_diff as drv  # noqa: E402

TOLERANCE = 1.0e-6

#: our metric -> Gmsh measure name, for the measures qualityMeasures.cpp computes.
MAPPING: dict[str, dict[str, str]] = {
    "tri": {"gamma": "gamma", "eta": "eta"},
    "quad": {"eta": "eta"},
    "tet": {"gamma": "gamma", "eta": "eta", "volume": "volume"},
}

ARTIFACT = REPO / "artifacts" / "mesh_quality_external_diff" / "gmsh_source_diff.json"


def ours_for(kind: str, points, elems):
    """Per-element measures from the library (same functions the driver uses)."""
    fn = drv.METRIC_FN[kind]
    try:
        return fn(points, elems)
    except TypeError:
        return fn(elems, points)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--oracle", required=False, default=None)
    ap.add_argument("--require-oracle", action="store_true")
    ap.add_argument("--artifact", default=str(ARTIFACT))
    ap.add_argument("--tolerance", type=float, default=TOLERANCE)
    args = ap.parse_args(argv)

    tmpdir = Path(args.artifact).resolve().parent
    tmpdir.mkdir(parents=True, exist_ok=True)
    report: dict = {
        "driver": "tools/mesh_quality_gmsh_source_diff.py",
        "oracle": "gmsh_source_oracle (Gmsh qualityMeasures.cpp, source build)",
        "tolerance": args.tolerance,
        "meshes": {},
    }

    disagreed = 0
    compared_pairs = 0
    ran_any = False

    for name, (kind, points, elems) in sorted(drv.meshes().items()):
        if kind not in MAPPING:
            report["meshes"][name] = {"kind": kind, "status": "SKIPPED",
                                      "why": "kind not computed by qualityMeasures.cpp"}
            continue
        mesh_file = tmpdir / f"_{name}.mkm"
        drv.write_mesh_file(mesh_file, kind, points, elems)
        status, external = drv.run_oracle(args.oracle, mesh_file, "W")
        mesh_file.unlink(missing_ok=True)
        entry: dict = {"kind": kind, "elements": len(elems), "status": status}
        if status != "RAN" or kind not in external:
            entry["why"] = "oracle did not run"
            report["meshes"][name] = entry
            continue
        ran_any = True
        ours = ours_for(kind, points, elems)
        cmp = drv.compare(ours, external[kind], MAPPING[kind])
        entry["metrics"] = cmp
        worst = 0.0
        for metric, info in cmp.items():
            if info.get("status") != "COMPARED":
                continue
            compared_pairs += 1
            worst = max(worst, info.get("max_abs_diff_rel", 0.0) or 0.0)
            if (info.get("beyond_tolerance", 0) or 0) > 0:
                disagreed += 1
        entry["worst_rel_diff"] = worst
        report["meshes"][name] = entry

    report["pairs_compared"] = compared_pairs
    report["pairs_beyond_tolerance"] = disagreed
    report["oracle_ran"] = ran_any
    Path(args.artifact).parent.mkdir(parents=True, exist_ok=True)
    Path(args.artifact).write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")

    for name, entry in report["meshes"].items():
        if entry.get("status") != "RAN":
            print(f"{name:16s} {entry['status']}")
            continue
        for metric, info in entry.get("metrics", {}).items():
            print(f"{name:16s} {metric:10s} {info.get('status'):9s} "
                  f"worst_rel={info.get('max_abs_diff_rel', float('nan')):.3e} "
                  f"beyond_tol={info.get('beyond_tolerance', 0)}/{entry['elements']}")
    print(f"pairs compared {compared_pairs}, beyond {args.tolerance:g}: {disagreed}")

    if not ran_any:
        return 2 if args.require_oracle else 0
    return 1 if disagreed else 0


if __name__ == "__main__":
    raise SystemExit(main())