#!/usr/bin/env python3
"""``tensorcad-dtmb5415`` -- build a DTMB 5415 hull from a 27-variable design vector.

This is the solver-side entry point of :mod:`tensorcad.dtmb5415` (the parametric rebuild of
the Zenodo 9001 x 27 design space).  It exists so the model is reachable from the command
line and from other tools, not just from a research script:

.. code-block:: console

    # reference hull (v = 0), JSON record on stdout
    tensorcad-dtmb5415 --json

    # a stored dataset configuration (needs the 1.0 GB design database)
    tensorcad-dtmb5415 --config-index 430 --stl out/5415_430.stl

    # an explicit design vector, with the analytic dV/dv chain
    tensorcad-dtmb5415 --design-vector "0.1,-0.2,0,..." --gradient --out myrun

Exit codes follow ``tensorcad-hull``: 0 ok, 1 geometry quality gate failed, 2 invalid input.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Sequence

#: Quality gate used by ``--gate`` (the same thresholds as run_acceptance.py).
GATE = {
    "finite_vertices": lambda q: q["non_finite_vertices"] == 0,
    "no_degenerate_faces": lambda q: q["degenerate_faces"] == 0,
    "no_non_manifold_edges": lambda q: q["non_manifold_edges"] == 0,
    "min_angle_deg": lambda q: q["min_angle_deg"] >= 20.0,
    "frac_angle_below_5deg": lambda q: q["frac_angle_below_5deg"] <= 1e-3,
    "aspect_p99": lambda q: q["aspect_p99"] <= 10.0,
    "stations_monotone": lambda q: bool(q["stations_monotone"]),
    "half_breadth_non_negative": lambda q: q["half_breadth_min_m"] >= -1e-6,
    "edge_max_m": lambda q: q["edge_max_m"] <= 8.0,
}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="tensorcad-dtmb5415",
        description="Parametric DTMB 5415 hull: 27 design variables -> mesh + hydrostatics.",
    )
    src = parser.add_mutually_exclusive_group()
    src.add_argument("--design-vector", default=None,
                     help="comma separated 27 normalised design variables in [-1, 1]; "
                          "default is the reference hull (all zeros)")
    src.add_argument("--config-index", type=int, default=None,
                     help="row of the Zenodo design database to reproduce (reads the 1.0 GB "
                          "DTMB5415-database.mat; provenance is echoed in the record)")
    parser.add_argument("--mat", default=None,
                        help="override the design database path used by --config-index")
    parser.add_argument("--waterline", type=float, default=None,
                        help="waterline z in metres for the hydrostatics (default: the design "
                             "draft, 0.248 m * 142 / 5.720)")
    parser.add_argument("--out", default="tensorcad-output/dtmb5415",
                        help="directory for the JSON record (and the STL when asked)")
    parser.add_argument("--stl", default=None,
                        help="also write a binary STL of the hull shell to this path")
    parser.add_argument("--gradient", action="store_true",
                        help="include the analytic d volume / d v chain (27 values)")
    parser.add_argument("--gate", action="store_true",
                        help="exit 1 unless the mesh passes the quality gate")
    parser.add_argument("--json", action="store_true", help="print the full record")
    parser.add_argument("--no-write", action="store_true",
                        help="compute and print only, write no files")
    return parser


def _design_vector(args: argparse.Namespace, out: dict[str, Any]) -> Sequence[float]:
    import numpy as np

    if args.config_index is not None:
        import scipy.io as sio

        from .dtmb5415 import DATASET

        mat = args.mat or DATASET["source_file"]
        path = Path(mat)
        if not path.exists():
            raise ValueError(f"design database not found: {path}")
        d = sio.loadmat(str(path), variable_names=["variablesMatrix"])
        v = np.asarray(d["variablesMatrix"], dtype=float).T
        if not 0 <= args.config_index < v.shape[0]:
            raise ValueError(f"--config-index must be in [0, {v.shape[0] - 1}]")
        out["design_source"] = {"database": str(path),
                                "database_sha256": DATASET["source_sha256"],
                                "config_index": int(args.config_index)}
        return v[args.config_index]
    if args.design_vector is None:
        out["design_source"] = {"reference_hull": True}
        return []
    vals = [float(x) for x in args.design_vector.replace(" ", "").split(",") if x != ""]
    out["design_source"] = {"design_vector_explicit": True}
    return vals


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    from .dtmb5415 import DTMB5415Space

    def fail(stage: str, exc: Exception) -> int:
        print(f"tensorcad-dtmb5415: {stage} failed: {exc}", file=sys.stderr)
        if args.json:
            print(json.dumps({"error": str(exc), "stage": stage}, sort_keys=True))
        return 2

    try:
        space = DTMB5415Space()
        record: dict[str, Any] = {}
        v = _design_vector(args, record)
        v = space.zero() if len(v) == 0 else v
        record.update(space.to_dict(v))
        record["volume_m3"] = float(space.volume(v))
        if args.waterline is not None:
            record["hydrostatics"] = space.hydrostatics(v, waterline_z=float(args.waterline))
        if args.gradient:
            record["volume_gradient"] = [float(x) for x in space.volume_gradient(v)]
        gate = {name: bool(fn(record["quality"])) for name, fn in GATE.items()}
        record["quality_gate"] = {"checks": gate, "pass": all(gate.values())}
    except Exception as exc:                     # noqa: BLE001 - one error contract
        return fail("build", exc)

    written: dict[str, str] = {}
    if not args.no_write:
        out = Path(args.out)
        out.mkdir(parents=True, exist_ok=True)
        meta = out / "dtmb5415_record.json"
        meta.write_text(json.dumps(record, indent=2))
        written["record"] = str(meta)
        if args.stl:
            from .export import write_stl

            stl = Path(args.stl)
            stl.parent.mkdir(parents=True, exist_ok=True)
            write_stl(space.mesh(v), stl)
            written["stl"] = str(stl)
    record["artifacts"] = written

    if args.json:
        print(json.dumps(record, indent=2, sort_keys=True))
    else:
        h = record["hydrostatics"]
        print(json.dumps({
            "design_source": record["design_source"],
            "geometry": {"vertices": record["quality"]["n_vertices"],
                         "faces": record["quality"]["n_faces"],
                         "edge_median_m": round(record["quality"]["edge_median_m"], 4)},
            "hydrostatics": {k: (round(vv, 4) if isinstance(vv, float) else vv)
                             for k, vv in h.items()},
            "volume_m3": round(record["volume_m3"], 3),
            "quality_gate": record["quality_gate"]["pass"],
            "artifacts": written,
        }, indent=2))
    return 0 if record["quality_gate"]["pass"] else 1


if __name__ == "__main__":  # pragma: no cover - module entry point
    raise SystemExit(main(sys.argv[1:]))