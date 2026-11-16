"""Re-audit the 500-ship fleet under the hardened gate set, in resumable slices.

The round-3/4 acceptance statistic ("N ships, M gates, P % pass") was produced by a
gate set whose G4 only looked at the sign of the divergence-theorem volume. That sign
is blind to a locally inverted patch (see ``tests/test_shipd_gate_orientation.py``), so
the statistic was reported by an instrument that could not have failed it for that
defect. ``G4`` now also requires every shared edge to be traversed in opposite
directions by its two faces.

This script answers the question the fix raises -- does the previous number survive? --
by rebuilding the whole fleet with the hardened gate set and recording, per ship, the
winding-conflict count, all ten gate results and the nominal-field correlation inputs.

Runs in slices so it fits a short local window; each slice writes its own JSON and the
slices are merged by :func:`merge`.

    python3 tools/shipd_rebuild/reaudit_500.py --start 0 --stop 170
    python3 tools/shipd_rebuild/reaudit_500.py --merge
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tensorcad import mesh as meshlib  # noqa: E402
from tensorcad.cb_sections import (  # noqa: E402
    CalibratedHullParameters,
    build_calibrated_hull,
)
from tensorcad.shipd_population import (  # noqa: E402
    build_shipd_database,
    nominal_performance_fields,
)
from tensorcad.shipd_rebuild import rebuild_ship  # noqa: E402

OUT = ROOT / "runs" / "shipd_rebuild_r5"
SHARDS = OUT / "reaudit_shards"
N_SHIPS = 500
SEED = 42
NX, NZ = 81, 21


def _fleet_shape(ship: dict) -> CalibratedHullParameters:
    """The round-3 calibrated family for this parameter vector.

    ``rebuild_ship`` solves the same exponents internally; this duplicate call exists
    only to raise the shape's feasibility verdict for the failure ledger.
    """
    return CalibratedHullParameters(
        length=ship["L"], beam=ship["B"], draft=ship["T"],
        block_coefficient=ship["Cb"],
        midship_coefficient=ship["Cm"], waterplane_coefficient=ship["Cwp"],
        nx=NX, nz=NZ,
    )


def run_slice(start: int, stop: int) -> dict:
    ships = build_shipd_database(N_SHIPS, seed=SEED)[start:stop]
    rows = []
    t0 = time.time()
    for ship in ships:
        params = _fleet_shape(ship)
        problems = params.problems()
        nominal = nominal_performance_fields(ship)
        result = rebuild_ship(ship, nx=NX, nz=NZ, method_preference="calibrated")
        mesh_conflicts = (
            len(meshlib.edge_winding_conflicts(result.mesh))
            if result.mesh is not None else None
        )
        rows.append({
            "ship_id": ship["ship_id"],
            "ship_type": ship["ship_type"],
            "L": ship["L"], "B": ship["B"], "T": ship["T"],
            "Cb": ship["Cb"], "Cm": ship["Cm"], "Cwp": ship["Cwp"],
            "nominal_wetted_surface_m2": nominal["wetted_surface_m2"],
            "nominal_viscous_resistance_kN": nominal["viscous_resistance_kN"],
            "nominal_design_froude": nominal["design_froude"],
            "shape_feasible": not problems,
            "shape_problems": problems,
            "verdict": result.verdict,
            "strict_verdict": result.strict_verdict,
            "error_message": result.error_message,
            "volume_error": result.volume_error,
            "triangle_count": result.triangle_count,
            "winding_conflicts": mesh_conflicts,
            "signed_volume": result.volume_mesh,
            "wetted_surface_m2": result.wetted_surface_area,
            "realised_Cb": result.realised_Cb,
            "realised_Cm": result.realised_Cm,
            "realised_Cwp": result.realised_Cwp,
            "gates": [g.as_dict() for g in result.gates],
            "coefficient_gates": [g.as_dict() for g in result.coefficient_gates],
        })
    return {"start": start, "stop": stop, "n": len(rows),
            "wall_time_s": time.time() - t0, "rows": rows}


def merge() -> dict:
    rows: list[dict] = []
    for path in sorted(SHARDS.glob("slice_*.json")):
        rows.extend(json.loads(path.read_text())["rows"])
    rows.sort(key=lambda r: r["ship_id"])
    ids = [r["ship_id"] for r in rows]
    expected = list(range(1, N_SHIPS + 1))
    coverage = {"n_rows": len(rows), "complete": ids == expected,
                "missing": sorted(set(expected) - set(ids))[:20]}
    gate_names = ["G1_watertight", "G2_edge_manifold", "G3_euler_char",
                  "G4_outward_oriented", "G5_no_zero_area", "G6_volume_consistency",
                  "G7_dimension_agreement", "G8_triangle_count"]
    gate_pass = {g: 0 for g in gate_names}
    for row in rows:
        for gate in row["gates"]:
            if gate["name"] in gate_pass and gate["passed"]:
                gate_pass[gate["name"]] += 1
    n = len(rows)
    verified = [r for r in rows if r["verdict"] == "VERIFIED"]
    strict = [r for r in rows if r["strict_verdict"] == "VERIFIED"]
    vol = [r["volume_error"] for r in rows if r["volume_error"] is not None]
    dev = []
    for row in rows:
        if row["wetted_surface_m2"] and row["nominal_wetted_surface_m2"]:
            dev.append((row["wetted_surface_m2"] - row["nominal_wetted_surface_m2"])
                       / row["nominal_wetted_surface_m2"])
    conflicts = [r["winding_conflicts"] for r in rows if r["winding_conflicts"] is not None]
    summary = {
        "n_ships": n,
        "coverage": coverage,
        "pass_rate": len(verified) / n if n else 0.0,
        "strict_pass_rate": len(strict) / n if n else 0.0,
        "gate_pass_rates": {g: gate_pass[g] / n if n else 0.0 for g in gate_names},
        "volume_error_mean": sum(vol) / len(vol) if vol else None,
        "volume_error_max": max(vol) if vol else None,
        "winding_conflicts_max": max(conflicts) if conflicts else None,
        "winding_conflicts_total": sum(conflicts) if conflicts else None,
        "ships_with_conflicts": sum(1 for c in conflicts if c),
        "surface_deviation": {
            "n": len(dev),
            "mean": sum(dev) / len(dev) if dev else None,
            "mean_abs": sum(abs(d) for d in dev) / len(dev) if dev else None,
            "max_abs": max((abs(d) for d in dev), default=None),
            "n_over_3pct": sum(1 for d in dev if abs(d) > 0.03),
            "n_over_15pct": sum(1 for d in dev if abs(d) > 0.15),
        },
        "shape_infeasible": sum(1 for r in rows if not r["shape_feasible"]),
        "failure_ledger": _ledger(rows),
    }
    payload = {"summary": summary, "rows": rows}
    (OUT / "reaudit_500.json").write_text(json.dumps(payload, indent=1))
    return summary


def _ledger(rows: list[dict]) -> dict[str, list[int]]:
    """Classify every ship into exactly one bucket, worst first.

    A verdict-based ledger says "0 failures" and tells the reader nothing about the
    population that passes while still missing its target, so the buckets include the
    silent degradations: a mesh that passes all ten gates yet reproduces the nominal
    wetted surface poorly, and a parameter vector whose calibrated shape is infeasible.
    """
    buckets: dict[str, list[int]] = {}

    def put(key: str, ship_id: int) -> None:
        buckets.setdefault(key, []).append(ship_id)

    for row in rows:
        if row["verdict"] == "REFUSED":
            put("refused_by_pipeline", row["ship_id"])
            continue
        if row["verdict"] != "VERIFIED" or row["strict_verdict"] != "VERIFIED":
            failed = [g["name"] for g in row["gates"] if not g["passed"]]
            failed += [g["name"] for g in row["coefficient_gates"] if not g["passed"]]
            put("gate_failed:" + ",".join(failed) if failed else "verdict_failed", row["ship_id"])
            continue
        nominal = row["nominal_wetted_surface_m2"]
        measured = row["wetted_surface_m2"]
        if nominal and measured:
            dev = abs(measured - nominal) / nominal
            if dev > 0.03:
                put("passes_gates_surface_off_target_over_3pct", row["ship_id"])
                continue
        put("clean", row["ship_id"])
    return {key: sorted(v) for key, v in sorted(buckets.items())}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--stop", type=int, default=N_SHIPS)
    parser.add_argument("--merge", action="store_true")
    args = parser.parse_args(argv)
    OUT.mkdir(parents=True, exist_ok=True)
    SHARDS.mkdir(parents=True, exist_ok=True)
    if args.merge:
        print(json.dumps(merge(), indent=1))
        return 0
    payload = run_slice(args.start, args.stop)
    path = SHARDS / f"slice_{args.start:03d}_{args.stop:03d}.json"
    path.write_text(json.dumps(payload))
    print(f"slice {args.start}-{args.stop}: {payload['n']} ships "
          f"in {payload['wall_time_s']:.1f} s -> {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())