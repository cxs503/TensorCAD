#!/usr/bin/env python3
"""Is a surviving mutant *observable*? Two very different conclusions hide behind "survived".

The 2026-11-06 audit set the precedent (`HUL-4-midship-split`): a mutant the suite does not
notice is only a verification gap if it changes an answer a user can see. If the mutation
cannot move any output -- a tie-break on a value that is equal by construction, a stopping
tolerance the iteration overshoots anyway -- then no test can kill it and the honest record
is "structural no-op", not "gap".

This probe decides that question by measurement instead of by argument. For each surviving
mutant it patches the module exactly as the audit did, imports a fresh copy, evaluates a
fixed battery of quantities on representative cases, restores the source, and reports the
relative change against the unmutated run.

  * `observable`  -- some quantity moves; the survivor is a real hole in the suite.
  * `no-op`       -- nothing moves to round-off; no test could have caught it.
  * `verdict`     -- a *status* (pass / fail / indeterminate / not_assessed) flips, which is
                     the worst case: a rule check that silently answers a different question.

Usage::

    python3 tools/regulatory_gate_power/observability.py --out runs/2026-11-08-regulatory-observability.json
"""

from __future__ import annotations

import argparse
import importlib
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from tools.regulatory_gate_power.mutate import (  # noqa: E402
    FILES,
    MUTANTS,
    ROOT,
    clean_tree,
    sha256,
)

SURVIVORS = ("STA-2-flooding-angle-clamp", "DAM-2-permeability-waterplane", "DAM-3-moment-tolerance")

MODULES = {"stat": "tensorcad.statutory", "is": "tensorcad.intact_stability", "dmg": "tensorcad.damage"}


def _box_mesh(length: float, beam: float, draft: float, freeboard: float = 1.0):
    from tensorcad.geometry import Mesh

    xs, ys, zs = [0.0, length], [-0.5 * beam, 0.5 * beam], [-draft, freeboard]
    vertices: list[tuple[float, float, float]] = []
    index: dict[tuple[float, float, float], int] = {}
    faces: list[tuple[int, int, int]] = []

    def vid(point):
        if point not in index:
            index[point] = len(vertices)
            vertices.append(point)
        return index[point]

    def quad(a, b, c, d):
        faces.append((vid(a), vid(b), vid(c)))
        faces.append((vid(a), vid(c), vid(d)))

    quad((xs[0], ys[0], zs[0]), (xs[0], ys[1], zs[0]), (xs[1], ys[1], zs[0]), (xs[1], ys[0], zs[0]))
    quad((xs[0], ys[0], zs[1]), (xs[1], ys[0], zs[1]), (xs[1], ys[1], zs[1]), (xs[0], ys[1], zs[1]))
    quad((xs[0], ys[0], zs[0]), (xs[1], ys[0], zs[0]), (xs[1], ys[0], zs[1]), (xs[0], ys[0], zs[1]))
    quad((xs[1], ys[1], zs[0]), (xs[0], ys[1], zs[0]), (xs[0], ys[1], zs[1]), (xs[1], ys[1], zs[1]))
    quad((xs[0], ys[1], zs[0]), (xs[0], ys[0], zs[0]), (xs[0], ys[0], zs[1]), (xs[0], ys[1], zs[1]))
    quad((xs[1], ys[0], zs[0]), (xs[1], ys[1], zs[0]), (xs[1], ys[1], zs[1]), (xs[1], ys[0], zs[1]))
    return Mesh(vertices=vertices, faces=faces)


def statutory_metrics() -> dict:
    """A curve whose peak (45 deg) lies past the flooding angle (35 deg)."""
    from tensorcad import statutory

    def curve(angle_deg: float) -> float:
        return 0.8 * math.sin(math.radians(2.0 * angle_deg))

    assessment = statutory.check_intact_stability(curve, flooding_angle_deg=35.0)
    by_code = {check.code: check for check in assessment.checks}
    return {
        "2.2.3_status": by_code["2.2.3"].status,
        "2.2.3_attained": by_code["2.2.3"].attained,
        "angle_of_max_gz_deg": assessment.angle_of_max_gz_deg,
        "verdict": assessment.verdict,
    }


def damage_metrics() -> dict:
    from tensorcad.damage import Compartment, flooded_equilibrium

    L, B, T = 10.0, 2.0, 0.8
    displacement = L * B * T
    box = _box_mesh(L, B, T)
    out: dict = {}
    cases = {
        "midship_mu1.0": Compartment(4.0, 6.0),
        "midship_mu0.6": Compartment(4.0, 6.0, permeability=0.6),
        "aft_mu0.6": Compartment(0.0, 2.0, permeability=0.6),
        "port_half_mu0.6": Compartment(4.0, 6.0, y0=-1.0, y1=0.0, permeability=0.6),
    }
    for name, compartment in cases.items():
        state = flooded_equilibrium(box, compartment, displacement, kg=0.6)
        out[name] = {
            "draft_at_reference": state.draft_at_reference,
            "trim_angle_deg": state.trim_angle_deg,
            "heel_angle_deg": state.heel_angle_deg,
            "waterplane_area": state.waterplane_area,
            "transverse_second_moment": state.transverse_second_moment,
            "volume_net": state.volume_net,
            "tcb": state.tcb,
            "gm_transverse": state.gm_transverse,
        }
    return out


BUILDERS = {"stat": statutory_metrics, "dmg": damage_metrics}


def evaluate(file_key: str) -> dict:
    module_name = MODULES[file_key]
    for name in [n for n in sys.modules if n == module_name or n.startswith(module_name + ".")]:
        del sys.modules[name]
    if file_key == "stat":
        importlib.import_module("tensorcad")
    importlib.import_module(module_name)
    return BUILDERS[file_key]()


def flatten(prefix: str, value) -> dict:
    if callable(value) and not isinstance(value, type):
        value = value()
    if isinstance(value, dict):
        out: dict = {}
        for key, item in value.items():
            out.update(flatten(f"{prefix}.{key}" if prefix else str(key), item))
        return out
    return {prefix: value}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    out_path = ROOT / args.out
    out_path.parent.mkdir(parents=True, exist_ok=True)

    originals = {key: path.read_bytes() for key, path in FILES.items()}
    digests = {key: sha256(data) for key, data in originals.items()}

    def restore_all() -> None:
        for key, data in originals.items():
            FILES[key].write_bytes(data)
            assert sha256(FILES[key].read_bytes()) == digests[key], f"restore failed for {key}"

    if not clean_tree():
        print("[abort] tree dirty; refusing to patch")
        return 1

    baseline = {key: flatten("", evaluate(key)) for key in ("stat", "dmg")}
    report: dict = {"baseline": baseline, "mutants": {}}

    for mutant in MUTANTS:
        if mutant["id"] not in SURVIVORS:
            continue
        key = mutant["file"]
        text = originals[key].decode()
        FILES[key].write_text(text.replace(mutant["old"], mutant["new"]))
        try:
            mutated = flatten("", evaluate(key))
        finally:
            restore_all()
        deltas = {}
        flips = {}
        for name, value in baseline[key].items():
            other = mutated[name]
            if isinstance(value, str) or isinstance(other, str):
                if value != other:
                    flips[name] = {"before": value, "after": other}
                continue
            scale = max(abs(value), 1e-30)
            relative = abs(other - value) / scale
            if relative > 1e-12:
                deltas[name] = {"before": value, "after": other, "relative_change": relative}
        worst = max((item["relative_change"] for item in deltas.values()), default=0.0)
        if flips:
            kind = "verdict"
        elif worst > 1e-12:
            kind = "observable"
        else:
            kind = "no-op"
        report["mutants"][mutant["id"]] = {
            "classification": kind,
            "worst_relative_change": worst,
            "changed_quantities": deltas,
            "flipped_statuses": flips,
            "claim": mutant["claim"],
        }
        print(f"[{kind:10s}] {mutant['id']}: worst relative change {worst:.3e}, "
              f"{len(deltas)} quantity(ies) moved, {len(flips)} status flip(s)", flush=True)

    out_path.write_text(json.dumps(report, indent=2))
    print(f"wrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())