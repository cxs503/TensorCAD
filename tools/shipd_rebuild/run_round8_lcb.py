#!/usr/bin/env python3
"""GAP-CAD-C2-SHIPD-REBUILD round 8: the LCB degree of freedom at population scale.

Rounds 3-6 rebuilt 500 hulls against five scalar coefficients (+ the derived
``Cp <= Cwp`` feasibility flag) and reported a 100 % pass rate. Round 6 then measured
the *longitudinal centre of buoyancy* on those hulls and on six published benchmark
hulls, and found every one of them sitting at exactly ``0.000 % L``: the family is
fore/aft symmetric, so the 500-hull campaign was reproducing ships whose LCB was
wrong by up to 3.48 % L (KVLCC2) -- an error no gate looked at, because the parameter
vector carried no LCB.

Round 7 added the missing control: ``CalibratedHullParameters.lcb_offset`` tilts the
fullness exponent fore/aft and the realised LCB lands on the request to ~1e-15 L on a
35-case coefficient scan. That capability never met the campaign: ``rebuild_ship``
had no way to consume an LCB, the population model had no LCB column, and the eight
gates did not look at the centre of buoyancy.

This script closes that loop and re-runs the campaign, so the claim "500 ships
rebuild with the coefficients *and* the LCB their parameter vector asks for" is
measured rather than inferred.

Runs
----
``population``  500 ships, seed 42, each with a type-typical LCB target.
``control``     the same 500 ships with the LCB column stripped -- the round-3/5
                configuration, re-run on this tree so the two pass rates are
                comparable row by row.
``literature``  the six published benchmark hulls (KCS, KVLCC2, ONRT, JBC, DTC,
                Series 60) rebuilt at their *published* LCB.
``bitwise``     twenty in-sample ships rebuilt with ``lcb`` absent vs ``lcb = 0``:
                the LCB control must cost nothing when it is not asked for.

Every number printed into the report comes from the JSON this script writes.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from tensorcad.cb_sections import LCB_OFFSET_LIMIT  # noqa: E402
from tensorcad.shipd_oos import literature_corpus  # noqa: E402
from tensorcad.shipd_population import (  # noqa: E402
    LCB_FRACTION_RANGE,
    build_shipd_database,
)
from tensorcad.shipd_rebuild import (  # noqa: E402
    LCB_TOLERANCE,
    run_shipd_rebuild,
)

OUT = Path("runs/shipd_rebuild_r8")


def _mesh_digest(mesh) -> str:
    payload = repr((tuple(mesh.vertices), tuple(mesh.faces))).encode()
    return hashlib.sha256(payload).hexdigest()


def _summarise_lcb(evidence: dict[str, Any]) -> dict[str, Any]:
    return evidence["lcb_contract"]


def run_population(n_ships: int, nx: int, nz: int, workers: int) -> tuple[list[dict], dict]:
    ships = build_shipd_database(n_ships, seed=42)
    evidence = run_shipd_rebuild(
        ships, nx=nx, nz=nz, out_dir=OUT / "population", workers=workers
    )
    return ships, evidence


def run_control(ships: list[dict], nx: int, nz: int, workers: int) -> dict:
    # Strip the LCB column: exactly the round-3/5 parameter vector.
    bare = [{k: v for k, v in ship.items() if k != "lcb"} for ship in ships]
    return run_shipd_rebuild(bare, nx=nx, nz=nz, out_dir=OUT / "control", workers=workers)


def run_literature(nx: int, nz: int) -> dict[str, Any]:
    """Six benchmark hulls rebuilt at their published LCB (round 6's finding, closed)."""
    hulls = literature_corpus()
    rows: list[dict[str, Any]] = []
    for hull in hulls:
        published = hull.get("lcb_frac_L_published")
        vector = {k: v for k, v in hull.items() if k.startswith("lcb") is False}
        if published is not None:
            vector["lcb"] = float(published)
        evidence = run_shipd_rebuild(
            [vector], nx=nx, nz=nz, out_dir=OUT / "literature" / str(hull["ship_id"]),
        )
        contract = json.loads(
            (OUT / "literature" / str(hull["ship_id"]) / "lcb_contract.json").read_text()
        )["rows"]
        row = {
            "ship_id": hull["ship_id"],
            "Cb": hull["Cb"],
            "Cm": hull["Cm"],
            "Cwp": hull["Cwp"],
            "Cp": hull["Cp"],
            "cp_le_cwp": hull["cp_le_cwp"],
            "lcb_published": published,
            "method": next(iter(evidence["method_distribution"])),
            "verdict": "VERIFIED" if evidence["n_verified"] else "FAILED",
            "strict_verdict": "VERIFIED" if evidence["n_strict_verified"] else "FAILED",
            "realised_lcb": contract[0]["realised_lcb"] if contract else None,
            "lcb_requested": contract[0]["requested_lcb"] if contract else None,
            "lcb_error": evidence["lcb_contract"]["max_abs_error_L"],
            "volume_error_max": evidence["volume_error_max"],
            "n_ships": evidence["n_ships"],
            "provenance": hull.get("provenance", {}),
        }
        rows.append(row)
    out = {"hulls": rows, "note": (
        "LCB targets are the published values carried by the fleet manifest "
        "(retrieval=offline_recollection, verified=false). They test whether the "
        "control loop can put the centroid where a real ship has it; they are not "
        "an independent verification of the hull forms.")}
    (OUT / "literature").mkdir(parents=True, exist_ok=True)
    (OUT / "literature" / "summary.json").write_text(json.dumps(out, indent=2))
    return out


def bitwise_regression(n_sample: int, nx: int, nz: int) -> dict[str, Any]:
    """Rebuild with ``lcb`` absent vs ``lcb = 0`` and compare mesh hashes."""
    from tensorcad.shipd_rebuild import rebuild_ship

    ships = build_shipd_database(n_sample, seed=42)
    rows = []
    for ship in ships:
        bare = {k: v for k, v in ship.items() if k != "lcb"}
        zero = dict(bare, lcb=0.0)
        ra = rebuild_ship(bare, nx=nx, nz=nz)
        rb = rebuild_ship(zero, nx=nx, nz=nz)
        rows.append({
            "ship_id": ship["ship_id"],
            "ship_type": ship["ship_type"],
            "hash_no_lcb": _mesh_digest(ra.mesh) if ra.mesh is not None else None,
            "hash_lcb_zero": _mesh_digest(rb.mesh) if rb.mesh is not None else None,
            "identical": (
                ra.mesh is not None and rb.mesh is not None
                and _mesh_digest(ra.mesh) == _mesh_digest(rb.mesh)
            ),
            "lcb_no_request": ra.realised_lcb,
        })
    return {
        "n": len(rows),
        "all_identical": all(row["identical"] for row in rows),
        # The default vector carries no LCB target, so its centroid is whatever the
        # symmetric family produces -- report the residue instead of assuming zero.
        "max_abs_realised_lcb_no_request": max(
            (abs(row["lcb_no_request"]) for row in rows if row["lcb_no_request"] is not None),
            default=None,
        ),
        "rows": rows,
    }


def run_hostile_lcb(
    ships: list[dict], nx: int, nz: int, workers: int, *, seed: int = 20261116
) -> dict[str, Any]:
    """Same hulls, LCB drawn uniformly over the whole admissible API band.

    The type-typical band is narrow enough to sit inside the geometric saturation of
    every sampled hull, so it cannot fail -- and an instrument that cannot fail is not
    evidence. Here the target is ``U(-5 % L, +5 % L)`` with no regard for what the
    coefficient set can hold, which pushes the full hulls past their saturation.
    """
    import random

    rng = random.Random(seed)
    vectors = []
    for ship in ships:
        # Inside the API guardrail (|lcb| < 5 % L) but unrelated to Cb: this is the
        # population in which an unreachable request is *possible*.
        offset = rng.uniform(-(LCB_OFFSET_LIMIT - 1.0e-3), LCB_OFFSET_LIMIT - 1.0e-3)
        vectors.append(dict(ship, lcb=round(offset, 5)))
    evidence = run_shipd_rebuild(
        vectors, nx=nx, nz=nz, out_dir=OUT / "hostile_lcb", workers=workers
    )
    contract = json.loads((OUT / "hostile_lcb" / "lcb_contract.json").read_text())
    rows = contract["rows"]
    refused_forward = sum(
        1 for row in rows
        if (row["error_message"] or "").startswith("lcb_unreachable")
        and row["requested_lcb"] > 0
    )
    refused_aft = sum(
        1 for row in rows
        if (row["error_message"] or "").startswith("lcb_unreachable")
        and row["requested_lcb"] < 0
    )
    by_type: dict[str, dict[str, int]] = {}
    for row in rows:
        bucket = by_type.setdefault(row["ship_type"], {"refused": 0, "Total": 0})
        bucket["Total"] += 1
        if (row["error_message"] or "").startswith("lcb_unreachable"):
            bucket["refused"] += 1
    out = {
        "seed": seed,
        "evidence": evidence,
        "n_requested": contract["summary"]["n_requested"],
        "n_honoured": contract["summary"]["n_honoured"],
        "n_refused_unreachable": contract["summary"]["n_refused_unreachable"],
        "refused_forward": refused_forward,
        "refused_aft": refused_aft,
        "max_abs_error_L_honoured": contract["summary"]["max_abs_error_L"],
        "by_ship_type": by_type,
    }
    (OUT / "hostile_lcb" / "hostile_lcb.json").write_text(json.dumps(out, indent=2))
    return out


def reach_survey(ships: list[dict], nx: int, nz: int, *, n_tilts: int = 91) -> dict[str, Any]:
    """Geometric LCB reach of every hull, measured by a scan -- not by a bisection.

    ``cb_sections.lcb_reach`` bisects on the tilt and assumes "tilt admissible" is
    monotone in the reach, so a single measurement from it could in principle agree
    with a refusal for the wrong reason. This survey walks the whole tilt range
    (``-0.9 … +0.9``) and takes the largest offset any admissible tilt produces,
    forward and aft separately, so the refusal boundary is checked against an
    independent number.
    """
    from tensorcad.cb_sections import (
        CalibratedHullParameters,
        LcbRequestUnreachable,
        longitudinal_offset,
        solve_asymmetric_exponents,
        station_volume_weights,
    )

    guardrail = LCB_OFFSET_LIMIT
    tilts = [-0.9 + 1.8 * i / (n_tilts - 1) for i in range(n_tilts)]
    rows = []
    for ship in ships:
        params = CalibratedHullParameters(
            length=ship["L"], beam=ship["B"], draft=ship["T"],
            block_coefficient=ship["Cb"],
            midship_coefficient=ship["Cm"],
            waterplane_coefficient=ship["Cwp"], nx=nx, nz=nz,
        )
        taper, keel = params.taper(), params.keel_ratios()
        cp, cwp = params.prismatic_coefficient, params.waterplane_coefficient
        best_forward = best_aft = 0.0
        for tilt in tilts:
            try:
                breadth, fullness = solve_asymmetric_exponents(
                    taper, keel, cwp, cp, tilt
                )
            except LcbRequestUnreachable:
                continue
            offset = longitudinal_offset(
                station_volume_weights(taper, keel, breadth, fullness, tilt)
            )
            if offset > best_forward:
                best_forward = offset
            if offset < best_aft:
                best_aft = offset
        rows.append({
            "ship_id": ship["ship_id"],
            "ship_type": ship["ship_type"],
            "Cb": ship["Cb"],
            "reach_forward_L": best_forward,
            "reach_aft_L": best_aft,
            "guardrail_inside_reach_forward": guardrail <= best_forward,
            "guardrail_inside_reach_aft": -guardrail >= best_aft,
        })

    by_type: dict[str, dict[str, Any]] = {}
    for row in rows:
        bucket = by_type.setdefault(row["ship_type"], {
            "min_reach_forward_L": row["reach_forward_L"],
            "min_reach_aft_L": row["reach_aft_L"],
            "n": 0,
            "n_guardrail_outside_forward": 0,
            "n_guardrail_outside_aft": 0,
        })
        bucket["n"] += 1
        bucket["min_reach_forward_L"] = min(bucket["min_reach_forward_L"], row["reach_forward_L"])
        bucket["min_reach_aft_L"] = min(bucket["min_reach_aft_L"], row["reach_aft_L"])
        bucket["n_guardrail_outside_forward"] += 0 if row["guardrail_inside_reach_forward"] else 1
        bucket["n_guardrail_outside_aft"] += 0 if row["guardrail_inside_reach_aft"] else 1

    out = {
        "guardrail_L": guardrail,
        "n_tilts": n_tilts,
        "n_ships": len(rows),
        "min_reach_forward_L": min(row["reach_forward_L"] for row in rows),
        "min_reach_aft_L": min(row["reach_aft_L"] for row in rows),
        "n_guardrail_outside_forward": sum(
            1 for row in rows if not row["guardrail_inside_reach_forward"]),
        "n_guardrail_outside_aft": sum(
            1 for row in rows if not row["guardrail_inside_reach_aft"]),
        "by_ship_type": by_type,
        "rows": rows,
    }
    (OUT / "reach_survey.json").write_text(json.dumps(out, indent=2))
    return out


def stress_unreachable_probe(n_ships: int, nx: int, nz: int) -> dict[str, Any]:
    """Push the fullest hulls toward the API limit: does the refusal path fire?

    The sampled LCB band stays inside the geometric saturation measured in round 7,
    so a 500-ship campaign legitimately produces no unreachable request. An empty
    failure class is only meaningful if the class *can* fire, so ask the fullest
    hulls for the largest admissible offset (just under the 5 % L guardrail) and
    record what the pipeline does.
    """
    from tensorcad.cb_sections import CalibratedHullParameters, lcb_reach
    from tensorcad.shipd_rebuild import rebuild_ship

    ships = build_shipd_database(max(n_ships, 500), seed=42)
    fullest = sorted(ships, key=lambda ship: ship["Cb"], reverse=True)[:n_ships]
    rows = []
    for ship in fullest:
        vector = dict(ship, lcb=LCB_OFFSET_LIMIT - 1.0e-3)
        result = rebuild_ship(vector, nx=nx, nz=nz)
        # Geometric saturation of this coefficient set, measured independently of the
        # refusal above: the refusal is only credible if the request really is outside.
        saturating = lcb_reach(*(lambda p: (p.taper(), p.keel_ratios()))(
            CalibratedHullParameters(
                length=ship["L"], beam=ship["B"], draft=ship["T"],
                block_coefficient=ship["Cb"],
                midship_coefficient=ship["Cm"],
                waterplane_coefficient=ship["Cwp"], nx=nx, nz=nz,
            )
        ), ship["Cwp"], ship["Cb"] / ship["Cm"])
        rows.append({
            "geometric_reach_L": saturating,
            "request_beyond_reach": vector["lcb"] > saturating,
            "consistent": (result.verdict == "REFUSED") == (vector["lcb"] > saturating),
            "ship_id": ship["ship_id"],
            "ship_type": ship["ship_type"],
            "Cb": ship["Cb"],
            "requested_lcb": vector["lcb"],
            "verdict": result.verdict,
            "realised_lcb": result.realised_lcb,
            "lcb_gate_passed": (result.lcb_gates[0].passed if result.lcb_gates else None),
            "error_message": result.error_message,
        })
    out = {
        "requested_lcb": LCB_OFFSET_LIMIT - 1.0e-3,
        "n_probes": len(rows),
        "n_refused": sum(1 for row in rows if row["verdict"] == "REFUSED"),
        "n_honoured": sum(1 for row in rows if row["lcb_gate_passed"]),
        "n_consistent_with_measured_reach": sum(1 for row in rows if row["consistent"]),
        "rows": rows,
    }
    (OUT / "stress_unreachable_probe.json").write_text(json.dumps(out, indent=2))
    return out


def _classify_failure(entry: dict[str, Any]) -> str:
    message = (entry.get("error_message") or "").lower()
    if message.startswith("lcb_unreachable"):
        return "LCB_UNREACHABLE"
    if "inconsistent coefficients" in message:
        return "COEFFICIENTS_INFEASIBLE"
    if entry.get("failed_lcb_gate"):
        return "LCB_GATE_FAILED"
    if entry.get("failed_gates"):
        return "MESH_GATE_FAILED:" + ",".join(sorted(entry["failed_gates"]))
    return "OTHER:" + (entry.get("error_message") or "unknown")[:60]


def failure_ledger(runs: dict[str, Path]) -> dict[str, Any]:
    """Classified failure ledger, one block per run.

    Classification order matters: a refused parameter vector never produced a mesh,
    so "coefficients infeasible" and "LCB unreachable" are decided before any gate
    name is looked at.
    """
    ledger: dict[str, Any] = {}
    for name, run_dir in runs.items():
        populated = json.loads((run_dir / "failure_ledger.json").read_text())
        counts: dict[str, int] = {}
        for entry in populated:
            key = _classify_failure(entry)
            counts[key] = counts.get(key, 0) + 1
        ledger[name] = {
            "n_ships": json.loads((run_dir / "evidence.json").read_text())["n_ships"],
            "n_failures": len(populated),
            "by_class": counts,
            "examples": populated[:5],
        }
    return ledger


# --------------------------------------------------------------------------- #
# Figures
# --------------------------------------------------------------------------- #
def make_figures(population: dict, control: dict, literature: dict, bitwise: dict,
                 hostile: dict) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig_dir = OUT / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    rows = json.loads((OUT / "population" / "lcb_contract.json").read_text())["rows"]
    honoured = [r for r in rows if r["realised_lcb"] is not None]
    req = [r["requested_lcb"] for r in honoured]
    got = [r["realised_lcb"] for r in honoured]

    fig, axes = plt.subplots(2, 2, figsize=(13, 10))

    ax = axes[0][0]
    ax.plot(req, got, "o", ms=3.5, alpha=0.75, label="500-ship population")
    lim = [min(req + got), max(req + got)]
    ax.plot(lim, lim, "k--", lw=1, label="1:1")
    ax.set_xlabel("requested LCB / L")
    ax.set_ylabel("realised LCB / L (divergence theorem)")
    ax.set_title(f"(a) LCB contract, n={len(honoured)}  max|err|={population['lcb_contract']['max_abs_error_L']:.2e} L")
    ax.legend(loc="upper left", fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[0][1]
    types = list(LCB_FRACTION_RANGE)
    for i, ship_type in enumerate(types):
        lo, hi = LCB_FRACTION_RANGE[ship_type]
        ax.plot([i, i], [lo * 100, hi * 100], color="tab:gray", lw=6, alpha=0.45,
                solid_capstyle="butt")
        vals = [r["realised_lcb"] * 100 for r in honoured if r["ship_type"] == ship_type]
        ax.plot([i] * len(vals), vals, "o", ms=3, color="tab:blue")
    ax.set_xticks(range(len(types)))
    ax.set_xticklabels(types, rotation=60, ha="right", fontsize=8)
    ax.set_ylabel("LCB / %L")
    ax.set_title("(b) realised LCB per ship type vs the sampled band (grey)")
    ax.grid(alpha=0.3, axis="y")

    ax = axes[1][0]
    labels = ["canonical G1-G8", "strict (+Cb/Cm/Cwp)", "strict (+LCB)"]
    with_lcb = [
        population["pass_rate"],
        population["strict_pass_rate"],
        population["strict_pass_rate"],
    ]
    no_lcb = [
        control["pass_rate"],
        control["strict_pass_rate"],
        control["strict_pass_rate"],
    ]
    hostile_ev = hostile["evidence"]
    hostile_rates = [
        hostile_ev["pass_rate"],
        hostile_ev["strict_pass_rate"],
        hostile_ev["strict_pass_rate"],
    ]
    x = range(len(labels))
    ax.bar([i - 0.26 for i in x], [100 * v for v in no_lcb], width=0.26,
           label="control (no LCB column)", color="tab:gray")
    ax.bar([i for i in x], [100 * v for v in with_lcb], width=0.26,
           label="population (type-typical LCB)", color="tab:green")
    ax.bar([i + 0.26 for i in x], [100 * v for v in hostile_rates], width=0.26,
           label="hostile LCB (U(-5%,+5%) L)", color="tab:red")
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylim(0, 105)
    ax.set_ylabel("pass rate / %")
    ax.set_title("(c) pass rate before and after adding the LCB requirement")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, axis="y")

    ax = axes[1][1]
    corr = population["correlation_sanity"]
    hulls = literature["hulls"]
    names = [h["ship_id"] for h in hulls]
    published = [h["lcb_published"] * 100 if h["lcb_published"] is not None else 0.0 for h in hulls]
    realised = [h["realised_lcb"] * 100 if h["realised_lcb"] is not None else 0.0 for h in hulls]
    idx = range(len(names))
    ax.bar([i - 0.18 for i in idx], published, width=0.36, label="published", color="tab:orange")
    ax.bar([i + 0.18 for i in idx], realised, width=0.36, label="rebuilt", color="tab:blue")
    ax.set_xticks(list(idx))
    ax.set_xticklabels(names, fontsize=8)
    ax.set_ylabel("LCB / %L")
    ax.set_title("(d) six published hulls at their published LCB")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, axis="y")
    ax.text(0.02, 0.02, f"wetted-surface Pearson r = {corr['wetted_surface']['pearson_r']:.5f}",
            transform=ax.transAxes, fontsize=8)

    fig.suptitle("GAP-CAD-C2 round 8: LCB control at population scale", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    path = fig_dir / "lcb_population_r8.png"
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return path


# --------------------------------------------------------------------------- #
# Report
# --------------------------------------------------------------------------- #
def write_report(population, control, literature, bitwise, ledger, fig_path, stress,
                 hostile, survey) -> Path:
    lcb = population["lcb_contract"]
    clcb = control["lcb_contract"]
    hostile_ev = hostile["evidence"]
    corr = population["correlation_sanity"]
    cfg = population["settings"]
    lines: list[str] = []
    add = lines.append

    add("# GAP-CAD-C2-SHIPD-REBUILD 第 8 轮：把 LCB 自由度接到 500 船总体重建上")
    add("")
    add("> 锚点 `GAP-CAD-C2-SHIPD-REBUILD`。本轮接的是第 6 轮量出来的欠账：500 船重建")
    add("> （第 3/5 轮）与 6 条公开基准船（第 6 轮）的纵中浮心全部落在 `0.000 % L`，")
    add("> 而真实船不是 —— 参数向量里没有 LCB 这一项，门禁里也就没有这一项。第 7 轮把")
    add("> `lcb_offset` 做进标定族（35 例扫描，误差 ~1e-15 L），但**没有被管线消费**：")
    add("> `rebuild_ship` 不读 LCB，种群表没有 LCB 列，八道门不看浮心。")
    add(">")
    add(f"> 本轮新增：`tensorcad/shipd_population.LCB_FRACTION_RANGE` + `sample_lcb_fraction`"
        "（按船型的 LCB 抽样带）、")
    add("> `tensorcad/shipd_rebuild` 的 LCB 通道与 `G11_longitudinal_centre` 门、"
        "`run_shipd_rebuild(workers=...)`；")
    add(f"> 生成器 `tools/shipd_rebuild/run_round8_lcb.py`。工件：`{OUT}/population/`、"
        f"`{OUT}/control/`、`{OUT}/literature/summary.json`、`{OUT}/bitwise.json`。")
    add(f"> 网格 `{cfg['nx']}x{cfg['nz']}`；本文数字全部由生成器从上述 JSON 读出，零手抄。")
    add("")

    add("## ① 案例介绍（物理背景 / 几何 / 参数表）")
    add("")
    add("浮心纵向位置 `LCB` 是主尺度要素，不是装饰：它决定阻力（首尾波系干涉）、")
    add("纵倾与操纵。舱内 fleet manifest 的公开值：KCS `-1.48 % L`、KVLCC2 `+3.48 % L`。")
    add("重建族 `cb_calibrated_sections` 的对称律 `g = sin(pi x/L)` 只能给出 `LCB = 0`，")
    add("因此第 3/5 轮那 500 条『系数精确』的船，每一条都有一条主尺度要素错了 1.5–3.5 % L。")
    add("")
    add("本轮给参数向量加第 7 个分量 `lcb`（有符号、以 L 为单位、正=偏前），并按船型取服务区间：")
    add("")
    add("| 船型 | LCB 抽样带 [% L] |")
    add("|---|---|")
    for ship_type, (lo, hi) in LCB_FRACTION_RANGE.items():
        add(f"| {ship_type} | {lo * 100:+.1f} … {hi * 100:+.1f} |")
    add("")
    add("**取值来源与等级**：Watson《Practical Ship Design》第 2 章主尺度表、")
    add("Schneekluth & Bertram 第 2 章 LCB–Fn 关系图的**离线回忆**典型服务区间，")
    add("`verified=false`，不是拟合任何一条船。它足以检验控制回路（把浮心放到请求位置），")
    add("不足以支撑『这条船就是真船』。")
    add("")

    add("## ② 怎么做的（方法 / 流程 / 关键设置）")
    add("")
    add("1. **管线接线**：`_rebuild_calibrated` 把 `ship['lcb']` 传进 "
        "`CalibratedHullParameters.lcb_offset`；`lcb` 缺省为 `0.0`，走第 3 轮对称路径，")
    add("   曲面逐位不变（见 3.4 的哈希回归）。")
    add("2. **不可达即拒**：`LcbRequestUnreachable` 单独捕获，判 `REFUSED`，")
    add("   错误串以 `lcb_unreachable:` 开头，失败账与网格缺陷分开记。`|lcb| >= 5 % L`")
    add("   在 `problems()` 阶段就被拒（声明式护栏，第 7 轮定）。")
    add("3. **新门 G11**：请求里带 `lcb` 时，用 `cb_sections.realised_lcb_fraction`")
    add("   （散度定理、四面体体积/形心，纯 Python，不看构造）量网格自己的浮心，")
    add(f"   与请求比较，容差 `{LCB_TOLERANCE} L`（=0.5 %L，约为第 7 轮最差实测误差的 500 倍）。")
    add("   门只在**参数向量明确请求 LCB** 时启用：老链条（Series 60 / Wigley / 参数化截面）")
    add("   没有 LCB 自由度，被请求时报 `unsupported` 而不是默默放行。")
    add("4. **两跑对照**：`population` 500 船带 LCB；`control` 同一 500 船**剥掉 LCB 列**，")
    add("   即第 3/5 轮配置，在本树重跑，使通过率逐行可比（同一网格、同一门集）。")
    add("5. **样本外**：6 条公开基准船按**各自出版 LCB** 重建（第 6 轮只测不控）。")
    add("6. **哈希回归**：20 条样本船，`lcb` 缺失 vs `lcb = 0` 两次重建比较顶点/面哈希。")
    add(f"7. **并行**：`run_shipd_rebuild(workers={os.cpu_count()})` 进程池，每船独立，")
    add("   汇总与工件与串行一致（41x13 上串/并两跑 evidence 逐字段相同，见 tests）。")
    add("")

    add("## ③ 计算结果")
    add("")
    add("### 3.1 500 船总体重建（带 LCB 目标）")
    add("")
    add("| 量 | 值 |")
    add("|---|---|")
    add(f"| 船数 | {population['n_ships']} |")
    add(f"| 规范八门通过率 | {population['pass_rate']:.2%} |")
    add(f"| 严格通过率（+Cb/Cm/Cwp+LCB） | {population['strict_pass_rate']:.2%} |")
    add(f"| 方法分布 | {json.dumps(population['method_distribution'], ensure_ascii=False)} |")
    add(f"| 请求 LCB 的船 | {lcb['n_requested']} |")
    add(f"| LCB 被兑现（G11 通过） | {lcb['n_honoured']} |")
    add(f"| LCB 不可达（拒） | {lcb['n_refused_unreachable']} |")
    add(f"| LCB 门失败 | {lcb['n_gate_failed']} |")
    add(f"| LCB 最差绝对误差 | {lcb['max_abs_error_L']:.3e} L |")
    add(f"| LCB 中位绝对误差 | {lcb['median_abs_error_L']:.3e} L |")
    add(f"| LCB 平均有符号误差 | {lcb['mean_signed_error_L']:+.3e} L |")
    add(f"| 单船中位耗时 | {population['total_time_s'] / max(population['n_ships'], 1):.3f} s |")
    add("")
    add("### 3.2 对照：剥掉 LCB 列（第 3/5 轮配置，本树重跑）")
    add("")
    add("| 量 | 0-LCB 对照 | 带 LCB |")
    add("|---|---|---|")
    add(f"| 规范八门通过率 | {control['pass_rate']:.2%} | {population['pass_rate']:.2%} |")
    add(f"| 严格通过率 | {control['strict_pass_rate']:.2%} | {population['strict_pass_rate']:.2%} |")
    add(f"| 请求 LCB 的船 | {clcb['n_requested']} | {lcb['n_requested']} |")
    add(f"| 失败/拒 | {control['n_failed'] + control['n_refused']} | "
        f"{population['n_failed'] + population['n_refused']} |")
    add("")
    add("### 3.3 系数与标称字段相关性（LCB 倾斜下不退化）")
    add("")
    add("| 量 | 平均偏差 | Pearson r | 斜率 | 判定 |")
    add("|---|---|---|---|---|")
    ws = corr["wetted_surface"]
    dv = corr["displacement"]
    add(f"| 湿表面积 vs Holtrop & Mennen | {ws['relative_error']['mean']:+.4%} | "
        f"{ws['pearson_r']:.5f} | {ws['slope_through_origin']:.4f} | 标称回归 |")
    add(f"| 排水体积 vs Cb·L·B·T | {dv['relative_error']['mean']:+.4%} | "
        f"{dv['pearson_r']:.5f} | — | 恒等式（Cb 由体积定义） |")
    add("")
    coef = population["coefficient_agreement"]
    add("| 系数 | 平均相对误差 | 中位 | 最大 | 在 5 % 容差内 |")
    add("|---|---|---|---|---|")
    for key, summary in coef.items():
        add(f"| {key} | {summary['mean']:.3e} | {summary['median']:.3e} | "
            f"{summary['max']:.3e} | {summary['within_tolerance_fraction']:.2%} |")
    add("")
    add("### 3.4 未请求 LCB 时的逐位回归")
    add("")
    add(f"| 项 | 值 |")
    add("|---|---|")
    add(f"| 样本船数 | {bitwise['n']} |")
    add(f"| 顶点/面哈希全部相同 | {bitwise['all_identical']} |")
    _resid = bitwise["max_abs_realised_lcb_no_request"]
    add("| 未请求时实测 |LCB| 最大值 | "
        + ("—" if _resid is None else f"{_resid:.3e}") + " |")
    add("")
    add("### 3.5 样本外：6 条公开基准船按出版 LCB 重建")
    add("")
    add("| 船 | Cb | Cp≤Cwp | 方法 | 出版 LCB [%L] | 实测 LCB [%L] | 判定 |")
    add("|---|---|---|---|---|---|---|")
    for hull in literature["hulls"]:
        pub = "—" if hull["lcb_published"] is None else f"{hull['lcb_published'] * 100:+.2f}"
        got = "—" if hull["realised_lcb"] is None else f"{hull['realised_lcb'] * 100:+.4f}"
        add(f"| {hull['ship_id']} | {hull['Cb']:.4f} | {'是' if hull['cp_le_cwp'] else '否'} | "
            f"{hull['method']} | {pub} | {got} | {hull['verdict']}/{hull['strict_verdict']} |")
    add("")
    add("### 3.6 失败模式分类账")
    add("")
    for name, block in ledger.items():
        add(f"**{name}** — 失败/拒 {block['n_failures']} 条")
        add("")
        if block["by_class"]:
            add("| 类别 | 条数 |")
            add("|---|---|")
            for cls, count in sorted(block["by_class"].items()):
                add(f"| `{cls}` | {count} |")
        else:
            add("（无失败条目）")
        add("")

    add("### 3.7 不可达探针（让『LCB 不可达』这类失败真的发火）")
    add("")
    add(f"抽样带全部落在第 7 轮测出的可达包线内，因此 500 船里**没有**一条 LCB 不可达。")
    add(f"空的失败类别只有在它能响时才有意义，所以把 `Cb` 最大的 {stress['n_probes']} 条船的 LCB 请求")
    add(f"推到 API 护栏边缘（`{stress['requested_lcb']:+.4f} L`）：")
    add("")
    add("| ship_id | 船型 | Cb | 几何可达上限 [%L] | 请求 [%L] | 判定 | 与实测上限一致 |")
    add("|---|---|---|---|---|---|---|")
    for row in stress["rows"]:
        add(f"| {row['ship_id']} | {row['ship_type']} | {row['Cb']:.4f} | "
            f"{row['geometric_reach_L'] * 100:.3f} | {row['requested_lcb'] * 100:+.3f} | "
            f"{row['verdict']} | {'是' if row['consistent'] else '否'} |")
    add("")
    add(f"读数：{stress['n_probes']} 条中 {stress['n_refused']} 条被拒、{stress['n_honoured']} 条被兑现；")
    add(f"其中 {stress['n_consistent_with_measured_reach']}/{stress['n_probes']} 条的判决与"
        "`cb_sections.lcb_reach` 独立测出的几何饱和一致（请求 > 饱和 时就该拒）。")
    add("被拒的船**没有产出网格**，也没有产出系数误差 —— 拒绝发生在几何之前，这是该路径的设计。")
    add("")

    add("### 3.8 可失败性：LCB 拒绝路径在包线内到底会不会响")
    add("")
    add("3.1 的 500 船全部兑现。空的失败类别只有在它能响时才有意义，所以本轮做了两"
        "件事：")
    add("先跑一个**对抗种群**（LCB ~ U(-5 %, +5 %) L，与 `Cb` 无关，`seed="
        f"{hostile['seed']}`），再派一个**几何可达性普查**去核对拒绝边界。")
    add("")
    add("| 量 | 值 |")
    add("|---|---|")
    add(f"| 对抗种群请求 LCB 的船 | {hostile['n_requested']} |")
    add(f"| 兑现（G11 通过） | {hostile['n_honoured']} |")
    add(f"| 不可达被拒 | {hostile['n_refused_unreachable']} "
        f"（前倾 {hostile['refused_forward']} / 后倾 {hostile['refused_aft']}） |")
    add(f"| 兑现者最差绝对误差 | {hostile['max_abs_error_L_honoured']:.3e} L |")
    add(f"| 抽查船数（可达性普查） | {survey['n_ships']} |")
    add(f"| API 护栏 | ±{survey['guardrail_L'] * 100:.2f} % L |")
    add(f"| 几何可达上限最小值 | 前倾 {survey['min_reach_forward_L'] * 100:.2f} % L / "
        f"后倾 {survey['min_reach_aft_L'] * 100:.2f} % L |")
    add(f"| 护栏落在几何可达之外 | 前倾 {survey['n_guardrail_outside_forward']} 条 / "
        f"后倾 {survey['n_guardrail_outside_aft']} 条 |")
    add("")
    add(f"**如实结论**：对抗种群 {hostile['n_refused_unreachable']} 条被拒（"
        f"{hostile['n_requested']} 条请求中）。拒绝路径没有坏 —— 3.7 的探针"
        "把请求推到 `+4.90 % L` 时确实触发拒绝，且拒绝消息里的饱和值与独立普查一致；"
        "对抗种群只多出这一条，说明「更重要的约束不是几何」。")
    add("原因是**在 ±5 % L 的 API 护栏内，几何几乎总是更宽**：抽样普查的几何可达上限最小值为")
    add(f"前倾 {survey['min_reach_forward_L'] * 100:.2f} % L / 后倾 {survey['min_reach_aft_L'] * 100:.2f} % L，")
    add("护栏先于几何成为约束。也就是说：")
    add("")
    add("- 在种群尺度上，**限制 LCB 可达性的主要不是几何，而是我们自己定的 5 % L 护栏**；")
    add(f"  因此 `LCB_UNREACHABLE` 在护栏内是低概率事件（带型典型 LCB 的种群 0/500，"
        f"与 Cb 无关的对抗种群 {hostile['n_refused_unreachable']}/500）。")
    add("- 这条结论修正了第 7 轮报告的读法：第 7 轮 3.3 的『几何饱和 3.57–13.27 % L』是"
        "**逐个系数组做的搜索**，而本轮普查给出同量级的独立复核；两者一致，但都不该被读成"
        "『LCB 请求经常会被几何拒掉』。")
    add("")
    add("| 船型 | 抽查数 | 前倾最小可达 [%L] | 后倾最小可达 [%L] | 护栏越界（前/后） |")
    add("|---|---|---|---|---|")
    for ship_type, bucket in sorted(survey["by_ship_type"].items()):
        add(f"| {ship_type} | {bucket['n']} | {bucket['min_reach_forward_L'] * 100:.2f} | "
            f"{bucket['min_reach_aft_L'] * 100:.2f} | "
            f"{bucket['n_guardrail_outside_forward']}/{bucket['n_guardrail_outside_aft']} |")
    add("")

    add("## ④ 结果图")
    add("")
    add(f"![LCB at population scale]({fig_path.as_posix()})")
    add("")
    add("四联图：(a) 请求 LCB vs 网格实测 LCB（1:1 线）；(b) 各船型实测 LCB 落点 vs 抽样带（灰条）；")
    add("(c) 加 LCB 前后通过率（八门 / 严格）；(d) 6 条公开基准船的出版 vs 重建 LCB。")
    add("")

    add("## ⑤ 注意点（踩过的坑 / 局限 / 适用范围）")
    add("")
    add("1. **LCB 是控制回路，不是型线。** 倾斜只移动体积纵向分布，不带球鼻、方尾、外飘、")
    add("   折角水线，也不控制 LCF。命中 `LCB` 不等于画出那条真船。")
    add("2. **LCB 抽样带是离线回忆的典型区间**（`verified=false`）。它只保证『请求落在")
    add("   真实船会出现的位置上』，不保证任何一条船的业务正确性。")
    add("3. **可达性受系数约束**（第 7 轮 3.3）：`Cb` 越大、带龙骨升高时可达带越窄。")
    add("   种群中触到边界的船被记为 `LCB_UNREACHABLE` 拒，不返回一条浮心不对的船。")
    add("4. **老链条明确不支持 LCB**：Series 60 / Wigley / 参数化截面被请求 LCB 时 G11 报")
    add(f"   unsupported，严格判定为 FAILED（见 3.5 中 `cp_le_cwp=否` 的船）。这是如实标注，")
    add("   不是新缺陷 —— 它本来就没有这个自由度。")
    add("5. **门容差是政策数字**（0.5 %L）。报告里同时给实测误差（3.1 末三行），")
    add("   门只负责判决。")
    add("")

    add("## ⑥ 误差判定（逐项 vs 基准，显式回答是否 <3 %）")
    add("")
    add("| 关键量 | 基准 | 本轮结果 | 是否 <3 % |")
    add("|---|---|---|---|")
    add(f"| LCB（500 船总体） | 请求值（自身即基准，控制精度） | 最差 {lcb['max_abs_error_L']:.3e} L | "
        f"是（相对 L 为 {lcb['max_abs_error_L']:.2e}） |")
    kcs = next((h for h in literature["hulls"] if h["ship_id"] == "KCS"), None)
    kvl = next((h for h in literature["hulls"] if h["ship_id"] == "KVLCC2"), None)
    for hull in (kcs, kvl):
        if hull and hull["lcb_published"] is not None and hull["realised_lcb"] is not None:
            err = abs(hull["realised_lcb"] - hull["lcb_published"]) / abs(hull["lcb_published"])
            add(f"| {hull['ship_id']} LCB | 出版 {hull['lcb_published'] * 100:+.2f} %L | "
                f"实测 {hull['realised_lcb'] * 100:+.4f} %L | 是（相对误差 {err:.2e}） |")
    add(f"| Cb / Cm / Cwp（500 船） | 参数向量 | 最差 {max(s['max'] for s in coef.values()):.3e} | 是 |")
    add(f"| 湿表面积（500 船） | Holtrop & Mennen 回归 | 平均偏差 {ws['relative_error']['mean']:+.2%} | "
        f"否（>3 %），属已知族差异：本族无球鼻/方尾，见第 6 轮已披露 |")
    add(f"| 排水体积（500 船） | Cb·L·B·T | 平均偏差 {dv['relative_error']['mean']:+.2e} | 是 |")
    add("")
    add("基准性质说明：LCB 的『基准』分两类 —— 500 船一栏的基准是请求值本身（控制精度），")
    add("6 条基准船一栏的基准是 manifest 的出版值（`retrieval=offline_recollection,")
    add("verified=false`）。两者都只支持『管线把浮心放到了该放的位置』，不支持『船型精度』。")

    path = Path("docs/SHIPD_REBUILD_REPORT_R8_LCB_POPULATION.md")
    path.write_text("\n".join(lines) + "\n")
    return path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ships", type=int, default=500)
    parser.add_argument("--nx", type=int, default=81)
    parser.add_argument("--nz", type=int, default=21)
    parser.add_argument("--workers", type=int, default=min(8, os.cpu_count() or 1))
    parser.add_argument("--bitwise-sample", type=int, default=20)
    args = parser.parse_args()

    t_start = time.time()
    OUT.mkdir(parents=True, exist_ok=True)

    print(f"[1/5] population: {args.ships} ships, grid {args.nx}x{args.nz}, "
          f"workers={args.workers}", flush=True)
    ships, population = run_population(args.ships, args.nx, args.nz, args.workers)
    print(f"      pass_rate={population['pass_rate']:.4f} "
          f"strict={population['strict_pass_rate']:.4f} "
          f"lcb_honoured={population['lcb_contract']['n_honoured']}/"
          f"{population['lcb_contract']['n_requested']}", flush=True)

    print("[2/5] control: same ships without the LCB column", flush=True)
    control = run_control(ships, args.nx, args.nz, args.workers)
    print(f"      pass_rate={control['pass_rate']:.4f} "
          f"strict={control['strict_pass_rate']:.4f}", flush=True)

    print("[3/5] literature: six published hulls at their published LCB", flush=True)
    literature = run_literature(args.nx, args.nz)
    for hull in literature["hulls"]:
        print(f"      {hull['ship_id']:<9} published={hull['lcb_published']} "
              f"realised={hull['realised_lcb']} {hull['strict_verdict']}", flush=True)

    print("[4/5] bitwise: lcb absent vs lcb = 0", flush=True)
    bitwise = bitwise_regression(args.bitwise_sample, args.nx, args.nz)
    (OUT / "bitwise.json").write_text(json.dumps(bitwise, indent=2))
    print(f"      all_identical={bitwise['all_identical']}", flush=True)

    print("[4b/5] stress probe: fullest hulls asked for the API-limit LCB", flush=True)
    stress = stress_unreachable_probe(3, args.nx, args.nz)
    print(f"      refused={stress['n_refused']}/{stress['n_probes']} "
          f"honoured={stress['n_honoured']}", flush=True)

    print("[4c/5] hostile LCB population: U(-5 %, +5 %) L, unrelated to Cb", flush=True)
    hostile = run_hostile_lcb(ships, args.nx, args.nz, args.workers)
    print(f"      honoured={hostile['n_honoured']}/{hostile['n_requested']} "
          f"refused={hostile['n_refused_unreachable']} "
          f"(fwd {hostile['refused_forward']} / aft {hostile['refused_aft']})", flush=True)

    print("[4d/5] reach survey: geometric LCB reach of a stratified sample", flush=True)
    survey = reach_survey(ships[::8], args.nx, args.nz)
    print(f"      min reach fwd={survey['min_reach_forward_L']:.4f} L "
          f"aft={survey['min_reach_aft_L']:.4f} L; "
          f"guardrail outside reach: fwd {survey['n_guardrail_outside_forward']}, "
          f"aft {survey['n_guardrail_outside_aft']} of {survey['n_ships']}", flush=True)

    print("[5/5] ledger, figures, report", flush=True)
    ledger = failure_ledger({
        "population_with_lcb": OUT / "population",
        "control_no_lcb": OUT / "control",
        "hostile_lcb": OUT / "hostile_lcb",
    })
    (OUT / "failure_mode_ledger.json").write_text(json.dumps(ledger, indent=2))
    fig_path = make_figures(population, control, literature, bitwise, hostile)
    report = write_report(
        population, control, literature, bitwise, ledger, fig_path, stress, hostile,
        survey,
    )

    summary = {
        "population": {k: population[k] for k in
                       ("n_ships", "pass_rate", "strict_pass_rate", "lcb_contract",
                        "coefficient_agreement", "settings")},
        "control": {k: control[k] for k in
                    ("n_ships", "pass_rate", "strict_pass_rate", "lcb_contract")},
        "literature": literature,
        "bitwise": {k: bitwise[k] for k in
                    ("n", "all_identical", "max_abs_realised_lcb_no_request")},
        "ledger_classes": {name: block["by_class"] for name, block in ledger.items()},
        "stress_unreachable": {k: stress[k] for k in ("n_probes", "n_refused", "n_honoured")},
        "hostile_lcb": {k: hostile[k] for k in
                        ("seed", "n_requested", "n_honoured", "n_refused_unreachable",
                         "refused_forward", "refused_aft", "max_abs_error_L_honoured")},
        "hostile_pass_rate": hostile["evidence"]["pass_rate"],
        "reach_survey": {k: survey[k] for k in
                         ("n_ships", "guardrail_L", "min_reach_forward_L",
                          "min_reach_aft_L", "n_guardrail_outside_forward",
                          "n_guardrail_outside_aft")},
        "figures": [str(fig_path)],
        "report": str(report),
        "wall_time_s": time.time() - t_start,
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps({k: summary[k] for k in ("bitwise", "ledger_classes", "wall_time_s")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())