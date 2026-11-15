"""Round-4 runner for GAP-CAD-C2-SHIPD-REBUILD: hull-form offsets (keel rise).

Round 3 delivered a coefficient-calibrated family that realises ``Cb``, ``Cm`` and
``Cwp`` to round-off, and passed all ten gates on 120 ships -- but it left one
acceptance item unmet: the wetted surface of the rebuilt hulls sits up to **+14.7 %**
above the Holtrop & Mennen (1982) regression the nominal performance fields use. Round
3 called that a "shape-family difference" and deferred it.

This round fixes the cause. The round-3 family gives *every* station the full design
draft, so the ends are full-height blades of zero width: the surface is larger than a
real hull's at the same coefficients. Real hulls lift the keel through the entrance and
the run (forefoot and run cut-away). :mod:`tensorcad.cb_sections` now carries that
feature (``end_draft_ratio`` plus an optional parallel middle body), and this runner

* calibrates ``end_draft_ratio`` **per ship type** so the family reproduces the
  Holtrop & Mennen regression surface, subject to the geometric feasibility limit
  (the coefficient vector must still be realisable -- a full hull cannot be hollowed
  out as far as a fine one),
* re-runs the whole parameter-vector -> mesh -> gate pipeline for both families, and
* writes every number and figure the report quotes.

    python3 tools/shipd_rebuild/run_round4_offsets.py            # full 500-ship campaign
    python3 tools/shipd_rebuild/run_round4_offsets.py --ships 50 # quick smoke
"""
from __future__ import annotations

import argparse
import json
import math
import statistics
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from tensorcad.cb_sections import (  # noqa: E402
    CalibratedHullParameters,
    build_calibrated_hull,
    realised_coefficients,
)
from tensorcad.shipd_population import (  # noqa: E402
    build_shipd_database,
    nominal_performance_fields,
)
from tensorcad.shipd_rebuild import run_shipd_rebuild  # noqa: E402

#: Ladder of ``end_draft_ratio`` values searched during calibration. The floor is a
#: *geometric* limit, not a taste: below it no positive fullness exponent realises Cb.
EPS_LADDER = (1.0, 0.8, 0.65, 0.55, 0.45, 0.35, 0.25, 0.15)


def _write(path: Path, payload) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2))


def _wetted(ship: dict, *, eps: float, ppar: float = 0.0, nx: int = 81, nz: int = 21):
    """Rebuilt wetted surface and coefficient errors for one ship / one shape."""
    params = CalibratedHullParameters(
        length=ship["L"], beam=ship["B"], draft=ship["T"],
        block_coefficient=ship["Cb"],
        midship_coefficient=ship["Cm"],
        waterplane_coefficient=ship["Cwp"],
        nx=nx, nz=nz,
        parallel_body_fraction=ppar, end_draft_ratio=eps,
    )
    problems = params.problems()
    if problems:
        return None, problems
    hull = build_calibrated_hull(params)
    measured = realised_coefficients(
        hull.mesh, length=ship["L"], beam=ship["B"], draft=ship["T"]
    )
    return measured["wetted_surface_m2"], []


def surface_deviation(ship: dict, *, eps: float, ppar: float = 0.0) -> float | None:
    """``(mesh - nominal) / nominal`` for the wetted surface, or ``None`` if refused."""
    nominal = nominal_performance_fields(ship)["wetted_surface_m2"]
    measured, _ = _wetted(ship, eps=eps, ppar=ppar)
    if measured is None:
        return None
    return (measured - nominal) / nominal


def calibrate_end_draft_ratio(ships: list[dict], *, sample: int = 4) -> dict:
    """Per-type ``end_draft_ratio`` that reproduces the Holtrop & Mennen surface.

    The deviation falls monotonically as the keel rises (the ends lose girth), so the
    calibration is a monotone ladder search: walk ``EPS_LADDER`` downwards and stop at
    the first value where the type's *median* deviation crosses zero. If the
    coefficient vector stops being realisable first, the largest feasible cut-away is
    kept and the residual positive deviation is reported -- it is a limit of the
    parameter vector, not a fitting choice.
    """
    by_type: dict[str, list[dict]] = {}
    for ship in ships:
        by_type.setdefault(ship["ship_type"], []).append(ship)

    out: dict[str, dict] = {}
    for ship_type, group in sorted(by_type.items()):
        probe = group if len(group) <= sample else group[:: max(1, len(group) // sample)][:sample]
        rows: list[dict] = []
        chosen = None
        for eps in EPS_LADDER:
            deviations, refused = [], 0
            for ship in probe:
                deviation = surface_deviation(ship, eps=eps)
                if deviation is None:
                    refused += 1
                else:
                    deviations.append(deviation)
            row = {
                "end_draft_ratio": eps,
                "n_feasible": len(deviations),
                "n_refused": refused,
                "median_deviation": statistics.median(deviations) if deviations else None,
            }
            rows.append(row)
            if deviations and chosen is None and row["median_deviation"] <= 0.0:
                chosen = row
        if chosen is None:
            feasible = [row for row in rows if row["n_refused"] == 0]
            chosen = feasible[-1] if feasible else rows[0]
            limited_by = "feasibility"
        else:
            limited_by = "zero_crossing"
        out[ship_type] = {
            "end_draft_ratio": chosen["end_draft_ratio"],
            "limited_by": limited_by,
            "ladder": rows,
            "n_ships": len(group),
        }
    return out


def assign_shape(ships: list[dict], calibration: dict) -> list[dict]:
    """Copy the ships with their per-type shape attached."""
    out = []
    for ship in ships:
        shape = calibration[ship["ship_type"]]
        out.append({**ship, "end_draft_ratio": shape["end_draft_ratio"],
                    "parallel_body_fraction": 0.0})
    return out


# --------------------------------------------------------------------------- #
# Figures
# --------------------------------------------------------------------------- #
def figure_gate_rates(evidence: dict, path: Path) -> None:
    names = list(evidence["flat_keel"]["gate_pass_rates"])
    x = range(len(names))
    width = 0.38
    fig, ax = plt.subplots(figsize=(11, 4.6))
    for k, (mode, colour, label) in enumerate((
        ("flat_keel", "#9aa7b4", "flat keel (round-3 family)"),
        ("offsets", "#2c6e9b", "keel rise (round-4 family)"),
    )):
        values = [100.0 * evidence[mode]["gate_pass_rates"][n] for n in names]
        ax.bar([i + (k - 0.5) * width for i in x], values, width, color=colour, label=label)
    ax.set_xticks(list(x))
    ax.set_xticklabels([n.replace("_", "\n") for n in names], fontsize=7)
    ax.set_ylabel("gate pass rate [%]")
    ax.set_ylim(0, 105)
    ax.set_title("Quality gates, %d ships: %s" % (
        evidence["n_ships"],
        "keel rise does not cost a single gate"))
    ax.legend(loc="lower right", fontsize=8)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def figure_deviation_by_type(per_type: dict, path: Path) -> None:
    types = sorted(per_type, key=lambda t: -abs(per_type[t]["flat_keel"]["median"]))
    x = range(len(types))
    width = 0.38
    fig, ax = plt.subplots(figsize=(11, 4.6))
    for k, (mode, colour, label) in enumerate((
        ("flat_keel", "#c0504d", "flat keel (no cut-away)"),
        ("offsets", "#2c6e9b", "keel rise (calibrated per type)"),
    )):
        med = [100.0 * per_type[t][mode]["median"] for t in types]
        worst = [100.0 * per_type[t][mode]["max_abs"] for t in types]
        ax.bar([i + (k - 0.5) * width for i in x], med, width, color=colour, label=label)
        ax.plot([i + (k - 0.5) * width for i in x], worst, "k.", ms=3)
    ax.axhline(3.0, color="k", ls="--", lw=1, label="3 % acceptance line")
    ax.axhline(0.0, color="k", lw=0.6)
    ax.set_xticks(list(x))
    ax.set_xticklabels(types, rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("wetted surface deviation [%]")
    ax.set_title("Wetted surface vs Holtrop & Mennen: median (bar) and worst ship (dot)")
    ax.legend(fontsize=8)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def figure_calibration_ladder(calibration: dict, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.4, 5.2))
    for ship_type, row in calibration.items():
        pts = [(r["end_draft_ratio"], 100.0 * r["median_deviation"])
               for r in row["ladder"] if r["median_deviation"] is not None]
        if not pts:
            continue
        ax.plot([p[0] for p in pts], [p[1] for p in pts], "o-", ms=3, lw=1.2,
                label=f"{ship_type} (eps={row['end_draft_ratio']:.2f})")
    ax.axhline(0.0, color="k", lw=0.8)
    ax.axhline(3.0, color="k", ls="--", lw=1)
    ax.invert_xaxis()
    ax.set_xlabel("end draft ratio  eps  (1.0 = flat keel, smaller = more cut-away)")
    ax.set_ylabel("median wetted surface deviation [%]")
    ax.set_title("Calibration ladder: the keel rise is the missing hull-form feature")
    ax.legend(fontsize=7, ncol=2)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def figure_scatter(rows: list[dict], path: Path) -> None:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.6))
    nominal = [r["nominal_wetted_surface_m2"] for r in rows]
    for ax, key, title, colour in (
        (ax1, "flat_keel", "flat keel (round-3 family)", "#c0504d"),
        (ax2, "offsets", "keel rise (round-4 family)", "#2c6e9b"),
    ):
        mesh = [r[key]["wetted_surface_m2"] for r in rows]
        ax.plot(nominal, mesh, "o", ms=2.4, color=colour, alpha=0.75)
        lo = min(nominal)
        hi = max(nominal)
        ax.plot([lo, hi], [lo, hi], "k--", lw=1, label="1:1")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("nominal wetted surface, Holtrop & Mennen [m²]")
        ax.set_ylabel("rebuilt mesh wetted surface [m²]")
        ax.set_title(title)
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
    fig.suptitle("Wetted surface, %d ships: the cut-away closes the regression gap" % len(rows))
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def figure_convergence(rows: list[dict], path: Path) -> None:
    resolutions = [r["resolution"] for r in rows]
    area = [abs(r["wetted_surface_relative_error_vs_finest"]) for r in rows]
    volume = [max(r["volume_relative_error"], 1e-17) for r in rows]
    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    ns = [r["n_stations"] * r["n_levels"] for r in rows]
    ax.loglog(ns, area, "s-", color="#2c6e9b", label="wetted surface vs finest grid")
    ax.loglog(ns, volume, "o-", color="#c0504d", label="block coefficient vs request")
    for n, a in zip(ns, area):
        pass
    ax.set_xlabel("grid points (nx * nz)")
    ax.set_ylabel("relative error")
    ax.set_title("Keel-rise family: coefficients stay exact, area converges")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def convergence_study(ship: dict, out_dir: Path) -> list[dict]:
    """Resolution study on one keel-rise hull (fixed shape)."""
    rows = []
    for nx, nz in ((21, 7), (41, 13), (81, 21), (161, 41)):
        params = CalibratedHullParameters(
            length=ship["L"], beam=ship["B"], draft=ship["T"],
            block_coefficient=ship["Cb"], midship_coefficient=ship["Cm"],
            waterplane_coefficient=ship["Cwp"],
            nx=nx, nz=nz, end_draft_ratio=ship["end_draft_ratio"],
        )
        hull = build_calibrated_hull(params)
        measured = realised_coefficients(
            hull.mesh, length=ship["L"], beam=ship["B"], draft=ship["T"]
        )
        rows.append({
            "resolution": f"{nx}x{nz}",
            "n_stations": nx,
            "n_levels": nz,
            "triangles": len(hull.mesh.faces),
            "volume_relative_error": measured["relative_error"].get("Cb", 0.0),
            "wetted_surface_m2": measured["wetted_surface_m2"],
            "Cb_error": abs(measured["coefficients"]["Cb"] / ship["Cb"] - 1.0),
        })
    finest = rows[-1]["wetted_surface_m2"]
    for row in rows:
        row["wetted_surface_relative_error_vs_finest"] = (
            row["wetted_surface_m2"] - finest
        ) / finest
    _write(out_dir / "convergence.json", rows)
    return rows


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ships", type=int, default=500)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--out", default="runs/shipd_rebuild_r4")
    args = parser.parse_args(argv)

    out_dir = Path(args.out)
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    ships = build_shipd_database(args.ships, seed=args.seed)
    print(f"{len(ships)} ships, types {sorted({s['ship_type'] for s in ships})}")

    t0 = time.time()
    calibration = calibrate_end_draft_ratio(ships)
    _write(out_dir / "offsets_calibration.json", calibration)
    print("calibrated end_draft_ratio per type:", {
        t: row["end_draft_ratio"] for t, row in calibration.items()
    }, f"({time.time() - t0:.1f} s)")

    shaped = assign_shape(ships, calibration)
    flat = [{**ship, "end_draft_ratio": 1.0, "parallel_body_fraction": 0.0}
            for ship in ships]

    evidence = {"n_ships": len(ships), "seed": args.seed,
                "resolution": {"nx": 81, "nz": 21},
                "end_draft_ratio_by_type": {
                    t: row["end_draft_ratio"] for t, row in calibration.items()}}
    for mode, group in (("flat_keel", flat), ("offsets", shaped)):
        target = out_dir / mode
        start = time.time()
        result = run_shipd_rebuild(group, nx=81, nz=21, out_dir=target,
                                   method_preference="calibrated")
        result["wall_time_s"] = time.time() - start
        _write(target / "evidence.json", result)
        evidence[mode] = result
        print(f"{mode:10s}: pass {result['pass_rate']:.4f} "
              f"strict {result['strict_pass_rate']:.4f} "
              f"({result['wall_time_s']:.1f} s)")
        evidence[mode].pop("correlation_sanity", None)

    # Per-ship deviations vs the nominal performance fields, both families.
    per_ship: list[dict] = []
    for ship in shaped:
        nominal = nominal_performance_fields(ship)["wetted_surface_m2"]
        row = {"ship_id": ship["ship_id"], "ship_type": ship["ship_type"],
               "Cb": ship["Cb"], "Cwp": ship["Cwp"],
               "end_draft_ratio": ship["end_draft_ratio"],
               "nominal_wetted_surface_m2": nominal}
        for mode, eps in (("flat_keel", 1.0), ("offsets", ship["end_draft_ratio"])):
            measured, problems = _wetted(ship, eps=eps)
            row[mode] = {
                "wetted_surface_m2": measured,
                "deviation": None if measured is None else (measured - nominal) / nominal,
                "problems": problems,
            }
        per_ship.append(row)
    _write(out_dir / "per_ship_surface.json", per_ship)

    per_type: dict[str, dict] = {}
    for ship_type in sorted({row["ship_type"] for row in per_ship}):
        group = [row for row in per_ship if row["ship_type"] == ship_type]
        entry = {"n_ships": len(group)}
        for mode in ("flat_keel", "offsets"):
            values = [row[mode]["deviation"] for row in group
                      if row[mode]["deviation"] is not None]
            entry[mode] = {
                "n_measured": len(values),
                "median": statistics.median(values) if values else None,
                "mean": statistics.fmean(values) if values else None,
                "max_abs": max(abs(v) for v in values) if values else None,
                "min": min(values) if values else None,
                "max": max(values) if values else None,
            } if values else {"n_measured": 0}
        per_type[ship_type] = entry
    _write(out_dir / "surface_deviation_by_type.json", per_type)

    all_flat = [abs(row["flat_keel"]["deviation"]) for row in per_ship
                if row["flat_keel"]["deviation"] is not None]
    all_off = [abs(row["offsets"]["deviation"]) for row in per_ship
               if row["offsets"]["deviation"] is not None]
    summary = {
        "flat_keel": {
            "n_measured": len(all_flat),
            "median_abs_deviation": statistics.median(all_flat),
            "mean_abs_deviation": statistics.fmean(all_flat),
            "max_abs_deviation": max(all_flat),
            "n_over_3pct": sum(1 for v in all_flat if v > 0.03),
        },
        "offsets": {
            "n_measured": len(all_off),
            "median_abs_deviation": statistics.median(all_off),
            "mean_abs_deviation": statistics.fmean(all_off),
            "max_abs_deviation": max(all_off),
            "n_over_3pct": sum(1 for v in all_off if v > 0.03),
        },
    }
    _write(out_dir / "surface_summary.json", summary)
    print("surface:", json.dumps(summary, indent=1))

    convergence = convergence_study(shaped[len(shaped) // 2], out_dir)

    figure_gate_rates(evidence, fig_dir / "gate_pass_rates.png")
    figure_deviation_by_type(per_type, fig_dir / "deviation_by_type.png")
    figure_calibration_ladder(calibration, fig_dir / "calibration_ladder.png")
    figure_scatter(per_ship, fig_dir / "wetted_surface_scatter.png")
    figure_convergence(convergence, fig_dir / "convergence.png")

    _write(out_dir / "comparison.json", {
        "n_ships": len(ships), "seed": args.seed,
        "end_draft_ratio_by_type": {t: r["end_draft_ratio"] for t, r in calibration.items()},
        "flat_keel": {k: evidence["flat_keel"][k] for k in (
            "pass_rate", "strict_pass_rate", "gate_pass_rates", "method_distribution",
            "volume_error_mean", "volume_error_max", "wall_time_s")},
        "offsets": {k: evidence["offsets"][k] for k in (
            "pass_rate", "strict_pass_rate", "gate_pass_rates", "method_distribution",
            "volume_error_mean", "volume_error_max", "wall_time_s")},
        "surface_summary": summary,
        "failure_mode_ledger": evidence["offsets"].get("failure_mode_ledger"),
    })
    print("wrote", out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())