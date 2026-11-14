"""Evidence build for the buttock-branch fix (GAP-CAD-C3-LINES-FAIRING, mission m21077).

Recomputes the *previous* third-family metric in place (the historic contour is reproduced
here verbatim, because the code that produced the old numbers no longer exists) and
compares it with the branch-tracing contour that now ships.  Writes the figures and a
JSON of every number the report quotes, so nothing is transcribed by hand.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from tensorcad.lines_fairing import (
    buttock_coverage,
    buttock_curves,
    buttock_levels,
    polyline_sign_flips,
    perturb_table,
    fair_offset_table,
    _grid_and_flat,
)
from tensorcad.offset_forms import offset_table
from tensorcad.offset_loft import OffsetProvenance, ShipOffsetTable

OUT = Path("docs/lines_fairing")
FORMS = ("wigley_iii", "series60_cb060", "kvlcc2_class")
OLD_FRACTIONS = (0.2, 0.35, 0.5, 0.65, 0.8)


def old_buttock_curves(eta, xi, zeta, levels):
    """The historic contour, verbatim: endpoint bracketing + bisection on a monotone section."""
    curves = []
    for level in levels:
        points = []
        for i, x in enumerate(xi):
            profile = [eta[i][j] for j in range(len(zeta))]
            if level < profile[0] or level > profile[-1]:
                continue
            low, high = 0, len(zeta) - 1
            if profile[low] == level:
                points.append((x, zeta[low]))
                continue
            while high - low > 1:
                middle = (low + high) // 2
                if profile[middle] <= level:
                    low = middle
                else:
                    high = middle
            if profile[high] == profile[low]:
                continue
            weight = (level - profile[low]) / (profile[high] - profile[low])
            points.append((x, zeta[low] + weight * (zeta[high] - zeta[low])))
        if len(points) >= 3:
            curves.append(points)
    return curves


def tumblehome_table(stations=21, waterlines=17):
    length, beam, draft = 10.0, 1.0, 0.65
    xs = [length * i / (stations - 1) for i in range(stations)]
    zs = [draft * j / (waterlines - 1) for j in range(waterlines)]

    def shape(z):
        zeta = z / draft
        return (zeta / 0.65) ** 0.7 if zeta <= 0.65 else 1.0 - 0.55 * (zeta - 0.65) / 0.35

    def plan(x):
        return 0.5 * beam * (1.0 - (2.0 * x / length - 1.0) ** 6)

    rows = tuple(tuple(max(0.0, plan(x) * shape(z)) for z in zs) for x in xs)
    return ShipOffsetTable(
        name="tumblehome test hull",
        length=length,
        stations=tuple(xs),
        waterlines=tuple(zs),
        half_breadths=rows,
        design_draft=draft,
        provenance=OffsetProvenance(
            source="synthetic: half-breadth plan with tumblehome above z/T = 0.65",
            method="y(x, z) = 0.5 B (1 - (2x/L - 1)^6) * shape(z/T)",
        ),
    )


def main() -> None:
    evidence = {"forms": {}, "tumblehome": {}}
    for key in FORMS:
        table = offset_table(key, n_stations=21, n_waterlines=17)
        grid, _ = _grid_and_flat(table)
        peak = max(max(row) for row in grid.eta)
        old_levels = tuple(peak * fraction for fraction in OLD_FRACTIONS)
        old_curves = old_buttock_curves(grid.eta, grid.xi, grid.zeta, old_levels)
        new_levels = buttock_levels(grid.eta, grid.xi, grid.zeta)
        new_curves = buttock_curves(grid.eta, grid.xi, grid.zeta, levels=new_levels)
        coverage = buttock_coverage(grid.eta, grid.xi, grid.zeta, levels=new_levels)
        noisy = perturb_table(table, amplitude=0.004 * 0.5 * table.beam, seed=20261114)
        noisy_grid, _ = _grid_and_flat(noisy)
        old_unfair = old_buttock_curves(noisy_grid.eta, noisy_grid.xi, noisy_grid.zeta, old_levels)
        new_unfair = buttock_curves(
            noisy_grid.eta, noisy_grid.xi, noisy_grid.zeta,
            levels=buttock_levels(noisy_grid.eta, noisy_grid.xi, noisy_grid.zeta),
        )
        result = fair_offset_table(
            noisy, tolerance=0.016 * 0.5 * table.beam, iterations=3, reference=table
        )
        fair_grid, _ = _grid_and_flat(result.table_after)
        new_faired = buttock_curves(
            fair_grid.eta, fair_grid.xi, fair_grid.zeta,
            levels=buttock_levels(noisy_grid.eta, noisy_grid.xi, noisy_grid.zeta),
        )
        evidence["forms"][key] = {
            "old_levels_fraction_of_peak": OLD_FRACTIONS,
            "old_curves": len(old_curves),
            "old_samples": sorted(len(curve) for curve in old_curves),
            "old_flips_clean": sum(polyline_sign_flips(curve) for curve in old_curves),
            "old_flips_noisy": sum(polyline_sign_flips(curve) for curve in old_unfair),
            "new_levels_fraction_of_peak": [round(level / peak, 4) for level in new_levels],
            "new_curves": len(new_curves),
            "new_samples": sorted(len(curve) for curve in new_curves),
            "new_flips_clean": sum(polyline_sign_flips(curve) for curve in new_curves),
            "new_flips_noisy": sum(polyline_sign_flips(curve) for curve in new_unfair),
            "new_flips_faired": sum(polyline_sign_flips(curve) for curve in new_faired),
            "fairing_deviation_m": result.deviation_max,
            "fairing_tolerance_m": result.tolerance,
            "fairing_accepted": result.accepted(),
            "total_flips": [result.flips_before, result.flips_after],
            "coverage": coverage,
        }

    # --- tumblehome: the case the old contour could not see at all
    hull = tumblehome_table()
    grid, _ = _grid_and_flat(hull)
    peak = max(max(row) for row in grid.eta)
    level = round(0.5 * peak, 6)
    old = old_buttock_curves(grid.eta, grid.xi, grid.zeta, (level,))
    new = buttock_curves(grid.eta, grid.xi, grid.zeta, levels=(level,))
    noise = perturb_table(hull, amplitude=0.004 * 0.5 * hull.beam, seed=7)
    noisy_grid, _ = _grid_and_flat(noise)
    levels = buttock_coverage(noisy_grid.eta, noisy_grid.xi, noisy_grid.zeta)["levels"]
    noisy_old = old_buttock_curves(noisy_grid.eta, noisy_grid.xi, noisy_grid.zeta, levels)
    noisy_new = buttock_curves(noisy_grid.eta, noisy_grid.xi, noisy_grid.zeta, levels=levels)
    result = fair_offset_table(
        noise, tolerance=0.032 * 0.5 * hull.beam, iterations=3, reference=hull
    )
    fair_grid, _ = _grid_and_flat(result.table_after)
    faired_new = buttock_curves(fair_grid.eta, fair_grid.xi, fair_grid.zeta, levels=levels)
    # exactness of the inversion against the continuous analytic shape
    shape = [max(0.0, (z / hull.design_draft / 0.65) ** 0.7)
             if z / hull.design_draft <= 0.65
             else 1.0 - 0.55 * (z / hull.design_draft - 0.65) / 0.35 for z in grid.zeta]
    analytic = []
    for fraction in (0.25, 0.5, 0.75):
        shape_value = fraction
        below = [z for z, s in zip(grid.zeta, shape) if s <= shape_value]
        above = [z for z, s in zip(grid.zeta, shape) if s >= shape_value]
        analytic.append((shape_value, min(above) if above else None))
    traced_level = new[0][0][1]
    evidence["tumblehome"] = {
        "level": level,
        "old_curves": len(old),
        "old_samples": [len(curve) for curve in old],
        "new_curves": len(new),
        "new_zetas": sorted(round(curve[0][1], 6) for curve in new),
        "deck_ordinate": max(row[-1] for row in grid.eta),
        "old_flips_noisy": sum(polyline_sign_flips(curve) for curve in noisy_old),
        "new_flips_noisy": sum(polyline_sign_flips(curve) for curve in noisy_new),
        "new_flips_faired": sum(polyline_sign_flips(curve) for curve in faired_new),
        "fairing_deviation_m": result.deviation_max,
        "fairing_tolerance_m": result.tolerance,
        "fairing_accepted": result.accepted(),
        "total_flips": [result.flips_before, result.flips_after],
        "buttock_flips": [
            [metric.flips_before, metric.flips_after]
            for metric in result.iterations[0].metrics
            if metric.direction == "buttocks"
        ][0],
        "trace_probe_zeta": traced_level,
        "grid_waterline_of_level": analytic,
    }

    # --- figure 1: the branch the old contour dropped
    figure, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(11.0, 4.2))
    ax_a.plot(grid.zeta, grid.eta[len(grid.xi) // 2], "o-", ms=3, lw=1.2,
              label="section amidships: y(z)")
    ax_a.axhline(level, color="tab:red", lw=1.2, ls="--",
                 label=f"buttock level y = {level:.3f}")
    for z in evidence["tumblehome"]["new_zetas"]:
        ax_a.axvline(z, color="tab:green", lw=1.0, ls=":")
    ax_a.set_xlabel("z / T")
    ax_a.set_ylabel("y / (B/2)")
    ax_a.set_title("tumblehome section: one level, two crossings")
    ax_a.legend(fontsize=8)
    ax_a.grid(alpha=0.3)

    for curve in old:
        ax_b.plot([point[0] for point in curve], [point[1] for point in curve], "s",
                  ms=4, color="tab:gray", label="old contour (1 branch)")
    for index, curve in enumerate(new):
        ax_b.plot([point[0] for point in curve], [point[1] for point in curve], "-", lw=1.6,
                  label="traced buttock" if index == 0 else None)
    ax_b.set_xlabel("x / L")
    ax_b.set_ylabel("z / T")
    ax_b.set_title(f"buttock y = {level:.3f}: {len(old)} curve (old) vs {len(new)} (fixed)")
    ax_b.legend(fontsize=8)
    ax_b.grid(alpha=0.3)
    figure.tight_layout()
    figure.savefig(OUT / "buttock_branches.png", dpi=150)
    plt.close(figure)

    # --- figure 2: third-family coverage and flips, old vs new
    figure, (ax_c, ax_d) = plt.subplots(1, 2, figsize=(11.0, 4.2))
    labels = list(FORMS)
    x = range(len(labels))
    old_samples = [min(evidence["forms"][key]["old_samples"]) for key in labels]
    new_samples = [min(evidence["forms"][key]["new_samples"]) for key in labels]
    ax_c.bar([value - 0.2 for value in x], old_samples, width=0.4, label="old levels (0.2..0.8 peak)")
    ax_c.bar([value + 0.2 for value in x], new_samples, width=0.4, label="levels from table")
    ax_c.set_xticks(list(x))
    ax_c.set_xticklabels(labels, fontsize=8, rotation=10)
    ax_c.set_ylabel("samples on the shortest buttock")
    ax_c.set_title("third-family coverage (21 x 17 grid)")
    ax_c.legend(fontsize=8)
    ax_c.grid(alpha=0.3, axis="y")

    old_flips = [evidence["forms"][key]["old_flips_noisy"] for key in labels]
    new_flips = [evidence["forms"][key]["new_flips_noisy"] for key in labels]
    ax_d.bar([value - 0.2 for value in x], old_flips, width=0.4, label="old contour")
    ax_d.bar([value + 0.2 for value in x], new_flips, width=0.4, label="branch tracing")
    ax_d.set_xticks(list(x))
    ax_d.set_xticklabels(labels, fontsize=8, rotation=10)
    ax_d.set_ylabel("buttock curvature sign flips (noisy table)")
    ax_d.set_title("flips the metric could see before the fairing")
    ax_d.legend(fontsize=8)
    ax_d.grid(alpha=0.3, axis="y")
    figure.tight_layout()
    figure.savefig(OUT / "buttock_coverage.png", dpi=150)
    plt.close(figure)

    (OUT / "buttock_branch_evidence.json").write_text(
        json.dumps(evidence, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps({
        "forms": {key: {k: v for k, v in value.items() if k != "coverage"}
                  for key, value in evidence["forms"].items()},
        "tumblehome": evidence["tumblehome"],
    }, indent=1)[:4000])


if __name__ == "__main__":
    main()