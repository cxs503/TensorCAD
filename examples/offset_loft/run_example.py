#!/usr/bin/env python3
"""Offset table -> NURBS loft: the GAP-CAD-T1-OFFSET-LOFT acceptance study.

Every number in ``results.json`` is computed here -- nothing is copied from a report,
and every ordinate comes from a documented closed form (:mod:`tensorcad.offset_forms`),
so the tables can be re-derived from the cited definition and their ``sha256`` checked.

The study answers three questions:

1. **Does the loft agree with the table's own Simpson integrals?**  Three published
   hull forms, table baseline vs NURBS quadrature vs the tessellated loft mesh, judged
   against the ``0.5 %`` acceptance and against the analytic reference integrals of the
   hull itself.
2. **How sensitive is that to the number of stations?**  The same hull is tabulated at
   11 ... 81 stations; the error of the loft is tracked against the finest table and
   against the closed form.  This is what justifies each landed table's resolution
   instead of asserting it.
3. **Does the waterline grid matter?**  A power-law section ``y ~ zeta**(1/m)`` carries a
   large share of its area in the first cell above the keel, so the waterline clustering
   ``z_j = T (j/(n-1))**p`` is swept for the two power-law hulls.

Outputs (``--output``, default ``docs/offset_loft``): ``results.json``, ``report.md``,
``plot.png`` and the landed tables ``tables/<name>.csv`` + ``.provenance.json``.

Usage::

    python examples/offset_loft/run_example.py [--output docs/offset_loft]
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tensorcad.offset_forms import builtin_forms, offset_table, reference_for  # noqa: E402
from tensorcad.offset_loft import (  # noqa: E402
    compare_offset_table,
    loft_offset_table,
    read_offset_bundle,
    station_sensitivity,
    table_integrals,
    write_offset_bundle,
)

#: The landed tables: hull -> (stations, waterlines, clustering exponent).
#:
#: The resolutions are the *output* of the station/waterline ladders below, not a guess:
#: the two power-law sections need clustered waterlines to resolve the first cell above
#: the keel, and the full-scale tanker needs 41 stations to resolve the ends of a
#: ``Cb = 0.81`` plan form (11 ... 31 stations still disagree by 0.3 ... 0.5 %).
LANDED = {
    "wigley_iii": dict(n_stations=21, n_waterlines=17, waterline_exponent=1.0),
    "series60_cb060": dict(n_stations=21, n_waterlines=17, waterline_exponent=3.0),
    "kvlcc2_class": dict(n_stations=41, n_waterlines=17, waterline_exponent=3.0),
}

STATION_COUNTS = (11, 15, 21, 31, 41, 61, 81)
CLUSTER_EXPONENTS = (1.0, 2.0, 3.0, 4.0)
TOLERANCE = 5e-3
MESH = dict(n_u=129, n_v=49)


def _reference_subset(reference: dict) -> dict:
    return {
        "volume": reference["volume"],
        "waterplane_area": reference["waterplane_area"],
        "lcb": reference["lcb"],
        "vcb": reference["vcb"],
    }


def run_study(output: Path) -> dict:
    forms = builtin_forms()
    results: dict = {
        "schema": "tensorcad.offset_loft.study/1",
        "tolerance": TOLERANCE,
        "landed_tables": {},
        "comparisons": {},
        "station_sensitivity": {},
        "waterline_sensitivity": {},
        "fast_path_cross_check": {},
        "failures": [],
    }

    # -- 1. land the tables, and prove the bundle round-trips -----------------
    table_dir = output / "tables"
    for key, config in LANDED.items():
        table = offset_table(key, **config)
        csv_path, json_path = write_offset_bundle(table, table_dir)
        reloaded = read_offset_bundle(csv_path)
        record = json.loads(json_path.read_text(encoding="utf-8"))
        results["landed_tables"][key] = {
            "csv": str(csv_path.relative_to(REPO_ROOT)),
            "provenance_json": str(json_path.relative_to(REPO_ROOT)),
            "digest": table.digest(),
            "csv_sha256": record["csv_sha256"],
            "n_stations": table.n_stations,
            "n_waterlines": table.n_waterlines,
            "waterline_exponent": table.waterline_exponent,
            "length": table.length,
            "beam": table.beam,
            "design_draft": table.design_draft,
            "official": table.provenance.official,
            "source": table.provenance.source,
            "licence": table.provenance.licence,
            "round_trip_identical": reloaded.digest() == table.digest(),
        }
        if not results["landed_tables"][key]["round_trip_identical"]:
            results["failures"].append(f"{key}: bundle round trip changed the ordinates")

    # -- 2. the acceptance ---------------------------------------------------
    for key, config in LANDED.items():
        table = offset_table(key, **config)
        reference = reference_for(key)
        comparison = compare_offset_table(
            table, reference=_reference_subset(reference), tolerance=TOLERANCE, **MESH
        )
        payload = comparison.to_dict()
        payload["reference_full"] = reference
        results["comparisons"][key] = payload
        if not comparison.passed:
            results["failures"].append(
                f"{key}: volume/waterplane deviation above {TOLERANCE:.1%}"
            )

    # -- 3. station-count sensitivity ---------------------------------------
    for key, config in LANDED.items():
        reference = reference_for(key)
        waterlines = config["n_waterlines"]
        exponent = config["waterline_exponent"]

        def build(n_stations: int, *, _key: str = key, _w: int = waterlines, _p: float = exponent):
            return offset_table(
                _key, n_stations=n_stations, n_waterlines=_w, waterline_exponent=_p
            )

        finest = table_integrals(build(max(STATION_COUNTS)))
        rows = station_sensitivity(
            build,
            STATION_COUNTS,
            reference=_reference_subset(reference),
            denominator=finest,
        )
        results["station_sensitivity"][key] = list(rows)
        passing = [row for row in rows if row["relative_error_volume_vs_denominator"] <= TOLERANCE
                   and row["relative_error_waterplane_vs_denominator"] <= TOLERANCE]
        if not passing:
            results["failures"].append(f"{key}: no station count reached the tolerance")

    # -- 4. waterline clustering sensitivity --------------------------------
    for key in ("series60_cb060", "kvlcc2_class"):
        reference = reference_for(key)
        rows = []
        for exponent in CLUSTER_EXPONENTS:
            table = offset_table(
                key,
                n_stations=LANDED[key]["n_stations"],
                n_waterlines=17,
                waterline_exponent=exponent,
            )
            baseline = table_integrals(table)
            loft = loft_offset_table(table)
            volume = loft.volume()
            waterplane = loft.waterplane_area()
            rows.append(
                {
                    "waterline_exponent": exponent,
                    "first_waterline": table.waterlines[1],
                    "volume_table_simpson": baseline.volume,
                    "volume_nurbs": volume,
                    "waterplane_table_simpson": baseline.waterplane_area,
                    "waterplane_nurbs": waterplane,
                    "relative_error_volume_vs_table": abs(volume - baseline.volume) / baseline.volume,
                    "relative_error_waterplane_vs_table": abs(
                        waterplane - baseline.waterplane_area
                    )
                    / baseline.waterplane_area,
                    "relative_error_volume_vs_reference": abs(volume - reference["volume"])
                    / reference["volume"],
                    "relative_error_waterplane_vs_reference": abs(
                        waterplane - reference["waterplane_area"]
                    )
                    / reference["waterplane_area"],
                }
            )
        results["waterline_sensitivity"][key] = rows

    # -- 5. fast basis-matrix path vs the slow derivative path ---------------
    coarse = offset_table("wigley_iii", n_stations=13, n_waterlines=9)
    loft = loft_offset_table(coarse)
    fast_volume = loft.volume()
    slow_volume = loft.volume_direct()
    fast_waterplane = loft.waterplane_area()
    slow_waterplane = loft.waterplane_area_direct()
    results["fast_path_cross_check"] = {
        "table": "wigley_iii 13 x 9",
        "volume_fast": fast_volume,
        "volume_slow": slow_volume,
        "volume_relative_difference": abs(fast_volume - slow_volume) / slow_volume,
        "waterplane_fast": fast_waterplane,
        "waterplane_slow": slow_waterplane,
        "waterplane_relative_difference": abs(fast_waterplane - slow_waterplane) / slow_waterplane,
    }
    if results["fast_path_cross_check"]["volume_relative_difference"] > 1e-9:
        results["failures"].append("fast basis-matrix volume disagrees with the derivative path")
    if results["fast_path_cross_check"]["waterplane_relative_difference"] > 1e-9:
        results["failures"].append(
            "fast basis-matrix waterplane disagrees with the derivative path"
        )
    return results


def make_plot(results: dict, output: Path) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(2, 2, figsize=(12.5, 9.0))

    # (a) body plan of the landed tables -----------------------------------
    axis = axes[0][0]
    for offset, (key, config) in enumerate(LANDED.items()):
        table = offset_table(key, **config)
        scale = max(row_i for row in table.half_breadths for row_i in row)
        stations = table.stations[:: max(1, table.n_stations // 9)]
        for x in stations:
            index = table.stations.index(x)
            ys = [y / scale for y in table.half_breadths[index]]
            zs = [z / table.design_draft for z in table.waterlines[: len(ys)]]
            axis.plot(ys, zs, "-", linewidth=0.8)
    axis.set_xlabel("y / (B/2)")
    axis.set_ylabel("z / T")
    axis.set_title("(a) landed tables: body plan (9 stations each, normalised)")
    axis.grid(True, alpha=0.3)

    # (b) station sensitivity ----------------------------------------------
    axis = axes[0][1]
    for key, rows in results["station_sensitivity"].items():
        counts = [row["n_stations"] for row in rows]
        errors = [row["relative_error_volume_vs_denominator"] for row in rows]
        axis.loglog(counts, errors, "o-", label=f"{key} (volume)")
    axis.axhline(results["tolerance"], color="k", linestyle="--", linewidth=1, label="0.5 % target")
    axis.set_xlabel("stations in the table")
    axis.set_ylabel("relative volume error vs 81-station table")
    axis.set_title("(b) station-count sensitivity")
    axis.grid(True, which="both", alpha=0.3)
    axis.legend(fontsize=7)

    # (c) waterline clustering ---------------------------------------------
    axis = axes[1][0]
    for key, rows in results["waterline_sensitivity"].items():
        exponents = [row["waterline_exponent"] for row in rows]
        axis.semilogy(
            exponents,
            [row["relative_error_volume_vs_reference"] for row in rows],
            "o-",
            label=f"{key} (volume vs closed form)",
        )
    axis.axhline(results["tolerance"], color="k", linestyle="--", linewidth=1)
    axis.set_xlabel("waterline clustering exponent p  (z_j = T (j/(n-1))^p)")
    axis.set_ylabel("relative error vs analytic reference")
    axis.set_title("(c) waterline clustering, 17 waterlines")
    axis.grid(True, which="both", alpha=0.3)
    axis.legend(fontsize=7)

    # (d) the acceptance ---------------------------------------------------
    axis = axes[1][1]
    keys = list(results["comparisons"])
    labels = ["volume\n(S-NURBS)", "volume\n(mesh)", "waterplane\n(S-NURBS)", "waterplane\n(mesh)"]
    positions = range(len(keys))
    width = 0.2
    for series, label in enumerate(labels):
        values = []
        for key in keys:
            entry = results["comparisons"][key]
            block = "volume" if series < 2 else "waterplane"
            field = "relative_error_nurbs" if series % 2 == 0 else "relative_error_mesh"
            values.append(entry[block][field] * 100.0)
        axis.bar([p + (series - 1.5) * width for p in positions], values, width, label=label)
    axis.axhline(results["tolerance"] * 100.0, color="k", linestyle="--", linewidth=1)
    axis.set_xticks(list(positions))
    axis.set_xticklabels([key.replace("_", "\n") for key in keys], fontsize=8)
    axis.set_ylabel("relative deviation vs table Simpson  [%]")
    axis.set_title("(d) acceptance: loft vs table Simpson (dashed = 0.5 %)")
    axis.grid(True, axis="y", alpha=0.3)
    axis.legend(fontsize=7)

    figure.tight_layout()
    path = output / "plot.png"
    figure.savefig(path, dpi=130)
    plt.close(figure)
    return path


def write_report(results: dict, output: Path) -> Path:
    comparisons = results["comparisons"]
    lines = [
        "# Offset table -> NURBS loft: acceptance study",
        "",
        "Generated by `examples/offset_loft/run_example.py`; every number below is read",
        "back from `results.json`, which the same run computed.",
        "",
        f"Acceptance: loft displacement **and** waterplane area within "
        f"**{results['tolerance']:.1%}** of the table's own Simpson integral.",
        "",
        "## Landed tables",
        "",
        "| table | stations | waterlines | p | L [m] | B [m] | T [m] | digest |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for key, entry in results["landed_tables"].items():
        lines.append(
            f"| {key} | {entry['n_stations']} | {entry['n_waterlines']} | "
            f"{entry['waterline_exponent']:.0f} | {entry['length']:.4g} | "
            f"{entry['beam']:.4g} | {entry['design_draft']:.4g} | "
            f"`{entry['digest'][:16]}...` |"
        )
    lines += [
        "",
        "Each table is a CSV plus a provenance sidecar (source citation, method, URL,",
        "licence, `sha256` of the ordinates). Re-reading the bundle verifies both digests:",
        f"round trip identical for every table = "
        f"{all(e['round_trip_identical'] for e in results['landed_tables'].values())}.",
        "",
        "## Acceptance",
        "",
        "| hull | dV vs table | dV mesh | dV vs closed form | dA vs table | dA mesh | dA vs closed form | pass |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for key, entry in comparisons.items():
        volume, waterplane = entry["volume"], entry["waterplane"]
        lines.append(
            f"| {key} | {volume['relative_error_nurbs']:.4%} | {volume['relative_error_mesh']:.4%} | "
            f"{volume['relative_error_nurbs_vs_reference']:.4%} | "
            f"{waterplane['relative_error_nurbs']:.4%} | {waterplane['relative_error_mesh']:.4%} | "
            f"{waterplane['relative_error_nurbs_vs_reference']:.4%} | {entry['passed']} |"
        )
    lines += [
        "",
        "## Station-count sensitivity",
        "",
        "| hull | stations | dV vs 81-station table | dA vs 81-station table | dV vs closed form |",
        "|---|---|---|---|---|",
    ]
    for key, rows in results["station_sensitivity"].items():
        for row in rows:
            lines.append(
                f"| {key} | {row['n_stations']} | "
                f"{row['relative_error_volume_vs_denominator']:.4%} | "
                f"{row['relative_error_waterplane_vs_denominator']:.4%} | "
                f"{row.get('relative_error_volume_vs_reference', float('nan')):.4%} |"
            )
    lines += [
        "",
        "## Waterline clustering (17 waterlines)",
        "",
        "| hull | p | first waterline [m] | dV vs table | dV vs closed form | dA vs closed form |",
        "|---|---|---|---|---|---|",
    ]
    for key, rows in results["waterline_sensitivity"].items():
        for row in rows:
            lines.append(
                f"| {key} | {row['waterline_exponent']:.0f} | {row['first_waterline']:.4g} | "
                f"{row['relative_error_volume_vs_table']:.4%} | "
                f"{row['relative_error_volume_vs_reference']:.4%} | "
                f"{row['relative_error_waterplane_vs_reference']:.4%} |"
            )
    check = results["fast_path_cross_check"]
    lines += [
        "",
        "## Fast path cross-check",
        "",
        f"- {check['table']}: volume {check['volume_fast']:.12g} (basis matrices) vs "
        f"{check['volume_slow']:.12g} (surface derivatives), relative difference "
        f"{check['volume_relative_difference']:.2e}",
        f"- waterplane {check['waterplane_fast']:.12g} vs {check['waterplane_slow']:.12g}, "
        f"relative difference {check['waterplane_relative_difference']:.2e}",
        "",
        "## Verdict",
        "",
    ]
    if results["failures"]:
        lines.append("**FAIL**")
        lines += [f"- {failure}" for failure in results["failures"]]
    else:
        lines.append(
            "All three landed tables pass: displacement and waterplane area of the lofted "
            "patch match the table Simpson integrals within "
            f"{max(max(e['volume']['relative_error_nurbs'], e['waterplane']['relative_error_nurbs']) for e in comparisons.values()):.3%}"
            f" (NURBS quadrature) and "
            f"{max(max(e['volume']['relative_error_mesh'], e['waterplane']['relative_error_mesh']) for e in comparisons.values()):.3%}"
            " (loft mesh), against a 0.5 % requirement."
        )
    path = output / "report.md"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default=str(REPO_ROOT / "docs" / "offset_loft"))
    arguments = parser.parse_args(argv)
    output = Path(arguments.output)
    output.mkdir(parents=True, exist_ok=True)

    results = run_study(output)
    results["plot"] = str(make_plot(results, output).relative_to(REPO_ROOT))
    report = write_report(results, output)
    results["report"] = str(report.relative_to(REPO_ROOT))
    results_path = output / "results.json"
    results_path.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {results_path}, {report}, {results['plot']}")
    for failure in results["failures"]:
        print("FAIL:", failure)
    return 0 if not results["failures"] else 1


if __name__ == "__main__":
    raise SystemExit(main())