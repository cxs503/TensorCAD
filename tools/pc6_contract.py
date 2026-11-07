#!/usr/bin/env python3
"""Build the PC6 geometry contract and the figures the report cites.

Everything in ``artifacts/geom_pc6.json`` and in the figures comes from this
script -- no number is typed by hand::

    python tools/pc6_contract.py [--output artifacts]

What it does:

1. generates the PC6 hull at the campaign design point (L 200 / B 24 / T 11),
   gates it (topology, dimensions, analytic displacement, self-intersection) and
   writes ``pc6_hull.stl``;
2. regenerates it for the beam sweep 24 -> 28 m and records the per-beam numbers
   plus the STLs under ``artifacts/pc6_sweep/``;
3. builds the de Long -> Pevek ice channel, writes ``pc6_channel.stl`` and checks
   that the corridor plan is simple;
4. runs the acceptance suite (``tests/test_pc6.py``) and records the counts;
5. writes ``artifacts/geom_pc6.json`` (mesh + implicit field + parameter table +
   acceptance block + digest) and three figures.

The exit status is non-zero when an acceptance criterion is false, so the script
cannot quietly publish a failing contract.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import re
import subprocess
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from dataclasses import replace

import tensorcad.pc6 as pc6
from tensorcad.icechannel import ice_channel
from tensorcad.pc6_contract import build_contract, load_contract

SCHEMA = "tensorcad.pc6-contract-run/1"
ROOT = pathlib.Path(__file__).resolve().parents[1]


def station_ring(p: pc6.PC6Parameters, s: float, levels: int = 41) -> list[tuple[float, float, float]]:
    """The closed section ring at ``s``: starboard up, deck across, port down."""
    ring: list[tuple[float, float, float]] = []
    for index in range(levels):
        u = index / (levels - 1)
        z = -p.draft * (1.0 - u)
        x = s * pc6.stem_x(p, z)
        ring.append((x, pc6.half_breadth(p, s, u), z))
    for index in range(levels - 1, -1, -1):
        u = index / (levels - 1)
        z = -p.draft * (1.0 - u)
        x = s * pc6.stem_x(p, z)
        ring.append((x, -pc6.half_breadth(p, s, u), z))
    return ring


def convergence_table(beams=None) -> list[dict[str, float]]:
    """Displacement and waterplane error against the analytic form, per grid."""
    p = pc6.PC6Parameters()
    form = pc6.analytic_form(p)
    rows: list[dict[str, float]] = []
    for stations, levels in ((41, 11), (57, 13), (73, 15), (97, 19), (129, 25), (161, 31)):
        model = pc6.generate(
            replace(p, longitudinal_stations=stations, vertical_levels=levels),
            check_self_intersection=False,
            volume_tolerance=0.05,
        )
        hydro = model.hydrostatics
        rows.append(
            {
                "stations": stations,
                "levels": levels,
                "faces": len(model.mesh.faces),
                "volume_relative_error": hydro.displacement_volume / form.volume - 1.0,
                "waterplane_relative_error": hydro.waterplane_area / form.waterplane_area - 1.0,
                "lcb_relative_error": hydro.longitudinal_centre_buoyancy / form.lcb - 1.0,
            }
        )
    return rows


def section_area_curve(p: pc6.PC6Parameters, count: int = 61) -> list[tuple[float, float]]:
    """Section area against station, from the analytic half-breadth profile."""
    curve: list[tuple[float, float]] = []
    levels = 200
    for index in range(count):
        s = index / (count - 1)
        area = 0.0
        previous = pc6.half_breadth(p, s, 0.0)
        for level in range(1, levels + 1):
            u = level / levels
            breadth = pc6.half_breadth(p, s, u)
            area += 0.5 * (breadth + previous) * p.draft / levels
            previous = breadth
        curve.append((s * p.length, 2.0 * area))
    return curve


def run_pytest() -> dict[str, object]:
    """Run the acceptance suite and parse its summary line."""
    completed = subprocess.run(
        [sys.executable, "-m", "pytest", "tests/test_pc6.py", "-q"],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    tail = completed.stdout.strip().splitlines()[-1] if completed.stdout.strip() else ""
    passed = int(match.group(1)) if (match := re.search(r"(\d+) passed", tail)) else 0
    failed = int(match.group(1)) if (match := re.search(r"(\d+) failed", tail)) else 0
    errors = int(match.group(1)) if (match := re.search(r"(\d+) error", tail)) else 0
    return {
        "command": "python -m pytest tests/test_pc6.py -q",
        "passed": passed,
        "failed": failed + errors,
        "summary": tail,
        "returncode": completed.returncode,
    }


def figures(contract: dict, output: pathlib.Path) -> dict[str, pathlib.Path]:
    """Three figures: the form, the validation numbers, the ice channel."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    p = pc6.PC6Parameters()
    form = pc6.analytic_form(p)
    paths: dict[str, pathlib.Path] = {}

    fig, axes = plt.subplots(2, 2, figsize=(9.0, 6.6))
    # profile: stem/keel outline and a few waterlines
    axes[0, 0].set_title("profile (z-x), x measured from the AP")
    for u in (0.0, 0.25, 0.5, 0.75, 1.0):
        z = -p.draft * (1.0 - u)
        xs = [s / 20.0 * pc6.stem_x(p, z) for s in range(21)]
        axes[0, 0].plot(xs, [z] * len(xs), lw=0.8)
    stem = [(pc6.stem_x(p, -p.draft * (1.0 - u)), -p.draft * (1.0 - u)) for u in (0.0, 1.0)]
    axes[0, 0].plot([v[0] for v in stem], [v[1] for v in stem], "k-", lw=1.4, label="raked stem")
    axes[0, 0].plot([0, 0], [-p.draft, 0], "k-", lw=1.4, label="transom")
    axes[0, 0].set_xlabel("x [m]")
    axes[0, 0].set_ylabel("z [m]")
    axes[0, 0].legend(fontsize=7, loc="lower right")

    # body plan: section rings at four stations
    axes[0, 1].set_title("body plan (y-z), four stations")
    for s in (0.0, 0.25, p.midship_fraction, 0.85, 1.0):
        ring = station_ring(p, s)
        axes[0, 1].plot([point[1] for point in ring], [point[2] for point in ring], lw=0.9,
                        label=f"s={s:.2f}")
    axes[0, 1].set_xlabel("y [m]")
    axes[0, 1].set_ylabel("z [m]")
    axes[0, 1].legend(fontsize=7)

    # waterline plan
    axes[1, 0].set_title("waterlines (x-y)")
    for u in (0.0, 0.25, 0.5, 0.75, 1.0):
        z = -p.draft * (1.0 - u)
        lengths = [s / 40.0 * pc6.stem_x(p, z) for s in range(41)]
        breadths = [pc6.half_breadth(p, s / 40.0, u) for s in range(41)]
        axes[1, 0].plot(lengths, breadths, lw=0.9)
        axes[1, 0].plot(lengths, [-b for b in breadths], lw=0.9)
    axes[1, 0].set_xlabel("x [m]")
    axes[1, 0].set_ylabel("y [m]")

    # section area curve against the block coefficient rectangle
    curve = section_area_curve(p)
    axes[1, 1].set_title("section area curve (mean of area x L / (L B T) = Cb)")
    areas = [area for _, area in curve]
    axes[1, 1].plot([x for x, _ in curve], areas, lw=1.2, label="2 * section area")
    axes[1, 1].axhline(p.beam * p.draft, color="0.6", ls=":", lw=0.8, label="B*T")
    mean = sum(areas) / len(areas)
    axes[1, 1].axhline(mean, color="tab:red", ls="--", lw=0.8,
                       label=f"mean = {mean:.1f} m^2 (Cb = {form.block_coefficient:.3f})")
    axes[1, 1].set_xlabel("x [m]")
    axes[1, 1].set_ylabel("area [m^2]")
    axes[1, 1].legend(fontsize=7)
    fig.tight_layout()
    path = output / "pc6_hull_form.png"
    fig.savefig(path, dpi=105)
    plt.close(fig)
    paths["form"] = path

    convergence = contract["_convergence"]
    sweep = contract["parameter_sweep"]
    fig, axes = plt.subplots(1, 3, figsize=(10.0, 3.3))
    axes[0].set_title("displacement error vs grid")
    axes[0].plot([row["faces"] for row in convergence],
                 [-row["volume_relative_error"] * 100.0 for row in convergence],
                 "o-", lw=1.0, ms=3, label="volume")
    axes[0].plot([row["faces"] for row in convergence],
                 [-row["waterplane_relative_error"] * 100.0 for row in convergence],
                 "s--", lw=1.0, ms=3, label="waterplane")
    axes[0].axhline(contract["acceptance"]["volume_tolerance"] * 100.0, color="tab:red", ls=":", lw=0.8,
                    label="gate tolerance")
    axes[0].set_xscale("log")
    axes[0].set_xlabel("triangles")
    axes[0].set_ylabel("|error| [%]")
    axes[0].legend(fontsize=7)

    axes[1].set_title("beam sweep: mesh vs analytic")
    beams = [row["beam_m"] for row in sweep["rows"]]
    axes[1].plot(beams, [row["displacement_volume_m3"] for row in sweep["rows"]], "o", ms=4,
                 label="mesh volume")
    axes[1].plot(beams, [row["analytic_volume_m3"] for row in sweep["rows"]], "x", ms=5,
                 label="analytic volume")
    axes[1].set_xlabel("beam B [m]")
    axes[1].set_ylabel("displacement [m^3]")
    axes[1].legend(fontsize=7)

    axes[2].set_title("implicit field at midship (y-z)")
    levels = 41
    ys = [-(0.55 * p.beam) + (1.1 * p.beam) * i / (levels - 1) for i in range(levels)]
    zs = [-p.draft * 1.02 + (1.04 * p.draft) * j / (levels - 1) for j in range(levels)]
    x_mid = p.midship_fraction * p.length
    grid = [[pc6.field_value(p, x_mid, y, z) for y in ys] for z in zs]
    contour = axes[2].contourf(ys, zs, grid, levels=[-1e9, 0.0, 1e9],
                               colors=["tab:blue", "white"])
    axes[2].plot([pc6.half_breadth(p, p.midship_fraction, (z + p.draft) / p.draft) for z in zs], zs,
                 "k-", lw=0.8)
    axes[2].plot([-pc6.half_breadth(p, p.midship_fraction, (z + p.draft) / p.draft) for z in zs], zs,
                 "k-", lw=0.8)
    axes[2].set_xlabel("y [m]")
    axes[2].set_ylabel("z [m]")
    axes[2].text(0.02, 0.05, "blue = f<0 (inside)", transform=axes[2].transAxes, fontsize=7)
    del contour
    fig.tight_layout()
    path = output / "pc6_hull_validation.png"
    fig.savefig(path, dpi=105)
    plt.close(fig)
    paths["validation"] = path

    channel = ice_channel()
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    ax.plot(list(channel.longitude), list(channel.latitude), "o-", ms=3, lw=1.1,
            label="route (great-circle segments)")
    # corridor walls, converted back to degrees for the map view only
    scale = 111.32  # km per degree of latitude; the local plane is in metres
    origin_lat, origin_lon = channel.origin
    for wall, style in ((channel.left, "--"), (channel.right, "--")):
        lons = [origin_lon + (x / 1000.0) / (scale * __import__("math").cos(__import__("math").radians(origin_lat)))
                for x, _ in wall]
        lats = [origin_lat + (y / 1000.0) / scale for _, y in wall]
        ax.plot(lons, lats, style, lw=0.8, color="tab:orange", alpha=0.8)
    ax.set_title(
        f"de Long -> Pevek corridor: {channel.total_km:.0f} km "
        f"({channel.total_km / 1.852:.0f} nm), width {channel.width:.0f} m"
    )
    ax.set_xlabel("longitude [deg E]")
    ax.set_ylabel("latitude [deg N]")
    for name, lon, lat in zip(channel.names, channel.longitude, channel.latitude, strict=True):
        ax.annotate(name, (lon, lat), fontsize=6, xytext=(3, 3), textcoords="offset points")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7)
    fig.tight_layout()
    path = output / "pc6_ice_channel.png"
    fig.savefig(path, dpi=105)
    plt.close(fig)
    paths["channel"] = path
    return paths


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="artifacts", help="artifact directory")
    parser.add_argument("--skip-tests", action="store_true", help="do not run pytest")
    args = parser.parse_args()

    output = (ROOT / args.output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    tests = run_pytest() if not args.skip_tests else {"passed": 0, "failed": 0, "summary": "skipped"}
    contract = build_contract(output, test_summary=tests)
    contract["_convergence"] = convergence_table()
    figure_paths = figures(contract, output)
    convergence = contract.pop("_convergence")
    contract["convergence"] = convergence

    # re-write the contract so the document on disk carries the convergence table
    from tensorcad.export import write_metadata

    write_metadata(json.loads(json.dumps(contract, default=float)), output / "geom_pc6.json")
    geometry = load_contract(output / "geom_pc6.json")

    summary = {
        "schema": SCHEMA,
        "artifact": str(output / "geom_pc6.json"),
        "figures": {name: str(path) for name, path in figure_paths.items()},
        "acceptance": contract["acceptance"],
        "hull": {
            "faces": contract["geometry"]["hull"]["topology"]["faces"],
            "vertices": contract["geometry"]["hull"]["topology"]["vertices"],
            "watertight": contract["geometry"]["hull"]["topology"]["watertight"],
            "euler_characteristic": contract["geometry"]["hull"]["topology"]["euler_characteristic"],
            "displacement_m3": contract["geometry"]["hydrostatics"]["mesh"]["displacement_volume_m3"],
            "mass_t": contract["geometry"]["hydrostatics"]["mass_displacement_t"],
            "block_coefficient": contract["geometry"]["hydrostatics"]["mesh"]["block_coefficient"],
            "relative_errors": contract["geometry"]["hydrostatics"]["relative_errors"],
        },
        "channel": {
            "total_km": contract["ice_channel"]["total_km"],
            "width_m": contract["ice_channel"]["corridor_width_m"],
            "watertight": contract["geometry"]["channel"]["topology"]["watertight"],
            "plan_is_simple": contract["geometry"]["channel"]["topology"]["plan_is_simple"],
        },
        "sweep": contract["parameter_sweep"],
        "convergence": convergence,
        "negative_controls": contract["negative_controls"],
        "tests": tests,
        "digest": contract["digest"]["value"],
        "reloaded_volume_m3": geometry.hydrostatics["displacement_volume_m3"],
    }
    (output / "pc6_contract_run.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n",
                                                 encoding="utf-8")
    print(json.dumps(summary["acceptance"], indent=2))
    print(f"hull faces={summary['hull']['faces']} displacement={summary['hull']['displacement_m3']:.1f} m^3")
    ok = (
        contract["acceptance"]["watertight_closed_manifold"]
        and contract["acceptance"]["beam_sweep_24_28_regenerates"]
        and contract["acceptance"]["implicit_field_agrees_with_mesh"]
        and contract["acceptance"]["negative_cases_detected"]
        and summary["tests"]["failed"] == 0
        and summary["tests"]["passed"] >= 6
    )
    print("ACCEPTANCE:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
