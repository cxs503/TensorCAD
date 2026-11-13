"""Evidence for docs/ring_waterline_regression.md: every number is computed here."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import matplotlib.pyplot as plt

from tensorcad import generate_hull, generate_wigley

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("trim_tests", ROOT / "tests" / "test_trim.py")
trim_tests = importlib.util.module_from_spec(spec)
spec.loader.exec_module(trim_tests)

from tensorcad.trim import SEAWATER_DENSITY, TrimError, hydrostatic_particulars  # noqa: E402


def ring_levels(mesh):
    return sorted({round(vertex[2], 12) for vertex in mesh.vertices})


def ring_table():
    rows = []
    for name, factory in (("hull", generate_hull), ("wigley", generate_wigley)):
        mesh = factory().mesh
        levels = ring_levels(mesh)
        keel, deck = levels[0], levels[-1]
        for z in levels[1:-1]:
            row = hydrostatic_particulars(mesh, z, station_count=21)
            jumps = {}
            for dz in (-1e-7, 1e-7):
                neighbour = hydrostatic_particulars(mesh, z + dz, station_count=21)
                jumps[f"{dz:+.0e}"] = {
                    "area": abs(neighbour.waterplane_area - row.waterplane_area) / row.waterplane_area,
                    "tpc": abs(neighbour.tpc - row.tpc) / row.tpc,
                    "volume": abs(neighbour.volume - row.volume) / row.volume,
                    "lcb": abs(neighbour.lcb - row.lcb),
                    "nudged": neighbour.waterline_nudged,
                }
            rows.append(
                {
                    "geometry": name,
                    "z": z,
                    "draft": row.draft,
                    "waterplane_area": row.waterplane_area,
                    "tpc": row.tpc,
                    "volume": row.volume,
                    "lcb": row.lcb,
                    "nudged": row.waterline_nudged,
                    "jumps": jumps,
                }
            )
        # boundaries
        rows.append({"geometry": name, "keel_ring": keel, "deck_ring": deck})
    return rows


def boundary_evidence():
    hull = generate_hull()
    mesh = hull.mesh
    levels = ring_levels(mesh)
    keel, deck = levels[0], levels[-1]
    refused = False
    message = ""
    try:
        hydrostatic_particulars(mesh, keel, station_count=21)
    except TrimError as exc:
        refused = True
        message = str(exc)
    on_ring = hydrostatic_particulars(mesh, deck, station_count=21)
    above = hydrostatic_particulars(mesh, deck + 1e-7, station_count=21)
    return {
        "keel_ring_z": keel,
        "keel_ring_refused": refused,
        "keel_ring_message": message,
        "deck_ring_z": deck,
        "deck_ring_on_ring_nudged": on_ring.waterline_nudged,
        "deck_ring_above_nudged": above.waterline_nudged,
        "deck_ring_area_on_ring": on_ring.waterplane_area,
        "deck_ring_area_above": above.waterplane_area,
        "deck_ring_area_rel_difference": abs(above.waterplane_area - on_ring.waterplane_area)
        / on_ring.waterplane_area,
    }


def prism_evidence():
    levels = (-2.0, -1.5, -1.0, -0.5, 0.0)
    length, root, taper = 10.0, 2.0, 0.5
    mesh = trim_tests._tapered_prism(levels, length=length, root_half_breadth=root, taper=taper)
    keel = levels[0]

    def half_breadth(z):
        return root + taper * z

    rows = []
    for z in levels[1:]:
        width = 2.0 * half_breadth(z)
        area = length * width
        volume = length * (z - keel) * (half_breadth(keel) + half_breadth(z))
        mct = SEAWATER_DENSITY / 1000.0 * (area * length * length / 12.0) / (100.0 * length)
        tpc = SEAWATER_DENSITY / 1000.0 * area / 100.0
        row = hydrostatic_particulars(mesh, z, station_count=21)
        rows.append(
            {
                "z": z,
                "measured": {
                    "waterplane_area": row.waterplane_area,
                    "lcf": row.lcf,
                    "volume": row.volume,
                    "tpc": row.tpc,
                    "mct": row.mct,
                },
                "closed_form": {"waterplane_area": area, "lcf": 0.5 * length, "volume": volume,
                                "tpc": tpc, "mct": mct},
                "rel_error": {
                    "waterplane_area": abs(row.waterplane_area - area) / area,
                    "lcf": abs(row.lcf - 0.5 * length) / (0.5 * length),
                    "volume": abs(row.volume - volume) / volume,
                    "tpc": abs(row.tpc - tpc) / tpc,
                    "mct": abs(row.mct - mct) / mct,
                },
                "nudged": row.waterline_nudged,
            }
        )
    return {"mesh_triangles": len(mesh), "levels": list(levels), "rows": rows}


def profile():
    """Waterplane area through a ring level, on a geometric ladder of offsets."""
    hull = generate_hull()
    z0 = -0.25 * hull.parameters.draft
    offsets = [0.0] + [(-1.0) ** k * 10.0 ** (-e) for e in range(2, 9) for k in (0,)]
    offsets = [0.0]
    for e in range(2, 9):
        offsets += [-10.0 ** (-e), 10.0 ** (-e)]
    samples = []
    for dz in offsets:
        row = hydrostatic_particulars(hull.mesh, z0 + dz, station_count=41)
        samples.append({"dz": dz, "area": row.waterplane_area, "nudged": row.waterline_nudged})
    base = next(s["area"] for s in samples if s["dz"] == 0.0)
    for s in samples:
        s["rel_difference"] = abs(s["area"] - base) / base
    return {"ring_z": z0, "station_count": 41, "samples": samples}


def figure(tables, prism, prof, path):
    figure, (left, right) = plt.subplots(1, 2, figsize=(12.0, 4.6))
    hull_rows = [r for r in tables if r.get("geometry") == "hull" and "z" in r]
    zs = [r["z"] for r in hull_rows]
    areas = [r["waterplane_area"] for r in hull_rows]
    left.plot(zs, areas, "o-", label="hydrostatic_particulars at ring levels", color="#1f4e79")
    left.set_xlabel("waterline z [m]")
    left.set_ylabel("waterplane area [m²]")
    left.set_title("Ring levels are ordinary waterlines (generated hull)")
    left.grid(alpha=0.3)

    offsets = [s["dz"] if s["dz"] != 0.0 else 1e-18 for s in prof["samples"]]
    diffs = [s["rel_difference"] for s in prof["samples"]]
    right.loglog([abs(o) for o in offsets], [max(d, 1e-18) for d in diffs], "s-", color="#a33")
    right.set_xlabel("|z - z_ring| [m]")
    right.set_ylabel("relative change of waterplane area")
    right.set_title(f"Continuity through the ring at z = {prof['ring_z']:.3f} m")
    right.grid(alpha=0.3, which="both")
    figure.tight_layout()
    figure.savefig(path, dpi=150)
    plt.close(figure)


def acceptance_figure(worst_exact, path):
    figure, axis = plt.subplots(figsize=(7.0, 4.0))
    names = list(worst_exact)
    values = [max(worst_exact[name], 1e-18) for name in names]
    axis.bar(names, values, color="#1f4e79")
    axis.set_yscale("log")
    axis.set_ylim(1e-18, 1.0)
    axis.axhline(3e-2, color="#a33", linestyle="--", label="3 % acceptance line")
    axis.axhline(1e-12, color="#777", linestyle=":", label="1e-12 (test tolerance floor)")
    axis.set_ylabel("relative error vs closed form")
    axis.set_title("Tapered prism: ring waterplane quantities vs the analytic answer")
    axis.legend(loc="upper right", fontsize=8)
    axis.grid(alpha=0.3, axis="y")
    for index, value in enumerate(values):
        axis.annotate(f"{value:.1e}", (index, value), ha="center", va="bottom", fontsize=8)
    figure.tight_layout()
    figure.savefig(path, dpi=150)
    plt.close(figure)


def main():
    tables = ring_table()
    prism = prism_evidence()
    prof = profile()
    boundaries = boundary_evidence()
    docs = ROOT / "docs"
    docs.mkdir(exist_ok=True)
    figure(tables, prism, prof, docs / "ring_waterline_regression.png")

    worst_jump = {"area": 0.0, "tpc": 0.0, "volume": 0.0, "lcb": 0.0}
    for row in tables:
        if "jumps" not in row:
            continue
        for jump in row["jumps"].values():
            for key in worst_jump:
                worst_jump[key] = max(worst_jump[key], jump[key])
            assert jump["nudged"] is False, row
    worst_exact = {
        key: max(r["rel_error"][key] for r in prism["rows"]) for key in prism["rows"][0]["rel_error"]
    }
    acceptance_figure(worst_exact, docs / "ring_waterline_accuracy.png")

    payload = {
        "schema": "tensorcad.trim.ring_waterline_evidence/1",
        "rings": tables,
        "worst_neighbour_jump": worst_jump,
        "prism": prism,
        "prism_worst_relative_error": worst_exact,
        "profile_through_ring": prof,
        "boundaries": boundaries,
    }
    (docs / "ring_waterline_regression.json").write_text(json.dumps(payload, indent=2) + "\n")

    lines = [
        "# The vertex-ring waterline: every level, and a shape that knows the answer",
        "",
        "Generated by `work/ring_waterline_evidence.py`; numbers in",
        "`docs/ring_waterline_regression.json` (no value below is typed by hand).",
        "",
        "## 1. Why this exists",
        "",
        "The parametric hull generator puts vertex rings on evenly spaced levels, so a",
        "\"round\" waterline such as `-T/2` lands exactly on a ring. m16609 established that",
        "such a cut is ordinary (`waterline_nudged is False`) at a single level. The two",
        "tests added here generalise that: every interior ring of both shipped geometries,",
        "and a geometry whose waterplane has a closed form.",
        "",
        "## 2. Continuity at every interior ring",
        "",
        "Worst relative change over all rings of both geometries for a 1e-7 m step across",
        f"the ring: area {worst_jump['area']:.2e}, TPC {worst_jump['tpc']:.2e},",
        f"volume {worst_jump['volume']:.2e}, LCB {worst_jump['lcb']:.2e} m.",
        "No cut at a ring needed the fallback (`waterline_nudged` False everywhere).",
        "",
        "| geometry | z [m] | draft [m] | area [m²] | TPC [t/cm] | volume [m³] | max rel jump |",
        "|---|---|---|---|---|---|---|",
    ]
    for row in tables:
        if "jumps" not in row:
            continue
        jump = max(
            max(j["area"], j["tpc"], j["volume"]) for j in row["jumps"].values()
        )
        lines.append(
            f"| {row['geometry']} | {row['z']:+.6f} | {row['draft']:.6f} | "
            f"{row['waterplane_area']:.9f} | {row['tpc']:.9f} | {row['volume']:.9f} | {jump:.2e} |"
        )
    lines += [
        "",
        "## 3. Exactness against a closed form (tapered prism)",
        "",
        "`half_breadth(z) = 2.0 + 0.5 z` independently of `x`, so the section is the",
        "rectangle `length x 2 half_breadth(z)` at every height, ring included, and Simpson",
        "over stations is exact for the constant section area. Worst relative error over the",
        f"four rings below the deck: area {worst_exact['waterplane_area']:.2e},",
        f"LCF {worst_exact['lcf']:.2e}, volume {worst_exact['volume']:.2e},",
        f"TPC {worst_exact['tpc']:.2e}, MCT {worst_exact['mct']:.2e}.",
        "",
        "| z [m] | area [m²] | area err | volume [m³] | volume err | MCT err |",
        "|---|---|---|---|---|---|",
    ]
    for row in prism["rows"]:
        lines.append(
            f"| {row['z']:+.2f} | {row['measured']['waterplane_area']:.9f} | "
            f"{row['rel_error']['waterplane_area']:.2e} | {row['measured']['volume']:.9f} | "
            f"{row['rel_error']['volume']:.2e} | {row['rel_error']['mct']:.2e} |"
        )
    lines += [
        "",
        "## 4. The two boundaries, stated rather than implied",
        "",
        f"* keel ring `z = {boundaries['keel_ring_z']:+.6f} m` is refused"
        f" (`keel_ring_refused = {boundaries['keel_ring_refused']}`), message:"
        f" `{boundaries['keel_ring_message']}`;",
        f"* deck ring `z = {boundaries['deck_ring_z']:+.6f} m`: the cut exactly on the ring is"
        f" ordinary (`nudged = {boundaries['deck_ring_on_ring_nudged']}`), while 1e-7 m above"
        f" it - above the capped body - the fallback fires and declares"
        f" (`nudged = {boundaries['deck_ring_above_nudged']}`), moving the area by"
        f" {boundaries['deck_ring_area_rel_difference']:.2e} relative, inside the documented"
        " O(1e-6) budget.",
        "",
        "## 5. Figure",
        "",
        "`docs/ring_waterline_regression.png`: left, waterplane area at every interior ring",
        "level of the generated hull; right, relative change of the area as the waterline",
        "approaches the `-T/4` ring from both sides on a geometric ladder of offsets.",
        "",
    ]
    (docs / "ring_waterline_regression.md").write_text("\n".join(lines))
    print(json.dumps({"worst_jump": worst_jump, "worst_exact": worst_exact,
                      "boundaries": boundaries}, indent=2))


if __name__ == "__main__":
    main()
