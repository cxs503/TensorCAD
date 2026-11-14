"""Evidence for the unified hull-family entry (``tensorcad.hull`` · ``family="series60"``).

Produces, from the public API only:

* ``artifacts/hull_family/hull_family.json`` -- five Series 60 parent forms through
  the unified entry, each against its own closed form and against the family
  module's own generator (bit-for-bit), the grid ladder, the transom and bulb
  pass-throughs, and the refusal matrix.
* ``docs/hull_family.png`` -- midship sections from the generated meshes, the
  volume-error ladder, and the breadth profile with and without a transom.

Run: ``python tools/hull_family/report.py`` (standard library + matplotlib for the
figure). Every number in ``docs/hull_family.md`` is injected from the JSON, not
transcribed.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from tensorcad.hull import HullParameters, generate_hull  # noqa: E402
from tensorcad.series60 import (  # noqa: E402
    Series60Parameters,
    analytic_series60_hydrostatics,
    generate_series60,
)

PARENTS = (0.60, 0.65, 0.70, 0.75, 0.80)
GRIDS = ((41, 13), (81, 21), (161, 41))
MAIN = dict(length=10.0, beam=1.5, draft=0.6, water_density=1000.0)
OUT_JSON = ROOT / "artifacts" / "hull_family" / "hull_family.json"
OUT_PNG = ROOT / "docs" / "hull_family.png"


def unified(**overrides) -> HullParameters:
    base = dict(MAIN, family="series60", block_coefficient=0.70, nx=41, nz=13)
    base.update(overrides)
    return HullParameters(**base)


def direct(**overrides) -> Series60Parameters:
    base = dict(MAIN, cb=0.70, nx=41, nz=13)
    base.update(overrides)
    return Series60Parameters(**base)


def midship_section(model, x_expect: float, tolerance: float = 1e-9):
    """Half-breadth against height, read off the generated mesh at one station."""
    points = [
        (vertex[2], abs(vertex[1]))
        for vertex in model.mesh.vertices
        if abs(vertex[0] - x_expect) <= tolerance
    ]
    points.sort()
    return points


def breadth_profile(model):
    return [width for _x, width in model.offsets]


def main() -> dict:
    evidence: dict = {"case": "unified hull-family entry: family='series60'", "families": ["sine", "series60"],
                      "principal_dimensions_m": {"L": MAIN["length"], "B": MAIN["beam"], "T": MAIN["draft"]},
                      "water_density_kg_m3": MAIN["water_density"], "parents": {}}

    # ---- five parents x three grids, against the closed form and the family module
    for cb in PARENTS:
        closed = analytic_series60_hydrostatics(direct(cb=cb))
        rows = []
        for nx, nz in GRIDS:
            model = generate_hull(unified(block_coefficient=cb, nx=nx, nz=nz))
            twin = generate_series60(direct(cb=cb, nx=nx, nz=nz))
            hydro = model.hydrostatics
            volume = hydro.displacement_volume
            rows.append({
                "nx": nx, "nz": nz,
                "triangles": len(model.mesh.faces),
                "displacement_volume_m3": volume,
                "volume_error_vs_closed_form": volume / closed["displacement_volume"] - 1.0,
                "waterplane_area_m2": hydro.waterplane_area,
                "waterplane_error_vs_closed_form": hydro.waterplane_area / closed["waterplane_area"] - 1.0,
                "block_coefficient_measured": hydro.block_coefficient,
                "block_coefficient_target": cb,
                "block_coefficient_error": hydro.block_coefficient / cb - 1.0,
                "vcb_z_m": hydro.vertical_centre_buoyancy,
                "vcb_error_vs_closed_form": hydro.vertical_centre_buoyancy
                / closed["vertical_centre_buoyancy_z"] - 1.0,
                "identical_to_family_module": (
                    model.mesh.vertices == twin.mesh.vertices
                    and model.mesh.faces == twin.mesh.faces
                    and volume == twin.hydrostatics.displacement_volume
                ),
                "max_station_width_m": max(breadth_profile(model)),
            })
        evidence["parents"][f"{cb:.2f}"] = {
            "closed_form": {key: closed[key] for key in (
                "displacement_volume", "waterplane_area", "maximum_section_area",
                "vertical_centre_buoyancy_z", "block_coefficient", "prismatic_coefficient",
                "midship_coefficient", "waterplane_coefficient",
                "longitudinal_centre_buoyancy")},
            "grids": rows,
        }

    # ---- transom and bulb pass-through
    plain = generate_hull(unified(nx=81, nz=21))
    transom = generate_hull(unified(nx=81, nz=21, transom_ratio=0.15))
    bulbed = generate_hull(unified(nx=81, nz=21, bulb_length_ratio=0.03, bulb_area_ratio=0.08))
    bulb = bulbed.metadata()["family_report"]["appendages"]["bulb"]
    evidence["transom"] = {
        "ratio": 0.15,
        "aft_station_width_m": transom.offsets[0][1],
        "expected_width_m": 0.15 * MAIN["beam"],
        "plain_aft_station_width_m": plain.offsets[0][1],
        "watertight": bool(transom.validation.errors == []),
        "displacement_volume_m3": transom.hydrostatics.displacement_volume,
        "plain_displacement_volume_m3": plain.hydrostatics.displacement_volume,
        "immersed_area_m2": transom.metadata()["family_report"]["transom"]["immersed_area_m2"],
        "immersed_area_ratio_of_mid_section":
            transom.metadata()["family_report"]["transom"]["immersed_area_ratio_of_mid_section"],
    }
    evidence["bulb"] = {
        "length_ratio": bulb["length_ratio"], "area_ratio": bulb["area_ratio"],
        "volume_m3": bulb["volume_m3"], "fp_area_m2": bulb["fp_area_m2"],
        "volume_share_of_hull": bulb["volume_m3"] / plain.hydrostatics.displacement_volume,
        "hull_displacement_unchanged": (
            bulbed.hydrostatics.displacement_volume == plain.hydrostatics.displacement_volume
        ),
        "appendage_meshes": bulbed.metadata()["appendages"]["bodies"],
    }

    # ---- refusals (the parameters a family cannot honour)
    refusals = []

    def expect_refusal(label, parameters, needle):
        try:
            parameters.validate()
        except ValueError as exc:
            refusals.append({"case": label, "refused": True,
                             "message_has_needle": needle in str(exc), "message": str(exc)})
        else:
            refusals.append({"case": label, "refused": False, "message_has_needle": False,
                             "message": "accepted -- this is a hole"})

    expect_refusal("series60 without block_coefficient", unified(block_coefficient=None),
                   "block_coefficient")
    expect_refusal("sine family with a Cb", HullParameters(block_coefficient=0.70),
                   "only apply to family='series60'")
    expect_refusal("sine family with a transom", HullParameters(transom_ratio=0.2),
                   "only apply to family='series60'")
    expect_refusal("series60 with a sine exponent", unified(fore_exponent=1.5),
                   "sine-family parameters")
    expect_refusal("unknown family", HullParameters(family="wigley"), "family must be one of")
    expect_refusal("Cb outside the ladder", unified(block_coefficient=0.95), "cb must lie in")
    evidence["refusals"] = refusals

    # ---- breadth profiles for the figure
    evidence["profiles"] = {
        "plain": breadth_profile(plain),
        "transom": breadth_profile(transom),
        "x": [x for x, _w in plain.offsets],
        "mid_sections": {
            f"{cb:.2f}": midship_section(
                generate_hull(unified(block_coefficient=cb, nx=161, nz=41)), 0.5 * MAIN["length"]
            ) for cb in PARENTS
        },
    }
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(evidence, indent=2, sort_keys=True))
    return evidence


def figure(evidence: dict) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.8))

    # (a) midship sections read off the generated meshes
    ax = axes[0]
    colours = plt.cm.viridis([i / (len(PARENTS) - 1) for i in range(len(PARENTS))])
    for colour, cb in zip(colours, PARENTS):
        points = evidence["profiles"]["mid_sections"][f"{cb:.2f}"]
        z = [p[0] for p in points]
        y = [p[1] for p in points]
        ax.plot(y, z, color=colour, lw=1.6, label=f"Cb = {cb:.2f}")
    ax.set_xlabel("half-breadth y [m]")
    ax.set_ylabel("z [m]  (0 = design waterline)")
    ax.set_title("(a) midship section, 161 x 41 mesh")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="lower right")

    # (b) volume error against the closed form, per grid, per parent
    ax = axes[1]
    for colour, cb in zip(colours, PARENTS):
        rows = evidence["parents"][f"{cb:.2f}"]["grids"]
        errors = [100.0 * row["volume_error_vs_closed_form"] for row in rows]
        ax.plot([f"{r['nx']}x{r['nz']}" for r in rows], errors, "o-", color=colour,
                lw=1.4, label=f"Cb = {cb:.2f}")
    ax.axhline(-2.0, color="crimson", ls="--", lw=1.0, label="2 % warning band")
    ax.set_ylabel("volume error vs closed form [%]")
    ax.set_xlabel("grid")
    ax.set_title("(b) loft under-fills the keel turn, 1st order")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)

    # (c) breadth profile, with and without a transom
    ax = axes[2]
    x = [xi / MAIN["length"] for xi in evidence["profiles"]["x"]]
    ax.plot(x, evidence["profiles"]["plain"], color="#1f77b4", lw=1.8, label="pointed stern")
    ax.plot(x, evidence["profiles"]["transom"], color="#d62728", lw=1.8, ls="--",
            label="transom ratio 0.15")
    ax.set_xlabel("x / L  (0 = aft end)")
    ax.set_ylabel("station width [m]")
    ax.set_title("(c) transom raises the aft deck breadth")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)

    fig.suptitle("tensorcad.hull family entry: Series 60 through the unified parameter set", fontsize=11)
    fig.tight_layout()
    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=150)
    print(f"wrote {OUT_PNG.relative_to(ROOT)}")


if __name__ == "__main__":
    data = main()
    print(f"wrote {OUT_JSON.relative_to(ROOT)}")
    worst = max(
        abs(row["volume_error_vs_closed_form"])
        for parent in data["parents"].values() for row in parent["grids"]
    )
    print(f"worst |volume error| over all parents and grids: {100.0 * worst:.3f} %")
    print("identical-to-family-module:", all(
        row["identical_to_family_module"]
        for parent in data["parents"].values() for row in parent["grids"]
    ))
    figure(data)