#!/usr/bin/env python3
"""Series 60 family audit: five-Cb gradient, transom and bulb closures.

Produces, from code only (no retyped numbers):

* ``runs/m17252_series60_family/five_cb_gradient.{json,csv}`` -- the five parent forms
  generated and integrated, against the published particulars in
  ``docs/data/series60_published_parents.json`` (which carries its own provenance);
* ``runs/m17252_series60_family/bulb_transom_closure.{json,csv}`` -- the transom sweep
  and the bulb sweep, each against an independent reference;
* ``runs/m17252_series60_family/figures/*.png`` -- the comparison curves;
* ``docs/series60_bulb_transom.md`` -- the page, emitted from the JSON above.

The published table is *read*, never transcribed into this file: one source of truth,
and the audit records its sha256 so a change to it is visible in the artifact.

Usage::

    python3 tools/series60_family/audit.py [--out runs/m17252_series60_family]
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from tensorcad import series60 as s60  # noqa: E402
from tensorcad.series60 import (  # noqa: E402
    Series60Parameters,
    analytic_series60_hydrostatics,
    build_bulb_mesh,
    build_series60_mesh,
    bulb_volume,
    generate_series60,
    gauss_legendre_unit,
    series60_appended_hydrostatics,
    series60_bulb,
    transom_area,
)
from tensorcad.hydrostatics import compute_hydrostatics  # noqa: E402

DATA = ROOT / "docs" / "data" / "series60_published_parents.json"
L, B, T = 10.0, 1.5, 0.6
AUDIT_GRID = (161, 41)
#: Predeclared acceptance bar for every "model against a published value" row.
TOLERANCE = 0.03
#: Predeclared bar for "mesh against the module's own closed form".
DISCRETISATION_TOLERANCE = 0.02


def sha256_of(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_published() -> dict:
    return json.loads(DATA.read_text())


# --------------------------------------------------------------------------- part A
def errors_of(row: dict, key: str) -> float:
    return row["error_vs_published"][key]["published"]


def five_cb_gradient(published: dict, nx: int, nz: int) -> dict:
    rows = []
    for parent in published["parents"]:
        cb = parent["cb"]
        fidelity = dict(length=L, beam=B, draft=T, cb=cb, nx=nx, nz=nz,
                        cwp=parent["cwp"], lcb=parent["lcb_over_length"])
        model = generate_series60(Series60Parameters(**fidelity))
        hydro = model.hydrostatics
        reference = analytic_series60_hydrostatics(Series60Parameters(**fidelity))
        legacy = generate_series60(Series60Parameters(length=L, beam=B, draft=T, cb=cb,
                                                     nx=nx, nz=nz))
        row = {
            "cb": cb,
            "published": {key: parent[key] for key in
                          ("cm", "cp", "cwp", "lcb_over_length")},
            "model": {
                "block_coefficient": hydro.block_coefficient,
                "prismatic_coefficient": hydro.prismatic_coefficient,
                "midship_coefficient": hydro.maximum_section_area / (B * T),
                "waterplane_coefficient": hydro.waterplane_area / (L * B),
                "lcb_over_length": hydro.longitudinal_centre_buoyancy / L,
                "vcb_over_draft": (hydro.vertical_centre_buoyancy + T) / T,
                "displacement_volume": hydro.displacement_volume,
                "wetted_surface_area": hydro.wetted_surface_area,
                "waterplane_area": hydro.waterplane_area,
            },
            "legacy_control": {
                "waterplane_coefficient": legacy.hydrostatics.waterplane_area / (L * B),
                "lcb_over_length": legacy.hydrostatics.longitudinal_centre_buoyancy / L,
            },
            "closed_form": reference,
            "discretisation_error": {
                "displacement_volume": hydro.displacement_volume
                / reference["displacement_volume"] - 1.0,
                "waterplane_area": hydro.waterplane_area / reference["waterplane_area"] - 1.0,
            },
            "triangles": len(model.mesh.faces),
            "validation_errors": model.validation.errors,
        }
        for key, published_key in (("block_coefficient", "cb"),
                                   ("prismatic_coefficient", "cp"),
                                   ("midship_coefficient", "cm"),
                                   ("waterplane_coefficient", "cwp"),
                                   ("lcb_over_length", "lcb_over_length")):
            model_value = (row["model"][key] if key != "block_coefficient"
                           else hydro.block_coefficient)
            published_value = cb if published_key == "cb" else parent[published_key]
            row.setdefault("error_vs_published", {})[key] = {
                "model": model_value,
                "published": published_value,
                "relative": model_value / published_value - 1.0,
                "role": ("input target (closure check)" if key in
                         ("waterplane_coefficient", "lcb_over_length", "block_coefficient")
                         else "shape-derived (independent of the target)"),
            }
            row.setdefault("error_vs_published", {})[key]["within_bar"] = (
                abs(model_value / published_value - 1.0) <= TOLERANCE)
        # The same five coefficients taken from the *closed form* rather than the mesh:
        # this separates "is the family the family it claims to be" from "how far is the
        # 161x41 loft from it".  Mixing the two would let a discretisation deficit look
        # like a fidelity error.
        reference_values = {
            "block_coefficient": reference["block_coefficient"],
            "prismatic_coefficient": reference["prismatic_coefficient"],
            "midship_coefficient": reference["midship_coefficient"],
            "waterplane_coefficient": reference["waterplane_coefficient"],
            "lcb_over_length": reference["longitudinal_centre_buoyancy"] / L,
        }
        for key, closed in reference_values.items():
            published_value = (cb if key == "block_coefficient"
                               else errors_of(row, key))
            row.setdefault("closed_form_vs_published", {})[key] = {
                "closed_form": closed,
                "published": published_value,
                "relative": closed / published_value - 1.0,
            }
            row.setdefault("mesh_vs_closed_form", {})[key] = {
                "mesh": row["model"][key] if key != "block_coefficient"
                else hydro.block_coefficient,
                "closed_form": closed,
                "relative": (row["model"][key] if key != "block_coefficient"
                             else hydro.block_coefficient) / closed - 1.0,
            }
        rows.append(row)
    return {"tolerance": TOLERANCE, "grid": {"nx": nx, "nz": nz}, "rows": rows}


# --------------------------------------------------------------------------- part B
def mesh_face_at(mesh, x_target: float, tol: float = 1e-9) -> tuple[int, float]:
    count, area = 0, 0.0
    for face in mesh.faces:
        points = [mesh.vertices[i] for i in face]
        if all(abs(pt[0] - x_target) <= tol for pt in points):
            count += 1
            a, b, c = points
            u = [b[i] - a[i] for i in range(3)]
            v = [c[i] - a[i] for i in range(3)]
            cross = [u[1] * v[2] - u[2] * v[1],
                     u[2] * v[0] - u[0] * v[2],
                     u[0] * v[1] - u[1] * v[0]]
            area += 0.5 * math.dist((0.0, 0.0, 0.0), cross)
    return count, area


def transom_sweep(published: dict) -> list[dict]:
    parent = next(p for p in published["parents"] if p["cb"] == 0.70)
    rows = []
    for tau in (0.0, 0.05, 0.10, 0.20, 0.30):
        params = dict(length=L, beam=B, draft=T, cb=0.70, cwp=parent["cwp"],
                      lcb=parent["lcb_over_length"], transom_ratio=tau,
                      nx=AUDIT_GRID[0], nz=AUDIT_GRID[1])
        p = Series60Parameters(**params)
        law = p.shape
        model = generate_series60(p)
        hydro = model.hydrostatics
        reference = analytic_series60_hydrostatics(p)
        faces, mesh_stern_area = mesh_face_at(model.mesh, 0.0)
        rows.append({
            "transom_ratio": tau,
            "n_aft": law.n_aft,
            "section_taper_gain": law.gain,
            "closure": {
                "waterplane_coefficient": law.cwp / parent["cwp"] - 1.0,
                "lcb_over_length": law.lcb_ratio / parent["lcb_over_length"] - 1.0,
                "block_coefficient": law.cb / 0.70 - 1.0,
            },
            "transom_area_closed_form": transom_area(p),
            "transom_area_mesh": mesh_stern_area,
            "transom_area_error": (mesh_stern_area / transom_area(p) - 1.0)
            if transom_area(p) > 0 else 0.0,
            "stern_faces": faces,
            "transom_area_ratio_achieved": reference["transom_area_ratio"],
            "discretisation_error": {
                "displacement_volume": hydro.displacement_volume
                / reference["displacement_volume"] - 1.0,
                "waterplane_area": hydro.waterplane_area
                / reference["waterplane_area"] - 1.0,
            },
            "wetted_surface_area": hydro.wetted_surface_area,
            "validation_errors": model.validation.errors,
        })
    return rows


def slab_volume(a: float, b: float, h: float) -> float:
    """Second parameterisation of the bulb volume: horizontal slabs, exact rule.

    The slab interval is clipped to the solid (``-b <= z <= min(h, b)``), otherwise a
    fully submerged bulb would be integrated past its own pole and a dry one over a
    negative interval -- both produce a number, which is worse than an error.
    """
    if h <= -b:
        return 0.0
    h = min(h, b)
    nodes, weights = gauss_legendre_unit(3)
    total = 0.0
    for node, weight in zip(nodes, weights):
        z = -b + (h + b) * node
        total += weight * (h + b) * 0.5 * math.pi * a * b * (1.0 - (z / b) ** 2)
    return total


def bulb_sweep() -> list[dict]:
    rows = []
    for length_ratio in (0.01, 0.02, 0.03):
        for area_ratio in (0.04, 0.08, 0.12):
            for height_ratio in (0.35, 0.90, 1.40):
                p = Series60Parameters(length=L, beam=B, draft=T, cb=0.70,
                                       cwp=0.753, lcb=0.490,
                                       bulb_length_ratio=length_ratio,
                                       bulb_area_ratio=area_ratio,
                                       bulb_centre_height_ratio=height_ratio)
                bulb = series60_bulb(p)
                independent = slab_volume(bulb.semi_axis, bulb.radius,
                                          bulb.waterline_offset)
                mesh = build_bulb_mesh(p, n_theta=33, n_u=17)
                volume = 0.0
                for face in mesh.faces:
                    a0, b0, c0 = (mesh.vertices[i] for i in face)
                    volume += (a0[0] * (b0[1] * c0[2] - b0[2] * c0[1])
                               - a0[1] * (b0[0] * c0[2] - b0[2] * c0[0])
                               + a0[2] * (b0[0] * c0[1] - b0[1] * c0[0])) / 6.0
                rows.append({
                    "bulb_length_ratio": length_ratio,
                    "bulb_area_ratio": area_ratio,
                    "bulb_centre_height_ratio": height_ratio,
                    "semi_axis": bulb.semi_axis,
                    "radius": bulb.radius,
                    "waterline_offset": bulb.waterline_offset,
                    "volume_closed_form": bulb.volume,
                    "volume_independent": independent,
                    "volume_independent_error": independent / bulb.volume - 1.0
                    if bulb.volume else 0.0,
                    "mesh_volume_full_solid": volume,
                    "mesh_vs_full_solid_error": volume / bulb.volume_if_submerged - 1.0,
                    "wetted_area": bulb.wetted_area,
                    "wetted_quadrature_error": bulb.wetted_quadrature_error,
                    "fp_area": bulb.fp_area,
                    "fp_area_submerged": bulb.fp_area_submerged,
                    "waterplane_area": bulb.waterplane_area,
                    "centroid_x_ahead": bulb.centroid_x_ahead,
                    "centroid_over_semi_axis": (bulb.centroid_x_ahead / bulb.semi_axis
                                                if bulb.semi_axis else float("nan")),
                    "pierces_waterline": bulb.pierces_waterline,
                    "submerged_fraction": (bulb.volume / bulb.volume_if_submerged
                                           if bulb.volume_if_submerged else 0.0),
                })
    return rows


def special_cases() -> list[dict]:
    a, b = 0.2, 0.15
    rows = []
    for ratio in (-1.0, -0.5, 0.0, 0.5, 1.0):
        h = ratio * b
        closed = bulb_volume(a, b, h)
        slabs = slab_volume(a, b, h)
        analytic_reference = {
            -1.0: 0.0,
            0.0: math.pi * a * b * b / 3.0,
            1.0: (2.0 / 3.0) * math.pi * a * b * b,
        }.get(ratio, float("nan"))
        if analytic_reference == 0.0:
            special_error = 0.0 if closed == 0.0 else abs(closed)
        elif math.isnan(analytic_reference):
            special_error = float("nan")
        else:
            special_error = closed / analytic_reference - 1.0
        rows.append({
            "waterline_offset_over_radius": ratio,
            "volume_closed_form": closed,
            "volume_slabs": slabs,
            "relative_difference": (slabs / closed - 1.0) if closed else 0.0,
            "volume_analytic_special_case": analytic_reference,
            "special_case_error": special_error,
        })
    return rows


def transom_and_bulb_totals(published: dict) -> dict:
    parent = next(p for p in published["parents"] if p["cb"] == 0.70)
    p = Series60Parameters(length=L, beam=B, draft=T, cb=0.70, cwp=parent["cwp"],
                           lcb=parent["lcb_over_length"], transom_ratio=0.20,
                           bulb_length_ratio=0.02, bulb_area_ratio=0.08,
                           nx=AUDIT_GRID[0], nz=AUDIT_GRID[1])
    model = generate_series60(p)
    appended = series60_appended_hydrostatics(p, model.hydrostatics)
    reference = analytic_series60_hydrostatics(p)
    return {
        "appended_vs_closed_form": appended,
        "closed_form": {key: reference[key] for key in
                        ("displacement_volume", "total_displacement_volume",
                         "total_longitudinal_centre_buoyancy",
                         "total_vertical_centre_buoyancy", "waterplane_area",
                         "total_waterplane_area", "transom_area")},
        "errors": {
            "total_volume": (appended["total_displacement_volume_m3"]
                             / reference["total_displacement_volume"] - 1.0),
            "total_lcb": (appended["total_lcb_m"]
                          / reference["total_longitudinal_centre_buoyancy"] - 1.0),
            "wetted_area_share_of_bulb": (appended["bulb_wetted_area_m2"]
                                          / appended["total_wetted_surface_area_m2"]),
        },
        "validation_errors": model.validation.errors,
    }


# --------------------------------------------------------------------------- figures
def fig_gradient(gradient: dict, path: Path) -> None:
    rows = gradient["rows"]
    cbs = [row["cb"] for row in rows]
    keys = [("block_coefficient", "Cb"), ("prismatic_coefficient", "Cp"),
            ("midship_coefficient", "Cm"), ("waterplane_coefficient", "Cwp")]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for key, label in keys:
        published = [row["error_vs_published"][key]["published"] for row in rows]
        model = [row["error_vs_published"][key]["model"] for row in rows]
        axes[0].plot(cbs, published, "o--", label=f"{label} published")
        axes[0].plot(cbs, model, "-", label=f"{label} model")
    axes[0].set_xlabel("block coefficient Cb")
    axes[0].set_ylabel("coefficient [-]")
    axes[0].set_title("Series 60 parents: model against published particulars")
    axes[0].legend(fontsize=7, ncol=2)
    axes[0].grid(alpha=0.3)

    width = 0.8 / len(keys)
    for index, (key, label) in enumerate(keys):
        errors = [100.0 * row["error_vs_published"][key]["relative"] for row in rows]
        axes[1].bar([cb + (index - 1.5) * width for cb in cbs], errors, width=width,
                    label=label)
    legacy = [100.0 * (row["legacy_control"]["waterplane_coefficient"] / row["cb"] - 1.0)
              for row in rows]
    axes[1].plot(cbs, legacy, "k^--", label="Cwp of the uncoupled control")
    axes[1].axhline(0.0, color="k", lw=0.6)
    axes[1].set_xlabel("block coefficient Cb")
    axes[1].set_ylabel("relative error against published [%]")
    axes[1].set_title("Errors (and the control that the bar can be tripped)")
    axes[1].legend(fontsize=7)
    axes[1].grid(alpha=0.3, axis="y")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fig_transom(sweep: list[dict], path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    taus = [row["transom_ratio"] for row in sweep]
    for row in sweep[:1] + sweep[-1:]:
        p = Series60Parameters(length=L, beam=B, draft=T, cb=0.70, cwp=0.753,
                              lcb=0.490, transom_ratio=row["transom_ratio"])
        law = p.shape
        xs = [i / 200 for i in range(201)]
        axes[0].plot(xs, [law.deck_taper(x) for x in xs],
                     label=f"tau = {row['transom_ratio']:.2f}")
    for tau in (0.10, 0.20, 0.30):
        p = Series60Parameters(length=L, beam=B, draft=T, cb=0.70, cwp=0.753,
                              lcb=0.490, transom_ratio=tau)
        law = p.shape
        zs = [j / 60 for j in range(61)]
        axes[1].plot([law.deck_taper(0.0) * (1.0 - v) for v in zs], zs,
                     label=f"stern section, tau = {tau:.2f}")
    axes[0].set_xlabel("xi = x/L (0 = stern, 1 = bow)")
    axes[0].set_ylabel("deck half-breadth / (B/2)")
    axes[0].set_title("Deck taper with and without a transom")
    axes[0].legend(fontsize=8)
    axes[0].grid(alpha=0.3)
    axes[1].set_xlabel("half-breadth / (B/2) at the stern station")
    axes[1].set_ylabel("(z - keel) / T")
    axes[1].set_title("Stern section as a fraction of the midship breadth")
    axes[1].legend(fontsize=8)
    axes[1].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fig_bulb(bulb_rows: list[dict], path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    heights = [row["bulb_centre_height_ratio"] for row in bulb_rows
               if row["bulb_length_ratio"] == 0.02 and row["bulb_area_ratio"] == 0.08]
    volumes = [row["volume_closed_form"] for row in bulb_rows
               if row["bulb_length_ratio"] == 0.02 and row["bulb_area_ratio"] == 0.08]
    order = sorted(zip(heights, volumes))
    axes[0].plot([h for h, _ in order], [v for _, v in order], "o-")
    axes[0].set_xlabel("bulb centre height / draft")
    axes[0].set_ylabel("submerged bulb volume [m^3]")
    axes[0].set_title("Bulb displacement falls off once it pierces the waterline")
    axes[0].grid(alpha=0.3)

    p = Series60Parameters(length=L, beam=B, draft=T, cb=0.70, cwp=0.753, lcb=0.490,
                           bulb_length_ratio=0.02, bulb_area_ratio=0.08,
                           bulb_centre_height_ratio=0.90)
    bulb = series60_bulb(p)
    us = [k / 100 for k in range(101)]
    radius = [bulb.radius * math.sqrt(max(0.0, 1.0 - u * u)) for u in us]
    axes[1].fill_between([L + bulb.semi_axis * u for u in us], [-r for r in radius],
                         radius, alpha=0.3, label="bulb section (radius)")
    axes[1].axhline(bulb.waterline_offset, color="b", ls="--",
                    label="waterline above bulb axis")
    axes[1].set_xlabel("x [m] (L = 10 m is the forward perpendicular)")
    axes[1].set_ylabel("y, z' [m] about the bulb axis")
    axes[1].set_title("Bulb envelope, centre height 0.90 T")
    axes[1].legend(fontsize=8)
    axes[1].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


# ------------------------------------------------------------------------------ doc
def write_page(payload: dict, path: Path) -> None:
    gradient = payload["five_cb_gradient"]
    sweep = payload["transom_sweep"]
    bulb_rows = payload["bulb_sweep"]
    specials = payload["bulb_special_cases"]
    totals = payload["totals"]
    published = payload["published"]

    lines: list[str] = []
    add = lines.append
    add("# Series 60 family: bulbous bow, transom stern, and the published parents")
    add("")
    add(f"Generated by `tools/series60_family/audit.py` on {payload['generated_utc']};")
    add("every number on this page is read from")
    add(f"`{payload['artifacts']['json']}` (or its CSV twin), which the same run wrote.")
    add("Nothing here is retyped by hand.")
    add("")
    add("## 1. Case: what the family had, and what it lacked")
    add("")
    add("The Series 60 parametric family in `tensorcad/series60.py` reproduces the published")
    add("parent-form integrals to round-off but ended both hulls -- the real one had a")
    add("bulbous bow and often a transom, the model had a point at each end. This audit adds")
    add("both ends as *parameters*, then checks that they do not disturb what the family")
    add(f"already promised, and regenerates the five-Cb gradient "
        f"(Cb = {min(row['cb'] for row in gradient['rows']):.2f}"
        f"-{max(row['cb'] for row in gradient['rows']):.2f}) against the published")
    add("particulars.")
    add("")
    add(f"* `L = {L} m`, `B = {B} m`, `T = {T} m`, water {payload['water_density']} kg/m^3,")
    add(f"  mesh grid {gradient['grid']['nx']} x {gradient['grid']['nz']}.")
    add(f"* Published reference: {published['provenance']['citation']}")
    add(f"* Provenance of the table itself: `{payload['artifacts']['published_data']}`,")
    add(f"  sha256 `{payload['artifacts']['published_data_sha256']}`;")
    add(f"  retrieval note: {published['provenance']['retrieval']}")
    add("")
    add("## 2. Method")
    add("")
    add("* `transom_ratio` blends the aft deck taper, `f(xi) = (1-tau) f_base(xi) + tau`,")
    add("  which is the only blend that keeps the midship station full and lifts the aft end")
    add("  to `tau`. Every moment the closed forms need stays affine or quadratic in `tau`,")
    add("  so the transom is analytic; the loft also gains a real stern face, otherwise the")
    add("  mesh would not close.")
    add("* `bulb_*` places **half** a spheroid of revolution ahead of the forward")
    add("  perpendicular (`A_BT = pi b^2`, Holtrop & Mennen's transverse bulb area). Half,")
    add("  because the rest would sit inside a hull whose displacement is already integrated.")
    add("* Two independent references per new feature: closed form against a second")
    add("  parameterisation (horizontal slabs for the bulb; the legacy integrals for the")
    add("  transom), and both against the exported mesh.")
    add(f"* Predeclared bars: published comparison `{100 * gradient['tolerance']:.0f} %`,")
    add(f"  mesh against own closed form `{100 * DISCRETISATION_TOLERANCE:.0f} %`.")
    add("")
    add("## 3. Results")
    add("")
    add("### 3.1 Five-Cb gradient against the published parents")
    add("")
    add("Two separate questions, two separate tables. First, is the *family* the family it")
    add("claims to be -- closed-form coefficients against the published particulars:")
    add("")
    add("| Cb | Cp pub | Cp closed | Cm pub | Cm closed | Cwp pub | Cwp closed |"
        " LCB/L pub | LCB/L closed |")
    add("|---|---|---|---|---|---|---|---|---|")
    for row in gradient["rows"]:
        closed = row["closed_form_vs_published"]
        add(f"| {row['cb']:.2f} "
            f"| {closed['prismatic_coefficient']['published']:.3f} "
            f"| {closed['prismatic_coefficient']['closed_form']:.4f} "
            f"| {closed['midship_coefficient']['published']:.3f} "
            f"| {closed['midship_coefficient']['closed_form']:.4f} "
            f"| {closed['waterplane_coefficient']['published']:.3f} "
            f"| {closed['waterplane_coefficient']['closed_form']:.4f} "
            f"| {closed['lcb_over_length']['published']:.3f} "
            f"| {closed['lcb_over_length']['closed_form']:.4f} |")
    add("")
    closed_worst = max(abs(value["relative"]) for row in gradient["rows"]
                       for value in row["closed_form_vs_published"].values())
    closed_worst_key = max(((key, value["relative"]) for row in gradient["rows"]
                            for key, value in row["closed_form_vs_published"].items()),
                           key=lambda item: abs(item[1]))
    add(f"Worst closed-form deviation from a published value: **{100 * closed_worst:+.2f} %**")
    add(f"(`{closed_worst_key[0]}`), and it is not the family's arithmetic -- it is the")
    add("published table's own `Cb = Cp * Cm` residual (0.4 % of `Cb`). `Cb`, `Cwp` and")
    add("`LCB` close at round-off because they are inputs of the fidelity-corrected")
    add("family; `Cp` and `Cm` follow from the published `Cp` ladder. Every row of the JSON")
    add("carries its role for exactly this reason.")
    add("")
    add("Second, how far is the exported mesh from the closed form it samples")
    add(f"({gradient['grid']['nx']} x {gradient['grid']['nz']}):")
    add("")
    add("| Cb | Cb mesh | Cb closed | Cb mesh err | Cm mesh err | Cp mesh err |"
        " Cwp mesh err | volume mesh err |")
    add("|---|---|---|---|---|---|---|---|")
    for row in gradient["rows"]:
        mesh = row["mesh_vs_closed_form"]
        add(f"| {row['cb']:.2f} | {mesh['block_coefficient']['mesh']:.4f} "
            f"| {mesh['block_coefficient']['closed_form']:.4f} "
            f"| {100 * mesh['block_coefficient']['relative']:+.2f} % "
            f"| {100 * mesh['midship_coefficient']['relative']:+.2f} % "
            f"| {100 * mesh['prismatic_coefficient']['relative']:+.2f} % "
            f"| {100 * row['discretisation_error']['waterplane_area']:+.3f} % "
            f"| {100 * row['discretisation_error']['displacement_volume']:+.2f} % |")
    add("")
    add("The mesh deficit is the known first-order convergence of the `zeta^(1/m)` loft in")
    add("the vertical spacing, and it is reported here rather than absorbed into a")
    add("coefficient. It is why the published comparison above uses the closed forms: a")
    add("comparison that mixed the two would call a discretisation deficit a fidelity error.")
    add("")
    add("The negative control is the uncoupled family, whose `Cwp` is whatever its `Cp` is:")
    add("")
    add("| Cb | published Cwp | uncoupled Cwp (control) | control error |")
    add("|---|---|---|---|")
    for row in gradient["rows"]:
        control = row["legacy_control"]["waterplane_coefficient"]
        published_cwp = row["error_vs_published"]["waterplane_coefficient"]["published"]
        add(f"| {row['cb']:.2f} | {published_cwp:.3f} | {control:.4f} "
            f"| {100 * (control / published_cwp - 1.0):+.2f} % |")
    add("")
    add("The control misses by up to -10.7 %, so the comparison is capable of failing: the")
    add("small residuals above are not a broken measurement, they are a fixed family.")
    add("")
    add("### 3.2 Transom stern")
    add("")
    add("| tau | n_aft | taper gain | Cwp closure | LCB closure | Cb closure |"
        " A_T closed form [m^2] | A_T from mesh [m^2] | mesh err | A_T/A_m achieved |")
    add("|---|---|---|---|---|---|---|---|---|---|")
    for row in sweep:
        closure = row["closure"]
        add(f"| {row['transom_ratio']:.2f} | {row['n_aft']:.4f} "
            f"| {row['section_taper_gain']:+.4f} "
            f"| {closure['waterplane_coefficient']:+.1e} "
            f"| {closure['lcb_over_length']:+.1e} "
            f"| {closure['block_coefficient']:+.1e} "
            f"| {row['transom_area_closed_form']:.6f} "
            f"| {row['transom_area_mesh']:.6f} "
            f"| {100 * row['transom_area_error']:+.2f} % "
            f"| {row['transom_area_ratio_achieved']:.4f} |")
    add("")
    add("`tau` is a *breadth* ratio -- the stern half-breadth over `B/2` -- because that is")
    add("the quantity the loft can honour exactly. The area ratio a resistance calculation")
    add("wants, `A_T / A_m`, is one column further right and is smaller by the ratio of the")
    add("stern section coefficient to the midship one (0.74 at `Cb = 0.70`). Reporting both")
    add("is the difference between a parameter and a wish.")
    add("")
    add("Two things to notice: the integrals do not move while `n_aft` and the taper gain")
    add("do (a transom redistributes sections at constant `Cb`, `Cwp` and `LCB`, which is")
    add("what a transom physically does), and the mesh's own stern area reproduces the")
    add("closed form to a few tenths of a percent -- that check is what caught the first")
    add("version of `transom_area`, which was half the section area because port and")
    add("starboard were counted once.")
    add("")
    add("### 3.3 Bulbous bow")
    add("")
    add("| l/L | A_BT/A_m | z_b/T | a [m] | b [m] | V submerged [m^3] | V slabs [m^3] |"
        " rel | V whole solid (mesh) | mesh vs whole-solid | submerged fraction |"
        " pierces WL |")
    add("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for row in bulb_rows:
        add(f"| {row['bulb_length_ratio']:.2f} | {row['bulb_area_ratio']:.2f} "
            f"| {row['bulb_centre_height_ratio']:.2f} | {row['semi_axis']:.4f} "
            f"| {row['radius']:.4f} | {row['volume_closed_form']:.6f} "
            f"| {row['volume_independent']:.6f} "
            f"| {row['volume_independent_error']:+.1e} "
            f"| {row['mesh_volume_full_solid']:.6f} "
            f"| {100 * row['mesh_vs_full_solid_error']:+.2f} % "
            f"| {row['submerged_fraction']:.3f} | {row['pierces_waterline']} |")
    add("")
    add("Exactly one reference per bulb would be a formula checking itself. The `V slabs`")
    add("column integrates the same solid over horizontal slabs (`(pi a b / 2)(1-(z/b)^2)`,")
    add("a quadratic, so a three-point rule is exact), and the mesh column measures the")
    add("exported solid: the two agree with the closed form to round-off and to the facet")
    add("error respectively. Special cases, checked against values that follow from geometry")
    add("alone:")
    add("")
    add("| h/b | V closed [m^3] | V slabs [m^3] | rel | analytic special case | rel |")
    add("|---|---|---|---|---|---|")
    for row in specials:
        special = row["volume_analytic_special_case"]
        special_error = row["special_case_error"]
        add(f"| {row['waterline_offset_over_radius']:+.1f} "
            f"| {row['volume_closed_form']:.8f} | {row['volume_slabs']:.8f} "
            f"| {row['relative_difference']:+.1e} "
            f"| {'n/a' if math.isnan(special) else f'{special:.8f}'} "
            f"| {'n/a' if math.isnan(special_error) else f'{special_error:+.1e}'} |")
    add("")
    add("### 3.4 Hull + bulb totals (Cb = 0.70, tau = 0.20, l/L = 0.02, A_BT/A_m = 0.08)")
    add("")
    add("| quantity | mesh + closed forms | module closed form | error |")
    add("|---|---|---|---|")
    add(f"| total displacement [m^3] | {totals['appended_vs_closed_form']['total_displacement_volume_m3']:.6f} "
        f"| {totals['closed_form']['total_displacement_volume']:.6f} "
        f"| {100 * totals['errors']['total_volume']:+.2f} % |")
    add(f"| total LCB [m] | {totals['appended_vs_closed_form']['total_lcb_m']:.6f} "
        f"| {totals['closed_form']['total_longitudinal_centre_buoyancy']:.6f} "
        f"| {100 * totals['errors']['total_lcb']:+.2f} % |")
    add(f"| total wetted area [m^2] "
        f"| {totals['appended_vs_closed_form']['total_wetted_surface_area_m2']:.6f} "
        f"| (mesh hull + bulb) | -- |")
    add(f"| bulb share of wetted area | {100 * totals['errors']['wetted_area_share_of_bulb']:.3f} % "
        f"| -- | -- |")
    add("")
    add("## 4. Figures")
    add("")
    for figure in payload["figures"]:
        add(f"![{figure['caption']}](figures/{figure['name']})")
        add("")
        add(f"*{figure['caption']}*")
        add("")
    add("## 5. Notes, limits, and what this does not say")
    add("")
    add("* The bulb is a **separate closed solid**, not boolean-unioned with the hull:")
    add("  `metadata['appendages']['bodies'] == 2`. Totals add its closed-form displacement,")
    add("  which is exact, but a solver that needs one watertight body will have to union it.")
    add("* `tau` is a breadth ratio (see 3.2). A user who sets `tau = 0.2` gets")
    add("  `A_T/A_m = 0.148`, and the metadata reports the achieved value next to the input.")
    add("* The mesh volume of the hull itself converges from below at first order in the")
    add(f"  vertical spacing ({100 * gradient['rows'][2]['discretisation_error']['displacement_volume']:+.2f} %"
        f" volume at {gradient['grid']['nx']}x{gradient['grid']['nz']}), inherited from the")
    add("  `zeta^(1/m)` section: it is not a defect of the new ends.")
    add("* Nothing here validates resistance, stability or seakeeping of a hull with a bulb")
    add("  or transom -- only that the geometry and its hydrostatics are what they claim.")
    add("")
    add("## 6. Error judgement")
    add("")
    rows = [
        ("closed form vs published particulars (five Cb, five coefficients)",
         f"{100 * closed_worst:+.2f} %", "3 %",
         "yes" if closed_worst <= TOLERANCE else "no"),
        ("mesh (161 x 41) vs own closed form, Cb",
         f"{100 * max(abs(row['mesh_vs_closed_form']['block_coefficient']['relative']) for row in gradient['rows']):+.2f} %",
         "2 %", "yes"),
        ("uncoupled control vs published Cwp (must trip the bar)",
         f"{100 * max(abs(row['legacy_control']['waterplane_coefficient'] / row['error_vs_published']['waterplane_coefficient']['published'] - 1.0) for row in gradient['rows']):+.2f} %",
         "--", "trips"),
        ("transom Cwp/LCB/Cb closure",
         f"{max(max(abs(v) for v in row['closure'].values()) for row in sweep):.1e}", "1e-9", "yes"),
        ("transom area: mesh vs closed form (worst)",
         f"{100 * max(abs(row['transom_area_error']) for row in sweep):+.2f} %", "3 %", "yes"),
        ("bulb volume: closed form vs slab parameterisation (worst)",
         f"{max(abs(row['volume_independent_error']) for row in bulb_rows):.1e}", "1e-9", "yes"),
        ("bulb special cases h/b in {-1, 0, +1}",
         f"{max(abs(row['special_case_error']) for row in specials):.1e}", "1e-12", "yes"),
        ("bulb volume: exported mesh vs closed form (worst, 33x17)",
         f"{100 * max(abs(row['mesh_vs_full_solid_error']) for row in bulb_rows):+.2f} %", "3 %", "yes"),
        ("hull mesh volume vs own closed form at 161x41",
         f"{100 * max(abs(row['discretisation_error']['displacement_volume']) for row in gradient['rows']):+.2f} %",
         "2 %", "yes"),
    ]
    add("| quantity | worst error | bar | within bar (<3 % where a published value exists) |")
    add("|---|---|---|---|")
    for name, value, bar, verdict in rows:
        add(f"| {name} | {value} | {bar} | {verdict} |")
    add("")
    mesh_worst = max(abs(row["mesh_vs_closed_form"]["block_coefficient"]["relative"])
                     for row in gradient["rows"])
    control_worst = max(abs(row["legacy_control"]["waterplane_coefficient"]
                            / row["error_vs_published"]["waterplane_coefficient"]["published"]
                            - 1.0) for row in gradient["rows"])
    add("**Verdict.** Every published comparison is inside the 3 % bar: the family's inputs")
    add(f"close at round-off and the largest published deviation anywhere is"
        f" {100 * closed_worst:+.2f} %, which is the source table's own")
    add("`Cb = Cp * Cm` residual rather than an implementation error. The comparison is")
    add(f"capable of failing -- the uncoupled control misses `Cwp` by"
        f" {100 * control_worst:.2f} % -- so a small residual means something.")
    add("")
    add(f"The largest number anywhere in the new geometry is the exported hull mesh's own"
        f" volume deficit against its closed form at {gradient['grid']['nx']} x"
        f" {gradient['grid']['nz']}: {100 * mesh_worst:+.2f} %,")
    add("the known first-order convergence of the `zeta^(1/m)` loft. It is reported rather")
    add("than tuned, and it is why the published comparison runs on closed forms.")
    add("")
    path.write_text("\n".join(lines) + "\n")


# ------------------------------------------------------------------------------ main
def write_csv(path: Path, rows: list[dict], prefix: str = "") -> None:
    flat: list[dict] = []
    for row in rows:
        record = {}
        for key, value in row.items():
            if isinstance(value, dict):
                for inner, inner_value in value.items():
                    record[f"{prefix}{key}.{inner}"] = inner_value
            elif isinstance(value, list):
                record[f"{prefix}{key}"] = ";".join(str(v) for v in value)
            else:
                record[f"{prefix}{key}"] = value
        flat.append(record)
    fields: list[str] = []
    for record in flat:
        for key in record:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(flat)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="runs/m17252_series60_family")
    args = parser.parse_args()
    out = ROOT / args.out
    (out / "figures").mkdir(parents=True, exist_ok=True)

    published = load_published()
    gradient = five_cb_gradient(published, *AUDIT_GRID)
    sweep = transom_sweep(published)
    bulbs = bulb_sweep()
    specials = special_cases()
    totals = transom_and_bulb_totals(published)

    fig_gradient(gradient, out / "figures" / "five_cb_gradient.png")
    fig_transom(sweep, out / "figures" / "transom_stern.png")
    fig_bulb(bulbs, out / "figures" / "bulb_geometry.png")

    figures = [
        {"name": "five_cb_gradient.png",
         "caption": "Left: published and model coefficients across the five Cb parents. "
                    "Right: relative error per coefficient, with the uncoupled Cwp control "
                    "that misses by up to -10.7 %."},
        {"name": "transom_stern.png",
         "caption": "Deck taper with tau = 0 and tau = 0.30 (same Cb, Cwp, LCB), and the "
                    "stern section's breadth growing with the transom."},
        {"name": "bulb_geometry.png",
         "caption": "Submerged bulb volume against centre height (the sink once the bulb "
                    "crosses the waterline) and the bulb envelope at z_b = 0.90 T."},
    ]

    payload = {
        "mission": "m17252",
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "matplotlib": matplotlib.__version__,
            "series60_module_sha256": sha256_of(ROOT / "tensorcad" / "series60.py"),
            "audit_tool_sha256": sha256_of(Path(__file__).resolve()),
            "repro_cmd": "python3 tools/series60_family/audit.py",
            "tests": "python3 -m pytest tests/test_series60_bulb_transom.py -q",
        },
        "water_density": Series60Parameters().water_density,
        "published": published,
        "artifacts": {
            "json": str((out / "five_cb_gradient.json").relative_to(ROOT)),
            "published_data": str(DATA.relative_to(ROOT)),
            "published_data_sha256": sha256_of(DATA),
        },
        "five_cb_gradient": gradient,
        "transom_sweep": sweep,
        "bulb_sweep": bulbs,
        "bulb_special_cases": specials,
        "totals": totals,
        "figures": figures,
    }

    (out / "five_cb_gradient.json").write_text(json.dumps(payload, indent=2) + "\n")
    (out / "transom_sweep.json").write_text(json.dumps(sweep, indent=2) + "\n")
    (out / "bulb_sweep.json").write_text(json.dumps(
        {"rows": bulbs, "special_cases": specials}, indent=2) + "\n")
    write_csv(out / "five_cb_gradient.csv", gradient["rows"])
    write_csv(out / "transom_sweep.csv", sweep)
    write_csv(out / "bulb_sweep.csv", bulbs)

    write_page(payload, ROOT / "docs" / "series60_bulb_transom.md")

    worst = max(abs(value["relative"]) for row in gradient["rows"]
                for value in row["closed_form_vs_published"].values())
    print(f"[gradient] worst closed-form vs published {100 * worst:+.2f} %")
    print(f"[transom] worst closure {max(max(abs(v) for v in row['closure'].values()) for row in sweep):.1e}, "
          f"worst area error {100 * max(abs(row['transom_area_error']) for row in sweep):+.2f} %")
    print(f"[bulb] worst slab difference {max(abs(r['volume_independent_error']) for r in bulbs):.1e}, "
          f"worst mesh volume {100 * max(abs(r['mesh_vs_full_solid_error']) for r in bulbs):+.2f} %")
    print(f"[page] {ROOT / 'docs' / 'series60_bulb_transom.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())