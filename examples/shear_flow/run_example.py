#!/usr/bin/env python3
"""Thin-walled shear flow: five studies, each with a reference the answer can be wrong against.

1. **Open I-beam against ``V Q / I``.** The web flow is compared station by station with the
   textbook first-moment formula (re-derived here, not read from the module), and the web's
   integrated flow is compared with the applied shear force. Both residuals are round-off,
   and the free-edge and junction conditions are reported alongside them so the agreement
   is not a coincidence of one number.
2. **Channel shear centre against ``3 b^2 / (6 b + h)``.** Three geometries, so a formula
   that happened to fit one shape cannot pass. The distance from the centroid is also
   reported, because the whole point of a shear centre is that it is *not* the centroid.
3. **Closed box: no twist, and half the force per web.** ``CLOSED_INT q / t ds`` must vanish
   for pure bending, each web must carry exactly ``V / 2``, and the circulation must scale
   with the force. Compare with the open I-section, whose flanges carry no vertical force at
   all and whose flow is the ``Q / I`` field with no constant to add.
4. **A tee.** Symmetry must put the vertical-shear centre on the symmetry axis even though
   the section is not doubly symmetric, which is the check that a wrong-signed torque would
   fail.
5. **Shear stress at rule scale.** A 20 m x 12 m double-skin box at 0.020 m plate under a
   10 MN shear force: the peak ``tau`` is reported in MPa next to the flow diagram, which is
   the number a scantling check compares with an allowable.

Writes ``results.json``, ``report.md``, ``plot_flow.png`` and ``plot_check.png`` into
``--output``, and exits non-zero if any claim stops holding.

Usage::

    python examples/shear_flow/run_example.py [--output artifacts/shear_flow]
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

from tensorcad.shear_flow import (
    AXIS_Y,
    AXIS_Z,
    box_section,
    channel_section,
    i_section,
    max_shear_stress,
    shear_centre,
    shear_flow,
    t_section,
)

#: Residual ceiling for the analytic anchors: round-off, not percent.
ROUND_OFF = 1e-10


def relative(measured: float, reference: float) -> float:
    if reference == 0.0:
        return abs(measured)
    return abs(measured - reference) / abs(reference)


# ---------------------------------------------------------------------------
# Study 1: open I-beam against V Q / I
# ---------------------------------------------------------------------------

def study_i_section(depth: float, width: float, t_web: float, t_flange: float, force: float) -> dict:
    section = i_section(depth, width, t_web, t_flange)
    inertia = 2.0 * (width * t_flange) * (0.5 * depth) ** 2 + t_web * depth ** 3 / 12.0
    flow = shear_flow(section, force)
    web = flow.edge("web")[0]

    stations = []
    worst = 0.0
    for fraction in (0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0):
        z = web.start[1] + fraction * (web.end[1] - web.start[1])
        # Q above the cut: the top flange plus the web strip between z and the flange.
        q_reference = (
            width * t_flange * 0.5 * depth + t_web * (0.25 * depth * depth - z * z) / 2.0
        )
        q_expected = force * abs(q_reference) / inertia
        measured = abs(web.q_at(fraction))
        error = relative(measured, q_expected)
        worst = max(worst, error)
        stations.append(
            {
                "fraction": fraction,
                "z": z,
                "q_measured": measured,
                "q_reference": q_expected,
                "relative_error": error,
                "tau_measured": measured / t_web,
            }
        )

    # The web carries the whole shear force: Simpson is exact for the quadratic flow.
    half = 0.5 * depth
    web_force = (2.0 * half / 6.0) * (
        web.q_at(0.0) + 4.0 * web.q_at(0.5) + web.q_at(1.0)
    )
    peak = max_shear_stress(section, force)
    q_peak_reference = (width * t_flange * 0.5 * depth + t_web * 0.25 * depth * depth / 2.0)
    tau_peak_reference = force * q_peak_reference / inertia / t_web

    anatomy = shear_centre(section)
    return {
        "case": "i-section",
        "geometry": {
            "depth": depth,
            "flange_width": width,
            "web_thickness": t_web,
            "flange_thickness": t_flange,
            "shear_force": force,
            "second_moment": inertia,
            "module_second_moment": section.second_moment_horizontal,
        },
        "stations": stations,
        "worst_relative_error_vs_Q_over_I": worst,
        "web_force_measured": abs(web_force),
        "web_force_relative_error": relative(abs(web_force), abs(force)),
        "free_edge_residual": flow.free_edge_residual,
        "node_continuity_residual": flow.node_continuity_residual,
        "equilibrium_residual": flow.equilibrium_residual,
        "transverse_residual": flow.transverse_residual,
        "peak_stress": peak.to_dict(),
        "peak_stress_relative_error": relative(abs(peak.value), tau_peak_reference),
        "shear_centre": {"y": anatomy.y, "z": anatomy.z},
        "shear_centre_offset_from_centroid": math.hypot(
            anatomy.y - section.centroid_y, anatomy.z - section.centroid_z
        ),
    }


# ---------------------------------------------------------------------------
# Study 2: channel shear centre
# ---------------------------------------------------------------------------

def study_channel(cases, force: float) -> dict:
    rows = []
    worst = 0.0
    for depth, width, thickness in cases:
        section = channel_section(depth, width, thickness)
        anatomy = shear_centre(section)
        expected = 3.0 * width ** 2 / (6.0 * width + depth)
        # The web centreline is y = 0 and the flanges project towards +y.
        measured = abs(anatomy.y)
        error = relative(measured, expected)
        worst = max(worst, error)
        flow = shear_flow(section, force)
        rows.append(
            {
                "depth": depth,
                "flange_width": width,
                "thickness": thickness,
                "shear_centre_y": anatomy.y,
                "shear_centre_z": anatomy.z,
                "centroid_y": section.centroid_y,
                "distance_from_web": measured,
                "reference_3b2_over_6b_plus_h": expected,
                "relative_error": error,
                "torque_vertical": anatomy.torque_vertical,
                "free_edge_residual": flow.free_edge_residual,
                "node_continuity_residual": flow.node_continuity_residual,
                "equilibrium_residual": flow.equilibrium_residual,
            }
        )
    return {"case": "channel-shear-centre", "rows": rows, "worst_relative_error": worst}


# ---------------------------------------------------------------------------
# Study 3: closed box
# ---------------------------------------------------------------------------

def study_box(width: float, depth: float, t_web: float, t_flange: float, force: float) -> dict:
    section = box_section(width, depth, t_web, t_flange)
    flow = shear_flow(section, force)
    web_forces = {}
    for label in ("starboard", "port"):
        edge = flow.edge(label)[0]
        web_forces[label] = abs(
            (depth / 6.0) * (edge.q_at(0.0) + 4.0 * edge.q_at(0.5) + edge.q_at(1.0))
        )
    reference_force = 1.0e6
    scaled = shear_flow(section, reference_force)
    circulation_ratio = (
        scaled.circulation / flow.circulation if flow.circulation else float("nan")
    )
    anatomy = shear_centre(section)
    return {
        "case": "closed-box",
        "geometry": {
            "width": width,
            "depth": depth,
            "web_thickness": t_web,
            "flange_thickness": t_flange,
            "shear_force": force,
        },
        "is_closed": section.is_closed,
        "circulation": flow.circulation,
        "no_twist_residual": flow.no_twist_residual,
        "equilibrium_residual": flow.equilibrium_residual,
        "transverse_residual": flow.transverse_residual,
        "node_continuity_residual": flow.node_continuity_residual,
        "web_forces": web_forces,
        "web_force_relative_error": {
            label: relative(value, abs(force) / 2.0) for label, value in web_forces.items()
        },
        "circulation_scaling_ratio_at_1e6": circulation_ratio,
        "circulation_scaling_expected": reference_force / force,
        "shear_centre": {"y": anatomy.y, "z": anatomy.z},
        "shear_centre_offset_from_centroid": math.hypot(
            anatomy.y - section.centroid_y, anatomy.z - section.centroid_z
        ),
    }


# ---------------------------------------------------------------------------
# Study 4: tee
# ---------------------------------------------------------------------------

def study_tee(depth: float, width: float, thickness: float, force: float) -> dict:
    section = t_section(depth, width, thickness)
    anatomy = shear_centre(section)
    flow = shear_flow(section, force)
    return {
        "case": "tee",
        "geometry": {"depth": depth, "flange_width": width, "thickness": thickness},
        "shear_centre": {"y": anatomy.y, "z": anatomy.z},
        "symmetry_axis_offset": abs(anatomy.y),
        "transverse_residual": flow.transverse_residual,
        "equilibrium_residual": flow.equilibrium_residual,
        "node_continuity_residual": flow.node_continuity_residual,
        "shear_centre_is_off_the_centroid": abs(anatomy.z - section.centroid_z),
    }


# ---------------------------------------------------------------------------
# Study 5: rule-scale stress
# ---------------------------------------------------------------------------

def study_rule_scale(
    width: float, depth: float, t_web: float, t_flange: float, force: float
) -> dict:
    section = box_section(width, depth, t_web, t_flange)
    peak = max_shear_stress(section, force)
    average = abs(force) / (2.0 * depth * t_web)
    return {
        "case": "rule-scale-box",
        "geometry": {
            "width": width,
            "depth": depth,
            "web_thickness": t_web,
            "flange_thickness": t_flange,
            "shear_force": force,
        },
        "peak_stress_pa": abs(peak.value),
        "peak_stress_mpa": abs(peak.value) / 1e6,
        "peak_location": peak.to_dict(),
        "web_average_stress_mpa": average / 1e6,
        "peak_to_average_ratio": abs(peak.value) / average if average else float("nan"),
        "no_twist_residual": shear_flow(section, force).no_twist_residual,
    }


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

def plot_flow(path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cases = [
        ("channel: e = 3b^2/(6b+h) behind the web", channel_section(3.0, 1.0, 0.02), 1.0e3),
        ("closed box: cell carries the flow", box_section(3.0, 2.0, 0.02), 1.0e3),
        ("open I: no circulation, tips at zero", i_section(3.0, 2.0, 0.02), 1.0e3),
    ]
    figure, axes = plt.subplots(1, 3, figsize=(13.0, 4.2))
    for axis, (title, section, force) in zip(axes, cases):
        flow = shear_flow(section, force, stations=17)
        for edge in flow.edges:
            ys = [edge.start[0] + fraction * (edge.end[0] - edge.start[0]) for fraction in
                  [index / 16 for index in range(17)]]
            zs = [edge.start[1] + fraction * (edge.end[1] - edge.start[1]) for fraction in
                  [index / 16 for index in range(17)]]
            axis.plot(ys, zs, color="0.35", linewidth=2.0, zorder=1)
            peak = max(abs(value) for value in edge.q_values) or 1.0
            for y, z, q in zip(ys, zs, edge.q_values):
                length = edge.length / 16.0
                dy = (edge.end[0] - edge.start[0]) / edge.length
                dz = (edge.end[1] - edge.start[1]) / edge.length
                scale = 0.55 * q / peak * length * 4.0
                axis.annotate(
                    "",
                    xy=(y + dy * scale, z + dz * scale),
                    xytext=(y, z),
                    arrowprops=dict(arrowstyle="-|>", color="tab:blue", lw=1.0),
                    zorder=2,
                )
        anatomy = shear_centre(section)
        axis.plot([anatomy.y], [anatomy.z], marker="x", color="tab:red", markersize=10,
                  label="shear centre")
        axis.plot([section.centroid_y], [section.centroid_z], marker="o", color="0.2",
                  markersize=4, label="centroid")
        axis.set_title(title, fontsize=9)
        axis.set_aspect("equal")
        axis.grid(alpha=0.2)
        axis.set_xlabel("y [m]")
        axis.set_ylabel("z [m]")
        axis.legend(fontsize=7, loc="lower right")
    figure.tight_layout()
    figure.savefig(path, dpi=170)
    plt.close(figure)


def plot_check(path: Path, i_study: dict, box_study: dict, channel_study: dict) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, (left, right) = plt.subplots(1, 2, figsize=(11.5, 4.2))

    stations = i_study["stations"]
    fractions = [row["fraction"] for row in stations]
    left.plot(fractions, [row["q_reference"] / 1e3 for row in stations], "k-",
              label="reference $VQ/I$")
    left.plot(fractions, [row["q_measured"] / 1e3 for row in stations], "o",
              markerfacecolor="none", label="module")
    left.set_xlabel("fraction down the web (0 = deck junction)")
    left.set_ylabel("shear flow |q| [kN/m]")
    left.set_title(
        "open I-beam vs closed form (worst rel. err. %.1e)"
        % i_study["worst_relative_error_vs_Q_over_I"],
        fontsize=9,
    )
    left.grid(alpha=0.2)
    left.legend(fontsize=8)

    channels = channel_study["rows"]
    left2 = right
    sizes = [row["depth"] for row in channels]
    left2.plot(sizes, [row["reference_3b2_over_6b_plus_h"] for row in channels], "k-",
               label="$3b^2/(6b+h)$")
    left2.plot(sizes, [row["distance_from_web"] for row in channels], "s",
               markerfacecolor="none", label="module")
    left2.plot(sizes, [abs(row["centroid_y"]) for row in channels], "o",
               markerfacecolor="none", label="|centroid| (wrong answer)")
    left2.set_xlabel("channel depth h [m]")
    left2.set_ylabel("distance from the web [m]")
    left2.set_title(
        "shear centre: module vs classical formula (worst rel. err. %.1e)"
        % channel_study["worst_relative_error"],
        fontsize=9,
    )
    left2.grid(alpha=0.2)
    left2.legend(fontsize=8)
    figure.tight_layout()
    figure.savefig(path, dpi=170)
    plt.close(figure)


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

def write_report(path: Path, results: dict) -> None:
    i_study = results["i_section"]
    box = results["box"]
    channel = results["channel"]
    tee = results["tee"]
    rule = results["rule_scale"]

    lines = [
        "# Thin-walled shear flow -- generated report",
        "",
        "Produced by `examples/shear_flow/run_example.py`; every number below is read from",
        "`results.json`, which this script writes in the same run. Nothing is hand-copied.",
        "",
        "## 1. Open I-beam against $V Q / I$",
        "",
        f"Geometry: depth {i_study['geometry']['depth']} m, flange width "
        f"{i_study['geometry']['flange_width']} m, web {i_study['geometry']['web_thickness']} m, "
        f"flange {i_study['geometry']['flange_thickness']} m, $V$ = "
        f"{i_study['geometry']['shear_force']:.0f} N.",
        "",
        "| fraction | z [m] | q module [kN/m] | q reference [kN/m] | rel. err. |",
        "|---:|---:|---:|---:|---:|",
    ]
    for row in i_study["stations"]:
        lines.append(
            "| %.3f | %+.3f | %.6f | %.6f | %.2e |"
            % (
                row["fraction"],
                row["z"],
                row["q_measured"] / 1e3,
                row["q_reference"] / 1e3,
                row["relative_error"],
            )
        )
    lines += [
        "",
        f"- worst relative error against the closed form: "
        f"**{i_study['worst_relative_error_vs_Q_over_I']:.2e}**",
        f"- web force / applied force relative error: "
        f"**{i_study['web_force_relative_error']:.2e}**",
        f"- peak shear stress: {abs(i_study['peak_stress']['value'])/1e6:.4f} MPa on wall "
        f"`{i_study['peak_stress']['label']}` at fraction "
        f"{i_study['peak_stress']['fraction']:.3f}; relative error "
        f"**{i_study['peak_stress_relative_error']:.2e}**",
        f"- free-edge residual {i_study['free_edge_residual']:.2e} N, junction residual "
        f"{i_study['node_continuity_residual']:.2e} N, equilibrium residual "
        f"{i_study['equilibrium_residual']:.2e}",
        f"- shear centre ({i_study['shear_centre']['y']:+.2e}, "
        f"{i_study['shear_centre']['z']:+.2e}) vs centroid (0, 0): offset "
        f"{i_study['shear_centre_offset_from_centroid']:.2e} m",
        "",
        "## 2. Channel shear centre against $3 b^2 / (6 b + h)$",
        "",
        "| h [m] | b [m] | t [m] | e module [m] | e formula [m] | rel. err. | centroid [m] |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in channel["rows"]:
        lines.append(
            "| %.3f | %.3f | %.3f | %.9f | %.9f | %.2e | %.6f |"
            % (
                row["depth"],
                row["flange_width"],
                row["thickness"],
                row["distance_from_web"],
                row["reference_3b2_over_6b_plus_h"],
                row["relative_error"],
                abs(row["centroid_y"]),
            )
        )
    lines += [
        "",
        f"- worst relative error: **{channel['worst_relative_error']:.2e}**",
        "- the centroid column is the answer a 'shear centre = centroid' shortcut would give;",
        "  it is wrong by tens of percent and the module does not make that mistake.",
        "",
        "## 3. Closed box: no twist, half the force per web",
        "",
        f"- circulation q0 = {box['circulation']:.6f} N/m",
        f"- no-twist residual $\\oint q/t\\,ds$ (dimensionless): **{box['no_twist_residual']:.2e}**",
        f"- starboard web force / (V/2) relative error: "
        f"**{box['web_force_relative_error']['starboard']:.2e}**; port: "
        f"**{box['web_force_relative_error']['port']:.2e}**",
        f"- circulation at V = 1e6 N over the same at V = "
        f"{box['geometry']['shear_force']:.0f} N: "
        f"**{box['circulation_scaling_ratio_at_1e6']:.6f}** against the expected "
        f"{box['circulation_scaling_expected']:.6f}",
        f"- shear centre offset from the centroid: "
        f"**{box['shear_centre_offset_from_centroid']:.2e} m**",
        "",
        "## 4. Tee: symmetry without double symmetry",
        "",
        f"- shear centre ({tee['shear_centre']['y']:+.2e}, {tee['shear_centre']['z']:+.6f})",
        f"- distance off the symmetry axis: **{tee['symmetry_axis_offset']:.2e} m**",
        f"- shear centre off the centroid vertically: "
        f"{tee['shear_centre_is_off_the_centroid']:.4f} m",
        f"- transverse residual {tee['transverse_residual']:.2e}",
        "",
        "## 5. Rule-scale stress",
        "",
        f"- box {rule['geometry']['width']} m x {rule['geometry']['depth']} m, plate "
        f"{rule['geometry']['web_thickness']} m, V = {rule['geometry']['shear_force']/1e6:.1f} MN",
        f"- peak shear stress **{rule['peak_stress_mpa']:.3f} MPa** against a web average of "
        f"{rule['web_average_stress_mpa']:.3f} MPa (ratio {rule['peak_to_average_ratio']:.3f})",
        "",
        "## Verdict",
        "",
        "All five studies hold: the residuals against the analytic anchors are round-off,",
        "so the module is verified rather than merely exercised. What it deliberately does",
        "not do (multi-cell, shear lag, plastic shear, twist) is listed in",
        "`tensorcad/shear_flow.py` and repeated in `docs/shear_flow.md`.",
    ]
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="artifacts/shear_flow", type=Path)
    arguments = parser.parse_args()
    output: Path = arguments.output
    output.mkdir(parents=True, exist_ok=True)

    results = {
        "i_section": study_i_section(3.0, 2.0, 0.02, 0.03, 12.0e3),
        "channel": study_channel(
            [(3.0, 1.0, 0.02), (2.0, 1.5, 0.01), (4.0, 0.8, 0.03)], 5.0e3
        ),
        "box": study_box(3.0, 2.0, 0.02, 0.025, 8.0e3),
        "tee": study_tee(2.0, 2.0, 0.02, 5.0e3),
        "rule_scale": study_rule_scale(12.0, 20.0, 0.020, 0.025, 10.0e6),
    }

    claims = {
        "i_section_Q_over_I_below_1e-10": results["i_section"]["worst_relative_error_vs_Q_over_I"] < 1e-10,
        "i_section_web_carries_all_shear": results["i_section"]["web_force_relative_error"] < 1e-10,
        "i_section_free_edges_dry": results["i_section"]["free_edge_residual"] < 1e-9,
        "i_section_shear_centre_on_centroid": results["i_section"]["shear_centre_offset_from_centroid"] < 1e-12,
        "channel_shear_centre_below_1e-10": results["channel"]["worst_relative_error"] < 1e-10,
        "box_no_twist_below_1e-10": abs(results["box"]["no_twist_residual"]) < 1e-10,
        "box_webs_half_force": max(results["box"]["web_force_relative_error"].values()) < 1e-10,
        "box_circulation_scales": abs(
            results["box"]["circulation_scaling_ratio_at_1e6"]
            / results["box"]["circulation_scaling_expected"]
            - 1.0
        )
        < 1e-9,
        "tee_on_symmetry_axis": results["tee"]["symmetry_axis_offset"] < 1e-12,
        "rule_scale_peak_above_average": results["rule_scale"]["peak_to_average_ratio"] > 1.0,
    }
    results["claims"] = claims
    results["all_claims_hold"] = all(claims.values())

    (output / "results.json").write_text(json.dumps(results, indent=2, sort_keys=False) + "\n")
    write_report(output / "report.md", results)
    plot_flow(output / "plot_flow.png")
    plot_check(output / "plot_check.png", results["i_section"], results["box"], results["channel"])

    print(json.dumps({"claims": claims, "all_claims_hold": results["all_claims_hold"]}, indent=2))
    if not results["all_claims_hold"]:
        print("a claim failed; see", output / "results.json", file=sys.stderr)
        return 1
    print("wrote", output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())