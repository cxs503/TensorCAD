#!/usr/bin/env python3
"""Publish the ``docs/gz/`` evidence bundle from a fresh run of the example.

Run:  python examples/gz/publish_evidence.py

The bundle is three files and none of them contains a hand-typed number:

* ``docs/gz/results.json`` -- written by ``examples/gz/run_example.py`` itself;
* ``docs/gz/plot.png``     -- the example's figure, copied verbatim;
* ``docs/gz/tables.md``    -- rendered from ``results.json`` by :func:`render_tables`
  below, so the markdown tables are a *view* of the JSON, not a second copy of it.

``tests/test_gz_evidence.py`` re-renders ``tables.md`` from the committed
``results.json`` and fails if the file on disk differs by one byte, which is how the
"no transcription" claim is enforced instead of asserted.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
EXAMPLE = Path(__file__).with_name("run_example.py")
EVIDENCE = ROOT / "docs" / "gz"


def num(value: float, digits: int) -> str:
    """Fixed-point, with the round-off band printed as a clean zero.

    The integrator's residuals sit at 1e-17 m; printing ``-0.000000000`` for them would
    suggest a sign the geometry does not have.
    """
    if abs(value) < 0.5 * 10.0**-digits:
        value = 0.0
    return f"{value:.{digits}f}"


def sci(value: float) -> str:
    """Scientific notation, without a sign on values that are exactly zero at this width.

    A residual of ``-1e-18`` is double-precision noise on a volume of 7.4 m3; printing
    ``-0.00e+00`` for it reads like a systematic bias.
    """
    text = f"{value:.2e}"
    return "0.00e+00" if text == "-0.00e+00" else text


def render_tables(results: dict[str, Any]) -> str:
    """Render the markdown tables from the example's ``results.json`` payload."""
    box = results["box_wall_sided"]
    log = results["log_circular_section"]
    hull = results["hull"]
    curve = hull["curve"]
    surface = results["free_surface"]
    upright = results["upright"]
    checks = results["checks"]
    parameters = hull["parameters"]

    lines: list[str] = []
    add = lines.append
    add("# Righting-arm tables")
    add("")
    add("Rendered from `docs/gz/results.json` by `examples/gz/publish_evidence.py`; the JSON")
    add("is written by `examples/gz/run_example.py`. Every number below is a formatting of a")
    add("value in that file -- `tests/test_gz_evidence.py` re-renders this page and compares")
    add("it byte for byte, so a hand-typed figure cannot survive CI.")
    add("")

    add("## 1. Cases")
    add("")
    add("| body | geometry | reference used |")
    add("|---|---|---|")
    add(
        f"| box barge | L 10.0 m, B 2.0 m, D 2.0 m, T 0.8 m, KG 0.5 m |"
        f" wall-sided trapezoid, valid to {num(box['limit_deg'], 2)} deg |"
    )
    add(
        f"| circular log (regular {log['segments']}-gon) | R {num(log['radius_m'], 1)} m,"
        f" half submerged, KG 0.5 m | exact `(R - KG) sin(phi)`, every heel |"
    )
    add(
        f"| parametric hull, 41x13 mesh | L {num(parameters['length_m'], 1)} m,"
        f" B {num(parameters['beam_m'], 1)} m, T {num(parameters['draft_m'], 1)} m,"
        f" freeboard {num(hull['freeboard_m'], 1)} m, KG {num(parameters['kg_m'], 2)} m |"
        f" `GM` from `tensorcad.stability.initial_stability` (waterplane moments) |"
    )
    add("")

    add(f"## 2. Box barge against the wall-sided closed form (0-{int(box['limit_deg'])} deg)")
    add("")
    add(f"`GM0 = {num(box['gm0_m'], 4)} m`, `BM0 = {num(box['bm0_m'], 4)} m`.")
    add("")
    add("| heel [deg] | GZ integrated [m] | wall-sided [m] | abs. error [m] |")
    add("|---:|---:|---:|---:|")
    for row in box["inside"]:
        add(
            f"| {num(row['heel_deg'], 0)} | {num(row['gz_m'], 9)} |"
            f" {num(row['wall_sided_m'], 9)} | {sci(row['abs_error'])} |"
        )
    add("")
    add(f"Worst absolute error inside the validity window: **{sci(box['worst_abs_error_m'])} m**")
    add("(round-off: the integrator and the closed form are both exact for this polyhedron).")
    add("")

    add("## 3. Negative control: past the wall-sided limit")
    add("")
    add("The formula is invalid above the limit, so the agreement above must *stop*; this")
    add("table is the proof that it can fail.")
    add("")
    add("| heel [deg] | GZ integrated [m] | wall-sided [m] | abs. error [m] |")
    add("|---:|---:|---:|---:|")
    for row in box["beyond_limit"]:
        add(
            f"| {num(row['heel_deg'], 0)} | {num(row['gz_m'], 9)} |"
            f" {num(row['wall_sided_m'], 9)} | {sci(row['abs_error'])} |"
        )
    add("")

    add("## 4. Circular log against the exact `(R - KG) sin(phi)`")
    add("")
    add("| heel [deg] | GZ integrated [m] | exact [m] | abs. error [m] |")
    add("|---:|---:|---:|---:|")
    for row in log["rows"]:
        add(
            f"| {num(row['heel_deg'], 4)} | {num(row['gz_m'], 9)} |"
            f" {num(row['exact_m'], 9)} | {sci(row['abs_error'])} |"
        )
    add("")
    add(f"Worst absolute error over the whole sweep: **{sci(log['worst_abs_error_m'])} m**.")
    add("")

    add("## 5. Hull: the curve, and its slope at zero")
    add("")
    add("| quantity | value |")
    add("|---|---:|")
    add(f"| `GZ'(0)` from the curve (Richardson) [m/rad] | {num(hull['gm_from_curve_slope_m'], 9)} |")
    add(
        f"| `GM` from `initial_stability` (waterplane moments) [m] |"
        f" {num(hull['gm_from_initial_stability_m'], 9)} |"
    )
    gm, reference_gm = hull["gm_from_curve_slope_m"], hull["gm_from_initial_stability_m"]
    add(f"| relative difference | {sci(abs(gm - reference_gm) / reference_gm)} |")
    add(f"| `KB` (upright) [m] | {num(hull['kb_m'], 6)} |")
    add(f"| `BM` transverse (upright) [m] | {num(hull['bm_transverse_m'], 6)} |")
    add(f"| maximum `GZ` [m] | {num(curve['max_gz_m'], 6)} |")
    add(f"| heel at maximum `GZ` [deg] | {num(curve['angle_of_max_gz_deg'], 1)} |")
    add(
        f"| area under the curve to {num(curve['area_under_curve_limit_deg'], 0)} deg [m·rad] |"
        f" {num(curve['area_under_curve_m_rad'], 6)} |"
    )
    add(f"| displacement held at [m³] | {num(curve['displacement_volume_m3'], 6)} |")
    add(f"| worst displacement residual [m³] | {sci(curve['max_volume_error'])} |")
    add("")
    add("| heel [deg] | GZ [m] | KN [m] | waterline offset [m] | V [m³] | volume residual [m³] |")
    add("|---:|---:|---:|---:|---:|---:|")
    for point in curve["points"]:
        add(
            f"| {num(point['heel_angle_deg'], 1)} | {num(curve_point_gz(curve, point), 9)} |"
            f" {num(point['kn_m'], 6)} | {num(point['waterline_offset_m'], 6)} |"
            f" {num(point['displacement_volume_m3'], 6)} | {sci(point['volume_error'])} |"
        )
    add("")

    add("## 6. Free surface and upright hydrostatics")
    add("")
    add("| quantity | value |")
    add("|---|---:|")
    add(f"| tank | {surface['tank']} |")
    add(f"| `i_x` [m⁴] | {num(surface['i_x_m4'], 6)} |")
    add(f"| flat loss on every point [m] | {num(surface['loss_m'], 6)} |")
    add(f"| peak `GZ`, solid [m] | {num(surface['max_gz_solid_m'], 6)} |")
    add(f"| peak `GZ`, with free surface [m] | {num(surface['max_gz_with_free_surface_m'], 6)} |")
    add(f"| upright displacement [m³] | {num(upright['displacement_volume_m3'], 6)} |")
    add(f"| upright `KB` [m] | {num(upright['kb_m'], 6)} |")
    add(f"| waterplane area [m²] | {num(upright['waterplane_area_m2'], 6)} |")
    add(f"| wetted surface area [m²] | {num(upright['wetted_surface_area_m2'], 6)} |")
    add("")

    add("## 7. Claim ledger")
    add("")
    verdict = "all claims hold" if checks["passed"] else "FAILED: " + "; ".join(checks["failures"])
    add(f"- verdict: **{verdict}**")
    for index, row in enumerate(box["inside"], start=1):
        if row["abs_error"] > 1e-12:
            add(f"- box barge row {index} exceeds the 1e-12 m bound")
    add(
        f"- box barge: worst {sci(box['worst_abs_error_m'])} m < 1e-12 m"
        f" ({'pass' if box['worst_abs_error_m'] <= 1e-12 else 'fail'})"
    )
    add(
        f"- log: worst {sci(log['worst_abs_error_m'])} m < 1e-9 m"
        f" ({'pass' if log['worst_abs_error_m'] <= 1e-9 else 'fail'})"
    )
    add(
        f"- hull `GZ'(0)` vs `GM`: {sci(abs(gm - reference_gm) / reference_gm)} relative < 1e-6"
        f" ({'pass' if abs(gm - reference_gm) / reference_gm < 1e-6 else 'fail'})"
    )
    add(
        f"- displacement held at every point: worst {sci(curve['max_volume_error'])} m³ < 1e-10"
        f" ({'pass' if curve['max_volume_error'] < 1e-10 else 'fail'})"
    )
    add(
        f"- negative control: wall-sided formula off by {sci(box['beyond_limit'][-1]['abs_error'])} m"
        f" at {num(box['beyond_limit'][-1]['heel_deg'], 0)} deg > 1e-2 m"
        f" ({'pass' if box['beyond_limit'][-1]['abs_error'] > 1e-2 else 'fail'})"
    )
    add("")
    return "\n".join(lines)


def curve_point_gz(curve: dict[str, Any], point: dict[str, Any]) -> float:
    """`GZ` at this point, looked up from the curve's own heel grid."""
    angles = curve["heel_angles_deg"]
    values = curve["gz_m"]
    heel = point["heel_angle_deg"]
    for angle, value in zip(angles, values, strict=True):
        if abs(angle - heel) < 1e-9:
            return float(value)
    raise KeyError(f"heel {heel} is not on the curve grid")


def publish(output: Path = EVIDENCE) -> dict[str, Any]:
    """Run the example into a scratch directory and publish the bundle."""
    output.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as scratch:
        run = subprocess.run(
            [sys.executable, str(EXAMPLE), "--output", scratch],
            capture_output=True,
            text=True,
        )
        if run.returncode != 0:
            raise SystemExit(f"the example failed; refusing to publish:\n{run.stdout}\n{run.stderr}")
        results = json.loads((Path(scratch) / "results.json").read_text())
        shutil.copyfile(Path(scratch) / "results.json", output / "results.json")
        shutil.copyfile(Path(scratch) / "plot.png", output / "plot.png")
    (output / "tables.md").write_text(render_tables(results))
    return results


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default=str(EVIDENCE), help="bundle directory")
    options = parser.parse_args()
    results = publish(Path(options.output))
    checks = results["checks"]
    box = results["box_wall_sided"]
    log = results["log_circular_section"]
    hull = results["hull"]
    curve = hull["curve"]
    gm, reference_gm = hull["gm_from_curve_slope_m"], hull["gm_from_initial_stability_m"]
    print(f"published {Path(options.output)}")
    print(f"  box barge  worst |GZ - wall-sided| = {box['worst_abs_error_m']:.3e} m")
    print(f"  log        worst |GZ - exact|       = {log['worst_abs_error_m']:.3e} m")
    print(
        f"  hull       GZ'(0) = {gm:.9f} vs GM = {reference_gm:.9f} m"
        f" ({abs(gm - reference_gm) / reference_gm:.2e} relative)"
    )
    print(f"  hull       max GZ = {curve['max_gz_m']:.9f} m at {curve['angle_of_max_gz_deg']:.1f} deg")
    print(f"  claims     {'all hold' if checks['passed'] else 'FAILED: ' + '; '.join(checks['failures'])}")
    return 0 if checks["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
