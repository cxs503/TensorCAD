"""Generate the elliptic-section evidence bundle: numbers, figure, document.

Writes ``artifacts/elliptic_section/{results.json,plot.png,report.md}`` and exits
non-zero if any claim made in the document stops holding. Every number in
``docs/elliptic_section.md`` is read from ``results.json`` -- none is typed by hand.

Run: ``python examples/elliptic_section/run_example.py``
"""

from __future__ import annotations

import json
import math
import pathlib
import sys

import tensorcad
from tensorcad.elliptic_section import (
    elliptic_gm0,
    elliptic_prism_mesh,
    elliptic_righting_arm,
    elliptic_righting_arm_curve,
    wall_sided_gm0,
)
from tensorcad.righting_arm import righting_arm_curve

OUT = pathlib.Path(__file__).resolve().parents[2] / "artifacts" / "elliptic_section"

SECTIONS = ((1.0, 1.0), (1.0, 0.6), (0.5, 1.0), (1.0, 0.35))
HEELS = (0.0, 5.0, 10.0, 20.0, 30.0, 40.0, 60.0, 80.0)
KG = 0.4
LENGTH = 4.0


def _signed_volume(mesh) -> float:
    total = 0.0
    for face in mesh.faces:
        a, b, c = (mesh.vertices[index] for index in face)
        total += (
            a[0] * (b[1] * c[2] - b[2] * c[1])
            - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0])
        ) / 6.0
    return total


def build_results() -> dict:
    cases = []
    for a, b in SECTIONS:
        mesh = elliptic_prism_mesh(a, b, LENGTH, ring=384)
        curve = righting_arm_curve(mesh, kg=KG, heels_deg=HEELS, relative_tolerance=1e-11)
        points = curve.to_dict()["points"]
        rows = []
        for point in points:
            analytic = elliptic_righting_arm(a, b, point["heel_deg"], kg=KG)
            rows.append(
                {
                    "heel_deg": point["heel_deg"],
                    "gz_mesh": point["righting_arm"],
                    "gz_analytic": analytic,
                    "abs_error": abs(point["righting_arm"] - analytic),
                }
            )
        worst = max(row["abs_error"] for row in rows if abs(row["gz_analytic"]) > 1e-3)
        scale = max(abs(row["gz_analytic"]) for row in rows)
        cases.append(
            {
                "a": a,
                "b": b,
                "ring": 384,
                "rows": rows,
                "worst_abs_error_m": worst,
                "worst_rel_error": worst / scale,
                "gm0_exact": elliptic_gm0(a, b, kg=KG),
                "gm0_wall_sided": wall_sided_gm0(a, b, kg=KG),
                "record": elliptic_righting_arm_curve(a, b, kg=KG),
            }
        )

    convergence = []
    for ring in (24, 48, 96, 192, 384):
        mesh = elliptic_prism_mesh(1.0, 0.6, LENGTH, ring=ring)
        curve = righting_arm_curve(mesh, kg=KG, heels_deg=HEELS, relative_tolerance=1e-11)
        error = max(
            abs(point["righting_arm"] - elliptic_righting_arm(1.0, 0.6, point["heel_deg"], kg=KG))
            for point in curve.to_dict()["points"]
        )
        convergence.append({"ring": ring, "worst_abs_error_m": error})

    circle = {}
    for radius in (0.25, 0.5, 1.0):
        exact = elliptic_gm0(radius, radius, kg=KG)
        wall = wall_sided_gm0(radius, radius, kg=KG)
        circle[str(radius)] = {
            "exact_gm0": exact,
            "wall_sided_gm0": wall,
            "exact_minus_wall": exact - wall,
            "relative_disagreement": abs(wall - exact) / abs(exact),
        }

    volume_check = []
    for a, b in SECTIONS:
        mesh = elliptic_prism_mesh(a, b, LENGTH, ring=384)
        volume_check.append(
            {
                "a": a,
                "b": b,
                "mesh_volume": _signed_volume(mesh),
                "analytic_volume": math.pi * a * b * LENGTH,
                "rel_error": abs(_signed_volume(mesh) - math.pi * a * b * LENGTH)
                / (math.pi * a * b * LENGTH),
            }
        )

    return {
        "schema": "tensorcad.elliptic_section.evidence/1",
        "tensorcad_version": getattr(tensorcad, "__version__", "unknown"),
        "kg": KG,
        "length": LENGTH,
        "heels_deg": HEELS,
        "cases": cases,
        "ring_convergence": convergence,
        "circle_limit": circle,
        "volume_check": volume_check,
    }


def checks(results: dict) -> list[tuple[str, bool, str]]:
    out = []
    for case in results["cases"]:
        out.append(
            (
                f"mesh GZ reproduces the closed form, a={case['a']} b={case['b']}",
                case["worst_rel_error"] < 3e-4,
                f"worst relative error {case['worst_rel_error']:.3e}",
            )
        )
    rings = [row["ring"] for row in results["ring_convergence"]]
    errors = [row["worst_abs_error_m"] for row in results["ring_convergence"]]
    out.append(
        (
            "discretisation error falls with ring count",
            all(earlier > later for earlier, later in zip(errors, errors[1:])),
            f"rings {rings} -> {['%.2e' % value for value in errors]}",
        )
    )
    for row in results["volume_check"]:
        out.append(
            (
                f"mesh volume == pi a b L, a={row['a']} b={row['b']}",
                row["rel_error"] < 2e-4,
                f"relative error {row['rel_error']:.2e}",
            )
        )
    for radius, row in results["circle_limit"].items():
        out.append(
            (
                f"circle limit R={radius}: exact GM0 = R - KG",
                abs(row["exact_gm0"] - (float(radius) - KG)) < 1e-12,
                f"{row['exact_gm0']:.9f} vs {float(radius) - KG:.9f}",
            )
        )
    return out


def make_figure(results: dict, path: pathlib.Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.4))
    for case in results["cases"]:
        label = f"a={case['a']:g}, b={case['b']:g}"
        heels = [row["heel_deg"] for row in case["rows"]]
        axes[0].plot(heels, [row["gz_analytic"] for row in case["rows"]], "-", label=f"{label} exact")
        axes[0].plot(
            heels,
            [row["gz_mesh"] for row in case["rows"]],
            "o",
            markersize=3.2,
            label=f"{label} mesh",
        )
    axes[0].set_xlabel("heel angle (deg)")
    axes[0].set_ylabel("GZ (m)")
    axes[0].set_title(f"Analytic vs mesh righting arm (KG = {KG} m)")
    axes[0].grid(alpha=0.3)
    axes[0].legend(fontsize=7, ncol=2)

    rings = [row["ring"] for row in results["ring_convergence"]]
    errors = [row["worst_abs_error_m"] for row in results["ring_convergence"]]
    axes[1].loglog(rings, errors, "s-", label="max |GZ_mesh - GZ_analytic| (a=1, b=0.6)")
    reference = [errors[0] * (rings[0] / ring) ** 2 for ring in rings]
    axes[1].loglog(rings, reference, "k--", label="O(h^2) reference")
    axes[1].set_xlabel("ellipse ring points")
    axes[1].set_ylabel("worst absolute error (m)")
    axes[1].set_title("Discretisation error of the mesh integrator")
    axes[1].grid(alpha=0.3, which="both")
    axes[1].legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def make_report(results: dict, verdicts: list[tuple[str, bool, str]]) -> str:
    lines = []
    lines.append("# Elliptic-section righting arm: analytic family vs mesh integrator\n")
    lines.append(
        "Generated by `python examples/elliptic_section/run_example.py`. Every number below is "
        "read from `artifacts/elliptic_section/results.json`.\n"
    )
    lines.append("## 1. Closed form\n")
    lines.append(
        "For a prism of elliptic cross-section (semi-axes `a` transverse, `b` vertical) floating at\n"
        "half immersion, the immersed section at every heel is an exact half-ellipse, so\n\n"
        "    KN(phi) = b sin(phi) + (4/(3 pi)) (a^2 - b^2) sin(phi) cos(phi) / sqrt(a^2 sin^2 + b^2 cos^2)\n\n"
        "    GZ(phi) = KN(phi) - KG sin(phi)\n\n"
        f"with `KG = {KG} m` from the keel and `L = {results['length']} m`.\n"
    )
    lines.append("## 2. Analytic GZ vs the mesh integrator\n")
    lines.append("| a | b | GM0 exact (m) | GM0 wall-sided (m) | worst |GZ_mesh - GZ_exact| | worst relative |")
    lines.append("|---:|---:|---:|---:|---:|---:|")
    for case in results["cases"]:
        lines.append(
            f"| {case['a']:g} | {case['b']:g} | {case['gm0_exact']:.6f} | "
            f"{case['gm0_wall_sided']:.6f} | {case['worst_abs_error_m']:.3e} | "
            f"{case['worst_rel_error']:.3e} |"
        )
    lines.append("")
    lines.append("### Full ladders\n")
    for case in results["cases"]:
        lines.append(f"**a = {case['a']:g}, b = {case['b']:g}**\n")
        lines.append("| heel (deg) | GZ mesh (m) | GZ exact (m) | |error| (m) |")
        lines.append("|---:|---:|---:|---:|")
        for row in case["rows"]:
            lines.append(
                f"| {row['heel_deg']:g} | {row['gz_mesh']:.9f} | {row['gz_analytic']:.9f} | "
                f"{row['abs_error']:.3e} |"
            )
        lines.append("")
    lines.append("## 3. Discretisation error under ring refinement (a = 1, b = 0.6)\n")
    lines.append("| ring points | worst absolute error (m) |")
    lines.append("|---:|---:|")
    for row in results["ring_convergence"]:
        lines.append(f"| {row['ring']} | {row['worst_abs_error_m']:.3e} |")
    lines.append("")
    lines.append("## 4. The circle limit, and why `KB + BM - KG` is not the answer\n")
    lines.append("| R | GM0 exact (m) | GM0 wall-sided (m) | exact - wall (m) | relative disagreement |")
    lines.append("|---:|---:|---:|---:|---:|")
    for radius, row in results["circle_limit"].items():
        lines.append(
            f"| {radius} | {row['exact_gm0']:.9f} | {row['wall_sided_gm0']:.9f} | "
            f"{row['exact_minus_wall']:.9f} | {row['relative_disagreement']:.2%} |"
        )
    lines.append("")
    lines.append("## 5. Gates\n")
    for label, ok, detail in verdicts:
        lines.append(f"- `{'PASS' if ok else 'FAIL'}` {label} -- {detail}")
    lines.append("")
    return "\n".join(lines)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    results = build_results()
    verdicts = checks(results)
    (OUT / "results.json").write_text(json.dumps(results, indent=2, sort_keys=True))
    make_figure(results, OUT / "plot.png")
    report = make_report(results, verdicts)
    (OUT / "report.md").write_text(report)
    (pathlib.Path(__file__).resolve().parents[2] / "docs" / "elliptic_section.md").write_text(report)
    failed = [label for label, ok, _ in verdicts if not ok]
    for label, ok, detail in verdicts:
        print(f"{'PASS' if ok else 'FAIL'}  {label}: {detail}")
    if failed:
        print(f"\n{len(failed)} gate(s) failed", file=sys.stderr)
        return 1
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())