"""Audit driver for the prism/pyramid metric set: four meshes, three independent checks.

Run from the repository root::

    python tools/prism_pyramid_audit.py

Writes ``artifacts/mesh_quality_prism_pyramid/audit.json`` (+ ``quality.png``,
``manifest.json``) and exits 1 if any criterion fails, so it is usable as a CI step as
well as a record. The three checks per mesh are:

1. **Closed form** -- the summed element volume against the exact geometric volume of the
   shape (a 1x2x2 plate prism stack, a unit cube filled by six pyramids), and, for the
   ideal elements, the analytic metric values.
2. **Second route** -- every metric against
   :func:`tensorcad.meshquality.reference.reference_prisms` /
   ``reference_pyramids`` (finite-difference Jacobian, Cayley-Menger volumes, ``atan2``
   angles). Tolerance 1e-6 (the mission acceptance); the measured maximum is reported.
3. **Gate** -- the shipped ``DEFAULT_GATE`` verdict, with the failing metric and element.
"""
from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tensorcad.meshquality import DEFAULT_GATE, measure  # noqa: E402
from tensorcad.meshquality.reference import reference_prisms, reference_pyramids  # noqa: E402

ARTIFACTS = ROOT / "artifacts" / "mesh_quality_prism_pyramid"
TOLERANCE = 1.0e-6


# --------------------------------------------------------------------------------------
# meshes
# --------------------------------------------------------------------------------------
def prism_plate(nx: int = 2, ny: int = 2, layers=(0.05, 0.10, 0.20)) -> dict:
    """A prism boundary layer of three geometrically growing layers over an nx by ny patch.

    Each column of the patch is split into two prisms along the diagonal, which is the
    arrangement a hybrid mesher produces above a quad surface mesh.
    """
    points: list[tuple[float, float, float]] = []
    prisms: list[tuple[int, ...]] = []
    zs = [0.0]
    for height in layers:
        zs.append(zs[-1] + height)
    for layer in range(len(layers)):
        z0, z1 = zs[layer], zs[layer + 1]
        for i in range(nx):
            for j in range(ny):
                x0, x1 = float(i), float(i + 1)
                y0, y1 = float(j), float(j + 1)
                base = len(points)
                points.extend(
                    [
                        (x0, y0, z0), (x1, y0, z0), (x0, y1, z0),
                        (x0, y0, z1), (x1, y0, z1), (x0, y1, z1),
                        (x1, y1, z0), (x1, y1, z1),
                    ]
                )
                prisms.append((base + 0, base + 1, base + 2, base + 3, base + 4, base + 5))
                prisms.append((base + 1, base + 6, base + 2, base + 4, base + 7, base + 5))
    return {
        "name": "prism_plate_layer",
        "kind": "prism",
        "points": points,
        "elements": prisms,
        "exact_volume": float(nx * ny * zs[-1]),
        "note": "3 growing prism layers, 2x2 columns, conforming diagonal split",
        "expected_gate": "PASS",
    }


def pyramid_cube() -> dict:
    """A unit cube filled by six pyramids with apex at the centre (the tet-hex transition)."""
    corners = [
        (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0), (0.0, 1.0, 1.0),
    ]
    points = corners + [(0.5, 0.5, 0.5)]
    apex = 8
    # The base quad must wind counter-clockwise seen from the apex (the right-hand
    # normal points *into* the pyramid), which is the reverse of the outward winding a
    # hexahedron face uses. Written the hex way the six cells come out at
    # scaled_jacobian = -1.0, which is exactly what ``inverted_pyramid_fan`` records.
    faces = (
        (0, 1, 2, 3),  # z = 0
        (4, 7, 6, 5),  # z = 1
        (0, 4, 5, 1),  # y = 0
        (1, 5, 6, 2),  # x = 1
        (2, 6, 7, 3),  # y = 1
        (3, 7, 4, 0),  # x = 0
    )
    pyramids = [tuple(face) + (apex,) for face in faces]
    return {
        "name": "pyramid_cube_fan",
        "kind": "pyramid",
        "points": points,
        "elements": pyramids,
        "exact_volume": 1.0,
        "note": "six pyramids, apex at the cube centre (hex-to-tet transition cell)",
        "expected_gate": "PASS",
    }


def sheared_prisms() -> dict:
    """Non-affine prisms: the top face is slid, so the Jacobian varies inside the cell."""
    points: list[tuple[float, float, float]] = []
    prisms: list[tuple[int, ...]] = []
    for layer, shear in enumerate((0.0, 0.4, 0.9)):
        z0, z1 = float(layer), float(layer + 1.0)
        base = len(points)
        points.extend(
            [
                (0.0, 0.0, z0), (1.0, 0.0, z0), (0.0, 1.0, z0),
                (shear, 0.0, z1), (1.0 + shear, 0.0, z1), (shear, 1.0, z1),
            ]
        )
        prisms.append((base, base + 1, base + 2, base + 3, base + 4, base + 5))
    return {
        "name": "sheared_prism_stack",
        "kind": "prism",
        "points": points,
        "elements": prisms,
        "exact_volume": 1.5,
        "note": "three prisms, top face sheared by 0, 0.4, 0.9 (non-affine map)",
        "expected_gate": "PASS",
    }


def collapsed_prism() -> dict:
    """A prism layer whose last cell is collapsed onto an edge: the failure the gate is for."""
    mesh = prism_plate()
    points = list(mesh["points"])
    points[3] = (0.0, 0.0, 0.05)
    points[4] = (1.0, 0.0, 0.05)
    points[5] = (0.5, 0.0, 0.05)
    mesh["points"] = points
    mesh["name"] = "collapsed_prism_layer"
    mesh["note"] = "prism_plate_layer with the first prism's top face collapsed (sliver)"
    mesh["exact_volume"] = float("nan")  # the collapsed cell has no volume to compare with
    mesh["expected_gate"] = "FAIL"
    return mesh


def inverted_pyramid_fan() -> dict:
    """The same six-pyramid fan with each base wound the hexahedron way (inverted cells).

    Kept as evidence that the sign of ``min_sicn`` / ``scaled_jacobian`` is load bearing:
    every other metric of this mesh is identical to ``pyramid_cube_fan``.
    """
    mesh = pyramid_cube()
    outward = (
        (0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)
    )
    mesh["elements"] = [tuple(face) + (8,) for face in outward]
    mesh["name"] = "inverted_pyramid_fan"
    mesh["note"] = "pyramid_cube_fan with the base wound outward: all six cells inverted"
    mesh["expected_gate"] = "FAIL"
    return mesh


MESHES = (prism_plate, pyramid_cube, sheared_prisms, collapsed_prism, inverted_pyramid_fan)


# --------------------------------------------------------------------------------------
# checks
# --------------------------------------------------------------------------------------
def _route(kind: str, points, elements):
    if kind == "prism":
        return reference_prisms(points, elements)
    return reference_pyramids(points, elements)


def audit_mesh(builder) -> dict:
    mesh = builder()
    kind = mesh["kind"]
    points = [tuple(map(float, p)) for p in mesh["points"]]
    elements = [tuple(int(i) for i in cell) for cell in mesh["elements"]]
    quality = measure(points, {kind: elements})[kind]
    reference = _route(kind, points, elements)

    per_metric: dict[str, dict[str, float]] = {}
    worst = 0.0
    for metric, values in reference.items():
        if metric not in quality.metrics:
            continue
        differences = [abs(a - b) for a, b in zip(quality[metric], values)]
        scale = max(1.0, max(abs(v) for v in values))
        worst = max(worst, max(differences))
        per_metric[metric] = {
            "max_abs_difference": max(differences),
            "relative": max(differences) / scale,
            "values": [round(v, 12) for v in values],
            "library_values": [round(v, 12) for v in quality[metric]],
        }

    summed_volume = math.fsum(quality["volume"])
    exact = mesh["exact_volume"]
    volume_error = None if not math.isfinite(exact) else abs(summed_volume - exact)

    report = DEFAULT_GATE.evaluate(measure(points, {kind: elements}))
    failures = [
        {
            "name": check.name,
            "metric": check.metric,
            "kind": check.kind,
            "value": check.value,
            "limit": check.limit,
            "detail": check.detail,
        }
        for check in report.failures
    ]

    expected_gate = mesh["expected_gate"]
    criteria = {
        "second_route_within_1e-6": worst <= TOLERANCE,
        "volume_within_1e-12_of_closed_form": (volume_error is None) or volume_error <= 1e-12,
        "gate_verdict_as_expected": ("PASS" if report.passed else "FAIL") == expected_gate,
    }
    return {
        "name": mesh["name"],
        "kind": kind,
        "note": mesh["note"],
        "elements": len(elements),
        "metrics": per_metric,
        "max_abs_difference_second_route": worst,
        "summed_volume": summed_volume,
        "exact_volume": exact,
        "volume_error": volume_error,
        "gate_passed": report.passed,
        "gate_expected": expected_gate,
        "gate_failures": failures,
        "criteria": criteria,
        "min_sicn_per_element": [round(v, 12) for v in quality["min_sicn"]],
        "skewness_per_element": [round(v, 12) for v in quality["skewness"]],
        "edge_ratio_per_element": [round(v, 12) for v in quality["edge_ratio"]],
    }


def ideal_element_records() -> list[dict]:
    """Analytic values of the ideal prism and the ideal pyramid, against the library."""
    from tensorcad.meshquality import prism_metrics, pyramid_metrics

    h = math.sqrt(3.0) / 2.0
    prism = ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.5, h, 0.0),
             (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (0.5, h, 1.0))
    pyramid = ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0), (0.5, 0.5, 1.0))
    prism_q = {k: v[0] for k, v in prism_metrics(list(prism), [(0, 1, 2, 3, 4, 5)]).metrics.items()}
    pyramid_q = {
        k: v[0] for k, v in pyramid_metrics(list(pyramid), [(0, 1, 2, 3, 4)]).metrics.items()
    }
    apex_angle = math.degrees(2.0 * math.asin(0.5 / math.sqrt(1.5)))
    return [
        {
            "element": "ideal prism (equilateral base, unit height)",
            "closed_form": {"volume": math.sqrt(3.0) / 4.0, "min_sicn": 1.0,
                            "scaled_jacobian": math.sqrt(3.0) / 2.0, "min_angle": 60.0,
                            "max_angle": 90.0, "skewness": 0.0, "edge_ratio": 1.0},
            "library": prism_q,
        },
        {
            "element": "ideal pyramid (unit square base, apex at unit height)",
            "closed_form": {"volume": 1.0 / 3.0, "min_sicn": 1.0,
                            "scaled_jacobian": math.sqrt(2.0 / 3.0), "min_angle": apex_angle,
                            "max_angle": 90.0, "skewness": (60.0 - apex_angle) / 60.0,
                            "edge_ratio": math.sqrt(1.5)},
            "library": pyramid_q,
        },
    ]


# --------------------------------------------------------------------------------------
# figure
# --------------------------------------------------------------------------------------
def make_figure(records: list[dict], ideals: list[dict], path: Path) -> None:
    names = [record["name"] for record in records]
    fig, axes = plt.subplots(1, 3, figsize=(16.0, 5.0))

    # (a) worst element of each mesh: min_sicn and skewness
    ax = axes[0]
    x = np.arange(len(records))
    sicn = [min(record["min_sicn_per_element"]) for record in records]
    skew = [max(record["skewness_per_element"]) for record in records]
    ax.bar(x - 0.2, sicn, width=0.4, label="min min_sicn (1 = ideal)", color="#2b6cb0")
    ax.bar(x + 0.2, skew, width=0.4, label="max equiangle skewness", color="#c05621")
    ax.axhline(0.1, color="#2b6cb0", linestyle="--", linewidth=1.0,
               label="gate: min_sicn >= 0.1")
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=18, ha="right", fontsize=8)
    ax.set_ylabel("worst element of the mesh")
    ax.set_title("(a) the numbers the gate acts on")
    ax.legend(fontsize=7)
    ax.grid(axis="y", alpha=0.3)

    # (b) library vs the second route, per metric, log scale
    ax = axes[1]
    metrics = list(records[0]["metrics"].keys())
    width = 0.8 / max(1, len(records))
    for index, record in enumerate(records):
        diffs = [max(record["metrics"][m]["max_abs_difference"], 1e-18) for m in metrics]
        ax.bar(np.arange(len(metrics)) + index * width, diffs, width=width,
               label=record["name"])
    ax.set_yscale("log")
    ax.axhline(TOLERANCE, color="black", linestyle="--", linewidth=1.0,
               label="acceptance 1e-6")
    ax.set_xticks(np.arange(len(metrics)) + 0.4 - width / 2)
    ax.set_xticklabels(metrics, rotation=40, ha="right", fontsize=7)
    ax.set_ylabel("max |library - second route|")
    ax.set_title("(b) second route (FD Jacobian / Cayley-Menger / atan2)")
    ax.legend(fontsize=6)
    ax.grid(axis="y", alpha=0.3)

    # (c) per-element min_sicn of the boundary-layer mesh and the collapsed one
    ax = axes[2]
    for record in (records[0], records[3]):
        values = record["min_sicn_per_element"]
        ax.plot(np.arange(len(values)), values, marker="o", markersize=3,
                label=record["name"], linewidth=1.2)
    ax.axhline(0.1, color="black", linestyle="--", linewidth=1.0, label="gate limit 0.1")
    ax.set_yscale("symlog", linthresh=0.1)
    ax.set_xlabel("element index")
    ax.set_ylabel("min_sicn")
    ax.set_title("(c) a collapsed prism is not a rounding difference")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)

    fig.suptitle("TensorCAD mesh quality: prism (wedge6) and pyramid metrics, "
                 "cross-checked against an independent route", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(path, dpi=160)
    plt.close(fig)


# --------------------------------------------------------------------------------------
def main() -> int:
    ARTIFACTS.mkdir(parents=True, exist_ok=True)
    records = [audit_mesh(builder) for builder in MESHES]
    ideals = ideal_element_records()

    ideal_failures = [
        {
            "element": record["element"],
            "metric": metric,
            "closed_form": value,
            "library": record["library"][metric],
            "difference": abs(value - record["library"][metric]),
        }
        for record in ideals
        for metric, value in record["closed_form"].items()
        if abs(value - record["library"][metric]) > 1e-9
    ]

    payload = {
        "domain": "mesh quality metrics, prism + pyramid",
        "tolerance_second_route": TOLERANCE,
        "meshes": records,
        "ideal_elements": ideals,
        "ideal_element_failures": ideal_failures,
        "node_table": {"prism": 6, "pyramid": 5},
    }
    audit_path = ARTIFACTS / "audit.json"
    audit_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    figure_path = ARTIFACTS / "quality.png"
    make_figure(records, ideals, figure_path)

    failed = [
        record["name"] for record in records if not all(record["criteria"].values())
    ]
    manifest = {
        "files": {
            path.name: {
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                "bytes": path.stat().st_size,
            }
            for path in (audit_path, figure_path)
        },
        "criteria_failed": failed,
        "ideal_element_failures": ideal_failures,
    }
    (ARTIFACTS / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    for record in records:
        volume_error = record["volume_error"]
        volume_text = "n/a" if volume_error is None else format(volume_error, ".3e")
        failing = [item["name"] for item in record["gate_failures"]]
        verdict = "PASS" if record["gate_passed"] else "FAIL"
        print(
            "{name:>24s}  n={n:>3d}  max|lib-ref|={diff:.3e}  vol_err={volume}  "
            "gate={verdict}  failing={failing}".format(
                name=record["name"],
                n=record["elements"],
                diff=record["max_abs_difference_second_route"],
                volume=volume_text,
                verdict=verdict,
                failing=failing,
            )
        )
    print(f"\nartifacts: {audit_path.relative_to(ROOT)}, {figure_path.relative_to(ROOT)}")
    if failed or ideal_failures:
        print(f"FAILED criteria: {failed} ideal={ideal_failures}")
        return 1
    print("all criteria satisfied")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())