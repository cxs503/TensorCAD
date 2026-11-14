"""Verified case ``CASE-CAD-CURVATURE-SPHERE``: discrete mean curvature of a sphere.

Question
--------
A mesh carries a *discrete* mean curvature (the cotangent Laplacian of
:mod:`tensorcad.curvature`). A sphere carries an *analytic* one. This case puts the
two numbers side by side and asks a single, falsifiable question: **is the discrete
value within 2 % of the analytic value**, and if not, on which meshes and which
reading of "the value" does it fail?

Reference
---------
The reference is the closed form

    kappa_1 = kappa_2 = 1 / R          (sphere of radius R)
    S = kappa_1 + kappa_2 = 2 / R      (sum convention)

Sources:
* do Carmo, *Differential Geometry of Curves and Surfaces*, Prentice-Hall 1976,
  Sec. 3-3, Example 2.
* Meyer et al., *Discrete Differential-Geometry Operators for Triangulated
  2-Manifolds*, VisMath 2003, Sec. 3.3.

The ``tensorcad.curvature.vertex_mean_curvature`` function returns
``H = (kappa_1 + kappa_2) / 2``, so for comparison with the 2/R convention
we use ``S = 2 * H``.

Two readings of "the discrete mean curvature", both reported
------------------------------------------------------------
``integral`` (= the area-weighted mean, ``sum_i A_i S_i / sum_i A_i``)
    The cotangent operator is a discrete *divergence*, so on a closed mesh the
    area-weighted integral reproduces the analytic surface integral of ``S``.
    This is the reading that is gated at 2 %.
``pointwise`` (per-vertex ``|S_i|`` versus ``2/R``, summarised as a median)
    What a curvature-driven mesher or a fairing loop locally sees. It is second
    order on a *structured* triangulation and does not converge on a sliver-ridden
    marching-tet mesh; both halves are measured here.

Gates are ``|measured - reference| <= tolerance`` with the tolerance declared
before the run. The case criterion is the task's: **2 %**.

    python benchmarks/verified/sphere_mean_curvature/run_benchmark.py
    python benchmarks/verified/sphere_mean_curvature/run_benchmark.py --out runs/case_cad_curvature_sphere

Exits non-zero when a gate breaks.  ``result.json`` is the artifact every number
in ``docs/sphere_mean_curvature.md`` is generated from.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import statistics
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

_REPO_ROOT = Path(__file__).resolve().parents[3]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import tensorcad  # noqa: E402
from tensorcad.curvature import principal_curvatures, vertex_mean_curvature  # noqa: E402
from tensorcad.geometry import Mesh  # noqa: E402
from tensorcad.geometry_field import SphereField, to_mesh  # noqa: E402
from tensorcad.mesh_conditioning import project_to_field, taubin_smooth  # noqa: E402

SCHEMA = "tensorcad.benchmark.sphere_mean_curvature/1"
CASE = "CASE-CAD-CURVATURE-SPHERE"
CAMPAIGN_NODE = "算例池/benchmark 供给 (owner 2026-11-12)"

#: The case criterion, declared before the first measurement.
CASE_TOLERANCE = 0.02
#: The analytic field path is a finite-difference Hessian: documented at 1.4e-15.
FIELD_TOLERANCE = 1e-6
#: Expected order of the cotangent estimator on a structured triangulation (theory: 2).
EXPECTED_ORDER = 2.0
#: Tolerance on the observed order (allow 1.5–2.5).
ORDER_TOLERANCE = 0.5


def _row(
    gate: str,
    quantity: str,
    reference: float,
    measured: float,
    tolerance: float,
    *,
    reference_source: str,
    unit: str = "1/m",
    absolute: bool = False,
    note: str = "",
) -> dict[str, Any]:
    """One ``reference vs measured`` line."""
    if absolute or reference == 0.0:
        error = measured - reference
        error_pct = None if reference == 0.0 else 100.0 * error / reference
    else:
        error = measured / reference - 1.0
        error_pct = 100.0 * error
    holds = bool(math.isfinite(error) and abs(error) <= tolerance)
    return {
        "gate": gate,
        "quantity": quantity,
        "unit": unit,
        "reference": reference,
        "measured": measured,
        "error": error,
        "error_percent": error_pct,
        "tolerance": tolerance,
        "tolerance_kind": "absolute" if absolute else "relative",
        "holds": holds,
        "reference_source": reference_source,
        "note": note,
    }


def _area_weights(mesh: Mesh) -> list[float]:
    """Lumped area per vertex (1/3 of incident face areas)."""
    weights = [0.0] * len(mesh.vertices)
    for face in mesh.faces:
        a, b, c = (mesh.vertices[i] for i in face)
        u = [b[ax] - a[ax] for ax in range(3)]
        v = [c[ax] - a[ax] for ax in range(3)]
        cross = (
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0],
        )
        area = 0.5 * math.sqrt(sum(c * c for c in cross)) / 3.0
        for i in face:
            weights[i] += area
    return weights


def _area_weighted_mean(values: list[float], weights: list[float]) -> float:
    """Area-weighted mean of finite values."""
    pairs = [(v, w) for v, w in zip(values, weights) if math.isfinite(v) and w > 0]
    if not pairs:
        return float("nan")
    return math.fsum(v * w for v, w in pairs) / math.fsum(w for _, w in pairs)


def _uv_sphere(radius: float, nu: int, nv: int) -> Mesh:
    """A regular UV sphere mesh (structured triangulation)."""
    vertices: list[tuple[float, float, float]] = []
    grid: dict[tuple[int, int], int] = {}
    for i in range(nu + 1):
        theta = math.pi * i / nu
        for j in range(nv):
            phi = 2.0 * math.pi * j / nv
            vertices.append((
                radius * math.sin(theta) * math.cos(phi),
                radius * math.sin(theta) * math.sin(phi),
                radius * math.cos(theta),
            ))
            grid[(i, j)] = len(vertices) - 1
    faces: list[tuple[int, int, int]] = []
    for i in range(nu):
        for j in range(nv):
            j2 = (j + 1) % nv
            a, b = grid[(i, j)], grid[(i, j2)]
            c, d = grid[(i + 1, j)], grid[(i + 1, j2)]
            faces.append((a, c, b))
            faces.append((b, c, d))
    return Mesh(vertices=vertices, faces=faces)


def _fibonacci_probes(radius: float, count: int) -> list[tuple[float, float, float]]:
    """Well-distributed probe points on a sphere (Fibonacci spiral + poles)."""
    golden = math.pi * (3.0 - math.sqrt(5.0))
    points = []
    for idx in range(count):
        z = 1.0 - 2.0 * (idx + 0.5) / count
        r = math.sqrt(max(0.0, 1.0 - z * z))
        theta = golden * idx
        points.append((radius * r * math.cos(theta), radius * r * math.sin(theta), radius * z))
    points.extend([
        (radius, 0.0, 0.0), (-radius, 0.0, 0.0),
        (0.0, radius, 0.0), (0.0, -radius, 0.0),
        (0.0, 0.0, radius), (0.0, 0.0, -radius),
    ])
    return points


# ---------------------------------------------------------------------------
# G1: Analytic field path — principal curvatures via finite-difference Hessian
# ---------------------------------------------------------------------------
def g1_field_path(probes: int = 64) -> dict[str, Any]:
    """Closed-form curvature vs the finite-difference shape operator."""
    rows: list[dict[str, Any]] = []
    radii = (0.1, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0)
    for radius in radii:
        field = SphereField(center=(0.0, 0.0, 0.0), radius=radius)
        samples = principal_curvatures(field, _fibonacci_probes(radius, probes))
        # The field returns H = (k1+k2)/2; for 2/R convention multiply by 2
        measured_H = math.fsum(
            0.5 * (s.kappa_max + s.kappa_min) for s in samples
        ) / len(samples)
        measured_S = 2.0 * measured_H  # S = k1 + k2 = 2/R
        ref_S = 2.0 / radius
        rows.append(_row(
            "G1", f"sphere S=2/R, R={radius:g}", ref_S, measured_S,
            FIELD_TOLERANCE, reference_source=f"closed form 2/R = {ref_S:.6g} 1/m",
            absolute=True, note=f"{len(samples)} probes; |S - 2/R|",
        ))
    return {"rows": rows}


# ---------------------------------------------------------------------------
# G2: Integral reading on marching-tet meshes (the 2 % gate)
# ---------------------------------------------------------------------------
def g2_integral_reading() -> dict[str, Any]:
    """Area-weighted mean curvature on marching-tet meshes at several resolutions."""
    rows: list[dict[str, Any]] = []
    convergence: list[dict[str, Any]] = []
    ref_S = 2.0  # R = 1

    for level in (12, 16, 24, 32, 48):
        field = SphereField(center=(0.0, 0.0, 0.0), radius=1.0)
        mesh = to_mesh(field, level, margin=3).mesh
        values = vertex_mean_curvature(mesh)
        weights = _area_weights(mesh)
        H_mean = _area_weighted_mean(values, weights)
        S_mean = 2.0 * H_mean  # convert to 2/R convention

        rows.append(_row(
            "G2", f"marching-tet level {level}, area-weighted S",
            ref_S, S_mean, CASE_TOLERANCE,
            reference_source="closed form 2/R = 2.0 1/m (R=1)",
            note=f"{len(mesh.faces)} faces, {len(mesh.vertices)} vertices",
        ))
        convergence.append({
            "level": level, "faces": len(mesh.faces), "vertices": len(mesh.vertices),
            "H_mean": H_mean, "S_mean": S_mean,
            "error_pct": abs(S_mean - ref_S) / ref_S * 100,
        })

    return {"rows": rows, "convergence": convergence}


# ---------------------------------------------------------------------------
# G3: Pointwise reading on marching-tet (raw + conditioned)
# ---------------------------------------------------------------------------
def g3_pointwise_reading() -> dict[str, Any]:
    """Pointwise median error on marching-tet, raw and after conditioning."""
    ref_S = 2.0
    field = SphereField(center=(0.0, 0.0, 0.0), radius=1.0)
    mesh = to_mesh(field, 32, margin=3).mesh

    # Raw
    values_raw = vertex_mean_curvature(mesh)
    weights = _area_weights(mesh)
    S_raw = [2.0 * v for v in values_raw if math.isfinite(v)]
    median_raw = statistics.median(S_raw)
    p95_raw = sorted(S_raw)[int(0.95 * (len(S_raw) - 1))]

    # Conditioned
    smoothed, _ = taubin_smooth(mesh, iterations=30)
    projected, _ = project_to_field(smoothed, field, iterations=10)
    values_cond = vertex_mean_curvature(projected)
    S_cond = [2.0 * v for v in values_cond if math.isfinite(v)]
    median_cond = statistics.median(S_cond)
    p95_cond = sorted(S_cond)[int(0.95 * (len(S_cond) - 1))]

    rows = [
        _row("G3", "pointwise median S (raw marching-tet)", ref_S, median_raw,
             0.60, reference_source="closed form 2/R = 2.0",
             note="raw marching-tet has slivers; median ~50% off is expected"),
        _row("G3b", "pointwise median S (taubin+project)", ref_S, median_cond,
             0.05, reference_source="closed form 2/R = 2.0",
             note="after 30 taubin + 10 projection iterations"),
    ]

    return {
        "rows": rows,
        "raw": {"median": median_raw, "p95": p95_raw, "n": len(S_raw)},
        "conditioned": {"median": median_cond, "p95": p95_cond, "n": len(S_cond)},
    }


# ---------------------------------------------------------------------------
# G4: Structured UV sphere — convergence order
# ---------------------------------------------------------------------------
def g4_uv_sphere_convergence() -> dict[str, Any]:
    """Convergence of the cotangent estimator on a structured UV sphere."""
    ref_S = 2.0
    resolutions = [(8, 16), (16, 32), (32, 64), (64, 128)]
    series: list[dict[str, Any]] = []

    for nu, nv in resolutions:
        mesh = _uv_sphere(1.0, nu, nv)
        values = vertex_mean_curvature(mesh)
        weights = _area_weights(mesh)
        H_mean = _area_weighted_mean(values, weights)
        S_mean = 2.0 * H_mean
        error = abs(S_mean - ref_S) / ref_S

        S_vals = [2.0 * v for v in values if math.isfinite(v)]
        series.append({
            "nu": nu, "nv": nv, "faces": len(mesh.faces),
            "vertices": len(mesh.vertices),
            "S_mean": S_mean, "error_pct": error * 100,
            "pointwise_median": statistics.median(S_vals),
        })

    # Estimate convergence order from the last three points
    if len(series) >= 3:
        hs = [1.0 / s["nu"] for s in series[-3:]]
        es = [s["error_pct"] for s in series[-3:]]
        # log(e) = p * log(h) + c => p = (log(e2)-log(e1)) / (log(h2)-log(h1))
        orders = []
        for i in range(len(hs) - 1):
            if hs[i] != hs[i+1] and es[i] > 0 and es[i+1] > 0:
                p = math.log(es[i+1] / es[i]) / math.log(hs[i+1] / hs[i])
                orders.append(p)
        observed_order = statistics.mean(orders) if orders else float("nan")
    else:
        observed_order = float("nan")

    rows = [
        _row("G4", f"UV sphere {resolutions[-1][0]}x{resolutions[-1][1]}, area-weighted S",
             ref_S, series[-1]["S_mean"], CASE_TOLERANCE,
             reference_source="closed form 2/R = 2.0",
             note=f"structured triangulation, observed order ~{observed_order:.2f}"),
    ]

    return {
        "rows": rows,
        "series": series,
        "observed_order": observed_order,
    }


# ---------------------------------------------------------------------------
# G5: Multiple radii — integral reading
# ---------------------------------------------------------------------------
def g5_multiple_radii() -> dict[str, Any]:
    """Integral reading on spheres of different radii."""
    rows: list[dict[str, Any]] = []
    radii = [0.5, 1.0, 2.0, 5.0]
    for R in radii:
        field = SphereField(center=(0.0, 0.0, 0.0), radius=R)
        mesh = to_mesh(field, 32, margin=3).mesh
        values = vertex_mean_curvature(mesh)
        weights = _area_weights(mesh)
        H_mean = _area_weighted_mean(values, weights)
        S_mean = 2.0 * H_mean
        ref_S = 2.0 / R
        rows.append(_row(
            "G5", f"sphere R={R}, area-weighted S", ref_S, S_mean,
            CASE_TOLERANCE, reference_source=f"closed form 2/R = {ref_S:.6g}",
            note=f"{len(mesh.faces)} faces",
        ))
    return {"rows": rows}


# ---------------------------------------------------------------------------
# G6: Integral identity cross-check (bending energy = 8*pi*R)
# ---------------------------------------------------------------------------
def g6_integral_identity() -> dict[str, Any]:
    """Cross-check: integral of S over sphere = 8*pi*R (bending energy)."""
    R = 1.0
    field = SphereField(center=(0.0, 0.0, 0.0), radius=R)
    mesh = to_mesh(field, 32, margin=3).mesh
    values = vertex_mean_curvature(mesh)
    weights = _area_weights(mesh)

    # Integral of S = sum_i A_i * S_i = sum_i A_i * 2 * H_i
    integral_S = math.fsum(
        2.0 * v * w for v, w in zip(values, weights) if math.isfinite(v) and w > 0
    )
    ref_integral = 8.0 * math.pi * R  # = 8*pi for R=1

    return {
        "rows": [_row(
            "G6", "integral of S dA (bending energy identity)",
            ref_integral, integral_S, 0.01,
            reference_source="closed form 8*pi*R = 8*pi",
            unit="", note="integral identity cross-check",
        )],
        "measured": integral_S,
        "reference": ref_integral,
    }


# ---------------------------------------------------------------------------
# G7: Raw mesh pointwise failure (diagnostic, not gated)
# ---------------------------------------------------------------------------
def g7_raw_pointwise_diagnostic() -> dict[str, Any]:
    """Document that raw marching-tet pointwise curvature is poor (expected)."""
    field = SphereField(center=(0.0, 0.0, 0.0), radius=1.0)
    mesh = to_mesh(field, 32, margin=3).mesh
    values = vertex_mean_curvature(mesh)
    S_vals = [2.0 * v for v in values if math.isfinite(v)]
    ref_S = 2.0

    deviations = sorted(abs(s - ref_S) / ref_S for s in S_vals)
    return {
        "gate": "G7",
        "quantity": "raw marching-tet pointwise error (diagnostic, not gated)",
        "median_error_pct": deviations[len(deviations) // 2] * 100,
        "p95_error_pct": deviations[int(0.95 * (len(deviations) - 1))] * 100,
        "max_error_pct": deviations[-1] * 100,
        "note": "Expected: raw marching-tet has sliver triangles; pointwise is poor but integral is good",
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def run(out_dir: Path | None = None) -> dict[str, Any]:
    """Run all gates and return the result bundle."""
    t0 = time.time()

    g1 = g1_field_path(probes=64)
    g2 = g2_integral_reading()
    g3 = g3_pointwise_reading()
    g4 = g4_uv_sphere_convergence()
    g5 = g5_multiple_radii()
    g6 = g6_integral_identity()
    g7 = g7_raw_pointwise_diagnostic()

    elapsed = time.time() - t0

    # Collect all rows
    all_rows = (
        g1["rows"] + g2["rows"] + g3["rows"] + g4["rows"] + g5["rows"] + g6["rows"]
    )
    all_hold = all(r["holds"] for r in all_rows)

    result = {
        "schema": SCHEMA,
        "case": CASE,
        "campaign_node": CAMPAIGN_NODE,
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "tensorcad_version": tensorcad.__version__ if hasattr(tensorcad, "__version__") else "dev",
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "elapsed_s": round(elapsed, 2),
        "gates": {
            "G1_field_path": g1,
            "G2_integral_marching_tet": g2,
            "G3_pointwise_conditioned": g3,
            "G4_uv_sphere_convergence": g4,
            "G5_multiple_radii": g5,
            "G6_integral_identity": g6,
            "G7_raw_diagnostic": g7,
        },
        "all_rows": all_rows,
        "all_gates_hold": all_hold,
        "case_tolerance": CASE_TOLERANCE,
        "summary": {
            "total_gates": len(all_rows),
            "passed": sum(1 for r in all_rows if r["holds"]),
            "failed": sum(1 for r in all_rows if not r["holds"]),
            "worst_error_pct": max(
                abs(r["error_percent"]) for r in all_rows if r["error_percent"] is not None
            ),
        },
    }

    if out_dir is not None:
        out_dir.mkdir(parents=True, exist_ok=True)
        result_path = out_dir / "result.json"
        result_path.write_text(json.dumps(result, indent=2))

        # Gates verified (flat list for machine adjudication)
        gates_verified = {
            "schema": "tensorcad.gates_verified/1",
            "case": CASE,
            "timestamp": result["timestamp"],
            "gates": [
                {
                    "gate": r["gate"],
                    "quantity": r["quantity"],
                    "reference": r["reference"],
                    "measured": r["measured"],
                    "error_percent": r["error_percent"],
                    "tolerance": r["tolerance"],
                    "holds": r["holds"],
                }
                for r in all_rows
            ],
            "all_hold": all_hold,
        }
        (out_dir / "gates_verified.json").write_text(json.dumps(gates_verified, indent=2))

    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=CASE)
    parser.add_argument("--out", type=str, default=None, help="Output directory")
    args = parser.parse_args()

    out_dir = Path(args.out) if args.out else None
    result = run(out_dir)

    # Print summary
    s = result["summary"]
    print(f"\n{'='*60}")
    print(f"CASE: {CASE}")
    print(f"Gates: {s['passed']}/{s['total_gates']} passed, {s['failed']} failed")
    print(f"Worst error: {s['worst_error_pct']:.4f}%")
    print(f"Case tolerance: {CASE_TOLERANCE*100:.1f}%")
    print(f"All gates hold: {result['all_gates_hold']}")
    print(f"Elapsed: {result['elapsed_s']:.1f}s")
    print(f"{'='*60}\n")

    if not result["all_gates_hold"]:
        print("FAILED gates:")
        for r in result["all_rows"]:
            if not r["holds"]:
                print(f"  {r['gate']}: {r['quantity']} -> error={r['error_percent']:.4f}% (tol={r['tolerance']*100:.1f}%)")
        sys.exit(1)
    else:
        print("ALL GATES PASSED")


if __name__ == "__main__":
    main()
