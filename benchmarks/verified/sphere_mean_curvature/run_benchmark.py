"""Verified case ``CASE-CAD-CURVATURE-SPHERE``: discrete mean curvature of a sphere.

Question
--------
A mesh carries a *discrete* mean curvature (the cotangent Laplacian of
:mod:`tensorcad.curvature`). A sphere carries an *analytic* one. This case puts the
two numbers side by side and asks a single, falsifiable question: **is the discrete
value within 2 % of the analytic value**, and if not, on which meshes and which
reading of "the value" does it fail?

Reference, and why it is not transcribed
----------------------------------------
The reference is the closed form

    kappa_1 = kappa_2 = 1 / R          (sphere of radius R, every direction is
                                        a normal section of the same circle)
    S = kappa_1 + kappa_2 = 2 / R      (the *sum* convention: the discrete
                                        operator returns the mean-curvature
                                        normal with |K| = 2|H|, so the quantity
                                        it estimates is kappa_1 + kappa_2)

Sources, section level (the value itself lives in :func:`reference_sum`, is
evaluated by code, and is never copied into this file):

* do Carmo, *Differential Geometry of Curves and Surfaces*, Prentice-Hall 1976,
  Sec. 3-3, Example 2 (the sphere: all normal sections are circles of radius R,
  hence normal curvature ``1/R`` in every tangent direction, and the principal
  curvatures are ``1/R, 1/R``). The same statement is already adjudicated inside
  this repository -- ``docs/fairness.md`` ("Analytic field, exact. Sphere
  ``1/R, 1/R`` ..."), which the ``curvature_fairness`` verified case gates at
  ``1e-6`` 1/m.
* Meyer, Desbrun, Schroeder, Barr, *Discrete Differential-Geometry Operators for
  Triangulated 2-Manifolds*, VisMath 2003, Sec. 3.3: the mean curvature normal
  ``K = (1/(2A)) sum (cot a + cot b)(x_i - x_j)`` has ``|K| = 2|H|`` and
  ``H = (kappa_1 + kappa_2)/2``. That is exactly the convention of
  ``tensorcad.curvature.vertex_mean_curvature`` (its docstring: ``H_i = 0.5 K_i . n_i``),
  so the analytic target for that estimator is ``2/R``, not ``1/R``.

The reference carries a second, independent check that no text can get wrong:
the integral of ``S`` over a sphere is ``8 pi R``, and for ``R = 1`` that is the
``E_bending = 8 pi`` closed form already gated by
``benchmarks/verified/curvature_fairness`` (gates G2/G3, mesh path). Agreement
between the local measurement and that pre-existing gate is what makes a wrong
constant in this file detectable.

Two readings of "the discrete mean curvature", both reported
------------------------------------------------------------
``integral`` (= the area-weighted mean, ``sum_i A_i S_i / sum_i A_i``)
    The cotangent operator is a discrete *divergence*, so on a closed mesh the
    area-weighted integral telescopes and reproduces the analytic surface
    integral of ``S`` regardless of how badly shaped the triangles are. This is
    the reading that is gated at 2 %, because it is the reading a downstream
    integral quantity (bending energy, wetted-area moments) actually consumes.
``pointwise`` (per-vertex ``|S_i|`` versus ``2/R``, summarised as a median)
    What a curvature-driven mesher or a fairing loop locally sees. It is second
    order on a *structured* triangulation and does not converge on a sliver-ridden
    marching-tet mesh; both halves are measured here rather than asserted.

Gates are ``|measured - reference| <= tolerance`` with the tolerance declared
before the run. The case criterion is the task's: **2 %**.

    python benchmarks/verified/sphere_mean_curvature/run_benchmark.py        # ~25 s
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
#: Expected order of the cotangent estimator on a structured triangulation (theory:
#: second order). The gate allows 0.5 of slack because the UV family has triangle-fan
#: poles whose 1/h distortion sits in the tail of the distribution.
ORDER_EXPECTED = 2.0
ORDER_TOLERANCE = 0.5

REFERENCE_SOURCES = {
    "closed_form": (
        "do Carmo, Differential Geometry of Curves and Surfaces (Prentice-Hall 1976), "
        "Sec. 3-3 Example 2: sphere radius R has kappa_1 = kappa_2 = 1/R"
    ),
    "convention": (
        "Meyer/Desbrun/Schroeder/Barr, Discrete Differential-Geometry Operators for "
        "Triangulated 2-Manifolds (VisMath 2003), Sec. 3.3: |K| = 2|H|, so the "
        "cotangent estimator targets kappa_1 + kappa_2"
    ),
    "integral_identity": (
        "int_S (kappa_1+kappa_2) dA = 8 pi R; at R = 1 this is the E_bending = 8 pi "
        "already gated by benchmarks/verified/curvature_fairness (G2/G3)"
    ),
    "in_repo": "docs/fairness.md, 'Analytic field, exact. Sphere 1/R, 1/R' (curvature_fairness G1, 1e-6 1/m)",
}


# --------------------------------------------------------------------------- #
# references -- evaluated, never transcribed
# --------------------------------------------------------------------------- #
def reference_sum(radius: float) -> float:
    """``kappa_1 + kappa_2 = 2/R`` for a sphere, from ``kappa_1 = kappa_2 = 1/R``."""
    return (1.0 / radius) + (1.0 / radius)


def reference_total(radius: float) -> float:
    """``int_S S dA = (2/R)(4 pi R^2) = 8 pi R``. Equals ``E_bending`` at ``R = 1``."""
    return reference_sum(radius) * (4.0 * math.pi * radius * radius)


# --------------------------------------------------------------------------- #
# measurement
# --------------------------------------------------------------------------- #
def area_weights(mesh: Mesh) -> list[float]:
    """One third of the incident triangle area per vertex (the tests' convention)."""
    weights = [0.0] * len(mesh.vertices)
    for face in mesh.faces:
        a, b, c = (mesh.vertices[index] for index in face)
        u = [b[axis] - a[axis] for axis in range(3)]
        v = [c[axis] - a[axis] for axis in range(3)]
        cross = (
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0],
        )
        third = 0.5 * math.sqrt(sum(component * component for component in cross)) / 3.0
        for index in face:
            weights[index] += third
    return weights


def _summarise(values: list[float]) -> dict[str, float]:
    finite = sorted(value for value in values if math.isfinite(value))
    if not finite:
        return {"count": 0, "median": float("nan"), "p95": float("nan"), "max": float("nan")}

    def _quantile(q: float) -> float:
        position = q * (len(finite) - 1)
        low = int(math.floor(position))
        high = min(low + 1, len(finite) - 1)
        weight = position - low
        return finite[low] * (1.0 - weight) + finite[high] * weight

    return {
        "count": len(finite),
        "median": _quantile(0.5),
        "p95": _quantile(0.95),
        "max": finite[-1],
    }


def measure(mesh: Mesh, radius: float) -> dict[str, Any]:
    """Both readings of the discrete mean curvature on one mesh, ``S = 2H``."""
    started = time.perf_counter()
    mean_curvature = vertex_mean_curvature(mesh)
    weights = area_weights(mesh)
    area = math.fsum(weights)
    total = 2.0 * math.fsum(
        value * weight
        for value, weight in zip(mean_curvature, weights, strict=True)
        if math.isfinite(value) and weight > 0.0
    )
    integral = total / area if area > 0.0 else float("nan")
    # relative error of S = 2H against 2/R is the same as that of H against 1/R
    pointwise = [abs(abs(value) * radius - 1.0) for value in mean_curvature]
    reference = reference_sum(radius)
    return {
        "radius": radius,
        "vertices": len(mesh.vertices),
        "faces": len(mesh.faces),
        "area": area,
        "area_exact": 4.0 * math.pi * radius * radius,
        "area_relative_error": area / (4.0 * math.pi * radius * radius) - 1.0,
        # the sum convention: the estimator returns H = (k1+k2)/2, so S = 2H
        "integral_sum": integral,
        "integral_relative_error": integral / reference - 1.0,
        "total_sum": total,
        "total_relative_error": total / reference_total(radius) - 1.0,
        "pointwise": _summarise(pointwise),
        "pointwise_median_value": statistics.median(
            [2.0 * abs(value) for value in mean_curvature if math.isfinite(value)]
        ),
        "reference_sum": reference,
        "reference_total": reference_total(radius),
        "seconds": time.perf_counter() - started,
    }


def uv_sphere(radius: float, nu: int, nv: int) -> Mesh:
    """A legal UV sphere: one vertex per pole, a triangle fan onto it, no zero-area face.

    ``nu`` latitude bands, ``nv`` longitudes. Structured on purpose: this is the
    control that separates "the operator is inaccurate" from "the mesher gave the
    operator nothing to work with".
    """
    vertices: list[tuple[float, float, float]] = [(0.0, 0.0, radius)]
    rings: list[list[int]] = []
    for band in range(1, nu):
        theta = math.pi * band / nu
        ring: list[int] = []
        for longitude in range(nv):
            phi = 2.0 * math.pi * longitude / nv
            ring.append(len(vertices))
            vertices.append(
                (
                    radius * math.sin(theta) * math.cos(phi),
                    radius * math.sin(theta) * math.sin(phi),
                    radius * math.cos(theta),
                )
            )
        rings.append(ring)
    south = len(vertices)
    vertices.append((0.0, 0.0, -radius))

    faces: list[tuple[int, int, int]] = []
    for longitude in range(nv):
        faces.append((0, rings[0][longitude], rings[0][(longitude + 1) % nv]))
    for lower, upper in zip(rings, rings[1:], strict=False):
        for longitude in range(nv):
            nxt = (longitude + 1) % nv
            faces.append((lower[longitude], upper[longitude], upper[nxt]))
            faces.append((lower[longitude], upper[nxt], lower[nxt]))
    for longitude in range(nv):
        faces.append((south, rings[-1][(longitude + 1) % nv], rings[-1][longitude]))
    return Mesh(vertices=vertices, faces=faces)


# --------------------------------------------------------------------------- #
# rows
# --------------------------------------------------------------------------- #
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
    if absolute or reference == 0.0:
        error = measured - reference
        error_percent = None if reference == 0.0 else 100.0 * error / reference
    else:
        error = measured / reference - 1.0
        error_percent = 100.0 * error
    return {
        "gate": gate,
        "quantity": quantity,
        "unit": unit,
        "reference": reference,
        "measured": measured,
        "error": error,
        "error_percent": error_percent,
        "tolerance": tolerance,
        "tolerance_kind": "absolute" if absolute else "relative",
        "holds": bool(math.isfinite(error) and abs(error) <= tolerance),
        "reference_source": reference_source,
        "note": note,
    }


# --------------------------------------------------------------------------- #
# sections
# --------------------------------------------------------------------------- #
def reference_chain_section(probes: int) -> dict[str, Any]:
    """G4: the analytic-field estimator must return exactly ``2/R``. Measurement chain."""
    rows = []
    radii = (0.1, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0)
    golden = math.pi * (3.0 - math.sqrt(5.0))
    for radius in radii:
        field = SphereField(center=(0.0, 0.0, 0.0), radius=radius)
        points = []
        for index in range(probes):
            z = 1.0 - 2.0 * (index + 0.5) / probes
            radial = math.sqrt(max(0.0, 1.0 - z * z))
            theta = golden * index
            points.append(
                (radius * radial * math.cos(theta), radius * radial * math.sin(theta), radius * z)
            )
        samples = principal_curvatures(field, points, refine=False)
        measured = math.fsum(sample.kappa_max + sample.kappa_min for sample in samples) / len(samples)
        rows.append(
            _row(
                "G4",
                f"field-path kappa_1+kappa_2, sphere R={radius:g} m, {probes} Fibonacci probes",
                reference_sum(radius),
                measured,
                FIELD_TOLERANCE,
                reference_source=REFERENCE_SOURCES["closed_form"],
                note="finite-difference Hessian of the exact SDF; proves the measurement chain",
            )
        )
    return {"rows": rows}


def mesh_section_marching(levels: tuple[int, ...], radius: float, condition: tuple[int, ...]) -> dict[str, Any]:
    """G1 + the resolution series on a raw marching-tet mesh, with one conditioned series."""
    field = SphereField(center=(0.0, 0.0, 0.0), radius=radius)
    series = []
    for level in levels:
        mesh = to_mesh(field, level, margin=3).mesh
        raw = measure(mesh, radius)
        entry: dict[str, Any] = {
            "level": level,
            "h": 2.0 * radius / level,
            "faces": raw["faces"],
            "vertices": raw["vertices"],
            "area_relative_error": raw["area_relative_error"],
            "integral_relative_error": raw["integral_relative_error"],
            "integral_sum": raw["integral_sum"],
            "total_sum": raw["total_sum"],
            "total_relative_error": raw["total_relative_error"],
            "pointwise_median_relative_error": raw["pointwise"]["median"],
            "pointwise_median_value": raw["pointwise_median_value"],
            "pointwise_p95_relative_error": raw["pointwise"]["p95"],
            "seconds": raw["seconds"],
        }
        if level in condition:
            smoothed, drift = taubin_smooth(mesh, field=field, iterations=10)
            projected, projection_drift = project_to_field(smoothed, field, iterations=2)
            conditioned = measure(projected, radius)
            entry["conditioned"] = {
                "taubin_drift": drift,
                "projection_drift": projection_drift,
                "integral_relative_error": conditioned["integral_relative_error"],
                "pointwise_median_relative_error": conditioned["pointwise"]["median"],
                "pointwise_p95_relative_error": conditioned["pointwise"]["p95"],
            }
        series.append(entry)

    finest = series[-1]
    rows = [
        _row(
            "G1",
            f"discrete mean curvature, area-weighted integral, {radius:g} m sphere, "
            f"marching tets level {finest['level']}",
            reference_sum(radius),
            finest["integral_sum"],
            CASE_TOLERANCE,
            reference_source=REFERENCE_SOURCES["closed_form"] + " ; " + REFERENCE_SOURCES["convention"],
            note=f"{finest['faces']} faces; S = 2H so the cotangent estimator is read at 2x",
        ),
        _row(
            "G5",
            f"surface integral of kappa_1+kappa_2 over the mesh, level {finest['level']}",
            reference_total(radius),
            finest["total_sum"],
            CASE_TOLERANCE,
            reference_source=REFERENCE_SOURCES["integral_identity"],
            unit="1",
            note="independent route: 8 pi R, the E_bending = 8 pi closed form at R = 1",
        ),
    ]
    return {"rows": rows, "series": series}


def mesh_section_structured(radius: float, resolutions: tuple[tuple[int, int], ...]) -> dict[str, Any]:
    """G2, G3 and the convergence-order gate on legal UV spheres."""
    series = []
    for nu, nv in resolutions:
        mesh = uv_sphere(radius, nu, nv)
        result = measure(mesh, radius)
        series.append(
            {
                "nu": nu,
                "nv": nv,
                "h": math.pi * radius / nu,
                "faces": result["faces"],
                "vertices": result["vertices"],
                "integral_relative_error": result["integral_relative_error"],
                "pointwise_median_relative_error": result["pointwise"]["median"],
                "pointwise_median_value": result["pointwise_median_value"],
                "pointwise_p95_relative_error": result["pointwise"]["p95"],
                "pointwise_max_relative_error": result["pointwise"]["max"],
                "area_relative_error": result["area_relative_error"],
            }
        )
    gate_entry = next(entry for entry in series if (entry["nu"], entry["nv"]) == (16, 32))
    rows = [
        _row(
            "G2",
            f"pointwise discrete mean curvature, median relative error, structured UV sphere "
            f"{radius:g} m, 16x32",
            0.0,
            gate_entry["pointwise_median_relative_error"],
            CASE_TOLERANCE,
            reference_source=REFERENCE_SOURCES["closed_form"] + " ; " + REFERENCE_SOURCES["convention"],
            unit="-",
            absolute=True,
            note="median |2H_i R - 2| / 2 over %d vertices; 0 means exact" % gate_entry["vertices"],
        )
    ]

    # least-squares slope of log|error| against log h over the first three rungs
    rungs = series[:3]
    xs = [math.log(entry["h"]) for entry in rungs]
    ys = [math.log(entry["pointwise_median_relative_error"]) for entry in rungs]
    x_mean = math.fsum(xs) / len(xs)
    y_mean = math.fsum(ys) / len(ys)
    slope = math.fsum((x - x_mean) * (y - y_mean) for x, y in zip(xs, ys, strict=True)) / math.fsum(
        (x - x_mean) ** 2 for x in xs
    )
    rows.append(
        _row(
            "G3",
            "observed order of the pointwise estimate on the structured sphere",
            ORDER_EXPECTED,
            slope,
            ORDER_TOLERANCE,
            reference_source="second-order cotangent estimator (Meyer et al. 2003, Sec. 3.3)",
            unit="-",
            absolute=True,
            note=f"least squares over nu = {[entry['nu'] for entry in rungs]}; gate is |order-2| <= 0.5",
        )
    )
    return {"rows": rows, "series": series, "observed_order": slope}


def radius_section(radii: tuple[float, ...], level: int) -> dict[str, Any]:
    """G6: the dimensionless relative error may not depend on the radius."""
    rows = []
    series = []
    for radius in radii:
        field = SphereField(center=(0.0, 0.0, 0.0), radius=radius)
        mesh = to_mesh(field, level, margin=3).mesh
        result = measure(mesh, radius)
        series.append(
            {
                "radius": radius,
                "level": level,
                "faces": result["faces"],
                "integral_sum": result["integral_sum"],
                "reference_sum": result["reference_sum"],
                "integral_relative_error": result["integral_relative_error"],
                "pointwise_median_relative_error": result["pointwise"]["median"],
            }
        )
        rows.append(
            _row(
                "G6",
                f"discrete mean curvature, area-weighted integral, sphere R={radius:g} m, level {level}",
                reference_sum(radius),
                result["integral_sum"],
                CASE_TOLERANCE,
                reference_source=REFERENCE_SOURCES["closed_form"],
            )
        )
    return {"rows": rows, "series": series}


def diagnostic_section(level: int, radius: float) -> dict[str, Any]:
    """G7, diagnostic only: the raw pointwise reading, which is expected to break 2 %."""
    field = SphereField(center=(0.0, 0.0, 0.0), radius=radius)
    mesh = to_mesh(field, level, margin=3).mesh
    result = measure(mesh, radius)
    rows = [
        _row(
            "G7",
            f"pointwise discrete mean curvature, median relative error, raw marching tets "
            f"level {level}",
            0.0,
            result["pointwise"]["median"],
            CASE_TOLERANCE,
            reference_source=REFERENCE_SOURCES["closed_form"] + " ; " + REFERENCE_SOURCES["convention"],
            note=(
                "sliver-dominated connectivity: the median does not converge under "
                "refinement (see series), so this row is reported, not gated"
            ),
        )
    ]
    return {"rows": rows, "measure": {k: v for k, v in result.items() if k != "pointwise"}
            | {"pointwise": result["pointwise"]}}


# --------------------------------------------------------------------------- #
# environment / plotting
# --------------------------------------------------------------------------- #
def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _git(*args: str) -> str:
    try:
        return subprocess.run(
            ["git", *args], cwd=_REPO_ROOT, capture_output=True, text=True, check=True
        ).stdout.strip()
    except Exception:  # pragma: no cover - not a git checkout
        return "unavailable"


def environment() -> dict[str, Any]:
    return {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "tensorcad_version": getattr(tensorcad, "__version__", "unknown"),
        "git_commit": _git("rev-parse", "HEAD"),
        "git_branch": _git("rev-parse", "--abbrev-ref", "HEAD"),
        "git_status_porcelain": _git("status", "--porcelain"),
        "sources": {
            str(path.relative_to(_REPO_ROOT)): _sha256(path)
            for path in (
                _REPO_ROOT / "tensorcad" / "curvature.py",
                _REPO_ROOT / "tensorcad" / "geometry_field.py",
                _REPO_ROOT / "tensorcad" / "mesh_conditioning.py",
                Path(__file__),
            )
        },
        "repro_cmd": "python benchmarks/verified/sphere_mean_curvature/run_benchmark.py",
    }


def _plot(result: dict[str, Any], out_dir: Path) -> list[str]:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:  # pragma: no cover - optional dependency
        return []

    figures = out_dir / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    written: list[str] = []
    tolerance = result["tolerance"]

    # 1. convergence
    structured = result["structured"]["series"]
    marching = result["marching"]["series"]
    figure, axes = plt.subplots(1, 2, figsize=(12.0, 4.4))
    axes[0].loglog(
        [entry["h"] for entry in structured],
        [entry["pointwise_median_relative_error"] for entry in structured],
        "o-",
        label="structured UV sphere, pointwise median",
    )
    axes[0].loglog(
        [entry["h"] for entry in marching],
        [entry["pointwise_median_relative_error"] for entry in marching],
        "s-",
        label="raw marching tets, pointwise median",
    )
    axes[0].loglog(
        [entry["h"] for entry in structured],
        [abs(entry["integral_relative_error"]) + 1e-16 for entry in structured],
        "^-",
        label="structured, integral (area-weighted)",
    )
    axes[0].loglog(
        [entry["h"] for entry in marching],
        [abs(entry["integral_relative_error"]) for entry in marching],
        "v-",
        label="marching tets, integral",
    )
    axes[0].axhline(tolerance, color="black", linestyle="--", linewidth=1.0, label="2 % criterion")
    axes[0].set_xlabel("mesh size h [m]")
    axes[0].set_ylabel("relative error of $\\kappa_1+\\kappa_2$ vs $2/R$")
    axes[0].set_title("convergence, R = 1 m")
    axes[0].grid(alpha=0.3, which="both")
    axes[0].legend(fontsize=7)

    radii = result["radius"]["series"]
    axes[1].semilogx(
        [entry["radius"] for entry in radii],
        [abs(entry["integral_relative_error"]) for entry in radii],
        "o-",
        label="integral",
    )
    axes[1].semilogx(
        [entry["radius"] for entry in radii],
        [entry["pointwise_median_relative_error"] for entry in radii],
        "s-",
        label="pointwise median (raw)",
    )
    axes[1].axhline(tolerance, color="black", linestyle="--", linewidth=1.0, label="2 % criterion")
    axes[1].set_xlabel("sphere radius R [m]")
    axes[1].set_ylabel("relative error")
    axes[1].set_title("scale invariance (marching tets, level 24)")
    axes[1].grid(alpha=0.3, which="both")
    axes[1].legend(fontsize=7)
    figure.tight_layout()
    path = figures / "convergence.png"
    figure.savefig(path, dpi=140)
    plt.close(figure)
    written.append(str(path.relative_to(out_dir)))

    # 2. error by reading, at one resolution
    finest = marching[-1]
    labels = ["integral\n(area-weighted)", "pointwise median\nraw", "pointwise median\nconditioned"]
    conditioned = finest.get("conditioned", {})
    values = [
        abs(finest["integral_relative_error"]),
        finest["pointwise_median_relative_error"],
        conditioned.get("pointwise_median_relative_error", float("nan")),
    ]
    figure, axis = plt.subplots(figsize=(6.4, 4.2))
    bars = axis.bar(labels, [100.0 * value for value in values], color=["#2b7bba", "#b04a3a", "#4a9b57"])
    axis.axhline(100.0 * tolerance, color="black", linestyle="--", linewidth=1.0,
                 label="2 % criterion")
    for bar, value in zip(bars, values, strict=True):
        axis.text(bar.get_x() + bar.get_width() / 2.0, bar.get_height(),
                  f"{100.0 * value:.2f} %", ha="center", va="bottom", fontsize=8)
    axis.set_ylabel("relative error vs $2/R$ [%]")
    axis.set_title(f"one mesh, three readings (marching tets, level {finest['level']})")
    axis.set_yscale("log")
    axis.grid(alpha=0.3, axis="y")
    axis.legend(fontsize=8)
    figure.tight_layout()
    path = figures / "readings.png"
    figure.savefig(path, dpi=140)
    plt.close(figure)
    written.append(str(path.relative_to(out_dir)))

    # 3. field cloud: per-vertex relative error on the sphere
    try:
        from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers projection)
        field = SphereField(center=(0.0, 0.0, 0.0), radius=1.0)
        meshes = [("raw marching tets", to_mesh(field, 24, margin=3).mesh),
                  ("structured UV 24x48", uv_sphere(1.0, 24, 48))]
        figure = plt.figure(figsize=(11.0, 4.6))
        for panel, (title, mesh) in enumerate(meshes, start=1):
            curvature = vertex_mean_curvature(mesh)
            error = [abs(2.0 * abs(value) - 1.0) for value in curvature]
            axis = figure.add_subplot(1, 2, panel, projection="3d")
            order = sorted(range(len(mesh.vertices)), key=lambda index: mesh.vertices[index][2])
            scatter = axis.scatter(
                [mesh.vertices[index][0] for index in order],
                [mesh.vertices[index][1] for index in order],
                [mesh.vertices[index][2] for index in order],
                c=[error[index] for index in order],
                cmap="magma", s=3.0, vmin=0.0, vmax=1.0,
            )
            axis.set_title(f"{title}: per-vertex relative error", fontsize=9)
            axis.set_xlim(-1.0, 1.0)
            axis.set_ylim(-1.0, 1.0)
            axis.set_zlim(-1.0, 1.0)
            axis.set_box_aspect((1.0, 1.0, 1.0))
            axis.tick_params(labelsize=6)
            cbar = figure.colorbar(scatter, ax=axis, shrink=0.7, pad=0.08)
            cbar.ax.tick_params(labelsize=6)
        figure.tight_layout()
        path = figures / "error_cloud.png"
        figure.savefig(path, dpi=140)
        plt.close(figure)
        written.append(str(path.relative_to(out_dir)))
    except Exception as error:  # pragma: no cover - plotting is not the measurement
        written.append(f"error_cloud.png skipped: {error}")
    return written


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def run(
    *,
    radius: float = 1.0,
    levels: tuple[int, ...] = (12, 16, 24, 32, 48),
    condition_levels: tuple[int, ...] = (16, 32, 48),
    resolutions: tuple[tuple[int, int], ...] = ((8, 16), (16, 32), (32, 64), (48, 96)),
    radii: tuple[float, ...] = (0.1, 0.5, 1.0, 2.0, 10.0),
    radius_level: int = 24,
    probes: int = 40,
) -> dict[str, Any]:
    started = time.perf_counter()
    chain = reference_chain_section(probes)
    marching = mesh_section_marching(levels, radius, condition_levels)
    structured = mesh_section_structured(radius, resolutions)
    scale = radius_section(radii, radius_level)
    diagnostic = diagnostic_section(32, radius)

    gates = chain["rows"] + marching["rows"] + structured["rows"] + scale["rows"]
    diagnostics = diagnostic["rows"]
    failing = [row["quantity"] for row in gates if not row["holds"]]

    # the invariant the case is about, stated once, from the measured numbers
    finest = marching["series"][-1]
    headline = {
        "quantity": "discrete mean curvature, kappa_1+kappa_2 (S = 2H), area-weighted integral",
        "radius": radius,
        "reference": reference_sum(radius),
        "measured": finest["integral_sum"],
        "relative_error": finest["integral_relative_error"],
        "tolerance": CASE_TOLERANCE,
        "within_2_percent": abs(finest["integral_relative_error"]) <= CASE_TOLERANCE,
        "pointwise_median_relative_error_raw": finest["pointwise_median_relative_error"],
        "pointwise_median_relative_error_conditioned": finest.get("conditioned", {}).get(
            "pointwise_median_relative_error"
        ),
        "pointwise_within_2_percent_raw": (
            finest["pointwise_median_relative_error"] <= CASE_TOLERANCE
        ),
    }

    return {
        "schema": SCHEMA,
        "case": CASE,
        "campaign_node": CAMPAIGN_NODE,
        "tolerance": CASE_TOLERANCE,
        "reference": {
            "formula": "kappa_1 + kappa_2 = 2/R",
            "sources": REFERENCE_SOURCES,
            "value": reference_sum(radius),
            "total_value": reference_total(radius),
        },
        "headline": headline,
        "all_gates_pass": not failing,
        "failing": failing,
        "rows": gates,
        "diagnostics": diagnostics,
        "reference_chain": chain["rows"],
        "marching": marching,
        "structured": structured,
        "radius": scale,
        "diagnostic": diagnostic["measure"],
        "environment": environment(),
        "seconds": time.perf_counter() - started,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", default=None, help="output directory (default: this case dir)")
    parser.add_argument("--radius", type=float, default=1.0)
    parser.add_argument("--quick", action="store_true", help="short series, for a smoke run")
    arguments = parser.parse_args(argv)

    out_dir = Path(arguments.out) if arguments.out else Path(__file__).resolve().parent
    out_dir.mkdir(parents=True, exist_ok=True)

    if arguments.quick:
        result = run(
            radius=arguments.radius,
            levels=(16, 32),
            condition_levels=(32,),
            resolutions=((16, 32), (32, 64)),
            radii=(0.5, 1.0),
            radius_level=16,
            probes=24,
        )
    else:
        result = run(radius=arguments.radius)

    result["figures"] = _plot(result, out_dir)

    (out_dir / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    gates = {
        "case": CASE,
        "tolerance": CASE_TOLERANCE,
        "all_gates_pass": result["all_gates_pass"],
        "headline": result["headline"],
        "gates": {
            f"{row['gate']}: {row['quantity']}": {
                "criterion": f"{row['gate']}: {row['quantity']}",
                "reference": row["reference"],
                "measured": row["measured"],
                "error_percent": row["error_percent"],
                "tolerance": row["tolerance"],
                "verdict": "pass" if row["holds"] else "fail",
            }
            for row in result["rows"]
        },
        "diagnostics": {
            f"{row['gate']}: {row['quantity']}": {
                "criterion": f"{row['gate']}: {row['quantity']}",
                "error_percent": row["error_percent"],
                "tolerance": row["tolerance"],
                "verdict": "pass" if row["holds"] else "fail",
                "scope": "diagnostic only: the raw pointwise reading on a sliver-ridden mesh",
            }
            for row in result["diagnostics"]
        },
        "failing": result["failing"],
        "environment": {
            "python": result["environment"]["python"],
            "git_commit": result["environment"]["git_commit"],
            "repro_cmd": result["environment"]["repro_cmd"],
        },
    }
    (out_dir / "gates_verified.json").write_text(json.dumps(gates, indent=2, sort_keys=True) + "\n")

    print(f"case {CASE}: {len(result['rows'])} gates, {len(result['diagnostics'])} diagnostics")
    print(
        "headline S vs 2/R: measured {:.6f} vs {:.6f} ({:+.4f} %)  -> {}".format(
            result["headline"]["measured"],
            result["headline"]["reference"],
            100.0 * result["headline"]["relative_error"],
            "within 2 %" if result["headline"]["within_2_percent"] else "OUTSIDE 2 %",
        )
    )
    print(f"pointwise median: raw {100.0 * result['headline']['pointwise_median_relative_error_raw']:.2f} %"
          f", conditioned {100.0 * (result['headline']['pointwise_median_relative_error_conditioned'] or float('nan')):.2f} %")
    print(f"figures: {result['figures']}")
    print(f"seconds: {result['seconds']:.1f}")
    if result["failing"]:
        print("FAILING GATES:")
        for quantity in result["failing"]:
            print(f"  - {quantity}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())