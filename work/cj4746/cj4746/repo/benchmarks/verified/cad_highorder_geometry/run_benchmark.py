#!/usr/bin/env python3
"""Verified benchmark: TensorCAD as a p-order geometry truth source.

Campaign node ``p127w2-n01-cad-highorder-geometry``. The script recomputes every number
that ``docs/highorder-geometry.md`` quotes and writes

* ``result.json``          -- measurements, environment, provenance and the raw rows,
* ``gates_verified.json``  -- J1..J7 with criterion, measured value and verdict,
* ``figures/*.png``        -- convergence, curvature error, cross-validation/sensitivity,
* ``artifacts/contract_*.npz|.h5|.json`` -- real element contracts plus their hashes.

Run from anywhere::

    python benchmarks/verified/cad_highorder_geometry/run_benchmark.py

Exit code 0 means every gate passed. A failing gate still writes its artifacts: an honest
FAIL is evidence, and a benchmark that hides its failures is worthless.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import platform
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

CASEDIR = Path(__file__).resolve().parent
ROOT = CASEDIR.parents[2]
sys.path.insert(0, str(ROOT))

import numpy as np

from tensorcad.geometry_field import EllipsoidField
from tensorcad.highorder import (
    contract,
    elliptic,
    quadrature,
)
from tensorcad.highorder import nurbs_geom as geom
from tensorcad.nurbs import NURBSSurface, sphere_surface
from tensorcad.voxel import to_voxel_grid

CASE = "cad_highorder_geometry"
SCHEMA = "tensorcad.highorder.benchmark/1.0"
GATES_SCHEMA = "tensorcad.gates/1.0"

#: The hull used for the hydrostatic cross-check: distinct semi-axes so that no symmetry
#: can hide a defect, and a draft/length ratio of the order of a real hull's.
ELLIPSOID = (1.0, 0.6, 0.42)
SPHERE_RADIUS = 1.0
CURVATURE_TOLERANCE = 1e-3
ROUTE_TOLERANCE = 5e-3


# ----------------------------------------------------------------------------------
# deterministic demo geometry
# ----------------------------------------------------------------------------------
def bezier_patch(degree: int) -> NURBSSurface:
    """A curved, non-rational, single-span Bezier patch of the requested degree."""
    p = degree
    net = []
    for i in range(p + 1):
        s = i / p if p else 0.0
        row = []
        for j in range(p + 1):
            t = j / p if p else 0.0
            z = 0.35 * math.sin(math.pi * s) * math.sin(math.pi * t) + 0.15 * s * t
            row.append((s, 0.7 * t, z))
        net.append(tuple(row))
    knots = (0.0,) * (p + 1) + (1.0,) * (p + 1)
    return NURBSSurface(degree=(p, p), knots=(knots, knots), control_net=tuple(net))


def patch_bounds(surface: NURBSSurface, samples: int = 65) -> dict[str, float]:
    """Axis-aligned bounds of the patch on a fixed parameter grid."""
    (u_low, u_high), (v_low, v_high) = surface.domain
    xs: list[float] = []
    ys: list[float] = []
    zs: list[float] = []
    for i in range(samples):
        u = u_low + (u_high - u_low) * i / (samples - 1)
        for j in range(samples):
            v = v_low + (v_high - v_low) * j / (samples - 1)
            x, y, z = surface.evaluate(u, v)
            xs.append(x)
            ys.append(y)
            zs.append(z)
    return {"L": max(xs) - min(xs), "B": max(ys) - min(ys), "T": max(zs) - min(zs)}


def curvature_and_wetted_area(axes: tuple[float, float, float]) -> tuple[float, float]:
    """Max principal curvature (closed form) and wetted surface (Legendre) for an ellipsoid."""
    surface = geom.ellipsoid_surface(*axes)
    kappa = 0.0
    for point in geom.curvature_scan(surface, grid=4, u_margin=0.05, v_margin=0.05):
        kappa = max(kappa, max(geom.ellipsoid_principal_curvatures(*axes, point.point)))
    return kappa, 0.5 * elliptic.triaxial_ellipsoid_area(*axes)


# ----------------------------------------------------------------------------------
# gates
# ----------------------------------------------------------------------------------
def gate_j1() -> dict[str, Any]:
    """Gauss exactness order 2p-1 plus curved-patch quadrature order >= 2p-1, p = 1..4."""
    rows = []
    passed = True
    for degree in (1, 2, 3, 4):
        surface = bezier_patch(degree)
        area = quadrature.surface_area(surface, order=32)
        study = quadrature.convergence_study(
            surface,
            lambda u, v, s=surface: quadrature.area_element(s, u, v),
            orders=(2, 3, 4, 6, 8),
            reference_order=32,
            exact=area,
        )
        probe = quadrature.gauss_exactness_probe(degree)
        target = 2 * degree - 1
        observed = study["fitted_order"]
        gate = bool(
            probe["exact_at_limit"]
            and probe["fails_beyond_limit"]
            and observed is not None
            and observed >= target
        )
        passed = passed and gate
        rows.append(
            {
                "degree": degree,
                "target_order": target,
                "exactness_degree": probe["degree_exact"],
                "exactness_error_at_limit": probe["error_at_degree_exact"],
                "exactness_error_beyond_limit": probe["error_one_degree_beyond"],
                "observed_order_fitted": observed,
                "pairwise_orders": study["pairwise_orders"],
                "gauss_orders": [row["order"] for row in study["rows"]],
                "relative_errors": [row["relative_error"] for row in study["rows"]],
                "passed": gate,
            }
        )
    return {
        "criterion": "exactness order = 2p-1 (exact there, wrong one degree beyond) and fitted "
        "quadrature order >= 2p-1 for p = 1..4",
        "rows": rows,
        "worst_observed_order": max(row["observed_order_fitted"] or 0.0 for row in rows),
        "tightest_margin": min(
            (row["observed_order_fitted"] or 0.0) - row["target_order"] for row in rows
        ),
        "passed": bool(passed),
    }


def gate_j2() -> dict[str, Any]:
    """Normals and principal curvatures vs the closed form and vs finite differences."""
    a, b, c = ELLIPSOID
    surface = geom.ellipsoid_surface(a, b, c)
    sphere = sphere_surface(1.5)
    ellipsoid_points = geom.curvature_scan(surface, grid=7)
    sphere_points = geom.curvature_scan(sphere, grid=5)

    kappa_max = max(
        max(geom.ellipsoid_principal_curvatures(a, b, c, point.point)) for point in ellipsoid_points
    )
    worst_curvature = 0.0
    worst_normal = 0.0
    worst_fd = 0.0
    scans: list[dict[str, Any]] = []
    for point in ellipsoid_points:
        closed = sorted(geom.ellipsoid_principal_curvatures(a, b, c, point.point), reverse=True)
        measured = sorted((abs(value) for value in point.principal_curvatures), reverse=True)
        curvature_error = max(abs(l - r) for l, r in zip(measured, closed, strict=True)) / kappa_max
        exact_normal = geom.ellipsoid_normal(a, b, c, point.point)
        normal_error = math.sqrt(
            sum((l - r) ** 2 for l, r in zip(point.normal, exact_normal, strict=True))
        )
        difference = geom.finite_difference_geometry(surface, point.u, point.v)
        fd_measured = sorted((abs(value) for value in difference.principal_curvatures), reverse=True)
        fd_error = max(abs(l - r) for l, r in zip(fd_measured, closed, strict=True)) / kappa_max
        worst_curvature = max(worst_curvature, curvature_error)
        worst_normal = max(worst_normal, normal_error)
        worst_fd = max(worst_fd, fd_error)
        scans.append(
            {
                "u": point.u,
                "v": point.v,
                "kappa_max_closed": closed[0],
                "kappa_min_closed": closed[1],
                "curvature_error": curvature_error,
                "normal_error": normal_error,
                "finite_difference_error": fd_error,
            }
        )
    sphere_kappa = 1.0 / 1.5
    sphere_curvature = max(
        abs(abs(value) - sphere_kappa) / sphere_kappa
        for point in sphere_points
        for value in point.principal_curvatures
    )
    sphere_normal = max(
        math.sqrt(
            sum(
                (l - r / math.hypot(*point.point)) ** 2
                for l, r in zip(point.normal, point.point, strict=True)
            )
        )
        for point in sphere_points
    )
    passed = (
        worst_curvature <= CURVATURE_TOLERANCE
        and worst_normal <= CURVATURE_TOLERANCE
        and worst_fd <= CURVATURE_TOLERANCE
        and sphere_curvature <= CURVATURE_TOLERANCE
        and sphere_normal <= CURVATURE_TOLERANCE
    )
    return {
        "criterion": "L_inf curvature error / kappa_max <= 1e-3 and unit-normal deviation <= 1e-3, "
        "against the closed-form ellipsoid/sphere and against a finite-difference route",
        "ellipsoid": {
            "semi_axes": [a, b, c],
            "samples": len(ellipsoid_points),
            "kappa_max": kappa_max,
            "worst_curvature_error": worst_curvature,
            "worst_normal_error": worst_normal,
            "worst_finite_difference_error": worst_fd,
        },
        "sphere": {
            "radius": 1.5,
            "samples": len(sphere_points),
            "worst_curvature_error": sphere_curvature,
            "worst_normal_error": sphere_normal,
        },
        "scans": scans,
        "passed": bool(passed),
    }


def gate_j3() -> dict[str, Any]:
    """Volume of the isoparametric element set vs analytic, relative error <= 1e-6."""
    sphere = sphere_surface(SPHERE_RADIUS)
    ellipsoid = geom.ellipsoid_surface(*ELLIPSOID)
    cases = []
    passed = True
    for name, surface, exact in (
        ("sphere_r1", sphere, 4.0 / 3.0 * math.pi * SPHERE_RADIUS**3),
        (
            "ellipsoid_1_0.6_0.42",
            ellipsoid,
            4.0 / 3.0 * math.pi * ELLIPSOID[0] * ELLIPSOID[1] * ELLIPSOID[2],
        ),
    ):
        exact_area = quadrature.surface_area(surface, order=16)
        rows = []
        for order in (2, 3, 4, 6):
            elements = quadrature.curved_elements(surface, order=order, subdivisions=(2, 2))
            volume = math.fsum(element.volume_contribution for element in elements)
            area = math.fsum(element.area for element in elements)
            rows.append(
                {
                    "order": order,
                    "elements": len(elements),
                    "volume": volume,
                    "volume_relative_error": abs(volume - exact) / exact,
                    "area_relative_error": abs(area - exact_area) / exact_area,
                }
            )
        order4 = next(row for row in rows if row["order"] == 4)
        gate = order4["volume_relative_error"] <= 1e-6
        passed = passed and gate
        cases.append(
            {
                "case": name,
                "exact_volume": exact,
                "exact_area": exact_area,
                "degree": list(surface.degree),
                "subdivisions": [2, 2],
                "rows": rows,
                "order4_volume_relative_error": order4["volume_relative_error"],
                "passed": gate,
            }
        )
    area_convergence = []
    exact_sphere_area = 4.0 * math.pi * SPHERE_RADIUS**2
    for order in (2, 3, 4, 6, 8):
        value = quadrature.surface_area(sphere, order=order)
        area_convergence.append(
            {
                "order": order,
                "area": value,
                "relative_error": abs(value - exact_sphere_area) / exact_sphere_area,
            }
        )
    return {
        "criterion": "order-4 curved elements on 2x2 subdivisions per knot span: volume relative "
        "error <= 1e-6 against 4/3 pi R^3 and 4/3 pi a b c",
        "cases": cases,
        "sphere_area_convergence": area_convergence,
        "passed": bool(passed),
    }


def gate_j4() -> dict[str, Any]:
    """Hull Cb and wetted surface: NURBS / STL / voxel vs the AGM-elliptic truth."""
    a, b, c = ELLIPSOID
    hull = geom.half_ellipsoid_patch(a, b, c)
    reference = geom.hydrostatics_reference(a, b, c)

    # Geometry gate first: a wrong L/B/T/Cb makes every downstream number meaningless.
    bounds = patch_bounds(hull)
    geometry_gate: dict[str, Any] = {
        key: {
            "measured": bounds[key],
            "expected": expected,
            "relative_error": abs(bounds[key] - expected) / expected,
        }
        for key, expected in (("L", 2 * a), ("B", 2 * b), ("T", c))
    }
    geometry_gate["Cb"] = {
        "measured": reference["Cb"],
        "expected": math.pi / 6.0,
        "relative_error": abs(reference["Cb"] - math.pi / 6.0) / (math.pi / 6.0),
    }
    if any(entry["relative_error"] >= 0.01 for entry in geometry_gate.values()):
        return {
            "criterion": "geometry gate Cb/L/B/T mutual check < 1% before any hydrostatic claim",
            "geometry_gate": geometry_gate,
            "refused": "geometry gate failed; hydrostatic claims refused",
            "passed": False,
        }

    nurbs_area = quadrature.surface_area(hull, order=16)
    nurbs_volume = quadrature.surface_volume(hull, order=16)
    nurbs_cb = nurbs_volume / (2 * a * 2 * b * c)

    stl_rows = []
    for samples in (32, 64, 96, 128):
        area, _ = quadrature.faceted_area_and_volume(hull, samples)
        stl_rows.append(
            {
                "samples": samples,
                "facets": len(hull.tessellate(samples, samples).mesh.faces),
                "wetted_surface": area,
                "wetted_surface_relative_error": abs(area - reference["wetted_surface"])
                / reference["wetted_surface"],
            }
        )

    voxel_rows = []
    for resolution in (48, 64, 96, 128):
        grid = to_voxel_grid(
            EllipsoidField((a, b, c), (0.0, 0.0, 0.0)), resolution=resolution, padding=2
        )
        half_volume = 0.5 * grid.volume()  # the grid is symmetric about z = 0
        cb = half_volume / (2 * a * 2 * b * c)
        cb_hard = 0.5 * grid.occupied_volume() / (2 * a * 2 * b * c)
        voxel_rows.append(
            {
                "resolution": resolution,
                "voxel_volume_full": grid.volume(),
                "half_hull_volume": half_volume,
                "Cb": cb,
                "Cb_relative_error": abs(cb - reference["Cb"]) / reference["Cb"],
                "Cb_hard_count": cb_hard,
                "Cb_hard_count_relative_error": abs(cb_hard - reference["Cb"]) / reference["Cb"],
            }
        )

    base_kappa, base_area = curvature_and_wetted_area((a, b, c))
    sensitivity = []
    for label, axes in (
        ("isotropic +1%", (1.01 * a, 1.01 * b, 1.01 * c)),
        ("isotropic -1%", (0.99 * a, 0.99 * b, 0.99 * c)),
        ("beam +1%", (a, 1.01 * b, c)),
        ("beam -1%", (a, 0.99 * b, c)),
        ("draft +1%", (a, b, 1.01 * c)),
        ("draft -1%", (a, b, 0.99 * c)),
    ):
        kappa, area = curvature_and_wetted_area(axes)
        dk = abs(kappa - base_kappa) / base_kappa
        ds = abs(area - base_area) / base_area
        sensitivity.append(
            {
                "perturbation": label,
                "curvature_relative_change": dk,
                "wetted_surface_relative_change": ds,
                "ratio_area_to_curvature": (ds / dk) if dk > 0.0 else None,
            }
        )

    kappas = (0.0, 0.2, 0.5, 0.8, 0.95)
    passed = bool(
        abs(nurbs_area - reference["wetted_surface"]) / reference["wetted_surface"]
        <= ROUTE_TOLERANCE
        and abs(nurbs_cb - reference["Cb"]) / reference["Cb"] <= ROUTE_TOLERANCE
        and all(row["wetted_surface_relative_error"] <= ROUTE_TOLERANCE for row in stl_rows)
        and all(row["Cb_relative_error"] <= ROUTE_TOLERANCE for row in voxel_rows)
    )
    return {
        "criterion": "wetted surface S and block coefficient Cb from NURBS quadrature, STL "
        "tessellation and voxel grid all within 0.5% of the independent AGM/Carlson "
        "elliptic-integral value",
        "geometry_gate": geometry_gate,
        "reference": reference,
        "nurbs": {
            "wetted_surface": nurbs_area,
            "wetted_surface_relative_error": abs(nurbs_area - reference["wetted_surface"])
            / reference["wetted_surface"],
            "volume": nurbs_volume,
            "volume_relative_error": abs(nurbs_volume - reference["volume"]) / reference["volume"],
            "Cb": nurbs_cb,
            "Cb_relative_error": abs(nurbs_cb - reference["Cb"]) / reference["Cb"],
        },
        "stl": stl_rows,
        "voxel": voxel_rows,
        "sensitivity": sensitivity,
        "elliptic_routes": {
            "triaxial_legendre": elliptic.triaxial_ellipsoid_area(a, b, c),
            "ellipe_vs_agm_max_difference": max(
                abs(elliptic.ellipe(k) - elliptic.ellipe_agm(k)) for k in kappas
            ),
            "ellipk_vs_agm_max_difference": max(
                abs(elliptic.ellipk(k) - elliptic.ellipk_agm(k)) for k in kappas
            ),
        },
        "passed": passed,
    }


def gate_j5() -> dict[str, Any]:
    """Contract write / read / recompute in three formats: bit identity and hash stability."""
    artifacts = CASEDIR / "artifacts"
    artifacts.mkdir(exist_ok=True)
    reports = []
    passed = True
    for name, surface in (
        ("sphere_r1", sphere_surface(SPHERE_RADIUS)),
        ("ellipsoid_1_0.6_0.42", geom.ellipsoid_surface(*ELLIPSOID)),
    ):
        elements = quadrature.curved_elements(surface, order=4, subdivisions=(2, 2))
        provenance = {
            "case": name,
            "generator": "benchmarks/verified/cad_highorder_geometry/run_benchmark.py",
            "degree": list(surface.degree),
            "quadrature_order": 4,
            "subdivisions": [2, 2],
        }
        payload = contract.build_payload(elements, provenance=provenance)
        formats = [".npz", ".json"]
        if importlib.util.find_spec("h5py") is not None:
            formats.insert(1, ".h5")
        else:
            reports.append(
                {
                    "format": ".h5",
                    "case": name,
                    "skipped": "h5py is not installed on this host; the writer is optional",
                    "passed": None,
                }
            )
        for suffix in formats:
            target = artifacts / f"contract_{name}{suffix}"
            contract.write_contract(target, elements, provenance=provenance)
            report = contract.verify_roundtrip(target, payload)
            report["format"] = suffix
            report["case"] = name
            report["bytes"] = target.stat().st_size
            reports.append(report)
            passed = bool(
                passed
                and report["all_bit_identical"]
                and report["hash_match"]
                and report["recomputed_max_abs_difference"] <= 1e-14
            )
        primary = artifacts / f"contract_{name}.npz"
        reports.append(
            {
                "format": "npz_as_written",
                "case": name,
                "rows": payload["header"]["rows"],
                "elements": payload["header"]["elements"],
                "content_hash": payload["header"]["content_hash"],
                "sha256_file": hashlib.sha256(primary.read_bytes()).hexdigest(),
                "recomputed": contract.recompute_integrals(payload["arrays"]),
            }
        )
    return {
        "criterion": "write-read-recompute is bit identical (recomputed integrals agree to "
        "<= 1e-14) and the content hash survives NPZ, HDF5 and JSON",
        "reports": reports,
        "worst_field_difference": max(
            entry["max_abs_difference"]
            for report in reports
            if "fields" in report
            for entry in report["fields"].values()
        ),
        "passed": bool(passed),
    }


def gate_j6() -> dict[str, Any]:
    """Negative control: the faceted patch must lose the p-order and miss the volume gate."""
    surface = sphere_surface(SPHERE_RADIUS)
    study = quadrature.tessellation_convergence(
        surface,
        (8, 16, 32, 64),
        exact_area=4.0 * math.pi * SPHERE_RADIUS**2,
        exact_volume=4.0 / 3.0 * math.pi * SPHERE_RADIUS**3,
    )
    area_order = study["area_fitted_order"]
    volume_order = study["volume_fitted_order"]
    coarse = [row["volume_relative_error"] for row in study["rows"][:2]]
    passed = bool(
        area_order is not None
        and area_order <= 2.05
        and volume_order is not None
        and volume_order <= 2.05
        and area_order < 3.0
        and all(error > 1e-3 for error in coarse)
    )
    return {
        "criterion": "faceted patch: observed order <= 2 for area and volume, and volume relative "
        "error > 1e-3 at 8 and 16 facets per direction (i.e. J1/J3 cannot be met by facets)",
        "rows": study["rows"],
        "area_fitted_order": area_order,
        "volume_fitted_order": volume_order,
        "area_pairwise_orders": study["area_pairwise_orders"],
        "volume_pairwise_orders": study["volume_pairwise_orders"],
        "passed": passed,
    }


def gate_j7() -> dict[str, Any]:
    """The gate module must run inside the 120 s CI budget (measured in a subprocess)."""
    test_file = ROOT / "tests" / "test_highorder_geometry.py"
    started = time.perf_counter()
    completed = subprocess.run(
        [sys.executable, "-m", "pytest", str(test_file), "-q"],
        capture_output=True,
        text=True,
        cwd=str(ROOT),
        check=False,
    )
    elapsed = time.perf_counter() - started
    tail = [line for line in completed.stdout.strip().splitlines() if line.strip()][-1:]
    return {
        "criterion": "pytest tests/test_highorder_geometry.py wall time <= 120 s and all green",
        "elapsed_seconds": elapsed,
        "returncode": completed.returncode,
        "summary": tail[0] if tail else "",
        "passed": bool(completed.returncode == 0 and elapsed <= 120.0),
    }


# ----------------------------------------------------------------------------------
# figures
# ----------------------------------------------------------------------------------
def make_figures(gates: dict[str, Any]) -> list[dict[str, Any]]:
    """Draw the report figures; a host without matplotlib keeps its gate evidence.

    The gates are the machine criterion and they do not need a plotting stack, so a missing
    matplotlib degrades to a recorded skip instead of failing a benchmark whose numbers are
    already computed and written.
    """
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return [
            {
                "name": "none",
                "path": "",
                "caption": "matplotlib is not installed on this host; the figures were not drawn. "
                "The gate evidence in result.json is complete and independent of plotting.",
            }
        ]

    figures_dir = CASEDIR / "figures"
    figures_dir.mkdir(exist_ok=True)
    figures: list[dict[str, Any]] = []

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    for row in gates["J1"]["rows"]:
        pairs = [
            (order, error)
            for order, error in zip(row["gauss_orders"], row["relative_errors"], strict=True)
            if error > 1e-14
        ]
        axes[0].loglog(
            [pair[0] for pair in pairs],
            [pair[1] for pair in pairs],
            "o-",
            label=f"p = {row['degree']} (fitted {row['observed_order_fitted']:.1f})",
        )
    axes[0].set_xlabel("Gauss points per direction")
    axes[0].set_ylabel("relative area error")
    axes[0].set_title("NURBS patch: spectral convergence")
    axes[0].grid(True, which="both", alpha=0.3)
    axes[0].legend(fontsize=8)
    faceted = gates["J6"]["rows"]
    axes[1].loglog(
        [row["samples"] for row in faceted],
        [row["area_relative_error"] for row in faceted],
        "s-",
        label="area (faceted)",
    )
    axes[1].loglog(
        [row["samples"] for row in faceted],
        [row["volume_relative_error"] for row in faceted],
        "^-",
        label="volume (faceted)",
    )
    guide = [row["samples"] for row in faceted]
    reference_error = faceted[0]["area_relative_error"]
    axes[1].loglog(
        guide,
        [reference_error * (guide[0] / value) ** 2 for value in guide],
        "k--",
        alpha=0.6,
        label="slope -2 guide",
    )
    axes[1].axhline(1e-3, color="crimson", linestyle=":", label="J3 volume gate 1e-3")
    axes[1].set_xlabel("facets per direction")
    axes[1].set_ylabel("relative error")
    axes[1].set_title("Negative control: linearized patch is second order")
    axes[1].grid(True, which="both", alpha=0.3)
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    path = figures_dir / "convergence.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    figures.append(
        {
            "name": "convergence.png",
            "path": str(path),
            "caption": "Left: Gauss-Legendre area error against points per direction for "
            "degree-1..4 Bezier patches (spectral decay; fitted order in the legend). Right: the "
            "same patch linearized -- second order only, and above the 1e-3 volume gate at every "
            "affordable facet count.",
        }
    )

    scans = gates["J2"]["scans"]
    grid = round(math.sqrt(len(scans)))
    curvature_map = np.zeros((grid, grid))
    normal_map = np.zeros((grid, grid))
    for index, entry in enumerate(scans):
        i, j = divmod(index, grid)
        curvature_map[i, j] = entry["curvature_error"]
        normal_map[i, j] = entry["normal_error"]
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.0))
    for axis, data, title in (
        (axes[0], curvature_map, "|kappa| - kappa_closed, over kappa_max"),
        (axes[1], normal_map, "|n - n_closed|"),
    ):
        image = axis.imshow(data, origin="lower", cmap="viridis", aspect="auto")
        axis.set_title(f"{title}\nmax {data.max():.2e} (gate 1e-3)", fontsize=9)
        axis.set_xlabel("v index")
        axis.set_ylabel("u index")
        fig.colorbar(image, ax=axis, fraction=0.046, pad=0.02)
    fig.suptitle("Ellipsoid 1.0 x 0.6 x 0.42: NURBS differential geometry vs closed form", fontsize=10)
    fig.tight_layout()
    path = figures_dir / "curvature_error.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    figures.append(
        {
            "name": "curvature_error.png",
            "path": str(path),
            "caption": "Curvature and unit-normal error over the parameter square of the "
            "ellipsoid, both against the closed-form implicit-surface value; several orders of "
            "magnitude below the 1e-3 gate away from the poles.",
        }
    )

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    stl = gates["J4"]["stl"]
    axes[0].semilogy(
        [row["samples"] for row in stl],
        [row["wetted_surface_relative_error"] for row in stl],
        "s-",
        label="STL facets: wetted surface",
    )
    voxel = gates["J4"]["voxel"]
    axes[0].semilogy(
        [row["resolution"] for row in voxel],
        [row["Cb_relative_error"] for row in voxel],
        "^-",
        label="voxel grid: Cb",
    )
    axes[0].axhline(ROUTE_TOLERANCE, color="crimson", linestyle=":", label="J4 gate 0.5%")
    axes[0].axhline(
        gates["J4"]["nurbs"]["wetted_surface_relative_error"],
        color="black",
        linestyle="--",
        label="NURBS quadrature (round-off)",
    )
    axes[0].set_xlabel("discretization (facets or voxels per direction)")
    axes[0].set_ylabel("relative error vs elliptic-integral truth")
    axes[0].set_title("Three routes to S and Cb")
    axes[0].grid(True, which="both", alpha=0.3)
    axes[0].legend(fontsize=8)
    sensitivity = gates["J4"]["sensitivity"]
    axes[1].plot(
        [100.0 * row["curvature_relative_change"] for row in sensitivity],
        [100.0 * row["wetted_surface_relative_change"] for row in sensitivity],
        "o",
    )
    for row in sensitivity:
        axes[1].annotate(
            row["perturbation"],
            (100.0 * row["curvature_relative_change"], 100.0 * row["wetted_surface_relative_change"]),
            fontsize=7,
            xytext=(3, 3),
            textcoords="offset points",
        )
    top = max(100.0 * row["wetted_surface_relative_change"] for row in sensitivity) * 1.2
    axes[1].plot([0.0, top], [0.0, 2.0 * top], "k--", alpha=0.5, label="isotropic scaling: 2:1")
    axes[1].set_xlabel("max principal curvature change [%]")
    axes[1].set_ylabel("wetted surface change [%]")
    axes[1].set_title("Sensitivity of wetted area to curvature error")
    axes[1].grid(True, alpha=0.3)
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    path = figures_dir / "cross_validation.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    figures.append(
        {
            "name": "cross_validation.png",
            "path": str(path),
            "caption": "Left: STL wetted surface and voxel Cb against the AGM/Carlson elliptic-"
            "integral truth (the NURBS route sits at round-off). Right: wetted-area change against "
            "curvature change for isotropic and shape perturbations of the same hull -- the ratio "
            "is metric dependent, so a curvature guarantee does not replace an area gate.",
        }
    )
    return figures


# ----------------------------------------------------------------------------------
# orchestration
# ----------------------------------------------------------------------------------
def git_revision() -> dict[str, str]:
    def run(args: list[str]) -> str:
        try:
            completed = subprocess.run(
                args, capture_output=True, text=True, cwd=str(ROOT), check=False
            )
            return completed.stdout.strip()
        except OSError:  # pragma: no cover - git is present in every CI job here
            return ""

    # Only the *source* is reported dirty: the benchmark rewrites its own evidence files on
    # every run (timestamps), so a whole-tree check would always say "dirty" and mean nothing.
    return {
        "commit": run(["git", "rev-parse", "HEAD"]),
        "branch": run(["git", "rev-parse", "--abbrev-ref", "HEAD"]),
        "sources_dirty": "yes"
        if run(["git", "status", "--porcelain", "--", "tensorcad", "tests", "pyproject.toml"])
        else "no",
    }


def module_digest() -> str:
    digest = hashlib.sha256()
    for path in sorted((ROOT / "tensorcad" / "highorder").glob("*.py")):
        digest.update(path.name.encode("ascii"))
        digest.update(path.read_bytes())
    return digest.hexdigest()


def main() -> int:
    started = time.perf_counter()
    gates = {
        "J1": gate_j1(),
        "J2": gate_j2(),
        "J3": gate_j3(),
        "J4": gate_j4(),
        "J5": gate_j5(),
        "J6": gate_j6(),
        "J7": gate_j7(),
    }
    figures = make_figures(gates)
    elapsed = time.perf_counter() - started

    result: dict[str, Any] = {
        "schema": SCHEMA,
        "case": CASE,
        "campaign_node": "p127w2-n01-cad-highorder-geometry",
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "machine": platform.machine(),
            "numpy": np.__version__,
            "processor": platform.processor(),
        },
        "provenance": {
            **git_revision(),
            "module_sha256": module_digest(),
            "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "test_file_sha256": hashlib.sha256(
                (ROOT / "tests" / "test_highorder_geometry.py").read_bytes()
            ).hexdigest(),
        },
        "geometry": {
            "sphere_radius": SPHERE_RADIUS,
            "ellipsoid_semi_axes": list(ELLIPSOID),
            "curvature_tolerance": CURVATURE_TOLERANCE,
            "route_tolerance": ROUTE_TOLERANCE,
        },
        "gates": gates,
        "figures": figures,
        "elapsed_seconds": elapsed,
        "passed": all(gates[key]["passed"] for key in gates),
    }
    (CASEDIR / "result.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    (CASEDIR / "gates_verified.json").write_text(
        json.dumps(
            {
                "schema": GATES_SCHEMA,
                "case": CASE,
                "campaign_node": result["campaign_node"],
                "generated_utc": result["generated_utc"],
                "gates": {
                    key: {
                        "criterion": gates[key]["criterion"],
                        "passed": gates[key]["passed"],
                        "evidence": "result.json#/gates/" + key,
                    }
                    for key in gates
                },
                "passed": result["passed"],
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    for key in gates:
        mark = "PASS" if gates[key]["passed"] else "FAIL"
        print(f"{key}: {mark}")
    print(f"figures: {len(figures)}   elapsed: {elapsed:.1f} s")
    print(f"overall: {'PASS' if result['passed'] else 'FAIL'}")
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
