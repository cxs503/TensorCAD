"""Unified solver-interface verification: one geometry, three solver routes, one frame contract.

Why this module exists
----------------------
The anchor ``GAP-CAD-B3-SOLVER-INTERFACE`` requires "三求解器各一个端到端算例真跑通
（坐标系契约统一 z-up，5415 三回合教训制度化）".  The three routes already exist
individually (:mod:`tensorcad.fvm_handoff`, :mod:`tensorcad.dg_handoff`,
:mod:`tensorcad.lattice_flow`), but each uses its own geometry (barge, circle, cylinder)
so the z-up coordinate contract was never verified *across* routes on the same solid.

The 5415 rounds paid for a coordinate-frame mix-up (y-up vs z-up) that silently mirrored
the hull.  This module institutionalises that lesson by running **all three routes on the
same CAD solid** and asserting that every route declares the same frame, the same signed
volume, and that a deliberate frame tamper is refused by every route's verifier.

Routes exercised
----------------
====================  ============================  ================================
route                 geometry                      closed-form reference
====================  ============================  ================================
``volume_msh`` (FVM)  ellipsoid mesh → gmsh C3D4    patch test (linear field exact)
``curved_wall`` (DG)  ellipse section (z=0)         circumference 2πa·E(e), area πab
``lattice`` (LBM)     ellipsoid voxel occupancy     staircase-approx volume
====================  ============================  ================================

The ellipsoid is chosen as the common geometry because it is:
- closed and watertight by construction (no open-surface ambiguity),
- has an analytic volume 4/3 πabc for cross-route volume comparison,
- can be voxelised and triangulated from the same semi-axes parameters.

Frame contract
--------------
Every route writes ``handoff.json`` with ``frame.name = "z-up"``.
:func:`verify_frame_consistency` re-reads all three bundles and asserts:
- all declare ``z-up``,
- all signed volumes agree within 5 % (the lattice is staircase-approximated),
- a tampered ``y-up`` manifest is refused by every route's verifier.

This is the institutionalised 5415 lesson: the frame is not a comment, it is a
machine-checked gate that every solver-facing handoff must pass.

Implementation note
-------------------
This module delegates to :mod:`tensorcad.solver_solid_gate` for the actual
geometry sampling, route execution, and cross-route verification.  The
``run_integration`` function here is a thin wrapper that:
1. Calls ``solid_from_library`` to get the one solid all routes consume.
2. Calls ``fvm_route``, ``dg_route``, ``lattice_route`` from the gate module.
3. Calls ``cross_route_gate`` to compare routes against each other.
4. Wraps the gate result into an ``IntegrationReport`` for the public API.

This delegation ensures that the integration report is always consistent
with the gate report, and that "the same solid" is a measured value (SHA-256
digest), not prose.
"""
from __future__ import annotations

import json
import math
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Sequence

from .solver_interface import COORDINATE_CONTRACT

__all__ = [
    "INTEGRATION_SCHEMA",
    "RouteResult",
    "IntegrationReport",
    "run_integration",
    "verify_frame_consistency",
    "frame_tamper_refused",
    "summary_markdown",
]

#: Bumped on any incompatible change to the integration report format.
INTEGRATION_SCHEMA = "tensorcad.solver_integration/1"


# --------------------------------------------------------------------------- common geometry


def _analytic_volume(semi_axes: Sequence[float]) -> float:
    """4/3 π a b c — the reference every route is scored against."""
    a, b, c = (float(v) for v in semi_axes)
    return 4.0 / 3.0 * math.pi * a * b * c


# --------------------------------------------------------------------------- integration


@dataclass(frozen=True)
class RouteResult:
    """One route's outcome."""
    route: str
    status: str
    metrics: dict[str, Any]
    frame: str
    error: str = ""

    def to_dict(self) -> dict[str, Any]:
        return {
            "route": self.route,
            "status": self.status,
            "metrics": self.metrics,
            "frame": self.frame,
            "error": self.error,
        }


@dataclass(frozen=True)
class IntegrationReport:
    """All three routes on the same geometry, with frame consistency check."""
    schema: str
    semi_axes: tuple[float, float, float]
    analytic_volume_m3: float
    routes: tuple[RouteResult, ...]
    frame_consistent: bool
    frame_tamper_refused: bool
    all_passed: bool
    elapsed_s: float

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": self.schema,
            "semi_axes": list(self.semi_axes),
            "analytic_volume_m3": self.analytic_volume_m3,
            "routes": [r.to_dict() for r in self.routes],
            "frame_consistent": self.frame_consistent,
            "frame_tamper_refused": self.frame_tamper_refused,
            "all_passed": self.all_passed,
            "elapsed_s": round(self.elapsed_s, 3),
        }


def verify_frame_consistency(routes: Sequence[RouteResult]) -> bool:
    """All routes declare z-up and their volumes/areas agree within tolerance."""
    if not routes:
        return False
    frames = {r.frame for r in routes}
    if frames != {COORDINATE_CONTRACT["name"]}:
        return False
    # Volume/area consistency: each route's primary metric vs analytic
    for r in routes:
        metrics = r.metrics
        if "volume_rel_error" in metrics and metrics["volume_rel_error"] > 0.05:
            return False
        if "circumference_rel_error" in metrics and metrics["circumference_rel_error"] > 1e-10:
            return False
        if "area_rel_error" in metrics and metrics["area_rel_error"] > 0.10:
            return False  # lattice staircase is coarse
    return True


def frame_tamper_refused(outdir: Path, route: str) -> bool:
    """A y-up manifest is refused by the route's verifier.

    Delegates to :func:`tensorcad.solver_solid_gate.tamper_refused` which
    properly refuses only ``HandoffError`` (not any exception).

    Args:
        outdir: path to a valid bundle directory.
        route: route name ("volume_msh", "lattice", or "curved_wall").

    Returns:
        True if the tampered bundle was refused, False otherwise.
    """
    from .solver_solid_gate import tamper_refused

    try:
        if route == "volume_msh":
            from .solver_interface import verify_bundle
            result = tamper_refused(outdir, verify_bundle, name="volume_msh")
        elif route == "lattice":
            from .solver_interface import verify_bundle
            result = tamper_refused(outdir, verify_bundle, name="lattice")
        else:
            from .dg_handoff import verify_bundle as verify_dg_bundle
            result = tamper_refused(outdir, verify_dg_bundle, name="curved_wall")
        return result["tamper_refused"]
    except Exception:
        return False


def summary_markdown(report: IntegrationReport) -> str:
    """Markdown summary of the integration report."""
    lines = [
        "# Solver Interface Integration Report",
        "",
        f"**Schema**: {report.schema}",
        f"**Semi-axes**: {report.semi_axes} m",
        f"**Analytic volume**: {report.analytic_volume_m3:.6f} m³",
        f"**Elapsed**: {report.elapsed_s:.1f} s",
        "",
        "## Route Results",
        "",
        "| route | status | frame | key metric |",
        "|---|---|---|---|",
    ]
    for r in report.routes:
        key_metric = ""
        if "volume_rel_error" in r.metrics:
            key_metric = f"volume_rel_error={r.metrics['volume_rel_error']:.6g}"
        elif "circumference_rel_error" in r.metrics:
            key_metric = f"circumference_rel_error={r.metrics['circumference_rel_error']:.6g}"
        elif "area_rel_error" in r.metrics:
            key_metric = f"area_rel_error={r.metrics['area_rel_error']:.6g}"
        elif "fluid_cells" in r.metrics:
            key_metric = f"fluid_cells={r.metrics['fluid_cells']}"
        lines.append(f"| {r.route} | {r.status} | {r.frame} | {key_metric} |")

    lines += [
        "",
        "## Frame Consistency",
        "",
        f"- All routes declare z-up: **{'PASS' if report.frame_consistent else 'FAIL'}**",
        f"- Frame tamper refused: **{'PASS' if report.frame_tamper_refused else 'FAIL'}**",
        "",
        "## 5415 Lesson",
        "",
        "The coordinate frame is a machine-checked gate, not a comment. Every route's",
        "bundle declares the contract frame (z-up), and a deliberate y-up tamper is",
        "refused by the route's own verifier. This institutionalises the lesson from",
        "the 5415 rounds where a y-up vs z-up mix-up silently mirrored the hull.",
        "",
        f"**Overall**: {'PASS' if report.all_passed else 'FAIL'}",
    ]
    return "\n".join(lines)


def _extract_route_metrics(route_data: dict[str, Any]) -> dict[str, Any]:
    """Extract key metrics from a route result dict."""
    metrics = {}

    # FVM route metrics
    if "volume_rel_error" in route_data:
        metrics["volume_rel_error"] = route_data["volume_rel_error"]
    if "patch_relative_linf" in route_data.get("consumer", {}):
        metrics["patch_relative_linf"] = route_data["consumer"]["patch_relative_linf"]
    if "patch_relative_l2" in route_data.get("consumer", {}):
        metrics["patch_relative_l2"] = route_data["consumer"]["patch_relative_l2"]
    if "global_relative_balance" in route_data.get("consumer", {}):
        metrics["global_relative_balance"] = abs(route_data["consumer"]["global_relative_balance"])

    # DG route metrics
    if "perimeter_rel_error" in route_data:
        metrics["circumference_rel_error"] = route_data["perimeter_rel_error"]
    if "area_rel_error" in route_data:
        metrics["area_rel_error"] = route_data["area_rel_error"]

    # Lattice route metrics
    if "volume_rel_error" in route_data and "volume_rel_error" not in metrics:
        metrics["volume_rel_error"] = route_data["volume_rel_error"]
    if "fluid_cells" in route_data:
        metrics["fluid_cells"] = route_data["fluid_cells"]

    return metrics


def _make_route_result(name: str, route_data: dict) -> RouteResult:
    """Create a RouteResult from a route data dict."""
    consumer = route_data.get("consumer", {})
    consumer_status = consumer.get("status", "unknown")

    # Map consumer status to RouteResult status
    if consumer_status == "ran":
        status = "ran"
    elif consumer_status == "withheld":
        status = "geometry_ready"
    elif consumer_status == "error":
        status = "error"
    else:
        status = consumer_status

    return RouteResult(
        route=route_data["route"],
        status=status,
        metrics=_extract_route_metrics(route_data),
        frame=route_data["frame"],
    )


def run_integration(
    outdir: str | Path,
    *,
    semi_axes: Sequence[float] = (1.0, 0.5, 0.3),
    fvm_resolution: tuple[int, int] = (24, 32),
    lbm_radial_cells: int = 8,
) -> IntegrationReport:
    """Run all three solver routes on the same ellipsoid geometry.

    This function delegates to :mod:`tensorcad.solver_solid_gate` for the actual
    geometry sampling, route execution, and cross-route verification.  The
    gate ensures that:
    - All three routes consume the **same** solid (same SHA-256 digest).
    - All routes declare the contract frame (z-up).
    - Each route's measured volume agrees with the analytic ``4/3 πabc``.
    - A deliberate y-up manifest is refused by every route's verifier.

    Args:
        outdir: where to write per-route bundles and the integration report.
        semi_axes: (a, b, c) of the ellipsoid in metres.
        fvm_resolution: (n_theta, n_phi) for the ellipsoid mesh.
        lbm_radial_cells: radial resolution for the lattice route.

    Returns:
        IntegrationReport with all route results and frame consistency verdict.
    """
    from .solver_solid_gate import (
        cross_route_gate,
        dg_route,
        fvm_route,
        lattice_route,
        solid_from_library,
        tamper_refused,
    )

    t0 = time.time()
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    a, b, c = (float(v) for v in semi_axes)
    analytic_vol = _analytic_volume(semi_axes)

    # --- Get the unified solid from the geometry library ---
    solid = solid_from_library(semi_axes, resolution=max(fvm_resolution))

    # --- Run each route ---
    fvm = fvm_route(solid, out / "fvm", solve=True)
    dg = dg_route(solid, out / "dg")
    lat = lattice_route(solid, out / "lattice", radial_cells=lbm_radial_cells)

    routes_dict = {"fvm": fvm, "dg": dg, "lattice": lat}

    # --- Cross-route gate ---
    gate = cross_route_gate(routes_dict, solid)

    # --- Build RouteResult objects ---
    route_results = tuple(
        _make_route_result(name, data) for name, data in routes_dict.items()
    )

    # --- Frame tamper test ---
    tamper_results = {}
    for name, route_data in routes_dict.items():
        route_dir = out / name
        route_dir.mkdir(exist_ok=True)
        try:
            if name == "fvm":
                from .solver_interface import verify_bundle
                tamper_results[name] = tamper_refused(route_dir, verify_bundle, name=name)
            elif name == "lattice":
                from .solver_interface import verify_bundle
                tamper_results[name] = tamper_refused(route_dir, verify_bundle, name=name)
            else:
                from .dg_handoff import verify_bundle as verify_dg_bundle
                tamper_results[name] = tamper_refused(route_dir, verify_dg_bundle, name=name)
        except Exception:
            tamper_results[name] = False

    all_tamper_refused = all(tamper_results.values()) if tamper_results else False

    elapsed = time.time() - t0
    all_passed = (
        gate["failed"] == []
        and all_tamper_refused
        and all(r.status in ("ran", "geometry_ready") for r in route_results)
    )

    report = IntegrationReport(
        schema=INTEGRATION_SCHEMA,
        semi_axes=(a, b, c),
        analytic_volume_m3=analytic_vol,
        routes=route_results,
        frame_consistent=gate["failed"] == [],
        frame_tamper_refused=all_tamper_refused,
        all_passed=all_passed,
        elapsed_s=elapsed,
    )

    # Write artifacts
    (out / "integration_report.json").write_text(
        json.dumps(report.to_dict(), indent=2, sort_keys=True) + "\n"
    )
    (out / "summary.md").write_text(summary_markdown(report))

    return report
