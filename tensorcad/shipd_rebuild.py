"""ShipD rebuild pipeline: parameter vectors → hull mesh → quality gate → report.

This module implements the GAP-CAD-C2-SHIPD-REBUILD capability:
- Takes a database of ship parameter vectors (L, B, T, Cb, Cm, Cwp, …)
- Reconstructs hull geometry using the best available parametric method
- Runs quality gates (watertight, manifold, volume consistency, Cb/L/B sanity)
- Produces a pass-rate report, correlation analysis, and failure-mode ledger

Reconstruction methods (selected automatically per ship):
1. **series60_form** — Series 60 section law (Cb ∈ [0.55, 0.85]); validated parent forms.
2. **wigley_form** — Wigley mathematical hull (any Cb, but simplified geometry).
3. **parametric_sections** — Generalized section loft with independent Cwp/Cb (Cb ∈ [0.45, 0.90]).
4. **refused** — Parameters outside all method windows; no geometry produced.

Quality gates (all must pass for VERIFIED):
- G1: Watertight (every edge used by exactly 2 faces)
- G2: Edge-manifold (no edge used by >2 faces)
- G3: Euler characteristic V−E+F = 2 (single closed surface)
- G4: Outward-oriented (signed volume > 0)
- G5: No zero-area faces
- G6: Volume consistency |V_mesh − Cb·L·B·T| / (Cb·L·B·T) < 3%
- G7: Realised L/B/T agree with requested within 2%
- G8: Triangle count ≥ 1 000 (not a placeholder)

External anchor: Ship parameter distributions from naval architecture literature
(Schneekluth & Bertram 1998; Watson 1998; IMO ship type statistics).
The 500-ship synthetic database uses these distributions as the population model.
"""
from __future__ import annotations

import csv
import json
import math
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Sequence

from tensorcad import mesh as meshlib
from tensorcad.cb_sections import (
    LcbRequestUnreachable,
    CalibratedHullParameters,
    build_calibrated_hull,
    realised_coefficients,
    relative_errors,
)
from tensorcad.geometry import Mesh
from tensorcad.hull import HullParameters, generate_hull
from tensorcad.series60 import Series60Parameters, build_series60_mesh
from tensorcad.wigley import WigleyParameters, generate_wigley
from tensorcad.hydrostatics import compute_hydrostatics
from tensorcad.shipd_population import correlation_sanity, nominal_performance_fields

SCHEMA = "tensorcad.shipd_rebuild/1.0"

# --------------------------------------------------------------------------- #
# Reconstruction method windows
# --------------------------------------------------------------------------- #
SERIES60_CB_WINDOW = (0.55, 0.85)
PARAMETRIC_CB_WINDOW = (0.45, 0.90)
WIGLEY_CB_WINDOW = (0.30, 0.95)

# --------------------------------------------------------------------------- #
# Quality gate thresholds
# --------------------------------------------------------------------------- #
VOLUME_TOLERANCE = 0.03       # 3% volume consistency
DIMENSION_TOLERANCE = 0.02    # 2% L/B/T agreement
MIN_TRIANGLES = 1_000         # Not a placeholder
DEFAULT_NX = 81               # Longitudinal stations
DEFAULT_NZ = 21               # Vertical levels

#: Block-coefficient window the coefficient-calibrated family claims.
CALIBRATED_CB_WINDOW = (0.30, 0.95)
#: Tolerance on the realised midship / waterplane coefficients (percent of the
#: requested value). Looser than the 3 % volume gate because these two are read
#: off the *sampled* station grid rather than integrated over the closed surface.
COEFFICIENT_TOLERANCE = 0.05
#: Longitudinal centre of buoyancy, as a fraction of L: a requested LCB is honoured
#: here or the hull is not the ship that was asked for. The band is a policy number
#: (~500x the worst realised error measured in round 7, 2.7e-05 L); the report
#: quotes the actual errors, the gate only decides pass/fail.
LCB_TOLERANCE = 0.005
#: Seawater density [kg/m^3] and kinematic viscosity [m^2/s] used by the nominal
#: performance block (ITTC fresh/seawater tables at 15 degC).
SEAWATER_DENSITY = 1025.0
SEAWATER_VISCOSITY = 1.1883e-6
GRAVITY = 9.80665


@dataclass(frozen=True)
class GateResult:
    name: str
    passed: bool
    detail: str
    value: float | None = None
    threshold: float | None = None

    def as_dict(self) -> dict[str, Any]:
        d = {"name": self.name, "passed": self.passed, "detail": self.detail}
        if self.value is not None:
            d["value"] = self.value
        if self.threshold is not None:
            d["threshold"] = self.threshold
        return d


@dataclass(frozen=True)
class RebuildResult:
    """Result of rebuilding one ship from its parameter vector."""
    ship_id: int
    ship_type: str
    method: str
    method_reason: str
    mesh: Mesh | None = None
    gates: list[GateResult] = field(default_factory=list)
    verdict: str = "PENDING"
    volume_mesh: float | None = None
    volume_parametric: float | None = None
    volume_error: float | None = None
    realised_L: float | None = None
    realised_B: float | None = None
    realised_T: float | None = None
    triangle_count: int | None = None
    error_message: str | None = None
    rebuild_time_s: float | None = None
    #: Coefficient agreement measured on the mesh (not assumed from the inputs).
    realised_Cb: float | None = None
    realised_Cm: float | None = None
    realised_Cwp: float | None = None
    coefficient_gates: list[GateResult] = field(default_factory=list)
    wetted_surface_area: float | None = None
    #: LCB contract: requested (fraction of L, signed, + forward), realised on the
    #: mesh by the divergence theorem, and their difference.
    requested_lcb: float | None = None
    realised_lcb: float | None = None
    lcb_error: float | None = None
    lcb_gates: list[GateResult] = field(default_factory=list)
    #: ``VERIFIED`` only when the eight canonical gates *and* the coefficient
    #: agreement gates pass. :attr:`verdict` keeps the canonical eight so round-2
    #: and round-3 pass rates stay directly comparable.
    strict_verdict: str = "PENDING"

    @property
    def all_gates_passed(self) -> bool:
        return all(g.passed for g in self.gates)

    def as_dict(self) -> dict[str, Any]:
        return {
            "ship_id": self.ship_id,
            "ship_type": self.ship_type,
            "method": self.method,
            "method_reason": self.method_reason,
            "verdict": self.verdict,
            "triangle_count": self.triangle_count,
            "volume_mesh": self.volume_mesh,
            "volume_parametric": self.volume_parametric,
            "volume_error": self.volume_error,
            "realised_L": self.realised_L,
            "realised_B": self.realised_B,
            "realised_T": self.realised_T,
            "gates": [g.as_dict() for g in self.gates],
            "error_message": self.error_message,
            "rebuild_time_s": self.rebuild_time_s,
            "realised_Cb": self.realised_Cb,
            "realised_Cm": self.realised_Cm,
            "realised_Cwp": self.realised_Cwp,
            "coefficient_gates": [g.as_dict() for g in self.coefficient_gates],
            "wetted_surface_area": self.wetted_surface_area,
            "requested_lcb": self.requested_lcb,
            "realised_lcb": self.realised_lcb,
            "lcb_error": self.lcb_error,
            "lcb_gates": [g.as_dict() for g in self.lcb_gates],
            "strict_verdict": self.strict_verdict,
        }


def _mesh_bbox(mesh: Mesh) -> tuple[float, float, float, float, float, float]:
    """Compute (xmin, ymin, zmin, xmax, ymax, zmax) for a mesh."""
    xs = [v[0] for v in mesh.vertices]
    ys = [v[1] for v in mesh.vertices]
    zs = [v[2] for v in mesh.vertices]
    return (min(xs), min(ys), min(zs), max(xs), max(ys), max(zs))


def _euler_characteristic(mesh: Mesh) -> int:
    """Compute Euler characteristic V - E + F for a mesh."""
    edge_set = set()
    for f in mesh.faces:
        for i in range(3):
            a, b = f[i], f[(i + 1) % 3]
            edge_set.add((min(a, b), max(a, b)))
    return len(mesh.vertices) - len(edge_set) + len(mesh.faces)


def _calibrated_feasible(cb: float, cm: float, cwp: float) -> tuple[bool, str]:
    """Can the coefficient-calibrated family realise this parameter triple?

    The family needs ``Cp = Cb/Cm <= Cwp`` (a waterplane cannot be finer than the
    prismatic coefficient) and three coefficients strictly inside ``(0, 1)``.
    """
    if not all(math.isfinite(v) and 0.0 < v < 1.0 for v in (cb, cm, cwp)):
        return False, f"coefficients outside (0,1): Cb={cb}, Cm={cm}, Cwp={cwp}"
    low, high = CALIBRATED_CB_WINDOW
    if not (low <= cb <= high):
        return False, f"Cb={cb:.4f} outside the calibrated window [{low}, {high}]"
    cp = cb / cm
    if cp > cwp:
        return False, f"Cp=Cb/Cm={cp:.4f} exceeds Cwp={cwp:.4f}"
    return True, (
        f"Cb={cb:.4f} calibrated family: Cp={cp:.4f} <= Cwp={cwp:.4f} "
        f"(realises Cb, Cm, Cwp exactly)"
    )


def _select_method(
    cb: float, cm: float, cwp: float, *, method_preference: str = "calibrated"
) -> tuple[str, str]:
    """Select the best reconstruction method for a ship's parameters.

    ``method_preference="calibrated"`` (default) puts the coefficient-calibrated
    family first: it is the only method that honours the requested ``Cb``, ``Cm``
    and ``Cwp`` instead of producing whatever the fixed breadth law gives.
    ``method_preference="legacy"`` reproduces the round-2 chain
    (Series 60 -> parametric sections -> Wigley -> refused) unchanged, which is
    what the round-2 evidence was collected with.
    """
    if method_preference == "calibrated":
        feasible, reason = _calibrated_feasible(cb, cm, cwp)
        if feasible:
            return "cb_calibrated_sections", reason

    low_s60, high_s60 = SERIES60_CB_WINDOW
    low_param, high_param = PARAMETRIC_CB_WINDOW
    low_wig, high_wig = WIGLEY_CB_WINDOW

    if low_s60 <= cb <= high_s60:
        if cm < 1.0:
            return "series60_form", f"Cb={cb:.4f} in Series 60 window [{low_s60}, {high_s60}]"
        else:
            return "parametric_sections", f"Cb={cb:.4f} in Series 60 window but Cm={cm:.4f}≈1"

    if low_param <= cb <= high_param:
        return "parametric_sections", f"Cb={cb:.4f} in parametric window [{low_param}, {high_param}]"

    if low_wig <= cb <= high_wig:
        return "wigley_form", f"Cb={cb:.4f} in Wigley window [{low_wig}, {high_wig}]"

    return "refused", f"Cb={cb:.4f} outside all method windows [{low_wig}, {high_wig}]"


def _rebuild_series60(ship: dict, nx: int, nz: int) -> RebuildResult:
    """Rebuild using Series 60 section law."""
    t0 = time.time()
    try:
        params = Series60Parameters(
            length=ship["L"], beam=ship["B"], draft=ship["T"],
            cb=ship["Cb"], nx=nx, nz=nz,
        )
        errors, warnings = params.problems()
        if errors:
            return RebuildResult(
                ship_id=ship["ship_id"], ship_type=ship["ship_type"],
                method="series60_form", method_reason="Series 60 selected but parameters rejected",
                verdict="FAILED", error_message="; ".join(errors),
                rebuild_time_s=time.time() - t0,
            )
        mesh, offsets = build_series60_mesh(params)
        elapsed = time.time() - t0
        return _evaluate_gates(ship, mesh, "series60_form",
                               f"Cb={ship['Cb']:.4f} in Series 60 window", elapsed)
    except Exception as e:
        return RebuildResult(
            ship_id=ship["ship_id"], ship_type=ship["ship_type"],
            method="series60_form", method_reason="Series 60 selected",
            verdict="FAILED", error_message=str(e),
            rebuild_time_s=time.time() - t0,
        )


def _rebuild_wigley(ship: dict, nx: int, nz: int) -> RebuildResult:
    """Rebuild using Wigley mathematical hull."""
    t0 = time.time()
    try:
        params = WigleyParameters(
            length=ship["L"], beam=ship["B"], draft=ship["T"],
            nx=nx, nz=nz,
        )
        model = generate_wigley(params)
        elapsed = time.time() - t0
        return _evaluate_gates(ship, model.mesh, "wigley_form",
                               f"Cb={ship['Cb']:.4f} in Wigley window", elapsed)
    except Exception as e:
        return RebuildResult(
            ship_id=ship["ship_id"], ship_type=ship["ship_type"],
            method="wigley_form", method_reason="Wigley selected as fallback",
            verdict="FAILED", error_message=str(e),
            rebuild_time_s=time.time() - t0,
        )


def _rebuild_parametric(ship: dict, nx: int, nz: int) -> RebuildResult:
    """Rebuild using generalized parametric sections."""
    t0 = time.time()
    try:
        params = HullParameters(
            length=ship["L"], beam=ship["B"], draft=ship["T"],
            nx=nx, nz=nz,
        )
        model = generate_hull(params)
        elapsed = time.time() - t0
        return _evaluate_gates(ship, model.mesh, "parametric_sections",
                               f"Cb={ship['Cb']:.4f} in parametric window", elapsed)
    except Exception as e:
        return RebuildResult(
            ship_id=ship["ship_id"], ship_type=ship["ship_type"],
            method="parametric_sections", method_reason="Parametric sections selected",
            verdict="FAILED", error_message=str(e),
            rebuild_time_s=time.time() - t0,
        )


def _rebuild_calibrated(ship: dict, nx: int, nz: int) -> RebuildResult:
    """Rebuild with the coefficient-calibrated family (realises Cb, Cm, Cwp)."""
    t0 = time.time()
    try:
        params = CalibratedHullParameters(
            length=ship["L"], beam=ship["B"], draft=ship["T"],
            block_coefficient=ship["Cb"],
            midship_coefficient=ship.get("Cm", 0.98),
            waterplane_coefficient=ship.get("Cwp", 0.80),
            nx=nx, nz=nz,
            # Optional hull-form features; absent -> the round-3 flat-keel family.
            parallel_body_fraction=float(ship.get("parallel_body_fraction", 0.0)),
            end_draft_ratio=float(ship.get("end_draft_ratio", 1.0)),
            # ``lcb`` absent -> 0.0 -> the round-3 symmetric family, bit for bit.
            lcb_offset=float(ship.get("lcb", 0.0) or 0.0),
        )
        problems = params.problems()
        if problems:
            return RebuildResult(
                ship_id=ship["ship_id"], ship_type=ship["ship_type"],
                method="cb_calibrated_sections",
                method_reason="calibrated family selected but parameters rejected",
                verdict="REFUSED", error_message="; ".join(problems),
                rebuild_time_s=time.time() - t0,
            )
        hull = build_calibrated_hull(params)
        elapsed = time.time() - t0
        return _evaluate_gates(
            ship, hull.mesh, "cb_calibrated_sections",
            f"Cb={ship['Cb']:.4f}, Cp<=Cwp -> calibrated family", elapsed,
        )
    except LcbRequestUnreachable as e:
        # A refused LCB is a parameter-vector finding, not a mesh defect: nothing was
        # built, so there is no surface to measure. Classified separately so the
        # failure ledger distinguishes "unreachable target" from "mesh broke a gate".
        return RebuildResult(
            ship_id=ship["ship_id"], ship_type=ship["ship_type"],
            method="cb_calibrated_sections",
            method_reason="calibrated family: requested LCB outside the reachable band",
            verdict="REFUSED", error_message=f"lcb_unreachable: {e}",
            requested_lcb=float(ship.get("lcb", 0.0) or 0.0),
            rebuild_time_s=time.time() - t0,
        )
    except Exception as e:  # pragma: no cover - defensive
        return RebuildResult(
            ship_id=ship["ship_id"], ship_type=ship["ship_type"],
            method="cb_calibrated_sections", method_reason="calibrated family",
            verdict="FAILED", error_message=str(e), rebuild_time_s=time.time() - t0,
        )


def _evaluate_gates(ship: dict, mesh: Mesh, method: str, reason: str, elapsed: float) -> RebuildResult:
    """Run all quality gates on a rebuilt mesh."""
    gates = []

    # G1: Watertight
    wt = meshlib.is_watertight(mesh)
    gates.append(GateResult("G1_watertight", wt,
                            "watertight" if wt else "not watertight"))

    # G2: Edge-manifold
    em = meshlib.is_edge_manifold(mesh)
    gates.append(GateResult("G2_edge_manifold", em,
                            "edge-manifold" if em else "non-manifold edges"))

    # G3: Euler characteristic
    euler = _euler_characteristic(mesh)
    euler_ok = euler == 2
    gates.append(GateResult("G3_euler_char", euler_ok,
                            f"V-E+F={euler} (expected 2)", value=float(euler), threshold=2.0))

    # G4: Outward-oriented -- and consistently wound, which the sign alone cannot
    # establish. Flipping one triangle moves the divergence-theorem sum by that one
    # tetrahedron: the total stays positive, the topology (G1-G3), the face areas
    # (G5) and the volume/moment tolerances (G6-G8) are all unchanged, so a mesh
    # with an inverted patch used to pass every gate. The directed-edge test in
    # ``mesh.edge_winding_conflicts`` sees it, so both are required here.
    vol = meshlib.total_signed_volume(mesh)
    conflicts = meshlib.edge_winding_conflicts(mesh)
    oriented = vol > 0 and not conflicts
    gates.append(GateResult("G4_outward_oriented", oriented,
                            f"signed volume={vol:.4f} m³, "
                            f"{len(conflicts)} winding-inconsistent edges",
                            value=vol, threshold=0.0))

    # G5: No zero-area faces
    zero_area = len(meshlib.zero_area_faces(mesh))
    no_zero = zero_area == 0
    gates.append(GateResult("G5_no_zero_area", no_zero,
                            f"{zero_area} zero-area faces", value=float(zero_area), threshold=0.0))

    # G6: Volume consistency
    vol_parametric = ship["Cb"] * ship["L"] * ship["B"] * ship["T"]
    vol_error = abs(vol - vol_parametric) / vol_parametric if vol_parametric > 0 else float("inf")
    vol_ok = vol_error < VOLUME_TOLERANCE
    gates.append(GateResult("G6_volume_consistency", vol_ok,
                            f"|V_mesh - V_param|/V_param = {vol_error:.4f}",
                            value=vol_error, threshold=VOLUME_TOLERANCE))

    # G7: Realised L/B/T agreement
    bbox = _mesh_bbox(mesh)
    rL = bbox[3] - bbox[0]
    rB = 2 * max(abs(bbox[1]), abs(bbox[4]))
    rT = bbox[5] - bbox[2]
    L_err = abs(rL - ship["L"]) / ship["L"] if ship["L"] > 0 else float("inf")
    B_err = abs(rB - ship["B"]) / ship["B"] if ship["B"] > 0 else float("inf")
    T_err = abs(rT - ship["T"]) / ship["T"] if ship["T"] > 0 else float("inf")
    dim_ok = L_err < DIMENSION_TOLERANCE and B_err < DIMENSION_TOLERANCE and T_err < DIMENSION_TOLERANCE
    gates.append(GateResult("G7_dimension_agreement", dim_ok,
                            f"L_err={L_err:.4f}, B_err={B_err:.4f}, T_err={T_err:.4f}",
                            value=max(L_err, B_err, T_err), threshold=DIMENSION_TOLERANCE))

    # G8: Triangle count
    tri_count = len(mesh.faces)
    tri_ok = tri_count >= MIN_TRIANGLES
    gates.append(GateResult("G8_triangle_count", tri_ok,
                            f"{tri_count} triangles (min={MIN_TRIANGLES})",
                            value=float(tri_count), threshold=float(MIN_TRIANGLES)))

    all_passed = all(g.passed for g in gates)
    verdict = "VERIFIED" if all_passed else "FAILED"

    # Coefficient agreement on the mesh itself. ``Cm`` is the largest sampled
    # station area and ``Cwp`` the waterplane flux; both are measured, never
    # copied from the request, so a surface that drifts off the parameter vector
    # shows up here even when the volume happens to be right.
    coefficient_gates: list[GateResult] = []
    realised_cb = realised_cm = realised_cwp = None
    realised_lcb = None
    wetted = None
    try:
        measured = realised_coefficients(
            mesh, length=ship["L"], beam=ship["B"], draft=ship["T"]
        )
        realised_cb = measured["coefficients"]["Cb"]
        realised_cm = measured["coefficients"]["Cm"]
        realised_cwp = measured["coefficients"]["Cwp"]
        wetted = measured["wetted_surface_m2"]
        # Measured by the independent divergence-theorem centroid in
        # ``cb_sections.realised_lcb_fraction``, never copied from the request.
        realised_lcb = measured.get("lcb_fraction")
        for name, got, want in (
            ("G9_midship_coefficient", realised_cm, ship.get("Cm", 0.98)),
            ("G10_waterplane_coefficient", realised_cwp, ship.get("Cwp", 0.80)),
        ):
            err = abs(got - want) / want if want else float("inf")
            coefficient_gates.append(GateResult(
                name, err < COEFFICIENT_TOLERANCE,
                f"realised={got:.5f} requested={want:.5f} err={err:.4%}",
                value=err, threshold=COEFFICIENT_TOLERANCE,
            ))
    except Exception as exc:  # pragma: no cover - defensive
        coefficient_gates.append(GateResult(
            "G9_midship_coefficient", False, f"measurement failed: {exc}",
        ))

    # LCB contract. Only a vector that *carries* an ``lcb`` target can be held to one:
    # the legacy Series 60 / Wigley / parametric chains have no LCB degree of freedom,
    # so asking them for one is reported as unsupported rather than silently passed.
    lcb_gates: list[GateResult] = []
    requested_lcb = ship.get("lcb")
    lcb_error = None
    if requested_lcb is not None:
        requested_lcb = float(requested_lcb)
        if method != "cb_calibrated_sections":
            lcb_gates.append(GateResult(
                "G11_longitudinal_centre", False,
                f"method {method} has no LCB control; requested "
                f"{requested_lcb:+.5f} L cannot be honoured",
            ))
        elif realised_lcb is None:
            lcb_gates.append(GateResult(
                "G11_longitudinal_centre", False, "LCB measurement failed",
            ))
        else:
            lcb_error = abs(realised_lcb - requested_lcb)
            lcb_gates.append(GateResult(
                "G11_longitudinal_centre", lcb_error < LCB_TOLERANCE,
                f"realised={realised_lcb:+.6f} L requested={requested_lcb:+.6f} L "
                f"err={lcb_error:.3e} L",
                value=lcb_error, threshold=LCB_TOLERANCE,
            ))

    strict_ok = (
        all_passed
        and all(g.passed for g in coefficient_gates)
        and all(g.passed for g in lcb_gates)
    )
    return RebuildResult(
        ship_id=ship["ship_id"], ship_type=ship["ship_type"],
        method=method, method_reason=reason,
        mesh=mesh, gates=gates, verdict=verdict,
        volume_mesh=vol, volume_parametric=vol_parametric,
        volume_error=vol_error,
        realised_L=rL, realised_B=rB, realised_T=rT,
        triangle_count=tri_count, rebuild_time_s=elapsed,
        realised_Cb=realised_cb, realised_Cm=realised_cm, realised_Cwp=realised_cwp,
        coefficient_gates=coefficient_gates, wetted_surface_area=wetted,
        requested_lcb=None if requested_lcb is None else float(requested_lcb),
        realised_lcb=realised_lcb, lcb_error=lcb_error, lcb_gates=lcb_gates,
        strict_verdict="VERIFIED" if strict_ok else "FAILED",
    )


def rebuild_ship(
    ship: dict,
    *,
    nx: int = DEFAULT_NX,
    nz: int = DEFAULT_NZ,
    method_preference: str = "calibrated",
) -> RebuildResult:
    """Rebuild one ship from its parameter vector."""
    method, reason = _select_method(
        ship["Cb"], ship.get("Cm", 0.98), ship.get("Cwp", 0.80),
        method_preference=method_preference,
    )

    if method == "refused":
        return RebuildResult(
            ship_id=ship["ship_id"], ship_type=ship["ship_type"],
            method="refused", method_reason=reason,
            verdict="REFUSED", error_message=reason,
        )

    if method == "cb_calibrated_sections":
        return _rebuild_calibrated(ship, nx, nz)
    if method == "series60_form":
        return _rebuild_series60(ship, nx, nz)
    elif method == "wigley_form":
        return _rebuild_wigley(ship, nx, nz)
    elif method == "parametric_sections":
        return _rebuild_parametric(ship, nx, nz)
    else:
        return RebuildResult(
            ship_id=ship["ship_id"], ship_type=ship["ship_type"],
            method="unknown", method_reason=f"Unknown method: {method}",
            verdict="REFUSED", error_message=f"Unknown method: {method}",
        )


def _rebuild_one(args: tuple[dict, int, int, str]) -> "RebuildResult":
    """Module-level trampoline so the rebuild loop can run in a process pool."""
    ship, nx, nz, method_preference = args
    return rebuild_ship(ship, nx=nx, nz=nz, method_preference=method_preference)


def run_shipd_rebuild(
    ships: Sequence[dict],
    *,
    nx: int = DEFAULT_NX,
    nz: int = DEFAULT_NZ,
    out_dir: str | Path = "runs/shipd_rebuild",
    method_preference: str = "calibrated",
    workers: int = 1,
) -> dict[str, Any]:
    """Run the full ShipD rebuild pipeline.

    ``method_preference`` selects the reconstruction chain: ``"calibrated"``
    (default) tries the coefficient-calibrated family first, ``"legacy"``
    reproduces the round-2 Series 60 / parametric / Wigley chain.

    ``workers > 1`` runs the per-ship rebuild in a process pool; the assembly and
    the artefacts are identical either way (each ship is independent).
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    results: list[RebuildResult] = []
    t_start = time.time()

    if workers and workers > 1 and len(ships) > 1:
        import multiprocessing

        with multiprocessing.Pool(workers) as pool:
            results = pool.map(
                _rebuild_one,
                [(ship, nx, nz, method_preference) for ship in ships],
                chunksize=max(1, len(ships) // (workers * 8)),
            )
    else:
        for ship in ships:
            result = rebuild_ship(ship, nx=nx, nz=nz, method_preference=method_preference)
            results.append(result)

    total_time = time.time() - t_start

    n_total = len(results)
    n_verified = sum(1 for r in results if r.verdict == "VERIFIED")
    n_failed = sum(1 for r in results if r.verdict == "FAILED")
    n_refused = sum(1 for r in results if r.verdict == "REFUSED")
    pass_rate = n_verified / n_total if n_total > 0 else 0.0

    method_counts: dict[str, int] = {}
    for r in results:
        method_counts[r.method] = method_counts.get(r.method, 0) + 1

    type_counts: dict[str, int] = {}
    for r in results:
        type_counts[r.ship_type] = type_counts.get(r.ship_type, 0) + 1

    gate_names = ["G1_watertight", "G2_edge_manifold", "G3_euler_char",
                  "G4_outward_oriented", "G5_no_zero_area", "G6_volume_consistency",
                  "G7_dimension_agreement", "G8_triangle_count"]
    gate_pass_counts = {g: 0 for g in gate_names}
    for r in results:
        for g in r.gates:
            if g.name in gate_pass_counts:
                gate_pass_counts[g.name] += 1 if g.passed else 0

    gate_pass_rates = {g: gate_pass_counts[g] / n_total if n_total > 0 else 0.0
                       for g in gate_names}

    failure_modes: dict[str, list[dict]] = {}
    for r in results:
        if r.verdict in ("FAILED", "REFUSED"):
            key = r.method + ":" + (r.error_message or "unknown")[:80]
            if key not in failure_modes:
                failure_modes[key] = []
            failure_modes[key].append({"ship_id": r.ship_id, "ship_type": r.ship_type, "error": r.error_message})

    n_strict = sum(1 for r in results if r.strict_verdict == "VERIFIED")
    strict_pass_rate = n_strict / n_total if n_total > 0 else 0.0

    # Realised-vs-requested coefficient errors, ship by ship.
    coefficient_errors: dict[str, list[float]] = {"Cb": [], "Cm": [], "Cwp": []}
    for ship, r in zip(ships, results):
        if r.realised_Cb is None:
            continue
        for key, got in (("Cb", r.realised_Cb), ("Cm", r.realised_Cm), ("Cwp", r.realised_Cwp)):
            want = ship[key]
            if want:
                coefficient_errors[key].append(abs(got - want) / want)

    def _coefficient_summary(values: list[float]) -> dict[str, float]:
        if not values:
            return {"mean": float("nan"), "max": float("nan"), "n": 0}
        ordered = sorted(values)
        return {
            "n": len(values),
            "mean": sum(values) / len(values),
            "median": ordered[len(ordered) // 2],
            "max": ordered[-1],
            "p95": ordered[int(0.95 * (len(ordered) - 1))],
            "within_tolerance_fraction": sum(
                1 for v in values if v < COEFFICIENT_TOLERANCE
            ) / len(values),
        }

    # LCB contract accounting: requested / honoured / refused / unsupported.
    lcb_rows: list[dict[str, Any]] = []
    for ship, r in zip(ships, results):
        if r.requested_lcb is None:
            continue
        lcb_rows.append({
            "ship_id": r.ship_id,
            "ship_type": r.ship_type,
            "Cb": ship["Cb"],
            "requested_lcb": r.requested_lcb,
            "realised_lcb": r.realised_lcb,
            "lcb_error": r.lcb_error,
            "verdict": r.verdict,
            "method": r.method,
            "lcb_gate_passed": (r.lcb_gates[0].passed if r.lcb_gates else None),
            "error_message": r.error_message,
        })
    lcb_errors = [row["lcb_error"] for row in lcb_rows if row["lcb_error"] is not None]
    lcb_summary = {
        "tolerance_L": LCB_TOLERANCE,
        "n_requested": len(lcb_rows),
        "n_honoured": sum(1 for row in lcb_rows if row["lcb_gate_passed"]),
        "n_refused_unreachable": sum(
            1 for row in lcb_rows
            if (row["error_message"] or "").startswith("lcb_unreachable")
        ),
        "n_unsupported_method": sum(
            1 for row in lcb_rows
            if row["method"] != "cb_calibrated_sections" and row["realised_lcb"] is not None
        ),
        "n_gate_failed": sum(1 for row in lcb_rows if row["lcb_gate_passed"] is False),
        "max_abs_error_L": max(lcb_errors) if lcb_errors else float("nan"),
        "median_abs_error_L": (
            sorted(lcb_errors)[len(lcb_errors) // 2] if lcb_errors else float("nan")
        ),
        "mean_signed_error_L": (
            sum(row["realised_lcb"] - row["requested_lcb"]
                for row in lcb_rows if row["realised_lcb"] is not None)
            / len([row for row in lcb_rows if row["realised_lcb"] is not None])
            if any(row["realised_lcb"] is not None for row in lcb_rows)
            else float("nan")
        ),
    }
    (out_dir / "lcb_contract.json").write_text(json.dumps(
        {"summary": lcb_summary, "rows": lcb_rows}, indent=2))

    nominal_fields = {
        ship["ship_id"]: nominal_performance_fields(ship) for ship in ships
    }
    (out_dir / "nominal_fields.json").write_text(json.dumps(nominal_fields, indent=2))
    (out_dir / "ship_database.json").write_text(json.dumps(list(ships), indent=2))

    correlation = correlation_sanity(ships, results)
    (out_dir / "correlation.json").write_text(json.dumps(correlation, indent=2))

    vol_errors = [r.volume_error for r in results if r.volume_error is not None and r.verdict == "VERIFIED"]
    vol_error_mean = sum(vol_errors) / len(vol_errors) if vol_errors else 0.0
    vol_error_max = max(vol_errors) if vol_errors else 0.0

    cb_errors = [(s["Cb"], r.volume_error) for s, r in zip(ships, results) if r.volume_error is not None]

    evidence = {
        "schema": SCHEMA,
        "n_ships": n_total,
        "n_verified": n_verified,
        "n_failed": n_failed,
        "n_refused": n_refused,
        "pass_rate": pass_rate,
        "strict_pass_rate": strict_pass_rate,
        "n_strict_verified": n_strict,
        "method_preference": method_preference,
        "coefficient_agreement": {
            key: _coefficient_summary(values) for key, values in coefficient_errors.items()
        },
        "coefficient_tolerance": COEFFICIENT_TOLERANCE,
        "lcb_contract": lcb_summary,
        "correlation_sanity": correlation,
        "method_distribution": method_counts,
        "ship_type_distribution": type_counts,
        "gate_pass_rates": gate_pass_rates,
        "volume_error_mean": vol_error_mean,
        "volume_error_max": vol_error_max,
        "total_time_s": total_time,
        "settings": {"nx": nx, "nz": nz, "volume_tolerance": VOLUME_TOLERANCE,
                     "dimension_tolerance": DIMENSION_TOLERANCE, "min_triangles": MIN_TRIANGLES},
        "failure_mode_ledger": {k: v[:10] for k, v in failure_modes.items()},
        "cb_volume_correlation": cb_errors[:100],
    }

    # Write artefacts
    ship_results = [r.as_dict() for r in results]
    (out_dir / "ship_results.json").write_text(json.dumps(ship_results, indent=2))
    (out_dir / "evidence.json").write_text(json.dumps(evidence, indent=2))

    csv_path = out_dir / "summary.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["ship_id", "ship_type", "L", "B", "T", "Cb", "method", "verdict",
                         "triangle_count", "volume_mesh", "volume_parametric", "volume_error",
                         "realised_L", "realised_B", "realised_T", "rebuild_time_s",
                         "requested_lcb", "realised_lcb", "lcb_error"])
        for r, s in zip(results, ships):
            writer.writerow([
                r.ship_id, r.ship_type, s["L"], s["B"], s["T"], s["Cb"],
                r.method, r.verdict, r.triangle_count, r.volume_mesh,
                r.volume_parametric, r.volume_error,
                r.realised_L, r.realised_B, r.realised_T, r.rebuild_time_s,
                r.requested_lcb, r.realised_lcb, r.lcb_error,
            ])

    failure_ledger = []
    for r in results:
        if r.verdict in ("FAILED", "REFUSED"):
            orig = next((s for s in ships if s["ship_id"] == r.ship_id), None)
            failure_ledger.append({
                "ship_id": r.ship_id, "ship_type": r.ship_type,
                "Cb": orig["Cb"] if orig else None,
                "method": r.method, "verdict": r.verdict,
                "error_message": r.error_message,
                "requested_lcb": r.requested_lcb,
                "realised_lcb": r.realised_lcb,
                "failed_gates": [g.name for g in r.gates if not g.passed]
                                + [g.name for g in r.lcb_gates if not g.passed],
                "failed_lcb_gate": any(not g.passed for g in r.lcb_gates),
            })
    (out_dir / "failure_ledger.json").write_text(json.dumps(failure_ledger, indent=2))

    return evidence
