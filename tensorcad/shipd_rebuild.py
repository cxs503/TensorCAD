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
from tensorcad.geometry import Mesh
from tensorcad.hull import HullParameters, generate_hull
from tensorcad.series60 import Series60Parameters, build_series60_mesh
from tensorcad.wigley import WigleyParameters, generate_wigley
from tensorcad.hydrostatics import compute_hydrostatics

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


def _select_method(cb: float, cm: float, cwp: float) -> tuple[str, str]:
    """Select the best reconstruction method for a ship's parameters."""
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

    # G4: Outward-oriented
    vol = meshlib.total_signed_volume(mesh)
    oriented = vol > 0
    gates.append(GateResult("G4_outward_oriented", oriented,
                            f"signed volume={vol:.4f} m³",
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

    return RebuildResult(
        ship_id=ship["ship_id"], ship_type=ship["ship_type"],
        method=method, method_reason=reason,
        mesh=mesh, gates=gates, verdict=verdict,
        volume_mesh=vol, volume_parametric=vol_parametric,
        volume_error=vol_error,
        realised_L=rL, realised_B=rB, realised_T=rT,
        triangle_count=tri_count, rebuild_time_s=elapsed,
    )


def rebuild_ship(ship: dict, *, nx: int = DEFAULT_NX, nz: int = DEFAULT_NZ) -> RebuildResult:
    """Rebuild one ship from its parameter vector."""
    method, reason = _select_method(ship["Cb"], ship.get("Cm", 0.98), ship.get("Cwp", 0.80))

    if method == "refused":
        return RebuildResult(
            ship_id=ship["ship_id"], ship_type=ship["ship_type"],
            method="refused", method_reason=reason,
            verdict="REFUSED", error_message=reason,
        )

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


def run_shipd_rebuild(
    ships: Sequence[dict],
    *,
    nx: int = DEFAULT_NX,
    nz: int = DEFAULT_NZ,
    out_dir: str | Path = "runs/shipd_rebuild",
) -> dict[str, Any]:
    """Run the full ShipD rebuild pipeline."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    results: list[RebuildResult] = []
    t_start = time.time()

    for ship in ships:
        result = rebuild_ship(ship, nx=nx, nz=nz)
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
                         "realised_L", "realised_B", "realised_T", "rebuild_time_s"])
        for r, s in zip(results, ships):
            writer.writerow([
                r.ship_id, r.ship_type, s["L"], s["B"], s["T"], s["Cb"],
                r.method, r.verdict, r.triangle_count, r.volume_mesh,
                r.volume_parametric, r.volume_error,
                r.realised_L, r.realised_B, r.realised_T, r.rebuild_time_s,
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
                "failed_gates": [g.name for g in r.gates if not g.passed],
            })
    (out_dir / "failure_ledger.json").write_text(json.dumps(failure_ledger, indent=2))

    return evidence
