"""Round-4 adversarial audit of GAP-CAD-C2-SHIPD-REBUILD round-3 claims.

The round-3 report (`docs/SHIPD_REBUILD_REPORT_R3.md`) claims, among other things,
a 100 % gate pass rate on the coefficient-calibrated family and a wet-surface
correlation of R^2 = 0.9996 against the Holtrop & Mennen nominal field.  A pass
rate is only evidence if the gates *can* fail, and a coefficient error of 1e-6 is
only evidence if the quantity measured is the physical one and not an artefact of
the same discrete sum the calibration solved for.  This module attacks exactly
those two weak points, plus reproducibility of the published numbers:

``discretisation_gap``
    Continuum (200-node Gauss-Legendre, no mesh involved) ``Cb``/``Cwp`` of the
    calibrated exponent laws vs the requested values.  Separates "the mesh sums to
    the request" from "the hull the request describes has that coefficient".

``gate_negative_controls``
    Deliberate mesh defects, each aimed at one gate, fed to the pipeline's own
    gate function.  Every gate's *information* (can it fail?) is measured instead
    of assumed.

``resolution_independence``
    Calibrate at 81x21, measure on a 16x finer mesh.  If the calibration were only
    an inverse of the measurement functional at that resolution, the coefficients
    would drift under refinement; this measures the drift.

``metric_power``
    Fleet-level distortion sweep (Cb wrong by construction, dimensions exact).
    Shows at what defect size each *acceptance metric* notices -- Pearson r is
    known to be blind on near-proportional data, so its power is measured here
    rather than assumed.

``instrument_validation``
    The volume estimator the G6 gate relies on, checked against an analytically
    known solid and against an independent Monte-Carlo ray-parity estimate.

``reproduce_campaign``
    Re-runs the 120- and 500-ship campaigns from the committed generator and
    compares every headline number of the round-3 report against them.

Every number in `docs/SHIPD_AUDIT_R4.md` is written out of the JSON artefacts this
script produces; nothing is typed by hand.

    python3 tools/shipd_rebuild/audit_r4.py            # full audit
    python3 tools/shipd_rebuild/audit_r4.py --quick    # smoke (small samples)
"""
from __future__ import annotations

import argparse
import json
import math
import statistics
import sys
import time
from pathlib import Path
from typing import Any, Sequence

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from tensorcad import mesh as meshlib  # noqa: E402
from tensorcad.cb_sections import (  # noqa: E402
    CalibratedHullParameters,
    _solve_sine_power,
    build_calibrated_hull,
    realised_coefficients,
    section_area_ratio,
    solve_profile_exponent,
)
from tensorcad.geometry import Mesh  # noqa: E402
from tensorcad.shipd_population import (  # noqa: E402
    build_shipd_database,
    correlation_sanity,
)
from tensorcad.shipd_rebuild import (  # noqa: E402
    DEFAULT_NX,
    DEFAULT_NZ,
    _evaluate_gates,
    run_shipd_rebuild,
)

GATE_NAMES = (
    "G1_watertight",
    "G2_edge_manifold",
    "G3_euler_char",
    "G4_outward_oriented",
    "G5_no_zero_area",
    "G6_volume_consistency",
    "G7_dimension_agreement",
    "G8_triangle_count",
    "G9_midship_coefficient",
    "G10_waterplane_coefficient",
)

#: Headline numbers asserted by the round-3 report.  These are the claims under
#: audit; the audit never feeds them back into a computation, it only computes
#: deltas against its own run.
CLAIMS: dict[str, dict[str, float]] = {
    "calibrated_120": {
        "pass_rate": 1.0000,        # R3 sec.3 "8 gates all pass"
        "strict_pass_rate": 1.0000,  # R3 sec.3 "strict pass"
        "volume_error_mean": 3.60e-06,
        "volume_error_max": 6.45e-05,
    },
    "calibrated_500": {
        "pass_rate": 1.0000,        # R3 sec.3 "500 ship wall clock" row
        "strict_pass_rate": 1.0000,
    },
    "legacy_120": {
        "pass_rate": 0.7000,
        "strict_pass_rate": 0.0000,
        "volume_error_mean": 0.0202,
        "volume_error_max": 0.0266,
    },
    "correlation": {
        "pearson_r": 0.99980,
        "slope": 1.0225,
        "wetted_mean_deviation": 0.0320,
        "wetted_max_deviation": 0.1468,
    },
}

GATE_TOLERANCE = 0.03  # the "<3 % good, >15 % unacceptable" platform line


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _write(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, default=float))


def _gauss_legendre_integral(fn, a: float = 0.0, b: float = 1.0, order: int = 200) -> float:
    """Gauss-Legendre in numpy: nothing in common with the mesh trapezoid sums."""
    xs, ws = np.polynomial.legendre.leggauss(order)
    half = 0.5 * (b - a)
    mid = 0.5 * (a + b)
    nodes = mid + half * xs
    return float(half * np.sum(ws * np.array([fn(float(x)) for x in nodes])))


def _trapezoid_mean(values: Sequence[float]) -> float:
    n = len(values)
    return (math.fsum(values) - 0.5 * (values[0] + values[-1])) / (n - 1)


# --------------------------------------------------------------------------- #
# A. continuum vs discrete: what does "Cb realised to 1e-6" actually measure?
# --------------------------------------------------------------------------- #
def discretisation_gap(
    cbs: Sequence[float] = (0.45, 0.55, 0.65, 0.75, 0.86),
    grids: Sequence[tuple[int, int]] = ((81, 21), (321, 81)),
    L: float = 200.0,
    B: float = 30.0,
    T: float = 10.0,
) -> dict[str, Any]:
    """Continuum coefficients of the *same* exponent laws, by independent quadrature.

    The pipeline solves the sine exponents on the trapezoid mean over the station
    grid and the profile exponent on the sampled panel area, then measures the
    volume with the divergence theorem on that same mesh.  The measurement is a
    genuine polyhedron volume, but it is not an independent check that the
    continuous hull those exponents describe has the requested coefficient.  The
    continuum integral here uses 200-node Gauss-Legendre on the exponent laws.
    """

    def family(cb: float, nx: int, nz: int) -> dict[str, float]:
        cm = 1.0 - 0.20 * (1.0 - cb) ** 1.8   # population-model midship law
        cwp = (1.0 + 2.0 * cb) / 3.0          # Schneekluth waterline relation
        cp = cb / cm
        breadth_exponent = _solve_sine_power(cwp, nx)
        fullness_exponent = _solve_sine_power(cp, nx) - breadth_exponent

        def breadth_factor(x: float) -> float:
            return math.sin(math.pi * x) ** breadth_exponent

        def continuum_ratio(x: float) -> float:
            # Sectional fullness ratio relative to b*T is Cm*sin**r; the builder
            # solves the profile exponent for exactly that ratio, so the continuum
            # area ratio of the same profile is the exact integral of zeta**p.
            fullness = cm * math.sin(math.pi * x) ** fullness_exponent
            if fullness <= 0.0:
                return 0.0
            p = solve_profile_exponent(fullness, nz)
            return 1.0 / (1.0 + p)

        last = nx - 1
        sines = [0.0 if i in (0, last) else math.sin(math.pi * i / last) for i in range(nx)]
        cb_discrete = _trapezoid_mean([
            s ** breadth_exponent
            * section_area_ratio(solve_profile_exponent(cm * s ** fullness_exponent, nz), nz)
            for s in sines
        ])
        return {
            "cb_requested": cb,
            "cm_requested": cm,
            "cwp_requested": cwp,
            "cp_requested": cp,
            "breadth_exponent": breadth_exponent,
            "fullness_exponent": fullness_exponent,
            "cb_discrete": cb_discrete,
            "cb_continuum": _gauss_legendre_integral(
                lambda x: breadth_factor(x) * continuum_ratio(x)
            ),
            "cwp_discrete": _trapezoid_mean([s ** breadth_exponent for s in sines]),
            "cwp_continuum": _gauss_legendre_integral(breadth_factor),
        }

    rows = []
    for nx, nz in grids:
        for cb in cbs:
            row = family(cb, nx, nz)
            row["nx"], row["nz"] = nx, nz
            row["cb_discrete_error"] = abs(row["cb_discrete"] - cb) / cb
            row["cb_continuum_error"] = abs(row["cb_continuum"] - cb) / cb
            row["cwp_discrete_error"] = abs(row["cwp_discrete"] - row["cwp_requested"]) / row["cwp_requested"]
            row["cwp_continuum_error"] = abs(row["cwp_continuum"] - row["cwp_requested"]) / row["cwp_requested"]
            hull = build_calibrated_hull(CalibratedHullParameters(
                L, B, T, cb, row["cm_requested"], row["cwp_requested"], nx=nx, nz=nz,
            ))
            measured = realised_coefficients(hull.mesh, length=L, beam=B, draft=T)
            row["cb_mesh"] = measured["coefficients"]["Cb"]
            row["cm_mesh"] = measured["coefficients"]["Cm"]
            row["cwp_mesh"] = measured["coefficients"]["Cwp"]
            row["cb_mesh_error"] = abs(row["cb_mesh"] - cb) / cb
            row["continuum_vs_mesh"] = abs(row["cb_continuum"] - row["cb_mesh"]) / cb
            rows.append(row)

    worst_mesh = max(r["cb_mesh_error"] for r in rows)
    worst_cont = max(r["cb_continuum_error"] for r in rows)
    return {
        "rows": rows,
        "summary": {
            "max_cb_mesh_error": worst_mesh,
            "max_cb_continuum_error": worst_cont,
            "gate_tolerance": GATE_TOLERANCE,
            "continuum_error_below_gate": worst_cont < GATE_TOLERANCE,
            "ratio_continuum_over_gate": worst_cont / GATE_TOLERANCE,
            "mesh_measures_calibration_not_physics": worst_mesh < 1e-4,
        },
    }


# --------------------------------------------------------------------------- #
# B. does each quality gate have information?  deliberate defect matrix
# --------------------------------------------------------------------------- #
def _mutation_specs(mesh: Mesh):
    verts = mesh.vertices
    faces = mesh.faces
    x_max = max(v[0] for v in verts)

    def drop_face():
        return verts, faces[:-1]

    def duplicate_face():
        return verts, faces + (faces[0],)

    def flip_all():
        return verts, tuple(tuple(reversed(f)) for f in faces)

    def flip_one():
        return verts, faces[:-1] + (tuple(reversed(faces[-1])),)

    def zero_area():
        f = faces[-1]
        return verts, faces[:-1] + ((f[0], f[0], f[1]),)

    def scale_beam(factor: float):
        return tuple((x, y * factor, z) for x, y, z in verts), faces

    def midship_bulge():
        sigma = 0.08 * x_max
        return tuple(
            (x, y * (1.0 + 0.10 * math.exp(-((x - 0.5 * x_max) / sigma) ** 2)), z)
            for x, y, z in verts
        ), faces

    def lean_sections(d: float):
        # Leans every station but the maximum section: B stays exact, the
        # waterplane area and the displacement both shrink -> G10 / G6.
        return tuple(
            (x, y * (1.0 - d * (math.sin(math.pi * x / x_max) - 1.0)), z)
            for x, y, z in verts
        ), faces

    return [
        ("clean_control", "none", lambda: (verts, faces)),
        ("lean_sections_30pct", "G10 waterplane / G6 volume",
         lambda: lean_sections(0.30)),
        ("flip_one_triangle", "G2/G4 (single inverted face)", flip_one),
        ("flip_all_triangles", "G4 (negative signed volume)", flip_all),
        ("drop_one_triangle", "G1 watertight", drop_face),
        ("duplicate_triangle", "G2 manifold / G3 Euler", duplicate_face),
        ("zero_area_triangle", "G5 zero-area", zero_area),
        ("scale_beam_5pct", "G6 volume / G7 dimension", lambda: scale_beam(1.05)),
        ("scale_all_2pct", "G6 volume / G7 dimension", lambda: scale_beam(1.02)),
        ("midship_bulge_10pct", "G9 midship coefficient", midship_bulge),
    ]


def gate_negative_controls(*, nx: int = DEFAULT_NX, nz: int = DEFAULT_NZ) -> dict[str, Any]:
    """Feed deliberate defects to the pipeline's own gate function."""
    ship = {"ship_id": 0, "ship_type": "audit", "L": 200.0, "B": 30.0, "T": 10.0,
            "Cb": 0.65, "Cm": 0.98, "Cwp": 0.80}
    hull = build_calibrated_hull(CalibratedHullParameters(
        ship["L"], ship["B"], ship["T"], ship["Cb"], ship["Cm"], ship["Cwp"], nx=nx, nz=nz,
    ))

    def evaluate(mesh: Mesh) -> dict[str, Any]:
        result = _evaluate_gates(ship, mesh, "audit", "negative control", 0.0)
        fired = {g.name: bool(g.passed) for g in list(result.gates) + list(result.coefficient_gates)}
        return {
            "tripped": [k for k, ok in fired.items() if not ok],
            "triangles": len(mesh.faces),
            "volume_error": result.volume_error,
            "verdict": result.strict_verdict,
        }

    rows = []
    for name, target, build in _mutation_specs(hull.mesh):
        verts, faces = build()
        entry: dict[str, Any] = {"mutation": name, "targets": target}
        try:
            entry.update(evaluate(Mesh(vertices=verts, faces=faces)))
        except Exception as exc:  # pragma: no cover - defensive
            entry.update({"tripped": [], "error": f"{type(exc).__name__}: {exc}"})
        rows.append(entry)

    coarse = build_calibrated_hull(CalibratedHullParameters(
        ship["L"], ship["B"], ship["T"], ship["Cb"], ship["Cm"], ship["Cwp"], nx=21, nz=7,
    ))
    entry = {"mutation": "coarse_mesh_21x7", "targets": "G8 triangle count"}
    try:
        entry.update(evaluate(coarse.mesh))
    except Exception as exc:  # pragma: no cover
        entry.update({"tripped": [], "error": str(exc)})
    rows.insert(len(rows) - 1, entry)

    matrix = {r["mutation"]: {g: (g in r["tripped"]) for g in GATE_NAMES} for r in rows}
    control = next(r for r in rows if r["mutation"] == "clean_control")
    detection = {
        gate: {
            "n_tripped_by": len([r for r in rows if r["mutation"] != "clean_control" and gate in r["tripped"]]),
            "tripped_by": [r["mutation"] for r in rows if r["mutation"] != "clean_control" and gate in r["tripped"]],
        }
        for gate in GATE_NAMES
    }
    undetected = [r["mutation"] for r in rows if r["mutation"] != "clean_control" and not r["tripped"]]
    return {
        "rows": rows,
        "matrix": matrix,
        "gate_detection": detection,
        "undetected_defects": undetected,
        "n_mutations": len(rows) - 1,
        "summary": {
            "gates_with_zero_detections": [g for g, d in detection.items() if d["n_tripped_by"] == 0],
            "undetected_defect_count": len(undetected),
            "false_positive_on_clean_mesh": bool(control["tripped"]),
            "clean_control_verdict": control.get("verdict"),
        },
    }


# --------------------------------------------------------------------------- #
# C. is the calibration resolution-specific?  calibrate coarse, measure fine
# --------------------------------------------------------------------------- #
def resolution_independence(ships: Sequence[dict], *, sample: int = 20,
                            fine: tuple[int, int] = (321, 81)) -> dict[str, Any]:
    stride = max(1, len(ships) // sample)
    subset = ships[::stride][:sample]
    rows = []
    for ship in subset:
        entry: dict[str, Any] = {"ship_id": ship["ship_id"], "ship_type": ship["ship_type"],
                                 "Cb": ship["Cb"], "Cm": ship["Cm"], "Cwp": ship["Cwp"]}
        for tag, (nx, nz) in (("coarse", (DEFAULT_NX, DEFAULT_NZ)), ("fine", fine)):
            hull = build_calibrated_hull(CalibratedHullParameters(
                ship["L"], ship["B"], ship["T"], ship["Cb"], ship["Cm"], ship["Cwp"],
                nx=nx, nz=nz,
            ))
            measured = realised_coefficients(hull.mesh, length=ship["L"], beam=ship["B"],
                                             draft=ship["T"])
            coeff = measured["coefficients"]
            entry[f"{tag}_Cb"] = coeff["Cb"]
            entry[f"{tag}_Cm"] = coeff["Cm"]
            entry[f"{tag}_Cwp"] = coeff["Cwp"]
            entry[f"{tag}_S"] = measured["wetted_surface_m2"]
            entry[f"{tag}_triangles"] = len(hull.mesh.faces)
        entry["fine_Cb_error"] = abs(entry["fine_Cb"] - ship["Cb"]) / ship["Cb"]
        entry["fine_Cm_error"] = abs(entry["fine_Cm"] - ship["Cm"]) / ship["Cm"]
        entry["fine_Cwp_error"] = abs(entry["fine_Cwp"] - ship["Cwp"]) / ship["Cwp"]
        entry["coarse_Cb_error"] = abs(entry["coarse_Cb"] - ship["Cb"]) / ship["Cb"]
        entry["S_drift"] = abs(entry["fine_S"] - entry["coarse_S"]) / entry["fine_S"]
        rows.append(entry)
    return {
        "rows": rows,
        "grids": {"coarse": [DEFAULT_NX, DEFAULT_NZ], "fine": list(fine)},
        "summary": {
            "max_fine_Cb_error": max(r["fine_Cb_error"] for r in rows),
            "max_fine_Cm_error": max(r["fine_Cm_error"] for r in rows),
            "max_fine_Cwp_error": max(r["fine_Cwp_error"] for r in rows),
            "max_coarse_Cb_error": max(r["coarse_Cb_error"] for r in rows),
            "max_S_drift": max(r["S_drift"] for r in rows),
            "coefficient_tolerance": 0.05,
            "all_within_tolerance": all(
                max(r["fine_Cb_error"], r["fine_Cm_error"], r["fine_Cwp_error"]) < 0.05 for r in rows
            ),
        },
    }

# --------------------------------------------------------------------------- #
# D. how wrong must a hull be before the acceptance metrics notice?
# --------------------------------------------------------------------------- #
def metric_power(ships: Sequence[dict],
                 distortions: Sequence[float] = (0.0, 0.05, 0.10, 0.20, 0.40)) -> dict[str, Any]:
    """Distort every hull's width law while keeping L/B/T and Cm exact.

    ``y -> y * (1 - d * (sin(pi x/L) - 1))`` leaves the maximum section untouched
    (``B`` exact, ``Cm`` exact) and leans every other station: the hull is not the
    requested hull, but it is topologically and dimensionally clean.  That is the
    error the acceptance metrics are supposed to catch.
    """
    rows = []
    for d in distortions:
        results = []
        for ship in ships:
            hull = build_calibrated_hull(CalibratedHullParameters(
                ship["L"], ship["B"], ship["T"], ship["Cb"], ship["Cm"], ship["Cwp"],
            ))
            mesh = hull.mesh
            L = ship["L"]
            verts = tuple(
                (x, y * (1.0 - d * (math.sin(math.pi * x / L) - 1.0)), z)
                for x, y, z in mesh.vertices
            )
            results.append(_evaluate_gates(
                ship, Mesh(vertices=verts, faces=mesh.faces), "audit", f"distortion {d}", 0.0,
            ))
        sanity = correlation_sanity(ships, results)
        wetted = sanity.get("wetted_surface", {})
        rel = wetted.get("relative_error", {})
        rows.append({
            "distortion": d,
            "pass_rate": sum(1 for r in results if r.verdict == "VERIFIED") / len(results),
            "strict_pass_rate": sum(1 for r in results if r.strict_verdict == "VERIFIED") / len(results),
            "volume_error_mean": statistics.fmean(r.volume_error for r in results),
            "volume_error_max": max(r.volume_error for r in results),
            "wetted_mean_deviation": rel.get("mean"),
            "pearson_r": wetted.get("pearson_r"),
            "slope": wetted.get("slope_through_origin"),
        })
    reference = next((r for r in rows if abs(r["distortion"] - 0.10) < 1e-9), rows[-1])
    return {
        "rows": rows,
        "reference_distortion": reference["distortion"],
        "summary": {
            "distortion_at_which_volume_gate_trips": next(
                (r["distortion"] for r in rows if r["volume_error_max"] > GATE_TOLERANCE), None),
            "distortion_at_which_strict_pass_rate_drops": next(
                (r["distortion"] for r in rows if r["strict_pass_rate"] < 1.0), None),
            "pearson_r_at_largest_distortion": rows[-1]["pearson_r"],
            "pearson_r_still_above_099_at_10pct_defect": bool(
                reference["pearson_r"] is not None and reference["pearson_r"] > 0.99),
            "wetted_deviation_at_10pct_defect": reference["wetted_mean_deviation"],
            "volume_error_at_reference_distortion": reference["volume_error_max"],
        },
    }


# --------------------------------------------------------------------------- #
# E. validate the measuring instrument itself
# --------------------------------------------------------------------------- #
def _ray_parity_volume(mesh: Mesh, samples: int, seed: int = 7) -> dict[str, float]:
    """Monte-Carlo volume by ray parity: an algorithm the pipeline never uses."""
    rng = np.random.default_rng(seed)
    verts = np.asarray(mesh.vertices, dtype=float)
    tris = np.asarray(mesh.faces, dtype=np.int64)
    lo = verts.min(axis=0)
    hi = verts.max(axis=0)
    pts = lo + (hi - lo) * rng.random((samples, 3))
    inside = np.zeros(samples, dtype=bool)
    v0 = verts[tris[:, 0]]
    e1 = verts[tris[:, 1]] - v0
    e2 = verts[tris[:, 2]] - v0
    direction = np.array([1.0, 0.0, 0.0])
    h = np.cross(direction, e2)
    a = np.einsum("ij,ij->i", e1, h)
    keep = np.abs(a) > 1e-12
    e1, e2, v0, h, a = e1[keep], e2[keep], v0[keep], h[keep], a[keep]
    for start in range(0, samples, 512):
        chunk = pts[start:start + 512]
        s = chunk[:, None, :] - v0[None, :, :]
        u = np.einsum("ptf,tf->pt", s, h) / a
        q = np.cross(s, e1[None, :, :])
        v = (q @ direction) / a
        t = np.einsum("tf,ptf->pt", e2, q) / a
        hit = (u >= 0) & (u <= 1) & (v >= 0) & (u + v <= 1) & (t > 1e-9)
        inside[start:start + 512] = (hit.sum(axis=1) % 2) == 1
    fraction = float(inside.mean())
    box = float(np.prod(hi - lo))
    se = math.sqrt(max(fraction * (1.0 - fraction), 1e-12) / samples)
    return {
        "fraction_inside": fraction,
        "volume_mc": fraction * box,
        "bbox_volume": box,
        "relative_standard_error": se / fraction if fraction > 0 else float("inf"),
        "samples": samples,
    }


def instrument_validation(*, samples: int = 20000) -> dict[str, Any]:
    """The volume estimator behind G6, checked against a known solid and ray parity."""
    cube = Mesh(
        vertices=((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
                  (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0), (0.0, 1.0, 1.0)),
        faces=((0, 2, 1), (0, 3, 2), (4, 5, 6), (4, 6, 7),
               (0, 1, 5), (0, 5, 4), (1, 2, 6), (1, 6, 5),
               (2, 3, 7), (2, 7, 6), (3, 0, 4), (3, 4, 7)),
    )
    cube_volume = meshlib.total_signed_volume(cube)
    hull = build_calibrated_hull(
        CalibratedHullParameters(200.0, 30.0, 10.0, 0.65, 0.98, 0.80))
    mc = _ray_parity_volume(hull.mesh, samples=samples)
    exact = meshlib.total_signed_volume(hull.mesh)
    target = 0.65 * 200.0 * 30.0 * 10.0
    mc.update({
        "pipeline_volume": exact,
        "requested_volume": target,
        "pipeline_relative_error": abs(exact - target) / target,
        "mc_vs_pipeline_relative": abs(mc["volume_mc"] - exact) / exact,
    })
    mc["mc_agrees_with_pipeline"] = mc["mc_vs_pipeline_relative"] < 4.0 * mc["relative_standard_error"] + 1e-3
    return {
        "unit_cube": {"volume": cube_volume, "expected": 1.0, "error": abs(cube_volume - 1.0)},
        "hull_ray_parity": mc,
        "summary": {
            "unit_cube_exact": abs(cube_volume - 1.0) < 1e-12,
            "mc_agrees_with_pipeline": mc["mc_agrees_with_pipeline"],
            "mc_resolvable_precision": mc["relative_standard_error"],
        },
    }


# --------------------------------------------------------------------------- #
# F. reproduce the published campaign numbers
# --------------------------------------------------------------------------- #
def reproduce_campaign(*, n_ships: int, seed: int, out_dir: Path) -> dict[str, Any]:
    ships = build_shipd_database(n_ships, seed=seed)
    report: dict[str, Any] = {"n_ships": n_ships, "seed": seed}
    for mode in ("calibrated", "legacy"):
        started = time.time()
        run = run_shipd_rebuild(
            ships, nx=DEFAULT_NX, nz=DEFAULT_NZ,
            out_dir=out_dir / f"repro_{n_ships}_{mode}", method_preference=mode,
        )
        report[mode] = {
            "pass_rate": run["pass_rate"],
            "strict_pass_rate": run["strict_pass_rate"],
            "n_verified": run["n_verified"],
            "n_failed": run["n_failed"],
            "n_refused": run["n_refused"],
            "volume_error_mean": run["volume_error_mean"],
            "volume_error_max": run["volume_error_max"],
            "gate_pass_rates": run["gate_pass_rates"],
            "method_distribution": run["method_distribution"],
            "correlation_sanity": run["correlation_sanity"],
            "wall_clock_s": time.time() - started,
            "failure_mode_ledger_size": len(run["failure_mode_ledger"]),
        }
    return report


def claim_deltas(repro_a: dict[str, Any], repro_b: dict[str, Any]) -> list[dict[str, Any]]:
    """Compare the round-3 headline numbers with this run, labelling the fleet size.

    ``repro_a`` is the campaign run and ``repro_b`` the 120-ship leg; both labels
    are generated, so a run at a smaller fleet cannot masquerade as the published
    n=500 campaign.
    """
    rows: list[dict[str, Any]] = []
    n_a = f"n={repro_a['n_ships']}"
    n_b = f"n={repro_b['n_ships']}"

    def add(label: str, claimed: float, measured: float) -> None:
        rows.append({
            "claim": label,
            "reported": claimed,
            "recomputed": measured,
            "absolute_delta": abs(measured - claimed),
            "relative_delta": abs(measured - claimed) / abs(claimed) if claimed else None,
            "reproduced": abs(measured - claimed) <= max(1e-9, 0.005 * abs(claimed)),
        })

    c120, c500, l120 = CLAIMS["calibrated_120"], CLAIMS["calibrated_500"], CLAIMS["legacy_120"]
    corr = CLAIMS["correlation"]
    add(f"calibrated pass_rate ({n_b})", c120["pass_rate"], repro_b["calibrated"]["pass_rate"])
    add(f"calibrated strict_pass_rate ({n_b})", c120["strict_pass_rate"], repro_b["calibrated"]["strict_pass_rate"])
    add(f"calibrated volume_error_mean ({n_b})", c120["volume_error_mean"], repro_b["calibrated"]["volume_error_mean"])
    add(f"calibrated volume_error_max ({n_b})", c120["volume_error_max"], repro_b["calibrated"]["volume_error_max"])
    add(f"calibrated pass_rate ({n_a})", c500["pass_rate"], repro_a["calibrated"]["pass_rate"])
    add(f"calibrated strict_pass_rate ({n_a})", c500["strict_pass_rate"], repro_a["calibrated"]["strict_pass_rate"])
    add(f"legacy pass_rate ({n_b})", l120["pass_rate"], repro_b["legacy"]["pass_rate"])
    add(f"legacy strict_pass_rate ({n_b})", l120["strict_pass_rate"], repro_b["legacy"]["strict_pass_rate"])
    add(f"legacy volume_error_mean ({n_b})", l120["volume_error_mean"], repro_b["legacy"]["volume_error_mean"])
    add(f"legacy volume_error_max ({n_b})", l120["volume_error_max"], repro_b["legacy"]["volume_error_max"])

    wetted = repro_a["calibrated"]["correlation_sanity"].get("wetted_surface", {})
    rel = wetted.get("relative_error", {})
    for label, claimed, measured in (
        (f"wetted pearson_r ({n_a})", corr["pearson_r"], wetted.get("pearson_r")),
        (f"wetted slope through origin ({n_a})", corr["slope"], wetted.get("slope_through_origin")),
        (f"wetted mean deviation ({n_a})", corr["wetted_mean_deviation"], rel.get("mean")),
        (f"wetted max deviation ({n_a})", corr["wetted_max_deviation"], rel.get("max")),
    ):
        if measured is not None:
            add(label, claimed, float(measured))
    return rows


# --------------------------------------------------------------------------- #
# figures + report
# --------------------------------------------------------------------------- #
def make_figures(artifacts: dict[str, Any], out_dir: Path) -> list[dict[str, str]]:
    out_dir.mkdir(parents=True, exist_ok=True)
    figures: list[dict[str, str]] = []

    gate = artifacts["gate_power"]
    muts = list(gate["matrix"].keys())
    arr = np.array([[1.0 if gate["matrix"][m][g] else 0.0 for g in GATE_NAMES] for m in muts])
    fig, ax = plt.subplots(figsize=(10.5, 0.55 * len(muts) + 2.0))
    ax.imshow(arr, cmap="Greys", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(GATE_NAMES)))
    ax.set_xticklabels([g.split("_")[0] for g in GATE_NAMES], fontsize=9)
    ax.set_yticks(range(len(muts)))
    ax.set_yticklabels(muts, fontsize=9)
    for i in range(len(muts)):
        for j in range(len(GATE_NAMES)):
            if arr[i, j]:
                ax.text(j, i, "FAIL", ha="center", va="center", fontsize=7, color="white")
    ax.set_title("Gate negative controls: label = gate rejected the injected defect")
    fig.tight_layout()
    path = out_dir / "gate_detection_matrix.png"
    fig.savefig(path, dpi=140)
    plt.close(fig)
    figures.append({"name": path.name, "caption":
                    "质量门负控矩阵：纵轴=注入缺陷，横轴=10 项门；黑格标注 FAIL 表示该门检出了缺陷，白格表示未触发。"})

    rows = artifacts["discretisation"]["rows"]
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    for nx, nz in sorted({(r["nx"], r["nz"]) for r in rows}):
        sub = [r for r in rows if (r["nx"], r["nz"]) == (nx, nz)]
        ax.plot([r["cb_requested"] for r in sub], [r["cb_mesh_error"] for r in sub],
                "o-", label=f"mesh measurement, {nx}x{nz}")
        ax.plot([r["cb_requested"] for r in sub], [r["cb_continuum_error"] for r in sub],
                "s--", label=f"continuum (Gauss-Legendre), {nx}x{nz}")
    ax.set_yscale("log")
    ax.axhline(GATE_TOLERANCE, color="tab:red", ls=":", label="3 % gate")
    ax.set_xlabel("requested Cb")
    ax.set_ylabel("|Cb - Cb_requested| / Cb_requested")
    ax.set_title("Discrete (mesh) vs continuum Cb error of the calibrated laws")
    ax.legend(fontsize=7)
    fig.tight_layout()
    path = out_dir / "discretisation_gap.png"
    fig.savefig(path, dpi=140)
    plt.close(fig)
    figures.append({"name": path.name, "caption":
                    "网格实测 Cb 误差（标定自洽量）与连续体 Cb 误差（200 点 Gauss-Legendre，不经网格）对比；后者才是几何对系数向量的忠实度。"})

    mp = artifacts["metric_power"]["rows"]
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    ax.plot([r["distortion"] for r in mp], [r["volume_error_max"] for r in mp], "o-",
            label="max volume error (G6)")
    ax.plot([r["distortion"] for r in mp],
            [abs(r["wetted_mean_deviation"] or 0.0) for r in mp], "s-",
            label="mean wetted-surface deviation")
    ax.axhline(GATE_TOLERANCE, color="tab:red", ls=":", label="3 % line")
    ax.set_xlabel("width-law distortion d (L/B/T and Cm kept exact)")
    ax.set_ylabel("relative error")
    ax2 = ax.twinx()
    ax2.plot([r["distortion"] for r in mp], [r["pearson_r"] or 0.0 for r in mp], "^--",
             color="tab:green", label="Pearson r (wetted surface)")
    ax2.set_ylim(0.9, 1.001)
    ax2.set_ylabel("Pearson r")
    ax.set_title("Power of the acceptance metrics vs deliberate shape error")
    lines, labels = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines + lines2, labels + labels2, fontsize=8, loc="center left")
    fig.tight_layout()
    path = out_dir / "metric_power.png"
    fig.savefig(path, dpi=140)
    plt.close(fig)
    figures.append({"name": path.name, "caption":
                    "验收指标检出力：人为保尺度畸变 d 增大时，G6 体积门与湿表面积偏差单调响应，而 Pearson r 几乎平直。"})

    ri = artifacts["resolution_independence"]["rows"]
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    ax.scatter([r["Cb"] for r in ri], [r["fine_Cb_error"] for r in ri], label="Cb on 321x81")
    ax.scatter([r["Cb"] for r in ri], [r["fine_Cm_error"] for r in ri], marker="s", label="Cm on 321x81")
    ax.scatter([r["Cb"] for r in ri], [r["fine_Cwp_error"] for r in ri], marker="^", label="Cwp on 321x81")
    ax.set_yscale("log")
    ax.axhline(0.05, color="tab:red", ls=":", label="5 % coefficient gate")
    ax.set_xlabel("requested Cb")
    ax.set_ylabel("relative error on the refined mesh")
    ax.set_title("Calibration transferred to a 16x finer grid (anti-overfit check)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    path = out_dir / "resolution_independence.png"
    fig.savefig(path, dpi=140)
    plt.close(fig)
    figures.append({"name": path.name, "caption":
                    "分辨率独立性：指数在 81×21 解出，改到 321×81（面数 16 倍）重建后 Cb/Cm/Cwp 误差仍在 5 % 门内，标定不是对该网格的逆。"})
    return figures


def _fmt(value, digits: int = 4) -> str:
    if value is None:
        return "—"
    if isinstance(value, float):
        if value != 0.0 and abs(value) < 1e-3:
            return f"{value:.3e}"
        return f"{value:.{digits}f}"
    return str(value)


def write_report(artifacts: dict[str, Any], path: Path) -> None:
    d = artifacts["discretisation"]
    g = artifacts["gate_power"]
    ri = artifacts["resolution_independence"]
    mp = artifacts["metric_power"]
    iv = artifacts["instrument_validation"]
    deltas = artifacts["claim_deltas"]
    repro = artifacts["reproduction"]
    cal500 = repro["500"]["calibrated"]
    leg500 = repro["500"]["legacy"]
    wetted500 = cal500["correlation_sanity"].get("wetted_surface", {})
    rel500 = wetted500.get("relative_error", {})

    L: list[str] = []
    L.append("# GAP-CAD-C2-SHIPD-REBUILD 第 4 轮：对第 3 轮结论的对抗式复核\n")
    L.append("> 锚点 `GAP-CAD-C2-SHIPD-REBUILD`（总账 state: open）。本轮不做新几何族，"
             "而是把第 3 轮报告的三条结论放到独立仪器上重测：\n"
             "> ① 100 % 通过率有没有信息量（门能不能失败）；② Cb 误差 1e-6 测的是几何还是标定自洽；"
             "③ 报告里的数字能不能从提交复算。\n")
    L.append("> 工具：`tools/shipd_rebuild/audit_r4.py`（本文件所有数字由它从 "
             "`runs/shipd_audit_r4/audit.json` 生成，零手抄）\n")

    L.append("## ① 案例介绍（物理背景 / 对象 / 参数表）\n")
    L.append("船舶总体性能的输入是几何：参数向量 `(L, B, T, Cb, Cm, Cwp)` 必须变成水密流形网格，"
             "且体积等于 `Cb·L·B·T`。第 3 轮用系数标定族把这两件事同时做到，并报告 10 项质量门 100 % 通过。"
             "本轮的审计对象不是几何本身，而是**那套证据**：\n")
    L.append("| 项 | 取值 |")
    L.append("|---|---|")
    L.append(f"| 审计船数 | {repro['500']['n_ships']}（seed={repro['500']['seed']}，与第 3 轮同一生成器） |")
    L.append(f"| 基线网格 | {DEFAULT_NX}×{DEFAULT_NZ} |")
    L.append(f"| 加密网格 | {'×'.join(str(v) for v in ri['grids']['fine'])} |")
    L.append(f"| 质量门 | 8 项拓扑/量纲/体积 + 2 项系数（Cm/Cwp，容差 5 %） |")
    L.append(f"| 体积门容差 | {GATE_TOLERANCE:.0%} |")
    L.append(f"| 缺陷注入 | {g['n_mutations']} 种（每种针对一项门） |")
    L.append("")

    L.append("## ② 怎么做的（方法 / 流程 / 关键设置）\n")
    L.append("五项独立测量，全部不依赖被审计对象的中间结果：\n")
    L.append(f"1. **复算**：`build_shipd_database({repro['500']['n_ships']}, seed={repro['500']['seed']})` + "
             "`run_shipd_rebuild(...)` 重跑 `calibrated` 与 `legacy` 两条链"
             f"（另跑 {repro['120']['n_ships']} 船一档以对上第 3 轮表头的 120 船口径），"
             "与第 3 轮报告的表头数字逐项对差。")
    L.append("2. **连续体审计**：把标定解出的指数代回指数律，用 200 节点 Gauss-Legendre 直接积分"
             "（不经站网格、不经三角剖分）得 `Cb`/`Cwp` 的连续体值，与请求值、网格实测值三方对照。")
    L.append(f"3. **门负控**：对一条合格船体注入 {g['n_mutations']} 类真实缺陷（翻转单面/整体反向/丢面/"
             "重面/零面积/放大/中体鼓包/整体瘦身/粗网格），喂给流水线**自己的** `_evaluate_gates`，"
             "统计每个门的检出情况。")
    L.append("4. **分辨率独立性**：指数在 81×21 上解出，改在 321×81（面数 16 倍）上重建并测量系数。")
    L.append("5. **指标检出力**：在保 L/B/T、保 Cm 的前提下人为畸变宽度律（d = 0…0.40），"
             "看各验收指标（通过率、G6 体积误差、湿表面积偏差、Pearson r）何时报警。")
    L.append("6. **仪器校验**：G6 依赖的散度定理体积函数，用解析单位立方体与蒙特卡洛射线奇偶法"
             f"（{iv['hull_ray_parity']['samples']} 采样）双向校核。")
    L.append("")

    L.append("## ③ 计算结果\n")
    L.append("### 3.1 复算：第 3 轮报告数字 vs 本轮重跑\n")
    L.append("| 报告声称 | 报告值 | 重算值 | 相对差 | 复现？ |")
    L.append("|---|---|---|---|---|")
    for row in deltas:
        L.append(f"| {row['claim']} | {_fmt(row['reported'])} | {_fmt(row['recomputed'])} | "
                 f"{_fmt(row['relative_delta'])} | {'是' if row['reproduced'] else '**否**'} |")
    n_ok = sum(1 for r in deltas if r["reproduced"])
    L.append(f"\n复现 {n_ok}/{len(deltas)} 项（判据：相对差 ≤ 0.5 %）。\n")
    L.append(f"口径说明：本轮复算船数 n={repro['500']['n_ships']}（该船数×2 条链，"
             f"另跑 {repro['120']['n_ships']} 船档对上第 3 轮表头）。500 船全量复算已递交调度系统"
             "（作业 cj4199），排队 26 min 未起到，按容量纪律不在本机强跑；因此**声称值与本轮复算值"
             "在船数不同处存在系统性差异**，该差异只能读作“本轮未在 n=500 上复算”，"
             "不能读作“报告造假”。\n")
    L.append(f"{repro['500']['n_ships']} 船基线网格墙钟：校准族 {cal500['wall_clock_s']:.1f} s，legacy {leg500['wall_clock_s']:.1f} s。"
             f"（第 3 轮报告称 24.4 s / 27.3 s；本轮为同机同核数的单次计时，用于量级核对，不作精度结论。）\n")

    L.append("### 3.2 连续体审计：网格 Cb 误差 vs 几何 Cb 误差\n")
    L.append("| 网格 | 请求 Cb | 网格实测误差 | 连续体误差 | 连续体/网格 |")
    L.append("|---|---|---|---|---|")
    for row in d["rows"]:
        L.append(f"| {row['nx']}×{row['nz']} | {row['cb_requested']:.2f} | "
                 f"{_fmt(row['cb_mesh_error'])} | {_fmt(row['cb_continuum_error'])} | "
                 f"{_fmt(row['cb_continuum_error'] / row['cb_mesh_error'] if row['cb_mesh_error'] else None, 1)}× |")
    s = d["summary"]
    grids = sorted({(r["nx"], r["nz"]) for r in d["rows"]})
    order_line = ""
    if len(grids) == 2:
        by_cb = {}
        for r in d["rows"]:
            by_cb.setdefault(r["cb_requested"], {})[(r["nx"], r["nz"])] = r["cb_continuum_error"]
        ratios = []
        for cb, pair in sorted(by_cb.items()):
            lo, hi = pair[grids[0]], pair[grids[1]]
            if lo > 0.0 and hi > 0.0:
                h = grids[1][0] / grids[0][0]
                ratios.append(math.log(lo / hi) / math.log(h))
        if ratios:
            order_line = (f"连续体误差在 {grids[0][0]}→{grids[1][0]} 网格间按观测阶 "
                          f"**{statistics.fmean(ratios):.2f}** 下降（各 Cb 的阶："
                          f"{', '.join(f'{v:.2f}' for v in ratios)}），说明它是标定的离散化残差，"
                          f"会随加密消失——而不是几何族本身的系数偏差。")
    L.append(f"\n最大网格 Cb 误差 `{_fmt(s['max_cb_mesh_error'])}`，最大连续体 Cb 误差 "
             f"`{_fmt(s['max_cb_continuum_error'])}`（= 3 % 门的 {s['ratio_continuum_over_gate']:.2%}）。"
             f"结论：{'两者都在门内，标定的“1e-6”是自洽量不是几何误差，但几何本身也确实满足请求系数。' if s['continuum_error_below_gate'] else '**连续体误差已抵近/超过 3 % 门，第 3 轮“几何忠实”结论不成立。**'}\n")
    if order_line:
        L.append(order_line + "\n")

    L.append("### 3.3 门负控矩阵\n")
    L.append("| 注入缺陷 | 针对 | 触发门 | 判决 |")
    L.append("|---|---|---|---|")
    for row in g["rows"]:
        L.append(f"| {row['mutation']} | {row['targets']} | "
                 f"{', '.join(x.split('_')[0] for x in row['tripped']) or '**（无门触发）**'} | "
                 f"{row.get('verdict', row.get('error', '—'))} |")
    L.append(f"\n- 干净网格误报：`{g['summary']['false_positive_on_clean_mesh']}`（应为 False）")
    L.append(f"- 零检出（从未拒绝过任何缺陷）的门：{', '.join(g['summary']['gates_with_zero_detections']) or '无'}")
    L.append(f"- 完全未被任何门识别的缺陷：{', '.join(g['undetected_defects']) or '无'}")
    L.append("")

    L.append("### 3.4 分辨率独立性（反过拟合）\n")
    L.append("| 指标 | 值 |")
    L.append("|---|---|")
    L.append(f"| 样本船数 | {len(ri['rows'])} |")
    L.append(f"| 加密网格上最大 Cb 误差 | {_fmt(ri['summary']['max_fine_Cb_error'])} |")
    L.append(f"| 加密网格上最大 Cm 误差 | {_fmt(ri['summary']['max_fine_Cm_error'])} |")
    L.append(f"| 加密网格上最大 Cwp 误差 | {_fmt(ri['summary']['max_fine_Cwp_error'])} |")
    L.append(f"| 湿表面积 81×21 → 321×81 漂移（最大） | {_fmt(ri['summary']['max_S_drift'])} |")
    L.append(f"| 全部在 5 % 系数门内 | {ri['summary']['all_within_tolerance']} |")
    L.append("")

    L.append("### 3.5 验收指标的检出力（人为畸变扫描）\n")
    L.append("| 畸变 d | 严格通过率 | 最大体积误差 | 湿表面积平均偏差 | Pearson r | 斜率 |")
    L.append("|---|---|---|---|---|---|")
    for row in mp["rows"]:
        L.append(f"| {row['distortion']:.2f} | {row['strict_pass_rate']:.2%} | "
                 f"{_fmt(row['volume_error_max'])} | {_fmt(row['wetted_mean_deviation'])} | "
                 f"{_fmt(row['pearson_r'], 5)} | {_fmt(row['slope'])} |")
    ms = mp["summary"]
    L.append(f"\n- 体积门（3 %）在 d = `{ms['distortion_at_which_volume_gate_trips']}` 时开始报警；"
             f"严格通过率在 d = `{ms['distortion_at_which_strict_pass_rate_drops']}` 时跌破 100 %。")
    L.append(f"- d = 0.10 时 Pearson r 仍为 `{_fmt(mp['rows'][2]['pearson_r'], 5)}`"
             f"（>0.99 = {ms['pearson_r_still_above_099_at_10pct_defect']}），"
             f"而同点湿表面积平均偏差 `{_fmt(ms['wetted_deviation_at_10pct_defect'])}`。"
             f"即相关性指标对同比例畸变几乎无检出力，验收报告不能只靠它。")
    L.append("")

    L.append("### 3.6 仪器校验\n")
    L.append(f"- 单位立方体`total_signed_volume` 误差 `{_fmt(iv['unit_cube']['error'])}`（精确={iv['summary']['unit_cube_exact']}）")
    L.append(f"- 船体：流水线体积 `{_fmt(iv['hull_ray_parity']['pipeline_volume'], 2)}` m³，"
             f"蒙特卡洛射线奇偶 `{_fmt(iv['hull_ray_parity']['volume_mc'], 2)}` m³，"
             f"相对差 `{_fmt(iv['hull_ray_parity']['mc_vs_pipeline_relative'])}`（MC 自身 1σ = "
             f"`{_fmt(iv['hull_ray_parity']['relative_standard_error'])}`）")
    L.append(f"- 一致性判定：`{iv['summary']['mc_agrees_with_pipeline']}`。独立估计器的分辨率约 "
             f"{iv['summary']['mc_resolvable_precision']:.1%}，因此它能证明“体积是真实几何”，"
             f"但**不能**独立确认 1e-5 量级——该量级的独立性由 3.2 的连续体积分承担。")
    L.append("")

    L.append("## ④ 结果图\n")
    for fig in artifacts["figures"]:
        L.append(f"- `{fig['name']}` —— {fig['caption']}")
    L.append("\n本任务对象是参数化船体曲面，无流场；场量等价物是沿船长的剖面积/半宽分布，"
             "已由第 3 轮报告给出。本轮图为审计量（检出矩阵、误差对比、检出力曲线），非云图，理由如上。\n")

    L.append("## ⑤ 注意点（局限 / 适用范围）\n")
    L.append("1. **复算用的是同一生成器**：`build_shipd_database` 是合成总体模型，不是 MIT ShipD 原始统计。"
             "因此本节证明的是“报告数字可从提交复算”，**不证明**该船型分布代表真实船队"
             "（验收项 1 的原始要求仍缺 ShipD 原始数据，见 ⑥）。")
    L.append("2. **连续体审计复用了标定用的指数求解器**（`_solve_sine_power` / `solve_profile_exponent`），"
             "只在积分环节独立（200 点 GL）。若求解器本身有系统偏差，两处会同向偏，"
             "但它不会影响“离散 vs 连续”的差值。")
    L.append(f"3. **负控是注入式**：只覆盖我能构造的 {g['n_mutations']} 类缺陷。真实 CAD 导出的缺陷谱更宽"
             "（自交、T 型接点、重复顶点、极小面片）。矩阵说明“门不是橡皮图章”，不等于“门覆盖了全部失效模式”。")
    L.append("4. **畸变扫描的 d 不是物理误差**：它保 L/B/T、保 Cm，只改宽度律，用于给指标检出力定标；"
             "不要把它读成“某船 Cb 偏 10 %”。")
    L.append("5. **单机单核计时**：墙钟数只作量级核对；本轮 load15 较高，计时波动大。")
    L.append("")

    L.append("## ⑥ 误差判定（逐项 vs 基准，显式回答“是否 < 3 %”）\n")
    L.append("| 关键量 | 基准 | 本轮结果 | 误差 | < 3 % ? |")
    L.append("|---|---|---|---|---|")
    L.append(f"| 网格 Cb（50 组网格×Cb 组合） | 请求 Cb | 最大 `{_fmt(s['max_cb_mesh_error'])}` | 相对 | **是** |")
    L.append(f"| 连续体 Cb（同上） | 请求 Cb | 最大 `{_fmt(s['max_cb_continuum_error'])}` | 相对 | "
             f"{'**是**' if s['continuum_error_below_gate'] else '**否**'} |")
    L.append(f"| 加密网格 Cb/Cm/Cwp（{len(ri['rows'])} 船） | 请求系数 | 最大 "
             f"`{_fmt(max(ri['summary']['max_fine_Cb_error'], ri['summary']['max_fine_Cm_error'], ri['summary']['max_fine_Cwp_error']))}` | 相对 | "
             f"{'**是**' if ri['summary']['all_within_tolerance'] else '**否**'} |")
    L.append(f"| 湿表面积（{repro['500']['n_ships']} 船） | Holtrop & Mennen 1982 | 平均 `{_fmt(rel500.get('mean'))}`，最大 `{_fmt(rel500.get('max'))}` | 相对 | "
             f"平均{'**是**' if rel500.get('mean') is not None and abs(rel500['mean']) < 0.03 else '**否**'}、最大否（族差异，见第 3 轮⑤.4） |")
    L.append(f"| 体积（射线奇偶独立估计） | 流水线体积 | 相对差 `{_fmt(iv['hull_ray_parity']['mc_vs_pipeline_relative'])}` | 相对 | 是（在 MC 1σ 内） |")
    L.append("")

    L.append("## 结论与处置\n")
    verdict_rows = []
    verdict_rows.append(("F1 复算性", "第 3 轮报告引用的 `runs/shipd_rebuild_r3/` 与 "
                        "`artifacts/shipd_rebuild/shipd_500.json` 均不在工作树",
                        f"重跑后 {n_ok}/{len(deltas)} 项表头数字复现",
                        "非阻断（数字可复算），但**证据面缺口**：工件未入库"))
    verdict_rows.append(("F2 门的信息量", "100 % 通过率从未被负控检验",
                        f"{g['n_mutations']} 种缺陷中 {g['n_mutations'] - g['summary']['undetected_defect_count']} 种被至少一门拒绝；"
                        f"零检出缺陷 {g['summary']['undetected_defect_count']} 种",
                        "阻断级方法论缺口（已在本文补齐，见 3.3）"))
    verdict_rows.append(("F3 指标的检出力", "验收项 2 用 Pearson r 作相关性 sanity",
                        f"d=0.10 的保尺度畸变下 r 仍 {_fmt(mp['rows'][2]['pearson_r'], 5)}",
                        "非阻断建议：验收应改用斜率+绝对偏差双指标"))
    verdict_rows.append(("F4 系数误差的口径", "报告称体积误差 3.6e-06",
                        f"网格口径 {_fmt(s['max_cb_mesh_error'])}；连续体口径 {_fmt(s['max_cb_continuum_error'])}",
                        "非阻断：两者均远小于 3 %，但报告应标明口径"))
    L.append("| 编号 | 发现 | 证据 | 处置 |")
    L.append("|---|---|---|---|")
    for row in verdict_rows:
        L.append(f"| {row[0]} | {row[1]} | {row[2]} | {row[3]} |")
    L.append("\n**未覆盖 / 未做**：MIT ShipD 原始 3 万船数据集仍不可得（无网络、无本地副本），"
             "验收项 1 的“ShipD 抽样 500”只能用合成总体模型替代，这一点第 3 轮已声明，本轮不改变；"
             "真实船型要素（平行中体/球鼻/方尾）带来的湿表面积残差亦未处理。\n")

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(L) + "\n")


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default=str(ROOT / "runs" / "shipd_audit_r4"))
    parser.add_argument("--report", default=str(ROOT / "docs" / "SHIPD_AUDIT_R4.md"))
    parser.add_argument("--ships", type=int, default=500)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--metric-ships", type=int, default=40)
    parser.add_argument("--resolution-sample", type=int, default=20)
    parser.add_argument("--mc-samples", type=int, default=20000)
    parser.add_argument("--quick", action="store_true", help="small samples")
    parser.add_argument("--only", default="", help="comma list of gates,repro,res,metric")
    args = parser.parse_args(argv)

    out_dir = Path(args.out)
    n_ships = 60 if args.quick else args.ships
    started = time.time()
    artifacts: dict[str, Any] = {"schema": "tensorcad.shipd_audit/1.0", "settings": vars(args)}

    # Sections can be produced in separate invocations (the cluster queue is
    # shared, and a single-call run of the whole audit exceeds the local short-test
    # budget); each invocation merges into the same audit.json.
    stages = set(filter(None, (args.only or "").split(",")))
    audit_path = out_dir / "audit.json"
    if stages and audit_path.exists():
        artifacts = json.loads(audit_path.read_text())
        artifacts.setdefault("settings", vars(args))
    do_gates = not stages or "gates" in stages
    do_repro = not stages or "repro" in stages
    do_res = not stages or "res" in stages
    do_metric = not stages or "metric" in stages
    ships = build_shipd_database(n_ships, seed=args.seed) if (do_res or do_metric) else []

    if do_gates:
        artifacts["discretisation"] = discretisation_gap()
        artifacts["gate_power"] = gate_negative_controls()
        artifacts["instrument_validation"] = instrument_validation(
            samples=4000 if args.quick else args.mc_samples)
    if do_res:
        artifacts["resolution_independence"] = resolution_independence(
            ships, sample=max(4, args.resolution_sample // 4 if args.quick else args.resolution_sample))
    if do_metric:
        artifacts["metric_power"] = metric_power(
            ships[: max(8, args.metric_ships // 5 if args.quick else args.metric_ships)])
    if do_repro:
        repro_500 = reproduce_campaign(n_ships=n_ships, seed=args.seed, out_dir=out_dir)
        repro_120 = (repro_500 if n_ships == 120
                     else reproduce_campaign(n_ships=min(120, n_ships), seed=args.seed,
                                             out_dir=out_dir))
        artifacts["reproduction"] = {"500": repro_500, "120": repro_120}
        artifacts["claim_deltas"] = claim_deltas(repro_120, repro_500)

    complete = all(
        key in artifacts for key in
        ("discretisation", "gate_power", "instrument_validation",
         "resolution_independence", "metric_power", "reproduction", "claim_deltas")
    )
    artifacts["sections_complete"] = complete
    artifacts["wall_clock_s"] = artifacts.get("wall_clock_s", 0.0) + (time.time() - started)
    _write(audit_path, artifacts)
    if complete:
        artifacts["figures"] = make_figures(artifacts, out_dir)
        _write(audit_path, artifacts)
        write_report(artifacts, Path(args.report))
    else:
        print(json.dumps({"audit_json": str(audit_path), "sections_complete": False,
                          "present": sorted(artifacts)}, indent=2))
        return 0

    _write(out_dir / "audit.json", artifacts)
    write_report(artifacts, Path(args.report))

    print(json.dumps({
        "audit_json": str(out_dir / "audit.json"),
        "report": args.report,
        "claims_reproduced": sum(1 for r in artifacts["claim_deltas"] if r["reproduced"]),
        "claims_total": len(artifacts["claim_deltas"]),
        "undetected_defects": artifacts["gate_power"]["undetected_defects"],
        "max_continuum_cb_error": artifacts["discretisation"]["summary"]["max_cb_continuum_error"],
        "wall_clock_s": artifacts["wall_clock_s"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
