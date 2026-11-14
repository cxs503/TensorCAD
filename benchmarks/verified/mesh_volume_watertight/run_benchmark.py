"""Verified case ``mesh_volume_watertight``: closed-body volume/watertightness vs analytic solids.

Anchor ``CASE-CAD-VOLUME-CLOSED`` (case pool, round 1).  Every reference here is a
closed form evaluated in code (:mod:`tensorcad.mesh_volume`), so nothing is
transcribed and a number can be *wrong* rather than merely *different*:

============================  ==============================================
quantity                      reference
============================  ==============================================
box volume / area             ``a b c`` / ``2(ab + bc + ca)``
sphere volume / area          ``4/3 pi r^3`` / ``4 pi r^2``
icosahedron volume            ``(5/12)(3+sqrt5)a^3``, ``a = 4R/sqrt(10+2sqrt5)``
                              -- the ``sub0`` mesh's own closed form
topology of a closed body     ``open_edges = 0``, ``euler = 2``, no non-manifold
                              edge, no duplicate facet, no winding conflict,
                              ``signed_volume > 0``, audit verdict ``passed``
============================  ==============================================

Thresholds are declared *before* the run and live in the module, not in this
script's argv, so they cannot be relaxed after a number is seen:

* :data:`tensorcad.mesh_volume.EXACT_TOLERANCE` (1e-12) for the box family: a
  triangle mesh represents a box exactly, only round-off is allowed;
* :data:`tensorcad.mesh_volume.RELATIVE_TOLERANCE` (1e-3, i.e. 0.1 %) for the
  curved body at the finest level of the ladder -- the case's acceptance line;
* the ladder itself is *reported*, not gated level by level: the observed order
  must land in ``[1.5, 2.5]`` (second order is what a piecewise-linear
  approximation of a smooth surface gives; the band is a sanity check on the
  fit, not a claim about the constant).

    python benchmarks/verified/mesh_volume_watertight/run_benchmark.py
    python benchmarks/verified/mesh_volume_watertight/run_benchmark.py --out runs/mesh_volume_watertight

Exits non-zero when a gate breaks.  ``result.json``, ``convergence.csv`` and the
two figures are generated here; every number in
``docs/mesh_volume_watertight_bench.md`` is read back from that artifact.
"""
from __future__ import annotations

import argparse
import csv
import importlib
import json
import math
import platform
import sys
import time
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import numpy as np

_REPO_ROOT = Path(__file__).resolve().parents[3]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import tensorcad
from tensorcad import mesh_volume as mv
from tensorcad.geometry_audit import (
    audit_watertightness,
    inject_duplicate_faces,
    inject_flip_faces,
)
from tensorcad.mesh_xval import TriMesh, icosphere, weld, write_stl

#: Box family: (name, a, b, c, rotation_deg, translation)
BoxCase = tuple[str, float, float, float, tuple[float, float, float], tuple[float, float, float]]

BOX_CASES: tuple[BoxCase, ...] = (
    ("box_2x2x2", 2.0, 2.0, 2.0, (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)),
    ("box_1.5x0.75x3", 1.5, 0.75, 3.0, (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)),
    ("box_2x2x2_rot_trans", 2.0, 2.0, 2.0, (31.0, -17.0, 73.0), (13.7, -42.1, 7.3)),
)

#: Sphere ladder, in subdivision levels of the icosphere (faces = 20 * 4**level).
SPHERE_RADIUS = 1.0
SPHERE_LEVELS = tuple(range(0, 7))
SPHERE_GATE_LEVEL = 6
#: Second sphere radius, to show the error is a dimensionless quantity.
SPHERE_SCALE_RADIUS = 2.5
SPHERE_SCALE_LEVEL = 6
#: The rung that is *supposed* to miss the 0.1 % line: if a coarse inscribed
#: polyhedron reported the analytic sphere volume, the measurement would be
#: disconnected from the geometry.
SPHERE_COARSE_LEVEL = 4
#: Box used for the STL round trip: half-sizes that are not exactly representable
#: in float32, so the format's own loss is visible instead of rounding away.
STL_BOX = (1.7, 1.1, 2.9)

ORDER_BAND = (1.5, 2.5)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _rotation(angles_deg: Sequence[float]) -> np.ndarray:
    """Deterministic Z-Y-X rotation matrix from three angles in degrees."""
    ax, ay, az = (math.radians(a) for a in angles_deg)
    rx = np.array([[1, 0, 0], [0, math.cos(ax), -math.sin(ax)], [0, math.sin(ax), math.cos(ax)]])
    ry = np.array([[math.cos(ay), 0, math.sin(ay)], [0, 1, 0], [-math.sin(ay), 0, math.cos(ay)]])
    rz = np.array([[math.cos(az), -math.sin(az), 0], [math.sin(az), math.cos(az), 0], [0, 0, 1]])
    return rz @ ry @ rx


def _transformed(mesh: TriMesh, angles_deg: Sequence[float], offset: Sequence[float]) -> TriMesh:
    points = mesh.points @ _rotation(angles_deg).T + np.asarray(offset, dtype=np.float64)
    return TriMesh(points, mesh.faces, mesh.name + "/rot_trans")


def _third_party_volume(mesh: TriMesh) -> float | None:
    """Volume from ``trimesh`` (an independent implementation), if importable."""
    try:
        trimesh = importlib.import_module("trimesh")
    except Exception:
        return None
    return float(trimesh.Trimesh(vertices=mesh.points, faces=mesh.faces, process=False).volume)


def _row(
    case: str,
    quantity: str,
    reference: float,
    measured: float,
    tolerance: float,
    *,
    absolute: bool = False,
    note: str = "",
) -> dict[str, Any]:
    error = measured - reference if absolute else (measured - reference) / reference
    return {
        "case": case,
        "quantity": quantity,
        "reference": reference,
        "measured": measured,
        "error": error,
        "error_kind": "absolute" if absolute else "relative",
        "tolerance": tolerance,
        "passed": abs(error) <= tolerance,
        "note": note,
    }


# ---------------------------------------------------------------------------
# the case
# ---------------------------------------------------------------------------
def run(outdir: Path) -> dict[str, Any]:
    started = time.time()
    outdir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    bodies: dict[str, dict[str, Any]] = {}

    # -- 1. box family: exactly representable, so only round-off is allowed ----
    for name, a, b, c, angles, offset in BOX_CASES:
        mesh = _transformed(mv.box_mesh(a, b, c), angles, offset)
        volume_ref = mv.analytic_box_volume(a, b, c)
        area_ref = mv.analytic_box_area(a, b, c)
        check = mv.check_closed_body(
            mesh,
            name=name,
            family="box",
            reference_volume=volume_ref,
            reference_area=area_ref,
            length=max(a, b, c),
        )
        bodies[name] = check.to_dict()
        rows.append(_row(name, "volume", check.reference_volume, check.measured_volume, mv.EXACT_TOLERANCE))
        rows.append(_row(name, "area", check.reference_area, check.measured_area, mv.EXACT_TOLERANCE))
        rows.append(
            _row(
                name,
                "route_agreement",
                check.measured_volume,
                check.tet_volume,
                mv.ROUTE_TOLERANCE,
                note="divergence theorem vs centroid-tet summation",
            )
        )
        rows.append(
            {
                "case": name,
                "quantity": "topology_closed",
                "reference": 1.0,
                "measured": 1.0 if check.topologically_closed else 0.0,
                "error": 0.0 if check.topologically_closed else 1.0,
                "error_kind": "flag",
                "tolerance": 0.0,
                "passed": bool(check.topologically_closed),
                "note": f"open={check.open_edges} euler={check.euler_characteristic} "
                f"defects={list(check.audit_defects)}",
            }
        )
        third = _third_party_volume(mesh)
        if third is not None:
            rows.append(
                _row(
                    name,
                    "trimesh_volume",
                    check.reference_volume,
                    third,
                    1e-9,
                    note="third-party implementation (trimesh), same triangles",
                )
            )

    # -- 2. sphere ladder: tessellation error, order measured ----------------
    ladder = mv.convergence_study(SPHERE_RADIUS, SPHERE_LEVELS)
    for check in ladder:
        bodies[check.name] = check.to_dict()
    icosa_ref = mv.analytic_icosahedron_volume(SPHERE_RADIUS)
    rows.append(
        _row(
            "sphere_r1_sub0",
            "icosahedron_closed_form",
            icosa_ref,
            ladder[0].measured_volume,
            mv.EXACT_TOLERANCE,
            note="(5/12)(3+sqrt5)a^3, a = 4R/sqrt(10+2sqrt5): the sub0 mesh's own "
            "closed form, independent of any mesh summation",
        )
    )
    order = mv.observed_order(ladder)
    order_area = mv.observed_order(ladder, metric="area_error")
    gate_row = ladder[[row.name for row in ladder].index(f"sphere_r{SPHERE_RADIUS:g}_sub{SPHERE_GATE_LEVEL}")]
    rows.append(
        _row(
            gate_row.name,
            "volume",
            gate_row.reference_volume,
            gate_row.measured_volume,
            mv.RELATIVE_TOLERANCE,
            note="finest level of the ladder (declared before the run)",
        )
    )
    rows.append(
        _row(
            gate_row.name,
            "area",
            gate_row.reference_area,
            gate_row.measured_area,
            mv.RELATIVE_TOLERANCE,
        )
    )
    rows.append(
        _row(
            gate_row.name,
            "route_agreement",
            gate_row.measured_volume,
            gate_row.tet_volume,
            mv.ROUTE_TOLERANCE,
        )
    )
    rows.append(
        {
            "case": gate_row.name,
            "quantity": "topology_closed",
            "reference": 1.0,
            "measured": 1.0 if gate_row.topologically_closed else 0.0,
            "error": 0.0 if gate_row.topologically_closed else 1.0,
            "error_kind": "flag",
            "tolerance": 0.0,
            "passed": bool(gate_row.topologically_closed),
            "note": f"open={gate_row.open_edges} euler={gate_row.euler_characteristic}",
        }
    )
    cylinder_route = mv.analytic_sphere_volume_via_cylinder(SPHERE_RADIUS)
    rows.append(
        _row(
            "reference_provenance",
            "sphere_volume_two_derivations",
            gate_row.reference_volume,
            cylinder_route,
            1e-15,
            note="4/3 pi r^3 vs (2/3)(pi r^2)(2r) (Euclid XII.18, circumscribed cylinder)",
        )
    )
    coarse = ladder[[row.name for row in ladder].index(f"sphere_r{SPHERE_RADIUS:g}_sub{SPHERE_COARSE_LEVEL}")]
    rows.append(
        {
            "case": coarse.name,
            "quantity": "coarse_must_exceed_tolerance",
            "reference": mv.RELATIVE_TOLERANCE,
            "measured": abs(coarse.volume_error),
            "error": 0.0,
            "error_kind": "flag",
            "tolerance": 0.0,
            "passed": abs(coarse.volume_error) > mv.RELATIVE_TOLERANCE,
            "note": "negative control on the ladder: an inscribed polyhedron at this "
            "resolution may not reproduce the sphere",
        }
    )
    rows.append(
        {
            "case": "sphere_ladder",
            "quantity": "observed_order_volume",
            "reference": 2.0,
            "measured": order,
            "error": order - 2.0,
            "error_kind": "flag",
            "tolerance": 0.0,
            "passed": ORDER_BAND[0] <= order <= ORDER_BAND[1],
            "note": f"least-squares slope of log|err| vs log h over sub0..sub{SPHERE_LEVELS[-1]}; "
            f"accepted band [{ORDER_BAND[0]:g}, {ORDER_BAND[1]:g}] around the theoretical 2",
        }
    )

    # -- 3. dimensionless check: the error does not depend on r --------------
    scaled = mv.convergence_study(SPHERE_SCALE_RADIUS, (SPHERE_SCALE_LEVEL,))[0]
    bodies[scaled.name] = scaled.to_dict()
    rows.append(
        _row(
            f"sphere_r{SPHERE_SCALE_RADIUS:g}_sub{SPHERE_SCALE_LEVEL}",
            "volume",
            scaled.reference_volume,
            scaled.measured_volume,
            mv.RELATIVE_TOLERANCE,
        )
    )
    rows.append(
        {
            "case": "scaling_invariance",
            "quantity": "volume_error_difference",
            "reference": ladder[SPHERE_SCALE_LEVEL].volume_error,
            "measured": scaled.volume_error,
            "error": scaled.volume_error - ladder[SPHERE_SCALE_LEVEL].volume_error,
            "error_kind": "absolute",
            "tolerance": 1e-12,
            "passed": abs(scaled.volume_error - ladder[SPHERE_SCALE_LEVEL].volume_error) <= 1e-12,
            "note": f"r={SPHERE_RADIUS:g} vs r={SPHERE_SCALE_RADIUS:g} at sub{SPHERE_SCALE_LEVEL}",
        }
    )

    # -- 4. file round trip: STL (float32) of the box -------------------------
    stl_path = outdir / "box_stl.stl"
    write_stl(mv.box_mesh(*STL_BOX), stl_path)
    stl_mesh = _read_stl_welded(stl_path)
    stl_check = mv.check_closed_body(
        stl_mesh,
        name="box_stl_roundtrip",
        family="box",
        reference_volume=mv.analytic_box_volume(*STL_BOX),
        reference_area=mv.analytic_box_area(*STL_BOX),
        length=max(STL_BOX),
    )
    bodies[stl_check.name] = stl_check.to_dict()
    rows.append(
        _row(
            stl_check.name,
            "volume",
            stl_check.reference_volume,
            stl_check.measured_volume,
            mv.FLOAT32_TOLERANCE,
            note="binary STL stores float32 coordinates, so the format's 1.2e-7 "
            "precision, not the 1e-12 in-memory box line, is the applicable tolerance",
        )
    )

    # -- 5. negative controls: a check that never fires proves nothing --------
    controls = mv.defect_detection(icosphere(SPHERE_RADIUS, 2), radius=SPHERE_RADIUS)
    for control in controls:
        rows.append(
            {
                "case": f"control:{control['name']}",
                "quantity": "defect_detected",
                "reference": 1.0,
                "measured": 1.0 if control["detected"] else 0.0,
                "error": 0.0 if control["detected"] else 1.0,
                "error_kind": "flag",
                "tolerance": 0.0,
                "passed": bool(control["detected"]),
                "note": f"{control['defect']} -> named {control['defects_named']}",
            }
        )

    # -- 6. blind spot: at the accepted resolution one defect is invisible to volume
    fine = icosphere(SPHERE_RADIUS, SPHERE_GATE_LEVEL)
    fine_surface = mv.as_geometry_mesh(fine)
    clean_volume = fine.volume()
    blind: dict[str, Any] = {"clean_volume": clean_volume, "level": SPHERE_GATE_LEVEL}
    for label, defect_mesh in (
        ("duplicate_single", inject_duplicate_faces(fine_surface, count=1)),
        ("flipped_single", inject_flip_faces(fine_surface, count=1)),
    ):
        audit = audit_watertightness(defect_mesh, length=2.0 * SPHERE_RADIUS)
        dirty = TriMesh(
            np.array(defect_mesh.vertices, dtype=np.float64),
            np.array(defect_mesh.faces, dtype=np.int64),
            label,
        )
        effect = (dirty.volume() - clean_volume) / clean_volume
        blind[label] = {
            "volume_effect": effect,
            "invisible_to_volume_gate": abs(effect) <= mv.RELATIVE_TOLERANCE,
            "detected": not audit.passed,
            "duplicate_faces": int(audit.duplicate_faces),
            "winding_conflicts": int(audit.winding_conflicts),
            "nonmanifold_edges": int(audit.nonmanifold_edges),
            "defects": [defect.kind for defect in audit.defects],
        }
        rows.append(
            {
                "case": f"blind_spot:{label}",
                "quantity": "defect_detected_by_audit",
                "reference": 1.0,
                "measured": 0.0 if audit.passed else 1.0,
                "error": 0.0 if not audit.passed else 1.0,
                "error_kind": "flag",
                "tolerance": 0.0,
                "passed": not audit.passed,
                "note": f"volume effect {effect:+.3e} at sub{SPHERE_GATE_LEVEL}; named "
                f"{[defect.kind for defect in audit.defects]}",
            }
        )

    # -- artifacts -----------------------------------------------------------
    passed = all(row["passed"] for row in rows)
    result = {
        "schema": "tensorcad.mesh_volume_watertight/1",
        "anchor": "CASE-CAD-VOLUME-CLOSED",
        "passed": passed,
        "tolerances": {
            "relative": mv.RELATIVE_TOLERANCE,
            "exact": mv.EXACT_TOLERANCE,
            "route": mv.ROUTE_TOLERANCE,
            "order_band": list(ORDER_BAND),
        },
        "gates": rows,
        "failures": [f"{r['case']}.{r['quantity']}" for r in rows if not r["passed"]],
        "bodies": bodies,
        "convergence": {
            "radius": SPHERE_RADIUS,
            "levels": list(SPHERE_LEVELS),
            "observed_order_volume": order,
            "observed_order_area": order_area,
            "first_level_within_tolerance": mv.first_level_within(ladder),
            "rows": [
                {
                    "name": row.name,
                    "faces": row.faces,
                    "h": row.mean_edge_length,
                    "volume": row.measured_volume,
                    "volume_error": row.volume_error,
                    "area_error": row.area_error,
                    "open_edges": row.open_edges,
                    "euler": row.euler_characteristic,
                }
                for row in ladder
            ],
        },
        "negative_controls": controls,
        "blind_spot": blind,
        "seconds": round(time.time() - started, 3),
        "environment": _environment(),
    }
    (outdir / "result.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    _write_convergence_csv(outdir / "convergence.csv", ladder)
    _plot(outdir, ladder, result["gates"], order)
    return result


def _read_stl_welded(path: Path) -> TriMesh:
    """Read our own binary STL back through an OSS reader and weld the triangle soup.

    Welding is what any consumer must do to recover topology from STL; doing it
    explicitly keeps "the file is a soup" (true) apart from "the file is open"
    (only true before welding).
    """
    trimesh = importlib.import_module("trimesh")
    loaded = trimesh.load(str(path), process=False)
    soup = TriMesh(
        np.asarray(loaded.vertices, dtype=np.float64),
        np.asarray(loaded.faces, dtype=np.int64),
        "stl",
    )
    welded, _merged = weld(soup)
    return welded


def _write_convergence_csv(path: Path, ladder: Sequence[Any]) -> Path:
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            ["level_name", "faces", "mean_edge_length", "volume", "volume_error",
             "area_error", "open_edges", "euler"]
        )
        for row in ladder:
            writer.writerow(
                [row.name, row.faces, f"{row.mean_edge_length:.12e}", f"{row.measured_volume:.15e}",
                 f"{row.volume_error:+.6e}", f"{row.area_error:+.6e}",
                 row.open_edges, row.euler_characteristic]
            )
    return path


def _plot(outdir: Path, ladder: Sequence[Any], gates: Sequence[dict[str, Any]], order: float) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    h = np.array([row.mean_edge_length for row in ladder])
    ev = np.array([abs(row.volume_error) for row in ladder])
    ea = np.array([abs(row.area_error) for row in ladder])

    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.loglog(h, ev, "o-", label="volume |error|")
    ax.loglog(h, ea, "s-", label="area |error|")
    ref = ev[0] * (h / h[0]) ** 2
    ax.loglog(h, ref, "k--", lw=1, label="slope 2 (reference)")
    ax.axhline(mv.RELATIVE_TOLERANCE, color="tab:red", ls=":", label="0.1 % acceptance line")
    ax.set_xlabel("mean edge length h (m)")
    ax.set_ylabel("relative error vs closed form")
    ax.set_title(f"inscribed icosphere r={SPHERE_RADIUS:g}: observed order {order:.2f}")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "convergence.png", dpi=140)
    plt.close(fig)

    labels = [f"{g['case']}.{g['quantity']}" for g in gates if g["error_kind"] != "flag"]
    values = [max(abs(g["error"]), 5e-17) for g in gates if g["error_kind"] != "flag"]
    colors = ["tab:green" if g["passed"] else "tab:red" for g in gates if g["error_kind"] != "flag"]
    fig, ax = plt.subplots(figsize=(7.2, 0.32 * len(labels) + 1.6))
    ax.barh(range(len(labels)), values, color=colors)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xscale("log")
    ax.axvline(mv.RELATIVE_TOLERANCE, color="tab:red", ls=":", label="0.1 %")
    ax.set_xlabel("|relative error| (floor 5e-17 for exact rows)")
    ax.set_title("closed-body gates: green pass, red fail")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "gates.png", dpi=140)
    plt.close(fig)


def _environment() -> dict[str, Any]:
    modules = {}
    for name in ("numpy", "trimesh", "matplotlib"):
        try:
            modules[name] = importlib.import_module(name).__version__
        except Exception:
            modules[name] = None
    return {
        "python": platform.python_version(),
        "tensorcad": getattr(tensorcad, "__version__", "unknown"),
        "modules": modules,
        "platform": platform.platform(),
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", default=str(Path(__file__).resolve().parent), help="artifact directory")
    args = parser.parse_args(argv)
    outdir = Path(args.out)
    result = run(outdir)
    print(f"{'case':32s} {'quantity':20s} {'reference':>18s} {'measured':>18s} {'error':>12s}  ok")
    for row in result["gates"]:
        print(
            f"{row['case']:32s} {row['quantity']:20s} {row['reference']:18.10g} "
            f"{row['measured']:18.10g} {row['error']:+12.3e}  {'✓' if row['passed'] else '✗'}"
        )
    conv = result["convergence"]
    print(
        f"\nobserved order (volume) {conv['observed_order_volume']:.3f}, (area) "
        f"{conv['observed_order_area']:.3f}; first level inside 0.1 %: "
        f"{conv['first_level_within_tolerance']}"
    )
    print(f"artifacts: {outdir}  [{result['seconds']} s]  passed={result['passed']}")
    if not result["passed"]:
        print("FAILURES: " + ", ".join(result["failures"]))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
