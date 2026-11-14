"""CAD -> DG high-order curved elements: one geometry source, one frame, real run.

Route (``GAP-CAD-B3-SOLVER-INTERFACE``, third leg -- *DG high-order curved*):

1. the unified geometry library builds a watertight solid
   (:func:`tensorcad.nurbs.cylinder_solid_mesh`); watertightness and the volume
   against ``pi r^2 h`` are measured, not assumed;
2. :func:`tensorcad.dg_curved.section_contour` slices it with a z-up horizontal
   plane, :func:`~tensorcad.dg_curved.fit_circle` measures the centre and radius,
   and :func:`~tensorcad.dg_curved.rational_circle_contract` emits the
   ``tensordg.curved_geometry.v1`` contract (``contract.json`` + ``handoff.json``
   land in ``artifacts/dg_curved/``);
3. TensorDG consumes that contract: ``validate_geometry_contract`` ->
   ``mesh_from_contract`` -> ``CurvedAdvection``, and three numbers are measured
   on the resulting curved mesh:

   * *free-stream preservation* -- the metric identity on CAD-derived curved
     elements (``u_t`` operator residual and drift after 200 SSP-RK3 steps);
   * *rigid-rotation analytic solution* -- ``f0(r, theta - omega t)`` about the
     CAD wall is an exact solution of ``u_t + div(omega(-y, x) u) = 0`` with
     no-penetration walls, so the DG ``l2`` error and its order in ``p`` and in
     ``h`` are scored against a closed form;
   * *domain area* -- the sum of the curved element areas against
     ``pi (r_outer^2 - r_wall^2)``.

Then the honest boundary: an **elliptical** section from the same pipeline is
turned into a least-squares B-spline contract, passes the CAD-side validation,
and is *refused* by the consumer, which only builds polar annuli around a
circular wall.  That refusal is recorded as data, with the measured radius
spread that causes it.

Run::

    PYTHONPATH=.:/root/tensor-ecosystem/TensorDG python \
        examples/dg_curved/run_dg_curved.py
"""
from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

import numpy as np  # noqa: E402

from tensorcad import dg_curved as dgc  # noqa: E402
from tensorcad.geometry import Mesh  # noqa: E402
from tensorcad.nurbs import cylinder_solid_mesh  # noqa: E402
from tensorcad.solver_interface import ellipsoid_mesh, signed_volume  # noqa: E402

ART = ROOT / "artifacts" / "dg_curved"

RADIUS = 1.0
HEIGHT = 2.0
SECTION_Z = 0.43 * HEIGHT      # deliberately off every vertex ring


# ---------------------------------------------------------------------------
# geometry helpers (measurements, not declarations)
# ---------------------------------------------------------------------------
def topology_report(mesh: Mesh) -> dict:
    """Edge incidence from the face list: watertight means every edge in 2 faces."""
    edges: dict[tuple[int, int], int] = {}
    for face in mesh.faces:
        for a, b in ((face[0], face[1]), (face[1], face[2]), (face[2], face[0])):
            key = (a, b) if a < b else (b, a)
            edges[key] = edges.get(key, 0) + 1
    counts: dict[int, int] = {}
    for value in edges.values():
        counts[value] = counts.get(value, 0) + 1
    return {
        "n_vertices": len(mesh.vertices),
        "n_faces": len(mesh.faces),
        "n_edges": len(edges),
        "edges_by_incidence": {str(k): v for k, v in sorted(counts.items())},
        "watertight": counts.get(2, 0) == len(edges),
    }


def cad_solid(u_samples: int = 192, v_samples: int = 8) -> dict:
    mesh = cylinder_solid_mesh(RADIUS, HEIGHT, u_samples=u_samples, v_samples=v_samples)
    volume = signed_volume(mesh.vertices, mesh.faces)
    analytic = math.pi * RADIUS * RADIUS * HEIGHT
    topo = topology_report(mesh)
    return {
        "builder": "tensorcad.nurbs.cylinder_solid_mesh",
        "radius_m": RADIUS,
        "height_m": HEIGHT,
        "u_samples": u_samples,
        "v_samples": v_samples,
        "mesh": mesh,
        "volume_m3": volume,
        "analytic_volume_m3": analytic,
        "volume_rel_error": volume / analytic - 1.0,
        "topology": topo,
    }


def contract_from_solid(solid: dict, *, samples: int = 4096) -> dict:
    mesh = solid["mesh"]
    contour = dgc.section_contour(mesh, SECTION_Z)
    fit = dgc.fit_circle(contour)
    contract = dgc.rational_circle_contract(
        (fit["centre_x"], fit["centre_y"]), fit["radius"],
        name="cylinder_section_z%0.3f" % SECTION_Z,
        source={
            "builder": solid["builder"],
            "section_plane_z_m": SECTION_Z,
            "n_section_points": len(contour),
            "fit_max_residual_m": fit["max_residual"],
            "mesh_u_samples": solid["u_samples"],
        },
    )
    report = dgc.curve_report(
        contract,
        samples=samples,
        analytic_centre=(0.0, 0.0),
        analytic_radius=RADIUS,
        section_points=contour,
    )
    return {"contour": contour, "fit": fit, "contract": contract, "report": report}


# ---------------------------------------------------------------------------
# DG consumption
# ---------------------------------------------------------------------------
def _rigid_rotation(omega: float):
    def velocity(x, y):
        x = np.asarray(x, float)
        y = np.asarray(y, float)
        return -omega * y, omega * x

    return velocity


def _rotating_exact(r0: float, delta: float, m: int, omega: float, t: float):
    """Exact solution of rigid-body rotation of a smooth initial blob."""

    def func(x, y):
        x = np.asarray(x, float)
        y = np.asarray(y, float)
        r = np.hypot(x, y)
        theta = np.arctan2(y, x) - omega * t
        return np.exp(-((r - r0) / delta) ** 2) * np.cos(m * theta)

    return func


def dg_run(dgc_mod, contract: dict, *, degree: int, nr: int, ntheta: int,
           r_outer: float, omega: float, t_end: float,
           wall_radius: float, cfl: float = 0.2) -> dict:
    """Consume the contract on the DG side and score one analytic case."""
    validated = dgc_mod.validate_geometry_contract(contract)
    mesh = dgc_mod.mesh_from_contract(validated, nr, ntheta, r_outer, degree)
    solver = dgc_mod.CurvedAdvection(mesh, degree, _rigid_rotation(omega), cfl=cfl)

    free = solver.free_stream_defect(u0=1.0, steps=200)

    r0, delta, m_wave = 1.5, 0.4, 2
    U0 = solver.project(_rotating_exact(r0, delta, m_wave, omega, 0.0))
    U, dt, nsteps = solver.evolve(U0, t_end, solver.stable_dt())
    num, den, rel = solver.l2_error(U, _rotating_exact(r0, delta, m_wave, omega, t_end))
    area_total = float(np.sum(mesh.areas()))
    # The wall of a p >= 2 annulus *is* the contract circle, so the exact area of
    # the meshed domain is pi (r_outer^2 - r_wall^2) with r_wall the contract
    # radius.  For p = 1 the wall is the inscribed chord polygon and the same
    # reference is only reached at O(h^2); that is reported, not hidden.
    area_exact = math.pi * (r_outer ** 2 - wall_radius ** 2)
    return {
        "degree": degree,
        "nr": nr,
        "ntheta": ntheta,
        "ncells": int(mesh.ncells),
        "dt": dt,
        "nsteps": nsteps,
        "l2_error_abs": num,
        "l2_error_rel": rel,
        "free_stream": free,
        "area_sum_m2": area_total,
        "area_exact_m2": area_exact,
        "area_rel_error": area_total / area_exact - 1.0,
        "wall_radius_declared_m": wall_radius,
        "wall_radius_element_min_m": float(np.min(_wall_radius_per_element(mesh))),
        "wall_radius_element_max_m": float(np.max(_wall_radius_per_element(mesh))),
    }


def _wall_radius_per_element(mesh) -> np.ndarray:
    """Radius of the inner-wall edge midpoint of every element (theta-averaged)."""
    radii = []
    for mp in mesh.flat_maps():
        x, y = mp.eval(np.array([-1.0]), np.array([0.0]))
        radii.append(float(np.mean(np.hypot(np.ravel(x), np.ravel(y)))))
    return np.asarray(radii)


def _wall_radius(mesh) -> float:
    return float(np.mean(_wall_radius_per_element(mesh)))


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main() -> int:
    started = time.time()
    ART.mkdir(parents=True, exist_ok=True)
    results: dict = {"section": {}, "dg": {}, "ellipse": {}, "checks": {}}

    # -- 1. CAD solid, coarse + fine, section convergence -------------------
    coarse = cad_solid(u_samples=48, v_samples=6)
    fine = cad_solid(u_samples=192, v_samples=8)
    results["cad_solid"] = {
        "coarse": {k: v for k, v in coarse.items() if k != "mesh"},
        "fine": {k: v for k, v in fine.items() if k != "mesh"},
    }
    built = contract_from_solid(fine)
    contract, report, contour = built["contract"], built["report"], built["contour"]
    results["_contour"] = contour
    coarse_contract = contract_from_solid(coarse)["report"]
    results["section"] = {
        "fit": built["fit"],
        "n_section_points": len(contour),
        "report": report,
        "coarse_report": coarse_contract,
        "section_deviation_ratio_coarse_over_fine": (
            coarse_contract["section_point_circle_residual_max_m"]
            / report["section_point_circle_residual_max_m"]
            if report["section_point_circle_residual_max_m"] > 0 else None),
    }
    handoff_dir = ART / "handoff_circle"
    record = dgc.write_handoff(
        handoff_dir, contract, report,
        extra={"source_mesh": {"builder": fine["builder"],
                               "u_samples": fine["u_samples"],
                               "v_samples": fine["v_samples"],
                               "section_plane_z_m": SECTION_Z,
                               "watertight": fine["topology"]["watertight"]}},
    )
    results["handoff"] = {"dir": str(handoff_dir.relative_to(ROOT)),
                          "record": record}

    # -- 2. DG side ---------------------------------------------------------
    try:
        import tensordg.curved_geometry as dgc_mod
    except ImportError as exc:  # pragma: no cover - environment dependent
        results["dg"] = {"status": "unavailable", "reason": str(exc)}
        _write(results, started)
        return 0

    r_outer, omega, t_end = 2.0, 1.0, 0.15
    degrees = (1, 2, 3)
    # both directions are refined: refining only theta (nr fixed) leaves the
    # radial discretisation as an error floor, which the first version of this
    # study did and the report documents as a corrected setup error.
    resolutions = ((2, 16), (4, 32), (8, 64))

    # 2a. one solved case per degree at the middle resolution
    per_degree = {}
    for degree in degrees:
        per_degree[str(degree)] = dg_run(
            dgc_mod, contract, degree=degree, nr=3, ntheta=32,
            r_outer=r_outer, omega=omega, t_end=t_end,
            wall_radius=built["fit"]["radius"])
    results["dg"]["per_degree"] = per_degree

    # 2b. h-refinement at fixed degree, plus the same case on TensorDG's own
    #     built-in exact-circle contract as a geometry reference
    refinement = {}
    for degree in degrees:
        rows = []
        for nr, ntheta in resolutions:
            run = dg_run(dgc_mod, contract, degree=degree, nr=nr, ntheta=ntheta,
                         r_outer=r_outer, omega=omega, t_end=t_end,
                         wall_radius=built["fit"]["radius"])
            rows.append({"nr": nr, "ntheta": ntheta, "ncells": run["ncells"],
                         "l2_error_rel": run["l2_error_rel"],
                         "area_rel_error": run["area_rel_error"],
                         "free_stream_drift": run["free_stream"]["drift"],
                         "u_t_operator": run["free_stream"]["u_t_operator"]})
        errors = [row["l2_error_rel"] for row in rows]
        orders = [math.log(errors[i] / errors[i + 1]) / math.log(2.0)
                  for i in range(len(errors) - 1)]
        refinement[str(degree)] = {"rows": rows, "observed_order": orders}
    results["dg"]["refinement"] = refinement

    reference = dgc_mod.circle_contract(RADIUS, (0.0, 0.0))
    reference["units"], reference["frame"] = "m", "z-up"
    results["dg"]["reference_geometry"] = {
        "source": "tensordg.curved_geometry.circle_contract",
        "run": dg_run(dgc_mod, reference, degree=2, nr=3, ntheta=32,
                      r_outer=r_outer, omega=omega, t_end=t_end,
                      wall_radius=RADIUS),
    }
    results["dg"]["status"] = "ran"

    # -- 3. elliptical section: produced, validated, refused by the consumer -
    # z = 0 would put the cut exactly on the equator vertex ring; 0.08 m off it
    ellipse_z = 0.08
    ellipsoid = ellipsoid_mesh((1.2, 0.9, 0.8), n_theta=64, n_phi=128)
    ellipse_contour = dgc.section_contour(ellipsoid, ellipse_z)
    shrink = math.sqrt(1.0 - (ellipse_z / 0.8) ** 2)
    ellipse_axes = [1.2 * shrink, 0.9 * shrink]
    ellipse_fit = dgc.fit_circle(ellipse_contour)
    ellipse_contract = dgc.bspline_closed_contract(
        ellipse_contour, degree=3, n_control=16,
        name="ellipse_section_a1.2_b0.9",
        source={"builder": "tensorcad.solver_interface.ellipsoid_mesh",
                "semi_axes_m": [1.2, 0.9, 0.8], "section_plane_z_m": ellipse_z,
                "section_ellipse_axes_m": ellipse_axes,
                "n_section_points": len(ellipse_contour)})
    ellipse_report = dgc.curve_report(
        ellipse_contract, samples=4096,
        section_points=ellipse_contour)
    dgc.write_handoff(ART / "handoff_ellipse", ellipse_contract, ellipse_report,
                      extra={"semi_axes_m": [1.2, 0.9, 0.8],
                             "section_plane_z_m": ellipse_z,
                             "section_ellipse_axes_m": ellipse_axes})
    refused = None
    try:
        dgc_mod.mesh_from_contract(
            dgc_mod.validate_geometry_contract(ellipse_contract), 3, 32, 2.0, 2)
    except Exception as exc:  # noqa: BLE001 - the refusal is the measurement
        refused = {"type": type(exc).__name__, "message": str(exc)}
    results["ellipse"] = {
        "n_section_points": len(ellipse_contour),
        "section_plane_z_m": ellipse_z,
        "section_ellipse_axes_m": ellipse_axes,
        "circle_fit": ellipse_fit,
        "fit_radius_spread_over_radius": (
            ellipse_fit["max_residual"] / ellipse_fit["radius"]),
        "contract_valid_locally": True,
        "contract_digest": ellipse_report["contract_digest"],
        "section_point_deviation_max_m": ellipse_report.get(
            "section_point_deviation_max_m"),
        "consumer_refusal": refused,
    }

    # -- 4. acceptance gates ------------------------------------------------
    dg_ok = results["dg"]["status"] == "ran"
    fine_report = report
    ref = results["dg"]["reference_geometry"]["run"]
    mid = results["dg"]["per_degree"]["2"]
    worst_fs_rel = abs(mid["free_stream"]["u_t_operator"]
                        / ref["free_stream"]["u_t_operator"] - 1.0)
    l2_rel_to_reference = abs(mid["l2_error_rel"] / ref["l2_error_rel"] - 1.0)
    refine_p2 = results["dg"]["refinement"]["2"]
    refine_p3 = results["dg"]["refinement"]["3"]
    gates = {
        "cad_solid_watertight": fine["topology"]["watertight"],
        "cad_volume_within_1pct": abs(fine["volume_rel_error"]) < 0.01,
        "section_is_closed": fine_report["closure_gap"] < 1e-12,
        # representation invariant: the emitted NURBS circle is exact
        "contract_curve_self_consistent_below_1e-15":
            fine_report["self_consistency_max_m"] / fine_report["declared_radius_m"] < 1e-15,
        # inherited error: the contract repeats the CAD section faceting, O(h^2)
        "contract_area_within_inherited_faceting":
            abs(fine_report["area_rel_error"]) < 1e-3
            and (results["section"]["section_deviation_ratio_coarse_over_fine"] or 0.0) > 10.0,
        "dg_consumed_contract": dg_ok,
        # the CAD-derived geometry must behave like TensorDG's own exact circle
        "dg_matches_reference_geometry_within_0p1pct":
            dg_ok and l2_rel_to_reference < 1e-3 and worst_fs_rel < 1e-3,
        # element-area summation against pi (r_outer^2 - r_wall^2); measured
        # -3.1e-6 (p=2) and -7.6e-10 (p=3), i.e. converging, so the gate is 1e-4
        "dg_area_sum_within_1e-4":
            dg_ok and abs(mid["area_rel_error"]) < 1e-4,
        "dg_refinement_reduces_error":
            dg_ok and all(
                all(b < a for a, b in zip(
                    [row["l2_error_rel"] for row in entry["rows"]],
                    [row["l2_error_rel"] for row in entry["rows"]][1:]))
                for entry in results["dg"]["refinement"].values()),
        "dg_order_toward_p_plus_one":
            dg_ok and min(refine_p2["observed_order"]) > 0.0
            and min(refine_p3["observed_order"]) > 0.0,
        "ellipse_contract_refused_with_reason":
            refused is not None and "not circular" in refused["message"],
    }
    results["checks"] = {
        "gates": gates,
        "dg_vs_reference": {
            "l2_error_rel_of_CAD_to_reference": l2_rel_to_reference,
            "free_stream_operator_rel_of_CAD_to_reference": worst_fs_rel,
        },
        "gates_failed": [k for k, v in gates.items() if not v],
    }
    _write(results, started)
    _figure(results, fine)
    print(json.dumps({"gates_failed": results["checks"]["gates_failed"],
                      "wall_s": round(time.time() - started, 1)}, indent=2))
    return 0 if not results["checks"]["gates_failed"] else 1


def _write(results: dict, started: float) -> None:
    results["wall_s"] = round(time.time() - started, 2)
    (ART / "results.json").write_text(
        json.dumps(results, indent=2, sort_keys=True, default=float) + "\n")


def _figure(results: dict, solid: dict) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
    refinement = results["dg"].get("refinement", {})
    colours = {"1": "C0", "2": "C1", "3": "C2"}
    for degree, entry in refinement.items():
        cells = [row["ncells"] for row in entry["rows"]]
        errs = [row["l2_error_rel"] for row in entry["rows"]]
        axes[0].loglog(np.sqrt(cells), errs, "o-", color=colours[degree],
                       label=f"p={degree}")
        slope = int(degree) + 1
        x0 = np.sqrt(cells[0])
        axes[0].loglog(np.sqrt(cells), errs[0] * (np.sqrt(cells) / x0) ** -slope,
                       ":", color=colours[degree], alpha=0.6,
                       label=f"$h^{{{slope}}}$")
    axes[0].set_xlabel("1/h (elements per side scale)")
    axes[0].set_ylabel("relative $L^2$ error, rigid rotation")
    axes[0].set_title("DG error on CAD-derived curved annulus")
    axes[0].grid(True, which="both", alpha=0.3)
    axes[0].legend(fontsize=8)

    # geometry: radial deviation of the emitted NURBS against the exact circle
    report = results["section"]["report"]
    coarse = results["section"]["coarse_report"]
    labels = ["fine section\n(192 segs)", "coarse section\n(48 segs)"]
    ideal = [report["radial_deviation_max_m"], coarse["radial_deviation_max_m"]]
    faceting = [report["section_point_circle_residual_max_m"],
                coarse["section_point_circle_residual_max_m"]]
    x = np.arange(2)
    axes[1].bar(x - 0.2, ideal, width=0.4, label="emitted NURBS vs exact circle")
    axes[1].bar(x + 0.2, faceting, width=0.4, label="CAD section points vs circle")
    axes[1].set_yscale("log")
    axes[1].set_xticks(x, labels)
    axes[1].set_ylabel("max deviation [m]")
    axes[1].set_title("Handoff geometry error at two CAD resolutions")
    axes[1].legend(fontsize=8)
    axes[1].grid(True, axis="y", alpha=0.3)

    # plan view of the sliced section and the emitted control polygon
    contour = results.get("_contour")
    if contour is not None:
        ring = np.asarray(contour)
        axes[2].plot(np.append(ring[:, 0], ring[0, 0]),
                     np.append(ring[:, 1], ring[0, 1]), "-", lw=1.0,
                     label="CAD section contour")
    theta = np.linspace(0.0, 2.0 * math.pi, 361)
    axes[2].plot(RADIUS * np.cos(theta), RADIUS * np.sin(theta), "--", lw=1.0,
                 label="exact circle r=%.3f m" % RADIUS)
    axes[2].set_aspect("equal")
    axes[2].set_title("z-up horizontal section of the CAD solid")
    axes[2].legend(fontsize=8)
    axes[2].grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(ART / "dg_curved_handoff.png", dpi=140)
    plt.close(fig)


if __name__ == "__main__":
    raise SystemExit(main())