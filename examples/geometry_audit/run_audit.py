"""Run the TensorCAD hull-geometry audit pipeline and write its evidence.

Everything on the report page and in ``artifacts/m5283_geometry_audit/gates_verified.json``
is produced here, from measured numbers -- no figure or number is transcribed by
hand.  Exits non-zero if any acceptance gate fails, so a stale page cannot pass CI.

Cases
-----
``C1``  watertightness recall on injected seam gaps (and a below-tolerance control)
``C2``  Gaussian-curvature-jump detection on injected bumps, with false-positive control
``C3``  offsets-table (lines plan) versus mesh integration of the same hull, < 1 %
``C4``  negative cases: gap, normal flip, duplicate shell, ``Cb`` offset -- must be rejected
``C5``  Wigley III and Series 60 hydrostatics against their closed forms, < 2 %
``C6``  the pytest suite and this gates file

Usage::

    python examples/geometry_audit/run_audit.py [--out artifacts/m5283_geometry_audit]
"""
from __future__ import annotations

import argparse
import json
import math
import random
import subprocess
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # noqa: E402  (must precede pyplot)
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tensorcad import mesh as meshlib  # noqa: E402
from tensorcad.exchange import (  # noqa: E402
    ExchangeError,
    ExchangeUnavailable,
    PatchGeometry,
    load_exchange,
    round_trip_facets,
    supported_routes,
    write_patchset,
    write_step_facets,
)
from tensorcad.geometry_audit import (  # noqa: E402
    CURVATURE_JUMP_THRESHOLD,
    GAP_TOLERANCE_FRACTION,
    HYDROSTATIC_TOLERANCE,
    HullSpec,
    audit_fairness,
    audit_hull_geometry,
    audit_watertightness,
    gaussian_curvature,
    inject_bump,
    inject_cb_offset,
    inject_duplicate_faces,
    inject_flip_faces,
    inject_seam_offset,
    integrate_offsets,
    offsets_from_function,
    offsets_from_mesh,
    write_gates_json,
)
from tensorcad.series60 import (  # noqa: E402
    Series60Parameters,
    analytic_series60_hydrostatics,
    build_series60_mesh,
)
from tensorcad.wigley import (  # noqa: E402
    WigleyParameters,
    analytic_block_coefficient,
    analytic_volume,
    analytic_wetted_surface_area,
    build_mesh,
    half_breadth,
    reference_parameters,
)

MISSION = "5283"
CASE_SEED = 20261106


# --------------------------------------------------------------------------- #
# cases
# --------------------------------------------------------------------------- #
def wigley_case() -> tuple[WigleyParameters, "object", object, object]:
    p = reference_parameters(length=10.0, beam=1.0, draft=0.625, nx=161, nz=41)
    mesh, _ = build_mesh(p)
    table = offsets_from_function(
        lambda x, z: half_breadth(p, x, p.draft - z),
        length=p.length,
        beam=p.beam,
        draft=p.draft,
        n_stations=201,
        n_levels=51,
        source="Wigley III analytic half-breadths",
    )
    spec = HullSpec(
        length=p.length,
        beam=p.beam,
        draft=p.draft,
        block_coefficient=analytic_block_coefficient(),
        name="Wigley III",
    )
    return p, mesh, table, spec


def series60_case(cb: float = 0.70) -> tuple[Series60Parameters, object, object, object]:
    s = Series60Parameters(length=10.0, beam=1.5, draft=0.6, nx=161, nz=41, cb=cb)
    mesh, _ = build_series60_mesh(s)
    n, m = s.exponents
    table = offsets_from_function(
        lambda x, z: 0.5 * s.beam * max(0.0, 1.0 - abs(2.0 * x / s.length - 1.0) ** n)
        * (max(0.0, min(z / s.draft, 1.0)) ** (1.0 / m)),
        length=s.length,
        beam=s.beam,
        draft=s.draft,
        n_stations=201,
        n_levels=51,
        source=f"Series 60 Cb={cb:.2f} parent-form half-breadths",
    )
    spec = HullSpec(
        length=s.length,
        beam=s.beam,
        draft=s.draft,
        block_coefficient=cb,
        name=f"Series 60 Cb={cb:.2f}",
    )
    return s, mesh, table, spec


# --------------------------------------------------------------------------- #
# C1: watertightness recall
# --------------------------------------------------------------------------- #
def c1_gap_recall(mesh, length: float, out_dir: Path) -> dict:
    rng = random.Random(CASE_SEED)
    rows = []
    for index in range(24):
        fraction = 0.001 + 0.019 * rng.random()  # 0.1 % .. 2 % of L
        delta = fraction * length
        x_split = length * (0.3 + 0.4 * rng.random())
        defected, patch_ids = inject_seam_offset(mesh, delta=delta, x_split=x_split)
        audit = audit_watertightness(defected, length=length)
        gap_defects = [d for d in audit.defects if d.kind in ("seam_gap", "open_boundary")]
        rows.append(
            {
                "case": index,
                "injected_gap_m": delta,
                "injected_gap_relative": fraction,
                "measured_gap_m": audit.max_gap,
                "measured_gap_relative": audit.max_gap / length,
                "detected": bool(gap_defects),
                "kinds": sorted({d.kind for d in audit.defects}),
                "patch_seams": audit.patch_seams,
                "seam_patches": len(set(patch_ids)),
            }
        )
    detected = sum(1 for row in rows if row["detected"])
    recall = detected / len(rows)

    # Below-tolerance control: a 0.05 % of L gap must not raise a gap defect, while
    # the mesh is still rejected as two disconnected shells -- both statements hold.
    control_mesh, _ = inject_seam_offset(mesh, delta=0.0005 * length, x_split=0.5 * length)
    control = audit_watertightness(control_mesh, length=length)
    control_gap_defects = [d for d in control.defects if d.kind in ("seam_gap", "open_boundary")]

    # Clean-control: the undamaged hulls must raise nothing.
    clean = audit_watertightness(mesh, length=length)
    (out_dir / "c1_gap_recall.json").write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")
    return {
        "id": "C1",
        "criterion": "seam gap > 0.001 L is detected (recall >= 0.95)",
        "injected_cases": len(rows),
        "detected": detected,
        "recall": recall,
        "recall_threshold": 0.95,
        "tolerance_relative": GAP_TOLERANCE_FRACTION,
        "measured_to_injected_ratio_median": float(
            np.median([row["measured_gap_m"] / row["injected_gap_m"] for row in rows])
        ),
        "below_tolerance_control": {
            "injected_gap_relative": 0.0005,
            "measured_gap_relative": control.max_gap / length,
            "gap_defects": [d.kind for d in control_gap_defects],
            "verdict": control.verdict if hasattr(control, "verdict") else (
                "accepted" if control.passed else "rejected"
            ),
            "note": (
                "the gap is below the 0.1 % tolerance so no gap defect is raised, but the two "
                "shells are disconnected and the topology gate rejects the model anyway"
            ),
        },
        "clean_control": {
            "defects": [d.kind for d in clean.defects],
            "passed": clean.passed,
        },
        "passed": recall >= 0.95 and clean.passed and not control_gap_defects,
        "evidence": "REAL (synthetic seams with exact ground truth on a measured mesh)",
    }


# --------------------------------------------------------------------------- #
# C2: curvature-jump detection
# --------------------------------------------------------------------------- #
def vertex_jump_ratios(field, base_field) -> list[float]:
    floor = 0.1 * base_field.mean_abs_curvature
    ratios = []
    for index, flag in enumerate(field.interior):
        if not flag:
            ratios.append(float("nan"))
            continue
        base = base_field.curvature[index]
        value = field.curvature[index]
        if not (math.isfinite(base) and math.isfinite(value)):
            ratios.append(float("nan"))
            continue
        ratios.append(abs(value - base) / max(abs(base), floor))
    return ratios


def spatial_jump_ratios(mesh, field, centre, radius: float, base_field) -> list[float]:
    """Recompute the spatial |dK|/K jump inside a bump footprint, independently.

    The audit's region list is *not* read here.  The comparison set is rebuilt from
    the mesh -- interior vertices, no designed creases and no vertex within
    ``CREASE_DILATION_RINGS`` of one -- and the criterion is the same relative jump
    the audit applies, so "this bump really does break C1 continuity" is measured
    rather than assumed.  The location ground truth stays the injected centre,
    which comes from the injector and not from the detector.
    """
    from tensorcad.geometry_audit import (
        CREASE_DILATION_RINGS,
        CREASE_DEGREES,
        _dilate,
        _edge_face_table,
        _unit,
        dihedral_angle_degrees,
    )

    floor = max(0.1 * base_field.mean_abs_curvature, 1e-12)
    normals = [_unit(meshlib.face_normal(mesh, face)) for face in mesh.faces]
    edge_faces = _edge_face_table(mesh)
    crease_vertices: set[int] = set()
    for edge, incident in edge_faces.items():
        if len(incident) != 2:
            continue
        if dihedral_angle_degrees(mesh, edge, normals=normals, incident=incident) > CREASE_DEGREES:
            crease_vertices.update(edge)
    excluded = _dilate(mesh, crease_vertices, CREASE_DILATION_RINGS)
    inside = {
        vertex
        for vertex, point in enumerate(mesh.vertices)
        if math.dist(point, centre) <= radius and vertex not in excluded
    }
    ratios: list[float] = []
    for edge, incident in edge_faces.items():
        if len(incident) != 2:
            continue
        u, v = edge
        if u not in inside or v not in inside:
            continue
        if not (field.interior[u] and field.interior[v]):
            continue
        ku, kv = field.curvature[u], field.curvature[v]
        if not (math.isfinite(ku) and math.isfinite(kv)):
            continue
        ratios.append(abs(ku - kv) / max(0.5 * (abs(ku) + abs(kv)), floor))
    return ratios


def c2_curvature_recall(mesh, out_dir: Path, seed: int = CASE_SEED) -> dict:
    base_field = gaussian_curvature(mesh)
    box = mesh.bounding_box()
    rng = random.Random(seed)
    dx = box.sizes[0] / 160.0
    sigma = 3.5 * dx
    rows = []
    for index in range(12):
        x = box.minimum[0] + box.sizes[0] * (0.2 + 0.6 * rng.random())
        z = box.minimum[2] + box.sizes[2] * (0.35 + 0.5 * rng.random())
        # put the bump centre on the surface: nearest vertex of the sampled hull
        target = None
        best = math.inf
        for vertex, point in enumerate(mesh.vertices):
            if point[1] <= 0.0:
                continue
            distance = (point[0] - x) ** 2 + (point[2] - z) ** 2
            if distance < best:
                best = distance
                target = vertex
        centre = mesh.vertices[target]
        # Amplitudes straddle the criterion: some bumps are far above the 50 % jump,
        # some are barely above it, some are below and must not be in the denominator.
        factor = 0.03 + 0.35 * rng.random()
        amplitude = factor * base_field.mean_abs_curvature * sigma * sigma
        bumped = inject_bump(mesh, centre=centre, amplitude=amplitude, sigma=sigma, field=base_field)
        field = gaussian_curvature(bumped)
        change = vertex_jump_ratios(field, base_field)  # |dK|/K at a vertex, pre -> post
        spatial = spatial_jump_ratios(bumped, field, centre, 2.0 * sigma, base_field)
        precondition = bool(spatial) and max(spatial) > CURVATURE_JUMP_THRESHOLD
        audit = audit_fairness(bumped, field=field)
        marked = any(
            region.count_as_defect and math.dist(region.centroid, centre) <= 3.0 * sigma
            for region in audit.regions
        )
        rows.append(
            {
                "case": index,
                "centre": [round(value, 6) for value in centre],
                "amplitude_m": amplitude,
                "amplitude_over_meanK_sigma2": factor,
                "sigma_m": sigma,
                "max_curvature_change_ratio": round(
                    max((value for value in change if math.isfinite(value)), default=float("nan")), 3
                ),
                "max_spatial_jump_ratio": round(max(spatial), 3) if spatial else float("nan"),
                "has_above_threshold_region": precondition,
                "regions_flagged": len([r for r in audit.regions if r.count_as_defect]),
                "flagged_area_fraction": audit.flagged_area_fraction,
                "marked": marked,
            }
        )
    evaluated = [row for row in rows if row["has_above_threshold_region"]]
    hits = sum(1 for row in evaluated if row["marked"])
    recall = hits / len(evaluated) if evaluated else 0.0

    # False-positive control on the undamaged hull, at its production resolution.
    clean = audit_fairness(mesh, field=base_field)
    (out_dir / "c2_curvature_recall.json").write_text(
        json.dumps(rows, indent=2) + "\n", encoding="utf-8"
    )
    return {
        "id": "C2",
        "criterion": "regions with |dK|/K > 50% are marked (recall >= 0.95)",
        "injected_cases": len(rows),
        "cases_above_threshold": len(evaluated),
        "cases_below_threshold": len(rows) - len(evaluated),
        "marked": hits,
        "recall": recall,
        "recall_threshold": 0.95,
        "threshold": CURVATURE_JUMP_THRESHOLD,
        "footprint_radius_sigma": 2.0,
        "min_spatial_jump_ratio": min(
            (row["max_spatial_jump_ratio"] for row in evaluated), default=float("nan")
        ),
        "max_spatial_jump_ratio": max(
            (row["max_spatial_jump_ratio"] for row in rows), default=float("nan")
        ),
        "clean_control": {
            "flagged_regions": len([r for r in clean.regions if r.count_as_defect]),
            "flagged_area_fraction": clean.flagged_area_fraction,
            "max_jump_ratio": clean.max_jump_ratio,
            "passed": clean.passed,
        },
        "passed": recall >= 0.95 and clean.passed,
        "evidence": "REAL (analytic base field, Gaussian bumps, per-case |dK|/K measured)",
    }


# --------------------------------------------------------------------------- #
# C3 / C5: hydrostatics
# --------------------------------------------------------------------------- #
def hydrostatic_gate(report, reference: dict[str, float], reference_label: str) -> dict:
    mesh_hydro = report.mesh_hydrostatics
    measured = {
        "displacement_volume": mesh_hydro.displacement_volume,
        "wetted_surface_area": mesh_hydro.wetted_surface_area,
        "block_coefficient": mesh_hydro.block_coefficient,
    }
    if report.consistency is not None:
        consistency = dict(report.consistency.differences)
    else:
        table_hydro = integrate_offsets(report.tables[0])
        consistency = {
            "displacement_volume": (table_hydro.displacement_volume - mesh_hydro.displacement_volume)
            / mesh_hydro.displacement_volume,
            "wetted_surface_area": (table_hydro.wetted_surface_area - mesh_hydro.wetted_surface_area)
            / mesh_hydro.wetted_surface_area,
            "block_coefficient": (table_hydro.block_coefficient - mesh_hydro.block_coefficient)
            / mesh_hydro.block_coefficient,
        }
    reference_errors = {
        key: (measured[key] - reference[key]) / reference[key]
        for key in measured
        if key in reference and reference[key]
    }
    return {
        "reference_label": reference_label,
        "measured": {key: round(value, 9) for key, value in measured.items()},
        "reference": {key: round(value, 9) for key, value in reference.items()},
        "reference_error": {key: round(value, 9) for key, value in reference_errors.items()},
        "table_vs_mesh": {key: round(value, 9) for key, value in consistency.items()},
        "table_source": report.tables[0].source,
        "wetted_area_reference_available": "wetted_surface_area" in reference,
    }


# --------------------------------------------------------------------------- #
# C4: negative cases
# --------------------------------------------------------------------------- #
def c4_negative_cases(mesh, spec, table, length: float, out_dir: Path) -> dict:
    cases = []

    gap_mesh, _ = inject_seam_offset(mesh, delta=0.01 * length, x_split=0.45 * length)
    flip_mesh = inject_flip_faces(mesh, count=40, seed=11)
    duplicate_mesh = inject_duplicate_faces(mesh, count=24, seed=12)
    cb_mesh = inject_cb_offset(mesh, factor=1.05)

    for name, candidate, expectation in (
        ("seam_gap_1pct_L", gap_mesh, "seam_gap"),
        ("normal_flip_40_faces", flip_mesh, "normal_inconsistency"),
        ("duplicate_shell_24_faces", duplicate_mesh, "duplicate_facets"),
        ("block_coefficient_plus_5pct", cb_mesh, "block_coefficient_mismatch"),
    ):
        report = audit_hull_geometry(candidate, spec=spec, offsets=table)
        kinds = sorted({defect.kind for defect in report.defects})
        deviations = {
            defect.kind: round(defect.relative, 9)
            for defect in report.defects
            if defect.relative is not None
        }
        details = [defect.detail for defect in report.defects]
        cases.append(
            {
                "case": name,
                "expected_defect": expectation,
                "verdict": report.verdict,
                "rejected": not report.accepted,
                "stages_run": list(report.stages_run),
                "defect_kinds": kinds,
                "expected_found": expectation in kinds,
                "relative_deviations": deviations,
                "details": details,
                "hydrostatics_computed": report.mesh_hydrostatics is not None,
            }
        )

    rows_ok = all(case["rejected"] and case["expected_found"] for case in cases)
    cb_case = next(case for case in cases if case["case"] == "block_coefficient_plus_5pct")
    (out_dir / "c4_negative_cases.json").write_text(json.dumps(cases, indent=2) + "\n",
                                                    encoding="utf-8")
    return {
        "id": "C4",
        "criterion": "a +5 % Cb offset (and every injected defect) is refused with the deviation named",
        "cases": cases,
        "cb_reported_deviation": cb_case["relative_deviations"].get("block_coefficient_mismatch"),
        "cb_reported_detail": next(
            (detail for detail in cb_case["details"] if detail.startswith("realised Cb")), ""
        ),
        "passed": rows_ok,
        "evidence": "REAL (defects injected into a measured mesh, audit re-run end to end)",
    }


# --------------------------------------------------------------------------- #
# exchange round trip
# --------------------------------------------------------------------------- #
def exchange_gate(mesh, length: float, out_dir: Path) -> dict:
    trip = round_trip_facets(mesh, directory=out_dir)
    defected, _ = inject_seam_offset(mesh, delta=0.01 * length, x_split=0.45 * length)
    step_path = out_dir / "defected_hull.step"
    write_step_facets(step_path, defected)
    reimported = load_exchange(step_path)
    write_patchset(out_dir / "defected_hull.patchset.json", PatchGeometry(
        mesh=defected, patch_ids=tuple("hull" for _ in defected.faces), source="injector"))
    from_step = audit_watertightness(reimported.mesh, length=length)
    from_patchset = audit_watertightness(
        load_exchange(out_dir / "defected_hull.patchset.json").mesh, length=length
    )
    sizes = {}
    for suffix in (".iges", ".igs"):
        probe = out_dir / f"probe{suffix}"
        probe.write_text("IGES probe\n", encoding="utf-8")
        try:
            load_exchange(probe)
            sizes[suffix] = "loaded"
        except ExchangeUnavailable as error:
            sizes[suffix] = f"ExchangeUnavailable: {error}"
        except ExchangeError as error:  # pragma: no cover
            sizes[suffix] = f"ExchangeError: {error}"
        probe.unlink()
    return {
        "round_trip": {key: round(value, 12) for key, value in trip.items()},
        "step_file_bytes": step_path.stat().st_size,
        "defected_step_reimport_detects_gap": not from_step.passed,
        "defected_step_reimport_max_gap_relative": from_step.max_gap / length,
        "patchset_reimport_detects_gap": not from_patchset.passed,
        "step_entity_metadata": reimported.metadata,
        "iges_route": sizes,
        "routes": supported_routes(),
        "passed": trip["vertex_deviation_m"] < 1e-12
        and not from_step.passed
        and not from_patchset.passed,
    }


# --------------------------------------------------------------------------- #
# figures
# --------------------------------------------------------------------------- #
def figure_gap_recall(c1: dict, out_dir: Path) -> Path:
    rows = json.loads((out_dir / "c1_gap_recall.json").read_text(encoding="utf-8"))
    injected = [row["injected_gap_relative"] for row in rows]
    measured = [row["measured_gap_relative"] for row in rows]
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2))
    ax = axes[0]
    ax.scatter(injected, measured, s=42, color="tab:blue", label="injected seam")
    ax.plot([0, 0.02], [0, 0.02], "k--", lw=1.0, label="measured = injected")
    ax.axhline(GAP_TOLERANCE_FRACTION, color="tab:red", ls=":", lw=1.4,
               label="gate tolerance 0.001 L")
    ax.set_xlabel("injected gap / L")
    ax.set_ylabel("measured gap / L")
    ax.set_title(f"C1  gap detection, recall {c1['detected']}/{c1['injected_cases']}")
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(alpha=0.3)
    ax = axes[1]
    ax.hist(
        [row["measured_gap_relative"] / row["injected_gap_relative"] for row in rows],
        bins=8, color="tab:green", alpha=0.8,
    )
    ax.set_xlabel("measured / injected gap")
    ax.set_ylabel("cases")
    ax.set_title("ground-truth agreement (1.0 = exact)")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    path = out_dir / "fig_c1_gap_recall.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def figure_curvature(mesh, out_dir: Path) -> Path:
    base = gaussian_curvature(mesh)
    box = mesh.bounding_box()
    centre_index = min(
        (
            vertex
            for vertex, point in enumerate(mesh.vertices)
            if point[1] > 0.0 and abs(point[0] - 7.0) < 0.2
        ),
        key=lambda vertex: abs(mesh.vertices[vertex][2] + 0.35),
    )
    centre = mesh.vertices[centre_index]
    sigma = 3.5 * box.sizes[0] / 160.0
    bumped = inject_bump(
        mesh, centre=centre, amplitude=12.0 * base.mean_abs_curvature * sigma * sigma,
        sigma=sigma, field=base,
    )
    field = gaussian_curvature(bumped)
    audit = audit_fairness(bumped, field=field)
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.4), sharey=True)
    for ax, values, title in (
        (axes[0], base, f"healthy hull: {len([r for r in audit_fairness(mesh, field=base).regions if r.count_as_defect])} regions"),
        (axes[1], field, f"bump injected: {len([r for r in audit.regions if r.count_as_defect])} regions"),
    ):
        xs = [point[0] for point in bumped.vertices]
        zs = [point[2] for point in bumped.vertices]
        ks = [min(abs(value), 20.0) if math.isfinite(value) else float("nan") for value in values.curvature]
        scatter = ax.scatter(xs, zs, c=ks, s=6, cmap="inferno_r", vmin=0.0, vmax=2.0)
        if ax is axes[1]:
            for region in audit.regions:
                if region.count_as_defect:
                    ax.plot([region.centroid[0]], [region.centroid[2]], marker="o",
                            mfc="none", mec="cyan", ms=16, mew=1.6)
        ax.set_xlabel("x [m]")
        ax.set_title(title, fontsize=9)
        ax.set_aspect("equal", adjustable="box")
    axes[0].set_ylabel("z [m]")
    fig.colorbar(scatter, ax=axes, label="|K| [1/m^2], clipped at 20", shrink=0.85)
    fig.suptitle("C2  discrete Gaussian curvature; cyan circles mark flagged jump regions", fontsize=10)
    path = out_dir / "fig_c2_curvature.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


def figure_hydrostatics(cases: dict, out_dir: Path) -> Path:
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    ax = axes[0]
    labels = []
    values = []
    for name, case in cases.items():
        for key, value in case["table_vs_mesh"].items():
            labels.append(f"{name}\n{key}")
            values.append(100.0 * value)
    positions = np.arange(len(values))
    colours = ["tab:green" if abs(value) <= 100.0 * HYDROSTATIC_TOLERANCE else "tab:red"
               for value in values]
    ax.bar(positions, values, color=colours)
    ax.axhline(100.0 * HYDROSTATIC_TOLERANCE, color="k", ls="--", lw=1.0)
    ax.axhline(-100.0 * HYDROSTATIC_TOLERANCE, color="k", ls="--", lw=1.0)
    ax.set_xticks(positions)
    ax.set_xticklabels(labels, fontsize=6, rotation=90)
    ax.set_ylabel("offsets-table vs mesh [%]")
    ax.set_title("C3  two independent integrations, 1 % band")
    ax.grid(alpha=0.3, axis="y")
    ax = axes[1]
    for name, case in cases.items():
        errors = case["reference_error"]
        keys = list(errors)
        ax.bar(np.arange(len(keys)) + (0.35 if name == "Series 60 Cb=0.70" else -0.05),
               100.0 * np.array([errors[key] for key in keys]),
               width=0.3, label=name)
        ax.set_xticks(np.arange(len(keys)))
        ax.set_xticklabels(keys, fontsize=7, rotation=20)
    ax.axhline(2.0, color="k", ls="--", lw=1.0)
    ax.axhline(-2.0, color="k", ls="--", lw=1.0)
    ax.set_ylabel("mesh vs closed form [%]")
    ax.set_title("C5  standard hulls against their closed forms, 2 % band")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    path = out_dir / "fig_c3_c5_hydrostatics.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default=str(REPO_ROOT / "artifacts" / "m5283_geometry_audit"))
    arguments = parser.parse_args()
    out_dir = Path(arguments.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    gates: list[dict] = []
    cases: dict[str, dict] = {}

    # --- Wigley III ------------------------------------------------------- #
    p, wigley_mesh, wigley_table, wigley_spec = wigley_case()
    wigley_report = audit_hull_geometry(wigley_mesh, spec=wigley_spec, offsets=wigley_table)
    wigley_report.to_json(out_dir / "audit_wigley_iii.json")
    wigley_reference = {
        "displacement_volume": analytic_volume(p.length, p.beam, p.draft),
        "wetted_surface_area": analytic_wetted_surface_area(p.length, p.beam, p.draft),
        "block_coefficient": analytic_block_coefficient(),
    }
    cases["Wigley III"] = hydrostatic_gate(
        wigley_report, wigley_reference, "Wigley III closed forms (integrated analytic surface)"
    )
    mesh_offsets = offsets_from_mesh(wigley_mesh, n_stations=81, n_levels=41)
    mesh_table_hydro = integrate_offsets(mesh_offsets)
    cases["Wigley III"]["mesh_digitised_table"] = {
        "source": mesh_offsets.source,
        "volume_relative_difference": (
            mesh_table_hydro.displacement_volume - wigley_report.mesh_hydrostatics.displacement_volume
        )
        / wigley_report.mesh_hydrostatics.displacement_volume,
        "wetted_relative_difference": (
            mesh_table_hydro.wetted_surface_area - wigley_report.mesh_hydrostatics.wetted_surface_area
        )
        / wigley_report.mesh_hydrostatics.wetted_surface_area,
    }

    # --- Series 60 -------------------------------------------------------- #
    s, s60_mesh, s60_table, s60_spec = series60_case(0.70)
    s60_report = audit_hull_geometry(s60_mesh, spec=s60_spec, offsets=s60_table)
    s60_report.to_json(out_dir / "audit_series60_cb070.json")
    s60_reference = analytic_series60_hydrostatics(s)
    reference = {
        "displacement_volume": s60_reference["displacement_volume"],
        "block_coefficient": s60_reference["block_coefficient"],
    }
    cases["Series 60 Cb=0.70"] = hydrostatic_gate(
        s60_report, reference, f"Series 60 Cb=0.70 parent form, Cp={s.exponents[0]:.4f}"
    )
    cases["Series 60 Cb=0.70"]["wetted_area_reference_available"] = False

    # --- C1 / C2 on the Wigley hull --------------------------------------- #
    c1 = c1_gap_recall(wigley_mesh, p.length, out_dir)
    c2 = c2_curvature_recall(wigley_mesh, out_dir)
    gates.extend([c1, c2])

    # --- C3 ---------------------------------------------------------------#
    worst = max(
        abs(value)
        for case in cases.values()
        for value in case["table_vs_mesh"].values()
    )
    gates.append(
        {
            "id": "C3",
            "criterion": "volume / wetted area / Cb from the offsets table and from the mesh agree within 1 %",
            "tolerance": HYDROSTATIC_TOLERANCE,
            "worst_relative_difference": worst,
            "cases": {name: case["table_vs_mesh"] for name, case in cases.items()},
            "passed": worst <= HYDROSTATIC_TOLERANCE,
            "evidence": "REAL (Simpson on a station x waterline table versus divergence theorem on triangles)",
        }
    )

    # --- C4 ---------------------------------------------------------------#
    gates.append(c4_negative_cases(wigley_mesh, wigley_spec, wigley_table, p.length, out_dir))

    # --- C5 ---------------------------------------------------------------#
    worst_reference = max(
        abs(value) for case in cases.values() for value in case["reference_error"].values()
    )
    gates.append(
        {
            "id": "C5",
            "criterion": "Wigley III and Series 60 hydrostatics within 2 % of their closed forms",
            "tolerance": 0.02,
            "worst_relative_error": worst_reference,
            "cases": {name: case["reference_error"] for name, case in cases.items()},
            "notes": {
                "Wigley III": "Cb = Cp*Cm = 4/9, Cp = Cm = Cwp = 2/3; wetted area from the analytic surface integral",
                "Series 60 Cb=0.70": "parent-form closed forms: Cb = n*m/((n+1)(m+1)); no closed-form wetted area exists, so C3 carries that check",
            },
            "passed": worst_reference <= 0.02,
            "evidence": "REAL (mesh divergence theorem against independent closed forms)",
        }
    )

    # --- exchange round trip ----------------------------------------------#
    exchange = exchange_gate(wigley_mesh, p.length, out_dir)

    # --- C6 ---------------------------------------------------------------#
    pytest_log = out_dir / "pytest.txt"
    # The artifact self-check inside that module is deselected here: it validates the
    # file this run is about to write, so it is run at the end (below) instead.
    process = subprocess.run(
        [
            sys.executable, "-m", "pytest", "-q", "tests/test_geometry_audit.py",
            "--deselect",
            "tests/test_geometry_audit.py::TestGatesArtifact::test_gates_file_is_complete",
        ],
        cwd=str(REPO_ROOT), capture_output=True, text=True,
    )
    pytest_log.write_text(process.stdout + process.stderr, encoding="utf-8")
    c6 = {
        "id": "C6",
        "criterion": "every criterion is a pytest case and gates_verified.json is produced",
        "test_module": "tests/test_geometry_audit.py",
        "pytest_returncode": process.returncode,
        "pytest_tail": process.stdout.strip().splitlines()[-1] if process.stdout.strip() else "",
        "pytest_log": str(pytest_log),
        "gates_file": str(out_dir / "gates_verified.json"),
        "passed": process.returncode == 0,
        "evidence": "REAL (the runner runs the pytest module and keeps its log)",
    }
    gates.append(c6)

    figures = {
        "fig_c1_gap_recall.png": figure_gap_recall(c1, out_dir),
        "fig_c2_curvature.png": figure_curvature(wigley_mesh, out_dir),
        "fig_c3_c5_hydrostatics.png": figure_hydrostatics(cases, out_dir),
    }

    summary = {
        "wigley_iii": wigley_report.to_dict(),
        "series60_cb070": s60_report.to_dict(),
        "series60_fairness_note": (
            "The Series 60 parent form has a zero-radius bilge (sections are zeta**(1/m) with "
            "1/m -> 0), so the curvature field is singular along the keel line. The fairness "
            "gate reports that as a curvature discontinuity; it is a property of the parent "
            "form, not of the mesh, and a production hull needs a finite bilge radius."
        ),
        "exchange": exchange,
        "figures": {name: str(path) for name, path in figures.items()},
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    gates_path = write_gates_json(
        out_dir / "gates_verified.json", gates, mission=MISSION,
        extra={"cases": cases, "exchange": exchange, "figures": summary["figures"]},
    )
    write_gates_json(
        gates_path, gates, mission=MISSION,
        extra={"cases": cases, "exchange": exchange, "figures": summary["figures"]},
    )
    payload = json.loads(gates_path.read_text(encoding="utf-8"))
    self_check = payload["all_passed"] and {gate["id"] for gate in payload["gates"]} == {
        "C1", "C2", "C3", "C4", "C5", "C6"
    }
    c6["gates_file_written"] = gates_path.exists()
    c6["gates_file_self_check"] = self_check
    c6["passed"] = bool(c6["passed"] and self_check and c6["gates_file_written"])
    write_gates_json(
        gates_path, gates, mission=MISSION,
        extra={"cases": cases, "exchange": exchange, "figures": summary["figures"]},
    )
    print(json.dumps(
        {
            "gates_file": str(gates_path),
            "all_passed": all(bool(gate["passed"]) for gate in gates),
            "gates": {
                gate["id"]: {
                    "passed": bool(gate["passed"]),
                    **{key: gate[key] for key in ("recall", "worst_relative_difference",
                                                 "worst_relative_error", "cb_reported_deviation")
                       if key in gate},
                }
                for gate in gates
            },
        },
        indent=2,
    ))
    return 0 if all(bool(gate["passed"]) for gate in gates) else 1


if __name__ == "__main__":
    raise SystemExit(main())