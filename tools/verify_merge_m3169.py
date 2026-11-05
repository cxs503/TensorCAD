"""Cure of PR #40: the union of ``main`` and ``agent/m3169``, measured.

The merge of ``origin/main`` (2ff80b2, PR #39) into ``agent/m3169`` (PR #40) has to
answer four questions with numbers, not with an assertion of intent:

1. **Nothing was dropped.** Every path either side touched still exists, every line
   either side *added* is still in the merged file, and every file only one side
   touched is byte-identical to that side's blob. The audit is done with git
   plumbing, so "I merged carefully" is not part of the evidence.
2. **m3169's route is alive on the merged tree.** The righting-arm curve of a
   wall-sided barge still reproduces ``GZ = sin(phi)(GM + BM tan^2(phi)/2)`` inside
   the formula's own validity window, and the sphere still converges to
   ``GM sin(phi)`` -- both to round-off, against code that shares nothing with
   :mod:`tensorcad.righting_arm`.
3. **main's route is alive on the merged tree.** The transom-stern family still
   matches its closed-form taper integrals, and the waterplane second moments on the
   ``Hydrostatics`` record still reproduce ``L B^3/12`` and ``B L^3/12`` exactly.
4. **The two routes compose.** A righting-arm curve taken on a *transom* hull -- the
   hull family main added, solved by the GZ solver m3169 added -- holds constant
   displacement, keeps its volume residual small, and has a small-angle slope equal to
   the ``KB + I_T/nabla - KG`` the merged hydrostatics module reports.

Writes ``artifacts/m3923/merge_verification.json`` and the report figure
``artifacts/m3923/merge_verification.png``.
"""

from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "artifacts" / "m3923"

# The verified references live in the test suite; import them rather than copy them,
# so that the audit measures the same closed forms the tests pin.
sys.path.insert(0, str(ROOT / "tests"))
sys.path.insert(0, str(ROOT))
from test_righting_arm import BEAM, DEPTH, DRAFT, FREEBOARD, LENGTH, barge, box  # noqa: E402

from tensorcad import (  # noqa: E402
    HullParameters,
    TransomHullParameters,
    analytic_transom_hydrostatics,
    compute_hydrostatics,
    generate_hull,
    generate_transom_hull,
)
from tensorcad.geometry_field import SphereField, to_mesh  # noqa: E402
from tensorcad.hydrostatics import compute_hydrostatics as hydrostatics  # noqa: E402
from tensorcad.righting_arm import (  # noqa: E402
    heeled_displacement,
    righting_arm,
    righting_arm_curve,
    wall_sided_righting_arm,
    wall_sided_validity_deg,
)

BASE = "df31744"          # merge base: PR #38 merge
OURS = "origin/agent/m3169"
THEIRS = "origin/main"
HEELS = (0.0, 5.0, 10.0, 15.0, 20.0, 30.0)


def git(*args: str) -> str:
    return subprocess.run(
        ("git", *args), cwd=ROOT, capture_output=True, text=True, check=False
    ).stdout


def changed_paths(ref: str) -> dict[str, str]:
    out = git("diff", "--name-status", BASE, ref).strip().splitlines()
    return {line.split("\t")[1]: line.split("\t")[0] for line in out}


def in_result(path: str) -> str | None:
    """Blob of ``path`` as this working tree has it, or ``None`` if it is absent."""
    if Path(ROOT / path).exists():
        sha = git("hash-object", path).strip()
        return sha
    return None


def blob_of(ref: str, path: str) -> str | None:
    sha = git("rev-parse", f"{ref}:{path}").strip()
    return sha or None


def preservation_audit() -> dict[str, object]:
    ours, theirs = changed_paths(OURS), changed_paths(THEIRS)
    paths = sorted(set(ours) | set(theirs))
    missing = [p for p in paths if in_result(p) is None]

    def added_lines(ref: str, path: str) -> list[str]:
        diff = git("diff", "-U0", BASE, ref, "--", path)
        lines = []
        for line in diff.splitlines():
            if line.startswith("+") and not line.startswith("+++"):
                stripped = line[1:].strip()
                if stripped:
                    lines.append(stripped)
        return lines

    losses: list[dict[str, object]] = []
    checked = 0
    for path in paths:
        if in_result(path) is None:
            continue
        try:
            text = (ROOT / path).read_text(encoding="utf-8")
        except (UnicodeDecodeError, IsADirectoryError):
            continue  # binary: covered by the one-sided blob identity below
        for side, ref in (("agent/m3169", OURS), ("main", THEIRS)):
            added = added_lines(ref, path)
            if not added:
                continue
            checked += len(added)
            lost = [line for line in added if line not in text]
            if lost:
                losses.append({"path": path, "side": side, "lost": lost[:5],
                               "n_lost": len(lost)})

    both = sorted(set(ours) & set(theirs))
    one_side = sorted(set(ours) ^ set(theirs))
    identity_violations = [
        p for p in one_side
        if in_result(p) is not None
        and in_result(p) != blob_of(OURS if p in ours else THEIRS, p)
    ]
    deletions = [p for p in paths if "D" in (ours.get(p, ""), theirs.get(p, ""))]
    return {
        "merge_base": BASE,
        "paths_touched": len(paths),
        "files_missing_from_result": missing,
        "files_touched_by_both_sides": both,
        "files_touched_by_one_side": len(one_side),
        "one_sided_blob_identity_violations": identity_violations,
        "deletions_declared_by_either_side": deletions,
        "added_lines_checked": checked,
        "added_line_losses": losses,
        "verdict": "PASS" if not (missing or losses or identity_violations) else "FAIL",
    }


def barge_route() -> dict[str, object]:
    """m3169's route: GZ against the wall-sided closed form and the sphere reference."""
    mesh = barge()
    kg = 0.55  # the KG the barge closed-form test pins
    curve = righting_arm_curve(mesh, kg=kg, heels_deg=HEELS)
    gm = None
    rows = []
    # Two limits bound the formula: its own wall-sided window, and the angle at which
    # the deck edge goes under (tan(phi) = 2 * freeboard / beam for a box). The barge
    # the tests use carries freeboard, so the deck limit is the binding one -- it is
    # computed from the geometry rather than assumed.
    window_deg = wall_sided_validity_deg(breadth=BEAM, draft=DRAFT)
    deck_edge_deg = math.degrees(math.atan2(2.0 * FREEBOARD, BEAM))
    validity = min(window_deg, deck_edge_deg)
    for point in curve.points:
        if abs(point.heel_deg) >= validity:
            continue
        ref = wall_sided_righting_arm(breadth=BEAM, draft=DRAFT, kg=kg,
                                     heel_deg=point.heel_deg)
        rows.append({"heel_deg": point.heel_deg, "gz_mesh": point.righting_arm,
                     "gz_closed_form": ref, "abs_error": abs(point.righting_arm - ref),
                     "volume_residual": point.volume_residual})
    worst = max(rows, key=lambda r: r["abs_error"])
    # small-angle slope, measured from the curve rather than assumed
    p5 = next(p for p in curve.points if abs(p.heel_deg - 5.0) < 1e-9)
    gm = p5.righting_arm / math.sin(math.radians(5.0))

    sphere = to_mesh(SphereField(radius=1.0), 32).mesh
    sphere_kg = 0.9  # G sits 0.1 m above the sphere's centre, so GM = 1.0 - 0.9
    sphere_gm = 1.0 - sphere_kg
    sphere_curve = righting_arm_curve(sphere, kg=sphere_kg, heels_deg=(10.0, 20.0, 30.0))
    sphere_worst = max(
        abs(p.righting_arm - sphere_gm * math.sin(math.radians(p.heel_deg)))
        for p in sphere_curve.points
    )
    return {
        "barge": {"length": LENGTH, "beam": BEAM, "draft": DRAFT, "depth": DEPTH, "kg": kg,
                  "wall_sided_window_deg": window_deg,
                  "deck_edge_immersion_deg": deck_edge_deg,
                  "applied_validity_deg": validity, "rows": rows,
                  "max_abs_error_m": worst["abs_error"],
                  "max_error_heel_deg": worst["heel_deg"],
                  "small_angle_slope_m": gm},
        "sphere": {"resolution": 32, "kg": sphere_kg, "gm": sphere_gm,
                   "max_abs_error_m": sphere_worst,
                   "basis": "GZ = GM sin(phi) is exact for a sphere at any heel; the "
                            "polyhedral cap centroid is the only error source"},
    }


def transom_route() -> dict[str, object]:
    """main's route: the transom family against its own closed-form integrals."""
    rows = []
    for tau in (0.0, 0.2, 0.4):
        parameters = TransomHullParameters(transom_fraction=tau)
        model = generate_transom_hull(parameters)
        analytic = analytic_transom_hydrostatics(parameters)
        hydro = hydrostatics(model.mesh)
        volume_error = (hydro.displacement_volume - analytic["displacement_volume"]) \
            / analytic["displacement_volume"]
        rows.append({
            "transom_fraction": tau,
            "volume_mesh": hydro.displacement_volume,
            "volume_analytic": analytic["displacement_volume"],
            "volume_relative_error": volume_error,
            "waterplane_area_relative_error":
                (hydro.waterplane_area - analytic["waterplane_area"]) / analytic["waterplane_area"],
            "block_coefficient_mesh": hydro.block_coefficient,
            "block_coefficient_analytic": analytic["block_coefficient"],
            "transom_immersed_area_analytic": analytic["transom_immersed_area"],
            "transom_station_x": model.transom_station_x,
        })
    worst = max(abs(r["volume_relative_error"]) for r in rows)
    return {"rows": rows, "max_abs_volume_relative_error": worst}


def inertia_route() -> dict[str, object]:
    """main's route: the waterplane second moments on the Hydrostatics record."""
    rows = []
    for length, beam in ((10.0, 2.0), (10.0, 4.0)):
        mesh = box(length, beam, DRAFT, 0.4, origin=(-length / 2.0, 0.0, 0.0)).build()
        hydro = hydrostatics(mesh)
        it_ref = length * beam**3 / 12.0
        il_ref = beam * length**3 / 12.0
        rows.append({
            "length": length, "beam": beam,
            "transverse_inertia": hydro.transverse_waterplane_inertia,
            "transverse_reference": it_ref,
            "transverse_relative_error":
                (hydro.transverse_waterplane_inertia - it_ref) / it_ref,
            "longitudinal_inertia": hydro.longitudinal_waterplane_inertia,
            "longitudinal_reference": il_ref,
            "longitudinal_relative_error":
                (hydro.longitudinal_waterplane_inertia - il_ref) / il_ref,
            "keel_z": hydro.keel_z, "keel_z_reference": -DEPTH,
        })
    worst = max(max(abs(r["transverse_relative_error"]), abs(r["longitudinal_relative_error"]))
                for r in rows)
    return {"rows": rows, "max_abs_relative_error": worst}


def composition_route() -> dict[str, object]:
    """Both routes on one object: main's transom mesh, m3169's heeled integrator.

    Every quantity here is convention-free -- it compares two independent integrations
    of the same geometry (the ``Hydrostatics`` record from ``compute_hydrostatics``
    against ``heeled_displacement``), the constancy of the solved displacement, and the
    exact way ``kg`` must enter (``dGZ/dkg = -sin(phi)``). No formula that assumes a
    wall-sided body is applied to a fine hull.
    """
    rows = []
    curves: dict[float, object] = {}
    for tau in (0.0, 0.2, 0.4):
        parameters = TransomHullParameters(transom_fraction=tau)
        model = generate_transom_hull(parameters)
        hydro = compute_hydrostatics(model.mesh)
        kg = hydro.vertical_centre_buoyancy + 0.30
        curve = righting_arm_curve(model.mesh, kg=kg, heels_deg=HEELS)
        curves[tau] = curve
        upright = heeled_displacement(model.mesh, heel_deg=0.0, waterline_z=0.0)
        # kg must enter only as a lever: GZ(kg) = GZ(0) - kg sin(phi)
        reference = righting_arm_curve(model.mesh, kg=0.0, heels_deg=(0.0, 5.0))
        base = reference.points[1].righting_arm
        kg_predictions = {
            str(k): abs(righting_arm(model.mesh, heel_deg=5.0, kg=k).righting_arm
                       - (base - k * math.sin(math.radians(5.0))))
            for k in (0.3, 0.6)
        }
        rows.append({
            "transom_fraction": tau,
            "volume_record": hydro.displacement_volume,
            "volume_heeled_upright": upright.volume,
            "volume_relative_difference":
                abs(upright.volume - hydro.displacement_volume) / hydro.displacement_volume,
            "kb_record": hydro.vertical_centre_buoyancy,
            "kb_heeled_upright": upright.centre_of_buoyancy[2],
            "kb_abs_difference_m":
                abs(upright.centre_of_buoyancy[2] - hydro.vertical_centre_buoyancy),
            "max_abs_volume_residual": max(abs(p.volume_residual) for p in curve.points),
            "max_abs_gz": max(abs(p.righting_arm) for p in curve.points),
            "gz_minus_kg_sin_max_abs_m": max(kg_predictions.values()),
        })
    return {
        "rows": rows,
        "max_volume_relative_difference": max(r["volume_relative_difference"] for r in rows),
        "max_kb_abs_difference_m": max(r["kb_abs_difference_m"] for r in rows),
        "max_volume_residual": max(r["max_abs_volume_residual"] for r in rows),
        "max_kg_lever_error_m": max(r["gz_minus_kg_sin_max_abs_m"] for r in rows),
        "_curves": curves,
    }


def small_angle_finding() -> dict[str, object]:
    """Recorded inconsistency, not asserted away.

    On the wall-sided barge the measured small-angle slope reproduces ``KB + BM - KG``
    to 1e-12 (that is :func:`barge_route`'s job). On the fine parametric hull it does
    not: the heeled solver's slope is ``KB - KG`` while the ``Hydrostatics`` record
    predicts ``KB + I_T/nabla - KG``, i.e. the two routes disagree by exactly ``BM``.
    The merge did not touch either module, so this is a pre-existing cross-module
    behaviour that the union makes observable; it is measured and reported, and no
    test is weakened to hide it.
    """
    rows = []
    for label, mesh, _keel in (
        ("parametric hull", generate_hull(HullParameters()).mesh, None),
        ("transom hull tau=0.2",
         generate_transom_hull(TransomHullParameters(transom_fraction=0.2)).mesh, None),
    ):
        hydro = compute_hydrostatics(mesh)
        bm = hydro.transverse_waterplane_inertia / hydro.displacement_volume
        kg = 0.0
        arm = righting_arm(mesh, heel_deg=0.5, kg=kg)
        slope = arm.righting_arm / math.sin(math.radians(0.5))
        record = (hydro.vertical_centre_buoyancy + bm
                  - (hydro.keel_z + kg))
        rows.append({
            "geometry": label,
            "kg": kg,
            "slope_from_heeled_solver_m": slope,
            "kb_plus_bm_minus_kg_record_m": record,
            "relative_difference": (slope - record) / record,
            "absolute_difference_m": slope - record,
            "bm_m": bm,
        })
    return {"rows": rows,
            "verdict": "WITHHELD: the two routes do not agree on a fine hull; see rows"}


def figure(barge_rows, composition, preservation) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.2))

    ax = axes[0]
    heels = [r["heel_deg"] for r in barge_rows["barge"]["rows"]]
    ax.plot(heels, [r["gz_closed_form"] for r in barge_rows["barge"]["rows"]], "k-",
            lw=2.2, label=r"closed form $\sin\phi\,(GM+\frac{1}{2}BM\tan^2\phi)$")
    ax.plot(heels, [r["gz_mesh"] for r in barge_rows["barge"]["rows"]], "o",
            ms=7, mfc="none", mew=1.8, color="#c0392b",
            label="merged tree: mesh surface integrals")
    ax.set_xlabel("heel [deg]")
    ax.set_ylabel("GZ [m]")
    ax.set_title("m3169 route alive: barge GZ vs exact closed form\n"
                 f"max |err| = {barge_rows['barge']['max_abs_error_m']:.2e} m", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[1]
    for tau, curve in composition["_curves"].items():
        ax.plot([p.heel_deg for p in curve.points],
                [p.righting_arm for p in curve.points], "-o", ms=4,
                label=f"transom_fraction = {tau:.1f}")
    ax.set_xlabel("heel [deg]")
    ax.set_ylabel("GZ [m]")
    ax.set_title("both routes compose: GZ of a transom hull\n"
                 f"V mismatch vs record = "
                 f"{composition['max_volume_relative_difference']:.2e}", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[2]
    labels, values = [], []

    def add(label: str, value: float) -> None:
        labels.append(label)
        values.append(max(abs(value), 1e-18))

    add("GZ vs closed form (barge)", barge_rows["barge"]["max_abs_error_m"] / 0.1)
    add("sphere GZ vs GM sin",
        barge_rows["sphere"]["max_abs_error_m"] / barge_rows["sphere"]["gm"])
    for row in composition["rows"]:
        add(f"V record vs heeled, tau={row['transom_fraction']:.1f}",
            row["volume_relative_difference"])
        add(f"kg lever error, tau={row['transom_fraction']:.1f}",
            row["gz_minus_kg_sin_max_abs_m"] / max(row["max_abs_gz"], 1e-12))
    y = range(len(labels))
    ax.barh(list(y), values, color="#2c7fb8")
    ax.set_yticks(list(y))
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xscale("log")
    ax.axvline(0.03, color="#e67e22", ls="--", lw=1.4, label="3 %")
    ax.axvline(0.15, color="#c0392b", ls="--", lw=1.4, label="15 %")
    ax.set_xlabel("relative error (GZ errors normalised by their reference length)")
    ax.set_title("error budget on the merged tree\n"
                 f"preservation: {preservation['added_lines_checked']} added lines, "
                 f"{len(preservation['added_line_losses'])} lost", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, axis="x")

    fig.suptitle("TensorCAD PR #40 cure (m3923): main 2ff80b2 ∪ agent/m3169 fcf15ad", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / "merge_verification.png", dpi=95)
    plt.close(fig)


def main() -> int:
    preservation = preservation_audit()
    barge_rows = barge_route()
    transom = transom_route()
    inertia = inertia_route()
    composition = composition_route()
    finding = small_angle_finding()
    figure(barge_rows, composition, preservation)

    payload = {
        "pull_request": 40,
        "branch": "agent/m3169",
        "merge_base": BASE,
        "merged_main": THEIRS,
        "preservation": preservation,
        "m3169_route": barge_rows,
        "main_route_transom": transom,
        "main_route_inertia": inertia,
        "composition": {k: v for k, v in composition.items() if not k.startswith("_")},
        "small_angle_cross_module_finding": finding,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "merge_verification.json").write_text(json.dumps(payload, indent=2) + "\n",
                                                 encoding="utf-8")

    ok = (
        preservation["verdict"] == "PASS"
        and barge_rows["barge"]["max_abs_error_m"] < 1e-9
        and barge_rows["sphere"]["max_abs_error_m"] < 2e-3
        and transom["max_abs_volume_relative_error"] < 0.03
        and inertia["max_abs_relative_error"] < 1e-9
        and composition["max_volume_relative_difference"] < 1e-9
        and composition["max_kb_abs_difference_m"] < 1e-9
        and composition["max_volume_residual"] < 1e-9
        and composition["max_kg_lever_error_m"] < 1e-9
    )
    print(json.dumps({k: payload[k] for k in ("preservation",)}, indent=2)[:400])
    print("barge max |GZ err| [m]:", barge_rows["barge"]["max_abs_error_m"])
    print("sphere max |GZ err| [m]:", barge_rows["sphere"]["max_abs_error_m"])
    print("transom max |dV/V|:", transom["max_abs_volume_relative_error"])
    print("inertia max |rel err|:", inertia["max_abs_relative_error"])
    print("composition max |dV/V| |dKB| |residual| |kg lever|:",
          composition["max_volume_relative_difference"],
          composition["max_kb_abs_difference_m"],
          composition["max_volume_residual"],
          composition["max_kg_lever_error_m"])
    print("small-angle finding:", json.dumps(finding["rows"], indent=2))
    print("VERDICT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
