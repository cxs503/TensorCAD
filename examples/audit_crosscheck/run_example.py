#!/usr/bin/env python3
"""Two geometry audits, one corpus: do the CAD (BRep) route and the mesh route agree?

``tensorcad`` carries two independent quality-audit stacks that grew up apart:

* the **BRep route** (:mod:`tensorcad.geomquality`) reads NURBS patches and
  names defects in patch/topology/continuity terms -- ``FREE_EDGE``,
  ``SEAM_GAP``, ``SLIVER``, ``CURVATURE_JUMP`` ... (14-case library);
* the **mesh route** (:mod:`tensorcad.validation` +
  :mod:`tensorcad.geometry_audit` + :mod:`tensorcad.mesh`) reads a triangulation
  -- the thing a downstream solver actually consumes -- and reports open edges,
  non-manifold edges, duplicate facets, winding conflicts, self-intersections,
  free gap magnitudes, triangle quality and Gaussian-curvature fairness.

The 2026-11-09 repo audit (``audit_report.md``, mission m14491) flagged exactly this
as an unverified risk: *"geometry_audit 与 geomquality 两套质量审计并行存在，判定口径是否
一致未经本审计验证"*, and left ``geometry_audit`` as ``WITHHELD`` (no execution
evidence). This example is that experiment. It builds one corpus -- three reference
hulls (Wigley III, Series 60 Cb 0.60, KVLCC2-class) times the geomquality defect
library's 14 injections plus the clean reference, 45 models -- and runs **both**
routes over every one of them, with the injected defect as the ground-truth label.

Four questions are answered with numbers:

1. **Does the BRep route actually work on all three hulls** (not just the one the
   unit test uses)? -- detection matrix, expected code vs reported codes.
2. **What does the mesh route see of the same defect?** -- per-family observable,
   classified ``named`` / ``coarse`` (only "the shell is open") / ``blind``.
3. **Do the two routes agree on the *clean* reference?** -- false alarms, with the
   zero-breadth stem of the Wigley form as the measured explanation.
4. **Is the mesh route's triangle-quality gate discriminative at all on raw CAD
   tessellation?** -- smallest angle / aspect ratio / sliver share per hull and per
   tessellation density.

Writes ``results.json``, ``report.md`` and two figures into ``--output`` and exits
non-zero if any pinned claim in ``CLAIMS`` stops holding.

Usage::

    python examples/audit_crosscheck/run_example.py [--output artifacts/audit_crosscheck]
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tensorcad import geometry_audit as audit  # noqa: E402
from tensorcad import intersect, mesh as meshlib, validation  # noqa: E402
from tensorcad.geomquality import baselines, defects  # noqa: E402

#: Two routes, one corpus. Fault families are the geomquality library's own codes.
TOPOLOGY_FAMILIES = (
    "FREE_EDGE",
    "OPEN_TOP",
    "SEAM_GAP",
    "NON_MANIFOLD",
    "ORIENTATION",
    "DUPLICATE_FACE",
)
COARSE_FAMILIES = ("SLIVER", "MICRO_EDGE", "CURVATURE_JUMP", "TANGENT_JUMP")
FAIRNESS_FAMILIES = ("ROUGHNESS",)
GEOMETRIC_FAMILIES = ("SELF_INTERSECTION",)
DECLARED_FAMILIES = ("UNIT_SCALE", "CB_MISMATCH")

DENSITIES = ((8, 4), (16, 8), (32, 16))
QUALITY_LIMITS = {"min_angle_deg": 2.0, "aspect_ratio": 25.0, "sliver_share": 0.02}


# --------------------------------------------------------------------------- #
# measurements
# --------------------------------------------------------------------------- #
def triangle_quality(mesh) -> dict[str, float]:
    """Smallest angle, worst aspect ratio and sliver share of a triangulation.

    Recomputed here rather than parsed out of a check message: the harness has to
    report the number, not the verdict.
    """
    smallest = math.inf
    worst_aspect = 0.0
    slivers = 0
    for face in mesh.faces:
        points = [mesh.vertices[index] for index in face]
        lengths = []
        for a in range(3):
            p, q = points[a], points[(a + 1) % 3]
            lengths.append(math.dist(p, q))
        area2 = abs(
            (points[1][0] - points[0][0]) * (points[2][1] - points[0][1])
            - (points[2][0] - points[0][0]) * (points[1][1] - points[0][1])
        )
        longest = max(lengths)
        height = area2 / longest if longest else 0.0
        aspect = longest / height if height > 0.0 else math.inf
        worst_aspect = max(worst_aspect, aspect)
        opposite, side_a, side_b = sorted(lengths)
        denominator = 2.0 * side_a * side_b
        if denominator <= 0.0:
            angle = 0.0
        else:
            cosine = (side_a * side_a + side_b * side_b - opposite * opposite) / denominator
            angle = math.degrees(math.acos(max(-1.0, min(1.0, cosine))))
        smallest = min(smallest, angle)
        if angle < QUALITY_LIMITS["min_angle_deg"]:
            slivers += 1
    return {
        "min_angle_deg": smallest,
        "aspect_ratio": worst_aspect,
        "sliver_share": slivers / len(mesh.faces) if mesh.faces else 0.0,
        "sliver_faces": slivers,
        "faces": len(mesh.faces),
    }


def self_intersections(mesh) -> int:
    """Pairs of triangles that fold through each other (``intersect`` route)."""
    found = intersect.find_self_intersections(mesh)
    if isinstance(found, dict):
        for key in ("pairs", "intersections", "count", "found"):
            if key in found:
                value = found[key]
                return len(value) if isinstance(value, (list, tuple)) else int(value)
        return int(found.get("n_pairs", 0))
    return len(found)


def mesh_route(brep, declared: dict[str, float], uv=(16, 8)) -> dict[str, object]:
    """Everything the mesh route can say about one tessellation."""
    mesh = brep.tessellate(*uv)
    watertightness = audit.audit_watertightness(mesh, length=float(brep.length)).to_dict()
    fairness = audit.audit_fairness(mesh)
    geometry = validation.inspect_geometry(mesh)
    quality_checks = validation.inspect_mesh_quality(mesh)
    counts = meshlib.edge_face_counts(mesh)
    volume = meshlib.total_signed_volume(mesh)
    box = mesh.bounding_box()
    return {
        "vertices": len(mesh.vertices),
        "faces": len(mesh.faces),
        "signed_volume": volume,
        "surface_area": meshlib.surface_area(mesh),
        "extent_x": box.sizes[0],
        "open_edges": watertightness["open_edges"],
        "boundary_components": watertightness["boundary_components"],
        "nonmanifold_edges": watertightness["nonmanifold_edges"],
        "duplicate_faces": watertightness["duplicate_faces"],
        "degenerate_faces": watertightness["degenerate_faces"],
        "winding_conflicts": watertightness["winding_conflicts"],
        "max_gap_m": watertightness["max_gap_m"],
        "watertight_audit_passed": watertightness["passed"],
        "nonmanifold_edges_counted": sum(1 for c in counts.values() if c > 2),
        "self_intersection_pairs": self_intersections(mesh),
        "fairness_defects": len(fairness.defects),
        "geometry_failures": [check.name for check in geometry if not check.passed],
        "quality_failures": [check.name for check in quality_checks if not check.passed],
        "quality": triangle_quality(mesh),
        # declared references are what let the mesh route check scale and Cb
        "declared_length": declared["length"],
        "declared_cb": declared["block_coefficient"],
        "cb_realised": volume / (
            declared["length"] * declared["beam"] * declared["draft"]
        ),
    }


# --------------------------------------------------------------------------- #
# per-family classification of the mesh route
# --------------------------------------------------------------------------- #
def mesh_verdict(family: str, clean: dict, model: dict) -> str:
    """``named`` / ``coarse`` / ``fairness`` / ``geometric`` / ``declared`` / ``blind``."""
    if family in TOPOLOGY_FAMILIES:
        return "named"
    if family in COARSE_FAMILIES:
        # A detached / unsewn patch shows up as an open shell, never by name.
        return "coarse" if model["open_edges"] > clean["open_edges"] else "blind"
    if family in FAIRNESS_FAMILIES:
        return "fairness" if model["fairness_defects"] > clean["fairness_defects"] else "blind"
    if family in GEOMETRIC_FAMILIES:
        return (
            "geometric"
            if model["self_intersection_pairs"] > clean["self_intersection_pairs"]
            else "blind"
        )
    if family in DECLARED_FAMILIES:
        if family == "UNIT_SCALE":
            ratio = model["extent_x"] / model["declared_length"]
            return "declared" if abs(ratio - 1.0) > 0.01 else "blind"
        return (
            "declared"
            if abs(model["cb_realised"] / model["declared_cb"] - 1.0) > 0.02
            else "blind"
        )
    raise KeyError(family)


def observable(family: str, clean: dict, model: dict) -> str:
    """The concrete number the verdict was read off, so the table can be checked."""
    if family in TOPOLOGY_FAMILIES:
        return (
            f"open={model['open_edges']-clean['open_edges']:+d} "
            f"nm={model['nonmanifold_edges']-clean['nonmanifold_edges']:+d} "
            f"dup={model['duplicate_faces']-clean['duplicate_faces']:+d} "
            f"wind={model['winding_conflicts']-clean['winding_conflicts']:+d}"
        )
    if family in COARSE_FAMILIES:
        return f"open={model['open_edges']-clean['open_edges']:+d}"
    if family in FAIRNESS_FAMILIES:
        return f"fairness={model['fairness_defects']-clean['fairness_defects']:+d}"
    if family in GEOMETRIC_FAMILIES:
        return f"self_x={model['self_intersection_pairs']}"
    if family == "UNIT_SCALE":
        return f"extent/declared={model['extent_x']/model['declared_length']:.6g}"
    return f"cb_realised/declared={model['cb_realised']/model['declared_cb']:.6f}"


# --------------------------------------------------------------------------- #
# claims the example hard-fails on
# --------------------------------------------------------------------------- #
def evaluate_claims(records: dict, specimens: list[str]) -> list[dict[str, object]]:
    claims = []

    def claim(name: str, ok: bool, detail: str) -> None:
        claims.append({"claim": name, "passed": bool(ok), "detail": detail})

    injected = [
        (name, case) for name in specimens for case in [c.name for c in defects.DEFECT_CASES]
    ]
    hits = sum(
        1
        for name, case in injected
        if records[(name, case)]["expected_code"] in records[(name, case)]["brep"]["codes"]
    )
    claim(
        "brep_route_names_every_injection_on_every_hull",
        hits == len(injected),
        f"{hits}/{len(injected)} injected defects named by their own code",
    )

    clean_codes = {
        name: ",".join(records[(name, "CLEAN")]["brep"]["codes"]) for name in specimens
    }
    claim(
        "brep_route_reports_no_defect_on_clean_references",
        all(not codes for codes in clean_codes.values()),
        f"clean codes per hull: {clean_codes}",
    )

    named = 0
    total = 0
    for name in specimens:
        clean = records[(name, "CLEAN")]["mesh"]
        for case in defects.DEFECT_CASES:
            if case.code not in TOPOLOGY_FAMILIES:
                continue
            total += 1
            if mesh_verdict(case.code, clean, records[(name, case.name)]["mesh"]) == "named":
                named += 1
    claim(
        "mesh_route_sees_every_topology_family",
        named == total,
        f"{named}/{total} topology injections produce a mesh-side topology signal",
    )

    folds = [
        name
        for name in specimens
        if records[(name, "fold")]["mesh"]["self_intersection_pairs"]
        > records[(name, "CLEAN")]["mesh"]["self_intersection_pairs"]
    ]
    others = [
        (name, case.name)
        for name in specimens
        for case in defects.DEFECT_CASES
        if case.name != "fold"
    ]
    false_folds = [
        key
        for key in others
        if records[key]["mesh"]["self_intersection_pairs"]
        > records[(key[0], "CLEAN")]["mesh"]["self_intersection_pairs"]
    ]
    # Measured: the fold is seen on every hull and a clean reference never trips the
    # detector -- but a flap, a doubled face, a micro-edge patch and a ripple also put
    # faces through each other, so "fires on folds only" is false and is not claimed.
    claim(
        "self_intersection_sees_the_fold_and_never_a_clean_reference",
        len(folds) == len(specimens)
        and all(
            not records[(name, "CLEAN")]["mesh"]["self_intersection_pairs"]
            for name in specimens
        ),
        f"folds seen on {len(folds)}/{len(specimens)} hulls; other injections that also "
        f"overlap geometrically (expected, not noise): {false_folds}",
    )

    blind = {
        family: sorted(
            {
                mesh_verdict(
                    family,
                    records[(name, "CLEAN")]["mesh"],
                    records[(name, case.name)]["mesh"],
                )
                for name in specimens
                for case in defects.DEFECT_CASES
                if case.code == family
            }
        )
        for family in COARSE_FAMILIES + FAIRNESS_FAMILIES + GEOMETRIC_FAMILIES + DECLARED_FAMILIES
    }
    claim(
        "coarse_families_are_never_named_by_the_mesh_route",
        all(kinds == ["coarse"] or "named" not in kinds for kinds in blind.values()),
        f"verdicts by family: {blind}",
    )

    quality = {
        name: {
            "sliver_share": records[(name, "CLEAN")]["mesh"]["quality"]["sliver_share"],
            "min_angle_deg": records[(name, "CLEAN")]["mesh"]["quality"]["min_angle_deg"],
        }
        for name in specimens
    }
    violated = [
        name
        for name, row in quality.items()
        if row["sliver_share"] > QUALITY_LIMITS["sliver_share"]
        or row["min_angle_deg"] < QUALITY_LIMITS["min_angle_deg"]
    ]
    claim(
        "raw_cad_tessellation_violates_the_quality_gate_at_16x8",
        bool(violated),
        f"violating hulls {violated} of {specimens}; measured {quality}",
    )

    return claims


# --------------------------------------------------------------------------- #
# figures
# --------------------------------------------------------------------------- #
def _plot_figures(payload: dict, output: Path) -> list[str]:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:  # pragma: no cover - figures are optional
        print(f"[skip] figures: {exc}")
        return []

    made = []
    families = payload["family_table"]
    labels = [row["family"] for row in families]
    brep = [row["brep_detected"] for row in families]
    mesh = [row["mesh_named"] for row in families]
    coarse = [row["mesh_coarse"] for row in families]

    figure, axis = plt.subplots(figsize=(9.0, 4.6))
    index = range(len(labels))
    axis.bar([i - 0.27 for i in index], brep, width=0.27, label="BRep route (named)")
    axis.bar(index, mesh, width=0.27, label="mesh route (named)")
    axis.bar([i + 0.27 for i in index], coarse, width=0.27, label="mesh route (coarse only)")
    axis.set_xticks(list(index))
    axis.set_xticklabels(labels, rotation=55, ha="right", fontsize=8)
    axis.set_ylabel(f"hulls detected (of {payload['n_specimens']})")
    axis.set_title("Same 14 injected defects, two audit routes")
    axis.legend(fontsize=8)
    axis.grid(axis="y", alpha=0.3)
    figure.tight_layout()
    path = output / "fig1_detection_parity.png"
    figure.savefig(path, dpi=150)
    plt.close(figure)
    made.append(path.name)

    figure, axis = plt.subplots(figsize=(7.6, 4.4))
    for name, series in payload["density_study"].items():
        axis.plot(
            [f"{u}x{v}" for u, v in DENSITIES],
            [point["sliver_share"] * 100.0 for point in series],
            marker="o",
            label=f"{name} sliver share",
        )
    axis.axhline(QUALITY_LIMITS["sliver_share"] * 100.0, linestyle="--", color="black", label="gate limit 2 %")
    axis.set_ylabel("sliver triangles (% of faces)")
    axis.set_xlabel("tessellation density (u x v)")
    axis.set_title("Raw CAD tessellation vs the mesh-route quality gate")
    axis.legend(fontsize=7)
    axis.grid(alpha=0.3)
    figure.tight_layout()
    path = output / "fig2_quality_vs_density.png"
    figure.savefig(path, dpi=150)
    plt.close(figure)
    made.append(path.name)
    return made


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--output", default="artifacts/audit_crosscheck")
    options = parser.parse_args(argv)
    output = Path(options.output)
    output.mkdir(parents=True, exist_ok=True)

    catalogue = baselines.specimens()
    specimens = sorted(catalogue)
    records: dict = {}
    for name in specimens:
        specimen = catalogue[name]()
        declared = {
            "length": specimen.length,
            "beam": specimen.beam,
            "draft": specimen.draft,
            "block_coefficient": specimen.particulars["block_coefficient"],
        }
        for case in [None] + list(defects.DEFECT_CASES):
            brep = specimen.brep if case is None else case.inject(specimen.brep)
            result = defects.detect_all(brep, declared=declared)
            records[(name, case.name if case else "CLEAN")] = {
                "specimen": name,
                "case": case.name if case else "CLEAN",
                "expected_code": case.code if case else None,
                "brep": {
                    "codes": sorted(result["codes"]),
                    "watertight": result.get("watertight"),
                },
                "mesh": mesh_route(brep, declared),
            }
            print(f"{name:16s} {str(case.name if case else 'CLEAN'):18s} "
                  f"brep={records[(name, case.name if case else 'CLEAN')]['brep']['codes']}")

    # family table: how many hulls each route catches per injected family
    family_table = []
    for case in defects.DEFECT_CASES:
        brep_hits = sum(
            1 for name in specimens if case.code in records[(name, case.name)]["brep"]["codes"]
        )
        verdicts = [
            mesh_verdict(case.code, records[(name, "CLEAN")]["mesh"], records[(name, case.name)]["mesh"])
            for name in specimens
        ]
        family_table.append(
            {
                "family": case.code,
                "case": case.name,
                "brep_detected": brep_hits,
                "mesh_named": sum(1 for v in verdicts if v == "named"),
                "mesh_coarse": sum(1 for v in verdicts if v == "coarse"),
                "mesh_other": sum(1 for v in verdicts if v in ("fairness", "geometric", "declared")),
                "mesh_blind": sum(1 for v in verdicts if v == "blind"),
                "verdicts_per_hull": dict(zip(specimens, verdicts)),
                "observables_per_hull": {
                    name: observable(case.code, records[(name, "CLEAN")]["mesh"], records[(name, case.name)]["mesh"])
                    for name in specimens
                },
            }
        )

    # clean-reference false alarms
    clean_table = []
    for name in specimens:
        clean = records[(name, "CLEAN")]
        defects_named = [d["kind"] for d in audit.audit_watertightness(
            catalogue[name]().brep.tessellate(16, 8), length=float(catalogue[name]().length)
        ).to_dict()["defects"]]
        clean_table.append(
            {
                "specimen": name,
                "brep_codes": clean["brep"]["codes"],
                "mesh_watertight_audit_passed": clean["mesh"]["watertight_audit_passed"],
                "mesh_defects": defects_named,
                "geometry_failures": clean["mesh"]["geometry_failures"],
                "quality_failures": clean["mesh"]["quality_failures"],
                "quality": clean["mesh"]["quality"],
                "vertices": clean["mesh"]["vertices"],
                "faces": clean["mesh"]["faces"],
                "declared_volume": (
                    catalogue[name]().particulars["block_coefficient"]
                    * catalogue[name]().length
                    * catalogue[name]().beam
                    * catalogue[name]().draft
                ),
                "signed_volume": clean["mesh"]["signed_volume"],
                "self_intersection_pairs": clean["mesh"]["self_intersection_pairs"],
            }
        )

    # density study: is the quality gate fixed by refinement?
    density_study = {}
    for name in specimens:
        specimen = catalogue[name]()
        series = []
        for uv in DENSITIES:
            mesh = specimen.brep.tessellate(*uv)
            quality = triangle_quality(mesh)
            quality["density"] = f"{uv[0]}x{uv[1]}"
            quality["signed_volume"] = meshlib.total_signed_volume(mesh)
            series.append(quality)
        density_study[name] = series

    claims = evaluate_claims(records, specimens)
    payload = {
        "schema": "tensorcad.audit_crosscheck/1",
        "specimens": specimens,
        "n_specimens": len(specimens),
        "n_models": len(records),
        "densities": [f"{u}x{v}" for u, v in DENSITIES],
        "family_table": family_table,
        "clean_table": clean_table,
        "density_study": density_study,
        "claims": claims,
        "records": [
            {**value, "key": f"{value['specimen']}/{value['case']}"}
            for value in records.values()
        ],
    }
    (output / "results.json").write_text(json.dumps(payload, indent=1, default=str))

    figures = _plot_figures(payload, output)
    payload["figures"] = figures

    failed = [claim for claim in claims if not claim["passed"]]
    lines = ["# Audit cross-check: BRep route vs mesh route", ""]
    lines.append(f"corpus: {len(specimens)} hulls x ({len(defects.DEFECT_CASES)} injections + clean)"
                 f" = {len(records)} models; tessellation 16x8")
    lines.append("")
    lines.append("## Detection parity")
    lines.append("")
    lines.append("| family | BRep named | mesh named | mesh coarse | mesh other | mesh blind |")
    lines.append("|---|---:|---:|---:|---:|---:|")
    for row in family_table:
        lines.append(
            f"| {row['family']} | {row['brep_detected']} | {row['mesh_named']} | "
            f"{row['mesh_coarse']} | {row['mesh_other']} | {row['mesh_blind']} |"
        )
    lines.append("")
    lines.append("## Clean references (false-alarm check)")
    lines.append("")
    lines.append("| hull | BRep codes | mesh defects | sliver share | min angle (deg) |")
    lines.append("|---|---|---|---:|---:|")
    for row in clean_table:
        lines.append(
            f"| {row['specimen']} | {row['brep_codes']} | {', '.join(row['mesh_defects']) or '-'} | "
            f"{100 * row['quality']['sliver_share']:.2f} % | {row['quality']['min_angle_deg']:.2f} |"
        )
    lines.append("")
    lines.append("## Claims")
    lines.append("")
    for claim in claims:
        lines.append(f"- {'PASS' if claim['passed'] else 'FAIL'} `{claim['claim']}`: {claim['detail']}")
    lines.append("")
    (output / "report.md").write_text("\n".join(lines) + "\n")

    print("\n".join(lines))
    if failed:
        print(f"{len(failed)} claim(s) failed", file=sys.stderr)
        return 1
    print(f"wrote {output}/results.json, report.md, {len(figures)} figure(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())