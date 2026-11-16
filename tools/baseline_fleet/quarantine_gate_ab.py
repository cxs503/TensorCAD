"""Measure what the patch-quarantine gate changes, and how often the hole is live.

Two measurements, both reproducible from this repository alone:

**A -- controlled A/B (synthetic barge).**  A four-sided barge (L = 40 m, B = 10 m,
T = 4 m, so 1600 m3) is described by three type-128 patches, one of which is
deliberately undecidable (a 2x2 control net with non-unit weights: the weight test is
silent and a 2x2 net has no interior control point, so both candidate readings score
0.0 and the reader must fail closed).  The same file is meshed and graded twice,
changing exactly one thing: whether the builder's patch census reaches
:func:`tensorcad.fleet.gates.evaluate_gates`.  The undecidable patch overlaps surface the
other two patches already carry, so the half-breadth field -- and therefore the volume --
is identical in both builds: G7 cannot see the loss, which is the fail-open G9 closes.

**B -- census of the committed fleet evidence.**  The quarantined-patch count recorded in
``benchmarks/verified/baseline_fleet/geometry_diagnostics.json`` is cross-tabulated with
the verdicts recorded in ``benchmarks/verified/baseline_fleet/result.json``, i.e. whether
any hull of the four-hull fleet was *passed* while carrying quarantined patches.  That is
the honest measure of how live the hole is on today's evidence, and it is reported as it
comes out.

Usage::

    python tools/baseline_fleet/quarantine_gate_ab.py --out runs/m22653

Writes ``quarantine_gate_ab.json`` and (if matplotlib is importable)
``quarantine_gate_ab.png`` into ``--out``.
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

import tensorcad.iges  # noqa: E402
from tensorcad.fleet.build import iges_field_mesh  # noqa: E402
from tensorcad.fleet.gates import evaluate_gates, measure, quarantine_census  # noqa: E402
from tensorcad.fleet.manifest import FleetManifest, Provenance, ShipRecord  # noqa: E402
from tensorcad.iges import BSplineSurface  # noqa: E402

L, B2, T = 40.0, 5.0, 4.0
BARGE_VOLUME = 2.0 * B2 * L * T


def _patch(x0: float, x1: float, weights: list[float]) -> BSplineSurface:
    coords = [(x0, B2, 0.0), (x1, B2, 0.0), (x0, B2, T), (x1, B2, T)]
    params: list[float] = [128.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                           0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0]
    params += list(weights)
    for point in coords:
        params.extend(point)
    return BSplineSurface.from_parameters(params)


class _FakeIges:
    def __init__(self, surfaces: list[BSplineSurface]) -> None:
        self._surfaces = surfaces

    def surfaces(self) -> list[BSplineSurface]:
        return list(self._surfaces)

    def trimmed_surfaces(self) -> list[tuple[object, object, list[object]]]:
        return []


def measure_ab(tmp: Path) -> dict[str, object]:
    surfaces_with = [_patch(0.0, 20.0, [1.0] * 4), _patch(0.0, 20.0, [1.0, 1.06, 0.94, 1.0]),
                     _patch(20.0, 40.0, [1.0] * 4)]
    surfaces_without = [surfaces_with[0], surfaces_with[2]]
    assert surfaces_with[1].layout_ambiguous is True
    assert all(s.layout_ambiguous is False for s in surfaces_without)

    record = ShipRecord(
        name="Barge", ship_class="synthetic", length=L, beam=2 * B2, draft=T,
        block_coefficient=BARGE_VOLUME / (L * 2 * B2 * T), displacement_basis="closure",
        displacement_volume=BARGE_VOLUME, geometry_route="iges_field",
        provenance=Provenance(source="synthetic", url="https://example.invalid/barge",
                              license="test", retrieval="test", verified=False),
    )
    path = tmp / "barge.igs"
    path.write_text("synthetic barge\n")
    original = tensorcad.iges.read_iges
    builds = {}
    for label, surfaces in (("with_undecidable_patch", surfaces_with), ("without_patch", surfaces_without)):
        tensorcad.iges.read_iges = lambda _p, _s=surfaces: _FakeIges(_s)  # type: ignore[assignment]
        built = iges_field_mesh(record, path, dx=2.0, dz=0.5, samples_per_metre=1.0)
        report = built.details["patch_quarantine"]
        builds[label] = {
            "patches_total": report["patches_total"],
            "patches_used": report["patches_used"],
            "patches_quarantined": len(report["patches_quarantined"]),
            "volume_m3": measure(built.mesh).signed_volume,
            "mesh": built,
        }
    tensorcad.iges.read_iges = original  # type: ignore[assignment]

    v_with = builds["with_undecidable_patch"]["volume_m3"]
    v_without = builds["without_patch"]["volume_m3"]
    mesh = builds["with_undecidable_patch"]["mesh"].mesh
    measured = measure(mesh).signed_volume
    graded = dataclasses.replace(record, displacement_volume=measured,
                                 block_coefficient=measured / (L * 2 * B2 * T))
    manifest = FleetManifest(ships=(graded,))
    dropped, total = quarantine_census(builds["with_undecidable_patch"]["mesh"].details)
    verdict_without = evaluate_gates(graded, mesh, manifest, volume_tolerance=0.01)
    verdict_with = evaluate_gates(graded, mesh, manifest, volume_tolerance=0.01,
                                  quarantine=dropped, patch_total=total)
    return {
        "barge": {"length_m": L, "beam_m": 2 * B2, "draft_m": T,
                  "analytic_volume_m3": BARGE_VOLUME},
        "builds": {
            label: {key: value for key, value in data.items() if key != "mesh"}
            for label, data in builds.items()
        },
        "volume_difference_m3": v_with - v_without,
        "volume_difference_relative": abs(v_with - v_without) / v_without,
        "verdict_percent_error_without_census": verdict_without["volume_relative_error"] * 100.0,
        "verdict_without_census": {
            "passed": verdict_without["passed"],
            "failed_gates": verdict_without["failed_gates"],
            "geometry_completeness": verdict_without["geometry_completeness"],
        },
        "verdict_with_census": {
            "passed": verdict_with["passed"],
            "failed_gates": verdict_with["failed_gates"],
            "geometry_completeness": verdict_with["geometry_completeness"],
            "quarantined": verdict_with["patches_quarantined"],
        },
    }


def measure_fleet_evidence() -> dict[str, object]:
    diagnostics = json.loads(
        (REPO / "benchmarks/verified/baseline_fleet/geometry_diagnostics.json").read_text()
    )["hulls"]
    results = json.loads(
        (REPO / "benchmarks/verified/baseline_fleet/result.json").read_text()
    )["hulls"]
    rows = []
    for hull, diag in diagnostics.items():
        entry = results.get(hull, {})
        verdict = entry.get("verdict")
        dropped = diag.get("ambiguous_patches")
        rows.append({
            "hull": hull,
            "patches_total": diag.get("patch_count"),
            "patches_quarantined": dropped,
            "quarantined_fraction": (None if not diag.get("patch_count") or dropped is None
                                     else dropped / diag["patch_count"]),
            "verdict": verdict,
            "passed": verdict == "VERIFIED",
            "reasons": entry.get("reasons", []),
        })
    live = [row["hull"] for row in rows if row["passed"] and (row["patches_quarantined"] or 0) > 0]
    return {
        "rows": rows,
        "hulls": len(rows),
        "hulls_with_quarantine": sum(1 for row in rows if (row["patches_quarantined"] or 0) > 0),
        "hulls_passed_with_quarantine": live,
        "note": (
            "Cross-tabulated from the committed evidence, which predates G9. A hull that "
            "carries quarantined patches *and* is VERIFIED is the fail-open; an empty list "
            "means the class is latent on this fleet, not that it cannot occur."
        ),
    }


def figure(payload: dict[str, object], out: Path) -> str | None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:  # pragma: no cover - optional
        return None
    rows = payload["fleet_evidence"]["rows"]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.0, 4.2))
    names = [row["hull"] for row in rows]
    fractions = [100.0 * (row["quarantined_fraction"] or 0.0) for row in rows]
    colours = ["#c0392b" if row["passed"] else "#7f8c8d" for row in rows]
    bars = ax1.bar(names, fractions, color=colours, width=0.55)
    for bar, row in zip(bars, rows):
        ax1.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.4,
                 f"{row['patches_quarantined']}/{row['patches_total']}\n{row['verdict']}",
                 ha="center", va="bottom", fontsize=8)
    ax1.set_ylabel("quarantined patches (% of file)")
    ax1.set_title("A. Committed fleet evidence: dropped patches vs verdict\n(red = VERIFIED: the fail-open class)")
    ax1.set_ylim(0, max(fractions) * 1.45 if max(fractions) else 5)

    ab = payload["ab"]
    labels = ["without census\n(G1-G8 only)", "with census\n(+G9)"]
    values = [ab["builds"]["with_undecidable_patch"]["volume_m3"],] * 2
    verdicts = [ab["verdict_without_census"], ab["verdict_with_census"]]
    colours2 = ["#27ae60" if v["passed"] else "#c0392b" for v in verdicts]
    bars2 = ax2.bar(labels, [100.0 * ab["verdict_percent_error_without_census"]] * 2,
                    color=colours2, width=0.5)
    for bar, v in zip(bars2, verdicts):
        ax2.text(bar.get_x() + bar.get_width() / 2, bar.get_height() * 0.5,
                 f"passed={v['passed']}\n{v['geometry_completeness']}",
                 ha="center", va="center", fontsize=9, color="white")
    ax2.axhline(1.0, color="black", linestyle="--", linewidth=1)
    ax2.text(1.02, 1.05, "G7 tolerance 1%", fontsize=8, transform=ax2.get_yaxis_transform())
    ax2.set_ylabel("|V_mesh - V_ref| / V_ref  (%)")
    ax2.set_title("B. Same mesh, same volume: only the census differs\n(synthetic barge, 1 of 3 patches undecidable)")
    fig.tight_layout()
    target = out / "patch_quarantine_ab.png"
    fig.savefig(target, dpi=150)
    plt.close(fig)
    return str(target)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="runs/m22653")
    parser.add_argument("--figdir", default="docs/figs",
                        help="where to write the figure (kept with the document that cites it)")
    parser.add_argument("--tmp", default=None, help="scratch dir for the synthetic IGES")
    args = parser.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    tmp = Path(args.tmp) if args.tmp else out / "scratch"
    tmp.mkdir(parents=True, exist_ok=True)

    payload = {
        "schema": "tensorcad.fleet.patch_quarantine_ab/1.0",
        "ab": measure_ab(tmp),
        "fleet_evidence": measure_fleet_evidence(),
    }
    figdir = Path(args.figdir)
    figdir.mkdir(parents=True, exist_ok=True)
    figure_path = figure(payload, figdir)
    if figure_path:
        payload["figure"] = str(Path(figure_path).relative_to(REPO)) if Path(figure_path).is_relative_to(REPO) else figure_path
    target = out / "quarantine_gate_ab.json"
    target.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({k: v for k, v in payload.items() if k != "fleet_evidence"}, indent=2)[:2200])
    print("wrote", target, "and", figure_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())