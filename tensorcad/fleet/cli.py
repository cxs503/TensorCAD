"""Command-line entry point for the baseline fleet gate.

    python -m tensorcad.fleet.cli --out-dir runs/m17335_fleet

Exit codes: ``0`` when every registered ship was admitted and every negative control
was rejected, ``1`` when at least one ship was refused or failed (the run still writes
its artefacts — a refusal is evidence, not a crash), and ``2`` when the manifest itself
is malformed.

This is the same entry point the main ``tensorcad`` CLI routes to when it is called
with ``fleet``; see :func:`tensorcad.cli.main`.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Sequence

from tensorcad.fleet.gates import VOLUME_TOLERANCE
from tensorcad.fleet.manifest import ManifestError
from tensorcad.fleet.run import run_fleet


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m tensorcad.fleet.cli",
        description="Build, gate and record the six-hull baseline fleet "
                    "(GAP-CAD-T1-BASELINE-FLEET).",
    )
    parser.add_argument("--out-dir", default="runs/m17335_fleet",
                        help="directory for evidence.json, ref_vs_got.csv, manifests/, figures")
    parser.add_argument("--manifest", default=None, help="fleet manifest JSON (default: packaged one)")
    parser.add_argument("--nx", type=int, default=481, help="stations for the geometry routes")
    parser.add_argument("--nz", type=int, default=97, help="waterline levels for the geometry routes")
    parser.add_argument("--volume-tolerance", type=float, default=VOLUME_TOLERANCE,
                        help="relative volume agreement required (mission: 0.01)")
    parser.add_argument("--no-plot", action="store_true", help="skip the summary figure")
    parser.add_argument("--allow-failures", action="store_true",
                        help="exit 0 even when a ship is refused (artefacts are still written)")
    return parser


def _plot(evidence: dict, path: Path) -> str | None:
    """Summary figure: volume agreement per ship and the negative-control verdicts."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:  # pragma: no cover - plotting is optional
        return f"matplotlib unavailable: {exc}"
    ships = [block for block in evidence["ships"] if block["topology"] is not None]
    refused = [block["ship"] for block in evidence["ships"] if block["topology"] is None]
    if not ships:
        return "no admitted geometry to plot"
    figure, (left, right) = plt.subplots(1, 2, figsize=(11.5, 4.4))
    names = [block["ship"] for block in ships]
    errors = [100.0 * block["volume_relative_error"] for block in ships]
    colours = ["#2b6cb0" if not block["reference_is_external"] else "#2f855a" for block in ships]
    bars = left.bar(names, errors, color=colours)
    left.axhline(1.0, color="#c53030", linestyle="--", linewidth=1.2, label="1 % mission tolerance")
    left.set_yscale("log")
    left.set_ylabel("|V_mesh - V_ref| / V_ref   [%]")
    left.set_title("Hull volume agreement\n(green: literature reference, blue: closure reference)")
    left.legend(loc="lower right", fontsize=8)
    for bar, value in zip(bars, errors):
        left.annotate(f"{value:.3g}%", (bar.get_x() + bar.get_width() / 2, value),
                      ha="center", va="bottom", fontsize=8)
    if refused:
        left.text(0.02, 0.02, "refused: " + ", ".join(refused), transform=left.transAxes, fontsize=8,
                  color="#c53030")
    controls = evidence["negative_controls"]
    labels = [control["control"].replace("_", "\n") for control in controls]
    verdicts = [0 if control["rejected"] else 1 for control in controls]
    right.bar(labels, verdicts, color=["#2f855a" if v == 0 else "#c53030" for v in verdicts])
    right.set_ylim(0, 1.5)
    right.set_yticks([0, 1])
    right.set_yticklabels(["rejected (correct)", "ACCEPTED (gate hole)"])
    right.set_title(f"Negative controls ({sum(1 for v in verdicts if v == 0)}/{len(verdicts)} rejected)")
    for axis in (left, right):
        axis.tick_params(axis="x", labelsize=8)
    figure.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=160)
    plt.close(figure)
    return None


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(list(argv) if argv is not None else None)
    try:
        evidence = run_fleet(
            args.out_dir,
            manifest_path=args.manifest,
            nx=args.nx,
            nz=args.nz,
            volume_tolerance=args.volume_tolerance,
        )
    except ManifestError as exc:
        print(f"manifest error: {exc}", file=sys.stderr)
        return 2
    out = Path(args.out_dir)
    plot_note = None
    if not args.no_plot:
        plot_note = _plot(evidence, out / "fleet_summary.png")
        if plot_note:
            print(f"figure skipped: {plot_note}", file=sys.stderr)
    summary = evidence["summary"]
    print(json.dumps(summary, indent=2))
    print(f"artefacts: {out / 'evidence.json'}, {out / 'ref_vs_got.csv'}, {out / 'manifests'}")
    if args.allow_failures or summary["controls_all_rejected"] and not summary["ships_refused_or_failed"]:
        return 0
    return 1


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())