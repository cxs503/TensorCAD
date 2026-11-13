"""Command-line quality gate: JSON mesh in, report out, non-zero exit on failure.

This is the CI entry point::

    python -m tensorcad.meshquality.cli --input mesh.json --json report.json
    python -m tensorcad.meshquality.cli --input hull.json --topology --max-genus 0

Exit status is 0 when every gated check passes, 1 when one fails, and 2 when the input
cannot be read -- so a pipeline can tell "the mesh is bad" from "the tool is broken".
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Sequence

from .io import load_spec
from .report import GateReport
from .topology import gate_non_orthogonality, gate_surface, non_orthogonality, surface_topology
from .report import DEFAULT_GATE, MetricLimit, QualityGate, measure


def _parse_limit(text: str) -> tuple[str, MetricLimit]:
    """``kind.metric>=value`` or ``kind.metric<=value`` -> (kind, MetricLimit)."""
    for operator in (">=", "<="):
        if operator in text:
            name, raw = text.split(operator, 1)
            kind, metric = name.split(".", 1)
            value = float(raw)
            if operator == ">=":
                return kind, MetricLimit(metric, minimum=value)
            return kind, MetricLimit(metric, maximum=value)
    raise argparse.ArgumentTypeError(f"expected kind.metric>=value or kind.metric<=value, got {text!r}")


def _gate_from_limits(overrides: Sequence[tuple[str, MetricLimit]]) -> QualityGate:
    limits = {kind: list(values) for kind, values in DEFAULT_GATE.limits.items()}
    for kind, limit in overrides:
        limits.setdefault(kind, [])
        limits[kind] = [existing for existing in limits[kind] if existing.metric != limit.metric]
        limits[kind].append(limit)
    return QualityGate({kind: tuple(values) for kind, values in limits.items()})


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m tensorcad.meshquality.cli",
        description="Compute element quality metrics and apply a gate (exit 1 on failure).",
    )
    parser.add_argument("--input", required=True, help="mesh spec JSON (points + elements)")
    parser.add_argument("--json", type=Path, default=None, help="write the full report here")
    parser.add_argument(
        "--limit",
        type=_parse_limit,
        action="append",
        default=[],
        metavar="KIND.METRIC>=V",
        help="override a threshold, e.g. hex.min_sicn>=0.2; repeatable",
    )
    parser.add_argument(
        "--topology",
        action="store_true",
        help="also gate the surface topology of the 'tri' elements (watertight, manifold, "
        "winding, outward); fails when the mesh has none",
    )
    parser.add_argument(
        "--max-non-orth",
        type=float,
        default=None,
        metavar="DEG",
        help="gate the worst internal-face non-orthogonality of 'hex'/'tet' cells",
    )
    parser.add_argument(
        "--max-genus",
        type=float,
        default=None,
        metavar="G",
        help="with --topology, also require the closed surface's genus to be <= G",
    )
    parser.add_argument("--no-gate", action="store_true", help="report metrics without gating")
    parser.add_argument("--quiet", action="store_true", help="print nothing on success")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        spec = load_spec(args.input)
        quality = measure(spec.points, spec.elements)
    except (OSError, ValueError, KeyError) as error:
        print(f"meshquality: cannot read {args.input}: {error}", file=sys.stderr)
        return 2

    checks = list(_gate_from_limits(args.limit).evaluate(quality).checks)
    payload: dict[str, object] = {
        "input": str(Path(args.input)),
        "metadata": spec.metadata,
        "quality": quality.to_dict(),
    }
    if args.topology:
        surface = spec.elements.get("tri", ())
        if not surface:
            print(
                "meshquality: --topology needs 'tri' surface elements (a volume-only mesh "
                "has no surface to gate)",
                file=sys.stderr,
            )
            return 2
        topology = surface_topology(spec.points, surface)
        payload["topology"] = topology.to_dict()
        checks.extend(gate_surface(topology, max_genus=args.max_genus))
    if args.max_non_orth is not None:
        kind = next((k for k in ("hex", "tet") if spec.elements.get(k)), None)
        if kind is None:
            print("meshquality: --max-non-orth needs 'hex' or 'tet' cells", file=sys.stderr)
            return 2
        faces = non_orthogonality(spec.points, spec.elements[kind], kind)
        payload["non_orthogonality"] = faces.to_dict()
        checks.extend(gate_non_orthogonality(faces, maximum=args.max_non_orth))
    report = GateReport(tuple(checks))
    payload["gate"] = report.to_dict() if not args.no_gate else None
    if args.json is not None:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if not args.quiet:
        print(json.dumps(payload["quality"], indent=2, sort_keys=True))
        for check in report.checks:
            print(f"[{'PASS' if check.passed else 'FAIL'}] {check.name}: {check.detail}")
    if not args.no_gate and not report.passed:
        print(f"meshquality: {len(report.failures)} check(s) failed", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":  # pragma: no cover - exercised through main()
    raise SystemExit(main())