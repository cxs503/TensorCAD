"""``tensorcad-neural``: fit, edit, inspect and extract neural geometry objects.

The CLI is the agent-facing surface for the neural geometry pipeline. Every
command prints measured numbers rather than a success banner, emits machine
readable JSON with ``--json``, and returns a non-zero status when an object fails
its own gates — a rolled-back edit exits 1 with the failed check names, not 0.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

from .neural_geometry import (
    NeuralGeometryFormatError,
    ObjectStore,
    TransactionError,
    dumps_canonical,
    read_manifest,
    sha256_file,
)

EXIT_OK = 0
EXIT_FAILED = 1
EXIT_ERROR = 2


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="tensorcad-neural",
        description="Manage geometry as a learned continuous field instead of a mesh file.",
    )
    parser.add_argument("--json", action="store_true", help="emit machine-readable JSON")
    subparsers = parser.add_subparsers(dest="command", required=True)

    fit = subparsers.add_parser("fit", help="fit a new object and publish revision 1")
    fit.add_argument("--output", required=True, help="bundle directory to create, e.g. results/sphere.ngo")
    fit.add_argument("--radius", type=float, default=0.5, help="sphere radius in metres")
    fit.add_argument("--steps", type=int, default=1800)
    fit.add_argument("--seed", type=int, default=503)
    fit.add_argument("--rays", type=int, default=2048, help="validation ray count")
    fit.add_argument("--json", action="store_true", default=argparse.SUPPRESS, help=argparse.SUPPRESS)

    edit = subparsers.add_parser("edit", help="apply a local dent as a transaction")
    edit.add_argument("--object", required=True, help="existing bundle directory")
    edit.add_argument("--depth", type=float, required=True, help="dent depth in metres")
    edit.add_argument("--sigma", type=float, default=0.22, help="dent angular width in radians")
    edit.add_argument("--steps", type=int, default=1000)
    edit.add_argument("--seed", type=int, default=503)
    edit.add_argument("--resolution", type=int, default=0, help="also extract a mesh at this resolution")
    edit.add_argument("--rays", type=int, default=4096, help="validation ray count")
    edit.add_argument("--json", action="store_true", default=argparse.SUPPRESS, help=argparse.SUPPRESS)

    extract = subparsers.add_parser("extract", help="derive a mesh (OBJ/STL) from the active revision")
    extract.add_argument("--object", required=True)
    extract.add_argument("--resolution", type=int, default=96, help="cells per axis")
    extract.add_argument("--format", default="obj,stl", help="comma separated: obj,stl")
    extract.add_argument("--output", default=None, help="directory for derived files")
    extract.add_argument("--json", action="store_true", default=argparse.SUPPRESS, help=argparse.SUPPRESS)

    info = subparsers.add_parser("info", help="print the manifest and the transaction history")
    info.add_argument("--object", required=True)
    info.add_argument("--json", action="store_true", default=argparse.SUPPRESS, help=argparse.SUPPRESS)

    verify = subparsers.add_parser("verify", help="re-open the object and check its integrity")
    verify.add_argument("--object", required=True)
    verify.add_argument("--json", action="store_true", default=argparse.SUPPRESS, help=argparse.SUPPRESS)
    return parser


def args_rest(args: argparse.Namespace) -> list[str]:
    """The argv of this command, for the provenance record (``main`` strips the program)."""
    return list(getattr(args, "_argv", []))


def _emit(payload: dict[str, Any], as_json: bool, lines: list[str] | None = None) -> None:
    if as_json:
        print(dumps_canonical(payload), end="")
    elif lines:
        print("\n".join(lines))
    else:
        print(json.dumps(payload, indent=2, sort_keys=True, default=str))


def _format_metrics(metrics: dict[str, Any], keys: list[tuple[str, str]]) -> list[str]:
    lines = []
    for key, label in keys:
        value = metrics.get(key)
        if isinstance(value, float):
            lines.append(f"  {label:<34s} {value:.6g}")
        else:
            lines.append(f"  {label:<34s} {value}")
    return lines


FIT_KEYS = [
    ("field_mae_m", "field MAE [m]"),
    ("field_rmse_m", "field RMSE [m]"),
    ("boundary_error_mean_m", "surface error mean [m]"),
    ("boundary_error_p95_m", "surface error P95 [m]"),
    ("boundary_error_max_m", "surface error max [m]"),
    ("ray_coverage", "ray coverage"),
    ("training_seconds", "training [s]"),
]

EDIT_KEYS = [
    ("ray_coverage", "ray coverage"),
    ("dent_peak_target_m", "dent peak target [m]"),
    ("dent_peak_measured_m", "dent peak measured [m]"),
    ("dent_peak_error_m", "dent peak error [m]"),
    ("dent_region_error_mean_m", "dent region error mean [m]"),
    ("protected_mean_displacement_m", "protected mean displacement [m]"),
    ("protected_p95_abs_displacement_m", "protected P95 displacement [m]"),
    ("protected_p99_abs_displacement_m", "protected P99 displacement [m]"),
    ("protected_max_abs_displacement_m", "protected max displacement [m]"),
]


def _cmd_fit(args: argparse.Namespace) -> int:
    from .neural.objects import NeuralObjectError, fit_sphere

    try:
        result = fit_sphere(
            args.output,
            radius=args.radius,
            steps=args.steps,
            seed=args.seed,
            rays=args.rays,
            command=["tensorcad-neural", "fit", *args_rest(args)],
        )
    except (NeuralObjectError, ValueError, TransactionError) as exc:
        _emit({"status": "error", "stage": "fit", "error": str(exc)}, args.json)
        return EXIT_ERROR
    payload = {
        "status": "committed",
        "object": str(result.store),
        "revision": result.revision,
        "weights_sha256": result.weights_sha256,
        "metrics": result.metrics,
        "checks": [check.to_dict() for check in result.checks],
    }
    _emit(
        payload,
        args.json,
        [
            f"fit  {result.store}  revision {result.revision}",
            f"  weights sha256                     {result.weights_sha256}",
            *_format_metrics(result.metrics, FIT_KEYS),
        ],
    )
    return EXIT_OK


def _cmd_edit(args: argparse.Namespace) -> int:
    from .neural.objects import NeuralObjectError, apply_dent_edit

    try:
        result = apply_dent_edit(
            args.object,
            depth=args.depth,
            sigma=args.sigma,
            steps=args.steps,
            seed=args.seed,
            extract_resolution=args.resolution or None,
            rays=args.rays,
            command=["tensorcad-neural", "edit", *args_rest(args)],
        )
    except (NeuralObjectError, NeuralGeometryFormatError, TransactionError, ValueError) as exc:
        _emit({"status": "error", "stage": "edit", "error": str(exc)}, args.json)
        return EXIT_ERROR
    payload = {
        "status": result.status,
        "parent_revision": result.parent_revision,
        "active_revision": result.active_revision,
        "weights_sha256": result.weights_sha256,
        "metrics": result.metrics,
        "checks": [check.to_dict() for check in result.checks],
        "failed_checks": list(result.failed_checks),
        "reason": result.reason,
        "derived": [record.get("path") for record in result.derived],
    }
    lines = [
        f"edit {args.object}: {result.status}",
        f"  parent revision                    {result.parent_revision}",
        f"  active revision                    {result.active_revision}",
        *_format_metrics(result.metrics, EDIT_KEYS),
    ]
    if result.reason:
        lines.append(f"  reason                             {result.reason}")
    _emit(payload, args.json, lines)
    return EXIT_OK if result.committed else EXIT_FAILED


def _cmd_extract(args: argparse.Namespace) -> int:
    from .neural.objects import NeuralObjectError, extract_object

    formats = tuple(part.strip() for part in args.format.split(",") if part.strip())
    try:
        _mesh, report, records = extract_object(
            args.object,
            resolution=args.resolution,
            formats=formats,
            output=args.output,
        )
    except (NeuralObjectError, NeuralGeometryFormatError, ValueError) as exc:
        _emit({"status": "error", "stage": "extract", "error": str(exc)}, args.json)
        return EXIT_ERROR
    payload = {
        "status": "derived",
        "object": str(Path(args.object)),
        "active_revision": ObjectStore(args.object).active_revision(),
        "resolution": args.resolution,
        "report": report.to_dict(),
        "artifacts": [
            {"format": record["format"], "path": record["path"], "sha256": record["sha256"]}
            for record in records
        ],
        "note": "derived mesh; the neural field remains the canonical geometry",
    }
    _emit(
        payload,
        args.json,
        [
            f"extract {args.object} at {args.resolution}^3 cells",
            f"  vertices                           {report.vertices}",
            f"  triangles                          {report.triangles}",
            f"  watertight                         {report.watertight}",
            f"  edge manifold                      {report.edge_manifold}",
            f"  euler characteristic               {report.euler_characteristic}",
            f"  signed volume [m^3]                {report.signed_volume:.6g}",
            f"  surface area [m^2]                 {report.surface_area:.6g}",
            *[
                f"  wrote                              {record['path']} ({record['format']})"
                for record in records
            ],
        ],
    )
    return EXIT_OK if report.watertight and report.edge_manifold else EXIT_FAILED


def _cmd_info(args: argparse.Namespace) -> int:
    store = ObjectStore(args.object)
    try:
        manifest = store.manifest()
    except (NeuralGeometryFormatError, OSError) as exc:
        _emit({"status": "error", "stage": "info", "error": str(exc)}, args.json)
        return EXIT_ERROR
    history = store.history()
    payload = {
        "manifest": manifest.summary(),
        "history": history,
    }
    lines = [
        f"object        {manifest.object_id}",
        f"bundle        {store.root}",
        f"revision      {manifest.revision} (parent {manifest.lineage.parent_revision})",
        f"operation     {manifest.lineage.operation}",
        f"format        {manifest.format_version}",
        f"architecture  {manifest.architecture_id}",
        f"field         {manifest.field.get('quantity')} [{manifest.units}], "
        f"tau={manifest.field.get('truncation_distance')}",
        f"validation    {manifest.validation_status}",
        f"weights       {manifest.weights.sha256}",
        f"derived       {', '.join(str(item.get('kind')) for item in manifest.derived) or 'none'}",
        "",
        "history:",
        *[
            f"  rev {entry.get('revision')}  {entry.get('status'):<12s} {entry.get('operation', ''):<10s}"
            f"  {entry.get('recorded_utc', '')}"
            + (f"  failed: {', '.join(entry.get('failed_checks', []))}" if entry.get("failed_checks") else "")
            for entry in history
        ],
    ]
    _emit(payload, args.json, lines)
    return EXIT_OK


def _cmd_verify(args: argparse.Namespace) -> int:
    """Re-open the object and check everything that does not need a training loop.

    The integrity checks (manifest, revision, weights hash, history, derived
    artifacts) are standard-library only, so they must keep working when the
    optional ``neural`` extra is absent. Only the live field query needs torch,
    and it is reported as skipped instead of taking the whole command down.
    """
    store = ObjectStore(args.object)
    checks: list[dict[str, Any]] = []

    def record(name: str, passed: bool, detail: str, severity: str = "error") -> None:
        checks.append({"name": name, "passed": passed, "detail": detail, "severity": severity})

    manifest = None
    try:
        manifest = read_manifest(store.root, verify_payload=False)
        record("manifest_parses", True, f"format {manifest.format_version}")
        record("revision_positive", manifest.revision >= 1, f"revision {manifest.revision}")
        actual = sha256_file(store.root / manifest.weights.path)
        record("weights_hash_matches", actual == manifest.weights.sha256, actual[:16])
        history = store.history()
        committed = [entry for entry in history if entry.get("status") == "committed"]
        record(
            "history_matches_revision",
            bool(committed) and int(committed[-1].get("revision", -1)) == manifest.revision,
            f"{len(history)} attempts, last committed revision "
            f"{committed[-1].get('revision') if committed else None}",
        )
    except Exception as exc:
        checks.append(
            {"name": "verify_ran", "passed": False,
             "detail": f"{type(exc).__name__}: {exc}", "severity": "error"}
        )

    if manifest is not None:
        for artifact in manifest.derived:
            path = artifact.get("path", "?")
            target = store.root / str(path)
            if not target.is_file():
                record(
                    f"derived_present:{path}",
                    False,
                    "recorded derived artifact is missing from the bundle",
                    severity="warning",
                )
                continue
            digest = sha256_file(target)
            record(
                f"derived_hash:{path}",
                digest == artifact.get("sha256"),
                f"{digest[:16]} vs recorded {str(artifact.get('sha256'))[:16]}",
            )
            expected_size = artifact.get("bytes")
            record(
                f"derived_size:{path}",
                expected_size is None or target.stat().st_size == expected_size,
                f"{target.stat().st_size} bytes",
            )

        try:
            from .neural.objects import build_field

            field = build_field(manifest, store.root / manifest.weights.path)
            import numpy as np

            points = np.array(
                [[0.0, 0.0, 0.0], [0.05, 0.0, 0.0], [0.9, 0.0, 0.0]], dtype=np.float64
            )
            values = field(points)
            record(
                "field_query",
                bool(np.all(np.isfinite(values))),
                "values at centre/inside/outside: " + ", ".join(f"{value:+.4f}" for value in values),
            )
        except ImportError as exc:
            record(
                "field_query",
                True,
                f"skipped: the neural extra is not installed ({exc})",
                severity="skipped",
            )
        except Exception as exc:
            record("field_query", False, f"{type(exc).__name__}: {exc}")

    failed = [check for check in checks if not check["passed"] and check["severity"] == "error"]
    warned = [check for check in checks if not check["passed"] and check["severity"] != "error"]
    passed = not failed
    payload = {
        "status": "verified" if passed else "failed",
        "checks": checks,
        "warnings": [check["name"] for check in warned],
    }
    lines = [
        *[
            f"  [{'PASS' if check['passed'] else ('WARN' if check['severity'] != 'error' else 'FAIL')}] "
            f"{check['name']}: {check['detail']}"
            for check in checks
        ],
        f"verify {args.object}: {'passed' if passed else 'FAILED'}"
        + (f" ({len(warned)} warning(s))" if warned else ""),
    ]
    _emit(payload, args.json, lines)
    return EXIT_OK if passed else EXIT_FAILED


COMMANDS = {
    "fit": _cmd_fit,
    "edit": _cmd_edit,
    "extract": _cmd_extract,
    "info": _cmd_info,
    "verify": _cmd_verify,
}


def main(argv: list[str] | None = None) -> int:
    arguments = list(sys.argv[1:] if argv is None else argv)
    args = build_parser().parse_args(arguments)
    args._argv = arguments
    handler = COMMANDS[args.command]
    return int(handler(args))


if __name__ == "__main__":  # pragma: no cover - module entry point
    raise SystemExit(main())
