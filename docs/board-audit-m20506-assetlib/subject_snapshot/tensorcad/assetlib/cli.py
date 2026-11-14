"""``python -m tensorcad.assetlib`` — the library's operator surface.

Sub-commands (all print JSON so a gate can read them):

``build-registry``  scan the store, re-hash every file, join the fleet specs and the
                    benchmark verdict, write the versioned registry + a history line;
``verify``          re-hash the registered bytes and report missing/mismatched ones;
``resolve``         resolve a ship name or asset URI to a verified path (exit 3 when
                    the byte does not match its fingerprint);
``growth``          differentiate the current registry against a previous snapshot;
``repair-specs``    re-point dead pointers in a fleet-spec file (``--dry-run`` first).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .registry import build_registry, latest_history, load_registry, write_registry
from .specs import repair_spec_file, scan_spec_pointers
from .store import AssetError, AssetStore

REPO = Path(__file__).resolve().parents[2]
DEFAULT_SPECS = REPO / "benchmarks" / "verified" / "baseline_fleet" / "fleet_specs.json"
DEFAULT_BENCHMARK = REPO / "benchmarks" / "verified" / "baseline_fleet" / "result.json"
DEFAULT_REGISTRY = REPO / "assets" / "geometry_registry.json"
DEFAULT_HISTORY = REPO / "assets" / "geometry_registry_history.jsonl"


def _url_maps() -> tuple[dict[str, str], dict[str, str]]:
    """Publisher URLs already recorded by the fetcher, keyed by page key / archive."""
    try:
        sys.path.insert(0, str(REPO / "tools" / "baseline_fleet"))
        from fetch_sources import PAGES, TARGETS  # type: ignore
    except Exception:
        return {}, {}
    return dict(PAGES), {spec[0]: spec[1] for spec in TARGETS.values()}


def _cmd_build(args: argparse.Namespace) -> int:
    store = AssetStore(args.store)
    registry = build_registry(
        store=store,
        specs_path=args.specs,
        benchmark_path=args.benchmark,
    )
    if not args.no_write:
        write_registry(registry, args.out, history=None if args.no_history else args.history)
    print(json.dumps({"snapshot": registry.snapshot, "written": None if args.no_write else str(args.out)}, indent=1))
    return 0 if not registry.snapshot.get("problems") else 1


def _cmd_verify(args: argparse.Namespace) -> int:
    registry = load_registry(args.registry)
    report = registry.verify(store=AssetStore(args.store or registry.store_root or None))
    print(json.dumps(report, indent=1))
    return 0 if report["ok"] else 1


def _cmd_resolve(args: argparse.Namespace) -> int:
    registry = load_registry(args.registry)
    try:
        path = registry.resolve(args.name, store=AssetStore(args.store or registry.store_root or None))
    except (AssetError, KeyError) as exc:
        reason = getattr(exc, "reason", "not_registered")
        print(json.dumps({"name": args.name, "resolved": False, "reason": reason, "detail": str(exc)}, indent=1))
        return 3
    print(json.dumps({"name": args.name, "resolved": True, "path": path}, indent=1))
    return 0


def _cmd_growth(args: argparse.Namespace) -> int:
    registry = load_registry(args.registry)
    previous = latest_history(args.history) if args.history else None
    print(json.dumps(registry.growth(previous), indent=1))
    return 0


def _cmd_scan(args: argparse.Namespace) -> int:
    registry = load_registry(args.registry) if Path(args.registry).is_file() else None
    report = scan_spec_pointers(
        args.specs, store=AssetStore(args.store or (registry.store_root if registry else "") or None), registry=registry
    )
    print(json.dumps(report.as_dict(), indent=1))
    return 0 if not report.pointers or report.dead or report.withheld else 0


def _cmd_repair(args: argparse.Namespace) -> int:
    registry = load_registry(args.registry) if Path(args.registry).is_file() else None
    url_map, archive_map = _url_maps()
    result = repair_spec_file(
        args.specs,
        out_path=args.out,
        store=AssetStore(args.store or (registry.store_root if registry else "") or None),
        registry=registry,
        repo_root=REPO,
        url_map=url_map,
        archive_url_map=archive_map,
        dry_run=args.dry_run,
    )
    print(json.dumps(result, indent=1))
    return 0 if not result["unresolved"] else 1


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="python -m tensorcad.assetlib", description=__doc__)
    parser.add_argument("--store", default=None, help="asset store root (default $TENSORCAD_ASSET_ROOT)")
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("build-registry", help="scan the store and write the registry")
    p.add_argument("--specs", default=str(DEFAULT_SPECS))
    p.add_argument("--benchmark", default=str(DEFAULT_BENCHMARK))
    p.add_argument("--out", default=str(DEFAULT_REGISTRY))
    p.add_argument("--history", default=str(DEFAULT_HISTORY))
    p.add_argument("--no-history", action="store_true")
    p.add_argument("--no-write", action="store_true")
    p.set_defaults(func=_cmd_build)

    p = sub.add_parser("verify", help="re-hash the registered bytes")
    p.add_argument("--registry", default=str(DEFAULT_REGISTRY))
    p.set_defaults(func=_cmd_verify)

    p = sub.add_parser("resolve", help="resolve a ship name or asset URI to a verified path")
    p.add_argument("name")
    p.add_argument("--registry", default=str(DEFAULT_REGISTRY))
    p.set_defaults(func=_cmd_resolve)

    p = sub.add_parser("growth", help="differentiate the registry against the last snapshot")
    p.add_argument("--registry", default=str(DEFAULT_REGISTRY))
    p.add_argument("--history", default=str(DEFAULT_HISTORY))
    p.set_defaults(func=_cmd_growth)

    p = sub.add_parser("scan-specs", help="classify the pointers in a fleet-spec file")
    p.add_argument("--specs", default=str(DEFAULT_SPECS))
    p.add_argument("--registry", default=str(DEFAULT_REGISTRY))
    p.set_defaults(func=_cmd_scan)

    p = sub.add_parser("repair-specs", help="re-point dead pointers in a fleet-spec file")
    p.add_argument("--specs", default=str(DEFAULT_SPECS))
    p.add_argument("--registry", default=str(DEFAULT_REGISTRY))
    p.add_argument("--out", default=None)
    p.add_argument("--dry-run", action="store_true")
    p.set_defaults(func=_cmd_repair)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())