"""Replace placeholder / duplicated STL files in a ship geometry library, with a ledger.

The defect this closes is the one the owner reported on 2026-11-13: the ship-platform
geometry library ships eight STL files of which six are 304-triangle placeholders and
two "different" ships are byte-identical. A gate that only says "rejected" leaves the
library broken — somebody still has to put a real hull in the file. This module is that
step, and it is deliberately narrow:

* it never fabricates a hull for a ship it has no record for; the file is reported
  ``withheld`` with the reason, and the ledger counts it, so the library cannot quietly
  look complete;
* it never substitutes a different ship's geometry; the replacement comes from the
  fleet manifest record that names the same ship, through the same geometry routes and
  the same gates the fleet run uses;
* every replacement is written to disk, read back, re-measured and gated, and the
  pre/post sha256 pair goes into the ledger, so "replaced" is an auditable claim about
  bytes rather than an intention.

Workflow (also exposed as ``python -m tensorcad.fleet.replace``)::

    python -m tensorcad.fleet.replace --source-dir <library> --out-dir <delivery>

The delivery directory keeps the source file names, so re-pointing the library is a
copy, and ``ledger.json`` says which file came from which record and why.
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Sequence

from tensorcad.fleet.build import (
    RouteRefusal,
    build_ship_geometry,
    read_stl,
    write_binary_stl,
)
from tensorcad.fleet.gates import VOLUME_TOLERANCE, evaluate_gates, measure, placeholder_reasons
from tensorcad.fleet.manifest import FleetManifest, ShipRecord, load_manifest, sha256_file
from tensorcad.fleet.sources import load_registry

SCHEMA = "tensorcad.fleet.replace.v1"
ANCHOR = "GAP-CAD-T1-BASELINE-FLEET"

#: Classification of one file in the source library.
PLACEHOLDER = "placeholder"
DUPLICATE = "duplicate"
VALID = "valid"


@dataclass
class FileReport:
    """What one file in the source library is, measured rather than assumed."""

    path: Path
    sha256: str
    bytes: int
    classification: str
    reasons: list[str] = field(default_factory=list)
    triangle_count: int | None = None
    signed_volume: float | None = None
    boundary_edge_count: int | None = None
    nonmanifold_edge_count: int | None = None
    error: str | None = None

    @property
    def name(self) -> str:
        return self.path.name

    @property
    def stem(self) -> str:
        return self.path.stem

    def as_dict(self) -> dict[str, Any]:
        return {
            "file": self.name,
            "stem": self.stem,
            "sha256": self.sha256,
            "bytes": self.bytes,
            "classification": self.classification,
            "reasons": list(self.reasons),
            "triangle_count": self.triangle_count,
            "signed_volume": self.signed_volume,
            "boundary_edge_count": self.boundary_edge_count,
            "nonmanifold_edge_count": self.nonmanifold_edge_count,
            "error": self.error,
        }


def _stl_files(directory: str | Path) -> list[Path]:
    return sorted(path for path in Path(directory).glob("*.stl") if path.is_file())


def audit_directory(directory: str | Path, *, max_placeholder_triangles: int = 304) -> list[FileReport]:
    """Classify every ``*.stl`` in ``directory`` as placeholder, duplicate or valid.

    ``duplicate`` is a *fleet-level* defect: it can only be seen by comparing files, so
    the report is produced for the whole directory and not per file in isolation. A file
    that cannot be parsed is reported with its error and classified ``placeholder``
    (it certainly is not usable geometry) instead of disappearing from the count.
    """
    reports: list[FileReport] = []
    for path in _stl_files(directory):
        sha = sha256_file(path)
        size = path.stat().st_size
        try:
            mesh = read_stl(path)
        except Exception as exc:  # noqa: BLE001 - a broken file is a finding, not a crash
            reports.append(
                FileReport(
                    path=path,
                    sha256=sha,
                    bytes=size,
                    classification=PLACEHOLDER,
                    reasons=[f"unreadable: {exc}"],
                    error=str(exc),
                )
            )
            continue
        facts = measure(mesh)
        reasons = placeholder_reasons(
            facts, max_placeholder_triangles=max_placeholder_triangles, mesh=mesh
        )
        reports.append(
            FileReport(
                path=path,
                sha256=sha,
                bytes=size,
                classification=PLACEHOLDER if reasons else VALID,
                reasons=reasons,
                triangle_count=facts.triangle_count,
                signed_volume=facts.signed_volume,
                boundary_edge_count=facts.boundary_edge_count,
                nonmanifold_edge_count=facts.nonmanifold_edge_count,
            )
        )
    by_sha: dict[str, list[FileReport]] = {}
    for report in reports:
        by_sha.setdefault(report.sha256, []).append(report)
    for sha, group in by_sha.items():
        if len(group) < 2:
            continue
        names = sorted(report.name for report in group)
        for report in group:
            report.reasons.append(
                f"byte_identical_geometry: {sha[:16]}... shared by {names}"
            )
            if report.classification == VALID:
                report.classification = DUPLICATE
    return reports


def plan_replacements(
    reports: Sequence[FileReport], manifest: FleetManifest
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Pair each defective file with the manifest record that can replace it.

    Returns ``(actionable, withheld)``. A file whose stem does not name a fleet record is
    withheld: fabricating geometry for an unregistered ship would be exactly the kind of
    silent invention the fleet anchor exists to prevent.
    """
    by_name = {record.name.upper(): record for record in manifest.ships}
    actionable: list[dict[str, Any]] = []
    withheld: list[dict[str, Any]] = []
    for report in reports:
        if report.classification == VALID:
            continue
        record = by_name.get(report.stem.upper())
        if record is None:
            withheld.append(
                {
                    **report.as_dict(),
                    "status": "WITHHELD",
                    "action": "withheld",
                    "reason": f"no fleet-manifest record named {report.stem!r}: geometry is not "
                              f"fabricated for a ship with no published particulars on file",
                }
            )
            continue
        actionable.append({**report.as_dict(), "action": "replace", "record": record.as_dict(),
                           "ship": record.name})
    return actionable, withheld


def _deliver_one(
    entry: dict[str, Any],
    record: ShipRecord,
    manifest: FleetManifest,
    out_dir: Path,
    *,
    nx: int,
    nz: int,
    volume_tolerance: float,
    registry,
) -> dict[str, Any]:
    """Build, gate, write and re-ingest one replacement; never raises on refusal."""
    block = dict(entry)
    try:
        built = build_ship_geometry(record, nx=nx, nz=nz)
    except RouteRefusal as exc:
        block.update({"action": "withheld", "reason": f"geometry route refused: {exc}",
                      "status": "REFUSED"})
        return block
    verdict = evaluate_gates(record, built.mesh, manifest, volume_tolerance=volume_tolerance,
                             registry=registry)
    target = out_dir / entry["file"]
    written_sha = write_binary_stl(built.mesh, target)
    stored_sha = sha256_file(target)
    reread = read_stl(target)
    facts = measure(reread)
    reference = record.reference_volume()
    volume_rel = abs(facts.signed_volume - reference) / reference
    reingest_ok = (
        stored_sha == written_sha
        and facts.boundary_edge_count == 0
        and facts.nonmanifold_edge_count == 0
        and facts.inconsistently_oriented_edge_count == 0
        and abs(facts.euler_characteristic - 2) == 0
        and abs(facts.signed_volume - verdict["topology"]["signed_volume"])
        / max(abs(verdict["topology"]["signed_volume"]), 1e-12)
        < 1e-6
    )
    block.update(
        {
            "status": "REPLACED" if verdict["passed"] and reingest_ok else "FAILED",
            "action": "replaced" if verdict["passed"] and reingest_ok else "withheld",
            "replacement_file": target.name,
            "replacement_sha256": stored_sha,
            "replacement_bytes": target.stat().st_size,
            "sha256_changed": stored_sha != entry["sha256"],
            "triangles": facts.triangle_count,
            "boundary_edge_count": facts.boundary_edge_count,
            "nonmanifold_edge_count": facts.nonmanifold_edge_count,
            "euler_characteristic": facts.euler_characteristic,
            "volume_mesh": facts.signed_volume,
            "volume_reference": reference,
            "volume_relative_error": volume_rel,
            "volume_tolerance": volume_tolerance,
            "reference_basis": record.reference_basis(),
            "gates_passed": [gate["name"] for gate in verdict["gates"] if gate["passed"]],
            "gates_failed": verdict["failed_gates"],
            "reingest_ok": reingest_ok,
            "verification_conclusion": (
                f"REPLACED {entry['file']}: sha256 {entry['sha256'][:16]}... -> {stored_sha[:16]}..., "
                f"{facts.triangle_count} triangles, watertight (boundary edges "
                f"{facts.boundary_edge_count}), manifold+oriented (non-manifold "
                f"{facts.nonmanifold_edge_count}), Euler {facts.euler_characteristic}; V_mesh "
                f"{facts.signed_volume:.6g} m^3 vs V_ref {reference:.6g} m^3 = "
                f"{volume_rel * 100:.4f}% (basis {record.reference_basis()}); gates "
                f"{len(verdict['gates']) - len(verdict['failed_gates'])}/{len(verdict['gates'])} passed; "
                f"re-ingest of the delivered file {'agrees' if reingest_ok else 'DISAGREES'}"
            ),
        }
    )
    if not (verdict["passed"] and reingest_ok):
        block["reason"] = (
            f"replacement did not pass: failed gates {verdict['failed_gates']}, "
            f"re-ingest {'ok' if reingest_ok else 'failed'}"
        )
    return block


def replace_library(
    source_dir: str | Path,
    out_dir: str | Path,
    *,
    manifest_path: str | Path | None = None,
    nx: int = 481,
    nz: int = 97,
    volume_tolerance: float = VOLUME_TOLERANCE,
    copy_unreplaced: bool = True,
    write: bool = True,
) -> dict[str, Any]:
    """Audit ``source_dir``, deliver replaced files to ``out_dir`` and return the ledger.

    Valid files are copied through by default, so ``out_dir`` is a drop-in replacement
    directory for the library rather than a partial patch.
    """
    source = Path(source_dir)
    delivery = Path(out_dir)
    manifest = load_manifest(manifest_path)
    registry = load_registry()
    reports = audit_directory(source)
    actionable, withheld = plan_replacements(reports, manifest)
    if write:
        delivery.mkdir(parents=True, exist_ok=True)
    delivered: list[dict[str, Any]] = []
    for entry in actionable:
        record = manifest.by_name(entry["ship"])
        delivered.append(
            _deliver_one(entry, record, manifest, delivery, nx=nx, nz=nz,
                         volume_tolerance=volume_tolerance, registry=registry)
        )
    copied: list[str] = []
    if copy_unreplaced and write:
        replaced = {block["file"] for block in delivered if block.get("replacement_file")}
        for report in reports:
            # everything that was not replaced is passed through unchanged: a drop-in
            # directory has to keep the library's file set, and a file the pipeline cannot
            # fix is a *reported* defect, not a deleted one
            if report.name in replaced:
                continue
            target = delivery / report.name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(report.path, target)
            copied.append(report.name)
    summary = {
        "files_scanned": len(reports),
        "placeholders": [report.name for report in reports if report.classification == PLACEHOLDER],
        "duplicates": [report.name for report in reports if report.classification == DUPLICATE],
        "valid": [report.name for report in reports if report.classification == VALID],
        "defective": len(actionable) + len(withheld),
        "replaced": [block["file"] for block in delivered if block["status"] == "REPLACED"],
        "withheld": len(withheld) + sum(1 for block in delivered if block["status"] != "REPLACED"),
        "passed_through": copied,
        "all_replacements_gated": all(
            block["status"] == "REPLACED" for block in delivered
        ) if delivered else True,
        "worst_replacement_volume_error": max(
            (block["volume_relative_error"] for block in delivered
             if block.get("volume_relative_error") is not None),
            default=None,
        ),
    }
    ledger = {
        "schema": SCHEMA,
        "anchor": ANCHOR,
        "settings": {
            "source_dir": str(source),
            "out_dir": str(delivery),
            "nx": nx,
            "nz": nz,
            "volume_tolerance": volume_tolerance,
            "manifest_sha256": manifest.sha256,
            "source_registry_sha256": registry.sha256,
        },
        "files": [block for block in delivered] + withheld,
        "summary": summary,
    }
    if write:
        (delivery / "ledger.json").write_text(json.dumps(ledger, indent=2, default=str) + "\n")
    return ledger


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m tensorcad.fleet.replace",
        description="Audit a ship geometry library and deliver gated replacements for "
                    "placeholder / duplicated STL files (GAP-CAD-T1-BASELINE-FLEET).",
    )
    parser.add_argument("--source-dir", required=True, help="the library directory to audit")
    parser.add_argument("--out-dir", required=True, help="delivery directory (drop-in replacement)")
    parser.add_argument("--manifest", default=None, help="fleet manifest JSON (default: packaged one)")
    parser.add_argument("--nx", type=int, default=481)
    parser.add_argument("--nz", type=int, default=97)
    parser.add_argument("--volume-tolerance", type=float, default=VOLUME_TOLERANCE)
    parser.add_argument("--no-copy-unreplaced", action="store_true",
                        help="do not pass unreplaced files through to the delivery directory")
    return parser


def main(argv: Sequence[str] | None = None) -> int:  # pragma: no cover - thin wrapper
    args = build_parser().parse_args(argv)
    ledger = replace_library(
        args.source_dir,
        args.out_dir,
        manifest_path=args.manifest,
        nx=args.nx,
        nz=args.nz,
        volume_tolerance=args.volume_tolerance,
        copy_unreplaced=not args.no_copy_unreplaced,
    )
    summary = ledger["summary"]
    print(
        f"scanned {summary['files_scanned']} files: {len(summary['placeholders'])} placeholder, "
        f"{len(summary['duplicates'])} duplicate, {len(summary['valid'])} valid -> "
        f"{len(summary['replaced'])} replaced, {summary['withheld']} withheld"
    )
    for block in ledger["files"]:
        print(f"  {block['file']:<16} {block['classification']:<12} {block['status']:<9} "
              f"{block.get('verification_conclusion') or block.get('reason')}")
    print(f"ledger: {Path(args.out_dir) / 'ledger.json'}")
    return 0 if summary["all_replacements_gated"] else 1


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())