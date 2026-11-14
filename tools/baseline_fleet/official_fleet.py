#!/usr/bin/env python3
"""Build the baseline fleet from the *official* published geometry.

    # 1. materialise the sources (network; digests pinned by tensorcad/fleet/official.py)
    python tools/baseline_fleet/fetch_sources.py --out artifacts/baseline_fleet

    # 2. read the published particulars out of the fetched pages
    python tools/baseline_fleet/official_fleet.py extract \\
        --pages artifacts/baseline_fleet/raw/pages \\
        --out artifacts/baseline_fleet_official

    # 3. build every hull from its official IGES through ``iges_field`` and gate it
    python tools/baseline_fleet/official_fleet.py build \\
        --src artifacts/baseline_fleet \\
        --out artifacts/baseline_fleet_official

What this closes
----------------
Earlier rounds admitted the fleet on *surrogate* geometry (Series-60 law / Beta
section-area law) because the official IGES archives could not be reached, and every
record carried ``retrieval="offline_recollection", verified=false``. This tool removes
both caveats for the four hulls that publish a hull IGES (KCS, KVLCC2, JBC, ONRT):
the digits are parsed from the publisher pages, the mesh is built from the publisher's
IGES, and each artefact carries the page sha256 and the archive member sha256.

What it deliberately does not do
--------------------------------
* It does not build DTC or Series 60 from a file: neither publishes a hull IGES, and
  inventing one is the defect this fleet exists to detect. They are reported ``WITHHELD``
  with the reason.
* It does not change the packaged manifest. The official fleet is written to
  ``artifacts/.../fleet_manifest_official.json``, so ``python -m tensorcad.fleet.cli``
  keeps working offline on the surrogate fleet while the official run is reproducible
  on a machine that has (or can fetch) the archives.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
import time
from dataclasses import replace
from pathlib import Path
from typing import Any, Sequence

if __package__ in (None, ""):  # allow ``python tools/baseline_fleet/official_fleet.py``
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from tensorcad.fleet.build import RouteRefusal, iges_field_mesh, write_binary_stl
from tensorcad.fleet.gates import VOLUME_TOLERANCE, evaluate_gates, measure
from tensorcad.fleet.manifest import (
    FleetManifest,
    ManifestError,
    Provenance,
    ShipRecord,
    load_manifest,
)
from tensorcad.fleet.official import (
    ARCHIVES,
    LICENCE,
    PAGES,
    PageExtraction,
    compare_with_record,
    consistency_check,
    parse_particulars,
    published_reference,
)
from tensorcad.fleet.sources import load_registry

SHIPS = ("KCS", "KVLCC2", "JBC", "ONRT")
#: Hulls with no official hull IGES in the published set, and why. Reported, never built.
WITHHELD = {
    "DTC": "the Duisburg Test Case page that used to carry the CAD (uni-due.de/ship/dtc.php) "
           "returns HTTP 404 as of this run and no DTC hull IGES is served by the workshop "
           "archives; the surrogate stays and is labelled as such",
    "Series60": "Series 60 is a published *family* (Todd 1963 sectional-area curves); no single "
                "official IGES exists for it, so its geometry stays the repository's Series-60 law",
}
#: Page file names as written by ``fetch_sources.py``.
PAGE_FILES = {
    "KCS": "t2015_kcs.html",
    "JBC": "t2015_jbc.html",
    "ONRT": "t2015_onrt.html",
    "KVLCC2": "simman_summary.html",
}


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


# --------------------------------------------------------------------------- #
# extraction
# --------------------------------------------------------------------------- #
def extract(args: argparse.Namespace) -> int:
    pages_dir = Path(args.pages)
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    packaged = {record.name: record for record in load_manifest().ships}
    payload: dict[str, Any] = {
        "schema": "tensorcad.fleet.official.extraction.v1",
        "pages_dir": str(pages_dir),
        "ships": {},
    }
    for ship in SHIPS:
        page = pages_dir / PAGE_FILES[ship]
        if not page.is_file():
            payload["ships"][ship] = {"ship": ship, "error": f"page not fetched: {page}"}
            print(f"{ship:<8} PAGE-MISSING {page}")
            continue
        raw = page.read_bytes()
        extraction = parse_particulars(ship, raw.decode("utf-8", "replace"), page_sha256=sha256_bytes(raw))
        fidelity = compare_with_record(packaged[ship], extraction) if ship in packaged else None
        block = {
            **extraction.as_dict(),
            "page_file": str(page),
            "page_sha256_matches_registry": extraction.page_sha256 == PAGES[ship]["sha256"],
            "consistency": consistency_check(extraction),
            "references": published_reference(ship, extraction),
            "compare_with_packaged_record": fidelity,
        }
        payload["ships"][ship] = block
        values = {name: item.value for name, item in extraction.values.items()}
        print(f"{ship:<8} {len(values):>2} digits  page_sha256={extraction.page_sha256[:16]}  "
              f"closure_rel={block['consistency'].get('relative_difference')}")
        if fidelity is not None:
            for item in fidelity["fields"]:
                flag = "ok  " if item["agrees"] else "DIFF"
                print(f"    {flag} {item['field']:<22} record={item['record']:<12g} "
                      f"published={item['published']:<12g} rel={item['relative_difference']:.2e}")
    (out_dir / "published_particulars.json").write_text(json.dumps(payload, indent=1) + "\n")
    print(f"wrote {out_dir / 'published_particulars.json'}")
    return 0


# --------------------------------------------------------------------------- #
# official records
# --------------------------------------------------------------------------- #
def official_record(ship: str, extraction: PageExtraction, *, geometry_member: str,
                    member_sha256: str) -> tuple[ShipRecord, dict[str, Any]]:
    """Turn an extraction into a fleet record whose provenance is fetched, not remembered.

    The only digit ever changed relative to the page is a *page-level decimal defect*, and it
    is changed only when the page's own closure test proves it: see :func:`resolve_digits`.
    """
    values, corrections = resolve_digits(ship, extraction)
    references = published_reference(ship, extraction)["candidates"]
    volume = next((item for item in references if item["basis"] != "closure"), None)
    basis = "literature_value_verified"
    if volume is None:
        basis = "closure"
    provenance = Provenance(
        source=f"{ARCHIVES[ship]['publisher']} workshop geometry and particulars",
        url=ARCHIVES[ship]["url"],
        license=LICENCE[ARCHIVES[ship]["publisher"]],
        retrieval="fetched",
        verified=True,
    )
    record = ShipRecord(
        name=ship,
        ship_class={
            "KCS": "KRISO Container Ship",
            "KVLCC2": "MOERI very large crude carrier (KRISO Tanker)",
            "JBC": "Japan Bulk Carrier",
            "ONRT": "ONR Tumblehome (ONR-5613-)",
        }[ship],
        length=values["length"],
        beam=values["beam"],
        draft=values["draft"],
        block_coefficient=values["block_coefficient"],
        displacement_basis=basis,
        provenance=provenance,
        displacement_volume=None if volume is None else float(volume["volume_m3"]),
        lcb_fraction=values.get("lcb_fraction"),
        geometry_route="iges_field",
        notes=(
            f"built from the publisher's IGES ({ARCHIVES[ship]['member']}); particulars parsed "
            f"from {extraction.page_url} (page sha256 {extraction.page_sha256[:16]}…)"
            + ("; corrections: " + "; ".join(corrections) if corrections else "")
        ),
        midship_coefficient=values.get("midship_coefficient"),
        geometry_sha256=member_sha256,
        geometry_path=geometry_member,
    )
    return record, {
        "corrections": corrections,
        "references": references,
        "values": values,
    }


def resolve_digits(ship: str, extraction: PageExtraction) -> tuple[dict[str, Any], list[str]]:
    """Read the record digits off the extraction, repairing a page-level decimal defect.

    A repair is applied only when (a) the page's own closure ``Cb*L*B*T`` misses its published
    displacement by more than 1 %, and (b) dividing the suspect digit by ten closes to within
    0.5 %. The corpus case is KVLCC2's ``Bwl = 580`` (the site drops the decimal point of
    58.0): 580 m gives 3.13e6 m^3 against a published 312 622 m^3 -- 900 % off -- while 58.0
    closes to 0.16 %. Any repair is recorded verbatim in the ship's ``corrections`` list so a
    reader sees the page digit and the digit admitted.
    """
    get = extraction.get
    length = get("lpp") or get("lwl")
    values: dict[str, Any] = {
        "length": length,
        "beam": get("bwl"),
        "draft": get("draft"),
        "block_coefficient": get("cb"),
        "midship_coefficient": get("cm"),
    }
    lcb_percent = get("lcb_percent")
    if lcb_percent is not None:
        values["lcb_fraction"] = lcb_percent / 100.0
    lcb_metres = get("lcb_metres_aft_fp")
    if lcb_metres is not None and length:
        # Published as metres aft of FP, positive aft: fwd-positive fraction is the negative.
        values["lcb_fraction"] = -lcb_metres / length
    corrections: list[str] = []
    check = consistency_check(extraction)
    published = get("displacement_volume")
    if published is None:
        mass = get("displacement_mass")
        published = None if mass is None else mass / 1.025
    if published and values["beam"] and check.get("relative_difference", 0.0) > 0.01:
        for factor in (0.1, 10.0):
            candidate = values["beam"] * factor
            closure = values["block_coefficient"] * length * candidate * values["draft"]
            if abs(closure - published) / published <= 0.005:
                corrections.append(
                    f"beam: page prints {extraction.get('bwl')!r} m; admitted {candidate:g} m because "
                    f"Cb*L*B*T closes to {closure:.1f} m^3 against the published {published:.1f} m^3 "
                    f"({abs(closure - published) / published * 100:.2f} %) while the page digit misses by "
                    f"{check['relative_difference'] * 100:.0f} %"
                )
                values["beam"] = candidate
                break
    return values, corrections


def official_manifest(records: Sequence[ShipRecord], *, placeholder: dict[str, Any]) -> dict[str, Any]:
    return {
        "schema": "tensorcad.fleet.baseline.v1",
        "anchor": "GAP-CAD-T1-BASELINE-FLEET",
        "generated_by": "tools/baseline_fleet/official_fleet.py",
        "purpose": "Baseline fleet built from the publishers' own hull IGES; every digit parsed from "
                   "the publisher page and every artefact pinned by sha256.",
        "honesty_note": "Only the four hulls that publish a hull IGES are in this manifest. DTC and "
                        "Series 60 are withheld with a stated reason rather than built from an invented "
                        "file. The IGES bytes are not redistributed (publishers state no licence); "
                        "tools/baseline_fleet/fetch_sources.py re-materialises them and re-hashes them "
                        "against the recorded digests.",
        "placeholder_signature": placeholder,
        "ships": [record_as_dict(record) for record in records],
    }


def record_as_dict(record: ShipRecord) -> dict[str, Any]:
    return {
        "name": record.name,
        "ship_class": record.ship_class,
        "length": record.length,
        "beam": record.beam,
        "draft": record.draft,
        "block_coefficient": record.block_coefficient,
        "displacement_volume": record.displacement_volume,
        "displacement_basis": record.displacement_basis,
        "lcb_fraction": record.lcb_fraction,
        "midship_coefficient": record.midship_coefficient,
        "geometry_route": record.geometry_route,
        "geometry_path": record.geometry_path,
        "geometry_sha256": record.geometry_sha256,
        "transom_ratio": record.transom_ratio,
        "provenance": record.provenance.as_dict(),
        "notes": record.notes,
    }


def iges_path_for(ship: str, src: Path) -> Path:
    """Locate the hull IGES under ``src`` in either cached layout.

    Earlier rounds staged the archives as ``raw/<SHIP>/<file>.igs``; ``fetch_sources.py``
    extracts them as ``iges/<member>``. Both are accepted, and the *first existing* path wins,
    so the tool never guesses a file it does not have (a missing source is a WITHHELD ship,
    not a crash).
    """
    member = ARCHIVES[ship]["member"]
    candidates = [src / "iges" / member, src / "raw" / ship / Path(member).name]
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    return candidates[0]


# --------------------------------------------------------------------------- #
# build + gate
# --------------------------------------------------------------------------- #
def build(args: argparse.Namespace) -> int:
    src = Path(args.src)
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    meshes = out_dir / "meshes"
    pages_dir = src / "raw" / "pages"
    packaged = load_manifest()
    registry = load_registry()

    records: list[ShipRecord] = []
    extractions: dict[str, PageExtraction] = {}
    per_ship: dict[str, Any] = {}
    for ship in SHIPS:
        page = pages_dir / PAGE_FILES[ship]
        iges = iges_path_for(ship, src)
        if not page.is_file() or not iges.is_file():
            per_ship[ship] = {"ship": ship, "status": "withheld",
                              "reason": f"source absent (page={page.is_file()}, iges={iges.is_file()})"}
            print(f"{ship:<8} WITHHELD (source absent)")
            continue
        extraction = parse_particulars(
            ship, page.read_text(encoding="utf-8", errors="replace"), page_sha256=sha256_file(page)
        )
        member_sha = sha256_file(iges)
        record, context = official_record(
            ship, extraction, geometry_member=str(iges.relative_to(src)), member_sha256=member_sha
        )
        records.append(record)
        extractions[ship] = extraction
    manifest = FleetManifest(
        ships=tuple(records), placeholder=dict(packaged.placeholder), path=None, sha256=None
    )
    (out_dir / "fleet_manifest_official.json").write_text(
        json.dumps(official_manifest(records, placeholder=dict(packaged.placeholder)), indent=1) + "\n"
    )

    for ship in args.ships.split(",") if args.ships else SHIPS:
        if ship not in {record.name for record in records}:
            continue
        record = next(item for item in records if item.name == ship)
        iges = iges_path_for(ship, src)
        started = time.time()
        try:
            built = iges_field_mesh(
                record,
                iges,
                dx=args.dx,
                dz=args.dz,
                ring_points=args.ring_points,
                samples_per_metre=args.samples_per_metre,
            )
        except Exception as exc:  # a route/parse failure is evidence, not a crash of the run
            per_ship[ship] = {
                "ship": ship,
                "status": "refused",
                "reason": f"{type(exc).__name__}: {exc}",
            }
            print(f"{ship:<8} REFUSED {type(exc).__name__}: {exc}")
            continue
        elapsed = time.time() - started
        facts = measure(built.mesh)
        stl_path = meshes / f"{ship}.stl"
        meshes.mkdir(parents=True, exist_ok=True)
        stl_sha = write_binary_stl(built.mesh, stl_path)
        verdict = evaluate_gates(record, built.mesh, manifest, registry=registry,
                                 volume_tolerance=args.volume_tolerance)
        # independent integration through a different code path
        from tensorcad.hydrostatics import compute_hydrostatics

        hydro = compute_hydrostatics(built.mesh)
        hydro_volume = getattr(hydro, "displacement_volume", None)
        per_ship[ship] = {
            "ship": ship,
            "status": "admitted" if verdict["passed"] else "failed_gates",
            "route": built.route,
            "record": record_as_dict(record),
            "context": context,
            "page_sha256": extractions[ship].page_sha256,
            "page_url": extractions[ship].page_url,
            "geometry_source": {
                "url": ARCHIVES[ship]["url"],
                "member": ARCHIVES[ship]["member"],
                "member_sha256": record.geometry_sha256,
                "bytes": iges.stat().st_size,
                "redistributed": False,
            },
            "build": {
                "dx": args.dx, "dz": args.dz, "ring_points": args.ring_points,
                "samples_per_metre": getattr(args, "samples_per_metre", None),
                "seconds": round(elapsed, 1),
                "details": {k: v for k, v in built.details.items() if k != "field"},
                "field": {k: v for k, v in built.details.get("field", {}).items() if k != "b"},
            },
            "topology": facts.as_dict(),
            "gates": verdict["gates"],
            "passed": verdict["passed"],
            "failed_gates": verdict["failed_gates"],
            "volume_relative_error": verdict["volume_relative_error"],
            "reference_basis": verdict["reference_basis"],
            "stl": {"path": str(stl_path), "sha256": stl_sha},
            "independent_volume": {
                "hydrostatics_displacement_volume": hydro_volume,
                "relative_difference": (
                    None if not hydro_volume else abs(hydro_volume - facts.volume) / facts.volume
                ),
            },
        }
        print(
            f"{ship:<8} {built.route:<11} tris={facts.triangle_count:>7} "
            f"V={facts.volume:>12.1f} ref={record.reference_volume():>12.1f} "
            f"err={1e2 * (verdict['volume_relative_error'] or 0):>7.3f}% "
            f"watertight={facts.boundary_edge_count == 0} failed={verdict['failed_gates']} "
            f"({elapsed:.0f}s)"
        )

    summary = {
        "schema": "tensorcad.fleet.official.build.v1",
        "anchor": "GAP-CAD-T1-BASELINE-FLEET",
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "settings": {
            "dx": args.dx, "dz": args.dz, "ring_points": args.ring_points,
            "samples_per_metre": getattr(args, "samples_per_metre", None),
            "volume_tolerance": args.volume_tolerance,
        },
        "ships": per_ship,
        "withheld": WITHHELD,
        "admitted": sorted(name for name, block in per_ship.items() if block.get("passed")),
        "failed": sorted(name for name, block in per_ship.items() if block.get("status") == "failed_gates"),
    }
    (out_dir / "official_build.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(f"wrote {out_dir / 'official_build.json'}")
    print(f"admitted: {summary['admitted']}  failed: {summary['failed']}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)

    p_extract = sub.add_parser("extract", help="parse published particulars out of fetched pages")
    p_extract.add_argument("--pages", default="artifacts/baseline_fleet/raw/pages")
    p_extract.add_argument("--out", default="artifacts/baseline_fleet_official")
    p_extract.set_defaults(func=extract)

    p_build = sub.add_parser("build", help="build + gate the official IGES fleet")
    p_build.add_argument("--src", default="artifacts/baseline_fleet",
                         help="output of tools/baseline_fleet/fetch_sources.py (iges/, raw/pages/)")
    p_build.add_argument("--out", default="artifacts/baseline_fleet_official")
    p_build.add_argument("--dx", type=float, default=0.5)
    p_build.add_argument("--dz", type=float, default=0.1)
    p_build.add_argument("--ring-points", type=int, default=240)
    p_build.add_argument("--samples-per-metre", type=float, default=8.0)
    p_build.add_argument("--volume-tolerance", type=float, default=VOLUME_TOLERANCE)
    p_build.add_argument("--ships", default="",
                         help="comma-separated subset to build (default: every fetched hull)")
    p_build.set_defaults(func=build)

    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except (ManifestError, RouteRefusal) as exc:
        print(f"error: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())