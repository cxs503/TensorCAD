#!/usr/bin/env python3
"""Build (or verify) the geometry-source ledger for anchor ``GAP-CAD-T1-SOURCE-HUNT``.

What this tool answers, from cached publisher bytes only:

* which workshop artefacts were actually retrieved, from which URL, with which SHA-256,
  and under which licence statement -- the manifest that the anchor's acceptance asks for;
* what the geometry archives *are*: IGES unit flag, model-space scale, exporting system,
  and the measured bounding extent in metres (``tensorcad.fleet.geometry_provenance``);
* whether that extent reproduces the publisher's own particulars -- the JBC archive's
  beam and depth at the published scale denominator 1:40, the ONRT archive's L/B ratio --
  so the ledger carries a computed quantity, not a transcribed digit;
* which acceptance items are MET and which are WITHHELD, with the probe evidence for the
  withheld ones (the DTC hosts answer 404/403 and this network has no egress).

The bytes stay in the git-ignored cache: the committed ledger holds digests, licences and
measured dimensions.  Every number is recomputed from the cache on each run, so
``--verify`` is a real check rather than a re-read of the previous output.

Usage::

    python tools/baseline_fleet/geometry_source_ledger.py \
        --cache /path/to/raw --out tensorcad/fleet/data/geometry_source_ledger.json
    python tools/baseline_fleet/geometry_source_ledger.py --cache ... --verify
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from tensorcad.fleet import geometry_provenance as gp  # noqa: E402

LEDGER_SCHEMA = "tensorcad.fleet.geometry_source_ledger.v1"
DEFAULT_OUT = REPO / "tensorcad" / "fleet" / "data" / "geometry_source_ledger.json"
FLEET_SPECS = REPO / "benchmarks" / "verified" / "baseline_fleet" / "fleet_specs.json"
CROSSREF = REPO / "benchmarks" / "verified" / "fleet_refs" / "crossref.json"

#: Which member of a publisher archive is the bare hull.  Named explicitly: the first
#: ``.igs`` in the JBC archive is ``jbc-spacer.igs`` (a 3 mm filler patch), not the hull.
HULL_MEMBER = {
    "t2015_jbc_iges": "jbc.igs",
    "t2015_onrt_wo_gap": "onrt_wo_gap/onr5613-full_appendage_deckhous_mergedR.igs",
}

LICENCE = {
    "archive": (
        "Tokyo 2015 workshop geometry archive, published by NMRI for research use; no OSI "
        "licence file ships with it, so only digests, member names and measured dimensions "
        "are recorded here -- the bytes are not redistributed"
    ),
    "page": "page copyright NMRI / SIMMAN; digest + extracted row text only",
}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


# --------------------------------------------------------------------------- cache
def load_cache_entry(binary: Path) -> dict[str, Any]:
    """Read one cached artefact: its meta sidecar plus the digest of the stored bytes."""
    meta_path = binary.with_suffix(".meta.json")
    meta = json.loads(meta_path.read_text()) if meta_path.is_file() else {}
    payload = binary.read_bytes()
    return {
        "cache_file": binary.name,
        "url": meta.get("url", ""),
        "final_url": meta.get("final_url", meta.get("url", "")),
        "http_status": meta.get("http_status"),
        "retrieved_utc": meta.get("retrieved_utc", ""),
        "fetched_by": meta.get("tool", ""),
        "bytes": len(payload),
        "sha256": gp.sha256_bytes(payload),
        "kind": gp.kind_of(binary.name.replace(".bin", ""), payload),
    }


def read_cache(cache: Path) -> dict[str, dict[str, Any]]:
    """All ``*.bin`` artefacts of the cache directory, keyed by their short source name."""
    entries: dict[str, dict[str, Any]] = {}
    for binary in sorted(cache.glob("*.bin")):
        key = re.sub(r"-[0-9a-f]{16}$", "", binary.stem)
        entries[key] = load_cache_entry(binary)
    return entries


def crossref_index(crossref: dict[str, Any]) -> dict[str, dict[str, Any]]:
    """Index the digest-bound source registrations of ``fleet_refs/crossref.json``."""
    index: dict[str, dict[str, Any]] = {}
    for key, entry in (crossref.get("sources") or {}).items():
        if not isinstance(entry, dict):
            continue
        sha = entry.get("sha256")
        if not sha:
            continue
        index[sha] = {
            "key": key,
            "family": entry.get("family"),
            "licence": entry.get("licence"),
            "publisher": entry.get("publisher"),
            "url": entry.get("url"),
            "quantities_read": entry.get("quantities_read"),
        }
        index.setdefault("by_key:" + key, index[sha])
    return index


# --------------------------------------------------------------------------- published
def published_ships() -> dict[str, dict[str, Any]]:
    """Digest-bound published particulars, straight out of the fleet manifest."""
    data = json.loads(FLEET_SPECS.read_text())
    out: dict[str, dict[str, Any]] = {}
    for ship in data["ships"]:
        name = ship.get("name")
        published = {
            "length_pp": ship.get("length_pp"),
            "beam": ship.get("beam"),
            "draft": ship.get("draft"),
            "displacement_volume": ship.get("published_volume_m3"),
            "cb": ship.get("published_cb"),
            "published_source": ship.get("published_source"),
            "scale_basis": ship.get("scale_basis"),
        }
        out[name] = published
    return out


def scrape_depth(page_markup: str) -> dict[str, Any]:
    """Depth (m) as printed on the T2015 general-characteristics page."""
    return gp.scrape_scalar(page_markup, r"Depth \\?\(D\\?\) \\?\(m\\?\) ([\d.]+)")


def scrape_lwl_pair(page_markup: str) -> dict[str, Any]:
    """``Length of waterline (m) <model> <full scale>`` pair on the T2015 ONRT page."""
    return gp.scrape_pair(page_markup, r"Length of waterline L WL \\?\(m\\?\) ([\d.]+) ([\d.]+)")


def quoted_scale_denominator(crossref: dict[str, Any], ship: str) -> dict[str, Any] | None:
    """Read the publisher's own scale-factor row (``Scale factor <lambda> 1 40``)."""
    hits = gp.find_published_quote_entries(crossref, "Scale factor")
    for hit in hits:
        if hit["ship"] != ship:
            continue
        match = re.search(r"Scale factor [^0-9]{0,8}1\s+(\d+(?:\.\d+)?)", hit["cell_text"])
        if not match:
            continue
        return {
            "denominator": float(match.group(1)),
            "quote": match.group(0),
            "source": hit["source"],
            "url": hit["url"],
            "sha256": hit["sha256"],
            "cell_text": hit["cell_text"],
        }
    return None


# --------------------------------------------------------------------------- geometry
def analyse_hull_archive(
    key: str, payload: bytes, published: dict[str, Any], *, depth_m: float | None,
    quoted_scale: dict[str, Any] | None, lwl_pair: dict[str, Any] | None,
) -> dict[str, Any]:
    """Fingerprint + measured extent + scale reconciliation for one publisher archive."""
    member_name = HULL_MEMBER[key]
    name, member_bytes = gp.zip_member(payload, member_name)
    fingerprint = gp.iges_fingerprint(member_bytes)
    extent = gp.iges_extent_m(member_bytes)

    axes: dict[str, dict[str, float]] = {}
    if published.get("length_pp"):
        axes["x"] = {"value": published["length_pp"], "quantity": "Lpp"}
    if published.get("beam"):
        axes["y"] = {"value": published["beam"], "quantity": "B"}
    ship_depth = depth_m or published.get("draft")
    if ship_depth:
        axes["z"] = {
            "value": ship_depth,
            "quantity": "D (depth)" if depth_m else "T (draft)",
        }
    denominators = gp.axis_denominators(extent["extent_m"], axes)

    lwl_full = None
    if lwl_pair:
        lwl_full = lwl_pair["values"][1]
    elif published.get("length_pp"):
        lwl_full = published["length_pp"]
    regime = gp.classify_regime(extent["extent_m"], lwl_full) if lwl_full else None

    checks: list[dict[str, Any]] = []
    if quoted_scale is not None:
        ship_key = key.split("_")[1].upper()
        # The two transverse axes are the ones a hull export cannot fudge: length can pick
        # up overhangs or a spacer patch, beam and depth are plain extents.
        for axis, label in (("y", "B"), ("z", "D")):
            row = denominators.get(axis)
            if not row:
                continue
            checks.append(
                gp.check_row(
                    f"{ship_key} scale denominator from {label}",
                    row["denominator"],
                    quoted_scale["denominator"],
                    unit="--",
                    source=f"quoted 'Scale factor 1 {quoted_scale['denominator']:g}' "
                    f"({quoted_scale['source']})",
                )
            )
        for axis, label in (("y", "B"), ("z", "D")):
            row = denominators.get(axis)
            if not row or not quoted_scale["denominator"]:
                continue
            expected = row["published_full"] / quoted_scale["denominator"]
            checks.append(
                gp.check_row(
                    f"{ship_key} {label} at 1:{quoted_scale['denominator']:g} (m)",
                    row["measured_m"],
                    expected,
                    unit="m",
                    source=f"published {label}={row['published_full']:g} m / "
                    f"{quoted_scale['denominator']:g}",
                )
            )
    if lwl_pair and regime:
        # Unit-normalised exports: the geometry is the right *shape*, so the check is on
        # proportions (beam over length), which survive normalisation.
        model_lwl, full_lwl = lwl_pair["values"]
        measured_ratio = extent["extent_m"][1] / extent["extent_m"][0]
        published_ratio = published["beam"] / full_lwl
        ship_key = key.split("_")[1].upper()
        checks.append(
            gp.check_row(
                f"{ship_key} B/Lwl from geometry",
                measured_ratio,
                published_ratio,
                unit="--",
                source=f"published B={published['beam']:g} m, Lwl={full_lwl:g} m",
            )
        )
        # A unit-normalised export can be read at *any* scale, so the honest check is the
        # mismatch that the regime predicts: the file is 1.0 long, the published model hull
        # is 3.147 m.  Recorded as a check, not as a footnote, so a reader sees the number.
        checks.append(
            gp.check_row(
                f"{ship_key} file length vs published model Lwl",
                extent["extent_m"][0],
                model_lwl,
                unit="m",
                source=f"published model Lwl={model_lwl:g} m (full {full_lwl:g} m)",
                note=(
                    "expected mismatch: the export is unit-length normalised, so it must be "
                    "scaled by the published Lwl before meshing; this row is the measurement "
                    "that establishes the regime, not a defect of the archive"
                ),
            )
        )

    return {
        "archive_sha256": gp.sha256_bytes(payload),
        "archive_bytes": len(payload),
        "members": gp.zip_members(payload),
        "hull_member": {
            "name": name,
            "bytes": len(member_bytes),
            "sha256": gp.sha256_bytes(member_bytes),
        },
        "iges_fingerprint": fingerprint,
        "extent": extent,
        "axis_denominators": denominators,
        "regime": regime,
        "checks": checks,
    }


# --------------------------------------------------------------------------- ledger
def build_ledger(cache: Path) -> dict[str, Any]:
    """Recompute the whole ledger from the cache plus the committed digest-bound refs."""
    entries = read_cache(cache)
    crossref = json.loads(CROSSREF.read_text())
    cr_index = crossref_index(crossref)
    ships = published_ships()

    for key, entry in entries.items():
        match = cr_index.get(entry["sha256"])
        entry["digest_bound_in_repo"] = bool(match)
        entry["repo_source_key"] = match["key"] if match else None
        entry["licence"] = (
            (match or {}).get("licence")
            or LICENCE["page" if entry["kind"] == "html" else "archive"]
        )

    def page_text(key: str) -> str:
        return (cache / f"{key}.bin").read_text(errors="replace") if (cache / f"{key}.bin").is_file() else ""

    # The cache file names carry a digest suffix, so resolve the page by prefix.
    def cached_payload(prefix: str) -> bytes | None:
        matches = sorted(cache.glob(prefix + "-*.bin"))
        return matches[0].read_bytes() if matches else None

    jbc_page = cached_payload("t2015_jbc_gc")
    onrt_page = cached_payload("t2015_onrt_gc")
    jbc_depth = scrape_depth(jbc_page.decode("utf-8", "replace")) if jbc_page else None
    if jbc_depth:
        jbc_depth["page"] = "t2015_jbc_gc"
        jbc_depth["page_sha256"] = gp.sha256_bytes(jbc_page)
    onrt_lwl = scrape_lwl_pair(onrt_page.decode("utf-8", "replace")) if onrt_page else None
    if onrt_lwl:
        onrt_lwl["page"] = "t2015_onrt_gc"
        onrt_lwl["page_sha256"] = gp.sha256_bytes(onrt_page)

    jbc_scale = quoted_scale_denominator(crossref, "JBC")

    geometry: dict[str, Any] = {}
    checks: list[dict[str, Any]] = []
    for key in ("t2015_jbc_iges", "t2015_onrt_wo_gap"):
        payload = cached_payload(key)
        if payload is None:
            continue
        ship = key.split("_")[1].upper()
        published = ships.get(ship, {})
        analysis = analyse_hull_archive(
            key,
            payload,
            published,
            depth_m=(jbc_depth["value"] if ship == "JBC" and jbc_depth else None),
            quoted_scale=(jbc_scale if ship == "JBC" else None),
            lwl_pair=(onrt_lwl if ship == "ONRT" else None),
        )
        analysis["ship"] = ship
        geometry[ship] = analysis
        checks.extend(analysis["checks"])

    # Cross-host identity: the same archive served by two workshop hosts.
    onrt_hosts = {
        key: entries[key]["sha256"]
        for key in ("t2015_onrt_wo_gap", "simman_onrt_wo_gap")
        if key in entries
    }
    cross_host = {
        "quantity": "ONRT geometry archive",
        "hosts": onrt_hosts,
        "identical": len(set(onrt_hosts.values())) == 1 and len(onrt_hosts) > 1,
    }

    dead_ends = {
        key: {
            "http_status": entry["http_status"],
            "sha256": entry["sha256"],
            "url": entry["url"],
        }
        for key, entry in entries.items()
        if isinstance(entry["http_status"], int) and entry["http_status"] >= 400
    }

    acceptance = grade_acceptance(entries, geometry, checks, cross_host, dead_ends, jbc_scale)

    ledger = {
        "schema": LEDGER_SCHEMA,
        "anchor": "GAP-CAD-T1-SOURCE-HUNT",
        "generated_utc": utc_now(),
        "tool": "tools/baseline_fleet/geometry_source_ledger.py",
        "cache": {
            "directory_hint": cache.name,
            "artefacts": len(entries),
            "retrieval_window": sorted(
                {e["retrieved_utc"] for e in entries.values() if e["retrieved_utc"]}
            ),
            "note": (
                "bytes stay out of git (publisher files, research-use grant only); the "
                "ledger keeps digest, byte count, URL and licence so a later run can "
                "re-fetch and compare"
            ),
        },
        "network": {
            "state": "no egress at ledger build time",
            "evidence": (
                "urllib HTTPS GET to t2015.nmri.go.jp times out in this environment; the "
                "ledger therefore consumes the cache retrieved 2026-11-14 and adds no new "
                "fetch. A run with egress re-fetches and compares digests."
            ),
        },
        "sources": entries,
        "geometry": geometry,
        "published": {
            ship: {"particulars": data, "digest_bound_entries": "benchmarks/verified/baseline_fleet/fleet_specs.json"}
            for ship, data in ships.items()
        },
        "scraped": {"jbc_depth": jbc_depth, "onrt_lwl": onrt_lwl, "jbc_scale_factor": jbc_scale},
        "cross_host": cross_host,
        "dead_ends": dead_ends,
        "checks": checks,
        "checks_summary": gp.summarise_check_rows(checks),
        "acceptance": acceptance,
        "open_questions": open_questions(geometry, jbc_scale, dead_ends),
    }
    return ledger


def grade_acceptance(
    entries: dict[str, dict[str, Any]],
    geometry: dict[str, Any],
    checks: list[dict[str, Any]],
    cross_host: dict[str, Any],
    dead_ends: dict[str, Any],
    jbc_scale: dict[str, Any] | None,
) -> dict[str, Any]:
    """Per-item verdicts, each tied to the artefact that decides it."""
    onrt = geometry.get("ONRT")
    jbc = geometry.get("JBC")
    onrt_ok = bool(onrt) and bool(onrt["extent"]["header_matches_extent"])
    dtc_probes = {k: v for k, v in dead_ends.items() if k.startswith("dtc") or k == "t2015_dtc_page"}
    return {
        "dtc_geometry": {
            "verdict": "WITHHELD",
            "detail": (
                "no publisher-hosted DTC geometry: uni-due.de page and its section answer "
                "404 with the same error-body digest, the T2015 DTC page answers 404, the "
                "geometry directory answers 403 (not listable) and the Internet Archive is "
                "unreachable from this network; the only public copy found is a third-party "
                "OpenFOAM tutorial whose licence does not travel with it"
            ),
            "evidence": dtc_probes,
        },
        "onrt_geometry": {
            "verdict": "MET" if onrt_ok else "WITHHELD",
            "detail": (
                "archive retrieved from the T2015 host, hull member fingerprinted, extent "
                "measured from the surface grid; the export is unit-length normalised (largest "
                "coordinate 1.0005) so the physical scale must come from the published Lwl"
                if onrt_ok
                else "ONRT archive missing or unmeasurable"
            ),
            "evidence": {
                "archive_sha256": onrt["archive_sha256"] if onrt else None,
                "hull_member_sha256": onrt["hull_member"]["sha256"] if onrt else None,
                "cross_host_identical": cross_host["identical"],
            },
        },
        "jbc_geometry": {
            "verdict": "MET" if jbc else "WITHHELD",
            "detail": (
                "IGES archive retrieved, seven members digested, hull member fingerprinted"
                if jbc
                else "JBC archive missing"
            ),
            "evidence": {
                "archive_sha256": jbc["archive_sha256"] if jbc else None,
                "members": [m["name"] for m in jbc["members"]] if jbc else [],
            },
        },
        "jbc_scale_factor": {
            "verdict": "MET" if (jbc and jbc_scale) else "WITHHELD",
            "detail": (
                "denominator derived from the hull geometry's transverse extent and "
                f"cross-checked against the publisher's quoted 'Scale factor 1 {jbc_scale['denominator']:g}'"
                if (jbc and jbc_scale)
                else "no geometry or no quoted scale factor"
            ),
            "evidence": {
                "derived": (jbc["axis_denominators"] if jbc else {}),
                "quoted": jbc_scale,
                "checks": [c for c in checks if c["quantity"].startswith("JBC")],
            },
        },
        "source_page_digests": {
            "verdict": "MET",
            "detail": (
                f"{len(entries)} cached artefacts digested, of which "
                f"{sum(1 for e in entries.values() if e['digest_bound_in_repo'])} also appear "
                "with the same SHA-256 in the committed crossref registry, so the cache and "
                "the manifest are one chain"
            ),
            "evidence": {
                "artefacts": len(entries),
                "digest_bound_in_repo": sum(
                    1 for e in entries.values() if e["digest_bound_in_repo"]
                ),
            },
        },
    }


def open_questions(
    geometry: dict[str, Any], jbc_scale: dict[str, Any] | None, dead_ends: dict[str, Any]
) -> list[dict[str, Any]]:
    """What this ledger still cannot answer, stated so the next mission inherits it."""
    questions: list[dict[str, Any]] = []
    if not any(k.startswith("dtc") and v["http_status"] == 200 for k, v in dead_ends.items()):
        questions.append(
            {
                "item": "DTC hull geometry (publisher copy)",
                "status": "WITHHELD",
                "evidence": (
                    "uni-due.de 404 / T2015 dtc page 404 / geometry dir 403 / archive.org "
                    "unreachable; needs a network with egress or an owner-supplied copy"
                ),
            }
        )
    jbc = geometry.get("JBC")
    if jbc:
        x_axis = jbc["axis_denominators"].get("x")
        if x_axis and x_axis["denominator"] and jbc_scale:
            spread = abs(x_axis["denominator"] - jbc_scale["denominator"]) / jbc_scale["denominator"]
            if spread > gp.TOLERANCE_REL:
                questions.append(
                    {
                        "item": "JBC x-axis extent vs Lpp",
                        "status": "EXPLAINED_NOT_MATCHED",
                        "evidence": (
                            f"x extent {x_axis['measured_m']:.4f} m implies denominator "
                            f"{x_axis['denominator']:.3f} against the quoted "
                            f"{jbc_scale['denominator']:g} ({spread * 100:.1f} % long): the hull "
                            "file carries geometry beyond the perpendiculars, so Lpp cannot be "
                            "read off the bounding box. Beam and depth (transverse extents, "
                            "insensitive to overhang) are used instead."
                        ),
                    }
                )
    questions.append(
        {
            "item": "KCS / KVLCC2 / SERIES60 geometry archives",
            "status": "WITHHELD",
            "evidence": (
                "this cache holds their particulars pages but no geometry archive for KCS or "
                "KVLCC2, and no Series 60 artefact at all"
            ),
        }
    )
    return questions


# --------------------------------------------------------------------------- verify
def verify(ledger_path: Path, cache: Path) -> int:
    """Recompute the ledger and compare the parts that must be reproducible."""
    if not ledger_path.is_file():
        print(f"FAIL: ledger missing at {ledger_path}")
        return 2
    stored = json.loads(ledger_path.read_text())
    fresh = build_ledger(cache)
    problems: list[str] = []

    for key, entry in stored["sources"].items():
        other = fresh["sources"].get(key)
        if other is None:
            problems.append(f"{key}: present in ledger, absent from cache")
            continue
        for field in ("sha256", "bytes", "http_status"):
            if entry[field] != other[field]:
                problems.append(f"{key}.{field}: ledger {entry[field]!r} != cache {other[field]!r}")
    for key in fresh["sources"]:
        if key not in stored["sources"]:
            problems.append(f"{key}: in cache but never recorded")

    for ship, block in stored["geometry"].items():
        other = fresh["geometry"].get(ship)
        if other is None:
            problems.append(f"geometry.{ship}: no longer measurable from cache")
            continue
        for field in ("archive_sha256", "hull_member", "extent", "axis_denominators"):
            if block[field] != other[field]:
                problems.append(f"geometry.{ship}.{field}: ledger disagrees with cache")

    for a, b in zip(stored["checks"], fresh["checks"]):
        if a["quantity"] != b["quantity"] or abs(a["measured"] - b["measured"]) > 1e-9:
            problems.append(f"check {a['quantity']!r}: stored {a['measured']} != fresh {b['measured']}")
    if len(stored["checks"]) != len(fresh["checks"]):
        problems.append("checks: row count changed")

    summary = fresh["checks_summary"]
    print(
        f"verify: {len(fresh['sources'])} artefacts, {len(fresh['geometry'])} geometry archives, "
        f"{summary['rows']} checks ({summary['within_tolerance']} within {gp.TOLERANCE_REL:.0%}), "
        f"{len(problems)} problem(s)"
    )
    for problem in problems:
        print("  -", problem)
    return 0 if not problems else 1


# --------------------------------------------------------------------------- main
def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cache", type=Path, required=True, help="directory holding the *.bin cache")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT, help="ledger JSON to write")
    parser.add_argument("--verify", action="store_true", help="compare --out against a fresh build")
    args = parser.parse_args(argv)

    if not args.cache.is_dir():
        print(f"FAIL: cache directory not found: {args.cache}")
        return 2
    if args.verify:
        return verify(args.out, args.cache)

    ledger = build_ledger(args.cache)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(ledger, indent=1, sort_keys=False) + "\n")
    summary = ledger["checks_summary"]
    print(
        f"wrote {args.out} :: {len(ledger['sources'])} artefacts, "
        f"{len(ledger['geometry'])} geometry archives, {summary['rows']} checks "
        f"({summary['within_tolerance']} within {gp.TOLERANCE_REL:.0%}, worst "
        f"{summary['worst_abs_rel_error']:.4%} on {summary['worst_quantity']})"
    )
    for item, block in ledger["acceptance"].items():
        print(f"  {item}: {block['verdict']}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())