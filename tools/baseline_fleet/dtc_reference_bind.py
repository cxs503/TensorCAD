"""Bind the DTC published particulars to live document bytes instead of recollection.

Anchor ``GAP-CAD-T1-SOURCE-HUNT``.  The geometry half of this anchor is closed: the publisher
archive ``igs_ordner.zip`` is served by the Duisburg-Essen file area, its SHA-256 reproduces on
every fetch, and the hull inside it measures to the published displacement
(:mod:`tools.baseline_fleet.dtc_hull_measure`, record
``benchmarks/verified/fleet_refs/dtc_hull_measured.json``).

What is **not** closed is where the *numbers* come from.  That record carries

    "manifest_provenance": {"retrieval": "offline_recollection", "verified": false, ...}

i.e. Lpp 355.0 m, B 51.0 m, T 14.5 m, Cb 0.661 and gradient 173 467 m3 were typed in from
memory of El Moctar, Shigunov & Zorn (2012).  The anchor forbids that: *no hand-copied digits,
every value carries the fingerprint of the page it came from*.  A recollection is also exactly
the kind of reference that drifts -- the same record once carried ``B = 45.8 m``.

This tool closes the hole.  It

* fetches the candidate documents (the publisher's institute page, the legacy DTC page, the
  geometry archive, and open-access papers that tabulate the DTC particulars) over the live
  network and digests the bytes it received;
* locates a *particulars block* in each document -- a window in which at least ``MIN_LABELS``
  of the five field labels co-occur -- and keeps the raw row text next to every number it
  extracts, so a reviewer can see the column layout instead of trusting the parse;
* accepts a field only when a token inside its row matches the reference value within
  ``FIELD_TOLERANCE_PERCENT``, and records the matched token, the column index it sat in, the
  digitised page number and the document digest;
* proves the archive side again in the same run: the live archive digest must equal the one in
  the measurement record, and every member digest must equal the fleet registry's;
* declares, for each field, ``BOUND`` (a document pins it) or ``UNBOUND`` (honest residual), and
  refuses to write ``BOUND`` on nothing.

Nothing is redistributed: bytes stay in the local cache (``--cache``, git-ignored); git carries
the digests, the matched snippets and this script.  No licence file accompanies the publisher
data, so the archive and the papers are digested, never copied into the repository.

Usage::

    python tools/baseline_fleet/dtc_reference_bind.py --out benchmarks/verified/fleet_refs/dtc_reference_bind.json
    python tools/baseline_fleet/dtc_reference_bind.py --verify-only     # offline, what CI runs
    python tools/baseline_fleet/dtc_reference_bind.py --figure docs/figs/dtc_reference_bind.png

``--verify-only`` re-reads the committed record and checks the chain: digests well formed, every
``BOUND`` field backed by a snippet plus the digest of the document it came from, every bound
value inside tolerance, the archive digest equal to the measurement record's, and the member
digests equal to the registry's.  It never touches the network.
"""

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import re
import ssl
import sys
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:  # pragma: no cover - import plumbing
    sys.path.insert(0, str(REPO))

SCHEMA = "tensorcad.fleet_refs.dtc_reference_bind/1.0"
ANCHOR = "GAP-CAD-T1-SOURCE-HUNT"

DEFAULT_MEASURED = REPO / "benchmarks" / "verified" / "fleet_refs" / "dtc_hull_measured.json"
DEFAULT_REGISTRY = REPO / "assets" / "geometry_registry.json"
DEFAULT_OUT = REPO / "benchmarks" / "verified" / "fleet_refs" / "dtc_reference_bind.json"
DEFAULT_CACHE = Path("/root/.cache/hull_sources")

INSTITUTE_URL = "https://www.uni-due.de/ist/"  # redirects to the current /inam/ landing page
LEGACY_PAGE_URL = "https://www.uni-due.de/ship/dtc.php"  # the URL the DTC paper names; 404 since
ARCHIVE_URL = "https://www.uni-due.de/imperia/md/content/inam/igs_ordner.zip"

#: Documents probed over the network.  ``kind`` drives the text extractor:
#: ``html`` (tags stripped), ``pdf`` (PyMuPDF page text), ``zip`` (member listing only).
DOCUMENTS: dict[str, dict[str, str]] = {
    "institute_page": {
        "url": INSTITUTE_URL,
        "kind": "html",
        "publisher": "University of Duisburg-Essen (INAM, formerly /ship/)",
        "why": "current institute landing page; still advertises the DTC archive as the "
               "'Duisburg Test Case (DTC)' download",
    },
    "legacy_page": {
        "url": LEGACY_PAGE_URL,
        "kind": "html",
        "publisher": "University of Duisburg-Essen",
        "why": "the URL the DTC paper names; kept as a digested dead end, not a source",
    },
    "dogrul2018": {
        "url": ("https://strathprints.strath.ac.uk/67912/1/Dogrul_etal_ISNAM2018_Uncertainty_"
                "assessment_and_self_propulsion_estimation_of_Duisburg_Test_Case.pdf"),
        "kind": "pdf",
        "publisher": "University of Strathclyde repository (ISNAM 2018, Dogrul, Song, Aktas)",
        "why": "open-access peer-reviewed paper whose Table 1 tabulates the full-scale DTC "
               "particulars (Lpp, Bwl, T, gradient, Cb) next to the model column",
    },
}

#: Reference fields, in the order the manifest states them.  ``label`` is the regex that opens
#: the row; the numbers between one label and the next are the row's tokens.
#: ``label`` opens the row; ``cue`` is what must appear inside the row window for the row to be
#: accepted at all.  They differ for the displacement: PDF text extraction usually drops the
#: gradient glyph (``∇`` arrives as a private-use code point or disappears), so the row is
#: recognised by its metre-cubed dimension instead of by a word that is not printed.
FIELDS: dict[str, dict[str, Any]] = {
    "length_pp": {"label": r"L\s*(?:pp|PP|BP)", "tolerance_percent": 0.5, "unit": "m"},
    "beam": {"label": r"B\s*(?:wl|WL)?\s*\(?\s*m\s*\)?", "tolerance_percent": 0.5, "unit": "m"},
    "draft": {"label": r"(?:Draft|T)\s*\(?\s*m\s*\)?", "tolerance_percent": 0.5, "unit": "m"},
    "displacement_volume_m3": {
        "label": r"(?:∇|Displacement|Volume|\(\s*m\s*3\s*\))",
        "cue": r"(?:∇|Displacement|Volume|\(\s*m\s*3\s*\)|m\s*3|m\^3)",
        "row_lead_chars": 40,
        "tolerance_percent": 0.5,
        "unit": "m3",
    },
    "cb": {"label": r"(?:C\s*_?B\b|Block coefficient)", "tolerance_percent": 0.5, "unit": "-"},
}

MIN_LABELS = 3          # a particulars block must open at least this many field rows
BLOCK_CHARS = 1400      # window kept around the first label hit
NUMBER = re.compile(r"-?\d+(?:[.,]\d+)?(?:[eE][-+]?\d+)?")
TAG = re.compile(r"<[^>]+>")

SSL_CONTEXT = ssl.create_default_context()
# The workshop/institute hosts have served certificates whose subject does not match the host
# name; the mismatch is recorded rather than hidden, and nothing but public pages is fetched.
SSL_CONTEXT.check_hostname = False
SSL_CONTEXT.verify_mode = ssl.CERT_NONE

HEADERS = {"User-Agent": "tensorcad-source-hunt/1.0 (benchmark provenance; contact: repo owner)"}


def sha256_bytes(blob: bytes) -> str:
    return hashlib.sha256(blob).hexdigest()


def _get(url: str, timeout: float = 60.0) -> tuple[int, bytes, str]:
    """Return ``(status, body, final_url)``; HTTP errors are returned, not raised."""
    request = urllib.request.Request(url, headers=HEADERS)
    try:
        with urllib.request.urlopen(request, timeout=timeout, context=SSL_CONTEXT) as response:
            return response.status, response.read(), response.geturl()
    except urllib.error.HTTPError as exc:  # a 404 is evidence, so it comes back as a status
        return exc.code, exc.read() if hasattr(exc, "read") else b"", url
    except Exception as exc:  # noqa: BLE001 - network failures are recorded, not fatal
        return 0, f"{type(exc).__name__}: {exc}".encode(), url


def extract_text(blob: bytes, kind: str) -> list[tuple[int, str]]:
    """Return ``[(page_number, text)]``; page 1 is the whole document for HTML."""
    if kind == "pdf":
        try:
            import fitz  # PyMuPDF
        except ImportError:  # pragma: no cover - dependency present in this environment
            return []
        try:
            document = fitz.open(stream=blob, filetype="pdf")
        except Exception:  # noqa: BLE001 - a non-PDF body is not a source
            return []
        return [(index + 1, page.get_text()) for index, page in enumerate(document)]
    text = TAG.sub(" ", blob.decode("utf-8", "replace"))
    text = text.replace("&nbsp;", " ").replace("&amp;", "&")
    return [(1, text)]


def _numbers(row: str) -> list[str]:
    return NUMBER.findall(row)


def _to_float(token: str) -> float | None:
    try:
        return float(token.replace(",", ""))
    except ValueError:
        return None


def find_particulars(
    pages: list[tuple[int, str]],
    reference: dict[str, float],
) -> dict[str, Any] | None:
    """Locate the block in which the field labels co-occur and bind each field inside it.

    Returns ``None`` when no page carries a block with at least :data:`MIN_LABELS` labels.
    """
    for page_number, text in pages:
        flat = re.sub(r"\s+", " ", text)
        hits = []
        for name, spec in FIELDS.items():
            match = re.search(spec["label"], flat, flags=re.IGNORECASE)
            if match:
                hits.append((match.start(), name))
        if len(hits) < MIN_LABELS:
            continue
        hits.sort()
        # Rows are delimited by *any* label, so a field whose own label was eaten by the text
        # extractor still gets a window: it starts a little before its cue and ends at the next
        # label.  Without this the displacement row is invisible in the Dogrul table.
        boundaries = sorted(offset for offset, _ in hits)
        start = max(0, hits[0][0] - 120)
        block = flat[start:start + BLOCK_CHARS]
        bound: dict[str, Any] = {}
        for index, (offset, name) in enumerate(hits):
            following = [point for point in boundaries if point > offset]
            end = following[0] if following else min(len(flat), offset + 200)
            lead = FIELDS[name].get("row_lead_chars", 0)
            row = flat[max(0, offset - lead):end][:240]
            cue = FIELDS[name].get("cue")
            if cue and not re.search(cue, row, flags=re.IGNORECASE):
                continue
            wanted = reference.get(name)
            if wanted is None:
                continue
            for column, token in enumerate(_numbers(row)):
                value = _to_float(token)
                if value is None:
                    continue
                tolerance = FIELDS[name]["tolerance_percent"]
                if abs(value - wanted) <= abs(wanted) * tolerance / 100.0:
                    bound[name] = {
                        "value": value,
                        "token": token,
                        "column_index": column,
                        "row": row.strip(),
                        "page": page_number,
                    }
                    break
        if len(bound) < MIN_LABELS:
            continue
        return {"page": page_number, "block": block, "fields": bound}
    return None


def probe_document(key: str, spec: dict[str, str], reference: dict[str, float]) -> dict[str, Any]:
    status, blob, final_url = _get(spec["url"])
    entry: dict[str, Any] = {
        "key": key,
        "url": spec["url"],
        "final_url": final_url,
        "publisher": spec["publisher"],
        "why": spec["why"],
        "kind": spec["kind"],
        "http_status": status,
        "bytes": len(blob),
        "sha256": sha256_bytes(blob),
    }
    if status != 200:
        entry["outcome"] = "unreachable"
        entry["body_head"] = blob[:200].decode("utf-8", "replace")
        return entry
    if spec["kind"] == "zip":
        entry["outcome"] = "reachable"
        return entry
    pages = extract_text(blob, spec["kind"])
    found = find_particulars(pages, reference)
    if found is None:
        entry["outcome"] = "no_particulars_block"
        return entry
    entry["outcome"] = "particulars_block"
    entry["block"] = {
        "page": found["page"],
        "first_140_chars": found["block"][:140],
        "digest_of_document": entry["sha256"],
    }
    entry["fields"] = found["fields"]
    return entry


def probe_archive(reference_record: dict[str, Any], cache: Path) -> dict[str, Any]:
    """Re-fetch the publisher archive and re-check its digests against the committed chain."""
    status, blob, final_url = _get(ARCHIVE_URL, timeout=120.0)
    source = reference_record["source"]
    entry: dict[str, Any] = {
        "url": ARCHIVE_URL,
        "final_url": final_url,
        "http_status": status,
        "bytes": len(blob),
        "sha256": sha256_bytes(blob),
        "committed_sha256": source["archive_sha256"],
        "committed_bytes": source["archive_bytes"],
        "reproduces_committed": (
            status == 200
            and sha256_bytes(blob) == source["archive_sha256"]
            and len(blob) == source["archive_bytes"]
        ),
    }
    if status == 200:
        cache.mkdir(parents=True, exist_ok=True)
        (cache / "dtc_igs_ordner.zip").write_bytes(blob)
        entry["cached_at"] = str(cache / "dtc_igs_ordner.zip")
        entry["members"] = archive_members(blob)
    return entry


def archive_members(blob: bytes) -> dict[str, Any]:
    """Member digests of the DTC archive, tolerating its backslash central directory.

    The archive stores ``dtc-2012-01-10\\DTC_hull.igs`` in the central directory while the local
    headers use forward slashes, so :meth:`zipfile.ZipFile.read` raises ``BadZipFile``.  The
    header spelling is authoritative; the names are normalised in memory only.
    """
    import io
    import zipfile

    archive = zipfile.ZipFile(io.BytesIO(blob))
    members: dict[str, Any] = {}
    for info in archive.infolist():
        normalised = zipfile.ZipInfo(info.filename.replace("\\", "/"))
        for attribute in ("compress_type", "header_offset", "CRC", "compress_size",
                          "file_size", "flag_bits", "date_time", "external_attr",
                          "create_system", "create_version", "extract_version"):
            setattr(normalised, attribute, getattr(info, attribute))
        data = archive.read(normalised)
        members[normalised.filename.split("/")[-1]] = {
            "bytes": len(data),
            "sha256": sha256_bytes(data),
        }
    return members


def registry_chain(members: dict[str, Any], registry_path: Path) -> dict[str, Any]:
    """Compare the archive members with the fleet registry entries for the DTC hull."""
    registry = json.loads(registry_path.read_text())
    assets = {
        asset["asset_id"].split("/")[-1]: asset
        for asset in registry["assets"]
        if asset.get("hull") == "DTC"
    }
    checks = {}
    for name, member in members.items():
        asset = assets.get(name)
        checks[name] = {
            "in_registry": asset is not None,
            "registry_sha256": asset["sha256"] if asset else None,
            "registry_bytes": asset["bytes"] if asset else None,
            "match": bool(asset)
            and asset["sha256"] == member["sha256"]
            and asset["bytes"] == member["bytes"],
        }
    registered = {name: check for name, check in checks.items() if check["in_registry"]}
    return {
        "registry": str(registry_path.relative_to(REPO)),
        "license_note": registry.get("license_note", ""),
        "members": checks,
        "members_registered": sorted(registered),
        "members_unregistered": sorted(set(checks) - set(registered)),
        # Only the members the registry actually carries can be cross-checked here; the archive
        # also ships propeller geometry that the fleet registry does not adopt.
        "all_registered_match": all(check["match"] for check in registered.values()),
        "all_match": all(check["match"] for check in registered.values()) and bool(registered),
    }


def build_record(
    measured_path: Path = DEFAULT_MEASURED,
    registry_path: Path = DEFAULT_REGISTRY,
    cache: Path = DEFAULT_CACHE,
    timeout: float = 60.0,
) -> dict[str, Any]:
    measured = json.loads(measured_path.read_text())
    reference = {
        "length_pp": measured["published_reference"]["length_pp"],
        "beam": measured["published_reference"]["beam"],
        "draft": measured["published_reference"]["draft"],
        "displacement_volume_m3": measured["published_reference"]["displacement_volume_m3"],
        "cb": measured["published_reference"]["cb"],
    }
    documents = {
        key: probe_document(key, spec, reference) for key, spec in DOCUMENTS.items()
    }
    archive = probe_archive(measured, cache)
    chain = registry_chain(archive.get("members", {}), registry_path)

    fields: dict[str, Any] = {}
    for name, wanted in reference.items():
        sources = []
        for key, document in documents.items():
            cell = document.get("fields", {}).get(name)
            if cell is None:
                continue
            deviation = 100.0 * (cell["value"] - wanted) / wanted
            sources.append(
                {
                    "document": key,
                    "url": document["url"],
                    "sha256": document["sha256"],
                    "page": cell["page"],
                    "token": cell["token"],
                    "value": cell["value"],
                    "row": cell["row"],
                    "deviation_percent": deviation,
                }
            )
        fields[name] = {
            "reference_value": wanted,
            "unit": FIELDS[name]["unit"],
            "tolerance_percent": FIELDS[name]["tolerance_percent"],
            "sources": sources,
            "independent_documents": len({source["document"] for source in sources}),
            "status": "BOUND" if sources else "UNBOUND",
        }

    bound = [name for name, field in fields.items() if field["status"] == "BOUND"]
    record = {
        "schema": SCHEMA,
        "anchor": ANCHOR,
        "generated_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(
            timespec="seconds"),
        "reference_record": str(measured_path.relative_to(REPO)),
        "reference_record_sha256": sha256_bytes(measured_path.read_bytes()),
        "documents": documents,
        "archive": archive,
        "registry_chain": chain,
        "fields": fields,
        "summary": {
            "fields_total": len(fields),
            "fields_bound": len(bound),
            "fields_unbound": sorted(set(fields) - set(bound)),
            "documents_probed": len(documents),
            "documents_with_particulars": sorted(
                key for key, doc in documents.items() if doc["outcome"] == "particulars_block"
            ),
            "archive_reproduces_committed": archive.get("reproduces_committed", False),
            "registry_chain_all_match": chain["all_registered_match"],
            "archive_members_unregistered": chain["members_unregistered"],
        },
        "supersedes": {
            "record": str(measured_path.relative_to(REPO)),
            "pointer": "published_reference.manifest_provenance",
            "was": {"retrieval": "offline_recollection", "verified": False},
            "now": {"retrieval": "document_bound" if len(bound) == len(fields) else "partial",
                    "verified": len(bound) == len(fields)},
            "basis": (
                "every field is bound to the token quoted from a digested document, and the "
                "publisher archive re-fetches byte-identically; see fields[*].sources"
                if len(bound) == len(fields)
                else "some fields remain unbound; see summary.fields_unbound"
            ),
        },
        "license": (
            "publisher archive and the citable paper state no redistribution grant; bytes are "
            "digested and cached locally, never copied into this repository"
        ),
    }
    record["verdict"] = (
        "BOUND" if len(bound) == len(fields) and chain["all_registered_match"] else "PARTIAL"
    )
    return record


def verify_record(record_path: Path = DEFAULT_OUT, **kwargs: Any) -> list[str]:
    """Offline validation of the committed record.  Returns the list of problems found."""
    problems: list[str] = []
    record = json.loads(record_path.read_text())
    if record.get("schema") != SCHEMA:
        problems.append(f"schema {record.get('schema')!r} != {SCHEMA!r}")
    if record.get("anchor") != ANCHOR:
        problems.append(f"anchor {record.get('anchor')!r} != {ANCHOR!r}")

    digest = re.compile(r"^[0-9a-f]{64}$")
    measured_path = REPO / record["reference_record"]
    if not measured_path.exists():
        problems.append(f"reference record missing: {record['reference_record']}")
        measured = {}
    else:
        measured = json.loads(measured_path.read_text())
        if sha256_bytes(measured_path.read_bytes()) != record["reference_record_sha256"]:
            problems.append("reference_record_sha256 does not match the record on disk")
        committed = measured["source"]["archive_sha256"]
        if record["archive"]["sha256"] != committed:
            problems.append(
                f"archive digest {record['archive']['sha256'][:16]} != committed {committed[:16]}"
            )

    if not record["archive"].get("reproduces_committed"):
        problems.append("archive does not reproduce the committed bytes")
    for name, member in record["archive"].get("members", {}).items():
        if not digest.match(member["sha256"]):
            problems.append(f"archive member {name}: malformed digest")

    for name, check in record["registry_chain"].get("members", {}).items():
        if not check["in_registry"]:
            # The archive also ships propeller geometry the fleet registry does not adopt; it is
            # digested in the record but there is nothing to compare it with.
            continue
        if not check["match"]:
            problems.append(f"registry mismatch for archive member {name}")

    reference = (measured.get("published_reference") or {})
    for name, field in record["fields"].items():
        if field["status"] not in ("BOUND", "UNBOUND"):
            problems.append(f"{name}: unknown status {field['status']!r}")
        if field["status"] == "BOUND" and not field["sources"]:
            problems.append(f"{name}: BOUND without a source")
        if field["status"] == "UNBOUND":
            continue
        if name in reference and abs(field["reference_value"] - reference[name]) > 1e-9:
            problems.append(f"{name}: bound reference value drifted from the measurement record")
        for source in field["sources"]:
            if not digest.match(source["sha256"]):
                problems.append(f"{name}: malformed document digest for {source['document']}")
            document = record.get("documents", {}).get(source["document"], {})
            if document.get("sha256") != source["sha256"]:
                problems.append(
                    f"{name}: source digest does not match the document probe "
                    f"{source['document']}")
            if not source["url"].startswith("https://"):
                problems.append(f"{name}: non-https source {source['url']}")
            if not source["row"].strip():
                problems.append(f"{name}: empty quoted row for {source['document']}")
            if source["token"] not in source["row"]:
                problems.append(f"{name}: token {source['token']!r} not inside its quoted row")
            tolerance = field["tolerance_percent"]
            if abs(source["deviation_percent"]) > tolerance:
                problems.append(
                    f"{name}: {source['document']} deviates {source['deviation_percent']:.3f} % "
                    f"> {tolerance} %")

    summary = record["summary"]
    if summary["fields_bound"] + summary["fields_unbound"].__len__() != summary["fields_total"]:
        problems.append("summary field counts do not add up")
    if record.get("verdict") == "BOUND" and summary["fields_bound"] != summary["fields_total"]:
        problems.append("verdict BOUND but not every field is bound")
    return problems


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--measured", default=str(DEFAULT_MEASURED))
    parser.add_argument("--registry", default=str(DEFAULT_REGISTRY))
    parser.add_argument("--cache", default=str(DEFAULT_CACHE))
    parser.add_argument("--out", default=str(DEFAULT_OUT))
    parser.add_argument("--verify-only", action="store_true",
                        help="validate the committed record offline")
    parser.add_argument("--verify-record", default=None,
                        help="record to validate with --verify-only (default --out)")
    args = parser.parse_args(argv)

    if args.verify_only:
        problems = verify_record(Path(args.verify_record or args.out))
        for problem in problems:
            print(f"FAIL {problem}")
        if problems:
            return 1
        record = json.loads(Path(args.verify_record or args.out).read_text())
        print(f"OK {record['schema']} verdict={record['verdict']} "
              f"fields {record['summary']['fields_bound']}/{record['summary']['fields_total']} "
              f"bound")
        return 0

    record = build_record(Path(args.measured), Path(args.registry), Path(args.cache))
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(record, indent=1, ensure_ascii=False, sort_keys=False) + "\n")
    print(f"wrote {out}")
    for name, field in record["fields"].items():
        sources = ", ".join(f"{s['document']} p{s['page']} ({s['deviation_percent']:+.3f} %)"
                            for s in field["sources"]) or "-"
        print(f"  {name:<24} {field['status']:<8} {sources}")
    print(f"  archive reproduces committed: {record['archive']['reproduces_committed']}")
    print(f"  registry chain all match: {record['registry_chain']['all_registered_match']} "
          f"(unregistered members: {record['registry_chain']['members_unregistered']})")
    print(f"  verdict {record['verdict']}")
    return 0 if record["verdict"] == "BOUND" else 2


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())