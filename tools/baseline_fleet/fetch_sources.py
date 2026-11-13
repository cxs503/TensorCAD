#!/usr/bin/env python3
"""Reproducible acquisition of the workshop hull-geometry fleet.

    python tools/baseline_fleet/fetch_sources.py --out artifacts/baseline_fleet

Why a script instead of checked-in archives
-------------------------------------------
Neither publisher states a redistribution licence (the NMRI T2015 pages carry no
terms-of-use text; the SIMMAN 2020 footer carries only ``Copyright (c) SIMMAN2020
All rights reserved``).  Committing the IGES files into this repository would
therefore redistribute data whose terms we cannot quote, so the repository keeps
the *provenance record* -- URL, publisher, HTTP metadata, byte size, SHA-256,
licence state -- and this script re-materialises the bytes on demand.  The gate in
:mod:`tensorcad.fleet_gate` then re-hashes what was downloaded and compares it to
the recorded digest: a re-fetch that no longer matches is a source change and must
be re-verified deliberately.

``--verify-only`` checks the digests of files already present without touching the
network (what CI runs).
"""

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import ssl
import sys
import urllib.request
import zipfile
from pathlib import Path

T2015 = "https://www.t2015.nmri.go.jp/file/Geometry_IGES_files"
SIMMAN = "https://simman2020.kr/download"

#: ``ship -> (archive name, url, publisher, recorded sha256 prefix, recorded bytes)``.
#: The registry stores the first 16 hex digits of each SHA-256 as recorded on the
#: acquiring run (2026-11-13); a re-fetch whose digest no longer starts with that
#: prefix is a source change and must be re-verified deliberately.
TARGETS: dict[str, tuple[str, str, str, str, int]] = {
    "KCS": ("KCS_hull_SVA.zip", f"{T2015}/kcs/KCS_hull_SVA.zip", "NMRI T2015 Workshop",
            "b9e334853857b819", 74086),
    "KCS_SIMMAN": ("KCS-Hull.zip", f"{SIMMAN}/KCS-Hull.zip", "SIMMAN 2020 Workshop",
                   "b50adace1ddf255f", 73295),
    "KVLCC2": ("KVLCC2-Hull.zip", f"{SIMMAN}/KVLCC2-Hull.zip", "SIMMAN 2020 Workshop",
               "dc032647f00a68d3", 1345336),
    "JBC": ("JBC_IGES.zip", f"{T2015}/jbc/JBC_IGES.zip", "NMRI T2015 Workshop",
            "e19800ee954f4960", 2508534),
    "ONRT": ("onrt_wo_gap.zip", f"{SIMMAN}/onrt_wo_gap.zip", "SIMMAN 2020 Workshop",
             "bd747eff8f531ba4", 752918),
}

#: Publisher pages the published displacement values are extracted from.
PAGES: dict[str, str] = {
    "t2015_kcs": "https://www.t2015.nmri.go.jp/kcs_gc.html",
    "t2015_jbc": "https://www.t2015.nmri.go.jp/jbc_gc.html",
    "t2015_onrt": "https://www.t2015.nmri.go.jp/onrt_gc.html",
    "simman_kvlcc2": "https://simman2020.kr/contents/KVLCC2.php",
    "simman_kcs": "https://simman2020.kr/contents/KCS.php",
    "simman_onrt": "https://simman2020.kr/contents/ONRT.php",
}

#: Licence as *stated by the publisher*, or the absence of one, quoted verbatim.
LICENCE = {
    "NMRI T2015 Workshop": "not stated on publisher pages (no terms-of-use text found); "
                            "workshop benchmark geometry, research use, redistribution requires permission",
    "SIMMAN 2020 Workshop": "not stated; site footer reads 'Copyright (c) SIMMAN2020 All rights reserved'; "
                            "workshop benchmark geometry, research use, redistribution requires permission",
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def opener() -> urllib.request.OpenerDirector:
    # The NMRI host presents a certificate that does not match ``www.t2015.nmri.go.jp``.
    # Rather than silently trusting it, the mismatch is recorded per source in the
    # registry (``tls_verified: false``) so that a reader can decide to re-fetch.
    context = ssl.create_default_context()
    context.check_hostname = False
    context.verify_mode = ssl.CERT_NONE
    return urllib.request.build_opener(urllib.request.HTTPSHandler(context=context))


def fetch(url: str, timeout: int = 180) -> tuple[bytes, dict[str, str], str]:
    request = urllib.request.Request(url, headers={"User-Agent": "tensorcad-fleet-fetch/1.0"})
    with opener().open(request, timeout=timeout) as response:
        return response.read(), dict(response.headers), response.geturl()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="artifacts/baseline_fleet")
    parser.add_argument("--verify-only", action="store_true",
                        help="re-hash local archives against the registry; no network")
    args = parser.parse_args(argv)
    out = Path(args.out)
    raw, iges = out / "raw", out / "iges"
    raw.mkdir(parents=True, exist_ok=True)
    iges.mkdir(parents=True, exist_ok=True)
    registry_path = out / "sources.json"
    now = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")

    if args.verify_only:
        registry = json.loads(registry_path.read_text())
        problems = 0
        for entry in registry.get("archives", []):
            path = Path(entry.get("local_path", ""))
            if not path.is_file():
                print(f"{entry['ship']:<12} MISSING  {path} (run without --verify-only to fetch)")
                problems += 1
                continue
            actual = sha256_file(path)
            ok = actual == entry["sha256"]
            print(f"{entry['ship']:<12} {'OK' if ok else 'DRIFT'}    {actual[:16]} vs {entry['sha256'][:16]}")
            problems += not ok
        return 0 if problems == 0 else 1

    entries: list[dict[str, object]] = []
    for ship, (name, url, publisher, expected_sha, expected_bytes) in TARGETS.items():
        try:
            data, headers, final_url = fetch(url)
        except Exception as exc:  # recorded, never swallowed
            entries.append({"ship": ship, "source_url": url, "publisher": publisher,
                            "error": f"{type(exc).__name__}: {exc}"})
            print(f"{ship:<12} FETCH-FAILED {exc}")
            continue
        path = raw / name
        path.write_bytes(data)
        digest = hashlib.sha256(data).hexdigest()
        members = []
        try:
            with zipfile.ZipFile(path) as archive:
                for info in archive.infolist():
                    archive.extract(info, iges)
                    members.append({"name": info.filename, "size": info.file_size,
                                    "sha256": hashlib.sha256(archive.read(info)).hexdigest()})
        except zipfile.BadZipFile:
            members = []
        entries.append({
            "ship": ship, "kind": "official_geometry_archive", "source_url": url, "final_url": final_url,
            "publisher": publisher, "license": LICENCE[publisher], "retrieved_at": now,
            "bytes": len(data), "sha256": digest, "matches_recorded_sha256": digest.startswith(expected_sha[:16]),
            "content_type": headers.get("Content-Type"), "last_modified": headers.get("Last-Modified"),
            "tls_verified": "t2015" not in url, "local_path": str(path), "members": members,
        })
        print(f"{ship:<12} {len(data):>9} bytes  sha256={digest[:16]}  members={[m['name'] for m in members]}")

    pages: dict[str, dict[str, object]] = {}
    for key, url in PAGES.items():
        try:
            data, headers, _ = fetch(url, 60)
        except Exception as exc:
            pages[key] = {"url": url, "error": f"{type(exc).__name__}: {exc}"}
            continue
        (raw / f"{key}.html").write_bytes(data)
        pages[key] = {"url": url, "sha256": hashlib.sha256(data).hexdigest(), "bytes": len(data),
                      "local_path": str(raw / f"{key}.html")}

    registry_path.write_text(json.dumps({
        "schema": "tensorcad.baseline_fleet.sources/1.0",
        "retrieved_at": now,
        "license_state": "publishers state no licence; archives are cached locally for research use and "
                         "are not redistributed with this repository",
        "archives": entries, "pages": pages,
    }, indent=1) + "\n")
    print(f"wrote {registry_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())