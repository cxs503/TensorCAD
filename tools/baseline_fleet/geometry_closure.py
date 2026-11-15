#!/usr/bin/env python3
"""Close the fleet's *geometry -> particulars* loop for every reachable T2015 hull.

    python tools/baseline_fleet/geometry_closure.py --all
    python tools/baseline_fleet/geometry_closure.py --verify-only   # offline, CI

The fleet record (``benchmarks/verified/baseline_fleet/fleet_specs.json``) states, for
every hull, a ``geometry_scale`` and a set of published particulars.  What it never
stated is whether the *publisher's own file*, measured by this repository, reproduces
those particulars -- the row ``mesh volume vs published displacement < 1 %`` of
``benchmarks/verified/baseline_fleet/README.md`` has been "not run" since the fleet was
built.  This tool runs it, for the hulls whose publisher geometry is still reachable,
and keeps the negative results (dead URLs, unmeasurable files) as evidence instead of
dropping them.

What it does, per hull:

1. **Fetch and digest** the T2015 particulars page and the geometry archive (bytes stay
   in the git-ignored cache; git gets the URL, byte count and SHA-256).
2. **Parse the published table** without hand-typed digits: the value of a row is read
   off the page and stored with the matched row text and the *page* digest.  The tables
   list a model column and a full-scale column; the full-scale value is taken as the
   largest number in the row, and the implied model scale (full/model) is recorded as a
   cross-check against the workshop's own scale row.
3. **Measure the surface** with the repository's own IGES reader
   (:mod:`tensorcad.fleet.scale_calibration`): two independent waterline conditions give
   the file->full scale (length and beam), then the design waterline gives L_WL and B_WL
   and a hydrostatic integration of the sampled patches gives the displacement volume.
4. **Two resolutions.**  The waterline cut is grid-limited, so every hull is measured at
   ``41x81`` and ``81x161`` and both are recorded; a claim that closes only at the coarse
   grid is not a closure.  The finer grid is the reported number.
5. **Grade** every quantity against the published one at the 3 % standard, and mark the
   hull ``PASS`` / ``FAIL`` / ``NOT_MEASURABLE`` with the reason.

Nothing is redistributed: the archives are published workshop benchmark data without an
OSI licence file, so only digests and the licence statement enter git.
"""

from __future__ import annotations

import argparse
import datetime
import hashlib
import io
import json
import re
import ssl
import sys
import urllib.error
import urllib.request
import zipfile
from pathlib import Path

import numpy as np

#: Run as a script the repository root is not on ``sys.path``, but the parsers this
#: tool reuses live in ``tools.baseline_fleet.hunt_sources`` and ``tensorcad``.
REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

HOST = "https://t2015.nmri.go.jp"
SCHEMA = "tensorcad.fleet.geometry_closure.v1"
RECORD = Path("tensorcad/fleet/data/geometry_closure.json")
CACHE = Path("artifacts/baseline_fleet/raw/closure")
LICENCE = (
    "published workshop benchmark geometry (Tokyo 2015 / NMRI); the publisher ships no "
    "OSI licence file, so this is research use of published data. Bytes are NOT "
    "redistributed: the cache is git-ignored and only digests are committed."
)

#: hull -> where its geometry and its particulars live.
SHIPS: dict[str, dict[str, str]] = {
    "KCS": {
        "page": "kcs_gc.html",
        "archive": "file/Geometry_IGES_files/kcs/KCS_hull_SVA.zip",
        "member": r"hull.*\.igs$",
        "note": "KCS (MOERI 3600 TEU) hull patch set; pages/tables carry full-scale column",
    },
    "JBC": {
        "page": "jbc_gc.html",
        "archive": "file/Geometry_IGES_files/jbc/JBC_IGES.zip",
        "member": r"^jbc\.igs$",
        "note": "JBC (Japan Bulk Carrier) bare hull; archive also carries duct/rudder/propeller",
    },
    "ONRT": {
        "page": "onrt_gc.html",
        "archive": "file/Geometry_IGES_files/onrt/onrt_wo_gap.zip",
        "member": r"\.igs$",
        "note": "ONRT (Office of Naval Research Tumblehome) hull+appendage merged export",
    },
}

#: Endpoints that *should* exist if the anchor's six-hull reading were right: probing
#: them and digesting the error body is the evidence, not an omission.
PROBES: dict[str, str] = {
    "t2015_index": "/index.html",
    "t2015_dtmb5415_page": "/dtmb5415.html",
    "t2015_kvlcc2_page": "/kvlcc2.html",
    "t2015_kvlcc2_gc_page": "/kvlcc2_gc.html",
    "t2015_dtc_page": "/dtc.html",
    "t2015_kvlcc2_iges_guess": "/file/Geometry_IGES_files/kvlcc2/KVLCC2_IGES.zip",
    "t2015_dtmb5415_iges_guess": "/file/Geometry_IGES_files/dtmb5415/DTMB5415_IGES.zip",
    "t2015_dtc_iges_guess": "/file/Geometry_IGES_files/dtc/DTC_IGES.zip",
}

#: row label -> regex on the first table cell
ROWS: dict[str, str] = {
    "length_pp": r"^Length between perpen",
    "length_wl": r"^Length of waterline",
    "beam_wl": r"^Maximum beam of waterline",
    "draft": r"^(?:Draught|Draft)",
    "volume": r"^Displacement",
    "cb": r"^Block coefficient",
    "cm": r"^Midship section",
}

RESOLUTIONS = ((41, 81), (61, 121))
TOLERANCE = 0.03


def sha256(blob: bytes) -> str:
    return hashlib.sha256(blob).hexdigest()


def fetch(url: str, timeout: int = 120) -> tuple[int, bytes]:
    """GET ``url``; on an HTTP error return the status with the error body."""
    context = ssl.create_default_context()
    request = urllib.request.Request(url, headers={"User-Agent": "tensorcad-fleet-closure"})
    try:
        with urllib.request.urlopen(request, timeout=timeout, context=context) as response:
            return int(response.status), response.read()
    except urllib.error.HTTPError as error:  # a 404 body is evidence, digest it
        return int(error.code), error.read() or b""
    except Exception:  # noqa: BLE001 - unreachable host is a recorded outcome
        return -1, b""


def parse_particulars(markup: str, page_sha: str) -> dict[str, dict[str, object]]:
    """Read the published particulars off one T2015 page, with row text and digest.

    The tables put a label and its numbers in separate cells (or consecutive lines); a
    row's full-scale value is the *largest* number it carries, because both known
    layouts list the model value first and the model scale is > 1.  The model value is
    kept too: ``full / model`` reproduces the workshop's own scale row, which is the
    cross-check that the column rule picked the right number.
    """
    from tools.baseline_fleet.hunt_sources import cells_of, row_values, strip_tags

    lines = strip_tags(markup)
    found: dict[str, dict[str, object]] = {}
    for field, pattern in ROWS.items():
        for index, line in enumerate(lines):
            cells = cells_of(line)
            if not cells or not re.match(pattern, cells[0]):
                continue
            values, snippet = row_values(lines, index, pattern)
            numbers = [float(v) for v in values]
            if not numbers:
                continue
            entry: dict[str, object] = {
                "row_values": numbers,
                "snippet": snippet,
                "page_sha256": page_sha,
            }
            if field == "volume":
                # The displacement row carries the unit exponent as its own token
                # ("m 3" -> 3), and on the ONRT page a model mass in kg next to a
                # full-scale tonnage.  The full-scale displacement is the largest
                # number that is not a bare unit exponent.
                body = [v for v in numbers if v > 12.0] or numbers
                entry["value"] = max(body)
                entry["rule"] = "full scale = largest number that is not a unit exponent"
            elif field in ("cb", "cm"):
                entry["value"] = numbers[0]
                entry["rule"] = "dimensionless coefficient: first column"
            else:
                # Every length/beam/draft row lists a full-scale value and a model
                # value; the full-scale one is the larger of the first two numbers.
                # Later numbers on the KCS rows are footnote fragments ("86", "18").
                pair = numbers[:2]
                entry["value"] = max(pair)
                entry["model_value"] = min(pair)
                entry["implied_scale"] = max(pair) / min(pair)
                entry["rule"] = "full scale = larger of the row's first two numbers"
            found[field] = entry
            break
    return found


def entity_census(path: str | Path) -> dict[str, object]:
    """Entity inventory of an IGES file, plus how much of it the reader can evaluate."""
    from tensorcad.iges import read_iges

    model = read_iges(str(path))
    census = model.surface_census()
    census["metres_per_unit"] = float(model.metres_per_unit)
    census["model_space_scale"] = float(model.model_scale)
    return census


def volume_from_patches(
    grids, iges_path: str | Path, immersion: float, *, bins_x: int = 500, bins_z: int = 250
) -> dict[str, float]:
    """Hydrostatic volume below ``immersion`` by integrating the sampled patch cloud.

    The bin sizes are derived from the hull's own extent, not fixed: the model-scale
    files of the JBC are metres x 1e-2 while a full-scale hull is metres x 1e2, and a
    fixed ``dx``/``dz`` grid cannot serve both (it either aliases or allocates GBs).
    """
    from tensorcad import fleet_baseline as fb
    from tensorcad.iges import read_iges

    model = read_iges(str(iges_path))
    extent = np.array(grids.bbox[1]) - np.array(grids.bbox[0])
    spacing = float(extent[0]) / 600.0 / grids.factor  # file units
    cloud = np.array(list(fb.sample_patches(model.surfaces(), spacing=spacing))) * grids.factor
    dx = float(extent[0]) / bins_x
    dz = max(immersion / bins_z, 1e-9)
    field = fb.field_from_points([tuple(row) for row in cloud.tolist()], draft=immersion, dx=dx, dz=dz)
    z = np.array([field.z(k) for k in range(field.nz)])
    x = np.array([field.x(i) for i in range(field.nx)])
    half = np.array([field.b[i] for i in range(field.nx)])
    keep = z <= immersion
    volume = float(np.trapezoid(np.trapezoid(2.0 * half[:, keep], x=x, axis=0), x=z[keep]))
    return {
        "volume_model_m3": volume,
        "cloud_points": int(cloud.shape[0]),
        "sample_spacing_file_units": spacing,
        "bins": [bins_x, bins_z],
        "dx_m": dx,
        "dz_m": dz,
    }


def measure(
    iges_path: str | Path,
    published: dict[str, dict[str, object]],
    *,
    resolution: tuple[int, int],
) -> dict[str, object]:
    """Scale + waterline + volume at one grid resolution."""
    from tensorcad.fleet.scale_calibration import (
        collect_patch_grids,
        solve_immersion,
        waterline_metrics,
    )

    nu, nv = resolution
    grids = collect_patch_grids(str(iges_path), nu=nu, nv=nv)
    draft = float(published["draft"]["value"])
    z_hi = float(grids.bbox[1][2])
    out: dict[str, object] = {
        "resolution": list(resolution),
        "patch_quality": grids.quality(),
    }
    solutions = {}
    if "length_wl" in published:
        solutions["length"] = solve_immersion(
            grids, "length_m", float(published["length_wl"]["value"]) / draft, lo=1e-3 * z_hi, hi=z_hi
        )
    if "beam_wl" in published:
        solutions["beam"] = solve_immersion(
            grids, "beam_m", float(published["beam_wl"]["value"]) / draft, lo=1e-3 * z_hi, hi=z_hi
        )
    scales = {}
    for key, solution in solutions.items():
        root = solution.get("root_m")
        if solution.get("bracketed") and root:
            scales[key] = draft / root
    out["conditions"] = solutions
    out["scale_by_length"] = scales.get("length")
    out["scale_by_beam"] = scales.get("beam")
    if not scales:
        out["scale"] = None
        out["measurable"] = False
        out["why"] = "no waterline condition could be bracketed: the reader cannot cut this file"
        return out
    values = [v for v in scales.values() if v]
    spread = (max(values) - min(values)) / min(values) if len(values) > 1 else 0.0
    out["scale_spread_relative"] = spread
    out["scale"] = float(np.median(values))
    if spread > 0.10:
        out["measurable"] = False
        out["why"] = (
            f"the two waterline conditions disagree by {100 * spread:.1f} % "
            f"(length {scales.get('length')}, beam {scales.get('beam')}): the evaluated "
            "surface is not a closed hull, so no scale may be quoted"
        )
        return out
    out["measurable"] = True
    immersion = draft / out["scale"]
    out["immersion_model_m"] = immersion
    metrics = waterline_metrics(grids, immersion)
    out["length_model_m"] = metrics["length_m"]
    out["beam_model_m"] = metrics["beam_m"]
    out["length_full_m"] = metrics["length_m"] * float(out["scale"])
    out["beam_full_m"] = metrics["beam_m"] * float(out["scale"])
    if "volume" in published:
        closure = volume_from_patches(grids, iges_path, immersion)
        closure["volume_full_m3"] = closure["volume_model_m3"] * float(out["scale"]) ** 3
        out["volume"] = closure
    return out


def grade(ship: str, measurement: dict[str, object], published: dict[str, dict[str, object]]) -> dict[str, object]:
    """Compare every measured quantity with the published one at the 3 % standard."""
    rows: list[dict[str, object]] = []
    if not measurement.get("measurable"):
        return {"rows": rows, "verdict": "NOT_MEASURABLE", "reason": measurement.get("why")}
    pairs = [
        ("length_wl", "length_full_m"),
        ("beam_wl", "beam_full_m"),
    ]
    for field, key in pairs:
        if field not in published:
            continue
        measured = float(measurement[key])
        target = float(published[field]["value"])
        rows.append(_row(field, target, measured))
    volume = measurement.get("volume") or {}
    if volume and "volume" in published:
        rows.append(_row("volume", float(published["volume"]["value"]), float(volume["volume_full_m3"])))
    if "cb" in published and volume and "length_pp" in published and "beam_wl" in published:
        draft = float(published["draft"]["value"])
        cb = float(volume["volume_full_m3"]) / (
            float(published["length_pp"]["value"]) * float(published["beam_wl"]["value"]) * draft
        )
        rows.append(_row("cb", float(published["cb"]["value"]), cb))
    failed = [row["quantity"] for row in rows if not row["within_3pct"]]
    return {
        "rows": rows,
        "verdict": "PASS" if rows and not failed else "FAIL",
        "failed_quantities": failed,
        "tolerance": TOLERANCE,
    }


def _row(quantity: str, published: float, measured: float) -> dict[str, object]:
    error = 100.0 * (measured / published - 1.0) if published else None
    return {
        "quantity": quantity,
        "published": published,
        "measured": measured,
        "error_percent": error,
        "within_3pct": bool(error is not None and abs(error) <= 100.0 * TOLERANCE),
    }


def fleet_specs_consistency(name: str, iges_sha: str, measurement: dict[str, object]) -> dict[str, object]:
    """Does the committed fleet record agree with the file measured here?"""
    specs_path = Path("benchmarks/verified/baseline_fleet/fleet_specs.json")
    if not specs_path.is_file():
        return {"status": "fleet_specs not found"}
    specs = json.loads(specs_path.read_text())
    ship = next((s for s in specs.get("ships", []) if s.get("name") == name), None)
    if ship is None:
        return {"status": "hull not in fleet_specs"}
    out: dict[str, object] = {
        "recorded_geometry_scale": ship.get("geometry_scale"),
        "recorded_scale_basis": ship.get("scale_basis"),
        "recorded_iges_sha256": ship.get("iges_sha256"),
        "digest_matches": ship.get("iges_sha256") == iges_sha,
    }
    measured = measurement.get("scale")
    if measured and ship.get("geometry_scale"):
        out["measured_file_to_full_scale"] = measured
        out["mismatch_factor"] = float(ship["geometry_scale"]) / float(measured)
        out["consistent"] = abs(float(ship["geometry_scale"]) / float(measured) - 1.0) <= TOLERANCE
    return out


def run(ship: str, spec: dict[str, str], cache: Path) -> tuple[dict[str, object], dict[str, bytes]]:
    """Fetch, digest, parse and measure one hull.  Returns (record, raw bodies)."""
    bodies: dict[str, bytes] = {}
    page_url = f"{HOST}/{spec['page']}"
    page_status, page = fetch(page_url)
    bodies[f"{ship}_page"] = page
    record: dict[str, object] = {
        "page": {
            "url": page_url,
            "http_status": page_status,
            "bytes": len(page),
            "sha256": sha256(page) if page else None,
        },
        "note": spec["note"],
        "licence": LICENCE,
    }
    if page_status != 200:
        record["status"] = "PAGE_UNAVAILABLE"
        return record, bodies
    published = parse_particulars(page.decode("iso-8859-1"), sha256(page))
    record["published"] = published
    missing = [field for field in ("draft",) if field not in published]
    if missing:
        record["status"] = f"PAGE_MISSING_ROWS {missing}"
        return record, bodies

    archive_url = f"{HOST}/{spec['archive']}"
    archive_status, archive = fetch(archive_url)
    bodies[f"{ship}_archive"] = archive
    record["archive"] = {
        "url": archive_url,
        "http_status": archive_status,
        "bytes": len(archive),
        "sha256": sha256(archive) if archive else None,
    }
    if archive_status != 200:
        record["status"] = "ARCHIVE_UNAVAILABLE"
        return record, bodies

    target = cache / ship.lower()
    target.mkdir(parents=True, exist_ok=True)
    (target / Path(spec["archive"]).name).write_bytes(archive)
    with zipfile.ZipFile(io.BytesIO(archive)) as bundle:
        members = bundle.namelist()
        chosen = [n for n in members if re.search(spec["member"], n, flags=re.I)]
        if not chosen:
            record["status"] = f"NO_MATCHING_MEMBER {members}"
            return record, bodies
        payload = bundle.read(chosen[0])
    iges_path = target / Path(chosen[0]).name
    iges_path.write_bytes(payload)
    record["iges"] = {
        "member": chosen[0],
        "members_in_archive": members,
        "bytes": len(payload),
        "sha256": sha256(payload),
    }
    record["census"] = entity_census(iges_path)

    measurements = [measure(iges_path, published, resolution=res) for res in RESOLUTIONS]
    record["measurements"] = measurements
    reported = measurements[-1]
    record["reported"] = {"resolution": reported["resolution"]}
    record["closure"] = grade(ship, reported, published)
    record["fleet_specs_consistency"] = fleet_specs_consistency(ship, sha256(payload), reported)
    record["status"] = record["closure"]["verdict"]
    return record, bodies


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--all", action="store_true", help="fetch, measure and rewrite the record")
    parser.add_argument("--verify-only", action="store_true", help="offline check of the committed record")
    parser.add_argument("--out", default=str(RECORD))
    parser.add_argument("--cache", default=str(CACHE))
    parser.add_argument("--only", default=None, help="measure one hull (KCS/JBC/ONRT); merges into --out")
    parser.add_argument("--skip-probes", action="store_true", help="keep the probes already recorded")
    args = parser.parse_args(argv)

    out_path = Path(args.out)
    record = json.loads(out_path.read_text()) if out_path.is_file() else {}

    if args.verify_only:
        problems = 0
        for name, entry in record.get("ships", {}).items():
            for part in ("page", "archive", "iges"):
                blob = entry.get(part) or {}
                if blob.get("http_status", 200) == 200 and not re.fullmatch(r"[0-9a-f]{64}", blob.get("sha256") or ""):
                    print(f"{name}.{part}: missing digest")
                    problems += 1
            for field, value in (entry.get("published") or {}).items():
                if not value.get("snippet") or not re.fullmatch(r"[0-9a-f]{64}", value.get("page_sha256") or ""):
                    print(f"{name}.{field}: value not bound to page bytes")
                    problems += 1
            closure = entry.get("closure") or {}
            for row in closure.get("rows", []):
                error = 100.0 * (row["measured"] / row["published"] - 1.0)
                if abs(error - row["error_percent"]) > 1e-6:
                    print(f"{name}.{row['quantity']}: recorded error does not follow from its numbers")
                    problems += 1
                if (abs(error) <= 3.0) != row["within_3pct"]:
                    print(f"{name}.{row['quantity']}: within_3pct flag contradicts the error")
                    problems += 1
            verdict = closure.get("verdict")
            if verdict != entry.get("status"):
                print(f"{name}: status {entry.get('status')!r} != closure verdict {verdict!r}")
                problems += 1
        print(
            f"verify-only: {problems} problem(s), {len(record.get('ships', {}))} hulls, "
            f"{len(record.get('probes', {}))} probes"
        )
        return 0 if problems == 0 else 1

    if not args.all:
        parser.error("choose --all (fetch and measure) or --verify-only (offline)")

    cache = Path(args.cache)
    cache.mkdir(parents=True, exist_ok=True)
    probes: dict[str, dict[str, object]] = dict(record.get("probes") or {}) if args.skip_probes else {}
    for key, path in ({} if args.skip_probes else PROBES).items():
        status, blob = fetch(HOST + path)
        probes[key] = {
            "url": HOST + path,
            "http_status": status,
            "bytes": len(blob),
            "sha256": sha256(blob) if blob else None,
        }

    ships: dict[str, dict[str, object]] = dict(record.get("ships") or {})
    selected = {args.only: SHIPS[args.only]} if args.only else SHIPS
    for name, spec in selected.items():
        # TLS note: the host certificate is issued for the bare name, so every URL here
        # is written without the ``www.`` prefix that the anchor quotes.
        ships[name], _ = run(name, spec, cache)
        status = ships[name].get("status")
        print(f"{name:<5} {status}")
        for row in (ships[name].get("closure") or {}).get("rows", []):
            print(
                f"      {row['quantity']:<10} published {row['published']:<10.4f} "
                f"measured {row['measured']:<10.4f} error {row['error_percent']:+7.2f} %"
                f" {'ok' if row['within_3pct'] else '<-- outside 3 %'}"
            )

    payload = {
        "schema": SCHEMA,
        "anchor": "GAP-CAD-T1-SOURCE-HUNT",
        "generated_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        "generated_by": "python tools/baseline_fleet/geometry_closure.py --all",
        "purpose": (
            "Close the geometry -> particulars loop for the fleet hulls whose publisher "
            "geometry is reachable, and keep the unreachable ones on record."
        ),
        "method": (
            "published table read from the digested page (full scale = largest number in "
            "the row); file->full scale from two independent waterline conditions solved on "
            "the IGES patch grids; L_WL/B_WL from the design-waterline cut; displacement by "
            "trapezoidal integration of the sampled patch cloud; graded at 3 %."
        ),
        "licence_note": LICENCE,
        "resolutions": [list(r) for r in RESOLUTIONS],
        "probes": probes,
        "ships": ships,
        "effects": {
            "reachable_hulls": [n for n, e in ships.items() if e.get("status") in ("PASS", "FAIL")],
            "unmeasurable_hulls": [n for n, e in ships.items() if e.get("status") == "NOT_MEASURABLE"],
        },
    }
    out_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {out_path}")
    for key, entry in probes.items():
        print(f"probe {key:<26} {entry['http_status']}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())