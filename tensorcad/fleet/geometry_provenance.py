"""Provenance of the standard-hull *geometry bytes* (anchor ``GAP-CAD-T1-SOURCE-HUNT``).

The fleet manifest (``benchmarks/verified/baseline_fleet/fleet_specs.json``) already
binds every published particular to the SHA-256 of the page it was read from.  What it
could not say is what the *geometry archives* themselves are: at which scale they were
authored, in which unit, and whether their measured dimensions reproduce the publisher's
particulars.  That question is answered here, from the bytes, with no network and no
hand-typed digit:

* :func:`iges_fingerprint` reads the IGES global section (unit flag, model-space scale,
  header maximum-coordinate field, exporting system, generation stamp);
* :func:`iges_extent_m` samples every B-spline surface and returns the bounding extent
  in metres -- sampling, not the control net, because control poles overshoot the
  surface at the ends of a pole-heavy patch (measured: 0.5 % on a hull's beam, 7.7 % on
  its length, see ``docs/gap_cad_t1_geometry_ledger.md``);
* :func:`axis_denominators` turns that extent plus the published full-scale dimensions
  into a per-axis scale denominator, so a single blundered axis cannot masquerade as a
  scale factor;
* :func:`classify_regime` says whether the archive is authored *at* the published model
  scale denominator or *normalised to unit length* -- two archives can both be
  "publisher geometry" and still need different treatment before they are meshed.

Everything is pure reading: no file is written outside the caller's ``--out``.
"""

from __future__ import annotations

import hashlib
import io
import math
import re
import zipfile
from pathlib import Path
from typing import Any, Iterable, Sequence

SCHEMA = "tensorcad.fleet.geometry_provenance.v1"

#: IGES global-section field indices (spec §2.2.3); the reader hands us the raw tuple.
_IGES_UNIT_FLAG = 13
_IGES_UNIT_NAME = 14
_IGES_GEN_STAMP = 17
_IGES_MIN_RESOLUTION = 18
_IGES_MAX_COORD = 19
_IGES_AUTHOR = 20
_IGES_ORGANISATION = 21

def _unit_table() -> dict[int, float]:
    """Reuse the reader's own unit table as the single source of truth."""
    from tensorcad.iges import UNIT_SCALE

    return dict(UNIT_SCALE)

#: A check is only reported as passing when the relative error is below this.
TOLERANCE_REL = 0.03


# --------------------------------------------------------------------------- digests
def sha256_bytes(payload: bytes) -> str:
    """SHA-256 of a byte string, lowercase hex (the only digest form used here)."""
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: str | Path, chunk: int = 1 << 20) -> str:
    """SHA-256 of a file, read in chunks so a 4 MB archive does not sit in memory twice."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(chunk), b""):
            digest.update(block)
    return digest.hexdigest()


def kind_of(name: str, payload: bytes) -> str:
    """Classify a cached artefact from its name and first bytes (html/zip/iges/pdf/other)."""
    lowered = name.lower()
    if payload[:2] == b"PK":
        return "zip"
    if payload[:4] == b"%PDF":
        return "pdf"
    head = payload[:4096].lstrip().lower()
    if head.startswith(b"<!doctype html") or head.startswith(b"<html") or b"<html" in head:
        return "html"
    if lowered.endswith(".igs") or lowered.endswith(".iges"):
        return "iges"
    return "other"


# --------------------------------------------------------------------------- archives
def zip_members(payload: bytes) -> list[dict[str, Any]]:
    """Census of a zip archive: one row per member with its own digest.

    A publisher archive is the unit of provenance, but a consumer needs the member, so
    the member digest is recorded as well -- that is what lets a future run say "the
    hull file changed, the propeller file did not".
    """
    rows: list[dict[str, Any]] = []
    with zipfile.ZipFile(io.BytesIO(payload)) as archive:
        for info in archive.infolist():
            if info.is_dir():
                continue
            data = archive.read(info.filename)
            rows.append(
                {
                    "name": info.filename,
                    "bytes": len(data),
                    "sha256": sha256_bytes(data),
                    "kind": kind_of(info.filename, data),
                }
            )
    return rows


def zip_member(payload: bytes, suffix: str) -> tuple[str, bytes]:
    """Return ``(name, bytes)`` of the first member whose name ends with ``suffix``."""
    with zipfile.ZipFile(io.BytesIO(payload)) as archive:
        for info in archive.infolist():
            if info.is_dir() or not info.filename.lower().endswith(suffix.lower()):
                continue
            return info.filename, archive.read(info.filename)
    raise KeyError(f"no member ending in {suffix!r}")


# --------------------------------------------------------------------------- IGES
def _float(field: str | None) -> float | None:
    if field is None:
        return None
    text = field.strip()
    if not text:
        return None
    try:
        return float(text.replace("D", "E").replace("d", "e"))
    except ValueError:
        return None


def _iges_file(payload: bytes):
    """Parse IGES bytes with the repository reader (imported lazily: keeps this cheap).

    ``tensorcad.iges.read_iges`` is path-based, so the bytes are staged through a
    temporary ``.igs`` file; nothing is written next to the caller's data.
    """
    import tempfile

    from tensorcad.iges import read_iges  # local import: heavy module

    with tempfile.NamedTemporaryFile("wb", suffix=".igs", delete=False) as handle:
        handle.write(payload)
        staged = Path(handle.name)
    try:
        return read_iges(str(staged))
    finally:
        staged.unlink(missing_ok=True)


def iges_fingerprint(payload: bytes) -> dict[str, Any]:
    """Global-section fingerprint of an IGES file: unit, scale, header consistency.

    ``header_max_coordinate`` is the file's own claim about its largest coordinate; it is
    returned so the caller can check it against the measured extent instead of trusting
    either one alone.
    """
    parsed = _iges_file(payload)
    fields: Sequence[str] = tuple(parsed.global_fields)
    unit_flag = None
    raw_flag = fields[_IGES_UNIT_FLAG] if len(fields) > _IGES_UNIT_FLAG else None
    try:
        unit_flag = int(str(raw_flag).strip())
    except (TypeError, ValueError):
        unit_flag = None
    table = _unit_table()
    metres_per_unit = table.get(unit_flag) if unit_flag is not None else None
    if metres_per_unit is None:
        # Fall back to the reader's own resolution, but say so.
        metres_per_unit = parsed.metres_per_unit
    return {
        "units": parsed.units,
        "unit_flag": unit_flag,
        "unit_name": (fields[_IGES_UNIT_NAME].strip() if len(fields) > _IGES_UNIT_NAME else ""),
        "metres_per_unit": metres_per_unit,
        "model_space_scale": parsed.model_scale,
        "header_max_coordinate": _float(fields[_IGES_MAX_COORD] if len(fields) > _IGES_MAX_COORD else None),
        "min_resolution": _float(
            fields[_IGES_MIN_RESOLUTION] if len(fields) > _IGES_MIN_RESOLUTION else None
        ),
        "generation_stamp": (fields[_IGES_GEN_STAMP].strip() if len(fields) > _IGES_GEN_STAMP else ""),
        "author": (fields[_IGES_AUTHOR].strip() if len(fields) > _IGES_AUTHOR else ""),
        "organisation": (
            fields[_IGES_ORGANISATION].strip() if len(fields) > _IGES_ORGANISATION else ""
        ),
        "entity_census": parsed.stats(),
        "metres_resolved_from": "unit_flag" if unit_flag in table else "reader",
    }


def iges_extent_m(payload: bytes, samples: int = 15) -> dict[str, Any]:
    """Bounding extent of every B-spline surface in an IGES file, in metres.

    The extent is measured by evaluating each surface on a ``samples x samples`` grid
    inside its own parameter domain.  Control-net extents are *not* used: a pole-heavy
    patch at the bow overshoots by several percent, which is exactly the size of error
    this function exists to detect.
    """
    parsed = _iges_file(payload)
    fingerprint = iges_fingerprint(payload)
    scale = float(fingerprint["metres_per_unit"]) * float(parsed.model_scale or 1.0)
    low = [math.inf] * 3
    high = [-math.inf] * 3
    counted = 0
    for surface in parsed.surfaces():
        counted += 1
        for point in surface.sample(samples, samples):
            for axis in range(3):
                value = point[axis] * scale
                if value < low[axis]:
                    low[axis] = value
                if value > high[axis]:
                    high[axis] = value
    if counted == 0:
        raise ValueError("IGES file contains no B-spline surface to measure")
    extent = [high[i] - low[i] for i in range(3)]
    header_max = fingerprint["header_max_coordinate"]
    # The header maximum-coordinate field is in file units, so compare in metres.
    header_max_m = header_max * scale if header_max else None
    measured_extreme = max(abs(low[i]) for i in range(3)) if any(
        math.isfinite(v) for v in low
    ) else 0.0
    measured_extreme = max(measured_extreme, max(abs(high[i]) for i in range(3)))
    return {
        "extent_m": [round(v, 6) for v in extent],
        "low_m": [round(v, 6) for v in low],
        "high_m": [round(v, 6) for v in high],
        "surfaces_measured": counted,
        "grid": [samples, samples],
        "method": "surface evaluation on a parameter grid (not the control net)",
        "header_max_coordinate_m": (round(header_max_m, 6) if header_max_m else None),
        "measured_extreme_m": round(measured_extreme, 6),
        "header_matches_extent": (
            None
            if not header_max_m
            else bool(abs(header_max_m - measured_extreme) <= 0.01 * max(header_max_m, 1e-9))
        ),
    }


# --------------------------------------------------------------------------- scale
def axis_denominators(
    extent_m: Sequence[float], published_full_scale: dict[str, float]
) -> dict[str, dict[str, float]]:
    """Per-axis scale denominator ``published_full / measured``.

    ``published_full_scale`` maps an axis to ``{"value": .., "quantity": ..}``.  A hull
    file does not necessarily carry Lpp along x (overhangs, spacer patches), so each axis
    is reported on its own and the caller decides which axes are trustworthy.
    """
    names = "xyz"
    out: dict[str, dict[str, float]] = {}
    for axis, index in zip(names, range(3)):
        entry = published_full_scale.get(axis)
        if not entry:
            continue
        measured = float(extent_m[index])
        if measured <= 0.0:
            raise ValueError(f"non-positive extent on axis {axis}")
        out[axis] = {
            "quantity": entry["quantity"],
            "published_full": float(entry["value"]),
            "measured_m": measured,
            "denominator": float(entry["value"]) / measured,
        }
    return out


def classify_regime(extent_m: Sequence[float], published_lwl: float) -> dict[str, Any]:
    """Decide whether an archive is at a published denominator or normalised to unit length.

    Unit-normalised publisher exports exist (the T2015 ONRT archive is one): their largest
    coordinate is exactly 1.0 and the hull proportions are right, but nothing in the file
    says "3.147 m".  Telling those two regimes apart matters, because a mesher that
    assumes model scale will build a hull 3x too small.
    """
    length = float(max(extent_m))
    if length <= 0.0:
        raise ValueError("degenerate extent")
    normalised = abs(length - 1.0) <= 0.02
    return {
        "extent_max_m": round(length, 6),
        "published_lwl_m": published_lwl,
        "regime": "unit_length_normalised" if normalised else "finite_length",
        "implied_denominator": (published_lwl / length) if length else None,
        "note": (
            "largest extent is 1.0 +-2 %, i.e. the export is normalised by a hull length; "
            "the physical scale must come from the published particulars, not from the file"
            if normalised
            else "the export carries a physical length; a denominator can be read from it"
        ),
    }


def relative_error(measured: float, reference: float) -> float:
    """Absolute relative error, with the reference guarded against zero."""
    if reference == 0.0:
        raise ValueError("zero reference")
    return abs(measured - reference) / abs(reference)


def check_row(
    quantity: str,
    measured: float,
    reference: float,
    *,
    unit: str,
    source: str,
    note: str = "",
) -> dict[str, Any]:
    """One row of the validation table, with the <3 % verdict spelled out."""
    error = relative_error(measured, reference)
    return {
        "quantity": quantity,
        "measured": measured,
        "reference": reference,
        "unit": unit,
        "reference_source": source,
        "abs_rel_error": error,
        "within_tolerance": bool(error < TOLERANCE_REL),
        "tolerance": TOLERANCE_REL,
        "note": note,
    }


# --------------------------------------------------------------------------- pages
def strip_html(markup: str) -> str:
    """HTML -> single-line text, so a table row becomes one whitespace-normalised string."""
    import html as _html

    text = re.sub(r"<script.*?</script>", " ", markup, flags=re.S | re.I)
    text = re.sub(r"<style.*?</style>", " ", text, flags=re.S | re.I)
    text = re.sub(r"<[^>]+>", " ", text)
    return re.sub(r"\s+", " ", _html.unescape(text))


def scrape_scalar(markup: str, pattern: str, *, cast: type = float) -> dict[str, Any]:
    """Read one scalar out of a publisher page, keeping the matched row text as evidence.

    The snippet is what makes the number auditable: a reader can search the page for it
    instead of trusting a transcribed digit.
    """
    text = strip_html(markup)
    match = re.search(pattern, text)
    if not match:
        raise ValueError(f"pattern not found in page text: {pattern}")
    raw = match.group(1)
    return {
        "value": cast(raw),
        "raw": raw,
        "snippet": text[max(0, match.start() - 40) : match.end() + 40].strip(),
    }


def scrape_pair(markup: str, pattern: str) -> dict[str, Any]:
    """Read a ``model | full scale`` pair from a publisher page (two numbers in one row)."""
    text = strip_html(markup)
    match = re.search(pattern, text)
    if not match:
        raise ValueError(f"pattern not found in page text: {pattern}")
    return {
        "values": [float(match.group(1)), float(match.group(2))],
        "raw": [match.group(1), match.group(2)],
        "snippet": text[max(0, match.start() - 40) : match.end() + 40].strip(),
    }


def find_published_quote_entries(
    crossref: dict[str, Any], needle: str
) -> list[dict[str, Any]]:
    """Find crossref entries whose preserved cell text contains ``needle``.

    Used for the JBC scale factor, which the publisher prints as ``Scale factor 1 40``
    inside a table cell; the crossref keeps the cell text verbatim, so the quote -- not a
    transcription -- is what the ledger cites.
    """
    hits: list[dict[str, Any]] = []
    for ship, block in (crossref.get("report") or {}).items():
        for quantity, body in (block.get("quantities") or {}).items():
            for entry in body.get("entries") or []:
                cell = entry.get("cell") or {}
                label = " ".join(str(v) for v in cell.values() if isinstance(v, str))
                if needle in label:
                    hits.append(
                        {
                            "ship": ship,
                            "quantity": quantity,
                            "cell_text": label,
                            "source": entry.get("source") or cell.get("source"),
                            "url": entry.get("url") or cell.get("url"),
                            "sha256": entry.get("sha256") or cell.get("sha256"),
                        }
                    )
    return hits


def summarise_check_rows(rows: Iterable[dict[str, Any]]) -> dict[str, Any]:
    """Count the validation rows and return the worst relative error, for the headline."""
    rows = list(rows)
    if not rows:
        return {"rows": 0, "within_tolerance": 0, "worst_abs_rel_error": None, "all_within": None}
    worst = max(rows, key=lambda r: r["abs_rel_error"])
    return {
        "rows": len(rows),
        "within_tolerance": sum(1 for r in rows if r["within_tolerance"]),
        "worst_abs_rel_error": worst["abs_rel_error"],
        "worst_quantity": worst["quantity"],
        "all_within": all(r["within_tolerance"] for r in rows),
    }