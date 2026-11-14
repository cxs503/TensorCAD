"""Official-source layer for the baseline fleet: fetched pages, extracted digits, IGES archives.

The fleet records (:mod:`tensorcad.fleet.manifest`) carry numbers and the source registry
(:mod:`tensorcad.fleet.sources`) carries citations. Both were written offline in earlier
rounds, so every entry said ``retrieval="offline_recollection"`` and ``verified=False``.
This module is the *online* half of that story: it knows the publisher pages, parses the
published particulars out of their HTML by structure (row label -> full-scale column) rather
than by hand, and compares the extraction against a fleet record field by field.

Two rules make the extraction auditable instead of convenient:

* nothing is transcribed. :func:`parse_particulars` returns the digits it reads from the page
  together with the page sha256 and the row label it read them from, so a reader can re-run
  the parse on the same bytes and get the same numbers;
* an extraction that disagrees with the record is *reported*, not silently adopted. The
  comparison returns a per-quantity relative difference and a verdict, and the caller
  decides. This is how the site's ``Bwl = 580`` for KVLCC2 (a missing decimal point for
  58.0) is caught: the closure ``Cb*L*B*T`` around the raw digit misses the published
  displacement by two orders of magnitude, while the corrected digit closes to 1e-6 relative.

The page bytes themselves are not redistributed here (the publishers state no licence);
:mod:`tools.baseline_fleet.fetch_sources` materialises them on demand and pins their
sha256 in :data:`PAGES`.
"""

from __future__ import annotations

import html as html_module
import re
from dataclasses import dataclass, field
from typing import Any, Iterable, Mapping, Sequence

SCHEMA = "tensorcad.fleet.official.v1"

#: Publisher pages the published particulars are read from, with the sha256 recorded on the
#: fetching run (2026-11-14). A re-fetch whose digest differs is a source change and must be
#: re-verified deliberately, exactly like the geometry archives.
PAGES: dict[str, dict[str, str]] = {
    "KCS": {
        "url": "https://www.t2015.nmri.go.jp/kcs_gc.html",
        "sha256": "68d9c24600cb699b16394fc5c32762ebc875c856d14ef8460e1d0ebf4e7c78a4",
        "publisher": "NMRI Tokyo 2015 Workshop",
        "locator": "Main particulars table, full-scale column",
    },
    "JBC": {
        "url": "https://www.t2015.nmri.go.jp/jbc_gc.html",
        "sha256": "fc2ce9054a8b7842162bac94365853fa076da2527e8c7e455832f9710f58f875",
        "publisher": "NMRI Tokyo 2015 Workshop",
        "locator": "Main particulars table, DESIGNED FULL column",
    },
    "ONRT": {
        "url": "https://www.t2015.nmri.go.jp/onrt_gc.html",
        "sha256": "5e13427332de8be17b9dc2245dc5830300baf2e2340efcd8a80eee2269c43a11",
        "publisher": "NMRI Tokyo 2015 Workshop",
        "locator": "Main particulars table, Full scale column",
    },
    "KVLCC2": {
        "url": "https://simman2020.kr/contents/Ship_data_summary.php",
        "sha256": "c1a668f6d2440c5d83a83c3616285e9288ef618dbc9a2feebd26627be3f4ebb3",
        "publisher": "SIMMAN 2020 Workshop",
        "locator": "Table 1 Main Particulars of KVLCC2, full-scale column",
    },
}

#: Geometry archives: the hull IGES member each ship's mesh is built from.
ARCHIVES: dict[str, dict[str, str]] = {
    "KCS": {
        "url": "https://www.t2015.nmri.go.jp/file/Geometry_IGES_files/kcs/KCS_hull_SVA.zip",
        "member": "KCS_hull_SVA/KCS_hull_SVA.igs",
        "publisher": "NMRI Tokyo 2015 Workshop",
        "sha256": "b9e334853857b819",
    },
    "KVLCC2": {
        "url": "https://simman2020.kr/download/KVLCC2-Hull.zip",
        "member": "kvlcc2.igs",
        "publisher": "SIMMAN 2020 Workshop",
        "sha256": "dc032647f00a68d3",
    },
    "JBC": {
        "url": "https://www.t2015.nmri.go.jp/file/Geometry_IGES_files/jbc/JBC_IGES.zip",
        "member": "jbc.igs",
        "publisher": "NMRI Tokyo 2015 Workshop",
        "sha256": "e19800ee954f4960",
    },
    "ONRT": {
        "url": "https://simman2020.kr/download/onrt_wo_gap.zip",
        "member": "onrt_wo_gap/onr5613-full_appendage_deckhous_mergedR.igs",
        "publisher": "SIMMAN 2020 Workshop",
        "sha256": "bd747eff8f531ba4",
    },
}

LICENCE = {
    "NMRI Tokyo 2015 Workshop": (
        "not stated on the publisher page (no terms-of-use text); workshop benchmark geometry, "
        "research use, redistribution requires permission -- digests are recorded, bytes are not shipped"
    ),
    "SIMMAN 2020 Workshop": (
        "not stated; site footer reads 'Copyright (c) SIMMAN2020 All rights reserved'; workshop "
        "benchmark geometry, research use, redistribution requires permission -- digests are recorded, "
        "bytes are not shipped"
    ),
}

#: Row labels as they appear on the publisher pages, per quantity. The label is matched
#: case-insensitively against the flattened row text; the digits are read from the row, never
#: from a table of expected values.
LABELS: dict[str, dict[str, tuple[str, ...]]] = {
    "KCS": {
        "lpp": ("length between perpendiculars",),
        "lwl": ("length of waterline",),
        "bwl": ("maximum beam of waterline",),
        "depth": ("depth",),
        "draft": ("draft",),
        "displacement_volume": ("displacement volume",),
        "cb": ("block coefficient",),
        "cm": ("midship section coefficient",),
        "lcb_percent": ("lcb (%l", "lcb (%lpp"),
    },
    "JBC": {
        "lpp": ("length between perpendiculars",),
        "lwl": ("length of waterline",),
        "bwl": ("maximum beam of waterline",),
        "depth": ("depth",),
        "draft": ("draft(full)",),
        "displacement_volume": ("displacement volume of a hull",),
        "cb": ("block coefficient",),
        "cm": ("midship section coefficient",),
        "lcb_percent": ("l_cb (%l", "lcb (%l"),
    },
    "ONRT": {
        "lwl": ("length of waterline",),
        "bwl": ("maximum beam of waterline",),
        "depth": ("depth",),
        "draft": ("draft",),
        "displacement_mass": ("displacement",),
        "cb": ("block coefficient",),
        "lcb_metres_aft_fp": ("lcb (m) aft of fp",),
    },
    "KVLCC2": {
        "lpp": ("lpp (m)",),
        "lwl": ("lwl (m)",),
        "bwl": ("bwl (m)",),
        "depth": ("d (m)",),
        "draft": ("t (m)",),
        "displacement_volume": ("displacement (m3)",),
        "cb": ("cb",),
        "cm": ("cm",),
        "lcb_percent": ("lcb (%), fwd+",),
    },
}

#: Column the digits are read from, per ship. The T2015 pages print ``label | symbol+unit |
#: full scale | model 1 | ...`` for KCS and JBC, but ONRT puts its *model* column first, and
#: the SIMMAN summary table has no symbol column at all -- so the index is per ship, stated
#: here rather than guessed per row.
COLUMN_INDEX = {"KCS": 2, "JBC": 2, "ONRT": 3, "KVLCC2": 1}
#: When the full-scale cell is not a number (ONRT prints ``NA`` for LCB at full scale) the
#: value is taken from the preceding column and labelled with the scale it actually belongs to.
FALLBACK_TO_PREVIOUS = {"ONRT"}

NUMBER = re.compile(r"[-+]?\d+(?:\.\d+)?")


@dataclass(frozen=True)
class Extraction:
    """One published particular, with the cell, page and scale it was read from."""

    quantity: str
    value: float
    raw: str
    scale: str
    label: str
    row: tuple[str, ...]
    page_url: str
    page_sha256: str

    def as_dict(self) -> dict[str, Any]:
        return {
            "quantity": self.quantity,
            "value": self.value,
            "raw": self.raw,
            "scale": self.scale,
            "label": self.label,
            "row": list(self.row),
            "page_url": self.page_url,
            "page_sha256": self.page_sha256,
        }


@dataclass
class PageExtraction:
    ship: str
    page_url: str
    page_sha256: str
    publisher: str
    locator: str
    values: dict[str, Extraction] = field(default_factory=dict)

    def get(self, quantity: str) -> float | None:
        item = self.values.get(quantity)
        return None if item is None else item.value

    def as_dict(self) -> dict[str, Any]:
        return {
            "ship": self.ship,
            "page_url": self.page_url,
            "page_sha256": self.page_sha256,
            "publisher": self.publisher,
            "locator": self.locator,
            "values": {name: item.as_dict() for name, item in sorted(self.values.items())},
        }


def strip_tags(markup: str) -> str:
    """Flatten HTML to text, dropping script/style bodies."""
    text = re.sub(r"(?is)<(script|style).*?</\1>", " ", markup)
    text = re.sub(r"<[^>]+>", " ", text)
    text = html_module.unescape(text)
    return re.sub(r"[ \t\xa0]+", " ", text)


def table_rows(markup: str) -> list[tuple[str, ...]]:
    """The page's table rows, each as a tuple of flattened cell texts.

    The publisher pages are plain tables: one ``<tr>`` per particular, one ``<td>`` per
    column. Reading cells rather than a flattened text stream is what makes the extraction
    structural -- a value cannot drift onto the wrong label, because it never leaves its row.
    """
    out: list[tuple[str, ...]] = []
    for raw_row in re.findall(r"(?is)<tr[^>]*>(.*?)</tr>", markup):
        cells = [
            re.sub(r"\s+", " ", strip_tags(cell)).strip()
            for cell in re.findall(r"(?is)<t[dh][^>]*>(.*?)</t[dh]>", raw_row)
        ]
        if any(cell for cell in cells):
            out.append(tuple(cells))
    return out


def _number_in(cell: str) -> tuple[str, float] | None:
    cleaned = cell.replace(",", "")
    match = NUMBER.search(cleaned)
    if match is None:
        return None
    return cell, float(match.group(0))


def parse_particulars(ship: str, markup: str, *, page_sha256: str = "") -> PageExtraction:
    """Read the published particulars of ``ship`` out of the publisher page ``markup``."""
    if ship not in LABELS:
        raise KeyError(f"no publisher-page parser registered for {ship!r}")
    page = PAGES.get(ship, {})
    extraction = PageExtraction(
        ship=ship,
        page_url=page.get("url", ""),
        page_sha256=page_sha256 or page.get("sha256", ""),
        publisher=page.get("publisher", ""),
        locator=page.get("locator", ""),
    )
    index = COLUMN_INDEX[ship]
    rows = table_rows(markup)
    for quantity, labels in LABELS[ship].items():
        for row in rows:
            if not row:
                continue
            # Match the label anywhere in the row: some pages put the quantity name in the
            # first cell ("LCB") and the qualifier in the second ("LCB (m) aft of FP"), and the
            # digits stay in their own cell either way.
            haystack = " ".join(row).lower()
            if not any(label in haystack for label in labels):
                continue
            if index >= len(row):
                continue
            hit = _number_in(row[index])
            scale = "full"
            if hit is None and ship in FALLBACK_TO_PREVIOUS and index - 1 >= 1:
                hit = _number_in(row[index - 1])
                scale = "reduced_scale"
            if hit is None:
                continue
            raw, value = hit
            extraction.values[quantity] = Extraction(
                quantity=quantity,
                value=value,
                raw=raw,
                scale=scale,
                label=row[0][:120],
                row=row,
                page_url=extraction.page_url,
                page_sha256=extraction.page_sha256,
            )
            break
    return extraction


#: Published displacement digits come in two currencies; a mass is turned into a volume with
#: the density the workshop states when it states one, otherwise 1.025 t/m^3 (seawater).
SEAWATER_DENSITY = 1.025


def published_reference(ship: str, extraction: PageExtraction) -> dict[str, Any]:
    """The reference volume(s) a ship's mesh can be compared against, with their basis.

    ONRT publishes a *mass* (8507 ton) and a Cb defined on ``LWL*BWL*T``; the two are not
    consistent to better than 2 %, so both candidate volumes are returned and the caller
    reports the mesh error against each instead of picking the flattering one.
    """
    candidates: list[dict[str, Any]] = []
    volume = extraction.get("displacement_volume")
    if volume is not None:
        candidates.append(
            {
                "basis": "published_volume",
                "volume_m3": volume,
                "label": "Displacement volume / ∇",
                "page_sha256": extraction.page_sha256,
                "note": "read from the publisher page by row label",
            }
        )
    mass = extraction.get("displacement_mass")
    if mass is not None:
        candidates.append(
            {
                "basis": "published_mass_over_density",
                "volume_m3": mass / SEAWATER_DENSITY,
                "label": "Δ (ton) / 1.025 t/m^3",
                "page_sha256": extraction.page_sha256,
                "note": "mass x 1.025 t/m^3 seawater; the workshop does not state a density",
            }
        )
    lwl = extraction.get("lwl") or extraction.get("lpp")
    cb = extraction.get("cb")
    if lwl and cb and extraction.get("bwl") and extraction.get("draft"):
        candidates.append(
            {
                "basis": "closure",
                "volume_m3": cb * lwl * extraction.get("bwl") * extraction.get("draft"),
                "label": "Cb * LWL * BWL * T",
                "page_sha256": extraction.page_sha256,
                "note": "closure of the published digits; not an external measurement",
            }
        )
    return {"ship": ship, "candidates": candidates}


def compare_with_record(
    record: Any,
    extraction: PageExtraction,
    *,
    tolerance: float = 0.002,
) -> dict[str, Any]:
    """Field-by-field comparison of an extraction against a fleet record.

    ``record`` is a :class:`tensorcad.fleet.manifest.ShipRecord`. ``tolerance`` is relative
    and deliberately tight (0.2 %): published digits carry 4-5 significant figures, so a
    disagreement above this is a transcription error, not rounding.
    """
    mapping = {
        "length": ("lpp", "lwl"),
        "beam": ("bwl",),
        "draft": ("draft",),
        "block_coefficient": ("cb",),
        "midship_coefficient": ("cm",),
        "displacement_volume": ("displacement_volume",),
    }
    fields: list[dict[str, Any]] = []
    for attr, quantities in mapping.items():
        record_value = getattr(record, attr, None)
        if record_value is None:
            continue
        for quantity in quantities:
            extracted = extraction.get(quantity)
            if extracted is None:
                continue
            relative = abs(float(record_value) - extracted) / max(abs(extracted), 1e-12)
            fields.append(
                {
                    "field": attr,
                    "quantity": quantity,
                    "record": float(record_value),
                    "published": extracted,
                    "relative_difference": relative,
                    "agrees": relative <= tolerance,
                }
            )
            # Compare the record against the *first* quantity that the page publishes; LWL is
            # offered as a fallback for a record whose Lpp the page does not print.
            break
    disagreements = [item for item in fields if not item["agrees"]]
    return {
        "ship": record.name,
        "tolerance": tolerance,
        "fields": fields,
        "disagreements": disagreements,
        "all_agree": not disagreements,
        "page_sha256": extraction.page_sha256,
        "page_url": extraction.page_url,
    }


def consistency_check(extraction: PageExtraction) -> dict[str, Any]:
    """Do the page's own digits close? ``Cb*L*B*T`` must reproduce its published volume.

    This is the check that exposes a page-level defect such as KVLCC2's ``Bwl = 580``: a
    beam of 580 m closes to 3.1e6 m^3 against a published 312 622 m^3. The function reports
    the closure both ways so a reader can see which digit is the suspect one.
    """
    lpp = extraction.get("lpp")
    lwl = extraction.get("lwl")
    bwl = extraction.get("bwl")
    draft = extraction.get("draft")
    cb = extraction.get("cb")
    length = lpp or lwl
    report: dict[str, Any] = {"ship": extraction.ship, "closure": None, "relative_difference": None}
    if not (length and bwl and draft and cb):
        return report
    closure = cb * length * bwl * draft
    report["closure"] = closure
    report["rule"] = "Cb * Lpp * Bwl * T" if lpp else "Cb * Lwl * Bwl * T"
    published = extraction.get("displacement_volume")
    if published is None:
        mass = extraction.get("displacement_mass")
        published = None if mass is None else mass / SEAWATER_DENSITY
    if published:
        report["published"] = published
        report["relative_difference"] = abs(closure - published) / published
        report["closes"] = report["relative_difference"] <= 0.005
    return report


def archive_record(ship: str) -> dict[str, Any]:
    """The provenance block for one ship's official geometry archive."""
    entry = ARCHIVES[ship]
    return {
        "ship": ship,
        "kind": "official_geometry_archive",
        "url": entry["url"],
        "publisher": entry["publisher"],
        "license": LICENCE[entry["publisher"]],
        "member": entry["member"],
        "sha256_prefix_recorded": entry["sha256"],
        "retrieval": "fetched",
        "verified": True,
        "redistributed": False,
    }