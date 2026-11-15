"""Longitudinal reference frames for LCB / LCG readings.

Why this module exists
---------------------
A longitudinal particular is not a number, it is *a number in a frame*:

* the reference length it is a percentage of -- ``Lpp`` or ``Lwl``;
* the point it is measured from -- AP, FP, the ``Lpp`` mid-point or the ``Lwl``
  mid-point (which sits forward of the ``Lpp`` mid-point by half the bulb overhang
  whenever a bulb protrudes ahead of the forward perpendicular);
* the sign convention -- forward-positive or aft-positive.

Two publishers can print ``1.68`` and ``1.625`` for the same physical location and
neither is wrong, and one publisher can print ``+1.48`` and ``-1.48`` for the same
location as another because the sign convention differs.  Comparing the printed
numerals directly therefore produces two kinds of false result: a *false conflict*
(a real disagreement invented by the frame) and a *false agreement* (a real
disagreement cancelled by the frame).  On KCS the frame term is not academic: the
bulb protrudes 2.5 m ahead of the FP over a 230 m length, so the ``Lwl`` mid-point
sits 0.543 % ``Lpp`` forward of the ``Lpp`` mid-point.

This module converts a published longitudinal reading into the library's canonical
quantity ``lcb_percent`` -- **percent of Lpp, positive forward, measured from the
Lpp mid-point** -- and, more importantly, classifies *how much of that conversion
the document itself states*.  The classification is the point: a frame conversion
that only changes the address of a full-scale numeral preserves the publisher's
information, whereas multiplying a model-scale numeral by a scale factor, or
substituting LCG for LCB, adds arithmetic of our own.  Only the first may vote on a
published value; the others are recorded as corroboration (see
``tensorcad.fleet.crossref.votes`` and :mod:`tensorcad.fleet.corroboration`).

The rule set is deliberately small and each rule is exercised in
``tests/test_fleet_lcb_frames.py``.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

#: The canonical frame every reading is converted into.
CANONICAL = "percent_Lpp_fwd_from_lpp_mid"

#: Provenance classes, from strongest to weakest.  Only ``EXPLICIT_FRAME`` votes.
EXPLICIT_FRAME = "explicit_frame"          # datum + reference length stated, numeral full scale
SCALE_DERIVED = "scale_derived"            # model-scale numeral, document states the scale
IDENTITY_DERIVED = "identity_derived"      # LCG used for LCB, document states even keel
COMPUTED_NOT_PUBLISHED = "computed_not_published"  # the source's own computation, not a particular
UNSTATED = "unstated"                      # datum or reference length not stated -> unusable

PROVENANCE_CLASSES = (
    EXPLICIT_FRAME,
    SCALE_DERIVED,
    IDENTITY_DERIVED,
    COMPUTED_NOT_PUBLISHED,
    UNSTATED,
)

#: Which provenance classes are allowed to vote on a published value.  A conversion
#: that needs arithmetic of ours is corroboration, never a second publisher.
VOTING_CLASSES = frozenset({EXPLICIT_FRAME})

#: What a class has to satisfy beyond its own definition before it can be trusted even
#: as corroboration.  ``IDENTITY_DERIVED`` requires the document to state even keel; a
#: trimmed condition makes LCG and LCB differ by the trimming moment, so the identity
#: would be wrong rather than merely imprecise.
_DATUM_ALIASES = {
    "ap": "ap",
    "after perpendicular": "ap",
    "aft perpendicular": "ap",
    "fp": "fp",
    "forward perpendicular": "fp",
    "lpp_mid": "lpp_mid",
    "lpp mid": "lpp_mid",
    "midship": "lpp_mid",
    "lwl_mid": "lwl_mid",
    "lwl mid": "lwl_mid",
    "waterline mid-point": "lwl_mid",
    "waterline midpoint": "lwl_mid",
}


class FrameError(ValueError):
    """A frame that cannot be interpreted, or a datum the converter refuses to guess."""


@dataclass(frozen=True)
class Frame:
    """Where a longitudinal numeral is measured from, and in what scale.

    ``datum``    one of ``ap``, ``fp``, ``lpp_mid``, ``lwl_mid``.
    ``sign``     ``fwd`` (positive forward) or ``aft`` (positive aft).
    ``reference`` ``Lpp`` or ``Lwl`` -- the length a published percentage is taken of.
    ``scale``    model-to-full-scale factor of the column the numeral was read from
                 (1.0 for a full-scale column).  It is kept here, not applied silently,
                 because applying it is exactly what demotes a cell to corroboration.
    """

    datum: str
    sign: str
    reference: str = "Lpp"
    scale: float = 1.0

    def __post_init__(self) -> None:
        datum = _DATUM_ALIASES.get(str(self.datum).strip().lower(), str(self.datum).strip().lower())
        if datum not in ("ap", "fp", "lpp_mid", "lwl_mid"):
            raise FrameError(f"uninterpretable datum {self.datum!r}")
        object.__setattr__(self, "datum", datum)
        if self.sign not in ("fwd", "aft"):
            raise FrameError(f"uninterpretable sign convention {self.sign!r}")
        if self.reference not in ("Lpp", "Lwl"):
            raise FrameError(f"uninterpretable reference length {self.reference!r}")
        if not self.scale > 0:
            raise FrameError(f"scale must be positive, got {self.scale!r}")


def _lwl_mid_offset(lpp: float, lwl: float | None, stated: float | None) -> float:
    """Distance of the ``Lwl`` mid-point forward of the ``Lpp`` mid-point [m].

    A hull whose waterline runs *aft* of the AP only (no forward overhang) has the two
    mid-points coincident; a bulb that protrudes forward moves the ``Lwl`` mid-point
    forward.  The offset is therefore taken from the document when it states one, and
    otherwise from ``(Lwl - Lpp) / 2`` -- valid only when the whole difference is forward
    overhang, which is the case for KCS (bulb) and is *not* assumed silently: the caller
    records which of the two it used.
    """
    if stated is not None:
        return float(stated)
    if lwl is None:
        raise FrameError("a Lwl datum needs either the document's offset or a Lwl value")
    return (float(lwl) - float(lpp)) / 2.0


def to_x_from_ap(
    value: float,
    frame: Frame,
    lpp: float,
    lwl: float | None = None,
    lwl_mid_offset: float | None = None,
) -> float:
    """Abscissa of the reading from the AP, in full-scale metres.

    ``value`` is in the scale of the column it was read from; ``lpp`` (and ``lwl``) are
    full scale, so a model-scale column is converted with ``frame.scale`` first.
    """
    if not lpp > 0:
        raise FrameError(f"Lpp must be positive, got {lpp!r}")
    v = float(value) * frame.scale
    if frame.datum == "ap":
        x = v if frame.sign == "fwd" else lpp - v
    elif frame.datum == "fp":
        x = lpp - v if frame.sign == "aft" else v
    elif frame.datum == "lpp_mid":
        x = lpp / 2.0 + (v if frame.sign == "fwd" else -v)
    else:  # lwl_mid
        offset = _lwl_mid_offset(lpp, lwl, lwl_mid_offset)
        x = lpp / 2.0 + offset + (-v if frame.sign == "aft" else v)
    return x


def to_percent_lpp(
    value: float,
    frame: Frame,
    lpp: float,
    lwl: float | None = None,
    lwl_mid_offset: float | None = None,
) -> float:
    """Canonical reading: percent of Lpp, forward-positive, from the Lpp mid-point."""
    x = to_x_from_ap(value, frame, lpp, lwl=lwl, lwl_mid_offset=lwl_mid_offset)
    return (x - float(lpp) / 2.0) / float(lpp) * 100.0


def canonicalise(record: Mapping[str, Any]) -> dict[str, Any]:
    """Convert one registry record and classify what the conversion cost.

    The record is a published *reading* (see ``benchmarks/verified/fleet_refs/lcb_frames.json``);
    the result keeps every numeral it used so the reader can redo the arithmetic:

    ``value_canonical`` ``percent_Lpp`` fwd+ from the Lpp mid-point;
    ``x_from_ap_m``     the same location as a distance from the AP;
    ``provenance``      one of :data:`PROVENANCE_CLASSES`;
    ``votes``           whether the cell may vote on the published quantity;
    ``derivation``      human-readable arithmetic, stored with any cell built from it.
    """
    frame = Frame(
        datum=record["frame"]["datum"],
        sign=record["frame"]["sign"],
        reference=record.get("frame", {}).get("reference", "Lpp"),
        scale=float(record["frame"].get("scale", 1.0)),
    )
    lpp = float(record["lpp_m"])
    lwl = record.get("lwl_m")
    offset = record.get("lwl_mid_offset_m")
    provenance = record["provenance"]
    if provenance not in PROVENANCE_CLASSES:
        raise FrameError(f"unknown provenance class {provenance!r}")

    if provenance == UNSTATED:
        raise FrameError("a reading whose datum is not stated cannot be converted")
    if provenance == IDENTITY_DERIVED and not record.get("even_keel_stated"):
        raise FrameError(
            "LCG may only stand in for LCB when the document states even keel "
            "(TFP == TAP); a trimmed condition moves LCB away from LCG"
        )

    value = float(record["value"])
    x = to_x_from_ap(value, frame, lpp, lwl=lwl, lwl_mid_offset=offset)
    canonical = (x - lpp / 2.0) / lpp * 100.0
    scale_note = f" x{frame.scale:g}" if frame.scale != 1.0 else ""
    arithmetic = (
        f"{record.get('quantity_published', 'LCB')} {value:g}{scale_note} m, "
        f"{frame.datum} datum, {frame.sign}+ -> x={x:.4f} m from AP; "
        f"({x:.4f} - {lpp / 2:.3f}) / {lpp:.3f} = {canonical:+.4f} %Lpp"
    )
    if provenance == IDENTITY_DERIVED:
        arithmetic += " [LCG=LCB under the stated even keel]"
    return {
        "id": record["id"],
        "ship": record["ship"],
        "source": record["source"],
        "quantity_published": record.get("quantity_published", "LCB"),
        "value_published": value,
        "unit": record.get("unit", "m"),
        "provenance": provenance,
        "votes": provenance in VOTING_CLASSES,
        "maps_to": record.get("maps_to", "lcb_percent"),
        "value_canonical": round(canonical, 6),
        "x_from_ap_m": round(x, 4),
        "lpp_m": lpp,
        "lwl_m": lwl,
        "lwl_mid_offset_used_m": (
            _lwl_mid_offset(lpp, lwl, offset) if frame.datum == "lwl_mid" else None
        ),
        "derivation": arithmetic,
        "sha256": record.get("sha256", ""),
        "url": record.get("url", ""),
        "page": record.get("page"),
        "quote": record.get("quote", ""),
        "frame_statement": record.get("frame_statement", ""),
    }


def residual(converted: Mapping[str, Any], reference_value: float) -> dict[str, Any]:
    """How far a converted reading sits from a reference, in %Lpp and relative to it.

    A residual is reported in the canonical unit (percentage points of Lpp) *and*
    relative to the reference, because a relative figure alone lets a small reading hide
    a large absolute difference and vice versa.
    """
    delta = float(converted["value_canonical"]) - float(reference_value)
    rel = abs(delta) / abs(float(reference_value)) if reference_value else float("inf")
    return {
        "reference_value": float(reference_value),
        "delta_pct_lpp": round(delta, 6),
        "delta_relative": round(rel, 6),
    }


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_digest(record: Mapping[str, Any], search_dirs: list[Path]) -> dict[str, Any]:
    """Check a record's pinned digest against local bytes, when those bytes still exist.

    A document whose bytes were never cached -- or are no longer on this host -- reports
    ``bytes_absent``; that is *unknown*, not *verified*, and the caller must say so.
    """
    pinned = record.get("sha256") or ""
    name = record.get("file")
    if name:
        for directory in search_dirs:
            candidate = Path(directory) / name
            if candidate.exists():
                actual = sha256_file(candidate)
                return {
                    "status": "verified" if actual == pinned else "digest_mismatch",
                    "file": str(candidate),
                    "sha256_actual": actual,
                    "sha256_pinned": pinned,
                }
    return {"status": "bytes_absent", "file": None, "sha256_pinned": pinned}