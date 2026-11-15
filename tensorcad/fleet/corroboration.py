"""Arithmetic corroboration of a published cell — a second check that is *not* a vote.

Anchor ``REPRO-GEOM-FLEET-REFS``.  The cross-source library
(:mod:`tensorcad.fleet.crossref`) counts a quantity ``CROSS_VERIFIED`` only when two
independent **publishers** state it.  Some documents state the same physical quantity in a
different *form* instead: a model-scale length where the library keeps a full-scale
percentage of Lpp (LCB), a mass where the library keeps a volume.  Two separate defects
live in that gap:

* refusing to use such a document throws away an independent measurement — the model
  basin printed the number it measured, it just printed it in metres at 1/49 scale;
* promoting it to a vote would be worse: the full-scale number would then be *this tool's
  multiplication*, so a transcription slip in the scale or in Lpp would be laundered into
  "two sources agree".

This module therefore produces a third kind of statement, kept apart from the census:

    ``CORROBORATES`` -- an independently *measured* numeral, converted with the document's
    own stated scale and the library's own Lpp, lands within tolerance of the published
    percentage.  It is recorded as corroboration, never merged into the family census, and
    never allowed to change a verdict (``votes`` is always ``False``).

The conversion is invertible and stated explicitly so a reader can redo it on paper::

    LCB [% Lpp, fwd+]  ->  distance aft of FP = Lpp/2 + pct/100 * Lpp
    model length * scale -> full-scale length

Only quantities where that algebra is exact (a length, a percentage of a stated Lpp) are
handled here; volume/mass conversions stay where they are, in ``crossref.resolve_basis``,
because there a *closure residual* must decide the reading and arithmetic alone cannot.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Any

__all__ = [
    "LCB_TOLERANCE",
    "LcbCorroboration",
    "lcb_aft_fp_from_percent",
    "lcb_percent_from_aft_fp",
    "corroborate_lcb_by_scale",
]

#: Relative tolerance for calling an arithmetic corroboration agreement.  Deliberately the
#: same 0.5 % the library uses between publishers: a conversion that cannot reach the same
#: band is not corroboration, it is a contradiction worth looking at.
LCB_TOLERANCE = 0.005


def lcb_aft_fp_from_percent(percent_fwd_positive: float, lpp: float) -> float:
    """Distance of the LCB aft of the forward perpendicular, from a percentage of Lpp.

    A positive percentage means the LCB is forward of midship, i.e. *less* than Lpp/2 behind
    the FP; a negative percentage puts it further aft.
    """
    if lpp <= 0:
        raise ValueError(f"Lpp must be positive, got {lpp!r}")
    return lpp / 2.0 - (percent_fwd_positive / 100.0) * lpp


def lcb_percent_from_aft_fp(distance_aft_fp: float, lpp: float) -> float:
    """Inverse of :func:`lcb_aft_fp_from_percent` (percentage of Lpp, forward positive)."""
    if lpp <= 0:
        raise ValueError(f"Lpp must be positive, got {lpp!r}")
    return ((lpp / 2.0) - distance_aft_fp) / lpp * 100.0


@dataclass(frozen=True)
class LcbCorroboration:
    """One model-scale LCB re-expressed in the library's form, plus its residual."""

    ship: str
    model_length_aft_fp_m: float
    scale: float
    full_scale_aft_fp_m: float
    lpp_m: float
    published_percent_fwd_positive: float
    published_aft_fp_m: float
    residual_m: float
    residual_relative: float
    tolerance: float
    verdict: str
    votes: bool = False

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)


def corroborate_lcb_by_scale(
    *,
    ship: str,
    model_length_aft_fp_m: float,
    scale: float,
    lpp_m: float,
    published_percent_fwd_positive: float,
    tolerance: float = LCB_TOLERANCE,
) -> LcbCorroboration:
    """Convert a model-scale ``LCB (m aft of FP)`` and compare it with a published percentage.

    ``scale`` is the model scale *denominator* the document states (1/49 -> ``49.0``).  The
    result never votes: it is a cross-check of one published number, not a second publisher.
    """
    for name, value in (("model_length_aft_fp_m", model_length_aft_fp_m), ("scale", scale)):
        if not float(value) > 0:
            raise ValueError(f"{name} must be positive, got {value!r}")
    full_scale = model_length_aft_fp_m * scale
    published_m = lcb_aft_fp_from_percent(published_percent_fwd_positive, lpp_m)
    residual = full_scale - published_m
    relative = abs(residual) / published_m
    return LcbCorroboration(
        ship=ship,
        model_length_aft_fp_m=model_length_aft_fp_m,
        scale=scale,
        full_scale_aft_fp_m=round(full_scale, 6),
        lpp_m=lpp_m,
        published_percent_fwd_positive=published_percent_fwd_positive,
        published_aft_fp_m=round(published_m, 6),
        residual_m=round(residual, 6),
        residual_relative=relative,
        tolerance=tolerance,
        verdict="CORROBORATES" if relative <= tolerance else "INCONSISTENT",
        votes=False,
    )