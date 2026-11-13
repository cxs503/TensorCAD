"""Offset tables of three published hull forms, with provenance and closed forms.

Why this exists
---------------
:mod:`tensorcad.offset_loft` needs *tables*.  The honest source of a table is the
ship's issuing body; the honest fallback is a hull whose geometry is published in
closed form, because then the ordinates are **evaluated**, not transcribed, and the
hydrostatics have an independent reference.  This module provides the fallback in a
form that keeps the distinction visible:

* the ordinates are computed at run time from the analytic twin in
  :mod:`tensorcad.geomquality.baselines` -- there is no hand-typed number anywhere in
  the pipeline, and the ``sha256`` of the resulting table is recorded in its provenance
  sidecar so a landed CSV can be checked against the definition that made it;
* every table carries a :class:`~tensorcad.offset_loft.OffsetProvenance` with the
  citation, the method, the URL of the defining code and ``official=False``: these are
  *documented reconstructions* of published hulls, not the issuing body's tables;
* every form carries the analytic reference integrals (``volume``, ``waterplane_area``,
  ``lcb``, ``vcb``) computed by
  :func:`tensorcad.geomquality.baselines.reference_integrals` from the half-breadth
  function itself, so the loft can be judged against the hull and not only against the
  table it came from.

The three forms
---------------
``wigley_iii``
    Parabolic hull ``y = (B/2)(1-(2x/L)^2)(1-(z/T-1)^2)``; ``Cb = 4/9``,
    ``Cwp = 2/3``, ``VCB = 5T/8`` exactly.  Polynomial in ``x`` and ``z``, so both the
    table sampling and a cubic loft are exact up to quadrature: this is the case that
    pins the *method*.
``series60_cb060``
    Series 60 parent form ``y = (B/2)(1-|2xi-1|^n) zeta^(1/m)`` with ``Cp = 0.625``,
    ``n = 5/3``, ``m = 24``.  The ``zeta^(1/24)`` section is nearly rectangular, the
    hardest interpolation case of the three: it is what makes the *waterline* count
    matter, and the study reports that.
``kvlcc2_class``
    Full-scale tanker (``L = 320 m``, ``B = 58 m``, ``T = 20.8 m``) calibrated to the
    published KVLCC2 particulars ``Cb = 0.8098`` and ``LCB = -3.48 % L``.  The official
    KRISO offset table is not in this repository (see ``docs/geomquality-pipeline.md``),
    so this is a class hull: the ``Cb``/``LCB`` agreement is a consistency check of the
    representation, not an independent KVLCC2 validation.  It contributes what the other
    two cannot: a full-scale hull with a fine, transom-free run and ``L/B = 5.5``.

Usage
-----
>>> from tensorcad.offset_forms import offset_table, reference_for
>>> table = offset_table("wigley_iii", n_stations=21, n_waterlines=17)
>>> reference = reference_for("wigley_iii")
>>> reference["volume"] / (table.length * table.beam * table.design_draft)
0.444...
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Mapping

from .geomquality.baselines import HullSpecimen, kvlcc2_class, series60, wigley_iii
from .offset_loft import OffsetProvenance, OffsetTableError, ShipOffsetTable

__all__ = [
    "HullForm",
    "builtin_forms",
    "offset_table",
    "reference_for",
    "sample_offset_table",
]


@dataclass(frozen=True)
class HullForm:
    """A published hull whose half-breadths are a documented function of ``(x, z)``."""

    key: str
    specimen_factory: Callable[[], HullSpecimen]
    provenance: OffsetProvenance

    def specimen(self) -> HullSpecimen:
        return self.specimen_factory()

    def reference(self) -> dict[str, Any]:
        """Analytic integrals of the hull itself, keyed like the table integrals.

        ``vcb`` is converted from the specimen's "negative depth" convention
        (``-z`` from the waterline, which is how the closed forms are written) to height
        above the keel, and ``lcb`` from midship to the aft perpendicular, so the numbers
        line up with :class:`~tensorcad.offset_loft.TableIntegrals` without anybody
        re-deriving the sign by hand.
        """
        specimen = self.specimen()
        reference = dict(specimen.reference)
        draft = specimen.draft
        out: dict[str, Any] = {
            "volume": reference["volume"],
            "waterplane_area": reference["waterplane_area"],
            "lcb": 0.5 * specimen.length + reference["lcb_x"],
            "vcb": draft + reference["vcb_z"],
            "method": reference.get("method", "analytic half-breadth + composite Gauss"),
        }
        if "closed_form" in reference:
            out["closed_form"] = dict(reference["closed_form"])
        return out

    def table(
        self,
        *,
        n_stations: int = 21,
        n_waterlines: int = 17,
        waterline_exponent: float = 1.0,
        name: str | None = None,
    ) -> ShipOffsetTable:
        return sample_offset_table(
            self.specimen(),
            provenance=self.provenance,
            n_stations=n_stations,
            n_waterlines=n_waterlines,
            waterline_exponent=waterline_exponent,
            name=name or self.key,
        )


def sample_offset_table(
    specimen: HullSpecimen,
    *,
    provenance: OffsetProvenance,
    n_stations: int = 21,
    n_waterlines: int = 17,
    waterline_exponent: float = 1.0,
    name: str | None = None,
) -> ShipOffsetTable:
    """Evaluate ``specimen.y_of_x_zeta`` on a uniform station/waterline grid.

    Stations run from the aft perpendicular (``x = 0``) to the fore perpendicular
    (``x = L``); the model coordinates used by the analytic twin are centred, so
    ``x_model = x - L/2``.  Ordinates are clamped at ``+0`` only to fold ``-1e-17``
    round-off at the perpendiculars.

    ``waterline_exponent`` clusters the waterlines towards the keel,
    ``z_j = T (j / (n_waterlines - 1)) ** p`` -- which is what a real offset table does
    and what a power-law section needs: for ``y ~ zeta**(1/m)`` the near-keel cell carries
    a large share of the section area, and a uniform grid under-resolves it (measured:
    ``p = 1`` leaves 1.7 % displacement error on the Series 60 section, ``p = 3`` leaves
    0.03 %; see ``docs/offset_loft/results.json``).
    """
    if n_stations < 3:
        raise OffsetTableError(f"n_stations must be >= 3, got {n_stations}")
    if n_waterlines < 3:
        raise OffsetTableError(f"n_waterlines must be >= 3, got {n_waterlines}")
    length = float(specimen.length)
    draft = float(specimen.draft)
    stations = tuple(length * i / (n_stations - 1) for i in range(n_stations))
    waterlines = tuple(
        draft * (j / (n_waterlines - 1)) ** waterline_exponent for j in range(n_waterlines)
    )
    rows: list[tuple[float, ...]] = []
    for x in stations:
        x_model = x - 0.5 * length
        row = tuple(max(0.0, float(specimen.y_of_x_zeta(x_model, z / draft))) for z in waterlines)
        rows.append(row)
    return ShipOffsetTable(
        name=name or specimen.name,
        length=length,
        stations=stations,
        waterlines=waterlines,
        half_breadths=tuple(rows),
        design_draft=draft,
        provenance=provenance,
        units="m",
        waterline_exponent=waterline_exponent,
    )


_CAD_CODE_URL = "https://github.com/cxs503/TensorCAD/blob/main/tensorcad/geomquality/baselines.py"

_FORMS: dict[str, HullForm] = {
    "wigley_iii": HullForm(
        key="wigley_iii",
        specimen_factory=wigley_iii,
        provenance=OffsetProvenance(
            source=(
                "Wigley III parabolic hull form y = (B/2)(1-(2x/L)^2)(1-(z/T-1)^2); "
                "the standard resistance/seakeeping verification hull (Wigley 1934; "
                "Journee & Massie 2001, 'Offshore Hydromechanics', section 3)"
            ),
            method=(
                "ordinates evaluated at run time from the closed form in "
                "tensorcad.geomquality.baselines.wigley_iii -- no hand-entered numbers"
            ),
            url=_CAD_CODE_URL,
            licence="n/a (mathematical definition)",
            official=False,
            notes=(
                "analytic definition, not an issuing body's tabulated ship table; "
                "used to pin the method because Cb = 4/9, Cwp = 2/3 and VCB = 5T/8 are exact"
            ),
        ),
    ),
    "series60_cb060": HullForm(
        key="series60_cb060",
        specimen_factory=series60,
        provenance=OffsetProvenance(
            source=(
                "Series 60 parent form, Todd (1951) DTMB Report 1447; "
                "Cp = 0.625 and n = 5/3 from the published parent-form table encoded in "
                "tensorcad.geomquality.baselines._SERIES60_PARENTS"
            ),
            method=(
                "ordinates evaluated at run time from "
                "y = (B/2)(1-|2xi-1|^n) zeta^(1/m) in "
                "tensorcad.geomquality.baselines.series60"
            ),
            url=_CAD_CODE_URL,
            licence="n/a (published parametric hull family)",
            official=False,
            notes=(
                "the family's Cwp equals Cp by construction (docs/series60_fidelity.md); "
                "the study therefore judges volume/waterplane against the same half-breadth "
                "function, not against measured Series 60 ship data"
            ),
        ),
    ),
    "kvlcc2_class": HullForm(
        key="kvlcc2_class",
        specimen_factory=kvlcc2_class,
        provenance=OffsetProvenance(
            source=(
                "KVLCC2-class tanker: published particulars L = 320 m, B = 58 m, T = 20.8 m, "
                "Cb = 0.8098, LCB = -3.48 % L (KRISO KVLCC2; ITTC 2010/2015 workshop benchmark)"
            ),
            method=(
                "TensorCAD class hull calibrated to those particulars "
                "(tensorcad.geomquality.baselines.kvlcc2_class): plan form "
                "1-(1-2xi)^a_aft / 1-(2xi-1)^a_fore with a zeta^(1/10) section; "
                "ordinates evaluated at run time, no hand-entered numbers"
            ),
            url=_CAD_CODE_URL,
            licence="n/a (reconstruction of published particulars)",
            official=False,
            notes=(
                "NOT the official KRISO offset table, which is not in this repository "
                "(docs/geomquality-pipeline.md, limitations) and was not reachable from the "
                "build sandbox; Cb/LCB agreement is a representation consistency check"
            ),
        ),
    ),
}


def builtin_forms() -> dict[str, HullForm]:
    """The three published forms this module can tabulate."""
    return dict(_FORMS)


def offset_table(
    key: str,
    *,
    n_stations: int = 21,
    n_waterlines: int = 17,
    waterline_exponent: float = 1.0,
    name: str | None = None,
) -> ShipOffsetTable:
    """Convenience wrapper: ``offset_table('wigley_iii', n_stations=21)``."""
    try:
        form = _FORMS[key]
    except KeyError as error:
        raise OffsetTableError(
            f"unknown hull form {key!r}; available: {sorted(_FORMS)}"
        ) from error
    return form.table(
        n_stations=n_stations,
        n_waterlines=n_waterlines,
        waterline_exponent=waterline_exponent,
        name=name,
    )


def reference_for(key: str) -> dict[str, Any]:
    """Analytic reference integrals (volume, waterplane, LCB, VCB) of a builtin form."""
    try:
        form = _FORMS[key]
    except KeyError as error:
        raise OffsetTableError(
            f"unknown hull form {key!r}; available: {sorted(_FORMS)}"
        ) from error
    return form.reference()