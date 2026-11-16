"""DTMB 5415 offset table from the parametric reconstruction.

Why this exists
---------------
:mod:`tensorcad.offset_forms` provides three analytic hull forms (Wigley III, Series 60,
KVLCC2-class) whose ordinates are evaluated from closed-form half-breadth functions.
Those forms are invaluable for pinning the *method*, but they are not real ships: their
sections are monotone, they have no bulbous bow, and their curvature is everywhere smooth.

The DTMB 5415 is the ITTC's standard validation hull for resistance and seakeeping
(Simoni et al. 2010; ITTC 2010/2015 workshop benchmark).  Its geometry comes from a
9001-configuration design database (Serani et al. 2025, Zenodo) and carries the features
a real ship has: a bulbous bow, a transom stern, flare, and a knuckle at the deck edge.
This module samples the baseline DTMB 5415 (design vector ``v = 0``) into a
:class:`~tensorcad.offset_loft.ShipOffsetTable`, giving the fairing and loft pipelines
a *real-ship* test case that the analytic forms cannot provide.

How the table is built
----------------------
The DTMB 5415 parametric model (:mod:`tensorcad.dtmb5415`) provides a structured
``(patch, girth, station)`` grid of the baseline hull.  Half-breadths are extracted by
interpolating this grid at the requested station positions and reading the starboard
``y`` coordinate at each waterline ``z``.  Where a section is not monotone in ``z``
(the bulbous bow overhangs), the outermost intersection is taken -- the same convention
a naval architect uses when reading an offset table from a lines plan.

The table carries full provenance: the dataset citation, the method (grid interpolation
at ``v = 0``), the SHA-256 of the ordinates, and ``official=False`` because this is a
reconstruction from the parametric model, not the original DTMB 5415 lines plan.

Limitations
-----------
* The table is sampled from a *mesh*, not from an analytic half-breadth function, so
  there is no independent closed-form reference for the integrals.  The hydrostatics
  are compared against the DTMB 5415 module's own volume computation (analytic shell
  volume, not displacement below the waterline).
* The bulbous bow creates non-monotone sections near the fore perpendicular; the
  outermost intersection is taken, which is the standard offset-table convention but
  loses the overhang detail.  A full B-rep loft would be needed to capture it.
* The transom stern is vertical at ``x = 0``; the station at ``x = 0`` has zero
  half-breadth at all waterlines above the transom bottom.

NumPy is imported inside the sampling function (same rule as :mod:`tensorcad.dtmb5415`).
"""
from __future__ import annotations

import math
from typing import Any

from .offset_loft import OffsetProvenance, OffsetTableError, ShipOffsetTable

__all__ = [
    "dtmb5415_offset_table",
    "dtmb5415_provenance",
    "DTMB5415_OFFSET_SCHEMA",
]

DTMB5415_OFFSET_SCHEMA = "tensorcad.dtmb5415_offsets.table/1"


def dtmb5415_provenance() -> OffsetProvenance:
    """Provenance record for the DTMB 5415 offset table."""
    return OffsetProvenance(
        source=(
            "DTMB 5415 surface combatant: ITTC standard validation hull for resistance "
            "and seakeeping (Simonsen et al. 2010; ITTC 2010/2015 workshop benchmark); "
            "Lpp = 142 m, B = 10.2695 m, T = 6.1566 m, Cb ~ 0.511"
        ),
        method=(
            "half-breadths sampled from the baseline (v=0) DTMB 5415 parametric mesh "
            "(tensorcad.dtmb5415) by interpolating the structured station grid at each "
            "requested station; outermost starboard intersection taken at each waterline "
            "(bulbous bow convention)"
        ),
        url=(
            "https://github.com/cxs503/TensorCAD/blob/main/tensorcad/dtmb5415_offsets.py"
        ),
        licence="n/a (reconstruction of published benchmark hull)",
        official=False,
        notes=(
            "NOT the original DTMB 5415 lines plan; sampled from the Zenodo design-space "
            "database parametric reconstruction (Serani et al. 2025); bulbous bow sections "
            "are non-monotone, outermost intersection taken; no independent closed-form "
            "reference for integrals -- compared against the parametric model's own volume"
        ),
    )


def _build_section_lookup(grid):
    """Build a lookup table from the DTMB 5415 structured grid.

    The grid has shape (2 patches, 49 girth, 133 station, 3 coords).
    Returns a list of (patch, station_index, x, y_values, z_values) tuples.
    """
    n_patches = grid.shape[0]
    n_girth = grid.shape[1]
    n_station = grid.shape[2]

    sections = []
    for p in range(n_patches):
        for s in range(n_station):
            x = float(grid[p, 0, s, 0])
            y_vals = [float(grid[p, g, s, 1]) for g in range(n_girth)]
            z_vals = [float(grid[p, g, s, 2]) for g in range(n_girth)]
            sections.append((p, s, x, y_vals, z_vals))

    # Sort by x for efficient lookup
    sections.sort(key=lambda item: item[2])
    return sections


def _interpolate_section_at_x(sections, x_target):
    """Interpolate the section curve (y, z) at longitudinal position x_target.

    Returns (y_values, z_values) for the girth points, or ([], []) if outside range.
    """
    if not sections:
        return [], []

    # Find bracketing sections
    left = None
    right = None
    for i, (p, s, x, y_vals, z_vals) in enumerate(sections):
        if x <= x_target:
            left = (p, s, x, y_vals, z_vals)
        if x >= x_target and right is None:
            right = (p, s, x, y_vals, z_vals)
            break

    if left is None:
        # x_target is before all sections
        return [], []
    if right is None:
        # x_target is after all sections
        return [], []

    _, _, x_left, y_left, z_left = left
    _, _, x_right, y_right, z_right = right

    if x_right == x_left:
        t = 0.0
    else:
        t = (x_target - x_left) / (x_right - x_left)

    n = len(y_left)
    y_out = [(1 - t) * y_left[g] + t * y_right[g] for g in range(n)]
    z_out = [(1 - t) * z_left[g] + t * z_right[g] for g in range(n)]

    return y_out, z_out


def _sample_half_breadths(
    n_stations: int,
    n_waterlines: int,
    waterline_exponent: float,
) -> tuple[tuple[float, ...], tuple[float, ...], tuple[tuple[float, ...], ...]]:
    """Sample half-breadths from the DTMB 5415 baseline structured grid.

    Returns (stations, waterlines, half_breadths).
    """
    import numpy as np

    from .dtmb5415 import DTMB5415Space, DRAFT_M, LPP_M

    space = DTMB5415Space()
    v = np.zeros(27)
    grid = space.station_grid(v)  # shape: (2, 49, 133, 3)

    length = LPP_M
    draft = DRAFT_M

    # Build section lookup
    sections = _build_section_lookup(grid)

    stations = tuple(length * i / (n_stations - 1) for i in range(n_stations))
    waterlines = tuple(
        draft * (j / (n_waterlines - 1)) ** waterline_exponent for j in range(n_waterlines)
    )

    half_breadths: list[tuple[float, ...]] = []

    for x_station in stations:
        # Interpolate section at this station
        y_vals, z_vals = _interpolate_section_at_x(sections, x_station)

        if not y_vals:
            # Station outside the hull -- all zeros
            half_breadths.append(tuple(0.0 for _ in waterlines))
            continue

        # For each waterline, find the outermost (maximum y) intersection
        row: list[float] = []
        for z_wl in waterlines:
            max_y = 0.0
            for g in range(len(y_vals) - 1):
                z0, z1 = z_vals[g], z_vals[g + 1]
                y0, y1 = y_vals[g], y_vals[g + 1]
                # Check if waterline crosses this girth segment
                if (z0 - z_wl) * (z1 - z_wl) <= 0:
                    # Linear interpolation
                    if abs(z1 - z0) > 1e-12:
                        t = (z_wl - z0) / (z1 - z0)
                        y_interp = (1 - t) * y0 + t * y1
                    else:
                        y_interp = max(y0, y1)
                    if y_interp > max_y:
                        max_y = y_interp
            row.append(max(0.0, max_y))
        half_breadths.append(tuple(row))

    return stations, waterlines, tuple(half_breadths)


def dtmb5415_offset_table(
    *,
    n_stations: int = 21,
    n_waterlines: int = 17,
    waterline_exponent: float = 1.0,
    name: str = "dtmb5415",
) -> ShipOffsetTable:
    """Create a ShipOffsetTable from the DTMB 5415 baseline parametric mesh.

    The table is sampled by interpolating the baseline (v=0) structured station grid
    at the requested station positions and reading the outermost starboard half-breadth
    at each waterline.  This is the standard offset-table convention for hulls with
    bulbous bows (where sections can be non-monotone in z).

    Parameters
    ----------
    n_stations : int
        Number of longitudinal stations (>= 3).
    n_waterlines : int
        Number of vertical waterlines (>= 3).
    waterline_exponent : float
        Clustering exponent for waterlines towards the keel.
        ``z_j = T * (j / (n-1)) ** p``.  Default 1.0 (uniform).
    name : str
        Name for the table.

    Returns
    -------
    ShipOffsetTable
        The offset table with full provenance.
    """
    if n_stations < 3:
        raise OffsetTableError(f"n_stations must be >= 3, got {n_stations}")
    if n_waterlines < 3:
        raise OffsetTableError(f"n_waterlines must be >= 3, got {n_waterlines}")

    from .dtmb5415 import DRAFT_M, LPP_M

    stations, waterlines, half_breadths = _sample_half_breadths(
        n_stations, n_waterlines, waterline_exponent
    )

    return ShipOffsetTable(
        name=name,
        length=LPP_M,
        stations=stations,
        waterlines=waterlines,
        half_breadths=half_breadths,
        design_draft=DRAFT_M,
        provenance=dtmb5415_provenance(),
        units="m",
        waterline_exponent=waterline_exponent,
    )
