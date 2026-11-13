"""Ship offset table -> NURBS loft, and the three ways of integrating the result.

Why this exists
---------------
A ship's geometry reaches a design office as an **offset table**: half-breadths ``y``
tabulated on a grid of stations ``x`` (longitudinal) and waterlines ``z`` (vertical).
Everything downstream -- hydrostatics, mesh generation, resistance -- needs a *surface*,
so the table has to be lofted.  ``GAP-CAD-T1-OFFSET-LOFT`` asks for exactly that
pipeline with a machine-decidable acceptance:

* the lofted surface reproduces the table's **displacement** and **waterplane area**
  to better than ``0.5 %`` of the table's own Simpson integral, and
* how that agreement depends on the number of stations is *measured*, not asserted.

What is here
------------
:class:`ShipOffsetTable`
    The table itself: stations, waterlines, half-breadths, design draft, and an
    :class:`OffsetProvenance` record (source, method, URL, licence, ``sha256`` of the
    ordinates).  Reading a table back verifies the digest, so a silently edited CSV
    cannot pass as the table it claims to be.
:func:`table_integrals`
    The **baseline**: sectional areas by composite Simpson over the waterlines, then
    volume/LCB/VCB by composite Simpson over the stations, waterplane area by Simpson
    over ``2 y(x, T)``.  Only the tabulated numbers are used -- no surface.
:func:`loft_offset_table` / :class:`OffsetLoft`
    The loft: every station becomes a ``(y, z)`` section curve and
    :func:`tensorcad.parametric_hull.loft_nurbs_surface` interpolates them into one
    tensor-product NURBS patch (cubic in both directions).  The patch is integrated
    with Gauss-Legendre quadrature per knot span on the surface's own derivatives.
    The half-section area is the Green's-theorem integral

    ``A(u) = integral y dz = integral y(u, v) * dz/dv (u, v) dv``

    (positive for the starboard half-section traversed keel -> deck and closed along the
    centre line), the displacement is ``V = 2 * integral A(u) dx/du du``, and the
    waterplane is ``Aw = 2 * integral y(u, v*(u)) dx/du du`` with ``v*(u)`` the parameter
    where the section crosses the design draft.
:meth:`OffsetLoft.mesh_integrals`
    The third route, and the one the acceptance literally names: tessellate the patch
    (mirrored to a closed body with a waterplane lid) and integrate with
    :func:`tensorcad.parametric_hull.compute_hydrostatics`.
:func:`compare_offset_table` / :func:`station_sensitivity`
    The acceptance: both surface routes against the table baseline and, when a closed
    form exists, against the analytic reference; then the same comparison over a ladder
    of station counts so the convergence is visible instead of assumed.

What is *not* done
------------------
No fairing or smoothing: the loft interpolates the tabulated ordinates, so table noise
stays in the surface (that is the point -- a fit would hide it).  No bulbous bow or
overhanging flare: a station whose section is not monotone in ``z`` is rejected with a
named reason instead of being silently mis-integrated by the waterplane root find.  No
superstructure, appendages, trim or heel: the model is the bare hull at zero trim.  The
integration assumes ``x`` increases with the loft parameter and that ``v = 0`` / ``v = 1``
are the keel and the deck edge of every station.
"""
from __future__ import annotations

import hashlib
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping, Sequence

from .nurbs import NURBSSurface, gauss_legendre
from .parametric_hull import (
    CrossSection,
    compute_hydrostatics,
    loft_nurbs_surface,
    tessellate_surface,
)

__all__ = [
    "LoftComparison",
    "OffsetLoft",
    "OffsetProvenance",
    "OffsetTableError",
    "ShipOffsetTable",
    "TableIntegrals",
    "compare_offset_table",
    "loft_offset_table",
    "read_offset_bundle",
    "read_offset_csv",
    "simpson_integrate",
    "station_sensitivity",
    "table_integrals",
    "write_offset_bundle",
    "write_offset_csv",
]


class OffsetTableError(ValueError):
    """Raised when a table, a bundle or a section cannot be used (named reason)."""


# ---------------------------------------------------------------------------
# composite Simpson
# ---------------------------------------------------------------------------
def _uniform_step(coords: Sequence[float], label: str) -> float:
    if len(coords) < 2:
        raise OffsetTableError(f"{label}: need >= 2 coordinates, got {len(coords)}")
    step = coords[1] - coords[0]
    if step <= 0.0:
        raise OffsetTableError(f"{label}: coordinates must increase, got step {step!r}")
    for a, b in zip(coords, coords[1:]):
        if abs((b - a) - step) > 1e-9 * max(abs(step), 1e-12):
            raise OffsetTableError(
                f"{label}: Simpson needs a uniform grid; {a!r} -> {b!r} is "
                f"{b - a!r} against the first step {step!r}"
            )
    return step


def _monotone_grid(coords: Sequence[float], label: str, *, minimum: int = 3) -> None:
    """Strictly increasing grid of at least ``minimum`` points (uniformity optional)."""
    if len(coords) < minimum:
        raise OffsetTableError(f"{label}: need >= {minimum} points, got {len(coords)}")
    for a, b in zip(coords, coords[1:]):
        if not (b > a):
            raise OffsetTableError(f"{label}: coordinates must increase, got {a!r} -> {b!r}")


def _is_uniform(coords: Sequence[float], *, rtol: float = 1e-9) -> bool:
    step = coords[1] - coords[0]
    for a, b in zip(coords, coords[1:]):
        if abs((b - a) - step) > rtol * max(abs(step), 1e-12):
            return False
    return True


def _poly_mul(a: Sequence[float], b: Sequence[float]) -> list[float]:
    out = [0.0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai == 0.0:
            continue
        for j, bj in enumerate(b):
            out[i + j] += ai * bj
    return out


def _poly_integral(a: Sequence[float], low: float, high: float) -> float:
    return sum(ai * (high ** (i + 1) - low ** (i + 1)) / (i + 1) for i, ai in enumerate(a))


def _quadratic_panel_weights(x0: float, x1: float, x2: float) -> tuple[float, float, float]:
    """Exact integrals of the three Lagrange quadratics on ``[x0, x2]``.

    ``sum_k w_k f_k`` integrates the parabola through ``(x0, f0), (x1, f1), (x2, f2)``;
    the weights are the integrated Lagrange basis in closed form, so the rule is exact
    for quadratics on *any* grid -- the unequal-spacing Simpson rule.
    """
    nodes = (x0, x1, x2)
    weights: list[float] = []
    for skip in range(3):
        poly = [1.0]
        denominator = 1.0
        for j in range(3):
            if j == skip:
                continue
            poly = _poly_mul(poly, [-nodes[j], 1.0])
            denominator *= nodes[skip] - nodes[j]
        weights.append(_poly_integral(poly, x0, x2) / denominator)
    return weights[0], weights[1], weights[2]


def simpson_integrate(values: Sequence[float], coords: Sequence[float]) -> float:
    """Composite Simpson-type integration on a monotone grid.

    * **Uniform grid.** Even interval counts use the 1/3 rule throughout; odd counts use
      one 3/8 panel on the first three intervals and 1/3 for the rest, which keeps the
      order at ``O(h^4)`` without dropping a point (the classic hybrid, e.g.
      Chapra & Canale, *Numerical Methods for Engineers*, section 17.4).
    * **Non-uniform grid** (what a real offset table looks like -- waterlines cluster
      towards the bilge).  Consecutive point *triples* are integrated as the parabola
      through them (the unequal-spacing Simpson rule of :func:`_quadratic_panel_weights`),
      exact for quadratics; an odd number of intervals is rejected with a named reason
      rather than silently degraded to the trapezoid.

    A uniform-grid reduction of the non-uniform branch is checked in the test suite
    (``test_unequal_spacing_rule_matches_classical_simpson``).
    """
    if len(values) != len(coords):
        raise OffsetTableError(
            f"Simpson: {len(values)} values against {len(coords)} coordinates"
        )
    n = len(values) - 1
    if n < 2:
        raise OffsetTableError("Simpson: need >= 3 points")
    f = [float(v) for v in values]
    if _is_uniform(coords):
        h = _uniform_step(coords, "Simpson")
        start = 0
        total = 0.0
        if n % 2:
            total += 3.0 * h / 8.0 * (f[0] + 3.0 * f[1] + 3.0 * f[2] + f[3])
            start = 3
            intervals = n - 3
        else:
            intervals = n
        if intervals >= 2:
            odd = sum(f[start + 1 : start + intervals : 2])
            even = sum(f[start + 2 : start + intervals : 2])
            total += h / 3.0 * (f[start] + 4.0 * odd + 2.0 * even + f[start + intervals])
        return total
    _monotone_grid(coords, "Simpson")
    if n % 2:
        raise OffsetTableError(
            f"unequal-spacing Simpson: {n} intervals is odd; supply an odd number of grid "
            "points so the quadratic panels cover the grid"
        )
    total = 0.0
    for i in range(0, n, 2):
        w0, w1, w2 = _quadratic_panel_weights(coords[i], coords[i + 1], coords[i + 2])
        total += w0 * f[i] + w1 * f[i + 1] + w2 * f[i + 2]
    return total


# ---------------------------------------------------------------------------
# provenance + table
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class OffsetProvenance:
    """Where the ordinates come from, in enough detail to re-derive or contest them.

    ``licence`` and ``url`` are mandatory in spirit: an offset table that cannot be
    traced is not evidence.  ``official`` records whether the ordinates are the issuing
    body's own table (``True``) or a documented reconstruction of the hull (``False``) --
    the two are *not* interchangeable and the comparison report says which is which.
    """

    source: str
    method: str
    url: str = ""
    licence: str = ""
    official: bool = False
    notes: str = ""

    def __post_init__(self) -> None:
        if not self.source.strip():
            raise OffsetTableError("provenance: source citation is empty")
        if not self.method.strip():
            raise OffsetTableError("provenance: method is empty")

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _fmt(value: float) -> str:
    """Deterministic, round-trip-exact float text (``repr`` is shortest round-trip)."""
    return repr(float(value))


@dataclass(frozen=True)
class ShipOffsetTable:
    """Half-breadths on a ``(station, waterline)`` grid, with provenance.

    ``stations`` and ``waterlines`` are metres on an increasing, uniform grid (``x = 0``
    at the aft perpendicular, ``x = length`` at the fore perpendicular; ``z = 0`` at the
    keel).  ``half_breadths[i][j]`` is the starboard half-breadth at station ``i``,
    waterline ``j``; the table is a *function of z* at every station, which is what the
    loft and the waterplane root find assume.
    """

    name: str
    length: float
    stations: tuple[float, ...]
    waterlines: tuple[float, ...]
    half_breadths: tuple[tuple[float, ...], ...]
    design_draft: float
    provenance: OffsetProvenance
    units: str = "m"
    #: clustering exponent ``p`` of the waterline grid, ``z_j = T (j / (n-1))**p``;
    #: ``1.0`` is uniform.  Recorded because the grid is part of the representation.
    waterline_exponent: float = 1.0

    def __post_init__(self) -> None:
        if not self.name.strip():
            raise OffsetTableError("table: name is empty")
        if not (self.length > 0.0):
            raise OffsetTableError(f"table: length must be positive, got {self.length!r}")
        if len(self.stations) < 3:
            raise OffsetTableError(
                f"table: need >= 3 stations for Simpson, got {len(self.stations)}"
            )
        if len(self.waterlines) < 3:
            raise OffsetTableError(
                f"table: need >= 3 waterlines for Simpson, got {len(self.waterlines)}"
            )
        _uniform_step(self.stations, "table stations")
        _monotone_grid(self.waterlines, "table waterlines")
        if not _is_uniform(self.waterlines) and (self.n_waterlines - 1) % 2:
            raise OffsetTableError(
                "table: a clustered (non-uniform) waterline grid needs an odd number of "
                f"points so the unequal-spacing Simpson rule covers it; got {self.n_waterlines}"
            )
        if abs(self.stations[0]) > 1e-9 * self.length or abs(
            self.stations[-1] - self.length
        ) > 1e-9 * self.length:
            raise OffsetTableError(
                "table: stations must span [0, length]; got "
                f"[{self.stations[0]!r}, {self.stations[-1]!r}] against length {self.length!r}"
            )
        if len(self.half_breadths) != len(self.stations):
            raise OffsetTableError(
                f"table: {len(self.half_breadths)} rows against {len(self.stations)} stations"
            )
        for index, row in enumerate(self.half_breadths):
            if len(row) != len(self.waterlines):
                raise OffsetTableError(
                    f"table: row {index} has {len(row)} ordinates against "
                    f"{len(self.waterlines)} waterlines"
                )
            for value in row:
                if not math.isfinite(value) or value < -1e-12:
                    raise OffsetTableError(
                        f"table: row {index} has a non-finite or negative half-breadth {value!r}"
                    )
        if self.design_waterline_index is None:
            raise OffsetTableError(
                f"table: design draft {self.design_draft!r} is not a tabulated waterline "
                f"{self.waterlines!r}; the Simpson waterplane needs that row"
            )

    # -- geometry of the grid ------------------------------------------------
    @property
    def n_stations(self) -> int:
        return len(self.stations)

    @property
    def n_waterlines(self) -> int:
        return len(self.waterlines)

    @property
    def design_waterline_index(self) -> int | None:
        tol = 1e-9 * max(abs(self.design_draft), 1.0)
        for index, z in enumerate(self.waterlines):
            if abs(z - self.design_draft) <= tol:
                return index
        return None

    @property
    def waterline_clustered(self) -> bool:
        return not _is_uniform(self.waterlines)

    @property
    def beam(self) -> float:
        return 2.0 * max(max(row) for row in self.half_breadths)

    @property
    def station_step(self) -> float:
        return self.stations[1] - self.stations[0]

    # -- serialisation -------------------------------------------------------
    def payload(self) -> dict[str, Any]:
        """Canonical, formatting-independent description of the numbers."""
        return {
            "name": self.name,
            "units": self.units,
            "length": self.length,
            "design_draft": self.design_draft,
            "waterline_exponent": self.waterline_exponent,
            "stations": list(self.stations),
            "waterlines": list(self.waterlines),
            "half_breadths": [list(row) for row in self.half_breadths],
        }

    def payload_text(self) -> str:
        return json.dumps(self.payload(), sort_keys=True, separators=(",", ":"))

    def digest(self) -> str:
        """``sha256`` of the canonical payload: the identity of these ordinates."""
        return hashlib.sha256(self.payload_text().encode("utf-8")).hexdigest()

    def to_csv(self) -> str:
        """CSV with a provenance header: an ``x\\z`` header row, then one row per station."""
        lines = [
            f"# tensorcad.offset_loft table {self.name}",
            f"# grid units={self.units} length={_fmt(self.length)} design_draft={_fmt(self.design_draft)} "
            f"waterline_exponent={_fmt(self.waterline_exponent)}",
            f"# provenance {json.dumps(self.provenance.to_dict(), sort_keys=True)}",
            f"# digest {self.digest()}",
            "x\\z," + ",".join(_fmt(z) for z in self.waterlines),
        ]
        for station, row in zip(self.stations, self.half_breadths):
            lines.append(_fmt(station) + "," + ",".join(_fmt(y) for y in row))
        return "\n".join(lines) + "\n"

    def sections(self) -> tuple[CrossSection, ...]:
        """One :class:`CrossSection` per station, ordered keel -> design waterline."""
        limit = self.design_waterline_index
        assert limit is not None  # validated in __post_init__
        drafts = self.waterlines[: limit + 1]
        return tuple(
            CrossSection(x=station, points=tuple((row[j], drafts[j]) for j in range(len(drafts))))
            for station, row in zip(self.stations, self.half_breadths)
        )

    def with_stations(self, indices: Sequence[int], *, name: str | None = None) -> "ShipOffsetTable":
        """Sub-table on the given station indices (used by the station-count study)."""
        stations = tuple(self.stations[i] for i in indices)
        rows = tuple(self.half_breadths[i] for i in indices)
        first = stations[0]
        return ShipOffsetTable(
            name=name or f"{self.name} [{len(indices)} stations]",
            length=stations[-1] - first,
            stations=tuple(x - first for x in stations),
            waterlines=self.waterlines,
            half_breadths=rows,
            design_draft=self.design_draft,
            provenance=self.provenance,
            units=self.units,
            waterline_exponent=self.waterline_exponent,
        )


def write_offset_csv(table: ShipOffsetTable, path: str | Path) -> Path:
    """Write ``table`` to ``path`` (UTF-8, deterministic text)."""
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(table.to_csv(), encoding="utf-8")
    return target


def read_offset_csv(
    path: str | Path,
    *,
    name: str | None = None,
    provenance: OffsetProvenance | None = None,
) -> ShipOffsetTable:
    """Read a table written by :meth:`ShipOffsetTable.to_csv`.

    The ``# provenance`` header is used unless an explicit record is supplied, and a
    ``# digest`` header is *verified*: ordinates edited after writing are rejected.
    """
    text = Path(path).read_text(encoding="utf-8")
    meta: dict[str, str] = {}
    body: list[str] = []
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("#"):
            key, _, rest = stripped[1:].strip().partition(" ")
            meta[key] = rest.strip()
            continue
        body.append(stripped)
    if len(body) < 4:
        raise OffsetTableError(f"{path}: no offset grid found")
    header = body[0].split(",")
    if header[0].strip() != "x\\z":
        raise OffsetTableError(f"{path}: first column header must be 'x\\z', got {header[0]!r}")
    waterlines = tuple(float(v) for v in header[1:])
    stations: list[float] = []
    rows: list[tuple[float, ...]] = []
    for line in body[1:]:
        cells = line.split(",")
        stations.append(float(cells[0]))
        rows.append(tuple(float(v) for v in cells[1:]))
    design_draft = waterlines[-1]
    waterline_exponent = 1.0
    grid = meta.get("grid", "")
    for token in grid.split():
        if token.startswith("design_draft="):
            design_draft = float(token.split("=", 1)[1])
        elif token.startswith("waterline_exponent="):
            waterline_exponent = float(token.split("=", 1)[1])
    if provenance is None:
        if "provenance" not in meta:
            raise OffsetTableError(f"{path}: no '# provenance' header; provenance is required")
        provenance = OffsetProvenance(**json.loads(meta["provenance"]))
    table = ShipOffsetTable(
        name=name or meta.get("tensorcad.offset_loft", Path(path).stem),
        length=stations[-1] - stations[0],
        stations=tuple(stations),
        waterlines=waterlines,
        half_breadths=tuple(rows),
        design_draft=design_draft,
        provenance=provenance,
        waterline_exponent=waterline_exponent,
    )
    recorded = meta.get("digest")
    if recorded is not None and recorded != table.digest():
        raise OffsetTableError(
            f"{path}: digest mismatch -- recorded {recorded}, recomputed {table.digest()}; "
            "the ordinates were edited after the table was written"
        )
    return table


def write_offset_bundle(table: ShipOffsetTable, directory: str | Path) -> tuple[Path, Path]:
    """Write ``<name>.csv`` plus ``<name>.provenance.json`` (digest + source + licence)."""
    folder = Path(directory)
    folder.mkdir(parents=True, exist_ok=True)
    stem = table.name.replace(" ", "_").replace("/", "-")
    csv_path = write_offset_csv(table, folder / f"{stem}.csv")
    record: dict[str, Any] = dict(table.provenance.to_dict())
    record.update(
        {
            "table": table.name,
            "units": table.units,
            "length": table.length,
            "design_draft": table.design_draft,
            "waterline_exponent": table.waterline_exponent,
            "n_stations": table.n_stations,
            "n_waterlines": table.n_waterlines,
            "digest": table.digest(),
            "csv_sha256": hashlib.sha256(csv_path.read_bytes()).hexdigest(),
            "csv": csv_path.name,
        }
    )
    json_path = folder / f"{stem}.provenance.json"
    json_path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return csv_path, json_path


def read_offset_bundle(path: str | Path) -> ShipOffsetTable:
    """Read a ``.csv`` + ``.provenance.json`` pair and verify both digests."""
    csv_path = Path(path)
    if csv_path.suffix != ".csv":
        raise OffsetTableError(f"{path}: expected a .csv file")
    json_path = csv_path.with_suffix(".provenance.json")
    if not json_path.exists():
        raise OffsetTableError(f"{json_path}: provenance sidecar missing")
    record: Mapping[str, Any] = json.loads(json_path.read_text(encoding="utf-8"))
    fields = {
        key: record[key]
        for key in ("source", "method", "url", "licence", "official", "notes")
        if key in record
    }
    table = read_offset_csv(
        csv_path, name=record.get("table"), provenance=OffsetProvenance(**fields)
    )
    if record.get("digest") != table.digest():
        raise OffsetTableError(
            f"{csv_path}: payload digest mismatch -- sidecar {record.get('digest')} "
            f"against recomputed {table.digest()}"
        )
    csv_digest = hashlib.sha256(csv_path.read_bytes()).hexdigest()
    if record.get("csv_sha256") != csv_digest:
        raise OffsetTableError(
            f"{csv_path}: csv_sha256 mismatch -- sidecar {record.get('csv_sha256')} "
            f"against file {csv_digest}"
        )
    return table


# ---------------------------------------------------------------------------
# baseline: Simpson on the table itself
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class TableIntegrals:
    """What the tabulated numbers alone say (composite Simpson, no surface)."""

    volume: float
    waterplane_area: float
    lcb: float
    vcb: float
    block_coefficient: float
    waterplane_coefficient: float
    prismatic_coefficient: float
    midship_coefficient: float
    section_areas: tuple[float, ...]
    method: str = "composite Simpson over waterlines, then over stations"

    def to_dict(self) -> dict[str, Any]:
        data = asdict(self)
        data["section_areas"] = [float(v) for v in self.section_areas]
        return data


def table_integrals(table: ShipOffsetTable) -> TableIntegrals:
    """Displacement, waterplane and coefficients by composite Simpson on the grid."""
    index = table.design_waterline_index
    assert index is not None  # validated in __post_init__
    drafts = table.waterlines[: index + 1]
    areas: list[float] = []
    vertical_moment: list[float] = []
    for row in table.half_breadths:
        half = row[: index + 1]
        areas.append(2.0 * simpson_integrate(half, drafts))
        vertical_moment.append(
            2.0 * simpson_integrate([y * z for y, z in zip(half, drafts)], drafts)
        )
    volume = simpson_integrate(areas, table.stations)
    lcb_moment = simpson_integrate([a * x for a, x in zip(areas, table.stations)], table.stations)
    vcb_moment = simpson_integrate(vertical_moment, table.stations)
    waterplane = simpson_integrate(
        [2.0 * row[index] for row in table.half_breadths], table.stations
    )
    length, beam, draft = table.length, table.beam, table.design_draft
    midship = max(areas)
    return TableIntegrals(
        volume=volume,
        waterplane_area=waterplane,
        lcb=lcb_moment / volume if volume else 0.0,
        vcb=vcb_moment / volume if volume else 0.0,
        block_coefficient=volume / (length * beam * draft),
        waterplane_coefficient=waterplane / (length * beam),
        prismatic_coefficient=volume / (midship * length) if midship else 0.0,
        midship_coefficient=midship / (beam * draft) if beam * draft else 0.0,
        section_areas=tuple(areas),
    )


# ---------------------------------------------------------------------------
# loft
# ---------------------------------------------------------------------------
def _distinct_spans(knots: Sequence[float], degree: int) -> tuple[float, ...]:
    low = knots[degree]
    high = knots[-1 - degree]
    interior = sorted({k for k in knots if low < k < high})
    return tuple([low, *interior, high])


def _panels(spans: Sequence[float], order: int) -> Iterable[tuple[float, float]]:
    nodes, weights = gauss_legendre(order)
    for a, b in zip(spans, spans[1:]):
        half = 0.5 * (b - a)
        mid = 0.5 * (a + b)
        for node, weight in zip(nodes, weights):
            yield mid + half * node, weight * half


def _find_span(knots: Sequence[float], degree: int, u: float) -> int:
    """Index of the knot span containing ``u`` (clamped at both ends)."""
    last = len(knots) - degree - 2
    if u >= knots[last + 1]:
        return last
    if u <= knots[degree]:
        return degree
    low, high = degree, last + 1
    while high - low > 1:
        mid = (low + high) // 2
        if u < knots[mid]:
            high = mid
        else:
            low = mid
    return low


def _basis(knots: Sequence[float], degree: int, span: int, u: float) -> list[float]:
    """The ``degree + 1`` non-zero B-spline basis values (Piegl & Tiller, A2.2)."""
    basis = [0.0] * (degree + 1)
    basis[0] = 1.0
    left = [0.0] * (degree + 1)
    right = [0.0] * (degree + 1)
    for j in range(1, degree + 1):
        left[j] = u - knots[span + 1 - j]
        right[j] = knots[span + j] - u
        saved = 0.0
        for r in range(j):
            denominator = right[r + 1] + left[j - r]
            temp = basis[r] / denominator if denominator else 0.0
            basis[r] = saved + right[r + 1] * temp
            saved = left[j - r] * temp
        basis[j] = saved
    return basis


def _basis_row(knots: Sequence[float], degree: int, u: float) -> list[float]:
    """All basis values of ``degree`` at ``u`` (length ``len(knots) - degree - 1``)."""
    count = len(knots) - degree - 1
    row = [0.0] * count
    if u >= knots[count]:
        row[count - 1] = 1.0
        return row
    span = _find_span(knots, degree, u)
    for index, value in enumerate(_basis(knots, degree, span, u)):
        row[span - degree + index] = value
    return row


def _basis_derivative_row(knots: Sequence[float], degree: int, u: float) -> list[float]:
    """First derivatives of all basis functions at ``u``.

    Uses ``dN_i/du = q (N_i^{q-1}/(U_{i+q}-U_i) - N_{i+1}^{q-1}/(U_{i+q+1}-U_{i+1}))``:
    the derivative of a degree-``q`` spline is a degree-``q-1`` spline on the same knot
    vector, so the derivative is exact and no separate derivative recursion is needed.
    """
    count = len(knots) - degree - 1
    if degree < 1:
        return [0.0] * count
    lower = _basis_row(knots, degree - 1, u)
    out = [0.0] * count
    for i in range(count):
        a = knots[i + degree] - knots[i]
        b = knots[i + degree + 1] - knots[i + 1]
        term_a = lower[i] / a if a else 0.0
        term_b = lower[i + 1] / b if b else 0.0
        out[i] = degree * (term_a - term_b)
    return out


@dataclass
class OffsetLoft:
    """A NURBS patch lofted through the table's stations, plus its integrals.

    The integrals do not call :meth:`NURBSSurface.derivatives` point by point (that is
    the honest but slow route -- ~60000 tensor evaluations per hull).  The
    tensor-product structure is used instead: with unit weights the patch is

    ``S(u, v) = sum_i sum_j B_i(u) B_j(v) P_ij``

    so the section at ``u`` is the ``v``-spline whose control points are
    ``c_j(u) = sum_i B_i(u) P_ij``, the half-section area is the bilinear form
    ``A(u) = sum_jk c_j^y c_k^z M1[j][k]`` with ``M1[j][k] = integral B_j B_k' dv``
    precomputed once over the ``v`` knot spans, and the displacement is
    ``V = 2 * integral A(u) dx/du du`` with ``dx/du`` from the exact derivative basis.
    Every factor is a polynomial on each knot span, so Gauss-Legendre with ``order`` 12
    per span is exact to round-off.  :meth:`volume_direct` keeps the slow
    ``NURBSSurface.derivatives`` route as an independent check of the fast path.
    """

    table: ShipOffsetTable
    surface: NURBSSurface
    degree_u: int
    degree_v: int

    def __post_init__(self) -> None:
        if self.surface.rational:
            raise OffsetTableError(
                "OffsetLoft: the fast integral path assumes unit weights; this surface is "
                "rational -- use volume_direct()/waterplane_area_direct() instead"
            )
        self._net = self.surface.control_net
        self._n_u = len(self._net)
        self._n_v = len(self._net[0])
        self._knots_u, self._knots_v = self.surface.knots
        self._p_u, self._p_v = self.surface.degree
        self._u_spans = _distinct_spans(self._knots_u, self._p_u)
        self._v_spans = _distinct_spans(self._knots_v, self._p_v)
        self._check_x_rows()
        self._m1: list[list[float]] | None = None
        self._m3: dict[int, list[list[list[float]]]] = {}

    # -- helpers -------------------------------------------------------------
    @property
    def u_spans(self) -> tuple[float, ...]:
        return self._u_spans

    @property
    def v_spans(self) -> tuple[float, ...]:
        return self._v_spans

    def _check_x_rows(self) -> float:
        """``x`` must be constant along ``v`` (each station is a planar section)."""
        spread = 0.0
        for row in self._net:
            xs = [point[0] for point in row]
            spread = max(spread, max(xs) - min(xs))
        if spread > 1e-9 * max(self.table.length, 1.0):
            raise OffsetTableError(
                f"OffsetLoft: x varies along the section by {spread!r} m; the station-based "
                "integrals assume planar sections"
            )
        return spread

    def _control_at(
        self, u: float
    ) -> tuple[list[tuple[float, float, float]], float, float]:
        """Section control points ``c_j(u)``, ``x(u)`` and ``dx/du`` at ``u``."""
        weights = _basis_row(self._knots_u, self._p_u, u)
        derivatives = _basis_derivative_row(self._knots_u, self._p_u, u)
        ctrl: list[tuple[float, float, float]] = []
        for j in range(self._n_v):
            x = y = z = 0.0
            for i, w in enumerate(weights):
                if w == 0.0:
                    continue
                point = self._net[i][j]
                x += w * point[0]
                y += w * point[1]
                z += w * point[2]
            ctrl.append((x, y, z))
        x_u = 0.0
        for i, d in enumerate(derivatives):
            if d:
                x_u += d * self._net[i][0][0]
        return ctrl, (ctrl[0][0] if ctrl else 0.0), x_u

    def _section_point(
        self, ctrl: Sequence[tuple[float, float, float]], v: float
    ) -> tuple[float, float, float]:
        weights = _basis_row(self._knots_v, self._p_v, v)
        x = y = z = 0.0
        for w, point in zip(weights, ctrl):
            if w == 0.0:
                continue
            x += w * point[0]
            y += w * point[1]
            z += w * point[2]
        return x, y, z

    def _v_matrices(self, order: int) -> list[list[float]]:
        """``M1[j][k] = integral B_j B_k' dv`` over the whole ``v`` domain."""
        if self._m1 is not None:
            return self._m1
        matrices = [[0.0] * self._n_v for _ in range(self._n_v)]
        for v, weight in _panels(self._v_spans, order):
            basis = _basis_row(self._knots_v, self._p_v, v)
            derivative = _basis_derivative_row(self._knots_v, self._p_v, v)
            for j, bj in enumerate(basis):
                if bj == 0.0:
                    continue
                row = matrices[j]
                for k, dk in enumerate(derivative):
                    if dk:
                        row[k] += weight * bj * dk
        self._m1 = matrices
        return matrices

    def _v_matrix_cubic(self, order: int) -> list[list[list[float]]]:
        """``M3[j][k][l] = integral B_j B_k B_l' dv`` (needed only for the VCB)."""
        cached = self._m3.get(order)
        if cached is not None:
            return cached
        count = self._n_v
        matrices = [[[0.0] * count for _ in range(count)] for _ in range(count)]
        for v, weight in _panels(self._v_spans, order):
            basis = _basis_row(self._knots_v, self._p_v, v)
            derivative = _basis_derivative_row(self._knots_v, self._p_v, v)
            for j, bj in enumerate(basis):
                if bj == 0.0:
                    continue
                plane = matrices[j]
                for k, bk in enumerate(basis):
                    if bk == 0.0:
                        continue
                    row = plane[k]
                    factor = weight * bj * bk
                    for l, dl in enumerate(derivative):
                        if dl:
                            row[l] += factor * dl
        self._m3[order] = matrices
        return matrices

    def u_at_x(self, target: float) -> float:
        """Loft parameter ``u`` of the station at ``x = target`` (invert the monotone ``x(u)``).

        The loft's ``u`` parametrisation is a mean chord-length, not a uniform station
        index, so section-wise queries must invert ``x(u)`` instead of assuming
        ``u = i / (n - 1)`` -- that assumption read a 5.9 m ordinate error on the
        full-scale tanker before it was removed.
        """
        low, high = self._u_spans[0], self._u_spans[-1]
        x_low, x_high = self._control_at(low)[1], self._control_at(high)[1]
        if not (x_low <= target <= x_high):
            raise OffsetTableError(
                f"station x={target!r} outside the lofted range [{x_low!r}, {x_high!r}]"
            )
        if abs(target - x_low) <= 1e-12 * max(self.table.length, 1.0):
            return low
        if abs(target - x_high) <= 1e-12 * max(self.table.length, 1.0):
            return high
        previous = self._control_at(low)[1]
        for _ in range(200):
            mid = 0.5 * (low + high)
            value = self._control_at(mid)[1]
            if abs(value - target) < 1e-14 * max(self.table.length, 1.0):
                return mid
            if value < target:
                low, previous = mid, value
            else:
                high = mid
            if high - low <= 1e-16:
                break
        return 0.5 * (low + high)

    def _v_at_height(self, u: float, height: float, *, tolerance: float = 1e-13) -> float:
        """Parameter ``v`` where the section at ``u`` crosses ``z = height``.

        Bisection on a *checked* monotone branch: the section is sampled first and a
        non-monotone ``z(v)`` (a bulbous bow, an overhanging flare) is rejected with a
        named reason instead of being silently answered with an endpoint.
        """
        ctrl, _, _ = self._control_at(u)
        low, high = self._v_spans[0], self._v_spans[-1]
        z_low = self._section_point(ctrl, low)[2]
        z_high = self._section_point(ctrl, high)[2]
        if height < z_low - tolerance or height > z_high + tolerance:
            raise OffsetTableError(
                f"waterplane {height!r} outside the section range [{z_low!r}, {z_high!r}] "
                f"at u={u!r}"
            )
        probes = [low + (high - low) * i / 32.0 for i in range(33)]
        samples = [self._section_point(ctrl, v)[2] for v in probes]
        scale = max(abs(z_high - z_low), 1e-12)
        for a, b in zip(samples, samples[1:]):
            if b - a < -1e-6 * scale:
                raise OffsetTableError(
                    f"section at u={u!r} is not monotone in z (z falls from {a!r} to {b!r}); "
                    "the waterplane root find needs a single-valued station curve"
                )
        if abs(z_low - height) <= tolerance:
            return low
        if abs(z_high - height) <= tolerance:
            return high
        for _ in range(200):
            mid = 0.5 * (low + high)
            if self._section_point(ctrl, mid)[2] < height:
                low = mid
            else:
                high = mid
            if high - low <= 1e-15 * max(abs(high), 1e-15):
                break
        return 0.5 * (low + high)

    # -- integrals -----------------------------------------------------------
    def _section_integrals(self, u: float, order: int) -> tuple[float, float]:
        """Half-section area and vertical moment ``integral z y dz`` at ``u``."""
        ctrl, _, _ = self._control_at(u)
        m1 = self._v_matrices(order)
        area = 0.0
        for j, cj in enumerate(ctrl):
            if cj[1] == 0.0:
                continue
            row = m1[j]
            for k, ck in enumerate(ctrl):
                if ck[2]:
                    area += cj[1] * ck[2] * row[k]
        m3 = self._v_matrix_cubic(order)
        moment = 0.0
        for j, cj in enumerate(ctrl):
            if cj[2] == 0.0:
                continue
            plane = m3[j]
            for k, ck in enumerate(ctrl):
                if ck[1] == 0.0:
                    continue
                row = plane[k]
                for l, cl in enumerate(ctrl):
                    if cl[2]:
                        moment += cj[2] * ck[1] * cl[2] * row[l]
        return area, moment

    def section_half_area(self, u: float, *, order: int = 12) -> float:
        """Starboard half-section area ``integral y dz`` at loft parameter ``u``."""
        return self._section_integrals(u, order)[0]

    def section_area(self, u: float, *, order: int = 12) -> float:
        return 2.0 * self.section_half_area(u, order=order)

    def section_areas_at_stations(self, *, order: int = 12) -> tuple[float, ...]:
        return tuple(
            self.section_area(self.u_at_x(x), order=order) for x in self.table.stations
        )

    def volume(self, *, order: int = 12) -> float:
        """Displacement ``2 * integral A(u) dx/du du`` (exact per knot span)."""
        total = 0.0
        for u, weight in _panels(self._u_spans, order):
            _, _, x_u = self._control_at(u)
            area, _ = self._section_integrals(u, order)
            total += weight * area * x_u
        return 2.0 * total

    def lcb(self, *, order: int = 12) -> float:
        volume_moment = 0.0
        volume = 0.0
        for u, weight in _panels(self._u_spans, order):
            _, x, x_u = self._control_at(u)
            area, _ = self._section_integrals(u, order)
            volume_moment += weight * x * area * x_u
            volume += weight * area * x_u
        return volume_moment / volume if volume else 0.0

    def vcb(self, *, order: int = 12) -> float:
        moment = 0.0
        volume = 0.0
        for u, weight in _panels(self._u_spans, order):
            _, _, x_u = self._control_at(u)
            area, z_moment = self._section_integrals(u, order)
            moment += weight * z_moment * x_u
            volume += weight * area * x_u
        return moment / volume if volume else 0.0

    def waterplane_half_breadth(self, u: float) -> float:
        ctrl, _, _ = self._control_at(u)
        v = self._v_at_height(u, self.table.design_draft)
        return self._section_point(ctrl, v)[1]

    def waterplane_area(self, *, order: int = 12) -> float:
        total = 0.0
        height = self.table.design_draft
        for u, weight in _panels(self._u_spans, order):
            ctrl, _, x_u = self._control_at(u)
            v = self._v_at_height(u, height)
            total += weight * 2.0 * self._section_point(ctrl, v)[1] * abs(x_u)
        return total

    def wetted_surface_area(self, *, order: int = 16) -> float:
        return 2.0 * self.surface.area(order)

    # -- independent slow routes (the tests use these to check the fast path) -
    def section_half_area_direct(self, u: float, *, order: int = 12) -> float:
        """``integral y dz/dv dv`` from :meth:`NURBSSurface.derivatives` (slow reference)."""
        total = 0.0
        for v, weight in _panels(self._v_spans, order):
            derivatives = self.surface.derivatives(u, v, 1)
            point, dv = derivatives[(0, 0)], derivatives[(0, 1)]
            total += weight * point[1] * dv[2]
        return total

    def volume_direct(self, *, order: int = 12) -> float:
        """``2 * integral integral y dz/dv dx/du`` from surface derivatives (slow)."""
        half = 0.0
        for u, weight_u in _panels(self._u_spans, order):
            for v, weight_v in _panels(self._v_spans, order):
                derivatives = self.surface.derivatives(u, v, 1)
                point, du, dv = derivatives[(0, 0)], derivatives[(1, 0)], derivatives[(0, 1)]
                half += weight_u * weight_v * point[1] * dv[2] * du[0]
        return 2.0 * half

    def waterplane_area_direct(self, *, order: int = 12) -> float:
        total = 0.0
        height = self.table.design_draft
        for u, weight in _panels(self._u_spans, order):
            v = self._v_at_height(u, height)
            derivatives = self.surface.derivatives(u, v, 1)
            point, du = derivatives[(0, 0)], derivatives[(1, 0)]
            total += weight * 2.0 * point[1] * abs(du[0])
        return total

    # -- mesh route ----------------------------------------------------------
    def mesh_integrals(self, *, n_u: int = 97, n_v: int = 33) -> dict[str, float]:
        """Tessellate the patch and integrate the closed body (the acceptance route)."""
        tessellated = tessellate_surface(
            self.surface, n_u=n_u, n_v=n_v, mirror_y=True, close_waterplane=True
        )
        result = compute_hydrostatics(
            tessellated.mesh,
            waterline_z=self.table.design_draft,
            length=self.table.length,
            beam=self.table.beam,
            draft=self.table.design_draft,
            n_hull_faces=tessellated.hull_face_count,
        )
        return {
            "volume": result.displaced_volume,
            "waterplane_area": result.waterplane_area,
            "lcb": result.lcb,
            "vcb": result.vcb,
            "wetted_surface_area": result.wetted_surface_area,
            "n_u": float(n_u),
            "n_v": float(n_v),
        }

    # -- deviation against the table ----------------------------------------
    def station_deviation(self) -> float:
        """Max ``|y_loft - y_table|`` at the tabulated ordinates (should be round-off)."""
        worst = 0.0
        limit = self.table.design_waterline_index
        assert limit is not None
        for station, row in zip(self.table.stations, self.table.half_breadths):
            u = self.u_at_x(station)
            ctrl, _, _ = self._control_at(u)
            for j in range(limit + 1):
                v = self._v_at_height(u, self.table.waterlines[j])
                worst = max(worst, abs(self._section_point(ctrl, v)[1] - row[j]))
        return worst


def loft_offset_table(
    table: ShipOffsetTable, *, degree_u: int = 3, degree_v: int = 3
) -> OffsetLoft:
    """Loft the table's stations into one interpolating NURBS patch."""
    surface = loft_nurbs_surface(
        table.sections(), degree_u=degree_u, degree_v=degree_v, interpolate=True
    )
    return OffsetLoft(
        table=table,
        surface=surface,
        degree_u=surface.degree[0],
        degree_v=surface.degree[1],
    )


# ---------------------------------------------------------------------------
# the acceptance
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class LoftComparison:
    """Table baseline, NURBS quadrature and loft mesh, side by side."""

    name: str
    tolerance: float
    volume_table: float
    volume_nurbs: float
    volume_mesh: float
    waterplane_table: float
    waterplane_nurbs: float
    waterplane_mesh: float
    lcb_table: float
    lcb_nurbs: float
    vcb_table: float
    vcb_nurbs: float
    station_deviation: float
    reference: dict[str, float] = None  # type: ignore[assignment]
    degree_u: int = 3
    degree_v: int = 3

    def __post_init__(self) -> None:
        if self.reference is None:
            object.__setattr__(self, "reference", {})

    @staticmethod
    def _rel(value: float, base: float) -> float:
        return abs(value - base) / abs(base) if base else float("inf")

    @property
    def volume_error_nurbs(self) -> float:
        return self._rel(self.volume_nurbs, self.volume_table)

    @property
    def volume_error_mesh(self) -> float:
        return self._rel(self.volume_mesh, self.volume_table)

    @property
    def waterplane_error_nurbs(self) -> float:
        return self._rel(self.waterplane_nurbs, self.waterplane_table)

    @property
    def waterplane_error_mesh(self) -> float:
        return self._rel(self.waterplane_mesh, self.waterplane_table)

    @property
    def volume_error_reference(self) -> float | None:
        value = self.reference.get("volume")
        return None if value is None else self._rel(self.volume_nurbs, value)

    @property
    def waterplane_error_reference(self) -> float | None:
        value = self.reference.get("waterplane_area")
        return None if value is None else self._rel(self.waterplane_nurbs, value)

    @property
    def passed(self) -> bool:
        return (
            self.volume_error_nurbs <= self.tolerance
            and self.volume_error_mesh <= self.tolerance
            and self.waterplane_error_nurbs <= self.tolerance
            and self.waterplane_error_mesh <= self.tolerance
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "tolerance": self.tolerance,
            "degrees": [self.degree_u, self.degree_v],
            "volume": {
                "table_simpson": self.volume_table,
                "nurbs_quadrature": self.volume_nurbs,
                "loft_mesh": self.volume_mesh,
                "relative_error_nurbs": self.volume_error_nurbs,
                "relative_error_mesh": self.volume_error_mesh,
                "relative_error_nurbs_vs_reference": self.volume_error_reference,
            },
            "waterplane": {
                "table_simpson": self.waterplane_table,
                "nurbs_quadrature": self.waterplane_nurbs,
                "loft_mesh": self.waterplane_mesh,
                "relative_error_nurbs": self.waterplane_error_nurbs,
                "relative_error_mesh": self.waterplane_error_mesh,
                "relative_error_nurbs_vs_reference": self.waterplane_error_reference,
            },
            "lcb": {
                "table_simpson": self.lcb_table,
                "nurbs_quadrature": self.lcb_nurbs,
            },
            "vcb": {
                "table_simpson": self.vcb_table,
                "nurbs_quadrature": self.vcb_nurbs,
            },
            "station_deviation": self.station_deviation,
            "reference": dict(self.reference),
            "passed": self.passed,
        }


def compare_offset_table(
    table: ShipOffsetTable,
    *,
    degree_u: int = 3,
    degree_v: int = 3,
    order: int = 12,
    tolerance: float = 5e-3,
    reference: Mapping[str, float] | None = None,
    mesh: bool = True,
    n_u: int = 97,
    n_v: int = 33,
) -> LoftComparison:
    """Compare the lofted surface against the table's own Simpson integrals."""
    baseline = table_integrals(table)
    loft = loft_offset_table(table, degree_u=degree_u, degree_v=degree_v)
    volume_nurbs = loft.volume(order=order)
    waterplane_nurbs = loft.waterplane_area(order=order)
    if mesh:
        mesh_result = loft.mesh_integrals(n_u=n_u, n_v=n_v)
        volume_mesh = mesh_result["volume"]
        waterplane_mesh = mesh_result["waterplane_area"]
    else:
        volume_mesh = volume_nurbs
        waterplane_mesh = waterplane_nurbs
    return LoftComparison(
        name=table.name,
        tolerance=tolerance,
        volume_table=baseline.volume,
        volume_nurbs=volume_nurbs,
        volume_mesh=volume_mesh,
        waterplane_table=baseline.waterplane_area,
        waterplane_nurbs=waterplane_nurbs,
        waterplane_mesh=waterplane_mesh,
        lcb_table=baseline.lcb,
        lcb_nurbs=loft.lcb(order=order),
        vcb_table=baseline.vcb,
        vcb_nurbs=loft.vcb(order=order),
        station_deviation=loft.station_deviation(),
        reference=dict(reference or {}),
        degree_u=loft.degree_u,
        degree_v=loft.degree_v,
    )


def station_sensitivity(
    build: Callable[[int], ShipOffsetTable],
    counts: Sequence[int],
    *,
    reference: Mapping[str, float] | None = None,
    denominator: TableIntegrals | None = None,
    degree_u: int = 3,
    degree_v: int = 3,
) -> tuple[dict[str, Any], ...]:
    """Loft the same hull at ``counts`` station counts and tabulate the convergence.

    ``build(n)`` must return the table sampled at ``n`` stations on the *same* hull, so
    the only variable is the longitudinal resolution.  Each row reports the loft's
    volume and waterplane area against (a) Simpson on that same table (the discretisation
    error of the sampling itself) and (b) ``denominator`` -- normally the finest table in
    the study, which is the best available stand-in for the hull -- plus (c) the closed
    form ``reference`` when the hull has one.
    """
    rows: list[dict[str, Any]] = []
    for n in counts:
        table = build(n)
        baseline = table_integrals(table)
        loft = loft_offset_table(table, degree_u=degree_u, degree_v=degree_v)
        volume = loft.volume()
        waterplane = loft.waterplane_area()
        row: dict[str, Any] = {
            "n_stations": table.n_stations,
            "station_step": table.station_step,
            "volume_table_simpson": baseline.volume,
            "volume_nurbs": volume,
            "waterplane_table_simpson": baseline.waterplane_area,
            "waterplane_nurbs": waterplane,
            "relative_error_volume_vs_own_table": abs(volume - baseline.volume) / baseline.volume,
            "relative_error_waterplane_vs_own_table": abs(
                waterplane - baseline.waterplane_area
            )
            / baseline.waterplane_area,
        }
        if denominator is not None:
            row["relative_error_volume_vs_denominator"] = (
                abs(volume - denominator.volume) / denominator.volume
            )
            row["relative_error_waterplane_vs_denominator"] = (
                abs(waterplane - denominator.waterplane_area) / denominator.waterplane_area
            )
        if reference:
            if "volume" in reference:
                row["relative_error_volume_vs_reference"] = (
                    abs(volume - reference["volume"]) / reference["volume"]
                )
            if "waterplane_area" in reference:
                row["relative_error_waterplane_vs_reference"] = (
                    abs(waterplane - reference["waterplane_area"])
                    / reference["waterplane_area"]
                )
        rows.append(row)
    return tuple(rows)
