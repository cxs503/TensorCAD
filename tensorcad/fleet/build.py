"""Geometry routes into the baseline fleet: form-family loft, offset table, STL.

Three routes, one contract: each returns a :class:`tensorcad.geometry.Mesh` for a
:class:`~tensorcad.fleet.manifest.ShipRecord`, or raises :class:`RouteRefusal` naming
the measured reason it will not build. A route never silently substitutes a box, a
coarser hull, or a clipped parameter for the geometry it was asked for.

* ``series60_form`` — this repository's Series-60 section law (``tensorcad.series60``),
  calibrated to the record's ``Cb``. Valid only inside the law's validated window
  ``Cb in [0.55, 0.85]``; outside it the route refuses.
* ``offset_table`` — an official offset table (stations x waterline half-breadths),
  lofted and closed at bow/stern. This is the route a downloaded NMRI/DTC/JBC offset
  file would take.
* ``imported_stl`` — a binary or ASCII STL, welded and oriented on the way in, so the
  same gates run on a vendor file as on a generated hull.
"""

from __future__ import annotations

import math
import struct
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

from tensorcad import mesh as meshlib
from tensorcad.geometry import Mesh, MeshBuilder
from tensorcad.fleet.manifest import ShipRecord

#: Window in which the Series-60 section law has been verified (``tensorcad.series60``).
SERIES60_CB_WINDOW = (0.55, 0.85)
SERIES60_UNITS = "m"
SERIES60_FRAME = "x aft->fore, y starboard, z up, waterline at z = 0"
STL_FRAME = "as authored: x,y,z of the file, waterline at z = z_waterline"


class RouteRefusal(ValueError):
    """The requested geometry route cannot represent this ship; the reason is the payload."""


@dataclass(frozen=True)
class BuiltGeometry:
    mesh: Mesh
    route: str
    details: dict[str, Any]

    def as_dict(self) -> dict[str, Any]:
        return {"route": self.route, "details": self.details}


def _series60_mesh(record: ShipRecord, *, nx: int, nz: int) -> BuiltGeometry:
    from tensorcad.series60 import Series60Parameters, build_series60_mesh, series60_section_law

    low, high = SERIES60_CB_WINDOW
    if not low <= record.block_coefficient <= high:
        raise RouteRefusal(
            f"series60_form refuses Cb = {record.block_coefficient}: the Series-60 section law is "
            f"validated only on Cb in [{low}, {high}]; clipping the block coefficient would build a "
            f"different ship and report it as this one"
        )
    parameters = Series60Parameters(
        length=record.length,
        beam=record.beam,
        draft=record.draft,
        nx=nx,
        nz=nz,
        cb=record.block_coefficient,
    )
    errors, warnings = parameters.problems()
    if errors:
        raise RouteRefusal("series60_form parameters rejected: " + "; ".join(errors))
    law = series60_section_law(parameters)
    mesh, offsets = build_series60_mesh(parameters)
    return BuiltGeometry(
        mesh=mesh,
        route="series60_form",
        details={
            "nx": nx,
            "nz": nz,
            "section_law": {
                "n_aft": law.n_aft,
                "n_fore": law.n_fore,
                "cm_mid": law.cm_mid,
                "gain": law.gain,
                "cb": law.cb,
                "cwp": law.cwp,
                "lcb_ratio": law.lcb_ratio,
                "transom_ratio": law.transom_ratio,
            },
            "offsets_points": len(offsets),
            "warnings": warnings,
            "units": SERIES60_UNITS,
            "frame": SERIES60_FRAME,
        },
    )


def loft_station_table(
    stations: Sequence[tuple[float, Sequence[float]]],
    levels: Sequence[float],
    *,
    units: str = "m",
    frame: str = SERIES60_FRAME,
) -> Mesh:
    """Loft half-breadth stations into a closed hull surface.

    ``stations`` is ``(x, (b_0..b_{n-1}))`` with ``b_k`` the half-breadth at
    ``levels[k]``, ordered keel -> deck (levels ascending). Each station becomes a ring
    ``+b_0 .. +b_{n-1}, -b_{n-2} .. -b_1`` so the keel and the deck edge are shared by
    both sides; consecutive rings are joined by quads and the two end stations are
    closed with a fan. A station whose half-breadths are all zero collapses to a single
    point (``MeshBuilder`` deduplicates identical positions), which is how the bow
    closes without a cap. That is what makes a lofted offset table watertight by
    construction instead of by repair.
    """
    levels = [float(z) for z in levels]
    if len(levels) < 2:
        raise RouteRefusal("an offset table needs at least two waterline levels")
    if any(b < 0.0 for _, breadth in stations for b in breadth):
        raise RouteRefusal("half-breadths must be non-negative")
    if any(len(breadth) != len(levels) for _, breadth in stations):
        raise RouteRefusal(
            f"every station must carry one half-breadth per waterline level ({len(levels)} expected)"
        )
    if len(stations) < 2:
        raise RouteRefusal("an offset table needs at least two stations")

    builder = MeshBuilder(units=units, coordinate_frame=frame)
    rings: list[list[int]] = []
    for x, breadth in stations:
        starboard = [builder.add_vertex((x, breadth[k], levels[k])) for k in range(len(levels))]
        last = len(levels) - 1
        # the port side runs deck -> keel. The keel point is shared with the starboard side
        # when the station closes to a point (b_keel == 0) and duplicated when it is flat.
        stop = 1 if float(breadth[0]) == 0.0 else 0
        port = [builder.add_vertex((x, -breadth[k], levels[k])) for k in range(last, stop - 1, -1)]
        rings.append(starboard + port)
    ring_size = len(rings[0])
    for i in range(len(rings) - 1):
        for j in range(ring_size):
            k = (j + 1) % ring_size
            builder.add_quad(rings[i][j], rings[i][k], rings[i + 1][k], rings[i + 1][j])
    for index, outward_sign in ((0, -1.0), (len(rings) - 1, 1.0)):
        _close_station(builder, rings[index], stations[index], outward_sign)
    return meshlib.orient_outward(builder.build()).mesh


def _close_station(builder: MeshBuilder, ring: Sequence[int],
                   station: tuple[float, Sequence[float]], outward_sign: float) -> None:
    """Fan a station ring shut, or leave it collapsed when the station has no breadth."""
    x, breadth = station
    if max(float(b) for b in breadth) <= 0.0:
        return
    levels = _station_levels(ring, builder)
    zs = [p[2] for p in levels]
    centre = builder.add_vertex((x, 0.0, sum(zs) / len(zs)))
    ring = list(ring)
    for j in range(len(ring)):
        k = (j + 1) % len(ring)
        if outward_sign < 0.0:
            builder.add_triangle(centre, ring[j], ring[k])
        else:
            builder.add_triangle(centre, ring[k], ring[j])


def _station_levels(ring: Sequence[int], builder: MeshBuilder) -> list[tuple[float, float, float]]:
    """Recover the positions of a ring's vertices (MeshBuilder stores them privately)."""
    positions = getattr(builder, "_vertices")
    return [positions[index] for index in ring]


def offset_table_mesh(
    record: ShipRecord,
    stations: Sequence[tuple[float, Sequence[float]]],
    levels: Sequence[float],
    *,
    source: str,
) -> BuiltGeometry:
    """Wrap :func:`loft_station_table` for a record and record where the table came from."""
    mesh = loft_station_table(stations, levels)
    return BuiltGeometry(
        mesh=mesh,
        route="offset_table",
        details={
            "source": source,
            "stations": len(stations),
            "levels": len(levels),
            "units": "m",
            "frame": SERIES60_FRAME,
        },
    )


def sample_offset_table(record: ShipRecord, *, stations: int = 41, levels: int = 13
                        ) -> tuple[list[tuple[float, list[float]]], list[float]]:
    """Sample the Series-60 law into offset-table form.

    This is how the ``offset_table`` route is exercised without a downloaded official
    table: it produces a table with exactly the shape an official one has, which lets
    the loft, the closure and the gates be proven, and lets the two routes be
    cross-checked against each other.
    """
    from tensorcad.series60 import Series60Parameters, series60_section_law

    low, high = SERIES60_CB_WINDOW
    if not low <= record.block_coefficient <= high:
        raise RouteRefusal(
            f"offset-table sampling needs the Series-60 law, which refuses Cb = "
            f"{record.block_coefficient} outside [{low}, {high}]"
        )
    parameters = Series60Parameters(
        length=record.length, beam=record.beam, draft=record.draft, nx=stations, nz=levels,
        cb=record.block_coefficient,
    )
    law = series60_section_law(parameters)
    level_z = [-record.draft * (1.0 - k / (levels - 1)) for k in range(levels)]
    table: list[tuple[float, list[float]]] = []
    for i in range(stations):
        xi = i / (stations - 1)
        x = record.length * xi
        taper = law.deck_taper(xi)
        m = law.exponent(xi) if not law.is_legacy else 1.0
        breadth = [0.5 * record.beam * taper * ((k / (levels - 1)) ** (1.0 / m)) for k in range(levels)]
        table.append((x, breadth))
    return table, level_z


def parametric_table(
    record: ShipRecord,
    *,
    stations: int = 481,
    levels: int = 97,
) -> tuple[list[tuple[float, list[float]]], list[float], dict[str, Any]]:
    """Offset table from a generalized hull-surface law solved to the record's particulars.

    The section-area distribution is a Beta shape ``g(xi) = k * xi^a * (1-xi)^b`` whose
    centroid fixes the longitudinal centre of buoyancy and whose normalisation fixes the
    block coefficient exactly::

        int_0^1 g dxi = Cb / Cm          (so the enclosed volume is exactly Cb*L*B*T)
        centroid(g)   = 0.5 + LCB/L      (LCB measured from midship, positive forward)

    The section shape is the power law ``b(z)/ (B/2) = g(xi) * (z/T)^(1/m)`` with
    ``m = Cm/(1-Cm)``, which reproduces the midship section coefficient by construction.
    Both conditions are solved in closed form (Beta normalisation via ``lgamma``), not by
    fitting the mesh, so the gate that compares mesh volume to the record is not
    comparing the mesh to itself.

    This is a *surrogate*: it carries the published ``Cb``, ``Cm`` and ``LCB`` and nothing
    else about the real hull (no bulb, no transom, no local fairing). The record's
    provenance says so, and the run reports ``reference_basis`` next to every number.
    """
    import math as _math

    cb = record.block_coefficient
    cm = record.midship_coefficient if record.midship_coefficient is not None else 0.98
    if not cb < cm <= 1.0:
        raise RouteRefusal(
            f"parametric_sections needs Cb < Cm <= 1, got Cb = {cb}, Cm = {cm}: a prismatic "
            f"coefficient above 1 is not a hull"
        )
    target_integral = cb / cm
    lcb = record.lcb_fraction or 0.0
    centroid = 0.5 + lcb
    if not 0.02 < centroid < 0.98:
        raise RouteRefusal(f"LCB fraction {lcb} puts the centroid at {centroid} of Lpp, outside (0.02, 0.98)")

    # Three conditions, two shape parameters plus a scale: the curve must peak at exactly
    # one midship section (otherwise the mesh cannot span the declared beam), its integral
    # must equal Cb/Cm, and its centroid must sit at the declared LCB. The centroid fixes b
    # in closed form for any a; the integral then fixes a by bisection, and the scale k is
    # pinned by the unit peak.
    a, b, k = _solve_section_area_shape(centroid, target_integral)
    m_vertical = cm / (1.0 - cm)
    x_peak = a / (a + b) if (a + b) > 0 else 0.5
    level_z = [-record.draft * (1.0 - j / (levels - 1)) for j in range(levels)]
    table: list[tuple[float, list[float]]] = []
    peak = 0.0
    for i in range(stations):
        xi = i / (stations - 1)
        if xi <= 0.0 or xi >= 1.0:
            g = 0.0  # fine ends: zero section area, so the hull closes without a cap
        else:
            g = k * (xi ** a) * ((1.0 - xi) ** b)
        peak = max(peak, g)
        breadth = [
            0.5 * record.beam * g * (((j / (levels - 1)) ** (1.0 / m_vertical)) if m_vertical > 0 else 1.0)
            for j in range(levels)
        ]
        table.append((record.length * xi, breadth))
    if peak > 1.0 + 1e-9:
        raise RouteRefusal(
            f"Beta section-area shape peaks at g = {peak:.4f} > 1: the requested Cb/Cm ratio cannot be "
            f"carried by a section-area distribution that stays under B*T*Cm"
        )
    details = {
        "law": "Beta section-area distribution x power-law section shape",
        "cb": cb,
        "cm": cm,
        "lcb_fraction": lcb,
        "beta_a": a,
        "beta_b": b,
        "beta_scale_k": k,
        "vertical_exponent": m_vertical,
        "normalisation_integral": target_integral,
        "centroid_x_over_l": centroid,
        "peak_section_area_ratio": peak,
        "peak_location_x_over_l": x_peak,
        "stations": stations,
        "levels": levels,
    }
    return table, level_z, details


def _beta_integral(a: float, b: float) -> float:
    """``B(a+1, b+1)``: the integral of ``xi^a (1-xi)^b`` over [0, 1]."""
    import math as _math

    return _math.exp(_math.lgamma(a + 1.0) + _math.lgamma(b + 1.0) - _math.lgamma(a + b + 2.0))


def _section_area_shape(centroid: float, a: float) -> tuple[float, float, float]:
    """``(b, k, integral)`` for a section-area curve with a unit peak and this centroid."""
    # centroid of the density proportional to xi^a (1-xi)^b is (a+1)/(a+b+2)
    b = (a + 1.0) / centroid - a - 2.0
    if b <= 1e-9:
        return b, float("nan"), float("nan")
    x_peak = a / (a + b)
    k = 1.0 / ((x_peak ** a) * ((1.0 - x_peak) ** b))
    return b, k, k * _beta_integral(a, b)


def _solve_section_area_shape(centroid: float, target_integral: float) -> tuple[float, float, float]:
    """Solve the section-area curve for the exact peak, integral and centroid demanded.

    The integral falls monotonically as ``a`` grows with the centroid held fixed, so a
    bisection in ``a`` lands on ``Cb/Cm``; ``b`` follows in closed form and ``k`` is pinned
    by ``g(x_peak) = 1``, which is what lets the mesh span the declared beam.
    """
    def integral_of(a: float) -> float:
        value = _section_area_shape(centroid, a)[2]
        return float("inf") if value != value else value

    # b = (a+1)/c - a - 2 must stay positive: that bounds a from below, because a curve whose
    # centroid sits aft of midship cannot be both very full at the stern and peaked amidships.
    a_floor = max(1e-9, (2.0 * centroid - 1.0) / (1.0 - centroid)) * (1.0 + 1e-6) + 1e-9
    low, high = a_floor, 60.0
    probe = integral_of(low)
    if probe != probe or probe < target_integral:
        raise RouteRefusal(
            f"Cb/Cm = {target_integral:.4f} cannot be carried by a section-area curve that peaks at "
            f"one midship section and has its centroid at {centroid:.4f} (max reachable "
            f"{probe if probe == probe else float('nan'):.4f})"
        )
    if integral_of(high) == integral_of(high) and integral_of(high) > target_integral:
        raise RouteRefusal(
            f"Cb/Cm = {target_integral:.4f} is below what a hull-shaped section-area curve can hold "
            f"(min {integral_of(high):.4f})"
        )
    for _ in range(200):
        mid = 0.5 * (low + high)
        value = integral_of(mid)
        if value > target_integral:
            low = mid
        else:
            high = mid
    a = 0.5 * (low + high)
    b, k, integral = _section_area_shape(centroid, a)
    if integral != integral:
        raise RouteRefusal("section-area solve failed to produce a finite shape")
    return a, b, k


def read_stl(path: str | Path) -> Mesh:
    """Read a binary or ASCII STL into a welded mesh.

    Binary is detected from the declared facet count matching the file length, which
    is the only reliable test: many exporters fill the 80-byte header with the word
    "solid", so a text sniff misclassifies them.
    """
    payload = Path(path).read_bytes()
    if len(payload) < 84:
        raise RouteRefusal(f"{path}: {len(payload)} bytes is shorter than an STL header")
    count = struct.unpack("<I", payload[80:84])[0]
    if len(payload) == 84 + 50 * count:
        return _binary_stl(payload, count, str(path))
    text = payload.decode("utf-8", "replace")
    if "facet" in text and "vertex" in text:
        return _ascii_stl(text, str(path))
    raise RouteRefusal(
        f"{path}: neither binary ({count} facets declared, file is {len(payload)} bytes) nor ASCII STL"
    )


def _binary_stl(payload: bytes, count: int, source: str) -> Mesh:
    builder = MeshBuilder(units="m", coordinate_frame=STL_FRAME)
    for index in range(count):
        base = 84 + 50 * index + 12  # skip the facet normal
        points = [struct.unpack_from("<3f", payload, base + 12 * k) for k in range(3)]
        builder.add_triangle(*(builder.add_vertex(point) for point in points))
    return builder.build()


def _ascii_stl(text: str, source: str) -> Mesh:
    builder = MeshBuilder(units="m", coordinate_frame=STL_FRAME)
    current: list[int] = []
    for line in text.splitlines():
        fields = line.split()
        if len(fields) == 4 and fields[0] == "vertex":
            current.append(builder.add_vertex((float(fields[1]), float(fields[2]), float(fields[3]))))
            if len(current) == 3:
                builder.add_triangle(*current)
                current = []
    if not builder.vertex_count:
        raise RouteRefusal(f"{source}: ASCII STL contains no vertices")
    return builder.build()


def stl_mesh(record: ShipRecord, path: str | Path, *, weld_tolerance: float = 1e-9) -> BuiltGeometry:
    """Import an STL through the gate: weld shared vertices, then orient outward."""
    mesh = read_stl(path)
    welded = meshlib.weld(mesh, tolerance=weld_tolerance)
    oriented = meshlib.orient_outward(welded).mesh
    return BuiltGeometry(
        mesh=oriented,
        route="imported_stl",
        details={
            "path": str(path),
            "weld_tolerance": weld_tolerance,
            "triangles_before_weld": len(mesh.faces),
            "triangles_after_weld": len(oriented.faces),
            "units": "m",
            "frame": STL_FRAME,
        },
    )


def build_ship_geometry(
    record: ShipRecord,
    *,
    nx: int = 121,
    nz: int = 25,
    stl_path: str | Path | None = None,
) -> BuiltGeometry:
    """Dispatch on ``record.geometry_route``; raise :class:`RouteRefusal` if unsupported."""
    route = record.geometry_route
    if route == "series60_form":
        return _series60_mesh(record, nx=nx, nz=nz)
    if route == "offset_table":
        if stl_path is not None:
            raise RouteRefusal("offset_table route does not take an STL path")
        table, levels = sample_offset_table(record, stations=nx, levels=nz)
        return offset_table_mesh(record, table, levels, source="series60 law sampled to table form")
    if route == "parametric_sections":
        table, levels, details = parametric_table(record, stations=nx, levels=nz)
        built = offset_table_mesh(record, table, levels,
                                  source="parametric sections solved to published Cb/Cm/LCB")
        built = BuiltGeometry(mesh=built.mesh, route="parametric_sections",
                              details={**built.details, **details})
        return built
    if route == "imported_stl":
        target = stl_path or record.geometry_path
        if target is None:
            raise RouteRefusal("imported_stl route needs a path (record.geometry_path or stl_path)")
        return stl_mesh(record, target)
    raise RouteRefusal(f"unknown geometry route {route!r}")

# --------------------------------------------------------------------------- #
# IGES field route (fourth route)
# --------------------------------------------------------------------------- #
IGES_FRAME = "x aft->fwd, y positive to starboard, z up from keel (m)"


def iges_field_mesh(
    record: ShipRecord,
    iges_path: str | Path,
    *,
    dx: float = 0.5,
    dz: float = 0.1,
    ring_points: int = 240,
    samples_per_metre: float = 8.0,
) -> BuiltGeometry:
    """Build a watertight hull mesh from an IGES file via the half-breadth field route.

    The IGES file is read as B-spline surfaces, sampled into a half-breadth field
    ``b(x, z) = max |y|``, and then lofted into a closed mesh by structured waterline
    slices. This guarantees watertightness by construction because topology is built
    from the field grid, not from patch seams.

    The cost is explicit: the ``max`` binning biases the surface outwards by
    O(cell^2 / radius), which is bounded by the reported ``field_cover`` metric.
    """
    from tensorcad.fleet_baseline import (
        build_solid,
        field_from_points,
        sample_patches,
    )
    from tensorcad.iges import read_iges

    iges_path = Path(iges_path)
    if not iges_path.exists():
        raise RouteRefusal(f"IGES file not found: {iges_path}")

    iges = read_iges(str(iges_path))
    surfaces = [s for s in iges.surfaces if hasattr(s, 'sample')]
    if not surfaces:
        raise RouteRefusal(f"no B-spline surfaces found in {iges_path}")

    spacing = 1.0 / samples_per_metre
    points = sample_patches(surfaces, spacing=spacing)

    field = field_from_points(
        points,
        draft=record.draft,
        dx=dx,
        dz=dz,
        min_x=0.0,
        max_x=record.length,
    )

    mesh, loft_info = build_solid(field, ring_points=ring_points, draft=record.draft)

    return BuiltGeometry(
        mesh=mesh,
        route="iges_field",
        details={
            "source_path": str(iges_path),
            "num_surfaces": len(surfaces),
            "num_sample_points": len(points),
            "field": field.to_dict(),
            "loft": loft_info,
            "dx": dx,
            "dz": dz,
            "ring_points": ring_points,
            "samples_per_metre": samples_per_metre,
            "units": "m",
            "frame": IGES_FRAME,
        },
    )


def build_ship_geometry(
    record: ShipRecord,
    *,
    nx: int = 121,
    nz: int = 25,
    stl_path: str | Path | None = None,
    iges_path: str | Path | None = None,
    dx: float = 0.5,
    dz: float = 0.1,
    ring_points: int = 240,
) -> BuiltGeometry:
    """Dispatch on ``record.geometry_route``; raise :class:`RouteRefusal` if unsupported.

    Parameters
    ----------
    record : ShipRecord
        The ship to build geometry for.
    nx, nz : int
        Resolution for the form-family routes (stations x waterline levels).
    stl_path : str or Path, optional
        Path to an STL file for the ``imported_stl`` route.
    iges_path : str or Path, optional
        Path to an IGES file for the ``iges_field`` route.
    dx, dz : float
        Field cell sizes for the ``iges_field`` route (metres).
    ring_points : int
        Vertices per waterline ring for the ``iges_field`` route.
    """
    route = record.geometry_route
    if route == "series60_form":
        return _series60_mesh(record, nx=nx, nz=nz)
    if route == "offset_table":
        if stl_path is not None:
            raise RouteRefusal("offset_table route does not take an STL path")
        table, levels = sample_offset_table(record, stations=nx, levels=nz)
        return offset_table_mesh(record, table, levels, source="series60 law sampled to table form")
    if route == "parametric_sections":
        table, levels, details = parametric_table(record, stations=nx, levels=nz)
        built = offset_table_mesh(record, table, levels,
                                  source="parametric sections solved to published Cb/Cm/LCB")
        built = BuiltGeometry(mesh=built.mesh, route="parametric_sections",
                              details={**built.details, **details})
        return built
    if route == "imported_stl":
        target = stl_path or record.geometry_path
        if target is None:
            raise RouteRefusal("imported_stl route needs a path (record.geometry_path or stl_path)")
        return stl_mesh(record, target)
    if route == "iges_field":
        target = iges_path or record.geometry_path
        if target is None:
            raise RouteRefusal("iges_field route needs a path (record.geometry_path or iges_path)")
        return iges_field_mesh(record, target, dx=dx, dz=dz, ring_points=ring_points)
    raise RouteRefusal(f"unknown geometry route {route!r}")
