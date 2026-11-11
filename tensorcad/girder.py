"""Hull-girder section properties: what a cross-section costs to bend.

Every other module in this package either *makes* a hull (:mod:`tensorcad.hull`,
:mod:`tensorcad.series60`, :mod:`tensorcad.wigley`) or asks what the hull does to the
water (:mod:`tensorcad.hydrostatics`, :mod:`tensorcad.stability`,
:mod:`tensorcad.resistance`). This one asks the structural question about the same
geometry: **given the transverse section that the shell encloses, how much bending
moment can it carry before the deck yields, and where does the beam first give way?**

It is the first structural module in the package, and it deliberately works on the
one input a geometry kernel actually owns: a simple closed polygon in the section
plane, ``(y, z)`` = (transverse, vertical) in metres. From that polygon it returns the
quantities a hull-girder check is written in:

* :class:`GirderSection` -- area, centroid, second moments about the centroidal axes,
  the product of inertia and the principal axes, the four elastic section moduli
  (deck, keel, port, starboard), the plastic neutral axis and the plastic section
  moduli, radii of gyration, extents and perimeter.
* :func:`girder_section` -- build one from any closed polygon.
* :func:`mirror_half_section` -- build the closed section from a half-breadth table
  ``(z, b(z))``, which is what an offset table or a loft gives you.
* :func:`section_modulus_extremes` -- the pair ``(I / z_max, I / z_min)`` a scanning
  loop actually wants.
* :func:`closed_cell_torsion_constant` / :func:`open_section_torsion_constant` -- the
  thin-walled torsion constant, Bredt for a single cell and the ``sum b t^3 / 3``
  rule for a strip assembly, which is the hull girder's other degree of freedom.

All of it is exact for a polygon: the area, first moments and second moments come from
the shoelace/Euler sums, so the only error is floating-point round-off (measured in
``docs/girder_section.md``: ``6.7e-16`` relative on a rectangle, ``1.1e-15`` on an
I-beam). Nothing is discretised and there is no convergence to study -- unlike the
mesh integrals in :mod:`tensorcad.hydrostatics`, this module never approximates the
geometry it is given. The plastic modulus is the one quantity that needs a search (the
neutral axis that halves the area); it is located by bisection on the exact clipped
area and converges to round-off in a few dozen steps.

What this module is **not**, stated so a consumer does not mistake a section property
for a strength result:

* Not a stress calculation. No bending moment is applied and no material is named; the
  section modulus is the *denominator* of ``sigma = M / Z``, and the numerator is the
  caller's load, which belongs to a rule set or a seakeeping calculation.
* Not an effective-breadth / shear-lag model. The section is taken as fully effective:
  a wide deck in a slender ship is not fully effective at the ends, and the standard
  corrections (Schade, von Karman) are outside the polygon's knowledge.
* Not a buckling or ultimate-strength check. The plastic modulus is the full-plastic
  first moment, not a collapse moment with plate slenderness and residual stress; ice
  class rules and CSR ultimate strength need those knobs and they are not invented here.
* Not a class-rule compliance engine. It reports numbers; comparing them to a rule
  minimum is a policy decision and lives with the rule, not with the geometry.
* Thin-walled torsion is the classical approximation, valid for ``t << depth``; it is
  documented as such and tested against its own closed form, not against 3-D elasticity.

The unit contract is metres and square metres throughout, and the sign contract is the
standard ship convention: ``y`` positive to starboard, ``z`` positive upward, the
polygon's orientation is normalised internally so both a clockwise and a
counter-clockwise offset loop give the same answer.
"""
from __future__ import annotations

import math
from collections.abc import Iterable, Sequence
from dataclasses import dataclass

Point = tuple[float, float]

#: Polygons longer than this skip the O(n^2) simplicity check; the check exists to
#: catch a mis-ordered offset table, not to police a 10 000-vertex shell section.
MAX_VERTICES_FOR_SIMPLICITY_CHECK = 4000

#: Bisection stops when the equal-area level moves by less than this fraction of the
#: section depth. 60 halvings of the initial bracket put it at round-off for a double.
_PLASTIC_TOLERANCE_RELATIVE = 1e-15
_PLASTIC_MAX_ITERATIONS = 200

#: The clipped-area consistency guard: below + above must rebuild the section area.
_SPLIT_TOLERANCE_RELATIVE = 1e-9


class GirderError(ValueError):
    """Raised when a section cannot be read as a simple closed polygon."""


# ---------------------------------------------------------------------------
# polygon primitives
# ---------------------------------------------------------------------------


def _orientation(a: Point, b: Point, c: Point) -> float:
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def _segments_cross(p1: Point, p2: Point, p3: Point, p4: Point) -> bool:
    d1 = _orientation(p3, p4, p1)
    d2 = _orientation(p3, p4, p2)
    d3 = _orientation(p1, p2, p3)
    d4 = _orientation(p1, p2, p4)
    return ((d1 > 0 > d2) or (d1 < 0 < d2)) and ((d3 > 0 > d4) or (d3 < 0 < d4))


def _clean_polygon(points: Iterable[Sequence[float]]) -> list[Point]:
    """Validate, de-duplicate, orient counter-clockwise and return the vertex list."""
    raw: list[Point] = []
    for index, point in enumerate(points):
        try:
            y, z = float(point[0]), float(point[1])
        except (TypeError, IndexError) as error:  # pragma: no cover - defensive
            raise GirderError(f"vertex {index} is not a (y, z) pair: {point!r}") from error
        if not (math.isfinite(y) and math.isfinite(z)):
            raise GirderError(f"vertex {index} is not finite: {point!r}")
        raw.append((y, z))

    cleaned: list[Point] = []
    for point in raw:
        if cleaned and cleaned[-1] == point:
            continue
        cleaned.append(point)
    if len(cleaned) > 1 and cleaned[0] == cleaned[-1]:
        cleaned.pop()
    if len(cleaned) < 3:
        raise GirderError(
            f"a section needs at least 3 distinct vertices, got {len(cleaned)}"
        )

    area = _signed_area(cleaned)
    if area == 0.0:
        raise GirderError("the section encloses no area (all vertices are collinear)")

    if len(cleaned) <= MAX_VERTICES_FOR_SIMPLICITY_CHECK:
        n = len(cleaned)
        for i in range(n):
            for j in range(i + 2, n):
                if i == 0 and j == n - 1:
                    continue
                if _segments_cross(cleaned[i], cleaned[(i + 1) % n], cleaned[j], cleaned[(j + 1) % n]):
                    raise GirderError(
                        f"the section polygon crosses itself between edges {i} and {j}; "
                        "pass the loops in offset order"
                    )

    if area < 0.0:
        cleaned.reverse()
    return cleaned


def _signed_area(points: Sequence[Point]) -> float:
    total = 0.0
    n = len(points)
    for index in range(n):
        y0, z0 = points[index]
        y1, z1 = points[(index + 1) % n]
        total += y0 * z1 - y1 * z0
    return 0.5 * total


def _moments(points: Sequence[Point]) -> tuple[float, float, float, float, float, float]:
    """Area, centroid and raw second moments of a counter-clockwise polygon.

    Returns ``(area, c_y, c_z, int_y2, int_z2, int_yz)`` where the integrals are about
    the section origin, ``int_y2 = INT y^2 dA`` and so on. The three second moments use
    the Euler sums, which are exact for a polygon: the boundary integral of a quadratic
    over a straight edge is the trapezoid rule, so no quadrature error is introduced.
    """
    area = 0.0
    sum_y = 0.0
    sum_z = 0.0
    int_y2 = 0.0
    int_z2 = 0.0
    int_yz = 0.0
    n = len(points)
    for index in range(n):
        y0, z0 = points[index]
        y1, z1 = points[(index + 1) % n]
        cross = y0 * z1 - y1 * z0
        area += cross
        sum_y += (y0 + y1) * cross
        sum_z += (z0 + z1) * cross
        int_y2 += (y0 * y0 + y0 * y1 + y1 * y1) * cross
        int_z2 += (z0 * z0 + z0 * z1 + z1 * z1) * cross
        int_yz += (y0 * z1 + 2.0 * y0 * z0 + 2.0 * y1 * z1 + y1 * z0) * cross
    area *= 0.5
    if area == 0.0:
        raise GirderError("the section encloses no area")
    c_y = sum_y / (6.0 * area)
    c_z = sum_z / (6.0 * area)
    return area, c_y, c_z, int_y2 / 12.0, int_z2 / 12.0, int_yz / 24.0


def _cumulative(points: Sequence[Point], level: float) -> tuple[float, float]:
    """Area and first moment of the part of the polygon with ``z <= level``.

    Both are boundary integrals, so no polygon has to be reconstructed: with
    ``sigma(z) = (level - z)^+`` the area is ``INT sigma dA = CLOSED_INT sigma dy``
    (Green's theorem, the ``dz`` term vanishing), and the first moment is
    ``INT z sigma dA = CLOSED_INT P dy`` with the continuous potential
    ``P(z) = -min(z, level)^2 / 2``, whose derivative is ``-z sigma``. Along a straight
    edge the area integrand is linear in the parameter and the moment integrand is
    quadratic, so three-point Simpson is exact for both -- provided the edge is split
    where it crosses the level, which is what the loop below does. Because it is a
    boundary integral the method is blind to holes: an annulus entered as one loop with
    a zero-width bridge subtracts its void without any nesting analysis, and vertices
    sitting exactly on the level need no special case.
    """
    area_total, _cy, cz, _iy, _iz, _iyz = _moments(points)
    moment_total = area_total * cz

    area_below = 0.0
    moment_below = 0.0
    n = len(points)
    for index in range(n):
        y0, z0 = points[index]
        y1, z1 = points[(index + 1) % n]
        pieces = [(y0, z0, y1, z1)]
        if (z0 > level) != (z1 > level):
            t = (level - z0) / (z1 - z0)
            y_cut = y0 + t * (y1 - y0)
            # Both halves matter: the area integrand is zero above the level, but the
            # moment potential is the constant -level^2 / 2 there, not zero.
            pieces = [(y0, z0, y_cut, level), (y_cut, level, y1, z1)]
        for ya, za, yb, zb in pieces:
            dy = yb - ya
            # Three-point Simpson on the parameter t in [0, 1]: exact for the linear
            # area integrand and the quadratic moment integrand.
            area_below += dy * (
                _sigma(za, level) + 4.0 * _sigma(0.5 * (za + zb), level) + _sigma(zb, level)
            ) / 6.0
            moment_below += dy * (
                _potential(za, level)
                + 4.0 * _potential(0.5 * (za + zb), level)
                + _potential(zb, level)
            ) / 6.0

    area_above = area_total - area_below
    moment_above = moment_total - moment_below
    tolerance = _SPLIT_TOLERANCE_RELATIVE * abs(area_total)
    rebuilder = abs(area_below + area_above - area_total)
    if area_below < -tolerance or area_above < -tolerance or rebuilder > tolerance:
        raise GirderError(
            "the half-plane cut did not reconstruct the section: "
            f"below {area_below:.12g} + above {area_above:.12g} != {area_total:.12g} "
            f"at level {level:.12g}"
        )
    return area_below, moment_below, area_above, moment_above


def _sigma(z: float, level: float) -> float:
    """The half-plane indicator's antiderivative in ``z``: ``(level - z)^+``."""
    return level - z if z < level else 0.0


def _potential(z: float, level: float) -> float:
    """``-min(z, level)^2 / 2``, continuous, with ``dP/dz = -z (level - z)^+``."""
    return -0.5 * (z if z < level else level) ** 2


def _plastic_axis(points: Sequence[Point], axis: str) -> tuple[float, float]:
    """Equal-area neutral axis and plastic modulus for bending about ``axis``."""
    if axis == "vertical":
        # Transposing is a reflection, so it flips the winding order; restore it before
        # the half-plane machinery (which is written for a counter-clockwise loop).
        transposed = [(z, y) for y, z in points]
        if _signed_area(transposed) < 0.0:
            transposed.reverse()
        level, modulus = _plastic_axis(transposed, "horizontal")
        return level, modulus

    zs = [z for _, z in points]
    z_min, z_max = min(zs), max(zs)
    area, _cy, _cz, _iy, _iz, _iyz = _moments(points)
    half = 0.5 * area
    low, high = z_min, z_max
    tolerance = _PLASTIC_TOLERANCE_RELATIVE * (z_max - z_min)
    level = 0.5 * (low + high)
    for _ in range(_PLASTIC_MAX_ITERATIONS):
        level = 0.5 * (low + high)
        residual = _cumulative(points, level)[0] - half
        if residual == 0.0:
            break
        if residual < 0.0:
            low = level
        else:
            high = level
        if high - low <= tolerance:
            break
    area_below, moment_below, area_above, moment_above = _cumulative(points, level)
    # Z_p = INT |z - z_p| dA = (moment above) - (moment below) - z_p (area above - area below)
    modulus = (moment_above - moment_below) - level * (area_above - area_below)
    return level, modulus


def _principal(
    i_horizontal: float, i_vertical: float, product: float
) -> tuple[float, float, float]:
    """Principal second moments and the angle (radians) of the major one."""
    mean = 0.5 * (i_horizontal + i_vertical)
    radius = math.hypot(0.5 * (i_horizontal - i_vertical), product)
    major = mean + radius
    minor = mean - radius
    # Rotating the section by this angle aligns the major principal axis with the
    # horizontal one: 2 theta = atan2(2 I_yz, I_hh - I_vv) is the root of
    # (I_vv - I_hh) sin(2 theta) / 2 + I_yz cos(2 theta) = 0, i.e. I_yz' = 0.
    angle = 0.5 * math.atan2(2.0 * product, i_horizontal - i_vertical)
    return major, minor, angle


# ---------------------------------------------------------------------------
# the section object
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class GirderSection:
    """Elastic and plastic properties of one transverse section of a hull girder.

    ``second_moment_horizontal`` is ``INT (z - z_c)^2 dA`` -- the one that appears in
    ``sigma = M z' / I`` for a vertical bending moment -- and
    ``second_moment_vertical`` is ``INT (y - y_c)^2 dA``, for a horizontal moment.
    ``product_of_inertia`` is ``INT (y - y_c)(z - z_c) dA`` and is zero exactly when the
    section is symmetric about one of the two centroidal axes.

    The four section moduli are the elastic ones: the extreme-fibre distance is measured
    from the centroid, so a section whose centroid is not at mid-depth has different
    deck and keel moduli, which is the normal case for a hull with a superstructure-free
    deck and a full bilge.
    """

    area: float
    perimeter: float
    centroid_y: float
    centroid_z: float
    second_moment_horizontal: float
    second_moment_vertical: float
    product_of_inertia: float
    principal_second_moments: tuple[float, float]
    principal_angle_deg: float
    section_modulus_deck: float
    section_modulus_keel: float
    section_modulus_starboard: float
    section_modulus_port: float
    plastic_neutral_axis_z: float
    plastic_modulus: float
    plastic_neutral_axis_y: float
    plastic_modulus_vertical: float
    radius_of_gyration_horizontal: float
    radius_of_gyration_vertical: float
    y_min: float
    y_max: float
    z_min: float
    z_max: float
    vertices: int
    points: tuple[Point, ...]
    basis: str = "polygon"

    @property
    def depth(self) -> float:
        """Section depth, keel to deck."""
        return self.z_max - self.z_min

    @property
    def beam(self) -> float:
        """Section breadth, port to starboard."""
        return self.y_max - self.y_min

    @property
    def plastic_to_elastic_ratio(self) -> float:
        """``min(Z_p / Z_elastic)``: how much shape factor the section has.

        This is the shape factor of the governing fibre. It is ``1.5`` for a rectangle
        and rises with the material concentrated away from the neutral axis, which is
        why a hull girder with a thick deck and a wide bottom carries more than
        ``sigma = M / Z`` suggests before it becomes a mechanism.
        """
        elastic = min(
            self.section_modulus_deck,
            self.section_modulus_keel,
            self.section_modulus_starboard,
            self.section_modulus_port,
        )
        return self.plastic_modulus / elastic if elastic > 0.0 else math.inf

    def first_moment(self, z: float) -> float:
        """Static moment ``Q(z) = INT_{z..}^{..} (zeta - z) dA`` above the level ``z``.

        ``Q`` is the numerator of the shear formula ``tau = V Q / (I t)``: it vanishes at
        the extreme fibres (above the deck there is no area) and is largest at the
        neutral axis, which is where a hull girder's shear stress peaks.
        """
        _area_below, _moment_below, area_above, moment_above = _cumulative(self.points, z)
        return moment_above - z * area_above

    def to_dict(self) -> dict[str, object]:
        return {
            "area": self.area,
            "perimeter": self.perimeter,
            "beam": self.beam,
            "depth": self.depth,
            "centroid_y": self.centroid_y,
            "centroid_z": self.centroid_z,
            "second_moment_horizontal": self.second_moment_horizontal,
            "second_moment_vertical": self.second_moment_vertical,
            "product_of_inertia": self.product_of_inertia,
            "principal_second_moments": list(self.principal_second_moments),
            "principal_angle_deg": self.principal_angle_deg,
            "section_modulus_deck": self.section_modulus_deck,
            "section_modulus_keel": self.section_modulus_keel,
            "section_modulus_starboard": self.section_modulus_starboard,
            "section_modulus_port": self.section_modulus_port,
            "plastic_neutral_axis_z": self.plastic_neutral_axis_z,
            "plastic_modulus": self.plastic_modulus,
            "plastic_neutral_axis_y": self.plastic_neutral_axis_y,
            "plastic_modulus_vertical": self.plastic_modulus_vertical,
            "radius_of_gyration_horizontal": self.radius_of_gyration_horizontal,
            "radius_of_gyration_vertical": self.radius_of_gyration_vertical,
            "plastic_to_elastic_ratio": self.plastic_to_elastic_ratio,
            "y_min": self.y_min,
            "y_max": self.y_max,
            "z_min": self.z_min,
            "z_max": self.z_max,
            "vertices": self.vertices,
            "basis": self.basis,
        }


def girder_section(points: Iterable[Sequence[float]], *, basis: str = "polygon") -> GirderSection:
    """Section properties of a simple closed polygon in the ``(y, z)`` plane.

    The polygon is validated (finite vertices, at least three of them, non-zero area, no
    self-intersection) and oriented internally, so an offset table can be handed over in
    either winding order. The result is exact for the polygon: ``docs/girder_section.md``
    measures the residual against closed forms at ``10^-15`` relative, which is round-off
    and not a discretisation error.
    """
    cleaned = _clean_polygon(points)
    area, c_y, c_z, int_y2, int_z2, int_yz = _moments(cleaned)

    i_y2 = int_y2 - area * c_y * c_y  # about the vertical centroidal axis
    i_z2 = int_z2 - area * c_z * c_z  # about the horizontal centroidal axis
    product = int_yz - area * c_y * c_z

    major, minor, angle = _principal(i_z2, i_y2, product)

    zs = [z for _, z in cleaned]
    ys = [y for y, _ in cleaned]
    z_min, z_max = min(zs), max(zs)
    y_min, y_max = min(ys), max(ys)

    deck = z_max - c_z
    keel = c_z - z_min
    starboard = y_max - c_y
    port = c_y - y_min

    plastic_z, plastic_modulus = _plastic_axis(cleaned, "horizontal")
    plastic_y, plastic_modulus_v = _plastic_axis(cleaned, "vertical")

    perimeter = 0.0
    n = len(cleaned)
    for index in range(n):
        y0, z0 = cleaned[index]
        y1, z1 = cleaned[(index + 1) % n]
        perimeter += math.hypot(y1 - y0, z1 - z0)

    return GirderSection(
        area=area,
        perimeter=perimeter,
        centroid_y=c_y,
        centroid_z=c_z,
        second_moment_horizontal=i_z2,
        second_moment_vertical=i_y2,
        product_of_inertia=product,
        principal_second_moments=(major, minor),
        principal_angle_deg=math.degrees(angle),
        section_modulus_deck=i_z2 / deck if deck > 0.0 else math.inf,
        section_modulus_keel=i_z2 / keel if keel > 0.0 else math.inf,
        section_modulus_starboard=i_y2 / starboard if starboard > 0.0 else math.inf,
        section_modulus_port=i_y2 / port if port > 0.0 else math.inf,
        plastic_neutral_axis_z=plastic_z,
        plastic_modulus=plastic_modulus,
        plastic_neutral_axis_y=plastic_y,
        plastic_modulus_vertical=plastic_modulus_v,
        radius_of_gyration_horizontal=math.sqrt(i_z2 / area),
        radius_of_gyration_vertical=math.sqrt(i_y2 / area),
        y_min=y_min,
        y_max=y_max,
        z_min=z_min,
        z_max=z_max,
        vertices=len(cleaned),
        points=tuple(cleaned),
        basis=basis,
    )


def section_modulus_extremes(
    points: Iterable[Sequence[float]], axis: str = "horizontal"
) -> tuple[float, float]:
    """``(I / distance_to_far_fibre, I / distance_to_near_fibre)`` for the given axis.

    This is the pair a scantling loop iterates over: it wants the governing modulus
    (the smaller one) and the one that would govern if the neutral axis moved.
    """
    section = girder_section(points)
    if axis == "horizontal":
        return section.section_modulus_deck, section.section_modulus_keel
    if axis == "vertical":
        return section.section_modulus_starboard, section.section_modulus_port
    raise GirderError(f"axis must be 'horizontal' or 'vertical', got {axis!r}")


def mirror_half_section(
    z_values: Sequence[float], half_breadths: Sequence[float]
) -> tuple[Point, ...]:
    """Close a half-breadth table into a full section, mirroring about the centreline.

    ``z_values`` must be ascending and ``half_breadths`` holds ``b(z) >= 0`` in metres.
    The starboard side is traced upward, the deck line across, the port side downward and
    the bottom closed, so the returned loop is the full section and can be handed
    straight to :func:`girder_section`. The table is what an offset file or a loft
    produces; the mirror is exact, so a section that is symmetric by construction keeps
    ``product_of_inertia == 0`` to round-off.
    """
    if len(z_values) != len(half_breadths):
        raise GirderError(
            f"z_values and half_breadths must have the same length, got "
            f"{len(z_values)} and {len(half_breadths)}"
        )
    if len(z_values) < 2:
        raise GirderError("a half-breadth table needs at least two heights")

    zs = [float(z) for z in z_values]
    bs = [float(b) for b in half_breadths]
    for index, (z, b) in enumerate(zip(zs, bs, strict=True)):
        if not (math.isfinite(z) and math.isfinite(b)):
            raise GirderError(f"entry {index} is not finite: ({z!r}, {b!r})")
        if b < 0.0:
            raise GirderError(f"half breadth at entry {index} is negative: {b!r}")
        if index and z <= zs[index - 1]:
            raise GirderError(
                f"z_values must be strictly ascending: entry {index - 1} is {zs[index - 1]!r} "
                f"and entry {index} is {z!r}"
            )

    points: list[Point] = [(b, z) for z, b in zip(zs, bs, strict=True)]
    points.extend((-b, z) for z, b in zip(reversed(zs), reversed(bs), strict=True))
    return tuple(points)


def closed_cell_torsion_constant(
    enclosed_area: float, wall_thickness: float, perimeter: float
) -> float:
    """Bredt torsion constant of a single-cell thin-walled tube, ``J = 4 A_m^2 t / P``.

    ``enclosed_area`` is the area enclosed by the wall *mid-line* and ``perimeter`` its
    mid-line length; for a rectangular box of outer size ``b x h`` and thickness ``t``
    these are ``(b - t)(h - t)`` and ``2((b - t) + (h - t))``. Valid for ``t`` small
    against the cell size; wall thickness variation and multi-cell decks are not modelled.
    """
    if not all(math.isfinite(value) for value in (enclosed_area, wall_thickness, perimeter)):
        raise GirderError("torsion inputs must be finite")
    if enclosed_area <= 0.0 or wall_thickness <= 0.0 or perimeter <= 0.0:
        raise GirderError(
            "torsion inputs must be positive, got "
            f"area={enclosed_area!r}, t={wall_thickness!r}, P={perimeter!r}"
        )
    return 4.0 * enclosed_area * enclosed_area * wall_thickness / perimeter


def open_section_torsion_constant(strips: Iterable[Sequence[float]]) -> float:
    """``J = SUM b t^3 / 3`` for a set of non-interacting thin strips.

    Each entry is ``(length, thickness)``. This is the open-section rule; it is the
    reason an open profile is orders of magnitude softer in torsion than a closed one of
    the same weight, and it is the model behind "a ship needs a closed deck" rules.
    """
    total = 0.0
    counted = 0
    for index, strip in enumerate(strips):
        try:
            length, thickness = float(strip[0]), float(strip[1])
        except (TypeError, IndexError) as error:  # pragma: no cover - defensive
            raise GirderError(f"strip {index} is not a (length, thickness) pair: {strip!r}") from error
        if not (math.isfinite(length) and math.isfinite(thickness)):
            raise GirderError(f"strip {index} is not finite: {strip!r}")
        if length <= 0.0 or thickness <= 0.0:
            raise GirderError(f"strip {index} must have positive length and thickness: {strip!r}")
        total += length * thickness**3 / 3.0
        counted += 1
    if counted == 0:
        raise GirderError("at least one strip is required")
    return total


def box_girder_section(
    beam: float, depth: float, *, thickness: float = 0.0, inner_beam: float = 0.0, inner_depth: float = 0.0
) -> tuple[Point, ...]:
    """Rectangular or double-walled box section, as a polygon, for anchors and tests.

    ``thickness=0`` gives the solid rectangle ``beam x depth``; a positive ``thickness``
    gives the ring between the outer box and the concentric void
    ``(beam - 2t) x (depth - 2t)``, or the explicitly sized void ``inner_beam x
    inner_depth`` when those are given.
    """
    if beam <= 0.0 or depth <= 0.0:
        raise GirderError("beam and depth must be positive")
    half_b = 0.5 * beam
    if inner_beam or inner_depth:
        b_i, d_i = 0.5 * inner_beam, 0.5 * inner_depth
        if b_i >= half_b or d_i >= 0.5 * depth:
            raise GirderError("the void must be smaller than the outer box")
    elif thickness > 0.0:
        b_i = half_b - thickness
        d_i = 0.5 * depth - thickness
        if b_i <= 0.0 or d_i <= 0.0:
            raise GirderError("thickness leaves no void in the box")
    else:
        b_i = d_i = 0.0
    outer = [(-half_b, 0.0), (half_b, 0.0), (half_b, depth), (-half_b, depth)]
    if b_i == 0.0 and d_i == 0.0:
        return tuple(outer)
    # One loop that traces the outer boundary and the void, joined by a zero-width
    # bridge: a polygon representation of an annulus, which the shoelace sums handle
    # exactly because the two crossings of the bridge cancel.
    half_d = 0.5 * depth
    return (
        (-half_b, 0.0),
        (half_b, 0.0),
        (half_b, depth),
        (-half_b, depth),
        (-half_b, half_d),
        (-b_i, half_d),
        (-b_i, half_d + d_i),
        (b_i, half_d + d_i),
        (b_i, half_d - d_i),
        (-b_i, half_d - d_i),
        (-b_i, half_d),
        (-half_b, half_d),
    )
