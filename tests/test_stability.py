"""Initial stability: the waterplane moments against closed forms.

Three references, in increasing order of how much they can go wrong:

* a **box**, where every quantity has an exact rational formula and the mesh integral
  should reproduce it to round-off (``A = L B``, ``I_T = L B^3 / 12``,
  ``BM = B^2 / (12 T)``, ``KB = T / 2``);
* a **hemisphere**, where the section is a circle (``A = pi R^2``,
  ``I_T = I_L = pi R^4 / 4``) and, pleasingly, ``KB = BM = 3 R / 8`` -- two different
  integrals that must land on the same number, which no shared bug can fake;
* the **parametric hull** of :mod:`tensorcad.hull`, whose waterplane is the sine
  planform ``|y| <= (B/2) sin(pi x / L)^0.7``: ``I_T`` is closed-form through the Beta
  integral in ``hydrostatics_analytic.taper_integral`` and ``I_L`` is a quadrature that
  carries its own self-convergence evidence.
"""

from __future__ import annotations

import math
import subprocess
import sys

import pytest

from tensorcad import mesh as meshlib
from tensorcad.geometry import Mesh, MeshBuilder
from tensorcad.hull import HullParameters, generate_hull
from tensorcad.hydrostatics import compute_hydrostatics
from tensorcad.hydrostatics_analytic import (
    analytic_hydrostatics,
    waterplane_inertia_longitudinal,
    waterplane_inertia_transverse,
)
from tensorcad.stability import (
    StabilityError,
    WaterplaneMoments,
    initial_stability,
    waterplane_moments,
)


# ---------------------------------------------------------------------------
# reference geometries
# ---------------------------------------------------------------------------
def box(length: float, beam: float, draft: float, *, x0: float = 0.0) -> Mesh:
    """Axis-aligned closed box, top face on ``z = 0``, keel on ``z = -draft``."""
    builder = MeshBuilder()
    corner = {}
    for x in (x0, x0 + length):
        for y in (-0.5 * beam, 0.5 * beam):
            for z in (-draft, 0.0):
                corner[(x, y, z)] = builder.add_vertex((x, y, z))

    def quad(p0, p1, p2, p3) -> None:
        a, b, c, d = (corner[point] for point in (p0, p1, p2, p3))
        builder.add_triangle(a, b, c)
        builder.add_triangle(a, c, d)

    left, right = x0, x0 + length
    lo, hi = -0.5 * beam, 0.5 * beam
    top, bottom = 0.0, -draft
    # outward-facing winding; the exact orientation does not matter to the section
    # integral (each loop is normalised), only to the volume's sign.
    quad((left, lo, bottom), (left, hi, bottom), (right, hi, bottom), (right, lo, bottom))
    quad((left, lo, top), (right, lo, top), (right, hi, top), (left, hi, top))
    quad((left, lo, bottom), (right, lo, bottom), (right, lo, top), (left, lo, top))
    quad((left, hi, bottom), (left, hi, top), (right, hi, top), (right, hi, bottom))
    quad((left, lo, bottom), (left, lo, top), (left, hi, top), (left, hi, bottom))
    quad((right, lo, bottom), (right, hi, bottom), (right, hi, top), (right, lo, top))
    return builder.build()


def hemisphere(radius: float, n_lat: int, n_lon: int) -> Mesh:
    """Closed lower half of a sphere: a UV surface plus a flat disc at ``z = 0``."""
    builder = MeshBuilder()
    rings: list[list[int]] = []
    for i in range(n_lat):
        theta = 0.5 * math.pi * i / n_lat  # 0 at the equator, below pi/2 at the keel
        ring = []
        for j in range(n_lon):
            phi = 2.0 * math.pi * j / n_lon
            ring.append(
                builder.add_vertex(
                    (
                        radius * math.cos(theta) * math.cos(phi),
                        radius * math.cos(theta) * math.sin(phi),
                        -radius * math.sin(theta),
                    )
                )
            )
        rings.append(ring)
    pole = builder.add_vertex((0.0, 0.0, -radius))
    for i in range(n_lat - 1):
        upper, lower = rings[i], rings[i + 1]
        for j in range(n_lon):
            j2 = (j + 1) % n_lon
            builder.add_triangle(upper[j], lower[j], lower[j2])
            builder.add_triangle(upper[j], lower[j2], upper[j2])
    last = rings[-1]
    for j in range(n_lon):
        j2 = (j + 1) % n_lon
        builder.add_triangle(last[j], pole, last[j2])
    # flat disc closing the hemisphere at the waterline
    equator = rings[0]
    centre = builder.add_vertex((0.0, 0.0, 0.0))
    for j in range(n_lon):
        j2 = (j + 1) % n_lon
        builder.add_triangle(centre, equator[j], equator[j2])
    mesh = builder.build()
    # The helper is a reference body, not a test of winding: orient it outward so the
    # signed volume (and therefore KB) is positive whatever the construction order.
    if meshlib.total_signed_volume(mesh) < 0.0:
        mesh = Mesh(
            vertices=mesh.vertices,
            faces=tuple((face[0], face[2], face[1]) for face in mesh.faces),
            units=mesh.units,
            coordinate_frame=mesh.coordinate_frame,
        )
    return mesh


def relative(value: float, reference: float) -> float:
    return value / reference - 1.0


# ---------------------------------------------------------------------------
# the box: every quantity is exact
# ---------------------------------------------------------------------------
def test_box_section_moments_are_exact() -> None:
    length, beam, draft = 10.0, 2.0, 0.8
    section = waterplane_moments(box(length, beam, draft), 0.0)
    assert isinstance(section, WaterplaneMoments)
    assert section.loops == 1
    assert relative(section.area, length * beam) == pytest.approx(0.0, abs=1e-14)
    assert relative(section.transverse_second_moment, length * beam**3 / 12.0) == (
        pytest.approx(0.0, abs=1e-14)
    )
    assert relative(section.longitudinal_second_moment, beam * length**3 / 12.0) == (
        pytest.approx(0.0, abs=1e-12)
    )
    assert section.longitudinal_centre_flotation == pytest.approx(0.5 * length, abs=1e-12)
    assert section.transverse_centre_flotation == pytest.approx(0.0, abs=1e-12)


def test_box_metacentric_height_reproduces_the_textbook_formula() -> None:
    length, beam, draft = 10.0, 2.0, 0.8
    kg = 0.55
    result = initial_stability(box(length, beam, draft), kg=kg)
    assert result.displacement_volume == pytest.approx(length * beam * draft, rel=1e-13)
    assert result.kb == pytest.approx(0.5 * draft, abs=1e-12)
    assert result.bm_transverse == pytest.approx(beam**2 / (12.0 * draft), rel=1e-13)
    # BM_T here is 0.4167 vs KB 0.4: a box's transverse metacentre sits above B.
    assert result.kb + result.bm_transverse > result.kb
    assert result.gm_transverse == pytest.approx(0.5 * draft + beam**2 / (12.0 * draft) - kg, abs=1e-12)
    # and the longitudinal metacentre is far higher: I_L / I_T = (L / B)^2
    assert result.gm_longitudinal == pytest.approx(0.5 * draft + length**2 / (12.0 * draft) - kg, abs=1e-10)
    assert result.bm_longitudinal / result.bm_transverse == pytest.approx((length / beam) ** 2, rel=1e-10)


def test_a_box_that_is_too_heavy_is_reported_unstable() -> None:
    length, beam, draft = 10.0, 2.0, 0.8
    # KM = 0.4 + 0.4167 = 0.8167, so KG above it is negative GM -- reported, not hidden.
    result = initial_stability(box(length, beam, draft), kg=1.2)
    assert result.gm_transverse < 0.0
    assert result.gm_transverse == pytest.approx(0.8166666667 - 1.2, abs=1e-6)


# ---------------------------------------------------------------------------
# the hemisphere: two independent integrals that must agree (KB = BM = 3R/8)
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("n_lat", [12, 24, 48])
def test_hemisphere_section_matches_the_circle(n_lat: int) -> None:
    radius = 1.0
    n_lon = 2 * n_lat
    section = waterplane_moments(hemisphere(radius, n_lat, n_lon), 0.0)
    # An inscribed n-gon under-counts a circle, and the deficit is O(1/n^2): measured
    # -1.14 %, -0.285 %, -0.0714 % of the area at n_lon = 24, 48, 96, i.e. one quarter
    # per doubling. The bounds below are that measured envelope, not a guessed one.
    assert relative(section.area, math.pi * radius**2) == pytest.approx(0.0, abs=0.02)
    assert relative(section.transverse_second_moment, math.pi * radius**4 / 4.0) == (
        pytest.approx(0.0, abs=0.03)
    )
    assert relative(section.longitudinal_second_moment, math.pi * radius**4 / 4.0) == (
        pytest.approx(0.0, abs=0.03)
    )
    # a circle's second moment is the same about every diameter through the centre
    assert section.transverse_second_moment == pytest.approx(section.longitudinal_second_moment, rel=1e-2)


def test_hemisphere_metacentre_sits_at_the_sphere_centre() -> None:
    """KB = 5R/8, BM = I/V = 3R/8, so KM = R exactly: two integrals, one identity.

    ``KB`` comes from the volume integration in :mod:`tensorcad.hydrostatics` and ``BM``
    from the waterplane moment in this module; neither code path knows about the other,
    so them summing to the sphere's radius is a real cross-check rather than a tautology.
    """
    radius = 1.0
    result = initial_stability(hemisphere(radius, 48, 96), kg=0.5)
    assert relative(result.kb, 5.0 * radius / 8.0) == pytest.approx(0.0, abs=0.01)
    assert relative(result.bm_transverse, 3.0 * radius / 8.0) == pytest.approx(0.0, abs=0.02)
    assert relative(result.km_transverse, radius) == pytest.approx(0.0, abs=0.01)
    assert result.gm_transverse == pytest.approx(result.km_transverse - 0.5, abs=1e-12)


def test_the_circle_error_falls_with_resolution() -> None:
    radius = 1.0
    coarse = waterplane_moments(hemisphere(radius, 8, 16), 0.0)
    fine = waterplane_moments(hemisphere(radius, 32, 64), 0.0)
    coarse_error = abs(relative(coarse.area, math.pi * radius**2))
    fine_error = abs(relative(fine.area, math.pi * radius**2))
    assert fine_error < coarse_error / 4.0


# ---------------------------------------------------------------------------
# the parametric hull against its analytic waterplane
# ---------------------------------------------------------------------------
def test_hull_section_agrees_with_the_analytic_planform() -> None:
    parameters = HullParameters(length=10.0, beam=2.0, draft=0.8, nx=41, nz=13)
    model = generate_hull(parameters)
    section = waterplane_moments(model.mesh, 0.0)
    analytic = analytic_hydrostatics(parameters)
    inertia_t = waterplane_inertia_transverse(parameters.length, parameters.beam)
    assert relative(section.area, analytic.waterplane_area) == pytest.approx(0.0, abs=0.03)
    assert relative(section.transverse_second_moment, inertia_t) == pytest.approx(0.0, abs=0.03)
    assert section.longitudinal_centre_flotation == pytest.approx(0.5 * parameters.length, abs=1e-9)
    assert section.transverse_centre_flotation == pytest.approx(0.0, abs=1e-12)


def test_hull_longitudinal_moment_agrees_with_quadrature() -> None:
    parameters = HullParameters(length=10.0, beam=2.0, draft=0.8, nx=41, nz=13)
    model = generate_hull(parameters)
    section = waterplane_moments(model.mesh, 0.0)
    inertia_l, evidence = waterplane_inertia_longitudinal(10.0, 2.0)
    assert evidence < 1e-6
    assert relative(section.longitudinal_second_moment, inertia_l) == pytest.approx(0.0, abs=0.03)


def test_hull_metacentric_height_against_the_analytic_surface() -> None:
    parameters = HullParameters(length=10.0, beam=2.0, draft=0.8, nx=161, nz=41)
    model = generate_hull(parameters)
    kg = 0.6
    result = initial_stability(model.mesh, kg=kg)
    analytic = analytic_hydrostatics(parameters)
    analytic_bm = (
        waterplane_inertia_transverse(parameters.length, parameters.beam) / analytic.displacement_volume
    )
    assert relative(result.displacement_volume, analytic.displacement_volume) == (
        pytest.approx(0.0, abs=0.005)
    )
    assert relative(result.bm_transverse, analytic_bm) == pytest.approx(0.0, abs=0.01)
    # KB from the mesh integration, against the closed form of the sqrt profile:
    # z-bar = -2T/5 measured from the waterline, i.e. 3T/5 above the keel.
    assert relative(result.kb, 0.6 * parameters.draft) == pytest.approx(0.0, abs=0.01)
    assert result.gm_transverse == pytest.approx(result.kb + result.bm_transverse - kg, abs=1e-12)


def test_hull_waterplane_area_agrees_with_the_flux_integral() -> None:
    """Two independent integrations of one number: Green on the loop, flux on the cap."""
    model = generate_hull(HullParameters(length=10.0, beam=2.0, draft=0.8, nx=81, nz=25))
    section = waterplane_moments(model.mesh, 0.0)
    hydrostatics = compute_hydrostatics(model.mesh, 0.0)
    assert relative(section.area, hydrostatics.waterplane_area) == pytest.approx(0.0, abs=1e-9)
    assert section.longitudinal_centre_flotation == pytest.approx(
        hydrostatics.longitudinal_centre_flotation, abs=1e-9
    )


def test_gm_falls_one_for_one_with_kg() -> None:
    model = generate_hull(HullParameters(length=10.0, beam=2.0, draft=0.8, nx=41, nz=13))
    light = initial_stability(model.mesh, kg=0.5)
    heavy = initial_stability(model.mesh, kg=0.7)
    assert light.gm_transverse - heavy.gm_transverse == pytest.approx(0.2, abs=1e-15)
    assert light.km_transverse == pytest.approx(heavy.km_transverse, rel=0.0, abs=0.0)


# ---------------------------------------------------------------------------
# more than one loop, and the failure modes
# ---------------------------------------------------------------------------
def test_two_disjoint_boxes_add_up() -> None:
    length, beam, draft, gap = 4.0, 1.0, 0.5, 6.0
    builder = MeshBuilder()
    meshes = [box(length, beam, draft, x0=0.0), box(length, beam, draft, x0=gap)]
    # merge the two boxes into one mesh by re-adding their triangles
    index = {}
    for mesh in meshes:
        for face in mesh.faces:
            triangle = []
            for vertex in face:
                point = mesh.vertices[vertex]
                triangle.append(index.setdefault(point, builder.add_vertex(point)))
            builder.add_triangle(*triangle)
    section = waterplane_moments(builder.build(), 0.0)
    assert section.loops == 2
    assert relative(section.area, 2.0 * length * beam) == pytest.approx(0.0, abs=1e-12)
    assert relative(section.transverse_second_moment, 2.0 * length * beam**3 / 12.0) == (
        pytest.approx(0.0, abs=1e-12)
    )
    offset = 0.5 * (length + gap) - 0.5 * length  # centre-to-centre distance
    expected_longitudinal = 2.0 * (beam * length**3 / 12.0 + length * beam * offset**2)
    assert relative(section.longitudinal_second_moment, expected_longitudinal) == (
        pytest.approx(0.0, abs=1e-12)
    )
    assert section.longitudinal_centre_flotation == pytest.approx(0.5 * (length + gap), abs=1e-12)


def test_a_waterline_that_misses_the_body_is_refused() -> None:
    with pytest.raises(StabilityError):
        waterplane_moments(box(10.0, 2.0, 0.8), 5.0)
    with pytest.raises(StabilityError):
        waterplane_moments(box(10.0, 2.0, 0.8), -5.0)


def test_a_non_finite_waterline_or_kg_is_refused() -> None:
    mesh = box(10.0, 2.0, 0.8)
    with pytest.raises(StabilityError):
        waterplane_moments(mesh, float("nan"))
    with pytest.raises(StabilityError):
        initial_stability(mesh, kg=float("inf"))


def test_a_body_with_no_displacement_is_refused() -> None:
    """A degenerate (zero-volume) surface must not produce an infinite GM."""
    builder = MeshBuilder()
    a = builder.add_vertex((0.0, 0.0, 0.0))
    b = builder.add_vertex((1.0, 0.0, 0.0))
    c = builder.add_vertex((0.0, 1.0, 0.0))
    builder.add_triangle(a, b, c)
    with pytest.raises(StabilityError):
        initial_stability(builder.build(), kg=1.0)


def test_the_report_serialises_to_rounded_json() -> None:
    payload = initial_stability(box(10.0, 2.0, 0.8), kg=0.55).to_dict()
    assert payload["gm_transverse_m"] == pytest.approx(0.2666666667, abs=1e-9)
    assert payload["waterplane"]["transverse_second_moment_m4"] == pytest.approx(10.0 * 8.0 / 12.0, abs=1e-9)
    assert payload["basis"].startswith("solid GM")


def test_the_module_imports_without_torch_or_numpy() -> None:
    code = """
import sys
class Blocker:
    def find_module(self, name, path=None):
        if name in ("torch", "numpy", "safetensors"):
            return self
    def load_module(self, name):
        raise ImportError("BLOCKED: " + name)
sys.meta_path.insert(0, Blocker())
import tensorcad.stability as stability
print(stability.SCHEMA_STABILITY)
"""
    completed = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=False)
    assert completed.returncode == 0, completed.stderr
    assert "tensorcad.stability.initial/1" in completed.stdout


# ===========================================================================
# m2925 suite: the waterplane-inertia and hydrostatics entry point, kept as it
# was written so the union does not quietly weaken either side's evidence.
#
# Three independent references are used, because each one closes a different hole:
#
# * A **rectangular barge** has closed-form waterplane inertias (``L*B**3/12`` and
#   ``B*L**3/12``), ``KB = T/2`` and ``BM = I/nabla``. The barge mesh is built here
#   vertex by vertex, so the answer is checked against arithmetic rather than against
#   another mesh, and it is checked on three tessellations so a lucky cancellation
#   cannot pass.
# * A **circular waterplane** (a vertical cylinder) has ``I_T = I_L = pi*R**4/4``; the
#   inscribed polygon used to mesh it converges to that from below, which tests the
#   transverse and longitudinal paths with a shape whose surface has no flat bottom.
# * The **analytic hull** of :mod:`tensorcad.hull` has the closed forms in
#   :mod:`tensorcad.hydrostatics_analytic` -- one of which is a gamma identity, the
#   other a tanh-sinh quadrature -- so the mesh integrator can be measured at the
#   resolutions users actually run.
#
# ``HullParameters``, ``analytic_hydrostatics``, ``generate_hull`` and
# ``initial_stability`` are the ones imported at the top of this file.
# ===========================================================================
from tensorcad.hydrostatics_analytic import (  # noqa: E402
    longitudinal_waterplane_inertia,
    tanh_sinh_quadrature,
    taper_integral,
    transverse_waterplane_inertia,
)


def box_mesh(
    length: float, beam: float, depth: float, nx: int = 1, ny: int = 1, nz: int = 1
) -> Mesh:
    """Outward-oriented barge: ``x`` aft->fore, ``y`` across, ``z`` keel (-depth) -> 0.

    The waterline is the top face at ``z = 0``, so the submerged body is the whole
    box and every quantity has a closed form. Vertices are shared by construction
    through :class:`MeshBuilder`, so the surface is watertight without any welding.
    """
    builder = MeshBuilder(units="m", coordinate_frame="box test frame")

    def corner(i: int, j: int, k: int) -> int:
        """Vertex index of a grid corner (``MeshBuilder`` deduplicates by position)."""
        return builder.add_vertex(
            (length * i / nx, -0.5 * beam + beam * j / ny, -depth + depth * k / nz)
        )

    def quad(a: int, b: int, c: int, d: int) -> None:
        builder.add_quad(a, b, c, d)

    for i in range(nx):
        for j in range(ny):
            # bottom (outward -z) and top (outward +z)
            quad(corner(i, j, 0), corner(i, j + 1, 0), corner(i + 1, j + 1, 0), corner(i + 1, j, 0))
            quad(
                corner(i, j, nz),
                corner(i + 1, j, nz),
                corner(i + 1, j + 1, nz),
                corner(i, j + 1, nz),
            )
    for j in range(ny):
        for k in range(nz):
            quad(corner(0, j, k), corner(0, j, k + 1), corner(0, j + 1, k + 1), corner(0, j + 1, k))
            quad(
                corner(nx, j, k),
                corner(nx, j + 1, k),
                corner(nx, j + 1, k + 1),
                corner(nx, j, k + 1),
            )
    for i in range(nx):
        for k in range(nz):
            quad(corner(i, 0, k), corner(i + 1, 0, k), corner(i + 1, 0, k + 1), corner(i, 0, k + 1))
            quad(
                corner(i, ny, k),
                corner(i, ny, k + 1),
                corner(i + 1, ny, k + 1),
                corner(i + 1, ny, k),
            )
    return builder.build()


def cylinder_mesh(radius: float, depth: float, segments: int) -> Mesh:
    """Vertical cylinder, closed, with the waterline as its top face at ``z = 0``."""
    builder = MeshBuilder(units="m", coordinate_frame="cylinder test frame")
    angles = [2.0 * math.pi * i / segments for i in range(segments)]

    def rim(theta: float, z: float) -> int:
        return builder.add_vertex(
            (radius * math.cos(theta), radius * math.sin(theta), z)
        )

    bottom = builder.add_vertex((0.0, 0.0, -depth))
    top = builder.add_vertex((0.0, 0.0, 0.0))
    for index in range(segments):
        a = angles[index]
        b = angles[(index + 1) % segments]
        builder.add_quad(rim(a, -depth), rim(b, -depth), rim(b, 0.0), rim(a, 0.0))
        builder.add_triangle(bottom, rim(b, -depth), rim(a, -depth))
        builder.add_triangle(top, rim(a, 0.0), rim(b, 0.0))
    return builder.build()


# --------------------------------------------------------------------------- barge


@pytest.mark.parametrize("tessellation", [(1, 1, 1), (4, 3, 2), (20, 6, 4)])
def test_rectangular_barge_is_exact(tessellation: tuple[int, int, int]) -> None:
    """A box waterplane is a rectangle: every quantity has an exact closed form."""
    length, beam, depth = 20.0, 6.0, 2.0
    mesh = box_mesh(length, beam, depth, *tessellation)
    hydro = compute_hydrostatics(mesh, 0.0, water_density=1000.0)

    transverse = length * beam**3 / 12.0
    longitudinal = beam * length**3 / 12.0
    volume = length * beam * depth

    assert hydro.basis == "closed surface integrals"
    assert hydro.displacement_volume == pytest.approx(volume, rel=1e-12)
    assert hydro.waterplane_area == pytest.approx(length * beam, rel=1e-12)
    assert hydro.transverse_waterplane_inertia == pytest.approx(transverse, rel=1e-12)
    assert hydro.longitudinal_waterplane_inertia == pytest.approx(longitudinal, rel=1e-12)
    assert hydro.keel_z == pytest.approx(-depth, abs=1e-12)

    stability = initial_stability(hydro, height_of_gravity=depth / 2.0)
    assert stability.height_of_centre_of_buoyancy == pytest.approx(depth / 2.0, rel=1e-12)
    assert stability.transverse_metacentric_radius == pytest.approx(transverse / volume, rel=1e-12)
    assert stability.longitudinal_metacentric_radius == pytest.approx(
        longitudinal / volume, rel=1e-12
    )
    # KG at half depth = KB, so GM is exactly BM.
    assert stability.transverse_metacentric_height == pytest.approx(transverse / volume, rel=1e-12)


def test_barge_inertia_does_not_depend_on_tessellation() -> None:
    """Refining a flat waterplane must not move the integrals at all."""
    coarse = compute_hydrostatics(box_mesh(12.0, 4.0, 1.5, 1, 1, 1), 0.0)
    fine = compute_hydrostatics(box_mesh(12.0, 4.0, 1.5, 6, 4, 3), 0.0)
    assert coarse.transverse_waterplane_inertia == pytest.approx(
        fine.transverse_waterplane_inertia, rel=1e-12
    )
    assert coarse.longitudinal_waterplane_inertia == pytest.approx(
        fine.longitudinal_waterplane_inertia, rel=1e-12
    )


# ------------------------------------------------------------------------ cylinder


def test_circular_waterplane_converges_to_the_closed_form() -> None:
    """An inscribed polygon approaches ``pi*R**4/4`` from below, and stays below."""
    radius, depth = 1.5, 2.0
    exact = math.pi * radius**4 / 4.0
    errors = []
    for segments in (24, 48, 96, 192):
        hydro = compute_hydrostatics(cylinder_mesh(radius, depth, segments), 0.0)
        transverse = hydro.transverse_waterplane_inertia
        longitudinal = hydro.longitudinal_waterplane_inertia
        # A waterplane that is symmetric in both axes must give the same two
        # moments, and an inscribed polygon can only be smaller than the circle.
        assert transverse == pytest.approx(longitudinal, rel=1e-9)
        assert transverse < exact
        errors.append(exact - transverse)
    # Inscribed regular polygon: the area deficit is O(1/n**2), so each doubling of
    # the segment count should cut the error by about four.
    ratios = [errors[i] / errors[i + 1] for i in range(len(errors) - 1)]
    assert all(3.5 < ratio < 4.5 for ratio in ratios), ratios
    assert errors[-1] / exact < 5e-4


# -------------------------------------------------------------------- analytic hull


def test_analytic_closed_forms_match_their_own_quadrature() -> None:
    """Transverse: gamma identity vs quadrature of ``int (2/3)*b**3 dx``.

    The closed form and the quadrature are derived separately (a gamma function and
    a double-exponential rule), so agreement is a real check rather than a
    restatement of the same algebra.
    """
    length, beam, exponent = 10.0, 2.0, 0.7
    closed = transverse_waterplane_inertia(length, beam, exponent)

    def integrand(xi: float) -> float:
        return (2.0 / 3.0) * (0.5 * beam) ** 3 * math.sin(math.pi * xi) ** (3.0 * exponent)

    numerical = length * tanh_sinh_quadrature(integrand, 0.0, 1.0, levels=7)
    assert closed == pytest.approx(numerical, rel=1e-12)


def test_longitudinal_inertia_quadrature_reports_its_own_error() -> None:
    """The reference carries its convergence evidence instead of asking to be believed."""
    evidence = longitudinal_waterplane_inertia(10.0, 2.0)
    assert evidence["value"] > 0.0
    assert evidence["relative_difference"] < 1e-10
    # Scaling laws: I_L goes with beam and with length cubed.
    doubled = longitudinal_waterplane_inertia(20.0, 2.0)["value"] / 8.0
    assert doubled == pytest.approx(evidence["value"], rel=1e-9)


@pytest.mark.parametrize("parameters", [(10.0, 2.0, 0.8), (60.0, 14.0, 6.0), (4.3561, 1.0, 0.4)])
def test_hull_mesh_inertia_matches_the_analytic_surface(
    parameters: tuple[float, float, float]
) -> None:
    """The sampled hull against the closed forms of the surface it samples."""
    length, beam, draft = parameters
    reference = analytic_hydrostatics((length, beam, draft))

    default = generate_hull(HullParameters(length=length, beam=beam, draft=draft))
    measured = default.hydrostatics
    assert measured.transverse_waterplane_inertia == pytest.approx(
        reference.transverse_waterplane_inertia, rel=0.03
    )
    assert measured.longitudinal_waterplane_inertia == pytest.approx(
        reference.longitudinal_waterplane_inertia, rel=0.03
    )

    fine = generate_hull(HullParameters(length=length, beam=beam, draft=draft, nx=161, nz=41))
    assert fine.hydrostatics.transverse_waterplane_inertia == pytest.approx(
        reference.transverse_waterplane_inertia, rel=0.005
    )
    assert fine.hydrostatics.longitudinal_waterplane_inertia == pytest.approx(
        reference.longitudinal_waterplane_inertia, rel=0.005
    )
    # Refining has to help, not just land inside the band by luck.
    assert abs(
        fine.hydrostatics.transverse_waterplane_inertia
        - reference.transverse_waterplane_inertia
    ) < abs(measured.transverse_waterplane_inertia - reference.transverse_waterplane_inertia)


def test_metacentric_heights_of_the_analytic_hull() -> None:
    """KB/BM/KM/GM of the default hull, against the closed-form reference."""
    model = generate_hull()
    reference = analytic_hydrostatics()
    stability = initial_stability(model.hydrostatics, height_of_gravity=0.6)

    # KB = 0.6*T for the sqrt profile (VCB = -0.4*T with the keel at -T).
    assert stability.height_of_centre_of_buoyancy == pytest.approx(0.6 * 0.8, rel=0.01)
    assert stability.height_of_centre_of_buoyancy == pytest.approx(
        reference.vertical_centre_buoyancy + 0.8, rel=0.01
    )
    assert stability.transverse_metacentric_radius == pytest.approx(
        reference.transverse_metacentric_radius, rel=0.03
    )
    assert stability.height_of_transverse_metacentre == pytest.approx(
        reference.height_of_metacentre_above_keel, rel=0.02
    )
    assert stability.transverse_metacentric_height == pytest.approx(
        stability.height_of_transverse_metacentre - 0.6, rel=1e-12
    )
    # A slender hull has a far larger longitudinal metacentric radius; that
    # ordering is a property of the form, not of the discretisation.
    assert stability.longitudinal_metacentric_radius > 10.0 * stability.transverse_metacentric_radius
    assert stability.is_stable_in_roll


def test_raising_the_centre_of_gravity_can_make_the_hull_unstable() -> None:
    """GM is reported through zero rather than clamped: an unstable hull is an answer."""
    hydro = generate_hull().hydrostatics
    stable = initial_stability(hydro, height_of_gravity=0.5)
    unstable = initial_stability(hydro, height_of_gravity=1.2)
    assert stable.transverse_metacentric_height > 0.0
    assert unstable.transverse_metacentric_height < 0.0
    assert not unstable.is_stable_in_roll
    # KG never changes KB/BM: only the subtraction moves.
    assert unstable.height_of_transverse_metacentre == pytest.approx(
        stable.height_of_transverse_metacentre, rel=1e-15
    )


# ------------------------------------------------------------------------ reporting


def test_stability_record_is_json_ready() -> None:
    stability = initial_stability(generate_hull().hydrostatics, height_of_gravity=0.55)
    record = stability.to_dict()
    for key in (
        "kb_m",
        "bm_transverse_m",
        "bm_longitudinal_m",
        "km_transverse_m",
        "kg_m",
        "gm_transverse_m",
        "is_stable_in_roll",
        "waterplane_inertia_transverse_m4",
        "waterplane_inertia_longitudinal_m4",
    ):
        assert key in record
    assert record["is_stable_in_roll"] is True
    assert record["gm_transverse_m"] == pytest.approx(
        record["km_transverse_m"] - record["kg_m"], abs=1e-8
    )


def test_hydrostatics_payload_carries_the_inertias() -> None:
    record = generate_hull().hydrostatics.to_dict()
    assert record["waterplane_inertia_transverse_m4"] > 0.0
    assert record["waterplane_inertia_longitudinal_m4"] > record["waterplane_inertia_transverse_m4"]
    assert record["keel_z_m"] == pytest.approx(-0.8, abs=1e-9)


# ------------------------------------------------------------------------ rejection


def test_height_of_gravity_must_be_a_finite_length() -> None:
    hydro = generate_hull().hydrostatics
    with pytest.raises(ValueError, match="finite length"):
        initial_stability(hydro, height_of_gravity=float("nan"))
    with pytest.raises(ValueError, match="cannot be negative"):
        initial_stability(hydro, height_of_gravity=-0.1)


def test_a_body_without_a_waterplane_is_rejected() -> None:
    """A closed surface entirely below the waterline has no free surface to float on."""
    mesh = box_mesh(4.0, 2.0, 1.0, 2, 2, 2)
    submerged = compute_hydrostatics(mesh, waterline_z=2.0)
    assert math.isnan(submerged.transverse_waterplane_inertia)
    with pytest.raises(ValueError, match="waterplane"):
        initial_stability(submerged, height_of_gravity=0.5)


def test_open_surface_inertia_is_not_claimed() -> None:
    """An unclosed hull shell gives no trustworthy moments; the value is nan."""
    mesh = generate_hull().mesh
    open_faces = tuple(face for face in mesh.faces if max(mesh.vertices[v][2] for v in face) < -0.2)
    shell = Mesh(vertices=mesh.vertices, faces=open_faces)
    hydro = compute_hydrostatics(shell, 0.0)
    assert hydro.basis == "open surface (volume not trustworthy)"
    assert math.isnan(hydro.transverse_waterplane_inertia)
    assert math.isnan(hydro.longitudinal_waterplane_inertia)
    with pytest.raises(ValueError, match="waterplane"):
        initial_stability(hydro, height_of_gravity=0.5)


def test_stability_module_is_standard_library_only() -> None:
    """The zero-dependency promise has to hold for this module too."""

    class Blocker:
        blocked: list[str] = []

        def find_module(self, name, path=None):  # noqa: ANN001, ARG002
            if name in ("torch", "numpy", "scipy"):
                self.blocked.append(name)
                return self
            return None

        def load_module(self, name):  # noqa: ANN001
            raise ImportError(f"blocked: {name}")

    blocker = Blocker()
    saved = {name: sys.modules.pop(name) for name in ("tensorcad.stability",) if name in sys.modules}
    sys.meta_path.insert(0, blocker)
    try:
        import tensorcad.stability as module  # noqa: F401,PLC0415
    finally:
        sys.meta_path.remove(blocker)
        sys.modules.update(saved)
    assert blocker.blocked == []