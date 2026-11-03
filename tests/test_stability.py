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
