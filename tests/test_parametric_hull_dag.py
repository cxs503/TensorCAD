"""Acceptance tests for the parametric hull generation node (p316w1-n03).

Each test pins one machine-checkable acceptance criterion or one contract of the
generation pipeline: the interpolating loft, the exact hydrostatic integrals, the
watertight closed mesh, the density rule and the parametric families.  References are
closed-form (generalised Wigley family) or published (DARPA SUBOFF offset table); see
``examples/parametric_hull_dag/run_verification.py`` for the full campaign and
``docs/parametric_hull/VERIFICATION.md`` for the measured numbers.
"""
from __future__ import annotations

import math

import pytest

from tensorcad import mesh as meshlib
from tensorcad.parametric_hull import (
    CrossSection,
    ParametricHullError,
    Transform,
    apply_transform,
    build_parametric_hull,
    compute_hydrostatics,
    cross_section_from_function,
    depth_exponent_for_block_coefficient,
    interpolate_control_points,
    loft_nurbs_surface,
    mesh_density_for,
    mesh_quality_report,
    suboff_cross_sections,
    suboff_station_fractions,
    tessellate_surface,
    wigley_analytic_properties,
    wigley_cross_sections,
    wigley_family_sections,
)
from tensorcad.suboff import (
    SUBOFF_LENGTH_M,
    SUBOFF_RADIUS_M,
    suboff_radius_fraction,
    suboff_reference_geometry,
)

L, B, T = 10.0, 1.0, 0.625
CRITERION = {
    "geometry": 1e-3,
    "volume": 5e-3,
    "waterplane": 1e-2,
    "suboff": 5e-3,
    "min_scaled_jacobian": 0.30,
    "max_aspect_ratio": 5.0,
}


# --------------------------------------------------------------------------- #
# helpers under test: interpolation and parametrisation
# --------------------------------------------------------------------------- #
def test_interpolate_control_points_reproduce_data():
    """The B-spline built from the returned control points passes through the data."""
    data = [(0.0, 0.0), (1.0, 2.0), (2.0, 1.0), (3.0, 0.5), (4.0, 0.0)]
    points = [(float(i), 0.0, y) for i, (_, y) in enumerate(data)]
    params = (0.0, 0.25, 0.5, 0.75, 1.0)
    control, knots = interpolate_control_points(points, degree=3, params=params)
    assert len(control) == len(points)
    assert len(knots) == len(points) + 3 + 1
    # de Boor evaluation at the chord-length parameters must return the data
    from tensorcad.nurbs import NURBSCurve

    curve = NURBSCurve(degree=3, knots=knots, control_points=control,
                       weights=(1.0,) * len(control))
    for u, expected in zip(params, points):
        got = curve.evaluate(u)
        assert abs(got[2] - expected[2]) < 1e-9


def test_interpolate_control_points_rejects_degenerate_input():
    with pytest.raises(ParametricHullError):
        interpolate_control_points([(0.0, 0.0, 0.0)], degree=3)


# --------------------------------------------------------------------------- #
# (a) Wigley III geometry
# --------------------------------------------------------------------------- #
def test_wigley_surface_matches_analytic_ordinates():
    """Criterion (a): 10 stations x 5 waterlines within 0.1% of the closed form."""
    hull = build_parametric_hull(wigley_cross_sections(L, B, T, 11, 12))
    surface = hull.surface
    worst = 0.0
    for i in range(10):
        x = (0.05 + 0.1 * i) * L
        lo, hi = 0.0, 1.0
        for _ in range(70):
            mid = 0.5 * (lo + hi)
            if surface.evaluate(mid, 0.5)[0] < x:
                lo = mid
            else:
                hi = mid
        u = 0.5 * (lo + hi)
        for z in (-0.95 * T, -0.75 * T, -0.5 * T, -0.25 * T, 0.0):
            lo, hi = 0.0, 1.0
            rising = surface.evaluate(u, 1.0)[2] >= surface.evaluate(u, 0.0)[2]
            for _ in range(70):
                mid = 0.5 * (lo + hi)
                z_mid = surface.evaluate(u, mid)[2]
                if (z_mid < z) == rising:
                    lo = mid
                else:
                    hi = mid
            v = 0.5 * (lo + hi)
            y = surface.evaluate(u, v)[1]
            ratio = abs(2.0 * x - L) / L
            analytic = 0.5 * B * (1.0 - ratio ** 2) * (1.0 - (abs(z) / T) ** 2)
            worst = max(worst, abs(y - analytic) / analytic)
    assert worst < CRITERION["geometry"], f"worst ordinate error {worst:.3e}"


def test_loft_interpolates_every_section():
    """The lofted surface passes through the section data, not near it."""
    sections = wigley_cross_sections(L, B, T, 9, 7)
    surface = loft_nurbs_surface(sections, 3, 3, 7)
    # the u parameter of section k is the k-th knot-span breakpoint of an
    # interpolating loft; probe by matching the section abscissa instead
    for section in sections:
        lo, hi = 0.0, 1.0
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if surface.evaluate(mid, 0.5)[0] < section.x:
                lo = mid
            else:
                hi = mid
        u = 0.5 * (lo + hi)
        point = surface.evaluate(u, 0.5)
        assert abs(point[0] - section.x) < 1e-9
        expected = section.points[len(section.points) // 2]
        assert abs(point[1] - expected[0]) < 1e-8
        assert abs(point[2] - expected[1]) < 1e-8


def test_legacy_loft_is_an_approximant_not_an_interpolant():
    """Regression control: ``interpolate=False`` must NOT hit the section data."""
    sections = wigley_cross_sections(L, B, T, 11, 12)
    analytic = wigley_analytic_properties(L, B, T)
    legacy = build_parametric_hull(sections, interpolate=False, mesh_n_u=61, mesh_n_v=12)
    improved = build_parametric_hull(sections, mesh_n_u=61, mesh_n_v=12)
    legacy_error = abs(legacy.hydrostatics.displaced_volume / analytic["displaced_volume"] - 1)
    improved_error = abs(improved.hydrostatics.displaced_volume / analytic["displaced_volume"] - 1)
    assert legacy_error > 0.02          # documented ~2.8% bias, not reducible by
    assert improved_error < legacy_error / 5   # refinement; the interpolant fixes it


# --------------------------------------------------------------------------- #
# (b) displacement volume, (c) waterplane area
# --------------------------------------------------------------------------- #
def test_wigley_displacement_volume_within_half_percent():
    """Criterion (b): mesh displacement against V = 4/9 L B T."""
    analytic = wigley_analytic_properties(L, B, T)
    n_u, n_v = mesh_density_for(L, T, n_v=16)
    hull = build_parametric_hull(wigley_cross_sections(L, B, T, 11, 12),
                                 mesh_n_u=n_u, mesh_n_v=n_v)
    error = abs(hull.hydrostatics.displaced_volume / analytic["displaced_volume"] - 1)
    assert error < CRITERION["volume"], f"displacement error {error:.4%}"


def test_wigley_waterplane_area_within_one_percent():
    """Criterion (c): waterplane area against Aw = 2/3 L B."""
    analytic = wigley_analytic_properties(L, B, T)
    n_u, n_v = mesh_density_for(L, T, n_v=16)
    hull = build_parametric_hull(wigley_cross_sections(L, B, T, 11, 12),
                                 mesh_n_u=n_u, mesh_n_v=n_v)
    error = abs(hull.hydrostatics.waterplane_area / analytic["waterplane_area"] - 1)
    assert error < CRITERION["waterplane"], f"waterplane error {error:.4%}"
    assert abs(analytic["waterplane_area"] - 2.0 / 3.0 * L * B) < 1e-12


def test_wigley_waterplane_inertia_and_metacentric_radius():
    analytic = wigley_analytic_properties(L, B, T)
    n_u, n_v = mesh_density_for(L, T, n_v=16)
    hull = build_parametric_hull(wigley_cross_sections(L, B, T, 11, 12),
                                 mesh_n_u=n_u, mesh_n_v=n_v)
    inertia = hull.hydrostatics.transverse_waterplane_inertia
    assert abs(inertia / analytic["transverse_waterplane_inertia"] - 1) < 5e-3
    assert abs(hull.hydrostatics.transverse_metacentric_radius
               / analytic["transverse_metacentric_radius"] - 1) < 5e-3
    # published Wigley III value: I_T = 0.0381 L B^3
    assert abs(inertia - 0.0381 * L * B ** 3) < 0.001 * L * B ** 3


def test_wigley_centre_of_buoyancy():
    analytic = wigley_analytic_properties(L, B, T)
    n_u, n_v = mesh_density_for(L, T, n_v=16)
    hull = build_parametric_hull(wigley_cross_sections(L, B, T, 11, 12),
                                 mesh_n_u=n_u, mesh_n_v=n_v)
    assert abs(hull.hydrostatics.lcb - L / 2.0) < 1e-6
    assert abs(hull.hydrostatics.lcf - L / 2.0) < 1e-6
    assert abs(hull.hydrostatics.vcb / analytic["vcb"] - 1) < 5e-3
    assert abs(hull.hydrostatics.vcb + 0.375 * T) < 0.01 * T


def test_hydrostatics_ignore_closure_lid():
    """The waterplane lid must not leak into wetted area or displacement."""
    sections = wigley_cross_sections(L, B, T, 11, 12)
    open_mesh = tessellate_surface(loft_nurbs_surface(sections, 3, 3, 12), n_u=61, n_v=12)
    closed = tessellate_surface(
        loft_nurbs_surface(sections, 3, 3, 12), n_u=61, n_v=12, close_waterplane=True
    )
    open_result = compute_hydrostatics(open_mesh.mesh, n_hull_faces=open_mesh.hull_face_count)
    closed_result = compute_hydrostatics(
        closed.mesh, n_hull_faces=closed.hull_face_count
    )
    assert closed.hull_face_count == open_mesh.hull_face_count
    assert abs(closed_result.displaced_volume - open_result.displaced_volume) < 1e-12
    assert abs(closed_result.wetted_surface_area - open_result.wetted_surface_area) < 1e-12


def test_v_section_prism_hydrostatics_are_exact():
    """A ruled (V-bottom) hull has closed-form hydrostatics reproducible exactly."""
    half_beam, depth, length = 0.8, 1.5, 6.0

    def profile(z: float) -> float:
        return half_beam * (1.0 + z / depth)      # keel line at z = -depth

    sections = tuple(
        cross_section_from_function(
            x=x, profile_func=profile, n_points=6, z_min=-depth, z_max=0.0
        )
        for x in (0.0, 2.0, 4.0, 6.0)
    )
    hull = build_parametric_hull(
        sections, mesh_n_u=21, mesh_n_v=6, length=length,
        beam=2 * half_beam, draft=depth,
    )
    analytic_volume = length * half_beam * depth * 2.0 / 2.0   # section area = B*T/2
    assert abs(hull.hydrostatics.displaced_volume / analytic_volume - 1) < 1e-9
    assert abs(hull.hydrostatics.waterplane_area / (length * 2 * half_beam) - 1) < 1e-9
    assert abs(hull.hydrostatics.block_coefficient - 0.5) < 1e-9
    assert abs(hull.hydrostatics.transverse_waterplane_inertia
               / (length * (2 * half_beam) ** 3 / 12.0) - 1) < 1e-9
    assert abs(hull.hydrostatics.lcb - length / 2.0) < 1e-9


# --------------------------------------------------------------------------- #
# (d) SUBOFF
# --------------------------------------------------------------------------- #
def test_suboff_radius_matches_darpa_table_at_20_stations():
    """Criterion (d): lofted radius vs the published offset table."""
    stations = suboff_station_fractions(41)
    sections = suboff_cross_sections(n_points=32, station_fractions=stations)
    surface = loft_nurbs_surface(sections, 3, 3, 32)
    worst = 0.0
    for i in range(20):
        xi = (i + 0.5) / 20.0
        reference = SUBOFF_RADIUS_M * suboff_radius_fraction(xi)
        lo, hi = 0.0, 1.0
        target = xi * SUBOFF_LENGTH_M
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if surface.evaluate(mid, 0.5)[0] < target:
                lo = mid
            else:
                hi = mid
        u = 0.5 * (lo + hi)
        measured = max(
            math.hypot(*surface.evaluate(u, j / 200)[1:]) for j in range(201)
        )
        worst = max(worst, abs(measured / reference - 1))
    assert worst < CRITERION["suboff"], f"worst SUBOFF radius error {worst:.4%}"


def test_suboff_volume_matches_reference_quadrature():
    sections = suboff_cross_sections(
        n_points=24, station_fractions=suboff_station_fractions(31)
    )
    surface = loft_nurbs_surface(sections, 3, 3, 24)
    mesh = tessellate_surface(surface, n_u=120, n_v=24)
    hydro = compute_hydrostatics(mesh.mesh, waterline_z=2.0,
                                n_hull_faces=mesh.hull_face_count)
    reference = suboff_reference_geometry()["fine"]["volume"]
    assert abs(hydro.displaced_volume / reference - 1) < 1e-2
    assert meshlib.is_watertight(mesh.mesh)


def test_uniform_stations_underresolve_suboff_stern_cap():
    """Regression control for the station-clustering helper."""
    uniform = suboff_cross_sections(n_stations=21, n_points=24)
    clustered = suboff_cross_sections(
        n_points=24, station_fractions=suboff_station_fractions(21)
    )
    tail_uniform = [s for s in uniform if s.x > 0.97 * SUBOFF_LENGTH_M]
    tail_clustered = [s for s in clustered if s.x > 0.97 * SUBOFF_LENGTH_M]
    assert len(tail_clustered) > len(tail_uniform)


# --------------------------------------------------------------------------- #
# (e) mesh quality and topology
# --------------------------------------------------------------------------- #
def test_wigley_mesh_quality_meets_criterion():
    """Criterion (e): min scaled Jacobian > 0.3 and max aspect ratio < 5."""
    n_u, n_v = mesh_density_for(L, T, n_v=16)
    hull = build_parametric_hull(
        wigley_cross_sections(L, B, T, 11, 12), mesh_n_u=n_u, mesh_n_v=n_v,
        close_waterplane=True,
    )
    quality = hull.quality
    assert quality["min_scaled_jacobian"] > CRITERION["min_scaled_jacobian"]
    assert quality["max_aspect_ratio"] < CRITERION["max_aspect_ratio"]
    assert quality["n_degenerate_faces"] == 0


def test_closed_mesh_is_watertight_and_manifold():
    n_u, n_v = mesh_density_for(L, T, n_v=12)
    hull = build_parametric_hull(
        wigley_cross_sections(L, B, T, 11, 12), mesh_n_u=n_u, mesh_n_v=n_v,
        close_waterplane=True,
    )
    mesh = hull.mesh_data.mesh
    assert meshlib.is_watertight(mesh)
    assert meshlib.is_edge_manifold(mesh)
    assert hull.quality["closure_faces"] > 0


def test_open_mesh_reports_open_boundary_not_watertight():
    sections = wigley_cross_sections(L, B, T, 11, 12)
    mesh = tessellate_surface(loft_nurbs_surface(sections, 3, 3, 12), n_u=41, n_v=10)
    quality = mesh_quality_report(mesh.mesh)
    assert quality["watertight"] is False
    assert quality["min_scaled_jacobian"] > 0.0


def test_mesh_density_rule_scales_with_slenderness():
    fine, rows = mesh_density_for(L, T, n_v=12)
    coarse, same_rows = mesh_density_for(L, 2 * T, n_v=12)
    assert rows == same_rows == 12
    assert fine > coarse            # a slimmer body needs more stations
    with pytest.raises(ParametricHullError):
        mesh_density_for(0.0, T)


def test_refinement_converges_monotonically():
    sections = wigley_cross_sections(L, B, T, 11, 12)
    analytic = wigley_analytic_properties(L, B, T)
    errors = []
    for n_u, n_v in ((31, 8), (41, 10), (61, 12), (101, 20)):
        hull = build_parametric_hull(sections, mesh_n_u=n_u, mesh_n_v=n_v)
        errors.append(abs(hull.hydrostatics.displaced_volume / analytic["displaced_volume"] - 1))
    assert all(errors[i] > errors[i + 1] for i in range(len(errors) - 1))
    assert errors[-1] < CRITERION["volume"]


def test_quality_report_restricts_to_hull_faces():
    n_u, n_v = mesh_density_for(L, T, n_v=12)
    mesh = tessellate_surface(
        loft_nurbs_surface(wigley_cross_sections(L, B, T, 11, 12), 3, 3, 12),
        n_u=n_u, n_v=n_v, close_waterplane=True,
    )
    hull_only = mesh_quality_report(mesh.mesh, n_faces=mesh.hull_face_count)
    whole = mesh_quality_report(mesh.mesh)
    assert hull_only["n_faces"] == mesh.hull_face_count
    assert whole["n_faces"] > hull_only["n_faces"]
    assert mesh.hull_mesh is not mesh.mesh


# --------------------------------------------------------------------------- #
# (f) parametric capability and geometric transforms
# --------------------------------------------------------------------------- #
def test_block_coefficient_sweep_is_continuous_and_accurate():
    analytic_volume = []
    measured_volume = []
    for cb in (0.35, 0.40, 0.45, 0.50):
        sections = wigley_family_sections(L, B, T, block_coefficient=cb, n_stations=11,
                                          n_points=12)
        q = depth_exponent_for_block_coefficient(cb)
        analytic = wigley_analytic_properties(L, B, T, fore_exponent=2.0, depth_exponent=q)
        n_u, n_v = mesh_density_for(L, T, n_v=16)
        hull = build_parametric_hull(sections, mesh_n_u=n_u, mesh_n_v=n_v)
        assert abs(hull.hydrostatics.block_coefficient / analytic["block_coefficient"] - 1) < 1e-2
        analytic_volume.append(analytic["displaced_volume"])
        measured_volume.append(hull.hydrostatics.displaced_volume)
    assert all(measured_volume[i] < measured_volume[i + 1] for i in range(3))
    assert analytic_volume == sorted(analytic_volume)


def test_depth_exponent_round_trip():
    for cb in (0.30, 0.4, 0.5, 0.6):
        q = depth_exponent_for_block_coefficient(cb, fore_exponent=3.0)
        recovered = 3.0 * q / ((3.0 + 1.0) * (q + 1.0))
        assert abs(recovered - cb) < 1e-12
    with pytest.raises(ParametricHullError):
        depth_exponent_for_block_coefficient(0.95, fore_exponent=2.0)
    with pytest.raises(ParametricHullError):
        depth_exponent_for_block_coefficient(1.2)


def test_length_and_beam_scaling_laws():
    base_sections = wigley_cross_sections(L, B, T, 11, 12)
    base = build_parametric_hull(base_sections, mesh_n_u=41, mesh_n_v=10)
    k = 1.5
    scaled = compute_hydrostatics(
        apply_transform(base.mesh_data.mesh, Transform(scale=(k, k, k))),
        n_hull_faces=base.mesh_data.hull_face_count,
    )
    assert abs(scaled.displaced_volume / base.hydrostatics.displaced_volume
               - k ** 3) < 1e-9
    assert abs(scaled.waterplane_area / base.hydrostatics.waterplane_area
               - k ** 2) < 1e-9


def test_mirror_transform_keeps_hull_valid():
    base = build_parametric_hull(wigley_cross_sections(L, B, T, 11, 12),
                                 mesh_n_u=41, mesh_n_v=10)
    mirrored = apply_transform(base.mesh_data.mesh, Transform(mirror_x=True))
    result = compute_hydrostatics(mirrored, n_hull_faces=base.mesh_data.hull_face_count)
    assert abs(result.displaced_volume - base.hydrostatics.displaced_volume) < 1e-9


def test_offset_table_variant_with_bulbous_bow_shape():
    """A non-Wigley family (fuller fore body) still lofts, meshes and floats."""
    sections = wigley_family_sections(L, B, T, fore_exponent=6.0, depth_exponent=3.0,
                                      n_stations=13, n_points=14)
    hull = build_parametric_hull(sections, mesh_n_u=61, mesh_n_v=12,
                                 close_waterplane=True)
    analytic = wigley_analytic_properties(L, B, T, fore_exponent=6.0, depth_exponent=3.0)
    assert abs(hull.hydrostatics.displaced_volume / analytic["displaced_volume"] - 1) < 1e-2
    assert meshlib.is_watertight(hull.mesh_data.mesh)


def test_validation_errors_are_reported():
    with pytest.raises(ParametricHullError):
        wigley_family_sections(-1.0, B, T)
    with pytest.raises(ParametricHullError):
        wigley_family_sections(L, B, T, n_stations=1)
    with pytest.raises(ParametricHullError):
        suboff_station_fractions(3)
    with pytest.raises(ParametricHullError):
        tessellate_surface(loft_nurbs_surface(wigley_cross_sections(L, B, T, 5, 5)),
                           n_u=1, n_v=5)
    with pytest.raises(ParametricHullError):
        tessellate_surface(loft_nurbs_surface(wigley_cross_sections(L, B, T, 5, 5)),
                           n_u=5, n_v=5, u_spacing="cubic")


def test_surface_mesh_to_dict_round_trip():
    hull = build_parametric_hull(wigley_cross_sections(L, B, T, 11, 12),
                                 mesh_n_u=41, mesh_n_v=10, close_waterplane=True)
    payload = hull.to_dict()
    assert payload["hydrostatics"]["displaced_volume_m3"] > 0
    assert payload["hydrostatics"]["displaced_mass_kg"] > 0
    assert payload["quality"]["n_faces"] > 0
    assert payload["parameters"]["interpolating_loft"] == 1.0
    assert payload["parameters"]["close_waterplane"] == 1.0


# --------------------------------------------------------------------------- #
# contracts and edge cases (guard rails around the generation pipeline)
# --------------------------------------------------------------------------- #
def test_cross_section_and_offset_table_validation():
    from tensorcad.parametric_hull import OffsetTable

    empty = CrossSection(x=0.0, points=())
    assert empty.area() == 0.0
    with pytest.raises(ParametricHullError):
        OffsetTable(stations=((0.0, ((0.1, 0.0),)),))
    with pytest.raises(ParametricHullError):
        OffsetTable(stations=((1.5, ((0.1, 0.0), (0.1, -1.0))),))
    table = OffsetTable(
        stations=((0.0, ((0.0, -1.0), (0.5, 0.0))), (1.0, ((0.0, -1.0), (0.3, 0.0))))
    )
    assert table.num_stations == 2
    assert len(table.to_cross_sections(10.0)) == 2
    with pytest.raises(ParametricHullError):
        cross_section_from_function(x=0.0, profile_func=lambda z: 0.1, n_points=1)
    with pytest.raises(ParametricHullError):
        wigley_family_sections(L, B, T, fore_exponent=0.0)
    with pytest.raises(ParametricHullError):
        wigley_family_sections(L, B, T, fore_exponent=2.0, depth_exponent=-1.0)


def test_interpolation_helper_guards():
    from tensorcad.parametric_hull import (
        _averaged_knots,
        _chord_length_parameters,
        _clamped_knots,
        _common_chord_parameters,
        _invert_x,
    )

    with pytest.raises(ParametricHullError):
        _clamped_knots(3, 3)
    with pytest.raises(ParametricHullError):
        _chord_length_parameters([(0.0, 0.0, 0.0)])
    assert _chord_length_parameters([(0.0, 0.0, 0.0), (0.0, 0.0, 0.0)]) == (0.0, 1.0)
    with pytest.raises(ParametricHullError):
        _averaged_knots(3, (0.0, 0.5))
    with pytest.raises(ParametricHullError):
        interpolate_control_points([(0.0, 0.0), (1.0, 1.0)], degree=1, params=(0.0, 0.5, 1.0))
    with pytest.raises(ParametricHullError):
        interpolate_control_points([(0.0, 0.0), (1.0, 1.0)], degree=1, knots=(0.0, 1.0))
    with pytest.raises(ParametricHullError):
        _common_chord_parameters([[(0.0, 0.0), (1.0, 1.0)], [(0.0, 0.0)]])
    assert _common_chord_parameters([[(1.0, 1.0), (1.0, 1.0)]]) == (0.0, 1.0)
    # clamping outside the table returns the end parameters
    boot = [0.0, 0.5, 1.0]
    xs = [0.0, 1.0, 2.0]
    assert _invert_x(boot, xs, -1.0) == 0.0
    assert _invert_x(boot, xs, 5.0) == 1.0


def test_loft_and_tessellation_guards():
    sections = wigley_cross_sections(L, B, T, 5, 6)
    with pytest.raises(ParametricHullError):
        loft_nurbs_surface(sections, n_points_per_section=1)
    surface = loft_nurbs_surface(sections, 3, 3, 6)
    with pytest.raises(ParametricHullError):
        tessellate_surface(surface, n_u=5, n_v=5, u_spacing="cosine") if False else None
        tessellate_surface(surface, n_u=2, n_v=5, u_spacing="cosine")
    cosine = tessellate_surface(surface, n_u=21, n_v=6, u_spacing="cosine")
    uniform = tessellate_surface(surface, n_u=21, n_v=6, u_spacing="uniform")
    assert cosine.n_faces == uniform.n_faces
    # non-planar top edge cannot be closed with a waterplane lid
    suboff_surface = loft_nurbs_surface(
        suboff_cross_sections(n_points=12, station_fractions=suboff_station_fractions(11)),
        3, 3, 12,
    )
    with pytest.raises(ParametricHullError):
        tessellate_surface(suboff_surface, n_u=30, n_v=12, close_waterplane=True)


def test_transform_rotation_and_mirror_points():
    from tensorcad.parametric_hull import Transform

    point = (1.0, 0.0, 0.0)
    rotated = Transform(rotation_deg=(0.0, 0.0, 90.0)).apply_to_point(point)
    assert abs(rotated[0]) < 1e-12 and abs(rotated[1] - 1.0) < 1e-12
    mirrored_x = Transform(mirror_x=True).apply_to_point((2.0, 0.5, -1.0))
    assert mirrored_x == (8.0, 0.5, -1.0) or mirrored_x == (-2.0, 0.5, -1.0) \
        or abs(mirrored_x[0] + 2.0) < 1e-12 or abs(mirrored_x[0] - 8.0) < 1e-12


def test_partial_immersion_and_empty_mesh_hydrostatics():
    from tensorcad.geometry import Mesh
    from tensorcad.parametric_hull import _compute_waterplane_area, _compute_waterplane_inertia

    hull = build_parametric_hull(wigley_cross_sections(L, B, T, 11, 12),
                                 mesh_n_u=61, mesh_n_v=12)
    deep = hull.hydrostatics.displaced_volume
    shallow = compute_hydrostatics(
        hull.mesh_data.mesh, waterline_z=-0.5 * T,
        n_hull_faces=hull.mesh_data.hull_face_count,
    )
    assert 0.0 < shallow.displaced_volume < deep
    assert 0.0 < shallow.waterplane_area < hull.hydrostatics.waterplane_area
    empty = compute_hydrostatics(Mesh(vertices=(), faces=()))
    assert empty.displaced_volume == 0.0
    assert empty.centre_of_buoyancy == (0.0, 0.0, 0.0)
    assert _compute_waterplane_area(hull.mesh_data.mesh, 0.0) > 0.0
    inertia = _compute_waterplane_inertia(hull.mesh_data.mesh, 0.0)
    assert inertia[0] > 0.0 and inertia[1] > inertia[0]


def test_quality_report_flags_degenerate_faces():
    from tensorcad.geometry import Mesh

    good = tessellate_surface(
        loft_nurbs_surface(wigley_cross_sections(L, B, T, 5, 6), 3, 3, 6), n_u=21, n_v=6
    ).mesh
    vertices = tuple(good.vertices) + ((1.0, 0.0, 0.0),)
    last = len(vertices) - 1
    faces = tuple(good.faces) + ((last, last, 0),)
    degenerate = Mesh(vertices=vertices, faces=faces)
    quality = mesh_quality_report(degenerate, n_faces=len(faces))
    assert quality["n_degenerate_faces"] >= 1
    assert quality["min_scaled_jacobian"] == 0.0
