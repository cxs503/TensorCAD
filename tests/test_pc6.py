"""Acceptance tests for the PC6 geometry service.

The suite is written as evidence, not as a smoke test: every assertion is a
machine-checkable criterion from the campaign ticket --

1. the hull is a closed manifold (watertight, edge-manifold, Euler 2, one shell);
2. the sampled solid reproduces L, B, T and the analytic displacement;
3. a beam sweep 24 -> 28 m regenerates a valid, watertight hull each time;
4. a non-positive draft is rejected;
5. a self-intersecting hull is detected by geometry, not only by the parameter
   range check;
6. the implicit field and the STL mesh describe the same solid;
7. the published contract round-trips through ``load_contract``;
8. the ice-channel corridor is a simple plan shape and a closed prism.
"""
from __future__ import annotations

import itertools
import json
import math
from dataclasses import replace
from pathlib import Path

import pytest

from tensorcad import mesh as meshlib
from tensorcad.icechannel import (
    DE_LONG_PEVEK_ROUTE,
    build_channel_mesh,
    channel_measurements,
    great_circle_km,
    ice_channel,
)
from tensorcad.intersect import find_self_intersections, triangles_intersect
from tensorcad.pc6 import (
    PC6Parameters,
    analytic_form,
    beam_sweep,
    build_mesh,
    field_value,
    generate,
    half_breadth,
    inside_by_winding_number,
    longitudinal_ratio,
    section_ratio,
    stem_x,
    topology,
)
from tensorcad.pc6_contract import build_contract, load_contract

# A deliberately coarse grid: every acceptance criterion below is either
# topological (grid-independent) or compared against the analytic reference, so a
# small mesh keeps the suite fast without weakening any assertion.
COARSE = replace(PC6Parameters(), longitudinal_stations=41, vertical_levels=11)
#: The coarse grid is a suite-speed choice, so the displacement check is given the
#: tolerance that grid actually meets; the shipped grid meets the shipped 0.5 %.
COARSE_TOLERANCE = 0.01


def test_hull_is_closed_manifold() -> None:
    model = generate(COARSE, check_self_intersection=False, volume_tolerance=COARSE_TOLERANCE)
    report = model.report["topology"]
    assert report["watertight"] is True
    assert report["edge_manifold"] is True
    assert report["boundary_edges"] == 0
    assert report["euler_characteristic"] == 2
    assert report["components"] == 1
    assert report["zero_area_faces"] == 0
    assert report["vertex_link_defects"] == 0
    # outward orientation: the divergence-theorem volume is positive and matches
    # the enclosed solid rather than its negative
    assert report["signed_volume_m3"] > 0.0
    assert model.validation.passed
    for check in model.validation.checks:
        assert check.passed, f"{check.name}: {check.detail}"


def test_principal_dimensions_and_analytic_hydrostatics() -> None:
    p = PC6Parameters()
    form = analytic_form(p)
    model = generate(p, check_self_intersection=False)
    sizes = model.report["principal_dimensions"]
    assert sizes["length"] == pytest.approx(p.length, rel=1e-12)
    assert sizes["beam"] == pytest.approx(p.beam, rel=1e-12)
    assert sizes["draft"] == pytest.approx(p.draft, rel=1e-12)
    # the analytic form is a plausible PC6 hull, not a degenerate one
    assert 0.70 < form.block_coefficient < 0.86
    assert form.midship_coefficient > 0.85
    assert form.prismatic_coefficient > form.block_coefficient
    assert form.waterplane_coefficient > form.block_coefficient
    # sampled integrals against the closed form (independent code path)
    hydro = model.hydrostatics
    assert hydro.displacement_volume == pytest.approx(form.volume, rel=5e-3)
    assert hydro.waterplane_area == pytest.approx(form.waterplane_area, rel=5e-3)
    assert hydro.longitudinal_centre_buoyancy == pytest.approx(form.lcb, rel=5e-3)
    assert hydro.vertical_centre_buoyancy == pytest.approx(form.vcb, rel=5e-3)
    assert hydro.block_coefficient == pytest.approx(form.block_coefficient, rel=5e-3)


def test_section_and_longitudinal_ratios_are_consistent() -> None:
    p = PC6Parameters()
    # g: transom value, unity at the maximum section, zero at the stem
    assert longitudinal_ratio(p, 0.0) == pytest.approx(p.transom_ratio)
    assert longitudinal_ratio(p, p.midship_fraction) == pytest.approx(1.0)
    assert longitudinal_ratio(p, 1.0) == 0.0
    # sec: flat of bottom at the keel, unity at the design waterline
    assert section_ratio(p, 1.0) == pytest.approx(1.0)
    assert section_ratio(p, 0.0) == pytest.approx(p.bilge_ratio / (1.0 + p.flare))
    assert half_breadth(p, p.midship_fraction, 1.0) == pytest.approx(0.5 * p.beam)
    # monotone increase of the half-breadth with height (no reentrant section)
    levels = [section_ratio(p, j / 40.0) for j in range(41)]
    assert all(b >= a - 1e-15 for a, b in itertools.pairwise(levels))
    # the stem is raked forward: it is further forward at the DWL than at the keel
    assert stem_x(p, 0.0) == pytest.approx(p.length)
    assert stem_x(p, -p.draft) == pytest.approx(p.length - p.stem_rake)


def test_beam_sweep_24_to_28_regenerates() -> None:
    models = beam_sweep((24.0, 26.0, 28.0), COARSE, volume_tolerance=COARSE_TOLERANCE)
    volumes = []
    for model in models:
        p = model.parameters
        assert model.report["topology"]["watertight"] is True
        assert model.report["topology"]["euler_characteristic"] == 2
        assert model.report["principal_dimensions"]["beam"] == pytest.approx(p.beam, rel=1e-12)
        assert model.validation.passed
        assert abs(model.report["reference_volume_error"]) < COARSE_TOLERANCE
        volumes.append(model.hydrostatics.displacement_volume)
    # a wider hull displaces strictly more at the same length and draft
    assert volumes[0] < volumes[1] < volumes[2]
    # the waterplane area scales with the beam: it is exactly B * L * int(g)
    waterplane_gain = (
        models[2].hydrostatics.waterplane_area / models[0].hydrostatics.waterplane_area
    )
    assert waterplane_gain == pytest.approx(models[2].parameters.beam / models[0].parameters.beam, rel=5e-3)


@pytest.mark.parametrize("draft", [0.0, -4.0])
def test_non_positive_draft_is_rejected(draft: float) -> None:
    with pytest.raises(ValueError) as excinfo:
        generate(replace(PC6Parameters(), draft=draft))
    assert "draft" in str(excinfo.value)


def test_degenerate_principal_ratios_are_rejected() -> None:
    with pytest.raises(ValueError, match="beam/draft"):
        generate(replace(PC6Parameters(), draft=2.0))  # B/T = 12: not a PC6 hull
    with pytest.raises(ValueError, match="length/beam"):
        generate(replace(PC6Parameters(), beam=60.0))  # L/B = 3.3: a barge
    with pytest.raises(ValueError, match="stem_rake"):
        generate(replace(PC6Parameters(), stem_rake=0.4 * PC6Parameters().length))


def test_self_intersecting_hull_is_detected() -> None:
    # the primitive itself
    flat = [(0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (0.0, 2.0, 0.0)]
    piercing = [(0.5, 0.5, -1.0), (0.5, 0.5, 1.0), (1.5, 1.5, 0.5)]
    assert triangles_intersect(flat, piercing)
    assert not triangles_intersect(flat, [(5.0, 5.0, 0.0), (6.0, 5.0, 0.0), (5.0, 6.0, 0.0)])
    # a clean hull has none
    assert find_self_intersections(build_mesh(COARSE)) == ()
    # a stem rake longer than the hull makes the loft fold through the transom;
    # the raw loft bypasses the parameter check so the *geometric* detector is
    # what has to catch it
    broken = replace(COARSE, stem_rake=1.25 * COARSE.length)
    raw = build_mesh(broken, validate=False)
    pairs = find_self_intersections(raw)
    assert pairs, "the self-intersection detector missed a folded loft"
    # and the published generator refuses the same parameters outright
    with pytest.raises(ValueError, match="stem_rake"):
        generate(broken)


def test_implicit_field_agrees_with_mesh() -> None:
    p = COARSE
    mesh = build_mesh(p)
    points = [
        (0.5 * p.length, 0.0, -0.5 * p.draft),  # clearly inside
        (0.5 * p.length, 0.0, 0.5 * p.draft),  # above the design waterline
        (0.5 * p.length, 0.0, -1.5 * p.draft),  # below the keel
        (-0.1 * p.length, 0.0, -0.5 * p.draft),  # aft of the transom
        (1.1 * p.length, 0.0, -0.5 * p.draft),  # forward of the stem
        (0.5 * p.length, 0.9 * p.beam, -0.5 * p.draft),  # outside the side shell
    ]
    for point in points:
        by_field = field_value(p, *point) < 0.0
        assert by_field == inside_by_winding_number(mesh, point), point
    # the sign changes across the shell: walking the half-breadth crosses it once
    s, u = 0.5, 0.5
    z = -p.draft * (1.0 - u)
    breadth = half_breadth(p, s, u)
    assert field_value(p, s * stem_x(p, z), breadth * 0.5, z) < 0.0
    assert field_value(p, s * stem_x(p, z), breadth * 1.5, z) > 0.0


def test_contract_round_trip(tmp_path: Path) -> None:
    contract = build_contract(
        tmp_path, COARSE, sweep_beams=(24.0, 28.0), field_samples=6,
        volume_tolerance=COARSE_TOLERANCE,
    )
    assert contract["acceptance"]["watertight_closed_manifold"] is True
    assert contract["acceptance"]["beam_sweep_24_28_regenerates"] is True
    assert contract["acceptance"]["negative_cases_detected"] is True
    assert contract["acceptance"]["implicit_field_agrees_with_mesh"] is True
    written = tmp_path / "geom_pc6.json"
    assert written.is_file()
    reloaded = json.loads(written.read_text(encoding="utf-8"))
    assert reloaded["schema"] == contract["schema"]
    assert reloaded["digest"]["value"] == contract["digest"]["value"]
    geometry = load_contract(written)
    assert geometry.parameters.beam == pytest.approx(COARSE.beam)
    assert meshlib.is_watertight(geometry.hull)
    assert meshlib.is_watertight(geometry.channel)
    recorded = contract["geometry"]["hydrostatics"]["mesh"]["displacement_volume_m3"]
    assert geometry.hydrostatics["displacement_volume_m3"] == pytest.approx(recorded, rel=1e-9)
    # the reloaded field still answers membership correctly
    assert geometry.inside(0.5 * COARSE.length, 0.0, -0.5 * COARSE.draft)
    assert not geometry.inside(0.5 * COARSE.length, 0.0, 0.5 * COARSE.draft)
    # the sampled field in the contract is the same function the code evaluates
    sampled = contract["implicit_field"]["hull"]
    nx, ny, nz = sampled["grid"]

    def axis(axis_name: str, position: int, count: int) -> float:
        low, high = sampled["bounds"][axis_name]
        return low + (high - low) * position / (count - 1)

    index = 0
    for iz in range(nz):
        for iy in range(ny):
            for ix in range(nx):
                point = (axis("x", ix, nx), axis("y", iy, ny), axis("z", iz, nz))
                assert sampled["values"][index] == pytest.approx(
                    field_value(COARSE, *point), abs=1e-3
                )
                index += 1


def test_contract_rejects_tampered_artifact(tmp_path: Path) -> None:
    contract = build_contract(
        tmp_path, COARSE, sweep_beams=(24.0,), field_samples=3,
        volume_tolerance=COARSE_TOLERANCE,
    )
    written = tmp_path / "geom_pc6.json"
    tampered = json.loads(written.read_text(encoding="utf-8"))
    tampered["geometry"]["hydrostatics"]["mesh"]["displacement_volume_m3"] *= 2.0
    written.write_text(json.dumps(tampered), encoding="utf-8")
    with pytest.raises(ValueError, match="deviates"):
        load_contract(written)
    # a missing STL is a loud failure, not a silent partial load
    (tmp_path / contract["artifacts"]["hull_stl"]).unlink()
    with pytest.raises(FileNotFoundError):
        load_contract(written)


def test_ice_channel_geometry() -> None:
    channel = ice_channel()
    assert len(channel.centre) == len(DE_LONG_PEVEK_ROUTE)
    # segment lengths are the great-circle distances, recomputed here directly
    for index, distance in enumerate(channel.segment_km):
        _, lat1, lon1 = DE_LONG_PEVEK_ROUTE[index]
        _, lat2, lon2 = DE_LONG_PEVEK_ROUTE[index + 1]
        assert distance == pytest.approx(great_circle_km(lat1, lon1, lat2, lon2), rel=1e-12)
    assert channel.total_km == pytest.approx(sum(channel.segment_km), rel=1e-12)
    assert 800.0 < channel.total_km < 1400.0  # de Long -> Pevek, campaign scale
    # the local plane preserves the along-track distance of the full route
    plan = sum(
        math.dist(channel.centre[i], channel.centre[i + 1]) for i in range(len(channel.centre) - 1)
    )
    assert plan / 1000.0 == pytest.approx(channel.total_km, rel=2e-3)
    # walls stay at half width on straight runs
    straight = [
        math.dist(channel.left[i], channel.right[i]) for i in (0, 1, 2, len(channel.centre) - 1)
    ]
    for width in straight:
        assert width == pytest.approx(channel.width, rel=1e-6)

    mesh = build_channel_mesh(channel)
    measurement = channel_measurements(mesh, channel)
    assert measurement["watertight"] is True
    assert measurement["edge_manifold"] is True
    assert measurement["euler_characteristic"] == 2
    assert measurement["signed_volume_m3"] > 0.0
    assert measurement["self_intersection_pairs"] == 0
    assert find_self_intersections(mesh) == ()
    # the trimmed prism volume equals plan length x width x depth
    assert measurement["volume_relative_error"] == pytest.approx(0.0, abs=1e-6)


def test_ice_channel_rejects_bad_input() -> None:
    with pytest.raises(ValueError):
        ice_channel(DE_LONG_PEVEK_ROUTE[:1])
    with pytest.raises(ValueError, match="width"):
        ice_channel(width=0.0)
    with pytest.raises(ValueError, match="depth"):
        ice_channel(depth=-1.0)
    with pytest.raises(ValueError, match="out-of-range"):
        ice_channel((("nowhere", 120.0, 0.0), ("pevek", 69.7, 170.3)))


def test_topology_only_reports_what_it_measures() -> None:
    mesh = build_mesh(COARSE)
    report = topology(mesh)
    assert report["faces"] == len(mesh.faces)
    assert report["edges"] == len(meshlib.edge_face_counts(mesh))
    assert report["euler_characteristic"] == (
        report["vertices"] - report["edges"] + report["faces"]
    )
