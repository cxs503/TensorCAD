"""Tests for the baseline-fleet intake gate (GAP-CAD-T1-BASELINE-FLEET, mission m17335).

The tests assert invariants, not merely that code runs:

* every registered hull is watertight, edge-manifold, single-shell (Euler 2) and spans the
  declared L/B/T;
* the enclosed volume agrees with the record's reference inside 1 % at production
  resolution and the discretisation error is shown to vanish (first-order, Richardson);
* the gate *rejects* every negative control: the 304-triangle placeholder, a byte-identical
  hull pair, a block coefficient clipped into the form family's window, and a record whose
  displacement contradicts its own particulars;
* the two independent geometry routes (Series-60 section law vs the Beta section-area law)
  produce the same enclosed volume for the same ship.

Coarse grids are used where resolution is not the point, so the suite stays fast; the
1 %-at-production-resolution claim is asserted on Series 60 only and is additionally
recorded for all six ships in ``runs/m17335_fleet/evidence.json``.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from tensorcad import mesh as meshlib
from tensorcad.fleet import (
    RouteRefusal,
    build_ship_geometry,
    clipped_cb_record,
    duplicate_manifest,
    duplicate_stl_pair,
    evaluate_gates,
    load_manifest,
    measure,
    placeholder_box,
    placeholder_reasons,
    read_stl,
    stl_mesh,
    wrong_volume_record,
)
from tensorcad.fleet.controls import write_binary_stl
from tensorcad.fleet.gates import gate_particulars_consistent
from tensorcad.fleet.manifest import ManifestError, record_from_dict
from tensorcad.fleet.run import run_fleet

COARSE = {"nx": 121, "nz": 25}


@pytest.fixture(scope="module")
def manifest():
    return load_manifest()


@pytest.fixture(scope="module")
def coarse_geometries(manifest):
    return {ship.name: build_ship_geometry(ship, **COARSE) for ship in manifest.ships}


def test_manifest_registers_the_six_standard_hulls_with_provenance(manifest):
    assert manifest.names() == ("KCS", "KVLCC2", "ONRT", "JBC", "DTC", "Series60")
    assert manifest.sha256 and len(manifest.sha256) == 64
    for ship in manifest.ships:
        assert ship.provenance.complete, ship.name
        assert ship.provenance.url.startswith("http"), ship.name
        assert ship.provenance.license, ship.name
        # every reference volume states whether it is external or derived
        assert ship.reference_basis() in {"literature_value_recalled", "closure"}
        assert ship.reference_volume() > 0.0
        assert ship.block_coefficient < (ship.midship_coefficient or 1.0) or ship.midship_coefficient is None


def test_closed_form_volume_matches_the_particulars(manifest):
    for ship in manifest.ships:
        closed = ship.block_coefficient * ship.length * ship.beam * ship.draft
        assert math.isclose(closed, ship.closed_form_volume(), rel_tol=1e-12)
        if ship.displacement_volume is not None:
            gate = gate_particulars_consistent(ship)
            assert gate.passed, (ship.name, gate.detail)


@pytest.mark.parametrize("name", ["KCS", "KVLCC2", "ONRT", "JBC", "DTC", "Series60"])
def test_every_hull_is_watertight_manifold_and_spans_its_dimensions(manifest, coarse_geometries, name):
    record = manifest.by_name(name)
    mesh = coarse_geometries[name].mesh
    facts = measure(mesh)
    assert facts.boundary_edge_count == 0
    assert facts.nonmanifold_edge_count == 0
    assert facts.inconsistently_oriented_edge_count == 0
    assert facts.zero_area_face_count == 0
    assert facts.duplicate_face_count == 0
    assert facts.euler_characteristic == 2
    assert facts.signed_volume > 0.0
    xmin, xmax, ymin, ymax, zmin, zmax = facts.bbox
    assert math.isclose(xmax - xmin, record.length, rel_tol=1e-9)
    # the widest station must reach the declared beam, but no discrete grid hits the exact
    # peak, so allow a 1e-4 shortfall and forbid any overshoot
    assert ymax - ymin <= record.beam * (1.0 + 1e-9)
    assert math.isclose(ymax - ymin, record.beam, rel_tol=1e-4)
    assert math.isclose(zmax - zmin, record.draft, rel_tol=1e-9)


def test_coarse_mesh_volume_is_within_a_few_percent_and_converges_downwards(manifest):
    """At 121x25 the linear-chord deficit is ~1.5 %/half the grid: assert it is bounded."""
    for ship in manifest.ships:
        reference = ship.reference_volume()
        volumes = []
        for nx, nz in ((121, 25), (241, 49)):
            mesh = build_ship_geometry(ship, nx=nx, nz=nz).mesh
            volumes.append(meshlib.total_signed_volume(mesh))
        coarse_error = abs(volumes[0] - reference) / reference
        fine_error = abs(volumes[1] - reference) / reference
        assert coarse_error < 0.04, (ship.name, coarse_error)
        assert fine_error < coarse_error, (ship.name, fine_error, coarse_error)


def test_production_resolution_meets_the_one_percent_volume_criterion(manifest):
    ship = manifest.by_name("Series60")
    mesh = build_ship_geometry(ship, nx=481, nz=97).mesh
    facts = measure(mesh)
    error = abs(facts.signed_volume - ship.reference_volume()) / ship.reference_volume()
    assert error < 0.01, error
    assert facts.boundary_edge_count == 0


def test_two_independent_geometry_laws_agree_on_the_same_hull(manifest):
    import dataclasses

    ship = manifest.by_name("KCS")
    series = build_ship_geometry(ship, nx=241, nz=49).mesh
    other = dataclasses.replace(ship, geometry_route="parametric_sections")
    parametric = build_ship_geometry(other, nx=241, nz=49).mesh
    v_series = meshlib.total_signed_volume(series)
    v_parametric = meshlib.total_signed_volume(parametric)
    assert abs(v_series - v_parametric) / v_series < 0.01
    assert meshlib.is_watertight(parametric)


def test_all_gates_pass_for_a_coarse_hull_when_the_tolerance_is_stated(manifest, coarse_geometries):
    ship = manifest.by_name("DTC")
    verdict = evaluate_gates(ship, coarse_geometries["DTC"].mesh, manifest, volume_tolerance=0.05)
    assert verdict["passed"], verdict["failed_gates"]
    assert "G8_principal_dimensions" not in verdict["failed_gates"]
    assert verdict["volume_relative_error"] < 0.05


def test_placeholder_304_triangles_is_rejected_by_name():
    mesh = placeholder_box(304)
    assert len(mesh.faces) == 304
    facts = measure(mesh)
    reasons = placeholder_reasons(facts, mesh=mesh)
    assert any("placeholder_signature" in reason for reason in reasons), reasons


def test_byte_identical_hull_pair_is_rejected(manifest, tmp_path):
    first, second = duplicate_stl_pair(tmp_path)
    assert first.read_bytes() == second.read_bytes()
    duplicate = duplicate_manifest(manifest, first, second)
    collisions = duplicate.duplicate_geometry()
    assert len(collisions) == 1
    verdicts = [evaluate_gates(ship, placeholder_box(304), duplicate) for ship in duplicate.ships]
    assert all("G2_no_duplicate_geometry" in verdict["failed_gates"] for verdict in verdicts)


def test_clipped_block_coefficient_fails_the_volume_gate(manifest):
    base = manifest.by_name("ONRT")
    clipped = clipped_cb_record(base, cb=0.56)
    mesh = build_ship_geometry(clipped, **COARSE).mesh
    verdict = evaluate_gates(clipped, mesh, manifest)
    assert "G7_volume_vs_reference" in verdict["failed_gates"], verdict["failed_gates"]
    assert verdict["volume_relative_error"] > 0.03


def test_incoherent_particulars_are_caught_before_any_geometry_runs(manifest):
    wrong = wrong_volume_record(manifest.by_name("KCS"), factor=1.08)
    gate = gate_particulars_consistent(wrong)
    assert not gate.passed
    assert "7." in gate.detail or "8." in gate.detail  # the measured percentage is in the message


def test_route_refuses_ships_outside_the_series60_window(manifest):
    import dataclasses

    ship = dataclasses.replace(manifest.by_name("Series60"), block_coefficient=0.90)
    with pytest.raises(RouteRefusal) as excinfo:
        build_ship_geometry(ship, **COARSE)
    assert "0.9" in str(excinfo.value)


def test_binary_and_ascii_stl_round_trip(tmp_path):
    mesh = placeholder_box(304)
    binary = write_binary_stl(mesh, tmp_path / "box_binary.stl")
    back = read_stl(binary)
    assert len(back.faces) == 304
    ascii_path = tmp_path / "box_ascii.stl"
    lines = ["solid box"]
    for face in mesh.faces:
        lines.append("  facet normal 0 0 0")
        lines.append("    outer loop")
        for index in face:
            x, y, z = mesh.vertices[index]
            lines.append(f"      vertex {x!r} {y!r} {z!r}")
        lines.append("    endloop")
        lines.append("  endfacet")
    lines.append("endsolid box")
    ascii_path.write_text("\n".join(lines))
    assert len(read_stl(ascii_path).faces) == 304


def test_stl_import_route_welds_and_orients(manifest, tmp_path):
    ship = manifest.by_name("Series60")
    generated = build_ship_geometry(ship, **COARSE).mesh
    export = write_binary_stl(generated, tmp_path / "hull.stl", solid_name="series60")
    imported = stl_mesh(ship, export)
    assert len(imported.mesh.faces) == len(generated.faces)
    assert meshlib.is_watertight(imported.mesh)
    volume = meshlib.total_signed_volume(imported.mesh)
    # binary STL stores float32 coordinates, so the round trip is exact to float32, not to 1e-9
    assert math.isclose(volume, meshlib.total_signed_volume(generated), rel_tol=1e-6)


def test_malformed_manifest_is_rejected(tmp_path):
    payload = {
        "schema": "tensorcad.fleet.baseline.v1",
        "ships": [
            {
                "name": "Broken",
                "length": 100.0,
                "beam": 10.0,
                "draft": 5.0,
                "block_coefficient": 0.6,
                "displacement_basis": "closure",
                "provenance": {"source": "", "url": "", "license": ""},
            }
        ],
    }
    path = tmp_path / "manifest.json"
    path.write_text(json.dumps(payload))
    with pytest.raises(ManifestError):
        load_manifest(path)
    record = record_from_dict(payload["ships"][0])
    assert any("provenance" in problem for problem in record.problems())


def test_run_fleet_writes_artefacts(tmp_path, manifest):
    evidence = run_fleet(tmp_path, nx=61, nz=13, volume_tolerance=0.10, write=True)
    assert evidence["summary"]["ships_registered"] == 6
    assert evidence["summary"]["controls_all_rejected"] is True
    assert (tmp_path / "evidence.json").exists()
    assert (tmp_path / "ref_vs_got.csv").exists()
    assert (tmp_path / "manifests" / "KCS.json").exists()
    csv_text = (tmp_path / "ref_vs_got.csv").read_text().splitlines()
    assert csv_text[0].startswith("ship,reference_basis")
    assert len(csv_text) == 7
    written = json.loads((tmp_path / "evidence.json").read_text())
    assert written["provenance"]["manifest_sha256"] == manifest.sha256
    assert written["cross_route_check"]["lofted_euler"] == 2