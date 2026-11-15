"""Contract tests for the official-geometry route of the baseline fleet.

The route's whole point is that a verdict is decidable from the artefact: a digest that
does not match is refused, a hull whose patches are trimmed is *named* as such instead of
being integrated as if it were whole, and an entry with no official file says so.  The
hull builds are skipped when the publisher bytes are not in the local cache, because the
bytes are deliberately not redistributed -- a skipped build is reported, never faked.
"""

from __future__ import annotations

import dataclasses

import pytest

from tensorcad.fleet import official_geometry as og
from tensorcad.geometry import Mesh


def test_series60_has_no_official_geometry():
    """A methodical series has offset tables, not a hull file: the verdict must say so."""
    record = og.verify_ship(og.SPECS["Series60"])
    assert record["verdict"] == "NO_OFFICIAL_GEOMETRY"
    assert "no official hull file exists" in record["detail"]


def test_euler_characteristic_on_closed_solids():
    """The topology check the gate uses must return 2 for a closed genus-0 surface."""
    tetra = Mesh(
        vertices=[(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)],
        faces=[(0, 2, 1), (0, 1, 3), (0, 3, 2), (1, 2, 3)],
    )
    assert og._euler_characteristic(tetra) == 2


def test_declared_digest_mismatch_is_refused():
    """A source change must not pass silently: a wrong pin is a named verdict."""
    ship = og.SPECS["KCS"]
    if og.materialise(ship) is None:
        pytest.skip("KCS publisher bytes not in the local cache")
    mis_pinned = dataclasses.replace(ship, iges_sha256="0" * 64)
    record = og.verify_ship(mis_pinned)
    assert record["verdict"] == "DIGEST_MISMATCH"
    assert record["provenance"]["geometry_sha256_measured"] != "0" * 64


def test_materialise_strips_the_cache_placeholder():
    ship = og.SPECS["KCS"]
    path = og.materialise(ship)
    if path is None:
        pytest.skip("KCS publisher bytes not in the local cache")
    assert "{cache}" not in str(path)
    assert path.is_file()


@pytest.mark.parametrize("name", ["KCS"])
def test_kcs_official_iges_is_verified(name):
    """The accepted hull: watertight, genus 0, volume and Cb inside the 1 % anchor band."""
    ship = og.SPECS[name]
    if og.materialise(ship) is None:
        pytest.skip(f"{name} publisher bytes not in the local cache")
    record = og.verify_ship(ship)
    assert record["verdict"] == "VERIFIED", record.get("detail")
    assert record["mesh"]["is_watertight"] is True
    assert record["mesh"]["is_edge_manifold"] is True
    assert record["mesh"]["boundary_edges"] == 0
    assert record["mesh"]["euler_characteristic"] == 2
    assert abs(record["volume_vs_published"]["rel_error_percent"]) < 100 * og.VOLUME_TOLERANCE
    assert abs(record["hydrostatics"]["cb_rel_error_percent"]) < 100 * og.CB_TOLERANCE
    # the verdict must carry provenance, not just a number
    assert record["provenance"]["geometry_sha256_measured"] == ship.iges_sha256


@pytest.mark.parametrize("name", ["JBC", "KVLCC2"])
def test_trim_or_scale_blocked_hulls_are_not_reported_as_verified(name):
    """Honest-FAIL guard: a hull we cannot build yet must never come back VERIFIED."""
    ship = og.SPECS[name]
    if og.materialise(ship) is None:
        pytest.skip(f"{name} publisher bytes not in the local cache")
    record = og.verify_ship(ship)
    assert record["verdict"] in {"WITHHELD", "TRIM_UNHANDLED", "FAILED", "DIGEST_MISMATCH"}
    assert record["verdict"] != "VERIFIED"
    # a refusal must name what was measured, not merely fail
    assert record.get("detail")
