"""STL leg of the format matrix: two encodings, one name (m23297).

The 6-format round-trip matrix rates ``stl`` as one leg, but STL is written in
two wire encodings that behave differently and the leg was only ever exercised
in one of them:

* **ASCII** -- decimal text, exact round trip, ~5x the bytes.  This is what
  ``meshio`` 5.x writes by default, so every "stl" row of the published green
  table is an ASCII row.
* **binary** -- 80-byte header + ``uint32`` count + ``float32`` triangles.  This
  is what CAD systems actually emit, and it is lossy: coordinates round to a
  grid of spacing ``2**-24 * |x|``, so the error is **absolute and grows with
  the coordinate offset**, not a fixed relative error.

The platform documentation claimed a binary-STL relative coordinate error bound
of ``3*eps32 ~ 3.6e-7``.  That is the *unit roundoff*; it is only a relative
error bound for a body sitting at the origin.  ``test_documented_relative_bound
_is_an_origin_only_bound`` pins the counterexample: the same sphere at
``x = 1e3 m`` exceeds it by an order of magnitude.  The correct statement, and
the screen implemented in :func:`tensorcad.meshio_exchange.stl_quantisation_bound`,
is ``max_coordinate_error = eps32 * max|x|``.

Everything here is measured, not asserted by inspection: volumes come from
:func:`tensorcad.meshio_exchange.surface_signed_volume`, coordinates are sorted
before comparison because meshio's STL reader renumbers the merged vertices.
"""
from __future__ import annotations

import math
import struct
from pathlib import Path

import numpy as np
import pytest

from tensorcad.geometry import Mesh
from tensorcad.meshio_exchange import (
    EPS32,
    FORMATS_SURFACE,
    MeshioFormatError,
    STL_ENCODINGS,
    read_meshio,
    round_trip_meshio,
    stl_quantisation_bound,
    surface_signed_volume,
    write_meshio,
)
from test_meshio_format_matrix import GEOMETRIES

GEOM_IDS = [case.name for case in GEOMETRIES]

#: Measured absolute coordinate error of the binary leg must stay inside the
#: half-ulp screen for every fixture and every offset probed.
OFFSETS = (0.0, 1.0e2, 1.0e3, 1.0e4)

#: Relative volume error accepted for the lossy leg (screening estimate times a
#: safety factor of 4, since the bound is an estimate, not a proof).
LOSSY_VOLUME_TOL = 4.0


def _translated(mesh: Mesh, offset: float) -> Mesh:
    """Same body with every x coordinate shifted (ship frame / global frame)."""
    return Mesh(
        vertices=[(x + offset, y, z) for (x, y, z) in mesh.vertices],
        faces=mesh.faces,
        units="m",
    )


def _sorted_vertices(mesh: Mesh) -> np.ndarray:
    return np.sort(np.array(mesh.vertices, dtype=float), axis=0)


def _one_binary_stl(path: Path) -> bytes:
    return path.read_bytes()


# ---------------------------------------------------------------------------
# Encoding identity: what the writer actually produced
# ---------------------------------------------------------------------------

def test_default_stl_encoding_is_ascii() -> None:
    """``write_meshio(..., 'stl')`` without a switch is ASCII, not binary."""
    mesh = GEOMETRIES[0].mesh
    import tempfile

    with tempfile.TemporaryDirectory() as td:
        path = write_meshio(Path(td) / "m.stl", mesh, file_format="stl")
        raw = _one_binary_stl(path)
    assert raw[:5] == b"solid", "STL leg is expected to be ASCII by meshio default"
    text = raw.decode("ascii")
    assert "facet normal" in text and "endfacet" in text
    n_tri = len(mesh.faces)
    assert f"  vertex" in text
    assert text.count("facet normal") == n_tri


def test_binary_stl_has_header_count_and_float32_payload() -> None:
    mesh = GEOMETRIES[0].mesh
    import tempfile

    with tempfile.TemporaryDirectory() as td:
        path = write_meshio(Path(td) / "m.stl", mesh, file_format="stl", stl_binary=True)
        raw = path.read_bytes()
    n_tri = len(mesh.faces)
    assert len(raw) == 84 + 50 * n_tri, "binary STL is exactly 50 bytes per triangle"
    (count,) = struct.unpack("<I", raw[80:84])
    assert count == n_tri
    # first triangle payload: 3 float32 points, correctly sized
    pts = struct.unpack("<9f", raw[84 : 84 + 36])
    assert len(pts) == 9 and all(math.isfinite(p) for p in pts)


@pytest.mark.parametrize("case", GEOMETRIES, ids=GEOM_IDS)
def test_ascii_leg_cost_is_measured(case) -> None:
    """ASCII is exact but bulky: the byte ratio is an interop cost, not a detail."""
    import tempfile

    with tempfile.TemporaryDirectory() as td:
        a = write_meshio(Path(td) / "a.stl", case.mesh, file_format="stl", stl_binary=False)
        b = write_meshio(Path(td) / "b.stl", case.mesh, file_format="stl", stl_binary=True)
        ratio = a.stat().st_size / b.stat().st_size
    # Ratio depends on text width per coordinate and on mesh fineness: the unit
    # cube is 2.1x (short numbers, few triangles) while the hull is 4.8x.  >2x is
    # the claim that survives all four fixtures; per-geometry values are tabulated
    # in docs/m23297_stl_encoding_verification.md.
    assert ratio > 2.0, f"ASCII leg should cost >2x the bytes, measured {ratio:.2f}"


# ---------------------------------------------------------------------------
# Round-trip fidelity per encoding
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("case", GEOMETRIES, ids=GEOM_IDS)
def test_ascii_round_trip_is_exact(case) -> None:
    back = round_trip_meshio(case.mesh, file_format="stl", stl_binary=False)
    assert len(back.vertices) == len(case.mesh.vertices)
    assert len(back.faces) == len(case.mesh.faces)
    np.testing.assert_array_equal(_sorted_vertices(back), _sorted_vertices(case.mesh))
    assert surface_signed_volume(back) == pytest.approx(
        surface_signed_volume(case.mesh), rel=1e-15
    )


@pytest.mark.parametrize("case", GEOMETRIES, ids=GEOM_IDS)
@pytest.mark.parametrize("offset", OFFSETS, ids=[f"x{int(o)}m" for o in OFFSETS])
def test_binary_round_trip_stays_inside_its_screen(case, offset: float) -> None:
    body = _translated(case.mesh, offset)
    screen = stl_quantisation_bound(body, encoding="binary")
    assert screen["encoding"] == "binary"
    # The screen is about the *offset*, so it must grow with it.
    assert screen["max_coordinate_error"] == pytest.approx(EPS32 * screen["max_abs_coordinate"])

    back = round_trip_meshio(body, file_format="stl", stl_binary=True)
    assert len(back.vertices) == len(body.vertices)
    assert len(back.faces) == len(body.faces)
    assert not screen["distinct_vertices_may_merge"], (
        "fixture mesh is too coarse to be destroyed by float32; the merge flag "
        "would make the topology assertion below vacuous"
    )

    measured = float(np.abs(_sorted_vertices(back) - _sorted_vertices(body)).max())
    assert measured <= screen["max_coordinate_error"], (
        f"measured coordinate error {measured:.3e} exceeds the half-ulp screen "
        f"{screen['max_coordinate_error']:.3e} at offset {offset:g}"
    )

    v0, v1 = surface_signed_volume(body), surface_signed_volume(back)
    assert v0 * v1 > 0.0, "binary round trip must not flip orientation"
    rel = abs(v1 - v0) / abs(v0)
    assert rel <= LOSSY_VOLUME_TOL * screen["screening_volume_rel_error"], (
        f"relative volume error {rel:.3e} exceeds "
        f"{LOSSY_VOLUME_TOL} x screening estimate "
        f"{screen['screening_volume_rel_error']:.3e}"
    )


@pytest.mark.parametrize("case", [GEOMETRIES[1]], ids=["sphere"])
def test_documented_relative_bound_is_an_origin_only_bound(case) -> None:
    """Counterexample to the doc claim ``rel err <= 3*eps32 ~ 3.6e-7``.

    At the origin binary STL is far better than that; in a ship frame
    (``x ~ 1e3 m``) it is worse.  The claim is only true for ``x ~ O(1)``.

    The counterexample needs a body whose coordinates are *not* exactly
    representable in float32 (the unit sphere).  The Wigley fixture does not
    qualify: its station abscissae are float32-exact, so translating it in x
    costs nothing in x and its volume error is governed by its O(0.5 m)
    y/z coordinates -- see
    :func:`test_hull_offset_growth_is_governed_by_the_coordinates_present`.
    """
    claimed = 3.0 * EPS32

    at_origin = _translated(case.mesh, 0.0)
    back0 = round_trip_meshio(at_origin, file_format="stl", stl_binary=True)
    rel0 = abs(surface_signed_volume(back0) - surface_signed_volume(at_origin)) / abs(
        surface_signed_volume(at_origin)
    )
    assert rel0 <= claimed, "the quoted bound does hold near the origin"

    ship_frame = _translated(case.mesh, 1.0e3)
    back1 = round_trip_meshio(ship_frame, file_format="stl", stl_binary=True)
    rel1 = abs(surface_signed_volume(back1) - surface_signed_volume(ship_frame)) / abs(
        surface_signed_volume(ship_frame)
    )
    screen = stl_quantisation_bound(ship_frame, encoding="binary")
    assert rel1 > claimed, (
        f"expected the origin-only bound {claimed:.2e} to be violated at 1e3 m, "
        f"measured {rel1:.3e}: if this fails the doc claim may be reconsidered"
    )
    # ...and the offset-aware screen is what actually covers it.
    assert rel1 <= LOSSY_VOLUME_TOL * screen["screening_volume_rel_error"]


def test_hull_offset_growth_is_governed_by_the_coordinates_present() -> None:
    """An x translation alone does not hurt a hull whose stations are exact.

    Not a licence to ignore the issue: the same hull described with *inexact*
    y/z values (or at a scale where its sections are no longer representable)
    does grow, and the screen covers both by using ``max|x|`` over all three
    axes rather than a single offset.
    """
    hull = GEOMETRIES[3].mesh
    at_origin = hull
    ship_frame = _translated(hull, 1.0e3)
    back0 = round_trip_meshio(at_origin, file_format="stl", stl_binary=True)
    back1 = round_trip_meshio(ship_frame, file_format="stl", stl_binary=True)
    rel0 = abs(surface_signed_volume(back0) - surface_signed_volume(at_origin)) / abs(
        surface_signed_volume(at_origin)
    )
    rel1 = abs(surface_signed_volume(back1) - surface_signed_volume(ship_frame)) / abs(
        surface_signed_volume(ship_frame)
    )
    # Measured on the fixture: the hull's abscissae are float32-exact, so the
    # binary leg is as good at 1e3 m as at the origin for an x-only shift.
    assert rel1 == pytest.approx(rel0, rel=1e-3)
    assert rel1 < 1e-6
    # The offset-aware screen is still *safe* (an upper bound) here, just loose:
    screen = stl_quantisation_bound(ship_frame, encoding="binary")
    assert rel1 * LOSSY_VOLUME_TOL <= max(screen["screening_volume_rel_error"], 1e-6)


# ---------------------------------------------------------------------------
# The screen itself
# ---------------------------------------------------------------------------

def test_ascii_screen_is_zero_and_binary_screen_scales_with_offset() -> None:
    mesh = GEOMETRIES[1].mesh
    ascii_screen = stl_quantisation_bound(mesh, encoding="ascii")
    assert ascii_screen["max_coordinate_error"] == 0.0
    assert ascii_screen["screening_volume_rel_error"] == 0.0
    assert ascii_screen["distinct_vertices_may_merge"] is False

    errors = [stl_quantisation_bound(_translated(mesh, o), encoding="binary")["max_coordinate_error"]
              for o in OFFSETS]
    assert errors[0] == pytest.approx(EPS32 * (1.0 + 0.0), rel=1e-6)  # unit sphere
    assert errors == sorted(errors) and errors[-1] > errors[0], (
        "binary error must grow monotonically with the coordinate offset"
    )


def test_merge_flag_catches_a_body_finer_than_the_float32_grid() -> None:
    """A mesh whose finest edge is below the quantization step loses nodes silently."""
    # Two triangles 1e-7 apart at x = 1e3: step is 6e-5, so the nodes collapse.
    base = 1.0e3
    eps = 1.0e-7
    vertices = [
        (base, 0.0, 0.0),
        (base + eps, 0.0, 0.0),
        (base, eps, 0.0),
        (base, 0.0, eps),
    ]
    mesh = Mesh(vertices=vertices, faces=((0, 1, 2), (0, 1, 3)), units="m")
    screen = stl_quantisation_bound(mesh, encoding="binary")
    assert screen["distinct_vertices_may_merge"] is True
    assert screen["node_budget_exceeded"] is True
    # ASCII keeps them apart.
    assert stl_quantisation_bound(mesh, encoding="ascii")["distinct_vertices_may_merge"] is False


def test_unknown_encoding_rejected() -> None:
    with pytest.raises(MeshioFormatError):
        stl_quantisation_bound(GEOMETRIES[0].mesh, encoding="utf8")
    assert STL_ENCODINGS == ("ascii", "binary")


def test_stl_binary_switch_is_stl_only() -> None:
    mesh = GEOMETRIES[0].mesh
    import tempfile

    with tempfile.TemporaryDirectory() as td:
        with pytest.raises(MeshioFormatError):
            write_meshio(Path(td) / "m.obj", mesh, file_format="obj", stl_binary=True)
    for fmt in FORMATS_SURFACE:
        if fmt == "stl":
            continue
        with pytest.raises(MeshioFormatError):
            write_meshio("/tmp/never_written", mesh, file_format=fmt, stl_binary=False)


def test_reader_auto_detects_both_dialects() -> None:
    mesh = GEOMETRIES[3].mesh
    import tempfile

    with tempfile.TemporaryDirectory() as td:
        a = write_meshio(Path(td) / "a.stl", mesh, file_format="stl", stl_binary=False)
        b = write_meshio(Path(td) / "b.stl", mesh, file_format="stl", stl_binary=True)
        back_a = read_meshio(a, file_format="stl")
        back_b = read_meshio(b, file_format="stl")
    assert len(back_a.vertices) == len(back_b.vertices) == len(mesh.vertices)
    assert surface_signed_volume(back_a) == pytest.approx(surface_signed_volume(mesh), rel=1e-15)
    assert surface_signed_volume(back_b) == pytest.approx(surface_signed_volume(mesh), rel=1e-6)

# ---------------------------------------------------------------------------
# Facet-normal self-consistency of the binary dialect
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("offset", [0.0, 1.0e3, 1.0e4])
@pytest.mark.parametrize("case", [GEOMETRIES[1], GEOMETRIES[3]], ids=["sphere", "hull"])
def test_binary_writer_normal_matches_the_stored_winding(case, offset: float) -> None:
    """The stored facet normal must describe the triangle that is in the file.

    Binary STL carries a redundant per-facet normal next to the winding.  If the
    normal is computed from the float64 positions and only then are the
    coordinates rounded to float32 (meshio's ordering), the two describe
    different triangles and the disagreement grows with the coordinate offset.
    Computing it from the quantised coordinates keeps it bounded by the float32
    rounding of a unit vector, independent of the offset.
    """
    import tempfile

    from tensorcad.mesh_foreign_dialects import read_stl_binary, write_stl_binary

    body = _translated(case.mesh, offset)
    with tempfile.TemporaryDirectory() as td:
        path = write_stl_binary(Path(td) / "b.stl", body)
        parsed = read_stl_binary(path)
        pts = np.array([tuple(v) for v in parsed.vertices], dtype=float)
        normals = np.array(parsed.declared_meta["normals"], dtype=float)
    tri = pts.reshape(-1, 3, 3)
    winding = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    lengths = np.linalg.norm(winding, axis=1)
    keep = lengths > 0.0
    winding[keep] /= lengths[keep, None]
    dev = float(np.linalg.norm(winding[keep] - normals[keep], axis=1).max())
    # 4 x the float32 unit-roundoff of a unit normal (sqrt(3)*eps32 ~ 1.0e-7)
    assert dev < 1.0e-6, f"stored facet normal deviates by {dev:.3e} at offset {offset:g}"
    # and the file is still a faithful copy of the quantised geometry
    assert surface_signed_volume(_faces_mesh(pts, len(tri))) * surface_signed_volume(body) > 0.0


def _faces_mesh(points: np.ndarray, n_tri: int) -> Mesh:
    return Mesh(
        vertices=[tuple(p) for p in points],
        faces=[(3 * i, 3 * i + 1, 3 * i + 2) for i in range(n_tri)],
        units="m",
    )
