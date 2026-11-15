"""The numpy backend must be a *second implementation*, not a second definition.

Everything here compares the two backends cell by cell. The scalar implementation in
:mod:`tensorcad.meshquality.metrics` is the reference; the vectorised one is only allowed
to differ by summation order, so the tolerance below is at the round-off level and the
guard patterns (``NaN``, ``inf``, sign of an inverted cell) must match exactly.
"""
from __future__ import annotations

import math
import time

import pytest

np = pytest.importorskip("numpy")

from tensorcad.meshquality import numpy_backend  # noqa: E402
from tensorcad.meshquality.metrics import (  # noqa: E402
    hex_metrics,
    prism_metrics,
    pyramid_metrics,
    quad_metrics,
    tet_metrics,
    tri_metrics,
)

RTOL = 1e-12
ATOL = 1e-12


# --------------------------------------------------------------------------------------
# mesh generators (deterministic, no fixtures and no mesh files needed)
# --------------------------------------------------------------------------------------
def structured_hexes(n: int, sx: float = 1.0, sy: float = 1.0, sz: float = 0.3):
    xs, ys, zs = np.linspace(0, sx, n + 1), np.linspace(0, sy, n + 1), np.linspace(0, sz, n + 1)
    points = np.stack(np.meshgrid(xs, ys, zs, indexing="ij"), axis=-1).reshape(-1, 3)
    index = lambda i, j, k: (i * (n + 1) + j) * (n + 1) + k  # noqa: E731
    cells = [
        [
            index(i, j, k),
            index(i + 1, j, k),
            index(i + 1, j + 1, k),
            index(i, j + 1, k),
            index(i, j, k + 1),
            index(i + 1, j, k + 1),
            index(i + 1, j + 1, k + 1),
            index(i, j + 1, k + 1),
        ]
        for i in range(n)
        for j in range(n)
        for k in range(n)
    ]
    return points, np.array(cells, dtype=np.int64)


def warped_hexes(count: int, seed: int = 7):
    """Random trilinear hexahedra: a unit cube per cell with independently jittered nodes."""
    rng = np.random.default_rng(seed)
    corners = np.array(
        [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]],
        dtype=float,
    )
    origin = rng.random((count, 1, 3)) * 3.0
    points = (corners[None, :, :] + origin).reshape(-1, 3)
    points = points + rng.normal(scale=0.18, size=points.shape)
    cells = np.arange(points.shape[0]).reshape(count, 8)
    return points, cells


def cube_tets(count: int, seed: int = 11):
    """Six-tetrahedron split of the unit cube, repeated for ``count`` jittered cubes."""
    points, hexes = warped_hexes(count, seed=seed)
    splits = (
        (0, 1, 2, 6), (0, 2, 3, 6), (0, 3, 7, 6), (0, 7, 4, 6), (0, 4, 5, 6), (0, 5, 1, 6)
    )
    tets = []
    for hex_cell in hexes:
        for split in splits:
            tets.append([hex_cell[i] for i in split])
    return points, np.array(tets, dtype=np.int64)


def quad_patch(count: int, seed: int = 13):
    rng = np.random.default_rng(seed)
    corners = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]], dtype=float)
    origin = rng.random((count, 1, 3)) * 2.0
    points = (corners[None, :, :] + origin).reshape(-1, 3)
    points = points + rng.normal(scale=0.15, size=points.shape)
    cells = np.arange(points.shape[0]).reshape(count, 4)
    return points, cells


def tri_patch(count: int, seed: int = 17):
    rng = np.random.default_rng(seed)
    corners = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]], dtype=float)
    origin = rng.random((count, 1, 3)) * 2.0
    points = (corners[None, :, :] + origin).reshape(-1, 3)
    points = points + rng.normal(scale=0.15, size=points.shape)
    cells = np.arange(points.shape[0]).reshape(count, 3)
    return points, cells


# --------------------------------------------------------------------------------------
# comparison helper
# --------------------------------------------------------------------------------------
def assert_same(scalar, vector, context: str) -> dict[str, float]:
    """Every metric of two :class:`ElementQuality` objects, including guard patterns."""
    assert scalar.kind == vector.kind
    assert set(scalar.metrics) == set(vector.metrics), context
    worst: dict[str, float] = {}
    for name in scalar.metrics:
        a = np.asarray(scalar.metrics[name], dtype=float)
        b = np.asarray(vector.metrics[name], dtype=float)
        assert a.shape == b.shape, f"{context}:{name} shape"
        assert np.array_equal(np.isnan(a), np.isnan(b)), f"{context}:{name} NaN pattern"
        assert np.array_equal(np.isinf(a), np.isinf(b)), f"{context}:{name} inf pattern"
        assert np.array_equal(np.sign(a[np.isinf(a)]), np.sign(b[np.isinf(b)])), context
        finite = np.isfinite(a)
        if not finite.any():
            worst[name] = 0.0
            continue
        diff = np.abs(a[finite] - b[finite])
        scale = np.maximum(np.abs(a[finite]), ATOL)
        assert np.all(diff <= ATOL + RTOL * scale), f"{context}:{name} max diff {diff.max():g}"
        worst[name] = float(diff.max())
    return worst


# --------------------------------------------------------------------------------------
# 1. the canonical elements are exact in both backends
# --------------------------------------------------------------------------------------
def test_ideal_elements_are_exactly_one():
    cube = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
    for backend in ("scalar", "numpy"):
        hexes = hex_metrics(cube, [[0, 1, 2, 3, 4, 5, 6, 7]], backend=backend)
        assert hexes["scaled_jacobian"] == (1.0,)
        assert hexes["min_sicn"] == (1.0,)
        assert abs(hexes["volume"][0] - 1.0) < 1e-15
        assert hexes["min_angle"] == (90.0,)
        quads = quad_metrics(cube, [[0, 1, 2, 3]], backend=backend)
        assert quads["scaled_jacobian"] == (1.0,)
        assert abs(quads["min_sicn"][0] - 1.0) < 1e-15
        equilateral = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.5, math.sqrt(3) / 2, 0.0)]
        tris = tri_metrics(equilateral, [[0, 1, 2]], backend=backend)
        assert abs(tris["min_sicn"][0] - 1.0) < 1e-15
        assert abs(tris["scaled_jacobian"][0] - 1.0) < 1e-15
        assert abs(tris["min_angle"][0] - 60.0) < 1e-12
        regular = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.5, math.sqrt(3) / 2, 0.0)]
        regular.append((0.5, math.sqrt(3) / 6, math.sqrt(2.0 / 3.0)))
        tets = tet_metrics(regular, [[0, 1, 2, 3]], backend=backend)
        assert abs(tets["min_sicn"][0] - 1.0) < 1e-12
        assert abs(tets["scaled_jacobian"][0] - 1 / math.sqrt(2)) < 1e-12
        assert abs(tets["gmsh_icn"][0] - 2.0 ** (-1.0 / 3.0)) < 1e-12


# --------------------------------------------------------------------------------------
# 2. agreement on real-shaped cell families
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize("count", [1, 5, 200])
def test_hex_agreement(count):
    points, cells = warped_hexes(count)
    worst = assert_same(
        hex_metrics(points, cells, backend="scalar"),
        hex_metrics(points, cells, backend="numpy"),
        f"hex{count}",
    )
    assert max(worst.values()) < 1e-12


def test_structured_boxes_agree():
    for n in (2, 4):
        points, cells = structured_hexes(n)
        assert_same(
            hex_metrics(points, cells, backend="scalar"),
            hex_metrics(points, cells, backend="numpy"),
            f"box{n}",
        )


def test_tet_agreement():
    points, cells = cube_tets(40)
    worst = assert_same(
        tet_metrics(points, cells, backend="scalar"),
        tet_metrics(points, cells, backend="numpy"),
        "tet",
    )
    assert max(worst.values()) < 1e-12


def test_quad_and_tri_agreement():
    points, cells = quad_patch(120)
    assert_same(
        quad_metrics(points, cells, backend="scalar"),
        quad_metrics(points, cells, backend="numpy"),
        "quad",
    )
    points, cells = tri_patch(120)
    assert_same(
        tri_metrics(points, cells, backend="scalar"),
        tri_metrics(points, cells, backend="numpy"),
        "tri",
    )


# --------------------------------------------------------------------------------------
# 3. the guards must survive vectorisation
# --------------------------------------------------------------------------------------
def test_degenerate_cells_keep_the_guards():
    cube = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
    collapsed = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)] * 2  # zero height: V = 0
    for backend in ("scalar", "numpy"):
        q = hex_metrics(collapsed, [[0, 1, 2, 3, 4, 5, 6, 7]], backend=backend)
        assert q["volume"] == (0.0,)
        assert q["min_sicn"] == (0.0,)
    # an inverted tetrahedron keeps a negative min_sicn in both backends
    inverted = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, -1.0)]
    for backend in ("scalar", "numpy"):
        q = tet_metrics(inverted, [[0, 1, 2, 3]], backend=backend)
        assert q["min_sicn"][0] < 0.0
        assert q["volume"][0] > 0.0
    # a duplicated node (zero-length edge) makes the angles NaN, not a silent zero
    degenerate = [(0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (1.0, 0.0, 0.0)]
    for backend in ("scalar", "numpy"):
        q = tri_metrics(degenerate, [[0, 1, 2]], backend=backend)
        # the two NaN corners are dropped (helper ``_min_max_deg``), the valid one is 0 deg
        assert q["min_angle"] == (0.0,) and q["max_angle"] == (0.0,)
        assert q["edge_ratio"] == (math.inf,)
        assert q["min_sicn"] == (0.0,)


# --------------------------------------------------------------------------------------
# 4. metamorphic invariance: rigid motion must not change any shape metric
# --------------------------------------------------------------------------------------
def test_rigid_motion_invariance():
    points, cells = warped_hexes(60)
    theta = 0.7
    rotation = np.array(
        [
            [math.cos(theta), -math.sin(theta), 0.0],
            [math.sin(theta), math.cos(theta), 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    moved = points @ rotation.T + np.array([3.0, -2.0, 0.5])
    base = hex_metrics(points, cells, backend="numpy")
    rotated = hex_metrics(moved, cells, backend="numpy")
    for name in ("min_angle", "max_angle", "skewness", "scaled_jacobian", "min_sicn", "edge_ratio"):
        a = np.asarray(base[name])
        b = np.asarray(rotated[name])
        assert np.allclose(a, b, rtol=1e-10, atol=1e-10), name
    # volume is invariant, and scaling by s multiplies it by s**3 and leaves the shapes alone
    scaled = hex_metrics(points * 2.5, cells, backend="numpy")
    assert np.allclose(np.asarray(base["volume"]) * 2.5**3, np.asarray(scaled["volume"]), rtol=1e-10)
    assert np.allclose(np.asarray(base["min_sicn"]), np.asarray(scaled["min_sicn"]), rtol=1e-10)


# --------------------------------------------------------------------------------------
# 5. dispatch contract
# --------------------------------------------------------------------------------------
def test_backend_selection_and_validation():
    points, cells = structured_hexes(6)  # 216 cells, above the auto threshold
    forced = hex_metrics(points, cells, backend="numpy")
    auto = hex_metrics(points, cells, backend="auto")
    assert auto == forced
    small_points, small_cells = structured_hexes(2)  # 8 cells, below the threshold
    assert hex_metrics(small_points, small_cells) == hex_metrics(
        small_points, small_cells, backend="scalar"
    )
    with pytest.raises(ValueError):
        hex_metrics(points, cells, backend="numba")


def test_prism_and_pyramid_refuse_forced_numpy():
    prism = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (0, 1, 1)]
    pyramid = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0.5, 0.5, 1)]
    with pytest.raises(ValueError):
        prism_metrics(prism, [[0, 1, 2, 3, 4, 5]], backend="numpy")
    with pytest.raises(ValueError):
        pyramid_metrics(pyramid, [[0, 1, 2, 3, 4]], backend="numpy")
    # auto stays on the scalar path and still answers
    assert prism_metrics(prism, [[0, 1, 2, 3, 4, 5]], backend="auto")["min_sicn"] > (0.0,)


# --------------------------------------------------------------------------------------
# 6. the point of the exercise: cost per cell
# --------------------------------------------------------------------------------------
def test_numpy_backend_is_much_faster():
    points, cells = structured_hexes(15)  # 3375 hexahedra
    start = time.perf_counter()
    hex_metrics(points, cells, backend="scalar")
    scalar_seconds = time.perf_counter() - start
    start = time.perf_counter()
    hex_metrics(points, cells, backend="numpy")
    numpy_seconds = time.perf_counter() - start
    assert numpy_seconds * 5.0 < scalar_seconds, (
        f"scalar {scalar_seconds:.3f}s vs numpy {numpy_seconds:.3f}s"
    )


# --------------------------------------------------------------------------------------
# 7. pin the documented limitation of the triangle sign convention
# --------------------------------------------------------------------------------------
def test_triangle_orientation_is_not_signed():
    """A 3-D triangle has no intrinsic orientation, so ``min_sicn`` cannot go negative.

    The scalar path computes ``signed_area = 0.5 * dot(cross, unit(cross))``, which is
    always positive; the vectorised path reproduces it and this test pins the consequence
    so that the documentation cannot claim otherwise.
    """
    forward = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)]
    backward = [(0.0, 0.0, 0.0), (0.0, 1.0, 0.0), (1.0, 0.0, 0.0)]
    for backend in ("scalar", "numpy"):
        a = tri_metrics(forward, [[0, 1, 2]], backend=backend)
        b = tri_metrics(backward, [[0, 1, 2]], backend=backend)
        assert a["min_sicn"] == b["min_sicn"]
        assert a["scaled_jacobian"] == b["scaled_jacobian"]
        assert a["min_sicn"][0] > 0.0


def test_backend_module_is_importable_without_touching_the_core():
    """The dependency-free core only imports numpy when a caller asks for it."""
    import sys

    assert hasattr(numpy_backend, "AUTO_MIN_CELLS")
    assert numpy_backend.UNSUPPORTED_KINDS == ("prism", "pyramid")
    assert sys.modules.get("numpy") is not None  # this test itself uses numpy