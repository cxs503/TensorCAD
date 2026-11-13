"""Differentiable DTMB 5415 parametric reconstruction from the 9001 x 27 design database.

The Zenodo design-space benchmark dataset (Serani, Palma, Wackers, Quagliarella, Gaggero,
Diez; see :data:`DATASET` for the provenance actually used here) publishes **9001
configurations of the DTMB 5415 hull**, each as 2250 surface points plus the **27
normalised design variables** that generated it.  What it does *not* publish is the
parametric model itself: the map ``design variables -> geometry`` has to be rebuilt from
the data.  This module is that rebuild, in the form TensorCAD needs: a callable whose
input is the 27-vector ``v`` and whose outputs are the point cloud, a triangulated surface
mesh, hydrostatics, and the **analytic derivative of each with respect to ``v``**.

How the reconstruction is parameterised
--------------------------------------

The 2250 points of one configuration are not an arbitrary cloud: they are a
``2 patches x 25 girth x 45 station`` structured grid (aft patch ``x in [-0.10, 80.13]``,
fore patch ``x in [81.79, 151.80]`` at full scale, interleaved in the published flattening
order ``index = 45 * (2 * j + patch) + i``).  The layout is *measured* off the baseline
column by :func:`point_grid`, which re-checks the invariants the layout implies, so a
re-ordered dataset cannot silently produce a scrambled surface.

The design space is then a **modal displacement field** over that grid::

    d(v)  = Phi . (C1 v + C2 (v (x) v))          # lagrangian displacement of the base grid
    P(v)  = P_ref + d(v)                          # 2250 points; v = 0 -> the reference hull

``Phi`` (6750 x 40) is the displacement PCA basis of the fit split, ``C1``/``C2`` the
least-squares coefficients on the 27 variables and their 378 distinct quadratic terms.  The
quadratic lift is there because the dataset README states the generating model is nonlinear;
``runs/m17352_5415_param/design_space_metrics.json`` records what the two families actually
achieve on a held-out split (linear 3.95 % of beam RMS, quadratic 2.45 %), which is the
honest statement of what a data-recovered surrogate can and cannot do.

Mesh
----

The surface mesh is a bicubic (natural) spline interpolation of the two structured patches in
their ``(girth, station)`` index space, refined to ``49 x 133`` per patch, mirrored to the
port side about ``y = 0`` and stitched across the patch seam.  The refinement and the
mirroring are **fixed linear operators** (the knots sit on the baseline index grid and never
move), so::

    vertices(v) = L . P(v) = L . P_ref + (L . Phi) . (C1 v + C2 vv)

is affine in the modal amplitudes and its Jacobian is exact.  The geometry is differentiable
by construction, not by finite-differencing a CAD kernel.

What the derivative is, and is not
----------------------------------

* :meth:`DTMB5415Space.vertex_jacobian` -- exact ``d vertex / d v``.
* :meth:`DTMB5415Space.volume` / :meth:`DTMB5415Space.volume_gradient` -- the enclosed volume
  of the hull shell, integrated as ``2 * integral A(x) dx`` over the half-section areas closed
  on the centreline; the gradient is analytic and matches a central difference of the same
  functional to round-off.
  The functional is deliberately **the shell volume to the deck edge, not the displacement
  below the waterline**: the waterline clip is piecewise smooth (its derivative jumps at
  every section that crosses the plane) and claiming a derivative through it would be a lie.
  Displacement below the waterline is available, non-differentiated, from
  :meth:`DTMB5415Space.hydrostatics`.
* The chain is therefore intact through geometry -> mesh -> integral, and bypassed only at the
  waterline clip; that gap is reported rather than hidden.

Data layout
-----------

``docs/data/dtmb5415/design_space.npz`` carries ``phi``, ``c1``, ``c2`` and ``base_cloud``
(~1.2 MB, float32 basis); ``design_space_metrics.json`` next to it carries the fit/holdout
residuals, the split rule and the sha256 of the 1.0 GB source ``.mat``.  Nothing here is
hand-typed: every number comes from ``tools/dtmb5415_param/fit_design_space.py``.

NumPy is imported inside the functions, never at module import, so ``import tensorcad`` still
works with the standard library alone (same rule as :mod:`tensorcad.voxel`).
"""
from __future__ import annotations

import json
import math
import os
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .geometry import Mesh

__all__ = [
    "DATASET",
    "GRID",
    "MESH_RESOLUTION",
    "REFERENCE",
    "SPACE_SCHEMA",
    "WATERLINE_Z",
    "DTMB5415Space",
    "point_grid",
    "spline_operator",
]

#: Schema tag for :meth:`DTMB5415Space.to_dict`; bump on incompatible change.
SPACE_SCHEMA = "tensorcad.dtmb5415.space/1"

#: Provenance of the design database this module was fitted from.  ``source_sha256`` is of
#: the 1.0 GB ``DTMB5415-database.mat`` actually read; matrix layout and the reference table
#: are quoted from the dataset's own README/metadata shipped inside the record.
DATASET: dict[str, Any] = {
    "title": "Design-Space Dimensionality Reduction Benchmark Dataset - DTMB 5415 Hull",
    "publisher": "Zenodo",
    "authors": ["Serani, A.", "Palma, G.", "Wackers, J.", "Quagliarella, D.",
                "Gaggero, S.", "Diez, M."],
    "reference": ("Serani et al. (2025) Extending parametric model embedding with physical "
                  "information for design-space dimensionality reduction in shape optimization, "
                  "Engineering with Computers"),
    "n_configurations": 9001,
    "n_design_variables": 27,
    "n_points_per_configuration": 2250,
    "geometry_scale": "full scale, Lpp = 142 m",
    "variables_matrix": "variablesMatrix (27 x 9001), normalised design variables in [-1, 1]",
    "geometry_matrices": "xMatrix/yMatrix/zMatrix (2250 x 9001), flattened hull coordinates",
    "column_0": "baseline (DTMB 5415 reference hull)",
    "readme_model_scale_table": {
        "Lpp_m": 5.720,
        "beam_m": 0.760,
        "draft_m": 0.248,
        "displacement_t": 0.549,
        "rho_kg_m3": 998.5,
        "LCG_m": 2.884,
        "VCG_m": 0.056,
        "note": "quoted from the dataset README; the geometry itself is full scale, Lpp = 142 m",
    },
    "source_file": "/root/ship-data-official/dtmb5415/DTMB5415-database.mat",
    "source_size_bytes": 1086721690,
    "source_sha256": "16b14ba598787eeed74eb6a10e9ea6fb701bd25496b5ccbf8a317ad6e364bee2",
}

#: The dense reference hull the config-0 acceptance check is measured against: the
#: platform's ``DTMB5415_Hull_exact.stl``, itself built from the dataset's baseline cloud
#: (same bounding box: x[-0.1038, 151.7966], y[0, 10.2695], z[-3.0245, 16.1708]).
REFERENCE: dict[str, Any] = {
    "path": ("/root/ship-performance-platform-incoming/ship-performance-platform/backend/"
             "data/geometry/ships/DTMB5415_Hull_exact.stl"),
    "sha256": "d40b17b1fc963e0d31ad8a30d2fe335758cb1e3c0324c968eafaf5fb9a3977fd",
    "n_triangles": 164926,
    "n_vertices": 83350,
    "units": "m, full scale",
}

#: Structured layout of one published configuration (measured, see :func:`point_grid`).
GRID = {"n_girth": 25, "n_station": 45, "n_patches": 2, "patch_names": ("fore", "aft")}

#: Refinement of the parametric surface used for the mesh (per patch).
MESH_RESOLUTION = {"n_station_out": 133, "n_girth_out": 49}

#: Scale factor of the shipped database: the record's README table is at model scale
#: (Lpp = 5.720 m), the geometry itself is full scale (Lpp = 142 m).
SCALE_LAMBDA = 142.0 / 5.720

#: Design draft used by :meth:`DTMB5415Space.hydrostatics`, measured from the baseline
#: (z = 0; the bulb hangs below, z_keel = -3.28 m).  Taken straight from the README table
#: scaled to full size: ``0.248 m * 142 / 5.720 = 6.1566 m``.  It is deliberately *not*
#: back-solved from a target displacement -- that made round 0 circular (the waterline was
#: tuned to a volume, then the volume was "validated" against it) and carried a t / (kg/m^3)
#: unit slip that made the target 1000x too small (see git history of this constant).
#: Geometry validation against the README table: ``tools/dtmb5415_param/verify_geometry.py``.
DRAFT_M = 0.248 * SCALE_LAMBDA
WATERLINE_Z = DRAFT_M

_DATA_DIR = Path(__file__).resolve().parent.parent / "docs" / "data" / "dtmb5415"


# ---------------------------------------------------------------------------
# structured grid of one published configuration
# ---------------------------------------------------------------------------
def point_grid(cloud: Any, check: bool = True) -> Any:
    """Reshape one 2250-point configuration into ``(2, n_girth, n_station, 3)``.

    The published flattening is ``index = 45 * (2 * j + patch) + i`` with ``patch = 0`` the
    fore body (bow -> midship, x descending) and ``patch = 1`` the aft body (midship ->
    transom, x descending).  The checks are exactly the invariants that layout implies, so a
    scrambled ordering raises instead of producing a silently broken surface:

    * every station column is near-planar -- the x spread inside a column is smaller than the
      longitudinal spacing between columns;
    * the girth runs from the deck edge (large ``y``) down to the centreline (``y == 0``).
    """
    import numpy as np

    cloud = np.asarray(cloud, dtype=np.float64)
    n_expected = GRID["n_patches"] * GRID["n_girth"] * GRID["n_station"]
    if cloud.shape != (n_expected, 3):
        raise ValueError(f"expected {n_expected} points, got {cloud.shape}")
    blk = cloud.reshape(2 * GRID["n_girth"], GRID["n_station"], 3)
    grid = np.stack([blk[0::2], blk[1::2]], axis=0)            # (2, girth, station, 3)
    if not check:
        return grid
    x_spread = float(np.abs(np.diff(grid[:, :, :, 0], axis=1)).max())
    x_extent = float(grid[:, :, :, 0].max() - grid[:, :, :, 0].min())
    if not x_spread < 0.05 * x_extent:
        raise ValueError(
            "station columns are not near-planar: within-column x spread "
            f"{x_spread:.3f} m >= 5% of the hull length {x_extent:.3f} m"
        )
    if not float(np.abs(grid[:, -1, :, 1]).max()) < 1e-6:
        raise ValueError("girth does not end on the centreline (y != 0 at the last girth row)")
    deck_wider = float(np.mean(grid[:, 0, :, 1] > grid[:, -1, :, 1] + 1e-9))
    # the stem station collapses to a point, so allow a small fraction of stations where the
    # deck edge is not wider than the centreline; a scrambled girth would fail almost all
    if not deck_wider > 0.9:
        raise ValueError(
            f"girth does not run deck edge -> centreline (only {100 * deck_wider:.1f}% of stations)")
    return grid


# ---------------------------------------------------------------------------
# fixed linear operator: natural cubic spline interpolation on the index grid
# ---------------------------------------------------------------------------
def spline_operator(n_in: int, n_out: int) -> Any:
    """Dense ``(n_out, n_in)`` matrix of natural-cubic interpolation on unit-spaced knots.

    The knots are the integer indices ``0 .. n_in - 1`` and never move, which is what makes
    the refinement a *fixed* linear operator: the constancy of the matrix across design
    variables is the differentiability guarantee, not an accident.
    """
    import numpy as np

    if n_in < 3 or n_out < 2:
        raise ValueError("need n_in >= 3 and n_out >= 2")
    h = 1.0
    n = n_in
    a = np.zeros((n, n))
    d2 = np.zeros((n, n))
    for i in range(1, n - 1):
        a[i, i - 1] = h
        a[i, i] = 4.0 * h
        a[i, i + 1] = h
        d2[i, i - 1] = 1.0 / h
        d2[i, i] = -2.0 / h
        d2[i, i + 1] = 1.0 / h
    a[0, 0] = 1.0
    a[-1, -1] = 1.0
    # natural spline: M solves  a M = 6 * (second difference of the data);  solving for all n
    # unit data vectors at once gives M[:, j] = d2 f / dx2 produced by a unit datum at node j
    m = np.linalg.solve(a, 6.0 * d2)
    op = np.zeros((n_out, n))
    for k, x in enumerate(np.linspace(0.0, n - 1, n_out)):
        i = min(int(math.floor(x)), n - 2)
        t = x - i
        h00 = 2 * t ** 3 - 3 * t ** 2 + 1
        h10 = t ** 3 - 2 * t ** 2 + t
        h01 = -2 * t ** 3 + 3 * t ** 2
        h11 = t ** 3 - t ** 2
        e_i = np.zeros(n)
        e_i[i] = 1.0
        e_j = np.zeros(n)
        e_j[i + 1] = 1.0
        op[k] = h00 * e_i + h10 * h * m[i] + h01 * e_j + h11 * h * m[i + 1]
    return op


@dataclass(frozen=True)
class _Operators:
    """Cached spline operators for the two index directions."""

    along_station: Any
    along_girth: Any


_OPS: _Operators | None = None


def _operators() -> _Operators:
    global _OPS
    if _OPS is None:
        _OPS = _Operators(
            along_station=spline_operator(GRID["n_station"], MESH_RESOLUTION["n_station_out"]),
            along_girth=spline_operator(GRID["n_girth"], MESH_RESOLUTION["n_girth_out"]),
        )
    return _OPS


def _refine(grid: Any) -> Any:
    """Bicubic interpolation of a ``(girth, station, ...)`` array onto the mesh grid."""
    import numpy as np

    ops = _operators()
    out = np.einsum("si,jic->jsc", ops.along_station, grid)
    return np.einsum("gj,jsc->gsc", ops.along_girth, out)


def _connectivity(n_girth: int, n_station: int) -> tuple[Any, Any, Any]:
    """Fixed triangle list for the starboard shell (two patches + seam bridge)."""
    import numpy as np

    fore = np.arange(n_girth * n_station).reshape(n_girth, n_station)
    aft = np.arange(n_girth * n_station, 2 * n_girth * n_station).reshape(n_girth, n_station)
    quads: list[tuple[int, int, int, int]] = []
    for patch in (fore, aft):
        for j in range(n_girth - 1):
            for i in range(n_station - 1):
                quads.append((patch[j, i], patch[j, i + 1],
                              patch[j + 1, i + 1], patch[j + 1, i]))
    # bridge across the fore/aft seam: the fore patch runs bow -> midship with descending x,
    # so the seam is its *last* station column and the aft patch's *first*.
    for j in range(n_girth - 1):
        quads.append((fore[j, -1], fore[j + 1, -1], aft[j + 1, 0], aft[j, 0]))
    tris = []
    for a, b, c, d in quads:
        tris.append((a, b, c))
        tris.append((a, c, d))
    return np.asarray(tris, dtype=np.int64), fore, aft

class DTMB5415Space:
    """The 27-variable DTMB 5415 design space: ``v -> cloud / mesh / hydrostatics / dV/dv``.

    Parameters
    ----------
    data_path :
        Directory holding ``design_space.npz`` (and ``design_space_metrics.json``); defaults
        to the copy shipped in ``docs/data/dtmb5415``.
    """

    def __init__(self, data_path: str | os.PathLike | None = None) -> None:
        import numpy as np

        directory = Path(data_path) if data_path is not None else _DATA_DIR
        npz = np.load(directory / "design_space.npz")
        self.phi = np.asarray(npz["phi"], dtype=np.float64)                 # (6750, n_modes)
        # the fit stores the coefficient matrices as (feature, mode); transpose once here
        self.c1 = np.asarray(npz["c1"], dtype=np.float64).T                  # (n_modes, 27)
        self.c2 = np.asarray(npz["c2"], dtype=np.float64).T                  # (n_modes, 405)
        self.base_cloud = np.asarray(npz["base_cloud"], dtype=np.float64)   # (2250, 3)
        # ``beam_m`` in the artifact is the HALF breadth of the base cloud (max |y|);
        # the moulded beam at the design waterline is 2x that (see ``beam_full_m``).
        self.beam_m = float(npz["beam_m"])
        self.beam_full_m = 2.0 * self.beam_m
        self.n_modes, self.n_variables = self.c1.shape
        metrics_path = directory / "design_space_metrics.json"
        self.metrics: dict[str, Any] = (
            json.loads(metrics_path.read_text()) if metrics_path.exists() else {}
        )
        self._tri_iu = np.triu_indices(self.n_variables)
        n_girth_out = MESH_RESOLUTION["n_girth_out"]
        n_station_out = MESH_RESOLUTION["n_station_out"]
        self._tris, self._fore_idx, self._aft_idx = _connectivity(n_girth_out, n_station_out)
        n_patch_vertices = 2 * n_girth_out * n_station_out
        self._n_patch_vertices = n_patch_vertices
        # the centreline is a topological fact of the layout (last girth row), verified once
        # against the base cloud so a layout change cannot pass unnoticed
        star = np.zeros(n_patch_vertices, dtype=bool)
        star.reshape(2, n_girth_out, n_station_out)[:, -1, :] = True
        base_star = point_grid(self.base_cloud)
        if float(np.abs(base_star[:, -1, :, 1]).max()) > 1e-6:
            raise ValueError("base cloud centreline row is not at y = 0")
        self._starboard_centre = star
        off = np.flatnonzero(~star)
        vmap = np.arange(n_patch_vertices)
        vmap[off] = n_patch_vertices + np.arange(off.size)
        mirrored = np.stack([vmap[self._tris[:, 0]], vmap[self._tris[:, 2]],
                             vmap[self._tris[:, 1]]], axis=1)
        self._all_tris = np.vstack([self._tris, mirrored])
        self._n_vertices = n_patch_vertices + off.size

    # -- design variables ------------------------------------------------------
    def zero(self) -> Any:
        """The reference design vector (all 27 variables zero)."""
        import numpy as np

        return np.zeros(self.n_variables)

    def _quadratic(self, v: Any) -> Any:
        import numpy as np

        return (v[:, None] * v[None, :])[self._tri_iu]

    def _amplitudes(self, v: Any) -> Any:
        """Modal amplitudes for a 27-vector ``v``, shape ``(n_modes,)``.

        Model: ``a(v) = C2 @ [v, vech(v v^T)]``.  ``C2`` was least-squares fitted on the
        full 405-column feature matrix, whose **first 27 columns are ``v`` itself** — so
        the linear part is already inside ``C2``.  ``C1`` (``design_space.npz: c1``) is the
        separate linear-only model M1 kept for diagnostics; adding ``C1 v`` here would
        count the linear term twice (round-0 defect: holdout RMS 0.833 m instead of
        0.252 m and 14/20 sampled configurations failing the geometry gate).
        """
        import numpy as np

        v = self._check(v)
        features = np.concatenate([v, self._quadratic(v)])       # 27 + 378 = 405, as fitted
        return self.c2 @ features

    def _d_amplitudes(self, v: Any) -> Any:
        """Exact ``d a / d v``, shape ``(n_modes, 27)``."""
        import numpy as np

        v = self._check(v)
        dq = np.zeros((self.n_variables * self.n_variables, self.n_variables))
        for k in range(self.n_variables):
            e = np.zeros(self.n_variables)
            e[k] = 1.0
            dq[:, k] = (np.outer(e, v) + np.outer(v, e)).ravel()
        iu, ju = self._tri_iu
        d_features = np.vstack([np.eye(self.n_variables), dq[iu * self.n_variables + ju]])
        # C2 already contains the linear block (see _amplitudes); C1 must not be added.
        return self.c2 @ d_features

    def _check(self, v: Any) -> Any:
        import numpy as np

        v = np.asarray(v, dtype=np.float64).reshape(-1) if not isinstance(v, np.ndarray) else v
        if v.size != self.n_variables:
            raise ValueError(f"design vector must have {self.n_variables} entries, got {v.size}")
        if not bool(np.all(np.isfinite(v))):
            raise ValueError("design vector contains non-finite entries")
        return v

    # -- geometry --------------------------------------------------------------
    def cloud(self, v: Sequence[float] | None = None) -> Any:
        """The 2250 reconstructed points, ``(2250, 3)``; ``v = 0`` is the reference hull."""
        import numpy as np

        v = self.zero() if v is None else self._check(v)
        return self.base_cloud + (self.phi @ self._amplitudes(v)).reshape(-1, 3)

    def cloud_jacobian(self, v: Sequence[float] | None = None) -> Any:
        """Exact ``d cloud / d v``, shape ``(2250, 3, 27)``."""
        v = self.zero() if v is None else v
        return (self.phi @ self._d_amplitudes(v)).reshape(-1, 3, self.n_variables)

    def station_grid(self, v: Sequence[float] | None = None) -> Any:
        """Refined starboard patch grid ``(2, n_girth_out, n_station_out, 3)``."""
        import numpy as np

        grid = point_grid(self.cloud(v))
        return np.stack([_refine(grid[k]) for k in range(2)], axis=0)

    def vertices(self, v: Sequence[float] | None = None) -> Any:
        """Mesh vertices ``(n_vertices, 3)``: starboard patches then their mirror image."""
        import numpy as np

        grid = self.station_grid(v)
        star = np.concatenate([grid[0].reshape(-1, 3), grid[1].reshape(-1, 3)], axis=0)
        off = np.flatnonzero(~self._starboard_centre)
        port = star[off].copy()
        port[:, 1] *= -1.0
        return np.vstack([star, port])

    def faces(self) -> Any:
        """Fixed triangle list ``(n_faces, 3)`` indexing :meth:`vertices`."""
        return self._all_tris

    def vertex_jacobian(self, v: Sequence[float] | None = None) -> Any:
        """Exact ``d vertices / d v``, shape ``(n_vertices, 3, 27)``."""
        import numpy as np

        if v is None:
            v = self.zero()
        dcloud = self.cloud_jacobian(v)                     # (2250, 3, 27)
        arr = np.stack([point_grid(dcloud[:, :, k], check=False)
                        for k in range(self.n_variables)],
                       axis=0)                              # (27, 2, girth, station, 3)
        arr = np.transpose(arr, (1, 2, 3, 4, 0)).reshape(
            GRID["n_patches"], GRID["n_girth"], GRID["n_station"], 3 * self.n_variables)
        ref = np.stack([_refine(arr[k]) for k in range(2)], axis=0).reshape(
            2, MESH_RESOLUTION["n_girth_out"], MESH_RESOLUTION["n_station_out"],
            3, self.n_variables)
        star = np.concatenate([ref[0].reshape(-1, 3, self.n_variables),
                               ref[1].reshape(-1, 3, self.n_variables)], axis=0)
        off = np.flatnonzero(~self._starboard_centre)
        port = star[off].copy()
        port[:, 1, :] *= -1.0
        return np.vstack([star, port])

    def mesh(self, v: Sequence[float] | None = None) -> Mesh:
        """The hull shell as a :class:`tensorcad.geometry.Mesh` (open at the deck edge)."""
        verts = self.vertices(v)
        return Mesh(
            vertices=tuple((float(a), float(b), float(c)) for a, b, c in verts),
            faces=tuple((int(a), int(b), int(c)) for a, b, c in self._all_tris),
            units="m",
            coordinate_frame="hull frame: x aft->fore, y starboard, z up; full scale",
        )

    # -- ordered station arrays (one longitudinal sweep, +x) --------------------
    def _stations(self, grid: Any) -> tuple[Any, Any, Any]:
        """``(x, y, z)`` station arrays ordered along +x: transom -> seam -> stem.

        ``x`` is ``(n_station_total,)`` (girth-mean station position), ``y``/``z`` are
        ``(n_station_total, n_girth)`` half-section contours (deck edge -> keel).
        """
        import numpy as np

        n_girth_out = grid.shape[1]
        n_station_out = grid.shape[2]
        sa = np.arange(n_station_out - 1, -1, -1)          # both patches run with x descending
        x = np.concatenate([grid[1][:, sa, 0].mean(axis=0), grid[0][:, sa, 0].mean(axis=0)])
        y = np.concatenate([grid[1][:, sa, 1].T, grid[0][:, sa, 1].T], axis=0)
        z = np.concatenate([grid[1][:, sa, 2].T, grid[0][:, sa, 2].T], axis=0)
        assert y.shape == (2 * n_station_out, n_girth_out)
        return x, y, z, sa

    @staticmethod
    def _section_area(y: Any, z: Any) -> Any:
        """Half-section area to the deck edge, closed along the centreline (shoelace)."""
        import numpy as np

        yc = np.concatenate([y, np.zeros((y.shape[0], 1))], axis=1)
        zc = np.concatenate([z, z[:, :1]], axis=1)
        dy = np.diff(yc, axis=1, append=yc[:, :1])
        dz = np.diff(zc, axis=1, append=zc[:, :1])
        return 0.5 * np.abs((yc * dz - zc * dy).sum(axis=1))

    @staticmethod
    def _section_area_partials(y: Any, z: Any) -> tuple[Any, Any]:
        """``d A / d y`` and ``d A / d z`` (shape ``(n_station, n_girth)``) of the shoelace area.

        The starboard half-section is ordered deck edge -> keel -> centreline chord, so the
        signed shoelace sum is positive on the design envelope; the sign is folded in so the
        derivative matches the absolute value returned by :meth:`_section_area`.
        """
        import numpy as np

        yc = np.concatenate([y, np.zeros((y.shape[0], 1))], axis=1)
        zc = np.concatenate([z, z[:, :1]], axis=1)
        n = yc.shape[1]
        dy = np.diff(yc, axis=1, append=yc[:, :1])
        dz = np.diff(zc, axis=1, append=zc[:, :1])
        signed = 0.5 * (yc * dz - zc * dy).sum(axis=1)
        sign = np.where(signed >= 0.0, 1.0, -1.0)[:, None]
        dA_dy = np.zeros_like(yc)
        dA_dz = np.zeros_like(zc)
        for k in range(n):
            kp = (k + 1) % n
            km = (k - 1) % n
            dA_dy[:, k] = 0.5 * (zc[:, kp] - zc[:, km])
            dA_dz[:, k] = 0.5 * (yc[:, km] - yc[:, kp])
        # the appended closing point is (y = 0, z = z_deck): its z is not independent, so its
        # derivative has to be folded back into the deck-edge point before dropping it
        dA_dz[:, 0] = dA_dz[:, 0] + dA_dz[:, -1]
        return sign * dA_dy[:, :-1], sign * dA_dz[:, :-1]

    # -- differentiable integral ----------------------------------------------
    def volume(self, v: Sequence[float] | None = None) -> float:
        """Enclosed shell volume ``2 * integral A(x) dx`` in m^3 (deck-edge closed)."""
        import numpy as np

        x, y, z, _ = self._stations(self.station_grid(v))
        area = self._section_area(y, z)
        dx = np.diff(x)
        if bool(np.any(dx <= 0.0)):
            raise ValueError("station x is not monotone -- design vector leaves the envelope")
        return float(2.0 * (0.5 * (area[:-1] + area[1:]) * dx).sum())

    def volume_gradient(self, v: Sequence[float] | None = None) -> Any:
        """Exact ``d volume / d v``, shape ``(27,)``.

        Chain: shoelace area partials -> trapezoid weights (including the ``dx`` weights,
        because the station x positions move with ``v`` as well) -> station vertices ->
        ``d vertices / d v``.
        """
        import numpy as np

        if v is None:
            v = self.zero()
        v = self._check(v)
        grid = self.station_grid(v)
        n_station_out = grid.shape[2]
        n_girth_out = grid.shape[1]
        x, y, z, sa = self._stations(grid)
        area = self._section_area(y, z)
        dx = np.diff(x)
        if bool(np.any(dx <= 0.0)):
            raise ValueError("station x is not monotone -- design vector leaves the envelope")
        dA_dy, dA_dz = self._section_area_partials(y, z)
        dV_dA = np.zeros(area.size)
        dV_dA[:-1] += dx
        dV_dA[1:] += dx
        dV_dX = np.zeros(area.size)
        dV_dX[:-1] -= 0.5 * (area[:-1] + area[1:])
        dV_dX[1:] += 0.5 * (area[:-1] + area[1:])
        G = np.zeros((area.size, n_girth_out, 3))
        # V = 2 * sum_i 0.5 (A_i + A_i+1) dx_i, so the mirror factor 2 is already inside the
        # trapezoid weights dV_dA = dx_(i-1) + dx_i; it still has to be applied to the x
        # weights dV_dX = +-0.5 (A_i + A_i+1), which were written with the 0.5 kept explicit.
        # (Round 0 applied `G *= 2.0` to everything: dV/dv came out exactly 2x the central
        # difference -- see tools/dtmb5415_param/run_acceptance.py, acceptance 3.)
        G[:, :, 0] = 2.0 * dV_dX[:, None] / float(n_girth_out)
        G[:, :, 1] = dV_dA[:, None] * dA_dy
        G[:, :, 2] = dV_dA[:, None] * dA_dz
        # scatter back to (patch, station) layout and then to the flattened vertex order
        G_patch = np.empty((2, n_station_out, n_girth_out, 3))
        G_patch[1][sa] = G[:n_station_out]
        G_patch[0][sa] = G[n_station_out:]
        flat = np.concatenate([G_patch[0].transpose(1, 0, 2).reshape(-1, 3),
                               G_patch[1].transpose(1, 0, 2).reshape(-1, 3)], axis=0)
        jac = self.vertex_jacobian(v)[: self._n_patch_vertices]
        return np.einsum("nc,nck->k", flat, jac)

    # -- physical (not differentiated) hydrostatics ---------------------------
    def hydrostatics(self, v: Sequence[float] | None = None,
                     waterline_z: float = WATERLINE_Z) -> dict[str, Any]:
        """Displacement, waterplane, beam and LCB at ``waterline_z`` (full scale, m).

        The waterline clip is piecewise smooth, so this is deliberately *not* differentiated:
        the analytic chain stops at :meth:`volume`.  Integrals use the refined station grid
        with the trapezoidal rule in ``x`` and an exact polygon clip in ``z`` per section.
        """
        import numpy as np

        x, y, z, _ = self._stations(self.station_grid(v))
        n_st = x.size
        a_sub = np.zeros(n_st)
        y_wl = np.zeros(n_st)
        for s in range(n_st):
            a_sub[s], y_wl[s] = _clip_section(y[s], z[s], waterline_z)
        dx = np.diff(x)
        if bool(np.any(dx <= 0.0)):
            raise ValueError("station x is not monotone -- design vector leaves the envelope")
        vol = float(2.0 * (0.5 * (a_sub[:-1] + a_sub[1:]) * dx).sum())
        awp = float(2.0 * (0.5 * (y_wl[:-1] + y_wl[1:]) * dx).sum())
        moment = float(2.0 * (0.5 * (x[:-1] * a_sub[:-1] + x[1:] * a_sub[1:]) * dx).sum())
        length = float(x[-1] - x[0])
        beam = float(2.0 * y_wl.max())
        t_wl = float(waterline_z)   # baseline (z = 0) is the moulded keel; the bulb
        # hangs below it (z_min = -3.28) and is not part of the moulded draft
        return {
            "waterline_z_m": float(waterline_z),
            "length_stations_m": length,
            "beam_waterline_m": beam,
            "draft_m": t_wl,
            "displacement_volume_m3": vol,
            "waterplane_area_m2": awp,
            "lcb_m": float(moment / vol) if vol > 0.0 else float("nan"),
            "cb": float(vol / (length * beam * t_wl)) if vol > 0.0 else float("nan"),
            "cwp": float(awp / (length * beam)) if vol > 0.0 else float("nan"),
            "keel_z_m": -3.279432,
        }

    # -- mesh quality (the numbers a gate needs) ------------------------------
    def quality(self, v: Sequence[float] | None = None) -> dict[str, Any]:
        """Topology/quality metrics of the mesh at ``v`` (all in the hull frame, metres)."""
        import numpy as np

        verts = self.vertices(v)
        tris = self._all_tris
        p = verts[tris]
        e0 = p[:, 1] - p[:, 0]
        e1 = p[:, 2] - p[:, 1]
        e2 = p[:, 0] - p[:, 2]
        area = 0.5 * np.linalg.norm(np.cross(e0, -e2), axis=1)
        lens = np.stack([np.linalg.norm(e0, axis=1), np.linalg.norm(e1, axis=1),
                         np.linalg.norm(e2, axis=1)], axis=1)
        perimeter = lens.sum(axis=1)
        # smallest angle via the two shortest edges of each triangle
        order = np.argsort(lens, axis=1)
        short = np.take_along_axis(lens, order[:, :2], axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            cos_min = ((short[:, 0] ** 2 + short[:, 1] ** 2 - np.take_along_axis(
                lens, order[:, 2:], axis=1)[:, 0] ** 2) / (2 * short[:, 0] * short[:, 1]))
        min_angle = np.degrees(np.arccos(np.clip(cos_min, -1.0, 1.0)))
        with np.errstate(divide="ignore", invalid="ignore"):
            aspect = perimeter ** 2 / (12.0 * math.sqrt(3.0) * area)
        edges = np.concatenate([tris[:, [0, 1]], tris[:, [1, 2]], tris[:, [2, 0]]], axis=0)
        edges = np.sort(edges, axis=1)
        _, counts = np.unique(edges, axis=0, return_counts=True)
        # station monotonicity per patch: x must decrease with station index
        grid = self.station_grid(v)
        dx_station = np.diff(grid[:, :, :, 0].mean(axis=1), axis=1)
        return {
            "n_vertices": int(verts.shape[0]),
            "n_faces": int(tris.shape[0]),
            "edge_median_m": float(np.median(lens.ravel())),
            "edge_min_m": float(lens.min()),
            "edge_max_m": float(lens.max()),
            "area_min_m2": float(area.min()),
            "degenerate_faces": int((area < 1e-9).sum()),
            "min_angle_deg": float(np.nanmin(min_angle)),
            "min_angle_p01_deg": float(np.nanpercentile(min_angle, 1)),
            "frac_angle_below_5deg": float((min_angle < 5.0).mean()),
            "aspect_p99": float(np.nanpercentile(aspect, 99)),
            "aspect_max": float(np.nanmax(aspect)),
            "non_manifold_edges": int((counts > 2).sum()),
            "boundary_edges": int((counts == 1).sum()),
            "non_finite_vertices": int((~np.isfinite(verts)).any(axis=1).sum()),
            "half_breadth_min_m": float(grid[:, :, :, 1].min()),
            "stations_monotone": bool(np.all(dx_station < 0.0)),
            "bbox": [[float(a), float(b)] for a, b in zip(verts.min(axis=0), verts.max(axis=0))],
        }

    def to_dict(self, v: Sequence[float] | None = None) -> dict[str, Any]:
        """A JSON-ready record: schema, provenance, design vector, geometry and metrics."""
        import numpy as np

        if v is None:
            v = self.zero()
        return {
            "schema": SPACE_SCHEMA,
            "dataset": {"title": DATASET["title"], "authors": DATASET["authors"],
                        "n_configurations": DATASET["n_configurations"],
                        "n_design_variables": DATASET["n_design_variables"],
                        "source_sha256": DATASET["source_sha256"]},
            "reference": {"path": REFERENCE["path"], "sha256": REFERENCE["sha256"]},
            "design_vector": [float(x) for x in np.asarray(v).reshape(-1)],
            "n_modes": int(self.n_modes),
            "metrics": self.metrics,
            "hydrostatics": self.hydrostatics(v),
            "quality": self.quality(v),
        }


def _clip_section(y: Any, z: Any, waterline_z: float) -> tuple[float, float]:
    """Area below ``waterline_z`` and half-breadth at that plane for one half-section.

    The half-section is closed along the centreline (deck edge -> keel -> straight chord back
    up the centreline), clipped to ``z <= waterline_z`` with Sutherland-Hodgman, then measured
    with the shoelace formula.  Returns ``(area, y_at_waterline)``.
    """
    import numpy as np

    poly = np.stack([np.concatenate([y, [0.0]]), np.concatenate([z, [z[0]]])], axis=1)
    out: list[tuple[float, float]] = []
    n = poly.shape[0]
    for i in range(n):
        cur = poly[i]
        nxt = poly[(i + 1) % n]
        cur_in = cur[1] <= waterline_z
        nxt_in = nxt[1] <= waterline_z
        if cur_in:
            out.append((float(cur[0]), float(cur[1])))
        if cur_in != nxt_in:
            t = (waterline_z - cur[1]) / (nxt[1] - cur[1])
            out.append((float(cur[0] + t * (nxt[0] - cur[0])), float(waterline_z)))
    if len(out) < 3:
        area = 0.0
    else:
        pts = np.asarray(out)
        yy, zz = pts[:, 0], pts[:, 1]
        dy = np.diff(yy, append=yy[:1])
        dz = np.diff(zz, append=zz[:1])
        area = float(abs(0.5 * (yy * dz - zz * dy).sum()))
    # half-breadth at the waterline: widest crossing along the girth
    y_wl = 0.0
    for i in range(n - 1):
        z0, z1 = poly[i, 1], poly[i + 1, 1]
        if (z0 - waterline_z) * (z1 - waterline_z) < 0.0:
            t = (waterline_z - z0) / (z1 - z0)
            y_wl = max(y_wl, float(poly[i, 0] + t * (poly[i + 1, 0] - poly[i, 0])))
    return area, y_wl
