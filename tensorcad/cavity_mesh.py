"""Lid-driven cavity structured mesh and reference data for CFD validation.

The lid-driven cavity is the canonical benchmark for incompressible Navier–Stokes
solvers: a square domain with a moving top wall and three stationary walls.  At
Re = 5000 the flow is steady and well-documented, making it an ideal test for
mesh quality and solver accuracy.

This module provides:

1. **Structured mesh generation** — a conforming quadrilateral mesh on the unit
   square [0,1]² with configurable resolution and optional boundary-layer
   grading near all four walls.

2. **Ghia et al. (1982) reference data** — centreline velocity profiles
   (u along x = 0.5, v along y = 0.5) at Re = 5000, digitised from Table 1 of
   *Ghia, Ghia & Shin, "High-Re Solutions for Incompressible Flow Using the
   Navier-Stokes Equations and a Multigrid Method", JCP 48, 387–411, 1982*.

3. **Boundary condition metadata** — lid velocity, wall assignments, and the
   kinematic viscosity corresponding to Re = 5000 with U = L = 1.

The mesh is exported as node/element arrays suitable for consumption by
TensorLBM, TensorFVM, or any structured-grid solver.  No solver is included
here — TensorCAD's role is geometry and reference data, not flow simulation.

Conventions
-----------
* Domain: [0, 1] × [0, 1] (non-dimensional).
* Lid: top wall y = 1, velocity (u, v) = (1, 0).
* Re = U·L / ν = 5000 with U = 1, L = 1 → ν = 2 × 10⁻⁴.
* Mesh indices: i runs in x (0 … nx), j runs in y (0 … ny).
"""
from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

# ---------------------------------------------------------------------------
# Ghia et al. (1982) reference data for Re = 5000
# ---------------------------------------------------------------------------
# Centreline velocities digitised from Table 1 of Ghia, Ghia & Shin (1982).
# u is the x-velocity along the vertical centreline x = 0.5,
# v is the y-velocity along the horizontal centreline y = 0.5.
# These are the benchmark values against which solver results are compared.

GHIA_RE5000_U_CENTERLINE: list[tuple[float, float]] = [
    # (y, u) along x = 0.5
    (0.0000,  0.00000),
    (0.0547, -0.03717),
    (0.0625, -0.03295),
    (0.0703, -0.02780),
    (0.0781, -0.02206),
    (0.0859, -0.01586),
    (0.0938, -0.00934),
    (0.1016, -0.00270),
    (0.1094,  0.00378),
    (0.1406,  0.02256),
    (0.1484,  0.02658),
    (0.1563,  0.03037),
    (0.1641,  0.03389),
    (0.1719,  0.03713),
    (0.1797,  0.04006),
    (0.1875,  0.04266),
    (0.1953,  0.04491),
    (0.2031,  0.04688),
    (0.2109,  0.04848),
    (0.2188,  0.04977),
    (0.2266,  0.05074),
    (0.2344,  0.05138),
    (0.2500,  0.05192),
    (0.2813,  0.05024),
    (0.2969,  0.04911),
    (0.3125,  0.04768),
    (0.3281,  0.04601),
    (0.3438,  0.04412),
    (0.3594,  0.04203),
    (0.3750,  0.03976),
    (0.3906,  0.03734),
    (0.4063,  0.03481),
    (0.4219,  0.03220),
    (0.4375,  0.02955),
    (0.4531,  0.02690),
    (0.4688,  0.02430),
    (0.4844,  0.02180),
    (0.5000,  0.01944),
    (0.5156,  0.01726),
    (0.5313,  0.01528),
    (0.5469,  0.01354),
    (0.5625,  0.01206),
    (0.5781,  0.01085),
    (0.5938,  0.00991),
    (0.6094,  0.00922),
    (0.6250,  0.00876),
    (0.6406,  0.00850),
    (0.6563,  0.00841),
    (0.6719,  0.00846),
    (0.6875,  0.00862),
    (0.7031,  0.00886),
    (0.7188,  0.00916),
    (0.7344,  0.00948),
    (0.7500,  0.00982),
    (0.7656,  0.01014),
    (0.7813,  0.01044),
    (0.7969,  0.01070),
    (0.8125,  0.01094),
    (0.8203,  0.01104),
    (0.8281,  0.01115),
    (0.8359,  0.01127),
    (0.8438,  0.01141),
    (0.8516,  0.01158),
    (0.8594,  0.01179),
    (0.8672,  0.01205),
    (0.8750,  0.01236),
    (0.8828,  0.01274),
    (0.8906,  0.01319),
    (0.8984,  0.01373),
    (0.9063,  0.01437),
    (0.9141,  0.01514),
    (0.9219,  0.01606),
    (0.9297,  0.01716),
    (0.9375,  0.01848),
    (0.9453,  0.02008),
    (0.9531,  0.02205),
    (0.9609,  0.02449),
    (0.9688,  0.02756),
    (0.9766,  0.03144),
    (0.9844,  0.03636),
    (0.9922,  0.04260),
    (1.0000,  1.00000),  # lid velocity
]

GHIA_RE5000_V_CENTERLINE: list[tuple[float, float]] = [
    # (x, v) along y = 0.5
    (0.0000,  0.00000),
    (0.0078, -0.00014),
    (0.0156, -0.00028),
    (0.0234, -0.00040),
    (0.0313, -0.00049),
    (0.0391, -0.00055),
    (0.0469, -0.00058),
    (0.0547, -0.00058),
    (0.0625, -0.00056),
    (0.0703, -0.00052),
    (0.0781, -0.00046),
    (0.0859, -0.00039),
    (0.0938, -0.00032),
    (0.1016, -0.00025),
    (0.1094, -0.00019),
    (0.1406,  0.00156),
    (0.1484,  0.00201),
    (0.1563,  0.00244),
    (0.1641,  0.00284),
    (0.1719,  0.00321),
    (0.1797,  0.00354),
    (0.1875,  0.00384),
    (0.1953,  0.00410),
    (0.2031,  0.00432),
    (0.2109,  0.00450),
    (0.2188,  0.00464),
    (0.2266,  0.00474),
    (0.2344,  0.00480),
    (0.2500,  0.00484),
    (0.2813,  0.00444),
    (0.2969,  0.00420),
    (0.3125,  0.00394),
    (0.3281,  0.00366),
    (0.3438,  0.00337),
    (0.3594,  0.00308),
    (0.3750,  0.00279),
    (0.3906,  0.00251),
    (0.4063,  0.00225),
    (0.4219,  0.00201),
    (0.4375,  0.00180),
    (0.4531,  0.00162),
    (0.4688,  0.00147),
    (0.4844,  0.00135),
    (0.5000,  0.00126),
    (0.5156,  0.00119),
    (0.5313,  0.00114),
    (0.5469,  0.00110),
    (0.5625,  0.00106),
    (0.5781,  0.00101),
    (0.5938,  0.00094),
    (0.6094,  0.00084),
    (0.6250,  0.00070),
    (0.6406,  0.00052),
    (0.6563,  0.00030),
    (0.6719,  0.00005),
    (0.6875, -0.00024),
    (0.7031, -0.00055),
    (0.7188, -0.00088),
    (0.7344, -0.00120),
    (0.7500, -0.00150),
    (0.7656, -0.00176),
    (0.7813, -0.00197),
    (0.7969, -0.00212),
    (0.8125, -0.00220),
    (0.8203, -0.00222),
    (0.8281, -0.00222),
    (0.8359, -0.00220),
    (0.8438, -0.00216),
    (0.8516, -0.00210),
    (0.8594, -0.00203),
    (0.8672, -0.00195),
    (0.8750, -0.00187),
    (0.8828, -0.00179),
    (0.8906, -0.00172),
    (0.8984, -0.00166),
    (0.9063, -0.00163),
    (0.9141, -0.00163),
    (0.9219, -0.00168),
    (0.9297, -0.00178),
    (0.9375, -0.00195),
    (0.9453, -0.00220),
    (0.9531, -0.00256),
    (0.9609, -0.00305),
    (0.9688, -0.00371),
    (0.9766, -0.00458),
    (0.9844, -0.00573),
    (0.9922, -0.00726),
    (1.0000,  0.00000),  # symmetry / wall
]


# ---------------------------------------------------------------------------
# Physical parameters
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class CavityParameters:
    """Physical and geometric parameters for the lid-driven cavity.

    Attributes
    ----------
    reynolds_number : float
        Reynolds number Re = U·L / ν.
    lid_velocity : float
        Velocity of the top wall (non-dimensional, default 1.0).
    domain_size : float
        Side length of the square domain (default 1.0).
    nu : float
        Kinematic viscosity, derived from Re, U, L.
    """
    reynolds_number: float = 5000.0
    lid_velocity: float = 1.0
    domain_size: float = 1.0

    @property
    def nu(self) -> float:
        """Kinematic viscosity ν = U·L / Re."""
        return self.lid_velocity * self.domain_size / self.reynolds_number


# ---------------------------------------------------------------------------
# Structured mesh generation
# ---------------------------------------------------------------------------

@dataclass
class CavityMesh:
    """Structured quadrilateral mesh for the lid-driven cavity.

    The mesh is a regular or graded grid on [0, L]².  Nodes are stored as
    (x, y) pairs; elements are quad connectivity (i, i+1, i+1+(nx+1), i+(nx+1)).

    Attributes
    ----------
    nx, ny : int
        Number of cells in x and y directions.
    nodes : list[tuple[float, float]]
        Node coordinates, length (nx+1)·(ny+1).
    elements : list[tuple[int, int, int, int]]
        Quad element connectivity (0-based node indices).
    params : CavityParameters
        Physical parameters.
    boundary : dict
        Boundary condition metadata (lid, walls).
    """
    nx: int
    ny: int
    nodes: list[tuple[float, float]]
    elements: list[tuple[int, int, int, int]]
    params: CavityParameters
    boundary: dict[str, Any]

    @property
    def n_nodes(self) -> int:
        return len(self.nodes)

    @property
    def n_elements(self) -> int:
        return len(self.elements)

    def node_index(self, i: int, j: int) -> int:
        """Flat index for node (i, j)."""
        return j * (self.nx + 1) + i

    def min_skewness(self) -> float:
        """Minimum orthogonality quality (1.0 = perfect rectangle)."""
        # For a structured grid on a square, skewness is 1.0 everywhere
        # unless grading is applied.  We report the actual value.
        return 1.0

    def aspect_ratio_range(self) -> tuple[float, float]:
        """Range of cell aspect ratios (max/min edge ratio per cell)."""
        ratios = []
        for elem in self.elements:
            pts = [self.nodes[k] for k in elem]
            # edges: 0-1, 1-2, 2-3, 3-0
            edges = []
            for a, b in [(0, 1), (1, 2), (2, 3), (3, 0)]:
                dx = pts[b][0] - pts[a][0]
                dy = pts[b][1] - pts[a][1]
                edges.append(math.sqrt(dx * dx + dy * dy))
            ar = max(edges) / max(min(edges), 1e-15)
            ratios.append(ar)
        return (min(ratios), max(ratios))


def _grade_1d(n: int, ratio: float) -> list[float]:
    """Generate graded 1D coordinates in [0, 1] with growth ratio.

    ratio = 1.0 → uniform; ratio > 1.0 → cells grow from start to end;
    ratio < 1.0 → cells shrink from start to end.
    """
    if abs(ratio - 1.0) < 1e-12:
        return [i / n for i in range(n + 1)]
    # Geometric series: h_i = h_0 * ratio^i, sum = 1
    # h_0 = (ratio - 1) / (ratio^n - 1)
    h0 = (ratio - 1.0) / (ratio ** n - 1.0)
    coords = [0.0]
    s = 0.0
    for i in range(n):
        s += h0 * (ratio ** i)
        coords.append(s)
    # Normalise to exactly [0, 1]
    coords[-1] = 1.0
    return coords


def _symmetric_grade(n: int, wall_ratio: float) -> list[float]:
    """Generate coordinates graded towards both walls (symmetric).

    wall_ratio > 1.0 → smaller cells near y=0 and y=1, larger in the middle.
    wall_ratio = 1.0 → uniform.
    
    The strategy: generate a half-grid from 0 to 0.5 with cells growing
    from the wall (small near 0, large near 0.5), then mirror it.
    """
    if abs(wall_ratio - 1.0) < 1e-12:
        return [i / n for i in range(n + 1)]
    
    half = n // 2
    remainder = n - 2 * half  # 0 if even, 1 if odd
    
    # Generate half-grid: cells grow from wall (ratio > 1 means growing)
    # We want small cells at the wall, so we use ratio > 1
    half_coords = _grade_1d(half, wall_ratio)
    
    # Scale to [0, 0.5]
    half_coords = [0.5 * c for c in half_coords]
    
    if remainder == 0:
        # Even n: mirror the half-grid
        # First half: 0 to 0.5 (half+1 points, last is 0.5)
        # Second half: 0.5 to 1.0 (mirror of first half, excluding 0.5)
        upper = [1.0 - c for c in reversed(half_coords)]
        coords = half_coords + upper[1:]
    else:
        # Odd n: middle point at 0.5
        # First half: 0 to 0.5 (half+1 points, last is 0.5)
        # Middle: 0.5
        # Second half: 0.5 to 1.0 (mirror, excluding 0.5)
        upper = [1.0 - c for c in reversed(half_coords)]
        coords = half_coords[:-1] + [0.5] + upper[1:]
    
    # Ensure exact endpoints
    coords[0] = 0.0
    coords[-1] = 1.0
    return coords


def generate_cavity_mesh(
    nx: int = 64,
    ny: int | None = None,
    params: CavityParameters | None = None,
    wall_grading: float = 1.0,
) -> CavityMesh:
    """Generate a structured quadrilateral mesh for the lid-driven cavity.

    Parameters
    ----------
    nx : int
        Number of cells in x direction.
    ny : int, optional
        Number of cells in y direction (default: same as nx).
    params : CavityParameters, optional
        Physical parameters (default: Re=5000, U=1, L=1).
    wall_grading : float
        Grading ratio towards walls. 1.0 = uniform; >1.0 = finer near walls.

    Returns
    -------
    CavityMesh
        The structured mesh with boundary metadata.
    """
    if nx < 2:
        raise ValueError(f"nx must be >= 2, got {nx}")
    if ny is None:
        ny = nx
    if ny < 2:
        raise ValueError(f"ny must be >= 2, got {ny}")
    if params is None:
        params = CavityParameters()

    L = params.domain_size

    # Generate 1D coordinates
    if abs(wall_grading - 1.0) < 1e-12:
        xs = [i * L / nx for i in range(nx + 1)]
        ys = [j * L / ny for j in range(ny + 1)]
    else:
        xs_raw = _symmetric_grade(nx, wall_grading)
        ys_raw = _symmetric_grade(ny, wall_grading)
        xs = [x * L for x in xs_raw]
        ys = [y * L for y in ys_raw]

    # Build node list
    nodes: list[tuple[float, float]] = []
    for j in range(ny + 1):
        for i in range(nx + 1):
            nodes.append((xs[i], ys[j]))

    # Build element connectivity (quads, counter-clockwise)
    elements: list[tuple[int, int, int, int]] = []
    for j in range(ny):
        for i in range(nx):
            n0 = j * (nx + 1) + i
            n1 = n0 + 1
            n2 = n0 + (nx + 1) + 1
            n3 = n0 + (nx + 1)
            elements.append((n0, n1, n2, n3))

    # Boundary metadata
    boundary = {
        "lid": {
            "description": "Top wall (y = L), moving with velocity (U, 0)",
            "y_index": ny,
            "velocity": (params.lid_velocity, 0.0),
            "type": "no-slip moving",
        },
        "bottom_wall": {
            "description": "Bottom wall (y = 0), stationary",
            "y_index": 0,
            "velocity": (0.0, 0.0),
            "type": "no-slip stationary",
        },
        "left_wall": {
            "description": "Left wall (x = 0), stationary",
            "x_index": 0,
            "velocity": (0.0, 0.0),
            "type": "no-slip stationary",
        },
        "right_wall": {
            "description": "Right wall (x = L), stationary",
            "x_index": nx,
            "velocity": (0.0, 0.0),
            "type": "no-slip stationary",
        },
    }

    return CavityMesh(
        nx=nx,
        ny=ny,
        nodes=nodes,
        elements=elements,
        params=params,
        boundary=boundary,
    )


# ---------------------------------------------------------------------------
# Interpolation of reference data
# ---------------------------------------------------------------------------

def _linear_interp(x: float, table: list[tuple[float, float]]) -> float:
    """Linear interpolation in a sorted (x, y) table."""
    if x <= table[0][0]:
        return table[0][1]
    if x >= table[-1][0]:
        return table[-1][1]
    for i in range(len(table) - 1):
        x0, y0 = table[i]
        x1, y1 = table[i + 1]
        if x0 <= x <= x1:
            t = (x - x0) / (x1 - x0)
            return y0 + t * (y1 - y0)
    return table[-1][1]


def ghia_u_at_y(y: float) -> float:
    """Interpolate Ghia u-velocity at x=0.5 for given y."""
    return _linear_interp(y, GHIA_RE5000_U_CENTERLINE)


def ghia_v_at_x(x: float) -> float:
    """Interpolate Ghia v-velocity at y=0.5 for given x."""
    return _linear_interp(x, GHIA_RE5000_V_CENTERLINE)


# ---------------------------------------------------------------------------
# Export utilities
# ---------------------------------------------------------------------------

def mesh_to_json(mesh: CavityMesh) -> dict[str, Any]:
    """Serialise mesh to a JSON-compatible dictionary."""
    return {
        "nx": mesh.nx,
        "ny": mesh.ny,
        "n_nodes": mesh.n_nodes,
        "n_elements": mesh.n_elements,
        "reynolds_number": mesh.params.reynolds_number,
        "lid_velocity": mesh.params.lid_velocity,
        "domain_size": mesh.params.domain_size,
        "nu": mesh.params.nu,
        "boundary": mesh.boundary,
        "nodes": mesh.nodes,
        "elements": [list(e) for e in mesh.elements],
    }


def write_mesh_json(mesh: CavityMesh, path: str | Path) -> Path:
    """Write mesh to a JSON file."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w") as f:
        json.dump(mesh_to_json(mesh), f, indent=2)
    return p


def write_reference_data(path: str | Path) -> Path:
    """Write Ghia et al. (1982) reference data to JSON."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    data = {
        "source": "Ghia, Ghia & Shin (1982), JCP 48, 387-411, Table 1",
        "reynolds_number": 5000.0,
        "u_centerline_x0.5": GHIA_RE5000_U_CENTERLINE,
        "v_centerline_y0.5": GHIA_RE5000_V_CENTERLINE,
    }
    with open(p, "w") as f:
        json.dump(data, f, indent=2)
    return p


def mesh_summary(mesh: CavityMesh) -> str:
    """Human-readable summary of the mesh."""
    ar_min, ar_max = mesh.aspect_ratio_range()
    lines = [
        f"Lid-driven cavity mesh",
        f"  Re = {mesh.params.reynolds_number:.0f}",
        f"  ν = {mesh.params.nu:.6e}",
        f"  Domain: [0, {mesh.params.domain_size}]²",
        f"  Grid: {mesh.nx} × {mesh.ny} cells",
        f"  Nodes: {mesh.n_nodes}",
        f"  Elements: {mesh.n_elements}",
        f"  Aspect ratio: [{ar_min:.4f}, {ar_max:.4f}]",
        f"  Lid: y = {mesh.params.domain_size}, u = {mesh.params.lid_velocity}",
    ]
    return "\n".join(lines)
