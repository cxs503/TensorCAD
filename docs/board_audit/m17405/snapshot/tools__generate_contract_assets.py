"""Generate 8 STL test assets for the geometry contract compliance table.

Assets
------
1. **wigley.stl** — existing valid Wigley hull (L=10 m, B=1 m, T=0.625 m).
2. **series60_cb60.stl** — valid Series 60 Cb=0.60 hull (L=20 m).
3. **series60_cb75.stl** — valid Series 60 Cb=0.75 hull (L=20 m).
4. **box_barge.stl** — valid rectangular barge (L=30 m, B=8 m, T=4 m).
5. **sphere_r1.stl** — valid unit sphere (r=1 m), closed, watertight.
6. **inverted_wigley.stl** — **red-team**: Wigley with inward normals (volume < 0).
7. **wrong_units.stl** — **red-team**: Wigley geometry labelled as ``mm`` (scale violation).
8. **left_handed.stl** — **red-team**: Wigley with left-handed frame declaration.

All files are written to ``artifacts/geom_contract_assets/``.
"""
from __future__ import annotations

import math
import shutil
from pathlib import Path

from tensorcad.geometry import Mesh, MeshBuilder
from tensorcad.export import write_stl
from tensorcad import wigley, series60, mesh as meshlib


def make_sphere(radius: float = 1.0, n_lat: int = 32, n_lon: int = 64) -> Mesh:
    """Build a watertight sphere mesh using a cube-to-sphere projection."""
    builder = MeshBuilder(units="m", coordinate_frame="x fore, y starboard, z up; right-handed")

    nx = n_lon // 4
    if nx < 1:
        nx = 1

    faces_def = [
        ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
        ((-1, 0, 0), (0, 0, 1), (0, 1, 0)),
        ((0, 1, 0), (0, 0, 1), (1, 0, 0)),
        ((0, -1, 0), (0, 0, -1), (1, 0, 0)),
        ((0, 0, 1), (1, 0, 0), (0, 1, 0)),
        ((0, 0, -1), (0, 1, 0), (1, 0, 0)),
    ]

    for normal, u_dir, v_dir in faces_def:
        for iu in range(nx):
            for iv in range(nx):
                corners_local = [
                    (-1 + 2 * iu / nx, -1 + 2 * iv / nx),
                    (-1 + 2 * (iu + 1) / nx, -1 + 2 * iv / nx),
                    (-1 + 2 * (iu + 1) / nx, -1 + 2 * (iv + 1) / nx),
                    (-1 + 2 * iu / nx, -1 + 2 * (iv + 1) / nx),
                ]
                corners_3d = []
                for u, v in corners_local:
                    x = normal[0] + u_dir[0] * u + v_dir[0] * v
                    y = normal[1] + u_dir[1] * u + v_dir[1] * v
                    z = normal[2] + u_dir[2] * u + v_dir[2] * v
                    length = math.sqrt(x * x + y * y + z * z)
                    corners_3d.append((radius * x / length, radius * y / length, radius * z / length))

                indices = [builder.add_vertex(c) for c in corners_3d]
                builder.add_quad(indices[0], indices[1], indices[2], indices[3])

    mesh = builder.build()
    oriented = meshlib.orient_outward(mesh)
    return oriented.mesh


def make_box_barge(length: float = 30.0, breadth: float = 8.0, draft: float = 4.0,
                   depth: float = 5.0, nx: int = 30, ny: int = 8, nz: int = 4) -> Mesh:
    """Build a closed rectangular barge mesh with shared vertices (watertight)."""
    builder = MeshBuilder(units="m", coordinate_frame="x aft->fore, y port->starboard, z up; right-handed")
    L, B, T, D = length, breadth, draft, depth

    # Build a grid of vertices on the box surface
    # The box has 6 faces. We create vertices at grid points and reference them by index.
    # To ensure watertightness, shared edges must use the same vertex indices.

    # We'll create a parametric mapping from (face, u, v) to vertex index.
    # But a simpler approach: create all vertices first, then create faces.

    # Grid points on each face (including edges)
    # Bottom: z = -T, x in [0, L], y in [-B/2, B/2]
    # Top: z = D-T, x in [0, L], y in [-B/2, B/2]
    # Starboard: y = B/2, x in [0, L], z in [-T, D-T]
    # Port: y = -B/2, x in [0, L], z in [-T, D-T]
    # Transom: x = L, y in [-B/2, B/2], z in [-T, D-T]
    # Stem: x = 0, y in [-B/2, B/2], z in [-T, D-T]

    # Create a 3D grid and only keep vertices on the surface
    # This ensures shared vertices at edges
    vertex_map: dict[tuple[int, int, int], int] = {}

    def get_vertex(ix: int, iy: int, iz: int) -> int:
        """Get or create vertex at grid position (ix, iy, iz)."""
        key = (ix, iy, iz)
        if key in vertex_map:
            return vertex_map[key]
        x = L * ix / nx
        y = -B / 2 + B * iy / ny
        z = -T + (D) * iz / nz  # z from -T to D-T
        idx = builder.add_vertex((x, y, z))
        vertex_map[key] = idx
        return idx

    # Bottom face (iz=0): z = -T
    for ix in range(nx):
        for iy in range(ny):
            i0 = get_vertex(ix, iy, 0)
            i1 = get_vertex(ix + 1, iy, 0)
            i2 = get_vertex(ix + 1, iy + 1, 0)
            i3 = get_vertex(ix, iy + 1, 0)
            # Normal pointing down (-z): CCW when viewed from below
            builder.add_triangle(i0, i3, i2)
            builder.add_triangle(i0, i2, i1)

    # Top face (iz=nz): z = D-T
    for ix in range(nx):
        for iy in range(ny):
            i0 = get_vertex(ix, iy, nz)
            i1 = get_vertex(ix + 1, iy, nz)
            i2 = get_vertex(ix + 1, iy + 1, nz)
            i3 = get_vertex(ix, iy + 1, nz)
            # Normal pointing up (+z): CCW when viewed from above
            builder.add_triangle(i0, i1, i2)
            builder.add_triangle(i0, i2, i3)

    # Starboard face (iy=ny): y = B/2
    for ix in range(nx):
        for iz in range(nz):
            i0 = get_vertex(ix, ny, iz)
            i1 = get_vertex(ix + 1, ny, iz)
            i2 = get_vertex(ix + 1, ny, iz + 1)
            i3 = get_vertex(ix, ny, iz + 1)
            # Normal pointing +y: CCW when viewed from starboard
            builder.add_triangle(i0, i3, i2)
            builder.add_triangle(i0, i2, i1)

    # Port face (iy=0): y = -B/2
    for ix in range(nx):
        for iz in range(nz):
            i0 = get_vertex(ix, 0, iz)
            i1 = get_vertex(ix + 1, 0, iz)
            i2 = get_vertex(ix + 1, 0, iz + 1)
            i3 = get_vertex(ix, 0, iz + 1)
            # Normal pointing -y: CCW when viewed from port
            builder.add_triangle(i0, i1, i2)
            builder.add_triangle(i0, i2, i3)

    # Transom face (ix=nx): x = L
    for iy in range(ny):
        for iz in range(nz):
            i0 = get_vertex(nx, iy, iz)
            i1 = get_vertex(nx, iy + 1, iz)
            i2 = get_vertex(nx, iy + 1, iz + 1)
            i3 = get_vertex(nx, iy, iz + 1)
            # Normal pointing +x: CCW when viewed from aft
            builder.add_triangle(i0, i3, i2)
            builder.add_triangle(i0, i2, i1)

    # Stem face (ix=0): x = 0
    for iy in range(ny):
        for iz in range(nz):
            i0 = get_vertex(0, iy, iz)
            i1 = get_vertex(0, iy + 1, iz)
            i2 = get_vertex(0, iy + 1, iz + 1)
            i3 = get_vertex(0, iy, iz + 1)
            # Normal pointing -x: CCW when viewed from forward
            builder.add_triangle(i0, i1, i2)
            builder.add_triangle(i0, i2, i3)

    mesh = builder.build()
    oriented = meshlib.orient_outward(mesh)
    return oriented.mesh


def invert_normals(mesh: Mesh) -> Mesh:
    """Return a copy of *mesh* with all face windings reversed."""
    flipped = tuple((f[2], f[1], f[0]) for f in mesh.faces)
    return Mesh(
        vertices=mesh.vertices,
        faces=flipped,
        units=mesh.units,
        coordinate_frame=mesh.coordinate_frame,
    )


def generate_all(output_dir: Path) -> dict[str, Path]:
    """Generate all 8 STL assets and return {name: path}."""
    output_dir.mkdir(parents=True, exist_ok=True)
    assets: dict[str, Path] = {}

    # 1. Wigley (existing — copy from docs/wigley)
    src = Path(__file__).resolve().parent.parent / "docs" / "wigley" / "wigley.stl"
    if src.exists():
        dst = output_dir / "wigley.stl"
        shutil.copy2(src, dst)
        assets["wigley"] = dst
    else:
        wig_params = wigley.WigleyParameters(length=10.0, beam=1.0, draft=0.625, nx=41, nz=13)
        wig_mesh, _ = wigley.build_mesh(wig_params)
        assets["wigley"] = write_stl(wig_mesh, output_dir / "wigley.stl")

    # 2. Series 60 Cb=0.60
    s60_params = series60.Series60Parameters(
        length=20.0, beam=3.0, draft=1.5, cb=0.60, nx=41, nz=13,
    )
    s60_mesh, _offsets = series60.build_series60_mesh(s60_params)
    assets["series60_cb60"] = write_stl(s60_mesh, output_dir / "series60_cb60.stl")

    # 3. Series 60 Cb=0.75
    s60_params_75 = series60.Series60Parameters(
        length=20.0, beam=3.0, draft=1.5, cb=0.75, nx=41, nz=13,
    )
    s60_mesh_75, _ = series60.build_series60_mesh(s60_params_75)
    assets["series60_cb75"] = write_stl(s60_mesh_75, output_dir / "series60_cb75.stl")

    # 4. Box barge
    barge = make_box_barge()
    assets["box_barge"] = write_stl(barge, output_dir / "box_barge.stl")

    # 5. Sphere r=1
    sphere = make_sphere()
    assets["sphere_r1"] = write_stl(sphere, output_dir / "sphere_r1.stl")

    # 6. Inverted Wigley (red-team: inward normals)
    wig_params = wigley.WigleyParameters(length=10.0, beam=1.0, draft=0.625, nx=41, nz=13)
    wig_mesh, _ = wigley.build_mesh(wig_params)
    wig_inverted = invert_normals(wig_mesh)
    assets["inverted_wigley"] = write_stl(wig_inverted, output_dir / "inverted_wigley.stl")

    # 7. Wrong units Wigley (red-team: geometry in metres but labelled mm)
    # The diagonal of the Wigley is ~10 m; labelled as mm it becomes 0.01 m → scale violation
    wrong_units_mesh = Mesh(
        vertices=wig_mesh.vertices,
        faces=wig_mesh.faces,
        units="mm",  # deliberately wrong
        coordinate_frame=wig_mesh.coordinate_frame,
    )
    assets["wrong_units"] = write_stl(wrong_units_mesh, output_dir / "wrong_units.stl")

    # 8. Left-handed frame Wigley (red-team: wrong coordinate system declaration)
    left_handed_mesh = Mesh(
        vertices=wig_mesh.vertices,
        faces=wig_mesh.faces,
        units="m",
        coordinate_frame="x fore, y starboard, z down; left-handed",
    )
    assets["left_handed"] = write_stl(left_handed_mesh, output_dir / "left_handed.stl")

    return assets


if __name__ == "__main__":
    out = Path(__file__).resolve().parent.parent / "artifacts" / "geom_contract_assets"
    generated = generate_all(out)
    for name, path in generated.items():
        print(f"  {name}: {path} ({path.stat().st_size} bytes)")
    print(f"\n{len(generated)} assets generated in {out}")
