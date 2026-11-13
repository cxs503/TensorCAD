"""Negative controls: geometries and records that *must* be rejected.

A gate that has never seen a bad input is decoration. Each control here reproduces
one of the defects the owner reported on 2026-11-13, or is a plausible way to fool a
gate, and the run harness records the verdict the fleet gate returns for it:

* :func:`placeholder_box` — a box tessellated to exactly 304 triangles, the signature
  count of the six stand-in STLs in the ship-platform geometry library.
* :func:`duplicate_stl_pair` — two "different" ships written from one mesh, i.e. the
  KVLCC1 == KVLCC2 byte-identical defect.
* :func:`clipped_cb_record` — a ship whose ``Cb`` was raised into the Series-60 window
  so the route would accept it; the volume gate must still fail.
* :func:`wrong_volume_record` — a record whose declared displacement disagrees with
  its own particulars, which must fail the particulars gate before any geometry runs.
* :func:`mis_transcribed_beam` — the defect the round-2 cross-check actually found: a
  published displacement quoted next to a beam that was transcribed wrongly (DTC carried
  B = 45.8 m against a cited 173467 m^3, 11.3 % off its own Cb*L*B*T), which must fail
  ``G6b_reference_consistency`` before any geometry is built.
"""

from __future__ import annotations

import dataclasses
import struct
from pathlib import Path
from typing import Iterable, Sequence

from tensorcad.geometry import Mesh, MeshBuilder
from tensorcad.fleet.manifest import FleetManifest, Provenance, ShipRecord

PLACEHOLDER_SOURCE = "negative control: 304-triangle placeholder box (owner 2026-11-13 defect)"


def _box_quads(builder: MeshBuilder, length: float, beam: float, height: float) -> None:
    """Six faces of a box, tessellated so the triangle budget can be hit exactly."""
    levels_x, levels_z, half = 30, 1, 1

    def point(x: float, y: float, z: float) -> int:
        return builder.add_vertex((x, y, z))

    # four sides, each a 30 x 1 quad strip; 4 * 30 = 120 quads
    for side in range(4):
        def at(u: float, v: float) -> int:
            if side == 0:
                return point(u * length, -beam / 2, v * height)
            if side == 1:
                return point(length, u * beam, v * height)
            if side == 2:
                return point(length - u * length, beam / 2, v * height)
            return point(0.0, beam / 2 - u * beam, v * height)

        for i in range(levels_x):
            builder.add_quad(
                at(i / levels_x, 0.0), at((i + 1) / levels_x, 0.0),
                at((i + 1) / levels_x, 1.0), at(i / levels_x, 1.0),
            )
    # top and bottom: 4 x 4 grids each, 16 quads per face -> 32 quads
    for z_value in (0.0, height):
        def at_plane(u: float, v: float, z_value: float = z_value) -> int:
            return point(u * length, -beam / 2 + v * beam, z_value)

        for i in range(4):
            for j in range(4):
                builder.add_quad(
                    at_plane(i / 4, j / 4), at_plane((i + 1) / 4, j / 4),
                    at_plane((i + 1) / 4, (j + 1) / 4), at_plane(i / 4, (j + 1) / 4),
                )


def placeholder_box(triangles: int = 304, *, length: float = 100.0, beam: float = 20.0,
                    height: float = 8.0) -> Mesh:
    """A box tessellated to ``triangles`` faces (default: the 304-triangle placeholder)."""
    if triangles != 304:
        raise ValueError("this control reproduces the 304-triangle signature only")
    builder = MeshBuilder(units="m", coordinate_frame="box control")
    _box_quads(builder, length, beam, height)
    mesh = builder.build()
    if len(mesh.faces) != triangles:
        raise AssertionError(f"placeholder control produced {len(mesh.faces)} triangles, expected {triangles}")
    return mesh


def write_binary_stl(mesh: Mesh, path: str | Path, *, solid_name: str = "placeholder") -> Path:
    """Write a binary STL, so the duplicate-file control is a real file-level defect."""
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    header = solid_name.encode("ascii", "replace")[:80].ljust(80, b" ")
    payload = bytearray(header)
    payload += struct.pack("<I", len(mesh.faces))
    for face in mesh.faces:
        payload += struct.pack("<3f", 0.0, 0.0, 0.0)
        for index in face:
            payload += struct.pack("<3f", *mesh.vertices[index])
        payload += struct.pack("<H", 0)
    target.write_bytes(bytes(payload))
    return target


def duplicate_stl_pair(directory: str | Path, *, triangles: int = 304) -> tuple[Path, Path]:
    """Two STL files with different names written from one mesh: byte-identical bytes."""
    root = Path(directory)
    mesh = placeholder_box(triangles)
    first = write_binary_stl(mesh, root / "KVLCC1.stl", solid_name="KVLCC1")
    second = write_binary_stl(mesh, root / "KVLCC2.stl", solid_name="KVLCC2")
    # the header differs by name; blank it so the files are byte-identical, as reported
    for path in (first, second):
        payload = bytearray(path.read_bytes())
        payload[:80] = b" " * 80
        path.write_bytes(bytes(payload))
    return first, second


def _clone(record: ShipRecord, **changes: object) -> ShipRecord:
    return dataclasses.replace(record, **changes)


def clipped_cb_record(base: ShipRecord, cb: float = 0.55) -> ShipRecord:
    """The same ship with ``Cb`` pushed into the form family's window but the volume kept."""
    return _clone(
        base,
        name=f"{base.name}_clipped_cb",
        block_coefficient=cb,
        notes=base.notes + " [negative control: Cb clipped into the Series-60 window]",
        provenance=Provenance(
            source="negative control derived from " + base.provenance.source,
            url=base.provenance.url,
            license=base.provenance.license,
            retrieval="derived_control",
            verified=False,
        ),
    )


def wrong_volume_record(base: ShipRecord, factor: float = 1.08) -> ShipRecord:
    """A record whose declared displacement contradicts its own ``Cb*L*B*T``."""
    volume = base.reference_volume() * factor
    return _clone(
        base,
        name=f"{base.name}_wrong_volume",
        displacement_volume=volume,
        displacement_basis="literature_value_recalled",
        notes=base.notes + f" [negative control: volume scaled by {factor}]",
        provenance=Provenance(
            source="negative control derived from " + base.provenance.source,
            url=base.provenance.url,
            license=base.provenance.license,
            retrieval="derived_control",
            verified=False,
        ),
    )


def duplicate_manifest(base: FleetManifest, first: Path, second: Path) -> FleetManifest:
    """A two-ship manifest whose geometry hashes collide, i.e. the reported defect."""
    from tensorcad.fleet.manifest import sha256_file

    sha = sha256_file(first)
    template = base.ships[0]
    ships = (
        _clone(template, name="KVLCC1", geometry_path=str(first), geometry_sha256=sha,
               geometry_route="imported_stl"),
        _clone(template, name="KVLCC2", geometry_path=str(second), geometry_sha256=sha,
               geometry_route="imported_stl"),
    )
    return FleetManifest(ships=ships, placeholder=dict(base.placeholder), path=base.path, sha256=base.sha256)

def mis_transcribed_beam(base: ShipRecord, *, beam: float = 45.8) -> ShipRecord:
    """The same record with one particular transcribed wrongly, volume left as cited.

    This is not a synthetic trap: it is the DTC entry as it stood before the round-2
    cross-check (B = 45.8 m next to a cited 173467 m^3), and it is the reason the gate
    exists. A record like this passed every pre-existing gate because it had no cited
    displacement and its reference was the closure, i.e. it was checked against itself.
    """
    return _clone(
        base,
        name=f"{base.name}_mis_transcribed",
        beam=beam,
        notes=base.notes + f" [negative control: beam mis-transcribed to {beam} m, volume kept as cited]",
        provenance=Provenance(
            source="negative control derived from " + base.provenance.source,
            url=base.provenance.url,
            license=base.provenance.license,
            retrieval="derived_control",
            verified=False,
        ),
    )


#: The six stand-in files and the byte-identical pair, as reported by the owner on
#: 2026-11-13. The partition is not reconstructable offline (the library itself is not in
#: this repository), so the fixture reproduces the two *defect classes* at their reported
#: sizes: six 304-triangle placeholders, and two files that share one sha256.
DEFECT_PLACEHOLDERS = ("KCS", "KVLCC2", "ONRT", "JBC", "DTC", "Series60")
DEFECT_DUPLICATE_PAIR = ("KVLCC1", "DTMB5415")


def defect_geometry_library(
    directory: str | Path,
    *,
    placeholders: Sequence[str] = DEFECT_PLACEHOLDERS,
    duplicate_pair: Sequence[str] = DEFECT_DUPLICATE_PAIR,
    duplicate_mesh=None,
) -> dict[str, list[Path]]:
    """Write an 8-file library carrying both reported defects; return the file paths.

    ``duplicate_mesh`` lets a caller supply the hull written twice (the default builds one
    with the Series-60 route at a coarse grid); both files must be *byte* identical, so the
    same mesh and the same header are used and the header is then blanked, as an exporter
    that ignores the solid name would produce.
    """
    root = Path(directory)
    root.mkdir(parents=True, exist_ok=True)
    box = placeholder_box(304)
    written: dict[str, list[Path]] = {"placeholders": [], "duplicates": []}
    for name in placeholders:
        written["placeholders"].append(write_binary_stl(box, root / f"{name}.stl", solid_name=name))
    if duplicate_mesh is None:
        from tensorcad.fleet.build import build_ship_geometry
        from tensorcad.fleet.manifest import load_manifest

        record = next(
            ship for ship in load_manifest().ships if ship.name == "KCS"
        )
        duplicate_mesh = build_ship_geometry(record, nx=81, nz=25).mesh
    for name in duplicate_pair:
        path = write_binary_stl(duplicate_mesh, root / f"{name}.stl", solid_name=name)
        payload = bytearray(path.read_bytes())
        payload[:80] = b" " * 80
        path.write_bytes(bytes(payload))
        written["duplicates"].append(path)
    return written
