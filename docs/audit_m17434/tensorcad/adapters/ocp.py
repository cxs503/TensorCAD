"""Optional OpenCASCADE (OCP) bridge: turn a validated mesh into a STEP file.

The core package never imports this module and never needs it: STEP is a B-rep format,
the core is a field protocol, and the seam between them is a triangle mesh plus the
handoff record that says what the mesh is. This module is that seam, and it is
deliberately one-directional and dumb:

    mesh  ->  sewn shell  ->  solid  ->  STEP

What it does **not** do is invent a B-rep. Sewing a triangle soup produces a faceted
solid -- thousands of planar faces, one per triangle -- which is exactly what a
solver-facing exchange file needs (a closed, valid solid of known volume) and exactly
what a designer does *not* want (no fillets, no analytic surfaces, no feature tree).
Both statements are in :mod:`tensorcad.adapters` and in the docs; conflating them is
how "we export STEP" turns into a claim nobody can check.

Everything is measured rather than assumed. :func:`write_step` returns the sewn shell's
free-edge count, whether the shell closed, whether OCP's own validity checker accepted
the solid, and the solid's volume; :func:`read_step` re-imports the file and reports the
volume it finds. A caller that wants the round trip is one call away
(:attr:`StepExport.round_trip_volume_error`).

Units: STEP is conventionally millimetres and OCP writes coordinates as given, so a
model in metres is scaled by ``length_scale`` (default 1000) on the way out and by the
same factor on the way back. The scale is recorded in the result, never assumed.

Requires the ``ocp`` extra (``pip install "tensorcad[ocp]"``).
"""
from __future__ import annotations

import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from ..geometry import Mesh

__all__ = [
    "OcpUnavailable",
    "StepExport",
    "StepImport",
    "have_ocp",
    "read_step",
    "write_step",
]


class OcpUnavailable(RuntimeError):
    """Raised when OCP is needed but not installed."""


def have_ocp() -> bool:
    """True when the optional OCP (OpenCASCADE) wheels are importable."""
    try:  # pragma: no cover - depends on the extra being installed
        import OCP.STEPControl  # noqa: F401
    except Exception:  # pragma: no cover - any import failure means "not usable"
        return False
    return True


def _ocp() -> Any:
    try:
        from OCP.BRepBuilderAPI import (
            BRepBuilderAPI_MakeFace,
            BRepBuilderAPI_MakePolygon,
            BRepBuilderAPI_MakeSolid,
            BRepBuilderAPI_Sewing,
        )
        from OCP.BRepCheck import BRepCheck_Analyzer
        from OCP.BRepGProp import BRepGProp
        from OCP.gp import gp_Pnt
        from OCP.GProp import GProp_GProps
        from OCP.IFSelect import IFSelect_RetDone
        from OCP.STEPControl import STEPControl_AsIs, STEPControl_Reader, STEPControl_Writer
        from OCP.TopoDS import TopoDS
    except Exception as exc:  # pragma: no cover - exercised only without the extra
        raise OcpUnavailable(
            "this call needs the optional OCP (OpenCASCADE) wheels: "
            'pip install "tensorcad[ocp]"'
        ) from exc
    return {
        "gp_Pnt": gp_Pnt,
        "MakePolygon": BRepBuilderAPI_MakePolygon,
        "MakeFace": BRepBuilderAPI_MakeFace,
        "MakeSolid": BRepBuilderAPI_MakeSolid,
        "Sewing": BRepBuilderAPI_Sewing,
        "Analyzer": BRepCheck_Analyzer,
        "VolumeProperties_s": BRepGProp.VolumeProperties_s,
        "SurfaceProperties_s": BRepGProp.SurfaceProperties_s,
        "GProps": GProp_GProps,
        "RetDone": IFSelect_RetDone,
        "AsIs": STEPControl_AsIs,
        "Writer": STEPControl_Writer,
        "Reader": STEPControl_Reader,
        "TopoDS": TopoDS,
    }


def _measure(ocp: dict[str, Any], shape: Any, properties: str) -> float:
    """Integrate a shape property, tolerating both OCP call conventions.

    ``BRepGProp.VolumeProperties_s(shape, props)`` fills ``props`` in place and returns
    ``None`` in the wheels this was measured against; some OCP builds return the
    properties instead. Reading ``props`` afterwards works either way.
    """
    props = ocp["GProps"]()
    returned = ocp[properties](shape, props)
    result = returned if returned is not None else props
    return float(result.Mass())


def _shells(ocp: dict[str, Any], shape: Any) -> list[Any]:
    """Every shell in a sewed shape, whether it sewed into one shell or a compound."""
    from OCP.TopAbs import TopAbs_SHELL
    from OCP.TopExp import TopExp_Explorer

    shells: list[Any] = []
    explorer = TopExp_Explorer(shape, TopAbs_SHELL)
    while explorer.More():
        shells.append(ocp["TopoDS"].Shell(explorer.Current()))
        explorer.Next()
    return shells


def _compound(ocp: dict[str, Any], solids: list[Any]) -> Any:
    """A compound of solids, for the case where sewing left more than one shell.

    Sewing never merges two closed shells that touch at a point or run through the same
    vertex, so a mesh with a detached island (a marching-tet artefact at coarse
    resolution, never at fine resolution) sews into one shell per island. Writing them
    as a compound keeps every island in the file and counts them, instead of silently
    dropping all but the largest -- which is exactly the bug a STEP export must not
    have, because the volume in the file would then disagree with the mesh.
    """
    from OCP.BRep import BRep_Builder
    from OCP.TopoDS import TopoDS_Compound

    builder = BRep_Builder()
    compound = TopoDS_Compound()
    builder.MakeCompound(compound)
    for solid in solids:
        builder.Add(compound, solid)
    return compound


@dataclass(frozen=True)
class StepExport:
    """What was written, and what OCP measured on the way out and back."""

    path: str
    faces: int
    free_edges: int
    multiple_edges: int
    shells: int
    solids: int
    shell_closed: bool
    solid_valid: bool
    volume_m3: float
    area_m2: float
    length_scale: float
    file_bytes: int
    seconds_sew: float
    seconds_write: float
    round_trip_volume_m3: float | None = None

    @property
    def round_trip_volume_error(self) -> float | None:
        """Relative volume error on re-import, or ``None`` if no round trip was done."""
        if self.round_trip_volume_m3 is None or self.volume_m3 == 0.0:
            return None
        return self.round_trip_volume_m3 / self.volume_m3 - 1.0

    @property
    def passed(self) -> bool:
        return (
            self.free_edges == 0
            and self.multiple_edges == 0
            and self.shells == self.solids
            and self.shell_closed
            and self.solid_valid
            and self.volume_m3 > 0.0
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "path": self.path,
            "faces": self.faces,
            "free_edges": self.free_edges,
            "multiple_edges": self.multiple_edges,
            "shells": self.shells,
            "solids": self.solids,
            "shell_closed": self.shell_closed,
            "solid_valid": self.solid_valid,
            "volume_m3": self.volume_m3,
            "area_m2": self.area_m2,
            "length_scale": self.length_scale,
            "file_bytes": self.file_bytes,
            "seconds_sew": self.seconds_sew,
            "seconds_write": self.seconds_write,
            "round_trip_volume_m3": self.round_trip_volume_m3,
            "round_trip_volume_error": self.round_trip_volume_error,
            "passed": self.passed,
        }


@dataclass(frozen=True)
class StepImport:
    """What a STEP file turned out to contain."""

    path: str
    volume_m3: float
    area_m2: float
    faces: int
    solids: int
    valid: bool
    length_scale: float

    def to_dict(self) -> dict[str, Any]:
        return {
            "path": self.path,
            "volume_m3": self.volume_m3,
            "area_m2": self.area_m2,
            "faces": self.faces,
            "solids": self.solids,
            "valid": self.valid,
            "length_scale": self.length_scale,
        }


def write_step(
    mesh: Mesh,
    path: str | Path,
    *,
    length_scale: float = 1000.0,
    sew_tolerance: float | None = None,
    round_trip: bool = False,
) -> StepExport:
    """Sew ``mesh`` into a solid and write it as STEP, reporting what happened.

    ``length_scale`` multiplies the coordinates on the way out: STEP is conventionally
    in millimetres and the models here are in metres, so the default is 1000.
    ``sew_tolerance`` defaults to ``1e-9 * length_scale``, i.e. a nanometre in model
    units and far below any cell size in this repository, so it can only join vertices
    that are meant to be identical.

    ``round_trip=True`` re-imports the file with :func:`read_step` and reports the
    volume it finds; that is the only honest way to say "the file is readable".
    """
    if not mesh.faces:
        raise ValueError("cannot write STEP for an empty mesh")
    if length_scale <= 0.0:
        raise ValueError(f"length_scale must be positive, got {length_scale!r}")
    ocp = _ocp()
    point = ocp["gp_Pnt"]
    tolerance = 1e-9 * length_scale if sew_tolerance is None else float(sew_tolerance)
    if tolerance <= 0.0:
        raise ValueError(f"sew_tolerance must be positive, got {sew_tolerance!r}")

    vertices = [(length_scale * x, length_scale * y, length_scale * z) for x, y, z in mesh.vertices]
    started = time.perf_counter()
    sewing = ocp["Sewing"](tolerance)
    for a, b, c in mesh.faces:
        polygon = ocp["MakePolygon"](point(*vertices[a]), point(*vertices[b]), point(*vertices[c]), True)
        face = ocp["MakeFace"](polygon.Wire())
        sewing.Add(face.Face())
    sewing.Perform()
    seconds_sew = time.perf_counter() - started
    shells = _shells(ocp, sewing.SewedShape())
    if not shells:
        raise RuntimeError("sewing produced no shell; the mesh has no closed surface")
    solids = [ocp["MakeSolid"](shell).Solid() for shell in shells]
    shape = solids[0] if len(solids) == 1 else _compound(ocp, solids)
    volume = sum(_measure(ocp, solid, "VolumeProperties_s") for solid in solids) / length_scale**3
    area = sum(_measure(ocp, solid, "SurfaceProperties_s") for solid in solids) / length_scale**2

    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    writer = ocp["Writer"]()
    writer.Transfer(shape, ocp["AsIs"])
    status = writer.Write(str(target))
    seconds_write = time.perf_counter() - started
    if status != ocp["RetDone"]:
        raise RuntimeError(f"STEP writer refused {target}: status {status!r}")

    record = StepExport(
        path=str(target),
        faces=len(mesh.faces),
        free_edges=int(sewing.NbFreeEdges()),
        multiple_edges=int(sewing.NbMultipleEdges()),
        shells=len(shells),
        solids=len(solids),
        shell_closed=all(bool(shell.Closed()) for shell in shells),
        solid_valid=all(bool(ocp["Analyzer"](solid).IsValid()) for solid in solids),
        volume_m3=float(volume),
        area_m2=float(area),
        length_scale=length_scale,
        file_bytes=target.stat().st_size,
        seconds_sew=seconds_sew,
        seconds_write=seconds_write,
    )
    if not round_trip:
        return record
    imported = read_step(target, length_scale=length_scale)
    return StepExport(
        **{**record.__dict__, "round_trip_volume_m3": imported.volume_m3}
    )


def read_step(path: str | Path, *, length_scale: float = 1000.0) -> StepImport:
    """Read a STEP file and measure it, so a round trip can be stated in numbers."""
    ocp = _ocp()
    target = Path(path)
    if not target.exists():
        raise FileNotFoundError(f"no STEP file at {target}")
    reader = ocp["Reader"]()
    status = reader.ReadFile(str(target))
    if status != ocp["RetDone"]:
        raise RuntimeError(f"STEP reader refused {target}: status {status!r}")
    reader.TransferRoots()
    shape = reader.OneShape()
    volume = _measure(ocp, shape, "VolumeProperties_s") / length_scale**3
    area = _measure(ocp, shape, "SurfaceProperties_s") / length_scale**2

    from OCP.TopAbs import TopAbs_FACE, TopAbs_SOLID
    from OCP.TopExp import TopExp_Explorer

    def count(kind: Any) -> int:
        explorer = TopExp_Explorer(shape, kind)
        found = 0
        while explorer.More():
            found += 1
            explorer.Next()
        return found

    return StepImport(
        path=str(target),
        volume_m3=float(volume),
        area_m2=float(area),
        faces=count(TopAbs_FACE),
        solids=count(TopAbs_SOLID),
        valid=bool(ocp["Analyzer"](shape).IsValid()),
        length_scale=length_scale,
    )
