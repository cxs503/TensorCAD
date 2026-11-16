"""Third-party mesh tool oracle (trimesh): execute a real foreign reader.

Why this module exists
======================
The independent readers in :mod:`tensorcad.mesh_interop_independent` and
:mod:`tensorcad.mesh_foreign_dialects` are written here, from the format
specifications.  That is much stronger than "meshio writes it and meshio reads
it back", but it is still *our* code judging *our* bytes: a systematic
misreading of a specification would pass the table and be invisible.

``docs/mesh_format_matrix_thirdparty/verification_report.md`` states exactly
this ceiling:

    No external CAD/solver tool was executed: gmsh, Salome, ParaView are not
    installed in this worktree, so "another tool can read it" is evidenced by
    spec-conformance of the bytes plus independent parsing, not by running that
    tool.  That is the honest ceiling of this verification.

``trimesh`` *is* an independent, third-party implementation (mikedh/trimesh,
5.x on this box, its own STL/OBJ loaders and its own mass-properties code) and
it is installed.  This module uses it in the two directions that matter:

* **consumer** -- files written by TensorCAD are loaded by trimesh and trimesh's
  own volume / watertightness / winding-consistency answers are compared with
  the reference mesh (:func:`trimesh_probe`);
* **producer** -- files written by *trimesh* (binary STL, OBJ) are read by the
  independent readers, so the dialect parser is validated against bytes it did
  not produce (:func:`trimesh_write`).

Coverage is honestly partial: trimesh reads STL and OBJ, not VTK/MSH/SU2/CGNS,
so there is no third-party oracle for four of the six formats.  Absence of a
check is reported as ``NOT_APPLICABLE``/``NOT_RUN`` -- never as a pass.

The module imports trimesh lazily and never imports meshio: it must not be able
to share a convention with the writer under test.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Sequence

import numpy as np

__all__ = [
    "TRIMESH_FORMATS",
    "trimesh_available",
    "trimesh_version",
    "trimesh_probe",
    "trimesh_write",
]

#: formats for which trimesh provides a real third-party reader
TRIMESH_FORMATS = ("stl", "obj")


def trimesh_available() -> bool:
    """True when trimesh can be imported in this environment."""
    try:
        import trimesh  # noqa: F401
    except Exception:  # pragma: no cover - environment dependent
        return False
    return True


def trimesh_version() -> str | None:
    try:
        import trimesh
    except Exception:  # pragma: no cover - environment dependent
        return None
    return getattr(trimesh, "__version__", "unknown")


def _require():
    try:
        import trimesh
    except Exception as exc:  # pragma: no cover - environment dependent
        raise RuntimeError("trimesh is required for the third-party oracle") from exc
    return trimesh


def trimesh_probe(path: str | Path) -> dict[str, Any]:
    """Load ``path`` with trimesh and report its own geometric verdicts.

    The returned volume is trimesh's divergence-theorem volume over *its* merged
    triangle soup -- computed by code that shares nothing with
    ``tensorcad.mesh_interop_independent.surface_signed_volume``.
    """
    trimesh = _require()
    p = Path(path)
    loaded = trimesh.load(str(p), force="mesh", process=True)
    mesh = loaded if hasattr(loaded, "faces") else trimesh.util.concatenate(
        list(loaded.geometry.values())
    )
    return {
        "tool": "trimesh",
        "tool_version": trimesh_version(),
        "file": p.name,
        "volume": float(mesh.volume),
        "watertight": bool(mesh.is_watertight),
        "winding_consistent": bool(mesh.is_winding_consistent),
        "faces": int(len(mesh.faces)),
        "vertices": int(len(mesh.vertices)),
        "euler_number": int(mesh.euler_number),
        "area": float(mesh.area),
    }


def trimesh_write(
    path: str | Path,
    vertices: np.ndarray,
    faces: Sequence[Sequence[int]],
    *,
    file_type: str,
) -> Path:
    """Write a mesh with trimesh's own exporter (binary STL by default for STL)."""
    trimesh = _require()
    if file_type not in TRIMESH_FORMATS:
        raise ValueError(f"trimesh oracle covers {TRIMESH_FORMATS}, not {file_type!r}")
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    mesh = trimesh.Trimesh(
        vertices=np.asarray(vertices, dtype=np.float64),
        faces=np.asarray(faces, dtype=np.int64),
        process=True,
    )
    payload = mesh.export(file_type=file_type)
    if isinstance(payload, str):
        p.write_text(payload)
    else:
        p.write_bytes(payload)
    return p