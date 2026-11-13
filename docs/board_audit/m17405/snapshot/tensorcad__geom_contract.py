"""Geometry contract: four conventions every TensorCAD mesh must honour.

The contract is deliberately small and machine-checkable.  A mesh that passes
all four checks is *contract-compliant*; downstream consumers (solvers, viewers,
export pipelines) can rely on the declared units, frame, scale and orientation
without guessing.

Four conventions
----------------
1. **Coordinate system** — right-handed, Cartesian.  The cross product of the
   first two basis vectors must point in the direction of the third.
2. **Units** — one of the SI length tokens ``m``, ``mm``, ``cm``, ``dm``, ``km``.
   The mesh declares its units in ``Mesh.units``; the contract requires the
   declaration to be present and valid.
3. **Scale** — the bounding-box diagonal must lie in a plausible range for the
   declared units.  A "ship" in metres should not have a 0.001 m diagonal
   (that is a model in millimetres mislabelled as metres).
4. **Orientation** — for a closed surface the signed volume must be positive
   (outward-facing normals).  For an open surface the majority of face normals
   must point into the same hemisphere (consistent winding).

The validator can be used as a library (:func:`validate_contract`) or as a CLI
(``python -m tensorcad.geom_contract <stl-file>``).
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Sequence

from . import mesh as meshlib
from .export import read_stl
from .geometry import Mesh

__all__ = [
    "ContractViolation",
    "ContractReport",
    "VALID_UNITS",
    "MIN_DIAGONAL_M",
    "MAX_DIAGONAL_M",
    "validate_contract",
    "validate_stl_file",
    "main",
]

#: Allowed SI length unit tokens.
VALID_UNITS: tuple[str, ...] = ("m", "mm", "cm", "dm", "km")

#: Plausible bounding-box diagonal range **in metres** for a ship-scale model.
#: Models whose diagonal falls outside this band are flagged as a scale
#: mismatch (e.g. a millimetre model mislabelled as metres).
MIN_DIAGONAL_M: float = 0.05   # 1 cm — smallest reasonable test model
MAX_DIAGONAL_M: float = 500.0  # 500 m — larger than any real ship


# --------------------------------------------------------------------------- #
# Violation and report types
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class ContractViolation:
    """One contract breach."""

    convention: str  # "coordinate_system" | "units" | "scale" | "orientation"
    severity: str    # "error" | "warning"
    message: str
    detail: str = ""

    def to_dict(self) -> dict[str, str]:
        d: dict[str, str] = {
            "convention": self.convention,
            "severity": self.severity,
            "message": self.message,
        }
        if self.detail:
            d["detail"] = self.detail
        return d


@dataclass(frozen=True)
class ContractReport:
    """Full contract-validation result for one mesh."""

    source: str
    passed: bool
    violations: tuple[ContractViolation, ...] = ()
    mesh_info: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {
            "source": self.source,
            "passed": self.passed,
            "violations": [v.to_dict() for v in self.violations],
            "mesh_info": self.mesh_info,
        }


# --------------------------------------------------------------------------- #
# Unit conversion helpers
# --------------------------------------------------------------------------- #
_UNIT_TO_METRE: dict[str, float] = {
    "m": 1.0,
    "mm": 1e-3,
    "cm": 1e-2,
    "dm": 1e-1,
    "km": 1e3,
}


def _diagonal_in_metres(mesh: Mesh) -> float:
    """Bounding-box diagonal expressed in metres."""
    bb = mesh.bounding_box()
    sx, sy, sz = bb.sizes
    factor = _UNIT_TO_METRE.get(mesh.units, 1.0)
    return math.sqrt(sx * sx + sy * sy + sz * sz) * factor


# --------------------------------------------------------------------------- #
# Individual convention checks
# --------------------------------------------------------------------------- #
def _check_coordinate_system(mesh: Mesh) -> list[ContractViolation]:
    """Verify the declared frame is right-handed.

    We cannot prove handedness from vertex positions alone, so we check the
    declaration string for the keyword ``right-handed`` (case-insensitive).
    If the frame is ``unspecified`` we emit a warning.
    """
    frame = mesh.coordinate_frame.lower()
    if "right-handed" in frame:
        return []
    if "left-handed" in frame:
        return [ContractViolation(
            convention="coordinate_system",
            severity="error",
            message="left-handed coordinate system is not supported",
            detail=f"coordinate_frame={mesh.coordinate_frame!r}",
        )]
    return [ContractViolation(
        convention="coordinate_system",
        severity="warning",
        message="coordinate frame is unspecified; right-handed is assumed",
        detail=f"coordinate_frame={mesh.coordinate_frame!r}",
    )]


def _check_units(mesh: Mesh) -> list[ContractViolation]:
    """Verify the unit declaration is a recognised SI length token."""
    if mesh.units in VALID_UNITS:
        return []
    return [ContractViolation(
        convention="units",
        severity="error",
        message=f"unrecognised unit {mesh.units!r}",
        detail=f"allowed units: {VALID_UNITS}",
    )]


def _check_scale(mesh: Mesh) -> list[ContractViolation]:
    """Verify the bounding-box diagonal is plausible for the declared units."""
    if not mesh.vertices:
        return [ContractViolation(
            convention="scale",
            severity="error",
            message="empty mesh has no scale",
        )]
    diag = _diagonal_in_metres(mesh)
    violations: list[ContractViolation] = []
    if diag < MIN_DIAGONAL_M:
        violations.append(ContractViolation(
            convention="scale",
            severity="error",
            message=f"bounding-box diagonal {diag:.6g} m is below minimum {MIN_DIAGONAL_M} m",
            detail=f"declared units={mesh.units!r}; the model may be in a different unit",
        ))
    elif diag > MAX_DIAGONAL_M:
        violations.append(ContractViolation(
            convention="scale",
            severity="error",
            message=f"bounding-box diagonal {diag:.6g} m exceeds maximum {MAX_DIAGONAL_M} m",
            detail=f"declared units={mesh.units!r}; the model may be in a different unit",
        ))
    return violations


def _check_orientation(mesh: Mesh) -> list[ContractViolation]:
    """Verify outward-facing normals for closed surfaces, consistent winding for open ones."""
    if not mesh.faces:
        return [ContractViolation(
            convention="orientation",
            severity="error",
            message="mesh has no faces",
        )]
    vol = meshlib.total_signed_volume(mesh)
    if meshlib.is_watertight(mesh):
        if vol <= 0.0:
            return [ContractViolation(
                convention="orientation",
                severity="error",
                message="closed mesh has non-positive signed volume (inward normals or inverted)",
                detail=f"signed_volume={vol:.9g}",
            )]
        return []
    # Open surface: check that the majority of normals point into the same hemisphere.
    normals = [meshlib.face_normal(mesh, f) for f in mesh.faces]
    # Use the first normal as the reference direction.
    ref = normals[0]
    ref_len = math.sqrt(ref[0] ** 2 + ref[1] ** 2 + ref[2] ** 2)
    if ref_len == 0.0:
        return [ContractViolation(
            convention="orientation",
            severity="error",
            message="first face has zero-area normal",
        )]
    ref = (ref[0] / ref_len, ref[1] / ref_len, ref[2] / ref_len)
    agreeing = 0
    for n in normals:
        n_len = math.sqrt(n[0] ** 2 + n[1] ** 2 + n[2] ** 2)
        if n_len == 0.0:
            continue
        dot = (n[0] * ref[0] + n[1] * ref[1] + n[2] * ref[2]) / n_len
        if dot >= 0.0:
            agreeing += 1
    ratio = agreeing / len(normals)
    if ratio < 0.5:
        return [ContractViolation(
            convention="orientation",
            severity="warning",
            message=f"only {ratio:.1%} of face normals agree with the majority direction",
            detail="open surface has inconsistent winding",
        )]
    return []


# --------------------------------------------------------------------------- #
# Public API
# --------------------------------------------------------------------------- #
def validate_contract(mesh: Mesh, *, source: str = "<unknown>") -> ContractReport:
    """Run all four convention checks on *mesh* and return a report."""
    violations: list[ContractViolation] = []
    violations.extend(_check_coordinate_system(mesh))
    violations.extend(_check_units(mesh))
    violations.extend(_check_scale(mesh))
    violations.extend(_check_orientation(mesh))
    errors = [v for v in violations if v.severity == "error"]
    info: dict[str, Any] = {
        "units": mesh.units,
        "coordinate_frame": mesh.coordinate_frame,
        "vertex_count": len(mesh.vertices),
        "triangle_count": len(mesh.faces),
        "watertight": meshlib.is_watertight(mesh),
    }
    if mesh.vertices:
        bb = mesh.bounding_box()
        info["bounding_box"] = bb.to_dict()
        info["diagonal_m"] = round(_diagonal_in_metres(mesh), 9)
        info["signed_volume"] = round(meshlib.total_signed_volume(mesh), 9)
    return ContractReport(
        source=source,
        passed=len(errors) == 0,
        violations=tuple(violations),
        mesh_info=info,
    )


def validate_stl_file(
    path: str | Path,
    *,
    units: str = "m",
    coordinate_frame: str = "unspecified right-handed frame",
) -> ContractReport:
    """Read a binary STL and validate it against the geometry contract."""
    mesh = read_stl(path, units=units, coordinate_frame=coordinate_frame)
    return validate_contract(mesh, source=str(path))


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Validate a binary STL against the TensorCAD geometry contract.",
    )
    parser.add_argument("stl", nargs="+", help="binary STL file(s) to validate")
    parser.add_argument(
        "--units",
        default="m",
        help="declared length unit (default: m)",
    )
    parser.add_argument(
        "--frame",
        default="unspecified right-handed frame",
        help="declared coordinate frame (default: 'unspecified right-handed frame')",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="output JSON instead of human-readable text",
    )
    args = parser.parse_args(argv)

    reports: list[ContractReport] = []
    for stl_path in args.stl:
        try:
            report = validate_stl_file(stl_path, units=args.units, coordinate_frame=args.frame)
        except Exception as exc:
            report = ContractReport(
                source=stl_path,
                passed=False,
                violations=(ContractViolation(
                    convention="read",
                    severity="error",
                    message=f"failed to read STL: {exc}",
                ),),
            )
        reports.append(report)

    if args.json:
        payload = [r.to_dict() for r in reports]
        print(json.dumps(payload, indent=2))
    else:
        for r in reports:
            status = "PASS" if r.passed else "FAIL"
            print(f"[{status}] {r.source}")
            for v in r.violations:
                tag = "!!" if v.severity == "error" else "!?"
                print(f"  {tag} [{v.convention}] {v.message}")
                if v.detail:
                    print(f"     {v.detail}")
            if r.mesh_info:
                print(f"  vertices={r.mesh_info.get('vertex_count', '?')}  "
                      f"triangles={r.mesh_info.get('triangle_count', '?')}  "
                      f"watertight={r.mesh_info.get('watertight', '?')}")
                if "diagonal_m" in r.mesh_info:
                    print(f"  diagonal={r.mesh_info['diagonal_m']:.6g} m  "
                          f"volume={r.mesh_info.get('signed_volume', '?')}")
            print()

    # Exit non-zero if any report has errors
    any_fail = any(not r.passed for r in reports)
    return 1 if any_fail else 0


if __name__ == "__main__":
    sys.exit(main())
