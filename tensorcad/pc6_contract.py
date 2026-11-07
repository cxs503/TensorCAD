"""The PC6 geometry contract: one JSON document the whole campaign consumes.

``artifacts/geom_pc6.json`` is the single source of geometric truth for the ice
legs. It carries both representations of the hull -- a binary STL triangulation
and an implicit field (analytic definition plus a coarse sample grid) -- together
with the parameter table, the topology and hydrostatic evidence, the ice channel
of the de Long -> Pevek route, the beam sweep 24 -> 28 m, the negative controls
that were run, and a content digest.

A consumer never reads a value back by hand::

    from tensorcad.pc6_contract import load_contract

    geometry = load_contract("artifacts/geom_pc6.json")
    geometry.hull                      # tensorcad.geometry.Mesh, already re-checked
    geometry.field(x, y, z)            # implicit field, negative inside
    geometry.channel                   # closed corridor solid of the ice channel

:func:`load_contract` re-runs the watertightness, topology and displacement
checks on the file it just read, so a truncated or hand-edited artifact fails
loudly instead of poisoning three downstream solvers.
"""
from __future__ import annotations

import hashlib
import json
import math
import re
from collections.abc import Sequence
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any

from .export import read_stl, write_metadata, write_stl
from .geometry import Mesh
from .hydrostatics import compute_hydrostatics
from .icechannel import (
    DEFAULT_DEPTH_M,
    DEFAULT_WIDTH_M,
    IceChannel,
    build_channel_mesh,
    channel_measurements,
    ice_channel,
)
from .intersect import find_self_intersections
from .pc6 import (
    COORDINATE_FRAME,
    UNITS,
    PC6Model,
    PC6Parameters,
    beam_sweep,
    build_mesh,
    field_value,
    generate,
    inside_by_winding_number,
    sampled_field,
    topology,
)

CONTRACT_SCHEMA = "tensorcad.pc6.geom/1.0"
CONTRACT_NAME = "geom_pc6.json"
HULL_STL = "pc6_hull.stl"
CHANNEL_STL = "pc6_channel.stl"
SWEEP_DIR = "pc6_sweep"

#: Deterministic pseudo-random points for the field/mesh agreement check: a fixed
#: linear congruential sequence, so the artifact is reproducible bit-for-bit.
_SAMPLE_SEED = 20261106


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


def _clean(value: Any) -> Any:
    """JSON-safe copy: non-finite floats become ``null`` rather than invalid JSON."""
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, dict):
        return {str(key): _clean(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_clean(item) for item in value]
    return value


def _lcg_points(count: int, bounds: dict[str, tuple[float, float]]) -> list[tuple[float, float, float]]:
    state = _SAMPLE_SEED
    points: list[tuple[float, float, float]] = []
    for _ in range(count):
        sample: list[float] = []
        for axis in ("x", "y", "z"):
            state = (1103515245 * state + 12345) % (2 ** 31)
            low, high = bounds[axis]
            sample.append(low + (high - low) * state / (2 ** 31))
        points.append((sample[0], sample[1], sample[2]))
    return points


def field_agreement(
    parameters: PC6Parameters, mesh: Mesh, *, samples: int = 24
) -> dict[str, Any]:
    """Compare the implicit field's sign against a solid-angle test on the mesh.

    The solid-angle (winding-number) test shares no code with the analytic
    membership test, so this is an independent check that the two representations
    in the contract describe the *same* solid.
    """
    p = parameters
    bounds = {
        "x": (-0.02 * p.length, 1.02 * p.length),
        "y": (-0.55 * p.beam, 0.55 * p.beam),
        "z": (-1.02 * p.draft, 0.02 * p.draft),
    }
    agree = 0
    inside = 0
    checked = 0
    for point in _lcg_points(samples, bounds):
        by_mesh = inside_by_winding_number(mesh, point)
        by_field = field_value(p, point[0], point[1], point[2]) < 0.0
        checked += 1
        inside += 1 if by_mesh else 0
        agree += 1 if by_mesh == by_field else 0
    return {
        "method": "sign of the implicit field vs winding number (solid angle) of the STL mesh",
        "samples": checked,
        "samples_inside_mesh": inside,
        "agreements": agree,
        "agreement_ratio": agree / checked if checked else float("nan"),
        "passed": agree == checked,
    }


def negative_controls(base: PC6Parameters | None = None) -> dict[str, Any]:
    """Run the two documented negative cases and record what happened.

    * a non-positive draft must be rejected by parameter validation;
    * a deliberately broken loft (stem rake larger than the hull is long) must be
      caught by the *geometric* self-intersection detector, not only by the
      parameter domain check. The raw loft is built with ``validate=False`` for
      exactly this purpose.
    """
    start = base if base is not None else PC6Parameters()
    draft_evidence: dict[str, Any] = {}
    for label, draft in (("zero", 0.0), ("negative", -4.0)):
        try:
            generate(replace(start, draft=draft))
        except ValueError as error:
            draft_evidence[label] = {
                "rejected": True,
                "error": str(error)[:240],
            }
        else:  # pragma: no cover - would be a real defect
            draft_evidence[label] = {"rejected": False, "error": None}

    domain_error: str | None
    broken = replace(start, stem_rake=1.25 * start.length)
    raw = build_mesh(broken, validate=False)
    pairs = find_self_intersections(raw)
    try:
        generate(broken)
    except ValueError as error:
        domain_rejected = True
        domain_error = str(error)[:240]
    else:  # pragma: no cover - would be a real defect
        domain_rejected = False
        domain_error = None
    return {
        "non_positive_draft": {
            "cases": draft_evidence,
            "passed": all(item["rejected"] for item in draft_evidence.values()),
        },
        "self_intersecting_hull": {
            "construction": "raw loft with stem_rake = 1.25 * L (stem curve crosses the transom)",
            "intersecting_face_pairs": len(pairs),
            "first_pairs": [list(pair) for pair in pairs[:5]],
            "detected": bool(pairs),
            "parameter_validation_also_rejects": domain_rejected,
            "parameter_error": domain_error,
            "passed": bool(pairs) and domain_rejected,
        },
    }


def _sweep_record(models: Sequence[PC6Model]) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    for model in models:
        p = model.parameters
        rows.append(
            {
                "beam_m": p.beam,
                "displacement_volume_m3": model.hydrostatics.displacement_volume,
                "analytic_volume_m3": model.form.volume,
                "volume_relative_error": model.report["reference_volume_error"],
                "block_coefficient": model.hydrostatics.block_coefficient,
                "analytic_block_coefficient": model.form.block_coefficient,
                "waterplane_area_m2": model.hydrostatics.waterplane_area,
                "analytic_waterplane_area_m2": model.form.waterplane_area,
                "wetted_surface_area_m2": model.hydrostatics.wetted_surface_area,
                "faces": len(model.mesh.faces),
                "watertight": bool(model.report["topology"]["watertight"]),
                "euler_characteristic": model.report["topology"]["euler_characteristic"],
                "validation_passed": model.validation.passed,
            }
        )
    return {
        "beams_m": [row["beam_m"] for row in rows],
        "rows": rows,
        "all_watertight": all(row["watertight"] for row in rows),
        "all_validation_passed": all(row["validation_passed"] for row in rows),
        "max_abs_volume_relative_error": max(
            (abs(row["volume_relative_error"]) for row in rows), default=float("nan")
        ),
        "monotone_volume": all(
            rows[i]["displacement_volume_m3"] < rows[i + 1]["displacement_volume_m3"]
            for i in range(len(rows) - 1)
        ),
    }


def build_contract(
    output_dir: str | Path,
    parameters: PC6Parameters | None = None,
    *,
    sweep_beams: Sequence[float] = (24.0, 26.0, 28.0),
    channel_width: float = DEFAULT_WIDTH_M,
    channel_depth: float = DEFAULT_DEPTH_M,
    field_samples: int = 24,
    volume_tolerance: float = 0.005,
    write_sweep_stl: bool = False,
    test_summary: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Generate everything, write the STL/JSON artifacts, return the contract."""
    target = Path(output_dir)
    target.mkdir(parents=True, exist_ok=True)
    p = parameters if parameters is not None else PC6Parameters()

    model = generate(p, volume_tolerance=volume_tolerance)
    hull_path = write_stl(model.mesh, target / HULL_STL, header="TensorCAD PC6 hull L200 B24 T11")
    sweep_models = beam_sweep(sweep_beams, p, volume_tolerance=volume_tolerance)
    # The sweep is evidence, not payload: its numbers and topology flags are in the
    # contract and the meshes regenerate with one call (``beam_sweep``), so the STLs
    # are opt-in to keep the artifact directory small.
    sweep_paths: list[str] = []
    if write_sweep_stl:
        sweep_dir = target / SWEEP_DIR
        sweep_dir.mkdir(parents=True, exist_ok=True)
        for swept in sweep_models:
            path = write_stl(
                swept.mesh,
                sweep_dir / f"pc6_hull_B{swept.parameters.beam:.0f}.stl",
                header=f"TensorCAD PC6 hull L{swept.parameters.length:.0f} "
                       f"B{swept.parameters.beam:.0f} T{swept.parameters.draft:.0f}",
            )
            sweep_paths.append(str(path.relative_to(target)))

    channel = ice_channel(width=channel_width, depth=channel_depth)
    channel_mesh = build_channel_mesh(channel)
    channel_path = write_stl(
        channel_mesh, target / CHANNEL_STL, header="TensorCAD PC6 de Long-Pevek ice channel"
    )
    channel_topology = channel_measurements(channel_mesh, channel)

    agreement = field_agreement(p, model.mesh, samples=field_samples)
    controls = negative_controls(p)

    hydro = model.hydrostatics
    form = model.form
    hydro_errors = {
        "displacement_volume": (hydro.displacement_volume - form.volume) / form.volume,
        "waterplane_area": (hydro.waterplane_area - form.waterplane_area) / form.waterplane_area,
        "lcb": (hydro.longitudinal_centre_buoyancy - form.lcb) / form.lcb,
        "lcf": (hydro.longitudinal_centre_flotation - form.lcf) / form.lcf,
        "vcb": (hydro.vertical_centre_buoyancy - form.vcb) / form.vcb,
        "block_coefficient": (hydro.block_coefficient - form.block_coefficient) / form.block_coefficient,
        "midship_section_area": (hydro.maximum_section_area - form.midship_section_area)
        / form.midship_section_area,
    }
    topology_block = dict(model.report["topology"])
    sweep = _sweep_record(sweep_models)

    contract: dict[str, Any] = {
        "schema": CONTRACT_SCHEMA,
        "kind": "pc6_ice_geometry_contract",
        "title": "PC6 parametric hull form and de Long -> Pevek ice channel",
        "units": UNITS,
        "coordinate_frame": COORDINATE_FRAME,
        "waterline_z_m": 0.0,
        "parameters": p.to_dict(),
        "geometry": {
            "hull": {
                "stl": HULL_STL,
                "stl_sha256": _sha256(hull_path),
                "stl_bytes": hull_path.stat().st_size,
                "format": "binary STL, outward-oriented triangles, metres",
                "topology": topology_block,
                "principal_dimensions": model.report["principal_dimensions"],
                "dimension_relative_errors": model.report["dimension_relative_errors"],
                "bounding_box": model.report["bounding_box"],
                "sampling": model.report["sampling"],
                "construction": "closed rings per station: two shell sheets + bottom + "
                                "waterplane deck + flat transom; starboard and port sheets "
                                "share the vertices on the raked stem curve",
            },
            "channel": {
                "stl": CHANNEL_STL,
                "stl_sha256": _sha256(channel_path),
                "stl_bytes": channel_path.stat().st_size,
                "format": "binary STL, outward-oriented triangles, metres",
                "topology": channel_topology,
            },
            "hydrostatics": {
                "mesh": hydro.to_dict(),
                "analytic": form.to_dict(),
                "relative_errors": hydro_errors,
                "mass_displacement_t": hydro.displacement_mass / 1000.0,
            },
        },
        "implicit_field": {
            "hull": sampled_field(p),
            "agreement_with_mesh": agreement,
        },
        "ice_channel": channel.to_dict(),
        "parameter_sweep": sweep,
        "negative_controls": controls,
        "consumers": {
            "loader": "tensorcad.pc6_contract.load_contract(path) -> PC6Geometry",
            "DEM hull contact": {
                "surface": HULL_STL,
                "note": "rigid contact surface; the implicit field gives a cheap "
                        "inside/outside test for particles near the shell",
            },
            "FEM ice belt": {
                "surface": CHANNEL_STL,
                "note": "closed prism of the corridor (width x depth below the ice "
                        "surface); import as the ice-belt domain, refine freely",
            },
            "LBM/FVM flow domain": {
                "field": "implicit_field.hull",
                "note": "classify voxels by the sign of the field over the recorded "
                        "bounding box; the STL is the fallback for body-fitted meshes",
            },
        },
        "acceptance": {
            "watertight_closed_manifold": bool(topology_block["watertight"])
            and bool(topology_block["edge_manifold"])
            and topology_block["euler_characteristic"] == 2
            and topology_block["components"] == 1,
            "beam_sweep_24_28_regenerates": bool(sweep["all_watertight"])
            and bool(sweep["all_validation_passed"]),
            "implicit_field_agrees_with_mesh": bool(agreement["passed"]),
            "negative_cases_detected": bool(controls["non_positive_draft"]["passed"])
            and bool(controls["self_intersecting_hull"]["passed"]),
            "volume_tolerance": volume_tolerance,
            "measured_volume_relative_error": model.report["reference_volume_error"],
        },
        "artifacts": {
            "contract": CONTRACT_NAME,
            "hull_stl": HULL_STL,
            "channel_stl": CHANNEL_STL,
            "sweep_stl": sweep_paths,
        },
    }
    if test_summary is not None:
        contract["acceptance"]["pytest"] = test_summary

    payload = json.dumps(_clean(contract), sort_keys=True, separators=(",", ":")).encode("utf-8")
    contract["digest"] = {
        "algorithm": "sha256 over this document without the digest block, canonical JSON",
        "value": hashlib.sha256(payload).hexdigest(),
    }
    write_metadata(_clean(contract), target / CONTRACT_NAME)
    return contract


@dataclass(frozen=True)
class PC6Geometry:
    """A loaded and re-verified contract."""

    contract: dict[str, Any]
    parameters: PC6Parameters
    hull: Mesh
    channel: Mesh
    hydrostatics: dict[str, Any]

    def field(self, x: float, y: float, z: float) -> float:
        return field_value(self.parameters, x, y, z)

    def inside(self, x: float, y: float, z: float) -> bool:
        return field_value(self.parameters, x, y, z) < 0.0

    @property
    def channel_corridor(self) -> IceChannel:
        return ice_channel(
            width=float(self.contract["ice_channel"]["corridor_width_m"]),
            depth=float(self.contract["ice_channel"]["corridor_depth_m"]),
        )


def _parameters_from_contract(contract: dict[str, Any]) -> PC6Parameters:
    table = contract["parameters"]
    values: dict[str, Any] = {}
    for name in PC6Parameters.__dataclass_fields__:
        if name in table and isinstance(table[name], dict) and "value" in table[name]:
            values[name] = table[name]["value"]
    return PC6Parameters(**values)


def load_contract(path: str | Path) -> PC6Geometry:
    """Load the contract, rebuild both solids and re-check them.

    Fails loudly when the STL is missing/corrupt, when either solid is not a closed
    manifold, or when the displacement of the rebuilt hull disagrees with the
    recorded value by more than the recorded tolerance.
    """
    source = Path(path)
    contract = json.loads(source.read_text(encoding="utf-8"))
    if contract.get("schema") != CONTRACT_SCHEMA:
        raise ValueError(
            f"{source}: schema {contract.get('schema')!r} is not {CONTRACT_SCHEMA!r}"
        )
    parameters = _parameters_from_contract(contract)
    parameters.require_valid()
    directory = source.parent
    hull = read_stl(
        directory / contract["artifacts"]["hull_stl"],
        units=UNITS,
        coordinate_frame=COORDINATE_FRAME,
    )
    channel = read_stl(
        directory / contract["artifacts"]["channel_stl"],
        units=UNITS,
        coordinate_frame=COORDINATE_FRAME,
    )
    for label, mesh in (("hull", hull), ("channel", channel)):
        topo = topology(mesh)
        if not topo["watertight"] or not topo["edge_manifold"]:
            raise ValueError(f"{source}: loaded {label} is not a closed manifold: {topo}")
        if topo["euler_characteristic"] != 2:
            raise ValueError(
                f"{source}: loaded {label} has Euler characteristic "
                f"{topo['euler_characteristic']}, expected 2"
            )
    recorded = contract["geometry"]["hydrostatics"]["mesh"]["displacement_volume_m3"]
    tolerance = float(contract["acceptance"]["volume_tolerance"])
    rebuilt = compute_hydrostatics(
        hull, 0.0, water_density=parameters.water_density, closed=True
    )
    error = (rebuilt.displacement_volume - recorded) / recorded
    if abs(error) > tolerance:
        raise ValueError(
            f"{source}: rebuilt hull displacement {rebuilt.displacement_volume:.3f} m^3 deviates "
            f"{error:+.2e} from the recorded {recorded:.3f} m^3 (tolerance {tolerance:g})"
        )
    return PC6Geometry(
        contract=contract,
        parameters=parameters,
        hull=hull,
        channel=channel,
        hydrostatics=rebuilt.to_dict(),
    )


def parse_pytest_summary(text: str) -> dict[str, Any]:
    """Extract ``failed``/``passed`` counts from a pytest ``-q`` output tail."""
    counts: dict[str, Any] = {}
    for label in ("passed", "failed", "error", "errors", "skipped"):
        match = re.search(rf"(\d+) {label}\b", text)
        if match:
            counts[label.rstrip("s")] = int(match.group(1))
    return {
        "passed": int(counts.get("passed", 0)),
        "failed": int(counts.get("failed", 0)) + int(counts.get("error", 0)),
        "source": "pytest -q tests/test_pc6.py",
    }


__all__ = [
    "CHANNEL_STL",
    "CONTRACT_NAME",
    "CONTRACT_SCHEMA",
    "HULL_STL",
    "PC6Geometry",
    "build_contract",
    "field_agreement",
    "load_contract",
    "negative_controls",
    "parse_pytest_summary",
]
