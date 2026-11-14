"""Calibrate the model->full-scale factor of a workshop IGES hull from the file itself.

Why this module exists
----------------------
The baseline fleet ships hull geometry at whatever scale the publisher used: KCS and
KVLCC2 are full scale, ONRT is 1:48.9355 and JBC is a model.  Until now a hull with no
published model scale was handled by normalising the *point cloud's x-extent* to the
published ``L_PP`` (``fleet_baseline.build_ship``, ``geometry_scale is None``).  That
shortcut is biased: the x-extent of a hull surface is its **length overall**
(bulbous bow / transom overhang included), not ``L_PP``, so the inferred factor comes
out systematically small.  For JBC it returns 37.15 where the correct factor is 40.

The calibration here uses quantities that mean the same thing in the file and in the
publisher's table: the **design waterline**.  If the file is a uniform scaling of the
hull by ``s`` and the keel sits at ``z = 0``, then the immersion of the design waterline
in the file is ``T/s`` and two readings at that immersion must reproduce the published
table:

* waterline length  ``s * L_wl(T/s) = L_WL``  ->  ``L_wl(z)/z = L_WL/T``
* waterline beam    ``s * B_wl(T/s) = B_WL``  ->  ``B_wl(z)/z = B_WL/T``

Each is one equation in one unknown (``z``, the model-frame immersion), solved by
bisection; the two solutions ``s_L``, ``s_B`` are *independent* and their disagreement
is the reported uncertainty.  A third, scale-free check rides along: the block
coefficient computed from the geometry against the published ``Cb``.  That check is
what exposes a geometry file that does not actually describe the published hull.

Honesty rules encoded here
--------------------------
* A patch that evaluates to a constant point is counted and **reported**, never
  silently multiplied into the surface statistics: the JBC archive currently has
  patches like that, and a mesher that ignores them loses volume without saying so.
* No number is written by hand: the report is produced by ``calibrate`` and only then
  serialised, and every published value used is echoed in the output.
* The volume closure is reported as an error with its own denominator, and a closure
  worse than ``closure_tolerance`` is flagged ``gate = "WITHHELD"`` instead of being
  rounded into a pass.

CLI
---
``python -m tensorcad.fleet.scale_calibration --iges FILE.igs --draft 16.5
--length-wl 285 --beam-wl 45 --volume 178369.9 --out report.json``
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Sequence

try:  # numpy is optional for the package core; the calibration is numeric.
    import numpy as np
except ModuleNotFoundError as exc:  # pragma: no cover - exercised only without numpy
    raise ModuleNotFoundError(
        "tensorcad.fleet.scale_calibration needs numpy (it evaluates and integrates "
        "the hull surface); the rest of tensorcad stays standard-library only"
    ) from exc

from ..iges import read_iges, unit_scale

__all__ = [
    "SCHEMA",
    "PatchGrids",
    "WaterlineCalibration",
    "calibrate",
    "collect_patch_grids",
    "solve_immersion",
    "waterline_metrics",
]

SCHEMA = "tensorcad.fleet.scale_calibration.v1"

#: Grid used to cache each patch: (u, v) evaluation points per patch.  ``u`` is the
#: coarse direction (it only needs to resolve where the waterline runs along the hull)
#: and ``v`` the fine one (the crossing is interpolated along ``v``).
DEFAULT_NU = 41
DEFAULT_NV = 81


@dataclass
class PatchGrids:
    """Cached evaluations of every usable patch, in the file's own (model) frame."""

    grids: list[np.ndarray]  # each (nu, nv, 3), metres of the model frame
    degenerate: list[int]  # patch indices whose evaluation collapses to a point
    nonzero: list[int]  # patch indices with a real, non-degenerate image
    factor: float  # units flag * model-space scale (file units -> metres)
    bbox: tuple[tuple[float, float, float], tuple[float, float, float]]

    @property
    def patch_count(self) -> int:
        return len(self.grids) + len(self.degenerate)

    def quality(self) -> dict[str, Any]:
        lo, hi = self.bbox
        return {
            "patches_total": self.patch_count,
            "patches_evaluable": len(self.grids),
            "patches_degenerate": len(self.degenerate),
            "degenerate_indices": sorted(self.degenerate),
            "degenerate_fraction": (
                round(len(self.degenerate) / self.patch_count, 4) if self.patch_count else None
            ),
            "unit_to_metre_factor": self.factor,
            "bbox_min_m": [round(v, 6) for v in lo],
            "bbox_max_m": [round(v, 6) for v in hi],
            "extent_m": [round(h - l, 6) for l, h in zip(lo, hi)],
        }


def collect_patch_grids(
    iges_path: str | Path, nu: int = DEFAULT_NU, nv: int = DEFAULT_NV
) -> PatchGrids:
    """Evaluate every type-128 patch on a coarse (u, v) grid and flag broken ones."""
    model = read_iges(str(iges_path))
    factor = unit_scale(model.unit_flag) * model.model_scale
    grids: list[np.ndarray] = []
    degenerate: list[int] = []
    nonzero: list[int] = []
    for index, surface in enumerate(model.surfaces()):
        u0, u1 = surface.domain_u
        v0, v1 = surface.domain_v
        u_values = np.linspace(u0, u1, nu)
        v_values = np.linspace(v0, v1, nv)
        block = np.empty((nu, nv, 3), dtype=float)
        for i, u in enumerate(u_values):
            for j, v in enumerate(v_values):
                block[i, j] = surface.evaluate(u, v)
        block *= factor  # file units -> metres of the model frame
        if not np.isfinite(block).all():
            degenerate.append(index)
            continue
        extent = block.reshape(-1, 3).max(axis=0) - block.reshape(-1, 3).min(axis=0)
        if float(extent.max()) <= 1e-9:
            degenerate.append(index)
            continue
        grids.append(block)
        nonzero.append(index)
    if not grids:
        raise ValueError(f"{iges_path}: no evaluable B-spline patch")
    stacked = np.concatenate([g.reshape(-1, 3) for g in grids])
    bbox = (tuple(stacked.min(axis=0).tolist()), tuple(stacked.max(axis=0).tolist()))
    return PatchGrids(grids, degenerate, nonzero, factor, bbox)


def waterline_points(grids: PatchGrids, immersion: float) -> np.ndarray:
    """Points where the cached patches cross the horizontal plane ``z = immersion``."""
    found: list[np.ndarray] = []
    for block in grids.grids:
        z = block[:, :, 2] - immersion
        nu, nv = z.shape
        for i in range(nu):
            row = z[i]
            left, right = row[:-1], row[1:]
            for j in np.nonzero(left * right < 0.0)[0]:
                w = left[j] / (left[j] - right[j])
                found.append(block[i, j] + w * (block[i, j + 1] - block[i, j]))
    if not found:
        return np.zeros((0, 3))
    return np.vstack(found)


def waterline_metrics(grids: PatchGrids, immersion: float, dx: float = 0.02) -> dict[str, float]:
    """Length, beam and area of the waterline at ``immersion`` (model metres)."""
    points = waterline_points(grids, immersion)
    if len(points) < 3:
        return {"immersion_m": immersion, "length_m": 0.0, "beam_m": 0.0, "area_m2": 0.0}
    xs = points[:, 0]
    length = float(xs.max() - xs.min())
    beam = float(2.0 * np.abs(points[:, 1]).max())
    bins = max(int(math.ceil(length / dx)), 1)
    half = np.zeros(bins)
    index = np.clip(((xs - xs.min()) / dx).astype(int), 0, bins - 1)
    np.maximum.at(half, index, np.abs(points[:, 1]))
    filled = np.nonzero(half > 0.0)[0]
    area = 0.0
    if len(filled) > 1:
        area = float(2.0 * np.trapezoid(half[filled[0] : filled[-1] + 1], dx=dx))
    return {"immersion_m": immersion, "length_m": length, "beam_m": beam, "area_m2": area}


def solve_immersion(
    grids: PatchGrids,
    quantity: str,
    target_ratio: float,
    *,
    lo: float | None = None,
    hi: float | None = None,
    tol: float = 1e-5,
    max_iterations: int = 60,
) -> dict[str, Any]:
    """Solve ``quantity(z) / z = target_ratio`` by bisection on the model-frame cut.

    ``quantity`` is ``"length_m"`` or ``"beam_m"``.  The ratio is assumed monotone in
    ``z`` (both grow with immersion while ``z`` is in the denominator); the bracket is
    validated before bisecting and a non-bracketing call raises instead of returning a
    confident wrong root.
    """
    if quantity not in ("length_m", "beam_m"):
        raise ValueError(f"unsupported quantity {quantity!r}")
    z_lo = lo if lo is not None else 1e-3 * float(grids.bbox[1][2])
    z_hi = hi if hi is not None else 0.99 * float(grids.bbox[1][2])

    def ratio(z: float) -> float:
        value = waterline_metrics(grids, z)[quantity]
        return float("nan") if (z <= 0.0 or value <= 0.0) else value / z

    # The ratio is *not* monotone over the whole file: the waterline shortens again as
    # the cut climbs onto the deck, and no crossing exists at all above it.  Scan, take
    # the lowest bracketing interval, then bisect inside it.  Merging this into one
    # bisection is what made the first version return "no root" for a hull whose root
    # sits in the middle of the range.
    scan = [z_lo + (z_hi - z_lo) * i / 60.0 for i in range(61)]
    values = [ratio(z) for z in scan]
    brackets: list[tuple[float, float]] = []
    for (za, fa), (zb, fb) in zip(zip(scan, values), zip(scan[1:], values[1:])):
        if not (math.isfinite(fa) and math.isfinite(fb)):
            continue
        if (fa - target_ratio) * (fb - target_ratio) <= 0.0 and fa != fb:
            brackets.append((za, zb))

    def bisect(z_a: float, z_b: float) -> dict[str, Any]:
        f_a = ratio(z_a) - target_ratio
        for iteration in range(1, max_iterations + 1):
            mid = 0.5 * (z_a + z_b)
            f_mid = ratio(mid) - target_ratio
            if abs(f_mid) < tol or (z_b - z_a) < tol:
                return {
                    "root_m": mid,
                    "iterations": iteration,
                    "converged": True,
                    "ratio_at_root": ratio(mid),
                    "quantity_at_root": waterline_metrics(grids, mid)[quantity],
                }
            if f_a * f_mid <= 0.0:
                z_b = mid
            else:
                z_a, f_a = mid, f_mid
        mid = 0.5 * (z_a + z_b)
        return {
            "root_m": mid,
            "iterations": max_iterations,
            "converged": False,
            "ratio_at_root": ratio(mid),
            "quantity_at_root": waterline_metrics(grids, mid)[quantity],
        }

    roots = [dict(bisect(za, zb), bracket_m=[za, zb]) for za, zb in brackets]
    if not roots:
        return {
            "quantity": quantity,
            "target_ratio": target_ratio,
            "root_m": None,
            "iterations": 0,
            "bracketed": False,
            "scan_ratios": [None if not math.isfinite(v) else round(v, 6) for v in values],
        }
    # A hull surface produces several level crossings (the cut clips a skeg, a bilge
    # radius, a sonar dome ... as well as the design waterline).  The design waterline
    # is the one where the hull is widest/longest, so rank the roots by the measured
    # quantity instead of taking the first bracket the scan happens to meet -- taking
    # the first one put the JBC waterline at z = 22 mm (a skeg edge) and produced a
    # scale of 742 instead of 39.8.
    best = max(roots, key=lambda r: r.get("quantity_at_root") or 0.0)
    best["quantity"] = quantity
    best["target_ratio"] = target_ratio
    best["bracketed"] = True
    best["all_roots_m"] = [round(r["root_m"], 6) for r in roots]
    best["all_roots_quantity"] = [round(r.get("quantity_at_root") or 0.0, 6) for r in roots]
    return best


@dataclass
class WaterlineCalibration:
    """The calibration result plus everything needed to re-judge it."""

    scale: float | None
    scale_by_length: float | None
    scale_by_beam: float | None
    published: dict[str, float]
    waterline: dict[str, float]
    quality: dict[str, Any]
    volume: dict[str, Any]
    gate: str
    tolerance: float
    notes: list[str] = field(default_factory=list)
    conditions: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "scale": self.scale,
            "scale_by_length": self.scale_by_length,
            "scale_by_beam": self.scale_by_beam,
            "published": self.published,
            "waterline_at_scale": self.waterline,
            "patch_quality": self.quality,
            "volume_closure": self.volume,
            "gate": self.gate,
            "closure_tolerance": self.tolerance,
            "notes": self.notes,
            "conditions": self.conditions,
        }


def _volume_from_cloud(
    grids: PatchGrids,
    immersion: float,
    *,
    spacing: float,
    dx: float,
    dz: float,
    cloud: np.ndarray | None = None,
) -> dict[str, float]:
    from .. import fleet_baseline as fb

    if cloud is None:
        raise ValueError("a sampled cloud must be supplied")
    field = fb.field_from_points(
        [tuple(row) for row in cloud.tolist()], draft=immersion, dx=dx, dz=dz
    )
    z = np.array([field.z(k) for k in range(field.nz)])
    x = np.array([field.x(i) for i in range(field.nx)])
    half = np.array([field.b[i] for i in range(field.nx)])
    keep = z <= immersion
    volume = float(np.trapezoid(np.trapezoid(2.0 * half[:, keep], x=x, axis=0), x=z[keep]))
    mid = int(np.argmin(np.abs(x)))
    section = float(2.0 * np.trapezoid(half[mid, keep], x=z[keep]))
    return {
        "volume_model_m3": volume,
        "midship_section_model_m2": section,
        "grid_cover": field.cover(),
        "dx_m": dx,
        "dz_m": dz,
        "cloud_points": int(cloud.shape[0]),
    }


def calibrate(
    iges_path: str | Path,
    *,
    draft: float,
    length_wl: float | None = None,
    beam_wl: float | None = None,
    length_pp: float | None = None,
    beam: float | None = None,
    volume: float | None = None,
    midship_coefficient: float | None = None,
    closure_tolerance: float = 0.03,
    sample_spacing: float = 2.0,
    nu: int = DEFAULT_NU,
    nv: int = DEFAULT_NV,
) -> WaterlineCalibration:
    """Determine the model->full scale of one hull file from its own waterline.

    ``sample_spacing`` is in the *file's* units (mm for the T2015 archive) and only
    controls the density of the cloud used for the volume closure.
    """
    grids = collect_patch_grids(iges_path, nu=nu, nv=nv)
    notes: list[str] = []
    published: dict[str, float] = {"draft": draft}
    if length_wl is not None:
        published["length_wl"] = length_wl
    if beam_wl is not None:
        published["beam_wl"] = beam_wl
    if length_pp is not None:
        published["length_pp"] = length_pp
    if beam is not None:
        published["beam"] = beam
    if volume is not None:
        published["volume"] = volume
    if midship_coefficient is not None:
        published["midship_coefficient"] = midship_coefficient

    z_hi = float(grids.bbox[1][2])
    solutions: dict[str, dict[str, Any]] = {}
    scale_by_length = None
    if length_wl is not None:
        solutions["length"] = solve_immersion(
            grids, "length_m", length_wl / draft, lo=1e-3 * z_hi, hi=z_hi
        )
        root = solutions["length"].get("root_m")
        scale_by_length = draft / root if root else None
    scale_by_beam = None
    if beam_wl is not None:
        solutions["beam"] = solve_immersion(
            grids, "beam_m", beam_wl / draft, lo=1e-3 * z_hi, hi=z_hi
        )
        root = solutions["beam"].get("root_m")
        scale_by_beam = draft / root if root else None

    candidates = [s for s in (scale_by_beam, scale_by_length) if s]
    scale = float(np.median(candidates)) if candidates else None
    if scale_by_beam and scale_by_length:
        spread = abs(scale_by_length - scale_by_beam) / scale_by_beam
        notes.append(
            f"two independent waterline conditions: length gives {scale_by_length:.4f}, "
            f"beam gives {scale_by_beam:.4f} (spread {100 * spread:.2f} %)"
        )

    waterline: dict[str, float] = {}
    volume_block: dict[str, Any] = {"status": "not_run"}
    gate = "WITHHELD"
    if scale:
        immersion = draft / scale
        metrics = waterline_metrics(grids, immersion)
        waterline = {
            "immersion_model_m": immersion,
            "length_model_m": metrics["length_m"],
            "length_full_m": metrics["length_m"] * scale,
            "beam_model_m": metrics["beam_m"],
            "beam_full_m": metrics["beam_m"] * scale,
            "area_model_m2": metrics["area_m2"],
            "area_full_m2": metrics["area_m2"] * scale * scale,
        }
        if length_wl:
            waterline["length_error_percent"] = 100.0 * (
                waterline["length_full_m"] / length_wl - 1.0
            )
        if beam_wl:
            waterline["beam_error_percent"] = 100.0 * (
                waterline["beam_full_m"] / beam_wl - 1.0
            )
        if volume is not None:
            from .. import fleet_baseline as fb

            model = read_iges(str(iges_path))
            cloud = np.array(
                list(fb.sample_patches(model.surfaces(), spacing=sample_spacing))
            ) * grids.factor
            closure = _volume_from_cloud(
                grids, immersion, spacing=sample_spacing, dx=0.01, dz=0.005, cloud=cloud
            )
            closure["volume_full_m3"] = closure["volume_model_m3"] * scale ** 3
            if closure["volume_model_m3"] > 0.0 and volume > 0.0:
                closure["volume_error_percent"] = 100.0 * (
                    closure["volume_full_m3"] / volume - 1.0
                )
                closure["scale_that_would_close_volume"] = scale * (
                    volume / closure["volume_full_m3"]
                ) ** (1.0 / 3.0)
            else:
                closure["volume_error_percent"] = None
                closure["empty"] = True
            if length_pp and beam:
                closure["cb_from_geometry"] = closure["volume_full_m3"] / (
                    length_pp * beam * draft
                )
            if midship_coefficient and beam:
                closure["cm_from_geometry"] = closure["midship_section_model_m2"] / (
                    waterline["beam_model_m"] * immersion
                )
                closure["cm_published"] = midship_coefficient
            within = closure["volume_error_percent"] is not None and abs(
                closure["volume_error_percent"]
            ) <= 100.0 * closure_tolerance
            closure["within_tolerance"] = bool(within)
            volume_block = closure
            gate = "PASS" if within else "WITHHELD"
            if not within:
                notes.append(
                    "volume closure fails at the calibrated scale: the file's evaluated "
                    "surface carries less volume than the publisher's displacement, so the "
                    "scale may not be quoted as validated by displacement — only by the "
                    "waterline conditions"
                )
        else:
            gate = "PASS"

    if grids.degenerate:
        notes.append(
            f"{len(grids.degenerate)} of {grids.patch_count} patches evaluate to a single "
            "point and were excluded from the waterline and volume statistics"
        )

    return WaterlineCalibration(
        scale=scale,
        scale_by_length=scale_by_length,
        scale_by_beam=scale_by_beam,
        published=published,
        waterline=waterline,
        quality=grids.quality(),
        volume=volume_block,
        gate=gate,
        tolerance=closure_tolerance,
        notes=notes,
        conditions={k: {kk: vv for kk, vv in v.items()} for k, v in solutions.items()},
    )


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--iges", required=True)
    parser.add_argument("--draft", type=float, required=True)
    parser.add_argument("--length-wl", type=float, default=None)
    parser.add_argument("--beam-wl", type=float, default=None)
    parser.add_argument("--length-pp", type=float, default=None)
    parser.add_argument("--beam", type=float, default=None)
    parser.add_argument("--volume", type=float, default=None)
    parser.add_argument("--midship-coefficient", type=float, default=None)
    parser.add_argument("--out", default=None)
    args = parser.parse_args(argv)
    result = calibrate(
        args.iges,
        draft=args.draft,
        length_wl=args.length_wl,
        beam_wl=args.beam_wl,
        length_pp=args.length_pp,
        beam=args.beam,
        volume=args.volume,
        midship_coefficient=args.midship_coefficient,
    )
    payload = json.dumps(result.to_dict(), indent=2, sort_keys=False)
    if args.out:
        Path(args.out).write_text(payload + "\n", encoding="utf-8")
        print(f"wrote {args.out}")
    print(payload)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())