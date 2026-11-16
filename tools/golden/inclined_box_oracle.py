"""Independent oracle for the hydrostatics of a rectangular barge at heel + trim.

Why this file exists
--------------------
:mod:`tensorcad.box_inclined` returns the exact hydrostatics of ``box ∩ halfspace`` by
clipping the six box faces and closing the cut with the waterplane polygon: a *boundary*
(triangulated polyhedron, divergence-theorem) route.  Its own test suite checks it against
one-axis closed forms (pure heel against the 2-D section oracle, pure trim against the
plane-cut moment), i.e. attitudes where an independent formula exists.  At combined heel
*and* trim no such formula exists, so the module's answer there has, until now, only been
compared with itself (``heel_and_trim_do_not_lose_volume``).

This oracle answers the same question from the *inside*: for every column ``(x, y)`` the
submerged set is the interval ``0 <= z <= z_w(x, y)``, so

    V   = ∫∫ z_w dA                       ∫ z dV = ∫∫ z_w²/2 dA
    x_B = ∫∫ x z_w dA / V                 y_B = ∫∫ y z_w dA / V
    A_wp = (1/n_z) · ∫∫ 1{0 < z_w < depth} dA        (the plane is tilted)

with ``z_w = clip((c - n_x x - n_y y)/n_z, 0, depth)``.  The depth integration is exact;
only the two horizontal directions are quadratured (composite Simpson), so the oracle has
a discretisation of its own and *no* shared step with the polyhedron route -- it never
builds a face, a cap polygon or a moment sum over triangles.

The offset ``c`` is the one the box floats at: it is solved here by bisection on
``V(c)`` against the requested displacement, exactly as the physical constraint reads.

Accuracy
--------
The integrands are piecewise smooth with kinks on the lines ``z_w = 0`` and
``z_w = depth``, so Simpson converges at second order.  :func:`convergence` measures that
order and the Richardson remainder on the very configuration the golden case uses; it is
the oracle's own error budget, produced on every run instead of claimed here.
"""
from __future__ import annotations

import math

import numpy as np

__all__ = [
    "oracle_hydrostatics",
    "waterline_normal",
    "convergence",
]


def waterline_normal(heel_deg: float, trim_deg: float) -> tuple[float, float, float]:
    """Up-normal of the waterline in body axes; the submerged side is ``n . p <= c``.

    Body axes as in ``docs/box_large_heel.md``: ``x`` forward, ``y`` starboard, ``z`` up
    from the keel, ``heel_deg > 0`` starboard side down, ``trim_deg > 0`` bow down.
    """
    phi = math.radians(heel_deg)
    tau = math.radians(trim_deg)
    n = (-math.sin(tau), -math.sin(phi) * math.cos(tau), math.cos(phi) * math.cos(tau))
    scale = math.sqrt(sum(component * component for component in n))
    return (n[0] / scale, n[1] / scale, n[2] / scale)


def _simpson_weights(count: int) -> np.ndarray:
    """Composite Simpson weights on ``count = n + 1`` equally spaced nodes (n even)."""
    if count < 3 or count % 2 == 0:
        raise ValueError("Simpson needs an odd number of nodes (n even)")
    weights = np.ones(count)
    weights[1:-1:2] = 4.0
    weights[2:-1:2] = 2.0
    return weights / 3.0


def _columns(
    n: tuple[float, float, float],
    c: float,
    length: float,
    beam: float,
    depth: float,
    resolution: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Submerged column height on the Simpson grid, with its axis weights and nodes."""
    dx = length / resolution
    dy = beam / resolution
    x = -0.5 * length + dx * np.arange(resolution + 1)
    y = -0.5 * beam + dy * np.arange(resolution + 1)
    wx = _simpson_weights(resolution + 1) * dx
    wy = _simpson_weights(resolution + 1) * dy
    nx, ny, nz = n
    plane = (c - nx * x[:, None] - ny * y[None, :]) / nz
    return np.clip(plane, 0.0, depth), x, y, wx, wy


def _volume_at(n, c, length, beam, depth, resolution) -> float:
    h, _, _, wx, wy = _columns(n, c, length, beam, depth, resolution)
    return float(wx @ h @ wy)


def _offset_for_volume(
    n: tuple[float, float, float],
    length: float,
    beam: float,
    depth: float,
    volume_m3: float,
    resolution: int,
    iterations: int = 200,
) -> float:
    """Bisect the waterline offset that displaces ``volume_m3`` (monotone in ``c``)."""
    nx, ny, nz = n
    corners = [
        nx * x + ny * y + nz * z
        for x in (-0.5 * length, 0.5 * length)
        for y in (-0.5 * beam, 0.5 * beam)
        for z in (0.0, depth)
    ]
    lo, hi = min(corners), max(corners)
    scale = max(length, beam, depth)
    # the bisection is run on the coarse grid first (cheap) and finished on the fine one
    coarse_grid = max(8, (resolution // 6) // 2 * 2)  # Simpson needs an even interval count
    for grid, steps in ((coarse_grid, iterations // 2), (resolution, iterations // 2)):
        for _ in range(steps):
            mid = 0.5 * (lo + hi)
            if _volume_at(n, mid, length, beam, depth, grid) < volume_m3:
                lo = mid
            else:
                hi = mid
            if hi - lo <= 1e-14 * scale:
                break
    return 0.5 * (lo + hi)


def oracle_hydrostatics(
    length: float,
    beam: float,
    depth: float,
    *,
    volume_m3: float,
    heel_deg: float = 0.0,
    trim_deg: float = 0.0,
    resolution: int = 1200,
    offset: float | None = None,
) -> dict[str, float]:
    """Volume, buoyancy centre and waterplane area of the barge, integrated from inside.

    ``resolution`` is the number of Simpson intervals per horizontal axis (even).  Pass
    ``offset`` only to pin the waterline for a convergence study; otherwise it is solved
    for the requested displacement.
    """
    if resolution % 2:
        raise ValueError("resolution must be even for composite Simpson")
    n = waterline_normal(heel_deg, trim_deg)
    if offset is None:
        offset = _offset_for_volume(n, length, beam, depth, volume_m3, resolution)
    h, x, y, wx, wy = _columns(n, offset, length, beam, depth, resolution)

    volume = float(wx @ h @ wy)
    if volume <= 0.0:
        raise ValueError("the waterline leaves no submerged volume")
    mx = float(wx @ (h * x[:, None]) @ wy)
    my = float(wx @ h @ (wy * y))
    mz = float(wx @ (0.5 * h * h) @ wy)
    inside = ((h > 0.0) & (h < depth)).astype(float)
    projected = float(wx @ inside @ wy)
    return {
        "displacement_volume_m3": volume,
        "lcb_x_m": mx / volume,
        "tcb_y_m": my / volume,
        "vcb_z_m": mz / volume,
        "waterplane_area_m2": projected / n[2],
        "waterline_offset_m": offset,
    }


def convergence(
    length: float,
    beam: float,
    depth: float,
    *,
    volume_m3: float,
    heel_deg: float,
    trim_deg: float,
    resolutions: tuple[int, ...] = (300, 600, 1200, 2400),
) -> dict[str, object]:
    """Self-convergence of the oracle on one configuration (its own error budget).

    The waterline offset is pinned to the finest run so that the only thing varying across
    the ladder is the horizontal quadrature: this measures the quadrature error, not the
    bisection.
    """
    runs = [
        oracle_hydrostatics(
            length,
            beam,
            depth,
            volume_m3=volume_m3,
            heel_deg=heel_deg,
            trim_deg=trim_deg,
            resolution=resolution,
            offset=(
                None
                if resolution != resolutions[-1]
                else oracle_hydrostatics(
                    length,
                    beam,
                    depth,
                    volume_m3=volume_m3,
                    heel_deg=heel_deg,
                    trim_deg=trim_deg,
                    resolution=resolutions[-1],
                )["waterline_offset_m"]
            ),
        )
        for resolution in resolutions
    ]
    finest = runs[-1]
    metrics = [key for key in finest if key != "waterline_offset_m"]
    table: dict[str, list[float]] = {}
    for key in metrics:
        reference = finest[key]
        table[key] = [
            abs(run[key] - reference) / abs(reference) if reference else abs(run[key] - reference)
            for run in runs
        ]
    # observed order between the last two refinements (both refine by 2)
    order: dict[str, float | None] = {}
    for key, errors in table.items():
        coarse, fine = errors[-3], errors[-2]
        order[key] = (
            math.log(coarse / fine) / math.log(2.0) if coarse > 0.0 and fine > 0.0 else None
        )
    return {
        "resolutions": list(resolutions),
        "config": {
            "length": length,
            "beam": beam,
            "depth": depth,
            "volume_m3": volume_m3,
            "heel_deg": heel_deg,
            "trim_deg": trim_deg,
        },
        "relative_deviation_from_finest": table,
        "observed_order": order,
        "finest": finest,
    }


def main(argv: list[str] | None = None) -> int:  # pragma: no cover - operator entry point
    import argparse
    import json

    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--length", type=float, default=10.0)
    parser.add_argument("--beam", type=float, default=4.0)
    parser.add_argument("--depth", type=float, default=3.0)
    parser.add_argument("--draft", type=float, default=1.5)
    parser.add_argument("--heel-deg", type=float, default=15.0)
    parser.add_argument("--trim-deg", type=float, default=3.0)
    parser.add_argument("--resolution", type=int, default=1200)
    parser.add_argument("--json", metavar="PATH")
    args = parser.parse_args(argv)

    volume = args.length * args.beam * args.draft
    hydro = oracle_hydrostatics(
        args.length,
        args.beam,
        args.depth,
        volume_m3=volume,
        heel_deg=args.heel_deg,
        trim_deg=args.trim_deg,
        resolution=args.resolution,
    )
    study = convergence(
        args.length,
        args.beam,
        args.depth,
        volume_m3=volume,
        heel_deg=args.heel_deg,
        trim_deg=args.trim_deg,
    )
    payload = {"hydrostatics": hydro, "convergence": study}
    if args.json:
        with open(args.json, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2)
            handle.write("\n")
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())