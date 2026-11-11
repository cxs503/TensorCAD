"""Apply the cap-adaptation change to tensorcad/parametric_hull.py (idempotent)."""
from __future__ import annotations

import re
import sys
from pathlib import Path

path = Path("tensorcad/parametric_hull.py")
src = path.read_text()

if "_adaptive_profile_rows" in src:
    print("already applied")
    sys.exit(0)

# 1. new helper, inserted just before tessellate_surface
helper = '''
def _adaptive_profile_rows(
    surface: "NURBSSurface",
    u_samples: Sequence[float],
    n_v: int,
    aspect_target: float,
) -> list[list[float]]:
    """Per-station profile samples that keep the *cap* cells isotropic.

    A tensor-product grid carries the same ``n_v`` profile points on every station.
    Near a tip, where the section radius ``r`` goes to zero, the circumferential cell
    size ``2*pi*r/n_theta`` goes to zero with it while the axial spacing ``dx`` does
    not: the last rings before the apex become fans of extremely thin triangles (the
    documented SUBOFF limitation: min scaled Jacobian 0.13, aspect ratio 7.9).  No
    amount of axial refinement fixes that, because refining ``dx`` only makes the
    meridional edge shorter while the circumferential edge stays short too.

    The fix is to drop columns: station ``i`` samples ``s_i + 1`` distinct profile
    parameters instead of ``n_v``, with ``s_i = round(L_i / (aspect_target * h_i))``
    clipped to ``[2, n_v - 1]``, where ``L_i`` is the polyline length of the profile
    and ``h_i`` the local axial spacing.  ``s_i`` therefore satisfies
    ``L_i / s_i ~ aspect_target * h_i``: the circumferential cell size matches the
    axial one.  The samples are placed on the *same* ``n_v`` index grid (index ``j``
    is snapped to ``round(j * s_i / (n_v - 1)) / s_i``), so coincident vertices appear
    in adjacent rows and the ordinary quad loop produces the zipper triangles; the
    downstream ``weld`` + ``remove_degenerate_faces`` (already in the pipeline for the
    degenerate stem column) collapses them.  The section polygon is coarsened only
    where the chord sagitta ``r * (1 - cos(L_i/(2 r s_i)))`` is already smaller than
    the axial spacing, and the volume effect is measured in
    ``docs/parametric_hull/CAP_ADAPTATION.md``.
    """
    if aspect_target <= 0.0:
        raise ParametricHullError("cap_aspect_target must be positive")
    xs = [surface.evaluate(u, 0.5)[0] for u in u_samples]
    rows: list[list[float]] = []
    for i, u in enumerate(u_samples):
        profile = [surface.evaluate(u, j / (n_v - 1)) for j in range(n_v)]
        length = sum(math.dist(profile[j], profile[j + 1]) for j in range(n_v - 1))
        if length <= 0.0:
            rows.append([j / (n_v - 1) for j in range(n_v)])
            continue
        if i == 0:
            h = abs(xs[1] - xs[0])
        elif i == len(u_samples) - 1:
            h = abs(xs[-1] - xs[-2])
        else:
            h = 0.5 * abs(xs[i + 1] - xs[i - 1])
        if h <= 0.0:
            rows.append([j / (n_v - 1) for j in range(n_v)])
            continue
        segments = int(round(length / (aspect_target * h)))
        segments = max(2, min(n_v - 1, segments))
        rows.append(
            [round(j * segments / (n_v - 1)) / segments for j in range(n_v)]
        )
    return rows


'''

anchor = "def tessellate_surface(\n"
assert anchor in src
src = src.replace(anchor, helper + anchor, 1)

# 2. signature: add the two keyword arguments
sig_start = src.index(anchor)
sig_end = src.index(") -> SurfaceMesh:", sig_start)
sig = src[sig_start:sig_end]
if "cap_adaptation" not in sig:
    src = src[:sig_end] + (
        "    cap_adaptation: bool = False,\n    cap_aspect_target: float = 1.0,\n"
    ) + src[sig_end:]

# 3. body: per-row profile samples
old_body = """    v_samples = [j / (n_v - 1) for j in range(n_v)]

    builder = MeshBuilder()
    grid: list[list[Vec3]] = [
        [surface.evaluate(u, v) for v in v_samples] for u in u_samples
    ]
"""
new_body = """    v_samples = [j / (n_v - 1) for j in range(n_v)]
    if cap_adaptation:
        rows_v = _adaptive_profile_rows(surface, u_samples, n_v, cap_aspect_target)
    else:
        rows_v = [v_samples] * len(u_samples)

    builder = MeshBuilder()
    grid: list[list[Vec3]] = [
        [surface.evaluate(u_samples[i], rows_v[i][j]) for j in range(n_v)]
        for i in range(len(u_samples))
    ]
"""
assert old_body in src
src = src.replace(old_body, new_body, 1)

path.write_text(src)
print(sig)
print("applied; cap_adaptation occurrences:", src.count("cap_adaptation"))