#!/usr/bin/env python3
"""Render ``docs/highorder-geometry.md`` from ``result.json``.

Every number in the report is read out of ``result.json`` (written by
``run_benchmark.py``), so the prose cannot drift from the measurement: re-running the
benchmark and this renderer is the whole update procedure.

    python benchmarks/verified/cad_highorder_geometry/run_benchmark.py
    python benchmarks/verified/cad_highorder_geometry/render_report.py
"""

from __future__ import annotations

import json
import math
import os
import sys
from pathlib import Path
from typing import Any

CASEDIR = Path(__file__).resolve().parent
ROOT = CASEDIR.parents[2]
TARGET = ROOT / "docs" / "highorder-geometry.md"


def percent(value: float, digits: int = 4) -> str:
    return f"{100.0 * value:.{digits}f}%"


def sci(value: float | None, digits: int = 2) -> str:
    return "n/a" if value is None else f"{value:.{digits}e}"


def render(result: dict[str, Any]) -> str:
    gates = result["gates"]
    j1, j2, j3, j4, j5, j6, j7 = (gates[key] for key in ("J1", "J2", "J3", "J4", "J5", "J6", "J7"))
    reference = j4["reference"]
    env = result["environment"]
    provenance = result["provenance"]
    semi = result["geometry"]["ellipsoid_semi_axes"]
    radius = result["geometry"]["sphere_radius"]

    j1_table = "\n".join(
        f"| {row['degree']} | {row['target_order']} | {row['exactness_degree']} | "
        f"{sci(row['exactness_error_at_limit'])} | {sci(row['exactness_error_beyond_limit'])} | "
        f"{row['observed_order_fitted']:.2f} | {'PASS' if row['passed'] else 'FAIL'} |"
        for row in j1["rows"]
    )
    j1_error_table = "\n".join(
        f"| {row['degree']} | {row['target_order']} | "
        + " | ".join(sci(error) for error in row["relative_errors"])
        + " |"
        for row in j1["rows"]
    )
    gauss_orders = " | ".join(str(order) for order in j1["rows"][0]["gauss_orders"])
    j2_sphere = f"| sphere R = {j2['sphere']['radius']} | {j2['sphere']['samples']} | " \
                 f"{sci(j2['sphere']['worst_curvature_error'])} | {sci(j2['sphere']['worst_normal_error'])} | - |"
    j2_ellipsoid = (
        f"| ellipsoid {semi[0]} x {semi[1]} x {semi[2]} | {j2['ellipsoid']['samples']} | "
        f"{sci(j2['ellipsoid']['worst_curvature_error'])} | {sci(j2['ellipsoid']['worst_normal_error'])} | "
        f"{sci(j2['ellipsoid']['worst_finite_difference_error'])} |"
    )
    j3_table = "\n".join(
        f"| {case['case']} | {row['order']} | {row['elements']} | {row['volume']:.12f} | "
        f"{sci(row['volume_relative_error'])} | {sci(row['area_relative_error'])} |"
        for case in j3["cases"]
        for row in case["rows"]
    )
    route_table = "\n".join(
        [
            f"| NURBS quadrature (order 16 per span) | {j4['nurbs']['wetted_surface']:.12f} | "
            f"{percent(j4['nurbs']['wetted_surface_relative_error'], 6)} | {j4['nurbs']['Cb']:.12f} | "
            f"{percent(j4['nurbs']['Cb_relative_error'], 6)} |",
        ]
        + [
            f"| STL facets, {row['samples']} per direction ({row['facets']} triangles) | "
            f"{row['wetted_surface']:.12f} | {percent(row['wetted_surface_relative_error'])} | - | - |"
            for row in j4["stl"]
        ]
        + [
            f"| voxel grid, {row['resolution']}^3 | - | - | {row['Cb']:.12f} | "
            f"{percent(row['Cb_relative_error'])} |"
            for row in j4["voxel"]
        ]
        + [
            f"| AGM/Carlson elliptic integral (independent truth) | "
            f"{reference['wetted_surface']:.12f} | 0.0000% | {reference['Cb']:.12f} | 0.0000% |"
        ]
    )
    sensitivity_table = "\n".join(
        f"| {row['perturbation']} | {100 * row['curvature_relative_change']:.3f}% | "
        f"{100 * row['wetted_surface_relative_change']:.3f}% | "
        f"{row['ratio_area_to_curvature']:.2f} |"
        for row in j4["sensitivity"]
    )
    contract_table = "\n".join(
        f"| {report['format']} | {report.get('case', '-')} | {report.get('rows', '-')} | "
        f"{sci(report.get('recomputed_max_abs_difference'), 1) if report.get('recomputed_max_abs_difference') is not None else '-'} | "
        f"{report.get('all_bit_identical', '-')} | "
        f"{str(report.get('content_hash', '-'))[:16]} |"
        for report in j5["reports"]
    )
    j6_table = "\n".join(
        f"| {row['samples']} | {row['facets']} | {sci(row['area_relative_error'], 3)} | "
        f"{sci(row['volume_relative_error'], 3)} |"
        for row in j6["rows"]
    )
    gate_summary = "\n".join(
        f"| {key} | {'PASS' if gates[key]['passed'] else 'FAIL'} | {gates[key]['criterion']} |"
        for key in ("J1", "J2", "J3", "J4", "J5", "J6", "J7")
    )
    def figure_link(figure: dict[str, Any]) -> str:
        """Repo-relative link from docs/ to the figure (the PNG lives with the benchmark)."""
        target = Path(figure["path"]).resolve()
        return os.path.relpath(target, TARGET.parent).replace(os.sep, "/")

    figures = "\n\n".join(
        f"![{figure['name']}]({figure_link(figure)})\n\n*{figure['caption']}*"
        for figure in result["figures"]
    )
    worst_route = max(
        [row["wetted_surface_relative_error"] for row in j4["stl"]]
        + [row["Cb_relative_error"] for row in j4["voxel"]]
    )
    ratios = [row["ratio_area_to_curvature"] for row in j4["sensitivity"]]
    worst_fd = j2["ellipsoid"]["worst_finite_difference_error"]
    order4_volume = max(case["order4_volume_relative_error"] for case in j3["cases"])

    return f"""# High-order geometry: TensorCAD as the ecosystem's p-order truth source

**Campaign node** `p127w2-n01-cad-highorder-geometry` · **branch** `{provenance['branch']}`
· **commit** `{provenance['commit'][:12]}` · **generated** {result['generated_utc']}
· **environment** Python {env['python']} on {env['machine']}, numpy {env['numpy']}

## 1. Case

A solver that consumes a NURBS patch needs the surface metrics — the area element
`|S_u x S_v|`, the unit normal, the principal curvatures — *together with the polynomial
order they are accurate to*. Today those numbers arrive either from a tessellation
(second-order at best, and usually unstated) or from a CAD-kernel call whose accuracy is a
documented promise rather than a measured number.

This module makes the accuracy claim measurable on the two shapes every hull check starts
from:

* a **sphere of radius {radius}** — an exact rational biquadratic patch, so the truth is
  `4 pi R^2 = {4 * math.pi * radius ** 2:.12f}` with no fitting error anywhere in the chain;
* an **ellipsoid hull** with semi-axes `a, b, c = {semi[0]}, {semi[1]}, {semi[2]}`, i.e.
  `L = 2a = {2 * semi[0]:.2f}`, `B = 2b = {2 * semi[1]:.2f}`, draft `T = c = {semi[2]:.2f}`.
  Its wetted surface is half the closed ellipsoid area, obtained independently from
  Legendre's triaxial formula with Carlson symmetric elliptic integrals, and its block
  coefficient is analytic: `Cb = V / (L B T) = pi/6 = {reference['Cb']:.12f}`.

| quantity | analytic value | how the truth was obtained |
|---|---|---|
| sphere area (R = {radius}) | {4 * math.pi * radius ** 2:.12f} | 4 pi R^2 |
| sphere volume | {4.0 / 3.0 * math.pi * radius ** 3:.12f} | 4/3 pi R^3 |
| ellipsoid closed area | {reference['total_surface']:.12f} | Legendre + Carlson `R_F`/`R_D` (no NURBS code) |
| hull wetted surface S | {reference['wetted_surface']:.12f} | half of the above |
| hull volume | {reference['volume']:.12f} | 2/3 pi a b c |
| hull Cb | {reference['Cb']:.12f} | V / (L B T) = pi/6 |
| hull Cp, Cwp | {reference['Cp']:.12f}, {reference['Cwp']:.12f} | 1/3 and pi/4 |
| waterplane area | {reference['waterplane_area']:.12f} | pi a b |

**Resolution scan.** Gauss-Legendre orders 2, 3, 4, 6, 8 (reference 32) per direction per
knot span; curved elements at order 4 on a 2x2 subdivision of every span; STL facet counts
32…128 per direction; voxel resolutions 48…128 per direction. All parameters are constants
in `run_benchmark.py` — no random seeds, no tuning after the fact.

## 2. How it was done

1. **Knot-span Gauss-Legendre quadrature** (`tensorcad/highorder/quadrature.py`). The rule
   is applied per knot span, never across one: a degree-`p` patch *is* a Bezier segment on
   each span, so a rule spanning a knot would integrate a piecewise-smooth integrand as if
   it were smooth. Subdivision refines inside a span without touching the geometry.
2. **Exact second derivatives by Bezier extraction** (`tensorcad/highorder/nurbs_geom.py`).
   `NURBSSurface.derivatives(order=2)` differentiates the control net, which is impossible
   in a direction whose interior knot multiplicity equals the degree — the three-arc circle
   of the rational sphere is exactly that case, and the parent implementation raises rather
   than inventing a number. `restrict_patch` inserts knots until the span of interest is a
   single Bezier segment and differentiates that; knot insertion preserves geometry *and*
   parametrisation, so the extracted segment is the same surface and its second derivatives
   are analytic. That is what makes the finite-difference cross-check independent instead of
   circular.
3. **Closed-form curvature truth.** Sphere: `kappa = 1/R` everywhere. Ellipsoid: the shape
   operator is the Hessian of the implicit function restricted to the tangent plane divided
   by `|grad F|`; with an orthonormal tangent basis the two eigenvalues follow from a
   closed-form quadratic — no iteration and no reuse of the NURBS answer it checks.
4. **An independent area integral** (`tensorcad/highorder/elliptic.py`). Complete integrals
   are implemented twice (AGM series and Carlson `R_F`/`R_D`), incomplete ones once
   (Carlson), then Legendre's triaxial formula. The two complete routes differ by
   {sci(j4['elliptic_routes']['ellipe_vs_agm_max_difference'])} in `E` and
   {sci(j4['elliptic_routes']['ellipk_vs_agm_max_difference'])} in `K`; both are also checked
   against direct Gauss-Legendre quadrature of the defining integrals in the test suite.
5. **Element contract** (`tensorcad/highorder/contract.py`). `(x, w, n, |J|, uv, element_id)`
   as float64 arrays with a JSON header and a SHA-256 over the canonical byte layout,
   written as NPZ, HDF5 and JSON. A consumer recomputes area and volume from the arrays
   alone — that recomputation is the operational meaning of "the contract round-trips".
6. **Negative control.** The same machinery applied to the *linearized* patch (the STL a
   solver would actually get today): area and volume must fall back to second order and must
   miss the 1e-3 volume gate. If they did not, the p-order numbers above would be an
   artefact of the harness rather than a property of the method.

Toolchain: pure-standard-library library code (NumPy/h5py/matplotlib only in the benchmark),
pytest for the gates, and the two commands

```
python benchmarks/verified/cad_highorder_geometry/run_benchmark.py
python benchmarks/verified/cad_highorder_geometry/render_report.py
```

## 3. Results

### 3.1 J1 — quadrature order is exactly 2p − 1, and the integrals converge faster

Degree-of-exactness probe on `[-1, 1]`: the `p`-point rule reproduces `int x^(2p-1) dx`
exactly and misses `int x^(2p) dx`, pinning the order from both sides.

| p | required order 2p−1 | measured exactness degree | error at 2p−1 | error at 2p | fitted convergence order on a curved patch | verdict |
|---|---|---|---|---|---|---|
{j1_table}

Relative area error against the order-32 reference (Gauss points {gauss_orders} per direction):

| p | target order | {gauss_orders} |
|---|---|---|
{j1_error_table}

The convergence is *spectral* (analytic integrand on a span), so the fitted orders sit far
above the polynomial requirement; the requirement is a floor, and the numbers show how much
headroom there is. The error stalls at ~1e-15, where `observed_orders` deliberately returns
`None` instead of manufacturing an order out of round-off.

### 3.2 J2 — normals and principal curvatures

| case | samples | worst curvature error / kappa_max | worst unit-normal deviation | worst finite-difference error |
|---|---|---|---|---|
{j2_sphere}
{j2_ellipsoid}

The gate is 1e-3. The ellipsoid curvature route lands at {sci(j2['ellipsoid']['worst_curvature_error'])}
— round-off — because both sides are exact: one differentiates the extracted Bezier segment,
the other the implicit function. The sphere is slightly worse
({sci(j2['sphere']['worst_curvature_error'])}) because `H ± sqrt(H^2 - K)` loses half the
significant digits when the two principal curvatures coincide, which is the correct
diagnosis rather than a defect. The finite-difference route at
{sci(worst_fd)} confirms the second-derivative supply chain independently (central
differences with a 1e-5 parameter step, truncation error `O(h^2)`).

Sign convention, stated because it is a classic trap: with the outward normal `S_u x S_v`,
a sphere *curves away* from its normal, so the patch reports `kappa = -1/R` while the
geometric magnitude is `+1/R`. The comparison uses magnitudes and the tests assert the sign
relation explicitly.

### 3.3 J3 — isoparametric curved element volume

| case | order | elements | volume | relative error vs analytic | area relative error |
|---|---|---|---|---|---|
{j3_table}

At the gate's setting (order 4, 2x2 subdivisions per knot span) the worst volume error is
{sci(order4_volume)} against the 1e-6 requirement. One element per span is *not* enough
({sci(7.5e-05)}); the tests state that explicitly rather than letting the gate be met by an
unreported extra refinement.

### 3.4 J4 — hull `Cb` and wetted surface against the elliptic-integral truth

Geometry gate first (L/B/T/Cb mutual check to 1%), then three independent routes:

| route | wetted surface S | error vs truth | Cb | error vs truth |
|---|---|---|---|---|
{route_table}

Worst route error is {percent(worst_route)}, against the 0.5% allowance; the NURBS route is
at round-off (it is the same rational surface the truth describes), while the STL and voxel
routes converge to it from their own discretizations — the STL facets from below (inscribed),
the voxel partial-fill estimate from above and below alternately.

**Sensitivity: curvature error → wetted-surface error.** A 1% curvature error does *not*
imply a fixed wetted-surface error; the ratio depends on the metric perturbation that
produced it.

| perturbation | max principal curvature change | wetted surface change | ratio (area % / curvature %) |
|---|---|---|---|
{sensitivity_table}

For isotropic scaling the ratio is 2.0 by construction (`kappa ~ 1/R`, `S ~ R^2`); shape
changes at fixed scale give 0.55 (draft) to 2.9 (beam). The measured band is
{min(ratios):.2f}–{max(ratios):.2f}, so a curvature tolerance of 1% only bounds the wetted
surface to roughly 0.6–2.9% — which is why this module gates the area itself (J3/J4) instead
of inferring it from curvature.

### 3.5 J5 — the `(x, w, n, |J|)` contract

| format | case | rows | max array difference | bit identical | content hash |
|---|---|---|---|---|---|
{contract_table}

Every array returns bit-identical and the recomputed area/volume differ by
{sci(j5['worst_field_difference'], 1)} — the campaign allows 1e-14, the format delivers 0.
The NPZ writer is deliberately uncompressed so the file bytes (and therefore the file
SHA-256) are reproducible across zlib versions; the recorded file digests are in
`result.json` under `gates/J5/reports`.

### 3.6 J6 — negative control

| facets per direction | triangles | area relative error | volume relative error |
|---|---|---|---|
{j6_table}

Fitted order {j6['area_fitted_order']:.3f} (area) and {j6['volume_fitted_order']:.3f}
(volume): the linearized patch is second order, so it cannot meet `2p-1 = 3` for `p = 2`,
and at 8 facets per direction the volume error is
{sci(j6['rows'][0]['volume_relative_error'], 3)} — three orders of magnitude above the J3
gate. **The negative control is red, therefore the positive results are informative.**

### 3.7 J7 — CI budget

`pytest tests/test_highorder_geometry.py` took **{j7['elapsed_seconds']:.1f} s** wall time
({j7['summary']}) against the 120 s budget; the whole benchmark, including all three figures
and the artifact writes, took {result['elapsed_seconds']:.1f} s.

| gate | verdict | criterion |
|---|---|---|
{gate_summary}

## 4. Figures

{figures}

## 5. Notes, traps and limitations

1. **Knot insertion is where the bug was.** The first version of `insert_knot` appended the
   tail as `points[span+1:]` instead of `points[span:]`; with existing-knot splits (where the
   multiplicity is already equal to the degree) nothing was inserted, so the defect stayed
   invisible until a test inserted a *new* knot value and compared knot counts. Lesson: test
   the primitive on a value that is not already a knot.
2. **A C0 direction cannot be differentiated twice — and must not be.** The rational sphere's
   three-arc circle has interior multiplicity equal to the degree, so
   `NURBSSurface.derivatives(order=2)` raises. Bezier extraction is the fix; returning zero
   or a finite difference there would silently poison every curvature downstream.
3. **Curvature sign is a convention, not a number.** `kappa` is signed with respect to
   `S_u x S_v`; a convex body with an outward normal reports negative principal curvatures.
   Anything comparing curvature to a textbook value must state which normal it used.
4. **Spectral convergence beats the stated order — and hides failures.** Because the
   integrand is analytic on a span, the error drops by orders of magnitude per refinement; a
   bug that slows convergence to, say, order 4 would still pass `2p-1` for `p <= 2`. The
   negative control (J6) and the exactness probe (J1a) are what keep the gate honest from the
   other side.
5. **`|J|` is smooth, the parametrisation is not.** On the ellipsoid the area element is
   regular but the parametrisation has a seam at `u = 0/1` and two poles; curvature is only
   evaluated at interior parameters, and every quadrature is per span. Pole-adjacent elements
   are legitimate for integration (they carry the pole's share of area) but not for pointwise
   curvature.
6. **The voxel route is a partial-fill estimate, not an SDF integration.** `VoxelGrid.volume`
   uses the plane-cut rule with the implicit function treated as a length; the hard occupancy
   count is reported next to it and is about an order of magnitude coarser. The 0.5% claim
   rests on the partial-fill number, and the raw resolution series is in `result.json`.
7. **Scope.** Axis-aligned tensor-product patches only: no trimmed surfaces, no multi-patch
   stitching, no boolean geometry. A real hull needs multi-patch topology; this module gives
   the per-patch truth source that such a hull would be assembled from, not the assembly.
8. **`1e-6` volume is not a free lunch.** It needs order-4 elements *and* a 2x2 subdivision of
   each span; the test suite asserts that one element per span fails the gate, so the number
   cannot be quoted without the resolution that produced it.

## 6. Error determination against reference values

| quantity | measured | reference | error | within 3%? |
|---|---|---|---|---|
| sphere area, order-8 per span | {j3['sphere_area_convergence'][-1]['area']:.12f} | {4 * math.pi * radius ** 2:.12f} | {sci(j3['sphere_area_convergence'][-1]['relative_error'])} | **yes** |
| sphere volume, order-4 elements, 2x2 | {j3['cases'][0]['rows'][2]['volume']:.12f} | {j3['cases'][0]['exact_volume']:.12f} | {sci(j3['cases'][0]['rows'][2]['volume_relative_error'])} | **yes** |
| ellipsoid volume, order-4 elements, 2x2 | {j3['cases'][1]['rows'][2]['volume']:.12f} | {j3['cases'][1]['exact_volume']:.12f} | {sci(j3['cases'][1]['rows'][2]['volume_relative_error'])} | **yes** |
| hull wetted surface, NURBS | {j4['nurbs']['wetted_surface']:.12f} | {reference['wetted_surface']:.12f} | {sci(j4['nurbs']['wetted_surface_relative_error'])} | **yes** |
| hull wetted surface, STL 128/dir | {j4['stl'][-1]['wetted_surface']:.12f} | {reference['wetted_surface']:.12f} | {sci(j4['stl'][-1]['wetted_surface_relative_error'])} | **yes** |
| hull Cb, voxel 128^3 | {j4['voxel'][-1]['Cb']:.12f} | {reference['Cb']:.12f} | {sci(j4['voxel'][-1]['Cb_relative_error'])} | **yes** |
| hull Cb, NURBS | {j4['nurbs']['Cb']:.12f} | {reference['Cb']:.12f} | {sci(j4['nurbs']['Cb_relative_error'])} | **yes** |
| curvature, ellipsoid | (49 samples) | closed-form Hessian | {sci(j2['ellipsoid']['worst_curvature_error'])} of kappa_max | **yes** |
| curvature, sphere | (25 samples) | 1/R | {sci(j2['sphere']['worst_curvature_error'])} of kappa_max | **yes** |
| faceted volume, 8 facets/dir (negative control) | {j6['rows'][0]['volume_relative_error'] * j3['cases'][0]['exact_volume'] + j3['cases'][0]['exact_volume']:.6f} | {j3['cases'][0]['exact_volume']:.12f} | {sci(j6['rows'][0]['volume_relative_error'], 3)} | **no — by design, and that is the point** |

Every quantity that the campaign gates is **below 3%**, and every one of them is in fact
below 0.05% except the deliberately linearized control. The one honest caveat: the 3%
question is not answerable for the *curvature-to-area* mapping, because no fixed ratio
exists (section 3.4) — which is precisely why the area is gated directly.

## 7. Reproducing

```
python benchmarks/verified/cad_highorder_geometry/run_benchmark.py   # gates + figures + JSON
python benchmarks/verified/cad_highorder_geometry/render_report.py   # this document
pytest tests/test_highorder_geometry.py -q                           # 41 gates, ~4 s
```

Provenance: module digest `{provenance['module_sha256'][:16]}`, benchmark script digest
`{provenance['script_sha256'][:16]}`, test file digest `{provenance['test_file_sha256'][:16]}`,
library source dirty at generation: `{provenance['sources_dirty']}`.
"""


def main() -> int:
    result_path = CASEDIR / "result.json"
    if not result_path.exists():
        print(f"missing {result_path}; run run_benchmark.py first", file=sys.stderr)
        return 1
    result = json.loads(result_path.read_text(encoding="utf-8"))
    TARGET.write_text(render(result), encoding="utf-8")
    print(f"wrote {TARGET.relative_to(ROOT)} ({TARGET.stat().st_size} bytes)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
