"""CLOUD field evidence for the intact-stability study (IS Code 2008 Part A).

Cases are the ones of report #5216, so every number can be cross-checked against a
published report instead of believed: a wall-sided barge (L=10, B=4, T=2, depth 4 m,
KG=1.20 m => GM0=0.466667 m) with its closed-form reference, the same barge at KG=1.55 m
(GM0=0.116667 m) as negative control, and a generated hull (L=12, B=2.4, T=0.9 + 0.6 m
freeboard, 81x25) with KG solved for GM0=0.300000 m, whose reference is a three-mesh
Richardson extrapolation (no closed form exists).

Outputs: ``cloud_fields.json`` (the field arrays), ``cloud_table.csv`` (table 1, long
form), ``cloud_results.json`` (enforced claims, IS Code table, uncertainty budget,
cross-check vs #5216), ``docs/stability_cloud.md``.

Run: python benchmarks/stability_cloud/run_cloud.py --output artifacts/m4801
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

from tensorcad import __version__
from tensorcad.bonjean import section_loops
from tensorcad.geometry import MeshBuilder
from tensorcad.gz import righting_arm_curve as gz_curve
from tensorcad.gz import with_freeboard
from tensorcad.hull import HullParameters, generate_hull
from tensorcad.righting_arm import (
    righting_arm_curve,
    wall_sided_righting_arm,
    wall_sided_validity_deg,
)

BARGE = dict(length=10.0, beam=4.0, draft=2.0, freeboard=2.0)
KG_BARGE, KG_HEAVY = 1.20, 1.55
HULL = dict(length=12.0, beam=2.4, draft=0.9, freeboard=0.6)
HULL_NX, HULL_NZ, GM0_TARGET = 81, 25, 0.300
LADDER = ((41, 13), (81, 25), (161, 41))
HEELS = tuple(float(a) for a in range(0, 91))
LADDER_HEELS = tuple(float(a) for a in range(0, 61, 5))
STEPS = (5.0, 2.0, 1.0, 0.5)
KG_GRID = tuple(round(0.60 + 0.025 * i, 6) for i in range(49))
#: IS Code 2008 Part A 2.2 criteria: (code, quantity, unit, required value).
CRITERIA = (("2.2.1", "A(0,30)", "m.rad", 0.055), ("2.2.2", "A(0,40)", "m.rad", 0.090),
            ("2.2.3", "A(30,40)", "m.rad", 0.030), ("2.2.4", "GZ(30)", "m", 0.20),
            ("2.2.5", "phi(GZmax)", "deg", 25.0), ("2.2.6", "GM0", "m", 0.15))
#: Published in report #5216 (mission 4753) for the same cases; recomputed here.
R5216 = {"barge": {"A(0,30)": 0.069437, "A(0,40)": 0.133014, "A(30,40)": 0.063576,
                   "GZ(30)": 0.288889, "GM0": 0.466667, "phi(GZmax)": 50.0},
         "barge_heavy": {"A(0,30)": 0.022548, "A(0,40)": 0.051131, "A(30,40)": 0.028584,
                         "GZ(30)": 0.113889, "GM0": 0.116667, "phi(GZmax)": 50.0},
         "hull": {"A(0,30)": 0.033962, "A(0,40)": 0.057996, "A(30,40)": 0.024035,
                  "GZ(30)": 0.124050, "GM0": 0.300000, "phi(GZmax)": 45.0}}


def box(length, beam, draft, freeboard):
    """Closed box: bottom z=0, waterline z=draft, deck z=draft+freeboard."""
    builder, corner = MeshBuilder(), {}
    top = draft + freeboard
    for x in (0.0, length):
        for y in (-0.5 * beam, 0.5 * beam):
            for z in (0.0, top):
                corner[(x, y, z)] = builder.add_vertex((x, y, z))
    lo, hi = -0.5 * beam, 0.5 * beam

    def quad(p0, p1, p2, p3):
        a, b, c, d = (corner[p] for p in (p0, p1, p2, p3))
        builder.add_triangle(a, b, c)
        builder.add_triangle(a, c, d)

    quad((0.0, lo, 0.0), (0.0, hi, 0.0), (length, hi, 0.0), (length, lo, 0.0))
    quad((0.0, lo, top), (length, lo, top), (length, hi, top), (0.0, hi, top))
    quad((0.0, lo, 0.0), (length, lo, 0.0), (length, lo, top), (0.0, lo, top))
    quad((0.0, hi, 0.0), (0.0, hi, top), (length, hi, top), (length, hi, 0.0))
    quad((0.0, lo, 0.0), (0.0, lo, top), (0.0, hi, top), (0.0, hi, 0.0))
    quad((length, lo, 0.0), (length, hi, 0.0), (length, hi, top), (length, lo, top))
    return builder.build()


def hull_model(nx, nz):
    return generate_hull(HullParameters(length=HULL["length"], beam=HULL["beam"],
                                        draft=HULL["draft"], nx=nx, nz=nz))


def hull_mesh(nx, nz):
    return with_freeboard(hull_model(nx, nz).mesh, HULL["freeboard"])


def hull_ref(nx, nz, kg=None):
    """KB, I_T, BM, KM of a mesh, and the KG that realises GM0 = 0.300 m."""
    h = hull_model(nx, nz).hydrostatics.to_dict()
    kb = h["vertical_centre_buoyancy_m"] - h["keel_z_m"]
    bm = h["waterplane_inertia_transverse_m4"] / h["displacement_volume_m3"]
    out = dict(nx=nx, nz=nz, volume_m3=h["displacement_volume_m3"],
               mass_kg=h["displacement_mass_kg"], cb=h["block_coefficient"],
               cp=h["prismatic_coefficient"], waterplane_area_m2=h["waterplane_area_m2"],
               kb_m=kb, it_m4=h["waterplane_inertia_transverse_m4"], bm_m=bm, km_m=kb + bm,
               kg_for_gm0_target_m=kb + bm - GM0_TARGET)
    out["gm0_m"] = None if kg is None else kb + bm - kg
    return out


def barge_ref(kg):
    kb, bm = 0.5 * BARGE["draft"], BARGE["beam"] ** 2 / (12.0 * BARGE["draft"])
    return dict(kb_m=kb, bm_m=bm, km_m=kb + bm, kg_m=kg, gm0_m=kb + bm - kg)


def ws_area(phi_deg, gm, bm):
    """int_0^phi GZ dphi of the wall-sided form: GM(1-cos)+(BM/2)(sec+cos-2)."""
    phi = math.radians(phi_deg)
    return gm * (1.0 - math.cos(phi)) + 0.5 * bm * (1.0 / math.cos(phi) + math.cos(phi) - 2.0)


def interp(heels, vals, x):
    for i in range(len(heels) - 1):
        if heels[i] <= x <= heels[i + 1]:
            span = heels[i + 1] - heels[i]
            w = 0.0 if span == 0.0 else (x - heels[i]) / span
            return vals[i] * (1.0 - w) + vals[i + 1] * w
    raise ValueError(f"heel {x} outside the sampled curve")


def trapz(heels, vals, a, b):
    """Trapezoid int GZ dphi in rad between two heels, linear at both ends."""
    nodes = [a, *[h for h in heels if a < h < b], b]
    return math.fsum((math.radians(n2) - math.radians(n1)) * 0.5
                     * (interp(heels, vals, n1) + interp(heels, vals, n2))
                     for n1, n2 in zip(nodes, nodes[1:]))


def peak(heels, vals):
    i = max(range(len(vals)), key=lambda k: vals[k])
    return heels[i], vals[i]


def crit_rows(gm0, heels, gz, theta_f=None):
    """Part A 2.2 checks on one measured curve; the comparison direction is explicit."""
    measured = {"2.2.1": trapz(heels, gz, 0.0, 30.0), "2.2.2": trapz(heels, gz, 0.0, 40.0),
                "2.2.3": trapz(heels, gz, 30.0, 40.0), "2.2.4": interp(heels, gz, 30.0),
                "2.2.5": peak(heels, gz)[0], "2.2.6": gm0}
    rows = [dict(criterion=c, quantity=q, unit=u, required=r, comparison=">=",
                 measured=measured[c], margin=measured[c] - r, passed=bool(measured[c] >= r))
            for c, q, u, r in CRITERIA]
    if theta_f is not None and theta_f < 40.0:
        alt = trapz(heels, gz, 0.0, theta_f)
        rows.append(dict(criterion="2.2.2(alt)", quantity=f"A(0,theta_f={theta_f:.3f})",
                         unit="m.rad", required=0.090, comparison=">=", measured=alt,
                         margin=alt - 0.090, passed=bool(alt >= 0.090),
                         note="theta_f = deck-edge immersion of the closed model; IS Code 2008 "
                              "2.2.2 lets it replace 40 deg when smaller"))
    return rows


def half_breadth(mesh, station, z):
    """Largest |y| of the section polygon at x = station on the line z = z."""
    best = None
    for loop in section_loops(mesh, station):
        for i in range(len(loop)):
            (y0, z0), (y1, z1) = loop[i], loop[(i + 1) % len(loop)]
            if (z0 - z) * (z1 - z) <= 0.0 and abs(z1 - z0) > 0.0:
                y = y0 + (z - z0) / (z1 - z0) * (y1 - y0)
                best = abs(y) if best is None or abs(y) > best else best
    return best


def spatial(mesh, xs, zs):
    """2y(x,z) over the wetted body plus the waterplane row y_hw(x) at z = 0."""
    rows = [[half_breadth(mesh, x, z) for z in zs] for x in xs]
    hw = [half_breadth(mesh, x, 0.0) or 0.0 for x in xs]
    contrib = [2.0 / 3.0 * h ** 3 for h in hw]
    it = math.fsum(0.5 * (xs[i + 1] - xs[i]) * (contrib[i] + contrib[i + 1])
                   for i in range(len(xs) - 1))
    return dict(x_m=list(xs), z_m=list(zs), half_breadth_rows_m=rows,
                half_breadth_at_waterline_m=hw, it_contribution_m4_per_m=contrib,
                it_from_stations_m4=it)


def density(xs, ys, hw):
    """rho(x,y) = y^2 masked to |y| <= y_hw(x); its double integral is I_T."""
    grid = [[float(y * y) if abs(y) <= h + 1e-12 else 0.0 for y in ys] for h in hw]
    it = 0.0
    for i in range(len(xs) - 1):
        for j in range(len(ys) - 1):
            it += (xs[i + 1] - xs[i]) * (ys[j + 1] - ys[j]) * 0.25 * (
                grid[i][j] + grid[i + 1][j] + grid[i][j + 1] + grid[i + 1][j + 1])
    return dict(x_m=list(xs), y_m=list(ys), density_y2_m2=grid, it_integral_m4=it)


def kn_of(heels, gz, kg):
    """KN(phi) = GZ(phi) + KG sin(phi): the loading-independent part of the arm."""
    return [v + kg * math.sin(math.radians(h)) for h, v in zip(heels, gz)]


def field_over_kg(heels, kn):
    """GZ(phi, KG) = KN(phi) - KG sin(phi): exact algebra on the measured curve."""
    return dict(kg_m=list(KG_GRID), heel_deg=list(heels),
                gz_m=[[k - kg * math.sin(math.radians(h)) for h, k in zip(heels, kn)]
                      for kg in KG_GRID])


def sweep_over_kg(heels, kn):
    """Criterion quantities as functions of the loading KG (integrated angular field)."""
    out = dict(kg_m=list(KG_GRID), A030=[], A040=[], A3040=[], GZ30=[], phi_max=[])
    for kg in KG_GRID:
        curve = [k - kg * math.sin(math.radians(h)) for h, k in zip(heels, kn)]
        out["A030"].append(trapz(heels, curve, 0.0, 30.0))
        out["A040"].append(trapz(heels, curve, 0.0, 40.0))
        out["A3040"].append(trapz(heels, curve, 30.0, 40.0))
        out["GZ30"].append(interp(heels, curve, 30.0))
        out["phi_max"].append(peak(heels, curve)[0])
    return out


def gci(fine, medium, coarse, ratio, safety=1.25):
    """Roache GCI on three levels, with an explicit fallback when it does not apply."""
    d21, d32 = fine - medium, medium - coarse
    res = dict(fine=fine, medium=medium, coarse=coarse, ratio=ratio,
               difference_fine_medium=d21, difference_medium_coarse=d32)
    if d21 == 0.0 or d32 == 0.0 or d21 * d32 <= 0.0:
        res.update(observed_order=None, extrapolated=None, gci95=None, u95=abs(d21),
                   fallback_reason="non-monotone or vanishing refinement differences: observed "
                                   "order undefined, no Richardson extrapolation; U95 falls back "
                                   "to |fine-medium|")
        return res
    order = abs(math.log(abs(d32 / d21)) / math.log(ratio))
    if not 0.1 <= order <= 5.0:
        res.update(observed_order=order, extrapolated=None, gci95=None, u95=abs(d21),
                   fallback_reason=f"observed order {order:.3f} outside [0.1, 5.0]: not in the "
                                   "asymptotic range, extrapolation withheld")
        return res
    res.update(observed_order=order,
               extrapolated=fine + (fine - medium) / (ratio ** order - 1.0),
               gci95=safety * abs(d21) / abs(fine) / (ratio ** order - 1.0),
               u95=abs(d21) / (ratio ** order - 1.0), fallback_reason=None)
    return res


def combine(**components):
    u = math.sqrt(math.fsum(v * v for v in components.values()))
    return dict(components=components, u=u, U95=2.0 * u)


def build():
    """Measure both cases and return (results, fields); every claim is enforced here."""
    claims, cases = [], {}
    fields = dict(schema="tensorcad.stability.cloud.fields/1", version=__version__)

    def claim(name, ok, detail):
        claims.append(dict(claim=name, passed=bool(ok), detail=detail))

    xs = [round(HULL["length"] * i / 80.0, 9) for i in range(81)]
    ys = [round(BARGE["beam"] * (i / 80.0 - 0.5), 9) for i in range(81)]
    zs = [round(1e-6 + (HULL["draft"] - 1e-6) * i / 60.0, 9) for i in range(61)]

    # ---- barge: measured curve, closed-form reference, angle ladder, criteria -----
    mesh = box(**BARGE)
    ref = barge_ref(KG_BARGE)
    curve = righting_arm_curve(mesh, kg=KG_BARGE, heels_deg=HEELS, waterline_z=BARGE["draft"])
    gz = [p.righting_arm for p in curve.points]
    gz2 = list(gz_curve(mesh, KG_BARGE, HEELS, waterline_z=BARGE["draft"]).gz)
    spread = max(abs(a - b) for a, b in zip(gz, gz2))
    volume = BARGE["length"] * BARGE["beam"] * BARGE["draft"]
    claim("barge geometry verified: L, B, T as stated, displacement L*B*T",
          abs(curve.displacement_volume - volume) / volume < 1e-12,
          f"V={curve.displacement_volume:.9f} m^3 vs L*B*T={volume:.9f} m^3")
    valid = wall_sided_validity_deg(breadth=BARGE["beam"], draft=BARGE["draft"],
                                    depth=BARGE["draft"] + BARGE["freeboard"])
    closed = [wall_sided_righting_arm(breadth=BARGE["beam"], draft=BARGE["draft"],
                                      kg=KG_BARGE, heel_deg=h) for h in HEELS]
    worst_in = max(abs(g - c) for h, g, c in zip(HEELS, gz, closed) if h <= valid)
    min_out = min(abs(g - c) for h, g, c in zip(HEELS, gz, closed) if h >= 60.0)
    claim("barge GZ field slice equals the wall-sided closed form inside its validity window",
          worst_in < 1e-9, f"worst |GZ-formula| = {worst_in:.2e} m over 0..{valid:.0f} deg")
    claim("the closed form stops being exact past its validity window (comparison can fail)",
          min_out > 1e-2, f"from 60 deg on, the smallest deviation is {min_out:.3e} m")
    claim("two independent waterline solvers agree on the whole barge curve", spread < 1e-9,
          f"max |righting_arm - gz| = {spread:.2e} m over 91 heels")
    anti = ws_area(30.0, ref["gm0_m"], ref["bm_m"])
    fine = [0.05 * i for i in range(601)]
    quad = trapz(fine, [wall_sided_righting_arm(breadth=BARGE["beam"], draft=BARGE["draft"],
                                                kg=KG_BARGE, heel_deg=h) for h in fine], 0.0, 30.0)
    claim("closed-form area antiderivative matches a fine trapezoid of the closed-form GZ",
          abs(anti - quad) / anti < 1e-8,
          f"A(0,30) = {anti:.9f} (antiderivative) vs {quad:.9f} (0.05 deg trapezoid)")
    ref30, ref40 = anti, ws_area(40.0, ref["gm0_m"], ref["bm_m"])
    ladder = []
    for step in STEPS:
        want = [round(step * i, 9) for i in range(int(HEELS[-1] / step) + 1)]
        pairs = [(h, v) for h, v in zip(HEELS, gz) if h in want]
        hh = [p[0] for p in pairs]
        a30, a40 = trapz(hh, [p[1] for p in pairs], 0.0, 30.0), trapz(hh, [p[1] for p in pairs], 0.0, 40.0)
        ladder.append(dict(step_deg=step, samples=len(hh), a030_measured=a30,
                           a030_reference=ref30, a030_relative_error=(a30 - ref30) / ref30,
                           a040_measured=a40, a040_reference=ref40,
                           a040_relative_error=(a40 - ref40) / ref40))
    tail = [r for r in ladder if r["step_deg"] <= 2.0]
    order = math.fsum(math.log(abs(a["a030_relative_error"]) / abs(b["a030_relative_error"]))
                      / math.log(a["step_deg"] / b["step_deg"])
                      for a, b in zip(tail, tail[1:])) / (len(tail) - 1)
    claim("the criterion integral converges at second order in the angle step",
          1.8 <= order <= 2.2, f"fitted order p = {order:.3f} over 2.0, 1.0, 0.5 deg")
    budget = combine(u_angle=abs(ladder[-2]["a030_measured"] - ladder[-1]["a030_measured"]) / 3.0,
                     u_solver=spread)
    crit = crit_rows(ref["gm0_m"], HEELS, gz, theta_f=valid)
    kn = kn_of(HEELS, gz, KG_BARGE)
    sp = spatial(mesh, xs, zs)
    fields["barge"] = dict(case="barge", kg_m=KG_BARGE, km_m=ref["km_m"], gm0_m=ref["gm0_m"],
                           heel_deg=list(HEELS), gz_m_at_kg=gz, gz_reference_m=closed, kn_m=kn,
                           gz_field=field_over_kg(HEELS, kn),
                           loading_sweep=sweep_over_kg(HEELS, kn), spatial=sp,
                           density=density(xs, ys, sp["half_breadth_at_waterline_m"]),
                           it_exact_m4=BARGE["length"] * BARGE["beam"] ** 3 / 12.0)
    cases["barge"] = dict(reference=dict(ref, validity_deg=valid), criteria=crit,
                          area_ladder=ladder, area_order=dict(order=order),
                          uncertainty=budget, solver_spread_m=spread,
                          curve_vs_closed_form=dict(worst_abs_m=worst_in,
                                                    min_deviation_beyond_validity_m=min_out))

    # ---- negative control: same hull, heavier loading -----------------------------
    hgz = [p.righting_arm for p in righting_arm_curve(
        mesh, kg=KG_HEAVY, heels_deg=HEELS, waterline_z=BARGE["draft"]).points]
    href = barge_ref(KG_HEAVY)
    hcrit = crit_rows(href["gm0_m"], HEELS, hgz, theta_f=valid)
    hkn = kn_of(HEELS, hgz, KG_HEAVY)
    fields["barge_heavy"] = dict(case="barge_heavy", kg_m=KG_HEAVY, km_m=href["km_m"],
                                 gm0_m=href["gm0_m"], heel_deg=list(HEELS), gz_m_at_kg=hgz,
                                 gz_field=field_over_kg(HEELS, hkn),
                                 loading_sweep=sweep_over_kg(HEELS, hkn))
    cases["barge_heavy"] = dict(reference=dict(href, validity_deg=valid), criteria=hcrit,
                             uncertainty=budget, solver_spread_m=spread)
    claim("the negative control flips the criteria (the comparison can fail)",
          len([r for r in hcrit if not r["passed"]]) >= 4
          and not [r for r in crit if not r["passed"]],
          "heavy KG fails " + str([r["criterion"] for r in hcrit if not r["passed"]])
          + "; design KG fails " + str([r["criterion"] for r in crit if not r["passed"]] or "none"))

    # ---- hull: KG solved for GM0 = 0.300 m, mesh ladder, GCI -----------------------
    kg = hull_ref(HULL_NX, HULL_NZ)["kg_for_gm0_target_m"]
    href = hull_ref(HULL_NX, HULL_NZ, kg)
    claim("hull geometry verified: Cb*L*B*T reproduces the hydrostatic displacement",
          abs(href["cb"] * HULL["length"] * HULL["beam"] * HULL["draft"] - href["volume_m3"]) < 1e-6,
          f"Cb={href['cb']:.6f}, V={href['volume_m3']:.6f} m^3, "
          f"L/B={HULL['length'] / HULL['beam']:.2f}, B/T={HULL['beam'] / HULL['draft']:.2f}")
    hmesh = hull_mesh(HULL_NX, HULL_NZ)
    hcurve = righting_arm_curve(hmesh, kg=kg, heels_deg=HEELS, waterline_z=0.0)
    hgz = [p.righting_arm for p in hcurve.points]
    probe = (0.0, 15.0, 30.0)
    alt = list(gz_curve(hmesh, kg, probe, waterline_z=0.0).gz)
    hspread = max(abs(interp(HEELS, hgz, h) - v) for h, v in zip(probe, alt))
    theta_f = math.degrees(math.atan(HULL["freeboard"] / (0.5 * HULL["beam"])))
    hcrit = crit_rows(href["gm0_m"], HEELS, hgz, theta_f=theta_f)
    meshes = []
    for nx, nz in LADDER:
        r = hull_ref(nx, nz, kg)
        c = righting_arm_curve(hull_mesh(nx, nz), kg=kg, heels_deg=LADDER_HEELS, waterline_z=0.0)
        arms = [p.righting_arm for p in c.points]
        meshes.append(dict(nx=nx, nz=nz, gm0_m=r["gm0_m"], kb_m=r["kb_m"], it_m4=r["it_m4"],
                           bm_m=r["bm_m"], a030=trapz(LADDER_HEELS, arms, 0.0, 30.0),
                           gz30=interp(LADDER_HEELS, arms, 30.0),
                           gz_by_heel={str(int(h)): v for h, v in zip(LADDER_HEELS, arms)}))
    ladder_gci = {name: gci(meshes[2][key], meshes[1][key], meshes[0][key], 2.0)
                  for name, key in (("gm0", "gm0_m"), ("a030", "a030"), ("gz30", "gz30"))}
    extrap = [dict(heel_deg=h,
                   value=gci(meshes[2]["gz_by_heel"][str(int(h))],
                             meshes[1]["gz_by_heel"][str(int(h))],
                             meshes[0]["gz_by_heel"][str(int(h))], 2.0)["extrapolated"])
              for h in LADDER_HEELS]
    hkn = kn_of(HEELS, hgz, kg)
    hsp = spatial(hull_model(HULL_NX, HULL_NZ).mesh, xs, zs)
    fields["hull"] = dict(case="hull", kg_m=kg, km_m=href["km_m"], gm0_m=href["gm0_m"],
                          heel_deg=list(HEELS), gz_m_at_kg=hgz, kn_m=hkn,
                          gz_field=field_over_kg(HEELS, hkn),
                          loading_sweep=sweep_over_kg(HEELS, hkn), spatial=hsp,
                          density=density(xs, ys, hsp["half_breadth_at_waterline_m"]),
                          it_exact_m4=href["it_m4"])
    hbudget = combine(u_solver=hspread, u_mesh=ladder_gci["a030"]["u95"])
    cases["hull"] = dict(reference=dict(href, kg_m=kg, theta_f_deg=theta_f), criteria=hcrit,
                         ladder=dict(meshes=meshes, gz_extrapolated=extrap, gci=ladder_gci),
                         uncertainty=hbudget, solver_spread_m=hspread)

    # ---- cross-check against the published numbers of report #5216 ----------------
    cross = {}
    for case in ("barge", "barge_heavy", "hull"):
        got = {r["quantity"]: r["measured"] for r in cases[case]["criteria"]
               if not r["criterion"].endswith("alt)")}
        cross[case] = {k: dict(this_run=got[k], report_5216=v,
                               relative_difference=None if not v else (got[k] - v) / v)
                       for k, v in R5216[case].items()}
    worst_cross = max(abs(row["relative_difference"]) for case in cross.values()
                      for row in case.values() if row["relative_difference"] is not None)
    claim("every criterion quantity reproduces report #5216 within 1 %", worst_cross < 0.01,
          f"worst |relative difference| = {100 * worst_cross:.3f} %")
    results = dict(schema="tensorcad.stability.cloud.results/1", version=__version__,
                   claims=claims, cases=cases, fields=fields,
                   cross_check_report_5216=cross,
                   uncertainty=dict(barge=budget, hull=hbudget),
                   provenance={
                       "barge": "box L=10, B=4, T=2, depth 4 m, waterline z=2, KG=1.20 m; "
                                "reference = wall_sided_righting_arm + elementary area integral",
                       "barge_heavy": "same box, KG=1.55 m (negative control)",
                       "hull": f"generate_hull L=12, B=2.4, T=0.9 + 0.6 m freeboard, "
                               f"{HULL_NX}x{HULL_NZ}, KG solved for GM0=0.300 m; reference = "
                               f"three-mesh Richardson extrapolation",
                       "report_cross_check": "report id 5216 (mission 4753) numbers, recomputed "
                                             "independently in this run"})
    return results, fields


def fmt(value, digits=6):
    if value is None:
        return "NO-DATA-TBD"
    if isinstance(value, bool):
        return "pass" if value else "FAIL"
    if isinstance(value, float):
        return f"{value:.{digits}f}" if abs(value) < 1e4 else f"{value:.3e}"
    return str(value)


def table_rows(results):
    """Table 1, long form: one row per level of one key quantity."""
    rows = []

    def add(block, case, kind, level, quantity, unit, value, reference=None, u95=None,
            required=None, passed=None, note=""):
        err = rel = None
        if reference is not None and value is not None:
            err = value - reference
            rel = err / reference if reference else None
        rows.append(dict(block=block, case=case, level_kind=kind, level=level,
                         quantity=quantity, unit=unit, value=value, reference=reference,
                         abs_error=err, rel_error=rel, u95=u95, required=required,
                         passed=passed, note=note))

    for case in ("barge", "barge_heavy", "hull"):
        data = results["cases"][case]
        field = results["fields"][case]
        if case in ("barge", "barge_heavy"):
            for h, g in zip(field["heel_deg"], field["gz_m_at_kg"]):
                inside = h <= data["reference"]["validity_deg"]
                ref = wall_sided_righting_arm(breadth=BARGE["beam"], draft=BARGE["draft"],
                                              kg=field["kg_m"], heel_deg=h) if inside else None
                add("heel", case, "heel_deg", h, "GZ", "m", g, ref, u95=1e-9,
                    note="closed form" if inside else
                    "beyond wall-sided validity window: reference not applicable")
        else:
            extrap = {row["heel_deg"]: row["value"] for row in data["ladder"]["gz_extrapolated"]}
            for h, g in zip(field["heel_deg"], field["gz_m_at_kg"]):
                add("heel", case, "heel_deg", h, "GZ", "m", g, extrap.get(h),
                    u95=2.0 * data["ladder"]["gci"]["gz30"]["u95"],
                    note="reference = Richardson extrapolation" if h in extrap
                    else "no extrapolation at this heel (outside the 0-60 deg ladder)")
        for row in data["criteria"]:
            add("criterion", case, "design_loading", f"KG={field['kg_m']:.4f}",
                row["quantity"], row["unit"], row["measured"], None,
                u95=data["uncertainty"]["U95"] if row["criterion"] in
                ("2.2.1", "2.2.2", "2.2.3") else None,
                required=row["required"], passed=row["passed"], note=row.get("note", ""))
        sweep = field["loading_sweep"]
        for i, k in enumerate(sweep["kg_m"]):
            for key, label, unit, req in (("A030", "A(0,30)", "m.rad", 0.055),
                                          ("A040", "A(0,40)", "m.rad", 0.090),
                                          ("A3040", "A(30,40)", "m.rad", 0.030),
                                          ("GZ30", "GZ(30)", "m", 0.20)):
                add("loading", case, "KG_m", k, label, unit, sweep[key][i],
                    required=req, passed=bool(sweep[key][i] >= req))
    for row in results["cases"]["hull"]["ladder"]["meshes"]:
        for label, key, unit in (("GM0", "gm0_m", "m"), ("A(0,30)", "a030", "m.rad"),
                                 ("GZ(30)", "gz30", "m")):
            add("mesh", "hull", "nx_x_nz", f"{row['nx']}x{row['nz']}", label, unit, row[key])
    return rows


def markdown(results):
    lines = ["# Stability CLOUD evidence: field distributions behind the IS Code 2008 checks",
             "",
             "Generated by `benchmarks/stability_cloud/run_cloud.py` (TensorCAD "
             f"{results['version']}). Every number below is read back from "
             "`artifacts/m4801/cloud_results.json` / `cloud_table.csv`; nothing is typed by "
             "hand. Field arrays: `artifacts/m4801/cloud_fields.json`.",
             "", "## Cases", ""]
    for name in ("barge", "barge_heavy", "hull"):
        field = results["fields"][name]
        lines.append(f"* **{name}** -- KG = {fmt(field['kg_m'])} m, GM0 = {fmt(field['gm0_m'])} m, "
                     f"KM = {fmt(field['km_m'])} m; {results['provenance'][name]}")
    lines += ["", "## Claims enforced by this run", "", "| claim | holds | detail |", "|---|---|---|"]
    for row in results["claims"]:
        lines.append(f"| {row['claim']} | {'PASS' if row['passed'] else 'FAIL'} | {row['detail']} |")
    lines += ["", "## IS Code 2008 Part A 2.2 criteria on the measured fields", "",
              "| case | criterion | quantity | required | measured | margin | verdict |",
              "|---|---|---|---|---|---|---|"]
    for name in ("barge", "barge_heavy", "hull"):
        for row in results["cases"][name]["criteria"]:
            lines.append(f"| {name} | {row['criterion']} | {row['quantity']} | {row['required']} | "
                         f"{fmt(row['measured'])} | {fmt(row['margin'])} | "
                         f"{'pass' if row['passed'] else '**fail**'} |")
    lines += ["", "## Cross-check against report #5216 (independent recomputation)", "",
              "| case | quantity | this run | #5216 | rel. difference |", "|---|---|---|---|---|"]
    for name, rows in results["cross_check_report_5216"].items():
        for quantity, row in rows.items():
            rel = row["relative_difference"]
            lines.append(f"| {name} | {quantity} | {fmt(row['this_run'])} | "
                         f"{fmt(row['report_5216'])} | "
                         f"{fmt(rel, 8) if rel is not None else 'n/a'} |")
    lines += ["", "## Field geometry (spatial evidence)", ""]
    for name in ("barge", "hull"):
        field = results["fields"][name]
        dens = field["density"]
        rel = (dens["it_integral_m4"] - field["it_exact_m4"]) / field["it_exact_m4"]
        lines.append(f"* {name}: I_T from the density field = {fmt(dens['it_integral_m4'])} m^4 vs "
                     f"mesh I_T = {fmt(field['it_exact_m4'])} m^4 ({100 * rel:+.3f} %); "
                     f"I_T from the station strip (2/3)y_hw^3 = "
                     f"{fmt(field['spatial']['it_from_stations_m4'])} m^4")
    lines += ["", "## Uncertainty budget", ""]
    for name, budget in results["uncertainty"].items():
        parts = ", ".join(f"u_{k} = {v:.3e}" for k, v in budget["components"].items())
        lines.append(f"* {name}: {parts} -> u_c = {budget['u']:.3e}, U95 = {budget['U95']:.3e}")
    for key, row in results["cases"]["hull"]["ladder"]["gci"].items():
        reason = row["fallback_reason"] or f"observed order {row['observed_order']:.2f}"
        lines.append(f"* hull mesh GCI [{key}]: fine = {fmt(row['fine'])}, medium = "
                     f"{fmt(row['medium'])}, coarse = {fmt(row['coarse'])}, U95 = "
                     f"{fmt(row['u95'], 8)} ({reason})")
    lines += ["* barge mesh term: not applicable -- the model is a polyhedron and the "
              "integral is exact for it (GCI needs a mesh family; fallback_reason recorded in "
              "the JSON), so U95 there carries only u_angle and u_solver.",
              "", "## Figures", "",
              "* `fig1_spatial_cloud.png` -- wetted breadth field 2y(x,z), I_T density field "
              "rho = y^2 and the station contribution.",
              "* `fig2_angular_cloud.png` -- GZ(phi, KG) field and criterion quantities vs KG.",
              "* `fig3_profiles_vs_reference.png` -- field slice vs closed form, angle-grid "
              "convergence, mesh ladder, criterion margins."]
    return "\n".join(lines) + "\n"


def figures(results, fields, out: Path):
    """Three figures; matplotlib is imported here so the JSON/CSV/Markdown never depend on it."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    hull, barge = fields["hull"], fields["barge"]
    xs = np.array(hull["spatial"]["x_m"])
    zs = np.array(hull["spatial"]["z_m"])

    fig, axes = plt.subplots(2, 2, figsize=(12.4, 8.4), constrained_layout=True)
    breadth = np.array(hull["spatial"]["half_breadth_rows_m"], dtype=float).T * 2.0
    img = axes[0, 0].pcolormesh(xs, zs, breadth, shading="auto", cmap="viridis")
    axes[0, 0].axhline(0.0, color="w", lw=1.4)
    axes[0, 0].text(0.4, -0.09, "waterline z = 0", color="w", fontsize=8)
    axes[0, 0].set_title(f"(a) hull wetted breadth field 2y(x,z) [m]\nL=12, B=2.4, T=0.9 m, "
                         f"{xs.size} x {zs.size} samples", fontsize=9)
    axes[0, 0].set_xlabel("x [m]")
    axes[0, 0].set_ylabel("z [m]  (0 = design waterline)")
    fig.colorbar(img, ax=axes[0, 0], label="2y [m]")
    for pos, label in enumerate(("hull", "barge")):
        case = fields[label]
        rho = np.array(case["density"]["density_y2_m2"], dtype=float).T
        ax = axes[0, 1] if pos == 0 else axes[1, 0]
        img = ax.pcolormesh(np.array(case["density"]["x_m"]), np.array(case["density"]["y_m"]),
                            rho, shading="auto", cmap="magma")
        rel = 100.0 * (case["density"]["it_integral_m4"] - case["it_exact_m4"]) / case["it_exact_m4"]
        ax.set_title(f"({'b' if pos == 0 else 'c'}) {label} waterplane density rho = y^2 [m^2]\n"
                     f"field integral {case['density']['it_integral_m4']:.5f} vs mesh I_T "
                     f"{case['it_exact_m4']:.5f} m^4 ({rel:+.3f} %)", fontsize=9)
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
        fig.colorbar(img, ax=ax, label="y^2 [m^2]")
    ax = axes[1, 1]
    for label, colour in (("hull", "#1f77b4"), ("barge", "#d62728")):
        curve = np.array(fields[label]["spatial"]["it_contribution_m4_per_m"], dtype=float)
        ax.plot(fields[label]["spatial"]["x_m"], curve / curve.max(), color=colour, label=label)
    ax.set_title("(d) station contribution (2/3) y_hw^3 [m^4/m], normalised\n"
                 "where the metacentric moment BM = I_T / V is carried", fontsize=9)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("normalised [-]")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.suptitle(f"CLOUD evidence 1 - spatial distribution behind BM = I_T/V "
                 f"(TensorCAD {__version__})", fontsize=11)
    fig.savefig(out / "fig1_spatial_cloud.png", dpi=110)
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(12.4, 8.4), constrained_layout=True)
    for row, (label, kg_d, kg_h) in enumerate((("barge", KG_BARGE, KG_HEAVY), ("hull", None, None))):
        case = fields[label]
        heels = np.array(case["gz_field"]["heel_deg"])
        kgs = np.array(case["gz_field"]["kg_m"])
        vals = np.array(case["gz_field"]["gz_m"], dtype=float)
        ax = axes[row, 0]
        img = ax.pcolormesh(heels, kgs, vals, shading="auto", cmap="RdBu_r", vmin=-0.2, vmax=0.6)
        cont = ax.contour(heels, kgs, vals, levels=[-0.05, 0.0, 0.05, 0.20, 0.40],
                          colors="k", linewidths=0.8)
        ax.clabel(cont, fmt="%.2f", fontsize=7)
        ax.axhline(case["km_m"], color="k", ls=":", lw=1.2)
        ax.text(55.0, case["km_m"] + 0.02, f"GM0=0 at KG=KM={case['km_m']:.4f} m", fontsize=7)
        if kg_d:
            ax.axhline(kg_d, color="#2ca02c", lw=1.3, label=f"design KG={kg_d:.2f} m")
            ax.axhline(kg_h, color="#ff7f0e", lw=1.3, ls="--",
                       label=f"heavy KG={kg_h:.2f} m (negative control)")
            ax.legend(fontsize=7, loc="lower right")
        ax.set_title(f"({'a' if row == 0 else 'c'}) {label}: GZ(phi,KG) [m]\n"
                     f"iso-lines include 0.20 m (2.2.4); KG {KG_GRID[0]:.3f}..{KG_GRID[-1]:.3f} m, "
                     f"heel step 1 deg", fontsize=9)
        ax.set_xlabel("heel phi [deg]")
        ax.set_ylabel("KG above keel [m]")
        fig.colorbar(img, ax=ax, label="GZ [m]")
        sweep = case["loading_sweep"]
        ax = axes[row, 1]
        ax.plot(sweep["kg_m"], sweep["A030"], color="#1f77b4", label="A(0,30) [m.rad]")
        ax.plot(sweep["kg_m"], sweep["A040"], color="#17becf", label="A(0,40) [m.rad]")
        ax.plot(sweep["kg_m"], sweep["A3040"], color="#9467bd", label="A(30,40) [m.rad]")
        for value, colour in ((0.055, "#1f77b4"), (0.090, "#17becf"), (0.030, "#9467bd")):
            ax.axhline(value, color=colour, ls=":", lw=1.0)
        ax.set_xlabel("KG above keel [m]")
        ax.set_ylabel("criterion area [m.rad]")
        twin = ax.twinx()
        twin.plot(sweep["kg_m"], sweep["GZ30"], color="#d62728", label="GZ(30) [m]")
        twin.axhline(0.20, color="#d62728", ls=":", lw=1.0)
        twin.set_ylabel("GZ(30) [m]")
        ax.axvline(case["kg_m"], color="k", lw=1.0, ls="-.")
        ax.set_title(f"({'b' if row == 0 else 'd'}) {label}: criterion quantities vs KG\n"
                     "dotted = 2.2.1/2.2.2/2.2.3/2.2.4 thresholds; dash-dot = evaluated "
                     f"KG={case['kg_m']:.4f} m", fontsize=9)
        handles = ax.get_lines()[:3] + twin.get_lines()[:1]
        ax.legend(handles, [h.get_label() for h in handles], fontsize=7, loc="upper right")
    fig.suptitle(f"CLOUD evidence 2 - angular/loading distribution of GZ and criterion areas "
                 f"(TensorCAD {__version__})", fontsize=11)
    fig.savefig(out / "fig2_angular_cloud.png", dpi=110)
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(12.4, 8.4), constrained_layout=True)
    bcase = results["cases"]["barge"]
    heels = np.array(barge["heel_deg"])
    gz = np.array(barge["gz_m_at_kg"], dtype=float)
    closed = np.array(barge["gz_reference_m"], dtype=float)
    keep = [i for i, h in enumerate(barge["heel_deg"]) if h <= bcase["reference"]["validity_deg"]]
    ax = axes[0, 0]
    ax.plot(heels, gz, color="#1f77b4", lw=1.6, label="mesh integral (righting_arm)")
    ax.plot(heels[keep], closed[keep], "k--", lw=1.0, label="wall-sided closed form")
    ax.axvline(bcase["reference"]["validity_deg"], color="grey", ls=":", lw=1.0)
    ax.fill_between(heels, 0.0, gz, where=heels <= 30.0, color="#1f77b4", alpha=0.12)
    ax.fill_between(heels, 0.0, gz, where=(heels >= 30) & (heels <= 40), color="#ff7f0e", alpha=0.18)
    ax.axhline(0.20, color="#d62728", ls=":", lw=1.0)
    ax.text(41.0, 0.205, "2.2.4: GZ(30) >= 0.20 m", fontsize=7, color="#d62728")
    ax.text(6.0, 0.02, "A(0,30)", fontsize=8)
    ax.text(32.0, 0.02, "A(30,40)", fontsize=8)
    ax.set_xlabel("heel phi [deg]")
    ax.set_ylabel("GZ [m]")
    ax.set_title(f"(a) barge KG=1.20 m: field slice vs closed form\nmax |difference| = "
                 f"{bcase['curve_vs_closed_form']['worst_abs_m']:.2e} m inside the validity window",
                 fontsize=9)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax = axes[0, 1]
    steps = np.array([r["step_deg"] for r in bcase["area_ladder"]])
    errs = np.array([abs(r["a030_relative_error"]) for r in bcase["area_ladder"]])
    ax.loglog(steps, errs, "o-", color="#1f77b4", label="|A(0,30) rel. error| vs angle step")
    ax.loglog(steps, errs[0] * (steps / steps[0]) ** 2, "k--", lw=1.0, label="h^2 reference")
    ax.set_xlabel("angle step h [deg]")
    ax.set_ylabel("relative error of A(0,30) [-]")
    ax.set_title(f"(b) angle-grid convergence of the criterion integral\nfitted order p = "
                 f"{bcase['area_order']['order']:.2f} (2 deg and finer)", fontsize=9)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")
    ax = axes[1, 0]
    hheels = np.array(hull["heel_deg"], dtype=float)
    ax.plot(hheels, np.array(hull["gz_m_at_kg"], dtype=float), color="#2ca02c", lw=1.6,
            label=f"hull {HULL_NX}x{HULL_NZ}")
    for row in results["cases"]["hull"]["ladder"]["meshes"]:
        hs = [float(k) for k in row["gz_by_heel"]]
        ax.plot(hs, [row["gz_by_heel"][k] for k in row["gz_by_heel"]], lw=0.9, alpha=0.85,
                label=f"mesh {row['nx']}x{row['nz']}")
    extrap = results["cases"]["hull"]["ladder"]["gz_extrapolated"]
    ax.plot([r["heel_deg"] for r in extrap], [r["value"] for r in extrap], "k+", ms=9,
            label="Richardson-extrapolated (3 meshes)")
    ax.set_xlabel("heel phi [deg]")
    ax.set_ylabel("GZ [m]")
    ax.set_title("(c) hull: mesh ladder of the GZ field and its extrapolation\n"
                 "no closed form exists for this hull, so the reference is the extrapolation",
                 fontsize=9)
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)
    ax = axes[1, 1]
    labels, ratios, colours = [], [], []
    for case_label in ("barge", "hull"):
        for row in results["cases"][case_label]["criteria"]:
            if row["criterion"].endswith("alt)"):
                continue
            labels.append(f"{case_label} {row['criterion']}")
            ratios.append(row["measured"] / row["required"])
            colours.append("#2ca02c" if row["passed"] else "#d62728")
    pos = np.arange(len(labels))
    ax.barh(pos, ratios, color=colours)
    ax.axvline(1.0, color="k", lw=1.2)
    ax.set_yticks(pos)
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("measured / required [-]")
    ax.set_title("(d) criterion margins (1.0 = the IS Code threshold)\n"
                 "green = pass, red = fail; the fields above explain each value", fontsize=9)
    ax.grid(alpha=0.3, axis="x")
    fig.suptitle(f"CLOUD evidence 3 - field slices against reference values "
                 f"(TensorCAD {__version__})", fontsize=11)
    fig.savefig(out / "fig3_profiles_vs_reference.png", dpi=110)
    plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("artifacts/m4801"))
    parser.add_argument("--docs", type=Path, default=Path("docs/stability_cloud.md"))
    parser.add_argument("--no-figures", action="store_true")
    args = parser.parse_args(argv)
    args.output.mkdir(parents=True, exist_ok=True)
    results, fields = build()
    (args.output / "cloud_fields.json").write_text(json.dumps(fields, indent=1))
    (args.output / "cloud_results.json").write_text(json.dumps(results, indent=1))
    rows = table_rows(results)
    with (args.output / "cloud_table.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    if not args.no_figures:
        figures(results, fields, args.output)
    args.docs.parent.mkdir(parents=True, exist_ok=True)
    args.docs.write_text(markdown(results))
    failed = [row["claim"] for row in results["claims"] if not row["passed"]]
    print(f"claims: {len(results['claims'])} total, {len(failed)} failed")
    print(f"table rows: {len(rows)}")
    for name in failed:
        print("FAILED CLAIM:", name)
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())