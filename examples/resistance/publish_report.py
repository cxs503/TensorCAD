#!/usr/bin/env python3
"""Build and publish the m13816 platform report (report content injected from artifacts)."""
import base64
import json
import os
import pathlib
import urllib.request

BASE = os.environ["TENSOR_API_BASE"]
USER = os.environ["TENSOR_AGENT_USER"]
SECRET = os.environ["TENSOR_AGENT_PASS"]


def post(path, payload, token=None):
    data = json.dumps(payload).encode()
    headers = {"Content-Type": "application/json"}
    if token:
        headers["Authorization"] = f"Bearer {token}"
    request = urllib.request.Request(BASE + path, data=data, headers=headers, method="POST")
    with urllib.request.urlopen(request, timeout=60) as response:
        return json.loads(response.read().decode())


auth = post("/auth/login", {"username": USER, "password": SECRET})
token = auth.get("access_token") or auth.get("token")

r = json.load(open("artifacts/resistance/results.json"))
fr, wa, pr, sc = r["friction"], r["wetted_area"], r["prohaska"], r["scaling"]
pm, ps = sc["prediction"]["model"], sc["prediction"]["ship"]
md = f"""# TensorCAD resistance: ITTC-57 / Prohaska / ITTC-78, measured

## 1. Case introduction
Geometry could not answer "how hard does this hull pull". `tensorcad/resistance.py` adds the
1966 ITTC first stage (friction line + Prohaska form factor) and the 1978 two-dimensional
model-to-ship scaling. Case: Wigley reference hull L=10 m, B=1 m, T=0.625 m (ITTC proportions),
sea water at 15 C (rho=1025.9 kg/m3, nu=1.1892e-6 m2/s); ship = geometrically similar L=100 m.
Reference values: Schoenherr (1932) implicit flat-plate line, the Prandtl-Schlichting line, and
`tensorcad.wigley.analytic_wetted_surface_area` (composite Gauss-Legendre on the analytic
surface, never on the mesh).

## 2. Method
- `ittc57_friction_coefficient`: CF = 0.075/(log10 Re - 2)^2.
- `schoenherr_friction_coefficient`: root of 0.242/sqrt(CF) = log10(Re CF) by bisection (tol 1e-13).
- `prohaska_fit`: CT/CF = (1+k) + y Fn^4/CF by least squares, with residuals and standard errors;
  needs >= 3 points; either friction coefficients or a model length, never both.
- `bare_hull_resistance`: R = 0.5 rho S V^2 x {{CF, (1+k)CF, CR, CA}}; CR and 1+k have no invented
  defaults (the defaults read as "not measured").
- `ship_from_model_test`: V_s = V_m sqrt(L_s/L_m); CR = CT_m - (1+k)CF_m; CT_s = (1+k)CFS + CR + CA.
- Wetted surface from the mesh: shell area minus the waterline deck patch, cross-checked against the
  analytic flank integral on 41x13 / 81x23 / 161x43 grids.

## 3. Results (read from artifacts/resistance/results.json)
ITTC-57 vs Schoenherr worst |diff| = {fr['worst_relative_difference_schoenherr_ship_range']*100:.2f} % for Re>=1e7
({fr['worst_relative_difference_schoenherr_full_range']*100:.2f} % over the whole 1e4..1e9 grid);
vs Prandtl-Schlichting {fr['worst_relative_difference_prandtl_schlichting_ship_range']*100:.2f} % for Re>=1e7.
Worst Schoenherr equation residual {fr['worst_absolute_schoenherr_residual']:.2e}.

| quantity | model (10 m) | ship (100 m) |
|---|---:|---:|
| Re | {pm['reynolds']:.4e} | {ps['reynolds']:.4e} |
| CF | {pm['friction_coefficient']:.4e} | {ps['friction_coefficient']:.4e} |
| (1+k)CF | {pm['viscous_coefficient']:.4e} | {ps['viscous_coefficient']:.4e} |
| CR | {sc['prediction']['residual_coefficient']:.4e} | {sc['prediction']['residual_coefficient']:.4e} |
| CT | {pm['total_coefficient']:.4e} | {ps['total_coefficient']:.4e} |
| R_T | {pm['total_resistance_n']:.2f} N | {ps['total_resistance_n']/1000:.2f} kN |
| P_E | -- | {ps['effective_power_w']/1e6:.3f} MW |

Wetted area (mesh vs analytic flank): finest grid relative error {wa['finest_relative_error']*100:.4f} %;
deck share of shell {wa['worst_deck_share']*100:.2f} % (using shell as wetted inflates friction 44.8 %).
Prohaska: injected 1+k = {pr['injected_form_factor']} recovered {pr['recovered_form_factor']:.12f};
injected y = {pr['injected_wave_coefficient']} recovered {pr['recovered_wave_coefficient']:.4e};
line residual {pr['max_relative_residual']:.1e}.
Scaling: R_s/R_m = {sc['prediction']['ship_resistance_ratio']:.1f} against the naive lambda^3 = {sc['naive_lambda_cubed_ratio']:.0f};
scale effect CT_m/CT_s = {ps and sc['prediction']['scale_effect']:.4f}; two independent routes to R_T agree to {sc['routes_agree_relative']:.1e} relative.

## 4. Figures
Left: the three flat-plate lines over Re = 1e4..1e9 (log-log). Right: the Prohaska plot CT/CF against
Fn^4/CF for the synthetic 10 m model test with the fitted line whose intercept is the form factor.

## 5. Notes / limitations
- Residual (wave) resistance is not modelled here; CR and CA are caller inputs, and a negative CR on a
  well-formed hull at low Fn is a real outcome that is reported rather than clipped.
- The 1957 line does not agree with Schoenherr below Re ~ 1e6 (16.1 % at 1e5, 43.4 % at 1e4); that is a
  property of the published lines and the tests assert the divergence instead of hiding it.
- A closed hull's shell includes the waterline deck patch (30.9 % here); `mesh_bare_hull_resistance` takes
  the mesh as given and documents which meshes need no correction.
- Prohaska assumes R_W ~ Fn^4, so the fitted 1+k is valid only below the hump (roughly Fn < 0.4).
- Fluid constants are conventional 15 C inputs, not verified against an ITTC document in this repo; every
  claim above is a ratio of like quantities and independent of them.
- Also fixed: the pinned public-surface contract test had been red at HEAD for ten exported names.

## 6. Error judgement vs reference values
| key quantity | reference | error | under 3 %? |
|---|---|---:|---|
| CF ITTC-57, Re>=1e7 | Schoenherr line | {fr['worst_relative_difference_schoenherr_ship_range']*100:.2f} % | yes |
| CF ITTC-57, Re>=1e7 | Prandtl-Schlichting line | {fr['worst_relative_difference_prandtl_schlichting_ship_range']*100:.2f} % | yes |
| CF ITTC-57, Re=1e5 | Schoenherr line | 16.07 % | no (stated: lines not interchangeable at low Re) |
| wetted area 161x43 | analytic flank integral | {abs(wa['finest_relative_error'])*100:.4f} % | yes |
| wetted area 41x13 | analytic flank integral | 0.0201 % | yes |
| 1+k (Prohaska) | injected 1.15 | {pr['relative_error_form_factor']*100:.1e} % | yes |
| ship R_T, model-to-ship route | direct construction with the same CR | {sc['routes_agree_relative']*100:.1e} % | yes |

Commit 34835af. Artifacts: artifacts/resistance/{{results.json, plot.png, report.md}}; docs/resistance.md;
tests/test_resistance.py (60 cases, all pass); examples/resistance/run_example.py exits non-zero when a
claim fails.
"""
png = base64.b64encode(pathlib.Path("artifacts/resistance/plot.png").read_bytes()).decode()
payload = {
    "title": "TensorCAD resistance: ITTC-57 friction, Prohaska form factor, ITTC-78 model-to-ship scaling (measured)",
    "summary": (
        "New tensorcad.resistance closing the documented resistance gap: ITTC-57 vs Schoenherr within 2.24 % "
        "in the ship range, Prohaska 1+k recovered to round-off, mesh wetted area within 0.0016 % of the "
        "analytic flank, 10 m -> 100 m Wigley scaling R_s/R_m = 615.5 rather than 1000."
    ),
    "scenario_key": "tensorcad_resistance_ittc",
    "mission_id": 13816,
    "tags": ["TensorCAD", "resistance", "ITTC-57", "Prohaska", "ITTC-78", "verification", "m13816"],
    "markdown": md,
    "figures": [
        {
            "name": "resistance_lines_prohaska.png",
            "caption": (
                "Left: ITTC-57, Schoenherr and Prandtl-Schlichting friction lines over Re = 1e4..1e9. "
                "Right: Prohaska plot of a synthetic 10 m model test with the fitted line; the intercept "
                "is the form factor."
            ),
            "png_base64": png,
        }
    ],
}
pathlib.Path("/tmp/report.json").write_text(json.dumps(payload))
answer = post("/tensor/reports", payload, token)
print(json.dumps(answer)[:800])