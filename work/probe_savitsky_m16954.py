"""Probe: does tensorcad.savitsky actually solve anything besides trim?

Read-only diagnostic used by the m16954 audit. Prints raw numbers, no claims.
"""
import math

from tensorcad.savitsky import (
    lift_coefficient,
    planing_prediction,
    resistance_components,
    running_trim,
)

RHO = 1025.5
G = 9.80665

print("A. running_trim over a (weight, speed, beam, deadrise) matrix")
print(f"{'weight_N':>10} {'V':>6} {'B':>5} {'beta':>5} {'tau_deg':>8} {'lambda':>8} {'CL_need':>10} {'CL_got':>10} {'resid%':>9}")
for w in (5e3, 1e4, 2e4, 4e4):
    for v in (6.0, 10.0, 15.0):
        for b in (2.0, 3.0):
            for beta in (10.0, 20.0, 30.0):
                try:
                    tau, lam = running_trim(w, v, b, beta, water_density=RHO)
                    cl_need = w / (0.5 * RHO * v * v * b * b)
                    cl_got = lift_coefficient(tau, lam, beta)
                    resid = 100.0 * (cl_got - cl_need) / cl_need
                    print(f"{w:10.0f} {v:6.1f} {b:5.1f} {beta:5.0f} {tau:8.3f} {lam:8.3f} "
                          f"{cl_need:10.4f} {cl_got:10.4f} {resid:9.2f}")
                except Exception as exc:  # noqa: BLE001
                    print(f"{w:10.0f} {v:6.1f} {b:5.1f} {beta:5.0f}  FAILED: {type(exc).__name__}: {exc}")

print()
print("B. planeing_prediction: does lambda move with weight at fixed speed/geometry?")
for w in (8e3, 1.2e4, 1.6e4, 2.4e4):
    try:
        p = planing_prediction(w, 12.0, 8.0, 2.5, 20.0)
        print(f"  W={w:8.0f} N  trim={p.trim_deg:7.3f} deg  lambda={p.wetted_length_ratio:6.3f} "
              f"Lw={p.wetted_length_m:6.3f} m  Rt={p.total_drag_n:9.1f} N  Pe={p.effective_power_w/1000:7.2f} kW")
    except Exception as exc:  # noqa: BLE001
        print(f"  W={w:8.0f} N  FAILED: {type(exc).__name__}: {exc}")

print()
print("C. lambda as a function of speed (same hull, same weight)")
for v in (6.0, 8.0, 10.0, 13.0, 16.0, 20.0):
    try:
        p = planing_prediction(1.5e4, v, 8.0, 2.5, 20.0)
        print(f"  V={v:5.1f} m/s  Fn={p.froude_number:5.3f}  trim={p.trim_deg:7.3f}  "
              f"lambda={p.wetted_length_ratio:6.3f}  CL={p.lift_coefficient:7.4f}  Rt={p.total_drag_n:9.1f} N")
    except Exception as exc:  # noqa: BLE001
        print(f"  V={v:5.1f} m/s  FAILED: {type(exc).__name__}: {exc}")

print()
print("D. lift vs displacement closure in planing_prediction (should be ~0%)")
for v in (8.0, 12.0, 16.0):
    p = planing_prediction(1.5e4, v, 8.0, 2.5, 20.0)
    q = 0.5 * RHO * v * v * 2.5 * 2.5
    print(f"  V={v:5.1f}  CL*0.5 rho V^2 B^2 = {p.lift_coefficient*q:9.1f} N vs W={p.weight_n:8.0f} N")

print()
print("E. pressure-drag term sanity: 0.008/sqrt(tau) at small tau")
for tau in (1.0, 2.0, 4.0, 8.0):
    print(f"  tau={tau:4.1f} deg  tan(tau)={math.tan(math.radians(tau)):8.5f}  "
          f"0.008/sqrt(tau)={0.008/math.sqrt(tau):8.5f}  "
          f"ratio={0.008/math.sqrt(tau)/math.tan(math.radians(tau)):8.4f}")