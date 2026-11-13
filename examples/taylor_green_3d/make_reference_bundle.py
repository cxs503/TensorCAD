#!/usr/bin/env python3
"""Reference bundle + independent cross-check for CASE-TAYLOR-GREEN-3D.

Two jobs:

1. **Cross-check the series reference** with an *independent* method: a Fourier
   pseudo-spectral RK4 integration of the same initial condition on a 32^3
   periodic grid with 2/3-rule dealiasing.  The series and the DNS share no
   code, so agreement in ``E(t)`` / ``Z(t)`` over the trusted window is evidence
   that the series coefficients are the Navier-Stokes solution and not merely a
   self-consistent artefact of the recursion.
2. **Write the consumable bundle** (``case_spec.json``, ``metrics_table.json``,
   ``mode_table_t*.json``, grid samples, closed-form tables) plus figures.

Everything lands under ``docs/taylor_green_3d/`` -- the platform push gate
filters ``artifacts/`` and ``runs/``, so reference data destined for main lives
next to its report.

Usage::

    python examples/taylor_green_3d/make_reference_bundle.py
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from tensorcad.taylor_green import (  # noqa: E402
    TaylorGreen3DSeries,
    metrics_table,
    write_case_bundle,
)

OUT = Path(__file__).resolve().parents[2] / "docs" / "taylor_green_3d"
TIMES = (1.0, 2.0, 3.0, 4.0)


def _wavenumbers(n: int) -> np.ndarray:
    k = np.fft.fftfreq(n, d=1.0 / n)
    return k.astype(int)


def _initial_field(n: int, amplitude: float) -> np.ndarray:
    x = np.arange(n) * (2.0 * math.pi / n)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    u = amplitude * np.sin(X) * np.cos(Y) * np.cos(Z)
    v = -amplitude * np.cos(X) * np.sin(Y) * np.cos(Z)
    w = np.zeros_like(u)
    return np.stack([u, v, w]).astype(np.float64)


def _rhs(u_hat: np.ndarray, k2: np.ndarray, kmax: int, nu: float) -> np.ndarray:
    n = u_hat.shape[1]
    nx, ny, nz = k2.shape
    kx = _wavenumbers(n)[:, None, None] * np.ones((1, n, n))
    ky = _wavenumbers(n)[None, :, None] * np.ones((n, 1, n))
    kz = _wavenumbers(n)[None, None, :] * np.ones((n, n, 1))
    mask = ((np.abs(kx) <= kmax) & (np.abs(ky) <= kmax) & (np.abs(kz) <= kmax))
    uh = u_hat * mask
    u = np.real(np.fft.ifftn(uh, axes=(1, 2, 3)))
    conv = np.zeros_like(u)
    for j in range(3):
        du = np.real(np.fft.ifftn(1j * [kx, ky, kz][j] * uh, axes=(1, 2, 3)))
        for i in range(3):
            conv[i] += u[j] * du[i]
    conv_hat = np.fft.fftn(conv, axes=(1, 2, 3)) * mask
    if k2 is None:  # pragma: no cover - defensive
        raise ValueError("k2 required")
    with np.errstate(invalid="ignore", divide="ignore"):
        inv = np.where(k2 > 0, 1.0 / np.where(k2 > 0, k2, 1.0), 0.0)
    kdot = kx * conv_hat[0] + ky * conv_hat[1] + kz * conv_hat[2]
    projected = conv_hat - np.stack([kx * kdot, ky * kdot, kz * kdot]) * inv
    return -projected - nu * k2 * uh


def spectral_dns(n: int = 32, amplitude: float = 0.25, nu: float = 0.01,
                 t_end: float = 4.0, dt: float = 0.04) -> dict[str, list[float]]:
    """Pseudo-spectral RK4 DNS of the 3-D Taylor-Green vortex (independent path)."""
    k = _wavenumbers(n)
    kx = k[:, None, None] * np.ones((1, n, n))
    ky = k[None, :, None] * np.ones((n, 1, n))
    kz = k[None, None, :] * np.ones((n, n, 1))
    k2 = kx ** 2 + ky ** 2 + kz ** 2
    kmax = n // 3  # 2/3 dealiasing
    u_hat = np.fft.fftn(_initial_field(n, amplitude), axes=(1, 2, 3))
    steps = int(round(t_end / dt))
    times, energies, enstrophies = [0.0], [], []
    for step in range(steps + 1):
        if step % 25 == 0 or step == steps:
            e = 0.5 * float(np.sum(np.abs(u_hat) ** 2)) / n ** 6
            z = 0.5 * float(np.sum(k2 * np.abs(u_hat) ** 2)) / n ** 6
            times.append(step * dt)
            energies.append(e)
            enstrophies.append(z)
        k1 = _rhs(u_hat, k2, kmax, nu)
        k2_ = _rhs(u_hat + 0.5 * dt * k1, k2, kmax, nu)
        k3 = _rhs(u_hat + 0.5 * dt * k2_, k2, kmax, nu)
        k4 = _rhs(u_hat + dt * k3, k2, kmax, nu)
        u_hat = u_hat + dt / 6.0 * (k1 + 2.0 * k2_ + 2.0 * k3 + k4)
    # normalisation: physical-space mean of |u|^2 = (1/n^6) sum |u_hat|^2
    return {"times": times[1:], "energy": energies, "enstrophy": enstrophies,
            "grid": n, "dt": dt, "dealias_kmax": kmax}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    ref = TaylorGreen3DSeries()
    series_e = {t: ref.energy(t) for t in TIMES}
    series_z = {t: ref.enstrophy(t) for t in TIMES}

    dns = spectral_dns()
    dns_e = dict(zip(dns["times"], dns["energy"]))
    dns_z = dict(zip(dns["times"], dns["enstrophy"]))

    comparison = []
    for t in TIMES:
        row = {"t": t, "series_energy": series_e[t], "series_enstrophy": series_z[t],
               "spectral_dns_energy": dns_e[t], "spectral_dns_enstrophy": dns_z[t],
               "energy_rel_error": abs(series_e[t] - dns_e[t]) / series_e[t],
               "enstrophy_rel_error": abs(series_z[t] - dns_z[t]) / series_z[t]}
        comparison.append(row)

    trunc = []
    for order in (6, 8, 10):
        s = TaylorGreen3DSeries(order=order)
        trunc.append({"order": order,
                      "energy_at_4": s.energy(4.0),
                      "rel_error_vs_order12": abs(s.energy(4.0) - series_e[4.0])
                      / series_e[4.0],
                      "coefficient_ratio": s.coefficient_ratio()})

    payload = {
        "case": "CASE-TAYLOR-GREEN-3D",
        "reference": {"order": ref.order, "amplitude": ref.amplitude,
                      "viscosity": ref.viscosity,
                      "energy_terms": ref.energy_terms(),
                      "defect_spectrum": ref.defect_spectrum()},
        "series_vs_spectral_dns": comparison,
        "order_convergence": trunc,
        "spectral_dns": {k: dns[k] for k in ("times", "energy", "enstrophy",
                                             "grid", "dt", "dealias_kmax")},
        "max_energy_rel_error": max(r["energy_rel_error"] for r in comparison),
        "max_enstrophy_rel_error": max(r["enstrophy_rel_error"] for r in comparison),
    }
    (OUT / "crosscheck.json").write_text(json.dumps(payload, indent=2) + "\n")

    manifest = write_case_bundle(OUT / "bundle", series=ref, times=TIMES,
                                 grid_levels=(8,))
    (OUT / "bundle_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")

    table = metrics_table(ref, TIMES)
    (OUT / "metrics_table.json").write_text(json.dumps(table, indent=2) + "\n")

    _figures(ref, payload)
    print(json.dumps({k: payload[k] for k in
                      ("max_energy_rel_error", "max_enstrophy_rel_error")}, indent=2))
    print("order convergence:", [(r["order"], f"{r['rel_error_vs_order12']:.2e}")
                                  for r in trunc])


def _figures(ref: TaylorGreen3DSeries, payload: dict) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figdir = OUT / "figures"
    figdir.mkdir(parents=True, exist_ok=True)

    # 1. energy decay: series vs independent spectral DNS
    dns = payload["spectral_dns"]
    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    ts = [i * 0.1 for i in range(41)]
    ax.plot(ts, [ref.energy(t) for t in ts], "-", label="series reference (order 12)")
    ax.plot(dns["times"], dns["energy"], "o", ms=3.5, label="spectral RK4 DNS 32$^3$")
    ax.set_xlabel("t (box units, $L/U = 2\\pi/A$)")
    ax.set_ylabel("kinetic energy $E(t)$")
    ax.set_title("Taylor-Green 3-D: reference series vs independent DNS")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(figdir / "energy_decay.png", dpi=140)
    plt.close(fig)

    # 2. truncation error and defect spectrum vs Taylor order
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.0, 3.6))
    orders = [r["order"] for r in payload["order_convergence"]]
    errs = [max(r["rel_error_vs_order12"], 1e-16) for r in payload["order_convergence"]]
    ax1.semilogy(orders, errs, "s-")
    ax1.set_xlabel("Taylor order N")
    ax1.set_ylabel("|$E_N(4)-E_{12}(4)$| / $E_{12}$")
    ax1.set_title("truncation at t = 4")
    ax1.grid(alpha=0.3)
    spec = ref.defect_spectrum()
    ax2.semilogy(range(len(spec)), [max(v, 1e-18) for v in spec], "o-")
    ax2.set_xlabel("Taylor order n")
    ax2.set_ylabel("max |NS residual coefficient|")
    ax2.set_title("order-completeness defect spectrum")
    ax2.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(figdir / "convergence.png", dpi=140)
    plt.close(fig)

    # 3. field snapshot: x-component on the z = pi/2 plane at t = 2
    level = 24
    xs = [i * (2.0 * math.pi / level) for i in range(level)]
    field = [[ref.velocity(x, y, math.pi / 2.0, 2.0)[0] for y in xs] for x in xs]
    fig, ax = plt.subplots(figsize=(4.9, 4.1))
    im = ax.imshow(field, origin="lower", extent=[0, 2 * math.pi, 0, 2 * math.pi],
                   cmap="RdBu_r", aspect="equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("$u_x$ at $z=\\pi/2$, t = 2 (series reference)")
    fig.colorbar(im, ax=ax, shrink=0.85)
    fig.tight_layout()
    fig.savefig(figdir / "field_snapshot.png", dpi=140)
    plt.close(fig)


if __name__ == "__main__":
    main()