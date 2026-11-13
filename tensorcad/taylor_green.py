"""Taylor-Green vortex reference solutions (incompressible Navier-Stokes).

Why this module exists
----------------------
The case-pool entry ``CASE-TAYLOR-GREEN-3D`` asks for a solver-vs-reference
decay check with ``ref: 级数精确解`` (series exact solution) inside an 8 % band.
TensorCAD's role in the ecosystem is geometry **and reference data** -- the same
role it plays for :mod:`tensorcad.kovasznay` (FVM order study) and
:mod:`tensorcad.cavity_mesh` (Ghia benchmark); the flow solvers live in
TensorLBM / TensorFVM / TensorDG and *consume* the artefacts written here.  This
module therefore ships the reference side of the Taylor-Green case in three
tiers of increasing difficulty:

1. :class:`TaylorGreenMode` -- the ``z``-invariant Taylor-Green mode.  Because
   ``lap(psi) = -2 psi`` the convective term is a pure gradient, so this is an
   **exact** unsteady Navier-Stokes solution in closed form
   (``u ~ exp(-2 nu t)``).  This is the tier a lattice-Boltzmann case can hit to
   a few tenths of a percent.
2. :class:`BeltramiABC` -- the ABC Beltrami flow with ``curl u = k u``.  Again
   exact (``u ~ exp(-k^2 nu t)``) but **z-dependent**, so it exercises the
   three-dimensional branch of a solver instead of a 2-D plane.
3. :class:`TaylorGreen3DSeries` -- the classic, genuinely three-dimensional
   Taylor-Green vortex.  Its convective term is *not* curl-free and therefore
   has **no closed form**; the reference is the time-Taylor series
   ``u(x,t) = sum_n u_n(x) t^n`` with coefficients generated from the projected
   Navier-Stokes recursion in a Fourier basis.  The truncation error is
   ``O(t^(order+1))`` and the coefficient decay is reported, so a consumer can
   pick a trustworthy ``t`` window.

Standard library only (``math``): TensorCAD keeps its zero-dependency guarantee
and every reference number is *recomputed* here rather than transcribed from a
table (zero hand-copying).

Verification surface (all machine-checkable)
--------------------------------------------
* :meth:`TaylorGreenMode.ns_residual` / :meth:`BeltramiABC.ns_residual` evaluate
  the momentum residual with analytic derivatives and give ``~1e-16``.
* :meth:`BeltramiABC.curl_eigen_residual` checks ``curl u = k u`` itself.
* :meth:`TaylorGreen3DSeries.divergence` is ``~1e-16``, and
  :meth:`TaylorGreen3DSeries.defect_order` shows the residual of the truncated
  series first appearing at ``t^order`` (maximal accuracy of the recursion).
* :func:`fit_decay_rate` recovers the analytic rate ``2 nu |k|^2`` from an
  energy history; :func:`compare_metrics` turns a solver field into relative
  errors against the reference -- exactly the quantities the case pool scores.

Reference (background, not transcribed)
---------------------------------------
Taylor, G. I. & Green, A. E. (1937). "Mechanism of the production of small
eddies from large ones." *Proc. R. Soc. Lond. A* 158(895), 499-521.
The time-Taylor-series route to a vorticity reference is the classical
Morf-Orszag-Frisch / Brachet family of papers; their tables are **not** copied
here -- the coefficients used below are generated and self-checked by code.
"""
from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Any, Mapping, Sequence

__all__ = [
    "TAYLOR_GREEN_PROVENANCE",
    "CASE_TAYLOR_GREEN_3D_ACCEPTANCE",
    "BeltramiABC",
    "TaylorGreen3DSeries",
    "TaylorGreenMode",
    "case_spec",
    "compare_metrics",
    "fit_decay_rate",
    "metrics_table",
    "mode_table",
    "sample_grid",
    "write_case_bundle",
]

TAYLOR_GREEN_PROVENANCE: dict[str, Any] = {
    "citation": (
        'Taylor, G. I. & Green, A. E. (1937). "Mechanism of the production of '
        'small eddies from large ones." Proc. R. Soc. Lond. A 158(895), 499-521.'
    ),
    "best_for": (
        "unsteady incompressible Navier-Stokes verification: exact z-invariant "
        "decay, exact z-dependent Beltrami decay, and a time-Taylor series for "
        "the full 3-D vortex (mode/energy/enstrophy decay references)"
    ),
    "series_method": (
        "time-Taylor series of the 3-D vortex in a projected Fourier basis "
        "(Morf-Orszag-Frisch family); coefficients generated and residual-"
        "checked here, no literature table transcribed"
    ),
}

CASE_TAYLOR_GREEN_3D_ACCEPTANCE: dict[str, Any] = {
    "case": "CASE-TAYLOR-GREEN-3D",
    "reference": "series exact solution / closed-form Fourier mode (this module)",
    "metrics": ["mode amplitude decay", "kinetic energy decay", "enstrophy"],
    "relative_tolerance": 0.08,
    "target_band": "<8% on every metric",
}

_TWO_PI = 2.0 * math.pi
Mode = tuple[int, int, int]
Vec = tuple[complex, complex, complex]


def _add(dst: dict[Mode, Vec], key: Mode, vec: Vec, scale: float = 1.0) -> None:
    cur = dst.get(key)
    if cur is None:
        dst[key] = (vec[0] * scale, vec[1] * scale, vec[2] * scale)
    else:
        dst[key] = (cur[0] + vec[0] * scale,
                    cur[1] + vec[1] * scale,
                    cur[2] + vec[2] * scale)


def _project_div_free(modes: Mapping[Mode, Vec]) -> dict[Mode, Vec]:
    """Leray projection onto the divergence-free subspace."""
    out: dict[Mode, Vec] = {}
    for key, vec in modes.items():
        kk = key[0] * key[0] + key[1] * key[1] + key[2] * key[2]
        if kk == 0:
            continue  # the mean mode carries no velocity in this family
        kd = key[0] * vec[0] + key[1] * vec[1] + key[2] * vec[2]
        out[key] = (vec[0] - key[0] * kd / kk,
                    vec[1] - key[1] * kd / kk,
                    vec[2] - key[2] * kd / kk)
    return out


def _as_list(modes: Mapping[Mode, Vec]) -> list[tuple[int, int, int, complex, complex, complex]]:
    """Flatten a mode dict once so the convolution loop stays tight."""
    return [(k[0], k[1], k[2], v[0], v[1], v[2]) for k, v in modes.items()]


def _convective_list(p: Sequence[tuple[int, int, int, complex, complex, complex]],
                     q: Sequence[tuple[int, int, int, complex, complex, complex]],
                     kmax: int) -> dict[Mode, Vec]:
    """Fourier coefficients of ``(p . grad) q`` cut at ``|k|_inf`` <= kmax."""
    out: dict[Mode, Vec] = {}
    for kx1, ky1, kz1, ax, ay, az in p:
        for kx2, ky2, kz2, bx, by, bz in q:
            kx = kx1 + kx2
            if kx > kmax or kx < -kmax:
                continue
            ky = ky1 + ky2
            if ky > kmax or ky < -kmax:
                continue
            kz = kz1 + kz2
            if kz > kmax or kz < -kmax:
                continue
            s = 1j * (ax * kx2 + ay * ky2 + az * kz2)
            key = (kx, ky, kz)
            cur = out.get(key)
            t0, t1, t2 = s * bx, s * by, s * bz
            if cur is None:
                out[key] = (t0, t1, t2)
            else:
                out[key] = (cur[0] + t0, cur[1] + t1, cur[2] + t2)
    return out


def _convective(p: Mapping[Mode, Vec], q: Mapping[Mode, Vec],
                kmax: int) -> dict[Mode, Vec]:
    """Convenience wrapper around :func:`_convective_list`."""
    return _convective_list(_as_list(p), _as_list(q), kmax)


def _eval_modes(modes: Mapping[Mode, Vec], x: float, y: float,
                z: float) -> tuple[float, float, float]:
    out = [0.0, 0.0, 0.0]
    for (kx, ky, kz), vec in modes.items():
        phase = kx * x + ky * y + kz * z
        c, s = math.cos(phase), math.sin(phase)
        for i in range(3):
            a, b = vec[i].real, vec[i].imag
            out[i] += a * c - b * s
    return (out[0], out[1], out[2])


def _scalar_field(factors: Sequence[tuple[str, int]]) -> dict[Mode, complex]:
    """Product of Cartesian ``sin``/``cos`` factors as complex Fourier modes."""
    out: dict[Mode, complex] = {(0, 0, 0): 1.0 + 0.0j}
    for kind, axis in factors:
        key_fwd, key_bwd = [0, 0, 0], [0, 0, 0]
        key_fwd[axis], key_bwd[axis] = 1, -1
        if kind == "sin":
            coeffs = ((tuple(key_fwd), -0.5j), (tuple(key_bwd), 0.5j))
        elif kind == "cos":
            coeffs = ((tuple(key_fwd), 0.5 + 0.0j), (tuple(key_bwd), 0.5 + 0.0j))
        else:  # pragma: no cover - guard
            raise ValueError(f"unknown factor {kind!r}")
        nxt: dict[Mode, complex] = {}
        for k1, v1 in out.items():
            for k2, v2 in coeffs:  # type: ignore[assignment]
                key = (k1[0] + k2[0], k1[1] + k2[1], k1[2] + k2[2])
                nxt[key] = nxt.get(key, 0.0j) + v1 * v2
        out = nxt
    return out# ---------------------------------------------------------------------------
# 1. z-invariant Taylor-Green mode: closed-form exact NS solution
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class TaylorGreenMode:
    """Exact decaying ``z``-invariant Taylor-Green mode on ``[0, 2*pi]^3``.

    ``u = A e^{-2 nu t} sin x cos y``, ``v = -A e^{-2 nu t} cos x sin y``,
    ``p = +(A^2/4) e^{-4 nu t}(cos 2x + cos 2y)``.
    """

    amplitude: float = 1.0
    viscosity: float = 0.01
    pressure_datum: float = 0.0

    def decay(self, t: float) -> float:
        return math.exp(-2.0 * self.viscosity * t)

    def velocity(self, x: float, y: float, z: float,
                 t: float) -> tuple[float, float, float]:
        d = self.amplitude * self.decay(t)
        return (d * math.sin(x) * math.cos(y),
                -d * math.cos(x) * math.sin(y),
                0.0)

    def pressure(self, x: float, y: float, z: float, t: float) -> float:
        quad = 0.25 * self.amplitude ** 2 * math.exp(-4.0 * self.viscosity * t)
        return self.pressure_datum + quad * (math.cos(2.0 * x) + math.cos(2.0 * y))

    def velocity_gradient(self, x: float, y: float, z: float,
                          t: float) -> tuple[tuple[float, float, float], ...]:
        """``g[j][i] = d u_i / d x_j``."""
        d = self.amplitude * self.decay(t)
        cx, sx, cy, sy = math.cos(x), math.sin(x), math.cos(y), math.sin(y)
        return (
            (d * cx * cy, d * sx * sy, 0.0),
            (-d * sx * sy, -d * cx * cy, 0.0),
            (0.0, 0.0, 0.0),
        )

    def vorticity(self, x: float, y: float, z: float,
                  t: float) -> tuple[float, float, float]:
        d = self.amplitude * self.decay(t)
        return (0.0, 0.0, 2.0 * d * math.sin(x) * math.sin(y))

    def ns_residual(self, x: float, y: float, z: float,
                    t: float) -> tuple[float, float, float]:
        """Momentum residual ``du/dt + (u.grad)u + grad p - nu lap u``."""
        nu = self.viscosity
        u = self.velocity(x, y, z, t)
        g = self.velocity_gradient(x, y, z, t)
        d = self.amplitude * self.decay(t)
        dudt = (-2.0 * nu * d * math.sin(x) * math.cos(y),
                2.0 * nu * d * math.cos(x) * math.sin(y),
                0.0)
        conv = tuple(sum(u[j] * g[j][i] for j in range(3)) for i in range(3))
        quad = 0.25 * self.amplitude ** 2 * math.exp(-4.0 * nu * t)
        gxp = quad * (-2.0 * math.sin(2.0 * x))
        gyp = quad * (-2.0 * math.sin(2.0 * y))
        return (dudt[0] + conv[0] + gxp - nu * (-2.0 * u[0]),
                dudt[1] + conv[1] + gyp - nu * (-2.0 * u[1]),
                dudt[2] + conv[2] - nu * (-2.0 * u[2]))

    def energy(self, t: float) -> float:
        return 0.25 * self.amplitude ** 2 * math.exp(-4.0 * self.viscosity * t)

    def enstrophy(self, t: float) -> float:
        return 0.5 * self.amplitude ** 2 * math.exp(-4.0 * self.viscosity * t)

    def decay_rate(self) -> float:
        """``2 nu |k|^2`` with ``|k|^2 = 2`` -> ``4 nu``."""
        return 4.0 * self.viscosity


# ---------------------------------------------------------------------------
# 2. ABC Beltrami flow: exact, z-dependent
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class BeltramiABC:
    """Exact Navier-Stokes solution with ``curl u = k u`` (Beltrami ABC flow).

    ``u = a(sin kz + cos ky, sin kx + cos kz, sin ky + cos kx)`` gives
    ``curl u = k u``, hence ``(u.grad)u = grad(|u|^2/2)`` and the momentum
    equation reduces to pure diffusion: ``u(t) = u(0) e^{-k^2 nu t}`` with
    ``p = -|u|^2/2 + p0``.
    """

    amplitude: float = 1.0
    wavenumber: float = 1.0
    viscosity: float = 0.01
    pressure_datum: float = 0.0

    def decay(self, t: float) -> float:
        return math.exp(-self.wavenumber ** 2 * self.viscosity * t)

    def velocity(self, x: float, y: float, z: float,
                 t: float) -> tuple[float, float, float]:
        a, k = self.amplitude * self.decay(t), self.wavenumber
        return (a * (math.sin(k * z) + math.cos(k * y)),
                a * (math.sin(k * x) + math.cos(k * z)),
                a * (math.sin(k * y) + math.cos(k * x)))

    def velocity_gradient(self, x: float, y: float, z: float,
                          t: float) -> tuple[tuple[float, float, float], ...]:
        """``g[j][i] = d u_i / d x_j``."""
        a, k = self.amplitude * self.decay(t), self.wavenumber
        sk = lambda v: math.sin(k * v)  # noqa: E731
        ck = lambda v: math.cos(k * v)  # noqa: E731
        return (
            (0.0, a * k * ck(x), -a * k * sk(x)),
            (-a * k * sk(y), 0.0, a * k * ck(y)),
            (a * k * ck(z), -a * k * sk(z), 0.0),
        )

    def vorticity(self, x: float, y: float, z: float,
                  t: float) -> tuple[float, float, float]:
        u = self.velocity(x, y, z, t)
        k = self.wavenumber
        return (k * u[0], k * u[1], k * u[2])

    def pressure(self, x: float, y: float, z: float, t: float) -> float:
        u = self.velocity(x, y, z, t)
        return self.pressure_datum - 0.5 * (u[0] ** 2 + u[1] ** 2 + u[2] ** 2)

    def curl_eigen_residual(self, x: float, y: float, z: float, t: float) -> float:
        """``max_i |(curl u)_i - k u_i|`` -- the Beltrami property itself."""
        u = self.velocity(x, y, z, t)
        w = self.vorticity(x, y, z, t)
        return max(abs(w[i] - self.wavenumber * u[i]) for i in range(3))

    def ns_residual(self, x: float, y: float, z: float,
                    t: float) -> tuple[float, float, float]:
        """Momentum residual ``du/dt + (u.grad)u + grad p - nu lap u``.

        ``grad p`` is obtained from ``p = -|u|^2/2`` analytically, i.e. as
        ``-sum_j u_j d_i u_j``, *independently* of the convective sum below, so
        the residual equals ``(u.grad)u - grad(|u|^2/2)`` and only vanishes if
        both the Beltrami identity and the decay rate ``k^2 nu`` are right.
        """
        u = self.velocity(x, y, z, t)
        g = self.velocity_gradient(x, y, z, t)
        nu, k = self.viscosity, self.wavenumber
        dudt = tuple(-k * k * nu * u[i] for i in range(3))
        conv = tuple(sum(u[j] * g[j][i] for j in range(3)) for i in range(3))
        grad_ke = tuple(sum(u[j] * g[i][j] for j in range(3)) for i in range(3))
        lap = tuple(-k * k * u[i] for i in range(3))
        return tuple(dudt[i] + conv[i] - grad_ke[i] - nu * lap[i] for i in range(3))

    def energy(self, t: float) -> float:
        return 1.5 * self.amplitude ** 2 * math.exp(
            -2.0 * self.wavenumber ** 2 * self.viscosity * t)

    def enstrophy(self, t: float) -> float:
        return self.wavenumber ** 2 * self.energy(t)

    def decay_rate(self) -> float:
        """``2 nu k^2``."""
        return 2.0 * self.wavenumber ** 2 * self.viscosity# ---------------------------------------------------------------------------
# 3. classic 3-D Taylor-Green vortex: time-Taylor series reference
# ---------------------------------------------------------------------------
def _initial_tg3d_modes(amplitude: float) -> dict[Mode, Vec]:
    """Fourier modes of ``(A sin x cos y cos z, -A cos x sin y cos z, 0)``."""
    ux = _scalar_field([("sin", 0), ("cos", 1), ("cos", 2)])
    uy = _scalar_field([("cos", 0), ("sin", 1), ("cos", 2)])
    modes: dict[Mode, Vec] = {}
    for key, val in ux.items():
        _add(modes, key, (complex(val), 0j, 0j), amplitude)
    for key, val in uy.items():
        _add(modes, key, (0j, complex(-val), 0j), amplitude)
    return _project_div_free(modes)


@lru_cache(maxsize=8)
def _build_tg3d_series_cached(order: int, amplitude: float,
                              viscosity: float) -> tuple[dict[Mode, Vec], ...]:
    """Taylor coefficients ``u_0 .. u_order`` of the 3-D Taylor-Green vortex."""
    series: list[dict[Mode, Vec]] = [_initial_tg3d_modes(amplitude)]
    kmax = order + 1
    for n in range(order):
        conv: dict[Mode, Vec] = {}
        for p in range(n + 1):
            for key, vec in _convective(series[p], series[n - p], kmax).items():
                _add(conv, key, vec)
        conv = _project_div_free(conv)
        nxt: dict[Mode, Vec] = {}
        denom = n + 1
        for key, vec in conv.items():
            _add(nxt, key, vec, -1.0 / denom)
        for key, vec in series[n].items():
            kk = key[0] ** 2 + key[1] ** 2 + key[2] ** 2
            _add(nxt, key, vec, -viscosity * kk / denom)
        series.append({k: v for k, v in nxt.items()
                       if max(abs(c) for c in v) > 1e-300})
    return tuple(series)


def _build_tg3d_series(order: int, amplitude: float,
                       viscosity: float) -> list[dict[Mode, Vec]]:
    """Cached entry point (the build is the expensive part of this module)."""
    return list(_build_tg3d_series_cached(order, amplitude, viscosity))


@dataclass(frozen=True)
class TaylorGreen3DSeries:
    """Time-Taylor series of the standard three-dimensional Taylor-Green vortex.

    Initial condition on the periodic box ``[0, 2*pi]^3``::

        u0 = A sin x cos y cos z,  v0 = -A cos x sin y cos z,  w0 = 0

    Its convective term is not curl-free, so no closed form exists; the
    reference is ``u(x,t) = sum_{n<=order} u_n(x) t^n`` with the coefficients
    from the projected Navier-Stokes recursion (see the module docstring).
    """

    order: int = 12
    amplitude: float = 0.25
    viscosity: float = 0.01
    _series: list[dict[Mode, Vec]] = field(default_factory=list, init=False,
                                          repr=False, compare=False)

    def __post_init__(self) -> None:
        if self.order < 1:
            raise ValueError("order must be >= 1")
        if self.amplitude <= 0.0:
            raise ValueError("amplitude must be positive")
        if self.viscosity < 0.0:
            raise ValueError("viscosity must be non-negative")
        object.__setattr__(self, "_series",
                           _build_tg3d_series(self.order, self.amplitude,
                                              self.viscosity))

    @property
    def series(self) -> list[dict[Mode, Vec]]:
        return self._series

    def _modes_at(self, t: float, upto: int | None = None) -> dict[Mode, Vec]:
        last = len(self._series) if upto is None else min(upto + 1, len(self._series))
        out: dict[Mode, Vec] = {}
        for n in range(last):
            factor = t ** n
            if factor == 0.0:
                continue
            for key, vec in self._series[n].items():
                _add(out, key, vec, factor)
        return out

    def initial_velocity(self, x: float, y: float, z: float) -> tuple[float, float, float]:
        a = self.amplitude
        return (a * math.sin(x) * math.cos(y) * math.cos(z),
                -a * math.cos(x) * math.sin(y) * math.cos(z),
                0.0)

    def velocity_modes(self, t: float) -> dict[Mode, Vec]:
        return self._modes_at(t)

    def velocity(self, x: float, y: float, z: float,
                 t: float) -> tuple[float, float, float]:
        return _eval_modes(self._modes_at(t), x, y, z)

    def divergence(self, x: float, y: float, z: float, t: float) -> float:
        total = 0.0
        for (kx, ky, kz), vec in self._modes_at(t).items():
            phase = kx * x + ky * y + kz * z
            c, s = math.cos(phase), math.sin(phase)
            kd = 1j * (kx * vec[0] + ky * vec[1] + kz * vec[2])
            total += kd.real * c - kd.imag * s
        return total

    @staticmethod
    def energy_modes(modes: Mapping[Mode, Vec]) -> float:
        return 0.5 * sum(sum(abs(c) ** 2 for c in vec) for vec in modes.values())

    def energy(self, t: float) -> float:
        return self.energy_modes(self._modes_at(t))

    def enstrophy(self, t: float) -> float:
        total = 0.0
        for key, vec in self._modes_at(t).items():
            kk = key[0] ** 2 + key[1] ** 2 + key[2] ** 2
            total += kk * sum(abs(c) ** 2 for c in vec)
        return 0.5 * total

    def dissipation_rate(self, t: float) -> float:
        """``dE/dt`` from the series (no finite differences involved)."""
        modes_t = self._modes_at(t)
        total = 0.0
        for n in range(1, len(self._series)):
            factor = n * t ** (n - 1)
            if factor == 0.0:
                continue
            for key, vec in self._series[n].items():
                for i in range(3):
                    total += factor * (vec[i] * modes_t.get(key, (0j, 0j, 0j))[i].conjugate()).real
        return total

    def defect_spectrum(self) -> list[float]:
        """Max-norm NS residual per Taylor order of the truncated series.

        Entry ``n`` is the coefficient of ``t^n`` in
        ``du/dt + P[(u.grad)u] - nu lap u`` evaluated for the *truncated*
        series.  Entries ``0 .. order-1`` must vanish to round-off (the
        recursion is order-complete) and entry ``order`` is the leading
        truncation term.
        """
        order = self.order
        kmax = order + 1
        conv: dict[int, dict[Mode, Vec]] = {}
        for p in range(order + 1):
            for q in range(order + 1 - p):
                acc = conv.setdefault(p + q, {})
                for key, vec in _convective(self._series[p], self._series[q],
                                            kmax).items():
                    _add(acc, key, vec)
        out: list[float] = []
        for n in range(order + 1):
            resid: dict[Mode, Vec] = {}
            if n < order:
                for key, vec in self._series[n + 1].items():
                    _add(resid, key, vec, n + 1)
            for key, vec in _project_div_free(conv.get(n, {})).items():
                _add(resid, key, vec)
            for key, vec in self._series[n].items():
                kk = key[0] ** 2 + key[1] ** 2 + key[2] ** 2
                _add(resid, key, vec, self.viscosity * kk)
            out.append(max((max(abs(c) for c in v) for v in resid.values()),
                           default=0.0))
        return out

    def coefficient_ratio(self) -> float:
        """``||u_order|| / ||u_{order-1}||`` -- proxy for the series radius."""
        def norm(u: Mapping[Mode, Vec]) -> float:
            return math.sqrt(2.0 * self.energy_modes(u))
        hi, lo = norm(self._series[-1]), norm(self._series[-2])
        return hi / lo if lo > 0.0 else float("nan")

    def energy_terms(self) -> list[float]:
        """Taylor coefficients of ``E(t)`` (``E(t) = sum_n e_n t^n``)."""
        order = self.order
        terms: list[float] = []
        for n in range(order + 1):
            total = 0.0
            for p in range(n + 1):
                q = n - p
                up, uq = self._series[p], self._series[q]
                for key, vp in up.items():
                    vq = uq.get(key)
                    if vq is None:
                        continue
                    total += sum((vp[i] * vq[i].conjugate()).real for i in range(3))
            terms.append(0.5 * total)
        return terms


# ---------------------------------------------------------------------------
# metrics / consumption helpers
# ---------------------------------------------------------------------------
def fit_decay_rate(times: Sequence[float], values: Sequence[float],
                   t_min: float = 0.0) -> dict[str, float]:
    """Least-squares exponential decay rate ``sigma`` from ``v(t) ~ v0 e^{-sigma t}``."""
    if len(times) != len(values) or len(times) < 2:
        raise ValueError("need matching sequences of length >= 2")
    pts = [(t, v) for t, v in zip(times, values) if t >= t_min and v > 0.0]
    if len(pts) < 2:
        raise ValueError("not enough positive samples in the window")
    n = float(len(pts))
    sx = sum(p[0] for p in pts)
    sy = sum(math.log(p[1]) for p in pts)
    sxx = sum(p[0] ** 2 for p in pts)
    sxy = sum(p[0] * math.log(p[1]) for p in pts)
    denom = n * sxx - sx * sx
    if denom == 0.0:
        raise ValueError("degenerate time samples")
    slope = (n * sxy - sx * sy) / denom
    intercept = (sy - slope * sx) / n
    resid = max(abs(math.log(v) - (intercept + slope * t)) for t, v in pts)
    return {"rate": -slope, "amplitude": math.exp(intercept),
            "max_log_residual": resid, "samples": n}


def compare_metrics(reference: Mapping[str, float],
                    measured: Mapping[str, float],
                    tolerance: float = 0.08) -> dict[str, Any]:
    """Relative errors of solver metrics against the reference values.

    ``relative_tolerance`` follows the case-pool band (8 %).
    """
    rows: list[dict[str, Any]] = []
    for name, ref in reference.items():
        if name not in measured:
            continue
        got = float(measured[name])
        err = got - ref
        rel = abs(err) / abs(ref) if ref != 0.0 else float("nan")
        rows.append({"metric": name, "reference": ref, "measured": got,
                     "abs_error": err, "rel_error": rel,
                     "pass": bool(rel <= tolerance)})
    worst = max((r["rel_error"] for r in rows), default=float("nan"))
    return {"tolerance": tolerance, "rows": rows,
            "max_rel_error": worst,
            "all_pass": bool(rows) and all(r["pass"] for r in rows)}


def metrics_table(series: TaylorGreen3DSeries,
                  times: Sequence[float]) -> dict[str, Any]:
    """Reference ``E(t)``, ``Z(t)``, ``dE/dt`` and the fitted decay rate."""
    ts = list(times)
    energies = [series.energy(t) for t in ts]
    rows = [{"t": t, "energy": e, "enstrophy": series.enstrophy(t),
             "dissipation": series.dissipation_rate(t),
             "two_nu_Z": -2.0 * series.viscosity * series.enstrophy(t)}
            for t, e in zip(ts, energies)]
    return {"case": "CASE-TAYLOR-GREEN-3D",
            "order": series.order, "amplitude": series.amplitude,
            "viscosity": series.viscosity, "rows": rows,
            "energy_terms": series.energy_terms(),
            "defect_spectrum": series.defect_spectrum(),
            "coefficient_ratio": series.coefficient_ratio()}


def mode_table(series: TaylorGreen3DSeries, t: float,
               max_modes: int = 32) -> list[dict[str, Any]]:
    """Exact velocity Fourier modes at time ``t`` (compact reference form)."""
    items = sorted(series.velocity_modes(t).items(),
                   key=lambda kv: (-sum(abs(c) ** 2 for c in kv[1]),
                                   kv[0][0], kv[0][1], kv[0][2]))
    out: list[dict[str, Any]] = []
    for (kx, ky, kz), vec in items[:max_modes]:
        if max(abs(c) for c in vec) <= 1e-14:
            continue
        out.append({"mode": [kx, ky, kz],
                    "amplitude": math.sqrt(sum(abs(c) ** 2 for c in vec)),
                    "re": [vec[0].real, vec[1].real, vec[2].real],
                    "im": [vec[0].imag, vec[1].imag, vec[2].imag]})
    return out


def sample_grid(field: Any, level: int, t: float) -> dict[str, Any]:
    """Velocity sample on a ``level^3`` uniform grid of ``[0, 2*pi)^3``.

    ``field`` is any object exposing ``velocity(x, y, z, t)``.
    """
    if level < 2:
        raise ValueError("level must be >= 2")
    h = _TWO_PI / level
    coords = [i * h for i in range(level)]
    u = [[[0.0] * level for _ in range(level)] for _ in range(level)]
    v = [[[0.0] * level for _ in range(level)] for _ in range(level)]
    w = [[[0.0] * level for _ in range(level)] for _ in range(level)]
    for i, x in enumerate(coords):
        for j, y in enumerate(coords):
            for k, z in enumerate(coords):
                vel = field.velocity(x, y, z, t)
                u[i][j][k], v[i][j][k], w[i][j][k] = vel
    return {"level": level, "time": t, "spacing": h, "coordinates": coords,
            "u": u, "v": v, "w": w}


def case_spec(series: TaylorGreen3DSeries, times: Sequence[float]) -> dict[str, Any]:
    """Solver-consumable specification of the case (dimensionless numbers)."""
    e0, z0 = series.energy(0.0), series.enstrophy(0.0)
    return {
        "case": "CASE-TAYLOR-GREEN-3D",
        "domain": {"box": [0.0, _TWO_PI, 0.0, _TWO_PI, 0.0, _TWO_PI],
                   "boundary": "periodic"},
        "reference": {
            "type": "time-Taylor series (this module, residual-checked)",
            "order": series.order,
            "coefficient_ratio": series.coefficient_ratio(),
            "trusted_t": [0.0, float(min(times))] if times else [0.0, 0.0],
        },
        "initial_condition": {
            "u": "A sin(x) cos(y) cos(z)",
            "v": "-A cos(x) sin(y) cos(z)",
            "w": "0",
            "amplitude": series.amplitude,
        },
        "parameters": {
            "amplitude_A": series.amplitude,
            "kinematic_viscosity_nu": series.viscosity,
            "reynolds_like_Re": (series.amplitude / series.viscosity
                                 if series.viscosity > 0 else float("inf")),
            "energy_E0": e0,
            "enstrophy_Z0": z0,
            "initial_dissipation_dE_dt": -2.0 * series.viscosity * z0,
        },
        "metrics": CASE_TAYLOR_GREEN_3D_ACCEPTANCE["metrics"],
        "relative_tolerance": CASE_TAYLOR_GREEN_3D_ACCEPTANCE["relative_tolerance"],
        "sample_times": list(times),
    }


def _write_json(path: Path, payload: Any, compact: bool = False) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    text = (json.dumps(payload, separators=(",", ":"), sort_keys=True)
            if compact else json.dumps(payload, indent=2, sort_keys=True))
    path.write_text(text + "\n", encoding="utf-8")
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def write_case_bundle(outdir: str | Path, series: TaylorGreen3DSeries | None = None,
                      times: Sequence[float] = (1.0, 2.0, 3.0, 4.0),
                      grid_levels: Sequence[int] = (8, 16),
                      closed_forms: bool = True) -> dict[str, Any]:
    """Write the reference bundle a solver repo can consume; returns a manifest."""
    series = series or TaylorGreen3DSeries()
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    manifest: dict[str, Any] = {"outdir": str(out), "case": "CASE-TAYLOR-GREEN-3D",
                                "files": {}}
    manifest["files"]["case_spec.json"] = _write_json(
        out / "case_spec.json", case_spec(series, times))
    manifest["files"]["metrics_table.json"] = _write_json(
        out / "metrics_table.json", metrics_table(series, times))
    for t in times:
        manifest["files"][f"mode_table_t{t:g}.json"] = _write_json(
            out / f"mode_table_t{t:g}.json", mode_table(series, t), compact=True)
    for level in grid_levels:
        manifest["files"][f"grid_{level}_t1.json"] = _write_json(
            out / f"grid_{level}_t1.json", sample_grid(series, level, 1.0),
            compact=True)
    if closed_forms:
        manifest["files"]["closed_forms.json"] = _write_json(
            out / "closed_forms.json", _closed_form_payload(times))
    return manifest


def _closed_form_payload(times: Sequence[float]) -> dict[str, Any]:
    mode = TaylorGreenMode()
    abc = BeltramiABC()
    rows = [{
        "t": t,
        "invariant_mode_energy": mode.energy(t),
        "invariant_mode_enstrophy": mode.enstrophy(t),
        "invariant_mode_analytic_rate": mode.decay_rate(),
        "abc_energy": abc.energy(t),
        "abc_enstrophy": abc.enstrophy(t),
        "abc_analytic_rate": abc.decay_rate(),
    } for t in times]
    return {"TaylorGreenMode": {"amplitude": mode.amplitude,
                                "viscosity": mode.viscosity},
            "BeltramiABC": {"amplitude": abc.amplitude,
                            "wavenumber": abc.wavenumber,
                            "viscosity": abc.viscosity},
            "rows": rows}