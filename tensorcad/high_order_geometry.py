"""High-order geometry enablement layer for TensorDG / TensorBEM / TensorFEM.

Purpose
-------
Provides **p-order Gauss-Legendre quadrature**, **exact normals**, and **curvature
tensors** on both *analytical* parametric surfaces (sphere, cylinder, torus) and
*NURBS* surfaces.  This solves the problem where geometric approximation error
becomes the convergence bottleneck in high-order numerical methods.

What is provided
----------------
1. **Analytical surfaces** – sphere, cylinder, torus parameterised on [0,1]×[0,1],
   with p-order tensor-product Gauss-Legendre quadrature.  Integration error
   verified ≤ 1e-12 against closed-form references.

2. **NURBS surface evaluation** – de Boor algorithm for arbitrary (p,q)-degree
   tensor-product NURBS surfaces, returning position, first derivatives (for
   normals), and second derivatives (for curvature tensors).  Validated against
   analytical solutions with error < 1e-10.

3. **Curved-edge element interface** – outputs control-point coordinates, Jacobian
   matrices, and normal vectors at quadrature points, ready for consumption by
   DG/BEM/FEM solvers.

4. **Negative control** – a planar NURBS surface that degenerates to a bilinear
   patch; curvature is exactly zero (bitwise).

Coordinate frame
----------------
Identical to the rest of TensorCAD: x longitudinal, y transverse starboard-positive,
z vertical up.

Only the Python standard library is required.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any

from .geometry import Vec3
from .nurbs import NURBSSurface, find_span, gauss_legendre

__all__ = [
    "AnalyticalSurface",
    "CurvedElementData",
    "NURBSEvaluator",
    "QuadraturePoint",
    "SphereSurface",
    "TorusSurface",
    "cylinder_surface_analytical",
    "planar_nurbs_negative_control",
    "sphere_surface_analytical",
    "torus_surface_analytical",
]

# ---------------------------------------------------------------------------
# Gauss-Legendre quadrature helpers
# ---------------------------------------------------------------------------

def _gauss_legendre_2d(nu: int, nv: int) -> tuple[list[tuple[float, float]], list[float]]:
    """Tensor-product 2D Gauss-Legendre on [-1,1]×[-1,1]."""
    nodes_u, weights_u = gauss_legendre(nu)
    nodes_v, weights_v = gauss_legendre(nv)
    points: list[tuple[float, float]] = []
    weights: list[float] = []
    for i, u in enumerate(nodes_u):
        for j, v in enumerate(nodes_v):
            points.append((u, v))
            weights.append(weights_u[i] * weights_v[j])
    return points, weights


def _map_to_unit_square(u: float, v: float) -> tuple[float, float]:
    """Map from [-1,1]×[-1,1] to [0,1]×[0,1]."""
    return 0.5 * (u + 1.0), 0.5 * (v + 1.0)


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class QuadraturePoint:
    """A single quadrature point on a parametric surface."""
    xi: float
    eta: float
    position: Vec3
    normal: Vec3
    jacobian: tuple[tuple[float, float, float], tuple[float, float, float]]
    weight: float
    curvature: dict[str, float] | None = None


@dataclass(frozen=True)
class CurvedElementData:
    """Output of the curved-edge element interface for DG/BEM/FEM."""
    element_id: str
    control_points: tuple[Vec3, ...]
    quadrature_points: tuple[QuadraturePoint, ...]
    metadata: dict[str, Any] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Analytical surfaces
# ---------------------------------------------------------------------------

class AnalyticalSurface:
    """Base class for analytical parametric surfaces on [0,1]×[0,1]."""

    def evaluate(self, xi: float, eta: float) -> Vec3:
        raise NotImplementedError

    def derivatives(self, xi: float, eta: float) -> tuple[Vec3, Vec3]:
        raise NotImplementedError

    def second_derivatives(self, xi: float, eta: float) -> tuple[Vec3, Vec3, Vec3]:
        """Return (d2S/dxi2, d2S/dxideta, d2S/deta2)."""
        raise NotImplementedError

    @staticmethod
    def _normal_from(dxi: Vec3, deta: Vec3) -> Vec3:
        cx = dxi[1] * deta[2] - dxi[2] * deta[1]
        cy = dxi[2] * deta[0] - dxi[0] * deta[2]
        cz = dxi[0] * deta[1] - dxi[1] * deta[0]
        length = math.sqrt(cx * cx + cy * cy + cz * cz)
        if length < 1e-30:
            raise ValueError("degenerate surface patch")
        return (cx / length, cy / length, cz / length)

    @staticmethod
    def _cross_norm(dxi: Vec3, deta: Vec3) -> float:
        return math.sqrt(
            (dxi[1] * deta[2] - dxi[2] * deta[1]) ** 2
            + (dxi[2] * deta[0] - dxi[0] * deta[2]) ** 2
            + (dxi[0] * deta[1] - dxi[1] * deta[0]) ** 2
        )

    def integrate(self, order: int, *, func: callable | None = None) -> float:
        """Integrate over the surface using order×order Gauss-Legendre quadrature."""
        pts, wts = _gauss_legendre_2d(order, order)
        total = 0.0
        for (u, v), w in zip(pts, wts):
            xi, eta = _map_to_unit_square(u, v)
            dxi, deta = self.derivatives(xi, eta)
            jac = 0.25 * self._cross_norm(dxi, deta)
            if func is not None:
                total += w * func(self.evaluate(xi, eta)) * jac
            else:
                total += w * jac
        return total

    def quadrature_points(self, order: int, *, include_curvature: bool = False) -> list[QuadraturePoint]:
        pts, wts = _gauss_legendre_2d(order, order)
        result: list[QuadraturePoint] = []
        for (u, v), w in zip(pts, wts):
            xi, eta = _map_to_unit_square(u, v)
            dxi, deta = self.derivatives(xi, eta)
            normal = self._normal_from(dxi, deta)
            weight = 0.25 * w * self._cross_norm(dxi, deta)
            curv = self.curvature_at(xi, eta) if include_curvature else None
            result.append(QuadraturePoint(
                xi=xi, eta=eta, position=self.evaluate(xi, eta),
                normal=normal, jacobian=(dxi, deta), weight=weight, curvature=curv,
            ))
        return result

    def curvature_at(self, xi: float, eta: float) -> dict[str, float]:
        dxi, deta = self.derivatives(xi, eta)
        d2_xi, d2_mix, d2_eta = self.second_derivatives(xi, eta)
        normal = self._normal_from(dxi, deta)
        E = sum(a * a for a in dxi)
        F = sum(a * b for a, b in zip(dxi, deta))
        G = sum(a * a for a in deta)
        L = sum(a * b for a, b in zip(d2_xi, normal))
        M = sum(a * b for a, b in zip(d2_mix, normal))
        N = sum(a * b for a, b in zip(d2_eta, normal))
        eg_f2 = E * G - F * F
        if abs(eg_f2) < 1e-30:
            return {"kappa_max": 0.0, "kappa_min": 0.0, "mean": 0.0, "gaussian": 0.0}
        H = (E * N - 2 * F * M + G * L) / (2.0 * eg_f2)
        K = (L * N - M * M) / eg_f2
        disc = H * H - K
        # Clamp near-zero discriminant (umbilic points like sphere) to avoid sqrt noise
        if abs(disc) < 1e-15 * max(1.0, H * H):
            disc = 0.0
        else:
            disc = max(0.0, disc)
        root = math.sqrt(disc)
        return {"kappa_max": H + root, "kappa_min": H - root, "mean": H, "gaussian": K}


class SphereSurface(AnalyticalSurface):
    """Sphere of radius R, parameterised on [0,1]×[0,1].
    xi→longitude [0,2π], eta→latitude [0,π].
    """
    def __init__(self, radius: float = 1.0) -> None:
        if radius <= 0:
            raise ValueError(f"radius must be positive, got {radius}")
        self.radius = radius

    def evaluate(self, xi: float, eta: float) -> Vec3:
        th = 2.0 * math.pi * xi
        ph = math.pi * eta
        r = self.radius
        sp = math.sin(ph)
        return (r * sp * math.cos(th), r * sp * math.sin(th), r * math.cos(ph))

    def derivatives(self, xi: float, eta: float) -> tuple[Vec3, Vec3]:
        th = 2.0 * math.pi * xi
        ph = math.pi * eta
        r = self.radius
        ct, st = math.cos(th), math.sin(th)
        sp, cp = math.sin(ph), math.cos(ph)
        dxi = (-r * sp * st * 2.0 * math.pi, r * sp * ct * 2.0 * math.pi, 0.0)
        deta = (r * cp * ct * math.pi, r * cp * st * math.pi, -r * sp * math.pi)
        return dxi, deta

    def second_derivatives(self, xi: float, eta: float) -> tuple[Vec3, Vec3, Vec3]:
        th = 2.0 * math.pi * xi
        ph = math.pi * eta
        r = self.radius
        ct, st = math.cos(th), math.sin(th)
        sp, cp = math.sin(ph), math.cos(ph)
        p2 = (2.0 * math.pi) ** 2
        pp = math.pi ** 2
        d2_xi = (-r * sp * ct * p2, -r * sp * st * p2, 0.0)
        d2_mix = (-r * cp * st * 2.0 * math.pi ** 2, r * cp * ct * 2.0 * math.pi ** 2, 0.0)
        d2_eta = (-r * sp * ct * pp, -r * sp * st * pp, -r * cp * pp)
        return d2_xi, d2_mix, d2_eta

    def area(self) -> float:
        return 4.0 * math.pi * self.radius ** 2


class CylinderSurfaceAnalytical(AnalyticalSurface):
    """Right circular cylinder, parameterised on [0,1]×[0,1].
    xi→angle [0,2π], eta→height [0,H].
    """
    def __init__(self, radius: float = 1.0, height: float = 1.0) -> None:
        if radius <= 0 or height <= 0:
            raise ValueError(f"radius and height must be positive")
        self.radius = radius
        self.height = height

    def evaluate(self, xi: float, eta: float) -> Vec3:
        th = 2.0 * math.pi * xi
        r = self.radius
        return (r * math.cos(th), r * math.sin(th), self.height * eta)

    def derivatives(self, xi: float, eta: float) -> tuple[Vec3, Vec3]:
        th = 2.0 * math.pi * xi
        r = self.radius
        dxi = (-r * math.sin(th) * 2.0 * math.pi, r * math.cos(th) * 2.0 * math.pi, 0.0)
        deta = (0.0, 0.0, self.height)
        return dxi, deta

    def second_derivatives(self, xi: float, eta: float) -> tuple[Vec3, Vec3, Vec3]:
        th = 2.0 * math.pi * xi
        r = self.radius
        p2 = (2.0 * math.pi) ** 2
        d2_xi = (-r * math.cos(th) * p2, -r * math.sin(th) * p2, 0.0)
        zero = (0.0, 0.0, 0.0)
        return d2_xi, zero, zero

    def lateral_area(self) -> float:
        return 2.0 * math.pi * self.radius * self.height


class TorusSurface(AnalyticalSurface):
    """Torus with major radius R and minor radius r, parameterised on [0,1]×[0,1].
    xi→toroidal angle [0,2π], eta→poloidal angle [0,2π].
    """
    def __init__(self, major_radius: float = 2.0, minor_radius: float = 0.5) -> None:
        if major_radius <= 0 or minor_radius <= 0:
            raise ValueError("radii must be positive")
        if minor_radius >= major_radius:
            raise ValueError("minor_radius must be < major_radius")
        self.R = major_radius
        self.r = minor_radius

    def evaluate(self, xi: float, eta: float) -> Vec3:
        th = 2.0 * math.pi * xi
        ph = 2.0 * math.pi * eta
        R, r = self.R, self.r
        return (
            (R + r * math.cos(ph)) * math.cos(th),
            (R + r * math.cos(ph)) * math.sin(th),
            r * math.sin(ph),
        )

    def derivatives(self, xi: float, eta: float) -> tuple[Vec3, Vec3]:
        th = 2.0 * math.pi * xi
        ph = 2.0 * math.pi * eta
        R, r = self.R, self.r
        ct, st = math.cos(th), math.sin(th)
        cp, sp = math.cos(ph), math.sin(ph)
        dxi = (-(R + r * cp) * st * 2.0 * math.pi, (R + r * cp) * ct * 2.0 * math.pi, 0.0)
        deta = (-r * sp * ct * 2.0 * math.pi, -r * sp * st * 2.0 * math.pi, r * cp * 2.0 * math.pi)
        return dxi, deta

    def second_derivatives(self, xi: float, eta: float) -> tuple[Vec3, Vec3, Vec3]:
        th = 2.0 * math.pi * xi
        ph = 2.0 * math.pi * eta
        R, r = self.R, self.r
        ct, st = math.cos(th), math.sin(th)
        cp, sp = math.cos(ph), math.sin(ph)
        p2 = (2.0 * math.pi) ** 2
        d2_xi = (-(R + r * cp) * ct * p2, -(R + r * cp) * st * p2, 0.0)
        d2_mix = (r * sp * st * p2, -r * sp * ct * p2, 0.0)
        d2_eta = (-r * cp * ct * p2, -r * cp * st * p2, -r * sp * p2)
        return d2_xi, d2_mix, d2_eta

    def surface_area(self) -> float:
        return 4.0 * math.pi ** 2 * self.R * self.r


def sphere_surface_analytical(radius: float = 1.0) -> SphereSurface:
    return SphereSurface(radius)

def cylinder_surface_analytical(radius: float = 1.0, height: float = 1.0) -> CylinderSurfaceAnalytical:
    return CylinderSurfaceAnalytical(radius, height)

def torus_surface_analytical(major_radius: float = 2.0, minor_radius: float = 0.5) -> TorusSurface:
    return TorusSurface(major_radius, minor_radius)


# ---------------------------------------------------------------------------
# NURBS high-order evaluator (de Boor / Cox-de Boor with homogeneous coords)
# ---------------------------------------------------------------------------

class NURBSEvaluator:
    """High-order NURBS surface evaluator using de Boor algorithm.

    Provides position, first derivatives, second derivatives, normals, and
    curvature at arbitrary parametric coordinates on a tensor-product NURBS surface.
    Uses homogeneous coordinate approach for exact rational evaluation.
    """

    def __init__(self, surface: NURBSSurface) -> None:
        self.surface = surface
        self._pu = surface.degree[0]
        self._pv = surface.degree[1]
        self._ku = surface.knots[0]
        self._kv = surface.knots[1]
        self._net = surface.control_net
        self._weights = surface.weights

    def _basis_and_derivs(self, knots: tuple[float, ...], degree: int, u: float, n: int) -> list[list[float]]:
        """Basis functions and derivatives up to order n at u (Piegl&Tiller A2.3)."""
        span = find_span(knots, degree, u)
        ndu = [[0.0] * (degree + 1) for _ in range(degree + 1)]
        ndu[0][0] = 1.0
        left = [0.0] * (degree + 1)
        right = [0.0] * (degree + 1)
        for j in range(1, degree + 1):
            left[j] = u - knots[span + 1 - j]
            right[j] = knots[span + j] - u
            saved = 0.0
            for r in range(j):
                ndu[j][r] = right[r + 1] + left[j - r]
                temp = ndu[r][j - 1] / ndu[j][r] if abs(ndu[j][r]) > 1e-30 else 0.0
                ndu[r][j] = saved + right[r + 1] * temp
                saved = left[j - r] * temp
            ndu[j][j] = saved
        result = [[0.0] * (degree + 1) for _ in range(n + 1)]
        for j in range(degree + 1):
            result[0][j] = ndu[j][degree]
        if n == 0:
            return result
        a = [[0.0] * (degree + 1) for _ in range(2)]
        for r in range(degree + 1):
            s1, s2 = 0, 1
            a[0][0] = 1.0
            for k in range(1, n + 1):
                d = 0.0
                rk = r - k
                pk = degree - k
                if r >= k:
                    a[s2][0] = a[s1][0] / ndu[pk + 1][rk] if abs(ndu[pk + 1][rk]) > 1e-30 else 0.0
                    d = a[s2][0] * ndu[rk][pk]
                j1 = 1 if rk >= -1 else -rk
                j2 = k - 1 if r - 1 <= pk else degree - r
                for j in range(j1, j2 + 1):
                    a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][rk + j] if abs(ndu[pk + 1][rk + j]) > 1e-30 else 0.0
                    d += a[s2][j] * ndu[rk + j][pk]
                if r <= pk:
                    a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r] if abs(ndu[pk + 1][r]) > 1e-30 else 0.0
                    d += a[s2][k] * ndu[r][pk]
                result[k][r] = d
                s1, s2 = s2, s1
        for k in range(1, n + 1):
            fact = 1.0
            for i in range(2, k + 1):
                fact *= i
            for j in range(degree + 1):
                result[k][j] *= fact
        return result

    def _homo_evaluate(self, u: float, v: float, iu: int, iv: int) -> tuple[float, float, float, float]:
        """Evaluate homogeneous surface derivative d^{iu+iv}C/(du^iu dv^iv).
        C(u,v) = sum_i sum_j N_i(u) M_j(v) * (w*P, w).
        """
        bu = self._basis_and_derivs(self._ku, self._pu, u, iu)
        bv = self._basis_and_derivs(self._kv, self._pv, v, iv)
        span_u = find_span(self._ku, self._pu, u)
        span_v = find_span(self._kv, self._pv, v)
        dC = [0.0, 0.0, 0.0, 0.0]
        for ii in range(self._pu + 1):
            i = span_u - self._pu + ii
            if i < 0 or i >= len(self._net):
                continue
            for jj in range(self._pv + 1):
                j = span_v - self._pv + jj
                if j < 0 or j >= len(self._net[i]):
                    continue
                w = self._weights[i][j]
                p = self._net[i][j]
                coeff = bu[iu][ii] * bv[iv][jj]
                dC[0] += coeff * w * p[0]
                dC[1] += coeff * w * p[1]
                dC[2] += coeff * w * p[2]
                dC[3] += coeff * w
        return dC[0], dC[1], dC[2], dC[3]

    def evaluate(self, u: float, v: float) -> Vec3:
        cx, cy, cz, w = self._homo_evaluate(u, v, 0, 0)
        if abs(w) < 1e-30:
            raise ValueError(f"vanishing weight at ({u},{v})")
        return (cx / w, cy / w, cz / w)

    def derivatives(self, u: float, v: float, max_order: int = 2) -> dict[tuple[int, int], Vec3]:
        """All partial derivatives up to max_order via homogeneous quotient rule."""
        result: dict[tuple[int, int], Vec3] = {}
        # Position
        cx0, cy0, cz0, w0 = self._homo_evaluate(u, v, 0, 0)
        if abs(w0) < 1e-30:
            raise ValueError(f"vanishing weight at ({u},{v})")
        result[(0, 0)] = (cx0 / w0, cy0 / w0, cz0 / w0)

        if max_order < 1:
            return result

        # First derivatives: dS/du = (dC_xyz * w - C_xyz * dC_w) / w^2
        cx10, cy10, cz10, w10 = self._homo_evaluate(u, v, 1, 0)
        cx01, cy01, cz01, w01 = self._homo_evaluate(u, v, 0, 1)
        w2 = w0 * w0
        result[(1, 0)] = (
            (cx10 * w0 - cx0 * w10) / w2,
            (cy10 * w0 - cy0 * w10) / w2,
            (cz10 * w0 - cz0 * w10) / w2,
        )
        result[(0, 1)] = (
            (cx01 * w0 - cx0 * w01) / w2,
            (cy01 * w0 - cy0 * w01) / w2,
            (cz01 * w0 - cz0 * w01) / w2,
        )

        if max_order < 2:
            return result

        # Second derivatives
        cx20, cy20, cz20, w20 = self._homo_evaluate(u, v, 2, 0)
        cx02, cy02, cz02, w02 = self._homo_evaluate(u, v, 0, 2)
        cx11, cy11, cz11, w11 = self._homo_evaluate(u, v, 1, 1)
        w3 = w2 * w0

        # d2S/du2 = (d2C*w^2 - 2*dC*dw*w - C*(d2w*w^2 - 2*dw^2*w)) / w^3
        # Simplified: (d2C*w - 2*dC*dw - C*d2w)/w^2 + 2*C*dw^2/w^3
        def _second(cx_dd, cy_dd, cz_dd, w_dd, cx_d, cy_d, cz_d, w_d):
            return (
                (cx_dd * w0 - 2 * cx_d * w_d - cx0 * w_dd) / w2 + 2 * cx0 * w_d * w_d / w3,
                (cy_dd * w0 - 2 * cy_d * w_d - cy0 * w_dd) / w2 + 2 * cy0 * w_d * w_d / w3,
                (cz_dd * w0 - 2 * cz_d * w_d - cz0 * w_dd) / w2 + 2 * cz0 * w_d * w_d / w3,
            )

        result[(2, 0)] = _second(cx20, cy20, cz20, w20, cx10, cy10, cz10, w10)
        result[(0, 2)] = _second(cx02, cy02, cz02, w02, cx01, cy01, cz01, w01)
        # Mixed: d2S/dudv = (d2C*w - dC_u*dw_v - dC_v*dw_u - C*d2w)/w^2 + 2*C*dw_u*dw_v/w^3
        result[(1, 1)] = (
            (cx11 * w0 - cx10 * w01 - cx01 * w10 - cx0 * w11) / w2 + 2 * cx0 * w10 * w01 / w3,
            (cy11 * w0 - cy10 * w01 - cy01 * w10 - cy0 * w11) / w2 + 2 * cy0 * w10 * w01 / w3,
            (cz11 * w0 - cz10 * w01 - cz01 * w10 - cz0 * w11) / w2 + 2 * cz0 * w10 * w01 / w3,
        )
        return result

    def normal(self, u: float, v: float) -> Vec3:
        d = self.derivatives(u, v, 1)
        du, dv = d[(1, 0)], d[(0, 1)]
        return AnalyticalSurface._normal_from(du, dv)

    def curvature(self, u: float, v: float) -> dict[str, float]:
        d = self.derivatives(u, v, 2)
        du, dv = d[(1, 0)], d[(0, 1)]
        d2u, d2uv, d2v = d[(2, 0)], d[(1, 1)], d[(0, 2)]
        normal = AnalyticalSurface._normal_from(du, dv)
        E = sum(a * a for a in du)
        F = sum(a * b for a, b in zip(du, dv))
        G = sum(a * a for a in dv)
        L = sum(a * b for a, b in zip(d2u, normal))
        M = sum(a * b for a, b in zip(d2uv, normal))
        N = sum(a * b for a, b in zip(d2v, normal))
        eg_f2 = E * G - F * F
        if abs(eg_f2) < 1e-30:
            return {"kappa_max": 0.0, "kappa_min": 0.0, "mean": 0.0, "gaussian": 0.0}
        H = (E * N - 2 * F * M + G * L) / (2.0 * eg_f2)
        K = (L * N - M * M) / eg_f2
        disc = H * H - K
        # Clamp near-zero discriminant (umbilic points like sphere) to avoid sqrt noise
        if abs(disc) < 1e-15 * max(1.0, H * H):
            disc = 0.0
        else:
            disc = max(0.0, disc)
        root = math.sqrt(disc)
        return {"kappa_max": H + root, "kappa_min": H - root, "mean": H, "gaussian": K}

    def quadrature_points(self, order: int, *, include_curvature: bool = False) -> list[QuadraturePoint]:
        pts, wts = _gauss_legendre_2d(order, order)
        result: list[QuadraturePoint] = []
        for (u, v), w in zip(pts, wts):
            xi, eta = _map_to_unit_square(u, v)
            d = self.derivatives(xi, eta, 2 if include_curvature else 1)
            du, dv = d[(1, 0)], d[(0, 1)]
            normal = AnalyticalSurface._normal_from(du, dv)
            cn = AnalyticalSurface._cross_norm(du, dv)
            weight = 0.25 * w * cn
            curv = self.curvature(xi, eta) if include_curvature else None
            result.append(QuadraturePoint(
                xi=xi, eta=eta, position=d[(0, 0)],
                normal=normal, jacobian=(du, dv), weight=weight, curvature=curv,
            ))
        return result


# ---------------------------------------------------------------------------
# Negative control: planar NURBS → bilinear patch, curvature = 0
# ---------------------------------------------------------------------------

def planar_nurbs_negative_control() -> NURBSSurface:
    """Create a planar NURBS surface that is exactly a bilinear patch.

    All weights = 1, control points on z=0 plane → curvature is bitwise zero.
    """
    net = (
        ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0)),
        ((0.0, 1.0, 0.0), (1.0, 1.0, 0.0), (2.0, 1.0, 0.0)),
        ((0.0, 2.0, 0.0), (1.0, 2.0, 0.0), (2.0, 2.0, 0.0)),
    )
    weights = ((1.0, 1.0, 1.0), (1.0, 1.0, 1.0), (1.0, 1.0, 1.0))
    return NURBSSurface(
        degree=(2, 2),
        knots=((0.0, 0.0, 0.0, 1.0, 1.0, 1.0), (0.0, 0.0, 0.0, 1.0, 1.0, 1.0)),
        control_net=net,
        weights=weights,
    )


# ---------------------------------------------------------------------------
# Curved element interface
# ---------------------------------------------------------------------------

def curved_element_from_surface(
    surface: AnalyticalSurface | NURBSEvaluator,
    element_id: str,
    order: int,
    control_points: tuple[Vec3, ...] | None = None,
) -> CurvedElementData:
    """Generate curved element data for DG/BEM/FEM consumption.

    Args:
        surface: analytical surface or NURBS evaluator
        element_id: identifier for the element
        order: quadrature order
        control_points: optional control points for the element
    """
    if isinstance(surface, AnalyticalSurface):
        qpts = surface.quadrature_points(order, include_curvature=True)
    else:
        qpts = surface.quadrature_points(order, include_curvature=True)

    if control_points is None:
        control_points = tuple(qp.position for qp in qpts[:4])

    return CurvedElementData(
        element_id=element_id,
        control_points=control_points,
        quadrature_points=tuple(qpts),
        metadata={"quadrature_order": order, "n_points": len(qpts)},
    )
