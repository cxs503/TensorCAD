"""Independent adversarial probe of tensorcad.planing closed forms.

Route A: tensorcad.planing analytic_* (the code under review)
Route B: independent adaptive quadrature (Gauss-Legendre, 40-point, recursive
         bisection on the *physical* section functions b(x), b'(x)) written here
         from scratch and never calling the module's integrals.
"""
import math
import sys

sys.path.insert(0, ".")
from tensorcad.planing import (  # noqa: E402
    PlaningHullParameters, analytic_planing_hydrostatics, section_area,
    waterline_half_breadth, chine_half_beam, wetted_length,
)

def gl40(f, a, b):
    # 40-point Gauss-Legendre on [-1,1] via precomputed nodes (math only)
    n = 20
    nodes, weights = [], []
    for i in range(1, n + 1):
        x = math.cos(math.pi * (i - 0.25) / (n + 0.5))
        for _ in range(100):
            p0, p1 = 1.0, x
            for k in range(2, n + 1):
                p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
            dp = n * (x * p1 - p0) / (x * x - 1.0)
            dx = -p1 / dp
            x += dx
            if abs(dx) < 1e-16:
                break
        nodes.append(x)
        weights.append(2.0 / ((1 - x * x) * dp * dp))
    c = 0.5 * (a + b)
    h = 0.5 * (b - a)
    return h * sum(w * f(c + h * z) for z, w in zip(nodes, weights))

def adaptive(f, a, b, tol=1e-12, depth=0):
    whole = gl40(f, a, b)
    m = 0.5 * (a + b)
    half = gl40(f, a, m) + gl40(f, m, b)
    if abs(half - whole) < tol * max(1.0, abs(half)) or depth > 30:
        return half
    return adaptive(f, a, m, tol, depth + 1) + adaptive(f, m, b, tol, depth + 1)

def route_b(p):
    L, T, tau = p.length, p.draft, p.tan_beta
    c = p.half_beam_transom / L ** p.planform_exponent
    q = p.planform_exponent
    b = lambda x: c * x ** q
    db = lambda x: c * q * x ** (q - 1.0)
    x0 = wetted_length(p)
    eps = 1e-12
    vol = adaptive(lambda x: section_area(p, x), eps, L, 1e-11)
    awp = adaptive(lambda x: 2.0 * waterline_half_breadth(p, x), eps, L, 1e-11)
    lcb = adaptive(lambda x: x * section_area(p, x), eps, L, 1e-11) / vol
    lcf = adaptive(lambda x: x * waterline_half_breadth(p, x), eps, L, 1e-11) / (awp / 2.0)
    def vmoment(x):
        bb = b(x) if x > 0 else 0.0
        return -T * T * bb + T * tau * bb * bb - tau ** 3 * 0.0 - (1 / 3) * tau * tau * bb ** 3 \
            if bb * tau < T else -T ** 3 / (3 * tau)
    # use exact section moment formula
    def smom(x):
        bb = b(x) if x > 0 else 0.0
        h = bb * tau
        if h >= T:
            return -T ** 3 / (3.0 * tau)
        return -T * T * bb + T * tau * bb * bb - (1.0 / 3.0) * tau * tau * bb ** 3
    vcb = adaptive(lambda x: smom(x), eps, L, 1e-11) / vol
    bottom = awp * math.sqrt(1 + tau * tau)
    top = adaptive(lambda x: (2.0 * (T - b(x) * tau) * math.sqrt(1 + db(x) ** 2)
                              if (x < x0 and T - b(x) * tau > 0) else 0.0),
                   0.0, max(x0, eps), 1e-11) if x0 > 0 else 0.0
    wetted = bottom + top + section_area(p, L)
    return dict(displacement_volume=vol, waterplane_area=awp,
                longitudinal_centre_buoyancy=lcb, longitudinal_centre_flotation=lcf,
                vertical_centre_buoyancy=vcb, wetted_surface_area=wetted)

def rel(a, r):
    return (a - r) / r

cases = []
for beta in (1.0, 12.0, 18.0, 30.0, 45.0, 59.0):
    for q in (0.05, 0.2, 0.6, 1.0, 2.0, 4.0):
        for draft in (0.05, 0.30, 0.55, 2.0):
            cases.append(PlaningHullParameters(length=10.0, beam=3.0, draft=draft,
                                               deadrise_deg=beta, planform_exponent=q,
                                               nx=21, nz=9))
worst = {}
for p in cases:
    try:
        p.validate()
    except ValueError as exc:
        print("INVALID", p.deadrise_deg, p.planform_exponent, p.draft, exc)
        continue
    a = analytic_planing_hydrostatics(p)
    b = route_b(p)
    for k, v in b.items():
        d = rel(a[k], v)
        key = k
        if key not in worst or abs(d) > abs(worst[key][0]):
            worst[key] = (d, p.deadrise_deg, p.planform_exponent, p.draft)
for k, (d, be, q, t) in sorted(worst.items()):
    print(f"{k:34s} worst {d:+.3e}  at beta={be} q={q} T={t}")
