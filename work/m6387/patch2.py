import pathlib
path = pathlib.Path("tensorcad/planing.py")
src = path.read_text()
start = src.index("def _wetted_topside_area(")
end = src.index("def analytic_planing_hydrostatics(")
new = '''#: 10-point Gauss-Legendre nodes and weights on ``[-1, 1]`` (exact to degree 19).
_GL10_X = (
    -0.9739065285171717, -0.8650633666889845, -0.6794095682990244,
    -0.4333953941292472, -0.1488743389816312, 0.1488743389816312,
    0.4333953941292472, 0.6794095682990244, 0.8650633666889845,
    0.9739065285171717,
)
_GL10_W = (
    0.0666713443086881, 0.1494513491505806, 0.2190863625159820,
    0.2692667193099963, 0.2955242247147529, 0.2955242247147529,
    0.2692667193099963, 0.2190863625159820, 0.1494513491505806,
    0.0666713443086881,
)
TOPSIDE_QUADRATURE_RTOL = 1e-12
TOPSIDE_QUADRATURE_MAX_DEPTH = 24


def _gauss10(f, lo: float, hi: float) -> float:
    centre = 0.5 * (lo + hi)
    half = 0.5 * (hi - lo)
    total = 0.0
    for node, weight in zip(_GL10_X, _GL10_W):
        total += weight * f(centre + half * node)
    return half * total


def _adaptive_gauss10(f, lo: float, hi: float, rtol: float, depth: int = 0) -> float:
    """Adaptive 10-point Gauss-Legendre on ``[lo, hi]``, relative tolerance ``rtol``.

    Bisection is driven by the difference between one panel and its two halves, which
    for a 10-point rule overestimates the error by orders of magnitude and so stops
    early on the smooth parts of the integrand while still resolving a boundary layer
    (``q = 0.05`` puts a ``t**39`` shape against ``t = 1``).
    """
    whole = _gauss10(f, lo, hi)
    middle = 0.5 * (lo + hi)
    split = _gauss10(f, lo, middle) + _gauss10(f, middle, hi)
    if depth >= TOPSIDE_QUADRATURE_MAX_DEPTH:
        return split
    if abs(split - whole) <= rtol * abs(split):
        return split
    return _adaptive_gauss10(f, lo, middle, rtol, depth + 1) + _adaptive_gauss10(
        f, middle, hi, rtol, depth + 1
    )


def _wetted_topside_area(p: PlaningHullParameters, x0: float,
                         rtol: float = TOPSIDE_QUADRATURE_RTOL) -> float:
    """``2 int_0^x0 (T - b tan(beta)) sqrt(1 + b'^2) dx``, adaptively quadratured.

    Evaluated in the substituted variable ``x = x0 t**e`` so that neither the integrand
    nor its factors may leave double range. Three numerical traps are avoided:

    * ``b'(x) = c q x**(q-1)`` diverges at a knife bow (``q < 1``). ``e = 2/q`` turns
      ``b`` into ``b(x0) t**2`` and the arc element into
      ``sqrt((dx/dt)**2 + (db/dt)**2)`` with ``db/dt = 2 b(x0) t``, i.e. a smooth
      integrand. For ``q >= 2`` that substitution would *introduce* a ``t**(2/q-1)``
      singularity at ``t = 0``, so ``e = 1`` there instead (``b'`` already vanishes at
      the stem for ``q > 1``).
    * ``sqrt(1 + b'**2) * dx`` is formed as ``hypot(dx/dt, db/dt)`` rather than as
      ``sqrt(1 + slope**2) * dx/dt``. At ``q = 0.05``, ``beta >= 30 deg`` and a light
      draft the slope reaches ``~1e21``; squaring it overflowed to ``inf`` and the
      whole wetted area came back ``inf``.
    * the integral is adaptive instead of a fixed 4000-interval Simpson rule: with a
      small ``q`` the arc element varies by ``e**(2/q)`` across ``[0, 1]``, so the
      uniform grid was 0.06 % short at ``q = 0.05`` and 0.85 % short at ``q >= 2``
      (see ``docs/planing_numerical_robustness.md``).
    """
    if x0 <= 0.0:
        return 0.0
    q = p.planform_exponent
    tau = p.tan_beta
    exponent = 2.0 / q if q < 2.0 else 1.0
    b_x0 = chine_half_beam(p, x0)
    power = exponent * q

    def integrand(t: float) -> float:
        height = p.draft - tau * b_x0 * t**power
        if height <= 0.0:
            return 0.0
        dx_dt = x0 * exponent * (1.0 if exponent == 1.0 else t ** (exponent - 1.0))
        db_dt = b_x0 * power * (1.0 if power == 1.0 else t ** (power - 1.0))
        return 2.0 * height * math.hypot(dx_dt, db_dt)

    if not rtol > 0.0:
        raise ValueError(f"rtol must be positive, got {rtol!r}")
    return _adaptive_gauss10(integrand, 0.0, 1.0, rtol)


'''
src = src[:start] + new + src[end:]
src = src.replace(
    "def analytic_wetted_area(p: PlaningHullParameters, samples: int = 4000) -> float:",
    "def analytic_wetted_area(p: PlaningHullParameters,\n                         rtol: float = TOPSIDE_QUADRATURE_RTOL) -> float:")
src = src.replace("    return bottom + _wetted_topside_area(p, x0, samples) + section_area(p, p.length)",
                  "    return bottom + _wetted_topside_area(p, x0, rtol) + section_area(p, p.length)")
src = src.replace(
    """    elementary closed form once ``q < 1`` (``b'`` diverges at a knife bow), so it is
    evaluated by composite Simpson in the substituted variable ``x = x0 t**e`` with
    ``e = 2/q`` (``q < 2``) or ``e = 1`` (``q >= 2``), which removes the endpoint
    singularity without creating one at ``t = 0``; see
    :func:`_wetted_topside_area`. The transom is closed form.""",
    """    elementary closed form once ``q < 1`` (``b'`` diverges at a knife bow), so it is
    evaluated by adaptive Gauss-Legendre quadrature in the substituted variable
    ``x = x0 t**e`` with ``e = 2/q`` (``q < 2``) or ``e = 1`` (``q >= 2``), which removes
    the endpoint singularity without creating one at ``t = 0``; see
    :func:`_wetted_topside_area`. The transom is closed form.""",
)
path.write_text(src)
print("ok")
