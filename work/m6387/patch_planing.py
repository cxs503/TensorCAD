import re, pathlib
path = pathlib.Path("tensorcad/planing.py")
src = path.read_text()

old_start = src.index("def _wetted_topside_area(")
old_end = src.index("def analytic_planing_hydrostatics(")
new_block = '''def _wetted_topside_area(p: PlaningHullParameters, x0: float, samples: int = 4000) -> float:
    """``2 int_0^x0 (T - b tan(beta)) sqrt(1 + b'^2) dx`` by Simpson in ``t``.

    Evaluated in the substituted variable ``x = x0 t**e`` so that neither the
    integrand nor its factors may leave double range. Two numerical traps are avoided:

    * ``b'(x) = c q x**(q-1)`` diverges at a knife bow (``q < 1``). ``e = 2/q`` turns
      ``b`` into ``b(x0) t**2`` and the arc element into
      ``sqrt((dx/dt)**2 + (db/dt)**2)`` with ``db/dt = 2 b(x0) t``, i.e. a smooth
      integrand. For ``q >= 2`` that substitution would *introduce* a ``t**(2/q-1)``
      singularity at ``t = 0``, so ``e = 1`` there instead (``b'`` already vanishes at
      the stem for ``q > 1``).
    * ``sqrt(1 + b'**2) * dx`` is formed as ``hypot(dx/dt, db/dt)`` rather than as
      ``sqrt(1 + slope**2) * dx/dt``. At ``q = 0.05``, ``beta >= 30 deg`` and a light
      draft the slope reaches ``~1e21``; squaring it overflowed to ``inf`` and the
      whole wetted area came back ``inf`` (see
      ``docs/planing_numerical_robustness.md``). ``hypot`` cannot overflow because it
      never forms the square.
    """
    if x0 <= 0.0:
        return 0.0
    if samples % 2:
        samples += 1
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

    h = 1.0 / samples
    total = integrand(0.0) + integrand(1.0)
    for i in range(1, samples):
        total += (4.0 if i % 2 else 2.0) * integrand(i * h)
    return total * h / 3.0


'''
src = src[:old_start] + new_block + src[old_end:]

# analytic_wetted_area docstring: substitution is no longer always 2/q
src = src.replace(
    """    elementary closed form once ``q < 1`` (``b'`` diverges at a knife bow), so it is
    evaluated by composite Simpson in the substituted variable ``x = x0 t**(2/q)``,
    which removes the endpoint singularity. The transom is closed form.""",
    """    elementary closed form once ``q < 1`` (``b'`` diverges at a knife bow), so it is
    evaluated by composite Simpson in the substituted variable ``x = x0 t**e`` with
    ``e = 2/q`` (``q < 2``) or ``e = 1`` (``q >= 2``), which removes the endpoint
    singularity without creating one at ``t = 0``; see
    :func:`_wetted_topside_area`. The transom is closed form.""",
)

# non-finite / zero analytic reference must be an ERROR, not a warning
old_check = '''        ref = reference[key]
        if name == "planing_lcb":'''
new_check = '''        ref = reference[key]
        if not math.isfinite(ref) or ref == 0.0:
            checks.append(
                Check(
                    name,
                    False,
                    f"closed form is not a usable reference ({ref!r}); the mesh value "
                    f"{measured:.6f} cannot be validated against it",
                    severity=ERROR,
                )
            )
            continue
        if name == "planing_lcb":'''
assert old_check in src
src = src.replace(old_check, new_check)

src = src.replace(
    """from .validation import (
    WARNING,""",
    """from .validation import (
    ERROR,
    WARNING,""",
)
path.write_text(src)
print("patched", len(src))
