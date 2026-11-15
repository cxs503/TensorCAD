"""Fuse value+partials in _gauss_newton (idempotent)."""
from pathlib import Path

path = Path("/root/tp-agent/work/m22464/tensorcad/geomquality/trim_mesh.py")
text = path.read_text()

old = """    for _ in range(iterations):
        base = surface.evaluate(u, v)
        residual = _sub(base, point)
        error = _norm(residual)
        if error <= tolerance:
            break
        exact_partials = getattr(surface, "partials", None) if exact else None
        if exact_partials is not None:
            _, du, dv = exact_partials(u, v)
        else:
            du = _scale(_sub(surface.evaluate(min(u + step_u, u1), v), base), 1.0 / step_u)
            dv = _scale(_sub(surface.evaluate(u, min(v + step_v, v1)), base), 1.0 / step_v)
"""
new = """    exact_partials = getattr(surface, "partials", None) if exact else None
    for _ in range(iterations):
        if exact_partials is not None:
            # One derivative-net evaluation yields the point *and* both tangents; asking
            # for ``evaluate`` and then ``partials`` doubled the cost of every iteration.
            base, du, dv = exact_partials(u, v)
        else:
            base = surface.evaluate(u, v)
        residual = _sub(base, point)
        error = _norm(residual)
        if error <= tolerance:
            break
        if exact_partials is None:
            du = _scale(_sub(surface.evaluate(min(u + step_u, u1), v), base), 1.0 / step_u)
            dv = _scale(_sub(surface.evaluate(u, min(v + step_v, v1)), base), 1.0 / step_v)
"""
assert text.count(old) == 1
path.write_text(text.replace(old, new))
print("patched", path)