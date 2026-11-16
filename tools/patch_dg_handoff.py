"""One-shot, reviewable patch: ellipse contracts + explicit frame in the DG handoff.

Run once from the repository root; the patch is idempotent (it refuses to apply
twice instead of duplicating code).
"""
from __future__ import annotations

from pathlib import Path

PATH = Path("tensorcad/dg_handoff.py")
text = PATH.read_text(encoding="utf-8")

if "ellipse_contract" in text:
    raise SystemExit("already patched")

# ---------------------------------------------------------------- exports
text = text.replace(
    '    "circle_contract",\n',
    '    "circle_contract",\n    "ellipse_contract",\n    "ellipse_perimeter",\n',
    1,
)

# ---------------------------------------------------------------- frame import
text = text.replace(
    "from dataclasses import dataclass\nfrom pathlib import Path",
    "from dataclasses import dataclass\nfrom pathlib import Path",
    1,
)

# ------------------------------------------------- ellipse contract + perimeter
ANCHOR = "# ----------------------------------------------------------------- measurement"
ELLIPSE = '''# ----------------------------------------------------------------- conics


def ellipse_perimeter(a: float, b: float, *, order: int = 24) -> float:
    """Perimeter of the ellipse with semi-axes ``a``, ``b`` by ``4 a E(e)``.

    The reference is deliberately derived in a **different parametrisation** from
    the rational-quadratic curve this module hands to TensorDG: here the ellipse
    is ``(A cos t, B sin t)`` with ``A = max(a, b)``, so the speed is
    ``A sqrt(1 - e^2 sin^2 t)`` and the perimeter is ``4 A E(e)``, the complete
    elliptic integral of the second kind evaluated by Gauss-Legendre quadrature
    on ``t in [0, pi/2]`` (a smooth integrand: 24 nodes reach machine precision).
    If the NURBS construction were wrong, this number would not notice.
    """
    A = max(float(a), float(b))
    B = min(float(a), float(b))
    if not (A > 0.0 and math.isfinite(A) and math.isfinite(B)):
        raise HandoffError("semi-axes must be positive and finite")
    e2 = 1.0 - (B / A) ** 2
    nodes, weights = _gauss_legendre(int(order))
    integral = 0.0
    for x, w in zip(nodes, weights):
        theta = 0.25 * math.pi * (x + 1.0)
        s = math.sin(theta)
        integral += w * math.sqrt(max(0.0, 1.0 - e2 * s * s))
    return 4.0 * A * 0.25 * math.pi * integral


def ellipse_contract(
    a: float, b: float, centre: Sequence[float] = (0.0, 0.0), *,
    segments: int = 8, plane: str = "xy", offset_m: float = 0.0,
    name: str = "cad-ellipse",
) -> dict[str, Any]:
    """Rational-quadratic NURBS **ellipse**, exact for any ``a``, ``b``.

    Built from :func:`circle_contract` by an affine map of the control points.
    That is not an approximation: a NURBS curve is affinely invariant (Piegl &
    Tiller 4.2.2), so scaling the unit-circle control polygon by ``(a, b)`` and
    keeping the tangent-intersection weights ``cos(dtheta/2)`` produces the
    ellipse exactly, which is what lets a CAD section that is not circular reach
    TensorDG without being re-fitted.

    A circle is the special case ``a == b``; ``radius_m`` is kept in the payload
    as the **enclosing** radius ``max(a, b)`` (what a ring-shaped DG domain must
    clear), not as a semi-axis, and ``semi_axes`` is the authoritative geometry.
    """
    ax, by = float(a), float(b)
    if not (ax > 0.0 and by > 0.0 and math.isfinite(ax) and math.isfinite(by)):
        raise HandoffError("semi-axes must be positive and finite")
    cx, cy = float(centre[0]), float(centre[1])
    base = circle_contract(1.0, (0.0, 0.0), segments=segments, plane=plane,
                           offset_m=offset_m, name=name)
    contract = dict(base)
    contract["control_points"] = [[cx + ax * p[0], cy + by * p[1]]
                                  for p in base["control_points"]]
    contract["centre"] = [cx, cy]
    contract["semi_axes"] = [ax, by]
    contract["radius_m"] = max(ax, by)
    contract["radius_kind"] = "enclosing_max_semi_axis"
    contract["producer"] = (
        "tensorcad.dg_handoff.ellipse_contract (affine map of the tangent-intersection circle)")
    return contract


''' + ANCHOR
text = text.replace(ANCHOR, ELLIPSE, 1)

# ---------------------------------------------------- report: ellipse references
OLD_START = "    degree, knots, cps, weights = _check_contract(contract)\n    radius = float(contract.get(\"radius_m\", 0.0))"
NEW_START = '''    degree, knots, cps, weights = _check_contract(contract)
    semi = contract.get("semi_axes")
    ellipse: tuple[float, float] | None = None
    if semi is not None:
        if len(semi) != 2:
            raise HandoffError("semi_axes must hold two positive numbers")
        ax, ay = float(semi[0]), float(semi[1])
        if not (ax > 0.0 and ay > 0.0 and math.isfinite(ax) and math.isfinite(ay)):
            raise HandoffError("semi_axes must be positive and finite")
        if not contract.get("centre"):
            raise HandoffError("an ellipse contract must carry its centre")
        ellipse = (ax, ay)
        radius = max(ax, ay)  # enclosing radius, only used for the outer-domain check
        cx, cy = float(contract["centre"][0]), float(contract["centre"][1])
    radius = float(contract.get("radius_m", 0.0)) if ellipse is None else radius'''
text = text.replace(OLD_START, NEW_START, 1)

# the pre-existing radius/centre resolution block must not override the ellipse
OLD_BLOCK = '''    centre = contract.get("centre")
    if radius <= 0.0:'''
NEW_BLOCK = '''    centre = contract.get("centre")
    if ellipse is not None:
        pass
    elif radius <= 0.0:'''
text = text.replace(OLD_BLOCK, NEW_BLOCK, 1)

# measured loop: ellipse residual instead of radial distance
OLD_RADIAL = '''            worst_radial = max(worst_radial, abs(math.hypot(x - cx, y - cy) - radius))'''
NEW_RADIAL = '''            if ellipse is None:
                worst_radial = max(worst_radial, abs(math.hypot(x - cx, y - cy) - radius))
            else:
                # implicit residual of the ellipse, dimensionless: 0 on the curve
                rho = math.hypot((x - cx) / ellipse[0], (y - cy) / ellipse[1])
                worst_radial = max(worst_radial, abs(rho - 1.0))'''
text = text.replace(OLD_RADIAL, NEW_RADIAL, 1)

OLD_REFS = '''    analytic_circumference = 2.0 * math.pi * radius
    analytic_area = math.pi * radius * radius'''
NEW_REFS = '''    if ellipse is None:
        analytic_circumference = 2.0 * math.pi * radius
        analytic_area = math.pi * radius * radius
        deviation_m = worst_radial
        deviation_rel = worst_radial / radius
    else:
        analytic_circumference = ellipse_perimeter(*ellipse)
        analytic_area = math.pi * ellipse[0] * ellipse[1]
        deviation_m = worst_radial * min(ellipse)
        deviation_rel = worst_radial'''
text = text.replace(OLD_REFS, NEW_REFS, 1)

OLD_DICT = '''        "radius_m": radius,
        "centre_m": [cx, cy],
        "closure_gap_m": closure_gap,
        "max_radial_deviation_m": worst_radial,
        "max_radial_deviation_rel": worst_radial / radius,'''
NEW_DICT = '''        "radius_m": radius,
        "radius_kind": "enclosing_max_semi_axis" if ellipse is not None else "radius_m",
        "semi_axes": list(ellipse) if ellipse is not None else None,
        "centre_m": [cx, cy],
        "closure_gap_m": closure_gap,
        "max_radial_deviation_m": deviation_m,
        "max_radial_deviation_rel": deviation_rel,'''
text = text.replace(OLD_DICT, NEW_DICT, 1)

# ------------------------------------------------------------ frame declaration
OLD_MANIFEST = '''        "plane": plane_map.to_dict(),'''
NEW_MANIFEST = '''        "frame": dict(COORDINATE_CONTRACT),
        "plane": plane_map.to_dict(),'''
text = text.replace(OLD_MANIFEST, NEW_MANIFEST, 1)

OLD_CHECK = '''    digest = _contract_digest(contract)
    if digest != manifest["geometry"]["digest"]:'''
NEW_CHECK = '''    frame = manifest.get("frame")
    if frame is None:
        # pre-frame bundles exist; they are readable but cannot enter a cross-route
        # frame gate, and the reader is told so rather than guessing a frame
        raise HandoffError(
            "bundle declares no frame; a solver-facing handoff must carry the "
            f"{COORDINATE_CONTRACT['name']!r} contract frame")
    if frame.get("name") != COORDINATE_CONTRACT["name"]:
        raise HandoffError(
            f"bundle frame {frame.get('name')!r} is not the contract frame "
            f"{COORDINATE_CONTRACT['name']!r}")

    digest = _contract_digest(contract)
    if digest != manifest["geometry"]["digest"]:'''
text = text.replace(OLD_CHECK, NEW_CHECK, 1)

# ---------------------------------------------------------------- frame import
if "COORDINATE_CONTRACT" not in text.split("class HandoffError")[0]:
    text = text.replace(
        "from typing import Any, Iterable, Sequence\n",
        "from typing import Any, Iterable, Sequence\n\n"
        "from .solver_interface import COORDINATE_CONTRACT\n",
        1,
    )

PATH.write_text(text, encoding="utf-8")
print("patched", PATH)