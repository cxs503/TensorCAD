"""Series 60 hull family: the systematic displacement-hull benchmark.

The Series 60 (Todd 1951) is the most widely cited parametric hull family in naval
architecture. Five parent forms cover Cb = 0.60, 0.65, 0.70, 0.75, 0.80 with
closed-form hydrostatic coefficients for verification.

Geometry: y(xi, zeta) = (B/2) * (1 - |2*xi-1|^n) * zeta**(1/m)
where n, m are derived from target Cb via published Cp values.

Closed forms: Cp=n/(n+1), Cm=m/(m+1), Cb=Cp*Cm, Cwp=Cp, VCB/T=(m+1)/(2m+1)
"""
from __future__ import annotations

import math
from dataclasses import asdict, dataclass, field
from typing import Any

from . import mesh as meshlib
from .geometry import Mesh, MeshBuilder
from .hull import COORDINATE_FRAME, UNITS
from .hydrostatics import Hydrostatics, compute_hydrostatics
from .validation import (
    WARNING,
    Check,
    DimensionTarget,
    ValidationReport,
    check_dimensions,
    inspect_geometry,
    merge,
)

_PARENT_FORMS: dict[float, tuple[float, float]] = {
    # (Cp, n) where n = Cp/(1-Cp), m = Cb/(Cp-Cb)
    # Cm = Cb/Cp increases with Cb (U-sections get fuller)
    0.60: (0.625, 5.0 / 3.0),  # Cp=0.625, n=1.667, m=24.0
    0.65: (0.670, 0.670/(1-0.670)),  # Cp=0.670, n=2.030, m=32.5
    0.70: (0.715, 0.715/(1-0.715)),  # Cp=0.715, n=2.509, m=46.7
    0.75: (0.760, 0.760/(1-0.760)),  # Cp=0.760, n=3.167, m=75.0
    0.80: (0.805, 0.805/(1-0.805)),  # Cp=0.805, n=4.128, m=160.0
}

SAMPLING_ERROR_WARNING = 0.02
MAX_STATIONS = 4001
MAX_LEVELS = 513


def _interpolate_exponents(cb: float) -> tuple[float, float]:
    if not math.isfinite(cb) or not 0.55 <= cb <= 0.85:
        raise ValueError(f"Cb must lie in [0.55, 0.85], got {cb!r}")
    sorted_cb = sorted(_PARENT_FORMS.keys())
    if cb <= sorted_cb[0]:
        cp, n = _PARENT_FORMS[sorted_cb[0]]
    elif cb >= sorted_cb[-1]:
        cp, n = _PARENT_FORMS[sorted_cb[-1]]
    else:
        for lo, hi in zip(sorted_cb, sorted_cb[1:]):
            if lo <= cb <= hi:
                frac = (cb - lo) / (hi - lo)
                cp = _PARENT_FORMS[lo][0] + frac * (_PARENT_FORMS[hi][0] - _PARENT_FORMS[lo][0])
                n = _PARENT_FORMS[lo][1] + frac * (_PARENT_FORMS[hi][1] - _PARENT_FORMS[lo][1])
                break
    if cp <= cb:
        raise ValueError(f"Cp ({cp:.6f}) must exceed Cb ({cb:.6f})")
    m = cb / (cp - cb)
    return n, m


def analytic_volume(L, B, T, n, m):
    return L * B * T * n * m / ((n + 1.0) * (m + 1.0))

def analytic_waterplane_area(L, B, n):
    return L * B * n / (n + 1.0)

def analytic_maximum_section_area(B, T, m):
    return B * T * m / (m + 1.0)

def analytic_vcb_above_keel(T, m):
    """Vertical centre of buoyancy above the keel for the ``zeta**(1/m)`` section.

    ``int zeta*y dzeta / int y dzeta`` with ``y ~ zeta**(1/m)`` gives
    ``(m+1)/(2m+1)``, not ``(m+1)/(2(m+2))``: the previous expression was off by
    6.5 % at m=24 (0.481 T against 0.510 T).  Both are inside the "0 < VCB < T"
    sanity band, which is why only an independent recomputation catches it.
    """
    return T * (m + 1.0) / (2.0 * m + 1.0)

def analytic_vcb_z(T, m):
    return analytic_vcb_above_keel(T, m) - T

def analytic_cb(n, m):
    return n * m / ((n + 1.0) * (m + 1.0))

def analytic_cp(n):
    return n / (n + 1.0)

def analytic_cm(m):
    return m / (m + 1.0)

def analytic_cwp(n):
    return n / (n + 1.0)


@dataclass(frozen=True)
class Series60Parameters:
    length: float = 10.0
    beam: float = 1.5
    draft: float = 0.6
    nx: int = 41
    nz: int = 13
    water_density: float = 1025.0
    cb: float = 0.60

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)

    @property
    def exponents(self) -> tuple[float, float]:
        return _interpolate_exponents(self.cb)

    @property
    def n(self) -> float:
        return self.exponents[0]

    @property
    def m(self) -> float:
        return self.exponents[1]

    def problems(self) -> tuple[list[str], list[str]]:
        errors: list[str] = []
        for name, mx in (("length", 1.0e4), ("beam", 1.0e3), ("draft", 1.0e3)):
            v = float(getattr(self, name))
            if not math.isfinite(v) or v <= 0.0 or v > mx:
                errors.append(f"{name} must be finite and in (0, {mx:g}], got {v!r}")
        if not 3 <= int(self.nx) <= MAX_STATIONS:
            errors.append(f"nx must be in [3, {MAX_STATIONS}], got {self.nx!r}")
        if not 3 <= int(self.nz) <= MAX_LEVELS:
            errors.append(f"nz must be in [3, {MAX_LEVELS}], got {self.nz!r}")
        if not math.isfinite(self.water_density) or not 300.0 <= self.water_density <= 2000.0:
            errors.append(f"water_density must be in [300, 2000], got {self.water_density!r}")
        if not math.isfinite(self.cb) or not 0.55 <= self.cb <= 0.85:
            errors.append(f"cb must lie in [0.55, 0.85], got {self.cb!r}")
        warnings: list[str] = []
        if not errors:
            lb = self.length / self.beam
            bt = self.beam / self.draft
            if not 1.5 <= lb <= 15.0:
                warnings.append(f"length_to_beam = {lb:.3g} outside [1.5, 15]")
            if not 1.0 <= bt <= 12.0:
                warnings.append(f"beam_to_draft = {bt:.3g} outside [1.0, 12]")
        return errors, warnings

    def validate(self) -> None:
        errors, _ = self.problems()
        if errors:
            raise ValueError("; ".join(errors))

    @property
    def grid_is_symmetric(self) -> bool:
        return self.nx % 2 == 1


def _longitudinal_taper(p: Series60Parameters, i: int) -> float:
    if i == 0 or i == p.nx - 1:
        return 0.0
    xi = i / (p.nx - 1)
    return 1.0 - math.pow(abs(2.0 * xi - 1.0), p.n)


def _vertical_taper(p: Series60Parameters, level: int) -> float:
    """Section fullness at ``level`` (0 at the keel, ``nz-1`` at the waterline).

    ``zeta ** (1/m)``, not ``zeta ** m``.  The closed forms in this module --
    ``Cm = m/(m+1)`` and ``Cb = Cp*Cm``, which are what pins the exponents to the
    published Series 60 parent forms -- follow from ``int_0^1 zeta**(1/m) dzeta =
    m/(m+1)``.  With the exponent inverted the lofted sections collapse to a sliver
    over the bottom 90 % of the draft and the mesh integrates 21 times less volume
    than its own analytic reference; the error does not shrink under refinement
    (measured -95.8 % at 241x101), which is how it was found during the geometry
    audit (docs/series60.md records the measurement).
    """
    if level >= p.nz - 1:
        return 1.0
    if level <= 0:
        return 0.0
    return math.pow(level / (p.nz - 1), 1.0 / p.m)


def _section_point(p: Series60Parameters, i: int, j: int) -> tuple[float, float, float]:
    x = p.length * i / (p.nx - 1)
    if j <= p.nz - 1:
        level = p.nz - 1 - j
        b = 0.5 * p.beam * _longitudinal_taper(p, i) * _vertical_taper(p, level)
        z = -p.draft * (1.0 - level / (p.nz - 1))
        return (x, b, z)
    else:
        level = j - (p.nz - 1)
        b = 0.5 * p.beam * _longitudinal_taper(p, i) * _vertical_taper(p, level)
        z = -p.draft * (1.0 - level / (p.nz - 1))
        return (x, -b, z)


def build_series60_mesh(p: Series60Parameters) -> tuple[Mesh, tuple[tuple[float, float], ...]]:
    p.validate()
    ring_size = 2 * p.nz - 1
    builder = MeshBuilder(units=UNITS, coordinate_frame=COORDINATE_FRAME)
    rings: list[list[int]] = []
    for i in range(p.nx):
        rings.append([builder.add_vertex(_section_point(p, i, j)) for j in range(ring_size)])
    for i in range(p.nx - 1):
        for j in range(ring_size - 1):
            builder.add_quad(rings[i][j], rings[i + 1][j], rings[i + 1][j + 1], rings[i][j + 1])
        builder.add_quad(rings[i][ring_size - 1], rings[i + 1][ring_size - 1],
                         rings[i + 1][0], rings[i][0])
    mesh = meshlib.orient_outward(builder.build()).mesh
    offsets = tuple((p.length * i / (p.nx - 1), 0.5 * p.beam * _longitudinal_taper(p, i))
                    for i in range(p.nx))
    return mesh, offsets


def _dimension_target(p: Series60Parameters) -> DimensionTarget:
    return DimensionTarget(length=p.length, beam=p.beam, draft=p.draft)


def _sampling_check(p: Series60Parameters, hydro: Hydrostatics) -> Check:
    n, m = p.exponents
    ref_vol = analytic_volume(p.length, p.beam, p.draft, n, m)
    rel = (hydro.displacement_volume - ref_vol) / ref_vol if ref_vol > 0 else float("nan")
    ok = abs(rel) <= SAMPLING_ERROR_WARNING
    return Check("series60_volume_convergence", ok,
                 f"volume {hydro.displacement_volume:.6f} vs analytic {ref_vol:.6f} = {rel:+.4%}",
                 severity=WARNING)


@dataclass(frozen=True)
class Series60Model:
    parameters: Series60Parameters
    mesh: Mesh
    offsets: tuple[tuple[float, float], ...]
    validation: ValidationReport = field(default_factory=ValidationReport)
    hydrostatics: Hydrostatics | None = None

    @property
    def vertices(self): return self.mesh.vertices
    @property
    def faces(self): return self.mesh.faces
    @property
    def units(self): return self.mesh.units
    @property
    def coordinate_frame(self): return self.mesh.coordinate_frame

    def metadata(self, *, with_discretisation: bool = False) -> dict[str, Any]:
        n, m = self.parameters.exponents
        data: dict[str, Any] = {
            "schema_version": "2.0", "model_type": "series_60_parametric_displacement_hull",
            "parameters": self.parameters.as_dict(), "exponents": {"n": n, "m": m},
            "units": self.units, "coordinate_frame": self.coordinate_frame,
            "generator": {"kind": "series_60_lofted_surface",
                          "longitudinal_taper": f"1 - |2*xi - 1|^{n:.4f}",
                          "vertical_taper": f"zeta^{m:.4f}",
                          "topology": "one closed section ring per station, capped at WL",
                          "backend": "tensorcad.geometry (standard library)"},
            "mesh": self.mesh.metadata(), "offset_count": len(self.offsets),
            "validation": self.validation.to_dict(),
            "engineering_status": "geometry and hydrostatics only",
        }
        if self.hydrostatics is not None:
            data["hydrostatics"] = self.hydrostatics.to_dict()
        return data


def generate_series60(parameters: Series60Parameters | None = None,
                      **overrides: Any) -> Series60Model:
    if parameters is None:
        p = Series60Parameters(**overrides)
    elif overrides:
        p = Series60Parameters(**{**parameters.as_dict(), **overrides})
    else:
        p = parameters
    p.validate()
    mesh, offsets = build_series60_mesh(p)
    report = merge(inspect_geometry(mesh), check_dimensions(mesh, _dimension_target(p)))
    if report.errors:
        raise ValueError("generated Series 60 hull failed validation: " + "; ".join(report.errors))
    hydrostatics = compute_hydrostatics(mesh, 0.0, water_density=p.water_density, closed=True)
    report = merge(list(report.checks), [_sampling_check(p, hydrostatics)])
    return Series60Model(parameters=p, mesh=mesh, offsets=offsets,
                         validation=report, hydrostatics=hydrostatics)


def analytic_series60_hydrostatics(p: Series60Parameters) -> dict[str, float]:
    n, m = p.exponents
    return {
        "displacement_volume": analytic_volume(p.length, p.beam, p.draft, n, m),
        "waterplane_area": analytic_waterplane_area(p.length, p.beam, n),
        "maximum_section_area": analytic_maximum_section_area(p.beam, p.draft, m),
        "vertical_centre_buoyancy": analytic_vcb_above_keel(p.draft, m),
        "vertical_centre_buoyancy_z": analytic_vcb_z(p.draft, m),
        "block_coefficient": analytic_cb(n, m),
        "prismatic_coefficient": analytic_cp(n),
        "midship_coefficient": analytic_cm(m),
        "waterplane_coefficient": analytic_cwp(n),
        "longitudinal_centre_buoyancy": p.length / 2.0,
        "longitudinal_centre_flotation": p.length / 2.0,
    }


def series60_analytic_relative_errors(p: Series60Parameters, hydro: Hydrostatics) -> dict[str, float]:
    ref = analytic_series60_hydrostatics(p)
    errors: dict[str, float] = {}
    for key in ("displacement_volume", "waterplane_area", "maximum_section_area",
                "vertical_centre_buoyancy", "longitudinal_centre_buoyancy",
                "longitudinal_centre_flotation", "block_coefficient",
                "prismatic_coefficient", "midship_coefficient", "waterplane_coefficient"):
        rv = ref.get(key)
        mv = getattr(hydro, key, None) if hydro else None
        if rv is not None and mv is not None and abs(rv) > 1e-15:
            errors[key] = (mv - rv) / rv
        else:
            errors[key] = float("nan")
    return errors


def estimate_series60_discretisation_error(p: Series60Parameters,
                                           max_triangles: int = 2_000_000) -> dict[str, Any]:
    ref = analytic_series60_hydrostatics(p)
    ref_vol, ref_awp = ref["displacement_volume"], ref["waterplane_area"]
    results: list[dict[str, Any]] = []
    nx, nz = p.nx, p.nz
    while True:
        pp = Series60Parameters(length=p.length, beam=p.beam, draft=p.draft,
                                nx=nx, nz=nz, water_density=p.water_density, cb=p.cb)
        model = generate_series60(pp)
        tc = len(model.mesh.faces)
        ve = (model.hydrostatics.displacement_volume - ref_vol) / ref_vol if ref_vol else float("nan")
        ae = (model.hydrostatics.waterplane_area - ref_awp) / ref_awp if ref_awp else float("nan")
        results.append({"nx": nx, "nz": nz, "triangles": tc,
                        "volume_m3": model.hydrostatics.displacement_volume,
                        "volume_error": ve, "waterplane_area_m2": model.hydrostatics.waterplane_area,
                        "waterplane_error": ae})
        if tc > max_triangles or nx > MAX_STATIONS // 2:
            break
        nx = 2 * nx - 1
        nz = 2 * nz - 1
    return {"reference": ref, "convergence": results}
