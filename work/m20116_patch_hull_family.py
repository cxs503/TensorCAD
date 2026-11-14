"""Apply the hull-family dispatch edits to tensorcad/hull.py (one-shot patch).

Every replacement is exact-string with a uniqueness assertion, so a shifted
context refuses rather than silently patching the wrong place.
"""
from pathlib import Path

HULL = Path("tensorcad/hull.py")
text = HULL.read_text()
original = text


def sub(old: str, new: str, count: int = 1) -> None:
    global text
    found = text.count(old)
    assert found == count, f"pattern found {found} times (expected {count}): {old[:70]!r}"
    text = text.replace(old, new, count)


# ---------------------------------------------------------------- 1. constants
sub(
    """MAX_LONGITUDINAL_STATIONS = 4001""",
    '''#: Hull families this module can generate. ``sine`` is the shipped analytic
#: surface (``sin^p`` longitudinal taper, ``sqrt`` vertical profile, one exponent
#: per half). ``series60`` delegates to :mod:`tensorcad.series60` for the Todd
#: (1951) station law with its Cb ladder, flat transom and bulbous bow, so a
#: caller sees one parameter set, one validation report and one hydrostatics
#: block whichever family it asks for.
FAMILIES: tuple[str, ...] = ("sine", "series60")

MAX_LONGITUDINAL_STATIONS = 4001''',
)

# ------------------------------------------------------------ 2. HullParameters
sub(
    """    fore_exponent: float = BREADTH_EXPONENT
    aft_exponent: float = BREADTH_EXPONENT
""",
    """    fore_exponent: float = BREADTH_EXPONENT
    aft_exponent: float = BREADTH_EXPONENT
    #: Which station law to loft. ``"sine"`` keeps the shipped surface bit for
    #: bit; ``"series60"`` routes the same principal dimensions through
    #: :mod:`tensorcad.series60`.
    family: str = "sine"
    #: --- Series 60 only (refused, not ignored, on ``family="sine"``) ---
    #: Target block coefficient of the Todd (1951) parent form, 0.55-0.85.
    block_coefficient: float | None = None
    #: Waterplane coefficient. ``None`` keeps ``Cwp == Cp`` (the legacy single
    #: exponent, which is -10.7 % on Cwp at Cb = 0.60 against the parents); a value
    #: decouples the deck plan from the section-area curve.
    waterplane_coefficient: float | None = None
    #: LCB over length from the aft end (0.5 = midship). Requires
    #: ``waterplane_coefficient``.
    lcb_ratio: float | None = None
    #: Aft-end deck half-breadth over ``B/2``: how wide the immersed flat transom is.
    transom_ratio: float = 0.0
    #: Bulb semi-axis ahead of the forward perpendicular over the length.
    bulb_length_ratio: float = 0.0
    #: Bulb transverse area at the forward perpendicular over the midship section
    #: area: the Holtrop & Mennen (1984) ``A_BT``.
    bulb_area_ratio: float = 0.0
    #: Bulb section centre height above the keel over the draft.
    bulb_centre_height_ratio: float = 0.35
""",
)

# ----------------------------------------------------- 3. docstring of the class
sub(
    """    ``fore_exponent`` and ``aft_exponent`` are the sine-taper exponents forward of
    and abaft midship.""",
    """    ``family`` selects the station law: ``"sine"`` (default) is the surface this
    module has always lofted, ``"series60"`` is the Todd (1951) Series 60 family in
    :mod:`tensorcad.series60`, parameterised by ``block_coefficient``, the optional
    ``waterplane_coefficient``/``lcb_ratio`` coupling, ``transom_ratio`` and the
    bulbous-bow ratios. The Series 60 parameters are refused on the sine family
    instead of being quietly dropped: a caller that asked for a Cb and silently got
    an exponent has been lied to.

    ``fore_exponent`` and ``aft_exponent`` are the sine-taper exponents forward of
    and abaft midship.""",
)

# -------------------------------------------------------------- 4. problems()
sub(
    """        errors = validate(PARAMETER_SPECS, values)
        warnings: list[str] = []
        if not errors:
            warnings = check_ratios(RATIO_SPECS, {k: float(v) for k, v in values.items()})
        return errors, warnings""",
    """        errors = validate(PARAMETER_SPECS, values)
        warnings: list[str] = []
        if not errors:
            warnings = check_ratios(RATIO_SPECS, {k: float(v) for k, v in values.items()})
        errors.extend(_family_problems(self))
        return errors, warnings""",
)

# ------------------------------------------- 5. family helpers after the class
sub(
    """def _breadth_exponent(parameters: HullParameters, index: int) -> float:""",
    '''#: Knobs that only mean something to the Series 60 station law.
_SERIES60_ONLY_KNOBS: tuple[str, ...] = (
    "block_coefficient",
    "waterplane_coefficient",
    "lcb_ratio",
    "transom_ratio",
    "bulb_length_ratio",
    "bulb_area_ratio",
)


def _series60_parameters(parameters: HullParameters):
    """Map the unified parameters onto :class:`tensorcad.series60.Series60Parameters`.

    Imported lazily: :mod:`tensorcad.series60` imports this module for the shared
    coordinate frame, so a module-level import would be circular.
    """
    if parameters.family != "series60":
        raise ValueError(
            f"the Series 60 parameter mapping needs family='series60', got {parameters.family!r}"
        )
    if parameters.block_coefficient is None:
        raise ValueError(
            "family='series60' needs block_coefficient: the target Cb (0.55-0.85) is what "
            "picks the parent form, so there is no sensible default to fall back on"
        )
    from .series60 import Series60Parameters

    return Series60Parameters(
        length=parameters.length,
        beam=parameters.beam,
        draft=parameters.draft,
        nx=parameters.nx,
        nz=parameters.nz,
        water_density=parameters.water_density,
        cb=float(parameters.block_coefficient),
        cwp=parameters.waterplane_coefficient,
        lcb=parameters.lcb_ratio,
        transom_ratio=float(parameters.transom_ratio),
        bulb_length_ratio=float(parameters.bulb_length_ratio),
        bulb_area_ratio=float(parameters.bulb_area_ratio),
        bulb_centre_height_ratio=float(parameters.bulb_centre_height_ratio),
    )


def _family_problems(parameters: HullParameters) -> list[str]:
    """Family-specific validation: unknown families, crossed knobs, delegated checks."""
    if parameters.family not in FAMILIES:
        return [f"family must be one of {FAMILIES}, got {parameters.family!r}"]
    if parameters.family == "sine":
        touched = [
            name
            for name in _SERIES60_ONLY_KNOBS
            if getattr(parameters, name) not in (None, 0.0)
        ]
        if touched:
            return [
                f"{', '.join(touched)} only apply to family='series60'; set family='series60' "
                "or drop them (a silently ignored Cb is worse than a refusal)"
            ]
        return []
    errors: list[str] = []
    if (
        parameters.fore_exponent != BREADTH_EXPONENT
        or parameters.aft_exponent != BREADTH_EXPONENT
    ):
        errors.append(
            "fore_exponent/aft_exponent are sine-family parameters: the Series 60 station "
            "law derives its own exponents from the target Cb"
        )
    try:
        grid = _series60_parameters(parameters)
    except ValueError as exc:  # missing block_coefficient
        return errors + [str(exc)]
    delegated_errors, _delegated_warnings = grid.problems()
    errors.extend(f"series60: {message}" for message in delegated_errors)
    return errors


def _breadth_exponent(parameters: HullParameters, index: int) -> float:''',
)

# ------------------------------------------------------- 6. HullModel fields
sub(
    """    validation: ValidationReport = field(default_factory=ValidationReport)
    hydrostatics: Hydrostatics | None = None
""",
    """    validation: ValidationReport = field(default_factory=ValidationReport)
    hydrostatics: Hydrostatics | None = None
    #: Closed mesh of the bulb appendage: a separate solid, not boolean-unioned
    #: with the hull (the displacement totals add its closed form).
    bulb_mesh: Mesh | None = None
    #: The delegated family's own report (station law, exponents, transom, bulb)
    #: when the surface came from another family module.
    family_report: dict[str, Any] = field(default_factory=dict)
""",
)

# ------------------------------------------------------------- 7. metadata
sub(
    """            "model_type": "closed_sampled_parametric_displacement_hull",
            "parameters": self.parameters.as_dict(),""",
    """            "model_type": "closed_sampled_parametric_displacement_hull",
            "family": self.parameters.family,
            "parameters": self.parameters.as_dict(),""",
)
sub(
    """        if self.hydrostatics is not None:
            data["hydrostatics"] = self.hydrostatics.to_dict()
        if with_discretisation:""",
    """        if self.hydrostatics is not None:
            data["hydrostatics"] = self.hydrostatics.to_dict()
        if self.family_report:
            data["family_report"] = self.family_report
        if self.bulb_mesh is not None:
            data["appendages"] = {
                "bodies": 2,
                "bulb_mesh": self.bulb_mesh.metadata(),
                "note": "the bulb is a separate closed solid; displacement totals add "
                        "its closed form (see family_report.appendages.bulb)",
            }
        if with_discretisation:""",
)

# ------------------------------------------------------- 8. build_mesh branch
sub(
    """    parameters.validate()
    p = parameters
    ring_size = 2 * p.nz - 1""",
    """    parameters.validate()
    p = parameters
    if p.family != "sine":
        return _build_family_mesh(p)
    ring_size = 2 * p.nz - 1""",
)

# ----------------------------------- 9. family mesh + model builders + branch
sub(
    """def generate_hull(parameters: HullParameters | None = None) -> HullModel:""",
    '''def _build_family_mesh(
    parameters: HullParameters,
) -> tuple[Mesh, tuple[tuple[float, float], ...]]:
    """Loft a non-sine family and convert its offsets to this module's convention.

    :mod:`tensorcad.series60` reports half-breadths; ``HullModel.offsets`` is the
    full girth width, so the factor of two is applied here rather than being
    duplicated (or forgotten) by every caller.
    """
    from .series60 import build_series60_mesh

    mesh, half_breadths = build_series60_mesh(_series60_parameters(parameters))
    return mesh, tuple((x, 2.0 * half) for x, half in half_breadths)


def _generate_family_hull(parameters: HullParameters) -> HullModel:
    """Generate a delegated family, then present it through this module's API."""
    from .series60 import generate_series60

    p = parameters
    grid = _series60_parameters(p)
    model = generate_series60(grid)
    hydrostatics = model.hydrostatics
    checks = list(model.validation.checks)
    realised_cb = (
        hydrostatics.block_coefficient
        if hydrostatics is not None and hydrostatics.block_coefficient
        else float("nan")
    )
    checks.append(
        Check(
            "family_dispatch",
            True,
            f"family='series60' routed to tensorcad.series60 with target Cb={grid.cb:g}",
        )
    )
    checks.append(
        Check(
            "series60_block_coefficient_measured",
            math.isfinite(realised_cb) and abs(realised_cb - grid.cb) <= SAMPLING_ERROR_WARNING * grid.cb,
            f"mesh Cb {realised_cb:.6f} against target {grid.cb:g} "
            f"({100.0 * (realised_cb / grid.cb - 1.0):+.2f} %); the deficit is the loft's "
            "first-order resolution of the steep section near the keel, see "
            "estimate_discretisation_error",
            severity=WARNING,
        )
    )
    report = merge(checks, [])
    metadata = model.metadata()
    family_report = {
        key: metadata[key]
        for key in ("station_law", "exponents", "transom", "appendages")
        if key in metadata
    }
    return HullModel(
        parameters=p,
        mesh=model.mesh,
        offsets=tuple((x, 2.0 * half) for x, half in model.offsets),
        validation=report,
        hydrostatics=hydrostatics,
        bulb_mesh=model.bulb_mesh,
        family_report=family_report,
    )


def generate_hull(parameters: HullParameters | None = None) -> HullModel:''',
)

sub(
    """    p = parameters if parameters is not None else HullParameters()
    mesh, offsets = build_mesh(p)""",
    """    p = parameters if parameters is not None else HullParameters()
    p.validate()
    if p.family != "sine":
        return _generate_family_hull(p)
    mesh, offsets = build_mesh(p)""",
)

# ------------------------------- 10. relative errors: family-aware reference
sub(
    """    p = parameters
    mesh, _ = build_mesh(p)
    waterline = mesh.bounding_box().maximum[2]
    measured = compute_hydrostatics(mesh, waterline, water_density=p.water_density, closed=True)
    reference = analytic_hydrostatics(""",
    """    p = parameters
    if p.family != "sine":
        return _family_relative_errors(p)
    mesh, _ = build_mesh(p)
    waterline = mesh.bounding_box().maximum[2]
    measured = compute_hydrostatics(mesh, waterline, water_density=p.water_density, closed=True)
    reference = analytic_hydrostatics(""",
)

sub(
    """def _sampling_check(parameters: HullParameters, hydrostatics: Hydrostatics) -> Check:""",
    '''def _family_relative_errors(parameters: HullParameters) -> dict[str, Any]:
    """Same signed table as :func:`analytic_relative_errors`, Series 60 reference.

    The reference is the resolved station law's closed form (the Todd ladder that
    fixed the exponents), so it is a genuine mesh-vs-analytic comparison of the
    loft, not a repeat of the law's own arithmetic.
    """
    from .series60 import analytic_series60_hydrostatics

    p = parameters
    mesh, _ = build_mesh(p)
    waterline = mesh.bounding_box().maximum[2]
    measured = compute_hydrostatics(mesh, waterline, water_density=p.water_density, closed=True)
    reference = analytic_series60_hydrostatics(_series60_parameters(p))
    # ``vertical_centre_buoyancy`` in the reference dict is measured up from the
    # keel; Hydrostatics reports z relative to the design waterline.
    names = (
        "displacement_volume",
        "waterplane_area",
        "maximum_section_area",
        "longitudinal_centre_buoyancy",
        "longitudinal_centre_flotation",
        "vertical_centre_buoyancy",
        "block_coefficient",
        "prismatic_coefficient",
    )
    errors: dict[str, float] = {}
    for name in names:
        key = "vertical_centre_buoyancy_z" if name == "vertical_centre_buoyancy" else name
        exact = float(reference[key])
        value = float(getattr(measured, name))
        errors[name] = (value - exact) / abs(exact) if exact else float("nan")
    estimate = estimate_discretisation_error(p)
    reported = estimate.get("relative_change", {}) if estimate.get("status") == "measured" else {}
    outside = {
        name: {"measured_relative_error": errors[name], "reported_relative_change": reported[name]}
        for name in reported
        if math.isfinite(errors.get(name, float("nan")))
        and abs(errors[name]) > reported[name]
    }
    return {
        "status": "measured",
        "family": p.family,
        "method": "mesh integrals against the closed forms of the resolved Series 60 station law",
        "grid": {"nx": p.nx, "nz": p.nz},
        "measured": {name: float(getattr(measured, name)) for name in names},
        "reference": {name: float(reference[
            "vertical_centre_buoyancy_z" if name == "vertical_centre_buoyancy" else name
        ]) for name in names},
        "relative_error": errors,
        "max_abs_relative_error": max(abs(value) for value in errors.values()),
        "outside_discretisation_estimate": outside,
        "reference_basis": "closed forms of the Todd (1951) Series 60 parent form "
                           "(station law resolved from the target Cb)",
    }


def _sampling_check(parameters: HullParameters, hydrostatics: Hydrostatics) -> Check:''',
)

# ------------------------- 11. discretisation estimate: keep every parameter
sub(
    """    refined = HullParameters(
        p.length,
        p.beam,
        p.draft,
        refined_nx,
        refined_nz,
        p.water_density,
        fore_exponent=p.fore_exponent,
        aft_exponent=p.aft_exponent,
    )""",
    """    # ``replace`` carries the family and its knobs, so the refinement measures the
    # same surface at a finer grid instead of a different form.
    refined = replace(p, nx=refined_nx, nz=refined_nz)""",
)

# ------------------------------- 12. resolve_resolution: family-aware reference
sub(
    """        model = generate_hull(current)
        exact = analytic_volume(
            current.length,
            current.beam,
            current.draft,
            fore_exponent=current.fore_exponent,
            aft_exponent=current.aft_exponent,
        )""",
    """        model = generate_hull(current)
        if current.family == "sine":
            exact = analytic_volume(
                current.length,
                current.beam,
                current.draft,
                fore_exponent=current.fore_exponent,
                aft_exponent=current.aft_exponent,
            )
        else:
            from .series60 import analytic_series60_hydrostatics

            exact = analytic_series60_hydrostatics(_series60_parameters(current))[
                "displacement_volume"
            ]""",
)

assert text != original
HULL.write_text(text)
print("patched", HULL, len(text.splitlines()), "lines")