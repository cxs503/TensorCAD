"""The repository's own surface ``Mesh`` into the canonical quality library.

Why this module exists
----------------------
:mod:`tensorcad.meshquality` is the canonical metric library (Verdict/Gmsh definitions,
six cell kinds, a gate) and :mod:`tensorcad.gate` is the production geometry gate every
hull export runs through. Until this module existed the two had no shared type: the gate
spoke a *fourth* aspect-ratio dialect, ``l_max**2 / A`` (equilateral = 2.3094), computed in
its own loop over ``mesh.faces``, and nothing in the pipeline could ask the canonical
library a question. A library nobody calls is a document, not a gate.

The four dialects this tree actually computes, all called "aspect ratio":

=================  ==========================  ============  ===================================
dialect            expression                  equilateral   produced by
=================  ==========================  ============  ===================================
``gate``           ``l_max**2 / A``            2.309401      :func:`tensorcad.gate._facet_aspect`
``legacy``         ``l_max**2 / (2 A)``        1.154701      ``validation``/``mesh_conditioning``
``edge_form``      ``l_max p / (4 sqrt3 A)``   1.0           ``geomquality.quality``
``verdict``        ``R / (2 r)``               1.0           ``meshquality.metrics`` (VTK/Verdict)
=================  ==========================  ============  ===================================

The ``gate`` dialect is exactly ``2 x legacy`` -- an identity, not a bound -- because the
only difference between the two expressions is the factor 2 in the denominator. That is
worth stating in code rather than in prose: a 2.0 difference between two numbers that both
call themselves an aspect ratio has already caused one threshold to be copied across
dialects in this repository (finding F1, ``docs/mesh-quality-lib-audit-m17404.md``).

What the bridge gives you
-------------------------
* :func:`qualify` -- a :class:`CanonicalMeshQuality` for any object with ``.vertices`` and
  ``.faces`` (i.e. :class:`tensorcad.geometry.Mesh`): every canonical metric per face, the
  four aspect dialects with their transfer factors, the surface topology block, and the
  index of the worst face so a CI failure names a location and not just a number.
* :func:`canonical_stats` -- the handful of scalars :func:`tensorcad.gate.measure_gated`
  records, so the production gate reports canonical, comparable numbers.
* :func:`legacy_to_canonical_limit` / :func:`gate_to_canonical_limit` -- the limit a
  canonical gate needs so that it can never pass a mesh the old limit rejects. The proof is
  :func:`tensorcad.meshquality.dialects.conservative_limit` (``verdict/legacy >= 1/2`` for
  every triangle); this module only composes it with the ``gate = 2 x legacy`` identity.

Deliberate limitations, stated rather than hidden: faces must be triangles (the repository's
``Mesh`` enforces that in ``__post_init__``), so quadrilaterals raise instead of being
silently fanned; a zero-area face is reported with ``area == 0`` and is not dropped from the
statistics; and the external oracles (Gmsh, VTK/meshmagick) are still absent on this host,
which changes nothing here because the bridge produces the *same* canonical numbers the rest
of the library already reports.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Mapping, Sequence

from .dialects import aspect_dialects, conservative_limit
from .metrics import ElementQuality, tri_metrics
from .report import SummaryStats, summary_stats
from .topology import SurfaceTopology, surface_topology

Vec3 = tuple[float, float, float]

#: The per-face metrics the bridge summarises. Anything else in ``ElementQuality`` is
#: reachable through :attr:`CanonicalMeshQuality.elements` but is not gated by default.
CANONICAL_METRICS: tuple[str, ...] = (
    "aspect_ratio",
    "min_sicn",
    "scaled_jacobian",
    "min_angle",
    "max_angle",
    "skewness",
    "gamma",
    "eta",
    "edge_ratio",
    "area",
)

#: The aspect dialects, in the order the table in the module docstring lists them.
ASPECT_DIALECTS_IN_BRIDGE: tuple[str, ...] = ("gate", "legacy", "edge_form", "verdict")


def gate_facet_aspect(area: float, longest: float) -> float:
    """``l_max**2 / A`` -- the expression :func:`tensorcad.gate._facet_aspect` maximises.

    Written with the same operation order as the gate's own loop (``longest / area`` where
    ``longest`` is already squared) so a delegation is bit-for-bit identical rather than
    merely close, and ``inf`` for a zero-area facet, which is what the gate's
    ``if area <= 0.0: continue`` skip means for the worst-case value.
    """
    if area <= 0.0:
        return math.inf
    return longest / area


def legacy_to_canonical_limit(limit: float, *, source: str = "legacy") -> float:
    """A ``verdict``-dialect limit no weaker than ``limit`` in ``source`` dialect.

    ``verdict / legacy >= 1/2`` for every triangle (proved in
    :mod:`tensorcad.meshquality.dialects`: ``verdict = legacy * alpha * beta`` with
    ``alpha * beta >= 1/2``), so ``verdict <= legacy_limit / 2`` implies
    ``legacy <= legacy_limit``. ``source="gate"`` first divides by the exact factor 2
    relating the gate and legacy expressions.
    """
    if source == "gate":
        limit = limit / 2.0
        source = "legacy"
    if source not in ("legacy", "edge_form", "verdict"):
        raise ValueError(f"unknown source dialect {source!r}; use 'gate', 'legacy', 'edge_form' or 'verdict'")
    return conservative_limit(limit, target="verdict", source=source)


def gate_to_canonical_limit(limit: float) -> float:
    """Shorthand for :func:`legacy_to_canonical_limit` with ``source="gate"``."""
    return legacy_to_canonical_limit(limit, source="gate")


@dataclass(frozen=True)
class CanonicalMeshQuality:
    """The canonical library's view of one surface mesh, in one record."""

    kind: str
    n_points: int
    elements: ElementQuality
    dialects: Mapping[str, tuple[float, ...]]
    topology: SurfaceTopology
    units: str = "m"

    # -- summaries -------------------------------------------------------------------
    def statistics(self) -> dict[str, SummaryStats]:
        """Outlier-resistant statistics per canonical metric, in face order."""
        return {metric: summary_stats(metric, self.elements[metric]) for metric in CANONICAL_METRICS}

    def stats_dict(self) -> dict[str, dict[str, float]]:
        return {metric: stats.to_dict() for metric, stats in self.statistics().items()}

    def scalars(self) -> dict[str, float]:
        """The scalars a gate compares: worst value of each metric, in its own direction."""
        stats = self.statistics()
        return {
            "max_aspect_ratio": stats["aspect_ratio"].maximum,
            "max_edge_ratio": stats["edge_ratio"].maximum,
            "max_skewness": stats["skewness"].maximum,
            "max_max_angle": stats["max_angle"].maximum,
            "min_angle": stats["min_angle"].minimum,
            "min_sicn": stats["min_sicn"].minimum,
            "min_scaled_jacobian": stats["scaled_jacobian"].minimum,
            "min_gamma": stats["gamma"].minimum,
            "min_eta": stats["eta"].minimum,
            "non_positive_area": float(stats["area"].non_positive),
        }

    def dialect_max(self) -> dict[str, float]:
        """Worst value of each aspect dialect over the faces (``inf`` if any is ``inf``)."""
        out: dict[str, float] = {}
        for dialect in ASPECT_DIALECTS_IN_BRIDGE:
            values = self.dialects[dialect]
            out[dialect] = max(values) if values else 0.0
        return out

    # -- localisation ----------------------------------------------------------------
    def worst(self, metric: str, *, maximum: bool = True) -> dict[str, Any]:
        """Index, value and centroid of the worst face for one metric.

        A gate that reports only an extremum cannot be acted on; the point of this method
        is that a CI failure can name the triangle to look at.
        """
        values = self.elements[metric]
        if not values:
            return {"metric": metric, "index": None, "value": math.nan, "centroid": None}
        target = math.inf if maximum else -math.inf
        index = None
        for position, value in enumerate(values):
            if value != value:  # nan never wins: an undefined metric is not the worst one
                continue
            if index is None or (value > target if maximum else value < target):
                index, target = position, value
        if index is None:
            return {"metric": metric, "index": None, "value": math.nan, "centroid": None}
        return {
            "metric": metric,
            "index": index,
            "value": target,
            "centroid": list(self.centroids[index]),
        }

    #: Set by :func:`qualify`; kept out of ``__init__`` so the record stays serialisable.
    centroids: tuple[Vec3, ...] = ()

    # -- serialisation ---------------------------------------------------------------
    def to_dict(self) -> dict[str, Any]:
        return {
            "kind": self.kind,
            "units": self.units,
            "points": self.n_points,
            "faces": len(self.elements),
            "statistics": self.stats_dict(),
            "aspect_dialects": {
                "maximum": self.dialect_max(),
                "identity_gate_over_legacy": 2.0,
            },
            "scalars": self.scalars(),
            "topology": self.topology.to_dict(),
            "worst": {
                "aspect_ratio": self.worst("aspect_ratio"),
                "min_sicn": self.worst("min_sicn", maximum=False),
                "min_angle": self.worst("min_angle", maximum=False),
                "skewness": self.worst("skewness"),
            },
        }


def _as_points(vertices: Sequence[Sequence[float]]) -> tuple[Vec3, ...]:
    return tuple((float(v[0]), float(v[1]), float(v[2])) for v in vertices)


def _as_triangles(faces: Sequence[Sequence[int]]) -> tuple[tuple[int, int, int], ...]:
    out: list[tuple[int, int, int]] = []
    for index, face in enumerate(faces):
        if len(face) != 3:
            raise ValueError(
                f"face {index} has {len(face)} vertices; the bridge converts the repository's "
                "triangle Mesh only -- fan a quad before handing it over, so the triangulation "
                "choice stays visible instead of being made silently here"
            )
        out.append((int(face[0]), int(face[1]), int(face[2])))
    return tuple(out)


def qualify(mesh: Any, *, backend: str = "auto") -> CanonicalMeshQuality:
    """Canonical metrics, aspect dialects and topology of one ``Mesh``-like object.

    ``mesh`` needs ``.vertices`` and ``.faces`` only, so a
    :class:`tensorcad.geometry.Mesh`, a mesh that came back through an importer, or a test
    fixture all qualify. ``backend`` is forwarded to :func:`tri_metrics`
    (``"auto"``/``"scalar"``/``"numpy"``).
    """
    points = _as_points(mesh.vertices)
    triangles = _as_triangles(mesh.faces)
    elements = tri_metrics(points, triangles, backend=backend)
    dialects = aspect_dialects(points, triangles)
    legacy = dialects["legacy"]
    gate_dialect = tuple(
        gate_facet_aspect(elements["area"][i], elements["max_edge"][i] ** 2) for i in range(len(elements))
    )
    merged = {**dict(dialects), "gate": gate_dialect}
    centroids: list[Vec3] = []
    for a, b, c in triangles:
        pa, pb, pc = points[a], points[b], points[c]
        centroids.append(
            ((pa[0] + pb[0] + pc[0]) / 3.0, (pa[1] + pb[1] + pc[1]) / 3.0, (pa[2] + pb[2] + pc[2]) / 3.0)
        )
    _ = legacy  # kept as a named local so the identity below is not an accidental re-read
    return CanonicalMeshQuality(
        kind="tri",
        n_points=len(points),
        elements=elements,
        dialects=merged,
        topology=surface_topology(points, triangles),
        units=str(getattr(mesh, "units", "m")),
        centroids=tuple(centroids),
    )


def canonical_stats(mesh: Any, *, backend: str = "auto") -> dict[str, float]:
    """:meth:`CanonicalMeshQuality.scalars` for a mesh, without keeping the record."""
    return qualify(mesh, backend=backend).scalars()


def dialect_identities(mesh: Any) -> dict[str, float]:
    """The exact relations the four dialects satisfy on this mesh, as worst-case deviations.

    ``gate = 2 * legacy`` is an identity of the expressions; ``edge_form = legacy * alpha``
    and ``verdict = legacy * alpha * beta`` are identities too, with the per-face ``alpha``,
    ``beta`` from :func:`tensorcad.meshquality.dialects.aspect_dialects`. The returned
    numbers are the maximum absolute deviations, so a test can assert they are at rounding
    level rather than "small".
    """
    quality = qualify(mesh, backend="scalar")
    dialects = quality.dialects
    worst_gate = 0.0
    worst_edge = 0.0
    worst_verdict = 0.0
    for i in range(len(quality.elements)):
        legacy = dialects["legacy"][i]
        alpha = dialects["alpha"][i]
        beta = dialects["beta"][i]
        worst_gate = max(worst_gate, abs(dialects["gate"][i] - 2.0 * legacy))
        if math.isfinite(legacy) and math.isfinite(alpha) and math.isfinite(dialects["edge_form"][i]):
            worst_edge = max(worst_edge, abs(dialects["edge_form"][i] - legacy * alpha))
        if (
            math.isfinite(legacy)
            and math.isfinite(alpha * beta)
            and math.isfinite(dialects["verdict"][i])
        ):
            worst_verdict = max(worst_verdict, abs(dialects["verdict"][i] - legacy * alpha * beta))
    return {
        "gate_minus_2_legacy": worst_gate,
        "edge_form_minus_legacy_alpha": worst_edge,
        "verdict_minus_legacy_alpha_beta": worst_verdict,
    }


__all__ = [
    "ASPECT_DIALECTS_IN_BRIDGE",
    "CANONICAL_METRICS",
    "CanonicalMeshQuality",
    "canonical_stats",
    "dialect_identities",
    "gate_facet_aspect",
    "gate_to_canonical_limit",
    "legacy_to_canonical_limit",
    "qualify",
]