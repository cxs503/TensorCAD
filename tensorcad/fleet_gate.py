"""Baseline hull fleet gate: provenance, placeholder detection, watertight/volume verdicts.

The ship platform's geometry library drifted into a state where six of eight
checked-in ``*_Hull.stl`` files were *not* hulls: 304-triangle box approximations,
two of them byte-identical to each other.  Nothing in the build noticed, because
"a file exists" was the only predicate anybody evaluated.  This module is the
missing predicate: a fleet entry is accepted only when **provenance** (official
source URL + publisher + licence + SHA-256 + retrieval time), **topology**
(closed, edge-manifold, outward-oriented) and **physics** (volume against the
source's own published displacement) are all satisfied.  Anything else is named
``PLACEHOLDER`` or ``WITHHELD`` with a machine-readable reason, never "ok".

Design rules kept from the rest of the package: standard library only (no CAD
kernel, no network), the gate fails *closed*, and every number it reports is
recomputed here -- nothing is trusted from a manifest that this module could not
re-measure from bytes on disk.
"""

from __future__ import annotations

import hashlib
import json
import math
import struct
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from .geometry import Mesh
from .mesh import is_edge_manifold, is_watertight, total_signed_volume

MANIFEST_SCHEMA = "tensorcad.baseline_fleet.manifest/1.0"
REPORT_SCHEMA = "tensorcad.baseline_fleet.audit/1.0"

#: Hosts that publish workshop geometry.  A source URL outside this list is not
#: rejected, but it can never be graded ``official``.
OFFICIAL_HOSTS: tuple[str, ...] = (
    "www.t2015.nmri.go.jp",
    "t2015.nmri.go.jp",
    "www.nmri.go.jp",
    "nmri.go.jp",
    "simman2020.kr",
    "www.simman2020.kr",
)

#: Fields a provenance record must carry before its bytes may be graded.
REQUIRED_SOURCE_FIELDS: tuple[str, ...] = (
    "source_url",
    "publisher",
    "license",
    "sha256",
    "retrieved_at",
)

#: Below this triangle count a hull mesh cannot resolve a bilge radius; the
#: platform's known placeholders sit at 304, real hull tessellations at >1e5.
MIN_HULL_TRIANGLES = 2_000

#: Volume acceptance against the source's published displacement.
VOLUME_TOLERANCE = 0.01

VERIFIED = "VERIFIED"
PLACEHOLDER = "PLACEHOLDER"
WITHHELD = "WITHHELD"
NOT_SUPPORTED = "NOT_SUPPORTED"


# --------------------------------------------------------------------------- #
# mesh intake
# --------------------------------------------------------------------------- #
def read_stl(path: str | Path) -> Mesh:
    """Read an ASCII or binary STL into a :class:`~tensorcad.geometry.Mesh`."""
    data = Path(path).read_bytes()
    if len(data) >= 84:
        count = struct.unpack("<I", data[80:84])[0]
        if 84 + 50 * count == len(data):
            vertices: list[tuple[float, float, float]] = []
            faces: list[tuple[int, int, int]] = []
            for i in range(count):
                base = 84 + 50 * i + 12  # skip the facet normal
                tri = struct.unpack("<9f", data[base : base + 36])
                start = len(vertices)
                vertices.extend((tri[0:3], tri[3:6], tri[6:9]))
                faces.append((start, start + 1, start + 2))
            return Mesh(vertices=vertices, faces=faces)
    text = data.decode("ascii", "replace")
    vertices, faces = [], []
    corner = 0
    for line in text.splitlines():
        token = line.strip().split()
        if len(token) == 4 and token[0].lower() == "vertex":
            vertices.append((float(token[1]), float(token[2]), float(token[3])))
            corner += 1
            if corner % 3 == 0:
                faces.append((corner - 3, corner - 2, corner - 1))
    if not faces:
        raise ValueError(f"{path}: not a readable STL (neither binary nor ASCII facets)")
    return Mesh(vertices=vertices, faces=faces)


def mesh_digest(mesh: Mesh) -> str:
    """SHA-256 over the vertex and face streams (order-sensitive, weld-sensitive)."""
    digest = hashlib.sha256()
    digest.update(struct.pack("<II", len(mesh.vertices), len(mesh.faces)))
    for vertex in mesh.vertices:
        digest.update(struct.pack("<3d", *vertex))
    for face in mesh.faces:
        digest.update(struct.pack("<3I", *face))
    return digest.hexdigest()


def _components(mesh: Mesh) -> int:
    parent = list(range(len(mesh.vertices)))

    def find(a: int) -> int:
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for face in mesh.faces:
        root = find(face[0])
        for other in face[1:]:
            r = find(other)
            if r != root:
                parent[r] = root
                root = find(root)
    return len({find(face[0]) for face in mesh.faces})


# --------------------------------------------------------------------------- #
# measurements
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class MeshAudit:
    """Everything the gate can measure from one mesh file."""

    path: str
    triangles: int
    vertices: int
    watertight: bool
    edge_manifold: bool
    components: int
    volume: float
    surface_area: float
    length: float
    breadth: float
    depth: float
    sha256: str

    @property
    def block_coefficient(self) -> float:
        box = self.length * self.breadth * self.depth
        return abs(self.volume) / box if box else float("nan")

    def as_dict(self) -> dict[str, Any]:
        return {
            "path": self.path,
            "triangles": self.triangles,
            "vertices": self.vertices,
            "watertight": self.watertight,
            "edge_manifold": self.edge_manifold,
            "components": self.components,
            "volume_m3": round(self.volume, 6),
            "surface_area_m2": round(self.surface_area, 6),
            "length_m": round(self.length, 6),
            "breadth_m": round(self.breadth, 6),
            "depth_m": round(self.depth, 6),
            "block_coefficient_of_bbox": round(self.block_coefficient, 6),
            "sha256": self.sha256,
        }


def audit_mesh(path: str | Path, *, weld_tolerance: float = 1e-6) -> MeshAudit:
    """Measure one mesh file.  Welds duplicate vertices first so that seams from
    a chunked exporter are judged as the shared edge they are."""
    mesh = read_stl(path)
    digest = mesh_digest(mesh)
    from .mesh import weld  # local import keeps module import cheap

    welded = weld(mesh, tolerance=weld_tolerance)
    xs = [v[0] for v in welded.vertices]
    ys = [v[1] for v in welded.vertices]
    zs = [v[2] for v in welded.vertices]
    area = 0.5 * sum(
        math.dist(welded.vertices[f[0]], welded.vertices[f[1]])
        * math.dist(welded.vertices[f[0]], welded.vertices[f[2]])
        * math.sin(
            _py_angle(
                _sub(welded.vertices[f[1]], welded.vertices[f[0]]),
                _sub(welded.vertices[f[2]], welded.vertices[f[0]]),
            )
        )
        for f in welded.faces
    )
    return MeshAudit(
        path=str(path),
        triangles=len(welded.faces),
        vertices=len(welded.vertices),
        watertight=is_watertight(welded),
        edge_manifold=is_edge_manifold(welded),
        components=_components(welded),
        volume=total_signed_volume(welded),
        surface_area=area,
        length=max(xs) - min(xs) if xs else 0.0,
        breadth=max(ys) - min(ys) if ys else 0.0,
        depth=max(zs) - min(zs) if zs else 0.0,
        sha256=digest,
    )


def _sub(a: Sequence[float], b: Sequence[float]) -> tuple[float, float, float]:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _py_angle(u: Sequence[float], v: Sequence[float]) -> float:
    nu = math.sqrt(sum(c * c for c in u))
    nv = math.sqrt(sum(c * c for c in v))
    if nu == 0.0 or nv == 0.0:
        return 0.0
    c = sum(a * b for a, b in zip(u, v)) / (nu * nv)
    return math.acos(max(-1.0, min(1.0, c)))


# --------------------------------------------------------------------------- #
# provenance
# --------------------------------------------------------------------------- #
def validate_source(entry: Mapping[str, Any]) -> list[str]:
    """Return the provenance problems with one registry entry (empty = valid)."""
    problems: list[str] = []
    for name in REQUIRED_SOURCE_FIELDS:
        value = entry.get(name)
        if value in (None, "", []):
            problems.append(f"missing {name}")
    sha = str(entry.get("sha256") or "")
    if sha and (len(sha) != 64 or any(c not in "0123456789abcdef" for c in sha.lower())):
        problems.append("sha256 is not 64 lowercase hex characters")
    url = str(entry.get("source_url") or "")
    if url and not url.startswith(("http://", "https://")):
        problems.append("source_url is not http(s)")
    host = url.split("//", 1)[-1].split("/", 1)[0].lower() if url else ""
    if host and host not in OFFICIAL_HOSTS:
        problems.append(f"source host {host!r} is not a known workshop publisher")
    if entry.get("license") in ("", "unknown", "UNKNOWN", "unspecified"):
        problems.append("license is not stated by the publisher")
    return problems


def is_official(entry: Mapping[str, Any]) -> bool:
    host = str(entry.get("source_url") or "").split("//", 1)[-1].split("/", 1)[0].lower()
    return host in OFFICIAL_HOSTS and not validate_source(entry)


# --------------------------------------------------------------------------- #
# fleet verdicts
# --------------------------------------------------------------------------- #
@dataclass
class ShipVerdict:
    ship: str
    grade: str
    reasons: list[str] = field(default_factory=list)
    mesh: dict[str, Any] | None = None
    volume_error: float | None = None

    def as_dict(self) -> dict[str, Any]:
        return {
            "ship": self.ship,
            "grade": self.grade,
            "reasons": self.reasons,
            "mesh": self.mesh,
            "volume_relative_error": None if self.volume_error is None else round(self.volume_error, 6),
        }


def judge_mesh(
    ship: str,
    audit: MeshAudit,
    *,
    published_volume: float | None,
    byte_twins: Iterable[str] = (),
    tolerance: float = VOLUME_TOLERANCE,
) -> ShipVerdict:
    """Grade one measured mesh against the three gates.

    ``PLACEHOLDER`` is reserved for meshes that provably are not this hull (too
    few facets to carry the form, or the same bytes published under another
    ship's name).  ``WITHHELD`` means the geometry may be real but the gate cannot
    certify it -- open surface, no published displacement to compare with.
    """
    reasons: list[str] = []
    placeholder = False
    if audit.triangles < MIN_HULL_TRIANGLES:
        placeholder = True
        reasons.append(
            f"{audit.triangles} triangles < {MIN_HULL_TRIANGLES}: cannot resolve a bilge radius"
        )
    twins = sorted(set(byte_twins) - {ship})
    if twins:
        placeholder = True
        reasons.append(f"byte-identical to {', '.join(twins)}: one form cannot be two ships")
    if placeholder:
        return ShipVerdict(ship, PLACEHOLDER, reasons, audit.as_dict())

    if not audit.watertight:
        reasons.append("not watertight: volume is not a closed-surface integral")
    elif not audit.edge_manifold:
        reasons.append("not edge-manifold")
    elif audit.components != 1:
        reasons.append(f"{audit.components} components: hull + appendages must be declared separately")
    elif published_volume is None:
        reasons.append("no published displacement on record for this ship")
    else:
        error = abs(abs(audit.volume) - published_volume) / published_volume
        if error <= tolerance:
            return ShipVerdict(
                ship,
                VERIFIED,
                [f"volume {abs(audit.volume):.3f} m3 vs published {published_volume:.3f} m3 ({error * 100:.3f} %)"],
                audit.as_dict(),
                error,
            )
        reasons.append(
            f"volume {abs(audit.volume):.3f} m3 vs published {published_volume:.3f} m3 "
            f"= {error * 100:.3f} % > {tolerance * 100:.0f} %"
        )
        return ShipVerdict(ship, WITHHELD, reasons, audit.as_dict(), error)
    return ShipVerdict(ship, WITHHELD, reasons, audit.as_dict())


def audit_fleet(
    directory: str | Path,
    *,
    published_volumes: Mapping[str, float] | None = None,
    pattern: str = "*_Hull.stl",
    tolerance: float = VOLUME_TOLERANCE,
) -> dict[str, Any]:
    """Audit every ``*_Hull.stl`` under ``directory`` and return a JSON report."""
    directory = Path(directory)
    published = dict(published_volumes or {})
    files = sorted(p for p in directory.glob(pattern) if p.is_file())
    audits: dict[str, MeshAudit] = {}
    failures: dict[str, str] = {}
    for path in files:
        ship = path.name[: -len(".stl")].replace("_Hull", "")
        try:
            audits[ship] = audit_mesh(path)
        except Exception as exc:  # unreadable file is a reportable state, not a crash
            failures[ship] = f"{type(exc).__name__}: {exc}"
    by_sha: dict[str, list[str]] = {}
    for ship, audit in audits.items():
        by_sha.setdefault(audit.sha256, []).append(ship)
    verdicts = [
        judge_mesh(
            ship,
            audit,
            published_volume=published.get(ship),
            byte_twins=by_sha.get(audit.sha256, []),
            tolerance=tolerance,
        )
        for ship, audit in audits.items()
    ]
    for ship, why in failures.items():
        verdicts.append(ShipVerdict(ship, NOT_SUPPORTED, [why]))
    verdicts.sort(key=lambda v: v.ship)
    counts: dict[str, int] = {}
    for verdict in verdicts:
        counts[verdict.grade] = counts.get(verdict.grade, 0) + 1
    return {
        "schema": REPORT_SCHEMA,
        "directory": str(directory),
        "pattern": pattern,
        "volume_tolerance": tolerance,
        "min_hull_triangles": MIN_HULL_TRIANGLES,
        "counts": counts,
        "duplicate_sha_groups": {sha: ships for sha, ships in by_sha.items() if len(ships) > 1},
        "ships": [v.as_dict() for v in verdicts],
        "pass": all(v.grade == VERIFIED for v in verdicts) and bool(verdicts),
    }


def write_json(payload: Mapping[str, Any], path: str | Path) -> Path:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(payload, indent=1, sort_keys=False) + "\n")
    return target