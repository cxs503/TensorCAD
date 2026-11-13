"""Knowledge governance for TensorCAD: one ledger, live probes, two exports.

The problem this module exists for is not "we lack documentation".  TensorCAD
already carries hundreds of verified numbers -- closed forms in ``wigley.py``,
exact section areas in ``box_gz.py``, benchmark gates under
``benchmarks/verified/`` and audit measurements in ``validation_report.json``.
What it did *not* have is a single machine-checkable index saying **which claim
is backed by which number, and is that number still true after today's commit**.
Documentation rots silently; a merged refactor can move a number and every
document that quotes it stays green because nothing compares them.

So the ledger (``knowledge/sediment_line.json``) stores knowledge as *sediment*:
each entry is a claim plus the machines that can contradict it today.

* ``probes`` -- a live recomputation through a public callable.  If the code
  changes, the probe value changes and the guard notices.
* ``artifacts`` -- a value at a JSON pointer inside a checked-in artifact, plus
  the artifact's SHA-256.  If the artifact is regenerated, the digest drifts and
  the entry must be re-verified deliberately (never silently).
* ``docs`` -- a verbatim phrase the prose must still contain.  If the sentence
  quoting the number is edited away, the entry fails.

Two intake mouths and two export mouths (the "double-in / double-out" shape):

* **passive intake** -- ``load_ledger`` reads what verification already produced;
* **active intake** -- :func:`harvest_candidates` walks the benchmark/audit
  artifacts and lists numeric leaves that nobody has pinned yet, i.e. the next
  entries worth reviewing;
* **agent outlet** -- :func:`export_agent_index` writes one JSONL line per claim
  so an agent can retrieve "what is verified, to what tolerance, how to
  reproduce" without reading prose;
* **training/evaluation outlet** -- :func:`export_training_corpus` writes a
  deterministic geometry -> hydrostatics corpus (and its held-out eval split)
  from the closed forms, usable as physical-AI training data with exact labels.

Design rules kept from the rest of the core package: standard library only, no
network, no clock (artifacts are byte-identical across runs), and a guard that
fails *closed* -- an unreadable pointer or an unimportable probe is a failure,
not a skip.
"""

from __future__ import annotations

import hashlib
import importlib
import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

#: Schema tag written into exports and expected in the ledger file.
SCHEMA = "tensorcad.sediment/1.0"
#: Export schema tag.
EXPORT_SCHEMA = "tensorcad.knowledge.export/1.0"
#: Entry kinds.  Kept short and closed so a typo cannot create a new class.
KINDS: tuple[str, ...] = (
    "analytic_identity",
    "exact_section",
    "benchmark_gate",
    "audit_measurement",
)
#: Probes may only reach into the package itself; the ledger is not a code-execution gadget.
PROBE_MODULE_PREFIX = "tensorcad."
#: Verdicts.  ``PASS`` is the only success.
PASS = "PASS"
DRIFT = "DRIFT"
MISMATCH = "MISMATCH"
MISSING = "MISSING"
BROKEN = "BROKEN"


class KnowledgeError(ValueError):
    """Raised when a ledger is malformed or a probe spec is not admissible."""


# --------------------------------------------------------------------------- #
# small helpers
# --------------------------------------------------------------------------- #
def sha256_file(path: str | Path) -> str:
    """ASCII-hex SHA-256 of a file, streamed (artifacts can be large)."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 16), b""):
            digest.update(chunk)
    return digest.hexdigest()


def resolve_pointer(document: Any, pointer: str) -> Any:
    """RFC 6901 JSON pointer (``""`` = whole document, ``/a/b/0`` = walk).

    Raises :class:`KnowledgeError` rather than returning ``None``: a missing
    pointer means the artifact schema moved, which the guard must report.
    """
    if pointer in ("", "/"):
        return document
    if not pointer.startswith("/"):
        raise KnowledgeError(f"JSON pointer must start with '/': {pointer!r}")
    node = document
    for raw in pointer.split("/")[1:]:
        token = raw.replace("~1", "/").replace("~0", "~")
        if isinstance(node, list):
            try:
                node = node[int(token)]
            except (ValueError, IndexError) as exc:
                raise KnowledgeError(f"pointer {pointer!r}: bad list index {token!r}") from exc
        elif isinstance(node, Mapping):
            if token not in node:
                raise KnowledgeError(f"pointer {pointer!r}: key {token!r} absent")
            node = node[token]
        else:
            raise KnowledgeError(f"pointer {pointer!r}: {token!r} descends into a scalar")
    return node


def _within(actual: Any, expected: Any, rel_tol: float, abs_tol: float) -> bool:
    """Numeric comparison for floats/ints; exact equality for anything else."""
    if isinstance(expected, bool) or isinstance(actual, bool):
        return actual == expected
    if isinstance(expected, (int, float)) and isinstance(actual, (int, float)):
        if math.isnan(float(actual)) or math.isnan(float(expected)):
            return False
        if math.isinf(float(expected)) or math.isinf(float(actual)):
            return float(actual) == float(expected)
        return math.isclose(float(actual), float(expected), rel_tol=rel_tol, abs_tol=abs_tol)
    return actual == expected


def _require_mapping(value: Any, what: str, entry_id: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise KnowledgeError(f"{entry_id}: {what} must be an object, got {type(value).__name__}")
    return value


def _require_keys(mapping: Mapping[str, Any], keys: Sequence[str], what: str, entry_id: str) -> None:
    missing = [key for key in keys if key not in mapping]
    if missing:
        raise KnowledgeError(f"{entry_id}: {what} missing key(s) {missing}")


# --------------------------------------------------------------------------- #
# ledger model
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Probe:
    """A live recomputation: call ``module:function`` and compare the result."""

    call: str
    expected: Any
    args: tuple[Any, ...] = ()
    kwargs: Mapping[str, Any] = field(default_factory=dict)
    rel_tol: float = 1e-12
    abs_tol: float = 0.0
    index: int | None = None
    attr: str | None = None
    note: str = ""


@dataclass(frozen=True)
class ArtifactRef:
    """A numeric leaf inside a checked-in artifact, pinned by digest."""

    path: str
    json_pointer: str
    expected: Any
    sha256: str
    rel_tol: float = 1e-12
    abs_tol: float = 0.0
    note: str = ""


@dataclass(frozen=True)
class DocRef:
    """A verbatim phrase the prose must still contain (doc drift detector)."""

    path: str
    contains: str
    note: str = ""


@dataclass(frozen=True)
class SedimentEntry:
    """One claim plus every machine that can contradict it."""

    id: str
    title: str
    claim: str
    kind: str
    acceptance: str
    probes: tuple[Probe, ...] = ()
    artifacts: tuple[ArtifactRef, ...] = ()
    docs: tuple[DocRef, ...] = ()
    tags: tuple[str, ...] = ()

    @property
    def checks(self) -> int:
        return len(self.probes) + len(self.artifacts) + len(self.docs)


@dataclass(frozen=True)
class Ledger:
    """Parsed sediment ledger with the digest of the file it came from."""

    entries: tuple[SedimentEntry, ...]
    sha256: str
    path: str = ""
    version: str = SCHEMA

    def by_id(self, entry_id: str) -> SedimentEntry:
        for entry in self.entries:
            if entry.id == entry_id:
                return entry
        raise KeyError(entry_id)

    def ids(self) -> tuple[str, ...]:
        return tuple(entry.id for entry in self.entries)


@dataclass
class Row:
    """Outcome of one check inside an entry."""

    check: str
    status: str
    detail: str

    def as_dict(self) -> dict[str, Any]:
        return {"check": self.check, "status": self.status, "detail": self.detail}


@dataclass
class GuardReport:
    """Per-entry verdicts; ``ok`` is true only when every check passed."""

    rows: list[dict[str, Any]] = field(default_factory=list)
    ledger_sha256: str = ""

    @property
    def ok(self) -> bool:
        return all(row["status"] == PASS for row in self.rows)

    def counts(self, status: str | None = None) -> int:
        if status is None:
            return len(self.rows)
        return sum(1 for row in self.rows if row["status"] == status)

    def failed(self) -> list[dict[str, Any]]:
        return [row for row in self.rows if row["status"] != PASS]

    def as_dict(self) -> dict[str, Any]:
        return {
            "schema": "tensorcad.sediment.guard/1.0",
            "ledger_sha256": self.ledger_sha256,
            "entries": len(self.rows),
            "checks": sum(len(row["checks"]) for row in self.rows),
            "counts": {
                status: self.counts(status) for status in (PASS, DRIFT, MISMATCH, MISSING, BROKEN)
            },
            "ok": self.ok,
            "rows": self.rows,
        }


# --------------------------------------------------------------------------- #
# loading / validation
# --------------------------------------------------------------------------- #
def _parse_probe(raw: Mapping[str, Any], entry_id: str) -> Probe:
    _require_keys(raw, ("call", "expected"), "probe", entry_id)
    call = raw["call"]
    if not isinstance(call, str) or ":" not in call:
        raise KnowledgeError(f"{entry_id}: probe call must look like 'module:function'")
    module, _, function = call.partition(":")
    if not module.startswith(PROBE_MODULE_PREFIX) or not function or function.startswith("_"):
        raise KnowledgeError(
            f"{entry_id}: probe {call!r} must be a public callable under {PROBE_MODULE_PREFIX!r}"
        )
    args = tuple(raw.get("args", ()))
    kwargs = dict(_require_mapping(raw.get("kwargs", {}), "probe kwargs", entry_id))
    rel_tol = float(raw.get("rel_tol", 1e-12))
    abs_tol = float(raw.get("abs_tol", 0.0))
    if rel_tol < 0 or abs_tol < 0:
        raise KnowledgeError(f"{entry_id}: tolerances must be non-negative")
    index = raw.get("index")
    if index is not None and not isinstance(index, int):
        raise KnowledgeError(f"{entry_id}: probe index must be an integer")
    attr = raw.get("attr")
    if attr is not None:
        if not isinstance(attr, str) or not attr.isidentifier() or attr.startswith("_"):
            raise KnowledgeError(f"{entry_id}: probe attr must be a public attribute name")
    return Probe(
        call=call,
        expected=raw["expected"],
        args=args,
        kwargs=kwargs,
        rel_tol=rel_tol,
        abs_tol=abs_tol,
        index=index,
        attr=attr,
        note=str(raw.get("note", "")),
    )


def _parse_artifact(raw: Mapping[str, Any], entry_id: str) -> ArtifactRef:
    _require_keys(raw, ("path", "json_pointer", "expected", "sha256"), "artifact", entry_id)
    digest = str(raw["sha256"]).strip().lower()
    if len(digest) != 64 or any(ch not in "0123456789abcdef" for ch in digest):
        raise KnowledgeError(f"{entry_id}: artifact sha256 must be 64 hex characters")
    return ArtifactRef(
        path=str(raw["path"]),
        json_pointer=str(raw["json_pointer"]),
        expected=raw["expected"],
        sha256=digest,
        rel_tol=float(raw.get("rel_tol", 1e-12)),
        abs_tol=float(raw.get("abs_tol", 0.0)),
        note=str(raw.get("note", "")),
    )


def _parse_doc(raw: Mapping[str, Any], entry_id: str) -> DocRef:
    _require_keys(raw, ("path", "contains"), "doc", entry_id)
    phrase = str(raw["contains"])
    if len(phrase.strip()) < 4:
        raise KnowledgeError(f"{entry_id}: doc phrase {phrase!r} is too short to be evidence")
    return DocRef(path=str(raw["path"]), contains=phrase, note=str(raw.get("note", "")))


def _parse_entry(raw: Mapping[str, Any]) -> SedimentEntry:
    _require_keys(raw, ("id", "title", "claim", "kind", "acceptance"), "entry", "<unknown>")
    entry_id = str(raw["id"])
    if not entry_id or " " in entry_id:
        raise KnowledgeError(f"{entry_id!r}: entry id must be a non-empty slug without spaces")
    kind = str(raw["kind"])
    if kind not in KINDS:
        raise KnowledgeError(f"{entry_id}: kind {kind!r} not in {KINDS}")
    probes = tuple(_parse_probe(p, entry_id) for p in raw.get("probes", ()))
    artifacts = tuple(_parse_artifact(a, entry_id) for a in raw.get("artifacts", ()))
    docs = tuple(_parse_doc(d, entry_id) for d in raw.get("docs", ()))
    if not (probes or artifacts or docs):
        raise KnowledgeError(f"{entry_id}: an entry with no probe/artifact/doc is not checkable")
    return SedimentEntry(
        id=entry_id,
        title=str(raw["title"]),
        claim=str(raw["claim"]),
        kind=kind,
        acceptance=str(raw["acceptance"]),
        probes=probes,
        artifacts=artifacts,
        docs=docs,
        tags=tuple(str(tag) for tag in raw.get("tags", ())),
    )


def parse_ledger(document: Mapping[str, Any], *, digest: str = "", path: str = "") -> Ledger:
    """Validate a ledger document in memory (no filesystem access)."""
    _require_keys(document, ("schema", "entries"), "ledger", "<ledger>")
    if document["schema"] != SCHEMA:
        raise KnowledgeError(f"ledger schema {document['schema']!r} != {SCHEMA!r}")
    raw_entries = document["entries"]
    if not isinstance(raw_entries, list) or not raw_entries:
        raise KnowledgeError("ledger.entries must be a non-empty list")
    entries = tuple(_parse_entry(raw) for raw in raw_entries)
    seen: set[str] = set()
    for entry in entries:
        if entry.id in seen:
            raise KnowledgeError(f"duplicate entry id {entry.id!r}")
        seen.add(entry.id)
    return Ledger(entries=entries, sha256=digest, path=path, version=SCHEMA)


def load_ledger(path: str | Path) -> Ledger:
    """Read + validate the ledger file, recording its digest for provenance."""
    path = Path(path)
    raw = path.read_bytes()
    document = json.loads(raw.decode("utf-8"))
    return parse_ledger(
        document,
        digest=hashlib.sha256(raw).hexdigest(),
        path=str(path),
    )


# --------------------------------------------------------------------------- #
# verification
# --------------------------------------------------------------------------- #
def _resolve_callable(call: str):
    module_name, _, function_name = call.partition(":")
    if not module_name.startswith(PROBE_MODULE_PREFIX) or function_name.startswith("_"):
        raise KnowledgeError(f"probe {call!r} is outside the admissible namespace")
    try:
        module = importlib.import_module(module_name)
    except ImportError as exc:  # pragma: no cover - only on a broken install
        raise KnowledgeError(f"probe {call!r}: cannot import {module_name!r}: {exc}") from exc
    try:
        function = getattr(module, function_name)
    except AttributeError as exc:
        raise KnowledgeError(f"probe {call!r}: {module_name!r} has no {function_name!r}") from exc
    if not callable(function):
        raise KnowledgeError(f"probe {call!r}: attribute is not callable")
    return function


def run_probe(probe: Probe) -> Any:
    """Evaluate one probe and return the observed value."""
    function = _resolve_callable(probe.call)
    value = function(*probe.args, **dict(probe.kwargs))
    if probe.index is not None:
        try:
            value = value[probe.index]
        except (IndexError, KeyError, TypeError) as exc:
            raise KnowledgeError(
                f"probe {probe.call!r}: cannot take index {probe.index}: {exc}"
            ) from exc
    if probe.attr is not None:
        try:
            value = getattr(value, probe.attr)
        except AttributeError as exc:
            raise KnowledgeError(
                f"probe {probe.call!r}: result has no attribute {probe.attr!r}"
            ) from exc
    return value


def _check_probe(probe: Probe) -> Row:
    suffix = ""
    if probe.index is not None:
        suffix += f"[{probe.index}]"
    if probe.attr is not None:
        suffix += f".{probe.attr}"
    label = f"probe:{probe.call}{suffix}"
    try:
        observed = run_probe(probe)
    except KnowledgeError as exc:
        return Row(label, BROKEN, str(exc))
    except Exception as exc:  # noqa: BLE001 - a raising probe is a broken claim, not a skip
        return Row(label, BROKEN, f"{type(exc).__name__}: {exc}")
    if _within(observed, probe.expected, probe.rel_tol, probe.abs_tol):
        return Row(label, PASS, f"observed {observed!r} == expected {probe.expected!r}")
    return Row(
        label,
        MISMATCH,
        f"observed {observed!r} != expected {probe.expected!r} "
        f"(rel_tol={probe.rel_tol:g}, abs_tol={probe.abs_tol:g})",
    )


def _check_artifact(ref: ArtifactRef, root: Path) -> list[Row]:
    path = root / ref.path
    rows: list[Row] = []
    if not path.is_file():
        return [Row(f"artifact:{ref.path}", MISSING, "file does not exist")]
    digest = sha256_file(path)
    if digest != ref.sha256:
        rows.append(
            Row(
                f"artifact-digest:{ref.path}",
                DRIFT,
                f"sha256 {digest} != pinned {ref.sha256} "
                "(artifact regenerated: re-verify the value, then update the ledger)",
            )
        )
        # Keep going: a changed artifact should also report whether the value still holds.
    try:
        document = json.loads(path.read_text(encoding="utf-8"))
        observed = resolve_pointer(document, ref.json_pointer)
    except (KnowledgeError, ValueError) as exc:
        rows.append(Row(f"artifact:{ref.path}#{ref.json_pointer}", BROKEN, str(exc)))
        return rows
    if _within(observed, ref.expected, ref.rel_tol, ref.abs_tol):
        rows.append(
            Row(
                f"artifact:{ref.path}#{ref.json_pointer}",
                PASS,
                f"observed {observed!r} == expected {ref.expected!r}",
            )
        )
    else:
        rows.append(
            Row(
                f"artifact:{ref.path}#{ref.json_pointer}",
                MISMATCH,
                f"observed {observed!r} != expected {ref.expected!r} "
                f"(rel_tol={ref.rel_tol:g}, abs_tol={ref.abs_tol:g})",
            )
        )
    return rows


def _check_doc(ref: DocRef, root: Path) -> Row:
    path = root / ref.path
    label = f"doc:{ref.path}"
    if not path.is_file():
        return Row(label, MISSING, "file does not exist")
    text = path.read_text(encoding="utf-8")
    if ref.contains in text:
        return Row(label, PASS, f"contains {ref.contains!r}")
    return Row(label, MISSING, f"phrase {ref.contains!r} no longer appears")


def verify_entry(entry: SedimentEntry, root: str | Path = ".") -> dict[str, Any]:
    """Run every check of one entry; never short-circuits."""
    root_path = Path(root)
    rows: list[Row] = []
    for probe in entry.probes:
        rows.append(_check_probe(probe))
    for ref in entry.artifacts:
        rows.extend(_check_artifact(ref, root_path))
    for doc in entry.docs:
        rows.append(_check_doc(doc, root_path))
    status = PASS
    for row in rows:
        if row.status != PASS:
            status = row.status
            break
    return {
        "id": entry.id,
        "kind": entry.kind,
        "status": status,
        "checks": [row.as_dict() for row in rows],
    }


def verify(ledger: Ledger, root: str | Path = ".") -> GuardReport:
    """Verify the whole ledger against the working tree at ``root``."""
    report = GuardReport(ledger_sha256=ledger.sha256)
    report.rows = [verify_entry(entry, root) for entry in ledger.entries]
    return report


# --------------------------------------------------------------------------- #
# active intake: harvest un-pinned numeric leaves
# --------------------------------------------------------------------------- #
#: Artifact globs the active intake walks, relative to the repository root.
HARVEST_GLOBS: tuple[str, ...] = (
    "benchmarks/verified/**/result.json",
    "benchmarks/verified/**/gates_verified.json",
    "validation_report.json",
    "artifacts/knowledge/*/manifest.json",
)


def _walk_numbers(node: Any, pointer: str = "", out: list[tuple[str, float]] | None = None):
    out = [] if out is None else out
    if isinstance(node, Mapping):
        for key in node:
            _walk_numbers(node[key], f"{pointer}/{key}", out)
    elif isinstance(node, list):
        for index, item in enumerate(node):
            _walk_numbers(item, f"{pointer}/{index}", out)
    elif isinstance(node, bool):
        pass
    elif isinstance(node, (int, float)):
        out.append((pointer, float(node)))
    return out


def harvest_candidates(
    root: str | Path = ".",
    ledger: Ledger | None = None,
    *,
    limit: int = 40,
) -> list[dict[str, Any]]:
    """List numeric leaves in benchmark artifacts that no ledger entry pins yet.

    This is the *active* intake mouth: it turns "we produced a number somewhere"
    into a queue of claims a reviewer can promote into the ledger.
    """
    root_path = Path(root)
    pinned: set[tuple[str, str]] = set()
    if ledger is not None:
        for entry in ledger.entries:
            for ref in entry.artifacts:
                pinned.add((ref.path, ref.json_pointer))
    candidates: list[dict[str, Any]] = []
    for pattern in HARVEST_GLOBS:
        for path in sorted(root_path.glob(pattern)):
            rel = path.relative_to(root_path).as_posix()
            try:
                document = json.loads(path.read_text(encoding="utf-8"))
            except ValueError:
                continue
            for pointer, value in _walk_numbers(document):
                if (rel, pointer) in pinned:
                    continue
                candidates.append(
                    {
                        "path": rel,
                        "json_pointer": pointer,
                        "value": value,
                        "sha256": sha256_file(path),
                    }
                )
    candidates.sort(key=lambda item: (item["path"], item["json_pointer"]))
    return candidates[:limit]


# --------------------------------------------------------------------------- #
# outlet 1: agent retrieval index
# --------------------------------------------------------------------------- #
def agent_index_records(ledger: Ledger) -> list[dict[str, Any]]:
    """One record per claim: what is verified, how tightly, how to re-derive it."""
    records = []
    for entry in sorted(ledger.entries, key=lambda item: item.id):
        records.append(
            {
                "id": entry.id,
                "title": entry.title,
                "claim": entry.claim,
                "kind": entry.kind,
                "acceptance": entry.acceptance,
                "tags": list(entry.tags),
                "reproduce": {
                    "probes": [
                        {
                            "call": probe.call,
                            "args": list(probe.args),
                            "kwargs": dict(probe.kwargs),
                            "expected": probe.expected,
                            "rel_tol": probe.rel_tol,
                            "abs_tol": probe.abs_tol,
                            **({"index": probe.index} if probe.index is not None else {}),
                            **({"attr": probe.attr} if probe.attr is not None else {}),
                        }
                        for probe in entry.probes
                    ],
                    "artifacts": [
                        {
                            "path": ref.path,
                            "json_pointer": ref.json_pointer,
                            "expected": ref.expected,
                            "sha256": ref.sha256,
                        }
                        for ref in entry.artifacts
                    ],
                    "docs": [{"path": doc.path, "contains": doc.contains} for doc in entry.docs],
                },
            }
        )
    return records


def _write_jsonl(path: Path, records: Iterable[Mapping[str, Any]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    payload = []
    for record in records:
        payload.append(json.dumps(record, sort_keys=True, ensure_ascii=False))
        count += 1
    path.write_text("\n".join(payload) + "\n" if payload else "", encoding="utf-8", newline="\n")
    return count


def export_agent_index(ledger: Ledger, path: str | Path) -> int:
    """Write the agent-facing retrieval index (JSONL, deterministic order)."""
    return _write_jsonl(Path(path), agent_index_records(ledger))


# --------------------------------------------------------------------------- #
# outlet 2: training corpus + held-out eval split
# --------------------------------------------------------------------------- #
#: Deterministic geometry sweep for the closed-form corpus (all combinations).
CORPUS_LENGTHS: tuple[float, ...] = (6.0, 8.0, 10.0, 12.0)
CORPUS_BEAMS: tuple[float, ...] = (0.8, 1.0, 1.2)
CORPUS_DRAFT_RATIO: float = 0.625
CORPUS_LONGITUDINAL: tuple[float, ...] = (1.5, 2.0, 2.5)
CORPUS_VERTICAL: tuple[float, ...] = (1.5, 2.0)
#: Held-out fraction expressed as "one hash residue in N" (deterministic, not random).
EVAL_MODULUS = 5
#: Exact closed forms are compared at this relative tolerance by consumers.
CORPUS_REL_TOL = 1e-9


def _corpus_targets(length: float, beam: float, draft: float, n: float, m: float) -> dict[str, dict[str, Any]]:
    """Exact Wigley hydrostatics for one geometry, straight from the closed forms."""
    from . import wigley

    values = {
        "displacement_volume_m3": wigley.analytic_volume(length, beam, draft, n, m),
        "waterplane_area_m2": wigley.analytic_waterplane_area(length, beam, n),
        "vertical_centre_buoyancy_m": wigley.analytic_vertical_centre_buoyancy(draft, m),
        "block_coefficient": wigley.analytic_block_coefficient(n, m),
        "prismatic_coefficient": wigley.analytic_prismatic_coefficient(n, m),
        "midship_coefficient": wigley.analytic_midship_coefficient(m),
        "waterplane_coefficient": wigley.analytic_waterplane_coefficient(n),
    }
    return {name: {"value": value, "rel_tol": CORPUS_REL_TOL} for name, value in values.items()}


def corpus_records(ledger: Ledger | None = None) -> list[dict[str, Any]]:
    """Build the deterministic corpus (geometry -> exact hydrostatics)."""
    ledger_sha = ledger.sha256 if ledger is not None else ""
    records = []
    for length in CORPUS_LENGTHS:
        for beam in CORPUS_BEAMS:
            draft = round(beam * CORPUS_DRAFT_RATIO, 12)
            for n in CORPUS_LONGITUDINAL:
                for m in CORPUS_VERTICAL:
                    record_id = (
                        f"wigley-L{length:g}-B{beam:g}-T{draft:g}-n{n:g}-m{m:g}"
                    )
                    records.append(
                        {
                            "id": record_id,
                            "split": _split_of(record_id),
                            "domain": "ship_hydrostatics",
                            "geometry_family": "wigley",
                            "inputs": {
                                "length_m": length,
                                "beam_m": beam,
                                "draft_m": draft,
                                "longitudinal_exponent": n,
                                "vertical_exponent": m,
                            },
                            "targets": _corpus_targets(length, beam, draft, n, m),
                            "provenance": {
                                "generator": "tensorcad.knowledge/export_training_corpus/1.0",
                                "source": "closed form: tensorcad.wigley analytic_* functions",
                                "reference": "Wigley parabolic hull form, n=m=2 is the classic case",
                                "ledger_sha256": ledger_sha,
                            },
                        }
                    )
    records.sort(key=lambda item: item["id"])
    return records


def _split_of(record_id: str) -> str:
    """Deterministic 80/20 assignment from the record id (no RNG, no clock)."""
    residue = int(hashlib.sha256(record_id.encode("utf-8")).hexdigest()[:8], 16) % EVAL_MODULUS
    return "eval" if residue == 0 else "train"


def export_training_corpus(
    outdir: str | Path,
    ledger: Ledger | None = None,
    *,
    write_agent_index: bool = True,
) -> dict[str, Any]:
    """Write corpus/eval/manifest (and the agent index) and return the manifest."""
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    records = corpus_records(ledger)
    train = [record for record in records if record["split"] == "train"]
    held_out = [record for record in records if record["split"] == "eval"]
    files: dict[str, Any] = {}

    for name, subset in (("training_corpus.jsonl", train), ("eval_set.jsonl", held_out)):
        count = _write_jsonl(out / name, subset)
        files[name] = {
            "records": count,
            "sha256": sha256_file(out / name),
            "bytes": (out / name).stat().st_size,
        }
    if write_agent_index and ledger is not None:
        count = export_agent_index(ledger, out / "agent_index.jsonl")
        files["agent_index.jsonl"] = {
            "records": count,
            "sha256": sha256_file(out / "agent_index.jsonl"),
            "bytes": (out / "agent_index.jsonl").stat().st_size,
        }

    manifest = {
        "schema": EXPORT_SCHEMA,
        "generator": "tensorcad.knowledge/export_training_corpus/1.0",
        "deterministic": True,
        "counts": {
            "records_total": len(records),
            "train": len(train),
            "eval": len(held_out),
            "provenance": {
                "family": "wigley",
                "exact": True,
                "rel_tol": CORPUS_REL_TOL,
                "targets": sorted(records[0]["targets"]) if records else [],
                "sweep": {
                    "length_m": list(CORPUS_LENGTHS),
                    "beam_m": list(CORPUS_BEAMS),
                    "draft_ratio": CORPUS_DRAFT_RATIO,
                    "longitudinal_exponent": list(CORPUS_LONGITUDINAL),
                    "vertical_exponent": list(CORPUS_VERTICAL),
                },
            },
        },
        "ledger": (
            {"path": ledger.path, "sha256": ledger.sha256, "entries": len(ledger.entries)}
            if ledger is not None
            else None
        ),
        "files": files,
    }
    (out / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n"
    )
    return manifest


def read_jsonl(path: str | Path) -> list[dict[str, Any]]:
    """Read a JSONL export back (used by tests and downstream consumers)."""
    return [
        json.loads(line)
        for line in Path(path).read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]