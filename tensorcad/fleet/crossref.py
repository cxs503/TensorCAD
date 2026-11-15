"""Cross-source reference library for the six benchmark hulls.

Anchor ``REPRO-GEOM-FLEET-REFS``.  The fleet benchmark
(:mod:`tensorcad.fleet.official`, ``benchmarks/verified/baseline_fleet``) compares a
reconstructed displacement against **one** published number per hull.  A single number
cannot tell a reconstruction error from a source typo, so this module keeps a *library* of
published principal particulars -- displacement volume ∇, block coefficient Cb, LCB (% Lpp,
forward positive) and length on the waterline Lwl, plus the Lpp/Bwl/T the coefficients are
defined on -- and asks, per quantity, whether at least two **independent publishers** state
it and whether they agree.

Design rules that the tests enforce
-----------------------------------
* Nothing is transcribed by hand.  :func:`parse_source` reads a *fetched byte string* whose
  SHA-256 is pinned in :data:`SOURCES`; every cell carries the raw cell text, the row it
  came from and the source digest, so a value can always be traced back to the bytes.
* Independence is a property of the *publisher family*, not of the URL.  Two pages served
  by the same workshop committee are one family; the same workshop series in two years is
  still one family unless the pages are demonstrably different editions (recorded in
  ``family``).  A quantity backed by one family is ``SINGLE_SOURCE`` even when three URLs
  repeat it -- that is the failure mode this anchor exists to expose.
* A disagreement between families is reported as ``CONFLICT`` with per-family values; it is
  never averaged away.
"""

from __future__ import annotations

import re
from dataclasses import asdict, dataclass, field
from typing import Any, Iterable, Mapping, Sequence

from .official import table_rows

#: Quantities the library tracks.  ``displacement_volume`` is always a volume in m^3;
#: a publisher that states a mass is converted at :data:`SEAWATER_DENSITY` and the raw
#: mass is kept in ``Cell.raw`` together with ``Cell.unit_basis == "mass"``.
QUANTITIES: tuple[str, ...] = (
    "displacement_volume",
    "cb",
    "lcb_percent",
    "lwl",
    "length_pp",
    "beam",
    "draft",
)

#: Human labels used in the report, in report order.
QUANTITY_LABELS: dict[str, str] = {
    "displacement_volume": "∇ (m³)",
    "cb": "Cb",
    "lcb_percent": "LCB (% Lpp, fwd+)",
    "lwl": "Lwl (m)",
    "length_pp": "Lpp (m)",
    "beam": "Bwl (m)",
    "draft": "T (m)",
}

#: Seawater density used only when a publisher states a *mass* (t) instead of a volume.
SEAWATER_DENSITY = 1.025

#: Agreement tolerance per quantity.  Published reproductions of the same principal
#: particular differ in the last printed digit; 0.5 % covers rounding, not a different ship.
TOLERANCE: dict[str, float] = {
    "displacement_volume": 0.005,
    "cb": 0.005,
    "lcb_percent": 0.02,  # percent-of-Lpp is a small number; 2 % of it is 0.03 %Lpp
    "lwl": 0.005,
    "length_pp": 0.005,
    "beam": 0.005,
    "draft": 0.005,
}

_NUMBER = re.compile(r"[-+]?\d+(?:\.\d+)?")

#: Row-label patterns.  Matched case-insensitively against the *whole* flattened row so a
#: page that prints "Maximum beam of waterline B WL (m)" and one that prints "Bwl (m)" are
#: both understood.  Order matters: the first pattern that matches wins, so LCB is tested
#: before Cb and Lwl before Lpp.
LABEL_PATTERNS: tuple[tuple[str, str], ...] = (
    ("lcb_percent", r"\blcb\b|longitudinal cent(?:er|re) of buoyancy"),
    ("lwl", r"\bl\s?wl\b|waterline length|length (?:of|on) (?:the )?water\s?line|length at (?:the )?waterline"),
    ("length_pp", r"\bl\s?pp\b|length between perpendiculars"),
    ("beam", r"\bb\s?wl\b|maximum beam|breadth"),
    ("draft", r"\bdraft\b|\bdraught\b|^\s*t\s*(?:\(|\b)"),
    ("displacement_volume", r"displacement"),
    ("cb", r"\bc\s?b\b|block coefficient"),
)


@dataclass(frozen=True)
class SourceSpec:
    """One fetched publisher document and how to read it."""

    key: str
    ship: str
    publisher: str
    #: Independence family.  Two sources from the same family never count as two sources.
    family: str
    url: str
    kind: str  # "html" | "pdf"
    #: SHA-256 of the fetched bytes (pinned after the first fetch, verified on every re-fetch).
    sha256: str
    licence: str
    note: str = ""
    #: Column holding the *full-scale* value; ``None`` = first numeric cell after the label.
    column: int | None = None
    #: For multi-ship pages: the sub-table heading that starts this ship's block.
    section: str | None = None
    #: Origin of the numerals *as the document itself declares it*; ``None`` = the document's
    #: own publisher.  Two documents that reproduce one upstream table are one lineage.
    lineage: str | None = None
    #: For ``kind="pdf"`` documents: ``(quantity, regex, basis)`` triples.  ``regex`` must have
    #: exactly one capture group around the numeral; ``basis`` is ``"volume" | "mass"``.
    quotes: tuple[tuple[str, str, str], ...] = ()


#: The document set.  Digests are filled in by ``tools/baseline_fleet/fetch_crossref.py --pin``.
SOURCES: tuple[SourceSpec, ...] = (
    SourceSpec(
        key="t2015_kcs",
        ship="KCS",
        publisher="NMRI Tokyo 2015 workshop (G2015 geometry pages)",
        family="NMRI-T2015",
        url="https://www.t2015.nmri.go.jp/kcs_gc.html",
        kind="html",
        sha256="68d9c24600cb699b16394fc5c32762ebc875c856d14ef8460e1d0ebf4e7c78a4",
        licence="page copyright NMRI; digests + extracted cells only, page bytes not redistributed",
        column=2,
    ),
    SourceSpec(
        key="t2015_jbc",
        ship="JBC",
        publisher="NMRI Tokyo 2015 workshop (G2015 geometry pages)",
        family="NMRI-T2015",
        url="https://www.t2015.nmri.go.jp/jbc_gc.html",
        kind="html",
        sha256="fc2ce9054a8b7842162bac94365853fa076da2527e8c7e455832f9710f58f875",
        licence="page copyright NMRI; digests + extracted cells only",
        column=2,
    ),
    SourceSpec(
        key="t2015_onrt",
        ship="ONRT",
        publisher="NMRI Tokyo 2015 workshop (G2015 geometry pages)",
        family="NMRI-T2015",
        url="https://www.t2015.nmri.go.jp/onrt_gc.html",
        kind="html",
        sha256="5e13427332de8be17b9dc2245dc5830300baf2e2340efcd8a80eee2269c43a11",
        licence="page copyright NMRI; digests + extracted cells only",
        column=3,
    ),
    SourceSpec(
        key="simman2020_summary",
        ship="KCS|KVLCC2|ONRT",
        publisher="SIMMAN 2020 workshop (KRISO/SNU) — Table 1/2/3 main particulars",
        family="SIMMAN-2020",
        url="https://simman2020.kr/contents/Ship_data_summary.php",
        kind="html",
        sha256="c1a668f6d2440c5d83a83c3616285e9288ef618dbc9a2feebd26627be3f4ebb3",
        licence="page copyright SIMMAN2020; digests + extracted cells only",
    ),
    SourceSpec(
        key="simman2014_kvlcc2",
        ship="KVLCC2",
        publisher="SIMMAN 2014 workshop (FORCE Technology) — KVLCC2 geometry page",
        family="SIMMAN-2014",
        url="https://simman2014.dk/ship-data/moeri-kvlcc2-tanker/geometry-and-conditions-moeri-kvlcc2-tanker/",
        kind="html",
        sha256="bfc9101e3a4afa4fbe02750de3f0ec830ccc948074265ad4e63184f5b67d310b",
        licence="page copyright FORCE Technology; digests + extracted cells only",
    ),
    SourceSpec(
        key="simman2014_kcs",
        ship="KCS",
        publisher="SIMMAN 2014 workshop (FORCE Technology) — MOERI container ship geometry page",
        family="SIMMAN-2014",
        url="https://simman2014.dk/ship-data/moeri-container-ship/geometry-and-conditions-moeri-container-ship/",
        kind="html",
        sha256="3cbc3a22daf3d494ffb26ce81a1bc951cd5dca6cb04c27dd2fd478f9a5dc9bf1",
        licence="page copyright FORCE Technology; digests + extracted cells only",
    ),
    SourceSpec(
        key="ittc_benchmarks",
        ship="KCS|KVLCC2|JBC|ONRT|Series60|DTC",
        publisher="ITTC Benchmark Repository — 'List of Benchmarks' (ITTC, Zürich)",
        family="ITTC",
        url="https://ittc.info/media/11250/list-of-benchmarks-2.pdf",
        kind="pdf",
        sha256="da101e9a352108e4cc73eb5f871fbe27d585f8e7bd4fb21af4261edd85fbaca1",
        licence="ITTC publication; bibliographic endorsement of the hull set (no numerals read)",
    ),
)


#: Second-round documents: open-access full texts obtained from third-party repositories when
#: the publishers' own pages were unreachable (see ``PROBES``).  Each carries the digit-level
#: quote it was read from, the digest of the fetched bytes, and the *lineage* it declares when
#: the text says it is reproducing another body's table.
_OA_SOURCES: tuple[SourceSpec, ...] = (
    SourceSpec(
        key="kvlcc2_ugent",
        ship="KVLCC2",
        publisher="Journal of Ship Research 61(3) 2017 (Ghent University) - Table 1",
        family="JSR-2017",
        lineage="JSR-2017",
        url="https://biblio.ugent.be/publication/8563318/file/8563319.pdf",
        kind="pdf",
        sha256="172f78d49f36be7455c20295fc51de7ba70ecf0293cc96455e59e0e0bb8f0568",
        licence="INSTITUTIONAL REPOSITORY; digests + quoted cells only, bytes not redistributed",
        note="Table 1 columns: KVLCC2 | DTC heavy ballast | DTC design draught | DTC light ballast; displacement printed in tonnes",
        quotes=(
            ("length_pp", r"Lpp\s*\(\s*m\s*\)\s*([0-9.]+)\s*[0-9.]+", "volume"),
            ("beam", r"\bB\s*\(\s*m\s*\)\s*([0-9.]+)\s*[0-9.]+", "volume"),
            ("draft", r"TFP\s*/\s*TAP\s*\(\s*m\s*\)\s*([0-9.]+)\s*/", "volume"),
            ("displacement_volume", r"\(t\)\s*([0-9,]{6,})", "mass"),
            ("cb", r"CB\s*\[\s*-\s*\]\s*([0-9.]+)\s*[0-9.]+", "volume"),
        ),
    ),
    SourceSpec(
        key="dtc_ugent",
        ship="DTC",
        publisher="Journal of Ship Research 61(3) 2017 (Ghent University) - Table 1 (design-draught column)",
        family="JSR-2017",
        lineage="JSR-2017",
        url="https://biblio.ugent.be/publication/8563318/file/8563319.pdf",
        kind="pdf",
        sha256="172f78d49f36be7455c20295fc51de7ba70ecf0293cc96455e59e0e0bb8f0568",
        licence="INSTITUTIONAL REPOSITORY; digests + quoted cells only",
        note="the column is labelled '(t)' but the numeral 173468 equals the volume in m3 (Cb closure 173467.4), so the basis is volume; a mass reading would be 2.4 % off",
        quotes=(
            ("length_pp", r"Lpp\s*\(\s*m\s*\)\s*[0-9.]+\s*([0-9.]+)", "volume"),
            ("beam", r"\bB\s*\(\s*m\s*\)\s*[0-9.]+\s*([0-9.]+)", "volume"),
            ("draft", r"TFP\s*/\s*TAP\s*\(\s*m\s*\)\s*[0-9.]+/[0-9.]+\s*[0-9.]+/[0-9.]+\s*([0-9.]+)/", "volume"),
            ("displacement_volume", r"\(t\)\s*[0-9,]+\s*[0-9,]+\s*([0-9,]{6,})", "volume"),
            ("cb", r"CB\s*\[\s*-\s*\]\s*[0-9.]+\s*([0-9.]+)", "volume"),
        ),
    ),
    SourceSpec(
        key="kvlcc2_hrcak",
        ship="KVLCC2",
        publisher="Brodogradnja 68(4) 2017, Table 1 'Principal dimensions of KVLCC2 (SIMMAN2008)'",
        family="BRODOGRADNJA-2017",
        lineage="SIMMAN-SERIES",
        url="https://hrcak.srce.hr/file/274163",
        kind="pdf",
        sha256="c20acf9e8eae047a0cb4070bd55aa48d8eaddf3a92d27f2ded6adf4499af5f16",
        licence="OA journal (CC BY-NC); digests + quoted cells only",
        note="the table caption names SIMMAN 2008, so this document reproduces the SIMMAN table and never counts as an independent origin",
        quotes=(
            ("length_pp", r"pp\s*L\s*([0-9.]+) m", "volume"),
            ("beam", r"\nB\s*\n([0-9.]+) m", "volume"),
            ("draft", r"\nT\s*\n([0-9.]+) m", "volume"),
            ("displacement_volume", r"\n(312622)", "volume"),
            ("cb", r"\nB\s*\nC\s*\n([0-9.]+)", "volume"),
        ),
    ),
    SourceSpec(
        key="kvlcc2_jstage",
        ship="KVLCC2",
        publisher="JASNAOE 23 (2016) - RANS KVLCC2 with body-force propeller, Table 1",
        family="JASNAOE-2016",
        lineage="JASNAOE-2016",
        url="https://www.jstage.jst.go.jp/article/jjasnaoe/23/0/23_1/_pdf",
        kind="pdf",
        sha256="c97b1d3e62b27e19e8a0bbbd9593411b9250b1141d7cb789c60e64dac8d2a5aa",
        licence="J-STAGE OA; digests + quoted cells only",
        note="states Lpp and draft only (propeller-focused paper)",
        quotes=(
            ("length_pp", r"Lpp\s*([0-9.]+) m", "volume"),
            ("draft", r"Draft\s*\n?([0-9.]+) m", "volume"),
        ),
    ),
    SourceSpec(
        key="jbc_cfdletters",
        ship="JBC",
        publisher="CFD Letters 15(4) 2023 - Table 1 General Parameters of the JBC at full scale",
        family="CFDLETTERS-2023",
        lineage="NMRI-T2015",
        url="https://semarakilmu.com.my/journals/index.php/CFD_Letters/article/download/897/940",
        kind="pdf",
        sha256="dd79da84862148d8115a221d638dcb9ab1fa9e10ba48ab3c02ae41e168dfec2a",
        licence="OA journal; digests + quoted cells only",
        note="workshop reproduction: the paper states the case comes from the Tokyo 2015 workshop",
        quotes=(
            ("length_pp", r"perpendiculars\s+LPP\s*\(\s*m\s*\)\s*([0-9.]+)", "volume"),
            ("lwl", r"waterline\s+LWL\s*\(\s*m\s*\)\s*([0-9.]+)", "volume"),
            ("beam", r"waterline\s+BWL\s*\(\s*m\s*\)\s*([0-9.]+)", "volume"),
            ("draft", r"Draft\s*\nT\s*\(\s*m\s*\)\s*([0-9.]+)", "volume"),
            ("displacement_volume", r"volume\s*.{1,4}\s*\(\s*m3\s*\)\s*([0-9.]+)", "volume"),
            ("cb", r"coefficient\s*\nCb\s*\n([0-9.]+)", "volume"),
        ),
    ),
    SourceSpec(
        key="jbc_scse",
        ship="JBC",
        publisher="Sustainable Chemical Engineering 3 (2025) - Table 2 Principal Particulars of JBC Hull",
        family="SCSE-2025",
        lineage="NMRI-T2015",
        url="https://scienpg.com/scs/index.php/engineering/article/download/scse.2025.3.22/pdf",
        kind="pdf",
        sha256="1b7503d516d7d995eca03bd1470bde8dc1f3a24e6c4697256be2133bb7b2d5fa",
        licence="OA journal; digests + quoted cells only",
        note="workshop reproduction carrying LCB and Cm columns",
        quotes=(
            ("length_pp", r"perpendicular\s+\S+\s+([0-9.]+)\s+[0-9.]+", "volume"),
            ("lwl", r"Length\s+waterline\s+\S+\s+([0-9.]+)\s+[0-9.]+", "volume"),
            ("beam", r"Beam\s+\S+\s+([0-9.]+)\s+[0-9.]+", "volume"),
            ("draft", r"Draft\s+\S+\s+([0-9.]+)\s+[0-9.]+", "volume"),
            ("displacement_volume", r"Displacement\s+volume\s+\S+\s+([0-9.]+)\s+[0-9.]+", "volume"),
            ("cb", r"efficient\s+\S+\s+([0-9.]+)\s+[0-9.]+", "volume"),
            ("lcb_percent", r"(?s)center of\s+buoy.{0,80}?([0-9]+\.[0-9]+)", "volume"),
        ),
    ),
    SourceSpec(
        key="jbc_memorial",
        ship="JBC",
        publisher="Memorial University of Newfoundland thesis 2025 - Table 3 JBC hull data",
        family="MEMORIAL-2025",
        lineage="NMRI-T2015",
        url="https://storage.googleapis.com/dissertations/2025/MEM_2025_KABONGE.pdf",
        kind="pdf",
        sha256="870e66a568e6dd1f56107cd7d88b1681d8eee40b184ed9acb3291734054f56ac",
        licence="INSTITUTIONAL REPOSITORY; digests + quoted cells only",
        note="the table's own source line cites NMRI (2015a) / Nobuyuki (2015)",
        quotes=(
            ("length_pp", r"3\.\s*LPP\s+([0-9.]+)", "volume"),
            ("lwl", r"4\.\s*LWL\s+([0-9.]+)", "volume"),
            ("beam", r"2\.\s*Beam,\s*B\s+([0-9.]+)", "volume"),
            ("draft", r"5\.\s*Draught\s+([0-9.]+)", "volume"),
            ("displacement_volume", r"6\.\s*Displacement\s+([0-9.]+)", "volume"),
        ),
    ),
    SourceSpec(
        key="dtc_strath",
        ship="DTC",
        publisher="ISNAM 2018 (Strathclyde repository) - Table 1 Main particulars of DTC hull",
        family="STRATH-ISNAM2018",
        lineage="DTC-PRIMARY",
        url="https://strathprints.strath.ac.uk/67912/1/Dogrul_etal_ISNAM2018_Uncertainty_assessment_and_self_propulsion_estimation_of_Duisburg_Test_Case.pdf",
        kind="pdf",
        sha256="550de4394a313431f162a71f5570cbad7c5e79869bea5e8c197f5c98a5e0cd33",
        licence="INSTITUTIONAL REPOSITORY; digests + quoted cells only",
        note="full-scale column of the DTC particulars table",
        quotes=(
            ("length_pp", r"Lpp\s*\(\s*m\s*\)\s*[0-9.]+\s*([0-9.]+)", "volume"),
            ("beam", r"Bwl\s*\(\s*m\s*\)\s*[0-9.]+\s*([0-9.]+)", "volume"),
            ("draft", r"T\s*\(\s*m\s*\)\s*[0-9.]+\s*([0-9.]+)", "volume"),
            ("displacement_volume", r"\(m3\)\s*[0-9.]+\s*([0-9.]+)", "volume"),
            ("cb", r"CB\s+[0-9.]+\s+([0-9.]+)", "volume"),
        ),
    ),
)

SOURCES = SOURCES + _OA_SOURCES


@dataclass
class Cell:
    """One quantity as stated by one source."""

    quantity: str
    value: float
    raw: str
    label: str
    source: str
    url: str
    sha256: str
    column: int
    unit_basis: str = "volume"
    note: str = ""

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)


def split_ships(spec: SourceSpec, payload: bytes | str) -> list[tuple[str, str]]:
    """Split a multi-ship page into one markup block per ship.

    The workshop summary page prints one numbered table per hull ("Table 2 Main Particulars
    of KCS"); the block a ship owns runs from its own heading to the next heading.  A
    single-ship page is returned whole.
    """
    markup = payload.decode("utf-8", "replace") if isinstance(payload, bytes) else payload
    ships = spec.ship.split("|")
    if len(ships) == 1:
        return [(ships[0], markup)]
    hits = [
        (m.start(), m.group(1).upper())
        for m in re.finditer(r"Table\s*\d+\s*Main Particulars of\s+([A-Za-z0-9]+)", markup, flags=re.I)
    ]
    out: list[tuple[str, str]] = []
    for index, (start, name) in enumerate(hits):
        end = hits[index + 1][0] if index + 1 < len(hits) else len(markup)
        out.append((name, markup[start:end]))
    return out


def _split_section(markup: str, spec: SourceSpec) -> str:
    """Restrict a page to the block that belongs to ``spec.ship``."""
    if not spec.section:
        return markup
    start = markup.find(spec.section)
    if start < 0:
        return ""
    nxt = re.search(r"<table", markup[start + len(spec.section):], flags=re.I)
    end = start + len(spec.section) + (nxt.start() if nxt else len(markup) - start)
    return markup[start:end]


def _pick(cells: Sequence[str], column: int | None) -> tuple[float, str, int, str] | None:
    """The numeric cell that carries the full-scale value, with its index and basis."""
    order: Iterable[int]
    if column is not None:
        order = [column] + [i for i in range(1, len(cells)) if i != column]
    else:
        order = range(1, len(cells))
    for index in order:
        if index >= len(cells):
            continue
        text = cells[index]
        match = _NUMBER.search(text.replace(",", ""))
        if match is None:
            continue
        basis = "volume"
        if re.search(r"(?i)\btonns?\b|\bton\b|\bt\s*$|\(t\)", text):
            basis = "mass"
        return float(match.group(0)), text, index, basis
    return None


def parse_source(spec: SourceSpec, payload: bytes | str) -> dict[str, Cell]:
    """Read the tracked quantities of ``spec.ship`` out of a fetched document."""
    if spec.kind != "html":
        return {}
    markup = payload.decode("utf-8", "replace") if isinstance(payload, bytes) else payload
    block = _split_section(markup, spec)
    out: dict[str, Cell] = {}
    for row in table_rows(block):
        if not row:
            continue
        # The *label* lives in the first cell.  Reading it from the flattened row instead
        # lets derivative rows win ("izz/Lpp 0.25" is not Lpp), which is how a page can
        # silently hand back the wrong quantity.
        label = row[0]
        if re.search(r"(?i)izz|ixx|hub|\bic\b", label):
            continue
        quantity = None
        for name, pattern in LABEL_PATTERNS:
            if re.search(pattern, label, flags=re.I):
                quantity = name
                break
        if quantity is None or quantity in out:
            continue
        picked = _pick(row, spec.column)
        if picked is None:
            continue
        value, raw, index, basis = picked
        if spec.column is not None and index != spec.column:
            # The configured full-scale column was empty, so a neighbour was read.  Its scale
            # is unknown: keep it for the audit trail, never let it vote.
            basis = "fallback_column"
        out[quantity] = Cell(
            quantity=quantity,
            value=value,
            raw=raw,
            label=row[0][:120],
            source=spec.key,
            url=spec.url,
            sha256=spec.sha256,
            column=index,
            unit_basis=basis,
        )
    if "displacement_volume" in out:
        cell = out["displacement_volume"]
        if cell.unit_basis == "mass":
            cell.value = round(cell.value / SEAWATER_DENSITY, 1)
            cell.note = f"mass {cell.raw} converted at {SEAWATER_DENSITY} t/m3"
    return out


def parse_pdf_text(spec: SourceSpec, text: str) -> dict[str, Cell]:
    """Read ``spec.quotes`` out of a PDF's extracted text (page breaks marked by ``\f``).

    A PDF has no cells to address, so the *quote* is the address: each pattern must match the
    table caption and its numeral together, and the matched text is stored verbatim as the
    cell's ``raw`` value together with the 1-based page it was found on.  A pattern that does
    not match yields no cell -- the document simply does not state that quantity.
    """
    out: dict[str, Cell] = {}
    if spec.kind != "pdf":
        return out
    pages = text.split("\f")
    for quantity, pattern, basis in spec.quotes:
        found = None
        for index, page in enumerate(pages, start=1):
            match = re.search(pattern, page)
            if match is not None:
                found = (index, match)
                break
        if found is None:
            continue
        page_no, match = found
        if match.re.groups < 1:  # a quote must capture the numeral, never accept the whole match
            continue
        raw = match.group(1)
        value = float(raw.replace(",", ""))
        note = ""
        if basis == "mass":
            value = value / SEAWATER_DENSITY
            note = f"mass {raw} t converted at {SEAWATER_DENSITY} t/m3"
        out[quantity] = Cell(
            quantity=quantity,
            value=round(value, 4),
            raw=raw,
            label=_quote_context(pages[page_no - 1], match).strip(),
            source=spec.key,
            url=spec.url,
            sha256=spec.sha256,
            column=0,
            unit_basis=basis,
            note=f"page {page_no}; {note}".strip("; "),
        )
    return out


def _quote_context(page: str, match: "re.Match[str]", width: int = 90) -> str:
    """Whitespace-collapsed text around a quote, so a reader can see the table it came from."""
    start = max(0, match.start() - width)
    end = min(len(page), match.end() + width)
    return re.sub(r"\s+", " ", page[start:end])

def spread(values: Sequence[float]) -> tuple[float, float]:
    """Absolute range and range relative to the mean (guard for a zero mean)."""
    lo, hi = min(values), max(values)
    mean = sum(values) / len(values)
    rel = abs(hi - lo) / abs(mean) if mean else float("inf")
    return hi - lo, rel


def cross_reference(extractions: Mapping[str, Mapping[str, Mapping[str, Cell]]]) -> dict[str, Any]:
    """Per ship, per quantity: the independent families that state it and whether they agree.

    ``extractions`` maps ``ship -> source_key -> quantity -> Cell``.
    """
    report: dict[str, Any] = {}
    for ship in sorted(extractions):
        per_quantity: dict[str, Any] = {}
        for quantity in QUANTITIES:
            entries = []
            for source_key, cells in sorted(extractions[ship].items()):
                cell = cells.get(quantity)
                if cell is None:
                    continue
                family = LINEAGE.get(source_key, FAMILIES.get(source_key, source_key))
                entries.append({"family": family, "cell": cell.as_dict()})
                if cell.unit_basis.startswith("fallback_column"):
                    continue
            families: dict[str, list[float]] = {}
            for entry in entries:
                # A cell read from a fallback column has an unknown scale: it stays in
                # ``entries`` for the audit trail but never votes.  The original loop put the
                # guard at the *end* of the collection loop, where it excluded nothing, so a
                # neighbouring-column reading could outvote (or silently corroborate) the real
                # full-scale value -- two tests had been failing on exactly that.
                if entry["cell"]["unit_basis"].startswith("fallback_column"):
                    continue
                families.setdefault(entry["family"], []).append(entry["cell"]["value"])
            family_values = {fam: sum(vals) / len(vals) for fam, vals in sorted(families.items())}
            if len(family_values) >= 2:
                delta, rel = spread(list(family_values.values()))
                tol = TOLERANCE[quantity]
                verdict = "CROSS_VERIFIED" if rel <= tol else "CONFLICT"
            else:
                delta, rel = 0.0, 0.0
                verdict = "SINGLE_SOURCE" if family_values else "MISSING"
            ordered = [
                {"family": fam, "value": value}
                for fam, value in sorted(family_values.items())
            ]
            per_quantity[quantity] = {
                "verdict": verdict,
                # A disagreement above the tolerance is *disputed*, not averaged: both
                # family values stay in ``ordered`` so a consumer must choose explicitly.
                "disputed": verdict == "CONFLICT",
                "published_values": ordered,
                "independent_families": len(family_values),
                "family_values": family_values,
                "range": delta,
                "range_relative": rel,
                "tolerance": TOLERANCE[quantity],
                "entries": entries,
            }
        # Same-family repeats are counted separately: they are corroboration, not independence.
        repeated = sum(
            1
            for quantity in QUANTITIES
            for src in [per_quantity[quantity]]
            if len(src["entries"]) > src["independent_families"] > 0
        )
        report[ship] = {
            "quantities": per_quantity,
            "cross_verified": sum(1 for q in per_quantity.values() if q["verdict"] == "CROSS_VERIFIED"),
            "single_source": sum(1 for q in per_quantity.values() if q["verdict"] == "SINGLE_SOURCE"),
            "conflict": sum(1 for q in per_quantity.values() if q["verdict"] == "CONFLICT"),
            "missing": sum(1 for q in per_quantity.values() if q["verdict"] == "MISSING"),
            "corroborated_by_repeat_family": repeated,
            "disputed_quantities": [
                q for q in QUANTITIES if per_quantity[q]["verdict"] == "CONFLICT"
            ],
        }
    return report


#: ``source key -> independence family``, derived from :data:`SOURCES`.
FAMILIES: dict[str, str] = {spec.key: spec.family for spec in SOURCES}

#: Independence *lineage*: the unit that may count as one independent source.  It is
#: stricter than :data:`FAMILIES`, because a later edition of the same workshop reproduces
#: the same upstream table: SIMMAN 2014 and SIMMAN 2020 state KVLCC2 as ``Bwl 580`` where the
#: ship is 58.0 m wide -- the identical defect in both editions is the evidence that one
#: table was copied from the other, so their agreement is *not* corroboration.  KCS is
#: backed by two lineages (the NMRI Tokyo 2015 page and the SIMMAN series), KVLCC2 by one.
LINEAGE: dict[str, str] = {
    "t2015_kcs": "NMRI-T2015",
    "t2015_jbc": "NMRI-T2015",
    "t2015_onrt": "NMRI-T2015",
    "simman2020_summary": "SIMMAN-SERIES",
    "simman2014_kvlcc2": "SIMMAN-SERIES",
    "simman2014_kcs": "SIMMAN-SERIES",
    "ittc_benchmarks": "ITTC",
}

LINEAGE.update(
    {
        "kvlcc2_ugent": "JSR-2017",
        "dtc_ugent": "JSR-2017",
        "kvlcc2_hrcak": "SIMMAN-SERIES",
        "kvlcc2_jstage": "JASNAOE-2016",
        "jbc_cfdletters": "NMRI-T2015",
        "jbc_scse": "NMRI-T2015",
        "jbc_memorial": "NMRI-T2015",
        "dtc_strath": "DTC-PRIMARY",
    }
)


def lineage_of(key: str) -> str:
    """Independence lineage of a source key (falls back to its publisher family)."""
    return LINEAGE.get(key, FAMILIES.get(key, key))


def closure_check(report: Mapping[str, Any], ships: Sequence[str] | None = None) -> dict[str, Any]:
    """Does each source's own ∇ match its own ``Cb·Lpp·Bwl·T``?

    A publisher table is internally consistent or it is not; the check needs no external
    truth.  When it fails, the *implied* breadth is returned, which is what exposes the
    KVLCC2 ``Bwl 580`` defect (the implied breadth is 58.0 m, i.e. the printed 580 is a
    factor-of-ten typo, and it travels with the table into every later edition).
    """
    out: dict[str, Any] = {}
    for ship in ships or sorted(report):
        cell = report[ship]["quantities"]
        volume = cell["displacement_volume"]["family_values"]
        cb = cell["cb"]["family_values"]
        lpp = cell["length_pp"]["family_values"]
        beam = cell["beam"]["family_values"]
        draft = cell["draft"]["family_values"]
        rows: dict[str, Any] = {}
        for family in sorted(volume):
            if not all(family in d for d in (cb, lpp, beam, draft)):
                continue
            product = cb[family] * lpp[family] * beam[family] * draft[family]
            rel = (product - volume[family]) / volume[family] if volume[family] else float("inf")
            implied = volume[family] / (cb[family] * lpp[family] * draft[family]) if cb[family] * lpp[family] * draft[family] else float("nan")
            rows[family] = {
                "cb_lpp_bwl_t_m3": product,
                "published_volume_m3": volume[family],
                "rel_error": rel,
                "consistent": abs(rel) <= 0.01,
                "implied_beam_m": implied,
            }
        out[ship] = rows
    return out


#: Closure identity tolerance.  The identity ``∇ = Cb·Lpp·Bwl·T`` is exact for a
#: wall-sided hull and holds to the printed rounding of a publisher table; 1 % is far
#: tighter than the 2.4 % error the ONRT cell shows when its mass/volume basis is wrong,
#: and far looser than the 0.08 % left by rounding.
CLOSURE_TOLERANCE = 0.01

#: Power-of-ten slips seen in printed tables (KVLCC2 ``Bwl 580`` for 58.0 m).  The scan
#: divides the printed number by each factor and asks whether the identity then closes.
SCALE_FACTORS: tuple[float, ...] = (10.0, 100.0, 0.1, 0.01)

#: The quantities the identity can attribute a slip to, in the order the scan tries them.
_SCALE_SCAN_ORDER: tuple[str, ...] = ("beam", "draft", "length_pp", "cb", "displacement_volume")


def _family_tables(report: Mapping[str, Any], ship: str) -> dict[str, dict[str, float]]:
    """``quantity -> {family: value}`` for one ship (the voting values only)."""
    return {q: dict(report[ship]["quantities"][q]["family_values"]) for q in QUANTITIES}


def _identity_residual(
    volume: float, cb: float, lpp: float, beam: float, draft: float
) -> float:
    """``Cb·Lpp·Bwl·T / ∇ - 1``; ``inf`` when a factor is missing or ∇ is zero."""
    product = cb * lpp * beam * draft
    if not volume:
        return float("inf")
    return product / volume - 1.0


def _scale_slip(
    tables: Mapping[str, dict[str, float]], family: str, raw_residual: float
) -> dict[str, Any] | None:
    """Find a single power-of-ten slip that makes the publication close on itself.

    The scan is deliberately narrow: exactly *one* quantity of *one* publisher is divided
    by 10, 100, 0.1 or 0.01, the cheapest factor first (10 before 100).  A table that needs
    two simultaneous corrections is not a digit slip and is left alone.  The printed number
    is never rewritten -- the caller records the factor next to it.
    """
    if abs(raw_residual) <= CLOSURE_TOLERANCE:
        return None
    for quantity in _SCALE_SCAN_ORDER:
        table = tables[quantity]
        if family not in table:
            continue
        printed = table[family]
        if not printed:
            continue
        for factor in SCALE_FACTORS:
            corrected = printed / factor
            trial = {q: dict(step) for q, step in tables.items()}
            trial[quantity][family] = corrected
            basis = next(
                (q for q in ("length_pp", "lwl") if family in trial[q]), None
            )
            if basis is None:
                continue
            residual = _identity_residual(
                trial["displacement_volume"][family],
                trial["cb"][family],
                trial[basis][family],
                trial["beam"][family],
                trial["draft"][family],
            )
            if abs(residual) <= CLOSURE_TOLERANCE:
                return {
                    "quantity": quantity,
                    "printed": printed,
                    "divide_by": factor,
                    "corrected": corrected,
                    "raw_residual": raw_residual,
                    "corrected_residual": residual,
                }
    return None


def adjudicate(report: Mapping[str, Any]) -> dict[str, Any]:
    """Decide *between* disagreeing publishers using only their own numbers.

    Two questions, both answered without an external reference value:

    **Is a published table internally consistent?**  Every publisher that states ∇ together
    with Cb, Lpp, Bwl and T implies a displacement of its own.  A family whose implied ∇
    misses its printed ∇ by more than :data:`CLOSURE_TOLERANCE` has a defect *in the table*,
    whatever the reconstruction does.

    **Which side of a dispute survives that test?**  For every ``CONFLICT`` cell the closure
    residual of each involved family is reported; when exactly one family closes, it is named
    as the ``preferred_family`` and the conflict is marked ``resolved_by_closure``.  When
    both close (a genuine difference of definition, e.g. LCB on different references) or
    neither does, the conflict stays ``unresolved`` and no value is picked: a library that
    guesses is worse than a library that shows both numbers.

    A failed closure is also scanned for a single power-of-ten slip (:func:`_scale_slip`),
    which is how the KVLCC2 ``Bwl 580`` cell is attributed to a printed typo rather than to a
    hull 10× as wide.  The slip is recorded on the cell itself
    (``cell['closure_defect']``) so a consumer can see why a family did not vote.
    """
    ships: dict[str, Any] = {}
    scale_typos: list[dict[str, Any]] = []
    conflicts: dict[str, dict[str, Any]] = {}
    for ship in sorted(report):
        tables = _family_tables(report, ship)
        volume = tables["displacement_volume"]
        rows: dict[str, Any] = {}
        for family in sorted(volume):
            needed = ("cb", "beam", "draft")
            if any(family not in tables[q] for q in needed):
                continue
            # Publishers define Cb on either length and print only one of them: NMRI's ONRT
            # row states ``CB = ∇/(Lwl·Bwl·T)`` and no Lpp.  Use Lpp when the table states
            # it and fall back to Lwl, recording which basis closed the identity.
            basis = next(
                (q for q in ("length_pp", "lwl") if family in tables[q]), None
            )
            if basis is None:
                continue
            length = tables[basis][family]
            product = (
                tables["cb"][family] * length * tables["beam"][family] * tables["draft"][family]
            )
            residual = _identity_residual(
                volume[family],
                tables["cb"][family],
                length,
                tables["beam"][family],
                tables["draft"][family],
            )
            slip = _scale_slip(tables, family, residual)
            rows[family] = {
                "length_basis": basis,
                "implied_volume_m3": product,
                "published_volume_m3": volume[family],
                "rel_error": residual,
                "consistent": abs(residual) <= CLOSURE_TOLERANCE,
                "explained_by_scale_slip": slip,
            }
            if slip is not None:
                scale_typos.append({"ship": ship, "family": family, **slip})
                # Bind the finding to the votes it came from: the entry keeps its printed
                # value and gains the reason it is not corroboration.
                for entry in report[ship]["quantities"][slip["quantity"]]["entries"]:
                    cell = entry["cell"]
                    if entry["family"] == family and cell["value"] == slip["printed"]:
                        cell["closure_defect"] = {
                            "quantity": slip["quantity"],
                            "printed": slip["printed"],
                            "divide_by": slip["divide_by"],
                            "corrected": slip["corrected"],
                            "corrected_residual": slip["corrected_residual"],
                        }
        ship_conflicts: dict[str, Any] = {}
        for quantity in QUANTITIES:
            row = report[ship]["quantities"][quantity]
            if row["verdict"] != "CONFLICT":
                continue
            involved = sorted(row["family_values"])
            residuals = {
                family: (rows[family]["rel_error"] if family in rows else None)
                for family in involved
            }
            consistent = [
                family
                for family in involved
                if residuals[family] is not None and abs(residuals[family]) <= CLOSURE_TOLERANCE
            ]
            preferred = consistent[0] if len(consistent) == 1 else None
            ship_conflicts[quantity] = {
                "families": involved,
                "family_values": dict(row["family_values"]),
                "closure_residual_by_family": residuals,
                "consistent_families": consistent,
                "preferred_family": preferred,
                "resolution": "resolved_by_closure" if preferred else "unresolved",
                "basis": (
                    f"{preferred} closes on its own Cb·Lpp·Bwl·T within "
                    f"{CLOSURE_TOLERANCE:.0%}; the other family does not"
                    if preferred
                    else "closure cannot separate the two readings (both close, or neither does)"
                ),
            }
        ships[ship] = {"closure": rows, "conflicts": ship_conflicts}
        conflicts.update({f"{ship}/{q}": row for q, row in ship_conflicts.items()})
    return {
        "closure_tolerance": CLOSURE_TOLERANCE,
        "scale_factors": list(SCALE_FACTORS),
        "ships": ships,
        "scale_typos": scale_typos,
        "conflicts": conflicts,
        "conflicts_resolved": sum(
            1 for row in conflicts.values() if row["resolution"] == "resolved_by_closure"
        ),
        "conflicts_unresolved": sum(
            1 for row in conflicts.values() if row["resolution"] == "unresolved"
        ),
    }


def coverage_summary(report: Mapping[str, Any]) -> dict[str, Any]:
    """Fleet-level tallies used by the gate and by the report."""
    ships = sorted(report)
    cells = [report[s]["quantities"][q] for s in ships for q in QUANTITIES]
    return {
        "ships": ships,
        "cells": len(cells),
        "cross_verified": sum(1 for c in cells if c["verdict"] == "CROSS_VERIFIED"),
        "single_source": sum(1 for c in cells if c["verdict"] == "SINGLE_SOURCE"),
        "conflict": sum(1 for c in cells if c["verdict"] == "CONFLICT"),
        "missing": sum(1 for c in cells if c["verdict"] == "MISSING"),
        "ships_fully_cross_verified": [
            s
            for s in ships
            if all(report[s]["quantities"][q]["verdict"] == "CROSS_VERIFIED" for q in QUANTITIES)
        ],
        "ships_failing_min_two_sources": [
            s
            for s in ships
            if any(
                report[s]["quantities"][q]["independent_families"] < 2
                for q in ("displacement_volume", "cb", "lcb_percent", "lwl")
            )
        ],
    }

#: A pinned digest is a SHA-256 hex string, nothing else.  ``pin_audit`` exists because a
#: digest that was typed rather than fetched still *looks* like a fingerprint: two of the
#: seven source digests in the first version of this registry shared their leading eight
#: hex digits with the real bytes and then diverged, which is what a hand-written
#: continuation looks like.  The audit recomputes every digest from the bytes and says so.
DIGEST_PATTERN = re.compile(r"^[0-9a-f]{64}$")


def pin_audit(fetched: Mapping[str, bytes | str]) -> dict[str, dict[str, Any]]:
    """Per source: is the pinned digest the digest of the bytes that were actually read?

    ``fetched`` maps source key -> payload.  A key that is not in ``fetched`` is reported as
    ``not_fetched``; it never counts as verified.  ``hand_written`` marks a pinned digest
    that is not even 64 lowercase hex digits, and ``drift`` marks a fetch whose digest no
    longer matches the pin -- a source change that must be re-read deliberately.
    """
    import hashlib

    out: dict[str, dict[str, Any]] = {}
    for spec in SOURCES:
        payload = fetched.get(spec.key)
        row: dict[str, Any] = {"url": spec.url, "family": spec.family}
        if payload is None:
            row["verdict"] = "not_fetched"
            row["pinned"] = spec.sha256
            out[spec.key] = row
            continue
        data = payload.encode("utf-8") if isinstance(payload, str) else payload
        digest = hashlib.sha256(data).hexdigest()
        row.update(pinned=spec.sha256, fetched=digest, bytes=len(data))
        if not spec.sha256:
            row["verdict"] = "unpinned"
        elif not DIGEST_PATTERN.match(spec.sha256):
            row["verdict"] = "hand_written"
        elif spec.sha256 == digest:
            row["verdict"] = "match"
        else:
            row["verdict"] = "drift"
        out[spec.key] = row
    return out


@dataclass(frozen=True)
class Reference:
    """One published quantity, resolved to a single value *and* to its provenance."""

    ship: str
    quantity: str
    value: float
    disputed: bool
    families: tuple[str, ...]
    #: Every family value, so a disputed reference can be quoted both ways.
    values: tuple[tuple[str, float], ...]
    #: ``(source_key, sha256)`` for every entry that backs the value.
    provenance: tuple[tuple[str, str], ...]

    def as_dict(self) -> dict[str, Any]:
        return {
            "ship": self.ship,
            "quantity": self.quantity,
            "value": self.value,
            "disputed": self.disputed,
            "families": list(self.families),
            "values": [{"family": fam, "value": value} for fam, value in self.values],
            "provenance": [{"source": key, "sha256": sha} for key, sha in self.provenance],
        }


def resolve_reference(
    report: Mapping[str, Any],
    ship: str,
    quantity: str = "displacement_volume",
    *,
    preferred_family: str | None = None,
) -> Reference | None:
    """Resolve one published quantity to a value plus provenance, never to an average.

    * agreed sources -> the agreed value;
    * one family -> that family's value, ``disputed=False`` (honest single source);
    * several families that disagree -> ``disputed=True`` with *all* values kept.  The
      returned ``value`` is the preferred family's, else the family with the most entries
      (ties broken alphabetically) -- a deterministic rule, not a mean, and the caller can
      see from ``disputed`` that the number is one of two.
    * nothing published -> ``None`` (the caller must fall back to a closure, visibly).
    """
    quantity_row = report[ship]["quantities"][quantity]
    families = quantity_row["family_values"]
    if not families:
        return None
    entries = quantity_row["entries"]
    counts: dict[str, int] = {}
    provenance: list[tuple[str, str]] = []
    for entry in entries:
        fam = entry["family"]
        counts[fam] = counts.get(fam, 0) + 1
        provenance.append((entry["cell"]["source"], entry["cell"]["sha256"]))
    if preferred_family is not None:
        if preferred_family not in families:
            raise KeyError(f"{preferred_family!r} states no {quantity} for {ship}")
        chosen = preferred_family
    elif len(families) == 1:
        chosen = next(iter(families))
    else:
        chosen = sorted(families, key=lambda fam: (-counts[fam], fam))[0]
    return Reference(
        ship=ship,
        quantity=quantity,
        value=families[chosen],
        disputed=quantity_row["verdict"] == "CONFLICT",
        families=tuple(sorted(families)),
        values=tuple(sorted(families.items())),
        provenance=tuple(sorted(set(provenance))),
    )
