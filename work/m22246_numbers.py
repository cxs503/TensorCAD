#!/usr/bin/env python3
"""Render the m22246 numbers section of ``docs/volume_mesh_format_bridge.md`` from the JSON.

Every figure in section 7 of that document is written by this script, read out of
``runs/m19950_volume_bridge/bridge_matrix.json`` (itself produced by
``examples/meshio_volume_bridge/run_example.py``).  Nothing is transcribed by hand, so the
document cannot drift away from the run it claims to describe -- the failure mode this whole
mission exists to close.

Usage::

    python examples/meshio_volume_bridge/run_example.py     # writes the JSON
    python work/m22246_numbers.py                           # rewrites section 7's tables
"""
from __future__ import annotations

import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
RUNS = REPO / "runs" / "m19950_volume_bridge" / "bridge_matrix.json"
DOC = REPO / "docs" / "volume_mesh_format_bridge.md"
BEGIN = "<!-- BEGIN m22246 numbers -->"
END = "<!-- END m22246 numbers -->"


def fmt(value: float, digits: int = 3) -> str:
    """Scientific with a fixed number of significant digits, and exact zeros kept as 0."""
    if value == 0.0:
        return "0"
    return f"{value:.{digits}e}"


def main() -> int:
    payload = json.loads(RUNS.read_text())
    rows = payload["gmsh"]
    lines: list[str] = []

    lines.append("| body | groups | verts | tets | boundary tri | V file [m³] | V tets [m³] | rel | watertight |")
    lines.append("|---|---|---|---|---|---|---|---|---|")
    for row in rows:
        lines.append(
            f"| {row['body']} | {' '.join(row['names'])} | {row['vertices']} | {row['tets']} | "
            f"{row['triangles']} | {row['volume_file']:.12f} | {row['volume_tets']:.12f} | "
            f"{fmt(row['rel_volume'])} | {'yes' if row['watertight'] else 'no'} |"
        )

    lines.append("")
    lines.append("Readers of the same file (each an independent implementation):")
    lines.append("")
    lines.append("| body | this module | meshio | FEM handoff (`solver_interface`) | cells seen by handoff |")
    lines.append("|---|---|---|---|---|")
    for row in rows:
        lines.append(
            f"| {row['body']} | {fmt(row['rel_volume'])} | {fmt(row['meshio_rel_volume'])} | "
            f"{fmt(row['handoff_rel_volume'])} | {row['handoff_tets']} tet / "
            f"{row['handoff_triangles']} tri |"
        )

    lines.append("")
    lines.append("Physical-group areas against the closed-form / quadrature anchors:")
    lines.append("")
    lines.append("| body | patch | area [m²] | anchor [m²] | rel |")
    lines.append("|---|---|---|---|---|")
    for row in rows:
        for patch, value in row["area_rel"].items():
            lines.append(
                f"| {row['body']} | `{patch}` | {row['areas'][patch]:.9f} | "
                f"{row['area_reference'][patch]:.9f} | {fmt(value)} |"
            )

    lines.append("")
    lines.append(
        "`all_pass` of the run that produced these numbers: "
        f"**{payload['all_pass']}** (criteria in the JSON: `gmsh_volume_rel` < "
        f"{fmt(payload['criteria']['gmsh_volume_rel'], 0)}, `gmsh_area_rel` < "
        f"{payload['criteria']['gmsh_area_rel']})."
    )

    rendered = "\n".join(lines)
    (REPO / "work" / "m22246_numbers.md").write_text(rendered + "\n", encoding="utf-8")

    text = DOC.read_text(encoding="utf-8")
    if BEGIN in text and END in text:
        head, rest = text.split(BEGIN, 1)
        _stale, tail = rest.split(END, 1)
        DOC.write_text(f"{head}{BEGIN}\n{rendered}\n{END}{tail}", encoding="utf-8")
        print(f"section 7 tables rewritten from {RUNS.relative_to(REPO)}")
    else:
        print(f"markers not found in {DOC}; wrote work/m22246_numbers.md instead")
        print(rendered)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())