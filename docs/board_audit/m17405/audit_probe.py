#!/usr/bin/env python3
"""Independent audit probe for m17357 (GAP-CAD-D3-GEOM-CONTRACT).

Does NOT reuse the subject's test fixtures: every mesh below is built by this
script, and every reference number is derived from an analytic formula or from
the declared generator parameters, so nothing is hand-copied from the subject's
report.  Target worktree is passed on the command line (subject checkout).

Checks
------
A. Analytic anchors      : STL signed volume vs closed-form volume for the
                           box barge, the unit sphere and the Wigley hull
                           -> quantifies whether the validator measures real
                           geometry and reports honest numbers.
B. Red-team detection    : 5 violations built independently here -> does
                           validate_contract flag the right convention?
C. False-negative probes : cases where the contract deliberately stays silent
                           (blind windows) -> bounds the claim "units/orientation
                           violations are auto-detected".
D. CLI contract          : exit code 0 on pass / 1 on fail, JSON parseable.
E. Reproducibility       : rerun tools/run_contract_compliance.py, compare the
                           regenerated compliance_report.json byte-for-byte.

Output: artifacts/geom_contract_audit/probe_results.json  (programmatic, no
hand-typed numbers anywhere in the report).
"""
from __future__ import annotations

import hashlib
import json
import math
import subprocess
import sys
from pathlib import Path

SUBJECT = Path(sys.argv[1] if len(sys.argv) > 1 else "/root/tp-agent/work/m17357")
OUT = Path(__file__).resolve().parents[3] / "artifacts" / "geom_contract_audit"
sys.path.insert(0, str(SUBJECT))

from tensorcad import mesh as meshlib                       # noqa: E402
from tensorcad.export import read_stl                       # noqa: E402
from tensorcad.geometry import Mesh, MeshBuilder            # noqa: E402
from tensorcad.geom_contract import validate_contract        # noqa: E402

ASSETS = SUBJECT / "artifacts" / "geom_contract_assets"


# --------------------------------------------------------------------------- #
# A. analytic anchors
# --------------------------------------------------------------------------- #
def anchor_checks() -> list[dict]:
    """got = signed volume recomputed from the shipped STL; ref = closed form."""
    rows = []
    # closed box: L x B x D  (generator: make_box_barge(length=30, breadth=8, draft=4, depth=5))
    box = read_stl(ASSETS / "box_barge.stl", units="m")
    ref_box = 30.0 * 8.0 * 5.0
    got_box = meshlib.total_signed_volume(box)
    rows.append(dict(case="box_barge closed box", ref=ref_box, got=got_box,
                     formula="V = L*B*D = 30*8*5"))

    # sphere r = 1 (generator: make_sphere(radius=1.0)); faceted, so expect a
    # small negative bias from the inscribed tessellation.
    sph = read_stl(ASSETS / "sphere_r1.stl", units="m")
    ref_sph = 4.0 / 3.0 * math.pi
    got_sph = meshlib.total_signed_volume(sph)
    rows.append(dict(case="sphere r=1 faceted (cube-to-sphere, 32x64)", ref=ref_sph,
                     got=got_sph, formula="V = 4/3*pi*r^3"))

    # Wigley parabolic hull (generator: L=10, B=1, T=0.625)
    #   y(x,z) = (B/2)(1-(2x/L)^2)(1-(z/T)^2)
    #   V = B * [L - (4/L^2)L^3/12] * [T - T/3] = (4/9) L B T
    wig = read_stl(ASSETS / "wigley.stl", units="m")
    L, B, T = 10.0, 1.0, 0.625
    ref_wig = (4.0 / 9.0) * L * B * T
    got_wig = meshlib.total_signed_volume(wig)
    rows.append(dict(case="Wigley parabolic hull L=10 B=1 T=0.625", ref=ref_wig,
                     got=got_wig, formula="V = (4/9)*L*B*T"))

    for r in rows:
        r["rel_err_pct"] = (r["got"] - r["ref"]) / r["ref"] * 100.0
        r["abs_err_pct"] = abs(r["rel_err_pct"])
        r["within_3pct"] = r["abs_err_pct"] < 3.0
        r["within_1pct"] = r["abs_err_pct"] < 1.0
    return rows


def prismatic_observation() -> list[dict]:
    """Out-of-scope observation: Series 60 STL volume vs Cb*L*B*T.

    The validator is only reporting what is in the mesh; a deviation here is an
    asset-generator property, recorded so the board knows it is not a validator
    artifact.
    """
    rows = []
    for name, cb in (("series60_cb60", 0.60), ("series60_cb75", 0.75)):
        m = read_stl(ASSETS / f"{name}.stl", units="m")
        ref = cb * 20.0 * 3.0 * 1.5
        got = meshlib.total_signed_volume(m)
        rows.append(dict(case=f"{name} vs Cb*L*B*T (L=20 B=3 T=1.5)", ref=ref, got=got,
                         rel_err_pct=(got - ref) / ref * 100.0,
                         implied_cb=got / (20.0 * 3.0 * 1.5), declared_cb=cb))
    return rows


# --------------------------------------------------------------------------- #
# B/C. in-memory probes (fixtures written here, not the subject's)
# --------------------------------------------------------------------------- #
def tetra(sign: float = 1.0, units: str = "m", frame: str = "right-handed",
          scale: float = 1.0) -> Mesh:
    """Closed tetrahedron; sign=-1 reverses the winding (inward normals)."""
    b = MeshBuilder(units=units, coordinate_frame=frame)
    p = [(0.0, 0.0, 0.0), (scale, 0.0, 0.0), (0.0, scale, 0.0), (0.0, 0.0, scale)]
    v = [b.add_vertex(pt) for pt in p]
    tris = [(0, 2, 1), (0, 1, 3), (0, 3, 2), (1, 2, 3)]
    for a, c, d in tris:
        if sign > 0:
            b.add_triangle(v[a], v[c], v[d])
        else:
            b.add_triangle(v[a], v[d], v[c])
    return b.build()


def open_pair(inconsistent: bool) -> Mesh:
    """Two triangles sharing no edge: a deliberately open surface."""
    b = MeshBuilder(units="m", coordinate_frame="right-handed")
    a0 = b.add_vertex((0.0, 0.0, 0.0)); a1 = b.add_vertex((1.0, 0.0, 0.0))
    a2 = b.add_vertex((0.0, 1.0, 0.0))
    b.add_triangle(a0, a1, a2)                       # normal +z
    c0 = b.add_vertex((5.0, 0.0, 0.0)); c1 = b.add_vertex((6.0, 0.0, 0.0))
    c2 = b.add_vertex((5.0, 1.0, 0.0))
    if inconsistent:
        b.add_triangle(c0, c2, c1)                   # normal -z
    else:
        b.add_triangle(c0, c1, c2)
    return b.build()


def red_team() -> list[dict]:
    cases = {
        "inward normals (closed tetra, reversed winding)": tetra(sign=-1.0),
        "left-handed declared frame": tetra(frame="x, y, z down; left-handed"),
        "unrecognised unit token": tetra(units="inch"),
        "below-minimum scale (1 mm tetra declared m)": tetra(scale=1e-3),
        "above-maximum scale (1000 m tetra)": tetra(scale=1000.0),
    }
    out = []
    for label, m in cases.items():
        rep = validate_contract(m, source=label)
        out.append(dict(case=label, passed=rep.passed,
                        conventions=[v.convention for v in rep.violations],
                        severities=[v.severity for v in rep.violations],
                        messages=[v.message for v in rep.violations]))
    return out


def false_negatives() -> list[dict]:
    """Counterexamples: violations inside the contract's blind windows."""
    out = []

    # 1. 100x unit mislabel that lands inside the diagonal band.
    m = tetra(units="cm", scale=100.0)   # 100 "cm" == 1 m true size
    rep = validate_contract(m, source="unit mislabel cm/m, 1 m object")
    out.append(dict(probe="unit mislabelled 100x but diagonal still inside "
                          "[MIN_DIAGONAL_M, MAX_DIAGONAL_M]", passed=rep.passed,
                    violations=[v.convention for v in rep.violations],
                    mesh_info_diag_m=rep.mesh_info.get("diagonal_m"),
                    note="contract detects a unit error only when the declared "
                         "diagonal leaves the plausible band"))

    # 2. exactly 50% inconsistent winding on an open surface -> no violation.
    rep = validate_contract(open_pair(True), source="open surface, 50% reversed")
    out.append(dict(probe="open surface with 50% of faces reversed",
                    passed=rep.passed,
                    violations=[v.convention for v in rep.violations],
                    note="_check_orientation flags only ratio < 0.5, so an exact "
                         "half-reversed open surface passes clean"))

    # 3. left-handed *geometry* declared right-handed, open surface.
    b = MeshBuilder(units="m", coordinate_frame="right-handed")
    p = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)]
    v = [b.add_vertex(pt) for pt in p]
    b.add_triangle(v[0], v[1], v[2])
    rep = validate_contract(b.build(), source="open surface, frame string lies")
    out.append(dict(probe="handedness cannot be verified geometrically; a mesh "
                          "labelled 'right-handed' is accepted on the string alone",
                    passed=rep.passed,
                    violations=[v.convention for v in rep.violations],
                    note="_check_coordinate_system is a declaration-string check "
                         "(documented in the source docstring)"))
    return out


# --------------------------------------------------------------------------- #
# D/E. CLI and reproducibility
# --------------------------------------------------------------------------- #
def cli_checks() -> list[dict]:
    py = sys.executable
    out = []
    # NOTE: an STL carries no unit metadata, so the CLI can only judge what the
    # caller declares.  Each row therefore states the declared units explicitly.
    for label, asset, declared, expect in (
        ("box_barge declared m -> exit 0", "box_barge.stl", "m", 0),
        ("inverted_wigley declared m -> exit 1", "inverted_wigley.stl", "m", 1),
        ("wrong_units (metres of geometry) declared mm -> exit 1",
         "wrong_units.stl", "mm", 1),
        ("wrong_units declared m -> exit 0 (declaration matches geometry)",
         "wrong_units.stl", "m", 0),
    ):
        r = subprocess.run([py, "-m", "tensorcad.geom_contract",
                            str(ASSETS / asset), "--units", declared, "--json"],
                           cwd=SUBJECT, capture_output=True, text=True, timeout=300)
        try:
            parsed = json.loads(r.stdout)
            n = len(parsed); ok_json = True
        except Exception:
            n = None; ok_json = False
        out.append(dict(case=label, declared_units=declared, expected_exit=expect,
                        exit=r.returncode,
                        exit_ok=r.returncode == expect, json_ok=ok_json,
                        reports=n, stderr=r.stderr[-200:]))
    return out


def repro_check() -> dict:
    path = ASSETS / "compliance_report.json"
    before = hashlib.sha256(path.read_bytes()).hexdigest()
    r = subprocess.run([sys.executable, "tools/run_contract_compliance.py"],
                       cwd=SUBJECT, capture_output=True, text=True, timeout=600)
    after = hashlib.sha256(path.read_bytes()).hexdigest()
    data = json.loads(path.read_text())
    return dict(cmd="python3 tools/run_contract_compliance.py", exit=r.returncode,
                sha256_before=before, sha256_after=after,
                identical=before == after, n_reports=len(data),
                n_pass=sum(1 for d in data if d["passed"]),
                n_fail=sum(1 for d in data if not d["passed"]))


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    result = dict(
        subject_worktree=str(SUBJECT),
        anchors=anchor_checks(),
        prismatic_observation_out_of_scope=prismatic_observation(),
        red_team=red_team(),
        false_negatives=false_negatives(),
        cli=cli_checks(),
        repro=repro_check(),
    )
    (OUT / "probe_results.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())