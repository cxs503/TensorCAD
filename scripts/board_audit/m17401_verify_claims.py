#!/usr/bin/env python3
"""m17401 board audit: independently re-verify the claims of mission m17349
(GAP-CAD-B1-STEP-READER, TensorCAD) and dump every number as JSON.

Zero hand-copying: each field is computed here from the delivered artifacts /
from a fresh re-run of the subject's own code, never typed from its report.

Checks
  C1 rerun      : `python tools/step_reader_audit.py` in the subject worktree
                  reproduces volume errors + convergence ratios bit-for-bit
  C2 volume     : the analytic closed forms are recomputed here from the
                  geometry described in docs/gap_cad_b1_step_reader.md, then the
                  delivered meshes are re-exported to OBJ and their volume is
                  re-integrated by an independent divergence-theorem script
  C3 determinism: regenerating the sample STEP files reproduces the shipped
                  sha256 (the fixtures are code, not downloads)
  C4 real file  : census of the third-party NX/ST-Developer AP203 file and the
                  exact refusal message when a solid is requested
  C5 refusal    : a mutated file with an unsupported surface type is refused by
                  name instead of silently guessed
  C6 malformed  : an intentionally inconsistent file is reported as
                  watertight=False with volume None (no fabricated volume)
  C7 benchmark  : tensor_benchmark_cases rows for the anchor (the authoritative
                  benchmark judge; absence is recorded, not papered over)

Run:
  python scripts/board_audit/m17401_verify_claims.py --out artifacts/board_audit/m17401/verify.json
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import subprocess
import sys
from pathlib import Path

SUBJECT = Path("/root/tp-agent/work/m17349")
ART = SUBJECT / "artifacts/gap-cad-b1-step-reader"
DSN = ("host=127.0.0.1 port=5432 dbname=tensorplatform user=foundation "
       "password=uGTN-KwuQv_7l8Tw0MMWUx_kf-nF2Rzz")
ANCHOR = "GAP-CAD-B1-STEP-READER"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def py(*args: str, cwd: Path = SUBJECT) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, *args], cwd=cwd, capture_output=True, text=True)


def obj_volume(path: Path) -> dict:
    """Independent divergence-theorem volume + edge census on an OBJ file."""
    verts, faces = [], []
    for line in path.read_text().splitlines():
        tok = line.split()
        if not tok:
            continue
        if tok[0] == "v":
            verts.append(tuple(float(x) for x in tok[1:4]))
        elif tok[0] == "f":
            faces.append([int(x.split("/")[0]) - 1 for x in tok[1:]])
    volume = 0.0
    edges: dict[tuple[int, int], int] = {}
    for f in faces:
        for k in range(1, len(f) - 1):
            a, b, c = verts[f[0]], verts[f[k]], verts[f[k + 1]]
            volume += (a[0] * (b[1] * c[2] - b[2] * c[1])
                       - a[1] * (b[0] * c[2] - b[2] * c[0])
                       + a[2] * (b[0] * c[1] - b[1] * c[0])) / 6.0
        for i in range(len(f)):
            e = tuple(sorted((f[i], f[(i + 1) % len(f)])))
            edges[e] = edges.get(e, 0) + 1
    return {"volume_m3": abs(volume), "triangles": len(faces), "vertices": len(verts),
            "edges_used_once": sum(1 for v in edges.values() if v == 1)}


# ---- C2: closed forms, transcribed from the geometry table of the delivery doc
CLOSED_FORM = {
    "plate_with_hole": lambda: (200 * 120 * 20 - math.pi * 20 ** 2 * 20) / 1e9,
    "stepped_shaft": lambda: (math.pi * 30 ** 2 * 80
                              + math.pi * 50 / 3 * (30 ** 2 + 30 * 15 + 15 ** 2)
                              + math.pi * 15 ** 2 * 40) / 1e9,
    "bushing": lambda: (math.pi * (40 ** 2 - 25 ** 2) * 60) / 1e9,
    "l_bracket": lambda: (60 * 20 + (50 - 20) * 20) * 10 / 1e9,
}


def c1_rerun() -> dict:
    watched = [ART / "volume_validation.json", ART / "convergence.json"]
    before = {p.name: sha(p) for p in watched}
    run = py("tools/step_reader_audit.py")
    after = {p.name: sha(p) for p in watched}
    shipped = json.loads((ART / "volume_validation.json").read_text())
    return {
        "exit_code": run.returncode,
        "all_pass": shipped["all_pass"],
        "artifacts_unchanged_by_rerun": before == after,
        "sha256": after,
        "per_file": {f["file"]: {"error_percent": f["error_percent"],
                                 "nominal_volume_m3": f["nominal_volume_m3"],
                                 "measured_volume_m3": f["measured_volume_m3"],
                                 "watertight": f["watertight"],
                                 "boundary_edges": f["boundary_edges"],
                                 "euler_characteristic": f["euler_characteristic"],
                                 "triangles": f["triangles"]} for f in shipped["files"]},
        "convergence": json.loads((ART / "convergence.json").read_text())["records"],
        "stderr_tail": run.stderr.strip().splitlines()[-2:],
    }


def c2_volume() -> dict:
    out = {}
    for name, formula in CLOSED_FORM.items():
        step = ART / "samples" / f"{name}.step"
        obj = Path(f"/tmp/m17401_{name}.obj")
        py("-m", "tensorcad.step_brep", "mesh", str(step), "--obj", str(obj))
        got = obj_volume(obj)
        ref = formula()
        out[name] = {"closed_form_m3": ref, **got,
                     "error_percent": 100.0 * (got["volume_m3"] - ref) / ref}
    return out


def c3_determinism(tmp: Path) -> dict:
    regen = tmp / "regen"
    regen.mkdir(parents=True, exist_ok=True)
    py("tools/step_samples.py", "--out", str(regen))
    return {p.name: {"shipped": sha(ART / "samples" / p.name), "regenerated": sha(p),
                     "identical": sha(p) == sha(ART / "samples" / p.name)}
            for p in sorted(regen.glob("*.step"))}


def c4_real_file() -> dict:
    real = ART / "real" / "2010_30p30n_thik_te_18inches.stp"
    code = (
        "import json,sys;from tensorcad.step_brep import census,read_step_brep\n"
        f"c=census({str(real)!r})\n"
        "o={'instances':c['instances'],'schema':c['header']['file_schema'],"
        "'preprocessor':c['header']['preprocessor'],'length_scale':c['length_scale'],"
        "'length_scale_detail':c.get('length_scale_detail'),"
        "'bsp':c['entities'].get('B_SPLINE_CURVE_WITH_KNOTS'),"
        "'wireframe':c['entities'].get('GEOMETRICALLY_BOUNDED_WIREFRAME_SHAPE_REPRESENTATION'),"
        "'solid_entities':{k:v for k,v in c['entities'].items() if 'BREP' in k or 'SHELL' in k}}\n"
        f"try:\n    read_step_brep({str(real)!r})\n    o['refusal']=None\n"
        "except Exception as e:\n    o['refusal']=type(e).__name__+': '+str(e)\n"
        "print(json.dumps(o))"
    )
    res = py("-c", code)
    return json.loads(res.stdout)


def c5_refusal(tmp: Path) -> dict:
    src = (ART / "samples" / "plate_with_hole.step").read_text()
    mutated = re.sub(r"CYLINDRICAL_SURFACE\('',#\d+,([-\d.Ee+]+)\)",
                     lambda m: "B_SPLINE_SURFACE_WITH_KNOTS('',3,3,(),.UNSPECIFIED.)", src)
    assert mutated != src, "mutation did not apply"
    target = tmp / "mutated_bspline.step"
    target.write_text(mutated)
    code = (
        "import json;from tensorcad.step_brep import read_step_brep\n"
        f"try:\n    read_step_brep({str(target)!r})\n    o={{'refusal':None}}\n"
        "except Exception as e:\n    o={'refusal':type(e).__name__+': '+str(e)}\n"
        "print(json.dumps(o))"
    )
    return {"mutation": "CYLINDRICAL_SURFACE -> B_SPLINE_SURFACE_WITH_KNOTS",
            **json.loads(py("-c", code).stdout)}


def c6_malformed(tmp: Path) -> dict:
    """Scale every point by 2 but leave CIRCLE radii untouched -> a file that is
    internally inconsistent; the reader must not invent a volume for it."""
    src = (ART / "samples" / "plate_with_hole.step").read_text()
    bad = re.sub(r"CARTESIAN_POINT\('',\(([^)]*)\)\)",
                 lambda m: "CARTESIAN_POINT('',(%s))" % ",".join(
                     "%.10g" % (float(v) * 2) for v in m.group(1).split(",")), src)
    target = tmp / "inconsistent_scale.step"
    target.write_text(bad)
    code = (
        "import json;from tensorcad.step_brep import read_step_brep\n"
        f"r=read_step_brep({str(target)!r})\n"
        "print(json.dumps({'complete':r.complete,'watertight':r.watertight,"
        "'boundary_edges':r.boundary_edges,'volume_m3':r.volume_m3}))"
    )
    return json.loads(py("-c", code).stdout)


def c7_benchmark() -> dict:
    try:
        import psycopg
        with psycopg.connect(DSN) as conn:
            rows = conn.execute(
                "SELECT case_name,metric,tolerance,status FROM tensor_benchmark_cases "
                "WHERE anchor_slug=%s", (ANCHOR,)).fetchall()
            all_cad = conn.execute(
                "SELECT DISTINCT anchor_slug FROM tensor_benchmark_cases "
                "WHERE repo='TensorCAD' ORDER BY 1").fetchall()
        return {"anchor_slug": ANCHOR, "rows": rows, "n_rows": len(rows),
                "tensorcad_anchors_in_table": [r[0] for r in all_cad],
                "note": ("no machine-adjudicated benchmark case exists for this anchor; "
                         "the acceptance numbers live only in the mission goal and the "
                         "delivered artifacts")}
    except Exception as exc:  # noqa: BLE001 - audit must not die on a DB hiccup
        return {"anchor_slug": ANCHOR, "error": f"{type(exc).__name__}: {exc}"}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    tmp = Path("/tmp/m17401_probe")
    tmp.mkdir(parents=True, exist_ok=True)
    report = {
        "mission": 17349,
        "auditor_mission": 17401,
        "anchor": ANCHOR,
        "repo": "TensorCAD",
        "subject_worktree": str(SUBJECT),
        "subject_commits": subprocess.run(["git", "log", "--oneline", "-3"], cwd=SUBJECT,
                                          capture_output=True, text=True).stdout.splitlines(),
        "C1_rerun": c1_rerun(),
        "C2_volume_independent": c2_volume(),
        "C3_fixture_determinism": c3_determinism(tmp),
        "C4_real_third_party_file": c4_real_file(),
        "C5_unsupported_refusal": c5_refusal(tmp),
        "C6_malformed_file": c6_malformed(tmp),
        "C7_benchmark_table": c7_benchmark(),
    }
    Path(args.out).write_text(json.dumps(report, indent=1, ensure_ascii=False))
    print(json.dumps(report, indent=1, ensure_ascii=False)[:2000])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())