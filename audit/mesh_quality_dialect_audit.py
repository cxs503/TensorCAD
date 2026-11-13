#!/usr/bin/env python3
"""m17437 交差核验：GAP-CAD-D1-MESH-QUALITY-LIB 指标方言一致性审计。

核验对象：m17404 交付的 tensorcad.geomquality.quality（网格质量指标库）。

本脚本回答一个报告没回答的问题：新库的 aspect_ratio 与仓内既有的
mesh_conditioning._aspect_ratio / validation.measure_quality 是**同一个量的不同写法**，
还是**不同的量**？若是后者，两处共用的阈值 25 就不可迁移，门禁会静默换判。

运行（不写 m17404 工作树，只读）：
    PYTHONPATH=/root/tp-agent/work/m17404 python audit/mesh_quality_dialect_audit.py
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

from tensorcad import mesh_conditioning as mc
from tensorcad.geomquality.quality import mesh_quality, triangle_quality


def legacy_aspect(tri) -> float:
    """仓内既有方言：l_max^2 / (2A)（mesh_conditioning._aspect_ratio）。"""
    longest = max(math.dist(tri[0], tri[1]), math.dist(tri[1], tri[2]), math.dist(tri[2], tri[0]))
    ux, uy, uz = (tri[1][0] - tri[0][0], tri[1][1] - tri[0][1], tri[1][2] - tri[0][2])
    vx, vy, vz = (tri[2][0] - tri[0][0], tri[2][1] - tri[0][1], tri[2][2] - tri[0][2])
    nx, ny, nz = (uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx)
    area2 = math.sqrt(nx * nx + ny * ny + nz * nz)
    if area2 <= 0.0:
        return math.inf
    return longest * longest / area2


def h(banner: str) -> None:
    print("\n" + banner)


def closed_form_ratio_report() -> dict:
    """两种方言之比的闭式：geomquality/legacy = p/(2*sqrt(3)*l_max)，取值 (1/sqrt3, sqrt3/2]。"""
    out = {}
    cases = {
        "equilateral_a1": ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.5, math.sqrt(3) / 2, 0.0)),
        "right_3_4_5": ((0.0, 0.0, 0.0), (4.0, 0.0, 0.0), (0.0, 3.0, 0.0)),
        "isosceles_2_2_0.02": ((0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (1.0, 0.02, 0.0)),
    }
    for name, tri in cases.items():
        g = triangle_quality(*tri)
        leg = legacy_aspect(tri)
        a, b, c = (math.dist(tri[0], tri[1]), math.dist(tri[1], tri[2]), math.dist(tri[2], tri[0]))
        p = a + b + c
        lmax = max(a, b, c)
        out[name] = {
            "geomquality_aspect": g.aspect_ratio,
            "legacy_aspect": leg,
            "ratio": g.aspect_ratio / leg,
            "closed_form_p_over_2sqrt3_lmax": p / (2.0 * math.sqrt(3.0) * lmax),
        }
    return out


def single_element_counterexample() -> dict:
    """One triangle, two gates, opposite verdicts — no mesh, no threshold tuning.

    Family: base 2, height h.  legacy = l_max^2/(2A) = 4/(2h), canonical = l_max*p/(4sqrt3 A).
    h = 0.0741 puts legacy at 26.99 (> 25: the shipped decimation guard rejects it) while
    the canonical aspect is 15.60 (< 25: the new default gate accepts it), and every other
    canonical measure is healthy (min angle 4.24 deg, gamma ~ 0.29, edge ratio 1.99).
    """
    h = 0.0741
    tri = ((0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (1.0, h, 0.0))
    g = triangle_quality(*tri)
    leg = legacy_aspect(tri)
    scaled = tuple((3.7 * x, 3.7 * y, 3.7 * z) for x, y, z in tri)
    gs = triangle_quality(*scaled)
    return {
        "triangle": tri,
        "legacy_aspect_lmax2_over_2A": leg,
        "canonical_aspect": g.aspect_ratio,
        "ratio": g.aspect_ratio / leg,
        "min_angle_deg": g.min_angle_deg,
        "gamma": g.gamma,
        "edge_ratio": g.edge_ratio,
        "verdict_legacy_gate_limit25": "FAIL" if leg > 25.0 else "PASS",
        "verdict_canonical_gate_limit25": "FAIL" if g.aspect_ratio > 25.0 else "PASS",
        "scale_invariance_rel_diff": [
            abs(legacy_aspect(scaled) - leg) / leg,
            abs(gs.aspect_ratio - g.aspect_ratio) / g.aspect_ratio,
        ],
    }


def mesh_threshold_report() -> dict:
    """在真船网格上量：同一个 25 阈值，两种方言判出多少不同的元素。"""
    sys.path.insert(0, "/root/tp-agent/work/m17404/tools")
    import mesh_quality_crosscheck as xc  # noqa: E402  复用它建网格的同一路径

    meshes = xc._hull_meshes()
    out = {}
    for name, mesh in meshes.items():
        verts = mesh.vertices
        faces = mesh.faces
        n_disagree = 0
        worst_ratio = math.inf
        max_new = 0.0
        max_legacy = 0.0
        band = 0
        only_aspect = 0
        for f in faces:
            idx = list(f)[:3]
            tri = [tuple(float(x) for x in verts[i]) for i in idx]
            g = triangle_quality(*tri)
            leg = legacy_aspect(tri)
            max_new = max(max_new, g.aspect_ratio)
            max_legacy = max(max_legacy, leg)
            new_fail = g.aspect_ratio > 25.0
            leg_fail = leg > 25.0
            if new_fail != leg_fail:
                n_disagree += 1
            if leg > 25.0 and g.aspect_ratio <= 25.0:
                band += 1
                if leg > 0.0:
                    worst_ratio = min(worst_ratio, leg / g.aspect_ratio)
                # Does this element pass every OTHER library check?  If so the aspect
                # dialect alone decides its verdict.
                if (g.min_angle_deg >= 2.0 and g.gamma >= 0.05 and g.edge_ratio <= 100.0):
                    only_aspect += 1
        out[name] = {
            "triangles": len(faces),
            "max_aspect_geomquality": max_new,
            "max_aspect_legacy_mesh_conditioning": max_legacy,
            "elements_where_verdict_differs_at_limit_25": n_disagree,
            "elements_failing_legacy_but_passing_new": band,
            "worst_legacy_over_new_ratio": (None if worst_ratio is math.inf else worst_ratio),
            "of_which_pass_every_other_library_check": only_aspect,
        }
    return out


def main() -> int:
    out_dir = Path(__file__).resolve().parent.parent / "artifacts" / "mesh_quality_audit"
    out_dir.mkdir(parents=True, exist_ok=True)
    h("== 1. 闭式对比：geomquality.aspect_ratio 与仓内既有 aspect_ratio 是否同一个量 ==")
    closed = closed_form_ratio_report()
    print(json.dumps(closed, indent=1))

    h("== 2. 单元素反例（无网格、无阈值调参）==")
    single = single_element_counterexample()
    print(json.dumps(single, indent=1))

    h("== 3. 真船网格上，共用阈值 25 的判决分歧 ==")
    meas = mesh_threshold_report()
    print(json.dumps(meas, indent=1))

    h("== 4. 结论 ==")
    ratios = [v["ratio"] for v in closed.values()]
    verdict = {
        "closed_form_ratio_range": [min(ratios), max(ratios)],
        "identity_holds_below_1e-12": all(abs(r - v["closed_form_p_over_2sqrt3_lmax"]) < 1e-12
                                          for r, v in zip(ratios, closed.values())),
        "same_quantity": False,
        "single_element_counterexample": single,
        "note": (
            "legacy = l_max^2/(2A)；geomquality = l_max*p/(4*sqrt(3)*A) = legacy * p/(2*sqrt(3)*l_max)，"
            "比值域 (1/sqrt3, sqrt3/2] = (0.577, 0.866]，恒 < 1。两者是不同归一化的量，"
            "所以同一个阈值 25 不可迁移：旧方言判超限的元素可被新库判合规。"
        ),
    }
    print(json.dumps(verdict, indent=1))

    artifact = {
        "audit_of": "m17404 / GAP-CAD-D1-MESH-QUALITY-LIB",
        "audit_by": "m17437 board-audit",
        "closed_form": closed,
        "mesh_measurements": meas,
        "verdict": verdict,
        "repro_cmd": "PYTHONPATH=/root/tp-agent/work/m17404 python3 audit/mesh_quality_dialect_audit.py",
    }
    out = out_dir / "dialect.json"
    out.write_text(json.dumps(artifact, indent=1, ensure_ascii=False) + "\n")
    print(f"artifact: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())