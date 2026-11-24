import json
import base64
import os

api_base = os.environ.get("TENSOR_API_BASE", "http://127.0.0.1:8010/api/v1")
user = os.environ.get("TENSOR_AGENT_USER")
password = os.environ.get("TENSOR_AGENT_PASS")

fig_dir = "/root/tp-agent/work/m25806/artifacts/m23297_stl_encoding"
figures = []
for f in ["stl_encoding_envelope.png", "stl_ascii_cost.png", "stl_normal_consistency.png"]:
    path = os.path.join(fig_dir, f)
    if os.path.exists(path):
        with open(path, "rb") as fh:
            figures.append({
                "name": f,
                "caption": f,
                "png_base64": base64.b64encode(fh.read()).decode("utf-8")
            })

report = {
    "title": "STL 腿的双编码矩阵与 float32 量化包络（m23297）",
    "summary": "TensorCAD 完成 STL 格式双编码（ASCII/binary）的往返保真测量与 float32 量化包络验证。纠正了文档中关于 binary STL 相对误差上限的错误陈述，实测证明其为绝对误差且随坐标偏移增长。新增 stl_quantisation_bound 预检、write_stl_binary 规范写器（修正面法向缺陷），并通过 40 个组合的独立解析器复算验证。",
    "scenario_key": "GAP-CAD-B2-MESH-FORMAT-MATRIX",
    "mission_id": "25806",
    "tags": ["TensorCAD", "meshio", "STL", "float32", "quantisation"],
    "markdown": open("/root/tp-agent/work/m25806/docs/m23297_stl_encoding_verification.md").read(),
    "figures": figures
}

with open("/root/tp-agent/work/m25806/runs/report.json", "w") as f:
    json.dump(report, f)
print("Report JSON written to runs/report.json")
