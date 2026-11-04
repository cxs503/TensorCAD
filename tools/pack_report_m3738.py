import base64, io, json, pathlib
from PIL import Image

root = pathlib.Path("/root/tp-agent/work/m3738")
tok = json.loads(pathlib.Path("/tmp/login_m3738.json").read_text())["access_token"]
caps = {
    "test_comparison.png": "逐块 pytest 通过数：合并树(m3738) vs main 基线(8389c7e)；新增 96 艉板 + 30 Wigley 测试全绿，失败集合与基线逐条相同。",
    "line_preservation.png": "4 个冲突文件的解析构成：并集新增行 / 重排连接行 / 被并集发布号取代的版本字面量；两侧丢失行数均为 0。",
}
figs = []
for name, cap in caps.items():
    im = Image.open(root / "docs/m3738" / name).convert("RGB")
    im.thumbnail((1000, 1000))
    buf = io.BytesIO()
    im.save(buf, "PNG", optimize=True)
    figs.append({"name": name, "caption": cap,
                 "png_base64": base64.b64encode(buf.getvalue()).decode()})
rep = json.loads(pathlib.Path("/tmp/report_m3738.json").read_text())
rep["figures"] = figs
out = pathlib.Path("/tmp/report_m3738_small.json")
out.write_text(json.dumps(rep))
print("payload_KB", round(out.stat().st_size / 1024, 1), "token_len", len(tok))
pathlib.Path("/tmp/token_m3738.txt").write_text(tok)