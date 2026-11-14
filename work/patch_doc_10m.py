import shutil
from pathlib import Path

root = Path('/srv/tensor/agent/work/m20965')
src = Path('/tmp/cj3081x/cj3081/repo/artifacts/perf_scale/ladder_fleet_10m.json')
for rel in ('benchmarks/perf_scale/ladder_fleet_10m.json', 'artifacts/perf_scale/ladder_fleet_10m.json'):
    dst = root / rel
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(src, dst)
    print('copied ->', rel)

p = root / 'docs/perf_scale_10m.md'
s = p.read_text(encoding='utf-8')
start = s.index('### 8.3 ')
new = '''### 8.3 10⁷ 带实测：已拿到（`cj3081`，n330, rc=0）

上一版本的这一节是一份 honest FAIL，并把根因写成「取件口径：作业写在 worktree 里的新文件不会被 fetch 打包」。
**这个诊断是错的**，本单已把它推翻（见 §8.4 第 1 条）：两次失败作业的日志原文是
`bash: line 1: cd: /srv/tensor/agent/work/m20813: Permission denied`（rc=1）——作业根本没启动，
被测工作树属于另一张单（m20813）且节点侧无读权限。取件侧一切正常：本轮 `cj3081` 的
`artifacts/perf_scale/ladder_fleet_10m.json` 原样出现在 fetch 包里。

命令（经 `run_compute` 递交，cores=4，作业墙钟 523 s）：

```
python scripts/perf_scale_10m.py --base 2001:501 --copies 2,3 --n-queries 2000 \\
    --output artifacts/perf_scale/ladder_fleet_10m.json
```

节点 n330（LoongArch64, 128 核），实测时 load1≈49/128。

| 子任务 | 体数 | 三角面 | 秒 | 吞吐 | 建索引 s | 索引 MB | 峰值 RSS MB |
|---|---:|---:|---:|---:|---:|---:|---:|
| `load`（生成单体） | 1 | 4,003,998 | 321.14 | 12,468 tri/s | – | – | 5,770 |
| `fleet_assemble` | 2 | 8,007,996 | 6.52 | 1.23 M tri/s | – | – | 5,770 |
| `point_locate_accel`（2000 射线） | 2 | 8,007,996 | 1.285 | 1,557 query/s | 68.31 | 637 | 5,770 |
| `volume` | 2 | 8,007,996 | 4.01 | 2.00 M tri/s | – | – | 5,770 |
| `fleet_assemble` | 3 | 12,011,994 | 9.50 | 1.26 M tri/s | – | – | 5,770 |
| `point_locate_accel`（2000 射线） | 3 | **12,011,994** | **1.247** | **1,604 query/s** | 100.26 | 954 | 5,770 |
| `volume` | 3 | **12,011,994** | **6.05** | 1.99 M tri/s | – | – | 5,770 |

索引规模：x2 档 2,030,644 cells、entries/tri=1.37；x3 档 3,030,075 cells、entries/tri=1.36
（近立方单元选址在细长多体包围盒上未爆桶）。原始数据：`benchmarks/perf_scale/ladder_fleet_10m.json`
（同一份也留在 `artifacts/perf_scale/`），作业 stdout 另有 gzip+base64 payload（§8.5）。

**精度（关键量 vs 解析基准，逐项）**

| 量 | 基准（`copies × analytic_volume(10,2,0.8)`） | 实测 | 相对误差 | <3%? |
|---|---|---:|---:|---|
| 体积 @8,007,996 tri（2 体） | 15.025396 m³ | 15.024947 m³ | −0.0030% | 是 |
| 体积 @12,011,994 tri（3 体） | 22.538095 m³ | 22.537421 m³ | −0.0030% | 是 |

**10⁷ 带的判定（取代 §5 的外推）**：12.0e6 三角面上，2000 条射线经网格加速只需 **1.25 s**。
按 §3.1 实测的暴力吞吐（7.9 query/s @411k，近似 ∝1/N）外推，同一查询集需
`2000 / (7.9 × 411k/12.0M) ≈ 7.4e3 s ≈ 2.1 h`，即加速结构在 10⁷ 带把点定位成本压低 **≈5.9e3 倍**。
结论更强：**10⁷ 三角面的点定位不需要集群**，单节点 1.25 s。真正的单机瓶颈是
(a) 参数化生成器（4.0e6 用 321 s，10⁷ 约 800 s）与 (b) 一次性建索引（12.0e6 用 100 s / 954 MB），
而不是查询本身。

**回归门（诚实口径）**：fleet 阶梯只有一个基础档 × 两个体数，每个 `scaling_*` 门都只有 1 个规模点，
因此全部 `not_run`，本次运行 `verdict=not_run`——**不是 pass**。§8.2 的 26 门全绿来自 quick 阶梯
（4 档、指数可拟合），两套数据不可混用。

### 8.4 本单修正的两个结论

1. **「10⁷ 取不回来」的根因在提交侧，不在取件侧**：`cj2967`/`cj2976` 的作业日志均为
   `cd: /srv/tensor/agent/work/m20813: Permission denied`、rc=1，作业一行都没执行。教训：递交前核对
   `stage_dir` 是否本工单工作树；失败要读 `jobs/cjN.log` 的 `rc` 与 `output_tail`，
   不能用未经证据支持的归因（原「取件口径」说法已被本轮 fetch 包直接证伪）。
2. **`not_run` 门曾被计入「全绿」**：`GateVerdict(passed=True, status=not_run)` 让
   `PerfBenchmarkResult.passed` 为 True、CLI 打印 "All gates passed."，而实际一个门都没算。
   现新增 `gates_checked` / `gates_not_run` / `verdict`（`pass`/`fail`/`not_run` 三态）并进 JSON；
   CLI 在无门可算时打印 "This is NOT a pass."。反证测试：
   `tests/test_perf_benchmark.py::test_verdict_distinguishes_not_checked_from_passed`、
   `::test_cli_does_not_claim_pass_when_no_gate_was_checked`（后者跑单档阶梯，断言输出含 "NOT a pass"）。

### 8.5 结果传输：三条独立通道（`tensorcad/perf_scale_artifact.py`）

一次 5–10 分钟的重计算不该死于打包约定，结果走三条互不依赖的通道：

1. `--output` JSON 落在**受版本控制**的 `benchmarks/perf_scale/ladder_fleet_10m.json`
   （`artifacts/` 部分子路径会被发布门禁滤掉，故基准产物不放那里）；
2. `PERF10M_JSON_BEGIN/END` 之间的 gzip+base64 payload（按 900 字符分块），可用
   `python -m tensorcad.perf_scale_artifact <job.log> -o recovered.json` 从日志还原；
3. 人读表格（子任务/体数/三角面/秒/建索引/索引 MB/RSS）。

`parse_emitted` 对**截断**报错（缺块=错误，绝不返回部分结果）；`tests/test_perf_scale_artifact.py`
（8 项）覆盖往返无损、分块、丢块、缺标记、坏 base64、日志恢复与失败退出码。

**局限（诚实记录）**
- `peak_rss_mb` 在 x2/x3 档与 `load` 档完全相同（5,770 MB），而 12.0e6 三角面的 Python 元组
  按体积推断不该停在这个量级；该字段取进程 `ru_maxrss` 高水位，在 loongarch64/n330 上疑似
  报告粒度或释放时机问题。**本单不基于 RSS 给出 10⁷ 内存结论**，需单独标定。
- 体数只有 2 和 3（无第三个规模点），所以复杂度门 `not_run`。要拿到 10⁷ 带的指数门，
  需再补一档，例如 `--base 1281:321 --copies 2,4,8`（1.6e6/3.3e6/6.6e6/1.3e7）。
'''
s = s[:start] + new
p.write_text(s, encoding='utf-8')
print('docs updated', len(s))