# Series 60 family, round 2: the capability on a product entry

Anchor `GAP-CAD-T2-PARAM-SERIES` (`TensorCAD`) · mission `m19694`. Every number on this page is rendered from `runs/m19694_series60_family/round2_evidence.json` by `tools/series60_family/round2_page.py`; nothing here is transcribed.

**What this round adds.** Two previous rounds, and the board audits of both, landed on the same sentence: the kernel is real and the artifacts reproduce, but the family sits behind an import, so no consumer can reach it, and the only "published comparison" was a closure labelled as a verification. Round 2 fixes both: a product entry (`tensorcad-hull series60 parent|family`), exports wired into the package namespace, and a comparison whose role labels say which rows are closures. None of that is an external anchor for the hull's *shape* -- that is still absent, and §6 says so.

## 1. Case and reference

The Series 60 (Todd 1951) is the standard systematic single-screw displacement-hull family; five parent forms carry the published particulars `Cb`, `Cm`, `Cp`, `Cwp`, `LCB/L`. The generator lofts

```
y(xi, zeta) = (B/2) * f_deck(xi) * zeta**(1/m(xi)),   xi = x/L, zeta = z/T
f_deck(xi)  = 1 - |2 xi - 1|^n(side)   (aft/fore exponents split by side)
```

| quantity | value |
|---|---|
| length `L` | 10.0 m |
| beam `B` | 1.5 m |
| draft `T` | 0.6 m |
| `L/B`, `B/T` | 6.667, 2.500 |
| `Cb` parents | 0.60, 0.65, 0.70, 0.75, 0.80 |
| grid ladder | 81x21 -> 161x41 -> 321x81 |
| second parameter axis | transom breadth ratio `tau` in 0.00, 0.15 at `Cb = 0.70` |

Reference table: `docs/data/series60_published_parents.json` (sha256 `16ac564989a87d2b...`). Citation carried by that file:

> Todd, F.H. -- Series 60 methodical series of single-screw merchant ship models (DTMB); parent-form particulars as tabulated in the resistance/design literature and cited in this repository as 'Todd 1951'

Provenance, stated plainly: the table is documentary -- it was inherited from this repository's own docs (no network in the sandbox, no fresh fetch), and it is self-consistent to within 0.4 % of `Cb` on `Cb = Cp * Cm`. Its sha256 is recorded in the evidence file so any edit is visible.

## 2. Method and implementation

| layer | where | what it does |
|---|---|---|
| kernel | `tensorcad/series60.py` | station law, inverse solve for `n_aft/n_fore`, closed forms, lofter, bulb/transom ends |
| **entry (new)** | `tensorcad/series60_cli.py` | `parent` (generate -> export -> gate) and `family` (five parents vs the published table); no constants of its own |
| dispatch (new) | `tensorcad/cli.py` | `tensorcad-hull series60 ...` reaches the above |
| exports (new) | `tensorcad/__init__.py` | `PARENT_CB`, `Series60Parameters`, `generate_series60`, `analytic_series60_hydrostatics`, ... |
| evidence | `tools/series60_family/round2_evidence.py` | drives the entry by subprocess, writes JSON/CSV/PNG |
| gate power | `tools/hull_gate_power/mutate.py` | one defect at a time, exact-string patch, tests that guard the claim |

The closed forms that make the comparison decidable:

```
Cp = n/(n+1),  Cm = m/(m+1),  Cb = Cp*Cm,  VCB/T = (m+1)/(2m+1),
Cwp = 2 * int_0^1 (1 - |2 xi - 1|^n) dxi       (transom term added separately)
```

Three checks keep the numbers honest:

1. the entry hands the *published* `Cwp` and `LCB` to the generator as targets and then re-measures them off the mesh -- that is an inverse-problem closure, and it is labelled `input target (closure check)` in every artifact;
2. the mesh is integrated on its own (`compute_hydrostatics`, closed surface) and compared with the closed form -- a discretisation check, bar 2 %;
3. the gate's reference volume comes from the analytic volume, never from the mesh (mutant `S60-17` breaks exactly that and is killed).

## 3. Parameter table (dimensionless)

Ladder grid `321x81`; coefficients are the measured mesh values, `closed` is the module's analytic value.

| `Cb` target | `Cb` closed | `Cb` mesh | `Cp` | `Cm` | `Cwp` | `LCB/L` | `VCB/T` |
|---|---|---|---|---|---|---|---|
| 0.60 | 0.600000 | 0.597959 | 0.626021 | 0.955175 | 0.699992 | 0.475041 | -0.462732 |
| 0.65 | 0.650000 | 0.647449 | 0.670940 | 0.964987 | 0.726991 | 0.483033 | -0.471657 |
| 0.70 | 0.700000 | 0.696836 | 0.715768 | 0.973550 | 0.752990 | 0.490021 | -0.480158 |
| 0.75 | 0.750000 | 0.746110 | 0.760491 | 0.981089 | 0.779988 | 0.496007 | -0.487768 |
| 0.80 | 0.800000 | 0.795209 | 0.805049 | 0.987777 | 0.806986 | 0.501999 | -0.494833 |

## 4. Comparison: published vs got, per criterion

Grid `321x81`. The `role` column is the part previous rounds got wrong: only `Cb` is imposed and re-measured; `Cwp` and `LCB/L` are imposed targets, so agreement is closure; `Cp` is the published ladder read to pin the exponents.

| `Cb` | quantity | role | published | got (mesh) | error | within ±3 % |
|---|---|---|---|---|---|---|
| 0.60 | block_coefficient | derived | 0.600000 | 0.597959 | -0.340 % | yes |
| 0.60 | lcb_over_length | input target (closure check) | 0.475000 | 0.475041 | +0.009 % | yes |
| 0.60 | midship_coefficient | published ladder (input, not a prediction) | 0.960000 | 0.955175 | -0.503 % | yes |
| 0.60 | prismatic_coefficient | published ladder (input, not a prediction) | 0.625000 | 0.626021 | +0.163 % | yes |
| 0.60 | waterplane_coefficient | input target (closure check) | 0.700000 | 0.699992 | -0.001 % | yes |
| 0.65 | block_coefficient | derived | 0.650000 | 0.647449 | -0.393 % | yes |
| 0.65 | lcb_over_length | input target (closure check) | 0.483000 | 0.483033 | +0.007 % | yes |
| 0.65 | midship_coefficient | published ladder (input, not a prediction) | 0.970000 | 0.964987 | -0.517 % | yes |
| 0.65 | prismatic_coefficient | published ladder (input, not a prediction) | 0.670000 | 0.670940 | +0.140 % | yes |
| 0.65 | waterplane_coefficient | input target (closure check) | 0.727000 | 0.726991 | -0.001 % | yes |
| 0.70 | block_coefficient | derived | 0.700000 | 0.696836 | -0.452 % | yes |
| 0.70 | lcb_over_length | input target (closure check) | 0.490000 | 0.490021 | +0.004 % | yes |
| 0.70 | midship_coefficient | published ladder (input, not a prediction) | 0.980000 | 0.973550 | -0.658 % | yes |
| 0.70 | prismatic_coefficient | published ladder (input, not a prediction) | 0.715000 | 0.715768 | +0.107 % | yes |
| 0.70 | waterplane_coefficient | input target (closure check) | 0.753000 | 0.752990 | -0.001 % | yes |
| 0.75 | block_coefficient | derived | 0.750000 | 0.746110 | -0.519 % | yes |
| 0.75 | lcb_over_length | input target (closure check) | 0.496000 | 0.496007 | +0.001 % | yes |
| 0.75 | midship_coefficient | published ladder (input, not a prediction) | 0.985000 | 0.981089 | -0.397 % | yes |
| 0.75 | prismatic_coefficient | published ladder (input, not a prediction) | 0.760000 | 0.760491 | +0.065 % | yes |
| 0.75 | waterplane_coefficient | input target (closure check) | 0.780000 | 0.779988 | -0.001 % | yes |
| 0.80 | block_coefficient | derived | 0.800000 | 0.795209 | -0.599 % | yes |
| 0.80 | lcb_over_length | input target (closure check) | 0.502000 | 0.501999 | -0.000 % | yes |
| 0.80 | midship_coefficient | published ladder (input, not a prediction) | 0.990000 | 0.987777 | -0.225 % | yes |
| 0.80 | prismatic_coefficient | published ladder (input, not a prediction) | 0.805000 | 0.805049 | +0.006 % | yes |
| 0.80 | waterplane_coefficient | input target (closure check) | 0.807000 | 0.806986 | -0.002 % | yes |

Criterion ledger (verdict = the number, not the prose):

| criterion | number | bar | verdict |
|---|---|---|---|
| five parents generated through the entry | 5 parents | 5 | PASS |
| measured vs published particulars | -0.658 % (Cb=0.70 measured midship_coefficient) | ±3 % | PASS |
| mesh vs closed form at `321x81` | +0.605 % | 2 % | PASS |
| second parameter axis (transom) | tau=0.00: Cwp=0.7530, A_T/A_m=0.0000, tau=0.15: Cwp=0.7530, A_T/A_m=0.1076 | gate passes | PASS |
| family reachable from `tensorcad-hull` | 5 entry runs | all exit 0/1 | PASS |

![five Cb gradient](../runs/m19694_series60_family/figures/five_cb_gradient.png)

*Bars: coefficients generated through the entry; black dashes: the published value. Right: residuals, with the two closure rows separated by the legend labels.*

## 5. Convergence and parameter sensitivity

### 5.1 Three-grid ladder, five parents

| grid | worst |mesh − closed form| | worst vs published |
|---|---|---|
| 81x21 | +2.443 % | -2.443 % |
| 161x41 | +1.215 % | -1.233 % |
| 321x81 | +0.605 % | -0.658 % |

The ladder is monotone in the discretisation error: **True** (+0.605 % at `321x81`). The residual is the `zeta**(1/m)` keel singularity: `m` runs from 24 (`Cb = 0.60`) to 160 (`Cb = 0.80`), so the section area sits in a layer of thickness `T/m` under the waterline and a coarse vertical grid clips it. This is why the entry refuses a 3-level grid rather than reporting it.

![convergence](../runs/m19694_series60_family/figures/convergence_ladder.png)

### 5.2 Second parameter axis: transom breadth ratio

| `tau` | `Cb` | `Cwp` | immersed transom area / midship area | gate |
|---|---|---|---|---|
| 0.00 | 0.700000 | 0.753000 | 0.000000 | pass |
| 0.15 | 0.700000 | 0.753000 | 0.107625 | pass |

`Cwp` and `Cb` are held to machine precision while `tau` moves: the stern fills by reallocating the aft taper, not by changing the volume. The transverse transom *area* ratio is not the breadth ratio (`tau = 0.15` gives about 0.11 of the midship area) because the aft section fullness is not the midship fullness -- the artifact reports both, and conflating them is a known way to mis-report `A_T/A_m`.

![transom sweep](../runs/m19694_series60_family/figures/transom_sweep.png)

## 6. Gate power, failure modes, artifacts

### 6.1 Gate power (mutation testing)

This round re-measured the four **entry** mutants: `4/4` killed = **100.0 %** against the 80 % bar (`PASS`). The family has 18 single-defect mutants defined in `tools/hull_gate_power/mutate.py`; the 14 kernel mutants (`S60-1..S60-14`, transom blend and bulb closed forms) carry their last recorded verdict in `docs/series60_gate_power.md` and were **not** re-run in this round -- the full set is submitted as a scheduler job, so until that artifact is fetched their status here is *unchanged*, not re-confirmed.

| mutant | claim it breaks | verdict | wall |
|---|---|---|---|
| `S60-15-entry-ignores-cb` | the entry generates the block coefficient the caller asked for | KILLED | 7.1 s |
| `S60-16-entry-vcb-frame` | the measured VCB is compared in the frame the closed form states (mesh z-frame against z-frame), not against the above-keel value | KILLED | 5.5 s |
| `S60-17-entry-self-referential-gate` | the gate's reference volume comes from outside the mesh; a mesh that only agrees with itself is not evidence | KILLED | 5.8 s |
| `S60-18-entry-closure-as-prediction` | Cwp and LCB are generator targets, and the record says so | KILLED | 6.8 s |

### 6.2 Failure modes this round found or avoided

| failure mode | symptom | guard now in place |
|---|---|---|
| CAP/FILL: kernel with no consumer | three audits in a row: "not exported, not in `cli.py`, not in `hull.py`" | `tensorcad-hull series60` dispatch + package exports; mutants S60-15..18 |
| a closure reported as a verification | `role` field said `independent of the target` for `Cp`/`Cm` | roles fixed in `series60_cli.py` and `tools/series60_family/audit.py`; mutant S60-18 |
| frame mismatch | z-frame `VCB` compared with the above-keel closed form reads −193 % on every hull | the pair is explicit; mutant S60-16 |
| gate that cannot refuse | a failed gate that still writes files leaves a consumer with an untrustworthy artifact | gate is evaluated before export, a refusal writes nothing |
| self-referential evidence | reference volume taken from the mesh | reference comes from the analytic volume; mutant S60-17 |

### 6.3 Artifacts and reproduction

```
python3 tools/series60_family/round2_evidence.py --out runs/m19694_series60_family
python3 -m tensorcad.cli series60 family --published docs/data/series60_published_parents.json --nx 321 --nz 81 --out /tmp/family
python3 tools/hull_gate_power/mutate.py --out runs/m19694_series60_family/gate_power_series60_round2.json --only <S60 ids>
python3 -m pytest tests/test_series60_cli.py -q          # 7 passed
```

| artifact | content |
|---|---|
| `runs/m19694_series60_family/round2_evidence.json` | criteria, per-grid records, transom points, published sha256 |
| `runs/m19694_series60_family/figures/convergence_ladder.png` | figure |
| `runs/m19694_series60_family/figures/five_cb_gradient.png` | figure |
| `runs/m19694_series60_family/figures/transom_sweep.png` | figure |
| `runs/m19694_series60_family/five_cb_gradient_ladder.csv` | 3 grids x 5 parents x 2 sources x 5 coefficients, long format |
| `runs/m19694_series60_family/gate_power_series60_round2.json` | mutation run, per mutant |

### 6.4 What is *not* established

- **external anchor for the shape: WITHHELD.** a Todd station-offset table for the five parents: absent in-repo and no network in this sandbox, so the family's *shape* is not externally anchored. Consequence: the coefficient closure and the discretisation ladder are verified; the sectional-area distribution is NOT.
- The published table is documentary (inherited, not freshly fetched). Its internal consistency is audited (`Cb` vs `Cp·Cm`, ≤ 0.4 % of `Cb`) but it is not an issuing-body offset table.
- `tensor_benchmark_cases` holds no registered case for this anchor, so the 3 %/2 % bars are code-declared constants, not platform-registered tolerances. The evidence file names both bars and the number that decides each one.
- Only the hull is integrated: a bulb, when requested, is a separate closed solid and is *added*, not boolean-unioned, so the gate covers the hull mesh alone.

Overall: `all_passed = True`.
