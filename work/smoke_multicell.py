import math, sys
sys.path.insert(0, ".")
from tensorcad.shear_flow import Wall, ThinWalledSection, shear_flow, box_section, i_section, channel_section
from tensorcad.shear_flow_multicell import (MultiCellSection, multi_cell_box_section,
                                            multi_cell_shear_flow, multi_cell_shear_centre)

# 1. two-cell box: topology
s = multi_cell_box_section([2.0, 2.0], 3.0, 0.02)
print("cells", s.cell_count, "cycle", s.cycle_rank, "areas", [s.cell_area(i) for i in range(s.cell_count)])
f = multi_cell_shear_flow(s, 8000.0)
print("equil", f.equilibrium_residual, "transv", f.transverse_residual, "cont", f.continuity_residual,
      "twist", f.twist_residual, "resultant", f.resultant)
sc = multi_cell_shear_centre(s)
print("shear centre", sc.y, sc.z, "centroid", s.centroid_y, s.centroid_z)

# 2. limit: inner web vanishing -> module single-cell box of width 4, height 3
t_in = 1e-9
s2 = multi_cell_box_section([2.0, 2.0], 3.0, 0.02, inner_thickness=t_in)
f2 = multi_cell_shear_flow(s2, 8000.0)
ref = shear_flow(ThinWalledSection([Wall(((0.0,0.0),(4.0,0.0)),0.02,"b"),
                                    Wall(((4.0,0.0),(4.0,3.0)),0.02,"s"),
                                    Wall(((4.0,3.0),(0.0,3.0)),0.02,"t"),
                                    Wall(((0.0,3.0),(0.0,0.0)),0.02,"p")]), 8000.0)
print("module single-cell box corners:", [(e.label, e.q_at(0.0), e.q_at(1.0)) for e in ref.edges])
print("multicell corner flows:", [(e.label, round(e.q_at(0.0),6), round(e.q_at(1.0),6)) for e in f2.edge_flows])

# 3. open tree: I-beam against the merged module
i_sec = i_section(depth=3.0, flange_width=2.0, web_thickness=0.02, flange_thickness=0.03)
mc = MultiCellSection(i_sec.walls)
print("I-beam cells", mc.cell_count)
a = multi_cell_shear_flow(mc, 12000.0); b = shear_flow(i_sec, 12000.0)
print("module I-beam resultant", b.resultant, "mc", a.resultant)
print("module first edges", [(e.label, round(e.q_at(0.5),6)) for e in b.edges])
print("mc first edges    ", [(e.label, round(e.q_at(0.5),6)) for e in a.edge_flows])

# 4. channel shear centre
ch = channel_section(depth=3.0, flange_width=1.0, thickness=0.02)
mcc = MultiCellSection(ch.walls)
scm = multi_cell_shear_centre(mcc); scb = None
from tensorcad.shear_flow import shear_centre
scb = shear_centre(ch)
print("channel shear centre module", scb.y, scb.z, "multicell", scm.y, scm.z)
print("formula", 3*1.0**2/(6*1.0+3.0), "centroid", ch.centroid_y)

# 5. module refuses the two-cell graph
try:
    ThinWalledSection(s.walls)
except Exception as exc:
    print("module refuses 2-cell:", type(exc).__name__, str(exc)[:60])
