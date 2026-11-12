import sys; sys.path.insert(0, "tests")
from tensorcad.shear_flow import Wall, ThinWalledSection, shear_flow
from tensorcad.shear_flow_multicell import multi_cell_box_section, multi_cell_shear_flow
from shear_flow_reduction import PerimeterOracle, field_at, points_along

V=8000.0
ref = ThinWalledSection([Wall(((0.0,0.0),(4.0,0.0)),0.02,"bottom"),
                         Wall(((4.0,0.0),(4.0,3.0)),0.02,"side_s"),
                         Wall(((4.0,3.0),(0.0,3.0)),0.02,"top"),
                         Wall(((0.0,3.0),(0.0,0.0)),0.02,"side_p")])
merged = shear_flow(ref, V)
chain = [Wall(((0.0,0.0),(4.0,0.0)),0.02,"bottom"),
         Wall(((4.0,0.0),(4.0,3.0)),0.02,"side_s"),
         Wall(((4.0,3.0),(0.0,3.0)),0.02,"top"),
         Wall(((0.0,3.0),(0.0,0.0)),0.02,"side_p")]
o = PerimeterOracle(chain, shear_force=V, coordinate=lambda p: p[1], cut_wall=2, cut_fraction=0.5)
for p in points_along(chain, count=3):
    print(f"{p}  merged={field_at(merged,p):14.6f} oracle={o.q_at(p):14.6f}")
