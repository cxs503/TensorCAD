import sys; sys.path.insert(0, "tests")
from tensorcad.shear_flow import Wall, ThinWalledSection, shear_flow, AXIS_Z
from tensorcad.shear_flow_multicell import multi_cell_box_section, multi_cell_shear_flow
from shear_flow_reduction import PerimeterOracle, field_at, points_along, difference_on_points, flow_scale

V = 8000.0
walls_ccw = [
    Wall(((0.0,0.0),(4.0,0.0)),0.02,"bottom"),
    Wall(((4.0,0.0),(4.0,3.0)),0.02,"side_s"),
    Wall(((0.0,3.0),(4.0,3.0)),0.02,"top"),     # traversed -y by the loop
    Wall(((0.0,3.0),(0.0,0.0)),0.02,"side_p"),
]
ref = ThinWalledSection([Wall(((0.0,0.0),(4.0,0.0)),0.02,"bottom"),
                         Wall(((4.0,0.0),(4.0,3.0)),0.02,"side_s"),
                         Wall(((4.0,3.0),(0.0,3.0)),0.02,"top"),
                         Wall(((0.0,3.0),(0.0,0.0)),0.02,"side_p")])
graph = multi_cell_box_section([4.0], 3.0, 0.02)
merged = shear_flow(ref, V)
gflow = multi_cell_shear_flow(graph, V)

# oracle chain must be continuous; walk CCW: bottom (0,0)->(4,0), side_s up, top y:4->0, side_p down
chain = [
    Wall(((0.0,0.0),(4.0,0.0)),0.02,"bottom"),
    Wall(((4.0,0.0),(4.0,3.0)),0.02,"side_s"),
    Wall(((4.0,3.0),(0.0,3.0)),0.02,"top"),
    Wall(((0.0,3.0),(0.0,0.0)),0.02,"side_p"),
]
for direction in (1,-1):
    for cut_wall, cut_f in ((2,0.5), (0,0.25)):
        o = PerimeterOracle(chain, shear_force=V, coordinate=lambda p: p[1], cut_wall=cut_wall,
                            cut_fraction=cut_f, direction=direction)
        pts = points_along(chain, count=7)
        print(f"dir={direction} cut={cut_wall}@{cut_f}: I={o.second_moment:.12f} zc={o.centroid} "
              f"circ={o.circulation:.6f} res={o.resultant} merged-vs-oracle={difference_on_points(merged, o, pts):.3e} "
              f"graph-vs-oracle={difference_on_points(gflow, o, pts):.3e} scale={flow_scale(merged):.1f}")
print("hand formula: corner", V*0.02*2*1.5/0.45, " web mid", V*(0.06+0.02*(1.5*1.5-0.5*1.5**2))/0.45)
print("merged corner", field_at(merged,(4.0,3.0)), "graph corner", field_at(gflow,(4.0,3.0)))
print("merged web mid", field_at(merged,(4.0,1.5)), "graph", field_at(gflow,(4.0,1.5)))
print("flange mid merged", field_at(merged,(2.0,3.0)), field_at(merged,(2.0,0.0)))
print("I merged", merged.second_moment, "graph", graph.second_moment_horizontal)
