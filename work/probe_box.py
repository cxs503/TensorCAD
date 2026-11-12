from tensorcad.shear_flow import Wall, ThinWalledSection, shear_flow
from tensorcad.shear_flow_multicell import multi_cell_box_section, multi_cell_shear_flow

V = 8000.0
def reference_box():
    return ThinWalledSection([
        Wall(((0.0,0.0),(4.0,0.0)),0.02,"bottom"),
        Wall(((4.0,0.0),(4.0,3.0)),0.02,"side_s"),
        Wall(((4.0,3.0),(0.0,3.0)),0.02,"top"),
        Wall(((0.0,3.0),(0.0,0.0)),0.02,"side_p"),
    ])
section = multi_cell_box_section([4.0], 3.0, 0.02)
print("graph nodes/edges", [(e.label,e.a,e.b) for e in section._segments])
print("nodes", section._nodes)
print("cells", section.cells)
graph = multi_cell_shear_flow(section, V)
mod = shear_flow(reference_box(), V)
for e in mod.edges:
    print("merged", e.label, e.start, e.end, e.q_at(0), e.q_at(0.5), e.q_at(1))
for e in graph.edge_flows:
    print("graph ", e.label, e.start, e.end, e.q_at(0), e.q_at(0.5), e.q_at(1))
print("merged resultant", mod.vertical_resultant, mod.horizontal_resultant)
print("graph resultant", graph.resultant)
print("merged centre", end=" ")
from tensorcad.shear_flow import shear_centre
print(shear_centre(reference_box()).to_dict())
