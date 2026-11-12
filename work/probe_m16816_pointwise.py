from tensorcad.shear_flow import i_section, shear_flow
from tensorcad.shear_flow_multicell import MultiCellSection, multi_cell_shear_flow

sec = i_section(depth=3.0, flange_width=2.0, web_thickness=0.02, flange_thickness=0.03)
graph = MultiCellSection(sec.walls)
flow = multi_cell_shear_flow(graph, 12000.0)
merged = shear_flow(sec, 12000.0)
print("walls merged:", [(w.label, w.points) for w in sec.walls])
print("edges graph:", [(e.label, e.a, e.b) for e in graph._segments])
print("n merged edges", len(merged.edges), "n graph", len(flow.edge_flows))
for mine, edge in zip(flow.edge_flows, merged.edges):
    print(mine.label, edge.label)
    for f in (0.0, 0.25, 0.5, 0.75, 1.0):
        print(f"   s={f}: graph={mine.q_at(f):12.6f} merged={edge.q_at(f):12.6f}")
print("resultant graph", flow.resultant, "merged", merged.vertical_resultant, merged.horizontal_resultant)
