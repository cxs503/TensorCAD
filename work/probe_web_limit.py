from tensorcad.shear_flow import Wall, ThinWalledSection, shear_flow
from tensorcad.shear_flow_multicell import multi_cell_box_section, multi_cell_shear_flow

V=8000.0
def ref():
    return ThinWalledSection([
        Wall(((0.0,0.0),(4.0,0.0)),0.02,"bottom"),
        Wall(((4.0,0.0),(4.0,3.0)),0.02,"side_s"),
        Wall(((4.0,3.0),(0.0,3.0)),0.02,"top"),
        Wall(((0.0,3.0),(0.0,0.0)),0.02,"side_p"),
    ])

mod = shear_flow(ref(), V)
ref_side = {e.label: e.q_at(0.5) for e in mod.edges}
ref_flange_mid_top = [e.q_at(0.5) for e in mod.edges if e.label=="top"][0]
print("merged side q(0.5):", ref_side, " top q(0.5)", ref_flange_mid_top)
for inner in (0.02, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10):
    sec = multi_cell_box_section([2.0,2.0], 3.0, 0.02, inner_thickness=inner)
    f = multi_cell_shear_flow(sec, V)
    ss = f.edge("box.side_s")[0].q_at(0.5)
    sp = f.edge("box.side_p")[0].q_at(0.5)
    web = f.edge("box.web1")[0].q_at(0.5)
    tp = f.edge("box.top0")[0]
    print(f"inner={inner:8.0e} side_s={ss:14.8f} ({abs(ss-ref_side['side_s']):11.3e}) side_p={sp:14.8f} ({abs(abs(sp)-abs(ref_side['side_p'])):11.3e}) web={web:12.6f} top0(s=0)={tp.q_at(0):12.6f} resultant={f.resultant[1]:12.6f}")
