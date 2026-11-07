"""Exchange tests: STEP / IGES / native round trips, plus a literal STEP fixture.

The fixture matters: a round trip proves that the reader and the writer agree, not that
either of them agrees with STEP.  :data:`STEP_FIXTURE` is a hand-written AP203 fragment
(one planar B-spline face, four boundary edges) that no code in the package produced.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from tensorcad.geomquality import baselines as baselines_module
from tensorcad.geomquality.brep import audit_topology
from tensorcad.geomquality.hydro import patch_hydrostatics
from tensorcad.geomquality.io import (
    ExchangeError,
    load_cad,
    read_iges,
    read_step,
    round_trip,
    save_cad,
)

STEP_FIXTURE = """ISO-10303-21;
HEADER;
FILE_DESCRIPTION(('hand written fixture'),'2;1');
FILE_NAME('fixture.step','2026-11-07T00:00:00',(''),(''),'','','');
FILE_SCHEMA(('CONFIG_CONTROL_DESIGN'));
ENDSEC;
DATA;
#1=CARTESIAN_POINT('',(0.,0.,0.));
#2=CARTESIAN_POINT('',(0.,0.,1.));
#3=CARTESIAN_POINT('',(1.,0.,0.));
#4=CARTESIAN_POINT('',(1.,0.,1.));
#5=VERTEX_POINT('',#1);
#6=VERTEX_POINT('',#2);
#7=VERTEX_POINT('',#3);
#8=VERTEX_POINT('',#4);
#10=B_SPLINE_CURVE_WITH_KNOTS('',1,(#1,#3),.UNSPECIFIED.,.F.,.F.,(2,2),(0.,1.),.UNSPECIFIED.);
#11=B_SPLINE_CURVE_WITH_KNOTS('',1,(#3,#4),.UNSPECIFIED.,.F.,.F.,(2,2),(0.,1.),.UNSPECIFIED.);
#12=B_SPLINE_CURVE_WITH_KNOTS('',1,(#4,#2),.UNSPECIFIED.,.F.,.F.,(2,2),(0.,1.),.UNSPECIFIED.);
#13=B_SPLINE_CURVE_WITH_KNOTS('',1,(#2,#1),.UNSPECIFIED.,.F.,.F.,(2,2),(0.,1.),.UNSPECIFIED.);
#20=EDGE_CURVE('',#5,#7,#10,.T.);
#21=EDGE_CURVE('',#7,#8,#11,.T.);
#22=EDGE_CURVE('',#8,#6,#12,.T.);
#23=EDGE_CURVE('',#6,#5,#13,.T.);
#30=ORIENTED_EDGE('',*,*,#20,.T.);
#31=ORIENTED_EDGE('',*,*,#21,.T.);
#32=ORIENTED_EDGE('',*,*,#22,.T.);
#33=ORIENTED_EDGE('',*,*,#23,.T.);
#40=EDGE_LOOP('',(#30,#31,#32,#33));
#50=FACE_OUTER_BOUND('',#40,.T.);
#60=B_SPLINE_SURFACE_WITH_KNOTS('',1,1,((#1,#2),(#3,#4)),.UNSPECIFIED.,.F.,.F.,.F.,(2,2),(2,2),(0.,1.),(0.,1.),.UNSPECIFIED.);
#70=ADVANCED_FACE('fixture_face',(#50),#60,.T.);
ENDSEC;
END-ISO-10303-21;
"""


def test_step_fixture_parses_without_the_writer(tmp_path: Path) -> None:
    path = tmp_path / "fixture.step"
    path.write_text(STEP_FIXTURE)
    brep = read_step(path)
    assert len(brep.faces) == 1
    assert len(brep.edges) == 4
    face = brep.faces[0]
    assert face.name == "fixture_face"
    assert face.domain == ((0.0, 1.0), (0.0, 1.0))
    # A unit square in the plane y = 0, parametrised as x = u, z = v.
    assert face.surface.evaluate(0.5, 0.5) == pytest.approx((0.5, 0.0, 0.5), abs=1e-12)
    assert face.surface.area(order=4) == pytest.approx(1.0, rel=1e-12)
    audit = audit_topology(brep)
    # One face on its own is an open shell: four free edges, watertightness 0.
    assert audit.watertight is False
    assert len(audit.free_edges) == 4
    assert {defect.code for defect in audit.defects} == {"FREE_EDGE"}


def test_inconsistent_loop_winding_is_reported(tmp_path: Path) -> None:
    """Negative control for the winding check: declare ``.F.``, keep the loop.

    The loop still winds anticlockwise about +n, so declaring the face reversed is an
    internal contradiction -- and that is precisely the "inside-out face" a kernel has to
    catch, because a viewer shows nothing.
    """
    flipped = STEP_FIXTURE.replace("#70=ADVANCED_FACE('fixture_face',(#50),#60,.T.);",
                                  "#70=ADVANCED_FACE('fixture_face',(#50),#60,.F.);")
    path = tmp_path / "flipped.step"
    path.write_text(flipped)
    brep = read_step(path)
    codes = {defect.code for defect in audit_topology(brep).defects}
    assert "ORIENTATION" in codes


@pytest.mark.parametrize("name", ["wigley_iii", "series60_cb0.60"])
def test_round_trip_preserves_geometry_and_topology(name: str, tmp_path: Path) -> None:
    specimen = baselines_module.specimens()[name]()
    results = round_trip(specimen.brep, tmp_path / name, stem=name)
    for suffix, result in results.items():
        assert result["watertight_before"] is True
        assert result["watertight_after"] is True, (suffix, result["codes_after"])
        # The polyline boundary curves are written as degree-1 splines, so the surfaces
        # come back to round-off and the topology comes back exactly.
        assert result["max_deviation"] < 1e-9, (suffix, result["max_deviation"])
        assert result["faces"] == len(specimen.brep.faces)
        assert result["edges"] == len(specimen.brep.edges)


def test_wetted_flag_survives_every_format(tmp_path: Path) -> None:
    specimen = baselines_module.wigley_iii()
    for suffix in (".brep.json", ".step", ".iges"):
        path = tmp_path / f"model{suffix}"
        save_cad(path, specimen.brep)
        reloaded = load_cad(path)
        flags = [face.wetted for face in reloaded.faces]
        assert any(flags) and not all(flags), (suffix, flags)
        # The waterplane stays excluded, so A_wp is the hull's, not zero.
        hydro = patch_hydrostatics(
            reloaded,
            length=specimen.length,
            beam=specimen.beam,
            draft=specimen.draft,
        )
        assert hydro.waterplane_area == pytest.approx(
            specimen.reference["waterplane_area"], rel=1e-9
        ), suffix
        assert hydro.volume == pytest.approx(specimen.reference["volume"], rel=1e-9)


def test_iges_reader_has_to_orient_the_shell(tmp_path: Path) -> None:
    """IGES carries no orientation; the reader must derive it or the volume is negative."""
    specimen = baselines_module.wigley_iii()
    path = tmp_path / "model.iges"
    save_cad(path, specimen.brep)
    reloaded = read_iges(path)
    audit = audit_topology(reloaded)
    assert audit.watertight is True
    assert audit.signed_volume > 0
    assert "ORIENTATION" not in {defect.code for defect in audit.defects}


def test_unsupported_inputs_raise(tmp_path: Path) -> None:
    with pytest.raises(ExchangeError):
        load_cad(tmp_path / "model.stl")
    # A type-144 (trimmed surface) record is out of scope, and must say so rather than
    # silently producing a wrong patch.
    specimen = baselines_module.wigley_iii()
    iges = tmp_path / "trimmed.iges"
    save_cad(iges, specimen.brep)
    iges.write_text(iges.read_text().replace("     128", "     144", 1))
    with pytest.raises(ExchangeError, match="144"):
        read_iges(iges)


def test_step_without_brep_entities_raises(tmp_path: Path) -> None:
    path = tmp_path / "facets.step"
    path.write_text(
        "ISO-10303-21;\nHEADER;\nENDSEC;\nDATA;\n#1=CARTESIAN_POINT('',(0.,0.,0.));\n"
        "ENDSEC;\nEND-ISO-10303-21;\n"
    )
    with pytest.raises(ExchangeError):
        read_step(path)


def test_mirror_welds_the_symmetry_plane() -> None:
    """A mirrored half model must be closed: that is what makes the 0/1 score meaningful."""
    specimen = baselines_module.wigley_iii()
    audit = audit_topology(specimen.brep)
    assert audit.watertight is True
    assert len(specimen.brep.faces) == 4
    # Two hull patches + two deck patches, welded along keel, stem, stern and waterline.
    assert len(specimen.brep.edges) == 8
