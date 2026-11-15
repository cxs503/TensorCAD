// Gmsh *source* oracle: Gmsh's own quality measures, compiled from Gmsh's source tree.
//
// Why this exists
// ---------------
// `gmsh_oracle.cpp` calls `gmsh::model::mesh::getElementQualities` through
// `libgmsh.so`; on this image that library cannot be loaded (its OpenCASCADE /
// FLTK / METIS / CGNS closure is absent), so the Gmsh half of the anchor's
// differential stayed NOT RUN. The functions `getElementQualities` calls are
// `MElement::gammaShapeMeasure()` / `etaShapeMeasure()` / `minSICNShapeMeasure()`,
// which dispatch into `src/mesh/qualityMeasures.cpp` (`qmTriangle::gamma`,
// `qmTriangle::eta`, `qmQuadrangle::eta`, `qmTetrahedron::qm(..., QMTET_GAMMA|ETA)`).
//
// This harness compiles *those Gmsh sources verbatim* -- nothing is re-typed and
// no formula is transcribed -- together with the element classes they operate on
// (`MTriangle`, `MQuadrangle`, `MTetrahedron`, `MHexahedron`, `MElement`,
// `MVertex`, `MFace`) and the numeric helpers (`Numeric.cpp`,
// `robustPredicates.cpp`, `ElementType.cpp`, `GaussIntegration.cpp`,
// `GaussQuadrature*.cpp`). So the numbers below are Gmsh's arithmetic, executed,
// not Gmsh's arithmetic read and re-implemented.
//
// What is deliberately *not* Gmsh: the mesh file reader and the `main` loop
// (below), and two symbols that only the display/vtable paths reference
// (`Msg::Error`, `CTX::instance`) which are stubbed in `gmsh_source_shims.cpp`
// because Gmsh's own definitions drag the whole Options/ColorTable graph. The
// quality measures never call them.
//
// Protocol (identical to `verdict_oracle.cpp` / `gmsh_oracle.cpp`):
//
//     W <kind> <metric> <index> <value>
//
// one line per (element, measure), so `tools/mesh_quality_external_diff.py`
// reads all three oracles with the same parser.
//
// Build: tools/mesh_quality_oracles/build_gmsh_source_oracle.sh (see README).
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// Gmsh headers (verbatim, from the Gmsh source tree).
#include "MElement.h"
#include "MHexahedron.h"
#include "MQuadrangle.h"
#include "MTetrahedron.h"
#include "MTriangle.h"
#include "MVertex.h"

namespace {

std::vector<std::string> split(const std::string &s)
{
  std::vector<std::string> out;
  std::istringstream is(s);
  std::string tok;
  while (is >> tok) out.push_back(tok);
  return out;
}

// Mesh file reader: "mq-mesh v1" (see tools/mesh_quality_oracles/README.md).
struct Mesh {
  std::vector<double> xyz;
  std::vector<std::pair<std::string, std::vector<std::vector<int>>>> blocks;
};

bool read_mesh(const char *path, Mesh &m)
{
  std::ifstream f(path);
  if(!f) return false;
  std::string line;
  while(std::getline(f, line)) {
    if(line.empty() || line[0] == '#') continue;
    std::vector<std::string> tok = split(line);
    if(tok.empty()) continue;
    if(tok[0] == "nodes") {
      std::size_t n = std::strtoul(tok[1].c_str(), nullptr, 10);
      for(std::size_t i = 0; i < n; ++i) {
        if(!std::getline(f, line)) return false;
        std::vector<std::string> c = split(line);
        if(c.size() < 3) return false;
        for(int k = 0; k < 3; ++k)
          m.xyz.push_back(std::strtod(c[k].c_str(), nullptr));
      }
    }
    else {
      std::size_t n = std::strtoul(tok[1].c_str(), nullptr, 10);
      std::vector<std::vector<int>> elems;
      elems.reserve(n);
      for(std::size_t i = 0; i < n; ++i) {
        if(!std::getline(f, line)) return false;
        std::vector<std::string> c = split(line);
        std::vector<int> e;
        e.reserve(c.size());
        for(std::size_t k = 0; k < c.size(); ++k)
          e.push_back(std::atoi(c[k].c_str()));
        elems.push_back(e);
      }
      m.blocks.push_back(std::make_pair(tok[0], elems));
    }
  }
  return true;
}

void emit(const char *kind, const char *metric, int idx, double value)
{
  std::printf("W %s %s %d %.17g\n", kind, metric, idx, value);
}

} // namespace

int main(int argc, char **argv)
{
  if(argc < 2) {
    std::fprintf(stderr, "usage: %s <mq-mesh-file>\n", argv[0]);
    return 2;
  }
  Mesh m;
  if(!read_mesh(argv[1], m)) {
    std::fprintf(stderr, "cannot read mesh file %s\n", argv[1]);
    return 2;
  }
  const std::size_t nv = m.xyz.size() / 3;
  std::vector<MVertex *> verts(nv);
  for(std::size_t i = 0; i < nv; ++i)
    verts[i] = new MVertex(m.xyz[3 * i], m.xyz[3 * i + 1], m.xyz[3 * i + 2],
                           nullptr, i);

  for(std::size_t b = 0; b < m.blocks.size(); ++b) {
    const std::string &kind = m.blocks[b].first;
    const std::vector<std::vector<int>> &elems = m.blocks[b].second;
    for(std::size_t e = 0; e < elems.size(); ++e) {
      const std::vector<int> &n = elems[e];
      const int idx = (int)e;
      if(kind == "tri" && n.size() >= 3) {
        MTriangle el(verts[n[0]], verts[n[1]], verts[n[2]], idx, 0);
        emit("tri", "gamma", idx, el.gammaShapeMeasure());
        emit("tri", "eta", idx, el.etaShapeMeasure());
      }
      else if(kind == "quad" && n.size() >= 4) {
        MQuadrangle el(verts[n[0]], verts[n[1]], verts[n[2]], verts[n[3]], idx,
                       0);
        emit("quad", "eta", idx, el.etaShapeMeasure());
      }
      else if(kind == "tet" && n.size() >= 4) {
        MTetrahedron el(verts[n[0]], verts[n[1]], verts[n[2]], verts[n[3]], idx,
                        0);
        emit("tet", "gamma", idx, el.gammaShapeMeasure());
        emit("tet", "eta", idx, el.etaShapeMeasure());
        emit("tet", "volume", idx, el.getVolume());
      }
      else {
        std::fprintf(stderr, "unsupported block kind '%s' (%zu nodes)\n",
                     kind.c_str(), n.size());
        return 3;
      }
    }
  }
  for(std::size_t i = 0; i < nv; ++i) delete verts[i];
  return 0;
}