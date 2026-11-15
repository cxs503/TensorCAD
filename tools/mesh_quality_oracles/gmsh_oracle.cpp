/* Gmsh oracle for the mesh-quality library (TensorCAD, GAP-CAD-D1-MESH-QUALITY-LIB).
 *
 * Reads a mesh in the line format documented in
 * `tools/mesh_quality_oracles/README.md`, pushes it into a *discrete* Gmsh model and
 * asks Gmsh itself for the per-element quality measures through
 * `gmsh::model::mesh::getElementQualities` (the same call the Gmsh Python API exposes).
 * Nothing is recomputed here: every number printed is Gmsh's own.
 *
 * Build (see build_oracles.sh):
 *   g++ -O2 -std=c++11 gmsh_oracle.cpp -I<prefix>/usr/include \
 *       -L<prefix>/usr/lib/<triplet> -lgmsh -Wl,-rpath,<prefix>/usr/lib/<triplet> -o gmsh_oracle
 *
 * Usage:  gmsh_oracle <mesh.txt>
 * Output: one line per (measure, element):  V <measure> <index> <value>
 */
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <gmsh.h>

namespace {

struct Mesh {
  std::vector<double> coords;                  // 3*N
  // kind -> list of element node-index lists
  std::vector<std::pair<std::string, std::vector<std::vector<int>>>> blocks;
};

std::vector<std::string> split(const std::string &s) {
  std::vector<std::string> out;
  std::istringstream is(s);
  std::string tok;
  while (is >> tok) out.push_back(tok);
  return out;
}

bool read_mesh(const char *path, Mesh &m) {
  std::ifstream f(path);
  if (!f) return false;
  std::string line;
  while (std::getline(f, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::vector<std::string> tok = split(line);
    if (tok.empty()) continue;
    if (tok[0] == "nodes") {
      std::size_t n = std::strtoul(tok[1].c_str(), nullptr, 10);
      m.coords.reserve(3 * n);
      for (std::size_t i = 0; i < n; ++i) {
        if (!std::getline(f, line)) return false;
        std::vector<std::string> c = split(line);
        if (c.size() < 3) return false;
        for (int k = 0; k < 3; ++k) m.coords.push_back(std::strtod(c[k].c_str(), nullptr));
      }
    } else {
      // "<kind> <count>" followed by count lines of node indices
      std::size_t n = strtoul(tok[1].c_str(), nullptr, 10);
      std::vector<std::vector<int>> elems;
      elems.reserve(n);
      for (std::size_t i = 0; i < n; ++i) {
        if (!std::getline(f, line)) return false;
        std::vector<std::string> c = split(line);
        std::vector<int> e;
        e.reserve(c.size());
        for (std::size_t k = 0; k < c.size(); ++k) e.push_back(atoi(c[k].c_str()));
        elems.push_back(e);
      }
      m.blocks.push_back(std::make_pair(tok[0], elems));
    }
  }
  return true;
}

int gmsh_type(const std::string &kind, std::size_t nodes) {
  if (kind == "tri" && nodes == 3) return 2;
  if (kind == "quad" && nodes == 4) return 3;
  if (kind == "tet" && nodes == 4) return 4;
  if (kind == "hex" && nodes == 8) return 5;
  if (kind == "prism" && nodes == 6) return 6;
  if (kind == "pyramid" && nodes == 5) return 7;
  return -1;
}

int dim_of(const std::string &kind) {
  return (kind == "tri" || kind == "quad") ? 2 : 3;
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cerr << "usage: gmsh_oracle <mesh.txt>\n";
    return 2;
  }
  Mesh mesh;
  if (!read_mesh(argv[1], mesh)) {
    std::cerr << "gmsh_oracle: cannot read " << argv[1] << "\n";
    return 2;
  }

  gmsh::initialize(argc, argv, false);
  gmsh::option::setNumber("General.Verbosity", 0);
  gmsh::option::setNumber("General.Terminal", 0);
  gmsh::model::add("mq");

  // one discrete entity per dimension, tags 1 (surface) and 2 (volume)
  gmsh::model::addDiscreteEntity(2, 1);
  gmsh::model::addDiscreteEntity(3, 2);

  const std::size_t numNodes = mesh.coords.size() / 3;
  std::vector<std::size_t> nodeTags(numNodes);
  for (std::size_t i = 0; i < numNodes; ++i) nodeTags[i] = i + 1;
  // Gmsh needs the nodes of an entity; register them on every entity we create.
  gmsh::model::mesh::addNodes(2, 1, nodeTags, mesh.coords);
  gmsh::model::mesh::addNodes(3, 2, nodeTags, mesh.coords);

  std::vector<std::string> measures = {"minSICN",   "gamma",    "eta",      "volume",
                                       "minEdge",   "maxEdge",  "outerRadius",
                                       "innerRadius"};

  for (std::size_t b = 0; b < mesh.blocks.size(); ++b) {
    const std::string &kind = mesh.blocks[b].first;
    const std::vector<std::vector<int>> &elems = mesh.blocks[b].second;
    if (elems.empty()) continue;
    const std::size_t npe = elems.front().size();
    const int etype = gmsh_type(kind, npe);
    if (etype < 0) continue;
    const int dim = dim_of(kind);
    const int entity = (dim == 2) ? 1 : 2;

    std::vector<int> elementTypes(1, etype);
    std::vector<std::vector<std::size_t>> elementTags(1);
    std::vector<std::vector<std::size_t>> elementNodes(1);
    elementTags[0].reserve(elems.size());
    elementNodes[0].reserve(elems.size() * npe);
    for (std::size_t i = 0; i < elems.size(); ++i) {
      elementTags[0].push_back(i + 1);  // per-kind tags; only used inside this loop
      for (std::size_t k = 0; k < npe; ++k) elementNodes[0].push_back((std::size_t)elems[i][k] + 1);
    }
    gmsh::model::mesh::addElements(dim, entity, elementTypes, elementTags, elementNodes);
    std::vector<std::size_t> etags; for (std::size_t i = 0; i < elems.size(); ++i) etags.push_back(i + 1);
    elementTags[0] = etags;

    for (std::size_t mi = 0; mi < measures.size(); ++mi) {
      std::vector<double> values;
      try {
        gmsh::model::mesh::getElementQualities(elementTags[0], values, measures[mi]);
      } catch (const std::exception &e) {
        std::cerr << "gmsh_oracle: " << measures[mi] << " failed: " << e.what() << "\n";
        continue;
      }
      if (values.size() != elems.size()) {
        std::cerr << "gmsh_oracle: " << measures[mi] << " returned " << values.size()
                  << " values for " << elems.size() << " elements\n";
        continue;
      }
      for (std::size_t i = 0; i < values.size(); ++i) {
        std::printf("V %s %s %zu %.17g\n", kind.c_str(), measures[mi].c_str(), i, values[i]);
      }
    }

    // remove the elements again so the next block's tags start clean
    
  }

  std::printf("V __gmsh__ __version__ 0 %.17g\n", 0.0);
  gmsh::finalize();
  return 0;
}