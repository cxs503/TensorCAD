/* Verdict oracle for the mesh-quality library (TensorCAD, GAP-CAD-D1-MESH-QUALITY-LIB).
 *
 * Verdict is the Sandia quality library whose definitions VTK's `vtkMeshQuality` calls,
 * and `meshmagick.print_quality()` delegates to `vtk.vtkMeshQuality`. Running Verdict
 * itself is therefore the same arithmetic meshmagick reports, without needing VTK or a
 * GPL runtime on this host.
 *
 * Reads the mesh line format documented in `tools/mesh_quality_oracles/README.md` and
 * prints one line per (measure, element):
 *
 *     W <kind> <measure> <index> <value>
 *
 * Build (see build_oracles.sh): compile this file together with the Verdict sources
 * (V_TriMetric.cpp, V_QuadMetric.cpp, V_TetMetric.cpp, V_HexMetric.cpp,
 * VerdictVector.cpp, V_KnifeMetric.cpp, V_WedgeMetric.cpp, V_PyramidMetric.cpp,
 * V_GaussIntegration.cpp) and `-I<verdict-master>`.
 */
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "verdict.h"

namespace {

struct Mesh {
  std::vector<double> coords;
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
      for (std::size_t i = 0; i < n; ++i) {
        if (!std::getline(f, line)) return false;
        std::vector<std::string> c = split(line);
        if (c.size() < 3) return false;
        for (int k = 0; k < 3; ++k) m.coords.push_back(std::strtod(c[k].c_str(), nullptr));
      }
    } else {
      std::size_t n = std::strtoul(tok[1].c_str(), nullptr, 10);
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

typedef double (*Metric)(int, const double[][3]);

double tri_area(int n, const double c[][3]) { return verdict::tri_area(n, c); }
double tri_aspect_ratio(int n, const double c[][3]) { return verdict::tri_aspect_ratio(n, c); }
double tri_aspect_frobenius(int n, const double c[][3]) { return verdict::tri_aspect_frobenius(n, c); }
double tri_condition(int n, const double c[][3]) { return verdict::tri_condition(n, c); }
double tri_edge_ratio(int n, const double c[][3]) { return verdict::tri_edge_ratio(n, c); }
double tri_maximum_angle(int n, const double c[][3]) { return verdict::tri_maximum_angle(n, c); }
double tri_minimum_angle(int n, const double c[][3]) { return verdict::tri_minimum_angle(n, c); }
double tri_radius_ratio(int n, const double c[][3]) { return verdict::tri_radius_ratio(n, c); }
double tri_scaled_jacobian(int n, const double c[][3]) { return verdict::tri_scaled_jacobian(n, c); }
double tri_equiangle_skew(int n, const double c[][3]) { return verdict::tri_equiangle_skew(n, c); }

double quad_aspect_ratio(int n, const double c[][3]) { return verdict::quad_aspect_ratio(n, c); }
double quad_condition(int n, const double c[][3]) { return verdict::quad_condition(n, c); }
double quad_edge_ratio(int n, const double c[][3]) { return verdict::quad_edge_ratio(n, c); }
double quad_maximum_angle(int n, const double c[][3]) { return verdict::quad_maximum_angle(n, c); }
double quad_minimum_angle(int n, const double c[][3]) { return verdict::quad_minimum_angle(n, c); }
double quad_radius_ratio(int n, const double c[][3]) { return verdict::quad_radius_ratio(n, c); }
double quad_scaled_jacobian(int n, const double c[][3]) { return verdict::quad_scaled_jacobian(n, c); }
double quad_equiangle_skew(int n, const double c[][3]) { return verdict::quad_equiangle_skew(n, c); }

double tet_aspect_ratio(int n, const double c[][3]) { return verdict::tet_aspect_ratio(n, c); }
double tet_aspect_gamma(int n, const double c[][3]) { return verdict::tet_aspect_gamma(n, c); }
double tet_condition(int n, const double c[][3]) { return verdict::tet_condition(n, c); }
double tet_edge_ratio(int n, const double c[][3]) { return verdict::tet_edge_ratio(n, c); }
double tet_minimum_angle(int n, const double c[][3]) { return verdict::tet_minimum_angle(n, c); }
double tet_radius_ratio(int n, const double c[][3]) { return verdict::tet_radius_ratio(n, c); }
double tet_scaled_jacobian(int n, const double c[][3]) { return verdict::tet_scaled_jacobian(n, c); }
double tet_volume(int n, const double c[][3]) { return verdict::tet_volume(n, c); }
double tet_equiangle_skew(int n, const double c[][3]) { return verdict::tet_equiangle_skew(n, c); }

double hex_condition(int n, const double c[][3]) { return verdict::hex_condition(n, c); }
double hex_jacobian(int n, const double c[][3]) { return verdict::hex_jacobian(n, c); }
double hex_scaled_jacobian(int n, const double c[][3]) { return verdict::hex_scaled_jacobian(n, c); }
double hex_equiangle_skew(int n, const double c[][3]) { return verdict::hex_equiangle_skew(n, c); }
double hex_distortion(int n, const double c[][3]) { return verdict::hex_distortion(n, c); }
double hex_max_aspect_frobenius(int n, const double c[][3]) { return verdict::hex_max_aspect_frobenius(n, c); }

double prism_volume(int n, const double c[][3]) { return verdict::wedge_volume(n, c); }
double prism_scaled_jacobian(int n, const double c[][3]) { return verdict::wedge_scaled_jacobian(n, c); }
double prism_condition(int n, const double c[][3]) { return verdict::wedge_condition(n, c); }
double prism_equiangle_skew(int n, const double c[][3]) { return verdict::wedge_equiangle_skew(n, c); }

double pyramid_volume(int n, const double c[][3]) { return verdict::pyramid_volume(n, c); }
double pyramid_scaled_jacobian(int n, const double c[][3]) { return verdict::pyramid_scaled_jacobian(n, c); }
double pyramid_equiangle_skew(int n, const double c[][3]) { return verdict::pyramid_equiangle_skew(n, c); }

struct Entry {
  const char *kind;
  const char *name;
  Metric fn;
};

const Entry ENTRIES[] = {
    {"tri", "area", tri_area},
    {"tri", "aspect_ratio", tri_aspect_ratio},
    {"tri", "aspect_frobenius", tri_aspect_frobenius},
    {"tri", "condition", tri_condition},
    {"tri", "edge_ratio", tri_edge_ratio},
    {"tri", "max_angle", tri_maximum_angle},
    {"tri", "min_angle", tri_minimum_angle},
    {"tri", "radius_ratio", tri_radius_ratio},
    {"tri", "scaled_jacobian", tri_scaled_jacobian},
    {"tri", "equiangle_skew", tri_equiangle_skew},
    {"quad", "aspect_ratio", quad_aspect_ratio},
    {"quad", "condition", quad_condition},
    {"quad", "edge_ratio", quad_edge_ratio},
    {"quad", "max_angle", quad_maximum_angle},
    {"quad", "min_angle", quad_minimum_angle},
    {"quad", "radius_ratio", quad_radius_ratio},
    {"quad", "scaled_jacobian", quad_scaled_jacobian},
    {"quad", "equiangle_skew", quad_equiangle_skew},
    {"tet", "aspect_ratio", tet_aspect_ratio},
    {"tet", "aspect_gamma", tet_aspect_gamma},
    {"tet", "condition", tet_condition},
    {"tet", "edge_ratio", tet_edge_ratio},
    {"tet", "minimum_angle", tet_minimum_angle},
    {"tet", "radius_ratio", tet_radius_ratio},
    {"tet", "scaled_jacobian", tet_scaled_jacobian},
    {"tet", "volume", tet_volume},
    {"tet", "equiangle_skew", tet_equiangle_skew},
    {"hex", "condition", hex_condition},
    {"hex", "jacobian", hex_jacobian},
    {"hex", "scaled_jacobian", hex_scaled_jacobian},
    {"hex", "equiangle_skew", hex_equiangle_skew},
    {"hex", "distortion", hex_distortion},
    {"hex", "max_aspect_frobenius", hex_max_aspect_frobenius},
    {"prism", "volume", prism_volume},
    {"prism", "scaled_jacobian", prism_scaled_jacobian},
    {"prism", "condition", prism_condition},
    {"prism", "equiangle_skew", prism_equiangle_skew},
    {"pyramid", "volume", pyramid_volume},
    {"pyramid", "scaled_jacobian", pyramid_scaled_jacobian},
    {"pyramid", "equiangle_skew", pyramid_equiangle_skew},
};

}  // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cerr << "usage: verdict_oracle <mesh.txt>\n";
    return 2;
  }
  Mesh mesh;
  if (!read_mesh(argv[1], mesh)) {
    std::cerr << "verdict_oracle: cannot read " << argv[1] << "\n";
    return 2;
  }

  for (std::size_t b = 0; b < mesh.blocks.size(); ++b) {
    const std::string &kind = mesh.blocks[b].first;
    const std::vector<std::vector<int>> &elems = mesh.blocks[b].second;
    for (std::size_t e = 0; e < sizeof(ENTRIES) / sizeof(ENTRIES[0]); ++e) {
      if (kind != ENTRIES[e].kind) continue;
      const std::size_t npe = elems.empty() ? 0 : elems.front().size();
      for (std::size_t i = 0; i < elems.size(); ++i) {
        std::vector<double> c(3 * npe);
        for (std::size_t k = 0; k < npe; ++k) {
          const std::size_t n = (std::size_t)elems[i][k];
          c[3 * k + 0] = mesh.coords[3 * n + 0];
          c[3 * k + 1] = mesh.coords[3 * n + 1];
          c[3 * k + 2] = mesh.coords[3 * n + 2];
        }
        const double v = ENTRIES[e].fn((int)npe, (const double(*)[3])c.data());
        std::printf("W %s %s %zu %.17g\n", kind.c_str(), ENTRIES[e].name, i, v);
      }
    }
  }
  return 0;
}