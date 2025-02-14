/// \ingroup base
/// \author Mattéo Clémot <matteo.clemot@univ-lyon1.fr>
/// \date February 2025.
///
/// \brief
///
///
///

#pragma once

#include <Debug.h>
#include <RipsPersistenceDiagramUtils.h>

namespace ttk::rpd {

  class PairCellsWithOracle : virtual public Debug {
  public:
    PairCellsWithOracle(const std::vector<std::vector<double>> &points,
                        rpd::MultidimensionalDiagram const& oracle,
                        bool distanceMatrix=false,
                        bool parallelSort=false);

    void run();

    void getGenerators(std::vector<Generator> &generators) const;

    void getCascades(std::vector<Cascade> &cascades, std::vector<EdgeSet> &critical) const;
    void getCascades(EdgeSet &cascades, std::vector<EdgeSet> &critical) const;

    static void callOracle(const std::vector<std::vector<double>> &points,
                           rpd::MultidimensionalDiagram &oracle,
                           double threshold = inf);

  private:
    int N{};
    std::vector<double> compressedDM {};
    const bool parallelSort_;

    rpd::MultidimensionalDiagram const& oracle_;
    double bound {0.};

    std::map<Edge, id_t> edgeToIndex {};
    std::vector<std::set<id_t, std::greater<id_t>>> graph {}; //std::greater to match the reverse co-lexicographic order used in Ripser
    std::vector<FiltratedEdge> edges {};
    std::vector<FiltratedTriangle> triangles {};
    std::vector<id_t> edgesIndices {};
    std::vector<id_t> edgesOrder {};

    std::vector<id_t> edgesPartner {};
    std::vector<id_t> trianglesPartner {};
    std::vector<std::vector<id_t>> boundaries {};
    int n_edges {};

    std::vector<std::vector<id_t>> cascadeEdges {};

    /**
     * distance matrix access function
     */
    inline double& DM(unsigned i, unsigned j) {
      return compressedDM[std::max(i,j) * (std::max(i,j) - 1) / 2 + std::min(i,j)];
    }

    const std::function<bool(id_t,id_t)> compEdges = [this](id_t i, id_t j) -> bool {
      if (edges[i].d == edges[j].d)
        return edges[i].e < edges[j].e;
      else
        return edges[i].d < edges[j].d;
    };

    const std::function<bool(id_t,id_t)> quickCompEdges = [this](id_t i, id_t j) -> bool {
      return edgesOrder[i] < edgesOrder[j];
    };

    void initializeWithBound();
    void pairCellsWithOracle();
    void eliminateBoundaryWithOracle(id_t t_id, id_t e_id);

    void fillRNG(std::vector<EdgeSet> &critical) const;
  };
}