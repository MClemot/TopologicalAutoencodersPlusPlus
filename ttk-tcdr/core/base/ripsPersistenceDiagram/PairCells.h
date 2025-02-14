/// \ingroup base
/// \author Mattéo Clémot <matteo.clemot@univ-lyon1.fr>
/// \date September 2024.
///
/// \brief TTK base class that executes the PairCells persistence algorithm
/// on a Rips complex.
///
/// This module defines the %PairCells class that takes a point cloud or a
/// distance matrix and can compute the persistence diagram, the persistent
/// generators and the persistent cascades of its Rips complex.
///

#pragma once

#include <Debug.h>
#include <RipsPersistenceDiagramUtils.h>

#ifdef TTK_ENABLE_CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#endif

namespace ttk::rpd {

  class PairCells : virtual public Debug {
  public:
#ifdef TTK_ENABLE_CGAL
    PairCells(const std::vector<CGAL::Epick::Point_2> &points,
              double upperBound=inf,
              bool parallelSort=false,
              bool parallelMatrixConstruction=false);
#endif
    PairCells(const std::vector<std::vector<double>> &points,
              bool distanceMatrix=false,
              double upperBound=inf,
              bool parallelSort=false,
              bool parallelMatrixConstruction=false);
    PairCells(float* data, int n, int dim,
              double upperBound=inf,
              bool parallelSort=false,
              bool parallelMatrixConstruction=false);

    void run();

    void getDiagram(MultidimensionalDiagram &diagrams) const;
    void getDiagramAndGenerators(MultidimensionalDiagram &diagrams, std::vector<Generator> &generators) const;

    void getCascades(std::vector<Cascade> &cascades, std::vector<EdgeSet> &critical) const;
    void getCascades(EdgeSet &cascades, std::vector<EdgeSet> &critical) const;
    void enrichCascades(std::set<Edge> &cascadeSet, std::vector<EdgeSet> &critical, std::vector<int> const& globalIndices) const;

  private:
    int N{};
    std::vector<double> compressedDM {};
    const double bound;
    const bool parallelSort_;
    const bool parallelMatrixConstruction_;

    std::vector<FiltratedEdge> edges {};
    std::vector<FiltratedTriangle> triangles {};
    std::vector<id_t> edgesIndices {};
    std::vector<id_t> edgesOrder {};
    std::vector<id_t> trianglesIndices {};

    std::vector<id_t> edgesPartner {};
    std::vector<id_t> trianglesPartner {};
    std::vector<std::vector<id_t>> boundaries {};
    int n_edges {};
    int n_pairedEdges {0};
    int n_triangles {};

    std::vector<std::vector<id_t>> cascadeEdges {};

    /**
     * distance matrix access function
     * \pre i < j
     */
    inline double& DM(unsigned i, unsigned j) {
      return compressedDM[j * (j - 1) / 2 + i];
    }

    const std::function<bool(id_t,id_t)> compTriangles = [this](id_t i, id_t j) -> bool {
      if (triangles[i].d == triangles[j].d)
        return triangles[i].t < triangles[j].t;
      else
        return triangles[i].d < triangles[j].d;
    };

    const std::function<bool(id_t,id_t)> compEdges = [this](id_t i, id_t j) -> bool {
      if (edges[i].d == edges[j].d)
        return edges[i].e < edges[j].e;
      else
        return edges[i].d < edges[j].d;
    };

    const std::function<bool(id_t,id_t)> quickCompEdges = [this](id_t i, id_t j) -> bool {
      return edgesOrder[i] < edgesOrder[j];
    };

    void initialize();
    void initializeWithBound();

    void executeKruskal();
    void symbolicPerturbation(double eps = std::numeric_limits<float>::epsilon());
    void apparentPairs();

    void pairCells();
    int eliminateBoundaries(id_t s);
  };
}