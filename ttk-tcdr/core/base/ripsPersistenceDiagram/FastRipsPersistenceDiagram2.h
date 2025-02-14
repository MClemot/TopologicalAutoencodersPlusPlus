/// \ingroup base
/// \author Mattéo Clémot <matteo.clemot@univ-lyon1.fr>
/// \date January 2024.

#pragma once

#ifdef TTK_ENABLE_CGAL

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_id_2.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_sphere.h>

#include <PairCells.h>
#include <RipsPersistenceDiagramUtils.h>

namespace ttk::rpd {

  class FastRipsPersistenceDiagram2 : virtual public Debug {
    using K = CGAL::Exact_predicates_inexact_constructions_kernel;
    using Vb = CGAL::Triangulation_vertex_base_with_info_2<unsigned int, K>;
    using Fb = CGAL::Triangulation_face_base_with_id_2<K>;
    using Tds = CGAL::Triangulation_data_structure_2<Vb, Fb>;
    using Delaunay = CGAL::Delaunay_triangulation_2<K, Tds>;
    using Point = Delaunay::Point_2;
    using Traits = CGAL::Search_traits_2<K>;
    using Tree = CGAL::Kd_tree<Traits>;
    using Fuzzy_sphere = CGAL::Fuzzy_sphere<Traits>;

  public:
    explicit FastRipsPersistenceDiagram2(const std::vector<std::vector<double>> &_points);
    explicit FastRipsPersistenceDiagram2(float* data, int n);

    template <typename T>
    void compute0Persistence(T &ph0, bool parallelSort = false);

    template <typename T>
    void computeDelaunayRips0And1Persistence(T &ph, bool parallelSort = false);

    template <typename T>
    void computeRips0And1Persistence(T &ph, bool parallelSort = false, bool parallelMML = false);

    void computeRipsCascades(EdgeSet &cascades, std::vector<EdgeSet> &criticalEdges, bool parallelSort = false, bool parallelMML = false);

    void exportRips1Generators(std::vector<Generator> &generators);

  private:
    Timer tm {};
    const unsigned N_p;
    unsigned N_f;
    std::vector<std::pair<Point,unsigned>> points;
    Delaunay del;

    std::vector<FiltratedQuadEdge> urquhart;
    std::vector<FiltratedQuadEdge> rng;
    std::vector<FiltratedEdge> death_poly;
    std::vector<double> birth_poly;

    //common to Delaunay-Rips and Rips
    void computeDelaunay();
    void computeUrquhart(UnionFind& UF, std::vector<FiltratedEdge>& max_delaunay, bool parallelSort);
    void compute1PH(std::vector<FiltratedQuadEdge> const& critical, UnionFind &UF, MultidimensionalDiagram &ph);

    //specific to Rips
    [[nodiscard]] bool isLensEmpty(Point const& p1, Point const& p2, Tree const& tree, double const& d) const;
    [[nodiscard]] bool isRightSemiLensEmpty(Point const& p1, Point const& p2, Tree const& tree) const;
    void reindexPolygons(UnionFind const& UF, std::vector<FiltratedEdge> const& max_delaunay, std::vector<int>& index_polys);
    void computePolygonRipsDeath(bool parallel, UnionFind &UF, std::vector<int> const& index_polys);
    void pComputePolygonRipsDeath(UnionFind &UF, std::vector<int> const& index_polys);

    void executePolygonPairCells(bool parallel, UnionFind &UF, std::vector<int> const& index_polys, EdgeSet &cascades, std::vector<EdgeSet> &critical);
  };

}

#endif