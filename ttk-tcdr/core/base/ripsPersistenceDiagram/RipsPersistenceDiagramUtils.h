/// \ingroup base
/// \author Mattéo Clémot <matteo.clemot@univ-lyon1.fr>
/// \date January 2024.

#pragma once

#include <vector>

#include <boost/dynamic_bitset.hpp>

namespace ttk::rpd {
  using id_t = int;
  using value_t = double;
  constexpr value_t inf = std::numeric_limits<value_t>::infinity();

  using Simplex = std::vector<id_t>;
  using FiltratedSimplex = std::pair<Simplex, value_t>;
  using PersistencePair= std::pair<FiltratedSimplex, FiltratedSimplex>;
  using Diagram = std::vector<PersistencePair>;
  using MultidimensionalDiagram = std::vector<Diagram>;

  using Edge = std::pair<id_t, id_t>;
  using EdgeSet = std::vector<Edge>;
  using EdgeSetSet = std::vector<EdgeSet>;
  using Cascade = EdgeSet;

  using Generator = std::pair<EdgeSet, std::pair<value_t,value_t>>;

  struct FiltratedEdge {
    std::pair<id_t, id_t> e;
    value_t d;
  };
  FiltratedEdge max(FiltratedEdge a, FiltratedEdge b);

  struct FiltratedQuadEdge {
    std::pair<id_t, id_t> e;
    int f1;
    int f2;
    value_t d;
  };

  struct FiltratedTriangle {
    std::tuple<id_t, id_t, id_t> t;
    value_t d;
  };

  class UnionFind {
  private:
    std::vector<int> parent, rank;
  public:
    explicit UnionFind(unsigned n);
    int find(int x);
    void merge(int x, int y);
    int merge_ret(int x, int y);
    [[nodiscard]] bool isRoot(int x) const;
  };

  class BoundaryContainer {
  public:
    explicit BoundaryContainer(std::vector<id_t> &simplices, unsigned size)
      : ids(simplices) {
      mask.resize(size, false);
      for(id_t const &id : ids)
        mask[id] = true;
    }
    void exclusiveAddBoundary(std::vector<id_t> const &boundary) {
      for(id_t const &id : boundary) {
        if(!mask[id])
          ids.emplace_back(id);
        else
          ids.erase(std::find(ids.begin(), ids.end(), id));
        mask.flip(id);
      }
    }

  private:
    std::vector<id_t> &ids;
    boost::dynamic_bitset<> mask;
  };

}