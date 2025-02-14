#include <PairCells.h>

#ifdef TTK_ENABLE_CGAL
ttk::rpd::PairCells::PairCells(const std::vector<CGAL::Epick::Point_2> &points, double upperBound, bool parallelSort, bool parallelMatrixConstruction) : bound(upperBound), parallelSort_(parallelSort), parallelMatrixConstruction_(parallelMatrixConstruction) {
  // inherited from Debug: prefix will be printed at the beginning of every msg
  this->setDebugMsgPrefix("PairCells");

  N = points.size();
  for (int i = 1; i < N; ++i) {
    for (int j = 0; j < i; ++j) {
      compressedDM.push_back(sqrt(CGAL::squared_distance(points[i], points[j])));
    }
  }
}
#endif

ttk::rpd::PairCells::PairCells(const std::vector<std::vector<double>> &points, bool distanceMatrix, double upperBound, bool parallelSort, bool parallelMatrixConstruction) : bound(upperBound), parallelSort_(parallelSort), parallelMatrixConstruction_(parallelMatrixConstruction) {
  // inherited from Debug: prefix will be printed at the beginning of every msg
  this->setDebugMsgPrefix("PairCells");

  N = points.size();
  if (distanceMatrix)
    compressedDM = points[0];
  else {
    const unsigned dim = points[0].size();
    for (int i = 1; i < N; ++i) {
      for (int j = 0; j < i; ++j) {
        double s = 0.;
        for (unsigned d = 0; d < dim; ++d)
          s += (points[i][d] - points[j][d]) * (points[i][d] - points[j][d]);
        compressedDM.push_back(sqrt(s));
      }
    }
  }
}

ttk::rpd::PairCells::PairCells(float* data, int n, int dim, double upperBound, bool parallelSort, bool parallelMatrixConstruction) : bound(upperBound), parallelSort_(parallelSort), parallelMatrixConstruction_(parallelMatrixConstruction) {
  // inherited from Debug: prefix will be printed at the beginning of every msg
  this->setDebugMsgPrefix("PairCells");

  N = n;
  for (int i = 1; i < N; ++i) {
    for (int j = 0; j < i; ++j) {
      double s = 0.;
      for (int d = 0; d < dim; ++d)
        s += (data[dim*i+d] - data[dim*j+d]) * (data[dim*i+d] - data[dim*j+d]);
      compressedDM.push_back(sqrt(s));
    }
  }
}

void ttk::rpd::PairCells::run() {
  Timer tm{};
  if (bound == inf || bound < 0)
    initialize();
  else
    initializeWithBound();
  printMsg("Initialized (#p=" + std::to_string(N) + ", #e=" + std::to_string(n_edges) + ", #t=" + std::to_string(n_triangles) + ")", 0, tm.getElapsedTime());
  pairCells();
  printMsg("Complete", 1, tm.getElapsedTime());
}

void ttk::rpd::PairCells::initialize() {
  //edges
  for (id_t i=0; i < N; ++i) {
    for (id_t j=i+1; j < N; ++j)
      edges.push_back({{i, j}, DM(i,j)});
  }
  n_edges = N*(N-1) / 2;
  executeKruskal();

  //triangles
  n_triangles = N*(N-1)*(N-2) / 6;
  triangles.resize(n_triangles);
  boundaries.resize(n_triangles);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for if(parallelMatrixConstruction_)
#endif // TTK_ENABLE_OPENMP
  for (id_t i=0; i < N-2; ++i) {
    const unsigned index_i = i * (i*i - 3*i*(N-1) + 3*N*N - 6*N + 2) / 6;
    for (id_t j=i+1; j < N-1; ++j) {
      const unsigned index_j = (i-j+1) * (i+j-2*N+2) / 2 + index_i;
      for (id_t k=j+1; k < N; ++k) {
        const double diameter = std::max(DM(i,j), std::max(DM(i,k), DM(j,k)));
        triangles[index_j + (k-j-1)] = {{i,j,k}, diameter};
        boundaries[index_j + (k-j-1)] = {n_edges-1 - (N-i)*(N-i-1)/2 + (j-i),
                                         n_edges-1 - (N-i)*(N-i-1)/2 + (k-i),
                                         n_edges-1 - (N-j)*(N-j-1)/2 + (k-j)};
      }
    }
  }

  apparentPairs();
}

void ttk::rpd::PairCells::initializeWithBound() {
  //edges
  std::map<Edge, id_t> edgeToIndex;
  std::vector<std::set<id_t>> graph (N); //not ok with std::unordered_set, why ?
  for (id_t i=0; i < N; ++i) {
    for (id_t j=i+1; j < N; ++j) {
      if (DM(i,j) <= bound) {
        edgeToIndex[{i,j}] = edges.size();
        edges.push_back({{i, j}, DM(i,j)});
        graph[i].insert(j); //only directed graph with edges (ij), i<j
      }
    }
  }
  n_edges = edges.size();
  executeKruskal();

  //triangles
  for (const auto& [e, f] : edges) {
    const auto& [i,j] = e;
    for (auto const& k : graph[j]) {
      if (graph[i].find(k) != graph[i].end()) { // c++20 -> std::unordered_set::contains
        triangles.push_back({{i,j,k}, std::max(DM(i,j), std::max(DM(i,k), DM(j,k)))});
        boundaries.push_back({edgeToIndex[{i, j}],
                              edgeToIndex[{i, k}],
                              edgeToIndex[{j, k}]});
      }
    }
  }
  n_triangles = triangles.size();

  apparentPairs();
}

void ttk::rpd::PairCells::executeKruskal() {
  UnionFind UF(N);
  edgesIndices.resize(n_edges);
  edgesOrder.resize(n_edges);
  edgesPartner.resize(n_edges, -1);
  std::iota(edgesIndices.begin(), edgesIndices.end(), 0);
  if (parallelSort_)
    TTK_PSORT(globalThreadNumber_, edgesIndices.begin(), edgesIndices.end(), compEdges)
  else
    std::sort(edgesIndices.begin(), edgesIndices.end(), compEdges);
  symbolicPerturbation();
  for (unsigned i=0; i<edgesIndices.size(); ++i)
    edgesOrder[edgesIndices[i]] = i;
  for (const id_t &id : edgesIndices) {
    const auto& [s,f] = edges[id];
    if (UF.find(s.first) != UF.find(s.second)) {
      UF.merge(s.first, s.second);
      edgesPartner[id] = -2;
      if (++n_pairedEdges == N-1)
        break;
    }
  }
}

void ttk::rpd::PairCells::symbolicPerturbation(double eps) { //todo to be improved (many calls when lots of collisions)
  bool generic;
  do {
    generic = true;
    for (unsigned i=0; i+1 < edgesIndices.size(); ++i) {
      if (edges[edgesIndices[i]].d >= edges[edgesIndices[i+1]].d) {
        generic = false;
        FiltratedEdge &e = edges[edgesIndices[i+1]];
        e.d += eps;
        DM(e.e.first, e.e.second) += eps;
      }
    }
    if (!generic) {
      if (parallelSort_)
        TTK_PSORT(globalThreadNumber_, edgesIndices.begin(), edgesIndices.end(), compEdges)
      else
        std::sort(edgesIndices.begin(), edgesIndices.end(), compEdges);
    }
  } while (!generic);
}

void ttk::rpd::PairCells::apparentPairs() {
  //arrays initialisations
  trianglesIndices.resize(n_triangles);
  trianglesPartner.resize(n_triangles, -1);
  cascadeEdges.resize(n_triangles);

  //full order
  std::iota(trianglesIndices.begin(), trianglesIndices.end(), 0);
  if (parallelSort_)
    TTK_PSORT(globalThreadNumber_, trianglesIndices.begin(), trianglesIndices.end(), compTriangles)
  else
    std::sort(trianglesIndices.begin(), trianglesIndices.end(), compTriangles);

  //apparent pairs
  id_t t = 0;
  for(const id_t &e : edgesIndices) {
    const double f = edges[e].d;
    while (triangles[trianglesIndices[t]].d < f)
      ++t;
    if (triangles[trianglesIndices[t]].d == f) {
      edgesPartner[e] = trianglesIndices[t];
      trianglesPartner[trianglesIndices[t]] = e;
      n_pairedEdges++;
    }
  }
}

void ttk::rpd::PairCells::pairCells() {
  for (const id_t &s : trianglesIndices) {
    if (trianglesPartner[s] == -1) {
      const int e = eliminateBoundaries(s);
      if (e != -1) {
        trianglesPartner[s] = e;
        edgesPartner[e] = s;
        if (++n_pairedEdges == n_edges)
          break;
      }
    }
  }
}

int ttk::rpd::PairCells::eliminateBoundaries(id_t s) {
  BoundaryContainer bc (boundaries[s], n_edges);
  while (!boundaries[s].empty()) {
    const id_t e = *std::max_element(boundaries[s].begin(), boundaries[s].end(), quickCompEdges);
    cascadeEdges[s].push_back(e);
    if (edgesPartner[e] == -1)
      return int(e);
    else
      bc.exclusiveAddBoundary(boundaries[edgesPartner[e]]);
  }
  return -1;
}

void ttk::rpd::PairCells::getDiagram(MultidimensionalDiagram &diagrams) const {
  diagrams.resize(2);
  for (const id_t &e : edgesIndices) {
    if (edgesPartner[e] == -2)
      diagrams[0].emplace_back(FiltratedSimplex{{0}, 0.}, FiltratedSimplex{{(int)edges[e].e.first, (int)edges[e].e.second}, edges[e].d});
    else if (edgesPartner[e] > 0 && edges[e].d < triangles[edgesPartner[e]].d) {
      const FiltratedSimplex birth ({(int)edges[e].e.first, (int)edges[e].e.second}, edges[e].d);
      const FiltratedSimplex death ({(int)std::get<0>(triangles[edgesPartner[e]].t),
                                  (int)std::get<1>(triangles[edgesPartner[e]].t),
                                  (int)std::get<2>(triangles[edgesPartner[e]].t)}, triangles[edgesPartner[e]].d);
      diagrams[1].emplace_back(birth, death);
    }
  }
  diagrams[0].emplace_back(FiltratedSimplex{{0}, 0.}, FiltratedSimplex{{-1}, inf});
}

void ttk::rpd::PairCells::getDiagramAndGenerators(MultidimensionalDiagram &diagrams, std::vector<Generator> &generators) const {
  diagrams.resize(2);
  for (const id_t &e : edgesIndices) {
    if (edgesPartner[e] == -2)
      diagrams[0].emplace_back(FiltratedSimplex{{0}, 0.}, FiltratedSimplex{{(int)edges[e].e.first, (int)edges[e].e.second}, edges[e].d});
    else if (edgesPartner[e] > 0 && edges[e].d < triangles[edgesPartner[e]].d) {
      const FiltratedSimplex birth ({(int)edges[e].e.first, (int)edges[e].e.second}, edges[e].d);
      const FiltratedSimplex death ({(int)std::get<0>(triangles[edgesPartner[e]].t),
                                     (int)std::get<1>(triangles[edgesPartner[e]].t),
                                     (int)std::get<2>(triangles[edgesPartner[e]].t)}, triangles[edgesPartner[e]].d);
      diagrams[1].emplace_back(birth, death);

      EdgeSet boundary;
      for (const id_t &edge : boundaries[edgesPartner[e]])
        boundary.emplace_back(edges[edge].e);
      generators.emplace_back(boundary, std::make_pair(edges[e].d, triangles[edgesPartner[e]].d));
    }
  }
  diagrams[0].emplace_back(FiltratedSimplex{{0}, 0.}, FiltratedSimplex{{-1}, inf});
}

void ttk::rpd::PairCells::getCascades(std::vector<Cascade> &cascades, std::vector<EdgeSet> &critical) const {
  critical = std::vector<EdgeSet>(3);
  for (const id_t &e : edgesIndices) {
    if (edgesPartner[e] == -2) //MST
      critical[0].emplace_back(edges[e].e);
    else if (edgesPartner[e] > 0 && edges[e].d < triangles[edgesPartner[e]].d) {
      const Edge &killerEdge= edges[cascadeEdges[edgesPartner[e]][0]].e;
      critical[1].emplace_back(edges[e].e); //RNG edge
      critical[2].emplace_back(killerEdge); //MML edge
      Cascade cascade = {killerEdge};
      for (unsigned i=1; i<cascadeEdges[edgesPartner[e]].size()-1; ++i)
        cascade.emplace_back(edges[cascadeEdges[edgesPartner[e]][i]].e);
      cascades.emplace_back(cascade);
    }
  }
}

void ttk::rpd::PairCells::getCascades(EdgeSet &cascades, std::vector<EdgeSet> &critical) const {
  critical = std::vector<EdgeSet>(3);
  std::set<id_t> cascadeSet;
  for (const id_t &e : edgesIndices) {
    if (edgesPartner[e] == -2) //MST
      critical[0].emplace_back(edges[e].e);
    else if (edgesPartner[e] > 0 && edges[e].d < triangles[edgesPartner[e]].d) {
      const Edge &killerEdge= edges[cascadeEdges[edgesPartner[e]][0]].e;
      critical[1].emplace_back(edges[e].e); //RNG edge
      critical[2].emplace_back(killerEdge); //MML edge
      for (unsigned i=1; i<cascadeEdges[edgesPartner[e]].size()-1; ++i)
        cascadeSet.insert(cascadeEdges[edgesPartner[e]][i]);
    }
  }
  for (const id_t &e : cascadeSet)
    cascades.emplace_back(edges[e].e);
}

void ttk::rpd::PairCells::enrichCascades(std::set<Edge> &cascadeSet, std::vector<EdgeSet> &critical, std::vector<int> const& globalIndices) const {
  for (const id_t &e : edgesIndices) {
    if (edgesPartner[e] > 0 && edges[e].d < triangles[edgesPartner[e]].d) {
      const Edge &killerEdge = edges[cascadeEdges[edgesPartner[e]][0]].e;
      critical[2].emplace_back(globalIndices[killerEdge.first], globalIndices[killerEdge.second]); //MML edge
      for (unsigned i=1; i<cascadeEdges[edgesPartner[e]].size()-1; ++i) {
        const Edge &edge = edges[cascadeEdges[edgesPartner[e]][i]].e;
        cascadeSet.emplace(globalIndices[edge.first], globalIndices[edge.second]);
      }
    }
  }
}