#include <PairCellsWithOracle.h>
#include <ripser.h>

ttk::rpd::PairCellsWithOracle::PairCellsWithOracle(const std::vector<std::vector<double>> &points, rpd::MultidimensionalDiagram const& oracle, bool distanceMatrix, bool parallelSort) : parallelSort_(parallelSort), oracle_(oracle) {
  // inherited from Debug: prefix will be printed at the beginning of every msg
  this->setDebugMsgPrefix("PairCellsWithOracle");

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

void ttk::rpd::PairCellsWithOracle::callOracle(const std::vector<std::vector<double>> &points, rpd::MultidimensionalDiagram &oracle, double threshold) {
  ripser::ripser(points, oracle, threshold, 1, false, false);

  std::sort(oracle[1].begin(), oracle[1].end(),
            [](const rpd::PersistencePair & x, const rpd::PersistencePair & y){
              if (x.second.second == y.second.second) { //reverse co-lexicographic order to match Ripser's choice
                Simplex const& t1 = x.second.first;
                Simplex const& t2 = y.second.first;
                return ! std::lexicographical_compare(t1.rbegin(), t1.rend(), t2.rbegin(), t2.rend());
              }
              else
                return x.second.second < y.second.second;
            });
}

void ttk::rpd::PairCellsWithOracle::run() {
  Timer tm{};
  for (auto const& [b,d] : oracle_[1]) {
    if (d.second < inf)
      bound = std::max(bound, d.second);
  }
  printMsg("Upper bound from the PD oracle: " + std::to_string(bound));
  initializeWithBound();
  printMsg("Initialized (#p=" + std::to_string(N) + ", #e=" + std::to_string(n_edges) + ")", 0, tm.getElapsedTime());
  pairCellsWithOracle();
  printMsg("Complete", 1, tm.getElapsedTime());
}

void ttk::rpd::PairCellsWithOracle::initializeWithBound() {
  //construct edges
  graph.resize(N);
  for (id_t i=0; i < N; ++i) {
    for (id_t j=i+1; j < N; ++j) {
      if (DM(i,j) <= bound) {
        edgeToIndex[{i,j}] = edges.size();
        edges.push_back({{i, j}, DM(i,j)});
        graph[i].insert(j);
        graph[j].insert(i);
      }
    }
  }
  n_edges = edges.size();

  //sort edges
  edgesIndices.resize(n_edges);
  edgesOrder.resize(n_edges);
  edgesPartner.resize(n_edges, -1);
  std::iota(edgesIndices.begin(), edgesIndices.end(), 0);
  if (parallelSort_)
    TTK_PSORT(globalThreadNumber_, edgesIndices.begin(), edgesIndices.end(), compEdges)
  else
    std::sort(edgesIndices.begin(), edgesIndices.end(), compEdges);
  for (unsigned i=0; i<edgesIndices.size(); ++i)
    edgesOrder[edgesIndices[i]] = i;
}

void ttk::rpd::PairCellsWithOracle::pairCellsWithOracle() {
  for (const auto &[birth,death] : oracle_[1]) {
    if (death.second < inf) {
      const id_t e_id = edgeToIndex[{birth.first[0], birth.first[1]}];
      const id_t t_id = triangles.size();
      triangles.push_back({{death.first[0], death.first[1], death.first[2]}, death.second});
      boundaries.push_back({edgeToIndex[{death.first[0], death.first[1]}],
                            edgeToIndex[{death.first[1], death.first[2]}],
                            edgeToIndex[{death.first[0], death.first[2]}]});
      cascadeEdges.emplace_back();
      trianglesPartner.push_back(e_id);
      edgesPartner[e_id] = t_id;
      eliminateBoundaryWithOracle(t_id, e_id);
    }
  }
}

void ttk::rpd::PairCellsWithOracle::eliminateBoundaryWithOracle(id_t t_id, id_t e_id) {
  BoundaryContainer bc (boundaries[t_id], n_edges);
  id_t youngest_id = -1;
  while (youngest_id != e_id) {
    youngest_id = *std::max_element(boundaries[t_id].begin(), boundaries[t_id].end(), quickCompEdges);
    if (youngest_id == e_id)
      return;
    cascadeEdges[t_id].push_back(youngest_id);
    if (edgesPartner[youngest_id] != -1)
      bc.exclusiveAddBoundary(boundaries[edgesPartner[youngest_id]]);
    else {
      const auto& [i,j] = edges[youngest_id].e;
      const double d = edges[youngest_id].d;
      for (id_t const& k : graph[i]) {
        if (DM(i,k) < d && DM(j,k) < d) {
          bc.exclusiveAddBoundary({youngest_id,
                                   edgeToIndex[{std::min(i,k),std::max(i,k)}],
                                   edgeToIndex[{std::min(j,k),std::max(j,k)}]});
          break;
        }
      }
    }
  }
}

void ttk::rpd::PairCellsWithOracle::getGenerators(std::vector<Generator> &generators) const {
  for (unsigned i=0; i<triangles.size(); ++i) {
    EdgeSet boundary;
    for (id_t const& e_id : boundaries[i])
      boundary.emplace_back(edges[e_id].e);
    generators.emplace_back(boundary, std::make_pair(edges[trianglesPartner[i]].d, triangles[i].d));
  }
}

void ttk::rpd::PairCellsWithOracle::getCascades(std::vector<Cascade> &cascades, std::vector<EdgeSet> &critical) const {
  fillRNG(critical);
  for (auto const& c : cascadeEdges) {
    critical[2].emplace_back(edges[c[0]].e);
    Cascade cascade;
    for (id_t const& e_id : c)
      cascade.emplace_back(edges[e_id].e);
    cascades.push_back(cascade);
  }
}

void ttk::rpd::PairCellsWithOracle::getCascades(EdgeSet &cascades, std::vector<EdgeSet> &critical) const {
  fillRNG(critical);
  std::set<id_t> cascadeSet;
  for (auto const& c : cascadeEdges) {
    critical[2].emplace_back(edges[c[0]].e);
    for (unsigned i=1; i<c.size(); ++i)
      cascadeSet.insert(c[i]);
  }
  for (const id_t &e_id : cascadeSet)
    cascades.emplace_back(edges[e_id].e);
}

void ttk::rpd::PairCellsWithOracle::fillRNG(std::vector<EdgeSet> &critical) const {
  critical = std::vector<EdgeSet>(3);
  for (auto const& [birth, death] : oracle_[0]) {
    if (death.second < rpd::inf)
      critical[0].emplace_back(death.first[0], death.first[1]);
  }
  for (auto const& [birth, death] : oracle_[1])
    critical[1].emplace_back(birth.first[0], birth.first[1]);
}