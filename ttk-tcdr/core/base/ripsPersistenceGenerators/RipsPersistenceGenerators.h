/// \ingroup base
/// \author Mattéo Clémot <matteo.clemot@univ-lyon1.fr>
/// \date June 2024.

#pragma once

// ttk common includes
#include <Debug.h>

#include <ripser.h>
#include <PairCells.h>
#include <PairCellsWithOracle.h>

namespace ttk {

  class RipsPersistenceGenerators : virtual public Debug {
  public:
    RipsPersistenceGenerators();
    void execute(const std::vector<std::vector<double>> &points,
                 std::vector<rpd::Diagram> &diagrams,
                 std::vector<rpd::Generator> &generators) {
      if (SimplexMaximumDiameter == 42.) {
        rpd::PairCellsWithOracle::callOracle(points, diagrams);
        rpd::PairCellsWithOracle pcwo(points, diagrams, false, false);
        pcwo.setDebugLevel(debugLevel_);
        pcwo.run();
        if (!OutputCascade)
          pcwo.getGenerators(generators);
        else {
          rpd::EdgeSet cascades;
          rpd::EdgeSetSet critical;
          pcwo.getCascades(cascades, critical);
          generators.emplace_back(critical[0], std::make_pair(0.,0.)); //MST
          generators.emplace_back(critical[1], std::make_pair(1.,1.)); //RNG
          generators.emplace_back(critical[2], std::make_pair(2.,2.)); //MML
          generators.emplace_back(cascades,    std::make_pair(3.,3.)); //cascade
        }
      }
      else {
        rpd::PairCells pc (points, false, SimplexMaximumDiameter);
        pc.setDebugLevel(this->debugLevel_);
        pc.run();
        if (!OutputCascade)
          pc.getDiagramAndGenerators(diagrams, generators);
        else {
          pc.getDiagram(diagrams);
          rpd::EdgeSet cascades;
          rpd::EdgeSetSet critical;
          pc.getCascades(cascades, critical);
          generators.emplace_back(critical[0], std::make_pair(0.,0.)); //MST
          generators.emplace_back(critical[1], std::make_pair(1.,1.)); //RNG
          generators.emplace_back(critical[2], std::make_pair(2.,2.)); //MML
          generators.emplace_back(cascades,    std::make_pair(3.,3.)); //cascade
        }
      }
    }
  protected:
    double SimplexMaximumDiameter {rpd::inf};
    bool OutputCascade{false};
  };

}