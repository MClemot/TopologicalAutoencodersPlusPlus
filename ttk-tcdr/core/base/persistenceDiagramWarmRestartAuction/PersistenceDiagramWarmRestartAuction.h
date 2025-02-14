#pragma once

#include <PersistenceDiagramAuction.h>

using ValuesPair = std::pair<double,double>;
using RipsPersistencePair = std::pair<std::pair<std::vector<int>, double>,
                                      std::pair<std::vector<int>, double>>;

using KDT = ttk::PersistenceDiagramAuction::KDT; //it could be better to instead using KDT = ttk::KDTree<double, std::array<double,2>>;

namespace ttk {
  template<typename T>
  class PersistenceDiagramWarmRestartAuction : public Debug {
  public:
    PersistenceDiagramWarmRestartAuction(const std::vector<T> &goodDiagram) {
      setDebugMsgPrefix("PersistenceDiagramWarmRestartAuction");

      std::vector<double> coordinates(0);
      std::vector<std::vector<double>> weights(1);
      for(unsigned i=0; i<goodDiagram.size(); ++i) {
        const ValuesPair &g = getPair(goodDiagram[i]);
        goods.emplace_back(g.first, g.second, false, i);
        coordinates.push_back(g.first);
        coordinates.push_back(g.second);
        weights[0].push_back(0.);
      }

      kdt = std::make_unique<KDT>(true, wasserstein_);
      correspondence_kdt_map = kdt->build(coordinates.data(), goodDiagram.size(), 2, weights, 1);
    }

    void setNewBidder(const std::vector<T> &bidderDiagram) {
      bidders.resize(0);

      for(unsigned i=0; i<bidderDiagram.size(); ++i) {
        const ValuesPair &b = getPair(bidderDiagram[i]);
        Bidder bidder (b.first, b.second, false, i);
        bidder.setPositionInAuction(i);
        bidders.emplace_back(bidder);
      }
    }

    void reinitializeGoodsPrice() {
      for (Good &g : goods)
        g.setPrice(0.);
    }

    double runAuction(std::vector<MatchingType> &matchings) {
      PersistenceDiagramAuction auction(bidders, goods, wasserstein_,
                                        1., 1., delta_,
                                        *kdt, correspondence_kdt_map);
      Timer t;

      matchings.resize(0);
      reinitializeGoodsPrice();

      const double w = auction.run(matchings);

      printMsg("Auction completed", 1.0, t.getElapsedTime());

      return w;
    }

    inline void setWasserstein(double wasserstein) {
      wasserstein_ = wasserstein;
    }

    inline void setDelta(double delta) {
      delta_ = delta;
    }

  private:
    double wasserstein_ {2.};
    double delta_ {0.01};

    std::unique_ptr<KDT> kdt;
    std::vector<KDT*> correspondence_kdt_map;

    GoodDiagram goods;
    BidderDiagram bidders;

    static inline ValuesPair getPair(const T& p);
  };

  template <>
  inline ValuesPair PersistenceDiagramWarmRestartAuction<RipsPersistencePair>::getPair(const RipsPersistencePair& p) {
    return {p.first.second, p.second.second};
  }

  template <>
  inline ValuesPair PersistenceDiagramWarmRestartAuction<PersistencePair>::getPair(const PersistencePair& p) {
    return {p.birth.sfValue, p.death.sfValue};
  }
}