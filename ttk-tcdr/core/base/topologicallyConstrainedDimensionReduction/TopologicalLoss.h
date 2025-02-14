/// \ingroup base
/// \class ttk::TopologicalLoss
/// \author Mattéo Clémot <matteo.clemot@univ-lyon1.fr>
/// \date July 2024.
///
/// This file defines the %TopologicalLoss class. It can compute the following
/// differentiable torch::Tensor quantities that can be used as loss function
/// terms for topological regularization in an autoencoder-based dimension
/// reduction technique :
/// (1) : "Topological Autoencoder" loss
/// (2) : 1-dimensional-extended "Topological Autoencoder" loss
/// (3) : cascade-extended "Topological Autoencoder" loss
/// (4) : asymmetric cascade-extended "Topological Autoencoder" loss
/// (5) : Wasserstein distance between 1-dimensional persistence diagrams
///

#pragma once

#ifdef TTK_ENABLE_TORCH

#include <PersistenceDiagramWarmRestartAuction.h>
#include <RipsPersistenceDiagram.h>
#include <PairCellsWithOracle.h>

#include <torch/torch.h>

namespace ttk {
  class TopologicalLoss {
  public:
    enum class REGUL : std::uint8_t {
      /** no regularization (returns 0) */
      NO_REGUL,
      /** "Topological Autoencoder" loss */
      TOPOAE,
      /** 1-dimensional-extended "Topological Autoencoder" loss */
      TOPOAE_DIM1,
      /** cascade-extended "Topological Autoencoder" loss */
      CASCADE,
      /** asymmetric cascade-extended "Topological Autoencoder" loss */
      ASYMMETRIC_CASCADE,
      /** Wasserstein distance between 1-dimensional persistence diagrams */
      W_DIM1,
    };

    TopologicalLoss(const torch::Tensor &input, std::vector<std::vector<double>> const& points, REGUL regul);

    torch::Tensor computeLoss(const torch::Tensor &latent);

  private:
    const torch::Tensor input_;
    const std::vector<std::vector<double>> &points_;
    const REGUL regul_;
    const torch::Reduction::Reduction reduction_ {torch::Reduction::Mean};
    const torch::DeviceType device {torch::kCPU};
    torch::Tensor latent_;
    int latentDimension;

    /* persistence containers */
    rpd::MultidimensionalDiagram inputPD;
    std::vector<torch::Tensor> inputCriticalPairIndices; // [0] is MST, [1] is RNG-MST, [2] is MML, [3] is strict cascade if required
    std::unique_ptr<PersistenceDiagramWarmRestartAuction<rpd::PersistencePair>> auction {nullptr};

    /* persistence computation methods */
    void precomputeInputPersistence();
    void computeLatent0Persistence(rpd::EdgeSet &latent0PD) const;
    template <typename PersistenceType>
    void computeLatent0And1Persistence(PersistenceType &latentPD) const;
    void computeLatentCascades(rpd::EdgeSet &latentCascades, std::vector<rpd::EdgeSet> &latentCritical) const;

    /* tensor tools */
    inline torch::Tensor pairsToTorch(const rpd::EdgeSet &edges) const;
    inline torch::Tensor diffDistances(const torch::Tensor &data, const torch::Tensor &indices) const;
    inline torch::Tensor diffEdgeSetMSE(const torch::Tensor &indices) const;
    torch::Tensor diffPD(const torch::Tensor &points, const rpd::Diagram &PD, const std::vector<unsigned> &indices) const;

    /* TopoAE-like distances */
    inline torch::Tensor diffRNGMML(const rpd::EdgeSetSet &latentCritical) const;
    torch::Tensor diffTopoAELoss() const;
    torch::Tensor diffTopoAELossDim1() const;
    torch::Tensor diffCascadeAELoss() const;
    torch::Tensor diffAsymmetricCascadeAELoss() const;

    /* Wasserstein distances */
    void performAuction(const rpd::Diagram &latentPD,
                        std::vector<unsigned> &directMatchingLatent,
                        std::vector<unsigned> &directMatchingInput,
                        std::vector<unsigned> &diagonalMatchingLatent);
    torch::Tensor diffW1();
  };
}

#endif