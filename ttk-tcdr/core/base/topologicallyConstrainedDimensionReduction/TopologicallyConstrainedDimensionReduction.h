/// \ingroup base
/// \class ttk::TCDR
/// \author Mattéo Clémot <matteo.clemot@univ-lyon1.fr>
/// \date 2024.
///
/// This module defines the %TCDR class
/// that embeds points into 2D, under topological constraints
///

#pragma once

// ttk includes
#include <Debug.h>
#include <DimensionReductionModel.h>
#include <TopologicalLoss.h>

namespace ttk {

  /**
   * The TCDR class provides a backend for dimension reduction using
   * autoencoders, with possible constraints on the preservation of the
   * topology of the input high dimensional point cloud when projecting in low
   * dimension
   */
  class TCDR : virtual public Debug {

  public:
    enum class OPTIMIZER : std::uint8_t {
      /** Adaptive Moment Estimation */
      ADAM = 0,
      /** Stochastic Gradient Descent */
      SGD = 1,
      /** Limited-memory Broyden–Fletcher–Goldfarb–Shanno */
      LBFGS = 2,
    };

    enum class MODEL : std::uint8_t {
      /** AutoEncoder architecture */
      AUTOENCODER = 0,
      /** AutoDecoder architecture */
      AUTODECODER = 1,
      /** Direct optimization*/
      DIRECT = 2,
    };

    using REGUL = TopologicalLoss::REGUL;

#ifdef TTK_ENABLE_TORCH

    TCDR();
    TCDR(bool useCUDA, int numberOfComponents, int epochs, double learningRate,
         OPTIMIZER optimizer, REGUL method, MODEL modelType, const std::string &architecture,
         const std::string &activation, int batchSize, bool batchNormalization,
         double regCoefficient, bool inputIsImages);

    /**
     * @brief Computes the projection with an AutoEncoder
     *
     * @param[out] outputEmbedding the final coordinates of the points
     *
     * @param[in] inputMatrix the high-dimension coordinates of the points
     *
     * @param[in] n the number of input points

  * @return 0 in case of success.
    */
    int execute(std::vector<std::vector<double>> &outputEmbedding,
                const std::vector<double> &inputMatrix,
                size_t n);

  protected:
    int NumberOfComponents{2};
    int Epochs{100};
    double LearningRate{0.01};
    OPTIMIZER Optimizer{OPTIMIZER::ADAM};
    REGUL Method{REGUL::NO_REGUL};
    MODEL ModelType{MODEL::AUTOENCODER};
    bool InputIsImages;
    std::string Architecture;
    std::string Activation;
    int BatchSize;
    bool BatchNormalization;
    double RegCoefficient{0.01};

  private:
    torch::DeviceType device {torch::kCPU};
    std::shared_ptr<DimensionReductionModel> model {nullptr};
    std::shared_ptr<torch::optim::Optimizer> torchOptimizer {nullptr};
    std::shared_ptr<TopologicalLoss> topologicalLossContainer {nullptr};

    void optimize(const torch::Tensor &input) const;
    void optimizeSimple(const torch::Tensor &input) const;

    inline void printLoss(int epoch, const torch::Tensor &loss) const;

#endif

  }; // TCDR class

}