#pragma once

// VTK Module
#include <ttkDimensionReductionMetricsModule.h>

// VTK Includes
#include <ttkAlgorithm.h>

// TTK Base Includes
#include <DimensionReductionMetrics.h>

class TTKDIMENSIONREDUCTIONMETRICS_EXPORT ttkDimensionReductionMetrics
  : public ttkAlgorithm // we inherit from the generic ttkAlgorithm class
  , protected ttk::DimensionReductionMetrics // and we inherit from the base class
{
private:
  bool SelectInputFieldsWithRegexp{false};
  std::string InputRegexpString{".*"};
  std::vector<std::string> InputScalarFields{};

  bool SelectRepresentationFieldsWithRegexp{false};
  std::string RepresentationRegexpString{".*"};
  std::vector<std::string> RepresentationScalarFields{};

public:
  static ttkDimensionReductionMetrics *New();
  vtkTypeMacro(ttkDimensionReductionMetrics, ttkAlgorithm);

  void SetInputScalarFields(const std::string &s) {
    InputScalarFields.push_back(s);
    Modified();
  }

  void ClearInputScalarFields() {
    InputScalarFields.clear();
    Modified();
  }

  void SetRepresentationScalarFields(const std::string &s) {
    RepresentationScalarFields.push_back(s);
    Modified();
  }

  void ClearRepresentationScalarFields() {
    RepresentationScalarFields.clear();
    Modified();
  }

  vtkSetMacro(SelectInputFieldsWithRegexp, bool);
  vtkGetMacro(SelectInputFieldsWithRegexp, bool);

  vtkSetMacro(InputRegexpString, const std::string &);
  vtkGetMacro(InputRegexpString, std::string);

  vtkSetMacro(SelectRepresentationFieldsWithRegexp, bool);
  vtkGetMacro(SelectRepresentationFieldsWithRegexp, bool);

  vtkSetMacro(RepresentationRegexpString, const std::string &);
  vtkGetMacro(RepresentationRegexpString, std::string);

  vtkSetMacro(Wasserstein, double);
  vtkGetMacro(Wasserstein, double);

  vtkSetMacro(SampleSize, int);
  vtkGetMacro(SampleSize, int);

  vtkSetMacro(NeighborhoodSize, unsigned);
  vtkGetMacro(NeighborhoodSize, unsigned);

protected:
  ttkDimensionReductionMetrics();
  ~ttkDimensionReductionMetrics() override = default;

  int FillInputPortInformation(int port, vtkInformation *info) override;

  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
