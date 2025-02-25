from Compute import *

metricId = dict()
metricId["W0"] = 0
metricId["W1"] = 1
metricId["TAE0"] = 2
metricId["TAE1"] = 3
metricId["TA"] = 4
metricId["LC"] = 5
metricId["Trust"] = 6
metricId["Cont"] = 7
metricId["LCMC"] = 8
metricId["MRREi"] = 9
metricId["MRREl"] = 10
metricId["RMSE"] = 11

def computeMetrics(file, method, metricNames):

    if not method in ttkmethods and not method in skmethods:
        return [-1 for _ in metricNames]

    input = CSVReader(registrationName='input',
                      FileName=[inputpath + file + ".csv"])
    input.HaveHeaders = header(file)

    output = CSVReader(registrationName='output',
                       FileName=[savedatapath + "{}_{}.csv".format(file, method)])

    tTKDimensionReductionMetrics = TTKDimensionReductionMetrics(registrationName='TTKDimensionReductionMetrics1',
                                                                Input=input,
                                                                Representation=output)
    setColumns(file, tTKDimensionReductionMetrics)
    tTKDimensionReductionMetrics.RepresentationColumns = ['Component_0', 'Component_1']
    tTKDimensionReductionMetrics.UpdatePipeline()
    metrics = servermanager.Fetch(tTKDimensionReductionMetrics).GetRowData().GetArray(1)

    return [metrics.GetValue(metricId[n]) for n in metricNames]


def export(csvFile, files, methods, metrics):
    f = open(csvFile, 'w')
    f.write("File, Method")
    for method in methods:
        f.write(", {}".format(method))
    f.write("\n")

    for file in files:
        f.write(file)
        metricsValues = dict()
        for method in methods:
            metricsValues[method] = computeMetrics(file, method, metrics)
        for i, metric in enumerate(metrics):
            f.write(", {}".format(metric))
            for method in methods:
                if metric in ["LC", "TA", "Trust", "Cont"]:
                    f.write(", {:#.2f}".format(metricsValues[method][i]))
                else:
                    f.write(", {:#.1e}".format(metricsValues[method][i]))
            f.write("\n")

    f.close()

