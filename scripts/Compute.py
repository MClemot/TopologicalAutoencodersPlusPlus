import paraview
import numpy as np
from sklearn.decomposition import PCA
from sklearn.manifold import MDS, Isomap, LocallyLinearEmbedding, TSNE
from umap import UMAP
import time
import os

paraview.compatibility.major = 5
paraview.compatibility.minor = 13

#### import the simple module from the paraview
from paraview.simple import *

methods = dict()
methods["TopoMap"] = "TopoMap (IEEE VIS 2020)"
methods["TopoAE"] = "Autoencoder"
methods["TopoAE1"] = "Autoencoder"
methods["TopoAE++"] = "Autoencoder"
methods["Carriere"] = "Autoencoder"

skmethods = dict()
skmethods["PCA"] = PCA(n_components=2)
skmethods["MDS"] = MDS(n_components=2)
skmethods["Isomap"] = Isomap(n_components=2, n_neighbors=8)
skmethods["LLE"] = LocallyLinearEmbedding(n_components=2)
skmethods["tSNE"] = TSNE(n_components=2)
skmethods["UMAP"] = UMAP(n_components=2)

dir_path = os.path.dirname(os.path.realpath(__file__))
inputpath = dir_path+"/../data/"
savedatapath = dir_path+"/../scripts_results/data/"
savefigurepath = dir_path+"/../scripts_results/figure/"

try:
    os.mkdir(dir_path+"/../scripts_results/")
    os.mkdir(savedatapath)
    os.mkdir(savefigurepath)
except FileExistsError:
    pass


def header(file):
    return file in ["3Clusters", "COIL20-1", "MoCap", "SingleCell"]

def threshold_v(file):
    if file in ["SingleCell", "MoCap"]:
        return 10
    return .5

def tries(method):
    if method in ["TopoAE", "TopoAE1", "TopoAE++", "Carriere"]:
        return 10
    else:
        return 1

def setColumns(file, object):
    if file[:6] == "COIL20":
        try:
            object.SelectFieldswithaRegexp = 1
            object.Regexp = '([1-9]|[1-9][0-9]{1,2}|10[01][0-9]|102[0-4])'
        except AttributeError:
            object.SelectInputFieldswithaRegexp = 1
            object.InputRegexp = '([1-9]|[1-9][0-9]{1,2}|10[01][0-9]|102[0-4])'
    elif file == "SingleCell":
        try:
            object.SelectFieldswithaRegexp = 1
            object.Regexp = 'expression.*'
        except AttributeError:
            object.SelectInputFieldswithaRegexp = 1
            object.InputRegexp = 'expression.*'
    elif file == "MoCap":
        try:
            object.SelectFieldswithaRegexp = 1
            object.Regexp = '[a-z]*[0-9]'
        except AttributeError:
            object.SelectInputFieldswithaRegexp = 1
            object.InputRegexp = '[a-z]*[0-9]'
    if file == "3Clusters":
        object.InputColumns = ['x', 'y', 'z']
    else:
        object.InputColumns = ['Field 0', 'Field 1', 'Field 2']

def setMethod(file, method, object):
    object.Method = methods[method]
    if method in ["TopoAE", "TopoAE1", "TopoAE++", "Carriere"]:
        object.Numberofepochs = 1000

        if file == "SingleCell":
            object.Regularizationcoefficient = 0.0001

        object.Hiddenlayers = '128 32'
        if file == "MoCap":
            object.Hiddenlayers = '32 32'

        if method == "TopoAE":
            object.Lossfunction = 'Topological Autoencoder (ICML 2020)'
        elif method == "TopoAE1":
            object.Lossfunction = 'Topological Autoencoder with dimension 1'
        elif method == "TopoAE++":
            object.Lossfunction = 'Asymmetric Cascade Autoencoder'
        elif method == "Carriere":
            object.Lossfunction = 'W1 regularization'
            if file == "MoCap":
                object.Regularizationcoefficient = 10
            if file == "SingleCell":
                object.Regularizationcoefficient = .01


def compute(file, method):

    ### methods handled by sklearn / umap
    if method in skmethods:
        if file == "3Clusters":
            X = np.genfromtxt(inputpath + file + ".csv", delimiter=',')[:,1:] #remove class number
        else:
            X = np.genfromtxt(inputpath + file + ".csv", delimiter=',')
        t = time.time()
        if header(file):
            Y = skmethods[method].fit_transform(X[1:])
        else:
            Y = skmethods[method].fit_transform(X)
        t = time.time() - t
        if file == "3Clusters":
            np.savetxt(savedatapath + "{}_{}.csv".format(file, method),
                       np.concatenate([np.genfromtxt(inputpath + file + ".csv", delimiter=',')[1:,:1], Y], axis=1),
                       delimiter=',', header="\"ClusterId\", \"Component_0\",\"Component_1\"")
        else:
            np.savetxt(savedatapath + "{}_{}.csv".format(file, method), Y, delimiter=',', header="\"Component_0\",\"Component_1\"")

    ### methods handled by TTK
    else:
        # create a new 'CSV Reader'
        csv = CSVReader(registrationName='csv', FileName=[inputpath + file + ".csv"])
        csv.HaveHeaders = header(file)

        # create a new 'TTK DimensionReduction'
        runs = []
        times = []
        best_va = (float('inf'), float('inf'))
        best_id = -1
        csv.UpdatePipeline()
        for k in range(tries(method)):
            runs.append(TTKDimensionReduction(registrationName='TTKDimensionReduction1', Input=csv, ModulePath='default'))
            setColumns(file, runs[-1])
            setMethod(file, method, runs[-1])

            t = time.time()
            runs[-1].UpdatePipeline()
            times.append(time.time() - t)

            tTKDimensionReductionMetrics = TTKDimensionReductionMetrics(registrationName='TTKDimensionReductionMetrics1',
                                                                         Input=runs[-1],
                                                                         Representation=runs[-1])
            setColumns(file, tTKDimensionReductionMetrics)
            tTKDimensionReductionMetrics.RepresentationColumns = ['Component_0', 'Component_1']
            tTKDimensionReductionMetrics.UpdatePipeline()
            metrics = servermanager.Fetch(tTKDimensionReductionMetrics).GetRowData().GetArray(1)
            w1 = metrics.GetValue(1)
            disto = metrics.GetValue(11)
            if w1 < best_va[0]:
                best_va = w1, disto
                best_id = k

        SaveData(savedatapath + "{}_{}.csv".format(file, method), runs[best_id])
        t = times[best_id]

    return t