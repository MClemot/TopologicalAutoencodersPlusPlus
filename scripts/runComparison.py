import argparse

import Compute
import Figure
import Metrics

available_files = ["3Clusters", "Twist", "K4", "K5", "COIL20-1", "MoCap", "SingleCell"]
available_methods = ["PCA", "MDS", "Isomap", "tSNE", "UMAP", "TopoMap", "TopoAE", "TopoAE+W1", "TopoAE++"]
available_metrics = ["W0", "W1", "RMSE", "LC", "TA", "Trust", "Cont"]

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--files', nargs='*')
    parser.add_argument('-m', '--methods', nargs='*')
    parser.add_argument('-is', '--imagesize', nargs=1, default=500)
    args = parser.parse_args()

    if args.files == ["*"] or args.files is None:
        args.files = available_files
    if args.methods == ["*"] or args.methods is None:
        args.methods = available_methods

    for file in args.files:
        for method in args.methods:
            Compute.compute(file, method)
            Figure.figure(file, method, args.imagesize)

    Metrics.export("./scripts_results/metrics.csv", args.files, args.methods, available_metrics)