# Topological Autoencoders++: Enhanced Cycle-Aware Dimensionality Reduction

This repository contains the code used for the manuscript referenced below.

*Topological Autoencoders++: Enhanced Cycle-Aware Dimensionality Reduction*

## Installation

Tested on Ubuntu 22.04.5 LTS. 

Choose between the automatic installation with a script or the step-by-step installation.

### Automatic installation

Run the `install.sh` script (change by `install+CUDA.sh` if you own a CUDA-compatible device and want to use it):
```
chmod +x install.sh
./install.sh
```

### Step-by-step installation 

1. Install dependencies

   ```
   sudo apt install -y cmake-qt-gui ninja-build libboost-system-dev libopengl-dev libxcursor-dev
   sudo apt install -y qtbase5-dev qtchooser qt5-qmake qtbase5-dev-tools qttools5-dev qtxmlpatterns5-dev-tools libqt5x11extras5-dev libqt5svg5-dev qttools5-dev libqt5x11extras5-dev libqt5svg5-dev qtxmlpatterns5-dev-tools 
   sudo apt install -y libcgal-dev
   sudo apt install -y python3-sklearn 
   sudo apt install -y git
   ```
   
   To enable the comparison with UMAP (and PCA, MDS, Isomap, t-SNE), do:
   ```
   pip install umap-learn
   ```

2. Download Torch

   First, go to the root of this repository.
   Choose a version of LibTorch to download on https://pytorch.org/get-started/locally/. If you own a CUDA-compatible device and want to use it, download a CUDA version. Otherwise, download the CPU version. Note that we tested the installation with LibTorch 2.4.0 with CUDA 12.1.
   
   Example for the CPU version:
   ```
   wget https://download.pytorch.org/libtorch/cpu/libtorch-cxx11-abi-shared-with-deps-2.4.0%2Bcpu.zip
   unzip libtorch-cxx11-abi-shared-with-deps-2.4.0%2Bcpu.zip
   rm libtorch-cxx11-abi-shared-with-deps-2.4.0%2Bcpu.zip
   ```
   
   Example for the CUDA version:
   ```
   wget https://download.pytorch.org/libtorch/cu121/libtorch-cxx11-abi-shared-with-deps-2.4.0%2Bcu121.zip
   unzip libtorch-cxx11-abi-shared-with-deps-2.4.0+cu121.zip
   rm libtorch-cxx11-abi-shared-with-deps-2.4.0+cu121.zip
   ```

3. Install ParaView

   ```
   git clone https://github.com/topology-tool-kit/ttk-paraview.git
   cd ttk-paraview
   mkdir build
   cd build
   cmake -G Ninja -DCMAKE_BUILD_TYPE=Release -DPARAVIEW_USE_PYTHON=ON -DPARAVIEW_INSTALL_DEVELOPMENT_FILES=ON -DCMAKE_INSTALL_PREFIX=../install ..
   ninja install
   ```
   Set the following environment variable (replace `3.11` by your version of Python)
   ```
   PV_PREFIX=`pwd`/../install
   export PATH=$PATH:$PV_PREFIX/bin
   export LD_LIBRARY_PATH=$PV_PREFIX/lib:$LD_LIBRARY_PATH
   export PYTHONPATH=$PV_PREFIX/lib/python3.11/site-packages
   ```

4. Install TTK
   
   ```
   cd ../../ttk-tcdr
   mkdir build
   cd build
   paraviewPath=`pwd`/../../ttk-paraview/install/lib/cmake/paraview-5.13
   torchPath=`pwd`/../../libtorch/share/cmake/Torch/
   cmake -G Ninja -DCMAKE_INSTALL_PREFIX=../install -DParaView_DIR=$paraviewPath -DTorch_DIR=$torchPath ..
   ninja install
   ```
   Set the following environment variable (replace `3.11` by your version of Python)
   ```
   TTK_PREFIX=`pwd`/../install
   export PV_PLUGIN_PATH=$TTK_PREFIX/bin/plugins/TopologyToolKit
   export LD_LIBRARY_PATH=$TTK_PREFIX/lib:$LD_LIBRARY_PATH
   export PYTHONPATH=$PYTHONPATH:$TTK_PREFIX/lib/python3.11/site-packages
   ```

## Reproduce experiments

1. To reproduce the comparison tables in the manuscript (figures 1, 13 and 14), run:

   ```
   python scripts/runComparison.py
   ```
   
   You can specify the files to run with the option `-f` or `--files` and the methods to compare with the option `-m` or `--methods`, e.g.:
   ```
   python scripts/runComparison.py --files Twist COIL20-1 -m TopoAE TopoAE++
   ```
   
   The file names to use are `3Clusters`, `Twist`, `K4`, `K5`, `COIL20-1`, `MoCap` and `SingleCell` and the available methods are listed below:

   | Method name to use | Implementation |
   |--------------------|----------------|
   | PCA                | scikit-learn   |
   | MDS                | scikit-learn   |
   | Isomap             | scikit-learn   |
   | LLE                | scikit-learn   |
   | tSNE               | scikit-learn   |
   | UMAP               | umap-learn     |
   | TopoMap            | TTK            |
   | TopoAE             | TTK            |
   | TopoAE+W1          | TTK            |
   | TopoAE++           | TTK            |

2. Run a ParaView state.
   ```
   paraview states/SingleCell.pvsm
   ```