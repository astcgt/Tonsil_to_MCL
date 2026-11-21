#!/usr/bin/env bash
# 建立分析環境 (conda env: tonsilmcl)  —— 不動 base
set -u
export MAMBA_NO_BANNER=1
ENVDIR=/raid/jupyter-alz.ee09/Tonsil_to_MCL/env/tonsilmcl
mamba create -y -p "$ENVDIR" -c conda-forge -c bioconda \
  python=3.11 \
  r-base=4.3 r-seurat r-seuratobject r-matrix r-hdf5r r-remotes r-biocmanager \
  bioconductor-singlecellexperiment bioconductor-zellkonverter bioconductor-basilisk \
  scanpy anndata h5py leidenalg python-igraph scikit-misc harmonypy \
  numba pyarrow fastparquet statsmodels seaborn matplotlib \
  pytables openpyxl tqdm joblib \
  2>&1 | tail -20
echo "=== conda env built ==="
"$ENVDIR/bin/python" -c "import scanpy, anndata; print('scanpy', scanpy.__version__)"
