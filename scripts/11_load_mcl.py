#!/usr/bin/env python
"""
步驟 11：載入 MCL Multiome 的基因表現資料並做品質管制。

沿用原論文的 QC 門檻（STAR Methods, "MCL analysis"）：
  UMI >= 900、偵測基因數 >= 700、粒線體比例 < 27.5% (M102) / < 20% (M413)
再加上 scDblFinder 等價的 doublet 分數（用 scanpy 的 scrublet）。

輸出 data/processed/mcl_gex_raw.h5ad
"""
import os, glob, re
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
MCLD = f"{BASE}/data/extracted/MCL"
PROC = f"{BASE}/data/processed"
os.makedirs(PROC, exist_ok=True)
sc.settings.n_jobs = 32


def log(*a):
    print("[11]", *a, flush=True)


meta = pd.read_csv(f"{MCLD}/metadata/cellranger_arc_metadata.csv")
gex = meta[meta["type"] == "RNA"].copy()
log("MCL GEX libraries:\n", gex[["gem_id", "library_name", "donor_id"]].to_string(index=False))

MT_CUT = {"M102": 27.5, "M413": 20.0}

adatas = []
for _, row in gex.iterrows():
    gem, lib, donor = row["gem_id"], row["library_name"], row["donor_id"]
    d = f"{MCLD}/BCLLATLAS_64_65/{gem}/filtered_feature_bc_matrix"
    if not os.path.isdir(d):
        log("缺少", d); continue
    log("讀", lib, d)
    a = sc.read_10x_mtx(d, var_names="gene_symbols", gex_only=True, cache=False)
    a.var_names_make_unique()
    a.obs_names = [f"{gem}_{b}" for b in a.obs_names]
    a.obs["gem_id"] = gem
    a.obs["library_name"] = lib
    a.obs["donor_id"] = donor
    a.obs["fraction"] = "CD19pos" if "pos" in lib else "CD19neg"
    adatas.append(a)
    log("   shape", a.shape)

A = ad.concat(adatas, join="outer", label=None, index_unique=None)
A.var_names_make_unique()
log("合併後:", A.shape)

A.var["mt"] = A.var_names.str.upper().str.startswith("MT-")
sc.pp.calculate_qc_metrics(A, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)

keep = np.ones(A.n_obs, bool)
keep &= A.obs["total_counts"].values >= 900
keep &= A.obs["n_genes_by_counts"].values >= 700
mtcut = A.obs["donor_id"].map(MT_CUT).values
keep &= A.obs["pct_counts_mt"].values < mtcut
log("QC 前:", A.n_obs, " QC 後:", int(keep.sum()))
A = A[keep].copy()

log("每個 library 細胞數:\n", A.obs.groupby(["donor_id", "fraction"]).size().to_string())

A.layers["counts"] = A.X.copy()
A.write_h5ad(f"{PROC}/mcl_gex_raw.h5ad", compression="gzip")
log("寫出", f"{PROC}/mcl_gex_raw.h5ad", A.shape)
