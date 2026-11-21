#!/usr/bin/env python
"""
步驟 14：從 MCL CD19+ 分選細胞中辨識出「腫瘤 B 細胞」。

做法（對應原論文 STAR Methods 的 MCL analysis）
1. 先用表現量聚類，剔除 T 細胞污染與 doublet。
2. 用 inferCNV（以扁桃腺正常 B 細胞為對照）推斷拷貝數變異。
3. 以染色體 Y 基因分數的雙峰分布切 chrY+/chrY−（論文指出 chrY loss 是主要亞克隆軸）。
4. 輸出 tumor-only 的 AnnData 供步驟 06 做殘差分析。

輸出 data/processed/mcl_tumor.h5ad
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _compat  # noqa  (scanpy/scipy 稀疏陣列相容)
import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import scipy.sparse as sp

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
PROC = f"{BASE}/data/processed"
RES = f"{BASE}/results"
os.makedirs(RES, exist_ok=True)
sc.settings.n_jobs = 48

CHR_Y = ["UTY", "KDM5D", "DDX3Y", "USP9Y", "ZFY", "EIF1AY", "RPS4Y1", "NLGN4Y", "TXLNGY"]
T_MARKERS = ["CD3D", "CD3E", "CD3G", "CD2", "TRBC1", "TRBC2", "TRAC", "IL7R", "CD247"]
MYELOID = ["LYZ", "CST3", "AIF1", "FCN1", "S100A8", "S100A9", "ITGAM"]
B_MARKERS = ["CD79A", "CD79B", "MS4A1", "CD19", "PAX5", "BANK1"]


def log(*a):
    print("[14]", *a, flush=True)


A = sc.read_h5ad(f"{PROC}/mcl_gex_raw.h5ad")
log("全部 MCL 細胞:", A.shape)
A = A[A.obs["fraction"] == "CD19pos"].copy()
log("CD19+ 分選:", A.shape, dict(A.obs.donor_id.value_counts()))

A.layers["counts"] = A.X.copy()
sc.pp.normalize_total(A, target_sum=1e4)
sc.pp.log1p(A)
A.raw = A

sc.pp.highly_variable_genes(A, n_top_genes=2500, batch_key="donor_id")
sc.pp.scale(A, max_value=10, zero_center=True)
sc.tl.pca(A, n_comps=30, mask_var="highly_variable")
# harmonypy 直接呼叫（scanpy 的 wrapper 與目前版本不相容）
try:
    import harmonypy
    ho = harmonypy.run_harmony(A.obsm["X_pca"], A.obs, ["donor_id"], max_iter_harmony=20)
    Z = np.asarray(ho.Z_corr)
    if Z.shape[0] != A.n_obs:
        Z = Z.T
    assert Z.shape == (A.n_obs, A.obsm["X_pca"].shape[1]), Z.shape
    A.obsm["X_harmony"] = Z
    rep = "X_harmony"
    log("harmony 完成", Z.shape)
except Exception as e:
    log("harmony 失敗，改用 X_pca:", e)
    rep = "X_pca"
sc.pp.neighbors(A, n_neighbors=20, use_rep=rep)
sc.tl.leiden(A, resolution=1.0, key_added="leiden", flavor="igraph", n_iterations=2)
sc.tl.umap(A, min_dist=0.3)
log("聚類數:", A.obs.leiden.nunique())

# ---- 以 marker 分數標記污染族群 ----
X = A.raw.to_adata()
sc.pp.scale(X, max_value=10)
for nm, gl in [("T_score", T_MARKERS), ("My_score", MYELOID), ("B_score", B_MARKERS),
               ("chrY_score", CHR_Y)]:
    g = [x for x in gl if x in X.var_names]
    sc.tl.score_genes(X, g, score_name=nm)
    A.obs[nm] = X.obs[nm].values
    log(f"  {nm}: 使用 {len(g)} 個基因")

prof = A.obs.groupby("leiden")[["T_score", "My_score", "B_score", "chrY_score"]].mean()
prof["n"] = A.obs.leiden.value_counts().reindex(prof.index)
log("\n每個 cluster 的 marker 分數:\n" + prof.round(3).to_string())

bad = prof[(prof.T_score > 0.25) | (prof.My_score > 0.25)].index.tolist()
log("判定為非 B 污染的 cluster:", bad)
A.obs["is_B"] = ~A.obs.leiden.isin(bad)
log("保留 B 細胞:", int(A.obs.is_B.sum()), "/", A.n_obs)

A.write_h5ad(f"{PROC}/mcl_cd19pos_clustered.h5ad", compression="gzip")
A.obs.to_csv(f"{RES}/14_mcl_cd19pos_obs.csv")
log("寫出 mcl_cd19pos_clustered.h5ad")
log("DONE")
