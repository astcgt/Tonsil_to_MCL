#!/usr/bin/env python
"""
步驟 13：把 R 匯出的 Multiome 矩陣組裝成 h5ad，並切出「正常 B 細胞參考」。

輸出
----
data/processed/tonsil_multiome_rna_all.h5ad   全部 68,749 細胞的 RNA
data/processed/ref_normalB_multiome.h5ad      46,216 個正常 B 細胞（NBC_MBC/GCBC/PC）
data/processed/tonsil_multiome_atac_B.h5ad    同一批 B 細胞的 ATAC peak 矩陣
"""
import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import scipy.io as sio
import scipy.sparse as sp

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
EXP = f"{BASE}/data/processed/multiome_export"
PROC = f"{BASE}/data/processed"
B_LEVELS = ["NBC_MBC", "GCBC", "PC"]


def log(*a):
    print("[13]", *a, flush=True)


def sanitize(obs):
    """h5ad 不接受混合型別欄位：非數值一律轉字串。"""
    for c in obs.columns:
        s = obs[c]
        if pd.api.types.is_numeric_dtype(s):
            continue
        obs[c] = s.astype(str).fillna("NA")
    return obs


md = pd.read_csv(f"{EXP}/multiome_metadata.csv", low_memory=False)
md = md.set_index("cell_barcode")
md = sanitize(md)
log("metadata", md.shape)

# ---------------- RNA ----------------
log("讀 rna_counts.mtx ...")
X = sio.mmread(f"{EXP}/rna_counts.mtx").T.tocsr().astype(np.float32)   # genes x cells -> cells x genes
genes = [l.strip() for l in open(f"{EXP}/rna_genes.txt")]
bcs = [l.strip() for l in open(f"{EXP}/rna_barcodes.txt")]
log("RNA", X.shape, len(genes), len(bcs))
A = ad.AnnData(X=X, obs=md.loc[bcs].copy(), var=pd.DataFrame(index=genes))
A.var_names_make_unique()
A.layers["counts"] = A.X.copy()

for r in ["umap", "harmony_rna", "pca"]:
    f = f"{EXP}/reduction_{r}.csv"
    if os.path.exists(f):
        e = pd.read_csv(f, index_col=0).loc[A.obs_names]
        A.obsm[f"X_{r}"] = e.values.astype(np.float32)
        log("加入降維", r, e.shape)

A.write_h5ad(f"{PROC}/tonsil_multiome_rna_all.h5ad", compression="gzip")
log("寫出 tonsil_multiome_rna_all.h5ad", A.shape)

# ---------------- 正常 B 參考 ----------------
isB = A.obs["annotation_level_1"].isin(B_LEVELS).values
ref = A[isB].copy()
ref.obs["annotation"] = ref.obs["annotation_20230508"].astype(str)
ref.obs["compartment"] = ref.obs["annotation_level_1"].astype(str)
log("正常 B 參考", ref.shape)
log(ref.obs["compartment"].value_counts().to_string())
ref.write_h5ad(f"{PROC}/ref_normalB_multiome.h5ad", compression="gzip")
log("寫出 ref_normalB_multiome.h5ad")

# ---------------- ATAC（同一批 B 細胞） ----------------
fa = f"{EXP}/atac_counts.mtx"
if os.path.exists(fa):
    log("讀 atac_counts.mtx ...")
    Xa = sio.mmread(fa).T.tocsr().astype(np.float32)
    peaks = [l.strip() for l in open(f"{EXP}/atac_peaks.txt")]
    abcs = [l.strip() for l in open(f"{EXP}/atac_barcodes.txt")]
    log("ATAC", Xa.shape, len(peaks))
    At = ad.AnnData(X=Xa, obs=md.loc[abcs].copy(), var=pd.DataFrame(index=peaks))
    keep = At.obs["annotation_level_1"].isin(B_LEVELS).values
    At = At[keep].copy()
    At.obs["annotation"] = At.obs["annotation_20230508"].astype(str)
    # 把 peak 名稱 (chr-start-end) 拆成座標欄位
    pr = pd.Series(At.var_names).str.rsplit("-", n=2, expand=True)
    At.var["chrom"] = pr[0].values
    At.var["start"] = pd.to_numeric(pr[1], errors="coerce").values
    At.var["end"] = pd.to_numeric(pr[2], errors="coerce").values
    At.write_h5ad(f"{PROC}/tonsil_multiome_atac_B.h5ad", compression="gzip")
    log("寫出 tonsil_multiome_atac_B.h5ad", At.shape)
    # 同時輸出 peak BED（之後用來量化 MCL fragments）
    bed = At.var[["chrom", "start", "end"]].dropna().copy()
    bed["start"] = bed["start"].astype(int); bed["end"] = bed["end"].astype(int)
    bed.to_csv(f"{PROC}/tonsil_atlas_peaks.bed", sep="\t", header=False, index=False)
    log("寫出 tonsil_atlas_peaks.bed", bed.shape)
log("DONE")
