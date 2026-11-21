#!/usr/bin/env python
"""
步驟 17：載入 MCL 的 ATAC peak 矩陣（與 GEX 同一個 filtered_feature_bc_matrix）。

CellRanger-ARC 的矩陣同時包含 "Gene Expression" 與 "Peaks" 兩種 feature。
步驟 11 只取了 GEX，這裡取 Peaks。

輸出 data/processed/mcl_atac_raw.h5ad
"""
import os, gzip
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import scipy.io as sio
import scipy.sparse as sp

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
MCLD = f"{BASE}/data/extracted/MCL"
PROC = f"{BASE}/data/processed"


def log(*a):
    print("[17]", *a, flush=True)


meta = pd.read_csv(f"{MCLD}/metadata/cellranger_arc_metadata.csv")
gex = meta[(meta["type"] == "RNA") & (meta["library_name"].str.contains("pos"))]
log("CD19+ 樣本:\n", gex[["gem_id", "library_name", "donor_id"]].to_string(index=False))

adatas = []
for _, row in gex.iterrows():
    gem, lib, donor = row["gem_id"], row["library_name"], row["donor_id"]
    d = f"{MCLD}/BCLLATLAS_64_65/{gem}/filtered_feature_bc_matrix"
    log("讀", lib)
    feats = pd.read_csv(f"{d}/features.tsv.gz", sep="\t", header=None,
                        names=["id", "name", "type", "chrom", "start", "end"])
    bcs = [l.strip() for l in gzip.open(f"{d}/barcodes.tsv.gz", "rt")]
    X = sio.mmread(f"{d}/matrix.mtx.gz").T.tocsr()      # cells x features
    ispeak = (feats["type"] == "Peaks").values
    log("   total features", X.shape, " peaks", int(ispeak.sum()))
    Xp = X[:, ispeak]
    vp = feats[ispeak].reset_index(drop=True)
    a = ad.AnnData(X=Xp.astype(np.float32),
                   obs=pd.DataFrame(index=[f"{gem}_{b}" for b in bcs]),
                   var=pd.DataFrame({"chrom": vp["chrom"].values,
                                     "start": vp["start"].values,
                                     "end": vp["end"].values},
                                    index=vp["name"].values))
    a.obs["gem_id"] = gem; a.obs["library_name"] = lib; a.obs["donor_id"] = donor
    adatas.append(a)
    log("   peaks AnnData", a.shape)

# 兩位病人的 peak set 不同 → 先各自存，之後用座標重疊對應
for a in adatas:
    donor = a.obs["donor_id"].iloc[0]
    f = f"{PROC}/mcl_atac_{donor}.h5ad"
    a.write_h5ad(f, compression="gzip")
    log("寫出", f, a.shape)
log("DONE")
