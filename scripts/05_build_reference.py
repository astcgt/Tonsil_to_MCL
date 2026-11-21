#!/usr/bin/env python
"""
步驟 05：建立「正常扁桃腺 B 細胞分化參考」與「MCL 查詢集」。

設計重點
--------
* 參考 (reference)：Tonsil Atlas 的 **Multiome** B 細胞（NBC/MBC、GCBC、PC）。
  刻意選 Multiome 而非 scRNA-seq，因為 MCL 查詢資料也是 Multiome（單核），
  同實驗室、同化學、同單核製備 → 把技術差異降到最低。
* 查詢 (query)：MCL 兩位病人（M102、M413）的 CD19+ 分選細胞。
* 兩者取基因交集，輸出成兩個 h5ad 給後續步驟。

輸出
----
data/processed/ref_normalB_multiome.h5ad
data/processed/query_mcl.h5ad
"""
import os, sys, glob, json, argparse
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import scipy.sparse as sp

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
PROC = f"{BASE}/data/processed"
os.makedirs(PROC, exist_ok=True)
sc.settings.n_jobs = 32


def log(*a):
    print("[05]", *a, flush=True)


# ----------------------------------------------------------------------------
def find_annotation_col(obs):
    """自動找出最細緻的細胞型別註解欄位。"""
    cands = [c for c in obs.columns
             if "annotation" in c.lower() or "cell_type" in c.lower() or c.lower() == "celltype"]
    # 以「唯一值數量」排序，取最細緻但不過碎的
    scored = []
    for c in cands:
        try:
            n = obs[c].astype(str).nunique()
        except Exception:
            continue
        if 2 <= n <= 200:
            scored.append((n, c))
    scored.sort(reverse=True)
    return scored, cands


# ----------------------------------------------------------------------------
def load_reference(h5ad_path):
    log("讀取參考:", h5ad_path)
    A = sc.read_h5ad(h5ad_path)
    log("  shape =", A.shape)
    scored, cands = find_annotation_col(A.obs)
    log("  候選註解欄位:", cands)
    log("  (n_unique, col):", scored[:8])
    return A


# ----------------------------------------------------------------------------
def load_mcl(mcl_dir):
    """讀 CellRanger-ARC 的 filtered_feature_bc_matrix（每個樣本一個資料夾）。"""
    adatas = {}
    # 找所有 filtered_feature_bc_matrix.h5 或資料夾
    h5s = sorted(glob.glob(os.path.join(mcl_dir, "**", "filtered_feature_bc_matrix.h5"),
                           recursive=True))
    mtxdirs = sorted(glob.glob(os.path.join(mcl_dir, "**", "filtered_feature_bc_matrix"),
                               recursive=True))
    log("找到 h5:", len(h5s), " mtx dir:", len(mtxdirs))
    for p in h5s:
        sample = _sample_name(p)
        log("  讀", sample, p)
        a = sc.read_10x_h5(p, gex_only=True)
        a.var_names_make_unique()
        adatas[sample] = a
    if not adatas:
        for p in mtxdirs:
            sample = _sample_name(p)
            log("  讀", sample, p)
            a = sc.read_10x_mtx(p, var_names="gene_symbols", gex_only=True)
            a.var_names_make_unique()
            adatas[sample] = a
    return adatas


def _sample_name(p):
    parts = [x for x in p.split(os.sep) if x not in ("outs", "filtered_feature_bc_matrix",
                                                     "filtered_feature_bc_matrix.h5")]
    return parts[-1]


# ----------------------------------------------------------------------------
if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref_h5ad", required=True)
    ap.add_argument("--mcl_dir", required=True)
    ap.add_argument("--mode", default="inspect", choices=["inspect", "build"])
    args = ap.parse_args()

    A = load_reference(args.ref_h5ad)
    print(A.obs.head(3).to_string())
    print("obs columns:", list(A.obs.columns))
    if args.mode == "inspect":
        M = load_mcl(args.mcl_dir)
        for k, v in M.items():
            log(k, v.shape)
        sys.exit(0)
