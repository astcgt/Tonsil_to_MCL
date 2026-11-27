#!/usr/bin/env python
"""步驟 120：分化時序（pseudotime）與 MCL 細胞在分化流形上的位置。

回答的問題
----------
1. 那 48 個狀態是哪裡來的？ → 它們是 Human Tonsil Atlas (Massoni-Badosa 2024) 的
   公開註釋，不是我們算的。本步驟把它們排成一條「時序」，讓它可以被看見。
2. MCL 細胞落在這條時序的哪裡？

方法（刻意和步驟 15 用不同的演算法，當作獨立檢查）
--------------------------------------------------
A. 時序軸：在圖譜自己的 Harmony 校正空間上做 diffusion map + DPT，
   根細胞取 NBC（未活化初始 B 細胞）狀態的中位細胞。DPT 是擴散距離，
   不是真實時間 —— 它只排序，不給單位。
B. MCL 投影：用步驟 15 的同一組「對位基因」（已排除 Ig 可變區、粒線體/核糖體、
   已知 MCL driver、CNA 區域），每顆細胞對基因做 z-score，
   以參考細胞擬合的 PCA(50) 投影 MCL，再取 k=15 最近參考鄰居的
   UMAP 座標中位數 / DPT 中位數 / 狀態多數票。
   → 這是 kNN 標籤轉移，和步驟 15 的「相關係數 softmax」完全不同的演算法。
   兩者的狀態一致率就是一個獨立的穩健性檢查。

輸出
----
results/120_ref_cells.csv.gz     參考細胞：UMAP、DPT、狀態
results/120_state_order.csv      48 個狀態按 DPT 中位數排序
results/120_mcl_cells.csv.gz     MCL 細胞：投影 UMAP、DPT、kNN 狀態、softmax 狀態
results/120_summary.json         一致率、分布、每位病人的 DPT
"""
import os
import sys
import json

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "8"

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _compat  # noqa

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from sklearn.decomposition import PCA

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
PROC = f"{BASE}/data/processed"
RES = f"{BASE}/results"
sc.settings.n_jobs = 8
K = 15
ROOT_STATE = "NBC"


def log(*a):
    print("[120]", *a, flush=True)


def dense(X):
    return X.toarray() if sp.issparse(X) else np.asarray(X)


def zrow(M):
    M = np.asarray(M, dtype=np.float32)
    return (M - M.mean(1, keepdims=True)) / (M.std(1, keepdims=True) + 1e-8)


def main():
    # ---------- 參考 ----------
    ref = sc.read_h5ad(f"{PROC}/ref_normalB_multiome.h5ad")
    log("參考:", ref.shape)
    sc.pp.normalize_total(ref, target_sum=1e4)
    sc.pp.log1p(ref)

    hv = json.load(open(f"{RES}/15_matching_genes.json"))
    hv = [g for g in hv if g in set(ref.var_names)]
    log("對位基因:", len(hv))

    # ---------- A. 時序軸（圖譜自己的 Harmony 空間） ----------
    cached = f"{RES}/120_ref_cells.csv.gz"
    if os.path.exists(cached) and os.environ.get("FORCE_DPT") != "1":
        cdf = pd.read_csv(cached).set_index("cell").reindex(ref.obs_names)
        if cdf["dpt"].notna().all():
            log("沿用已算好的 DPT:", cached, "（FORCE_DPT=1 可重算）")
            ann = ref.obs["annotation"].astype(str).values
            dpt = cdf["dpt"].to_numpy(dtype=float)
            umap = cdf[["umap1", "umap2"]].to_numpy(dtype=float)
            ref_df = cdf.reset_index()
            return stage_b(ref, hv, ann, dpt, umap, ref_df)
    log("neighbors on X_harmony_rna ...")
    sc.pp.neighbors(ref, n_neighbors=15, use_rep="X_harmony_rna")
    log("diffmap ...")
    sc.tl.diffmap(ref, n_comps=15)
    ann = ref.obs["annotation"].astype(str).values
    root_idx = np.where(ann == ROOT_STATE)[0]
    if len(root_idx) == 0:
        raise SystemExit(f"找不到根狀態 {ROOT_STATE}")
    DM = ref.obsm["X_diffmap"][:, 1:6]
    c = DM[root_idx].mean(0)
    ref.uns["iroot"] = int(root_idx[np.argmin(((DM[root_idx] - c) ** 2).sum(1))])
    log("root cell:", ref.obs_names[ref.uns["iroot"]], ROOT_STATE)
    log("dpt ...")
    sc.tl.dpt(ref)
    dpt = np.asarray(ref.obs["dpt_pseudotime"], dtype=float)
    bad = ~np.isfinite(dpt)
    if bad.any():
        log("非有限 DPT:", int(bad.sum()), "→ 設為該狀態中位數")
        for st in np.unique(ann[bad]):
            m = (ann == st) & ~bad
            dpt[(ann == st) & bad] = np.median(dpt[m]) if m.any() else np.nanmedian(dpt[~bad])
    umap = np.asarray(ref.obsm["X_umap"], dtype=float)

    ref_df = pd.DataFrame({
        "cell": ref.obs_names,
        "donor_id": ref.obs["donor_id"].astype(str).values,
        "state": ann,
        "umap1": umap[:, 0], "umap2": umap[:, 1],
        "dpt": dpt,
    })
    ref_df.to_csv(f"{RES}/120_ref_cells.csv.gz", index=False, compression="gzip")
    return stage_b(ref, hv, ann, dpt, umap, ref_df)


def stage_b(ref, hv, ann, dpt, umap, ref_df):
    # ---------- 48 狀態排序 + 每位捐贈者細胞數 ----------
    g = ref_df.groupby("state", observed=True)
    order = pd.DataFrame({
        "n_cells": g.size(),
        "n_donors_ge10": g["donor_id"].apply(lambda s: int((s.value_counts() >= 10).sum())),
        "dpt_median": g["dpt"].median(),
        "dpt_q25": g["dpt"].quantile(0.25),
        "dpt_q75": g["dpt"].quantile(0.75),
        "umap1": g["umap1"].median(),
        "umap2": g["umap2"].median(),
    }).sort_values("dpt_median")
    used15 = set(pd.read_csv(f"{RES}/15_normal_state_pseudobulk.csv",
                             index_col=0, usecols=[0]).index.astype(str))
    order["used_in_step15"] = [st in used15 for st in order.index.astype(str)]
    order.index.name = "state"
    order.to_csv(f"{RES}/120_state_order.csv")
    log("\n狀態時序（前 12 / 後 12）:\n",
        order[["n_cells", "n_donors_ge10", "dpt_median", "used_in_step15"]].head(12).round(3).to_string())
    log("\n", order[["n_cells", "n_donors_ge10", "dpt_median", "used_in_step15"]].tail(12).round(3).to_string())
    log("用於步驟 15 的狀態數:", int(order["used_in_step15"].sum()), "/", len(order))

    # ---------- B. MCL 投影 ----------
    R = zrow(dense(ref[:, hv].X))
    pca = PCA(n_components=50, random_state=0).fit(R)
    Rp = pca.transform(R).astype(np.float32)
    log("參考 PCA 完成，解釋變異:", round(float(pca.explained_variance_ratio_.sum()), 3))
    del R

    s15 = pd.read_csv(f"{RES}/15_mcl_obs.csv", index_col=0)
    mcl = sc.read_h5ad(f"{PROC}/mcl_cd19pos_clustered.h5ad")
    if mcl.raw is not None:
        mobs = mcl.obs.copy()
        mcl = mcl.raw.to_adata()
        mcl.obs = mobs
    # is_tumor 只存在於步驟 15 的 obs 表（步驟 14 判定），這裡補回來
    for c in ("is_tumor", "matched_normal_state", "match_corr"):
        mcl.obs[c] = s15.reindex(mcl.obs_names)[c].values
    log("MCL:", mcl.shape, "腫瘤細胞:", int((mcl.obs.is_tumor.astype(str) == "True").sum()))
    miss = [h for h in hv if h not in set(mcl.var_names)]
    if miss:
        raise SystemExit(f"MCL 缺少 {len(miss)} 個對位基因，與步驟 15 不一致")
    Q = zrow(dense(mcl[:, hv].X))
    Qp = pca.transform(Q).astype(np.float32)
    del Q

    # 分塊 kNN（避免 sklearn 的 BLAS 執行緒池在這台機器上崩潰）
    rn2 = (Rp ** 2).sum(1)
    ind = np.empty((Qp.shape[0], K), dtype=np.int64)
    dist = np.empty((Qp.shape[0], K), dtype=np.float32)
    CH = 1024
    for a in range(0, Qp.shape[0], CH):
        b = min(a + CH, Qp.shape[0])
        d2 = rn2[None, :] - 2.0 * (Qp[a:b] @ Rp.T) + (Qp[a:b] ** 2).sum(1)[:, None]
        pt = np.argpartition(d2, K, axis=1)[:, :K]
        dd = np.take_along_axis(d2, pt, axis=1)
        o = np.argsort(dd, axis=1)
        ind[a:b] = np.take_along_axis(pt, o, axis=1)
        dist[a:b] = np.sqrt(np.maximum(np.take_along_axis(dd, o, axis=1), 0))
    log("kNN 完成")

    nb_state = ann[ind]
    vote, votefrac = [], []
    for row in nb_state:
        v, cnt = np.unique(row, return_counts=True)
        j = int(np.argmax(cnt))
        vote.append(v[j]); votefrac.append(cnt[j] / K)
    mcl_df = pd.DataFrame({
        "cell": mcl.obs_names,
        "donor_id": mcl.obs["donor_id"].astype(str).values,
        "is_tumor": mcl.obs["is_tumor"].astype(str).values,
        "umap1": np.median(umap[ind, 0], axis=1),
        "umap2": np.median(umap[ind, 1], axis=1),
        "dpt": np.median(dpt[ind], axis=1),
        "knn_state": vote,
        "knn_state_frac": votefrac,
        "knn_dist_mean": dist.mean(1),
    })
    mcl_df["softmax_state"] = s15.reindex(mcl_df["cell"])["matched_normal_state"].values
    mcl_df["match_corr"] = s15.reindex(mcl_df["cell"])["match_corr"].values
    mcl_df.to_csv(f"{RES}/120_mcl_cells.csv.gz", index=False, compression="gzip")

    tum = mcl_df[mcl_df.is_tumor == "True"]
    non = mcl_df[mcl_df.is_tumor == "False"]
    agree = float((tum.knn_state == tum.softmax_state).mean())
    log("\n腫瘤細胞 kNN 狀態分布:\n", tum.knn_state.value_counts().head(10).to_string())
    log("\n兩種演算法狀態一致率（腫瘤）: %.3f" % agree)

    def blk(df):
        if len(df) == 0:
            return None
        return {
            "n": int(len(df)),
            "dpt_median": float(df.dpt.median()),
            "dpt_q25": float(df.dpt.quantile(.25)),
            "dpt_q75": float(df.dpt.quantile(.75)),
            "knn_state_top": df.knn_state.value_counts().head(6).to_dict(),
            "knn_state_frac_median": float(df.knn_state_frac.median()),
            "umap_median": [float(df.umap1.median()), float(df.umap2.median())],
        }

    summary = {
        "n_states_annotated": int(len(order)),
        "n_states_used_step15": int(order["used_in_step15"].sum()),
        "root_state": ROOT_STATE,
        "k": K,
        "pca_var_explained": float(pca.explained_variance_ratio_.sum()),
        "state_agreement_knn_vs_softmax_tumour": agree,
        "tumour": blk(tum),
        "nontumour": blk(non),
        "per_donor_tumour": {d: blk(x) for d, x in tum.groupby("donor_id")},
        "ref_dpt_by_state": {k: float(v) for k, v in order["dpt_median"].items()},
        "ref_state_n": {k: int(v) for k, v in order["n_cells"].items()},
        "ref_state_umap": {k: [float(a), float(b)]
                           for k, a, b in zip(order.index, order.umap1, order.umap2)},
    }
    json.dump(summary, open(f"{RES}/120_summary.json", "w"), indent=2, ensure_ascii=False)
    log("DONE → results/120_*")


if __name__ == "__main__":
    main()
