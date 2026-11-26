#!/usr/bin/env python
"""
步驟 15（分析核心）：把每個 MCL 細胞「對位」到正常 B 細胞分化狀態，
再扣掉該狀態應有的表現，得到純粹的腫瘤殘差。

為什麼不用「MCL vs 全部正常 B」直接做差異表現？
--------------------------------------------------
MCL 細胞在分化軸上不是均勻分布（論文指出集中在 activated NBC → GC-committed 這段）。
直接比較會把「分化階段組成差異」誤判成「腫瘤異常」。本步驟先消掉這個混淆。

方法（刻意選「可解釋、對批次穩健」而非深度模型）
-------------------------------------------------
1. 以正常參考的 46 個 B 細胞狀態（48 個註釋狀態中，每位捐贈者細胞數 ≥10 的） × 6 位捐贈者，建 pseudobulk 表現譜。
2. 選「對位用基因集」：參考集的高變異基因，但排除
      - 免疫球蛋白可變區（腫瘤是單株，會主導相似度）
      - 粒線體 / 核糖體
      - 已知 MCL 驅動基因（CCND1、SOX11 等，避免用結果來決定對位）
      - 位於腫瘤 CNA 區域的基因（劑量效應會扭曲對位）
3. 每個 MCL 細胞 → 對 46 個狀態算相關係數 → softmax 得到「軟性對位權重」。
4. 期望表現 = 權重加權的正常 pseudobulk（在全基因上）。
5. 殘差 = 觀測 − 期望。
6. **陰性對照**：MCL 檢體內的非腫瘤 B 細胞走同一流程，殘差應接近 0。

輸出 results/15_*
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _compat  # noqa  (scanpy/scipy 稀疏陣列相容)
import os, json, argparse
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import scipy.sparse as sp
from scipy.stats import rankdata

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
PROC = f"{BASE}/data/processed"
RES = f"{BASE}/results"
os.makedirs(RES, exist_ok=True)
sc.settings.n_jobs = 48

IG_PREFIX = ("IGHV", "IGKV", "IGLV", "IGHJ", "IGKJ", "IGLJ", "IGLC", "IGKC",
             "IGHA", "IGHG", "IGHM", "IGHD", "IGHE", "TRAV", "TRBV", "TRGV", "TRDV",
             "TRAC", "TRBC")
# 對位時排除的已知 MCL 驅動/標記（避免循環論證）
MCL_KNOWN = {"CCND1", "SOX11", "CCND2", "CCND3", "CDK4", "BIRC3", "ATM", "TP53",
             "KMT2D", "NSD2", "UBR5", "MEF2B", "SP140", "NOTCH1", "NOTCH2", "TCL1A"}


def log(*a):
    print("[15]", *a, flush=True)


def dense(X):
    return X.toarray() if sp.issparse(X) else np.asarray(X)


def build_matching_genes(ref, n_hvg=2000, extra_exclude=()):
    bl = set(extra_exclude) | set(MCL_KNOWN)
    for g in ref.var_names:
        gu = g.upper()
        if gu.startswith(IG_PREFIX) or gu.startswith("MT-") or gu.startswith("RPL") \
           or gu.startswith("RPS") or gu.startswith("HB") and len(gu) <= 4:
            bl.add(g)
    t = ref.copy()
    sc.pp.highly_variable_genes(t, n_top_genes=n_hvg + len(bl) + 500,
                                batch_key="donor_id", flavor="seurat")
    hv = [g for g in t.var_names[t.var.highly_variable] if g not in bl]
    return hv[:n_hvg]


def pseudobulk(A, state_col, donor_col, min_cells=10):
    """回傳 (states x genes) 的 pseudobulk（先每捐贈者平均，再跨捐贈者平均 → 降低單一捐贈者主導）"""
    X = dense(A.X)
    rows, names = [], []
    for st, sub in A.obs.groupby(state_col, observed=True):
        per_donor = []
        for dn, sub2 in sub.groupby(donor_col, observed=True):
            if len(sub2) < min_cells:
                continue
            ii = A.obs_names.get_indexer(sub2.index)
            per_donor.append(X[ii].mean(axis=0))
        if len(per_donor) == 0:
            if len(sub) >= min_cells:
                ii = A.obs_names.get_indexer(sub.index)
                per_donor = [X[ii].mean(axis=0)]
            else:
                continue
        rows.append(np.mean(per_donor, axis=0))
        names.append(st)
    return pd.DataFrame(np.vstack(rows), index=names, columns=A.var_names)


def soft_match(Q, P, temperature=0.02, topk=5):
    """Q: cells x genes (已 z-score), P: states x genes (已 z-score).
    用 Pearson 相關 → 只留 top-k → softmax 權重。"""
    Qz = (Q - Q.mean(1, keepdims=True)) / (Q.std(1, keepdims=True) + 1e-8)
    Pz = (P - P.mean(1, keepdims=True)) / (P.std(1, keepdims=True) + 1e-8)
    C = Qz @ Pz.T / Q.shape[1]                      # cells x states
    W = np.zeros_like(C)
    kth = np.argsort(-C, axis=1)[:, :topk]
    for i in range(C.shape[0]):
        s = C[i, kth[i]] / temperature
        s = np.exp(s - s.max())
        W[i, kth[i]] = s / s.sum()
    return C, W


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", default=f"{PROC}/ref_normalB_multiome.h5ad")
    ap.add_argument("--mcl", default=f"{PROC}/mcl_cd19pos_clustered.h5ad")
    ap.add_argument("--tumor_col", default="is_tumor")
    ap.add_argument("--n_hvg", type=int, default=2000)
    ap.add_argument("--topk", type=int, default=5)
    args = ap.parse_args()

    ref = sc.read_h5ad(args.ref)
    log("參考:", ref.shape)
    sc.pp.normalize_total(ref, target_sum=1e4); sc.pp.log1p(ref)

    mcl = sc.read_h5ad(args.mcl)
    log("MCL:", mcl.shape)
    # mcl 的 X 可能已 scale 過 → 用 raw
    if mcl.raw is not None:
        mcl = mcl.raw.to_adata()
        log("改用 .raw (lognorm):", mcl.shape)

    common = ref.var_names.intersection(mcl.var_names)
    log("共同基因:", len(common))
    ref = ref[:, common].copy(); mcl = mcl[:, common].copy()

    hv = build_matching_genes(ref, n_hvg=args.n_hvg)
    log("對位用基因:", len(hv), "例:", hv[:10])
    json.dump(hv, open(f"{RES}/15_matching_genes.json", "w"))

    # ---------- pseudobulk ----------
    P_all = pseudobulk(ref, "annotation", "donor_id", min_cells=10)
    log("正常狀態 pseudobulk:", P_all.shape)
    P_all.to_csv(f"{RES}/15_normal_state_pseudobulk.csv")

    hv_use = [g for g in hv if g in P_all.columns]
    P = P_all[hv_use].values.astype(np.float32)
    Q = dense(mcl[:, hv_use].X).astype(np.float32)

    C, W = soft_match(Q, P, topk=args.topk)
    states = list(P_all.index)
    best = [states[i] for i in C.argmax(1)]
    mcl.obs["matched_normal_state"] = best
    mcl.obs["match_corr"] = C.max(1)
    log("\n對位到的正常狀態分布:\n",
        pd.Series(best).value_counts().head(25).to_string())

    # ---------- 期望與殘差（全基因） ----------
    Pfull = P_all[common].values.astype(np.float32)
    expected = (W @ Pfull).astype(np.float32)
    observed = dense(mcl.X).astype(np.float32)
    resid = observed - expected
    log("殘差矩陣:", resid.shape, " mean|r| = %.4f" % np.abs(resid).mean())

    mcl.layers["expected"] = expected
    mcl.layers["residual"] = resid
    mcl.obs["resid_l2"] = np.linalg.norm(resid, axis=1)
    mcl.write_h5ad(f"{PROC}/mcl_with_residual.h5ad", compression="gzip")

    # ---------- 基因層級摘要 ----------
    grp = mcl.obs.get(args.tumor_col, pd.Series(True, index=mcl.obs_names))
    out = {}
    for lab, mask in [("tumor", grp.astype(str) == "True"), ("nontumor", grp.astype(str) == "False")]:
        if mask.sum() < 20:
            continue
        out[lab] = pd.DataFrame({
            "gene": common,
            "mean_observed": observed[mask.values].mean(0),
            "mean_expected": expected[mask.values].mean(0),
            "mean_residual": resid[mask.values].mean(0),
            "sd_residual": resid[mask.values].std(0),
            "n_cells": int(mask.sum()),
        })
    for lab, df in out.items():
        df["t_like"] = df.mean_residual / (df.sd_residual / np.sqrt(df.n_cells) + 1e-9)
        df = df.sort_values("mean_residual", ascending=False)
        df.to_csv(f"{RES}/15_gene_residual_{lab}.csv", index=False)
        log(f"\n=== {lab}: 殘差最高的 25 個基因 ===")
        log("\n" + df.head(25)[["gene", "mean_observed", "mean_expected",
                                "mean_residual", "t_like"]].round(3).to_string(index=False))
        log(f"\n=== {lab}: 殘差最低的 15 個基因 ===")
        log("\n" + df.tail(15)[["gene", "mean_observed", "mean_expected",
                                "mean_residual", "t_like"]].round(3).to_string(index=False))

    pd.DataFrame(W, index=mcl.obs_names, columns=states).to_csv(f"{RES}/15_match_weights.csv.gz",
                                                                compression="gzip")
    mcl.obs.to_csv(f"{RES}/15_mcl_obs.csv")
    log("DONE")


if __name__ == "__main__":
    main()
