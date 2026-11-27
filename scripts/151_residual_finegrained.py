#!/usr/bin/env python
"""步驟 151：用細分後的參考重算殘差 —— 粗解析度製造了多少假陽性？

邏輯
----
步驟 15 把 ncsMBC 當成**一格**（6,788 顆細胞的平均）。
步驟 140 顯示那一格裡至少有 4 個可再現的次狀態，而腫瘤只對應其中 3 個
（substate 3 / 5 / 4），完全避開第 4 個（substate 6，ZEB2+/FCRL5+）。

後果：**任何在 ncsMBC 內部有梯度的基因，都會拿到假的殘差。**
  若某基因在 substate 3 很高、在 substate 6 很低，
  那麼「ncsMBC 整格平均」會低估腫瘤該有的期望值 → 殘差被灌水。

做法
----
把參考從 46 個狀態改成「記憶區間用 9 個次狀態 + 其餘狀態照舊」，
其他一切（對位基因、softmax 溫度、top-k、期望值算法）與步驟 15 完全相同，
然後比較兩套殘差。只有解析度變了，所以差異完全歸因於解析度。

輸出 results/151_*
"""
import os
import sys
import json

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "8"

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _compat  # noqa

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.stats import spearmanr

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
PROC = f"{BASE}/data/processed"
RES = f"{BASE}/results"
TOPK = 5
TEMP = 0.02
MIN_CELLS = 10


def log(*a):
    print("[151]", *a, flush=True)


def dense(X):
    return X.toarray() if sp.issparse(X) else np.asarray(X)


def soft_match(Q, P, temperature=TEMP, topk=TOPK):
    Qz = (Q - Q.mean(1, keepdims=True)) / (Q.std(1, keepdims=True) + 1e-8)
    Pz = (P - P.mean(1, keepdims=True)) / (P.std(1, keepdims=True) + 1e-8)
    C = Qz @ Pz.T / Q.shape[1]
    W = np.zeros_like(C)
    kth = np.argsort(-C, axis=1)[:, :topk]
    for i in range(C.shape[0]):
        s = C[i, kth[i]] / temperature
        s = np.exp(s - s.max())
        W[i, kth[i]] = s / s.sum()
    return C, W


def pseudobulk(X, labels, donors, names, min_cells=MIN_CELLS):
    rows, keep = [], []
    lab = pd.Series(labels)
    for st in pd.unique(lab):
        idx = np.where(lab.values == st)[0]
        per = []
        for dn in pd.unique(donors[idx]):
            ii = idx[donors[idx] == dn]
            if len(ii) >= min_cells:
                per.append(X[ii].mean(0))
        if not per and len(idx) >= min_cells:
            per = [X[idx].mean(0)]
        if per:
            rows.append(np.mean(per, 0)); keep.append(st)
    return pd.DataFrame(np.vstack(rows), index=keep, columns=names)


def main():
    ref = sc.read_h5ad(f"{PROC}/ref_normalB_multiome.h5ad")
    sc.pp.normalize_total(ref, target_sum=1e4); sc.pp.log1p(ref)
    ref.obs["state"] = ref.obs["annotation"].astype(str)

    sub = pd.read_csv(f"{RES}/140_ref_substates.csv.gz", index_col=0)
    fine = ref.obs["state"].astype(str).copy()
    hit = sub.index.intersection(ref.obs_names)
    s = sub.loc[hit, "substate"].astype(str)
    s = s[s != "dropped"]
    fine.loc[s.index] = "mem_" + s.values
    ref.obs["fine"] = fine
    log("粗解析度狀態數: %d → 細解析度: %d"
        % (ref.obs.state.nunique(), ref.obs.fine.nunique()))

    mcl = sc.read_h5ad(f"{PROC}/mcl_cd19pos_clustered.h5ad")
    if mcl.raw is not None:
        mobs = mcl.obs.copy(); mcl = mcl.raw.to_adata(); mcl.obs = mobs
    s15 = pd.read_csv(f"{RES}/15_mcl_obs.csv", index_col=0)
    mcl.obs["is_tumor"] = s15.reindex(mcl.obs_names)["is_tumor"].astype(str).values

    common = ref.var_names.intersection(mcl.var_names)
    ref = ref[:, common].copy(); mcl = mcl[:, common].copy()
    hv = [g for g in json.load(open(f"{RES}/15_matching_genes.json")) if g in set(common)]
    log("共同基因 %d，對位基因 %d" % (len(common), len(hv)))

    Xr = dense(ref.X).astype(np.float32)
    donors = ref.obs["donor_id"].astype(str).to_numpy()
    out = {}
    resid = {}
    for tag, col in (("coarse", "state"), ("fine", "fine")):
        P = pseudobulk(Xr, ref.obs[col].astype(str).to_numpy(), donors, list(common))
        log("%s 參考 pseudobulk: %s" % (tag, P.shape))
        gi = [list(common).index(g) for g in hv]
        Q = dense(mcl.X).astype(np.float32)
        C, W = soft_match(Q[:, gi], P[hv].to_numpy(np.float32))
        exp = (W @ P.to_numpy(np.float32)).astype(np.float32)
        r = Q - exp
        tum = (mcl.obs.is_tumor.astype(str) == "True").to_numpy()
        non = (mcl.obs.is_tumor.astype(str) == "False").to_numpy()
        d = pd.DataFrame({"gene": common,
                          "mean_observed": Q[tum].mean(0),
                          "mean_expected": exp[tum].mean(0),
                          "resid_tum": r[tum].mean(0),
                          "resid_non": r[non].mean(0)})
        d["tumour_specific"] = d.resid_tum - d.resid_non
        d["match_corr"] = np.nan
        resid[tag] = d.set_index("gene")
        out[tag] = {"n_states": int(P.shape[0]),
                    "mean_match_corr": float(C.max(1)[tum].mean()),
                    "mean_abs_resid_tumour": float(np.abs(r[tum]).mean()),
                    "mean_abs_resid_nontumour": float(np.abs(r[non]).mean())}
        log("  腫瘤對位相關係數中位 %.3f | mean|resid| 腫瘤 %.4f, 非腫瘤 %.4f, 比值 %.2f"
            % (np.median(C.max(1)[tum]), out[tag]["mean_abs_resid_tumour"],
               out[tag]["mean_abs_resid_nontumour"],
               out[tag]["mean_abs_resid_tumour"] / out[tag]["mean_abs_resid_nontumour"]))

    A, B = resid["coarse"], resid["fine"]
    M = pd.DataFrame({"coarse": A.tumour_specific, "fine": B.tumour_specific,
                      "expected_coarse": A.mean_expected,
                      "expected_fine": B.mean_expected,
                      "observed": A.mean_observed})
    M["delta"] = M.fine - M.coarse
    M["shrink_pct"] = 100 * (1 - M.fine.abs() / M.coarse.abs().clip(lower=1e-6))
    M.round(4).to_csv(f"{RES}/151_residual_coarse_vs_fine.csv")

    rho = spearmanr(M.coarse, M.fine)
    log("\n粗 vs 細的腫瘤特異殘差 Spearman = %.3f" % rho.statistic)

    top_c = M.nlargest(100, "coarse")
    log("\n=== 粗解析度的 top-100 基因，改用細解析度後縮水多少 ===")
    log("  中位縮水 %.1f%%；%d / 100 個縮水超過一半；%d 個掉出前 500 名"
        % (top_c.shrink_pct.median(), int((top_c.shrink_pct > 50).sum()),
           int((M.fine.rank(ascending=False).reindex(top_c.index) > 500).sum())))

    PANEL = ["CCND1", "SOX11", "SLC7A11", "SLC3A2", "MAFG", "BACH1", "NFE2L2",
             "CCDC88A", "RGS12", "ADAMTS6", "TCF4", "ZNF804A", "FOXP1",
             "GAB1", "GAB2", "ARHGAP24", "MAST4"]
    log("\n=== 重點基因：粗 vs 細 ===")
    log("%-10s %10s %10s %8s %9s %9s %7s" %
        ("gene", "粗殘差", "細殘差", "縮水%", "粗排名", "細排名", "觀察值"))
    rc = M.coarse.rank(ascending=False); rf = M.fine.rank(ascending=False)
    pan = {}
    for g in PANEL:
        if g not in M.index:
            continue
        log("%-10s %10.3f %10.3f %7.0f%% %9d %9d %7.3f"
            % (g, M.loc[g, "coarse"], M.loc[g, "fine"], M.loc[g, "shrink_pct"],
               int(rc[g]), int(rf[g]), M.loc[g, "observed"]))
        pan[g] = {"coarse": float(M.loc[g, "coarse"]), "fine": float(M.loc[g, "fine"]),
                  "shrink_pct": float(M.loc[g, "shrink_pct"]),
                  "rank_coarse": int(rc[g]), "rank_fine": int(rf[g])}

    log("\n=== 縮水最多的 20 個基因（粗解析度前 300 名裡）===")
    cand = M[rc <= 300].nlargest(20, "shrink_pct")
    log("\n" + cand[["coarse", "fine", "shrink_pct"]].round(3).to_string())

    out["spearman_coarse_vs_fine"] = float(rho.statistic)
    out["top100_median_shrink_pct"] = float(top_c.shrink_pct.median())
    out["top100_n_shrink_over_50pct"] = int((top_c.shrink_pct > 50).sum())
    out["panel"] = pan
    out["most_shrunk"] = {g: float(v) for g, v in cand.shrink_pct.items()}
    json.dump(out, open(f"{RES}/151_finegrained.json", "w"), indent=2, ensure_ascii=False)
    log("\nDONE → results/151_*")


if __name__ == "__main__":
    main()
