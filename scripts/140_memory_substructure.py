#!/usr/bin/env python
"""步驟 140：把解析度拉進 ncsMBC 裡面。

問題
----
步驟 120 的時序圖上，95.8% 的腫瘤細胞塌縮成一格（ncsMBC），而且落在擴散時序的
最前端。那不是腫瘤真的「很早期」，是**參考的解析度不夠**：
48 個狀態是整個 B 細胞區間的切法，ncsMBC 這一格裡塞了 6,788 顆細胞，
它內部的結構被整個抹掉了。

做法
----
A. 只拿記憶區間的參考細胞重新嵌入（ncsMBC、csMBC、FCRL4/5+、MBC FCRL5+、
   Early MBC、Precursor MBCs、Reactivated proliferative MBCs），
   重選高變異基因、重做 PCA、用 Harmony 校正捐贈者，再 Leiden 分群。
B. **捐贈者再現性過濾**：一個次狀態若不是在 ≥5/6 位捐贈者都出現（各 ≥10 顆），
   就當成批次假象丟掉。這是整個分析最重要的一道關卡 ——
   細分很容易細分出批次效應。
C. 在這個子空間裡重算擴散擬時序（根取整體擬時序最小的那顆記憶細胞）。
D. 用同一組基因、同一個 PCA，把 19,975 顆腫瘤細胞 kNN 投影進來。
   問三件事：
     ① 腫瘤還是塌縮成一格嗎，還是散開了？
     ② 兩位病人落在同一個次狀態嗎？
     ③ 腫瘤自己在這個軸上有沒有內部梯度（＝腫瘤保留了分化方向嗎）？

排除基因：Ig 可變區（腫瘤是單株，IGHV 會主導一切）、粒線體、核糖體、
          已知 MCL driver —— 與步驟 15 同一套邏輯。

輸出 results/140_*
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
from sklearn.decomposition import PCA

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
PROC = f"{BASE}/data/processed"
RES = f"{BASE}/results"
sc.settings.n_jobs = 8

MEMORY = ["ncsMBC", "csMBC", "ncsMBC FCRL4/5+", "csMBC FCRL4/5+", "MBC FCRL5+",
          "Early MBC", "Precursor MBCs", "Reactivated proliferative MBCs"]
FOCUS = "ncsMBC"
K = 15
N_HVG = 2000
MIN_DONORS = 5
MIN_CELLS_PER_DONOR = 10
IG_PREFIX = ("IGHV", "IGKV", "IGLV", "IGHJ", "IGKJ", "IGLJ", "IGLC", "IGKC",
             "IGHA", "IGHG", "IGHM", "IGHD", "IGHE", "TRAV", "TRBV", "TRGV",
             "TRDV", "TRAC", "TRBC")
MCL_KNOWN = {"CCND1", "SOX11", "CCND2", "CCND3", "CDK4", "BIRC3", "ATM", "TP53",
             "KMT2D", "NSD2", "UBR5", "MEF2B", "SP140", "NOTCH1", "NOTCH2", "TCL1A"}


def log(*a):
    print("[140]", *a, flush=True)


def dense(X):
    return X.toarray() if sp.issparse(X) else np.asarray(X)


def zrow(M):
    M = np.asarray(M, np.float32)
    return (M - M.mean(1, keepdims=True)) / (M.std(1, keepdims=True) + 1e-8)


def main():
    ref = sc.read_h5ad(f"{PROC}/ref_normalB_multiome.h5ad")
    ref.obs["state"] = ref.obs["annotation"].astype(str)
    sub = ref[ref.obs.state.isin(MEMORY)].copy()
    del ref
    log("記憶區間參考細胞:", sub.shape)
    log("\n" + sub.obs.state.value_counts().to_string())
    sc.pp.normalize_total(sub, target_sum=1e4)
    sc.pp.log1p(sub)

    bl = set(MCL_KNOWN)
    for g in sub.var_names:
        gu = g.upper()
        if gu.startswith(IG_PREFIX) or gu.startswith("MT-") or \
           gu.startswith("RPL") or gu.startswith("RPS"):
            bl.add(g)
    t = sub.copy()
    sc.pp.highly_variable_genes(t, n_top_genes=N_HVG + len(bl) + 500,
                                batch_key="donor_id", flavor="seurat")
    hv = [g for g in t.var_names[t.var.highly_variable] if g not in bl][:N_HVG]
    del t
    log("子空間用的高變異基因:", len(hv))
    json.dump(hv, open(f"{RES}/140_subspace_genes.json", "w"))

    R = zrow(dense(sub[:, hv].X))
    pca = PCA(n_components=50, random_state=0).fit(R)
    sub.obsm["X_pca_sub"] = pca.transform(R).astype(np.float32)
    log("子空間 PCA 解釋變異: %.3f" % float(pca.explained_variance_ratio_.sum()))

    # scanpy 的 harmony 包裝在這個版本會把矩陣轉錯方向，所以直接呼叫 harmonypy
    import harmonypy
    ho = harmonypy.run_harmony(sub.obsm["X_pca_sub"], sub.obs, ["donor_id"],
                               max_iter_harmony=20)
    Zc = np.asarray(ho.Z_corr)
    if Zc.shape == (sub.n_obs, 50):
        pass
    elif Zc.shape == (50, sub.n_obs):
        Zc = Zc.T
    else:
        raise SystemExit(f"harmony 輸出形狀看不懂: {Zc.shape}")
    sub.obsm["X_harmony_sub"] = np.ascontiguousarray(Zc, dtype=np.float32)
    log("Harmony 校正完成:", sub.obsm["X_harmony_sub"].shape)

    sc.pp.neighbors(sub, n_neighbors=15, use_rep="X_harmony_sub")
    best = None
    for res in (0.3, 0.5, 0.8, 1.2):
        key = f"leiden_{res}"
        sc.tl.leiden(sub, resolution=res, key_added=key,
                     flavor="igraph", n_iterations=2, directed=False)
        ct = pd.crosstab(sub.obs[key], sub.obs.donor_id)
        ok = ((ct >= MIN_CELLS_PER_DONOR).sum(1) >= MIN_DONORS)
        log("resolution %.1f → %2d 群，其中 %2d 群通過 %d/6 捐贈者再現性"
            % (res, ct.shape[0], int(ok.sum()), MIN_DONORS))
        if best is None or int(ok.sum()) > best[1]:
            best = (key, int(ok.sum()), res)
    key, n_ok, res = best
    log("採用 resolution %.1f（%d 個可再現的次狀態）" % (res, n_ok))

    ct = pd.crosstab(sub.obs[key], sub.obs.donor_id)
    keep = ct.index[(ct >= MIN_CELLS_PER_DONOR).sum(1) >= MIN_DONORS]
    sub.obs["substate"] = sub.obs[key].astype(str)
    sub.obs.loc[~sub.obs[key].isin(keep), "substate"] = "dropped"
    log("丟掉 %d 群（捐贈者再現性不足），涉及 %d 顆細胞"
        % (ct.shape[0] - len(keep), int((sub.obs.substate == "dropped").sum())))

    # ---- 子空間內的擬時序 ----
    sc.tl.diffmap(sub, n_comps=15)
    prev = pd.read_csv(f"{RES}/120_ref_cells.csv.gz").set_index("cell")
    g_dpt = prev["dpt"].reindex(sub.obs_names).to_numpy(float)
    root = int(np.nanargmin(np.where(sub.obs.state.values == "Early MBC",
                                     g_dpt, np.inf)))
    if not np.isfinite(g_dpt[root]):
        root = int(np.nanargmin(g_dpt))
    sub.uns["iroot"] = root
    log("子空間根細胞:", sub.obs_names[root], sub.obs.state.iloc[root])
    sc.tl.dpt(sub)
    d = np.asarray(sub.obs["dpt_pseudotime"], float)
    bad = ~np.isfinite(d)
    if bad.any():
        d[bad] = np.nanmedian(d[~bad])
    sub.obs["sub_dpt"] = d

    prof = sub.obs.groupby("substate", observed=True).agg(
        n=("sub_dpt", "size"), dpt=("sub_dpt", "median"),
        dpt_q25=("sub_dpt", lambda s: s.quantile(.25)),
        dpt_q75=("sub_dpt", lambda s: s.quantile(.75)))
    comp = pd.crosstab(sub.obs.substate, sub.obs.state, normalize="index")
    prof["top_state"] = comp.idxmax(1)
    prof["top_state_frac"] = comp.max(1)
    prof["n_donors"] = pd.crosstab(sub.obs.substate, sub.obs.donor_id).ge(
        MIN_CELLS_PER_DONOR).sum(1)
    prof = prof.sort_values("dpt")
    log("\n=== 記憶區間的次狀態（按子空間擬時序排序） ===\n", prof.round(3).to_string())

    # ---- 每個次狀態的標記基因 ----
    use = sub[sub.obs.substate != "dropped"].copy()
    sc.tl.rank_genes_groups(use, "substate", method="wilcoxon", n_genes=12)
    marks = {}
    for g in use.obs.substate.cat.categories if hasattr(
            use.obs.substate, "cat") else sorted(use.obs.substate.unique()):
        try:
            marks[str(g)] = [str(x) for x in
                             use.uns["rank_genes_groups"]["names"][g][:12]]
        except Exception:
            pass
    log("\n=== 各次狀態的標記基因 ===")
    for g in prof.index:
        if g in marks:
            log("  %-10s (n=%5d, dpt %.3f, 主要來自 %s) : %s"
                % (g, prof.loc[g, "n"], prof.loc[g, "dpt"],
                   prof.loc[g, "top_state"], ", ".join(marks[g][:8])))

    sub.obs[["donor_id", "state", "substate", "sub_dpt"]].to_csv(
        f"{RES}/140_ref_substates.csv.gz", compression="gzip")

    # ---- 投影腫瘤細胞 ----
    s15 = pd.read_csv(f"{RES}/15_mcl_obs.csv", index_col=0)
    mcl = sc.read_h5ad(f"{PROC}/mcl_cd19pos_clustered.h5ad")
    if mcl.raw is not None:
        mobs = mcl.obs.copy(); mcl = mcl.raw.to_adata(); mcl.obs = mobs
    for c in ("is_tumor",):
        mcl.obs[c] = s15.reindex(mcl.obs_names)[c].values
    hv_use = [g for g in hv if g in set(mcl.var_names)]
    log("\n投影用基因: %d / %d" % (len(hv_use), len(hv)))
    Q = zrow(dense(mcl[:, hv_use].X))
    Rs = zrow(dense(sub[:, hv_use].X))
    p2 = PCA(n_components=50, random_state=0).fit(Rs)
    Rp = p2.transform(Rs).astype(np.float32)
    Qp = p2.transform(Q).astype(np.float32)

    rn2 = (Rp ** 2).sum(1)
    ind = np.empty((Qp.shape[0], K), np.int64)
    for a in range(0, Qp.shape[0], 1024):
        b = min(a + 1024, Qp.shape[0])
        d2 = rn2[None, :] - 2.0 * (Qp[a:b] @ Rp.T) + (Qp[a:b] ** 2).sum(1)[:, None]
        pt = np.argpartition(d2, K, axis=1)[:, :K]
        ind[a:b] = np.take_along_axis(pt, np.argsort(
            np.take_along_axis(d2, pt, axis=1), axis=1), axis=1)

    ss = sub.obs.substate.to_numpy().astype(str)
    sd = sub.obs.sub_dpt.to_numpy()
    vote, frac = [], []
    for row in ss[ind]:
        v, c = np.unique(row, return_counts=True)
        j = int(np.argmax(c)); vote.append(v[j]); frac.append(c[j] / K)
    out = pd.DataFrame({"cell": mcl.obs_names,
                        "donor_id": mcl.obs["donor_id"].astype(str).values,
                        "is_tumor": mcl.obs["is_tumor"].astype(str).values,
                        "substate": vote, "substate_frac": frac,
                        "sub_dpt": np.median(sd[ind], axis=1)})
    out.to_csv(f"{RES}/140_mcl_substates.csv.gz", index=False, compression="gzip")

    tum = out[out.is_tumor == "True"]
    non = out[out.is_tumor == "False"]
    vc = tum.substate.value_counts()
    log("\n=== 腫瘤細胞落在哪些次狀態 ===\n",
        (vc.head(12) / len(tum) * 100).round(1).to_string())
    log("\n最大一格佔 %.1f%%（原本 48 狀態下是 95.8%%）" % (100 * vc.iloc[0] / len(tum)))
    log("腫瘤子空間擬時序: 中位 %.3f, IQR %.3f–%.3f, 全距 %.3f"
        % (tum.sub_dpt.median(), tum.sub_dpt.quantile(.25),
           tum.sub_dpt.quantile(.75), tum.sub_dpt.max() - tum.sub_dpt.min()))
    log("同檢體正常 B: 中位 %.3f, IQR %.3f–%.3f"
        % (non.sub_dpt.median(), non.sub_dpt.quantile(.25), non.sub_dpt.quantile(.75)))
    log("\n=== 兩位病人分別落在哪 ===")
    for d_, x in tum.groupby("donor_id"):
        v = x.substate.value_counts()
        log("  %s (n=%d): %s | 擬時序中位 %.3f"
            % (d_, len(x), ", ".join("%s %.0f%%" % (k, 100 * v_ / len(x))
                                     for k, v_ in v.head(4).items()),
               x.sub_dpt.median()))

    summary = {
        "resolution": res, "n_substates_kept": int(len(keep)),
        "n_substates_dropped": int(ct.shape[0] - len(keep)),
        "donor_reproducibility_rule": f">= {MIN_DONORS}/6 donors with >= {MIN_CELLS_PER_DONOR} cells",
        "substate_profile": json.loads(prof.round(4).to_json(orient="index")),
        "markers": marks,
        "tumour": {"n": int(len(tum)),
                   "top_substates_pct": {k: round(100 * v / len(tum), 2)
                                         for k, v in vc.head(10).items()},
                   "largest_pct": float(100 * vc.iloc[0] / len(tum)),
                   "sub_dpt_median": float(tum.sub_dpt.median()),
                   "sub_dpt_iqr": [float(tum.sub_dpt.quantile(.25)),
                                   float(tum.sub_dpt.quantile(.75))]},
        "nontumour": {"n": int(len(non)),
                      "sub_dpt_median": float(non.sub_dpt.median()),
                      "sub_dpt_iqr": [float(non.sub_dpt.quantile(.25)),
                                      float(non.sub_dpt.quantile(.75))]},
        "per_donor": {d_: {"n": int(len(x)),
                           "top": x.substate.value_counts().head(4).to_dict(),
                           "sub_dpt_median": float(x.sub_dpt.median())}
                      for d_, x in tum.groupby("donor_id")},
    }
    json.dump(summary, open(f"{RES}/140_substructure.json", "w"),
              indent=2, ensure_ascii=False)
    log("\nDONE → results/140_*")


if __name__ == "__main__":
    main()
