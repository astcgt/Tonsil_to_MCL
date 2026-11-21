#!/usr/bin/env python
"""
步驟 22：MCL 腫瘤 vs 正常 B 細胞的染色質可及性差異。

難點：兩邊的 peak 集合是各自 call 出來的（圖譜 215,204 個 vs MCL 各自 ~90–109k），
      不能直接相比。做法是把 MCL 的 peak 以座標重疊投影到「圖譜 peak 集」上
      （一個 MCL peak 若重疊多個圖譜 peak，依重疊長度分配），
      之後兩邊就在同一個座標系上比較。

輸出
----
results/22_atac_pseudobulk.csv.gz   圖譜 peak × (正常狀態 / MCL 病人) 的 CPM
results/22_atac_diff.csv            MCL vs 正常的差異可及性
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _compat  # noqa
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import scipy.sparse as sp

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
PROC = f"{BASE}/data/processed"
RES = f"{BASE}/results"
os.makedirs(RES, exist_ok=True)


def log(*a):
    print("[22]", *a, flush=True)


def parse_peaks(names):
    """圖譜是 chr-start-end；MCL 是 chr:start-end。兩種都處理。"""
    s = pd.Series(names).astype(str).str.replace(":", "-", regex=False)
    p = s.str.rsplit("-", n=2, expand=True)
    return pd.DataFrame({"Chromosome": p[0].values,
                         "Start": pd.to_numeric(p[1], errors="coerce").values,
                         "End": pd.to_numeric(p[2], errors="coerce").values})


def build_overlap_map(src_names, dst_names):
    """回傳稀疏矩陣 M (n_src x n_dst)，M[i,j] = 重疊比例（以 src peak 長度為分母）。"""
    import pyranges as pr
    s = parse_peaks(src_names); s["idx_src"] = np.arange(len(s))
    d = parse_peaks(dst_names); d["idx_dst"] = np.arange(len(d))
    s = s.dropna(); d = d.dropna()
    s[["Start", "End"]] = s[["Start", "End"]].astype(int)
    d[["Start", "End"]] = d[["Start", "End"]].astype(int)
    gs = pr.PyRanges(s); gd = pr.PyRanges(d)
    j = gs.join(gd, suffix="_dst").df
    if len(j) == 0:
        raise RuntimeError("peak 無任何重疊，檢查染色體命名")
    ov = (np.minimum(j["End"], j["End_dst"]) - np.maximum(j["Start"], j["Start_dst"])).clip(lower=0)
    srclen = (j["End"] - j["Start"]).clip(lower=1)
    w = (ov / srclen).values
    M = sp.coo_matrix((w, (j["idx_src"].values, j["idx_dst"].values)),
                      shape=(len(src_names), len(dst_names))).tocsr()
    log(f"  重疊映射: {len(j)} 組配對, 有對應的 src peak = {int((M.sum(1) > 0).sum())}/{len(src_names)}")
    return M


# ---------------- 正常 B ATAC ----------------
log("讀正常 B ATAC ...")
N = sc.read_h5ad(f"{PROC}/tonsil_multiome_atac_B.h5ad")
log("  ", N.shape)
atlas_peaks = list(N.var_names)

# ---------------- MCL ATAC 投影 ----------------
mcl_pb = {}
for donor in ["M102", "M413"]:
    f = f"{PROC}/mcl_atac_{donor}.h5ad"
    if not os.path.exists(f):
        log("缺少", f); continue
    A = sc.read_h5ad(f)
    log(f"MCL {donor} ATAC:", A.shape)
    # 只保留被判定為腫瘤的細胞
    lab = f"{RES}/20_mcl_obs_labeled.csv"
    if os.path.exists(lab):
        ob = pd.read_csv(lab, index_col=0)
        tum = ob.index[ob["is_tumor"].astype(str) == "True"]
        keep = A.obs_names.isin(tum)
        log(f"  與腫瘤標記交集: {int(keep.sum())}/{A.n_obs}")
        if keep.sum() > 100:
            A = A[keep].copy()
    M = build_overlap_map(list(A.var_names), atlas_peaks)
    v = np.asarray(A.X.sum(axis=0)).ravel()          # 每個 MCL peak 的總 count
    proj = np.asarray(M.T @ v).ravel()               # 投影到圖譜 peak
    mcl_pb[f"MCL_{donor}"] = proj
    log(f"  投影完成, 總 count {v.sum():.3g} -> {proj.sum():.3g}")
    del A

# ---------------- 正常 pseudobulk（依分化狀態） ----------------
log("建正常 pseudobulk ...")
X = N.X.tocsr()
norm_pb = {}
for st, sub in N.obs.groupby("annotation", observed=True):
    if len(sub) < 50:
        continue
    ii = N.obs_names.get_indexer(sub.index)
    norm_pb[f"normal::{st}"] = np.asarray(X[ii].sum(axis=0)).ravel()
norm_pb["normal::ALL"] = np.asarray(X.sum(axis=0)).ravel()
log("  正常狀態數:", len(norm_pb) - 1)

PB = pd.DataFrame({**norm_pb, **mcl_pb}, index=atlas_peaks)
CPM = PB / PB.sum(axis=0) * 1e6
CPM.to_csv(f"{RES}/22_atac_pseudobulk.csv.gz", compression="gzip")
log("pseudobulk 寫出:", CPM.shape)

# ---------------- 差異 ----------------
mcl_cols = [c for c in CPM.columns if c.startswith("MCL_")]
if mcl_cols:
    lm = np.log2(CPM[mcl_cols] + 1)
    ln = np.log2(CPM[["normal::ALL"]] + 1)
    d = pd.DataFrame({
        "peak": CPM.index,
        "mcl_mean_log2cpm": lm.mean(axis=1).values,
        "normal_log2cpm": ln["normal::ALL"].values,
    })
    for c in mcl_cols:
        d[f"lfc_{c}"] = (lm[c] - ln["normal::ALL"]).values
    d["lfc_mean"] = d[[f"lfc_{c}" for c in mcl_cols]].mean(axis=1)
    d["consistent"] = (np.sign(d[[f"lfc_{c}" for c in mcl_cols]]).nunique(axis=1) == 1)
    pk = parse_peaks(d["peak"])
    d["chrom"] = pk["Chromosome"].values; d["start"] = pk["Start"].values; d["end"] = pk["End"].values
    d = d.sort_values("lfc_mean", ascending=False)
    d.to_csv(f"{RES}/22_atac_diff.csv", index=False)
    log("\n=== MCL 中可及性最高（相對正常 B）的 20 個 peak ===")
    log("\n" + d[d.consistent].head(20).round(3).to_string(index=False))
    log("\n=== 最低的 10 個 ===")
    log("\n" + d[d.consistent].tail(10).round(3).to_string(index=False))
log("DONE")
