#!/usr/bin/env python
"""
步驟 19：在「扣掉正常分化」的殘差上推論轉錄因子活性，找出最異常的調控子。

兩個層次互相驗證
----------------
A. 觀測 vs 期望的 TF 活性差（單細胞層級）
   對 observed 與 expected 兩個矩陣分別跑 decoupler ULM + CollecTRI，
   取差值 → ΔTF activity。這是「這個細胞的 TF 活性，比它應該有的狀態多/少多少」。
B. 直接在殘差矩陣上跑 ULM（pseudobulk 層級）
   殘差像是一個 contrast/logFC 矩陣，ULM 本來就是為這種輸入設計的。

兩個層次的排名若一致，結論較可信。另以 MCL 檢體內的非腫瘤 B 細胞作陰性對照。

輸出 results/19_*
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _compat  # noqa  (scanpy/scipy 稀疏陣列相容)
import os, json
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import decoupler as dc
import scipy.sparse as sp

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
PROC = f"{BASE}/data/processed"
RES = f"{BASE}/results"
os.makedirs(RES, exist_ok=True)
sc.settings.n_jobs = 48


def log(*a):
    print("[19]", *a, flush=True)


def dense(X):
    return X.toarray() if sp.issparse(X) else np.asarray(X)


# ---------------- TF-target 先驗 ----------------
net = dc.op.collectri(organism="human")
log("CollecTRI:", net.shape, "| TF 數:", net["source"].nunique())

A = sc.read_h5ad(f"{PROC}/mcl_with_residual.h5ad")
log("MCL（含殘差）:", A.shape)
log("欄位:", [c for c in A.obs.columns if c in
              ("donor_id", "leiden", "is_B", "is_tumor", "matched_normal_state", "match_corr")])

obs = A.obs.copy()
genes = A.var_names

# ---------------- A. 單細胞 ΔTF 活性 ----------------
log("跑 ULM on observed ...")
Aobs = ad.AnnData(X=A.X.copy(), obs=obs, var=A.var.copy())
dc.mt.ulm(Aobs, net=net, tmin=5, verbose=False)
act_obs = Aobs.obsm["score_ulm"]

log("跑 ULM on expected ...")
Aexp = ad.AnnData(X=A.layers["expected"].copy(), obs=obs, var=A.var.copy())
dc.mt.ulm(Aexp, net=net, tmin=5, verbose=False)
act_exp = Aexp.obsm["score_ulm"]

common_tf = act_obs.columns.intersection(act_exp.columns)
d_act = act_obs[common_tf] - act_exp[common_tf]
d_act.index = A.obs_names
log("ΔTF activity:", d_act.shape)
d_act.to_csv(f"{RES}/19_delta_tf_activity_cells.csv.gz", compression="gzip")

# 以「病人 × 腫瘤/非腫瘤」匯總
grp_cols = [c for c in ["donor_id", "is_tumor"] if c in obs.columns]
if "is_tumor" not in obs.columns:
    obs["is_tumor"] = obs.get("is_B", True)
    grp_cols = ["donor_id", "is_tumor"]

summ = d_act.groupby([obs[c].astype(str).values for c in grp_cols]).mean().T
summ.columns = ["_".join(map(str, c)) if isinstance(c, tuple) else str(c) for c in summ.columns]
summ.to_csv(f"{RES}/19_delta_tf_by_group.csv")

tumor_cols = [c for c in summ.columns if c.endswith("True")]
if tumor_cols:
    summ["tumor_mean"] = summ[tumor_cols].mean(axis=1)
    summ["tumor_min_abs"] = summ[tumor_cols].abs().min(axis=1)
    # 兩位病人方向一致才算數
    summ["consistent"] = (np.sign(summ[tumor_cols]).nunique(axis=1) == 1)
    rank = summ.sort_values("tumor_mean", ascending=False)
    log("\n=== 腫瘤中活性「上升」最多的 TF（兩位病人方向一致者） ===")
    log("\n" + rank[rank.consistent].head(25).round(3).to_string())
    log("\n=== 腫瘤中活性「下降」最多的 TF ===")
    log("\n" + rank[rank.consistent].tail(20).round(3).to_string())
    rank.to_csv(f"{RES}/19_tf_ranking_singlecell.csv")

# ---------------- B. pseudobulk 殘差上的 ULM ----------------
log("\n建 pseudobulk 殘差 ...")
key_cols = [c for c in ["donor_id", "leiden"] if c in obs.columns]
R = A.layers["residual"]
rows, idx = [], []
for k, sub in obs.groupby(key_cols, observed=True):
    if len(sub) < 30:
        continue
    ii = A.obs_names.get_indexer(sub.index)
    rows.append(np.asarray(R[ii]).mean(axis=0))
    idx.append("_".join(map(str, k)) if isinstance(k, tuple) else str(k))
PB = pd.DataFrame(np.vstack(rows), index=idx, columns=genes)
log("pseudobulk 殘差:", PB.shape)
PB.to_csv(f"{RES}/19_pseudobulk_residual.csv.gz", compression="gzip")

pb_ad = ad.AnnData(X=PB.values.astype(np.float32),
                   obs=pd.DataFrame(index=PB.index),
                   var=pd.DataFrame(index=PB.columns))
dc.mt.ulm(pb_ad, net=net, tmin=5, verbose=False)
sc_pb = pb_ad.obsm["score_ulm"]
pv_pb = pb_ad.obsm.get("padj_ulm", None)
sc_pb.to_csv(f"{RES}/19_pseudobulk_tf_score.csv")
res = pd.DataFrame({"tf": sc_pb.columns,
                    "mean_score": sc_pb.mean(0).values,
                    "min_abs": sc_pb.abs().min(0).values,
                    "n_pos": (sc_pb > 0).sum(0).values,
                    "n_neg": (sc_pb < 0).sum(0).values}).sort_values("mean_score", ascending=False)
res.to_csv(f"{RES}/19_tf_ranking_pseudobulk.csv", index=False)
log("\n=== pseudobulk 殘差上 TF 活性最高 25 ===")
log("\n" + res.head(25).round(3).to_string(index=False))
log("\n=== 最低 15 ===")
log("\n" + res.tail(15).round(3).to_string(index=False))
log("DONE")
