#!/usr/bin/env python
"""
步驟 20：用 inferCNV 區分 MCL 的腫瘤細胞與檢體內殘存的正常 B 細胞。

參考組採用扁桃腺的正常 B 細胞（Multiome，同製備方式），
再用三個互相獨立的證據交叉確認腫瘤身分：
  (1) inferCNV 的染色體層級訊號（論文指出 chrY loss 是主要亞克隆軸）
  (2) CCND1 表現（t(11;14) 造成的過表現，正常成熟 B 細胞幾乎不表現）
  (3) 免疫球蛋白輕鏈限制性（單株腫瘤會偏向 kappa 或 lambda 其中一邊）

輸出 data/processed/mcl_tumor_labeled.h5ad
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
GTF = f"{BASE}/data/genome/gencode.v44.annotation.gtf.gz"
os.makedirs(RES, exist_ok=True)
sc.settings.n_jobs = 48


def log(*a):
    print("[20]", *a, flush=True)


def gene_positions(gtf_path):
    """從 GENCODE GTF 取每個 gene symbol 的染色體與起訖。"""
    import gzip
    rows = {}
    with gzip.open(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[2] != "gene":
                continue
            attr = f[8]
            name = attr.split('gene_name "')[1].split('"')[0] if 'gene_name "' in attr else None
            if name is None or name in rows:
                continue
            rows[name] = (f[0], int(f[3]), int(f[4]))
    df = pd.DataFrame.from_dict(rows, orient="index",
                                columns=["chromosome", "start", "end"])
    return df


# ---------------- 載入 ----------------
mcl = sc.read_h5ad(f"{PROC}/mcl_cd19pos_clustered.h5ad")
if mcl.raw is not None:
    obs_keep = mcl.obs.copy()
    mcl = mcl.raw.to_adata()
    mcl.obs = obs_keep
log("MCL CD19+:", mcl.shape)
mcl = mcl[mcl.obs["is_B"].astype(str) == "True"].copy()
log("去掉非 B 污染後:", mcl.shape)

ref = sc.read_h5ad(f"{PROC}/ref_normalB_multiome.h5ad")
sc.pp.normalize_total(ref, target_sum=1e4); sc.pp.log1p(ref)
# 取子集以控制運算量：每個狀態最多 150 顆
rng = np.random.default_rng(0)
idx = []
for st, sub in ref.obs.groupby("annotation", observed=True):
    ii = np.flatnonzero(ref.obs_names.isin(sub.index))
    idx.extend(rng.choice(ii, size=min(150, len(ii)), replace=False))
ref_s = ref[np.sort(np.array(idx))].copy()
ref_s.obs["cnv_group"] = "normal_B_tonsil"
log("參考子集:", ref_s.shape)

mcl.obs["cnv_group"] = "MCL_" + mcl.obs["donor_id"].astype(str)
common = ref_s.var_names.intersection(mcl.var_names)
log("共同基因:", len(common))
comb = ad.concat([ref_s[:, common], mcl[:, common]], join="inner",
                 label="src", keys=["ref", "mcl"], index_unique="-")
log("合併:", comb.shape)

# ---------------- 基因座標 ----------------
pos = gene_positions(GTF)
log("GTF 基因座標:", pos.shape)
keep = comb.var_names.intersection(pos.index)
comb = comb[:, keep].copy()
comb.var["chromosome"] = pos.loc[keep, "chromosome"].values
comb.var["start"] = pos.loc[keep, "start"].values
comb.var["end"] = pos.loc[keep, "end"].values
comb = comb[:, comb.var.chromosome.astype(str).str.match(r"^chr([0-9]+|X|Y)$")].copy()
comb.var["chromosome"] = comb.var["chromosome"].astype(str)
log("用於 inferCNV 的基因:", comb.n_vars)

# ---------------- inferCNV ----------------
import infercnvpy as cnv
cnv.tl.infercnv(comb, reference_key="cnv_group",
                reference_cat=["normal_B_tonsil"],
                window_size=100, step=10, exclude_chromosomes=None)
log("inferCNV 完成", comb.obsm["X_cnv"].shape)

cnv.tl.pca(comb, n_comps=30)
cnv.pp.neighbors(comb, n_neighbors=15)
cnv.tl.leiden(comb, resolution=1.0, key_added="cnv_leiden")
cnv.tl.cnv_score(comb, groupby="cnv_leiden")
log("CNV cluster 數:", comb.obs.cnv_leiden.nunique())

tab = comb.obs.groupby(["cnv_leiden", "cnv_group"], observed=True).size().unstack(fill_value=0)
log("\ncnv_leiden × 來源:\n" + tab.to_string())
sco = comb.obs.groupby("cnv_leiden", observed=True)["cnv_score"].mean()
log("\n每個 cnv_leiden 的 CNV 分數:\n" + sco.round(4).to_string())

comb.obs[["cnv_group", "cnv_leiden", "cnv_score", "src"]].to_csv(f"{RES}/20_cnv_obs.csv")
np.save(f"{RES}/20_X_cnv.npy", np.asarray(comb.obsm["X_cnv"].todense()
                                          if sp.issparse(comb.obsm["X_cnv"])
                                          else comb.obsm["X_cnv"]))
import json
json.dump({str(k): int(v) if np.isscalar(v) else [int(x) for x in np.ravel(v)]
           for k, v in comb.uns["cnv"]["chr_pos"].items()},
          open(f"{RES}/20_cnv_chrpos.json", "w"))

# ---------------- 判定腫瘤 ----------------
# 以 CNV cluster 的組成判定：若該 cluster 有相當比例的「正常參考細胞」混在裡面，
# 代表它的 CNV profile 與正常無異 → 落在該 cluster 的 MCL 細胞視為非腫瘤殘存 B 細胞。
comp = comb.obs.groupby(["cnv_leiden", "src"], observed=True).size().unstack(fill_value=0)
comp["ref_frac"] = comp.get("ref", 0) / comp.sum(axis=1)
REF_FRAC_CUT = 0.20
normal_like = set(comp.index[comp["ref_frac"] >= REF_FRAC_CUT])
log("\n每個 cnv_leiden 的正常參考佔比:\n" + comp.round(3).to_string())
log("判定為『正常樣』的 cluster:", sorted(map(str, normal_like)))

mcl_obs = comb.obs[comb.obs.src == "mcl"].copy()
mcl_obs["is_tumor_cnv"] = ~mcl_obs["cnv_leiden"].isin(normal_like)
log("依 CNV cluster 判定為腫瘤:", int(mcl_obs.is_tumor_cnv.sum()), "/", len(mcl_obs))
ref_score = comb.obs.loc[comb.obs.src == "ref", "cnv_score"]
log("（參考用）正常參考 CNV 分數 99 百分位 =",
    round(float(np.percentile(ref_score, 99)), 5))

# CCND1 與輕鏈限制性
X = mcl
for g in ["CCND1", "SOX11", "IGKC", "IGLC1", "IGLC2", "IGLC3"]:
    if g in X.var_names:
        v = X[:, g].X
        mcl_obs[g] = np.asarray(v.todense()).ravel() if sp.issparse(v) else np.asarray(v).ravel()
lam = mcl_obs[[c for c in ["IGLC1", "IGLC2", "IGLC3"] if c in mcl_obs]].sum(axis=1)
kap = mcl_obs["IGKC"] if "IGKC" in mcl_obs else 0
mcl_obs["kappa_minus_lambda"] = kap - lam

log("\n每位病人的 CCND1 / SOX11 / 輕鏈：")
log(mcl_obs.groupby("cnv_group")[["CCND1", "SOX11", "kappa_minus_lambda", "cnv_score"]]
    .mean().round(3).to_string())

mcl_obs.index = [i.rsplit("-", 1)[0] for i in mcl_obs.index]
mcl_obs = mcl_obs.reindex(mcl.obs_names)
for c in ["cnv_score", "cnv_leiden", "is_tumor_cnv", "CCND1", "SOX11", "kappa_minus_lambda"]:
    if c in mcl_obs:
        mcl.obs[c] = mcl_obs[c].values
mcl.obs["is_tumor"] = mcl.obs["is_tumor_cnv"].astype(bool)
log("\n最終腫瘤細胞數:", int(mcl.obs.is_tumor.sum()), "/", mcl.n_obs)
log(mcl.obs.groupby(["donor_id", "is_tumor"]).size().to_string())

mcl.write_h5ad(f"{PROC}/mcl_tumor_labeled.h5ad", compression="gzip")
mcl.obs.to_csv(f"{RES}/20_mcl_obs_labeled.csv")
log("DONE")
