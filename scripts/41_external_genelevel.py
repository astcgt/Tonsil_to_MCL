#!/usr/bin/env python
"""
步驟 41：外部世代的「逐基因」配對驗證。

步驟 40 用模組分數發現：
  * 八聚體／B 細胞身分模組 → 17/19 樣本一致下降（可推廣）
  * 第二類抗原呈現模組     → 無一致趨勢（不可推廣）
但模組分數會被高表現基因（CD74、HLA-DRA）主導。
這一步改為逐基因檢定，看清楚到底哪些基因可以推廣、哪些不行。
"""
import os, sys, glob, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _compat  # noqa
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.stats import wilcoxon

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
EXT = f"{BASE}/data/external_mcl"
RES = f"{BASE}/results"
sc.settings.n_jobs = 48

GENES = ["CIITA", "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1",
         "CD74", "HLA-DMB", "NLRC5", "HLA-A", "HLA-B", "B2M",
         "POU2F2", "POU2AF1", "MEF2C", "MEF2A", "PAX5", "EBF1", "SPIB", "BANK1",
         "IRF4", "JUN", "JUNB", "FOS", "FOSL2", "BATF", "RELB", "NFKB2",
         "CCND1", "SOX11", "ROR1", "MSR1", "BCL2", "MYC", "TCL1A", "CR1", "FCRL1", "ATM"]
BCELL = ["CD79A", "CD79B", "MS4A1", "CD19", "BANK1", "PAX5"]
TCELL = ["CD3D", "CD3E", "CD2", "TRAC", "IL7R", "CD247", "TRBC2"]
MYE = ["LYZ", "CST3", "AIF1", "FCN1", "S100A8", "S100A9", "ITGAM", "CD14"]


def log(*a):
    print("[41]", *a, flush=True)


rows = []
for f in sorted(glob.glob(f"{EXT}/GSE303064/*.h5")):
    s = re.sub(r".*_(S\d+)\.filtered.*", r"\1", os.path.basename(f))
    A = sc.read_10x_h5(f)
    A.var_names_make_unique()
    A.var["mt"] = A.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(A, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)
    A = A[(A.obs.total_counts >= 800) & (A.obs.n_genes_by_counts >= 400)
          & (A.obs.pct_counts_mt < 20)].copy()
    if A.n_obs < 300:
        continue
    sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
    for nm, gl in [("B", BCELL), ("T", TCELL), ("My", MYE)]:
        gg = [x for x in gl if x in A.var_names]
        sc.tl.score_genes(A, gg, score_name=nm + "_s")
    B = A[(A.obs.B_s > 0.15) & (A.obs.T_s < 0.1) & (A.obs.My_s < 0.1)].copy()
    if B.n_obs < 150 or "CCND1" not in B.var_names:
        continue
    X = B.X.toarray() if sp.issparse(B.X) else np.asarray(B.X)
    gi = {k: i for i, k in enumerate(B.var_names)}
    cc = X[:, gi["CCND1"]]
    mal = cc > 0.5
    nor = cc == 0                     # 嚴格：CCND1 完全不表現才算非惡性
    if mal.sum() < 100 or nor.sum() < 50:
        log(f"  {s}: 惡性 {int(mal.sum())} / 嚴格非惡性 {int(nor.sum())} → 跳過")
        continue
    rec = dict(sample=s, n_mal=int(mal.sum()), n_nor=int(nor.sum()))
    for g in GENES:
        if g in gi:
            rec[g] = float(X[mal, gi[g]].mean() - X[nor, gi[g]].mean())
        else:
            rec[g] = np.nan
    rows.append(rec)
    log(f"  {s}: 惡性 {rec['n_mal']}, 嚴格非惡性 {rec['n_nor']}")
    del A, B, X

df = pd.DataFrame(rows)
df.to_csv(f"{RES}/41_external_genelevel.csv", index=False)
log(f"\n可用樣本 {len(df)} 個")

res = []
for g in GENES:
    if g not in df:
        continue
    v = df[g].dropna().values
    if len(v) < 5:
        continue
    try:
        p = wilcoxon(v).pvalue
    except Exception:
        p = np.nan
    res.append(dict(gene=g, median_delta=np.median(v),
                    n_neg=int((v < 0).sum()), n=len(v), pvalue=p))
r = pd.DataFrame(res)
r["padj"] = np.minimum(1, r.pvalue * len(r))
r = r.sort_values("median_delta")
r.to_csv(f"{RES}/41_external_gene_stats.csv", index=False)
pd.set_option("display.width", 200)
log("\n=== 外部 11 位病人：惡性 B − 同檢體內嚴格非惡性 B（逐基因）===")
log("（median_delta 為負 = 腫瘤中被壓制；n_neg/n = 幾個樣本方向一致）")
log("\n" + r.round(4).to_string(index=False))
log("DONE")
