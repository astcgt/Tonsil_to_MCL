#!/usr/bin/env python
"""
步驟 40：在外部公開世代上驗證「抗原呈現被關掉」是不是被套細胞淋巴瘤的普遍特徵。

為什麼需要這一步
----------------
主要分析只有兩位病人。外部世代雖然沒有染色質資料（公開資料中沒有任何
病人來源的被套細胞淋巴瘤多體學），但有足夠多的病人可以回答一個關鍵問題：
**CIITA / 第二類抗原呈現的壓制，是這兩位病人的特殊情況，還是普遍現象？**

設計（同檢體內配對，避開批次效應）
-----------------------------------
GSE303064 沒有做 CD19 分選，所以同一個檢體裡同時有惡性 B 細胞與正常 B 細胞。
在每個檢體內部：
  1. 分出 B 細胞（CD79A/CD79B/MS4A1 陽性，排除 T/NK/髓系）
  2. 在 B 細胞中以 CCND1 把惡性與非惡性分開（t(11;14) 讓惡性細胞 CCND1 高表現）
  3. 比較兩群的「CIITA–MHC-II 模組分數」
如果壓制是普遍的，每個檢體內惡性 B 都該低於非惡性 B。

輸出 results/40_external_validation.csv
"""
import os, sys, glob, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _compat  # noqa
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
EXT = f"{BASE}/data/external_mcl"
RES = f"{BASE}/results"
sc.settings.n_jobs = 48

MHC2 = ["CIITA", "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1",
        "HLA-DQA1", "HLA-DQB1", "CD74", "HLA-DMA", "HLA-DMB"]
MHC1 = ["HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "NLRC5"]
BCELL = ["CD79A", "CD79B", "MS4A1", "CD19", "BANK1", "PAX5"]
TCELL = ["CD3D", "CD3E", "CD2", "TRAC", "IL7R", "CD247", "TRBC2"]
MYE = ["LYZ", "CST3", "AIF1", "FCN1", "S100A8", "S100A9", "ITGAM", "CD14"]
OCT = ["POU2F2", "POU2AF1", "MEF2C", "PAX5", "EBF1"]


def log(*a):
    print("[40]", *a, flush=True)


def score(A, genes, name):
    g = [x for x in genes if x in A.var_names]
    if len(g) < 2:
        A.obs[name] = np.nan
        return 0
    sc.tl.score_genes(A, g, score_name=name, ctrl_size=min(50, max(10, len(g) * 5)))
    return len(g)


rows = []

# ---------------- GSE303064（11 位病人、未分選，含正常 B 細胞） ----------------
files = sorted(glob.glob(f"{EXT}/GSE303064/*.h5"))
log(f"GSE303064: {len(files)} 個樣本")
for f in files:
    sample = re.sub(r".*_(S\d+)\.filtered.*", r"\1", os.path.basename(f))
    try:
        A = sc.read_10x_h5(f)
    except Exception as e:
        log(f"  {sample} 讀取失敗: {e}"); continue
    A.var_names_make_unique()
    A.var["mt"] = A.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(A, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)
    keep = (A.obs.total_counts >= 800) & (A.obs.n_genes_by_counts >= 400) & (A.obs.pct_counts_mt < 20)
    A = A[keep].copy()
    if A.n_obs < 300:
        log(f"  {sample}: QC 後僅 {A.n_obs} 顆，跳過"); continue
    sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
    for nm, gl in [("B", BCELL), ("T", TCELL), ("My", MYE)]:
        score(A, gl, nm + "_score")
    isB = (A.obs.B_score > 0.15) & (A.obs.T_score < 0.1) & (A.obs.My_score < 0.1)
    B = A[isB].copy()
    if B.n_obs < 150:
        log(f"  {sample}: B 細胞僅 {B.n_obs} 顆，跳過"); continue

    cc = np.asarray(B[:, "CCND1"].X.todense()).ravel() if "CCND1" in B.var_names else np.zeros(B.n_obs)
    B.obs["CCND1"] = cc
    mal = cc > 0.5
    n_mal, n_nor = int(mal.sum()), int((~mal).sum())
    score(B, MHC2, "MHC2"); score(B, MHC1, "MHC1"); score(B, OCT, "OCT")
    if n_mal < 100 or n_nor < 30:
        log(f"  {sample}: 惡性 {n_mal} / 非惡性 {n_nor}（數量不足，仍記錄）")
    rows.append(dict(cohort="GSE303064", sample=sample, n_cells=A.n_obs, n_B=B.n_obs,
                     n_malignant=n_mal, n_normalB=n_nor,
                     CCND1_mal=float(cc[mal].mean()) if n_mal else np.nan,
                     CCND1_nor=float(cc[~mal].mean()) if n_nor else np.nan,
                     MHC2_mal=float(B.obs.MHC2[mal].mean()) if n_mal else np.nan,
                     MHC2_nor=float(B.obs.MHC2[~mal].mean()) if n_nor else np.nan,
                     MHC1_mal=float(B.obs.MHC1[mal].mean()) if n_mal else np.nan,
                     MHC1_nor=float(B.obs.MHC1[~mal].mean()) if n_nor else np.nan,
                     OCT_mal=float(B.obs.OCT[mal].mean()) if n_mal else np.nan,
                     OCT_nor=float(B.obs.OCT[~mal].mean()) if n_nor else np.nan))
    log(f"  {sample}: 細胞 {A.n_obs}, B {B.n_obs}, 惡性 {n_mal}, 非惡性 {n_nor}, "
        f"ΔMHC2 = {rows[-1]['MHC2_mal'] - rows[-1]['MHC2_nor']:+.3f}"
        if n_mal and n_nor else f"  {sample}: 記錄")
    del A, B

df = pd.DataFrame(rows)
if len(df):
    df["dMHC2"] = df.MHC2_mal - df.MHC2_nor
    df["dMHC1"] = df.MHC1_mal - df.MHC1_nor
    df["dOCT"] = df.OCT_mal - df.OCT_nor
    df.to_csv(f"{RES}/40_external_validation.csv", index=False)
    ok = df.dropna(subset=["dMHC2"])
    ok = ok[(ok.n_malignant >= 100) & (ok.n_normalB >= 30)]
    log(f"\n=== 可用於配對比較的樣本：{len(ok)} 個 ===")
    pd.set_option("display.width", 220)
    log("\n" + ok[["sample", "n_B", "n_malignant", "n_normalB",
                   "CCND1_mal", "CCND1_nor", "dMHC2", "dMHC1", "dOCT"]].round(3).to_string(index=False))
    if len(ok) >= 3:
        from scipy.stats import wilcoxon, ttest_1samp
        for col, nm in [("dMHC2", "第二類抗原呈現"), ("dMHC1", "第一類抗原呈現"), ("dOCT", "八聚體/B 細胞身分模組")]:
            v = ok[col].dropna().values
            if len(v) >= 3:
                try:
                    w = wilcoxon(v).pvalue
                except Exception:
                    w = np.nan
                t = ttest_1samp(v, 0).pvalue
                log(f"\n{nm}（惡性 − 非惡性）: 中位數 {np.median(v):+.3f}, "
                    f"{int((v<0).sum())}/{len(v)} 個樣本為負, Wilcoxon p={w:.4g}, t-test p={t:.4g}")
log("DONE")
