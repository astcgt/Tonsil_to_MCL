#!/usr/bin/env python
"""步驟 231：步驟 41／230 的「惡性 vs 嚴格非惡性」切法帶有定序深度的混淆。

為什麼需要這一支
----------------
步驟 230 的事先註冊放了一個無所不在的陰性對照 `RPL13A`（核糖體蛋白，每顆細胞都要）。
它**失敗了**：median delta +0.428、19/19 個樣本同向、padj 6.1×10⁻⁵。
一個每顆細胞都在表現的基因不應該在「惡性 vs 非惡性」之間系統性地差這麼多，
所以差的不是生物，是切法。

被檢驗的機制
------------
`nor = (CCND1 == 0)` 這個條件不只挑出非惡性細胞，**也順便挑出低深度／高 dropout 的細胞**
—— 因為在一顆讀得淺的細胞裡，任何中低表現基因都更可能剛好是 0，CCND1 也是。
CP10K 正規化**修不掉這一點**：縮放改變尺度，但不會把 dropout 造成的零補回來，
所以淺細胞的 log1p(CP10K) 對每個基因都系統性偏低。

可否證的預測
------------
P1 惡性組的 UMI 中位數 **高於**嚴格非惡性組，且在**每一個樣本**都如此。
P2 跨樣本，`UMI 比值` 與 `RPL13A delta` **正相關**（深度差越大，假訊號越大）。
P1 與 P2 都成立 → 這個切法的任何**正向** delta 都不可直接解讀。

對既有結論的影響（事先寫明，避免事後挑選）
-----------------------------------------
偏誤方向是**讓惡性組的每個基因都偏高**。所以：
  * 步驟 41 的「B 細胞身分／八聚體模組在惡性細胞中下降（17/19）」是**逆著偏誤**成立的，
    屬於保守，結論不受威脅；
  * 任何**正向** delta（包含步驟 230 盤裡的 GPX4 +0.187、SLC3A2 +0.062）**不可採信**；
  * 步驟 41 的「第二類抗原呈現無一致趨勢」有可能是真實下降被偏誤抵銷，應標為未決。

輸出 results/231_depth_confound.json
"""
import os
import sys
import glob
import json
import re

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _compat  # noqa

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.stats import spearmanr, wilcoxon

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
EXT = f"{BASE}/data/external_mcl/GSE303064"
RES = f"{BASE}/results"
sc.settings.n_jobs = 16

BCELL = ["CD79A", "CD79B", "MS4A1", "CD19", "BANK1", "PAX5"]
TCELL = ["CD3D", "CD3E", "CD2", "TRAC", "IL7R", "CD247", "TRBC2"]
MYE = ["LYZ", "CST3", "AIF1", "FCN1", "S100A8", "S100A9", "ITGAM", "CD14"]
UBIQ = ["RPL13A", "RPS6", "ACTB", "GAPDH", "TMSB4X"]   # 無所不在的對照基因


def log(*a):
    print("[231]", *a, flush=True)


def main():
    rows = []
    for f in sorted(glob.glob(f"{EXT}/*.h5")):
        s = re.sub(r".*_(S\d+)\.filtered.*", r"\1", os.path.basename(f))
        A = sc.read_10x_h5(f)
        A.var_names_make_unique()
        A.var["mt"] = A.var_names.str.upper().str.startswith("MT-")
        sc.pp.calculate_qc_metrics(A, qc_vars=["mt"], inplace=True,
                                   percent_top=None, log1p=False)
        A = A[(A.obs.total_counts >= 800) & (A.obs.n_genes_by_counts >= 400)
              & (A.obs.pct_counts_mt < 20)].copy()
        if A.n_obs < 300:
            continue
        tot = A.obs.total_counts.to_numpy().copy()
        ng = A.obs.n_genes_by_counts.to_numpy().copy()
        sc.pp.normalize_total(A, target_sum=1e4)
        sc.pp.log1p(A)
        for nm, gl in [("B", BCELL), ("T", TCELL), ("My", MYE)]:
            sc.tl.score_genes(A, [x for x in gl if x in A.var_names],
                              score_name=nm + "_s")
        keep = ((A.obs.B_s > 0.15) & (A.obs.T_s < 0.1) & (A.obs.My_s < 0.1)).to_numpy()
        B = A[keep].copy()
        tot, ng = tot[keep], ng[keep]
        if B.n_obs < 150 or "CCND1" not in B.var_names:
            continue
        X = B.X.toarray() if sp.issparse(B.X) else np.asarray(B.X)
        gi = {k: i for i, k in enumerate(B.var_names)}
        cc = X[:, gi["CCND1"]]
        mal, nor = cc > 0.5, cc == 0
        if mal.sum() < 100 or nor.sum() < 50:
            continue
        rec = dict(sample=s, n_mal=int(mal.sum()), n_nor=int(nor.sum()),
                   umi_mal=float(np.median(tot[mal])),
                   umi_nor=float(np.median(tot[nor])),
                   genes_mal=float(np.median(ng[mal])),
                   genes_nor=float(np.median(ng[nor])))
        for g in UBIQ:
            rec["delta_" + g] = (float(X[mal, gi[g]].mean() - X[nor, gi[g]].mean())
                                 if g in gi else np.nan)
        rows.append(rec)
        log("  %s：UMI 惡性 %.0f vs 非惡性 %.0f（比值 %.2f）| RPL13A delta %+.4f"
            % (s, rec["umi_mal"], rec["umi_nor"],
               rec["umi_mal"] / rec["umi_nor"], rec["delta_RPL13A"]))
        del A, B, X

    d = pd.DataFrame(rows)
    d["umi_ratio"] = d.umi_mal / d.umi_nor
    d["genes_ratio"] = d.genes_mal / d.genes_nor
    d.to_csv(f"{RES}/231_depth_confound.csv", index=False)
    n = len(d)

    p1 = bool((d.umi_ratio > 1).all())
    rs = spearmanr(d.umi_ratio, d.delta_RPL13A)
    p2 = bool(rs.statistic > 0 and rs.pvalue < 0.05)

    out = {"cohort": "GSE303064", "n_samples": int(n),
           "P1_malignant_deeper": {
               "rule": "每一個樣本的惡性組 UMI 中位數都高於嚴格非惡性組",
               "umi_ratio_median": float(d.umi_ratio.median()),
               "umi_ratio_range": [float(d.umi_ratio.min()), float(d.umi_ratio.max())],
               "genes_ratio_median": float(d.genes_ratio.median()),
               "n_samples_above_1": int((d.umi_ratio > 1).sum()),
               "pass": p1},
           "P2_confound_scales_with_depth_gap": {
               "rule": "跨樣本 UMI 比值與 RPL13A delta 正相關且 p < 0.05",
               "spearman_r": float(rs.statistic), "pvalue": float(rs.pvalue),
               "pass": p2},
           "ubiquitous_gene_deltas": {
               g: {"median": float(d["delta_" + g].median()),
                   "n_positive": int((d["delta_" + g] > 0).sum()), "n": int(n),
                   "wilcoxon_p": float(wilcoxon(d["delta_" + g].dropna()).pvalue)}
               for g in UBIQ if "delta_" + g in d and d["delta_" + g].notna().sum() >= 5},
           "consequences": {
               "bias_direction": "讓惡性組的每個基因都偏高",
               "step41_identity_module_down": "逆著偏誤成立，屬保守，結論不受威脅",
               "any_positive_delta": "不可直接解讀，含步驟 230 的 GPX4 與 SLC3A2",
               "step41_MHCII_no_trend": "可能是真實下降被偏誤抵銷，標為未決"}}
    out["verdict"] = ("confound_confirmed" if (p1 and p2) else
                      "confound_not_supported")

    log("\n=== P1：惡性組是不是讀得比較深 ===")
    log("  UMI 比值中位 %.2f（範圍 %.2f–%.2f）| %d/%d 個樣本 > 1"
        % (d.umi_ratio.median(), d.umi_ratio.min(), d.umi_ratio.max(),
           int((d.umi_ratio > 1).sum()), n))
    log("  偵測基因數比值中位 %.2f" % d.genes_ratio.median())
    log("  P1 %s" % ("通過" if p1 else "未通過"))
    log("\n=== P2：深度差越大，假訊號越大？ ===")
    log("  Spearman r = %.3f, p = %.3g → P2 %s"
        % (rs.statistic, rs.pvalue, "通過" if p2 else "未通過"))
    log("\n=== 無所不在的對照基因，在這個切法下的 delta ===")
    for g, v in out["ubiquitous_gene_deltas"].items():
        log("  %-8s median %+8.4f | %2d/%2d 同向 | Wilcoxon p %.3g"
            % (g, v["median"], v["n_positive"], v["n"], v["wilcoxon_p"]))
    log("\n*** 判定：%s ***" % out["verdict"])
    if out["verdict"] == "confound_confirmed":
        log("    `CCND1 == 0` 當非惡性組會順便挑出低深度細胞，CP10K 修不掉。")
        log("    這個切法的任何正向 delta 都不可直接解讀。")

    json.dump(out, open(f"{RES}/231_depth_confound.json", "w"),
              indent=2, ensure_ascii=False)
    log("\n→ results/231_depth_confound.json、231_depth_confound.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
