#!/usr/bin/env python
"""
步驟 101：MAFG 矛盾的最後一塊 —— 「協同性」有獨立的染色質對應物嗎？

要解的問題
----------
MAFG 的元件在 MCL 裡**獲得相關性卻失去可及性**（−0.2 ~ −2.0）。
2.7 節 B 已經部分解釋了：這些元件是 CpG 島性質、本來就開著（偵測率 0.11、
甲基化 0.032），所以不可能再明顯「新開」。
但「協同性」這個說法到目前為止只有 RNA 層面的證據（peak–gene 相關性）。
這一步問：**不看 RNA，只看染色質，元件與 MAFG 啟動子的共開放度有變強嗎？**

為什麼這個設計避開了批次混淆
----------------------------
先前的 peak–gene 分析是「扁桃體圖譜的正常細胞」對「MCL 的腫瘤細胞」，
兩邊是不同資料集、不同批次。
這一步改成**同一個檢體內**：M413 的 790 顆非腫瘤 B 細胞 vs 同一管的腫瘤細胞。
同一次 10x 反應、同一個 ATAC 定量流程，所以批次、定量、peak 座標全部一致。
M102 只有 29 顆非腫瘤細胞，統計力不足，照實報告但不下結論。

兩個讀出
--------
A. **共開放度（log odds ratio）**：元件與 MAFG 啟動子在單細胞層面是否一起被偵測到。
   刻意不做 metacell —— metacell 的數量會直接影響相關係數大小，
   而兩組細胞數差 12 倍（9,974 vs 790），很容易做出假的「腫瘤相關性更強」。
   改用 odds ratio，並把腫瘤細胞**降抽到和正常組一樣的 790 顆**（重複 30 次取平均）。

B. **MAFG 表現分層**：多體學是同一顆細胞核同時量 RNA 與 ATAC，
   所以可以直接問「MAFG 表現高的細胞，這些元件是否更開」。
   每顆細胞的定序深度會同時影響兩邊，所以用**同窗格未重接線 peak 當對照**
   （它們承受完全相同的深度混淆），並額外報告偏相關（扣掉 log 片段總數）。
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _compat  # noqa
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
RES = f"{BASE}/results"
RNG = np.random.default_rng(0)
TSS = {"MAFG": ("chr17", 81927735), "SLC7A11": ("chr4", 138242349)}
NBOOT = 30


def log(*a):
    print("[101]", *a, flush=True)


# ------------------------------------------------------------------ 載入
lab = ad.read_h5ad(f"{BASE}/data/processed/mcl_tumor_labeled.h5ad")
rw = pd.read_csv(f"{RES}/52_rewiring_annotated.csv.gz")
pk = pd.read_csv(f"{BASE}/data/processed/atlas_peaks_sorted.csv")

rows_a, rows_b = [], []
for donor in ["M413", "M102"]:
    bc = pd.read_csv(f"{BASE}/data/processed/mcl_atac_atlaspeaks_{donor}_cells.csv",
                     header=None)[0].astype(str).values
    Z = sp.load_npz(f"{BASE}/data/processed/mcl_atac_atlaspeaks_{donor}.npz").tocsr()
    Z = Z[:len(bc)] if Z.shape[0] > len(bc) else Z
    Z = (Z > 0).astype(np.int8)                       # 二元偵測
    o = lab.obs.reindex(bc)
    keep = o.is_tumor.notna().values
    Z, bc, o = Z[keep], bc[keep], o[keep]
    tum = np.flatnonzero(o.is_tumor.astype(bool).values)
    nor = np.flatnonzero(~o.is_tumor.astype(bool).values)
    depth = np.asarray(Z.sum(axis=1)).ravel()
    log(f"\n{'='*70}\n{donor}: ATAC {Z.shape} | 腫瘤 {len(tum):,} / 非腫瘤 {len(nor):,}"
        f" | 片段深度中位 腫瘤 {np.median(depth[tum]):.0f} / 正常 {np.median(depth[nor]):.0f}")
    if len(nor) < 50:
        log("  非腫瘤細胞不足 50 顆，只報告描述性數字，不做檢定")

    for gene, (chrom, tss) in TSS.items():
        g = rw[(rw.gene == gene) & (rw.chrom == chrom)]
        dn = g[g.gained]
        cl = g[(~g.gained) & (~g.lost)]
        # 啟動子 peak：包含 TSS 的，否則最近的
        cand = pk[(pk.chrom == chrom) & (pk.start <= tss) & (pk.end >= tss)]
        if len(cand) == 0:
            sub = pk[pk.chrom == chrom].copy()
            sub["d"] = np.abs((sub.start + sub.end) // 2 - tss)
            cand = sub.nsmallest(1, "d")
        ppid = int(cand.pid.iloc[0])
        log(f"\n--- {gene}: 啟動子 peak pid={ppid} | 新生 {len(dn)} / 對照 {len(cl)} ---")

        P = np.asarray(Z[:, ppid].todense()).ravel()

        def logor(idx, a, b):
            """Haldane 校正後的 log odds ratio。"""
            n11 = float(((a[idx] == 1) & (b[idx] == 1)).sum()) + .5
            n10 = float(((a[idx] == 1) & (b[idx] == 0)).sum()) + .5
            n01 = float(((a[idx] == 0) & (b[idx] == 1)).sum()) + .5
            n00 = float(((a[idx] == 0) & (b[idx] == 0)).sum()) + .5
            return np.log((n11 * n00) / (n10 * n01))

        # ---------------- A. 共開放度
        for kind, sub in [("de_novo", dn), ("control", cl)]:
            for _, r in sub.iterrows():
                j = int(r.peak_id)
                if j == ppid:
                    continue
                E = np.asarray(Z[:, j].todense()).ravel()
                lo_n = logor(nor, E, P) if len(nor) >= 20 else np.nan
                # 腫瘤降抽到和正常組一樣的細胞數，重複取平均
                if len(nor) >= 20:
                    k = min(len(nor), len(tum))
                    vals = [logor(RNG.choice(tum, k, replace=False), E, P)
                            for _ in range(NBOOT)]
                    lo_t, lo_sd = float(np.mean(vals)), float(np.std(vals))
                else:
                    lo_t, lo_sd = logor(tum, E, P), np.nan
                rows_a.append(dict(
                    donor=donor, gene=gene, kind=kind, peak_id=j, dist=int(r.dist),
                    det_normal=round(float(E[nor].mean()), 5) if len(nor) else np.nan,
                    det_tumour=round(float(E[tum].mean()), 5),
                    logor_normal=round(lo_n, 4), logor_tumour=round(lo_t, 4),
                    logor_tumour_sd=round(lo_sd, 4) if lo_sd == lo_sd else np.nan,
                    delta=round(lo_t - lo_n, 4) if lo_n == lo_n else np.nan))

        # ---------------- B. 表現分層（只用腫瘤細胞，同一顆細胞核）
        if gene in lab.var_names:
            x = lab[:, gene].X
            x = np.asarray(x.todense()).ravel() if sp.issparse(x) else np.asarray(x).ravel()
            expr = pd.Series(x, index=lab.obs_names).reindex(bc).values
            tot = o.total_counts.values.astype(float)
            cpm = expr / np.maximum(tot, 1) * 1e4
            from scipy.stats import spearmanr, rankdata
            import statsmodels.api as sm
            ld = np.log1p(depth)

            # ⚠️ 踩過的陷阱（保留紀錄）：一開始用「秩殘差化」做偏相關是錯的。
            # 被偵測與否是二元變數，秩轉換後只有兩個值；把 log 深度迴歸掉之後，
            # 殘差幾乎只剩深度本身的結構，而兩個變數共用它，
            # 於是被 induce 出 +0.71 ~ +0.99 的假相關 —— 連對照組都一樣高、p 全是 0。
            # 正確做法是對二元結果直接做邏輯迴歸，讀表現項的係數。
            def logit_coef(y, x, c):
                """logit(P(peak 被偵測)) ~ 表現 + log 深度，回傳表現項的係數與 Wald p。"""
                X = sm.add_constant(np.c_[(x - x.mean()) / (x.std() + 1e-12),
                                          (c - c.mean()) / (c.std() + 1e-12)])
                try:
                    m = sm.Logit(y, X).fit(disp=0, maxiter=200)
                    return float(m.params[1]), float(m.pvalues[1])
                except Exception:
                    return np.nan, np.nan

            for kind, sub in [("de_novo", dn), ("control", cl)]:
                for _, r in sub.iterrows():
                    j = int(r.peak_id)
                    E = np.asarray(Z[:, j].todense()).ravel()[tum].astype(float)
                    if E.sum() < 10:
                        continue
                    s = spearmanr(E, cpm[tum])
                    b, pv = logit_coef(E, cpm[tum], ld[tum])
                    rows_b.append(dict(donor=donor, gene=gene, kind=kind, peak_id=j,
                                       dist=int(r.dist), n_pos=int(E.sum()),
                                       rho=round(float(s.statistic), 4), p=float(s.pvalue),
                                       beta_expr=round(b, 4) if b == b else np.nan,
                                       p_expr=pv))

A = pd.DataFrame(rows_a); A.to_csv(f"{RES}/101_coaccess.csv", index=False)
B = pd.DataFrame(rows_b); B.to_csv(f"{RES}/101_expr_strat.csv", index=False)

from scipy.stats import mannwhitneyu, wilcoxon
pd.set_option("display.width", 200)

log("\n" + "=" * 72)
log("=== A. 共開放度：元件與啟動子的 log odds ratio（腫瘤已降抽到同細胞數）===")
for (donor, gene), s in A.groupby(["donor", "gene"]):
    dn = s[(s.kind == "de_novo") & s.delta.notna()]
    cl = s[(s.kind == "control") & s.delta.notna()]
    log(f"\n{donor} / {gene}")
    log(f"  偵測率 新生 正常 {s[s.kind=='de_novo'].det_normal.median():.4f} → 腫瘤 "
        f"{s[s.kind=='de_novo'].det_tumour.median():.4f} | 對照 "
        f"{s[s.kind=='control'].det_normal.median():.4f} → "
        f"{s[s.kind=='control'].det_tumour.median():.4f}")
    if len(dn) < 3 or len(cl) < 5:
        log("  非腫瘤細胞不足，跳過檢定"); continue
    u = mannwhitneyu(dn.delta, cl.delta, alternative="greater")
    try:
        w = wilcoxon(dn.delta).pvalue
    except Exception:
        w = np.nan
    log(f"  Δ log-OR 新生 {dn.delta.median():+.3f} (n={len(dn)}) vs 對照 "
        f"{cl.delta.median():+.3f} (n={len(cl)}) | Mann-Whitney p = {u.pvalue:.4g}"
        f" | 新生自身 vs 0: Wilcoxon p = {w:.4g}")
    log(dn[["peak_id", "dist", "det_normal", "det_tumour", "logor_normal",
            "logor_tumour", "delta"]].sort_values("dist").to_string(index=False))

log("\n" + "=" * 72)
log("=== B. MAFG／SLC7A11 表現分層：元件可及性是否隨基因表現上升（腫瘤細胞內）===")
for (donor, gene), s in B.groupby(["donor", "gene"]):
    dn = s[(s.kind == "de_novo") & s.beta_expr.notna()]
    cl = s[(s.kind == "control") & s.beta_expr.notna()]
    if len(dn) < 3 or len(cl) < 5:
        log(f"\n{donor} / {gene}: 可用 peak 不足（新生 {len(dn)}、對照 {len(cl)}）"); continue
    u = mannwhitneyu(dn.beta_expr, cl.beta_expr, alternative="greater")
    log(f"\n{donor} / {gene}: 表現項係數中位 新生 {dn.beta_expr.median():+.4f} (n={len(dn)})"
        f" vs 對照 {cl.beta_expr.median():+.4f} (n={len(cl)}) | Mann-Whitney p = {u.pvalue:.4g}")
    log(f"  新生中 p<0.05 的有 {int((dn.p_expr<0.05).sum())}/{len(dn)}；"
        f"對照 {int((cl.p_expr<0.05).sum())}/{len(cl)}")
    log(dn[["peak_id", "dist", "n_pos", "rho", "beta_expr", "p_expr"]]
        .sort_values("beta_expr", ascending=False).head(8)
        .to_string(index=False, float_format=lambda x: f"{x:.4g}"))
log("\nDONE")
