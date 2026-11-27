#!/usr/bin/env python
"""步驟 134：能不能在資料裡「模擬不富足的環境」？

問題
----
步驟 121 顯示 SLC7A11 在 1,178 株的中位 gene effect 是 **+0.086** —— 敲掉反而略好。
一個合理的解釋是：培養基（RPMI 約 200 µM 胱胺酸）遠比體內（40–60 µM）充裕，
而且培養皿裡沒有脂質過氧化壓力，所以這條抗鐵死亡的保險絲根本沒被點火。
如果這個解釋對，那在「胱胺酸相對不足」的細胞裡，SLC7A11 的依賴性應該要回來。

DepMap 沒有不同培養基的版本，但有一個**自然實驗**
------------------------------------------------
細胞自己有第二條取得半胱胺酸的路：硫轉移路徑（CBS → CTH，由甲硫胺酸合成半胱胺酸）。
硫轉移能力強的細胞，等於「自己會做」，對外部胱胺酸的依賴就低；
硫轉移能力弱的細胞，等於長期處在相對不足的狀態。
→ 用 CBS/CTH 表現把 1,178 株分層，看 SLC7A11 的依賴性是否隨之變化。
這是現有資料裡最接近「換培養基」的設計。

另外三個不需要新資料的讀出
--------------------------
A. 共依賴性（co-dependency）：跟 MAFG / SLC7A11 的依賴性輪廓最相關的基因是誰？
   這是 DepMap 判斷「這個基因實際上在哪條路上」的標準做法。
B. GPX4 依賴性：GPX4 是鐵死亡的主開關。如果 MCL 株對 GPX4 有依賴，
   代表鐵死亡這條軸在 MCL 裡即使在富足培養基中也是帶負載的。
C. AIFM2(FSP1) 表現 vs GPX4 依賴性：文獻上 FSP1 高的細胞對 GPX4 抑制劑有抗性。
   能不能在這份資料裡重現，是這套讀出有沒有效的陽性對照。

輸出 results/134_*
"""
import os
import json

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "8"

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
DM = f"{BASE}/data/depmap"
RES = f"{BASE}/results"

FOCUS = ["MAFG", "SLC7A11", "GPX4", "AIFM2", "SLC3A2", "GCLM", "GCLC"]
TRANSSULF = ["CBS", "CTH", "MAT2A"]          # 自己合成半胱胺酸的能力
EXPR_WANT = FOCUS + TRANSSULF + ["MAFF", "MAFK", "NFE2L2", "ACSL4", "LPCAT3"]
DEP_CUT = -0.5


def log(*a):
    print("[134]", *a, flush=True)


def load(path, genes=None):
    head = pd.read_csv(path, nrows=0)
    idx = head.columns[0]
    if genes is None:
        df = pd.read_csv(path, index_col=0)
    else:
        want = {c.split(" (")[0]: c for c in head.columns[1:]
                if c.split(" (")[0] in set(genes)}
        df = pd.read_csv(path, usecols=[idx] + list(want.values()), index_col=0)
    df.columns = [c.split(" (")[0] for c in df.columns]
    return df


def main():
    out = {}
    log("讀取 CRISPRGeneEffect 全表 …")
    EFF = load(f"{DM}/CRISPRGeneEffect.csv").astype(np.float32)
    log("gene effect:", EFF.shape)
    mod = pd.read_csv(f"{DM}/Model.csv", index_col=0)
    lin = mod["OncotreeLineage"].reindex(EFF.index)
    mcl_ids = [i for i in pd.read_csv(f"{DM}/MCL_cell_lines_DepMap_24Q4.csv")["ModelID"]
               if i in EFF.index]
    mcl_names = list(mod.loc[mcl_ids, "StrippedCellLineName"])
    log("MCL 株:", mcl_names)
    out["mcl_lines"] = mcl_names

    # ---------------- A. 共依賴性 ----------------
    log("\n=== A. 共依賴性：依賴性輪廓最相關的基因 ===")
    M = EFF.to_numpy()
    ok = np.isfinite(M).mean(0) > 0.95
    Mf = pd.DataFrame(M[:, ok], index=EFF.index, columns=EFF.columns[ok])
    Mf = Mf.fillna(Mf.mean())
    Z = (Mf - Mf.mean()) / (Mf.std() + 1e-9)
    Zv = Z.to_numpy(np.float32)
    n = Zv.shape[0]
    codep = {}
    for g in FOCUS:
        if g not in Z.columns:
            continue
        v = Z[g].to_numpy(np.float32)
        r = (Zv.T @ v) / n
        s = pd.Series(r, index=Z.columns).drop(index=[g]).sort_values(ascending=False)
        codep[g] = {"top_positive": s.head(12).round(3).to_dict(),
                    "top_negative": s.tail(6).round(3).to_dict()}
        log("  %-8s 最相關: %s" % (g, ", ".join("%s %.2f" % (k, v_)
                                                for k, v_ in list(s.head(8).items()))))
    out["codependency"] = codep

    # ---------------- B. GPX4 / 鐵死亡軸 ----------------
    log("\n=== B. 鐵死亡主開關 GPX4 ===")
    gp = EFF["GPX4"].dropna()
    rest = gp.drop(index=[i for i in mcl_ids if i in gp.index], errors="ignore")
    gm = gp.reindex(mcl_ids).dropna()
    lymph = gp[[i for i in gp.index if lin.get(i) == "Lymphoid" and i not in mcl_ids]]
    log("  全部 %d 株   中位 %+.3f   依賴比例 %.1f%%"
        % (len(gp), gp.median(), 100 * (gp < DEP_CUT).mean()))
    log("  MCL  %d 株   中位 %+.3f   依賴比例 %.0f%%   個別: %s"
        % (len(gm), gm.median(), 100 * (gm < DEP_CUT).mean(),
           ", ".join("%s %.2f" % (nm, v) for nm, v in
                     zip(mod.loc[gm.index, "StrippedCellLineName"], gm))))
    log("  非 MCL 淋巴系 %d 株 中位 %+.3f" % (len(lymph), lymph.median()))
    out["gpx4"] = {"n_all": int(len(gp)), "median_all": float(gp.median()),
                   "frac_dep_all": float((gp < DEP_CUT).mean()),
                   "mcl": {nm: float(v) for nm, v in
                           zip(mod.loc[gm.index, "StrippedCellLineName"], gm)},
                   "median_mcl": float(gm.median()),
                   "median_lymphoid_nonmcl": float(lymph.median()),
                   "mwu_p_mcl_more_dependent": float(
                       mannwhitneyu(gm, rest, alternative="less").pvalue) if len(gm) >= 3 else None}

    # ---------------- 表現量 ----------------
    EXP = load(f"{DM}/OmicsExpressionProteinCodingGenesTPMLogp1.csv", EXPR_WANT)
    log("\n表現量矩陣:", EXP.shape)
    common = EFF.index.intersection(EXP.index)
    log("同時有 CRISPR 與表現的株:", len(common))

    def corr(dep_gene, expr_gene):
        a = EFF.loc[common, dep_gene]
        b = EXP.loc[common, expr_gene]
        m = a.notna() & b.notna()
        r, p = spearmanr(a[m], b[m])
        return float(r), float(p), int(m.sum())

    # C. FSP1 高 → 對 GPX4 抑制較不敏感（陽性對照）
    r, p, nn = corr("GPX4", "AIFM2")
    log("\n=== C. 陽性對照：AIFM2(FSP1) 表現 vs GPX4 依賴性 ===")
    log("  Spearman rho = %+.3f  p = %.3g  (n=%d)" % (r, p, nn))
    log("  正值 = FSP1 越高的株，GPX4 敲除傷害越小 → 與文獻方向一致")
    out["fsp1_buffers_gpx4"] = {"rho": r, "p": p, "n": nn}

    # ---------------- D. 自然實驗：硫轉移能力分層 ----------------
    log("\n=== D. 自然實驗：用硫轉移能力（CBS+CTH）當『內生胱胺酸充裕度』 ===")
    ts = EXP.loc[common, [g for g in TRANSSULF if g in EXP.columns]]
    score = ts[["CBS", "CTH"]].mean(1) if {"CBS", "CTH"} <= set(ts.columns) else ts.mean(1)
    score = score.dropna()
    res = {}
    for dep_gene in ["SLC7A11", "SLC3A2", "GPX4", "MAFG", "GCLM"]:
        d = EFF.loc[score.index, dep_gene].dropna()
        s = score.reindex(d.index)
        r, p = spearmanr(s, d)
        q = s.quantile([1 / 3, 2 / 3])
        lo = d[s <= q.iloc[0]]          # 硫轉移弱 = 相對不足
        hi = d[s >= q.iloc[1]]          # 硫轉移強 = 自己會做
        u = mannwhitneyu(lo, hi, alternative="less")
        res[dep_gene] = {"spearman_rho": float(r), "p": float(p), "n": int(len(d)),
                         "median_low_transsulfuration": float(lo.median()),
                         "median_high_transsulfuration": float(hi.median()),
                         "delta": float(lo.median() - hi.median()),
                         "mwu_p_low_more_dependent": float(u.pvalue),
                         "frac_dep_low": float((lo < DEP_CUT).mean()),
                         "frac_dep_high": float((hi < DEP_CUT).mean())}
        log("  %-8s rho %+.3f (p %.2g) | 硫轉移弱 %+.3f vs 強 %+.3f → Δ %+.3f, p %.2g"
            % (dep_gene, r, p, lo.median(), hi.median(),
               lo.median() - hi.median(), u.pvalue))
        log("           依賴比例: 弱 %.1f%% vs 強 %.1f%%"
            % (100 * (lo < DEP_CUT).mean(), 100 * (hi < DEP_CUT).mean()))
    out["transsulfuration_stratification"] = res
    log("  → 負的 Δ 代表：硫轉移弱（＝內生胱胺酸相對不足）的細胞更依賴這個基因。")

    # 三株 MCL 落在哪
    mm = [i for i in mcl_ids if i in score.index]
    if mm:
        pct = [float((score < score[i]).mean()) for i in mm]
        log("\n  三株 MCL 的硫轉移分數百分位: %s"
            % ", ".join("%s %.0f%%" % (nm, 100 * p_) for nm, p_ in
                        zip(mod.loc[mm, "StrippedCellLineName"], pct)))
        out["mcl_transsulfuration_percentile"] = dict(
            zip(mod.loc[mm, "StrippedCellLineName"], pct))

    # SLC7A11 表現 vs 自身依賴（表現越高越上癮？）
    r, p, nn = corr("SLC7A11", "SLC7A11")
    log("\n  SLC7A11 表現 vs 自身依賴性: rho %+.3f, p %.3g (n=%d)" % (r, p, nn))
    out["slc7a11_expr_vs_dep"] = {"rho": r, "p": p, "n": nn}

    # 哪些株真的依賴 SLC7A11
    s7 = EFF["SLC7A11"].dropna()
    dep_lines = s7[s7 < DEP_CUT]
    byl = lin.reindex(dep_lines.index).value_counts()
    log("\n  真正依賴 SLC7A11 的株: %d / %d" % (len(dep_lines), len(s7)))
    log("  譜系分布 (前 6): %s" % byl.head(6).to_dict())
    out["slc7a11_dependent_lines"] = {"n": int(len(dep_lines)), "n_total": int(len(s7)),
                                      "by_lineage": {k: int(v) for k, v in byl.head(10).items()}}

    json.dump(out, open(f"{RES}/134_ferroptosis_context.json", "w"),
              indent=2, ensure_ascii=False)
    log("\nDONE → results/134_ferroptosis_context.json")


if __name__ == "__main__":
    main()
