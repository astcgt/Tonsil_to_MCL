#!/usr/bin/env python
"""步驟 203：用 DepMap 否證一個機制預測。

預測從哪裡來
------------
步驟 202 的化學計量模型（Human-GEM，四道驗證門檻全過）給出一個很硬的結果：
把轉硫途徑（甲硫胺酸 → 半胱胺酸）關掉之後，敲掉 `SLC7A11` 會讓穀胱甘肽產能
**在每一個胱胺酸濃度下都歸零**——xCT 變成唯一入口，沒有旁路。

而本專案自己的表現資料說，轉硫途徑在 B 細胞與 MCL 就是關著的：
  CBS   MCL 0.0003 / 正常 B 0.0002（46 個正常狀態中 0 個 > 0.1）
  CBSL  0.0001 / 0.0000（0/46）
  CTH   0.0090 / 0.0112（1/46）

這剛好可以解釋步驟 135 那個還沒有機制的觀察：依賴 `SLC7A11` 的 18 株裡
12 株是淋巴系、而且每一株都是 B 細胞腫瘤（OR 26.6、經驗 p = 0.021）。
如果機制是對的，那不是「B 細胞剛好這樣」，而是「沒有轉硫途徑的細胞才會這樣」。

可否證的預測
------------
P1  全庫 1,178 株裡，轉硫途徑表現越低，`SLC7A11` 的 gene effect 越負。
    也就是 corr(轉硫表現, SLC7A11 gene effect) 應該是**正的**。
P2  這個關係必須在**非淋巴系**的株裡也成立。否則它只是譜系的代名詞，
    不是機制 —— 淋巴系同時低 CBS 又高度依賴 xCT，光看全庫無法分辨這兩件事。
P3  效應量要大於「隨便挑一個基因的表現也會有的相關」（經驗虛無分佈）。
P4  陰性對照：轉硫表現不該同樣地預測一個無關的必需基因
    （這裡用 `RPL13A`，核糖體蛋白，所有株都依賴）。

任何一項不過都要寫進結論，不是只報過的那幾項。

輸出 results/203_*
"""
import os
import json
import warnings

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "8"
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from scipy import stats

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
RES = f"{BASE}/results"
DM = f"{BASE}/data/depmap"

TS_GENES = ["CBS", "CBSL", "CTH"]
TARGET = "SLC7A11"
NEG_CTRL = "RPL13A"
EXTRA = ["SLC3A2", "GCLC", "GCLM", "MAT2A", "SLC7A11"]
N_NULL = 2000
RNG = np.random.default_rng(0)


def log(*a):
    print("[203]", *a, flush=True)


def gene_cols(path, genes):
    """只讀需要的欄位：DepMap 的表頭是 `SYMBOL (ENTREZ)`。"""
    hdr = pd.read_csv(path, nrows=0)
    first = hdr.columns[0]
    want = {}
    for c in hdr.columns[1:]:
        sym = c.split(" (")[0]
        if sym in genes:
            want[sym] = c
    return first, want, list(hdr.columns)


def main():
    eff_path = f"{DM}/CRISPRGeneEffect.csv"
    exp_path = f"{DM}/OmicsExpressionProteinCodingGenesTPMLogp1.csv"

    need_eff = set([TARGET, NEG_CTRL] + EXTRA)
    f1, eff_cols, _ = gene_cols(eff_path, need_eff)
    EFF = pd.read_csv(eff_path, usecols=[f1] + list(eff_cols.values()))
    EFF = EFF.rename(columns={v: k for k, v in eff_cols.items()}).set_index(f1)
    log("CRISPR gene effect：%d 株 × %d 基因" % EFF.shape)

    f2, exp_cols, all_exp = gene_cols(exp_path, set(TS_GENES + EXTRA))
    EXP = pd.read_csv(exp_path, usecols=[f2] + list(exp_cols.values()))
    EXP = EXP.rename(columns={v: k for k, v in exp_cols.items()}).set_index(f2)
    log("表現：%d 株 × %d 基因（找到 %s）"
        % (EXP.shape[0], EXP.shape[1], ", ".join(sorted(exp_cols))))

    mod = pd.read_csv(f"{DM}/Model.csv").set_index("ModelID")
    lin = mod["OncotreeLineage"] if "OncotreeLineage" in mod.columns else mod["lineage"]

    ts_avail = [g for g in TS_GENES if g in EXP.columns]
    df = pd.DataFrame({
        "slc7a11_eff": EFF[TARGET],
        "negctrl_eff": EFF[NEG_CTRL] if NEG_CTRL in EFF.columns else np.nan,
    })
    for g in ts_avail:
        df[f"exp_{g}"] = EXP[g]
    df["ts_score"] = df[[f"exp_{g}" for g in ts_avail]].mean(axis=1)
    df["lineage"] = lin.reindex(df.index)
    df["lymphoid"] = df.lineage.astype(str).str.contains("Lymph", case=False)
    df = df.dropna(subset=["slc7a11_eff", "ts_score"])
    log("可用於分析的株數：%d（淋巴系 %d）"
        % (len(df), int(df.lymphoid.sum())))
    df.to_csv(f"{RES}/203_lines.csv")

    out = {"n_lines": int(len(df)), "n_lymphoid": int(df.lymphoid.sum()),
           "transsulfuration_genes_found": ts_avail}

    # P1：全庫相關
    r, p = stats.spearmanr(df.ts_score, df.slc7a11_eff)
    out["P1_all_lines"] = {"spearman_r": float(r), "p": float(p),
                           "direction_expected": "正（表現越高越不依賴）",
                           "pass": bool(r > 0 and p < 0.05)}
    log("P1 全庫：Spearman r = %+.3f，p = %.2e → %s"
        % (r, p, "過" if out["P1_all_lines"]["pass"] else "不過"))

    # P2：只看非淋巴系
    nl = df[~df.lymphoid]
    r2, p2 = stats.spearmanr(nl.ts_score, nl.slc7a11_eff)
    out["P2_non_lymphoid_only"] = {"n": int(len(nl)), "spearman_r": float(r2),
                                   "p": float(p2),
                                   "pass": bool(r2 > 0 and p2 < 0.05)}
    log("P2 非淋巴系（n = %d）：r = %+.3f，p = %.2e → %s"
        % (len(nl), r2, p2, "過" if out["P2_non_lymphoid_only"]["pass"] else "不過"))
    ly = df[df.lymphoid]
    r3, p3 = stats.spearmanr(ly.ts_score, ly.slc7a11_eff)
    out["P2b_lymphoid_only"] = {"n": int(len(ly)), "spearman_r": float(r3),
                                "p": float(p3)}
    log("   （淋巴系內部 n = %d：r = %+.3f，p = %.3f）" % (len(ly), r3, p3))
    log("   轉硫分數：淋巴系中位 %.3f vs 非淋巴系 %.3f（Mann-Whitney p = %.2e）"
        % (ly.ts_score.median(), nl.ts_score.median(),
           stats.mannwhitneyu(ly.ts_score, nl.ts_score,
                              alternative="less").pvalue))
    out["ts_score_lymphoid_vs_other"] = {
        "median_lymphoid": float(ly.ts_score.median()),
        "median_other": float(nl.ts_score.median()),
        "mannwhitney_p_lymphoid_lower": float(
            stats.mannwhitneyu(ly.ts_score, nl.ts_score,
                               alternative="less").pvalue)}

    # P3：經驗虛無 —— 隨機基因的表現對 SLC7A11 依賴性的相關
    hdr = pd.read_csv(exp_path, nrows=0)
    pool = [c for c in hdr.columns[1:]]
    pick = RNG.choice(len(pool), size=min(N_NULL, len(pool)), replace=False)
    cols = [pool[i] for i in pick]
    log("P3 讀 %d 個隨機基因的表現建虛無分佈（分批）…" % len(cols))
    null_r = []
    step = 250
    for i in range(0, len(cols), step):
        sub = pd.read_csv(exp_path, usecols=[hdr.columns[0]] + cols[i:i + step])
        sub = sub.set_index(hdr.columns[0]).reindex(df.index)
        for c in sub.columns:
            v = sub[c]
            ok = v.notna() & df.slc7a11_eff.notna()
            if ok.sum() > 200:
                rr, _ = stats.spearmanr(v[ok], df.slc7a11_eff[ok])
                if rr == rr:
                    null_r.append(float(rr))
    null_r = np.array(null_r)
    emp_p = float((null_r >= r).mean())
    out["P3_empirical_null"] = {
        "n_random_genes": int(len(null_r)),
        "null_mean_r": float(null_r.mean()), "null_sd_r": float(null_r.std()),
        "null_q95": float(np.quantile(null_r, 0.95)),
        "observed_r": float(r), "empirical_p": emp_p,
        "pass": bool(emp_p < 0.05)}
    log("P3 虛無：%d 個隨機基因 r = %.3f ± %.3f（95 分位 %.3f）；"
        "觀察到的 %.3f → 經驗 p = %.4f → %s"
        % (len(null_r), null_r.mean(), null_r.std(),
           np.quantile(null_r, 0.95), r, emp_p,
           "過" if out["P3_empirical_null"]["pass"] else "不過"))

    # P4：陰性對照
    if df.negctrl_eff.notna().sum() > 100:
        r4, p4 = stats.spearmanr(df.ts_score, df.negctrl_eff,
                                 nan_policy="omit")
        out["P4_negative_control"] = {
            "gene": NEG_CTRL, "spearman_r": float(r4), "p": float(p4),
            "pass": bool(abs(r4) < abs(r))}
        log("P4 陰性對照 %s：r = %+.3f，p = %.2e → %s"
            % (NEG_CTRL, r4, p4,
               "過（效應明顯小於標的）" if out["P4_negative_control"]["pass"]
               else "不過"))

    # 附帶：18 株依賴者 vs 其餘的轉硫分數
    dep = df[df.slc7a11_eff < -0.5]
    oth = df[df.slc7a11_eff >= -0.5]
    u = stats.mannwhitneyu(dep.ts_score, oth.ts_score, alternative="less")
    out["dependent_vs_rest"] = {
        "n_dependent": int(len(dep)), "n_rest": int(len(oth)),
        "median_ts_dependent": float(dep.ts_score.median()),
        "median_ts_rest": float(oth.ts_score.median()),
        "mannwhitney_p": float(u.pvalue)}
    log("依賴者（gene effect < −0.5，n = %d）轉硫分數中位 %.3f "
        "vs 其餘（n = %d）%.3f，Mann-Whitney p = %.2e"
        % (len(dep), dep.ts_score.median(), len(oth),
           oth.ts_score.median(), u.pvalue))

    passed = [k for k, v in out.items()
              if isinstance(v, dict) and v.get("pass") is True]
    failed = [k for k, v in out.items()
              if isinstance(v, dict) and v.get("pass") is False]
    out["summary"] = {"passed": passed, "failed": failed}
    log("\n通過：%s" % (", ".join(passed) or "（無）"))
    log("未通過：%s" % (", ".join(failed) or "（無）"))
    json.dump(out, open(f"{RES}/203_transsulfuration.json", "w"),
              indent=2, ensure_ascii=False)
    log("DONE → results/203_*")


if __name__ == "__main__":
    main()
