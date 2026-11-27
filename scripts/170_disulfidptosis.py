#!/usr/bin/env python
"""步驟 170（B4）：disulfidptosis 的假說，在 DepMap 泛癌資料上站得住嗎？

機制（Liu et al., Nat Cell Biol 2023）
-------------------------------------
SLC7A11 把胱胺酸運進細胞，還原成半胱胺酸需要 **NADPH**，而 NADPH 主要來自
葡萄糖經磷酸戊糖途徑（PPP）。所以：
  SLC7A11 越高 → 需要還原的胱胺酸越多 → 對葡萄糖／PPP／NADPH 供應越依賴。
葡萄糖一斷，胱胺酸還原不完，細胞內堆積二硫鍵，**肌動蛋白骨架（WAVE 調節複合體）
被交聯而崩塌** → disulfidptosis。

可檢驗的預測（DepMap 是葡萄糖充裕的條件，所以效應應該**小但存在**）
------------------------------------------------------------------
P1  SLC7A11 表現越高的株，對葡萄糖轉運（SLC2A1/3）與 PPP（G6PD、PGD、TKT、
    TALDO1、RPIA、RPE）與 NADPH 供應（IDH1、ME1、NNT、PGD）**越依賴**
    → gene effect 與 SLC7A11 表現應為**負相關**。
P2  這個關係要比「隨機基因與 SLC7A11 表現的相關」強（用經驗 null 檢定），
    否則只是「表現高的株普遍比較脆弱」這種全域效應。
P3  **方向性對照**：與 SLC7A11 無關的代謝路徑（脂肪酸合成 FASN/ACACA、
    嘌呤合成 PPAT/GART）不該有同樣的關係。
P4  **陽性對照**：SLC3A2（SLC7A11 的必要伴侶）與 SLC7A11 表現的關係，
    方向應與 P1 相同且更強。
P5  WAVE 調節複合體（NCKAP1、WASF2、CYFIP1、ABI1/2、BRK1、RAC1）是執行者。
    文獻說敲掉它會**保護**細胞不死於 disulfidptosis。在葡萄糖充裕下不該看到，
    但如果 SLC7A11 高的株有基線二硫壓力，敲掉 WRC 的 gene effect 應該
    **比較不負**（＝正相關）。這是一個與 P1 **方向相反**的預測，
    所以它同時是機制特異性的檢定：若 P1 與 P5 都成立且方向相反，
    就很難用「全域脆弱度」解釋。
P6  譜系：這個關係在淋巴系是否更強？三株 MCL 落在哪？

輸出 results/170_*
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
RNG = np.random.default_rng(0)

SETS = {
    "glucose_transport": ["SLC2A1", "SLC2A3"],
    "PPP_oxidative": ["G6PD", "PGD", "TALDO1", "TKT", "RPIA", "RPE"],
    "NADPH_other": ["IDH1", "ME1", "NNT", "GSR"],
    "WRC_actin": ["NCKAP1", "WASF2", "CYFIP1", "ABI1", "ABI2", "BRK1", "RAC1"],
    "ctrl_fatty_acid": ["FASN", "ACACA", "SCD"],
    "ctrl_purine": ["PPAT", "GART", "PAICS", "ATIC"],
    "ctrl_panessential": ["RAN", "RPL13A", "POLR2A"],
    "positive_ctrl_partner": ["SLC3A2", "GCLC", "GCLM"],
}
DRIVER = "SLC7A11"
DEP_CUT = -0.5


def log(*a):
    print("[170]", *a, flush=True)


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
    log("讀取 DepMap …")
    EFF = pd.read_csv(f"{DM}/CRISPRGeneEffect.csv", index_col=0).astype(np.float32)
    EFF.columns = [c.split(" (")[0] for c in EFF.columns]
    EXP = load(f"{DM}/OmicsExpressionProteinCodingGenesTPMLogp1.csv",
               [DRIVER] + sum(SETS.values(), []) + ["SLC2A1"])
    mod = pd.read_csv(f"{DM}/Model.csv", index_col=0)
    lin = mod["OncotreeLineage"].reindex(EFF.index).fillna("Unknown")
    mcl_ids = [i for i in pd.read_csv(f"{DM}/MCL_cell_lines_DepMap_24Q4.csv")["ModelID"]
               if i in EFF.index]
    common = EFF.index.intersection(EXP.index)
    drv = EXP.loc[common, DRIVER].dropna()
    log("可用細胞株: %d | SLC7A11 表現中位 %.2f (log2 TPM+1)" % (len(drv), drv.median()))
    out = {"n_lines": int(len(drv)), "driver": DRIVER}

    # ---------- 經驗 null：隨機基因表現 vs 各依賴性的相關分布 ----------
    good = EFF.columns[(EFF.notna().mean(0) > 0.95) & (EFF.std(0) > 0.05)]
    log("可用於 null 的基因: %d" % len(good))

    def rho_dep_vs_driver(dep_gene, drv_vec):
        a = EFF.loc[drv_vec.index, dep_gene]
        m = a.notna()
        if m.sum() < 100:
            return np.nan
        return spearmanr(drv_vec[m], a[m]).statistic

    # null = 把 SLC7A11 表現換成 200 個隨機基因的表現，重算同一組相關
    EXPALL = None
    log("建立經驗 null（200 個隨機表現基因 × 同一組依賴性）…")
    cols = pd.read_csv(f"{DM}/OmicsExpressionProteinCodingGenesTPMLogp1.csv",
                       nrows=0).columns[1:]
    pick = RNG.choice(len(cols), size=260, replace=False)
    EXPALL = pd.read_csv(f"{DM}/OmicsExpressionProteinCodingGenesTPMLogp1.csv",
                         usecols=[0] + [int(i) + 1 for i in pick], index_col=0)
    EXPALL.columns = [c.split(" (")[0] for c in EXPALL.columns]
    EXPALL = EXPALL.loc[[i for i in common if i in EXPALL.index]]
    # 只留與 SLC7A11 表現分布相近（變異量相近）的基因，避免 null 太鬆
    sd0 = drv.std()
    keep = [c for c in EXPALL.columns
            if np.isfinite(EXPALL[c].std()) and 0.5 * sd0 < EXPALL[c].std() < 2 * sd0]
    keep = keep[:200]
    log("  null 基因（變異量與 SLC7A11 相近）: %d" % len(keep))

    rows = []
    for setname, genes in SETS.items():
        for g in genes:
            if g not in EFF.columns:
                continue
            r = rho_dep_vs_driver(g, drv)
            null = np.array([rho_dep_vs_driver(g, EXPALL[c].dropna()) for c in keep])
            null = null[np.isfinite(null)]
            # 雙尾經驗 p：|觀察| 在 null 的 |值| 分布中的位置
            emp_p = float((np.abs(null) >= abs(r)).mean()) if len(null) else np.nan
            z = (r - null.mean()) / (null.std() + 1e-9) if len(null) else np.nan
            rows.append(dict(set=setname, gene=g, rho=r, emp_p=emp_p, z_vs_null=z,
                             null_mean=float(null.mean()) if len(null) else np.nan,
                             dep_median=float(EFF[g].median()),
                             frac_dep=float((EFF[g] < DEP_CUT).mean())))
    T = pd.DataFrame(rows)
    T.to_csv(f"{RES}/170_slc7a11_vs_dependency.csv", index=False)

    pd.set_option("display.width", 220)
    log("\n=== SLC7A11 表現 vs 各基因依賴性的 Spearman ===")
    log("（負值 = SLC7A11 越高越依賴該基因，符合 disulfidptosis 預測）")
    log("%-22s %-10s %7s %8s %9s %10s %8s" %
        ("set", "gene", "rho", "emp_p", "z_vs_null", "dep_median", "frac_dep"))
    for s_ in SETS:
        sub = T[T.set == s_]
        for r in sub.itertuples():
            flag = " *" if (r.emp_p < 0.05) else ""
            log("%-22s %-10s %+7.3f %8.3f %+9.2f %10.3f %7.1f%%%s"
                % (r.set, r.gene, r.rho, r.emp_p, r.z_vs_null, r.dep_median,
                   100 * r.frac_dep, flag))
        if len(sub):
            log("%-22s %-10s %+7.3f  （該組平均）" % ("", "— mean —", sub.rho.mean()))

    # ---------- 組層級判讀 ----------
    verdict = {}
    for s_ in SETS:
        sub = T[T.set == s_].dropna(subset=["rho"])
        if not len(sub):
            continue
        verdict[s_] = {"n_genes": int(len(sub)), "mean_rho": float(sub.rho.mean()),
                       "n_sig": int((sub.emp_p < 0.05).sum()),
                       "mean_z_vs_null": float(sub.z_vs_null.mean())}
    out["by_set"] = verdict
    log("\n=== 組層級 ===")
    log("%-22s %7s %9s %12s" % ("set", "平均rho", "顯著數", "平均z vs null"))
    for k, v in verdict.items():
        log("%-22s %+7.3f %6d/%d %+12.2f"
            % (k, v["mean_rho"], v["n_sig"], v["n_genes"], v["mean_z_vs_null"]))

    log("\n判讀規則：")
    log("  P1 成立 = glucose_transport / PPP_oxidative / NADPH_other 平均 rho 明顯為負")
    log("  P3 成立 = ctrl_fatty_acid / ctrl_purine / ctrl_panessential 接近 0")
    log("  P5 成立 = WRC_actin 為**正**（方向相反）→ 這是機制特異性最強的證據")

    # ---------- P6 譜系 ----------
    log("\n=== P6. 淋巴系 vs 其他：葡萄糖/PPP 依賴與 SLC7A11 的關係 ===")
    key = [g for g in SETS["glucose_transport"] + SETS["PPP_oxidative"]
           if g in EFF.columns]
    lym = [i for i in drv.index if lin.get(i) == "Lymphoid"]
    oth = [i for i in drv.index if lin.get(i) != "Lymphoid"]
    byl = {}
    for g in key:
        rl = spearmanr(drv[lym], EFF.loc[lym, g]).statistic if len(lym) > 20 else np.nan
        ro = spearmanr(drv[oth], EFF.loc[oth, g]).statistic
        byl[g] = {"lymphoid": float(rl), "other": float(ro)}
        log("  %-8s 淋巴系 rho %+.3f (n=%d) | 其他 rho %+.3f (n=%d)"
            % (g, rl, len(lym), ro, len(oth)))
    out["lineage"] = byl

    log("\n=== 三株 MCL 的位置 ===")
    mm = [i for i in mcl_ids if i in drv.index]
    for i in mm:
        pct = float((drv < drv[i]).mean())
        g1 = EFF.loc[i, "SLC2A1"] if "SLC2A1" in EFF.columns else np.nan
        g6 = EFF.loc[i, "G6PD"] if "G6PD" in EFF.columns else np.nan
        log("  %-8s SLC7A11 表現 %.2f (第 %.0f 百分位) | SLC2A1 依賴 %+.2f | G6PD 依賴 %+.2f"
            % (mod.loc[i, "StrippedCellLineName"], drv[i], 100 * pct, g1, g6))
    out["mcl"] = {str(mod.loc[i, "StrippedCellLineName"]): {
        "slc7a11_expr": float(drv[i]),
        "slc7a11_expr_pct": float((drv < drv[i]).mean()),
        "SLC2A1_dep": float(EFF.loc[i, "SLC2A1"]) if "SLC2A1" in EFF.columns else None,
        "G6PD_dep": float(EFF.loc[i, "G6PD"]) if "G6PD" in EFF.columns else None}
        for i in mm}

    # ---------- 分層：SLC7A11 最高/最低三分位 ----------
    log("\n=== 分層檢定：SLC7A11 表現最高 vs 最低三分位 ===")
    q = drv.quantile([1 / 3, 2 / 3])
    hi = drv.index[drv >= q.iloc[1]]; lo = drv.index[drv <= q.iloc[0]]
    strat = {}
    for s_ in ("glucose_transport", "PPP_oxidative", "WRC_actin", "ctrl_purine"):
        for g in SETS[s_]:
            if g not in EFF.columns:
                continue
            a, b = EFF.loc[hi, g].dropna(), EFF.loc[lo, g].dropna()
            u = mannwhitneyu(a, b, alternative="less")
            strat[g] = {"set": s_, "median_high_SLC7A11": float(a.median()),
                        "median_low_SLC7A11": float(b.median()),
                        "delta": float(a.median() - b.median()),
                        "p_high_more_dependent": float(u.pvalue)}
            log("  %-22s %-8s 高 %+.3f vs 低 %+.3f → Δ %+.3f (p=%.3g)"
                % (s_, g, a.median(), b.median(), a.median() - b.median(), u.pvalue))
    out["stratified"] = strat

    json.dump(out, open(f"{RES}/170_disulfidptosis.json", "w"),
              indent=2, ensure_ascii=False)
    log("\nDONE → results/170_*")


if __name__ == "__main__":
    main()
