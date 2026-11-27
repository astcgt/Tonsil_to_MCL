#!/usr/bin/env python
"""步驟 135：步驟 134 冒出來的兩件事，各給一個嚴謹的檢定。

觀察一：全部 1,178 株裡只有 18 株真的依賴 SLC7A11，而其中 12 株是淋巴系。
        → 這不是「沒有細胞需要它」，而是「只有淋巴系需要它」。
        檢定：Fisher 精確檢定，並且用**依賴率配對的隨機基因**當背景，
        排除「淋巴系本來就對什麼都比較敏感」這個解釋。

觀察二：SLC7A11 / SLC3A2 / GCLC / GCLM / AIFM2 五個基因的依賴性輪廓互相高度相關，
        而且都和 B 細胞身分基因（CD79A、EBF1、PAX5、FOXO1）相關。
        → 檢定：這五個基因的兩兩相關平均值，是否高於隨機五基因組的分布。

輸出 results/135_*
"""
import os
import json

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "8"

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
DM = f"{BASE}/data/depmap"
RES = f"{BASE}/results"
DEP_CUT = -0.5
MODULE = ["SLC7A11", "SLC3A2", "GCLC", "GCLM", "AIFM2"]
BCELL = ["CD79A", "CD79B", "EBF1", "PAX5", "FOXO1", "SYK", "BLNK", "IRF4"]
RNG = np.random.default_rng(0)


def log(*a):
    print("[135]", *a, flush=True)


def main():
    EFF = pd.read_csv(f"{DM}/CRISPRGeneEffect.csv", index_col=0).astype(np.float32)
    EFF.columns = [c.split(" (")[0] for c in EFF.columns]
    mod = pd.read_csv(f"{DM}/Model.csv", index_col=0)
    lin = mod["OncotreeLineage"].reindex(EFF.index).fillna("Unknown")
    sub = mod["OncotreeSubtype"].reindex(EFF.index).fillna("")
    is_lymph = (lin == "Lymphoid").to_numpy()
    log("株數 %d | 淋巴系 %d" % (len(EFF), int(is_lymph.sum())))
    out = {"n_lines": int(len(EFF)), "n_lymphoid": int(is_lymph.sum())}

    # ---------- 觀察一：譜系富集 ----------
    dep_frac_all = (EFF < DEP_CUT).mean(0)
    res = {}
    log("\n=== 依賴這個基因的細胞株裡，淋巴系佔多少？ ===")
    log("%-9s %6s %6s %7s %7s %8s %10s %10s" %
        ("gene", "n_dep", "n_lym", "lym%", "OR", "p", "emp_p", "dep_frac"))
    for g in MODULE + ["MAFG", "GPX4", "NFE2L2", "BACH1"]:
        if g not in EFF.columns:
            continue
        v = EFF[g]
        d = (v < DEP_CUT).to_numpy() & v.notna().to_numpy()
        a = int((d & is_lymph).sum())
        b = int((d & ~is_lymph).sum())
        c = int((~d & is_lymph).sum())
        e = int((~d & ~is_lymph).sum())
        orr, pv = fisher_exact([[a, b], [c, e]], alternative="greater")

        # 依賴率配對的隨機基因背景：抽 2,000 個依賴率相近的基因，
        # 看它們的「依賴株中淋巴系比例」分布
        f = float(d.mean())
        near = dep_frac_all[(dep_frac_all > f * 0.6) & (dep_frac_all < f * 1.6)
                            & (dep_frac_all > 0)].index
        near = [x for x in near if x != g]
        emp = []
        if len(near) >= 50:
            for gg in RNG.choice(near, size=min(2000, len(near)), replace=False):
                dd = (EFF[gg] < DEP_CUT).to_numpy() & EFF[gg].notna().to_numpy()
                if dd.sum() >= 3:
                    emp.append((dd & is_lymph).sum() / dd.sum())
        obs = a / max(a + b, 1)
        emp_p = float((np.array(emp) >= obs).mean()) if emp else np.nan
        res[g] = {"n_dependent": a + b, "n_lymphoid_dependent": a,
                  "frac_lymphoid_among_dependent": obs,
                  "odds_ratio": float(orr), "fisher_p": float(pv),
                  "empirical_p_vs_rate_matched_genes": emp_p,
                  "n_matched_genes": len(emp), "dep_frac": f}
        log("%-9s %6d %6d %6.1f%% %7.1f %8.2g %10s %10.3f"
            % (g, a + b, a, 100 * obs, orr, pv,
               ("%.4f" % emp_p) if emp else "-", f))
    out["lineage_enrichment"] = res

    # 淋巴系內部：B 系 vs 其他
    s7 = EFF["SLC7A11"]
    dl = s7[(s7 < DEP_CUT)].index
    log("\n=== 依賴 SLC7A11 的 %d 株是誰 ===" % len(dl))
    for i in dl:
        log("  %-14s %-10s %-34s %+.2f"
            % (mod.loc[i, "StrippedCellLineName"], lin[i], sub[i][:34], s7[i]))
    out["slc7a11_dependent_lines"] = [
        {"line": str(mod.loc[i, "StrippedCellLineName"]), "lineage": str(lin[i]),
         "subtype": str(sub[i]), "gene_effect": float(s7[i])} for i in dl]

    # ---------- 觀察二：模組內聚性 ----------
    log("\n=== 五個路徑基因的依賴性輪廓，互相有多像？ ===")
    M = EFF[[g for g in MODULE if g in EFF.columns]].dropna()
    C = M.corr()
    iu = np.triu_indices(len(C), 1)
    obs = float(C.to_numpy()[iu].mean())
    log("\n" + C.round(3).to_string())
    log("\n模組內兩兩相關平均 = %.3f" % obs)

    good = EFF.columns[(EFF.notna().mean(0) > 0.95) & (EFF.std(0) > 0.05)]
    E = EFF[good].fillna(EFF[good].mean())
    null = []
    for _ in range(2000):
        gs = RNG.choice(len(good), size=len(MODULE), replace=False)
        c = np.corrcoef(E.iloc[:, gs].to_numpy().T)
        null.append(c[np.triu_indices(len(MODULE), 1)].mean())
    null = np.array(null)
    pv = float((null >= obs).mean())
    log("隨機五基因組的分布: 平均 %.3f, 95 百分位 %.3f  → 經驗 p = %.4f"
        % (null.mean(), np.percentile(null, 95), pv))
    out["module_coherence"] = {"genes": MODULE, "mean_pairwise_r": obs,
                               "null_mean": float(null.mean()),
                               "null_p95": float(np.percentile(null, 95)),
                               "empirical_p": pv,
                               "matrix": json.loads(C.round(4).to_json())}

    # 與 B 細胞身分基因的相關
    log("\n=== 路徑基因 vs B 細胞身分基因的依賴性相關 ===")
    bc = [g for g in BCELL if g in EFF.columns]
    X = EFF[[g for g in MODULE if g in EFF.columns] + bc].dropna()
    CB = X.corr().loc[[g for g in MODULE if g in X.columns], bc]
    log("\n" + CB.round(3).to_string())
    out["module_vs_bcell_identity"] = json.loads(CB.round(4).to_json())
    log("\n模組 × B 細胞身分 的平均相關 = %.3f" % float(CB.to_numpy().mean()))
    out["module_vs_bcell_mean_r"] = float(CB.to_numpy().mean())

    json.dump(out, open(f"{RES}/135_lineage_enrichment.json", "w"),
              indent=2, ensure_ascii=False)
    log("\nDONE → results/135_lineage_enrichment.json")


if __name__ == "__main__":
    main()
