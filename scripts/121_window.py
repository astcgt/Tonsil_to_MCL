#!/usr/bin/env python
"""步驟 121：治療窗口 —— 「正常細胞碰到這個藥會怎樣」能用現有資料回答到什麼程度。

三個可以真的算的代理指標
------------------------
A. 泛譜系 CRISPR 依賴性（DepMap 24Q4，約 1,180 株）
   把 MAFG（蛋白質藥物的標的）與下游五個基因，放到所有譜系裡看。
   敲掉這個基因，MCL 株掉多少、其他譜系掉多少。
   ⚠️ DepMap 全是癌細胞株，不是正常細胞。它回答的是
      「哪些**增生中的**細胞離不開這個基因」，是治療窗口的下界而非正常組織毒性。
B. 旁系同源冗餘（MAFF / MAFK）
   如果正常細胞靠 MAFF/MAFK 補位，單打 MafG 的窗口就寬 —— 但藥效也可能被補掉。
   檢定：各株的 MAFG 依賴性 vs MAFF/MAFK 表現（相關係數）。
C. 真正的正常人類細胞：扁桃腺圖譜 46 個 B 細胞狀態的表現量
   這是唯一一組真正的正常原始細胞。回答「標的在正常對應細胞裡開著嗎」。

輸出 results/121_*
"""
import os
import sys
import json

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "8"

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
DM = f"{BASE}/data/depmap"
RES = f"{BASE}/results"

TARGET = "MAFG"
PARALOG = ["MAFF", "MAFK"]
PARTNER = ["NFE2L2", "BACH1"]
PATHWAY = ["SLC7A11", "SLC3A2", "GCLM", "GCLC", "AIFM2"]
CTRL_ESS = ["RPL13A", "POLR2A", "RAN"]          # 泛必需，當「沒有窗口」的參照
CTRL_MCL = ["CCND1", "BTK", "BCL2"]             # MCL 臨床相關，當「已知藥物標的」參照
GENES = [TARGET] + PARALOG + PARTNER + PATHWAY + CTRL_ESS + CTRL_MCL
DEP_CUT = -0.5                                   # DepMap 慣用的「依賴」門檻


def log(*a):
    print("[121]", *a, flush=True)


def pick_cols(path, genes):
    head = pd.read_csv(path, nrows=0)
    want = {}
    for c in head.columns[1:]:
        sym = c.split(" (")[0]
        if sym in genes:
            want[sym] = c
    return want


def load_matrix(path, genes, label):
    want = pick_cols(path, genes)
    missing = sorted(set(genes) - set(want))
    if missing:
        log(f"{label}: 找不到 {missing}")
    df = pd.read_csv(path, usecols=[head_index(path)] + list(want.values()), index_col=0)
    df.columns = [c.split(" (")[0] for c in df.columns]
    return df


def head_index(path):
    return pd.read_csv(path, nrows=0).columns[0]


def main():
    out = {}

    # ---------- A. 泛譜系依賴性 ----------
    eff = load_matrix(f"{DM}/CRISPRGeneEffect.csv", GENES, "GeneEffect")
    log("GeneEffect:", eff.shape)
    mod = pd.read_csv(f"{DM}/Model.csv", index_col=0)
    lin = mod["OncotreeLineage"].reindex(eff.index)
    sub = mod["OncotreeSubtype"].reindex(eff.index)
    mcl_ids = [i for i in pd.read_csv(f"{DM}/MCL_cell_lines_DepMap_24Q4.csv")["ModelID"]
               if i in eff.index]
    log("有 CRISPR 資料的 MCL 株:", len(mcl_ids),
        list(mod.loc[mcl_ids, "StrippedCellLineName"]))

    rows = []
    for g in eff.columns:
        v = eff[g].dropna()
        m = v.reindex(mcl_ids).dropna()
        rest = v.drop(index=m.index, errors="ignore")
        p = (mannwhitneyu(m, rest, alternative="less").pvalue
             if len(m) >= 3 else np.nan)
        rows.append({
            "gene": g,
            "n_lines": len(v),
            "median_all": v.median(),
            "mean_all": v.mean(),
            "frac_dependent_all": float((v < DEP_CUT).mean()),
            "n_mcl": len(m),
            "median_mcl": m.median() if len(m) else np.nan,
            "frac_dependent_mcl": float((m < DEP_CUT).mean()) if len(m) else np.nan,
            "delta_mcl_minus_rest": (m.median() - rest.median()) if len(m) else np.nan,
            "mwu_p_mcl_more_dependent": p,
            "min_lineage_median": np.nan, "min_lineage": "",
        })
    dep = pd.DataFrame(rows).set_index("gene")

    # 每譜系中位數 → 哪一類細胞最受不了
    linmed = {}
    for g in eff.columns:
        s = eff[g].groupby(lin).median().dropna()
        n = eff[g].groupby(lin).count()
        s = s[n >= 10]
        linmed[g] = s
        if len(s):
            dep.loc[g, "min_lineage_median"] = float(s.min())
            dep.loc[g, "min_lineage"] = str(s.idxmin())
    L = pd.DataFrame(linmed)
    L.index.name = "OncotreeLineage"
    L.round(4).to_csv(f"{RES}/121_depmap_lineage_median.csv")
    dep.round(4).to_csv(f"{RES}/121_depmap_summary.csv")
    log("\n=== 泛譜系依賴性（越負越必需；< -0.5 算依賴） ===\n",
        dep[["n_lines", "median_all", "frac_dependent_all", "n_mcl", "median_mcl",
             "frac_dependent_mcl", "min_lineage", "min_lineage_median"]].round(3).to_string())

    out["depmap"] = json.loads(dep.round(5).to_json(orient="index"))
    out["lineage_median"] = json.loads(L.round(4).to_json(orient="index"))
    out["mcl_lines"] = list(mod.loc[mcl_ids, "StrippedCellLineName"])

    # 其他淋巴系（非 MCL）—— 最接近「正常 B 細胞」的可得代理
    lymph = [i for i in eff.index if lin.get(i) == "Lymphoid" and i not in mcl_ids]
    bcell = [i for i in lymph if isinstance(sub.get(i), str)
             and ("B-Cell" in sub.get(i) or "Burkitt" in sub.get(i)
                  or "Lymphoma" in sub.get(i))]
    out["lymphoid_nonmcl"] = {
        g: {"n": int(eff.loc[lymph, g].notna().sum()),
            "median": float(eff.loc[lymph, g].median()),
            "frac_dependent": float((eff.loc[lymph, g] < DEP_CUT).mean())}
        for g in eff.columns}
    out["bcell_nonmcl"] = {
        "n_lines": len(bcell),
        **{g: {"median": float(eff.loc[bcell, g].median()),
               "frac_dependent": float((eff.loc[bcell, g] < DEP_CUT).mean())}
           for g in eff.columns}}
    log("\n非 MCL 的 B 細胞腫瘤株 (n=%d) 中位依賴性:" % len(bcell))
    log("\n".join("  %-8s %+.3f" % (g, out["bcell_nonmcl"][g]["median"]) for g in eff.columns))

    # ---------- B. 旁系同源冗餘 ----------
    expr = load_matrix(f"{DM}/OmicsExpressionProteinCodingGenesTPMLogp1.csv",
                       [TARGET] + PARALOG + PATHWAY, "Expression")
    log("Expression:", expr.shape)
    red = {}
    common = eff.index.intersection(expr.index)
    for g in PARALOG + [TARGET]:
        if g not in expr.columns:
            continue
        a = eff.loc[common, TARGET]
        b = expr.loc[common, g]
        ok = a.notna() & b.notna()
        r, p = spearmanr(a[ok], b[ok])
        red[g] = {"n": int(ok.sum()), "spearman_MAFG_dependency_vs_expr": float(r),
                  "p": float(p),
                  "expr_median_log2tpm1": float(b[ok].median())}
    out["paralog_redundancy"] = red
    log("\n=== MAFG 依賴性 vs 旁系同源表現（正相關 = 旁系多的株反而更不依賴 MAFG） ===")
    for g, d in red.items():
        log("  %-6s rho = %+.3f  p = %.3g  (中位表現 %.2f)"
            % (g, d["spearman_MAFG_dependency_vs_expr"], d["p"], d["expr_median_log2tpm1"]))

    # ---------- C. 真正的正常細胞：46 個扁桃腺 B 狀態 ----------
    pb = pd.read_csv(f"{RES}/15_normal_state_pseudobulk.csv", index_col=0)
    tum = pd.read_csv(f"{RES}/15_gene_residual_tumor.csv").set_index("gene")
    norm = {}
    for g in [TARGET] + PARALOG + PARTNER + PATHWAY:
        if g not in pb.columns:
            norm[g] = None
            continue
        col = pb[g]
        rec = {
            "normal_states_median": float(col.median()),
            "normal_states_max": float(col.max()),
            "normal_states_argmax": str(col.idxmax()),
            "frac_states_expressed_gt0.1": float((col > 0.1).mean()),
            "ncsMBC": float(col["ncsMBC"]) if "ncsMBC" in col.index else None,
        }
        if g in tum.index:
            rec.update({
                "tumour_observed": float(tum.loc[g, "mean_observed"]),
                "tumour_expected": float(tum.loc[g, "mean_expected"]),
                "tumour_residual": float(tum.loc[g, "mean_residual"]),
            })
        norm[g] = rec
    out["normal_tonsil_B"] = norm
    log("\n=== 正常扁桃腺 B 細胞（46 狀態 lognorm 表現） vs 腫瘤實測 ===")
    log("%-9s %9s %9s %9s %9s %9s" % ("gene", "正常中位", "正常最大", "ncsMBC", "腫瘤實測", "腫瘤期望"))
    for g, d in norm.items():
        if not d:
            continue
        log("%-9s %9.3f %9.3f %9s %9s %9s" % (
            g, d["normal_states_median"], d["normal_states_max"],
            "%.3f" % d["ncsMBC"] if d["ncsMBC"] is not None else "-",
            "%.3f" % d["tumour_observed"] if "tumour_observed" in d else "-",
            "%.3f" % d["tumour_expected"] if "tumour_expected" in d else "-"))

    out["notes"] = {
        "dep_cut": DEP_CUT,
        "caveat_depmap": "DepMap 皆為癌細胞株，非正常細胞；此處只能當『增生細胞的依賴性』代理",
        "caveat_normal": "唯一真正的正常人類細胞是扁桃腺 B 細胞；本專案未取得任何正常非 B 組織資料",
    }
    json.dump(out, open(f"{RES}/121_window.json", "w"), indent=2, ensure_ascii=False)
    log("\nDONE → results/121_window.json, 121_depmap_summary.csv, 121_depmap_lineage_median.csv")


if __name__ == "__main__":
    main()
