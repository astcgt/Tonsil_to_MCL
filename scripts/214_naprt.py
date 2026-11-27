#!/usr/bin/env python
"""步驟 214：`NAPRT` 是不是 `NAMPT` 的 B 細胞選擇性的機制與病人篩選標記？

這一題是被一個「失敗」指出來的
------------------------------
步驟 213 的 P2 不但沒過，方向還是反的：`NAMPT` 表現越高，依賴性越**弱**
（非淋巴系 1,012 株 r = +0.148、p = 2.4×10⁻⁶）。
那不是自我成癮，而是 NAMPT 抑制劑文獻裡的標準圖像 ——
細胞做 NAD⁺ 有兩條路：
  補救途徑（salvage）        菸鹼醯胺 → NMN → NAD⁺，限速酵素 **NAMPT**
  Preiss-Handler 途徑        菸鹼酸  → NaMN → NAD⁺，限速酵素 **NAPRT**
只有**失去 NAPRT** 的細胞才會真的被 NAMPT 抑制劑殺死。
所以如果 B 細胞腫瘤普遍 NAPRT 低，那它就同時是
（a）步驟 213 那個 92% vs 27% 的機制，（b）一個現成的病人篩選標記。

事前指定
--------
P1  `NAPRT` 表現在 B 細胞腫瘤株顯著低於非淋巴系株。
P2  `NAPRT` 表現與 `NAMPT` 依賴性正相關（NAPRT 高 → gene effect 較不負 → 較不依賴）。
P3  **關鍵的一關**：P2 必須在**非淋巴系**株內也成立。
    否則 `NAPRT` 只是「是不是淋巴系」的同義詞，解釋不了任何事。
P4（陰性對照）`NAPRT` 表現不該同樣地預測一個無關必需基因（`RPL13A`）的依賴性。
P5  病人層面的窗口：在本專案的分化配對殘差裡，
    `NAPRT` 在 MCL 腫瘤與分化配對的正常 B 細胞各是多少？
    如果兩者都低，那窗口不來自 NAPRT 本身（要照實說）。

輸出 results/214_*
"""
import os
import json

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "8"

import numpy as np
import pandas as pd
from scipy import stats

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
RES = f"{BASE}/results"
DM = f"{BASE}/data/depmap"
GENES_EXP = ["NAPRT", "NAMPT", "NMRK1", "NMRK2", "NADSYN1", "QPRT"]
GENES_EFF = ["NAMPT", "NAPRT", "RPL13A"]


def log(*a):
    print("[214]", *a, flush=True)


def load(path, genes):
    hdr = pd.read_csv(path, nrows=0)
    f = hdr.columns[0]
    m = {}
    for c in hdr.columns[1:]:
        m.setdefault(c.split(" (")[0], c)
    have = [g for g in genes if g in m]
    D = pd.read_csv(path, usecols=[f] + [m[g] for g in have]).set_index(f)
    D.columns = [c.split(" (")[0] for c in D.columns]
    return D, have


def main():
    X, hx = load(f"{DM}/OmicsExpressionProteinCodingGenesTPMLogp1.csv", GENES_EXP)
    E, he = load(f"{DM}/CRISPRGeneEffect.csv", GENES_EFF)
    log("表現 %s｜依賴 %s" % (hx, he))
    mod = pd.read_csv(f"{DM}/Model.csv").set_index("ModelID")
    lin = mod["OncotreeLineage"].astype(str)
    dis = mod["OncotreePrimaryDisease"].astype(str)

    idx = X.index.intersection(E.index)
    d = pd.DataFrame(index=idx)
    for g in hx:
        d["exp_" + g] = X[g].reindex(idx)
    for g in he:
        d["eff_" + g] = E[g].reindex(idx)
    d["lineage"] = lin.reindex(idx)
    d["disease"] = dis.reindex(idx)
    d["lymph"] = d.lineage.astype(str).str.contains("Lymph", case=False)
    d["Bcell"] = d.lymph & ~d.disease.astype(str).str.contains(
        "T-Cell|T Cell|NK", case=False, regex=True)
    d = d.dropna(subset=["exp_NAPRT", "eff_NAMPT"])
    d.to_csv(f"{RES}/214_naprt_lines.csv")
    log("可用株數 %d（B 細胞腫瘤 %d、非淋巴系 %d）"
        % (len(d), int(d.Bcell.sum()), int((~d.lymph).sum())))

    out = {"n_lines": int(len(d)), "n_Bcell": int(d.Bcell.sum()),
           "n_non_lymphoid": int((~d.lymph).sum())}

    B = d[d.Bcell]; O = d[~d.lymph]
    u = stats.mannwhitneyu(B.exp_NAPRT, O.exp_NAPRT, alternative="less")
    out["P1_NAPRT_low_in_Bcell"] = {
        "median_Bcell": float(B.exp_NAPRT.median()),
        "median_non_lymphoid": float(O.exp_NAPRT.median()),
        "mannwhitney_p": float(u.pvalue),
        "frac_Bcell_below_1": float((B.exp_NAPRT < 1).mean()),
        "frac_other_below_1": float((O.exp_NAPRT < 1).mean()),
        "pass": bool(u.pvalue < 0.05
                     and B.exp_NAPRT.median() < O.exp_NAPRT.median())}
    log("\nP1 `NAPRT` 表現：B 細胞腫瘤中位 %.3f vs 非淋巴系 %.3f，單尾 p = %.2e"
        % (B.exp_NAPRT.median(), O.exp_NAPRT.median(), u.pvalue))
    log("   低於 1 (log2 TPM+1) 的比例：B 細胞 %.0f%% vs 其他 %.0f%% → %s"
        % (100 * (B.exp_NAPRT < 1).mean(), 100 * (O.exp_NAPRT < 1).mean(),
           "過" if out["P1_NAPRT_low_in_Bcell"]["pass"] else "不過"))

    for nm, sub, key in [("P2_all_lines", d, "P2"), ("P3_non_lymphoid", O, "P3")]:
        r, p = stats.spearmanr(sub.exp_NAPRT, sub.eff_NAMPT)
        out[nm] = {"n": int(len(sub)), "spearman_r": float(r), "p": float(p),
                   "expected_sign": "正（NAPRT 高 → 較不依賴 NAMPT）",
                   "pass": bool(r > 0 and p < 0.05)}
        log("%s `NAPRT` 表現 vs `NAMPT` 依賴（n=%d）：r = %+.3f，p = %.2e → %s"
            % (key, len(sub), r, p, "過" if out[nm]["pass"] else "不過"))

    r4, p4 = stats.spearmanr(d.exp_NAPRT, d.eff_RPL13A, nan_policy="omit")
    r2 = out["P2_all_lines"]["spearman_r"]
    out["P4_negative_control"] = {"gene": "RPL13A", "spearman_r": float(r4),
                                 "p": float(p4),
                                 "pass": bool(abs(r4) < abs(r2))}
    log("P4 陰性對照 `RPL13A`：r = %+.3f，p = %.2e（標的 %+.3f）→ %s"
        % (r4, p4, r2, "過" if out["P4_negative_control"]["pass"] else "不過"))

    # 分層：NAPRT 低 vs 高的株，NAMPT 依賴性差多少
    thr = 1.0
    lo = d[d.exp_NAPRT < thr]; hi = d[d.exp_NAPRT >= thr]
    uu = stats.mannwhitneyu(lo.eff_NAMPT, hi.eff_NAMPT, alternative="less")
    out["stratified"] = {
        "threshold_log2tpm1": thr,
        "n_naprt_low": int(len(lo)), "n_naprt_high": int(len(hi)),
        "median_NAMPT_effect_low": float(lo.eff_NAMPT.median()),
        "median_NAMPT_effect_high": float(hi.eff_NAMPT.median()),
        "frac_dependent_low": float((lo.eff_NAMPT < -0.5).mean()),
        "frac_dependent_high": float((hi.eff_NAMPT < -0.5).mean()),
        "mannwhitney_p": float(uu.pvalue)}
    log("\n以 `NAPRT` < %.0f 分層：低 %d 株 NAMPT effect 中位 %.3f（%.0f%% 依賴）"
        " vs 高 %d 株 %.3f（%.0f%% 依賴），p = %.2e"
        % (thr, len(lo), lo.eff_NAMPT.median(),
           100 * (lo.eff_NAMPT < -0.5).mean(), len(hi), hi.eff_NAMPT.median(),
           100 * (hi.eff_NAMPT < -0.5).mean(), uu.pvalue))

    # P5：病人層面
    R = pd.read_csv(f"{RES}/15_gene_residual_tumor.csv", index_col=0)
    pat = {}
    for g in ["NAPRT", "NAMPT", "NMRK1", "QPRT", "SLC7A11"]:
        if g in R.index:
            pat[g] = {"mcl_observed": float(R.loc[g, "mean_observed"]),
                      "normal_b_expected": float(R.loc[g, "mean_expected"]),
                      "residual": float(R.loc[g, "mean_residual"])}
    out["P5_patient_level"] = pat
    log("\nP5 病人層面（log1p CP10K）：")
    for g, v in pat.items():
        log("   %-8s MCL %.4f｜分化配對正常 B %.4f｜殘差 %+.4f"
            % (g, v["mcl_observed"], v["normal_b_expected"], v["residual"]))

    passed = [k for k, v in out.items() if isinstance(v, dict) and v.get("pass") is True]
    failed = [k for k, v in out.items() if isinstance(v, dict) and v.get("pass") is False]
    out["summary"] = {"passed": passed, "failed": failed}
    log("\n通過：%s｜未通過：%s"
        % (", ".join(passed) or "（無）", ", ".join(failed) or "（無）"))
    json.dump(out, open(f"{RES}/214_naprt.json", "w"), indent=2, ensure_ascii=False)
    log("DONE → results/214_*")


if __name__ == "__main__":
    main()
