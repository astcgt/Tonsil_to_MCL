#!/usr/bin/env python
"""步驟 223：那個被減數到底是「非腫瘤 B 細胞」還是「M413 這一位捐贈者」？

背景
----
`tumour_specific = mean_tumour − mean_nontumour`，而被減數來自 **819 個非腫瘤 B 細胞，
其中 790 個（96.5%）來自捐贈者 M413**，另一位 M102 只有 29 個。
步驟 220／221 已經證明（11 個種子、含 6 個保留集種子）這個被減數會**改變哪些方法彼此一致**，
而且 `ContrastiveVI` 的分數與它相關 ρ = 0.967、與自己的腫瘤項只有 −0.012。

但那些都還沒回答一個更直接的問題：**換掉捐贈者，基因排名會不會變？**
`130 --group-by-donor` 把非腫瘤組按捐贈者拆開輸出，這一支用它來回答。

事先寫死的預測（本檔在任何 `_ds*` 結果落地之前寫成）
----------------------------------------------------
D1 **正確性自檢（守門，失敗就作廢整支）**：`_ds{N}` 的合併欄位
   （`mean_tumour`、`mean_nontumour`、`tumour_specific`）必須與 `_s{N}` **逐位元相同**。
   `--group-by-donor` 只增加分組、不改計算，所以任何差異都代表補丁改了行為 → 作廢。
D2 **被減數其實是一位捐贈者**：把被減數換成「只用 M413」與「只用 M102」各重算一次排名。
   預測：合併 vs 只用 M413 的 Spearman **> 0.99**（因為 790/819 = 96.5%），
   而合併 vs 只用 M102 **明顯更低（< 0.90）**。
   成立 → 「非腫瘤比較組」在數值上等於 M413 一個人。
D3 **換捐贈者對各方法的衝擊，與被減數佔比成正比**：
   以各方法自己的雜訊底為單位，報告 `SLC7A11` 在三種被減數下的值。
   預測：變動幅度 **C > D > B > A**（對應 220 量到的被減數變異佔比
   C 1.04 > D 1.71?…… 注意 D 的佔比最高但其分數與腫瘤項仍正相關，
   所以這裡只事先排序 **C 最大、A 最小**，中間兩者不預測）。
D4 **2 位捐贈者估不出這個被減數**：每個基因在兩位捐贈者之間的差距
   （|m_M413 − m_M102|）的中位數，對照該方法 `mean_nontumour` 的跨基因標準差。
   預測：**比值 > 1**，即「換一位捐贈者造成的改變」大於「基因之間的真實差異」
   → 這一項在 n=2 位捐贈者下不可估計，與跑幾個種子完全無關。

輸出 results/223_*
"""
import glob
import json
import re
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
RES = f"{BASE}/results"
METHODS = ["A_gene_space", "B_scvi_counterfactual", "C_contrastivevi", "D_latent_arithmetic"]
FOCUS = "SLC7A11"
POS = ["CCND1", "SOX11"]
POOLED_COLS = ["mean_tumour", "mean_nontumour", "tumour_specific"]


def log(*a):
    print("[223]", *a, flush=True)


def main():
    ds = {}
    for p in sorted(glob.glob(f"{RES}/130_residual_by_method_ds*.csv.gz")):
        m = re.search(r"_ds(\d+)\.csv\.gz$", p)
        if m:
            ds[int(m.group(1))] = p
    if not ds:
        raise SystemExit("還沒有 130_residual_by_method_ds*.csv.gz")
    seeds = sorted(ds)
    log("捐贈者拆分的種子:", seeds)

    out = {"seeds": seeds, "gates": {}}

    # ----------------------------------------------------------- D1 守門 ---
    d1 = {}
    for s in seeds:
        a = pd.read_csv(f"{RES}/130_residual_by_method_s{s}.csv.gz")
        b = pd.read_csv(ds[s])
        A = {k: v.set_index("gene") for k, v in a.groupby("method")}
        B = {k: v.set_index("gene") for k, v in b.groupby("method")}
        per = {}
        for mth in METHODS:
            same = all(np.array_equal(A[mth][c].to_numpy(), B[mth].loc[A[mth].index, c].to_numpy())
                       for c in POOLED_COLS)
            mad = max(float(np.max(np.abs(A[mth][c].to_numpy()
                                         - B[mth].loc[A[mth].index, c].to_numpy())))
                      for c in POOLED_COLS)
            per[mth] = {"bit_identical": bool(same), "max_abs_diff": mad}
        d1[s] = per
    ok = all(v["bit_identical"] for s in d1 for v in d1[s].values())
    out["gates"]["D1_patch_correctness"] = {
        "requirement": "_ds 的合併欄位必須與 _s 逐位元相同",
        "per_seed": d1, "pass": bool(ok)}
    log("\n=== D1：補丁正確性（合併欄位是否逐位元重現）===")
    for s in seeds:
        log("  種子 %d: %s | 最大絕對差 %.3e" % (
            s, all(v["bit_identical"] for v in d1[s].values()),
            max(v["max_abs_diff"] for v in d1[s].values())))
    if not ok:
        out["verdict"] = "VOID_patch_changed_behaviour"
        json.dump(out, open(f"{RES}/223_donor_leverage.json", "w"), indent=2, ensure_ascii=False)
        log("D1 失敗 → `--group-by-donor` 改變了既有計算，整支作廢，不下任何結論。")
        return 1

    # 取出 per-donor 欄位
    b0 = pd.read_csv(ds[seeds[0]])
    dcols = [c for c in b0.columns if c.startswith("mean_nontumour_")]
    donors = [c.replace("mean_nontumour_", "") for c in dcols]
    ncell = {d: int(b0[f"n_nontumour_{d}"].iloc[0]) for d in donors}
    log("\n捐贈者欄位: %s | 細胞數 %s" % (donors, ncell))
    out["donors"] = ncell
    big = max(ncell, key=ncell.get)
    small = min(ncell, key=ncell.get)

    # --------------------------------------------- D2 / D3 / D4 ---
    d2, d3, d4 = {}, {}, {}
    for mth in METHODS:
        rho_big, rho_small, rho_cross = [], [], []
        foc, d4v = {k: [] for k in ["pooled", big, small]}, []
        for s in seeds:
            b = pd.read_csv(ds[s])
            g = b[b.method == mth].set_index("gene")
            mt = g["mean_tumour"]
            pooled = mt - g["mean_nontumour"]
            v_big = mt - g[f"mean_nontumour_{big}"]
            v_small = mt - g[f"mean_nontumour_{small}"]
            rho_big.append(spearmanr(pooled, v_big).statistic)
            rho_small.append(spearmanr(pooled, v_small).statistic)
            rho_cross.append(spearmanr(v_big, v_small).statistic)
            fl = float(g["mean_nontumour"].abs().mean())
            foc["pooled"].append(float(pooled[FOCUS] / fl))
            foc[big].append(float(v_big[FOCUS] / fl))
            foc[small].append(float(v_small[FOCUS] / fl))
            gap = (g[f"mean_nontumour_{big}"] - g[f"mean_nontumour_{small}"]).abs()
            d4v.append(float(np.median(gap) / np.std(g["mean_nontumour"], ddof=1)))
        d2[mth] = {f"rho_pooled_vs_{big}_only": float(np.mean(rho_big)),
                   f"rho_pooled_vs_{small}_only": float(np.mean(rho_small)),
                   "rho_between_the_two_donors": float(np.mean(rho_cross))}
        d3[mth] = {k: {"mean": float(np.mean(v)), "by_seed": v} for k, v in foc.items()}
        d3[mth]["swing_big_minus_small"] = float(np.mean(foc[big]) - np.mean(foc[small]))
        d4[mth] = float(np.mean(d4v))

    p_d2 = bool(all(d2[m][f"rho_pooled_vs_{big}_only"] > 0.99 for m in METHODS)
                and all(d2[m][f"rho_pooled_vs_{small}_only"] < 0.90 for m in METHODS))
    sw = {m: abs(d3[m]["swing_big_minus_small"]) for m in METHODS}
    p_d3 = bool(max(sw, key=sw.get) == "C_contrastivevi"
                and min(sw, key=sw.get) == "A_gene_space")
    p_d4 = bool(all(v > 1.0 for v in d4.values()))

    out["gates"]["D2_subtrahend_is_one_donor"] = {
        "prediction": f"合併 vs 只用 {big} 的 ρ > 0.99，合併 vs 只用 {small} 的 ρ < 0.90",
        "values": d2, "pass": p_d2}
    out["gates"]["D3_swing_by_method"] = {
        "prediction": "換捐贈者的衝擊 C 最大、A 最小",
        "focus_gene": FOCUS, "values": d3,
        "abs_swing": sw, "pass": p_d3}
    out["gates"]["D4_not_estimable_from_two_donors"] = {
        "prediction": "兩位捐贈者之間的逐基因差距中位數 / 跨基因標準差 > 1",
        "values": d4, "pass": p_d4}

    log("\n=== D2：換掉被減數的來源，基因排名還一樣嗎（Spearman）===")
    log("%-24s %18s %18s %18s" % ("method", f"合併 vs 只用{big}", f"合併 vs 只用{small}",
                                  "兩位捐贈者之間"))
    for m in METHODS:
        log("%-24s %18.4f %18.4f %18.4f"
            % (m, d2[m][f"rho_pooled_vs_{big}_only"], d2[m][f"rho_pooled_vs_{small}_only"],
               d2[m]["rho_between_the_two_donors"]))
    log("D2 預測成立 = %s" % p_d2)

    log("\n=== D3：%s 在三種被減數下的值（單位＝各方法自己的雜訊底）===" % FOCUS)
    log("%-24s %12s %12s %12s %12s" % ("method", "合併", f"只用{big}", f"只用{small}", "差距"))
    for m in METHODS:
        log("%-24s %12.2f %12.2f %12.2f %12.2f"
            % (m, d3[m]["pooled"]["mean"], d3[m][big]["mean"], d3[m][small]["mean"],
               d3[m]["swing_big_minus_small"]))
    log("D3 預測成立 = %s（最大 %s、最小 %s）"
        % (p_d3, max(sw, key=sw.get), min(sw, key=sw.get)))

    log("\n=== D4：兩位捐贈者的差距 vs 基因之間的真實差異 ===")
    for m in METHODS:
        log("  %-24s 比值 %.3f%s" % (m, d4[m], "   ← >1，估不出來" if d4[m] > 1 else ""))
    log("D4 預測成立 = %s" % p_d4)

    out["verdict"] = {"D2": p_d2, "D3": p_d3, "D4": p_d4}
    json.dump(out, open(f"{RES}/223_donor_leverage.json", "w"), indent=2, ensure_ascii=False)
    log("\n→ results/223_donor_leverage.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
