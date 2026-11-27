#!/usr/bin/env python
"""步驟 221：把陽性對照改放在「腫瘤項」上，並檢查 D 的自我抵銷機制。

為什麼需要這一支
----------------
步驟 196 的陽性對照（CCND1 與 SOX11 要進前 50）是放在 `tumour_specific` 上判定的，
而 `tumour_specific = mean_tumour − mean_nontumour`，被減數只有 819 個非腫瘤 B 細胞、
其中 790 個（96.5%）來自單一捐贈者 M413。外部審查指出一個會把結論弄反的後果：
**通過陽性對照不等於那個方法看見了生物訊號**，有可能只是它的被減數剛好把 CCND1 排很高。

步驟 220 已經用我自己的資料證實了配對層級的反轉（平均看不見）：
  B–D 0.823 → 0.340、A–D 0.284 → 0.062（崩塌）
  A–C 0.052 → 0.409（相減之後才「一致」）
所以被減數不是加雜訊，而是**改變了哪些方法彼此一致**。

這一支要回答三個可直接否證的問題（全部不需要 GPU，只用 130 已經寫出的檔）
--------------------------------------------------------------------------
T1 **排名反轉**：每個方法、每個種子，CCND1 與 SOX11 在 `mean_tumour` 上的排名
   對照在 `tumour_specific` 上的排名。
   事先預測（來自審查，待本支驗證）：C 在 `mean_tumour` 上把 CCND1 排在倒數
   （> #3900 / 4012），但在 `tumour_specific` 上排進前 5。
   成立 → C 的陽性對照是被減數通過的，不是腫瘤訊號通過的。
T2 **以雜訊底為尺度的陽性對照**（取代排名門檻，因為排名門檻的判定會隨門檻翻邊：
   審查實測 top-20 時連 A 都失敗、top-201 時 D 只失敗 1/10）。
   規則：陽性對照必須在**每一個種子**都超過該方法自己雜訊底的 3 倍。
   事先預測：D 在 CCND1 上失敗（有種子低到約 1.4×），A、B、C 通過。
T3 **D 的自我抵銷機制**。D 用一個全域的 20 維 `delta = mean(z_腫瘤) − mean(z_的正常鄰居)`，
   而 CCND1 是 t(11;14) 的產物、是資料裡最大且最普遍的腫瘤位移，
   因此它**主導了 delta 的方向**。減掉 delta 等於沿著 CCND1 定義的軸移動，
   所以 `decode(z) − decode(z−delta)` 在 CCND1 上必然偏小 ——
   這個構造會系統性地漏掉「定義它自己方向」的那些基因。
   不需要跑 decoder 的檢驗：把 D 的每基因分數對 A 的每基因分數回歸，
   **限制在 A 的前 100 名**。若有飽和／自我抵銷，斜率應為**負**；
   作為對照，在全基因組上斜率應為正。

輸出 results/221_*
"""
import glob
import json
import re
import sys

import numpy as np
import pandas as pd
from scipy.stats import linregress

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
RES = f"{BASE}/results"
DET = "A_gene_space"
METHODS = [DET, "B_scvi_counterfactual", "C_contrastivevi", "D_latent_arithmetic"]
POS = ["CCND1", "SOX11"]
FOCUS = "SLC7A11"
FLOOR_MULT = 3.0            # T2 的事先門檻
TOPN_SAT = 100              # T3 限制在 A 的前 100 名


def log(*a):
    print("[221]", *a, flush=True)


def load():
    out = {}
    for p in sorted(glob.glob(f"{RES}/130_residual_by_method_s*.csv.gz")):
        m = re.search(r"_s(\d+)\.csv\.gz$", p)
        if not m:
            continue
        d = pd.read_csv(p)
        if set(d.method.unique()) != set(METHODS):
            continue
        out[int(m.group(1))] = {k: v.set_index("gene") for k, v in d.groupby("method")}
    return out


def rank_desc(s):
    """1 = 最大，與 130 的 (x > x[g]).sum() + 1 等價（無 ties 時）。"""
    return s.rank(ascending=False, method="min").astype(int)


def main():
    R = load()
    if not R:
        raise SystemExit("找不到 130_residual_by_method_s*.csv.gz")
    seeds = sorted(R)
    log("種子:", seeds)
    nG = len(R[seeds[0]][DET])
    out = {"seeds": seeds, "n_genes": int(nG)}

    # ------------------------------------------------------------- T1 ---
    rows = []
    for g in POS + [FOCUS]:
        for m in METHODS:
            rt, rs = [], []
            for s in seeds:
                d = R[s][m]
                rt.append(int(rank_desc(d.mean_tumour)[g]))
                rs.append(int(rank_desc(d.tumour_specific)[g]))
            rows.append({"gene": g, "method": m,
                         "rank_mean_tumour": rt, "rank_tumour_specific": rs,
                         "median_rank_tumour_term": int(np.median(rt)),
                         "median_rank_contrast": int(np.median(rs)),
                         "inversion": int(np.median(rt) - np.median(rs))})
    t1 = pd.DataFrame(rows)
    c_ccnd1 = t1[(t1.gene == "CCND1") & (t1.method == "C_contrastivevi")].iloc[0]
    t1_pred = bool(c_ccnd1["median_rank_tumour_term"] > 3900
                   and c_ccnd1["median_rank_contrast"] <= 5)
    out["T1"] = {"table": t1.to_dict(orient="records"),
                 "prediction": "C 把 CCND1 在腫瘤項排 >#3900、在相減後排 ≤#5",
                 "C_CCND1_tumour_term_median_rank": int(c_ccnd1["median_rank_tumour_term"]),
                 "C_CCND1_contrast_median_rank": int(c_ccnd1["median_rank_contrast"]),
                 "confirmed": t1_pred}

    log("\n=== T1：陽性對照在「腫瘤項」vs「相減後」的排名（共 %d 個基因）===" % nG)
    log("%-9s %-24s %22s %22s %10s" % ("gene", "method", "腫瘤項（各種子）",
                                       "相減後（各種子）", "反轉"))
    for r in rows:
        log("%-9s %-24s %22s %22s %+10d"
            % (r["gene"], r["method"],
               str(r["rank_mean_tumour"]), str(r["rank_tumour_specific"]), r["inversion"]))
    log("T1 預測成立 = %s（C 的 CCND1：腫瘤項 #%d → 相減後 #%d）"
        % (t1_pred, c_ccnd1["median_rank_tumour_term"], c_ccnd1["median_rank_contrast"]))

    # ------------------------------------------------------------- T2 ---
    t2, fails = {}, {}
    for m in METHODS:
        t2[m] = {}
        for g in POS + [FOCUS]:
            v = []
            for s in seeds:
                d = R[s][m]
                fl = float(d["mean_nontumour"].abs().mean())
                v.append(float(d.loc[g, "tumour_specific"] / fl))
            t2[m][g] = {"by_seed": v, "mean": float(np.mean(v)),
                        "min": float(np.min(v)),
                        "cv": float(np.std(v, ddof=1) / abs(np.mean(v))) if len(v) > 1 else 0.0,
                        "passes_3x_every_seed": bool(np.min(v) >= FLOOR_MULT)}
        fails[m] = [g for g in POS if not t2[m][g]["passes_3x_every_seed"]]
    out["T2"] = {"threshold": FLOOR_MULT, "values": t2, "failures": fails,
                 "prediction": "D 在 CCND1 上失敗；A、B、C 通過",
                 "confirmed": bool(fails["D_latent_arithmetic"] == ["CCND1"]
                                   and all(not fails[m] for m in METHODS
                                           if m != "D_latent_arithmetic"))}
    # 門檻敏感度：2× 與 3× 是否給出同一個答案
    out["T2"]["sensitivity"] = {
        f"{k}x": {m: [g for g in POS if t2[m][g]["min"] < k] for m in METHODS}
        for k in (2.0, 3.0, 4.0)}

    log("\n=== T2：陽性對照 ÷ 各方法自己的雜訊底（事先門檻：每個種子都要 ≥ %.1f×）===" % FLOOR_MULT)
    log("%-24s %18s %18s %18s" % ("method", "CCND1", "SOX11", FOCUS))
    for m in METHODS:
        log("%-24s" % m + "".join(
            " %8.2f (min %5.2f)" % (t2[m][g]["mean"], t2[m][g]["min"])
            for g in POS + [FOCUS]))
    for m in METHODS:
        log("  %-24s 未通過: %s" % (m, fails[m] or "—"))
    log("T2 預測成立 = %s" % out["T2"]["confirmed"])
    log("  門檻敏感度 2×/3×/4×: %s" % json.dumps(out["T2"]["sensitivity"], ensure_ascii=False))

    # ------------------------------------------------------------- T3 ---
    sat = {}
    for s in seeds:
        a = R[s][DET]["mean_tumour"]
        topA = a.sort_values(ascending=False).head(TOPN_SAT).index
        for m in METHODS[1:]:
            b = R[s][m]["mean_tumour"]
            g_all = linregress(a.to_numpy(), b.to_numpy())
            g_top = linregress(a.loc[topA].to_numpy(), b.loc[topA].to_numpy())
            sat.setdefault(m, {"slope_all": [], "slope_topA": []})
            sat[m]["slope_all"].append(float(g_all.slope))
            sat[m]["slope_topA"].append(float(g_top.slope))
    for m in sat:
        sat[m] = {k: {"by_seed": v, "mean": float(np.mean(v))} for k, v in sat[m].items()}
        sat[m]["negative_slope_in_topA"] = bool(sat[m]["slope_topA"]["mean"] < 0)
    # A 的前 50 名在各方法的中位數排名（反轉的另一個看法）
    med_top50 = {}
    for s in seeds:
        a = R[s][DET]["mean_tumour"]
        t50 = a.sort_values(ascending=False).head(50).index
        for view in ("mean_tumour", "tumour_specific"):
            for m in METHODS:
                r = rank_desc(R[s][m][view])
                med_top50.setdefault(view, {}).setdefault(m, []).append(
                    float(np.median(r.loc[t50])))
    med_top50 = {v: {m: float(np.median(x)) for m, x in d.items()}
                 for v, d in med_top50.items()}
    out["T3"] = {"regression_on_A": sat,
                 "median_rank_of_A_top50": med_top50,
                 "prediction": "D 在 A 的前 100 名內斜率為負，全基因組為正",
                 "confirmed": bool(sat["D_latent_arithmetic"]["negative_slope_in_topA"]
                                   and sat["D_latent_arithmetic"]["slope_all"]["mean"] > 0)}

    log("\n=== T3：各方法的腫瘤項對 A 的腫瘤項回歸（檢驗自我抵銷）===")
    log("%-24s %14s %22s" % ("method", "全基因組斜率", "限 A 前 %d 名的斜率" % TOPN_SAT))
    for m in METHODS[1:]:
        log("%-24s %14.4f %22.4f%s"
            % (m, sat[m]["slope_all"]["mean"], sat[m]["slope_topA"]["mean"],
               "   ← 負，符合自我抵銷" if sat[m]["negative_slope_in_topA"] else ""))
    log("\nA 的前 50 名在各方法的中位數排名（共 %d 個基因）：" % nG)
    for view, d in med_top50.items():
        log("  %-18s %s" % (view, {m[0]: int(x) for m, x in d.items()}))
    log("T3 預測成立 = %s" % out["T3"]["confirmed"])

    json.dump(out, open(f"{RES}/221_tumour_term_controls.json", "w"),
              indent=2, ensure_ascii=False)
    t1.to_csv(f"{RES}/221_positive_control_ranks.csv", index=False)
    log("\n→ results/221_tumour_term_controls.json、221_positive_control_ranks.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
