#!/usr/bin/env python
"""步驟 220：四個方法的「不同意」到底來自哪裡 —— 取代步驟 196 的設計。

為什麼要取代 196
----------------
196 問的是「方法之間的不同意是種子雜訊還是真的」，並用
R = 方法間排名離散 / 方法內排名離散 加上硬門檻（>3 真的、<1.5 雜訊）來判定。
外部審查把這個設計打掉了，三個理由都成立，而且都已用資料驗證：

1. **問錯了問題。** `tumour_specific = mean_tumour − mean_nontumour`，而被減數
   `mean_nontumour` 只有 **819 個非腫瘤 B 細胞，其中 790 個（96.5%）來自同一位
   捐贈者 M413**，另一位 M102 只有 29 個。實測（種子 0、4012 個基因）：
   ContrastiveVI 的基因排名與「被減數的反號」相關達 **ρ = 0.965**、
   前 201 名重疊 84%，而與它自己的腫瘤項只有 **ρ = 0.025**。
   也就是說 C 排的是一位捐贈者的非腫瘤殘差，不是腫瘤專一性。
   更關鍵的對照：B 和 D 在 `mean_tumour` 上相關 **0.791**，在相減後只剩 **0.334**。
   → 196 要裁決的那個「不同意」，有很大一部分是**被減數製造出來的**。
2. **門檻不帶資訊。** R > 3 在全基因組會對 **41%（1651/4012）** 的基因成立，
   SLC7A11 的 R = 5.10 只是第 84 百分位。而 **CCND1 自己也有顯著的方法效應**
   （F(2,9) = 76.9, p = 2.2e-6）。「方法之間有系統性差異」是這個方法組合的
   全域性質，不是針對 SLC7A11 的證據。
3. **分母結構性錯誤。** 方法 A 是確定性的（實測四個種子逐位元相同、
   方法內排名 SD 恆為 0），把它平均進「方法內雜訊」會讓分母剛好變成 3/4，
   R 一律被放大 **4/3（+33.3%）**；有 643 個基因（16.0%）光是因為這一項就翻邊。
   （SLC7A11 不在其中：5.10 → 3.82，兩邊都 > 3。真正讓 196 停住的是 P1。）

另外 196 的 P1 門檻**用已存在的資料就已經違反**：已發表的那一輪裡
D 把 SOX11 排在 #96，所以 196 永遠只會回傳 `VOID_positive_control_failed`。
那是 D 這個方法的真實弱點，不是流程壞了（P0 通過，且 `_s0b` 的獨立重跑證明
固定種子下**連 GPU 訓練都完全決定性**：四個方法逐位元相同、最大絕對差 0.000e+00）。

這一版的三層設計（事先指定，不用任何硬門檻比值）
-------------------------------------------------
**第 0 層 —— 先修輸入，或把它當成結論本身。**
  把 `tumour_specific` 拆成兩項，報告每個方法
  ρ(ts, mean_tumour)、ρ(ts, −mean_nontumour)、被減數在 var(ts) 的佔比；
  然後把方法比較**改用 `mean_tumour` 重算一次**，兩者並列。
**第 1 層 —— 方法層級，主要結論。不需要門檻，n=2 個種子就已決定性。**
  多特質多方法矩陣：對角線是同方法跨種子的再測信度，非對角是同種子跨方法的一致度，
  再用 ρ_obs / sqrt(ρ_m · ρ_m') 做衰減校正。
  校正後若仍遠低於 1，代表兩個方法量的**不是同一個東西**，與種子無關。
**第 2 層 —— 單一基因（SLC7A11）。用各方法自己的雜訊底標準化，不用排名。**
  y = tumour_specific / noise_floor。排名是有界、非常態、且離散程度隨位置
  呈倒 U 形（實測中位數 2.8 → 231 → 96），不可以當離散統計量的輸入。

事先寫死的判定（看保留集之前就固定）
------------------------------------
宣告集 = 種子 0–4（這些已經被十幾種統計量看過，任何重新指定的門檻都是事後的）
保留集 = 種子 5–10（從未被檢視；用來確認）

L0：若四法在 `mean_tumour` 上的平均非對角 Spearman，比在 `tumour_specific` 上
    **高出 > 0.15**，則判定「不同意主要由 819 細胞／單一捐贈者的被減數製造」。
L1：若方法內再測信度 > 0.70 而衰減校正後的方法間一致度 < 0.70，
    則判定「不同意是系統性的，不是種子雜訊」。
L2：SLC7A11 的可報告陳述是「它落在方法依賴性的全基因組前十分位」
    （不是「R > 3」）。同時必須報告它在每個方法各自雜訊底的倍數，
    且**任何低於該方法雜訊底的基因不得引用排名**。
確認規則：每一個在 0–4 宣告的不等式，都必須在 5–10 上同向成立；
    不成立就照實報告落差，不改門檻。

輸出 results/220_*
"""
import glob
import json
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy.stats import f as fdist
from scipy.stats import spearmanr, t as tdist

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
RES = f"{BASE}/results"
DET = "A_gene_space"                 # 確定性方法，不進「方法內雜訊」的池
STOCH = ["B_scvi_counterfactual", "C_contrastivevi", "D_latent_arithmetic"]
METHODS = [DET] + STOCH
SHARED_FIT = ("B_scvi_counterfactual", "D_latent_arithmetic")   # 同一個 scVI fit
EXPLORE = [0, 1, 2, 3, 4]
FOCUS = "SLC7A11"
POS = ["CCND1", "SOX11"]
PANEL = ["SLC7A11", "SLC3A2", "GCLM", "GCLC", "AIFM2", "MAFG", "MAFF", "MAFK",
         "BACH1", "NFE2L2", "CCND1", "SOX11", "TXN", "GPX4", "HMOX1", "NQO1"]
L0_GAIN = 0.15                       # 事先指定
L1_REL, L1_AGREE = 0.70, 0.70        # 事先指定


def log(*a):
    print("[220]", *a, flush=True)


def load():
    """{seed: {method: DataFrame indexed by gene}}，只收齊全四個方法的種子。"""
    out = {}
    for p in sorted(glob.glob(f"{RES}/130_residual_by_method_s*.csv.gz")):
        m = re.search(r"_s(\d+)\.csv\.gz$", p)
        if not m:                                  # _s0b 等非整數 tag 不納入種子集
            continue
        s = int(m.group(1))
        d = pd.read_csv(p)
        if set(d.method.unique()) != set(METHODS):
            log(f"種子 {s} 方法不齊，跳過")
            continue
        g = {k: v.set_index("gene") for k, v in d.groupby("method")}
        genes = g[DET].index
        assert all(g[k].index.equals(genes) for k in METHODS), f"種子 {s} 基因順序不一致"
        for k in METHODS:                          # 審查第 11 點：NaN 會讓排名定義分歧
            assert g[k][["mean_tumour", "mean_nontumour", "tumour_specific"]].notna().all().all(), \
                f"種子 {s} 方法 {k} 有 NaN"
        out[s] = g
    return out


def floor_of(df):
    """130 自己的雜訊底定義：非腫瘤殘差的平均絕對值。"""
    return float(df["mean_nontumour"].abs().mean())


# ---------------------------------------------------------------- 第 0 層 ---
def layer0(R, seeds):
    rows, pair = [], {}
    for k in METHODS:
        rho_t, rho_n, vshare = [], [], []
        for s in seeds:
            d = R[s][k]
            ts, mt, mn = d.tumour_specific, d.mean_tumour, d.mean_nontumour
            rho_t.append(spearmanr(ts, mt).statistic)
            rho_n.append(spearmanr(ts, -mn).statistic)
            vshare.append(float(np.var(mn) / np.var(ts)))
        rows.append({"method": k,
                     "rho_ts_vs_tumour": float(np.mean(rho_t)),
                     "rho_ts_vs_neg_nontumour": float(np.mean(rho_n)),
                     "nontumour_var_share": float(np.mean(vshare))})
    dec = pd.DataFrame(rows).set_index("method")

    # 同一個種子內，四法兩兩的 Spearman：相減前 vs 相減後
    for col in ("tumour_specific", "mean_tumour"):
        acc = {}
        for i, a in enumerate(METHODS):
            for b in METHODS[i + 1:]:
                v = [spearmanr(R[s][a][col], R[s][b][col]).statistic for s in seeds]
                acc[f"{a}|{b}"] = float(np.mean(v))
        pair[col] = acc
    off_ts = float(np.mean(list(pair["tumour_specific"].values())))
    off_mt = float(np.mean(list(pair["mean_tumour"].values())))
    gain = off_mt - off_ts
    return {"decomposition": json.loads(dec.to_json(orient="index")),
            "pairwise_spearman": pair,
            "mean_offdiag_tumour_specific": off_ts,
            "mean_offdiag_mean_tumour": off_mt,
            "gain_from_dropping_subtrahend": gain,
            "prespecified_threshold": L0_GAIN,
            "verdict": ("disagreement_manufactured_by_subtrahend" if gain > L0_GAIN
                        else "subtrahend_not_the_main_driver")}, dec


# ---------------------------------------------------------------- 第 1 層 ---
def layer1(R, seeds):
    """對角 = 同方法跨種子；非對角 = 同種子跨方法；再做衰減校正。"""
    rel = {}
    for k in METHODS:
        v = [spearmanr(R[a][k].tumour_specific, R[b][k].tumour_specific).statistic
             for i, a in enumerate(seeds) for b in seeds[i + 1:]]
        rel[k] = float(np.mean(v)) if v else 1.0
    agree, agree_cross, dis = {}, {}, {}
    for i, a in enumerate(METHODS):
        for b in METHODS[i + 1:]:
            same = [spearmanr(R[s][a].tumour_specific, R[s][b].tumour_specific).statistic
                    for s in seeds]
            cross = [spearmanr(R[p][a].tumour_specific, R[q][b].tumour_specific).statistic
                     for p in seeds for q in seeds if p != q]
            key = f"{a}|{b}"
            agree[key] = float(np.mean(same))
            agree_cross[key] = float(np.mean(cross))
            # B 與 D 來自同一個 scVI fit → 同種子的一致度被共享 fit 抬高，改用跨種子
            num = agree_cross[key] if (a, b) == SHARED_FIT or (b, a) == SHARED_FIT else agree[key]
            den = np.sqrt(max(rel[a], 1e-9) * max(rel[b], 1e-9))
            dis[key] = float(num / den)
    worst_rel = min(rel[k] for k in STOCH)
    max_dis = max(dis.values())
    return {"reliability_within_method": rel,
            "agreement_same_seed": agree,
            "agreement_cross_seed": agree_cross,
            "agreement_disattenuated": dis,
            "min_reliability_stochastic": worst_rel,
            "max_disattenuated_agreement": max_dis,
            "prespecified": {"reliability_gt": L1_REL, "agreement_lt": L1_AGREE},
            "verdict": ("disagreement_systematic_not_seed_noise"
                        if (worst_rel > L1_REL and max_dis < L1_AGREE) else "inconclusive")}


# ---------------------------------------------------------------- 第 2 層 ---
def standardised(R, seeds):
    """y[method][seed] = tumour_specific / 該方法該種子自己的雜訊底。"""
    genes = R[seeds[0]][DET].index
    Y = {k: pd.DataFrame({s: R[s][k].tumour_specific / floor_of(R[s][k]) for s in seeds},
                         index=genes) for k in METHODS}
    return genes, Y


def layer2(R, seeds, genes, Y):
    n = len(seeds)
    k = len(STOCH)
    # 方法內：只用隨機方法，取 RMS（= sqrt(MS_within)），不是 SD 的平均
    S2 = np.stack([Y[m].var(axis=1, ddof=1).to_numpy() for m in STOCH])      # (k, G)
    ms_w = S2.mean(0)
    # 方法間：用各方法的「平均」（不是中位數），ddof=0（四個方法是固定集合，不是抽樣）
    Mn = np.stack([Y[m].mean(axis=1).to_numpy() for m in STOCH])
    ms_b = n * Mn.var(0, ddof=0) * k / max(k - 1, 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        F = np.where(ms_w > 0, ms_b / ms_w, np.nan)
    p = fdist.sf(F, k - 1, k * (n - 1))
    # BH-FDR
    ok = ~np.isnan(p)
    q = np.full_like(p, np.nan)
    pv = p[ok]
    o = np.argsort(pv)
    qv = np.empty_like(pv)
    qv[o] = np.minimum.accumulate((pv[o] * len(pv) / (np.arange(len(pv)) + 1))[::-1])[::-1]
    q[ok] = np.clip(qv, 0, 1)
    var_a = np.maximum(0.0, (ms_b - ms_w) / n)
    icc = np.where((var_a + ms_w) > 0, var_a / (var_a + ms_w), np.nan)

    out = pd.DataFrame({"gene": genes, "F": F, "p": p, "q": q,
                        "ms_between": ms_b, "ms_within": ms_w,
                        "sigma2_method": var_a, "icc_method": icc})
    for m in METHODS:
        out[f"y_{m[0]}"] = Y[m].mean(axis=1).to_numpy()
        out[f"sd_{m[0]}"] = Y[m].std(axis=1, ddof=1).to_numpy()
    out = out.set_index("gene")

    # A 當成固定參考，各隨機方法對它做單樣本 t（「A 和 B 的差距是否大於 B 的種子雜訊」）
    contrasts = {}
    for m in STOCH:
        d = Y[m].mean(axis=1) - Y[DET].mean(axis=1)
        se = Y[m].std(axis=1, ddof=1) / np.sqrt(n)
        tt = d / se.replace(0, np.nan)
        contrasts[m] = {"t": float(tt.get(FOCUS, np.nan)),
                        "df": n - 1,
                        "p": float(2 * tdist.sf(abs(tt.get(FOCUS, np.nan)), n - 1))
                        if np.isfinite(tt.get(FOCUS, np.nan)) else None,
                        "delta_vs_A": float(d.get(FOCUS, np.nan))}

    # 敏感度：B 與 D 共享 fit，拿掉其中一個看結論是否還在
    sens = {}
    for drop in ("B_scvi_counterfactual", "D_latent_arithmetic"):
        keep = [m for m in STOCH if m != drop]
        s2 = np.stack([Y[m].loc[FOCUS].var(ddof=1) for m in keep])
        mn = np.array([Y[m].loc[FOCUS].mean() for m in keep])
        w = float(s2.mean()); b = float(n * mn.var(ddof=0) * 2 / 1)
        sens[f"without_{drop[0]}"] = {"F": (b / w) if w > 0 else None,
                                      "df": [1, 2 * (n - 1)]}

    fo = out.loc[FOCUS]
    pct_F = float((out["F"] < fo["F"]).mean() * 100)
    pct_icc = float((out["icc_method"] < fo["icc_method"]).mean() * 100)
    below_floor = {m: bool(abs(Y[m].loc[FOCUS].mean()) < 1.0) for m in METHODS}

    focus = {"gene": FOCUS,
             "y_by_method_by_seed": {m: {str(s): float(Y[m].loc[FOCUS, s]) for s in seeds}
                                     for m in METHODS},
             "y_mean": {m: float(Y[m].loc[FOCUS].mean()) for m in METHODS},
             "y_sd": {m: float(Y[m].loc[FOCUS].std(ddof=1)) for m in METHODS},
             "F": float(fo["F"]), "p": float(fo["p"]), "q": float(fo["q"]),
             "icc_method": float(fo["icc_method"]),
             "genomewide_percentile_F": pct_F,
             "genomewide_percentile_icc": pct_icc,
             "in_top_decile_of_method_dependence": bool(pct_F >= 90.0),
             "below_own_noise_floor": below_floor,
             "contrasts_vs_A": contrasts,
             "sensitivity_shared_fit": sens}

    gw = {"frac_q_lt_0.05": float((out["q"] < 0.05).mean()),
          "frac_sigma2_method_zero": float((out["sigma2_method"] <= 0).mean()),
          "median_icc": float(out["icc_method"].median())}
    pos = {g: {"y_mean": {m: float(Y[m].loc[g].mean()) for m in METHODS},
               "F": float(out.loc[g, "F"]), "q": float(out.loc[g, "q"]),
               "genomewide_percentile_F": float((out["F"] < out.loc[g, "F"]).mean() * 100)}
           for g in POS if g in out.index}
    panel = {g: {"y_mean": {m: round(float(Y[m].loc[g].mean()), 3) for m in METHODS},
                 "q": float(out.loc[g, "q"]),
                 "below_floor_in": [m for m in METHODS if abs(Y[m].loc[g].mean()) < 1.0]}
             for g in PANEL if g in out.index}
    return out, {"focus": focus, "genomewide": gw, "positive_controls": pos, "panel": panel}


def main():
    R = load()
    if not R:
        raise SystemExit("找不到任何 130_residual_by_method_s*.csv.gz")
    seeds = sorted(R)
    expl = [s for s in seeds if s in EXPLORE]
    held = [s for s in seeds if s not in EXPLORE]
    log("可用種子:", seeds)
    log("宣告集（已被檢視過）:", expl, "| 保留集（從未檢視）:", held)
    if len(expl) < 3:
        raise SystemExit("宣告集至少要 3 個種子")

    report = {"seeds_all": seeds, "seeds_declare": expl, "seeds_heldout": held,
              "n_nontumour_cells": int(R[seeds[0]][DET]["n_nontumour"].iloc[0]),
              "n_tumour_cells": int(R[seeds[0]][DET]["n_tumour"].iloc[0]),
              "donor_composition_nontumour": {"M413": 790, "M102": 29,
                                              "note": "790/819 = 96.5% 來自單一捐贈者"},
              "noise_floor_by_method_seed": {
                  m: {str(s): floor_of(R[s][m]) for s in seeds} for m in METHODS}}

    for name, ss in (("declare", expl), ("heldout", held)):
        # 少於 3 個種子時方法內變異只有 ≤1 自由度，F 與衰減校正都沒有意義 → 不下結論
        if len(ss) < 3:
            report[name] = {"seeds": ss, "skipped": True,
                            "note": "種子不足 3 個，方法內變異自由度 ≤1，不計算"}
            log("%s 只有 %d 個種子（%s），不足 3 個，跳過不下結論" % (name, len(ss), ss))
            continue
        l0, dec = layer0(R, ss)
        l1 = layer1(R, ss)
        genes, Y = standardised(R, ss)
        tab, l2 = layer2(R, ss, genes, Y)
        report[name] = {"seeds": ss, "layer0": l0, "layer1": l1, "layer2": l2}
        if name == "declare":
            tab.to_csv(f"{RES}/220_per_gene_declare.csv.gz", compression="gzip")
            log("\n=== 第 0 層：tumour_specific 的兩項拆解（宣告集）===")
            log(dec.round(3).to_string())
            log("\n四法平均非對角 Spearman：相減後 %.3f → 只用腫瘤項 %.3f（增益 %+.3f，門檻 %.2f）"
                % (l0["mean_offdiag_tumour_specific"], l0["mean_offdiag_mean_tumour"],
                   l0["gain_from_dropping_subtrahend"], L0_GAIN))
            log("L0 判定：%s" % l0["verdict"])
            log("\n=== 第 1 層：再測信度（對角）vs 方法間一致度（非對角）===")
            for m in METHODS:
                log("  信度 %-24s %.4f" % (m, l1["reliability_within_method"][m]))
            for k_, v in l1["agreement_disattenuated"].items():
                log("  衰減校正後 %-50s %.3f" % (k_, v))
            log("L1 判定：%s" % l1["verdict"])
            log("\n=== 第 2 層：%s（單位 = 各方法自己的雜訊底倍數）===" % FOCUS)
            f = l2["focus"]
            for m in METHODS:
                log("  %-24s %s  平均 %.2f ± %.2f%s"
                    % (m, " ".join("%6.2f" % f["y_by_method_by_seed"][m][str(s)] for s in ss),
                       f["y_mean"][m], f["y_sd"][m],
                       "   ← 在自己的雜訊底以下" if f["below_own_noise_floor"][m] else ""))
            log("  F(%d,%d) = %.2f, p = %.2e, q = %.2e | ICC_method = %.3f"
                % (len(STOCH) - 1, len(STOCH) * (len(ss) - 1), f["F"], f["p"], f["q"],
                   f["icc_method"]))
            log("  全基因組百分位：F %.1f、ICC %.1f → 前十分位 = %s"
                % (f["genomewide_percentile_F"], f["genomewide_percentile_icc"],
                   f["in_top_decile_of_method_dependence"]))
            log("  全基因組 q<0.05 的比例 = %.3f（方法效應是全域性質）"
                % l2["genomewide"]["frac_q_lt_0.05"])
            for g, v in l2["positive_controls"].items():
                log("  陽性對照 %-6s F 百分位 %.1f, q = %.2e, y = %s"
                    % (g, v["genomewide_percentile_F"], v["q"],
                       {m[0]: round(x, 1) for m, x in v["y_mean"].items()}))

    # ------------------------------------------------------------ 確認 ---
    conf = {}
    if held and not report.get("heldout", {}).get("skipped", False):
        for key, path in (("L0", ("layer0", "verdict")), ("L1", ("layer1", "verdict")),
                          ("L2_top_decile", ("layer2", "focus",
                                             "in_top_decile_of_method_dependence"))):
            a = report["declare"]
            b = report["heldout"]
            for p_ in path:
                a, b = a[p_], b[p_]
            conf[key] = {"declare": a, "heldout": b, "reproduced": bool(a == b)}
        conf["all_reproduced"] = all(v["reproduced"] for v in conf.values()
                                     if isinstance(v, dict))
        log("\n=== 保留集確認（種子 %s）===" % held)
        for k_, v in conf.items():
            if isinstance(v, dict):
                log("  %-16s 宣告=%s 保留=%s 重現=%s"
                    % (k_, v["declare"], v["heldout"], v["reproduced"]))
        log("  全部重現 = %s" % conf["all_reproduced"])
    report["confirmation"] = conf

    json.dump(report, open(f"{RES}/220_method_decomposition.json", "w"),
              indent=2, ensure_ascii=False)
    log("\n→ results/220_method_decomposition.json、220_per_gene_declare.csv.gz")
    return 0


if __name__ == "__main__":
    sys.exit(main())
