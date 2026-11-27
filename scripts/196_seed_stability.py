#!/usr/bin/env python
"""步驟 196：四種殘差方法「互不同意」到底是隨機種子的雜訊，還是真的？

為什麼要做這一題
----------------
步驟 130 用四種做法重算同一件事，結果彼此相當不同意：
top-100 的 Jaccard 只有 0.081–0.282，201 個基因裡只有 25 個四法全中。
SLC7A11 在四法的排名是 A #158 / B #150 / C #414 / D #53（3/4 進前 5%）。
但 130 只跑了**一個種子**（scvi.settings.seed = 0），所以那份「不同意」有兩種
完全不同的解釋，而它們對已發表的宣稱有相反的意涵：

  解釋一：不同意是**種子雜訊**。生成模型每次訓練都不一樣，四個數字只是同一個
          分布的四次抽樣。→ 方法之間的差異不該被當成證據，A（確定性、步驟 15
          的做法）的結果可以當代表。
  解釋二：不同意是**真的**。每個方法自己很穩，但它們系統性地指向不同基因。
          → 「MCL 專一地開啟 SLC7A11」是**方法依賴**的宣稱，報告裡必須寫明
          它站在哪一個方法上，而不是說「四種方法都做了」。

可判定的量
----------
同一個基因、同一個方法、不同種子 → 排名的離散程度 = **方法內雜訊**
同一個基因、不同方法（各取其種子中位數） → 排名的離散程度 = **方法間差異**
比值 R = 方法間 / 方法內。這兩個量單位相同（排名），所以可以直接比。

事先寫死的判定（跑完之前就固定，避免看到數字再挑規則）
-------------------------------------------------------
P0 決定性：種子 0 重跑 vs 已發表的 130_method_comparison.json。
   方法 A 沒有任何 RNG，**必須完全一致**；不一致代表流程本身有問題，整個測試作廢。
   B/C/D 在固定種子下若 |Δ排名| > 25，代表連固定種子都不決定性（cuDNN 等），
   這件事本身要報告。
P1 陽性對照：CCND1 與 SOX11 必須在**每一個方法 × 每一個種子**都進前 50。
   這兩個是 t(11;14) 易位產物與 MCL 的定義性標記，測不出來就是跑壞了，不下結論。
P2 本題：SLC7A11 的 R = 方法間離散 / 方法內離散。
      R > 3.0   → 不同意是真的，宣稱是方法依賴的，報告必須改寫
      R < 1.5   → 不同意是種子雜訊，A 的結果可當代表
      1.5–3.0   → 判不出來，照實寫「判不出來」
P3 陰性對照：取在種子 0 時四個方法都落在 |tumour_specific| 最低十分位的基因。
   這些是雜訊基因，它們的排名**應該**在種子之間亂跑。若它們的方法內離散
   沒有明顯大於 CCND1/SOX11（要求中位數至少 5 倍），代表這個排名指標被
   表現量之類的東西主導、而不是被訊號主導 → P2 的比值不可解釋，不下結論。
P4 已發表措辭的穩健性：SLC7A11「3/4 方法同意」（排名進前 5% = 前 201 名）
   在各種子下的同意數。若跨種子變動超過 ±1，報告裡的「3 of 4」必須改成區間。

輸出 results/196_seed_stability.json
"""
import json
import sys

import numpy as np
import pandas as pd

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
RES = f"{BASE}/results"
SEEDS = [0, 1, 2, 3, 4]
METHODS = ["A_gene_space", "B_scvi_counterfactual", "C_contrastivevi", "D_latent_arithmetic"]
POS = ["CCND1", "SOX11"]          # P1 陽性對照
FOCUS = "SLC7A11"                 # P2 的主角
TOP_FRAC = 0.05                   # P4：前 5%
PANEL = ["SLC7A11", "SLC3A2", "GCLM", "GCLC", "AIFM2", "MAFG", "MAFF", "MAFK",
         "BACH1", "NFE2L2", "CCND1", "SOX11", "TXN", "GPX4", "HMOX1", "NQO1"]


def log(*a):
    print("[196]", *a, flush=True)


def ranks_for(seed):
    """回傳 genes × methods 的排名表（1 = tumour_specific 最高），與 130 的定義相同。"""
    p = f"{RES}/130_tumour_specific_by_method_s{seed}.csv"
    piv = pd.read_csv(p, index_col=0)
    miss = [m for m in METHODS if m not in piv.columns]
    if miss:
        raise SystemExit(f"種子 {seed} 缺方法 {miss}")
    # 130 用 (piv[m] > piv[m][g]).sum() + 1 ⇒ 等同 method="min" 的降序排名
    return piv[METHODS].rank(ascending=False, method="min").astype(int), piv[METHODS]


def main():
    R, V = {}, {}
    for s in SEEDS:
        try:
            R[s], V[s] = ranks_for(s)
        except FileNotFoundError:
            log(f"種子 {s} 的輸出還沒有，跳過")
    if len(R) < 3:
        raise SystemExit(f"只有 {len(R)} 個種子完成，至少要 3 個才有離散度可談")
    seeds = sorted(R)
    genes = R[seeds[0]].index
    for s in seeds[1:]:
        if not R[s].index.equals(genes):
            raise SystemExit("不同種子的基因集不一致，無法比較")
    log(f"完成的種子：{seeds} | 基因數 {len(genes)}")

    out = {"seeds_done": seeds, "n_genes": int(len(genes)), "gates": {}, "verdict": None}

    # ---------------------------------------------------------- P0 決定性 ---
    pub = json.load(open(f"{RES}/130_method_comparison.json"))
    p0 = {}
    if 0 in R:
        for m in METHODS:
            pr = pub["snr"][m]["panel_rank"]
            d = {g: int(R[0].loc[g, m]) - int(pr[g]) for g in PANEL if g in R[0].index and g in pr}
            p0[m] = {"max_abs_delta_rank": int(max(abs(v) for v in d.values())),
                     "delta_by_gene": d}
        a_exact = p0["A_gene_space"]["max_abs_delta_rank"] == 0
        nondet = {m: p0[m]["max_abs_delta_rank"] > 25 for m in METHODS if m != "A_gene_space"}
        out["gates"]["P0_determinism"] = {
            "A_reproduces_published_exactly": bool(a_exact),
            "max_abs_delta_rank": {m: p0[m]["max_abs_delta_rank"] for m in METHODS},
            "nondeterministic_at_fixed_seed": nondet,
            "pass": bool(a_exact),
            "note": ("方法 A 沒有 RNG，完全重現 = 流程可信"
                     if a_exact else "方法 A 竟然和已發表的不同 → 流程有問題，測試作廢"),
        }
        if not a_exact:
            out["verdict"] = "VOID_pipeline_not_reproducible"
            json.dump(out, open(f"{RES}/196_seed_stability.json", "w"),
                      indent=2, ensure_ascii=False)
            log("P0 失敗：方法 A 沒有重現已發表的排名。停止，不下任何結論。")
            log(json.dumps(p0["A_gene_space"]["delta_by_gene"], ensure_ascii=False))
            return 1
    else:
        out["gates"]["P0_determinism"] = {"pass": None, "note": "種子 0 沒跑，無法檢查決定性"}

    # ------------------------------------------------------- P1 陽性對照 ---
    bad = []
    for g in POS:
        for m in METHODS:
            for s in seeds:
                r = int(R[s].loc[g, m])
                if r > 50:
                    bad.append({"gene": g, "method": m, "seed": s, "rank": r})
    out["gates"]["P1_positive_control"] = {
        "requirement": "CCND1 與 SOX11 在每個方法 × 每個種子都要進前 50",
        "violations": bad,
        "pass": len(bad) == 0,
    }
    if bad:
        out["verdict"] = "VOID_positive_control_failed"
        json.dump(out, open(f"{RES}/196_seed_stability.json", "w"), indent=2, ensure_ascii=False)
        log(f"P1 失敗：{len(bad)} 次違規。跑壞了，不下結論。")
        for b in bad[:10]:
            log("  ", b)
        return 1

    # 方法內離散（跨種子）／方法間離散（各方法取種子中位數）
    def spreads(g):
        within = {m: float(np.std([R[s].loc[g, m] for s in seeds], ddof=1)) for m in METHODS}
        med = {m: float(np.median([R[s].loc[g, m] for s in seeds])) for m in METHODS}
        between = float(np.std(list(med.values()), ddof=1))
        wbar = float(np.mean(list(within.values())))
        return within, med, between, wbar

    # ------------------------------------------------------- P3 陰性對照 ---
    s0 = seeds[0]
    absv = V[s0][METHODS].abs()
    dec = absv.quantile(0.10)
    null_genes = [g for g in genes if all(absv.loc[g, m] <= dec[m] for m in METHODS)]
    rng = np.random.default_rng(0)
    if len(null_genes) > 50:
        null_genes = list(rng.choice(null_genes, 50, replace=False))
    log(f"陰性對照基因 {len(null_genes)} 個（四法皆落在 |殘差| 最低十分位）")
    nw = [np.mean(list(spreads(g)[0].values())) for g in null_genes] if null_genes else []
    pw = [np.mean(list(spreads(g)[0].values())) for g in POS]
    ratio_np = float(np.median(nw) / max(np.mean(pw), 1e-9)) if nw else float("nan")
    out["gates"]["P3_negative_control"] = {
        "requirement": "雜訊基因的方法內排名離散，中位數要 ≥ 陽性對照的 5 倍",
        "n_null_genes": len(null_genes),
        "median_within_sd_null": float(np.median(nw)) if nw else None,
        "mean_within_sd_positive": float(np.mean(pw)),
        "ratio": ratio_np,
        "pass": bool(ratio_np >= 5.0),
        "note": ("排名指標被訊號主導，比值可解釋" if ratio_np >= 5.0 else
                 "雜訊基因的排名和陽性對照一樣穩 → 指標被訊號以外的東西主導，P2 不可解釋"),
    }

    # ------------------------------------------------------------ P2 本題 ---
    within, med, between, wbar = spreads(FOCUS)
    Rratio = float(between / max(wbar, 1e-9))
    if Rratio > 3.0:
        p2 = "disagreement_is_real"
    elif Rratio < 1.5:
        p2 = "disagreement_is_seed_noise"
    else:
        p2 = "undetermined"
    out["gates"]["P2_main"] = {
        "gene": FOCUS,
        "rank_by_seed": {m: {str(s): int(R[s].loc[FOCUS, m]) for s in seeds} for m in METHODS},
        "within_method_sd": within,
        "mean_within_method_sd": wbar,
        "median_rank_per_method": med,
        "between_method_sd": between,
        "R_between_over_within": Rratio,
        "thresholds": {">3.0": "disagreement_is_real", "<1.5": "seed_noise"},
        "result": p2,
    }

    # ----------------------------------------------- P4 「3/4 同意」穩健性 ---
    top_n = int(round(TOP_FRAC * len(genes)))
    agree = {str(s): int(sum(int(R[s].loc[FOCUS, m]) <= top_n for m in METHODS)) for s in seeds}
    vals = list(agree.values())
    out["gates"]["P4_agreement_wording"] = {
        "top_n": top_n,
        "n_methods_agreeing_by_seed": agree,
        "published_claim": 3,
        "range": [min(vals), max(vals)],
        "stable": bool(max(vals) - min(vals) <= 1),
        "note": ("已發表的「3/4」措辭站得住" if max(vals) - min(vals) <= 1
                 else "跨種子變動超過 ±1 → 報告必須改成區間"),
    }

    # 整個 panel 的同一組量，供報告引用
    out["panel"] = {}
    for g in PANEL:
        if g not in genes:
            continue
        w, md, b, wb = spreads(g)
        out["panel"][g] = {"within_method_sd": w, "median_rank_per_method": md,
                           "between_method_sd": b, "mean_within_method_sd": wb,
                           "R": float(b / max(wb, 1e-9))}

    # 最終判定：只有在 P3 通過時 P2 才算話
    if not out["gates"]["P3_negative_control"]["pass"]:
        out["verdict"] = "P2_uninterpretable_metric_fails_negative_control"
    else:
        out["verdict"] = p2

    json.dump(out, open(f"{RES}/196_seed_stability.json", "w"), indent=2, ensure_ascii=False)

    log("\n=== SLC7A11 跨種子排名 ===")
    log("%-26s" % "method" + "".join("%9s" % f"s{s}" for s in seeds) + "%12s%12s"
        % ("種子 SD", "中位數"))
    for m in METHODS:
        log("%-26s" % m + "".join("%9d" % R[s].loc[FOCUS, m] for s in seeds)
            + "%12.1f%12.0f" % (within[m], med[m]))
    log("\n方法內離散平均 = %.1f 名 | 方法間離散 = %.1f 名 | R = %.2f" % (wbar, between, Rratio))
    log("P0 方法 A 完全重現 = %s" % out["gates"]["P0_determinism"].get(
        "A_reproduces_published_exactly"))
    log("P1 陽性對照 = 通過")
    log("P3 陰性對照 比值 %.1f（要 ≥5）= %s" % (ratio_np, out["gates"]["P3_negative_control"]["pass"]))
    log("P4 同意數 %s = %s" % (agree, out["gates"]["P4_agreement_wording"]["note"]))
    log("\n判定：%s" % out["verdict"])
    log("→ results/196_seed_stability.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
