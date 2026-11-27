#!/usr/bin/env python
"""步驟 237：small Maf 的冗餘是「部分」還是「完全」？用互惠性來判。

為什麼 234 還不夠
-----------------
步驟 234 S3 只測了一個方向：`MAFF`+`MAFK` 表現高 → 較不依賴 `MAFG`
（r = +0.0776、經驗 p = 0.0190）。但**單向相關不足以區分兩種情況**：

  (a) **互相緩衝（接近完全冗餘）**：三個成員可以彼此替代 →
      反方向也該成立（`MAFG` 表現高 → 較不依賴 `MAFK`）。
      → 只打一個沒用，必須全打。
  (b) **單向／不對稱**：只有 `MAFG` 會被別人補，反之不然 →
      冗餘的程度與方向都受限，打 `MafG` 仍可能有效。

這個差別決定整條路線生死，而 DepMap 的現有資料就能判 —— 需要的是**互惠性檢定**。

事先寫死的檢定與判準（在讀任何數字之前寫成）
--------------------------------------------
B1 **正向（重現 234，當守門）**：`MAFF`+`MAFK` 表現 vs `MAFG` gene effect，
   Spearman r > 0 且參數 p < 0.05。不成立 → 與 234 矛盾，整支作廢。
B2 **反向（本支的新檢定）**：`MAFG` 表現 vs `MAFK` gene effect，以及 vs `MAFF` gene effect。
   規則：r > 0、p < 0.05、且相對 1,000 個隨機基因的經驗 p < 0.05。
B3 **互惠性判定**：
   B1 與 B2 **都**成立 → `mutual_buffering`（接近完全冗餘，必須全打）
   只有 B1 成立       → `asymmetric_MAFG_is_buffered`（不對稱，打 MafG 仍可能有效）
   只有 B2 成立       → `asymmetric_reverse_only`
   都不成立           → `no_reciprocity`
B4 **分層效應量（不只看相關係數）**：把株按 `MAFK` 表現分三分位，
   比較最低 vs 最高三分位的 `MAFG` gene effect（Mann-Whitney + 中位數差）。
   反向同做。相關係數很小的時候，效應量才看得出實際意義。
B5 **共依賴性**：三個基因的 gene effect 輪廓彼此的 Spearman。
   同一個複合體 → 正相關；互相緩衝 → 不一定。這一項**描述性**，不當判準。
B6 **陰性對照**：隨機挑 20 對無關基因做同一套互惠性檢定，
   通過率不得明顯高於 5%（否則這個檢定本身會亂報）。

事先預測（記入誠實清單）
------------------------
- B1 會過（234 已經測過）。
- B2 我預測**會過但更弱** —— `MAFG` 是三者中表現最高、最廣的，
  所以它緩衝別人的能力應該比別人緩衝它更強，反向訊號會存在但小。
- B3 我預測 `mutual_buffering`。
- B4 我預測效應量很小（中位數差 < 0.05），也就是「統計上有、生物上小」。

輸出 results/237_*
"""
import json
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
DM = f"{BASE}/data/depmap"
RES = f"{BASE}/results"
FAM = ["MAFG", "MAFF", "MAFK"]
N_NULL = 1000
N_NEGPAIR = 20
SEED = 0


def log(*a):
    print("[237]", *a, flush=True)


def colmap(hdr, genes):
    s2c = {}
    for c in hdr.columns[1:]:
        s2c.setdefault(c.split(" (")[0], c)
    return {g: s2c[g] for g in genes if g in s2c}


def load(fn, genes):
    hdr = pd.read_csv(fn, nrows=0)
    cm = colmap(hdr, genes)
    D = pd.read_csv(fn, usecols=[hdr.columns[0]] + list(cm.values())).set_index(hdr.columns[0])
    D.columns = [c.split(" (")[0] for c in D.columns]
    return D


def rho(x, y):
    m = x.notna() & y.notna()
    if m.sum() < 100:
        return np.nan, np.nan, int(m.sum())
    r = spearmanr(x[m], y[m])
    return float(r.statistic), float(r.pvalue), int(m.sum())


def main():
    E = load(f"{DM}/CRISPRGeneEffect.csv", FAM)
    X = load(f"{DM}/OmicsExpressionProteinCodingGenesTPMLogp1.csv", FAM)
    common = E.index.intersection(X.index)
    E, X = E.loc[common], X.loc[common]
    log("共同細胞株 %d" % len(common))
    out = {"n_lines": int(len(common))}

    # 經驗虛無用的隨機基因池（一次抽，B2 與 B6 共用）
    allx = pd.read_csv(f"{DM}/OmicsExpressionProteinCodingGenesTPMLogp1.csv", nrows=0)
    rng = np.random.default_rng(SEED)
    pool = [c for c in allx.columns[1:] if c.split(" (")[0] not in set(FAM)]
    sel = list(rng.choice(pool, size=min(N_NULL, len(pool)), replace=False))
    log("抽 %d 個隨機基因的表現量…" % len(sel))
    blocks = []
    for a in range(0, len(sel), 200):
        blocks.append(pd.read_csv(f"{DM}/OmicsExpressionProteinCodingGenesTPMLogp1.csv",
                                  usecols=[allx.columns[0]] + sel[a:a + 200]
                                  ).set_index(allx.columns[0]).reindex(common))
    NX = pd.concat(blocks, axis=1)
    log("虛無表現矩陣 %s" % (NX.shape,))

    def emp_p(dep, obs_r):
        vals = []
        for c in NX.columns:
            r, _, n = rho(NX[c], dep)
            if np.isfinite(r):
                vals.append(r)
        vals = np.asarray(vals)
        return float((vals >= obs_r).mean()), vals

    # ------------------------------------------------------------ B1 ---
    partner = X[["MAFF", "MAFK"]].sum(1)
    r1, p1, n1 = rho(partner, E["MAFG"])
    b1 = bool(r1 > 0 and p1 < 0.05)
    out["B1_forward"] = {"rule": "MAFF+MAFK 表現 vs MAFG gene effect，r > 0 且 p < 0.05",
                         "spearman_r": r1, "pvalue": p1, "n": n1, "pass": b1}
    log("\n=== B1 正向（守門，重現 234）===")
    log("  MAFF+MAFK 表現 vs MAFG gene effect：r = %+.4f, p = %.3g, n = %d → %s"
        % (r1, p1, n1, "過" if b1 else "不過"))
    if not b1:
        out["verdict"] = "VOID_contradicts_234"
        json.dump(out, open(f"{RES}/237_paralog_buffering.json", "w"),
                  indent=2, ensure_ascii=False)
        return 1

    # ------------------------------------------------------------ B2 ---
    log("\n=== B2 反向（新檢定）：MAFG 表現 vs 另外兩個的 gene effect ===")
    rev = {}
    for tgt in ["MAFK", "MAFF"]:
        r, p, n = rho(X["MAFG"], E[tgt])
        ep, _ = emp_p(E[tgt], r)
        ok = bool(r > 0 and p < 0.05 and ep < 0.05)
        rev[tgt] = {"spearman_r": r, "pvalue": p, "n": n,
                    "empirical_p": ep, "pass": ok}
        log("  MAFG 表現 vs %s gene effect：r = %+.4f, p = %.3g, 經驗 p = %.4f → %s"
            % (tgt, r, p, ep, "過" if ok else "不過"))
    b2 = bool(any(v["pass"] for v in rev.values()))
    out["B2_reverse"] = {"rule": "MAFG 表現 vs MAFK/MAFF gene effect，"
                                 "r > 0、p < 0.05、經驗 p < 0.05",
                         "per_target": rev, "pass_any": b2}

    # ------------------------------------------------------------ B3 ---
    if b1 and b2:
        v = "mutual_buffering"
    elif b1:
        v = "asymmetric_MAFG_is_buffered"
    elif b2:
        v = "asymmetric_reverse_only"
    else:
        v = "no_reciprocity"
    out["B3_reciprocity"] = {"call": v}
    log("\n=== B3 互惠性判定：%s ===" % v)
    if v == "mutual_buffering":
        log("  → 接近完全冗餘：三個成員互相替代，**只打一個不會有用，必須全打**。")
    elif v == "asymmetric_MAFG_is_buffered":
        log("  → 不對稱：`MAFG` 被別人補，但反之不成立。")
        log("     冗餘的方向受限，打 `MafG` 仍可能有效，但要同時壓住補位的那一個。")

    # ------------------------------------------------------------ B4 ---
    log("\n=== B4 分層效應量（相關係數小的時候，看實際差多少）===")
    strat = {}
    for expr_g, dep_g in [("MAFK", "MAFG"), ("MAFF", "MAFG"), ("MAFG", "MAFK")]:
        v_ = X[expr_g]
        d_ = E[dep_g]
        m = v_.notna() & d_.notna()
        lo_t, hi_t = v_[m].quantile(1 / 3), v_[m].quantile(2 / 3)
        lo = d_[m & (v_ <= lo_t)]
        hi = d_[m & (v_ >= hi_t)]
        u = mannwhitneyu(lo, hi, alternative="less")
        strat["%s_expr_vs_%s_dep" % (expr_g, dep_g)] = {
            "n_low": int(len(lo)), "n_high": int(len(hi)),
            "median_dep_low_expr": float(lo.median()),
            "median_dep_high_expr": float(hi.median()),
            "median_diff": float(hi.median() - lo.median()),
            "mannwhitney_p": float(u.pvalue)}
        log("  %-6s 表現低 vs 高 → %-6s gene effect 中位 %+.4f vs %+.4f"
            "（差 %+.4f，MW p = %.3g）"
            % (expr_g, dep_g, lo.median(), hi.median(),
               hi.median() - lo.median(), u.pvalue))
    out["B4_stratified_effect"] = strat

    # ------------------------------------------------------------ B5 ---
    log("\n=== B5 共依賴性（描述性，不當判準）===")
    co = {}
    for i, a in enumerate(FAM):
        for b in FAM[i + 1:]:
            r, p, n = rho(E[a], E[b])
            co["%s_vs_%s" % (a, b)] = {"spearman_r": r, "pvalue": p, "n": n}
            log("  %s vs %s 的 gene effect 輪廓：r = %+.4f, p = %.3g" % (a, b, r, p))
    out["B5_codependency"] = co

    # ------------------------------------------------------------ B6 ---
    log("\n=== B6 陰性對照：%d 對隨機基因跑同一套互惠性檢定 ===" % N_NEGPAIR)
    he = pd.read_csv(f"{DM}/CRISPRGeneEffect.csv", nrows=0)
    dep_syms = [c.split(" (")[0] for c in he.columns[1:]]
    expr_syms = set(c.split(" (")[0] for c in allx.columns[1:])
    cand = [g for g in dep_syms if g in expr_syms and g not in set(FAM)]
    pairs = [(rng.choice(cand), rng.choice(cand)) for _ in range(N_NEGPAIR * 3)]
    pairs = [(a, b) for a, b in pairs if a != b][:N_NEGPAIR]
    npass = 0
    details = []
    for a, b in pairs:
        try:
            Ei = load(f"{DM}/CRISPRGeneEffect.csv", [a, b]).reindex(common)
            Xi = load(f"{DM}/OmicsExpressionProteinCodingGenesTPMLogp1.csv",
                      [a, b]).reindex(common)
            rf, pf, _ = rho(Xi[b], Ei[a])
            rr, pr, _ = rho(Xi[a], Ei[b])
            ok = bool(rf > 0 and pf < 0.05 and rr > 0 and pr < 0.05)
            npass += ok
            details.append({"pair": [a, b], "r_fwd": rf, "p_fwd": pf,
                            "r_rev": rr, "p_rev": pr, "reciprocal": ok})
        except Exception:
            continue
    rate = npass / max(len(details), 1)
    b6 = bool(rate <= 0.20)
    out["B6_negative_control"] = {
        "rule": "隨機基因對的互惠通過率不得明顯高於 5%（容忍到 20%）",
        "n_pairs": len(details), "n_reciprocal": npass,
        "rate": rate, "pass": b6, "details": details[:8]}
    log("  %d / %d 對隨機基因也呈互惠（%.0f%%）→ %s"
        % (npass, len(details), rate * 100,
           "檢定可信" if b6 else "**檢定本身會亂報，互惠結論不可信**"))

    out["verdict"] = {"reciprocity": v,
                      "negative_control_ok": b6,
                      "call": v if b6 else "VOID_test_not_specific"}
    log("\n*** 最終判定：%s ***" % out["verdict"]["call"])
    json.dump(out, open(f"{RES}/237_paralog_buffering.json", "w"),
              indent=2, ensure_ascii=False)
    log("\n→ results/237_paralog_buffering.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
