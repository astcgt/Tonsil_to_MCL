#!/usr/bin/env python
"""步驟 216：「沒有備援手臂」—— `SLC7A11` 的 B 細胞限制終於有機制了嗎？

三條獨立線索收斂到同一個地方
----------------------------
鐵死亡有兩條互相獨立的防線：
  **GPX4 手臂**  需要兩種燃料：穀胱甘肽（來自 xCT → `SLC7A11`）與**硒**（GPX4 是硒蛋白）
  **FSP1 手臂**  `AIFM2`，走 CoQ10，完全不需要穀胱甘肽

今天三個互相獨立的分析指向同一件事：
  (1) 步驟 204（子代理）：18 個候選特徵裡，`exp_GPX4` 是唯一通過**最嚴分層**
      （B 系非漿細胞）的生物學特徵
  (2) 步驟 205（子代理）：和 `SLC7A11` 共依賴的倖存者首位是 **`SCLY`**
      （硒半胱胺酸裂解酶 —— 硒的來源，也就是 GPX4 的另一種燃料）。
      陽性對照有效：`GCLM` 排第 2、`GCLC` 排第 4
  (3) 步驟 214（本人）：`AIFM2`/FSP1 在 MCL（0.0277）與分化配對正常 B（0.0337）
      **都近乎不存在**

合起來的假說很簡單：**B 細胞沒有 FSP1 備援，所以只剩 GPX4 一條路，
於是對 GPX4 的兩種燃料都無可取代地依賴。**
這和今天被否證的轉硫途徑假說不同 —— 它不預測「B 細胞比較特別」，
而預測「**任何**缺乏 FSP1 的細胞都會這樣」，所以可以在非淋巴系株裡檢定。

事前指定
--------
P1  `AIFM2` 表現在 B 細胞腫瘤株顯著低於非淋巴系株。
P2  **關鍵的一關**：在**非淋巴系**株內，`AIFM2` 表現越低，`SLC7A11` 越依賴。
    這一關決定它是機制還是只是「是淋巴系」的同義詞。
P3  效應量要大於隨機基因的虛無分佈（1,000 個隨機基因表現對
    `SLC7A11` 依賴性的相關）。
P4（陰性對照）`AIFM2` 表現不該同樣預測一個無關必需基因（`RPL13A`）的依賴性。
P5  複合指標：同時 `AIFM2` 低且 `GPX4` 依賴的株，`SLC7A11` 依賴率應最高。
    用 2×2 分層看，而且必須在非淋巴系株內也成立。

輸出 results/216_*
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
N_NULL = 1000
RNG = np.random.default_rng(0)


def log(*a):
    print("[216]", *a, flush=True)


def take(path, genes):
    hdr = pd.read_csv(path, nrows=0)
    f = hdr.columns[0]
    m = {}
    for c in hdr.columns[1:]:
        m.setdefault(c.split(" (")[0], c)
    have = [g for g in genes if g in m]
    D = pd.read_csv(path, usecols=[f] + [m[g] for g in have]).set_index(f)
    D.columns = [c.split(" (")[0] for c in D.columns]
    return D, have, hdr, f


def main():
    EFFG = ["SLC7A11", "GPX4", "AIFM2", "SCLY", "RPL13A"]
    EXPG = ["AIFM2", "GPX4", "SCLY", "SLC7A11"]
    E, he, ehdr, ef = take(f"{DM}/CRISPRGeneEffect.csv", EFFG)
    X, hx, xhdr, xf = take(f"{DM}/OmicsExpressionProteinCodingGenesTPMLogp1.csv", EXPG)
    M = pd.read_csv(f"{DM}/Model.csv").set_index("ModelID")
    idx = E.index.intersection(X.index)
    d = pd.DataFrame(index=idx)
    for g in he:
        d["eff_" + g] = E[g].reindex(idx)
    for g in hx:
        d["exp_" + g] = X[g].reindex(idx)
    d["lin"] = M["OncotreeLineage"].astype(str).reindex(idx)
    d["dis"] = M["OncotreePrimaryDisease"].astype(str).reindex(idx)
    d["lymph"] = d.lin.str.contains("Lymph", case=False)
    d["B"] = d.lymph & ~d.dis.str.contains("T-Cell|T Cell|NK|T-Lympho",
                                           case=False, regex=True)
    d = d.dropna(subset=["eff_SLC7A11", "exp_AIFM2"])
    d.to_csv(f"{RES}/216_lines.csv")
    nl = d[~d.lymph]
    log("可用 %d 株（B 細胞腫瘤 %d、非淋巴系 %d）"
        % (len(d), int(d.B.sum()), len(nl)))
    out = {"n_lines": int(len(d)), "n_B": int(d.B.sum()),
           "n_non_lymphoid": int(len(nl))}

    u = stats.mannwhitneyu(d[d.B].exp_AIFM2, nl.exp_AIFM2, alternative="less")
    out["P1_AIFM2_low_in_B"] = {
        "median_B": float(d[d.B].exp_AIFM2.median()),
        "median_non_lymphoid": float(nl.exp_AIFM2.median()),
        "p": float(u.pvalue),
        "pass": bool(u.pvalue < 0.05
                     and d[d.B].exp_AIFM2.median() < nl.exp_AIFM2.median())}
    log("\nP1 `AIFM2` 表現：B 細胞中位 %.3f vs 非淋巴系 %.3f，單尾 p = %.2e → %s"
        % (d[d.B].exp_AIFM2.median(), nl.exp_AIFM2.median(), u.pvalue,
           "過" if out["P1_AIFM2_low_in_B"]["pass"] else "不過"))

    r_all, p_all = stats.spearmanr(d.exp_AIFM2, d.eff_SLC7A11)
    r_nl, p_nl = stats.spearmanr(nl.exp_AIFM2, nl.eff_SLC7A11)
    out["P2_within_non_lymphoid"] = {
        "all": {"r": float(r_all), "p": float(p_all), "n": int(len(d))},
        "non_lymphoid": {"r": float(r_nl), "p": float(p_nl), "n": int(len(nl))},
        "expected_sign": "正（AIFM2 高 → gene effect 較不負 → 較不依賴）",
        "pass": bool(r_nl > 0 and p_nl < 0.05)}
    log("P2 `AIFM2` 表現 vs `SLC7A11` 依賴：全庫 r = %+.3f (p=%.1e)｜"
        "**非淋巴系 r = %+.3f (p=%.1e)** → %s"
        % (r_all, p_all, r_nl, p_nl,
           "過" if out["P2_within_non_lymphoid"]["pass"] else "不過"))

    log("P3 建虛無（%d 個隨機基因的表現對 `SLC7A11` 依賴性，非淋巴系內）…" % N_NULL)
    pool = [c for c in xhdr.columns[1:]]
    pick = RNG.permutation(len(pool))[:N_NULL]
    nullr, step = [], 250
    cols = [pool[i] for i in pick]
    for i in range(0, len(cols), step):
        sub = pd.read_csv(f"{DM}/OmicsExpressionProteinCodingGenesTPMLogp1.csv",
                          usecols=[xf] + cols[i:i + step]).set_index(xf)
        sub = sub.reindex(nl.index)
        for c in sub.columns:
            v = sub[c]
            ok = v.notna()
            if ok.sum() > 300:
                rr = stats.spearmanr(v[ok], nl.eff_SLC7A11[ok]).statistic
                if rr == rr:
                    nullr.append(float(rr))
    nullr = np.array(nullr)
    emp = float((nullr >= r_nl).mean())
    out["P3_empirical_null"] = {
        "n_random_genes": int(len(nullr)), "null_mean": float(nullr.mean()),
        "null_sd": float(nullr.std()),
        "null_q95": float(np.quantile(nullr, 0.95)),
        "observed": float(r_nl), "empirical_p": emp, "pass": bool(emp < 0.05)}
    log("   虛無 r = %.3f ± %.3f（95 分位 %.3f）；觀察 %+.3f → 經驗 p = %.4f → %s"
        % (nullr.mean(), nullr.std(), np.quantile(nullr, 0.95), r_nl, emp,
           "過" if out["P3_empirical_null"]["pass"] else "不過"))

    r4 = stats.spearmanr(nl.exp_AIFM2, nl.eff_RPL13A, nan_policy="omit")
    out["P4_negative_control"] = {"r": float(r4.statistic), "p": float(r4.pvalue),
                                  "pass": bool(abs(r4.statistic) < abs(r_nl))}
    log("P4 陰性對照 `RPL13A`：r = %+.3f (p=%.2e)，標的 %+.3f → %s"
        % (r4.statistic, r4.pvalue, r_nl,
           "過" if out["P4_negative_control"]["pass"] else "不過"))

    # P5：2×2（AIFM2 低/高 × GPX4 依賴/不依賴），非淋巴系內
    q = nl.exp_AIFM2.quantile(0.25)
    cells = {}
    for alab, amask in (("AIFM2低", nl.exp_AIFM2 <= q), ("AIFM2高", nl.exp_AIFM2 > q)):
        for glab, gmask in (("GPX4依賴", nl.eff_GPX4 < -0.5),
                            ("GPX4不依賴", nl.eff_GPX4 >= -0.5)):
            sub = nl[amask & gmask]
            if len(sub) == 0:
                continue
            cells["%s×%s" % (alab, glab)] = {
                "n": int(len(sub)),
                "frac_SLC7A11_dependent": float((sub.eff_SLC7A11 < -0.5).mean()),
                "median_SLC7A11_effect": float(sub.eff_SLC7A11.median())}
    out["P5_two_by_two_non_lymphoid"] = cells
    log("\nP5 非淋巴系內的 2×2（`AIFM2` 下四分位、`GPX4` 依賴門檻 −0.5）：")
    for k, vv in cells.items():
        log("   %-22s n=%4d｜`SLC7A11` 依賴 %.1f%%｜effect 中位 %+.3f"
            % (k, vv["n"], 100 * vv["frac_SLC7A11_dependent"],
               vv["median_SLC7A11_effect"]))
    lo = cells.get("AIFM2低×GPX4依賴", {}).get("median_SLC7A11_effect")
    hi = cells.get("AIFM2高×GPX4不依賴", {}).get("median_SLC7A11_effect")
    out["P5_pass"] = {"low_backup_median": lo, "high_backup_median": hi,
                      "pass": bool(lo is not None and hi is not None and lo < hi)}
    log("   → 無備援（AIFM2低×GPX4依賴）%.3f vs 有備援 %.3f → %s"
        % (lo if lo is not None else float("nan"),
           hi if hi is not None else float("nan"),
           "過" if out["P5_pass"]["pass"] else "不過"))

    passed = [k for k, v in out.items() if isinstance(v, dict) and v.get("pass") is True]
    failed = [k for k, v in out.items() if isinstance(v, dict) and v.get("pass") is False]
    out["summary"] = {"passed": passed, "failed": failed}
    log("\n通過：%s｜未通過：%s"
        % (", ".join(passed) or "（無）", ", ".join(failed) or "（無）"))
    json.dump(out, open(f"{RES}/216_no_backup_arm.json", "w"), indent=2,
              ensure_ascii=False)
    log("DONE → results/216_*")


if __name__ == "__main__":
    main()
