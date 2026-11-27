#!/usr/bin/env python
"""步驟 217：「表現越高越不依賴」是發現，還是 DepMap 的背景性質？

這一題是三個互相獨立的「失敗」逼出來的
--------------------------------------
今天有三處事前指定的檢定，都因為**方向相反**而不過：
  步驟 213 P2  `NAMPT` 表現越高 → 越不依賴（非淋巴系 1,012 株 r = +0.148）
  步驟 214 P2  `NAPRT` 表現 vs `NAMPT` 依賴，方向相反
  步驟 205 P8（子代理）`SLC7A11` 自己 r_exp = +0.357，而且共依賴訊號住在
               `SLC7A11` **低**表現的株裡 —— 和病人腫瘤（高於配對正常 200 倍）相反

三個不同的基因、三個不同的分析，撞出同一個方向。這通常不是三個發現，
而是一個背景性質。如果「自己的表現越高 → 自己越不被需要」在 DepMap 裡
**普遍成立**，那上面三處就都不是訊號，我今天的解讀必須全部降級。

事前指定
--------
H   隨機抽基因，算 corr(自己的表現, 自己的 gene effect) 的分佈。
    若中位數顯著大於 0，則「方向相反」是背景，不是發現。
P1  報告分佈的中位數、四分位與正值比例。
P2  `SLC7A11`、`NAMPT`、`NAPRT` 落在這個分佈的哪個百分位。
P3（陰性對照）把表現與依賴的株標籤打亂，中位數應該回到 0 ——
    確認流程本身不會製造出偏移。

⚠️ 這一步不產生新的生物學，它只決定今天三個結論該怎麼寫。

輸出 results/217_*
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
N_GENES = 600
FOCUS = ["SLC7A11", "NAMPT", "NAPRT", "GPX4", "AIFM2", "SCLY", "GCLC", "GCLM"]
RNG = np.random.default_rng(0)


def log(*a):
    print("[217]", *a, flush=True)


def main():
    eh = pd.read_csv(f"{DM}/CRISPRGeneEffect.csv", nrows=0)
    xh = pd.read_csv(f"{DM}/OmicsExpressionProteinCodingGenesTPMLogp1.csv", nrows=0)
    ef, xf = eh.columns[0], xh.columns[0]
    emap, xmap = {}, {}
    for c in eh.columns[1:]:
        emap.setdefault(c.split(" (")[0], c)
    for c in xh.columns[1:]:
        xmap.setdefault(c.split(" (")[0], c)
    shared = sorted(set(emap) & set(xmap))
    log("表現與依賴都有的基因 %d 個" % len(shared))
    pick = list(RNG.choice(shared, size=min(N_GENES, len(shared)), replace=False))
    genes = sorted(set(pick) | set(g for g in FOCUS if g in emap and g in xmap))
    log("抽樣 %d 個（含 %d 個關注基因）" % (len(genes), len(FOCUS)))

    E = pd.read_csv(f"{DM}/CRISPRGeneEffect.csv",
                    usecols=[ef] + [emap[g] for g in genes]).set_index(ef)
    E.columns = [c.split(" (")[0] for c in E.columns]
    X = pd.read_csv(f"{DM}/OmicsExpressionProteinCodingGenesTPMLogp1.csv",
                    usecols=[xf] + [xmap[g] for g in genes]).set_index(xf)
    X.columns = [c.split(" (")[0] for c in X.columns]
    idx = E.index.intersection(X.index)
    E, X = E.loc[idx], X.loc[idx]
    log("共有 %d 株" % len(idx))

    rows = []
    for g in genes:
        a, b = X[g], E[g]
        ok = a.notna() & b.notna()
        if ok.sum() < 300 or a[ok].std() < 1e-6:
            continue
        rows.append((g, float(stats.spearmanr(a[ok], b[ok]).statistic),
                     int(ok.sum())))
    D = pd.DataFrame(rows, columns=["gene", "r_self", "n"]).dropna()
    D.to_csv(f"{RES}/217_selfcorr.csv", index=False)
    rand = D[~D.gene.isin(FOCUS)].r_self.values
    log("\n=== 隨機基因的 corr(自己的表現, 自己的 gene effect) ===")
    log("  n = %d 個基因｜中位數 %+.4f｜平均 %+.4f｜四分位 %+.3f ~ %+.3f"
        % (len(rand), np.median(rand), rand.mean(),
           np.quantile(rand, .25), np.quantile(rand, .75)))
    log("  正值比例 %.1f%%" % (100 * (rand > 0).mean()))
    w = stats.wilcoxon(rand)
    log("  Wilcoxon 檢定中位數 ≠ 0：p = %.3e" % w.pvalue)

    out = {"n_random_genes": int(len(rand)),
           "median": float(np.median(rand)), "mean": float(rand.mean()),
           "q25": float(np.quantile(rand, .25)),
           "q75": float(np.quantile(rand, .75)),
           "frac_positive": float((rand > 0).mean()),
           "wilcoxon_p": float(w.pvalue),
           "background_is_positive": bool(np.median(rand) > 0 and w.pvalue < 0.05)}

    log("\n=== 關注基因落在這個分佈的哪裡 ===")
    foc = {}
    for g in FOCUS:
        if g not in set(D.gene):
            continue
        r = float(D.loc[D.gene == g, "r_self"].iloc[0])
        pct = float((rand < r).mean())
        foc[g] = {"r_self": r, "percentile_vs_random": pct}
        log("  %-8s r = %+.3f → 落在隨機分佈的第 %.0f 百分位 %s"
            % (g, r, 100 * pct,
               "（不特別）" if 0.05 < pct < 0.95 else "（極端）"))
    out["focus_genes"] = foc

    # P3：打亂株標籤
    sh = []
    for g in D.gene.sample(min(200, len(D)), random_state=1):
        a, b = X[g], E[g]
        ok = a.notna() & b.notna()
        if ok.sum() < 300:
            continue
        bb = RNG.permutation(b[ok].values)
        sh.append(float(stats.spearmanr(a[ok], bb).statistic))
    sh = np.array(sh)
    out["P3_shuffled_control"] = {
        "n": int(len(sh)), "median": float(np.median(sh)),
        "pass": bool(abs(np.median(sh)) < 0.02)}
    log("\nP3 打亂株標籤後中位數 %+.4f（應接近 0）→ %s"
        % (np.median(sh), "過" if out["P3_shuffled_control"]["pass"] else "不過"))

    out["interpretation"] = (
        "背景為正" if out["background_is_positive"] else "背景不為正")
    log("\n結論：%s" % (
        "「表現越高越不依賴」是 DepMap 的**背景性質**，"
        "今天三處『方向相反』的檢定因此不是獨立的發現，解讀必須降級"
        if out["background_is_positive"] else
        "背景不為正，所以那三處的方向相反確實需要個別解釋"))
    json.dump(out, open(f"{RES}/217_selfcorr_null.json", "w"), indent=2,
              ensure_ascii=False)
    log("DONE → results/217_*")


if __name__ == "__main__":
    main()
