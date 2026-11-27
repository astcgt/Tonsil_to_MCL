#!/usr/bin/env python
"""步驟 201（B3，修正版）：只控制「半胱胺酸從哪裡來」，而不是重建整個培養基。

第一版錯在哪
------------
步驟 200 只調了胱胺酸交換反應的上限，但 Human-GEM 的預設培養基是全開的
（最大生長速率 138.5，等於什麼都能從外面拿）。細胞可以直接攝取半胱胺酸、
穀胱甘肽等等，所以把胱胺酸關到 0 也毫無影響，全部敲除效應都是 0。

正確的設計
----------
真正要問的不是「RPMI 長什麼樣」，而是：
**當半胱胺酸只能從胱胺酸來的時候，xCT 會不會變成必需？**

所以只做一件最小的介入：關掉所有**直接提供半胱胺酸的旁路**
（半胱胺酸本身、穀胱甘肽、高半胱胺酸/高胱胺酸、半胱胺醯甘胺酸…），
保留甲硫胺酸（硫轉移途徑 CBS→CTH 仍可用，這是生物學上正確的替代路線，
也正是步驟 134 用表現量當代理去測的那條）。

然後做 2 × 2：胱胺酸 {充裕, 匱乏} × 甲硫胺酸 {充裕, 匱乏}
每一格做單基因敲除，看 SLC7A11 何時才變成瓶頸。

這個設計會自己回答步驟 134 留下的問題：
  若「胱胺酸匱乏」單獨不夠、要「胱胺酸＋甲硫胺酸都匱乏」才讓 xCT 必需，
  那就說明硫轉移途徑是真正的緩衝，而 DepMap 測不到是因為培養基兩者都充裕。

⚠️ 限制：FBA 是穩態化學計量，沒有動力學、氧化壓力或鐵死亡。
   它回答「這條路在化學計量上是不是瓶頸」，不回答「細胞會不會死」。

輸出 results/201_*
"""
import os
import json
import warnings

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "4"
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
RES = f"{BASE}/results"
GEM = f"{BASE}/data/aux/gem/Human-GEM.xml"

KO_GENES = ["SLC7A11", "SLC3A2", "GCLC", "GCLM", "GSS", "G6PD", "PGD",
            "SLC2A1", "GPX4", "CTH", "CBS", "MAT2A", "FASN"]
# 直接或間接提供半胱胺酸的旁路，要關掉（胱胺酸本身除外，它是我們的自變數）
BYPASS = ["cysteine", "glutathione", "homocysteine", "homocystine",
          "cysteinylglycine", "cystathionine", "taurine", "hypotaurine",
          "gamma-glutamylcysteine", "n-acetyl", "thiosulfate", "sulfide"]
RICH, POOR = 1.0, 0.01     # mmol/gDW/h，相對的充裕 / 匱乏


def log(*a):
    print("[201]", *a, flush=True)


def ex_by_name(model, needle):
    out = []
    for r in model.exchanges:
        nm = (r.name or "").lower()
        mn = " ".join((m.name or "").lower() for m in r.metabolites)
        if needle in nm or needle in mn:
            out.append(r)
    return out


def main():
    from cobra.io import read_sbml_model
    log("載入 Human-GEM …")
    m = read_sbml_model(GEM)
    log("模型：%d 反應、%d 基因" % (len(m.reactions), len(m.genes)))
    wt0 = m.slim_optimize()
    log("未介入的最大生長速率 = %.3f" % wt0)

    cys_ex = [r for r in m.exchanges if (r.name or "").lower() == "exchange of cystine"]
    if not cys_ex:
        cys_ex = ex_by_name(m, "cystine")
    CYS = cys_ex[0]
    met_ex = [r for r in m.exchanges
              if (r.name or "").lower() in ("exchange of methionine",
                                            "exchange of l-methionine")]
    if not met_ex:
        met_ex = [r for r in ex_by_name(m, "methionine")
                  if "homo" not in (r.name or "").lower()]
    MET = met_ex[0]
    log("胱胺酸交換: %s (%s)" % (CYS.id, CYS.name))
    log("甲硫胺酸交換: %s (%s)" % (MET.id, MET.name))

    closed = []
    for needle in BYPASS:
        for r in ex_by_name(m, needle):
            if r.id in (CYS.id, MET.id):
                continue
            if r.lower_bound < 0:
                r.lower_bound = 0.0
                closed.append((r.id, r.name))
    log("關掉 %d 個直接提供半胱胺酸的旁路交換反應" % len(closed))
    for rid, nm in closed[:10]:
        log("    %s  %s" % (rid, nm))
    wt1 = m.slim_optimize()
    log("關掉旁路後的生長速率 = %.3f" % wt1)

    sym2rxn = {}
    for g in KO_GENES:
        hits = [gg for gg in m.genes if (gg.name or "").upper() == g
                or gg.id.upper() == g]
        sym2rxn[g] = sorted({r.id for gg in hits for r in gg.reactions})

    rows = []
    for cys_ub in (RICH, POOR):
        for met_ub in (RICH, POOR):
            with m:
                CYS.lower_bound = -cys_ub
                MET.lower_bound = -met_ub
                s0 = m.slim_optimize()
                s0 = 0.0 if (s0 is None or not np.isfinite(s0)) else s0
                cond = ("胱胺酸%s / 甲硫胺酸%s"
                        % ("充裕" if cys_ub == RICH else "匱乏",
                           "充裕" if met_ub == RICH else "匱乏"))
                log("\n--- %s → 野生型生長 %.4f ---" % (cond, s0))
                for g, rx in sym2rxn.items():
                    if not rx:
                        continue
                    with m:
                        for rid in rx:
                            m.reactions.get_by_id(rid).bounds = (0, 0)
                        s1 = m.slim_optimize()
                        s1 = 0.0 if (s1 is None or not np.isfinite(s1)) else s1
                    eff = (1 - s1 / s0) if s0 > 1e-9 else np.nan
                    rows.append(dict(cystine=cys_ub, methionine=met_ub,
                                     condition=cond, gene=g, growth_wt=s0,
                                     growth_ko=s1, ko_effect=eff))
                    if np.isfinite(eff) and eff > 0.001:
                        log("    %-8s 敲除效應 %.3f" % (g, eff))

    D = pd.DataFrame(rows)
    D.to_csv(f"{RES}/201_cystine_fba.csv", index=False)
    piv = D.pivot(index="gene", columns="condition", values="ko_effect")
    pd.set_option("display.width", 240)
    log("\n=== 敲除效應（1 − 相對生長），四種營養條件 ===")
    log("\n" + piv.round(4).to_string())

    out = {"model": "Human-GEM", "cystine_exchange": CYS.id,
           "methionine_exchange": MET.id,
           "n_bypass_closed": len(closed),
           "bypass_closed": [{"id": a, "name": b} for a, b in closed],
           "growth_unconstrained": float(wt0),
           "ko_effect": json.loads(piv.round(6).to_json()),
           "caveat": "FBA 為穩態化學計量；無動力學、氧化壓力與鐵死亡"}
    if "SLC7A11" in piv.index:
        v = piv.loc["SLC7A11"]
        log("\n★ SLC7A11 在四種條件下的敲除效應：")
        for c in piv.columns:
            log("   %-26s %.4f" % (c, v[c]))
        out["slc7a11"] = {str(c): float(v[c]) for c in piv.columns}
    json.dump(out, open(f"{RES}/201_cystine_fba.json", "w"), indent=2,
              ensure_ascii=False)
    log("\nDONE → results/201_*")


if __name__ == "__main__":
    main()
