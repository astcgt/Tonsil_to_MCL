#!/usr/bin/env python
"""步驟 200（B3）：字面意義上「把培養基的胱胺酸調低」。

為什麼
------
步驟 121/134 顯示 SLC7A11 在 1,178 株的中位 gene effect 是 +0.086（敲掉反而略好）。
最常見的解釋是「培養基的胱胺酸（RPMI 約 200 µM）遠比體內（40–60 µM）充裕，
所以這條保險絲沒被點火」。步驟 134 用硫轉移能力當代理去檢驗，結果是中性。
那個代理很弱（mRNA ≠ 通量）。

基因組尺度代謝模型可以**直接**做這件事：把胱胺酸交換反應的攝取上限
從培養基濃度調到體內濃度，看細胞的最大生長速率與必需反應集合怎麼變。
這是唯一能真正「換培養基」的計算方法。

做什麼
------
1. 載入 Human-GEM（SysBio Chalmers），設成 RPMI-1640 類的培養基。
2. 胱胺酸攝取上限從 RPMI 水準掃到體內水準再到 0，記錄最大生長速率。
3. 在每個濃度下，對 SLC7A11（xCT 反應）、GCLC、GCLM、GSS、G6PD、SLC2A1
   做單反應敲除，看生長下降多少 → 得到「依賴性 vs 胱胺酸供應」的曲線。
4. **關鍵讀出**：SLC7A11 的敲除效應是否隨胱胺酸下降而放大？
   若是，就從代謝化學層面支持「DepMap 測不到是因為培養基太營養」。

⚠️ 限制（必須與結果一起報）
   FBA 假設穩態與最大化生長，沒有動力學、沒有氧化壓力、沒有鐵死亡。
   它能說的是「化學計量上這條路是不是瓶頸」，不能說「細胞會不會死」。

輸出 results/200_*
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

# 要測的基因（用基因符號找對應反應）
KO_GENES = ["SLC7A11", "SLC3A2", "GCLC", "GCLM", "GSS", "G6PD", "PGD",
            "SLC2A1", "GPX4", "CTH", "CBS", "FASN"]
# 胱胺酸攝取上限掃描（mmol/gDW/h）。RPMI 的胱胺酸約 200 µM，血漿約 40–60 µM。
SCAN = [1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.0]


def log(*a):
    print("[200]", *a, flush=True)


def main():
    import cobra
    from cobra.io import read_sbml_model
    log("載入 Human-GEM（約 50 MB，需要一兩分鐘）…")
    m = read_sbml_model(GEM)
    log("模型：%d 反應、%d 代謝物、%d 基因"
        % (len(m.reactions), len(m.metabolites), len(m.genes)))

    # 找胱胺酸的交換反應
    cys_ex = [r for r in m.exchanges
              if "cystine" in (r.name or "").lower()
              or any("cystine" in (met.name or "").lower() for met in r.metabolites)]
    log("胱胺酸交換反應候選: %s" % [(r.id, r.name) for r in cys_ex][:4])
    if not cys_ex:
        log("找不到胱胺酸交換反應，改找 cysteine"); 
        cys_ex = [r for r in m.exchanges
                  if "cysteine" in (r.name or "").lower()]
        log("  cysteine 候選: %s" % [(r.id, r.name) for r in cys_ex][:4])
    if not cys_ex:
        log("都找不到，結束"); return
    EX = cys_ex[0]
    log("使用交換反應: %s (%s)" % (EX.id, EX.name))

    bio = m.objective.expression
    base = m.optimize()
    log("預設培養基下最大生長速率 = %.4f (status %s)" % (base.objective_value, base.status))

    # gene → reactions
    sym2rxn = {}
    for g in KO_GENES:
        hits = [gg for gg in m.genes
                if (gg.name or "").upper() == g or gg.id.upper() == g]
        rx = sorted({r.id for gg in hits for r in gg.reactions})
        sym2rxn[g] = rx
        log("  %-8s → %d 個反應 %s" % (g, len(rx), rx[:3]))

    rows = []
    for ub in SCAN:
        with m:
            EX.lower_bound = -ub            # 攝取為負值
            s0 = m.slim_optimize()
            if s0 is None or not np.isfinite(s0):
                s0 = 0.0
            for g, rx in sym2rxn.items():
                if not rx:
                    continue
                with m:
                    for rid in rx:
                        m.reactions.get_by_id(rid).bounds = (0, 0)
                    s1 = m.slim_optimize()
                    s1 = 0.0 if (s1 is None or not np.isfinite(s1)) else s1
                rows.append(dict(cystine_uptake_ub=ub, gene=g,
                                 growth_wt=s0, growth_ko=s1,
                                 rel_growth=(s1 / s0) if s0 > 1e-9 else np.nan,
                                 ko_effect=(1 - s1 / s0) if s0 > 1e-9 else np.nan))
        log("胱胺酸上限 %.4f → 野生型生長 %.4f" % (ub, s0))

    D = pd.DataFrame(rows)
    D.to_csv(f"{RES}/200_cystine_fba.csv", index=False)
    piv = D.pivot(index="gene", columns="cystine_uptake_ub", values="ko_effect")
    piv = piv[sorted(piv.columns, reverse=True)]
    pd.set_option("display.width", 240)
    log("\n=== 敲除效應（1 − 相對生長）隨胱胺酸供應下降的變化 ===")
    log("（欄＝胱胺酸攝取上限，由充裕到零；值越大＝越依賴）")
    log("\n" + piv.round(3).to_string())

    out = {"model": "Human-GEM", "exchange": EX.id,
           "scan": SCAN,
           "ko_effect": json.loads(piv.round(5).to_json()),
           "gene_to_reactions": sym2rxn,
           "caveat": ("FBA 假設穩態與最大化生長，沒有動力學、氧化壓力與鐵死亡；"
                      "它回答的是化學計量上的瓶頸，不是細胞會不會死")}
    # 關鍵判讀：SLC7A11 的敲除效應有沒有隨胱胺酸下降而放大
    if "SLC7A11" in piv.index:
        v = piv.loc["SLC7A11"].dropna()
        if len(v) >= 3:
            rich, poor = float(v.iloc[0]), float(v.iloc[-1])
            out["slc7a11"] = {"ko_effect_rich_medium": rich,
                              "ko_effect_poor_medium": poor,
                              "amplification": poor - rich}
            log("\n★ SLC7A11 敲除效應：充裕 %.3f → 匱乏 %.3f（放大 %+.3f）"
                % (rich, poor, poor - rich))
            log("   若放大量明顯為正，就從化學計量層面支持『培養基太營養所以測不到』")
    json.dump(out, open(f"{RES}/200_cystine_fba.json", "w"), indent=2, ensure_ascii=False)
    log("\nDONE → results/200_*")


if __name__ == "__main__":
    main()
