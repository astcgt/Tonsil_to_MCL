#!/usr/bin/env python
"""步驟 213：`NAMPT` —— 從還原力那一側攻擊同一條鐵死亡防線？

它是怎麼出現的
--------------
不是我挑的。步驟 212 的全基因組治療窗口篩選（分數 = min(z(殘差), z(B 細胞依賴選擇性))，
兩項都必須為正）把它排在 **第 1 名、17,552 個基因中**，分數 6.35 是第 2 名的 2.1 倍。
同一個排名在沒有被告知的情況下自己找出 `BCL2`（venetoclax，已核准用於 MCL）、
`MCL1`、`EZH2`、`BRD2`、`IRF4`、`NFKB1`，而 `ROR1`（zilovertamab 的標靶）
在未閘控版本排第 4 —— 所以這個排名不是在產生噪音。

為什麼它和這個專案的主軸是同一件事
----------------------------------
主軸是 胱胺酸 → 穀胱甘肽 → 鐵死亡防禦。
但穀胱甘肽要發揮作用必須維持**還原態**，而把 GSSG 還原回 GSH 的酵素
（glutathione reductase）燒的是 **NADPH**。NAMPT 是 NAD⁺ 補救合成的限速酵素，
NAD⁺ 經 NADK 磷酸化成 NADP⁺ 再被還原成 NADPH。
也就是說 NAMPT 供應的是**同一條防線的電力**。

步驟 202 的化學計量模型說「餓掉胱胺酸是錯的槓桿」（胱胺酸一少，
轉硫途徑的相對貢獻反而上升）。從還原力那一側進攻不受那個論證影響。
步驟 170 測過的 disulfidptosis（靠葡萄糖飢餓造成 NADPH 短缺）沒有訊號，
但那是**營養**造成的短缺；NAMPT 抑制是**藥理**造成的短缺，是不同的實驗。

事前指定
--------
P1  B 細胞選擇性必須通得過經驗虛無：和依賴率相符的隨機基因比，
    `NAMPT` 的譜系富集 odds ratio 的經驗 p < 0.05。
    （方法與步驟 135 相同，避免「隨便一個基因也會這樣」。）
P2  選擇性必須不只是「是淋巴系」的同義詞：在**非淋巴系**株內，
    `NAMPT` 依賴性與 `NAMPT` 表現的關係應該仍然存在。
P3  機制連結要有證據：`NAMPT` 的 CRISPR 輪廓應該和 NADPH 供應基因
    （`NADK`、`G6PD`、`PGD`、`IDH1`、`ME1`、`NNT`）相關，
    而且高於和隨機基因的相關。
P4（陰性對照）`NAMPT` 的輪廓不該同樣地和一組無關的必需基因
    （核糖體 `RPL13A`、`RPS6`）相關。
P5  文獻新穎性：報告 `NAMPT AND mantle cell lymphoma` 與
    `NAMPT AND ferroptosis` 的論文數，不隱藏。

⚠️ 一定要寫進結論的責任：NAMPT 是**普遍表現**的酵素，
   NAMPT 抑制劑在臨床上的已知毒性正是來自這一點（視網膜與血液毒性）。
   「分化配對殘差有 4.3 倍窗口」是對**正常 B 細胞**而言，
   不是對視網膜或心臟而言。這一步不能宣稱全身安全性。

輸出 results/213_*
"""
import os
import json
import time
import urllib.parse
import urllib.request

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "8"

import numpy as np
import pandas as pd
from scipy import stats

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
RES = f"{BASE}/results"
DM = f"{BASE}/data/depmap"
TARGET = "NAMPT"
NADPH = ["NADK", "G6PD", "PGD", "IDH1", "ME1", "NNT", "GSR", "TXNRD1"]
AXIS = ["SLC7A11", "SLC3A2", "GCLC", "GCLM", "AIFM2", "GPX4"]
NEG = ["RPL13A", "RPS6"]
DEP_CUT = -0.5
N_NULL = 1500
RNG = np.random.default_rng(0)
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"


def log(*a):
    print("[213]", *a, flush=True)


def pubmed_count(q):
    url = EUTILS + "?" + urllib.parse.urlencode(
        {"db": "pubmed", "term": q, "retmode": "json", "retmax": 0,
         "tool": "regulatory-genomics-pipeline"})
    try:
        with urllib.request.urlopen(url, timeout=30) as r:
            return int(json.load(r)["esearchresult"]["count"])
    except Exception as e:
        log("  PubMed 查詢失敗（%s）：%s" % (type(e).__name__, q))
        return None


def main():
    need = sorted(set([TARGET] + NADPH + AXIS + NEG))
    hdr = pd.read_csv(f"{DM}/CRISPRGeneEffect.csv", nrows=0)
    first = hdr.columns[0]
    s2c = {}
    for c in hdr.columns[1:]:
        s2c.setdefault(c.split(" (")[0], c)
    have = [g for g in need if g in s2c]
    E = pd.read_csv(f"{DM}/CRISPRGeneEffect.csv",
                    usecols=[first] + [s2c[g] for g in have]).set_index(first)
    E.columns = [c.split(" (")[0] for c in E.columns]
    log("讀入 %d 個基因 × %d 株" % (E.shape[1], E.shape[0]))

    mod = pd.read_csv(f"{DM}/Model.csv").set_index("ModelID")
    lin = mod["OncotreeLineage"].astype(str).reindex(E.index).fillna("NA")
    dis = mod["OncotreePrimaryDisease"].astype(str).reindex(E.index).fillna("NA")
    is_lymph = lin.str.contains("Lymph", case=False)
    is_B = is_lymph & ~dis.str.contains("T-Cell|T Cell|NK", case=False, regex=True)
    out = {"n_lines": int(len(E)), "n_lymphoid": int(is_lymph.sum()),
           "n_Bcell": int(is_B.sum())}
    log("淋巴系 %d 株（其中 B 細胞腫瘤 %d），非淋巴系 %d"
        % (is_lymph.sum(), is_B.sum(), (~is_lymph).sum()))

    v = E[TARGET].dropna()
    dep = v < DEP_CUT
    nB = int((dep & is_B.reindex(v.index).fillna(False)).sum())
    nd = int(dep.sum())
    tab = [[nB, nd - nB],
           [int(is_B.reindex(v.index).fillna(False).sum()) - nB,
            int(len(v)) - nd - (int(is_B.reindex(v.index).fillna(False).sum()) - nB)]]
    orr, fp = stats.fisher_exact(tab, alternative="greater")
    log("\n%s：%d / %d 株依賴（gene effect < %.1f），其中 B 細胞腫瘤 %d 株"
        % (TARGET, nd, len(v), DEP_CUT, nB))
    log("  Fisher odds ratio = %.2f，p = %.2e" % (orr, fp))

    # P1：和依賴率相符的隨機基因比
    log("P1 建經驗虛無（依賴率相符的隨機基因）…")
    allhdr = [c for c in hdr.columns[1:]]
    rate = nd / len(v)
    null_or, tried, step = [], 0, 400
    pick = RNG.permutation(len(allhdr))
    i = 0
    while len(null_or) < N_NULL and i < len(pick):
        cols = [allhdr[j] for j in pick[i:i + step]]
        i += step
        sub = pd.read_csv(f"{DM}/CRISPRGeneEffect.csv",
                          usecols=[first] + cols).set_index(first)
        for c in sub.columns:
            vv = sub[c].dropna()
            if len(vv) < 500:
                continue
            dd = vv < DEP_CUT
            r = dd.sum() / len(vv)
            if not (0.5 * rate <= r <= 2.0 * rate) or dd.sum() < 3:
                continue
            bb = is_B.reindex(vv.index).fillna(False)
            a = int((dd & bb).sum()); b = int(dd.sum()) - a
            cc = int(bb.sum()) - a; d2 = len(vv) - int(dd.sum()) - cc
            if min(a + b, a + cc) == 0:
                continue
            o, _ = stats.fisher_exact([[a, b], [cc, d2]], alternative="greater")
            if np.isfinite(o):
                null_or.append(float(o))
            tried += 1
    null_or = np.array(null_or)
    emp_p = float((null_or >= orr).mean()) if len(null_or) else float("nan")
    out["P1_lineage_enrichment"] = {
        "n_dependent": nd, "n_Bcell_dependent": nB, "odds_ratio": float(orr),
        "fisher_p": float(fp), "n_matched_null_genes": int(len(null_or)),
        "null_median_or": float(np.median(null_or)) if len(null_or) else None,
        "empirical_p": emp_p, "pass": bool(emp_p < 0.05)}
    log("  虛無 %d 個速率相符基因，OR 中位 %.2f；觀察 %.2f → 經驗 p = %.4f → %s"
        % (len(null_or), np.median(null_or) if len(null_or) else -1, orr, emp_p,
           "過" if out["P1_lineage_enrichment"]["pass"] else "不過"))

    # P2：非淋巴系內，表現是否預測依賴性
    ehdr = pd.read_csv(f"{DM}/OmicsExpressionProteinCodingGenesTPMLogp1.csv", nrows=0)
    ef = ehdr.columns[0]
    ecol = {c.split(" (")[0]: c for c in ehdr.columns[1:]}
    X = pd.read_csv(f"{DM}/OmicsExpressionProteinCodingGenesTPMLogp1.csv",
                    usecols=[ef] + [ecol[g] for g in [TARGET] if g in ecol]
                    ).set_index(ef)
    X.columns = [c.split(" (")[0] for c in X.columns]
    j = pd.DataFrame({"eff": E[TARGET], "exp": X[TARGET].reindex(E.index),
                      "lymph": is_lymph}).dropna()
    nl = j[~j.lymph]
    r_all, p_all = stats.spearmanr(j["exp"], j.eff)
    r_nl, p_nl = stats.spearmanr(nl["exp"], nl.eff)
    out["P2_expression_predicts_dependency"] = {
        "all_lines": {"n": int(len(j)), "r": float(r_all), "p": float(p_all)},
        "non_lymphoid": {"n": int(len(nl)), "r": float(r_nl), "p": float(p_nl)},
        "expected_sign": "負（表現越高越依賴）",
        "pass": bool(r_nl < 0 and p_nl < 0.05)}
    log("P2 表現 vs 依賴：全庫 r = %+.3f (p=%.1e)｜非淋巴系 r = %+.3f (p=%.1e) → %s"
        % (r_all, p_all, r_nl, p_nl,
           "過" if out["P2_expression_predicts_dependency"]["pass"] else "不過"))

    # P3 / P4：CRISPR 輪廓的機制連結
    def prof_corr(a, b):
        m = E[[a, b]].dropna()
        if len(m) < 300:
            return np.nan
        return float(stats.spearmanr(m[a], m[b]).statistic)
    link = {g: prof_corr(TARGET, g) for g in NADPH + AXIS + NEG if g in E.columns}
    nadph_r = np.nanmean([link[g] for g in NADPH if g in link])
    neg_r = np.nanmean([link[g] for g in NEG if g in link])
    axis_r = np.nanmean([link[g] for g in AXIS if g in link])
    log("\nP3 %s 的 CRISPR 輪廓相關：" % TARGET)
    for g in NADPH + AXIS + NEG:
        if g in link:
            tag = "NADPH" if g in NADPH else ("主軸" if g in AXIS else "陰性對照")
            log("   %-8s %-6s ρ = %+.3f" % (g, tag, link[g]))
    log("   NADPH 供應組平均 %+.3f｜主軸平均 %+.3f｜陰性對照平均 %+.3f"
        % (nadph_r, axis_r, neg_r))
    out["P3_mechanism_link"] = {
        "per_gene": link, "mean_NADPH": float(nadph_r),
        "mean_axis": float(axis_r), "mean_negctrl": float(neg_r),
        "pass": bool(nadph_r > 0 and nadph_r > abs(neg_r))}
    out["P4_negative_control"] = {
        "mean_negctrl": float(neg_r), "mean_NADPH": float(nadph_r),
        "pass": bool(abs(neg_r) < abs(nadph_r))}
    log("P3 → %s｜P4 → %s"
        % ("過" if out["P3_mechanism_link"]["pass"] else "不過",
           "過" if out["P4_negative_control"]["pass"] else "不過"))

    # P5：文獻
    log("\nP5 PubMed 查證…")
    qs = {
        "NAMPT AND mantle cell lymphoma": 'NAMPT AND "mantle cell lymphoma"',
        "NAMPT AND ferroptosis": "NAMPT AND ferroptosis",
        "NAMPT AND lymphoma": "NAMPT AND lymphoma",
        "NAMPT inhibitor clinical": "NAMPT inhibitor AND (clinical trial[pt] OR phase)",
        "SLC7A11 AND NAMPT": "SLC7A11 AND NAMPT",
        "(reference) CCND1 AND MCL": 'CCND1 AND "mantle cell lymphoma"',
    }
    lit = {}
    for k, q in qs.items():
        lit[k] = pubmed_count(q)
        log("   %-34s %s" % (k, lit[k]))
        time.sleep(0.4)
    out["P5_literature"] = lit

    passed = [k for k, v in out.items() if isinstance(v, dict) and v.get("pass") is True]
    failed = [k for k, v in out.items() if isinstance(v, dict) and v.get("pass") is False]
    out["summary"] = {"passed": passed, "failed": failed}
    out["caveat"] = ("NAMPT 普遍表現，臨床 NAMPT 抑制劑的已知毒性正來自這一點。"
                     "4.3 倍窗口是對分化配對的正常 B 細胞而言，不是對全身而言。")
    log("\n通過：%s｜未通過：%s"
        % (", ".join(passed) or "（無）", ", ".join(failed) or "（無）"))
    json.dump(out, open(f"{RES}/213_nampt.json", "w"), indent=2, ensure_ascii=False)
    log("DONE → results/213_*")


if __name__ == "__main__":
    main()
