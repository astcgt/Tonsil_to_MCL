#!/usr/bin/env python
"""步驟 234：`MAFG` 的標的檔案 —— 我們要打的那個蛋白，它自己安全嗎？冗餘嗎？

為什麼這一支該在設計更多結合劑之前做
------------------------------------
這個專案對 `SLC7A11` 做了大量 DepMap 分析，但對**真正要打的那個蛋白 `MAFG`**
從來沒有建過標的檔案。這是標準標的評估的第一步，而它可能直接否決整條路線：
如果 `MAFG` 是廣泛必需基因或高度受限基因，「治療窗口」的論述就垮了。

而且 DepMap 可以**免費檢定這個專案設計理念的核心前提**。報告的論述是：
「AP-1／BACH1 的讀取因子彼此冗餘，但它們共用一個要求 —— 在 MARE 上 `MafG` 是必要夥伴，
所以 `MafG` 是整組蛋白收斂到的唯一一點。」
但 small Maf 家族有三個成員（`MAFG`／`MAFF`／`MAFK`），它們本身也可能彼此冗餘。
**如果 `MAFF`／`MAFK` 能替代 `MAFG`，那「收斂到唯一一點」就是錯的，只打 `MafG` 不會有用。**

事先寫死的判定（在看到任何數字之前寫成）
----------------------------------------
S1 **族群遺傳學的受限程度**（gnomAD）：報 `MAFG` 的 pLI 與 LOEUF，附對照基因。
   規則：**LOEUF < 0.35 且 pLI > 0.9** → 高度受限，治療窗口的論述必須加上安全性限縮。
S2 **是否廣泛必需**（DepMap CRISPR）：gene effect < −0.5 的細胞株比例。
   規則：**> 50%** → `MAFG` 是廣泛必需 → 治療窗口垮掉，必須在報告裡寫明。
S3 **家族冗餘（本支最重要的一關，直接檢定設計理念）**：
   在細胞株層級，`MAFF` + `MAFK` 的表現量是否預測對 `MAFG` 的依賴性較弱？
   方向：若可替代，則 `MAFF`/`MAFK` 高 → 不那麼需要 `MAFG` → gene effect 較不負
   → **表現量與 gene effect 正相關**。
   規則：Spearman **r > 0**、參數 p < 0.05、且相對 1,000 個隨機基因對的經驗 p < 0.05
   → **冗餘成立**，「`MafG` 是唯一必要夥伴」必須降級為「`MafG` 是三個可能夥伴之一」。
   陰性對照：同樣的檢定換成與 `MAFG` 無關的隨機基因，不得通過。
S4 **表現廣度**：`MAFG` 在多少個 OncotreeLineage 的中位表現 > 1 (log2 TPM+1)。
   規則：> 80% 的譜系 → 普遍表現，和 `NAMPT` 同一類風險（全身毒性）。

事先預測（會記入誠實清單，不論對錯）
------------------------------------
- S1：`MAFG` 的 LOEUF 我猜 **0.3–0.7**（中度受限，不會進最受限那一層）。
- S2：**不是**廣泛必需，我猜 < 10% 的株依賴。
- S3：**這是我最擔心的一關。我預測冗餘會成立**（r > 0 且顯著），
  因為 small Maf 三個成員的 DNA 結合域高度相似、文獻上在 NRF2 訊號裡本來就被當作可互換。
  如果成立，報告第四節的設計理念必須改寫 —— 這比任何結合劑分數都重要。
- S4：普遍表現，> 80% 的譜系。

輸出 results/234_mafg_dossier.json
"""
import json
import os
import ssl
import sys
import urllib.request

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
DM = f"{BASE}/data/depmap"
RES = f"{BASE}/results"

TARGET = "MAFG"
FAMILY = ["MAFF", "MAFK"]
# 對照基因：已上市藥物的標的（窗口存在）、已知廣泛必需、已知高度受限
COMPARATORS = ["BTK", "BCL2", "CRBN", "IKZF1", "RPL13A", "POLR2A", "SOX11", "SLC7A11"]
DEP_CUT = -0.5
S2_FRAC_CUT = 0.50
N_NULL = 1000
SEED = 0

_SYS_CA = "/etc/ssl/certs/ca-certificates.crt"
CTX = (ssl.create_default_context(cafile=_SYS_CA)
       if os.path.exists(_SYS_CA) else ssl.create_default_context())


def log(*a):
    print("[234]", *a, flush=True)


def gnomad(sym):
    q = ('{gene(gene_symbol:"%s",reference_genome:GRCh38)'
         '{gene_id symbol gnomad_constraint'
         '{pLI oe_lof oe_lof_lower oe_lof_upper mis_z lof_z}}}' % sym)
    req = urllib.request.Request(
        "https://gnomad.broadinstitute.org/api",
        data=json.dumps({"query": q}).encode(),
        headers={"Content-Type": "application/json",
                 "User-Agent": "regulatory-genomics-pipeline"})
    try:
        with urllib.request.urlopen(req, timeout=45, context=CTX) as r:
            g = json.load(r)["data"]["gene"]
        return g.get("gnomad_constraint") or {}
    except Exception as e:
        log("  gnomAD 查詢失敗 %s：%s" % (sym, e))
        return {}


def pick(hdr, genes):
    s2c = {}
    for c in hdr.columns[1:]:
        s2c.setdefault(c.split(" (")[0], c)
    return {g: s2c[g] for g in genes if g in s2c}


def main():
    out = {"target": TARGET}

    # ------------------------------------------------------------ S1 ---
    log("=== S1：gnomAD 族群遺傳學受限程度 ===")
    con = {}
    for g in [TARGET] + FAMILY + COMPARATORS:
        con[g] = gnomad(g)
    log("%-10s %8s %8s %10s %8s" % ("gene", "pLI", "LOEUF", "oe_lof", "mis_z"))
    for g, c in con.items():
        if c:
            log("%-10s %8.3f %8.3f %10.3f %8.2f"
                % (g, c.get("pLI") or np.nan, c.get("oe_lof_upper") or np.nan,
                   c.get("oe_lof") or np.nan, c.get("mis_z") or np.nan))
        else:
            log("%-10s （無資料）" % g)
    tc = con.get(TARGET, {})
    pli, loeuf = tc.get("pLI"), tc.get("oe_lof_upper")
    s1_constrained = bool(loeuf is not None and pli is not None
                          and loeuf < 0.35 and pli > 0.9)
    out["S1_constraint"] = {
        "rule": "LOEUF < 0.35 且 pLI > 0.9 → 高度受限，窗口論述需加安全性限縮",
        "values": con, "target_pLI": pli, "target_LOEUF": loeuf,
        "highly_constrained": s1_constrained}
    log("  %s pLI %.3f、LOEUF %.3f → %s"
        % (TARGET, pli or np.nan, loeuf or np.nan,
           "高度受限" if s1_constrained else "非最受限層（中度或以下）"))

    # ------------------------------------------------------- 讀 DepMap ---
    need = [TARGET] + FAMILY + COMPARATORS
    he = pd.read_csv(f"{DM}/CRISPRGeneEffect.csv", nrows=0)
    ce = pick(he, need)
    E = pd.read_csv(f"{DM}/CRISPRGeneEffect.csv",
                    usecols=[he.columns[0]] + list(ce.values())).set_index(he.columns[0])
    E.columns = [c.split(" (")[0] for c in E.columns]
    hx = pd.read_csv(f"{DM}/OmicsExpressionProteinCodingGenesTPMLogp1.csv", nrows=0)
    cx = pick(hx, need)
    X = pd.read_csv(f"{DM}/OmicsExpressionProteinCodingGenesTPMLogp1.csv",
                    usecols=[hx.columns[0]] + list(cx.values())).set_index(hx.columns[0])
    X.columns = [c.split(" (")[0] for c in X.columns]
    mod = pd.read_csv(f"{DM}/Model.csv").set_index("ModelID")
    log("\nCRISPR %d 株 × %d 基因 | 表現 %d 株" % (E.shape[0], E.shape[1], X.shape[0]))

    # ------------------------------------------------------------ S2 ---
    log("\n=== S2：是否廣泛必需 ===")
    s2 = {}
    for g in E.columns:
        v = E[g].dropna()
        s2[g] = {"n": int(len(v)), "median": float(v.median()),
                 "frac_dependent": float((v < DEP_CUT).mean())}
    log("%-10s %6s %10s %14s" % ("gene", "n", "中位效應", "依賴株比例"))
    for g in [TARGET] + FAMILY + [c for c in COMPARATORS if c in s2]:
        if g in s2:
            log("%-10s %6d %10.3f %13.1f%%"
                % (g, s2[g]["n"], s2[g]["median"], s2[g]["frac_dependent"] * 100))
    s2_broad = bool(s2[TARGET]["frac_dependent"] > S2_FRAC_CUT)
    out["S2_common_essential"] = {
        "rule": "gene effect < %.1f 的株比例 > %.0f%% → 廣泛必需，窗口垮掉"
                % (DEP_CUT, S2_FRAC_CUT * 100),
        "values": s2, "broadly_essential": s2_broad}
    log("  %s 依賴株比例 %.1f%% → %s"
        % (TARGET, s2[TARGET]["frac_dependent"] * 100,
           "廣泛必需" if s2_broad else "非廣泛必需"))

    # ------------------------------------------------------------ S3 ---
    log("\n=== S3：small Maf 家族冗餘（直接檢定設計理念）===")
    common = E.index.intersection(X.index)
    dep = E.loc[common, TARGET]
    fam_present = [g for g in FAMILY if g in X.columns]
    partner = X.loc[common, fam_present].sum(1)
    ok = dep.notna() & partner.notna()
    rs = spearmanr(partner[ok], dep[ok])
    log("  %s 表現總和 vs %s gene effect：r = %+.4f, p = %.3g, n = %d"
        % ("+".join(fam_present), TARGET, rs.statistic, rs.pvalue, int(ok.sum())))

    # 經驗虛無：隨機基因的表現量對同一個 gene effect
    rng = np.random.default_rng(SEED)
    allx = pd.read_csv(f"{DM}/OmicsExpressionProteinCodingGenesTPMLogp1.csv", nrows=0)
    pool = [c for c in allx.columns[1:]
            if c.split(" (")[0] not in set([TARGET] + FAMILY)]
    sel = list(rng.choice(pool, size=min(N_NULL, len(pool)), replace=False))
    log("  抽 %d 個隨機基因建經驗虛無…" % len(sel))
    null = []
    CH = 200
    for a in range(0, len(sel), CH):
        blk = pd.read_csv(f"{DM}/OmicsExpressionProteinCodingGenesTPMLogp1.csv",
                          usecols=[allx.columns[0]] + sel[a:a + CH]
                          ).set_index(allx.columns[0]).reindex(common)
        for c in blk.columns:
            v = blk[c]
            m = dep.notna() & v.notna()
            if m.sum() >= 100:
                null.append(spearmanr(v[m], dep[m]).statistic)
    null = np.asarray([x for x in null if np.isfinite(x)])
    emp_p = float((null >= rs.statistic).mean())
    s3 = bool(rs.statistic > 0 and rs.pvalue < 0.05 and emp_p < 0.05)

    # 陰性對照：隨機基因對（不應通過）
    neg_r = float(np.median(null))
    out["S3_family_redundancy"] = {
        "rule": "MAFF+MAFK 表現與 MAFG gene effect 正相關、p < 0.05、經驗 p < 0.05 "
                "→ 冗餘成立，「MafG 是唯一必要夥伴」必須降級",
        "partners_used": fam_present,
        "spearman_r": float(rs.statistic), "pvalue": float(rs.pvalue),
        "n_lines": int(ok.sum()),
        "null_n": int(len(null)), "null_median": neg_r,
        "null_q95": float(np.quantile(null, 0.95)),
        "empirical_p": emp_p,
        "redundancy_supported": s3}
    log("  經驗虛無：n = %d、中位 %+.4f、95 分位 %+.4f → 經驗 p = %.4f"
        % (len(null), neg_r, np.quantile(null, 0.95), emp_p))
    log("  S3 冗餘 %s" % ("成立 → 設計理念必須改寫" if s3 else "不成立 → 設計理念撐得住"))

    # 補充：逐一看 MAFF 與 MAFK
    per = {}
    for g in fam_present:
        v = X.loc[common, g]
        m = dep.notna() & v.notna()
        r2 = spearmanr(v[m], dep[m])
        per[g] = {"spearman_r": float(r2.statistic), "pvalue": float(r2.pvalue)}
        log("    單獨 %-5s r = %+.4f, p = %.3g" % (g, r2.statistic, r2.pvalue))
    out["S3_family_redundancy"]["per_partner"] = per

    # ------------------------------------------------------------ S4 ---
    log("\n=== S4：表現廣度 ===")
    lin = mod["OncotreeLineage"].astype(str).reindex(X.index)
    med = X.groupby(lin)[TARGET].median().dropna().sort_values()
    frac = float((med > 1.0).mean())
    out["S4_expression_breadth"] = {
        "rule": "> 80% 的譜系中位表現 > 1 (log2 TPM+1) → 普遍表現，全身毒性風險",
        "n_lineages": int(len(med)),
        "frac_lineages_expressed": frac,
        "by_lineage": {str(k): float(v) for k, v in med.items()},
        "ubiquitous": bool(frac > 0.80)}
    log("  %d 個譜系，中位表現 > 1 的佔 %.1f%% → %s"
        % (len(med), frac * 100, "普遍表現" if frac > 0.80 else "有組織選擇性"))
    log("  最低 3 個：%s" % {str(k): round(float(v), 2) for k, v in med.head(3).items()})
    log("  最高 3 個：%s" % {str(k): round(float(v), 2) for k, v in med.tail(3).items()})

    # ------------------------------------------------------------ 總結 ---
    blockers = []
    if s1_constrained:
        blockers.append("S1_highly_constrained")
    if s2_broad:
        blockers.append("S2_broadly_essential")
    if s3:
        blockers.append("S3_family_redundant")
    out["verdict"] = {
        "blockers": blockers,
        "call": ("design_rationale_must_change" if s3 else
                 "no_blocker_found" if not blockers else "qualified")}
    log("\n*** 判定：%s ***" % out["verdict"]["call"])
    if s3:
        log("    `MAFF`／`MAFK` 看起來可以替代 `MAFG` ——")
        log("    「讀取因子冗餘但收斂到唯一一個必要夥伴」這個論述必須改寫。")
        log("    這直接影響第四節的設計理念，比任何結合劑分數都重要。")
    elif not blockers:
        log("    三道關卡都沒有否決 `MAFG` 當標的。")

    json.dump(out, open(f"{RES}/234_mafg_dossier.json", "w"),
              indent=2, ensure_ascii=False)
    log("\n→ results/234_mafg_dossier.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
