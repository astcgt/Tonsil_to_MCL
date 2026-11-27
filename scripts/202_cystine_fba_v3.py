#!/usr/bin/env python
"""步驟 202（B3 第三版）：半胱胺酸要匱乏到什麼程度，xCT 才會變成瓶頸？

這個問題從哪裡來
----------------
DepMap 的 `SLC7A11` 依賴性只在 18 株（1.5%）出現，三株 MCL 是 0/3。
但 DepMap 的培養基（RPMI-1640）胱胺酸是 208 µM，遠高於人血漿的游離胱胺酸
（約 40–60 µM），而腫瘤微環境更低。所以「測不到依賴」可能只是
「在那個培養基裡測不到」。這一步用化學計量去定出那條界線在哪。

前兩版為什麼失敗（都是我的設計錯誤，不是程式 bug）
--------------------------------------------------
v1（步驟 200）只調了胱胺酸交換反應的上限。但 Human-GEM 的 SBML 把 **1,654 條
交換反應全部設成 (−1000, 1000)**，等於所有營養無限供應（成長率 138.5，
四種營養條件完全一樣）。關一條沒有意義。

v2（步驟 201）把 46 條「直接提供半胱胺酸」的交換關掉，但沒有關掉其餘 1,600 條，
所以四格成長率還是 138.5。全部敲除效應都是 0。

v3 的四個修正
-------------
1. **先關掉全部 1,654 條交換的攝取，再只開一個定義過的 RPMI-1640 培養基。**
2. **目標函數不能是「最大化 GSH 合成反應」。** γ-穀胺醯循環會讓 GSH 分解後
   再生出半胱胺酸、麩胺酸、甘胺酸，又合成回 GSH —— 淨投入為零的空轉，
   於是合成通量可以無限大（實測：關掉所有硫醇來源後仍有 249.8）。
   改成加一條**不可逆的需求反應** `GSH[c] -->` 當目標。
   因為 GSH 被不可逆移除，循環無法供應它，每一單位淨產出都需要一單位淨半胱胺酸。
3. **三道一定要先過的驗證門檻**（任何一道不過就 raise，不准往下解讀）：
   G1 陽性對照的**化學計量**：胱胺酸 0.208 + 甲硫胺酸 0.101 下，
      最大 GSH 淨產出必須等於 2×0.208 + 0.101 = 0.517
      （一個胱胺酸給兩個半胱胺酸；甲硫胺酸走轉硫途徑給一個）。
   G2 fail-stop：關掉所有硫醇來源與甲硫胺酸，產出必須是 0。
   G3 硫酸鹽不得成為半胱胺酸來源（人類沒有硫酸鹽同化途徑）。
      關掉硫酸鹽後產出不應改變。
   G4 **敲掉 `SLC7A11` 必須真的有效果**。這一道是第一次跑完才發現要加的：
      前三道全過，但敲掉 `SLC7A11` 在任何胱胺酸濃度下效果都是 0。
      原因是 Human-GEM 是**泛人類**模型，裡面還有三條 xCT 以外的胱胺酸入口
      （MAR08687、MAR06891、MAR06892），GPR 全都是 `SLC7A9 and SLC3A1`
      ——那是腎臟與腸道的 b⁰,⁺ 系統（rBAT / b⁰,⁺AT）。

4. **用本專案自己的表現資料把模型限定到 B 細胞情境**，而不是憑印象關反應：
      基因        MCL 實測   正常 B 期望     （log1p CP10K，`results/15_gene_residual_tumor.csv`）
      SLC7A11     0.2061    0.0009
      SLC3A2      0.2724    0.0946
      SLC3A1      0.0172    0.0273   ← 比 SLC7A11 低 12 倍
      SLC7A9      0.0043    0.0023   ← 比 SLC7A11 低 48 倍
   判定規則寫死在程式裡：兩個情境的平均表現都 < 0.05 就視為「這個細胞沒有這個基因」，
   把它敲掉。只對「把胞外胱胺酸/半胱胺酸送進胞內」的反應做這件事，
   不做全基因組的情境化（那會動到整個網絡，風險遠大於收益）。
   GPR 含 `or` 的反應一律保守處理：只報告、不動。

⚠️ FBA 是穩態化學計量：沒有動力學、沒有氧化壓力、沒有鐵死亡。
   它回答「這條路在化學計量上何時成為瓶頸」，不回答「細胞會不會死」。

輸出 results/202_*
"""
import os
import json
import warnings

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "4"
warnings.filterwarnings("ignore")

import re
import numpy as np
import pandas as pd

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
RES = f"{BASE}/results"
GEM = f"{BASE}/data/aux/gem/Human-GEM.xml"

# RPMI-1640 的濃度 (mM) 直接當攝取上限（交換反應的 lower_bound 取負值）
RPMI = {
    "MAR09034": 11.1,     # glucose
    "MAR09066": 1.15,     # arginine
    "MAR09062": 0.379,    # asparagine
    "MAR09070": 0.150,    # aspartate
    "MAR09363": 0.208,    # cystine   ← 要掃的變數之一
    "MAR09071": 0.136,    # glutamate
    "MAR09063": 2.05,     # glutamine
    "MAR09067": 0.133,    # glycine
    "MAR09038": 0.0968,   # histidine
    "MAR09039": 0.382,    # isoleucine
    "MAR09040": 0.382,    # leucine
    "MAR09041": 0.274,    # lysine
    "MAR09042": 0.101,    # methionine ← 要掃的變數之一
    "MAR09043": 0.0909,   # phenylalanine
    "MAR09068": 0.174,    # proline
    "MAR09069": 0.286,    # serine
    "MAR09044": 0.168,    # threonine
    "MAR09045": 0.0245,   # tryptophan
    "MAR09064": 0.111,    # tyrosine
    "MAR09046": 0.171,    # valine
}
FREE = ["MAR09048",  # O2
        "MAR09047",  # H2O
        "MAR09079",  # H+
        "MAR09072",  # Pi
        "MAR09081",  # K+
        "MAR09077",  # Na+
        "MAR09082",  # Ca2+
        "MAR13072",  # Mg2+
        "MAR09078",  # HCO3-
        "MAR09074",  # sulfate
        "MAR09076"]  # Fe2+
CYSTINE = "MAR09363"
METHIONINE = "MAR09042"
SULFATE = "MAR09074"
THIOL_PAT = r"cystein|cystin|cystyl|glutathion|taurine|hypotaurine"

PANEL = ["SLC7A11", "SLC3A2", "GCLC", "GCLM", "GSS", "GPX4", "AIFM2",
         "CBS", "CTH", "MAT2A", "G6PD", "PGD", "SLC2A1"]
# 轉硫途徑：甲硫胺酸 → 高半胱胺酸 → 胱硫醚 → 半胱胺酸。
# 同一套表現判定（兩情境平均 < 0.05）說它在 B 細胞與 MCL 都等於不存在：
#   CBS  0.0003 / 0.0002（46 個正常狀態中 0 個 > 0.1）
#   CBSL 0.0001 / 0.0000（0/46）
#   CTH  0.0090 / 0.0112（1/46）
# 所以「泛人類」模型與「B 細胞」模型在這條路上是不同的，兩個都跑、分開報。
TRANSSULF = ["CBS", "CBSL", "CTH"]
CYSTEINE_EX = None   # 執行時填入游離半胱胺酸的交換反應 id
ABSENT_CUT = 0.05      # log1p CP10K；兩個情境都低於這個值才算「沒有」


def bcell_expression():
    """回傳 {基因: (MCL 實測, 正常 B 期望)}，單位 log1p CP10K。"""
    d = pd.read_csv(f"{RES}/15_gene_residual_tumor.csv", index_col=0)
    return {g: (float(d.loc[g, "mean_observed"]),
                float(d.loc[g, "mean_expected"])) for g in d.index}


def import_routes(m, names=("cystine",)):
    """所有把胞外**胱胺酸**送進胞內的反應（交換反應不算）。

    一開始也掃了 cysteine，結果錯得很具體：Human-GEM 裡有大量
    `cysteine + proline <=> cysteine + proline` 這種反應 —— 那是脯胺酸轉運蛋白，
    半胱胺酸只是逆向受質。於是 SLC1A5、SLC6A14、SLC38A5、CLTRN 這些
    和胱胺酸入口無關的基因也被敲掉，動到了與問題無關的網絡。
    游離半胱胺酸的攝取本來就被培養基關掉（它在 thiol 清單裡），
    所以這裡只需要、也只應該處理胱胺酸。
    """
    out = []
    for nm in names:
        me = [x for x in m.metabolites if x.name == nm and x.compartment == "e"]
        mc = [x for x in m.metabolites if x.name == nm and x.compartment == "c"]
        if not me or not mc:
            continue
        for r in me[0].reactions:
            if r in mc[0].reactions:
                out.append((nm, r))
    return out


def contextualise(m, expr):
    """把 B 細胞不表現的轉運次單元敲掉，回傳一份可稽核的紀錄。"""
    rec, dropped = [], set()
    for nm, r in import_routes(m):
        rule = r.gene_reaction_rule
        gs = sorted({g.name for g in r.genes})
        ev = {g: expr.get(g, (float("nan"),) * 2) for g in gs}
        pure_and = (" or " not in rule.lower()) and bool(gs)
        absent = [g for g, (t, n) in ev.items()
                  if t == t and t < ABSENT_CUT and n < ABSENT_CUT]
        act = "—"
        if pure_and and absent:
            act = "敲掉 " + "/".join(absent)
            dropped.update(absent)
        elif absent:
            act = "GPR 含 or，保守不動"
        rec.append({"metabolite": nm, "reaction": r.id,
                    "gpr_symbols": " and ".join(gs) if pure_and else rule,
                    "pure_and": pure_and,
                    "expression": {g: {"mcl": ev[g][0], "normal_b": ev[g][1]}
                                   for g in gs},
                    "absent_genes": absent, "action": act})
    for g in sorted(dropped):
        gid = gid_of(m, g)
        if gid:
            m.genes.get_by_id(gid).knock_out()
    return rec, sorted(dropped)


def log(*a):
    print("[202]", *a, flush=True)


def build():
    import cobra
    m = cobra.io.read_sbml_model(GEM)
    thiol = [r.id for r in m.exchanges
             if re.search(THIOL_PAT, r.name.lower())]
    g = [x for x in m.metabolites if x.name == "GSH" and x.compartment == "c"][0]
    dm = cobra.Reaction("DM_GSH_c")
    dm.lower_bound, dm.upper_bound = 0.0, 1000.0
    dm.add_metabolites({g: -1.0})
    m.add_reactions([dm])
    return m, thiol


def cysteine_exchange(m):
    for r in m.exchanges:
        if r.name.strip().lower() == "exchange of cysteine":
            return r.id
    return None


def set_medium(m, thiol, cystine, methionine, close_all_thiol=False,
               cysteine=0.0):
    """關掉全部交換的攝取，只開 RPMI；非胱胺酸的硫醇來源一律關閉。"""
    for r in m.exchanges:
        r.lower_bound = 0.0
    for rid, c in RPMI.items():
        m.reactions.get_by_id(rid).lower_bound = -float(c)
    for rid in FREE:
        m.reactions.get_by_id(rid).lower_bound = -1000.0
    for rid in thiol:                       # 其他硫醇旁路永遠關著
        if rid != CYSTINE:
            m.reactions.get_by_id(rid).lower_bound = 0.0
    m.reactions.get_by_id(CYSTINE).lower_bound = -float(cystine)
    m.reactions.get_by_id(METHIONINE).lower_bound = -float(methionine)
    if cysteine and CYSTEINE_EX:
        m.reactions.get_by_id(CYSTEINE_EX).lower_bound = -float(cysteine)
    if close_all_thiol:
        for rid in thiol:
            m.reactions.get_by_id(rid).lower_bound = 0.0


def solve(m):
    s = m.optimize()
    return float(s.objective_value) if s.status == "optimal" else float("nan")


def gid_of(m, sym):
    hit = [x.id for x in m.genes if x.name == sym]
    return hit[0] if hit else None


def main():
    global CYSTEINE_EX
    log("載入 Human-GEM …")
    m, thiol = build()
    CYSTEINE_EX = cysteine_exchange(m)
    log("游離半胱胺酸的交換反應：%s" % CYSTEINE_EX)
    m.objective = "DM_GSH_c"
    log("反應 %d、代謝物 %d、基因 %d；交換反應 %d 條；含硫醇的交換 %d 條"
        % (len(m.reactions), len(m.metabolites), len(m.genes),
           len(m.exchanges), len(thiol)))

    expr = bcell_expression()
    routes, dropped = contextualise(m, expr)
    log("\n把模型限定到 B 細胞情境（只動胱胺酸/半胱胺酸的入口反應）：")
    for r in routes:
        log("  %-10s %-34s %s" % (r["reaction"], r["gpr_symbols"][:34], r["action"]))
    log("  → 敲掉的基因：%s" % (", ".join(dropped) or "（無）"))

    gates = {}

    # G1：陽性對照的化學計量
    with m:
        set_medium(m, thiol, 0.208, 0.101)
        v = solve(m)
    expect = 2 * 0.208 + 0.101
    gates["G1_stoichiometry"] = {
        "observed": v, "expected": expect,
        "reasoning": "一個胱胺酸給兩個半胱胺酸，甲硫胺酸經轉硫途徑給一個",
        "pass": bool(abs(v - expect) < 1e-3)}
    log("G1 化學計量：實測 %.4f vs 預期 %.4f（2×0.208+0.101）→ %s"
        % (v, expect, "過" if gates["G1_stoichiometry"]["pass"] else "不過"))

    # G2：fail-stop
    with m:
        set_medium(m, thiol, 0.0, 0.0, close_all_thiol=True)
        v = solve(m)
    gates["G2_failstop"] = {"observed": v, "expected": 0.0,
                            "pass": bool(abs(v) < 1e-6)}
    log("G2 fail-stop（所有硫醇來源＋甲硫胺酸關閉）：%.6f → %s"
        % (v, "過" if gates["G2_failstop"]["pass"] else "不過"))

    # G3：硫酸鹽不能當半胱胺酸來源
    with m:
        set_medium(m, thiol, 0.208, 0.101)
        base = solve(m)
        m.reactions.get_by_id(SULFATE).lower_bound = 0.0
        nosulf = solve(m)
    gates["G3_no_sulfate_assimilation"] = {
        "with_sulfate": base, "without_sulfate": nosulf,
        "pass": bool(abs(base - nosulf) < 1e-6)}
    log("G3 硫酸鹽：有 %.4f / 無 %.4f → %s"
        % (base, nosulf, "過" if gates["G3_no_sulfate_assimilation"]["pass"] else "不過"))

    # G4：敲掉 SLC7A11 必須真的有效果（只有胱胺酸、沒有甲硫胺酸的條件下）
    with m:
        set_medium(m, thiol, 0.208, 0.0)
        wt = solve(m)
        with m:
            m.genes.get_by_id(gid_of(m, "SLC7A11")).knock_out()
            ko = solve(m)
    gates["G4_slc7a11_is_required"] = {
        "wt": wt, "slc7a11_ko": ko,
        "reasoning": "胱胺酸是唯一硫來源、甲硫胺酸關閉時，xCT 應為唯一入口",
        "pass": bool(wt > 1e-6 and abs(ko) < 1e-6)}
    log("G4 敲掉 SLC7A11（只有胱胺酸）：野生型 %.4f → 敲除 %.6f → %s"
        % (wt, ko, "過" if gates["G4_slc7a11_is_required"]["pass"] else "不過"))

    if not all(g["pass"] for g in gates.values()):
        json.dump({"gates": gates}, open(f"{RES}/202_cystine_fba.json", "w"),
                  indent=2, ensure_ascii=False)
        raise SystemExit("驗證門檻未通過，拒絕輸出任何敲除結果")
    log("三道門檻全部通過，開始敲除分析\n")

    # --- 2×2：胱胺酸 {充裕, 零} × 甲硫胺酸 {充裕, 零} ---------------------
    rows = []
    grid = [(0.208, 0.101, "RPMI（DepMap 的條件）"),
            (0.208, 0.0, "胱胺酸充裕 / 無甲硫胺酸"),
            (0.0, 0.101, "無胱胺酸 / 甲硫胺酸充裕"),
            (0.0, 0.0, "兩者皆無")]
    for cy, me_, label in grid:
        with m:
            set_medium(m, thiol, cy, me_)
            wt = solve(m)
            for sym in PANEL:
                g = gid_of(m, sym)
                if g is None:
                    continue
                with m:
                    m.genes.get_by_id(g).knock_out()
                    ko = solve(m)
                rows.append(dict(cystine_mM=cy, methionine_mM=me_,
                                 condition=label, gene=sym, ensembl=g,
                                 gsh_wt=wt, gsh_ko=ko,
                                 ko_effect=ko - wt,
                                 ko_frac=(ko - wt) / wt if wt > 1e-9 else np.nan))
        log("%-26s 野生型 GSH 淨產出 = %.4f" % (label, wt))
    D = pd.DataFrame(rows)
    D.to_csv(f"{RES}/202_fba_knockout.csv", index=False)

    piv = D.pivot_table(index="gene", columns="condition",
                        values="ko_frac").reindex(PANEL)
    log("\n敲除造成的 GSH 淨產出相對變化（−1.00 ＝ 完全歸零）：")
    log("\n" + piv.round(3).to_string())

    # --- 掃描：胱胺酸濃度多低，xCT 才變成瓶頸 ----------------------------
    scan = []
    cy_grid = [0.208, 0.15, 0.10, 0.06, 0.04, 0.02, 0.01, 0.0]
    me_grid = [0.101, 0.05, 0.02, 0.0]
    for cy in cy_grid:
        for me_ in me_grid:
            with m:
                set_medium(m, thiol, cy, me_)
                wt = solve(m)
                g = gid_of(m, "SLC7A11")
                with m:
                    m.genes.get_by_id(g).knock_out()
                    ko = solve(m)
            scan.append(dict(cystine_mM=cy, methionine_mM=me_,
                             gsh_wt=wt, gsh_ko_slc7a11=ko,
                             frac_lost=(wt - ko) / wt if wt > 1e-9 else np.nan))
    S = pd.DataFrame(scan)
    S.to_csv(f"{RES}/202_fba_scan.csv", index=False)
    sp = S.pivot_table(index="cystine_mM", columns="methionine_mM",
                       values="frac_lost").sort_index(ascending=False)
    log("\n敲掉 SLC7A11 之後失去的 GSH 產能比例"
        "（列＝胱胺酸 mM、欄＝甲硫胺酸 mM）：")
    log("\n" + sp.round(3).to_string())

    # --- 情境 B：把轉硫途徑也按表現關掉（B 細胞 / MCL 的真實情況）----------
    ts_expr = {g: expr.get(g, (float("nan"),) * 2) for g in TRANSSULF}
    ts_absent = [g for g, (t, n) in ts_expr.items()
                 if t == t and t < ABSENT_CUT and n < ABSENT_CUT]
    log("\n情境 B：轉硫途徑的表現判定")
    for g in TRANSSULF:
        t, n = ts_expr[g]
        log("  %-5s MCL %.4f / 正常 B %.4f → %s"
            % (g, t, n, "視為不存在" if g in ts_absent else "保留"))
    scanB = []
    with m:
        for g in ts_absent:
            gi = gid_of(m, g)
            if gi:
                m.genes.get_by_id(gi).knock_out()
        for cy in [0.208, 0.10, 0.04, 0.01]:
            for me_ in [0.101, 0.0]:
                with m:
                    set_medium(m, thiol, cy, me_)
                    wt = solve(m)
                    with m:
                        m.genes.get_by_id(gid_of(m, "SLC7A11")).knock_out()
                        ko = solve(m)
                scanB.append(dict(cystine_mM=cy, methionine_mM=me_, gsh_wt=wt,
                                  gsh_ko_slc7a11=ko,
                                  frac_lost=(wt - ko) / wt if wt > 1e-9 else np.nan))
    B = pd.DataFrame(scanB)
    B.to_csv(f"{RES}/202_fba_bcell_context.csv", index=False)
    log("\n情境 B（無轉硫途徑）敲掉 SLC7A11 失去的 GSH 產能比例：")
    log("\n" + B.pivot_table(index="cystine_mM", columns="methionine_mM",
                             values="frac_lost").sort_index(ascending=False)
        .round(3).to_string())

    # --- 情境 C：培養基裡有游離半胱胺酸（血清、細胞交互餵養、胞外還原）------
    scanC = []
    if CYSTEINE_EX:
        with m:
            for g in ts_absent:
                gi = gid_of(m, g)
                if gi:
                    m.genes.get_by_id(gi).knock_out()
            for cysfree in [0.0, 0.005, 0.02, 0.05, 0.2]:
                with m:
                    set_medium(m, thiol, 0.208, 0.0, cysteine=cysfree)
                    wt = solve(m)
                    with m:
                        m.genes.get_by_id(gid_of(m, "SLC7A11")).knock_out()
                        ko = solve(m)
                scanC.append(dict(free_cysteine_mM=cysfree, gsh_wt=wt,
                                  gsh_ko_slc7a11=ko,
                                  frac_lost=(wt - ko) / wt if wt > 1e-9 else np.nan))
        C = pd.DataFrame(scanC)
        C.to_csv(f"{RES}/202_fba_free_cysteine.csv", index=False)
        log("\n情境 C（無轉硫途徑、胱胺酸 0.208）游離半胱胺酸對 xCT 依賴性的影響：")
        log("\n" + C.round(4).to_string(index=False))

    out = {"model": "Human-GEM", "objective": "DM_GSH_c（不可逆 GSH 需求）",
           "medium": "RPMI-1640（20 個成分）＋ 11 個自由無機物/氣體",
           "n_exchanges_closed": len(m.exchanges) - len(RPMI) - len(FREE),
           "n_thiol_exchanges_blocked": len(thiol) - 1,
           "context_specific_routes": routes,
           "context_dropped_genes": dropped,
           "absent_threshold_log1p_cp10k": ABSENT_CUT,
           "gates": gates,
           "transsulfuration_expression": {g: {"mcl": ts_expr[g][0],
                                               "normal_b": ts_expr[g][1]}
                                           for g in TRANSSULF},
           "transsulfuration_absent": ts_absent,
           "scenario_B_bcell_no_transsulfuration": scanB,
           "scenario_C_free_cysteine": scanC,
           "ko_frac_by_condition": json.loads(piv.round(4).to_json()),
           "slc7a11_scan": json.loads(sp.round(4).to_json()),
           "caveat": ("FBA 為穩態化學計量；無動力學、氧化壓力與鐵死亡。"
                      "回答的是「何時成為化學計量瓶頸」，不是「細胞會不會死」。")}
    json.dump(out, open(f"{RES}/202_cystine_fba.json", "w"),
              indent=2, ensure_ascii=False)
    log("\nDONE → results/202_*")


if __name__ == "__main__":
    main()
