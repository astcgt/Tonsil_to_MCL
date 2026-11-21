#!/usr/bin/env python
"""
步驟 96：PIONEER2 介面資料分析。

我們要回答的問題
----------------
Q1. MafG 的交互作用夥伴有哪些？PIONEER2 預測的介面殘基，和我們從
    晶體結構 7X5E 取出的介面是否一致？（外部的第三方介面定義）
Q2. 我們設計的結合劑所瞄準的殘基，是不是 MafG 用來接「很多夥伴」的
    共用介面？若是，一個結合劑就能擋掉多條路徑；若不是，專一性較好但覆蓋少。
Q3. 軸線上哪些蛋白質的交互作用最多？哪些介面被多個夥伴共用（樞紐介面）？
Q4. FSP1 的失敗能不能從介面資料得到第二個解釋？

資料格式
--------
source=PDB      → protein1_interface 是二元介面殘基編號（逗號分隔，1-based）
source=PIONEER2 → 每殘基預測機率，格式 "殘基:機率;殘基:機率;..."
"""
import os, sys, gzip, json
from collections import defaultdict, Counter
import numpy as np
import pandas as pd

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
RES = f"{BASE}/results"
PROB_CUT = 0.50          # PIONEER2 機率門檻（作者網站預設觀念：>0.5 視為介面）

ids = pd.read_csv(f"{RES}/94_uniprot_ids.csv")
U2G = dict(zip(ids.uniprot, ids.gene))
AXIS = set(ids.uniprot)

# 我們自己從晶體結構定出的 MafG 介面（來自 scripts/56_mafg_target_def.py）
MAFG_DIMER = [75, 76, 79, 80, 82, 83, 86, 87, 89, 90, 93, 94, 97, 100, 101, 104, 107, 108, 111]
MAFG_DNA = [34, 53, 56, 57, 58, 60, 61, 62, 64, 65, 67, 68, 69, 71, 72]
# 結合劑實際接觸的熱點（來自 77_interface_compare.csv 的最佳設計）
BINDER_TARGET = MAFG_DNA


def log(*a):
    print("[96]", *a, flush=True)


def parse_iface(s, kind):
    """回傳介面殘基編號集合。"""
    if not s or s == "-":
        return set()
    if kind == "PDB":
        out = set()
        for x in s.split(","):
            x = x.strip()
            if x.isdigit():
                out.add(int(x))
        return out
    out = set()
    for pair in s.split(";"):
        if ":" not in pair:
            continue
        r, p = pair.split(":", 1)
        try:
            if float(p) >= PROB_CUT:
                out.add(int(r))
        except ValueError:
            continue
    return out


rows = []
with gzip.open(f"{RES}/95_pioneer2_axis.tsv.gz", "rt") as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) < 5:
            continue
        p1, p2, src, i1, i2 = f[0], f[1], f[2], f[3], f[4]
        for a, b, ia in [(p1, p2, i1), (p2, p1, i2)]:
            if a not in AXIS:
                continue
            res = parse_iface(ia, src)
            rows.append(dict(gene=U2G[a], self_uniprot=a, partner=b,
                             partner_gene=U2G.get(b, ""), source=src,
                             n_iface=len(res), residues=sorted(res)))
d = pd.DataFrame(rows)
log(f"展開成 {len(d):,} 筆「軸線蛋白質 → 夥伴」記錄")
log(f"來源分佈: {d.source.value_counts().to_dict()}")

# ---------------------------------------------------------------- Q3
log("\n" + "=" * 78)
log("=== Q3：軸線蛋白質的交互作用數量 ===")
g = (d.groupby("gene")
     .agg(n_partners=("partner", "nunique"),
          n_pdb=("source", lambda s: (s == "PDB").sum()),
          n_pred=("source", lambda s: (s == "PIONEER2").sum()),
          median_iface=("n_iface", "median"))
     .sort_values("n_partners", ascending=False))
log("\n" + g.to_string())
g.to_csv(f"{RES}/96_partner_counts.csv")

# ---------------------------------------------------------------- Q1
log("\n" + "=" * 78)
log("=== Q1：MafG 的介面，PIONEER2 vs 我們的晶體結構定義 ===")
m = d[d.gene == "MAFG"]
log(f"MafG 有 {m.partner.nunique()} 個夥伴（PDB {(m.source=='PDB').sum()}、預測 {(m.source=='PIONEER2').sum()}）")
if len(m):
    cnt = Counter()
    for _, r in m.iterrows():
        cnt.update(r.residues)
    tot = m.partner.nunique()
    freq = pd.DataFrame([dict(pos=p, n_partners=c, frac=round(c / tot, 3),
                              in_dimer=p in MAFG_DIMER, in_dna=p in MAFG_DNA)
                         for p, c in cnt.items()]).sort_values("n_partners", ascending=False)
    freq.to_csv(f"{RES}/96_mafg_residue_frequency.csv", index=False)
    log(f"\nMafG 被最多夥伴使用的前 20 個殘基：")
    log(freq.head(20).to_string(index=False))

    from scipy.stats import mannwhitneyu
    allpos = set(freq.pos)
    dim = freq[freq.in_dimer].n_partners
    dna = freq[freq.in_dna].n_partners
    oth = freq[~freq.in_dimer & ~freq.in_dna].n_partners
    log(f"\n各區域被夥伴使用的次數（中位數）：")
    log(f"  二聚面（我們的定義，{len(dim)} 個殘基在資料中）: {dim.median() if len(dim) else float('nan')}")
    log(f"  DNA 結合面（結合劑標的，{len(dna)} 個）:        {dna.median() if len(dna) else float('nan')}")
    log(f"  其餘位置（{len(oth)} 個）:                      {oth.median() if len(oth) else float('nan')}")
    for nm, v in [("二聚面", dim), ("DNA 面", dna)]:
        if len(v) >= 3 and len(oth) >= 3:
            u = mannwhitneyu(v, oth, alternative="greater")
            log(f"  {nm} vs 其餘  Mann-Whitney p = {u.pvalue:.4g}")

    # ---------------------------------------------------------------- Q2
    log("\n" + "=" * 78)
    log("=== Q2：結合劑標的是共用介面還是專一介面？ ===")
    bt = freq[freq.pos.isin(BINDER_TARGET)]
    log(f"結合劑瞄準的 {len(BINDER_TARGET)} 個殘基中，{len(bt)} 個出現在 PIONEER2 資料裡")
    if len(bt):
        log(f"  這些殘基平均被 {bt.n_partners.mean():.1f} 個夥伴使用 "
            f"（佔全部 {tot} 個夥伴的 {100*bt.frac.mean():.1f}%）")
        log(f"  最常被共用的: " + ", ".join(
            f"{int(r.pos)}({int(r.n_partners)})" for _, r in bt.nlargest(6, "n_partners").iterrows()))
    # 哪些夥伴會被結合劑擋到（夥伴介面與結合劑標的重疊）
    hit = []
    for _, r in m.iterrows():
        ov = set(r.residues) & set(BINDER_TARGET)
        if ov:
            hit.append(dict(partner=r.partner, partner_gene=r.partner_gene, source=r.source,
                            n_overlap=len(ov), overlap=sorted(ov),
                            frac_of_partner_iface=round(len(ov) / max(len(r.residues), 1), 3)))
    h = pd.DataFrame(hit).sort_values("n_overlap", ascending=False)
    h.to_csv(f"{RES}/96_binder_blocked_partners.csv", index=False)
    log(f"\n結合劑標的與 {len(h)} 個夥伴的介面重疊（共 {tot} 個夥伴，{100*len(h)/tot:.0f}%）")
    if len(h):
        log("\n重疊最多的前 12 個夥伴：")
        show = h.head(12).copy()
        show["partner_gene"] = show.partner_gene.replace("", "-")
        log(show[["partner", "partner_gene", "source", "n_overlap",
                  "frac_of_partner_iface"]].to_string(index=False))

# ---------------------------------------------------------------- Q4
log("\n" + "=" * 78)
log("=== Q4：FSP1（AIFM2）的介面資料 ===")
fs = d[d.gene == "AIFM2"]
log(f"FSP1 有 {fs.partner.nunique()} 個夥伴")
if len(fs):
    c2 = Counter()
    for _, r in fs.iterrows():
        c2.update(r.residues)
    f2 = pd.DataFrame([dict(pos=p, n_partners=c) for p, c in c2.items()]).sort_values(
        "n_partners", ascending=False)
    f2.to_csv(f"{RES}/96_fsp1_residue_frequency.csv", index=False)
    log(f"被最多夥伴使用的前 12 個殘基: " + ", ".join(
        f"{int(r.pos)}({int(r.n_partners)})" for _, r in f2.head(12).iterrows()))
    # 和我們的輔因子口袋比較
    try:
        pocket = [int(x) for x in open(f"{RES}/62_fsp1_hotspots.txt").read().split()]
    except Exception:
        pocket = []
    if pocket:
        inp = f2[f2.pos.isin(pocket)]
        log(f"\n我們的輔因子口袋 {len(pocket)} 個殘基中，{len(inp)} 個出現在 PIONEER2 資料")
        log(f"  口袋殘基平均被 {inp.n_partners.mean() if len(inp) else 0:.1f} 個夥伴使用")
        log(f"  非口袋殘基平均被 {f2[~f2.pos.isin(pocket)].n_partners.mean():.1f} 個夥伴使用")

d.drop(columns=["residues"]).to_csv(f"{RES}/96_axis_interactions.csv", index=False)
log(f"\n寫出 96_axis_interactions.csv、96_partner_counts.csv、"
    f"96_mafg_residue_frequency.csv、96_binder_blocked_partners.csv")
log("DONE")
