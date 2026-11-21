#!/usr/bin/env python
"""
步驟 73：把「演化約束（ESM）」與「結構介面（共結晶）」疊在一起。

這一步要回答兩個對藥物設計至關重要的問題：
  Q1 **哪些殘基同時是「結構上的介面」又是「演化上不能動的」？**
     → 這些就是最該打的位置：藥物結合在這裡，腫瘤很難靠突變逃脫
       （因為突變這些位置會把蛋白質自己弄壞）
  Q2 **介面上有哪些位置是可以容忍突變的？**
     → 這些是潛在的抗藥性逃逸口，設計時要避免只依賴它們

輸出 results/73_druggable_residues.csv
"""
import os, sys, json
import numpy as np, pandas as pd

BASE="/raid/jupyter-alz.ee09/Tonsil_to_MCL"; RES=f"{BASE}/results"
def log(*a): print("[73]",*a,flush=True)

cons=pd.read_csv(f"{RES}/71_esm_position_constraint.csv")
log("ESM 位置約束表:", cons.shape, "基因:", sorted(cons.gene.unique()))

ifaces={}
# MAFG：兩個介面（來自 7X5E）
rows=[l.strip() for l in open(f"{RES}/56_mafg_hotspots.txt") if l.strip() and not l.startswith("#")]
sec=None; mg=[]
for l in open(f"{RES}/56_mafg_hotspots.txt"):
    l=l.strip()
    if l.startswith("## "): sec="dimer" if "二聚" in l else "dna"; continue
    if not l or l.startswith("#") or "=" in l: continue
    p=l.split(",")
    if len(p)>=4: mg.append((int(p[0]), p[1], float(p[2]), p[3]))
ifaces["MAFG"]=pd.DataFrame(mg, columns=["pos","aa","dist","iface"])
# FSP1：輔因子口袋（來自 8YOQ）
fs=[]
for l in open(f"{RES}/62_fsp1_hotspots.txt"):
    l=l.strip()
    if not l or l.startswith("#") or "=" in l: continue
    p=l.split(",")
    if len(p)>=3: fs.append((int(p[0]), p[1], float(p[2]), "cofactor_pocket"))
ifaces["AIFM2"]=pd.DataFrame(fs, columns=["pos","aa","dist","iface"])

out=[]
for g,IF in ifaces.items():
    c=cons[cons.gene==g].copy()
    if not len(c): log(f"{g}: 無 ESM 資料，跳過"); continue
    c["pct"]=c.mean_llr.rank(pct=True)*100        # 百分位，越小越受約束
    m=IF.merge(c[["pos","wt","mean_llr","pct"]], on="pos", how="left")
    chk=(m.aa==m.wt)
    log(f"\n=== {g}：介面殘基 {len(m)} 個，與 ESM 序列對得上的 {int(chk.sum())} 個 ===")
    if chk.sum()<len(m)*0.8:
        log(f"  ⚠ 殘基編號可能有偏移（共結晶 vs UniProt），不一致範例: "
            + ", ".join(f"{int(r.pos)}:{r.aa}vs{r.wt}" for _,r in m[~chk].head(5).iterrows()))
    m["gene"]=g
    m["constrained"]=m.pct<=25        # 全蛋白質中最受約束的前 25%
    m["tolerant"]=m.pct>=60
    out.append(m)
    base=c.pct.describe()
    log(f"  介面殘基的約束百分位中位數 = {m.pct.median():.1f}（全蛋白質基準 50.0）")
    log(f"  其中「高約束（前 25%）」的介面殘基：{int(m.constrained.sum())} 個")
    log("    " + " ".join(f"{r.aa}{int(r.pos)}" for _,r in m[m.constrained].sort_values('pct').iterrows()))
    log(f"  「可容忍（後 40%）」的介面殘基（潛在抗藥性逃逸口）：{int(m.tolerant.sum())} 個")
    log("    " + " ".join(f"{r.aa}{int(r.pos)}" for _,r in m[m.tolerant].sort_values('pct',ascending=False).iterrows()))
    # 統計檢定：介面殘基是否比隨機位置更受約束
    from scipy.stats import mannwhitneyu
    other=c[~c.pos.isin(m.pos)]
    if len(other)>10:
        u,p=mannwhitneyu(m.mean_llr.dropna(), other.mean_llr, alternative="less")
        log(f"  介面 vs 非介面的約束差異：Mann-Whitney p = {p:.3g}"
            f"（介面中位 {m.mean_llr.median():.2f} vs 非介面 {other.mean_llr.median():.2f}）")

if out:
    A=pd.concat(out, ignore_index=True)
    A.to_csv(f"{RES}/73_druggable_residues.csv", index=False)
    log(f"\n寫出 {RES}/73_druggable_residues.csv")
    log("\n=== 最優先的藥物接觸點（介面 + 高演化約束，依約束度排序）===")
    pd.set_option("display.width",200)
    top=A[A.constrained].sort_values(["gene","pct"])
    log("\n"+top[["gene","pos","aa","iface","dist","mean_llr","pct"]].round(2).to_string(index=False))
log("DONE")
