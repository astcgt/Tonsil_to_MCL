#!/usr/bin/env python
"""
步驟 56：定義 MAFG 的兩個可成藥表面，並產出對應的藥物候選。

為什麼標 MAFG
-------------
它是本研究以 RNA+ATAC 重接線分析找出的**樞紐**：
在 chr17q25 獲得 19 個全新的調控元件（正常相關性 ≈0 → 腫瘤 0.28–0.53，0 條斷裂），
而它所驅動的抗氧化／鐵死亡防禦下游（SLC7A11、SLC3A2、GCLM、GCLC、AIFM2）
正好是 DepMap 中 MCL 選擇性依賴最強的一群。
文獻查證：MAFG + mantle cell lymphoma = 0 篇，MAFG + lymphoma 全部只有 2 篇。

兩個表面（取自 PDB 7X5E，NRF2–MafG–DNA 共結晶 2.3 Å）
-----------------------------------------------------
介面 1：MafG ↔ NRF2 的白胺酸拉鍊二聚體介面
        → 擋住它，NRF2/BACH1 就無法與 MafG 組成有功能的異二聚體
介面 2：MafG 的鹼性區（直接接觸 DNA）
        → 擋住它，即使二聚體形成也無法結合 MARE/ARE

另外直接從共結晶取出 MARE 序列，做成誘餌雙股寡核苷酸。
"""
import os, sys
from collections import defaultdict
import numpy as np, pandas as pd

BASE="/raid/jupyter-alz.ee09/Tonsil_to_MCL"; ST=f"{BASE}/data/structures"; RES=f"{BASE}/results"
PDB=f"{ST}/7x5e.pdb"; CUT=4.5
AA3={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I',
     'LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
DNA={"DA","DT","DG","DC"}
COMP=str.maketrans("ACGT","TGCA")
def log(*a): print("[56]", *a, flush=True)

atoms=defaultdict(list); lines=defaultdict(list)
for l in open(PDB):
    if not l.startswith("ATOM"): continue
    r=l[17:20].strip(); c=l[21]; n=int(l[22:26])
    xyz=(float(l[30:38]),float(l[38:46]),float(l[46:54]))
    kind = "prot" if r in AA3 else ("dna" if r in DNA else None)
    if kind: atoms[c].append((n,r,xyz,kind)); lines[c].append(l)

def contacts(cA, cB, cut=CUT):
    A=[a for a in atoms[cA]]; B=[a for a in atoms[cB]]
    PA=np.array([a[2] for a in A]); PB=np.array([a[2] for a in B])
    out={}
    for i0 in range(0,len(PA),3000):
        blk=PA[i0:i0+3000]
        d=np.sqrt(((blk[:,None,:]-PB[None,:,:])**2).sum(-1)).min(axis=1)
        for j,v in enumerate(d):
            n,r,_,_=A[i0+j]
            if n not in out or v<out[n][1]: out[n]=(r,v)
    return {n:(r,v) for n,(r,v) in out.items() if v<=cut}

log("MafG = 鏈 A（殘基 23–121）；NRF2 = 鏈 B；DNA = 鏈 C/D")
iface_nrf2 = contacts("A","B")
iface_dna  = {}
for d in ("C","D"):
    for n,(r,v) in contacts("A",d).items():
        if n not in iface_dna or v<iface_dna[n][1]: iface_dna[n]=(r,v)

def fmt(dd): return " ".join(f"{AA3.get(r,'X')}{n}" for n,(r,v) in sorted(dd.items()))
log(f"\n=== 介面 1：MafG ↔ NRF2 二聚體介面（{len(iface_nrf2)} 個殘基）===")
log("  " + fmt(iface_nrf2))
log(f"\n=== 介面 2：MafG ↔ DNA 鹼性區（{len(iface_dna)} 個殘基）===")
log("  " + fmt(iface_dna))
core1=[n for n,(r,v) in iface_nrf2.items() if v<=3.8]
core2=[n for n,(r,v) in iface_dna.items() if v<=3.5]
log(f"\n核心熱點：二聚體介面 {len(core1)} 個、DNA 介面 {len(core2)} 個")

with open(f"{RES}/56_mafg_hotspots.txt","w") as fh:
    fh.write("# MAFG 可成藥表面（來源 PDB 7X5E，鏈 A = MafG）\n")
    fh.write("## 介面1 MafG-NRF2 二聚體\n")
    for n,(r,v) in sorted(iface_nrf2.items()): fh.write(f"{n},{AA3.get(r,'X')},{v:.2f},dimer\n")
    fh.write("## 介面2 MafG-DNA\n")
    for n,(r,v) in sorted(iface_dna.items()): fh.write(f"{n},{AA3.get(r,'X')},{v:.2f},dna\n")
    fh.write("\n# RFdiffusion 熱點\n")
    fh.write("dimer_hotspots=" + ",".join(f"A{n}" for n in sorted(core1)) + "\n")
    fh.write("dna_hotspots="   + ",".join(f"A{n}" for n in sorted(core2)) + "\n")

# 只保留 MafG 單鏈作為設計標的
with open(f"{RES}/56_mafg_target.pdb","w") as fh:
    for l in lines["A"]: fh.write(l)
    fh.write("END\n")
log(f"寫出 {RES}/56_mafg_target.pdb（MafG 單體 {len(lines['A'])} 原子）")

# ---------------- MARE 誘餌寡核苷酸 ----------------
mare = "TTGCTGACTCATCAT"     # 7X5E 共結晶中的 CsMBE
rows=[]
def add(nm, top, note):
    rows.append(dict(name=nm, top_strand=top, bottom_strand=top.translate(COMP)[::-1],
                     length=len(top), gc=round((top.count("G")+top.count("C"))/len(top),2), note=note))
add("MARE_from_7X5E", mare, "直接取自 NRF2–MafG 共結晶所結合的 DNA")
add("MARE_x3_tandem", mare*3, "三串聯，提高吸附容量")
add("ARE_canonical", "TGACTCAGCA", "典型抗氧化反應元件核心 (TGACnnnGC)")
add("ARE_x4_tandem", "TGACTCAGCA"*4, "四串聯 ARE")
add("MARE_hairpin", mare + "TTCG" + mare.translate(COMP)[::-1],
    "髮夾式單股（TTCG 四環），穩定度高於雙股")
add("scramble_control", "TACTGTCATGCTACT", "陰性對照：組成相近但無 TGACTCA 核心")
df=pd.DataFrame(rows); df.to_csv(f"{RES}/56_mafg_decoys.csv", index=False)
pd.set_option("display.width",220)
log(f"\n=== MARE 誘餌寡核苷酸（{len(df)} 條）===")
log("\n"+df[["name","top_strand","length","gc","note"]].to_string(index=False))
log("DONE")
