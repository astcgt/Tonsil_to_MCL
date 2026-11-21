#!/usr/bin/env python
"""
步驟 62：定義 FSP1 (AIFM2) 的可成藥表面。

為什麼是 FSP1
-------------
在本研究發現的「胱胺酸–穀胱甘肽–鐵死亡防禦」軸線上，FSP1 是最適合做蛋白質藥物的節點：
  * DepMap MCL 選擇性依賴 z = −2.60（三株 MCL 細胞株都必需）
  * 文獻：AIFM2/FSP1 + mantle cell lymphoma = **0 篇**；+ lymphoma 只有 5 篇
  * 它是**可溶性酵素**（不像 SLC7A11 是 12 次跨膜蛋白），有 1.79 Å 晶體結構
  * 它走的是 GPX4 之外的獨立路徑，所以和既有的 GPX4 抑制劑不重疊

結構：PDB 8YOQ（hFSP1 + 6-OH-FAD + NADP，1.79 Å）
可成藥表面 = FAD 與 NAD(P) 的輔因子口袋周邊（催化核心）
"""
import os, sys
from collections import defaultdict
import numpy as np, pandas as pd

BASE="/raid/jupyter-alz.ee09/Tonsil_to_MCL"; ST=f"{BASE}/data/structures"; RES=f"{BASE}/results"
PDB=f"{ST}/8yoq.pdb"; CUT=4.5
AA3={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I',
     'LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
def log(*a): print("[62]",*a,flush=True)

prot=defaultdict(list); het=defaultdict(list); keep=[]
for l in open(PDB):
    if l.startswith("ATOM"):
        r=l[17:20].strip()
        if r in AA3:
            prot[l[21]].append((int(l[22:26]), r,
                                (float(l[30:38]),float(l[38:46]),float(l[46:54]))))
            keep.append(l)
    elif l.startswith("HETATM"):
        r=l[17:20].strip()
        if r not in ("HOH","WAT"):
            het[r].append((l[21], int(l[22:26]),
                           (float(l[30:38]),float(l[38:46]),float(l[46:54]))))
chains=sorted(prot, key=lambda c: -len(prot[c]))
log("蛋白質鏈:", {c: len({x[0] for x in prot[c]}) for c in chains})
log("配體:", {k: len(v) for k,v in het.items()})

main=chains[0]
P=np.array([a[2] for a in prot[main]])
res=[(a[0],a[1]) for a in prot[main]]
pockets={}
for lig,atoms in het.items():
    L=np.array([a[2] for a in atoms if a[0]==main])
    if len(L)==0: continue
    d=np.sqrt(((P[:,None,:]-L[None,:,:])**2).sum(-1)).min(axis=1)
    near={}
    for i,v in enumerate(d):
        n,r=res[i]
        if v<=CUT and (n not in near or v<near[n][1]): near[n]=(r,v)
    pockets[lig]=near
    log(f"\n=== {lig} 口袋（{len(near)} 個殘基，鏈 {main}）===")
    log("  " + " ".join(f"{AA3[r]}{n}" for n,(r,v) in sorted(near.items())))

allp={}
for lig,near in pockets.items():
    for n,(r,v) in near.items():
        if n not in allp or v<allp[n][1]: allp[n]=(r,v)
core=[n for n,(r,v) in allp.items() if v<=3.8]
log(f"\n=== 合併後的催化口袋：{len(allp)} 個殘基，核心熱點 {len(core)} 個 ===")
log("  核心: " + " ".join(f"{AA3[allp[n][0]]}{n}" for n in sorted(core)))

with open(f"{RES}/62_fsp1_hotspots.txt","w") as fh:
    fh.write(f"# FSP1 (AIFM2) 輔因子口袋，來源 PDB 8YOQ 鏈 {main}\n")
    for n,(r,v) in sorted(allp.items()): fh.write(f"{n},{AA3[r]},{v:.2f}\n")
    fh.write("\nhotspots=" + ",".join(f"{main}{n}" for n in sorted(core)) + "\n")
with open(f"{RES}/62_fsp1_target.pdb","w") as fh:
    for l in keep:
        if l[21]==main: fh.write(l)
    fh.write("END\n")
log(f"\n寫出 {RES}/62_fsp1_target.pdb（鏈 {main}）與熱點清單")
seq="".join(AA3[r] for n,r in sorted({(n,r) for n,r,_ in prot[main]}))
log(f"FSP1 鏈長 {len(seq)} 殘基")
log("DONE")
