#!/usr/bin/env python
"""
步驟 36：定義 IRF4 的可成藥表面，並據此產出兩類藥物候選。

標的依據
--------
IRF4 是 DepMap 中對被套細胞淋巴瘤最具選擇性的依賴基因（gene effect −1.32，z = −3.19），
在我們的腫瘤殘差中也是上升的。它是轉錄因子、沒有傳統小分子口袋，
但有一個清楚的功能性表面：**DNA 結合域（殘基 21–129，色胺酸五聯重複）**。

用 PDB 7OGS（IRF4 DNA 結合域 + DNA，2.37 Å）定出實際接觸 DNA 的殘基，
這些殘基就是「阻斷型結合劑」要蓋住的熱點。

產出
----
1. `results/36_irf4_hotspots.txt`        —— RFdiffusion3 用的熱點殘基清單
2. `results/36_irf4_target.pdb`          —— 只留單條蛋白質鏈的標的結構
3. `results/36_decoy_oligos.csv`         —— 誘餌雙股寡核苷酸（直接吸附 IRF4）
"""
import os, sys, itertools
from collections import defaultdict
import numpy as np
import pandas as pd

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
ST = f"{BASE}/data/structures"
RES = f"{BASE}/results"
PDB = f"{ST}/7ogs.pdb"
CUT = 4.5          # Å
COMP = str.maketrans("ACGT", "TGCA")

AA3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H',
       'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
       'TYR':'Y','VAL':'V'}
DNA = {"DA", "DT", "DG", "DC"}


def log(*a):
    print("[36]", *a, flush=True)


def parse(pdb):
    prot, dna, lines = [], [], []
    for l in open(pdb):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        res = l[17:20].strip(); ch = l[21]; num = int(l[22:26])
        xyz = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
        if res in AA3:
            prot.append((ch, num, res, xyz)); lines.append(l)
        elif res in DNA:
            dna.append((ch, num, res, xyz))
    return prot, dna, lines


prot, dna, lines = parse(PDB)
log(f"蛋白質原子 {len(prot)}，DNA 原子 {len(dna)}")

# 只用 A 鏈（單體）與它實際接觸的 DNA（C/D 鏈）
pA = [p for p in prot if p[0] == "A"]
dCD = [d for d in dna if d[0] in ("C", "D")]
P = np.array([p[3] for p in pA]); D = np.array([d[3] for d in dCD])
log(f"A 鏈原子 {len(P)}，C/D 鏈 DNA 原子 {len(D)}")

# 兩兩距離（分塊）
contact = defaultdict(float)
for i0 in range(0, len(P), 2000):
    blk = P[i0:i0 + 2000]
    d = np.sqrt(((blk[:, None, :] - D[None, :, :]) ** 2).sum(-1))
    mn = d.min(axis=1)
    for j, v in enumerate(mn):
        ch, num, res, _ = pA[i0 + j]
        key = (num, res)
        if key not in contact or v < contact[key]:
            contact[key] = v

hot = sorted([(n, r, v) for (n, r), v in contact.items() if v <= CUT])
log(f"\n=== IRF4 DNA 結合域中直接接觸 DNA 的殘基（≤{CUT} Å，共 {len(hot)} 個）===")
s = " ".join(f"{AA3[r]}{n}" for n, r, v in hot)
log("  " + s)

# 取最核心的（≤3.5 Å）作為 RFdiffusion 熱點
core = [(n, r, v) for n, r, v in hot if v <= 3.5]
log(f"\n=== 核心熱點（≤3.5 Å，共 {len(core)} 個）===")
log("  " + " ".join(f"{AA3[r]}{n}" for n, r, v in core))
with open(f"{RES}/36_irf4_hotspots.txt", "w") as fh:
    fh.write("# IRF4 DNA 結合域的 DNA 接觸殘基（來源 PDB 7OGS A 鏈）\n")
    fh.write("# 格式: 殘基編號,單字母,最短距離(Å)\n")
    for n, r, v in hot:
        fh.write(f"{n},{AA3[r]},{v:.2f}\n")
    fh.write("\n# RFdiffusion 熱點字串（A 鏈）\n")
    fh.write("hotspots=[" + ",".join(f"'A{n}'" for n, r, v in core) + "]\n")

# 輸出單鏈標的結構
with open(f"{RES}/36_irf4_target.pdb", "w") as fh:
    for l in lines:
        if l[21] == "A":
            fh.write(l)
    fh.write("END\n")
log(f"\n寫出 {RES}/36_irf4_target.pdb（僅 A 鏈，{len(pA)} 原子）")

# ---------------- 誘餌雙股寡核苷酸 ----------------
# 7OGS 中 IRF4 結合的是 ISRE 類元件；AICE（與 BATF-JUN 共用）與 EICE（與 PU.1/SPIB 共用）
# 是 IRF4 在 B 細胞中的兩種主要複合位點。
rows = []
def add(name, top, note):
    bot = top.translate(COMP)[::-1]
    rows.append(dict(name=name, top_strand_5to3=top, bottom_strand_5to3=bot,
                     length=len(top), gc=round((top.count("G")+top.count("C"))/len(top), 2),
                     note=note))

add("ISRE_from_7OGS", "TGTACTTTCGGTTTCAGTTA",
    "直接取自 7OGS 共結晶的 DNA（IRF4 雙體實際結合的序列）")
add("ISRE_high_affinity_x2", "TTTCACTTTCACTTTCACTTTCAC",
    "四個串聯的 AANTGAAA 反向核心（TTTC 重複），提高吸附容量")
add("AICE1_BATF_IRF4", "TGAGTCATTTC",
    "AP-1 核心(TGAGTCA)+IRF 核心(TTTC)，BATF-IRF4 複合體位點；與本研究中被打開的 AP-1/BATF 染色質一致")
add("AICE2_BATF_IRF4", "TTTCTGAGTCA", "AICE 的另一種排列")
add("EICE_PU1_IRF4", "GGAAGTGAAAC",
    "ETS 核心(GGAA)+IRF 核心(GAAA)，PU.1/SPIB-IRF4 複合體位點")
add("AICE_x3_tandem", "TGAGTCATTTCTGAGTCATTTCTGAGTCATTTC",
    "三串聯 AICE，作為高親和力誘餌")
add("scrambled_control", "TGCATAGCATGCATAGCATG",
    "陰性對照：鹼基組成相近但無 IRF 核心")

df = pd.DataFrame(rows)
df.to_csv(f"{RES}/36_decoy_oligos.csv", index=False)
log(f"\n=== 誘餌雙股寡核苷酸候選（{len(df)} 條）===")
log("\n" + df[["name", "top_strand_5to3", "length", "gc", "note"]].to_string(index=False))
log("\n建議化學修飾：磷硫酸酯骨架（抗核酸酶）＋ 3' 端反向 dT 封端；")
log("            或做成髮夾式(hairpin)單股以提高穩定度。")
log("DONE")
