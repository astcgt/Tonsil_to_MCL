#!/usr/bin/env python
"""步驟 136（第 7 階的可行版本）：用**不是共摺疊信心分數**的方法，再判一次特異性。

為什麼不是 AF3 / Chai-1
-----------------------
第 7 階原本想用另一個共摺疊模型（AlphaFold3 或 Chai-1）重跑同一批設計。
AF3 的權重是申請制；Chai-1、ESMFold（需要 openfold + omegaconf）都要安裝新套件，
本機不允許安裝。所以改用兩個**不需要任何新套件**的正交軸：

A. 預測結構的物理性質 —— 埋藏表面積（Shrake-Rupley SASA）、界面接觸數、
   界面疏水比例。這些量完全不看模型的信心分數，只看它建出來的那個結構長什麼樣。
B. StaB-ddG（third_party，ICML 2025）—— 反向摺疊對數機率 + 熱力學循環，
   算式裡沒有任何共摺疊信心分數。
   用法很直接：**把標的從 MafG「突變成」MafK / MafF**（兩者長度相同、可一對一對齊，
   分別差 15 / 24 個殘基），看結合 ΔΔG 變多少。
   ΔΔG ≈ 0 → 沒有選擇性；ΔΔG 明顯為正 → 換成旁系同源後結合變差，也就是有選擇性。

輸出 results/136_*
"""
import os
import sys
import json
import glob
import shutil
import argparse
import subprocess
from collections import defaultdict

import numpy as np
import pandas as pd

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
RES = f"{BASE}/results"
AUX = f"{BASE}/data/aux/uniprot"
PRED = f"{RES}/132_specificity/boltz_results_inputs/predictions"
REPO = f"{BASE}/third_party/StaB-ddG"
PY = f"{BASE}/env/protdes/bin/python"
WORK = f"{RES}/136_orthogonal"
HYDRO = set("AVLIMFWCY")


def log(*a):
    print("[136]", *a, flush=True)


# --------------------------------------------------------- A. 界面幾何 ---
def interface_geometry(cif):
    from Bio.PDB import MMCIFParser
    from Bio.PDB.SASA import ShrakeRupley
    import warnings
    warnings.filterwarnings("ignore")
    st = MMCIFParser(QUIET=True).get_structure("x", cif)[0]
    chains = {c.id: c for c in st}
    if not {"A", "B"} <= set(chains):
        return None
    sr = ShrakeRupley()

    def sasa(entity):
        sr.compute(entity, level="R")
        return sum(r.sasa for r in entity.get_residues())

    import copy
    a_only = copy.deepcopy(chains["A"])
    b_only = copy.deepcopy(chains["B"])
    sa, sb = sasa(a_only), sasa(b_only)
    sr.compute(st, level="R")
    sab = sum(r.sasa for ch in ("A", "B") for r in chains[ch])
    bsa = sa + sb - sab

    # 界面殘基與接觸
    ca = {r.id[1]: r for r in chains["A"]}
    cb = {r.id[1]: r for r in chains["B"]}
    pa = np.array([[a.coord for a in r if a.element != "H"] for r in chains["A"]],
                  dtype=object)
    contacts, ifa, ifb = 0, set(), set()
    A = [(r.id[1], np.array([at.coord for at in r])) for r in chains["A"]]
    B = [(r.id[1], np.array([at.coord for at in r])) for r in chains["B"]]
    for ia, xa in A:
        for ib, xb in B:
            d = np.sqrt(((xa[:, None, :] - xb[None, :, :]) ** 2).sum(-1))
            n = int((d < 5.0).sum())
            if n:
                contacts += n
                ifa.add(ia); ifb.add(ib)
    hyd = sum(1 for i in ifa
              if ca[i].get_resname() in
              {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "CYS", "TYR"})
    return dict(bsa=float(bsa), contacts=contacts, n_iface_A=len(ifa),
                n_iface_B=len(ifb),
                hydrophobic_frac_A=hyd / max(len(ifa), 1))


def panel_a():
    rows = []
    cifs = sorted(glob.glob(f"{PRED}/**/*_model_*.cif", recursive=True))
    log("找到 %d 個預測結構" % len(cifs))
    for i, c in enumerate(cifs):
        nm = os.path.basename(c).rsplit("_model_", 1)[0]
        g = interface_geometry(c)
        if g:
            rows.append(dict(name=nm, file=os.path.basename(c), **g))
        if (i + 1) % 25 == 0:
            log("  %d/%d" % (i + 1, len(cifs)))
    D = pd.DataFrame(rows)
    D["target"] = D.name.str.split("__").str[0]
    D["binder"] = D.name.str.split("__").str[1]
    D.to_csv(f"{RES}/136_interface_geometry_raw.csv", index=False)
    R = D.groupby(["binder", "target"]).agg(
        n=("bsa", "size"), bsa_mean=("bsa", "mean"), bsa_sd=("bsa", "std"),
        contacts_mean=("contacts", "mean"),
        iface_res_A=("n_iface_A", "mean"),
        hydrophobic=("hydrophobic_frac_A", "mean")).reset_index()
    R.to_csv(f"{RES}/136_interface_geometry.csv", index=False)
    pd.set_option("display.width", 240)
    order = ["MAFG", "MAFF", "MAFK", "BACH1", "NFE2L2", "JUN", "FOS"]
    for m in ["bsa_mean", "contacts_mean"]:
        p = R.pivot(index="binder", columns="target", values=m)
        p = p[[c for c in order if c in p.columns]]
        log("\n=== %s（Å² 或接觸原子對數；列 = 結合劑） ===\n" % m, p.round(0).to_string())
    return R


# ------------------------------------------------------- B. StaB-ddG ---
def frag(gene):
    d = json.load(open(f"{AUX}/{gene}.json"))
    return d["sequence"]["value"][22:121]


def cif_to_pdb(cif, out):
    from Bio.PDB import MMCIFParser, PDBIO
    import warnings
    warnings.filterwarnings("ignore")
    st = MMCIFParser(QUIET=True).get_structure("x", cif)
    io = PDBIO(); io.set_structure(st); io.save(out)


def panel_b(gpu="4", trials=3, mc=20):
    mg, mf, mk = frag("MAFG"), frag("MAFF"), frag("MAFK")
    subs = {"to_MAFK": [(i + 1, mg[i], mk[i]) for i in range(99) if mg[i] != mk[i]],
            "to_MAFF": [(i + 1, mg[i], mf[i]) for i in range(99) if mg[i] != mf[i]]}
    log("MafG→MafK %d 個取代；MafG→MafF %d 個" % (len(subs["to_MAFK"]), len(subs["to_MAFF"])))

    cifs = sorted(glob.glob(f"{PRED}/**/MAFG__D1_best_model_*.cif", recursive=True))
    if not cifs:
        log("找不到 MAFG__D1_best 的預測結構"); return None
    os.makedirs(f"{WORK}/pdb", exist_ok=True)
    pdb = f"{WORK}/pdb/complex.pdb"
    cif_to_pdb(cifs[0], pdb)
    log("用結構:", os.path.basename(cifs[0]))

    # 檢查鏈 A 的序列真的是 MafG 片段
    from Bio.PDB import PDBParser
    from Bio.PDB.Polypeptide import protein_letters_3to1 as t31
    st = PDBParser(QUIET=True).get_structure("x", pdb)[0]
    seqA = "".join(t31.get(r.get_resname(), "X") for r in st["A"])
    if seqA != mg:
        log("⚠️ 鏈 A 序列與 MafG 片段不符:", seqA[:40], "…")
        log("   MafG 片段:", mg[:40], "…")
        return None
    log("鏈 A 序列與 MafG 片段相符（99 aa）")

    # 熱點正對照：七肽登記上的白胺酸改成丙胺酸
    hot = [p for p in range(1, 100) if mg[p - 1] == "L" and 47 <= p <= 85]
    entries = []
    for k, v in subs.items():
        entries.append((k, ",".join("%sA%d%s" % (w, p, m) for p, w, m in v)))
    for p, w, m in subs["to_MAFK"]:
        entries.append(("single_K_%s%d%s" % (w, p, m), "%sA%d%s" % (w, p, m)))
    for p in hot:
        entries.append(("ctrl_L%dA" % p, "LA%dA" % p))
    log("StaB-ddG 項目數: %d（含 %d 個白胺酸對照）" % (len(entries), len(hot)))

    csv = f"{WORK}/mutations.csv"
    # StaB-ddG 的格式：#Pdb = "<檔名>_<鏈組1>_<鏈組2>"，再一欄 mutation
    pd.DataFrame([{"#Pdb": "complex_A_B", "mutation": m} for _, m in entries]).to_csv(
        csv, index=False)
    pd.DataFrame([{"label": n, "mutation": m} for n, m in entries]).to_csv(
        f"{WORK}/mutation_labels.csv", index=False)

    out = f"{WORK}/stabddg_out"
    os.makedirs(out, exist_ok=True)
    env = dict(os.environ, CUDA_VISIBLE_DEVICES=str(gpu))
    cmd = [PY, "run_stabddg.py", "--pdb_dir", f"{WORK}/pdb", "--csv_path", csv,
           "--output_dir", out, "--trials", str(trials), "--mc_samples", str(mc),
           "--seed", "0", "--run_name", "spec"]
    log("跑 StaB-ddG …")
    p = subprocess.run(cmd, cwd=REPO, env=env, capture_output=True, text=True)
    log("return code:", p.returncode)
    if p.returncode != 0:
        log(p.stdout[-1200:]); log(p.stderr[-2500:])
        return None
    res = sorted(glob.glob(f"{out}/**/*.csv", recursive=True))
    log("輸出:", [os.path.basename(x) for x in res])
    if not res:
        return None
    D = pd.concat([pd.read_csv(x) for x in res], ignore_index=True)
    D.to_csv(f"{RES}/136_stabddg_specificity.csv", index=False)
    log("\n" + D.head(40).to_string())
    return D


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gpu", default="4")
    ap.add_argument("--skip_a", action="store_true")
    ap.add_argument("--skip_b", action="store_true")
    a = ap.parse_args()
    os.makedirs(WORK, exist_ok=True)
    out = {}
    if not a.skip_a:
        log("=== A. 界面幾何（正交軸一：結構的物理性質） ===")
        R = panel_a()
        out["geometry"] = json.loads(R.round(3).to_json(orient="records"))
    if not a.skip_b:
        log("\n=== B. StaB-ddG（正交軸二：結合 ΔΔG） ===")
        D = panel_b(gpu=a.gpu)
        if D is not None:
            out["stabddg"] = json.loads(D.round(4).to_json(orient="records"))
    json.dump(out, open(f"{RES}/136_orthogonal_specificity.json", "w"),
              indent=2, ensure_ascii=False)
    log("\nDONE → results/136_*")


if __name__ == "__main__":
    main()
