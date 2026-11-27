#!/usr/bin/env python
"""步驟 143：「BACH1 專屬表位」是真的嗎？用 StaB-ddG 正交檢定。

步驟 141 找到 MafG 上 9 個位置（UniProt 61,64,65,68,71,115,117,118,121），
BACH1 預測會埋住它們（10–70 Å²），而 NRF2 不會。
其中 NRF2 那一半是**實驗事實**：晶體 7X5E 上這 9 個位置的埋藏全部是 0.0 Å²。
BACH1 那一半是預測（沒有 MafG–BACH1 的實驗結構）。

檢定
----
把這 9 個位置全部換成丙胺酸，分別算在兩個複合體上的結合 ΔΔG：
  (a) 晶體 MafG–NRF2（7X5E）     → 預期 ΔΔG ≈ 0（NRF2 根本沒碰這些位置）
  (b) 預測 MafG–BACH1            → 若表位是真的，ΔΔG 應明顯為正
對照：同樣打掉 9 個**共用核心**位置，兩個複合體都應該明顯為正。

這個設計自己帶陰性對照與陽性對照，所以結果可解讀。

輸出 results/143_*
"""
import os
import sys
import ast
import json
import glob
import shutil
import argparse
import subprocess
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
RES = f"{BASE}/results"
REPO = f"{BASE}/third_party/StaB-ddG"
PY = f"{BASE}/env/protdes/bin/python"
WORK = f"{RES}/143_epitope"
XTAL = f"{RES}/102_stabddg/pdbs_fwd/mafg_nrf2_AB.pdb"
P141 = f"{RES}/141_partners"

BACH1_ONLY = [39, 42, 43, 46, 49, 93, 95, 96, 99]      # 片段編號（UniProt = +22）


def log(*a):
    print("[143]", *a, flush=True)


def cif_to_pdb(cif, out):
    from Bio.PDB import MMCIFParser, PDBIO
    st = MMCIFParser(QUIET=True).get_structure("x", cif)
    io = PDBIO(); io.set_structure(st); io.save(out)


def chainA_seq(pdb):
    from Bio.PDB import PDBParser
    from Bio.PDB.Polypeptide import protein_letters_3to1 as t31
    st = PDBParser(QUIET=True).get_structure("x", pdb)[0]
    return {r.id[1]: t31.get(r.get_resname(), "X") for r in st["A"]}


def run(pdb_dir, csv, out, gpu, trials=3, mc=20, tag="run"):
    env = dict(os.environ, CUDA_VISIBLE_DEVICES=str(gpu))
    cmd = [PY, "run_stabddg.py", "--pdb_dir", pdb_dir, "--csv_path", csv,
           "--output_dir", out, "--trials", str(trials), "--mc_samples", str(mc),
           "--seed", "0", "--run_name", tag]
    p = subprocess.run(cmd, cwd=REPO, env=env, capture_output=True, text=True)
    if p.returncode != 0:
        log("StaB-ddG 失敗:", p.stderr[-1500:]); return None
    f = f"{out}/{tag}.csv"
    if not os.path.exists(f):
        log("找不到", f); return None
    D = pd.read_csv(f, index_col=0)
    D["mut"] = D["Mutation"].apply(lambda s: ",".join(ast.literal_eval(s)))
    pc = [c for c in D.columns if c.startswith("pred_")]
    D["ddg"] = D[pc].mean(1); D["sd"] = D[pc].std(1)
    return D[["mut", "ddg", "sd"]]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gpu", default="4")
    a = ap.parse_args()
    os.makedirs(WORK, exist_ok=True)

    P = pd.read_csv(f"{RES}/141_per_residue_bsa.csv", index_col=0)
    shared = P[(P.MAFG_x_BACH1 > 10) & (P.MAFG_x_NRF2 > 10)]
    shared9 = shared.nlargest(9, "MAFG_x_NRF2")
    log("BACH1 專屬表位 9 個位置 (UniProt):", [int(i) + 22 for i in BACH1_ONLY])
    log("共用核心對照 9 個位置 (UniProt):", [int(i) + 22 for i in shared9.index])

    targets = {}
    # (a) 晶體 MafG–NRF2
    os.makedirs(f"{WORK}/xtal", exist_ok=True)
    shutil.copy(XTAL, f"{WORK}/xtal/mafgnrf2.pdb")
    targets["NRF2_crystal"] = (f"{WORK}/xtal", "mafgnrf2_A_B", f"{WORK}/xtal/mafgnrf2.pdb")
    # (b) 預測 MafG–BACH1
    cifs = sorted(glob.glob(f"{P141}/**/MAFG_x_BACH1_model_*.cif", recursive=True))
    if not cifs:
        log("找不到 MafG–BACH1 預測結構"); return
    os.makedirs(f"{WORK}/bach1", exist_ok=True)
    cif_to_pdb(cifs[0], f"{WORK}/bach1/mafgbach1.pdb")
    targets["BACH1_predicted"] = (f"{WORK}/bach1", "mafgbach1_A_B",
                                  f"{WORK}/bach1/mafgbach1.pdb")
    # (c) 預測 MafG–MafG（同源二聚，也是腫瘤側）
    cifs2 = sorted(glob.glob(f"{P141}/**/MAFG_x_MAFG_model_*.cif", recursive=True))
    if cifs2:
        os.makedirs(f"{WORK}/homo", exist_ok=True)
        cif_to_pdb(cifs2[0], f"{WORK}/homo/mafghomo.pdb")
        targets["MAFG_homodimer_predicted"] = (f"{WORK}/homo", "mafghomo_A_B",
                                               f"{WORK}/homo/mafghomo.pdb")

    out = {}
    for name, (pdir, pid, pdb) in targets.items():
        seqA = chainA_seq(pdb)
        def muts(idx):
            m = []
            for i in idx:
                i = int(i)
                if i in seqA and seqA[i] != "A":
                    m.append("%sA%dA" % (seqA[i], i))
            return m
        mb, ms = muts(BACH1_ONLY), muts(shared9.index)
        entries = [("bach1_only_9A", ",".join(mb)),
                   ("shared_core_9A", ",".join(ms))]
        entries += [("single_bo_%s" % x, x) for x in mb]
        csv = f"{WORK}/{name}.csv"
        pd.DataFrame([{"#Pdb": pid, "mutation": m} for _, m in entries]).to_csv(
            csv, index=False)
        lab = {m: n for n, m in entries}
        log("\n跑 %s（鏈 A 解析 %d 個殘基；表位 %d 個、共用核心 %d 個可突變）"
            % (name, len(seqA), len(mb), len(ms)))
        D = run(pdir, csv, f"{WORK}/out_{name}", a.gpu, tag=name)
        if D is None:
            continue
        D["label"] = D["mut"].map(lab)
        D.to_csv(f"{RES}/143_{name}.csv", index=False)
        row = {}
        for k in ("bach1_only_9A", "shared_core_9A"):
            r = D[D.label == k]
            row[k] = {"ddg": float(r.ddg.iloc[0]), "sd": float(r.sd.iloc[0])} if len(r) else None
        row["singles"] = {r.label: float(r.ddg) for r in
                          D[D.label.astype(str).str.startswith("single_")].itertuples()}
        out[name] = row
        log("  表位 9 殘基全打掉 : ΔΔG = %s" %
            ("%.2f ± %.2f" % (row['bach1_only_9A']['ddg'], row['bach1_only_9A']['sd'])
             if row.get('bach1_only_9A') else "-"))
        log("  共用核心 9 殘基     : ΔΔG = %s" %
            ("%.2f ± %.2f" % (row['shared_core_9A']['ddg'], row['shared_core_9A']['sd'])
             if row.get('shared_core_9A') else "-"))

    log("\n=== 判讀表 ===")
    log("%-26s %18s %18s" % ("複合體", "打掉 BACH1 專屬 9 個", "打掉共用核心 9 個"))
    for k, v in out.items():
        log("%-26s %18s %18s"
            % (k,
               "%.2f" % v["bach1_only_9A"]["ddg"] if v.get("bach1_only_9A") else "-",
               "%.2f" % v["shared_core_9A"]["ddg"] if v.get("shared_core_9A") else "-"))
    log("\n預期：NRF2_crystal 的第一欄應接近 0（晶體上這 9 個位置埋藏就是 0），")
    log("      BACH1_predicted 的第一欄應明顯為正；兩者的第二欄都應明顯為正。")

    json.dump({"bach1_only_positions_uniprot": [int(i) + 22 for i in BACH1_ONLY],
               "shared_core_positions_uniprot": [int(i) + 22 for i in shared9.index],
               "results": out},
              open(f"{RES}/143_epitope_test.json", "w"), indent=2, ensure_ascii=False)
    log("\nDONE → results/143_*")


if __name__ == "__main__":
    main()
