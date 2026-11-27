#!/usr/bin/env python
"""步驟 141：MafG 的三個夥伴各用到它哪幾個殘基？有沒有「只有 BACH1 用」的表位？

治療邏輯
--------
MafG 不單獨工作。它至少有三種複合體：
  MafG–NRF2   抗氧化反應，正常細胞需要      → 想保留
  MafG–BACH1  抑制型；腫瘤裡 BACH1 富集 1.38× → 想擋掉
  MafG–MafG   同源二聚，抑制型
現在的設計品瞄準 DNA 結合鹼性區 + 白胺酸拉鍊，那是**三者共用**的面，所以會全擋。
如果 BACH1 用到一些 NRF2 不用的 MafG 殘基，那就存在一個「BACH1 選擇性表位」，
擋住它可以只關掉腫瘤那一半，留下正常細胞的抗氧化反應。

這正是 BCL-2 vs BCL-xL（venetoclax）走過的路：同一個摺疊家族，
靠界面上少數不同的殘基做出選擇性。

做法
----
1. 用 Boltz-2 預測三個複合體（都用 bZIP 同源片段，長度可比）：
   MAFG+NFE2L2、MAFG+BACH1、MAFG+MAFG，各 5 次取樣。
2. 對每個複合體算 MafG 鏈的**逐殘基埋藏表面積**（Shrake-Rupley）。
3. 定義：
   BACH1 專屬表位 = BACH1 埋得多（ΔSASA > 10 Å²）但 NRF2 埋得少（< 5 Å²）的位置
   共用核心     = 三個都埋
4. 把設計品 D1 實際接觸的殘基疊上去（取自步驟 132 的預測），
   看它現在蓋的是共用核心還是 BACH1 專屬表位。

輸出 results/141_*
"""
import os
import sys
import json
import glob
import copy
import argparse
import subprocess
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
RES = f"{BASE}/results"
AUX = f"{BASE}/data/aux/uniprot"
OUT = f"{RES}/141_partners"
D132 = f"{RES}/132_specificity/boltz_results_inputs/predictions"
UP, DOWN = 28, 7
NRF2 = ("HLTRDELRAKALHIPFPVEKIINLPVVDFNEMMSKEQFNEAQLALIRDIRRRRGKNKVAAQNCRKRK"
        "LENNIVELEQDLDHLKDEKEKLLKEKGENDKSLHLLKKQLS")


def log(*a):
    print("[141]", *a, flush=True)


def frag(gene):
    d = json.load(open(f"{AUX}/{gene}.json"))
    s = d["sequence"]["value"]
    dom = [f for f in d.get("features", [])
           if f["type"] == "Domain" and f.get("description") == "bZIP"][0]
    a = dom["location"]["start"]["value"]; b = dom["location"]["end"]["value"]
    return s[max(0, a - UP - 1):min(len(s), b + DOWN)]


def per_residue_bsa(cif, chain="A"):
    """鏈 chain 上每個殘基，從單獨存在到形成複合體，埋掉多少 SASA。"""
    from Bio.PDB import MMCIFParser
    from Bio.PDB.SASA import ShrakeRupley
    st = MMCIFParser(QUIET=True).get_structure("x", cif)[0]
    if chain not in [c.id for c in st]:
        return None
    sr = ShrakeRupley()
    alone = copy.deepcopy(st[chain])
    sr.compute(alone, level="R")
    free = {r.id[1]: r.sasa for r in alone}
    sr.compute(st, level="R")
    bound = {r.id[1]: r.sasa for r in st[chain]}
    return pd.Series({i: free[i] - bound.get(i, 0.0) for i in free}).sort_index()


def main():
    global OUT
    ap = argparse.ArgumentParser()
    ap.add_argument("--gpu", default="4")
    ap.add_argument("--samples", type=int, default=5)
    ap.add_argument("--tag", default="",
                    help="輸出加後綴（例如 _n20）。提高取樣數重跑時用，\n                         以免覆蓋已發表的 n=5 結果與 MD 的起始結構。")
    ap.add_argument("--score_only", action="store_true")
    a = ap.parse_args()
    OUT = f"{RES}/141_partners{a.tag}"
    os.makedirs(OUT, exist_ok=True)
    yml = f"{OUT}/inputs"; os.makedirs(yml, exist_ok=True)

    MAFG = frag("MAFG")
    BACH1 = frag("BACH1")
    log("MAFG %d aa | BACH1 bZIP %d aa | NRF2 %d aa" % (len(MAFG), len(BACH1), len(NRF2)))

    pairs = {"MAFG_x_NRF2": NRF2, "MAFG_x_BACH1": BACH1, "MAFG_x_MAFG": MAFG}
    if not a.score_only:
        for f in glob.glob(f"{yml}/*.yaml"):
            os.remove(f)
        for nm, partner in pairs.items():
            with open(f"{yml}/{nm}.yaml", "w") as fh:
                fh.write("version: 1\nsequences:\n")
                fh.write(f"  - protein:\n      id: A\n      sequence: {MAFG}\n      msa: empty\n")
                fh.write(f"  - protein:\n      id: B\n      sequence: {partner}\n      msa: empty\n")
        env = dict(os.environ, CUDA_VISIBLE_DEVICES=a.gpu)
        cmd = [f"{BASE}/env/boltz/bin/python", "-c",
               "import sys; from boltz.main import cli; sys.argv[0]='boltz'; cli()",
               "predict", yml, "--out_dir", OUT, "--cache", f"{BASE}/models/boltz",
               "--output_format", "mmcif", "--diffusion_samples", str(a.samples),
               "--recycling_steps", "3", "--devices", "1", "--accelerator", "gpu",
               "--no_kernels"]
        log("跑 Boltz-2（3 個複合體 × %d 次取樣）…" % a.samples)
        p = subprocess.run(cmd, env=env, capture_output=True, text=True)
        log("return code:", p.returncode)
        if p.returncode != 0:
            log(p.stderr[-2500:]); return

    # ---------- 逐殘基埋藏 ----------
    prof = {}
    for nm in pairs:
        cifs = sorted(glob.glob(f"{OUT}/**/{nm}_model_*.cif", recursive=True))
        if not cifs:
            log("缺少", nm); continue
        S = [per_residue_bsa(c) for c in cifs]
        S = [s for s in S if s is not None]
        prof[nm] = pd.concat(S, axis=1).mean(1)
        log("%s：%d 個模型，MafG 側總埋藏 %.0f Å²"
            % (nm, len(S), prof[nm].sum()))

    # 設計品 D1 在步驟 132 的預測裡用到哪些 MafG 殘基
    cifs = sorted(glob.glob(f"{D132}/**/MAFG__D1_best_model_*.cif", recursive=True))
    if cifs:
        S = [per_residue_bsa(c) for c in cifs]
        prof["MAFG_x_D1design"] = pd.concat([s for s in S if s is not None],
                                            axis=1).mean(1)
        log("MAFG_x_D1design：%d 個模型，MafG 側總埋藏 %.0f Å²"
            % (len(S), prof["MAFG_x_D1design"].sum()))

    P = pd.DataFrame(prof).fillna(0.0)
    P.index.name = "mafg_pos_frag"
    P["uniprot_pos"] = P.index + 22
    seq = frag("MAFG")
    P["aa"] = [seq[i - 1] for i in P.index]
    P.round(2).to_csv(f"{RES}/141_per_residue_bsa{a.tag}.csv")

    nrf2, bach1 = P["MAFG_x_NRF2"], P["MAFG_x_BACH1"]
    homo = P.get("MAFG_x_MAFG", pd.Series(0.0, index=P.index))
    d1 = P.get("MAFG_x_D1design", pd.Series(0.0, index=P.index))

    bach_only = P[(bach1 > 10) & (nrf2 < 5)]
    nrf2_only = P[(nrf2 > 10) & (bach1 < 5)]
    shared = P[(bach1 > 10) & (nrf2 > 10)]
    log("\n=== BACH1 專屬表位（BACH1 埋 >10 Å²、NRF2 埋 <5 Å²）===")
    log("共 %d 個位置" % len(bach_only))
    if len(bach_only):
        log("\n" + bach_only[["aa", "uniprot_pos", "MAFG_x_BACH1", "MAFG_x_NRF2",
                              "MAFG_x_MAFG", "MAFG_x_D1design"]].round(1).to_string())
    log("\n=== NRF2 專屬（想保留的那一面）：%d 個位置 ===" % len(nrf2_only))
    if len(nrf2_only):
        log("\n" + nrf2_only[["aa", "uniprot_pos", "MAFG_x_NRF2", "MAFG_x_BACH1"]]
            .round(1).to_string())
    log("\n=== 共用核心：%d 個位置 ===" % len(shared))

    # 設計品現在蓋的是哪一類
    d1_cov = P[d1 > 10]
    log("\n=== 設計品 D1 埋住的 MafG 位置：%d 個 ===" % len(d1_cov))
    if len(d1_cov):
        n_sh = int(((d1 > 10) & (bach1 > 10) & (nrf2 > 10)).sum())
        n_bo = int(((d1 > 10) & (bach1 > 10) & (nrf2 < 5)).sum())
        n_no = int(((d1 > 10) & (nrf2 > 10) & (bach1 < 5)).sum())
        log("  其中落在共用核心 %d 個、BACH1 專屬 %d 個、NRF2 專屬 %d 個"
            % (n_sh, n_bo, n_no))
        log("  → 若「共用核心」佔絕大多數，代表目前的設計無法區分兩種複合體。")

    corr = P[[c for c in P.columns if c.startswith("MAFG_x_")]].corr(method="spearman")
    log("\n=== 四種夥伴的逐殘基埋藏輪廓相關（越高＝用的面越像） ===\n",
        corr.round(3).to_string())

    json.dump({"bach1_only_positions": [
                   {"frag": int(i), "uniprot": int(r.uniprot_pos), "aa": r.aa,
                    "bsa_bach1": float(r.MAFG_x_BACH1), "bsa_nrf2": float(r.MAFG_x_NRF2),
                    "bsa_d1": float(r.get("MAFG_x_D1design", 0))}
                   for i, r in bach_only.iterrows()],
               "nrf2_only_positions": [
                   {"frag": int(i), "uniprot": int(r.uniprot_pos), "aa": r.aa,
                    "bsa_nrf2": float(r.MAFG_x_NRF2), "bsa_bach1": float(r.MAFG_x_BACH1)}
                   for i, r in nrf2_only.iterrows()],
               "n_shared_core": int(len(shared)),
               "profile_correlation": json.loads(corr.round(4).to_json()),
               "total_bsa": {k: float(v.sum()) for k, v in prof.items()}},
              open(f"{RES}/141_partner_interfaces{a.tag}.json", "w"),
              indent=2, ensure_ascii=False)
    log("\nDONE → results/141_*")


if __name__ == "__main__":
    main()
