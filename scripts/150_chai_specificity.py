#!/usr/bin/env python
"""步驟 150（第 7 階，真正的版本）：用 Chai-1 重跑整個特異性盤。

為什麼需要
----------
步驟 132（Boltz-2）說設計品對 MafG / MafF / MafK 一視同仁（ipSAE 0.633 / 0.593 / 0.653）。
步驟 136–137（StaB-ddG）說有選擇性，而且比天然 NRF2 更挑（ΔΔG +7.3 vs +0.95）。
兩個方法矛盾。Chai-1 是**第三個、架構與訓練資料都不同**的共摺疊模型，
用它當裁判。

設計與步驟 132 完全一致：7 個標的 × 4 個結合劑 × 5 次 diffusion 取樣。
標的片段、結合劑序列都直接沿用，所以兩個模型的結果可以逐格比較。

輸出 results/150_*
"""
import os
import sys
import json
import glob
import random
import argparse
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
os.environ.setdefault("DISABLE_PANDERA_IMPORT_WARNING", "True")

import numpy as np
import pandas as pd

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
RES = f"{BASE}/results"
AUX = f"{BASE}/data/aux/uniprot"
OUT = f"{RES}/150_chai"
UP, DOWN = 28, 7
TARGETS = ["MAFG", "MAFF", "MAFK", "BACH1", "NFE2L2", "JUN", "FOS"]
NRF2 = ("HLTRDELRAKALHIPFPVEKIINLPVVDFNEMMSKEQFNEAQLALIRDIRRRRGKNKVAAQNCRKRK"
        "LENNIVELEQDLDHLKDEKEKLLKEKGENDKSLHLLKKQLS")


def log(*a):
    print("[150]", *a, flush=True)


def frag(gene):
    d = json.load(open(f"{AUX}/{gene}.json"))
    s = d["sequence"]["value"]
    dom = [f for f in d.get("features", [])
           if f["type"] == "Domain" and f.get("description") == "bZIP"][0]
    a = dom["location"]["start"]["value"]; b = dom["location"]["end"]["value"]
    return s[max(0, a - UP - 1):min(len(s), b + DOWN)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gpu", default="4")
    ap.add_argument("--samples", type=int, default=5)
    ap.add_argument("--tag", default="",
                    help="輸出加後綴（例如 _n20）。提高取樣數重跑時用，\n                         以免覆蓋已發表的 n=5 結果與 MD 的起始結構。")
    ap.add_argument("--score_only", action="store_true")
    a = ap.parse_args()
    global OUT
    OUT = f"{RES}/150_chai{a.tag}"
    os.environ["CUDA_VISIBLE_DEVICES"] = a.gpu
    os.makedirs(OUT, exist_ok=True)

    tseq = {g: frag(g) for g in TARGETS}
    meta87 = json.load(open(f"{RES}/87_meta.json"))
    D1 = meta87["R1_mafg_dna_block_90_2_model_1_m4"]["seq"]
    D2 = meta87["R1_mafg_dna_block_90_1_model_0_m1"]["seq"]
    rng = random.Random(1)
    sc = list(D1); rng.shuffle(sc)
    BINDERS = {"D1_best": D1, "D2_second": D2, "CTRLPOS_NRF2": NRF2,
               "CTRLNEG_scram": "".join(sc)}
    combos = {f"{t}__{b}": (tseq[t], bs) for t in TARGETS for b, bs in BINDERS.items()}
    log("%d 個複合體 × %d 次取樣（與步驟 132 完全相同的輸入）" % (len(combos), a.samples))

    if not a.score_only:
        from chai_lab.chai1 import run_inference
        import torch
        for i, (nm, (ts, bs)) in enumerate(sorted(combos.items()), 1):
            od = Path(f"{OUT}/{nm}")
            if (od / "scores.model_idx_0.npz").exists():
                log("[%2d/%d] %s 已有結果，跳過" % (i, len(combos), nm)); continue
            od.mkdir(parents=True, exist_ok=True)
            fa = Path(f"{OUT}/{nm}.fasta")
            fa.write_text(f">protein|name=A\n{ts}\n>protein|name=B\n{bs}\n")
            log("[%2d/%d] %s …" % (i, len(combos), nm))
            try:
                run_inference(fasta_file=fa, output_dir=od,
                              num_trunk_recycles=3, num_diffn_timesteps=200,
                              num_diffn_samples=a.samples, seed=0,
                              use_esm_embeddings=True, use_msa_server=False,
                              device="cuda:0", low_memory=True)
            except Exception as e:
                log("  失敗:", type(e).__name__, str(e)[:300])
            torch.cuda.empty_cache()

    # ---------------- 收分數 ----------------
    rows = []
    for nm in combos:
        for f in sorted(glob.glob(f"{OUT}/{nm}/scores.model_idx_*.npz")):
            z = np.load(f)
            k = {x: z[x] for x in z.files}
            def sc1(name):
                v = k.get(name)
                return float(np.asarray(v).ravel()[0]) if v is not None else np.nan
            rows.append(dict(name=nm,
                             model=os.path.basename(f).split("idx_")[1].split(".")[0],
                             iptm=sc1("iptm"), ptm=sc1("ptm"),
                             aggregate_score=sc1("aggregate_score"),
                             chain_chain_clashes=sc1("total_clashes"),
                             has_inter_chain_clashes=sc1("has_inter_chain_clashes")))
    if not rows:
        log("沒有任何結果"); return
    D = pd.DataFrame(rows)
    D["target"] = D.name.str.split("__").str[0]
    D["binder"] = D.name.str.split("__").str[1]
    D.to_csv(f"{RES}/150_chai_raw{a.tag}.csv", index=False)

    R = D.groupby(["binder", "target"]).agg(
        n=("iptm", "size"), iptm_mean=("iptm", "mean"), iptm_sd=("iptm", "std"),
        agg_mean=("aggregate_score", "mean"), ptm_mean=("ptm", "mean")).reset_index()
    R.to_csv(f"{RES}/150_chai_summary{a.tag}.csv", index=False)
    pd.set_option("display.width", 240)
    order = [t for t in TARGETS]
    for m in ["iptm_mean", "agg_mean"]:
        p = R.pivot(index="binder", columns="target", values=m)
        p = p[[c for c in order if c in p.columns]]
        log("\n=== Chai-1 %s（列 = 結合劑） ===\n" % m, p.round(3).to_string())

    # 與 Boltz-2 並排
    try:
        B = pd.read_csv(f"{RES}/132_specificity_summary.csv")
        bp = B.pivot(index="binder", columns="target", values="iptm_mean")
        cp = R.pivot(index="binder", columns="target", values="iptm_mean")
        common = [c for c in order if c in bp.columns and c in cp.columns]
        log("\n=== 小 Maf 家族內的選擇性：兩個模型並排（ipTM） ===")
        log("%-15s %22s %22s" % ("結合劑", "Boltz-2 (MAFG/MAFF/MAFK)",
                                 "Chai-1 (MAFG/MAFF/MAFK)"))
        sm = [c for c in ("MAFG", "MAFF", "MAFK") if c in common]
        for b in cp.index:
            if b not in bp.index:
                continue
            log("%-15s %22s %22s"
                % (b, " / ".join("%.3f" % bp.loc[b, c] for c in sm),
                   " / ".join("%.3f" % cp.loc[b, c] for c in sm)))
        from scipy.stats import spearmanr
        x = bp[common].stack(); y = cp[common].stack()
        ix = x.index.intersection(y.index)
        rr = spearmanr(x[ix], y[ix])
        log("\n兩個模型在 %d 格上的 Spearman = %.3f (p = %.3g)"
            % (len(ix), rr.statistic, rr.pvalue))
    except Exception as e:
        log("並排比較失敗:", e)

    json.dump({"summary": json.loads(R.round(4).to_json(orient="records"))},
              open(f"{RES}/150_chai_specificity{a.tag}.json", "w"), indent=2, ensure_ascii=False)
    log("\nDONE → results/150_*")


if __name__ == "__main__":
    main()
