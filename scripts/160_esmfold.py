#!/usr/bin/env python
"""步驟 160：ESMFold —— 做兩件 Boltz-2 與 Chai-1 做不到的事，再加一個第三意見。

ESMFold 的歸納偏差和前兩者不同：**單序列、無 MSA、以蛋白質語言模型為主幹**，
而且它是在**單體**上訓練的。這讓它特別適合回答兩個共摺疊模型答不了的問題。

A. 設計品自己摺得起來嗎？（共摺疊分數答不了的問題）
   一個好的迷你結合劑應該**在沒有標的的時候就是摺好的**。
   如果它只有在結合時才成形，entropy 代價高、容易聚集、難表現。
   這是 de novo binder 領域的標準檢查項，而共摺疊分數無法區分
   「自己會摺」和「被標的帶著摺」。

B. 標的自己摺得起來嗎？
   步驟 132 裡 D1 對 NFE2L2 的 ipSAE 只有 0.008，但那一格的 target pLDDT 只有 64.6。
   所以「不結合」和「模型不會摺這個標的」分不開。ESMFold 可以單獨回答。

C. 第三個共摺疊意見（用 glycine linker 串成單鏈的標準作法）。
   ⚠️ ESMFold 是單體模型，用 linker 做複合體是公認較弱的用法，
      所以 C 的權重低於 A 和 B，只當佐證。

輸出 results/160_*
"""
import os
import sys
import json
import argparse
import warnings

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import torch

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
RES = f"{BASE}/results"
AUX = f"{BASE}/data/aux/uniprot"
UP, DOWN = 28, 7
TARGETS = ["MAFG", "MAFF", "MAFK", "BACH1", "NFE2L2", "JUN", "FOS"]
NRF2 = ("HLTRDELRAKALHIPFPVEKIINLPVVDFNEMMSKEQFNEAQLALIRDIRRRRGKNKVAAQNCRKRK"
        "LENNIVELEQDLDHLKDEKEKLLKEKGENDKSLHLLKKQLS")


def log(*a):
    print("[160]", *a, flush=True)


def frag(gene):
    d = json.load(open(f"{AUX}/{gene}.json"))
    s = d["sequence"]["value"]
    dom = [f for f in d.get("features", [])
           if f["type"] == "Domain" and f.get("description") == "bZIP"][0]
    a = dom["location"]["start"]["value"]; b = dom["location"]["end"]["value"]
    return s[max(0, a - UP - 1):min(len(s), b + DOWN)]


def calc_d0(L):
    L = np.maximum(26.0, np.asarray(L, float))
    return np.maximum(1.0, 1.24 * (L - 15.0) ** (1.0 / 3.0) - 1.8)


def ipsae(pae, ch, a, b, cut=10.0):
    """與步驟 93 同一條公式（Dunbrack ipSAE），回傳 min(A→B, B→A)。"""
    out = {}
    for c1, c2 in ((a, b), (b, a)):
        m = np.outer(ch == c1, ch == c2) & (pae < cut)
        n0 = m.sum(1)
        d0 = calc_d0(n0)
        rows = np.where((ch == c1) & (n0 > 0))[0]
        best = 0.0
        for i in rows:
            v = 1.0 / (1.0 + (pae[i, m[i]] / d0[i]) ** 2)
            best = max(best, float(v.mean()))
        out[(c1, c2)] = best
    return min(out[(a, b)], out[(b, a)])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gpu", default="3")
    ap.add_argument("--recycles", type=int, default=4)
    a = ap.parse_args()
    os.environ["CUDA_VISIBLE_DEVICES"] = a.gpu

    import esm
    log("載入 ESMFold（首次執行會下載約 2.6 GB 權重）…")
    model = esm.pretrained.esmfold_v1().eval().cuda()
    model.set_chunk_size(128)
    log("模型就緒")

    tseq = {g: frag(g) for g in TARGETS}
    meta = json.load(open(f"{RES}/87_meta.json"))
    import random
    D1 = meta["R1_mafg_dna_block_90_2_model_1_m4"]["seq"]
    D2 = meta["R1_mafg_dna_block_90_1_model_0_m1"]["seq"]
    rng = random.Random(1); sc = list(D1); rng.shuffle(sc)
    BINDERS = {"D1_best": D1, "D2_second": D2, "CTRLPOS_NRF2": NRF2,
               "CTRLNEG_scram": "".join(sc)}
    # 另外兩個設計品，看「自己會不會摺」有沒有區分力
    for k in ("R2_mafg_zeroescape_90_0_model_0_m4",
              "R1_mafg_dimer_block_70_2_model_3_m1",
              "R1_mafg_dna_block_90_0_model_3_m5"):
        if k in meta and "seq" in meta[k]:
            BINDERS[k] = meta[k]["seq"]

    # ---------------- A + B：單體 ----------------
    mono = {**{f"target_{k}": v for k, v in tseq.items()},
            **{f"binder_{k}": v for k, v in BINDERS.items()}}
    rows = []
    with torch.no_grad():
        for nm, s in mono.items():
            o = model.infer([s], num_recycles=a.recycles)
            pl = float(o["mean_plddt"][0]) / (100.0 if float(o["mean_plddt"][0]) > 1.5 else 1.0)
            ptm = float(o["ptm"][0])
            plr = o["plddt"][0].mean(-1).cpu().numpy()
            plr = plr / (100.0 if plr.max() > 1.5 else 1.0)
            rows.append(dict(name=nm, kind=nm.split("_")[0], length=len(s),
                             mean_plddt=pl, ptm=ptm,
                             frac_res_plddt_gt70=float((plr > 0.70).mean()),
                             frac_res_plddt_gt90=float((plr > 0.90).mean())))
            log("  %-34s len %3d  pLDDT %.3f  pTM %.3f  (>0.7 的殘基 %.0f%%)"
                % (nm, len(s), pl, ptm, 100 * rows[-1]["frac_res_plddt_gt70"]))
            torch.cuda.empty_cache()
    M = pd.DataFrame(rows)
    M.to_csv(f"{RES}/160_esmfold_monomer.csv", index=False)

    log("\n=== A. 設計品自己摺得起來嗎 ===")
    b = M[M.kind == "binder"].sort_values("mean_plddt", ascending=False)
    log("\n" + b[["name", "length", "mean_plddt", "ptm",
                  "frac_res_plddt_gt70"]].round(3).to_string(index=False))
    log("\n=== B. 標的自己摺得起來嗎（解釋步驟 132 的低分格） ===")
    t = M[M.kind == "target"].sort_values("mean_plddt", ascending=False)
    log("\n" + t[["name", "length", "mean_plddt", "ptm"]].round(3).to_string(index=False))

    # ---------------- C：複合體（第三意見） ----------------
    log("\n=== C. 複合體（glycine linker；ESMFold 的弱項，只當佐證） ===")
    crows = []
    core = ["D1_best", "D2_second", "CTRLPOS_NRF2", "CTRLNEG_scram"]
    with torch.no_grad():
        for tg in TARGETS:
            for bd in core:
                s = tseq[tg] + ":" + BINDERS[bd]
                o = model.infer([s], num_recycles=a.recycles)
                pae = o["predicted_aligned_error"][0].cpu().numpy()
                ci = o["chain_index"][0].cpu().numpy()
                pl = o["plddt"][0].mean(-1).cpu().numpy()
                pl = pl / (100.0 if pl.max() > 1.5 else 1.0)
                ch = np.where(ci == 0, "A", "B")
                cross = np.outer(ch == "A", ch == "B") | np.outer(ch == "B", ch == "A")
                crows.append(dict(
                    target=tg, binder=bd,
                    ptm=float(o["ptm"][0]),
                    ipsae_min=float(ipsae(pae, ch, "A", "B")),
                    ipae=float(pae[cross].mean()),
                    target_plddt=float(pl[ch == "A"].mean()),
                    binder_plddt=float(pl[ch == "B"].mean())))
                log("  %-7s x %-14s ipSAE %.3f  i-PAE %5.2f  pTM %.3f"
                    % (tg, bd, crows[-1]["ipsae_min"], crows[-1]["ipae"],
                       crows[-1]["ptm"]))
                torch.cuda.empty_cache()
    C = pd.DataFrame(crows)
    C.to_csv(f"{RES}/160_esmfold_complex.csv", index=False)
    pd.set_option("display.width", 240)
    for m in ("ipsae_min", "ipae"):
        p = C.pivot(index="binder", columns="target", values=m)[TARGETS]
        log("\n=== ESMFold %s ===\n" % m, p.round(3).to_string())

    out = {"monomer": json.loads(M.round(4).to_json(orient="records")),
           "complex": json.loads(C.round(4).to_json(orient="records"))}
    json.dump(out, open(f"{RES}/160_esmfold.json", "w"), indent=2, ensure_ascii=False)
    log("\nDONE → results/160_*")


if __name__ == "__main__":
    main()
