#!/usr/bin/env python
"""步驟 132（第 8 階）：設計品的**特異性**——它只黏 MafG，還是誰都黏？

為什麼這是最大的計算風險
------------------------
設計瞄準的是 MafG 的 DNA 結合鹼性區，而這正是 small Maf 家族最保守的一段。
MAFF / MAFK 與 MAFG 的 bZIP 結構域座標完全相同（UniProt 都是 51–114）。
如果設計品對三個 small Maf 一視同仁，那它就不是「MafG 抑制劑」，
而是「small Maf 家族抑制劑」—— 治療邏輯與毒性評估都要重寫。
先前所有評估都只跑過 MafG 一個標的，所以這件事完全沒被測過。

設計
----
標的（A 鏈，7 個）：都取「bZIP 結構域往前 28、往後 7」的同源片段，所以長度可比
  MAFG   on-target
  MAFF   旁系同源（最嚴格的特異性考驗）
  MAFK   旁系同源
  BACH1  共用同一個 TGACTCA 位點的夥伴
  NFE2L2 天然夥伴（NRF2）
  JUN    AP-1，同樣讀 TGACTCA
  FOS    AP-1
結合劑（B 鏈，4 個）：最佳設計、次佳設計、天然 NRF2（陽性對照）、打散序列（陰性對照）

7 × 4 = 28 個複合體 × 5 次 diffusion 取樣。

判讀
----
特異性 = on-target(MAFG) 分數 − 最高的 off-target 分數。
同時報告標的自己的 pLDDT：分數低如果只是「模型不會摺這個標的」，那不是特異性。

輸出 results/132_*
"""
import os
import sys
import json
import glob
import random
import argparse
import subprocess

import numpy as np
import pandas as pd

PAE_CUT = 10.0          # 與 Overath / Bennett 一致

# --- 以下四個函式逐字取自 scripts/93_ipsae.py（該檔沒有 __main__ 保護，
#     直接 import 會把它整份分析跑一遍，所以在這裡內聯） ---
def calc_d0(L):
    L = np.maximum(26.0, np.asarray(L, dtype=float))
    return np.maximum(1.0, 1.24 * (L - 15.0) ** (1.0 / 3.0) - 1.8)


def atoms_from_cif(path):
    """取每個殘基的鏈代號、CB 座標（無 CB 用 CA）與 pLDDT。

    pDockQ / pDockQ2 需要座標與 pLDDT；ipSAE 只需要鏈代號。
    """
    hdr, rows = [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith("_atom_site."):
                hdr.append(line.strip().split(".", 1)[1])
            elif line.startswith(("ATOM", "HETATM")):
                q = line.split()
                if len(q) >= len(hdr):
                    rows.append(q[:len(hdr)])
    if not rows:
        return None
    df = pd.DataFrame(rows, columns=hdr)
    an = "label_atom_id" if "label_atom_id" in df else "auth_atom_id"
    ch = "label_asym_id" if "label_asym_id" in df else "auth_asym_id"
    rn = "label_seq_id" if "label_seq_id" in df else "auth_seq_id"
    bf = "B_iso_or_equiv" if "B_iso_or_equiv" in df else None
    for c in ("Cartn_x", "Cartn_y", "Cartn_z"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df[rn] = pd.to_numeric(df[rn], errors="coerce")
    if bf:
        df[bf] = pd.to_numeric(df[bf], errors="coerce")
    # 每個殘基優先取 CB，無 CB（甘胺酸）取 CA
    sel = df[df[an].isin(["CA", "CB"])].copy()
    sel["pref"] = (sel[an] == "CB").astype(int)
    sel = sel.sort_values([ch, rn, "pref"]).drop_duplicates(subset=[ch, rn], keep="last")
    xyz = sel[["Cartn_x", "Cartn_y", "Cartn_z"]].values
    pl = sel[bf].values if bf else np.full(len(sel), np.nan)
    return sel[ch].values, sel[rn].values, xyz, pl


def chains_from_cif(path):
    """取每個殘基的鏈代號（只用 CA 原子，一殘基一筆）。"""
    hdr, rows = [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith("_atom_site."):
                hdr.append(line.strip().split(".", 1)[1])
            elif line.startswith(("ATOM", "HETATM")):
                p = line.split()
                if len(p) >= len(hdr):
                    rows.append(p[:len(hdr)])
    if not rows:
        return None
    df = pd.DataFrame(rows, columns=hdr)
    an = "label_atom_id" if "label_atom_id" in df else "auth_atom_id"
    ch = "label_asym_id" if "label_asym_id" in df else "auth_asym_id"
    rn = "label_seq_id" if "label_seq_id" in df else "auth_seq_id"
    df = df[df[an] == "CA"]
    df = df.drop_duplicates(subset=[ch, rn])
    return df[ch].values, pd.to_numeric(df[rn], errors="coerce").values


def ipsae_pair(pae, chains, a, b, cut=PAE_CUT):
    """回傳 (asym A→B, asym B→A, min, max)。"""
    out = {}
    for c1, c2 in [(a, b), (b, a)]:
        m = np.outer(chains == c1, chains == c2) & (pae < cut)
        n0 = m.sum(axis=1)
        d0 = calc_d0(n0)
        best = 0.0
        rows = np.where((chains == c1) & (n0 > 0))[0]
        for i in rows:
            v = 1.0 / (1.0 + (pae[i, m[i]] / d0[i]) ** 2)
            best = max(best, float(v.mean()))
        out[(c1, c2)] = best
    x, y = out[(a, b)], out[(b, a)]
    return x, y, min(x, y), max(x, y)


def interface_pae(pae, chains, a, b):
    """介面 PAE：跨鏈 PAE 的平均（兩個方向都算）。BindCraft 的 i_pAE 對應量。"""
    m = np.outer(chains == a, chains == b) | np.outer(chains == b, chains == a)
    return float(pae[m].mean()) if m.any() else np.nan

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
RES = f"{BASE}/results"
AUX = f"{BASE}/data/aux/uniprot"
OUT = f"{RES}/132_specificity"
UP, DOWN = 28, 7          # bZIP 結構域前後取多少殘基


def log(*a):
    print("[132]", *a, flush=True)


def frag(gene):
    d = json.load(open(f"{AUX}/{gene}.json"))
    s = d["sequence"]["value"]
    dom = [f for f in d.get("features", [])
           if f["type"] == "Domain" and f.get("description") == "bZIP"]
    if not dom:
        raise SystemExit(f"{gene}: 找不到 bZIP 結構域註釋")
    a = dom[0]["location"]["start"]["value"]
    b = dom[0]["location"]["end"]["value"]
    lo, hi = max(1, a - UP), min(len(s), b + DOWN)
    return s[lo - 1:hi], lo, hi, d.get("primaryAccession", "")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gpu", default="4")
    ap.add_argument("--samples", type=int, default=5)
    ap.add_argument("--score_only", action="store_true")
    ap.add_argument("--tag", default="",
                    help="輸出加後綴（例如 _n20），用來在不覆蓋已發表的 n=5 "
                         "結果與 MD 起始結構的前提下，用更多取樣重算一次。")
    a = ap.parse_args()

    global OUT
    OUT = f"{RES}/132_specificity{a.tag}"
    os.makedirs(OUT, exist_ok=True)
    yml = f"{OUT}/inputs"
    os.makedirs(yml, exist_ok=True)

    # ---------- 標的 ----------
    TARGETS = ["MAFG", "MAFF", "MAFK", "BACH1", "NFE2L2", "JUN", "FOS"]
    tseq = {}
    log("標的片段（bZIP 結構域 −%d / +%d）:" % (UP, DOWN))
    for g in TARGETS:
        s, lo, hi, acc = frag(g)
        tseq[g] = s
        log("  %-7s %s  %3d–%-3d  %3d aa" % (g, acc, lo, hi, len(s)))

    # ---------- 結合劑 ----------
    meta87 = json.load(open(f"{RES}/87_meta.json"))
    NRF2 = ("HLTRDELRAKALHIPFPVEKIINLPVVDFNEMMSKEQFNEAQLALIRDIRRRRGKNKVAAQNCRKRK"
            "LENNIVELEQDLDHLKDEKEKLLKEKGENDKSLHLLKKQLS")
    D1 = meta87["R1_mafg_dna_block_90_2_model_1_m4"]["seq"]
    D2 = meta87["R1_mafg_dna_block_90_1_model_0_m1"]["seq"]
    rng = random.Random(1)
    sc = list(D1); rng.shuffle(sc); SCRAM = "".join(sc)
    BINDERS = {"D1_best": D1, "D2_second": D2, "CTRLPOS_NRF2": NRF2, "CTRLNEG_scram": SCRAM}
    for k, v in BINDERS.items():
        log("結合劑 %-14s %3d aa" % (k, len(v)))

    combos = {}
    for t in TARGETS:
        for b, bs in BINDERS.items():
            combos[f"{t}__{b}"] = (tseq[t], bs)

    if not a.score_only:
        for f in glob.glob(f"{yml}/*.yaml"):
            os.remove(f)
        for nm, (ts, bs) in combos.items():
            with open(f"{yml}/{nm}.yaml", "w") as fh:
                fh.write("version: 1\nsequences:\n")
                fh.write(f"  - protein:\n      id: A\n      sequence: {ts}\n      msa: empty\n")
                fh.write(f"  - protein:\n      id: B\n      sequence: {bs}\n      msa: empty\n")
        log("%d 個複合體 × %d 次取樣 → 跑 Boltz-2（GPU %s）" % (len(combos), a.samples, a.gpu))
        env = dict(os.environ, CUDA_VISIBLE_DEVICES=a.gpu)
        # env/boltz/bin/boltz 的 shebang 還指向舊的專案路徑（Tonsil_MCL_diff），
        # 所以用模組進入點呼叫，不依賴那個 console script
        cmd = [f"{BASE}/env/boltz/bin/python", "-c",
               "import sys; from boltz.main import cli; sys.argv[0]='boltz'; cli()",
               "predict", yml, "--out_dir", OUT,
               "--cache", f"{BASE}/models/boltz", "--output_format", "mmcif",
               "--diffusion_samples", str(a.samples), "--recycling_steps", "3",
               "--devices", "1", "--accelerator", "gpu", "--no_kernels"]
        p = subprocess.run(cmd, env=env, capture_output=True, text=True)
        log("boltz return code:", p.returncode)
        if p.returncode != 0:
            log(p.stdout[-1500:]); log(p.stderr[-3000:])

    # ---------- 打分：confidence + ipSAE ----------
    rows = []
    for f in sorted(glob.glob(f"{OUT}/**/confidence_*.json", recursive=True)):
        d = json.load(open(f))
        base = os.path.basename(f)[len("confidence_"):].rsplit(".json", 1)[0]
        nm, _, mdl = base.rpartition("_model_")
        rec = dict(name=nm, model=mdl, iptm=d.get("iptm"),
                   iplddt=d.get("complex_iplddt"), ipde=d.get("complex_ipde"),
                   plddt=d.get("complex_plddt"))
        cif = f.replace("confidence_", "").replace(".json", ".cif")
        pae = os.path.join(os.path.dirname(f), f"pae_{base}.npz")
        if os.path.exists(cif) and os.path.exists(pae):
            try:
                P = np.load(pae)["pae"]
                ch, rn, xyz, pl = atoms_from_cif(cif)
                _ab, _ba, ipmin, _mx = ipsae_pair(P, ch, "A", "B")
                rec["ipsae_min"] = float(ipmin)
                rec["ipae"] = float(interface_pae(P, ch, "A", "B"))
                rec["target_plddt"] = float(np.mean([v for c, v in zip(ch, pl) if c == "A"]))
                rec["binder_plddt"] = float(np.mean([v for c, v in zip(ch, pl) if c == "B"]))
            except Exception as e:  # noqa
                log("ipSAE 失敗", base, type(e).__name__, e)
        rows.append(rec)

    if not rows:
        log("沒有結果，結束"); return
    D = pd.DataFrame(rows)
    D["target"] = D.name.str.split("__").str[0]
    D["binder"] = D.name.str.split("__").str[1]
    D.to_csv(f"{RES}/132_specificity_raw{a.tag}.csv", index=False)

    agg = {c: ["mean", "std"] for c in ["iptm", "ipsae_min", "ipde", "ipae",
                                        "target_plddt", "binder_plddt"] if c in D}
    R = D.groupby(["binder", "target"]).agg(n=("iptm", "size"), **{
        f"{c}_{s}": (c, s) for c, ss in agg.items() for s in ss}).reset_index()
    R.to_csv(f"{RES}/132_specificity_summary{a.tag}.csv", index=False)

    pd.set_option("display.width", 260)
    for metric in ["iptm_mean", "ipsae_min_mean", "ipde_mean"]:
        if metric not in R:
            continue
        piv = R.pivot(index="binder", columns="target", values=metric)
        piv = piv[[c for c in ["MAFG", "MAFF", "MAFK", "BACH1", "NFE2L2", "JUN", "FOS"]
                   if c in piv.columns]]
        log("\n=== %s（列 = 結合劑，欄 = 標的） ===\n" % metric, piv.round(3).to_string())

    tp = R.pivot(index="binder", columns="target", values="target_plddt_mean")
    log("\n=== 標的自己的 pLDDT（低 = 模型不會摺這個標的，分數低不算特異） ===\n",
        tp.round(3).to_string())

    out = {"targets": {g: {"seq": tseq[g], "len": len(tseq[g])} for g in TARGETS},
           "binders": {k: {"len": len(v)} for k, v in BINDERS.items()},
           "samples": a.samples, "specificity": {}}
    for metric, better_high in [("iptm_mean", True), ("ipsae_min_mean", True),
                                ("ipde_mean", False)]:
        if metric not in R:
            continue
        piv = R.pivot(index="binder", columns="target", values=metric)
        for b in piv.index:
            off = piv.loc[b, [c for c in piv.columns if c != "MAFG"]]
            on = piv.loc[b, "MAFG"]
            best_off = off.max() if better_high else off.min()
            gap = (on - best_off) if better_high else (best_off - on)
            out["specificity"].setdefault(str(b), {})[metric] = {
                "on_target_MAFG": float(on),
                "best_off_target": str(off.idxmax() if better_high else off.idxmin()),
                "best_off_value": float(best_off),
                "gap_on_minus_best_off": float(gap),
            }
    log("\n=== 特異性（on-target 減掉最強的 off-target） ===")
    for b, dd in out["specificity"].items():
        for metric, v in dd.items():
            log("  %-14s %-16s MAFG %.3f vs %-6s %.3f  → 差距 %+.3f"
                % (b, metric, v["on_target_MAFG"], v["best_off_target"],
                   v["best_off_value"], v["gap_on_minus_best_off"]))
    json.dump(out, open(f"{RES}/132_specificity{a.tag}.json", "w"), indent=2, ensure_ascii=False)
    log("\nDONE → results/132_*")


if __name__ == "__main__":
    main()
