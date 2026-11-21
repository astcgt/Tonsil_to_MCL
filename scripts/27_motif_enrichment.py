#!/usr/bin/env python
"""
步驟 27：在「MCL 喪失 / 新增」的開放染色質區域中做轉錄因子 motif 富集。

目的
----
把「哪些 TF 的結合位點，系統性地出現在被關掉（或被打開）的調控元件裡」找出來。
這能回答「誰是上游驅動者」——例如若 SOX 家族 motif 顯著富集在 MCL 喪失的 peak 中，
就支持「異位表現的 SOX11 直接壓制這些元件」這個假說。

方法
----
1. 前景：MCL 一致喪失（或新增）的 peak；背景：GC 含量與長度配對的其他 peak。
2. 用 JASPAR2024 CORE vertebrates 非冗餘 PFM，轉成 log-odds PWM，
   以 0.8 × 最高分 為門檻掃描（兩股）。
3. 以超幾何檢定比較前景 vs 背景的命中率。
"""
import os, sys, argparse, math
from collections import defaultdict
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from multiprocessing import Pool

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
RES = f"{BASE}/results"
GEN = f"{BASE}/data/genome"
FASTA = f"{GEN}/GRCh38.primary_assembly.genome.fa"
PFM = f"{GEN}/JASPAR2024_CORE_vertebrates_nr_pfms.txt"
IDX = {"A": 0, "C": 1, "G": 2, "T": 3}


def log(*a):
    print("[27]", *a, flush=True)


def read_jaspar(path):
    """JASPAR 原始格式：>MA0001.1 AGL3 / A [ ... ] 四行"""
    motifs = {}
    with open(path) as fh:
        lines = [l.rstrip("\n") for l in fh]
    i = 0
    while i < len(lines):
        if not lines[i].startswith(">"):
            i += 1; continue
        hdr = lines[i][1:].split()
        mid, name = hdr[0], hdr[1] if len(hdr) > 1 else hdr[0]
        mat = []
        for j in range(1, 5):
            row = lines[i + j]
            row = row[row.find("[") + 1: row.rfind("]")] if "[" in row else row.split(None, 1)[1]
            mat.append([float(x) for x in row.split()])
        motifs[f"{mid}_{name}"] = np.array(mat)       # 4 x L
        i += 5
    return motifs


def pfm_to_pwm(pfm, bg=(0.25, 0.25, 0.25, 0.25), pseudo=0.25):
    c = pfm + pseudo
    p = c / c.sum(axis=0, keepdims=True)
    return np.log2(p / np.array(bg)[:, None])


def seq_to_idx(s):
    a = np.full(len(s), -1, dtype=np.int8)
    for ch, k in IDX.items():
        a[np.frombuffer(s.encode(), dtype="S1") == ch.encode()] = k
    return a


def scan_max(idx, pwm):
    """回傳序列（單股）對該 PWM 的最高分。"""
    L = pwm.shape[1]
    n = len(idx) - L + 1
    if n <= 0:
        return -1e9
    best = -1e9
    # 向量化：建立 (n, L) 的視窗索引
    win = np.lib.stride_tricks.sliding_window_view(idx, L)
    valid = (win >= 0).all(axis=1)
    if not valid.any():
        return -1e9
    w = win[valid]
    sc = pwm[w, np.arange(L)].sum(axis=1)
    return float(sc.max())


def rc_pwm(pwm):
    return pwm[::-1, ::-1]


_G = {}


def _init(seqs, motifs, thr):
    _G["seqs"] = seqs; _G["motifs"] = motifs; _G["thr"] = thr


def _work(mkey):
    pwm = _G["motifs"][mkey]
    pwmr = rc_pwm(pwm)
    tmax = pwm.max(axis=0).sum()
    cut = _G["thr"] * tmax
    hits = np.zeros(len(_G["seqs"]), dtype=bool)
    for i, idx in enumerate(_G["seqs"]):
        s1 = scan_max(idx, pwm)
        if s1 >= cut:
            hits[i] = True; continue
        if scan_max(idx, pwmr) >= cut:
            hits[i] = True
    return mkey, hits


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n_fg", type=int, default=1500)
    ap.add_argument("--n_bg", type=int, default=6000)
    ap.add_argument("--thr", type=float, default=0.80)
    ap.add_argument("--procs", type=int, default=96)
    args = ap.parse_args()

    import pyfaidx
    fa = pyfaidx.Fasta(FASTA, as_raw=True, sequence_always_upper=True)
    d = pd.read_csv(f"{RES}/22_atac_diff.csv")
    d = d[d.consistent].copy()
    d["len"] = d.end - d.start
    d = d[(d.len > 150) & (d.len < 3000)]
    log("可用 peak:", len(d))

    lost = d[(d.normal_log2cpm > 2)].sort_values("lfc_mean").head(args.n_fg)
    gain = d[(d.mcl_mean_log2cpm > 2)].sort_values("lfc_mean", ascending=False).head(args.n_fg)
    log("前景 lost:", len(lost), " gain:", len(gain))

    def grab(df):
        out = []
        for _, r in df.iterrows():
            try:
                s = str(fa[r.chrom][int(r.start):int(r.end)])
            except Exception:
                continue
            if len(s) < 50:
                continue
            out.append(s)
        return out

    # ---- 長度 + GC 配對背景（關鍵：不配對的話長 peak 會讓所有 motif 都「富集」）----
    rng = np.random.default_rng(0)
    rest = d[~d.peak.isin(set(lost.peak) | set(gain.peak))].copy()

    def gc_of(df):
        g = []
        for _, r in df.iterrows():
            try:
                s = str(fa[r.chrom][int(r.start):int(r.end)])
                g.append((s.count("G") + s.count("C")) / max(1, len(s)))
            except Exception:
                g.append(np.nan)
        return np.array(g)

    log("計算背景池的 GC ...")
    pool = rest.sample(n=min(60000, len(rest)), random_state=1).reset_index(drop=True)
    pool["gcfrac"] = gc_of(pool)
    pool = pool.dropna(subset=["gcfrac"])
    log("背景池:", len(pool))

    def matched_bg(fg, k=3):
        """每個前景 peak 配 k 個長度(±15%)與 GC(±0.03)相近的背景 peak。"""
        fg = fg.copy(); fg["gcfrac"] = gc_of(fg)
        used, picks = set(), []
        L = pool["len"].values; G = pool["gcfrac"].values
        for _, r in fg.iterrows():
            if np.isnan(r.gcfrac):
                continue
            ok = (np.abs(L - r.len) <= 0.15 * r.len) & (np.abs(G - r.gcfrac) <= 0.03)
            cand = np.flatnonzero(ok)
            cand = [c for c in cand if c not in used]
            if len(cand) == 0:                      # 放寬
                ok = (np.abs(L - r.len) <= 0.35 * r.len) & (np.abs(G - r.gcfrac) <= 0.06)
                cand = [c for c in np.flatnonzero(ok) if c not in used]
            if len(cand) == 0:
                continue
            take = rng.choice(cand, size=min(k, len(cand)), replace=False)
            used.update(take.tolist()); picks.extend(take.tolist())
        return pool.iloc[picks]

    bg_lost = matched_bg(lost)
    bg_gain = matched_bg(gain)
    log("配對背景: lost->%d, gain->%d" % (len(bg_lost), len(bg_gain)))

    sets = {"lost": grab(lost), "gain": grab(gain),
            "bg_lost": grab(bg_lost), "bg_gain": grab(bg_gain)}
    for k, v in sets.items():
        log(f"  {k}: {len(v)} 條序列，平均長度 {np.mean([len(x) for x in v]):.0f}")
    gcs = {k: np.mean([(s.count("G") + s.count("C")) / len(s) for s in v]) for k, v in sets.items()}
    log("  GC 含量:", {k: round(v, 3) for k, v in gcs.items()})

    motifs = {k: pfm_to_pwm(v) for k, v in read_jaspar(PFM).items()}
    log("JASPAR motif 數:", len(motifs))

    results = {}
    for key in ["lost", "gain", "bg_lost", "bg_gain"]:
        idxs = [seq_to_idx(s) for s in sets[key]]
        with Pool(args.procs, initializer=_init, initargs=(idxs, motifs, args.thr)) as p:
            res = dict(p.map(_work, list(motifs.keys()), chunksize=4))
        results[key] = {m: int(h.sum()) for m, h in res.items()}
        log(f"  {key} 掃描完成")

    rows = []
    for fgkey in ["lost", "gain"]:
        nfg = len(sets[fgkey]); nbg = len(sets["bg_" + fgkey])
        for m in motifs:
            a = results[fgkey][m]; b = results["bg_" + fgkey][m]
            if a + b == 0:
                continue
            # 超幾何：在 (nfg+nbg) 個區域中有 a+b 個命中，抽 nfg 個，觀察到 a 個
            pv = hypergeom.sf(a - 1, nfg + nbg, a + b, nfg)
            fg_rate = a / nfg; bg_rate = b / nbg
            rows.append(dict(set=fgkey, motif=m,
                             fg_hits=a, fg_n=nfg, fg_rate=round(fg_rate, 4),
                             bg_hits=b, bg_n=nbg, bg_rate=round(bg_rate, 4),
                             fold=round(fg_rate / max(bg_rate, 1e-6), 3),
                             pvalue=pv))
    df = pd.DataFrame(rows)
    df["padj"] = df.groupby("set")["pvalue"].transform(
        lambda p: np.minimum(1, p * len(p)))
    df = df.sort_values(["set", "pvalue"])
    df.to_csv(f"{RES}/27_motif_enrichment.csv", index=False)
    for k in ["lost", "gain"]:
        s = df[(df.set == k) & (df.fold > 1)].head(30)
        log(f"\n=== {k}: 富集最顯著的 30 個 motif ===")
        log("\n" + s[["motif", "fg_rate", "bg_rate", "fold", "pvalue", "padj"]].to_string(index=False))
    log("DONE")


if __name__ == "__main__":
    main()
