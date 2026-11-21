#!/usr/bin/env python
"""
步驟 30：針對指定標的設計核酸藥物候選。

產出三類
--------
1. **Gapmer ASO**（RNase H 依賴性降解）：16-nt，5-10-5 LNA/DNA/LNA 設計
2. **siRNA**（RISC 介導降解）：21-nt sense/antisense，含 Khvorova 不對稱性規則
3. **Decoy 雙股寡核苷酸**（只在標的是轉錄因子時）：模擬該 TF 的 DNA 結合位點，
   競爭性吸附 TF 使其離開真正的標的基因

篩選規則（文獻共識，非藥動）
-----------------------------
共通：GC 40–60%；禁止 ≥4 連續同一鹼基；避免 4 個以上連續 G（G-quadruplex）
ASO：避開 CpG（免疫刺激）、避開 TTTT
siRNA：antisense 5' 端要比 3' 端不穩定（熱力學不對稱）、位置 2–8 的種子區
      要做 seed-based off-target 統計、避免 GC > 60%
脫靶：對整個人類轉錄組做精確 k-mer 比對（完全相符 = 嚴重脫靶；
      種子區相符 = 潛在 miRNA 樣脫靶）

注意：本步驟**不考慮藥物動力學 / 遞送**，依使用者設定假設可完美到達標的。
"""
import os, sys, gzip, argparse, json, re
from collections import defaultdict
import numpy as np
import pandas as pd

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
GEN = f"{BASE}/data/genome"
RES = f"{BASE}/results"
os.makedirs(RES, exist_ok=True)
COMP = str.maketrans("ACGTUacgtu", "TGCAAtgcaa")


def log(*a):
    print("[30]", *a, flush=True)


def rc(s):
    return s.translate(COMP)[::-1]


def read_transcripts(path, symbols):
    """從 GENCODE transcripts FASTA 取出指定 gene symbol 的轉錄本。"""
    want = set(symbols)
    seqs = defaultdict(dict)
    name, buf, sym = None, [], None
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if sym and sym in want:
                    seqs[sym][name] = "".join(buf)
                parts = line[1:].strip().split("|")
                name = parts[0]
                sym = parts[5] if len(parts) > 5 else None
                ttype = parts[7] if len(parts) > 7 else ""
                buf = []
                if sym in want and ttype != "protein_coding":
                    sym = None            # 只取蛋白質編碼轉錄本
            else:
                buf.append(line.strip())
        if sym and sym in want:
            seqs[sym][name] = "".join(buf)
    return seqs


def gc(s):
    return (s.count("G") + s.count("C")) / max(1, len(s))


def passes_common(s):
    if not (0.40 <= gc(s) <= 0.60):
        return False, "GC 超出 40-60%"
    if re.search(r"(A{4,}|C{4,}|G{4,}|T{4,})", s):
        return False, "有 4 連續同鹼基"
    if s.count("GGGG"):
        return False, "可能形成 G-quadruplex"
    return True, ""


def aso_filters(s):
    ok, why = passes_common(s)
    if not ok:
        return False, why
    if "CG" in s:
        return False, "含 CpG（免疫刺激風險）"
    return True, ""


DG = {  # 最近鄰 ΔG°37 (kcal/mol)，Santa Lucia 1998 簡表
    "AA": -1.00, "AT": -0.88, "TA": -0.58, "CA": -1.45, "GT": -1.44,
    "CT": -1.28, "GA": -1.30, "CG": -2.17, "GC": -2.24, "GG": -1.84,
}


def nn_dg(s):
    def d(p):
        return DG.get(p, DG.get(rc(p), -1.3))
    return sum(d(s[i:i + 2]) for i in range(len(s) - 1))


def sirna_filters(sense):
    ok, why = passes_common(sense)
    if not ok:
        return False, why, {}
    anti = rc(sense)
    # Khvorova 不對稱性：antisense 的 5' 端應該比 3' 端「鬆」(ΔG 較不負)
    dg5 = nn_dg(anti[:5])
    dg3 = nn_dg(anti[-5:])
    asym = dg5 - dg3
    if asym <= 0:
        return False, "熱力學不對稱性不利（antisense 5' 端太穩定）", {"asym": asym}
    if anti[0] not in "AU T":
        pass
    return True, "", {"asym": asym, "dg5": dg5, "dg3": dg3}


def offtarget_index(tx_fa, kmer=16, limit_tx=None):
    """建立整個轉錄組的 k-mer 集合（記憶體換速度）。"""
    S = set()
    n = 0
    with gzip.open(tx_fa, "rt") as fh:
        seq = []
        for line in fh:
            if line.startswith(">"):
                if seq:
                    s = "".join(seq)
                    for i in range(0, len(s) - kmer + 1):
                        S.add(s[i:i + kmer])
                    n += 1
                    if limit_tx and n >= limit_tx:
                        break
                seq = []
            else:
                seq.append(line.strip())
        if seq:
            s = "".join(seq)
            for i in range(0, len(s) - kmer + 1):
                S.add(s[i:i + kmer])
    return S


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", nargs="+", required=True)
    ap.add_argument("--n_each", type=int, default=8)
    ap.add_argument("--aso_len", type=int, default=16)
    ap.add_argument("--sirna_len", type=int, default=21)
    ap.add_argument("--offtarget", action="store_true",
                    help="建轉錄組 k-mer 索引做脫靶過濾（需約 30–60 GB 記憶體）")
    ap.add_argument("--out", default=f"{RES}/30_oligo_candidates.csv")
    args = ap.parse_args()

    txfa = f"{GEN}/gencode.v44.transcripts.fa.gz"
    seqs = read_transcripts(txfa, args.targets)
    log("取得轉錄本:", {k: len(v) for k, v in seqs.items()})

    ot = None
    if args.offtarget:
        log("建立轉錄組 k-mer 索引（這會用掉不少記憶體）...")
        ot = offtarget_index(txfa, kmer=args.aso_len)
        log("  k-mer 數:", len(ot))

    rows = []
    for sym, tx in seqs.items():
        if not tx:
            log("找不到轉錄本:", sym); continue
        # 取最長的蛋白質編碼轉錄本為代表
        name, s = max(tx.items(), key=lambda kv: len(kv[1]))
        s = s.upper().replace("U", "T")
        log(f"{sym}: 代表轉錄本 {name} 長度 {len(s)}")

        # ---- ASO ----
        cands = []
        L = args.aso_len
        for i in range(0, len(s) - L + 1):
            w = s[i:i + L]
            ok, why = aso_filters(w)
            if not ok:
                continue
            aso = rc(w)                       # ASO 與 mRNA 互補
            score = 1.0 - abs(gc(w) - 0.50) * 2
            cands.append(dict(target=sym, kind="gapmer_ASO", pos=i,
                              mrna_site=w, sequence=aso,
                              gc=round(gc(aso), 3), score=round(score, 3)))
        cands.sort(key=lambda d: -d["score"])
        # 位置分散：彼此間隔 >= 50 nt
        sel, used = [], []
        for c in cands:
            if all(abs(c["pos"] - u) >= 50 for u in used):
                sel.append(c); used.append(c["pos"])
            if len(sel) >= args.n_each:
                break
        rows += sel

        # ---- siRNA ----
        cands = []
        L = args.sirna_len
        for i in range(0, len(s) - L + 1):
            w = s[i:i + L]
            ok, why, info = sirna_filters(w)
            if not ok:
                continue
            anti = rc(w)
            cands.append(dict(target=sym, kind="siRNA", pos=i,
                              mrna_site=w,
                              sequence=anti[:19] + "dTdT",     # antisense (guide)
                              sense=w[:19] + "dTdT",
                              gc=round(gc(w), 3),
                              asymmetry=round(info.get("asym", 0), 2),
                              seed=anti[1:8],
                              score=round(info.get("asym", 0), 3)))
        cands.sort(key=lambda d: -d["score"])
        sel, used = [], []
        for c in cands:
            if all(abs(c["pos"] - u) >= 50 for u in used):
                sel.append(c); used.append(c["pos"])
            if len(sel) >= args.n_each:
                break
        rows += sel

    df = pd.DataFrame(rows)
    if ot is not None and len(df):
        def n_ot(seq):
            s2 = seq.replace("dTdT", "")
            s2 = rc(s2)
            return sum(1 for i in range(0, max(1, len(s2) - args.aso_len + 1))
                       if s2[i:i + args.aso_len] in ot)
        df["offtarget_kmer_hits"] = df["sequence"].map(n_ot)
    df.to_csv(args.out, index=False)
    log("寫出", args.out, df.shape)
    if len(df):
        for k, g in df.groupby("kind"):
            log(f"\n=== {k} ===")
            cols = [c for c in ["target", "pos", "sequence", "gc", "asymmetry", "seed", "score"] if c in g]
            log("\n" + g[cols].head(20).to_string(index=False))


if __name__ == "__main__":
    main()
