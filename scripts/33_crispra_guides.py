#!/usr/bin/env python
"""
步驟 33：設計「重新打開 CIITA 啟動子 III」用的 CRISPR 活化（CRISPRa）導引 RNA。

為什麼是這個模態
----------------
我們的證據顯示 CIITA 在 MCL 中被關掉**不是因為基因被刪除或突變**
（chr16 沒有拷貝數缺失、細胞株中也沒有 CIITA 突變），
而是啟動子 III 這段序列失去染色質可及性。
這種「序列還在、只是被關掉」的狀況，最直接的反轉方式就是把一個
序列特異的轉錄活化因子（dCas9-VP64 / dCas9-SAM / dCas9-p300）
導引到該啟動子上 —— 導引 RNA 本身就是一個 RNA 藥物。

設計規則（依 CRISPRa 文獻共識）
-------------------------------
* 最佳視窗：距離轉錄起始點 **-400 ~ -50 bp**（dCas9-VP64/SAM）
* PAM：SpCas9 的 NGG，protospacer 20 nt
* GC 含量 40–80%
* 避開 TTTT（Pol III 終止子）
* 避開在基因組其他地方有完全相符的 protospacer+PAM（脫靶）

同時也設計對照：鎖定 SOX11 啟動子的 CRISPR 干擾（CRISPRi）導引，
用來關掉異位活化的 SOX11。
"""
import os, sys, re, argparse
import numpy as np
import pandas as pd

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
GEN = f"{BASE}/data/genome"
RES = f"{BASE}/results"
FASTA = f"{GEN}/GRCh38.primary_assembly.genome.fa"
COMP = str.maketrans("ACGT", "TGCA")


def log(*a):
    print("[33]", *a, flush=True)


def rc(s):
    return s.translate(COMP)[::-1]


TARGETS = [
    # (名稱, 染色體, TSS, 股, 模態, 視窗起, 視窗迄  —— 相對 TSS)
    ("CIITA_pIII", "chr16", 10877202, "+", "CRISPRa", -400, -50),
    ("CIITA_pIII_wide", "chr16", 10877202, "+", "CRISPRa", -600, +100),
    ("SOX11_prom", "chr2", 5692384, "+", "CRISPRi", -150, +250),
    ("IRF4_prom", "chr6", 391739, "+", "CRISPRi", -150, +250),
]


def find_guides(seq, offset, chrom, tss, strand):
    out = []
    for m in re.finditer(r"(?=([ACGT]{21}GG))", seq):
        w = m.group(1)
        proto, pam = w[1:21], w[21:]
        g = _check(proto)
        if g is None:
            continue
        pos = offset + m.start() + 1
        out.append(dict(chrom=chrom, genomic_start=pos, strand_of_guide="+",
                        protospacer=proto, pam="NGG",
                        dist_to_tss=(pos - tss) if strand == "+" else (tss - pos),
                        gc=round((proto.count("G") + proto.count("C")) / 20, 2)))
    rseq = rc(seq)
    L = len(seq)
    for m in re.finditer(r"(?=([ACGT]{21}GG))", rseq):
        w = m.group(1)
        proto = w[1:21]
        g = _check(proto)
        if g is None:
            continue
        pos = offset + (L - (m.start() + 23))
        out.append(dict(chrom=chrom, genomic_start=pos, strand_of_guide="-",
                        protospacer=proto, pam="NGG",
                        dist_to_tss=(pos - tss) if strand == "+" else (tss - pos),
                        gc=round((proto.count("G") + proto.count("C")) / 20, 2)))
    return out


def _check(p):
    gc = (p.count("G") + p.count("C")) / len(p)
    if not (0.40 <= gc <= 0.80):
        return None
    if "TTTT" in p:
        return None
    if re.search(r"(A{5,}|C{5,}|G{5,}|T{5,})", p):
        return None
    return p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=f"{RES}/33_crispr_guides.csv")
    ap.add_argument("--offtarget", action="store_true", default=True)
    args = ap.parse_args()

    import pyfaidx
    fa = pyfaidx.Fasta(FASTA, as_raw=True, sequence_always_upper=True)

    rows = []
    for name, chrom, tss, strand, mode, w0, w1 in TARGETS:
        s = tss + w0 if strand == "+" else tss - w1
        e = tss + w1 if strand == "+" else tss - w0
        seq = str(fa[chrom][s:e])
        gs = find_guides(seq, s, chrom, tss, strand)
        for g in gs:
            g["target"] = name
            g["mode"] = mode
        log(f"{name} ({chrom}:{s}-{e}, {mode}): 找到 {len(gs)} 條合格導引")
        rows += gs
    df = pd.DataFrame(rows)
    if len(df) == 0:
        log("沒有找到合格導引"); return

    # --- 全基因組脫靶：protospacer 的 20-mer 在基因組出現幾次 ---
    if args.offtarget:
        log("建立基因組 20-mer 索引做脫靶檢查（只檢查完全相符）...")
        MAPB = {"A": 0, "C": 1, "G": 2, "T": 3}

        def enc(s):
            v = 0
            for ch in s:
                if ch not in MAPB:
                    return None
                v = v * 4 + MAPB[ch]
            return v

        want = {}
        for p in df.protospacer.unique():
            for q in (p, rc(p)):
                v = enc(q)
                if v is not None:
                    want.setdefault(v, set()).add(p)
        counts = {p: 0 for p in df.protospacer.unique()}
        K = 20
        mult = (np.uint64(4) ** np.arange(K - 1, -1, -1, dtype=np.uint64))
        wanted = np.array(sorted(want.keys()), dtype=np.uint64)
        for chrom in [f"chr{i}" for i in list(range(1, 23))] + ["chrX", "chrY"]:
            if chrom not in fa:
                continue
            seqc = str(fa[chrom][:])
            a = np.frombuffer(seqc.encode(), dtype=np.uint8)
            code = np.full(a.shape, 255, dtype=np.uint8)
            for ch, v in MAPB.items():
                code[a == ord(ch)] = v
            del a, seqc
            win = np.lib.stride_tricks.sliding_window_view(code, K)
            step = 10_000_000          # 分塊，控制記憶體
            for i0 in range(0, win.shape[0], step):
                w = win[i0:i0 + step]
                good = (w != 255).all(axis=1)
                if not good.any():
                    continue
                vals = (w[good].astype(np.uint64) * mult).sum(axis=1)
                pos = np.searchsorted(wanted, vals)
                pos = np.clip(pos, 0, len(wanted) - 1)
                hit = wanted[pos] == vals
                if not hit.any():
                    continue
                for v, c in zip(*np.unique(vals[hit], return_counts=True)):
                    for p in want[int(v)]:
                        counts[p] += int(c)
            del code, win
            log(f"  {chrom} 掃描完成")
        df["genome_exact_hits"] = df.protospacer.map(counts)
        df["unique_in_genome"] = df.genome_exact_hits <= 1

    df = df.sort_values(["target", "dist_to_tss"])
    df.to_csv(args.out, index=False)
    log("寫出", args.out, df.shape)
    for t, g in df.groupby("target"):
        ok = g[g.get("unique_in_genome", True)] if "unique_in_genome" in g else g
        log(f"\n=== {t}（{g['mode'].iloc[0]}）通過脫靶篩選 {len(ok)}/{len(g)} ===")
        cols = [c for c in ["genomic_start", "strand_of_guide", "protospacer",
                            "dist_to_tss", "gc", "genome_exact_hits"] if c in ok.columns]
        log("\n" + ok[cols].head(15).to_string(index=False))


if __name__ == "__main__":
    main()
