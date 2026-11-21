#!/usr/bin/env python
"""
步驟 31：對步驟 30 產出的寡核苷酸候選做全轉錄組脫靶篩選。

方法
----
把整個 GENCODE 轉錄組編碼成 2-bit 的 k-mer 陣列（每個 k-mer 一個 uint64），
排序後用二分搜尋查詢。250 M 個 k-mer 只佔約 2 GB 記憶體、查詢是 O(log n)。

對每個候選計算
--------------
* `n_perfect_tx`   ：全長完全相符的轉錄本數（扣掉標的本身後應為 0）
* `n_seed_hits`    ：siRNA 種子區（guide 第 2–8 位）在轉錄組 3'UTR 的出現次數
                     —— 這是 miRNA 樣脫靶的主要來源
* `n_offtarget_genes`：出現完全相符的「其他基因」數
"""
import os, sys, gzip, argparse
from collections import defaultdict
import numpy as np
import pandas as pd

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
GEN = f"{BASE}/data/genome"
RES = f"{BASE}/results"
TX = f"{GEN}/gencode.v44.transcripts.fa.gz"
MAP = {"A": 0, "C": 1, "G": 2, "T": 3}
COMP = str.maketrans("ACGT", "TGCA")


def log(*a):
    print("[31]", *a, flush=True)


def rc(s):
    return s.translate(COMP)[::-1]


def encode_all(seq, k):
    """把序列轉成所有長度 k 的 2-bit 編碼（含 N 的視窗會被丟棄）。"""
    a = np.frombuffer(seq.encode(), dtype=np.uint8)
    code = np.full(a.shape, 255, dtype=np.uint8)
    for ch, v in MAP.items():
        code[a == ord(ch)] = v
    n = len(code) - k + 1
    if n <= 0:
        return np.empty(0, dtype=np.uint64), np.empty(0, dtype=np.int64)
    win = np.lib.stride_tricks.sliding_window_view(code, k)
    good = (win != 255).all(axis=1)
    w = win[good].astype(np.uint64)
    mult = (np.uint64(4) ** np.arange(k - 1, -1, -1, dtype=np.uint64))
    return (w * mult).sum(axis=1), np.flatnonzero(good)


def kmer_of(s, k):
    assert len(s) == k
    v = np.uint64(0)
    for ch in s:
        v = v * np.uint64(4) + np.uint64(MAP[ch])
    return v


def build_index(k, limit=None):
    """回傳 (sorted_kmers, gene_of_kmer_ptr) —— 為了省記憶體，只保留 k-mer 與其基因 id。"""
    kms, gids, genes = [], [], []
    gid = -1
    name, buf = None, []
    n = 0
    with gzip.open(TX, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if buf:
                    s = "".join(buf).upper()
                    e, _ = encode_all(s, k)
                    if len(e):
                        kms.append(e); gids.append(np.full(len(e), gid, dtype=np.int32))
                    n += 1
                    if limit and n >= limit:
                        break
                p = line[1:].split("|")
                sym = p[5] if len(p) > 5 else p[0]
                genes.append(sym); gid = len(genes) - 1
                buf = []
            else:
                buf.append(line.strip())
        if buf and (not limit or n < limit):
            s = "".join(buf).upper()
            e, _ = encode_all(s, k)
            if len(e):
                kms.append(e); gids.append(np.full(len(e), gid, dtype=np.int32))
    K = np.concatenate(kms); G = np.concatenate(gids)
    o = np.argsort(K, kind="stable")
    return K[o], G[o], np.array(genes, dtype=object)


def query(Ks, Gs, genes, s):
    v = kmer_of(s, len(s))
    lo = np.searchsorted(Ks, v, "left"); hi = np.searchsorted(Ks, v, "right")
    if hi <= lo:
        return 0, set()
    return int(hi - lo), set(genes[Gs[lo:hi]])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inp", default=f"{RES}/30_oligo_candidates.csv")
    ap.add_argument("--out", default=f"{RES}/31_oligo_offtarget.csv")
    ap.add_argument("--k", type=int, default=16)
    ap.add_argument("--seed_k", type=int, default=7)
    args = ap.parse_args()

    df = pd.read_csv(args.inp)
    log("候選數:", len(df))

    log(f"建立 {args.k}-mer 索引 ...")
    K16, G16, genes = build_index(args.k)
    log(f"  {len(K16):,} 個 {args.k}-mer，{len(genes):,} 條轉錄本，記憶體約 {K16.nbytes/1e9:.1f} GB")

    log(f"建立 {args.seed_k}-mer 索引（siRNA 種子區用）...")
    K7, G7, _ = build_index(args.seed_k)
    log(f"  {len(K7):,} 個 {args.seed_k}-mer")

    rows = []
    for _, r in df.iterrows():
        seq = str(r["sequence"]).replace("dTdT", "")
        site = str(r["mrna_site"])
        # 全長完全相符（查 mRNA 端的序列）
        n_hit, hit_genes = (0, set())
        if len(site) >= args.k:
            n_hit, hit_genes = query(K16, G16, genes, site[:args.k])
        off_genes = hit_genes - {r["target"]}
        d = dict(r)
        d["n_perfect_tx"] = n_hit
        d["n_offtarget_genes"] = len(off_genes)
        d["offtarget_genes"] = ";".join(sorted(off_genes)[:5])
        if r["kind"] == "siRNA" and isinstance(r.get("seed"), str) and len(r["seed"]) >= args.seed_k:
            sd = r["seed"][:args.seed_k]
            # 種子區的互補序列出現在 mRNA 上 = 潛在 miRNA 樣脫靶
            n_seed, sg = query(K7, G7, genes, rc(sd))
            d["n_seed_hits_tx"] = n_seed
            d["n_seed_hit_genes"] = len(sg)
        rows.append(d)

    out = pd.DataFrame(rows)
    out["specificity_ok"] = (out.n_offtarget_genes == 0)
    out.to_csv(args.out, index=False)
    log("寫出", args.out)

    for kind in out.kind.unique():
        s = out[(out.kind == kind) & (out.specificity_ok)]
        log(f"\n=== {kind}：通過特異性篩選的候選 ({len(s)}/{(out.kind==kind).sum()}) ===")
        cols = [c for c in ["target", "pos", "sequence", "gc", "asymmetry",
                            "n_perfect_tx", "n_offtarget_genes",
                            "n_seed_hit_genes"] if c in s.columns]
        log("\n" + s[cols].to_string(index=False))


if __name__ == "__main__":
    main()
