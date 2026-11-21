#!/usr/bin/env python
"""
步驟 24：在指定基因座上，畫出 MCL 腫瘤細胞 vs 正常扁桃腺 B 細胞的
         ATAC 片段覆蓋度剖面（對等比較，兩邊都從 fragments 重新量化）。

這一步的目的是把「某個調控子異常」落實到**具體的基因組座標**上：
到底是啟動子、哪一段增強子、還是一整片新生的可及性區域。

用法
----
python 24_locus_profile.py --genes CCND1 SOX11 IRF4 --flank 400000 --bin 200
"""
import os, sys, gzip, argparse, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
PROC = f"{BASE}/data/processed"
RES = f"{BASE}/results"
GTF = f"{BASE}/data/genome/gencode.v44.annotation.gtf.gz"
MCLD = f"{BASE}/data/extracted/MCL/BCLLATLAS_64_65"
FRAGD = f"{BASE}/data/extracted/fragments"
os.makedirs(RES, exist_ok=True)


def log(*a):
    print("[24]", *a, flush=True)


def gene_locus(gtf, symbols):
    want = set(symbols)
    out = {}
    with gzip.open(gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if f[2] != "gene":
                continue
            a = f[8]
            if 'gene_name "' not in a:
                continue
            nm = a.split('gene_name "')[1].split('"')[0]
            if nm in want and nm not in out:
                out[nm] = dict(chrom=f[0], start=int(f[3]), end=int(f[4]), strand=f[6])
                if len(out) == len(want):
                    break
    return out


def coverage_from_fragments(frag_path, chrom, start, end, binsize,
                            keep_barcodes=None, prefix=""):
    """用 tabix 取區間內片段，回傳每個 bin 的片段數與使用到的細胞數。"""
    import pysam
    nbin = (end - start) // binsize + 1
    cov = np.zeros(nbin, dtype=np.float64)
    seen = set()
    try:
        tbx = pysam.TabixFile(frag_path)
    except Exception as e:
        log("  無法開啟", frag_path, e)
        return cov, 0
    contigs = set(tbx.contigs)
    c = chrom if chrom in contigs else chrom.replace("chr", "")
    if c not in contigs:
        log("  找不到染色體", chrom, "可用前幾個:", list(contigs)[:5])
        return cov, 0
    for row in tbx.fetch(c, max(0, start), end, parse=None):
        f = row.split("\t")
        bc = prefix + f[3]
        if keep_barcodes is not None and bc not in keep_barcodes:
            continue
        s, e = int(f[1]), int(f[2])
        b0 = max(0, (s - start) // binsize)
        b1 = min(nbin - 1, (e - start) // binsize)
        cov[b0:b1 + 1] += 1
        seen.add(bc)
    tbx.close()
    return cov, len(seen)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes", nargs="+", default=["CCND1", "SOX11", "IRF4"])
    ap.add_argument("--flank", type=int, default=400000)
    ap.add_argument("--bin", type=int, default=200)
    ap.add_argument("--out", default=f"{RES}/24_locus_profiles.npz")
    args = ap.parse_args()

    loci = gene_locus(GTF, args.genes)
    log("基因座標:", json.dumps(loci, indent=1))

    # ---- MCL 腫瘤細胞條碼 ----
    obs = pd.read_csv(f"{RES}/20_mcl_obs_labeled.csv", index_col=0)
    tum = obs.index[obs["is_tumor"].astype(str) == "True"]
    meta = pd.read_csv(f"{BASE}/data/extracted/MCL/metadata/cellranger_arc_metadata.csv")
    pos = meta[(meta["type"] == "RNA") & meta["library_name"].str.contains("pos")]
    mcl_frag = {r["donor_id"]: (f"{MCLD}/{r['gem_id']}/atac_fragments.tsv.gz", r["gem_id"] + "_")
                for _, r in pos.iterrows()}
    log("MCL fragments:", {k: v[0].split('/')[-2] for k, v in mcl_frag.items()})

    # ---- 正常 B 細胞條碼（依 gem_id 分檔） ----
    nmd = pd.read_csv(f"{PROC}/multiome_export/multiome_metadata.csv", low_memory=False)
    nB = nmd[nmd["annotation_level_1"].isin(["NBC_MBC", "GCBC", "PC"])]
    log("正常 B 細胞:", len(nB), "| gem_id 數:", nB.gem_id.nunique())
    normal_frag = []
    if os.path.isdir(FRAGD):
        for gem in sorted(nB.gem_id.unique()):
            cands = [os.path.join(dp, f) for dp, _, fs in os.walk(FRAGD)
                     for f in fs if f.endswith(".tsv.gz") and gem in os.path.join(dp, f)]
            if cands:
                normal_frag.append((gem, cands[0]))
        log("找到正常 fragments 檔:", len(normal_frag), "/", nB.gem_id.nunique())
    else:
        log("⚠ 尚未解壓正常 fragments（", FRAGD, "不存在）→ 只輸出 MCL 剖面")

    out = {}
    for g, L in loci.items():
        chrom = L["chrom"]
        tss = L["start"] if L["strand"] == "+" else L["end"]
        s, e = max(0, tss - args.flank), tss + args.flank
        nbin = (e - s) // args.bin + 1
        log(f"\n=== {g} {chrom}:{s}-{e} ({nbin} bins) ===")

        for donor, (fp, pref) in mcl_frag.items():
            keep = set(tum[tum.str.startswith(pref)])
            keep_bare = {k[len(pref):] for k in keep}
            cov, ncell = coverage_from_fragments(fp, chrom, s, e, args.bin,
                                                 keep_barcodes=keep_bare)
            out[f"{g}|MCL_{donor}"] = cov
            out[f"{g}|MCL_{donor}|ncell"] = np.array([ncell])
            log(f"  MCL {donor}: {ncell} 顆腫瘤細胞, 總片段 {cov.sum():.0f}")

        if normal_frag:
            tot = np.zeros(nbin); ncell_all = 0
            for gem, fp in normal_frag:
                sub = nB[nB.gem_id == gem]
                keep_bare = {b.split(gem + "_")[-1] for b in sub.cell_barcode}
                cov, nc = coverage_from_fragments(fp, chrom, s, e, args.bin,
                                                  keep_barcodes=keep_bare)
                tot += cov; ncell_all += nc
            out[f"{g}|normalB"] = tot
            out[f"{g}|normalB|ncell"] = np.array([ncell_all])
            log(f"  正常 B: {ncell_all} 顆, 總片段 {tot.sum():.0f}")

        out[f"{g}|coords"] = np.array([s, e, args.bin])
        out[f"{g}|tss"] = np.array([tss])
        out[f"{g}|chrom"] = np.array([chrom], dtype=object)

    np.savez_compressed(args.out, **out)
    log("\n寫出", args.out)


if __name__ == "__main__":
    main()
