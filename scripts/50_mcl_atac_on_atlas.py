#!/usr/bin/env python
"""
步驟 50：把 MCL 的 ATAC fragments 重新量化到「扁桃腺圖譜的 peak 座標」上（單細胞層級）。

為什麼必須做
------------
MCL 與正常各自 call 出不同的 peak 集合（90k / 109k vs 215k），
特徵空間不同就無法做聯合的調控網絡比較。
這裡直接從 fragments 重新計數，讓兩邊落在同一套座標上。

做法：每條染色體各自處理，用 searchsorted 把 fragment 的中點對到 peak 區間，
以 (cell, peak) 累加成稀疏矩陣。以染色體平行化。
"""
import os, sys, gzip, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd, scipy.sparse as sp
from multiprocessing import Pool

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
PROC = f"{BASE}/data/processed"
MCLD = f"{BASE}/data/extracted/MCL/BCLLATLAS_64_65"
PEAKS = f"{PROC}/tonsil_atlas_peaks.bed"

def log(*a): print("[50]", *a, flush=True)

pk = pd.read_csv(PEAKS, sep="\t", header=None, names=["chrom","start","end"])
pk["pid"] = np.arange(len(pk))
pk = pk.sort_values(["chrom","start"]).reset_index(drop=True)
log("圖譜 peak:", len(pk))
BYCHR = {c: g.reset_index(drop=True) for c, g in pk.groupby("chrom")}

def one_chrom(args):
    frag_path, chrom, bc_index = args
    g = BYCHR.get(chrom)
    if g is None: return None
    starts = g.start.values; ends = g.end.values; pids = g.pid.values
    import pysam
    try:
        tbx = pysam.TabixFile(frag_path)
    except Exception:
        return None
    if chrom not in set(tbx.contigs):
        tbx.close(); return None
    rows, cols = [], []
    for line in tbx.fetch(chrom):
        f = line.split("\t")
        bi = bc_index.get(f[3])
        if bi is None: continue
        mid = (int(f[1]) + int(f[2])) // 2
        j = np.searchsorted(starts, mid, side="right") - 1
        if j >= 0 and mid < ends[j]:
            rows.append(bi); cols.append(pids[j])
    tbx.close()
    return (np.array(rows, dtype=np.int32), np.array(cols, dtype=np.int32))

if __name__ == "__main__":
    meta = pd.read_csv(f"{BASE}/data/extracted/MCL/metadata/cellranger_arc_metadata.csv")
    pos = meta[(meta["type"]=="RNA") & meta["library_name"].str.contains("pos")]
    obs_lab = pd.read_csv(f"{BASE}/results/20_mcl_obs_labeled.csv", index_col=0)
    chroms = sorted(BYCHR.keys())
    out = {}
    for _, r in pos.iterrows():
        gem, donor = r["gem_id"], r["donor_id"]
        pref = gem + "_"
        cells = [c for c in obs_lab.index if c.startswith(pref)]
        bare = [c[len(pref):] for c in cells]
        bc_index = {b: i for i, b in enumerate(bare)}
        fp = f"{MCLD}/{gem}/atac_fragments.tsv.gz"
        log(f"{donor}: {len(cells)} 顆細胞, {fp}")
        with Pool(24) as p:
            res = p.map(one_chrom, [(fp, c, bc_index) for c in chroms])
        R = np.concatenate([x[0] for x in res if x is not None])
        C = np.concatenate([x[1] for x in res if x is not None])
        M = sp.coo_matrix((np.ones(len(R), np.float32), (R, C)),
                          shape=(len(cells), len(pk))).tocsr()
        M.sum_duplicates()
        log(f"  矩陣 {M.shape}, 非零 {M.nnz:,}, 總 fragment {M.sum():,.0f}, "
            f"每顆細胞中位數 {np.median(np.asarray(M.sum(1)).ravel()):.0f}")
        sp.save_npz(f"{PROC}/mcl_atac_atlaspeaks_{donor}.npz", M)
        pd.Series(cells).to_csv(f"{PROC}/mcl_atac_atlaspeaks_{donor}_cells.csv", index=False, header=False)
        out[donor] = M.shape
    pk.to_csv(f"{PROC}/atlas_peaks_sorted.csv", index=False)
    log("完成:", out)
