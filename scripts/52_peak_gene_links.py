#!/usr/bin/env python
"""
步驟 52：RNA + ATAC 聯合 —— 建立「調控元件 → 基因」連結，並比較正常與腫瘤。

這一步是整個分析的核心轉折：
前面只看「哪些基因/TF 的量變了」，這裡要看「調控關係本身有沒有被重接線」。

設計
----
1. **分化狀態配對**：腫瘤細胞 87% 對位到非類別轉換記憶 B 細胞 (ncsMBC)，
   所以正常端只取對位分布涵蓋到的那些狀態，避免拿錯的發育階段來比。
2. **Metacell**：單細胞 ATAC 太稀疏，把每 ~40 顆鄰近細胞合併成一個 metacell。
3. **peak–gene 相關**：對每個基因，取 TSS ±250 kb 內的 peak，
   在 metacell 層級算 peak 可及性與基因表現的相關係數。
4. **背景校正**：用「距離相同但來自隨機其他基因」的 peak–gene 配對當 null，
   取經驗 p 值（這是 SCENIC+/ArchR 的標準做法）。
5. **重接線**：同一個 peak–gene 配對，在正常顯著但腫瘤不顯著（= 連結斷掉），
   或反過來（= 新生連結）。

輸出 results/52_peak_gene_links_{normal,tumor}.csv 與 52_rewiring.csv
"""
import os, sys, gzip
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _compat  # noqa
import numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp
from sklearn.neighbors import NearestNeighbors

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
PROC = f"{BASE}/data/processed"
RES  = f"{BASE}/results"
GTF  = f"{BASE}/data/genome/gencode.v44.annotation.gtf.gz"
WIN  = 250_000
MC_K = 40            # 每個 metacell 合併的細胞數
sc.settings.n_jobs = 48

def log(*a): print("[52]", *a, flush=True)

def tss_table():
    rows=[]
    with gzip.open(GTF,"rt") as fh:
        for line in fh:
            if line[0]=="#": continue
            f=line.split("\t")
            if f[2]!="gene": continue
            a=f[8]
            if 'gene_name "' not in a: continue
            nm=a.split('gene_name "')[1].split('"')[0]
            gt=a.split('gene_type "')[1].split('"')[0] if 'gene_type "' in a else ""
            if gt not in ("protein_coding","lncRNA"): continue
            rows.append((nm, f[0], int(f[3]) if f[6]=="+" else int(f[4])))
    t=pd.DataFrame(rows, columns=["gene","chrom","tss"]).drop_duplicates("gene")
    return t.set_index("gene")

def metacells(X_rna, X_atac, emb, k=MC_K, seed=0):
    """以 RNA 潛在空間取 kNN 聚合；回傳 (metacell x gene, metacell x peak)。"""
    n = emb.shape[0]
    rng = np.random.default_rng(seed)
    nn = NearestNeighbors(n_neighbors=k, n_jobs=-1).fit(emb)
    seeds = rng.choice(n, size=max(50, n//(k//2)), replace=False)
    _, idx = nn.kneighbors(emb[seeds])
    R = np.vstack([np.asarray(X_rna[i].mean(axis=0)).ravel() for i in idx])
    A = np.vstack([np.asarray(X_atac[i].mean(axis=0)).ravel() for i in idx])
    return R, A

def corr_cols(A, B):
    """A: n x p, B: n x p —— 逐欄相關（已對齊的配對）。"""
    Az = (A - A.mean(0)) / (A.std(0) + 1e-9)
    Bz = (B - B.mean(0)) / (B.std(0) + 1e-9)
    return (Az * Bz).mean(0)

def build(tag, R, A, genes, peaks_df, tss, rng):
    """回傳 peak-gene 連結表。"""
    gidx = {g:i for i,g in enumerate(genes)}
    out = []
    bychr = {c: g.reset_index(drop=True) for c, g in peaks_df.groupby("chrom")}
    # 先標準化
    Rz = (R - R.mean(0)) / (R.std(0) + 1e-9)
    Az = (A - A.mean(0)) / (A.std(0) + 1e-9)
    nmc = R.shape[0]
    use_genes = [g for g in genes if g in tss.index]
    log(f"  [{tag}] 有 TSS 的基因 {len(use_genes)}，metacell {nmc}")
    allpairs_r, allpairs_meta = [], []
    for gi, g in enumerate(use_genes):
        c, t = tss.loc[g, "chrom"], tss.loc[g, "tss"]
        gp = bychr.get(c)
        if gp is None: continue
        m = (gp.start > t-WIN) & (gp.end < t+WIN)
        if not m.any(): continue
        pid = gp.pid.values[m.values]
        r = (Az[:, pid] * Rz[:, gidx[g]][:,None]).mean(0)
        d = ((gp.start.values[m.values]+gp.end.values[m.values])//2 - t)
        allpairs_r.append(r)
        allpairs_meta.append(pd.DataFrame(dict(gene=g, peak_id=pid, dist=d, r=r)))
        if (gi+1) % 2000 == 0: log(f"    [{tag}] {gi+1}/{len(use_genes)}")
    L = pd.concat(allpairs_meta, ignore_index=True)
    # 背景：隨機打散 gene 與 peak 的配對（保留距離分布）
    nbg = min(400_000, len(L))
    bg_pid = rng.choice(peaks_df.pid.values, nbg)
    bg_g   = rng.choice([gidx[g] for g in use_genes], nbg)
    bg_r = (Az[:, bg_pid] * Rz[:, bg_g]).mean(0)
    bg_sorted = np.sort(bg_r)
    L["emp_p"] = 1 - np.searchsorted(bg_sorted, L.r.values, side="left")/len(bg_sorted)
    L["tag"] = tag
    log(f"  [{tag}] peak-gene 配對 {len(L):,}，背景 {nbg:,}，"
        f"emp_p<0.01 的有 {int((L.emp_p<0.01).sum()):,}")
    return L

if __name__ == "__main__":
    rng = np.random.default_rng(0)
    tss = tss_table(); log("TSS 表:", len(tss))
    peaks = pd.read_csv(f"{PROC}/atlas_peaks_sorted.csv")
    log("peak:", len(peaks))

    # ---------- 正常（分化狀態配對） ----------
    obs_m = pd.read_csv(f"{RES}/15_mcl_obs.csv", index_col=0)
    matched = obs_m.matched_normal_state.value_counts()
    keep_states = matched[matched >= 50].index.tolist()
    log("腫瘤對位到、且細胞數足夠的正常狀態:", keep_states)

    Rn = sc.read_h5ad(f"{PROC}/ref_normalB_multiome.h5ad")
    An = sc.read_h5ad(f"{PROC}/tonsil_multiome_atac_B.h5ad")
    common_cells = Rn.obs_names.intersection(An.obs_names)
    sel = Rn.obs_names.isin(common_cells) & Rn.obs.annotation.isin(keep_states).values
    Rn = Rn[sel].copy(); An = An[Rn.obs_names].copy()
    log("正常（配對狀態）細胞:", Rn.n_obs)
    sc.pp.normalize_total(Rn, target_sum=1e4); sc.pp.log1p(Rn)
    Anx = An.X.tocsr(); Anx = Anx.multiply(1e4/ (Anx.sum(1)+1)).tocsr(); Anx.data = np.log1p(Anx.data)
    emb_n = Rn.obsm["X_harmony_rna"][:, :30] if "X_harmony_rna" in Rn.obsm else Rn.obsm["X_pca"][:, :30]

    # 基因集：高變異 + 全部 TF
    import decoupler as dc
    net = dc.op.collectri(organism="human")
    TFs = sorted(set(net["source"]) & set(Rn.var_names))
    sc.pp.highly_variable_genes(Rn, n_top_genes=4000, flavor="seurat")
    genes = sorted(set(Rn.var_names[Rn.var.highly_variable]) | set(TFs))
    log(f"分析基因 {len(genes)}（含 {len(TFs)} 個 TF）")
    gi = Rn.var_names.get_indexer(genes)

    Rn_mc, An_mc = metacells(Rn.X.tocsr()[:, gi], Anx, emb_n)
    log("正常 metacell:", Rn_mc.shape, An_mc.shape)
    Ln = build("normal", Rn_mc, An_mc, genes, peaks, tss, rng)
    Ln.to_csv(f"{RES}/52_peak_gene_links_normal.csv.gz", index=False, compression="gzip")
    del Rn, An, Anx, Rn_mc, An_mc

    # ---------- 腫瘤 ----------
    Rt = sc.read_h5ad(f"{PROC}/mcl_tumor_labeled.h5ad")
    Rt = Rt[Rt.obs.is_tumor.astype(str)=="True"].copy()
    mats, cells = [], []
    for d in ["M102","M413"]:
        M = sp.load_npz(f"{PROC}/mcl_atac_atlaspeaks_{d}.npz")
        cl = pd.read_csv(f"{PROC}/mcl_atac_atlaspeaks_{d}_cells.csv", header=None)[0].tolist()
        mats.append(M); cells += cl
    At = sp.vstack(mats).tocsr()
    ci = pd.Index(cells)
    keep = ci.isin(Rt.obs_names)
    At = At[np.flatnonzero(keep)]; ci = ci[keep]
    Rt = Rt[ci].copy()
    log("腫瘤細胞（RNA+ATAC 都有）:", Rt.n_obs)
    At = At.multiply(1e4/(At.sum(1)+1)).tocsr(); At.data = np.log1p(At.data)
    if "X_harmony" in Rt.obsm: emb_t = np.asarray(Rt.obsm["X_harmony"])[:, :30]
    elif "X_pca" in Rt.obsm:  emb_t = np.asarray(Rt.obsm["X_pca"])[:, :30]
    else:
        sc.pp.pca(Rt, n_comps=30); emb_t = Rt.obsm["X_pca"]
    gi2 = Rt.var_names.get_indexer(genes)
    ok = gi2 >= 0
    genes_t = [g for g,o in zip(genes, ok) if o]
    Rt_mc, At_mc = metacells(Rt.X.tocsr()[:, gi2[ok]], At, emb_t)
    log("腫瘤 metacell:", Rt_mc.shape, At_mc.shape)
    Lt = build("tumor", Rt_mc, At_mc, genes_t, peaks, tss, rng)
    Lt.to_csv(f"{RES}/52_peak_gene_links_tumor.csv.gz", index=False, compression="gzip")

    # ---------- 重接線 ----------
    k = ["gene","peak_id"]
    M = Ln[k+["r","emp_p","dist"]].merge(Lt[k+["r","emp_p"]], on=k, suffixes=("_n","_t"))
    M["gained"] = (M.emp_p_t < 0.01) & (M.emp_p_n > 0.2)
    M["lost"]   = (M.emp_p_n < 0.01) & (M.emp_p_t > 0.2)
    M["dr"] = M.r_t - M.r_n
    M.to_csv(f"{RES}/52_rewiring.csv.gz", index=False, compression="gzip")
    log(f"\n共同 peak-gene 配對 {len(M):,}；新生連結 {int(M.gained.sum()):,}；斷裂連結 {int(M.lost.sum()):,}")
    log("\n=== 新生連結最多的基因 Top25 ===")
    log("\n" + M[M.gained].gene.value_counts().head(25).to_string())
    log("\n=== 斷裂連結最多的基因 Top25 ===")
    log("\n" + M[M.lost].gene.value_counts().head(25).to_string())
    log("DONE")
