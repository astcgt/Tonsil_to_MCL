#!/usr/bin/env python
"""
步驟 79：建立 eRegulon（TF → 調控元件 → 目標基因），完成聯合 GRN 的最後一塊。

前面已經有：
  * peak–基因連結（正常 / 腫瘤各一套，共 217,705 組）
  * 哪些連結是新生、哪些是斷裂
這一步補上「是誰坐在那個元件上」：
  對每個新生 / 斷裂的 peak 掃 879 個 JASPAR motif，
  再要求 TF 本身在腫瘤中有表現，組成三元組 TF → peak → gene。

輸出：每個 TF 在腫瘤中「獲得 / 失去」多少目標基因 = 調控網絡重接線的 TF 層級摘要。
"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
from multiprocessing import Pool

BASE="/raid/jupyter-alz.ee09/Tonsil_to_MCL"; RES=f"{BASE}/results"; GEN=f"{BASE}/data/genome"
FASTA=f"{GEN}/GRCh38.primary_assembly.genome.fa"; PFM=f"{GEN}/JASPAR2024_CORE_vertebrates_nr_pfms.txt"
IDX={"A":0,"C":1,"G":2,"T":3}
def log(*a): print("[79]",*a,flush=True)

def read_jaspar(p):
    mo={}; L=[l.rstrip("\n") for l in open(p)]; i=0
    while i<len(L):
        if not L[i].startswith(">"): i+=1; continue
        h=L[i][1:].split(); mid,name=h[0],(h[1] if len(h)>1 else h[0]); mat=[]
        for j in range(1,5):
            r=L[i+j]; r=r[r.find("[")+1:r.rfind("]")] if "[" in r else r.split(None,1)[1]
            mat.append([float(x) for x in r.split()])
        mo[name.upper()]=np.array(mat); i+=5
    return mo
def pwm(p,ps=.25):
    c=p+ps; q=c/c.sum(0,keepdims=True); return np.log2(q/0.25)
def enc(s):
    a=np.full(len(s),-1,np.int8); b=np.frombuffer(s.encode(),dtype="S1")
    for ch,k in IDX.items(): a[b==ch.encode()]=k
    return a
def hit(idx,P,cut):
    L=P.shape[1]
    if len(idx)<L: return False
    w=np.lib.stride_tricks.sliding_window_view(idx,L)
    v=(w>=0).all(1)
    if not v.any(): return False
    return bool(P[w[v],np.arange(L)].sum(1).max()>=cut)
_G={}
def _init(seqs,mos,thr): _G.update(seqs=seqs,mos=mos,thr=thr)
def _work(k):
    P=_G["mos"][k]; Pr=P[::-1,::-1]; cut=_G["thr"]*P.max(0).sum()
    return k, np.array([hit(ix,P,cut) or hit(ix,Pr,cut) for ix in _G["seqs"]])

if __name__=="__main__":
    import pyfaidx
    fa=pyfaidx.Fasta(FASTA, as_raw=True, sequence_always_upper=True)
    M=pd.read_csv(f"{RES}/52_rewiring_annotated.csv.gz")
    sel=M[M.gained|M.lost].copy()
    peaks=sel[["peak_id","chrom","start","end"]].drop_duplicates().reset_index(drop=True)
    log(f"要掃 motif 的 peak：{len(peaks)}（新生 {M.gained.sum()} + 斷裂 {M.lost.sum()} 組連結）")
    seqs=[]
    keep=[]
    for i,r in peaks.iterrows():
        try: s=str(fa[r.chrom][int(r.start):int(r.end)])
        except Exception: continue
        if len(s)>=60: seqs.append(enc(s)); keep.append(i)
    peaks=peaks.iloc[keep].reset_index(drop=True)
    log(f"實際可用 {len(seqs)} 個 peak")
    mos={k:pwm(v) for k,v in read_jaspar(PFM).items()}
    log(f"motif {len(mos)} 個，開始掃描 ...")
    with Pool(96, initializer=_init, initargs=(seqs,mos,0.85)) as p:
        res=dict(p.map(_work, list(mos.keys()), chunksize=4))
    log("掃描完成")
    H=pd.DataFrame(res, index=peaks.peak_id.values)
    H.to_csv(f"{RES}/79_peak_motif_hits.csv.gz", compression="gzip")

    # TF 必須在腫瘤中有表現
    ex=pd.read_csv(f"{RES}/15_gene_residual_tumor_specific.csv").set_index("gene")
    expressed=set(ex[ex.mean_observed>0.05].index)
    log(f"腫瘤中有表現的基因 {len(expressed)}")

    rows=[]
    for tf in H.columns:
        base=tf.split("::")[0].upper()
        if base not in expressed: continue
        pk=set(H.index[H[tf].values])
        g=sel[sel.peak_id.isin(pk)]
        ng=int(g.gained.sum()); nl=int(g.lost.sum())
        if ng+nl==0: continue
        rows.append(dict(tf=tf, gene=base, n_gained_links=ng, n_lost_links=nl,
                         n_gained_genes=g[g.gained].gene.nunique(),
                         n_lost_genes=g[g.lost].gene.nunique(),
                         net=ng-nl,
                         tf_resid=float(ex.loc[base,"resid_tumor_specific"]) if base in ex.index else np.nan))
    R=pd.DataFrame(rows)
    # 全體基準：每個 TF 的 peak 中有多少比例是 gained
    tot_g=int(sel.gained.sum()); tot_l=int(sel.lost.sum())
    R["gained_frac"]=R.n_gained_links/(R.n_gained_links+R.n_lost_links)
    base_frac=tot_g/(tot_g+tot_l)
    from scipy.stats import binomtest
    R["p_enrich"]=[binomtest(int(r.n_gained_links), int(r.n_gained_links+r.n_lost_links),
                             base_frac, alternative="greater").pvalue for _,r in R.iterrows()]
    R["padj"]=np.minimum(1, R.p_enrich*len(R))
    R=R.sort_values("p_enrich")
    R.to_csv(f"{RES}/79_eregulon_rewiring.csv", index=False)
    pd.set_option("display.width",230)
    log(f"\n全體基準：新生佔 {100*base_frac:.1f}%")
    log("\n=== 其結合位點最偏向『新生連結』的 TF（= 在腫瘤中擴張地盤）Top 25 ===")
    log("\n"+R.head(25)[["tf","n_gained_links","n_lost_links","gained_frac","n_gained_genes","tf_resid","padj"]].round(4).to_string(index=False))
    log("\n=== 最偏向『斷裂連結』的 TF（= 在腫瘤中失守）Top 15 ===")
    log("\n"+R.sort_values("gained_frac").head(15)[["tf","n_gained_links","n_lost_links","gained_frac","n_lost_genes","tf_resid"]].round(4).to_string(index=False))
    log("DONE")
