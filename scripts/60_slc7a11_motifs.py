#!/usr/bin/env python
"""
步驟 60：SLC7A11 新生調控域的 motif 組成 —— 是誰把它打開的？

把 SLC7A11 的 5 個新生元件當前景，
用「長度 + GC 配對」的背景（同樣在 MCL 中可及、但沒有新生連結的 peak）做富集檢定。
前景只有 5 個序列、統計力有限，所以同時列出每個元件各自的最佳命中，
讓證據以「哪個 TF 在每個元件都出現」的形式呈現，而不是只看 p 值。
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
from scipy.stats import hypergeom
from multiprocessing import Pool
sys.path.insert(0, f"{os.path.dirname(os.path.abspath(__file__))}")
from importlib import import_module
m27 = import_module("27_motif_enrichment") if False else None

BASE="/raid/jupyter-alz.ee09/Tonsil_to_MCL"; RES=f"{BASE}/results"; GEN=f"{BASE}/data/genome"
FASTA=f"{GEN}/GRCh38.primary_assembly.genome.fa"; PFM=f"{GEN}/JASPAR2024_CORE_vertebrates_nr_pfms.txt"
IDX={"A":0,"C":1,"G":2,"T":3}
def log(*a): print("[60]",*a,flush=True)

def read_jaspar(path):
    mo={}; L=[l.rstrip("\n") for l in open(path)]; i=0
    while i<len(L):
        if not L[i].startswith(">"): i+=1; continue
        h=L[i][1:].split(); mid,name=h[0],(h[1] if len(h)>1 else h[0]); mat=[]
        for j in range(1,5):
            r=L[i+j]; r=r[r.find("[")+1:r.rfind("]")] if "[" in r else r.split(None,1)[1]
            mat.append([float(x) for x in r.split()])
        mo[f"{mid}_{name}"]=np.array(mat); i+=5
    return mo
def pwm(p,ps=.25):
    c=p+ps; q=c/c.sum(0,keepdims=True); return np.log2(q/0.25)
def enc(s):
    a=np.full(len(s),-1,np.int8)
    b=np.frombuffer(s.encode(),dtype="S1")
    for ch,k in IDX.items(): a[b==ch.encode()]=k
    return a
def best(idx,P):
    L=P.shape[1]
    if len(idx)<L: return -1e9,-1
    w=np.lib.stride_tricks.sliding_window_view(idx,L)
    v=(w>=0).all(1)
    if not v.any(): return -1e9,-1
    sc=P[w[v],np.arange(L)].sum(1)
    k=int(sc.argmax()); return float(sc[k]), int(np.flatnonzero(v)[k])
_G={}
def _init(seqs,mos,thr): _G.update(seqs=seqs,mos=mos,thr=thr)
def _work(k):
    P=_G["mos"][k]; Pr=P[::-1,::-1]; cut=_G["thr"]*P.max(0).sum()
    hits=np.zeros(len(_G["seqs"]),bool); pos=[-1]*len(_G["seqs"])
    for i,ix in enumerate(_G["seqs"]):
        s1,p1=best(ix,P); s2,p2=best(ix,Pr)
        if max(s1,s2)>=cut: hits[i]=True; pos[i]=p1 if s1>=s2 else p2
    return k,hits,pos

if __name__=="__main__":
    import pyfaidx
    fa=pyfaidx.Fasta(FASTA,as_raw=True,sequence_always_upper=True)
    M=pd.read_csv(f"{RES}/52_rewiring_annotated.csv.gz")
    atac=pd.read_csv(f"{RES}/22_atac_diff.csv")
    fg=M[(M.gained)&(M.gene=="SLC7A11")][["chrom","start","end","dist","r_t"]].copy()
    fg["name"]=fg.apply(lambda r:f"{r.chrom}:{int(r.start)}-{int(r.end)} ({int(r.dist):+,} bp)",axis=1)
    log(f"前景：SLC7A11 的 {len(fg)} 個新生元件")
    # 背景：在 MCL 中可及、未被判為新生/斷裂、長度與 GC 配對
    used=set(M.loc[M.gained|M.lost,"peak_id"])
    cand=M[(~M.peak_id.isin(used))][["chrom","start","end"]].drop_duplicates()
    cand["len"]=cand.end-cand.start
    rng=np.random.default_rng(0)
    def gc(s): return (s.count("G")+s.count("C"))/max(1,len(s))
    def grab(df):
        out=[]
        for _,r in df.iterrows():
            try: out.append(str(fa[r.chrom][int(r.start):int(r.end)]))
            except Exception: out.append("")
        return out
    fgs=grab(fg); fg_len=[len(s) for s in fgs]; fg_gc=[gc(s) for s in fgs]
    pool=cand.sample(n=min(40000,len(cand)),random_state=1).reset_index(drop=True)
    pool["seq"]=grab(pool); pool=pool[pool.seq.str.len()>50]
    pool["gc"]=pool.seq.map(gc)
    bg_idx=[]
    for L,G in zip(fg_len,fg_gc):
        ok=(np.abs(pool["len"].values-L)<=0.3*L)&(np.abs(pool.gc.values-G)<=0.05)
        c=[i for i in np.flatnonzero(ok) if i not in bg_idx]
        if len(c)<200:
            ok=(np.abs(pool["len"].values-L)<=0.6*L)&(np.abs(pool.gc.values-G)<=0.10)
            c=[i for i in np.flatnonzero(ok) if i not in bg_idx]
        bg_idx += list(rng.choice(c,size=min(300,len(c)),replace=False))
    bgs=pool.seq.iloc[bg_idx].tolist()
    log(f"配對背景 {len(bgs)} 條（前景長度 {fg_len}，GC {[round(x,2) for x in fg_gc]}）")

    mos={k:pwm(v) for k,v in read_jaspar(PFM).items()}
    log("motif 數:",len(mos))
    res={}
    for tag,seqs in [("fg",fgs),("bg",bgs)]:
        ix=[enc(s) for s in seqs]
        with Pool(96,initializer=_init,initargs=(ix,mos,0.85)) as p:
            r=p.map(_work,list(mos.keys()),chunksize=4)
        res[tag]={k:(h,po) for k,h,po in r}
        log(f"  {tag} 掃描完成")
    rows=[]
    nf,nb=len(fgs),len(bgs)
    for k in mos:
        a=int(res["fg"][k][0].sum()); b=int(res["bg"][k][0].sum())
        if a==0: continue
        pv=hypergeom.sf(a-1,nf+nb,a+b,nf)
        rows.append(dict(motif=k,fg_hits=a,fg_n=nf,bg_rate=round(b/nb,4),
                         fold=round((a/nf)/max(b/nb,1e-6),2),pvalue=pv,
                         elements="".join("1" if x else "0" for x in res["fg"][k][0])))
    R=pd.DataFrame(rows).sort_values(["fg_hits","pvalue"],ascending=[False,True])
    R.to_csv(f"{RES}/60_slc7a11_motifs.csv",index=False)
    pd.set_option("display.width",220)
    log("\n=== 在全部 5 個 SLC7A11 新生元件都出現的 motif ===")
    log("\n"+R[R.fg_hits==5].head(30)[["motif","bg_rate","fold","pvalue","elements"]].to_string(index=False))
    log("\n=== 在 4 個元件出現、且富集最強的 ===")
    log("\n"+R[R.fg_hits==4].sort_values("pvalue").head(20)[["motif","bg_rate","fold","pvalue","elements"]].to_string(index=False))
    log("\n（elements 欄位的 5 位元對應下列元件順序）")
    for i,n in enumerate(fg.name): log(f"   位元{i+1}: {n}")
    log("DONE")
