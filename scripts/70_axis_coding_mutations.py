#!/usr/bin/env python
"""
步驟 70：找出抗氧化／鐵死亡防禦軸線上的「編碼區」突變。

三個來源
--------
1. 22 位 MCL 病人的全基因組定序（Zenodo 16098829，249,321 個體細胞 SNV）
2. DepMap 八株 MCL 細胞株的體細胞突變
3. cBioPortal mcl_idibips_2013（29 例）

把 WGS 的 SNV 對到 GENCODE 的 CDS 上，翻譯出胺基酸變化，
作為後續 ESM / 結構擾動實驗的輸入。
特別關注 KEAP1 / NFE2L2 / CUL3 —— 若 KEAP1 失能或 NRF2 的 degron 突變，
就能一舉解釋整條軸線為何被打開。
"""
import os, sys, gzip, json
from collections import defaultdict
import numpy as np, pandas as pd

BASE="/raid/jupyter-alz.ee09/Tonsil_to_MCL"; RES=f"{BASE}/results"; GEN=f"{BASE}/data/genome"
GTF=f"{GEN}/gencode.v44.annotation.gtf.gz"; FASTA=f"{GEN}/GRCh38.primary_assembly.genome.fa"
AXIS=["SLC7A11","SLC3A2","GCLM","GCLC","AIFM2","GPX4","MAFG","MAFF","MAFK",
      "BACH1","BACH2","NFE2L2","NFE2L1","KEAP1","CUL3","RBX1","TXNRD1","GSR","OXR1","FTH1","SLC40A1"]
CODON={}
B="TCAG"
AAs="FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
i=0
for a in B:
    for b in B:
        for c in B:
            CODON[a+b+c]=AAs[i]; i+=1
def log(*a): print("[70]",*a,flush=True)

def load_cds(gtf, genes):
    want=set(genes); out=defaultdict(list); tx_of={}
    with gzip.open(gtf,"rt") as fh:
        for line in fh:
            if line[0]=="#": continue
            f=line.rstrip("\n").split("\t")
            if f[2] not in ("CDS",): continue
            a=f[8]
            if 'gene_name "' not in a: continue
            g=a.split('gene_name "')[1].split('"')[0]
            if g not in want: continue
            if 'tag "MANE_Select"' not in a and 'tag "Ensembl_canonical"' not in a: continue
            t=a.split('transcript_id "')[1].split('"')[0]
            out[(g,t)].append((f[0], int(f[3]), int(f[4]), f[6], int(f[7])))
            tx_of[g]=t
    # 每個基因只留一個轉錄本（優先 MANE）
    cds={}
    for (g,t),v in out.items():
        if g in cds and len(cds[g][1])>=len(v): continue
        cds[g]=(t, sorted(v, key=lambda x:x[1]))
    return cds

cds=load_cds(GTF, AXIS)
log("取得 CDS 的基因:", sorted(cds))
import pyfaidx
fa=pyfaidx.Fasta(FASTA, as_raw=True, sequence_always_upper=True)
COMP=str.maketrans("ACGTN","TGCAN")

def build_prot(g):
    t,ex=cds[g]
    strand=ex[0][3]; chrom=ex[0][0]
    pieces=[]; coords=[]
    for (c,s,e,st,ph) in ex:
        seq=str(fa[c][s-1:e]); pieces.append(seq); coords += list(range(s,e+1))
    nt="".join(pieces)
    if strand=="-":
        nt=nt.translate(COMP)[::-1]; coords=coords[::-1]
    prot="".join(CODON.get(nt[i:i+3],"X") for i in range(0,len(nt)-2,3))
    return dict(tx=t, chrom=chrom, strand=strand, nt=nt, coords=coords, prot=prot)

P={g:build_prot(g) for g in cds}
for g,d in sorted(P.items()):
    log(f"  {g:9s} {d['tx']}  {len(d['prot'])} aa  {'(含終止碼)' if d['prot'].endswith('*') else ''}")
json.dump({g:{k:v for k,v in d.items() if k in ("tx","chrom","strand","prot")} for g,d in P.items()},
          open(f"{RES}/70_axis_proteins.json","w"), indent=1)

# ---------- WGS ----------
mut=pd.read_csv(f"{BASE}/data/external_mcl/wgs/MCL.somaticMutations.tsv", sep="\t")
mut.columns=["chrom","start","end","ref","alt","sample"]
mut["patient"]=mut["sample"].str.replace(r"[PR]\d*$","",regex=True)
rows=[]
for g,d in P.items():
    pos2i={p:i for i,p in enumerate(d["coords"])}
    sub=mut[(mut.chrom==d["chrom"]) & mut.start.isin(pos2i)]
    for _,r in sub.iterrows():
        i=pos2i[r.start]
        ref=r.ref; alt=r.alt
        if d["strand"]=="-":
            ref=ref.translate(COMP); alt=alt.translate(COMP)
        if d["nt"][i]!=ref:
            continue
        ci=i//3; off=i%3
        cod=d["nt"][ci*3:ci*3+3]
        if len(cod)<3: continue
        new=cod[:off]+alt+cod[off+1:]
        wt=CODON.get(cod,"X"); mt=CODON.get(new,"X")
        rows.append(dict(source="WGS_22pt", gene=g, patient=r.patient, sample=r["sample"],
                         chrom=r.chrom, pos=int(r.start), ref=r.ref, alt=r.alt,
                         aa_pos=ci+1, aa_ref=wt, aa_alt=mt,
                         kind="synonymous" if wt==mt else ("nonsense" if mt=="*" else "missense")))
W=pd.DataFrame(rows)
log(f"\n=== WGS：軸線基因編碼區突變 {len(W)} 個 ===")
if len(W):
    pd.set_option("display.width",220)
    log("\n"+W[["gene","patient","chrom","pos","ref","alt","aa_pos","aa_ref","aa_alt","kind"]].to_string(index=False))
    log("\n依基因統計:\n"+W.groupby(["gene","kind"]).size().to_string())

# ---------- DepMap ----------
dm=pd.read_csv(f"{BASE}/results/09_depmap_mcl_mutations.csv", low_memory=False)
mo=pd.read_csv(f"{BASE}/data/depmap/Model.csv", low_memory=False).set_index("ModelID")
d2=dm[dm.HugoSymbol.isin(AXIS)].copy()
if len(d2): d2["line"]=d2.ModelID.map(mo.CellLineName)
log(f"\n=== DepMap 八株 MCL 細胞株：軸線基因突變 {len(d2)} 個 ===")
if len(d2):
    cols=[c for c in ["line","HugoSymbol","ProteinChange","VariantInfo","LikelyLoF"] if c in d2.columns]
    log("\n"+d2[cols].to_string(index=False))

# ---------- cBioPortal ----------
try:
    cb=json.load(open(f"{BASE}/data/aux/cbio_mcl_mut.json"))
    c2=[x for x in cb if (x.get("gene") or {}).get("hugoGeneSymbol") in AXIS]
    log(f"\n=== cBioPortal MCL (29 例)：軸線基因突變 {len(c2)} 個 ===")
    for x in c2:
        log(f"  {x['sampleId']}  {(x.get('gene') or {}).get('hugoGeneSymbol')}  {x.get('proteinChange')}  {x.get('mutationType')}")
except Exception as e:
    log("cBioPortal 讀取失敗:", e)

allm=[]
if len(W): allm.append(W)
if len(d2):
    allm.append(pd.DataFrame(dict(source="DepMap_MCL_lines", gene=d2.HugoSymbol,
                                  patient=d2.get("line"), aa_change=d2.get("ProteinChange"),
                                  kind=d2.get("VariantInfo"))))
if allm:
    pd.concat(allm, ignore_index=True).to_csv(f"{RES}/70_axis_coding_mutations.csv", index=False)
    log(f"\n寫出 {RES}/70_axis_coding_mutations.csv")
log("DONE")
