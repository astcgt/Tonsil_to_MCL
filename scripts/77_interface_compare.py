#!/usr/bin/env python
"""
步驟 77：比較「設計結合劑」與「天然 NRF2」各自接觸 MafG 的哪些殘基。

這決定了設計品是不是真的打在「零逃逸」的位置上。
從 Boltz-2 預測出來的複合體結構直接算接觸。
"""
import os, sys, gzip, glob, json
from collections import defaultdict
import numpy as np, pandas as pd

BASE="/raid/jupyter-alz.ee09/Tonsil_to_MCL"; RES=f"{BASE}/results"
PRED=f"{RES}/72_boltz_calibrated/boltz_results_inputs/predictions"
CUT=5.0; MAFG_OFFSET=23       # 鏈 A 的第 1 個殘基 = MafG 的 23
AA3={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I',
     'LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
def log(*a): print("[77]",*a,flush=True)

def read_cif(p):
    op=gzip.open if p.endswith(".gz") else open
    hdr,rows=[],[]
    with op(p,"rt") as fh:
        for l in fh:
            if l.startswith("_atom_site."): hdr.append(l.strip().split(".",1)[1])
            elif l.startswith(("ATOM","HETATM")):
                q=l.split()
                if len(q)>=len(hdr): rows.append(q[:len(hdr)])
    if not rows: return None
    d=pd.DataFrame(rows, columns=hdr)
    for c in ("Cartn_x","Cartn_y","Cartn_z"): d[c]=pd.to_numeric(d[c],errors="coerce")
    d["auth_seq_id"]=pd.to_numeric(d["auth_seq_id"],errors="coerce")
    return d.dropna(subset=["Cartn_x","auth_seq_id"])

esc=pd.read_csv(f"{RES}/75_escape_map_MAFG.csv")
zero=set(esc[esc.n_tolerated==0].pos.astype(int))
tol =set(esc[esc.n_tolerated>0].pos.astype(int))
log(f"MafG 介面：零逃逸殘基 {len(zero)} 個、可容忍殘基 {len(tol)} 個")

rows=[]
for d in sorted(glob.glob(f"{PRED}/*")):
    nm=os.path.basename(d)
    cifs=glob.glob(os.path.join(d,"*.cif"))+glob.glob(os.path.join(d,"*.cif.gz"))
    if not cifs: continue
    df=read_cif(sorted(cifs)[0])
    if df is None: continue
    ch="auth_asym_id" if "auth_asym_id" in df else "label_asym_id"
    chains=df[ch].unique()
    if len(chains)<2: continue
    A=df[df[ch]==chains[0]]; B=df[df[ch]==chains[1]]
    PA=A[["Cartn_x","Cartn_y","Cartn_z"]].values; PB=B[["Cartn_x","Cartn_y","Cartn_z"]].values
    dmat=np.sqrt(((PA[:,None,:]-PB[None,:,:])**2).sum(-1))
    near=dmat.min(axis=1)
    res=A["auth_seq_id"].values.astype(int)
    contact={}
    for i,v in enumerate(near):
        if v<=CUT:
            p=res[i]+MAFG_OFFSET-1
            if p not in contact or v<contact[p]: contact[p]=v
    nz=len(set(contact)&zero); nt=len(set(contact)&tol)
    rows.append(dict(name=nm, n_contact=len(contact),
                     zero_escape_hit=nz, tolerant_hit=nt,
                     zero_frac=round(nz/max(1,len(contact)),3),
                     residues=" ".join(f"{p}" for p in sorted(contact))))
R=pd.DataFrame(rows).sort_values("zero_escape_hit", ascending=False)
R.to_csv(f"{RES}/77_interface_compare.csv", index=False)
pd.set_option("display.width",250); pd.set_option("display.max_colwidth",70)
log("\n=== 各複合體接觸到的 MafG 殘基（零逃逸 vs 可容忍）===")
log("\n"+R[["name","n_contact","zero_escape_hit","tolerant_hit","zero_frac"]].to_string(index=False))
log("\n=== 接觸殘基明細 ===")
for _,r in R.iterrows():
    log(f"  {r['name']}\n     {r.residues}")
nat=R[R.name.str.contains("CTRL_POS")]
if len(nat):
    log(f"\n天然 NRF2 打中 {int(nat.zero_escape_hit.iloc[0])} 個零逃逸殘基")
    des=R[~R.name.str.contains("CTRL")]
    log(f"設計品打中零逃逸殘基的範圍：{des.zero_escape_hit.min()}–{des.zero_escape_hit.max()}；"
        f"超過天然的有 {int((des.zero_escape_hit>nat.zero_escape_hit.iloc[0]).sum())}/{len(des)} 個")
log("DONE")
