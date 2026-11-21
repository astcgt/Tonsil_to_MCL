#!/usr/bin/env python
"""步驟 82：通用的 RFdiffusion3 結合劑評分器（支援多個標的）。"""
import os, sys, gzip, glob, json, argparse
import numpy as np, pandas as pd
BASE="/raid/jupyter-alz.ee09/Tonsil_to_MCL"; RES=f"{BASE}/results"
AA3={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I',
     'LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
HYDRO=set("AVLIMFWYC"); CUT=5.0
def log(*a): print("[82]",*a,flush=True)

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

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--cifdir", required=True); ap.add_argument("--n_target", type=int, required=True)
    ap.add_argument("--pdb_first_res", type=int, required=True)
    ap.add_argument("--hotspots", required=True, help="逗號分隔的 PDB 殘基編號")
    ap.add_argument("--prefix", default=None); ap.add_argument("--out", required=True)
    a=ap.parse_args()
    HOT=set(int(x) for x in a.hotspots.split(","))
    hot_renum={h-a.pdb_first_res+1 for h in HOT}
    rows=[]
    pat=f"{a.cifdir}/*{a.prefix}*.cif.gz" if a.prefix else f"{a.cifdir}/*.cif.gz"
    for cif in sorted(glob.glob(pat)):
        j=cif.replace(".cif.gz",".json")
        meta=json.load(open(j)) if os.path.exists(j) else {}
        df=read_cif(cif)
        if df is None: continue
        rs=df["auth_seq_id"].values
        rn="auth_comp_id" if "auth_comp_id" in df else "label_comp_id"
        T=df[rs<=a.n_target]; B=df[rs>a.n_target]
        if len(T)==0 or len(B)==0: continue
        Tx=T[["Cartn_x","Cartn_y","Cartn_z"]].values; Bx=B[["Cartn_x","Cartn_y","Cartn_z"]].values
        d=np.sqrt(((Tx[:,None,:]-Bx[None,:,:])**2).sum(-1)); near=d<=CUT
        tr=T["auth_seq_id"].values.astype(int); br=B["auth_seq_id"].values.astype(int)
        ct=set(tr[near.any(1)]); cb=set(br[near.any(0)])
        ii,jj=np.where(near)
        pairs=len({(int(x),int(y)) for x,y in zip(tr[ii],br[jj])})
        hit=len(ct&hot_renum)
        sm=(B.drop_duplicates(subset=["auth_seq_id"]).sort_values("auth_seq_id")[["auth_seq_id",rn]].values)
        seq="".join(AA3.get(r,"X") for _,r in sm)
        bres={int(n):AA3.get(r,"X") for n,r in sm}
        ifa=[bres.get(i,"X") for i in cb]
        hyd=np.mean([x in HYDRO for x in ifa]) if ifa else np.nan
        m=meta.get("metrics",{})
        rows.append(dict(design=os.path.basename(cif).replace(".cif.gz","").replace("input_",""),
            binder_len=len(seq), hotspot_covered=hit, hotspot_total=len(hot_renum),
            hotspot_frac=round(hit/max(1,len(hot_renum)),3),
            iface_pairs=pairs, iface_binder_res=len(cb),
            iface_hydrophobic_frac=round(float(hyd),3) if hyd==hyd else np.nan,
            clashes_sc=m.get("n_clashing.interresidue_clashes_w_sidechain"),
            clashes_bb=m.get("n_clashing.interresidue_clashes_w_backbone"),
            helix=m.get("helix_fraction"), sheet=m.get("sheet_fraction"),
            ala_content=m.get("alanine_content"), sequence=seq))
    if not rows: log("沒有可評分的設計"); return
    R=pd.DataFrame(rows)
    R["score"]=(R.hotspot_frac*100 + R.iface_pairs.clip(upper=140)*0.25
                - R.clashes_sc.fillna(0)*3 - R.clashes_bb.fillna(0)*6
                - (R.ala_content.fillna(0).clip(lower=0.15)-0.15)*200)
    R=R.sort_values("score",ascending=False)
    R.to_csv(a.out, index=False)
    pd.set_option("display.width",240)
    log(f"\n=== {a.prefix or '全部'}：{len(R)} 個設計，前 10 ===")
    log("\n"+R.head(10)[["design","binder_len","hotspot_covered","hotspot_total","iface_pairs",
                          "iface_hydrophobic_frac","clashes_sc","ala_content","score"]].round(3).to_string(index=False))
    log(f"寫出 {a.out}")
main()
