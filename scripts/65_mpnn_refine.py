#!/usr/bin/env python
"""
步驟 65：把 RFdiffusion3 的骨架用 ProteinMPNN 重新設計序列。

RFdiffusion3 自帶的 sequence head 會給出富含丙胺酸的低複雜度序列，
不適合直接合成。ProteinMPNN 以骨架幾何為條件重新取樣序列，
產出的胺基酸組成與天然蛋白接近。

做法：標的鏈的序列固定不動，只重新設計結合劑那一段。
"""
import os, sys, gzip, glob, subprocess, argparse, json
import numpy as np, pandas as pd

BASE="/raid/jupyter-alz.ee09/Tonsil_to_MCL"; RES=f"{BASE}/results"
def log(*a): print("[65]",*a,flush=True)

def cif_to_pdb(cif, pdb, n_target):
    """把 RFdiffusion3 的 mmCIF 轉成兩鏈 PDB：A = 標的，B = 結合劑。"""
    op = gzip.open if cif.endswith(".gz") else open
    hdr, rows = [], []
    with op(cif,"rt") as fh:
        for line in fh:
            if line.startswith("_atom_site."): hdr.append(line.strip().split(".",1)[1])
            elif line.startswith(("ATOM","HETATM")):
                p=line.split()
                if len(p)>=len(hdr): rows.append(p[:len(hdr)])
    if not rows: return 0
    df=pd.DataFrame(rows, columns=hdr)
    ci={c:i for i,c in enumerate(hdr)}
    n=0
    with open(pdb,"w") as fh:
        serial=0
        for _,r in df.iterrows():
            try: rs=int(float(r["auth_seq_id"]))
            except Exception: continue
            ch = "A" if rs<=n_target else "B"
            num = rs if rs<=n_target else rs-n_target
            el = r["type_symbol"]; an=r["auth_atom_id"].strip('"'); rn=r["auth_comp_id"]
            try: x,y,z=float(r["Cartn_x"]),float(r["Cartn_y"]),float(r["Cartn_z"])
            except Exception: continue
            serial+=1; n+=1
            name = an if len(an)>=4 else f" {an:<3}"
            fh.write(f"ATOM  {serial:5d} {name:4s}{'':1s}{rn:>3s} {ch}{num:4d}{'':4s}"
                     f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{0.00:6.2f}{'':10s}{el:>2s}\n")
        fh.write("TER\nEND\n")
    return n

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--ranking", required=True)
    ap.add_argument("--cifdir", required=True)
    ap.add_argument("--n_target", type=int, required=True)
    ap.add_argument("--top", type=int, default=8)
    ap.add_argument("--nseq", type=int, default=8)
    ap.add_argument("--gpu", default="3")
    ap.add_argument("--tag", required=True)
    a=ap.parse_args()
    out=f"{RES}/65_mpnn_{a.tag}"; os.makedirs(out, exist_ok=True)
    pdbdir=os.path.join(out,"pdb"); os.makedirs(pdbdir, exist_ok=True)
    env=dict(os.environ, CUDA_VISIBLE_DEVICES=a.gpu)

    df=pd.read_csv(a.ranking).head(a.top)
    log(f"處理前 {len(df)} 個設計（標的 {a.n_target} 殘基）")
    recs=[]
    for _,r in df.iterrows():
        nm=str(r.design)
        cif=os.path.join(a.cifdir, f"input_{nm}.cif.gz")
        if not os.path.exists(cif):
            cand=glob.glob(os.path.join(a.cifdir, f"*{nm}*.cif.gz"))
            if not cand: log(f"  找不到 {nm}"); continue
            cif=cand[0]
        pdb=os.path.join(pdbdir, f"{nm}.pdb")
        n=cif_to_pdb(cif, pdb, a.n_target)
        if n==0: log(f"  {nm} 轉檔失敗"); continue
        cmd=[f"{BASE}/env/protdes/bin/mpnn","--model_type","protein_mpnn",
             "--checkpoint_path",f"{BASE}/models/foundry/proteinmpnn_v_48_020.pt",
             "--is_legacy_weights","True","--structure_path",pdb,
             "--out_directory",out,"--name",nm,
             "--designed_chains","B",
             "--temperature","0.1","--batch_size",str(a.nseq),
             "--number_of_batches","1","--write_fasta","True","--seed","0"]
        p=subprocess.run(cmd, env=env, capture_output=True, text=True)
        fas=glob.glob(os.path.join(out, f"**/{nm}*.fa*"), recursive=True)
        if not fas:
            log(f"  {nm}: MPNN 無輸出（rc={p.returncode}）"); 
            if p.returncode!=0: log("    "+p.stderr[-400:].replace("\n"," | "))
            continue
        seqs=[]
        for f in fas:
            cur=None
            for line in open(f):
                line=line.strip()
                if line.startswith(">"): cur=line
                elif line and cur: seqs.append(line); cur=None
        for i,s in enumerate(seqs[:a.nseq]):
            recs.append(dict(design=nm, mpnn_idx=i, sequence=s, length=len(s),
                             hotspot_covered=r.get("hotspot_covered"),
                             iface_pairs=r.get("iface_pairs"),
                             rfd3_sequence=r.get("sequence")))
        log(f"  {nm}: 取得 {len(seqs)} 條序列")
    if recs:
        R=pd.DataFrame(recs)
        def comp(s):
            return {aa: round(s.count(aa)/len(s),3) for aa in "AGLVIMFWYKRDEQNSTHCP"}
        R["ala_frac"]=R.sequence.map(lambda s: round(s.count("A")/len(s),3))
        R["charge"]=R.sequence.map(lambda s: s.count("K")+s.count("R")-s.count("D")-s.count("E"))
        R["hydrophobic_frac"]=R.sequence.map(lambda s: round(sum(c in "AVLIMFWYC" for c in s)/len(s),3))
        R.to_csv(f"{RES}/65_mpnn_{a.tag}.csv", index=False)
        pd.set_option("display.width",240); pd.set_option("display.max_colwidth",120)
        log(f"\n=== ProteinMPNN 精修後（{a.tag}），每個骨架取第 1 條 ===")
        log("\n"+R[R.mpnn_idx==0][["design","length","ala_frac","hydrophobic_frac","charge","sequence"]].to_string(index=False))
        log(f"\n原始 RFdiffusion 序列的丙胺酸比例平均 "
            f"{np.mean([s.count('A')/len(s) for s in R.rfd3_sequence.dropna().unique()]):.3f}"
            f" → 精修後 {R.ala_frac.mean():.3f}")
        log(f"寫出 {RES}/65_mpnn_{a.tag}.csv（共 {len(R)} 條序列）")
main()
