#!/usr/bin/env python
"""
步驟 88：用「晶體結構當模板」重新評估 FSP1 結合劑。

為什麼需要
----------
前一次評估中，FSP1（373 aa 酵素）在無 MSA 的單序列模式下自己就摺不好
（pLDDT 0.46），所以複合體的 ipTM（0.13–0.20）量的是「模型不會摺 FSP1」，
而不是「結合劑不好」。這是技術性失敗，不是設計失敗。

解法：Boltz-2 支援 `templates:` 欄位。把 PDB 8YOQ 的 FSP1 鏈當模板餵進去，
標的摺疊就被固定住，模型只需要把結合劑擺上去 —— 這才是我們要問的問題。

同時加入對照：打散序列的結合劑（同樣給模板）。
"""
import os, sys, json, glob, subprocess, argparse, random, gzip
import numpy as np, pandas as pd
BASE="/raid/jupyter-alz.ee09/Tonsil_to_MCL"; RES=f"{BASE}/results"; ST=f"{BASE}/data/structures"
FSP1=("MGSQVSVESGALHVVIVGGGFGGIAAASQLQALNVPFMLVDMKDSFHHNVAALRASVETGFAKKTFISYSVTFKDNFRQGLVVGIDLKNQMVLLQGG"
      "EALPFSHLILATGSTGPFPGKFNEVSSQQAAIQAYEDMVRQVQRSRFIVVVGGGSAGVEMAAEIKTEYPEKEVTLIHSQVALADKELLPSVRQEVKEI"
      "LLRKGVQLLLSERVSNLEELPLNEYREYIKVQTDKGTEVATNLVILCTGIKINSSAYRKAFESRLASSGALRVNEHLQVEGHSNVYAIGDCADVRTPK"
      "MAYLAGLHANIAVANIVNSVKQRPLQAYKPGALTFLLSMGRNDGVGQISGFYVGRLMVRLTKSRDLFVSKSWKTMRQSPP")
def log(*a): print("[88]",*a,flush=True)

def pdb_to_cif(pdb, out, chain="A"):
    """把單鏈 PDB 轉成最小 mmCIF（Boltz 模板用）。"""
    rows=[]; n=0
    for l in open(pdb):
        if not l.startswith("ATOM"): continue
        if l[21]!=chain: continue
        n+=1
        rows.append(dict(group="ATOM", id=n, type_symbol=l[76:78].strip() or l[12:16].strip()[0],
            label_atom_id=l[12:16].strip(), label_comp_id=l[17:20].strip(),
            label_asym_id="A", label_seq_id=int(l[22:26]),
            x=float(l[30:38]), y=float(l[38:46]), z=float(l[46:54])))
    with open(out,"w") as fh:
        fh.write("data_template\n#\nloop_\n")
        for c in ["_atom_site.group_PDB","_atom_site.id","_atom_site.type_symbol",
                  "_atom_site.label_atom_id","_atom_site.label_alt_id","_atom_site.label_comp_id",
                  "_atom_site.label_asym_id","_atom_site.label_entity_id","_atom_site.label_seq_id",
                  "_atom_site.pdbx_PDB_ins_code","_atom_site.Cartn_x","_atom_site.Cartn_y",
                  "_atom_site.Cartn_z","_atom_site.occupancy","_atom_site.B_iso_or_equiv",
                  "_atom_site.auth_seq_id","_atom_site.auth_comp_id","_atom_site.auth_asym_id",
                  "_atom_site.auth_atom_id","_atom_site.pdbx_PDB_model_num"]:
            fh.write(c+"\n")
        for r in rows:
            fh.write(f"ATOM {r['id']} {r['type_symbol']} {r['label_atom_id']} . {r['label_comp_id']} A 1 "
                     f"{r['label_seq_id']} ? {r['x']:.3f} {r['y']:.3f} {r['z']:.3f} 1.00 0.00 "
                     f"{r['label_seq_id']} {r['label_comp_id']} A {r['label_atom_id']} 1\n")
        fh.write("#\n")
    return n

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--gpu", default="3"); ap.add_argument("--samples", type=int, default=3)
    a=ap.parse_args()
    out=f"{RES}/88_fsp1_template"; os.makedirs(out, exist_ok=True)
    yml=os.path.join(out,"inputs"); os.makedirs(yml, exist_ok=True)
    for f in glob.glob(os.path.join(yml,"*.yaml")): os.remove(f)
    tcif=os.path.join(out,"fsp1_template.cif")
    n=pdb_to_cif(f"{RES}/62_fsp1_target.pdb", tcif)
    log(f"模板 CIF 寫出：{n} 個原子 → {tcif}")

    B=pd.concat([pd.read_csv(f"{RES}/85_binders_fsp1nadp.csv"),
                 pd.read_csv(f"{RES}/85_binders_fsp1mouth.csv")], ignore_index=True)
    B=B.sort_values("score2",ascending=False).drop_duplicates("design").head(4)
    meta={}
    def write(nm, b):
        with open(os.path.join(yml,f"{nm}.yaml"),"w") as fh:
            fh.write("version: 1\nsequences:\n")
            fh.write(f"  - protein:\n      id: A\n      sequence: {FSP1}\n      msa: empty\n")
            fh.write(f"  - protein:\n      id: B\n      sequence: {b}\n      msa: empty\n")
            fh.write(f"templates:\n  - pdb: {RES}/88_fsp1_template_renum.pdb\n    chain_id: [A]\n")
    for _,r in B.iterrows():
        nm=f"{r.design}_m{int(r.mpnn_idx)}"
        write(nm, r.binder)
        meta[nm]=dict(kind=("NADP 口袋" if "nadp" in r.design else "口袋入口"),
                      blen=int(r.blen), hotspot=int(r.hotspot_covered), seq=r.binder)
    rng=random.Random(2)
    for i,(_,r) in enumerate(B.head(2).iterrows()):
        s=list(r.binder); rng.shuffle(s); s="".join(s)
        write(f"CTRL_NEG_{i+1}", s); meta[f"CTRL_NEG_{i+1}"]=dict(kind="陰性對照(打散)", blen=len(s))
    log(f"{len(meta)} 個複合體 × {a.samples} 次取樣（含模板）")
    env=dict(os.environ, CUDA_VISIBLE_DEVICES=a.gpu)
    cmd=[f"{BASE}/env/boltz/bin/boltz","predict",yml,"--out_dir",out,"--cache",f"{BASE}/models/boltz",
         "--output_format","mmcif","--diffusion_samples",str(a.samples),"--recycling_steps","3",
         "--devices","1","--accelerator","gpu","--no_kernels"]
    p=subprocess.run(cmd, env=env, capture_output=True, text=True)
    log("return code:", p.returncode)
    if p.returncode!=0: log(p.stderr[-2500:])
    rows=[]
    for f in glob.glob(os.path.join(out,"**","confidence_*.json"), recursive=True):
        d=json.load(open(f)); nm=os.path.basename(f)[len("confidence_"):].rsplit("_model_",1)[0]
        rows.append(dict(name=nm, iptm=d.get("iptm"), iplddt=d.get("complex_iplddt"),
                         plddt=d.get("complex_plddt"), ipde=d.get("complex_ipde")))
    if not rows: log("無結果"); return
    D=pd.DataFrame(rows)
    R=D.groupby("name").agg(n=("iptm","size"), iptm_mean=("iptm","mean"), iptm_sd=("iptm","std"),
                            iptm_max=("iptm","max"), plddt_mean=("plddt","mean"),
                            iplddt_mean=("iplddt","mean"), ipde_mean=("ipde","mean")).reset_index()
    R["kind"]=R.name.map(lambda n: meta.get(n,{}).get("kind","?"))
    R["hotspot"]=R.name.map(lambda n: meta.get(n,{}).get("hotspot"))
    R=R.sort_values("iptm_mean",ascending=False)
    R.to_csv(f"{RES}/88_fsp1_template_scores.csv", index=False)
    pd.set_option("display.width",230)
    log("\n=== FSP1 結合劑（給晶體結構當模板）===")
    log("\n"+R[["name","kind","hotspot","iptm_mean","iptm_sd","iptm_max","plddt_mean","iplddt_mean","ipde_mean"]].round(3).to_string(index=False))
    log(f"\n標的 pLDDT（給模板後）：{R.plddt_mean.mean():.3f}　"
        f"（前一次無模板時只有 0.46–0.47）")
    log("DONE")
main()
