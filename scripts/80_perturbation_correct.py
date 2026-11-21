#!/usr/bin/env python
"""
步驟 80：修正版的蛋白質層級擾動實驗 —— 把突變對到「正確的結合夥伴」上。

前一版（步驟 74）的錯誤：
  我把 DNA 結合鹼性區的突變拿去測 MafG–NRF2 複合體，
  但那兩個表面不重疊，所以測不到任何東西（ΔipTM 只有 ±0.04）。

修正後的配對：
  A 組｜DNA 結合區突變 (K53/R56/R57/T58/K60/N61/R62/Y64/R69/K71/R72)
        → 測 **MafG + MARE DNA**（Boltz-2 支援 DNA 鏈）
  B 組｜二聚體介面突變 (Q75/L79/E80/Q82/K83/L86/Q87/E89/L93/N97/M100)
        → 測 **MafG + NRF2**
  C 組｜磷酸化擬態，同樣各自對到正確夥伴
  D 組｜對照：在「不相關表面」上的突變，驗證訊號不是雜訊

MARE 序列取自共結晶 7X5E：TTGCTGACTCATCAT
"""
import os, sys, json, glob, subprocess, argparse
import numpy as np, pandas as pd

BASE="/raid/jupyter-alz.ee09/Tonsil_to_MCL"; RES=f"{BASE}/results"
MAFG="SLTDEELVTMSVRELNQHLRGLSKEEIVQLKQRRRTLKNRGYAASCRVKRVTQKEELEKQKAELQQEVEKLASENASMKLELDALRSKYEALQTFARTV"
OFF=23
NRF2="HLTRDELRAKALHIPFPVEKIINLPVVDFNEMMSKEQFNEAQLALIRDIRRRRGKNKVAAQNCRKRKLENNIVELEQDLDHLKDEKEKLLKEKGENDKSLHLLKKQLS"
MARE_T="TTGCTGACTCATCAT"; MARE_B="ATGATGAGTCAGCAA"
DNA_RES=[53,56,57,58,60,61,62,64,69,71,72]
DIM_RES=[75,79,80,82,83,86,87,89,93,97,100]
CTRL_RES=[30,45,115]      # 既不在 DNA 面也不在二聚面
def log(*a): print("[80]",*a,flush=True)
def mut(seq,p,new):
    i=p-OFF; return seq[:i]+new+seq[i+1:], seq[i]

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--gpu", default="4"); ap.add_argument("--samples", type=int, default=1)
    a=ap.parse_args()
    out=f"{RES}/80_perturb_correct"; os.makedirs(out, exist_ok=True)
    yml=os.path.join(out,"inputs"); os.makedirs(yml, exist_ok=True)
    for f in glob.glob(os.path.join(yml,"*.yaml")): os.remove(f)
    dms=pd.read_csv(f"{RES}/71_esm_dms_MAFG.csv.gz")
    meta={}
    def write(name, prot, with_dna):
        with open(os.path.join(yml,f"{name}.yaml"),"w") as fh:
            fh.write("version: 1\nsequences:\n")
            fh.write(f"  - protein:\n      id: A\n      sequence: {prot}\n      msa: empty\n")
            if with_dna:
                fh.write(f"  - dna:\n      id: C\n      sequence: {MARE_T}\n")
                fh.write(f"  - dna:\n      id: D\n      sequence: {MARE_B}\n")
            else:
                fh.write(f"  - protein:\n      id: B\n      sequence: {NRF2}\n      msa: empty\n")

    write("WT_MafG_DNA", MAFG, True);  meta["WT_MafG_DNA"]=dict(assay="MafG+DNA", grp="野生型", desc="野生型 MafG + MARE DNA")
    write("WT_MafG_NRF2", MAFG, False);meta["WT_MafG_NRF2"]=dict(assay="MafG+NRF2", grp="野生型", desc="野生型 MafG + NRF2")

    def add(p, new, assay, grp, tag=""):
        s,wt=mut(MAFG,p,new)
        nm=f"{'D' if assay=='MafG+DNA' else 'P'}_{wt}{p}{new}{tag}"
        write(nm, s, assay=="MafG+DNA")
        llr=dms[(dms.pos==p)&(dms.mt==new)].llr
        meta[nm]=dict(assay=assay, grp=grp, desc=f"MafG {wt}{p}{new}", pos=p, wt=wt, mt=new,
                      esm_llr=float(llr.iloc[0]) if len(llr) else None)

    for p in DNA_RES:   add(p,"A","MafG+DNA","A組 DNA面突變→測DNA")
    for p in DIM_RES:   add(p,"A","MafG+NRF2","B組 二聚面突變→測NRF2")
    for p in CTRL_RES:
        add(p,"A","MafG+DNA","D組 無關表面對照"); add(p,"A","MafG+NRF2","D組 無關表面對照")
    # 磷酸化擬態：各自對到正確夥伴
    for p in [58,64,67]:
        for new in ["D","E"]: add(p,new,"MafG+DNA","C組 磷酸化擬態(DNA面)")
    for p in [111]:
        for new in ["D","E"]: add(p,new,"MafG+NRF2","C組 磷酸化擬態(二聚面)")
    # 電荷反轉（比 Ala 更激烈，測上限）
    for p in [53,60,56]: add(p,"D","MafG+DNA","A組 電荷反轉", "_chg")

    log(f"共 {len(meta)} 個複合體")
    env=dict(os.environ, CUDA_VISIBLE_DEVICES=a.gpu)
    cmd=[f"{BASE}/env/boltz/bin/boltz","predict",yml,"--out_dir",out,
         "--cache",f"{BASE}/models/boltz","--output_format","mmcif",
         "--diffusion_samples",str(a.samples),"--recycling_steps","3",
         "--devices","1","--accelerator","gpu","--no_kernels"]
    p=subprocess.run(cmd, env=env, capture_output=True, text=True)
    log("return code:", p.returncode)
    if p.returncode!=0: log(p.stderr[-2000:])
    rows=[]
    for f in glob.glob(os.path.join(out,"**","confidence_*.json"), recursive=True):
        d=json.load(open(f))
        nm=os.path.basename(f)[len("confidence_"):].rsplit("_model_",1)[0]
        rows.append(dict(name=nm, iptm=d.get("iptm"), ptm=d.get("ptm"),
                         iplddt=d.get("complex_iplddt"), plddt=d.get("complex_plddt"),
                         ipde=d.get("complex_ipde")))
    if not rows: log("無結果"); return
    R=pd.DataFrame(rows)
    for k in ["assay","grp","desc","pos","wt","mt","esm_llr"]:
        R[k]=R.name.map(lambda n: meta.get(n,{}).get(k))
    bd=float(R[R.name=="WT_MafG_DNA"].iptm.iloc[0]); bn=float(R[R.name=="WT_MafG_NRF2"].iptm.iloc[0])
    R["base"]=np.where(R.assay=="MafG+DNA", bd, bn)
    R["d_iptm"]=R.iptm-R.base
    R=R.sort_values(["grp","d_iptm"])
    R.to_csv(f"{RES}/80_perturb_correct.csv", index=False)
    pd.set_option("display.width",230)
    log(f"\n野生型基準：MafG+DNA ipTM = {bd:.3f}　｜　MafG+NRF2 ipTM = {bn:.3f}")
    for g,s in R.groupby("grp"):
        log(f"\n=== {g} ===")
        log("\n"+s[["desc","assay","iptm","d_iptm","iplddt","ipde","esm_llr"]].round(3).to_string(index=False))
    log("DONE")
main()
