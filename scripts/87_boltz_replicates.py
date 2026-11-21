#!/usr/bin/env python
"""
步驟 87：最終評估 —— 多次 diffusion 取樣，報告分布而非單一數字。

發現的問題：Boltz-2 單次取樣的 ipTM 變異很大
（同一個設計兩次執行分別得到 0.759 與 0.320）。
所以最終排名必須用多次取樣的平均 ± 標準差，並與對照組一起跑。

同一批內同時含：陽性對照（天然 NRF2）、陰性對照（兩種打散）、全部候選設計。
"""
import os, sys, json, glob, subprocess, argparse, random
import numpy as np, pandas as pd
BASE="/raid/jupyter-alz.ee09/Tonsil_to_MCL"; RES=f"{BASE}/results"
MAFG="SLTDEELVTMSVRELNQHLRGLSKEEIVQLKQRRRTLKNRGYAASCRVKRVTQKEELEKQKAELQQEVEKLASENASMKLELDALRSKYEALQTFARTV"
NRF2="HLTRDELRAKALHIPFPVEKIINLPVVDFNEMMSKEQFNEAQLALIRDIRRRRGKNKVAAQNCRKRKLENNIVELEQDLDHLKDEKEKLLKEKGENDKSLHLLKKQLS"
def log(*a): print("[87]",*a,flush=True)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--gpu", default="4"); ap.add_argument("--samples", type=int, default=5)
    a=ap.parse_args()
    out=f"{RES}/87_boltz_replicates"; os.makedirs(out, exist_ok=True)
    yml=os.path.join(out,"inputs"); os.makedirs(yml, exist_ok=True)
    for f in glob.glob(os.path.join(yml,"*.yaml")): os.remove(f)
    def write(nm, b):
        with open(os.path.join(yml,f"{nm}.yaml"),"w") as fh:
            fh.write("version: 1\nsequences:\n")
            fh.write(f"  - protein:\n      id: A\n      sequence: {MAFG}\n      msa: empty\n")
            fh.write(f"  - protein:\n      id: B\n      sequence: {b}\n      msa: empty\n")
    meta={}
    write("CTRL_POS_NRF2", NRF2); meta["CTRL_POS_NRF2"]=dict(kind="陽性對照", blen=len(NRF2))
    cands=[]
    B1=pd.read_csv(f"{RES}/65_mpnn_mafg_binders.csv").sort_values("score",ascending=False).drop_duplicates("design")
    for _,r in B1.head(4).iterrows():
        cands.append((f"R1_{r.design}_m{int(r.mpnn_idx)}", r.binder,
                      "第一輪·"+("阻斷二聚" if "dimer" in r.design else "阻斷DNA結合"),
                      int(r.hotspot_covered)))
    Z=pd.read_csv(f"{RES}/85_binders_mafgZE.csv").sort_values("score2",ascending=False).drop_duplicates("design")
    for _,r in Z.head(2).iterrows():
        cands.append((f"R2_{r.design}_m{int(r.mpnn_idx)}", r.binder, "第二輪·零逃逸熱點", int(r.hotspot_covered)))
    for nm,b,k,h in cands:
        write(nm,b); meta[nm]=dict(kind=k, blen=len(b), hotspot=h, seq=b)
    rng=random.Random(1)
    for i,(nm,b,k,h) in enumerate(cands[:2]):
        s=list(b); rng.shuffle(s); s="".join(s)
        write(f"CTRL_NEG_{i+1}", s); meta[f"CTRL_NEG_{i+1}"]=dict(kind="陰性對照(打散)", blen=len(s))
    log(f"{len(meta)} 個複合體 × {a.samples} 次取樣")
    env=dict(os.environ, CUDA_VISIBLE_DEVICES=a.gpu)
    cmd=[f"{BASE}/env/boltz/bin/boltz","predict",yml,"--out_dir",out,"--cache",f"{BASE}/models/boltz",
         "--output_format","mmcif","--diffusion_samples",str(a.samples),"--recycling_steps","3",
         "--devices","1","--accelerator","gpu","--no_kernels"]
    p=subprocess.run(cmd, env=env, capture_output=True, text=True)
    log("return code:", p.returncode)
    if p.returncode!=0: log(p.stderr[-2000:])
    rows=[]
    for f in glob.glob(os.path.join(out,"**","confidence_*.json"), recursive=True):
        d=json.load(open(f))
        base=os.path.basename(f)[len("confidence_"):]
        nm=base.rsplit("_model_",1)[0]
        rows.append(dict(name=nm, iptm=d.get("iptm"), iplddt=d.get("complex_iplddt"),
                         ipde=d.get("complex_ipde"), plddt=d.get("complex_plddt")))
    if not rows: log("無結果"); return
    D=pd.DataFrame(rows)
    R=D.groupby("name").agg(n=("iptm","size"),
        iptm_mean=("iptm","mean"), iptm_sd=("iptm","std"), iptm_max=("iptm","max"), iptm_min=("iptm","min"),
        ipde_mean=("ipde","mean"), ipde_sd=("ipde","std"),
        iplddt_mean=("iplddt","mean")).reset_index()
    R["kind"]=R.name.map(lambda n: meta.get(n,{}).get("kind","?"))
    R["blen"]=R.name.map(lambda n: meta.get(n,{}).get("blen"))
    R["hotspot"]=R.name.map(lambda n: meta.get(n,{}).get("hotspot"))
    R=R.sort_values("iptm_mean", ascending=False)
    R.to_csv(f"{RES}/87_boltz_replicates.csv", index=False)
    pd.set_option("display.width",250)
    log(f"\n=== 最終評估：每個複合體 {a.samples} 次 diffusion 取樣 ===")
    log("\n"+R[["name","kind","blen","hotspot","iptm_mean","iptm_sd","iptm_min","iptm_max",
                "ipde_mean","ipde_sd","iplddt_mean"]].round(3).to_string(index=False))
    pos=R[R.kind=="陽性對照"]; neg=R[R.kind.str.contains("陰性")]
    des=R[~R.kind.str.contains("對照")]
    if len(pos) and len(neg):
        log(f"\n陽性對照（天然 NRF2）ipTM = {pos.iptm_mean.iloc[0]:.3f} ± {pos.iptm_sd.iloc[0]:.3f}")
        log(f"陰性對照（打散）ipTM = {neg.iptm_mean.mean():.3f} ± {neg.iptm_sd.mean():.3f}")
        log(f"設計品 ipTM 範圍 {des.iptm_mean.min():.3f} – {des.iptm_mean.max():.3f}")
        from scipy.stats import mannwhitneyu
        if len(des)>2 and len(neg)>=2:
            dd=D[D.name.isin(des.name)]; nn=D[D.name.isin(neg.name)]
            u,pv=mannwhitneyu(dd.iptm, nn.iptm, alternative="greater")
            log(f"設計品 vs 陰性對照（所有取樣點）：Mann-Whitney p = {pv:.4g}"
                f"（n={len(dd)} vs {len(nn)}）")
    json.dump(meta, open(f"{RES}/87_meta.json","w"), ensure_ascii=False, indent=1)
    log("DONE")
main()
