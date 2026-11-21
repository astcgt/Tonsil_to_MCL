#!/usr/bin/env python
"""
步驟 72：Boltz-2 共摺疊評估（含校準對照）。

單看一個 ipTM 數字沒有意義，必須有參照點。所以同時跑：
  * **陽性對照**：MafG–NRF2 天然異二聚體（有 2.3 Å 共結晶，PDB 7X5E）
                  這是「真的會結合」的上限參考
  * **陰性對照**：把設計序列隨機打散（保留胺基酸組成）
                  這是「不會結合」的下限參考
  * 待測：RFdiffusion3 + ProteinMPNN 設計出來的結合劑

評估指標
  ipTM           界面預測模板建模分數（binder 設計最常用的篩選指標）
  complex_iplddt 界面每殘基信心
  complex_ipde   界面距離誤差（越小越好）

隱私：single-sequence 模式（空 MSA），設計序列不離開本機。
"""
import os, sys, json, glob, subprocess, argparse, random
import numpy as np, pandas as pd

BASE="/raid/jupyter-alz.ee09/Tonsil_to_MCL"; RES=f"{BASE}/results"
MAFG="SLTDEELVTMSVRELNQHLRGLSKEEIVQLKQRRRTLKNRGYAASCRVKRVTQKEELEKQKAELQQEVEKLASENASMKLELDALRSKYEALQTFARTV"
# NRF2 (NFE2L2) 的 Neh1 CNC-bZIP 區，取自 7X5E 鏈 B
NRF2="HLTRDELRAKALHIPFPVEKIINLPVVDFNEMMSKEQFNEAQLALIRDIRRRRGKNKVAAQNCRKRKLENNIVELEQDLDHLKDEKEKLLKEKGENDKSLHLLKKQLS"
# BACH1 的 bZIP 區（MafG 的另一個夥伴，UniProt O14867 約 566-650）
BACH1="EQLSPGSRAHSEAHSPSSSVHSPIGSNSNSGDSSHSSSHSSQHSQHSAHSQHSHSHSHSQ"
def log(*a): print("[72]",*a,flush=True)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--gpu", default="4"); ap.add_argument("--top", type=int, default=6)
    ap.add_argument("--samples", type=int, default=2)
    a=ap.parse_args()
    out=f"{RES}/72_boltz_calibrated"; os.makedirs(out, exist_ok=True)
    yml=os.path.join(out,"inputs"); os.makedirs(yml, exist_ok=True)
    for f in glob.glob(os.path.join(yml,"*.yaml")): os.remove(f)

    def write(name, chainB, chainA=MAFG):
        with open(os.path.join(yml,f"{name}.yaml"),"w") as fh:
            fh.write("version: 1\nsequences:\n")
            fh.write(f"  - protein:\n      id: A\n      sequence: {chainA}\n      msa: empty\n")
            fh.write(f"  - protein:\n      id: B\n      sequence: {chainB}\n      msa: empty\n")

    meta={}
    # 陽性對照
    write("CTRL_POS_MafG_NRF2", NRF2); meta["CTRL_POS_MafG_NRF2"]=dict(kind="陽性對照(天然 NRF2)", blen=len(NRF2))
    # 設計
    R=pd.read_csv(f"{RES}/65_mpnn_mafg_binders.csv")
    top=R.sort_values("score",ascending=False).drop_duplicates("design").head(a.top)
    rng=random.Random(0)
    for _,r in top.iterrows():
        nm=f"{r.design}_m{int(r.mpnn_idx)}"
        write(nm, r.binder); meta[nm]=dict(kind=("阻斷二聚" if "dimer" in r.design else "阻斷 DNA 結合"),
                                           blen=int(r.blen), hotspot=int(r.hotspot_covered),
                                           iface=int(r.iface_pairs), seq=r.binder)
    # 陰性對照：打散前兩名設計的序列
    for i,(_,r) in enumerate(top.head(2).iterrows()):
        s=list(r.binder); rng.shuffle(s); s="".join(s)
        nm=f"CTRL_NEG_scram{i+1}"
        write(nm, s); meta[nm]=dict(kind="陰性對照(打散)", blen=len(s), seq=s)
    log(f"共 {len(meta)} 個複合體要預測")

    env=dict(os.environ, CUDA_VISIBLE_DEVICES=a.gpu)
    cmd=[f"{BASE}/env/boltz/bin/boltz","predict",yml,"--out_dir",out,
         "--cache",f"{BASE}/models/boltz","--output_format","mmcif",
         "--diffusion_samples",str(a.samples),"--recycling_steps","3",
         "--devices","1","--accelerator","gpu","--no_kernels"]
    log("執行 Boltz-2 ...")
    p=subprocess.run(cmd, env=env, capture_output=True, text=True)
    log("return code:", p.returncode)
    if p.returncode!=0: log(p.stderr[-2000:])

    rows=[]
    for f in glob.glob(os.path.join(out,"**","confidence_*.json"), recursive=True):
        d=json.load(open(f))
        base=os.path.basename(f)[len("confidence_"):].rsplit("_model_",1)[0]
        rows.append(dict(name=base, iptm=d.get("iptm"), ptm=d.get("ptm"),
                         conf=d.get("confidence_score"),
                         iplddt=d.get("complex_iplddt"), plddt=d.get("complex_plddt"),
                         ipde=d.get("complex_ipde")))
    if not rows:
        log("沒有結果"); return
    R2=pd.DataFrame(rows).groupby("name").agg(
        iptm=("iptm","max"), ptm=("ptm","max"), conf=("conf","max"),
        iplddt=("iplddt","max"), plddt=("plddt","max"), ipde=("ipde","min"),
        n_samples=("iptm","size")).reset_index()
    R2["kind"]=R2.name.map(lambda n: meta.get(n,{}).get("kind","?"))
    R2["binder_len"]=R2.name.map(lambda n: meta.get(n,{}).get("blen"))
    R2["hotspot"]=R2.name.map(lambda n: meta.get(n,{}).get("hotspot"))
    R2=R2.sort_values("iptm", ascending=False)
    R2.to_csv(f"{RES}/72_boltz_calibrated.csv", index=False)
    pd.set_option("display.width",230)
    log("\n=== Boltz-2 共摺疊評估（含校準對照）===")
    log("\n"+R2[["name","kind","binder_len","hotspot","iptm","iplddt","ipde","ptm","plddt"]].round(3).to_string(index=False))
    pos=R2[R2.kind.str.contains("陽性")].iptm.max()
    neg=R2[R2.kind.str.contains("陰性")].iptm.max()
    des=R2[~R2.kind.str.contains("對照")]
    log(f"\n校準：天然 MafG–NRF2 ipTM = {pos:.3f}；打散對照最高 = {neg:.3f}")
    log(f"設計品 ipTM 範圍 {des.iptm.min():.3f} – {des.iptm.max():.3f}，"
        f"超過陰性對照的有 {int((des.iptm>neg).sum())}/{len(des)} 個，"
        f"達到或超過天然複合體的有 {int((des.iptm>=pos).sum())}/{len(des)} 個")
    json.dump(meta, open(f"{RES}/72_boltz_meta.json","w"), ensure_ascii=False, indent=1)
    log("DONE")
main()
