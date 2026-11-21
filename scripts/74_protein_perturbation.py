#!/usr/bin/env python
"""
步驟 74：蛋白質層級的擾動實驗組（四類操作）。

對應使用者要求的四種操作，全部在「折疊 / 結構 / 交互作用」層面評估：

  (1) 點突變（等同 SNP）
      對 MafG 的 DNA 結合區與二聚體介面上，ESM 認定最受約束的殘基，
      做丙胺酸取代與「ESM 最偏好的替代胺基酸」取代，
      預測 MafG–NRF2 複合體的 ipTM 變化。
      → 回答「若這個位置被突變，交互作用會不會垮」

  (2) 增加序列：降解標籤融合（biodegrader）
      把設計好的結合劑接上 E3 連接酶招募序列，
      預測「結合劑融合體 + MafG」是否還結合，以及是否能同時抓住 E3。
      → 回答「加一段序列能不能把阻斷劑變成降解劑」

  (3) 增加序列：核定位訊號與連接子
      binder 前後各加 SV40 NLS 與不同長度 GS linker，
      看預測結構是否仍摺疊、界面是否保留。

  (4) 增加訊號：磷酸化擬態
      介面上的 Ser/Thr/Tyr 換成 Asp/Glu（磷酸化擬態），
      看交互作用是否被破壞 → 這等於問「細胞用磷酸化能不能自己關掉這個介面」

所有評估都用 Boltz-2 共摺疊（single-sequence，不外送序列）。
"""
import os, sys, json, glob, subprocess, argparse
import numpy as np, pandas as pd

BASE="/raid/jupyter-alz.ee09/Tonsil_to_MCL"; RES=f"{BASE}/results"
# 7X5E 鏈 A = MafG 殘基 23–121；UniProt O15525 編號與之一致
MAFG_23_121="SLTDEELVTMSVRELNQHLRGLSKEEIVQLKQRRRTLKNRGYAASCRVKRVTQKEELEKQKAELQQEVEKLASENASMKLELDALRSKYEALQTFARTV"
MAFG_OFFSET=23
NRF2="HLTRDELRAKALHIPFPVEKIINLPVVDFNEMMSKEQFNEAQLALIRDIRRRRGKNKVAAQNCRKRKLENNIVELEQDLDHLKDEKEKLLKEKGENDKSLHLLKKQLS"
# E3 招募序列
DEGRONS={
  # VHL 的 HIF1a ODD 核心（需脯胺酸羥化，實務上要配合羥化酶；此處作為序列層面的可行性測試）
  "VHL_ODD":"DALTLLAPAAGDTIISLDF",
  # SPOP 的 SBC 共識序列（Ser/Thr-rich，不需 PTM）
  "SPOP_SBC":"ADSSTSS",
  # KEAP1 的 NRF2 ETGE 高親和 degron —— 在 NRF2 高活性的細胞裡 KEAP1 仍在，可借用
  "KEAP1_ETGE":"LQLDEETGEFLPIQ",
}
NLS="PKKKRKV"
def log(*a): print("[74]",*a,flush=True)

def mutate(seq, offset, pos, new):
    i=pos-offset
    assert 0<=i<len(seq), (pos,i,len(seq))
    return seq[:i]+new+seq[i+1:], seq[i]

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--gpu", default="4"); ap.add_argument("--n_res", type=int, default=8)
    ap.add_argument("--samples", type=int, default=1)
    ap.add_argument("--experiments", nargs="+",
                    default=["point","degron","nls","phospho"])
    a=ap.parse_args()
    out=f"{RES}/74_protein_perturbation"; os.makedirs(out, exist_ok=True)
    yml=os.path.join(out,"inputs"); os.makedirs(yml, exist_ok=True)
    for f in glob.glob(os.path.join(yml,"*.yaml")): os.remove(f)
    meta={}

    def write(name, chains):
        """chains = [(id, seq), ...]"""
        with open(os.path.join(yml,f"{name}.yaml"),"w") as fh:
            fh.write("version: 1\nsequences:\n")
            for cid,s in chains:
                fh.write(f"  - protein:\n      id: {cid}\n      sequence: {s}\n      msa: empty\n")

    # 讀 ESM 約束與介面
    dr=pd.read_csv(f"{RES}/73_druggable_residues.csv") if os.path.exists(f"{RES}/73_druggable_residues.csv") else None
    cons=pd.read_csv(f"{RES}/71_esm_position_constraint.csv")
    mg=cons[cons.gene=="MAFG"].set_index("pos")
    dms=pd.read_csv(f"{RES}/71_esm_dms_MAFG.csv.gz")

    # --- 基準 ---
    write("BASE_MafG_NRF2", [("A",MAFG_23_121),("B",NRF2)])
    meta["BASE_MafG_NRF2"]=dict(exp="基準", desc="野生型 MafG–NRF2")

    # --- (1) 點突變 ---
    if "point" in a.experiments:
        if dr is not None:
            cand=dr[(dr.gene=="MAFG")&(dr.constrained)].sort_values("pct")
        else:
            cand=pd.DataFrame()
        if not len(cand):
            # 退而求其次：取 ESM 最受約束且落在 23–121 的位置
            c=mg[(mg.index>=MAFG_OFFSET)&(mg.index<=121)].sort_values("mean_llr")
            cand=pd.DataFrame(dict(pos=c.index[:a.n_res], aa=c.wt.values[:a.n_res], iface="esm_top"))
        cand=cand.head(a.n_res)
        log(f"(1) 點突變：選了 {len(cand)} 個受約束的介面殘基")
        for _,r in cand.iterrows():
            p=int(r.pos)
            if not (MAFG_OFFSET<=p<=121): continue
            for new,tag in [("A","Ala"),]:
                s,wt=mutate(MAFG_23_121, MAFG_OFFSET, p, new)
                nm=f"PT_{wt}{p}{new}"
                write(nm, [("A",s),("B",NRF2)])
                meta[nm]=dict(exp="點突變", desc=f"MafG {wt}{p}{new}（{r.get('iface','')}）",
                              pos=p, wt=wt, mt=new,
                              esm_llr=float(dms[(dms.pos==p)&(dms.mt==new)].llr.iloc[0])
                                        if len(dms[(dms.pos==p)&(dms.mt==new)]) else None)
            # ESM 最偏好的替代
            sub=dms[dms.pos==p].sort_values("llr", ascending=False)
            if len(sub):
                best=sub.iloc[0]
                s,wt=mutate(MAFG_23_121, MAFG_OFFSET, p, best.mt)
                nm=f"PT_{wt}{p}{best.mt}_esmbest"
                write(nm, [("A",s),("B",NRF2)])
                meta[nm]=dict(exp="點突變", desc=f"MafG {wt}{p}{best.mt}（ESM 最偏好的替代）",
                              pos=p, wt=wt, mt=best.mt, esm_llr=float(best.llr))

    # --- (4) 磷酸化擬態 ---
    if "phospho" in a.experiments:
        sty=[(i+MAFG_OFFSET, c) for i,c in enumerate(MAFG_23_121) if c in "STY"]
        if dr is not None:
            ifpos=set(dr[dr.gene=="MAFG"].pos.astype(int))
            sty=[(p,c) for p,c in sty if p in ifpos]
        log(f"(4) 磷酸化擬態：介面上的 Ser/Thr/Tyr 共 {len(sty)} 個")
        for p,c in sty[:6]:
            for new in ["D","E"]:
                s,wt=mutate(MAFG_23_121, MAFG_OFFSET, p, new)
                nm=f"PH_{wt}{p}{new}"
                write(nm, [("A",s),("B",NRF2)])
                meta[nm]=dict(exp="磷酸化擬態", desc=f"MafG {wt}{p}{new}（模擬磷酸化）", pos=p, wt=wt, mt=new)

    # --- (2)(3) 增加序列：以最佳結合劑為骨幹 ---
    if ("degron" in a.experiments) or ("nls" in a.experiments):
        B=pd.read_csv(f"{RES}/65_mpnn_mafg_binders.csv").sort_values("score",ascending=False)
        best=B.iloc[0]
        binder=best.binder
        log(f"(2)(3) 以 {best.design}(m{int(best.mpnn_idx)}) 為骨幹，結合劑 {len(binder)} aa")
        write("ADD_binder_WT", [("A",MAFG_23_121),("B",binder)])
        meta["ADD_binder_WT"]=dict(exp="增加序列", desc="原始結合劑（基準）", blen=len(binder))
        if "degron" in a.experiments:
            for dn,ds in DEGRONS.items():
                for where in ["C","N"]:
                    fused = binder+"GGSGGSGGS"+ds if where=="C" else ds+"GGSGGSGGS"+binder
                    nm=f"ADD_degron_{dn}_{where}term"
                    write(nm, [("A",MAFG_23_121),("B",fused)])
                    meta[nm]=dict(exp="增加序列-降解標籤",
                                  desc=f"結合劑 + {dn}（接在 {where} 端）", blen=len(fused))
        if "nls" in a.experiments:
            for L in [0,4,12]:
                link="GGSG"*(L//4)
                fused=NLS+link+binder
                nm=f"ADD_NLS_link{L}"
                write(nm, [("A",MAFG_23_121),("B",fused)])
                meta[nm]=dict(exp="增加序列-核定位訊號",
                              desc=f"SV40 NLS + {L} aa linker + 結合劑", blen=len(fused))

    log(f"\n共 {len(meta)} 個複合體要預測")
    env=dict(os.environ, CUDA_VISIBLE_DEVICES=a.gpu)
    cmd=[f"{BASE}/env/boltz/bin/boltz","predict",yml,"--out_dir",out,
         "--cache",f"{BASE}/models/boltz","--output_format","mmcif",
         "--diffusion_samples",str(a.samples),"--recycling_steps","3",
         "--devices","1","--accelerator","gpu","--no_kernels"]
    p=subprocess.run(cmd, env=env, capture_output=True, text=True)
    log("Boltz return code:", p.returncode)
    if p.returncode!=0: log(p.stderr[-2500:])

    rows=[]
    for f in glob.glob(os.path.join(out,"**","confidence_*.json"), recursive=True):
        d=json.load(open(f))
        nm=os.path.basename(f)[len("confidence_"):].rsplit("_model_",1)[0]
        rows.append(dict(name=nm, iptm=d.get("iptm"), ptm=d.get("ptm"),
                         iplddt=d.get("complex_iplddt"), plddt=d.get("complex_plddt"),
                         ipde=d.get("complex_ipde"), conf=d.get("confidence_score")))
    if not rows: log("無結果"); return
    R=pd.DataFrame(rows).groupby("name").agg(
        iptm=("iptm","max"), ptm=("ptm","max"), iplddt=("iplddt","max"),
        plddt=("plddt","max"), ipde=("ipde","min")).reset_index()
    for k in ["exp","desc","esm_llr","blen","pos","wt","mt"]:
        R[k]=R.name.map(lambda n: meta.get(n,{}).get(k))
    base=R[R.name=="BASE_MafG_NRF2"].iptm
    b=float(base.iloc[0]) if len(base) else np.nan
    bb=R[R.name=="ADD_binder_WT"].iptm
    bb=float(bb.iloc[0]) if len(bb) else np.nan
    R["d_iptm"]=np.where(R.exp.astype(str).str.startswith("增加序列"), R.iptm-bb, R.iptm-b)
    R=R.sort_values(["exp","d_iptm"])
    R.to_csv(f"{RES}/74_protein_perturbation.csv", index=False)
    pd.set_option("display.width",250); pd.set_option("display.max_colwidth",46)
    log(f"\n基準：野生型 MafG–NRF2 ipTM = {b:.3f}；原始結合劑 ipTM = {bb:.3f}")
    for e,sub in R.groupby("exp"):
        log(f"\n=== {e} ===")
        log("\n"+sub[["desc","iptm","d_iptm","iplddt","ipde","esm_llr"]].round(3).to_string(index=False))
    json.dump(meta, open(f"{RES}/74_perturbation_meta.json","w"), ensure_ascii=False, indent=1)
    log("DONE")
main()
