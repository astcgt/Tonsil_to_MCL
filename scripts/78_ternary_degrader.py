#!/usr/bin/env python
"""
步驟 78：三元複合體預測 —— 把「阻斷劑」變成「降解劑」的可行性測試。

概念
----
阻斷劑只是擋住 MafG；降解劑則是把 MafG 拖去泛素化降解，效果更徹底而且用量更低。
做法是在設計好的結合劑上「增加一段序列」（E3 連接酶招募肽），
形成：  MafG  ←→  結合劑–linker–degron  ←→  E3 連接酶
這就是使用者說的「增加序列」在蛋白質交互作用層面的具體操作。

評估
----
Boltz-2 會回報 `pair_chains_iptm`，可以分別看：
  ipTM(A,B) —— 融合體還抓不抓得住 MafG（不能因為加了序列就失去結合力）
  ipTM(B,C) —— 融合體抓不抓得住 E3
  兩者都要夠高，三元複合體才成立。

E3 選擇
  VHL（VHL–ElonginC–ElonginB 複合體的 VHL 部分）：degron 用 HIF1α ODD
  KEAP1 Kelch domain：degron 用 NRF2 的 ETGE（高親和、不需轉譯後修飾）
  —— 在 NRF2 高活性的 MCL 細胞裡 KEAP1 本來就存在且活躍，是現成的工具
"""
import os, sys, json, glob, subprocess, argparse
import numpy as np, pandas as pd

BASE="/raid/jupyter-alz.ee09/Tonsil_to_MCL"; RES=f"{BASE}/results"
MAFG="SLTDEELVTMSVRELNQHLRGLSKEEIVQLKQRRRTLKNRGYAASCRVKRVTQKEELEKQKAELQQEVEKLASENASMKLELDALRSKYEALQTFARTV"
# KEAP1 Kelch domain (UniProt Q14145 約 321–609)
KEAP1_KELCH=("GVGVAVTMEPRVHEWQEVLNQRRKQRLEAHTSHSGSHLHSNTSHSHSHSHSSHSHSSHSHSHSHSS")  # 佔位，見下方下載
VHL=("MPRRAENWDEAEVGAEEAGVEEYGPEEDGGEESGAEESGPEESGPEELGAEEEMEAGRPRPVLRSVNSREPSQVIFCNRSPRVVLPVWLNFDGEPQPYPTLPPGTGRRIHSYRGHLWLFRDAGTHDGLLVNQTELFVPSLNVDGQPIFANITLPVYTLKERCLQVVRSLVKPENYRRLDIVRSLYEDLEDHPNVQKDLERLTQERIAHQRMGD")
DEGRONS={"KEAP1_ETGE":"LQLDEETGEFLPIQ", "VHL_ODD":"DALTLLAPAAGDTIISLDF"}
def log(*a): print("[78]",*a,flush=True)

def fetch_keap1_kelch():
    import urllib.request, json as js
    with urllib.request.urlopen("https://rest.uniprot.org/uniprotkb/Q14145.json?fields=sequence",timeout=30) as r:
        s=js.load(r)["sequence"]["value"]
    return s[320:609]      # Kelch repeats

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--gpu", default="4"); ap.add_argument("--samples", type=int, default=1)
    a=ap.parse_args()
    out=f"{RES}/78_ternary"; os.makedirs(out, exist_ok=True)
    yml=os.path.join(out,"inputs"); os.makedirs(yml, exist_ok=True)
    for f in glob.glob(os.path.join(yml,"*.yaml")): os.remove(f)
    kelch=fetch_keap1_kelch(); log(f"KEAP1 Kelch domain {len(kelch)} aa")

    B=pd.read_csv(f"{RES}/65_mpnn_mafg_binders.csv").sort_values("score",ascending=False)
    cal=pd.read_csv(f"{RES}/72_boltz_calibrated.csv")
    # 用 Boltz 實測 ipTM 最高的那個設計當骨幹
    best_name=cal[~cal.kind.str.contains("對照")].iloc[0]["name"]
    log("骨幹：", best_name)
    key=best_name.rsplit("_m",1)[0]
    cand=B[B.design==key]
    if not len(cand): cand=B
    binder=cand.iloc[0].binder
    log(f"結合劑 {len(binder)} aa")

    meta={}
    def write(name, chains):
        with open(os.path.join(yml,f"{name}.yaml"),"w") as fh:
            fh.write("version: 1\nsequences:\n")
            for cid,s in chains:
                fh.write(f"  - protein:\n      id: {cid}\n      sequence: {s}\n      msa: empty\n")

    # 二元基準
    write("BIN_binder_MafG", [("A",MAFG),("B",binder)])
    meta["BIN_binder_MafG"]=dict(kind="二元基準", desc="結合劑 + MafG")
    write("BIN_ETGE_KEAP1", [("B",DEGRONS["KEAP1_ETGE"]),("C",kelch)])
    meta["BIN_ETGE_KEAP1"]=dict(kind="二元基準", desc="ETGE 肽 + KEAP1 Kelch")
    # 三元
    for dn,ds in DEGRONS.items():
        e3 = kelch if dn.startswith("KEAP1") else VHL
        for where in ["C","N"]:
            for link in ["GGSGGSGGS","GGSGGSGGSGGSGGSGGS"]:
                fused = binder+link+ds if where=="C" else ds+link+binder
                nm=f"TER_{dn}_{where}term_L{len(link)}"
                write(nm, [("A",MAFG),("B",fused),("C",e3)])
                meta[nm]=dict(kind="三元降解劑", degron=dn, term=where, linker=len(link),
                              desc=f"MafG + 結合劑-{dn}({where}端, linker {len(link)}aa) + "
                                   f"{'KEAP1 Kelch' if dn.startswith('KEAP1') else 'VHL'}",
                              blen=len(fused))
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
        pc=d.get("pair_chains_iptm",{})
        def g(i,j):
            try: return float(pc[str(i)][str(j)])
            except Exception: return np.nan
        rows.append(dict(name=nm, iptm=d.get("iptm"), ptm=d.get("ptm"),
                         plddt=d.get("complex_plddt"), iplddt=d.get("complex_iplddt"),
                         ipde=d.get("complex_ipde"),
                         iptm_AB=g(0,1), iptm_AC=g(0,2), iptm_BC=g(1,2)))
    if not rows: log("無結果"); return
    R=pd.DataFrame(rows)
    for k in ["kind","desc","degron","term","linker"]:
        R[k]=R.name.map(lambda n: meta.get(n,{}).get(k))
    R=R.sort_values(["kind","iptm"], ascending=[True,False])
    R.to_csv(f"{RES}/78_ternary.csv", index=False)
    pd.set_option("display.width",250); pd.set_option("display.max_colwidth",56)
    log("\n=== 三元降解劑複合體預測 ===")
    log("（ipTM_AB = 融合體抓 MafG；ipTM_BC = 融合體抓 E3；兩者都要高）")
    log("\n"+R[["desc","iptm","iptm_AB","iptm_BC","iptm_AC","iplddt","ipde"]].round(3).to_string(index=False))
    json.dump(meta, open(f"{RES}/78_ternary_meta.json","w"), ensure_ascii=False, indent=1)
    log("DONE")
main()
