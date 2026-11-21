#!/usr/bin/env python
"""
步驟 64：用 Boltz-2 共摺疊評估設計出來的結合劑。

讀出什麼
--------
Boltz-2 對「標的 + 結合劑」兩條鏈一起摺疊，回報：
  * ipTM  —— 界面預測模板建模分數，是 de novo binder 最常用的篩選指標
  * pTM   —— 整體摺疊信心
  * complex_plddt —— 每殘基信心的平均
這些共同反映「模型認為這個複合體長得出來、而且界面是對的」的程度。

隱私考量
--------
刻意用 **single-sequence 模式（空 MSA）**，不呼叫外部 MSA 伺服器，
所以設計出來的序列不會離開這台機器。
代價是標的摺疊精度略低；但對「同一標的下比較不同結合劑」這個用途，
相對排序仍然有效。
"""
import os, sys, json, glob, argparse, subprocess, shutil
import numpy as np, pandas as pd

BASE="/raid/jupyter-alz.ee09/Tonsil_to_MCL"; RES=f"{BASE}/results"
def log(*a): print("[64]",*a,flush=True)

TARGETS = {
  "FSP1": "MGSQVSVESGALHVVIVGGGFGGIAAASQLQALNVPFMLVDMKDSFHHNVAALRASVETGFAKKTFISYSVTFKDNFRQGLVVGIDLKNQMVLLQGGEALPFSHLILATGSTGPFPGKFNEVSSQQAAIQAYEDMVRQVQRSRFIVVVGGGSAGVEMAAEIKTEYPEKEVTLIHSQVALADKELLPSVRQEVKEILLRKGVQLLLSERVSNLEELPLNEYREYIKVQTDKGTEVATNLVILCTGIKINSSAYRKAFESRLASSGALRVNEHLQVEGHSNVYAIGDCADVRTPKMAYLAGLHANIAVANIVNSVKQRPLQAYKPGALTFLLSMGRNDGVGQISGFYVGRLMVRLTKSRDLFVSKSWKTMRQSPP",
  # 標的序列取自共結晶結構的該鏈
  "MAFG": "SLTDEELVTMSVRELNQHLRGLSKEEIVQLKQRRRTLKNRGYAASCRVKRVTQKEELEKQKAELQQEVEKLASENASMKLELDALRSKYEALQTFARTV",
  "IRF4": "GNGKLRQWLIDQIDSGKYPGLVWENEEKSIFRIPWKHAGKQDYNREEDAALFKAWALFKGKFREGIDKPDPPTWKTRLRCALNKSNDFEELVERSQLDISDPYKVYRIVPEGAK",
}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--ranking", required=True)
    ap.add_argument("--target", required=True)
    ap.add_argument("--top", type=int, default=8)
    ap.add_argument("--gpu", default="4")
    ap.add_argument("--out", default=None)
    a=ap.parse_args()
    out = a.out or f"{RES}/64_boltz_{a.target.lower()}"
    os.makedirs(out, exist_ok=True)
    env=dict(os.environ, CUDA_VISIBLE_DEVICES=a.gpu)

    df=pd.read_csv(a.ranking).head(a.top)
    log(f"評估 {a.target} 的前 {len(df)} 個設計")
    yml_dir=os.path.join(out,"inputs"); os.makedirs(yml_dir, exist_ok=True)
    for _,r in df.iterrows():
        seq=str(r.sequence).replace("X","G")
        nm=str(r.design).replace("/","_")
        with open(os.path.join(yml_dir,f"{nm}.yaml"),"w") as fh:
            fh.write("version: 1\nsequences:\n")
            fh.write(f"  - protein:\n      id: A\n      sequence: {TARGETS[a.target]}\n      msa: empty\n")
            fh.write(f"  - protein:\n      id: B\n      sequence: {seq}\n      msa: empty\n")
    log(f"寫出 {len(df)} 個 YAML 到 {yml_dir}")

    cmd=[f"{BASE}/env/boltz/bin/boltz","predict",yml_dir,
         "--out_dir",out,"--cache",f"{BASE}/models/boltz",
         "--output_format","mmcif","--diffusion_samples","1",
         "--recycling_steps","3","--devices","1","--accelerator","gpu","--no_kernels"]
    log("執行:", " ".join(cmd))
    p=subprocess.run(cmd, env=env, capture_output=True, text=True)
    log("return code:", p.returncode)
    if p.returncode!=0:
        log("stderr 末段:\n"+p.stderr[-2500:])
    rows=[]
    for f in sorted(glob.glob(os.path.join(out,"**","confidence_*.json"), recursive=True)):
        d=json.load(open(f))
        nm=os.path.basename(f).replace("confidence_","").replace("_model_0.json","")
        rows.append(dict(design=nm,
                         confidence=d.get("confidence_score"),
                         ptm=d.get("ptm"), iptm=d.get("iptm"),
                         complex_plddt=d.get("complex_plddt"),
                         complex_iplddt=d.get("complex_iplddt"),
                         complex_pde=d.get("complex_pde"),
                         complex_ipde=d.get("complex_ipde")))
    if rows:
        R=pd.DataFrame(rows).sort_values("iptm", ascending=False)
        R=R.merge(df[["design","binder_len","hotspot_covered","iface_pairs","sequence"]],
                  on="design", how="left")
        R.to_csv(f"{RES}/64_boltz_{a.target.lower()}_scores.csv", index=False)
        pd.set_option("display.width",230)
        log(f"\n=== Boltz-2 共摺疊評估（{a.target}）===")
        log("\n"+R[["design","iptm","ptm","complex_plddt","binder_len","hotspot_covered","iface_pairs"]].round(3).to_string(index=False))
        log(f"\n寫出 {RES}/64_boltz_{a.target.lower()}_scores.csv")
    else:
        log("沒有找到 confidence JSON —— 檢查上面的 stderr")
main()
