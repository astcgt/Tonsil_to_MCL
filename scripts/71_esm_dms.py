#!/usr/bin/env python
"""
步驟 71：用蛋白質語言模型對軸線標的做「完整 in-silico 飽和突變掃描」。

為什麼要做
----------
22 位病人的全基因組定序顯示，這條抗氧化軸線上的 21 個基因**一個編碼突變都沒有**
（同一份資料在已知 driver 基因上抓到 50 個突變、參考鹼基 100% 相符，所以不是流程問題）。
沒有天然突變可用，就改用模型把「每個位置換成每一種胺基酸」全部算一遍，
得到這個蛋白質的**適應度地景 (fitness landscape)**：
  * 分數越低 = 該位置越不能動 = 功能核心 = 藥物該打的地方
  * 分數高的位置 = 可以容忍突變 = 未來抗藥性可能從這裡來

方法：ESM-2 650M 的 masked-marginals
  對位置 i 蓋上 <mask> → 前向傳遞 → 取該位置的 log-softmax
  LLR = log p(突變型) − log p(野生型)
（Meier et al. 2021 的標準做法；650M 在吞吐與精度上是最佳折衷，
 瓶頸是「位點數 × 序列長度」而非模型大小。）
"""
import os, sys, json, argparse
import numpy as np, pandas as pd

BASE="/raid/jupyter-alz.ee09/Tonsil_to_MCL"; RES=f"{BASE}/results"
AA="ACDEFGHIKLMNPQRSTVWY"
def log(*a): print("[71]",*a,flush=True)

# 標的序列（取自 UniProt 正規序列 / 共結晶）
TARGETS={}

def load_targets():
    import urllib.request, json as js
    spec={"MAFG":"O15525","AIFM2":"Q9BRQ8","SLC7A11":"Q9UPY5","SLC3A2":"P08195",
          "GCLM":"P48507","BACH1":"O14867","NFE2L2":"Q16236","KEAP1":"Q14145"}
    for g,u in spec.items():
        try:
            with urllib.request.urlopen(f"https://rest.uniprot.org/uniprotkb/{u}.json?fields=sequence,ft_domain,ft_region,ft_binding",timeout=30) as r:
                d=js.load(r)
            TARGETS[g]=dict(uniprot=u, seq=d["sequence"]["value"],
                            feats=[(f["type"], f["location"]["start"]["value"], f["location"]["end"]["value"],
                                    f.get("description","")) for f in d.get("features",[])])
            log(f"  {g} ({u}) {len(TARGETS[g]['seq'])} aa")
        except Exception as e:
            log(f"  {g} 取得失敗: {e}")

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--genes", nargs="+", default=["MAFG","AIFM2","SLC3A2","GCLM","BACH1","NFE2L2","KEAP1"])
    ap.add_argument("--model", default="facebook/esm2_t33_650M_UR50D")
    ap.add_argument("--gpu", default="3")
    ap.add_argument("--batch", type=int, default=8)
    a=ap.parse_args()
    os.environ["CUDA_VISIBLE_DEVICES"]=a.gpu
    import torch
    from transformers import AutoTokenizer, AutoModelForMaskedLM
    dev="cuda" if torch.cuda.is_available() else "cpu"
    log("device", dev, "model", a.model)
    tok=AutoTokenizer.from_pretrained(a.model)
    model=AutoModelForMaskedLM.from_pretrained(a.model, torch_dtype=torch.float32).to(dev).eval()

    load_targets()
    allrows=[]
    for g in a.genes:
        if g not in TARGETS: continue
        seq=TARGETS[g]["seq"]; L=len(seq)
        log(f"\n=== {g}: {L} aa，{L*19} 個單點突變 ===")
        ids=tok(seq, return_tensors="pt")["input_ids"].to(dev)   # 1 x (L+2)
        mask_id=tok.mask_token_id
        logprobs=np.zeros((L, len(AA)), dtype=np.float32)
        aa_idx={aa: tok.convert_tokens_to_ids(aa) for aa in AA}
        with torch.no_grad():
            for s in range(0, L, a.batch):
                e=min(s+a.batch, L)
                batch=ids.repeat(e-s,1).clone()
                for k,i in enumerate(range(s,e)):
                    batch[k, i+1]=mask_id
                out=model(input_ids=batch).logits        # B x (L+2) x V
                lp=torch.log_softmax(out, dim=-1)
                for k,i in enumerate(range(s,e)):
                    row=lp[k, i+1]
                    logprobs[i]=np.array([row[aa_idx[aa]].item() for aa in AA])
                if (s//a.batch) % 20 == 0:
                    log(f"   {e}/{L}")
        rows=[]
        for i,wt in enumerate(seq):
            if wt not in aa_idx: continue
            wtlp=logprobs[i, AA.index(wt)]
            for j,mt in enumerate(AA):
                if mt==wt: continue
                rows.append(dict(gene=g, pos=i+1, wt=wt, mt=mt,
                                 llr=float(logprobs[i,j]-wtlp)))
        df=pd.DataFrame(rows)
        df.to_csv(f"{RES}/71_esm_dms_{g}.csv.gz", index=False, compression="gzip")
        # 位置層級的約束度 = 該位置所有突變 LLR 的平均（越負越不能動）
        pos=df.groupby(["pos","wt"]).llr.mean().reset_index().rename(columns={"llr":"mean_llr"})
        pos["gene"]=g
        pos["constraint_rank"]=pos.mean_llr.rank()
        allrows.append(pos)
        log(f"  最受約束的 20 個位置：")
        t=pos.sort_values("mean_llr").head(20)
        log("   " + " ".join(f"{r.wt}{r.pos}({r.mean_llr:.1f})" for _,r in t.iterrows()))
        log(f"  最可容忍的 10 個位置：")
        t=pos.sort_values("mean_llr").tail(10)
        log("   " + " ".join(f"{r.wt}{r.pos}({r.mean_llr:.1f})" for _,r in t.iterrows()))
    if allrows:
        P=pd.concat(allrows, ignore_index=True)
        P.to_csv(f"{RES}/71_esm_position_constraint.csv", index=False)
        log(f"\n寫出 {RES}/71_esm_position_constraint.csv（{len(P)} 個位置）")
    json.dump({g:{"uniprot":v["uniprot"],"seq":v["seq"]} for g,v in TARGETS.items()},
              open(f"{RES}/71_target_sequences.json","w"), indent=1)
    log("DONE")

main()
