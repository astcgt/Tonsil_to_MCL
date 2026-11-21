#!/usr/bin/env python
"""
步驟 55：對新生調控元件做「掃窗式 in-silico 突變」(window-scan ISM)。

為什麼用掃窗而不是單鹼基
------------------------
單鹼基 ISM 在 2 kb 元件上要跑 6,000 次前向傳遞，太慢。
改成把元件切成 20 bp 的窗，逐窗做保留組成的打散，
找出「哪一段序列真的在驅動」，再回頭在那幾段上掃 motif。
這能定位到 20 bp 解析度，足以落在單一 TF 結合位點上。

標的
----
MAFG   chr17 新生調控域（19 個元件，正常相關性 ≈0 → 腫瘤 0.28–0.53）
SLC7A11 chr4 新生元件（5 個）
"""
import os, sys, json, argparse
import numpy as np, pandas as pd

BASE="/raid/jupyter-alz.ee09/Tonsil_to_MCL"; RES=f"{BASE}/results"; GEN=f"{BASE}/data/genome"
FASTA=f"{GEN}/GRCh38.primary_assembly.genome.fa"
SEQ_LEN=524_288; HALF=SEQ_LEN//2; BIN=32; OUT_BINS=6144; OUT_SPAN=OUT_BINS*BIN

TARGETS = {
 "MAFG":   dict(chrom="chr17", tss=81927000, strand="-",
                elements=[(81831639,81833881),(81890554,81892436),(81926271,81926994),
                          (81927413,81928897),(81901674,81903076),(81859752,81861678),
                          (82022368,82024013),(82050834,82052437)]),
 "SLC7A11":dict(chrom="chr4",  tss=138242000, strand="+",
                elements=[(138341498,138342345),(138392601,138393029),
                          (138222944,138223474),(138348035,138348490),(138406192,138406449)]),
}

def log(*a): print("[55]", *a, flush=True)

def di_shuffle(seq, rng):
    pairs=[seq[i:i+2] for i in range(0,len(seq)-1,2)]
    rng.shuffle(pairs)
    out="".join(pairs)
    return (out+seq[len(out):])[:len(seq)]

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--gpu", default="4"); ap.add_argument("--win", type=int, default=20)
    ap.add_argument("--out", default=f"{RES}/55_ism_results.csv")
    a=ap.parse_args()
    os.environ["CUDA_VISIBLE_DEVICES"]=a.gpu
    import torch, pyfaidx, grelu.resources
    torch.backends.cudnn.enabled=False
    log("torch", torch.__version__, "cuda", torch.cuda.is_available())
    model=grelu.resources.load_model(repo_id="Genentech/borzoi-model", filename="human_rep0.ckpt")
    model=model.eval().to("cuda")
    tasks=pd.DataFrame(model.data_params["tasks"])
    txt=tasks.astype(str).agg(" ".join, axis=1)
    sel={"RNA_GM12878": np.flatnonzero(txt.str.contains("RNA:GM12878",case=False).values),
         "ATAC_B": np.flatnonzero(txt.str.contains("ATAC:Memory B|ATAC:Plasma B",case=False,regex=True).values),
         "DNASE_GM12878": np.flatnonzero(txt.str.contains("DNASE:GM12878",case=False).values)}
    fa=pyfaidx.Fasta(FASTA, as_raw=True, sequence_always_upper=True)
    rng=np.random.default_rng(0)

    def undo(x):
        cs=384.0; y=np.array(x,dtype=np.float64); m=y>cs
        y[m]=(y[m]-cs)**2+cs; return y**2
    def predict(seq):
        with torch.no_grad():
            p=model.predict_on_seqs(seq, device=0)
        return undo(np.asarray(p).squeeze())

    rows=[]
    for gene,cfg in TARGETS.items():
        chrom, tss = cfg["chrom"], cfg["tss"]
        ws = tss-HALF
        wt = str(fa[chrom][ws:ws+SEQ_LEN]); assert len(wt)==SEQ_LEN
        out_start = ws+(SEQ_LEN-OUT_SPAN)//2
        tb = (tss-out_start)//BIN
        pw = predict(wt)
        def sc(P):
            lo,hi=max(0,tb-64),min(OUT_BINS,tb+64)
            return {k: float(P[v,lo:hi].sum()) if len(v) else np.nan for k,v in sel.items()}
        base=sc(pw)
        log(f"\n=== {gene} {chrom}:{tss} WT: " + " ".join(f"{k}={v:.3g}" for k,v in base.items()))
        for (es,ee) in cfg["elements"]:
            if not (ws <= es and ee <= ws+SEQ_LEN):
                log(f"  元件 {es}-{ee} 超出視窗，跳過"); continue
            # 先整段刪
            i0,i1=es-ws, ee-ws
            mut=wt[:i0]+di_shuffle(wt[i0:i1],rng)+wt[i1:]
            s=sc(predict(mut))
            d_whole={k: np.log2((s[k]+1)/(base[k]+1)) for k in s}
            log(f"  [{gene}] 整段刪 {chrom}:{es}-{ee} ({ee-es} bp) → " +
                " ".join(f"{k}={d_whole[k]:+.3f}" for k in d_whole))
            rows.append(dict(gene=gene, chrom=chrom, el_start=es, el_end=ee,
                             win_start=np.nan, win_end=np.nan, mode="whole_element", **d_whole))
            # 掃窗（只對有明顯效應的元件做，省算力）
            if abs(d_whole["RNA_GM12878"]) < 0.08:
                continue
            log(f"    → 效應明顯，開始 {a.win} bp 掃窗（{(ee-es)//a.win} 個窗）")
            for p in range(es, ee-a.win+1, a.win):
                j0,j1=p-ws, p+a.win-ws
                m2=wt[:j0]+di_shuffle(wt[j0:j1],rng)+wt[j1:]
                s2=sc(predict(m2))
                d2={k: np.log2((s2[k]+1)/(base[k]+1)) for k in s2}
                rows.append(dict(gene=gene, chrom=chrom, el_start=es, el_end=ee,
                                 win_start=p, win_end=p+a.win, mode="window", **d2))
            df=pd.DataFrame(rows); df.to_csv(a.out, index=False)
    df=pd.DataFrame(rows); df.to_csv(a.out, index=False)
    log(f"\n寫出 {a.out}，共 {len(df)} 列")
    w=df[df["mode"]=="window"].copy()
    if len(w):
        pd.set_option("display.width",200)
        for g,sub in w.groupby("gene"):
            log(f"\n=== {g}: 影響最大的 15 個 20 bp 窗 ===")
            log("\n"+sub.reindex(sub.RNA_GM12878.abs().sort_values(ascending=False).index)
                .head(15)[["chrom","win_start","win_end","RNA_GM12878","ATAC_B","DNASE_GM12878"]]
                .round(4).to_string(index=False))
    log("DONE")

main()
