#!/usr/bin/env python
"""
步驟 100：甲基化讀出的修正版 —— 只在 CpG 位點取值。

步驟 97 的錯誤（我自己的，已在此修正）
--------------------------------------
METHYLATION 輸出是 **1 bp 解析度，而且只在 CpG 位點有值，其餘位置為 0**。
步驟 97 對整個 501 bp 視窗取平均，於是被非 CpG 的零稀釋了約 500/3 倍：
隨機背景算出 0.0070，看起來「幾乎沒有甲基化」，實際還原回去接近 1.0。
若照步驟 97 的數字下結論，會得到完全相反的解讀。

另外步驟 97 把兩個軌道平均了，但它們是兩種不同實驗
（WGBS 與 snmC-seq3），這一步分開報告。

修正後的讀出
------------
從參考基因組取出視窗內所有 CpG 位置，只在這些位置取預測值的平均。
這樣數字就是可解讀的甲基化比例（0–1）。

校準錨點（必須先過關，否則讀出不可信）
  CpG 島啟動子（GAPDH / ACTB）→ 應該接近 0（未甲基化）
  隨機基因體背景           → 應該接近 0.8（高甲基化）

要回答的問題
------------
SLC7A11 那五個新生元件，在正常 B 細胞裡是不是「被甲基化關著」？
若是 → 去甲基化藥物（azacitidine / decitabine，血液腫瘤常用）
反而可能打開這條鐵死亡防禦軸線，是個臨床警訊。
"""
import sys, os
import numpy as np
import pandas as pd

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
RES = f"{BASE}/results"
_KEY = open(f"{BASE}/.secrets/alphagenome_key").read().strip()

from alphagenome.data import genome
from alphagenome.models.v2 import dna_client, dna_output
import pyfaidx

RNG = np.random.default_rng(0)
TSS = {"SLC7A11": ("chr4", 138242349), "MAFG": ("chr17", 81927735)}
ANCHOR = {"GAPDH_promoter": ("chr12", 6534517), "ACTB_promoter": ("chr7", 5530601)}
COND = {"blood B cell": "UBERON:0000178_CL:0000236",
        "blood memory B cell": "UBERON:0000178_CL:0000787"}
SEQLEN = 524288
HALF = 250

FA = pyfaidx.Fasta(f"{BASE}/data/genome/GRCh38.primary_assembly.genome.fa",
                   as_raw=True, sequence_always_upper=True)


def log(*a):
    print("[100]", *a, flush=True)


def features(gene, chrom, tss, half):
    rw = pd.read_csv(f"{RES}/52_rewiring_annotated.csv.gz")
    g = rw[(rw.gene == gene) & (rw.chrom == chrom)]
    lo, hi = tss - half, tss + half
    out = []
    for kind, sub in [("de_novo", g[g.gained]), ("control", g[(~g.gained) & (~g.lost)])]:
        for _, r in sub.iterrows():
            mid = (int(r.start) + int(r.end)) // 2
            if lo + 5000 < mid < hi - 5000:
                out.append(dict(label=f"{kind}_{int(r.peak_id)}", kind=kind, mid=mid,
                                dist=int(mid - tss)))
    return out


model = dna_client.create(_KEY)
log("連線建立（金鑰未輸出）")

rows = []
for gene, (chrom, tss) in {**TSS, **ANCHOR}.items():
    iv = genome.Interval(chrom, tss - SEQLEN // 2, tss + SEQLEN // 2)
    if gene in TSS:
        feats = features(gene, chrom, tss, SEQLEN // 2)
    else:
        feats = [dict(label=f"anchor_{gene}", kind="anchor_cpg_island", mid=tss, dist=0)]
    for m in RNG.integers(iv.start + 2000, iv.end - 2000, 80):
        feats.append(dict(label=f"random_{m}", kind="random_background",
                          mid=int(m), dist=int(m - tss)))

    # 視窗內所有 CpG 的位置（以 C 的座標為準）
    seq = str(FA[chrom][iv.start:iv.end])
    cpg = np.frombuffer(seq.encode(), dtype="S1")
    is_cg = (cpg[:-1] == b"C") & (cpg[1:] == b"G")
    cpg_idx = np.flatnonzero(is_cg)
    log(f"\n=== {gene} | 視窗 CpG {len(cpg_idx):,} 個 | 特徵 {len(feats)} 個 ===")

    for cname, cid in COND.items():
        try:
            pr = model.predict_interval(interval=iv, condition=cid,
                                        requested_outputs=[dna_output.OutputType.METHYLATION])
        except Exception as e:
            log(f"  {cname} 失敗: {type(e).__name__}: {str(e)[:150]}"); continue
        me = getattr(pr, "methylation", None)
        if me is None:
            log(f"  {cname}: 無 methylation"); continue
        V = np.asarray(me.values)
        names = list(me.metadata["name"]) if hasattr(me, "metadata") else \
            [f"track{i}" for i in range(V.shape[-1])]
        assays = [n.split("_")[-1] for n in names]
        log(f"  {cname}: shape {V.shape}, 解析度 {me.resolution} bp, 軌道 {assays}")

        for ti, assay in enumerate(assays):
            v = V[:, ti]
            for f in feats:
                c0 = f["mid"] - HALF - iv.start
                c1 = f["mid"] + HALF - iv.start
                sel = cpg_idx[(cpg_idx >= max(c0, 0)) & (cpg_idx <= min(c1, len(v) - 1))]
                if len(sel) == 0:
                    val, n_cpg = np.nan, 0
                else:
                    s = v[sel]
                    s = s[np.isfinite(s)]
                    val, n_cpg = (float(s.mean()) if len(s) else np.nan), len(sel)
                rows.append(dict(gene=gene, condition=cname, assay=assay,
                                 label=f["label"], kind=f["kind"], dist=f["dist"],
                                 n_cpg=n_cpg, methyl=val))

df = pd.DataFrame(rows)
df.to_csv(f"{RES}/100_methylation_cpg.csv", index=False)

from scipy.stats import mannwhitneyu
pd.set_option("display.width", 200)

log("\n" + "=" * 76)
log("=== 校準：讀出可信嗎？（CpG 島啟動子應接近 0、隨機背景應接近 0.8）===")
for (cond, assay), s in df.groupby(["condition", "assay"]):
    an = s[s.kind == "anchor_cpg_island"]
    bg = s[s.kind == "random_background"]
    log(f"{cond:22s} {assay:10s} | CpG 島錨點 {an.methyl.mean():.3f} "
        f"(n={len(an)}) | 隨機背景中位 {bg.methyl.median():.3f} (n={len(bg)}, "
        f"CpG 中位 {bg.n_cpg.median():.0f})")

log("\n" + "=" * 76)
log("=== 問題：新生元件在正常 B 細胞裡被甲基化關著嗎？===")
for (gene, cond, assay), s in df[df.gene.isin(TSS)].groupby(["gene", "condition", "assay"]):
    dn = s[(s.kind == "de_novo") & s.methyl.notna()]
    cl = s[(s.kind == "control") & s.methyl.notna()]
    bg = s[s.kind == "random_background"]
    if len(dn) < 3 or len(cl) < 5:
        continue
    u = mannwhitneyu(dn.methyl, cl.methyl, alternative="greater")
    log(f"\n{gene} / {cond} / {assay}")
    log(f"  新生元件 {dn.methyl.median():.3f} (n={len(dn)}, CpG 中位 {dn.n_cpg.median():.0f})"
        f" vs 未重接線對照 {cl.methyl.median():.3f} (n={len(cl)})"
        f" | 隨機背景 {bg.methyl.median():.3f} | p = {u.pvalue:.3g}")
    log(dn[["label", "dist", "n_cpg", "methyl"]].sort_values("dist").round(3).to_string(index=False))
log("\nDONE")
