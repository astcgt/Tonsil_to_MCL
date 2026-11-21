#!/usr/bin/env python
"""
步驟 91：用 AlphaGenome 2 對「MCL 新生 TRE/MARE 元件」做獨立驗證與 in-silico 擾動。

為什麼這一步是必要的
--------------------
第 2 節的新生元件（SLC7A11 五個、MAFG 十九個）完全是從我們自己的兩位病人資料
推出來的。AlphaGenome 2 是在正常組織上訓練的外部模型，它提供兩個我們自己資料
做不到的檢定：

  Part 1（獨立性檢定）：這些元件在**正常 B 細胞**裡被預測成開放還是關閉？
      若預測為關閉 → 支持「它們是癌症特異的」（我們自己的資料也這麼說，但這是外部證據）。

  Part 2（因果檢定）：把 TRE/MARE 核心 (TGACTCA) 打掉，模型預測的
      SLC7A11 / MAFG 表現與可及性會不會掉？用 score_ism_variants 做全飽和突變，
      **同一個窗格內的非 motif 位置就是內部對照**。

  Part 3（讀出有效性檢定）：對 HMOX1 的經典 ARE 做同一套流程。
      HMOX1 是教科書等級的 NRF2/small-Maf 標的。如果連它都測不到，
      那 Part 2 的陰性結果只代表「讀出無效」，不代表生物學陰性。
      ——這是從 Borzoi 那次學到的教訓：沒有陽性對照的陰性結果沒有價值。

金鑰處理：只從 .secrets/alphagenome_key 讀取，不寫入任何輸出。
"""
import os, sys, json, re
import numpy as np
import pandas as pd
import pysam

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
RES = f"{BASE}/results"
FA = f"{BASE}/data/genome/GRCh38.primary_assembly.genome.fa"

with open(f"{BASE}/.secrets/alphagenome_key") as fh:
    _KEY = fh.read().strip()

from alphagenome.data import genome
from alphagenome.models.v2 import dna_client, variant_scorers, dna_output

SEQLEN = 524288                      # 500 KB 輸入窗（涵蓋 ±262 kb，夠蓋所有元件）
ISM_HALF = 16                        # MARE 核心兩側各 16 bp → 33 bp 窗
# 三個 B 細胞條件：blood B cell 的 ATAC 與 RNA-seq 都是 observed（非 imputed）
CONDS = {
    "blood B cell":          "UBERON:0000178_CL:0000236",
    "blood naive B cell":    "UBERON:0000178_CL:0000788",
    "blood memory B cell":   "UBERON:0000178_CL:0000787",
}
TSS = {"SLC7A11": ("chr4", 138242349, "-"),
       "MAFG":    ("chr17", 81927735, "-"),
       "HMOX1":   ("chr22", 35380361, "+")}


def log(*a):
    print("[91]", *a, flush=True)


def rc(s):
    return s[::-1].translate(str.maketrans("ACGTacgt", "TGCAtgca"))


model = dna_client.create(_KEY)
fa = pysam.FastaFile(FA)
log("AlphaGenome 2 連線建立（金鑰未輸出）")


# ---------------------------------------------------------------- Part 1
def part1():
    """新生元件在正常 B 細胞裡被預測成開放還是關閉？"""
    d = pd.read_csv(f"{RES}/52_rewiring_annotated.csv.gz")
    out_rows = []
    for gene in ["SLC7A11", "MAFG"]:
        chrom, tss, _ = TSS[gene]
        iv = genome.Interval(chrom, tss - SEQLEN // 2, tss + SEQLEN // 2)
        els = d[(d.gene == gene) & d.gained].sort_values("dr", ascending=False)
        els = els[(els.start > iv.start) & (els.end < iv.end)]
        # 同一基因、同窗格內「沒有重接線」的 peak 作為對照
        ctl = d[(d.gene == gene) & (~d.gained) & (~d.lost)
                & (d.start > iv.start) & (d.end < iv.end)]
        for cname, curie in CONDS.items():
            pred = model.predict_interval(
                interval=iv, requested_outputs=[dna_output.OutputType.ATAC],
                condition=curie)
            at = pred.atac
            # 只取 observed 的軌道（避免用 imputed 下結論）
            im = at.metadata.get("imputation_mode")
            obs = (np.where(np.asarray(im) == "observed")[0]
                   if im is not None else np.arange(at.values.shape[1]))
            if len(obs) == 0:
                log(f"  {gene} / {cname}: 無 observed ATAC 軌道，跳過")
                continue
            vals = at.values[:, obs].mean(axis=1)
            res = at.resolution
            def bin_mean(s, e):
                i0 = max(0, (s - iv.start) // res); i1 = min(len(vals), (e - iv.start) // res + 1)
                return float(vals[i0:i1].mean()) if i1 > i0 else np.nan

            def fixed501(s, e):
                """⚠️ 關鍵修正：以 peak 中點為中心取固定 501 bp。

                第一版用 bin_mean（peak 全長平均）得到 MAFG 20.95×／p=1.1e-5，
                但那是假象 —— peak 寬度與預測值的 Spearman 是 +0.73，而新生元件
                比對照寬 3.5 倍（1,402 vs 406 bp），所以「平均 over 更寬的區域」
                本身就會拿到更高的值。固定寬度讓寬度在設計上不可能成為混淆，
                而且 501 bp 正好對齊 AlphaGenome 自己 CenterMaskScorer 的寬度。
                """
                mid = (s + e) // 2
                return bin_mean(mid - 250, mid + 251)

            gen_bg = float(np.median(vals))
            for kind, sub in [("新生元件", els), ("同窗格未重接線對照", ctl.head(40))]:
                for _, r in sub.iterrows():
                    out_rows.append(dict(
                        gene=gene, cond=cname, kind=kind,
                        peak_id=int(r.peak_id), dr=round(r.dr, 3), dist=int(r.dist),
                        width=int(r.end) - int(r.start),
                        pred_atac=bin_mean(int(r.start), int(r.end)),      # 舊（有寬度混淆）
                        pred_atac_501=fixed501(int(r.start), int(r.end)),  # 新（固定寬度）
                        window_median=gen_bg))
            log(f"  {gene} / {cname}: 新生 {len(els)} 個、對照 {min(len(ctl),40)} 個")
    p1 = pd.DataFrame(out_rows)
    p1.to_csv(f"{RES}/91_ag_part1_baseline.csv", index=False)
    from scipy.stats import mannwhitneyu, spearmanr
    summ = []
    for col, tag in [("pred_atac", "peak 全長平均（舊，有寬度混淆）"),
                     ("pred_atac_501", "固定 501 bp（修正版）")]:
        log(f"\n=== Part 1 讀出：{tag} ===")
        for (g, c), s in p1.groupby(["gene", "cond"]):
            a = s[s.kind == "新生元件"]
            b = s[s.kind == "同窗格未重接線對照"]
            va, vb = a[col].dropna(), b[col].dropna()
            if len(va) < 2 or len(vb) < 2:
                continue
            u = mannwhitneyu(va, vb, alternative="two-sided")
            # 寬度是否還是混淆？
            rw = spearmanr(s.width, s[col], nan_policy="omit")
            # 最保守對照：取最寬的 10 個
            btop = b.nlargest(10, "width")[col].dropna()
            u2 = mannwhitneyu(va, btop, alternative="two-sided") if len(btop) >= 3 else None
            log(f"  {g:9s} {c:20s} 新生 {va.median():.4f} (n={len(va)}) vs 對照 {vb.median():.4f} "
                f"(n={len(vb)}) = {va.median()/max(vb.median(),1e-9):5.2f}×  p={u.pvalue:.3g}"
                f"  | 寬度相關 ρ={rw.statistic:+.2f}"
                + (f"  | 最寬10對照 {btop.median():.4f} = {va.median()/max(btop.median(),1e-9):.2f}× "
                   f"p={u2.pvalue:.3g}" if u2 else ""))
            summ.append(dict(readout=tag, gene=g, cond=c,
                             novel=round(float(va.median()), 5), ctrl=round(float(vb.median()), 5),
                             ratio=round(float(va.median() / max(vb.median(), 1e-9)), 3),
                             p=float(u.pvalue), width_rho=round(float(rw.statistic), 3),
                             ctrl_widest10=round(float(btop.median()), 5) if u2 else np.nan,
                             ratio_widest10=round(float(va.median() / max(btop.median(), 1e-9)), 3) if u2 else np.nan,
                             p_widest10=float(u2.pvalue) if u2 else np.nan,
                             n_novel=len(va), n_ctrl=len(vb)))
    pd.DataFrame(summ).to_csv(f"{RES}/91_ag_part1_summary.csv", index=False)
    return p1


# ---------------------------------------------------------------- Part 2/3
def find_mare(chrom, lo, hi):
    """在區間內找 TRE/MARE 核心（兩股），回傳正股座標。"""
    seq = fa.fetch(chrom, lo, hi).upper()
    hits = []
    for strand, s in [("+", seq), ("-", rc(seq))]:
        for mm in re.finditer("TGACTCA", s):
            off = mm.start() if strand == "+" else len(seq) - mm.end()
            hits.append((lo + off, strand))
    return sorted(set(hits))


def ism_site(label, gene, chrom, center, strand_note):
    """對一個 MARE 位點做全飽和 in-silico 突變。"""
    g_chrom, tss, _ = TSS[gene]
    iv = genome.Interval(g_chrom, tss - SEQLEN // 2, tss + SEQLEN // 2)
    ism_iv = genome.Interval(chrom, center - ISM_HALF, center + 7 + ISM_HALF)
    sc = [variant_scorers.RECOMMENDED_VARIANT_SCORERS["RNA_SEQ"],
          variant_scorers.RECOMMENDED_VARIANT_SCORERS["ATAC"]]
    log(f"  ISM {label}: {chrom}:{ism_iv.start}-{ism_iv.end} ({ism_iv.width} bp) "
        f"在 {gene} 的 500 kb 窗內")
    rows = []
    for cname, curie in CONDS.items():
        try:
            res = model.score_ism_variants(interval=iv, ism_interval=ism_iv,
                                           condition=curie, variant_scorers=sc,
                                           progress_bar=False, max_workers=5)
        except Exception as e:
            log(f"    {cname} 失敗: {type(e).__name__}: {e}")
            continue
        for per_variant in res:
            for ad in per_variant:
                # ⚠️ score_ism_variants 回傳的 AnnData 在 uns 裡沒有 variant_scorer，
                # 所以不能靠字串判斷讀出型別（第一版因此把 RNA 分數誤標成 ATAC）。
                # 改用結構特徵：GeneMask（RNA）結果的 obs 有 gene_name，CenterMask（ATAC）沒有。
                is_gene = "gene_name" in ad.obs.columns
                vstr = str(ad.uns.get("variant", ""))
                mpos = int(re.search(r":(\d+)", vstr).group(1)) if re.search(r":(\d+)", vstr) else np.nan
                im = ad.var.get("imputation_mode")
                obs = (np.where(np.asarray(im) == "observed")[0]
                       if im is not None else np.arange(ad.shape[1]))
                if len(obs) == 0:
                    continue
                X = np.asarray(ad.X)
                if is_gene:
                    gi = [i for i, o in enumerate(ad.obs.gene_name.values) if o == gene]
                    if not gi:
                        continue
                    val = float(np.nanmean(X[np.ix_(gi, obs)]))
                    kind = "RNA"
                else:
                    val = float(np.nanmean(X[:, obs]))
                    kind = "ATAC"
                rows.append(dict(site=label, gene=gene, cond=cname, readout=kind,
                                 variant=vstr, pos=mpos, score=val,
                                 in_motif=bool(center <= mpos < center + 7)))
    return rows


def part23():
    sites = []
    # SLC7A11：兩個新生元件中的 TRE 核心
    for pos, st in [(138406332, "-"), (138348373, "+")]:
        sites.append((f"SLC7A11_el_{pos}", "SLC7A11", "chr4", pos, st))
    # MAFG：啟動子區（−769 bp）的 TGACTCAGC（TRE + ARE 重疊）
    sites.append(("MAFG_promoter_81926966", "MAFG", "chr17", 81926966, "+"))
    # Part 3 陽性對照：HMOX1 上游的經典 ARE
    hm = find_mare("chr22", 35380361 - 12000, 35380361 + 2000)
    log(f"HMOX1 上游 ±12kb 找到 {len(hm)} 個 TRE 核心: {hm}")
    for pos, st in hm[:3]:
        sites.append((f"HMOX1_ARE_{pos}", "HMOX1", "chr22", pos, st))

    allrows = []
    for label, gene, chrom, pos, st in sites:
        allrows += ism_site(label, gene, chrom, pos, st)
    p2 = pd.DataFrame(allrows)
    if p2.empty:
        log("ISM 無結果"); return p2
    p2.to_csv(f"{RES}/91_ag_ism.csv.gz", index=False, compression="gzip")
    log(f"\n=== Part 2/3：ISM 結果（motif 內 vs 同窗格 motif 外）===")
    from scipy.stats import mannwhitneyu
    summ = []
    for (site, cond, ro), s in p2.groupby(["site", "cond", "readout"]):
        a = s[s.in_motif].score.dropna()
        b = s[~s.in_motif].score.dropna()
        if len(a) < 3 or len(b) < 3:
            continue
        u = mannwhitneyu(a.abs(), b.abs(), alternative="greater")
        summ.append(dict(site=site, cond=cond, readout=ro,
                         n_motif=len(a), n_flank=len(b),
                         motif_med_abs=round(float(a.abs().median()), 5),
                         flank_med_abs=round(float(b.abs().median()), 5),
                         ratio=round(float(a.abs().median() / max(b.abs().median(), 1e-12)), 2),
                         motif_min=round(float(a.min()), 5),
                         p=float(u.pvalue)))
    sm = pd.DataFrame(summ)
    # 18–36 個檢定一起做 Bonferroni；未校正的 p 不足以下結論
    sm["p_bonf"] = np.minimum(1.0, sm.p * len(sm))
    sm = sm.sort_values(["readout", "site", "cond"])
    sm.to_csv(f"{RES}/91_ag_ism_summary.csv", index=False)
    pd.set_option("display.width", 220)
    log("\n" + sm.to_string(index=False))
    return p2


if __name__ == "__main__":
    if os.environ.get("SKIP_P1") != "1":
        log("=== Part 1：新生元件在正常 B 細胞的預測可及性 ===")
        try:
            part1()
        except Exception as e:
            log(f"Part 1 失敗: {type(e).__name__}: {e}")
    log("\n=== Part 2/3：MARE 核心的 in-silico 飽和突變 ===")
    try:
        part23()
    except Exception as e:
        import traceback; traceback.print_exc()
        log(f"Part 2/3 失敗: {type(e).__name__}: {e}")
    log("DONE")
