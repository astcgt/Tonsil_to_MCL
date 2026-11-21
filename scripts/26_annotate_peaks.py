#!/usr/bin/env python
"""
步驟 26：把差異可及性的 peak 標註到最近的基因，並與表現量殘差整合。

輸出
----
results/26_peak_annotated.csv      每個 peak 的最近基因與距離
results/26_gene_level_evidence.csv 基因層級的整合證據表
"""
import os, sys, gzip
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
RES = f"{BASE}/results"
GTF = f"{BASE}/data/genome/gencode.v44.annotation.gtf.gz"


def log(*a):
    print("[26]", *a, flush=True)


def load_tss(gtf):
    rows = []
    with gzip.open(gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if f[2] != "gene":
                continue
            a = f[8]
            if 'gene_name "' not in a:
                continue
            nm = a.split('gene_name "')[1].split('"')[0]
            gt = a.split('gene_type "')[1].split('"')[0] if 'gene_type "' in a else ""
            tss = int(f[3]) if f[6] == "+" else int(f[4])
            rows.append((f[0], tss, nm, gt, int(f[3]), int(f[4])))
    df = pd.DataFrame(rows, columns=["chrom", "tss", "gene", "gene_type", "gstart", "gend"])
    return df


tss = load_tss(GTF)
log("GTF 基因:", len(tss))
tssp = tss[tss.gene_type.isin(["protein_coding", "lncRNA"])].copy()
log("蛋白質編碼 + lncRNA:", len(tssp))

d = pd.read_csv(f"{RES}/22_atac_diff.csv")
log("peaks:", len(d))
d["mid"] = ((d.start + d.end) / 2).astype(int)

# 每條染色體各自做最近鄰搜尋
out = []
for c, sub in d.groupby("chrom"):
    g = tssp[tssp.chrom == c].sort_values("tss")
    if len(g) == 0:
        sub = sub.copy(); sub["nearest_gene"] = None; sub["dist_to_tss"] = np.nan
        out.append(sub); continue
    pos = g.tss.values
    idx = np.searchsorted(pos, sub["mid"].values)
    idx = np.clip(idx, 1, len(pos) - 1)
    left, right = pos[idx - 1], pos[idx]
    pick = np.where(np.abs(sub["mid"].values - left) <= np.abs(right - sub["mid"].values),
                    idx - 1, idx)
    sub = sub.copy()
    sub["nearest_gene"] = g.gene.values[pick]
    sub["dist_to_tss"] = sub["mid"].values - pos[pick]
    out.append(sub)
ann = pd.concat(out, ignore_index=True)
ann.to_csv(f"{RES}/26_peak_annotated.csv", index=False)
log("寫出 26_peak_annotated.csv")

sig = ann[(ann.consistent) & (ann.mcl_mean_log2cpm > 3)].sort_values("lfc_mean", ascending=False)
log("\n=== 新增可及性最強的 40 個 peak（含最近基因）===")
log("\n" + sig.head(40)[["peak", "nearest_gene", "dist_to_tss",
                         "mcl_mean_log2cpm", "normal_log2cpm", "lfc_mean"]]
    .round(2).to_string(index=False))

# ---- 基因層級：把 ±100 kb 內的 peak 匯總 ----
near = ann[(ann.dist_to_tss.abs() < 100000)].copy()
gl = near.groupby("nearest_gene").agg(
    n_peaks=("peak", "size"),
    atac_lfc_mean=("lfc_mean", "mean"),
    atac_lfc_max=("lfc_mean", "max"),
    mcl_log2cpm_max=("mcl_mean_log2cpm", "max"),
    n_consistent=("consistent", "sum")).reset_index().rename(columns={"nearest_gene": "gene"})

ex = pd.read_csv(f"{RES}/15_gene_residual_tumor_specific.csv")
merged = ex.merge(gl, on="gene", how="left")

dep = pd.read_csv(f"{RES}/09_depmap_mcl_dependency.csv")
merged = merged.merge(dep[["gene", "mcl_mean", "mcl_max", "other_mean",
                           "delta_vs_all", "z_vs_all"]]
                      .rename(columns={"mcl_mean": "depmap_mcl_geneeffect",
                                       "mcl_max": "depmap_mcl_max",
                                       "other_mean": "depmap_other_mean",
                                       "delta_vs_all": "depmap_delta",
                                       "z_vs_all": "depmap_z"}),
                      on="gene", how="left")
merged.to_csv(f"{RES}/26_gene_level_evidence.csv", index=False)
log("\n寫出 26_gene_level_evidence.csv", merged.shape)

top = merged.sort_values("resid_tumor_specific", ascending=False).head(30)
log("\n=== 表現殘差 Top30 + 其染色質與依賴性證據 ===")
log("\n" + top[["gene", "resid_tumor_specific", "mean_observed", "mean_expected",
                "atac_lfc_mean", "atac_lfc_max", "n_peaks",
                "depmap_mcl_geneeffect", "depmap_z"]].round(3).to_string(index=False))
log("DONE")
