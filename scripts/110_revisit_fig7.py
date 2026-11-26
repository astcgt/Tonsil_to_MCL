#!/usr/bin/env python
"""
步驟 110：重現並重新檢視 Massoni-Badosa et al., Immunity 2024 的 MCL 分析（原文 Figure 7G–I）。

原文對 MCL 的主張（Results 最後一段、Figure 7G–I、Discussion）：
  (G) 兩位病人都有一個帶 chrY loss 的主要 cluster，可把腫瘤細胞分成 chrY+/chrY−（c1/c2）；
      inferCNV 以檢體內的非腫瘤 B 細胞（c3）為參考。
  (H) 每個 cluster 的細胞數，以及正常／腫瘤 B 細胞標記的平均表現。
  (I) 兩個 cluster 有升高的 metallothionein 表現。
  Discussion：「MCL cells are not frozen in a particular maturation stage,
      but they retain a certain differentiation potential of normal B cells.」

本腳本做四件事：
  A. 重現 chrY 亞克隆結構（原文 Fig 7G 上半 + S7K–S7N）
  B. 重現 inferCNV 的腫瘤／非腫瘤分離（原文 Fig 7G 下半）
  C. 重現 metallothionein cluster（原文 Fig 7I），並檢查「兩個 cluster」是否成立
  D. 把原文的定性結論（「不被凍結在某一個成熟階段」）換成定量檢定：
     每顆腫瘤細胞對 46 個正常狀態的軟性對位權重，算 Shannon 熵，
     並以同一套流程對檢體內殘存的正常 B 細胞做同樣計算當尺規。

輸入
  data/processed/mcl_tumor_labeled.h5ad   （20,794 顆 CD19+ 細胞，含 chrY_score/cnv_score/leiden）
  results/15_match_weights.csv.gz         （每顆細胞 × 46 個正常狀態的軟性對位權重）
  results/15_mcl_obs.csv                  （含 matched_normal_state）
輸出
  results/110_revisit_fig7.json           （所有數字）
  results/110_cluster_profile.csv         （每個 cluster 的細胞數、chrY、MT、標記表現）
  results/110_state_assignment.csv        （每個正常狀態被對位到的細胞數與比例）
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
PROC = f"{BASE}/data/processed"
RES = f"{BASE}/results"

sc.settings.n_jobs = 16


def log(*a):
    print(*a, flush=True)


# 原文 Fig 7H 的「正常與腫瘤 B 細胞標記」。分組依原子論文對 NBC/MBC/GC 區間的定義。
MARKERS = {
    "MCL (t(11;14))": ["CCND1", "SOX11"],
    "Naive": ["TCL1A", "IGHD", "FCER2", "SELL"],
    "NBC activation": ["CD69", "EGR2", "EGR3", "MYC", "CCL3", "CCL4"],
    "GC-committed": ["MEF2B", "RGS13", "CCND2", "MIR155HG", "TRAF1", "BCL6"],
    "GC": ["AICDA", "CXCR4", "CD83", "MKI67"],
    "Memory": ["CD27", "TNFRSF13B", "AIM2"],
    "Plasma cell": ["PRDM1", "XBP1"],
}
MT_GENES = ["MT1X", "MT2A", "MT1E", "MT1F", "MT1G", "MT1M", "MT1H", "MT1A"]

# 原文把 MCL 狀態描述為「from activated NBCs to GC-committed cells」。
# 下面是參考註釋裡落在那個區間的狀態名稱（annotation_20230508 的字串）。
INTERVAL_ACTIVATED_TO_GCCOMMITTED = [
    "NBC early activation",
    "NBC IFN-activated",
    "NBC CD229+",
    "Proliferative NBC",
    "Early GC-commited NBC",
    "GC-commited NBC",
    "preGC",
]
MEMORY_STATES = [
    "ncsMBC",
    "csMBC",
    "ncsMBC FCRL4/5+",
    "csMBC FCRL4/5+",
    "MBC FCRL5+",
    "Early MBC",
    "Precursor MBCs",
    "Reactivated proliferative MBCs",
]

out = {}

# ---------------------------------------------------------------- 載入
log("載入 mcl_tumor_labeled.h5ad ...")
A = sc.read_h5ad(f"{PROC}/mcl_tumor_labeled.h5ad")
log("  ", A.shape)
A.obs["leiden"] = A.obs["leiden"].astype(str)

out["n_cells_cd19pos"] = int(A.n_obs)
out["n_tumor"] = int(A.obs.is_tumor.sum())
out["n_nontumor"] = int((~A.obs.is_tumor).sum())
out["per_donor"] = {
    str(d): {
        "n": int((A.obs.donor_id == d).sum()),
        "n_tumor": int(((A.obs.donor_id == d) & A.obs.is_tumor).sum()),
        "n_nontumor": int(((A.obs.donor_id == d) & ~A.obs.is_tumor).sum()),
    }
    for d in sorted(A.obs.donor_id.unique())
}
log("  細胞數", out["per_donor"])

# ------------------------------------------------- A. chrY 亞克隆（Fig 7G 上）
log("\n[A] chrY 亞克隆結構")
chry = {}
for d in sorted(A.obs.donor_id.unique()):
    m = (A.obs.donor_id == d) & A.obs.is_tumor
    s = A.obs.loc[m, "chrY_score"].to_numpy()
    # 雙峰性檢定：單一高斯 vs 兩個高斯的 BIC 差
    from sklearn.mixture import GaussianMixture

    x = s.reshape(-1, 1)
    bic1 = GaussianMixture(1, random_state=0).fit(x).bic(x)
    g2 = GaussianMixture(2, random_state=0).fit(x)
    bic2 = g2.bic(x)
    lab = g2.predict(x)
    mu = g2.means_.ravel()
    lo = int(np.argmin(mu))  # chrY 低 = chrY loss
    frac_loss = float((lab == lo).mean())
    # 兩群的分離度（Cohen's d）
    a, b = s[lab == lo], s[lab != lo]
    d_coh = float(
        (b.mean() - a.mean())
        / np.sqrt(((a.var(ddof=1) * (len(a) - 1)) + (b.var(ddof=1) * (len(b) - 1))) / (len(a) + len(b) - 2))
    )
    chry[str(d)] = {
        "n_tumor": int(m.sum()),
        "bic_1component": float(bic1),
        "bic_2component": float(bic2),
        "delta_bic": float(bic1 - bic2),
        "mean_chrY_loss": float(mu[lo]),
        "mean_chrY_retained": float(mu[1 - lo]),
        "frac_chrY_loss": frac_loss,
        "cohens_d": d_coh,
    }
    log(f"  {d}: chrY− 占 {frac_loss:.1%}  ΔBIC={bic1-bic2:.0f}  d={d_coh:.2f}")
out["chrY"] = chry

# chrY 是否與 cluster 結構對齊（原文說「a major cluster with chromosome Y loss」）
cl_chry = (
    A.obs[A.obs.is_tumor]
    .groupby(["donor_id", "leiden"], observed=True)["chrY_score"]
    .agg(["mean", "size"])
    .reset_index()
)
cl_chry = cl_chry[cl_chry["size"] >= 50]
out["chrY_cluster_spread"] = {
    str(d): {
        "n_clusters": int((cl_chry.donor_id == d).sum()),
        "min_cluster_mean": float(cl_chry.loc[cl_chry.donor_id == d, "mean"].min()),
        "max_cluster_mean": float(cl_chry.loc[cl_chry.donor_id == d, "mean"].max()),
    }
    for d in sorted(cl_chry.donor_id.unique())
}
log("  cluster 層級 chrY 範圍", out["chrY_cluster_spread"])

# ------------------------------------------- B. inferCNV 分離（Fig 7G 下）
log("\n[B] inferCNV 腫瘤／非腫瘤分離")
cnv = {}
for d in sorted(A.obs.donor_id.unique()):
    t = A.obs.loc[(A.obs.donor_id == d) & A.obs.is_tumor, "cnv_score"].to_numpy()
    n = A.obs.loc[(A.obs.donor_id == d) & ~A.obs.is_tumor, "cnv_score"].to_numpy()
    if len(n) < 10:
        cnv[str(d)] = {"n_nontumor": int(len(n)), "note": "too few in-sample normal B"}
        continue
    u = stats.mannwhitneyu(t, n, alternative="greater")
    cnv[str(d)] = {
        "n_tumor": int(len(t)),
        "n_nontumor": int(len(n)),
        "median_tumor": float(np.median(t)),
        "median_nontumor": float(np.median(n)),
        "fold": float(np.median(t) / np.median(n)),
        "mwu_p": float(u.pvalue),
        "auc": float(u.statistic / (len(t) * len(n))),
    }
    log(f"  {d}: 腫瘤 {np.median(t):.4f} vs 非腫瘤 {np.median(n):.4f}  AUC={cnv[str(d)]['auc']:.3f}  p={u.pvalue:.2e}")
out["infercnv"] = cnv

# 檢體內非腫瘤 B 與扁桃腺正常 B 參考的比較（原文用 c3 當 reference）
cnv_obs = pd.read_csv(f"{RES}/20_cnv_obs.csv", index_col=0)
ref_sc = cnv_obs.loc[cnv_obs.src == "ref", "cnv_score"]
p99 = float(ref_sc.quantile(0.99))
out["infercnv_ref"] = {
    "n_ref_normal_B": int(len(ref_sc)),
    "median_ref": float(ref_sc.median()),
    "p99_ref": p99,
}
# 檢體內殘存的正常 B 太少（M102 只有 29 顆），小樣本 AUC 不穩。
# 另外用 5,298 顆扁桃腺正常 B 當外部參考分布，報告超過其 99 百分位的腫瘤細胞比例。
for d in sorted(A.obs.donor_id.unique()):
    t = A.obs.loc[(A.obs.donor_id == d) & A.obs.is_tumor, "cnv_score"]
    n = A.obs.loc[(A.obs.donor_id == d) & ~A.obs.is_tumor, "cnv_score"]
    out["infercnv"][str(d)]["frac_tumor_above_ref_p99"] = float((t > p99).mean())
    out["infercnv"][str(d)]["frac_nontumor_above_ref_p99"] = float((n > p99).mean()) if len(n) else None
    log(f"  {d}: 腫瘤超過參考 p99 = {(t > p99).mean():.1%}"
        + (f"；檢體內非腫瘤 B = {(n > p99).mean():.1%}（n={len(n)}）" if len(n) else ""))
log("  扁桃腺正常 B 參考:", out["infercnv_ref"])

# -------------------------------------- C. metallothionein cluster（Fig 7I）
log("\n[C] metallothionein")
present = {k: [g for g in v if g in A.var_names] for k, v in MARKERS.items()}
mt_present = [g for g in MT_GENES if g in A.var_names]
out["mt_genes_used"] = mt_present
sc.tl.score_genes(A, mt_present, score_name="MT_score")
for k, gs in present.items():
    if gs:
        sc.tl.score_genes(A, gs, score_name=f"sc_{k}")

tum = A.obs[A.obs.is_tumor].copy()
prof = (
    tum.groupby(["donor_id", "leiden"], observed=True)
    .agg(n=("MT_score", "size"), MT_score=("MT_score", "mean"), chrY_score=("chrY_score", "mean"),
         cnv_score=("cnv_score", "mean"), CCND1=("CCND1", "mean"), SOX11=("SOX11", "mean"))
    .reset_index()
)
for k in present:
    if f"sc_{k}" in tum:
        prof[k] = (
            tum.groupby(["donor_id", "leiden"], observed=True)[f"sc_{k}"].mean().reset_index(drop=True).to_numpy()
        )
prof = prof[prof.n >= 50].sort_values(["donor_id", "MT_score"], ascending=[True, False])
prof.to_csv(f"{RES}/110_cluster_profile.csv", index=False)

# 「升高」的定義：逐 cluster 對同病人其餘腫瘤細胞做 Mann-Whitney，
# Bonferroni 校正（除以該病人的 cluster 數），再要求 rank-biserial 效應量 > 0.3。
# 這個定義不依賴我們的 cluster 分法要和原文一致，只問「有沒有 cluster 的 MT 明顯較高」。
# （不用離群值門檻：cluster 只有 8–9 個時，mean+2SD 的 z 上限是 (n-1)/sqrt(n) ≈ 2.5，
#   而 MT 在 cluster 之間是連續梯度而非兩個離群點，離群值檢定在此無法成立。）
mt_hits = {}
for d in sorted(prof.donor_id.unique()):
    p = prof[prof.donor_id == d]
    m = (A.obs.donor_id == d) & A.obs.is_tumor
    kept = set(p.leiden)
    per_cl = {}
    for cl in p.leiden:
        sel = (m & (A.obs.leiden == cl)).to_numpy()
        rest = (m & (A.obs.leiden != cl) & A.obs.leiden.isin(kept)).to_numpy()
        a_ = A.obs.loc[sel, "MT_score"].to_numpy()
        b_ = A.obs.loc[rest, "MT_score"].to_numpy()
        u = stats.mannwhitneyu(a_, b_, alternative="greater")
        rb = 2 * u.statistic / (len(a_) * len(b_)) - 1  # rank-biserial, -1..1
        per_cl[str(cl)] = {
            "n": int(len(a_)),
            "MT_score": float(a_.mean()),
            "mwu_p": float(u.pvalue),
            "rank_biserial": float(rb),
        }
    alpha = 0.05 / len(p)
    hit = [c for c, v in per_cl.items() if v["mwu_p"] < alpha and v["rank_biserial"] > 0.3]
    # Kruskal-Wallis：MT 是否被 cluster 結構解釋（而非均勻散布）
    groups = [A.obs.loc[(m & (A.obs.leiden == cl)).to_numpy(), "MT_score"].to_numpy() for cl in p.leiden]
    kw = stats.kruskal(*groups)
    n_tot = sum(len(g) for g in groups)
    eta2 = float((kw.statistic - len(groups) + 1) / (n_tot - len(groups)))
    mt_hits[str(d)] = {
        "n_clusters_tested": int(len(p)),
        "bonferroni_alpha": float(alpha),
        "n_elevated": len(hit),
        "elevated_clusters": hit,
        "elevated_n_cells": [per_cl[c]["n"] for c in hit],
        "elevated_MT_score": [per_cl[c]["MT_score"] for c in hit],
        "elevated_rank_biserial": [per_cl[c]["rank_biserial"] for c in hit],
        "max_MT_score": float(p.MT_score.max()),
        "min_MT_score": float(p.MT_score.min()),
        "median_MT_score": float(p.MT_score.median()),
        "cluster_spread": float(p.MT_score.max() - p.MT_score.min()),
        "kruskal_p": float(kw.pvalue),
        "kruskal_eta2": eta2,
        "all_cluster_scores": {str(r.leiden): float(r.MT_score) for r in p.itertuples()},
        "per_cluster": per_cl,
    }
    log(f"  {d}: {len(hit)}/{len(p)} 個 cluster MT 顯著升高 -> {hit}；"
        f"cluster 間全距 {p.MT_score.max()-p.MT_score.min():.2f}，Kruskal η²={eta2:.3f}")
out["metallothionein"] = mt_hits

# 病人層級：這個程式是否為病人特異（原文 Fig 7I 只展示 M102）
_mt_by_donor = {}
for d in sorted(A.obs.donor_id.unique()):
    v = A.obs.loc[((A.obs.donor_id == d) & A.obs.is_tumor).to_numpy(), "MT_score"].to_numpy()
    _mt_by_donor[str(d)] = v
if len(_mt_by_donor) == 2:
    (d1, v1), (d2, v2) = list(_mt_by_donor.items())
    u = stats.mannwhitneyu(v1, v2, alternative="two-sided")
    out["metallothionein_donor_contrast"] = {
        "donors": [d1, d2],
        "median": [float(np.median(v1)), float(np.median(v2))],
        "mwu_p": float(u.pvalue),
        "higher_in": d1 if np.median(v1) > np.median(v2) else d2,
    }
    log(f"  病人間對比：{d1} 中位數 {np.median(v1):.3f} vs {d2} {np.median(v2):.3f}, p={u.pvalue:.2e}")

# MT 高 cluster 的逐基因倍率（相對同病人其餘腫瘤細胞），用來確認不是單一基因撐起來的
import scipy.sparse as _sp

for d in sorted(prof.donor_id.unique()):
    hits = mt_hits[str(d)]["elevated_clusters"]
    if not hits:
        continue
    m = (A.obs.donor_id == d) & A.obs.is_tumor
    sel = m & A.obs.leiden.isin(hits)
    rest = m & ~A.obs.leiden.isin(hits)
    per_gene = {}
    for g in mt_present:
        v = A[:, g].X
        v = np.asarray(v.todense()).ravel() if _sp.issparse(v) else np.asarray(v).ravel()
        a_, b_ = v[sel.to_numpy()], v[rest.to_numpy()]
        per_gene[g] = {
            "mean_elevated": float(a_.mean()),
            "mean_rest": float(b_.mean()),
            "log2fc": float(np.log2((a_.mean() + 0.01) / (b_.mean() + 0.01))),
            "mwu_p": float(stats.mannwhitneyu(a_, b_, alternative="greater").pvalue),
        }
    mt_hits[str(d)]["per_gene"] = per_gene
    nsig = sum(1 for v in per_gene.values() if v["mwu_p"] < 1e-10)
    mt_hits[str(d)]["n_genes_sig"] = nsig
    log(f"     逐基因：{nsig}/{len(mt_present)} 個 MT 基因 p < 1e-10")

# MT 高 cluster 的細胞 vs 其餘，檢定
for d in sorted(prof.donor_id.unique()):
    hits = mt_hits[str(d)]["elevated_clusters"]
    if not hits:
        continue
    m = tum.donor_id == d
    a = tum.loc[m & tum.leiden.isin(hits), "MT_score"]
    b = tum.loc[m & ~tum.leiden.isin(hits), "MT_score"]
    u = stats.mannwhitneyu(a, b, alternative="greater")
    out["metallothionein"][str(d)]["mwu_p_vs_rest"] = float(u.pvalue)
    out["metallothionein"][str(d)]["n_cells_elevated"] = int(len(a))

# ---------------------------- D. 定量化原文的定性結論（軟性對位 + 熵）
log("\n[D] 定量檢定：腫瘤細胞是否被凍結在單一成熟階段")
W = pd.read_csv(f"{RES}/15_match_weights.csv.gz", index_col=0)
obs15 = pd.read_csv(f"{RES}/15_mcl_obs.csv", index_col=0)
states = list(W.columns)
out["n_states_matched"] = len(states)
out["n_states_annotated_in_reference"] = 48
out["state_filter"] = "min_cells=10 per state per donor (15_residual_core.py)"

# 硬性指派（最大權重）
best = W.idxmax(axis=1)
common = obs15.index.intersection(W.index)
tum_idx = obs15.index[obs15.is_tumor.astype(bool)].intersection(W.index)
non_idx = obs15.index[~obs15.is_tumor.astype(bool)].intersection(W.index)

vc = best.loc[tum_idx].value_counts()
tab = pd.DataFrame({"n_cells": vc, "frac": vc / vc.sum()})
tab.index.name = "matched_normal_state"
tab.to_csv(f"{RES}/110_state_assignment.csv")
out["top_states"] = [
    {"state": str(s), "n": int(tab.loc[s, "n_cells"]), "frac": float(tab.loc[s, "frac"])} for s in tab.index[:8]
]
log("  前幾名對位狀態:")
for r in out["top_states"][:5]:
    log(f"    {r['state']:<34s} {r['n']:>6d}  {r['frac']:.1%}")

iv = [s for s in INTERVAL_ACTIVATED_TO_GCCOMMITTED if s in states]
mm = [s for s in MEMORY_STATES if s in states]
out["interval_states_used"] = iv
out["memory_states_used"] = mm
out["frac_hard_in_paper_interval"] = float(best.loc[tum_idx].isin(iv).mean())
out["frac_hard_in_memory"] = float(best.loc[tum_idx].isin(mm).mean())
# 軟性質量（不只看最大權重，而是看權重總和）
out["soft_mass_paper_interval"] = float(W.loc[tum_idx, iv].sum(axis=1).mean())
out["soft_mass_memory"] = float(W.loc[tum_idx, mm].sum(axis=1).mean())
log(f"  硬性指派落在原文區間（activated NBC → GC-committed）: {out['frac_hard_in_paper_interval']:.1%}")
log(f"  硬性指派落在記憶 B 區間:                              {out['frac_hard_in_memory']:.1%}")
log(f"  軟性權重質量 原文區間 / 記憶區間: {out['soft_mass_paper_interval']:.3f} / {out['soft_mass_memory']:.3f}")

# 每顆細胞的 Shannon 熵（以狀態數取對數歸一化）：測「是否凍結在單一狀態」
Wn = W.to_numpy()
Wn = np.clip(Wn, 1e-12, None)
H = -(Wn * np.log(Wn)).sum(axis=1) / np.log(Wn.shape[1])
H = pd.Series(H, index=W.index)
eff = pd.Series(np.exp(-(Wn * np.log(Wn)).sum(axis=1)), index=W.index)  # 有效狀態數
out["entropy"] = {
    # 15_residual_core.py 的 soft_match 只保留相關係數最高的 5 個狀態再做 softmax
    # （temperature = 0.02），所以「有效狀態數」的上限是 5，不是 46。腫瘤與正常細胞
    # 走的是同一套流程、同一個 k，因此兩者可以直接相比。
    "softmax_topk": 5,
    "softmax_temperature": 0.02,
    "max_possible_effective_states": 5,
    "tumor_median": float(H.loc[tum_idx].median()),
    "nontumor_median": float(H.loc[non_idx].median()) if len(non_idx) > 10 else None,
    "tumor_effective_states_median": float(eff.loc[tum_idx].median()),
    "nontumor_effective_states_median": float(eff.loc[non_idx].median()) if len(non_idx) > 10 else None,
    "n_nontumor_used": int(len(non_idx)),
}
if len(non_idx) > 10:
    u = stats.mannwhitneyu(H.loc[tum_idx], H.loc[non_idx], alternative="two-sided")
    out["entropy"]["mwu_p"] = float(u.pvalue)
    out["entropy"]["direction"] = (
        "tumour_more_spread" if H.loc[tum_idx].median() > H.loc[non_idx].median() else "tumour_more_focused"
    )
log("  熵 / 有效狀態數:", out["entropy"])

# 每位病人分開（原文只在 Fig 7 展示 M102）
out["per_donor_state"] = {}
for d in sorted(obs15.donor_id.unique()):
    idx = obs15.index[(obs15.donor_id == d) & obs15.is_tumor.astype(bool)].intersection(W.index)
    b = best.loc[idx].value_counts()
    out["per_donor_state"][str(d)] = {
        "n": int(len(idx)),
        "top_state": str(b.index[0]),
        "top_frac": float(b.iloc[0] / b.sum()),
        "frac_paper_interval": float(best.loc[idx].isin(iv).mean()),
        "entropy_median": float(H.loc[idx].median()),
        "effective_states_median": float(eff.loc[idx].median()),
        "n_states_with_1pct": int((b / b.sum() >= 0.01).sum()),
    }
    log(f"  {d}: 主狀態 {b.index[0]} ({b.iloc[0]/b.sum():.1%}), "
        f"≥1% 的狀態數 {out['per_donor_state'][str(d)]['n_states_with_1pct']}")

# 治療窗口的那個數字（SLC7A11），從既有殘差表取出以便 Section 2 引用
try:
    gr = pd.read_csv(f"{RES}/15_gene_residual_tumor_specific.csv", index_col=0)
    for g in ["SLC7A11", "MAFG"]:
        if g in gr.index:
            out.setdefault("residual_examples", {})[g] = {
                c: float(gr.loc[g, c]) for c in gr.columns if gr[c].dtype.kind == "f"
            }
except Exception as e:
    log("  (殘差表讀取跳過)", e)

with open(f"{RES}/110_revisit_fig7.json", "w") as f:
    json.dump(out, f, indent=2, ensure_ascii=False)
log("\n寫出 results/110_revisit_fig7.json")


# ------------------------------------------------ 圖用的直方圖（給 111 畫圖）
log("\n[E] 產生圖用的分布摘要")
fig = {}


def hist(v, lo, hi, nb=46):
    c, e = np.histogram(np.asarray(v, float), bins=nb, range=(lo, hi))
    return {"counts": [int(x) for x in c], "lo": float(lo), "hi": float(hi), "nbins": int(nb)}


# a. chrY 分數（腫瘤細胞，逐病人），含兩成分分界
fig["chrY"] = {}
for d in sorted(A.obs.donor_id.unique()):
    v = A.obs.loc[((A.obs.donor_id == d) & A.obs.is_tumor).to_numpy(), "chrY_score"]
    fig["chrY"][str(d)] = hist(v, -1.0, 1.0, 48)
    fig["chrY"][str(d)]["mean_loss"] = out["chrY"][str(d)]["mean_chrY_loss"]
    fig["chrY"][str(d)]["mean_ret"] = out["chrY"][str(d)]["mean_chrY_retained"]
    fig["chrY"][str(d)]["frac_loss"] = out["chrY"][str(d)]["frac_chrY_loss"]

# b. inferCNV 分數：扁桃腺正常 B 參考 / 檢體內非腫瘤 B / 腫瘤
CLO, CHI = 0.006, 0.020
fig["cnv"] = {
    "ref": hist(ref_sc, CLO, CHI, 40),
    "nontumor": hist(A.obs.loc[(~A.obs.is_tumor).to_numpy(), "cnv_score"], CLO, CHI, 40),
    "tumor": hist(A.obs.loc[A.obs.is_tumor.to_numpy(), "cnv_score"], CLO, CHI, 40),
    "p99_ref": p99,
}

# c. metallothionein：逐 cluster（已在 out 裡），加上顯著者標記
fig["mt"] = {
    str(d): {
        "clusters": out["metallothionein"][str(d)]["all_cluster_scores"],
        "n": {str(r.leiden): int(r.n) for r in prof[prof.donor_id == d].itertuples()},
        "elevated": out["metallothionein"][str(d)]["elevated_clusters"],
    }
    for d in sorted(prof.donor_id.unique())
}

# d. 軟性對位：狀態排名 + 有效狀態數分布
fig["states"] = out["top_states"]
fig["interval_states"] = iv
fig["eff"] = {
    "tumor": hist(eff.loc[tum_idx], 1.0, 5.0, 40),
    "nontumor": hist(eff.loc[non_idx], 1.0, 5.0, 40),
    "tumor_median": out["entropy"]["tumor_effective_states_median"],
    "nontumor_median": out["entropy"]["nontumor_effective_states_median"],
}
out["fig"] = fig
with open(f"{RES}/110_revisit_fig7.json", "w") as f:
    json.dump(out, f, indent=2, ensure_ascii=False)
log("已把圖用摘要寫回 results/110_revisit_fig7.json")
