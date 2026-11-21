#!/usr/bin/env python
"""
步驟 51：對所有候選基因做文獻新穎性篩選。

使用者要求「找到的 pathway 和藥物是目前沒有其他文獻指出的」，
所以先用 PubMed 計數把已被充分研究的濾掉。

對每個基因查三個層次：
  n_mcl        : <gene> AND (mantle cell lymphoma)
  n_lymphoma   : <gene> AND (lymphoma OR leukemia)
  n_total      : <gene> 全部（當作「這個基因有多紅」的基準）
新穎性定義：n_mcl 很小（<=3）但在我們的資料中訊號很強。
"""
import os, sys, time, json, urllib.parse, urllib.request
import pandas as pd, numpy as np

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
RES = f"{BASE}/results"
EUT = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

def log(*a): print("[51]", *a, flush=True)

def count(term, retry=3):
    q = urllib.parse.urlencode({"db":"pubmed","term":term,"retmode":"json","rettype":"count"})
    for _ in range(retry):
        try:
            with urllib.request.urlopen(f"{EUT}?{q}", timeout=30) as r:
                d = json.load(r)
            return int(d["esearchresult"]["count"])
        except Exception:
            time.sleep(1.5)
    return np.nan

# ---- 候選清單：來自本研究各層分析 ----
cands = {}
ex = pd.read_csv(f"{RES}/15_gene_residual_tumor_specific.csv")
for g in ex.sort_values("resid_tumor_specific", ascending=False).head(25).gene:
    cands.setdefault(g, []).append("表現殘差上升 Top25")
for g in ex.sort_values("resid_tumor_specific").head(20).gene:
    cands.setdefault(g, []).append("表現殘差下降 Top20")
tf = pd.read_csv(f"{RES}/19_tf_ranking_paired.csv", index_col=0)
s = tf[tf.same_dir]
for g in s.sort_values("paired_mean", ascending=False).head(20).index:
    cands.setdefault(g, []).append("TF 活性上升 Top20")
for g in s.sort_values("paired_mean").head(15).index:
    cands.setdefault(g, []).append("TF 活性下降 Top15")
dep = pd.read_csv(f"{RES}/09_depmap_mcl_dependency.csv")
for g in dep.sort_values("z_vs_all").head(25).gene:
    cands.setdefault(g, []).append("DepMap MCL 選擇性依賴 Top25")

log(f"候選基因 {len(cands)} 個，開始查 PubMed ...")
rows = []
for i, (g, src) in enumerate(sorted(cands.items())):
    n_mcl = count(f'"{g}"[All Fields] AND ("mantle cell lymphoma"[All Fields])')
    time.sleep(0.34)
    n_lym = count(f'"{g}"[All Fields] AND (lymphoma[All Fields] OR leukemia[All Fields])')
    time.sleep(0.34)
    n_all = count(f'"{g}"[All Fields]')
    time.sleep(0.34)
    rows.append(dict(gene=g, n_mcl=n_mcl, n_lymphoma=n_lym, n_total=n_all,
                     sources=" ; ".join(src)))
    if (i+1) % 20 == 0: log(f"  已查 {i+1}/{len(cands)}")
df = pd.DataFrame(rows)

# 加回量化證據
df = df.merge(ex[["gene","resid_tumor_specific","mean_observed","mean_expected"]], on="gene", how="left")
df = df.merge(dep[["gene","mcl_mean","z_vs_all"]].rename(
    columns={"mcl_mean":"depmap_effect","z_vs_all":"depmap_z"}), on="gene", how="left")
tfp = tf.reset_index().rename(columns={"index":"gene"})[["gene","paired_mean","same_dir"]]
df = df.merge(tfp.rename(columns={"paired_mean":"tf_paired","same_dir":"tf_consistent"}), on="gene", how="left")
df.to_csv(f"{RES}/51_novelty_screen.csv", index=False)
pd.set_option("display.width", 250)

log("\n=== 新穎候選：MCL 文獻 <= 3 篇，但本研究訊號強 ===")
nov = df[(df.n_mcl <= 3)].copy()
nov["evidence"] = (nov.resid_tumor_specific.abs().fillna(0)
                   + nov.depmap_z.abs().fillna(0)/3
                   + nov.tf_paired.abs().fillna(0))
nov = nov.sort_values("evidence", ascending=False)
log("\n" + nov.head(30)[["gene","n_mcl","n_lymphoma","n_total","resid_tumor_specific",
                          "depmap_z","tf_paired","sources"]].round(3).to_string(index=False))
log("\n=== 對照：文獻最多的（已充分研究，排除）===")
log("\n" + df.sort_values("n_mcl", ascending=False).head(15)[["gene","n_mcl","n_lymphoma","n_total"]].to_string(index=False))
log("DONE")
