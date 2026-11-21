#!/usr/bin/env python
"""
步驟 99：接觸圖結果的第三層混淆排除 —— 把 peak 的活性強度也扣掉。

要排除的另一種解釋
------------------
「新生元件比未重接線 peak 更常接觸啟動子」這個結果，有一個不涉及生物學的解釋：
**強 peak 既容易被判為相關、也容易被預測成高接觸**。
換句話說 peak 的活性本身可能同時驅動自變數與因變數。

步驟 98 已經把寬度（活性的代理）與距離迴歸掉。這一步直接用**實測可及性**：
從扁桃體圖譜的 46,216 顆正常 B 細胞算每個 peak 的偵測率與平均訊號，
當成第三個共變數再迴歸一次。

若效應在扣掉「寬度 + 距離 + 實測可及性」之後仍然存在，
那麼剩下的解釋就只有「這些元件真的位在接觸優勢位置」。
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _compat  # noqa
import numpy as np
import pandas as pd
import anndata as ad

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
RES = f"{BASE}/results"


def log(*a):
    print("[99]", *a, flush=True)


df = pd.read_csv(f"{RES}/98_contact_multi.csv")
df["pid"] = df.label.str.rsplit("_", n=1).str[-1].astype(int)
pids = np.sort(df.pid.unique())
log(f"需要可及性的 peak: {len(pids)} 個")

log("載入扁桃體正常 B 細胞 ATAC ...")
A = ad.read_h5ad(f"{BASE}/data/processed/tonsil_multiome_atac_B.h5ad")
X = A.X.tocsc() if hasattr(A.X, "tocsc") else A.X
n = A.shape[0]
det = np.asarray((X[:, pids] > 0).sum(axis=0)).ravel() / n
sig = np.asarray(X[:, pids].sum(axis=0)).ravel() / n
acc = pd.DataFrame(dict(pid=pids, acc_det=det, acc_mean=sig))
log(f"偵測率範圍 {det.min():.4f} – {det.max():.4f}（中位 {np.median(det):.4f}）")
acc.to_csv(f"{RES}/99_peak_accessibility.csv", index=False)
del A, X

df = df.merge(acc, on="pid", how="left")

# ---- 先確認這個混淆是否真的存在
from scipy.stats import mannwhitneyu, spearmanr, combine_pvalues, rankdata

dn = df[df.kind == "de_novo"].drop_duplicates("pid")
cl = df[df.kind == "control"].drop_duplicates("pid")
log("\n=== 混淆存在嗎？===")
for gene in df.gene.unique():
    a = dn[dn.gene == gene]; b = cl[cl.gene == gene]
    u = mannwhitneyu(a.acc_det, b.acc_det, alternative="greater")
    log(f"{gene:8s} 新生元件偵測率中位 {a.acc_det.median():.4f} (n={len(a)}) vs "
        f"對照 {b.acc_det.median():.4f} (n={len(b)})  p = {u.pvalue:.3g}")
    s = df[df.gene == gene]
    log(f"         O/E ~ 可及性 rho = {spearmanr(s.acc_det, s.oe).statistic:+.3f}"
        f" | O/E ~ 寬度 rho = {spearmanr(s.width, s.oe).statistic:+.3f}")


def adj(s, cols):
    y = rankdata(s.oe)
    X = [np.ones(len(s))] + [rankdata(np.log(s[c] + 1e-9)) for c in cols]
    X = np.column_stack(X)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return y - X @ beta


rows = []
for (gene, grp, cond), s in df.groupby(["gene", "group", "condition"]):
    if (s.kind == "de_novo").sum() < 3 or (s.kind == "control").sum() < 5:
        continue
    m = (s.kind == "de_novo").values
    out = dict(gene=gene, group=grp, condition=cond)
    for name, cols in [("p_wd", ["width", "bin_dist"]),
                       ("p_wda", ["width", "bin_dist", "acc_det"]),
                       ("p_a", ["acc_det"])]:
        r = adj(s, cols)
        out[name] = float(mannwhitneyu(r[m], r[~m], alternative="greater").pvalue)
        out[name.replace("p_", "d_")] = round(float(np.median(r[m]) - np.median(r[~m])), 2)
    rows.append(out)
R = pd.DataFrame(rows)
R.to_csv(f"{RES}/99_contact_stats_accadj.csv", index=False)

pd.set_option("display.width", 200)
log("\n" + "=" * 78)
log("=== 三種校正下的結果（p_wd=寬度+距離；p_wda=再加實測可及性；p_a=只校正可及性）===")
for (gene, grp), s in R.groupby(["gene", "group"]):
    st_wd = combine_pvalues(s.p_wd.clip(1e-300, 1), method="stouffer").pvalue
    st_wda = combine_pvalues(s.p_wda.clip(1e-300, 1), method="stouffer").pvalue
    log(f"{gene:8s} {grp:18s} n={len(s):2d} | 顯著 {int((s.p_wd<0.05).sum())}/{len(s)} → "
        f"{int((s.p_wda<0.05).sum())}/{len(s)} (加可及性) | "
        f"Stouffer {st_wd:.3g} → {st_wda:.3g} | 殘差差 {s.d_wd.median():+.1f} → {s.d_wda.median():+.1f}")
log("\nDONE")
