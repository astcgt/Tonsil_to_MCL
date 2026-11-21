#!/usr/bin/env python
"""
步驟 92：AlphaGenome 2 的 cell_conditioning —— 用我們自己的細胞定義條件。

為什麼這個設計比步驟 91 好
--------------------------
步驟 91(A) 失敗在 peak 寬度混淆：我拿「新生元件」和「其他 peak」比，
兩組的區間長度不同，所以比較本身就不公平。

這一步把比較改成**同一段 DNA、只換細胞狀態**：
  條件 A = 腫瘤細胞（同一位病人）
  條件 B = 同一份檢體裡的嚴格非腫瘤 B 細胞
兩邊讀的是**完全相同的座標**，所以寬度、GC、距離全部自動消掉。
而且兩組細胞來自同一次 10x 反應，donor 與批次也配對。

AlphaGenome 的文件明確警告 encoder 不校正批次效應，
所以我們**按病人分開跑**，絕不跨病人比較。

M413：790 顆嚴格非腫瘤 B 細胞 → 主要分析（各取 256 顆，達 API 上限）
M102：只有 29 顆 → 次要分析（統計力不足，照實報告）

金鑰只從 .secrets/alphagenome_key 讀取，不寫入任何輸出。
"""
import os, sys, re, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _compat  # noqa
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
RES = f"{BASE}/results"
with open(f"{BASE}/.secrets/alphagenome_key") as fh:
    _KEY = fh.read().strip()

from alphagenome.data import genome
from alphagenome.models.v2 import dna_client, dna_output

SEQLEN = 524288
MAXC = dna_client.MAX_CONDITIONING_CELLS     # 256
RNG = np.random.default_rng(0)
TSS = {"SLC7A11": ("chr4", 138242349), "MAFG": ("chr17", 81927735)}


def log(*a):
    print("[92]", *a, flush=True)


# ------------------------------------------------- 基因符號 → Ensembl ID
def symbol2ensg():
    import gzip
    mp = {}
    with gzip.open(f"{BASE}/data/genome/gencode.v44.annotation.gtf.gz", "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if f[2] != "gene":
                continue
            g = re.search(r'gene_id "([^".]+)', f[8])
            n = re.search(r'gene_name "([^"]+)"', f[8])
            if g and n and n.group(1) not in mp:
                mp[n.group(1)] = g.group(1)
    return mp


def build(A, idx, mp):
    """取出指定細胞，轉成 AlphaGenome 要的格式：原始整數 UMI + Ensembl ID。"""
    S = A[idx].copy()
    keep = [i for i, s in enumerate(S.var_names) if s in mp]
    S = S[:, keep].copy()
    S.var["ensg"] = [mp[s] for s in S.var_names]
    S = S[:, ~S.var.ensg.duplicated().values].copy()
    S.var_names = S.var.ensg.values
    S.var_names_make_unique()
    X = S.X.toarray() if sp.issparse(S.X) else np.asarray(S.X)
    S.X = sp.csr_matrix(np.rint(X).astype(np.float32))      # 確保整數
    S.obs = S.obs[["donor_id", "is_tumor"]].copy()
    return S


model = dna_client.create(_KEY)
log("連線建立（金鑰未輸出）")

log("建立 基因符號 → Ensembl 對照表 ...")
MP = symbol2ensg()
log(f"  {len(MP):,} 個基因")

lab = ad.read_h5ad(f"{BASE}/data/processed/mcl_tumor_labeled.h5ad")
raw = ad.read_h5ad(f"{BASE}/data/processed/mcl_gex_raw.h5ad")
raw = raw[raw.obs_names.isin(lab.obs_names)].copy()
raw.obs = lab.obs.loc[raw.obs_names, ["donor_id", "is_tumor", "is_tumor_cnv"]].copy()
log(f"原始計數矩陣對齊後: {raw.shape}")

# 元件座標
rw = pd.read_csv(f"{RES}/52_rewiring_annotated.csv.gz")
MARE = pd.read_csv(f"{RES}/90_mare_sites.csv")

rows, diag = [], []
for donor in ["M413", "M102"]:
    d = raw[raw.obs.donor_id == donor]
    tum = np.where(d.obs.is_tumor.astype(bool).values)[0]
    nor = np.where((~d.obs.is_tumor.astype(bool).values)
                   & (~d.obs.is_tumor_cnv.astype(bool).values))[0]
    n = min(len(tum), len(nor), MAXC)
    if n < 10:
        log(f"{donor}: 配對細胞不足（腫瘤 {len(tum)} / 正常 {len(nor)}），跳過"); continue
    tsel = RNG.choice(tum, n, replace=False)
    nsel = RNG.choice(nor, n, replace=False)
    log(f"\n=== {donor}: 每組 {n} 顆（腫瘤池 {len(tum)}、正常池 {len(nor)}）===")
    Atum = build(d, tsel, MP)
    Anor = build(d, nsel, MP)
    log(f"  條件矩陣: {Atum.shape} / {Anor.shape}")

    # --- 文件建議的診斷：兩種狀態在 embedding 空間分得開嗎？
    try:
        Etum = model.encode_cells(Atum, progress_bar=False)
        Enor = model.encode_cells(Anor, progress_bar=False)
        ct, cn = Etum.mean(0), Enor.mean(0)
        wt = np.linalg.norm(Etum - ct, axis=1).mean()
        wn = np.linalg.norm(Enor - cn, axis=1).mean()
        between = float(np.linalg.norm(ct - cn))
        sep = between / max((wt + wn) / 2, 1e-9)
        # 分離度比值會被組內散度稀釋，所以再做兩個嚴謹檢定
        E = np.vstack([Etum, Enor]); y = np.r_[np.ones(len(Etum)), np.zeros(len(Enor))]
        # ① 置換檢定：組間中心距離是否大於隨機分組
        null = []
        for _ in range(2000):
            yp = RNG.permutation(y)
            null.append(np.linalg.norm(E[yp == 1].mean(0) - E[yp == 0].mean(0)))
        p_perm = float((np.sum(np.array(null) >= between) + 1) / (len(null) + 1))
        # ② 交叉驗證 AUC：單顆細胞能不能被分類
        auc = np.nan
        try:
            from sklearn.linear_model import LogisticRegression
            from sklearn.model_selection import cross_val_score, StratifiedKFold
            auc = float(np.mean(cross_val_score(
                LogisticRegression(max_iter=2000, C=1.0), E, y,
                cv=StratifiedKFold(5, shuffle=True, random_state=0), scoring="roc_auc")))
        except Exception as e:
            log(f"    AUC 計算略過: {type(e).__name__}")
        log(f"  embedding 診斷: 組間中心距離 {between:.3f} / 組內平均散度 "
            f"{(wt+wn)/2:.3f} → 比值 {sep:.2f}")
        log(f"    置換檢定 p = {p_perm:.4g}（2000 次）｜ 單細胞分類 5-fold AUC = {auc:.3f}")
        diag.append(dict(donor=donor, n=n, latent_dim=int(Etum.shape[1]),
                         between=round(between, 4), within=round(float((wt + wn) / 2), 4),
                         separation=round(sep, 3), p_perm=p_perm,
                         cv_auc=round(auc, 4) if auc == auc else np.nan))
    except Exception as e:
        log(f"  embedding 診斷失敗: {type(e).__name__}: {e}")

    # --- 同一段 DNA、只換細胞狀態
    for gene, (chrom, tss) in TSS.items():
        iv = genome.Interval(chrom, tss - SEQLEN // 2, tss + SEQLEN // 2)
        preds = {}
        for nm, cond in [("tumour", Atum), ("normal", Anor)]:
            try:
                # ⚠️ API 限制：condition 與 cell_conditioning 不能同時給。
                # 用 cell_conditioning 時，細胞狀態完全由我們自己的細胞定義。
                preds[nm] = model.predict_interval(
                    interval=iv,
                    requested_outputs=[dna_output.OutputType.SC_ATAC,
                                       dna_output.OutputType.SC_GENE_EXPRESSION_DENOISED],
                    cell_conditioning=cond, cell_condition_name=nm)
            except Exception as e:
                log(f"    {gene}/{nm} 失敗: {type(e).__name__}: {str(e)[:150]}")
        if len(preds) < 2:
            continue

        els = rw[(rw.gene == gene) & rw.gained]
        els = els[(els.start > iv.start) & (els.end < iv.end)]
        mare = MARE[MARE.gene == gene]
        # ⭐ 決定性對照：同窗格內「未重接線」的 peak。
        # 因為讀出是同一座標的腫瘤/正常比值，寬度與 GC 在比值中完全抵銷，
        # 所以這組對照可以直接回答「效應是否專一於新生元件」。
        ctl = rw[(rw.gene == gene) & (~rw.gained) & (~rw.lost)
                 & (rw.start > iv.start) & (rw.end < iv.end)]

        # --- 基因層級讀出：去噪單細胞表現（每個基因一個值，不是分箱軌道）
        g0 = getattr(preds["normal"], "sc_gene_expression_denoised", None)
        g1 = getattr(preds["tumour"], "sc_gene_expression_denoised", None)
        if g0 is not None and g1 is not None:
            gm = g0.gene_metadata
            col = "gene_name" if "gene_name" in gm.columns else "gene_id"
            hit = np.where(gm[col].astype(str).values == gene)[0]
            if len(hit):
                a = float(np.nanmean(g0.values[hit, :]))
                b = float(np.nanmean(g1.values[hit, :]))
                rows.append(dict(donor=donor, gene=gene, track="基因表現(去噪, CPM)",
                                 feature=f"【{gene} 本身的預測表現】", dist=0, dr=np.nan,
                                 normal=a, tumour=b,
                                 log2fc=float(np.log2((b + 1e-6) / (a + 1e-6)))))
                # 同窗格其他基因作對照，確認不是全域縮放
                other = np.setdiff1d(np.arange(len(gm)), hit)
                if len(other):
                    oa = np.nanmean(g0.values[other, :], axis=1)
                    ob = np.nanmean(g1.values[other, :], axis=1)
                    ol = np.log2((ob + 1e-6) / (oa + 1e-6))
                    rows.append(dict(donor=donor, gene=gene, track="基因表現(去噪, CPM)",
                                     feature=f"【同窗格其他 {len(other)} 個基因對照】",
                                     dist=0, dr=np.nan,
                                     normal=float(np.nanmedian(oa)), tumour=float(np.nanmedian(ob)),
                                     log2fc=float(np.nanmedian(ol))))
                log(f"    {gene} 預測表現: 正常 {a:.4f} → 腫瘤 {b:.4f} CPM "
                    f"(log2FC {np.log2((b+1e-6)/(a+1e-6)):+.3f})")

        # --- 位置層級讀出：單細胞 ATAC（分箱軌道）
        for track_attr in ["sc_atac"]:
            t0 = getattr(preds["normal"], track_attr, None)
            t1 = getattr(preds["tumour"], track_attr, None)
            if t0 is None or t1 is None:
                continue
            res = t0.resolution
            v0 = t0.values.mean(axis=1)
            v1 = t1.values.mean(axis=1)

            def seg(v, s, e):
                i0 = max(0, (s - iv.start) // res); i1 = min(len(v), (e - iv.start) // res + 1)
                return float(v[i0:i1].mean()) if i1 > i0 else np.nan

            for _, r in els.iterrows():
                mid = (int(r.start) + int(r.end)) // 2
                a = seg(v0, mid - 250, mid + 251)
                b = seg(v1, mid - 250, mid + 251)
                rows.append(dict(donor=donor, gene=gene, track=track_attr,
                                 feature=f"新生元件 {int(r.peak_id)}",
                                 dist=int(r.dist), dr=round(float(r.dr), 3),
                                 normal=a, tumour=b,
                                 log2fc=float(np.log2((b + 1e-6) / (a + 1e-6)))))
            for _, r in mare.iterrows():
                a = seg(v0, int(r.pos) - 250, int(r.pos) + 251)
                b = seg(v1, int(r.pos) - 250, int(r.pos) + 251)
                rows.append(dict(donor=donor, gene=gene, track=track_attr,
                                 feature=f"MARE {int(r.pos)} ({r.motif})",
                                 dist=int(r.dist), dr=round(float(r.dr), 3),
                                 normal=a, tumour=b,
                                 log2fc=float(np.log2((b + 1e-6) / (a + 1e-6)))))
            for _, r in ctl.iterrows():
                mid = (int(r.start) + int(r.end)) // 2
                a = seg(v0, mid - 250, mid + 251)
                b = seg(v1, mid - 250, mid + 251)
                rows.append(dict(donor=donor, gene=gene, track=track_attr,
                                 feature=f"〔未重接線對照 {int(r.peak_id)}〕",
                                 dist=int(r.dist), dr=round(float(r.dr), 3),
                                 normal=a, tumour=b,
                                 log2fc=float(np.log2((b + 1e-6) / (a + 1e-6)))))
            # 全窗格對照：隨機 200 個位置，確認差異不是全域偏移
            mids = RNG.integers(iv.start + 1000, iv.end - 1000, 200)
            gl = [np.log2((seg(v1, m - 250, m + 251) + 1e-6) / (seg(v0, m - 250, m + 251) + 1e-6))
                  for m in mids]
            rows.append(dict(donor=donor, gene=gene, track=track_attr,
                             feature="【全窗格隨機對照 n=200】", dist=0, dr=np.nan,
                             normal=float(np.median(v0)), tumour=float(np.median(v1)),
                             log2fc=float(np.median(gl))))
            log(f"    {gene} / {track_attr}: 元件中位 log2FC "
                f"{np.median([r['log2fc'] for r in rows if r['gene']==gene and r['track']==track_attr and '元件' in r['feature'] and r['donor']==donor]):+.3f}"
                f" | 全域對照 {np.median(gl):+.3f}")

df = pd.DataFrame(rows)
if not df.empty:
    df.to_csv(f"{RES}/92_ag_cellcond.csv", index=False)
pd.DataFrame(diag).to_csv(f"{RES}/92_ag_encode_diag.csv", index=False)
if df.empty:
    log("序列預測全部失敗，只有 embedding 診斷可用"); sys.exit(0)

log("\n" + "=" * 70)
log("=== 結果：同一段 DNA，只換細胞狀態（腫瘤 vs 同檢體配對正常）===")
from scipy.stats import wilcoxon
pd.set_option("display.width", 220)
from scipy.stats import mannwhitneyu
for (donor, gene, track), s in df.groupby(["donor", "gene", "track"]):
    el = s[s.feature.str.startswith("新生元件")]
    pk = s[s.feature.str.startswith("〔未重接線對照")]
    if len(el) >= 3 and len(pk) >= 3:
        u = mannwhitneyu(el.log2fc, pk.log2fc, alternative="greater")
        log(f"\n>>> {donor}/{gene}/{track} 新生元件 vs 同窗格未重接線 peak（決定性對照）")
        log(f"    新生 {el.log2fc.median():+.3f} (n={len(el)}) vs 未重接線 "
            f"{pk.log2fc.median():+.3f} (n={len(pk)})  Mann-Whitney p = {u.pvalue:.4g}")
    el = s[s.feature.str.contains("元件|MARE") & ~s.feature.str.startswith("〔")]
    bg = s[s.feature.str.contains("全窗格")]
    if el.empty:
        continue
    bgv = float(bg.log2fc.iloc[0]) if len(bg) else 0.0
    try:
        w = wilcoxon(el.log2fc - bgv).pvalue if len(el) >= 5 else np.nan
    except Exception:
        w = np.nan
    log(f"\n{donor} / {gene} / {track}")
    log(f"  元件 log2FC 中位 {el.log2fc.median():+.3f} (n={len(el)})"
        f" | 全窗格對照 {bgv:+.3f} | 扣除對照後 Wilcoxon p={w:.4g}")
    log("\n" + el[["feature", "dist", "dr", "normal", "tumour", "log2fc"]]
        .sort_values("log2fc", ascending=False).round(4).to_string(index=False))
log("\nDONE")
