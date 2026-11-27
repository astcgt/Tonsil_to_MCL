#!/usr/bin/env python
"""步驟 230：xCT／鐵死亡模組在 GSE303064 的 19 位獨立 MCL 病人身上複製得出來嗎？

規則全部來自 `results/230_PREREGISTRATION.md`（在任何數字落地之前寫完並提交，
commit 8c17124）。本檔不新增任何統計量：細胞 QC、B 細胞門檻、惡性／嚴格非惡性的
切法、跨樣本檢定與多重校正，全部沿用 `scripts/41_external_genelevel.py`。
新增的只有兩件事：
  1. 基因盤換成鐵死亡模組（步驟 41 那 39 個基因裡一個都沒有）
  2. G0 可偵測性門檻 —— 低表現基因在 3' scRNA 裡 dropout 嚴重，
     測不到就判 `no_power`，不得寫成陰性結果

輸出 results/230_*
"""
import os
import sys
import glob
import json
import re

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _compat  # noqa

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.stats import wilcoxon

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
EXT = f"{BASE}/data/external_mcl"
RES = f"{BASE}/results"
sc.settings.n_jobs = 16

# --- 事先固定的基因盤（230_PREREGISTRATION.md）---
TARGET = ["SLC7A11"]
MODULE = ["MAFG", "SLC3A2", "GCLC", "GCLM", "GSS", "GPX4", "AIFM2", "SCLY"]
LEVER = ["NAMPT", "NAPRT"]
POS = ["CCND1", "SOX11"]
NEG = ["RPL13A", "CD3E", "LYZ"]
GENES = TARGET + MODULE + LEVER + POS + NEG

# --- 步驟 41 的既有門檻，不改 ---
BCELL = ["CD79A", "CD79B", "MS4A1", "CD19", "BANK1", "PAX5"]
TCELL = ["CD3D", "CD3E", "CD2", "TRAC", "IL7R", "CD247", "TRBC2"]
MYE = ["LYZ", "CST3", "AIF1", "FCN1", "S100A8", "S100A9", "ITGAM", "CD14"]

# --- G0 門檻 ---
G0_MIN_FRAC = 0.05      # 惡性細胞中偵測率
G0_MIN_SAMPLES = 10     # 至少幾個樣本達標
R3_MIN_POSITIVE = 15    # R3 要求的同向樣本數


def log(*a):
    print("[230]", *a, flush=True)


def main():
    rows, det = [], []
    for f in sorted(glob.glob(f"{EXT}/GSE303064/*.h5")):
        s = re.sub(r".*_(S\d+)\.filtered.*", r"\1", os.path.basename(f))
        A = sc.read_10x_h5(f)
        A.var_names_make_unique()
        A.var["mt"] = A.var_names.str.upper().str.startswith("MT-")
        sc.pp.calculate_qc_metrics(A, qc_vars=["mt"], inplace=True,
                                   percent_top=None, log1p=False)
        A = A[(A.obs.total_counts >= 800) & (A.obs.n_genes_by_counts >= 400)
              & (A.obs.pct_counts_mt < 20)].copy()
        if A.n_obs < 300:
            log("  %s：QC 後只剩 %d 顆 → 跳過" % (s, A.n_obs))
            continue
        A.layers["counts"] = A.X.copy()          # G0 需要原始 count
        sc.pp.normalize_total(A, target_sum=1e4)
        sc.pp.log1p(A)
        for nm, gl in [("B", BCELL), ("T", TCELL), ("My", MYE)]:
            gg = [x for x in gl if x in A.var_names]
            sc.tl.score_genes(A, gg, score_name=nm + "_s")
        B = A[(A.obs.B_s > 0.15) & (A.obs.T_s < 0.1) & (A.obs.My_s < 0.1)].copy()
        if B.n_obs < 150 or "CCND1" not in B.var_names:
            log("  %s：B 細胞 %d 顆 → 跳過" % (s, B.n_obs))
            continue
        X = B.X.toarray() if sp.issparse(B.X) else np.asarray(B.X)
        C = B.layers["counts"]
        C = C.toarray() if sp.issparse(C) else np.asarray(C)
        gi = {k: i for i, k in enumerate(B.var_names)}
        cc = X[:, gi["CCND1"]]
        mal = cc > 0.5
        nor = cc == 0
        if mal.sum() < 100 or nor.sum() < 50:
            log("  %s：惡性 %d / 嚴格非惡性 %d → 跳過"
                % (s, int(mal.sum()), int(nor.sum())))
            continue
        rec = dict(sample=s, n_mal=int(mal.sum()), n_nor=int(nor.sum()))
        drec = dict(sample=s)
        for g in GENES:
            if g in gi:
                rec[g] = float(X[mal, gi[g]].mean() - X[nor, gi[g]].mean())
                drec[g] = float((C[mal, gi[g]] > 0).mean())
            else:
                rec[g] = np.nan
                drec[g] = np.nan
        rows.append(rec)
        det.append(drec)
        log("  %s：惡性 %d、嚴格非惡性 %d | SLC7A11 偵測率 %.3f、delta %+.4f"
            % (s, rec["n_mal"], rec["n_nor"], drec["SLC7A11"], rec["SLC7A11"]))
        del A, B, X, C

    df = pd.DataFrame(rows)
    dd = pd.DataFrame(det)
    df.to_csv(f"{RES}/230_external_module_delta.csv", index=False)
    dd.to_csv(f"{RES}/230_external_detection.csv", index=False)
    n = len(df)
    log("\n可用樣本 %d 個" % n)
    if n < 5:
        json.dump({"verdict": "VOID_too_few_samples", "n_samples": n},
                  open(f"{RES}/230_external_module.json", "w"),
                  indent=2, ensure_ascii=False)
        return 1

    # ------------------------------------------------- 跨樣本檢定 ---
    res = []
    for g in GENES:
        if g not in df:
            continue
        v = df[g].dropna().to_numpy()
        if len(v) < 5:
            continue
        try:
            p = float(wilcoxon(v).pvalue)
        except Exception:
            p = float("nan")
        res.append(dict(gene=g, median_delta=float(np.median(v)),
                        n_pos=int((v > 0).sum()), n_neg=int((v < 0).sum()),
                        n=len(v), pvalue=p,
                        det_median=float(np.median(dd[g].dropna()))
                        if g in dd else float("nan")))
    r = pd.DataFrame(res)
    r["padj"] = np.minimum(1.0, r.pvalue * len(r))     # 步驟 41 的保守校正
    r.to_csv(f"{RES}/230_external_gene_stats.csv", index=False)

    def row(g):
        m = r[r.gene == g]
        return None if m.empty else m.iloc[0]

    out = {"cohort": "GSE303064", "n_samples": n,
           "n_genes_tested": int(len(r)),
           "prereg": "results/230_PREREGISTRATION.md",
           "stats": r.to_dict(orient="records")}

    # ------------------------------------------------- G0 守門 ---
    g0_ok_samples = int((dd["SLC7A11"] >= G0_MIN_FRAC).sum()) if "SLC7A11" in dd else 0
    g0 = bool(g0_ok_samples >= G0_MIN_SAMPLES)
    out["G0_detectability"] = {
        "rule": "SLC7A11 原始 count > 0 的惡性細胞比例 ≥ %.2f，且達標樣本數 ≥ %d"
                % (G0_MIN_FRAC, G0_MIN_SAMPLES),
        "n_samples_passing": g0_ok_samples,
        "per_sample_detection": dd.set_index("sample")["SLC7A11"].round(4).to_dict()
                                 if "SLC7A11" in dd else {},
        "pass": g0}

    log("\n=== G0 可偵測性（事先門檻：≥ %.0f%% 惡性細胞、≥ %d 個樣本）==="
        % (G0_MIN_FRAC * 100, G0_MIN_SAMPLES))
    log("  SLC7A11 達標樣本 %d / %d | 偵測率中位數 %.4f"
        % (g0_ok_samples, n, float(np.median(dd["SLC7A11"].dropna()))))
    log("  G0 %s" % ("通過" if g0 else "未通過 → 本世代判 no_power，不下陰性結論"))

    # ------------------------------------------------- R1 / R2 守門 ---
    r1 = all((row(g) is not None and row(g)["median_delta"] > 0
              and row(g)["padj"] < 0.05) for g in POS)
    r2_bad = [g for g in NEG if row(g) is not None
              and row(g)["median_delta"] > 0 and row(g)["padj"] < 0.05]
    out["R1_positive_controls"] = {
        "rule": "CCND1 與 SOX11 的 median delta > 0 且 padj < 0.05",
        "values": {g: (None if row(g) is None else
                       {"median_delta": float(row(g)["median_delta"]),
                        "n_pos": int(row(g)["n_pos"]), "n": int(row(g)["n"]),
                        "padj": float(row(g)["padj"])}) for g in POS},
        "pass": bool(r1)}
    out["R2_negative_controls"] = {
        "rule": "RPL13A／CD3E／LYZ 不得有 padj < 0.05 的正向差異",
        "violations": r2_bad,
        "values": {g: (None if row(g) is None else
                       {"median_delta": float(row(g)["median_delta"]),
                        "padj": float(row(g)["padj"])}) for g in NEG},
        "pass": bool(not r2_bad)}

    log("\n=== R1 陽性對照（守門）===")
    for g in POS:
        q = row(g)
        log("  %-8s median delta %+8.4f | %2d/%2d 同向 | padj %.3g"
            % (g, q["median_delta"], q["n_pos"], q["n"], q["padj"]) if q is not None
            else "  %-8s 不在資料裡" % g)
    log("  R1 %s" % ("通過" if r1 else "未通過 → 整支作廢"))
    log("\n=== R2 陰性對照（守門）===")
    for g in NEG:
        q = row(g)
        log("  %-8s median delta %+8.4f | padj %.3g"
            % (g, q["median_delta"], q["padj"]) if q is not None
            else "  %-8s 不在資料裡" % g)
    log("  R2 %s" % ("通過" if not r2_bad else "未通過（%s）→ 整支作廢" % r2_bad))

    # ------------------------------------------------- 全盤 ---
    log("\n=== 全基因盤：惡性 B − 同檢體內嚴格非惡性 B（log1p CP10K）===")
    log("%-9s %12s %9s %10s %12s" % ("gene", "median_delta", "同向", "padj", "偵測率中位"))
    for _, q in r.sort_values("median_delta", ascending=False).iterrows():
        tag = ("★標的" if q.gene in TARGET else
               "模組" if q.gene in MODULE else
               "槓桿" if q.gene in LEVER else
               "陽對" if q.gene in POS else "陰對")
        log("%-9s %+12.4f %5d/%-3d %10.3g %12.4f  %s"
            % (q.gene, q.median_delta, q.n_pos, q.n, q.padj, q.det_median, tag))

    # ------------------------------------------------- R3 主要判定 ---
    q = row("SLC7A11")
    r3 = bool(q is not None and q["median_delta"] > 0
              and q["n_pos"] >= R3_MIN_POSITIVE and q["padj"] < 0.05)
    out["R3_primary"] = {
        "rule": "SLC7A11 median delta > 0、≥ %d/%d 樣本同向、padj < 0.05"
                % (R3_MIN_POSITIVE, n),
        "values": None if q is None else
                  {"median_delta": float(q["median_delta"]),
                   "n_pos": int(q["n_pos"]), "n": int(q["n"]),
                   "pvalue": float(q["pvalue"]), "padj": float(q["padj"])},
        "pass": r3}

    # ------------------------------------------------- R4 模組 vs 單基因 ---
    core = ["GCLC", "GCLM", "SLC3A2"]
    core_sig = [g for g in core if row(g) is not None and row(g)["padj"] < 0.05]
    out["R4_module_vs_single_gene"] = {
        "rule": "R3 成立但 GCLC／GCLM／SLC3A2 全不顯著 → 模組說法降級為單一轉運蛋白",
        "core_significant": core_sig,
        "downgrade_required": bool(r3 and not core_sig)}

    # ------------------------------------------------- 判定 ---
    if not (r1 and not r2_bad):
        verdict = "VOID_controls_failed"
    elif not g0:
        verdict = "no_power"
    elif r3:
        verdict = "replicates"
    else:
        verdict = "does_not_replicate"
    out["verdict"] = verdict

    log("\n=== R3 主要判定 ===")
    if q is not None:
        log("  SLC7A11 median delta %+.4f | %d/%d 同向（需 ≥ %d）| padj %.3g"
            % (q["median_delta"], q["n_pos"], q["n"], R3_MIN_POSITIVE, q["padj"]))
    log("  R3 %s" % ("通過" if r3 else "未通過"))
    log("\n*** 判定：%s ***" % verdict)
    if verdict == "no_power":
        log("    這個世代答不了這個問題。不得寫成陰性結果。")
    elif verdict == "does_not_replicate":
        log("    第二個指向反面的獨立世代，而且亞型是對的 ——")
        log("    依預先註冊的 R5，不得再用「亞型不同」解釋掉。")
    elif verdict == "replicates":
        log("    19 位獨立病人、正確亞型上複製成立。")
        if out["R4_module_vs_single_gene"]["downgrade_required"]:
            log("    但 R4 要求把「模組」降級為「單一轉運蛋白上調」。")

    json.dump(out, open(f"{RES}/230_external_module.json", "w"),
              indent=2, ensure_ascii=False)
    log("\n→ results/230_external_module.json、230_external_gene_stats.csv、"
        "230_external_module_delta.csv、230_external_detection.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
