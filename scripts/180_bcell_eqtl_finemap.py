#!/usr/bin/env python
"""步驟 180（B2）：分選 B 細胞的 eQTL 精細定位 —— 五個新生元件裡有因果變異嗎？

為什麼這比步驟 133 的 GTEx 版本好三個等級
------------------------------------------
1. **細胞類型對了**：GTEx 只有全血（混合）與 EBV 轉化的 LCL。
   eQTL Catalogue 有**分選的記憶 B 細胞**（OneK1K n=726/849、IBDverse n=227–313…），
   那正是腫瘤對應的正常狀態（ncsMBC）。
2. **解決了連鎖不平衡**：用 SuSiE 的 credible set（精細定位後的因果候選），
   不是 LD 拖著一整串的顯著變異。落在元件裡的 credible-set 變異是可解讀的。
3. **多個獨立世代可以合併**，檢定力比單一 GTEx 組織高得多。

⚠️ 一個重要情報：eQTL Catalogue **已經登記了記憶 B 細胞的 caQTL**
   （TenK10K QTD100034，n=913；FG_multiome QTD00068，n=1108），
   但目前**尚未在 FTP 上釋出**（全部 404）。那是未來最強的驗證資料，
   等它釋出應該立刻重跑。本步驟先用表現層的 eQTL。

設計與步驟 133 相同（前景＝五個新生元件，背景＝同窗其他 peak，
並做寬度配對），所以兩者可以對照。

輸出 results/180_*
"""
import os
import io
import gzip
import json
import subprocess

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
RES = f"{BASE}/results"
PROC = f"{BASE}/data/processed"
AUX = f"{BASE}/data/aux/eqtlcat"
CACHE = f"{AUX}/cs"
FTP = "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie"
FTP_BETA = "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/r8_beta/susie"

GENES = {"SLC7A11": ("ENSG00000151012", "chr4", 138242349),
         "MAFG": ("ENSG00000197063", "chr17", 81927735),
         "SLC3A2": ("ENSG00000168003", "chr11", 62856004),
         "GCLM": ("ENSG00000023909", "chr1", 93886831),
         "AIFM2": ("ENSG00000042286", "chr10", 70087919)}
WINDOW = 250_000
PIP_CUT = 0.01          # credible set 內的成員，PIP 很低的也算（保守：擴大前景機會）

# 已釋出、與 B 細胞最相關的資料集（study, dataset, label, n, 類別）
DATASETS = [
    ("QTS000038", "QTD000607", "OneK1K B_memory", 726, "memory_B"),
    ("QTS000048", "QTD000885", "OneK1K_re Memory_B", 849, "memory_B"),
    ("QTS000048", "QTD000880", "OneK1K_re Age-assoc_B", 132, "memory_B"),
    ("QTS000045", "QTD000830", "IBDverse Mem B CD27lo", 227, "memory_B"),
    ("QTS000045", "QTD000831", "IBDverse Mem B PDE4D", 301, "memory_B"),
    ("QTS000045", "QTD000832", "IBDverse Mem B PTPRJ", 313, "memory_B"),
    ("QTS000047", "QTD000870", "Perez_re Memory_B", 147, "memory_B"),
    ("QTS000046", "QTD000847", "Randolph_re Mem_B_NI", 79, "memory_B"),
    ("QTS000038", "QTD000608", "OneK1K B_naive", 876, "naive_B"),
    ("QTS000048", "QTD000887", "OneK1K_re Naive_B", 889, "naive_B"),
    ("QTS000045", "QTD000834", "IBDverse Naive B", 307, "naive_B"),
    ("QTS000026", "QTD000474", "Schmiedel B-cell_naive", 91, "naive_B"),
    ("QTS000038", "QTD000606", "OneK1K B_intermediate", 618, "other_B"),
    ("QTS000037", "QTD000597", "Perez_2022 B", 163, "other_B"),
    ("QTS000013", "QTD000110", "GEUVADIS LCL", 445, "LCL"),
    ("QTS000044", "QTD000764", "MAGE LCL", 675, "LCL"),
    ("QTS000029", "QTD000539", "TwinsUK LCL", 418, "LCL"),
    ("QTS000015", "QTD000221", "GTEx_v10 LCL", 326, "LCL"),
]

ELEMENTS = {"SLC7A11": [("chr4", 138222944, 138223474, "-19.1 kb", True),
                        ("chr4", 138392601, 138393029, "+150.5 kb", True),
                        ("chr4", 138406192, 138406449, "+164.0 kb", True),
                        ("chr4", 138341498, 138342345, "+99.6 kb", False),
                        ("chr4", 138348035, 138348490, "+105.9 kb", False)]}


def log(*a):
    print("[180]", *a, flush=True)


def fetch(qts, qtd):
    """舊資料集在 susie/ 下是 .tsv.gz，新的在 r8_beta/susie/ 下是 .parquet。"""
    os.makedirs(CACHE, exist_ok=True)
    for ext in (".tsv.gz", ".parquet"):
        f0 = f"{CACHE}/{qtd}.credible_sets{ext}"
        if os.path.exists(f0) and os.path.getsize(f0) > 1000:
            return f0
    for root, ext in ((FTP, ".tsv.gz"), (FTP_BETA, ".parquet"),
                      (FTP_BETA, ".tsv.gz"), (FTP, ".parquet")):
        fn = f"{CACHE}/{qtd}.credible_sets{ext}"
        url = f"{root}/{qts}/{qtd}/{qtd}.credible_sets{ext}"
        p = subprocess.run(["curl", "-sfL", "--retry", "2", "-m", "300", "-o", fn, url],
                           capture_output=True)
        if p.returncode == 0 and os.path.exists(fn) and os.path.getsize(fn) > 1000:
            return fn
        if os.path.exists(fn):
            os.remove(fn)
    log("  下載失敗（四種路徑都試過）:", qtd)
    return None


def main():
    mare = pd.read_csv(f"{RES}/90_mare_sites.csv")
    mg = mare[mare.gene == "MAFG"].drop_duplicates("peak_id")
    ELEMENTS["MAFG"] = [(r.chrom, int(r.el_start), int(r.el_end),
                         "%+.1f kb" % (r.dist / 1000), True) for _, r in mg.iterrows()]
    peaks = pd.read_csv(f"{PROC}/atlas_peaks_sorted.csv")

    want = {v[0]: k for k, v in GENES.items()}
    rows = []
    log("下載 %d 個資料集的 credible sets（每個約 0.8 MB）…" % len(DATASETS))
    for qts, qtd, label, n, cat in DATASETS:
        fn = fetch(qts, qtd)
        if not fn:
            continue
        got = 0
        if fn.endswith(".parquet"):
            df = pd.read_parquet(fn)
        else:
            df = pd.read_csv(fn, sep="\t")
        df["gid"] = df["gene_id"].astype(str).str.split(".").str[0]
        df = df[df.gid.isin(want)]
        for _, f in df.iterrows():
            v = str(f["variant"])
            try:
                chrom, pos = v.split("_")[0], int(v.split("_")[1])
            except Exception:
                continue
            rows.append(dict(dataset=label, category=cat, n=n,
                             gene=want[f["gid"]], chrom=chrom, pos=pos,
                             variant=v, rsid=str(f.get("rsid", "")),
                             pip=float(f["pip"]), cs_size=int(f["cs_size"]),
                             pvalue=float(f["pvalue"]), beta=float(f["beta"])))
            got += 1
        log("  %-26s n=%-4d → %4d 個 credible-set 變異（我們的 5 個基因）" % (label, n, got))

    if not rows:
        log("沒有抓到任何 credible set"); return
    CS = pd.DataFrame(rows)
    CS.to_csv(f"{RES}/180_credible_sets.csv.gz", index=False, compression="gzip")
    log("\n總計 %d 筆 credible-set 變異，涵蓋 %d 個資料集"
        % (len(CS), CS.dataset.nunique()))
    log("\n每個基因有 credible set 的資料集數:")
    log("\n" + CS.groupby(["gene", "category"]).dataset.nunique().to_string())

    out = {"datasets": len(DATASETS), "n_cs_variants": int(len(CS)),
           # 2026-11-27 更新：這個註解已過期。TenK10K 的記憶 B caQTL 不在
           # eQTL Catalogue 上，作者自己發布在 Hugging Face
           # （anglixue/TenK10K_multiome，CC-BY-4.0）。六個元件全部有對應的
           # 受測 peak、q 全部 ≥ 0.983，但 chr4 的顯著基礎率只有 2.01%
           # （301/14,998），六個元件的期望命中數 = 0.125 —— 觀察到 0 是基礎率。
           # 詳見 results/206_external_validation_scouting.md。
           "caqtl_note": ("記憶 B 細胞 caQTL 已登記但尚未釋出："
                          "TenK10K QTD100034 n=913、FG_multiome QTD00068 n=1108"),
           "genes": {}}

    for gene, (gid, chrom, tss) in GENES.items():
        sub = CS[(CS.gene == gene) & (CS.chrom == chrom)]
        fg = ELEMENTS.get(gene, [])
        fg_iv = [(a, b) for (_c, a, b, _d, _k) in fg]
        lo, hi = tss - WINDOW, tss + WINDOW
        bg = peaks[(peaks.chrom == chrom) & (peaks.end > lo) & (peaks.start < hi)].copy()
        if fg_iv:
            bg = bg[~bg.apply(lambda r: any(r.start < b and r.end > a
                                            for a, b in fg_iv), axis=1)]
        bg["width"] = bg.end - bg.start
        log("\n=== %s | credible-set 變異 %d 個 | 前景 %d 元件 | 背景 %d peak ==="
            % (gene, len(sub), len(fg_iv), len(bg)))
        if not len(sub):
            out["genes"][gene] = {"n_cs_variants": 0}
            log("  這個基因在所有 B 細胞資料集裡都沒有 credible set → 無法檢定")
            continue

        uniq = sub.drop_duplicates("pos")
        def count(ivs, df):
            n = sum(int(((df.pos >= a) & (df.pos < b)).sum()) for a, b in ivs)
            bp = sum(b - a for a, b in ivs)
            return n, bp
        n_fg, bp_fg = count(fg_iv, uniq)
        bg_iv = list(zip(bg.start, bg.end))
        n_bg, bp_bg = count(bg_iv, uniq)
        wide = bg.nlargest(max(2 * max(len(fg_iv), 1), 5), "width")
        n_w, bp_w = count(list(zip(wide.start, wide.end)), uniq)

        def rate(n, bp):
            return 1000.0 * n / bp if bp else np.nan
        orr = pv = orw = pvw = np.nan
        if bp_fg and bp_bg:
            orr, pv = fisher_exact([[n_fg, max(bp_fg - n_fg, 0)],
                                    [n_bg, max(bp_bg - n_bg, 0)]], alternative="greater")
        if bp_fg and bp_w:
            orw, pvw = fisher_exact([[n_fg, max(bp_fg - n_fg, 0)],
                                     [n_w, max(bp_w - n_w, 0)]], alternative="greater")
        exp_fg = rate(n_bg, bp_bg) * bp_fg / 1000.0
        log("  唯一位置 %d 個；前景 %d bp → %d 個（期望 %.2f）；背景 %d bp → %d 個 (%.3f/kb)"
            % (len(uniq), bp_fg, n_fg, exp_fg, bp_bg, n_bg, rate(n_bg, bp_bg)))
        log("  Fisher OR %.2f, p = %.3g | 寬度配對 OR %.2f, p = %.3g"
            % (orr, pv, orw, pvw))

        hits = []
        for (c, a, b, lab, consistent) in fg:
            h = uniq[(uniq.pos >= a) & (uniq.pos < b)]
            for _, r in h.iterrows():
                allcs = sub[sub.pos == r.pos]
                hits.append(dict(element=lab, two_patient_consistent=consistent,
                                 chrom=c, pos=int(r.pos), rsid=r.rsid,
                                 max_pip=float(allcs.pip.max()),
                                 n_datasets=int(allcs.dataset.nunique()),
                                 datasets=sorted(allcs.dataset.unique())[:6],
                                 best_p=float(allcs.pvalue.min())))
        if hits:
            log("  ★ 落在新生元件裡的 credible-set 變異:")
            for h in hits:
                log("    %-10s %s:%d %s  max PIP %.3f  出現在 %d 個資料集  最佳 p %.2g"
                    % (h["element"], h["chrom"], h["pos"], h["rsid"],
                       h["max_pip"], h["n_datasets"], h["best_p"]))
        else:
            log("  元件內沒有 credible-set 變異")

        out["genes"][gene] = {
            "n_cs_variants": int(len(sub)), "n_unique_positions": int(len(uniq)),
            "foreground": {"bp": int(bp_fg), "n": int(n_fg),
                           "expected": float(exp_fg)},
            "background": {"bp": int(bp_bg), "n": int(n_bg),
                           "per_kb": rate(n_bg, bp_bg),
                           "odds_ratio": float(orr), "p": float(pv)},
            "width_matched": {"bp": int(bp_w), "n": int(n_w),
                              "odds_ratio": float(orw), "p": float(pvw)},
            "hits": hits,
            "datasets_with_cs": sorted(sub.dataset.unique()),
        }

    json.dump(out, open(f"{RES}/180_bcell_eqtl.json", "w"), indent=2, ensure_ascii=False)
    log("\nDONE → results/180_*")


if __name__ == "__main__":
    main()
