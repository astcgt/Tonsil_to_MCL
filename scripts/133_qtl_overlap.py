#!/usr/bin/env python
"""步驟 133：人類族群遺傳學的外部證據 —— 元件裡有沒有 eQTL？

為什麼這條線值得做
------------------
整份報告的弱點是「只有兩位病人」。eQTL 不需要第三位病人：它問的是
**在幾百上千個健康人身上，這段序列裡的常見變異會不會改變下游基因的表現**。
如果會，這段序列就是真人身上的功能性調控元件 —— 這比任何序列模型的預測都硬。

設計（刻意抄 2.6 節 A 的教訓：前景與背景必須配對）
---------------------------------------------------
前景：SLC7A11 的 5 個新生元件 / MAFG 的 MARE 元件
背景：同一個 ±250 kb 視窗裡、圖譜叫出來的**其他所有 peak**
檢定：每 kb 的顯著 eQTL 變異密度，Fisher 精確檢定
      再加一個**寬度配對**版本（背景只取最寬的 N 個），因為寬度本身就是
      「這個調控區有多強」的代理，上次就是這裡出事的。

資料來源：GTEx v10 single-tissue eQTL（REST API，不需下載整份資料）
  Whole_Blood、Cells_EBV-transformed_lymphocytes、Spleen 三個與 B 細胞最相關的組織

⚠️ 已知限制（必須跟結果一起報）
  1. GTEx 的顯著 eQTL 變異是**連鎖不平衡下的一整群**，落在元件裡不等於它就是因果變異。
     要判因果需要 fine-mapping（SuSiE / DAP-G）或 coloc，本步驟只做富集，不宣稱因果。
  2. GTEx 沒有扁桃腺，也沒有 MCL。全血與 LCL 是最接近的替代。
  3. **caQTL（染色質可及性 QTL）目前沒有公開 API**，eQTL Catalogue 的 API 已下線。
     要做 caQTL 必須下載 tabix 索引檔，而本機沒有 tabix。見報告的「待辦」。

輸出 results/133_*
"""
import os
import json
import time
import subprocess
import itertools

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
RES = f"{BASE}/results"
PROC = f"{BASE}/data/processed"
CACHE = f"{BASE}/data/aux/gtex"
os.makedirs(CACHE, exist_ok=True)

TISSUES = ["Whole_Blood", "Cells_EBV-transformed_lymphocytes", "Spleen"]
GENES = {
    "SLC7A11": ("ENSG00000151012.13", "chr4", 138242349),
    "MAFG": ("ENSG00000197063.11", "chr17", 81927735),
    "SLC3A2": ("ENSG00000168003.18", "chr11", 62856004),
}
WINDOW = 250_000
PSIG = 1e-4          # GTEx API 回傳的就是該組織的顯著 eQTL；再加一道保守門檻

ELEMENTS = {
    "SLC7A11": [("chr4", 138222944, 138223474, "-19.1 kb"),
                ("chr4", 138392601, 138393029, "+150.5 kb"),
                ("chr4", 138406192, 138406449, "+164.0 kb"),
                ("chr4", 138341498, 138342345, "+99.6 kb"),
                ("chr4", 138348035, 138348490, "+105.9 kb")],
}


def log(*a):
    print("[133]", *a, flush=True)


def gtex(gencode, tissue):
    """GTEx v10 single-tissue eQTL，分頁抓完，存快取。"""
    fn = f"{CACHE}/{gencode}_{tissue}.json"
    if os.path.exists(fn):
        return json.load(open(fn))
    rows, page = [], 0
    while True:
        url = ("https://gtexportal.org/api/v2/association/singleTissueEqtl"
               f"?gencodeId={gencode}&tissueSiteDetailId={tissue}"
               f"&datasetId=gtex_v10&page={page}&itemsPerPage=1000")
        p = subprocess.run(["curl", "-s", "-m", "60", url],
                           capture_output=True, text=True)
        try:
            d = json.loads(p.stdout)
        except Exception:
            log("  解析失敗", tissue, "page", page, p.stdout[:200]); break
        rows += d.get("data", [])
        pi = d.get("paging_info", {})
        if page + 1 >= pi.get("numberOfPages", 0):
            break
        page += 1
        time.sleep(0.4)
    json.dump(rows, open(fn, "w"))
    return rows


def main():
    peaks = pd.read_csv(f"{PROC}/atlas_peaks_sorted.csv")
    mare = pd.read_csv(f"{RES}/90_mare_sites.csv")
    out = {"tissues": TISSUES, "window_bp": WINDOW, "p_cut": PSIG, "genes": {}}

    # MAFG 的前景元件：90_mare_sites 裡屬於 MAFG 的那些 peak
    mg = mare[mare.gene == "MAFG"].drop_duplicates("peak_id")
    ELEMENTS["MAFG"] = [(r.chrom, int(r.el_start), int(r.el_end),
                         "%+.1f kb" % (r.dist / 1000)) for _, r in mg.iterrows()]
    log("MAFG 前景元件:", len(ELEMENTS["MAFG"]))

    allrows = []
    for gene, (gid, chrom, tss) in GENES.items():
        lo, hi = tss - WINDOW, tss + WINDOW
        fg = ELEMENTS.get(gene, [])
        fg_iv = [(a, b) for (_c, a, b, _d) in fg]
        bg = peaks[(peaks.chrom == chrom) & (peaks.end > lo) & (peaks.start < hi)].copy()
        # 背景排除任何與前景重疊的 peak
        def overlaps_fg(r):
            return any(r.start < b and r.end > a for a, b in fg_iv)
        bg["is_fg"] = bg.apply(overlaps_fg, axis=1) if len(fg_iv) else False
        bgc = bg[~bg.is_fg].copy()
        bgc["width"] = bgc.end - bgc.start
        log("\n=== %s  %s:%d±%dkb | 前景 %d 個元件 | 背景 %d 個 peak ==="
            % (gene, chrom, tss, WINDOW // 1000, len(fg_iv), len(bgc)))

        ev = []
        for t in TISSUES:
            r = gtex(gid, t)
            log("  GTEx %-34s %5d 個顯著 eQTL 變異" % (t, len(r)))
            for x in r:
                if x.get("pValue", 1) <= PSIG:
                    ev.append((t, int(x["pos"]), x["snpId"], float(x["pValue"]),
                               float(x.get("nes", np.nan))))
        E = pd.DataFrame(ev, columns=["tissue", "pos", "rsid", "p", "nes"]).drop_duplicates()
        if E.empty:
            log("  這個基因在這三個組織沒有顯著 eQTL，跳過")
            out["genes"][gene] = {"n_evariants": 0}
            continue

        def count_in(ivs):
            n = sum(int(((E.pos >= a) & (E.pos < b)).sum()) for a, b in ivs)
            bp = sum(b - a for a, b in ivs)
            return n, bp

        n_fg, bp_fg = count_in(fg_iv)
        bg_iv = list(zip(bgc.start, bgc.end))
        n_bg, bp_bg = count_in(bg_iv)
        # 寬度配對：背景取最寬的 N 個（N = 前景元件數 × 2）
        wide = bgc.nlargest(max(2 * len(fg_iv), 5), "width")
        n_w, bp_w = count_in(list(zip(wide.start, wide.end)))

        def rate(n, bp):
            return 1000.0 * n / bp if bp else np.nan

        tab = [[n_fg, max(bp_fg - n_fg, 0)], [n_bg, max(bp_bg - n_bg, 0)]]
        orr, pv = fisher_exact(tab, alternative="greater")
        tabw = [[n_fg, max(bp_fg - n_fg, 0)], [n_w, max(bp_w - n_w, 0)]]
        orw, pvw = fisher_exact(tabw, alternative="greater")

        hits = []
        for a, b in fg_iv:
            h = E[(E.pos >= a) & (E.pos < b)]
            for _, r in h.iterrows():
                hits.append(dict(gene=gene, element=f"{chrom}:{a}-{b}", **r.to_dict()))

        log("  前景 %6d bp, %3d 個 eQTL 變異 → %.3f /kb" % (bp_fg, n_fg, rate(n_fg, bp_fg)))
        log("  背景 %6d bp, %3d 個                → %.3f /kb  | OR %.2f, p = %.3g"
            % (bp_bg, n_bg, rate(n_bg, bp_bg), orr, pv))
        log("  寬度配對背景 %6d bp, %3d 個        → %.3f /kb  | OR %.2f, p = %.3g"
            % (bp_w, n_w, rate(n_w, bp_w), orw, pvw))
        if hits:
            log("  落在元件裡的變異:")
            for h in hits[:15]:
                log("    %s  %s  %s  p=%.2g  nes=%+.3f"
                    % (h["element"], h["rsid"], h["tissue"], h["p"], h["nes"]))
        allrows += hits

        out["genes"][gene] = {
            "n_evariants_total": int(len(E)),
            "foreground": {"n_elements": len(fg_iv), "bp": int(bp_fg),
                           "n_variants": int(n_fg), "per_kb": rate(n_fg, bp_fg)},
            "background_all": {"n_peaks": int(len(bgc)), "bp": int(bp_bg),
                               "n_variants": int(n_bg), "per_kb": rate(n_bg, bp_bg),
                               "odds_ratio": float(orr), "p": float(pv)},
            "background_width_matched": {"n_peaks": int(len(wide)), "bp": int(bp_w),
                                         "n_variants": int(n_w), "per_kb": rate(n_w, bp_w),
                                         "odds_ratio": float(orw), "p": float(pvw),
                                         "median_width_fg": float(np.median(
                                             [b - a for a, b in fg_iv])) if fg_iv else None,
                                         "median_width_bg": float(wide.width.median())},
            "elements": [{"chrom": c, "start": a, "end": b, "label": d}
                         for (c, a, b, d) in fg],
        }

    if allrows:
        pd.DataFrame(allrows).to_csv(f"{RES}/133_qtl_hits.csv", index=False)
    json.dump(out, open(f"{RES}/133_qtl_overlap.json", "w"), indent=2, ensure_ascii=False)
    log("\nDONE → results/133_qtl_overlap.json")


if __name__ == "__main__":
    main()
