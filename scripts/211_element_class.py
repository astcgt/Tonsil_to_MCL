#!/usr/bin/env python
"""步驟 211（新方向，淺層）：腫瘤新開的元件，和 NRF2 用的標準 ARE 是同一類嗎？

這一題從哪裡來
--------------
MafG 的蛋白質結合劑路線在步驟 210 之後三個軸全部關閉：
  C5 `MAFG` 在任何細胞都不是遺傳依賴（DepMap 中位 −0.157、MCL 0/3、OR 0.86）
  C6 設計品分不開 MafG / MafF / MafK（三個共摺疊模型差距 −0.011 ~ +0.051；MD ΔQ = −0.016）
  D5 夥伴選擇性在幾何上不存在 —— **NRF2 自己就遮住 BACH1 足跡的 82.1%**，
     因為兩者都用同一段白胺酸拉鍊和 small Maf 異二聚化，天生互斥；
     那 9 個「BACH1 專屬」位置只是共用介面的邊緣差異（而且家族內完全保守）

所以選擇性不可能來自蛋白質。剩下唯一的地方是 **DNA 元件本身**：
如果腫瘤新開的那些 MARE 在序列類別上就和 NRF2 用的標準 ARE 不同，
那麼誘餌寡核苷酸或序列專一的藥物就有可能只打腫瘤的程式、放過正常的細胞保護程式。

元件類別怎麼分
--------------
small Maf 複合體讀兩種相關但不同的元件：
  **TRE（AP-1 位點）** `TGACTCA` —— 7 bp，不需要側翼 GC
  **ARE / 完整 MARE** `TGACnnnGC` —— small Maf 的 extended homology region（EHR）
      讀的就是那個側翼 **GC**；NRF2–smallMaf 的標準結合位點是這一種
關鍵的區別不在核心而在**側翼有沒有 GC**。

事前指定（雙尾，因為我沒有理由預設方向）
----------------------------------------
H   腫瘤新開的 peak 與穩定開著的 peak 之間，`TGAC` 核心裡「完整 ARE」的比例不同。
    報告方向，但不預設方向。
P1  差異必須在**寬度與 GC 含量配對**之後仍然存在。
    這個專案已經被同一個陷阱騙過三次（motif 長度、GC 含量、peak 寬度），
    每一次都產生漂亮而錯誤的結果，所以配對不是加分項而是必要條件。
P2  差異必須大於「把前景／背景標籤在配對層內打亂」的經驗虛無分佈。
P3（陰性對照）組成相同但不是已知 TF 位點的樣式 `GTCAnnnCG`
    不該顯示同樣的差異。
P4（健全性檢查）前景與背景配對後的寬度中位數與 GC 中位數必須幾乎相同，
    否則配對失敗，結果不准解讀。

輸出 results/211_*
"""
import os
import json
import re

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "8"

import numpy as np
import pandas as pd

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
RES = f"{BASE}/results"
FASTA = f"{BASE}/data/genome/GRCh38.primary_assembly.genome.fa"
N_PERM = 5000
RNG = np.random.default_rng(0)

# TGAC 核心；完整 ARE 需要 +7/+8 的 GC
CORE = "TGAC"
ARE = re.compile(r"TGAC[ACGT]{3}GC")
ARE_RC = re.compile(r"GC[ACGT]{3}GTCA")
CTRL = re.compile(r"GTCA[ACGT]{3}CG")          # 組成相同、非已知位點
CTRL_RC = re.compile(r"CG[ACGT]{3}TGAC")


def log(*a):
    print("[211]", *a, flush=True)


def revcomp(s):
    return s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def count_pat(seq, pat, pat_rc):
    return len(pat.findall(seq)) + len(pat_rc.findall(seq))


def main():
    import pyfaidx
    fa = pyfaidx.Fasta(FASTA, as_raw=True, sequence_always_upper=True)

    d = pd.read_csv(f"{RES}/52_rewiring_annotated.csv.gz")
    log("重接線表 %d 列" % len(d))
    # 以 peak 為單位去重（一個 peak 可以連到多個基因）
    pk = (d.groupby("pid")
            .agg(chrom=("chrom", "first"), start=("start", "first"),
                 end=("end", "first"), gained=("gained", "max"),
                 lost=("lost", "max"))
            .reset_index())
    pk["width"] = pk.end - pk.start
    fg = pk[(pk.gained == 1) | (pk.gained is True)]
    fg = pk[pk.gained.astype(bool)]
    bg = pk[(~pk.gained.astype(bool)) & (~pk.lost.astype(bool))]
    log("前景（腫瘤新開）%d 個 peak｜背景（穩定開著）%d 個" % (len(fg), len(bg)))

    def seqs_of(df):
        out = []
        for c, s, e in zip(df.chrom, df.start, df.end):
            try:
                out.append(str(fa[c][int(s):int(e)]))
            except Exception:
                out.append("")
        return out

    def gc(s):
        if not s:
            return np.nan
        n = s.count("G") + s.count("C")
        return n / len(s)

    log("取序列…")
    fg = fg.copy(); bg = bg.copy()
    fg["seq"] = seqs_of(fg); bg["seq"] = seqs_of(bg)
    fg = fg[fg.seq.str.len() > 50]; bg = bg[bg.seq.str.len() > 50]
    fg["gc"] = [gc(s) for s in fg.seq]; bg["gc"] = [gc(s) for s in bg.seq]

    # --- 寬度 + GC 配對 ---------------------------------------------------
    log("做寬度 × GC 配對…")
    bg = bg.dropna(subset=["gc"])
    bg_sorted = bg.sort_values(["width", "gc"]).reset_index(drop=True)
    wq = np.quantile(bg_sorted.width, np.linspace(0, 1, 21))
    gq = np.quantile(bg_sorted.gc, np.linspace(0, 1, 21))
    bg_sorted["wbin"] = np.digitize(bg_sorted.width, wq)
    bg_sorted["gbin"] = np.digitize(bg_sorted.gc, gq)
    fgm = fg.dropna(subset=["gc"]).copy()
    fgm["wbin"] = np.digitize(fgm.width, wq)
    fgm["gbin"] = np.digitize(fgm.gc, gq)
    pool = {k: v.index.values for k, v in bg_sorted.groupby(["wbin", "gbin"])}
    pick, keep = [], []
    for i, r in fgm.iterrows():
        cand = pool.get((r.wbin, r.gbin))
        if cand is None or len(cand) == 0:
            continue
        pick.append(int(RNG.choice(cand)))
        keep.append(i)
    fgm = fgm.loc[keep]
    bgm = bg_sorted.loc[pick]
    log("配對後：前景 %d、背景 %d" % (len(fgm), len(bgm)))

    sane = {"fg_width_median": float(fgm.width.median()),
            "bg_width_median": float(bgm.width.median()),
            "fg_gc_median": float(fgm.gc.median()),
            "bg_gc_median": float(bgm.gc.median())}
    sane["pass"] = bool(abs(sane["fg_width_median"] - sane["bg_width_median"])
                        / sane["fg_width_median"] < 0.05
                        and abs(sane["fg_gc_median"] - sane["bg_gc_median"]) < 0.01)
    log("P4 配對健全性：寬度中位 %.0f vs %.0f｜GC 中位 %.3f vs %.3f → %s"
        % (sane["fg_width_median"], sane["bg_width_median"],
           sane["fg_gc_median"], sane["bg_gc_median"],
           "過" if sane["pass"] else "不過"))

    def stats_of(df, pat, pat_rc):
        ncore = np.array([s.count(CORE) + revcomp(s).count(CORE) for s in df.seq])
        nare = np.array([count_pat(s, pat, pat_rc) for s in df.seq])
        tot_core, tot_are = int(ncore.sum()), int(nare.sum())
        return tot_core, tot_are, (tot_are / tot_core if tot_core else np.nan)

    out = {"n_fg": int(len(fgm)), "n_bg": int(len(bgm)),
           "P4_matching_sanity": sane}
    results = {}
    for nm, pat, pat_rc in [("ARE_full", ARE, ARE_RC),
                            ("control_pattern", CTRL, CTRL_RC)]:
        cf, af, rf = stats_of(fgm, pat, pat_rc)
        cb, ab, rb = stats_of(bgm, pat, pat_rc)
        results[nm] = {"fg_cores": cf, "fg_hits": af, "fg_ratio": rf,
                       "bg_cores": cb, "bg_hits": ab, "bg_ratio": rb,
                       "difference": rf - rb,
                       "ratio_of_ratios": (rf / rb) if rb else np.nan}
        log("%-16s 前景 %d/%d = %.4f｜背景 %d/%d = %.4f｜差 %+.4f（倍率 %.3f）"
            % (nm, af, cf, rf, ab, cb, rb, rf - rb,
               results[nm]["ratio_of_ratios"]))
    out["patterns"] = results

    # --- P2：配對層內打亂標籤 --------------------------------------------
    both = pd.concat([fgm.assign(grp=1), bgm.assign(grp=0)])
    obs = results["ARE_full"]["difference"]
    cores = np.array([s.count(CORE) + revcomp(s).count(CORE) for s in both.seq])
    hits = np.array([count_pat(s, ARE, ARE_RC) for s in both.seq])
    g = both.grp.values.copy()
    null = np.empty(N_PERM)
    for i in range(N_PERM):
        gp = RNG.permutation(g)
        c1, h1 = cores[gp == 1].sum(), hits[gp == 1].sum()
        c0, h0 = cores[gp == 0].sum(), hits[gp == 0].sum()
        null[i] = (h1 / c1 if c1 else np.nan) - (h0 / c0 if c0 else np.nan)
    emp_p = float((np.abs(null) >= abs(obs)).mean())
    out["P2_label_permutation"] = {
        "n_perm": N_PERM, "null_mean": float(np.nanmean(null)),
        "null_sd": float(np.nanstd(null)), "observed": float(obs),
        "empirical_p_two_sided": emp_p, "pass": bool(emp_p < 0.05)}
    log("P2 打亂標籤：虛無 %.5f ± %.5f；觀察 %+.5f → 雙尾經驗 p = %.4f → %s"
        % (np.nanmean(null), np.nanstd(null), obs, emp_p,
           "過" if out["P2_label_permutation"]["pass"] else "不過"))

    ctrl_diff = results["control_pattern"]["difference"]
    out["P3_negative_control"] = {
        "control_difference": float(ctrl_diff),
        "target_difference": float(obs),
        "pass": bool(abs(ctrl_diff) < abs(obs))}
    log("P3 陰性對照差異 %+.5f vs 標的 %+.5f → %s"
        % (ctrl_diff, obs, "過" if out["P3_negative_control"]["pass"] else "不過"))

    passed = [k for k, v in out.items() if isinstance(v, dict) and v.get("pass") is True]
    failed = [k for k, v in out.items() if isinstance(v, dict) and v.get("pass") is False]
    out["summary"] = {"passed": passed, "failed": failed,
                      "verdict": ("worth_deepening"
                                  if (sane["pass"] and emp_p < 0.05
                                      and abs(ctrl_diff) < abs(obs))
                                  else "not_worth_deepening")}
    log("\n通過：%s｜未通過：%s" % (", ".join(passed) or "（無）",
                                   ", ".join(failed) or "（無）"))
    log("判定：%s" % out["summary"]["verdict"])
    json.dump(out, open(f"{RES}/211_element_class.json", "w"),
              indent=2, ensure_ascii=False)
    log("DONE → results/211_*")


if __name__ == "__main__":
    main()
