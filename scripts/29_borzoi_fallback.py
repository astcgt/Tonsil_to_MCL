#!/usr/bin/env python
"""
步驟 29（AlphaGenome 2 的本機替代）：用 Borzoi 做 in-silico 元件刪除，
檢驗「MCL 中喪失/新增的那些調控元件，在序列層面是否足以解釋基因表現的改變」。

為什麼用 Borzoi 而不是 Enformer
-------------------------------
Enformer 的 5,313 條 human track 裡 RNA-seq 是 0 條，而且幾乎沒有 B 細胞 ATAC。
Borzoi 有 `RNA:GM12878`（淋巴母細胞）28 條、`ATAC:Memory B`、`ATAC:Plasma B`、
`DNASE:GM12878`，是開源模型裡最接近 B 細胞情境的。

重要技術細節
------------
* 輸入 524,288 bp，但**實際輸出只有中央 196,608 bp**（6,144 bins × 32 bp）。
* Borzoi 的輸出經過 squashed scaling，必須用 `undo_squashed_scale()` 還原
  才能算正確的 fold change（方向不受影響，但幅度會失真）。
* WT 與 mutant 必須讓目標基因落在**完全相同的 offset**，否則 bin 不對齊。

跑法：env/genomic/bin/python scripts/29_borzoi_fallback.py --gpu 4
"""
import os, sys, json, argparse
import numpy as np
import pandas as pd

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
RES = f"{BASE}/results"
GEN = f"{BASE}/data/genome"
FASTA = f"{GEN}/GRCh38.primary_assembly.genome.fa"

SEQ_LEN = 524_288
HALF = SEQ_LEN // 2
BIN = 32
OUT_BINS = 6144
OUT_SPAN = OUT_BINS * BIN              # 196,608

ELEMENTS = {
    # CIITA 有多重啟動子：pI 在 chr16:10,866,222（髓系）、pIII 在 chr16:10,877,200（B 細胞專用，
    # 承載 CIITA-201/202 等主要轉錄本）。評分以 pIII 為準。
    "CIITA": dict(chrom="chr16", tss=10877202, direction="lost", elements=[
        ("chr16", 10876429, 10877713, -1.23, "pIII 啟動子本體"),
        ("chr16", 10864696, 10865195, -2.30, "pI 上游 -1.5 kb（喪失最強）"),
        ("chr16", 10868192, 10869312, -1.12, "pI/pIII 之間"),
        ("chr16", 10855506, 10856375, -1.54, "遠端上游 -21 kb"),
        ("chr16", 10846567, 10847741, -1.24, "遠端上游 -30 kb"),
        ("chr16", 10852148, 10852628, -1.21, "遠端上游 -25 kb"),
        ("chr16", 10891088, 10892072, -1.09, "基因內 +14 kb"),
        ("chr16", 10911702, 10913211, -1.84, "基因內 +35 kb"),
    ]),
    "SOX11": dict(chrom="chr2", tss=5692384, direction="gained", elements=[
        ("chr2", 5695637, 5697306, +3.66, "基因內新生元件 +3.3 kb"),
        ("chr2", 5725875, 5726268, +2.64, "下游增強子 +33.5 kb"),
        ("chr2", 5692044, 5692374, +2.53, "啟動子 -0.3 kb"),
    ]),
}


def log(*a):
    print("[29]", *a, flush=True)


def di_shuffle(seq, rng):
    """保留雙核苷酸組成的打散（Altschul-Erikson），作為『破壞元件但不改變鹼基組成』的對照。"""
    s = list(seq)
    # 簡化版：以 2-mer 為單位重排
    pairs = [s[i:i + 2] for i in range(0, len(s) - 1, 2)]
    rng.shuffle(pairs)
    out = "".join("".join(p) for p in pairs)
    return (out + seq[len(out):])[:len(seq)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gpu", default="4")
    ap.add_argument("--repo", default="Genentech/borzoi-model")
    ap.add_argument("--ckpt", default="human_rep0.ckpt")
    ap.add_argument("--out", default=f"{RES}/29_borzoi_results.json")
    args = ap.parse_args()
    os.environ["CUDA_VISIBLE_DEVICES"] = args.gpu

    import torch, pyfaidx
    torch.backends.cudnn.enabled = False  # 本機 cuDNN 與 torch 2.4.1 的 conv1d 不相容
    import grelu.resources
    log("torch", torch.__version__, "cuda", torch.cuda.is_available())

    model = grelu.resources.load_model(repo_id=args.repo, filename=args.ckpt)
    model = model.eval().to("cuda" if torch.cuda.is_available() else "cpu")
    tasks = pd.DataFrame(model.data_params["tasks"])
    log("Borzoi tracks:", len(tasks), "| 欄位:", list(tasks.columns)[:6])

    # 挑 B 細胞相關 track
    txt = tasks.astype(str).agg(" ".join, axis=1)
    sel = {}
    sel["RNA_GM12878"] = np.flatnonzero(txt.str.contains("RNA:GM12878", case=False).values)
    sel["ATAC_B"] = np.flatnonzero(txt.str.contains("ATAC:Memory B|ATAC:Plasma B|ATAC:Naive B",
                                                    case=False, regex=True).values)
    sel["DNASE_GM12878"] = np.flatnonzero(txt.str.contains("DNASE:GM12878", case=False).values)
    sel["CAGE_Bcell"] = np.flatnonzero(txt.str.contains("CAGE:CD19|CAGE:tonsil", case=False).values)
    for k, v in sel.items():
        log(f"  track 群組 {k}: {len(v)} 條")
    json.dump({k: v.tolist() for k, v in sel.items()},
              open(f"{RES}/29_borzoi_track_idx.json", "w"))

    fa = pyfaidx.Fasta(FASTA, as_raw=True, sequence_always_upper=True)
    rng = np.random.default_rng(0)

    def predict(seq):
        import grelu.sequence.format as gsf
        with torch.no_grad():
            p = model.predict_on_seqs(seq, device=0 if torch.cuda.is_available() else "cpu")
        return np.asarray(p).squeeze()        # (n_tracks, 6144)

    def undo(x):
        """反轉 Borzoi 的 squashed scale（sum_sqrt + clip_soft=384）。"""
        try:
            from borzoi_pytorch.pytorch_borzoi_utils import undo_squashed_scale
            return np.asarray(undo_squashed_scale(torch.tensor(x))).astype(np.float64)
        except Exception:
            cs = 384.0
            y = np.array(x, dtype=np.float64)
            m = y > cs
            y[m] = (y[m] - cs) ** 2 + cs
            return y ** 2

    results = {}
    for gene, cfg in ELEMENTS.items():
        chrom, tss = cfg["chrom"], cfg["tss"]
        win_start = tss - HALF
        wt = str(fa[chrom][win_start:win_start + SEQ_LEN])
        assert len(wt) == SEQ_LEN, len(wt)
        out_start = win_start + (SEQ_LEN - OUT_SPAN) // 2
        tss_bin = (tss - out_start) // BIN
        log(f"\n=== {gene} {chrom}:{tss} | 輸出視窗 {out_start}-{out_start+OUT_SPAN} | TSS bin {tss_bin} ===")

        pw = undo(predict(wt))
        def score(P):
            lo, hi = max(0, tss_bin - 64), min(OUT_BINS, tss_bin + 64)   # TSS ±2 kb
            return {k: float(P[v, lo:hi].sum()) if len(v) else np.nan for k, v in sel.items()}
        base = score(pw)
        results[f"{gene}|WT"] = base
        log("  WT:", {k: round(v, 1) for k, v in base.items()})

        for (c, s, e, lfc, label) in cfg["elements"]:
            i0, i1 = s - win_start, e - win_start
            if not (0 <= i0 < i1 <= SEQ_LEN):
                log(f"  [skip] {label} 超出輸入視窗"); continue
            in_out_window = (out_start <= s) and (e <= out_start + OUT_SPAN)
            mut = wt[:i0] + di_shuffle(wt[i0:i1], rng) + wt[i1:]
            assert len(mut) == SEQ_LEN
            pm = undo(predict(mut))
            sc = score(pm)
            d = {k: (np.log2((sc[k] + 1) / (base[k] + 1)) if base.get(k) is not None else np.nan)
                 for k in sc}
            tag = f"{gene}|{label}"
            results[tag] = dict(abs_scores=sc, log2FC_vs_WT=d,
                                observed_atac_lfc=lfc, element=f"{c}:{s}-{e}",
                                element_in_output_window=bool(in_out_window))
            log(f"  刪除 {label:22s} → log2FC: " +
                " ".join(f"{k}={d[k]:+.3f}" for k in d if not np.isnan(d[k])))

        # 全部元件同時刪除
        mut = wt
        for (c, s, e, lfc, label) in cfg["elements"]:
            i0, i1 = s - win_start, e - win_start
            if 0 <= i0 < i1 <= SEQ_LEN:
                mut = mut[:i0] + di_shuffle(wt[i0:i1], rng) + mut[i1:]
        pm = undo(predict(mut))
        sc = score(pm)
        d = {k: np.log2((sc[k] + 1) / (base[k] + 1)) for k in sc}
        results[f"{gene}|ALL_elements"] = dict(abs_scores=sc, log2FC_vs_WT=d)
        log(f"  刪除【全部元件】 → log2FC: " +
            " ".join(f"{k}={d[k]:+.3f}" for k in d if not np.isnan(d[k])))

    json.dump(results, open(args.out, "w"), indent=1, default=float)
    log("\n寫出", args.out)


if __name__ == "__main__":
    main()
