#!/usr/bin/env python
"""
步驟 28：用 AlphaGenome 2 驗證「cis 調控異常」的因果可能性。

要回答的兩個問題
----------------
Q1 (CIITA 沉默)：MCL 中喪失的那幾段 CIITA 上游元件，若在序列層面被「破壞」，
   模型是否預測 CIITA 在 B 細胞中的表現會下降？下降幅度有多大？
   → 用 in-silico 刪除（把元件換成隨機/打散序列），比較 WT vs 刪除後的 RNA_SEQ 訊號。

Q2 (SOX11 活化)：MCL 中新生的那三段 SOX11 元件，是否足以驅動 SOX11 表現？
   → 反向操作：在正常 B 細胞 condition 下，比較帶有 / 不帶這些元件的序列。

需要 ALPHA_GENOME_API_KEY 環境變數。沒有 key 時本腳本會明確報錯並建議改跑
29_borzoi_fallback.py（本機 Borzoi，不需 key）。
"""
import os, sys, json, argparse
import numpy as np
import pandas as pd

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
RES = f"{BASE}/results"
GEN = f"{BASE}/data/genome"
FASTA = f"{GEN}/GRCh38.primary_assembly.genome.fa"

# 由步驟 22/26 的實測結果決定的元件座標
ELEMENTS = {
    "CIITA": {
        "chrom": "chr16", "tss": 10866222, "strand": "+",
        "lost_elements": [
            ("chr16", 10864696, 10865195, -2.30),   # 最強喪失，距 TSS -1.5 kb
            ("chr16", 10855506, 10856375, -1.54),   # -10.7 kb
            ("chr16", 10846567, 10847741, -1.24),   # -19.7 kb
            ("chr16", 10852148, 10852628, -1.21),   # -14.1 kb
            ("chr16", 10868192, 10869312, -1.12),   # +2.0 kb (gene body)
            ("chr16", 10876429, 10877713, -1.23),   # +10.2 kb
        ],
    },
    "SOX11": {
        "chrom": "chr2", "tss": 5692384, "strand": "+",
        "gained_elements": [
            ("chr2", 5695637, 5697306, +3.66),      # 最強新生，+3.3 kb
            ("chr2", 5725875, 5726268, +2.64),      # +33.5 kb
            ("chr2", 5692044, 5692374, +2.53),      # 啟動子 -0.3 kb
        ],
    },
}


def log(*a):
    print("[28]", *a, flush=True)


def shuffle_seq(s, rng):
    """保留雙核苷酸組成的打散（比單純隨機更嚴謹的陰性對照）。"""
    a = list(s)
    rng.shuffle(a)
    return "".join(a)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--condition", default=None,
                    help="AlphaGenome 2 的 composite ontology CURIE（B 細胞相關）")
    ap.add_argument("--list_conditions", action="store_true")
    ap.add_argument("--out", default=f"{RES}/28_alphagenome_results.json")
    args = ap.parse_args()

    key = os.environ.get("ALPHA_GENOME_API_KEY")
    if not key:
        print("""
========================================================================
需要 AlphaGenome 2 API key 才能執行本步驟。

取得方式（免費、非商業研究用）：
  1. 開 https://alphagenome.google/api ，用 Google 帳號登入後申請
  2. 設定環境變數：
       export ALPHA_GENOME_API_KEY='你的key'
  3. 重新執行本腳本

另外（有時效）：11/30 之前可以在 https://forms.gle/ftBfHGH7NUrCZEjN6
送出 use-case 表單爭取額外配額。本專案的題目（MCL 的 CIITA 沉默與
SOX11 異位活化之 cis 調控機制）很適合送。

沒有 key 時，請改跑不需要 key 的本機版本：
  env/genomic/bin/python scripts/29_borzoi_fallback.py
========================================================================
""")
        sys.exit(2)

    from alphagenome.data import genome
    from alphagenome.models.v2 import dna_client
    import pyfaidx

    model = dna_client.create(key)
    if args.list_conditions:
        conds = model.supported_conditions()
        df = pd.DataFrame(conds) if not isinstance(conds, pd.DataFrame) else conds
        df.to_csv(f"{RES}/28_supported_conditions.csv", index=False)
        log("可用 condition 數:", len(df))
        txt = df.astype(str).apply(lambda r: " ".join(r.values), axis=1)
        for kw in ["B cell", "lymphocyte", "tonsil", "lymph", "naive B", "memory B"]:
            hit = df[txt.str.contains(kw, case=False, na=False)]
            log(f"  含 '{kw}': {len(hit)}")
            if len(hit):
                log("\n" + hit.head(10).to_string())
        return

    assert args.condition, "請用 --list_conditions 先挑一個 B 細胞 condition，再用 --condition 指定"
    fa = pyfaidx.Fasta(FASTA, as_raw=True, sequence_always_upper=True)
    SEQLEN = dna_client.SEQUENCE_LENGTH_1MB
    HALF = SEQLEN // 2
    rng = np.random.default_rng(0)
    out = {}

    for gene, cfg in ELEMENTS.items():
        chrom, tss = cfg["chrom"], cfg["tss"]
        start = tss - HALF
        end = start + SEQLEN
        wt = str(fa[chrom][start:end])
        assert len(wt) == SEQLEN, len(wt)
        log(f"\n=== {gene} {chrom}:{start}-{end} ===")

        base = model.predict_sequence(
            sequence=wt, condition=args.condition,
            requested_outputs=[dna_client.OutputType.RNA_SEQ,
                               dna_client.OutputType.SC_ATAC,
                               dna_client.OutputType.CHIP_HISTONE],
            interval=genome.Interval(chrom, start, end))
        out[f"{gene}|WT"] = _summarize(base, tss, start)

        key_list = "lost_elements" if "lost_elements" in cfg else "gained_elements"
        for (c, s, e, lfc) in cfg[key_list]:
            assert c == chrom
            i0, i1 = s - start, e - start
            if not (0 <= i0 < i1 <= SEQLEN):
                log("  元件超出視窗，跳過", c, s, e); continue
            mut = wt[:i0] + shuffle_seq(wt[i0:i1], rng) + wt[i1:]
            assert len(mut) == SEQLEN
            pred = model.predict_sequence(
                sequence=mut, condition=args.condition,
                requested_outputs=[dna_client.OutputType.RNA_SEQ,
                                   dna_client.OutputType.SC_ATAC,
                                   dna_client.OutputType.CHIP_HISTONE],
                interval=genome.Interval(chrom, start, end))
            tag = f"{gene}|del_{c}:{s}-{e}"
            out[tag] = _summarize(pred, tss, start)
            out[tag]["observed_atac_lfc"] = lfc
            log(f"  {tag}: 完成")

    json.dump(out, open(args.out, "w"), indent=1, default=float)
    log("寫出", args.out)


def _summarize(pred, tss, win_start):
    """把 AlphaGenome 的輸出壓成『TSS ±2 kb 的訊號總和』。"""
    d = {}
    for name in ["rna_seq", "sc_atac", "chip_histone"]:
        o = getattr(pred, name, None)
        if o is None:
            continue
        v = np.asarray(o.values)
        res = max(1, int(round((o.interval.width if hasattr(o, "interval") else v.shape[0]) / v.shape[0])))
        c = (tss - win_start) // res
        lo, hi = max(0, c - 2000 // res), min(v.shape[0], c + 2000 // res)
        d[name + "_tss2kb"] = float(np.nansum(v[lo:hi]))
        d[name + "_total"] = float(np.nansum(v))
    return d


if __name__ == "__main__":
    main()
