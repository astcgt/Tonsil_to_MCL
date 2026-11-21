#!/usr/bin/env python
"""
步驟 97：AlphaGenome 2 的兩個未做測試 —— 甲基化 與 染色質接觸。

這一步要回答報告第七節問題 1 剩下的兩項，以及補上 SLC7A11 故事裡最大的一個洞。

Q1（接觸圖）：SLC7A11 的遠端新生元件「碰得到」啟動子嗎？
---------------------------------------------------------
目前「元件 → SLC7A11」這個連結完全建立在**相關性**上（peak-gene correlation）。
最遠的元件在 +164 kb。如果它和啟動子中間有一道 TAD 邊界，相關性就很可能是假的。
接觸圖可以直接檢驗這件事。

讀出設計（距離配對是關鍵）：
  接觸強度隨距離單調衰減，所以絕對值毫無意義。
  對每個元件，取「同一條對角線（相同 bin 距離）」的所有 bin 對當虛無分布，
  算出 observed/expected 與百分位。這樣距離就被完全抵銷。
  決定性對照 = 同窗格內「未重接線」的 peak，用完全相同的方式計算。

另外計算絕緣分數（insulation score）剖面，檢查元件與啟動子之間有沒有邊界。

Q2（甲基化）：那五個元件在正常 B 細胞裡是不是「被甲基化關著」？
--------------------------------------------------------------
若是 → 多出一條治療邏輯，而且是**臨床警訊**：
去甲基化藥物（azacitidine / decitabine，血液腫瘤常用）反而可能打開這條
鐵死亡防禦軸線。這是可檢驗、且文獻上沒人提的預測。

軌道等級（這決定結論能說多重）
------------------------------
`supported_conditions()` 把每個條件 × 輸出型別標成 observed / imputed_supported /
imputed_exploratory 三級。我們只用 **observed** 級的淋巴系條件：
  - blood B cell (UBERON:0000178_CL:0000236)：methylation 與 contact_maps 都是 observed
  - GM12878 (EFO:0002784)：contact_maps observed（B 淋巴母細胞，第二個獨立軌道）
  - blood memory B cell (UBERON:0000178_CL:0000787)：methylation observed（第二個獨立軌道）
B 細胞的 contact_maps 在「組織+細胞型」條件下只有 imputed_exploratory，
所以**刻意改用 blood B cell 這個 observed 條件**，而不是 tonsil naive B cell。

⚠️ 必須在報告裡寫明的限制：這些軌道都是**正常 B 細胞**，不是 MCL。
所以這一步問的是「這段距離在 B 細胞的三維結構裡搭不搭得起來」，
不是「MCL 裡接觸有沒有改變」。前者才是我們需要的 sanity check。

金鑰只從 .secrets/alphagenome_key 讀取，不寫入任何輸出。
"""
import os, sys, json
import numpy as np
import pandas as pd

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
RES = f"{BASE}/results"
_KEY = open(f"{BASE}/.secrets/alphagenome_key").read().strip()

from alphagenome.data import genome
from alphagenome.models.v2 import dna_client, dna_output

RNG = np.random.default_rng(0)
TSS = {"SLC7A11": ("chr4", 138242349), "MAFG": ("chr17", 81927735)}
# 校準錨點：CpG 島啟動子（應該幾乎沒有甲基化）
ANCHOR = {"GAPDH_promoter": ("chr12", 6534517), "ACTB_promoter": ("chr7", 5530601)}

CONTACT_COND = {"blood B cell": "UBERON:0000178_CL:0000236",
                "GM12878": "EFO:0002784"}
METHYL_COND = {"blood B cell": "UBERON:0000178_CL:0000236",
               "blood memory B cell": "UBERON:0000178_CL:0000787"}

LEN_CONTACT = 1048576
LEN_METHYL = 524288


def log(*a):
    print("[97]", *a, flush=True)


# --------------------------------------------------------------- 元件座標
def features(gene, chrom, tss, half):
    """回傳 (label, kind, midpoint, dist) 清單：新生元件 + 未重接線對照。"""
    rw = pd.read_csv(f"{RES}/52_rewiring_annotated.csv.gz")
    g = rw[(rw.gene == gene) & (rw.chrom == chrom)]
    lo, hi = tss - half, tss + half
    out = []
    for kind, sub in [("de_novo", g[g.gained]),
                      ("control", g[(~g.gained) & (~g.lost)])]:
        for _, r in sub.iterrows():
            mid = (int(r.start) + int(r.end)) // 2
            if not (lo + 5000 < mid < hi - 5000):
                continue
            out.append(dict(label=f"{kind}_{int(r.peak_id)}", kind=kind, mid=mid,
                            dist=int(mid - tss), width=int(r.end) - int(r.start),
                            dr=float(r.dr)))
    return out


# --------------------------------------------------------------- CpG 密度
def cpg_counts(chrom, mids, half=250):
    import pyfaidx
    fa = pyfaidx.Fasta(f"{BASE}/data/genome/GRCh38.primary_assembly.genome.fa",
                       as_raw=True, sequence_always_upper=True)
    out = []
    for m in mids:
        s = str(fa[chrom][m - half:m + half + 1])
        out.append(dict(n_cpg=s.count("CG"), gc=(s.count("G") + s.count("C")) / max(len(s), 1)))
    return out


# =============================================================== PART 1 接觸圖
def diag_null(M, d, i, j, guard=2):
    """同一條對角線（bin 距離 d）的所有 bin 對，排除與 i/j 重疊的。"""
    n = M.shape[0]
    a = np.arange(0, n - d)
    b = a + d
    ok = (np.abs(a - i) > guard) & (np.abs(a - j) > guard) & \
         (np.abs(b - i) > guard) & (np.abs(b - j) > guard)
    v = M[a[ok], b[ok]]
    return v[np.isfinite(v)]


def insulation(M, w=10):
    """對角菱形絕緣分數：越低 = 越像邊界。"""
    n = M.shape[0]
    out = np.full(n, np.nan)
    for k in range(w, n - w):
        blk = M[k - w:k, k + 1:k + w + 1]
        out[k] = np.nanmean(blk)
    med = np.nanmedian(out)
    return np.log2((out + 1e-9) / (med + 1e-9))


def run_contact(model):
    rows, ins_rows = [], []
    for gene, (chrom, tss) in TSS.items():
        iv = genome.Interval(chrom, tss - LEN_CONTACT // 2, tss + LEN_CONTACT // 2)
        feats = features(gene, chrom, tss, LEN_CONTACT // 2)
        log(f"\n=== 接觸圖 {gene} ({chrom}:{tss:,}) | 新生 "
            f"{sum(f['kind']=='de_novo' for f in feats)} / 對照 "
            f"{sum(f['kind']=='control' for f in feats)} ===")
        for cname, cid in CONTACT_COND.items():
            try:
                pr = model.predict_interval(interval=iv, condition=cid,
                                            requested_outputs=[dna_output.OutputType.CONTACT_MAPS])
            except Exception as e:
                log(f"  {cname} 失敗: {type(e).__name__}: {str(e)[:200]}")
                continue
            cm = getattr(pr, "contact_maps", None)
            if cm is None:
                log(f"  {cname}: 回傳沒有 contact_maps"); continue
            V = np.asarray(cm.values)
            res = cm.resolution
            log(f"  {cname}: shape {V.shape}, 解析度 {res} bp, "
                f"軌道數 {V.shape[-1] if V.ndim == 3 else 1}")
            M = np.nanmean(V, axis=-1) if V.ndim == 3 else V
            M = (M + M.T) / 2.0
            n = M.shape[0]
            ipro = (tss - iv.start) // res
            if not (0 <= ipro < n):
                log("  啟動子落在範圍外，跳過"); continue

            ins = insulation(M)
            for k, val in enumerate(ins):
                if np.isfinite(val):
                    ins_rows.append(dict(gene=gene, condition=cname,
                                         pos=int(iv.start + k * res), insulation=round(float(val), 4)))

            for f in feats:
                j = (f["mid"] - iv.start) // res
                if not (0 <= j < n) or abs(j - ipro) < 2:
                    continue
                d = abs(j - ipro)
                v = float(M[ipro, j])
                nl = diag_null(M, d, ipro, j)
                if len(nl) < 20 or not np.isfinite(v):
                    continue
                lo, hi = (min(ipro, j), max(ipro, j))
                seg = ins[lo:hi + 1]
                seg = seg[np.isfinite(seg)]
                rows.append(dict(
                    gene=gene, condition=cname, label=f["label"], kind=f["kind"],
                    dist=f["dist"], width=f["width"], dr=round(f["dr"], 3),
                    bin_dist=int(d), contact=round(v, 6),
                    exp_median=round(float(np.median(nl)), 6),
                    oe=round(float(v / max(np.median(nl), 1e-12)), 4),
                    pct=round(float(np.mean(nl < v)), 4), n_null=len(nl),
                    min_insulation_between=round(float(seg.min()), 4) if len(seg) else np.nan))
    return pd.DataFrame(rows), pd.DataFrame(ins_rows)


# =============================================================== PART 2 甲基化
def run_methyl(model):
    rows = []
    targets = {**{g: v for g, v in TSS.items()}, **ANCHOR}
    for gene, (chrom, tss) in targets.items():
        iv = genome.Interval(chrom, tss - LEN_METHYL // 2, tss + LEN_METHYL // 2)
        if gene in TSS:
            feats = features(gene, chrom, tss, LEN_METHYL // 2)
        else:
            feats = [dict(label=f"anchor_{gene}", kind="anchor_cpg_island",
                          mid=tss, dist=0, width=501, dr=np.nan)]
        # 全窗格隨機對照（代表基因體背景，預期高甲基化）
        for m in RNG.integers(iv.start + 2000, iv.end - 2000, 60):
            feats.append(dict(label=f"random_{m}", kind="random_background",
                              mid=int(m), dist=int(m - tss), width=501, dr=np.nan))
        cg = cpg_counts(chrom, [f["mid"] for f in feats])
        for f, c in zip(feats, cg):
            f.update(c)
        log(f"\n=== 甲基化 {gene} | 特徵 {len(feats)} 個 ===")
        for cname, cid in METHYL_COND.items():
            try:
                pr = model.predict_interval(interval=iv, condition=cid,
                                            requested_outputs=[dna_output.OutputType.METHYLATION])
            except Exception as e:
                log(f"  {cname} 失敗: {type(e).__name__}: {str(e)[:200]}")
                continue
            me = getattr(pr, "methylation", None)
            if me is None:
                log(f"  {cname}: 回傳沒有 methylation"); continue
            V = np.asarray(me.values)
            res = me.resolution
            log(f"  {cname}: shape {V.shape}, 解析度 {res} bp")
            v = np.nanmean(V, axis=-1) if V.ndim == 2 else V

            def seg(mid, half=250):
                i0 = max(0, (mid - half - iv.start) // res)
                i1 = min(len(v), (mid + half - iv.start) // res + 1)
                s = v[i0:i1]
                s = s[np.isfinite(s)]
                return float(s.mean()) if len(s) else np.nan

            for f in feats:
                rows.append(dict(gene=gene, condition=cname, label=f["label"],
                                 kind=f["kind"], dist=f["dist"], width=f["width"],
                                 dr=f["dr"], n_cpg=f["n_cpg"], gc=round(f["gc"], 3),
                                 methyl=round(seg(f["mid"]), 5)))
    return pd.DataFrame(rows)


# =============================================================== 主程式
model = dna_client.create(_KEY)
log("連線建立（金鑰未輸出）")

ct, ins = run_contact(model)
if not ct.empty:
    ct.to_csv(f"{RES}/97_contact_oe.csv", index=False)
    ins.to_csv(f"{RES}/97_insulation.csv.gz", index=False)

me = run_methyl(model)
if not me.empty:
    me.to_csv(f"{RES}/97_methylation.csv", index=False)

# --------------------------------------------------------------- 統計
from scipy.stats import mannwhitneyu
summary = {}
log("\n" + "=" * 72)
log("=== Q1：遠端新生元件碰得到啟動子嗎？（距離配對後的 observed/expected）===")
for (gene, cond), s in ct.groupby(["gene", "condition"]) if not ct.empty else []:
    dn = s[s.kind == "de_novo"]
    cl = s[s.kind == "control"]
    if len(dn) < 2 or len(cl) < 3:
        log(f"\n{gene} / {cond}: 樣本不足（新生 {len(dn)}、對照 {len(cl)}）"); continue
    u = mannwhitneyu(dn.oe, cl.oe, alternative="greater")
    k = f"{gene}|{cond}"
    summary[k] = dict(n_de_novo=len(dn), n_control=len(cl),
                      oe_de_novo=round(float(dn.oe.median()), 3),
                      oe_control=round(float(cl.oe.median()), 3),
                      pct_de_novo=round(float(dn.pct.median()), 3),
                      p=float(u.pvalue))
    log(f"\n{gene} / {cond}")
    log(f"  新生元件 O/E 中位 {dn.oe.median():.3f} (n={len(dn)}, 距離配對百分位 "
        f"{dn.pct.median():.2f}) vs 未重接線對照 {cl.oe.median():.3f} (n={len(cl)})"
        f"  Mann-Whitney p = {u.pvalue:.4g}")
    log(dn[["label", "dist", "bin_dist", "oe", "pct", "min_insulation_between"]]
        .sort_values("dist").round(3).to_string(index=False))

log("\n" + "=" * 72)
log("=== Q2：正常 B 細胞裡，這些元件是被甲基化關著的嗎？===")
for (gene, cond), s in me.groupby(["gene", "condition"]) if not me.empty else []:
    dn = s[s.kind == "de_novo"]
    cl = s[s.kind == "control"]
    bg = s[s.kind == "random_background"]
    an = s[s.kind == "anchor_cpg_island"]
    if len(an):
        log(f"\n{gene} / {cond}: CpG 島錨點 {an.methyl.iloc[0]:.3f}"
            f"（預期接近 0）| 隨機背景中位 {bg.methyl.median():.3f}（預期高）")
        continue
    if len(dn) < 2 or len(cl) < 3:
        continue
    u = mannwhitneyu(dn.methyl, cl.methyl, alternative="greater")
    k = f"methyl|{gene}|{cond}"
    summary[k] = dict(n_de_novo=len(dn), methyl_de_novo=round(float(dn.methyl.median()), 4),
                      methyl_control=round(float(cl.methyl.median()), 4),
                      methyl_background=round(float(bg.methyl.median()), 4),
                      cpg_de_novo=int(dn.n_cpg.median()), p=float(u.pvalue))
    log(f"\n{gene} / {cond}")
    log(f"  新生元件甲基化中位 {dn.methyl.median():.4f} (n={len(dn)}, CpG 中位 "
        f"{dn.n_cpg.median():.0f}) vs 未重接線對照 {cl.methyl.median():.4f}"
        f" | 隨機背景 {bg.methyl.median():.4f}   p = {u.pvalue:.4g}")
    log(dn[["label", "dist", "n_cpg", "gc", "methyl"]].sort_values("dist").round(4).to_string(index=False))

json.dump(summary, open(f"{RES}/97_summary.json", "w"), indent=2, ensure_ascii=False)
log("\nDONE")
