#!/usr/bin/env python
"""
步驟 98：接觸圖的跨條件複製 + 寬度／距離混淆校正。

為什麼需要這一步
----------------
步驟 97 給出兩個互相矛盾的 observed 級結果：
  SLC7A11：GM12878 顯著（O/E 1.382, p=0.023），blood B cell 不顯著（1.008, p=0.74）
  MAFG：兩個條件都顯著，但**新生元件寬 3.4 倍（1402 vs 412 bp）且 O/E 與寬度正相關**
        → 和第六節第 15/16 列同一個寬度陷阱，結果不可信

所以這一步做兩件事：
1. **跨 21 個 observed 條件複製**，分成三組來判斷訊號是細胞型專一還是結構性：
     LCL 組（6）：GM12878 + 5 株 HapMap B 淋巴母細胞株（與 GM12878 同細胞型）
     原發血液組（9）：blood B cell、各種 T 細胞、單核球、K562
     非血液組（6）：HepG2、IMR-90、MCF-7、HCT116、Panc1、Caco-2（專一性對照）
2. **把寬度與距離從讀出裡扣掉**：在每個（基因 × 條件）內部，
   把 O/E、log 寬度、log 距離都轉成秩，用最小平方把後兩者迴歸掉，
   再比較殘差。這樣寬度與距離就不可能是結果的來源。

虛無分布仍然是「同一條對角線」，所以距離在第一層就已經抵銷；
秩迴歸是第二層保險，針對的是 O/E 本身隨距離／寬度系統性漂移。

金鑰只從 .secrets/alphagenome_key 讀取，不寫入任何輸出。
"""
import os, sys, time, json
import numpy as np
import pandas as pd

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
RES = f"{BASE}/results"
_KEY = open(f"{BASE}/.secrets/alphagenome_key").read().strip()

from alphagenome.data import genome
from alphagenome.models.v2 import dna_client, dna_output

TSS = {"SLC7A11": ("chr4", 138242349), "MAFG": ("chr17", 81927735)}
SEQLEN = 1048576

GROUPS = {
    "LCL(B淋巴母細胞)": {
        "GM12878": "EFO:0002784", "GM18519": "CLO:0024886", "GM18508": "CLO:0024919",
        "GM18498": "CLO:0024937", "GM18907": "CLO:0026238", "GM18861": "CLO:0026256"},
    "原發血液細胞": {
        "blood B cell": "UBERON:0000178_CL:0000236",
        "blood Th1": "UBERON:0000178_CL:0000545",
        "blood Th2": "UBERON:0000178_CL:0000546",
        "blood CD4 T": "UBERON:0000178_CL:0000624",
        "blood CD8 T": "UBERON:0000178_CL:0000625",
        "blood naive CD4 T": "UBERON:0000178_CL:0000895",
        "blood memory CD4 T": "UBERON:0000178_CL:0000897",
        "blood monocyte": "UBERON:0000178_CL:0001054",
        "K562": "EFO:0002067"},
    "非血液(專一性對照)": {
        "HepG2": "EFO:0001187", "IMR-90": "EFO:0001196", "MCF-7": "EFO:0001203",
        "HCT116": "EFO:0002824", "Panc1": "EFO:0002713", "Caco-2": "EFO:0001099"},
}
OUT = f"{RES}/98_contact_multi.csv"


def log(*a):
    print("[98]", *a, flush=True)


def features(gene, chrom, tss, half):
    rw = pd.read_csv(f"{RES}/52_rewiring_annotated.csv.gz")
    g = rw[(rw.gene == gene) & (rw.chrom == chrom)]
    lo, hi = tss - half, tss + half
    out = []
    for kind, sub in [("de_novo", g[g.gained]), ("control", g[(~g.gained) & (~g.lost)])]:
        for _, r in sub.iterrows():
            mid = (int(r.start) + int(r.end)) // 2
            if not (lo + 5000 < mid < hi - 5000):
                continue
            out.append(dict(label=f"{kind}_{int(r.peak_id)}", kind=kind, mid=mid,
                            dist=int(mid - tss), width=int(r.end) - int(r.start),
                            dr=float(r.dr)))
    return out


def diag_null(M, d, i, j, guard=2):
    n = M.shape[0]
    a = np.arange(0, n - d); b = a + d
    ok = (np.abs(a - i) > guard) & (np.abs(a - j) > guard) & \
         (np.abs(b - i) > guard) & (np.abs(b - j) > guard)
    v = M[a[ok], b[ok]]
    return v[np.isfinite(v)]


def predict_retry(model, iv, cid, tries=4):
    for k in range(tries):
        try:
            return model.predict_interval(
                interval=iv, condition=cid,
                requested_outputs=[dna_output.OutputType.CONTACT_MAPS])
        except Exception as e:
            if k == tries - 1:
                raise
            time.sleep(2 ** k * 3)
    return None


model = dna_client.create(_KEY)
log("連線建立（金鑰未輸出）")

rows = []
for gene, (chrom, tss) in TSS.items():
    iv = genome.Interval(chrom, tss - SEQLEN // 2, tss + SEQLEN // 2)
    feats = features(gene, chrom, tss, SEQLEN // 2)
    log(f"\n=== {gene}: 新生 {sum(f['kind']=='de_novo' for f in feats)} / "
        f"對照 {sum(f['kind']=='control' for f in feats)} ===")
    for grp, conds in GROUPS.items():
        for cname, cid in conds.items():
            try:
                pr = predict_retry(model, iv, cid)
            except Exception as e:
                log(f"  {grp}/{cname} 失敗: {type(e).__name__}: {str(e)[:120]}")
                continue
            cm = getattr(pr, "contact_maps", None)
            if cm is None:
                log(f"  {grp}/{cname}: 無 contact_maps"); continue
            V = np.asarray(cm.values); res = cm.resolution
            M = np.nanmean(V, axis=-1) if V.ndim == 3 else V
            M = (M + M.T) / 2.0
            n = M.shape[0]
            ipro = (tss - iv.start) // res
            got = 0
            for f in feats:
                j = (f["mid"] - iv.start) // res
                if not (0 <= j < n) or abs(j - ipro) < 2:
                    continue
                d = abs(j - ipro)
                v = float(M[ipro, j])
                nl = diag_null(M, d, ipro, j)
                if len(nl) < 20 or not np.isfinite(v):
                    continue
                rows.append(dict(gene=gene, group=grp, condition=cname,
                                 label=f["label"], kind=f["kind"], dist=f["dist"],
                                 width=f["width"], dr=round(f["dr"], 3), bin_dist=int(d),
                                 oe=round(float(v / max(np.median(nl), 1e-12)), 5),
                                 pct=round(float(np.mean(nl < v)), 4)))
                got += 1
            log(f"  {grp}/{cname}: {got} 個特徵")
            pd.DataFrame(rows).to_csv(OUT, index=False)     # 增量存檔

df = pd.DataFrame(rows)
df.to_csv(OUT, index=False)
log(f"\n共 {len(df):,} 筆 ({df.condition.nunique()} 個條件)")

# ------------------------------------------------- 寬度／距離秩迴歸校正
from scipy.stats import rankdata, mannwhitneyu, combine_pvalues


def adj_resid(s):
    """把 log 寬度與 log 距離（秩）從 O/E（秩）裡迴歸掉。"""
    y = rankdata(s.oe)
    X = np.c_[np.ones(len(s)), rankdata(np.log(s.width)), rankdata(np.log(s.bin_dist))]
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return y - X @ beta


res_rows = []
for (gene, grp, cond), s in df.groupby(["gene", "group", "condition"]):
    dn = s[s.kind == "de_novo"]; cl = s[s.kind == "control"]
    if len(dn) < 3 or len(cl) < 5:
        continue
    r = adj_resid(s)
    m = s.kind.values == "de_novo"
    u_raw = mannwhitneyu(dn.oe, cl.oe, alternative="greater")
    u_adj = mannwhitneyu(r[m], r[~m], alternative="greater")
    res_rows.append(dict(gene=gene, group=grp, condition=cond,
                         n_de_novo=len(dn), n_control=len(cl),
                         oe_de_novo=round(float(dn.oe.median()), 3),
                         oe_control=round(float(cl.oe.median()), 3),
                         p_raw=float(u_raw.pvalue), p_adj=float(u_adj.pvalue),
                         resid_delta=round(float(np.median(r[m]) - np.median(r[~m])), 2)))
R = pd.DataFrame(res_rows)
R.to_csv(f"{RES}/98_contact_stats.csv", index=False)

pd.set_option("display.width", 200)
log("\n" + "=" * 78)
log("=== 每個條件：新生元件 vs 未重接線對照（p_adj 已扣掉寬度與距離）===")
for gene, s in R.groupby("gene"):
    log(f"\n--- {gene} ---")
    log(s[["group", "condition", "oe_de_novo", "oe_control", "p_raw", "p_adj", "resid_delta"]]
        .sort_values(["group", "p_adj"]).to_string(index=False, float_format=lambda x: f"{x:.4g}"))

log("\n" + "=" * 78)
log("=== 按組彙總（Stouffer 合併 p_adj，同組條件視為獨立複製）===")
for (gene, grp), s in R.groupby(["gene", "group"]):
    nsig = int((s.p_adj < 0.05).sum())
    st = combine_pvalues(s.p_adj.clip(1e-300, 1), method="stouffer")
    log(f"{gene:8s} {grp:18s} n={len(s):2d} 條件 | 顯著 {nsig}/{len(s)} | "
        f"殘差差中位 {s.resid_delta.median():+.1f} | Stouffer p = {st.pvalue:.4g}")

log("\n" + "=" * 78)
log("=== 逐元件跨條件一致性（距離配對百分位的中位數）===")
for gene, s in df[df.kind == "de_novo"].groupby("gene"):
    t = s.pivot_table(index=["label", "dist"], columns="group", values="pct", aggfunc="median")
    log(f"\n--- {gene} ---")
    log(t.round(3).sort_index(level="dist").to_string())

log("\nDONE")
