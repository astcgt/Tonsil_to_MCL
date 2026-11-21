#!/usr/bin/env python
"""
步驟 54：MCL 全基因組體細胞突變 vs 重接線的調控元件。

問題：腫瘤專有的那些「新生調控元件」裡，有沒有復發性的非編碼體細胞突變？
資料：22 位 MCL 病人 / 35 個樣本的 WGS（249,321 個 SNV，Zenodo 16098829）

做法
----
1. 把突變對到圖譜 peak 上
2. 比較三組 peak 的突變密度：新生連結(gained)、斷裂連結(lost)、其餘背景
   （以 peak 總長度標準化，並用 bootstrap 給信賴區間）
3. 列出被最多病人重複打中的 peak（復發性非編碼突變熱點）
"""
import os, sys
import numpy as np, pandas as pd

BASE="/raid/jupyter-alz.ee09/Tonsil_to_MCL"; RES=f"{BASE}/results"
def log(*a): print("[54]", *a, flush=True)

mut = pd.read_csv(f"{BASE}/data/external_mcl/wgs/MCL.somaticMutations.tsv", sep="\t")
mut.columns = ["chrom","start","end","ref","alt","sample"]
mut["patient"] = mut["sample"].str.replace(r"[PR]\d*$","",regex=True)
log(f"突變 {len(mut):,}，樣本 {mut['sample'].nunique()}，病人 {mut.patient.nunique()}")

pk = pd.read_csv(f"{BASE}/data/processed/atlas_peaks_sorted.csv")
M  = pd.read_csv(f"{RES}/52_rewiring_annotated.csv.gz")
gained = set(M.loc[M.gained,"peak_id"]); lost = set(M.loc[M.lost,"peak_id"])
tested = set(M.peak_id)
log(f"peak 總數 {len(pk):,}；被檢定過 {len(tested):,}；新生 {len(gained):,}；斷裂 {len(lost):,}")

# 把突變對到 peak
hits=[]
for c, g in pk.groupby("chrom"):
    mg = mut[mut.chrom==c]
    if len(mg)==0: continue
    g = g.sort_values("start")
    s=g.start.values; e=g.end.values; pid=g.pid.values
    j = np.searchsorted(s, mg.start.values, side="right")-1
    ok = (j>=0) & (mg.start.values < e[np.clip(j,0,len(e)-1)])
    if ok.any():
        hits.append(pd.DataFrame(dict(peak_id=pid[j[ok]], patient=mg.patient.values[ok],
                                      sample=mg["sample"].values[ok], pos=mg.start.values[ok],
                                      ref=mg.ref.values[ok], alt=mg.alt.values[ok], chrom=c)))
H = pd.concat(hits, ignore_index=True)
log(f"落在 peak 內的突變 {len(H):,}（{100*len(H)/len(mut):.1f}%）")
H.to_csv(f"{RES}/54_mutations_in_peaks.csv.gz", index=False, compression="gzip")

pk["len"] = pk.end - pk.start
plen = pk.set_index("pid")["len"]
cnt  = H.peak_id.value_counts()
def dens(ids, label):
    ids = [i for i in ids if i in plen.index]
    L = plen.loc[ids].sum(); n = cnt.reindex(ids).fillna(0).sum()
    npat = H[H.peak_id.isin(ids)].groupby("peak_id").patient.nunique()
    return dict(group=label, n_peaks=len(ids), total_bp=int(L), n_mut=int(n),
                mut_per_Mb=round(1e6*n/max(L,1),2),
                peaks_with_mut=int((cnt.reindex(ids).fillna(0)>0).sum()),
                peaks_multi_patient=int((npat>=2).sum()))
bg = tested - gained - lost
rows=[dens(gained,"新生連結 peak"), dens(lost,"斷裂連結 peak"),
      dens(bg,"其餘被檢定的 peak"), dens(set(pk.pid)-tested,"未被檢定的 peak")]
R=pd.DataFrame(rows)
log("\n=== 各組 peak 的體細胞突變密度 ===")
log("\n"+R.to_string(index=False))
r0 = R.loc[R.group=="其餘被檢定的 peak","mut_per_Mb"].iloc[0]
for _,x in R.iterrows():
    log(f"  {x.group}: 相對背景 {x.mut_per_Mb/r0:.2f}×")

log("\n=== 被最多『病人』重複打中的 peak（復發性非編碼突變熱點）Top 30 ===")
npat = H.groupby("peak_id").agg(n_patients=("patient","nunique"), n_mut=("patient","size")).reset_index()
npat = npat.merge(pk[["pid","chrom","start","end"]], left_on="peak_id", right_on="pid")
gm = M.groupby("peak_id").agg(genes=("gene", lambda s: ",".join(sorted(set(s))[:3])),
                              any_gained=("gained","max"), any_lost=("lost","max")).reset_index()
npat = npat.merge(gm, on="peak_id", how="left")
top = npat.sort_values(["n_patients","n_mut"], ascending=False).head(30)
pd.set_option("display.width",230)
log("\n"+top[["chrom","start","end","n_patients","n_mut","genes","any_gained","any_lost"]].to_string(index=False))
npat.to_csv(f"{RES}/54_peak_mutation_recurrence.csv", index=False)

log("\n=== 抗氧化軸線元件上的突變 ===")
AX=["SLC7A11","SLC3A2","GCLM","GCLC","MAFG","MAFF","MAFK","BACH1","NFE2L2","OXR1","FTH1"]
ax_pk = set(M.loc[M.gene.isin(AX),"peak_id"])
sub = npat[npat.peak_id.isin(ax_pk)].sort_values("n_patients", ascending=False)
log("\n"+sub.head(25)[["chrom","start","end","n_patients","n_mut","genes","any_gained"]].to_string(index=False))
log("DONE")
