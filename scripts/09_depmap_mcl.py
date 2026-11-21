#!/usr/bin/env python
"""
步驟 09：DepMap 中 MCL 細胞株的「選擇性依賴」分析。

輸出
----
results/09_depmap_mcl_dependency.csv   每個基因在 MCL 株 vs 其他株的 gene effect 差異
results/09_depmap_mcl_expression.csv   MCL 株的表現量
results/09_depmap_mcl_mutations.csv    MCL 株的體細胞突變
"""
import os, sys
import numpy as np
import pandas as pd

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
D = f"{BASE}/data/depmap"
R = f"{BASE}/results"
os.makedirs(R, exist_ok=True)


def log(*a):
    print("[09]", *a, flush=True)


model = pd.read_csv(f"{D}/Model.csv", low_memory=False)
mcl = model[model["OncotreeCode"].astype(str) == "MCL"]
log("MCL 細胞株:", len(mcl))
log(mcl[["ModelID", "CellLineName", "OncotreeSubtype", "OncotreeLineage"]].to_string(index=False))
mcl_ids = set(mcl["ModelID"])

# 近緣對照：其他成熟 B 細胞腫瘤（DLBCL/BL/CLL/FL/MM）
lymphoid = model[model["OncotreeLineage"].astype(str) == "Lymphoid"]
bcell_codes = {"DLBCLNOS", "BL", "CLLSLL", "FL", "MM", "PCM", "MZL", "LPL", "HGBCL", "PBL"}
bcell_other = set(lymphoid[lymphoid["OncotreeCode"].astype(str).isin(bcell_codes)]["ModelID"]) - mcl_ids
log("其他 B 系腫瘤株:", len(bcell_other))

# ---------------- CRISPR gene effect ----------------
log("讀 CRISPRGeneEffect.csv ...")
ge = pd.read_csv(f"{D}/CRISPRGeneEffect.csv", index_col=0, low_memory=False)
log("shape", ge.shape)
ge.columns = [c.split(" (")[0] for c in ge.columns]

have_mcl = sorted(mcl_ids & set(ge.index))
have_b = sorted(bcell_other & set(ge.index))
have_all = [i for i in ge.index if i not in mcl_ids]
log("有 CRISPR 的 MCL 株:", have_mcl)
log("有 CRISPR 的其他 B 系株:", len(have_b), " | 全部非-MCL 株:", len(have_all))

if len(have_mcl) == 0:
    log("警告：沒有任何 MCL 株有 CRISPR 資料，跳過依賴性分析")
    sys.exit(0)

m_mcl = ge.loc[have_mcl]
m_oth = ge.loc[have_all]
m_b = ge.loc[have_b] if have_b else None

res = pd.DataFrame({
    "gene": ge.columns,
    "mcl_mean": m_mcl.mean(axis=0).values,
    "mcl_min": m_mcl.min(axis=0).values,
    "mcl_max": m_mcl.max(axis=0).values,
    "other_mean": m_oth.mean(axis=0).values,
    "other_sd": m_oth.std(axis=0).values,
})
if m_b is not None:
    res["bcell_other_mean"] = m_b.mean(axis=0).values
for i, mid in enumerate(have_mcl):
    nm = model.set_index("ModelID").loc[mid, "CellLineName"]
    res[f"{nm}"] = m_mcl.loc[mid].values

res["delta_vs_all"] = res["mcl_mean"] - res["other_mean"]
res["z_vs_all"] = res["delta_vs_all"] / res["other_sd"].replace(0, np.nan)
if m_b is not None:
    res["delta_vs_bcell"] = res["mcl_mean"] - res["bcell_other_mean"]

# 「選擇性必需」= MCL 本身很必需（<-0.5）且比其他株更必需
res["selective"] = (res["mcl_max"] < -0.3) & (res["delta_vs_all"] < -0.2)
res = res.sort_values("z_vs_all")
res.to_csv(f"{R}/09_depmap_mcl_dependency.csv", index=False)
log("寫出", f"{R}/09_depmap_mcl_dependency.csv")

log("\n=== MCL 選擇性依賴 Top 30 (z_vs_all 最負) ===")
cols = ["gene", "mcl_mean", "mcl_max", "other_mean", "delta_vs_all", "z_vs_all"]
print(res.head(30)[cols].to_string(index=False))

log("\n=== MCL 全株皆強必需 (mcl_max < -0.8) 且非泛必需 (other_mean > -0.4) ===")
sel = res[(res.mcl_max < -0.8) & (res.other_mean > -0.4)].sort_values("delta_vs_all")
print(sel.head(40)[cols].to_string(index=False))
sel.to_csv(f"{R}/09_depmap_mcl_selective_essential.csv", index=False)

# ---------------- 表現量 ----------------
log("\n讀表現量 ...")
ex = pd.read_csv(f"{D}/OmicsExpressionProteinCodingGenesTPMLogp1.csv", index_col=0, low_memory=False)
ex.columns = [c.split(" (")[0] for c in ex.columns]
ids = [i for i in mcl_ids if i in ex.index]
sub = ex.loc[ids]
sub.index = [model.set_index("ModelID").loc[i, "CellLineName"] for i in ids]
sub.T.to_csv(f"{R}/09_depmap_mcl_expression.csv")
log("MCL 株表現量 shape:", sub.shape)
for g in ["CCND1", "SOX11", "CDK4", "CDK6", "BCL2", "MYC", "ATM", "TP53", "BTK", "MCL1"]:
    if g in sub.columns:
        log(f"  {g:7s}", " ".join(f"{n}={v:.2f}" for n, v in sub[g].items()))

# ---------------- 突變 ----------------
log("\n讀突變 ...")
mu = pd.read_csv(f"{D}/OmicsSomaticMutations.csv", low_memory=False)
mm = mu[mu["ModelID"].isin(mcl_ids)]
mm.to_csv(f"{R}/09_depmap_mcl_mutations.csv", index=False)
log("MCL 株突變筆數:", len(mm))
if len(mm):
    top = mm[mm.get("LikelyLoF", False) == True]["HugoSymbol"].value_counts().head(20) \
        if "LikelyLoF" in mm.columns else mm["HugoSymbol"].value_counts().head(20)
    log("常見突變基因:\n", top.to_string())
log("DONE")
