#!/usr/bin/env python
"""步驟 130：四種「扣掉分化」的做法，在同一組基因上互相比較。

為什麼要做
----------
步驟 15 的殘差是**字面上的逐基因算術相減**（log1p 空間）。它可解釋、對批次穩健，
但沒有 per-gene 不確定度，也沒有 null 分布。本步驟把同一個問題用另外三種做法重算，
再看四者是否指向同一群基因 —— 一致 = 結論穩；不一致 = 知道線性假設在哪裡壞掉。

四種做法（全部輸出 cells × genes 的殘差，單位一致：log1p(1e4 × 正規化表現)）
---------------------------------------------------------------------------
A  gene-space      期望 = Σ_s w_s · 正常狀態 pseudobulk        （步驟 15 的做法）
B  scVI 反事實     把腫瘤細胞在 latent space 投到正常流形上（k=15 NN），
                   再用**該細胞自己的批次與文庫大小**解碼回表現量
D  latent 算術     scGen 的做法：δ = mean(z_腫瘤) − mean(z_其正常鄰居)，
                   z_cf = z − δ，解碼 → 這才是字面上的「向量空間中減去」
C  ContrastiveVI   不做減法：直接學一個「共享 latent（＝分化）」和一個
                   「顯著 latent（＝只有腫瘤才有的變異）」，
                   期望 = 只用共享 latent 解碼

⚠️ 一個會把結論弄反的技術細節（第一版踩到了）
---------------------------------------------
殘差兩邊的「平滑程度」必須一致。
單一細胞的 log1p(稀疏計數) 因為 Jensen 不等式，系統性低於 log1p(平滑速率)。
所以若「實測」用原始稀疏計數、「期望」用模型解碼出的平滑速率，
每個表現基因都會拿到系統性的負殘差 —— 跟生物無關。
  A：兩邊都是 log1p(稀疏) 的平均  → 一致
  C：兩邊都是解碼速率            → 一致
  B、D：第一版「實測」用原始計數  → 不一致，已改成用模型對該細胞自己的重建
       （＝解碼它自己的 latent），這樣兩邊都是平滑的，偏差相消。
B、D 因此量的是「把細胞移到正常流形上，模型預測會變什麼」，
這正是反事實該有的定義。

共同的雜訊底
------------
同一管檢體裡的非腫瘤 B 細胞走完全相同的流程。它們的殘差就是經驗 null。
每個方法各自算自己的訊雜比，所以方法之間可以比較「誰把訊號從雜訊裡拉得更開」。

輸出 results/130_*
"""
import os
import sys
import json

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _compat  # noqa

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import scipy.sparse as sp
import torch
from scipy.stats import spearmanr

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
PROC = f"{BASE}/data/processed"
RES = f"{BASE}/results"
N_HVG = 4000
N_LATENT = 20
K = 15
PANEL = ["SLC7A11", "SLC3A2", "GCLM", "GCLC", "AIFM2", "MAFG", "MAFF", "MAFK",
         "BACH1", "NFE2L2", "CCND1", "SOX11", "TXN", "GPX4", "HMOX1", "NQO1"]
import argparse

_ap = argparse.ArgumentParser()
_ap.add_argument("--seed", type=int, default=0,
                 help="torch / numpy / scvi 的隨機種子。0 = 已發表的那一輪")
_ap.add_argument("--tag", default="",
                 help="輸出檔名後綴。空字串 = 覆寫已發表的輸出，請在重跑時務必帶 tag")
_ap.add_argument("--group-by-donor", action="store_true",
                 help="把非腫瘤比較組按捐贈者拆開，額外輸出 per-donor 的欄位。"
                      "預設關閉；開啟時 tumour_specific 仍用合併後的 nontumour 計算，"
                      "所以輸出是既有輸出的嚴格超集，合併欄位必須逐位元重現")
_args, _ = _ap.parse_known_args()
SEED = _args.seed
TAG = _args.tag
SPLIT_DONOR = _args.group_by_donor

IG_PREFIX = ("IGHV", "IGKV", "IGLV", "IGHJ", "IGKJ", "IGLJ", "IGLC", "IGKC",
             "IGHA", "IGHG", "IGHM", "IGHD", "IGHE", "TRAV", "TRBV", "TRGV", "TRDV",
             "TRAC", "TRBC")


def log(*a):
    print("[130]", *a, flush=True)


def dense(X):
    return X.toarray() if sp.issparse(X) else np.asarray(X)


def lognorm(M, target=1e4):
    """counts (cells x genes) → log1p(target-normalised)."""
    s = np.asarray(M.sum(1)).ravel()
    s[s == 0] = 1.0
    return np.log1p(M / s[:, None] * target)


def rate_to_lognorm(R, target=1e4):
    """scVI 的 rate（每個細胞加總為 1）→ 同一尺度。"""
    s = R.sum(1, keepdims=True)
    s[s == 0] = 1.0
    return np.log1p(R / s * target)


# --------------------------------------------------------------- 資料 ---
def build():
    ref = sc.read_h5ad(f"{PROC}/ref_normalB_multiome.h5ad")
    ref.X = ref.layers["counts"].copy()
    ref.obs["source"] = "reference"
    ref.obs["batch"] = "ref_" + ref.obs["donor_id"].astype(str)
    ref.obs["state"] = ref.obs["annotation"].astype(str)
    ref = ref[:, ~ref.var_names.duplicated()].copy()

    s15 = pd.read_csv(f"{RES}/15_mcl_obs.csv", index_col=0)
    raw = sc.read_h5ad(f"{PROC}/mcl_gex_raw.h5ad")
    keep = [c for c in s15.index if c in set(raw.obs_names)]
    log("MCL 原始 counts 對得上的細胞:", len(keep), "/", len(s15))
    mcl = raw[keep].copy()
    del raw
    mcl.obs["source"] = "mcl"
    mcl.obs["batch"] = "mcl_" + mcl.obs["donor_id"].astype(str)
    mcl.obs["is_tumor"] = s15.reindex(mcl.obs_names)["is_tumor"].astype(str).values
    mcl.obs["state"] = s15.reindex(mcl.obs_names)["matched_normal_state"].astype(str).values
    mcl = mcl[:, ~mcl.var_names.duplicated()].copy()

    common = ref.var_names.intersection(mcl.var_names)
    ref, mcl = ref[:, common].copy(), mcl[:, common].copy()
    log("共同基因:", len(common))

    # 技術性黑名單（Ig 可變區／粒線體／核糖體）—— 不排除已知 driver，
    # 因為生成模型的工作本來就是把腫瘤變異放進 salient latent
    bl = set()
    for g in common:
        gu = g.upper()
        if gu.startswith(IG_PREFIX) or gu.startswith("MT-") or \
           gu.startswith("RPL") or gu.startswith("RPS"):
            bl.add(g)
    log("技術黑名單:", len(bl))

    tmp = ref.copy()
    sc.pp.highly_variable_genes(tmp, n_top_genes=N_HVG + len(bl) + 500,
                                flavor="seurat_v3", batch_key="donor_id")
    hv = [g for g in tmp.var_names[tmp.var.highly_variable] if g not in bl][:N_HVG]
    del tmp
    G = list(dict.fromkeys(hv + [g for g in PANEL if g in set(common)]))
    log("建模基因:", len(G), "| panel 全在內:",
        all(g in set(G) for g in PANEL if g in set(common)))

    A = ad.concat([ref[:, G], mcl[:, G]], label=None, index_unique=None)
    A.obs["source"] = list(ref.obs["source"]) + list(mcl.obs["source"])
    A.obs["batch"] = list(ref.obs["batch"]) + list(mcl.obs["batch"])
    A.obs["state"] = list(ref.obs["state"]) + list(mcl.obs["state"])
    A.obs["is_tumor"] = ["ref"] * ref.n_obs + list(mcl.obs["is_tumor"])
    A.X = sp.csr_matrix(A.X)
    log("合併:", A.shape, A.obs.source.value_counts().to_dict())
    return A, G


# ------------------------------------------------- 分組平均的累加器 ---
class GroupMean:
    """分塊累加 per-gene 平均，避免存整個 cells×genes 矩陣。"""

    def __init__(self, groups, n_genes):
        self.g = {k: (np.zeros(n_genes, np.float64), np.zeros(n_genes, np.float64), 0)
                  for k in groups}

    def add(self, key, M):
        s, s2, n = self.g[key]
        self.g[key] = (s + M.sum(0), s2 + (M ** 2).sum(0), n + M.shape[0])

    def mean(self, key):
        s, s2, n = self.g[key]
        return s / max(n, 1)

    def sd(self, key):
        s, s2, n = self.g[key]
        m = s / max(n, 1)
        return np.sqrt(np.maximum(s2 / max(n, 1) - m ** 2, 0))

    def n(self, key):
        return self.g[key][2]


def summarise(name, resid_fn, A, G, idx_groups, chunk=4096):
    """resid_fn(rows) -> (observed_lognorm, expected_lognorm) for those rows."""
    gm = GroupMean(list(idx_groups), len(G))
    for key, idx in idx_groups.items():
        for a in range(0, len(idx), chunk):
            rows = idx[a:a + chunk]
            obs, exp = resid_fn(rows)
            gm.add(key, obs - exp)
    out = pd.DataFrame({"gene": G})
    for key in idx_groups:
        out[f"mean_{key}"] = gm.mean(key)
        out[f"sd_{key}"] = gm.sd(key)
        out[f"n_{key}"] = gm.n(key)
    out["tumour_specific"] = out["mean_tumour"] - out["mean_nontumour"]
    out["method"] = name
    return out


def main():
    torch.manual_seed(SEED)
    np.random.seed(SEED)
    log("seed = %d | tag = %r" % (SEED, TAG))
    A, G = build()
    gi = {g: i for i, g in enumerate(G)}

    is_ref = (A.obs.source == "reference").values
    i_ref = np.where(is_ref)[0]
    i_mcl = np.where(~is_ref)[0]
    i_tum = np.where((A.obs.is_tumor == "True").values)[0]
    i_non = np.where((A.obs.is_tumor == "False").values)[0]
    log("ref %d | mcl %d (tumour %d, non-tumour %d)"
        % (len(i_ref), len(i_mcl), len(i_tum), len(i_non)))
    groups = {"tumour": i_tum, "nontumour": i_non}
    if SPLIT_DONOR:
        # 為什麼要這一段：非腫瘤比較組只有 819 個細胞，其中 790 個（96.5%）來自
        # 單一捐贈者 M413。tumour_specific 是「腫瘤 − 非腫瘤」，所以那個被減數
        # 實際上等於一位捐贈者。拆開來才能量「換掉捐贈者後基因排名會不會變」。
        # 合併的 nontumour 保留不動 → tumour_specific 的定義不變、輸出可逐位元對照。
        bt = A.obs["batch"].to_numpy()
        for dn in sorted(set(bt[i_non])):
            sub = i_non[bt[i_non] == dn]
            if len(sub) >= 5:
                groups["nontumour_" + str(dn)] = sub
        log("非腫瘤按捐贈者拆開:",
            {k: int(len(v)) for k, v in groups.items() if k.startswith("nontumour")})

    X = A.X
    OBS_CACHE = {}

    def observed(rows):
        key = (rows[0], len(rows))
        if key not in OBS_CACHE:
            OBS_CACHE.clear()
            OBS_CACHE[key] = lognorm(dense(X[rows]))
        return OBS_CACHE[key]

    results = []

    # ============================================================ A ===
    log("\n=== A：gene-space 減法（步驟 15 的做法） ===")
    W = pd.read_csv(f"{RES}/15_match_weights.csv.gz", index_col=0)
    names = A.obs_names.to_numpy()
    Wal = W.reindex(names[i_mcl])
    cov = float(Wal.notna().all(1).mean())
    log("match weights 對得上的比例: %.3f" % cov)
    # 正常狀態 pseudobulk（在 G 上，每捐贈者先平均）
    refl = lognorm(dense(X[i_ref]))
    rdf = A.obs.iloc[i_ref]
    rows, rn = [], []
    for st, sub in rdf.groupby("state", observed=True):
        per = []
        for dn, s2 in sub.groupby("batch", observed=True):
            if len(s2) < 10:
                continue
            per.append(refl[rdf.index.get_indexer(s2.index)].mean(0))
        if not per and len(sub) >= 10:
            per = [refl[rdf.index.get_indexer(sub.index)].mean(0)]
        if per:
            rows.append(np.mean(per, 0)); rn.append(st)
    P = pd.DataFrame(np.vstack(rows), index=rn, columns=G)
    log("正常狀態 pseudobulk:", P.shape)
    del refl
    Wm = Wal[[c for c in W.columns if c in P.index]].fillna(0).to_numpy(np.float32)
    Pm = P.loc[[c for c in W.columns if c in P.index]].to_numpy(np.float32)
    pos = {c: k for k, c in enumerate(names[i_mcl])}

    def resid_A(rows):
        obs = observed(rows)
        w = Wm[[pos[n] for n in names[rows]]]
        w = w / np.maximum(w.sum(1, keepdims=True), 1e-9)
        return obs, (w @ Pm)

    results.append(summarise("A_gene_space", resid_A, A, G, groups))
    log("A 完成")

    # ============================================================ scVI ===
    import scvi
    scvi.settings.seed = SEED
    log("\n=== 訓練 scVI ===")
    B = A.copy()
    scvi.model.SCVI.setup_anndata(B, batch_key="batch")
    m = scvi.model.SCVI(B, n_latent=N_LATENT, n_layers=2, gene_likelihood="nb")
    m.train(max_epochs=200, early_stopping=True, batch_size=512,
            accelerator="gpu", devices=1)
    m.save(f"{RES}/130_scvi_model{TAG}", overwrite=True, save_anndata=False)
    Z = m.get_latent_representation()
    log("latent:", Z.shape)
    np.save(f"{RES}/130_scvi_latent{TAG}.npy", Z)

    mod = m.module.eval()
    dev = next(mod.parameters()).device
    bkey = m.adata_manager.get_state_registry("batch").categorical_mapping
    bidx = pd.Categorical(A.obs["batch"], categories=list(bkey)).codes.astype(np.int64)
    lib = np.log(np.asarray(X.sum(1)).ravel().clip(1))

    def decode(z, rows):
        with torch.no_grad():
            zt = torch.tensor(z, dtype=torch.float32, device=dev)
            lt = torch.tensor(lib[rows], dtype=torch.float32, device=dev).reshape(-1, 1)
            bt = torch.tensor(bidx[rows], device=dev).reshape(-1, 1)
            g = mod.generative(z=zt, library=lt, batch_index=bt)
            return g["px"].scale.cpu().numpy()

    # ---- B：投到正常流形（k=15 NN） ----
    log("\n=== B：scVI 反事實（投到正常流形） ===")
    Zr = Z[i_ref].astype(np.float32)
    rn2 = (Zr ** 2).sum(1)

    def nn_ref(zq, k=K):
        d2 = rn2[None, :] - 2.0 * (zq @ Zr.T) + (zq ** 2).sum(1)[:, None]
        pt = np.argpartition(d2, k, axis=1)[:, :k]
        return Zr[pt].mean(1)

    def resid_B(rows):
        z = Z[rows].astype(np.float32)
        zcf = nn_ref(z)
        # 兩邊都走同一個解碼器 → Jensen 偏差相消
        return rate_to_lognorm(decode(z, rows)), rate_to_lognorm(decode(zcf, rows))

    results.append(summarise("B_scvi_counterfactual", resid_B, A, G, groups))
    log("B 完成")

    # ---- D：latent 算術（scGen） ----
    log("\n=== D：latent space 向量減法（scGen 的做法） ===")
    delta = Z[i_tum].mean(0) - nn_ref(Z[i_tum].astype(np.float32)).mean(0)
    log("δ 長度 = %.3f（latent 20 維）" % float(np.linalg.norm(delta)))
    json.dump({"delta": delta.tolist()}, open(f"{RES}/130_latent_delta{TAG}.json", "w"))

    def resid_D(rows):
        z = Z[rows].astype(np.float32)
        return rate_to_lognorm(decode(z, rows)), rate_to_lognorm(decode(z - delta, rows))

    results.append(summarise("D_latent_arithmetic", resid_D, A, G, groups))
    log("D 完成")
    del m, mod
    torch.cuda.empty_cache()

    # ====================================================== ContrastiveVI ===
    log("\n=== C：ContrastiveVI（學出來，不是減出來） ===")
    from scvi.external import ContrastiveVI
    C = A.copy()
    ContrastiveVI.setup_anndata(C, batch_key="batch")
    cm = ContrastiveVI(C, n_background_latent=N_LATENT, n_salient_latent=10, n_layers=2)
    cm.train(background_indices=list(i_ref), target_indices=list(i_mcl),
             max_epochs=150, batch_size=512, accelerator="gpu", devices=1)
    log("ContrastiveVI 訓練完成")
    bg = cm.get_specific_normalized_expression(C, expression_type="background",
                                               library_size=1.0, return_numpy=True)
    fl = cm.get_salient_normalized_expression(C, library_size=1.0, return_numpy=True)
    BG, FL = rate_to_lognorm(np.asarray(bg)), rate_to_lognorm(np.asarray(fl))
    np.save(f"{RES}/130_cvi_salient_latent{TAG}.npy", cm.get_latent_representation(
        C, representation_kind="salient"))

    def resid_C(rows):
        return FL[rows], BG[rows]

    results.append(summarise("C_contrastivevi", resid_C, A, G, groups))
    log("C 完成")

    # ============================================================ 比較 ===
    R = pd.concat(results, ignore_index=True)
    R.to_csv(f"{RES}/130_residual_by_method{TAG}.csv.gz", index=False, compression="gzip")

    piv = R.pivot(index="gene", columns="method", values="tumour_specific")
    meth = list(piv.columns)
    rho = pd.DataFrame(index=meth, columns=meth, dtype=float)
    jac = pd.DataFrame(index=meth, columns=meth, dtype=float)
    tops = {m: set(piv[m].sort_values(ascending=False).head(100).index) for m in meth}
    for a in meth:
        for b in meth:
            rho.loc[a, b] = spearmanr(piv[a], piv[b]).statistic
            jac.loc[a, b] = len(tops[a] & tops[b]) / len(tops[a] | tops[b])
    log("\n=== 四法之間的 Spearman（tumour-specific 殘差） ===\n", rho.round(3).to_string())
    log("\n=== top-100 基因的 Jaccard 重疊 ===\n", jac.round(3).to_string())

    snr = {}
    for m_ in meth:
        d = R[R.method == m_].set_index("gene")
        floor = d["mean_nontumour"].abs().mean()
        snr[m_] = {
            "noise_floor_mean_abs_nontumour": float(floor),
            "tumour_mean_abs": float(d["mean_tumour"].abs().mean()),
            "global_ratio": float(d["mean_tumour"].abs().mean() / max(floor, 1e-9)),
            "panel_z_over_floor": {g: float(d.loc[g, "tumour_specific"] / max(floor, 1e-9))
                                   for g in PANEL if g in d.index},
            "panel_rank": {g: int((piv[m_] > piv[m_].get(g, -9e9)).sum() + 1)
                           for g in PANEL if g in piv.index},
        }
    log("\n=== 每個方法的訊雜比 ===")
    for m_ in meth:
        log("  %-24s 雜訊底 %.4f | 腫瘤 %.4f | 比值 %.2f"
            % (m_, snr[m_]["noise_floor_mean_abs_nontumour"],
               snr[m_]["tumour_mean_abs"], snr[m_]["global_ratio"]))
    log("\n=== panel 基因在各方法的排名（共 %d 個基因） ===" % len(piv))
    hdr = "%-9s" % "gene" + "".join("%26s" % m_ for m_ in meth)
    log(hdr)
    for g in PANEL:
        if g not in piv.index:
            continue
        log("%-9s" % g + "".join("%26s" % ("#%d  (%+.3f)" % (snr[m_]["panel_rank"][g],
                                                             piv.loc[g, m_]))
                                 for m_ in meth))

    json.dump({"spearman": json.loads(rho.to_json()),
               "jaccard_top100": json.loads(jac.to_json()),
               "snr": snr, "n_genes": int(len(piv)),
               "n_tumour": int(len(i_tum)), "n_nontumour": int(len(i_non)),
               "n_reference": int(len(i_ref))},
              open(f"{RES}/130_method_comparison{TAG}.json", "w"), indent=2, ensure_ascii=False)
    piv.to_csv(f"{RES}/130_tumour_specific_by_method{TAG}.csv")
    log("\nDONE → results/130_*")


if __name__ == "__main__":
    main()
