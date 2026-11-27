#!/usr/bin/env python
"""步驟 233：步驟 120 的狀態指派，準確率是多少？「腫瘤很集中」是生物還是平滑假象？

為什麼需要這一支
----------------
步驟 120 報的 92.4% 是 **kNN 與 softmax 兩個方法彼此一致**，不是準確率。
兩者用同一個參考系、同一組對位基因、同一個鄰域結構 —— **共同偏誤會表現為高一致率**。
整個管線從來沒有量過：把一顆**已知狀態**的細胞丟進去，它會不會回到自己的狀態。

而且從既有檔案就能看出一個更尖的問題。腫瘤的投影 dpt 是
**15 個鄰居的中位數**，所以天生被平滑：
    腫瘤投影 dpt      IQR 0.0431–0.0515（寬 0.0084）
    參考真 ncsMBC dpt IQR 0.0281–0.0497（寬 0.0215）→ 窄了 2.6 倍
我上一輪把「腫瘤在分化軸上幾乎沒有展幅」當成生物發現來講，
但它有多少是 kNN 中位數造成的，**當時沒有量**。這一支就是去量它。
（另外：腫瘤的 dpt 中位 0.0481 其實**落在參考 ncsMBC 自己的 IQR 之內**，
所以我上一輪說的「+0.0096 位移」本來就不該叫位移 —— 這一點已經可以直接更正。）

方法：對參考細胞做留出法，走**完全相同**的管線
-----------------------------------------------
每個狀態分層抽出 15%（至少 10 顆）當留出集，剩下 85% 擬合 PCA(50)，
再用同一套分塊 kNN（K=15）把留出細胞投影回去，取鄰居的狀態多數票與 dpt 中位數。
留出細胞的真實狀態與真實 dpt 已知，所以準確率與平滑程度都可以直接量。

事先寫死的判定（在跑之前寫完）
------------------------------
C0 **健全性（守門）**：留出細胞的「投影 dpt」與「真實 dpt」Spearman > 0.80。
   不過 → 投影本身壞了，整支作廢。
C1 **狀態回收率**：報每個狀態的回收率。腫瘤被指派的那個狀態 `ncsMBC`
   必須 ≥ **0.80**，「腫瘤是 ncsMBC」才撐得起它目前的說法。
   事先預測：ncsMBC ≥ 0.85（它有 6,788 顆，是最大的狀態之一）。
C2 **平滑假象的大小（本支的主要問題）**：留出 `ncsMBC` 細胞的**投影** dpt 的 IQR 寬度。
   規則：若該寬度 ≤ 腫瘤 IQR 寬度（0.0084）的 **2 倍**，
   則「腫瘤在分化軸上特別集中」**必須撤回** —— 因為正常細胞走同一套管線也一樣集中。
   事先預測：會 ≤ 2 倍，也就是**我上一輪的「緊」那個重構站不住**。
C3 **中位數位移**：留出 `ncsMBC` 的投影 dpt 中位數 − 真實 dpt 中位數。
   若它與腫瘤的 +0.0096 同向且量級相當（0.5–2 倍之間），則該位移判定為管線假象。
C4 **鄰居距離校準（描述性，不當門檻）**：留出細胞的 knn_dist_mean，
   對照腫瘤 15.78 與同檢體非腫瘤 15.19。
   **刻意不當門檻**：留出細胞與參考系同批次，距離一定偏小，這個比較被批次混淆。
   腫瘤真正公平的對照是同檢體非腫瘤細胞，而那個數字已經有了。

輸出 results/233_projection_control.json
"""
import os
import sys
import json

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "8"

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _compat  # noqa

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.stats import spearmanr
from sklearn.decomposition import PCA

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
PROC = f"{BASE}/data/processed"
RES = f"{BASE}/results"
K = 15
HOLDOUT_FRAC = 0.15
MIN_PER_STATE = 10
SEED = 0
FOCUS = "ncsMBC"
C1_MIN = 0.80
C2_MULT = 2.0
TUMOUR_IQR_W = 0.0084          # 來自 results/120_mcl_cells.csv.gz，寫死以便判定可重現
TUMOUR_DPT_MED = 0.0481
TUMOUR_KNN_DIST = 15.783
NONTUM_KNN_DIST = 15.190


def log(*a):
    print("[233]", *a, flush=True)


def dense(X):
    return X.toarray() if sp.issparse(X) else np.asarray(X)


def zrow(M):
    M = np.asarray(M, dtype=np.float32)
    return (M - M.mean(1, keepdims=True)) / (M.std(1, keepdims=True) + 1e-8)


def knn_chunked(Q, R, k):
    rn2 = (R ** 2).sum(1)
    ind = np.empty((Q.shape[0], k), dtype=np.int64)
    dist = np.empty((Q.shape[0], k), dtype=np.float32)
    CH = 1024
    for a in range(0, Q.shape[0], CH):
        b = min(a + CH, Q.shape[0])
        d2 = rn2[None, :] - 2.0 * (Q[a:b] @ R.T) + (Q[a:b] ** 2).sum(1)[:, None]
        pt = np.argpartition(d2, k, axis=1)[:, :k]
        dd = np.take_along_axis(d2, pt, axis=1)
        o = np.argsort(dd, axis=1)
        ind[a:b] = np.take_along_axis(pt, o, axis=1)
        dist[a:b] = np.sqrt(np.maximum(np.take_along_axis(dd, o, axis=1), 0))
    return ind, dist


def main():
    ref = sc.read_h5ad(f"{PROC}/ref_normalB_multiome.h5ad")
    log("參考:", ref.shape)
    sc.pp.normalize_total(ref, target_sum=1e4)
    sc.pp.log1p(ref)
    hv = json.load(open(f"{RES}/15_matching_genes.json"))
    hv = [g for g in hv if g in set(ref.var_names)]
    log("對位基因:", len(hv))

    cdf = pd.read_csv(f"{RES}/120_ref_cells.csv.gz").set_index("cell").reindex(ref.obs_names)
    if not cdf["dpt"].notna().all():
        raise SystemExit("120_ref_cells.csv.gz 的 dpt 不完整，先跑步驟 120")
    dpt_true = cdf["dpt"].to_numpy(float)
    state_true = cdf["state"].astype(str).to_numpy()

    # ---------- 分層留出 ----------
    rng = np.random.default_rng(SEED)
    held = np.zeros(ref.n_obs, dtype=bool)
    for st in np.unique(state_true):
        idx = np.where(state_true == st)[0]
        n_h = max(min(MIN_PER_STATE, len(idx) // 2), int(round(len(idx) * HOLDOUT_FRAC)))
        if len(idx) < 4:
            continue
        held[rng.choice(idx, size=min(n_h, len(idx) - 2), replace=False)] = True
    log("留出 %d / %d 顆（%d 個狀態）"
        % (held.sum(), ref.n_obs, len(np.unique(state_true[held]))))

    Xz = zrow(dense(ref[:, hv].X))
    pca = PCA(n_components=50, random_state=0).fit(Xz[~held])
    Rp = pca.transform(Xz[~held]).astype(np.float32)
    Qp = pca.transform(Xz[held]).astype(np.float32)
    log("PCA 擬合於保留的 %d 顆，解釋變異 %.3f"
        % ((~held).sum(), float(pca.explained_variance_ratio_.sum())))
    del Xz

    ind, dist = knn_chunked(Qp, Rp, K)
    log("kNN 完成")

    ann_keep = state_true[~held]
    dpt_keep = dpt_true[~held]
    nb = ann_keep[ind]
    vote, votefrac = [], []
    for row in nb:
        v, c = np.unique(row, return_counts=True)
        j = int(np.argmax(c))
        vote.append(v[j]); votefrac.append(c[j] / K)
    H = pd.DataFrame({
        "true_state": state_true[held],
        "pred_state": vote,
        "vote_frac": votefrac,
        "dpt_true": dpt_true[held],
        "dpt_proj": np.median(dpt_keep[ind], axis=1),
        "knn_dist_mean": dist.mean(1)})
    H.to_csv(f"{RES}/233_holdout_cells.csv.gz", index=False, compression="gzip")

    out = {"n_reference": int(ref.n_obs), "n_holdout": int(held.sum()),
           "k": K, "holdout_frac": HOLDOUT_FRAC, "seed": SEED,
           "pca_var_explained": float(pca.explained_variance_ratio_.sum())}

    # ---------------------------------------------------------- C0 ---
    rs = spearmanr(H.dpt_true, H.dpt_proj)
    c0 = bool(rs.statistic > 0.80)
    out["C0_projection_sane"] = {"rule": "投影 dpt 與真實 dpt 的 Spearman > 0.80",
                                 "spearman_r": float(rs.statistic),
                                 "pvalue": float(rs.pvalue), "pass": c0}
    log("\n=== C0 健全性 ===")
    log("  投影 dpt vs 真實 dpt：Spearman r = %.4f (p = %.3g) → %s"
        % (rs.statistic, rs.pvalue, "通過" if c0 else "未通過，整支作廢"))
    if not c0:
        out["verdict"] = "VOID_projection_broken"
        json.dump(out, open(f"{RES}/233_projection_control.json", "w"),
                  indent=2, ensure_ascii=False)
        return 1

    # ---------------------------------------------------------- C1 ---
    H["correct"] = H.true_state == H.pred_state
    acc = H.correct.mean()
    _rows = []
    for st, g in H.groupby("true_state"):
        wrong = g.loc[~g.correct, "pred_state"]
        _rows.append({"true_state": st, "n": len(g),
                      "recovery": float(g.correct.mean()),
                      "top_confusion": (wrong.value_counts().index[0]
                                        if len(wrong) else "—")})
    per = pd.DataFrame(_rows).set_index("true_state").sort_values("recovery")
    per.to_csv(f"{RES}/233_state_recovery.csv")
    foc = float(per.loc[FOCUS, "recovery"]) if FOCUS in per.index else float("nan")
    c1 = bool(foc >= C1_MIN)
    out["C1_state_recovery"] = {
        "rule": "%s 的回收率 ≥ %.2f" % (FOCUS, C1_MIN),
        "overall_accuracy": float(acc),
        "focus_state": FOCUS, "focus_recovery": foc,
        "n_states": int(len(per)),
        "per_state": {str(k): {"n": int(v["n"]), "recovery": float(v["recovery"]),
                               "top_confusion": str(v["top_confusion"])}
                      for k, v in per.iterrows()},
        "pass": c1}
    log("\n=== C1 狀態回收率（留出細胞，走完全相同的管線）===")
    log("  整體準確率 %.3f（%d 個狀態）" % (acc, len(per)))
    log("  %s 回收率 %.3f（門檻 %.2f）→ %s" % (FOCUS, foc, C1_MIN, "通過" if c1 else "未通過"))
    log("\n  最差的 8 個狀態：")
    for k, v in per.head(8).iterrows():
        log("    %-34s n=%5d 回收 %.3f  最常誤判成 → %s"
            % (k[:34], int(v["n"]), v["recovery"], v["top_confusion"]))
    log("  最好的 5 個：")
    for k, v in per.tail(5).iterrows():
        log("    %-34s n=%5d 回收 %.3f" % (k[:34], int(v["n"]), v["recovery"]))

    # ---------------------------------------------------------- C2 ---
    hf = H[H.true_state == FOCUS]
    w_proj = float(hf.dpt_proj.quantile(.75) - hf.dpt_proj.quantile(.25))
    w_true = float(hf.dpt_true.quantile(.75) - hf.dpt_true.quantile(.25))
    c2_artefact = bool(w_proj <= C2_MULT * TUMOUR_IQR_W)
    out["C2_smoothing_artefact"] = {
        "rule": "留出 %s 的投影 dpt IQR 寬度 ≤ %.1f × 腫瘤的 %.4f "
                "→ 「腫瘤特別集中」必須撤回" % (FOCUS, C2_MULT, TUMOUR_IQR_W),
        "holdout_proj_iqr_width": w_proj,
        "holdout_true_iqr_width": w_true,
        "tumour_iqr_width": TUMOUR_IQR_W,
        "smoothing_factor_in_normal_cells": (w_true / w_proj) if w_proj > 0 else None,
        "tightness_claim_must_be_withdrawn": c2_artefact}
    log("\n=== C2 「腫瘤很集中」是生物還是平滑假象 ===")
    log("  留出 %s：真實 dpt IQR 寬 %.4f → 投影後 %.4f（被平滑 %.1f 倍）"
        % (FOCUS, w_true, w_proj, (w_true / w_proj) if w_proj > 0 else float("nan")))
    log("  腫瘤投影 dpt IQR 寬 %.4f" % TUMOUR_IQR_W)
    log("  → %s" % ("正常細胞走同一管線也一樣集中 → **「腫瘤特別集中」必須撤回**"
                    if c2_artefact else
                    "正常細胞投影後仍明顯較寬 → 腫瘤的集中程度有生物成分"))

    # ---------------------------------------------------------- C3 ---
    shift_h = float(hf.dpt_proj.median() - hf.dpt_true.median())
    shift_t = float(TUMOUR_DPT_MED - hf.dpt_true.median())
    ratio = (shift_h / shift_t) if shift_t != 0 else None
    c3_artefact = bool(ratio is not None and 0.5 <= ratio <= 2.0)
    out["C3_median_displacement"] = {
        "rule": "留出細胞的位移與腫瘤的位移同向且在 0.5–2 倍之間 → 判定為管線假象",
        "holdout_shift": shift_h, "tumour_shift": shift_t,
        "ratio": ratio, "displacement_is_artefact": c3_artefact}
    log("\n=== C3 中位數位移 ===")
    log("  留出 %s：真實中位 %.4f → 投影中位 %.4f（位移 %+.4f）"
        % (FOCUS, hf.dpt_true.median(), hf.dpt_proj.median(), shift_h))
    log("  腫瘤：投影中位 %.4f，對同一真實中位的位移 %+.4f" % (TUMOUR_DPT_MED, shift_t))
    log("  比值 %s → %s" % (("%.2f" % ratio) if ratio is not None else "—",
                            "位移是管線假象" if c3_artefact else "位移不能單以管線解釋"))

    # ---------------------------------------------------------- C4 ---
    out["C4_neighbour_distance"] = {
        "note": "描述性。留出細胞與參考系同批次，距離必然偏小，此比較被批次混淆；"
                "腫瘤公平的對照是同檢體非腫瘤細胞（15.190）。",
        "holdout_median": float(H.knn_dist_mean.median()),
        "holdout_q25_q75": [float(H.knn_dist_mean.quantile(.25)),
                            float(H.knn_dist_mean.quantile(.75))],
        "tumour_median": TUMOUR_KNN_DIST,
        "intra_sample_nontumour_median": NONTUM_KNN_DIST}
    log("\n=== C4 鄰居距離（描述性）===")
    log("  留出參考細胞 %.3f | 同檢體非腫瘤 %.3f | 腫瘤 %.3f"
        % (H.knn_dist_mean.median(), NONTUM_KNN_DIST, TUMOUR_KNN_DIST))
    log("  （留出細胞同批次，數字必然偏小，不可直接拿來說腫瘤是外插）")

    out["verdict"] = {"C0": c0, "C1_focus_pass": c1,
                      "C2_tightness_claim_withdrawn": c2_artefact,
                      "C3_displacement_is_artefact": c3_artefact}
    json.dump(out, open(f"{RES}/233_projection_control.json", "w"),
              indent=2, ensure_ascii=False)
    log("\n→ results/233_projection_control.json、233_state_recovery.csv、"
        "233_holdout_cells.csv.gz")
    return 0


if __name__ == "__main__":
    sys.exit(main())
