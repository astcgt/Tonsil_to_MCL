#!/usr/bin/env python
"""步驟 212：全基因組的「治療窗口」篩選 —— 把這個專案的邏輯一次做到底。

為什麼要做這一步
----------------
到今天為止，所有想用「選擇性結合」取得治療窗口的路線都關閉了：
  C5 `MAFG` 在任何細胞都不是遺傳依賴
  C6 設計品分不開 MafG / MafF / MafK（三個共摺疊模型 + MD 都同意）
  D5 夥伴選擇性在幾何上不存在（NRF2 自己就遮住 BACH1 足跡的 82.1%）
  211 元件層面的序列選擇性不存在（完整 ARE 比例倍率 0.990，經驗 p = 0.80）

但這個專案真正的工具從頭到尾不是結合劑，是**分化配對**。
窗口不必來自分子的選擇性，可以來自**細胞狀態的差異**：
一個基因如果（a）在腫瘤裡遠高於「分化程度相同的那顆正常 B 細胞」，
而且（b）是 B 淋巴系腫瘤特有的遺傳依賴 —— 那打它就天生有窗口，
不需要藥物自己會分辨。

這一步把那個邏輯對**全基因組**跑一次。它同時是一個方法驗證：
如果 `SLC7A11` 在沒有被告知的情況下自己浮到前面，那說明這條路徑不是
事後挑出來的；如果有別的基因排得更前面，那就是一個更好的標的。

事前指定
--------
分數 = min( z(殘差), z(B 淋巴系依賴選擇性) )，且兩項都必須 > 0
     （第一版用相加，被 P4 抓到缺陷：一項強就能把另一項為負的基因拖進前段。
       治療窗口的定義是「兩個條件同時成立」，所以取兩項的較小者。）
  殘差       = `results/15_gene_residual_tumor.csv` 的 mean_residual（腫瘤 − 分化配對期望）
  依賴選擇性 = mean(非淋巴系 gene effect) − mean(B 淋巴系 gene effect)
               （正值＝B 淋巴系更依賴）

P1  `SLC7A11` 必須落在前 1%。這是方法驗證：它是用完全獨立的證據鏈找到的，
    如果這個無監督排名把它丟到後面，那排名本身有問題。
P2  泛必需基因（兩組都強烈依賴）不該霸榜 —— 它們沒有選擇性，
    必須被依賴選擇性那一項自動壓下去。用核糖體基因當檢查。
P3  經驗虛無：把兩份資料之間的基因標籤打亂 2,000 次，
    `SLC7A11` 的名次必須明顯優於虛無分佈。
P4（陰性對照）在腫瘤裡**下降**的基因不該出現在前段。

⚠️ 限制：依賴性來自細胞株，殘差來自病人。這兩者是不同的系統，
   而且本專案已經證明細胞株是這條軸線的壞模型
   （三株 MCL 的 `SLC7A11` 都低於中位，病人卻高 200 倍）。
   所以這是一個**產生假說**的排名，不是一個效力預測。

輸出 results/212_*
"""
import os
import json

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "8"

import numpy as np
import pandas as pd
from scipy import stats

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
RES = f"{BASE}/results"
DM = f"{BASE}/data/depmap"
N_PERM = 2000
RNG = np.random.default_rng(0)
RIBO_PREFIX = ("RPL", "RPS")


def log(*a):
    print("[212]", *a, flush=True)


def main():
    R = pd.read_csv(f"{RES}/15_gene_residual_tumor.csv", index_col=0)
    log("殘差表 %d 個基因" % len(R))

    hdr = pd.read_csv(f"{DM}/CRISPRGeneEffect.csv", nrows=0)
    first = hdr.columns[0]
    sym2col = {}
    for c in hdr.columns[1:]:
        sym2col.setdefault(c.split(" (")[0], c)
    shared = [g for g in R.index if g in sym2col]
    log("與 DepMap 共有 %d 個基因，分批讀取 gene effect…" % len(shared))

    mod = pd.read_csv(f"{DM}/Model.csv").set_index("ModelID")
    lin = mod["OncotreeLineage"].astype(str)
    dis = mod["OncotreePrimaryDisease"].astype(str) if "OncotreePrimaryDisease" in mod.columns else lin

    cols = [sym2col[g] for g in shared]
    chunks, step = [], 2000
    for i in range(0, len(cols), step):
        sub = pd.read_csv(f"{DM}/CRISPRGeneEffect.csv",
                          usecols=[first] + cols[i:i + step]).set_index(first)
        chunks.append(sub)
    E = pd.concat(chunks, axis=1)
    E.columns = [c.split(" (")[0] for c in E.columns]
    log("gene effect：%d 株 × %d 基因" % E.shape)

    ln = lin.reindex(E.index).fillna("NA")
    dz = dis.reindex(E.index).fillna("NA")
    is_lymph = ln.str.contains("Lymph", case=False)
    # B 細胞腫瘤：淋巴系且疾病名稱不含 T-cell / NK
    is_B = is_lymph & ~dz.str.contains("T-Cell|T Cell|NK", case=False, regex=True)
    log("譜系：淋巴系 %d 株，其中判為 B 細胞腫瘤 %d 株，非淋巴系 %d 株"
        % (int(is_lymph.sum()), int(is_B.sum()), int((~is_lymph).sum())))

    depB = E[is_B.values].mean()
    depO = E[(~is_lymph).values].mean()
    sel = depO - depB                      # 正值＝B 細胞腫瘤更依賴

    df = pd.DataFrame({
        "residual": R["mean_residual"].reindex(sel.index),
        "observed": R["mean_observed"].reindex(sel.index),
        "expected": R["mean_expected"].reindex(sel.index),
        "t_like": R["t_like"].reindex(sel.index),
        "dep_Bcell": depB, "dep_other": depO, "dep_selectivity": sel,
    }).dropna(subset=["residual", "dep_selectivity"])
    log("可評分的基因 %d 個" % len(df))

    df["z_residual"] = stats.zscore(df.residual, nan_policy="omit")
    df["z_selectivity"] = stats.zscore(df.dep_selectivity, nan_policy="omit")
    # 相加的分數有一個缺陷，是 P4 自己抓到的：一項很強就能把另一項為負的基因
    # 拖進前段（NMNAT1 殘差 −0.031、UMPS −0.026、ADSL −0.025 都是這樣進來的）。
    # 治療窗口的定義是「兩項都成立」，不是「加起來很大」，所以改成閘控：
    # 必須殘差 > 0 且依賴選擇性 > 0，否則不計分。
    df["score_sum"] = df.z_residual + df.z_selectivity
    gate = (df.residual > 0) & (df.dep_selectivity > 0)
    df["score"] = np.where(gate, np.minimum(df.z_residual, df.z_selectivity),
                           -np.inf)
    df["gated"] = gate
    df["rank"] = df.score.rank(ascending=False).astype(int)
    df = df.sort_values("score", ascending=False)
    df.to_csv(f"{RES}/212_window_screen.csv")

    n = len(df)
    out = {"n_genes": int(n), "n_B_lines": int(is_B.sum()),
           "n_other_lines": int((~is_lymph).sum())}

    log("\n=== 前 20 名 ===")
    cols_show = ["rank", "observed", "expected", "residual",
                 "dep_Bcell", "dep_other", "dep_selectivity", "score",
                 "score_sum"]
    log("\n" + df.head(20)[cols_show].round(3).to_string())

    r_slc = int(df.loc["SLC7A11", "rank"]) if "SLC7A11" in df.index else -1
    out["P1_SLC7A11_rank"] = {"rank": r_slc, "n": int(n),
                              "percentile": r_slc / n if r_slc > 0 else None,
                              "pass": bool(0 < r_slc <= 0.01 * n)}
    log("\nP1 `SLC7A11` 名次 %d / %d（前 %.2f%%）→ %s"
        % (r_slc, n, 100 * r_slc / n if r_slc > 0 else -1,
           "過" if out["P1_SLC7A11_rank"]["pass"] else "不過"))

    ribo = [g for g in df.index if g.startswith(RIBO_PREFIX)]
    top1 = df.head(max(1, int(0.01 * n))).index
    n_ribo_top = len([g for g in ribo if g in top1])
    out["P2_pan_essentials_not_on_top"] = {
        "n_ribosomal_genes": len(ribo), "n_in_top_1pct": n_ribo_top,
        "median_rank_ribosomal": float(df.loc[ribo, "rank"].median()) if ribo else None,
        "pass": bool(n_ribo_top <= max(1, 0.01 * len(ribo)))}
    log("P2 核糖體基因 %d 個，中位名次 %.0f，前 1%% 裡有 %d 個 → %s"
        % (len(ribo), df.loc[ribo, "rank"].median() if ribo else -1, n_ribo_top,
           "過" if out["P2_pan_essentials_not_on_top"]["pass"] else "不過"))

    # P3：打亂兩份資料之間的基因對應
    zr = df.z_residual.values
    zs = df.z_selectivity.values
    res_v = df.residual.values
    idx_slc = list(df.index).index("SLC7A11") if "SLC7A11" in df.index else None
    null_rank = np.empty(N_PERM)
    for i in range(N_PERM):
        perm = RNG.permutation(len(zs))
        sp, selp = zs[perm], df.dep_selectivity.values[perm]
        gp = (res_v > 0) & (selp > 0)
        sc = np.where(gp, np.minimum(zr, sp), -np.inf)
        null_rank[i] = (sc > sc[idx_slc]).sum() + 1
    emp_p = float((null_rank <= r_slc).mean())
    out["P3_label_permutation"] = {
        "n_perm": N_PERM, "null_median_rank": float(np.median(null_rank)),
        "observed_rank": r_slc, "empirical_p": emp_p,
        "pass": bool(emp_p < 0.05)}
    log("P3 打亂對應：虛無名次中位 %.0f；觀察 %d → 經驗 p = %.4f → %s"
        % (np.median(null_rank), r_slc, emp_p,
           "過" if out["P3_label_permutation"]["pass"] else "不過"))

    down = df[df.residual < -0.05]
    out["P4_downregulated_not_on_top"] = {
        "n_down": int(len(down)),
        "n_down_in_top_1pct": int(len([g for g in down.index if g in top1])),
        "pass": bool(len([g for g in down.index if g in top1]) == 0)}
    log("P4 腫瘤中下降的基因 %d 個，前 1%% 裡有 %d 個 → %s"
        % (len(down), out["P4_downregulated_not_on_top"]["n_down_in_top_1pct"],
           "過" if out["P4_downregulated_not_on_top"]["pass"] else "不過"))

    out["top30"] = json.loads(df.head(30)[cols_show].round(4).to_json(orient="index"))
    out["axis_genes"] = json.loads(
        df.reindex([g for g in ["SLC7A11", "SLC3A2", "GCLC", "GCLM", "AIFM2",
                                "GPX4", "MAFG", "CCND1", "SOX11"]
                    if g in df.index])[cols_show].round(4).to_json(orient="index"))
    passed = [k for k, v in out.items() if isinstance(v, dict) and v.get("pass") is True]
    failed = [k for k, v in out.items() if isinstance(v, dict) and v.get("pass") is False]
    out["summary"] = {"passed": passed, "failed": failed}
    log("\n通過：%s｜未通過：%s"
        % (", ".join(passed) or "（無）", ", ".join(failed) or "（無）"))
    json.dump(out, open(f"{RES}/212_window_screen.json", "w"),
              indent=2, ensure_ascii=False)
    log("DONE → results/212_*")


if __name__ == "__main__":
    main()
