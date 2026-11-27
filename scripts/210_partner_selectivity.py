#!/usr/bin/env python
"""步驟 210（D5b，淺層）：現有設計品擋住的是 BACH1 還是 NRF2？零 GPU。

為什麼要先問這一題
------------------
步驟 141 在 MafG 上找到 9 個「BACH1 埋 >10 Å²、NRF2 埋 <5 Å²」的專屬位置
（NRF2 在共結晶 7X5E 裡對它們的埋藏是 0.0 Å²）。原本的計畫是拿這 9 個位置當
熱點重新設計結合劑，去取回 **MafG 專一性**。

但先做了一個便宜的檢查，結論是那條路走不通：
這 9 個位置在 MafG 與 MafK 之間 **9/9 完全相同**，與 MafF **8/9 相同**
（第 9 個 V121 只是 C 端比對缺口）。表位是整個 small Maf 家族共有的，
所以瞄準它不可能分開三個旁系。

但它換來另一個東西：**夥伴選擇性** —— 擋住 BACH1、保留 NRF2。
在這個框架下，「family binder」不再是缺陷而是預期行為。
所以在花 30+ GPU 小時重新設計之前，先問現成的設計品 D1 有沒有這個性質。

#### 必須處理的混淆
BACH1 的逐殘基埋藏輪廓與 MafG 同二聚體的相關是 **0.970**，與 NRF2 只有 0.772。
也就是說 BACH1 和同二聚體用的幾乎是同一面（標準的白胺酸拉鍊面），而 NRF2 用的面
不太一樣。於是**任何**黏在拉鍊面上的東西都會「比較擋得到 BACH1」。
所以單看 D1 偏好 BACH1 毫無意義，必須和一個「泛拉鍊面結合者」比 ——
MafG 同二聚體就是現成的那個比較基準。

事前指定的讀數與預測
--------------------
遮蔽比例 occl(P) = Σᵢ BSA_P(i)·1[BSA_D1(i) > 5 Å²] / Σᵢ BSA_P(i)
  （夥伴 P 在 MafG 上的埋藏面積，有多少比例落在設計品也埋住的位置上）
選擇性 S(X) = occl_X(BACH1) − occl_X(NRF2)，X ∈ {D1 設計品, MafG 同二聚體}

  P1  S(D1) > 0                      —— 設計品偏好遮蔽 BACH1
  P2  S(D1) > S(同二聚體)            —— 而且比泛拉鍊面結合者更偏好（關鍵的一關）
  P3  S(D1) 大於隨機足跡的 95 分位   —— 不是足跡大小造成的
  P4（陰性對照）用 NRF2 自己當「設計品」時，S 必須是負的 —— 檢查讀數方向沒反

⚠️ 這是幾何重疊，不是競爭結合自由能。它回答「設計品坐在誰的位子上」，
   不回答「它能不能把誰推開」。後者是步驟 211 的三鏈共摺疊。

輸出 results/210_*
"""
import os
import json

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "8"

import numpy as np
import pandas as pd

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
RES = f"{BASE}/results"
BSA_CUT = 5.0          # Å²，算「這個位置被佔住了」的門檻
N_PERM = 20000
RNG = np.random.default_rng(0)

BACH1_SPEC = [61, 64, 65, 68, 71, 115, 117, 118, 121]
NRF2_SPEC = [50]


def log(*a):
    print("[210]", *a, flush=True)


def occl(partner_bsa, blocker_bsa, cut=BSA_CUT):
    tot = partner_bsa.sum()
    if tot <= 0:
        return float("nan")
    return float(partner_bsa[blocker_bsa > cut].sum() / tot)


def main():
    d = pd.read_csv(f"{RES}/141_per_residue_bsa.csv", index_col=0)
    N = d["MAFG_x_NRF2"].values
    B = d["MAFG_x_BACH1"].values
    H = d["MAFG_x_MAFG"].values
    D = d["MAFG_x_D1design"].values
    pos = d["uniprot_pos"].values

    out = {"bsa_cut_A2": BSA_CUT, "n_residues": int(len(d))}

    log("足跡大小（BSA > %.0f Å² 的位置數）：NRF2 %d、BACH1 %d、同二聚體 %d、D1 %d"
        % (BSA_CUT, (N > BSA_CUT).sum(), (B > BSA_CUT).sum(),
           (H > BSA_CUT).sum(), (D > BSA_CUT).sum()))
    out["footprint_sizes"] = {"NRF2": int((N > BSA_CUT).sum()),
                              "BACH1": int((B > BSA_CUT).sum()),
                              "MAFG_homodimer": int((H > BSA_CUT).sum()),
                              "D1_design": int((D > BSA_CUT).sum())}

    blockers = {"D1_design": D, "MAFG_homodimer": H, "NRF2_as_blocker": N}
    res = {}
    for nm, X in blockers.items():
        oB, oN = occl(B, X), occl(N, X)
        res[nm] = {"occl_BACH1": oB, "occl_NRF2": oN, "selectivity": oB - oN}
        log("%-16s 遮蔽 BACH1 %.3f｜遮蔽 NRF2 %.3f｜選擇性 S = %+.3f"
            % (nm, oB, oN, oB - oN))
    out["occlusion"] = res

    S_d1 = res["D1_design"]["selectivity"]
    S_hd = res["MAFG_homodimer"]["selectivity"]
    S_nr = res["NRF2_as_blocker"]["selectivity"]

    # P3：隨機足跡的虛無分佈（保留足跡大小，打亂位置）
    null = np.empty(N_PERM)
    for i in range(N_PERM):
        Xp = RNG.permutation(D)
        null[i] = occl(B, Xp) - occl(N, Xp)
    q95 = float(np.quantile(null, 0.95))
    emp_p = float((null >= S_d1).mean())
    out["P3_permutation_null"] = {
        "n_perm": N_PERM, "null_mean": float(null.mean()),
        "null_sd": float(null.std()), "null_q95": q95,
        "observed": S_d1, "empirical_p": emp_p, "pass": bool(emp_p < 0.05)}

    out["P1_design_prefers_BACH1"] = {"S": S_d1, "pass": bool(S_d1 > 0)}
    out["P2_better_than_generic_zipper_binder"] = {
        "S_design": S_d1, "S_homodimer": S_hd,
        "difference": S_d1 - S_hd, "pass": bool(S_d1 > S_hd)}
    out["P4_negative_control_direction"] = {
        "S_NRF2_as_blocker": S_nr, "pass": bool(S_nr < 0)}

    log("")
    log("P1 S(D1) = %+.3f > 0 → %s" % (S_d1, "過" if out["P1_design_prefers_BACH1"]["pass"] else "不過"))
    log("P2 S(D1) %+.3f vs S(同二聚體) %+.3f（差 %+.3f）→ %s"
        % (S_d1, S_hd, S_d1 - S_hd,
           "過" if out["P2_better_than_generic_zipper_binder"]["pass"] else "不過"))
    log("P3 隨機足跡 %.3f ± %.3f（95 分位 %.3f）；觀察 %+.3f → 經驗 p = %.4f → %s"
        % (null.mean(), null.std(), q95, S_d1, emp_p,
           "過" if out["P3_permutation_null"]["pass"] else "不過"))
    log("P4 用 NRF2 當遮蔽者 S = %+.3f（應為負）→ %s"
        % (S_nr, "過" if out["P4_negative_control_direction"]["pass"] else "不過"))

    # 專屬位置的覆蓋
    ix_b = np.isin(pos, BACH1_SPEC)
    ix_n = np.isin(pos, NRF2_SPEC)
    cov = {"BACH1_specific_positions": BACH1_SPEC,
           "n_BACH1_specific_covered_by_D1": int((D[ix_b] > BSA_CUT).sum()),
           "n_BACH1_specific": int(ix_b.sum()),
           "n_NRF2_specific_covered_by_D1": int((D[ix_n] > BSA_CUT).sum()),
           "n_NRF2_specific": int(ix_n.sum()),
           "D1_bsa_on_BACH1_specific_A2": float(D[ix_b].sum()),
           "D1_bsa_on_NRF2_specific_A2": float(D[ix_n].sum())}
    out["epitope_coverage"] = cov
    log("\nD1 覆蓋 BACH1 專屬位置 %d/%d（共 %.0f Å²）、NRF2 專屬 %d/%d（%.0f Å²）"
        % (cov["n_BACH1_specific_covered_by_D1"], cov["n_BACH1_specific"],
           cov["D1_bsa_on_BACH1_specific_A2"],
           cov["n_NRF2_specific_covered_by_D1"], cov["n_NRF2_specific"],
           cov["D1_bsa_on_NRF2_specific_A2"]))

    passed = [k for k, v in out.items() if isinstance(v, dict) and v.get("pass") is True]
    failed = [k for k, v in out.items() if isinstance(v, dict) and v.get("pass") is False]
    out["summary"] = {"passed": passed, "failed": failed,
                      "verdict": ("worth_deepening" if (S_d1 > S_hd and emp_p < 0.05)
                                  else "not_worth_deepening")}
    log("\n通過：%s" % (", ".join(passed) or "（無）"))
    log("未通過：%s" % (", ".join(failed) or "（無）"))
    log("判定：%s" % out["summary"]["verdict"])
    json.dump(out, open(f"{RES}/210_partner_selectivity.json", "w"),
              indent=2, ensure_ascii=False)
    log("DONE → results/210_*")


if __name__ == "__main__":
    main()
