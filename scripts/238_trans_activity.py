#!/usr/bin/env python
"""步驟 238：既然因果在反式那一側，就該量反式因子的**活性**，不只表現量。

為什麼
------
步驟 235 證明腫瘤新開的可及性**不寫在局部序列裡**（語言模型 AUC 0.5071，
陽性對照 0.8906），四條證據收斂到「序列是必要的，但不是區別所在」。
因果鏈因此指向反式：**哪些因子的活性在腫瘤裡變了。**

步驟 19 已經算過逐捐贈者、腫瘤 vs 非腫瘤的 TF 活性（regulon 分數），
但從來沒有針對這條軸線的家族結構去讀它。這一支補上。

**n = 2 位捐贈者這件事必須先處理，而不是最後才提**
--------------------------------------------------
只有 M102 與 M413 兩位。所以「兩位都同方向」這件事**本身不是證據** ——
純機率下就有 50% 會同方向。因此本支的主要統計量**不是方向一致性**，而是：
    每位捐贈者的 |delta| 落在**全部 TF 的 |delta| 分布**裡的百分位。
方向一致性只當附加條件，不當主要證據。

事先寫死的檢定與判準
--------------------
A1 **家族層級**：small Maf（MAFG/MAFF/MAFK）、AP-1（JUN/JUNB/FOS/FOSL1/FOSL2）、
   BACH（BACH1/BACH2）、NFE2L2 —— 各家族成員的 delta 百分位。
   規則：一個家族要算「活性有變」，必須至少一個成員在**兩位捐贈者**都落在
   全 TF |delta| 分布的 **第 80 百分位以上**，且兩位同方向。
A2 **軸線的關鍵問題**：`MAFG` 自己的活性在腫瘤裡有變嗎？
   若 `MAFG` 的 delta 百分位很低 → 它是**恆定的必要條件，不是訊號源** →
   因果應歸給 AP-1／BACH 那一側，而 `MafG` 只是它們共用的接口。
A3 **全域對照**：列出所有 TF 裡 delta 最大且兩位同方向的前 15 名，
   檢查這條軸線的成員是否真的名列前茅，或者只是被我們挑出來看而已。
   （這是防「只看自己想看的基因」的對照。）
A4 **陰性對照**：一組與這條軸線無關的 TF（CTCF、YY1、SP1、NRF1、E2F1、MYC）
   不該整組都排在前面。若它們也都在前面 → 這個分數在量的是「腫瘤 vs 正常的
   全域差異」而不是特定調控，判定降級為不可解讀。

事先預測（記入誠實清單）
------------------------
- A2 我預測 `MAFG` 的活性變化**很小**（它是普遍表現的接口蛋白），
  也就是 `MafG` 是恆定要求而不是訊號源。
- A1 我預測 **AP-1 家族會是變化最大的那一組**（報告的核心發現是 TGACTCA 位點富集）。
- A4 我預測陰性對照不會整組排在前面。

輸出 results/238_*
"""
import json
import sys

import numpy as np
import pandas as pd

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
RES = f"{BASE}/results"

FAMILIES = {
    "small_Maf": ["MAFG", "MAFF", "MAFK"],
    "AP-1": ["JUN", "JUNB", "JUND", "FOS", "FOSB", "FOSL1", "FOSL2"],
    "BACH": ["BACH1", "BACH2"],
    "NFE2_family": ["NFE2L2", "NFE2", "NFE2L1", "NFE2L3"],
}
NEG_CTRL = ["CTCF", "YY1", "SP1", "NRF1", "E2F1", "MYC"]
PCT_CUT = 80.0


def log(*a):
    print("[238]", *a, flush=True)


def main():
    d = pd.read_csv(f"{RES}/19_tf_ranking_paired.csv", index_col=0)
    log("TF 表 %s｜欄位 %s" % (d.shape, list(d.columns)[:8]))
    need = ["d_M102", "d_M413", "same_dir"]
    for c in need:
        if c not in d.columns:
            raise SystemExit("缺欄位 %s" % c)
    d = d.dropna(subset=["d_M102", "d_M413"])
    n_tf = len(d)
    log("可用 TF %d 個" % n_tf)

    # 全 TF 的 |delta| 分布 → 百分位（兩位捐贈者各自算）
    for dn in ("M102", "M413"):
        a = d["d_" + dn].abs()
        d["pct_" + dn] = a.rank(pct=True) * 100
    d["min_pct"] = d[["pct_M102", "pct_M413"]].min(1)
    d["same_dir_bool"] = np.sign(d.d_M102) == np.sign(d.d_M413)

    out = {"n_tf": int(n_tf), "pct_cut": PCT_CUT,
           "n_donors": 2,
           "caveat": "只有 2 位捐贈者。『兩位同方向』純機率下有 50% 會發生，"
                     "所以主要統計量是 |delta| 在全 TF 分布裡的百分位，"
                     "方向一致性只當附加條件。"}

    # ------------------------------------------------------------ A1 ---
    log("\n=== A1 家族層級：各成員的 delta 與百分位（%d 個 TF 為分母）===" % n_tf)
    log("%-10s %-9s %9s %9s %7s %7s %8s" %
        ("family", "TF", "d_M102", "d_M413", "pct102", "pct413", "同方向"))
    fam_out = {}
    for fam, members in FAMILIES.items():
        rows = []
        for g in members:
            if g not in d.index:
                continue
            r = d.loc[g]
            rows.append({"tf": g, "d_M102": float(r.d_M102), "d_M413": float(r.d_M413),
                         "pct_M102": float(r.pct_M102), "pct_M413": float(r.pct_M413),
                         "min_pct": float(r.min_pct),
                         "same_dir": bool(r.same_dir_bool)})
            log("%-10s %-9s %+9.3f %+9.3f %7.1f %7.1f %8s"
                % (fam, g, r.d_M102, r.d_M413, r.pct_M102, r.pct_M413,
                   "是" if r.same_dir_bool else "否"))
        qual = [x for x in rows if x["min_pct"] >= PCT_CUT and x["same_dir"]]
        fam_out[fam] = {"members": rows, "qualifying": [x["tf"] for x in qual],
                        "changed": bool(qual),
                        "best_min_pct": max([x["min_pct"] for x in rows], default=np.nan)}
    out["A1_by_family"] = fam_out
    log("\n  家族判定（至少一個成員兩位都 ≥ 第 %.0f 百分位且同方向）：" % PCT_CUT)
    for fam, v in fam_out.items():
        log("    %-12s %s  %s" % (fam, "有變" if v["changed"] else "沒有變",
                                  v["qualifying"] or ""))

    # ------------------------------------------------------------ A2 ---
    log("\n=== A2 軸線的關鍵問題：`MAFG` 自己的活性有變嗎 ===")
    if "MAFG" in d.index:
        m = d.loc["MAFG"]
        a2_const = bool(m.min_pct < 50)
        out["A2_mafg_is_constitutive"] = {
            "d_M102": float(m.d_M102), "d_M413": float(m.d_M413),
            "pct_M102": float(m.pct_M102), "pct_M413": float(m.pct_M413),
            "min_pct": float(m.min_pct),
            "interpretation": ("恆定的必要接口，不是訊號源" if a2_const
                               else "活性本身也變了，可能參與驅動"),
            "is_constitutive": a2_const}
        log("  MAFG：delta %+.4f / %+.4f｜百分位 %.1f / %.1f → **%s**"
            % (m.d_M102, m.d_M413, m.pct_M102, m.pct_M413,
               out["A2_mafg_is_constitutive"]["interpretation"]))
        if a2_const:
            log("  → 因果應歸給 AP-1／BACH 那一側；`MafG` 是它們共用的接口，")
            log("     這也正是它當標的的理由，但它不是變化的來源。")

    # ------------------------------------------------------------ A3 ---
    log("\n=== A3 全域對照：所有 TF 裡變化最大且兩位同方向的前 15 名 ===")
    top = (d[d.same_dir_bool]
           .assign(absmin=lambda x: x[["d_M102", "d_M413"]].abs().min(1))
           .sort_values("absmin", ascending=False).head(15))
    axis_all = set(sum(FAMILIES.values(), []))
    log("%-4s %-10s %9s %9s %7s %s" % ("#", "TF", "d_M102", "d_M413", "min_pct", "屬軸線"))
    for i, (g, r) in enumerate(top.iterrows(), 1):
        log("%-4d %-10s %+9.3f %+9.3f %7.1f %s"
            % (i, g, r.d_M102, r.d_M413, r.min_pct, "★" if g in axis_all else ""))
    n_axis_top = sum(1 for g in top.index if g in axis_all)
    out["A3_global_top15"] = {
        "top": [{"tf": g, "d_M102": float(r.d_M102), "d_M413": float(r.d_M413),
                 "min_pct": float(r.min_pct), "on_axis": bool(g in axis_all)}
                for g, r in top.iterrows()],
        "n_axis_members_in_top15": int(n_axis_top)}
    log("  軸線成員進前 15 名的個數：%d" % n_axis_top)

    # ------------------------------------------------------------ A4 ---
    log("\n=== A4 陰性對照：與這條軸線無關的 TF ===")
    nc = []
    for g in NEG_CTRL:
        if g in d.index:
            r = d.loc[g]
            nc.append({"tf": g, "min_pct": float(r.min_pct),
                       "same_dir": bool(r.same_dir_bool)})
            log("  %-8s 百分位 %.1f / %.1f｜同方向 %s"
                % (g, r.pct_M102, r.pct_M413, "是" if r.same_dir_bool else "否"))
    n_hi = sum(1 for x in nc if x["min_pct"] >= PCT_CUT and x["same_dir"])
    a4 = bool(len(nc) == 0 or n_hi < max(1, len(nc) // 2))
    out["A4_negative_control"] = {
        "rule": "無關 TF 不得有一半以上同時達標",
        "values": nc, "n_qualifying": n_hi, "pass": a4}
    log("  %d / %d 個無關 TF 也達標 → %s"
        % (n_hi, len(nc), "對照通過" if a4 else "**分數在量全域差異，判定降級**"))

    changed = [f for f, v in fam_out.items() if v["changed"]]
    out["verdict"] = {
        "families_with_changed_activity": changed,
        "negative_control_ok": a4,
        "call": ("trans_activity_localised_to_" + "+".join(changed)
                 if changed and a4 else
                 "no_family_passes" if a4 else "VOID_score_measures_global_shift")}
    log("\n*** 判定：%s ***" % out["verdict"]["call"])

    d_out = d[["d_M102", "d_M413", "pct_M102", "pct_M413", "min_pct", "same_dir_bool"]]
    d_out.to_csv(f"{RES}/238_tf_activity_pct.csv")
    json.dump(out, open(f"{RES}/238_trans_activity.json", "w"),
              indent=2, ensure_ascii=False)
    log("\n→ results/238_trans_activity.json、238_tf_activity_pct.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
