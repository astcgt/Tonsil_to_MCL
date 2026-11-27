#!/usr/bin/env python
"""步驟 236：既然交叉反應是必要條件，現有設計品要重新評價一次。

為什麼要重讀同一批資料
----------------------
步驟 234 S3 證明 small Maf 家族部分冗餘（`MAFF`+`MAFK` 表現 vs `MAFG` gene effect
r = +0.0776、經驗 p = 0.0190）。這把設計目標整個反過來：

  舊判準：設計品應該**只**黏 `MafG`，黏到 `MafF`／`MafK` 是 off-target。
  新判準：只打掉 `MafG` 會被另外兩個補上 → **必須三個都黏**才可能有效。
          真正的 off-target 變成**非 small Maf 的那些蛋白**（BACH1、FOS、JUN、NFE2L2）。

資料一行都不用重跑：步驟 132（Boltz-2）與步驟 150（Chai-1）已經有
4 個 binder × 7 個 target × n=20 取樣。這一支只換判準重讀。

**用哪個模型判定**：只用 **Chai-1**。這是步驟 132 的既有結論 ——
Boltz-2 的陰性對照在 `MafK` 上拿 0.646（n=20），和真設計無法區分，判定不可用；
ESMFold 複合體模式連陽性對照都不過。Chai-1 是唯一兩個對照都過的。
Boltz-2 的數字仍然列出以示透明，但**不參與判定**。

尺度：各 target 的絕對 iptm 不可互相比較（不同蛋白的基線不同），
所以一律換算成**該 target 自己的打亂序列陰性對照**為單位：
    z = (score_binder − score_scram) / sd_scram
這與步驟 220／221 用雜訊底標準化的做法一致。

事先寫死的判準（在讀任何數字之前寫成）
--------------------------------------
R1 **泛 Maf 覆蓋（新的主要判準）**：設計品在 `MAFG`、`MAFF`、`MAFK` **三個**上
   都必須 z ≥ **2.0**。三個都過 → `pan_maf_coverage`；只過 `MAFG` → `mafg_only`（舊目標，
   但依 234 是無效的）；一個都不過 → `no_binding`。
R2 **真正的 off-target（判準不變，只換成員）**：在 `BACH1`、`FOS`、`JUN`、`NFE2L2`
   上必須 z < 2.0。任何一個超過 → 記為 off-target 命中。
R3 **陽性對照**：天然 `NRF2`（CTRLPOS_NRF2）對 `MAFG` 必須 z ≥ 2.0，
   否則這個尺度本身不可用，整支作廢。
R4 **陰性對照自檢**：打亂序列對自己當基線，z 恆為 0 —— 用來確認換算沒寫錯。

事先預測（記入誠實清單）
------------------------
- R3 我預測**會過**（Chai-1 兩個對照都過是步驟 132 的既有結論）。
- R1 我預測**會是 `pan_maf_coverage`**：舊報告把「設計品也黏 MafK」當成缺點寫，
  所以那個結合是存在的。若成立，**現有設計品在新判準下比舊判準下好** ——
  這會是這批資料第一次被正確評價。
- R2 我預測 `BACH1` 可能命中（它和 small Maf 共用 MARE 結合模式），其餘不命中。

輸出 results/236_*
"""
import json
import sys

import numpy as np
import pandas as pd

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
RES = f"{BASE}/results"

SMALL_MAF = ["MAFG", "MAFF", "MAFK"]
OFF_TARGET = ["BACH1", "FOS", "JUN", "NFE2L2"]
NEG = "CTRLNEG_scram"
POS = "CTRLPOS_NRF2"
DESIGNS = ["D1_best", "D2_second"]
Z_CUT = 2.0


def log(*a):
    print("[236]", *a, flush=True)


def load(path, score_col, sd_col):
    d = pd.read_csv(path)
    d = d.rename(columns={score_col: "score", sd_col: "sd"})
    return d[["binder", "target", "n", "score", "sd"]]


def zscale(d):
    """換算成各 target 自己的陰性對照為單位。"""
    neg = d[d.binder == NEG].set_index("target")
    rows = []
    for _, r in d.iterrows():
        nb = neg.loc[r.target] if r.target in neg.index else None
        if nb is None or not np.isfinite(nb["sd"]) or nb["sd"] <= 0:
            z = np.nan
        else:
            z = (r.score - nb["score"]) / nb["sd"]
        rows.append(dict(binder=r.binder, target=r.target, n=int(r.n),
                         score=float(r.score),
                         neg_score=float(nb["score"]) if nb is not None else np.nan,
                         neg_sd=float(nb["sd"]) if nb is not None else np.nan,
                         z=float(z)))
    return pd.DataFrame(rows)


def main():
    out = {}
    tabs = {}
    for name, path, sc, sd in [
            ("chai", f"{RES}/150_chai_summary_n20.csv", "iptm_mean", "iptm_sd"),
            ("boltz", f"{RES}/132_specificity_summary_n20.csv", "iptm_mean", "iptm_std")]:
        try:
            z = zscale(load(path, sc, sd))
            tabs[name] = z
            z.to_csv(f"{RES}/236_z_{name}.csv", index=False)
        except Exception as e:
            log("讀 %s 失敗：%s" % (path, e))

    if "chai" not in tabs:
        raise SystemExit("缺 Chai-1 資料，無法判定")
    Z = tabs["chai"]

    def zv(binder, target):
        m = Z[(Z.binder == binder) & (Z.target == target)]
        return float(m.z.iloc[0]) if len(m) else np.nan

    # ------------------------------------------------- R4 自檢 ---
    neg_self = [zv(NEG, t) for t in SMALL_MAF + OFF_TARGET]
    r4 = bool(np.nanmax(np.abs(neg_self)) < 1e-9)
    out["R4_scale_sanity"] = {"rule": "打亂序列對自己的 z 必須恆為 0",
                              "values": neg_self, "pass": r4}
    log("=== R4 尺度自檢：打亂序列自己的 z = %s → %s"
        % ([round(x, 6) for x in neg_self], "過" if r4 else "不過"))

    # ------------------------------------------------- R3 陽性對照 ---
    zpos = zv(POS, "MAFG")
    r3 = bool(np.isfinite(zpos) and zpos >= Z_CUT)
    out["R3_positive_control"] = {
        "rule": "天然 NRF2 對 MAFG 的 z ≥ %.1f" % Z_CUT,
        "z": zpos, "pass": r3}
    log("\n=== R3 陽性對照：天然 NRF2 → MAFG 的 z = %.2f（門檻 %.1f）→ %s"
        % (zpos, Z_CUT, "過" if r3 else "不過"))
    if not r3:
        out["verdict"] = "VOID_scale_unusable"
        log("\n*** 陽性對照不過 → 這個尺度不可用，不下任何結論 ***")
        json.dump(out, open(f"{RES}/236_pan_maf.json", "w"), indent=2, ensure_ascii=False)
        return 1

    # ------------------------------------------------- 主表 ---
    log("\n=== Chai-1，以各 target 自己的打亂序列陰性對照為單位（z）===")
    log("%-14s %8s %8s %8s | %8s %8s %8s %8s"
        % ("binder", *SMALL_MAF, *OFF_TARGET))
    for b in [POS] + DESIGNS + [NEG]:
        log("%-14s %8.2f %8.2f %8.2f | %8.2f %8.2f %8.2f %8.2f"
            % (b, *[zv(b, t) for t in SMALL_MAF], *[zv(b, t) for t in OFF_TARGET]))
    log("（左半＝新判準要求全部 ≥ %.1f；右半＝真正的 off-target，要求全部 < %.1f）"
        % (Z_CUT, Z_CUT))

    # ------------------------------------------------- R1 / R2 ---
    res = {}
    for b in DESIGNS:
        maf = {t: zv(b, t) for t in SMALL_MAF}
        off = {t: zv(b, t) for t in OFF_TARGET}
        n_maf = sum(1 for v in maf.values() if np.isfinite(v) and v >= Z_CUT)
        hits = [t for t, v in off.items() if np.isfinite(v) and v >= Z_CUT]
        call = ("pan_maf_coverage" if n_maf == 3 else
                "partial_maf_%d_of_3" % n_maf if n_maf > 0 else "no_binding")
        if n_maf == 1 and np.isfinite(maf["MAFG"]) and maf["MAFG"] >= Z_CUT:
            call = "mafg_only"
        res[b] = {"small_maf_z": maf, "off_target_z": off,
                  "n_small_maf_passing": n_maf,
                  "off_target_hits": hits, "call": call}
        log("\n  %s：small Maf 過 %d/3 → **%s**" % (b, n_maf, call))
        log("    off-target 命中：%s" % (hits or "無"))
    out["R1_pan_maf_coverage"] = {
        "rule": "MAFG/MAFF/MAFK 三個都 z ≥ %.1f" % Z_CUT, "per_design": res}
    out["R2_off_target"] = {
        "rule": "BACH1/FOS/JUN/NFE2L2 都必須 z < %.1f" % Z_CUT,
        "hits": {b: res[b]["off_target_hits"] for b in DESIGNS}}

    # ------------------------------------------------- Boltz 透明度 ---
    if "boltz" in tabs:
        B = tabs["boltz"]
        def bz(b, t):
            m = B[(B.binder == b) & (B.target == t)]
            return float(m.z.iloc[0]) if len(m) else np.nan
        out["boltz_for_transparency"] = {
            "note": "Boltz-2 的陰性對照在 MafK 上與真設計無法區分（步驟 132，n=20 時 0.646），"
                    "判定不可用，此處僅列出不參與判定",
            "z": {b: {t: bz(b, t) for t in SMALL_MAF + OFF_TARGET}
                  for b in DESIGNS + [POS]}}
        log("\n=== Boltz-2（僅供透明，不參與判定）===")
        for b in DESIGNS + [POS]:
            log("%-14s %s" % (b, {t: round(bz(b, t), 2) for t in SMALL_MAF}))

    # --------------------------------- 設計品之間的取捨（新判準下的真正問題）---
    log("\n=== 新判準下的取捨：覆蓋率 vs off-target ===")
    log("%-14s %-22s %-18s %s" % ("binder", "small Maf 覆蓋", "最大 off-target", "對照：天然 NRF2"))
    zpos_off = {t: zv(POS, t) for t in OFF_TARGET}
    trade = {}
    for b in DESIGNS:
        hits = res[b]["off_target_hits"]
        worst = max(OFF_TARGET, key=lambda t: (zv(b, t) if np.isfinite(zv(b, t)) else -99))
        wz = zv(b, worst)
        exceeds_pos = bool(np.isfinite(wz) and np.isfinite(zpos_off.get(worst, np.nan))
                           and wz > zpos_off[worst])
        trade[b] = {"n_small_maf": res[b]["n_small_maf_passing"],
                    "worst_off_target": worst, "worst_off_target_z": wz,
                    "natural_NRF2_z_same_target": zpos_off.get(worst),
                    "exceeds_natural_partner": exceeds_pos}
        log("%-14s %-22s %-18s %s"
            % (b, "%d/3" % res[b]["n_small_maf_passing"],
               "%s z=%.2f" % (worst, wz),
               "z=%.2f%s" % (zpos_off.get(worst, np.nan),
                             "  ← 設計品超過天然夥伴" if exceeds_pos else "")))
    out["trade_off"] = trade
    log("\n  判讀：off-target 若**超過天然 NRF2 在同一個 target 上的分數**，")
    log("  就不能用『共摺疊模型對 bZIP 普遍給高分』解釋掉，那是設計品特有的結合。")

    clean = [b for b in DESIGNS
             if trade[b]["n_small_maf"] == 3 and not res[b]["off_target_hits"]]
    best = res.get("D1_best", {}).get("call")
    out["verdict"] = {
        "any_design_meets_corrected_spec": clean,
        "D1_best": best,
        "D2_second": res.get("D2_second", {}).get("call"),
        "reading": ("舊判準把『也黏 MafF/MafK』當缺點；在 234 的冗餘結果之後，"
                    "pan_maf_coverage 才是有效的形狀。")}
    log("\n*** 判定（D1_best）：%s ***" % best)
    if best == "pan_maf_coverage":
        log("    三個 small Maf 都黏 → 在新判準下，現有設計品的形狀是**對的**。")
        log("    舊報告把這件事寫成缺點，是因為判準建立在未經檢定的『唯一必要夥伴』前提上。")
    elif best and best.startswith("partial"):
        log("    只覆蓋部分 small Maf → 依 234，未覆蓋的那些會補上來，效果會被削弱。")
    elif best == "mafg_only":
        log("    只黏 MafG → 這正是 234 說會失敗的形狀。")

    json.dump(out, open(f"{RES}/236_pan_maf.json", "w"), indent=2, ensure_ascii=False)
    log("\n→ results/236_pan_maf.json、236_z_chai.csv、236_z_boltz.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
