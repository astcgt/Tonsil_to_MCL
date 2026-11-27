#!/usr/bin/env python
"""步驟 204：`SLC7A11` 依賴性為什麼只出現在 B 系？逐一檢定可能的解釋。

要解釋的事實（步驟 135，不重算）
--------------------------------
1,178 株裡只有 18 株依賴 `SLC7A11`（gene effect < −0.5），其中 12 株是淋巴系，
而且 12 株全部是 B 細胞腫瘤（OR 26.6、Fisher p = 3.7e−10、對 1,320 個依賴率
配對隨機基因的經驗 p = 0.021）。GCLC / GCLM / AIFM2 同方向，SLC3A2 沒有選擇性，
MAFG 完全沒有富集。五個基因的 CRISPR 輪廓構成一個真實模組（r = 0.278，經驗 p < 0.0005）。

已經被否證的解釋（步驟 203，不再重試）
--------------------------------------
「B 細胞沒有轉硫途徑（CBS/CBSL/CTH）所以只能靠 xCT 進胱胺酸」四項預測全不過：
全庫 r = −0.030（方向還是反的）、非淋巴系 r = −0.054、對 1,999 個隨機基因的
經驗 p = 0.789、陰性對照 RPL13A 的相關還更強。唯一成立的組間差異完全可以用
譜系組成解釋。

這一步的判準（先寫下來再跑）
----------------------------
一個候選解釋要算「機制」而不是「『是不是 B 細胞』的代名詞」，必須通過
兩類彼此獨立的譜系控制，兩類都報，缺一不可：

  (A) 分層：在**非淋巴系**（n ≈ 1,085）或**B 系內部**（n ≈ 72，其中 12 株依賴）
      這兩個分層裡，特徵仍然預測 `SLC7A11` 的 gene effect。
      非淋巴系只有 6 株依賴者，所以這個分層對二元對比本來就沒有力；
      每個分層都會印出 80% 檢定力下可偵測的最小 r（≈ 2.8/sqrt(n−3)），
      判決因此分三種：pass / fail / underpowered，不把沒力氣當成否證。
  (B) 衰減：把 `effect ~ is_B` 的 B 係數，和 `effect ~ is_B + feature` 的
      B 係數相比，看加進特徵之後 B 的效應縮掉多少。虛無分佈取自
      **B 系特異度相當的隨機基因**（|B vs 非 B 標準化均值差| 配對 ±25%），
      所以「只要夠 B 系特異就能衰減」會被扣掉。

結果變數（主要 + 次要 + 陰性對照），全部先指定
----------------------------------------------
  主要   `SLC7A11` gene effect。
  次要   `GCLC` gene effect（42 株依賴，其中 9 株非淋巴系，比主要變數有力），
         以及模組分數＝SLC7A11/GCLC/GCLM/AIFM2 gene effect 的 z 分數平均。
         真的是穀胱甘肽臂的機制，應該連這兩個一起預測到。
  陰性對照之一  `RPL13A`（泛必需基因）：任何特徵都不該預測它。
  陰性對照之二  B 身分依賴分數（CD79A/PAX5/EBF1 的 gene effect z 平均）。
         這一個是用來抓代名詞的：一個特徵如果把 `SLC7A11` 和「B 身分基因本身
         的依賴性」預測得一樣好，那它測到的是 B 細胞身分，不是胱胺酸代謝。

虛無分佈（每一個效應量都有，兩種各自意義不同）
----------------------------------------------
  Null A  隨機基因（或同樣大小的隨機基因組，z 標準化後平均）的表現 vs 同一個
          結果變數 → 回答「這個基因的表現特別嗎」。面板特徵配對面板大小。
  Null B  同一個特徵 vs **依賴率配對的隨機基因**的 gene effect
          → 回答「`SLC7A11` 在依賴率相近的基因裡特別嗎」（步驟 135 的做法）。
          只對 `SLC7A11` 與 `GCLC` 這兩個有明確依賴率的結果變數建。

檢定與各自的預設預測
--------------------
T1 收斂範圍：依賴性到底是「淋巴系」還是更窄。
    預測：B 系 > T/NK。這一題的關鍵是**漿細胞骨髓瘤**——它也是 B 系但已經
    離開 B 細胞程式；如果 19 株骨髓瘤裡 0 株依賴，那範圍比「B 細胞」更窄，
    而且這個對比完全在淋巴系內部，不可能被譜系混淆。
    陰性對照：`RPL13A` 不該有任何組間結構。虛無：依賴率配對隨機基因的
    「依賴株中 B（非漿細胞）比例」。
T2 自身表現：依賴的株是不是本來就表現比較多 `SLC7A11`。
    預測：表現越高越依賴（Spearman r < 0），且要在分層裡成立。
    注意有個下限效應：完全不表現的基因不可能成為依賴，所以全庫的相關
    一定有一部分來自這個下限，判決要看分層。
    陰性對照特徵：ACTB/RPL13A/GAPDH 的表現（持家基因面板）。
T3 NRF2 軸：文獻裡 xCT 依賴的標準解釋是 KEAP1/NFE2L2 改變把 xCT 開到底。
    預測（刻意不對稱）：NRF2 活性應該解釋**非淋巴系**那幾株，但 B 系
    不該是 NRF2 高的——步驟 135 已經顯示 `NFE2L2` 依賴性在淋巴系是反富集
    （OR 0.19）。如果這樣，NRF2 是另一條路，不是 B 系的答案。
    NRF2 標的分數刻意不含 `SLC7A11` 以免循環。
T4 備援臂：B 系是不是少了 GPX4 以外的那條防線（FSP1/AIFM2–CoQ）。
    預測：AIFM2 表現高應該有保護作用（r > 0），B 系的 AIFM2 應該偏低。
    反向證據是 `AIFM2` 自己的依賴性也在淋巴系富集（OR 74.8）——兩條臂
    都非冗餘地被需要，這和「缺備援」的預測其實是衝突的，要照實寫。
T5 氧化還原系統的專一性地圖：把依賴性富集攤到五條臂上
    （胱胺酸–GSH / 硫氧還蛋白 / GPX–硒 / CoQ–FSP1 / NADPH / PUFA 脂質），
    每個基因都用依賴率配對隨機基因建經驗 p。
    預測：如果機制是「B 細胞只能靠 GSH」，富集應該**只**出現在胱胺酸–GSH 臂，
    而硫氧還蛋白臂不該富集。這一題是用來定位機制，不是用來證實。
T6 技術性解釋（如果成立必須照實報告，這是合格的結論）：
    (a) 推估生長速率（`CRISPRInferredModelGrowthRate`，三個 library 欄位
        大多互斥，所以各欄分別 z 標準化，不平均）；
    (b) 篩選 library 本身（Avana / Humagne / KY）；
    (c) 懸浮 vs 貼附；
    (d) 培養液 FormulationID；關鍵檢定是**同一種培養液內**（RPMI + 10% FBS，
        393 株，其中 56 株淋巴系）B 系與非淋巴系的依賴率對比，
        外加 T 細胞與髓系株當「培養條件相同的對照組」；
    (e) 每株的全基因組平均 gene effect（株層級的整體篩選敏感度）——如果淋巴系
        對什麼都比較敏感，富集有一部分是株層級假象。
    Model.csv 讀不到硫醇類添加物，所以「培養條件排除」只涵蓋生長型態與配方。
T7 胱胺酸以外的半胱胺酸來源、麩胺酸輸出、ISR/ATF4 驅動。
    預測：替代來源（SLC1A4/SLC1A5/MPST/TST/SQOR/GGT1）越弱越依賴。
    面板裡逐基因檢定時用 max 統計量的虛無（k 個隨機基因的最大 |r|），
    否則挑最好的那個基因等於偷看。
T8 無偏掃描（**探索性，結論只能當假說，不算通過**）：
    (a) B 系內部 72 株：全轉錄體表現 vs `SLC7A11` gene effect；
    (b) 非淋巴系：同上；
    (c) 共依賴：全基因組 gene effect vs `SLC7A11` gene effect，分層算；
    (d) B 系內部 12 株依賴者 vs 其餘 B 系株的差異表現。
    (a) 和 (d) 天生就是譜系控制過的，所以是這一步最有價值的發現管道。

第一次跑完之後才加進來的東西（為了可追溯，照實記下）
----------------------------------------------------
這些不是原始的預設預測，解讀時要當成「事後加的防呆」而不是「事先的預測」：
  * 分層 `B_nonPlasmaCell`（n = 53，12 株依賴）。原因：B 系分層裡那 19 株
    漿細胞骨髓瘤**全部不依賴**，所以任何「非漿細胞 B 細胞偏高或偏低」的特徵
    都會在 B 系分層裡過關——它只是「是不是骨髓瘤」的代名詞。
    必須把骨髓瘤去掉，看特徵能不能分開 12 株依賴者與 41 株不依賴的 B 細胞株。
    （附帶說明：19 株骨髓瘤 0 株依賴這件事，是寫 T1 之前先探過資料才知道的，
    不是 T1 跑出來才發現。）
  * 分層 `non_haematopoietic`（同時排掉淋巴與髓系）。原因：髓系株的 AIFM2
    表現比 B 系還低，卻只有 1/42 依賴，所以 T4 需要一個排掉髓系的對照。
  * 面板逐基因分解 + max 統計量虛無，否則「面板裡挑最好的基因」等於偷看。
  * underpowered 的判準放寬：效應量小於 80% 檢定力下的最小可偵測 r 就一律
    標 underpowered（原本只在方向也相反時才標），不把沒力氣當成否證。
  * 衰減檢定補上「方向是否與預設預測一致」的標記。衰減量本身不看方向，
    關係跑反了也會讓 is_B 的係數縮小，所以方向相反的一律不算支持。
  * 三個頭條對比（淋巴系內部、B 系內部、同培養液）各自補上經驗虛無。

已知的限制（讀結論前要先知道）
------------------------------
  * `RPL13A` 在 CRISPRGeneEffect 裡只有 349/1,178 株有值（不是每個 library
    都篩它），所以它當結果變數的那個陰性對照檢定力偏低。
    組間結構那一關（Kruskal）仍然可用，因為那是用現有的株算的。
  * T5 把 `PRDX6` 放在硫氧還蛋白臂其實是放錯了：PRDX6 是 1-Cys
    peroxiredoxin，用的是穀胱甘肽而不是硫氧還蛋白。它富集（經驗 p 0.015）
    不是硫氧還蛋白臂的反例，反而和「只有穀胱甘肽這條路被需要」一致。
    這個歸類是跑之前就寫死的，不事後搬動，照實記在這裡。
  * T5 的 `SLC7A11`/`GCLC`/`GCLM`/`AIFM2` 是步驟 135 的輸入，富集不算獨立確認；
    `GSS`、`PRDX6`、`ACSL4`、`SCLY`、`SEPHS2`、`PDSS1`、`FDFT1` 才是新的。
  * 培養條件只能查到生長型態與 FormulationID；Model.csv 沒有硫醇類添加物
    （2-ME、麩胱甘肽、N-乙醯半胱胺酸）與血清批次，那部分無法排除。

輸出 results/204_*
"""
import os
import json
import warnings

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "8"
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from scipy import stats

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
RES = f"{BASE}/results"
DM = f"{BASE}/data/depmap"

DEP_CUT = -0.5
N_NULL = 2000
RNG = np.random.default_rng(204)

# ---- 結果變數 ----------------------------------------------------------
TARGET = "SLC7A11"
MODULE = ["SLC7A11", "GCLC", "GCLM", "AIFM2"]
BIDENT = ["CD79A", "PAX5", "EBF1"]
PAN_ESSENTIAL = "RPL13A"

# ---- 候選特徵：表現面板 ------------------------------------------------
# 方向 pred 是「與 gene effect 的 Spearman r 的預期符號」；
# gene effect 越負＝越依賴，所以「越高越依賴」＝ r < 0＝"neg"。
EXPR_PANELS = {
    "exp_SLC7A11":          (["SLC7A11"], "neg",
                             "T2 自身表現：越高越依賴"),
    "exp_SLC3A2":           (["SLC3A2"], "neg",
                             "T2 伴隨次單元表現"),
    "nrf2_target_score":    (["NQO1", "TXNRD1", "SRXN1", "GCLM", "ABCC1",
                              "G6PD", "GSR"], "neg",
                             "T3 NRF2 標的活性（不含 SLC7A11）"),
    "exp_AIFM2":            (["AIFM2"], "pos",
                             "T4 FSP1 備援臂：越高越不依賴"),
    "exp_GPX4":             (["GPX4"], "pos",
                             "T4 GPX4 表現"),
    "txn_arm_score":        (["TXN", "TXNRD1", "PRDX1", "PRDX2", "PRDX6"],
                             "pos", "T4 硫氧還蛋白臂（與 NRF2 面板共用 TXNRD1）"),
    "gsh_synth_exp":        (["GCLC", "GCLM", "GSS"], "neg",
                             "T4 GSH 合成能力"),
    "isr_atf4_score":       (["DDIT3", "TRIB3", "CHAC1", "ASNS", "SESN2",
                              "PHGDH", "PSAT1"], "neg",
                             "T7 ISR/ATF4 驅動（不含 SLC7A11）"),
    "alt_cys_score":        (["SLC1A4", "SLC1A5", "MPST", "TST", "SQOR",
                              "GGT1"], "pos",
                             "T7 胱胺酸以外的半胱胺酸來源：越強越不依賴"),
    "glu_export_score":     (["GLS", "GLUD1", "SLC1A5"], "neg",
                             "T7 麩胺酸供應（xCT 的交換對象）"),
    "pufa_lipid_score":     (["ACSL4", "LPCAT3", "MBOAT7", "ELOVL5",
                              "FADS2"], "neg",
                             "T7 PUFA 磷脂預備（鐵死亡易感性）"),
    "mufa_score":           (["SCD", "SCD5"], "pos",
                             "T7 單元不飽和脂肪酸（保護）"),
    "bcell_identity_score": (["CD79A", "CD79B", "PAX5", "EBF1", "MS4A1",
                              "CD19"], "neg",
                             "基準：純粹的 B 身分，代名詞的對照上限"),
    "housekeeping_ctrl":    (["ACTB", "RPL13A", "GAPDH"], "none",
                             "陰性對照特徵：不該有任何關係"),
}

# ---- T5 氧化還原各臂（檢定的是 gene effect 的 B 系富集）---------------
REDOX_ARMS = {
    "cystine_GSH": ["SLC7A11", "SLC3A2", "GCLC", "GCLM", "GSS", "GSR",
                    "CHAC1", "GGT1", "OPLAH"],
    "thioredoxin": ["TXN", "TXNRD1", "PRDX1", "PRDX2", "PRDX6", "TXN2",
                    "TXNRD2", "SRXN1"],
    "GPX_selenium": ["GPX4", "GPX1", "SEPHS2", "SEPSECS", "PSTK", "EEFSEC",
                     "SCLY"],
    "CoQ_FSP1": ["AIFM2", "COQ2", "COQ3", "COQ5", "COQ6", "PDSS1", "PDSS2",
                 "FDFT1"],
    "NADPH": ["G6PD", "PGD", "TKT", "IDH1", "ME1", "NNT", "TALDO1"],
    "PUFA_lipid": ["ACSL4", "LPCAT3", "SCD", "FADS2", "ELOVL5", "MBOAT7"],
}


def log(*a):
    print("[204]", *a, flush=True)


# =======================================================================
# 讀資料
# =======================================================================
def load():
    log("讀 CRISPRGeneEffect.csv …")
    EFF = pd.read_csv(f"{DM}/CRISPRGeneEffect.csv", index_col=0)
    EFF.columns = [c.split(" (")[0] for c in EFF.columns]
    EFF = EFF.astype(np.float32)
    log("  gene effect：%d 株 × %d 基因" % EFF.shape)

    log("讀 OmicsExpressionProteinCodingGenesTPMLogp1.csv …")
    EXP = pd.read_csv(f"{DM}/OmicsExpressionProteinCodingGenesTPMLogp1.csv",
                      index_col=0)
    EXP.columns = [c.split(" (")[0] for c in EXP.columns]
    EXP = EXP.loc[:, ~EXP.columns.duplicated()].astype(np.float32)
    log("  表現：%d 株 × %d 基因" % EXP.shape)

    mod = pd.read_csv(f"{DM}/Model.csv", index_col=0)
    gr = pd.read_csv(f"{DM}/CRISPRInferredModelGrowthRate.csv", index_col=0)
    return EFF, EXP, mod, gr


def nrf2_alterations():
    """KEAP1 / CUL3 失能，或 NFE2L2 的非同義變異 → NRF2 軸改變（二元）。"""
    want = {"KEAP1", "NFE2L2", "CUL3"}
    keep = []
    cols = ["ModelID", "HugoSymbol", "ProteinChange", "VariantInfo",
            "VepImpact", "LikelyLoF", "Hotspot", "TumorSuppressorHighImpact"]
    log("掃 OmicsSomaticMutations.csv 取 KEAP1/NFE2L2/CUL3 …")
    for ch in pd.read_csv(f"{DM}/OmicsSomaticMutations.csv", usecols=cols,
                          chunksize=400_000, low_memory=False):
        keep.append(ch[ch.HugoSymbol.isin(want)])
    mut = pd.concat(keep, ignore_index=True)
    log("  命中 %d 個變異" % len(mut))

    def truthy(s):
        return s.astype(str).str.lower().isin({"true", "1", "yes", "t"})

    nonsyn = ~mut.VariantInfo.astype(str).str.contains(
        "SILENT|INTRON|UTR|FLANK|IGR|RNA|SPLICE_REGION", case=False, na=False)
    lof = truthy(mut.LikelyLoF) | truthy(mut.TumorSuppressorHighImpact) | \
        mut.VepImpact.astype(str).str.upper().eq("HIGH")
    ts_lof = mut.HugoSymbol.isin({"KEAP1", "CUL3"}) & lof
    nrf2_gof = mut.HugoSymbol.eq("NFE2L2") & nonsyn & \
        mut.VepImpact.astype(str).str.upper().isin({"HIGH", "MODERATE"})
    alt = mut[ts_lof | nrf2_gof]
    return (set(alt.ModelID.dropna().astype(str)),
            set(alt[alt.HugoSymbol.eq("NFE2L2")].ModelID.dropna().astype(str)),
            set(alt[alt.HugoSymbol.isin({"KEAP1", "CUL3"})]
                .ModelID.dropna().astype(str)))


# =======================================================================
# 小工具
# =======================================================================
def zrow(df):
    """逐欄 z 標準化（跨細胞株）。"""
    return (df - df.mean()) / df.std(ddof=0).replace(0, np.nan)


def panel_score(EXP, genes):
    """面板分數：各基因先跨株 z 標準化再平均，避免高表現基因主宰。"""
    g = [x for x in genes if x in EXP.columns]
    if not g:
        return None, []
    return zrow(EXP[g]).mean(axis=1), g


def spearman_cols(X, y):
    """X 的每一欄對 y 的 Spearman r（向量化，X、y 都不含 NaN）。"""
    R = X.rank()
    Rz = (R - R.mean()) / R.std(ddof=0).replace(0, np.nan)
    yz = stats.rankdata(y)
    yz = (yz - yz.mean()) / yz.std()
    return Rz.mul(yz, axis=0).mean(axis=0)


def min_detectable_r(n):
    return 2.8 / np.sqrt(max(n - 3, 1)) if n > 4 else np.nan


def verdict(obs, pred, null, n, para_p):
    """三態判決：方向對 + 超過虛無的 5/95 分位 + 參數 p < 0.05 → pass；
    這個分層有力但沒過 → fail；力氣不夠分辨 → underpowered。"""
    if pred == "none":
        return "control"
    if null is None or len(null) < 50:
        return "no_null"
    null = np.asarray(null, dtype=float)
    null = null[np.isfinite(null)]
    ok_dir = (obs < 0) if pred == "neg" else (obs > 0)
    beat = (obs <= np.quantile(null, 0.05)) if pred == "neg" \
        else (obs >= np.quantile(null, 0.95))
    if ok_dir and beat and para_p < 0.05:
        return "pass"
    # 沒過的時候要先問「這個分層分辨得出來嗎」：效應量小於 80% 檢定力下
    # 的最小可偵測 r，就只能說沒力氣，不能說否證。
    if n < 40 or abs(obs) < min_detectable_r(n):
        return "underpowered"
    return "fail"


def bh(p):
    """Benjamini-Hochberg q 值。"""
    p = np.asarray(p, dtype=float)
    o = np.argsort(p)
    q = np.empty_like(p)
    m = len(p)
    q[o] = np.minimum.accumulate((p[o] * m / np.arange(1, m + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)


def emp_p(null, obs, pred):
    null = np.asarray(null, dtype=float)
    null = null[np.isfinite(null)]
    if len(null) < 50:
        return np.nan
    if pred == "neg":
        return float((null <= obs).mean())
    if pred == "pos":
        return float((null >= obs).mean())
    return float((np.abs(null) >= abs(obs)).mean())


def fisher2(a, b, c, d):
    orr, p = stats.fisher_exact([[a, b], [c, d]], alternative="greater")
    return float(orr), float(p)


# =======================================================================
# 主程式
# =======================================================================
def main():
    EFF, EXP, mod, gr = load()
    out = {}

    # ---------- 分組：和步驟 135 一樣用 == "Lymphoid"，並檢核既有數字 ----
    lin = mod["OncotreeLineage"].reindex(EFF.index).fillna("Unknown")
    dis = mod["OncotreePrimaryDisease"].reindex(EFF.index).fillna("")
    sub = mod["OncotreeSubtype"].reindex(EFF.index).fillna("")
    form = mod["FormulationID"].reindex(EFF.index).fillna("")
    gp = mod["GrowthPattern"].reindex(EFF.index)
    is_lymph = lin.eq("Lymphoid")
    s7 = EFF[TARGET]
    dep = s7 < DEP_CUT
    assert len(EFF) == 1178, len(EFF)
    assert int(is_lymph.sum()) == 93, int(is_lymph.sum())
    assert int(dep.sum()) == 18, int(dep.sum())
    assert int((dep & is_lymph).sum()) == 12, int((dep & is_lymph).sum())
    log("既有數字檢核通過：1,178 株 / 93 淋巴系 / 18 依賴 / 12 淋巴系依賴")

    # 免疫系細分組：B 系要把漿細胞骨髓瘤分出來（它已離開 B 細胞程式）
    b_dis = {"Mature B-Cell Neoplasms", "B-Lymphoblastic Leukemia/Lymphoma"}
    is_pc = sub.eq("Plasma Cell Myeloma")
    is_bcell = is_lymph & (dis.isin(b_dis) | sub.isin(b_dis))
    is_t = is_lymph & (dis.str.contains("T and NK|T-Lymphoblastic") |
                       sub.str.contains("T and NK|T-Lymphoblastic"))
    # 注意：不能用 contains("Hodgkin")，那會把 "Non-Hodgkin Lymphoma" 也抓進來
    is_hodgkin = is_lymph & dis.eq("Hodgkin Lymphoma")
    is_bnp = is_bcell & ~is_pc                       # B 系但非漿細胞
    grp = pd.Series("NonHaematopoietic", index=EFF.index)
    grp[lin.eq("Myeloid")] = "Myeloid"
    grp[is_t] = "T_NK"
    grp[is_hodgkin] = "Hodgkin"
    grp[is_bcell & is_pc] = "B_PlasmaCell"
    grp[is_bnp] = "B_nonPlasmaCell"
    grp[is_lymph & grp.eq("NonHaematopoietic")] = "Lymphoid_other"

    # ---------- T1 收斂範圍 ---------------------------------------------
    log("\n=== T1 依賴性到底落在哪一群 ===")
    rows = []
    for g, idx in grp.groupby(grp).groups.items():
        v = s7.reindex(idx)
        rows.append({"group": g, "n": len(idx), "n_dependent": int((v < DEP_CUT).sum()),
                     "dep_rate": float((v < DEP_CUT).mean()),
                     "median_effect": float(v.median()),
                     "median_effect_RPL13A": float(
                         EFF[PAN_ESSENTIAL].reindex(idx).median())})
    gt = pd.DataFrame(rows).sort_values("dep_rate", ascending=False)
    log("\n" + gt.to_string(index=False))
    gt.to_csv(f"{RES}/204_group_summary.csv", index=False)

    t1 = {}
    # 關鍵一：完全在淋巴系內部的對比（不可能被譜系混淆）
    a = int((dep & is_bnp).sum()); b = int((~dep & is_bnp).sum())
    other_ly = is_lymph & ~is_bnp
    c = int((dep & other_ly).sum()); d = int((~dep & other_ly).sum())
    orr, p = fisher2(a, c, b, d)
    t1["within_lymphoid_Bnp_vs_otherLymphoid"] = {
        "n_Bnp": a + b, "dep_Bnp": a, "n_other_lymphoid": c + d,
        "dep_other_lymphoid": c, "odds_ratio": orr, "fisher_p": p}
    log("淋巴系內部：B（非漿細胞）%d/%d 依賴 vs 其他淋巴系 %d/%d"
        "  OR = %s，Fisher p = %.2e" % (a, a + b, c, c + d,
                                        "inf" if np.isinf(orr) else "%.1f" % orr, p))
    # 關鍵二：B 系內部，漿細胞 vs 非漿細胞
    a2 = a; b2 = b
    c2 = int((dep & is_bcell & is_pc).sum()); d2 = int((~dep & is_bcell & is_pc).sum())
    orr2, p2 = fisher2(a2, c2, b2, d2)
    t1["within_Bcell_nonPC_vs_plasmaCell"] = {
        "dep_nonPC": a2, "n_nonPC": a2 + b2, "dep_PC": c2, "n_PC": c2 + d2,
        "odds_ratio": orr2, "fisher_p": p2,
        "median_effect_nonPC": float(s7[is_bnp].median()),
        "median_effect_PC": float(s7[is_bcell & is_pc].median())}
    log("B 系內部：非漿細胞 %d/%d vs 漿細胞骨髓瘤 %d/%d  Fisher p = %.4f"
        % (a2, a2 + b2, c2, c2 + d2, p2))
    # T/NK、Hodgkin、髓系（培養條件接近的對照組）
    for nm, m in (("T_NK", is_t), ("Hodgkin", is_hodgkin),
                  ("Myeloid", lin.eq("Myeloid"))):
        t1[f"control_group_{nm}"] = {
            "n": int(m.sum()), "n_dependent": int((dep & m).sum()),
            "median_effect": float(s7[m].median())}
    # 虛無：依賴率配對隨機基因的「依賴株中 B（非漿細胞）比例」
    # D 是「這株依賴這個基因」的布林矩陣，預先算好，後面所有配對虛無共用
    D = ((EFF < DEP_CUT) & EFF.notna()).to_numpy()
    gpos = {g: i for i, g in enumerate(EFF.columns)}
    dep_frac_all = pd.Series(D.mean(0), index=EFF.columns)
    bnp_arr = is_bnp.to_numpy()

    def matched_frac_null(gene, n_draw=N_NULL, rows=None, group=None,
                          min_dep=3):
        """依賴率配對隨機基因的「依賴株中屬於某一組的比例」分佈。

        rows 限定只看哪些株（例如只看淋巴系、或只看同一種培養液的株），
        配對依賴率也在同一個子集裡算，所以子集內的觀察值與虛無是同一個尺度。
        回傳 (虛無分佈, 觀察值)。"""
        if group is None:
            group = bnp_arr
        rmask = np.ones(D.shape[0], dtype=bool) if rows is None \
            else np.asarray(rows, dtype=bool)
        Dr = D[rmask]
        grp_r = np.asarray(group, dtype=bool)[rmask]
        frac_sub = Dr.mean(0)
        f = float(frac_sub[gpos[gene]])
        obs_d = Dr[:, gpos[gene]]
        obs = float((obs_d & grp_r).sum() / max(obs_d.sum(), 1))
        cand = np.where((frac_sub > f * 0.6) & (frac_sub < f * 1.6)
                        & (frac_sub > 0))[0]
        cand = cand[cand != gpos[gene]]
        if len(cand) < 50:
            return np.array([]), obs
        sel = RNG.choice(cand, size=min(n_draw, len(cand)), replace=False)
        sub = Dr[:, sel]
        tot = sub.sum(0)
        hit = (sub & grp_r[:, None]).sum(0)
        keep = tot >= min_dep
        return hit[keep] / tot[keep], obs

    nullfrac, obsfrac = matched_frac_null(TARGET)
    t1["empirical_null_frac_Bnp_among_dependents"] = {
        "observed": obsfrac, "n_matched_genes": len(nullfrac),
        "null_mean": float(np.mean(nullfrac)),
        "null_q95": float(np.quantile(nullfrac, 0.95)),
        "empirical_p": float((np.array(nullfrac) >= obsfrac).mean())}
    log("依賴株中 B（非漿細胞）比例 = %.3f；%d 個依賴率配對基因 %.3f ± %.3f"
        "（95 分位 %.3f）→ 經驗 p = %.4f"
        % (obsfrac, len(nullfrac), np.mean(nullfrac), np.std(nullfrac),
           np.quantile(nullfrac, 0.95),
           t1["empirical_null_frac_Bnp_among_dependents"]["empirical_p"]))
    # 三個頭條對比，各自要有自己的經驗虛無（Fisher p 不算）：
    # 配對依賴率都在同一個子集裡算，所以觀察值與虛無同尺度。
    log("\n三個頭條對比的經驗虛無（依賴率配對隨機基因，子集內計算）：")
    for nm, rows_m, desc in (
            ("within_lymphoid", is_lymph.to_numpy(),
             "只看 93 株淋巴系：依賴者中 B（非漿細胞）的比例"),
            ("within_Bcell", is_bcell.to_numpy(),
             "只看 72 株 B 系：依賴者中非漿細胞的比例"),
            ("rpmi10_matched", (form.eq("RPMI + 10% FBS") &
                                (is_bnp | ~is_lymph)).to_numpy(),
             "只看 RPMI+10%FBS 的 B（非漿細胞）與非淋巴系株")):
        nf, ob = matched_frac_null(TARGET, rows=rows_m)
        ep = float((nf >= ob).mean()) if len(nf) else np.nan
        t1[f"empirical_null_{nm}"] = {
            "description": desc, "observed_frac": ob,
            "n_matched_genes": int(len(nf)),
            "null_mean": float(np.mean(nf)) if len(nf) else np.nan,
            "null_q95": float(np.quantile(nf, 0.95)) if len(nf) else np.nan,
            "empirical_p": ep}
        log("  %-16s 觀察 %.3f vs 虛無 %.3f（95 分位 %.3f，%d 個基因）→ 經驗 p = %s"
            % (nm, ob, np.mean(nf) if len(nf) else np.nan,
               np.quantile(nf, 0.95) if len(nf) else np.nan, len(nf),
               "%.4f" % ep if ep == ep else "-"))

    # 陰性對照：RPL13A 不該有組間結構
    r13 = EFF[PAN_ESSENTIAL]

    def group_kruskal(v):
        """組間差異；只納入該基因有 ≥ 3 個非缺值的組（Hodgkin 的 RPL13A 全缺）。"""
        gs = [v[grp == g].dropna() for g in sorted(grp.unique())]
        gs = [x for x in gs if len(x) >= 3]
        return stats.kruskal(*gs), len(gs)

    kw, n13 = group_kruskal(r13)
    kw7, n7 = group_kruskal(s7)
    t1["negative_control_RPL13A_group_kruskal_p"] = float(kw.pvalue)
    t1["negative_control_RPL13A_n_groups"] = n13
    t1["SLC7A11_group_kruskal_p"] = float(kw7.pvalue)
    log("組間差異 Kruskal：SLC7A11 p = %.2e；陰性對照 RPL13A p = %.3f"
        % (kw7.pvalue, kw.pvalue))
    # 那 6 株非淋巴系依賴者有沒有 B 標記（可能是標錯譜系）
    nl_dep = EFF.index[(dep & ~is_lymph).to_numpy()]
    bmk = [g for g in ["CD79A", "CD19", "PAX5", "MS4A1"] if g in EXP.columns]
    t1["non_lymphoid_dependents_Bmarkers"] = []
    log("\n非淋巴系的 %d 株依賴者，B 標記表現（全庫 90 百分位當基準）：" % len(nl_dep))
    q90 = EXP[bmk].quantile(0.90)
    for i in nl_dep:
        e = EXP[bmk].reindex([i]).iloc[0] if i in EXP.index else None
        rec = {"model": str(i),
               "line": str(mod.loc[i, "StrippedCellLineName"]) if i in mod.index else "",
               "lineage": str(lin[i]), "effect": float(s7[i])}
        if e is not None and e.notna().any():
            rec.update({f"exp_{g}": float(e[g]) for g in bmk})
            rec["any_B_marker_above_q90"] = bool((e > q90).any())
        t1["non_lymphoid_dependents_Bmarkers"].append(rec)
        log("  %-14s %-24s effect %+.2f  %s" % (
            rec["line"], rec["lineage"], rec["effect"],
            " ".join("%s=%.2f" % (g, rec.get("exp_" + g, np.nan)) for g in bmk)))
    out["T1_scope"] = t1

    # ---------- 建特徵表 ------------------------------------------------
    log("\n=== 組裝細胞株層級特徵 ===")
    shared = EFF.index.intersection(EXP.index)
    log("gene effect ∩ 表現 = %d 株（18 株依賴者保留 %d 株，12 株 B 系依賴者保留 %d 株）"
        % (len(shared), int(dep.reindex(shared).sum()),
           int((dep & is_bnp).reindex(shared).sum())))
    F = pd.DataFrame(index=EFF.index)
    panel_used = {}
    for k, (genes, pred, desc) in EXPR_PANELS.items():
        sc, used = panel_score(EXP, genes)
        if sc is None:
            log("  !! %s 一個基因都找不到，跳過" % k)
            continue
        F[k] = sc.reindex(EFF.index)
        panel_used[k] = used
        if len(used) != len(genes):
            log("  %s：%d/%d 個基因有表現資料（缺 %s）"
                % (k, len(used), len(genes),
                   ",".join(sorted(set(genes) - set(used)))))

    # T6 的技術特徵
    grz = pd.DataFrame(index=EFF.index)
    for c in gr.columns:                     # 三個 library 欄位互斥，各自 z
        v = gr[c].reindex(EFF.index)
        grz[c] = (v - v.mean()) / v.std(ddof=0)
    F["growth_rate_z"] = grz.mean(axis=1)    # 每株幾乎只有一欄非 NaN
    lib = gr.reindex(EFF.index).notna().idxmax(axis=1).where(
        gr.reindex(EFF.index).notna().any(axis=1))
    F["screen_lib_avana"] = lib.eq("Achilles-Avana-2D").astype(float)
    F["is_suspension"] = gp.eq("Suspension").astype(float)
    F["is_rpmi"] = form.str.contains("RPMI").astype(float)
    F["line_mean_effect"] = EFF.mean(axis=1)
    nrf2_any, nrf2_gof, keap_lof = nrf2_alterations()
    F["nrf2_pathway_alt"] = pd.Series(
        [1.0 if str(i) in nrf2_any else 0.0 for i in EFF.index], index=EFF.index)
    LINE_FEATS = {"growth_rate_z": ("neg", "T6a 生長速率：越快越依賴"),
                  "screen_lib_avana": ("none", "T6b 篩選 library（技術）"),
                  "is_suspension": ("neg", "T6c 懸浮培養"),
                  "is_rpmi": ("neg", "T6d RPMI 培養液"),
                  "line_mean_effect": ("pos", "T6e 株層級整體篩選敏感度"),
                  "nrf2_pathway_alt": ("neg", "T3 KEAP1/NFE2L2/CUL3 改變")}
    log("NRF2 軸有改變的株數：%d（NFE2L2 %d / KEAP1·CUL3 %d）"
        % (int(F.nrf2_pathway_alt.sum()), len(nrf2_gof), len(keap_lof)))

    # 結果變數
    Y = pd.DataFrame(index=EFF.index)
    Y["SLC7A11"] = s7
    Y["GCLC"] = EFF["GCLC"]
    Y["module_z"] = zrow(EFF[[g for g in MODULE if g in EFF.columns]]).mean(axis=1)
    Y["RPL13A_negctrl"] = EFF[PAN_ESSENTIAL]
    Y["Bidentity_negctrl"] = zrow(EFF[[g for g in BIDENT
                                       if g in EFF.columns]]).mean(axis=1)

    # 第四個分層 B_nonPlasmaCell 是必要的：B 系內部通過的候選，有可能只是
    # 「是不是漿細胞骨髓瘤」的代名詞（骨髓瘤 0/19 依賴），把骨髓瘤去掉才知道
    # 特徵是不是真的在 12 株依賴者與 41 株非依賴的 B 細胞株之間有差別。
    STRATA = {"all": pd.Series(True, index=EFF.index),
              "non_lymphoid": ~is_lymph,
              "B_lineage": is_bcell,
              "B_nonPlasmaCell": is_bnp,
              # 髓系株的 AIFM2 表現比 B 系還低卻幾乎不依賴，所以非造血系
              # 這個分層（同時排掉淋巴與髓系）是 T4 必要的對照
              "non_haematopoietic": grp.eq("NonHaematopoietic")}
    for k, m in STRATA.items():
        log("分層 %-13s n = %-5d 依賴者 %-3d  80%% 檢定力最小可偵測 r ≈ %.3f"
            % (k, int(m.sum()), int((dep & m).sum()),
               min_detectable_r(int(m.sum()))))
    out["strata"] = {k: {"n": int(m.sum()), "n_dependent": int((dep & m).sum()),
                         "min_detectable_r_80pct": float(min_detectable_r(int(m.sum())))}
                     for k, m in STRATA.items()}

    # ---------- 虛無分佈 A：隨機基因 / 隨機面板的表現 -------------------
    log("\n=== 建 Null A（隨機基團表現）與 Null B（依賴率配對基因）===")
    expr_pool = EXP.loc[:, EXP.notna().mean() > 0.95]
    expr_pool = expr_pool.loc[:, expr_pool.std() > 0.1]
    Zpool = zrow(expr_pool).fillna(0.0)      # 已跨株 z 化，缺值補 0 ＝ 補平均
    A_pool = Zpool.to_numpy()
    sizes = sorted({1} | {len(v) for v in panel_used.values()})
    nullA = {}        # (outcome, stratum, k) -> array
    for st, mask in STRATA.items():
        idx = EFF.index[mask.to_numpy()].intersection(Zpool.index)
        pos = Zpool.index.get_indexer(idx)
        for k in sizes:
            cols = RNG.integers(0, A_pool.shape[1], size=(N_NULL, k))
            S = np.empty((len(idx), N_NULL), dtype=np.float32)
            Asub = A_pool[pos]
            for j in range(N_NULL):
                S[:, j] = Asub[:, cols[j]].mean(axis=1)
            Sdf = pd.DataFrame(S, index=idx)
            for yname in Y.columns:
                y = Y[yname].reindex(idx)
                ok = y.notna().to_numpy()
                nullA[(yname, st, k)] = spearman_cols(
                    Sdf[ok], y[ok]).to_numpy()
        log("  Null A 完成分層 %s（n = %d）" % (st, len(idx)))

    # ---------- 虛無分佈 B：依賴率配對基因的 gene effect ----------------
    nullB_mat = {}
    for yname, gene in (("SLC7A11", "SLC7A11"), ("GCLC", "GCLC")):
        f = float((EFF[gene] < DEP_CUT).mean())
        nr = dep_frac_all[(dep_frac_all > f * 0.6) & (dep_frac_all < f * 1.6)
                          & (dep_frac_all > 0)].index
        nr = [x for x in nr if x != gene]
        pick = list(RNG.choice(nr, size=min(N_NULL, len(nr)), replace=False))
        M = EFF[pick]
        nullB_mat[yname] = M.fillna(M.mean())
        log("  Null B：%s 的依賴率 %.4f，配對到 %d 個基因"
            % (yname, f, len(pick)))

    # ---------- 逐候選 × 逐分層 × 逐結果變數 ----------------------------
    log("\n=== 候選解釋逐項檢定 ===")
    ALL_FEATS = {k: (EXPR_PANELS[k][1], EXPR_PANELS[k][2]) for k in panel_used}
    ALL_FEATS.update(LINE_FEATS)
    rows = []
    for feat, (pred, desc) in ALL_FEATS.items():
        k = len(panel_used.get(feat, [1]))
        for st, mask in STRATA.items():
            idx = EFF.index[mask.to_numpy()]
            for yname in Y.columns:
                xs = F[feat].reindex(idx)
                ys = Y[yname].reindex(idx)
                ok = xs.notna() & ys.notna()
                n = int(ok.sum())
                if n < 20 or xs[ok].std() == 0:
                    continue
                r, p = stats.spearmanr(xs[ok], ys[ok])
                nA = nullA.get((yname, st, k))
                pA = emp_p(nA, r, pred) if nA is not None else np.nan
                # Null B：同一特徵 vs 依賴率配對基因的 effect
                pB, nb = np.nan, None
                if yname in nullB_mat:
                    M = nullB_mat[yname].reindex(idx)[ok.to_numpy()]
                    nb = spearman_cols(M, xs[ok]).to_numpy()
                    pB = emp_p(nb, r, pred)
                # 表現類特徵問的是「這個基因的表現特別嗎」→ Null A；
                # 培養、生長速率、突變這些株層級特徵沒有「隨機基因的表現」
                # 可以比，要問的是「SLC7A11 在依賴率相近的基因裡特別嗎」→ Null B。
                nuse = nA if feat in panel_used else nb
                vd = verdict(r, pred, nuse, n, p)
                rows.append({"feature": feat, "description": desc,
                             "predicted_sign": pred, "outcome": yname,
                             "stratum": st, "n": n, "spearman_r": float(r),
                             "parametric_p": float(p),
                             "empirical_p_nullA_random_genes": pA,
                             "empirical_p_nullB_ratematched_genes": pB,
                             "null_used_for_verdict": "A" if feat in panel_used
                             else "B",
                             "min_detectable_r_80pct": float(min_detectable_r(n)),
                             "verdict": vd})
        pr = [x for x in rows if x["feature"] == feat
              and x["outcome"] == "SLC7A11"]

        def _f(v):
            return "%.3f" % v if v == v else "-"
        log("%-22s " % feat + " | ".join(
            "%s r=%+.3f pA=%s pB=%s %s"
            % (x["stratum"][:9], x["spearman_r"],
               _f(x["empirical_p_nullA_random_genes"]),
               _f(x["empirical_p_nullB_ratematched_genes"]),
               x["verdict"]) for x in pr))
    CAND = pd.DataFrame(rows)
    CAND.to_csv(f"{RES}/204_candidate_ranking.csv", index=False)

    # ---------- 面板逐基因分解：哪一個基因在推動面板分數 ----------------
    # 面板裡挑最好的那個基因等於偷看，所以虛無用 max 統計量：
    # 從「隨機單基因」的虛無分佈抽 k 個，取預測方向上最極端的那個。
    log("\n=== 面板逐基因分解（最佳基因對 max 統計量虛無）===")
    drows = []
    for feat, (genes, pred, desc) in EXPR_PANELS.items():
        used = panel_used.get(feat, [])
        if len(used) < 2 or pred == "none":
            continue
        for st, mask in STRATA.items():
            idx = EFF.index[mask.to_numpy()]
            y = Y[TARGET].reindex(idx)
            base = nullA.get((TARGET, st, 1))
            if base is None:
                continue
            draws = RNG.choice(np.asarray(base, dtype=float),
                               size=(N_NULL, len(used)))
            ext = draws.min(1) if pred == "neg" else draws.max(1)
            rs = {}
            for g in used:
                x = EXP[g].reindex(idx)
                o = x.notna() & y.notna()
                if o.sum() < 20 or x[o].std() == 0:
                    continue
                r, p = stats.spearmanr(x[o], y[o])
                rs[g] = r
                drows.append({"panel": feat, "gene": g, "stratum": st,
                              "n": int(o.sum()), "spearman_r": float(r),
                              "parametric_p": float(p),
                              "empirical_p_single_gene_null": emp_p(base, r, pred),
                              "is_panel_best": False})
            if not rs:
                continue
            bg = min(rs, key=rs.get) if pred == "neg" else max(rs, key=rs.get)
            pm = float((ext <= rs[bg]).mean()) if pred == "neg" \
                else float((ext >= rs[bg]).mean())
            for d in drows:
                if d["panel"] == feat and d["stratum"] == st and d["gene"] == bg:
                    d["is_panel_best"] = True
                    d["empirical_p_maxstat_null"] = pm
            log("  %-20s %-16s 最佳基因 %-8s r = %+.3f  max 統計量經驗 p = %.4f"
                % (feat, st, bg, rs[bg], pm))
    DEC = pd.DataFrame(drows)
    DEC.to_csv(f"{RES}/204_panel_decomposition.csv", index=False)
    out["panel_decomposition_best"] = DEC[DEC.is_panel_best].to_dict("records")

    # ---------- (B) 衰減檢定 -------------------------------------------
    log("\n=== 衰減檢定：加入特徵後，is_B 的係數縮掉多少 ===")
    isB = is_bnp.astype(float)        # 用「B 系非漿細胞」當要被解釋的指標
    yv = Y["SLC7A11"]

    att = {}
    bB_arr = isB.to_numpy()
    # 衰減量本身不看方向：關係跑反了也會讓 is_B 的係數縮小。
    # 所以每一列都要標記「方向和預設預測一致嗎」，方向反的不能算通過。
    rall = {r0["feature"]: (r0["spearman_r"], r0["predicted_sign"])
            for r0 in rows if r0["outcome"] == "SLC7A11"
            and r0["stratum"] == "all"}
    for feat in ALL_FEATS:
        x = F[feat]
        # 純 is_B 的係數與加了特徵之後的係數，都限制在同一組株上才可比
        ok = yv.notna() & x.notna()
        if ok.sum() < 50 or x[ok].std() == 0:
            continue
        yk = yv[ok].to_numpy()
        bk = isB[ok].to_numpy()
        v = x[ok].to_numpy()
        v = (v - v.mean()) / (v.std() or 1)
        rawm = float(np.linalg.lstsq(
            np.vstack([np.ones(len(yk)), bk]).T, yk, rcond=None)[0][1])
        adj = float(np.linalg.lstsq(
            np.vstack([np.ones(len(yk)), bk, v]).T, yk, rcond=None)[0][1])
        shrink = 1 - adj / rawm if rawm != 0 else np.nan
        # 虛無：B 系特異度相當的隨機基因（標準化均值差配對 ±25%）
        pe = np.nan
        okb = bB_arr[ok.to_numpy()] == 1
        dobs = abs(v[okb].mean() - v[~okb].mean())      # v 已標準化
        idxo = Zpool.index.intersection(EFF.index[ok.to_numpy()])
        yo = yv.reindex(idxo).to_numpy()
        bo = isB.reindex(idxo).to_numpy()
        A = A_pool[Zpool.index.get_indexer(idxo)]
        d_all = np.abs(A[bo == 1].mean(axis=0) - A[bo == 0].mean(axis=0)) / \
            (A.std(axis=0) + 1e-9)
        cand = np.where((d_all > 0.75 * dobs) & (d_all < 1.25 * dobs))[0]
        if len(cand) >= 30:
            pick = RNG.choice(cand, size=min(600, len(cand)), replace=False)
            raw0 = np.linalg.lstsq(np.vstack([np.ones(len(yo)), bo]).T, yo,
                                   rcond=None)[0][1]
            sh = []
            for j in pick:
                w = A[:, j]
                w = (w - w.mean()) / (w.std() + 1e-9)
                bj = np.linalg.lstsq(np.vstack([np.ones(len(yo)), bo, w]).T,
                                     yo, rcond=None)[0][1]
                sh.append(1 - bj / raw0)
            pe = float((np.array(sh) >= shrink).mean())
        robs, pr_sign = rall.get(feat, (np.nan, "none"))
        wrong = bool(pr_sign in ("neg", "pos") and robs == robs and
                     ((pr_sign == "neg" and robs > 0) or
                      (pr_sign == "pos" and robs < 0)))
        att[feat] = {"beta_isB_raw": rawm, "beta_isB_adjusted": adj,
                     "shrinkage": float(shrink), "n": int(ok.sum()),
                     "B_specificity_std_diff": float(dobs),
                     "n_matched_null_genes": int(len(cand)),
                     "null_median_shrinkage": float(np.median(sh))
                     if pe == pe else np.nan,
                     "empirical_p": pe,
                     "all_lines_r": float(robs) if robs == robs else None,
                     "attenuation_wrong_sign_all_lines": wrong,
                     "wrong_sign_vs_prediction": wrong}
        log("%-22s beta_isB %+.4f → %+.4f（縮 %5.1f%%）經驗 p = %s%s"
            % (feat, rawm, adj, 100 * shrink,
               "%.3f" % pe if pe == pe else "-",
               "   ← 方向與預測相反，不算支持" if wrong else ""))
    out["attenuation_of_isB"] = att

    # ---------- T5 氧化還原各臂的 B 系富集 ------------------------------
    log("\n=== T5 依賴性富集攤到各條氧化還原臂 ===")
    rows = []
    for arm, genes in REDOX_ARMS.items():
        for g in genes:
            if g not in EFF.columns:
                continue
            d = D[:, gpos[g]]
            if d.sum() < 1:
                rows.append({"arm": arm, "gene": g, "n_dependent": 0,
                             "n_dep_Bnp": 0, "frac_Bnp": np.nan,
                             "odds_ratio": np.nan, "fisher_p": np.nan,
                             "empirical_p": np.nan, "dep_frac": 0.0})
                continue
            aa = int((d & bnp_arr).sum()); bb = int((d & ~bnp_arr).sum())
            cc = int((~d & bnp_arr).sum()); dd2 = int((~d & ~bnp_arr).sum())
            orr, pv = fisher2(aa, bb, cc, dd2)
            e, _ = matched_frac_null(g, 1000)
            obs = aa / max(aa + bb, 1)
            rows.append({"arm": arm, "gene": g, "n_dependent": aa + bb,
                         "n_dep_Bnp": aa, "frac_Bnp": obs,
                         "odds_ratio": orr, "fisher_p": pv,
                         "empirical_p": float((e >= obs).mean())
                         if len(e) else np.nan,
                         "dep_frac": float(dep_frac_all[g])})
    RX = pd.DataFrame(rows)
    RX.to_csv(f"{RES}/204_redox_dependency_map.csv", index=False)
    armsum = pd.DataFrame([{
        "arm": arm,
        "n_genes": len(x),
        "n_genes_with_dependents": int((x.n_dependent > 0).sum()),
        "n_genes_emp_p_lt_005": int((x.empirical_p < 0.05).sum()),
        "median_frac_Bnp": float(x.frac_Bnp.median()),
        "genes_enriched": ",".join(x.gene[x.empirical_p < 0.05]),
    } for arm, x in RX.groupby("arm")]).set_index("arm")
    log("\n" + armsum.to_string())
    out["T5_redox_arms"] = json.loads(armsum.to_json(orient="index"))
    out["limitations"] = {
        "RPL13A_outcome_coverage": int(EFF[PAN_ESSENTIAL].notna().sum()),
        "RPL13A_note": "RPL13A 只有這麼多株有 gene effect（全庫 1,178），"
                       "所以它當結果變數的陰性對照檢定力偏低。",
        "PRDX6_arm_misassignment": "PRDX6 是用穀胱甘肽的 1-Cys peroxiredoxin，"
                                   "放在 thioredoxin 臂是跑前就寫死的歸類錯誤；"
                                   "它富集與『只有 GSH 這條路被需要』一致，"
                                   "不是硫氧還蛋白臂的反例。",
        "culture_not_fully_excluded": "Model.csv 沒有硫醇類添加物與血清批次。",
        "T5_circularity": "SLC7A11/GCLC/GCLM/AIFM2 是步驟 135 的輸入；"
                          "GSS/PRDX6/ACSL4/SCLY/SEPHS2/PDSS1/FDFT1 是新的。"}
    out["T5_caveat"] = (
        "SLC7A11/GCLC/GCLM/AIFM2 是步驟 135 的輸入，它們在這裡富集不是獨立確認；"
        "共檢定 %d 個基因，未校正下期望約 %.1f 個偽陽性。"
        "『硫氧還蛋白臂與 NADPH 臂沒有富集』只對有足夠依賴株可檢定的基因成立，"
        "每個基因的 n_dependent 都列在 204_redox_dependency_map.csv。"
        % (len(RX), 0.05 * RX.empirical_p.notna().sum()))
    out["T5_genes_enriched"] = RX[RX.empirical_p < 0.05][
        ["arm", "gene", "n_dependent", "n_dep_Bnp", "frac_Bnp",
         "odds_ratio", "empirical_p"]].to_dict("records")

    # ---------- T6 培養條件：同配方內的對比 -----------------------------
    log("\n=== T6 同一種培養液內的對比（RPMI + 10% FBS）===")
    t6 = {}
    rp = form.eq("RPMI + 10% FBS")
    for nm, m in (("B_nonPlasmaCell", is_bnp), ("B_PlasmaCell", is_bcell & is_pc),
                  ("T_NK", is_t), ("Myeloid", lin.eq("Myeloid")),
                  ("NonLymphoid", ~is_lymph)):
        mm = rp & m
        t6[f"rpmi10_{nm}"] = {"n": int(mm.sum()),
                              "n_dependent": int((dep & mm).sum()),
                              "median_effect": float(s7[mm].median())
                              if mm.sum() else np.nan}
        log("  RPMI+10%%FBS 內 %-18s n = %-4d 依賴 %-3d 中位 effect %s"
            % (nm, int(mm.sum()), int((dep & mm).sum()),
               "%+.3f" % s7[mm].median() if mm.sum() else "-"))
    aa = int((dep & rp & is_bnp).sum()); bb = int((~dep & rp & is_bnp).sum())
    cc = int((dep & rp & ~is_lymph).sum()); dd2 = int((~dep & rp & ~is_lymph).sum())
    orr, pv = fisher2(aa, cc, bb, dd2)
    t6["rpmi_matched_Bnp_vs_nonlymphoid"] = {
        "dep_Bnp": aa, "n_Bnp": aa + bb, "dep_nonlymphoid": cc,
        "n_nonlymphoid": cc + dd2, "odds_ratio": orr, "fisher_p": pv}
    log("  同配方對比：B（非漿細胞）%d/%d vs 非淋巴系 %d/%d  OR %.1f，p = %.2e"
        % (aa, aa + bb, cc, cc + dd2, orr, pv))
    # 懸浮 vs 貼附，限非淋巴系（技術性解釋的直接檢定）
    susp = gp.eq("Suspension")
    for nm, m in (("non_lymphoid_suspension", (~is_lymph) & susp),
                  ("non_lymphoid_adherent", (~is_lymph) & gp.eq("Adherent"))):
        t6[nm] = {"n": int(m.sum()), "n_dependent": int((dep & m).sum()),
                  "dep_rate": float(dep[m].mean()) if m.sum() else np.nan}
        log("  %-26s n = %-4d 依賴 %-2d（%.2f%%）"
            % (nm, int(m.sum()), int((dep & m).sum()), 100 * dep[m].mean()))
    t6["B_nonPlasmaCell_dep_rate_for_reference"] = float(dep[is_bnp].mean())
    t6["caveat"] = ("Model.csv 讀不到硫醇類添加物（2-ME、麩胱甘肽、"
                    "N-乙醯半胱胺酸）與血清批次，所以「培養條件排除」"
                    "只涵蓋生長型態與 FormulationID 這兩項。")
    out["T6_culture"] = t6

    # ---------- T3 NRF2 的不對稱預測 ------------------------------------
    log("\n=== T3 NRF2 軸 ===")
    t3 = {}
    na = F["nrf2_pathway_alt"].eq(1)
    aa = int((dep & na).sum()); bb = int((~dep & na).sum())
    cc = int((dep & ~na).sum()); dd2 = int((~dep & ~na).sum())
    orr, pv = fisher2(aa, cc, bb, dd2)
    t3["alteration_vs_dependency_all"] = {"n_altered": int(na.sum()),
                                          "dep_altered": aa, "dep_wt": cc,
                                          "odds_ratio": orr, "fisher_p": pv}
    nlm = ~is_lymph
    aa2 = int((dep & na & nlm).sum()); bb2 = int((~dep & na & nlm).sum())
    cc2 = int((dep & ~na & nlm).sum()); dd3 = int((~dep & ~na & nlm).sum())
    orr2, pv2 = fisher2(aa2, cc2, bb2, dd3)
    t3["alteration_vs_dependency_non_lymphoid"] = {
        "n_altered": int((na & nlm).sum()), "dep_altered": aa2, "dep_wt": cc2,
        "odds_ratio": orr2, "fisher_p": pv2}
    t3["n_altered_in_Bnp"] = int((na & is_bnp).sum())
    t3["nrf2_score_Bnp_vs_nonlymphoid"] = {
        "median_Bnp": float(F.nrf2_target_score[is_bnp].median()),
        "median_non_lymphoid": float(F.nrf2_target_score[nlm].median()),
        "mannwhitney_p_Bnp_higher": float(stats.mannwhitneyu(
            F.nrf2_target_score[is_bnp].dropna(),
            F.nrf2_target_score[nlm].dropna(), alternative="greater").pvalue)}
    log("  NRF2 軸改變：全庫 OR %.1f（p = %.3g）；非淋巴系 OR %.1f（p = %.3g）；"
        "B 系裡有改變的 %d 株" % (orr, pv, orr2, pv2, t3["n_altered_in_Bnp"]))
    log("  NRF2 標的分數中位：B（非漿細胞）%.3f vs 非淋巴系 %.3f（B 更高 p = %.3f）"
        % (t3["nrf2_score_Bnp_vs_nonlymphoid"]["median_Bnp"],
           t3["nrf2_score_Bnp_vs_nonlymphoid"]["median_non_lymphoid"],
           t3["nrf2_score_Bnp_vs_nonlymphoid"]["mannwhitney_p_Bnp_higher"]))
    t3["caveat"] = ("KEAP1 的旗標只抓 LikelyLoF / 高衝擊 / 腫瘤抑制高衝擊，"
                    "錯過一般錯義突變；而非淋巴系只有 6 株依賴者，"
                    "所以『NRF2 連非淋巴系那幾株也解釋不了』這半邊是"
                    "檢定力不足，不是否證。B 系那半邊（53 株裡只有 1 株有改變、"
                    "NRF2 標的分數顯著偏低）才是成立的否證。")
    out["T3_nrf2"] = t3

    # ---------- 各特徵在分組間的分佈（判斷「有沒有把 B 放在極端」）-----
    log("\n=== 每個候選特徵在各組的中位數（要當機制，B 必須在極端）===")
    gmed = F.groupby(grp).median().T
    log("\n" + gmed.round(3).to_string())
    gmed.to_csv(f"{RES}/204_feature_by_group.csv")
    out["feature_group_medians"] = json.loads(gmed.round(4).to_json())

    # 面板成員基因各自的組間中位數（z 分數），用來看是哪個基因把 B 放到極端
    pg = sorted({g for v in panel_used.values() for g in v} |
                {g for v in REDOX_ARMS.values() for g in v if g in EXP.columns})
    GZ = zrow(EXP[[g for g in pg if g in EXP.columns]]).reindex(EFF.index)
    gene_med = GZ.groupby(grp).median().T
    gene_med.to_csv(f"{RES}/204_gene_expression_by_group.csv")
    log("\n（面板成員基因的組間中位數已寫入 204_gene_expression_by_group.csv）")

    # ---------- T8 無偏掃描（探索性）-----------------------------------
    log("\n=== T8 無偏掃描（探索性，不計入通過）===")
    t8 = {"NOTE": "exploratory: hypotheses only, not passes"}
    scr = {}
    for st, mask in (("B_lineage", is_bcell), ("B_nonPlasmaCell", is_bnp),
                     ("non_lymphoid", ~is_lymph)):
        idx = EFF.index[mask.to_numpy()].intersection(Zpool.index)
        y = Y["SLC7A11"].reindex(idx)
        ok = y.notna()
        r = spearman_cols(Zpool.loc[idx][ok.to_numpy()], y[ok])
        n = int(ok.sum())
        tt = r * np.sqrt((n - 2) / np.clip(1 - r ** 2, 1e-12, None))
        pv = pd.Series(2 * stats.t.sf(np.abs(tt), n - 2), index=r.index)
        d = pd.DataFrame({"gene": r.index, "spearman_r": r.values,
                          "p": pv.values, "n": n})
        d["bh_q"] = bh(d.p.values)
        d = d.sort_values("spearman_r")
        scr[st] = d
        d.head(40).to_csv(f"{RES}/204_screen_expr_{st}_top_dep.csv", index=False)
        d.tail(40).iloc[::-1].to_csv(
            f"{RES}/204_screen_expr_{st}_top_protect.csv", index=False)
        log("  表現掃描 %s（n = %d）最「越高越依賴」的前 12 個："
            % (st, n))
        log("    " + ", ".join("%s(%+.2f)" % (a, b) for a, b in
                               zip(d.gene.head(12), d.spearman_r.head(12))))
        t8[f"expr_screen_{st}_top_dependency"] = d.head(25).to_dict("records")
        t8[f"expr_screen_{st}_top_protective"] = \
            d.tail(25).iloc[::-1].to_dict("records")

    # 共依賴
    for st, mask in (("B_lineage", is_bcell), ("non_lymphoid", ~is_lymph)):
        idx = EFF.index[mask.to_numpy()]
        E = EFF.loc[idx]
        E = E.loc[:, E.notna().mean() > 0.9]
        E = E.loc[:, E.std() > 0.03].fillna(E.mean())
        c = E.corrwith(E[TARGET], method="spearman").drop(TARGET, errors="ignore")
        c = c.sort_values(ascending=False)
        pd.DataFrame({"gene": c.index[:60], "r": c.values[:60]}).to_csv(
            f"{RES}/204_codependency_{st}.csv", index=False)
        log("  共依賴 %s 前 12：%s" % (st, ", ".join(
            "%s(%.2f)" % (a, b) for a, b in zip(c.index[:12], c.values[:12]))))
        t8[f"codependency_{st}_top"] = [
            {"gene": str(a), "r": float(b)} for a, b in
            zip(c.index[:25], c.values[:25])]

    # B 系內部：12 株依賴者 vs 其餘
    idx = EFF.index[is_bcell.to_numpy()].intersection(Zpool.index)
    dd = dep.reindex(idx)
    A = Zpool.loc[idx].to_numpy()
    g1, g0 = A[dd.to_numpy()], A[~dd.to_numpy()]
    u = np.array([stats.mannwhitneyu(g1[:, j], g0[:, j]).pvalue
                  for j in range(A.shape[1])])
    dlt = np.nanmean(g1, axis=0) - np.nanmean(g0, axis=0)
    de = pd.DataFrame({"gene": Zpool.columns, "delta_z": dlt, "p": u}) \
        .sort_values("p")
    de.head(60).to_csv(f"{RES}/204_de_dependent_vs_rest_within_B.csv",
                       index=False)
    log("  B 系內部 %d 株依賴者 vs %d 株其餘，差異表現前 12："
        % (int(dd.sum()), int((~dd).sum())))
    log("    " + ", ".join("%s(%+.2f,p=%.1e)" % (a, b, c) for a, b, c in
                           zip(de.gene.head(12), de.delta_z.head(12),
                               de.p.head(12))))
    t8["de_dependent_vs_rest_within_B"] = de.head(25).to_dict("records")
    out["T8_exploratory"] = t8

    # ---------- 排名與總結 ---------------------------------------------
    prim = CAND[(CAND.outcome == "SLC7A11") &
                (CAND.predicted_sign != "none")].copy()
    tab = []
    for feat, sg in prim.groupby("feature"):
        row = {"feature": feat,
               "attenuation_shrinkage": att.get(feat, {}).get("shrinkage"),
               "attenuation_emp_p": att.get(feat, {}).get("empirical_p"),
               "wrong_sign": att.get(feat, {}).get("wrong_sign_vs_prediction")}
        for st in STRATA:
            s = sg[sg.stratum == st]
            if len(s):
                s = s.iloc[0]
                row[f"r_{st}"] = round(float(s.spearman_r), 3)
                row[f"emp_p_{st}"] = (
                    round(float(s.empirical_p_nullA_random_genes), 4)
                    if s.null_used_for_verdict == "A" and
                    s.empirical_p_nullA_random_genes ==
                    s.empirical_p_nullA_random_genes
                    else round(float(s.empirical_p_nullB_ratematched_genes), 4)
                    if s.empirical_p_nullB_ratematched_genes ==
                    s.empirical_p_nullB_ratematched_genes else None)
                row[f"verdict_{st}"] = s.verdict
        tab.append(row)
    TAB = pd.DataFrame(tab).sort_values("r_B_nonPlasmaCell")
    TAB.to_csv(f"{RES}/204_primary_summary.csv", index=False)
    log("\n=== 主要結果變數（SLC7A11）各分層判決 ===")
    log("\n" + TAB.to_string(index=False))
    out["primary_summary"] = TAB.to_dict("records")
    # B_lineage 分層通過不算：那個分層裡 19 株漿細胞骨髓瘤全都不依賴，
    # 任何「非漿細胞 B 細胞比較高/低」的特徵都會過。真正的譜系控制是
    # non_lymphoid（或非造血系）與 B_nonPlasmaCell。
    # 另外，方向與預設預測相反的一律不算通過。
    ctrl_strata = ["non_lymphoid", "non_haematopoietic", "B_nonPlasmaCell"]
    # 通過與否只看各分層自己的判決（verdict 本身已經要求方向正確）。
    # 衰減檢定的「方向相反」旗標不拿來否決分層判決：像 exp_GPX4 全庫 r 幾乎是 0
    # （−0.006），那個旗標只是反映全庫沒有關係，不代表通過的那一層方向錯。
    survived = sorted(set(prim.feature[prim.stratum.isin(ctrl_strata) &
                                       (prim.verdict == "pass")]))
    # 但是跨分層變號要照實標出來：只在一個小分層通過、在大分層方向相反的，
    # 穩定性很差，只能當假說。
    flips = {}
    for feat, sg in prim.groupby("feature"):
        rs = {r0.stratum: r0.spearman_r for r0 in sg.itertuples()}
        ps = [s for s in ctrl_strata
              if s in rs and sg[sg.stratum == s].verdict.iloc[0] == "pass"]
        if not ps:
            continue
        big = [s for s in ("non_lymphoid", "non_haematopoietic", "all")
               if s in rs]
        if big and any(np.sign(rs[ps[0]]) != np.sign(rs[s]) and abs(rs[s]) > 0.03
                       for s in big):
            flips[feat] = {"passing_stratum": ps[0],
                           "r_by_stratum": {k: round(float(v), 3)
                                            for k, v in rs.items()}}
    # 最嚴格的那一關：B 系去掉漿細胞之後還成立的（12 株依賴 vs 41 株不依賴）
    strict = sorted(set(prim.feature[(prim.stratum == "B_nonPlasmaCell") &
                                     (prim.verdict == "pass")]))
    failed = sorted(set(prim.feature) - set(survived))

    # 陰性對照之二的把關：一個特徵如果在它通過的那個分層裡，把「B 身分基因
    # 自己的依賴性」也預測得一樣好，那它測到的是 B 細胞狀態，不是胱胺酸代謝。
    bid_fail = {}
    allc = CAND[CAND.outcome.isin([TARGET, "Bidentity_negctrl"])]
    for feat in survived:
        for st in ctrl_strata:
            v = prim[(prim.feature == feat) & (prim.stratum == st)]
            if not len(v) or v.verdict.iloc[0] != "pass":
                continue
            b = allc[(allc.feature == feat) & (allc.stratum == st) &
                     (allc.outcome == "Bidentity_negctrl")]
            if len(b) and b.verdict.iloc[0] == "pass":
                bid_fail[f"{feat}@{st}"] = {
                    "r_SLC7A11": round(float(v.spearman_r.iloc[0]), 3),
                    "r_Bidentity_control": round(float(b.spearman_r.iloc[0]), 3),
                    "emp_p_Bidentity_control": round(
                        float(b.empirical_p_nullA_random_genes.iloc[0]), 4),
                    "note": "這個分層裡連 B 身分基因的依賴性也一起預測到了，"
                            "所以在這個分層它是 B 細胞狀態的指標，"
                            "不是 SLC7A11 專一的機制"}
    for k, v in bid_fail.items():
        log("  ！陰性對照沒守住：%s 對 SLC7A11 r = %+.3f，"
            "但對 B 身分依賴也有 r = %+.3f（經驗 p = %.3f）"
            % (k, v["r_SLC7A11"], v["r_Bidentity_control"],
               v["emp_p_Bidentity_control"]))
    out["summary"] = {
        "survived_lineage_control_any_stratum": survived,
        "survived_strictest_stratum_B_nonPlasmaCell": strict,
        "sign_flips_across_strata_treat_as_hypothesis": flips,
        "attenuation_wrong_sign_all_lines": sorted(
            k for k, v in att.items()
            if v.get("attenuation_wrong_sign_all_lines")),
        "did_not_survive": failed,
        "failed_Bidentity_negative_control": bid_fail,
        "n_candidates_tested": int(prim.feature.nunique()),
        "decision_rule": "pass 必須：方向正確 + 超過虛無分佈的 5/95 分位 "
                         "+ 參數 p < 0.05。算『撐過譜系控制』的分層只有 "
                         "non_lymphoid、non_haematopoietic、B_nonPlasmaCell；"
                         "B_lineage 不算，因為那層裡 19 株漿細胞骨髓瘤全不依賴，"
                         "任何區分漿細胞與非漿細胞的特徵都會過。",
        "null_used": "表現類特徵用 Null A（隨機基因/隨機面板的表現），"
                     "株層級特徵用 Null B（依賴率配對隨機基因的 gene effect）；"
                     "兩種的經驗 p 都記在 204_candidate_ranking.csv。",
    }
    log("最嚴格分層（B 系去掉漿細胞，12 依賴 vs 41 不依賴）還通過的：%s"
        % (", ".join(strict) or "（無）"))
    for k, v in flips.items():
        log("  跨分層變號（只能當假說）：%s 在 %s 通過，但 %s"
            % (k, v["passing_stratum"],
               " / ".join("%s r=%+.3f" % (a, b)
                          for a, b in v["r_by_stratum"].items())))
    log("\n在譜系控制後仍然通過的：%s" % (", ".join(survived) or "（無）"))
    log("沒有通過的：%s" % ", ".join(failed))

    # ---------- 追加：T1 那幾個組間對比各自的經驗虛無 -------------------
    # 放在最後面是刻意的：前面所有抽樣都已經做完，追加在這裡不會動到
    # 任何已經報出去的數字。
    # 另外要講清楚一件事：淋巴系內部那個對比的比較組（T/NK + Hodgkin + 漿細胞）
    # 是**看到 19 株骨髓瘤 0 株依賴之後才這樣分的**，所以下面 B vs T/NK、
    # B vs 髓系這兩個事前就該做的對比，必須各自補上虛無。
    log("\n=== 追加：B vs T/NK、B vs 髓系的經驗虛無 ===")
    extra = {}
    for nm, rows_m, group_m, desc in (
            ("allB_vs_other_lymphoid", is_lymph, is_bcell,
             "只看淋巴系：依賴者中屬於 B 系（含漿細胞）的比例"),
            ("Bnp_vs_TNK", is_bnp | is_t, is_bnp,
             "只看 B（非漿細胞）與 T/NK：依賴者中屬於 B 的比例"),
            ("Bnp_vs_myeloid", is_bnp | lin.eq("Myeloid"), is_bnp,
             "只看 B（非漿細胞）與髓系：依賴者中屬於 B 的比例")):
        nf, ob = matched_frac_null(TARGET, rows=rows_m.to_numpy(),
                                   group=group_m.to_numpy())
        ep = float((nf >= ob).mean()) if len(nf) else np.nan
        q95 = float(np.quantile(nf, 0.95)) if len(nf) else np.nan
        n_t = int((dep & rows_m).sum())
        extra[nm] = {
            "description": desc, "n_lines_in_subset": int(rows_m.sum()),
            "n_dependent_in_subset": n_t, "observed_frac": ob,
            "n_matched_genes": int(len(nf)),
            "null_mean": float(np.mean(nf)) if len(nf) else np.nan,
            "null_q95": q95, "empirical_p": ep,
            # 虛無的 95 分位就是 1.0 時，沒有任何觀察值可能顯著，
            # 這種情況只能說檢定力不足，不能說支持或否證。
            "null_saturated_underpowered": bool(q95 >= 1.0)}
        log("  %-24s 觀察 %.3f（%d 株依賴）vs 虛無 %.3f（95 分位 %.3f）"
            "→ 經驗 p = %s%s"
            % (nm, ob, n_t, np.mean(nf) if len(nf) else np.nan, q95,
               "%.4f" % ep if ep == ep else "-",
               "   ← 虛無 95 分位已達 1.0，檢定力不足" if q95 >= 1.0 else ""))
    # 淋巴系內部那個對比也標上同樣的警告
    wl = out["T1_scope"].get("empirical_null_within_Bcell", {})
    if wl.get("null_q95", 0) >= 1.0:
        wl["null_saturated_underpowered"] = True
        wl["interpretation"] = ("虛無的 95 分位就是 1.0，所以漿細胞那個分界"
                                "再怎麼乾淨都不可能達到顯著；"
                                "經驗 p = %.3f 只能算提示，不是已確立的事實。"
                                % wl.get("empirical_p", float("nan")))
    out["T1_extra_group_contrasts"] = extra

    # 總結要照自己訂的規則來：撐過分層(A) 與 衰減(B) 兩關都過的才算
    both = sorted(f for f in out["summary"]["survived_lineage_control_any_stratum"]
                  if (att.get(f, {}).get("empirical_p") or 1) < 0.05
                  and not att.get(f, {}).get("attenuation_wrong_sign_all_lines"))
    out["summary"]["survived_both_A_strata_and_B_attenuation"] = both
    out["summary"]["survived_A_only_failed_B_attenuation"] = sorted(
        set(out["summary"]["survived_lineage_control_any_stratum"]) - set(both))
    n_ctrl_tests = int(len(prim[prim.stratum.isin(ctrl_strata)]))
    n_ctrl_pass = int(len(prim[prim.stratum.isin(ctrl_strata) &
                               (prim.verdict == "pass")]))
    out["summary"]["multiplicity"] = {
        "n_tests_in_control_strata": n_ctrl_tests,
        "n_passes": n_ctrl_pass,
        "expected_passes_by_chance_at_0.05": round(0.05 * n_ctrl_tests, 2)}
    log("\n兩關（分層 + 衰減）都過的：%s" % (", ".join(both) or "（無）"))
    log("只過分層、衰減沒過的：%s"
        % (", ".join(out["summary"]["survived_A_only_failed_B_attenuation"])
           or "（無）"))
    log("多重檢定：譜系控制分層裡共 %d 個檢定、%d 個通過，"
        "0.05 下期望約 %.1f 個偶然通過"
        % (n_ctrl_tests, n_ctrl_pass, 0.05 * n_ctrl_tests))

    json.dump(out, open(f"{RES}/204_lineage_mechanism.json", "w"),
              indent=2, ensure_ascii=False, default=str)
    log("\nDONE → results/204_*")


if __name__ == "__main__":
    main()
