#!/usr/bin/env python
"""步驟 205：在 B 系淋巴腫瘤裡找「下一個槓桿」——與 xCT 軸共依賴的可成藥基因。

為什麼要找新槓桿
----------------
這條軸的三個明顯施力點都已經被自己的資料否掉了，不要再回去試：
  (1) 餓胱胺酸：步驟 202 的化學計量模型說降低胱胺酸會讓 xCT **變得比較不必需**
      （轉硫途徑的相對份額上升）；而 B 細胞本來就沒有轉硫途徑（CBS 0/46、
      CBSL 0/46、CTH 1/46 個正常 B 狀態 > 0.1），所以敲掉 SLC7A11 在**每一個**
      濃度下都抽掉 100% 的穀胱甘肽產能——濃度不是可以調的旋鈕。
  (2) 雙硫死亡（disulfidptosis）：步驟 170 完全沒訊號，每個基因集對虛無的平均
      z 落在 −0.97 到 +0.50 之間，連陽性對照組都只有 +0.50，所以那個測試同時
      也沒有檢定力。
  (3) 打 MAFG 蛋白：MAFG 在任何地方都不是遺傳依賴（DepMap 中位 gene effect
      −0.157、只有 2.5% 的株低於 −0.5、MCL 0/3、選擇性 z = +0.36 方向還相反）。

所以這一步改問一個不同的問題：**在沒有旁路的 B 系細胞裡，還有哪個基因的喪失
會跟這條軸一起致死？** 用 DepMap 的 CRISPR gene effect 做共依賴掃描。

名詞要誠實
----------
DepMap 是單基因敲除，算出來的是**共依賴**（co-dependency，兩個基因在同一批株
裡一起變得必需）與**表現→依賴耦合**，不是實驗證明的合成致死。全文與輸出一律
用「共依賴」，不寫「合成致死」。

細胞株分群（在看任何結果之前就定好）
------------------------------------
DepMap 24Q4 有 CRISPR 的 1,178 株裡，Lymphoid 譜系 93 株，但那 93 株並不同質：
18 株 T/NK、3 株霍奇金、19 株漿細胞骨髓瘤（PCM）。而步驟 135 那 12 株依賴
SLC7A11 的淋巴株全是成熟或前驅 B 細胞腫瘤，沒有一株 PCM。

  主要組 B53  ：成熟 ＋ 前驅 B 細胞腫瘤，**不含 PCM**（n ≈ 53）。
                理由是先驗的、不是看結果挑的：本專案整條軸都建立在
                「與正常扁桃體成熟 B 細胞做分化配對」上，漿細胞是終末分化、
                另一個分化區間。
  敏感性組 B72：B53 ＋ PCM（n ≈ 72）。命中必須在這一組維持同號。
  外組 1 HAEM ：T/NK ＋ 霍奇金 ＋ Myeloid 譜系。用來把「B 細胞」跟
                「血液系／懸浮培養」這兩件事分開。
  外組 2 SOLID：其餘所有非淋巴非髓系的實體瘤株。

B53 內部的粗層（控制譜系時用，只有 3 層，才不會把自由度吃光）：
  mature_GC    DLBCLNOS, GCB, ABC, BL, PEL, MBN
  mature_nonGC CLLSLL, MCL, SMZL, BPLL
  precursor_B  BLL, BLLKMT2A, BLLTCF3PBX1, BLLETV6RUNX1, BLLHYPO

錨點
----
主錨 SLC7A11 單基因。次錨 = SLC7A11, GCLC, GCLM, AIFM2 四個基因在 B53 內
z 標準化後的平均（**故意不含 SLC3A2**：它有 45% 的株是依賴，而且同時是 LAT1
的重鏈，放進去會把 SLC7A5／mTORC1 整串拖進來）。

候選基因宇宙（先過濾，才不會被泛必需基因與沒人依賴的基因淹掉）
----------------------------------------------------------------
  U1 在 B53 內資料完整（無缺值）。
  U2 在 B53 內 gene effect 的標準差 ≥ 0.10（平的基因沒有資訊）。
  U3 在 B53 內至少 10% 的株是依賴（gene effect < −0.5）。
     「相關但沒人依賴」就是 MAFG 那個失敗的翻版。
  U4 非泛必需：全 1,178 株裡依賴比例 < 0.90（泛必需基因會霸佔任何共依賴榜，
     而且本來就是爛標的）。泛必需旗標一律照實報。

可否證的預測
------------
Gate0 模組自身的共依賴（步驟 135 全庫 mean pairwise r = 0.278）必須在 B53
      **內部**還活著。如果模組在 B 系內部散掉，那整個掃描就沒有可追的訊號，
      要照實說。
P1  主統計量 = B53 內的 Spearman（rank 相關）。至少要有一個非模組基因的
    r_B 超過「層內置換」的最大統計量虛無的 95 分位（FWER 0.05）。
P2  選擇效應要校正。拿一萬多個基因掃出來的 top k 去跟「隨機夥伴基因」比是
    循環論證，因為掃描本身就是那個分佈。所以用**換錨虛無**：把 SLC7A11 換成
    **整個**在全庫依賴比例上與它配對（rate-matched、且在 B53 內資料完整且
    不平）的錨點池，整個掃描重跑，比較觀察到的第 k 名與虛無的第 k 名。
    三個會讓這個檢定假性通過的陷阱，事先堵掉：
      (a) SLC7A11 自己也在候選宇宙裡，它跟自己的 r = 1.0 會霸佔第一名 →
          觀察端與虛無端都**把模組五個基因的欄位剔掉**。
      (b) 換錨基因若本身也在宇宙裡，它跟自己的 r = 1.0 會污染虛無 →
          把那一格設成 −inf。
      (c) 經驗 p 用 (b+1)/(m+1)，不讓它變成 0；FDR 才不會假性為 0。
    每個虛無錨點都保留它自己真正的夥伴基因，所以這個檢定是**保守**的。
P3  譜系控制，三個版本都要報：
    P3a 粗層：在 B53 的 3 個粗層內各自去均值之後，r 要保留 ≥ 70%。
    P3b B 細胞身分分數：步驟 135 顯示 PAX5/EBF1/CD79A/CD79B/SYK 跟模組的
        全庫 r = 0.242，跟模組自身的 0.278 幾乎一樣高——也就是全庫共依賴
        大半只是「是 B 細胞」。所以把這五個基因 gene effect 的 z 平均當成
        一個「B 細胞身分依賴分數」偏迴歸掉，r 還要保留 ≥ 70%。
    P3c 全控制：粗層虛擬變數 ＋ 篩選深度 ＋ B 細胞身分分數一起迴歸掉，
        r 要同號且保留 ≥ 50%。
P4  非 B 系要弱。r 在 SOLID 與 HAEM 兩個外組都要明顯小於 r_B
    （預設門檻 r_out ≤ 0.15）；在 HAEM 也強的話就是「血液系」不是「B 系」。
P5  篩選品質共變項。每株的全基因平均 gene effect（篩選深度／生長速率的代理）
    偏迴歸掉之後，r 要保留 ≥ 70%。
P6  Jackknife。把 B53 裡最依賴 SLC7A11 的 3 株丟掉重算，r 要同號且保留 ≥ 60%。
P7  B72 敏感性組同號。
P8  表現→依賴耦合（第二條獨立的證據線，不是必要條件）：病人的腫瘤是
    SLC7A11 **高**（比配對正常 B 高 200 倍），而 DepMap 那 3 株 MCL 全部低於
    中位（第 44／17／12 百分位），所以只能用 B 系整組。預測：命中基因的依賴性
    會被 B53 內的 SLC7A11 表現預測（r < 0，表現越高越依賴），而且效應量要
    超過 1,000 個隨機基因表現當預測子的虛無。
P9  機制導向小組（事先寫死，所以多重比較負擔很小、檢定力比全基因掃描高）。
    既然 B 細胞沒有轉硫途徑、GSH 完全押在 xCT 上，機制上可能的夥伴只有
    幾類：GSH 以外的硫醇還原系統、硒半胱胺酸／GPX4 成熟、BH4 與 CoQ 這些
    不經 GSH 的脂質 ROS 防線、NADPH 供應、鐵的處理、半胱胺酸的其他來源與
    消耗、以及這條軸的轉錄調控。小組名單見 PANEL，對整個換錨池算每個基因的
    經驗 p（同樣 (b+1)/(m+1)、同樣剔掉自我相關），再在小組內做 BH FDR。
    預測：至少一個**非模組**基因 FDR < 0.10。SLC3A2 本身在小組名單裡，但它是
    模組成員，所以它過了只能算陽性對照、不算發現。

對照（在看結果之前就宣告，判準全部寫成數字）
--------------------------------------------
陽性對照（它們是對照，**不會**進入推薦名單）：
  GCLC、GCLM、AIFM2 **三個分開報**，每個都必須落在候選宇宙的前 5%。
  如果三個裡有任何一個落在 5% 之外，或因為 U3 被濾掉（AIFM2 全庫只有 0.6%
  的株是依賴，很可能根本進不了宇宙），就要照實寫「全基因掃描在這個樣本數下
  檢定力不足」，不能只報過的那幾個。
陰性對照：
  (a) B 細胞身分基因 PAX5、EBF1、CD79A、CD79B、SYK。判準（數字，先寫死）：
      它們在 P3a 粗層去均值之後的 r 必須**低於 P1 的 FWER 門檻**，而且不得
      被判為存活。原始 r 接近 0 時「保留比例」沒有意義，所以不用比例當判準。
      如果有任何一個沒塌，就是 3 層粗層不夠力，結論必須寫出來，並以 P3b
      （偏迴歸掉 B 細胞身分分數）為主要的譜系控制。
  (b) RPL13A（核糖體蛋白、泛必需），跟步驟 203 用同一個陰性對照。
  (c) MAFG、NFE2L2：已經被否掉的兩個，不應該冒出來。

泛必需不是二分的
----------------
U4 的 0.90 只是硬下限。U3（B53 內至少 10% 的株依賴）本身就會把廣泛必需的基因
富集進來，所以每個候選都必須把 depfrac_all 與 depfrac_B53 兩個數字報出來；
在 60–89% 的株都是依賴的基因不能當成乾淨的標的，要明寫。

排名規則（先寫死）
------------------
先按 r_B（主錨）排序，然後必須通過 P3a、P3b、P3c、P5、P6、P7 與 U4 才算
「存活」；
P1／P2 是全名單層級的檢定，照實報。最後名單最多 5 個，再加上
PubMed 的新穎性與抑制劑查詢（排名凍結之後才查）。

任何一項不過都要寫進結論，不是只報過的那幾項。「什麼都沒活下來」是可以接受
而且必須照實寫的答案。

輸出 results/205_*
"""
import os
import json
import time
import urllib.parse
import urllib.request
import ssl
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
EUT = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
TOOLNAME = "codep-screen"          # 一般性工具名；不放任何個人資訊

ANCHOR = "SLC7A11"
MODULE = ["SLC7A11", "SLC3A2", "GCLC", "GCLM", "AIFM2"]
MOD_ANCHOR = ["SLC7A11", "GCLC", "GCLM", "AIFM2"]   # 次錨：刻意不含 SLC3A2
POS_CTRL = ["GCLC", "GCLM", "AIFM2", "SLC3A2"]
NEG_CTRL_BID = ["PAX5", "EBF1", "CD79A", "CD79B", "SYK"]
NEG_CTRL_ESS = ["RPL13A"]
NEG_CTRL_DEAD = ["MAFG", "NFE2L2"]

# P9 機制導向小組（事先寫死，不得在看結果後增刪）
PANEL = {
    "thiol_non_GSH": ["TXNRD1", "TXNRD2", "TXN", "TXN2", "PRDX1", "PRDX2",
                      "PRDX6", "SRXN1", "GLRX", "GSR"],
    "selenocysteine_GPX4": ["GPX4", "SEPHS2", "PSTK", "SEPSECS", "EEFSEC",
                            "SCLY", "SELENOI", "SEPHS1"],
    "BH4_CoQ_nonGSH_lipid": ["GCH1", "PTS", "SPR", "QDPR", "COQ2", "COQ3",
                             "FDFT1", "HMGCR", "POR", "NQO1"],
    "NADPH": ["G6PD", "PGD", "IDH1", "ME1", "NNT", "NADK", "TALDO1"],
    "iron": ["FTH1", "FTL", "NCOA4", "SLC40A1", "STEAP3", "TFRC", "IREB2",
             "ACO1"],
    "cysteine_supply_use": ["CHAC1", "GGCT", "GSS", "OPLAH", "CTH", "CBS",
                            "MPST", "NFS1", "SLC1A5", "SLC3A2", "SLC7A5"],
    "lipid_PUFA_substrate": ["ACSL4", "LPCAT3", "ALOX5", "ALOX15", "ELOVL5",
                             "SCD", "FADS2"],
    "axis_transcription": ["ATF4", "NFE2L2", "KEAP1", "MAFG", "MAFK", "MAFF",
                           "BACH1", "CEBPB", "EIF2AK4", "DDIT3"],
}

B_CODES = {"PCM", "DLBCLNOS", "BL", "BLL", "CLLSLL", "BLLKMT2A", "MCL", "PEL",
           "BLLTCF3PBX1", "BLLETV6RUNX1", "GCB", "MBN", "SMZL", "BLLHYPO",
           "BPLL", "ABC", "LPL", "MZL", "FL", "HGBCL", "PBL", "MM"}
PLASMA_CODES = {"PCM", "MM"}
TNK_CODES = {"ALCLALKP", "ALCLALKN", "TLL", "ATLL", "MTNN", "PTCL", "TPLL",
             "TMN", "NKCL"}
HL_CODES = {"HL", "CHL", "NLPHL"}
STRATA = {
    "mature_GC": {"DLBCLNOS", "GCB", "ABC", "BL", "PEL", "MBN", "FL", "HGBCL",
                  "PBL"},
    "mature_nonGC": {"CLLSLL", "MCL", "SMZL", "BPLL", "LPL", "MZL"},
    "precursor_B": {"BLL", "BLLKMT2A", "BLLTCF3PBX1", "BLLETV6RUNX1",
                    "BLLHYPO"},
    "plasma": {"PCM", "MM"},
}

DEP_CUT = -0.5          # 步驟 135 的慣例
MIN_SD_B = 0.10         # U2
MIN_DEPFRAC_B = 0.10    # U3
PANESS_FRAC = 0.90      # U4
R_OUT_MAX = 0.15        # P4
KEEP_FRAC = 0.70        # P3／P5
KEEP_FRAC_JK = 0.60     # P6
N_PERM = 1000           # P1 層內置換
N_ANCHOR_MAX = 2000     # P2 換錨池上限（原則上用整個配對池）
KEEP_FRAC_FULL = 0.50   # P3c
N_EXP_NULL = 1000       # P8 隨機表現預測子
TOP_K_REPORT = 40
RNG = np.random.default_rng(205)


def log(*a):
    print("[205]", *a, flush=True)


# ---------------------------------------------------------------- 小工具
def zrank(M):
    """對每一欄做 rank 再 z 標準化；rank 上的 Pearson 就是 Spearman。"""
    R = stats.rankdata(M, axis=0).astype(np.float64)
    R -= R.mean(axis=0, keepdims=True)
    sd = R.std(axis=0, keepdims=True)
    sd[sd == 0] = np.inf
    return R / sd


def zrank_vec(v):
    r = stats.rankdata(v).astype(np.float64)
    r -= r.mean()
    s = r.std()
    return r / (s if s > 0 else np.inf)


def center_within(M, labels):
    """在每個層內把欄位去均值（控制粗層），之後重新單位化。"""
    out = np.array(M, dtype=np.float64, copy=True)
    for lab in np.unique(labels):
        m = labels == lab
        out[m] = out[m] - out[m].mean(axis=0, keepdims=True)
    nrm = np.linalg.norm(out, axis=0, keepdims=True)
    nrm[nrm == 0] = np.inf
    return out / nrm


def resid_on(M, cov):
    """把共變項 cov 從每一欄迴歸掉，再單位化。"""
    M = np.array(M, dtype=np.float64, copy=True)
    c = cov / np.linalg.norm(cov)
    out = M - np.outer(c, c @ M)
    nrm = np.linalg.norm(out, axis=0, keepdims=True)
    nrm[nrm == 0] = np.inf
    return out / nrm


def resid_on_design(M, X):
    """把設計矩陣 X（會自動加截距）從 M 的每一欄迴歸掉，再單位化。"""
    M = np.array(M, dtype=np.float64, copy=True)
    X = np.column_stack([np.ones(len(M)), np.asarray(X, dtype=np.float64)])
    Q, _ = np.linalg.qr(X)
    out = M - Q @ (Q.T @ M)
    nrm = np.linalg.norm(out, axis=0, keepdims=True)
    nrm[nrm == 0] = np.inf
    return out / nrm


def bh_fdr(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    o = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        prev = min(prev, p[o[i]] * n / (i + 1))
        q[o[i]] = prev
    return q


def _ssl_ctx():
    """這台機器的 HTTPS 走會換憑證的代理，Python 自帶的 CA bundle 驗不過
    （SSL: CERTIFICATE_VERIFY_FAILED），所以改用系統的 CA bundle。"""
    for ca in ("/etc/ssl/certs/ca-certificates.crt",
               "/etc/pki/tls/certs/ca-bundle.crt"):
        if os.path.exists(ca):
            try:
                return ssl.create_default_context(cafile=ca)
            except Exception:
                pass
    return None


SSLCTX = _ssl_ctx()


def pubmed_count(term, retry=3):
    q = urllib.parse.urlencode({"db": "pubmed", "term": term,
                                "retmode": "json", "rettype": "count",
                                "tool": TOOLNAME})
    last = None
    for _ in range(retry):
        try:
            with urllib.request.urlopen(f"{EUT}?{q}", timeout=30,
                                        context=SSLCTX) as r:
                d = json.load(r)
            return int(d["esearchresult"]["count"])
        except Exception as e:
            last = "%s: %s" % (type(e).__name__, e)
            time.sleep(2.0)
    log("   PubMed 查詢失敗（%s）：%s" % (last, term))
    return None


# ---------------------------------------------------------------- 主程式
def main():
    out = {"design": {
        "anchor": ANCHOR, "module_anchor": MOD_ANCHOR,
        "dep_cut": DEP_CUT, "min_sd_B": MIN_SD_B,
        "min_depfrac_B": MIN_DEPFRAC_B, "paness_frac": PANESS_FRAC,
        "r_out_max": R_OUT_MAX, "keep_frac": KEEP_FRAC,
        "keep_frac_jackknife": KEEP_FRAC_JK, "n_perm": N_PERM,
        "n_anchor_max": N_ANCHOR_MAX, "n_exp_null": N_EXP_NULL,
        "label": "co-dependency and expression-dependency coupling; "
                 "NOT demonstrated synthetic lethality"}}

    # ---------- 讀資料 ----------
    log("讀 Model.csv …")
    mod = pd.read_csv(f"{DM}/Model.csv", low_memory=False).set_index("ModelID")
    log("讀 CRISPRGeneEffect.csv（全矩陣）…")
    t0 = time.time()
    EFF = pd.read_csv(f"{DM}/CRISPRGeneEffect.csv", index_col=0,
                      low_memory=False)
    EFF.columns = [c.split(" (")[0] for c in EFF.columns]
    EFF = EFF.loc[:, ~EFF.columns.duplicated()]
    log("  %d 株 × %d 基因（%.0f 秒）" % (EFF.shape[0], EFF.shape[1],
                                          time.time() - t0))

    code = mod["OncotreeCode"].astype(str).reindex(EFF.index)
    lin = mod["OncotreeLineage"].astype(str).reindex(EFF.index)
    sub = mod["OncotreeSubtype"].astype(str).reindex(EFF.index)

    is_b = code.isin(B_CODES) & (lin == "Lymphoid")
    is_plasma = code.isin(PLASMA_CODES)
    is_tnk = code.isin(TNK_CODES) & (lin == "Lymphoid")
    is_hl = code.isin(HL_CODES) & (lin == "Lymphoid")
    is_mye = lin == "Myeloid"

    B53 = EFF.index[is_b & ~is_plasma]
    B72 = EFF.index[is_b]
    HAEM = EFF.index[is_tnk | is_hl | is_mye]
    SOLID = EFF.index[~(lin == "Lymphoid") & ~is_mye]

    strat = pd.Series("other", index=EFF.index)
    for k, v in STRATA.items():
        strat[code.isin(v)] = k

    log("B53 n=%d｜B72 n=%d｜HAEM n=%d（T/NK %d＋HL %d＋Myeloid %d）｜SOLID n=%d"
        % (len(B53), len(B72), len(HAEM), int(is_tnk.sum()), int(is_hl.sum()),
           int(is_mye.sum()), len(SOLID)))
    log("B53 粗層：%s" % strat[B53].value_counts().to_dict())
    out["groups"] = {
        "n_B53": len(B53), "n_B72": len(B72), "n_HAEM": len(HAEM),
        "n_SOLID": len(SOLID),
        "B53_strata": {k: int(v) for k, v in
                       strat[B53].value_counts().items()},
        "B53_subtypes": {k: int(v) for k, v in
                         code[B53].value_counts().items()},
        "B53_lines": sorted(B53.tolist())}

    # ---------- 候選宇宙 ----------
    EB = EFF.loc[B53]
    depfrac_all = (EFF < DEP_CUT).sum(axis=0) / EFF.notna().sum(axis=0)
    depfrac_B = (EB < DEP_CUT).sum(axis=0) / EB.notna().sum(axis=0)
    sd_B = EB.std(axis=0)
    complete_B = EB.notna().all(axis=0)
    paness = depfrac_all >= PANESS_FRAC

    univ = complete_B & (sd_B >= MIN_SD_B) & (depfrac_B >= MIN_DEPFRAC_B) \
        & ~paness
    genes = EB.columns[univ].tolist()
    gidx = {g: i for i, g in enumerate(genes)}
    log("宇宙過濾：完整 %d → 加 sd≥%.2f %d → 加 B 依賴比例≥%.2f %d → 去泛必需 %d"
        % (int(complete_B.sum()),
           MIN_SD_B, int((complete_B & (sd_B >= MIN_SD_B)).sum()),
           MIN_DEPFRAC_B,
           int((complete_B & (sd_B >= MIN_SD_B)
                & (depfrac_B >= MIN_DEPFRAC_B)).sum()),
           len(genes)))
    log("  泛必需（全庫依賴比例 ≥ %.2f）共 %d 個基因" % (PANESS_FRAC,
                                                        int(paness.sum())))
    out["universe"] = {
        "n_complete_in_B53": int(complete_B.sum()),
        "n_universe": len(genes),
        "n_panessential": int(paness.sum()),
        "n_total_genes": int(EFF.shape[1])}

    n_b = len(B53)
    EBU = EB[genes].to_numpy(dtype=np.float64)
    Zuniv = zrank(EBU) / np.sqrt(n_b)

    # ---------- Gate0：模組在 B53 內還有沒有共依賴 ----------
    modg = [g for g in MODULE if g in EB.columns]
    Zmod = zrank(EB[modg].to_numpy(dtype=np.float64)) / np.sqrt(n_b)
    Cm = Zmod.T @ Zmod
    iu = np.triu_indices(len(modg), 1)
    mod_r_B = float(Cm[iu].mean())
    null_mod = np.empty(2000)
    for i in range(2000):
        pick = RNG.choice(len(genes), size=len(modg), replace=False)
        S = Zuniv[:, pick]
        null_mod[i] = (S.T @ S)[iu].mean()
    g0p = float((null_mod >= mod_r_B).mean())
    pw = {}
    for i, a in enumerate(modg):
        for b in modg[i + 1:]:
            pw["%s~%s" % (a, b)] = round(float(Cm[i, modg.index(b)]), 4)
    out["Gate0_module_coherence_in_B53"] = {
        "genes": modg, "mean_pairwise_r": mod_r_B, "pairwise": pw,
        "null_mean": float(null_mod.mean()),
        "null_q95": float(np.quantile(null_mod, 0.95)),
        "empirical_p": g0p,
        "reference_all_lines_r_from_135": 0.2781638299743254,
        "pass": bool(g0p < 0.05)}
    log("Gate0 模組在 B53 內 mean pairwise r = %+.3f（虛無 %.3f，95 分位 %.3f）"
        "→ 經驗 p = %.4f → %s"
        % (mod_r_B, null_mod.mean(), np.quantile(null_mod, 0.95), g0p,
           "過" if g0p < 0.05 else "不過"))
    log("   兩兩：%s" % pw)

    # ---------- 錨點 ----------
    za_B = zrank_vec(EB[ANCHOR].to_numpy(dtype=np.float64)) / np.sqrt(n_b)
    ma = [g for g in MOD_ANCHOR if g in EB.columns]
    modz = zrank(EB[ma].to_numpy(dtype=np.float64))
    zm_B = zrank_vec(modz.mean(axis=1)) / np.sqrt(n_b)

    linemean_B = EBU.mean(axis=1)
    zlm_B = zrank_vec(linemean_B)
    r_lm_anchor = float(np.dot(zlm_B / np.sqrt(n_b), za_B))
    log("篩選深度共變項：B53 每株宇宙平均 gene effect %.3f ± %.3f；"
        "與主錨的 Spearman r = %+.3f"
        % (linemean_B.mean(), linemean_B.std(), r_lm_anchor))
    out["line_mean_covariate"] = {"mean": float(linemean_B.mean()),
                                  "sd": float(linemean_B.std()),
                                  "spearman_r_with_anchor": r_lm_anchor}

    # ---------- 主掃描 ----------
    rB = Zuniv.T @ za_B
    rB_mod = Zuniv.T @ zm_B
    strat_B = strat[B53].to_numpy()
    Zc = center_within(Zuniv * np.sqrt(n_b), strat_B)
    ac = center_within((za_B * np.sqrt(n_b)).reshape(-1, 1), strat_B).ravel()
    rB_strat = Zc.T @ ac
    Zp = resid_on(Zuniv * np.sqrt(n_b), zlm_B)
    ap = resid_on((za_B * np.sqrt(n_b)).reshape(-1, 1), zlm_B).ravel()
    rB_partial = Zp.T @ ap

    # P3b：B 細胞身分依賴分數（PAX5/EBF1/CD79A/CD79B/SYK 的 z 平均）
    bid_g = [g for g in NEG_CTRL_BID if g in EB.columns
             and not EB[g].isna().any()]
    bid_raw = zrank(EB[bid_g].to_numpy(dtype=np.float64)).mean(axis=1)
    zbid = zrank_vec(bid_raw)
    r_bid_anchor = float(np.dot(zbid / np.sqrt(n_b), za_B))
    log("P3b B 細胞身分分數用 %s；它與主錨的 Spearman r = %+.3f"
        % (",".join(bid_g), r_bid_anchor))
    Zb = resid_on(Zuniv * np.sqrt(n_b), zbid)
    ab = resid_on((za_B * np.sqrt(n_b)).reshape(-1, 1), zbid).ravel()
    rB_partial_bid = Zb.T @ ab
    out["bcell_identity_score"] = {
        "genes": bid_g, "spearman_r_with_anchor": r_bid_anchor,
        "reference_all_lines_module_vs_bidentity_r_from_135":
            0.24178731257089617}

    # P3c：粗層虛擬變數 ＋ 篩選深度 ＋ B 身分分數一起迴歸掉
    D = np.column_stack(
        [(strat_B == s).astype(float) for s in np.unique(strat_B)[1:]]
        + [zlm_B, zbid])
    Zf = resid_on_design(Zuniv * np.sqrt(n_b), D)
    af = resid_on_design((za_B * np.sqrt(n_b)).reshape(-1, 1), D).ravel()
    rB_full = Zf.T @ af

    res = pd.DataFrame({
        "gene": genes, "r_B53": rB, "r_B53_module_anchor": rB_mod,
        "r_B53_strat_resid": rB_strat, "r_B53_partial_linemean": rB_partial,
        "r_B53_partial_bidentity": rB_partial_bid,
        "r_B53_full_control": rB_full,
        "depfrac_B53": depfrac_B[genes].to_numpy(),
        "depfrac_all": depfrac_all[genes].to_numpy(),
        "sd_B53": sd_B[genes].to_numpy(),
        "mean_eff_B53": EBU.mean(axis=0),
        "mean_eff_all": EFF[genes].mean(axis=0).to_numpy()})

    # ---------- 外組與敏感性組 ----------
    def group_r(idx, label):
        ii = [i for i in idx if pd.notna(EFF.at[i, ANCHOR])]
        n = len(ii)
        a = zrank_vec(EFF.loc[ii, ANCHOR].to_numpy(dtype=np.float64)) \
            / np.sqrt(n)
        G = EFF.loc[ii, genes]
        ok = (G.isna().mean(axis=0) <= 0.05).to_numpy()
        Gn = G.loc[:, ok]
        Gn = Gn.fillna(Gn.median(axis=0))
        Z = zrank(Gn.to_numpy(dtype=np.float64)) / np.sqrt(n)
        log("  %s：n=%d（主錨有值），%d/%d 個基因可算（其餘該組缺值 > 5%%）"
            % (label, n, int(ok.sum()), len(genes)))
        return pd.Series(Z.T @ a,
                         index=[g for g, k in zip(genes, ok) if k]), n

    for gname, idx in (("SOLID", SOLID), ("HAEM", HAEM), ("B72", B72)):
        s, nn = group_r(idx, gname)
        res["r_%s" % gname] = res.gene.map(s)
        out["groups"]["n_%s_anchor_ok" % gname] = nn
    res["delta_B53_minus_SOLID"] = res.r_B53 - res.r_SOLID

    # ---------- Jackknife ----------
    drop3 = EB[ANCHOR].nsmallest(3)
    keep = [i for i in B53 if i not in set(drop3.index)]
    aj = zrank_vec(EFF.loc[keep, ANCHOR].to_numpy(dtype=np.float64)) \
        / np.sqrt(len(keep))
    ZJ = zrank(EFF.loc[keep, genes].to_numpy(dtype=np.float64)) \
        / np.sqrt(len(keep))
    res["r_B53_jackknife3"] = ZJ.T @ aj
    log("Jackknife 丟掉：%s" % ", ".join(
        "%s (%s, %.2f)" % (i, sub.get(i, "?"), v) for i, v in drop3.items()))
    out["jackknife_dropped"] = [
        {"line": i, "subtype": str(sub.get(i, "")), "slc7a11_eff": float(v)}
        for i, v in drop3.items()]

    # ---------- P1：層內置換的最大統計量虛無 ----------
    log("P1 層內置換 %d 次，取最大統計量…" % N_PERM)
    a_raw = EB[ANCHOR].to_numpy(dtype=np.float64)
    idx_by_strat = [np.where(strat_B == s)[0] for s in np.unique(strat_B)]
    perm_max = np.empty(N_PERM)
    perm_k = np.empty((N_PERM, TOP_K_REPORT))
    for i in range(N_PERM):
        p = a_raw.copy()
        for ii in idx_by_strat:
            p[ii] = RNG.permutation(p[ii])
        rr = Zuniv.T @ (zrank_vec(p) / np.sqrt(n_b))
        perm_max[i] = rr.max()
        perm_k[i] = np.sort(rr)[::-1][:TOP_K_REPORT]
    thr_fwer = float(np.quantile(perm_max, 0.95))
    nonmod = res[~res.gene.isin(MODULE)]
    n_pass_p1 = int((nonmod.r_B53 > thr_fwer).sum())
    out["P1_within_stratum_permutation"] = {
        "n_perm": N_PERM, "null_max_r_mean": float(perm_max.mean()),
        "null_max_r_q95_FWER05": thr_fwer,
        "observed_max_r_nonmodule": float(nonmod.r_B53.max()),
        "observed_max_gene": str(nonmod.loc[nonmod.r_B53.idxmax(), "gene"]),
        "n_nonmodule_genes_above_FWER": n_pass_p1,
        "pass": bool(n_pass_p1 >= 1)}
    log("P1 虛無最大 r 的 95 分位（FWER 0.05）= %.3f；"
        "觀察到的非模組最大 r = %.3f（%s）→ 超過門檻的基因 %d 個 → %s"
        % (thr_fwer, nonmod.r_B53.max(),
           nonmod.loc[nonmod.r_B53.idxmax(), "gene"], n_pass_p1,
           "過" if n_pass_p1 >= 1 else "不過"))
    # 每個基因的置換 p（未校正，純參考）
    res["perm_p_uncorrected"] = np.nan

    # ---------- P2：換錨虛無（校正選擇效應） ----------
    da = float(depfrac_all[ANCHOR])
    msk = ((depfrac_all - da).abs() <= 0.01) \
        & complete_B.reindex(depfrac_all.index).fillna(False) \
        & (sd_B.reindex(depfrac_all.index) >= MIN_SD_B) \
        & ~pd.Series(depfrac_all.index.isin(MODULE),
                     index=depfrac_all.index)
    cand_anchor = depfrac_all.index[msk.fillna(False)].tolist()
    if len(cand_anchor) > N_ANCHOR_MAX:
        pick = RNG.choice(len(cand_anchor), size=N_ANCHOR_MAX, replace=False)
        swap_genes = [cand_anchor[i] for i in pick]
    else:
        swap_genes = list(cand_anchor)
    m_sw = len(swap_genes)
    log("P2 換錨：SLC7A11 全庫依賴比例 %.4f，配對池 %d 個基因，用 %d 個"
        % (da, len(cand_anchor), m_sw))
    SW = zrank(EB[swap_genes].to_numpy(dtype=np.float64)) / np.sqrt(n_b)
    swap_idx = {g: i for i, g in enumerate(swap_genes)}
    null_r = (Zuniv.T @ SW).T           # (n_anchor, n_gene)
    # (b) 換錨基因本身若在宇宙裡，它跟自己的 r = 1.0 要剔掉
    n_self = 0
    for i, g in enumerate(swap_genes):
        if g in gidx:
            null_r[i, gidx[g]] = -np.inf
            n_self += 1
    # (a) 觀察端與虛無端都把模組五個基因的欄位剔掉
    nonmod_col = np.array([g not in MODULE for g in genes])
    null_sorted = np.sort(null_r[:, nonmod_col], axis=1)[:, ::-1]
    obs_sorted = np.sort(rB[nonmod_col])[::-1]
    p2 = {"n_anchors": m_sw, "n_anchor_self_masked": n_self,
          "n_module_columns_removed": int((~nonmod_col).sum()),
          "note": "觀察端與虛無端都已剔除模組基因欄位與換錨基因的自我相關；"
                  "每個虛無錨點保留自己真正的夥伴基因，所以本檢定偏保守"}
    for k in (1, 5, 10, 25):
        b = int((null_sorted[:, k - 1] >= obs_sorted[k - 1]).sum())
        p2["top%d" % k] = {
            "observed_r": float(obs_sorted[k - 1]),
            "null_mean": float(null_sorted[:, k - 1].mean()),
            "null_q95": float(np.quantile(null_sorted[:, k - 1], 0.95)),
            "empirical_p": (b + 1) / (m_sw + 1)}
        log("   非模組第 %2d 名：觀察 %.3f vs 換錨虛無 %.3f（95 分位 %.3f）"
            "→ 經驗 p = %.4f"
            % (k, obs_sorted[k - 1], null_sorted[:, k - 1].mean(),
               np.quantile(null_sorted[:, k - 1], 0.95),
               p2["top%d" % k]["empirical_p"]))
    p2["pass"] = bool(p2["top5"]["empirical_p"] < 0.05)
    out["P2_anchor_swap_null"] = p2
    res["emp_p_anchor_swap"] = \
        ((null_r >= rB[None, :]).sum(axis=0) + 1) / (m_sw + 1)

    # P2b（事後）：把 B 細胞身分分數從觀察端與虛無端都偏迴歸掉之後再比 top k。
    # 為什麼要補這一項：未校正的非模組 top 20 裡有 CD79A、CD79B、SYK、CD19、
    # FOXO1 五個 B 身分基因，而身分分數與主錨的 r = +0.498；配對錨點沒有這種
    # 耦合，所以原始的 P2 有一部分是「是 B 細胞」撐起來的。
    SWb = resid_on(SW * np.sqrt(n_b), zbid)
    null_bid = (Zb.T @ SWb).T
    for i, g in enumerate(swap_genes):
        if g in gidx:
            null_bid[i, gidx[g]] = -np.inf
    nb_sorted = np.sort(null_bid[:, nonmod_col], axis=1)[:, ::-1]
    ob_sorted = np.sort(rB_partial_bid[nonmod_col])[::-1]
    p2b = {"note": "post hoc；觀察端與虛無端都偏迴歸掉 B 細胞身分分數"}
    for k in (1, 5, 10, 25):
        b = int((nb_sorted[:, k - 1] >= ob_sorted[k - 1]).sum())
        p2b["top%d" % k] = {
            "observed_r": float(ob_sorted[k - 1]),
            "null_mean": float(nb_sorted[:, k - 1].mean()),
            "null_q95": float(np.quantile(nb_sorted[:, k - 1], 0.95)),
            "empirical_p": (b + 1) / (m_sw + 1)}
        log("   P2b（去 B 身分）第 %2d 名：觀察 %.3f vs 虛無 %.3f"
            "（95 分位 %.3f）→ 經驗 p = %.4f"
            % (k, ob_sorted[k - 1], nb_sorted[:, k - 1].mean(),
               np.quantile(nb_sorted[:, k - 1], 0.95),
               p2b["top%d" % k]["empirical_p"]))
    p2b["pass"] = bool(p2b["top5"]["empirical_p"] < 0.05)
    out["P2b_anchor_swap_bidentity_partialled"] = p2b

    # ---------- P8：表現 → 依賴耦合 ----------
    log("讀表現矩陣（全矩陣）…")
    t0 = time.time()
    EXP = pd.read_csv(f"{DM}/OmicsExpressionProteinCodingGenesTPMLogp1.csv",
                      index_col=0, low_memory=False)
    EXP.columns = [c.split(" (")[0] for c in EXP.columns]
    EXP = EXP.loc[:, ~EXP.columns.duplicated()]
    log("  %d 株 × %d 基因（%.0f 秒）" % (EXP.shape[0], EXP.shape[1],
                                          time.time() - t0))
    bx = [i for i in B53 if i in EXP.index and pd.notna(EXP.at[i, ANCHOR])]
    nx = len(bx)
    log("P8 B53 裡有表現資料的 %d 株；SLC7A11 表現中位 %.3f（全庫中位 %.3f）"
        % (nx, float(EXP.loc[bx, ANCHOR].median()),
           float(EXP[ANCHOR].median())))
    ZE = zrank(EFF.loc[bx, genes].to_numpy(dtype=np.float64)) / np.sqrt(nx)
    ze_anchor = zrank_vec(EXP.loc[bx, ANCHOR].to_numpy(dtype=np.float64)) \
        / np.sqrt(nx)
    r_exp = ZE.T @ ze_anchor
    res["r_exp_slc7a11_B53"] = r_exp

    exp_pool = [c for c in EXP.columns if c not in MODULE]
    pickx = RNG.choice(len(exp_pool), size=min(N_EXP_NULL, len(exp_pool)),
                       replace=False)
    XP = EXP.loc[bx, [exp_pool[i] for i in pickx]].to_numpy(dtype=np.float64)
    keepx = (~np.isnan(XP).any(axis=0)) & (np.nanstd(XP, axis=0) > 0)
    ZX = zrank(XP[:, keepx]) / np.sqrt(nx)
    null_exp = (ZE.T @ ZX).T            # (n_pred, n_gene)
    res["emp_p_exp_null"] = (null_exp <= r_exp[None, :]).mean(axis=0)
    thr_exp = float(np.quantile(null_exp.min(axis=1), 0.05))
    nm_exp = res[~res.gene.isin(MODULE)].r_exp_slc7a11_B53
    out["P8_expression_coupling"] = {
        "n_lines_with_expression": nx,
        "n_null_predictors": int(ZX.shape[1]),
        "null_min_r_q05_FWER05": thr_exp,
        "observed_min_r_nonmodule": float(nm_exp.min()),
        "n_nonmodule_below_FWER": int((nm_exp < thr_exp).sum()),
        "pass": bool((nm_exp < thr_exp).sum() >= 1)}
    log("P8 虛無最小 r 的 5 分位（FWER 0.05）= %.3f；觀察到的非模組最小 r = %.3f"
        " → 超過門檻 %d 個"
        % (thr_exp, nm_exp.min(),
           out["P8_expression_coupling"]["n_nonmodule_below_FWER"]))
    med = float(EXP.loc[bx, ANCHOR].median())
    hi = [i for i in bx if EXP.at[i, ANCHOR] > med]
    lo = [i for i in bx if EXP.at[i, ANCHOR] <= med]
    res["eff_hi_slc7a11"] = res.gene.map(EFF.loc[hi, genes].mean())
    res["eff_lo_slc7a11"] = res.gene.map(EFF.loc[lo, genes].mean())
    res["delta_hi_minus_lo"] = res.eff_hi_slc7a11 - res.eff_lo_slc7a11
    out["P8_expression_coupling"].update({"n_hi": len(hi), "n_lo": len(lo),
                                          "median_log_tpm_split": med})

    # ---------- 存活判定 ----------
    res["panessential"] = res.depfrac_all >= PANESS_FRAC
    res["broadly_essential_60_89"] = (res.depfrac_all >= 0.60) \
        & (res.depfrac_all < PANESS_FRAC)
    res["P3a_strat"] = (res.r_B53_strat_resid >= KEEP_FRAC * res.r_B53) \
        & (res.r_B53_strat_resid > 0)
    res["P3b_bidentity"] = \
        (res.r_B53_partial_bidentity >= KEEP_FRAC * res.r_B53) \
        & (res.r_B53_partial_bidentity > 0)
    res["P3c_full_control"] = \
        (res.r_B53_full_control >= KEEP_FRAC_FULL * res.r_B53) \
        & (res.r_B53_full_control > 0)
    res["P4_outgroup"] = (res.r_SOLID <= R_OUT_MAX) \
        & (res.r_HAEM <= R_OUT_MAX)
    res["P5_partial"] = (res.r_B53_partial_linemean >= KEEP_FRAC * res.r_B53) \
        & (res.r_B53_partial_linemean > 0)
    res["P6_jackknife"] = (res.r_B53_jackknife3 >= KEEP_FRAC_JK * res.r_B53) \
        & (res.r_B53_jackknife3 > 0)
    res["P7_B72"] = res.r_B72 > 0
    res["above_FWER"] = res.r_B53 > thr_fwer
    res["survivor"] = res.P3a_strat & res.P3b_bidentity \
        & res.P3c_full_control & res.P5_partial & res.P6_jackknife \
        & res.P7_B72 & ~res.panessential

    # ---------- P9：機制導向小組 ----------
    log("\nP9 機制導向小組（事先寫死）…")
    prows = []
    for grp, gs in PANEL.items():
        for g in gs:
            if g not in EB.columns:
                prows.append({"panel": grp, "gene": g, "status": "不在矩陣"})
                continue
            v = EB[g]
            if v.isna().any():
                ok = v.notna().to_numpy()
                zg = np.zeros(n_b)
                zg[ok] = zrank_vec(v[ok].to_numpy(dtype=np.float64)) \
                    / np.sqrt(ok.sum())
                r = float(stats.spearmanr(v, EB[ANCHOR],
                                          nan_policy="omit").statistic)
            else:
                zg = zrank_vec(v.to_numpy(dtype=np.float64)) / np.sqrt(n_b)
                r = float(np.dot(zg, za_B))
            # 換錨虛無：同一個基因對整個配對錨點池的 r（剔掉自我相關）
            nullv = SW.T @ zg
            if g in swap_idx:
                nullv[swap_idx[g]] = -np.inf
            pe = (int((nullv >= r).sum()) + 1) / (m_sw + 1)
            # 同一件事，但觀察端與虛無端都偏迴歸掉 B 細胞身分分數。
            # P2b 顯示無偏掃描的 top k 去掉 B 身分之後就不顯著了，所以小組
            # 這一層一定要有同樣的版本，否則分不出「硒」與「是 B 細胞」。
            if not v.isna().any():
                rgb = resid_on(zg.reshape(-1, 1) * np.sqrt(n_b),
                               zbid).ravel()
                r_b = float(np.dot(rgb, ab))
                nvb = SWb.T @ rgb
                if g in swap_idx:
                    nvb[swap_idx[g]] = -np.inf
                pe_b = (int((nvb >= r_b).sum()) + 1) / (m_sw + 1)
            else:
                r_b, pe_b = np.nan, np.nan
            why = []
            if not bool(complete_B[g]):
                why.append("B53 有缺值")
            if float(sd_B[g]) < MIN_SD_B:
                why.append("sd %.3f" % sd_B[g])
            if float(depfrac_B[g]) < MIN_DEPFRAC_B:
                why.append("B53 依賴比例 %.3f" % depfrac_B[g])
            if bool(paness[g]):
                why.append("泛必需 %.3f" % depfrac_all[g])
            row = {"panel": grp, "gene": g,
                   "status": "宇宙內" if g in gidx
                   else "宇宙外（%s）" % "；".join(why),
                   "r_B53": r, "emp_p_anchor_swap": pe,
                   "r_B53_bidentity_partialled": r_b,
                   "emp_p_anchor_swap_bidentity": pe_b,
                   "depfrac_B53": float(depfrac_B[g]),
                   "depfrac_all": float(depfrac_all[g]),
                   "panessential": bool(paness[g])}
            if g in gidx:
                rr = res[res.gene == g].iloc[0]
                row.update({
                    "r_B53_strat_resid": float(rr.r_B53_strat_resid),
                    "r_SOLID": float(rr.r_SOLID)
                    if pd.notna(rr.r_SOLID) else np.nan,
                    "r_HAEM": float(rr.r_HAEM)
                    if pd.notna(rr.r_HAEM) else np.nan,
                    "r_exp_slc7a11_B53": float(rr.r_exp_slc7a11_B53),
                    "survivor": bool(rr.survivor)})
            prows.append(row)
    panel = pd.DataFrame(prows)
    inm = panel.emp_p_anchor_swap.notna()
    panel.loc[inm, "fdr_within_panel"] = bh_fdr(
        panel.loc[inm, "emp_p_anchor_swap"].to_numpy())
    inmb = panel.emp_p_anchor_swap_bidentity.notna()
    panel.loc[inmb, "fdr_within_panel_bidentity"] = bh_fdr(
        panel.loc[inmb, "emp_p_anchor_swap_bidentity"].to_numpy())
    inu = panel.status == "宇宙內"
    panel = panel.sort_values("r_B53", ascending=False)
    panel.to_csv(f"{RES}/205_mechanism_panel.csv", index=False,
                 float_format="%.5f")
    nonmod_panel = ~panel.gene.isin(MODULE)
    nsig = int((panel.fdr_within_panel < 0.10).sum())
    nsig_nm = int(((panel.fdr_within_panel < 0.10) & nonmod_panel).sum())
    nsig_nm_b = int(((panel.fdr_within_panel_bidentity < 0.10)
                     & nonmod_panel).sum())
    out["P9_mechanism_panel"] = {
        "n_panel_genes": int(len(panel)),
        "n_in_universe": int(inu.sum()),
        "n_fdr_below_010_all": nsig,
        "n_fdr_below_010_nonmodule": nsig_nm,
        "n_fdr_below_010_nonmodule_bidentity_partialled": nsig_nm_b,
        "pass_bidentity_partialled": bool(nsig_nm_b >= 1),
        "bidentity_note": "P2b 顯示無偏掃描去掉 B 身分之後就沒訊號了，所以"
                          "以 *_bidentity_partialled 這一組為主要判準",
        "module_members_in_panel": sorted(
            panel.gene[panel.gene.isin(MODULE)].tolist()),
        "pass": bool(nsig_nm_b >= 1),
        "top": panel[inu].head(12)[
            ["panel", "gene", "r_B53", "emp_p_anchor_swap",
             "fdr_within_panel", "r_SOLID", "r_HAEM", "depfrac_B53",
             "panessential"]].round(4).to_dict("records")}
    log("P9：小組 %d 個基因、宇宙內 %d 個；FDR < 0.10 共 %d 個，其中非模組 %d 個"
        " → %s" % (len(panel), int(inu.sum()), nsig, nsig_nm,
                   "過" if nsig_nm >= 1 else "不過"))
    log("   去 B 身分之後 FDR < 0.10 的非模組基因：%d 個" % nsig_nm_b)
    log("\n" + panel[inu].head(15)[
        ["panel", "gene", "r_B53", "emp_p_anchor_swap", "fdr_within_panel",
         "r_B53_bidentity_partialled", "emp_p_anchor_swap_bidentity",
         "fdr_within_panel_bidentity", "r_SOLID", "r_HAEM", "depfrac_B53"]]
        .round(3).to_string(index=False))

    # ---------- 輸出掃描表 ----------
    res = res.sort_values("r_B53", ascending=False)
    res.drop(columns=["perm_p_uncorrected"]).to_csv(
        f"{RES}/205_codependency_scan.csv.gz", index=False,
        float_format="%.5f")
    log("\n掃描結果寫出：results/205_codependency_scan.csv.gz（%d 個基因）"
        % len(res))

    nm = res[~res.gene.isin(MODULE)]
    cols = ["gene", "r_B53", "r_B53_strat_resid", "r_B53_partial_bidentity",
            "r_B53_full_control", "r_B53_partial_linemean",
            "r_B53_jackknife3", "r_SOLID", "r_HAEM", "depfrac_B53",
            "depfrac_all", "broadly_essential_60_89",
            "emp_p_anchor_swap", "survivor", "above_FWER"]
    log("\n=== 主錨 SLC7A11、B53 內共依賴 Top 20（非模組）===")
    log("\n" + nm.head(20)[cols].round(3).to_string(index=False))

    # ---------- 對照 ----------
    ctrl = {}
    ranks = {g: int(i + 1) for i, g in enumerate(res.gene)}
    for g in POS_CTRL + NEG_CTRL_BID + NEG_CTRL_ESS + NEG_CTRL_DEAD:
        if g in res.gene.values:
            r0 = res[res.gene == g].iloc[0]
            ctrl[g] = {
                "in_universe": True, "rank": ranks[g],
                "pct_of_universe": round(ranks[g] / len(res), 4),
                "r_B53": float(r0.r_B53),
                "r_B53_strat_resid": float(r0.r_B53_strat_resid),
                "retained_after_strat": round(
                    float(r0.r_B53_strat_resid / r0.r_B53), 3)
                if r0.r_B53 != 0 else None,
                "r_SOLID": float(r0.r_SOLID) if pd.notna(r0.r_SOLID) else None,
                "r_HAEM": float(r0.r_HAEM) if pd.notna(r0.r_HAEM) else None,
                "depfrac_all": float(r0.depfrac_all),
                "panessential": bool(r0.panessential)}
        elif g in EB.columns:
            why = []
            if not bool(complete_B[g]):
                why.append("B53 有缺值")
            if float(sd_B[g]) < MIN_SD_B:
                why.append("sd %.3f < %.2f" % (sd_B[g], MIN_SD_B))
            if float(depfrac_B[g]) < MIN_DEPFRAC_B:
                why.append("B53 依賴比例 %.3f < %.2f"
                           % (depfrac_B[g], MIN_DEPFRAC_B))
            if bool(paness[g]):
                why.append("泛必需 %.3f" % depfrac_all[g])
            zg = zrank_vec(EB[g].to_numpy(dtype=np.float64)) / np.sqrt(n_b)
            ctrl[g] = {"in_universe": False,
                       "excluded_because": "；".join(why) or "?",
                       "r_B53": float(np.dot(zg, za_B)),
                       "depfrac_B53": float(depfrac_B[g]),
                       "depfrac_all": float(depfrac_all[g]),
                       "panessential": bool(paness[g])}
        else:
            ctrl[g] = {"in_universe": False,
                       "excluded_because": "不在 CRISPR 矩陣"}
    out["controls"] = {
        "positive_module_members": {g: ctrl[g] for g in POS_CTRL},
        "negative_bcell_identity": {g: ctrl[g] for g in NEG_CTRL_BID},
        "negative_panessential": {g: ctrl[g] for g in NEG_CTRL_ESS},
        "negative_ruled_out": {g: ctrl[g] for g in NEG_CTRL_DEAD}}
    log("\n=== 對照 ===")
    for lab, gs in (("陽性（模組成員，不進名單）", POS_CTRL),
                    ("陰性 B 細胞身分", NEG_CTRL_BID),
                    ("陰性泛必需", NEG_CTRL_ESS),
                    ("陰性已否證", NEG_CTRL_DEAD)):
        for g in gs:
            c = ctrl[g]
            if c.get("in_universe"):
                log("  %-22s %-7s r_B53 = %+.3f｜層內去均值後 %+.3f"
                    "（保留 %s）｜r_SOLID %+.3f｜排名 %d/%d（前 %.1f%%）"
                    % (lab, g, c["r_B53"], c["r_B53_strat_resid"],
                       c["retained_after_strat"], c["r_SOLID"] or 0,
                       c["rank"], len(res), 100 * c["pct_of_universe"]))
            else:
                log("  %-22s %-7s 不在宇宙（%s）｜r_B53 = %s"
                    % (lab, g, c.get("excluded_because", "?"),
                       ("%+.3f" % c["r_B53"]) if "r_B53" in c else "NA"))
    # 陽性對照：GCLC／GCLM／AIFM2 三個分開判（事先寫死：各自都要在前 5%）
    pos_req = ["GCLC", "GCLM", "AIFM2"]
    pos_each = {}
    for g in pos_req:
        c = ctrl[g]
        if c.get("in_universe"):
            pos_each[g] = {"in_universe": True,
                           "pct_of_universe": c["pct_of_universe"],
                           "r_B53": c["r_B53"],
                           "in_top5pct": bool(c["pct_of_universe"] <= 0.05)}
        else:
            pos_each[g] = {"in_universe": False,
                           "excluded_because": c.get("excluded_because"),
                           "r_B53": c.get("r_B53"), "in_top5pct": False}
    n_pos_ok = sum(1 for v in pos_each.values() if v["in_top5pct"])
    out["controls"]["positive_each"] = pos_each
    out["controls"]["n_positive_in_top5pct"] = n_pos_ok
    out["controls"]["positive_all_three_pass"] = bool(n_pos_ok == 3)
    out["controls"]["positive_verdict"] = (
        "GCLC／GCLM／AIFM2 中有 %d／3 落在候選宇宙前 5%%。"
        "未達 3／3 就必須寫成『全基因掃描在 n=%d 下檢定力不足』。"
        % (n_pos_ok, n_b))
    log("  陽性對照判定：%d/3 落在前 5%%（%s）"
        % (n_pos_ok, {g: v.get("pct_of_universe") for g, v in
                      pos_each.items()}))

    # 陰性對照（B 細胞身分）：數字判準 = P3a 殘差 r 低於 FWER 門檻且非存活
    bid_ver = {}
    for g in NEG_CTRL_BID:
        c = ctrl[g]
        if c.get("in_universe"):
            r0 = res[res.gene == g].iloc[0]
            bid_ver[g] = {
                "r_B53": float(r0.r_B53),
                "r_B53_strat_resid": float(r0.r_B53_strat_resid),
                "below_FWER_threshold": bool(r0.r_B53_strat_resid < thr_fwer),
                "survivor": bool(r0.survivor),
                "collapsed": bool(r0.r_B53_strat_resid < thr_fwer
                                  and not r0.survivor)}
        else:
            bid_ver[g] = {"in_universe": False,
                          "excluded_because": c.get("excluded_because"),
                          "r_B53": c.get("r_B53"), "collapsed": None}
        out["controls"]["negative_bcell_identity"][g].update(
            {"collapsed": bid_ver[g]["collapsed"]})
    n_bid_bad = sum(1 for v in bid_ver.values() if v["collapsed"] is False)
    # 宣告的判準（P3a 殘差 < FWER 門檻）在沒有任何基因破 0.525 的情況下
    # 幾乎是自動通過的，所以它本身沒有鑑別力。真正要看的是保留比例：
    retain = {g: (bid_ver[g].get("r_B53_strat_resid", float("nan"))
                  / bid_ver[g]["r_B53"])
              if bid_ver[g].get("r_B53") else None
              for g in NEG_CTRL_BID}
    p3a_bad = [g for g, v in retain.items() if v is not None and v > 0.9]
    out["controls"]["bcell_identity_verdict"] = {
        "detail": bid_ver, "n_not_collapsed": n_bid_bad,
        "pass_declared_rule": bool(n_bid_bad == 0),
        "declared_rule_is_vacuous": True,
        "declared_rule_why_vacuous":
            "判準是『P3a 殘差 r < FWER 門檻 %.3f』，但全掃描沒有任何基因破這個"
            "門檻，所以這條規則幾乎自動通過、沒有鑑別力。" % thr_fwer,
        "P3a_retained_ratio": {g: (round(v, 3) if v is not None else None)
                               for g, v in retain.items()},
        "P3a_inadequate_genes": p3a_bad,
        "P3a_inadequate": bool(len(p3a_bad) > 0),
        "P3b_residual_r": {g: round(float(
            res[res.gene == g].iloc[0].r_B53_partial_bidentity), 3)
            for g in NEG_CTRL_BID if g in res.gene.values},
        "verdict": "3 層粗層（P3a）對 B 細胞身分基因幾乎沒有作用——%s 的 r 保留"
                   "超過 90%%。只有 P3b（偏迴歸掉 B 身分依賴分數）會讓它們塌掉，"
                   "所以譜系控制必須以 P3b／P3c 為主，P3a 不足以當證據。"
                   % "、".join(p3a_bad)}
    log("  P3a 對 B 身分基因的保留比例：%s → P3a %s"
        % ({g: (round(v, 2) if v is not None else None)
            for g, v in retain.items()},
           "不夠力" if p3a_bad else "夠力"))
    log("  P3b 之後的殘差 r：%s"
        % out["controls"]["bcell_identity_verdict"]["P3b_residual_r"])
    log("  陰性 B 身分判定：%d/%d 沒塌（判準：P3a 殘差 r < FWER 門檻 %.3f "
        "且非存活）" % (n_bid_bad, len(bid_ver), thr_fwer))

    # ---------- 候選名單 ----------
    cand = nm[nm.survivor & nm.P4_outgroup].copy()
    log("\n存活（P3a＋P3b＋P3c＋P5＋P6＋P7＋非泛必需）：%d／%d；"
        "再加 P4 外組弱：%d"
        % (int(nm.survivor.sum()), len(nm),
           int((nm.survivor & nm.P4_outgroup).sum())))
    level = "survivor+P4"
    short = cand.head(5)
    if len(short) == 0:
        level = "survivor"
        short = nm[nm.survivor].head(5)
        log("沒有基因同時通過所有過濾 → 退一步報存活者前 5 名")
    if len(short) == 0:
        level = "unfiltered"
        short = nm.head(5)
        log("連存活者都沒有 → 報未過濾前 5 名，並明寫沒通過哪一項")
    shortlist = short.gene.tolist()
    out["shortlist_level"] = level
    out["shortlist_genes"] = shortlist
    log("\n=== 候選名單（層級 %s）===" % level)
    log("\n" + short[cols + ["r_exp_slc7a11_B53", "emp_p_exp_null",
                             "delta_hi_minus_lo"]].round(3)
        .to_string(index=False))
    out["shortlist_detail"] = short[
        cols + ["r_B53_module_anchor", "r_B72", "r_exp_slc7a11_B53",
                "emp_p_exp_null", "delta_hi_minus_lo", "mean_eff_B53",
                "mean_eff_all", "sd_B53", "P3a_strat", "P3b_bidentity",
                "P3c_full_control", "P4_outgroup", "P5_partial",
                "P6_jackknife", "P7_B72"]
    ].round(4).to_dict("records")
    keepg = set(shortlist) | set(MODULE) | set(NEG_CTRL_BID) \
        | set(NEG_CTRL_ESS) | set(NEG_CTRL_DEAD)
    res[res.gene.isin(keepg)].to_csv(f"{RES}/205_shortlist_detail.csv",
                                     index=False, float_format="%.5f")
    nm.head(TOP_K_REPORT).to_csv(f"{RES}/205_top%d.csv" % TOP_K_REPORT,
                                 index=False, float_format="%.5f")
    json.dump(out, open(f"{RES}/205_next_lever.json", "w"), indent=2,
              ensure_ascii=False)
    log("中途結果已寫出（PubMed 之前）")

    # ---------- 事後分析（明確標示：以下不是預先指定的檢定） ----------
    # 預先指定的部分到這裡結束。P9 把訊號指向硒半胱胺酸／GPX4 那一臂之後，
    # 下面幾項是為了把那個觀察講清楚而補做的，所有 p 值都要當成事後、
    # 不可用來宣告顯著性。
    log("\n=== 事後分析（非預先指定）===")
    post = {"disclaimer": "post hoc, exploratory; p-values are not "
                          "pre-specified and must not be read as "
                          "significance claims"}

    # (1) 小組各臂的平均 r（看訊號集中在哪一臂）
    arm = panel[panel.emp_p_anchor_swap.notna()].groupby("panel").agg(
        n=("gene", "size"), mean_r=("r_B53", "mean"),
        max_r=("r_B53", "max"),
        n_fdr010=("fdr_within_panel", lambda s: int((s < 0.10).sum())))
    arm = arm.sort_values("mean_r", ascending=False)
    log("\n各臂平均 r_B53：\n" + arm.round(3).to_string())
    post["panel_arms"] = arm.round(4).reset_index().to_dict("records")

    # (2) 2×2 共依賴：SLC7A11 依賴的 B 株裡有幾株也依賴候選基因
    #     虛無：在 B53 內依賴比例配對（±0.03）的基因
    dep_anchor = set(EB.index[EB[ANCHOR] < DEP_CUT])
    log("\nB53 裡依賴 SLC7A11 的株：%d" % len(dep_anchor))
    # 硒相關基因組：事先寫死的小組（硒半胱胺酸合成／GPX4）＋ 無偏掃描自己
    # 撞出來的硒「攝取」路徑。LRP8（ApoER2）是 selenoprotein P 的受體，
    # LDLRAP1（ARH）是 LDL 受體家族內吞所需的接合子——這兩個不在預先指定的
    # 小組裡，是無偏掃描的 top 10 自己帶出來的，所以純屬事後觀察。
    sec_arm = ["SCLY", "SEPSECS", "SEPHS2", "PSTK", "EEFSEC", "GPX4",
               "LRP8", "LDLRAP1"]
    post["selenium_uptake_note"] = (
        "LRP8 (ApoER2) 與 LDLRAP1 (ARH) 不在預先指定的 PANEL 裡，卻自己出現在"
        "無偏掃描的非模組 top 10；它們是 selenoprotein P 攝取路徑。這是事後"
        "觀察，不是預先指定的檢定，但它表示無偏掃描與機制小組各自獨立地指向"
        "同一條硒的軸線。")
    ov = []
    for g in list(dict.fromkeys(shortlist + sec_arm)):
        if g not in EB.columns:
            continue
        dg = set(EB.index[EB[g] < DEP_CUT])
        k = len(dep_anchor & dg)
        tab = [[k, len(dg) - k],
               [len(dep_anchor) - k,
                n_b - len(dg) - len(dep_anchor) + k]]
        orv, fp = stats.fisher_exact(tab)
        fb = float(depfrac_B[g])
        pool = [h for h in genes
                if abs(float(depfrac_B[h]) - fb) <= 0.03 and h != g]
        nullk = []
        for h in pool:
            nullk.append(len(dep_anchor & set(EB.index[EB[h] < DEP_CUT])))
        nullk = np.array(nullk) if pool else np.array([np.nan])
        ep = float((nullk >= k).mean()) if pool else float("nan")
        ov.append({"gene": g, "n_dep_B53": len(dg),
                   "n_overlap_with_anchor_dep": k,
                   "odds_ratio": float(orv), "fisher_p": float(fp),
                   "n_rate_matched": len(pool),
                   "null_mean_overlap": float(np.nanmean(nullk)),
                   "empirical_p_vs_rate_matched": ep})
        log("  %-8s 自身依賴 %2d 株，與 SLC7A11 依賴者重疊 %2d 株"
            "（OR %.1f、Fisher p %.1e；配對虛無重疊 %.1f，經驗 p %.3f）"
            % (g, len(dg), k, orv, fp, np.nanmean(nullk), ep))
    post["overlap_2x2"] = ov

    # (3) 留一粗層：r 是不是只靠某一層撐起來的
    loo = []
    for g in list(dict.fromkeys(shortlist + sec_arm)):
        if g not in gidx:
            continue
        row = {"gene": g, "r_B53": float(rB[gidx[g]])}
        for s in np.unique(strat_B):
            m = strat_B != s
            row["drop_%s" % s] = float(
                np.dot(zrank_vec(EB[g].to_numpy(dtype=np.float64)[m]),
                       zrank_vec(a_raw[m])) / m.sum())
        loo.append(row)
    post["leave_one_stratum_out"] = loo
    log("\n留一粗層：\n" + pd.DataFrame(loo).round(3).to_string(index=False))

    # (3b) 把 GPX4 的依賴輪廓偏迴歸掉：候選基因是不是只是 GPX4 的影子？
    #      GPX4 用 GSH 當還原劑，所以「xCT→GSH」跟「做得出 GPX4」本來就同一條
    #      路；這一項是要看候選基因有沒有在 GPX4 之外還加上東西。
    if "GPX4" in EB.columns:
        zgpx = zrank_vec(EB["GPX4"].to_numpy(dtype=np.float64))
        ag = resid_on((za_B * np.sqrt(n_b)).reshape(-1, 1), zgpx).ravel()
        pg = []
        for g in list(dict.fromkeys(shortlist + sec_arm)):
            if g not in EB.columns or g == "GPX4":
                continue
            zg1 = zrank_vec(EB[g].to_numpy(dtype=np.float64))
            rg = resid_on(zg1.reshape(-1, 1), zgpx).ravel()
            r_pg = float(np.dot(rg, ag))
            r0 = float(np.dot(zg1 / np.sqrt(n_b), za_B))
            r_g_gpx = float(np.dot(zg1 / np.sqrt(n_b),
                                   zgpx / np.sqrt(n_b)))
            pg.append({"gene": g, "r_B53": r0,
                       "r_with_GPX4_profile": r_g_gpx,
                       "r_B53_partial_GPX4": r_pg,
                       "retained": round(r_pg / r0, 3) if r0 else None})
        post["partial_out_GPX4"] = pg
        log("\n偏迴歸掉 GPX4 依賴輪廓之後：\n"
            + pd.DataFrame(pg).round(3).to_string(index=False))

    # (4) 候選基因在 B 株裡到底有沒有表現（標的必須被表現）
    ex = []
    for g in list(dict.fromkeys(shortlist + sec_arm + [ANCHOR])):
        if g in EXP.columns:
            ex.append({"gene": g,
                       "median_logTPM_B53": float(EXP.loc[bx, g].median()),
                       "median_logTPM_all": float(EXP[g].median()),
                       "frac_B53_above_1": float(
                           (EXP.loc[bx, g] > 1).mean())})
    post["expression_of_candidates"] = ex
    log("\n候選基因表現：\n" + pd.DataFrame(ex).round(3).to_string(index=False))
    out["post_hoc"] = post

    # ---------- PubMed（排名凍結之後才查） ----------
    log("\nPubMed 查詢（排名已凍結）…")
    panel_top = panel[inu].head(8).gene.tolist()
    lit_genes = list(dict.fromkeys(shortlist + panel_top + MODULE))
    lit = {}
    for g in lit_genes:
        q = {
            "n_lymphoma": '"%s"[All Fields] AND lymphoma[All Fields]' % g,
            "n_lymphoma_or_leukemia":
                '"%s"[All Fields] AND (lymphoma[All Fields] OR '
                'leukemia[All Fields])' % g,
            "n_mcl": '"%s"[All Fields] AND ("mantle cell lymphoma"'
                     '[All Fields])' % g,
            "n_inhibitor": '"%s"[All Fields] AND (inhibitor[All Fields] OR '
                           'antagonist[All Fields])' % g,
            "n_inhibitor_smallmol":
                '"%s"[All Fields] AND inhibitor[All Fields] AND '
                '("small molecule"[All Fields] OR drug[All Fields])' % g,
            "n_ferroptosis": '"%s"[All Fields] AND ferroptosis[All Fields]'
                             % g,
            "n_total": '"%s"[All Fields]' % g}
        d = {}
        for k, term in q.items():
            d[k] = pubmed_count(term)
            time.sleep(0.40)
        lit[g] = d
        log("  %-9s lymphoma %s｜lym/leu %s｜MCL %s｜inhibitor %s｜"
            "小分子 %s｜ferroptosis %s｜total %s"
            % (g, d["n_lymphoma"], d["n_lymphoma_or_leukemia"], d["n_mcl"],
               d["n_inhibitor"], d["n_inhibitor_smallmol"],
               d["n_ferroptosis"], d["n_total"]))
    out["pubmed"] = lit
    out["pubmed_note"] = (
        "計數法與步驟 51 相同。n_inhibitor 只是『文獻裡有人寫過抑制劑／拮抗劑』"
        "的代理，不等於有臨床級化合物；沒有查 DrugBank，所以凡是沒有另外手動"
        "確認的化合物都要標為未驗證。")

    # ---------- 總結 ----------
    checks = {k: v["pass"] for k, v in out.items()
              if isinstance(v, dict) and "pass" in v}
    checks["positive_controls_all_three_in_top5pct"] = \
        out["controls"]["positive_all_three_pass"]
    checks["bcell_identity_P3a_adequate"] = not out["controls"][
        "bcell_identity_verdict"]["P3a_inadequate"]
    out["amendment_log"] = [
        "第一次執行在讀到 Gate0（模組在 B53 內 r = 0.382）與候選宇宙大小"
        "（1,613）之後中止。當時還沒有看過任何基因層級的結果（top 名單、"
        "對照、候選名單都沒看）。",
        "中止後才加進去的：P2 把模組欄位與換錨基因的自我相關剔掉、經驗 p 改"
        "成 (b+1)/(m+1)、換錨池從抽 300 改成用整池 2,000；P3b（偏迴歸掉 B "
        "細胞身分分數）與 P3c（全控制）兩個譜系控制，並寫進存活判定；陽性與"
        "陰性對照的數字判準。這些都是在看到基因層級結果之前寫定的。",
        "完全事後（看過結果才做，不可用來宣稱顯著性）：P2b、第 (1)–(4) 項"
        "事後分析、以及把 LRP8／LDLRAP1 併進硒組。"]
    out["summary"] = {
        "checks": checks,
        "n_universe": len(genes),
        "n_survivors": int(nm.survivor.sum()),
        "n_survivors_with_outgroup_filter": int(
            (nm.survivor & nm.P4_outgroup).sum()),
        "n_above_FWER_nonmodule": n_pass_p1,
        "FWER05_threshold_r": thr_fwer,
        "shortlist": shortlist, "shortlist_level": level}
    log("\n通過：%s" % (", ".join(k for k, v in checks.items() if v) or "（無）"))
    log("未通過：%s" % (", ".join(k for k, v in checks.items()
                                 if v is False) or "（無）"))
    json.dump(out, open(f"{RES}/205_next_lever.json", "w"), indent=2,
              ensure_ascii=False)
    log("DONE → results/205_*")


if __name__ == "__main__":
    main()
