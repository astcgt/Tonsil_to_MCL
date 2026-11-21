# 分析計畫 v1 — 從正常扁桃腺 B 細胞分化軸「扣除」後，找出 MCL 的異常調控網絡並設計藥物

## 0. 目標鏈（使用者需求拆解）
1. 用 HCA Tonsil Atlas 建立**正常 B 細胞分化連續軸**（NBC → 活化 → GC-committed → DZ/LZ GCBC → MBC / PC）。
2. 把 **MCL 腫瘤細胞投影**到這條正常軸上，得到每個腫瘤細胞的「正常對應狀態」。
3. **扣除正常分化成分**：殘差 = 觀測 − 由正常對應狀態預測的期望值。
4. 在殘差上做**基因調控網絡 (GRN) 推論**，排序「最異常的調控子 (regulon)」。
5. 判斷最異常者的**異常型態**（cis 調控/結構變異 vs 編碼突變 vs 純表觀劑量）。
6. 依型態選工具驗證：
   - cis 調控 → **AlphaGenome 2**（需 API key）
   - 編碼突變 → **ESM-2/ESM-C + Boltz-2/OpenFold3 (+PIONEER2 介面)**；AlphaFold3 權重受限
   - 依賴性驗證 → **DepMap CRISPR**（MCL 細胞株）
7. 找弱點 → 生成**蛋白質 binder / ASO / siRNA / decoy 寡核苷酸**候選（暫不考慮藥動）。

## 1. 資料來源
| 用途 | 來源 | 大小 | 取得方式 |
|---|---|---|---|
| 正常扁桃腺 scRNA (121 類註解) | Zenodo 8373756 `TonsilAtlasSeuratRNA.tar.gz` | 10 GB | 直接下載 ✅ |
| 正常扁桃腺 Multiome (RNA+ATAC) | Zenodo 8373756 `TonsilAtlasSeuratMultiome.tar.gz` | 5.7 GB | 直接下載 ✅ |
| **MCL 兩病人 Multiome (M102/M413)** | Zenodo 10373041 `MCL.zip` | 8.9 GB | 直接下載 ✅ |
| MCL 復發/原發 scRNA 擴充世代 | researchdata.se 2024-497；Zenodo 14626636 | 待查 | 補強 n 數 |
| MCL 體細胞突變 | cBioPortal `mcl_idibips_2013` (29 例) | 小 | API ✅ |
| MCL 細胞株依賴性 | DepMap CRISPR (figshare 鏡像) | ~0.5 GB | 待試 |
| TF–target 先驗 | CollecTRI / DoRothEA / JASPAR2024 | 小 | 直接下載 |

## 2. 方法設計（關鍵：如何「扣除正常分化」）
三種互補做法，取交集以避免單一方法的偽陽性：

- **(A) kNN 參考匹配殘差**：把 MCL 細胞映射到正常參考的潛在空間，取 k 個最近的正常鄰居，期望表現 = 鄰居平均；殘差 = 觀測 − 期望。
- **(B) 迴歸殘差**：`expr ~ s(正常分化 pseudotime) + 細胞週期 + 病人批次`，取 Pearson/Deviance 殘差。
- **(C) 潛在空間解耦**：scVI/scArches 參考映射，用 decoder 重建「若此細胞是正常的」表現。

GRN 層：
- 以 **decoupler + CollecTRI** 算每細胞 TF 活性（快、穩健），再用 **pySCENIC/GRNBoost2** 在腫瘤細胞內重建 de novo 調控網絡做交叉驗證。
- 「異常度」定義 = (MCL TF 活性 − 匹配正常期望 TF 活性) 的效應量，並以病人為單位做混合模型避免單病人驅動。
- 以 Multiome ATAC 檢查該 TF 的 **motif 可及性** 與 **peak–gene link**，判斷是 cis 還是 trans。

## 3. 預期（需用資料證實，不可先射箭）
MCL 的定義性病變是 **t(11;14)(q13;q32) → IGH 增強子劫持 CCND1**，屬 cis 調控異常；另 **SOX11** 在傳統型 MCL 異常高表現且不屬正常成熟 B 細胞程式。若殘差分析排名前列出現 CCND1/E2F、SOX11，即與已知生物學自洽，可作為方法學的 positive control；真正有價值的是**排名次高但未被充分研究**的調控子。

## 4. 運算資源規劃
- CPU：AMD EPYC 7742 ×2（256 執行緒），分析主要以 CPU 平行。
- GPU：**只用 GPU 3 與 GPU 4**（GPU 0–2 依指示避開；GPU 5/6 已被佔滿）。
- 記憶體：系統 1.0 TB，可用 ~950 GB；單步驟上限設 300 GB。
- 儲存：`/raid` 剩 2.2 TB（已 84%）→ 本專案硬上限 **150 GB**；暫存走 scratchpad。

## 5. 已知阻礙（需使用者決策）
1. **AlphaGenome 2**：2026-11-11 發布，僅 API（非商業免費），需使用者申請 API key。無 key 則改用 Enformer/Borzoi（開源權重）做 cis 效應預測。
2. **AlphaFold3 權重**：需向 DeepMind 申請（數週）。替代：Boltz-2（MIT，開源權重）、OpenFold3、Chai-1。
3. **PIONEER2**：僅網頁伺服器 (pioneer2.yulab.org)，無批次 API → 手動查詢少量殘基。
4. **MCL n=2**：樣本數不足，需補外部世代。
