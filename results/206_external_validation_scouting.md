# 206：外部驗證資料盤點

**日期**：2026-11-27
**要驗的主張**：`SLC7A11` 附近五個候選增強子元件在 MCL 腫瘤細胞中新生開啟、在分化配對的正常 B 細胞中不開啟；
`MAFG` 啟動子的 MARE 元件具功能性。
**弱點**：主要證據只來自兩名病人。
**前提**：步驟 180 已證明表現量 eQTL 在此處**沒有檢定力**（不是負結果），不再重做。

---

## 零、一段話結論

**我們在等的記憶 B 細胞 ATAC caQTL 已經釋出了，不在 eQTL Catalogue，而在作者自己的 Hugging Face 倉庫，
CC-BY-4.0，今天可下載。本盤點已經把它查完了：結果是沒有檢定力的空值，和步驟 180 同一類。**
記憶 B 細胞 chr4 上只有 2.01% 的 peak 有顯著 caQTL，五個元件的期望命中數是 **0.100**，
實際命中 0——這是基準率，不是證據。TenK10K 到此為止，不必再等、也不必再查。

**真正該去拿的是 BLUEPRINT 的 5 名原發 MCL 病人**（六種組織蛋白修飾、公開、GRCh38、peak 檔每個 1–4 MB，
合計約 60 MB，不需申請），加上同一條管線的分化配對正常記憶 B 細胞。
那是唯一能把「兩名病人」擴充成「五名病人對配對正常」的東西。

**沒有任何現存資料能回答「這五個元件是否調控 `SLC7A11`」。**
ENCODE 的 Flow-FISH CRISPR 篩選確實存在（353 個實驗），B 系細胞也做過
（GM12878 6 個、BJAB 12 個），但這 **18 個 B 系篩選全部指向 `PPIF`**；
以 `SLC7A11` 或 `MAFG` 為 examined locus 的實驗數都是 **0**。

---

## 一、先修正三個座標事實

這三點會讓後續任何「有沒有覆蓋」的判讀出錯，先釘住。

**（一）`SLC7A11` 在負股上。** Ensembl `ENSG00000151012`：chr4:138,164,097–138,242,541，
**strand = −1**（GRCh38）；MANE `ENST00000280612.9` 的 TSS 在高座標端 138,242,349。因此：

- 標為 `-19.1 kb` 的元件實際座標是 **chr4:138,222,944–138,223,474**，它落在**基因體內部（內含子）**，
  是從 TSS 往轉錄方向走 19.1 kb，**不是上游**。
- 四個標為 `+99.6` 至 `+164.0 kb` 的元件（chr4:138,341,498–138,406,449）座標高於 TSS，
  以轉錄方向論**才是上游**。
- 要檢查的 chr4 視窗是 **138.20–138.42 Mb**。若用 138.3–138.5 Mb 會**漏掉內含子那個元件**。

**（二）`MAFG` 也在負股上。** Ensembl `ENSG00000197063`：chr17:81,918,270–81,927,769，**strand = −1**；
MANE `ENST00000357736.9` 的 TSS 在 81,927,735。
MARE 位點 chr17:81,926,966 的座標低於 TSS 769 bp，所以它是 TSS **下游** 769 bp，
落在基因 5' 端內部，**不是上游**。承載它的元件是 chr17:81,926,271–81,926,994。
「啟動子近端」成立（769 bp），「上游」不成立。
另有替代轉錄本 `XM_047436069.1` 起點在 81,926,855，相對於**它**，MARE 才算上游 111 bp。
此外 `MILIP` 的基因體從 81,927,829 開始，就在這個元件上方。

**（三）hg19 對應座標**（Ensembl REST GRCh38→GRCh37，供 EpiMap 等舊建構使用）：

| 元件 | hg38 | hg19 |
|---|---|---|
| 內含子（−19.1 kb） | chr4:138,222,944–138,223,474 | chr4:139,144,098–139,144,628 |
| +99.6 kb | chr4:138,341,498–138,342,345 | chr4:139,262,652–139,263,499 |
| +105.9 kb | chr4:138,348,035–138,348,490 | chr4:139,269,189–139,269,644 |
| +150.5 kb | chr4:138,392,601–138,393,029 | chr4:139,313,755–139,314,183 |
| +164.0 kb | chr4:138,406,192–138,406,449 | chr4:139,327,346–139,327,603 |
| MAFG MARE 元件 | chr17:81,926,271–81,926,994 | chr17:79,884,147–79,884,870 |
| MARE motif 位點 | chr17:81,926,966 | chr17:79,884,842 |

---

## 二、任務一：等待中的兩個染色質 QTL 資料集

### 2.1 eQTL Catalogue：仍然完全沒有釋出

登記本身是真的，而且比步驟 180 的註記更精確。這兩個資料集不在主要命名空間裡，
而在另一套 `QTS1000xx` / `QTD1000xx` 染色質 QTL 命名空間，只列在 GitHub 上的
`eQTL-Catalogue/eQTL-Catalogue-resources` → `data_tables/chromatin_QTLs/chromatin_dataset_id_map.tsv`：

| dataset | study | 細胞 | 檢測 | n | 基因型 |
|---|---|---|---|---|---|
| QTD100034 | QTS100012 TenK10K | memory B cell | **ATAC** | 913 | wgs |
| QTD00068 | QTS100014 FG_multiome_L2 | memory B cell (L2) | **ATAC** | 1108 | finngen_imputed |

這張表目前 8,021 bytes，與本專案 `data/aux/eqtlcat/chromatin_map.tsv` **位元組完全相同**，
也就是今天早上到現在沒有變動。最後一次改動是 **2026-11-25**（「Add tissue labels to the dataset_id_map」），
所以這個分支是活的、正在維護，只是還沒有檔案。

**FTP 查核（全部 404）**：`/chromatin/`、`/chromatin_QTLs/`、`/caQTL/`、`/hQTL/`、`/r8_beta/chromatin/`、
`/susie/QTS100012/`、`/sumstats/QTS100012/`、`/r8_beta/susie/QTS100012/`、`/r8_beta/sumstats/QTS100012/`、
`/r8_beta/susie/QTS100014/`。整個 `QTS1000xx` 命名空間在 FTP 上不存在。

已釋出的部分：`/susie/`（r7）、`/r8_beta/`（16 個 study，QTS000015–000054，2026-09-09／10 更新）、
`/imported/`（只有 GTEx_V8、MetaLCL）。`tabix_ftp_paths.tsv` 共 758 個資料集，
量化方法為 ge 280、exon 113、tx 113、txrev 113、leafcutter 113、microarray 25、aptamer 1。
`Data_access` 頁另外連到 **Zenodo 22846842**（r8-beta 合併 credible sets），
檔案為 ge 779 MB、majiq 607 MB、txrev 325 MB、exon 294 MB、tx 284 MB、leafcutter 229 MB
——**全部是表現量與剪接，沒有任何 ATAC／caQTL 檔案**。
（前兩個超過 500 MB 上限，只記 URL 與大小，不下載。）

另外兩件事：

1. eQTL Catalogue 的 **RESTful API 已經停用**（回傳「This API is no longer available」）。以後只能走 FTP 與 GitHub 的表。
2. TenK10K 與 FinnGen 的**表現量**分支（QTS000064 phase 1、QTS000066）列在
   `dataset_metadata_upcoming.tsv`，標為 `ingested`，但 FTP 上同樣 404。
   它們的量化方法是 `ge`，**即使釋出也只是把步驟 180 用更大的 n 重做一次，不是染色質層的檢定**。

### 2.2 TenK10K 的 caQTL 已由作者自行釋出，而且**已經查完了**

**位置與授權**

- **Hugging Face `anglixue/TenK10K_multiome`**，公開、**CC-BY-4.0**、3,979 個檔案、
  最後更新 **2026-11-26**。
  作者的 GitHub（`powellgenomicslab/tenk10k_phase1_multiome`）明確指認這裡是正式發布位置：
  「We released all key summary results on the Hugging Face platform」。
- **論文**：Nature，DOI `10.1038/s41586-026-11173-4`（2026-10-30 線上）。
- **原始資料**另走 EGA 管制存取：`EGAS50000001653`（單細胞）、`EGAS50000001654`（WGS）。
- **內容**：922 捐贈者、28 種免疫細胞類型的 caQTL summary statistics、colocalization、
  fine-mapping、ChromBPNet 變異效應分數、cell-state 依賴 caQTL、GRN。
- **README 的三則更新**：2026-06-27 修正 caQTL×eQTL coloc 的 allele flipping；
  2026-10-30 註明 `scDeepID`／`CASTIE` 仍在開發；2026-11-16 更新 cell-state interaction 結果。
  引用 coloc 或 cell-state 結果必須用最新版。

**相關檔案與實測大小**

| 檔案 | 大小 | 說明 |
|---|---|---|
| `caQTL_common_TensorQTL/B_memory/TenK10K.sig_cis_qtl_pairs.chr4.csv` | 3,843,273 B | chr4 每個受測 peak 的 permutation 結果（含 qval）；14,998 個 peak |
| 同上 `.chr17.csv` | 3,949,717 B | chr17，15,389 個 peak |
| `Fine_mapping/all_celltypes_caQTL_susie_annotated.zip` | 194,764,002 B | 28 細胞類型的 SuSiE credible sets |
| `SMR/smr_sig_results_caQTL2eQTL.csv.zip` | 6,563,516 B | caQTL→eQTL 的 SMR 顯著結果，81,550 列 |
| `caQTL_common_TensorQTL/TenK10K.caQTL.1Mb.final_summary_qval005_feb02_2026_freeze.csv` | 66,246,863 B | 跨細胞類型摘要 |
| `caQTL_common_TensorQTL/B_memory/TenK10K.cis_qtl_pairs.chr4.parquet` | 2,328,809,720 B | 全部 nominal pairs，**超過 500 MB 上限，未下載**；parquet 支援 HTTP range |
| `ChromBPNet/B_memory/.../chrombpnet_nobias.h5` | — | 記憶 B 細胞 ATAC 訓練的序列→可及性模型，CPU 可跑 |

**查核結果（已完成）**

六個元件在記憶 B 細胞中都有對應的受測 peak：

| 元件 | 記憶 B 細胞中的受測 peak | caQTL qval |
|---|---|---|
| 內含子 −19.1 kb | chr4:138,222,946–138,223,873 | 0.9996 |
| +99.6 kb | chr4:138,341,549–138,342,501 | 0.983 |
| +105.9 kb | chr4:138,347,756–138,348,706 | 0.9996 |
| +150.5 kb | chr4:138,392,329–138,393,367 | 0.9996 |
| +164.0 kb | chr4:138,405,189–138,407,318 | 0.9996 |
| MAFG MARE 元件 | chr17:81,926,271–81,927,283 | 0.9999 |

**判決：沒有檢定力的空值，與步驟 180 同一類。**

- chr4 上 14,998 個受測 peak 裡，qval < 0.05 的有 **301 個（2.01%）**。
  五個元件的**期望命中數是 0.100**。實際命中 0——這完全在基準率之內，不構成任何證據。
  chr17 的基準率是 374／15,389（2.43%），MARE 元件同樣是空值。
- `SMR/smr_sig_results_caQTL2eQTL.csv` 共 81,550 列，
  **`SLC7A11`／`ENSG00000151012`／`MAFG`／`ENSG00000197063` 的命中數都是 0**。
  也就是 TenK10K 自己的 caQTL→eQTL SMR 在 28 種細胞類型裡都沒有把任何 peak 連到這兩個基因。
- 因此 **`Fine_mapping/all_celltypes_caQTL_susie_annotated.zip` 不是步驟 180 的替換件**：
  qval 全部 ≥ 0.98，這五個 peak 在 `B_memory` 裡根本沒有 credible set 可供檢定。
  下載它不會產生新資訊。

**關於「受測 peak」能推論什麼——不能推太多**

- peak 集合**不是記憶 B 細胞專屬的**。五個元件的 peak 區間在 `CD14_Mono`（單核球，21,500 個 chr4 peak）
  裡**以完全相同的邊界出現**。所以正確的說法是
  「這五段通過了 TenK10K 在記憶 B 細胞的量化與檢定門檻」，
  **不是**「它們在記憶 B 細胞中被呼叫為新生 peak」。
- peak 區間也比我們的元件寬得多（例：+164.0 kb 的 peak 2,129 bp 對我們的 258 bp），
  這是合併式共識 peak 集的特徵，不是元件層解析度。
- TenK10K 是**血液**記憶 B 細胞的 pseudobulk，不是本專案用的扁桃體狀態配對參考。
- 這些觀察**不要求**修改本專案的主張。本專案的主張本來就是量的
  （`normal_log2cpm` 0.183–1.497、`ratio` 2.21–13.0），不是有／無。
  只是在對外敘述時，「新生開啟」必須明確寫成量的增加。

### 2.3 FinnGen multiome（n = 1,108）：已發表，取得條件**未解決**

- **論文**：Nature，DOI `10.1038/s41586-026-11078-2`（2026-09-30 線上）。
  1,108 名芬蘭人、約 1,000 萬 PBMC 配對 snATAC + snRNA，
  得到 **338,100 個 cis-caQTL、涵蓋 210,584 個 peak**，以及 51,083 個 cis-eQTL。
- **瀏覽器**：`https://cascade.finngen.fi/`——線上、公開可達（HTTP 200）、不需登入，
  但是 JavaScript 應用，純 HTTP 取不到區間查詢結果。
- **狀態：未解決。** 找不到 summary statistics 的大量下載位置，也找不到明文取得條款；
  medRxiv 預印本（`10.1101/2025.11.25.25340489`）拒絕自動抓取（HTTP 403）。
  也未確認 memory B 是否在其 Level 2 標註中單獨列出。
- **下一步（約二十分鐘）**：用瀏覽器開 CASCADE 查 `chr4:138,222,944-138,406,449`
  與 `chr17:81,926,271-81,926,994`，並讀 Nature 論文的 Data Availability。
  這是 n 最大的那一份，值得確認——但考慮 TenK10K 在 n = 913 下已是基準率空值，
  期望值不高。

---

## 三、任務二：覆蓋座標、依決定性排序

### 1. BLUEPRINT 原發 MCL 對配對正常 B 細胞的組織蛋白 peak（公開，GRCh38）

**判決：唯一直接攻擊「只有兩名病人」這個弱點的資料。第一順位。**

- **5 名原發 MCL 病人**：`007MCL`、`2392MCL`、`6902MCL`、`817-03MCL`、`8761MCL`
  （`BIOMATERIAL_TYPE = Primary Tissue`，venous blood）。
  每人都有 **H3K27ac、H3K4me1、H3K4me3、H3K27me3、H3K9me3、H3K36me3、Input**，
  外加 total-RNA-Seq 與 WGBS。另有 MCL 細胞株 JVM-2、Z-138 的同一組標記（不要用，見 README 的細胞株保留）。
- **分化配對的正常 B 細胞**（同一條管線、同一個建構）有 H3K27ac narrow peak 的樣本數：
  naive B 8、CD38-negative naive B 4、**class switched memory B 3**、**unswitched memory B 1**、
  germinal center B 3，共 19 個。
- **索引**：`https://ftp.ebi.ac.uk/pub/databases/blueprint/releases/current_release/homo_sapiens/20160816.data.index`（9.3 MB，11,277 列）
- **檔案路徑**：`https://ftp.ebi.ac.uk/pub/databases/blueprint/data/homo_sapiens/GRCh38/venous_blood/<donor>/Mantle_Cell_Lymphoma/ChIP-Seq/NCMLS/`
- **大小**：MACS2 peak BED 每個 **0.8–4 MB**，全部加起來約 60 MB。實測可公開取得：
  `S01FJZH1.ERX1122580.H3K27ac.bwa.GRCh38.20160126.bed.gz` → HTTP 200、1,654,269 bytes；
  `S015BHH1.ERX1007488.H3K27ac.bwa.GRCh38.20151026.bed.gz`（csMBC pool 2）→ HTTP 200、1,227,338 bytes。
- **取得方式**：直接 HTTP，不需申請。
- **授權**：processed data 依 BLUEPRINT `README_data_reuse_20140721` 公開；raw BAM 受管制（EGA）。
- **保留**（必須寫進任何結論）：
  - H3K27ac／H3K4me1 量的是增強子**活性**，不是**可及性**。不是與本專案 ATAC 同一個檢測。
  - MCL 樣本取自**靜脈血**（白血病性 MCL），本專案的腫瘤來自組織。
  - **開放 FTP 上沒有 MCL 的 ATAC／DNase**：`007MCL/Mantle_Cell_Lymphoma/` 底下只有
    `Bisulfite-Seq/`、`ChIP-Seq/`、`RNA-Seq/` 三個目錄。可及性只能走下面第 2 項的管制申請。
  - peak call 是二元的，所以檢定是「每個元件在 5 個腫瘤對 4 個記憶 B 裡有沒有 peak」，檢定力低但不為零。
    **做之前先算基準率**（全基因體 peak 覆蓋率 × 5），否則會重蹈 2.2 節的覆轍。
  - **必須先用同一批捐贈者的 RNA-Seq（gene quant 檔就在索引裡）確認這 5 名病人表現 `SLC7A11`**。
    若不表現，任何 null 都無法解讀——這正是 README 記錄過的 DepMap 細胞株陷阱。
  - **獨立性**：本專案兩例（M102、M413）與 BLUEPRINT 這 5 名捐贈者**是否重疊尚未查核**。
    在寫「獨立複製」之前要先確認。
  - bigWig 訊號檔存在但每個 **0.7–1.6 GB**（例：Z-138 H3K27ac 1,230 MB、csMBC pool 8 H3K27ac 1,602 MB），
    **超過 500 MB 上限，不要下載**；需要定量訊號時用 HTTP range 讀取區間。

### 2. BLUEPRINT MCL ATAC-seq，EGAD00001002918（管制存取）

**判決：檢測完全對位的同一個檢定，但要申請。與第 1 項並行啟動。**

- `EGAD00001002918`：「BLUEPRINT September 2016, ATAC-seq Mantle Cell Lymphoma from venous blood, on Genome GRCh38」，
  **5 個樣本**、5 runs、5 alignments，釋出日 2016-11-21。
- 另有 `EGAD00001002655`：「BLUEPRINT ChIP-Seq from two mantle cell lymphoma patients」。
- 相關 study：`EGAS00001001637`（MCL 450K 甲基化，對應 Queirós 2016 Cancer Cell，82–86 例）。
- **取得方式**：BLUEPRINT Data Access Committee 申請（`blueprint-dac@ebi.ac.uk`）。
  EGA 的資料集紀錄未公布檔案大小。
- **工作量**：申請數週；核准後要自己從 BAM 跑 peak calling。
- **限制**：受管制、有同意書條款，不能把個體層資料放進公開產出。

### 3. ENCODE-rE2G，記憶 B 細胞（公開，GRCh38）——**已查完**

**判決：座標覆蓋。五個 `SLC7A11` 元件全部 0 條連結；MARE 元件有部分真實支持。**

- 資料集：`ENCSR300MTO`（記憶 B，捐贈者 ENCDO352FHC）、`ENCSR228RAU`（ENCDO512FEI），
  Jesse Engreitz lab，released。另有 naive B、activated B、B cell、GM12878/12891/12892/12864/12865。
- **大小**：thresholded element-gene links **3.7–6.5 MB**；candidate element 參考 BED 1.3 MB；
  完整預測 379–485 MB；bigInteract 255–267 MB（可區間讀取）。
- **五個 `SLC7A11` 元件**：在三個記憶 B 檔案（`ENCFF754SST`、`ENCFF324LSG`、`ENCFF714FCF`）
  與 GM12878（`ENCFF162ADG`）中**全部 0 條 thresholded 連結**。
- **而且連候選元件都不是**：記憶 B 細胞的 candidate element 參考檔（`ENCFF096OFJ`，1.3 MB）裡，
  五個元件各自重疊 **0 個**候選元件；只有 MARE 元件重疊 1 個。
  同時 `SLC7A11` 本身是被建模的目標——在記憶 B 細胞中它連到自己的啟動子
  （chr4:138,241,658–138,242,599，score 1.000）與 +935 kb、+934 kb、+1.05 Mb、+206 kb 四個遠端元件。
  所以方法找得到 `SLC7A11` 的調控元件，只是我們這五個連候選名單都沒進。
- **與 TenK10K 的差異要照實寫**：ENCODE 說這五段在記憶 B 細胞裡不夠開放、進不了候選名單；
  TenK10K（922 捐贈者匯總）把它們量化並檢定了。
  最可能的解釋是**閾值與樣本數差異**——
  也就是「開」或「關」在這裡是**閾值依賴**的。不應該選一邊講。
- **MARE 元件的結果分兩個捐贈者，要分開報**：
  - `ENCFF714FCF`（捐贈者 ENCDO512FEI）：元件 chr17:81,926,539–81,927,039
    （**含 MARE motif 81,926,966**），`class = genic`、**`isSelfPromoter = FALSE`**、
    → `MAFG` **score 0.657**，另連 `SIRT7` 0.342。
    這是一條**非自身啟動子**的元件—基因連結，是 MARE 元件調控 `MAFG` 的真實外部支持。
  - `ENCFF754SST`（捐贈者 ENCDO352FHC）：元件較寬（chr17:81,926,534–81,928,781），
    `class = promoter`、**`isSelfPromoter = TRUE`**、score 0.99999。
    self-promoter 連結的分數**是按定義給的**，這一條不能當證據。
  - GM12878（`ENCFF162ADG`）：同一元件連到 `MAFG` 0.920、`MAFG-AS1` 0.907、
    `PYCR1` 0.414、`SIRT7` 0.312、`PCYT2` 0.202（不含 `MILIP`）——目標基因指派不唯一。
  - **結論**：兩個記憶 B 捐贈者中有一個給出 0.657 的非自身啟動子 `MAFG` 連結。
    這支持「MARE 元件是 `MAFG` 的調控元件」，但它是在**正常**記憶 B 細胞裡成立，
    而且另一個捐贈者的檔案無法用。支持是部分的、捐贈者依賴的。

### 4. EpiMap（Boix 2021，hg19，公開）——**已查完**

**判決：只在內含子元件上有弱訊號，而且與第 3 項矛盾。MARE 的訊號不成立。**

- 檔案：`https://personal.broadinstitute.org/cboix/epimap/links/pergroup/links_by_group.hsc_b_cell.tsv.gz`（18.8 MB）、
  `links_by_group.lymphoblastoid.tsv.gz`（16.5 MB）。833 個 epigenome，建議使用 `links_corr_only`。

| 元件 | HSC & B-cell | Lymphoblastoid |
|---|---|---|
| 內含子 −19.1 kb | → `SLC7A11` 0.158、0.114 | → `SLC7A11` 0.409、0.342 |
| +99.6 / +105.9 / +150.5 / +164.0 kb | 全部 0 | 全部 0 |

- EpiMap 在該視窗中指向 `SLC7A11` 的連結有 70 條（B cell）／37 條（LCL），所以不是覆蓋缺口。
- **內含子元件的訊號很弱，而且與 ENCODE-rE2G 相反**（rE2G 在三個記憶 B 檔案都給 0）。
  EpiMap 的分數 0.11–0.16（LCL 0.34–0.41），而同檔第一列的分數是 0.017；
  分數定義（link 在該組 biosample 中出現的比例）未讀原始定義文件。
  `hsc_b_cell` 是 **HSC 加 B 細胞的混合組，不是正常記憶 B 細胞**。
  只能寫成「EpiMap 在一個混合的 B 系組別中給了內含子元件一個弱連結，而 rE2G 不同意」。
- **MARE 的 EpiMap 訊號不成立**：該連結區間是 hg19 chr17:79,883,736–79,884,156，
  而 MARE motif 位點對應 hg19 **79,884,842**。
  該區間與我們 723 bp 的元件只重疊**底部 9 bp**（79,884,147–79,884,156），
  **完全不含 motif**，比 motif 低 686 bp。不可引用。

### 5. ReMap 2022 轉錄因子 ChIP-seq 圖譜（hg38，公開，可免下載區間查詢）

**判決：六個元件全部覆蓋。提供「我們提名的因子會結合某些元件」的外部證據，但特異性未經檢定。**

- 取得方式：UCSC REST `getData/track?genome=hg38;track=ReMapTFs`，逐區間查詢，不必下載。
- 各元件 peak 數：內含子 151（80 個因子）／+99.6 kb 134（73）／+105.9 kb 324（140）／
  +150.5 kb 130（67）／+164.0 kb 365（124）／MAFG MARE 元件 1,937（**546**）。
- **MARE 家族相關命中**：
  - **+164.0 kb：`BACH1`，細胞為 GM12878**（B 系）。另有 JUNB（GM12878）、ATF3／FOSL1／NFE2（K-562）。
  - 內含子：`MAFK`（OCI-Ly7，B 細胞淋巴瘤株）、`NFE2L2`（A-549）、JUN／JUNB／JUND。
  - +150.5 kb：`JUNB`（GM12878）、`NFE2`（K-562）、FOS（MV4-11）。
  - +105.9 kb：FOS／FOSL1／FOSL2／JUN／JUNB／JUND／ATF3／`NFE2L2`，但**全部是非 B 系細胞**。
  - +99.6 kb：**沒有任何 MARE 家族因子**。
  - MAFG MARE 元件：`BACH1`（GM12878、K-562、HepG2、A-549…）、`BACH2`（B-cell、DOHH2、OCI-Ly7）、
    `MAFF`（GM12878）、`MAFG`（K-562）、`MAFK`（OCI-Ly7、HeLa、HepG2…）、`NFE2L2`、`JUNB`（GM12878）。
    B 系 biotype 共 14 種。
- **BACH1 的位置兩個資源不一致**：ReMap 把 GM12878 的 BACH1 peak 放在 **+164.0 kb**；
  UCSC `encRegTfbsClustered`（ENCODE TF ChIP-seq clusters）放在
  **+150.5 kb**（chr4:138,392,629–138,392,985，score 220）。兩者都要報。
- **必須一起寫出來的保留**：MARE 那個視窗有 1,937 條 peak、546 個不同因子。
  那是一個 ChIP-seq 熱點／高可及性啟動子，在那裡找到 MARE 家族因子**幾乎是必然的**。
  要主張特異性，必須先做寬度與可及性配對的背景檢定——又是 README 記錄過三次的同一個陷阱。
- **ReMap 裡沒有任何 MCL 細胞株**。出現的 B 系 biotype 是 GM128xx（LCL）、OCI-Ly1/3/7、
  SU-DHL-2/4、DOHH2、Raji，全是 DLBCL／Burkitt／LCL。

### 6. ENCODE4 cCRE 登記（hg38，公開，免下載）

**判決：六個元件全部落在已登記的 cCRE 上。只回答「這在人體某處是候選元件」。**

- 例：+150.5 kb = `EH38E3605851`，class **Distal enhancer**，
  DNase maxZ 4.97、H3K27ac maxZ 2.80、H3K4me3 1.52、CTCF 1.54。
  +164.0 kb = `EH38E3605860`；內含子 = `EH38E3605794`／`EH38E3605795`；
  MARE 元件 = `EH38E3251121`／`EH38E1892757`。
- 取得方式：UCSC REST `track=cCREregistry`。
- 與第 3 項的記憶 B 候選元件檔不衝突——cCRE 登記是跨全部生物樣本的聯集，
  記憶 B 的候選名單是細胞類型專屬的。
- 用途：只適合當成「這些座標不是隨機序列」的背景陳述，不能當驗證。

### 7. 4D Nucleome，B 系 3D 基因體（公開）

**判決：覆蓋座標（全基因體），但回答更弱的問題。**

- GM12878 共 89 個資料集：in situ Hi-C 27、ChIP-seq 12、CUT&RUN 11、in situ ChIA-PET 7、
  TCC 7、PLAC-seq 7、Dilution Hi-C 5、ChIA-Drop 4、sci-Hi-C 3。
  另有 lymphoblastoid 56 個（Dilution Hi-C 28、in situ Hi-C 25）。
- 能回答：四個遠端元件在 B 系淋巴母細胞核內是否與 `SLC7A11` 啟動子接觸。
- 不能回答：接觸不等於調控；GM12878 是 EBV 轉化的 LCL，既不是記憶 B 也不是 MCL。
- 取得方式：`data.4dnucleome.org`，公開。接觸矩陣很大，應以區間查詢取用。

### 8. GSE97541（Notch B 細胞 regulome）——第二套原發 MCL H3K27ac

**判決：真的是原發 MCL，但格式與建構都有代價。排在 BLUEPRINT 之後。**

- **3 例原發 MCL**：`MCL-M06`、`MCL-M08`、`MCL-M10`，各有 H3K27ac 與配對 WCE input。
  另有 3 例原發 CLL、MCL 細胞株 HBL2／JVM2／Maver-1／Mino／Rec-1／SP-49／Z138 的 H3K27ac，
  以及 Rec-1 的 4C-seq。
- **processed data 只有 `.tdf`**（IGV 二進位格式），**每個約 208 MB**
  （實測 `GSM2760324_MCL-M06_H3K27ac.tdf` = 208,413,333 bytes）。
  六個原發 MCL 檔案合計約 1.2 GB——**不要全下，一次一個**。
- raw data 在 **dbGaP `phs000939`**，受管制。
- 其他代價：`.tdf` 要轉換才能用；序列中**沒有配對的正常 B 細胞**；建構需先確認（很可能是 hg19）。

---

## 四、明確不覆蓋座標、或根本不存在的

這些查過了，列出來是為了不必再查。

- **ENCODE 的增強子擾動實驗存在，但從來沒有指向我們的基因。**
  物件型別是 `type=FunctionalCharacterizationExperiment`（共 869 個），
  assay 包括 **Flow-FISH CRISPR screen 353**、pooled clone sequencing 134、MPRA 125、
  proliferation CRISPR screen 78、Control MPRA 54、STARR-seq 42、FACS CRISPR screen 25、
  perturbation followed by snATAC-seq 8、CRISPR screen 5。
  但 `examined_loci.gene.symbol=SLC7A11` → **0**，`=MAFG` → **0**。
  B 系細胞確實做過：GM12878 共 17 個（MPRA 9、**Flow-FISH CRISPR 6**、snATAC 擾動 1、pooled clone 1）、
  **BJAB 12 個（全部 Flow-FISH CRISPR）**——
  但這 **18 個 B 系 CRISPR 篩選全部指向 `PPIF`**。
  GM12878 的 9 個 MPRA 沒有標註 examined_loci 或 examined_regions，
  要確認是否碰到我們的座標需逐個查其 `elements_references`（**尚未做**，是唯一剩下的缺口）。
  **準確的說法是：能力與細胞類型都存在，只是沒有人把它指向這個基因座。**
- **ENCODE 沒有任何 mantle cell lymphoma 生物樣本**（搜尋「mantle cell lymphoma」→ total 0）。
- **GeneHancer Double Elite**：內含子元件重疊 `GH04J138215` → `SLC7A11`（score 5.8）；
  **四個遠端元件什麼都不重疊**；MARE 元件重疊 `GH17J082095` → **`HGS`，不是 `MAFG`**。
  授權：GeneCards 僅限學術使用、不可再散布，不適合放進公開產出。
- **GEO 裡的 MCL ATAC-seq 全是細胞株**：GSE122510（Z-138，2 樣本）、GSE305144（MINO、CCMCL1）、
  GSE173869（Granta-519）、GSE230394（其實是 DLBCL）。**GEO 沒有原發 MCL 的 ATAC-seq。**
- **Pirtobrutinib MCL 單細胞多體學**（Yan 等，Am J Hematol 2026;101(7):1632–1649，PMID 41998849）
  在縱向 MCL 病人樣本上做了 scATAC-seq，但**找不到寄存編號**：
  GEO 搜尋「pirtobrutinib」只回傳 GSE267421（一個 DLBCL 系列），文章不在 PMC，
  Wiley 拒絕自動抓取（HTTP 403）。最便宜的路是寄一封信給通訊作者。
- **ENCODE 的 CLL enhancer-gene 連結檔（`ENCSR022QAQ` / `ENCFF488BCT`）本盤點不採用。**
  它不是 rE2G，而是 Guigó lab 的 ABC 式輸出，欄位結構不同（第 10 欄才是 Score、
  只含 `class=enhancer` 的連結、不含 self-promoter），要用必須依它自己的 header 重新解析。
- **eQTL Catalogue 的染色質分支**：見第二節。

**順帶一提（當文獻讀，不當驗證資料）**：GSE305144
「Transcriptional reprogramming in mantle cell lymphoma promotes BTK inhibitor resistance」（63 樣本）
與本專案的軸線直接相撞——**BRG1 ATPase domain 突變透過限制 ROS 與游離鐵，
保護 MCL 細胞免於 BTKi 誘發的鐵死亡**，機制走 ATF4／UPR。
但它的 ATAC-seq 與 BRG1／H3K27ac CUT&RUN 只做在 **MINO 與 CCMCL1 細胞株**上，
原發 MCL 細胞只做 RNA-seq。考慮 README 已記錄「DepMap 三株 MCL 都把 `SLC7A11` 表現在中位數以下」，
不適合拿來驗證開放性主張。

---

## 五、建議：先打哪一個

**第一順位：BLUEPRINT 的 5 名原發 MCL 對配對正常記憶 B 細胞的 H3K27ac／H3K4me1 peak。**

1. **它攻擊真正的弱點。** 其他每一項都在回答比較弱的問題（預測、結合、接觸、基因型）。
   只有這一項把「兩名病人」換成「五名病人，對上分化配對的正常記憶 B 細胞」。
2. **建構已經對了。** GRCh38，不需要 liftOver，省掉一個出錯來源。
3. **成本接近零。** peak 檔合計約 60 MB，公開 HTTP，不必申請，不碰 GPU。
4. **順序**：先用索引裡的 RNA-Seq gene quant 確認這 5 名病人表現 `SLC7A11`；
   再算全基因體 peak 覆蓋率當基準率；最後才做元件層比較。

**同時並行（兩件）**：

- 寄出 `EGAD00001002918`（MCL ATAC-seq、5 名病人、GRCh38）的 BLUEPRINT DAC 申請。
  那是上面那個檢定的「同檢測」版本，申請本身要數週。
- 用瀏覽器開 `cascade.finngen.fi` 查我們的區間，並讀 FinnGen multiome 論文
  （`10.1038/s41586-026-11078-2`）的 Data Availability。二十分鐘。
  期望值不高（見 2.2 節），但那是 n = 1,108 的那一份。

**第三順位（只在前面有正向結果時才值得）**：GSE97541 的三例原發 MCL H3K27ac 當第二個獨立隊列；
4DN 的 GM12878 Hi-C／PLAC-seq 看遠端元件是否接觸 `SLC7A11` 啟動子；
TenK10K 附的 `ChromBPNet/B_memory` 模型對六個元件打序列層的可及性分數。

**不要期待的事。** 今天沒有任何資料能回答「這五個元件是否調控 `SLC7A11`」。
那需要 B 系細胞裡的元件擾動實驗。ENCODE 有這個能力、也在 GM12878 與 BJAB 上做過
Flow-FISH CRISPR，但**從來沒有指向這個基因座**（`SLC7A11` 與 `MAFG` 各 0 個實驗）。

---

## 六、該寫進報告限制一節的四件事

1. **記憶 B 細胞 caQTL 已經到了，也已經查了，結果是沒有檢定力的空值。**
   TenK10K Phase 1（922 捐贈者、n = 913 記憶 B）在六個元件上的 caQTL qval 全部 ≥ 0.98；
   chr4 的顯著基準率是 2.01%，五個元件的期望命中數 0.100；
   作者自己的 caQTL→eQTL SMR（81,550 列）對 `SLC7A11` 與 `MAFG` 的命中數都是 0。
   **這與步驟 180 同一類：無檢定力，不是負結果。**
   README 目前寫「那些是正確的檢定」——檢測與細胞類型確實對了，但檢定力不在那裡。
2. **四個遠端元件（+99.6、+105.9、+150.5、+164.0 kb）沒有被任何增強子—基因資源提名**
   （ENCODE-rE2G 記憶 B 0 條、且連候選元件都不是；EpiMap 0 條；GeneHancer Double Elite 0 個），
   而這兩個資源對 `SLC7A11` 都保有檢定力（rE2G 找到 5 個其他元件、EpiMap 找到 70／37 條）。
   這與「它們在正常狀態下不是 `SLC7A11` 的增強子」相容，
   但**不構成對「它們在 MCL 中調控 `SLC7A11`」的任何支持**。
   唯一的正向外部線索是 GM12878 的 `BACH1` peak，位置在 +150.5 kb 或 +164.0 kb
   （ENCODE TF clusters 與 ReMap 不一致）。
3. **`MAFG` MARE 的外部支持是部分的。** 兩個 ENCODE-rE2G 記憶 B 捐贈者中，
   一個（ENCDO512FEI）給出含 motif 的元件 → `MAFG`、`class=genic`、`isSelfPromoter=FALSE`、
   score **0.657**；另一個（ENCDO352FHC）給出的是 self-promoter 連結，分數按定義給的、不可用。
   EpiMap 的對應連結與 motif 差 686 bp，不可用。
   ReMap 在該區找到 BACH1／BACH2／MAFF／MAFG／MAFK／NFE2L2，但那是 546 因子的 ChIP-seq 熱點，
   未經背景配對不能當特異性證據。
4. **「在正常 B 細胞中不開啟」在對外敘述上要寫成量的差異。**
   本專案自己的數字本來就是量的（`normal_log2cpm` 0.183–1.497、`ratio` 2.21–13.0）。
   外部資源對「開／關」的判定彼此相反：ENCODE-rE2G 的記憶 B 候選元件名單裡五個元件一個都沒進，
   TenK10K 卻把它們全部量化並檢定了（但其 peak 集合是跨細胞類型共享的共識集，
   相同邊界也出現在單核球，且 peak 比元件寬得多）。
   這個差異來自閾值與樣本數，不是來自事實，應該照實報告而不是選一邊。
