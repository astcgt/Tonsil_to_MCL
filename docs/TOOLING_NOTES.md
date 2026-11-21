# 工具現況備忘（2026-11-21 查證）

## ⚠️ 需要使用者行動（有時效）

### 1. AlphaGenome 2 API key —— 免費，但要你的 Google 帳號
- 申請頁：<https://alphagenome.google/api>（會轉到 `deepmind.google.com/science/alphagenome/api`）
- 非商業研究**免費**。官方文件沒有提到審核等待期。
- 拿到後設環境變數：`export ALPHA_GENOME_API_KEY=<你的key>`
- 條款：<https://alphagenome.google/terms>（研究用途、輸出不得拿去訓練其他模型）

### 2. Use-case 表單 —— **截止 2026-11-30**
- <https://forms.gle/ftBfHGH7NUrCZEjN6>
- 這**不是**拿 key 的必要步驟，是「爭取額外配額與 credits」。
- 我們這個題目（MCL t(11;14) enhancer hijacking 對 B 細胞 CCND1 的影響）很適合送。

---

## AlphaGenome 2 規格重點
| 項目 | 內容 |
|---|---|
| 套件 | `pip install alphagenome`（2.0.0），**需 Python ≥3.11**（本專案 env 是 3.11 ✅） |
| 子模組 | `from alphagenome.models.v2 import dna_client` |
| 輸入長度 | 16 KB / 128 KB / 512 KB / **1 MB**（與 v1 相同，沒變大） |
| 條件化 | v2 **一次只能指定一個 condition**（composite ontology CURIE），要比多組織得迴圈 |
| 關鍵能力 | `SC_GENE_EXPRESSION`(1,943 biosamples)、`SC_ATAC`(704)、`RNA_SEQ`、`CHIP_HISTONE`、`CONTACT_MAPS`、`CELL_EMBEDDINGS` |
| **結構變異** | **`predict_sequence()` 接受任意 ACGTN 序列** → 可自組 der(14) 嵌合序列。官方 FAQ 有 `arbitrary-sequences` 專章，但明言**大型 SV 未經 benchmark** → 當假說生成而非定量預測 |
| 變異分數 | v2 **只有 `raw_score`，沒有 `quantile_score`** |
| 硬限制 | 每 request ≤23 scorers；expression conditioning ≤256 cells；ISM 區間 ≤10 bp；gRPC timeout 60 s；`RESOURCE_EXHAUSTED` 可重試（代表有配額節流，數字未公布） |
| Atlas | 90 億 SNV 預計算分數，`from alphagenome.atlas import atlas`。**只能逐一查詢（client 端 10 threads），不能整批下載；只有 SNV，沒有 SV** |
| 官方範例 | `colabs/v2/tal1_oncogenic_mutations.ipynb`（致癌非編碼突變造出新 enhancer）、`enhancer_to_gene_bcl11a.ipynb`、`haplotypes.ipynb` |

## Plan B：Borzoi via gReLU（現在就能跑，不用等 key）
- `pip install gReLU`（MIT，Python ≥3.8）→ `grelu.resources.load_model("Genentech/borzoi-model", "human_rep0.ckpt")`
- **必須用 `predict_on_seqs()`**，不能用 `predict_variant_effects()`（官方明文不支援 indel，更不用說重排）
- 輸入 524,288 bp，但**實際輸出只有中央 196,608 bp** → enhancer–gene 最大可模擬距離 ≈ 360 kb
- B 細胞相關 track：`RNA:GM12878`(28)、`ATAC:Memory B`、`ATAC:Plasma B`、`DNASE:GM12878`；**tonsil 只有 2 條 CAGE，沒有 RNA/ATAC**
- 兩個一定要做對：(1) REF/ALT 讓目標基因落在**相同 offset**；(2) 用 `undo_squashed_scale()` 反轉尺度，否則 fold change 幅度是錯的（方向仍對）
- 腳本骨架已存於 `scripts/mcl_t11_14_borzoi.py`（**未實際執行過**）

### 其他開源選項
- **scooby**（MIT，gagneurlab）：Borzoi-based，**同時給 scRNA + scATAC 單細胞解析度**，最接近 AG2 的開源方案。但 `-scooby-2` 權重繼承 AlphaGenome 非商業條款。
- **Enformer 不適用**：5,313 條 track 裡 **RNA-seq = 0 條**，且 naive/memory B ATAC = 0 條，視窗也太小。
- **Evo2 不適用**：只輸出 likelihood/embedding，沒有細胞型別特異的表現/可及性輸出頭。

## 蛋白質端
| 工具 | 授權 | 取得 | A100 40GB |
|---|---|---|---|
| AlphaFold3 | 程式碼 Apache-2.0，**權重禁止一切商業用途（含商業研究）、禁止分享、禁止拿 output 訓練模型** | 已不需申請表：`https://storage.googleapis.com/alphafold3/af3.bin.zst` | 可，但需調設定（官方只正式支援 80GB） |
| **Boltz-2** | **MIT（學術+商業皆可）** | `pip install boltz[cuda]` | 可。**只支援單一小分子親和力，不支援 protein-protein affinity** |
| OpenFold3 | Apache-2.0 | repo `aqlaboratory/openfold-3`，需 Python ≥3.12 | 可 |
| Chai-1 | Apache-2.0（含權重） | `pip install chai_lab` | 中小型 complex 可 |
| Chai-2 | 封閉申請制 | 拿不到 | — |

### ESM —— 已移交 Chan Zuckerberg Biohub 並改為 MIT
`facebook/*`、`EvolutionaryScale/*` 全部 redirect 到 `biohub/*`，授權 MIT、非 gated。
**點突變 fitness/LLR 建議用 `facebook/esm2_t33_650M_UR50D`**（bf16 僅 ~1.3 GB，吞吐遠勝 3B/15B，瓶頸是位點數×序列長度而非模型大小）。

### Binder design
- `RFdiffusion → ProteinMPNN → nrbennet/dl_binder_design`（注意是單 t 的 `nrbennet`）。RFdiffusion 的 conda env 是 CUDA 11.1，**建議直接用 Docker image**。
- 新平台 **RosettaCommons Foundry**：`pip install "rc-foundry[all]"`，含 RFdiffusion3/ProteinMPNN/LigandMPNN，BSD-3 含權重，安裝體驗好很多。
- **BindCraft 時間是最大瓶頸**：第三方估計單張 A100 上每個通過 filter 的設計約 34 小時，官方建議產出 ≥100 個 → 單卡數週等級。**本專案不採用全量 BindCraft**，改走 RFdiffusion 小批量 + 嚴格 in-silico 篩選。

## 第三方資料的授權狀態（決定能不能轉載）

| 來源 | 授權 | 我們怎麼處理 |
|---|---|---|
| AlphaFold 3 文件 `docs/input.md` | **Apache-2.0**（第 1 條把 documentation 算在 Work 內；第 4 條要求附授權副本與保留出處） | 逐字轉載於 `af3_input.md`，檔首加出處標註，授權全文放 `licenses/` |
| **PIONEER2 介面資料集** | **沒有授權聲明**（2026-11-22 查證） | **不轉載。** 只公開腳本 `scripts/95_pioneer2.py`，由使用者自行從原始 URL 重建 |

### PIONEER2 的查證紀錄（2026-11-22）

檔案來源：`https://pioneer2.yulab.org/static/PIONEER2_interactome_wide_interfaces.txt`
（由 `scripts/95_pioneer2.py` 串流下載並過濾，衍生子集 `results/95_pioneer2_axis.tsv.gz`，37 MB）

讀過的頁面與結果：

| 頁面 | 狀態 | 有授權嗎 |
|---|---|---|
| `pioneer2.yulab.org/` | 200 | 無。導覽列只有 Home / Download / About / Help，沒有 Terms 或 License 連結 |
| `/downloads` | 200 | 無。只描述檔案內容（352,124 個實驗偵測到的二元交互作用），未附任何條款 |
| `/about` | 200 | 無。只有「Version 1, May 2026」、聯絡信箱 `weill-pioneer@cornell.edu`、頁尾 `© Yu Lab` |
| `/help` | 200 | 無 |
| `/terms`、`/download` | 404 | 不存在 |
| 資料檔本身 | — | 衍生檔表頭只有欄位名，**沒有嵌入的授權或引用要求**；原始檔已不在本機 |

也找不到論文、Zenodo 或 GitHub 等有條款的其他發布點，網站亦未要求特定引用格式。

**結論**：沒有授權聲明 → 預設著作權（保留所有權利）→ 在 CC BY 4.0 的公開倉庫裡轉載衍生子集
需要權利人許可。因此**不放進任何公開倉庫**。
要取得該檔，執行 `scripts/95_pioneer2.py`，它會從上述原始 URL 重新下載並重建同一個子集。

若將來要公開那個檔案，正確做法是先寫信給 `weill-pioneer@cornell.edu` 取得書面許可或他們的授權條款。
