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
