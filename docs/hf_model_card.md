---
license: cc-by-4.0
language:
- en
- zh
tags:
- biology
- genomics
- single-cell
- chromatin
- protein-design
- lymphoma
pretty_name: "Tonsil to MCL — rewiring tables, structure predictions and target surfaces"
---

# Tonsil → MCL: data and structure predictions

Companion artifact repository for a study of **differentiation-matched regulatory rewiring in
mantle cell lymphoma (MCL)** and protein-level drug design against the circuit it opens.

**Code, methods and the full report live in the code repository:
<https://github.com/astcgt/Tonsil_to_MCL>**
This repository holds the files that are too large or too numerous for that repository, so that the
numbers in the report can be checked without re-running the pipeline.

---

## What the study asks

Almost every MCL study compares tumour cells with a normal B cell. MCL cells resemble early B cells,
so that comparison mostly measures development rather than cancer. This study compares each tumour
cell with **the normal tonsil B cell at the same stage of development**, and asks what remains.

What remains is not a known cancer gene but a seven-letter DNA site, `TGACTCA`, the landing site of
the AP-1 and small Maf transcription factor families. It leads to the cystine, glutathione and
ferroptosis-defence pathway through `SLC7A11`.

## What is in this repository

| Folder | Contents | Size |
|---|---|---|
| `tables/` | Peak-level and pseudobulk tables that exceed the code repository's size budget | 110 MB |
| `structures/` | Target surfaces used for binder design, as prepared (renumbered, single chain where needed) | 0.6 MB |
| `predictions/` | Raw co-folding and design outputs: 15 pipeline steps, 2,061 files (`.cif`, `.npz`, `.json`, `.yaml`) | 149 MB |

### `tables/`

| File | Rows are | Columns of interest |
|---|---|---|
| `22_atac_diff.csv` | one per chromatin peak (215,204) | mean accessibility in tumour and normal, log fold change per donor, consistency flag, coordinates |
| `26_peak_annotated.csv` | the same peaks | plus nearest gene and distance to its transcription start site |
| `22_atac_pseudobulk.csv.gz` | one per peak | accessibility per normal differentiation state and per tumour sample |
| `15_normal_state_pseudobulk.csv` | one per normal differentiation state (46) | expression of every gene, used as the expectation that tumour cells are compared against |

All four are aggregate: per-peak or per-state summaries. Per-cell tables, per-patient variant calls
and raw data are deliberately not published in either repository.

One further file is held back for a different reason. The study also used the PIONEER2 protein
interface dataset, and a filtered subset of it sits in the pipeline's outputs. That dataset states
no licence on any page of its web server (checked 2026-11-22), so it defaults to all rights
reserved and we do not republish it. The code repository's `scripts/95_pioneer2.py` rebuilds the
subset from the original source, and `docs/TOOLING_NOTES.md` records what was checked.

### `structures/`

`mafg_target.pdb` and `fsp1_target.pdb` are the prepared target surfaces, derived from PDB entries
`7X5E` (MafG) and `8YOQ` (FSP1). `fsp1_template_renumbered.pdb` is the same FSP1 structure with
residues renumbered to start at 1 and `SEQRES` records added, which co-folding with a template
requires. `irf4_target.pdb` belongs to an earlier line of work that the study later set aside.

### `predictions/`

The outputs of Boltz-2 co-folding and of RFdiffusion3 / ProteinMPNN design, one folder per pipeline
step, named by the step number used in the code repository. Each prediction folder carries the
structure (`.cif`), the predicted error matrices (`pae_*.npz`, `pde_*.npz`) and the confidence
summary (`confidence_*.json`). These are the inputs to every interface metric reported in the study,
so a reader can recompute ipSAE, ipTM, pDockQ, pDockQ2 and LIS from them directly.

## Headline results these files support

- MCL opens `TGACTCA` sites genome-wide: AP-1 enriched 1.28–1.42×, BACH1 1.38×, against a background
  matched for **both length and GC content**. Without that matching, 503 motifs look significant; with
  it, 25 do.
- `SLC7A11` gains five regulatory elements absent from normal B cells. Expected expression at the
  matched normal state is 0.001; measured in tumour cells it is 0.206.
- Five consecutive genes of the pathway are MCL-selective dependencies in DepMap
  (`SLC7A11` z = −2.89, `GCLM` −2.74, `AIFM2` −2.60, `SLC3A2` −2.18, `GCLC` −1.65).
- The pathway carries no driver mutation: 1.10× background in non-coding regions, and 0 coding
  mutations across 21 pathway genes against 50 in known drivers from the same pipeline.
- The five elements are already wired and already unlocked, yet closed: they sit in
  contact-privileged positions relative to the promoter across 21 measured contact maps, they are
  relatively unmethylated in normal B cells (0.21 against 0.67 for neighbouring peaks and 0.90 for
  genomic background), and their accessibility is as low as peaks that do nothing. The missing part
  is the transcription factor, not chromatin licensing.
- A 90-residue mini-protein against MafG passes all four wet-lab-calibrated filters and exceeds the
  native NRF2 partner on three of four metrics.
- The MafG–NRF2 interface carries a defined epitope: a binding-energy model, calibrated against
  measured energies on a homology-clean split, marks 14 hotspots on MafG spanning residues 69–107 and
  11 on NRF2 spanning 514–553, including every leucine of the zipper register. The same run
  overturned our own earlier per-mutation ranking, which came from a co-folding confidence score.

## Limitations

- **Two patients.** No third patient-derived MCL multiome exists in public data.
- **Predictions, not measurements.** Every structure here is a model. No binding affinity,
  no crystal structure, no cell assay.
- **One arm is weaker than the other.** The `SLC7A11` element claim has clean internal evidence and
  two independent external supports. The claim that MafG is the druggable hub is weaker: its contact
  signal was retracted after peak accessibility explained it, and its gained coupling failed two
  chromatin tests. Of the four interface residues we first called critical, only **L79** is ranked
  critical by two methods.
- The study's own report lists 75 hypotheses that were tested and rejected, including several of the
  authors' own errors. Read that section before building on any single number here.

## Provenance of the primary data

Derived from public data: the human tonsil atlas (Massoni-Badosa et al., *Immunity* 2024) and its
Zenodo deposits, MCL samples from the same atlas project, DepMap, JASPAR2024, GENCODE v44 and the
PDB. The code repository's README lists every accession and URL.

## Citation

If you use these files, cite the code repository and the primary data sources it lists:

```
Tonsil_to_MCL: differentiation-matched regulatory rewiring in mantle cell lymphoma,
with protein-level drug design. https://github.com/astcgt/Tonsil_to_MCL
```

Licensed **CC BY 4.0**: use freely, including commercially, with attribution.

---

## 中文說明：扁桃體 → MCL 的資料與結構預測

這是一項研究的**成果檔案倉庫**，研究主題為**套細胞淋巴瘤（MCL）中經分化配對後的調控重接線**，
以及針對它所開啟的那條迴路進行的蛋白質層級藥物設計。

**程式碼、方法與完整報告均置於程式碼倉庫：
<https://github.com/astcgt/Tonsil_to_MCL>**
本倉庫存放的是該處容量無法容納、或檔案數量過多的檔案，目的在於使報告中的數字**無須重跑整條流程即可查核**。

---

### 本研究的問題

幾乎所有 MCL 研究都是將腫瘤細胞與正常 B 細胞相比。但 MCL 細胞的表型近似早期 B 細胞，因此這類比較
所測得的大多是發育差異，而非癌化差異。本研究改為將每個腫瘤細胞與**分化程度完全相同的正常扁桃體
B 細胞**相比，再追問其殘餘為何。

殘餘下來的並非任何已知癌症基因，而是一段七個字母的 DNA 位點 `TGACTCA`，即 AP-1 與 small Maf
轉錄因子家族的結合位點，並經由 `SLC7A11` 通往胱胺酸、穀胱甘肽與鐵死亡防禦路徑。

### 本倉庫的內容

| 資料夾 | 內容 | 大小 |
|---|---|---|
| `tables/` | 超出程式碼倉庫體積上限的 peak 層級與 pseudobulk 彙總表 | 110 MB |
| `structures/` | 結合劑設計所用的標的表面，為已處理版本（已重新編號，必要時取單鏈） | 0.6 MB |
| `predictions/` | 共摺疊與設計的原始輸出：15 個流程步驟、2,061 個檔案（`.cif`、`.npz`、`.json`、`.yaml`） | 149 MB |

#### `tables/`

| 檔案 | 每一列代表 | 主要欄位 |
|---|---|---|
| `22_atac_diff.csv` | 一個染色質 peak（共 215,204 個） | 腫瘤與正常的平均可及性、各捐贈者的 log fold change、一致性旗標、座標 |
| `26_peak_annotated.csv` | 同一批 peak | 另附最近的基因及其轉錄起始位點的距離 |
| `22_atac_pseudobulk.csv.gz` | 一個 peak | 各正常分化狀態與各腫瘤樣本的可及性 |
| `15_normal_state_pseudobulk.csv` | 一個正常分化狀態（共 46 個） | 每個基因的表現量，作為腫瘤細胞比較時的期望值 |

這四份均屬彙總層級，即以 peak 或以狀態為單位的摘要。逐細胞表格、各病人的變異位點與原始資料，
兩個倉庫均刻意不予公開。

另有一份檔案因其他原因未予公開。本研究亦使用了 **PIONEER2** 蛋白質介面資料集，其篩選後的子集存在
於流程輸出之中。該資料集的網站任何一頁均未標示授權（查核日期 2026-11-22），依預設即為保留一切
權利，因此我們不予重新散布。程式碼倉庫的 `scripts/95_pioneer2.py` 可自原始來源重建該子集，
`docs/TOOLING_NOTES.md` 則記錄了查核的過程。

#### `structures/`

`mafg_target.pdb` 與 `fsp1_target.pdb` 為已處理完成的標的表面，分別取自 PDB 條目 `7X5E`（MafG）
與 `8YOQ`（FSP1）。`fsp1_template_renumbered.pdb` 為同一個 FSP1 結構，其殘基重新編號為自 1 起算並
補上 `SEQRES` 記錄，此為帶模板共摺疊的必要條件。`irf4_target.pdb` 屬於本研究後來擱置的一條早期路線。

#### `predictions/`

Boltz-2 共摺疊與 RFdiffusion3／ProteinMPNN 設計的輸出，一個流程步驟對應一個資料夾，並以程式碼倉庫
所使用的步驟編號命名。每個預測資料夾均含結構（`.cif`）、預測誤差矩陣（`pae_*.npz`、`pde_*.npz`）
與信賴度摘要（`confidence_*.json`）。這些是本研究所有介面指標的輸入，因此讀者可直接由其重新計算
ipSAE、ipTM、pDockQ、pDockQ2 與 LIS。

### 這些檔案所支持的主要結果

- MCL 在全基因組範圍開放 `TGACTCA` 位點：AP-1 富集 1.28–1.42×、BACH1 1.38×，背景**同時就長度與
  GC 含量配對**。未做此層配對時有 503 個 motif 達到顯著；配對之後僅餘 25 個。
- `SLC7A11` 新生出五個正常 B 細胞所無的調控元件。在分化配對的正常狀態下，其期望表現量為 0.001，
  而在腫瘤細胞中實測為 0.206。
- 這條路徑連續五個基因在 DepMap 中均為 MCL 選擇性依賴
  （`SLC7A11` z = −2.89、`GCLM` −2.74、`AIFM2` −2.60、`SLC3A2` −2.18、`GCLC` −1.65）。
- 這條路徑不帶驅動突變：非編碼區僅為背景的 1.10×；21 個路徑基因共 0 個編碼突變，而同一流程在已知
  driver 上檢出 50 個。
- 那五個元件線路已接通、也已解鎖，卻仍處於關閉狀態：在 21 張實測接觸圖中，它們相對於啟動子均位於
  接觸上占優勢的位置；在正常 B 細胞中甲基化程度偏低（0.21，鄰近 peak 為 0.67、基因組背景為 0.90）；
  而其可及性低至與無任何活性的 peak 相當。所缺的是轉錄因子，而非染色質層面的許可。
- 一個針對 MafG 的 90 殘基迷你蛋白質通過全部四項濕實驗校準門檻，並在四個指標中有三個優於天然的
  NRF2 夥伴。
- MafG–NRF2 介面帶有一個明確的抗原決定位：一個以實測結合能、並在無同源性洩漏的切分上校準過的模型，
  在 MafG 上標出 14 個熱點（殘基 69–107），在 NRF2 上標出 11 個（514–553），其中包含白胺酸拉鍊
  七肽重複位相（register）上的每一個白胺酸。同一次計算也推翻了我們先前的單點突變排名，該排名來自
  共摺疊的信賴度分數。

### 限制

- **僅有兩位病人。** 公開資料中不存在第三例病人來源的 MCL 多體學資料。
- **此處為預測，而非量測。** 本倉庫每一個結構都是模型，不含任何結合親和力、晶體結構或細胞實驗。
- **兩條線的強度並不相同。** `SLC7A11` 元件的主張具備明確的內部證據與兩項獨立的外部支持。
  「MafG 為可成藥樞紐」這一主張則較弱：其接觸訊號在被 peak 可及性解釋之後已撤回，
  其獲得的協同性在兩項染色質檢定中均未通過。我們最初認定為關鍵的四個介面殘基中，
  僅 **L79** 被兩種方法同時列為關鍵。
- 本研究的報告列出 75 個經檢定後推翻的假設，其中包含作者自身的錯誤。在引用本倉庫任何單一數字之前，
  建議先閱讀該節。

### 原始資料來源

全部衍生自公開資料：人類扁桃體圖譜（Massoni-Badosa 等人，*Immunity* 2024）及其 Zenodo 存放檔、
同一圖譜專案的 MCL 樣本、DepMap、JASPAR2024、GENCODE v44 與 PDB。
每一筆取用編號與網址均列於程式碼倉庫的 README。

### 引用

若使用這些檔案，請引用程式碼倉庫以及其所列出的原始資料來源：

```
Tonsil_to_MCL: differentiation-matched regulatory rewiring in mantle cell lymphoma,
with protein-level drug design. https://github.com/astcgt/Tonsil_to_MCL
```

授權為 **CC BY 4.0**：可自由使用，包含商業用途，但須標注來源。
