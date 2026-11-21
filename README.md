# Tonsil → MCL

**Differentiation-matched regulatory rewiring in mantle cell lymphoma, and protein-level drug design against the hub.**

從扁桃體正常 B 細胞分化出發，找出被套細胞淋巴瘤打開的調控迴路，再對迴路樞紐設計蛋白質藥物。

<p align="center">
  <a href="#english">🇬🇧 <b>English</b></a> &nbsp;·&nbsp;
  <a href="#中文">🇹🇼 <b>中文</b></a> &nbsp;·&nbsp;
  <a href="https://huggingface.co/ast-atp/Tonsil_to_MCL">🤗 <b>Data &amp; structures</b></a> &nbsp;·&nbsp;
  <a href="#data-sources--資料來源">📚 <b>Data sources / 資料來源</b></a> &nbsp;·&nbsp;
  <a href="#ethics--data-use--倫理與資料使用">⚖️ <b>Ethics / 倫理</b></a>
</p>

---

<a name="english"></a>
## English

### The question

Almost every study asks how mantle cell lymphoma (MCL) differs from a normal B cell.
MCL cells look like early B cells, so that comparison mostly measures development, not cancer.

This project asks a different question: **how does an MCL cell differ from the normal B cell at the same stage of development?**

### The method

```
1. Build a reference of normal B-cell development
   46,216 tonsil B cells · 47 states · RNA and chromatin in the same nucleus

2. Give each tumour cell a soft position on those 47 states
   The model then predicts what a normal cell at that position would express

3. Subtract the prediction from the measurement
   What remains is the part that development cannot explain

4. Rebuild the control network from RNA and chromatin together
   Compare the tumour network with the normal network

5. Test with external models, then design a drug
```

### What we found

The strongest remaining signal is not a known cancer gene. It is a seven-letter piece of DNA: **TGACTCA** — the landing site for the AP-1 and small Maf protein families.

| Finding | Evidence |
|---|---|
| MCL opens TGACTCA sites genome-wide | AP-1 enriched 1.28–1.42×, BACH1 1.38× (background matched for length and GC) |
| B-cell identity factors lose their sites | POU family 0.70–0.85× |
| `SLC7A11` gains five elements absent from normal cells | Expected expression 0.001, measured 0.206 |
| The tumour needs five consecutive pathway genes | `SLC7A11` z = −2.89, `GCLM` −2.74, `FSP1` −2.60, `SLC3A2` −2.18, `GCLC` −1.65 |
| The pathway is not mutated | 1.10× background in non-coding regions; **0** coding mutations in 21 pathway genes, against **50** in known drivers (same pipeline) |
| The literature is almost empty | `"system xc−" + MCL` = **0** papers; `MAFG` 0; `SLC3A2` 0; `AIFM2/FSP1` 0; `SLC7A11` 1 |

The pathway runs in one line: cystine enters the cell → the cell builds glutathione → glutathione plus FSP1 blocks ferroptosis.
**The tumour does not keep this pathway by accident. It lives on it.**

### External validation

All of the above comes from the same two patients, which is a weakness. We ran three tests with AlphaGenome 2 (trained on normal tissue only). Two failed and one succeeded; all three are reported.

| Test | Result |
|---|---|
| Baseline accessibility of the new elements | ❌ **Failed and retracted.** A peak-width artefact. Width correlates with prediction at +0.74 and our elements are 3.5× wider. With a width-matched control the effect reverses (0.55×, p = 0.57) |
| Saturation mutagenesis of the MARE site | ✅ **MAFG promoter confirmed.** Breaking the site changes predicted `MAFG` expression 3.1× more than neighbouring bases (Bonferroni p ≈ 1×10⁻⁴). Positive control `HMOX1` works, and is selective — only 1 of 3 TGACTCA copies scores |
| Cell-state conditioning (our own cells) | ✅ **SLC7A11 elements confirmed, in two patients.** Same DNA, only the cell state changes, so width and composition cancel exactly. Elements open ~2× under the tumour state against ~1.2× for non-rewired peaks in the same window (p = 7.9×10⁻⁴ and 1.9×10⁻⁵) |

### Drug design

The AP-1 and BACH1 families are redundant, but they all need the same partner at the site: **MafG**.

Evolutionary constraint (a protein language model) and the crystal interface converge on the same residues — two completely different kinds of evidence (`MAFG` p = 8.3×10⁻¹¹). Those positions also leave little room to escape: the interface accepts 10.5% of substitutions against 42.4% for the whole protein.

The best design is a **90-residue mini-protein** that passes **all four wet-lab-calibrated filters**, thresholds derived from 3,766 binders tested in the laboratory:

| | ipSAE | ipTM | Interface error | pLDDT | Filters |
|---|---|---|---|---|---|
| **Best design** | **0.631** | 0.790 | 6.97 Å | 0.902 | **4 / 4** |
| Native NRF2 (positive control) | 0.539 | 0.713 | 9.16 Å | 0.937 | 3 / 4 |
| Scrambled (negative control) | 0.013–0.076 | 0.39–0.47 | 17–20 Å | 0.58–0.65 | 0 / 4 |

```
LAAGGSVDTEKLAELQKKLEEMELELAKMQQELAELQKELEKLQKELAEAEAEMKAAGIPAKLSEEERARLLHWAALVEAPLHRVRRRPL
```

### What we rejected

17 hypotheses were tested and dropped. Three were our own first conclusions. The full table is in the report and the showcase page.

Three separate artefacts came from the same cause: a foreground group and a background group that differ in length or base composition — motif length, then GC content, then peak width. Each one produced a beautiful and false result.

### Honest limits

- **Two patients.** No third patient-derived MCL multiome exists in public data.
- **`MAFG` elements gain coupling but lose openness.** We have no explanation. Our first explanation came from the width artefact and was retracted.
- **Delivery is out of scope by design.** A 90-residue protein does not cross a membrane by itself. Every design here assumes it reaches the nucleus.
- **The structure model cannot score single mutations on a protein–DNA interface.** 0 of 24 mutations passed the noise floor; the baseline is saturated.

---

<a name="中文"></a>
## 中文

### 問題

幾乎所有研究都在問「被套細胞淋巴瘤和正常 B 細胞差在哪」。但 MCL 細胞看起來就像早期 B 細胞，所以那種比較主要量到的是分化，不是癌症。

這個專案改問：**MCL 細胞和「分化程度完全相同的那顆正常 B 細胞」差在哪？**

### 方法

```
1. 建立正常 B 細胞分化參考
   46,216 顆扁桃體 B 細胞 · 47 個分化狀態 · RNA 與染色質同一顆細胞核

2. 讓每顆腫瘤細胞在這 47 個狀態上取得軟性位置
   再由模型推算「若這顆細胞是正常的，它該表現多少」

3. 把實測減掉推算值
   剩下的就是分化解釋不了的部分

4. 用 RNA 與染色質共同重建調控網絡
   比較腫瘤網絡與正常網絡

5. 用外部模型驗證，然後設計藥物
```

### 發現

最強的殘餘訊號不是任何已知癌症基因，而是一段 7 個字母的 DNA：**TGACTCA** —— AP-1 與 small Maf 家族的落腳點。

| 發現 | 證據 |
|---|---|
| MCL 全基因組打開 TGACTCA 位點 | AP-1 富集 1.28–1.42×、BACH1 1.38×（背景已配對長度與 GC） |
| B 細胞身分因子失去位點 | POU 家族 0.70–0.85× |
| `SLC7A11` 長出五個正常細胞沒有的元件 | 期望表現 0.001，實測 0.206 |
| 腫瘤需要連續五個路徑基因 | `SLC7A11` z = −2.89、`GCLM` −2.74、`FSP1` −2.60、`SLC3A2` −2.18、`GCLC` −1.65 |
| 這條路徑沒有突變驅動 | 非編碼區僅背景 1.10×；21 個路徑基因 **0** 個編碼突變，而同一流程在已知 driver 抓到 **50** 個 |
| 文獻幾乎是空的 | `"system xc−" + MCL` = **0** 篇；`MAFG` 0；`SLC3A2` 0；`AIFM2/FSP1` 0；`SLC7A11` 1 |

路徑是一條線：胱胺酸進入細胞 → 細胞合成穀胱甘肽 → 穀胱甘肽加上 FSP1 阻擋鐵死亡。
**腫瘤不是順便開著這條路，它靠這條路活著。**

### 外部驗證

以上全部來自同樣兩位病人，這是弱點。我們用 AlphaGenome 2（只學過正常組織）做了三個測試，兩個失敗、一個成功，三個都報告。

| 測試 | 結果 |
|---|---|
| 新生元件的基線可及性 | ❌ **失敗並撤回。** peak 寬度假象。寬度與預測值相關 +0.74，而我們的元件寬 3.5 倍。改用寬度配對對照後效應反轉（0.55×、p = 0.57） |
| MARE 位點的飽和突變 | ✅ **MAFG 啟動子獲確認。** 破壞該位點使預測的 `MAFG` 表現改變量為鄰近鹼基的 3.1 倍（Bonferroni p ≈ 1×10⁻⁴）。陽性對照 `HMOX1` 有效，且具選擇性 —— 三個 TGACTCA 只有一個顯著 |
| 細胞狀態條件化（用我們自己的細胞） | ✅ **SLC7A11 元件獲確認，兩位病人一致。** 同一段 DNA、只換細胞狀態，所以寬度與組成完全抵銷。元件在腫瘤狀態下開啟約 2 倍，同窗格未重接線的 peak 僅約 1.2 倍（p = 7.9×10⁻⁴ 與 1.9×10⁻⁵） |

### 藥物設計

AP-1 與 BACH1 家族成員彼此冗餘，但都需要同一個夥伴：**MafG**。

演化約束（蛋白質語言模型）與晶體介面收斂到同一組殘基 —— 兩種完全不同的證據（`MAFG` p = 8.3×10⁻¹¹）。這些位置也幾乎沒有逃逸空間：介面只接受 10.5% 的取代，全蛋白質平均是 42.4%。

最佳設計是一個 **90 殘基的迷你蛋白質**，通過**全部四個濕實驗校準門檻**（門檻取自 3,766 個實驗測試過的結合劑）：

| | ipSAE | ipTM | 介面誤差 | pLDDT | 通過 |
|---|---|---|---|---|---|
| **最佳設計** | **0.631** | 0.790 | 6.97 Å | 0.902 | **4 / 4** |
| 天然 NRF2（陽性對照） | 0.539 | 0.713 | 9.16 Å | 0.937 | 3 / 4 |
| 亂序（陰性對照） | 0.013–0.076 | 0.39–0.47 | 17–20 Å | 0.58–0.65 | 0 / 4 |

### 被推翻的假設

17 個假設被檢定後放棄，其中三個是我們自己最初的結論。完整表格在報告與展示頁。

三個獨立的假象來自同一個原因：前景組與背景組在長度或鹼基組成上不同 —— 先是 motif 長度，再是 GC 含量，最後是 peak 寬度。每一次都產生漂亮但錯誤的結果。

### 誠實的限制

- **只有兩位病人。** 公開資料中不存在第三例病人來源的 MCL 多體學。
- **`MAFG` 元件獲得協同性卻失去開放度。** 目前沒有解釋。我們最初的解釋來自寬度假象，已撤回。
- **藥物動力學依設計排除在範圍外。** 90 殘基的蛋白質無法自行穿過細胞膜，所有設計都預設可到達細胞核。
- **結構模型無法評分蛋白質–DNA 介面上的單點突變。** 24 個突變 0 個超過雜訊地板，基線已飽和。

---

<a name="data-sources--資料來源"></a>
## Data sources / 資料來源

All primary data is public. Nothing in this repository was generated by us at the bench.
所有原始資料皆為公開資料，本倉庫沒有任何我們自己實驗產生的資料。

### Primary single-cell data / 主要單細胞資料

| Dataset | Source | Link |
|---|---|---|
| **Human Tonsil Atlas** (normal B-cell reference, 46,216 cells, 47 states; and the MCL multiome samples M102 / M413) | Massoni-Badosa et al., *Immunity* 2024, "An atlas of cells in the human tonsil" | [DOI 10.1016/j.immuni.2024.11.006](https://doi.org/10.1016/j.immuni.2024.11.006) · Zenodo [8373756](https://zenodo.org/records/8373756) |
| **MCL external cohort** (11 patients, scRNA, used for replication) | GEO | [GSE303064](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE303064) |
| **MCL subclones** (5 patients, SingleCellExperiment with inferCNV labels) | Zenodo | [18200523](https://zenodo.org/records/18200523) |
| **MCL somatic mutations** (22 patients, whole genome) | Zenodo | [16098829](https://zenodo.org/records/16098829) |
| Additional GEO series used during exploration | GEO | [GSE184031](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184031) · [GSE239497](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE239497) |

### Reference genome and annotation / 參考基因組與註釋

| Resource | Version | Link |
|---|---|---|
| Human genome | GRCh38 primary assembly | [GENCODE release 44](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz) |
| Gene annotation | GENCODE v44 | [GTF](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.annotation.gtf.gz) |
| Transcription factor motifs | JASPAR2024 CORE vertebrates, non-redundant | [download](https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt) |
| Protein sequences and accessions | UniProt (REST) | [rest.uniprot.org](https://rest.uniprot.org/) |

### Structures / 結構

| PDB | Content | Used for |
|---|---|---|
| [7X5E](https://www.rcsb.org/structure/7X5E) | MafG with DNA | MafG dimer and DNA-binding interfaces; binder target definition |
| [8YOQ](https://www.rcsb.org/structure/8YOQ) | FSP1 (AIFM2), 1.79 Å, with 6-OH-FAD and NADP | FSP1 cofactor pocket; template for co-folding |
| [7OGS](https://www.rcsb.org/structure/7OGS) | IRF4 with DNA | Earlier target iteration (superseded) |

### Functional genomics and screens / 功能基因體與篩選

| Resource | Content | Link |
|---|---|---|
| **DepMap** | CRISPR gene-effect scores; MCL lines JVM-2, MAVER-1, Z-138 | [depmap.org](https://depmap.org/portal/) |
| **PIONEER2** | Interactome-wide interface predictions (9.3 GB bulk file; we stream-filtered it to 35 proteins) | [pioneer2.yulab.org/downloads](https://pioneer2.yulab.org/downloads) · preprint [bioRxiv 2025.06.09.658393](https://doi.org/10.1101/2025.06.09.658393) |
| **PubMed** | Novelty screening via E-utilities | [eutils](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/) |

### Models / 模型

| Model | Role | Link |
|---|---|---|
| **AlphaGenome 2** (Google DeepMind) | DNA sequence → chromatin and expression; saturation mutagenesis; cell-state conditioning | [alphagenome.google/api](https://alphagenome.google/api) |
| **ESM-2 650M** (Meta) | Protein language model; saturation mutagenesis of 7 proteins, 60,000+ variants | [github.com/facebookresearch/esm](https://github.com/facebookresearch/esm) |
| **Boltz-2** (MIT) | Co-folding; 110 complexes / 159 predictions. Open equivalent of AlphaFold 3 | [github.com/jwohlwend/boltz](https://github.com/jwohlwend/boltz) |
| **RFdiffusion3** (Baker lab) | Generative backbone design | [github.com/RosettaCommons](https://github.com/RosettaCommons) |
| **ProteinMPNN** (Baker lab) | Sequence design on a backbone | [github.com/dauparas/ProteinMPNN](https://github.com/dauparas/ProteinMPNN) |
| **Borzoi** (Calico) | Sequence → RNA coverage; in-silico deletion. Successor to Enformer | [github.com/calico/borzoi](https://github.com/calico/borzoi) · via [gReLU](https://genentech.github.io/gReLU/) |
| **IPSAE** (Dunbrack lab) | Interface scoring. We reimplemented the formulas and reconciled exactly against the official tool | [github.com/DunbrackLab/IPSAE](https://github.com/DunbrackLab/IPSAE) · [bioRxiv 2025.02.10.637595](https://doi.org/10.1101/2025.02.10.637595) |

**AlphaFold 3 was not used.** Its weights require an application to DeepMind. Boltz-2 is the open model of the same class and is what we ran.
**沒有使用 AlphaFold 3** —— 其權重需向 DeepMind 申請。Boltz-2 是同級的開源模型，我們用的是它。

### Calibration thresholds / 校準門檻

The binder filters are not ours. They come from published benchmarks:

| Threshold | Source |
|---|---|
| `ipSAE > 0.61` | Overath et al., 3,766 experimentally tested binders / 15 targets, [bioRxiv 2025.08.14.670059](https://doi.org/10.1101/2025.08.14.670059) |
| `ipTM ≥ 0.50`, `interface pAE ≤ 10.85 Å`, `pLDDT ≥ 0.80` | BindCraft default filters |
| `pae_interaction < 10 Å` | Bennett et al., *Nature Communications* 2023 — an independent laboratory converging on the same number |

---

<a name="ethics--data-use--倫理與資料使用"></a>
## Ethics and data use / 倫理與資料使用

This repository contains **code and aggregate results only**.

**Deliberately excluded from this repository:**

| Excluded | Why |
|---|---|
| Per-cell tables with cell barcodes and donor identifiers | Patient-derived. Redistribution is outside the scope of the original data use terms |
| Per-patient somatic variant calls (`chrom`/`pos`/`ref`/`alt`) | Individual-level genotype data |
| Raw sequencing data, count matrices, genome FASTA | Size, and the same considerations above |
| The 9.3 GB PIONEER2 bulk file and the 38 MB filtered extract | Redistribution of another group's dataset |
| Model weights | Size, and each model's own licence |
| API keys | Never committed. Stored outside the repository with `600` permissions |

**What is included:** every analysis script, the reports, the interactive showcase page, and results at gene level, peak level, design level or as statistical summaries.

To reproduce the analysis, download the primary data from the links above and run the scripts in numerical order. Every script states its inputs and outputs in its docstring.

**倉庫只包含程式碼與彙總層級結果。** 刻意排除：含細胞條碼與捐贈者編號的逐細胞表格、每位病人的體細胞變異位點、原始定序資料、PIONEER2 的原始與過濾檔、模型權重、API 金鑰。要重現分析，請自上方連結下載原始資料，再依編號順序執行腳本。

### Citation

If this work is useful, please cite the primary data sources above — above all the Human Tonsil Atlas (Massoni-Badosa et al., *Immunity* 2024), without which none of this exists.

---

## Repository layout / 倉庫結構

```
scripts/     73 analysis scripts, numbered in execution order
             73 個分析腳本，依執行順序編號
results/     Aggregate results: CSV and JSON (57 MB)
             彙總結果：CSV 與 JSON
docs/        Tool status notes, API application steps
             工具現況備忘、API 申請步驟
licenses/    Licenses of third-party files redistributed here
             本倉庫轉載的第三方檔案的授權副本
REPORT_v3.md Full report with every rejected hypothesis (Chinese)
             完整報告，含所有被推翻的假設（中文）
af3_input.md THIRD-PARTY: verbatim copy of AlphaFold 3 docs/input.md
             (Google DeepMind, Apache-2.0) — see the header in that file
             第三方：AlphaFold 3 docs/input.md 的逐字複本，出處見該檔檔首
```

**Large files / 大型檔案**：the tables and structure predictions that exceed this repository's
size budget live at <https://huggingface.co/ast-atp/Tonsil_to_MCL> — peak-level and pseudobulk
tables, the prepared target surfaces, and every co-folding and design output, which is enough to
recompute each interface metric reported here without re-running the pipeline.
超出本倉庫體積預算的彙總表與結構預測放在上面那個 Hugging Face 倉庫，可用來重算報告裡的每個界面指標。

**Not redistributed / 刻意未轉載**：`results/95_pioneer2_axis.tsv.gz` is a derived subset of the
PIONEER2 interface dataset. That dataset states **no licence** on any of its pages (checked
2026-11-22, see `docs/TOOLING_NOTES.md`), so it defaults to all rights reserved and we do not
republish it here or on Hugging Face. Run `scripts/95_pioneer2.py` to rebuild it from the original
source. PIONEER2 的介面資料集**沒有授權聲明**，所以我們不轉載衍生檔；用該腳本可自行重建。

**License / 授權**：CC BY 4.0 (`LICENSE`). Use freely, including commercially, with attribution.
可自由使用（含商業用途），但須標注來源。

**Third-party content / 第三方內容**：`af3_input.md` is not our work. It is a
verbatim copy of `docs/input.md` from
[google-deepmind/alphafold3](https://github.com/google-deepmind/alphafold3),
© Google DeepMind, redistributed under the Apache License 2.0 (copy in
`licenses/`). Everything else in this repository is our own.
`af3_input.md` 不是本專案的產出，是上述檔案的逐字複本，依 Apache-2.0 轉載。

### Key scripts / 關鍵腳本

| Script | What it does |
|---|---|
| `13_assemble_reference.py` | Build the normal B-cell multiome reference |
| `15_residual_core.py` | Differentiation matching and residual calculation |
| `50_mcl_atac_on_atlas.py` | Requantify MCL ATAC onto the atlas peak coordinates |
| `52_peak_gene_links.py` | Peak–gene links with an empirical null |
| `27_motif_enrichment.py` | Motif enrichment with a length- and GC-matched background |
| `71_esm_dms.py` | ESM-2 saturation mutagenesis |
| `91_alphagenome.py` | AlphaGenome 2 baseline and saturation mutagenesis |
| `92_ag_cellcond.py` | AlphaGenome 2 cell-state conditioning |
| `93_ipsae.py` | ipSAE, pDockQ, pDockQ2, LIS, and the calibrated filters |
| `95_pioneer2.py` | Stream-filter the 9.3 GB PIONEER2 file without storing it |

---

## Status / 狀態

This is **exploratory computational work**. Nothing here has been tested in a laboratory.
Every claim is a hypothesis with a stated confidence level and, where we found one, a stated failure mode.

這是**探索性的計算工作**，沒有任何結果經過實驗驗證。每一項主張都是帶有明確信心水準的假說，並在我們找到時標明其失效模式。
