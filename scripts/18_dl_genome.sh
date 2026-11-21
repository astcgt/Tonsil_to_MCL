#!/usr/bin/env bash
# 下載 hg38 基因體與註解（motif 掃描、AlphaGenome/Borzoi 序列輸入都需要）
set -u
OUT=/raid/jupyter-alz.ee09/Tonsil_to_MCL/data/genome
mkdir -p "$OUT"; cd "$OUT"

get(){ [ -s "$2" ] && { echo "[skip] $2"; return; }; echo "[get ] $2";
  curl -sL --fail --retry 4 -o "$2" "$1" && echo "[done] $2 $(du -h "$2"|cut -f1)" || echo "[FAIL] $2"; }

# hg38 (與 CellRanger refdata GRCh38-2020-A 同一組 primary assembly)
get "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz" \
    "GRCh38.primary_assembly.genome.fa.gz"
get "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.annotation.gtf.gz" \
    "gencode.v44.annotation.gtf.gz"

# JASPAR 2024 CORE vertebrates (非冗餘) PFM
get "https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt" \
    "JASPAR2024_CORE_vertebrates_nr_pfms.txt"

echo "=== 解壓基因體 ==="
if [ -s GRCh38.primary_assembly.genome.fa.gz ] && [ ! -s GRCh38.primary_assembly.genome.fa ]; then
  gunzip -k GRCh38.primary_assembly.genome.fa.gz && echo "[gunzip] fa ok"
fi
ls -lh "$OUT"
echo "=== GENOME DONE ==="
