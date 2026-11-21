#!/usr/bin/env bash
# 下載外部公開 MCL 單細胞世代（擴充樣本數用）
set -u
OUT=/raid/jupyter-alz.ee09/Tonsil_to_MCL/data/external_mcl
mkdir -p "$OUT"/{GSE303064,GSE239497,GSE184031,zenodo18200523,wgs}
cd "$OUT"

get(){ # $1=url $2=dest
  if [ -s "$2" ]; then echo "[skip] $2"; return; fi
  echo "[get ] $2"
  curl -sL --fail --retry 4 --retry-delay 8 -o "$2" "$1" \
    && echo "[done] $2 $(du -h "$2"|cut -f1)" || echo "[FAIL] $2"
}

# ① GSE303064 — 11 病人 / 20 樣本 / ~99k 細胞（含非惡性 B 與完整微環境）
get "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE303nnn/GSE303064/suppl/GSE303064_RAW.tar" \
    "GSE303064/GSE303064_RAW.tar"

# ② GSE239497 — 6 病人 / 7 樣本，CD19+ 磁珠富集（與 tonsil MCL CD19+ 最可比）
get "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE239nnn/GSE239497/suppl/GSE239497_RAW.tar" \
    "GSE239497/GSE239497_RAW.tar"

# ③ GSE184031 — 4 病人，輕量
get "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184031/suppl/GSE184031_RAW.tar" \
    "GSE184031/GSE184031_RAW.tar"

# ④ Zenodo 18200523 — 5 病人，SingleCellExperiment（已含 cell_type 與 inferCNV subclone）
get "https://zenodo.org/api/records/18200523/files/Supplemental%20Data%20Object%201.rds/content" \
    "zenodo18200523/MCL_subclones_SCE.rds"

# ⑤ 配套 WGS 體細胞突變（22 病人 / 35 樣本）
get "https://zenodo.org/api/records/16098829/files/MCL.somaticMutations.tsv/content" \
    "wgs/MCL.somaticMutations.tsv"

echo "=== 解開 tar ==="
for g in GSE303064 GSE239497 GSE184031; do
  if [ -s "$g/${g}_RAW.tar" ]; then
    ( cd "$g" && tar xf "${g}_RAW.tar" && echo "[untar] $g -> $(ls | wc -l) files" )
  fi
done
echo "=== EXTERNAL MCL DONE ==="
du -sh "$OUT"/*
