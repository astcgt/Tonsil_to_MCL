#!/usr/bin/env bash
# 下載 Tonsil Atlas + MCL 公開資料 (Zenodo)
set -u
OUT=/raid/jupyter-alz.ee09/Tonsil_to_MCL/data/raw_zenodo
cd "$OUT"
dl () { # $1=record  $2=filename
  if [ -s "$2" ]; then echo "[skip] $2 已存在"; return; fi
  echo "[get ] $2"
  curl -L --fail --retry 5 --retry-delay 10 -C - \
    -o "$2" "https://zenodo.org/records/$1/files/$2?download=1" \
    && echo "[done] $2 $(du -h "$2" | cut -f1)" || echo "[FAIL] $2"
}
dl 8373756  tonsil_atlas_donor_metadata.csv
dl 8373756  README.md
dl 10373041 MCL.zip
dl 8373756  TonsilAtlasSeuratRNA.tar.gz
dl 8373756  TonsilAtlasSeuratMultiome.tar.gz
echo "=== ALL DOWNLOADS FINISHED ==="
ls -lh "$OUT"
