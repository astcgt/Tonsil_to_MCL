#!/usr/bin/env bash
set -u
OUT=/raid/jupyter-alz.ee09/Tonsil_to_MCL/data/depmap
B=https://ndownloader.figshare.com/files
get(){ [ -s "$OUT/$2" ] && { echo "[skip] $2"; return; }; curl -sL --retry 4 -o "$OUT/$2" "$B/$1" && echo "[done] $2 $(du -h $OUT/$2|cut -f1)"; }
get 51064667 CRISPRGeneEffect.csv
get 51064631 CRISPRGeneDependency.csv
get 51065489 OmicsExpressionProteinCodingGenesTPMLogp1.csv
get 51065732 OmicsSomaticMutations.csv
get 51065795 depmap_README.txt
echo "=== DEPMAP DONE ==="
