#!/usr/bin/env bash
# 解壓並盤點下載到的資料結構
set -u
BASE=/raid/jupyter-alz.ee09/Tonsil_to_MCL
RAW=$BASE/data/raw_zenodo
EXT=$BASE/data/extracted
mkdir -p "$EXT"

echo "############ 1. MCL.zip 內容 ############"
unzip -l "$RAW/MCL.zip" | head -60
echo "...(total)"; unzip -l "$RAW/MCL.zip" | tail -3

echo
echo "############ 2. TonsilAtlasSeuratRNA.tar.gz 內容 ############"
tar tzf "$RAW/TonsilAtlasSeuratRNA.tar.gz" | head -40
echo "...(n files)"; tar tzf "$RAW/TonsilAtlasSeuratRNA.tar.gz" | wc -l

echo
echo "############ 3. TonsilAtlasSeuratMultiome.tar.gz 內容 ############"
tar tzf "$RAW/TonsilAtlasSeuratMultiome.tar.gz" | head -40
echo "...(n files)"; tar tzf "$RAW/TonsilAtlasSeuratMultiome.tar.gz" | wc -l
