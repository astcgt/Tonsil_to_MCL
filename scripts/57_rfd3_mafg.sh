#!/usr/bin/env bash
# 步驟 57：設計阻斷 MAFG 的迷你蛋白質結合劑（兩種機制各設計一批）
set -u
BASE=/raid/jupyter-alz.ee09/Tonsil_to_MCL
OUT=$BASE/results/57_rfd3_mafg; mkdir -p "$OUT"
export CUDA_VISIBLE_DEVICES=${GPU:-3}
cat > "$OUT/input.json" <<'JSON'
{
  "mafg_dimer_block_70": {
    "input": "/raid/jupyter-alz.ee09/Tonsil_to_MCL/results/56_mafg_target.pdb",
    "contig": "A23-121,70",
    "select_hotspots": "A75,A76,A79,A80,A82,A83,A86,A87,A89,A90,A93,A94,A97,A100"
  },
  "mafg_dimer_block_90": {
    "input": "/raid/jupyter-alz.ee09/Tonsil_to_MCL/results/56_mafg_target.pdb",
    "contig": "A23-121,90",
    "select_hotspots": "A75,A76,A79,A80,A82,A83,A86,A87,A89,A90,A93,A94,A97,A100"
  },
  "mafg_dna_block_70": {
    "input": "/raid/jupyter-alz.ee09/Tonsil_to_MCL/results/56_mafg_target.pdb",
    "contig": "A23-121,70",
    "select_hotspots": "A53,A56,A57,A60,A61,A62,A69,A71,A72"
  },
  "mafg_dna_block_90": {
    "input": "/raid/jupyter-alz.ee09/Tonsil_to_MCL/results/56_mafg_target.pdb",
    "contig": "A23-121,90",
    "select_hotspots": "A53,A56,A57,A60,A61,A62,A69,A71,A72"
  }
}
JSON
echo "=== 開始擴散取樣（GPU $CUDA_VISIBLE_DEVICES）==="
"$BASE/env/protdes/bin/rfd3" out_dir="$OUT" inputs="$OUT/input.json" \
  diffusion_batch_size=4 n_batches=4 2>&1 | tail -25
echo "=== 產出 $(ls "$OUT"/*.cif.gz 2>/dev/null | wc -l) 個設計 ==="
