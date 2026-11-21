#!/usr/bin/env bash
# 步驟 76：第二輪結合劑設計 —— 熱點改用「ESM 認定零逃逸」的殘基
# 這樣設計出來的結合劑，其接觸點在演化上完全不能突變，理論上最難產生抗藥性
set -u
BASE=/raid/jupyter-alz.ee09/Tonsil_to_MCL
OUT=$BASE/results/76_rfd3_v2; mkdir -p "$OUT"
export CUDA_VISIBLE_DEVICES=${GPU:-3}
cat > "$OUT/input.json" <<'JSON'
{
  "mafg_zeroescape_70": {
    "input": "/raid/jupyter-alz.ee09/Tonsil_to_MCL/results/56_mafg_target.pdb",
    "contig": "A23-121,70",
    "select_hotspots": "A53,A56,A57,A58,A60,A61,A62,A64,A65,A68,A69,A71,A72"
  },
  "mafg_dual_90": {
    "input": "/raid/jupyter-alz.ee09/Tonsil_to_MCL/results/56_mafg_target.pdb",
    "contig": "A23-121,90",
    "select_hotspots": "A53,A56,A57,A60,A61,A64,A69,A75,A79,A80,A93,A97"
  },
  "fsp1_pocketmouth_80": {
    "input": "/raid/jupyter-alz.ee09/Tonsil_to_MCL/results/62_fsp1_target.pdb",
    "contig": "A7-373,80",
    "select_hotspots": "A18,A19,A20,A22,A23,A41,A43,A49,A53,A69,A80,A108,A110"
  },
  "fsp1_nadp_80": {
    "input": "/raid/jupyter-alz.ee09/Tonsil_to_MCL/results/62_fsp1_target.pdb",
    "contig": "A7-373,80",
    "select_hotspots": "A149,A150,A151,A154,A156,A174,A175,A181,A207,A212,A243,A244,A293"
  }
}
JSON
echo "=== 開始第二輪擴散取樣（GPU $CUDA_VISIBLE_DEVICES）==="
"$BASE/env/protdes/bin/rfd3" out_dir="$OUT" inputs="$OUT/input.json" \
  diffusion_batch_size=4 n_batches=4 2>&1 | tail -20
echo "=== 產出 $(ls "$OUT"/*.cif.gz 2>/dev/null | wc -l) 個設計 ==="
