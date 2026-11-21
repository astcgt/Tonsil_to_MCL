#!/usr/bin/env bash
# 步驟 37：用 RFdiffusion3 設計阻斷 IRF4 DNA 結合面的迷你蛋白質結合劑
# 標的：PDB 7OGS 的 A 鏈（IRF4 DNA 結合域，殘基 20–133）
# 熱點：由步驟 36 從共結晶結構算出的 DNA 接觸殘基
set -u
BASE=/raid/jupyter-alz.ee09/Tonsil_to_MCL
ENVP=$BASE/env/protdes
OUT=$BASE/results/37_rfd3_binders
mkdir -p "$OUT"
export CUDA_VISIBLE_DEVICES=${GPU:-3}
export FOUNDRY_CHECKPOINT_DIR=$BASE/models/foundry

cat > "$OUT/input.json" <<'JSON'
{
  "irf4_dbd_blocker_60": {
    "input": "/raid/jupyter-alz.ee09/Tonsil_to_MCL/results/36_irf4_target.pdb",
    "contig": "A20-133,60",
    "select_hotspots": "A23,A24,A54,A56,A57,A74,A78,A91,A94,A95,A96,A98,A99,A102,A103"
  },
  "irf4_dbd_blocker_80": {
    "input": "/raid/jupyter-alz.ee09/Tonsil_to_MCL/results/36_irf4_target.pdb",
    "contig": "A20-133,80",
    "select_hotspots": "A23,A24,A54,A56,A57,A74,A78,A91,A94,A95,A96,A98,A99,A102,A103"
  },
  "irf4_dbd_blocker_100": {
    "input": "/raid/jupyter-alz.ee09/Tonsil_to_MCL/results/36_irf4_target.pdb",
    "contig": "A20-133,100",
    "select_hotspots": "A23,A24,A54,A56,A57,A74,A78,A91,A94,A95,A96,A98,A99,A102,A103"
  }
}
JSON

echo "=== 輸入規格 ==="; cat "$OUT/input.json"
echo "=== 開始擴散取樣（GPU $CUDA_VISIBLE_DEVICES）==="
"$ENVP/bin/rfd3" \
  out_dir="$OUT" \
  inputs="$OUT/input.json" \
  diffusion_batch_size=4 \
  n_batches=3 \
  2>&1 | tail -60

echo "=== 產出 ==="
ls -la "$OUT" | head -30
