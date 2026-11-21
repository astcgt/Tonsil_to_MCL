#!/usr/bin/env bash
# 步驟 38：用 ProteinMPNN 為 RFdiffusion3 產出的骨架設計序列
# 只設計 binder 鏈（B 鏈），標的 IRF4（A 鏈）保持原序列不動
set -u
BASE=/raid/jupyter-alz.ee09/Tonsil_to_MCL
ENVP=$BASE/env/protdes
IN=$BASE/results/37_rfd3_binders
OUT=$BASE/results/38_mpnn_seqs
mkdir -p "$OUT"
export CUDA_VISIBLE_DEVICES=${GPU:-3}

N=$(ls "$IN"/*.pdb 2>/dev/null | wc -l)
echo "找到 $N 個骨架"
if [ "$N" -eq 0 ]; then echo "沒有骨架可處理，先跑 37_run_rfd3_binder.sh"; exit 1; fi

i=0
for pdb in "$IN"/*.pdb; do
  i=$((i+1))
  nm=$(basename "$pdb" .pdb)
  "$ENVP/bin/mpnn" \
    --model_type protein_mpnn \
    --checkpoint_path "$BASE/models/foundry/proteinmpnn_v_48_020.pt" \
    --is_legacy_weights True \
    --structure_path "$pdb" \
    --out_directory "$OUT" \
    --name "$nm" \
    --fixed_chains A \
    --designed_chains B \
    --temperature 0.1 \
    --batch_size 8 \
    --number_of_batches 1 \
    --write_fasta True \
    --seed 0 \
    > "$OUT/${nm}.log" 2>&1 || echo "  [warn] $nm 失敗（見 $OUT/${nm}.log）"
  [ $((i % 10)) -eq 0 ] && echo "  已處理 $i / $N"
done
echo "=== 完成，產出 ==="
ls "$OUT"/*.fa 2>/dev/null | wc -l
echo "=== MPNN DONE ==="
