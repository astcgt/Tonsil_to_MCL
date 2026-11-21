#!/usr/bin/env bash
# 建立「基因體序列模型」環境（gReLU + Borzoi），與主分析環境隔離以免 torch 版本衝突
set -u
export MAMBA_NO_BANNER=1
export CONDA_PKGS_DIRS=/raid/jupyter-alz.ee09/Tonsil_to_MCL/env/pkgs
ENVDIR=/raid/jupyter-alz.ee09/Tonsil_to_MCL/env/genomic

mamba create -y -p "$ENVDIR" -c conda-forge python=3.11 pip 2>&1 | tail -5
"$ENVDIR/bin/pip" install --no-cache-dir -q \
  "torch==2.4.1" --index-url https://download.pytorch.org/whl/cu121 2>&1 | tail -3
"$ENVDIR/bin/pip" install --no-cache-dir -q gReLU borzoi-pytorch pyfaidx 2>&1 | tail -10

echo "=== 檢查 ==="
"$ENVDIR/bin/python" - <<'EOF'
import torch
print("torch", torch.__version__, "cuda", torch.cuda.is_available())
import grelu, grelu.resources, grelu.variant, grelu.interpret.score
print("gReLU OK")
try:
    import borzoi_pytorch; print("borzoi-pytorch OK")
except Exception as e:
    print("borzoi-pytorch:", e)
EOF
echo "=== GENOMIC ENV DONE ==="
