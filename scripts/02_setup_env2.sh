#!/usr/bin/env bash
# 第二階段：安裝分析用 Python 套件（GRN、參考映射、CNV、富集）
set -u
ENVDIR=/raid/jupyter-alz.ee09/Tonsil_to_MCL/env/tonsilmcl
PY="$ENVDIR/bin/python"
PIP="$ENVDIR/bin/pip"

$PIP install --no-cache-dir -q \
  "decoupler" "omnipath" \
  "scvi-tools" \
  "infercnvpy" \
  "gseapy" "celltypist" \
  "pyranges" "pybiomart" \
  "arboreto" "ctxcore" "pyscenic" \
  "adjustText" "upsetplot" "scikit-network" \
  2>&1 | tail -25

echo "=== 套件檢查 ==="
$PY - <<'EOF'
import importlib
for m in ["scanpy","anndata","decoupler","scvi","infercnvpy","gseapy",
          "celltypist","pyranges","pyscenic","arboreto","harmonypy","torch"]:
    try:
        mod = importlib.import_module(m)
        v = getattr(mod, "__version__", "?")
        print(f"{m:14s} OK   {v}")
    except Exception as e:
        print(f"{m:14s} FAIL {type(e).__name__}: {str(e)[:80]}")
EOF
echo "=== DONE ==="
