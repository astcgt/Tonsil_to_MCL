#!/usr/bin/env bash
# 步驟 224：捐贈者拆分（種子 0–4）到齊後自動跑 223 的判定。純 CPU，不佔 GPU worker。
set -u
cd /raid/jupyter-alz.ee09/Tonsil_to_MCL
NEED=5; MAX_MIN=${1:-90}; t0=$(date +%s)
count() { ls results/130_residual_by_method_ds[0-9].csv.gz 2>/dev/null | wc -l; }
echo "[224] $(date +%H:%M:%S) 等 $NEED 個捐贈者拆分結果，上限 ${MAX_MIN} 分鐘"
while :; do
  n=$(count); el=$(( ($(date +%s) - t0) / 60 ))
  if [ "$n" -ge "$NEED" ]; then
    echo "[224] $(date +%H:%M:%S) $n/$NEED 到齊（${el} 分），跑 223"
    env/tonsilmcl/bin/python -I scripts/223_donor_leverage.py 2>&1 | tee logs/223_donor.log
    echo "[224] $(date +%H:%M:%S) DONOR_ANALYSED rc=${PIPESTATUS[0]}"; exit 0
  fi
  if [ "$(pgrep -cf 'group-by-'donor)" -eq 0 ] && [ "$n" -lt "$NEED" ]; then
    echo "[224] $(date +%H:%M:%S) JOBS_GONE：沒有工作在跑，但只有 $n/$NEED 個檔案"
    for s in 0 1 2 3 4; do
      [ -f "results/130_residual_by_method_ds$s.csv.gz" ] || \
        { echo "--- ds$s log 尾段:"; tail -4 "logs/130_donorsplit_s$s.log" 2>/dev/null; }
    done
    exit 1
  fi
  [ "$el" -ge "$MAX_MIN" ] && { echo "[224] TIMEOUT：只到 $n/$NEED"; exit 1; }
  [ $((el % 10)) -eq 0 ] && echo "[224] $(date +%H:%M:%S) $n/$NEED（${el} 分）"
  sleep 60
done
