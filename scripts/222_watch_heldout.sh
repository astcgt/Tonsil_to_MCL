#!/usr/bin/env bash
# 步驟 222：保留集（種子 5–10）到齊後，自動跑 220 與 221 的確認。
#
# 為什麼獨立一支：220/221 完全不用 GPU，塞進 GPU 佇列會白佔一張卡空等。
# 不用 pgrep（會匹配到自己），改用「檔案數」當完成訊號，並明確涵蓋失敗路徑。
set -u
cd /raid/jupyter-alz.ee09/Tonsil_to_MCL
NEED=6
MAX_MIN=${1:-90}
t0=$(date +%s)
PY="env/tonsilmcl/bin/python -I"
count() { ls results/130_residual_by_method_s[5-9].csv.gz \
             results/130_residual_by_method_s10.csv.gz 2>/dev/null | wc -l; }
echo "[222] $(date +%H:%M:%S) 等保留集 $NEED 個種子，上限 ${MAX_MIN} 分鐘"
while :; do
  n=$(count); el=$(( ($(date +%s) - t0) / 60 ))
  if [ "$n" -ge "$NEED" ]; then
    echo "[222] $(date +%H:%M:%S) 保留集 $n/$NEED 到齊（${el} 分），跑確認"
    $PY scripts/220_method_decomposition.py 2>&1 | tee logs/220_confirm.log
    r1=${PIPESTATUS[0]}
    $PY scripts/221_tumour_term_controls.py 2>&1 | tee logs/221_confirm.log
    r2=${PIPESTATUS[0]}
    echo "[222] $(date +%H:%M:%S) HELDOUT_CONFIRMED rc=$r1/$r2"
    exit 0
  fi
  # 跑中的 130 行程歸零但檔案沒到齊 → 有工作死了，不要無聲等到 timeout
  if [ "$(pgrep -cf '130_residual_mo''dels.py')" -eq 0 ] && [ "$n" -lt "$NEED" ]; then
    echo "[222] $(date +%H:%M:%S) JOBS_GONE：沒有 130 在跑，但只有 $n/$NEED 個檔案"
    for s in 5 6 7 8 9 10; do
      [ -f "results/130_residual_by_method_s$s.csv.gz" ] || \
        { echo "--- 種子 $s 的 log 尾段:"; tail -5 "logs/130_seed_s$s.log" 2>/dev/null; }
    done
    exit 1
  fi
  if [ "$el" -ge "$MAX_MIN" ]; then
    echo "[222] $(date +%H:%M:%S) TIMEOUT：${MAX_MIN} 分鐘只到 $n/$NEED"; exit 1
  fi
  [ $((el % 10)) -eq 0 ] && echo "[222] $(date +%H:%M:%S) $n/$NEED（${el} 分）"
  sleep 60
done
