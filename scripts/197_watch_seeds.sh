#!/usr/bin/env bash
# 步驟 197：等 130 的五個種子都落地，立刻跑 196 的判定。
#
# 為什麼不塞進 195 的佇列：196 必須等**全部**種子完成。若當成佇列裡的一項工作，
# 會有一個 GPU worker 被它卡住空等（兩張卡、五個工作 → 最後一輪一定有人空轉）。
# 196 完全不用 GPU，所以獨立起一支純 CPU 的看守，兩張卡都不浪費。
#
# 這支腳本自己的檔名不會被 pgrep 整串匹配（190/194 踩過那個坑），因為它不用 pgrep。
set -u
cd /raid/jupyter-alz.ee09/Tonsil_to_MCL
NEED=5
MAX_MIN=${1:-150}
t0=$(date +%s)
echo "[197] $(date +%H:%M:%S) 等 $NEED 個種子，上限 ${MAX_MIN} 分鐘"
while :; do
  n=$(ls results/130_tumour_specific_by_method_s[0-9].csv 2>/dev/null | wc -l)
  el=$(( ($(date +%s) - t0) / 60 ))
  if [ "$n" -ge "$NEED" ]; then
    echo "[197] $(date +%H:%M:%S) $n/$NEED 到齊（${el} 分），跑 196"
    env/tonsilmcl/bin/python scripts/196_seed_stability.py 2>&1 | tee logs/196_seed_stability.log
    echo "[197] $(date +%H:%M:%S) SEEDS_ANALYSED rc=${PIPESTATUS[0]}"
    exit 0
  fi
  if [ "$el" -ge "$MAX_MIN" ]; then
    echo "[197] $(date +%H:%M:%S) TIMEOUT：${MAX_MIN} 分鐘只到 $n/$NEED"
    # 至少 3 個種子就還有離散度可談，照樣跑，讓 196 自己決定夠不夠
    if [ "$n" -ge 3 ]; then
      echo "[197] 用現有的 $n 個種子跑 196（196 自己會檢查夠不夠）"
      env/tonsilmcl/bin/python scripts/196_seed_stability.py 2>&1 | tee logs/196_seed_stability.log
      echo "[197] $(date +%H:%M:%S) SEEDS_ANALYSED_PARTIAL rc=${PIPESTATUS[0]}"
    fi
    exit 1
  fi
  [ $((el % 10)) -eq 0 ] && echo "[197] $(date +%H:%M:%S) $n/$NEED（${el} 分）"
  sleep 60
done
