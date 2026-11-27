#!/usr/bin/env bash
# 步驟 194：等 MD tier1 的 15 段全部完成就退出，讓上層立刻接手做門檻判定。
#
# 寫成獨立檔案而不是 inline 字串，是因為 pgrep -f 會匹配到自己的指令列：
# 如果把 "190_md_interface.py" 寫在 -c 的字串裡，這支監看程式自己就會被數進去，
# 於是「還有沒有 MD 在跑」永遠是 true。檔案化之後指令列只有腳本路徑。
#
# 退出條件有三個，沉默不算正常：
#   1. 15/15 完成        → TIER1_COMPLETE
#   2. 沒有 MD 程序在跑但佇列裡還有 MD 工作 → STALL（worker 掛了）
#   3. 超過上限時間      → TIMEOUT
set -u
cd /raid/jupyter-alz.ee09/Tonsil_to_MCL
PAT=190_md_inter'face.py'          # 拆開寫，連這支腳本自己的檔案內容都不會整串匹配
MAXMIN=${1:-480}
t0=$(date +%s)
while :; do
  n=$(ls results/190_md/*/10ns_rep*/done.txt 2>/dev/null | wc -l)
  if [ "$n" -ge 15 ]; then
    echo "TIER1_COMPLETE $n/15 —— 15 段 10 ns 全部跑完，可以做門檻判定了"
    exit 0
  fi
  procs=$(ps -eo pid,args | grep -c "[/ ]$PAT" || true)
  qmd=$(grep -c "$PAT" results/SCHED_queue.txt 2>/dev/null || echo 0)
  if [ "$procs" -eq 0 ] && [ "$qmd" -gt 0 ]; then
    echo "STALL 完成 $n/15，但沒有 MD 程序在跑、佇列裡還有 $qmd 個 MD 工作 —— worker 可能掛了"
    exit 2
  fi
  if [ $(( ($(date +%s) - t0) / 60 )) -ge "$MAXMIN" ]; then
    echo "TIMEOUT 等了 $MAXMIN 分鐘，只完成 $n/15"
    exit 3
  fi
  sleep 60
done
