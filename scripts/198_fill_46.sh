#!/usr/bin/env bash
# 步驟 198：192 收工後接棒守住它放掉的卡（通常是 4 與 6）。
#
# 為什麼需要：192 的 worker 在佇列見底時 break 退出，192 主行程跑完階段 2／3 就結束，
# 它佔的卡從此沒有任何排程器看守。這支在 192 真正結束後，把空出來的卡交給 199 打包守衛。
# 用精確 PID 而非 pgrep 模式比對（pgrep -f 會匹配到本腳本自己的命令列）。
set -u
cd /raid/jupyter-alz.ee09/Tonsil_to_MCL
PID192=${1:-3311792}
CARDS=${2:-"4 6"}
PER_GPU=${3:-2}
MAX_MIN=${4:-240}
t0=$(date +%s)
log() { echo "[198] $(date +%H:%M:%S) $*" | tee -a results/DAEMON_status.txt; }

log "等 192（pid $PID192）結束，之後接管【$CARDS】，上限 ${MAX_MIN} 分鐘"
while kill -0 "$PID192" 2>/dev/null; do
  el=$(( ($(date +%s) - t0) / 60 ))
  if [ "$el" -ge "$MAX_MIN" ]; then log "TIMEOUT：192 仍在跑（${el} 分），不接棒"; exit 1; fi
  sleep 30
done
log "192 已結束（等了 $(( ($(date +%s) - t0) / 60 )) 分）"

V=$(env/protdes/bin/python -c "import json;print(json.load(open('results/190_md_tier1.json'))['gate'].get('verdict','?'))" 2>/dev/null || echo "?")
log "tier1 門檻判定 = $V"
grep -E '門檻判定|判定為' results/SCHED_status.txt 2>/dev/null | tail -3

TAKE=""
for g in $CARDS; do
  u=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits -i "$g" 2>/dev/null)
  if [ -n "$u" ] && [ "$u" -le 2000 ]; then TAKE="$TAKE $g"
  else log "GPU $g 仍被佔用（${u} MiB），不接管"; fi
done
TAKE=$(echo $TAKE)
[ -z "$TAKE" ] && { log "沒有卡可接管"; exit 1; }
log "接管【$TAKE】× 每卡 $PER_GPU 個 worker"
exec bash scripts/199_pack_daemon.sh "$TAKE" "$PER_GPU" 16000
