#!/usr/bin/env bash
# 步驟 199：同卡打包排程器 —— 一張卡上同時跑多個輕量工作
#
# 為什麼要另開一支而不改 195
# --------------------------
# 195 正在執行中，而 bash 是邊讀邊執行的，改動正在跑的 .sh 會讓它跳到錯位的位置。
# 所以 195 一行都不動，這支是它的「可打包」版本，兩者共用同一個加鎖佇列，
# 可以同時存在（flock 保證只有一個 worker 拿到同一行工作）。
#
# 為什麼可以打包（有實測，不是推測）
# ----------------------------------
# 2026-11-27 量測：130_residual_models.py 一個行程只用
#   顯存 962 MiB / 40960、功耗 86 W / 400、記憶體頻寬 10%、RSS 2.7 GB、CPU 213%
# 也就是一張 A100 有**四分之三的時間在等資料**。
# 併行安全性用實驗證實過，不是假設：把 seed 0 重跑成 _s0b 並刻意與 seed 2 共用 GPU 3，
# 四個方法的輸出與獨占整張卡時**逐位元相同**（最大絕對差 0.000e+00、排名完全一致、
# 前 100 Jaccard 1.0000）。所以打包只改變總時長，不改變任何一個工作的結果。
#
# 和 195 的唯一差別
# -----------------
# 195 的開工條件是「這張卡 memory.used ≤ 1000 MiB」—— 等於禁止任何共用，
# 第二個工作上去之後第三個就永遠等不到。
# 這支改成**顯存預算**：memory.used + 單一工作估計值 ≤ 預算才開工。
# 預算刻意設得保守（預設 16 GiB / 40 GiB），萬一別人也落到這張卡上，還留 24 GiB 給他。
#
# 用法
#   bash scripts/199_pack_daemon.sh "3 5" 2 16000     # 3 和 5 各加 2 個 worker
#   touch results/DAEMON_stop                          # 停（跑完手上的才退）
set -u
BASE=/raid/jupyter-alz.ee09/Tonsil_to_MCL
cd "$BASE"
Q=results/DAEMON_queue.txt
LOCK=results/DAEMON_queue.lock
STATUS=results/DAEMON_status.txt
STOP=results/DAEMON_stop

GPUS="${1:-3 5}"
PER_GPU="${2:-2}"
BUDGET_MIB="${3:-16000}"
EST_MIB="${EST_MIB:-1400}"     # 單一工作的顯存估計（實測 962 MiB + 餘裕）
IDLE_SLEEP=20
IDLE_REPORT=15

touch "$Q" "$LOCK"
say() { echo "[199] $(date +%H:%M:%S) $*" | tee -a "$STATUS"; }

# 不自行移除停止哨兵：若還掛著就拒絕啟動，避免把別人剛下的「停」悄悄取消
if [ -f "$STOP" ]; then
  echo "[199] 停止哨兵 $STOP 還在。要重新開工請先 rm 掉它。" >&2
  exit 1
fi

take_next() {
  ( flock 9
    s=$(head -1 "$Q" 2>/dev/null)
    [ -n "$s" ] && sed -i '1d' "$Q"
    printf '%s' "$s"
  ) 9<"$LOCK"
}

requeue_front() { ( flock 9; sed -i "1i $1" "$Q" ) 9<"$LOCK"; }

budget_ok() {           # 0 = 這張卡還塞得下一個工作
  u=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits -i "$1" 2>/dev/null)
  [ -n "$u" ] || return 1
  [ $((u + EST_MIB)) -le "$BUDGET_MIB" ]
}

worker() {
  g=$1; slot=$2; idle=0
  while :; do
    if [ -f "$STOP" ]; then say "GPU $g/$slot 收到停止哨兵，收工"; return; fi
    job=$(take_next)
    if [ -z "$job" ]; then
      idle=$((idle + 1))
      [ $((idle % IDLE_REPORT)) -eq 1 ] && say "GPU $g/$slot 待命中（佇列空，已等 $((idle * IDLE_SLEEP)) 秒）"
      sleep "$IDLE_SLEEP"; continue
    fi
    idle=0
    waited=0
    while ! budget_ok "$g"; do
      [ $((waited % 300)) -eq 0 ] && say "GPU $g/$slot 顯存預算不足（已用 $(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits -i "$g" 2>/dev/null) MiB，預算 ${BUDGET_MIB}），等待 ${waited}s"
      sleep 30; waited=$((waited + 30))
      if [ "$waited" -ge 3600 ]; then
        say "GPU $g/$slot 等超過 1 小時，把工作退回佇列前端"
        requeue_front "$job"; job=""; break
      fi
    done
    [ -z "$job" ] && continue
    cmd=${job//\{GPU\}/$g}
    say "GPU $g/$slot ← ${job:0:110}"
    t0=$(date +%s)
    bash -c "$cmd"; rc=$?
    say "GPU $g/$slot ✓ rc=$rc（$(( ($(date +%s) - t0) / 60 )) 分）${job:0:60}"
  done
}

say "打包守衛啟動：卡【$GPUS】× 每卡 $PER_GPU 個 worker，顯存預算 ${BUDGET_MIB} MiB／工作估 ${EST_MIB} MiB（佇列現有 $(wc -l < "$Q") 項）"
for g in $GPUS; do
  for s in $(seq 1 "$PER_GPU"); do
    worker "$g" "$s" &
    sleep 10          # 錯開啟動，避免同時讀同一批 h5ad 把 I/O 打爆
  done
done
wait
say "打包守衛全部 worker 收工"
