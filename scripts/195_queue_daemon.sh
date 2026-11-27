#!/usr/bin/env bash
# 步驟 195：常駐佇列守衛（修掉 192 的「抽乾就收工」）
#
# 192 的真正缺口
# --------------
# 192 的 worker 是 `job=$(take_next); [ -z "$job" ] && break`。
# 佇列一旦抽乾，worker 就**永久退出**。15 個 tier1 + 3 個獨立工作被 4 張卡吃完之後，
# GPU 3 的 worker 在 15:42:32 收工、GPU 5 在 15:40:25 收工，而 GPU 4、6 還在跑長工作。
# 此時再往佇列追加工作，3 和 5 已經沒有 worker 會去拉 → 兩張卡空著。
# 「無縫銜接」在佇列有料時完全正常（證據：15:14:53 GPU 5 做完 MD 立刻接 132、
# 15:36:52 接 141，中間 0 秒），壞的是**佇列見底之後沒有人守著**。
#
# 這一版的差別
# ------------
#   * worker 拉到空佇列時 **睡 20 秒再拉**，不退出
#   * 只有看到停止哨兵檔 results/DAEMON_stop 才收工 → 可以隨時追加工作
#   * 每次啟動工作前先確認那張卡是空的（別人佔住就等，不硬擠上去）
#   * 佇列檔**獨立於 192**（192 的階段 3 會 `: > $Q` 清空它自己的佇列，
#     共用同一個檔案會被它清掉）
#
# 用法
#   bash scripts/195_queue_daemon.sh 3 5          # 守 GPU 3 與 5
#   bash scripts/195_enqueue.sh "<指令含 {GPU}>"  # 隨時追加
#   touch results/DAEMON_stop                     # 收工（跑完手上的才退）
set -u
BASE=/raid/jupyter-alz.ee09/Tonsil_to_MCL
cd "$BASE"
Q=results/DAEMON_queue.txt
LOCK=results/DAEMON_queue.lock
STATUS=results/DAEMON_status.txt
STOP=results/DAEMON_stop
IDLE_SLEEP=20
IDLE_REPORT=15          # 每 15 次空轉（約 5 分鐘）才報一次，避免 log 爆掉

GPUS="${*:-3 5}"
touch "$Q" "$LOCK"
rm -f "$STOP"

say() { echo "[195] $(date +%H:%M:%S) $*" | tee -a "$STATUS"; }

take_next() {
  ( flock 9
    s=$(head -1 "$Q" 2>/dev/null)
    [ -n "$s" ] && sed -i '1d' "$Q"
    printf '%s' "$s"
  ) 9<"$LOCK"
}

gpu_free() {            # 回傳 0 = 這張卡可用
  u=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits -i "$1" 2>/dev/null)
  [ -n "$u" ] && [ "$u" -le 1000 ]
}

worker() {
  g=$1
  idle=0
  while :; do
    if [ -f "$STOP" ]; then say "GPU $g 收到停止哨兵，收工"; return; fi
    job=$(take_next)
    if [ -z "$job" ]; then
      idle=$((idle + 1))
      [ $((idle % IDLE_REPORT)) -eq 1 ] && say "GPU $g 待命中（佇列空，已等 $((idle * IDLE_SLEEP)) 秒）"
      sleep "$IDLE_SLEEP"
      continue
    fi
    idle=0
    # 卡被別人佔住就先把工作放回佇列前端，等它空出來，不要硬擠
    waited=0
    while ! gpu_free "$g"; do
      [ $((waited % 300)) -eq 0 ] && say "GPU $g 被佔用，等待中（${waited}s）"
      sleep 30; waited=$((waited + 30))
      if [ "$waited" -ge 3600 ]; then
        say "GPU $g 佔用超過 1 小時，把工作退回佇列"
        ( flock 9; sed -i "1i $job" "$Q" ) 9<"$LOCK"
        job=""; break
      fi
    done
    [ -z "$job" ] && continue
    cmd=${job//\{GPU\}/$g}
    say "GPU $g ← ${job:0:120}"
    t0=$(date +%s)
    bash -c "$cmd"; rc=$?
    say "GPU $g ✓ rc=$rc（$(( ($(date +%s) - t0) / 60 )) 分）${job:0:70}"
  done
}

say "守衛啟動，卡：$GPUS（佇列現有 $(wc -l < "$Q") 項）"
for g in $GPUS; do worker "$g" & sleep 5; done
wait
say "全部 worker 收工"
