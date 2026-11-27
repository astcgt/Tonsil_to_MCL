#!/usr/bin/env bash
# 步驟 192：GPU 工作排程器（取代 191_md_autopilot.sh）
#
# 為什麼不用 191：191 把「一個系統」綁在「一張卡」上，又用 pgrep 判斷前一份
# 工作有沒有結束（pgrep -f 會匹配到等待迴圈自己）。這一版改成：
#   * 工作單位是一行可執行的指令，裡面用 {GPU} 當佔位符
#   * 4 個 worker 各綁一張卡，從同一個加鎖佇列拉工作
#   * 佇列是一個「活的」檔案：執行中可以追加工作，worker 會接著拉
#
# 階段：
#   1. tier1（10 ns × 3 重複 × 5 系統，含正負對照）—— 已完成的 rep 會秒回跳過
#   2. 分析 + 門檻判定。陰性對照若沒散掉就停，不往下燒 30 ns
#   3. 判定為續做 → 排入 tier2（30 ns）；否則只排獨立工作
set -u
BASE=/raid/jupyter-alz.ee09/Tonsil_to_MCL
cd "$BASE"
PY=env/protdes/bin/python
Q=results/SCHED_queue.txt
LOCK=results/SCHED_queue.lock
STATUS=results/SCHED_status.txt
GPUS_WANTED="3 4 5 6"
SYSTEMS="MAFG__D1_best MAFK__D1_best MAFF__D1_best MAFG__CTRLPOS_NRF2 MAFG__CTRLNEG_scram"

say() { echo "[192] $(date +%H:%M:%S) $*"; echo "[192] $(date +%H:%M:%S) $*" >> "$STATUS"; }

: > "$LOCK"
: > "$STATUS"

take_next() {
  ( flock 9
    s=$(head -1 "$Q" 2>/dev/null)
    [ -n "$s" ] && sed -i '1d' "$Q"
    printf '%s' "$s"
  ) 9<"$LOCK"
}

enqueue() {            # 可在執行中追加，worker 會接著拉
  ( flock 9; printf '%s\n' "$1" >> "$Q" ) 9<"$LOCK"
}

pick_gpus() {
  GPUS=""
  for g in $GPUS_WANTED; do
    u=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits -i "$g" 2>/dev/null)
    if [ -z "$u" ]; then say "  GPU $g 查不到，跳過"
    elif [ "$u" -gt 1000 ]; then say "  GPU $g 被佔用（${u} MiB），跳過"
    else GPUS="$GPUS $g"; fi
  done
  GPUS=$(echo $GPUS)
}

worker() {
  g=$1
  while :; do
    job=$(take_next)
    [ -z "$job" ] && break
    cmd=${job//\{GPU\}/$g}
    say "  GPU $g ← ${job:0:110}"
    bash -c "$cmd"
    say "  GPU $g ✓ rc=$? ${job:0:70}"
  done
  say "  GPU $g worker 收工"
}

drain() {              # 起 worker 把佇列清空
  pick_gpus
  if [ -z "$GPUS" ]; then say "沒有空的 GPU，中止"; exit 1; fi
  say "可用的卡：$GPUS（佇列 $(wc -l < "$Q") 項）"
  for g in $GPUS; do worker "$g" & sleep 15; done
  wait
}

# ---------------- 階段 1：tier1 把五系統補齊到 10 ns -----------------------
say "階段 1：tier1（10 ns × 3 重複 × 5 系統）"
: > "$Q"
for r in 1 2 3; do
  for s in $SYSTEMS; do
    enqueue "$PY scripts/190_md_interface.py --tier 1 --gpu {GPU} --only $s --rep $r >> logs/190_t1_${s}_rep${r}.log 2>&1"
  done
done
drain

# ---------------- 階段 2：分析與門檻判定 -----------------------------------
say "階段 2：tier1 分析"
$PY scripts/190_md_interface.py --tier 1 --analyze 2>&1 | tail -40 | tee -a "$STATUS"
V=$($PY -c "import json;print(json.load(open('results/190_md_tier1.json'))['gate'].get('verdict','?'))" 2>/dev/null)
say "門檻判定 = $V"

# ---------------- 階段 3：依判定排 tier2，獨立工作一律排在後面 -------------
: > "$Q"
case "$V" in
  A4_both_dissociated|A5_assay_cannot_discriminate|no_data|need_two_systems)
    say "判定為「停止」：不做 tier2。原因：$V"
    echo "$V" > results/190_md_stopped.txt
    ;;
  *)
    say "判定為「續做」：排入 tier2（30 ns × 3 重複 × 5 系統）"
    for r in 1 2 3; do
      for s in $SYSTEMS; do
        enqueue "$PY scripts/190_md_interface.py --tier 2 --gpu {GPU} --only $s --rep $r >> logs/190_t2_${s}_rep${r}.log 2>&1"
      done
    done
    ;;
esac

if [ -s results/SCHED_aux_queue.txt ]; then
  say "追加獨立工作 $(wc -l < results/SCHED_aux_queue.txt) 項"
  while IFS= read -r l; do [ -n "$l" ] && enqueue "$l"; done < results/SCHED_aux_queue.txt
fi

if [ -s "$Q" ]; then
  drain
  case "$V" in
    A4_both_dissociated|A5_assay_cannot_discriminate|no_data|need_two_systems) ;;
    *) say "tier2 總分析"
       $PY scripts/190_md_interface.py --tier 2 --analyze 2>&1 | tail -40 | tee -a "$STATUS" ;;
  esac
fi
say "DONE"
