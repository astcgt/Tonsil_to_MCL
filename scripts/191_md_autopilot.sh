#!/usr/bin/env bash
# 步驟 191：MD 的自動接續。
# 等 tier1 的 6 段跑完 → 自動分析 → 依門檻判定決定要不要起 tier2。
# 判定規則寫在 scripts/190_md_interface.py 的 gate()，這裡只負責執行與配卡。
#
# v2 的三個修正：
#   1. 原本用 `pgrep -f 'tier 2 --gpu 3'` 等前一份工作結束，但等待迴圈自己的
#      指令列就含有這個字串，pgrep 永遠找得到自己 → 參考系統會無限等待。
#      改成 4 個 worker 共用一個加鎖的工作佇列，不再用 pgrep 判斷。
#   2. 卡從 2 張（3,4）擴到 4 張（3,4,5,6）。啟動前逐卡確認真的空著，
#      別人佔用的卡就跳過，worker 數量自動縮減。
#   3. tier2 系統表加入陰性對照（打亂序列），原本只排在跑不完的 200 ns tier3。
set -u
BASE=/raid/jupyter-alz.ee09/Tonsil_to_MCL
cd "$BASE"
PY=env/protdes/bin/python
STATE=results/190_md
GPUS_WANTED="3 4 5 6"
# 佇列順序＝配卡順序。主比較的兩個系統排最前面，確保各自一張卡；
# 第五個（陰性對照）等最先跑完的那張卡。
SYSTEMS="MAFG__D1_best MAFK__D1_best MAFF__D1_best MAFG__CTRLPOS_NRF2 MAFG__CTRLNEG_scram"

say() { echo "[191] $*"; }

say "等 tier1 的 6 段完成…"
while [ "$(ls $STATE/MAF{G,K}__D1_best/rep*/done.txt 2>/dev/null | wc -l)" -lt 6 ]; do
  sleep 120
done
say "tier1 全部完成，開始分析"

$PY scripts/190_md_interface.py --tier 1 --analyze 2>&1 | tail -40

V=$($PY -c "import json;print(json.load(open('results/190_md_tier1.json'))['gate'].get('verdict','?'))" 2>/dev/null)
say "門檻判定 = $V"

case "$V" in
  A4_both_dissociated|no_data|need_two_systems)
    say "判定為「停止」：不起 tier2。原因：$V"
    echo "$V" > results/190_md_autopilot_stopped.txt
    say "DONE"
    exit 0
    ;;
esac
say "判定為「續做」：啟動 tier2（5 系統 × 3 重複 × 30 ns）"

# --- 選卡：只用真的空著的卡，別人的工作不碰 -------------------------------
# tier1 的 python 寫完 done.txt 之後還會自己跑一次分析才退出，那段時間
# GPU 3/4 的記憶體還沒釋放。先等自己的卡放手，最多等 40 分鐘。
for g in 3 4; do
  for _ in $(seq 40); do
    u=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits -i "$g" 2>/dev/null)
    [ -n "$u" ] && [ "$u" -le 1000 ] && break
    sleep 60
  done
done

GPUS=""
for g in $GPUS_WANTED; do
  used=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits -i "$g" 2>/dev/null)
  if [ -z "$used" ]; then
    say "  GPU $g 查不到，跳過"
  elif [ "$used" -gt 1000 ]; then
    say "  GPU $g 已被佔用（${used} MiB），跳過"
  else
    GPUS="$GPUS $g"
  fi
done
GPUS=$(echo $GPUS)
if [ -z "$GPUS" ]; then
  say "沒有空的 GPU，中止"
  echo "no_free_gpu" > results/190_md_autopilot_stopped.txt
  exit 1
fi
say "  可用的卡：$GPUS"

# --- 加鎖的工作佇列 --------------------------------------------------------
Q="$STATE/t2_queue.txt"
LOCK="$STATE/t2_queue.lock"
printf '%s\n' $SYSTEMS > "$Q"
: > "$LOCK"

take_next() {
  (
    flock 9
    s=$(head -1 "$Q")
    [ -n "$s" ] && sed -i '1d' "$Q"
    printf '%s' "$s"
  ) 9<"$LOCK"
}

worker() {
  g=$1
  while :; do
    s=$(take_next)
    [ -z "$s" ] && break
    say "  GPU $g ← $s"
    $PY scripts/190_md_interface.py --tier 2 --gpu "$g" --only "$s" \
        > "logs/190_md_t2_${s}.log" 2>&1
    say "  GPU $g 完成 $s（return $?）"
  done
  say "  GPU $g 佇列空了，worker 收工"
}

for g in $GPUS; do
  worker "$g" &
  sleep 20          # 錯開建系統的時間，避免同時溶劑化吃爆記憶體
done
wait
say "tier2 全部跑完，做總分析"
$PY scripts/190_md_interface.py --tier 2 --analyze 2>&1 | tail -40
say "DONE"
