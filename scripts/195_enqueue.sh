#!/usr/bin/env bash
# 把一行工作追加到 195 守衛的活佇列。指令裡用 {GPU} 當卡號佔位符。
set -u
cd /raid/jupyter-alz.ee09/Tonsil_to_MCL
LOCK=results/DAEMON_queue.lock; Q=results/DAEMON_queue.txt
touch "$Q" "$LOCK"
( flock 9; printf '%s\n' "$1" >> "$Q" ) 9<"$LOCK"
echo "已排入（佇列共 $(wc -l < "$Q") 項）：${1:0:100}"
