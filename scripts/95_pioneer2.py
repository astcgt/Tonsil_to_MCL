#!/usr/bin/env python
"""
步驟 95：抓取 PIONEER2 全互作體介面檔，但**不落地完整 9.3 GB**。

為什麼不直接下載整個檔
----------------------
專案已佔 122 GB，自訂上限 150 GB。完整檔 9,991,839,205 bytes（9.3 GiB）
落地後還要再讀一次，對磁碟不划算。而我們只需要 35 個蛋白質相關的行。

做法：分段下載 → 立即過濾 → 刪掉原始分段
-----------------------------------------
* 160 個分段，每段約 62 MB，同時最多 10 個 worker
  → 峰值磁碟約 620 MB，而不是 9.3 GB。
* 行邊界處理（這是唯一容易出錯的地方）：
  每段多抓 OVERLAP bytes，然後
    - 除了第 0 段，丟掉第一個換行符之前的殘行；
    - 只保留到「超過本段原始終點的第一個換行符」為止。
  這樣每一行剛好被一個分段完整收下，不重複也不遺漏。
* 每段可獨立重試；中斷後重跑只會補未完成的段。

輸出
----
results/95_pioneer2_axis.tsv.gz  只含我們 35 個蛋白質參與的交互作用
results/95_pioneer2_stats.csv    逐段統計（可核對覆蓋率）
"""
import os, sys, gzip, io, json, subprocess, tempfile, time
from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
RES = f"{BASE}/results"
# 分段下載的暫存目錄。這些分段加起來約 9.3 GB，所以放在暫存區而非倉庫裡；
# 要指到別的磁碟時用環境變數 P2_CHUNK_DIR 覆寫。
TMP = os.environ.get("P2_CHUNK_DIR") or os.path.join(tempfile.gettempdir(), "p2chunks")
URL = "https://pioneer2.yulab.org/static/PIONEER2_interactome_wide_interfaces.txt"
TOTAL = 9991839205
NCHUNK = 160
WORKERS = 10
OVERLAP = 262144          # 256 KB，遠大於任何一行（最長的蛋白質序列約 35 kB）

os.makedirs(TMP, exist_ok=True)
ids = pd.read_csv(f"{RES}/94_uniprot_ids.csv")
KEEP = set(ids.uniprot)
G = dict(zip(ids.uniprot, ids.gene))


def log(*a):
    print("[95]", *a, flush=True)


def fetch(lo, hi):
    """抓 [lo, hi] 位元組，失敗重試 3 次。"""
    for attempt in range(3):
        r = subprocess.run(
            ["curl", "-sL", "--fail", "--max-time", "900",
             "-r", f"{lo}-{hi}", URL],
            stdout=subprocess.PIPE)
        if r.returncode == 0 and r.stdout:
            return r.stdout
        time.sleep(3 * (attempt + 1))
    return None


def do_chunk(i):
    part = f"{TMP}/hit_{i:04d}.tsv"
    meta = f"{TMP}/hit_{i:04d}.json"
    if os.path.exists(meta):                      # 已完成的段直接跳過（可續跑）
        m = json.load(open(meta))
        return i, part, m["lines"], m["bytes"]
    lo = i * (TOTAL // NCHUNK)
    hi = TOTAL - 1 if i == NCHUNK - 1 else (i + 1) * (TOTAL // NCHUNK) - 1
    end_guard = hi - lo                     # 本段原始終點（相對位置）
    raw = fetch(lo, min(hi + OVERLAP, TOTAL - 1))
    if raw is None:
        return i, None, 0, 0
    # 行邊界：丟掉開頭殘行
    start = 0
    if i > 0:
        nl = raw.find(b"\n")
        if nl < 0:
            return i, [], 0, len(raw)
        start = nl + 1
    # 行邊界：保留到超過原始終點的第一個換行符
    nl = raw.find(b"\n", end_guard)
    stop = len(raw) if nl < 0 else nl + 1
    body = raw[start:stop]

    nbytes = len(body)
    nlines = 0
    with open(part, "wb") as fh:
        for line in io.BytesIO(body):             # 逐行迭代，不一次展開成 list
            line = line.rstrip(b"\n")
            if not line:
                continue
            nlines += 1
            f = line.split(b"\t", 3)
            if len(f) < 3:
                continue
            if f[0].decode("ascii", "ignore") in KEEP or f[1].decode("ascii", "ignore") in KEEP:
                fh.write(line + b"\n")
    del raw, body
    json.dump(dict(lines=nlines, bytes=nbytes), open(meta, "w"))
    return i, part, nlines, nbytes


log(f"檔案 {TOTAL:,} bytes，切成 {NCHUNK} 段（每段約 {TOTAL//NCHUNK/1e6:.0f} MB），"
    f"{WORKERS} 個 worker，峰值磁碟約 {WORKERS*(TOTAL//NCHUNK+OVERLAP)/1e6:.0f} MB")
log(f"過濾對象：{len(KEEP)} 個 UniProt ID")

out_path = f"{RES}/95_pioneer2_axis.tsv.gz"
stats, parts, nhits_total = [], [], [0]
t0 = time.time()
done = 0
with ThreadPoolExecutor(max_workers=WORKERS) as ex:
    futs = {ex.submit(do_chunk, i): i for i in range(NCHUNK)}
    for fu in as_completed(futs):
        i, part, nlines, nbytes = fu.result()
        done += 1
        if part is None:
            log(f"  段 {i} 失敗（3 次重試後）")
            stats.append(dict(chunk=i, ok=False, lines=0, bytes=0, hits=0))
            continue
        nh = sum(1 for _ in open(part, "rb")) if os.path.exists(part) else 0
        parts.append(part)
        nhits_total[0] += nh
        stats.append(dict(chunk=i, ok=True, lines=nlines, bytes=nbytes, hits=nh))
        if done % 20 == 0 or done == NCHUNK:
            el = time.time() - t0
            tb = sum(s["bytes"] for s in stats)
            log(f"  {done}/{NCHUNK} 段 | 已掃 {tb/1e9:.2f} GB | "
                f"命中 {nhits_total[0]:,} 行 | {tb/1e6/max(el,1):.0f} MB/s | {el/60:.1f} 分")

st = pd.DataFrame(stats).sort_values("chunk")
st.to_csv(f"{RES}/95_pioneer2_stats.csv", index=False)
cov = st.bytes.sum()
log(f"\n覆蓋率檢查：掃過 {cov:,} / {TOTAL:,} bytes = {100*cov/TOTAL:.3f}%")
log(f"失敗段數：{(~st.ok).sum()}　總行數：{st.lines.sum():,}　命中行數：{nhits_total[0]:,}")

HDR = b"protein1\tprotein2\tsource\tprotein1_interface\tprotein2_interface\tprotein1_sequence\tprotein2_sequence"
n_written = 0
with gzip.open(out_path, "wb") as fh:
    fh.write(HDR + b"\n")
    for part in sorted(parts):
        for line in open(part, "rb"):
            if not line.startswith(b"protein1\t"):
                fh.write(line); n_written += 1
log(f"寫出 {n_written:,} 行")
log(f"寫出 {out_path}（{os.path.getsize(out_path)/1e6:.1f} MB）")
log("DONE")
