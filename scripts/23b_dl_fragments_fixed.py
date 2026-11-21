#!/usr/bin/env python3
"""
下載正常扁桃腺 fragments（修正版）。

前一版的 bug：同時傳 `-r {lo+have}-{hi}` 與 `-C -`，curl 會把本地檔案大小
再加到 range 起點上（double counting），結果附加了重複的位元組，分塊過大。

修正：續傳時**不用 `-C -`**，改用「指定 range + 以 append 模式寫入」。
並在合併前逐一檢查每個分塊的大小是否完全正確。
"""
import os, sys, subprocess, concurrent.futures as cf, time

OUT = "/raid/jupyter-alz.ee09/Tonsil_to_MCL/data/raw_zenodo"
PART = os.path.join(OUT, "_parts")
REC, NAME = "8373756", "fragments_files.tar.gz"
NCHUNK = 12
URL = f"https://zenodo.org/records/{REC}/files/{NAME}?download=1"
os.makedirs(PART, exist_ok=True)


def log(*a):
    print("[23b]", *a, flush=True)


def size_of(u):
    r = subprocess.run(["curl", "-sI", "-L", "--max-time", "90", u],
                       capture_output=True, text=True)
    for line in r.stdout.splitlines():
        if line.lower().startswith("content-length:"):
            return int(line.split(":")[1].strip())
    return None


def get_chunk(a):
    u, path, lo, hi = a
    want = hi - lo + 1
    for attempt in range(12):
        have = os.path.getsize(path) if os.path.exists(path) else 0
        if have == want:
            return (path, True, "ok")
        if have > want:                       # 之前被寫壞 → 砍掉重來
            os.remove(path); have = 0
        # 關鍵：只用 range，不用 -C -，並以 append 模式寫入
        r = subprocess.run(
            ["curl", "-sL", "--fail", "--retry", "2", "--retry-delay", "5",
             "--max-time", "10800", "-r", f"{lo+have}-{hi}", u],
            stdout=open(path, "ab"))
        if os.path.getsize(path) == want:
            return (path, True, "ok")
        time.sleep(4)
    return (path, False, f"size={os.path.getsize(path)} want={want}")


def main():
    total = size_of(URL)
    final = os.path.join(OUT, NAME)
    log(f"目標大小 {total:,} bytes ({total/1e9:.2f} GB)")
    if os.path.exists(final) and os.path.getsize(final) == total:
        log("已完整"); return
    step = total // NCHUNK + 1
    parts, jobs = [], []
    for i in range(NCHUNK):
        lo, hi = i * step, min(total - 1, (i + 1) * step - 1)
        if lo > hi:
            break
        p = os.path.join(PART, f"{NAME}.part{i:02d}")
        want = hi - lo + 1
        cur = os.path.getsize(p) if os.path.exists(p) else 0
        status = "完整" if cur == want else ("過大→重抓" if cur > want else f"缺 {want-cur:,}")
        log(f"  part{i:02d}  目前 {cur:>13,} / 應為 {want:>13,}  [{status}]")
        parts.append((p, lo, hi)); jobs.append((URL, p, lo, hi))

    t0 = time.time()
    with cf.ThreadPoolExecutor(max_workers=NCHUNK) as ex:
        res = list(ex.map(get_chunk, jobs))
    bad = [(p, m) for p, ok, m in res if not ok]
    if bad:
        log("[ERR] 失敗:", bad); sys.exit(1)
    log(f"所有分塊就緒（{(time.time()-t0)/60:.1f} 分鐘），開始合併並驗證")

    with open(final, "wb") as fo:
        for p, lo, hi in parts:
            want = hi - lo + 1
            assert os.path.getsize(p) == want, (p, os.path.getsize(p), want)
            with open(p, "rb") as fi:
                while True:
                    b = fi.read(1 << 24)
                    if not b:
                        break
                    fo.write(b)
    got = os.path.getsize(final)
    ok = got == total
    log(f"合併完成 {got:,} / {total:,}  {'OK' if ok else 'MISMATCH'}")
    if ok:
        for p, _, _ in parts:
            os.remove(p)
        # gzip 完整性檢查
        r = subprocess.run(["gzip", "-t", final], capture_output=True, text=True)
        log("gzip 完整性:", "通過" if r.returncode == 0 else f"失敗 {r.stderr[:200]}")
    log("=== FRAGMENTS DONE ===")


main()
