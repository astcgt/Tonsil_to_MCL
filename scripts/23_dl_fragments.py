#!/usr/bin/env python3
"""下載正常扁桃腺 Multiome/scATAC 的 fragments 檔（供與 MCL 做對等的覆蓋度比較）。"""
import os, sys, subprocess, concurrent.futures as cf, time

OUT = "/raid/jupyter-alz.ee09/Tonsil_to_MCL/data/raw_zenodo"
PART = os.path.join(OUT, "_parts")
REC, NAME = "8373756", "fragments_files.tar.gz"
NCHUNK = 12
os.makedirs(PART, exist_ok=True)
URL = f"https://zenodo.org/records/{REC}/files/{NAME}?download=1"


def size_of(u):
    r = subprocess.run(["curl", "-sI", "-L", "--max-time", "90", u],
                       capture_output=True, text=True)
    for line in r.stdout.splitlines():
        if line.lower().startswith("content-length:"):
            return int(line.split(":")[1].strip())
    return None


def need(path, lo, hi):
    have = os.path.getsize(path) if os.path.exists(path) else 0
    return (hi - lo + 1) - have


def get_chunk(a):
    u, path, lo, hi = a
    for _ in range(10):
        rem = need(path, lo, hi)
        if rem <= 0:
            return (path, True)
        have = (hi - lo + 1) - rem
        cmd = ["curl", "-sL", "--fail", "--retry", "3", "--retry-delay", "5",
               "--max-time", "10800", "-r", f"{lo+have}-{hi}"]
        if have:
            cmd += ["-C", "-"]
        cmd += ["-o", path, u]
        subprocess.run(cmd)
        if need(path, lo, hi) <= 0:
            return (path, True)
        time.sleep(5)
    return (path, False)


final = os.path.join(OUT, NAME)
total = size_of(URL)
print(f"[plan] {NAME} {total/1e9:.2f} GB", flush=True)
if os.path.exists(final) and os.path.getsize(final) == total:
    print("[skip] 已完整", flush=True); sys.exit(0)
if os.path.exists(final):
    os.remove(final)
step = total // NCHUNK + 1
parts, jobs = [], []
for i in range(NCHUNK):
    lo, hi = i * step, min(total - 1, (i + 1) * step - 1)
    if lo > hi:
        break
    p = os.path.join(PART, f"{NAME}.part{i:02d}")
    parts.append((p, lo, hi)); jobs.append((URL, p, lo, hi))

t0 = time.time()
with cf.ThreadPoolExecutor(max_workers=len(jobs)) as ex:
    res = list(ex.map(get_chunk, jobs))
if any(not ok for _, ok in res):
    print("[ERR] 有區塊失敗", flush=True); sys.exit(1)
print(f"[info] 下載完成 {(time.time()-t0)/60:.1f} min，合併中", flush=True)
with open(final, "wb") as fo:
    for p, lo, hi in parts:
        with open(p, "rb") as fi:
            while True:
                b = fi.read(1 << 24)
                if not b:
                    break
                fo.write(b)
got = os.path.getsize(final)
print(f"[done] {NAME} {got/1e9:.2f} GB {'OK' if got==total else 'MISMATCH'}", flush=True)
if got == total:
    for p, _, _ in parts:
        os.remove(p)
print("=== FRAGMENTS DOWNLOADED ===", flush=True)
