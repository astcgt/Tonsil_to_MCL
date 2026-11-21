#!/usr/bin/env python3
"""Zenodo 多連線分塊下載器（安靜輸出、可續傳、合併後驗證大小）。"""
import os, sys, subprocess, concurrent.futures as cf, time

OUT = "/raid/jupyter-alz.ee09/Tonsil_to_MCL/data/raw_zenodo"
PART = os.path.join(OUT, "_parts")
FILES = [
    ("10373041", "MCL.zip"),
    ("8373756",  "TonsilAtlasSeuratRNA.tar.gz"),
    ("8373756",  "TonsilAtlasSeuratMultiome.tar.gz"),
]
NCHUNK = 10
os.makedirs(PART, exist_ok=True)


def url(rec, name):
    return "https://zenodo.org/records/%s/files/%s?download=1" % (rec, name)


def size_of(u):
    r = subprocess.run(["curl", "-sI", "-L", "--max-time", "90", u],
                       capture_output=True, text=True)
    for line in r.stdout.splitlines():
        if line.lower().startswith("content-length:"):
            return int(line.split(":")[1].strip())
    return None


def need(path, lo, hi):
    want = hi - lo + 1
    have = os.path.getsize(path) if os.path.exists(path) else 0
    return want - have


def get_chunk(args):
    u, path, lo, hi = args
    for _ in range(10):
        rem = need(path, lo, hi)
        if rem <= 0:
            return (path, True)
        have = (hi - lo + 1) - rem
        cmd = ["curl", "-sL", "--fail", "--retry", "3", "--retry-delay", "5",
               "--max-time", "7200", "-r", "%d-%d" % (lo + have, hi)]
        if have:
            cmd += ["-C", "-"]
        cmd += ["-o", path, u]
        subprocess.run(cmd)
        if need(path, lo, hi) <= 0:
            return (path, True)
        time.sleep(5)
    return (path, False)


def main():
    jobs, plan = [], []
    for rec, name in FILES:
        final = os.path.join(OUT, name)
        total = size_of(url(rec, name))
        if total is None:
            print("[ERR] size unknown:", name, flush=True)
            continue
        if os.path.exists(final) and os.path.getsize(final) == total:
            print("[skip] %s already complete" % name, flush=True)
            continue
        if os.path.exists(final):
            os.remove(final)          # 舊的部分檔案改走分塊流程
        step = total // NCHUNK + 1
        parts = []
        for i in range(NCHUNK):
            lo, hi = i * step, min(total - 1, (i + 1) * step - 1)
            if lo > hi:
                break
            p = os.path.join(PART, "%s.part%02d" % (name, i))
            parts.append((p, lo, hi))
            jobs.append((url(rec, name), p, lo, hi))
        plan.append((name, final, total, parts))
        print("[plan] %-36s %6.2f GB -> %d chunks" % (name, total / 1e9, len(parts)), flush=True)

    if not jobs:
        print("=== NOTHING TO DO ===", flush=True)
        return
    t0 = time.time()
    with cf.ThreadPoolExecutor(max_workers=len(jobs)) as ex:
        res = list(ex.map(get_chunk, jobs))
    bad = [p for p, ok in res if not ok]
    if bad:
        print("[ERR] failed chunks:", bad, flush=True)
        sys.exit(1)
    print("[info] chunks done in %.1f min" % ((time.time() - t0) / 60), flush=True)

    for name, final, total, parts in plan:
        with open(final, "wb") as fo:
            for p, lo, hi in parts:
                with open(p, "rb") as fi:
                    while True:
                        b = fi.read(1 << 24)
                        if not b:
                            break
                        fo.write(b)
        got = os.path.getsize(final)
        ok = (got == total)
        print("[done] %-36s %6.2f GB  %s" % (name, got / 1e9, "OK" if ok else "SIZE MISMATCH"), flush=True)
        if ok:
            for p, _, _ in parts:
                os.remove(p)
    print("=== ALL DOWNLOADS FINISHED ===", flush=True)


main()
