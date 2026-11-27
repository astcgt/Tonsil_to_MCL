#!/usr/bin/env python
"""步驟 159：讓 ESMFold 在這台機器上跑起來（可重現、可還原的環境修補）。

本機有三個障礙，都和模型本身無關：
  1. openfold 不在環境裡，而且它的 setup.py 要編 CUDA 擴充 —— 本機沒有 nvcc。
     → 直接取原始碼（third_party/openfold），不安裝、不編譯。
  2. openfold/utils/kernel/attention_core.py 在 import 時就要載入編譯好的
     attn_core_inplace_cuda。ESMFold 的預設推論路徑不會呼叫這個 kernel。
     → third_party/openfold_shim/attn_core_inplace_cuda.py 放一個替身，
       被呼叫時會立刻拋錯（不會悄悄算出錯的數值）。
  3. fair-esm 的 trunk.py 有一個 Python 3.11 不接受的可變 dataclass 預設值。
     → 本腳本把那一行改成 field(default_factory=...)，並留下 .orig 備份。

執行一次即可。之後跑 ESMFold 要帶：
  PYTHONPATH=third_party/openfold:third_party/openfold_shim env/genomic/bin/python ...
還原：把 .orig 備份覆蓋回去。
"""
import os
import re
import shutil
import sys

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
V1 = f"{BASE}/env/genomic/lib/python3.11/site-packages/esm/esmfold/v1"
FILES = ["trunk.py", "esmfold.py"]
# 任何 "name: SomeConfig = SomeConfig()" 都要改成 default_factory
PAT = re.compile(r"^(\s*)(\w+): ([\w.\[\], ]+?) = (\w+)\(\)\s*$", re.M)


def log(*a):
    print("[159]", *a, flush=True)


def main():
    for fn in FILES:
        path = os.path.join(V1, fn)
        src = open(path).read()
        hits = PAT.findall(src)
        if not hits:
            log("%s：沒有可變預設值需要處理" % fn); continue
        if not os.path.exists(path + ".orig"):
            shutil.copy(path, path + ".orig")
            log("已備份 →", path + ".orig")
        if "from dataclasses import dataclass, field" not in src:
            src = src.replace("from dataclasses import dataclass",
                              "from dataclasses import dataclass, field", 1)
        src = PAT.sub(
            lambda m: ("%s%s: %s = field(default_factory=%s)"
                       "  # patched: Python 3.11 不接受可變 dataclass 預設值"
                       % (m.group(1), m.group(2), m.group(3), m.group(4))), src)
        open(path, "w").write(src)
        log("%s：修補了 %d 處 %s" % (fn, len(hits), [h[3] for h in hits]))

    for p in (f"{BASE}/third_party/openfold/openfold",
              f"{BASE}/third_party/openfold_shim/attn_core_inplace_cuda.py"):
        log(("有" if os.path.exists(p) else "缺") + ":", p)

    log("\n驗證 import …")
    os.environ["PYTHONPATH"] = (f"{BASE}/third_party/openfold:"
                                f"{BASE}/third_party/openfold_shim")
    import subprocess
    r = subprocess.run([f"{BASE}/env/genomic/bin/python", "-c",
                        "import torch; from esm.esmfold.v1.pretrained import esmfold_v1;"
                        " print('ESMFold import OK | torch', torch.__version__,"
                        " '| CUDA', torch.cuda.is_available())"],
                       env=dict(os.environ), capture_output=True, text=True)
    print(r.stdout.strip() or r.stderr.strip()[-800:])
    sys.exit(0 if "ESMFold import OK" in r.stdout else 1)


if __name__ == "__main__":
    main()
