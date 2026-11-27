#!/usr/bin/env python
"""步驟 235b：只做一件事 —— 把序列用 DNA 語言模型嵌成向量。

為什麼單獨切出來
----------------
嵌入需要 torch ≥ 2.6（transformers 4.57 因 CVE-2025-32434 拒絕用舊版 torch 載入
`.bin` 權重，而 Nucleotide Transformer 的所有 repo 都只提供 `.bin`、沒有 safetensors）。
符合條件的環境是 `env/esmfold`（torch 2.7.1+cu126），但它沒有 pandas／sklearn／pyfaidx；
而有那些套件的 `env/tonsilmcl` 只有 torch 2.4.1。
**刻意不去動任何一個能用的環境**（裝 pandas 可能牽動 numpy 版本，那是會弄壞既有環境的那種改動），
改成把流程切成三段，各自在已經具備條件的環境裡跑：
    235（env/tonsilmcl）  取序列 → 235b（env/esmfold）嵌入 → 235（env/tonsilmcl）檢定
本檔只依賴 numpy + torch + transformers。

用法
----
    env/esmfold/bin/python -I scripts/235b_embed.py <輸入.txt> <輸出.npy> [batch]

輸入是一行一條序列的純文字檔（已經裁切成固定長度），輸出是 float16 的 (n, hidden) 矩陣。
**不傳 trust_remote_code**：v1 的 config 沒有 auto_map，用 transformers 內建的
EsmForMaskedLM 就能載入，不執行任何從 HuggingFace 下載的程式碼。
"""
import os
import sys

os.environ["HF_HOME"] = "/var/tmp/tonsilmcl_hf"
for _v in ("REQUESTS_CA_BUNDLE", "CURL_CA_BUNDLE", "SSL_CERT_FILE"):
    os.environ[_v] = "/etc/ssl/certs/ca-certificates.crt"

import numpy as np

MODEL = "InstaDeepAI/nucleotide-transformer-500m-human-ref"


def log(*a):
    print("[235b]", *a, flush=True)


def main():
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    src, dst = sys.argv[1], sys.argv[2]
    batch = int(sys.argv[3]) if len(sys.argv) > 3 else 16

    seqs = [l.strip() for l in open(src) if l.strip()]
    log("讀入 %d 條序列，長度中位 %d" % (len(seqs), int(np.median([len(s) for s in seqs]))))

    import torch
    from transformers import AutoTokenizer, EsmConfig, EsmForMaskedLM
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    log("torch %s | device %s" % (torch.__version__, dev))
    tok = AutoTokenizer.from_pretrained(MODEL)
    cfg = EsmConfig.from_pretrained(MODEL)
    mod = EsmForMaskedLM.from_pretrained(MODEL)
    npar = sum(p.numel() for p in mod.parameters())
    log("參數 %.0f M | hidden %d | 層 %d" % (npar / 1e6, cfg.hidden_size,
                                             cfg.num_hidden_layers))
    if npar < 100e6:
        raise SystemExit("參數量不符預期，權重可能沒載進來，中止")
    mod = mod.to(dev).eval()

    out = []
    for a in range(0, len(seqs), batch):
        blk = seqs[a:a + batch]
        enc = tok(blk, return_tensors="pt", padding=True, truncation=True)
        enc = {k: v.to(dev) for k, v in enc.items()}
        with torch.no_grad():
            h = mod(**enc, output_hidden_states=True).hidden_states[-1]
            m = enc["attention_mask"].unsqueeze(-1).to(h.dtype)
            v = (h * m).sum(1) / m.sum(1).clamp(min=1)
        out.append(v.float().cpu().numpy())
        if (a // batch) % 25 == 0:
            log("  %d / %d" % (a, len(seqs)))
    Z = np.vstack(out).astype(np.float16)
    np.save(dst, Z)
    log("寫出 %s  shape %s" % (dst, Z.shape))
    return 0


if __name__ == "__main__":
    sys.exit(main())
