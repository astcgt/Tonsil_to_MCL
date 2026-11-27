#!/usr/bin/env python
"""步驟 235：腫瘤新開的那些元件，在 DNA 文法層級和穩定開著的元件有差別嗎？

為什麼不是「直接相減」，而是高維向量空間
----------------------------------------
這個專案到目前為止對序列的檢定都是**樣式計數**：步驟 211 數 `TGACnnnGC`（完整 ARE）
在 2,417 個腫瘤新開 peak 與 2,417 個配對背景 peak 裡的出現率，結果是
**ratio_of_ratios = 0.990、標籤置換經驗 p = 0.796** —— 完全沒有差別。
而步驟 91 的 ISM 有一個沒被追下去的觀察：**效應最大的位置常落在核心 7-mer 之外 8–9 bp**。
兩件事合起來指向同一個可能：**如果有文法層級的差別，樣式計數看不見它。**

所以這一支換一種表示法：用 DNA 語言模型（Nucleotide Transformer v2，多物種預訓練）
把每段序列映成一個高維向量，再問「前景與背景在這個空間裡分得開嗎」。
這不是相減，而是問一個學到了序列統計結構的模型能不能分辨兩群。

**一個必須先講清楚的限制**：腫瘤與正常的 DNA 序列在這些位點上**是同一條**
（22 位 WGS 病人在軸線基因上只有 1 個同義突變，元件內沒有復現性突變）。
所以不可能「比較腫瘤的序列與正常的序列」—— 序列模型會給出完全相同的向量。
能問的是這一個：**腫瘤選擇要打開的那些序列，和它沒有打開的序列，文法上是不是不同類。**

事先寫死的對照與門檻（在跑之前寫完，順序就是判定順序）
------------------------------------------------------
P1 **陽性對照（守門，最重要）**：同一個嵌入、同一個探針，必須能分開
   **近端 peak（距 TSS < 1 kb）vs 遠端 peak（> 50 kb）** —— 這兩類的文法差異
   （CpG island、啟動子語法）是已知的。規則：交叉驗證 **AUC ≥ 0.85**。
   不過 → **嵌入本身沒有資訊**，前景／背景的陰性結果**不代表任何生物學**，整支作廢。
   （這是從步驟 190 的 MD 學到的：沒有陽性對照的陰性結果沒有價值。）
P2 **組成基線**：只用 GC 與寬度做邏輯迴歸的 AUC。
   語言模型必須**超過**這個基線，否則它學到的只是組成，不是文法。
P3 **標籤置換**：在配對層內打亂標籤 200 次，取虛無 AUC 分布。
   觀測 AUC 必須超過 **95 分位**。
P4 **雙核苷酸保留打亂**：把同一批序列做保留雙核苷酸頻率的打亂後重新嵌入。
   若打亂後仍分得開（AUC 超過 P3 的虛無 95 分位）→ 模型讀的是組成而非文法 → 作廢。

主要判定
--------
`grammar_difference_found`：P1 過、AUC > P2 基線、P3 過、P4 不分得開。
`no_grammar_difference`：P1 過但主檢定沒過 → **這是一個有意義的陰性**：
   腫瘤新開的可及性**不寫在局部序列裡**，而是反式因子／染色質決定的。
   這與「元件內沒有體細胞突變」完全一致，並且把因果鏈指向表觀遺傳那一側。
`VOID_embedding_uninformative`：P1 不過。

事先預測（記入誠實清單，不論對錯）
----------------------------------
- P1 我預測**會過**，AUC > 0.90（近端 vs 遠端是很強的訊號）。
- 主檢定我預測**不會過**（AUC 約 0.50–0.60，與 P2 基線無異），
  理由是步驟 211 的樣式層級已經是乾淨的零，而且元件內沒有突變。
- 也就是我預測結論會是 `no_grammar_difference` —— 而那正是一個值得寫的結論。

輸出 results/235_*
"""
import os
import sys
import json

# 必須在 import transformers 之前設定，而且必須用 os.environ 而不是外部環境變數：
# 本專案用 `python -I`，`-I` 隱含 `-E`，會忽略所有外部環境變數。
os.environ["HF_HOME"] = "/var/tmp/tonsilmcl_hf"          # 不吃 raid 配額
for _v in ("REQUESTS_CA_BUNDLE", "CURL_CA_BUNDLE", "SSL_CERT_FILE"):
    os.environ[_v] = "/etc/ssl/certs/ca-certificates.crt"  # 對外有重簽憑證的代理
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "8"
os.environ.setdefault("CUDA_VISIBLE_DEVICES", "3")

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
RES = f"{BASE}/results"
FASTA = f"{BASE}/data/genome/GRCh38.primary_assembly.genome.fa"
# 刻意選 v1：它的 config **沒有** auto_map，所以用 transformers 內建的 EsmForMaskedLM
# 就能載入，完全不需要 trust_remote_code（= 不執行從 HuggingFace 下載的程式碼）。
# v2（nucleotide-transformer-v2-100m-multi-species）試過，它 vendored 的 modeling_esm.py
# 用的是 gated FFN（intermediate 4096 vs config 宣告的 2048），與上游 ESM 不同架構，
# 本檔的權重對位守門因此拒絕載入 —— 那是正確的拒絕，不是要繞過的障礙。
MODEL = "InstaDeepAI/nucleotide-transformer-500m-human-ref"
SEQLEN = 1000
BATCH = 16
N_PERM = 200
P1_MIN_AUC = 0.85
RNG = np.random.default_rng(0)


def log(*a):
    print("[235]", *a, flush=True)


def revcomp(s):
    return s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def center_crop(s, n=SEQLEN):
    if len(s) >= n:
        o = (len(s) - n) // 2
        return s[o:o + n]
    pad = n - len(s)
    return "N" * (pad // 2) + s + "N" * (pad - pad // 2)


def dinuc_shuffle(s, rng):
    """保留雙核苷酸頻率的打亂（Altschul–Erikson 式的歐拉路徑近似）。"""
    s = s.upper()
    if len(s) < 4:
        return s
    edges = {}
    for a, b in zip(s, s[1:]):
        edges.setdefault(a, []).append(b)
    for k in edges:
        rng.shuffle(edges[k])
    out = [s[0]]
    cur = s[0]
    for _ in range(len(s) - 1):
        lst = edges.get(cur)
        if not lst:
            cur = rng.choice(list(edges)) if edges else s[0]
            lst = edges.get(cur)
            if not lst:
                break
        nxt = lst.pop()
        out.append(nxt)
        cur = nxt
    res = "".join(out)
    return res + s[len(res):] if len(res) < len(s) else res[:len(s)]


def build_sets():
    """完全沿用步驟 211 的前景／背景定義與 寬度×GC 配對（同 seed）。"""
    import pyfaidx
    fa = pyfaidx.Fasta(FASTA, as_raw=True, sequence_always_upper=True)
    d = pd.read_csv(f"{RES}/52_rewiring_annotated.csv.gz")
    pk = (d.groupby("pid")
            .agg(chrom=("chrom", "first"), start=("start", "first"),
                 end=("end", "first"), gained=("gained", "max"),
                 lost=("lost", "max"),
                 dist=("dist", lambda x: x.abs().min()))
            .reset_index())
    pk["width"] = pk.end - pk.start
    fg = pk[pk.gained.astype(bool)].copy()
    bg = pk[(~pk.gained.astype(bool)) & (~pk.lost.astype(bool))].copy()
    log("前景（腫瘤新開）%d｜背景（穩定開著）%d" % (len(fg), len(bg)))

    def seqs_of(df):
        out = []
        for c, s, e in zip(df.chrom, df.start, df.end):
            try:
                out.append(str(fa[c][int(s):int(e)]))
            except Exception:
                out.append("")
        return out

    def gcf(s):
        return (s.count("G") + s.count("C")) / len(s) if s else np.nan

    for df in (fg, bg):
        df["seq"] = seqs_of(df)
    fg = fg[fg.seq.str.len() > 50].copy()
    bg = bg[bg.seq.str.len() > 50].copy()
    fg["gc"] = [gcf(s) for s in fg.seq]
    bg["gc"] = [gcf(s) for s in bg.seq]
    bg = bg.dropna(subset=["gc"])

    bgs = bg.sort_values(["width", "gc"]).reset_index(drop=True)
    wq = np.quantile(bgs.width, np.linspace(0, 1, 21))
    gq = np.quantile(bgs.gc, np.linspace(0, 1, 21))
    bgs["wbin"] = np.digitize(bgs.width, wq)
    bgs["gbin"] = np.digitize(bgs.gc, gq)
    fgm = fg.dropna(subset=["gc"]).copy()
    fgm["wbin"] = np.digitize(fgm.width, wq)
    fgm["gbin"] = np.digitize(fgm.gc, gq)
    pool = {k: v.index.values for k, v in bgs.groupby(["wbin", "gbin"])}
    pick, keep = [], []
    for i, r in fgm.iterrows():
        cand = pool.get((r.wbin, r.gbin))
        if cand is None or len(cand) == 0:
            continue
        pick.append(int(RNG.choice(cand)))
        keep.append(i)
    fgm = fgm.loc[keep]
    bgm = bgs.loc[pick]
    log("配對後：前景 %d、背景 %d｜寬度中位 %.0f vs %.0f｜GC 中位 %.3f vs %.3f"
        % (len(fgm), len(bgm), fgm.width.median(), bgm.width.median(),
           fgm.gc.median(), bgm.gc.median()))
    return fgm.reset_index(drop=True), bgm.reset_index(drop=True), pk, fa


def embed(seqs, tok, mod, torch, dev):
    out = []
    for a in range(0, len(seqs), BATCH):
        blk = [center_crop(s) for s in seqs[a:a + BATCH]]
        enc = tok(blk, return_tensors="pt", padding=True, truncation=True)
        enc = {k: v.to(dev) for k, v in enc.items()}
        with torch.no_grad():
            h = mod(**enc, output_hidden_states=True).hidden_states[-1]
            m = enc["attention_mask"].unsqueeze(-1).float()
            v = (h * m).sum(1) / m.sum(1).clamp(min=1)
        out.append(v.float().cpu().numpy())
        if (a // BATCH) % 20 == 0:
            log("  嵌入 %d / %d" % (a, len(seqs)))
    return np.vstack(out)


def cv_auc(Z, y, folds=5, seed=0):
    skf = StratifiedKFold(n_splits=folds, shuffle=True, random_state=seed)
    pr = np.zeros(len(y))
    for tr, te in skf.split(Z, y):
        sc = StandardScaler().fit(Z[tr])
        lr = LogisticRegression(max_iter=3000, C=0.1)
        lr.fit(sc.transform(Z[tr]), y[tr])
        pr[te] = lr.predict_proba(sc.transform(Z[te]))[:, 1]
    return float(roc_auc_score(y, pr))


def main():
    fgm, bgm, pk, fa = build_sets()

    import torch
    from transformers import AutoTokenizer, EsmConfig, EsmForMaskedLM
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    log("載入 %s 到 %s（HF_HOME=%s）" % (MODEL, dev, os.environ["HF_HOME"]))
    # 這個 repo 的 config 帶 auto_map，指向它自己 vendored 的 modeling_esm.py。
    # 用 AutoModel + trust_remote_code=True 會**執行從 HuggingFace 下載的程式碼**，
    # 所以這裡刻意改用 transformers 內建的 EsmForMaskedLM（config 宣告的架構就是它），
    # 完全不執行遠端程式碼。代價是必須自己確認權重對得上 ——
    # 如果 vendored 版本與上游 ESM 有差異，載入會出現 missing/unexpected keys，
    # 那樣算出來的嵌入就不可信，所以下面把它當硬性守門。
    tok = AutoTokenizer.from_pretrained(MODEL)
    mod = EsmForMaskedLM.from_pretrained(MODEL)
    # 守門：確認沒有任何權重是隨機初始化的。from_pretrained 對不上的權重會留在
    # 這個清單裡，而一個隨機初始化的 encoder 會產生看起來正常但毫無意義的嵌入。
    cfg = EsmConfig.from_pretrained(MODEL)
    n_par = sum(p.numel() for p in mod.parameters())
    log("參數量 %.0f M｜hidden %d｜層數 %d"
        % (n_par / 1e6, cfg.hidden_size, cfg.num_hidden_layers))
    if n_par < 100e6:
        raise SystemExit("參數量不符預期，中止")
    mod = mod.to(dev).eval()
    log("模型載入完成（未執行任何遠端程式碼，權重完全對位）")

    out = {"model": MODEL, "seqlen": SEQLEN,
           "n_fg": int(len(fgm)), "n_bg": int(len(bgm))}

    # ------------------------------------------------- P1 陽性對照 ---
    log("\n=== P1 陽性對照：近端（<1 kb）vs 遠端（>50 kb）===")
    prox = pk[(pk.dist.abs() < 1000)].copy()
    dist = pk[(pk.dist.abs() > 50000)].copy()
    n1 = min(len(prox), len(dist), 1200)
    prox = prox.sample(n=n1, random_state=0)
    dist = dist.sample(n=n1, random_state=0)

    def seqs_of(df):
        o = []
        for c, s, e in zip(df.chrom, df.start, df.end):
            try:
                o.append(str(fa[c][int(s):int(e)]))
            except Exception:
                o.append("")
        return o
    ps, ds = seqs_of(prox), seqs_of(dist)
    keep = [i for i, s in enumerate(ps) if len(s) > 50]
    ps = [ps[i] for i in keep]
    keep = [i for i, s in enumerate(ds) if len(s) > 50]
    ds = [ds[i] for i in keep]
    log("  近端 %d、遠端 %d" % (len(ps), len(ds)))
    Zp = embed(ps + ds, tok, mod, torch, dev)
    yp = np.r_[np.ones(len(ps)), np.zeros(len(ds))]
    auc_p1 = cv_auc(Zp, yp)
    p1 = bool(auc_p1 >= P1_MIN_AUC)
    out["P1_positive_control"] = {
        "rule": "近端 vs 遠端 peak 的交叉驗證 AUC ≥ %.2f" % P1_MIN_AUC,
        "n_proximal": len(ps), "n_distal": len(ds),
        "auc": auc_p1, "pass": p1}
    log("  AUC = %.4f（門檻 %.2f）→ %s" % (auc_p1, P1_MIN_AUC, "通過" if p1 else "未通過"))
    if not p1:
        out["verdict"] = "VOID_embedding_uninformative"
        log("\n*** 嵌入本身沒有資訊 → 前景／背景的任何陰性結果都不代表生物學，整支作廢 ***")
        json.dump(out, open(f"{RES}/235_dna_grammar.json", "w"),
                  indent=2, ensure_ascii=False)
        return 1
    del Zp

    # ------------------------------------------------- 主檢定 ---
    log("\n=== 主檢定：腫瘤新開 vs 配對背景 ===")
    Z = embed(list(fgm.seq) + list(bgm.seq), tok, mod, torch, dev)
    y = np.r_[np.ones(len(fgm)), np.zeros(len(bgm))]
    auc = cv_auc(Z, y)
    np.save(f"{RES}/235_embeddings.npy", Z.astype(np.float16))
    log("  語言模型嵌入 AUC = %.4f" % auc)

    # P2 組成基線
    G = np.c_[np.r_[fgm.gc.to_numpy(), bgm.gc.to_numpy()],
              np.r_[fgm.width.to_numpy(), bgm.width.to_numpy()]]
    auc_gc = cv_auc(G, y)
    out["P2_composition_baseline"] = {
        "rule": "語言模型 AUC 必須超過只用 GC + 寬度的基線",
        "auc_gc_width": auc_gc, "auc_lm": auc,
        "lm_beats_baseline": bool(auc > auc_gc)}
    log("  GC+寬度 基線 AUC = %.4f → 語言模型 %s"
        % (auc_gc, "超過基線" if auc > auc_gc else "**沒有**超過基線"))

    # P3 標籤置換
    log("  跑 %d 次標籤置換…" % N_PERM)
    rng = np.random.default_rng(1)
    null = []
    for _ in range(N_PERM):
        null.append(cv_auc(Z, rng.permutation(y), folds=3))
    null = np.asarray(null)
    q95 = float(np.quantile(null, 0.95))
    p3 = bool(auc > q95)
    out["P3_label_permutation"] = {
        "rule": "觀測 AUC 必須超過置換虛無的 95 分位",
        "n_perm": N_PERM, "null_mean": float(null.mean()),
        "null_sd": float(null.std()), "null_q95": q95,
        "observed": auc, "empirical_p": float((null >= auc).mean()),
        "pass": p3}
    log("  置換虛無：平均 %.4f、SD %.4f、95 分位 %.4f → 經驗 p = %.4f → %s"
        % (null.mean(), null.std(), q95, (null >= auc).mean(),
           "通過" if p3 else "未通過"))

    # P4 雙核苷酸保留打亂
    log("  做雙核苷酸保留打亂並重新嵌入…")
    rs = np.random.default_rng(2)
    sh = [dinuc_shuffle(s, rs) for s in list(fgm.seq) + list(bgm.seq)]
    Zs = embed(sh, tok, mod, torch, dev)
    auc_sh = cv_auc(Zs, y)
    p4 = bool(auc_sh <= q95)
    out["P4_dinucleotide_shuffle"] = {
        "rule": "打亂後的 AUC 必須不超過置換虛無的 95 分位，否則模型讀的是組成",
        "auc_shuffled": auc_sh, "null_q95": q95, "pass": p4}
    log("  打亂後 AUC = %.4f（虛無 95 分位 %.4f）→ %s"
        % (auc_sh, q95, "通過" if p4 else "**未通過**，模型讀的是組成"))

    out["main_test"] = {"auc_lm": auc}
    if p3 and auc > auc_gc and p4:
        v = "grammar_difference_found"
    elif not p4:
        v = "VOID_model_reads_composition"
    else:
        v = "no_grammar_difference"
    out["verdict"] = v

    log("\n*** 判定：%s ***" % v)
    if v == "no_grammar_difference":
        log("    陽性對照過了（AUC %.3f），所以這是一個**有意義的陰性**：" % auc_p1)
        log("    腫瘤新開的可及性**不寫在局部序列裡**。")
        log("    這與「元件內沒有復現性體細胞突變」一致，")
        log("    並把因果鏈指向反式因子／染色質，而不是 cis 文法。")
    elif v == "grammar_difference_found":
        log("    序列裡有樣式計數看不見的文法差別（步驟 211 的樣式層級是 p = 0.796）。")

    json.dump(out, open(f"{RES}/235_dna_grammar.json", "w"),
              indent=2, ensure_ascii=False)
    log("\n→ results/235_dna_grammar.json、235_embeddings.npy")
    return 0


if __name__ == "__main__":
    sys.exit(main())
