#!/usr/bin/env python
"""步驟 218（階段一）：BLUEPRINT 的 5 位原發 MCL 病人能不能用？

為什麼這是今天最有價值的剩餘一步
--------------------------------
這個專案最大的弱點一直是 **N = 2**（公開資料裡沒有第三例病人來源的 MCL 多體學）。
外部資料偵查（`results/206_external_validation_scouting.md`）找到 BLUEPRINT
在 GRCh38 上有 **5 位原發 MCL 病人**（007MCL / 2392MCL / 6902MCL / 817-03MCL / 8761MCL），
每人 6 種組蛋白標記 + DNA 甲基化 + total-RNA-seq，配對 3 位 class-switched 記憶 B
加 1 位 unswitched，**同一條流程、開放下載、不需申請**。
如果可用，N 就從 2 變成 7。

但在比較任何 peak 之前，有一個把關必須先過
------------------------------------------
這正是今天在 DepMap 上踩到的陷阱：**三株 MCL 細胞株的 `SLC7A11` 都低於中位**
（第 44、17、12 百分位），而病人腫瘤高於配對正常 B 細胞 200 倍 ——
在那種細胞裡做鐵死亡實驗會得到一個無意義的陰性。
同一個邏輯適用於這裡：**如果這 5 位病人不表現 `SLC7A11`，
那他們在那五個元件上的任何陰性結果都無法解讀。**

事前指定（階段一只做把關，不碰 peak）
------------------------------------
G1  至少 3 / 5 位 MCL 病人的 `SLC7A11` 必須明確表現（FPKM > 1）。
    不過就停止，並把這件事當成結論寫出來（「這個外部世代測不了這條軸線」）。
G2  陽性對照：`CCND1` 必須在 5/5 位病人高表現（t(11;14) 是 MCL 的定義性易位）。
    不過就代表我抓錯樣本或檔案，整個資源不可信。
G3  陰性對照：`CD3E`（T 細胞標記）應該低 —— 確認樣本不是被 T 細胞污染的。

階段二（只有 G1–G3 都過才做）會下載 H3K27ac / H3K4me1 的 peak call，
而且必須先算全基因組的 peak 基礎率再比較那五個元件 ——
今天的 caQTL 分析就是因為沒先算基礎率，差點把「0 個命中」當成陰性結果
（實際期望命中數只有 0.125）。

⚠️ 兩個要帶著走的限制：H3K27ac 是**活性**不是可及性；
   BLUEPRINT 的 MCL 樣本來自**靜脈血**，本專案的是組織。
⚠️ 與 M102 / M413 的捐贈者是否重疊**未經確認**，寫「獨立」之前必須先查。

輸出 results/218_*
"""
import os
import json
import subprocess

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "8"

import pandas as pd

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
RES = f"{BASE}/results"
# 下載的暫存目錄：可用 TONSILMCL_SCRATCH 覆寫，預設放在專案外的 /var/tmp
SCRATCH = os.environ.get("TONSILMCL_SCRATCH", "/var/tmp/tonsilmcl_scratch")
DL = f"{SCRATCH}/bp_rna"
IDX = f"{SCRATCH}/bp_index/idx.tsv"
FTP = "https://ftp.ebi.ac.uk/pub/databases/"
MCL = ["007MCL", "2392MCL", "6902MCL", "817-03MCL", "8761MCL"]
GENES = ["SLC7A11", "CCND1", "SOX11", "MAFG", "SLC3A2", "GCLC", "GCLM",
         "AIFM2", "GPX4", "SCLY", "GSS", "NAMPT", "CD3E", "CD19", "MS4A1"]


def log(*a):
    print("[218]", *a, flush=True)


def main():
    os.makedirs(DL, exist_ok=True)
    if not os.path.exists(IDX):
        raise SystemExit(
            "缺少 BLUEPRINT 的檔案索引：%s\n"
            "這個索引不在本倉庫裡（它是 EBI 的清單檔，不是我們的產出）。\n"
            "取得方式：到 %sblueprint/ 的 index 目錄下載清單檔，\n"
            "需要的欄位是 DONOR_ID、FILE_TYPE、FILE、FILE_SIZE，存成上面那個路徑。\n"
            "已有一份的話，設 TONSILMCL_SCRATCH 指向它的上層目錄即可。" % (IDX, FTP))
    d = pd.read_csv(IDX, sep="\t", low_memory=False)
    q = d[(d.DONOR_ID.isin(MCL)) & (d.FILE_TYPE == "RNA_GENE_QUANT_STAR_CRG")]
    log("索引裡 MCL 的基因定量檔 %d 個" % len(q))
    # 每位捐贈者取一個檔
    one = q.sort_values("FILE").groupby("DONOR_ID").first().reset_index()
    log("每位捐贈者取一個檔，共 %d 位" % len(one))

    rows, got = [], {}
    for _, r in one.iterrows():
        url = FTP + r.FILE
        dest = f"{DL}/{r.DONOR_ID}.results"
        if not os.path.exists(dest) or os.path.getsize(dest) < 1e6:
            log("下載 %s（%.1f MB）…" % (r.DONOR_ID, float(r.FILE_SIZE) / 1e6))
            p = subprocess.run(["curl", "-sS", "--max-time", "600",
                                "-o", dest, url], capture_output=True, text=True)
            if p.returncode != 0:
                log("  失敗：%s" % p.stderr[-200:]); continue
        try:
            t = pd.read_csv(dest, sep="\t")
        except Exception as e:
            log("  讀取失敗 %s：%s" % (r.DONOR_ID, e)); continue
        got[r.DONOR_ID] = t
        log("  %s：%d 列，欄位 %s" % (r.DONOR_ID, len(t), list(t.columns)[:6]))

    if not got:
        raise SystemExit("沒有成功取得任何檔案")

    # RSEM 的 gene_id 通常是 ENSG，需要用符號對照；先看有沒有 gene_name 欄
    sample = next(iter(got.values()))
    idcol = "gene_id" if "gene_id" in sample.columns else sample.columns[0]
    valcol = ("FPKM" if "FPKM" in sample.columns
              else ("TPM" if "TPM" in sample.columns else sample.columns[-1]))
    log("\n用 %s 當識別欄、%s 當表現量" % (idcol, valcol))

    # 從 GENCODE 取 ENSG ↔ 符號
    gtf = f"{BASE}/data/genome/gencode.v44.annotation.gtf.gz"
    import gzip
    import re
    want = set(GENES)
    sym2ens = {}
    with gzip.open(gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if f[2] != "gene":
                continue
            m = re.search(r'gene_name "([^"]+)"', f[8])
            if m and m.group(1) in want:
                e = re.search(r'gene_id "([^"]+)"', f[8]).group(1)
                sym2ens[m.group(1)] = e.split(".")[0]
    log("對到 %d / %d 個基因的 ENSG" % (len(sym2ens), len(GENES)))

    for donor, t in got.items():
        ids = t[idcol].astype(str).str.split(".").str[0]
        t = t.assign(_ens=ids)
        for g in GENES:
            e = sym2ens.get(g)
            if e is None:
                continue
            m = t[t._ens == e]
            rows.append(dict(donor=donor, gene=g,
                             value=float(m[valcol].iloc[0]) if len(m) else float("nan")))
    D = pd.DataFrame(rows)
    P = D.pivot(index="gene", columns="donor", values="value").reindex(GENES)
    P.to_csv(f"{RES}/218_blueprint_mcl_expression.csv")
    log("\n=== 5 位原發 MCL 病人的表現量（%s）===" % valcol)
    log("\n" + P.round(2).to_string())

    out = {"n_donors": int(len(got)), "value_column": valcol,
           "donors": sorted(got), "expression": json.loads(P.round(3).to_json())}
    sl = P.loc["SLC7A11"].dropna()
    out["G1_SLC7A11_expressed"] = {
        "n_above_1": int((sl > 1).sum()), "n_donors": int(len(sl)),
        "values": {k: float(v) for k, v in sl.items()},
        "pass": bool((sl > 1).sum() >= 3)}
    cc = P.loc["CCND1"].dropna()
    out["G2_CCND1_positive_control"] = {
        "n_above_10": int((cc > 10).sum()), "n_donors": int(len(cc)),
        "values": {k: float(v) for k, v in cc.items()},
        "pass": bool((cc > 10).sum() == len(cc) and len(cc) > 0)}
    cd = P.loc["CD3E"].dropna()
    out["G3_CD3E_negative_control"] = {
        "median": float(cd.median()) if len(cd) else None,
        "pass": bool(len(cd) > 0 and cd.median() < 5)}
    for k in ("G1_SLC7A11_expressed", "G2_CCND1_positive_control",
              "G3_CD3E_negative_control"):
        log("%s → %s" % (k, "過" if out[k]["pass"] else "不過"))
    out["stage2_authorised"] = bool(all(out[k]["pass"] for k in (
        "G1_SLC7A11_expressed", "G2_CCND1_positive_control",
        "G3_CD3E_negative_control")))
    log("\n階段二（下載 peak call）是否放行：%s" % out["stage2_authorised"])
    json.dump(out, open(f"{RES}/218_blueprint_gate.json", "w"), indent=2,
              ensure_ascii=False)
    log("DONE → results/218_*")


if __name__ == "__main__":
    main()
