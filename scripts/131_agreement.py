#!/usr/bin/env python
"""步驟 131：步驟 130 四種方法的「一致性」統計 —— 門檻定義寫死在這裡。

稽核時發現 results/130_agreement.json 是用臨時指令產生的，門檻沒有留下記錄。
這個檔案把定義固定下來，讓它可以被重現與被質疑。

定義（唯一的一個，不要再用別的）
--------------------------------
「某基因被某方法認可」＝ 該方法的 tumour_specific 殘差排名進入**前 5%**
（4,012 個基因 → 前 200 名）。
「n/4 一致」＝ 有 n 個方法認可它。

另外報告一個不同的量，**不要和上面混用**：
  z = tumour_specific / sd(mean_nontumour across genes)
  這是以「非腫瘤殘差在基因間的標準差」為尺度的訊噪比，用來比較方法的靈敏度，
  不是用來定義一致性。SLC7A11 在前 5% 的定義下是 3/4，在 z>3 的定義下是 1/4 ——
  兩者都對，只是問的問題不同。
"""
import json
import pandas as pd

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
RES = f"{BASE}/results"
TOP_FRAC = 0.05
PANEL = ["SLC7A11", "SLC3A2", "GCLM", "GCLC", "AIFM2", "MAFG", "MAFF", "MAFK",
         "BACH1", "NFE2L2", "CCND1", "SOX11"]
TRUTH = ["CCND1", "SOX11"]


def main():
    piv = pd.read_csv(f"{RES}/130_tumour_specific_by_method.csv", index_col=0)
    R = pd.read_csv(f"{RES}/130_residual_by_method.csv.gz")
    meth = list(piv.columns)
    cut = int(round(len(piv) * TOP_FRAC))
    top = pd.DataFrame({m: piv[m].rank(ascending=False) <= cut for m in meth})
    cnt = top.sum(1)

    z = {}
    for m, d in R.groupby("method"):
        d = d.set_index("gene")
        sd = d["mean_nontumour"].std()
        z[m] = {"noise_sd": float(sd),
                "n_z_gt3": int(((d["tumour_specific"] / sd) > 3).sum()),
                "panel_z": {g: float(d.loc[g, "tumour_specific"] / sd)
                            for g in PANEL if g in d.index},
                "truth_z": {g: float(d.loc[g, "tumour_specific"] / sd)
                            for g in TRUTH if g in d.index}}

    out = {"definition": {
               "agreement": f"tumour_specific 排名進入前 {TOP_FRAC:.0%}"
                            f"（{len(piv)} 個基因取前 {cut} 名）",
               "z": "tumour_specific / sd(mean_nontumour across genes)；"
                    "與 agreement 是兩個不同的量，勿混用"},
           "n_genes": int(len(piv)), "top_n": cut, "methods": meth,
           "n_agree_4of4": int((cnt == 4).sum()),
           "n_agree_3of4": int((cnt == 3).sum()),
           "n_agree_1of4": int((cnt == 1).sum()),
           "agree_4of4_genes": list(piv[cnt == 4]
                                    .sort_values(meth[0], ascending=False).index),
           "panel_agreement": {g: int(cnt[g]) for g in PANEL if g in cnt.index},
           "panel_rank_per_method": {
               g: {m: int(piv[m].rank(ascending=False)[g]) for m in meth}
               for g in PANEL if g in piv.index},
           "sensitivity_z": z}
    json.dump(out, open(f"{RES}/130_agreement.json", "w"), indent=2,
              ensure_ascii=False)
    print("[131] 門檻：前 %d 名 / %d 個基因" % (cut, len(piv)))
    print("[131] 4/4 一致 %d 個、3/4 %d 個、僅 1 法 %d 個"
          % (out["n_agree_4of4"], out["n_agree_3of4"], out["n_agree_1of4"]))
    print("[131] panel:", out["panel_agreement"])
    print("[131] DONE → results/130_agreement.json（門檻已寫入檔案）")


if __name__ == "__main__":
    main()
