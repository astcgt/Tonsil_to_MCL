#!/usr/bin/env python
"""
步驟 06：把 MCL 細胞投影到「正常 B 細胞分化流形」上，並計算扣除正常分化後的殘差。

核心邏輯
--------
1. 只用**正常** B 細胞訓練 scVI → 潛在空間 Z 只編碼「正常分化變異 + 捐贈者批次」。
2. 用 scArches surgery 把 MCL 細胞投進 Z（編碼器凍結，只學新批次向量）。
   因為編碼器沒看過腫瘤變異，MCL 細胞會被投到「最像它的正常狀態」位置。
3. 對每個 MCL 細胞取 k 個最近的正常鄰居 → 期望表現 μ（= 若它是正常細胞該長的樣子）。
4. 殘差 r = 觀測 − 期望。這就是「扣掉正常分化後剩下的腫瘤訊號」。
5. 另算 novelty score = 到參考流形的距離，標記「正常程式解釋不了」的細胞。

為什麼不直接做 MCL vs 正常 B 的差異表現？
  因為 MCL 細胞的組成偏向某個分化窗口（論文說是 activated NBC → GC-committed），
  直接比會把「分化階段差異」誤判成「腫瘤異常」。本步驟就是要消掉這個混淆。
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _compat  # noqa  (scanpy/scipy 稀疏陣列相容)
import os, sys, json, argparse, warnings
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import scipy.sparse as sp
from sklearn.neighbors import NearestNeighbors

warnings.filterwarnings("ignore")
BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
PROC = f"{BASE}/data/processed"
RES = f"{BASE}/results"
os.makedirs(RES, exist_ok=True)
sc.settings.n_jobs = 32


def log(*a):
    print("[06]", *a, flush=True)


IG_PREFIX = ("IGHV", "IGKV", "IGLV", "IGHJ", "IGKJ", "IGLJ", "IGHD", "IGLC", "IGKC", "IGHA",
             "IGHG", "IGHM", "IGHE", "TRAV", "TRBV", "TRGV", "TRDV")


def gene_blocklist(var_names):
    """做映射時要排除的基因：免疫球蛋白可變區（腫瘤是單株的，會主導訊號）、
    粒線體、核糖體。注意：這些只在『建立潛在空間/HVG』時排除，殘差仍會算全基因。"""
    bl = set()
    for g in var_names:
        gu = g.upper()
        if gu.startswith(IG_PREFIX):
            bl.add(g)
        if gu.startswith("MT-") or gu.startswith("MT."):
            bl.add(g)
        if gu.startswith("RPL") or gu.startswith("RPS"):
            bl.add(g)
    return bl


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", default=f"{PROC}/ref_normalB_multiome.h5ad")
    ap.add_argument("--query", default=f"{PROC}/query_mcl.h5ad")
    ap.add_argument("--ref_batch", default="donor_id")
    ap.add_argument("--query_batch", default="sample")
    ap.add_argument("--annot", default="annotation")
    ap.add_argument("--n_latent", type=int, default=30)
    ap.add_argument("--n_hvg", type=int, default=3000)
    ap.add_argument("--k", type=int, default=50)
    ap.add_argument("--epochs_ref", type=int, default=200)
    ap.add_argument("--epochs_query", type=int, default=150)
    ap.add_argument("--gpu", default="3")
    ap.add_argument("--out_prefix", default=f"{RES}/06")
    args = ap.parse_args()

    os.environ["CUDA_VISIBLE_DEVICES"] = args.gpu
    import torch, scvi
    scvi.settings.seed = 0
    log("torch", torch.__version__, "cuda", torch.cuda.is_available(),
        "device", torch.cuda.get_device_name(0) if torch.cuda.is_available() else "cpu")

    ref = sc.read_h5ad(args.ref)
    qry = sc.read_h5ad(args.query)
    log("ref", ref.shape, "| query", qry.shape)

    # ---------- 共同基因 ----------
    common = ref.var_names.intersection(qry.var_names)
    log("共同基因:", len(common))
    ref = ref[:, common].copy()
    qry = qry[:, common].copy()

    # ---------- 原始 counts ----------
    for A in (ref, qry):
        if "counts" not in A.layers:
            A.layers["counts"] = A.X.copy()

    # ---------- HVG（只在參考上選，排除 Ig/MT/RP） ----------
    bl = gene_blocklist(common)
    log("blocklist 基因數:", len(bl))
    tmp = ref.copy()
    sc.pp.normalize_total(tmp, target_sum=1e4)
    sc.pp.log1p(tmp)
    sc.pp.highly_variable_genes(tmp, n_top_genes=args.n_hvg + len(bl),
                                batch_key=args.ref_batch if args.ref_batch in tmp.obs else None,
                                flavor="seurat")
    hvg = [g for g in tmp.var_names[tmp.var.highly_variable] if g not in bl][: args.n_hvg]
    log("HVG 使用數:", len(hvg))
    del tmp

    ref_h = ref[:, hvg].copy()
    qry_h = qry[:, hvg].copy()
    ref_h.X = ref_h.layers["counts"]
    qry_h.X = qry_h.layers["counts"]

    # ---------- 訓練 scVI（只用正常細胞） ----------
    bkey = args.ref_batch if args.ref_batch in ref_h.obs else None
    scvi.model.SCVI.setup_anndata(ref_h, layer="counts", batch_key=bkey)
    model = scvi.model.SCVI(ref_h, n_latent=args.n_latent, n_layers=2,
                            gene_likelihood="nb", dispersion="gene-batch" if bkey else "gene")
    log("訓練參考 scVI ...")
    model.train(max_epochs=args.epochs_ref, early_stopping=True,
                plan_kwargs={"weight_decay": 0.0})
    model.save(f"{args.out_prefix}_scvi_ref", overwrite=True)
    Zr = model.get_latent_representation()
    log("Z_ref", Zr.shape)

    # ---------- scArches surgery：把 MCL 投進去 ----------
    if bkey:
        qry_h.obs[bkey] = qry_h.obs[args.query_batch].astype(str).values
    scvi.model.SCVI.prepare_query_anndata(qry_h, f"{args.out_prefix}_scvi_ref")
    qmodel = scvi.model.SCVI.load_query_data(qry_h, f"{args.out_prefix}_scvi_ref")
    log("訓練 query surgery ...")
    qmodel.train(max_epochs=args.epochs_query, plan_kwargs={"weight_decay": 0.0})
    Zq = qmodel.get_latent_representation()
    log("Z_query", Zq.shape)

    np.save(f"{args.out_prefix}_Zref.npy", Zr)
    np.save(f"{args.out_prefix}_Zqry.npy", Zq)

    # ---------- kNN：找每個 MCL 細胞的正常對應狀態 ----------
    nn = NearestNeighbors(n_neighbors=args.k, algorithm="auto", n_jobs=-1).fit(Zr)
    dist, idx = nn.kneighbors(Zq)
    log("kNN 完成", idx.shape)

    # 參考自身的 kNN 距離分布 → 用來標準化 novelty
    d_self, _ = nn.kneighbors(Zr, n_neighbors=args.k + 1)
    base = np.median(d_self[:, 1:].mean(axis=1))
    novelty = dist.mean(axis=1) / base
    log("novelty score: median=%.2f  p90=%.2f" % (np.median(novelty), np.percentile(novelty, 90)))

    # 匹配到的正常細胞型別（多數決）
    if args.annot in ref.obs:
        ann = ref.obs[args.annot].astype(str).values
        matched = []
        for row in idx:
            vals, cnt = np.unique(ann[row], return_counts=True)
            matched.append(vals[np.argmax(cnt)])
        qry.obs["matched_normal_state"] = matched
        log("匹配到的正常狀態分布:\n",
            pd.Series(matched).value_counts().head(20).to_string())
    qry.obs["novelty"] = novelty

    # ---------- 期望表現與殘差 ----------
    log("計算期望表現 ...")
    refn = ref.copy()
    sc.pp.normalize_total(refn, target_sum=1e4); sc.pp.log1p(refn)
    qn = qry.copy()
    sc.pp.normalize_total(qn, target_sum=1e4); sc.pp.log1p(qn)

    Xr = refn.X.toarray() if sp.issparse(refn.X) else np.asarray(refn.X)
    Xq = qn.X.toarray() if sp.issparse(qn.X) else np.asarray(qn.X)
    # 鄰居平均（分批以控制記憶體）
    exp = np.zeros_like(Xq, dtype=np.float32)
    step = 2000
    for s in range(0, idx.shape[0], step):
        e = min(s + step, idx.shape[0])
        exp[s:e] = Xr[idx[s:e]].mean(axis=1)
    resid = Xq - exp
    log("residual", resid.shape, "mean|r| = %.4f" % np.abs(resid).mean())

    np.save(f"{args.out_prefix}_expected.npy", exp.astype(np.float32))
    np.save(f"{args.out_prefix}_residual.npy", resid.astype(np.float32))
    qn.layers["expected"] = exp
    qn.layers["residual"] = resid
    qn.write_h5ad(f"{args.out_prefix}_mcl_with_residual.h5ad", compression="gzip")
    refn.write_h5ad(f"{args.out_prefix}_ref_lognorm.h5ad", compression="gzip")
    log("完成，輸出於", args.out_prefix)


if __name__ == "__main__":
    main()
