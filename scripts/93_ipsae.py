#!/usr/bin/env python
"""
步驟 93：為既有的 Boltz-2 預測加算 ipSAE，並套用濕實驗校準過的判定門檻。

為什麼要做這一步
----------------
步驟 87 用 ipTM 排序設計。但三份獨立 benchmark 都指出 ipTM 不是最佳判別指標：

  ① Overath et al., bioRxiv 10.1101/2025.08.14.670059
     3,766 個實驗測試過的 binder / 15 個標的，統一 pipeline 重跑。
     結論：以介面為焦點的指標（最顯著是 ipSAE）勝過 ipTM 與 ipAE，
     平均精確度（average precision）提高 1.4 倍。
     作者建議門檻：ipSAE_min > 0.61。
  ② BindCraft 的預設 filter（直接取自其 default_filters.json）
     i_pTM ≥ 0.50、i_pAE ≤ 0.35（= 原始 PAE 10.85 Å）、pLDDT ≥ 0.80。
  ③ Bennett et al., Nat Commun 2023（獨立實驗室、獨立 pipeline）
     pae_interaction < 10 Å。
  ②③ 收斂到幾乎相同的數字（10.85 vs 10），這是該領域最穩固的單一結論。

  ⚠️ 反面證據：Ramon et al. 2026 實測 ipSAE（0.71–0.81）**並未**勝過 ipTM（0.77–0.82）。
     所以本腳本**同時**報告兩者，不用其中一個取代另一個。

ipSAE 定義（Dunbrack, bioRxiv 10.1101/2025.02.10.637595）
---------------------------------------------------------
對鏈對 (A→B)，只取跨鏈且 PAE < 門檻的殘基對：
  n0res[i] = 符合條件的 j 的個數
  d0res[i] = max(1.0, 1.24 * (max(26, n0res[i]) - 15)^(1/3) - 1.8)
  ipsae_byres[i] = mean_j 1/(1 + (PAE[i,j] / d0res[i])^2)
  ipsae_asym(A→B) = max_i ipsae_byres[i]
  ipSAE_min = min(ipsae_asym(A→B), ipsae_asym(B→A))   ← benchmark 用的那個

本腳本自行實作上述公式（公式來源為官方 ipsae.py，已逐行核對），
不執行下載來的第三方程式。
"""
import os, sys, glob, json, math, re
import numpy as np
import pandas as pd

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
RES = f"{BASE}/results"
PRED = f"{RES}/87_boltz_replicates/boltz_results_inputs/predictions"
PAE_CUT = 10.0          # 與 Overath / Bennett 一致
DIST_CUT = 10.0

# 濕實驗校準門檻
THR = dict(ipsae_min=0.61,      # Overath et al.
           iptm=0.50,           # BindCraft
           ipae=10.85,          # BindCraft 0.35 × 31；Bennett 用 10.0
           plddt=0.80)          # BindCraft


def log(*a):
    print("[93]", *a, flush=True)


def calc_d0(L):
    L = np.maximum(26.0, np.asarray(L, dtype=float))
    return np.maximum(1.0, 1.24 * (L - 15.0) ** (1.0 / 3.0) - 1.8)


def atoms_from_cif(path):
    """取每個殘基的鏈代號、CB 座標（無 CB 用 CA）與 pLDDT。

    pDockQ / pDockQ2 需要座標與 pLDDT；ipSAE 只需要鏈代號。
    """
    hdr, rows = [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith("_atom_site."):
                hdr.append(line.strip().split(".", 1)[1])
            elif line.startswith(("ATOM", "HETATM")):
                q = line.split()
                if len(q) >= len(hdr):
                    rows.append(q[:len(hdr)])
    if not rows:
        return None
    df = pd.DataFrame(rows, columns=hdr)
    an = "label_atom_id" if "label_atom_id" in df else "auth_atom_id"
    ch = "label_asym_id" if "label_asym_id" in df else "auth_asym_id"
    rn = "label_seq_id" if "label_seq_id" in df else "auth_seq_id"
    bf = "B_iso_or_equiv" if "B_iso_or_equiv" in df else None
    for c in ("Cartn_x", "Cartn_y", "Cartn_z"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df[rn] = pd.to_numeric(df[rn], errors="coerce")
    if bf:
        df[bf] = pd.to_numeric(df[bf], errors="coerce")
    # 每個殘基優先取 CB，無 CB（甘胺酸）取 CA
    sel = df[df[an].isin(["CA", "CB"])].copy()
    sel["pref"] = (sel[an] == "CB").astype(int)
    sel = sel.sort_values([ch, rn, "pref"]).drop_duplicates(subset=[ch, rn], keep="last")
    xyz = sel[["Cartn_x", "Cartn_y", "Cartn_z"]].values
    pl = sel[bf].values if bf else np.full(len(sel), np.nan)
    return sel[ch].values, sel[rn].values, xyz, pl


def chains_from_cif(path):
    """取每個殘基的鏈代號（只用 CA 原子，一殘基一筆）。"""
    hdr, rows = [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith("_atom_site."):
                hdr.append(line.strip().split(".", 1)[1])
            elif line.startswith(("ATOM", "HETATM")):
                p = line.split()
                if len(p) >= len(hdr):
                    rows.append(p[:len(hdr)])
    if not rows:
        return None
    df = pd.DataFrame(rows, columns=hdr)
    an = "label_atom_id" if "label_atom_id" in df else "auth_atom_id"
    ch = "label_asym_id" if "label_asym_id" in df else "auth_asym_id"
    rn = "label_seq_id" if "label_seq_id" in df else "auth_seq_id"
    df = df[df[an] == "CA"]
    df = df.drop_duplicates(subset=[ch, rn])
    return df[ch].values, pd.to_numeric(df[rn], errors="coerce").values


def ipsae_pair(pae, chains, a, b, cut=PAE_CUT):
    """回傳 (asym A→B, asym B→A, min, max)。"""
    out = {}
    for c1, c2 in [(a, b), (b, a)]:
        m = np.outer(chains == c1, chains == c2) & (pae < cut)
        n0 = m.sum(axis=1)
        d0 = calc_d0(n0)
        best = 0.0
        rows = np.where((chains == c1) & (n0 > 0))[0]
        for i in rows:
            v = 1.0 / (1.0 + (pae[i, m[i]] / d0[i]) ** 2)
            best = max(best, float(v.mean()))
        out[(c1, c2)] = best
    x, y = out[(a, b)], out[(b, a)]
    return x, y, min(x, y), max(x, y)


def interface_pae(pae, chains, a, b):
    """介面 PAE：跨鏈 PAE 的平均（兩個方向都算）。BindCraft 的 i_pAE 對應量。"""
    m = np.outer(chains == a, chains == b) | np.outer(chains == b, chains == a)
    return float(pae[m].mean()) if m.any() else np.nan


def extra_metrics(pae, ch, xyz, pl, a, b):
    """pDockQ、pDockQ2、LIS —— 公式逐行核對官方 ipsae.py。

    pDockQ  : x = mean_plddt(介面殘基) * log10(npairs)
              0.724 / (1 + exp(-0.052*(x-152.611))) + 0.018
    pDockQ2 : x = mean_plddt(介面殘基) * mean(ptm(PAE, d0=10))
              1.31 / (1 + exp(-0.075*(x-84.733))) + 0.005
    LIS     : PAE<12 的跨鏈對，取 mean((12-PAE)/12)
    """
    CUT = 8.0
    dist = np.sqrt(((xyz[:, None, :] - xyz[None, :, :]) ** 2).sum(-1))
    ma, mb = (ch == a), (ch == b)
    cross = np.outer(ma, mb)
    near = cross & (dist <= CUT)
    npairs = int(near.sum())
    if npairs == 0:
        return 0.0, 0.0, 0.0
    involved = np.unique(np.r_[np.where(near.any(1))[0], np.where(near.any(0))[0]])
    mean_pl = float(np.nanmean(pl[involved]))
    x1 = mean_pl * math.log10(npairs)
    pdockq = 0.724 / (1 + math.exp(-0.052 * (x1 - 152.611))) + 0.018
    ptm_vals = 1.0 / (1.0 + (pae[near] / 10.0) ** 2)
    x2 = mean_pl * float(ptm_vals.mean())
    pdockq2 = 1.31 / (1 + math.exp(-0.075 * (x2 - 84.733))) + 0.005
    sel = pae[cross]
    v = sel[sel < 12]
    lis = float(((12 - v) / 12).mean()) if v.size else 0.0
    return round(pdockq, 4), round(pdockq2, 4), round(lis, 4)


rows = []
dirs = sorted([d for d in glob.glob(f"{PRED}/*") if os.path.isdir(d)])
log(f"找到 {len(dirs)} 個預測目錄")
for d in dirs:
    name = os.path.basename(d)
    for cif in sorted(glob.glob(f"{d}/*_model_*.cif")):
        stem = os.path.basename(cif).replace(".cif", "")
        npz = f"{d}/pae_{stem}.npz"
        cjs = f"{d}/confidence_{stem}.json"
        if not os.path.exists(npz):
            continue
        z = np.load(npz)
        key = [k for k in z.files if "pae" in k.lower()] or list(z.files)
        pae = np.asarray(z[key[0]], dtype=float)
        if pae.ndim == 3:
            pae = pae[0]
        got = atoms_from_cif(cif)
        if got is None:
            continue
        ch, rn, xyz, pl = got
        if len(ch) != pae.shape[0]:
            log(f"  {stem}: 殘基數不符（cif {len(ch)} vs pae {pae.shape[0]}），跳過")
            continue
        uc = list(pd.unique(ch))
        if len(uc) < 2:
            continue
        a, b = uc[0], uc[1]
        x, y, mn, mx = ipsae_pair(pae, ch, a, b)
        ipae = interface_pae(pae, ch, a, b)
        pdq, pdq2, lis = extra_metrics(pae, ch, xyz, pl, a, b)
        conf = json.load(open(cjs)) if os.path.exists(cjs) else {}
        rows.append(dict(
            design=name, model=stem.split("_model_")[-1],
            n_res=len(ch), chain_a=a, chain_b=b, len_a=int((ch == a).sum()),
            len_b=int((ch == b).sum()),
            ipsae_a2b=round(x, 4), ipsae_b2a=round(y, 4),
            ipsae_min=round(mn, 4), ipsae_max=round(mx, 4),
            iface_pae=round(ipae, 3),
            pdockq=pdq, pdockq2=pdq2, lis=lis,
            iptm=round(float(conf.get("iptm", np.nan)), 4),
            protein_iptm=round(float(conf.get("protein_iptm", np.nan)), 4),
            complex_plddt=round(float(conf.get("complex_plddt", np.nan)), 4),
            complex_iplddt=round(float(conf.get("complex_iplddt", np.nan)), 4),
            complex_ipde=round(float(conf.get("complex_ipde", np.nan)), 4)))

df = pd.DataFrame(rows)
if df.empty:
    log("沒有可計算的預測"); sys.exit(1)
df.to_csv(f"{RES}/93_ipsae_per_model.csv", index=False)
log(f"逐 model 計算完成：{len(df)} 筆")

KIND = {"CTRL_POS_NRF2": "陽性對照（天然 NRF2）",
        "CTRL_NEG_1": "陰性對照（亂序）", "CTRL_NEG_2": "陰性對照（亂序）"}
df["kind"] = df.design.map(lambda s: KIND.get(s, "設計"))

# benchmark 做法：每個設計取多 seed 的 ipSAE_min 最大值
g = df.groupby(["design", "kind"]).agg(
    n=("ipsae_min", "size"),
    ipsae_min_max=("ipsae_min", "max"),       # Anthropic 報告的做法
    ipsae_min_mean=("ipsae_min", "mean"),
    ipsae_min_sd=("ipsae_min", "std"),
    iptm_mean=("iptm", "mean"), iptm_sd=("iptm", "std"),
    iface_pae_mean=("iface_pae", "mean"),
    pdockq_mean=("pdockq", "mean"), pdockq2_mean=("pdockq2", "mean"),
    lis_mean=("lis", "mean"), lis_max=("lis", "max"),
    plddt_mean=("complex_plddt", "mean"),
    iplddt_mean=("complex_iplddt", "mean")).reset_index()

g["pass_ipsae"] = g.ipsae_min_max > THR["ipsae_min"]
g["pass_iptm"] = g.iptm_mean >= THR["iptm"]
g["pass_ipae"] = g.iface_pae_mean <= THR["ipae"]
g["pass_plddt"] = g.plddt_mean >= THR["plddt"]
g["n_pass"] = g[["pass_ipsae", "pass_iptm", "pass_ipae", "pass_plddt"]].sum(axis=1)
g["verdict"] = np.where(g.n_pass == 4, "全部通過",
                        np.where(g.n_pass >= 3, "通過 3/4", "未通過"))
g = g.sort_values("ipsae_min_max", ascending=False)
g.to_csv(f"{RES}/93_ipsae_summary.csv", index=False)

pd.set_option("display.width", 250)
log("\n" + "=" * 100)
log("=== ipSAE_min（多 seed 取最大）＋ 濕實驗校準門檻 ===")
log(f"門檻：ipSAE_min > {THR['ipsae_min']}（Overath）｜ipTM ≥ {THR['iptm']}、"
    f"介面 PAE ≤ {THR['ipae']} Å、pLDDT ≥ {THR['plddt']}（BindCraft/Bennett）")
cols = ["design", "kind", "n", "ipsae_min_max", "ipsae_min_mean",
        "iptm_mean", "iface_pae_mean", "pdockq_mean", "pdockq2_mean", "lis_mean",
        "plddt_mean", "n_pass", "verdict"]
log("\n" + g[cols].round(4).to_string(index=False))

# 兩個指標是否同意？
from scipy.stats import spearmanr, mannwhitneyu
r = spearmanr(df.ipsae_min, df.iptm)
log(f"\nipSAE_min vs ipTM 的 Spearman（逐 model, n={len(df)}）= {r.statistic:+.3f} (p={r.pvalue:.3g})")
des = df[df.kind == "設計"].ipsae_min
neg = df[df.kind.str.startswith("陰性")].ipsae_min
if len(des) > 2 and len(neg) > 2:
    u = mannwhitneyu(des, neg, alternative="greater")
    log(f"設計 vs 陰性對照（ipSAE_min, 全部 seed）: {des.median():.4f} vs {neg.median():.4f}, "
        f"Mann-Whitney p = {u.pvalue:.4g}")
    u2 = mannwhitneyu(df[df.kind == "設計"].iptm, df[df.kind.str.startswith("陰性")].iptm,
                      alternative="greater")
    log(f"同一批資料用 ipTM 做同樣檢定: p = {u2.pvalue:.4g}（供比較）")
log("\n寫出 93_ipsae_per_model.csv、93_ipsae_summary.csv")
log("DONE")
