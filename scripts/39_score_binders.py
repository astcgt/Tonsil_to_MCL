#!/usr/bin/env python
"""
步驟 39：評分並篩選 RFdiffusion3 產出的 IRF4 結合劑設計。

評分標準（都可從結構直接算，不需額外模型）
-------------------------------------------
1. **熱點覆蓋率**：15 個 DNA 接觸熱點中，有幾個被結合劑蓋住（任一原子 ≤5 Å）
   —— 這直接對應「能不能擋住 IRF4 結合 DNA」
2. **界面接觸數**：結合劑與標的之間 ≤5 Å 的殘基對數（界面大小的代理）
3. **結構品質**：RFdiffusion3 自己報的 clash 數、chainbreak、二級結構比例、
   迴旋半徑、丙胺酸含量（過高代表序列沒設計好）
4. **疏水界面比例**：界面殘基中疏水胺基酸的比例（太低不易結合，太高易聚集）

輸出 results/39_binder_ranking.csv 與前幾名的 PDB
"""
import os, sys, gzip, glob, json, re
from collections import defaultdict
import numpy as np
import pandas as pd

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
IN = f"{BASE}/results/37_rfd3_binders"
RES = f"{BASE}/results"
HOTSPOTS = [23, 24, 54, 56, 57, 74, 78, 91, 94, 95, 96, 98, 99, 102, 103]
CUT = 5.0
HYDRO = set("AVLIMFWYC")
AA3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H',
       'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
       'TYR':'Y','VAL':'V'}


def log(*a):
    print("[39]", *a, flush=True)


def read_cif(path):
    """最小 mmCIF atom_site 解析器（RFdiffusion3 把標的與結合劑放在同一條鏈，
    用 auth_seq_id 區分：1–114 是標的 IRF4，115 以後是新設計的結合劑）。"""
    op = gzip.open if path.endswith(".gz") else open
    hdr, rows = [], []
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith("_atom_site."):
                hdr.append(line.strip().split(".", 1)[1])
            elif line.startswith(("ATOM", "HETATM")):
                p = line.split()
                if len(p) >= len(hdr):
                    rows.append(p[:len(hdr)])
    if not rows or not hdr:
        return None
    df = pd.DataFrame(rows, columns=hdr)
    for c in ("Cartn_x", "Cartn_y", "Cartn_z"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["auth_seq_id"] = pd.to_numeric(df["auth_seq_id"], errors="coerce")
    return df.dropna(subset=["Cartn_x", "auth_seq_id"])


N_TARGET = 114       # 7OGS A 鏈 20–133，在輸出中被重新編號成 1–114


def analyse(cif, meta):
    df = read_cif(cif)
    if df is None or len(df) == 0:
        return None
    rs = df["auth_seq_id"].values
    rn_col = "auth_comp_id" if "auth_comp_id" in df else "label_comp_id"
    T = df[rs <= N_TARGET]
    B = df[rs > N_TARGET]
    if len(T) == 0 or len(B) == 0:
        return None
    Tx = T[["Cartn_x", "Cartn_y", "Cartn_z"]].values
    Bx = B[["Cartn_x", "Cartn_y", "Cartn_z"]].values
    d = np.sqrt(((Tx[:, None, :] - Bx[None, :, :]) ** 2).sum(-1))
    near = d <= CUT
    t_res = T["auth_seq_id"].values.astype(int)
    b_res = B["auth_seq_id"].values.astype(int)
    contact_t = set(t_res[near.any(axis=1)])
    contact_b = set(b_res[near.any(axis=0)])
    ii, jj = np.where(near)
    pairs = len({(int(a), int(b)) for a, b in zip(t_res[ii], b_res[jj])})
    # 輸出中的編號 1..114 對應 7OGS 原始殘基 20..133 → 熱點要平移
    hot_renum = {h - 19 for h in HOTSPOTS}
    hot_cov = len(contact_t & hot_renum)

    seqmap = (B.drop_duplicates(subset=["auth_seq_id"])
              .sort_values("auth_seq_id")[["auth_seq_id", rn_col]].values)
    seq = "".join(AA3.get(r, "X") for _, r in seqmap)
    bres = {int(n): AA3.get(r, "X") for n, r in seqmap}
    iface_aa = [bres.get(i, "X") for i in contact_b]
    hyd = np.mean([a in HYDRO for a in iface_aa]) if iface_aa else np.nan

    m = meta.get("metrics", {})
    return dict(
        design=os.path.basename(cif).replace(".cif.gz", "").replace("input_", ""),
        binder_len=len(seq),
        hotspot_covered=hot_cov, hotspot_frac=round(hot_cov / len(HOTSPOTS), 3),
        iface_target_res=len(contact_t), iface_binder_res=len(contact_b),
        iface_pairs=pairs,
        iface_hydrophobic_frac=round(float(hyd), 3) if hyd == hyd else np.nan,
        clashes_sc=m.get("n_clashing.interresidue_clashes_w_sidechain"),
        clashes_bb=m.get("n_clashing.interresidue_clashes_w_backbone"),
        chainbreaks=m.get("n_chainbreaks"),
        helix=round(m.get("helix_fraction", float("nan")), 2),
        sheet=round(m.get("sheet_fraction", float("nan")), 2),
        loop=round(m.get("loop_fraction", float("nan")), 2),
        radius_gyration=round(m.get("radius_of_gyration", float("nan")), 1),
        ala_content=round(m.get("alanine_content", float("nan")), 3),
        sequence=seq)


rows = []
cifs = sorted(glob.glob(f"{IN}/*.cif.gz")) + sorted(glob.glob(f"{IN}/*.cif"))
log(f"找到 {len(cifs)} 個設計")
for c in cifs:
    j = c.replace(".cif.gz", ".json").replace(".cif", ".json")
    meta = json.load(open(j)) if os.path.exists(j) else {}
    try:
        r = analyse(c, meta)
    except Exception as e:
        log(f"  {os.path.basename(c)} 解析失敗: {type(e).__name__}: {e}"); continue
    if r:
        rows.append(r)

if not rows:
    log("沒有可評分的設計"); sys.exit(0)
df = pd.DataFrame(rows)
# 綜合分數：熱點覆蓋為主，界面大小次之，扣掉 clash 與過高丙胺酸
df["score"] = (df.hotspot_frac * 100
               + df.iface_pairs.clip(upper=120) * 0.25
               - df.clashes_sc.fillna(0) * 3
               - df.clashes_bb.fillna(0) * 6
               - (df.ala_content.fillna(0).clip(lower=0.15) - 0.15) * 200)
df = df.sort_values("score", ascending=False)
df.to_csv(f"{RES}/39_binder_ranking.csv", index=False)
pd.set_option("display.width", 240)
log(f"\n=== 共 {len(df)} 個設計，前 15 名 ===")
cols = ["design", "binder_len", "hotspot_covered", "iface_pairs", "iface_binder_res",
        "iface_hydrophobic_frac", "clashes_sc", "helix", "sheet", "ala_content", "score"]
log("\n" + df.head(15)[cols].to_string(index=False))
log("\n=== 前 5 名的結合劑序列 ===")
for _, r in df.head(5).iterrows():
    log(f"  {r.design}  （{r.binder_len} aa，蓋住 {r.hotspot_covered}/15 個熱點）")
    log(f"    {r.sequence}")
log("\n寫出 " + f"{RES}/39_binder_ranking.csv")
