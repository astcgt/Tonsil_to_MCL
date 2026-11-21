#!/usr/bin/env python
"""
步驟 102：用 StaB-ddG 做正交的「結合 ΔΔG」裁決。

為什麼要這一步
--------------
3.4 與 80 步的讀出是 Boltz-2 的共摺疊信心分數（ipTM / ipde）。
那量的是「模型能不能建出這個界面」，不是「結合有多強」。
StaB-ddG（ICML 2025）的歸納偏差完全不同：
反向摺疊對數機率 + 熱力學循環，整個算式裡沒有任何共摺疊信心分數。
所以它是真正的第二意見，不是同一個量的另一種包裝。

ΔΔG 的符號約定與上游一致：**負 = 穩定化、正 = 去穩定化**。

做了什麼
--------
1. 結構用實驗結構 7X5E（NRF2–MafG 異源二聚體 + DNA），不用預測結構。
   鏈 A = MafG（UniProt O15525 編號 23–121）、鏈 B = NRF2（Q16236 452–560）。
   鏈 E/F 是同一個不對稱單位裡的**第二份拷貝**，當獨立的結構重複。
2. DNA 鏈（C/D/G/H）與水分子移除：StaB-ddG 是蛋白質–蛋白質的量，
   ProteinMPNN 的 parse_PDB 只讀主鏈 N/CA/C/O，核酸會被當成 gap。
   這是必須誠實交代的限制（見 REPORT）。
3. 殘基重新編號成每條鏈從 1 開始（StaB-ddG 的硬性要求），
   並保留 UniProt ↔ 新編號的對照表。
   位置對不上會被上游的 `assert mut_seq[mut_pos] == wt_aa` 當場抓到。
4. 突變集合取自 80 步的全部單點突變（二聚面、DNA 面、無關表面對照）。
   DNA 面與對照面的突變在這裡的作用是**空間專一性對照**：
   它們對 MafG–NRF2 的結合不該有大影響。
   同時算出每個突變殘基到夥伴鏈的最短重原子距離，當連續的共變數。
5. 反對稱性自我稽核：對每個突變另外做一份「突變型主鏈」——
   把該殘基的名稱換成突變後的胺基酸、只留主鏈原子（座標完全不動），
   再跑反向突變。StaB-ddG 的 ΔΔG 是同一個主鏈上兩條序列的對數機率差，
   所以反向的值在期望上必須正好是 −1 倍。
   殘差因此直接量到 Monte-Carlo 雜訊，也會抓到我自己改檔案時的任何錯誤。

輸出
----
results/102_stabddg.csv        每個突變的 ΔΔG（AB 拷貝、EF 拷貝、反向、重複試驗）
results/102_resmap.csv         UniProt ↔ 新編號對照
results/102_summary.json       相關性、專一性、反對稱性殘差
"""
import os, sys, json, shutil, subprocess, argparse
import numpy as np
import pandas as pd

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
RES = f"{BASE}/results"
REPO = f"{BASE}/third_party/StaB-ddG"
PY = f"{BASE}/env/protdes/bin/python"
PDB7X5E = f"{BASE}/data/structures/7x5e.pdb"
WORK = f"{RES}/102_stabddg"

AA3to1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M",
}
AA1to3 = {v: k for k, v in AA3to1.items() if k != "MSE"}
BACKBONE = {"N", "CA", "C", "O"}

# 7X5E 的兩份拷貝：(複合體名稱, MafG 鏈, NRF2 鏈)
COPIES = [("mafg_nrf2_AB", "A", "B"), ("mafg_nrf2_EF", "E", "F")]


# ---------------------------------------------------------------- PDB 讀寫
def read_protein(path):
    """只讀 ATOM 的標準胺基酸，回傳 {chain: [(reskey, resname, [lines])]}（保序）。"""
    chains = {}
    index = {}
    for ln in open(path):
        if ln[:4] != "ATOM" and not (ln[:6] == "HETATM" and ln[17:20] == "MSE"):
            continue
        altloc = ln[16]
        if altloc not in (" ", "A"):
            continue
        rname = ln[17:20].strip()
        if rname not in AA3to1:
            continue
        ch = ln[21]
        key = ln[22:27]  # resseq + icode
        chains.setdefault(ch, [])
        index.setdefault(ch, {})
        if key not in index[ch]:
            index[ch][key] = len(chains[ch])
            chains[ch].append([key, rname, []])
        chains[ch][index[ch][key]][2].append(ln)
    return chains


def atom_xyz(lines, heavy_only=True):
    out = []
    for ln in lines:
        name = ln[12:16].strip()
        if heavy_only and name.startswith("H"):
            continue
        out.append([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
    return np.asarray(out, float)


def write_complex(path, chains, order, mutate=None):
    """
    寫出重新編號（每鏈從 1 起算）的複合體。
    mutate = (chain, reskey, new_aa1)：把該殘基改名、只留主鏈原子（座標不動）。
    """
    with open(path, "w") as fh:
        serial = 1
        for ch in order:
            for newpos, (key, rname, lines) in enumerate(chains[ch], start=1):
                rn = rname
                keep = lines
                if mutate and mutate[0] == ch and mutate[1] == key:
                    rn = AA1to3[mutate[2]]
                    keep = [l for l in lines if l[12:16].strip() in BACKBONE]
                for ln in keep:
                    out = (ln[:6] + f"{serial:5d}" + ln[11:17] + f"{rn:>3s}" + ln[20:22]
                           + f"{newpos:4d}" + " " + ln[27:])
                    fh.write("ATOM  " + out[6:])
                    serial += 1
            fh.write(f"TER   {serial:5d}\n")
            serial += 1
        fh.write("END\n")


# ---------------------------------------------------------------- 準備
def prep():
    os.makedirs(f"{WORK}/pdbs_fwd", exist_ok=True)
    os.makedirs(f"{WORK}/pdbs_rev", exist_ok=True)
    chains = read_protein(PDB7X5E)
    for ch in chains:
        print(f"  鏈 {ch}: {len(chains[ch])} 個殘基 "
              f"({chains[ch][0][0].strip()}–{chains[ch][-1][0].strip()})")

    # ---- UniProt ↔ 新編號對照，並檢查鏈內連續
    rmap = []
    for name, cm, cn in COPIES:
        for role, ch in (("MafG", cm), ("NRF2", cn)):
            for newpos, (key, rname, _) in enumerate(chains[ch], start=1):
                rmap.append(dict(copy=name, role=role, chain=ch,
                                 uniprot_pos=int(key.strip()), new_pos=newpos,
                                 aa=AA3to1[rname]))
    rmap = pd.DataFrame(rmap)
    rmap.to_csv(f"{RES}/102_resmap.csv", index=False)

    # ---- 突變集合（80 步的全部單點突變，去重）
    pert = pd.read_csv(f"{RES}/80_perturb_correct.csv")
    pert = pert[pert["pos"].notna()].copy()
    pert["pos"] = pert["pos"].astype(int)

    def face(grp):
        if "二聚面" in grp:
            return "dimer"
        if "DNA面" in grp or "電荷反轉" in grp:
            return "dna"
        return "control"

    pert["face"] = pert["grp"].map(face)
    muts = (pert.groupby(["pos", "wt", "mt"], as_index=False)
                .agg(face=("face", lambda s: sorted(set(s))[0]),
                     esm_llr=("esm_llr", "first")))
    muts["mut_id"] = muts.wt + muts.pos.astype(str) + muts.mt
    print(f"  突變集合：{len(muts)} 個單點突變 "
          f"（{(muts.face=='dimer').sum()} 二聚面 / {(muts.face=='dna').sum()} DNA 面 "
          f"/ {(muts.face=='control').sum()} 對照）")

    # ---- 寫正向複合體，並檢查野生型胺基酸對得上
    rows = []
    for name, cm, cn in COPIES:
        write_complex(f"{WORK}/pdbs_fwd/{name}.pdb", chains, [cm, cn])
        keys = {int(k.strip()): (i + 1, AA3to1[r])
                for i, (k, r, _) in enumerate(chains[cm])}
        coords = {int(k.strip()): atom_xyz(l) for k, _, l in chains[cm]}
        # 到夥伴鏈的最短重原子距離
        partner = np.vstack([atom_xyz(l) for _, _, l in chains[cn]])
        for _, m in muts.iterrows():
            if m.pos not in keys:
                print(f"  ⚠️ {name}: 位置 {m.pos} 不在鏈 {cm} 的結構範圍內，跳過")
                continue
            newpos, aa = keys[m.pos]
            assert aa == m.wt, f"{name} {cm}{m.pos}: 結構是 {aa}，突變表寫 {m.wt}"
            xyz = coords[m.pos]
            d = float(np.sqrt(((xyz[:, None, :] - partner[None, :, :]) ** 2)
                              .sum(-1)).min())
            rows.append(dict(copy=name, mut_id=m.mut_id, face=m.face,
                             uniprot_pos=m.pos, wt=m.wt, mt=m.mt,
                             esm_llr=m.esm_llr, chain=cm, new_pos=newpos,
                             dist_partner=round(d, 2),
                             pdb=f"{name}_{cm}_{cn}",
                             mutation=f"{m.wt}{cm}{newpos}{m.mt}"))
    fwd = pd.DataFrame(rows)
    fwd[["pdb", "mutation"]].rename(columns={"pdb": "#Pdb"}).to_csv(
        f"{WORK}/mutations_fwd.csv", index=False)

    # ---- 反向：突變型主鏈 + 反向突變（只做 AB 拷貝）
    name, cm, cn = COPIES[0]
    keymap = {int(k.strip()): k for k, _, _ in chains[cm]}
    rev = []
    for _, m in fwd[fwd["copy"] == name].iterrows():
        tag = f"revAB{m.mut_id}"
        write_complex(f"{WORK}/pdbs_rev/{tag}.pdb", chains, [cm, cn],
                      mutate=(cm, keymap[m.uniprot_pos], m.mt))
        rev.append(dict(**{"#Pdb": f"{tag}_{cm}_{cn}"},
                        mutation=f"{m.mt}{cm}{m.new_pos}{m.wt}", mut_id=m.mut_id))
    rev = pd.DataFrame(rev)
    rev[["#Pdb", "mutation"]].to_csv(f"{WORK}/mutations_rev.csv", index=False)
    fwd.to_csv(f"{WORK}/fwd_meta.csv", index=False)
    rev.to_csv(f"{WORK}/rev_meta.csv", index=False)
    print(f"  正向 {len(fwd)} 列（兩份拷貝）、反向 {len(rev)} 列（AB 拷貝）")
    return fwd, rev


def prep_scan():
    """
    MafG–NRF2 界面的完整電腦丙胺酸掃描。
    這仍然在 StaB-ddG 的驗證範圍內（天然界面、實驗結構），
    而且全部突變共用一個結構，所以 190 個突變和 13 個一樣快。
    問的問題是：要擋住 MafG–NRF2，競爭性抑制劑該蓋住哪一塊表面。
    """
    chains = read_protein(PDB7X5E)
    name, cm, cn = COPIES[0]
    os.makedirs(f"{WORK}/pdbs_fwd", exist_ok=True)
    write_complex(f"{WORK}/pdbs_fwd/{name}.pdb", chains, [cm, cn])
    rows = []
    for role, ch, other in (("MafG", cm, cn), ("NRF2", cn, cm)):
        partner = np.vstack([atom_xyz(l) for _, _, l in chains[other]])
        for newpos, (key, rname, lines) in enumerate(chains[ch], start=1):
            aa = AA3to1[rname]
            if aa == "A":
                continue
            d = float(np.sqrt(((atom_xyz(lines)[:, None, :] - partner[None, :, :]) ** 2)
                              .sum(-1)).min())
            rows.append(dict(role=role, chain=ch, uniprot_pos=int(key.strip()),
                             new_pos=newpos, wt=aa, mt="A",
                             dist_partner=round(d, 2),
                             pdb=f"{name}_{cm}_{cn}",
                             mutation=f"{aa}{ch}{newpos}A"))
    scan = pd.DataFrame(rows)
    scan[["pdb", "mutation"]].rename(columns={"pdb": "#Pdb"}).to_csv(
        f"{WORK}/mutations_scan.csv", index=False)
    scan.to_csv(f"{WORK}/scan_meta.csv", index=False)
    print(f"  丙胺酸掃描：{len(scan)} 個位置"
          f"（MafG {(scan.role=='MafG').sum()}、NRF2 {(scan.role=='NRF2').sum()}）")
    return scan


def analyze_scan(trials):
    scan = pd.read_csv(f"{WORK}/scan_meta.csv")
    out = parse_out(pd.read_csv(f"{WORK}/out_scan/output.csv"), trials)
    t = scan.merge(out, left_on=["pdb", "mutation"], right_on=["pdb", "mutation"],
                   how="left", validate="one_to_one").drop(columns=["pdb"])
    name, cm, cn = COPIES[0]
    for role, ch in (("MafG", cm), ("NRF2", cn)):
        bur = interface_burial(f"{WORK}/pdbs_fwd/mafg_nrf2_AB.pdb", ch)
        m = t.role == role
        t.loc[m, "dsasa"] = [round(bur[p][0], 1) for p in t.loc[m, "new_pos"]]
    t = t.sort_values("ddg", ascending=False)
    t.to_csv(f"{RES}/102_alascan.csv", index=False)

    from scipy.stats import spearmanr
    hot = t[t.ddg >= 1.0]
    s = dict(n=int(len(t)),
             n_hotspot=int(len(hot)),
             hotspot_threshold_kcal=1.0,
             noise_floor=dict(
                 n_far=int((t.dist_partner > 10).sum()),
                 max_abs_ddg_far=float(t.loc[t.dist_partner > 10, "ddg"]
                                       .abs().max()),
                 trial_sd_median=float(np.nanmedian(t.ddg_sd))),
             vs_dsasa=dict(zip(("r", "p"), [float(x) for x in spearmanr(
                 t.ddg, t.dsasa)])),
             vs_dist=dict(zip(("r", "p"), [float(x) for x in spearmanr(
                 t.ddg, t.dist_partner)])))
    for role in ("MafG", "NRF2"):
        sub = t[t.role == role]
        s[role] = dict(
            n_hotspot=int((sub.ddg >= 1.0).sum()),
            top=[dict(res=f"{r.wt}{r.uniprot_pos}", ddg=round(float(r.ddg), 2),
                      dist=float(r.dist_partner), dsasa=float(r.dsasa))
                 for r in sub.head(12).itertuples()],
            hotspot_span=[int(sub.loc[sub.ddg >= 1.0, "uniprot_pos"].min()),
                          int(sub.loc[sub.ddg >= 1.0, "uniprot_pos"].max())])
    with open(f"{RES}/102_alascan_summary.json", "w") as fh:
        json.dump(s, fh, indent=2, ensure_ascii=False)
    print(json.dumps(s, indent=2, ensure_ascii=False))
    return t, s


# ---------------------------------------------------------------- 埋藏面積
# 殘基體積（Å³，Zamyatnin 1972），用來算側鏈縮小了多少
VOL = dict(A=88.6, R=173.4, N=114.1, D=111.1, C=108.5, Q=143.8, E=138.4,
           G=60.1, H=153.2, I=166.7, L=166.7, K=168.6, M=162.9, F=189.9,
           P=112.7, S=89.0, T=116.1, W=227.8, Y=193.6, V=140.0)


def interface_burial(pdb_path, mafg_chain):
    """
    每個 MafG 殘基在結合時被埋掉的溶劑可及表面積（Å²）。
    這是必須排除的混淆：反向摺疊類模型對「把埋在界面裡的大疏水殘基換成丙胺酸」
    幾乎一定給出去穩定化的答案，所以要問的是「ΔΔG 有沒有超出單純幾何埋藏」。
    """
    from Bio.PDB import PDBParser
    from Bio.PDB.SASA import ShrakeRupley
    sr = ShrakeRupley()
    st = PDBParser(QUIET=True).get_structure("c", pdb_path)[0]
    sr.compute(st, level="R")
    bound = {r.id[1]: r.sasa for r in st[mafg_chain]}
    st2 = PDBParser(QUIET=True).get_structure("c", pdb_path)[0]
    for ch in [c.id for c in st2 if c.id != mafg_chain]:
        st2.detach_child(ch)
    sr.compute(st2, level="R")
    free = {r.id[1]: r.sasa for r in st2[mafg_chain]}
    return {k: (free[k] - bound[k], free[k]) for k in bound}


# ---------------------------------------------------------------- 執行
def run_stab(pdb_dir, csv_path, out_dir, trials, mc, gpu, seed=0):
    env = dict(os.environ, PYTHONPATH=REPO, CUDA_VISIBLE_DEVICES=str(gpu))
    cmd = [PY, "run_stabddg.py", "--pdb_dir", pdb_dir, "--csv_path", csv_path,
           "--output_dir", out_dir, "--trials", str(trials),
           "--mc_samples", str(mc), "--seed", str(seed),
           "--run_name", "output"]
    print("  $", " ".join(cmd))
    p = subprocess.run(cmd, cwd=REPO, env=env, capture_output=True, text=True)
    if p.returncode != 0:
        print(p.stdout[-3000:]); print(p.stderr[-3000:])
        raise SystemExit(f"StaB-ddG 失敗（returncode {p.returncode}）")
    for line in p.stdout.splitlines():
        if "Finished loading" in line or "Saved predictions" in line:
            print("   ", line)
    return pd.read_csv(f"{out_dir}/output.csv")


def parse_out(df, trials):
    """output.csv 的 Mutation 欄是 "['QA82A']" 字串。取平均與重複試驗標準差。"""
    cols = [f"pred_{i+1}" for i in range(trials)]
    out = df.copy()
    out["mutation"] = out["Mutation"].astype(str).str.strip("[]'\" ")
    out["ddg"] = out[cols].mean(axis=1)
    out["ddg_sd"] = out[cols].std(axis=1, ddof=1) if trials > 1 else np.nan
    out["pdb"] = out["Name"]
    return out[["pdb", "mutation", "ddg", "ddg_sd"] + cols]


# ---------------------------------------------------------------- 主程式
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trials", type=int, default=3)
    ap.add_argument("--mc", type=int, default=20)
    ap.add_argument("--gpu", type=int, default=3)
    ap.add_argument("--skip-run", action="store_true")
    ap.add_argument("--scan-only", action="store_true",
                    help="只做界面丙胺酸掃描")
    a = ap.parse_args()

    if a.scan_only:
        print("== 界面丙胺酸掃描 ==")
        prep_scan()
        if not a.skip_run:
            run_stab(f"{WORK}/pdbs_fwd", f"{WORK}/mutations_scan.csv",
                     f"{WORK}/out_scan", a.trials, a.mc, a.gpu, seed=2)
        analyze_scan(a.trials)
        return

    print("== 1. 準備結構與突變表 ==")
    fwd, rev = prep()

    if not a.skip_run:
        print("== 2. 正向（野生型結構，兩份拷貝）==")
        run_stab(f"{WORK}/pdbs_fwd", f"{WORK}/mutations_fwd.csv",
                 f"{WORK}/out_fwd", a.trials, a.mc, a.gpu, seed=0)
        print("== 3. 反向（突變型主鏈，反對稱性稽核）==")
        run_stab(f"{WORK}/pdbs_rev", f"{WORK}/mutations_rev.csv",
                 f"{WORK}/out_rev", a.trials, a.mc, a.gpu, seed=1)

    print("== 4. 整理 ==")
    of = parse_out(pd.read_csv(f"{WORK}/out_fwd/output.csv"), a.trials)
    orv = parse_out(pd.read_csv(f"{WORK}/out_rev/output.csv"), a.trials)

    f = fwd.merge(of, on=["pdb", "mutation"], how="left", validate="one_to_one")
    ab = f[f["copy"] == "mafg_nrf2_AB"].copy()
    ef = f[f["copy"] == "mafg_nrf2_EF"][["mut_id", "ddg", "ddg_sd", "dist_partner"]]
    ef.columns = ["mut_id", "ddg_EF", "ddg_sd_EF", "dist_partner_EF"]
    tab = ab.merge(ef, on="mut_id", how="left")

    orv["mut_id"] = orv["pdb"].str.replace(r"^revAB(.*)_[A-Z]_[A-Z]$", r"\1", regex=True)
    tab = tab.merge(orv[["mut_id", "ddg", "ddg_sd"]].rename(
        columns={"ddg": "ddg_rev", "ddg_sd": "ddg_sd_rev"}), on="mut_id", how="left")
    tab["antisym_resid"] = tab.ddg + tab.ddg_rev

    # Boltz-2 的讀出接回來（同一個突變、MafG+NRF2 那一側）
    pert = pd.read_csv(f"{RES}/80_perturb_correct.csv")
    pert = pert[(pert.assay == "MafG+NRF2") & pert.pos.notna()].copy()
    pert["mut_id"] = pert.wt + pert.pos.astype(int).astype(str) + pert.mt
    wt = pd.read_csv(f"{RES}/80_perturb_correct.csv")
    wt = wt[wt.name == "WT_MafG_NRF2"].iloc[0]
    pert["d_ipde"] = pert["ipde"] - float(wt["ipde"])
    tab = tab.merge(pert[["mut_id", "d_iptm", "d_ipde"]], on="mut_id", how="left")

    # 幾何混淆：界面埋藏面積與側鏈縮小量
    bur = interface_burial(f"{WORK}/pdbs_fwd/mafg_nrf2_AB.pdb", "A")
    tab["dsasa"] = [round(bur[p][0], 1) for p in tab.new_pos]
    tab["sasa_free"] = [round(bur[p][1], 1) for p in tab.new_pos]
    tab["d_vol"] = [round(VOL[w] - VOL[m], 1) for w, m in zip(tab.wt, tab.mt)]

    tab = tab.drop(columns=["copy", "pdb"]).sort_values("ddg", ascending=False)

    # ---- 統計
    from scipy.stats import pearsonr, spearmanr, mannwhitneyu

    def safe(fn, x, y):
        m = np.isfinite(x) & np.isfinite(y)
        if m.sum() < 4:
            return (np.nan, np.nan, int(m.sum()))
        r, p = fn(x[m], y[m])
        return (float(r), float(p), int(m.sum()))

    s = {}
    s["n_mutations"] = int(len(tab))
    s["copy_replicate"] = dict(zip(("r", "p", "n"),
        safe(pearsonr, tab.ddg.to_numpy(float), tab.ddg_EF.to_numpy(float))))
    s["antisymmetry"] = dict(
        mean_resid=float(np.nanmean(tab.antisym_resid)),
        max_abs_resid=float(np.nanmax(np.abs(tab.antisym_resid))),
        rmse=float(np.sqrt(np.nanmean(tab.antisym_resid ** 2))),
        ddg_range=[float(tab.ddg.min()), float(tab.ddg.max())],
        trial_sd_median=float(np.nanmedian(tab.ddg_sd)))
    for lab, sub in (("all", tab), ("dimer", tab[tab.face == "dimer"])):
        s[f"vs_boltz_{lab}"] = dict(
            iptm=dict(zip(("r", "p", "n"), safe(spearmanr, sub.ddg.to_numpy(float),
                                                sub.d_iptm.to_numpy(float)))),
            ipde=dict(zip(("r", "p", "n"), safe(spearmanr, sub.ddg.to_numpy(float),
                                                sub.d_ipde.to_numpy(float)))),
            esm=dict(zip(("r", "p", "n"), safe(spearmanr, sub.ddg.to_numpy(float),
                                               sub.esm_llr.to_numpy(float)))))
    s["vs_distance"] = dict(zip(("r", "p", "n"),
        safe(spearmanr, tab.ddg.to_numpy(float), tab.dist_partner.to_numpy(float))))
    g = {k: v.ddg.to_numpy(float) for k, v in tab.groupby("face")}
    s["by_face"] = {k: dict(n=len(v), median=float(np.median(v)),
                            max=float(v.max())) for k, v in g.items()}
    if "dimer" in g and "dna" in g:
        u, p = mannwhitneyu(g["dimer"], g["dna"], alternative="greater")
        s["dimer_gt_dna"] = dict(U=float(u), p=float(p))

    # 幾何混淆：ΔΔG 有沒有超出「界面埋藏面積 + 側鏈縮小量」能解釋的範圍
    from scipy.stats import rankdata
    dim = tab[tab.face == "dimer"].copy()
    for lab, sub in (("all", tab), ("dimer", dim)):
        s[f"vs_geometry_{lab}"] = dict(
            dsasa=dict(zip(("r", "p", "n"), safe(spearmanr, sub.ddg.to_numpy(float),
                                                 sub.dsasa.to_numpy(float)))),
            d_vol=dict(zip(("r", "p", "n"), safe(spearmanr, sub.ddg.to_numpy(float),
                                                 sub.d_vol.to_numpy(float)))))
    y = rankdata(dim.ddg)
    X = np.c_[np.ones(len(dim)), rankdata(dim.dsasa), rankdata(dim.d_vol)]
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    dim["ddg_geom_resid"] = y - X @ beta
    s["geometry_adjusted_dimer"] = {
        r.mut_id: round(float(r.ddg_geom_resid), 2)
        for r in dim.sort_values("ddg_geom_resid", ascending=False).itertuples()}
    tab = tab.merge(dim[["mut_id", "ddg_geom_resid"]], on="mut_id", how="left")

    # 兩個方法的共識：各自的去穩定化排名
    dim2 = tab[(tab.face == "dimer") & tab.d_iptm.notna()].copy()
    dim2["rank_stab"] = rankdata(-dim2.ddg)          # 1 = 最去穩定
    dim2["rank_iptm"] = rankdata(dim2.d_iptm)        # 1 = ipTM 掉最多
    dim2["rank_ipde"] = rankdata(-dim2.d_ipde)       # 1 = 界面距離誤差升最多
    s["consensus_dimer"] = {
        r.mut_id: dict(ddg=round(float(r.ddg), 3), rank_stab=int(r.rank_stab),
                       d_ipde=round(float(r.d_ipde), 3), rank_ipde=int(r.rank_ipde),
                       d_iptm=round(float(r.d_iptm), 3), rank_iptm=int(r.rank_iptm))
        for r in dim2.sort_values("rank_stab").itertuples()}
    s["both_top4"] = dict(
        vs_iptm=sorted(dim2.loc[(dim2.rank_stab <= 4) & (dim2.rank_iptm <= 4), "mut_id"]),
        vs_ipde=sorted(dim2.loc[(dim2.rank_stab <= 4) & (dim2.rank_ipde <= 4), "mut_id"]))
    tab.to_csv(f"{RES}/102_stabddg.csv", index=False)

    with open(f"{RES}/102_summary.json", "w") as fh:
        json.dump(s, fh, indent=2, ensure_ascii=False)

    print(json.dumps(s, indent=2, ensure_ascii=False))
    print(tab[["mut_id", "face", "dist_partner", "ddg", "ddg_sd", "ddg_EF",
               "antisym_resid", "d_iptm", "esm_llr"]].to_string(index=False))


if __name__ == "__main__":
    main()
