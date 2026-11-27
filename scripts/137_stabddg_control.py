#!/usr/bin/env python
"""步驟 137：StaB-ddG 說「把 MafG 換成 MafK，結合掉 +7.3 kcal/mol」。這個數字可信嗎？

決定性的對照
------------
天然的 NRF2 是公認能和 MafG、MafF、MafK **三個都**形成異源二聚體的。
所以對 NRF2 來說，把夥伴從 MafG 換成 MafK，真實的 ΔΔG 應該接近 0。

→ 拿同一組 15 個取代，在**實驗晶體結構 7X5E 的 MafG–NRF2 界面**上跑一次。
  若也得到 +7 → StaB-ddG 只是在說「改 15 個胺基酸就會扣 7 分」，與選擇性無關，
    那它對我們設計品的 +7.3 就不能採信。
  若接近 0 → 那 +7.3 是真的訊號，設計品確實比 NRF2 更挑標的。

結構用 scripts/102 已經準備好的 results/102_stabddg/pdbs_fwd/mafg_nrf2_AB.pdb
（鏈 A = MafG，UniProt 23–121 重新編號為 1 起；鏈 B = NRF2；DNA 已移除）。

輸出 results/137_*
"""
import os
import json
import glob
import ast
import shutil
import argparse
import subprocess

import pandas as pd
import numpy as np

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
RES = f"{BASE}/results"
AUX = f"{BASE}/data/aux/uniprot"
REPO = f"{BASE}/third_party/StaB-ddG"
PY = f"{BASE}/env/protdes/bin/python"
SRC = f"{RES}/102_stabddg/pdbs_fwd/mafg_nrf2_AB.pdb"
WORK = f"{RES}/137_control"


def log(*a):
    print("[137]", *a, flush=True)


def frag(g):
    return json.load(open(f"{AUX}/{g}.json"))["sequence"]["value"][22:121]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gpu", default="4")
    ap.add_argument("--trials", type=int, default=3)
    ap.add_argument("--mc", type=int, default=20)
    a = ap.parse_args()
    os.makedirs(f"{WORK}/pdb", exist_ok=True)
    shutil.copy(SRC, f"{WORK}/pdb/xtal.pdb")

    from Bio.PDB import PDBParser
    from Bio.PDB.Polypeptide import protein_letters_3to1 as t31
    import warnings
    warnings.filterwarnings("ignore")
    st = PDBParser(QUIET=True).get_structure("x", SRC)[0]
    seqA = {r.id[1]: t31.get(r.get_resname(), "X") for r in st["A"]}
    log("晶體鏈 A（MafG）解析到 %d 個殘基，編號 %d–%d"
        % (len(seqA), min(seqA), max(seqA)))
    log("晶體鏈 B（NRF2）解析到 %d 個殘基" % len(list(st["B"])))

    mg, mf, mk = frag("MAFG"), frag("MAFF"), frag("MAFK")
    sets = {"to_MAFK": mk, "to_MAFF": mf}
    entries, report = [], {}
    for name, other in sets.items():
        muts, skipped = [], []
        for i in range(99):
            if mg[i] == other[i]:
                continue
            p = i + 1
            if p not in seqA:
                skipped.append(p); continue
            if seqA[p] != mg[i]:
                skipped.append(p); continue
            muts.append("%sA%d%s" % (mg[i], p, other[i]))
        report[name] = {"n_applied": len(muts), "n_skipped": len(skipped),
                        "skipped_positions": skipped, "mutations": muts}
        log("%s：%d 個取代可套用，%d 個因為晶體沒解析或序列不符而跳過 %s"
            % (name, len(muts), len(skipped), skipped))
        entries.append((name, ",".join(muts)))
        for m in muts:
            entries.append(("single_%s_%s" % (name[-4:], m), m))

    # 同樣的白金對照：七肽登記白胺酸 → 丙胺酸
    hot = [p for p in range(1, 100)
           if p in seqA and seqA[p] == "L" and 47 <= p <= 85]
    for p in hot:
        entries.append(("ctrl_L%dA" % p, "LA%dA" % p))
    log("加上 %d 個白胺酸對照，共 %d 個項目" % (len(hot), len(entries)))

    csv = f"{WORK}/mutations.csv"
    pd.DataFrame([{"#Pdb": "xtal_A_B", "mutation": m} for _, m in entries]).to_csv(
        csv, index=False)
    pd.DataFrame([{"label": n, "mutation": m} for n, m in entries]).to_csv(
        f"{WORK}/labels.csv", index=False)

    out = f"{WORK}/out"
    os.makedirs(out, exist_ok=True)
    env = dict(os.environ, CUDA_VISIBLE_DEVICES=str(a.gpu))
    cmd = [PY, "run_stabddg.py", "--pdb_dir", f"{WORK}/pdb", "--csv_path", csv,
           "--output_dir", out, "--trials", str(a.trials), "--mc_samples", str(a.mc),
           "--seed", "0", "--run_name", "ctrl"]
    log("跑 StaB-ddG（晶體 MafG–NRF2 界面）…")
    p = subprocess.run(cmd, cwd=REPO, env=env, capture_output=True, text=True)
    log("return code:", p.returncode)
    if p.returncode != 0:
        log(p.stdout[-1200:]); log(p.stderr[-2500:]); return

    D = pd.read_csv(f"{out}/ctrl.csv", index_col=0)
    D["mut"] = D["Mutation"].apply(lambda s: ",".join(ast.literal_eval(s)))
    pc = [c for c in D.columns if c.startswith("pred_")]
    D["ddg"] = D[pc].mean(1); D["sd"] = D[pc].std(1)
    L = pd.read_csv(f"{WORK}/labels.csv")
    M = L.merge(D[["mut", "ddg", "sd"]], left_on="mutation", right_on="mut", how="left")
    M.to_csv(f"{RES}/137_stabddg_control.csv", index=False)

    pd.set_option("display.width", 200)
    whole = M[M.label.str.startswith("to_")]
    ctrl = M[M.label.str.startswith("ctrl_")].sort_values("ddg", ascending=False)
    log("\n=== 晶體 MafG–NRF2 界面：把 MafG 換成旁系同源 ===\n",
        whole[["label", "ddg", "sd"]].round(3).to_string(index=False))
    log("\n=== 同一結構上的白胺酸熱點對照 ===\n",
        ctrl[["label", "ddg", "sd"]].round(3).to_string(index=False))

    # 與設計品的結果並排
    try:
        d1 = json.load(open(f"{RES}/136_stabddg_parsed.json"))
    except Exception:
        d1 = {}
    comp = {}
    for k in ("to_MAFK", "to_MAFF"):
        row = whole[whole.label == k]
        comp[k] = {"NRF2_crystal": float(row.ddg.iloc[0]) if len(row) else None,
                   "designed_binder_predicted": d1.get(k)}
    log("\n=== 並排比較（ΔΔG kcal/mol，正 = 結合變差） ===")
    log("%-10s %22s %26s" % ("轉換", "天然 NRF2（晶體）", "設計品 D1（預測結構）"))
    for k, v in comp.items():
        log("%-10s %22s %26s"
            % (k, "%.2f" % v["NRF2_crystal"] if v["NRF2_crystal"] is not None else "-",
               "%.2f" % v["designed_binder_predicted"]
               if v["designed_binder_predicted"] is not None else "-"))
    log("\n判讀：NRF2 在生物學上與三個 small Maf 都能結合，所以它那一欄應該接近 0。")

    json.dump({"applied": report, "comparison": comp,
               "leucine_controls": {r.label: float(r.ddg) for r in ctrl.itertuples()},
               "all": json.loads(M.round(4).to_json(orient="records"))},
              open(f"{RES}/137_stabddg_control.json", "w"), indent=2, ensure_ascii=False)
    log("\nDONE → results/137_*")


if __name__ == "__main__":
    main()
