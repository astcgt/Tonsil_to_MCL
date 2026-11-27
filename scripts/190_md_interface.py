#!/usr/bin/env python
"""步驟 190（B1）：分子動力學 —— 設計品在 MafG 上待得住，在 MafK 上待不住嗎？

為什麼要做
----------
共摺疊模型只回答「模型能不能把這個界面建出來」。MD 回答的是不同的問題：
**這個界面在水裡、在室溫下、隨時間推移，撐不撐得住。**
目前四個計算方法對「D1 對 MafG 有沒有選擇性」意見分歧
（Boltz 說沒有但陰性對照失效；Chai-1 說 MafG 高 0.051；StaB-ddG 說差很多）。
MD 是唯一一個不看任何「信心分數」的判準。

分階段，每階有明確的續做門檻
----------------------------
tier1（最小）：MAFG–D1、MAFK–D1，各 3 個重複 × 20 ns，顯溶劑
  讀出：原生接觸保留率 Q、界面 RMSD、埋藏面積、兩鏈質心距離
  ★ 門檻 A（是否續做 tier2）：
     A1 兩個系統的平均 Q 差距 > 0.10 且 3 個重複方向一致  → 有訊號，續做
     A2 兩個系統都穩定（Q > 0.8）且差距 < 0.05            → 無差異，續做一次確認
     A3 重複之間互相矛盾（Q 的標準差 > 0.15）             → 取樣不足，續做
     A4 兩個都在 20 ns 內解離（Q < 0.4）                  → 結構本身有問題，停止並回報
tier2（中等）：延長到 100 ns，加入 MAFF–D1 與晶體 MafG–NRF2 當參考
  ★ 門檻 B：tier1 的差距在 100 ns 是否維持且 p < 0.05（重複間 t 檢定）
tier3（完整）：5 個重複 × 200 ns，加上 steered MD 算解離功

跑法：
  scripts/190_md_interface.py --tier 1 --gpu 4
  scripts/190_md_interface.py --analyze           # 只分析並印出門檻判定
"""
import os
import sys
import json
import glob
import time
import argparse
import warnings

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
RES = f"{BASE}/results"
WORK = f"{RES}/190_md"
PRED132 = f"{RES}/132_specificity/boltz_results_inputs/predictions"

TIERS = {
    # tier1 擴成完整五系統：先用最便宜的 10 ns 問「這個讀數分辨得出東西嗎」。
    # 陰性對照（打亂序列）若在 10 ns 也維持高 Q，就代表 Q 沒有鑑別力，
    # 不值得再花 67 GPU 小時做 30 ns。
    1: dict(systems=["MAFG__D1_best", "MAFK__D1_best", "MAFF__D1_best",
                     "MAFG__CTRLPOS_NRF2", "MAFG__CTRLNEG_scram"],
            reps=3, ns=10),
    # tier2 的長度是依實測速度（OpenCL 150 ns/day、192k 原子）回推的：
    # 4 系統 × 3 重複 × 30 ns = 360 ns ≈ 2.4 GPU-day ≈ 兩張卡並行 1.2 天
    2: dict(systems=["MAFG__D1_best", "MAFK__D1_best", "MAFF__D1_best",
                     "MAFG__CTRLPOS_NRF2", "MAFG__CTRLNEG_scram"],
            reps=3, ns=30),
    3: dict(systems=["MAFG__D1_best", "MAFK__D1_best", "MAFF__D1_best",
                     "MAFG__CTRLPOS_NRF2", "MAFG__CTRLNEG_scram"], reps=5, ns=200),
}
DT_FS = 4.0            # 4 fs with hydrogen mass repartitioning
REPORT_PS = 50.0       # 每 50 ps 存一次


def log(*a):
    print("[190]", *a, flush=True)


def cif_to_pdb(cif, out):
    from Bio.PDB import MMCIFParser, PDBIO
    st = MMCIFParser(QUIET=True).get_structure("x", cif)
    io = PDBIO(); io.set_structure(st); io.save(out)


def build(system):
    """把 Boltz 預測的複合體做成可以跑 MD 的溶劑化系統。"""
    from openmm import app, unit, Platform
    import openmm as mm
    src = sorted(glob.glob(f"{PRED132}/**/{system}_model_0.cif", recursive=True))
    if not src:
        src = sorted(glob.glob(f"{PRED132}/**/{system}_model_*.cif", recursive=True))
    if not src:
        raise SystemExit(f"找不到 {system} 的預測結構")
    os.makedirs(f"{WORK}/{system}", exist_ok=True)
    raw = f"{WORK}/{system}/raw.pdb"
    cif_to_pdb(src[0], raw)

    # 預測結構沒有 C 端 OXT，力場模板對不上 → 用 PDBFixer 補齊端基與缺失原子
    from pdbfixer import PDBFixer
    fx = PDBFixer(filename=raw)
    fx.findMissingResidues()
    fx.missingResidues = {}          # 不補整段缺失殘基，只補原子
    fx.findNonstandardResidues()
    fx.replaceNonstandardResidues()
    fx.removeHeterogens(keepWater=False)
    fx.findMissingAtoms()
    fx.addMissingAtoms()
    fx.addMissingHydrogens(7.4)
    fixed = f"{WORK}/{system}/fixed.pdb"
    app.PDBFile.writeFile(fx.topology, fx.positions, open(fixed, "w"))

    pdb = app.PDBFile(fixed)
    ff = app.ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
    mod = app.Modeller(pdb.topology, pdb.positions)
    mod.deleteWater()
    # 延長型捲曲螺旋在立方盒裡大半是水 → 十二面體盒，原子數少約 30%
    try:
        mod.addSolvent(ff, model="tip3p", padding=0.9 * unit.nanometer,
                       boxShape="dodecahedron",
                       ionicStrength=0.15 * unit.molar, neutralize=True)
    except TypeError:
        mod.addSolvent(ff, model="tip3p", padding=0.9 * unit.nanometer,
                       ionicStrength=0.15 * unit.molar, neutralize=True)
    sysm = ff.createSystem(mod.topology, nonbondedMethod=app.PME,
                           nonbondedCutoff=1.0 * unit.nanometer,
                           constraints=app.HBonds, rigidWater=True,
                           hydrogenMass=4.0 * unit.amu)
    with open(f"{WORK}/{system}/system.xml", "w") as fh:
        fh.write(mm.XmlSerializer.serialize(sysm))
    app.PDBFile.writeFile(mod.topology, mod.positions,
                          open(f"{WORK}/{system}/solvated.pdb", "w"))
    n = mod.topology.getNumAtoms()
    log("  %s：溶劑化後 %d 個原子" % (system, n))
    return n


def run_rep(system, rep, ns, gpu):
    from openmm import app, unit, Platform
    import openmm as mm
    d = f"{WORK}/{system}"
    out = f"{d}/{ns}ns_rep{rep}"
    os.makedirs(out, exist_ok=True)
    if os.path.exists(f"{out}/traj.dcd") and os.path.exists(f"{out}/done.txt"):
        log("  %s rep%d 已完成，跳過" % (system, rep)); return
    pdb = app.PDBFile(f"{d}/solvated.pdb")
    sysm = mm.XmlSerializer.deserialize(open(f"{d}/system.xml").read())
    integ = mm.LangevinMiddleIntegrator(310 * unit.kelvin, 1.0 / unit.picosecond,
                                        DT_FS * unit.femtoseconds)
    integ.setRandomNumberSeed(1000 + rep)
    sysm.addForce(mm.MonteCarloBarostat(1 * unit.bar, 310 * unit.kelvin, 25))
    # 這台機器的驅動是 CUDA 12.0，而 openmm-cuda-12 的 wheel 用更新的 CUDA 編譯，
    # 載入 PTX 會失敗（CUDA_ERROR_UNSUPPORTED_PTX_VERSION）。
    # OpenCL 在 A100 上實測 173 ns/day，夠用，而且不必賭 wheel 版本。
    sim = None
    for pname in ("CUDA", "OpenCL"):
        try:
            plat = Platform.getPlatformByName(pname)
            sim = app.Simulation(pdb.topology, sysm, integ, plat,
                                 {"Precision": "mixed", "DeviceIndex": "0"})
            sim.context.setPositions(pdb.positions)
            sim.context.getState(getEnergy=True)      # 這一步才會真的編譯 kernel
            log("  平台:", pname)
            break
        except Exception as e:
            log("  %s 不可用（%s），改試下一個" % (pname, type(e).__name__))
            sim = None
    if sim is None:
        raise SystemExit("沒有可用的 GPU 平台")
    t0 = time.time()
    sim.minimizeEnergy(maxIterations=5000)
    sim.context.setVelocitiesToTemperature(310 * unit.kelvin, 2000 + rep)
    sim.step(int(200 * 1000 / DT_FS))           # 200 ps 平衡
    nsteps = int(ns * 1e6 / DT_FS)
    every = int(REPORT_PS * 1000 / DT_FS)
    # enforcePeriodicBox=False 很關鍵：預設是 True，OpenMM 會把「每個分子」
    # 各自折回盒內。設計品與標的是兩個分子，其中一個跨過盒面就整體位移一個
    # 盒向量，事後算距離會看到 3–4 nm 的假解離。關掉之後座標在時間上連續。
    sim.reporters.append(app.DCDReporter(f"{out}/traj.dcd", every,
                                        enforcePeriodicBox=False))
    sim.reporters.append(app.StateDataReporter(
        f"{out}/log.csv", every * 4, step=True, time=True, temperature=True,
        potentialEnergy=True, speed=True))
    sim.step(nsteps)
    app.PDBFile.writeFile(sim.topology,
                          sim.context.getState(getPositions=True).getPositions(),
                          open(f"{out}/final.pdb", "w"))
    open(f"{out}/done.txt", "w").write("%s %d ns %.1f min\n"
                                       % (system, ns, (time.time() - t0) / 60))
    log("  %s rep%d：%d ns 完成，耗時 %.1f 分鐘"
        % (system, rep, ns, (time.time() - t0) / 60))


# ----------------------------------------------------------------- 分析 ---
def analyse_openmm(systems, ns, cutoff=0.45):
    """不依賴 mdtraj：用 OpenMM 的 DCDFile 讀軌跡，自己算接觸。"""
    from openmm import app, unit
    import openmm as mm
    rows = []
    for system in systems:
        d = f"{WORK}/{system}"
        if not os.path.exists(f"{d}/solvated.pdb"):
            continue
        top = app.PDBFile(f"{d}/solvated.pdb").topology
        # 只取蛋白質重原子，並記錄鏈
        idx, ch = [], []
        for atom in top.atoms():
            if atom.residue.name in ("HOH", "NA", "CL"):
                continue
            if atom.element is None or atom.element.symbol == "H":
                continue
            idx.append(atom.index); ch.append(atom.residue.chain.index)
        idx = np.array(idx); ch = np.array(ch)
        cid = sorted(set(ch.tolist()))
        if len(cid) < 2:
            log("  %s：只有一條鏈，跳過" % system); continue
        A = np.where(ch == cid[0])[0]; B = np.where(ch == cid[1])[0]
        for rep_dir in sorted(glob.glob(f"{d}/{ns}ns_rep*")):
            dcd = f"{rep_dir}/traj.dcd"
            if not os.path.exists(dcd):
                continue
            # OpenMM 沒有讀 DCD 的公開 API，所以自己解析（見 read_dcd）
            frames, boxes = read_dcd(dcd)
            if frames is None or len(frames) < 3:
                continue
            P = frames[:, idx, :]            # (T, natoms_prot, 3) in nm
            Ms = [box_vectors(boxes[t]) if boxes and t < len(boxes) else None
                  for t in range(len(P))]
            Mi = [None if M is None else np.linalg.inv(M) for M in Ms]
            n_box = sum(M is not None for M in Ms)
            if n_box == 0:
                log("  %s %s：軌跡沒有晶胞資訊，距離未做週期修正"
                    % (system, os.path.basename(rep_dir)))
            d0 = mic_dist_mat(P[0][A], P[0][B], Ms[0], Mi[0])
            native = d0 < cutoff
            nn = int(native.sum())
            if nn == 0:
                log("  %s %s：第一幀沒有界面接觸" % (system, os.path.basename(rep_dir)))
                continue
            q, com = [], []
            for t in range(len(P)):
                dm = mic_dist_mat(P[t][A], P[t][B], Ms[t], Mi[t])
                q.append(float((dm[native] < cutoff * 1.5).mean()))
                com.append(mic_vec(P[t][A].mean(0), P[t][B].mean(0),
                                   Ms[t], Mi[t]))
            q = np.array(q); com = np.array(com)
            half = len(q) // 2
            rows.append(dict(system=system, rep=os.path.basename(rep_dir),
                             n_frames=len(q), n_native_contacts=nn,
                             Q_mean_last_half=float(q[half:].mean()),
                             Q_min=float(q.min()), Q_final=float(q[-1]),
                             com_start=float(com[0]),
                             com_mean_last_half=float(com[half:].mean()),
                             com_drift=float(com[half:].mean() - com[0])))
    return pd.DataFrame(rows)


def read_dcd(path):
    """最小的 DCD 讀取器（OpenMM 寫出的格式），回傳 (T, N, 3) 的 nm 座標。"""
    import struct
    with open(path, "rb") as fh:
        raw = fh.read()
    off = 0

    def blk():
        nonlocal off
        n = struct.unpack_from("<i", raw, off)[0]; off += 4
        b = raw[off:off + n]; off += n
        off += 4
        return b
    hdr = blk()
    nset = struct.unpack_from("<i", hdr, 4)[0]
    has_box = struct.unpack_from("<i", hdr, 44)[0]
    blk()                                   # title
    natom = struct.unpack_from("<i", blk(), 0)[0]
    out, boxes = [], []
    while off < len(raw):
        try:
            bx = None
            if has_box:
                b = blk()
                # OpenMM 寫的是 6 個 double：
                #   (a, sin(pi/2-gamma), b, sin(pi/2-beta), sin(pi/2-alpha), c)
                # 也就是 (a, cos gamma, b, cos beta, cos alpha, c)，長度單位 Å。
                bx = struct.unpack_from("<6d", b, 0)
            x = np.frombuffer(blk(), dtype="<f4")
            y = np.frombuffer(blk(), dtype="<f4")
            z = np.frombuffer(blk(), dtype="<f4")
        except Exception:
            break
        if len(x) != natom:
            break
        out.append(np.stack([x, y, z], 1))
        boxes.append(bx)
    if not out:
        return None, None
    return np.asarray(out) / 10.0, boxes    # Å → nm


def box_vectors(bx):
    """把 DCD 的 6 個數字還原成 nm 的 3×3 盒向量矩陣。"""
    import math
    from openmm.app.internal.unitcell import computePeriodicBoxVectors
    if bx is None:
        return None
    a, cg, b, cb, ca, c = bx
    try:
        v = computePeriodicBoxVectors(
            a / 10.0, b / 10.0, c / 10.0,
            math.acos(max(-1.0, min(1.0, ca))),
            math.acos(max(-1.0, min(1.0, cb))),
            math.acos(max(-1.0, min(1.0, cg))))
    except Exception:
        return None
    try:                                   # 可能帶單位，也可能已經是純數
        from openmm import unit as _u
        v = v.value_in_unit(_u.nanometer)
    except Exception:
        pass
    return np.array([[float(v[i][j]) for j in range(3)] for i in range(3)],
                    dtype=float)


def mic_dist_mat(a, b, M, Minv):
    """三斜盒的最小影像距離矩陣。

    前一版用「剛體解折疊」：假設 OpenMM 逐分子折回盒內後，整條鏈只差一個
    整數盒向量組合，在 n∈{-1,0,1}³ 裡搜。十二面體盒的約化三斜向量是斜的，
    ±1 不夠，於是留下一批沒修好的幀（症狀：質心相距 4 nm 但 Q_final 還有 0.86）。
    改成每一組原子對各自取最小影像，就不需要假設折疊方式：
    把位移換到分數座標、減掉最接近的整數、再換回來。
    對接觸距離（< 0.7 nm）這在任何折疊狀態下都是精確的。
    代價是質心距離會在「半個盒寬」飽和 —— 判斷結合/解離足夠，
    但不能當成真實位移量來報告。
    """
    d = b[None, :, :] - a[:, None, :]
    if M is None:
        return np.sqrt((d ** 2).sum(-1))
    f = d @ Minv
    f -= np.round(f)
    d = f @ M
    return np.sqrt((d ** 2).sum(-1))


def mic_vec(va, vb, M, Minv):
    """單一向量的最小影像距離（給質心用）。"""
    d = vb - va
    if M is None:
        return float(np.linalg.norm(d))
    f = d @ Minv
    f -= np.round(f)
    return float(np.linalg.norm(f @ M))


def gate(df, ns=0):
    """續做門檻判定。

    主比較固定是「同一個設計品在 MafG 與 MafK 上的 Q 差距」。
    之前用 sorted(index)[:2] 取前兩個系統，tier1 只有兩個系統時剛好正確，
    但 tier2 起有對照系統，字母序會把 MAFF 和 MAFG__CTRLNEG 配成一對，
    算出一個沒有意義的 ΔQ。這裡改成明確指定。
    對照系統只做自己的合格檢查，不進入主比較、也不參與解離判定。
    """
    if df.empty:
        return {"verdict": "no_data"}
    g = df.groupby("system").Q_mean_last_half
    m, s = g.mean(), g.std()
    sysm = sorted(m.index)
    v = {"per_system": {k: {"Q_mean": float(m[k]), "Q_sd": float(s[k]),
                            "n_reps": int(g.count()[k])} for k in sysm}}

    # 對照的事前預期：陽性（天然 NRF2–MafG 界面）要穩住，陰性（打亂序列）要散掉。
    ctl = {}
    for key, name, expect, ok in [
            ("positive", "MAFG__CTRLPOS_NRF2", "Q > 0.60", lambda q: q > 0.60),
            ("negative", "MAFG__CTRLNEG_scram", "Q < 0.40", lambda q: q < 0.40)]:
        if name in m.index:
            q = float(m[name])
            ctl[key] = {"system": name, "Q_mean": q, "expect": expect,
                        "pass": bool(ok(q))}
    if ctl:
        v["controls"] = ctl
        v["controls_all_pass"] = all(c["pass"] for c in ctl.values())

    # 對照沒有鑑別力就停：陰性對照該散掉卻沒散，代表這個 Q 讀數
    # 連 Boltz 的偽陽性都拒絕不了，再延長時間也不會變得有意義。
    if "negative" in ctl and not ctl["negative"]["pass"]:
        v["verdict"] = "A5_assay_cannot_discriminate"
        v["action"] = ("停止：打亂序列的陰性對照在 %d ns 仍維持 Q=%.3f（預期 < 0.40），"
                       "這個讀數沒有鑑別力，不值得延長"
                       % (ns, ctl["negative"]["Q_mean"]))
        return v

    MAIN = ["MAFG__D1_best", "MAFK__D1_best"]
    if all(k in m.index for k in MAIN):
        a, b = MAIN
    elif len(sysm) >= 2:
        a, b = sysm[0], sysm[1]
    else:
        v["verdict"] = "need_two_systems"
        return v

    mm, ss = m[[a, b]], s[[a, b]]
    diff = float(mm[a] - mm[b])
    v["difference"] = {"pair": [a, b], "delta_Q": diff}
    maxsd = float(max(ss.fillna(0)))
    if (mm < 0.4).all():
        v["verdict"] = "A4_both_dissociated"
        v["action"] = "停止：兩個系統都在 tier1 內解離，先檢查起始結構"
    elif maxsd > 0.15:
        v["verdict"] = "A3_replicates_disagree"
        v["action"] = "續做 tier2：重複之間矛盾，取樣不足"
    elif abs(diff) > 0.10:
        v["verdict"] = "A1_signal"
        v["action"] = "續做 tier2：兩系統差距 > 0.10 且重複一致"
    elif (mm > 0.8).all() and abs(diff) < 0.05:
        v["verdict"] = "A2_both_stable_no_difference"
        v["action"] = "續做 tier2 一次確認（延長到 100 ns）"
    else:
        v["verdict"] = "ambiguous"
        v["action"] = "續做 tier2"
    return v


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tier", type=int, default=1)
    ap.add_argument("--gpu", default="4")
    ap.add_argument("--analyze", action="store_true")
    ap.add_argument("--only", default="", help="只跑這個系統")
    ap.add_argument("--rep", type=int, default=0,
                    help="只跑這個重複（0＝全部）。排程以 系統×重複 為單位配卡時用。")
    a = ap.parse_args()
    os.makedirs(WORK, exist_ok=True)
    cfg = TIERS[a.tier]

    if not a.analyze:
        os.environ["CUDA_VISIBLE_DEVICES"] = a.gpu
        log("tier %d：%d 個系統 × %d 重複 × %d ns（GPU %s）"
            % (a.tier, len(cfg["systems"]), cfg["reps"], cfg["ns"], a.gpu))
        todo = [a.only] if a.only else cfg["systems"]
        for s_ in todo:
            if not os.path.exists(f"{WORK}/{s_}/system.xml"):
                log("建立系統 %s …" % s_)
                build(s_)
        reps = [a.rep] if a.rep else list(range(1, cfg["reps"] + 1))
        for s_ in todo:
            for r in reps:
                run_rep(s_, r, cfg["ns"], a.gpu)
        if a.rep:                       # 單一 rep 的工作不做分析，交給排程最後統一跑
            log("rep 模式完成，略過分析")
            return

    log("\n=== 分析 ===")
    df = analyse_openmm(cfg["systems"], cfg["ns"])
    if df.empty:
        log("沒有可分析的軌跡"); return
    df.to_csv(f"{RES}/190_md_tier{a.tier}_per_rep.csv", index=False)
    pd.set_option("display.width", 220)
    log("\n" + df.round(3).to_string(index=False))
    v = gate(df, cfg["ns"])
    log("\n=== 門檻判定 ===")
    for k, d in v.get("per_system", {}).items():
        log("  %-22s Q(後半段) = %.3f ± %.3f  (n=%d)"
            % (k, d["Q_mean"], d["Q_sd"] if np.isfinite(d["Q_sd"]) else 0, d["n_reps"]))
    for key, c in v.get("controls", {}).items():
        log("  對照（%s）%-22s Q = %.3f  預期 %s  →  %s"
            % ("陽性" if key == "positive" else "陰性", c["system"],
               c["Q_mean"], c["expect"], "通過" if c["pass"] else "未通過"))
    if "difference" in v:
        log("  差距 ΔQ = %+.3f  (%s − %s)"
            % (v["difference"]["delta_Q"], *v["difference"]["pair"]))
    log("  判定：%s" % v["verdict"])
    log("  下一步：%s" % v.get("action", "-"))
    json.dump({"tier": a.tier, "config": cfg, "gate": v,
               "per_rep": json.loads(df.round(4).to_json(orient="records"))},
              open(f"{RES}/190_md_tier{a.tier}.json", "w"),
              indent=2, ensure_ascii=False)
    log("\nDONE → results/190_md_tier%d*" % a.tier)


if __name__ == "__main__":
    main()
