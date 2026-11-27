#!/usr/bin/env python
"""步驟 215：不延長模擬，改讀數 —— 用 MM/GBSA 重新評分已經跑完的軌跡。

為什麼要做
----------
步驟 190 的門檻判定在 tier1 就停了，而且停得對：
打亂序列的陰性對照維持 **Q = 0.890 ± 0.092**，真正的設計品 0.875 ± 0.020，
無法區分 → `A5_assay_cannot_discriminate`，tier2 的 67 GPU 小時取消。

但壞掉的是**讀數**，不是模擬。Q（原生接觸比例）是用「預測的姿態」定義的，
所以它問的是「這個預測的姿態會不會散」，不是「這個序列黏不黏」。
共摺疊模型既然把打亂序列也擺成接觸，短程 MD 當然讓它留在那裡。

那 11 條軌跡已經花掉約 14 GPU 小時。換一個對序列敏感的讀數就可以不浪費：
單軌跡 MM/GBSA 的近似交互作用能
    ΔE = E(複合體) − E(鏈 A 單獨) − E(鏈 B 單獨)
全部用 amber14 力場 + GBn2 隱式溶劑，座標直接取自既有軌跡。

事前指定（同一個陰性對照當門檻，不換標準）
------------------------------------------
P1  打亂序列的陰性對照 ΔE 必須明顯比真正的設計品**差**（較不負）。
    不然 MM/GBSA 也沒有鑑別力，這條路就徹底關閉，不再投入。
P2  陽性對照（天然 NRF2–MafG 界面）應該和設計品同級或更好。
P3  只有 P1 過了，才可以去看 MafG vs MafK 的差異 —— 否則那個數字沒有意義。
P4  同一系統的重複之間標準差必須小於「設計品與陰性對照之差」，
    否則差異落在雜訊裡。

⚠️ 單軌跡 MM/GBSA 沒有熵項、沒有顯式溶劑校正、假設單獨鏈不鬆弛。
   它是一個排序啟發式，絕對值不可信。這裡只用它做
   「分不分得出打亂序列」這一件事，而那是它做得到的事。
⚠️ GPU 全部被佔用（本專案自己的工作），所以一律用 CPU 平台。
⚠️ **必須先解折疊，這是第一版的 bug。** 能量要的是真實座標，不是距離，
   所以步驟 190 用的「最小影像距離」在這裡沒有用。第一版沒處理，結果自相矛盾：
   `MAFK` 的交互作用能是 −4.2 ± 0.1 kJ/mol（埋藏 1,500 Å² 的複合體不可能這麼小），
   `MAFG` 三個重複是 −148.8 / −204.8 / −97.3、單一重複內標準差到 ±407。
   原因有兩層：(1) 沒有解折疊；(2) `MAFG`/`MAFK` 的軌跡是在 `enforcePeriodicBox`
   被改成 False **之前**跑的（會折疊），`MAFF`/`CTRLPOS`/`CTRLNEG` 是之後跑的
   （不折疊）—— 於是變成拿折疊的和沒折疊的互相比較。
   這一版對每一幀把 B 鏈剛體移到離 A 鏈最近的週期影像（n ∈ {-2..2}³，
   並在時間上保持連續），對本來就連續的軌跡位移量會是零，所以兩種都安全。

輸出 results/215_*
"""
import os
import json
import glob
import importlib.util

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "8"

import numpy as np
import pandas as pd

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
RES = f"{BASE}/results"
WORK = f"{RES}/190_md"
STRIDE = 10            # 每 10 幀取一幀（10 ns / 200 幀 → 20 個樣本）
LAST_HALF = True       # 只用後半段，前半段還在鬆弛


def log(*a):
    print("[215]", *a, flush=True)


def load_reader():
    spec = importlib.util.spec_from_file_location(
        "md", f"{BASE}/scripts/190_md_interface.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def unwrap_B(P, nA, boxes, md):
    """把 B 鏈剛體移到離 A 鏈最近的週期影像，並在時間上保持連續。

    十二面體的約化三斜向量是斜的，所以搜尋範圍要到 ±2（±1 不夠，
    步驟 190 的第一版就是因此留下一批沒修好的幀）。
    """
    ns = np.array([[i, j, k] for i in (-2, -1, 0, 1, 2)
                   for j in (-2, -1, 0, 1, 2)
                   for k in (-2, -1, 0, 1, 2)], dtype=float)
    out = P.copy()
    comA = P[:, :nA, :].mean(1)
    shift = np.zeros(3)
    prev = None
    nfix = 0
    for t in range(len(P)):
        M = md.box_vectors(boxes[t]) if boxes and t < len(boxes) else None
        if M is None:
            continue
        comB = P[t, nA:, :].mean(0)
        cand = ns @ M
        target = comA[t] if prev is None else prev
        dist = np.linalg.norm(comB + shift + cand - target, axis=1)
        best = cand[int(np.argmin(dist))]
        if np.linalg.norm(best) > 1e-9:
            shift = shift + best
            nfix += 1
        out[t, nA:, :] = P[t, nA:, :] + shift
        prev = comB + shift
    return out, nfix


def build_contexts(solvated):
    from openmm import app, unit, Platform
    import openmm as mm
    pdb = app.PDBFile(solvated)
    mod = app.Modeller(pdb.topology, pdb.positions)
    mod.delete([r for r in pdb.topology.residues()
                if r.name in ("HOH", "NA", "CL")])
    ptop, ppos = mod.topology, mod.positions
    nA = sum(1 for a in ptop.atoms() if a.residue.chain.index == 0)
    nB = sum(1 for a in ptop.atoms() if a.residue.chain.index == 1)
    ff = app.ForceField("amber14-all.xml", "implicit/gbn2.xml")
    plat = Platform.getPlatformByName("CPU")
    ctxs = {}
    for key, drop in (("complex", None), ("A", 1), ("B", 0)):
        m2 = app.Modeller(ptop, ppos)
        if drop is not None:
            m2.delete([r for r in ptop.residues() if r.chain.index == drop])
        sysm = ff.createSystem(m2.topology, nonbondedMethod=app.NoCutoff,
                               constraints=None, rigidWater=False)
        integ = mm.VerletIntegrator(1.0 * unit.femtosecond)
        ctxs[key] = mm.Context(sysm, integ, plat)
    return ctxs, nA, nB


def main():
    md = load_reader()
    from openmm import unit
    systems = sorted(d for d in os.listdir(WORK)
                     if os.path.isdir(f"{WORK}/{d}")
                     and os.path.exists(f"{WORK}/{d}/solvated.pdb"))
    rows = []
    for s in systems:
        traj_dirs = sorted(glob.glob(f"{WORK}/{s}/10ns_rep*"))
        traj_dirs = [t for t in traj_dirs
                     if os.path.exists(f"{t}/done.txt")
                     and os.path.exists(f"{t}/traj.dcd")]
        if not traj_dirs:
            continue
        ctxs, nA, nB = build_contexts(f"{WORK}/{s}/solvated.pdb")
        npro = nA + nB
        log("%-22s 蛋白質原子 %d（A %d + B %d），完成的重複 %d 個"
            % (s, npro, nA, nB, len(traj_dirs)))
        for td in traj_dirs:
            frames, boxes = md.read_dcd(f"{td}/traj.dcd")
            if frames is None or len(frames) < 10:
                continue
            P = frames[:, :npro, :]
            P, nfix = unwrap_B(P, nA, boxes, md)
            # 解折疊之後做一次健全性檢查：質心距離應該落在 0.5–4 nm，
            # 不然就是解折疊失敗，這一條軌跡不准用。
            com = np.linalg.norm(P[:, :nA, :].mean(1) - P[:, nA:, :].mean(1),
                                 axis=1)
            if not (0.3 < float(np.median(com)) < 5.0):
                log("   %s 解折疊後質心距離中位 %.2f nm，不合理 → 跳過"
                    % (os.path.basename(td), float(np.median(com))))
                continue
            if nfix:
                log("   %s 解折疊了 %d 幀（質心距離中位 %.2f nm）"
                    % (os.path.basename(td), nfix, float(np.median(com))))
            lo = len(P) // 2 if LAST_HALF else 0
            sel = range(lo, len(P), STRIDE)
            de = []
            for t in sel:
                xyz = P[t]
                ctxs["complex"].setPositions(xyz * unit.nanometer)
                ec = ctxs["complex"].getState(
                    getEnergy=True).getPotentialEnergy().value_in_unit(
                        unit.kilojoule_per_mole)
                ctxs["A"].setPositions(xyz[:nA] * unit.nanometer)
                ea = ctxs["A"].getState(
                    getEnergy=True).getPotentialEnergy().value_in_unit(
                        unit.kilojoule_per_mole)
                ctxs["B"].setPositions(xyz[nA:] * unit.nanometer)
                eb = ctxs["B"].getState(
                    getEnergy=True).getPotentialEnergy().value_in_unit(
                        unit.kilojoule_per_mole)
                de.append(ec - ea - eb)
            de = np.array(de)
            rows.append(dict(system=s, rep=os.path.basename(td),
                             n_frames_used=len(de),
                             dE_mean_kJ=float(de.mean()),
                             dE_sd_kJ=float(de.std()),
                             dE_mean_kcal=float(de.mean() / 4.184)))
            log("   %-12s %2d 幀｜ΔE = %+9.1f ± %5.1f kJ/mol （%+7.1f kcal/mol）"
                % (os.path.basename(td), len(de), de.mean(), de.std(),
                   de.mean() / 4.184))
    D = pd.DataFrame(rows)
    if D.empty:
        log("沒有可評分的軌跡"); return
    D.to_csv(f"{RES}/215_mmgbsa_per_rep.csv", index=False)
    g = D.groupby("system").dE_mean_kcal
    m, sd, cnt = g.mean(), g.std(), g.count()
    log("\n=== 每個系統（kcal/mol，越負越好）===")
    for k in m.index:
        log("  %-22s %+8.1f ± %5.1f  (n=%d)"
            % (k, m[k], sd[k] if np.isfinite(sd[k]) else 0, cnt[k]))

    out = {"method": "單軌跡 MM/GBSA（amber14 + GBn2，CPU 平台）",
           "stride": STRIDE, "last_half_only": LAST_HALF,
           "per_system_kcal": {k: {"mean": float(m[k]),
                                   "sd": float(sd[k]) if np.isfinite(sd[k]) else None,
                                   "n_reps": int(cnt[k])} for k in m.index}}
    D1 = "MAFG__D1_best"; NEG = "MAFG__CTRLNEG_scram"
    POS = "MAFG__CTRLPOS_NRF2"; MK = "MAFK__D1_best"
    if D1 in m.index and NEG in m.index:
        gap = float(m[NEG] - m[D1])          # 正值＝陰性對照較差（較不負）
        noise = float(max(sd.get(D1, 0) or 0, sd.get(NEG, 0) or 0))
        out["P1_negative_control_is_worse"] = {
            "design_kcal": float(m[D1]), "scrambled_kcal": float(m[NEG]),
            "gap_kcal": gap, "pass": bool(gap > 0)}
        out["P4_gap_exceeds_replicate_noise"] = {
            "gap_kcal": gap, "max_replicate_sd_kcal": noise,
            "pass": bool(gap > noise)}
        log("\nP1 陰性對照 %+.1f vs 設計品 %+.1f → 差距 %+.1f kcal/mol → %s"
            % (m[NEG], m[D1], gap, "過" if gap > 0 else "不過"))
        log("P4 差距 %+.1f vs 重複之間最大標準差 %.1f → %s"
            % (gap, noise, "過" if gap > noise else "不過"))
    if POS in m.index and D1 in m.index:
        out["P2_positive_control"] = {
            "natural_NRF2_kcal": float(m[POS]), "design_kcal": float(m[D1]),
            "pass": bool(m[POS] <= m[D1] + abs(0.25 * m[D1]))}
        log("P2 陽性對照（天然 NRF2–MafG）%+.1f vs 設計品 %+.1f → %s"
            % (m[POS], m[D1], "過" if out["P2_positive_control"]["pass"] else "不過"))
    p1ok = out.get("P1_negative_control_is_worse", {}).get("pass", False)
    p4ok = out.get("P4_gap_exceeds_replicate_noise", {}).get("pass", False)
    if MK in m.index and D1 in m.index:
        out["P3_mafg_vs_mafk"] = {
            "MAFG_kcal": float(m[D1]), "MAFK_kcal": float(m[MK]),
            "difference_kcal": float(m[D1] - m[MK]),
            "interpretable": bool(p1ok and p4ok),
            "note": ("P1/P4 未通過時這個差異沒有意義，僅列出數字"
                     if not (p1ok and p4ok) else "可解讀")}
        log("P3 MafG %+.1f vs MafK %+.1f（差 %+.1f）→ %s"
            % (m[D1], m[MK], m[D1] - m[MK],
               "可解讀" if (p1ok and p4ok) else "**不可解讀**（P1/P4 未通過）"))
    out["verdict"] = ("rescue_succeeded" if (p1ok and p4ok)
                      else "readout_also_cannot_discriminate")
    log("\n判定：%s" % out["verdict"])
    json.dump(out, open(f"{RES}/215_mmgbsa.json", "w"), indent=2,
              ensure_ascii=False)
    log("DONE → results/215_*")


if __name__ == "__main__":
    main()
