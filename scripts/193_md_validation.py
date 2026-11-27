#!/usr/bin/env python
"""步驟 193：MD 讀數的方法學檢驗。

為什麼需要這一步
----------------
步驟 190 原本的分析有一個會把結論整個翻轉的錯誤：`DCDReporter` 預設
`enforcePeriodicBox=True`，OpenMM 會把**每個分子**各自折回週期盒內。
設計品與標的是兩個分子，其中一個跨過盒面就整體位移一個盒向量，
而分析程式用的是普通歐氏距離。症狀是自相矛盾的數字：

    MAFG rep1  質心漂移 +3.55 nm（看起來解離）但 Q_final = 0.862
    三個重複的 Q_min 全部剛好 0.000
    MAFK 三個重複全部「第一幀沒有界面接觸」而被整批丟棄 → need_two_systems

修正後的讀數（三斜最小影像，MIC）會讓所有數字變得一致，但「讓數字變好看」
本身不是證據。這支腳本做四個檢驗，證明修正是對的而不只是方便：

  檢驗一  盒向量矩陣的慣例沒有轉置錯：把 B 鏈平移任意整數盒向量組合
          （a、b、c、-a+2b），MIC 距離必須完全不變。
  檢驗二  MIC 不會無中生有：把 B 鏈剛體平移「非盒向量」的 1.5–5.0 nm
          （40 個隨機方向 × 5 個距離），模擬真正的解離，
          MIC 報出的 Q 必須是 0。
  檢驗三  實際資料沒有別名風險：MIC 只在「真實距離 < 最小垂直寬度的一半」
          時精確。檢查每一幀的質心距離都遠低於這個上限。
  檢驗四  第 0 幀的 Q 必須是 1.000（原生接觸就是在第 0 幀定義的）。

輸出 results/193_md_validation.json
"""
import os
import json
import importlib.util

import numpy as np

BASE = "/raid/jupyter-alz.ee09/Tonsil_to_MCL"
RES = f"{BASE}/results"
WORK = f"{RES}/190_md"
CUTOFF = 0.45          # 與 190 一致
Q_CUT = CUTOFF * 1.5


def log(*a):
    print("[193]", *a, flush=True)


def load_md():
    spec = importlib.util.spec_from_file_location(
        "md", f"{BASE}/scripts/190_md_interface.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def chain_atoms(solvated):
    from openmm import app
    top = app.PDBFile(solvated).topology
    idx, ch = [], []
    for a in top.atoms():
        if a.residue.name in ("HOH", "NA", "CL"):
            continue
        if a.element is None or a.element.symbol == "H":
            continue
        idx.append(a.index)
        ch.append(a.residue.chain.index)
    idx, ch = np.array(idx), np.array(ch)
    cid = sorted(set(ch.tolist()))
    sizes = [sum(1 for r in c.residues()) for c in top.chains()]
    return idx, np.where(ch == cid[0])[0], np.where(ch == cid[1])[0], sizes


def perp_widths(M):
    """三斜盒的三個垂直寬度（MIC 的有效範圍是其中最小值的一半）。"""
    V = abs(np.linalg.det(M))
    return np.array([V / np.linalg.norm(np.cross(M[(i + 1) % 3], M[(i + 2) % 3]))
                     for i in range(3)])


def main():
    md = load_md()
    rng = np.random.default_rng(0)
    out = {"cutoff_nm": CUTOFF, "q_cut_nm": Q_CUT, "systems": {}}
    ok_all = True

    systems = sorted(d for d in os.listdir(WORK)
                     if os.path.isdir(f"{WORK}/{d}")
                     and os.path.exists(f"{WORK}/{d}/solvated.pdb"))
    for s in systems:
        traj = None
        for cand in sorted(os.listdir(f"{WORK}/{s}")):
            p = f"{WORK}/{s}/{cand}/traj.dcd"
            if os.path.exists(p):
                traj = p
                break
        if traj is None:
            continue
        idx, A, B, sizes = chain_atoms(f"{WORK}/{s}/solvated.pdb")
        frames, boxes = md.read_dcd(traj)
        if frames is None:
            continue
        P = frames[:, idx, :]
        M = md.box_vectors(boxes[0])
        Mi = np.linalg.inv(M)
        w = perp_widths(M)
        lim = float(w.min() / 2)

        d0 = md.mic_dist_mat(P[0][A], P[0][B], M, Mi)
        native = d0 < CUTOFF
        nn = int(native.sum())
        q0 = float((d0[native] < Q_CUT).mean())

        # 檢驗一：整數盒向量平移不變性。
        # 只在 MIC 聲稱精確的範圍內檢查，也就是距離 < 最小垂直寬度的一半。
        # round 式 MIC 對超出這個範圍的原子對本來就不保證找到真正的最近影像，
        # 所以把那些對一起檢查會測到一個沒有人依賴的性質。
        # （實測 MAFG__CTRLPOS_NRF2：724,686 對裡只有 1 對不變，距離 5.46 nm，
        #  高於上限 4.67 nm，且不屬於原生接觸集。）
        within = d0 < lim
        shifts = {"+a": M[0], "+b": M[1], "+c": M[2],
                  "-a+2b": -M[0] + 2 * M[1], "+a-b+c": M[0] - M[1] + M[2]}
        inv, inv_beyond = {}, {}
        for k, v in shifts.items():
            dm = md.mic_dist_mat(P[0][A], P[0][B] + v, M, Mi)
            dev = np.abs(d0 - dm)
            inv[k] = float(dev[within].max())
            inv_beyond[k] = float(dev[~within].max()) if (~within).any() else 0.0
        t1 = max(inv.values()) < 1e-5
        # 原生接觸集必須完全落在有效範圍內（否則 Q 本身就不可信）
        t1b = bool((d0[native] < lim).all())
        t1 = t1 and t1b

        # 檢驗二：非盒向量的真實分離必須報 Q = 0
        fab = {}
        for mag in (1.5, 2.0, 3.0, 4.0, 5.0):
            qs = []
            for _ in range(40):
                v = rng.normal(size=3)
                v = v / np.linalg.norm(v) * mag
                dm = md.mic_dist_mat(P[0][A], P[0][B] + v, M, Mi)
                qs.append(float((dm[native] < Q_CUT).mean()))
            fab["%.1fnm" % mag] = {"median": float(np.median(qs)),
                                   "max": float(np.max(qs))}
        t2 = all(v["max"] == 0.0 for v in fab.values())

        # 檢驗三：實際質心距離 vs MIC 上限
        step = max(1, len(P) // 20)
        coms = []
        for t in range(0, len(P), step):
            Mt = md.box_vectors(boxes[t]) if boxes and t < len(boxes) else M
            coms.append(md.mic_vec(P[t][A].mean(0), P[t][B].mean(0),
                                   Mt, np.linalg.inv(Mt)))
        t3 = max(coms) < lim

        t4 = abs(q0 - 1.0) < 1e-9
        ok = bool(t1 and t2 and t3 and t4)
        ok_all = ok_all and ok
        out["systems"][s] = {
            "trajectory": os.path.relpath(traj, BASE),
            "chain_residues": sizes[:2],
            "n_heavy_atoms": [int(len(A)), int(len(B))],
            "box_perp_widths_nm": [float(x) for x in w],
            "mic_valid_below_nm": lim,
            "n_native_contacts": nn,
            "test1_lattice_invariance_max_dev_nm_within_limit": inv,
            "test1_max_dev_nm_beyond_limit": inv_beyond,
            "test1_n_pairs_beyond_limit": int((~within).sum()),
            "test1_n_pairs_total": int(within.size),
            "test1_all_native_contacts_within_limit": t1b,
            "test1_pass": bool(t1),
            "test2_fabricated_Q_after_rigid_separation": fab,
            "test2_pass": bool(t2),
            "test3_com_range_nm": [float(min(coms)), float(max(coms))],
            "test3_pass": bool(t3),
            "test4_frame0_Q": q0,
            "test4_pass": bool(t4),
            "all_pass": ok,
        }
        log("%-22s 原生接觸 %4d｜檢驗一 %s（有效範圍內最大偏差 %.1e nm，"
            "範圍外 %d/%d 對）｜檢驗二 %s｜"
            "檢驗三 %s（質心 %.2f–%.2f < %.2f）｜檢驗四 Q0=%.3f %s"
            % (s, nn, "過" if t1 else "不過", max(inv.values()),
               int((~within).sum()), int(within.size),
               "過" if t2 else "不過", "過" if t3 else "不過",
               min(coms), max(coms), lim, q0, "過" if t4 else "不過"))

    out["all_systems_pass"] = bool(ok_all)
    json.dump(out, open(f"{RES}/193_md_validation.json", "w"),
              indent=2, ensure_ascii=False)
    log("\n全部通過：%s" % ok_all)
    log("DONE → results/193_md_validation.json")
    if not ok_all:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
