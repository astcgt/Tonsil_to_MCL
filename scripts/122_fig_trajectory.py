#!/usr/bin/env python3
"""分化時序圖（圖 6）。三個面板，全部畫自 results/120_*。

  a  正常 B 細胞分化流形（扁桃腺圖譜自己的 UMAP），48 個註釋狀態按分區上色
  b  同一張流形上，兩位病人的 MCL 腫瘤細胞被投影到哪裡（kNN 標籤轉移）
  c  擴散時序軸：48 個狀態按 DPT 中位數排序，MCL 的位置畫成同軸上的密度

風格沿用 scripts/103_figures.py / 111_fig_revisit.py：先寫成帶 CSS 變數的向量 SVG，
再按配色方案各壓平一次。

跑法：env/tonsilmcl/bin/python -I scripts/122_fig_trajectory.py
輸出 results/figs/trajectory.svg 與 trajectory.{light,dark}.svg
"""
import json
import os
import re
import gzip
import csv
import math
import random

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results")
FIGS = os.path.join(RES, "figs")
os.makedirs(FIGS, exist_ok=True)

C1 = "var(--series-1, #0072b2)"
C2 = "var(--series-2, #e69f00)"
C3 = "var(--series-3, #009e73)"
C4 = "var(--series-4, #cc79a7)"
C5 = "var(--series-5, #56b4e9)"
CN = "var(--neutral, #8c8c88)"
TXT = "var(--text-primary, #111111)"
TXT2 = "var(--text-secondary, #4a4a48)"
LINE = "var(--line, #d6d6d2)"

THEME = {
    "light": {"--series-1": "#0072b2", "--series-2": "#e69f00", "--series-3": "#009e73",
              "--series-4": "#cc79a7", "--series-5": "#56b4e9", "--neutral": "#8c8c88",
              "--surface-2": "#ffffff", "--text-primary": "#111111",
              "--text-secondary": "#4a4a48", "--line": "#d6d6d2"},
    "dark": {"--series-1": "#2787c9", "--series-2": "#c08400", "--series-3": "#099c72",
             "--series-4": "#c3729f", "--series-5": "#4e9fcf", "--neutral": "#8a8a84",
             "--surface-2": "#171716", "--text-primary": "#ededea",
             "--text-secondary": "#bcbcb6", "--line": "#343431"},
}

# 分區：五個發育區塊（顏色），不是五個狀態
GROUP = [
    ("Naive / activated", C1, ("NBC", "NBC early activation", "NBC CD229+",
                               "NBC IFN-activated", "Proliferative NBC")),
    ("GC-committed", C5, ("Early GC-commited NBC", "GC-commited NBC", "preGC")),
    ("Germinal centre", C3, ("LZ", "LZ proliferative", "LZ_DZ transition",
                             "LZ_DZ reentry commitment", "DZ_LZ transition",
                             "DZ non proliferative", "DZ cell cycle exit",
                             "DZ early Sphase", "DZ late Sphase", "DZ early G2Mphase",
                             "DZ late G2Mphase", "GC DZ Noproli",
                             "PC committed Light Zone GCBC", "DZ migratory PC precursor")),
    ("Memory", C2, ("ncsMBC", "csMBC", "ncsMBC FCRL4/5+", "csMBC FCRL4/5+",
                    "MBC FCRL5+", "Early MBC", "Precursor MBCs",
                    "Reactivated proliferative MBCs")),
    ("Plasma cell", C4, None),   # 其餘全歸這裡
]
GCOL, GNAME = {}, {}
for gi, (nm, col, members) in enumerate(GROUP):
    if members:
        for m in members:
            GCOL[m] = col
            GNAME[m] = nm


def group_of(state):
    return GCOL.get(state, GROUP[-1][1]), GNAME.get(state, GROUP[-1][0])


def resolve(svg, mode="light"):
    th = THEME[mode]
    return re.sub(r"var\((--[\w-]+), (#\w+)\)", lambda m: th.get(m.group(1), m.group(2)), svg)


def t(x, y, s, size=9.0, fill=None, anchor="start", weight="400"):
    return ('<text x="%.1f" y="%.1f" font-size="%s" font-weight="%s" fill="%s" '
            'text-anchor="%s" font-family="Helvetica,Arial,sans-serif">%s</text>'
            % (x, y, size, weight, fill or TXT2, anchor, s))


def axis(x0, y0, x1, y1, col=None, w=1.0, dash=""):
    return ('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="%s"%s/>'
            % (x0, y0, x1, y1, col or LINE, w, (' stroke-dasharray="%s"' % dash) if dash else ""))


def esc(s):
    return (s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


def read_gz(path):
    with gzip.open(path, "rt") as fh:
        return list(csv.DictReader(fh))


# ----------------------------------------------------------------- data ---
S = json.load(open(os.path.join(RES, "120_summary.json")))
REF = read_gz(os.path.join(RES, "120_ref_cells.csv.gz"))
MCL = [r for r in read_gz(os.path.join(RES, "120_mcl_cells.csv.gz"))
       if r["is_tumor"] == "True"]
ORDER = list(csv.DictReader(open(os.path.join(RES, "120_state_order.csv"))))

XS = [float(r["umap1"]) for r in REF]
YS = [float(r["umap2"]) for r in REF]
X0, X1 = min(XS), max(XS)
Y0, Y1 = min(YS), max(YS)
PAD = 8


def scale(w, h):
    """UMAP → panel, aspect preserved."""
    sx = (w - 2 * PAD) / (X1 - X0)
    sy = (h - 2 * PAD) / (Y1 - Y0)
    s = min(sx, sy)
    ox = PAD + ((w - 2 * PAD) - s * (X1 - X0)) / 2
    oy = PAD + ((h - 2 * PAD) - s * (Y1 - Y0)) / 2

    def f(x, y):
        return (ox + s * (float(x) - X0), h - oy - s * (float(y) - Y0))
    return f


def dots(rows, f, colf, r=0.85, op=0.5, n=7000, seed=0):
    """One dot per cell, as a zero-length subpath with a round cap —
    about a third of the bytes of a <circle>, and the same picture."""
    rnd = random.Random(seed)
    if len(rows) > n:
        rows = rnd.sample(rows, n)
    by = {}
    for row in rows:
        by.setdefault(colf(row), []).append(f(row["umap1"], row["umap2"]))
    out = []
    for c, pts in by.items():
        d = "".join("M%.1f %.1fh0" % (x, y) for x, y in pts)
        out.append('<path fill="none" stroke="%s" stroke-opacity="%.2f" '
                   'stroke-width="%.2f" stroke-linecap="round" d="%s"/>'
                   % (c, op, 2 * r, d))
    return "".join(out)


# --------------------------------------------------------------- panel a ---
PH = 300          # plot height shared by panels a and b
PW = 230


def legend_rows(x, y, items, size=8.0, pitch=10.5):
    o = []
    for i, (col, label, op) in enumerate(items):
        yy = y + i * pitch
        o.append('<circle cx="%.1f" cy="%.1f" r="2.9" fill="%s" fill-opacity="%.2f"/>'
                 % (x + 3, yy - 2.7, col, op))
        o.append(t(x + 9.5, yy, label, size, TXT2))
    return "".join(o)


def panel_a(W=PW, H=364):
    f = scale(W, PH)
    o = [dots(REF, f, lambda r: group_of(r["state"])[0], r=0.9, op=0.45, seed=1)]

    label = {"NBC": ("NBC", "end", 0, 3),
             "ncsMBC": ("ncsMBC", "end", 0, 3),
             "csMBC": ("csMBC", "start", 0, 3),
             "preGC": ("preGC", "end", 0, 3),
             "LZ": ("GC light zone", "start", 0, 3),
             "DZ late Sphase": ("GC dark zone", "end", 0, 3),
             "Mature IgG+ PC": ("plasma cell", "end", 0, -7)}
    for st, (txt, anc, ddx, ddy) in label.items():
        cx, cy = f(*S["ref_state_umap"][st])
        o.append('<circle cx="%.1f" cy="%.1f" r="2.7" fill="%s" stroke="%s" '
                 'stroke-width="1.2"/>'
                 % (cx, cy, group_of(st)[0], "var(--surface-2, #ffffff)"))
        dx = (-5 if anc == "end" else 5) + ddx
        o.append(t(cx + dx, cy + ddy, esc(txt), 8.4, TXT, anc, "600"))

    o.append(legend_rows(0, 18, [(col, nm, 1.0) for nm, col, _ in GROUP]))
    o.append(t(0, H - 12, "46,216 normal tonsil B cells from 6 donors,", 7.6, TXT2))
    o.append(t(0, H - 2, "48 annotated states · the atlas's own UMAP", 7.6, TXT2))
    return "".join(o), W, H


# --------------------------------------------------------------- panel b ---
def panel_b(W=PW, H=364):
    f = scale(W, PH)
    cols = {"M102": C4, "M413": C3}
    o = [dots(REF, f, lambda r: CN, r=0.9, op=0.16, n=7000, seed=1)]
    o.append(dots(MCL, f, lambda r: cols.get(r["donor_id"], C4),
                  r=0.95, op=0.40, n=7000, seed=2))

    cx, cy = f(*S["ref_state_umap"]["ncsMBC"])
    o.append('<circle cx="%.1f" cy="%.1f" r="15" fill="none" stroke="%s" '
             'stroke-width="1.1" stroke-dasharray="3 2.5"/>' % (cx, cy, TXT))
    tt = S["tumour"]
    top = max(tt["knn_state_top"].items(), key=lambda kv: kv[1])
    o.append(t(cx - 19, cy - 4, "ncsMBC", 8.4, TXT, "end", "600"))
    o.append(t(cx - 19, cy + 5.5, "%.0f%% of tumour cells" % (100 * top[1] / tt["n"]),
               7.6, TXT2, "end"))

    items = [(CN, "normal reference", 0.45)]
    for d, col in cols.items():
        n = sum(1 for r in MCL if r["donor_id"] == d)
        items.append((col, "%s · %s tumour cells" % (d, f"{n:,}"), 1.0))
    o.append(legend_rows(0, 18, items))
    o.append(t(0, H - 12, "k = 15 nearest-neighbour projection onto the", 7.6, TXT2))
    o.append(t(0, H - 2, "reference · 92% agreement with step 15's soft match", 7.6, TXT2))
    return "".join(o), W, H


# --------------------------------------------------------------- panel c ---
def panel_c(W=312, H=364):
    x0, w = 74, W - 84
    y0, rowh = 10, 6.0
    rows = ORDER                       # already sorted by dpt_median
    o = []

    # MCL tumour as a band across every state row
    tt = S["tumour"]
    ytop, ybot = y0 - 6, y0 + len(rows) * rowh - 2
    o.append('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="%s" '
             'fill-opacity="0.13"/>'
             % (x0 + w * tt["dpt_q25"], ytop, w * (tt["dpt_q75"] - tt["dpt_q25"]),
                ybot - ytop, TXT))
    o.append(axis(x0 + w * tt["dpt_median"], ytop, x0 + w * tt["dpt_median"], ybot,
                  TXT, 1.0, "2 2"))
    o.append(t(x0 + w * tt["dpt_median"] + 4, ytop + 7, "MCL tumour", 7.8, TXT,
               "start", "700"))
    o.append(t(x0 + w * tt["dpt_median"] + 4, ytop + 15.5,
               "19,975 cells, both patients", 7.0, TXT2))

    for i, r in enumerate(rows):
        st = r["state"]
        col, _ = group_of(st)
        yy = y0 + i * rowh
        q25, q75, med = float(r["dpt_q25"]), float(r["dpt_q75"]), float(r["dpt_median"])
        o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
                 'stroke-width="1.7" stroke-opacity="0.45" stroke-linecap="round"/>'
                 % (x0 + w * q25, yy, x0 + w * q75, yy, col))
        o.append('<circle cx="%.1f" cy="%.1f" r="1.9" fill="%s"/>' % (x0 + w * med, yy, col))
        used = r["used_in_step15"] in ("True", "true")
        nm = st if len(st) <= 22 else st[:21] + "\u2026"
        o.append(t(x0 - 5, yy + 2.0, esc(nm), 5.4, TXT2 if used else CN, "end"))
        if not used:
            o.append(t(x0 - 1.5, yy + 2.0, "\u2715", 5.0, CN, "middle"))

    ybase = y0 + len(rows) * rowh + 2
    o.append(axis(x0, ybase, x0 + w, ybase))
    for v in (0, 0.25, 0.5, 0.75, 1.0):
        o.append(axis(x0 + w * v, ybase, x0 + w * v, ybase + 3))
        o.append(t(x0 + w * v, ybase + 11, "%.2f" % v, 7.2, TXT2, "middle"))
    o.append(t(x0 + w / 2, ybase + 21, "diffusion pseudotime, rooted at NBC",
               7.8, TXT2, "middle"))

    # in-sample normal B, for contrast
    nt = S.get("nontumour")
    yd = ybase + 33
    o.append(t(x0 - 5, yd + 2.4, "in-sample normal B", 7.0, CN, "end"))
    o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
             'stroke-width="2.0" stroke-linecap="round"/>'
             % (x0 + w * nt["dpt_q25"], yd, x0 + w * nt["dpt_q75"], yd, CN))
    o.append('<circle cx="%.1f" cy="%.1f" r="2.1" fill="%s"/>'
             % (x0 + w * nt["dpt_median"], yd, CN))
    o.append(t(x0 + w * nt["dpt_q75"] + 5, yd + 2.4, "819 cells", 7.0, CN))

    o.append(t(0, H - 12, "dot = median, bar = interquartile range.  \u2715 = 2 of 48 states", 7.2, CN))
    o.append(t(0, H - 2, "dropped from step 15 (fewer than 10 cells).", 7.2, CN))
    return "".join(o), W, H


def main():
    gap, head = 16, 17
    items = [panel_a(), panel_b(), panel_c()]
    W = sum(w for _, w, _ in items) + gap * (len(items) - 1)
    H = max(h for _, _, h in items) + head
    out, x = [], 0
    for i, (inner, w, h) in enumerate(items):
        out.append('<text x="%.0f" y="11" font-size="12.5" font-weight="700" fill="%s" '
                   'font-family="Helvetica,Arial,sans-serif">%s</text>'
                   % (x, TXT, chr(97 + i)))
        out.append('<g transform="translate(%.1f,%d)">%s</g>' % (x, head, inner))
        x += w + gap
    svg = ('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" width="%d" '
           'height="%d" role="img">%s</svg>' % (W, H, W, H, "".join(out)))
    with open(os.path.join(FIGS, "trajectory.svg"), "w") as fh:
        fh.write(svg)
    print("trajectory.svg  %d x %d  %.1f kB" % (W, H, len(svg) / 1024))
    for mode in ("light", "dark"):
        flat = resolve(svg, mode)
        flat = flat.replace("<svg ", '<svg style="background:%s" ' % THEME[mode]["--surface-2"], 1)
        with open(os.path.join(FIGS, "trajectory.%s.svg" % mode), "w") as fh:
            fh.write(flat)
        print("trajectory.%s.svg  %.1f kB" % (mode, len(flat) / 1024))


if __name__ == "__main__":
    main()
