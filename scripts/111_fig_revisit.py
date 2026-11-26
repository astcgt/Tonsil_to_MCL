#!/usr/bin/env python3
"""Figure for Section 1: reproducing and re-examining Massoni-Badosa et al. Fig. 7G-I.

Four panels, all drawn from results/110_revisit_fig7.json:
  a  chromosome Y module score in tumour cells, per patient      (their Fig. 7G, top)
  b  inferCNV score: tonsil reference B, in-sample normal B, tumour (their Fig. 7G, bottom)
  c  metallothionein score per tumour cluster, both patients      (their Fig. 7I)
  d  our quantitative re-examination: which normal state each tumour cell matches,
     and how concentrated that match is

Style follows scripts/103_figures.py: vector SVG written against the page's own
colour variables, then flattened once per colour scheme.

Run:  env/tonsilmcl/bin/python -I scripts/111_fig_revisit.py
Writes results/figs/revisit.svg and revisit.{light,dark}.svg
"""
import json
import os
import re

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIGS = os.path.join(BASE, "results", "figs")
os.makedirs(FIGS, exist_ok=True)

D = json.load(open(os.path.join(BASE, "results", "110_revisit_fig7.json")))
F = D["fig"]

C1 = "var(--series-1, #0072b2)"
C2 = "var(--series-2, #e69f00)"
C3 = "var(--series-3, #009e73)"
C4 = "var(--series-4, #cc79a7)"
CN = "var(--neutral, #8c8c88)"
TXT = "var(--text-primary, #111111)"
TXT2 = "var(--text-secondary, #4a4a48)"
LINE = "var(--line, #d6d6d2)"

THEME = {
    "light": {"--series-1": "#0072b2", "--series-2": "#e69f00", "--series-3": "#009e73",
              "--series-4": "#cc79a7", "--neutral": "#8c8c88", "--surface-2": "#ffffff",
              "--text-primary": "#111111", "--text-secondary": "#4a4a48", "--line": "#d6d6d2"},
    "dark": {"--series-1": "#2787c9", "--series-2": "#c08400", "--series-3": "#099c72",
             "--series-4": "#c3729f", "--neutral": "#8a8a84", "--surface-2": "#171716",
             "--text-primary": "#ededea", "--text-secondary": "#bcbcb6", "--line": "#343431"},
}


def resolve(svg, mode="light"):
    th = THEME[mode]
    return re.sub(r"var\((--[\w-]+), (#\w+)\)", lambda m: th.get(m.group(1), m.group(2)), svg)


def t(x, y, s, size=9.5, fill=None, anchor="start", weight="400"):
    return ('<text x="%.1f" y="%.1f" font-size="%s" font-weight="%s" fill="%s" '
            'text-anchor="%s" font-family="Helvetica,Arial,sans-serif">%s</text>'
            % (x, y, size, weight, fill or TXT2, anchor, s))


def axis(x0, y0, x1, y1, col=None, w=1.0, dash=""):
    return ('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="%s"%s/>'
            % (x0, y0, x1, y1, col or LINE, w, (' stroke-dasharray="%s"' % dash) if dash else ""))


def smooth(c, k=2):
    """Light moving average so a 48-bin histogram reads as a density curve."""
    n = len(c)
    return [sum(c[max(0, i - k):min(n, i + k + 1)]) / len(c[max(0, i - k):min(n, i + k + 1)])
            for i in range(n)]


def density_path(h, x0, w, ybase, hgt, k=2, close=True):
    c = smooth(h["counts"], k)
    mx = max(c) or 1
    pts = []
    for i, v in enumerate(c):
        xx = x0 + w * (i + 0.5) / len(c)
        yy = ybase - hgt * v / mx
        pts.append((xx, yy))
    d = "M%.1f %.1f " % (x0, ybase) if close else "M%.1f %.1f " % pts[0]
    d += " ".join("L%.1f %.1f" % p for p in pts)
    if close:
        d += " L%.1f %.1f Z" % (x0 + w, ybase)
    return d


def xpos(h, val, x0, w):
    return x0 + w * (val - h["lo"]) / (h["hi"] - h["lo"])


# --------------------------------------------------------------- panel a ---
def panel_a(W=252, H=250):
    o = []
    o.append(t(0, 11, "chrY module score, tumour cells", 10, TXT, weight="600"))
    o.append(t(0, 24, "their Fig. 7G: one major clone has lost chrY", 9, TXT2))
    x0, w = 30, W - 44
    rows = [("M102", C1), ("M413", C2)]
    for i, (d, col) in enumerate(rows):
        h = F["chrY"][d]
        yb = 44 + i * 100 + 76
        o.append(axis(x0, yb, x0 + w, yb))
        o.append('<path d="%s" fill="%s" fill-opacity="0.26" stroke="%s" stroke-width="1.5"/>'
                 % (density_path(h, x0, w, yb, 68), col, col))
        # the two component means, and the split
        for val, lab in ((h["mean_loss"], "chrY−"), (h["mean_ret"], "chrY+")):
            xx = xpos(h, val, x0, w)
            o.append(axis(xx, yb - 72, xx, yb, TXT2, 0.9, "2 2"))
            o.append(t(xx, yb - 76, lab, 8.5, TXT2, "middle"))
        o.append(t(x0, 44 + i * 100 + 10, d, 9.5, TXT, weight="600"))
        o.append(t(x0 + 26, 44 + i * 100 + 10,
                   "%.0f%% chrY−  ·  d = %.2f" % (h["frac_loss"] * 100,
                                                       D["chrY"][d]["cohens_d"]), 9, TXT2))
        for v in (-1, 0, 1):
            xx = xpos(h, v, x0, w)
            o.append(axis(xx, yb, xx, yb + 3, LINE))
            o.append(t(xx, yb + 13, "%g" % v, 8.5, TXT2, "middle"))
    o.append(t(x0 + w / 2, H - 4, "score", 9, TXT2, "middle"))
    return "".join(o), W, H


# --------------------------------------------------------------- panel b ---
def panel_b(W=234, H=250):
    o = []
    o.append(t(0, 11, "inferCNV score", 10, TXT, weight="600"))
    o.append(t(0, 24, "their Fig. 7G: tumour separates from normal B", 9, TXT2))
    x0, w = 26, W - 40
    yb = 150
    series = [("tonsil reference B (n = 5,298)", "ref", CN),
              ("normal B in the same sample (n = 819)", "nontumor", C3),
              ("tumour (n = 19,975)", "tumor", C4)]
    for lab, key, col in series:
        h = F["cnv"][key]
        o.append('<path d="%s" fill="%s" fill-opacity="0.22" stroke="%s" stroke-width="1.5"/>'
                 % (density_path(h, x0, w, yb, 100), col, col))
    o.append(axis(x0, yb, x0 + w, yb))
    h = F["cnv"]["ref"]
    xx = xpos(h, F["cnv"]["p99_ref"], x0, w)
    o.append(axis(xx, yb - 104, xx, yb, TXT2, 0.9, "2 2"))
    o.append(t(xx + 3, yb - 96, "99th pct of", 8.5, TXT2))
    o.append(t(xx + 3, yb - 87, "reference", 8.5, TXT2))
    for v, lb in ((0.008, "0.008"), (0.012, "0.012"), (0.016, "0.016")):
        xp = xpos(h, v, x0, w)
        o.append(axis(xp, yb, xp, yb + 3, LINE))
        o.append(t(xp, yb + 13, lb, 8.5, TXT2, "middle"))
    y = yb + 28
    for lab, key, col in series:
        o.append('<rect x="%.0f" y="%.0f" width="9" height="9" fill="%s" fill-opacity="0.55" '
                 'stroke="%s" stroke-width="1.2"/>' % (x0 - 12, y - 8, col, col))
        o.append(t(x0 + 1, y, lab, 8.6, TXT2))
        y += 11.4
    o.append(t(x0 - 12, y + 6, "100%% (M102) and %.0f%% (M413) of tumour cells"
               % (D["infercnv"]["M413"]["frac_tumor_above_ref_p99"] * 100), 8.6, TXT2))
    o.append(t(x0 - 12, y + 16, "sit above the reference 99th pct; none of the", 8.6, TXT2))
    o.append(t(x0 - 12, y + 26, "normal B cells beside them do.", 8.6, TXT2))
    return "".join(o), W, H


# --------------------------------------------------------------- panel c ---
def panel_c(W=268, H=250):
    o = []
    o.append(t(0, 11, "Metallothionein score per tumour cluster", 10, TXT, weight="600"))
    o.append(t(0, 24, "their Fig. 7I, shown there for M102 only", 9, TXT2))
    x0 = 34
    colw = (W - x0 - 10)
    lo, hi = -0.32, 0.68
    yb0 = 44

    def xs(v):
        return x0 + colw * (v - lo) / (hi - lo)

    o.append(axis(xs(0), yb0, xs(0), 196, LINE, 1, "2 2"))
    y = yb0 + 12
    for d, col in (("M102", C1), ("M413", C2)):
        o.append(t(x0 - 32, y + 3, d, 9.5, TXT, weight="600"))
        cl = F["mt"][d]["clusters"]
        order = sorted(cl, key=lambda k: -cl[k])
        for k in order:
            v = cl[k]
            hit = k in F["mt"][d]["elevated"]
            x1_, x2_ = (xs(min(0, v)), xs(max(0, v)))
            o.append('<rect x="%.1f" y="%.1f" width="%.1f" height="7" fill="%s" '
                     'fill-opacity="%s"/>' % (x1_, y, max(1.0, x2_ - x1_), col,
                                              "0.95" if hit else "0.42"))
            if hit:
                o.append(t(x2_ + 4, y + 6.5, "cluster %s · %s cells" % (k, f"{F['mt'][d]['n'][k]:,}"),
                           8.5, TXT, weight="600"))
            y += 9.4
        y += 9
    o.append(axis(x0, 196, x0 + colw, 196))
    for v in (-0.2, 0, 0.2, 0.4, 0.6):
        o.append(axis(xs(v), 196, xs(v), 199, LINE))
        o.append(t(xs(v), 209, "%g" % v, 8.5, TXT2, "middle"))
    o.append(t(x0 + colw / 2, H - 27, "mean score", 9, TXT2, "middle"))
    o.append(t(x0 - 32, H - 11, "One M102 cluster carries the programme; in M413", 8.8, TXT2))
    o.append(t(x0 - 32, H - 1, "every cluster is flat and negative.", 8.8, TXT2))
    return "".join(o), W, H


# --------------------------------------------------------------- panel d ---
def panel_d(W=300, H=250):
    o = []
    o.append(t(0, 11, "Which normal state does each tumour cell match?", 10, TXT, weight="600"))
    o.append(t(0, 24, "our test of their verbal conclusion", 9, TXT2))
    x0 = 118
    bw = W - x0 - 30
    iv = set(F["interval_states"])
    y = 44
    mx = max(s["frac"] for s in F["states"])
    for s in F["states"][:6]:
        inint = s["state"] in iv
        col = C2 if inint else C1
        o.append(t(x0 - 5, y + 7.5, s["state"][:26], 8.8,
                   TXT if s["frac"] > 0.1 else TXT2, "end"))
        o.append('<rect x="%.1f" y="%.1f" width="%.1f" height="8.5" fill="%s" '
                 'fill-opacity="%s"/>' % (x0, y, max(0.8, bw * s["frac"] / mx), col,
                                          "0.9" if s["frac"] > 0.1 else "0.6"))
        o.append(t(x0 + bw * s["frac"] / mx + 4, y + 7.5,
                   "%.1f%%" % (s["frac"] * 100), 8.5, TXT2))
        y += 11.6
    o.append(t(x0 - 5, y + 10, "orange = the interval they name", 8.5, C2, "end"))
    o.append(t(x0 - 5, y + 20, "(activated NBC → GC-committed)", 8.5, TXT2, "end"))

    yb = 196
    o.append(t(0, yb - 46, "How concentrated is that match?", 9.5, TXT, weight="600"))
    o.append(t(0, yb - 35, "effective number of normal states per cell (max 5)", 8.8, TXT2))
    xx0, ww = 10, W - 30
    for key, col, lab in (("nontumor", C3, "normal B"), ("tumor", C4, "tumour")):
        h = F["eff"][key]
        o.append('<path d="%s" fill="%s" fill-opacity="0.24" stroke="%s" stroke-width="1.5"/>'
                 % (density_path(h, xx0, ww, yb, 28), col, col))
    o.append(axis(xx0, yb, xx0 + ww, yb))
    h = F["eff"]["tumor"]
    for key, col, lab in (("tumor", C4, "tumour %.2f" % F["eff"]["tumor_median"]),
                          ("nontumor", C3, "normal B %.2f" % F["eff"]["nontumor_median"])):
        xv = xpos(h, F["eff"][key + "_median"], xx0, ww)
        o.append(axis(xv, yb - 30, xv, yb, col, 1.1, "2 2"))
    for v in (1, 2, 3, 4, 5):
        xv = xpos(h, v, xx0, ww)
        o.append(axis(xv, yb, xv, yb + 3, LINE))
        o.append(t(xv, yb + 13, "%d" % v, 8.5, TXT2, "middle"))
    o.append('<rect x="%.0f" y="%.0f" width="9" height="9" fill="%s" fill-opacity="0.55" '
             'stroke="%s" stroke-width="1.2"/>' % (xx0, yb + 22, C3, C3))
    o.append(t(xx0 + 13, yb + 30, "normal B in the same sample  %.2f"
               % F["eff"]["nontumor_median"], 8.8, TXT2))
    o.append('<rect x="%.0f" y="%.0f" width="9" height="9" fill="%s" fill-opacity="0.55" '
             'stroke="%s" stroke-width="1.2"/>' % (xx0, yb + 35, C4, C4))
    o.append(t(xx0 + 13, yb + 43, "tumour  %.2f   p = 1.3 × 10⁻²⁷"
               % F["eff"]["tumor_median"], 8.8, TXT2))
    return "".join(o), W, H


def main():
    gap, head = 14, 17
    items = [panel_a(), panel_b(), panel_c(), panel_d()]
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
    with open(os.path.join(FIGS, "revisit.svg"), "w") as fh:
        fh.write(svg)
    print("revisit.svg  %d x %d  %.1f kB" % (W, H, len(svg) / 1024))
    for mode in ("light", "dark"):
        flat = resolve(svg, mode)
        flat = flat.replace("<svg ", '<svg style="background:%s" ' % THEME[mode]["--surface-2"], 1)
        with open(os.path.join(FIGS, "revisit.%s.svg" % mode), "w") as fh:
            fh.write(flat)


if __name__ == "__main__":
    main()
