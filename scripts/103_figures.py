#!/usr/bin/env python3
"""Structure and locus figures for results/showcase.html.

Everything is drawn from real coordinates or real tables, as vector SVG, with
the page's own colour variables so the figures follow light and dark mode.

  7X5E                      the MafG-NRF2 heterodimer on its DNA site
  Boltz-2 predictions       the designed binder on the same MafG surface
  results/102_alascan.csv   binding energy cost per position (StaB-ddG)
  results/26_peak_annotated.csv + 59_axis_elements_atac.csv   the locus tracks
  data/genome/JASPAR2024_CORE_vertebrates_nr_pfms.txt         the motif logo

Run:  env/protdes/bin/python -I scripts/103_figures.py
Writes results/figs/*.svg and, if cairosvg and Pillow are importable,
results/figs/*.gif as well.
"""
import csv
import gzip
import json
import math
import os
import re
import sys

import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIGS = os.path.join(BASE, "results", "figs")
os.makedirs(FIGS, exist_ok=True)

PDB7X5E = os.path.join(BASE, "data", "structures", "7x5e.pdb")
PRED = "/var/tmp/hf-tonsil/repo/predictions/72_boltz_calibrated/boltz_results_inputs/predictions"
DESIGN = os.path.join(PRED, "mafg_dna_block_90_2_model_1_m4",
                      "mafg_dna_block_90_2_model_1_m4_model_0.cif")

AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
       "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
       "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V"}
NUC = {"DA", "DC", "DG", "DT", "A", "C", "G", "T", "U"}


# ---------------------------------------------------------------- parsing ---
def parse_pdb(path, atoms=("CA", "P", "C1'")):
    """-> list of dicts, in file order."""
    out = []
    for ln in open(path):
        if not ln.startswith("ATOM"):
            continue
        name = ln[12:16].strip()
        if name not in atoms:
            continue
        out.append(dict(chain=ln[21], resi=int(ln[22:26]), resn=ln[17:20].strip(),
                        atom=name, xyz=np.array([float(ln[30:38]), float(ln[38:46]),
                                                 float(ln[46:54])])))
    return out


def parse_cif(path, atoms=("CA", "P", "C1'")):
    out = []
    for ln in open(path):
        if not ln.startswith("ATOM"):
            continue
        f = ln.split()
        if f[3] not in atoms:
            continue
        out.append(dict(chain=f[15], resi=int(f[7]), resn=f[5], atom=f[3],
                        xyz=np.array([float(f[10]), float(f[11]), float(f[12])])))
    return out


def chain_trace(atoms, chain, atom="CA"):
    sel = [a for a in atoms if a["chain"] == chain and a["atom"] == atom]
    sel.sort(key=lambda a: a["resi"])
    return sel


def seq_of(trace):
    return "".join(AA3.get(a["resn"], "X") for a in trace)


# ------------------------------------------------------------- superpose ---
def kabsch(P, Q):
    """Rotation+translation that puts P onto Q."""
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, qc - R @ pc


def match_by_sequence(tr_a, tr_b):
    """Pair residues of two traces of the same protein by their sequence."""
    from difflib import SequenceMatcher
    sa, sb = seq_of(tr_a), seq_of(tr_b)
    pairs = []
    for blk in SequenceMatcher(None, sa, sb, autojunk=False).get_matching_blocks():
        for k in range(blk.size):
            pairs.append((tr_a[blk.a + k], tr_b[blk.b + k]))
    return pairs


# ---------------------------------------------------------------- drawing ---
def catmull(points, n=5):
    """Smooth a polyline through its points."""
    p = np.asarray(points, float)
    if len(p) < 3:
        return p
    ext = np.vstack([p[0], p, p[-1]])
    out = []
    for i in range(len(ext) - 3):
        p0, p1, p2, p3 = ext[i:i + 4]
        for t in np.linspace(0, 1, n, endpoint=False):
            t2, t3 = t * t, t * t * t
            out.append(0.5 * ((2 * p1) + (-p0 + p2) * t +
                              (2 * p0 - 5 * p1 + 4 * p2 - p3) * t2 +
                              (-p0 + 3 * p1 - 3 * p2 + p3) * t3))
    out.append(p[-1])
    return np.asarray(out)


def rotation(ax, ay):
    ca, sa = math.cos(ay), math.sin(ay)
    cb, sb = math.cos(ax), math.sin(ax)
    Ry = np.array([[ca, 0, sa], [0, 1, 0], [-sa, 0, ca]])
    Rx = np.array([[1, 0, 0], [0, cb, -sb], [0, sb, cb]])
    return Rx @ Ry


class Scene:
    """Painter's algorithm over short, continuous spans of a tube.

    A span is drawn as one path with a halo in the surface colour, so a span
    in front hides what is behind it. Spans overlap by one point, so the join
    inside a tube stays closed."""

    SPAN = 10

    def __init__(self, W, H, surface="var(--surface-2, #ffffff)"):
        self.W, self.H, self.surface = W, H, surface
        self.items = []          # (depth, halo_d, core_d, color, width, halo, op, kind)

    def path(self, pts, color, width=3.4, halo=2.3, opacity=None):
        pts = np.asarray(pts, float)
        if len(pts) < 2:
            return
        def dstr(a):
            return "M%d %d" % (round(a[0][0]), round(a[0][1])) + "".join(
                "L%d %d" % (round(q[0]), round(q[1])) for q in a[1:])
        for i in range(0, len(pts) - 1, self.SPAN):
            halo_seg = pts[i:i + self.SPAN + 1]
            core_seg = pts[i:i + self.SPAN + 2]      # one segment of overlap
            if len(halo_seg) < 2:
                continue
            self.items.append((float(halo_seg[:, 2].mean()), dstr(halo_seg),
                               dstr(core_seg), color, width, halo, opacity, "p"))

    def ball(self, p, color, r=3.6, halo=1.8):
        self.items.append((float(p[2]), None, None, color, r, halo, 1.0, "c", p))

    def seg(self, a, b, color, width=1.6, halo=1.4):
        d = "M%d %d L%d %d" % (round(a[0]), round(a[1]), round(b[0]), round(b[1]))
        self.items.append((float((a[2] + b[2]) / 2), d, d, color, width, halo, None, "p"))

    def svg(self, extra="", cls="", title=""):
        zs = [it[0] for it in self.items] or [0]
        lo, hi = min(zs), max(zs)
        rng = (hi - lo) or 1.0
        out = []
        for it in sorted(self.items, key=lambda i: i[0]):
            z, d, dcore, color, w, halo, op, kind = it[:8]
            o = op if op is not None else round(0.55 + 0.45 * (z - lo) / rng, 2)
            if kind == "p":
                if halo:
                    out.append('<path d="%s" stroke="%s" stroke-width="%.1f" '
                               'stroke-linecap="butt"/>' % (d, self.surface, w + 2 * halo))
                out.append('<path d="%s" stroke="%s" stroke-width="%.1f" stroke-opacity="%s"/>'
                           % (dcore, color, w, o))
            else:
                p = it[8]
                if halo:
                    out.append('<circle cx="%.1f" cy="%.1f" r="%.1f" fill="%s" stroke="none"/>'
                               % (p[0], p[1], w + halo, self.surface))
                out.append('<circle cx="%.1f" cy="%.1f" r="%.1f" fill="%s" stroke="none"/>'
                           % (p[0], p[1], w, color))
        return ('<svg class="%s" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" '
                'width="%d" height="%d" role="img" fill="none" stroke-linecap="round" '
                'stroke-linejoin="round">%s%s%s</svg>'
                % (cls, self.W, self.H, self.W, self.H,
                   ("<title>%s</title>" % title) if title else "",
                   "".join(out), extra))


def project(xyz, R, scale, center, W, H, flip=True):
    p = (xyz - center) @ R.T * scale
    out = p.copy()
    out[:, 0] += W / 2
    out[:, 1] = (H / 2 - p[:, 1]) if flip else (H / 2 + p[:, 1])
    return out


def fit(all_xyz, W, H, pad=18):
    c = (all_xyz.max(0) + all_xyz.min(0)) / 2
    span = (all_xyz.max(0) - all_xyz.min(0))[:2].max()
    return c, min(W - 2 * pad, H - 2 * pad) / max(span, 1e-6)


def view_matrix(dna_xyz, prot_xyz):
    """Screen x along the DNA axis, screen y towards the protein."""
    d = dna_xyz - dna_xyz.mean(0)
    _, _, vt = np.linalg.svd(d, full_matrices=False)
    x = vt[0] / np.linalg.norm(vt[0])
    v = prot_xyz.mean(0) - dna_xyz.mean(0)
    y = v - x * (v @ x)
    y /= np.linalg.norm(y)
    return np.vstack([x, y, np.cross(x, y)])


C1, C2, C3, C4 = ("var(--series-1, #0072b2)", "var(--series-2, #e69f00)",
                  "var(--series-3, #009e73)", "var(--series-4, #cc79a7)")
CN = "var(--neutral, #8c8c88)"
TXT = "var(--text-primary, #111111)"
TXT2 = "var(--text-secondary, #4a4a48)"


def load_crystal():
    at = parse_pdb(PDB7X5E)
    mafg = chain_trace(at, "A")
    nrf2 = chain_trace(at, "B")
    dna = {c: [a for a in at if a["chain"] == c and a["atom"] in ("P", "C1'")]
           for c in "CD"}
    return at, mafg, nrf2, dna


def dna_parts(dna):
    """Backbone polylines and the base-pair rungs, found by distance."""
    bb = {c: np.array([a["xyz"] for a in sorted([x for x in dna[c] if x["atom"] == "P"],
                                                key=lambda a: a["resi"])]) for c in dna}
    c1 = {c: sorted([x for x in dna[c] if x["atom"] == "C1'"], key=lambda a: a["resi"])
          for c in dna}
    # Watson-Crick pairing: find the strand offset that makes the two
    # sequences complementary, rather than trusting the nearest C1' atom
    comp = {"A": "T", "T": "A", "G": "C", "C": "G"}
    sa = [a["resn"][-1] for a in c1["C"]]
    sb = [b["resn"][-1] for b in c1["D"]]
    best_s, best_n = None, -1
    for s in range(0, len(sa) + len(sb)):
        n = sum(1 for i, x in enumerate(sa)
                if 0 <= s - i < len(sb) and comp.get(x) == sb[s - i])
        if n > best_n:
            best_s, best_n = s, n
    rungs = []
    for i, a in enumerate(c1["C"]):
        j = best_s - i
        if 0 <= j < len(c1["D"]):
            b = c1["D"][j]
            rungs.append((a, b, float(np.linalg.norm(a["xyz"] - b["xyz"]))))
    return bb, c1, rungs


def motif_positions(c1):
    """The TGACTCA bases of the crystallised site, on both strands."""
    seq = "".join(a["resn"][-1] for a in c1["C"])
    i = seq.find("TGACTCA")
    return (set(a["resi"] for a in c1["C"][i:i + 7]) if i >= 0 else set()), i, seq


def ddg_table():
    """position -> binding energy cost, from the 191-position alanine scan."""
    out = {}
    with open(os.path.join(BASE, "results", "102_alascan.csv")) as fh:
        for r in csv.DictReader(fh):
            out[(r["chain"], int(r["uniprot_pos"]))] = float(r["ddg"])
    return out


RAMP_LIGHT = ["var(--ddg-0, #b6bcbf)", "var(--ddg-1, #8fc0dd)",
              "var(--ddg-2, #2e86bd)", "var(--ddg-3, #0a4f78)"]
RAMP_FLAT = ["#b6bcbf", "#8fc0dd", "#2e86bd", "#0a4f78"]
RAMP_DARK = RAMP_LIGHT
BINS = [0.5, 1.0, 2.0]


def ddg_bin(v):
    return sum(v >= b for b in BINS)


class Complex:
    """The crystal, plus any predicted partner superposed on its MafG chain."""

    def __init__(self):
        at, self.mafg, self.nrf2, dna = load_crystal()
        self.bb, self.c1, self.rungs = dna_parts(dna)
        self.motif, _, self.seq = motif_positions(self.c1)
        prot = np.array([a["xyz"] for a in self.mafg + self.nrf2])
        dnax = np.vstack([self.bb["C"], self.bb["D"]])
        self.V = view_matrix(dnax, prot)
        allp = np.vstack([prot, dnax])
        v = allp @ self.V.T
        # centre on the centroid in x and z (rotation changes them) and on the
        # middle of the bounding box in y, so the panel is not bottom-heavy
        cv = np.array([v[:, 0].mean(), (v[:, 1].min() + v[:, 1].max()) / 2, v[:, 2].mean()])
        self.center = self.V.T @ cv
        self.binder = None

    def add_design(self, path=DESIGN):
        at = parse_cif(path)
        chains = sorted({a["chain"] for a in at})
        traces = {c: chain_trace(at, c) for c in chains}
        # the chain that matches MafG by sequence is the target; the other is the binder
        target_seq = seq_of(self.mafg)
        best = max(chains, key=lambda c: len(
            [p for p in match_by_sequence(traces[c], self.mafg)]))
        pairs = match_by_sequence(traces[best], self.mafg)
        P = np.array([a["xyz"] for a, _ in pairs])
        Q = np.array([b["xyz"] for _, b in pairs])
        R, t = kabsch(P, Q)
        other = [c for c in chains if c != best]
        self.binder = [a["xyz"] @ R.T + t
                       for a in sorted(traces[other[0]], key=lambda a: a["resi"])]
        self.binder = np.array(self.binder)
        self.fit_rmsd = float(np.sqrt((((P @ R.T + t) - Q) ** 2).sum(1).mean()))
        self.n_match = len(pairs)

    def radius(self):
        pts = [np.array([a["xyz"] for a in self.mafg + self.nrf2]),
               self.bb["C"], self.bb["D"]]
        if self.binder is not None:
            pts.append(self.binder)
        p = (np.vstack(pts) - self.center) @ self.V.T
        return math.hypot(np.abs(p[:, 0]).max(), np.abs(p[:, 2]).max()), np.abs(p[:, 1]).max()

    def scene(self, W, H, angle=0.0, parts=("mafg", "nrf2", "dna"), colour="chain",
              ddg=None, ramp=RAMP_LIGHT, labels=(), pad=16):
        rx, ry = self.radius()
        scale = min((W / 2 - pad) / max(rx, 1e-6), (H / 2 - pad) / max(ry, 1e-6))
        R = rotation(0.0, math.radians(angle)) @ self.V
        sc = Scene(W, H)

        def proj(xyz):
            return project(np.atleast_2d(np.asarray(xyz, float)), R, scale,
                           self.center, W, H)

        if "dna" in parts:
            for c in "CD":
                sc.path(proj(catmull(self.bb[c], 5)), CN, 3.8, 2.4)
            for a, b, _ in self.rungs:
                pa, pb = proj(a["xyz"])[0], proj(b["xyz"])[0]
                hot = a["resi"] in self.motif
                sc.seg(pa, pb, C3 if hot else CN, 3.0 if hot else 1.6)
        K = 4
        for key, trace, col in (("mafg", self.mafg, C1), ("nrf2", self.nrf2, C2)):
            if key not in parts:
                continue
            xyz = np.array([a["xyz"] for a in trace])
            pts = proj(catmull(xyz, K))
            if colour == "ddg" and ddg is not None:
                ch = "A" if key == "mafg" else "B"
                bins = [ddg_bin(ddg.get((ch, a["resi"]), 0.0)) for a in trace]
                per = [bins[min(len(bins) - 1, i // K)] for i in range(len(pts))]
                i = 0
                while i < len(pts) - 1:
                    j = i
                    while j < len(pts) - 1 and per[j + 1] == per[i]:
                        j += 1
                    sc.path(pts[i:j + 2], ramp[per[i]], 4.0, 2.3)
                    i = j + 1
            else:
                sc.path(pts, col, 3.4, 2.4)
        if self.binder is not None and "binder" in parts:
            sc.path(proj(catmull(self.binder, 4)), C3, 3.4, 2.4)
        txt = []
        for (ch, resi, lab, col) in labels:
            trace = self.mafg if ch == "A" else self.nrf2
            hit = [a for a in trace if a["resi"] == resi]
            if not hit:
                continue
            p = proj(hit[0]["xyz"])[0]
            sc.ball(p, col, 3.8, 1.8)
            side = -1 if ch == "A" else 1
            txt.append('<text x="%.1f" y="%.1f" text-anchor="%s" font-size="9.5" '
                       'fill="%s" font-weight="700">%s</text>'
                       % (p[0] + side * 7, p[1] + 3.2, "end" if side < 0 else "start",
                          TXT, lab))
        return sc, "".join(txt)


# ------------------------------------------------------------ composition ---
def panels(items, gap=10, letters=True, head=16):
    """Lay out several scenes side by side with panel letters a, b, c."""
    W = sum(w for _, w, _ in items) + gap * (len(items) - 1)
    H = max(h for _, _, h in items) + head
    out = []
    x = 0
    for i, (inner, w, h) in enumerate(items):
        lab = ('<text x="0" y="11" font-size="12.5" font-weight="700" fill="%s">%s</text>'
               % (TXT, chr(97 + i))) if letters else ""
        out.append('<g transform="translate(%d,0)">%s<g transform="translate(0,%d)">%s</g></g>'
                   % (x, lab, head, inner))
        x += w + gap
    return W, H, "".join(out)


def strip_svg(svg):
    """Inner content of an <svg> produced by Scene.svg()."""
    return re.sub(r"^<svg[^>]*>|</svg>$", "", svg)


def legend_row(x, y, entries, size=10.5):
    out = []
    for (col, txt) in entries:
        out.append('<rect x="%.0f" y="%.0f" width="9" height="9" fill="%s"/>' % (x, y - 8, col))
        out.append('<text x="%.0f" y="%.0f" font-size="%s" fill="%s">%s</text>'
                   % (x + 13, y, size, TXT2, txt))
        x += 17 + 6.1 * len(txt)
    return "".join(out), x


# -------------------------------------------------------- figure builders ---
def fig_structures(cx, ddg, W=258, H=420):
    hot = sorted([(v, k) for k, v in ddg.items() if v >= 2.0], reverse=True)
    labels = [("A", 79, "L79", RAMP_LIGHT[3]), ("B", 525, "L525", RAMP_LIGHT[3]),
              ("A", 93, "L93", RAMP_LIGHT[3]), ("B", 532, "L532", RAMP_LIGHT[3])]
    a, ta = cx.scene(W, H, 0)
    b, tb = cx.scene(W, H, 0, colour="ddg", ddg=ddg, labels=labels)
    c, tc = cx.scene(W, H, 0, parts=("mafg", "binder", "dna"))
    items = [(strip_svg(s.svg(extra=t)), W, H) for s, t in ((a, ta), (b, tb), (c, tc))]
    Wt, Ht, body = panels(items)
    l1, _ = legend_row(0, Ht + 14, [(C1, "MafG"), (C2, "NRF2"), (CN, "DNA"),
                                    (C3, "TGACTCA site · designed binder")])
    l2 = ""
    x = 0
    for i, lab in enumerate(["&lt; 0.5", "0.5–1.0", "1.0–2.0", "≥ 2.0"]):
        l2 += ('<rect x="%.0f" y="%.0f" width="9" height="9" fill="%s"/>'
               '<text x="%.0f" y="%.0f" font-size="10.5" fill="%s">%s</text>'
               % (W + 10 + x, Ht + 29 - 8, RAMP_LIGHT[i], W + 10 + x + 12, Ht + 29, TXT2, lab))
        x += 15 + 6.1 * len(lab)
    l2 += ('<text x="%.0f" y="%.0f" font-size="10.5" fill="%s">kcal mol⁻¹</text>'
           % (W + 10 + x + 2, Ht + 29, TXT2))
    return ('<svg class="mol" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" '
            'width="%d" height="%d" role="img" fill="none" stroke-linecap="round" '
            'stroke-linejoin="round">%s%s%s</svg>'
            % (Wt, Ht + 38, Wt, Ht + 38, body, l1, l2))


def read_pfm(name):
    rows, grab = {}, False
    for ln in open(os.path.join(BASE, "data", "genome",
                                "JASPAR2024_CORE_vertebrates_nr_pfms.txt")):
        if ln.startswith(">"):
            grab = ln.split()[0][1:] == name
            continue
        if grab:
            b = ln.split()[0]
            rows[b] = [float(v) for v in re.findall(r"[\d.]+", ln.split("[")[1])]
            if b == "T":
                break
    return rows


def logo_svg(pfm, W, H, title, pad_l=26, pad_b=16, pad_t=12):
    bases = "ACGT"
    col = {"A": C3, "C": C1, "G": C2, "T": C4}
    n = len(pfm["A"])
    bw = (W - pad_l - 4) / n
    ph = H - pad_b - pad_t
    out = ['<text x="0" y="9" font-size="10.5" font-weight="700" fill="%s">%s</text>'
           % (TXT, title)]
    out.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="1"/>'
               % (pad_l, pad_t, pad_l, pad_t + ph, TXT))
    for b, lab in ((0, "0"), (1, "1"), (2, "2")):
        y = pad_t + ph * (1 - b / 2)
        out.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="1"/>'
                   % (pad_l - 4, y, pad_l, y, TXT))
        out.append('<text x="%.1f" y="%.1f" font-size="9.5" fill="%s" text-anchor="end">%s</text>'
                   % (pad_l - 7, y + 3, TXT2, lab))
    out.append('<text transform="translate(9,%.1f) rotate(-90)" font-size="9.5" fill="%s" '
               'text-anchor="middle">bits</text>' % (pad_t + ph / 2, TXT2))
    for i in range(n):
        tot = sum(pfm[b][i] for b in bases) or 1
        p = {b: pfm[b][i] / tot for b in bases}
        ic = 2 + sum(p[b] * math.log2(p[b]) for b in bases if p[b] > 0)
        y = pad_t + ph
        for b in sorted(bases, key=lambda b: p[b]):
            h = ph * p[b] * ic / 2
            if h < 0.6:
                continue
            out.append('<g transform="translate(%.2f,%.2f) scale(%.4f,%.4f)">'
                       '<text x="0" y="0" font-size="20" font-weight="700" fill="%s" '
                       'font-family="Helvetica,Arial,sans-serif" text-anchor="middle">%s</text></g>'
                       % (pad_l + bw * (i + 0.5), y, bw / 13.4, h / 14.4, col[b], b))
            y -= h
    return "".join(out)


def fig_logos(cx, W=700):
    """Two information-content logos, and the DNA that the crystal holds."""
    maf, ap1 = read_pfm("MA0089.3"), read_pfm("MA0099.4")
    half, LH = (W - 26) / 2, 104
    out = ['<text x="0" y="10" font-size="12.5" font-weight="700" fill="%s">a</text>' % TXT,
           '<text x="%.0f" y="10" font-size="12.5" font-weight="700" fill="%s">b</text>'
           % (half + 26, TXT),
           '<g transform="translate(0,16)">%s</g>'
           % logo_svg(ap1, half, LH, "AP-1  (FOS::JUN, JASPAR MA0099.4)"),
           '<g transform="translate(%.0f,16)">%s</g>'
           % (half + 26, logo_svg(maf, half, LH, "small Maf  (MAFG::NFE2L1, JASPAR MA0089.3)"))]
    # the site in the crystal, both strands
    y0 = LH + 40
    out.append('<text x="0" y="%.0f" font-size="12.5" font-weight="700" fill="%s">c</text>'
               % (y0, TXT))
    out.append('<text x="16" y="%.0f" font-size="10.5" fill="%s">The same site in the '
               'crystal structure 7X5E, chains C and D</text>' % (y0, TXT2))
    pair = {a["resi"]: b["resn"][-1] for a, b, _ in cx.rungs}
    top = [(a["resi"], a["resn"][-1]) for a in cx.c1["C"]]
    bw = 15.5
    x0 = 40
    for k, (resi, ch) in enumerate(top):
        on = resi in cx.motif
        x = x0 + bw * k
        for row, base in ((0, ch), (1, pair.get(resi, " "))):
            out.append('<text x="%.1f" y="%.1f" font-size="12" font-family="ui-monospace,'
                       'Menlo,monospace" font-weight="%s" fill="%s" text-anchor="middle">%s</text>'
                       % (x, y0 + 22 + row * 16, "700" if on else "400",
                          C3 if on else TXT2, base))
        if resi in pair:
            out.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
                       'stroke-width="1" stroke-opacity="%s"/>'
                       % (x, y0 + 25, x, y0 + 31, C3 if on else CN, "1" if on else "0.6"))
    for row, lab in ((0, "5'"), (1, "3'")):
        out.append('<text x="%.0f" y="%.1f" font-size="10" fill="%s" text-anchor="end">%s</text>'
                   % (x0 - 12, y0 + 22 + row * 16, TXT2, lab))
    xe = x0 + bw * len(top)
    out.append('<text x="%.0f" y="%.1f" font-size="10" fill="%s">%s</text>'
               % (xe + 6, y0 + 22, TXT2, "3'"))
    out.append('<text x="%.0f" y="%.1f" font-size="10" fill="%s">%s</text>'
               % (xe + 6, y0 + 38, TXT2, "5'"))
    out.append('<text x="%.0f" y="%.1f" font-size="10.5" font-weight="700" fill="%s">'
               'TGACTCA</text>' % (xe + 34, y0 + 30, C3))
    H = y0 + 48
    return ('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" width="%d" '
            'height="%d" role="img">%s</svg>' % (W, H, W, H, "".join(out)))


# ------------------------------------------------------------ locus view ---
# SLC7A11, GENCODE v44: chr4:138,164,097-138,242,349, minus strand.
GENE = dict(name="SLC7A11", chrom="chr4", start=138164097, end=138242349, strand="-")
TSS = GENE["end"]


def locus_data(win=(138150000, 138420000)):
    lo, hi = win
    peaks = []
    with open(os.path.join(BASE, "results", "26_peak_annotated.csv")) as fh:
        for r in csv.DictReader(fh):
            if r["chrom"] != GENE["chrom"]:
                continue
            s, e = int(r["start"]), int(r["end"])
            if e < lo or s > hi:
                continue
            peaks.append(dict(start=s, end=e, mid=(s + e) / 2,
                              mcl=float(r["mcl_mean_log2cpm"]),
                              normal=float(r["normal_log2cpm"]), peak=r["peak"]))
    el = []
    with open(os.path.join(BASE, "results", "59_axis_elements_atac.csv")) as fh:
        for r in csv.DictReader(fh):
            if r["gene"] == "SLC7A11":
                el.append(dict(peak=r["peak"], start=int(r["start"]), end=int(r["end"]),
                               dist=int(r["dist"]), mcl=float(r["mcl_mean_log2cpm"]),
                               normal=float(r["normal_log2cpm"])))
    pct = {}
    vals = {}
    with open(os.path.join(BASE, "results", "98_contact_multi.csv")) as fh:
        for r in csv.DictReader(fh):
            if r["gene"] == "SLC7A11" and r["kind"] == "de_novo" and r["group"].startswith("LCL"):
                vals.setdefault(int(r["dist"]), []).append(float(r["pct"]))
    for d, v in vals.items():
        v.sort()
        pct[d] = v[len(v) // 2] if len(v) % 2 else (v[len(v) // 2 - 1] + v[len(v) // 2]) / 2
    return peaks, el, pct, win


def fig_locus(W=720):
    peaks, el, pct, (lo, hi) = locus_data()
    eset = {e["peak"] for e in el}
    pad_l, pad_r = 104, 16
    pw = W - pad_l - pad_r
    def X(bp):
        return pad_l + (bp - lo) / (hi - lo) * pw
    arc_h, y_gene = 58, 92
    th = 44                                   # track height
    y_n, y_m = y_gene + 42, y_gene + 42 + th + 34
    yax = y_m + th + 15
    H = yax + 30
    mx = max([p["mcl"] for p in peaks] + [p["normal"] for p in peaks] + [1])
    o = []
    # contact arcs
    for i, e in enumerate(sorted(el, key=lambda e: e["dist"])):
        p = pct.get(e["dist"], 0)
        x0, x1 = X((e["start"] + e["end"]) / 2), X(TSS)
        h = arc_h * (0.25 + 0.75 * p)
        o.append('<path d="M%.1f %.1f Q %.1f %.1f %.1f %.1f" stroke="%s" stroke-width="%.1f" '
                 'fill="none" stroke-opacity="%.2f"/>'
                 % (x0, y_gene - 6, (x0 + x1) / 2, y_gene - 6 - h * 2, x1, y_gene - 6,
                    C1, 1.0 + 1.6 * p, 0.35 + 0.65 * p))
        o.append('<text x="%.1f" y="%.1f" font-size="9.5" fill="%s" text-anchor="middle" '
                 'font-weight="700">%.2f</text>'
                 % (x0, y_gene - 14 - h * 0.62 - 11 * (i % 2), C1, p))
    o.append('<text x="0" y="%.0f" font-size="10.5" font-weight="700" fill="%s">Contact rank'
             '</text>' % (y_gene - 44, TXT))
    o.append('<text x="0" y="%.0f" font-size="9.5" fill="%s">distance-matched,</text>'
             % (y_gene - 31, TXT2))
    o.append('<text x="0" y="%.0f" font-size="9.5" fill="%s">6 B-cell lines</text>'
             % (y_gene - 19, TXT2))
    # gene model
    gx0, gx1 = X(GENE["start"]), X(GENE["end"])
    o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="4"/>'
             % (gx0, y_gene, gx1, y_gene, TXT))
    for k in range(5):
        xa = gx1 - (gx1 - gx0) * (k + 0.5) / 5
        o.append('<path d="M%.1f %.1f l5 -3.4 l0 6.8 z" fill="%s"/>' % (xa, y_gene, TXT))
    o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="1.4"/>'
             % (gx1, y_gene - 14, gx1, y_gene + 6, TXT))
    o.append('<text x="%.1f" y="%.1f" font-size="10.5" font-weight="700" fill="%s" '
             'text-anchor="end">%s</text>' % (gx1 - 6, y_gene + 17, TXT, GENE["name"]))
    o.append('<text x="%.1f" y="%.1f" font-size="9.5" fill="%s">TSS</text>'
             % (gx1 + 4, y_gene + 17, TXT2))
    # accessibility tracks
    for (y0, key, lab) in ((y_n, "normal", "Normal B cells"), (y_m, "mcl", "MCL cells")):
        base = y0 + th
        o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="1"/>'
                 % (pad_l, base, W - pad_r, base, "var(--line, #d6d6d2)"))
        o.append('<text x="0" y="%.0f" font-size="10.5" font-weight="700" fill="%s">%s</text>'
                 % (y0 + 11, TXT, lab))
        o.append('<text x="0" y="%.0f" font-size="9.5" fill="%s">open chromatin</text>'
                 % (y0 + 24, TXT2))
        o.append('<text x="0" y="%.0f" font-size="9.5" fill="%s">log₂ CPM</text>'
                 % (y0 + 36, TXT2))
        for p in peaks:
            h = th * max(p[key], 0) / mx
            x0, x1 = X(p["start"]), X(p["end"])
            w = max(2.4, x1 - x0)
            hot = p["peak"] in eset
            o.append('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="%s"/>'
                     % (x0, base - h, w, h, C1 if hot else CN))
        o.append('<text x="%.1f" y="%.1f" font-size="9" fill="%s" text-anchor="end">%.0f</text>'
                 % (pad_l - 5, y0 + 8, TXT2, mx))
        o.append('<text x="%.1f" y="%.1f" font-size="9" fill="%s" text-anchor="end">0</text>'
                 % (pad_l - 5, base + 3, TXT2))
    # element flags
    for i, e in enumerate(sorted(el, key=lambda e: e["dist"])):
        x = X((e["start"] + e["end"]) / 2)
        o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="1" '
                 'stroke-dasharray="2 2"/>' % (x, y_gene + 6, x, y_n - 14, C1))
        o.append('<text x="%.1f" y="%.1f" font-size="9" fill="%s" text-anchor="middle" '
                 'font-weight="700">%s</text>'
                 % (x, y_n - 16 - 11 * (i % 2), C1, ("%+d kb" % round(e["dist"] / 1000))))
    # axis
    o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="1"/>'
             % (pad_l, yax, W - pad_r, yax, TXT))
    for kb in range(-100, 200, 50):
        bp = TSS + kb * 1000
        if not (lo <= bp <= hi):
            continue
        o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="1"/>'
                 % (X(bp), yax, X(bp), yax + 4, TXT))
        o.append('<text x="%.1f" y="%.1f" font-size="9.5" fill="%s" text-anchor="middle">%s</text>'
                 % (X(bp), yax + 15, TXT2, ("%+d" % kb) if kb else "0"))
    o.append('<text x="%.1f" y="%.1f" font-size="10" fill="%s" text-anchor="middle">'
             'Distance from the SLC7A11 promoter (kb)</text>' % (pad_l + pw / 2, H - 2, TXT))
    leg, _ = legend_row(0, H - 2, [(C1, "the five new elements"), (CN, "other peaks")])
    o.append(leg)
    return ('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" width="%d" height="%d" '
            'role="img">%s</svg>' % (W, H, W, H, "".join(o)))


# ------------------------------------------------------------- animation ---
THEME = {
    "light": {"--series-1": "#0072b2", "--series-2": "#e69f00", "--series-3": "#009e73",
              "--series-4": "#cc79a7", "--neutral": "#8c8c88", "--surface-2": "#ffffff",
              "--text-primary": "#111111", "--text-secondary": "#4a4a48",
              "--line": "#d6d6d2", "--ddg-0": "#b6bcbf", "--ddg-1": "#8fc0dd",
              "--ddg-2": "#2e86bd", "--ddg-3": "#0a4f78"},
    "dark": {"--series-1": "#2787c9", "--series-2": "#c08400", "--series-3": "#099c72",
             "--series-4": "#c3729f", "--neutral": "#8a8a84", "--surface-2": "#171716",
             "--text-primary": "#ededea", "--text-secondary": "#bcbcb6",
             "--line": "#343431", "--ddg-0": "#5d666b", "--ddg-1": "#2b7aa8",
             "--ddg-2": "#4ea7d8", "--ddg-3": "#9fd4f0"},
}


def resolve(svg, mode="light"):
    th = THEME[mode]
    return re.sub(r"var\((--[\w-]+), (#\w+)\)",
                  lambda m: th.get(m.group(1), m.group(2)), svg)


OUTLINE = "#4d4d49"      # halo colour for the transparent animation: a soft
                         # outline that stays visible on a white and on a dark page


def make_gif(cx, path, n=24, W=250, H=404, scale=1.35, ms=100, mode="light",
             swap=False, **kw):
    """One turn of the complex as a GIF.

    The background is transparent and the halo is a grey outline, so the same
    file sits on the light and on the dark page."""
    import cairosvg
    from PIL import Image
    import io
    frames = []
    for k in range(n):
        half = k < n / 2
        if swap:
            kw = dict(kw, parts=("mafg", "nrf2", "dna") if half
                      else ("mafg", "binder", "dna"))
        sc, txt = cx.scene(W, H, 360.0 * k / n, **kw)
        sc.surface = OUTLINE
        if swap:
            txt += ('<text x="%d" y="%d" font-size="10.5" fill="#8a8a86" '
                    'text-anchor="middle" font-weight="700">%s</text>'
                    % (W / 2, H - 4,
                       "natural partner NRF2 &#183; crystal" if half
                       else "designed binder &#183; model"))
        svg = resolve(sc.svg(extra=txt), mode)
        png = cairosvg.svg2png(bytestring=svg.encode(), scale=scale, background_color=None)
        im = Image.open(io.BytesIO(png)).convert("RGBA")
        bg = Image.new("RGBA", im.size, (255, 255, 255, 0))
        im = Image.alpha_composite(bg, im)
        a = im.getchannel("A").point(lambda v: 255 if v > 110 else 0)
        q = im.convert("RGB").quantize(colors=63, method=Image.MEDIANCUT)
        q.paste(63, (0, 0), a.point(lambda v: 255 - v))     # index 63 = transparent
        pal = q.getpalette()
        pal[63 * 3:63 * 3 + 3] = [255, 0, 255]
        q.putpalette(pal)
        q.info["transparency"] = 63
        frames.append(q)
    frames[0].save(path, save_all=True, append_images=frames[1:], duration=ms,
                   loop=0, optimize=True, disposal=2, transparency=63)
    return os.path.getsize(path)


# -------------------------------------------------------- study overview ---
def _t(x, y, s, size=10, fill=None, anchor="start", weight="400"):
    return ('<text x="%.1f" y="%.1f" font-size="%s" fill="%s" text-anchor="%s" '
            'font-weight="%s">%s</text>' % (x, y, size, fill or TXT2, anchor, weight, s))


def _arrow(x1, y1, x2, y2, col=None, w=1.4):
    col = col or TXT
    ang = math.atan2(y2 - y1, x2 - x1)
    hx, hy = x2 - 6 * math.cos(ang), y2 - 6 * math.sin(ang)
    px, py = 3.2 * math.sin(ang), -3.2 * math.cos(ang)
    return ('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="%s"/>'
            '<path d="M%.1f %.1f L%.1f %.1f L%.1f %.1f z" fill="%s"/>'
            % (x1, y1, hx, hy, col, w, x2, y2, hx + px, hy + py, hx - px, hy - py, col))


def _wave(x, y, w, amp, n, phase=0.0, step=3.0):
    pts = []
    k = 0.0
    while k <= w:
        pts.append("%s%.1f %.1f" % ("M" if not pts else "L", x + k,
                                    y + amp * math.sin(2 * math.pi * n * k / w + phase)))
        k += step
    return "".join(pts)


def _dna_icon(x, y, w, col=None, amp=7, n=2.4):
    col = col or CN
    o = ['<path d="%s" stroke="%s" stroke-width="2.4" fill="none"/>' % (_wave(x, y, w, amp, n), col),
         '<path d="%s" stroke="%s" stroke-width="2.4" fill="none"/>'
         % (_wave(x, y, w, amp, n, math.pi), col)]
    for k in range(1, 13):
        xx = x + w * k / 13.0
        y1 = y + amp * math.sin(2 * math.pi * n * (w * k / 13.0) / w)
        y2 = y + amp * math.sin(2 * math.pi * n * (w * k / 13.0) / w + math.pi)
        o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="1" '
                 'stroke-opacity=".7"/>' % (xx, y1, xx, y2, col))
    return "".join(o)


def _cell(x, y, r, col=None):
    col = col or TXT2
    return ('<circle cx="%.1f" cy="%.1f" r="%.1f" fill="none" stroke="%s" stroke-width="1.4"/>'
            '<circle cx="%.1f" cy="%.1f" r="%.1f" fill="%s" fill-opacity=".14" stroke="%s" '
            'stroke-width="1"/>' % (x, y, r, col, x, y, r * 0.52, col, col))


def fig_overview(cx, W=1100):
    H = 300
    o = []
    cols = [0, 272, 560, 846]
    heads = ["Measure", "Match", "Find", "Attack"]
    subs = ["one nucleus, two readouts", "against the same stage of development",
            "what development cannot explain", "the protein that runs the circuit"]
    for i, x in enumerate(cols):
        o.append(_t(x, 12, chr(97 + i), 12.5, TXT, weight="700"))
        o.append(_t(x + 16, 12, heads[i], 12, TXT, weight="700"))
        o.append(_t(x + 16, 25, subs[i], 9.5, TXT2))
        if i:
            o.append('<line x1="%.1f" y1="40" x2="%.1f" y2="%d" stroke="%s" stroke-width="1" '
                     'stroke-opacity=".45"/>' % (x - 16, x - 16, H - 10, "var(--line, #d6d6d2)"))

    # ---- a. multiome
    x0 = cols[0]
    o.append('<ellipse cx="%.1f" cy="86" rx="30" ry="24" fill="%s" fill-opacity=".13" '
             'stroke="%s" stroke-width="1.4"/>' % (x0 + 32, C4, C4))
    o.append('<circle cx="%.1f" cy="86" r="11" fill="%s" fill-opacity=".30" stroke="%s" '
             'stroke-width="1"/>' % (x0 + 32, C4, C4))
    o.append(_t(x0 + 32, 122, "tonsil", 9.5, TXT2, "middle"))
    o.append(_t(x0 + 32, 134, "46,216 B cells", 9.5, TXT2, "middle"))
    o.append(_arrow(x0 + 68, 86, x0 + 92, 86))
    o.append(_cell(x0 + 112, 86, 18, C1))
    o.append(_t(x0 + 112, 122, "one nucleus", 9.5, TXT2, "middle"))
    o.append(_arrow(x0 + 134, 78, x0 + 158, 62))
    o.append(_arrow(x0 + 134, 94, x0 + 158, 112))
    o.append('<path d="%s" stroke="%s" stroke-width="2.2" fill="none"/>'
             % (_wave(x0 + 162, 58, 64, 5, 2.2), C2))
    o.append(_t(x0 + 162, 42, "RNA", 10, TXT, weight="700"))
    o.append(_t(x0 + 196, 42, "22,345 genes", 9.5, TXT2))
    o.append(_dna_icon(x0 + 162, 112, 64, CN))
    for dx in (0, 46):
        o.append('<circle cx="%.1f" cy="112" r="7.5" fill="%s" fill-opacity=".85"/>'
                 % (x0 + 166 + dx, CN))
    o.append(_t(x0 + 162, 140, "open chromatin", 10, TXT, weight="700"))
    o.append(_t(x0 + 162, 152, "215,204 peaks", 9.5, TXT2))
    o.append(_t(x0, 182, "The same cell gives both readouts, so the two layers", 9.5, TXT2))
    o.append(_t(x0, 194, "can be compared without matching samples.", 9.5, TXT2))

    # ---- b. differentiation-matched comparison
    x0 = cols[1]
    arc = "M%.1f 150 Q %.1f 54 %.1f 150" % (x0 + 8, x0 + 112, x0 + 216)
    o.append('<path d="%s" stroke="%s" stroke-width="1.6" fill="none" stroke-opacity=".5"/>'
             % (arc, CN))
    for k in range(47):
        tt = k / 46.0
        px = (1 - tt) ** 2 * (x0 + 8) + 2 * (1 - tt) * tt * (x0 + 112) + tt ** 2 * (x0 + 216)
        py = (1 - tt) ** 2 * 150 + 2 * (1 - tt) * tt * 54 + tt ** 2 * 150
        o.append('<circle cx="%.1f" cy="%.1f" r="2.6" fill="%s" fill-opacity=".75"/>'
                 % (px, py, CN))
    o.append(_t(x0 + 8, 166, "naive", 9.5, TXT2))
    o.append(_t(x0 + 112, 166, "germinal centre", 9.5, TXT2, "middle"))
    o.append(_t(x0 + 216, 166, "memory", 9.5, TXT2, "end"))
    o.append(_t(x0, 182, "47 normal states", 9.5, TXT2))
    # the tumour cell and the normal state it is matched to
    tt = 0.33
    mx_ = (1 - tt) ** 2 * (x0 + 8) + 2 * (1 - tt) * tt * (x0 + 112) + tt ** 2 * (x0 + 216)
    my_ = (1 - tt) ** 2 * 150 + 2 * (1 - tt) * tt * 54 + tt ** 2 * 150
    tx, ty = mx_, 62
    o.append('<circle cx="%.1f" cy="%.1f" r="7" fill="%s"/>' % (tx, ty, C4))
    o.append(_t(tx + 12, ty + 3, "tumour cell", 9.5, C4, weight="700"))
    o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="1" '
             'stroke-dasharray="3 3"/>' % (tx, ty + 9, mx_, my_ - 7, C4))
    o.append('<circle cx="%.1f" cy="%.1f" r="4.6" fill="none" stroke="%s" stroke-width="1.8"/>'
             % (mx_, my_, C4))
    o.append(_t(x0 + 112, 136, "matched normal state", 9.5, C4, "middle"))
    # expected against measured
    bx = x0 + 8
    o.append(_t(bx, 212, "SLC7A11 at that state", 9.5, TXT, weight="700"))
    for k, (lab, v, col) in enumerate((("expected", 0.001, CN), ("measured", 0.206, C4))):
        yy = 226 + k * 18
        o.append(_t(bx, yy + 8, lab, 9.5, TXT2))
        o.append('<rect x="%.1f" y="%.1f" width="%.1f" height="9" fill="%s"/>'
                 % (bx + 52, yy, max(2.0, 150 * v / 0.206), col))
        o.append(_t(bx + 52 + max(2.0, 150 * v / 0.206) + 5, yy + 8, "%.3f" % v, 9.5, TXT2))
    o.append(_t(bx, 272, "the gap is the part development cannot explain", 9.5, TXT2))

    # ---- c. the circuit
    x0 = cols[2]
    o.append(_t(x0, 52, "strongest remaining signal", 9.5, TXT2))
    o.append('<rect x="%.1f" y="58" width="128" height="26" fill="%s" fill-opacity=".14" '
             'stroke="%s" stroke-width="1"/>' % (x0, C3, C3))
    o.append(_t(x0 + 64, 76, "TGACTCA", 14, C3, "middle", "700"))
    o.append(_t(x0 + 136, 70, "AP-1 and", 9.5, TXT2))
    o.append(_t(x0 + 136, 82, "small Maf", 9.5, TXT2))
    o.append(_arrow(x0 + 64, 90, x0 + 64, 106))
    o.append(_t(x0 + 72, 102, "genome-wide, 1.28 to 1.42&#215;", 9.5, TXT2))
    # locus sketch
    ly = 132
    o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="2.6"/>'
             % (x0, ly, x0 + 70, ly, TXT))
    o.append(_t(x0 + 2, ly - 7, "SLC7A11", 9.5, TXT, weight="700"))
    for k, dx in enumerate((92, 118, 150, 196, 214)):
        o.append('<rect x="%.1f" y="%.1f" width="5" height="12" fill="%s"/>'
                 % (x0 + dx, ly - 6, C1))
    o.append('<path d="M%.1f %.1f Q %.1f %.1f %.1f %.1f" stroke="%s" stroke-width="1.2" '
             'fill="none" stroke-opacity=".8"/>' % (x0 + 216, ly - 8, x0 + 150, ly - 46,
                                                    x0 + 70, ly - 8, C1))
    o.append(_t(x0 + 92, ly + 24, "five new elements", 9.5, C1, weight="700"))
    o.append(_t(x0 + 92, ly + 36, "up to 164 kb away", 9.5, TXT2))
    # pathway chain
    py0 = 212
    chain = [("SLC7A11", "cystine in"), ("GCLC/GCLM", "glutathione"), ("GSS", "glutathione"),
             ("AIFM2", "FSP1 defence")]
    bw2 = 58
    for k, (g, sub) in enumerate(chain):
        bx = x0 + k * (bw2 + 14)
        o.append('<rect x="%.1f" y="%.1f" width="%d" height="22" fill="%s" fill-opacity=".14" '
                 'stroke="%s" stroke-width="1"/>' % (bx, py0, bw2, C3, C3))
        o.append(_t(bx + bw2 / 2, py0 + 14, g, 9, TXT, "middle", "700"))
        o.append(_t(bx + bw2 / 2, py0 + 34, sub, 8.5, TXT2, "middle"))
        if k < len(chain) - 1:
            o.append(_arrow(bx + bw2 + 2, py0 + 11, bx + bw2 + 11, py0 + 11))
    o.append(_t(x0, py0 + 54, "five steps in a row are MCL-selective in DepMap", 9.5, TXT2))

    # ---- d. the drug
    x0 = cols[3]
    sm, tsm = cx.scene(118, 168, 0, pad=6)
    sb, tsb = cx.scene(118, 168, 0, parts=("mafg", "binder", "dna"), pad=6)
    o.append('<g transform="translate(%.1f,46)">%s</g>' % (x0, strip_svg(sm.svg(extra=tsm))))
    o.append('<g transform="translate(%.1f,46)">%s</g>' % (x0 + 134, strip_svg(sb.svg(extra=tsb))))
    o.append(_arrow(x0 + 120, 130, x0 + 132, 130))
    o.append(_t(x0, 226, "MafG with its natural", 9.5, TXT2))
    o.append(_t(x0, 238, "partner NRF2 on DNA", 9.5, TXT2))
    o.append(_t(x0 + 134, 226, "MafG with the designed", 9.5, TXT2))
    o.append(_t(x0 + 134, 238, "90-residue binder", 9.5, TXT2))
    o.append(_t(x0, 262, "25 hotspot positions mark the epitope.", 9.5, TXT2))
    o.append(_t(x0, 274, "Add a degradation tag and the cell removes MafG.", 9.5, TXT2))
    return ('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" width="%d" height="%d" '
            'role="img">%s</svg>' % (W, H, W, H, "".join(o)))


# ------------------------------------------------------------------ main ---
def main():
    cx = Complex()
    cx.add_design()
    ddg = ddg_table()
    out = {
        "overview.svg": fig_overview(cx),
        "motif.svg": fig_logos(cx),
        "locus_slc7a11.svg": fig_locus(),
        "structures.svg": fig_structures(cx, ddg),
    }
    for name, svg in out.items():
        with open(os.path.join(FIGS, name), "w") as fh:
            fh.write(svg)
        print("%-22s %6.1f kB" % (name, len(svg) / 1024))
        # Markdown cannot resolve the page's colour variables, so write one
        # file per mode for README.md and the report to point at.
        for mode in ("light", "dark"):
            flat = resolve(svg, mode)
            flat = flat.replace("<svg ", '<svg style="background:%s" '
                                % THEME[mode]["--surface-2"], 1)
            with open(os.path.join(FIGS, name.replace(".svg", ".%s.svg" % mode)), "w") as fh:
                fh.write(flat)
    print("MafG superposition: %d residues, RMSD %.2f A" % (cx.n_match, cx.fit_rmsd))
    try:
        import cairosvg  # noqa: F401
        from PIL import Image  # noqa: F401
    except ImportError:
        print("cairosvg or Pillow missing: the GIF was not written")
        return
    kb = make_gif(cx, os.path.join(FIGS, "mafg_rotation.gif"), swap=True)
    print("%-22s %6.1f kB" % ("mafg_rotation.gif", kb / 1024))


if __name__ == "__main__":
    main()
