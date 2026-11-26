#!/usr/bin/env python3
"""Estimate the printed height of each page of the report, per language.

There is no HTML renderer in this environment, so page overflow cannot be
checked by printing. This walks the template's block elements, applies the
CSS values actually used, and estimates height from text metrics:

  * average character advance 0.50 em for the serif body, 0.52 em for the
    sans used in tables and captions, 1.0 em for CJK
  * a block's height = ceil(chars / chars_per_line) * font_size * line_height
    + margins + padding + borders
  * a table row = max over its cells of the wrapped cell height + cell padding
  * an inline SVG = width * viewBox_height / viewBox_width

It is an estimate, not a renderer: treat anything within ~5% of the limit as
needing a human look. Run after editing docs/report_template.html.

Run:  env/tonsilmcl/bin/python -I scripts/113_estimate_pages.py
"""
import math
import os
import re
import sys

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TPL = os.path.join(BASE, "docs", "report_template.html")

MM_PER_PT = 25.4 / 72.0
PAGE_H_MM = 297.0
MARGIN_TOP_MM, MARGIN_BOT_MM, MARGIN_SIDE_MM = 13.0, 12.0, 14.0  # overwritten from @page
AVAIL_MM = PAGE_H_MM - MARGIN_TOP_MM - MARGIN_BOT_MM
WIDTH_MM = 210.0 - 2 * MARGIN_SIDE_MM

CJK = re.compile(r"[　-鿿＀-￯]")


def css_num(block, prop, default):
    """Read a px/pt value for a property out of the stylesheet text."""
    m = re.search(prop + r"\s*:\s*([\d.]+)(pt|px|mm)", block)
    if not m:
        return default
    v, u = float(m.group(1)), m.group(2)
    return v if u == "pt" else (v * 0.75 if u == "px" else v / MM_PER_PT)


def text_of(h):
    h = re.sub(r"<[^>]+>", " ", h)
    h = h.replace("&nbsp;", " ").replace("&amp;", "&").replace("&lt;", "<").replace("&gt;", ">")
    return re.sub(r"\s+", " ", h).strip()


def strip_other_lang(h, lang):
    other = "zh-Hant" if lang == "en" else "en"
    prev = None
    while prev != h:
        prev = h
        h = re.sub(r'<(\w+)\b[^>]*(?<![-\w])lang="%s"[^>]*>((?:(?!</?\1[\s>]).)*?)</\1>' % other,
                   "", h, flags=re.S)
    return h


def chars_per_line(text, width_mm, size_pt, sans=False):
    if not text:
        return 1
    cjk = len(CJK.findall(text))
    frac_cjk = cjk / max(1, len(text))
    adv_em = frac_cjk * 1.0 + (1 - frac_cjk) * (0.52 if sans else 0.50)
    adv_mm = adv_em * size_pt * MM_PER_PT
    return max(1.0, width_mm / adv_mm)


def wrapped_mm(text, width_mm, size_pt, lh, sans=False):
    if not text:
        return 0.0
    cpl = chars_per_line(text, width_mm, size_pt, sans)
    lines = max(1, math.ceil(len(text) / cpl))
    return lines * size_pt * lh * MM_PER_PT


def svg_mm(tag, width_mm):
    m = re.search(r'viewBox="0 0 ([\d.]+) ([\d.]+)"', tag)
    if not m:
        return 40.0
    w, h = float(m.group(1)), float(m.group(2))
    return width_mm * h / w


# ----------------------------------------------------------------- walk ---
def table_mm(tbl, width_mm, lang, base_pt, lh):  # noqa: ARG001
    """Tables use the sans stack at 8.5pt by default (see stylesheet)."""
    size = TABLE_PT
    cap = re.search(r"<caption>(.*?)</caption>", tbl, flags=re.S)
    h = 0.0
    if cap:
        h += wrapped_mm(text_of(strip_other_lang(cap.group(1), lang)),
                        width_mm, CAPTION_PT, 1.38, True) + 2.6 * MM_PER_PT
    rows = re.findall(r"<tr>(.*?)</tr>", tbl, flags=re.S)
    tight = "tight" in (re.search(r'<table[^>]*>', tbl).group(0) or "")
    pad = (2.3 if tight else 3.2) * 2 * MM_PER_PT
    for r in rows:
        cells = re.findall(r"<t[hd][^>]*>(.*?)</t[hd]>", r, flags=re.S)
        if not cells:
            continue
        widths = []
        for c in re.findall(r"<t[hd]([^>]*)>", r):
            m = re.search(r"width:\s*([\d.]+)%", c)
            widths.append(float(m.group(1)) / 100 if m else None)
        n = len(cells)
        fill = [w if w else (1.0 - sum(x for x in widths if x)) / max(1, sum(1 for x in widths if not x))
                for w in widths] or [1.0 / n] * n
        hr = 0.0
        for c, frac in zip(cells, fill + [1.0 / n] * n):
            cw = max(12.0, width_mm * frac) - 5 * MM_PER_PT
            hr = max(hr, wrapped_mm(text_of(strip_other_lang(c, lang)), cw, size, 1.38, True))
        h += hr + pad
    return h + 4 * MM_PER_PT + 7 * MM_PER_PT


BLOCK = re.compile(
    r"<(h1|h2|h3|p|ul|ol|table|figure|figcaption|div|section|svg)\b([^>]*)>", re.S)


def page_mm(sheet, lang, verbose=False):
    base = ZH_PT if lang == "zh" else BASE_PT
    total = 0.0
    parts = []
    s = strip_other_lang(sheet, lang)

    # masthead
    m = re.search(r'<div class="masthead">(.*?)</div>\s*</div>', s, flags=re.S)
    if m:
        total += 12.0
        parts.append(("masthead", 12.0))
        s = s.replace(m.group(0), "")

    def add(name, mm):
        nonlocal total
        total += mm
        parts.append((name, mm))

    # consume top-level blocks in order
    for m in re.finditer(r"<h1\b[^>]*>(.*?)</h1>", s, flags=re.S):
        add("h1", wrapped_mm(text_of(m.group(1)), WIDTH_MM * 0.89, H1_PT, 1.20, True)
            + 5 * MM_PER_PT)
    for m in re.finditer(r'<p class="byline">(.*?)</p>', s, flags=re.S):
        add("byline", wrapped_mm(text_of(m.group(1)), WIDTH_MM, 8.6, 1.4, True) + 2 * MM_PER_PT)
    for m in re.finditer(r'<p class="affil">(.*?)</p>', s, flags=re.S):
        add("affil", wrapped_mm(text_of(m.group(1)), WIDTH_MM, 7.6, 1.45, True) + 9 * MM_PER_PT)
    for m in re.finditer(r'<div class="abstract">(.*?)</div>', s, flags=re.S):
        inner = WIDTH_MM - (10 + 10) * MM_PER_PT - 2.4 * MM_PER_PT
        hh = 0.0
        for p in re.findall(r"<p\b[^>]*>(.*?)</p>", m.group(1), flags=re.S):
            hh += wrapped_mm(text_of(p), inner, ABSTRACT_PT, 1.50)
        add("abstract", hh + (8 + 8 + 10) * MM_PER_PT)
    for m in re.finditer(r"<h2\b[^>]*>(.*?)</h2>", s, flags=re.S):
        mt, mb = margins("h2", 6.5, 3.0)
        add("h2", wrapped_mm(text_of(m.group(1)), WIDTH_MM, H2_PT, 1.3, True)
            + (mt + mb + 2.0) * MM_PER_PT)
    for m in re.finditer(r"<h3\b[^>]*>(.*?)</h3>", s, flags=re.S):
        mt, mb = margins("h3", 7.0, 2.5)
        add("h3", wrapped_mm(text_of(m.group(1)), WIDTH_MM, H3_PT, 1.3, True)
            + (mt + mb) * MM_PER_PT)
    for m in re.finditer(r'<div class="flow">(.*?)</div>', s, flags=re.S):
        inner = WIDTH_MM - 12 * MM_PER_PT
        add("flow", wrapped_mm(text_of(m.group(1)), inner, 8.4, 1.75, True)
            + (5 + 5 + 12) * MM_PER_PT)
    for m in re.finditer(r'<div class="key[^"]*">(.*?)</div>', s, flags=re.S):
        inner = WIDTH_MM - (9 + 9) * MM_PER_PT - 2.4 * MM_PER_PT
        hh = wrapped_mm(text_of(re.search(r'class="h"[^>]*>(.*?)</span>', m.group(1),
                                          flags=re.S).group(1) if 'class="h"' in m.group(1) else "",),
                        inner, 7.8, 1.3, True) + 2.4 * MM_PER_PT
        for p in re.findall(r"<p\b[^>]*>(.*?)</p>", m.group(1), flags=re.S):
            hh += wrapped_mm(text_of(p), inner, KEY_PT, 1.45)
        add("key", hh + (6 + 6 + 6 + 7) * MM_PER_PT)
    for m in re.finditer(r'<div class="seq">(.*?)</div>', s, flags=re.S):
        add("seq", wrapped_mm(text_of(m.group(1)), WIDTH_MM - 12 * MM_PER_PT, 7.5, 1.5, True)
            + (4 + 4 + 8) * MM_PER_PT)
    for m in re.finditer(r'<div class="methods">(.*?)</div>\s*(?=<h2|<div class="pageno")',
                         s, flags=re.S):
        hh = 0.0
        for q in re.findall(r"<p\b[^>]*>(.*?)</p>", m.group(1), flags=re.S):
            hh += wrapped_mm(text_of(q), WIDTH_MM, METHODS_PT, 1.42) + 3 * MM_PER_PT
        add("methods", hh)
    for m in re.finditer(r'<div class="foot">(.*?)</div>', s, flags=re.S):
        add("foot", wrapped_mm(text_of(m.group(1)), WIDTH_MM, 7.3, 1.45, True)
            + (7 + 4) * MM_PER_PT)
    for m in re.finditer(r"<table\b.*?</table>", s, flags=re.S):
        add("table", table_mm(m.group(0), WIDTH_MM, lang, base, BODY_LH))
    for m in re.finditer(r'<ol class="refs">(.*?)</ol>', s, flags=re.S):
        hh = 0.0
        for li in re.findall(r"<li>(.*?)</li>", m.group(1), flags=re.S):
            hh += wrapped_mm(text_of(li), (WIDTH_MM - 7) / 2, 7.8, 1.42) + 2.6 * MM_PER_PT
        add("refs", hh / 2 + 3 * MM_PER_PT)
    # body paragraphs and list items that are not inside the structures above
    rest = s
    for pat in [r'<div class="abstract">.*?</div>', r'<div class="key[^"]*">.*?</div>',
                r'<div class="flow">.*?</div>', r'<div class="seq">.*?</div>',
                r'<div class="foot">.*?</div>', r"<table\b.*?</table>",
                r'<ol class="refs">.*?</ol>', r"<figcaption>.*?</figcaption>",
                r'<p class="byline">.*?</p>', r'<p class="affil">.*?</p>',
                r'<div class="methods">.*?</div>\s*(?=<h2|<div class="pageno")']:
        rest = re.sub(pat, "", rest, flags=re.S)
    for m in re.finditer(r'<p\b(?![^>]*class="(?:byline|affil)")[^>]*>(.*?)</p>', rest, flags=re.S):
        add("p", wrapped_mm(text_of(m.group(1)), WIDTH_MM, base, BODY_LH) + P_MB * MM_PER_PT)
    for m in re.finditer(r"<ul>(.*?)</ul>", rest, flags=re.S):
        hh = 0.0
        for li in re.findall(r"<li\b[^>]*>(.*?)</li>", m.group(1), flags=re.S):
            hh += wrapped_mm(text_of(li), WIDTH_MM - 11 * MM_PER_PT, base, BODY_LH) \
                + 1.8 * MM_PER_PT
        add("ul", hh + P_MB * MM_PER_PT)
    # figures: walk the figpage in document order so a figrow counts once,
    # at half width, with the row taking the height of its taller member
    fp = re.search(r'<div class="figpage">(.*)$', s, flags=re.S)
    if fp:
        body_fp = fp.group(1)
        rows = []
        for m in re.finditer(r'<div class="figrow">(.*?)</div>\s*(?=<figure|<div class="pageno")',
                             body_fp, flags=re.S):
            rows.append((m.start(), m.end(), m.group(1), "row"))
        taken = [(a, b) for a, b, _, _ in rows]
        for m in re.finditer(r"<figure([^>]*)>(.*?)</figure>", body_fp, flags=re.S):
            if any(a <= m.start() < b for a, b in taken):
                continue
            rows.append((m.start(), m.end(), m.group(0), "fig"))
        rows.sort()

        def one_fig(chunk, w):
            inner = re.search(r"<figure[^>]*>(.*?)</figure>", chunk, flags=re.S)
            inner = inner.group(1) if inner else chunk
            tok = re.search(r"\{\{(FIG_\w+)\}\}", inner)
            vb = FIGVB.get(tok.group(1)) if tok else None
            hh = (w * vb[1] / vb[0]) if vb else 40.0
            cap = re.search(r"<figcaption>(.*?)</figcaption>", inner, flags=re.S)
            if cap:
                hh += wrapped_mm(text_of(strip_other_lang(cap.group(1), lang)), w,
                                 FIGCAP_PT, 1.42, True) + 2.4 * MM_PER_PT
            return hh

        for _, _, chunk, kind in rows:
            if kind == "row":
                halves = re.findall(r"<figure[^>]*>.*?</figure>", chunk, flags=re.S)
                w = (174.0 - 6 * MM_PER_PT) / 2
                add("figrow", max(one_fig(h, w) for h in halves) + FIGGAP_PT * MM_PER_PT)
            else:
                w = 92.0 if "fignarrow" in chunk[:60] else 174.0
                add("figure", one_fig(chunk, w) + FIGGAP_PT * MM_PER_PT)
    if verbose:
        for n, v in parts:
            print("      %-10s %6.1f mm" % (n, v))
    return total


# figure aspect ratios, read from the rendered SVGs
FIGVB = {}
for tok, name in (("FIG_OVERVIEW", "overview"), ("FIG_REVISIT", "revisit"),
                  ("FIG_MOTIF", "motif"), ("FIG_LOCUS", "locus_slc7a11"),
                  ("FIG_STRUCTURES", "structures")):
    p = os.path.join(BASE, "results", "figs", name + ".light.svg")
    if os.path.exists(p):
        m = re.search(r'viewBox="0 0 ([\d.]+) ([\d.]+)"', open(p).read(600))
        if m:
            FIGVB[tok] = (float(m.group(1)), float(m.group(2)))

html = open(TPL).read()
style = re.search(r"<style>(.*?)</style>", html, flags=re.S).group(1)
STYLE = style
BASE_PT = css_num(re.search(r"\nbody\{font-family[^}]*\}", style).group(0), "font-size", 9.6)
ZH_PT = css_num(re.search(r'body\[data-lang="zh"\]\{[^}]*\}', style).group(0), "font-size", BASE_PT)
_pg = re.search(r"@page \{[^}]*margin:\s*([\d.]+)mm\s+([\d.]+)mm\s+([\d.]+)mm\s+([\d.]+)mm", style)
if _pg:
    MARGIN_TOP_MM, MARGIN_SIDE_MM, MARGIN_BOT_MM = (float(_pg.group(1)), float(_pg.group(2)),
                                                    float(_pg.group(3)))
    AVAIL_MM = PAGE_H_MM - MARGIN_TOP_MM - MARGIN_BOT_MM
    WIDTH_MM = 210.0 - 2 * MARGIN_SIDE_MM
H1_PT = css_num(re.search(r"h1\{[^}]*\}", style).group(0), "font-size", 16.5)
H2_PT = css_num(re.search(r"h2\{[^}]*\}", style).group(0), "font-size", 10.4)
H3_PT = css_num(re.search(r"h3\{[^}]*\}", style).group(0), "font-size", 9.1)
TABLE_PT = css_num(re.search(r"table\{[^}]*\}", style).group(0), "font-size", 8.5)
CAPTION_PT = css_num(re.search(r"caption\{[^}]*\}", style).group(0), "font-size", 8.4)
BODY_LH = 1.44
P_MB = 4.2


def margins(sel, default_top, default_bot):
    m = re.search(re.escape(sel) + r"\{[^}]*margin:\s*([\d.]+)pt\s+0\s+([\d.]+)pt", STYLE)
    return (float(m.group(1)), float(m.group(2))) if m else (default_top, default_bot)
ABSTRACT_PT = css_num(re.search(r"\.abstract\{[^}]*\}", style).group(0), "font-size", 9.3)
METHODS_PT = css_num(re.search(r"\.methods\{[^}]*\}", style).group(0), "font-size", 7.2)
FIGCAP_PT = css_num(re.search(r"figcaption\{[^}]*\}", style).group(0), "font-size", 7.5)
FIGGAP_PT = css_num(re.search(r"\.figpage\{[^}]*\}", style).group(0), "gap", 5.5)
KEY_PT = css_num(re.search(r"\.key\{[^}]*\}", style).group(0), "font-size", 9.0)

body = html[html.index("<body"):]
body = re.sub(r"<script.*?</script>", "", body, flags=re.S)
sheets = re.findall(r'<section class="sheet">(.*?)</section>', body, flags=re.S)

print("available height per page: %.0f mm  (A4 minus %g/%g mm margins)"
      % (AVAIL_MM, MARGIN_TOP_MM, MARGIN_BOT_MM))
print("body %.1f pt · tables %.1f pt · abstract %.1f pt\n" % (BASE_PT, TABLE_PT, ABSTRACT_PT))
worst = 0.0
verbose = "-v" in sys.argv
for lang in ("en", "zh"):
    print("--- %s ---" % lang)
    for i, sh in enumerate(sheets, 1):
        h = page_mm(sh, lang, verbose)
        worst = max(worst, h / AVAIL_MM)
        flag = "OVERFLOW" if h > AVAIL_MM else ("tight" if h > AVAIL_MM * 0.95 else "ok")
        print("  page %d: %6.1f mm  (%3.0f%% of page)  %s" % (i, h, 100 * h / AVAIL_MM, flag))
    print()
print("worst page fill: %.0f%%" % (100 * worst))
sys.exit(1 if worst > 1.0 else 0)
