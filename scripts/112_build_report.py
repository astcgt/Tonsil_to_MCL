#!/usr/bin/env python3
"""Assemble the assignment report into two self-contained single-language files.

docs/report_template.html carries both languages, every translated element
tagged with a lang attribute. This script cuts one language out of the markup
(not merely hides it with CSS), drops the language toggle and its script, pins
<html lang> and <body data-lang>, then substitutes each {{FIG_*}} token with the
light-mode SVG inlined. The two results need no network, no sibling directory
and no JavaScript:

  REPORT_assignment.html      English - the version to hand in
  REPORT_assignment.zh.html   中文    - the same report, for reading

Run:  env/tonsilmcl/bin/python -I scripts/112_build_report.py
"""
import os
import re
import sys
from html.parser import HTMLParser

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIGS = os.path.join(BASE, "results", "figs")
TPL = os.path.join(BASE, "docs", "report_template.html")

FIG = {
    "FIG_OVERVIEW": "overview.light.svg",
    "FIG_REVISIT": "revisit.light.svg",
    "FIG_MOTIF": "motif.light.svg",
    "FIG_LOCUS": "locus_slc7a11.light.svg",
    "FIG_STRUCTURES": "structures.light.svg",
}

# lang attribute in the template -> (output file, data-lang, <html lang>, <title>)
BUILDS = {
    "en": ("REPORT_assignment.html", "en", "en",
           "Differentiation-matched analysis of mantle cell lymphoma"),
    "zh-Hant": ("REPORT_assignment.zh.html", "zh", "zh-Hant",
                "套細胞淋巴瘤的分化配對分析"),
}

VOID = {"area", "base", "br", "col", "embed", "hr", "img", "input",
        "link", "meta", "param", "source", "track", "wbr"}
TAG = re.compile(r"<(/?)([A-Za-z][A-Za-z0-9]*)\b[^>]*?(/?)>", re.S)
CJK = re.compile(r"[\u2e80-\u9fff\uf900-\ufaff\ufe30-\ufe4f\uff00-\uffef]")


def close_of(html, name, i):
    """Index just past the tag that closes <name>, opened before index i.

    Counts nesting, so an element holding a tag of its own name - a
    <span lang="en"> wrapping a <span class="gene"> - is cut in the right place.
    """
    depth = 1
    low = name.lower()
    for t in TAG.finditer(html, i):
        if t.group(2).lower() != low:
            continue
        if t.group(1) == "/":
            depth -= 1
            if depth == 0:
                return t.end()
        elif t.group(3) != "/" and low not in VOID:
            depth += 1
    sys.exit("unclosed <%s> in the template near offset %d" % (name, i))


def cut_lang(html, other):
    """Delete every element carrying lang="<other>". The document root is kept.

    The lookbehind matters: a bare \\b would let lang="en" match the
    data-lang="en" on <body> and delete the whole document.
    """
    pat = re.compile(r'<([A-Za-z][A-Za-z0-9]*)\b[^>]*?(?<![-\w])lang="%s"[^>]*?>'
                     % re.escape(other))
    pos, cuts = 0, 0
    while True:
        m = pat.search(html, pos)
        if m is None:
            return html, cuts
        if m.group(1).lower() in ("html", "body"):
            pos = m.end()
            continue
        html = html[:m.start()] + html[close_of(html, m.group(1), m.end()):]
        pos, cuts = m.start(), cuts + 1


class Balance(HTMLParser):
    """Minimal well-formedness check: every non-void element is closed in order."""

    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.stack, self.errors = [], []

    def handle_starttag(self, tag, attrs):
        if tag not in VOID:
            self.stack.append((tag, self.getpos()[0]))

    def handle_endtag(self, tag):
        if tag in VOID:
            return
        if not self.stack:
            self.errors.append("line %d: </%s> with nothing open" % (self.getpos()[0], tag))
        elif self.stack[-1][0] != tag:
            self.errors.append("line %d: </%s> closes <%s> opened on line %d"
                               % (self.getpos()[0], tag, *self.stack[-1]))
            self.stack.pop()
        else:
            self.stack.pop()


def inline(path):
    svg = open(path).read().strip()
    # the per-mode files carry a hard background so they read on a dark page;
    # on paper the sheet is already white, so drop it.
    svg = re.sub(r'\s*style="background:#[0-9a-fA-F]+"', "", svg, count=1)
    # let CSS set the size; keep the viewBox so the aspect ratio is preserved
    svg = re.sub(r'^(<svg[^>]*?)\s+width="\d+"\s+height="\d+"', r"\1", svg, count=1)
    if "<svg" not in svg or "viewBox" not in svg.split(">", 1)[0]:
        sys.exit("not a sizeable SVG: %s" % path)
    return svg


def build(tpl, lang):
    out_name, data_lang, html_lang, title = BUILDS[lang]
    other = next(k for k in BUILDS if k != lang)

    html, cuts = cut_lang(tpl, other)

    # the toggle, its script and its stylesheet exist only to switch languages
    html, n = re.subn(r'\n<div class="toolbar".*?</div>\n', "\n", html, flags=re.S)
    if n != 1:
        sys.exit("language toolbar not found in the template")
    html, n = re.subn(r"\n<script>.*?</script>\n", "\n", html, flags=re.S)
    if n != 1:
        sys.exit("language script not found in the template")
    for dead in (r"\n/\*[^\n]*toolbar[^\n]*\*/\n(?:(?!\n)[^\n]*\n)+",
                 r"\n  \.toolbar\{display:none !important\}",
                 r'\nbody\[data-lang="[a-z]+"\] \[lang="[A-Za-z-]+"\]\{display:none\}'):
        html, n = re.subn(dead, "\n" if dead.endswith("+") else "", html)
        if not n:
            sys.exit("dead language CSS not found: %s" % dead)

    # pin what the script used to set, so the file needs no JavaScript and no
    # stale localStorage can hide the one language it contains
    html, n = re.subn(r'<html lang="[^"]*">', '<html lang="%s">' % html_lang, html, count=1)
    html, m = re.subn(r'<body data-lang="[^"]*">', '<body data-lang="%s">' % data_lang,
                      html, count=1)
    if n != 1 or m != 1:
        sys.exit("could not pin <html lang> / <body data-lang>")
    html = re.sub(r"<title>.*?</title>", "<title>%s</title>" % title, html, count=1, flags=re.S)

    bal = Balance()
    bal.feed(html)
    if bal.errors or bal.stack:
        sys.exit("%s: unbalanced markup after the language cut\n  %s"
                 % (out_name, "\n  ".join(bal.errors + ["line %d: <%s> never closed" % (l, t)
                                                        for t, l in bal.stack])))

    for token, name in FIG.items():
        p = os.path.join(FIGS, name)
        if not os.path.exists(p):
            sys.exit("missing figure: %s\n  run scripts/103_figures.py and "
                     "scripts/111_fig_revisit.py" % p)
        html = html.replace("{{%s}}" % token, inline(p))

    problems = []
    if re.findall(r"\{\{[A-Z_]+\}\}", html):
        problems.append("unsubstituted tokens: %s" % re.findall(r"\{\{[A-Z_]+\}\}", html))
    if 'lang="%s"' % other in html:
        problems.append('lang="%s" still present' % other)
    if lang == "en":
        leak = CJK.findall(html)
        if leak:
            problems.append("CJK characters in the English file: %r" % "".join(leak[:40]))
    n_sheet, n_svg = html.count('class="sheet"'), html.count("<svg")
    if n_sheet != 5:
        problems.append("%d sheets, expected 5" % n_sheet)
    if n_svg != len(FIG):
        problems.append("%d SVGs, expected %d" % (n_svg, len(FIG)))
    if problems:
        sys.exit("%s:\n  %s" % (out_name, "\n  ".join(problems)))

    out = os.path.join(BASE, out_name)
    with open(out, "w") as fh:
        fh.write(html)
    print("%-26s %4.0f kB  %d pages  %d figures  %d elements cut"
          % (out_name, len(html.encode()) / 1024, n_sheet, n_svg, cuts))


def main():
    tpl = open(TPL).read()
    for lang in BUILDS:
        build(tpl, lang)


if __name__ == "__main__":
    main()
