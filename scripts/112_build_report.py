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
REPORT = os.path.join(BASE, "REPORT_v3.md")

# Every row of section 6 carries one of these verdict classes. The counts the
# report quotes are derived from the tags, so no number is maintained twice.
# Boundary rule: a pre-registered prediction that fails is 否決, not 撤回 - a
# pre-registration failing is the system working, not a retraction.
CLASSES = ("撤回", "否決", "陷阱", "維持", "無檢定力")
CLASS_TOKEN = {"撤回": "N_RETRACTED", "否決": "N_REFUTED", "陷阱": "N_PITFALL",
               "維持": "N_UPHELD", "無檢定力": "N_NOPOWER"}

# Numbers quoted outside the table, as (file, pattern, what it should equal).
# "total" is every row; the other keys are verdict classes. Each pattern
# captures the figure so it can be compared with section 6, the only place any
# of them is maintained. A pattern that stops matching is an error rather than
# a skipped check: the sentence was reworded, and the gate has to be reworded
# with it instead of quietly stopping.
QUOTED = [
    ("README.md", r"\n(\d+) tests are recorded one per row", "total"),
    ("README.md", r"of which (\d+) are conclusions of our own", "撤回"),
    ("README.md", r"逐列記錄 (\d+) 項檢定", "total"),
    ("README.md", r"其中 (\d+) 項是我們自己發表後", "撤回"),
    ("README.md", r"with all (\d+) entries of the test record \(Chinese\)", "total"),
    ("README.md", r"含全部 (\d+) 項檢定紀錄（中文）", "total"),
    ("docs/hf_model_card.md", r"test record of (\d+) entries", "total"),
    ("docs/hf_model_card.md", r"(\d+) 項檢定紀錄", "total"),
]

ONES = ("zero one two three four five six seven eight nine ten eleven twelve "
        "thirteen fourteen fifteen sixteen seventeen eighteen nineteen").split()
TENS = ("", "", "twenty", "thirty", "forty", "fifty", "sixty", "seventy",
        "eighty", "ninety")

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


def spell(n):
    """English words for 0-99, capitalised: the sentence opens with it."""
    if not 0 <= n < 100:
        sys.exit("cannot spell %d in words; the report sentence assumes under a hundred" % n)
    w = ONES[n] if n < 20 else TENS[n // 10] + ("-" + ONES[n % 10] if n % 10 else "")
    return w[0].upper() + w[1:]


def count_verdicts():
    """Tally section 6 of REPORT_v3.md by verdict class: the one maintained count.

    Returns (total, {class: n}). Every numbered row must carry exactly one of
    CLASSES in its second cell, so a row added without a verdict fails the build
    rather than quietly falling out of the totals.
    """
    md = open(REPORT, encoding="utf-8").read()
    try:
        body = md[md.index("## 六、"):md.index("## 七、")]
    except ValueError:
        sys.exit("cannot locate section 6 of %s, which the report quotes" % REPORT)
    tally, untagged = dict.fromkeys(CLASSES, 0), []
    for line in body.split("\n"):
        m = re.match(r"^\|\s*(\d+)\s*\|\s*([^|]*?)\s*\|", line)
        if not m:
            continue
        if m.group(2) in tally:
            tally[m.group(2)] += 1
        else:
            untagged.append("row %s: %r" % (m.group(1), m.group(2)))
    if untagged:
        sys.exit("section 6 rows without a recognised verdict class (%s):\n  %s"
                 % ("／".join(CLASSES), "\n  ".join(untagged)))
    total = sum(tally.values())
    if total < 2:
        sys.exit("section 6 of %s has %d tagged rows, which cannot be right" % (REPORT, total))
    return total, tally


def check_quoted(total, tally):
    """Fail if any document states a figure that section 6 does not support."""
    want = dict(tally, total=total)
    bad = []
    for rel, pat, key in QUOTED:
        text = open(os.path.join(BASE, rel), encoding="utf-8").read()
        m = re.search(pat, text)
        if m is None:
            bad.append("%s: the sentence this gate watches has been reworded, so the "
                       "figure there is no longer checked. Update QUOTED in "
                       "scripts/112_build_report.py to match it: %s" % (rel, pat))
        elif int(m.group(1)) != want[key]:
            bad.append("%s:%d says %s for %s; section 6 of REPORT_v3.md has %d"
                       % (rel, text[:m.start(1)].count("\n") + 1, m.group(1), key, want[key]))
    if bad:
        sys.exit("figures quoted outside section 6 are out of step:\n  %s" % "\n  ".join(bad))


def main():
    tpl = open(TPL).read()
    total, tally = count_verdicts()
    check_quoted(total, tally)
    subs = [("N_ENTRIES_WORD", spell(total)), ("N_ENTRIES", str(total))]
    subs += [(CLASS_TOKEN[k], str(tally[k])) for k in CLASSES]
    for token, value in subs:
        tpl, k = re.subn(r"\{\{%s\}\}" % token, value, tpl)
        if not k:
            sys.exit("the template no longer carries {{%s}}; put it back or drop this "
                     "substitution" % token)
    for lang in BUILDS:
        build(tpl, lang)
    print("test record: %d entries in REPORT_v3.md section 6 (%s); README and the "
          "model card agree on the total"
          % (total, "、".join("%s %d" % (k, tally[k]) for k in CLASSES)))


if __name__ == "__main__":
    main()
