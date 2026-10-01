#!/usr/bin/env python3
r"""
Render a GitHub-flavoured markdown document with LaTeX math to PDF.

Written for WPMW's docs and used unchanged as the github-md-to-pdf skill's
script; it works on any repository's markdown.  The source file is never
modified: every rewrite below is applied to a temporary copy, which is then
rendered by Pandoc with the xelatex engine.

Math forms
----------
GitHub accepts six.  Pandoc reads ``$...$`` and ``$$...$$`` natively, and
``\(...\)`` / ``\[...\]`` with the ``tex_math_single_backslash`` extension,
which is enabled.  The two forms WPMW actually uses (see "Style guide for
math in WPMW markdown" in this directory's README) are rewritten:

- inline ``$`...`$`` (backtick-guarded, so ``_`` / ``*`` inside the LaTeX are
  not read as emphasis) becomes ``$...$``.  A span may wrap across a source
  line; a match that would cross a blank line (a paragraph break) is left
  alone, since that can only be an unclosed span, which the linter reports.
- fenced ```` ```math ```` blocks become ``$$...$$``.  A fence nested in a
  list item or blockquote carries that context's indentation and/or ``>``
  on every line; the prefix is captured, required on the closing fence, and
  reapplied to every line of the restored block, so the math stays nested.

Doubled-backslash escapes (``\\;``, ``\\{``, ``\\}``, ...) are GitHub's
workaround for its CommonMark preprocessor stripping a lone backslash before
punctuation inside ``$...$``/``$$...$$``.  Real LaTeX reads ``\\`` as a line
break, so they are reduced to single-backslash form -- but only inside
genuine ``$...$``/``$$...$$`` spans.  Content from ```` ```math ```` blocks
and ``$`...`$`` spans is already natural LaTeX and may hold a deliberate
``\\`` (e.g. in ``\begin{cases}``), so it is stashed first and left untouched.

Other GitHub-versus-Pandoc differences
--------------------------------------
- Lists.  A bullet list (or a list starting at 1.) directly under a
  paragraph line, with no blank line, is a list on GitHub; Pandoc would merge
  it into the paragraph.  A blank line is inserted.
- Mermaid.  ```` ```mermaid ```` blocks are drawn to PNG with mermaid-cli
  (``mmdc``, or ``npx -p @mermaid-js/mermaid-cli``), pointing Puppeteer at a
  local Chromium.  If that fails the block stays as code, with a warning.
  ``--no-mermaid`` skips the attempt.
- Images.  ``--image-map PREFIX=DIR`` points image URLs starting with PREFIX
  at local files -- for figures on a branch not yet pushed, or no network.
- Strikethrough.  Pandoc's template loads ``soul`` and emits ``\st{...}``,
  which corrupts or crashes on content containing math.  ``ulem``'s
  ``\sout`` does not, so it is loaded and ``\st`` aliased to it.
- Glyphs.  ✓ ✔ ✗ ✘ are missing from DejaVu Serif and would be silently
  dropped; they are mapped to math symbols with ``newunicodechar`` when that
  package is installed.

Pre-flight and post-flight
--------------------------
Before rendering, the source is scanned for control characters and for a
TAB that ate a LaTeX command.  Both are the signature of a backslash escape
(``\a``, ``\b``, ``\f``, ``\v``, ``\t``) consumed by a non-raw Python string
when the markdown was generated -- e.g. ``\approx`` -> BEL + ``pprox``.
GitHub shows garbage there too, so the render is refused (exit status 2) and
the lines are listed; fix the SOURCE.  ``--allow-control-chars`` overrides.

After rendering, ``pdftotext`` is scanned for leftover ``$```, ```` ```math ````
or raw LaTeX commands, which indicate math that was not typeset.

Requirements
------------
``pandoc``, ``xelatex`` and the TeX packages lmodern, soul, ulem, hyphenat
and microtype; newunicodechar is optional.  On Debian/Ubuntu::

    apt-get install -y pandoc texlive-xetex lmodern texlive-plain-generic \
                       texlive-latex-extra

Usage
-----
    python3 md_to_pdf.py docs/supplement/emission_and_absorption.md
    python3 md_to_pdf.py doc.md -o /tmp/out.pdf
    python3 md_to_pdf.py doc.md -o out.pdf \
        --image-map https://raw.githubusercontent.com/billpage/wpmw/output/figures/=./figs

With no ``-o``, inside the wpmw repository the PDF is written via
``output_path()`` (see "Output path convention" in this directory's README)
using the input's stem; a standalone copy writes ``<stem>.pdf`` beside the
input.

Known limitation
----------------
Equivalent in content to GitHub's preview, not pixel-identical: this uses
LaTeX's own typesetting, not GitHub's browser-side MathJax/KaTeX.  Pixel
fidelity would need a headless browser on the pushed page -- out of scope.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile

try:  # inside the wpmw repository: default output goes through output_path()
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from wpmwlib.wpmw_utils import output_path  # noqa: E402
except ImportError:  # standalone copy (e.g. the github-md-to-pdf skill)
    output_path = None

_HEADER_INCLUDES = r"""
\usepackage[htt]{hyphenat}
\usepackage{microtype}
\emergencystretch=3em
\sloppy
\usepackage[normalem]{ulem}
\let\st\sout
\usepackage{amssymb}
\IfFileExists{newunicodechar.sty}{%
  \usepackage{newunicodechar}%
  \newunicodechar{✓}{\ensuremath{\checkmark}}%
  \newunicodechar{✔}{\ensuremath{\checkmark}}%
  \newunicodechar{✗}{\ensuremath{\times}}%
  \newunicodechar{✘}{\ensuremath{\times}}%
}{}
"""

# Group 1 is whatever leading indentation and/or blockquote markers (">")
# precede the fence -- required to match identically on the closing fence
# via the \1 backreference, so a fence nested in a list item or blockquote
# doesn't let an unrelated later ``` anywhere in the document act as its
# closing marker (see module docstring).
_FENCED_MATH = re.compile(r"^([ \t>]*)```math\n(.*?)\n\1```", re.DOTALL | re.MULTILINE)
_BACKTICK_DOLLAR = re.compile(r"\$`([^`]+?)`\$", re.DOTALL)
_BLOCK_MATH = re.compile(r"\$\$(.+?)\$\$", re.DOTALL)
_INLINE_MATH = re.compile(
    r"(?<![\\$])\$(?![ \t\n$`])([^\n$]+?)(?<![ \t])\$(?![0-9$])"
)
# Two literal backslashes followed by one ASCII-punctuation character --
# GitHub's documented workaround for its own CommonMark backslash-strip.
_DOUBLED_ESCAPE = re.compile(r"\\\\([!\"#$%&'()*+,\-./:;<=>?@\[\]^_`{|}~])")


_MERMAID = re.compile(r"^([ \t>]*)```mermaid[ \t]*\n(.*?)\n\1```[ \t]*$",
                      re.DOTALL | re.MULTILINE)
# A bullet item, or an ordered item numbered 1 -- the only list starts that
# CommonMark lets interrupt a paragraph (GitHub renders them as lists; Pandoc
# needs a blank line first).
_LIST_START = re.compile(r"^\s*(?:[-*+]|1[.)])\s+\S")
_LIST_ITEM = re.compile(r"^\s*(?:[-*+]|\d+[.)])\s+\S")
_CONTROL = re.compile(r"[\x00-\x08\x0b\x0c\x0e-\x1f]")
_EATEN_TAB = re.compile(r"\t(?:hinspace|heta|imes|ilde|ext|extbf|frac|au)\b")


def _reduce_doubled_escapes(expr: str) -> str:
    """\\\\X -> \\X, for real LaTeX. See module docstring."""
    return _DOUBLED_ESCAPE.sub(r"\\\1", expr)


# ------------------------------------------------------------------ pre-flight
def find_control_chars(text: str) -> list[tuple[int, str]]:
    """Lines holding control characters or a TAB that ate a LaTeX command --
    the signature of a backslash escape (\\a \\b \\f \\v \\t) consumed by a
    non-raw Python string when the markdown was generated."""
    hits = []
    for n, line in enumerate(text.split("\n"), 1):
        if _CONTROL.search(line.replace("\t", "")) or _EATEN_TAB.search(line):
            hits.append((n, line.encode("unicode_escape").decode()[:120]))
    return hits


# ------------------------------------------------------------------ rewriting
def separate_lists(text: str) -> str:
    """Insert a blank line between a paragraph line and a list that starts
    directly under it.  GitHub renders such a list; Pandoc would merge it into
    the paragraph.  Skips code fences and multi-line $$ blocks."""
    out: list[str] = []
    in_fence = in_display = False
    for line in text.split("\n"):
        if line.lstrip(" \t>").startswith("```"):
            in_fence = not in_fence
        elif not in_fence and line.count("$$") % 2 == 1:
            in_display = not in_display
        if not in_fence and not in_display and out and _LIST_START.match(line):
            prev = out[-1]
            if (prev.strip() and not _LIST_ITEM.match(prev)
                    and not prev.startswith((" ", "\t", ">", "|", "#"))):
                out.append("")
        out.append(line)
    return "\n".join(out)


def map_images(text: str, mappings: list[tuple[str, str]]) -> str:
    """Point image links whose URL starts with PREFIX at files in DIR."""
    for prefix, local in mappings:
        def sub(m: re.Match) -> str:
            path = os.path.join(local, m.group(2)[len(prefix):])
            if not os.path.exists(path):
                print(f"WARNING: mapped image not found: {path}", file=sys.stderr)
            return f"{m.group(1)}({path}"
        text = re.sub(r"(!\[[^\]]*\])\((" + re.escape(prefix) + r"[^)\s]*)", sub, text)
    return text


def _mmdc_cmd() -> list[str] | None:
    if shutil.which("mmdc"):
        return ["mmdc"]
    if shutil.which("npx"):
        return ["npx", "-y", "-p", "@mermaid-js/mermaid-cli", "mmdc"]
    return None


def _chromium() -> str | None:
    cands = [os.path.join(os.environ.get("PLAYWRIGHT_BROWSERS_PATH", "/opt/pw-browsers"), "chromium")]
    cands += [shutil.which(n) for n in ("chromium", "chromium-browser", "google-chrome")]
    for c in cands:
        if c and os.path.exists(c):
            return c
    return None


def render_mermaid(text: str, workdir: str) -> str:
    """Replace ```mermaid blocks by PNGs drawn with mermaid-cli; on any
    failure leave the block as code and warn."""
    cmd = _mmdc_cmd()
    n = 0

    def sub(m: re.Match) -> str:
        nonlocal n
        n += 1
        pre, body = m.group(1), m.group(2)
        body = "\n".join(l[len(pre):] if l.startswith(pre) else l for l in body.split("\n"))
        if cmd is None:
            print("WARNING: mermaid-cli not available; diagram left as code", file=sys.stderr)
            return m.group(0)
        src = os.path.join(workdir, f"mermaid{n}.mmd")
        png = os.path.join(workdir, f"mermaid{n}.png")
        with open(src, "w", encoding="utf-8") as f:
            f.write(body + "\n")
        args = cmd + ["-i", src, "-o", png, "-s", "3", "-b", "white"]
        chrome = _chromium()
        if chrome:
            cfg = os.path.join(workdir, "puppeteer.json")
            with open(cfg, "w") as f:
                json.dump({"executablePath": chrome, "args": ["--no-sandbox"]}, f)
            args += ["-p", cfg]
        env = dict(os.environ, PUPPETEER_SKIP_DOWNLOAD="1")
        try:
            r = subprocess.run(args, capture_output=True, text=True, env=env, timeout=180)
            ok = r.returncode == 0 and os.path.exists(png)
            err = r.stderr
        except (OSError, subprocess.TimeoutExpired) as e:
            ok, err = False, str(e)
        if not ok:
            print(f"WARNING: mermaid render failed; left as code\n{err[-400:]}", file=sys.stderr)
            return m.group(0)
        return f"{pre}![Diagram {n}]({png})"

    return _MERMAID.sub(sub, text)


def convert_math_delimiters(text: str) -> str:
    """Rewrite the two GitHub math forms Pandoc doesn't parse natively, and
    reduce doubled-backslash escapes within genuine $...$/$$...$$ spans.

    ```` ```math ... ``` ```` fenced blocks and ``` $`...`$ ``` inline spans
    are stashed out first (their content is already natural single-backslash
    LaTeX and must not be touched -- it may contain a deliberate ``\\`` line
    break), then restored verbatim under Pandoc-native ``$$`` / ``$``
    delimiters. What's left after stashing is genuine ``$...$``/``$$...$$``
    content as WPMW's docs actually write it, where a doubled escape is
    reduced to single-backslash form before Pandoc sees it. Already-native
    ``\\(...\\)`` and ``\\[...\\]`` are left untouched (handled by Pandoc's
    ``tex_math_single_backslash`` extension instead).

    A fenced block nested inside a list item or blockquote carries that
    context's marker (indentation and/or ``>``) on every line, including
    the fence lines themselves. That shared prefix is captured and stripped
    from the content before stashing, then reapplied to every line of the
    restored ``$$...$$`` block, so the math stays correctly nested rather
    than being hoisted out of its list item or blockquote.

    A ``` $`...`$ ``` span may itself wrap across a source line break (a
    long expression hand-wrapped in prose); the non-greedy, backtick-free
    content class already stops at the first ``` `$ ```, so this is safe
    for the ordinary case. As a guard against the one degenerate case that
    isn't -- a genuinely unclosed opening ``` $` ``` with no backtick
    before its eventual, unrelated closing marker -- a match that would
    cross a blank line (a paragraph break, where CommonMark could not have
    intended a single span to continue) is left untouched rather than
    stashed, matching this project's own linter, which already flags an
    unclosed ``` $`...`$ ``` as a lint error rather than something this
    script should try to repair.
    """
    fenced_store: list[tuple[str, str]] = []

    def _stash_fenced(m: re.Match) -> str:
        prefix, content = m.group(1), m.group(2)
        if prefix:
            content = "\n".join(
                line[len(prefix):] if line.startswith(prefix) else line
                for line in content.split("\n")
            )
        fenced_store.append((prefix, content))
        return f"{prefix}\x00FENCED{len(fenced_store) - 1}\x00"

    text = _FENCED_MATH.sub(_stash_fenced, text)

    backtick_store: list[str] = []

    def _stash_backtick(m: re.Match) -> str:
        content = m.group(1)
        if re.search(r"\n[ \t]*\n", content):
            return m.group(0)
        backtick_store.append(content)
        return f"\x00BACKTICK{len(backtick_store) - 1}\x00"

    text = _BACKTICK_DOLLAR.sub(_stash_backtick, text)

    text = _BLOCK_MATH.sub(
        lambda m: f"$${_reduce_doubled_escapes(m.group(1))}$$", text)
    text = _INLINE_MATH.sub(
        lambda m: f"${_reduce_doubled_escapes(m.group(1))}$", text)

    for i, (prefix, content) in enumerate(fenced_store):
        restored = f"$$\n{content.strip(chr(10))}\n$$"
        if prefix:
            restored = "\n".join(prefix + line for line in restored.split("\n"))
        text = text.replace(f"{prefix}\x00FENCED{i}\x00", restored)
    for i, content in enumerate(backtick_store):
        text = text.replace(f"\x00BACKTICK{i}\x00", f"${content}$")

    return text


# ------------------------------------------------------------------ post-flight
_POST_PAT = re.compile(r"\$`|```math|\\(?:frac|hbar|thinspace|mathrm|sqrt|partial|begin|left|right)\b")


def _literal_code(src_text: str) -> str:
    """The source's code -- fenced blocks other than math/mermaid, and inline
    code spans other than $`...`$ math -- where the post-flight's patterns
    may legitimately appear as text (e.g. a style guide quoting the syntax)."""
    blocks = re.findall(r"^[ \t>]*(`{3,})(?!math|mermaid)[^\n]*\n(.*?)^[ \t>]*\1",
                        src_text, flags=re.DOTALL | re.MULTILINE)
    spans = re.findall(r"(?<![$`])(`+)(.+?)(?<!`)\1(?![`$])", src_text)
    return "\n".join(b for _, b in blocks) + "\n" + "\n".join(c for _, c in spans)


def post_check(pdf: str, src_text: str = "") -> list[str]:
    """Lines of the PDF's text that look like untypeset math.  A pattern is
    reported only if it occurs more often in the PDF than in the source's own
    code, so documents that quote the syntax are not flagged."""
    if not shutil.which("pdftotext"):
        return []
    txt = subprocess.run(["pdftotext", pdf, "-"], capture_output=True, text=True).stdout
    allowed: dict[str, int] = {}
    for m in _POST_PAT.finditer(_literal_code(src_text)):
        allowed[m.group(0)] = allowed.get(m.group(0), 0) + 1
    found: dict[str, list[str]] = {}
    for line in txt.split("\n"):
        for m in _POST_PAT.finditer(line):
            found.setdefault(m.group(0), []).append(line)
    out: list[str] = []
    for tok, lines in found.items():
        if len(lines) > allowed.get(tok, 0):
            out += lines[:3]
    return out[:10]


def prepare(text: str, workdir: str, image_map=(), mermaid: bool = True) -> str:
    """All source-to-Pandoc rewriting, in order."""
    if mermaid:
        text = render_mermaid(text, workdir)
    text = map_images(text, list(image_map))
    text = separate_lists(text)
    return convert_math_delimiters(text)


def render_pdf(src_path: str, dst_path: str, timeout: int = 180, image_map=(),
               toc: bool = True, mermaid: bool = True,
               allow_control_chars: bool = False) -> None:
    """Convert ``src_path`` (GitHub-flavoured math markdown) to a PDF at
    ``dst_path``.  Exits 2 if the pre-flight finds control characters."""
    with open(src_path, encoding="utf-8") as f:
        text = f.read()
    bad = find_control_chars(text)
    if bad:
        print("Control characters / eaten escapes in the source -- fix the SOURCE:",
              file=sys.stderr)
        for n, l in bad:
            print(f"  line {n}: {l}", file=sys.stderr)
        if not allow_control_chars:
            sys.exit(2)

    with tempfile.TemporaryDirectory() as td:
        prepped = prepare(text, td, image_map, mermaid)
        md_path = os.path.join(td, "prepped.md")
        header_path = os.path.join(td, "header.tex")
        with open(md_path, "w", encoding="utf-8") as f:
            f.write(prepped)
        with open(header_path, "w", encoding="utf-8") as f:
            f.write(_HEADER_INCLUDES)
        cmd = [
            "pandoc", md_path, "-o", dst_path,
            "-f", "markdown+tex_math_single_backslash",
            "--pdf-engine=xelatex",
            "-V", "mainfont=DejaVu Serif",
            "-V", "monofont=DejaVu Sans Mono",
            "-V", "geometry:margin=1in",
            "-V", "linkcolor=blue",
            "-H", header_path,
            "--resource-path", os.path.dirname(os.path.abspath(src_path)),
        ]
        if toc:
            cmd.append("--toc")
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        for l in r.stderr.splitlines():
            if any(k in l for k in ("Missing character", "Could not fetch", "rror")):
                print(l, file=sys.stderr)
        if r.returncode:
            miss = re.findall(r"File `([^']+\.sty)' not found", r.stderr)
            if miss:
                print(f"Missing TeX package(s): {miss} -- see the install line above",
                      file=sys.stderr)
            raise subprocess.CalledProcessError(r.returncode, cmd, r.stdout, r.stderr)

    left = post_check(dst_path, text)
    if left:
        print("WARNING: possible untypeset math in the PDF text:", file=sys.stderr)
        for l in left:
            print("  " + l[:120], file=sys.stderr)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    parser.add_argument("input", help="path to the source markdown file")
    parser.add_argument("-o", "--output",
                        help="destination PDF (default: output_path(<stem>.pdf) inside "
                             "wpmw, else <input stem>.pdf beside the input)")
    parser.add_argument("--image-map", action="append", default=[], metavar="PREFIX=DIR",
                        help="rewrite image URLs starting with PREFIX to files in DIR "
                             "(repeatable)")
    parser.add_argument("--no-toc", action="store_true", help="omit the table of contents")
    parser.add_argument("--no-mermaid", action="store_true",
                        help="leave ```mermaid blocks as code")
    parser.add_argument("--allow-control-chars", action="store_true",
                        help="render despite the control-character pre-flight")
    args = parser.parse_args()

    maps = [tuple(m.split("=", 1)) for m in args.image_map]
    stem = os.path.splitext(os.path.basename(args.input))[0]
    if args.output:
        dst = args.output
    elif output_path is not None:
        dst = output_path(f"{stem}.pdf")
    else:
        dst = os.path.splitext(args.input)[0] + ".pdf"

    render_pdf(args.input, dst, image_map=maps, toc=not args.no_toc,
               mermaid=not args.no_mermaid,
               allow_control_chars=args.allow_control_chars)
    print(f"Wrote {dst}")


if __name__ == "__main__":
    main()
