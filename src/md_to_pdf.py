r"""
Render a WPMW markdown document (with embedded LaTeX math) to PDF.

WPMW's markdown math conventions do not use plain CommonMark math syntax —
see "Style guide for math in WPMW markdown" in this directory's README.
Instead, docs use two GitHub-specific forms that survive CommonMark's
emphasis parser:

- inline:  ``$`...`$``          (backtick-guarded, so ``_`` / ``*`` inside
  the LaTeX are not read as markdown emphasis)
- display: fenced ```` ```math ```` blocks

Pandoc does not recognise either form as math on its own, so this script
rewrites both into Pandoc-native ``$...$`` / ``$$...$$`` in a temporary
copy of the file, then hands that copy to Pandoc with the xelatex engine.
The source document on disk is never modified.

Pandoc is also invoked with the ``tex_math_single_backslash`` extension, so
GitHub's other two math forms — inline ``\(...\)`` and display ``\[...\]``
— render correctly too, in case they ever appear (current WPMW docs use
only the two forms above; this is included for compatibility with GitHub's
full math syntax rather than because it's presently exercised).

This script also reduces GitHub's doubled-backslash escape workaround
(``\\;``, ``\\:``, ``\\{``, ``\\}``, ...) back to single-backslash form
within genuine ``$...$``/``$$...$$`` spans before handing off to Pandoc.
GitHub's CommonMark preprocessor silently strips a lone backslash from
``\X`` (X = ASCII punctuation) inside ``$...$``/``$$...$$`` before MathJax
ever sees it; doubling it (``\\X``) is the documented on-GitHub workaround
(see the style guide below). Real LaTeX (via Pandoc/xelatex) has no such
preprocessor — ``\\`` is its own native line-break command — so a doubled
escape that survives unreduced renders as a spurious line break plus a
literal punctuation character instead of the intended spacing command.
This reduction is scoped to content that was originally ``$...$``/``$$...$$``
in the source; already-native single-backslash LaTeX inside ```` ```math ````
fenced blocks or ``` $`...`$ ``` backtick-dollar spans (which were never
subject to GitHub's strip in the first place, and may contain a *genuine*
intentional ``\\`` line break, e.g. inside ``\begin{cases}...\end{cases}``)
is left completely untouched.

GFM strikethrough (``~~...~~``) containing inline math is also handled.
Pandoc's default LaTeX template auto-loads the ``soul`` package and emits
``\st{...}`` for strikeout content; ``soul`` builds ``\st`` by measuring
and re-kerning each syllable, which corrupts or crashes outright on
content containing a math-mode box (``Extra }, or forgotten $`` deep in
``soul``'s internal macros, from real documents, not a contrived example).
``ulem``'s ``\sout`` does not have this problem -- it draws a straight
strike rule under the content without trying to re-typeset it -- so this
script loads ``ulem`` (with ``normalem``, so it doesn't also redefine
``\emph`` to underline) after ``soul`` and aliases ``\let\st\sout``,
letting Pandoc's own output stand unmodified while replacing which macro
actually draws the strike.

Requires ``pandoc`` and a LaTeX engine (``xelatex``) on PATH, plus the
``lmodern`` font-metrics package. On Debian/Ubuntu::

    apt-get install pandoc texlive-xetex lmodern

Usage
-----
    python3 md_to_pdf.py docs/supplement/emission_and_absorption.md
    python3 md_to_pdf.py docs/supplement/emission_and_absorption.md -o /tmp/out.pdf

With no ``-o``, the PDF is written via ``output_path()`` (see the "Output
path convention" section of this directory's README) using the input
file's stem, e.g. ``emission_and_absorption.pdf``.

Known limitation
-----------------
Equivalent in content to GitHub's preview, not pixel-identical: this uses
LaTeX's own typesetting (fonts, line-breaking, no syntax highlighting),
not a copy of GitHub's browser-side MathJax/KaTeX rendering. For
pixel-fidelity to the GitHub web preview specifically, the file would need
to already be pushed and rendered via a headless browser instead — out of
scope for this script.
"""

from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(__file__))
from wpmwlib.wpmw_utils import output_path  # noqa: E402

_HEADER_INCLUDES = r"""
\usepackage[htt]{hyphenat}
\usepackage{microtype}
\emergencystretch=3em
\sloppy
\usepackage[normalem]{ulem}
\let\st\sout
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


def _reduce_doubled_escapes(expr: str) -> str:
    """\\\\X -> \\X, for real LaTeX. See module docstring."""
    return _DOUBLED_ESCAPE.sub(r"\\\1", expr)


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


def render_pdf(src_path: str, dst_path: str, timeout: int = 90) -> None:
    """Convert ``src_path`` (WPMW markdown) to a PDF at ``dst_path``."""
    with open(src_path, encoding="utf-8") as f:
        prepped = convert_math_delimiters(f.read())

    with tempfile.TemporaryDirectory() as td:
        md_path = os.path.join(td, "prepped.md")
        header_path = os.path.join(td, "header.tex")
        with open(md_path, "w", encoding="utf-8") as f:
            f.write(prepped)
        with open(header_path, "w", encoding="utf-8") as f:
            f.write(_HEADER_INCLUDES)

        subprocess.run(
            [
                "pandoc", md_path, "-o", dst_path,
                "-f", "markdown+tex_math_single_backslash",
                "--pdf-engine=xelatex",
                "-V", "mainfont=DejaVu Serif",
                "-V", "monofont=DejaVu Sans Mono",
                "-V", "geometry:margin=1in",
                "-V", "linkcolor=blue",
                "-H", header_path,
                "--toc",
            ],
            check=True,
            timeout=timeout,
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("input", help="path to the source markdown file")
    parser.add_argument(
        "-o", "--output",
        help="destination PDF path (default: output_path(<stem>.pdf))",
    )
    args = parser.parse_args()

    stem = os.path.splitext(os.path.basename(args.input))[0]
    dst = args.output or output_path(f"{stem}.pdf")

    render_pdf(args.input, dst)
    print(f"Wrote {dst}")


if __name__ == "__main__":
    main()
