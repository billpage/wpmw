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
"""

_FENCED_MATH = re.compile(r"```math\n(.*?)\n```", re.DOTALL)
_BACKTICK_DOLLAR = re.compile(r"\$`([^`\n]+?)`\$")
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
    """
    fenced_store: list[str] = []

    def _stash_fenced(m: re.Match) -> str:
        fenced_store.append(m.group(1))
        return f"\x00FENCED{len(fenced_store) - 1}\x00"

    text = _FENCED_MATH.sub(_stash_fenced, text)

    backtick_store: list[str] = []

    def _stash_backtick(m: re.Match) -> str:
        backtick_store.append(m.group(1))
        return f"\x00BACKTICK{len(backtick_store) - 1}\x00"

    text = _BACKTICK_DOLLAR.sub(_stash_backtick, text)

    text = _BLOCK_MATH.sub(
        lambda m: f"$${_reduce_doubled_escapes(m.group(1))}$$", text)
    text = _INLINE_MATH.sub(
        lambda m: f"${_reduce_doubled_escapes(m.group(1))}$", text)

    for i, content in enumerate(fenced_store):
        text = text.replace(f"\x00FENCED{i}\x00",
                             f"$$\n{content.strip(chr(10))}\n$$")
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
