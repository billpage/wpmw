"""
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


def convert_math_delimiters(text: str) -> str:
    """Rewrite WPMW's GitHub math conventions to Pandoc-native math.

    ```` ```math ... ``` ```` fenced blocks become ``$$ ... $$``, and
    ``` $`...`$ ``` inline spans become ``$...$``.
    """
    text = re.sub(
        r"```math\n(.*?)\n```",
        lambda m: f"$$\n{m.group(1).strip(chr(10))}\n$$",
        text,
        flags=re.DOTALL,
    )
    text = re.sub(r"\$`([^`]+?)`\$", lambda m: f"${m.group(1)}$", text)
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
