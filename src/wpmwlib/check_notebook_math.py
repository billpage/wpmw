"""
Jupyter-notebook LaTeX-math linter for the WPMW project.

``docs/*.md`` files use GitHub's backtick-dollar convention, ``$`...`$`` for
inline math and fenced ```` ```math ... ``` ```` blocks for display math --
adopted specifically to dodge CommonMark's backslash-stripping when GitHub
renders a *markdown* file (see the module docstring of
``wpmwlib.check_md_math`` for the details and the incident that motivated
it).

``docs/notebooks/*.ipynb`` files must NOT use that convention. GitHub's
notebook viewer renders each markdown cell through a plain-MathJax pipeline
that does not know about either GitHub-.md-specific extension: the
backticks inside ``$`...`$`` are parsed as a literal code span (rendering as
monospace text with visible backslashes) and a ```` ```math ```` fence is
just an unrecognised code-block language tag (rendering as a literal code
block, complete with the word "math" as its first line). The fix is the
plain form, ``$...$`` and ``$$...$$`` -- and, because a notebook's markdown-
cell source is stored as a raw JSON string with no CommonMark pass in
between, none of the backslash-stripping problems the backtick-dollar form
exists to dodge apply here in the first place. (First hit: the
``emission_and_absorption.ipynb`` notebook, built by copying the tutorial's
markdown verbatim into cells, 2026-09.)

This linter checks the opposite direction from ``check_md_math``: it flags
any ``$`...`$`` or ```` ```math ```` still present in a notebook's markdown
cells (a sign a doc excerpt was copied in without going through the
``convert_math_for_notebook`` conversion), plus a basic sanity check that
``$`` and ``$$`` delimiters are balanced within each cell. It does not need
``check_md_math``'s GFM-escape or render passes -- those exist entirely to
work around the .md-specific backslash-stripping and fenced-math handling
that notebooks don't go through.

Run as::

    PYTHONPATH=src python3 -m wpmwlib.check_notebook_math docs/notebooks
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

_BACKTICK_DOLLAR = re.compile(r"\$`[^`]*?`\$")
_FENCED_MATH = re.compile(r"```math\b")
_DISPLAY_MATH = re.compile(r"\$\$.*?\$\$", re.DOTALL)


@dataclass
class Issue:
    file: Path
    cell: int
    severity: str
    message: str
    preview: str


def _cell_source(cell: dict) -> str:
    src = cell.get("source", "")
    return "".join(src) if isinstance(src, list) else src


def _preview(text: str, at: int, width: int = 70) -> str:
    lo, hi = max(0, at - width // 2), min(len(text), at + width // 2)
    snippet = text[lo:hi].replace("\n", " ↵ ")
    return snippet


def scan_notebook(path: Path) -> list[Issue]:
    issues: list[Issue] = []
    try:
        nb = json.loads(path.read_text())
    except (json.JSONDecodeError, OSError) as exc:
        return [Issue(path, -1, "READ", f"could not parse notebook: {exc}", "")]

    for i, cell in enumerate(nb.get("cells", [])):
        if cell.get("cell_type") != "markdown":
            continue
        text = _cell_source(cell)

        m = _BACKTICK_DOLLAR.search(text)
        if m:
            issues.append(Issue(
                path, i, "WRONGCONV",
                "found $`...`$ (the .md-file convention) in a notebook "
                "markdown cell -- convert to $...$ before committing",
                _preview(text, m.start())))

        m = _FENCED_MATH.search(text)
        if m:
            issues.append(Issue(
                path, i, "WRONGCONV",
                "found a fenced ```math block (the .md-file convention) in "
                "a notebook markdown cell -- convert to $$...$$ before "
                "committing",
                _preview(text, m.start())))

        # Balance check: remove paired $$...$$ blocks, then every remaining
        # single $ should pair up within the same cell.
        stripped = _DISPLAY_MATH.sub("", text)
        n_single = stripped.count("$")
        if n_single % 2 != 0:
            issues.append(Issue(
                path, i, "UNBALANCED",
                f"odd number of single $ ({n_single}) after removing "
                "$$...$$ blocks -- an inline math span is probably left "
                "unclosed",
                _preview(text, 0)))

    return issues


def _walk_ipynb(paths: Iterable[Path]) -> list[Path]:
    out: list[Path] = []
    for p in paths:
        if p.is_file() and p.suffix == ".ipynb":
            out.append(p)
        elif p.is_dir():
            out.extend(sorted(p.rglob("*.ipynb")))
    return out


def _format_report(issues: list[Issue], all_files: list[Path]) -> str:
    lines: list[str] = []
    by_file: dict[Path, list[Issue]] = {}
    for it in issues:
        by_file.setdefault(it.file, []).append(it)
    for nb_path in all_files:
        rel = nb_path
        try:
            rel = nb_path.relative_to(Path.cwd())
        except ValueError:
            pass
        if nb_path not in by_file:
            lines.append(f"OK   {rel}")
            continue
        lines.append(f"=== {rel} ({len(by_file[nb_path])} issue(s)) ===")
        for it in by_file[nb_path]:
            lines.append(f"  cell {it.cell:3d} [{it.severity:10s}] {it.message}")
            lines.append(f"            {it.preview}")
        lines.append("")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m wpmwlib.check_notebook_math",
        description="Lint WPMW notebooks for the plain-MathJax math "
                    "delimiters GitHub's notebook renderer expects "
                    "(as opposed to the $`...`$ / ```math convention used "
                    "in docs/*.md).",
    )
    parser.add_argument("paths", nargs="*", default=["docs/notebooks"],
                        help="Files or directories to scan "
                             "(default: docs/notebooks).")
    args = parser.parse_args(argv)

    paths = [Path(p) for p in args.paths]
    missing = [p for p in paths if not p.exists()]
    if missing:
        print(f"error: paths not found: {', '.join(map(str, missing))}",
              file=sys.stderr)
        return 2

    all_files = _walk_ipynb(paths)
    issues: list[Issue] = []
    for nb_path in all_files:
        issues.extend(scan_notebook(nb_path))

    report = _format_report(issues, all_files)
    if report:
        print(report)
    print(f"Summary: {len(issues)} issue(s) across {len(all_files)} "
          f"notebook(s).")
    return 0 if not issues else 1


if __name__ == "__main__":
    sys.exit(main())
