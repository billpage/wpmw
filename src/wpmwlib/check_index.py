"""
Index-consistency linter for the WPMW project.

``docs/analysis/INDEX.md`` lists every labelled result (theorem, proposition,
lemma, corollary, definition, postulate, conjecture) and every labelled open
item in the analysis notes, with a one-line statement and a status for each.
That is a second copy of facts the notes already state, so it drifts the
moment a note gains a theorem and the index does not. This linter is the
guard, in the same spirit as ``check_abstracts`` for the README ladder: it
re-derives the set of labelled results and open items from the notes and
checks the index against it.

What it treats as a *definition* in a note (``docs/analysis/*.md`` except
``README.md`` and ``INDEX.md``):

* a paragraph that opens with a bold span of the form
  ``**Theorem C1.**``, ``**Corollary K1.1 (the lattice cannot ...).**``,
  ``**Definition (H).**``, ``**Postulate (S) -- sea carrier lock.**`` (also
  inside a blockquote). The span must end in a full stop, which is what
  separates a definition (``**Theorem C1.**``) from a mention at the start of a
  sentence (``**Theorem E7** shows ...``);
* the postulate shorthands ``**P0 (ontology).**`` and ``**(E) Existence.**``;
* ``**No-go lemma.**``;
* a heading of the form ``## 3. Theorem Z3: ...`` or ``## 2. The channel
  table (Lemma C1)``.

Open items are the bullets that start ``- **ID`` inside a section whose
heading contains "Open items", where ID looks like ``K-LS1``, ``S-SP3`` or
``CLS3``. Numbered open items in the older notes have no ID and are not
enforced; they may appear in the index but nothing requires it.

What it checks against the index:

* MISSING  -- a definition or open item in a note has no row in the index.
* BROKEN   -- any relative link in the index (in a row or anywhere else)
  points at a file that does not exist, or at an anchor that no heading in
  that file produces.
* UNMATCHED (warning) -- a row that matches no definition or open item, e.g. a
  hand-entered row for something the extractor cannot see, or a row left behind
  after a result was renamed.

Rows are recognised by their first cell, which must be a link:
``[Thm C1](file.md#anchor)`` for a result and ``[K-LS1](file.md#anchor)`` for an
open item. Type abbreviations: Thm, Prop, Lem, Cor, Def, Post, Conj.

Run as::

    PYTHONPATH=src python3 -m wpmwlib.check_index docs/analysis
    PYTHONPATH=src python3 -m wpmwlib.check_index docs/analysis --suggest

``--suggest`` prints a skeleton index row for everything MISSING, with the
section number and anchor filled in, so adding a note's entries is a matter of
completing the statement and standing columns.

Exits 0 if there are no errors (warnings do not count), 1 if any error is
found, 2 on bad input.
"""

from __future__ import annotations

import argparse
import re
import sys
import unicodedata
from dataclasses import dataclass
from pathlib import Path

INDEX_NAME = "INDEX.md"
SKIP = {"README.md", INDEX_NAME}

ABBR = {
    "Theorem": "Thm", "Proposition": "Prop", "Lemma": "Lem",
    "Corollary": "Cor", "Definition": "Def", "Postulate": "Post",
    "Conjecture": "Conj",
}
ABBR_INV = {v: k for k, v in ABBR.items()}
_TYPES = "|".join(ABBR)

_LABEL = r"(?:[A-Z]+[0-9]*|[0-9]+)(?:\.[0-9]+)?[a-z]?"
_SPAN = re.compile(
    rf"^(?P<type>{_TYPES})\s+(?:\((?P<name>[^)]+)\)|(?P<label>{_LABEL}))"
    r"(?P<tail>.*)$"
)
_TAIL = re.compile(r"(?:\s*\([^)]*\))?(?:\s*[\u2014\u2013,:-]\s*[^.]*)?\.")
_BOLD = re.compile(r"^(?:>\s*)?\*\*(?P<span>[^*]+?)\*\*")
_POST_P = re.compile(r"^P[0-9]+(?=\s+\()")
_POST_LETTER = re.compile(r"^\((?P<l>[A-Z])\)\s+[A-Z][a-z]+\.$")
_HEAD = re.compile(r"^(?P<h>#{1,6})\s+(?P<text>.+?)\s*$")
_HEAD_A = re.compile(
    rf"^(?P<type>{_TYPES})\s+(?P<label>{_LABEL})(?:\s+\([^)]*\))?\s*(?::|$)"
)
_HEAD_B = re.compile(
    r"\((?P<type>Theorem|Proposition|Lemma)\s+(?P<label>" + _LABEL + r")\)\s*$"
)
_ID = re.compile(r"^[A-Z]{1,4}(?:-[A-Z]{2})?[0-9]+[a-z]?$")
_ROW = re.compile(
    r"^\|\s*\[(?P<text>[^\]]+)\]\((?P<file>[^)#]+)(?:#(?P<anchor>[^)]+))?\)\s*\|"
)


@dataclass(frozen=True)
class Item:
    file: str
    kind: str          # "result" or "open"
    key: tuple         # ("Thm", "C1") for results, ("K-LS1",) for open items
    line: int
    section: str       # section number of the enclosing heading, or ""
    anchor: str | None  # GitHub anchor of the enclosing heading, if safe
    text: str          # first words of the statement, for --suggest

    def label(self) -> str:
        if self.kind == "result":
            return f"{self.key[0]} {self.key[1]}"
        return self.key[0]


@dataclass
class Issue:
    severity: str
    where: str
    message: str


# ---------------------------------------------------------------- slugs ---

def github_slug(text: str) -> str:
    """Approximate GitHub's heading-anchor rule: lower-case, drop everything
    except letters, digits, spaces, hyphens and underscores, spaces to
    hyphens."""
    out = []
    for ch in text.lower():
        cat = unicodedata.category(ch)
        if ch in " -_" or cat[0] in "LN" or cat == "Mn":
            out.append(ch)
    return "".join(out).replace(" ", "-")


def heading_slugs(lines: list[str]) -> dict[int, str]:
    """Map line index -> anchor slug, with GitHub's ``-1``, ``-2`` suffixes
    for repeated headings. Fenced code blocks are skipped."""
    seen: dict[str, int] = {}
    result: dict[int, str] = {}
    fence = False
    for i, line in enumerate(lines):
        if line.startswith("```"):
            fence = not fence
            continue
        if fence:
            continue
        m = _HEAD.match(line)
        if not m:
            continue
        slug = github_slug(m.group("text"))
        n = seen.get(slug, 0)
        seen[slug] = n + 1
        result[i] = slug if n == 0 else f"{slug}-{n}"
    return result


def _safe_heading(text: str) -> bool:
    """A heading whose anchor we are confident about: plain ASCII, no math or
    code spans, no emphasis markers."""
    return text.isascii() and not any(c in text for c in "`$*[]<>")


# ----------------------------------------------------------- extraction ---

def _paragraph(lines: list[str], i: int) -> str:
    """Join the paragraph starting at line i (until a blank line), removing
    blockquote markers from continuation lines."""
    parts = [lines[i]]
    j = i + 1
    while j < len(lines) and lines[j].strip():
        parts.append(re.sub(r"^>\s?", "", lines[j]))
        j += 1
    return " ".join(p.strip() for p in parts)


def _first_words(text: str, n: int = 14) -> str:
    text = re.sub(r"\*\*[^*]+?\*\*", "", text, count=1)
    text = re.sub(r"[`$*_>]", "", text)
    return " ".join(text.split()[:n])


def parse_note(path: Path) -> list[Item]:
    lines = path.read_text(encoding="utf8").splitlines()
    slugs = heading_slugs(lines)
    items: dict[tuple, Item] = {}
    section, anchor = "", None
    open_section = False
    fence = False
    for i, line in enumerate(lines):
        if line.startswith("```"):
            fence = not fence
            continue
        if fence:
            continue
        hm = _HEAD.match(line)
        if hm:
            text = hm.group("text")
            num = re.match(r"(\d+(?:\.\d+)*)\.?\s", text + " ")
            section = num.group(1) if num else ""
            anchor = slugs[i] if _safe_heading(text) else None
            if len(hm.group("h")) == 2:
                open_section = bool(re.search(r"open items", text, re.I))
            body = re.sub(r"^\d+(?:\.\d+)*\.?\s+", "", text)
            key = None
            m = _HEAD_A.match(body)
            if m:
                key = (ABBR[m.group("type")], m.group("label"))
            else:
                m = _HEAD_B.search(body)
                if m:
                    key = (ABBR[m.group("type")], m.group("label"))
            if key and ("result", key) not in items:
                items[("result", key)] = Item(
                    path.name, "result", key, i + 1, section, anchor, body)
            continue
        if open_section:
            om = re.match(r"^\s*[-*]\s+(?:~~)?\*\*(?P<id>[A-Za-z0-9-]+?)[.\s(]", line)
            if om and _ID.match(om.group("id")):
                key = (om.group("id"),)
                if ("open", key) not in items:
                    items[("open", key)] = Item(
                        path.name, "open", key, i + 1, section, anchor,
                        _first_words(line))
        if not (line.startswith("**") or line.startswith("> **")):
            continue
        para = _paragraph(lines, i)
        bm = _BOLD.match(para)
        if not bm:
            continue
        span = bm.group("span").strip()
        key = None
        if span == "No-go lemma.":
            key = ("Lem", "(no-go)")
        elif _POST_P.match(span) and span.endswith("."):
            key = ("Post", _POST_P.match(span).group(0))
        elif _POST_LETTER.match(span):
            key = ("Post", f"({_POST_LETTER.match(span).group('l')})")
        else:
            sm = _SPAN.match(span)
            if sm and _TAIL.fullmatch(sm.group("tail")):
                lab = sm.group("label") or f"({sm.group('name')})"
                key = (ABBR[sm.group("type")], lab)
        if key and ("result", key) not in items:
            items[("result", key)] = Item(
                path.name, "result", key, i + 1, section, anchor,
                _first_words(para))
    return list(items.values())


def parse_notes(directory: Path) -> list[Item]:
    items: list[Item] = []
    for path in sorted(directory.glob("*.md")):
        if path.name in SKIP:
            continue
        items.extend(parse_note(path))
    return items


# ---------------------------------------------------------------- index ---

@dataclass
class Row:
    line: int
    text: str
    file: str
    anchor: str | None

    def key(self) -> tuple | None:
        parts = self.text.split(None, 1)
        if len(parts) == 2 and parts[0] in ABBR_INV:
            return ("result", (parts[0], parts[1].strip()))
        if len(parts) == 1 and _ID.match(parts[0]):
            return ("open", (parts[0],))
        return None


def parse_index(path: Path) -> list[Row]:
    rows = []
    for i, line in enumerate(path.read_text(encoding="utf8").splitlines()):
        m = _ROW.match(line)
        if m:
            rows.append(Row(i + 1, m.group("text").strip(), m.group("file"),
                            m.group("anchor")))
    return rows


# ---------------------------------------------------------------- check ---

def check(directory: Path) -> tuple[list[Issue], list[Item]]:
    issues: list[Issue] = []
    index_path = directory / INDEX_NAME
    if not index_path.exists():
        return [Issue("ERROR", str(directory), f"no {INDEX_NAME}")], []

    items = parse_notes(directory)
    rows = parse_index(index_path)

    row_keys: dict[tuple, Row] = {}
    for r in rows:
        k = r.key()
        if k is not None:
            row_keys.setdefault((k, r.file), r)

    missing: list[Item] = []
    for it in items:
        if ((it.kind, it.key), it.file) not in row_keys:
            missing.append(it)
            issues.append(Issue(
                "ERROR", f"{it.file}:{it.line}",
                f"MISSING {it.label()} (section {it.section or '?'}) has no "
                f"row in {INDEX_NAME}"))

    known = {((it.kind, it.key), it.file) for it in items}
    for r in rows:
        k = r.key()
        if k is None:
            issues.append(Issue(
                "WARNING", f"{INDEX_NAME}:{r.line}",
                f"UNMATCHED row label {r.text!r}: not of the form 'Thm X' "
                "or an open-item ID"))
        elif (k, r.file) not in known:
            issues.append(Issue(
                "WARNING", f"{INDEX_NAME}:{r.line}",
                f"UNMATCHED row {r.text!r} in {r.file}: no such definition or "
                "open item found in the note (hand-entered, or renamed?)"))
    issues.extend(check_links(directory, index_path))
    return issues, missing


_LINK = re.compile(r"\]\((?P<target>[^)\s]+)\)")


def check_links(directory: Path, index_path: Path) -> list[Issue]:
    """Every relative link in the index must resolve: the file must exist,
    and a ``#anchor`` into a markdown file must match one of its headings."""
    issues: list[Issue] = []
    slug_cache: dict[Path, set[str]] = {}
    fence = False
    for i, line in enumerate(index_path.read_text(encoding="utf8").splitlines()):
        if line.startswith("```"):
            fence = not fence
            continue
        if fence:
            continue
        for m in _LINK.finditer(line):
            target = m.group("target")
            if re.match(r"^[a-z]+:", target) or target.startswith("#"):
                continue
            path_part, _, anchor = target.partition("#")
            file = (directory / path_part).resolve() if path_part else index_path
            if not file.exists():
                issues.append(Issue("ERROR", f"{INDEX_NAME}:{i + 1}",
                                    f"BROKEN link: {path_part} does not exist"))
                continue
            if anchor and file.suffix == ".md":
                if file not in slug_cache:
                    slug_cache[file] = set(heading_slugs(
                        file.read_text(encoding="utf8").splitlines()).values())
                if anchor not in slug_cache[file]:
                    issues.append(Issue(
                        "ERROR", f"{INDEX_NAME}:{i + 1}",
                        f"BROKEN anchor: #{anchor} is not a heading of "
                        f"{path_part}"))
    return issues


def suggest_row(it: Item) -> str:
    sec = it.section or "?"
    target = f"{it.file}#{it.anchor}" if it.anchor else it.file
    if it.kind == "result":
        return f"| [{it.label()}]({target}) | {sec} | {it.text} | -- |"
    return f"| [{it.label()}]({target}) | {it.text} | open |"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m wpmwlib.check_index",
        description="Check docs/analysis/INDEX.md against the labelled "
                    "results and open items in the analysis notes.")
    parser.add_argument("directory", nargs="?", default="docs/analysis",
                        help="directory holding the notes and INDEX.md "
                             "(default: docs/analysis)")
    parser.add_argument("--suggest", action="store_true",
                        help="print a skeleton row for every MISSING entry")
    args = parser.parse_args(argv)

    directory = Path(args.directory)
    if not directory.is_dir():
        print(f"error: not a directory: {directory}", file=sys.stderr)
        return 2

    issues, missing = check(directory)
    n_items = len(parse_notes(directory))
    errors = [i for i in issues if i.severity == "ERROR"]
    warnings = [i for i in issues if i.severity == "WARNING"]
    for it in issues:
        print(f"[{it.severity:7s}] {it.where}: {it.message}")
    if args.suggest and missing:
        print("\nSuggested rows (complete the statement and standing "
              "columns, then place under the note's heading):\n")
        for it in missing:
            print(suggest_row(it))
    print(f"Summary: {len(errors)} error(s), {len(warnings)} warning(s) "
          f"across {n_items} labelled results and open items.")
    return 0 if not errors else 1


if __name__ == "__main__":
    sys.exit(main())
