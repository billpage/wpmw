"""
Where is a term used?  A concordance over the project's markdown.

For each term it lists every file (``docs/**/*.md``, ``ORIENTATION.md`` and
the READMEs) that uses it, with a count, the heading each use sits under
(as a working GitHub link), and the first line of context.  It is meant for
the question a reader asks of a word they have just met: where is this used,
and where is it explained?  Uses in a heading, and sentences that look
definitional ("which X", "X is", "we call", "denotes", "means", ...), are
marked ``*`` and listed first, because those are the places that explain a
term rather than use it.

It reads and never writes.  ``GLOSSARY.md`` is the curated companion; this is
the uncurated one, and cannot drift.

Run as::

    PYTHONPATH=src python3 -m wpmwlib.where_used eikonal
    PYTHONPATH=src python3 -m wpmwlib.where_used "dark catalysis" --lines 3
    PYTHONPATH=src python3 -m wpmwlib.where_used eikonal --src     # also src/*.py
    PYTHONPATH=src python3 -m wpmwlib.where_used --candidates 3    # ranked word list

``--candidates N`` prints words that occur in at least N files, ranked by how
much more the notes use them than English does, with a mark for those the
glossary already has.  It needs ``wordfreq``.  Where pip refuses to touch the
system Python (``externally-managed-environment``, as on Linux Mint 22 and
Ubuntu 24.04), use a virtual environment::

    python3 -m venv ~/.venvs/wpmw && ~/.venvs/wpmw/bin/pip install wordfreq
    PYTHONPATH=src ~/.venvs/wpmw/bin/python -m wpmwlib.where_used --candidates 4  It is a menu to
prune by hand, not a verdict.  It cannot see a term used once, in a
specialised sense, in one tutorial: those are for readers to suggest.

Exits 0 if the term was found, 1 if not.
"""

from __future__ import annotations

import argparse
import collections
import re
import signal
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEFINING = re.compile(
    r"(which\s+\w+\.|\bwe call\b|\bcalled\b|\bdenotes?\b|\bmeans?\b|"
    r"\bis the\b|\bare the\b|\bis a\b|\bis an\b|\bnamed?\b|\bsense\b|"
    r"\bdefined?\b|\bdefinition\b)", re.I)


def slug(heading: str) -> str:
    """GitHub's heading anchor."""
    h = re.sub(r"\[([^\]]*)\]\([^)]*\)", r"\1", heading.strip().lower())
    h = re.sub(r"[`*_$\\]", "", h)
    h = re.sub(r"[^\w\- ]", "", h)
    return h.replace(" ", "-")


def markdown_files(src: bool) -> list[Path]:
    files = sorted((ROOT / "docs").rglob("*.md"))
    files += [p for p in (ROOT / "ORIENTATION.md", ROOT / "README.md",
                          ROOT / "src" / "README.md") if p.exists()]
    if src:
        files += sorted((ROOT / "src").glob("*.py"))
    return files


def _defining(line: str, pat) -> bool:
    """A use that explains the term: it is bold or italic there, or a defining
    phrase ('which X', 'we call', 'denotes', 'means', ...) sits within 60
    characters of it.  Abstract blockquotes are skipped: they summarise."""
    if line.lstrip().startswith(">"):
        return False
    for m in pat.finditer(line):
        a, b = m.span()
        if re.search(r"\*\*[^*]*$", line[:a]) and re.search(r"^[^*]*\*\*", line[b:]):
            return True
        if re.search(r"(?<!\*)\*[^*]*$", line[:a]) and re.search(r"^[^*]*\*(?!\*)", line[b:]):
            return True
        if DEFINING.search(line[max(0, a - 60):b + 60]):
            return True
    return False


def scan(term: str, files: list[Path]):
    pat = re.compile(r"(?<![\w-])" + re.escape(term).replace(r"\ ", r"[\s-]+")
                     + r"\w*", re.I)
    out = collections.OrderedDict()
    for path in files:
        heading, anchor, in_code = "", "", False
        lines = path.read_text(errors="replace").splitlines()
        for i, line in enumerate(lines, 1):
            if path.suffix == ".md":
                if line.startswith("```"):
                    in_code = not in_code
                m = None if in_code else re.match(r"^(#{1,6})\s+(.*)", line)
                if m:
                    heading, anchor = m.group(2), slug(m.group(2))
            if pat.search(line):
                is_head = path.suffix == ".md" and line.startswith("#")
                kind = "*" if (is_head or _defining(line, pat)) else " "
                out.setdefault(path, []).append(
                    (i, heading, anchor, kind, line.strip()))
    return out


def report(term: str, hits, n_lines: int) -> None:
    n_files = len(hits)
    total = sum(len(v) for v in hits.values())
    print(f"{term!r}: {total} uses in {n_files} file(s)\n")

    def key(item):
        path, rows = item
        return (-sum(r[3] == "*" for r in rows), -len(rows), str(path))
    for path, rows in sorted(hits.items(), key=key):
        rel = path.relative_to(ROOT)
        print(f"{rel}  ({len(rows)})")
        shown = sorted(rows, key=lambda r: (r[3] != "*", r[0]))[:n_lines]
        for ln, head, anchor, kind, text in sorted(shown):
            where = f"{rel}#{anchor}" if anchor else str(rel)
            print(f"  {kind} L{ln:<5} {where}")
            print(f"           {text[:150]}")
        print()


def candidates(min_files: int, top: int) -> int:
    """Words the notes use far more than English does, in >= min_files files.

    Keyness = (share of the notes' words) / (share of English words), so it
    finds both jargon ('eikonal', 'positon') and ordinary words put to
    special use ('sea', 'floor', 'reach', 'pair').  Needs ``wordfreq``
    (pip install wordfreq); a development aid, not a dependency of the model.
    """
    try:
        from wordfreq import zipf_frequency
    except ImportError:
        print("--candidates needs the wordfreq package: pip install wordfreq\n"
              "(in a virtual environment if pip says externally-managed-environment; see the module docstring)")
        return 2
    files = [p for p in markdown_files(False) if p.name != "INDEX.md"
             and "GLOSSARY" not in p.name]
    df, tf = collections.Counter(), collections.Counter()
    cap = collections.defaultdict(lambda: [0, 0])      # [lower, Capitalised]
    for path in files:
        text = path.read_text(errors="replace")
        text = re.sub(r"```.*?```", " ", text, flags=re.S)
        text = re.sub(r"\$`.*?`\$|\$[^$\n]*\$|`[^`\n]*`|\]\([^)]*\)|https?://\S+",
                      " ", text)
        text = re.sub(r"\\[A-Za-z]+", " ", text)          # stray LaTeX commands
        raw = re.findall(r"[A-Za-z][A-Za-z\-]{3,}", text)
        for w in raw:
            cap[w.lower()][w[0].isupper()] += 1
        words = [w.lower() for w in raw]
        df.update(set(words))
        tf.update(words)
    n_tok = sum(tf.values())
    gloss = ROOT / "docs" / "GLOSSARY.md"
    have = gloss.read_text().lower() if gloss.exists() else ""
    rows = []
    for w, n in df.items():
        if n < min_files:
            continue
        lo, up = cap[w]
        if up > 3 * lo:                    # a proper noun: bibliography's job
            continue
        eng = max(10 ** (zipf_frequency(w, "en") - 9), 1e-9)
        rows.append((tf[w] / n_tok / eng, w, n, tf[w],
                     zipf_frequency(w, "en")))
    rows.sort(reverse=True)
    print(f"{len(rows)} words in at least {min_files} of {len(files)} files; "
          f"top {top} by keyness (* = already a heading in GLOSSARY.md)\n")
    print(f"  {'word':24s} {'files':>5} {'uses':>6} {'zipf':>5}  keyness")
    for k, w, n, t, z in rows[:top]:
        mark = "*" if re.search(r"(?m)^#{3,4}\s+.*\b" + re.escape(w), have) else " "
        print(f"{mark} {w:24s} {n:5d} {t:6d} {z:5.1f}  {k:9.0f}")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("term", nargs="?")
    ap.add_argument("--lines", type=int, default=2,
                    help="lines of context per file (default 2)")
    ap.add_argument("--src", action="store_true", help="also search src/*.py")
    ap.add_argument("--candidates", type=int, metavar="N",
                    help="list words used in at least N files instead")
    ap.add_argument("--top", type=int, default=150,
                    help="length of the --candidates list (default 150)")
    a = ap.parse_args(argv)
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)      # quiet under | head
    if a.candidates:
        return candidates(a.candidates, a.top)
    if not a.term:
        ap.error("give a term, or --candidates N")
    files = [p for p in markdown_files(a.src) if p.name != "GLOSSARY.md"]
    gloss = ROOT / "docs" / "GLOSSARY.md"
    if gloss.exists():
        for line in gloss.read_text().splitlines():
            if re.match(r"^#{3,4}\s", line) and re.search(
                    r"(?<![\w-])" + re.escape(a.term), line, re.I):
                print(f"Glossary entry: docs/GLOSSARY.md#{slug(line.lstrip('# '))}\n")
                break
    hits = scan(a.term, files)
    if not hits:
        print(f"{a.term!r}: not found")
        return 1
    report(a.term, hits, a.lines)
    return 0


if __name__ == "__main__":
    sys.exit(main())
