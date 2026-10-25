# docs

Project documentation, split into four directories that play different
roles. See the top-level [`README.md`](../README.md) for repository-wide
conventions (cloning, output paths, the `output` branch), and
[`ORIENTATION.md`](../ORIENTATION.md) for the conceptual guide — what the
central idea is, which words mean what, and what has been retracted. This
page is a guide to what's in each subdirectory here and where to start
reading.

## The four directories

- **[`algorithm/`](algorithm)** — The canonical specifications.
  Precise enough to re-implement in another language, and the only
  directory you need if you just want to run or reproduce what's in `src/`.
  Start with `phase_space_crystal_lattice_algorithm.md`; it alone reproduces
  every result in `src/`.

- **[`analysis/`](analysis)** — The mathematical case for why the
  algorithm is what it is. A **derivation ladder** of notes, each taking as
  input something its predecessor postulated and ending with the open items
  that motivate the next. Read in order starting from
  `phase_space_crystal_lattice_review.md` if you want to follow the
  reasoning rather than just the result. To look up a theorem or open item by
  its label (Theorem K4, S-SP3, C2), use [`analysis/INDEX.md`](analysis/INDEX.md).

- **[`supplement/`](supplement)** — Background that supports, but
  isn't part of, the algorithm specs: redrafted source memos from David
  Cyganski, the interaction-diagram drawings, targeted test cases (e.g. the
  inverted pair barrier), a close reading of Takabayasi (1954), and the
  `figures/` cited by those documents.

- **[`notebooks/`](notebooks)** — Jupyter notebook companions to selected
  supplement tutorials, interleaving that tutorial's prose with the same
  code from its demo script and the figures it generates, for interactive
  reading. Outputs are stripped here for reviewable diffs; a rendered copy
  of each notebook, with outputs and figures intact, lives on the `output`
  branch (see the top-level README's "Notebooks" section). These are a
  companion rendering, not a replacement for the tutorial or the script,
  which remain canonical.

- **[`interactive/`](interactive)** — self-contained interactive web pages
  (one HTML file each, no external scripts) that let a reader vary the
  parameters of a figure in an analysis note. They are served by GitHub
  Pages (see the top-level README's "Interactive pages" section). Like the
  notebooks they add no claims: each page names the note and the demo whose
  computation it repeats. Slide decks prepared for calls with collaborators
  live here too, as one HTML file each with their speaker notes; they quote
  the notes and demos they summarise and add nothing to them. The first is
  [`worlds_and_particles.html`](interactive/worlds_and_particles.html), for
  the call of 27 October 2026.

## Reading aids

- **[`GLOSSARY.md`](GLOSSARY.md)** — the terms the notes in `algorithm/`,
  `analysis/` and `supplement/` use, with what each means *here*, where it is
  defined, and what not to confuse it with. For where a word is used, run
  `PYTHONPATH=src python3 -m wpmwlib.where_used WORD`.

## How they relate

`algorithm/` is self-contained — read it alone to implement or run the
model. `analysis/` and `supplement/` are where the algorithm's claims get
justified, corrected, and cross-checked against both the source material and
each other; `analysis/` carries the theorems and the ladder structure,
`supplement/` carries the background documents and worked test cases those
theorems draw on. `notebooks/` contains Jupyter notesbooks with illustrative
code and tutorial level explanations. `interactive/` contains interactive
HTML-based graphics and examples that are accessible through rendered github
pages. Each subdirectory's own `README.md` links every document it contains.
