# docs

Project documentation, split into four directories that play different
roles. See the top-level [`README.md`](../README.md) for repository-wide
conventions (cloning, output paths, the `output` branch); this page is a
guide to what's in each subdirectory here and where to start reading.

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
  reasoning rather than just the result.

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

## How they relate

`algorithm/` is self-contained — read it alone to implement or run the
model. `analysis/` and `supplement/` are where the algorithm's claims get
justified, corrected, and cross-checked against both the source material and
each other; `analysis/` carries the theorems and the ladder structure,
`supplement/` carries the background documents and worked test cases those
theorems draw on. `notebooks/` sits on top of `supplement/`, one level
further towards interactive reading rather than new claims. Each
subdirectory's own `README.md` links every document it contains.
