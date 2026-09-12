# WPMW

Notes, derivations, and Python implementations for a phase-space-crystal-lattice
interpretation of the Wigner equation, plus the underlying extended Fokker–Planck
(xFP) machinery that motivates it.

This repository is a working set of research notes and code. It is not a polished
library and the algorithms here are still being explored. Conventions may change.

## Getting started

For most collaborators who want to read the documentation or run scripts,
clone only the `main` branch:

```bash
git clone --single-branch --branch main https://github.com/billpage/wpmw.git
cd wpmw
```

This skips the `output` branch, which stores demo figures as binary blobs and
can grow over time. The figures are still viewable online — they are embedded
in the documentation via raw GitHub URLs and render directly in the GitHub
markdown viewer without cloning anything.

If you also need to commit new figures to the `output` branch (see
[Sharing figures](#sharing-figures-via-the-output-branch) below), add the
worktree after cloning:

```bash
git fetch origin output
git worktree add ../wpmw-output output
```

## Repository layout

```
docs/
  algorithm/    Core algorithm specifications
  supplement/   Extended analysis and discussion
  analysis/     Mathematical derivations and review notes
  notebooks/    Notebook companions to supplement tutorials (outputs stripped;
                rendered copies live on the `output` branch)
src/            Python implementations (runnable demos / scripts)
  wpmwlib/      Shared library modules imported by the scripts
references/     bibliography.md (links to papers; PDFs are NOT committed)
```

See [`docs/README.md`](docs/README.md) for a guide to what's in each
documentation subdirectory and where to start reading.

## Output path convention

All Python scripts in WPMW must write their outputs (PNGs, MP4s, CSVs, etc.)
through helpers in `src/wpmwlib/wpmw_utils.py`. Do not hardcode paths such as
`/home/claude/`, `/mnt/user-data/outputs/`, or `/kaggle/working/`.

**`output_path(filename)`** — runtime scratch output (Claude container, Kaggle,
local runs):

```python
from wpmwlib.wpmw_utils import output_path

fig.savefig(output_path("my_figure.png"), dpi=150, bbox_inches="tight")
```

**`docs_path(filename)`** — figures destined for the `output` branch so they can
be embedded in documentation and shared with collaborators (see below):

```python
from wpmwlib.wpmw_utils import output_path, docs_path

fig.savefig(output_path("my_figure.png"), dpi=150, bbox_inches="tight")
dp = docs_path("my_figure.png")
if dp:
    fig.savefig(dp, dpi=150, bbox_inches="tight")
```

Do not use `shutil.copy2` to mirror outputs to a second location — write
directly to each helper.

Environment variable summary:

| Context | `WPMW_OUTPUT` | `WPMW_DOCS` |
|---|---|---|
| Claude container | `/mnt/user-data/outputs` | (unset) |
| Kaggle | `/kaggle/working` | (unset) |
| Local (post-patch) | (unset → `./output`) | `~/wpmw-output` (worktree) |

## Sharing figures via the `output` branch

Demo figures are stored on an orphan `output` branch — separate from the source
history — and embedded in documentation via stable raw GitHub URLs. This keeps
`main` free of binary blobs while making figures accessible to anyone with a link.

**One-time worktree setup** (after the `output` branch exists on origin):

```bash
cd ~/wpmw
git worktree add ../wpmw-output output
```

**Running a script locally to produce keeper figures:**

```bash
export WPMW_OUTPUT=./output
export WPMW_DOCS=~/wpmw-output
python src/my_script.py
```

**Committing the figures from the worktree:**

```bash
cd ~/wpmw-output
git add figures/my_figure.png
git commit -m "add my_figure output from my_script.py"
git push origin output
```

**Embedding a figure in a markdown doc on `main`:**

```markdown
![My figure](https://raw.githubusercontent.com/billpage/wpmw/output/figures/my_figure.png)
```

The URL pattern is always:
`https://raw.githubusercontent.com/billpage/wpmw/output/<subdir>/<filename>`

## Notebooks

Some ladder steps also have a Jupyter notebook companion that interleaves a
supplement tutorial's prose with the same code from its demo script, plus
its generated figures — meant for interactive reading, not as a replacement
for the tutorial or the script, which remain canonical.

Notebooks follow the same two-copy split as figures, using the same
`docs_path()` helper with an explicit `subdir`:

```python
from wpmwlib.wpmw_utils import docs_path

dp = docs_path("my_notebook.ipynb", subdir="notebooks")
```

- **Rendered copy** (outputs and figures intact) — committed to the `output`
  branch worktree at `notebooks/<filename>.ipynb`, exactly parallel to
  `figures/`. This is the copy worth reading or running interactively.
- **Stripped copy** (`jupyter nbconvert --clear-output`, no outputs, no
  execution counts) — committed to `docs/notebooks/<filename>.ipynb` on
  `main`, so future patches to it stay reviewable as plain-text diffs
  instead of noisy JSON-with-embedded-images diffs.

**Linking to both from a markdown doc on `main`:**

```markdown
[Notebook (rendered, output branch)](https://github.com/billpage/wpmw/blob/output/notebooks/my_notebook.ipynb) ·
[Notebook (source, main)](https://github.com/billpage/wpmw/blob/main/docs/notebooks/my_notebook.ipynb)
```

Use the `/blob/` URL, not `/raw/` — GitHub renders `.ipynb` files (markdown,
code, and any saved outputs) at the `/blob/` path; `/raw/` would just serve
the JSON.

## File conventions

- Follow the existing directory structure for new files.
- Drop version suffixes from filenames — git tracks revisions.
- Do not commit copyrighted PDFs, personal files, or admin material.
- Reference papers go in `references/bibliography.md` as links, not as PDFs.

## License

This repository uses two licenses, split by content type:

- **Code** (`src/`, including `wpmwlib/`) — [GNU Affero General Public
  License v3.0 or later](LICENSE) (`AGPL-3.0-or-later`). Copyleft that
  extends to network use: a modified version served over a network must
  still offer its source to users interacting with it, not only to those
  who receive a distributed copy.
- **Documentation, derivations, and figures** (`docs/`, `references/`,
  and the images on the `output` branch) — [Creative Commons
  Attribution-ShareAlike 4.0 International](LICENSE-DOCS)
  (`CC-BY-SA-4.0`). Adaptations must credit the original authors and stay
  under the same license, whether or not they are formally published.

See [`LICENSE`](LICENSE) and [`LICENSE-DOCS`](LICENSE-DOCS) for the full
terms. If you build on this work, please cite it as described below.

## Citation

The primary content of this repository is the research notes, derivations,
and analysis under `docs/` — the Python code in `src/` is a secondary,
numerical-verification companion to those notes. Cite the notes as below;
a machine-readable version is in [`CITATION.cff`](CITATION.cff) as the
`preferred-citation`, which GitHub's "Cite this repository" option and
tools such as Zenodo read in preference to the plain software entry.

Plain text:

> Bill Page (2026). *WPMW: notes and derivations on a compensated
> Liouville algorithm for quantum dynamics via signed particles*
> [Research notes]. https://github.com/billpage/wpmw/tree/main/docs

BibTeX:

```bibtex
@techreport{wpmw2026notes,
  author = {Page, Bill},
  title  = {{WPMW}: notes and derivations on a compensated {Liouville}
             algorithm for quantum dynamics via signed particles},
  year   = {2026},
  url    = {https://github.com/billpage/wpmw/tree/main/docs},
  note   = {Accessed: insert date}
}
```

If you specifically need to cite the code rather than the notes (e.g. you
ran the simulations but aren't drawing on the derivations), use the
top-level entry in `CITATION.cff` instead — same author, but scoped to
`src/` and the AGPL-3.0-or-later license.
