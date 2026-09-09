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

If you use or build on this repository, please credit the authors. A
machine-readable citation is also available in
[`CITATION.cff`](CITATION.cff), which GitHub uses to generate a "Cite
this repository" option.

Plain text:

> Bill Page and contributors (2026). *WPMW: A compensated Liouville
> algorithm for quantum dynamics via signed particles* [Software and
> research notes]. GitHub. https://github.com/billpage/wpmw

BibTeX:

```bibtex
@misc{wpmw2026,
  author       = {Page, Bill and contributors},
  title        = {{WPMW}: A compensated {Liouville} algorithm for quantum
                   dynamics via signed particles},
  year         = {2026},
  howpublished = {GitHub repository},
  url          = {https://github.com/billpage/wpmw},
  note         = {Accessed: insert date}
}
```
