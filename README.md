# WPMW

Notes, derivations, and Python implementations for a phase-space-crystal-lattice
interpretation of the Wigner equation, plus the underlying extended Fokker–Planck
(xFP) machinery that motivates it.

This repository is a working set of research notes and code. It is not a polished
library and the algorithms here are still being explored. Conventions may change.

## Repository layout

```
docs/
  algorithm/    Core algorithm specifications
  supplement/   Extended analysis and discussion
  analysis/     Mathematical derivations and review notes
src/            Python implementations
references/     bibliography.md (links to papers; PDFs are NOT committed)
```

## Output path convention

All Python scripts in WPMW must write their outputs (PNGs, MP4s, CSVs, etc.)
through the `output_path()` helper in `src/wpmw_utils.py`. Do not hardcode paths
such as `/home/claude/`, `/mnt/user-data/outputs/`, or `/kaggle/working/`.

```python
from wpmw_utils import output_path

fig.savefig(output_path("my_figure.png"), dpi=150, bbox_inches="tight")
```

The output directory is controlled by the `WPMW_OUTPUT` environment variable.
If unset, files go to `./output`. Do not use `shutil.copy2` to mirror outputs to
a second location — write directly to `output_path(...)`.

When running scripts in the Claude container, set
`WPMW_OUTPUT=/mnt/user-data/outputs` before execution. The Kaggle runner notebook
sets `WPMW_OUTPUT=/kaggle/working` automatically.

## Contributing changes

Bill Page is the gatekeeper for all commits to this repository. The standard flow
for collaborator-proposed changes is:

1. Clone the repository and make changes on a local branch.
2. Generate a patch with `git format-patch` and share the resulting `.patch`
   file(s).
3. Bill applies the patch locally and pushes:
   ```bash
   cd ~/wpmw
   git am ~/Downloads/NNNN-<patch-name>.patch
   git push origin main
   ```

Notes on file conventions:

- Follow the existing directory structure for new files.
- Drop version suffixes from filenames — git tracks revisions.
- Do not commit copyrighted PDFs, personal files, or admin material.
- Reference papers go in `references/bibliography.md` as links, not as PDFs.

## License

TODO — license not yet selected. Until a license is added, default copyright
applies and the contents should be treated as all-rights-reserved by the author.
