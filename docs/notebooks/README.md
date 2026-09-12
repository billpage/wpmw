# docs/notebooks

Jupyter notebook companions to selected `docs/supplement/` tutorials. Each
notebook interleaves that tutorial's prose with the same code from its
companion demo script in `src/`, reorganised so the code runs where the
tutorial first discusses it, plus the figures the script generates. These
are a companion rendering for interactive reading, not a replacement — where
a notebook and its source tutorial or script disagree, the tutorial and the
script win.

Outputs are stripped from the copies here so patches stay reviewable as
plain-text diffs. A rendered copy of each notebook, with outputs and figures
intact, lives on the `output` branch at `notebooks/<filename>.ipynb`; see
the top-level README's "Notebooks" section for the convention and how to
link to both copies from a markdown doc.

## Contents

- **[`emission_and_absorption.ipynb`](emission_and_absorption.ipynb)** —
  Companion to
  [`../supplement/emission_and_absorption.md`](../supplement/emission_and_absorption.md)
  and [`../../src/demo_emission_and_absorption.py`](../../src/demo_emission_and_absorption.py).
  Rendered copy (outputs and figures intact):
  [`output` branch](https://github.com/billpage/wpmw/blob/output/notebooks/emission_and_absorption.ipynb).
