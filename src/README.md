# src

Python implementations for the WPMW project.

## Modules

- `wpmw_utils.py` — shared helpers (currently the `output_path()` router).
- `phase_space_crystal_lattice.py` — implementation of the algorithm in
  `docs/algorithm/phase_space_crystal_lattice_algorithm.md`. Exposes both
  the deterministic mesh-density form (spec §3c) and the Monte-Carlo particle
  form (spec §6), and both the Fourier-mode (§3b) and differential (§7b) jump
  forms.
- `wigner_split_fourier.py` — reference solver: Strang-split spectral Fourier
  on the Wigner equation. Specialized to QHO; for general V the force-kick
  kernel must be replaced by the full Wigner–Moyal kernel.
- `demo_qho_ground_state.py` — demo: QHO ground-state preservation; compares
  the crystal-lattice solver and the split-Fourier reference. Insensitive to
  the force-term sign because the ground state is rotationally symmetric.
- `demo_coherent_state.py` — demo: cosine-potential dynamics for a coherent
  Gaussian placed on the downhill slope. Sensitive to the sign convention;
  reproduces the textbook Newtonian centroid trajectory.
- `sign_convention_check.py` — regression test for the §6.3 sign correction
  in `docs/supplement/phase_space_crystal_lattice_supplement.md`. Compares
  three candidate discrete update rules (V2 general formula, V2 simplified /
  Python, and original spec §3c) on a coherent state, confirming that only
  the QLE-consistent form drives the centroid downhill.

## Output path convention

All scripts in this directory must write files through `output_path()` from
`wpmw_utils`, never via hardcoded absolute paths:

```python
from wpmw_utils import output_path

fig.savefig(output_path("my_figure.png"), dpi=150, bbox_inches="tight")
```

The output directory is controlled by the `WPMW_OUTPUT` environment variable.
If unset, files go to `./output`.

See the top-level `README.md` for the full convention.
