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
- `demo_qho_ground_state.py` — demonstration that exercises both solvers on
  the QHO ground state and writes a comparison figure via `output_path()`.

## Sign-convention finding (action item)

The spec's verbatim Fourier-mode update (§3c) implements
`dW/dt = -V'(x) dW/dp` in the continuum limit, which is the *opposite* sign
from the spec's own stated QLE (`dW/dt = +V'(x) dW/dp`). Empirically this
pushes a coherent state *toward* a potential hill instead of away from it.
For rotationally symmetric states (e.g. the QHO ground state) the sign is
irrelevant, but for general states it produces backward dynamics.

`PhaseSpaceCrystalLattice.step_jump_fourier` exposes a `qle_sign` parameter
to flip the sign back to the QLE convention. The differential form
(`step_jump_differential`) always uses the QLE-consistent sign. The spec
itself should be reviewed and corrected (or its sign convention for V or for
the QLE more clearly justified).

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
