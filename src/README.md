# src

Python implementations for the WPMW project.

## Modules

- `wpmw_utils.py` — shared helpers: `output_path()` (runtime scratch output)
  and `docs_path()` (figures for the `output` branch, shareable via URL).
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
- `demo_cat_state.py` — demo: free-particle (V = 0) evolution of a
  Schrödinger-cat state (two colliding Gaussian wave packets in
  superposition).  Compares the PSC solver to the closed-form Wigner
  function at four times and overlays sample classical trajectories from
  interference-node seeds, illustrating that the cat state's nodes are
  rigidly transported along free-particle characteristics.
- `demo_cat_state_microdynamics.py` — companion demo: the same cat-state
  problem, but as a Monte-Carlo crystal-lattice microdynamics
  simulation.  Samples ~2×10⁷ positons from the shifted distribution
  W' = W + 2/h (everywhere non-negative, since |W| ≤ 2/h for any pure
  state), streams each one ballistically (V=0 means no mediated jumps),
  and reconstructs W from the binned counts as ρ_emp − 2/h.
  Demonstrates that the Wigner function — including its negative
  interference fringes — is faithfully recovered from a strictly
  non-negative ensemble, and that individual positon trajectories are
  trivial classical horizontal lines regardless of whether their initial
  conditions happen to coincide with Wigner-function nodes.

  Sample output figures (committed on the `output` branch):

  ![Cat-state microdynamics evolution](https://raw.githubusercontent.com/billpage/wpmw/output/figures/cat_state_microdynamics_evolution.png)

  3×4 grid at default parameters (sampling grid 256², reconstruction grid
  128², N = 2×10⁷): MC reconstruction (top), exact closed-form Wigner
  (middle), and pointwise difference (bottom) at t = 0, t_c/2, t_c, 3 t_c/2.
  L² deviation between MC and exact stays at ~0.38 across all times (peak
  |W| ≈ 0.32), dominated by Poisson sampling of the 2/h background.

  ![Cat-state microdynamics marginals](https://raw.githubusercontent.com/billpage/wpmw/output/figures/cat_state_microdynamics_marginals.png)

  Position-space probability density |ψ(x, t)|² at the same times — the
  textbook double-slit-like fringe pattern at t = t_c is recovered cleanly
  from the all-positon ensemble.

  ![Cat-state microdynamics trajectories](https://raw.githubusercontent.com/billpage/wpmw/output/figures/cat_state_microdynamics_trajectories.png)

  Trajectory portrait of the six tagged test positons.  Three with t = 0
  initial coordinates near interference-node lines (red), three at non-node
  locations including the lobe centres (blue).  Every trajectory is a
  horizontal line at constant p; "near-node" status leaves no microdynamic
  signature.
- `sign_convention_check.py` — regression test for the §6.3 sign correction
  in `docs/supplement/phase_space_crystal_lattice_supplement.md`. Compares
  three candidate discrete update rules (V2 general formula, V2 simplified /
  Python, and original spec §3c) on a coherent state, confirming that only
  the QLE-consistent form drives the centroid downhill.

## Output path convention

All scripts in this directory must write files through helpers from
`wpmw_utils`, never via hardcoded absolute paths.

Use `output_path()` for all runtime scratch output:

```python
from wpmw_utils import output_path

fig.savefig(output_path("my_figure.png"), dpi=150, bbox_inches="tight")
```

Use `docs_path()` additionally for figures that should be committed to the
`output` branch and embedded in documentation:

```python
from wpmw_utils import output_path, docs_path

fig.savefig(output_path("my_figure.png"), dpi=150, bbox_inches="tight")
dp = docs_path("my_figure.png")
if dp:
    fig.savefig(dp, dpi=150, bbox_inches="tight")
```

`docs_path()` returns `None` when `WPMW_DOCS` is unset (cloud environments),
so the `if dp:` guard is always required.

See the top-level `README.md` for the full convention including the `output`
branch worktree setup and figure embedding instructions.
