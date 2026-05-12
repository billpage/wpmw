# src

Python implementations for the WPMW project.

The runnable demos and regression tests live directly in this directory.
Shared library code lives in the `wpmwlib/` subpackage (see below) and is
imported by the scripts here.

## Library modules (`wpmwlib/`)

- `wpmwlib/wpmw_utils.py` — shared helpers: `output_path()` (runtime scratch
  output) and `docs_path()` (figures for the `output` branch, shareable via
  URL).
- `wpmwlib/phase_space_crystal_lattice.py` — implementation of the algorithm
  in `docs/algorithm/phase_space_crystal_lattice_algorithm.md`. Exposes both
  the deterministic mesh-density form (spec §3c) and the Monte-Carlo particle
  form (spec §6), and both the Fourier-mode (§3b) and differential (§7b) jump
  forms.
- `wpmwlib/wigner_split_fourier.py` — reference solver: Strang-split spectral
  Fourier on the Wigner equation. Specialized to QHO; for general V the
  force-kick kernel must be replaced by the full Wigner–Moyal kernel.
- `wpmwlib/check_md_math.py` — markdown LaTeX-math linter (see "Markdown
  math linter" below).

## Runnable scripts

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
  rigidly transported along free-particle characteristics.  Also exposes
  shared constants and helpers (`HBAR`, `MASS`, `X0`, `P0`, `SIGMA`, `T_C`,
  `L`, `W_cat_initial`, `W_cat_exact`, `sample_node_seeds`) that
  `demo_cat_state_microdynamics.py` reuses via a sibling import.
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
- `demo_cosine_well_microdynamics.py` — companion to the cat-state
  microdynamics demo for a single-period **cosine well**
  V(x) = -V_p cos(2πx/L) (minimum at x = 0).  Where the cat-state demo
  runs at V = 0 (positon evolution is pure ballistic streaming and *is*
  the full QLE), this demo compares two distinct evolutions on the same
  grid:  (i) the full QLE on the PSC mesh via
  `PhaseSpaceCrystalLattice.strang_step_fourier`, exact for a single-mode
  cosine including all higher-derivative Moyal terms that produce
  Wigner-function negative regions; and (ii) classical-positon Monte
  Carlo — 5×10⁶ positons sampled from W' = W + 2/h at t = 0, evolved
  under Hamilton's equations alone (no quantum jumps), binned, and
  reconstructed as ρ_emp − 2/h.  This is exactly what the QLE would give
  in the ℏ → 0 limit (the leading Liouville term).  Initial state is a
  min-uncertainty squeezed-vacuum Gaussian at (0, 0) with
  σ_x = 2 σ_{x,gs} (energy ≈ 1.06 ℏω in the harmonic approximation),
  chosen so the well's quartic anharmonicity seeds visible Wigner
  negativity within a few classical periods.  The MC ensemble is
  integrated with a float32 in-place velocity-Verlet stepper that runs
  ≈ 4× faster than the fp64 version at this size; cumulative drift over
  800 steps is far below the binning resolution.

  Sample output figures (committed on the `output` branch):

  ![Cosine-well microdynamics evolution](https://raw.githubusercontent.com/billpage/wpmw/output/figures/cosine_well_microdynamics_evolution.png)

  3×4 grid at default parameters (mesh 128², reconstruction grid 48²,
  N = 5×10⁶, 800 steps over 4 T_period): full QLE on the mesh (top),
  classical-positon MC (middle), and pointwise difference MC − QLE
  (bottom) at t = 0, T_period, 2 T_period, 4 T_period.  The QLE row
  develops visible negative regions on the wavepacket flanks; the
  classical row stays strictly non-negative; the difference row evolves
  from pure shot noise at t = 0 to a clearly structured pattern by
  t = 4 T_period that marks where the QLE has put negative weight that
  the classical evolution misses.  Six tagged classical Hamilton orbits
  are overlaid on rows 1 and 2.

  ![Cosine-well microdynamics negativity](https://raw.githubusercontent.com/billpage/wpmw/output/figures/cosine_well_microdynamics_negativity.png)

  Wigner negativity ∫|min(W, 0)| dx dp over the run.  QLE accumulates
  real negativity from 0 to ≈ 0.22 over four periods.  Classical MC sits
  at a flat shot-noise floor of ≈ 1.5 (independent of t) — that floor is
  fundamental for finite-N empirical reconstruction of W, not evidence
  of physical negativity in the classical evolution.

  ![Cosine-well microdynamics marginals](https://raw.githubusercontent.com/billpage/wpmw/output/figures/cosine_well_microdynamics_marginals.png)

  Position-space probability density ρ(x) = ∫ W dp at the same four
  times.  Integrating over momentum collapses the negative regions, so
  QLE and MC agree well in the marginal even where they differ pointwise
  in W.

  ![Cosine-well microdynamics trajectories](https://raw.githubusercontent.com/billpage/wpmw/output/figures/cosine_well_microdynamics_trajectories.png)

  Phase-space portrait of the six classical Hamilton orbits used as
  overlays in the evolution figure.  Three "near-bottom" orbits inside
  ~σ_{x,gs} of the well minimum (red) are nearly elliptical — the
  harmonic regime.  Three "wider" orbits at amplitude ~ 2 σ_{x,gs}
  (blue) show faint deformation from cosine anharmonicity but remain
  well bound (turning points well below V_max = +V_p).  Tick marks at
  t = T_period, 2 T_period, 3 T_period sit close together on each orbit,
  confirming the period.
- `sign_convention_check.py` — regression test for the §6.3 sign correction
  in `docs/supplement/phase_space_crystal_lattice_supplement.md`. Compares
  three candidate discrete update rules (V2 general formula, V2 simplified /
  Python, and original spec §3c) on a coherent state, confirming that only
  the QLE-consistent form drives the centroid downhill.

## Output path convention

All scripts in this directory must write files through helpers from
`wpmwlib.wpmw_utils`, never via hardcoded absolute paths.

Use `output_path()` for all runtime scratch output:

```python
from wpmwlib.wpmw_utils import output_path

fig.savefig(output_path("my_figure.png"), dpi=150, bbox_inches="tight")
```

Use `docs_path()` additionally for figures that should be committed to the
`output` branch and embedded in documentation:

```python
from wpmwlib.wpmw_utils import output_path, docs_path

fig.savefig(output_path("my_figure.png"), dpi=150, bbox_inches="tight")
dp = docs_path("my_figure.png")
if dp:
    fig.savefig(dp, dpi=150, bbox_inches="tight")
```

`docs_path()` returns `None` when `WPMW_DOCS` is unset (cloud environments),
so the `if dp:` guard is always required.

See the top-level `README.md` for the full convention including the `output`
branch worktree setup and figure embedding instructions.

## Markdown math linter

`wpmwlib/check_md_math.py` lints the project's markdown files for math that
will not render correctly on GitHub. It catches the two failure modes we have
actually hit:

- macros that vanilla LaTeX accepts but GitHub's MathJax config blocks —
  notably `\operatorname`, `\bm`, `\href`, `\DeclareMathOperator`,
  `\newcommand`, `\definecolor`, `\colorbox`, `\label` / `\ref` / `\eqref`,
  `\tag`, `\intertext`, `\verb`, `\mathds`;
- multi-line `$$...$$` display blocks placed inside a list item, which
  GitHub silently re-parses as nested bullets.

Optionally it also feeds every expression to KaTeX (strict mode) and to
MathJax 3 to catch malformed LaTeX (mismatched delimiters, unknown macros,
etc.). Those passes need `node` plus `katex` and `mathjax-full` from npm; if
they are not available, the linter prints a notice and skips them, still
running the static and structural passes.

### Run locally

From the repository root:

```bash
# fast, no Node needed — static + structural passes only:
python -m wpmwlib.check_md_math --no-render

# full check — first install the npm packages once:
npm install --no-save katex mathjax-full
python -m wpmwlib.check_md_math
```

The default scan targets are `docs/` and `README.md`. Pass any file or
directory to override. Exit code is `0` on a clean scan, `1` if issues
were reported.

The same check runs in CI on every push and pull request via
`.github/workflows/check_md_math.yml`.

### Style guide for math in WPMW markdown

A short cheat sheet for keeping new docs lint-clean:

- **Avoid backslash-escaped TeX shortcuts inside math.** GitHub's markdown
  preprocessor strips the leading backslash from any ``\X`` where X is ASCII
  punctuation, *even inside math blocks*, before MathJax sees the content.
  This silently corrupts spacing and turns ``\bigl\{...\bigr\}`` into the
  hard "Missing or unrecognized delimiter for \\bigl" error. Replace with
  one of the safe forms below:

  | Don't write | Write instead |
  | --- | --- |
  | `\,` (thin space) | `\thinspace` (preferred) or `\\,` |
  | `\!` (negative thin space) | `\negthinspace` (preferred) or `\\!` |
  | `\;` (thick space) | `\\;` (no working letter-named form) |
  | `\:` (medium space) | `\\:` (no working letter-named form) |
  | `\{` (literal left brace) | `\lbrace` (preferred) or `\\{` |
  | `\}` (literal right brace) | `\rbrace` (preferred) or `\\}` |

  Note: `\thickspace` and `\medspace` look like the natural letter-named
  alternatives for `\;` and `\:`, but they are *not* defined in MathJax 3
  with only `base` and `ams` packages — GitHub's actual config — so on
  GitHub they render as raw text instead of math spacing. Use doubled
  backslash (`\\;` / `\\:`) for thick and medium spaces.

  The linter's GFM pass enforces this rule.

- For upright function names (erf, erfc, sgn, Tr, ...) use `\mathrm{...}`,
  not `\operatorname{...}` — same glyph, math-mode spacing, universally
  supported.
- For prose-like content inside math (subscripts like `_{\text{short-range}}`,
  unit labels, etc.) use `\text{...}`.
- Keep `$$...$$` display blocks on a **single source line** when they appear
  inside a numbered or bulleted list item. If you need visual line breaks,
  either use `$$\begin{aligned} ... \\ ... \end{aligned}$$` on one line, or
  switch the block to a ```` ```math ```` fenced code block (see below) —
  fenced code blocks are recognised inside list items even when split over
  multiple lines.
- Outside of list items, multi-line `$$...$$` is fine — preferred for long
  derivations.
- **Alternative display syntax: ```` ```math ```` fenced blocks.**
  GitHub also accepts a fenced-code form for display math:

  ````
  ```math
  \frac{\partial W}{\partial t} + \frac{p}{m}\frac{\partial W}{\partial x} = 0
  ```
  ````

  This is equivalent to `$$...$$` for display math, but is **exempt from
  the CommonMark backslash-strip pipeline** (verified empirically: see
  the experiment branch). Inside a fenced math block you can write
  `\,`, `\;`, `\!`, `\bigl\{`, `\bigr\}` directly — the backslashes
  reach MathJax intact.

  Two reasons to prefer the fenced form:

  1. *Heavy use of backslash-escapes* — equations with lots of TeX
     spacing or sized-delimiter braces are clearer with `\,` `\;`
     `\bigl\{` than with `\thinspace` `\\;` `\bigl\\{`. Switch to a
     fenced block and write the natural TeX.
  2. *Awkward markdown context* — fenced blocks survive list-item
     nesting, blockquote nesting, and `<details>` better than `$$...$$`.
     The structural pass already suggests this as one fix when a
     multi-line `$$` block is inside a list item.

  Trade-offs: extra surrounding syntax, no inline use (display-only),
  and visual diff churn if you switch a long-established `$$...$$`
  block.

  The linter applies the static pass (forbidden macros) and both
  render passes to fenced content, but skips the GFM pass — fenced
  math is exempt by design.
- Bold math: `\boldsymbol{x}` or `\mathbf{x}`, not `\bm{x}`.
- Inline math: write `$x$5` carefully — GitHub treats `$` adjacent to digits
  inconsistently. A space (`$x$ 5`) avoids the problem entirely.
- **Inline math with `}_{` (subscript right after a brace): wrap in
  backtick-dollar.** GitHub's markdown preprocessor treats a `_`
  preceded by `}` as the start of an italic span and eats the
  underscore before MathJax sees the math. So `$V^{(2)}_{\vec q}$`
  inside flowing prose renders as the literal string
  `$V^{(2)}{\vec q}$` (underscore gone, dollar signs visible).
  This is community-discussion
  [#65772](https://github.com/orgs/community/discussions/65772).

  The fix is GitHub's documented alternative inline-math syntax,
  `$`...`$` (backtick-dollar): the backticks make the content a code
  span as far as markdown is concerned, so the inline emphasis rule
  is skipped entirely. Rewrite:

  | Don't write | Write instead |
  | --- | --- |
  | `$V^{(2)}_{\vec q}$` | `` $`V^{(2)}_{\vec q}`$ `` |
  | `$\|\Gamma^{(2)}_{\vec q}\|/\hbar$` | `` $`\|\Gamma^{(2)}_{\vec q}\|/\hbar`$ `` |

  Inline math without the `}_{` pattern (e.g. `$\vec r_{ij}$`,
  `$V_2$`) is fine as plain `$...$`. Display math `$$...$$` is also
  not affected by this rule.

  The linter's GFM pass enforces this.
