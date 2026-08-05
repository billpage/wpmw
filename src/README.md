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
  forms.  Also implements the four-rule (Focus/Defocus/Right-Hop/Left-Hop)
  jump forms — `step_jump_four_rule` (mesh) and `step_jump_four_rule_mc`
  (particle MC) — specified in
  `docs/analysis/four_rule_microdynamics_equivalence.md`.
- `wpmwlib/sea_dressed_lattice.py` — `SeaDressedLattice`, a subclass of
  `PhaseSpaceCrystalLattice` implementing the sixteen-channel two-body
  microdynamics of `docs/analysis/sea_dressed_microdynamics.md`. Carries three
  integer fields per cell — unpaired positons `U+`, unpaired negatons `U-`,
  and ground-state neutral pairs `S` — with the observable excess
  `E = U+ - U-` and sea background `S_bar = B`, the integer image of the
  crystal shift 2/h. Exposes the Monte-Carlo channel step
  (`step_jump_sea_mc`), the mesh-form channel generator
  (`channel_generator_mesh`) used for the exactness check, a recombination
  step (`step_recombine`, which leaves `E` invariant and so lives entirely
  outside the QLE generator), and the diagnostics `sea_min_fraction`,
  `unpaired_total` and `worldline_invariants`.
- `wpmwlib/wigner_split_fourier.py` — reference solver: Strang-split spectral
  Fourier on the Wigner equation. Specialized to QHO; for general V the
  force-kick kernel must be replaced by the full Wigner–Moyal kernel.
- `wpmwlib/check_md_math.py` — markdown LaTeX-math linter (see "Markdown
  math linter" below).

## Runnable scripts

Grouped by role. Within each group the order is the order to read them in.

### Solver demos

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
### Microdynamics ladder

Each demo verifies one rung of the derivation ladder in `docs/analysis/`;
see that directory's README for the ladder itself.

- `demo_four_rule_equivalence.py` — verification companion to
  `docs/analysis/four_rule_microdynamics_equivalence.md`, which analyses the
  proposal (D. Cyganski, 2026) to replace the single mediated-jump rule with
  four two-body rules (Focus, Defocus, Right-Hop, Left-Hop) carrying signed,
  occupancy-dependent bias rates.  Part A checks the generator identity: on
  random mesh fields and for mode sets q = 1, 2, {1,2,3}, one Euler step of
  the independently-assembled four-rule mesh form agrees with the original
  single-rule stencil to machine precision (worst deviation ~1e-16).  Part B
  evolves a squeezed Gaussian in the cosine well three ways on the same
  64² lattice with an identical splitting sequence — deterministic mesh QLE,
  single-rule MC, four-rule MC (ν = 1.6×10⁶, two classical periods) — showing
  both MC runs track the mesh at their shot-noise floors, with the four-rule
  floor ≈ 5.6× lower at equal ν because its bias rates are built from
  background-subtracted excess counts.

  Sample output figures (committed on the `output` branch):

  ![Four-rule equivalence evolution](https://raw.githubusercontent.com/billpage/wpmw/output/figures/four_rule_equivalence_evolution.png)

  4×4 grid: mesh QLE, single-rule MC, four-rule MC, and pointwise
  (four-rule − mesh) at t = 0, T/2, T, 2T.  The single-rule row carries
  uniform background shot noise (its mediator count includes the 2/h
  background); the four-rule row is visibly cleaner.

  ![Four-rule equivalence metrics](https://raw.githubusercontent.com/billpage/wpmw/output/figures/four_rule_equivalence_metrics.png)

  Relative L² deviation from the mesh for both MC runs (left) and Wigner
  negativity for all three evolutions (right).
- `demo_four_action_foundations.py` — verification companion to
  `docs/supplement/four_action_foundations.md`, which audits the claim
  (D. Cyganski, 3 August 2026) that momentum and energy balance determine the
  four-action rates.  Six parts.  A: every member of the exact rate family
  conserves number, momentum and energy (worst residual 2.5e-14), so
  conservation is implied by exactness and cannot select a member.  B: on the
  conserving line b = -1/2 the QLE residual is 0, 1.04, 3.11 for a = 1/2, 0, 2
  while every conservation residual stays at machine epsilon — a = 0 is a
  fully conserving, Ehrenfest-exact model that is not the QLE.  C: the focus
  gain is invisible to moments 0-2 (1e-16) and appears at moment 3 (8e-4).
  D: with a free momentum grid step the model is exactly Moyal evolution at
  hbar_eff = 2*delta/k, matching to 1e-14 across a factor of 40 in delta and
  equalling the true QLE only at delta = hbar*k/2.  E: only G = Gamma/2 leaves
  both rates supported on the endpoint cells.  F: a two-body mass-action rate
  linearised about a uniform sea is endpoint-symmetric with exactly zero
  antisymmetric part, so it cannot produce the focus rate.  Figures:
  `four_action_uniqueness_map.png`, `four_action_hbar_effective.png`.

- `demo_sea_dressed_dynamics.py` — verification companion to
  `docs/analysis/sea_dressed_microdynamics.md`, which takes the step left
  open by the four-rule note's no-go lemma: the occupancy stencils are
  realised as sixteen local, momentum-conserving, two-body collision
  channels between world-particles, with a Dirac sea of positon–negaton
  pairs supplying the pinned-density reservoir that linearity requires.
  Part A assembles the sixteen-channel mean-field generator channel by
  channel — crossing conjugates included — on *independent* random `U+` and
  `U-` fields and compares it against the original single-rule stencil
  acting on `E = U+ - U-`; agreement to machine precision verifies the
  channel table's bookkeeping.  Part B evolves the same squeezed Gaussian in
  the cosine well as the four-rule demo six ways: deterministic mesh QLE,
  four-rule MC as the noise-floor reference, sea-dressed MC against a pinned
  sea, and sea-dressed MC with a live ledger at three recombination rates
  κ_rec.  With κ_rec = 0 the sea drains and the orphan load inflates —
  every capture, split and emission event orphans a partner — and increasing
  κ_rec converges the live run onto the pinned one.  The worldline
  invariants (`U+ + S` and `U- + S` totals) are asserted constant across
  every jump and recombination event.

  Sample output figures (committed on the `output` branch):

  ![Sea-dressed evolution](https://raw.githubusercontent.com/billpage/wpmw/output/figures/sea_dressed_evolution.png)

  4×4 grid at four snapshot times: mesh QLE, sea-dressed MC against a pinned
  sea, and live-sea MC at the smallest and largest κ_rec.  The suptitle
  carries the Part A generator-identity error alongside ν.

  ![Sea-dressed metrics](https://raw.githubusercontent.com/billpage/wpmw/output/figures/sea_dressed_metrics.png)

  Three panels over the run: relative L² deviation from the mesh (fidelity),
  worst-cell sea depletion min S/B, and the unpaired population `U+ + U-`.
  Read together they show the failure mode of an unreplenished sea and its
  repair by recombination.
- `demo_phase_resonance_rates.py` — verification companion to
  `docs/analysis/phase_resonance_microdynamics.md`, which derives the sea's
  polarisation instead of postulating it by giving world-particles a de
  Broglie phase (P0–P3).  Part A checks the beat kinematics: two legs split
  by 2q·dp beat as a full-contrast grating at exactly the mode-q wavelength,
  drifting at the mean velocity p̄/m (Proposition 1).  Part B checks Bragg
  selection: one split-operator kick of a plane wave scatters to exactly
  p₀ ± 2q·dp with amplitude i·V_p·dt/(2ħ) per sideband — the full quantum,
  linear in V_p, with the refractive factor i.  Part C checks the midpoint
  identity: the exact first-order Wigner change of the kicked plane wave
  equals the single-cosine QLE stencil at the midpoint sites n₀ ± q, so the
  stencil's half-quantum offsets are the interference midpoints of a
  full-quantum transfer.  Part D is a rate-table toy in which sea pairs carry
  particle-level data only, verifying that the per-channel rate, its
  quadrature, the factor γ/2 and the direction field σ are all *derived*
  (Theorem 3); that the pattern-phase spread across pairs pumped at random
  locations is zero to machine precision (Lemma 2 — coherence for free); and
  that time-averaging along the transition worldline dephases every row
  except the resonant one.

  Sample output figures (committed on the `output` branch):

  ![Phase grating of the sea](https://raw.githubusercontent.com/billpage/wpmw/output/figures/phase_grating_sea.png)

  One pair's beat, the coherence of the grating across pairs, and the
  quadrature relation to V(x).

  ![Phase grating in space-time](https://raw.githubusercontent.com/billpage/wpmw/output/figures/phase_grating_spacetime.png)

  Beat crests in the (x, t) plane and the refraction of a worldline crossing
  them.

  ![Phase-resonance rate law](https://raw.githubusercontent.com/billpage/wpmw/output/figures/phase_resonance_rate_law.png)

  The derived rate law and the row-selection rule.
- `demo_contact_vertex_reduction.py` — companion to §13 of
  `docs/analysis/phase_resonance_microdynamics.md`, which replaces the P5
  weight by a two-level unitary contact interaction
  `h = g0 + g1·C·exp(i δ)`, `P_flip = sin²(|h| τ_e)`, with δ the struck
  beat's pattern phase at the vertex.  No phase offset is fitted anywhere.
  R1: the note's rate law re-emerges with the analytically predicted
  coefficient and the correct quadrature sign, δ₀ = 0 being *derived* from
  hermiticity rather than matched.  R2: the residual scales as O(C³),
  because the O(C²) term of |h|² is direction-symmetric and cancels in the
  net, so linear response is exact at first order.  R3: g0 → 0 kills the
  response — "no noise, no force" — since the coefficient tracks
  sin(2 g0 τ_e), making the bare phase-blind exchange traffic the carrier of
  the quantum force rather than removable noise.  R4: with two pumped modes
  the gross rate acquires an O(C²) spatial modulation at the difference
  wavevector, which is the signature that discriminates the amplitude vertex
  from the bare affine P5 (which predicts exactly zero).

  Sample output figure (committed on the `output` branch):

  ![Contact-vertex concept](https://raw.githubusercontent.com/billpage/wpmw/output/figures/contact_vertex_concept.png)

  Anatomy of the contact, the phasor sum, the emerging rate law, and the
  consequence tree.
- `demo_phase_alignment.py` — verification companion to
  `docs/analysis/phase_alignment_microdynamics.md` and
  `docs/algorithm/phase_alignment_microdynamics_algorithm.md`.  The note
  eliminates the beat, the grating and the resonance condition in favour of
  one scalar carried by a pair, the misalignment μ of two transported clock
  phases; this demo checks that the replacement loses nothing.  Part A
  (Lemma 4): μ is invariant under a global phase shift and under
  re-referencing either worldline, and the pair amplitude depends on the
  partners only through μ.  Part B (Proposition 3): μ winds in space at
  ∂μ/∂x = Δp/ħ and along a path of velocity v at
  μ̇ = (Δp/ħ)(v − v̄_pair), so "beating" is nothing but Δp ≠ 0 and no
  propagating object is involved.  Part C is the exchange theorem —
  requiring μ to be stationary along the transition path, together with
  momentum conservation, forces the vertex to be a momentum *swap*, after
  which energy conservation is automatic rather than imposed.  Parts D–E
  check the locality of μ under the pump and the quadrature and linearity of
  the vertex bias; Part F checks that the resulting stencil preserves L¹ on
  the ring.

  Sample output figure (committed on the `output` branch):

  ![Phase-alignment contact interaction](https://raw.githubusercontent.com/billpage/wpmw/output/figures/phase_alignment_contact.png)

  Six panels: transporting two clocks to a common point to define μ; the
  three states of a pair; the pump making μ a function of place; the vertex
  as a momentum swap; the stationarity condition; and the ensemble limit.
- `demo_relational_pairing.py` — verification companion to
  `docs/analysis/relational_pairing_and_carrier_lock.md`, which removes the
  stored partner index of §2.2 of the phase-alignment specification in
  favour of a misalignment defined over *every* ordered pair of
  world-particle indices.  Part A shows μ_ij is antisymmetric and satisfies
  the cocycle identity exactly, so the two-index family is the coboundary of
  the one-index field Φ and a stored index carries no relational state.
  Part B exhibits the obstruction: a free per-pair carrier phase suppresses
  the all-pairs average by 1/B while leaving the partnered average
  untouched.  Part C imposes the carrier-lock postulate (S) — a common
  transported phase per (cell, row), the literal reading of "phase-space
  *crystal* lattice" — and recovers exact equivalence.  Part D is the
  factorisation: under (S) the double sum collapses onto a per-cell,
  per-row order parameter Z, turning the O(N_exc · B) vertex loop into
  O(N_exc + N_sea) with a measured speedup of ≈6850× at N = 8192.  Part E
  reassembles Γ_q(x) from Z alone with the carrier cancelling identically.
  Part F is the deciding experiment: under per-particle noise the two forms
  agree in the mean and differ by √N in variance, while under
  pair-correlated noise they differ *in the mean*.  Part G sizes the
  consumable-sea deficit for the cosine-well parameters, and Part H
  separates the two labels in (S) — an exact momentum row against a
  sweepable spatial cell.

  Sample output figures (committed on the `output` branch):

  ![Relational pairing](https://raw.githubusercontent.com/billpage/wpmw/output/figures/relational_pairing.png)

  The obstruction of Part B, the stencil rate assembled from Z alone against
  the analytic Γ_q(x), and the variance separation of Part F.

  ![Lattice cells and rows](https://raw.githubusercontent.com/billpage/wpmw/output/figures/lattice_cells_and_rows.png)

  Exact momentum rows against sweepable spatial cells, and the intra-cell
  pump spread over a full cosine-well run.
- `demo_pairing_resource_arithmetic.py` — verification companion to
  `docs/analysis/permanent_pairing_density_matrix.md`, which revisits
  permanent positon–negaton pairing under the density-matrix reading (pair =
  sampled element of ρ; positon = ket leg, negaton = bra leg).  Runs the
  exact split-operator Schrödinger arm on the canonical cosine well and
  compares two accountings of the sea's aligned-pair budget.  Check 1: the
  local load factor |W|·(h/2) never exceeds 1 — Wigner's bound, the same
  inequality that makes W' = W + 2/h non-negative — with a narrow-packet
  control saturating it to 6 digits.  Check 2: the global off-diagonal mass
  C(t) of ρ(P,P') stays bounded (peak 10.3) against one state-mass unit of
  aligned stock per momentum row, with half the coherence within 4 one-leg
  hops of the diagonal.  Check 3: the consumable model of the relational
  note's §8 exhausts the busiest cell at t* ≈ T/76 on accepted events alone,
  reproducing the ≈770× shortfall ν-independently.  Check 4: flow
  feasibility is the triangle inequality on W' ≥ 0.

  Sample output figure (committed on the `output` branch):

  ![Pairing resource arithmetic](https://raw.githubusercontent.com/billpage/wpmw/output/figures/pairing_resource_arithmetic.png)

  The load-factor map capped by the Wigner bound, bounded storage demand
  against linearly growing consumable demand, and the off-diagonal
  density-matrix mass the split pairs must carry.

- `demo_coherence_ladder.py` — verification companion to
  `docs/analysis/coherence_ladder.md`, which indexes ρ(P,P') by splitting
  rung k = (P−P')/(2q·dp) and derives the vertex channel table from
  stationarity.  Part A confirms the table: the exact channel (own winding,
  sinc = 1 at p_in = mate's row) plus the leg-local ladder and compound
  sideband channels, and nothing else.  Part B assembles the four
  leg-local channels — struck leg ket/bra × direction, sea strikers, phase
  continuity, refractive −i·e^{±iφ} with conjugate bra factors, one shared
  constant c = V_q/2ħ — and matches the exact commutator −(i/ħ)[V,ρ] on a
  random mixed state to 1.7×10⁻¹⁶ relative, Hermitian to 1.4×10⁻¹⁷, trace
  drift 2.4×10⁻¹⁸, with the exact V = 0 freeze and zero diagonal flux on a
  coherence-free state (Lemma 3 at the ladder level).  Part C draws the
  ladder diagram.

  Sample output figure (committed on the `output` branch):

  ![Coherence ladder](https://raw.githubusercontent.com/billpage/wpmw/output/figures/coherence_ladder.png)

  The density matrix as a ladder of splitting rungs with the four one-leg
  hops, and the striker roles: excess for the population boundary, the
  background sea for interior rungs.

- `demo_position_pair_ladder.py` — verification companion to
  `docs/analysis/position_pair_ladder.md`, the same construction in the
  position representation ρ(X,X') on a periodic lattice of M sites.  Part A
  checks Theorem P1: the von Neumann generator is exactly four one-leg hops
  of amplitude ±iJ/ħ (sign by species) plus a diagonal pump −i·ΔV/ħ, matching
  the commutator to 1.8×10⁻¹⁶ relative, exactly Hermitian and traceless, with
  every modulus and population exactly frozen at J = 0.  Part B checks
  Theorem P2 — momentum is the misalignment of a rung-1 pair, p̄ = ħμ/a, with
  j = (ħ/ma)|ρ₁|sin μ and exact lattice continuity to 1.1×10⁻¹⁶.  Part C
  recovers the Euler force term from the pump alone, second order in a over
  six refinements.  Part D exhibits the two obstructions: a diagonal state has
  exactly zero population flux, free spreading is ballistic
  (dlog⟨x²⟩/dlog t = 2.000), hop probability is O(δt²), and local gauge
  sweeps arg ρ uniformly around the circle.  Part E measures the ℓ¹ resource
  arithmetic against its bounds and tabulates the channel rate Λ = 4J/ħ
  against the lattice-independent momentum-side amplitude.  Part F runs the
  complex-weight pair ensemble: unbiased, 1/√N, noise amplified by exactly
  e^{Λt}.  Part G runs 200,000 self-conjugate walkers guided by sin μ, which
  track the exact probability density to the shot-noise floor.

  Sample output figure (committed on the `output` branch):

  ![Position-space coherence ladder](https://raw.githubusercontent.com/billpage/wpmw/output/figures/position_pair_ladder.png)

  Eight panels: the pair board, momentum from misalignment, the Euler force
  term, the amplitude-versus-probability obstruction, the gauge circle, the
  resource arithmetic, the pair-ensemble Monte Carlo, and the guided
  self-conjugate walkers.

### Figure generators and regression tests

- `gen_microdynamics_4d_figures.py` — generates the five schematic
  (blackboard-style) figures for
  `docs/supplement/phase_space_crystal_lattice_4d_supplement.md`, in the
  visual idiom of the Sozi slides reproduced in the crystal-lattice
  supplement: red positons, green negatons, blue mediator.  Contrasts two
  particles in 1 spatial dimension against one particle in 2 spatial
  dimensions, which share a 4-dimensional joint phase space but slice it
  differently.  No physics is computed here — the script draws, it does not
  verify.

  Output figures (committed on the `output` branch):

  ![4D layouts](https://raw.githubusercontent.com/billpage/wpmw/output/figures/microdynamics_4d_layouts.png)

  The two natural slicings of the shared 4-dimensional joint phase space.

  ![4D Fourier modes](https://raw.githubusercontent.com/billpage/wpmw/output/figures/microdynamics_4d_fourier_modes.png)

  Mode support: the 2p/1D pair potential lives on a 1-D line in joint
  wavevector space, while 1p/2D modes can fill all of Z².

  ![4D jumps](https://raw.githubusercontent.com/billpage/wpmw/output/figures/microdynamics_4d_jumps.png)

  Per-event momentum-space displacement: 2p/1D jumps are confined to the
  anti-diagonal; 1p/2D jumps point in arbitrary directions.

  ![4D centre-of-mass and relative coordinates](https://raw.githubusercontent.com/billpage/wpmw/output/figures/microdynamics_4d_com_relative.png)

  The same jump rule in the original (p₁, p₂) basis and in the rotated
  (P, p_rel) basis, where the separation is manifest — the algebraic shadow
  of "2p/1D is hidden 1+1-D".

  ![4D starburst](https://raw.githubusercontent.com/billpage/wpmw/output/figures/microdynamics_4d_starburst.png)

  Crystal-lattice mediator picture: a virtual quantum exchanged between the
  two particles in 2p/1D, against a single particle absorbing a vectorial
  kick in 1p/2D.
- `gen_phase_alignment_interaction_diagrams.py` — companion to
  `docs/supplement/phase_alignment_interaction_diagrams.md`.  Draws every
  elementary process of the phase-alignment layer twice, in the format of §4
  of the crystal-lattice supplement: once as a lattice picture in the (x, p)
  plane and once as trajectories in the (x, t) plane.  Three checks run
  before anything is drawn.  Part A confirms the vertex is a swap: the
  two-condition system returns p_in = p_b and p_out = p_a at zero residual,
  Σp and Σp² are conserved without either being imposed, and the union of
  the two worldlines is *identical* to two straight lines crossing — the
  exchange permutes labels, not trajectories.  Part B confirms the sinc
  dephasing envelope for off-stationary candidates against direct numerical
  averaging.  Part C builds a sea row as a ring amplitude and reads the
  misalignment off the Kapitza–Dirac sidebands, whose relative amplitude
  comes out at μ₁/2 as predicted; this is the check that exposes the sign
  discrepancy recorded in §7 of the document, where the kick sign of §6 of
  the algorithm specification reproduces −Γ_q(x) while the opposite sign
  reproduces +Γ_q(x) and Lemma 5 as printed.

  Output figures (committed on the `output` branch), one per process, each
  with a space–momentum panel and a space–time panel:

  ![Free leg](https://raw.githubusercontent.com/billpage/wpmw/output/figures/pa_int_0_free_leg.png)

  Process 0, the free leg: momentum is constant and only the carried clock
  advances, so the worldline is straight through any potential whatever.

  ![The pump](https://raw.githubusercontent.com/billpage/wpmw/output/figures/pa_int_1_pump.png)

  Process 1, the pump: the only place the potential touches a
  world-particle, and it touches phase alone.  Nothing moves in the (x, p)
  panel; the clock kick differs from place to place because V does.

  ![Exchange, s = +1](https://raw.githubusercontent.com/billpage/wpmw/output/figures/pa_int_2_exchange_up.png)

  Process 2, the exchange at s = +1: the excess particle hops up by 2q·dp
  and one sea leg hops the other way.

  ![Exchange, s = −1](https://raw.githubusercontent.com/billpage/wpmw/output/figures/pa_int_3_exchange_down.png)

  Process 3, the same vertex run downward, drawn with p_out = 0 so that the
  label-permutation reading is hard to miss.

  ![Suppressed channels](https://raw.githubusercontent.com/billpage/wpmw/output/figures/pa_int_4_suppressed.png)

  Process 4, the suppressed channels: nothing is forbidden, and misaligned
  traffic simply averages away.
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
will not render correctly on GitHub. It catches the failure modes we have
actually hit:

- inline `$...$` math that GitHub's preprocessor will not recognise as math
  at all: math inside a `*...*` emphasis span, and math whose opening `$` is
  glued to a hyphen or a quotation mark (`-$x$`, `"$x$`). In every such case
  the dollar signs survive to the rendered page;
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

The module is kept in lockstep with the standalone version at
[billpage/GitHubLinter](https://github.com/billpage/GitHubLinter), which
carries the same checks without the WPMW packaging. Fixes flow in both
directions; when they diverge, diff the two files before adding a check to
either.

**One negative result worth recording.** An earlier revision flagged
`` `$...$` `` — backtick *outside* the dollar signs — as broken, on the
theory that GitHub's math pipeline runs before inline-code processing and
would reinterpret the content as math anyway. It does not. Tested upstream
against GitHub's own renderer (`POST https://api.github.com/markdown`,
`mode=gfm`), every such instance came back as a plain `<code>` element with
no `<math-renderer>` wrapper, while the bare form without backticks did
produce one. A code span that forms at all removes its content from
consideration entirely, as CommonMark requires. The check is gone; writing
`` `$...$` `` in prose to *talk about* dollar-math is safe, which is what
the style guide below does throughout.

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
`.github/workflows/check_md_math.yml`, which passes `src/README.md`
explicitly in addition to the defaults — this file documents the math
syntax rules and quotes the broken forms as examples, so it is worth
keeping under the linter's eye. The module's own defaults are left
matching upstream `GitHubLinter` rather than being changed here.

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
- **Inline math with `}_` or `'_` (subscript right after a brace or prime): wrap in
  backtick-dollar.** GitHub's markdown preprocessor treats any `_`
  preceded by punctuation as the start of an italic span — regardless of what
  follows the `_`. All of `}_q` (letter), `}_0` (digit), `}_{`
  (brace), `}_\vec` (command), and `'_i` (prime) trigger the trap.
  The underscore is eaten, the whole `$...$` fails to render, and
  other inline math later in the same paragraph often cascades and
  breaks too.
  Reference: community discussion
  [#65772](https://github.com/orgs/community/discussions/65772).

  The fix is GitHub's documented alternative inline-math syntax,
  `$`...`$` (backtick-dollar). The backticks make the content a
  code span as far as markdown is concerned, so the inline emphasis
  rule is skipped entirely.

  | Don't write | Write instead |
  | --- | --- |
  | `$V^{(2)}_{\vec q}$` | `` $`V^{(2)}_{\vec q}`$ `` |
  | `$\|\Gamma^{(2)}_q(r)\|$` | `` $`\|\Gamma^{(2)}_q(r)\|`$ `` |
  | `$W^{(2)}_0(x, p)$` | `` $`W^{(2)}_0(x, p)`$ `` |
  | `$(X_i, X'_i)$` | `` $`(X_i, X'_i)`$ `` |

  **Important:** any doubled-backslash spacing such as `\\,` or `\\;`
  inside the expression must be simplified to `\,` / `\;` inside the
  backtick-dollar form, because the backticks bypass CommonMark's
  processing — the extra backslash is no longer needed.

  Inline math without a punctuation-then-underscore pattern (e.g.
  `$\vec r_{ij}$`, `$V_2$`) is fine as plain `$...$`. Display
  math `$$...$$` is also not affected by this rule.

  The linter's GFM pass enforces this.

- **Never put `$...$` inside a `*...*` or `_..._` emphasis span.** GitHub
  renders the markdown to HTML *first* and only then scans for `$...$`
  pairs to hand to MathJax. Math that has ended up inside the resulting
  `<em>` is not picked up, and the dollar signs are left on the page
  verbatim — in italics, which makes it look almost deliberate.

  | Don't write | Write instead |
  | --- | --- |
  | `*carries the same $\mu$.*` | `` *carries the same $`\mu`$.* `` |
  | `*linear in $N$*` | `` *linear in $`N`$* `` |

  This bites hardest in **figure captions**, which are italicised by
  house style and often mention symbols; one caption with three `$\mu$`
  in it produces six stray dollar signs.

  Doubled delimiters (`**strong**`, `__strong__`) are *not* known to have
  this problem — the repository has long-standing `**...$x$...**` run-in
  headers that render correctly — so the linter does not flag them. If a
  counterexample ever turns up, widen `_EMPH_SPAN` in the linter.

  The linter's emphasis-span pass enforces this.

- **Inline math with two `^*` (complex conjugate) in the same paragraph.**
  The `*` after `^` is left-flanking per CommonMark and can open
  emphasis. If TWO `^*` expressions appear in the same paragraph,
  the first opens an italic span and the second closes it, eating
  both `$...$` regions between them. The fix is the same:
  `$`...`$` for any expression containing `^*` when another
  such expression is nearby. The linter does not yet detect this
  automatically — watch for it manually.
