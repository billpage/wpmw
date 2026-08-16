# docs/analysis

Mathematical derivations and detailed review notes.

This directory holds the worked-out mathematics behind the algorithms. It
began as review notes on the source documents — term-by-term matches between
the extended Fokker–Planck (xFP) expansion and the Moyal series for the Wigner
Quantum Liouville Equation, and moment-problem analyses for jump densities
under polynomial potentials. Most of it is now a **derivation ladder**: a
sequence of notes, each taking as input something its predecessor postulated,
and each ending with the open items that motivate the next.

## The ladder

Read in this order. Each note states in its §0 what it inherits and what it
retracts.

1. **[`phase_space_crystal_lattice_review.md`](phase_space_crystal_lattice_review.md)** — Review of the two source
   documents (Cyganski's *Extended Fokker–Planck Eq. and the QLE V2* memo and
   the *Wigner Collisions Diagram* Sozi deck), cross-referenced at the
   equation-and-page level. The entry point for anyone tracing a claim back to
   its origin.
2. **[`four_rule_microdynamics_equivalence.md`](four_rule_microdynamics_equivalence.md)** — Analysis of Cyganski's
   proposal (Zoom, 2026) to replace the single mediated-jump rule with four
   two-body rules (Focus, Defocus, Right-Hop, Left-Hop). Proves exact
   equivalence at any particle number `ν`, shows the four-rule form is ≈5.6×
   quieter, and identifies the `G`-freedom — a family of exact rate
   assignments of which the single rule is the `G = 0` member. Ends on a no-go
   lemma: pairwise mass action among tracked particles is quadratic in
   occupancy while the QLE generator is linear, so a fully collision-based
   microdynamics needs a species whose density is *pinned*.
3. **[`sea_dressed_microdynamics.md`](sea_dressed_microdynamics.md)** — Takes the step that lemma leaves open.
   Realises the collision term as sixteen local, two-body,
   momentum-conserving channels against a pinned Dirac sea of positon–negaton
   pairs, exact at pinned sea. Postulates the sea's polarisation: the rate
   field `Γ_q(x)`, its sign structure, and the half-quantum stencil offsets
   all enter as assumptions.
4. **[`phase_resonance_microdynamics.md`](phase_resonance_microdynamics.md)** — Derives that polarisation rather
   than postulating it, by making phase a particle-level property (P0–P5,
   Theorems 1–3). Contains the parity result (fundamental particles occupy
   even momentum sites), the rate-table no-go (Theorem 2: phase-blind
   transition rules cannot reproduce linear rates), and the dark-sea lemma.
5. **[`phase_alignment_microdynamics.md`](phase_alignment_microdynamics.md)** — A change of variables on the
   predecessor: the beat, the grating and the resonance condition are replaced
   by a single scalar, the misalignment `μ` of two transported clock phases.
   No new postulates and no different predictions, but Theorem 4 is stronger
   than what it replaces: requiring `μ` to hold still through a vertex forces
   the vertex to be a **momentum swap**, from which energy conservation and
   the selection rule follow rather than being imposed.
6. **[`relational_pairing_and_carrier_lock.md`](relational_pairing_and_carrier_lock.md)** — Removes stored partnership
   from the algorithm specification (§2.2), at the cost of one postulate (S),
   the sea carrier lock. Proposition R1 shows a partner index carries no
   relational state; Theorem R4 factorises the vertex weight through a
   per-cell, per-row order parameter `Z_r`, cutting the encounter loop from
   `O(N_exc · B)` to `O(N_exc + N_sea)`; and §8 records a defect the indexed
   formulation concealed — under permanent partnership the sea is a
   consumable resource with no source, short by ≈770× for the cosine-well
   parameters.
7. **[`permanent_pairing_density_matrix.md`](permanent_pairing_density_matrix.md)** — Reinstates permanent pairing
   under a density-matrix reading: a pair is a sampled element of ρ (positon =
   ket leg, negaton = bra leg, μ = arg ρ), an excess particle a diagonal
   sample. Retracts two *inferences* of its predecessor while keeping its
   calculations: pairing adds nothing to the one-particle marginal but is the
   entire content of the two-point function; and the pump — whose sidebands
   *are* split-pair amplitudes — is the source Corollary R5.1 declared
   absent. Derives the stencil and the mediated counting from one-leg hops of
   bound pairs, shows storage feasibility is exactly the Wigner bound
   |W| ≤ 2/h (the same inequality as W' ≥ 0), makes postulate (S)
   unnecessary, and leaves one load-bearing theorem: split pairs mediate with
   the same vertex constant as pump-excited pairs.
8. **[`coherence_ladder.md`](coherence_ladder.md)** — Indexes ρ by splitting rung and derives the
   complete first-order channel table of the pairing vertex from
   stationarity: exact (K3), leg-local ladder, and compound classes. Proves
   the ladder theorem — four leg-local channels (struck leg ket/bra ×
   direction), sea strikers for interior edges, phase continuity, conjugate
   bra factors and one constant reproduce the commutator elementwise at
   machine precision, freezing exactly at V = 0. Corrects five earlier
   statements, including the spec's "sea–sea optional" (the sea is the
   interior engine) and its erase amplitude (all transfers carry μ₁, by the
   licensing argument). Leaves open: striker back-reaction neutrality and
   compound-channel cancellation.
9. **[`position_pair_ladder.md`](position_pair_ladder.md)** — The same construction in the *position*
   representation, ρ(X, X'), where each leg carries a place and a clock and
   nothing else. The potential and the kinetic operator swap jobs: the
   potential only winds μ, and all motion is four one-leg hops of amplitude
   ±iJ/ħ. Momentum turns out not to be carried at all — it is the
   misalignment of a nearest-neighbour conjugate pair, p̄ = ħμ/a — and the
   Euler force term follows from the pump alone. Corrects four statements of
   the density-matrix algorithm specification, in particular that the
   continuum's missing uniform bound is restored on the lattice, and that the
   spec's mean-field pair-Bohm machinery is unnecessary. Ends on two
   obstructions: hop *probability* is second order in δt, and local gauge
   sweeps arg ρ around the whole circle, so no positon-only sea exists here.
   The observable sector survives as a genuine positive-rate particle process
   guided by sin μ; the coherence sector does not.
10. **[`fourd_microdynamics.md`](fourd_microdynamics.md)** — The ladder pushed into four-dimensional
    phase space: two particles on a line, and one particle in the plane.
    Everything above generalises under one substitution — the shift operator
    becomes shift-by-`q` for a *vector* wavevector — so the exactness family,
    the participant-locality selection of the symmetric member, and endpoint
    locality all survive unchanged (Theorem A1). What does not survive is the
    4-D supplement's framing: Theorem A3 shows a momentum direction is
    conserved iff it is orthogonal to every active mode wavevector, so the
    2p/1D-versus-1p/2D distinction dissolves — the two are the same equation,
    verified to zero. Corrects the four-rule note's claim that the hop channel
    cannot be a particle–particle exchange (true for external modes only), and
    shows that Theorem 4 of the phase-alignment note is unique **only in one
    spatial dimension**: in d ≥ 2 its conditions leave a (d−1)-parameter
    family whose every member conserves energy and does not dephase, so the
    swap must be postulated (exchange-only). Establishes that the sea is not
    merely convenient but the *only* available collision partner once N > 1,
    that the crystal shift does not commute with products, and that peak
    |W|/(2/h)^d is the state purity — so entanglement is literally
    excess-to-background loss. Closes with the harmonic potential: exactly
    classical in the mean (the Moyal bracket truncates), and the
    microdynamics' worst case, since the per-mode injected momentum variance
    is 2mω²ℏ independent of q and the total therefore grows linearly in the
    mode cutoff.
11. **[`open_position_space.md`](open_position_space.md)** — What happens when position space is not
    closed. Separates the two jobs the periodic box has been doing: fixing
    the momentum quantum, and keeping worlds in view. Only the first is real,
    and it has three independent sources — ring circumference, coherence
    horizon, and periodicity of `V` — of which the last works on all of ℝ.
    Theorem O1: the Wigner kernel's modulus is independent of position for
    *every* potential, so a world far from a localised scatterer is struck at
    full rate and its free behaviour is a cancellation whose fringe frequency
    grows as 2x/ℏ; the escape problem is therefore not a boundary problem and
    no absorber can fix it. Theorem O5: for an `a`-periodic potential every
    world's `p mod (πℏ/a)` is exactly conserved with no box at all, so an
    L = na ring is not an approximation to open space but exactly `n` of its
    sectors — the modular companion to Theorem A3. Corrects the algorithm
    specifications on both points, retracts open item 4 of the inverted-pair
    supplement, and records that the sea-dressed layer does *not* generalise
    naively, its constant background having infinite total on non-compact
    phase space.

12. **[`interworld_coupling.md`](interworld_coupling.md)** — Why four rules, and why not more. Reads
    the potential as a coupling between the two legs of a position pair,
    `U = V(x₁) − V(x₂)`, in midpoint and *full* separation coordinates
    `X = (x₁+x₂)/2`, `Y = x₁−x₂`. Proposition I1: the coupling vanishes at
    coincidence, is antisymmetric under leg exchange, and vanishes
    identically for a free particle — the last being the sharpest available
    test of any proposed interworld force law, and the point of difference
    from many-interacting-worlds models, whose interworld potential is what
    makes a free packet spread. Proposition I2: for one cosine mode the
    coupling factorises into a midpoint amplitude (the classical force) and a
    separation grating of *twice* the potential's period, because each leg
    moves only `Y/2` — which derives the half-quantum offset that the
    sea-dressed note postulates. Theorem I3: the available momentum channels
    are exactly the Fourier spectrum of the coupling in `Y`, so they are
    discrete iff the coupling is periodic; hence `M` modes give `2M` shifts
    and `4M` rules, and a one-mode potential has exactly four. Theorem I4:
    the Moyal series *is* the odd-power expansion in `Y`, so a coupling
    linear in the separation is exactly classical — which is why the harmonic
    and inverted harmonic alike have no jump channel, and why the inverted
    pair barrier cannot test the four rules. Theorem I5: the coupling exerts
    no force and does no work; it is the winding rate of the misalignment,
    `dμ/dt = −U/ℏ`. Contradicts the hypothesis that prompted it — the four
    rules need no non-linear interworld force, the apparent four-wave mixing
    being the bilinearity of `ρ` rather than a `χ⁽³⁾` medium — and states
    plainly which of its claims are theorems, which are Fourier-dual
    restatements, and which is the interpretive postulate that reading `Y` as
    a physical separation requires.

13. **[`species_phase_duality.md`](species_phase_duality.md)** — What the
    relationship is between world-particle *species* (positon/negaton) and
    world-particle *phase*, and what it costs to add an annihilation process.
    Theorem D1: species and phase are one degree of freedom in conjugate
    bases — phase is a coordinate on the leg separation `Y`, sign is a value
    on the momentum `p` — so no carrier holds both sharply, and the sign is a
    functional of a whole `Y`-fibre rather than an attribute of a pair. The
    decisive case is the cat at rest, whose pair phase is identically zero
    everywhere and whose Wigner negativity is 0.29. Theorem D2: the crystal
    shift is the *identity operator*, `ρ → ρ + 2·1`, dark because `[H,1] = 0`
    for every Hamiltonian in every dimension — not because anything cancels;
    in the pair basis it is `δ(Y)`, so the sea sits at exactly zero leg
    separation, which is why its darkness depends only on relative position.
    Theorem D3 extends I2/I5 to arbitrarily many modes,
    `dμ/dt = −2 Σ_q Γ_q(X) sin(k_q Y/2)`, whose maximum is the `γ_max` of the
    annihilation burden — so the pathwise `L¹` growth rate is a phase
    precession rate — and Corollary D4.1 thereby *derives* the sea
    polarisation the sea-dressed note postulates. Theorem D5: annihilation
    exact in `(X, p)` is non-local over the *coherence length* in the leg
    positions and independent of `Δx`, and Corollary D5.1 licenses anonymous
    annihilation for all future time, correcting §7.4 of the representation-
    cost supplement. Theorem D6 prices the soft blob: a momentum kernel of
    width `σ_p` is exactly an imposed coherence length `ħ/σ_p`. §7 records
    that no unraveling linear in the ensemble is `L¹`-stationary (Proposition
    U1) and that the governing exponent is `ρ(|L|) = 2.341`, settling N1
    negatively. §8 specifies the annihilation substep for the algorithm notes
    and §9 measures it: 0.32 per cent wall-clock overhead against a benefit of
    `exp(2ρt)` in particle count, exactly unbiased in every moment, with no
    plateau over three decades of ensemble size.

## Companion code

Every note lists its verifying script in §0 and its numerical results in a
late section. In ladder order: `src/demo_four_rule_equivalence.py`,
`src/demo_sea_dressed_dynamics.py`, `src/demo_phase_resonance_rates.py`,
`src/demo_phase_alignment.py`, `src/demo_relational_pairing.py`,
`src/demo_pairing_resource_arithmetic.py`,
`src/demo_coherence_ladder.py`,
`src/demo_position_pair_ladder.py`,
`src/demo_fourd_microdynamics.py`,
`src/demo_open_position_space.py`,
`src/demo_interworld_coupling.py`,
`src/demo_species_phase_duality.py`. All
non-trivial claims in these notes are verified numerically before inclusion.

## Related

Drawings of the elementary processes of the phase-alignment layer are in
[`../supplement/phase_alignment_interaction_diagrams.md`](../supplement/phase_alignment_interaction_diagrams.md); the simulation
counterparts of the ladder are in [`../algorithm/`](../algorithm/README.md).
