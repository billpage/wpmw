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
    Definitions (H) and (R) in §3 fix the coherence horizon `L_c` and the
    **reach** `y_max = L_c/2`, the greatest half ket–bra separation a world
    instantiates and hence the greatest distance at which it consults the
    potential; the reach is the central parameter of step 14.
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

13. **[`species_sectors_and_annihilation.md`](species_sectors_and_annihilation.md)**
    — What world-particle *species* is, what sectors the excess population
    divides into, and what it takes to add an annihilation process. Two
    ensembles are separated first and must not be conflated: E1 draws carriers
    from the quasi-density `W(x, p)` with species the sign of `W` and structure
    group `Z2`, while E2 draws from the density matrix with positon and negaton
    naming the ket and bra legs and structure group `U(1)`. Theorem D0 shows
    the Weyl transform relates the represented objects but *not* the ensembles
    — the species censuses are anti-correlated, so no carrier-level map exists
    and "positon" in the two layers is a homonym. Theorem D1: species and phase
    are one degree of freedom in conjugate bases, the decisive case being the
    cat at rest whose pair phase is identically zero everywhere and whose
    negativity is 0.29. Theorem D2.1: an operator is dark under every
    Hamiltonian iff it is `c·1`, so the neutral sea (`c = 0`, invisible in the
    observable but live in the dynamics) and the crystal shift (`c = 2`,
    visible and provably inert) differ only in `c` — and since `2·1` commutes
    with everything it cannot be the pumped medium the sea-dressed layer needs.
    Theorem D3 extends I2/I5 to arbitrarily many modes. Theorems D8–D10 split
    the excess into an unpaired sector carrying the Born density (the diagonal
    of `rho`, never negative in any position column) and a column-balanced
    sector carrying all the coherence (the off-diagonal), and show the jump
    substep conserves every column sum while streaming conserves every row sum
    — `2(M + N)` free invariants. Theorem D15: the four actions *split and
    combine* bound sea pairs rather than creating them, so positon and negaton
    number are each exactly conserved. Theorem D16 gives the sizing floor in
    units of the Wigner capacity, and Theorem D17 an adaptive per-cell
    allocation that reproduces the accuracy at 0.7 to 1.8 per cent of the
    world-particle count. Section 10 specifies the annihilation substep for the
    algorithm notes, including the orphans-only requirement whose violation
    drives the sea negative while leaving the observable exactly right.

    *This note replaces an earlier step 13, `species_phase_duality.md`; its
    §0.3 records what changed.*

14. **[`compensated_liouville_splitting.md`](compensated_liouville_splitting.md)** — The classical force as
    deterministic acceleration. In the variable `s` conjugate to momentum the
    whole potential term is multiplication by
    `M(x,s) = (i/ℏ)[V(x+y) − V(x−y)]` with `y = ℏs/2` the half ket–bra
    separation — so a world consults V only within its own reach, and no
    Fourier decomposition of V is needed anywhere. Subtracting the part linear
    in `y` splits off the classical Liouville force with the *full* `V′`; both
    factors are diagonal in the same variables, so the factorisation carries
    no Trotter error within the potential substep (C1; §2.1 spells out in
    which representation, why the symbol's independence of `p` is what
    makes the two factors commute, and why the free/potential Strang error
    is untouched; §2.2 says what the operator actually does to worlds —
    a convolution in momentum whose kernel is real and odd, hence signed,
    hence never a one-body Markov jump generator), and the residual is
    exactly the odd part of the
    cubic Taylor remainder of V (C2). Theorem C3 is the point: restricted to a
    bounded coherence reach the residual kernel has zero zeroth *and* first
    moments, so it is a bounded signed jump measure that conserves worlds and
    carries no net momentum — a focus-and-hop that delivers no force, leaving
    the entire classical force in the deterministic step. This is the
    phase-space analogue of the Bohm–Nelson classical/quantum split, with the
    difference that the quantum part is an interaction rather than a force,
    and the price of that difference is a condition on the reach (C4): the
    split gains for `k·y_max ≪ π` and loses beyond `π/2`, and reach and
    momentum quantum are the same parameter, `Δp = πℏ/(2y_max)`. Theorem C5:
    the reach `y_max = L_c/2` of Definition (R) is the one parameter
    everything depends on, and
    a ring pins every world at `u = qπ` for every mode — the two ket–bra arms
    meet at the antipode, so the symbol vanishes and the residual exactly
    cancels the classical term — which is why the reorganisation looks empty
    when tested on a ring, and why a ring is not a valid testbed for the reach
    condition. Theorem C7 is the open-line payoff: if `V‴` vanishes on
    `[x−y_max, x+y_max]` then a world at `x` takes no events at all, and for a
    barrier the interaction region is the barrier profile translated outward
    by exactly the reach — the finite-reach refinement of Theorem O1, and the
    sense in which the coherence horizon restores locality in position. Closes
    with Coulomb, whose Moyal series is geometric and converges iff the reach
    misses the nucleus.

    *Promoted to a specification.*
    [`../algorithm/compensated_liouville_algorithm.md`](../algorithm/compensated_liouville_algorithm.md)
    turns this note into an implementable open-line algorithm and records what
    the continuum argument leaves out: the reach fixes the momentum grid
    outright, the compensation must be taken against the kernel's own first
    moment rather than `V′(x)`, the Nyquist rung must be zeroed, and — an
    erratum for §4 here — the total-variation figures tabulated in this note
    are functions of the rung count rather than absolute numbers, because
    under a *hard* coherence horizon the event rate diverges logarithmically
    and the momentum churn linearly. Read that specification's §4.4 before
    quoting any event budget from §4 below.

15. **[`eckart_barrier_compensated.md`](eckart_barrier_compensated.md)** — The
    first open-line test problem that actually exercises the hop channel
    against a closed form. Theorem E7 rules out polynomials (no jump measure
    without a reach) and Theorem I4 rules out the harmonic and inverted
    harmonic alike (no jump channel at all), which between them had left the
    project without one. The Eckart pair barrier `V₀ sech²(r/a)` is bounded,
    asymptotically free on both sides, has a non-vanishing third derivative,
    and has an exact transmission coefficient. Theorem K1 generalises C6: the
    reach ceiling is the distance to the nearest *complex* singularity of
    `V`, so Coulomb's real-axis pole is the special case, sech² has the
    uniform ceiling `y_max < πa/2`, and for a soft core the softening length
    **is** the ceiling — with the corollary that `Δp > ℏ/a`, so fewer than
    `β` rungs span the barrier's own momentum scale and a reach-limited
    lattice cannot resolve the packet at all. Theorem K2 gives the exact
    far-field law, the hyperbolic continuation of Lemma C0, and corrects §5.1
    of the splitting note: for an exponential tail the reach *rescales* the
    interaction profile rather than translating it. Theorem K3: weighted by
    the potential's own spectrum the budget ratio saturates at 1 from below,
    so on the open line compensation never loses — the opposite of the ring,
    where §6.3 pins it at 1. Theorems K4–K6 are the point. The classical
    outcome functional is exactly invariant under streaming plus
    deterministic acceleration, so the entire quantum correction to the
    transmission is delivered by the residual channel (0.044475 measured
    against a closed-form 0.044134); that correction arrives as a small
    imbalance between two large opposed flows of positon–negaton pairs across
    the classical separatrix, net/gross = 0.19; and the cancellation tightens
    as `1/β`, because a packet centred on the barrier kills the slope term
    and leaves only the Jacobian `dp/dE` across the tunnelling window — so
    resolving `T` to fixed relative accuracy costs `β²` particles, and the
    semiclassical limit is the expensive one.

16. **[`sea_population_equilibrium.md`](sea_population_equilibrium.md)** — What
    the signed residual channel costs the *ledger*, as opposed to the
    observable. Prices open item CLA3 of the compensated specification. The
    Moyal equation fixes `u+ - u-` and says nothing about `u+ + u-`, so the
    population is extra structure; S1 shows any `E`-preserving sink is
    bilinear, which is exactly the class Proposition U1 left open. Momentum
    conservation forces the consumed sea pair to sit on the parent's own row
    (S0), making "ionisation" a derivation rather than a metaphor. The spec's
    rate `R = sum_q |K_res|` silently selects the **emissive** realisation of
    every event, and that choice is ruinous: the sea is relocated, not
    consumed, with the worst-cell deficit growing without bound (S4), and
    throttling by a finite sea moves `W` by 40 per cent in the core while
    leaving norm and `<p>` exact (S5). The **absorptive** realisation —
    a deposition of `-1` realised by removing a positon rather than adding a
    negaton — is identical in the observable and opposite in the ledger, so
    debit and credit both land on row `p`. It is per event and not per leg
    (S6, body momentum), hence supply limited. Theorem S7 is the note's
    centre: with absorptive fraction `f`, `dN = 2(1-2f) n_ev` and
    `dS = (2f-1) n_ev`, so `f = 1/2` closes both ledgers at once. Measured
    `f = 0.434`, and absorption already restores QLE fidelity by four orders
    at fixed `dt` — the emissive gap being the sign problem showing up as
    integration error rather than variance. Leaves the prediction that a
    supra-minimal initial ensemble drives `f` to `1/2`.

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
`src/demo_species_sectors_and_annihilation.py`,
`src/demo_compensated_liouville_splitting.py`,
`src/demo_compensated_liouville_algorithm.py`,
`src/demo_eckart_barrier_compensated.py`. All
non-trivial claims in these notes are verified numerically before inclusion.

## Related

Drawings of the elementary processes of the phase-alignment layer are in
[`../supplement/phase_alignment_interaction_diagrams.md`](../supplement/phase_alignment_interaction_diagrams.md); the simulation
counterparts of the ladder are in [`../algorithm/`](../algorithm/README.md).
