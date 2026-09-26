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

11b. **[`reach_energy_coupling.md`](reach_energy_coupling.md)** — What the
    reach actually controls, prompted by the question of whether the
    coherence horizon of `open_position_space.md` §3 is necessary at all.
    Settles open items 1 and 4 of that note and corrects three of its
    statements. Theorem E1: the reach is a *period*, not an aperture — a
    momentum lattice exists iff the difference field `D_x(y) = V(x+y) −
    V(x−y)` is periodic in the half separation `y`, with `Δp = πℏ/period`,
    and the ring, the periodicity of `V` and the postulated horizon are one
    mechanism with three sources for that period; a window of finite support
    does not by itself deliver a lattice. Proposition E1.1: a sharp window
    is exact on the lattice when `2y_max` is a whole number of periods of
    the potential and fails otherwise, so commensuration, not sharpness, is
    what matters, and Theorem E2 makes the fold of a decaying potential's
    tail back into one period exact. Theorem E3: signed world number and
    energy are reach-independent for *every* even window, both being moment
    conditions that follow from oddness — which closes open item 4, and
    contradicts an earlier draft of the note that expected the reach to
    couple to the energy balance once `V` has more than one harmonic.
    Theorem E4 extends the four-action energy ledger to any number of modes:
    for the symmetric member the focus channel does no net work and the hop
    channel delivers the whole classical power. The coupling that does exist
    lives one order higher. Theorem E5 finds a `1/y_max²` contamination of
    the third moment when a horizon profile is applied to the *full* kernel,
    and Theorem E6 shows it vanishes exactly when the profile is applied to
    the compensated residual instead. Theorem E7: for a polynomial `V` the
    Moyal series terminates, so there is no jump measure on the open line
    until a reach is imposed. Theorem E8: under the compensated split the
    residual event budget grows without bound with the reach for every `V`
    with `V′(x) ≠ 0` — linearly for bounded `V`, cubically for the quartic
    double well.

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
    semiclassical limit is the expensive one. Extended after the fact with §8.4, which
    prices the classical-trajectory ontology on this barrier: Proposition K8
    shows momentum conservation forces the emitted pair to be an *ionised*
    neutral pair co-located with the parent, and that `Gamma` is a horizon
    parameter rather than a property of the potential. Theorem K9 gives the
    supply condition `min(u+, u-) >= (dt/2)(|K| * N)` for the ledger to close
    with every trajectory Newtonian; it fails only where the minority species
    is exactly zero, which for a positive initial `W` is everywhere at `t=0`,
    so the ensemble must bootstrap. A standing dressed population `rho` of 3
    to 10 at reach `4 pi a` is the usable window, and in world-particles that
    is three to ten unpaired bodies against a sea of exactly two pairs per
    Planck cell.

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
    integration error rather than variance. S9 then shows `f` is an
    *attractor*: `f > 1/2` drains bodies and so removes partners, `f < 1/2`
    does the reverse, and a twentyfold padding of the initial ensemble and the
    minimal one converge on the same value from opposite sides. Closure
    therefore needs no tuning — and, for the same reason, says nothing about
    the one dimensionless constant the recombination rate still carries.

17. **[`compensated_ontology.md`](compensated_ontology.md)** — The ladder's
    picture stated as a proposed ontology, and bounded. Four postulates: (E)
    the world is a signed counting measure on phase space; (A) only ensembles
    whose expectation is a Wigner function of some `rho >= 0` occur; (S) every
    world-particle streams on a Newtonian arc under the *full* classical
    force; (D) pairs are ionised from and recombined into the sea at rate
    `Gamma = sum_q |K_res|`, absorptive fraction `f = 1/2`. Theorem G1 is only
    an assembly of C1–C3, K4, K8 and S7, and the note is really about what
    G1 leaves out. **Theorem G2**: for quadratic `V` the demographic channel
    is *empty*, measured at 1.2e-15 at every reach — while the published
    field-less signed-particle formulation has a rate of O(100) for the same
    system, which is the sharpest available argument that the compensated
    split is ontology and not numerics. So creation and annihilation cannot be
    the whole of the quantum, and **Theorem G3** locates the rest: the
    residual generator and the admissibility constraint are independent
    functions of `hbar`, and a quartic sweep closes the first continuously
    while the second does not move. Proposition G3.1 exhibits four Gaussians
    the dynamics cannot tell apart, differing by a factor eight in phase-space
    area; the inequality that separates them is exactly the Wigner bound
    `|W| <= 2/h`. **Theorem G4** was not anticipated and is the note's own
    contribution: admissibility is preserved by (S)+(D) and *not* by (S)
    alone, the least eigenvalue of the reconstructed `rho` going from the 1e-8
    grid floor to -0.10 under classical carrier transport — so the residual
    channel is doing kinematic work, and the sea is not optional. **Theorem
    G5** records the outstanding defect: `Gamma` grows roughly linearly with
    the coherence reach while the generator converges to 1e-14, so how many
    worlds exist is a property of the regulator. Corrects the framing that
    prompted it — the compensated split eliminates quantum *force*, not
    quantum *kinematics* — and leaves the split gauge-invariance theorem
    (G-SP1) as the load-bearing open item. §8 delimits the scope: `N` bodies
    verbatim, spin not at all, measurement untouched, and states narrowly what
    survives a survey of the signed-particle, Bell-type-QFT, Dirac-sea and
    negative-probability literatures now recorded in
    [`../../references/bibliography.md`](../../references/bibliography.md).

18. **[`stochastic_ledger.md`](stochastic_ledger.md)** — The ledger run as an
    exact Markov jump process on integer counts, rather than as a mean field
    on a mesh. Answers the standing request for a formal stochastic treatment
    of the free-body and sea-pair populations, and corrects the law it was
    meant to confirm. Theorem N1 strengthens J1: `P = S + N/2` is conserved
    *pathwise*, in exact integers, on every trajectory and not merely in
    expectation — a deterministic invariant of the generator, not a
    martingale. Theorem N2 is why the ledger is invisible: both realisations
    of an event move `E` by `+1` at the upper daughter and `-1` at the lower,
    so the Bernoulli(`f`) choice between them lies entirely in the kernel of
    the observable map, and the only noise `W` ever sees is Poisson
    event-timing noise common to both branches. The framing consequence is
    worth stating plainly — this is not a stochastic mechanics in Nelson's
    sense: trajectories are exactly Newtonian by (S), the diffusion in the
    phase variables is identically zero, and the whole stochastic content is
    demographic. **Theorem N3 is the note's centre and a correction to S7.**
    Stationarity of the body count gives `Gamma_tot (1 - 2f) = R_sink` with no
    closure and nothing about the potential, where `R_sink` is the total rate
    of every *other* body-removing channel — so `f = 1/2` is the sinkless
    special case and not the law, and any `kappa > 0` forces `f < 1/2` by a
    computable amount, verified to between 0.01 and 0.9 per cent over a
    fortyfold range of `kappa`. That invalidates the argument S-SP3 used to
    rule out a floor — that a systematic leak cannot change sign — because
    S-SP7 identifies a sink *and* a source, and their difference changes sign
    freely. Theorems N4 and N5 close the availability question: `f` is the
    probability that both legs find a partner, which under an
    independent-occupancy closure is `(1 - e^-lam)^2`, so `f = 1/2` pins `lam*
    = -ln(1 - 2^-1/2) = 1.227947` bodies per species per cell, `f'(lam*) =
    sqrt2 - 1` exactly, `Var(Lambda)/M = 1/f' = 2.414214`, and a Fano factor
    `0.983028`, slightly sub-Poissonian — five numbers carrying no free
    constant, measured at 2.4487 against 2.4559 in the well-mixed limit and
    flat across `Q = 1` to `12`, as a closure with no channel index in it
    requires. In world-particles that is 2.456 bodies against exactly two sea
    pairs per Planck cell, the low end of K9's empirical window. Theorem N6
    reassigns a role the project had given to recombination: with transport
    off, per-cell occupancy is a reflected critical random walk whose spread
    grows without bound, `kappa` damps it only partly and drags `f` to 0.32,
    and turning streaming on instead holds the spread flat and leaves `f` at
    0.5005. Streaming, not recombination, is the local regulator — and by N3
    no sink could have been, since every sink moves `f` off one half.

19. **[`soft_core_coulomb.md`](soft_core_coulomb.md)** — The potential open
    item K-LS5 nominated, put through the geometry the Eckart note put `V0
    sech^2(r/a)` through. Theorem Z1: for `-Z/sqrt(r^2 + eps^2)` the third
    derivative is `(Z/eps^4) u (6u^2 - 9) (1+u^2)^(-7/2)` with `u = r/eps`, so
    it vanishes at `r = 0` and `r = +- eps sqrt(3/2)` and the soft core
    carries **four** emission lobes — the same count K7 found for sech², but
    with a quiet radius set by the softening length alone and carrying no `Z`
    at all. Pure Coulomb has `V''' = -6Z/r^4`, one sign on each half-line and
    no interior zero, so it carries two: the inner pair of lobes is
    manufactured by the softening and collapses onto the origin as `eps -> 0`,
    which makes the soft core a source of structure rather than a
    regularisation convenience. Theorem Z2: `Gamma(0) = 0` identically at
    every reach, so the nucleus is dark exactly as the sech² summit is, and
    the interior quiet ring survives the horizon, drifting outward by 0.4 per
    cent at `y_max = eps/4` and 6.1 per cent at `y_max = 0.99 eps`. **Theorem
    Z3 is what this potential can raise and sech² cannot**: by Corollary K1.2
    the reach ceiling is `R(x) = sqrt(x^2 + eps^2)`, which *varies with
    position*, so either the momentum quantum `dp = pi hbar / 2 y_max` varies
    with `x` — and the phase-space crystal is not uniform — or a single
    uniform lattice must take the infimum `y_max < eps`, the tightest ceiling
    in the problem and located at the one point where nothing is emitted
    anyway. Theorem Z4 prices the uniform choice: the well is harmonic for
    `eps >> a0` with `sigma_r = eps^(3/4)/sqrt2`, so resolving the ground
    state to `k` rungs per `sigma_p` needs `eps >= k^4 pi^4 / 4` in units `a0
    = hbar^2 / mu Z` — 24.35 `a0` for one rung, 389.6 for two, the cost
    quartic in resolution; measured crossings 34.6 and 432.9, the ratio to
    prediction falling from 1.42 to 1.11 as the anharmonic correction dies. An
    unsoftened atom therefore cannot live on a uniform reach-limited crystal
    at all. Theorem Z5 is the consolation and the reason this is the right
    vehicle for S-SP6: at threshold the horizon spans 2.8 `sigma_r` and
    `Gamma(sigma_r)/Gamma_max = 0.974`, so unlike the Eckart barrier — where
    K-LS2's ceiling and resolution conditions have no common ground — the
    window where the lattice works is *not* the window where the demographic
    channel is empty. Leaves the non-uniform lattice as Z-LS1 and the
    transmission calculation, which has no closed form and needs a
    split-operator reference, as Z-LS2.

20. **[`dark_sea_and_worldline_identity.md`](dark_sea_and_worldline_identity.md)** —
    Whether the phase machinery transplants from the phase-resonance
    formulation onto the compensated one, and what the sea's darkness is a
    statement about. Theorem Y1 answers the sharpest test first, negatively
    for the no-go: the compensated kernel is read from the potential
    directly, `K(lam V) = lam K(V)` to machine precision, so the phase-blind
    no-go (Theorem 2) does not reach it, and the phase it needs is not
    absent: the potential acts on the density matrix only by winding its
    phase, the ladder `mu = arg rho`, and the Wigner transform turns that
    winding into the real, signed kernel `K`, while `arg rho` itself
    survives as the sign structure of `W`. By Corollary Y1.1 a
    particle-level phase can therefore never be *necessary* for `E`, and any
    phase rule lives in the kernel of the observable map — where, by Theorem
    N2, the realisation choice already lives. Theorem Y2 asks what the sea's
    ineligibility for recombination is worth, since all four ways of
    settling the two legs move `E` identically: with the sea eligible there
    are realisations with `dN = dS = 0` whose effect on the counts is
    exactly the excluded hop, and which conserve momentum by shifting one
    aligned pair a row, while `P = S + N/2` stays conserved throughout — so
    ineligibility, not arithmetic, is what gives postulate (S) and
    Proposition K8 their content, and with it every event carries `|dS| =
    1`. Theorems Y3 and Y4 supply the phase account darkness was missing:
    the two bodies a catalysed recombination consumes are a winding pair
    whose midpoint is the parent's own row, with `|Psi| = 2|sin(mu/2)|` and
    an envelope drifting at the parent's velocity, and a momentum kink of
    `+-xi` carries a phase ramp pinned at the point where it happens, so the
    created pair's misalignment freezes at `mu(x_k)` and is dark exactly
    when the kink sits at a node — a node-located kink and a phase reset at
    the parent's position being one freedom seen twice, not two rules.
    **Theorem Y5 is the cost.** If world-particles are conserved and none
    ever changes momentum at an event, per-row species counts are
    event-invariant and the residual channel can do nothing; so the
    two-clock reading forces piecewise worldlines, and the alternative is
    birth and death. Theorem Y6 prices that choice with tagged positons at
    the Eckart barrier: the individual momentum walk is driftless (`sum_q
    xi_q K_q = 0`) with diffusion `D_p = (1/2) sum_q xi_q^2 |K_q|`, an
    *unsigned* moment `E` never feels; `Gamma(0) = 0`, so no kink ever
    happens at the summit and an individual crosses by flank activation and
    classical passage, never through; and the measured crossing of tagged
    positons — 0.327 against 0.213 for `E` and 0.069 classically at `E_0 =
    V_0/2` — tracks neither the classical nor the quantum transmission,
    falling *below* both above the barrier, and is insensitive to the sea's
    momentum profile. Tunnelling therefore survives as a demographic account
    of `E` and not as a statement about identity, which corrects
    `eckart_barrier_compensated.md` §8.3 and the ledger note's "the
    diffusion in the phase variables is identically zero". Section 8 records
    a defect found along the way: the population clamp in
    `demo_emission_and_absorption.py` costs 16 per cent of `E` in 250 steps
    and biases `f`.

21. **[`fourd_compensated_ledger.md`](fourd_compensated_ledger.md)** — The
    compensated world form taken out of 1+1 dimensions for the first time: two
    particles on a line, the ledger run on exact integer counts with bodies
    carrying their own sub-cell positions, and the one cell size nothing fixes
    — the bin area `A` in the centre-of-mass directions — scanned. Proposition
    M1: for a pair potential the residual symbol depends on the ket–bra
    separation only through `y1 - y2`, so every event moves relative momentum
    and leaves `P` untouched, and Theorem M2 collapses step 10's leak law:
    under (S) every event is momentum-neutral in every direction and the
    field's share goes wholly through the force. Corollary M3 turns step 10's
    worst case into a free one — the harmonic sector carries no events, so
    Theorem D and the ring-seam obstruction do not arise — while Proposition M4
    records that the residual is not always cheaper than the uncompensated
    kernel on the same lattice (0.83 for a cosine, 2.8 for the attractive
    Poeschl–Teller well at the ledger runs' reach, as Theorem C4 of step 14
    predicts). **Proposition M5 is the note's centre**: the equilibrium body
    count is proportional to the number of occupied joint cells, 2.4 to 4.6
    bodies per cell, so shrinking `A` from infinity to `h` and `h/4` multiplies
    it by 4.8 and 13.2, with fidelity error growing as its square root.
    Proposition M6: the cells filled are those escaped bodies reach, not the
    state's support, and doubling the window raises the count again.
    Proposition M7 restricts N4 and N5: under real streaming, partner
    availability is structural, 0.55 to 0.74 whether the mean occupancy of a
    requested cell is 1.3 or 10.3. Proposition M8 answers step 10's open item 1
    for this layer: at `A` infinite, partners drawn from the wrong
    centre-of-mass component move the conditional potential energy by 4 within
    a quarter time unit, so co-location in the centre of mass is not a gauge
    choice. Proposition M9 finds the same trade-off inside the 1D cell, which
    the mesh ledger cannot show. Corollary M10: the cost is a product over
    degrees of freedom — the exponential of step 10's Propositions B2 and B3,
    returned as body inflation now that the sea is no longer a partner.
22. **[`sea_phase_reference.md`](sea_phase_reference.md)** — Whether the
    residual-kernel weights, which the compensated algorithm reads from the
    potential as a precomputed field, could instead emerge from the
    world-particles themselves, and what the sea would have to be for that to
    work. Theorem L1 writes each weight as a sum over contacts, `K_q(x) = -(B
    dp/hbar) int dy U_res(x,y) w(y) sin(2 xi_q y/hbar)`, in which the
    prefactor `B = 1/(pi hbar)` is the sea's own density and the strip the
    kernel sweeps holds exactly one aligned pair, `B dp 2 y_max = 1`; but by
    Proposition L2 sampling those contacts impulsively costs about `10^5`
    contacts per unit time to match the exact kernel, because each weight is a
    signed interference sum whose cancellations are the point, while replacing
    the pinned density `B` by the actual ledger sea costs almost nothing.
    Theorem L3 shows why per-world phases are the natural carrier: at their
    midpoint two bodies' extended phases wind apart at their energy
    difference, which for partners in one momentum row is the kernel's own
    `U(x,y)`. Theorem L4 is the obstacle it meets: with species signs, every
    aligned pair cancels exactly in any pair sum, so the sea can serve only as
    a set of unsigned clocks — a phase reference, invisible as amplitude.
    Definition L0 makes the required relation precise: the **sea lock**, each
    aligned pair's clock sampling its row's plane wave, which classical flow
    carries as a smooth phase field `theta = S/hbar` over phase space that
    folds only when projected onto position. A gas of random phases never
    locks under the event rules; a locked sea largely persists. Theorem L5
    introduces **dark catalysis**, an interaction between aligned pairs at the
    kernel's own rate in which both members of a pair fire together and their
    deposits cancel exactly, so it moves phase and nothing else; by
    Proposition L6 it couples nearly the whole sea near the barrier, makes
    created pairs dark, and at reach scale maintains the lock against the
    damage events do — 0.470 against 0.461 with no events at all and 0.395
    with co-located coupling. **Theorem L7 is the negative result.** A freshly
    locked sea dephases at a rate that is a force-free pair sum linear in `V`,
    reproduced to correlation 0.999, but the force it subtracts is the
    partners' rather than the parent's: it is the trapezoid rule's error for
    `U = int V' dx'` over the chord where the compensated kernel is the
    midpoint rule's, so its third moment — the leading quantum correction — is
    about `-2` times the kernel's. **Theorem L8 removes the obstacle.**
    Comparing two clocks at a distance needs a lever momentum, and whoever
    supplies it is the reader: an inertial reader reads the full `U`, each
    partner's own lever reads the trapezoid residual, and a reader streaming
    with the parent under (S) reads `U_res` exactly at re-lock, to all orders
    in `y`. The compensated split is therefore the choice of the parent as
    reader, the sea at the vertex is a set of clocks at positions, the
    leg-force term is how fast the lever mismatch grows, not a postulate, and
    the reader's lever is P2 applied to the reader, not a new rule.
    **Proposition L9** caps the drift of that reading by the row for the
    partners a vertex actually reads, so in the particle model the limit is
    the lock itself: more than half of the ceiling near 0.46 came from the
    continuous momenta within each row of the initial lock, and a sea that
    keeps its rows, re-locked at the kernel's own rate, reads 0.82 against a
    control of 0.89 — but keeping its rows means its momenta ignore the force,
    which (S) does not allow. Sections 5, 8 and 9 carry the figures, and §5 an
    interactive page of the phase field.
23. **[`force_blind_sea.md`](force_blind_sea.md)** — Postulate (S) gives every
    world-particle the full classical force, the members of an aligned sea
    pair included, and nothing in the QLE asks for that: a pair is `W`-null,
    and the Moyal equation fixes `u+ - u-` and says nothing about `u+ + u-`.
    This note provisionally adopts **(S′)**: free bodies obey the full force,
    aligned pairs move inertially — advecting at `p/m` with no drift in `p` —
    while their clocks still wind with `V`, so the sea is force-blind but
    potential-sensitive; a motionless pair would need a Hamiltonian clock, and
    a potential-blind one would leave the kernel nothing to be read from.
    **Proposition Q1** organises everything: in every ledger of the chain the
    sea enters an event only as its source or its sink and is never read, so
    the body fields — and with them `E`, `N` and `f` — are independent of how
    it moves, verified bitwise; the exceptions are exactly where the sea is
    read, a throttled rate, a sea-weighted kernel, its phases and identity.
    **Proposition Q2** gives the kinematics: under (S′) a sea clock stays on
    its row's plane wave up to the eikonal phase and a same-row pair winds at
    exactly `U`, so Theorem L8's drift term vanishes; under (S) `d(hbar theta
    - p x) = -H dt - x dp`, a uniform force keeps the row lock, and the
    mismatch begins at `V''`, which is the row mixing step 22 measured.
    **Proposition Q3**: what changes is the sea itself — deficits stay in
    their rows, Theorem S8's early dip goes but the worst cell hovers near
    zero, and fast recombination repairs it where slow does not. **Proposition
    Q4**: with step 20's tags now conserved, tagged bodies cross the Eckart
    barrier inside dark pairs, through the summit — `T_tag` 0.89 against 0.34
    at `E0 = V0/2`, with only 9 per cent of the tag above the barrier energy —
    while `T_E` is unchanged, so individual crossing is still not tunnelling.
    Two demo defects are repaired: a species mask read at the destination row,
    the whole of step 22's drift of `Sum E` (L-SP8), and the tag bookkeeping
    of step 20's Part E. Section 8 lists the corrections to steps 15, 16, 17,
    20 and 22 that (S′) would require.

## Index

[`INDEX.md`](INDEX.md) lists every labelled theorem, proposition, lemma,
corollary, definition, postulate and open item in these notes, each with a
section reference, a one-line statement and a *standing* — whether a later
note corrected, restricted, superseded or retracted it, or which open item it
feeds. It also carries the prefix registry (which theorem letters are
reused), a ledger of what each note corrected, and the open items with their
equivalents in the algorithm and supplement folders. Use the ladder above to
read; use the index to look a label up.

The index is updated in the same patch as the note it indexes, alongside
that note's ladder entry above and the abstract blockquote under its title.
Both are checked:

```bash
PYTHONPATH=src python3 -m wpmwlib.check_abstracts docs/analysis
PYTHONPATH=src python3 -m wpmwlib.check_index docs/analysis
```

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
`src/demo_reach_energy_coupling.py`,
`src/demo_interworld_coupling.py`,
`src/demo_species_sectors_and_annihilation.py`,
`src/demo_compensated_liouville_splitting.py`,
`src/demo_compensated_liouville_algorithm.py`,
`src/demo_eckart_barrier_compensated.py`,
`src/demo_compensated_ontology.py`,
`src/demo_stochastic_ledger.py`,
`src/demo_soft_core_coulomb.py`. All
non-trivial claims in these notes are verified numerically before inclusion.

## Related

Drawings of the elementary processes of the phase-alignment layer are in
[`../supplement/phase_alignment_interaction_diagrams.md`](../supplement/phase_alignment_interaction_diagrams.md); the simulation
counterparts of the ladder are in [`../algorithm/`](../algorithm/README.md).
