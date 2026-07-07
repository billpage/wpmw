# Sea-Dressed Two-Body Microdynamics: A Fully World-Particle Realization of the Crystal-Lattice Collision Term

## 0. Status and provenance

This note is the sequel to
[`four_rule_microdynamics_equivalence.md`](four_rule_microdynamics_equivalence.md).
That analysis established that the four-rule scheme (Focus, Defocus,
Right-Hop, Left-Hop) reproduces the quantum Liouville equation (QLE) exactly,
but ended on a no-go lemma: rate laws built from pairwise mass action among
the tracked particles are quadratic in occupancy, while the QLE generator is
linear, so a fully collision-based microdynamics is impossible **unless the
rates are linearized by species whose density is pinned** (a reservoir).

The present note takes the step that lemma leaves open. The proposal
(B. Page, July 2026 conversation) is that the reservoir is not an external
bookkeeping device but the model's own Dirac sea: the uniform background
$W' = W + 2/h$ of the crystal shift, reinterpreted as a dense population of
bound positon–negaton pairs, whose hidden microdynamics (dissociation,
recombination, polarization) supplies the pairwise interactions that the
excess particles alone cannot. The result is a rate table of **sixteen
collision channels**, each a local, momentum-conserving, two-body event
between world-particles, whose mean field at pinned sea is *exactly* the
symmetric member of the four-rule family, and hence exactly the QLE.

Derivation and implementation were developed jointly with Claude (Anthropic),
July 2026. Companion code: `src/wpmwlib/sea_dressed_lattice.py` and
`src/demo_sea_dressed_dynamics.py`. All numerical claims in §10 are outputs
of that demo.

## 1. The question, and where the answer must live

The four-rule scheme removed the mediated rule's action-at-a-distance in
momentum, but its rate laws still consult two things no collision participant
carries: the mean-field mode amplitude $\Gamma_q(x)$, and (in the hop channel)
the occupancy of cells the hopping particle never touches. The no-go lemma
says these consultations cannot be eliminated among the *excess* particles
alone. Schematically:

```math
\text{mass action: rate} \propto N_a N_b \quad (\text{quadratic}), \qquad
\text{QLE generator: rate} \propto N_a \quad (\text{linear}).
```

The loophole is standard chemical kinetics: a bilinear rate $k \cdot B \cdot N$
is *linear in $N$* if the partner species $B$ is pinned at constant density
(a buffered reagent, a photon bath, a catalyst). The crystal lattice already
owns a canonical pinned density: the background $2/h$, which the algorithm
spec introduces as a passive shift. The move made here is to promote it to
an active medium:

- **ground sea pairs** at density $B$ per cell — the crystal shift, embodied;
- **polarized sea excitations** — coherent pair excitations carrying momentum
  $\pm q h/L$, whose local density realizes $\Gamma_q(x)$ as a *prepared state of
  the medium* rather than a consulted field;
- **unpaired positons and negatons** — the dynamical excess, including
  Wigner negativity as a literal negaton surplus.

Mediation is thereby not eliminated but **embodied**: every channel's rate
constant carries one factor of the local polarization amplitude. What was a
field consulted by a rule becomes a population collided with by a particle.
Correspondingly, the hop channel's momentum-conservation problem dissolves:
the sea partner or the exchanged excitation takes the recoil, and *every*
event in the table below is a momentum-conserving two-body collision.

## 2. Ontology: two species and a pairing relation

The model has exactly two fundamental species — **positons** and
**negatons** — plus a two-place pairing relation. The three "species" of the
implementation are *states*, not new particles:

| implementation species | ontological reading | contributes to $E$ |
|---|---|---|
| $U^+$ unpaired positon | positon not currently bound | $+1$ |
| $U^-$ unpaired negaton ("hole") | negaton not currently bound | $-1$ |
| $S$ ground sea pair | bound positon–negaton pair at one cell | $0$ |
| (polarized excitation) | excited bound pair carrying $\pm q h/L$ | $0$ |

The observable excess field is $E = U^+ - U^-$, and the Wigner distribution is
$W = E / (nu \cdot dx \cdot dp)$. Wigner negativity is a local negaton surplus —
$W < 0$ cells are cells where unpaired negatons outnumber unpaired positons.

Three consequences worth stating before the table:

1. **No worldline ever begins or ends.** Every creation/annihilation in the
   model is of *composites and states* (pairs form and dissociate, excitations
   are exchanged); the fundamental particles persist. §7 makes this a checked
   invariant.
2. **Crossing symmetry is built in.** Every channel acting on unpaired
   positons has a conjugate acting on unpaired negatons with the mirrored
   effect on $E$. The *signed* rates of the four-rule scheme decompose into
   pairs of individually non-negative channel rates; the sign of $E$ never
   enters a rate law.
3. **The sea is a non-equilibrium reservoir.** The polarization pattern
   $\Gamma_q(x)$ is a prepared, pumped state (the analogy is a laser medium,
   not a thermal bath). Detailed balance is deliberately broken: each channel
   runs one way, with the direction set by the local polarization sign. This
   is why linear one-way rates are legitimate: the microreversibility that
   would force quadratic reverse rates is pre-empted by the preparation of
   the medium (see §6).

## 3. State variables and reservoirs

Per lattice cell $(n, m)$ (momentum index $n$, position index $m$):

```math
U^+_{n,m} \in \mathbb{Z}_{\ge 0}, \qquad
U^-_{n,m} \in \mathbb{Z}_{\ge 0}, \qquad
S_{n,m} \in \mathbb{Z}_{\ge 0}, \qquad
E_{n,m} = U^+_{n,m} - U^-_{n,m}.
```

The sea is initialized virgin: $S = B := \mathrm{round}(nu \cdot (2/h) \cdot dx \cdot dp)$ in every
cell — the exact integer image of the crystal shift (25 000 pairs per cell at
the demo's parameters). The polarized-excitation populations are held pinned
at level 2 (below), entering only through the rate constants; §12 discusses
promoting them to tracked reservoirs.

Three levels of idealization organize what follows:

- **Level 1 (pinned sea).** $S/B \equiv 1$ in all rates, and recombination is
  instantaneous ($\kappa \to \infty$). This is the regime of the exactness theorem (§5).
- **Level 2 (live sea).** $S$ is a real ledger debited by every channel and
  credited by a finite-rate recombination channel. Pinning becomes a
  *derived* property — a timescale separation — rather than an axiom (§8).
- **Level 3 (self-consistent sea).** The polarization itself is dynamical,
  imprinted on the sea by the sources of $V(x)$. Open (§12).

## 4. The channel table

Fix a mode $q$ and a position column $x$ with signed rate field

```math
\Gamma_q(x) = -\frac{V_q}{\hbar}\,\sin\!\Big(\frac{2\pi q x}{L} + \phi_q\Big),
\qquad \sigma = \mathrm{sign}(\Gamma_q(x)), \qquad \gamma = |\Gamma_q(x)|.
```

Define the roles **hi** $= n + \sigma q$ and **lo** $= n - \sigma q$ relative
to a center $n$. (For $\sigma = -1$ the two blocks below are the mirror
images $q \to -q$; at most one block is active per column per mode.) Write
$s_c = S_c / B$ for the sea factor at cell $c$ ($s_c \equiv 1$ at level 1).

Eight channels per polarization block. Rates are per center $n$ per unit
time; effects list the integer updates per event.

| ch. | rate | effect on $U^+ / U^- / S$ | $\Delta E$ stencil | collision reading |
|---|---|---|---|---|
| **K1** | $`(\gamma /2) s_{lo} U^+_{hi}`$ | $`U^+_{hi}-1, U^+_n+2, U^-_{lo}+1, S_{lo}-1`$ | focus $(-1,+2,-1)$ | unpaired positon at *hi* and sea positon at *lo* scatter into $n$; the sea partner's negaton is orphaned at *lo* |
| **K1b** | $`(\gamma /2) s_{lo} U^-_{hi}`$ | $`U^-_{hi}-1, U^-_n+2, U^+_{lo}+1, S_{lo}-1`$ | defocus $(+1,-2,+1)$ | crossing conjugate of K1 (negaton capture) |
| **K2** | $`(\gamma /2) s_n^2 U^+_{lo}`$ | $`U^+_{hi}+1, U^+_{lo}+1, U^-_n+2, S_n-2`$ | defocus | two sea positons at $n$ scatter to *hi*/*lo*, Bose-stimulated by the unpaired-positon occupancy at *lo*; both negatons orphaned at $n$ |
| **K2b** | $`(\gamma /2) s_n^2 U^-_{lo}`$ | $`U^-_{hi}+1, U^-_{lo}+1, U^+_n+2, S_n-2`$ | focus | crossing conjugate of K2 |
| **K3** | $`(\gamma /2) U^+_{hi}`$ | $`U^+_{hi}-1, U^+_{lo}+1`$ | hop *hi*→*lo* $(-1,0,+1)$ | unpaired positon absorbs a polarized excitation carrying $`-\sigma \cdot 2q\,dp`$ |
| **K3b** | $`(\gamma /2) U^-_{hi}`$ | $`U^-_{hi}-1, U^-_{lo}+1`$ | hop *lo*→*hi* on $E$ | crossing conjugate of K3 |
| **K4** | $`(\gamma /2) s_{hi} U^+_{lo}`$ | $`U^+_{lo}+1, U^-_{hi}+1, S_{hi}-1`$ | hop *hi*→*lo* | a sea positon at *hi* emits an excitation and lands at *lo*, destination-stimulated by $`U^+_{lo}`$; its negaton is orphaned at *hi* |
| **K4b** | $`(\gamma /2) s_{hi} U^-_{lo}`$ | $`U^-_{lo}+1, U^+_{hi}+1, S_{hi}-1`$ | hop *lo*→*hi* on $E$ | crossing conjugate of K4 |

**Momentum bookkeeping (every row is a genuine two-body event):**

- K1/K1b: the two captured particles change momentum by $`-\sigma q\,dp`$ and
  $`+\sigma q\,dp`$; total $0$.
- K2/K2b: the two sea partners recoil to $`\pm \sigma q\,dp`$; total $0$.
- K3/K3b: the particle's $`-\sigma \cdot 2q\,dp`$ is supplied by the absorbed excitation.
- K4/K4b: the particle's $`-\sigma \cdot 2q\,dp`$ is carried off by the emitted excitation.

**Recombination channel** (the hidden sea microdynamics, level 2):

| ch. | rate (per cell) | effect | $\Delta E$ |
|---|---|---|---|
| **R** | $`\kappa \cdot U^+_n U^-_n / B`$ | $`U^+_n-1, U^-_n-1, S_n+1`$ | $0$ |

R is the only quadratic rate in the model, and it is allowed to be: it lies
entirely outside the QLE generator (it changes $E$ by exactly zero), so the
linearity constraint never touches it. Its role is stability, not dynamics
(§8).

## 5. Exactness at pinned sea

At level 1 ($s \equiv 1$), sum the signed $\Delta E$ contributions of the block. Per
center $n$ (writing the $\sigma$-block quantities and recalling
$E = U^+ - U^-$):

```math
f_n \;=\; \underbrace{\tfrac{\gamma}{2}U^+_{hi}}_{K1}
        - \underbrace{\tfrac{\gamma}{2}U^-_{hi}}_{K1b}
        - \underbrace{\tfrac{\gamma}{2}U^+_{lo}}_{K2}
        + \underbrace{\tfrac{\gamma}{2}U^-_{lo}}_{K2b}
 \;=\; \tfrac{\gamma}{2}\,\big(E_{hi} - E_{lo}\big),
```

```math
h^{hi\to lo}_n \;=\; \underbrace{\tfrac{\gamma}{2}U^+_{hi}}_{K3}
        - \underbrace{\tfrac{\gamma}{2}U^-_{hi}}_{K3b}
        + \underbrace{\tfrac{\gamma}{2}U^+_{lo}}_{K4}
        - \underbrace{\tfrac{\gamma}{2}U^-_{lo}}_{K4b}
 \;=\; \tfrac{\gamma}{2}\,\big(E_{hi} + E_{lo}\big).
```

Restoring the sign conventions of the four-rule note
($hi = n+q$, $lo = n-q$ when $\Gamma > 0$; a *hi*→*lo* hop is a left-hop,
i.e. negative $h_n$):

```math
f_n = \frac{\Gamma}{2}\big(E_{n+q} - E_{n-q}\big), \qquad
h_n = -\frac{\Gamma}{2}\big(E_{n+q} + E_{n-q}\big),
```

which is **precisely member (5), the symmetric member**, of the exact-rate
family of the four-rule analysis — and therefore, by the theorem proved
there, the sixteen-channel generator equals the single-rule generator equals
the discrete QLE stencil, identically in the state and at every lattice
resolution. Note the derivation never assumed $E \ge 0$: the crossing
conjugates carry the sign of $E$ channel-by-channel, with every individual
rate non-negative. §10 (Part A) verifies the identity numerically on
independent random $U^+$, $U^-$ fields to machine precision.

**Why the symmetric member, and only it, shows up.** Call a rate assignment
*participant-local* if every channel's rate is proportional to the occupancy
of a cell that the channel's event actually touches (a source or a
destination). The focus/defocus event touches $n$ and $n \pm q$, so a
participant-local focus bias must have the operator form
$F = a_+ A + a_- A^{-1} + a_0 I$ (in the shift-operator notation of the
four-rule note). Membership in the exact family requires
$F = (A - A^{-1})G$, whose Fourier symbol vanishes at $`q\theta \in \{0, \pi\}`$; imposing that on the symbol of $F$ forces $a_0 = 0$ and
$a_- = -a_+$. The hop event touches only $n \pm q$ — the hopped-over cell
$n$ is *not* a participant — so $H$ may contain no $I$-component beyond null
biases; the family formula $H = (2 - A - A^{-1})G - \Gamma I + H_0$ then
forces $`G = (\Gamma/2)\,I`$ up to the usual $\ker(A^{-1}-A)$ null terms. That
is member (5). In other words: **demanding that the sea-dressed channels be
local collisions among their own participants uniquely selects the symmetric
member** — the same member that minimizes gross traffic. The two-body
ontology and the variance-optimal bookkeeping point at the same rates, which
is a satisfying coincidence and possibly not a coincidence.

## 6. Collision readings, stimulated channels, and broken detailed balance

The channels divide into two classes.

**Clean mass action (K1, K1b, K3, K3b).** K1 is an ordinary bimolecular
capture: rate $\propto$ (density of unpaired positons at *hi*) $\times$ (density of sea
partners at *lo*), with the bilinearity rendered linear by the pinned sea
factor — the no-go lemma's loophole, used exactly as intended. K3 is
absorption from the polarized-excitation reservoir: rate $\propto$ (reservoir
density) $\times$ (absorber density), with the reservoir density $\propto \gamma$ pinned.
These channels are Einstein absorption kinetics with the sea polarization in
the role of the photon bath.

**Destination-stimulated (K2, K2b, K4, K4b).** K4's reading is spontaneous
emission by a *sea* positon, Bose-stimulated by the occupancy of the
destination cell: physical rate $`\propto B(1 + U^+_{lo}/\mathcal{N})`$ with
$\mathcal{N}$ the micro-mode multiplicity of a cell. The constant part is a
uniform hop rate — a **null bias** ($H_0$) with no mean-field effect — and
the linear part is exactly K4's table rate. Because the *emitter* is drawn
from the pinned sea, no second dynamic occupancy enters and the rate is
exactly linear: this is the trick that makes the stimulated channels honest.
K2 is the same structure for a sea–sea scattering event with the products at
$n \pm q$.

One residual idealization must be flagged. Full Bose kinetics would
stimulate K2 by *both* destinations, $\propto (1 + U_{lo}/\mathcal{N})(1 + U_{hi}/\mathcal{N})$, whose expansion contains a bilinear term
$U_{lo}U_{hi}/\mathcal{N}$. The constant term is null (a uniform defocus rate
telescopes to zero), the linear terms are the table's (a *Bose-symmetric
variant* of the table reweights K1 to rate $`\gamma s_{lo} U^+_{hi}`$ and K2 to
$`(\gamma /2)(U^+_{lo} + U^+_{hi})`$; its mean field is identical), and the bilinear
term is suppressed as $O(1/\mathcal{N})$, vanishing in the fine-mode-structure
limit $\mathcal{N} \to \infty$ — a limit exactly parallel to the reservoir
limit $B \to \infty$. At the level of this note the process is *defined* by
the linear rate laws; the $1/\mathcal{N}$ bookkeeping is recorded as an open
item (§12).

**Broken detailed balance is a feature, not an omission.** Each channel runs
one way; the reverse processes (e.g. the microreverse of K1) are not in the
table. In an equilibrium medium this would be inconsistent — microreversible
mass action would force quadratic reverse rates and resurrect the no-go
lemma. But the sea is *prepared*: its polarization is a pumped,
non-equilibrium state, and one-way stimulated kinetics in a pumped medium is
ordinary physics (gain in a laser). The arrow supplied by the preparation is
exactly the arrow $\Gamma_q(x)$ used to supply by fiat.

## 7. Worldline continuity

Inspect the effect columns of §4: every channel conserves both

```math
\mathcal{P} = \sum_{n,m}\big(U^+_{n,m} + S_{n,m}\big)
\quad\text{(total positons)}, \qquad
\mathcal{M} = \sum_{n,m}\big(U^-_{n,m} + S_{n,m}\big)
\quad\text{(total negatons)},
```

and so does R ($-1, -1, +1$ in $U^+, U^-, S$) and free streaming (a
permutation of cells). No fundamental worldline is ever created or
destroyed: what the occupancy picture records as birth and death of quanta
is, underneath, pairing changes and momentum exchange among persistent
particles. The demo asserts $\mathcal{P}$ and $\mathcal{M}$ exactly at every
run's end (§10). This is the substrate on which the labeled ("tagged
world-particle") unraveling of the model can be built — each particle's
history is piecewise free-streaming punctuated by discrete collisions — with
the identity assignment at collision events a gauge choice invisible to all
occupancy observables. That construction and its non-uniqueness are
discussed in the July 2026 conversation record and deferred here.

## 8. Back-reaction: the sea as a live ledger, and why recombination is not optional

At level 2 every capture, split, and emission debits the pair ledger $S$ and
credits the orphan populations. Two coupled degradation mechanisms follow:

1. **Depletion.** Rates carry $s = S/B$; as $S$ drains where the dynamics is
   busy, the effective $\Gamma$ weakens *inhomogeneously* and the evolution
   falls behind the QLE. (The channels drain the sea asymmetrically — K2
   fastest through $s^2$, K3 not at all — so depletion also *distorts*, not
   merely slows.)
2. **Orphan avalanche.** Six of the eight channels orphan a partner. The
   crossing-conjugate pairs (K1, K1b) etc. fire at rates proportional to
   $U^+$ and $U^-$ *separately* — only their means cancel. Orphans therefore
   beget events which beget orphans: without recombination the unpaired
   population, and with it the gross event rate and the shot noise, grows
   exponentially at rate $\sim |\Gamma|$. This was observed directly during
   development: an early "pinned" variant without instantaneous re-pairing
   destroyed the state (relative L2 error > 600) within two periods.

The recombination channel R is therefore not decorative bookkeeping but the
**stabilizer of the entire construction**: it drains the orphan load and
restores the ledger, at exactly zero cost to the observable dynamics
($\Delta E = 0$ per event). Level 1 is recovered as the $\kappa \to \infty$ limit — pinning is
the statement that recombination is fast compared to $|\Gamma|$, i.e. a
timescale separation, not an axiom. §10 (Part B) shows the convergence
$\kappa = 0 \to 20 \to 200 \to$ pinned quantitatively.

## 9. Algorithm

Implemented in `src/wpmwlib/sea_dressed_lattice.py`
(class `SeaDressedLattice`, extending `PhaseSpaceCrystalLattice`):

1. **Streaming.** Integer-roll advection applied to $U^+$, $U^-$, and $S$
   alike (spec §3a). A uniform sea is streaming-invariant; depleted spots
   ride their momentum row.
2. **Jump substep** (`step_jump_sea_mc`). Per mode and polarization block,
   the eight channels are tau-leaped in the fixed order K1, K1b, K2, K2b,
   K3, K3b, K4, K4b: draw $e \sim \mathrm{Poisson}(rate \cdot dt)$ lattice-wide, cap by the
   drained populations, apply. Within any single channel each drained cell is
   drained from exactly one center role, so the caps are safe under fully
   vectorized application; between channels the counts update sequentially.
   In pinned mode the step ends with instantaneous re-pairing
   $r = \min(U^+, U^-)$ per cell.
3. **Recombination substep** (`step_recombine`). $e \sim \mathrm{Poisson}(\kappa U^+ U^- dt / B)$ per cell, capped by $\min(U^+, U^-)$.
4. **Observables.** $`W = (U^+ - U^-)/(nu\,dx\,dp)`$; diagnostics
   $\min_c S_c/B$ (worst-cell depletion), $\sum(U^+ + U^-)$ (orphan load),
   and the worldline invariants of §7.

Cost: the sixteen channels are $O(NM)$ vectorized passes per mode per
substep; the full 64×64, 41-macro-step, 16-substep demo runs all four
sea-dressed solvers plus references in ~20 s (NumPy, single core).

## 10. Numerical verification

`src/demo_sea_dressed_dynamics.py`, seed 20260706. Parameters identical to
the four-rule demo Part B: cosine well $V = -V_p\cos(2\pi x/L)$, $L = 8$,
$V_p = 1.5$, 64×64 lattice, $nu = 1.6\times10^6$ (so $`B = 25\,000`$ pairs per
cell), squeezed Gaussian initial state, $`dt_{adv} = m\,dx/dp`$ (exact integer
advection), 16 jump substeps per macro step
($`|\Gamma|_{max}\,dt_{jump} = 0.030`$), duration two well periods.

**Part A — generator identity.** The sixteen-channel mean-field generator,
assembled channel by channel on *independent* random fields $U^+$ and $U^-$
(so the crossing structure is genuinely exercised, not hidden by $U^- = 0$),
against the single-rule QLE stencil acting on $E = U^+ - U^-$; mode sets
$q=1$, $q=2$, and $`q=\{1,2,3\}`$ mixed, three random trials each:

```
worst absolute deviation: 5.7e-17   (relative ~ 1e-15)
```

Machine precision — the table of §4 is an exact regrouping of the QLE.

**Part B — stochastic evolution with a live sea.** Six solvers, identical
splitting sequence: mesh QLE (target), four-rule MC (reference), sea-dressed
pinned, and live sea at $`\kappa \in \{0, 20, 200\}`$. Relative L2 deviation from the
mesh at $t = 2T$:

| solver | relL2 at 2T | final $\min S/B$ | final orphan load |
|---|---|---|---|
| four-rule MC (reference) | 0.074 | — | — |
| sea MC, pinned (level 1) | 0.108 | (frozen) | 1.8 M |
| sea MC, live, $\kappa = 200$ | 0.172 | 0.945 | 5.1 M |
| sea MC, live, $\kappa = 20$  | 0.352 | 0.796 | 29 M |
| sea MC, live, $\kappa = 0$   | 0.764 | 0.000 | 206 M |

Readings, in order of §8's predictions:

- **Pinned tracks the QLE at the shot-noise floor.** The 1.45× noise ratio
  to the four-rule MC is the price of drawing the crossing-conjugate
  channels independently (gross traffic) where the four-rule MC draws a
  single signed net — the same gross-versus-net tradeoff identified in the
  four-rule note, now appearing one level down.
- **$\kappa = 0$ exhibits both failure modes at once**: the worst cell drains the
  sea *completely* ($\min S/B = 0$), and the orphan load runs away to 206 M
  unpaired particles (an 8× multiple of the entire sea), washing out the
  interference structure (relL2 0.76).
- **Convergence in $\kappa$.** $0.764 \to 0.352 \to 0.172$ approaching the pinned
  0.108, with the sea held at 95% everywhere by $`\kappa = 200 \approx 130\,|\Gamma|_{max}`$. Pinning is visibly a derived, quantitative property
  of fast recombination.
- **Worldline invariants** $(\mathcal{P}, \mathcal{M})$ conserved exactly
  (integer equality) across all live runs — ~10⁸ collision events with no
  particle created or destroyed.

Figures (committed to the `output` branch by the demo when `WPMW_DOCS` is
set):

![evolution](https://raw.githubusercontent.com/billpage/wpmw/output/figures/sea_dressed_evolution.png)

![metrics](https://raw.githubusercontent.com/billpage/wpmw/output/figures/sea_dressed_metrics.png)

## 11. What has and has not been achieved

Achieved:

- The QLE collision term is generated, exactly, by individually non-negative,
  local, momentum-conserving, two-body collision channels among
  world-particles — no rate law consults the occupancy of a non-participant,
  and no event violates momentum conservation. The no-go lemma is satisfied,
  not evaded: its reservoir loophole is instantiated by the model's own
  background.
- The reservoir idealization is demoted from axiom to limit: level 2 makes
  the sea a live ledger and exhibits pinning as fast recombination, with the
  approach to the ideal measured.
- Fundamental-particle worldlines are exactly conserved, providing the
  substrate for trajectory/identity unravelings.

Not achieved, and stated plainly:

- **Mediation is embodied, not eliminated.** The polarization pattern that
  sets every rate constant is *prepared* to match $\Gamma_q(x)$; at levels 1
  and 2 the model does not explain the preparation, it assumes it — exactly
  as the original algorithm assumes $V(x)$. The gain is ontological
  (a consulted field becomes a collided-with population) and structural
  (momentum conservation, two-body locality), not a derivation of the
  potential.
- **The stimulated channels carry a residual $O(1/\mathcal{N})$
  idealization** (§6): exact linearity of the Bose enhancement is imposed
  rather than derived.
- **The excitation ledger is untracked**: K3 consumes and K4 emits
  polarization quanta of opposite sign; at level 2 the reservoir is pinned
  against this traffic.

## 12. Open items

1. **Level 3 (self-consistent sea).** Let the sources of $V(x)$ imprint the
   polarization dynamically — a wave equation or relaxation dynamics on the
   sea's excitation populations, with $\Gamma_q(x)$ emerging as its steady
   state. This would close the loop the four-rule note's "mediator" language
   gestures at, and would make the excitation ledger (item 3) mandatory.
2. **Quantify the Bose bilinear.** Implement the double-stimulated K2 with
   finite $\mathcal{N}$ and measure the departure from the QLE as
   $1/\mathcal{N}$ scaling; identify $\mathcal{N}$ physically (micro-modes
   per cell vs. $B$).
3. **Excitation momentum ledger.** Track $\rho_{\pm 2q}(x)$ explicitly;
   check whether K3/K4 traffic self-consistently maintains or degrades the
   polarization, and whether R has an excitation-emitting analogue.
4. **Variance structure.** The pinned run pays 1.45× over signed-net
   four-rule sampling. Explore variance reduction that preserves the
   two-body ontology (correlated draws of crossing conjugates; the
   Bose-symmetric variant of §6; the gauge freedom $G$ of the exact family).
5. **Streaming conventions for composites.** Pairs here ride the momentum
   row of their cell; alternatives (static sea frame, pair
   center-of-mass kinematics) change nothing at uniform $S$ but matter once
   depletion structures the ledger.
6. **Canonical identity unraveling.** With worldline conservation
   established, implement a tagged-particle demo (e.g. minimal-transport
   selection at collision events) and study the geometry of the resulting
   trajectory ensemble against Bohm and Nelson trajectories for the same
   states.
