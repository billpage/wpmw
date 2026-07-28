# Phase-Alignment Microdynamics Algorithm

**A specification for simulating the quantum Liouville equation from a live sea of world-particle pairs, with the misalignment of transported clock phases as the only relational state.**

---

## 0. Status of this specification

This document specifies the simulation counterpart of
`docs/analysis/phase_alignment_microdynamics.md`. Where
`docs/algorithm/phase_space_crystal_lattice_algorithm.md` takes the rate
field $\Gamma_q(x)$ as given and applies it as a stencil, and where
`docs/analysis/sea_dressed_microdynamics.md` realises that stencil as
channels against a *pinned* sea, this specification describes a **live
sea**: the pairs are dynamical, they carry phases, and $\Gamma_q(x)$ is
not an input but an output, assembled from vertex statistics.

The specification is therefore **more expensive and less accurate** than
the mesh algorithm for any problem the mesh algorithm can solve. Its
purpose is different: it is an executable statement of the ontology, so
that claims about what the model *is* can be tested rather than
asserted. Anyone wanting Wigner dynamics should use
`wpmwlib/phase_space_crystal_lattice.py`.

**Revision, July 2026.** §§2, 4, 5, 7 and 9 have been rewritten to remove
stored partnership indices, following
`docs/analysis/relational_pairing_and_carrier_lock.md`. Misalignment is
now defined over all ordered pairs of world-particle indices, mediation
is selected by cell and momentum row rather than by a label, and the sea
enters the vertex only through the order parameter of §2.5. The price is
one new postulate, (S) of §2.2; the returns are an exact factorisation
that removes the dominant cost, a $B$-independent vertex rule, and the
disappearance of a sea-consumption defect that the indexed formulation
concealed. No numerical result of the analysis note changes value.

Status of the pieces:

- §§1–5 (state, streaming, encounters, the vertex) are fully determined
  by the analysis note and implementable as written, given (S).
- §6 (the pump) is determined at first order; the steady state under
  continuous pumping is **[open]** — it is open item 2 of the analysis
  note, and the calibration constant of §8 is currently fixed by L1
  exactness rather than derived from a pump/drain balance.
- §7 (reconstruction) is standard.
- Items marked **[choice]** are implementation decisions not fixed by the
  physics; alternatives are noted.

Companion code: `src/demo_phase_alignment.py` verifies the kinematic and
vertex-level claims of §§2–5 in isolation, and
`src/demo_relational_pairing.py` verifies the coboundary, factorisation
and carrier-lock claims of §2.2, §2.5 and §4. A full live-sea integrator is
not yet in `wpmwlib`; §10 lists the validation ladder it must climb.

---

## 1. Domain, grids and constants

Periodic ring $x \in [0, L)$, one spatial degree of freedom, $\hbar$ and
$m$ as in the rest of the project. The momentum half-grid is

```math
dp \;=\; \frac{\pi\hbar}{L},
\qquad p_n \;=\; n \cdot dp,
\qquad n \in \mathbb{Z}.
```

By Theorem 1 of the phase-resonance note, **fundamental particles occupy
even $n$**; odd $n$ carries two-leg interference content only. The
simulation must therefore initialise all worldlines on even sites and
must never assign an odd momentum to a single leg. This is a hard
invariant and a good assertion to compile in.

The potential is a finite Fourier sum

```math
V(x) \;=\; \sum_{q \ge 1} V_q \cos\!\big(K_q x + \phi_q\big),
\qquad K_q = \frac{2\pi q}{L},
```

and each mode $q$ supplies a momentum quantum $2q \cdot dp$.

Cells for encounter detection: $M$ bins of width $\Delta x = L/M$. The
cell width is a numerical parameter and must be checked for convergence;
it sets the encounter rate and hence enters the calibration of §8.

---

## 2. State representation

### 2.1 Worldlines

The state is an ensemble of world-particles, each

```math
\big(x_j,\; p_j,\; \theta_j,\; \varepsilon_j\big),
\qquad
\varepsilon_j = \pm 1,
\qquad
\theta_j \in [0, 2\pi).
```

Two populations:

- **excess particles** — the positons sampled from the shifted Wigner
  function $W' = W + 2/h$, as in the existing microdynamics demos;
- **sea partners** — positon/negaton pairs, $B$ pairs per cell.

There is no third population, and nothing else is stored per particle.

### 2.2 Pairing

**No partnership is stored.** For any two world-particles the
misalignment $\mu_{ij}$ of §2.3 is the coboundary of the one-index field
$\Phi$, hence antisymmetric and additive,
$\mu_{ij} + \mu_{jk} = \mu_{ik}$, at every event. By Proposition R1 of
[`../analysis/relational_pairing_and_carrier_lock.md`](../analysis/relational_pairing_and_carrier_lock.md)
a stored index could therefore carry no state; it could only select. The
selection is recovered from kinematics: an ordered pair $(a, b)$ mediates
a $q$-mode vertex for an excess particle at row $n$ iff both legs occupy
the encounter cell, $a$ sits on row $n + q$, and $b$ sits on row $n - q$
(Corollary 4.3, with §1's momentum quantum $2q \cdot dp$).

**Postulate (S) — sea carrier lock.** *Sea particles sharing a cell and a
momentum row share a transported phase, up to the pumped misalignment.*

(S) is required, not optional. Lemma 1 of the phase-resonance note fixes
only the phase *difference* within a pair, so without (S) the residual
per-pair gauge group makes $\mu_{ij}$ gauge-variant for non-partners and
the all-pairs average vanishes as $1/B$ (Proposition R2). (S) is the
literal content of *phase-space **crystal** lattice*, and it is the
condition under which Lemma 5 holds: it is the sea-wide strengthening of
Lemma 1's per-pair darkness. It is testable and not a convention — see
rung 3a of §10 — and an implementation must instrument
$\lvert\hat Z_n(x)\rvert$ of §2.5 at runtime rather than assume it.

### 2.3 The derived quantity

No misalignment is stored. It is computed on demand from two worldlines'
leg data at the point where it is needed, for any pair of indices
$(i, j)$:

```math
\Phi_j(x,t) \;=\; \theta_j \;+\; \frac{p_j\,(x - x_j) - E_j\,t}{\hbar},
\qquad
\mu_{ij} \;=\; \Phi_i(x^{*},t^{*}) - \Phi_j(x^{*},t^{*}) \pmod{2\pi}.
```

Storing $\mu$ would be a bug: it is a function of the evaluation event,
and caching it invites evaluating a vertex with a stale value.

### 2.4 Sea depth

$B$ pairs per cell, with $B \approx 2\lvert W\rvert_{\max}\Delta x\thinspace dp/h$
scaled to the ensemble size. $B$ controls the amplitude bookkeeping of
the analysis note's open item 4 and must be swept: results that do not
converge in $B$ are not results.

### 2.5 The order parameter

Per cell $x_m$ and momentum row $n$, define

```math
Z_n(x_m, t) \;=\; \sum_{j \thinspace\in\thinspace (x_m,\thinspace n)} e^{i\Phi_j(x_m,t)},
\qquad
\hat Z_n \;=\; \frac{Z_n}{N_n},
```

the sum running over sea particles in that cell and row, with $N_n$ their
count. By Theorem R4 of the analysis note this is the *only* sea datum a
vertex needs. It is rebuilt from the $\Phi_j$ once per step and never
carried across steps, and $\lvert\hat Z_n\rvert$ is the runtime
diagnostic for (S): it equals $1$ under carrier lock and falls to
$O(N_n^{-1/2})$ for an incoherent sea.

---

## 3. Free evolution

Between vertices every world-particle streams and winds:

```math
\dot x_j = \frac{p_j}{m},
\qquad
\dot p_j = 0,
\qquad
\dot\theta_j = \frac{p_j^2/2m - V(x_j)}{\hbar}.
```

Momentum is constant on a leg — the potential acts only through the pump
(§6) and through vertices (§5), never as a classical force. This is the
sharpest structural difference from the classical-positon Monte Carlo in
`demo_cosine_well_microdynamics.py`, where $\dot p = -V'(x)$; here the
entire force emerges from vertex statistics.

Integrate with the same splitting used elsewhere in the project.
**[choice]** — exact free flight between vertex times is preferable to
fixed-step integration, since $p$ is piecewise constant and the phase
integral is then analytic on each leg.

---

## 4. Encounter detection

An **encounter** is the presence of an excess particle in a cell whose
sea rows $n \pm q$ are occupied, over a window $\tau_e$. There is no
pairwise loop: by Theorem R4 of the analysis note the sum over all
admissible mediating pairs factorises through the order parameter of
§2.5,

```math
\sum_{a,\thinspace b} \big[\thinspace w_0 + \kappa\thinspace\cos\mu_{ab}\thinspace\big]
\;=\; w_0\thinspace N_{n+q}\thinspace N_{n-q}
\;+\; \kappa\thinspace\mathrm{Re}\big(Z_{n+q}\thinspace\overline{Z_{n-q}}\big),
```

exactly and not as an approximation. The loop is therefore

```
for each cell:
    build Z_n, N_n for every occupied row      # one pass over the sea
    for each excess particle c in cell at row n:
        for s in (+1, -1):
            read Zhat_{n+sq}, Zhat_{n-sq}      # two numbers
            propose an exchange with window tau_e
```

Cost is $O(N_{\mathrm{exc}} + N_{\mathrm{sea}})$ per step rather than
$O(N_{\mathrm{exc}} \cdot B)$, and the sea-subsampling workaround of the
previous revision is no longer needed.

The window $\tau_e$ is a numerical parameter entering only through the
product with the coupling; see §8.

---

## 5. The vertex

This is the core of the specification. For an excess particle $c$ at
momentum $p_c = p_n$ in a cell, mediated by the sea rows $n + q$ and
$n - q$ of that cell:

### 5.1 Candidate exchange

By Theorem 4 the only exchange that survives averaging is the swap: the
excess particle arrives on the mate's momentum and leaves on the struck
partner's. With partnership removed this becomes a **row condition**. A
candidate exchange in direction $s = \pm 1$ is admissible iff both

```math
N_{n+sq} \;>\; 0
\qquad\text{and}\qquad
N_{n-sq} \;>\; 0
```

in the cell, in which case the out-state is

```math
p_c \;\mapsto\; p_{n + 2sq},
```

with one sea particle moved from row $n + sq$ to row $n - sq$ — the swap
of Theorem 4, read as a transfer of occupancy between the two mediating
rows. The excess particle's incoming row $n$ is the transition midpoint,
whose parity equals the parity of $q$ (Theorem 1); both mediating rows
and both particle rows remain even, so §1's invariant is preserved
automatically.

Implementations may either (i) test the condition and reject
non-matching encounters, or (ii) allow all exchanges and weight them by
the dephasing envelope

```math
\Big\lvert \mathrm{sinc}\Big(\tfrac{1}{2}\,\dot\mu\,\tau_e\Big) \Big\rvert,
\qquad
\dot\mu = \frac{\Delta p}{\hbar}\big(v_{\mathrm{exch}} - \bar v_{\mathrm{pair}}\big).
```

Route (ii) is the physically honest one and reduces to (i) as
$\tau_e \to \infty$; route (i) is much cheaper. **[choice]** — start
with (i), and use (ii) to check that the neglected off-stationary
traffic is genuinely negligible at the chosen $\tau_e$.

### 5.2 Reading the misalignment

Evaluate the order parameters of §2.5 at the encounter event
$(x^{*}, t^{*})$ from the current leg data of the sea particles in the
cell. Under (S) the argument of
$\hat Z_{n+sq}\thinspace\overline{\hat Z_{n-sq}}$ *is* the misalignment
$\mu$ of the predecessor note, and its modulus is $1$; departures of the
modulus from $1$ measure failure of (S) and must be logged, not clipped.
No search for any mate is required, and no particle refers to any other
particle by index.

### 5.3 Firing

Draw $u \sim \mathrm{Uniform}(0,1)$ and exchange if $u < P$, with either

```math
P \;=\; w_0 \;+\; \kappa\thinspace\mathrm{Re}\big(\hat Z_{n+sq}\thinspace\overline{\hat Z_{n-sq}}\big)
\qquad\text{(affine form)}
```

or the contact form

```math
P \;=\; \sin^2\!\big(\lvert h \rvert \tau_e \big),
\qquad
h = g_0 + g_1\thinspace\mu_1\thinspace e^{i\mu}.
```

The affine form is cheaper; the contact form is the one that derives
$\delta_0 = 0$ and saturates correctly at large coupling. **[choice]**.

Clamp: the affine form requires $\kappa \le \min(w_0, 1 - w_0)$. An
implementation must assert this rather than clipping silently, since a
clipped probability breaks the linear response the whole construction
exists to reproduce. The clamp is now the only condition needed, since
$\lvert\hat Z_{n+sq}\thinspace\overline{\hat Z_{n-sq}}\rvert \le 1$
holds identically (Corollary R4.1) and the rule is independent of the sea
depth $B$.

### 5.4 Phase update

On firing, **phase is continuous through the vertex for both legs**; only
momenta change.

- **Excess particle** — $\theta_c$ carries through. Its leg reference
  data are reset to $(x^{*}, p_{n+2sq}, \theta_c)$ at $t^{*}$.
- **Struck sea particle** — one particle is drawn uniformly from row
  $n + sq$ in the cell and moved to row $n - sq$, carrying its $\theta$
  through unchanged. Its contribution therefore migrates from
  $Z_{n+sq}$ to $Z_{n-sq}$: the vertex is a **transfer of coherence
  between rows**, which is the index-free reading of Corollary 4.2's
  "the pair exits aligned".
- There is no mate to leave untouched.

This replaces the previous revision's rule, in which the struck
partner's exit phase was solved from the condition $\mu = 0$ against a
stored mate. Under (S) the two agree at the level of the mean generator;
the equivalence is argued rather than proved, and rung 4 of §10 is the
test. See open item 2 of the analysis note.

On not firing, nothing changes; the encounter is discarded.

### 5.5 Invariants to assert

Every vertex must preserve, exactly and in floating point to rounding:

1. $\sum p$ over the two participating worldlines and the two mediating
   rows;
2. $\sum p^2/2m$ over the same (Corollary 4.1 — automatic, so a
   violation indicates a coding error, not a physics choice);
3. worldline count and species count;
4. the even-site invariant of §1 for every single leg;
5. $N_{n+sq} + N_{n-sq}$ in the cell, the occupancy transferred between
   the mediating rows summing to zero.

There is no partnership invariant, because there are no partnerships.

Asserting (2) is the single most valuable test in the whole
specification, because it fails loudly if the exchange has been
implemented as anything other than a swap.

---

## 6. The pump

The potential enters by acting on every sea particle through the same
phase vertex, over a pump interval $\tau_p$:

```math
\theta_j \;\mapsto\; \theta_j - \frac{V(x_j)\,\tau_p}{\hbar}.
```

At first order this displaces each pair's misalignment by

```math
\mu \;\mapsto\; \mu + \mu_1,
\qquad
\lvert\mu_1\rvert = \frac{V_q \tau_p}{\hbar} = C,
```

with the mode structure of Lemma 5: after the kick, $\mu$ at any event is
a function of the event alone. Two consequences the implementation should
exploit and test:

- **Populations are untouched at first order** (Lemma 3). A
  regression test should confirm that the $(x, p)$ histogram of the sea
  is unchanged by the pump to $O(V_p)$; if it is not, the pump has been
  implemented as a force rather than a phase.
- **$\mu$ becomes place-valued** (Lemma 5). Assert that
  $\lvert\hat Z_n(x)\rvert$ of §2.5 is $1$ to rounding immediately after
  a pump, in every occupied row of every cell. Under (S) this is the
  runtime form of the spread test of Part D of `demo_phase_alignment.py`,
  and it is strictly stronger: a spread test over stored pairs passes
  even when the pairs are mutually randomised, whereas
  $\lvert\hat Z_n\rvert$ does not (§7 of the analysis note, Model II).

**[open]** The steady state. Under continuous pumping, vertices drain
coherence — each firing migrates one sea particle between the mediating
rows — while the pump replenishes it. The fixed point of this balance
should supply the calibration constant of §8 from first principles.
It has not been derived, and until it is, the pump interval and the
vertex rate are tuned together, which remains the least satisfactory
part of this specification.

Note what is *no longer* part of this gap. In the previous revision the
drain had no matching source at all: by Corollaries 4.2 and 4.3 a struck
pair exits at $\Delta p = 0$ and can never mediate again, while by
Lemma 3 the pump cannot restore it, so the sea was a consumable resource
short by some three orders of magnitude for the cosine-well parameters
(Proposition R5 of the analysis note). Under §2.2 admissibility is
re-derived from row occupancy each step and no particle is ever spent, so
what remains open is the value of a constant, not the existence of a
steady state.

---

## 7. Reconstruction

The Wigner function is recovered exactly as in the existing
microdynamics demos: bin the excess particles on the $(x, p)$ grid and
subtract the background,

```math
W(x_m, p_n) \;=\; \rho_{\mathrm{emp}}(x_m, p_n) \;-\; \frac{2}{h}.
```

Sea particles are **not** binned into $W$: they are the medium, not the
state. The diagnostics worth accumulating separately are the per-cell,
per-row order parameter $\hat Z_n(x)$ of §2.5 — its modulus tracking (S)
and its argument the misalignment — and the quadrature
$\mathrm{Re}(\hat Z_{n+q}\thinspace\overline{\hat Z_{n-q}})$, whose
profile should reproduce $-V'(x)$ up to the calibration. This replaces
the misalignment distribution $f(x, \bar p, \Delta p, \mu)$ of the
analysis note's open item 2 with a two-index array of complex numbers,
which is both cheaper and, by Theorem R4, sufficient.

---

## 8. Calibration

The mean generator must reproduce

```math
\partial_t W_n \;=\; -\frac{p_n}{m}\,\partial_x W_n
\;+\; \Gamma_q(x)\,\big(W_{n+q} - W_{n-q}\big).
```

The advection term is exact by construction (§3). The collision term
carries one overall constant, the product of the encounter frequency,
the pump duty cycle and the vertex coupling. Fix it by requiring L1
exactness of the reconstructed evolution over one step, then verify —
do not re-fit — against $\Gamma_q(x) = -(V_q/\hbar)\sin(K_q x + \phi_q)$
with the corrected sign of
`docs/supplement/phase_space_crystal_lattice_supplement.md` §6.3.

The gross traffic is pure noise: the $w_0$ part of §5.3 produces
equal-and-opposite exchanges that cancel in the mean, and by Theorem R4
it enters only as $w_0 N_{n+q}N_{n-q}$, so it costs nothing to evaluate. This is the
$G$-freedom of `docs/analysis/four_rule_microdynamics_equivalence.md`,
and it means variance can be reduced by lowering $w_0$ without changing
any mean. **[choice]** — the smallest $w_0$ consistent with the clamp of
§5.3 is the variance-optimal one.

---

## 9. Reference loop

```
initialise:
    sample excess positons from W' = W + 2/h        (even momentum sites)
    populate B sea particles per cell per row       (carrier locked, (S))

for each step:
    pump:      theta_j -= V(x_j) * tau_p / hbar     for every sea particle
    stream:    advance x_j, theta_j on every leg    (p_j constant)
    assemble:  Z[m][n] <- sum exp(i Phi_j(x_m, t))  (one pass over the sea)
               N[m][n] <- count
               assert |Z[m][n]| / N[m][n] == 1      (postulate (S))
    encounter: for each excess particle c in cell m at row n:
                   for s in (+1, -1):
                       if N[m][n+s*q] == 0 or N[m][n-s*q] == 0:  skip
                       zz <- (Z/N)[m][n+s*q] * conj((Z/N)[m][n-s*q])
                       P  <- w0 + kappa * Re(zz)    (or the contact form)
                       if uniform() < P:
                           p_c <- p_c + 2*s*q*dp    (theta_c continuous)
                           move one sea particle    (theta continuous)
                               from row n+s*q to row n-s*q in cell m
                           assert momentum and energy unchanged
    diagnose:  bin excess particles -> W
               per-cell Re(Zhat[n+q] conj(Zhat[n-q])) -> compare with -V'(x)
```

---

## 10. Validation ladder

In order; each rung is a regression test, and none should be skipped.

1. **Kinematics.** Parts A–B of `demo_phase_alignment.py`: invariance
   and winding rates of $\mu$.
2. **Vertex algebra.** Part C: the swap solution, and conservation of
   both $\sum p$ and $\sum p^2$ without imposition.
3. **Pump.** Parts D–E: place-valued $\mu$; sea populations unchanged at
   $O(V_p)$; per-cell $\cos\mu$ profile proportional to $-V'(x)$.
3a. **Postulate (S), the deciding experiment.** Parts B–C and F of
   `demo_relational_pairing.py`. Build a sea whose pairs are individually
   aligned but mutually randomised and confirm that
   $\lvert\hat Z_n\rvert$ collapses to $O(N^{-1/2})$ while a
   partnered spread test still passes. This is the one rung that
   distinguishes (S) from the previous revision's convention, and it must
   be run before rung 4 is believed. Also confirm Theorem R4's
   factorisation identity to rounding, and the $\sqrt{N}$ variance
   advantage of Model I.
4. **Single step against the stencil.** One pump-plus-vertex step on a
   random $W$, compared with one Euler step of the mesh form in
   `wpmwlib/phase_space_crystal_lattice.py`. Agreement should be at the
   shot-noise floor and should improve as $N^{-1/2}$.
5. **Free particle.** $V = 0$: no pump, hence $\mu \equiv 0$, hence no
   exchanges fire beyond the cancelling gross traffic, and the evolution
   must reduce to the ballistic streaming already validated in
   `demo_cat_state_microdynamics.py`.
6. **Cosine well.** The problem of `demo_cosine_well_microdynamics.py`,
   run live, compared against the mesh QLE. This is the first test that
   exercises the steady state of §6 and is where the **[open]** item is
   expected to bite.
7. **Convergence sweeps.** In $B$, in $\Delta x$, in $\tau_e$, in
   $\tau_p$. A result that has not been swept in all four is not a
   result. The $B$ sweep is now a pure shot-noise check, since by
   Corollary R4.1 the vertex rule is exactly $B$-independent.

---

## 11. Known gaps

- **Steady state** (§6, **[open]**) — the calibration is fitted, not
  derived.
- **Row-index convention** — the analysis note's open item 1. The
  specification above is written entirely in momenta precisely to avoid
  depending on the unresolved convention; any translation to row indices
  must be checked first.
- **Multi-mode potentials** — §1 admits a Fourier sum, but first-order
  additivity of the mode contributions to $\mu$ is assumed and not
  verified.
- **Uniformly offset pairs** — a struck sea particle is generally
  separated from the others on its new row, so a nonuniform potential
  dephases it at rate $[V(x_a) - V(x_b)]/\hbar$ and degrades
  $\lvert\hat Z_n\rvert$. The specification does not model this drift;
  whether it matters over a run has not been estimated. In the present
  variables it is the question of whether (S) is dynamically preserved —
  open item 3 of the analysis note — and $\lvert\hat Z_n(x)\rvert$ is
  the instrument.
- **Phase continuity at the vertex** — §5.4's rule that both legs carry
  phase through the vertex replaces the previous revision's
  exit-alignment condition. Under (S) the two agree in the mean
  generator, but this is argued and not proved. Analysis note, open
  item 2.
- **Species bookkeeping** — §2.5 writes $Z_n$ without the factor
  $\varepsilon_j$. This is correct when $\varepsilon$ is constant within a
  row, which holds for the pumped sea, but the general bilinear form has
  not been fixed. Analysis note, open item 1.

Two entries of the previous revision have been closed rather than
carried forward. The **cost** entry — $O(N_{\mathrm{exc}} \cdot B)$ per
step — is superseded by Corollary R4.2. The **sea consumption** defect
diagnosed as Proposition R5 does not arise under §2.2.
