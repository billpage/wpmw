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

Status of the pieces:

- §§1–5 (state, streaming, encounters, the vertex) are fully determined
  by the analysis note and implementable as written.
- §6 (the pump) is determined at first order; the steady state under
  continuous pumping is **[open]** — it is open item 2 of the analysis
  note, and the calibration constant of §8 is currently fixed by L1
  exactness rather than derived from a pump/drain balance.
- §7 (reconstruction) is standard.
- Items marked **[choice]** are implementation decisions not fixed by the
  physics; alternatives are noted.

Companion code: `src/demo_phase_alignment.py` verifies the kinematic and
vertex-level claims of §§2–5 in isolation. A full live-sea integrator is
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

Each sea particle carries an index to its partner. Partnership is a
permanent property of the two worldlines: vertices change momenta and
phases, never partnerships. **[choice]** — an implementation may instead
re-derive pairing from co-location, but permanent indices make Lemma 4's
invariance testable and are strongly recommended.

### 2.3 The derived quantity

No misalignment is stored. It is computed on demand from the two
partners' data at the point where it is needed:

```math
\Phi_j(x,t) \;=\; \theta_j \;+\; \frac{p_j\,(x - x_j) - E_j\,t}{\hbar},
\qquad
\mu \;=\; \Phi_a(x^{*},t^{*}) - \Phi_b(x^{*},t^{*}) \pmod{2\pi}.
```

Storing $\mu$ would be a bug: it is a function of the evaluation event,
and caching it invites evaluating a vertex with a stale value.

### 2.4 Sea depth

$B$ pairs per cell, with $B \approx 2\lvert W\rvert_{\max}\Delta x\thinspace dp/h$
scaled to the ensemble size. $B$ controls the amplitude bookkeeping of
the analysis note's open item 4 and must be swept: results that do not
converge in $B$ are not results.

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

An **encounter** is a co-location of an excess particle with a sea
particle in the same cell during a window $\tau_e$.

```
for each cell:
    for each excess particle c in cell:
        for each sea particle a in cell:
            propose an encounter (c, a) with window tau_e
```

Cost is $O(N_{\mathrm{exc}} \cdot B)$ per step, which dominates the
simulation. **[choice]** — subsample the sea partners in a cell and
scale the rate correspondingly; this is statistically sound because the
vertex rule depends on the struck partner only through $p_a$, $p_b$ and
$\mu$, and by Lemma 5 the last of these is common to all pairs at a
given place.

The window $\tau_e$ is a numerical parameter entering only through the
product with the coupling; see §8.

---

## 5. The vertex

This is the core of the specification. For an encounter between excess
particle $c$ (momentum $p_c$) and sea partner $a$ (partner $b$):

### 5.1 Candidate exchange

By Theorem 4, the only exchange that survives averaging is the swap.
The candidate is therefore admissible **iff**

```math
p_c \;=\; p_b,
```

in which case the out-state is

```math
p_c \mapsto p_a,
\qquad
p_a \mapsto p_b,
\qquad
p_b \text{ unchanged}.
```

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

Evaluate $\mu$ at the encounter event $(x^{*}, t^{*})$ from §2.3, using
the partner's *current* leg data. No search for the mate's position is
required — only its stored leg parameters.

### 5.3 Firing

Draw $u \sim \mathrm{Uniform}(0,1)$ and exchange if $u < P$, with either

```math
P \;=\; w_0 \;+\; \kappa\thinspace\cos\mu
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
exists to reproduce.

### 5.4 Phase update

On firing:

- **Excess particle** — $\theta_c$ is continuous through the vertex. Its
  leg reference data are reset to $(x^{*}, p_a, \theta_c)$ at $t^{*}$.
- **Struck partner** — takes $p_b$, with exit phase chosen so the pair
  leaves aligned: solve $\mu = 0$ at $(x^{*}, t^{*})$ for
  $\theta_a^{\mathrm{out}}$. By Lemma 4 this is one equation in one
  unknown and always solvable from data present at the vertex.
- **Mate** — untouched in every respect.

On not firing, nothing changes; the encounter is discarded.

### 5.5 Invariants to assert

Every vertex must preserve, exactly and in floating point to rounding:

1. $\sum p$ over the three participating worldlines;
2. $\sum p^2/2m$ over the same (Corollary 4.1 — automatic, so a
   violation indicates a coding error, not a physics choice);
3. worldline count and species count;
4. the even-site invariant of §1 for every single leg;
5. partnership indices.

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
- **$\mu$ becomes place-valued** (Lemma 5). Assert that the spread of
  $\mu$ across pairs in a cell is at rounding level immediately after a
  pump. This is Part D of the companion demo, promoted to a runtime
  check.

**[open]** The steady state. Under continuous pumping, vertices drain
misalignment (each firing returns a pair to $\mu = 0$) while the pump
replenishes it. The fixed point of this balance should reproduce the
pinned $\Gamma_q(x)$ and should supply the calibration constant of §8
from first principles. Neither has been derived. Until it is, the pump
interval and the vertex rate are tuned together, which is the least
satisfactory part of this specification.

---

## 7. Reconstruction

The Wigner function is recovered exactly as in the existing
microdynamics demos: bin the excess particles on the $(x, p)$ grid and
subtract the background,

```math
W(x_m, p_n) \;=\; \rho_{\mathrm{emp}}(x_m, p_n) \;-\; \frac{2}{h}.
```

Sea particles are **not** binned into $W$: they are the medium, not the
state. Diagnostics worth binning separately are the misalignment
distribution $f(x, \bar p, \Delta p, \mu)$ of the analysis note's open
item 2, and the per-cell mean $\cos\mu$, whose profile should reproduce
$-V'(x)$ up to the calibration.

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
equal-and-opposite exchanges that cancel in the mean. This is the
$G$-freedom of `docs/analysis/four_rule_microdynamics_equivalence.md`,
and it means variance can be reduced by lowering $w_0$ without changing
any mean. **[choice]** — the smallest $w_0$ consistent with the clamp of
§5.3 is the variance-optimal one.

---

## 9. Reference loop

```
initialise:
    sample excess positons from W' = W + 2/h        (even momentum sites)
    populate B aligned pairs per cell               (mu = 0 everywhere)

for each step:
    pump:      theta_j -= V(x_j) * tau_p / hbar     for every sea particle
    stream:    advance x_j, theta_j on every leg    (p_j constant)
    encounter: for each (excess, sea partner) co-location:
                   if p_c != p_b:  skip             (or weight by sinc)
                   mu <- Phi_a(x*, t*) - Phi_b(x*, t*)
                   P  <- w0 + kappa * cos(mu)       (or the contact form)
                   if uniform() < P:
                       swap p_c <-> p_a
                       set theta_a so the pair leaves aligned
                       assert momentum and energy unchanged
    diagnose:  bin excess particles -> W
               bin per-cell <cos mu> -> compare with -V'(x)
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
   result.

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
- **Uniformly offset pairs** — a post-vertex pair has $\Delta p = 0$ but
  is generally separated, so a nonuniform potential dephases it at rate
  $[V(x_a) - V(x_b)]/\hbar$. The specification currently ignores this
  drift; whether it matters over a run has not been estimated.
- **Cost.** $O(N_{\mathrm{exc}} \cdot B)$ per step with $B$ in the
  thousands makes this orders of magnitude slower than the mesh form.
  The subsampling of §4 is the only mitigation currently identified.
