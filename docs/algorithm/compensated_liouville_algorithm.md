# Compensated Liouville Algorithm

**A complete specification for evolving the Wigner distribution as exact Newtonian flow plus a zero-mean signed mediated-jump channel.**

---

## 0. Status of this specification

This document promotes four analysis notes into an implementable
specification.

| note | step | theorems | what it supplies |
|---|---|---|---|
| [`../analysis/compensated_liouville_splitting.md`](../analysis/compensated_liouville_splitting.md) | 14 | C0–C7 | the split, its exactness, the reach condition, the quiet region |
| [`../analysis/eckart_barrier_compensated.md`](../analysis/eckart_barrier_compensated.md) | 15 | K1–K9 | the reach ceiling, the conserved classical outcome, a worked open-line problem with a closed form |
| [`../analysis/sea_population_equilibrium.md`](../analysis/sea_population_equilibrium.md) | 16 | S0–S9 | what an event *is*, and the rule that closes the ledger |
| [`../analysis/compensated_ontology.md`](../analysis/compensated_ontology.md) | 17 | G1–G5 | the four postulates this algorithm is a realisation of |

Companion verification: `src/demo_compensated_liouville_algorithm.py`, and for
§5 `src/demo_sea_population_equilibrium.py` and
`src/demo_emission_and_absorption.py`.

A tutorial introduction to §5, self-contained for a reader without the
background assumed here, is
[`../supplement/emission_and_absorption.md`](../supplement/emission_and_absorption.md).

**The one-sentence relation to the rest of `algorithm/`.** This is
[`phase_space_crystal_lattice_algorithm.md`](phase_space_crystal_lattice_algorithm.md)
with the rate field replaced and one step added: the mediated-jump rule of its
[§3b](phase_space_crystal_lattice_algorithm.md#3b-potential-driven-mediated-jumps) is used verbatim, the rate field $`\Gamma_q(x)`$ is replaced by a
compensated field $`\Gamma^{\mathrm{res}}_q(x)`$ computed from $`V`$ without
any Fourier decomposition of $`V`$, and a deterministic acceleration step is
placed in front of it.

**What the analysis note supplies** and this document does not re-derive: the
symbol, the exactness of the split (C1), its identification with the cubic
Taylor remainder (C2), the moment conditions at bounded reach (C3), the
crossover condition (C4), and the quiet region (C7).

**What is new here.** Five things, all of them consequences of putting the
split on a grid rather than on the real line, and all verified in the
companion demo:

1. §2 — the reach and the momentum cell are not merely "one parameter" in
   the loose sense of C4; the grid is *forced* by the reach, and the number of
   momentum cells is the ket–bra rung count. Neither is a convergence knob.
2. §3.2 — the deterministic step must be compensated against the **kernel's
   own first moment**, not against $`V'(x)`$. Using $`V'(x)`$ leaves up to
   $`4\times10^{-2}`$ of spurious force in the jump channel.
3. §4.1 — the Nyquist rung must be zeroed, or a real $`W`$ acquires an
   imaginary part at the $`10^{-2}`$ level every substep.
4. §4.4 — **the coherence horizon needs a profile, not a cutoff.** Under the
   hard cutoff of Definition (H) the total jump rate grows like $`\ln N_p`$ and
   the momentum churn grows like $`N_p`$, so neither has a limit and the
   discrete operator has no third moment — no semiclassical expansion at all.
   A soft profile fixes all three and is specified here as the default, with
   the hard cutoff as the special case $`w \equiv 1`$.
5. §4.1, §6 — the profile multiplies **the residual only**. Windowing the
   whole symbol is the obvious move and it silently destroys the Newtonian
   step.

Item 4 is also an erratum for the note's §4: the total-variation figures
tabulated there are quoted as absolute quantities and are in fact functions of
the grid used to compute them. The ratios are stable to about ten per cent and
the note's qualitative conclusion is unaffected.

Item 4 has a consequence for the ladder that this document cannot settle on its
own. Definition (H) of
[`../analysis/open_position_space.md`](../analysis/open_position_space.md) §3
states the horizon as a sharp fact — a pair beyond $`L_c`$ "is not
instantiated". What §4.4 shows is that the sharp version is the one choice of
profile under which the jump channel has an infinite event rate. If the
specified default here is right, Definition (H) wants restating as a **rung
occupancy** $`w(y)`$ rather than a cutoff, of which the indicator function is
one member and not a good one. That is a change to a postulate, and it is
flagged here rather than made.

**What the step 15–17 notes changed here.** Four things, of which the first is
a correction and not a clarification.

1. **§4.3's world-form rate was wrong by default.** It gave one rate
   $`R(x) = \sum_q|K_{\rm res}|`$ and one update rule. Read literally that
   prescribes the *emissive* unravelling — every event creates a pair — and
   S4 shows that relocates sea pairs without bound while S5 shows a finite sea
   then moves $`W`$ by 40 per cent in the core. There are two realisations of
   every event with identical observable content and opposite ledger content,
   and the mode is chosen **per event by supply**, not per channel by sign.
   §5 replaces §4.3's world form entirely.
2. **The world form needs state the mesh form does not.** $`W`$ is a
   *difference* of two populations, and the difference does not determine the
   populations. Three fields, not one; §5.1.
3. **The reach acquired a ceiling and a competing preference.** K1 caps
   $`y_{\max}`$ at the distance to the nearest complex singularity of $`V`$,
   while the ledger of §5 wants a reach several times larger. §2.4 states both
   and CLA7 records that they are not obviously compatible.
4. **The deterministic step acquired an exact invariant.** K4: the classical
   outcome functional is conserved by streaming plus the compensated force, so
   the entire quantum correction is delivered by the residual channel and can
   be measured as such. §3.4 makes it a required test rather than a result.

**The demography is now closed, and the vocabulary is fixed.** An event has
exactly two realisations, emissive and absorptive, chosen by supply. The two
body-moving realisations that would close the ledger exactly — Cyganski's
Right-Hop and Left-Hop — are excluded, because they return the discontinuous
momentum change postulate (S) exists to deny, and Proposition L1 in §5.7 shows
there is no variant of them that avoids it. So the residual ledger leak of one
to two per cent is not a gap in the specification but the measured price of the
ontology, and nothing in §9 offers a flag to turn it off.

Earlier drafts of this document called §4 "the hop channel". That was wrong in
the project's own terms: §4 is the crystal-lattice **mediated-jump** rule, and
a **hop** is the four-action rule that moves one body between momentum cells,
which is precisely the thing this algorithm does not do. §4.0 states both
definitions and flags the one note that still has them inverted.

Items marked **[choice]** are implementation decisions not fixed by the
underlying physics. Items marked **[open]** need further work. Items marked
**[normative]** are requirements: an implementation that violates one is
wrong, not merely different.

---

## 1. The split

### 1.1 The symbol

Let $`s`$ be conjugate to momentum, $`y = \hbar s/2`$ the half ket–bra
separation. Take the partial Fourier transform of $`W`$ in $`p`$ alone,
$`\hat W(x, s) = \int W(x,p)\thinspace e^{-ips}\thinspace dp`$. The entire
potential term of the QLE is then multiplication by

```math
M(x, s) \;=\; \frac{i}{\hbar}\bigl[V(x + y) - V(x - y)\bigr] ,
\qquad y = \frac{\hbar s}{2} .
```

Subtracting the part linear in $`y`$ gives the split this document
implements,

```math
M \;=\; \underbrace{i\thinspace V'(x)\thinspace s}_{M_{\mathrm{cl}}}
\;+\;
\underbrace{\frac{i}{\hbar}\bigl[V(x+y) - V(x-y) - 2y\thinspace V'(x)\bigr]}_{M_{\mathrm{res}}} .
```

### 1.2 Why the two steps can be applied in sequence

| | statement | where |
|---|---|---|
| C1 | $`e^{\tau M} = e^{\tau M_{\mathrm{cl}}}e^{\tau M_{\mathrm{res}}}`$ exactly, for every $`\tau`$ — both factors are multiplication operators in $`(x,s)`$ | note §2, §2.1 |
| C2 | $`M_{\mathrm{res}}`$ is the odd part of the cubic Taylor remainder of $`V`$ about $`x`$; it vanishes iff $`V''' \equiv 0`$ on the reach | note §2.2 |
| C3 | at bounded reach the residual kernel has zero zeroth **and** zero first moment: it conserves worlds and delivers no net momentum | note §3 |
| C7 | $`V''' \equiv 0`$ on $`[x-y_{\max}, x+y_{\max}]`$ implies no events at $`x`$ | note §5 |

C1 is what licenses the sequential application with **no Trotter error inside
the potential substep**. It says nothing about the Strang error between
streaming and the potential, which this algorithm still has; §6 is where that
is accounted for.

C3 is what makes the two steps non-overlapping: the whole classical force sits
in step 1 and step 2 does not apply it a second time.

### 1.3 What this algorithm is an implementation of

The split is a rearrangement of one operator, so as *numerics* it needs no
interpretation. As a **world-particle algorithm** it is the realisation of the
four postulates of
[`../analysis/compensated_ontology.md`](../analysis/compensated_ontology.md)
§1, and every normative requirement below traces to one of them:

| | postulate | what the code owes it |
|---|---|---|
| (E) | a world is a locally finite *signed* counting measure | two body fields, §5.1 |
| (A) | admissible ensembles are those whose mean is a Wigner transform of some $`\rho \succeq 0`$ | an initial-state check and a diagnostic; nothing in the step enforces it |
| (S) | between events every world streams under the **full** classical force, momentum continuous without exception | §3.1, §3.2, and §5.3 — jumps only, no hops: no code path may move a body between momentum rows |
| (D) | superimposed on (S), a birth-and-death process against a sea of neutral bound pairs | §4, §5 |

Theorem G1 assembles C1–C3, K4, K8 and S7 into the statement that (E), (S),
(D) reproduce the QLE in expectation. Two consequences are worth having in
front of an implementer.

**Theorem G2.** For quadratic $`V`$ the residual vanishes identically at every
reach, so the algorithm degenerates to exact Newtonian flow and takes **no
events at all**. Measured at $`1.2\times10^{-15}`$. This is validation rung 7
and it is the test the uncompensated scheme cannot pass.

**Theorem G4.** The residual channel is not only what makes the expectation
match the QLE; it is what keeps the ensemble inside the admissible set. Under
streaming alone in a quartic potential the least eigenvalue of the
reconstructed $`\rho`$ runs from the grid floor $`10^{-8}`$ to $`-0.10`$. So
switching the jump channel off is not a cheap approximation; it leaves the
state space.

---

## 2. Grids and parameters: the reach *is* the momentum grid

This section is the part of the specification most likely to be got wrong,
because it inverts the usual relation between a grid and an answer.

### 2.1 The primitives

Three numbers are chosen, and everything else follows.

| symbol | meaning | kind |
|---|---|---|
| $`y_{\max}`$ | the **reach**: half the coherence horizon $`L_c`$; the greatest distance from its own position at which a world consults $`V`$ | physical postulate |
| $`N_p`$ | the number of ket–bra **rungs** inside the horizon; even | physical postulate |
| $`M_x`$ | position cells | numerical |

$`y_{\max}`$ and $`N_p`$ are inherited from
[`../analysis/open_position_space.md`](../analysis/open_position_space.md) §3,
Definitions (H) and (R). They are statements about which ket–bra pairs exist.
Only $`M_x`$ is a resolution parameter in the ordinary sense.

### 2.2 What follows

```math
\Delta p \;=\; \frac{\pi\hbar}{2\thinspace y_{\max}} \;=\; \frac{\pi\hbar}{L_c},
\qquad
P_{\max} \;=\; \frac{N_p\thinspace\Delta p}{2},
\qquad
a \;=\; \frac{2 L_c}{N_p} ,
```

where $`a`$ is the rung spacing in the separation coordinate. The momentum
grid and the rung grid are exact Fourier duals of one another:

```
     rungs in separation                momentum transfers
     |y| <= y_max, spacing a/2          xi = q dp, |q| < N_p/2

  -y_max            +y_max                -P_max        +P_max
     |...............|          <--->        |.............|
     <-> a/2                                 <-> dp = pi hbar / L_c
     N_p rungs                               N_p cells

     resolution in y   <---->   extent in xi
     extent in y       <---->   resolution in xi
```

Read across: **the reach sets the momentum cell; the rung count sets the
momentum extent.** Verified in the demo (Part A): at $`y_{\max} = 1`$ and
$`N_p = 64`$, $`\Delta p = 1.5708`$, $`P_{\max} = 50.265`$, $`a = 0.0625`$,
and exactly $`64`$ rungs fall inside the horizon.

### 2.3 Refining the momentum grid is not a refinement **[normative]**

Halving $`\Delta p`$ doubles $`y_{\max}`$. It does not converge anything; it
asserts that worlds consult the potential twice as far away, which is a
different physical model with a different answer. In the limit
$`\Delta p \to 0`$ one recovers the untruncated Wigner equation, where the
residual is a differential operator rather than a jump measure and Theorem C3
fails outright.

The consequence for testing: **a grid-convergence study in $`\Delta p`$ is
meaningless for this algorithm.** The convergence parameters are $`\Delta t`$
and $`M_x`$, and those only.

The analogy that keeps this straight: the reach is an **aperture** in the
ket–bra separation, and the momentum-transfer spectrum is its diffraction
pattern. Narrowing the aperture coarsens the pattern — that is optics, not
error. §4.4 pushes the analogy one step further, and it turns out to predict
the event rate.

### 2.4 Choosing the primitives

- $`y_{\max}`$ from the physics if it is known; otherwise from C4, which wants
  $`u = k\thinspace y_{\max} \ll \pi`$ for the dominant scale $`k`$ on which
  $`V`$ varies. $`u \le \pi/4`$ is a reasonable working value: the residual
  is then about a tenth of the classical term. But see the ceiling and the
  conflict below.
- $`N_p`$ large enough that $`P_{\max}`$ covers the state's momentum support,
  and that $`\sigma_p`$ spans several cells. With a short reach $`\Delta p`$
  is large, so this is the binding constraint. Note the reach condition of C4
  reads $`\sigma_p \gg \hbar k/2`$ — a state too narrow in momentum to be
  resolved on the grid is also a state for which the split buys nothing.
- $`M_x`$ to resolve $`V'''`$, and to satisfy the CFL condition of §6.2.
- Position domain: the open line. See §2.5.

**There is a hard ceiling [normative].** Theorem K1: the expansion of $`M`$ in
odd powers of $`y`$ converges at $`x`$ if and only if

```math
y_{\max} \;<\; R(x) \;=\; \mathrm{dist}\bigl(0,\;
  \{\thinspace y : x \pm y \text{ is a singularity of } V\thinspace\}\bigr),
```

the distance to the nearest complex singularity of $`V`$ — not, as C6
suggested from the Coulomb case, to the nearest real one. For
$`V_0\thinspace\mathrm{sech}^2(x/a)`$ this is $`\sqrt{x^2 + (\pi a/2)^2}`$,
so $`\inf_x R = \pi a/2`$; for soft-core Coulomb it is the softening length
at the origin; for an entire $`V`$ such as a Gaussian barrier there is no
ceiling at all. Beyond the ceiling the residual is still computable but it no
longer has a semiclassical reading, and the reach is not a small parameter in
anything.

**And a competing preference, which is a genuine conflict [open].** The ledger
of §5 gets *better* with reach: the shortfall $`1/2 - f`$ measures
$`0.059, 0.020, 0.013`$ at $`y_{\max}/a = 2\pi, 4\pi, 8\pi`$ on the Eckart
barrier, with a cost-benefit knee near $`4\pi a`$. That is a factor of eight
above K1's ceiling of $`\pi a/2`$ for the same potential. Corollary K1.1
sharpens the bind: a reach-limited lattice resolves a packet only if
$`\sigma_r < a/2`$, narrower in position than half the barrier. The two
requirements are not jointly satisfiable on this potential and the
specification does not know which one yields; see CLA7. Until it is settled,
**quote the reach with every result** and do not treat a run inside the
ceiling and a run outside it as the same calculation.

### 2.5 The open line is the setting; the ring is a special case

**Nothing in this algorithm needs a position box.** The only transform is in
the rung direction, and the rung ladder is the coherence structure of a single
world, not a boundary condition. Position enters as a spectator parameter
throughout (note §2.1): the rate field at $`x`$ is built by sampling $`V`$ on
$`[x - y_{\max},\; x + y_{\max}]`$ and depends on nothing else.

So the specification is written for $`x \in \mathbb{R}`$, represented by a
window wide enough that

- $`W`$ is negligible at the edge, with margin $`P_{\max}\Delta t/m`$ so that
  no departure point of §3.3 leaves the window, and
- the window exceeds the active set of §7.4 by at least one reach.

Widening the window then adds only quiet cells: it costs storage and buys
nothing, which is the practical form of Theorem C7.

The ring enters in exactly three places, all of them special cases.

1. **Where a periodic $`V`$ quantises the jumps.** Theorem O5 gives a momentum
   quantum $`\pi\hbar/a`$ from the period of $`V`$, on all of $`\mathbb{R}`$,
   with no box. This is independent of the reach-induced $`\Delta p`$ of §2.2
   and the two coincide only when $`y_{\max} = a/2`$.
2. **Where the algorithm degenerates to the crystal lattice.** §4.2: at
   $`y_{\max} = L/2`$ on a ring of circumference $`L`$ with a single-mode
   $`V`$, the uncompensated rate field is exactly the two-atom stencil
   $`\Gamma_q`$.
3. **As a source of closed-form counterexamples**, per note §6. It is not a
   home: note §6.3 shows a ring pins every world at maximal reach, which is
   precisely the value of the one parameter the reach condition depends on at
   which the split buys nothing. **A ring is not a valid benchmark for this
   algorithm [normative].** §6's measurement uses a periodic window only
   because the test potential is a cosine, so that spectral advection is exact
   and the flow integrator does not contaminate the result.

---

## 3. Step 1 — deterministic acceleration

### 3.1 It must be a flow, not another split **[normative]**

Step 1 is the classical Liouville evolution

```math
\partial_t W \;=\; -\frac{p}{m}\thinspace\partial_x W \;+\; V'(x)\thinspace\partial_p W ,
```

and it must be integrated **as a flow** — by characteristics, or by a
symplectic integrator on worlds. Strang-splitting streaming against the
classical force inside step 1 would reinstate exactly the commutator term the
reordering was meant to remove, and the algorithm would gain nothing but an
extra substep.

This is the single structural requirement of the whole specification and it is
easy to violate by accident, because the classical force *is* available as a
multiplication operator in $`(x,s)`$ and applying it there is one line.

### 3.2 The force is the kernel's first moment **[normative]**

Theorem C3 says the residual delivers no net momentum. That is a statement
about the continuum. On the rung grid the first moment of the sampled kernel
is not $`V'(x)`$, so subtracting $`i V'(x) s`$ leaves a residue.

Instead, define the deterministic acceleration by the discrete first moment of
the full kernel — Theorem O2 of
[`../analysis/open_position_space.md`](../analysis/open_position_space.md)
read as a definition rather than as a result:

```math
V'_{\mathrm{eff}}(x) \;=\;
  \sum_{q} \xi_q\thinspace K(x, \xi_q) \;\Big/\; \sum_{q} \xi_q\thinspace B(\xi_q) ,
\qquad
a(x) \;=\; -\thinspace V'_{\mathrm{eff}}(x) ,
```

where $`K`$ is the kernel of $`M`$ on the rung grid, $`B`$ the kernel of the
reference symbol $`i s`$, both sums run over the whole grid, and $`a(x)`$ is
the acceleration the flow of §3.1 integrates. Then set
$`M_{\mathrm{cl}} = i\thinspace V'_{\mathrm{eff}}(x)\thinspace s`$ and
$`M_{\mathrm{res}} = M - M_{\mathrm{cl}}`$ as before.

With this definition the residual's discrete first moment vanishes to machine
precision at **every** $`N_p`$, and $`a(x) \to V'(x)`$ as $`N_p`$ grows.
Verified (Part C) at $`x = 1`$:

| potential | $`y_{\max}`$ | $`N_p`$ | residue, naive $`V'`$ | residue, discrete |
|---|---|---|---|---|
| cosine well | 1.00 | 32 | $`6.4\times10^{-3}`$ | $`2.3\times10^{-15}`$ |
| cosine well | 2.00 | 32 | $`2.3\times10^{-2}`$ | $`2.7\times10^{-15}`$ |
| Gaussian barrier | 2.00 | 32 | $`3.9\times10^{-2}`$ | $`3.1\times10^{-15}`$ |
| soft-core Coulomb | 2.00 | 32 | $`4.4\times10^{-2}`$ | $`2.2\times10^{-15}`$ |

against forces of order $`0.5`$ to $`0.8`$, so the naive compensation misses
by up to six per cent at a long reach and a small rung count. A spurious force
in a channel specified to carry none is not an accuracy question — it is a
violation of the conservation law the split exists to establish.

### 3.3 The two forms

- **Mesh.** Semi-Lagrangian transport: for each cell, integrate the
  Hamiltonian flow backwards over $`\Delta t`$ and interpolate $`W`$ there.
  **[choice]** of interpolant; cubic or better, since linear interpolation is
  diffusive enough to swamp the residual channel. Departure points must land
  inside the window (§2.5). Spectral advection is available only on a ring and
  is not the general form.
- **World ensemble.** A symplectic integrator per world — velocity Verlet is
  sufficient and is the recommended default **[choice]**. This is the natural
  open-line form: no interpolation, no boundary, and it is where the split
  pays, since it costs one force evaluation per world per step against a
  global transform for the jump channel.

### 3.4 The invariant step 1 conserves, and how to use it

**Theorem K4.** Let $`\Sigma`$ be any union of whole classical trajectories —
for a barrier, the classical outcome functional

```math
\Sigma \;=\; \bigl\{\, E > V_0,\; p > 0 \,\bigr\}
        \;\cup\; \bigl\{\, E < V_0,\; r > 0 \,\bigr\},
\qquad E = \frac{p^2}{2\mu} + V(r).
```

Then $`\int_\Sigma W`$ is exactly invariant under step 1, so its entire
motion is delivered by step 2.

This is the sharpest integration test the algorithm has, and it is cheap.
Step 1 is supposed to be the whole of the force; if $`\int_\Sigma W`$ drifts
under step 1 alone, either the flow integrator is wrong or §3.2's compensation
has leaked force into the wrong channel. Measured on the Eckart barrier over
the approach: $`0.500003, 0.499990, 0.499619, 0.499803, 0.500310`$ at
$`t = 0, 2, 4, 6, 8`$.

Read the other way it is the physical content of the whole reordering. The
classical answer is readable at $`t = 0`$ and never needs to be evolved; what
step 2 computes is the quantum correction and nothing else. On the Eckart
barrier that correction is $`0.0445`$ against a closed-form $`0.0441`$, and by
Theorem K5 it arrives as the small imbalance of two much larger opposed flows
across $`\partial\Sigma`$ — net/gross $`\approx 0.19`$ at $`\beta = 2`$,
falling as $`0.55/\beta`$. **The sampling cost of a transmission therefore
grows as $`\beta^2`$** (K6), which is worth knowing before a semiclassical
problem is attempted.

---

## 4. Step 2 — the mediated-jump channel

### 4.0 Two words this project uses precisely **[normative]**

They are one letter apart, they both describe momentum moving by
$`\pm\xi_q`$, and they are not the same thing. Getting them the wrong way
round produces code that is a correct solver and a wrong ontology.

| term | what it names | established in |
|---|---|---|
| **jump** | the *mediated* transfer: a world at $`(x, p)`$ mediates while momentum $`\pm\xi_q`$ is delivered at $`p \pm \xi_q`$. The mediator is unchanged. Under the compensated reading of §5 no body moves at all — a pair is ionised or bound. | [`phase_space_crystal_lattice_algorithm.md`](phase_space_crystal_lattice_algorithm.md) title and §3b, "single-rule mediated-jump process" |
| **hop** | one body moved from momentum cell $`n - q`$ to $`n + q`$. Right-Hop and Left-Hop, two of Cyganski's four actions. | [`../analysis/four_rule_microdynamics_equivalence.md`](../analysis/four_rule_microdynamics_equivalence.md) §2; used in this sense in [`../analysis/reach_energy_coupling.md`](../analysis/reach_energy_coupling.md) §4, where each hop changes $`T`$ by $`2p_n\xi_q/m`$ |

**This algorithm is a jump algorithm and contains no hops.** §4 is the jump
channel, reusing the crystal-lattice rule verbatim; §5.3 excludes the two hop
realisations because a hop moves a body between momentum rows, which postulate
(S) forbids. That the two coincide at the level of *mean-field occupancy* —
the mediated-jump rule is the pure-hop member of the four-rule family, §4 of
the four-rule note — is exactly why the distinction has to be made at the
world-particle level and cannot be read off the generator.

*Erratum flagged, not made.* §8.1 of
[`../analysis/eckart_barrier_compensated.md`](../analysis/eckart_barrier_compensated.md)
is titled "a hop is not a jump" and uses the two words the other way round
from the table above. Its argument is right and its labels are inverted; the
correction belongs in that note.

### 4.1 Building the rate field

For each position cell $`x_m`$:

1. Evaluate $`V`$ at the $`2N_p`$ points $`x_m \pm y_k`$, where
   $`y_k = \hbar s_k/2`$ are the rungs, and form $`M(x_m, s_k)`$.
2. **Zero the Nyquist rung [normative].** With $`N_p`$ even the grid holds
   $`y = -y_{\max}`$ but not $`+y_{\max}`$, so the odd symbol is sampled
   asymmetrically and one rung has no partner. The symmetric treatment of a
   maximal-separation rung visited from both sides is the average of the two
   endpoint values, and for an odd function that average is zero. Set
   $`M(x_m, s_{N_p/2}) = 0`$, and the same for the reference symbol $`is`$ used
   in §3.2.
3. Form $`M_{\mathrm{res}} = M - M_{\mathrm{cl}}`$ per §3.2.
4. **Apply the horizon profile to the residual, and only to the residual
   [normative]:** $`M_{\mathrm{res}} \leftarrow w(y/y_{\max})\thinspace
   M_{\mathrm{res}}`$, with $`w`$ as specified in §4.4. $`M_{\mathrm{cl}}`$ is
   left alone.
5. Inverse-transform in $`s`$ to get the signed transfer rates
   $`K_{\mathrm{res}}(x_m, \xi_q)`$.
6. Set $`\Gamma^{\mathrm{res}}_q(x_m) = -K_{\mathrm{res}}(x_m, \xi_q)`$.

Step 2 is not cosmetic. Verified (Part B), cosine well, $`y_{\max} = 1`$,
$`N_p = 64`$:

| | $`\max\lvert\mathrm{Im}\thinspace K\rvert/\max\lvert K\rvert`$ | oddness in $`q`$ | $`\max\lvert\mathrm{Im}\thinspace W\rvert`$ after one substep |
|---|---|---|---|
| as sampled | $`1.15\times10^{-1}`$ | $`2.30\times10^{-1}`$ | $`1.63\times10^{-2}`$ |
| Nyquist zeroed | $`1.09\times10^{-16}`$ | $`2.34\times10^{-16}`$ | $`5.6\times10^{-16}`$ |

Left alone, a single rung drives a real $`W`$ complex at the percent level
every substep, and the realness and oddness of the kernel — the properties
that give number conservation for free and force the sign structure — are lost
at the same order.

A soft horizon does not make step 2 redundant. The profile zeroes the residual
at that rung for free, since $`w(\pm 1) = 0`$, but $`V'_{\mathrm{eff}}`$ is
read off the *unwindowed* $`M`$ and stays contaminated.

**Why the profile goes on the residual only.** $`w(s)\thinspace i V'(x) s`$ is
not the symbol of a drift, so windowing the whole symbol takes the Newtonian
step apart. The test that shows it is validation rung 7: for an open-line
quadratic $`V`$ the residual is identically zero, so the algorithm must reduce
to exact Newtonian flow. Verified (Part H), $`V = m\omega^2x^2/2`$ with
$`\omega = 1.7`$, at $`x = 1`$:

| placement | $`\max\lvert M_{\mathrm{res}}\rvert`$ | error in the force |
|---|---|---|
| hard horizon | $`4.6\times10^{-13}`$ | $`2.8\times10^{-14}`$ |
| window the whole symbol | $`2.2\times10^{1}`$ | $`5.5\times10^{-2}`$ |
| window the residual | $`1.2\times10^{-13}`$ | $`2.8\times10^{-14}`$ |

Windowing the whole symbol manufactures a residual of size 22 in a potential
that has none, and shifts the classical force by two per cent.

Read physically the rule is not a trick: the horizon grades the **coherence**
channel. The classical force is the local, first-moment part of the kernel and
is not a coherence effect, so there is nothing there for a rung occupancy to
attenuate.

**No Fourier decomposition of $`V`$ appears anywhere.** The rate field is
built by sampling $`V`$ over each world's own reach. This is the note's §1
point as an implementation fact, and §6.2 there is the reason it matters: the
mode-by-mode bound on the residual diverges even where the residual is exactly
zero.

### 4.2 The update rule — unchanged

With $`\Gamma^{\mathrm{res}}_q`$ in hand, the rule is §3b/§3c of
[`phase_space_crystal_lattice_algorithm.md`](phase_space_crystal_lattice_algorithm.md)
verbatim:

```math
W(x_m, p_n) \;\mathrel{+}=\; \Delta t \sum_{q \ge 1}
  \Gamma^{\mathrm{res}}_q(x_m)\thinspace
  \bigl[\thinspace W(x_m, p_{n+q}) - W(x_m, p_{n-q})\thinspace\bigr] .
```

A world at cell $`(m,n)`$ mediates a transfer between cells $`n+q`$ and
$`n-q`$; the mediator is unchanged; no world is created or destroyed. Verified
(Part D) against the spectral action of $`M_{\mathrm{res}}`$ to
$`7.1\times10^{-16}`$ relative.

The correspondence is exact in the other direction too. Run the *uncompensated*
field at maximal reach $`y_{\max} = L/2`$ on a single-mode potential and the
kernel collapses onto the two atoms $`q = \pm 1`$ with
$`-K(x, \xi_1) = \Gamma_1(x)`$ to $`8.9\times10^{-16}`$, everything else below
$`2.2\times10^{-16}`$. So the crystal-lattice algorithm is the
$`\Gamma^{\mathrm{res}} \to \Gamma`$, $`y_{\max} \to L/2`$ corner of this one.

### 4.3 Mesh and world forms

- **Mesh, exact.** Apply $`\exp(\Delta t\thinspace M_{\mathrm{res}})`$ in the
  $`(x, s)`$ representation: FFT each momentum column, multiply, inverse FFT.
  Exact in $`\Delta t`$ for the substep, and this is the recommended default
  for the mesh form. Explicit Euler on the stencil of §4.2 amplifies
  negativity and must not be used for production runs **[normative]** — this
  is the standing project result on the momentum substep, unchanged here.
  The mesh form carries one field and needs nothing from §5.
- **World ensemble, sampled.** Total rate at $`x`$ is
  $`R(x) = \sum_{q \ne 0} \lvert K_{\mathrm{res}}(x, \xi_q)\rvert`$. Draw an
  event time from $`R(x)`$ and draw $`q`$ with probability
  $`\lvert K_q\rvert / R`$. **What to then do with the drawn event is §5, and
  is not obtainable from $`R`$ alone [normative].**

The sign of $`K_q`$ has to go somewhere, and there is no way to make it go
away. $`K_{\mathrm{res}}`$ is odd and non-zero, hence takes both signs, hence
is never the generator of a one-body Markov jump process — Proposition T3 of
[`../supplement/takabayasi_1954_stochastic_picture.md`](../supplement/takabayasi_1954_stochastic_picture.md),
unweakened by compensation. The two admissible unravelings are the project's
standing pair:

- **positon/negaton [normative for the ontology]** — the sign is a species
  label, and the ledger of §5 is what makes it work. This is postulate (D) and
  it matches
  [`../analysis/species_sectors_and_annihilation.md`](../analysis/species_sectors_and_annihilation.md).
- **signed weights [choice]** — weights multiply by $`\mathrm{sign}(K_q)`$ and
  every weight grows as $`e^{R t}`$ to compensate the depletion. Cheaper per
  step, with the sign problem relocated into variance. Legitimate as a solver;
  it carries no ledger and represents no ontology, so nothing in §5 applies to
  it and nothing in §5 can be tested with it.

Compensation shrinks $`R`$; it does nothing to the sign structure. That
distinction is worth keeping sharp, because "the quantum channel now carries
no force" invites the misreading that it has become classical.

**Erratum.** Until step 16 this section gave $`R`$ and the update rule and
stopped, which reads as an instruction to *create* a pair at every event. That
is the emissive unravelling, and Theorems S4 and S5 show it is the wrong one:
it drains the worst sea cell without bound, and a sea deep enough to survive
is a sea whose throttling moves $`W`$ by 40 per cent in the core. §5 is the
replacement.

### 4.4 The horizon profile, and the event budget

**The specified horizon is soft [normative].** The occupancy profile
$`w(y/y_{\max})`$ must be

| requirement | what it buys |
|---|---|
| even | both moment conditions of C3 survive **exactly** |
| $`w(0) = 1`$ | the leading $`\hbar^2`$ Moyal term is untouched |
| $`w(\pm 1) = 0`$ | no seam, so the event rate converges |
| $`w'(\pm 1) = 0`$ | $`1/q^3`$ tail, so the momentum churn converges |
| supported on the reach | the quiet region of C7 survives |

The raised cosine $`w(t) = \cos^2(\pi t/2)`$ meets all five and is the
default **[choice]**. The hard cutoff of Definition (H) is the special case
$`w \equiv 1`$, and it is the one choice that fails the last three.

**Why the hard cutoff fails.** It is a hard-edged aperture in the separation
coordinate. $`M_{\mathrm{res}}`$ is odd, so its periodic extension over the
rung window carries a jump of $`2\lvert M_{\mathrm{res}}(x,s_{\max})\rvert`$ at
the seam, and a jump gives coefficients falling as $`1/q`$:

```math
\lvert K_{\mathrm{res}}(x, \xi_q)\rvert \;\simeq\;
  \frac{\lvert M_{\mathrm{res}}(x, s_{\max})\rvert}{\pi\thinspace q}
\qquad\Longrightarrow\qquad
R(x) \;\simeq\; \frac{2}{\pi}\thinspace
  \bigl\lvert M_{\mathrm{res}}(x, s_{\max})\bigr\rvert \thinspace \ln N_p \;+\; O(1) .
```

Verified (Part E), fitted slope against predicted, at $`x = 1`$:

| potential | $`y_{\max}`$ | $`R`$ at $`N_p = 256`$ | $`R`$ at $`N_p = 16384`$ | fitted slope | $`(2/\pi)\lvert M_{\mathrm{res}}(s_{\max})\rvert`$ |
|---|---|---|---|---|---|
| cosine well | 0.25 | 0.00714 | 0.01429 | 0.00170 | 0.00170 |
| cosine well | 1.00 | 0.44541 | 0.89001 | 0.10598 | 0.10573 |
| cosine well | 2.00 | 3.27872 | 6.52013 | 0.77268 | 0.77085 |
| Gaussian barrier | 1.00 | 0.23750 | 0.45822 | 0.05292 | 0.05285 |
| Gaussian barrier | 2.00 | 5.66761 | 11.10204 | 1.29542 | 1.29236 |

Four figures of agreement. But the diverging event count is the mildest of the
three symptoms. Cosine well, $`y_{\max} = 1`$, measured tail exponent in
brackets:

| horizon | quantity | $`N_p = 256`$ | 1024 | 4096 | 16384 |
|---|---|---|---|---|---|
| hard $`(q^{-1.03})`$ | $`R = \sum_q\lvert K\rvert`$ | 0.44541 | 0.59539 | 0.74308 | 0.89001 |
| | $`\sum_q\lvert \xi K\rvert`$ | 14.35 | 58.56 | 235.4 | 942.7 |
| | $`\sum_q \xi^3 K`$ | 31.91 | 128.98 | 517.2 | 2070.3 |
| raised cosine $`(q^{-3.00})`$ | $`R`$ | 0.01631 | 0.01657 | 0.01663 | 0.01665 |
| | $`\sum_q\lvert \xi K\rvert`$ | 0.05339 | 0.05447 | 0.05475 | 0.05482 |
| | $`\sum_q \xi^3 K`$ | $`-0.13381`$ | $`-0.12981`$ | $`-0.12880`$ | $`-0.12847`$ |

The event rate diverges logarithmically, the **momentum churn diverges
linearly**, and the third moment diverges linearly too — so under a hard
horizon the discrete operator has no third moment, which is to say no
semiclassical expansion, term by term, at all.

Under the soft horizon all three converge, and the third moment converges to
$`-0.128465`$ — which is the analytic leading Moyal coefficient,
$`-0.128465`$, to six figures. **A soft horizon is what gives the discrete
operator a Moyal expansion.** That is the argument for making it the default,
and it is a stronger one than the numerical convenience: a hard cutoff makes
the model's own semiclassical limit ill-defined.

Two things the profile does not cost.

- **Conservation.** $`w`$ even keeps
  $`\sum_q K_{\mathrm{res}} = \sum_q \xi_q K_{\mathrm{res}} = 0`$ exactly, for
  any profile. The compensation of §3.2 is untouched.
- **The leading quantum correction.** With $`M_{\mathrm{res}}(0) =
  M_{\mathrm{res}}'(0) = M_{\mathrm{res}}''(0) = 0`$ and $`w`$ even,
  $`(w M_{\mathrm{res}})'''(0) = w(0)\thinspace M_{\mathrm{res}}'''(0)`$, so
  $`w(0)=1`$ preserves the $`\hbar^2`$ term exactly. The profile attenuates only
  the far rungs — which are exactly the ones producing the sidelobes.

**Quote $`N_p`$ with any event rate anyway.** Even under a soft horizon $`R`$
depends on $`N_p`$ until it has converged, and $`N_p`$ is a physical parameter
(the rung count), not a tolerance. This is the erratum of §0.

**The reach still controls the rate.** The seam value obeys the C2 bound
$`\lvert M_{\mathrm{res}}\rvert \le (2/\hbar)(y_{\max}^3/6)\sup\lvert V'''\rvert`$,
so short reach means cheap channel under either horizon. C4 survives intact.

![Event budget under a hard and a soft horizon](https://raw.githubusercontent.com/billpage/wpmw/output/figures/compensated_algorithm_budget.png)

## 5. The ledger — what the world form has to carry

The mesh form evolves one field and is complete. The world form is not:
it represents $`W`$ as a *difference* of two populations, and a difference
does not determine its terms. This section specifies the missing structure.
It has no counterpart in
[`phase_space_crystal_lattice_algorithm.md`](phase_space_crystal_lattice_algorithm.md)
and it is the most recently settled part of the algorithm.

### 5.1 Three fields, one of which is the observable

| field | meaning | what constrains it |
|---|---|---|
| $`E = u^+ - u^-`$ | signed density; **the observable**, $`E \equiv W`$ | the QLE, exactly |
| $`N = u^+ + u^-`$ | body density | the unravelling rule, and nothing else |
| $`S`$ | bound sea pairs per cell, background $`B = 2/h`$ | the same rule |

**Theorem S2, first half.** The $`E`$ equation is closed and is exactly the
QLE whatever $`N`$ and $`S`$ do. Two consequences, and the second is the
uncomfortable one. Nothing chosen in this section can corrupt the observable.
And nothing measured on the observable can validate anything chosen in this
section — the ledger has to be checked on its own terms, which is what
validation rungs 11 to 13 are for.

Two structural facts an implementation must maintain **[normative]**:

- $`u^\pm \ge 0`$ pointwise. Clamp after every allocation.
- $`N \ge \lvert E\rvert`$ pointwise, by definition. Spectral transport rings,
  so re-impose it after each transport substep; measured excursions are
  $`1.6\times10^{-3}`$ at $`\Delta t = 0.02`$ and $`6.7\times10^{-4}`$ at
  $`0.01`$, i.e. first order and not structural.

### 5.2 What an event is

$`K_{-q} = -K_q`$, so the $`(q, -q)`$ pair is **one event with two legs of
opposite sign**, not two channels. A parent at momentum row $`p`$ deposits
$`+1`$ at $`p + \xi_q`$ and $`-1`$ at $`p - \xi_q`$.

**Proposition S0 / K8 (co-location is forced).** The bound pair the event
consumes or produces sits at the parent's *own* momentum row.

*Proof.* The children carry $`(p + \xi_q) + (p - \xi_q) = 2p`$. A neutral pair
drawn from row $`p'`$ carries $`2p'`$, so $`p' = p`$. Nothing else can supply
the difference, because by postulate (S) no world-particle — the parent
included — changes its momentum discontinuously. $`\square`$

So the event is ionisation or recombination of a pair co-located with the
parent in *both* $`r`$ and $`p`$: the $`c = 0`$ member of the dark family of
[`../analysis/species_sectors_and_annihilation.md`](../analysis/species_sectors_and_annihilation.md),
not the $`c = 2`$ crystal shift. This is a derivation, not a postulate, and it
is why the sea has to be indexed by momentum row and not merely by cell.

### 5.3 Two realisations, and two that are excluded

> **Arriving here from a link?** This section is normative and terse, and it
> assumes the rest of the document. Three places to get the background:
>
> - [`../supplement/emission_and_absorption.md`](../supplement/emission_and_absorption.md)
>   is a tutorial treatment of exactly this table, written to stand alone. It
>   builds the signed ensemble, the sea, and the event from scratch, draws both
>   realisations in phase space *and* as space-time worldlines for all three
>   species, and gives the demographic reading of why the absorptive fraction
>   finds $`1/2`$ without being told to.
> - [`../analysis/compensated_liouville_splitting.md`](../analysis/compensated_liouville_splitting.md)
>   (theorems C0–C7) derives the split that produces this residual channel in
>   the first place — why there is anything left over once the classical force
>   has been removed exactly.
> - [`../analysis/eckart_barrier_compensated.md`](../analysis/eckart_barrier_compensated.md)
>   §8 shows what the same two realisations look like as a physical process:
>   tunnelling through a barrier, read as two large opposed flows of pairs
>   across the classical separatrix.

```
   momentum row                     deposit   emissive      absorptive
   ------------------------------   -------   -----------   -----------
   p + xi_q   daughter A              +1      create  +      remove  -
   p          parent (sea pair)        0      ionise it      bind it
   p - xi_q   daughter B              -1      create  -      remove  +
                                              dN = +2        dN = -2
                                              dS = -1        dS = +1
```

Both realisations are identical in $`E`$ and opposite in the ledger. **These
two are the whole of the algorithm's demography [normative].**

Two further realisations are arithmetically available and are excluded. They
are the **hops** of §4.0 — Cyganski's Right-Hop and Left-Hop, one body carried
from momentum cell $`n - q`$ to $`n + q`$:

| | at $`p+\xi_q`$ | at $`p-\xi_q`$ | $`\Delta N`$ | $`\Delta S`$ | what moves |
|---|---|---|---|---|---|
| $`M^+`$ | create $`+`$ | remove $`+`$ | 0 | 0 | one positon, by $`+2\xi_q`$ |
| $`M^-`$ | remove $`-`$ | create $`-`$ | 0 | 0 | one negaton, by $`-2\xi_q`$ |

Each is neutral in *both* ledgers, which is exactly why they are tempting; and
each moves a body between momentum rows, which is exactly why they are not
available under postulate (S). §5.7 records what they would have bought and why
the price is not payable.

**Theorem S6 (the mode is per event, not per leg).** A mixed realisation — one
leg created and the other removed — is one of $`M^\pm`$, and changes body
momentum by $`2\xi_q`$, which nothing can supply. Hence an event is wholly
emissive or wholly absorptive, and absorption needs a partner of the right
species at **both** daughters. Absorption is therefore supply limited; when
supply fails, the event falls back to emissive.

The distinction is the whole ontological content of the compensated split. A
jump changes *how many worlds there are*; a hop would change *where one of them
is*. Only the first is available.

### 5.4 The ledger identity, and why it needs no tuning

**Theorem S7.** With absorptive fraction $`f`$ over $`n_{\rm ev}`$ events,

```math
\Delta N = 2\,(1 - 2f)\, n_{\rm ev},
\qquad
\Delta S = (2f - 1)\, n_{\rm ev},
```

so $`f = 1/2`$ closes the body ledger and the sea ledger **simultaneously**.
Not two problems with two knobs: one problem with one number. Verified here to
$`2\times10^{-14}`$ and $`5\times10^{-12}`$ respectively.

**Theorem S9.** $`f`$ is an attractor near $`1/2`$, reached from both sides and
independent of the initial ensemble. The mechanism is S7 read as feedback: if
$`f > 1/2`$ bodies drain, partners grow scarce, and $`f`$ falls; if
$`f < 1/2`$, the reverse. Padding the initial ensemble with $`\pm`$ pairs to
$`\rho = N_0/\lvert E\rvert = 1`$ and $`20`$ — which leaves $`E`$ untouched —
gives $`f = 0.462`$ and $`0.836`$ at $`t = 0.5`$, agreeing to within $`0.015`$
by $`t = 8`$.

**No separate recombination rate is required for closure.** Absorption is
itself a pair-removal channel satisfying S1, and S9 says it suffices. A
bilinear recombination sink $`-\tfrac{\kappa}{2}(N^2 - E^2)`$ with
$`\kappa \propto \Gamma_{\rm tot}(x)`$ remains available as a population
knob **[choice]** — the proportionality is forced, or the ledger acquires a
reach dependence — but it is no longer load-bearing, and its one dimensionless
constant is still not fixed by anything (S-SP1).

**The attractor sits slightly below $`1/2`$, and the shortfall is real.** It is
flat in $`\Delta t`$ ($`0.4718, 0.4706, 0.4704`$ at
$`\Delta t = 0.02, 0.01, 0.005`$) and controlled by the reach
($`1/2 - f = 0.059, 0.020, 0.013`$ at $`y_{\max}/a = 2\pi, 4\pi, 8\pi`$). The
mechanism is the **conjunction cost** of S6: each leg finds its partner 67 to
76 per cent of the time, but the two availabilities are strongly
anti-correlated — mean disagreement $`\approx 0.5`$ against means of
$`0.7`$ — so typically one leg has a partner and the other does not. Requiring
both costs about a quarter of all events, and that cost barely moves with
reach. Whether the shortfall vanishes or floors near $`0.013`$ is S-SP3.

### 5.5 The allocation rule **[normative]**

Per position cell, per channel $`q \ge 1`$, per parent species, with demand
$`D = \lvert K_q(x)\rvert\thinspace u^{\rm parent}\thinspace\Delta t`$:

1. Read the partner supply at the two daughters, $`\mathrm{cap}_A`$ at
   $`p + \xi_q`$ and $`\mathrm{cap}_B`$ at $`p - \xi_q`$, each being the
   *opposite* species to the sign the event deposits there.
2. $`A = \min(D,\; \mathrm{cap}_A,\; \mathrm{cap}_B)`$ absorptive events;
   $`\mathrm{Em} = D - A`$ emissive.
3. Apply both: absorptive removes one body from each daughter and credits
   $`S`$ at the parent row; emissive adds one body to each daughter and debits
   $`S`$ at the parent row.
4. Clamp $`u^\pm \ge 0`$.

The sign that decides which species each leg wants is
$`\mathrm{sgn}\,K_q`$ times the parent's own sign. This is the point the
project has had to correct twice: **the mode is chosen per event by supply,
not per channel by sign.** A channel with $`K_q < 0`$ does not "run backwards";
it deposits the opposite species, and whether that deposition is realised by
adding or by removing depends on what is locally available.

Caps must be read **live**, against the populations as they stand when the
channel is reached, not against a snapshot taken at the top of the step;
otherwise the same partner is spent twice. That makes the tau-leap
order-dependent. Open item S-SP4 measures it: the observable moves by
$`3.0`$ to $`3.5\times10^{-3}`$ against a splitting error of
$`1.0\times10^{-2}`$, and $`f`$ from $`0.431`$ to $`0.442`$ — both inside the
error the method already carries. But $`N(T)`$ spreads by about 14 per cent
across orderings, so **every body count and sea-depth figure carries an
$`O(15\text{ per cent})`$ method uncertainty that the fidelity figures do not**.

### 5.6 The sea is a per-cell counter, not a particle list **[normative]**

$`B = 2/h`$ is a density of *pairs*, so $`Bh = 2`$: every phase-space cell of
area $`h`$ holds exactly two bound sea pairs, independent of any grid. That
number is not a parameter and not a discretisation artefact.

It is also why the sea must not be carried as particles. With a Monte-Carlo
multiplier $`\nu`$ the sea scales as $`4\nu`$ bodies per $`h`$-cell times the
phase-space volume of the window, so at $`\nu = 10^6`$ it is
$`\mathcal{O}(10^8)`$ bodies for the Eckart problem — and it *grows if the
window is merely widened*, which is not physics. The sea is uniform, dark,
inert, and analytically known everywhere, so it should be a scalar counter per
$`(x, p)`$ cell, materialised only when an event touches it. Only the unpaired
population needs an explicit list, and at $`\nu = 10^6`$, $`\rho = 3`$ that is
three million bodies — ordinary.

This also relocates the cost. It is not storage; it is the partner lookup in
two momentum rows that every event of §5.5 performs. The per-cell occupancy
representation matters more than the particle representation for exactly that
reason.

That the sea:unpaired ratio depends on how wide a window was drawn is
uncomfortable if the sea is read as real rather than as bookkeeping. The
presumed resolution — that sea outside the packet's support never participates,
because there is no demand — is not established by anything measured. See
CLA8.

### 5.7 Standing population, and why the ledger is allowed to leak

The minimal ensemble cannot absorb at all to begin with: for a positive $`W`$,
$`u^- = 0`$ everywhere at $`t = 0`$, so $`f`$ starts at zero and can rise only
as emissive events manufacture the partners absorption later consumes.
Absorption then starves at the moving frontier of the minority species.

**The specified repair is a standing supra-minimal population [choice].** Pad
the initial ensemble with $`\pm`$ pairs to $`\rho = N_0/\lvert E\rvert`$, which
leaves $`E`$ and therefore every observable untouched. On the Eckart barrier at
$`y_{\max} = 4\pi a`$:

| $`\rho`$ | $`1/2 - f`$ | starved demand | rel $`L^2`$ against the mesh QLE |
|---|---|---|---|
| 1 | 0.0325 | 11.6% | $`5.5\times10^{-2}`$ |
| 3 | 0.0247 | 5.1% | $`5.3\times10^{-2}`$ |
| 10 | 0.0114 | 2.0% | $`8.2\times10^{-2}`$ |
| 100 | 0.0021 | 0.7% | $`8.3\times10^{-1}`$ |

The shortfall falls without stopping, so no finite $`\rho`$ closes the ledger
exactly; and large $`\rho`$ costs fidelity, because $`E`$ becomes a small
difference of two large populations again and integration error is amplified
by $`N/\lVert E\rVert`$. **A ratio $`\rho \approx 3`$ to $`10`$ is the usable
window**:
1 to 2.5 per cent residual imbalance at essentially unchanged fidelity. In
world-particles, $`\rho = 3`$ is four positons and two negatons per
$`h`$-cell against the sea's two pairs — the entire ontological content of a
tunnelling event is a handful of bodies.

**The residual leak is accepted, and is not a defect of the implementation.**
The shortfall $`1/2 - f`$ never quite reaches zero, so the sea slowly loses
pairs where the minority species is scarce — 0.6 per cent globally over 24
time units, with a worst cell fluctuating around $`-1\thinspace B`$. That is
the price of postulate (S), and the following says the price cannot be
negotiated down.

**Proposition L1 (there is no exact-momentum third channel).** Let a batch of
demand $`D`$ be allocated as $`a`$ absorptive, $`e`$ emissive, and $`m_\pm`$
of the excluded $`M^\pm`$ realisations, with $`a = e`$ for ledger closure.

1. If $`m_+ = m_-`$, so that body momentum is conserved event by event, the
   partner demand on each leg is $`a + m_\pm/2 = D/2`$, **independent of
   $`m`$**. Momentum-balanced hops relieve the supply constraint of S6 by
   exactly nothing.
2. If $`m_+ \ne m_-`$, the feasibility condition becomes
   $`\mathrm{cap}_A + \mathrm{cap}_B \ge D`$ — a sum where S6 required a
   minimum, and with measured availabilities summing to $`\approx 1.45\,D`$
   this closes the ledger easily. The price is a net momentum injection
   $`2\xi_q(m_+ - m_-)`$ into the ensemble.

*Proof.* Substitute $`a = e = (D - m_+ - m_-)/2`$ into the leg demands
$`a + m_-`$ and $`a + m_+`$. In case 1 both equal $`D/2`$; in case 2 their sum
is $`D`$. $`\square`$

So the alternatives are a leaky ledger with every trajectory Newtonian, or a
closed ledger with body momentum conserved only in the mean. There is no exact
middle, and no amount of care in the allocation manufactures one.

**This specification takes the first, and does not implement the second
[normative].** Postulate (S) is the claim the compensated reading exists to
make: the entire content of
[`../analysis/eckart_barrier_compensated.md`](../analysis/eckart_barrier_compensated.md)
§8 is that tunnelling needs no world-particle to move discontinuously, and an
$`M^\pm`$ channel would give that back for a two per cent bookkeeping
improvement. There is no `hops` flag in the interface of §9 and no code path
for one. CLA9 records what the exploration measured, so that the option is not
reinvented from the ledger arithmetic alone.

### 5.8 What the implementation must expose

Cheap diagnostics do not see ledger failures — under the supply-throttled run
of S5 the norm is conserved to $`10^{-15}`$ and $`\langle p \rangle`$ to
$`10^{-7}`$ while $`W`$ is 40 per cent wrong in the core. So the following are
required, not optional **[normative]**: the absorptive fraction $`f`$; the
running counts $`n_{\rm abs}`$ and $`n_{\rm emi}`$; the ledger residuals
against S7; $`\min_x S/B`$ and $`\langle S\rangle/B`$; the starved-demand
fraction; and $`\min(N - \lvert E \rvert)`$.

One of these is a correctness check rather than a diagnostic. Because §5.3
admits only the emissive and absorptive realisations, **every event must change
$`N`$ by exactly $`\pm 2`$ and $`S`$ by exactly $`\mp 1`$**, and S7 is then an
identity rather than an approximation. An implementation whose ledger residual
is small but not at the floor has admitted a hop somewhere — most likely by
allocating the two legs of one event independently — and is running a different
ontology from the one specified. Rung 11 is the test; treat a failure as a bug,
not as noise.

---

## 6. Composition and step size

### 6.1 What the reordering buys

Strang, at macro step $`\Delta t`$:

```math
\Phi_{\mathrm{cl}}(\Delta t/2) \;\circ\;
\exp\bigl(\Delta t\thinspace M_{\mathrm{res}}\bigr) \;\circ\;
\Phi_{\mathrm{cl}}(\Delta t/2) ,
```

with $`\Phi_{\mathrm{cl}}`$ the classical flow of §3. Both this and the
uncompensated scheme split the **same** operator
$`M_{\mathrm{cl}} + M_{\mathrm{res}}`$ and both are second order; the
reordering shrinks the **constant** by the symbol ratio, because the surviving
commutator is with $`M_{\mathrm{res}}`$ rather than with the whole potential
term. Verified (Part F), cosine well, $`y_{\max} = 1`$ so $`u = \pi/4`$,
$`L^1`$ error at $`t = 1`$, hard horizon:

| $`\Delta t`$ | uncompensated | compensated | gain |
|---|---|---|---|
| 0.1250 | $`4.23\times10^{-1}`$ | $`4.40\times10^{-2}`$ | 9.6 |
| 0.0312 | $`3.50\times10^{-2}`$ | $`3.69\times10^{-3}`$ | 9.5 |
| 0.0156 | $`2.42\times10^{-3}`$ | $`2.51\times10^{-4}`$ | 9.7 |
| 0.0078 | $`2.28\times10^{-4}`$ | $`2.19\times10^{-5}`$ | 10.4 |

against a predicted
$`\lvert M_{\mathrm{res}}\rvert/\lvert M_{\mathrm{cl}}\rvert = 0.0904`$, i.e.
$`11.1\times`$, and the leading estimate $`u^2/6 = 0.103`$. The gain is flat
in $`\Delta t`$, which is the signature of a changed constant rather than a
changed order.

Under the specified soft horizon the residual is smaller again — the profile
attenuates the far rungs, where $`M_{\mathrm{res}}`$ is largest, giving a ratio
of $`0.0072`$ — and the measured gains run into the hundreds. Those figures are
lower bounds: the compensated error there sits close to the reference's own
floor, so the comparison stops being a measurement. The hard-horizon column is
the honest one to quote for the *mechanism*; the soft horizon simply has less
residual to split badly.

Be clear about what this is not. The scheme is still second order; the
Strang error between streaming and the jump channel is still there; and on a
mesh the classical flow costs about what the jump channel costs, so the gain is
in accuracy, not in time. In the world ensemble the gain is in both.

### 6.2 Choosing $`\Delta t`$

Three conditions, whichever binds:

- CFL on streaming: $`\Delta t \le m\thinspace\Delta x / P_{\max}`$.
- Event rate: $`\max_x R(x)\thinspace\Delta t \ll 1`$, with $`R`$ from §4.4.
  In the world form this is the mean events per world per step; in the mesh
  form the exact substep of §4.3 removes the condition, and only the Strang
  error remains. Note that $`R`$ grows with the reach without bound (G5), so
  the reach sets the timestep in the world form even though it does not in the
  mesh form.
- Ledger accuracy, world form only: the allocation of §5.5 is first order in
  $`\Delta t`$ and its error is amplified by $`N/\lVert E\rVert`$, so a run with
  a large standing population needs a smaller $`\Delta t`$ than the same run
  with $`\rho = 1`$ (§5.7). The absorptive fraction $`f`$ itself is *not*
  $`\Delta t`$-sensitive — $`0.4718, 0.4706, 0.4704`$ over a factor of four —
  so a converged $`f`$ is no evidence that $`\Delta t`$ is small enough.
- Splitting accuracy: from §6.1, and this is usually the loosest of the three
  once the reach is short.

---

## 7. Reference pseudocode

### 7.1 Setup (once, for static $`V`$)

```python
import numpy as np

dp   = np.pi * hbar / (2 * y_max)          # forced by the reach
p    = (np.arange(N_p) - N_p // 2) * dp
s    = 2 * np.pi * np.fft.fftfreq(N_p, d=dp)
y    = hbar * s / 2                        # the rungs
q    = np.round(np.fft.fftfreq(N_p, d=1.0) * N_p).astype(int)
xi   = q * dp                              # momentum transfer atoms
nyq  = N_p // 2

X, Y = np.meshgrid(x, y, indexing="ij")     # x is an OPEN-line window
M    = (1j / hbar) * (V(X + Y) - V(X - Y)) # sample V over each reach
S    = np.broadcast_to(s, M.shape).copy()
M[:, nyq] = 0.0                            # 4.1 [normative]
S[:, nyq] = 0.0

def first_moment(sym):                     # discrete first moment of a kernel
    return np.real((xi * np.fft.ifft(sym, axis=-1)).sum(axis=-1))

dV_eff = first_moment(M) / first_moment(1j * S)      # 3.2 [normative]
accel  = -dV_eff                           # what the flow of 3.1 integrates
M_cl   = 1j * dV_eff[:, None] * S
M_res  = M - M_cl

w      = np.cos(np.pi * y / (2 * y_max)) ** 2        # 4.4 horizon profile
M_res *= w[None, :]                        # the RESIDUAL only [normative]

K_res  = np.real(np.fft.ifft(M_res, axis=1))         # signed transfer rates
Gamma  = -K_res                                      # 4.2 orientation
rate   = np.abs(K_res).sum(axis=1)                   # R(x), see 4.4
active = rate > tol                                  # 7.4
```

### 7.2 The step

```python
def step(W, dt):
    W = classical_flow(W, dt / 2, accel)   # 3.1: a FLOW, not a split
    Wh = np.fft.fft(W[active], axis=1)     # 7.4: quiet cells skipped
    Wh *= np.exp(dt * M_res[active])       # 4.3, exact in dt
    W[active] = np.real(np.fft.ifft(Wh, axis=1))
    return classical_flow(W, dt / 2, accel)
```

The world-ensemble step replaces `classical_flow` by a per-world symplectic
integrator and the middle block by §7.3. It streams `u_plus`, `u_minus` and
`sea` alike — the sea is carried, not pinned.

### 7.3 The ledger step (world form)

The middle block of §7.2 replaced by the allocation rule of §5.5. `up`, `um`
and `sea` are $`(M_x, N_p)`$ arrays; `roll(f, k)` is `np.roll(f, k, axis=1)`.

```python
def step_events(up, um, sea, dt):
    """One potential substep in the world form.  Returns (n_abs, n_emi)."""
    n_abs = n_emi = 0.0
    for q in range(1, N_p // 2):               # (q, -q) is ONE event (5.2)
        lam = np.abs(K_res[:, q])[:, None]
        if lam.max() < 1e-14:                  # quiet channel (7.4)
            continue
        sg = np.sign(K_res[:, q])[:, None]
        for parent, sp in ((up, +1.0), (um, -1.0)):
            D = lam * parent * dt              # demand this channel makes
            t = sg * sp                        # sign deposited at p + xi_q
            # partners absorption needs: OPPOSITE species at each daughter,
            # read LIVE, not from a snapshot (5.5)
            capA = np.where(t > 0, roll(um, -q), roll(up, -q))   # at p + xi
            capB = np.where(t > 0, roll(up, +q), roll(um, +q))   # at p - xi
            A  = np.minimum(D, np.minimum(capA, capB))           # absorptive
            Em = D - A                                           # emissive
            n_abs += A.sum(); n_emi += Em.sum()

            aq, am = roll(A, q), roll(A, -q)   # absorptive: consume, then bind
            um -= np.where(t > 0, aq, 0.0); up -= np.where(t > 0, 0.0, aq)
            up -= np.where(t > 0, am, 0.0); um -= np.where(t > 0, 0.0, am)
            sea += A                           # credited at the PARENT row

            eq, em = roll(Em, q), roll(Em, -q)  # emissive: ionise, then place
            up += np.where(t > 0, eq, 0.0); um += np.where(t > 0, 0.0, eq)
            um += np.where(t > 0, em, 0.0); up += np.where(t > 0, 0.0, em)
            sea -= Em                          # debited at the PARENT row

            np.maximum(up, 0.0, out=up)        # 5.1 [normative]
            np.maximum(um, 0.0, out=um)
    return n_abs, n_emi
```

Both `sea` credits and debits land on the parent row, which is what removes
the transport asymmetry that sinks the emissive unravelling (S4). All three
fields are streamed by the classical flow of §3 — the sea is not static, it is
carried.

### 7.4 The active set

Theorem C7 is an optimisation as well as a theorem: where $`V'''`$ vanishes
over the whole reach the rate field is identically zero, and the jump channel —
and, in the world form, all the bookkeeping around it — can be skipped
outright.

Verified (Part G) on a harmonic trap plus a $`C^\infty`$ bump supported on
$`[-1,1]`$, over $`x \in [-6, 6]`$:

| reach $`y_{\max}`$ | 0.25 | 0.50 | 1.00 | 2.00 |
|---|---|---|---|---|
| predicted active region | $`\lvert x\rvert < 1.25`$ | $`< 1.50`$ | $`< 2.00`$ | $`< 3.00`$ |
| measured active fraction | 0.200 | 0.242 | 0.321 | 0.483 |

The measured edge tracks $`b + y_{\max}`$ with unit slope, falling silent two
to four position cells inside it because the bump's third derivative
approaches zero smoothly. **Mask on the measured rate field, not on an
analytic edge [normative]** — a potential without compact support has no exact
edge, and §5.1 of the note shows the interaction there is the potential's own
profile translated outward by one reach, not a cliff.

The crystal-lattice algorithm has no counterpart to this saving: its rate
field $`\Gamma_q(x) = -(V_q/\hbar)\sin(2\pi q x/L + \phi_q)`$ is non-zero at
almost every $`x`$ for every mode, and the cancellation that makes the
residual vanish is spread across the whole mode sum (note §6.2).

![Splitting error and the active set](https://raw.githubusercontent.com/billpage/wpmw/output/figures/compensated_algorithm_validation.png)

### 7.5 Cost

| | setup | per step |
|---|---|---|
| mesh, this algorithm | $`O(M_x N_p \log N_p)`$ | $`O(M_x^{\mathrm{act}} N_p \log N_p)`$ + flow |
| mesh, crystal lattice | $`O(Q M_x)`$ | $`O(Q M_x N_p)`$ |
| worlds, this algorithm | as above | $`O(N_w)`$ + events at rate $`R`$ |

$`M_x^{\mathrm{act}}`$ is the active-set size. The setup cost is paid once for
a static $`V`$ and repaid every step; for a time-dependent $`V`$ it is paid
every step and the comparison changes **[open]**.

---

## 8. Validation sequence

Rungs an implementation should pass, in order. Each is checked in the
companion demo at the section given.

1. **Split identity.** $`M - (M_{\mathrm{cl}} + M_{\mathrm{res}})`$ is
   identically zero on the grid. (Trivial, but catches Nyquist-handling
   mismatches between the two symbols.)
2. **Kernel structure.** $`K_{\mathrm{res}}`$ real and odd to $`10^{-16}`$;
   zeroth moment zero. §4.1.
3. **Momentum neutrality.** First moment of $`K_{\mathrm{res}}`$ zero to
   $`10^{-15}`$ at every $`N_p`$. §3.2.
4. **Stencil equals spectral.** The §4.2 rule reproduces
   $`\exp(\Delta t M_{\mathrm{res}})`$ to first order in $`\Delta t`$ and the
   generators agree to $`10^{-15}`$. §4.2.
5. **Crystal-lattice corner.** Uncompensated, single mode, maximal reach:
   two atoms, $`-K_1 = \Gamma_1`$. §4.2.
6. **Quiet region.** Harmonic trap plus compact bump: rate field below
   tolerance beyond $`b + y_{\max}`$, edge tracking the reach. §7.4.
7. **Harmonic exactness.** Open-line quadratic $`V`$: $`M_{\mathrm{res}}`$
   identically zero, so the algorithm reduces to the exact classical flow and
   the error is the flow integrator's alone. This is the sharpest single test;
   it is the one the crystal-lattice algorithm cannot pass, and it is also the
   test that catches a horizon profile applied in the wrong place (§4.1).
7b. **Open-line operation.** No wrap anywhere; the active set compact for a
   trap-plus-bump potential; widening the window changes nothing. §2.5, and
   Part I of the demo.
7c. **Sampler neutrality.** In the world form, the signed mean transfer over
   many drawn events is at the statistical floor — the sampler inherits C3.
   Measured: $`4.3\times10^{-3}`$ of the mean transfer scale over 3635 events,
   against a floor of $`1.7\times10^{-2}`$. Part I.
8. **Convergence against a reference.** Second order in $`\Delta t`$, with the
   constant below the uncompensated scheme's by the symbol ratio. §6.1.
9. **Agreement with `wigner_split_fourier`** on the QHO, and with the
   crystal-lattice solver on a cosine well at maximal reach.

Rungs 10 to 15 test the world form and the ledger. Nothing in 1 to 9 touches
them, and by Theorem S2 no observable does either.

10. **Quiet in the world form.** Quadratic $`V`$, world ensemble, a whole run:
    zero events, $`N`$ and $`S`$ unchanged to machine precision, every
    worldline a Newtonian arc. The mesh version of this is rung 7; the world
    version is the one that would catch a ledger that leaks on an *empty*
    channel.
11. **The ledger identity.** $`\Delta N - 2(1-2f)n_{\rm ev}`$ and
    $`\Delta S - (2f-1)n_{\rm ev}`$ at the floor. Measured
    $`2\times10^{-14}`$ and $`5\times10^{-12}`$. Fails loudly if a leg is
    applied to the wrong species or a sea credit lands on a daughter row.
12. **Absorptive against emissive.** Eckart barrier, $`T = 6`$,
    $`\Delta t = 0.02`$: $`N = 2.32`$ and rel $`L^2 = 2.0\times10^{-2}`$
    absorptive, against $`N = 1.0\times10^{5}`$ and $`1.8\times10^{2}`$
    emissive. Both integrate the same $`E`$; the difference is amplification
    by $`N/\lVert E\rVert`$. The error must **halve** with $`\Delta t`$ — if it
    saturates instead, something is throttling the rate and §5.1's promise is
    broken.
13. **The attractor.** $`f`$ from $`\rho = 1`$ and $`\rho = 20`$ agreeing to
    $`0.015`$ by $`t = 8`$, approached from below and above respectively.
14. **Open-line transmission.** Eckart, $`\beta = 2`$: $`\Sigma`$ at
    $`t = 0`$ gives the classical $`0.500`$, the run gives $`0.5445`$ against
    a closed-form $`0.5442`$, and the C1 check
    $`\lvert\text{full} - \text{compensated}\rvert = 4\times10^{-6}`$. This is
    the only rung with an exact answer from outside the project.
15. **Admissibility.** Quartic $`\lambda = 0.05`$ from
    $`(|0\rangle + |2\rangle)/\sqrt2`$: the least eigenvalue of the
    reconstructed $`\rho`$ stays at the grid floor under (S)+(D) and reaches
    $`-0.10`$ under (S) alone (G4). Switching off the jump channel must
    *visibly* break this, or the channel is not being applied.

---

## 9. The `wpmwlib` interface

Sketch for the module this specification is meant to become, following the
conventions of `wpmwlib/phase_space_crystal_lattice.py` — array shape
$`(M_x, N_p)`$, position axis first **[choice]**; note this transposes that
module's convention and the two must not be mixed silently.

```python
class CompensatedLiouville:
    """Wigner evolution as Newtonian flow plus a zero-mean signed jump channel.

    Parameters
    ----------
    y_max : float          reach; sets dp = pi*hbar/(2*y_max).  NOT a tolerance.
    N_p   : int            rung count (even); sets P_max = N_p*dp/2.
    M_x   : int            position cells.
    box   : float          position window.
    V, dV : callable       potential and its derivative (dV for diagnostics only).
    mass, hbar : float
    flow  : {'semi_lagrangian', 'symplectic'}
    horizon : {'raised_cosine', 'hard'}    see 4.4; soft is the default, and
                                          'hard' (w = 1) is the special case.
    """

    # --- setup -------------------------------------------------------
    def build_rate_field(self) -> None: ...     # 4.1; sets .Gamma, .rate, .active
    @property
    def accel(self) -> np.ndarray: ...          # 3.2, the kernel's first moment
    @property
    def dp(self) -> float: ...                  # forced; read-only

    # --- substeps ----------------------------------------------------
    def step_classical(self, dt) -> None: ...   # 3.1
    def step_residual(self, dt) -> None: ...    # 4.3, exact
    def step(self, dt) -> None: ...             # 6.1, Strang

    # --- diagnostics -------------------------------------------------
    def moments(self) -> tuple: ...             # zeroth and first, per x
    def event_budget(self) -> np.ndarray: ...   # R(x); meaningless without N_p
    def negativity(self) -> float: ...
```

Four interface requirements follow from the body of this document.
`dp` is derived and read-only, so that no caller can refine it (§2.3).
`accel` is derived from the kernel, not from `dV` (§3.2). `build_rate_field`
zeroes the Nyquist rung before anything else touches the symbol, and applies
the horizon profile to the residual only, after the force has been extracted
(§4.1). And there is no `L` or `periodic` parameter: `box` is a window on the
open line, not a circumference (§2.5).

The world form is a **second class**, not a mode flag on the first, because it
carries state the mesh form has no analogue of (§5.1):

```python
class SignedLedger:
    """World-form state: two body fields and a per-cell sea counter (5.1).

    Parameters
    ----------
    cla  : CompensatedLiouville   supplies K_res, the reach, and the grid
    rho  : float                  standing population N_0/|E| (5.7).  1 is the
                                  minimal ensemble; 3 to 10 is the usable
                                  window on a barrier; large rho costs
                                  fidelity, not bodies.
    kappa: float or None          optional bilinear recombination sink; must
                                  scale as Gamma_tot(x) if used at all (5.4).

    There is deliberately no option to enable the M+/M- realisations of 5.3.
    They would close the ledger and they move a body between momentum rows,
    which postulate (S) forbids; Proposition L1 shows there is no variant that
    does the first without the second.
    """

    # --- state -------------------------------------------------------
    u_plus:  np.ndarray     # (M_x, N_p), >= 0
    u_minus: np.ndarray     # (M_x, N_p), >= 0
    sea:     np.ndarray     # (M_x, N_p) scalar counter, background B = 2/h

    @property
    def E(self) -> np.ndarray: ...   # u_plus - u_minus; this IS W
    @property
    def N(self) -> np.ndarray: ...   # u_plus + u_minus; >= |E| [normative]

    # --- substeps ----------------------------------------------------
    def step_classical(self, dt) -> None: ...   # 3.1, on ALL THREE fields
    def step_events(self, dt) -> tuple: ...     # 7.3; returns (n_abs, n_emi)
    def step(self, dt) -> tuple: ...            # 6.1, Strang

    # --- diagnostics, all required by 5.8 ----------------------------
    def absorptive_fraction(self) -> float: ... # f; the S7/S9 number
    def ledger_residual(self) -> tuple: ...     # rung 11
    def sea_min_fraction(self) -> float: ...    # min S/B; < 0 is a deficit
    def unpaired_total(self) -> float: ...      # N integrated
    def starved_fraction(self) -> float: ...    # demand absorption could not
                                                # meet, per 5.4
```

Four requirements on this class. `sea` is a **counter array, never a particle
list** (§5.6). `step_classical` streams all three fields, since the sea is
carried by the flow and not pinned. `E` is a derived property with no setter: a
caller who writes to it has separated the observable from the population it is
supposed to be a difference of, which is the one way to break Theorem S2. And
`step_events` writes both legs of an event in one operation, never leg by leg,
because a per-leg write is by S6 an $`M^\pm`$ hop in disguise.

The diagnostic names `sea_min_fraction` and `unpaired_total` match
`wpmwlib/sea_dressed_lattice.py` deliberately; the two classes carry the same
three fields for different generators and should be comparable without
translation.

---

## 10. Open implementation questions

- **CLA1 (what sets the profile).** §4.4 makes the horizon's *shape* a free
  function of the model, not just its width, and shows that the sharp
  Definition (H) is the one shape under which the event rate, the momentum
  churn and the third moment all fail to exist. Two questions follow. What
  physically sets $`w`$ — is it a decoherence profile, a rung-occupancy
  statistic of the sea, or something else? And is there a distinguished
  profile, e.g. the one minimising $`\sum_q\lvert\xi_q K_q\rvert`$ at fixed
  $`w(0) = 1`$? This is CLS4 of the note ("what sets $`L_c`$?") with more
  structure to bite on, and it is now the load-bearing open item: the
  specification's default is a **[choice]** where a derivation is wanted.
- **CLA1b (Definition (H)).** If CLA1 resolves in favour of a profile, the
  coherence horizon should be restated in
  [`../analysis/open_position_space.md`](../analysis/open_position_space.md)
  §3 as a rung occupancy rather than a cutoff. Theorems O1–O5 need re-reading
  against that: O5's coset invariant is unaffected (it depends on the transfer
  lattice, not on occupancy), but Proposition O4's localisation edge becomes a
  profile edge, and C7 becomes exact only where $`w`$ has compact support —
  which the raised cosine does.
- **CLA2 (world-form validation).** Everything in §6.1 is measured on a mesh,
  and so is everything in §5 — the ledger runs are mean-field tau-leaps, not
  ensembles. The claimed cost advantage is a world-ensemble claim and is still
  unmeasured. The decisive test is a signed-ensemble run against the mesh
  reference at fixed variance, and it now has a prediction attached (S-SP5):
  the absorptive ensemble should show the emissive ensemble's variance reduced
  by roughly $`N_{\rm em}/N_{\rm abs}`$, which was $`4\times10^{4}`$ in the
  §5.7 runs.
- **CLA3 (annihilation load). Resolved, and the answer was not the expected
  one.** The question was how much annihilation the unravelling costs. The
  answer is that no separate annihilation rate is needed: absorption is itself
  a pair-removal channel satisfying S1, and by S9 the absorptive fraction
  regulates itself onto the closing value $`f = 1/2`$ from any preparation.
  What survives is smaller and different. The bilinear sink's one
  dimensionless constant is still unfixed — and, worse for pinning it, closure
  happens at *every* initial ratio, so closure carries no information about
  the constant (S-SP1). The pricing in
  [`../analysis/species_sectors_and_annihilation.md`](../analysis/species_sectors_and_annihilation.md)
  against nodal structure and this note's pricing against event rate have
  still not been compared on one problem.
- **CLA7 (the reach is over-determined).** §2.4. K1 caps $`y_{\max}`$ at the
  analyticity strip, $`\pi a/2`$ for the Eckart barrier; the ledger of §5 wants
  $`\approx 4\pi a`$, eight times larger; and Corollary K1.1 says a
  reach-limited lattice cannot resolve a packet of width $`\sigma_r \ge a/2`$
  at all. All three cannot hold. Either the ledger's reach preference is an
  artefact of the mean-field tau-leap, or K1's ceiling constrains only the
  Moyal *reading* and not the algorithm, or the Eckart barrier is a bad test
  case for the reach. Not decided, and it blocks a single recommended
  parameter set.
- **CLA8 (the box).** §5.6. The sea is uniform at two pairs per $`h`$-cell, so
  its total scales with the phase-space volume of the window while the
  unpaired population does not — 40 to 140 to one in the Eckart runs, and more
  if the window is widened. On a closed system that ratio would be physics; on
  the open line it is a boundary choice. The presumed resolution is that only
  the sea within the packet's support ever participates, but nothing measured
  establishes it, and the ontology of §1.3 needs it to be true.
- **CLA9 (the deprecated $`M^\pm`$ hop channel).** Not an open question; a
  record, so that the option is not reinvented from the ledger arithmetic.
  Admitting the hop realisations of §5.3 at lowest priority does close the
  ledger and does keep every cell out of sea deficit. Two things measured
  during that exploration are worth keeping. The minimum hop fraction that
  closes the ledger is exactly twice the two-channel shortfall,
  $`2(1/2 - f)`$ —
  $`0.123, 0.063, 0.025`$ at $`y_{\max}/a = 2\pi, 4\pi, 8\pi`$ — so the damage
  is a quantity already tabulated in §5.4. And the allocation has a free
  local/global balance constant $`c`$. Greedy absorption, $`c = 1`$, eats the
  partner supply and *increases* the hop count, while $`c = 1/2`$ measured
  best for QLE fidelity. None of this is in an analysis note, and under the
  present specification none of it is implemented. If (S) is ever revisited,
  start here.
- **CLA10 (does the census mean anything).** This is G-SP1, restated as an
  implementation question. $`\Gamma`$ grows without bound with the reach while
  the generator converges to $`10^{-14}`$, so the number of events — and hence
  the number of world-particles — is a property of the regulator. Until split
  gauge invariance is proved, $`f`$, $`N`$ and $`S`$ are quantities of the
  representation and not of the physics, and no result in §5 should be quoted
  without its reach.
- **CLA4 (time-dependent $`V`$).** The whole rate field is rebuilt each step,
  and §7.5's cost comparison inverts. Whether the active set still pays is
  unknown.
- **CLA5 (the sub-$`\Delta p`$ acceleration).** This is CLS3 of the note,
  partly answered. Step 1 moves momentum by an amount that is generically not
  a multiple of $`\Delta p`$, so it cannot be represented on the momentum
  lattice at all; §3.3's two forms both evade this by keeping momentum
  continuous within step 1 (semi-Lagrangian, or per-world). Whether a strictly
  lattice-resident version exists — and hence whether Theorem O5's coset
  invariant survives the compensated split as an exact statement — is open.
- **CLA6 (multi-dimensional).** The symbol construction is per-degree-of-
  freedom and the rung grid becomes a lattice in $`\mathbf{y}`$, so the setup
  cost scales as $`(M_x N_p)^d`$. The active set may be the only thing that
  makes $`d = 2`$ affordable; untested.

---

## 11. Sources

- [`../analysis/compensated_liouville_splitting.md`](../analysis/compensated_liouville_splitting.md) — theorems C0–C7, and the physical argument this specification implements.
- [`../analysis/eckart_barrier_compensated.md`](../analysis/eckart_barrier_compensated.md) — K1 (the reach ceiling, §2.4), K4–K6 (the conserved outcome functional and the $`1/\beta`$ cost law, §3.4), K8 (co-location, §5.2), K9 and §8.4 (the supply condition and the standing population, §5.7).
- [`../analysis/sea_population_equilibrium.md`](../analysis/sea_population_equilibrium.md) — S1–S9, the whole of §5, and the erratum to §4.3.
- [`../analysis/compensated_ontology.md`](../analysis/compensated_ontology.md) — postulates (E), (A), (S), (D) and theorems G1–G5; §1.3, and validation rungs 10 and 15.
- [`../analysis/reach_energy_coupling.md`](../analysis/reach_energy_coupling.md) — Theorem E1 ($`\Delta p = \pi\hbar/2y_{\max}`$, §2.2) and the compensated cancellation of the $`1/y_{\max}^2`$ contamination.
- [`phase_space_crystal_lattice_algorithm.md`](phase_space_crystal_lattice_algorithm.md) — §3b/§3c, the jump rule reused verbatim in §4.2.
- [`../analysis/open_position_space.md`](../analysis/open_position_space.md) — Definitions (H) and (R), Theorem O2 (the first moment is the force), Theorem O5 (the coset invariant).
- [`../supplement/takabayasi_1954_stochastic_picture.md`](../supplement/takabayasi_1954_stochastic_picture.md) — Proposition T3, on why the signed kernel admits no one-body Markov unraveling.
- [`../supplement/emission_and_absorption.md`](../supplement/emission_and_absorption.md) — tutorial treatment of §5.3 and §5.4; Theorem J1 (the pair count $`P = S + N/2`$ is conserved by the event channel at any absorptive fraction) and Theorem J2 (that identity is global only, and does not constrain the local dynamics). Open item J-SP2 proposes a change to the allocation loop of §5.5.
- [`../analysis/species_sectors_and_annihilation.md`](../analysis/species_sectors_and_annihilation.md) — the positon/negaton unraveling of §4.3, and the $`c = 0`$ dark pair of §5.2.
- David Cyganski, *Extended Fokker–Planck Eq. and the QLE V2* (project memo).
