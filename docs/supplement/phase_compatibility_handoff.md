# Phase in the compensated Liouville model: a hand-off

**Scoping document for an analysis note that has not been written. States the
question, what is already established on each side of it, the one measurement
that would settle it, and what fails if it cannot be settled. Written to be
the first thing read in a fresh conversation.**

---

## 0. What this is

Not a results note. A hand-off, written because the question below needs a
conversation of its own and because the material bearing on it is scattered
across six documents that were written for other purposes.

**Read this first, then §2's reading list, then nothing else until the
question in §1 has an answer.**

Nothing here is new. Everything is cited to a note that establishes it, except
where §4 marks a claim as unestablished.

---

## 1. The question

The compensated Liouville model and the phase-resonance line of work both have
a process called **absorption**, and they are not obviously the same process.

**Ledger absorption** (§4 of
[`emission_and_absorption.md`](emission_and_absorption.md), Theorem S6 of
[`../analysis/sea_population_equilibrium.md`](../analysis/sea_population_equilibrium.md)).
An event at parent row $`p`$ with channel $`q`$ consumes a negaton from row
$`p + \xi_q`$ and a positon from row $`p - \xi_q`$ and binds them into a
neutral pair at row $`p`$. The parent streams on undisturbed. **No body
changes momentum.**

**Phase-resonance absorption**, the K3 vertex (§4 of
[`../analysis/phase_resonance_microdynamics.md`](../analysis/phase_resonance_microdynamics.md)).
The excess particle's in-leg at $`p_{hi}`$ meets the excited partner of a
phase-matched beating pair; the particle exits at
$`p_{lo} = p_{hi} - \sigma\,2q\,\delta p`$, the struck partner moves to its
mate's momentum with gauge-matched phase, and the pair exits dark. **The
particle changes momentum.**

> **The question.** Does the phase machinery — the particle-level phase of
> P0–P5, the misalignment $`\mu`$, and darkness as interference — transplant
> onto the mediated-jump vertex, in which no body changes momentum? If it
> does, write the note. If it does not, the compensated model needs a
> different account of why the sea is dark, and §5 says what breaks.

---

## 2. Reading list, in order

1. §8.1 of
   [`../analysis/eckart_barrier_compensated.md`](../analysis/eckart_barrier_compensated.md)
   — "A jump is not a hop", the postulate (S) argument, and Proposition K8.
   Four pages. This is the constraint.
2. §§1–4 of
   [`../analysis/phase_resonance_microdynamics.md`](../analysis/phase_resonance_microdynamics.md)
   — phase as a particle property (P0–P5), darkness by interference, the
   vertices. This is the machinery.
3. §2 of
   [`../analysis/phase_alignment_microdynamics.md`](../analysis/phase_alignment_microdynamics.md)
   — the misalignment $`\mu`$, Lemma 4 (completeness), Proposition 3 (winding
   rates). This is the machinery in its cleanest form and is the version to
   work with.
4. [`../analysis/four_rule_microdynamics_equivalence.md`](../analysis/four_rule_microdynamics_equivalence.md)
   — the exact equivalence of the hop and jump formulations, and the
   $`G`$-freedom. §5 below is about the limits of that equivalence.
5. §14 of [`emission_and_absorption.md`](emission_and_absorption.md) — the
   bootstrap, and why the creation rule matters.
6. §1 of
   [`what_the_reach_is.md`](what_the_reach_is.md) — the two objects called
   $`\mu`$, since confusing them is the first thing that will happen.

---

## 3. What is established

### 3.1 On the compensated side

- **Postulate (S).** No world-particle ever changes its momentum
  discontinuously. Every body — excess positons and every member of the
  created sea alike — follows a genuine Newtonian worldline under the full
  classical force, for its entire life. Eckart §8.1.
- **Proposition K8.** If the event conserves body momentum and the parent
  streams on undisturbed, the pair it produces must come from a bound pair at
  the parent's own momentum row. Co-location is derived, not postulated.
- **Theorem N2** of
  [`../analysis/stochastic_ledger.md`](../analysis/stochastic_ledger.md). Both
  realisations of an event move $`E`$ identically, so the choice between them
  lies in the kernel of the observable map. The ledger carries a noise channel
  the Wigner function cannot see.
- **Theorem N1.** $`P = S + N/2`$ is conserved pathwise, in exact integers.
- A bound pair contributes zero to $`E`$ *arithmetically*, its two members
  carrying opposite sign. No phase is invoked anywhere in the compensated
  formulation — $`W`$ is real and there is nothing in the algorithm that looks
  like a phase.

### 3.2 On the phase side

- **P0–P5.** Every world-particle, sea or excess, carries a de Broglie phase
  $`\theta`$ winding along its worldline at the Lagrangian rate, and
  broadcasts it as a plane wave riding the worldline; the negaton broadcasts
  with an intrinsic sign flip.
- **Lemma 4 (completeness).** A pair's entire gauge-invariant relational
  content is $`\mu = \Phi_a - \Phi_b`$. The pair state is $`(p_a, p_b, \mu)`$
  and nothing else, and $`|\Psi| = 2|\sin(\mu/2)|`$.
- **Darkness is $`\mu \equiv 0`$**, by interference rather than by bookkeeping
  subtraction.
- **Proposition 3.** $`\partial\mu/\partial x = \Delta p/\hbar`$, so a
  co-located, co-moving pair has $`\Delta p = 0`$, its misalignment never
  winds, and darkness set at creation is preserved for all time under
  classical streaming. This is the stability result the creation rule needs.
- **Theorem 2 (no-go).** No phase-blind microdynamics can be linear in
  $`V_q`$. The phase variable is *necessary*, not decorative.

### 3.3 The numerical check worth having in hand

The two quantities called $`\mu`$ are structurally the same object — a
gauge-invariant relative phase of a two-ended thing, whose gradient is the
momentum conjugate to the ends' separation — and numerically different:

| | ladder $`\mu`$ | pair $`\mu`$ |
|---|---|---|
| between | two legs of a ket–bra pair, one body | two members of a sea pair, two bodies |
| winds at | $`\partial\mu/\partial X = p/\hbar`$ | $`\partial\mu/\partial x = \Delta p/\hbar`$ |
| $`\mu = 0`$ | no current: *still* | $`|\Psi| = 0`$: *dark* |

Verified in Part A of `src/demo_reach_tutorial.py`. Keeping them apart is a
precondition for thinking clearly about the question.

---

## 4. What is not established

- **That the two absorptions are the same event.** §1. This is the question.
- **That equivalence of generators settles it.** It does not, and this is the
  crux.
  [`../analysis/four_rule_microdynamics_equivalence.md`](../analysis/four_rule_microdynamics_equivalence.md)
  proves the hop and jump formulations exactly equivalent at any particle
  number $`\nu`$. But that is equivalence of the *generator*, hence of $`E`$,
  and the property at issue — that worldline identity survives the event — is
  not a property of the generator. Theorem N2 is the same structure seen from
  the ledger side: two realisations agreeing exactly on the observable and
  disagreeing on the inventory. **Equivalence in $`E`$ is precisely what
  cannot settle an identity question.** Anyone starting this work should be
  clear about that before reaching for the equivalence theorem, because it
  looks like it answers the question and does not.
- **Whether a mediated-jump vertex can carry a phase rule at all.** In the
  ledger vertex the two consumed bodies come from *different* momentum rows
  and the pair is created at a *third*. Their phases $`\Phi_a`$ and $`\Phi_b`$
  are evaluated at different momenta, so whatever $`\mu`$ they had is not
  obviously the $`\mu`$ of the created pair. Either the created pair's
  misalignment is *forced* by the consumed bodies' phases — in which case
  darkness at creation is a claim to be checked and may simply be false — or
  it is *free*, in which case the creation rule is a new postulate and should
  be counted as one. The project currently has it both ways.
- **Whether the inventories match.** A single world-particle at $`(x,p)`$ maps
  under the $`y`$-transform to $`e^{2ipy/\hbar}`$, spread across every rung of
  the position-pair ladder; a single rung maps to something spread across every
  momentum. So the two pictures share the $`\mu`$-structure but not the
  inventory, and neither set is a coarse-graining of the other. Open item
  R-SP8 of [`what_the_reach_is.md`](what_the_reach_is.md).

---

## 5. What breaks if the answer is no

Not catastrophic, but not nothing.

- **The bootstrap of §14 loses its justification.** "Absorption creates dark
  pairs" becomes a free postulate rather than a consequence of a vertex, and
  the pleasant argument that the sea is paired *by construction* weakens to a
  stipulation.
- **Darkness falls back on arithmetic.** A bound pair contributes zero to
  $`E`$ because its members carry opposite sign. That is true, sufficient for
  every numerical result the project has, and ontologically thin: it says the
  sea is invisible by cancellation in a ledger, not that anything about the
  two bodies makes them invisible. §1.2 of the phase-resonance note exists
  because that thinness was felt to be a defect.
- **Theorem 2 becomes an argument against the compensated formulation.** If
  no phase-blind microdynamics can be linear in $`V_q`$, and the compensated
  jump formulation is phase-blind, then either the no-go's hypotheses do not
  cover it — most likely, since it is not a rate table — or something is
  wrong. Checking which is the sharpest single test available and should
  probably be the new note's §1.

---

## 6. Suggested shape for the note

Working title: *Phase in the compensated model: what the jump vertex can
carry.*

1. The two vertices, side by side, in one table. §1 above.
2. Whether Theorem 2's no-go reaches the mediated-jump formulation. The
   sharpest test; do it first, because a negative answer here changes
   everything downstream.
3. Transport: is $`\mu`$ well defined for a pair whose members sit on
   different momentum rows between events? Proposition 3 assumes a pair; the
   ledger vertex's two consumed bodies are not one.
4. The creation rule: forced or free. §4 above.
5. If forced and consistent — the derivation of darkness for the ledger
   vertex, which is the result the note exists to produce.
6. If not — what replaces it, and whether §5's fallback is acceptable.
7. Consequences for J-SP3 and J-SP4 of
   [`emission_and_absorption.md`](emission_and_absorption.md), and for R-SP8
   of [`what_the_reach_is.md`](what_the_reach_is.md).

A companion demo is likely to be small: the phase content is analytic and the
numbers in §3.3 are already computed. The expensive part is the reading.

---

## 7. Sources

- de Broglie, L. — *Recherches sur la théorie des quanta*, Ann. Phys. **3**
  (1925) 22–128. The phase that winds along the worldline, which is what
  P0–P5 promote to a particle property.
- Bloch, F. — *Nuclear Induction*, Phys. Rev. **70** (1946) 460–474; and
  Fano, U. — *Description of States in Quantum Mechanics by Density Matrix and
  Operator Techniques*, Rev. Mod. Phys. **29** (1957) 74–93. Populations and
  coherences, the standard vocabulary for the diagonal/off-diagonal
  distinction that the ladder $`\mu`$ lives on.
- Madelung, E. — *Quantentheorie in hydrodynamischer Form*, Z. Phys. **40**
  (1927) 322–326. The probability current, whose vanishing is what the ladder
  $`\mu = 0`$ case asserts.
- Takabayasi, T. — *The formulation of quantum mechanics in terms of ensemble
  in phase space*, Prog. Theor. Phys. **11** (1954) 341–373. Already in
  [`../../references/bibliography.md`](../../references/bibliography.md) and
  discussed in
  [`takabayasi_1954_stochastic_picture.md`](takabayasi_1954_stochastic_picture.md);
  the closest historical precedent for attaching phase to a phase-space
  ensemble, and worth re-reading with this question in mind.
- The project notes listed in §2.
