# Darkness and identity: what the compensated vertex can carry

> Whether the phase machinery transplants from the phase-resonance formulation onto the compensated one, and what the sea's darkness is a statement about. Theorem Y1 answers the sharpest test first, negatively for the no-go: the compensated kernel is read from the potential directly, `K(lam V) = lam K(V)` to machine precision, so the phase-blind no-go (Theorem 2) does not reach it and the phase it needs is not absent but already integrated into the kernel as the ladder `mu`. By Corollary Y1.1 a particle-level phase can therefore never be *necessary* for `E`, and any phase rule lives in the kernel of the observable map — where, by Theorem N2, the realisation choice already lives. Theorem Y2 asks what the sea's ineligibility for recombination is worth, since all four ways of settling the two legs move `E` identically: with the sea eligible there are realisations with `dN = dS = 0` whose effect on the counts is exactly the excluded hop, and which conserve momentum by shifting one aligned pair a row, while `P = S + N/2` stays conserved throughout — so ineligibility, not arithmetic, is what gives postulate (S) and Proposition K8 their content, and with it every event carries `|dS| = 1`. Theorems Y3 and Y4 supply the phase account darkness was missing: the two bodies a catalysed recombination consumes are a winding pair whose midpoint is the parent's own row, with `|Psi| = 2|sin(mu/2)|` and an envelope drifting at the parent's velocity, and a momentum kink of `+-xi` carries a phase ramp pinned at the point where it happens, so the created pair's misalignment freezes at `mu(x_k)` and is dark exactly when the kink sits at a node — a node-located kink and a phase reset at the parent's position being one freedom seen twice, not two rules. **Theorem Y5 is the cost.** If world-particles are conserved and none ever changes momentum at an event, per-row species counts are event-invariant and the residual channel can do nothing; so the two-clock reading forces piecewise worldlines, and the alternative is birth and death. Theorem Y6 prices that choice with tagged positons at the Eckart barrier: the individual momentum walk is driftless (`sum_q xi_q K_q = 0`) with diffusion `D_p = (1/2) sum_q xi_q^2 |K_q|`, an *unsigned* moment `E` never feels; `Gamma(0) = 0`, so no kink ever happens at the summit and an individual crosses by flank activation and classical passage, never through; and the measured crossing of tagged positons — 0.327 against 0.213 for `E` and 0.069 classically at `E_0 = V_0/2` — tracks neither the classical nor the quantum transmission, falling *below* both above the barrier, and is insensitive to the sea's momentum profile. Tunnelling therefore survives as a demographic account of `E` and not as a statement about identity, which corrects `eckart_barrier_compensated.md` §8.3 and the ledger note's "the diffusion in the phase variables is identically zero". Section 8 records a defect found along the way: the population clamp in `demo_emission_and_absorption.py` costs 16 per cent of `E` in 250 steps and biases `f`.

*Ladder abstract — see the [full list](README.md#the-ladder).*

**Status.** Analysis note, step 20 of the ladder. Companion demo:
`src/demo_dark_sea_and_identity.py` (`--full` for the barrier table).
Prompted by the standing question of whether the two processes the project
calls absorption are the same event — open items J-SP3 and J-SP4 of
[`../supplement/emission_and_absorption.md`](../supplement/emission_and_absorption.md)
§12 — and by the thinness of the arithmetic account of darkness noted in
§1.2 of [`phase_resonance_microdynamics.md`](phase_resonance_microdynamics.md).

---

## 0. What this note settles, and what it corrects

The project has two accounts of why the sea is invisible. In the compensated
formulation an aligned pair contributes zero to $`E`$ because its members
carry opposite sign: darkness is arithmetic. In the phase formulation a pair
is dark because its two transported clocks are aligned, $`\mu \equiv 0`$, so
its broadcast amplitude vanishes: darkness is interference. The first is
sufficient for every numerical result the project has, and says nothing about
the two bodies; the second says something about them but was built on a vertex
in which a body changes momentum.

**Settles.**

- **Whether the no-go reaches the compensated model.** It does not (§2). The
  compensated vertex consults the potential itself, which is the exception
  Theorem 2 names. The phase the no-go demands is present, but as the ladder
  $`\mu`$ integrated into the kernel, not as an attribute of a body.
- **What darkness has to be a statement about.** Not $`E`$ (§3). Every way of
  settling a recombination's two legs moves $`E`$ identically, so the sea's
  exclusion is a claim about identity, and it is what keeps the excluded hop
  excluded.
- **That the phase machinery does transplant** (§4, §5), with one free datum:
  where the kink happens.
- **What the two-clock reading costs** (§6, §7): piecewise worldlines, and
  with them an individual-level barrier crossing that no observable sees.

**Corrects.**

- **The framing of the question.** It is not true that the phase-resonance
  vertex bends a worldline and the compensated one does not. In the
  phase-alignment form the vertex is a swap, and §3 of
  [`../supplement/phase_alignment_interaction_diagrams.md`](../supplement/phase_alignment_interaction_diagrams.md)
  already reads a swap as two straight lines with permuted labels. The
  difference between the two vertices is elsewhere (§1).
- **The status of dark creation in the phase formulation.** §6 of
  [`phase_alignment_microdynamics.md`](phase_alignment_microdynamics.md) fixes
  the struck partner's exit phase *by requiring* the pair to leave aligned.
  That is a creation rule of the same standing as §14.3 of the emission
  supplement, not a consequence.
- **[`eckart_barrier_compensated.md`](eckart_barrier_compensated.md) §8.3.**
  "A world does not pass through the barrier" survives in the narrow sense —
  no body is ever in a classically forbidden place — but under piecewise
  worldlines individual bodies cross the separatrix in both directions at a
  rate set by $`D_p`$ (§7).
- **[`stochastic_ledger.md`](stochastic_ledger.md) §1.** "The diffusion in the
  phase variables is identically zero" holds under birth and death. Under
  piecewise worldlines there is a momentum diffusion, still not Nelson's
  (which acts in position) and still invisible to $`W`$.
- **`demo_emission_and_absorption.py`.** Its population clamp breaks $`E`$
  (§8).

**Inherits.** Postulate (S) and Proposition K8 from the Eckart note;
Theorems N1–N3 from the ledger note; Lemma 1, P0–P5 and Theorem 2 from the
phase-resonance note; Lemma 4 and Proposition 3 from the phase-alignment note;
the reach from [`open_position_space.md`](open_position_space.md).

---

## 1. The two vertices, side by side

**Catalysed recombination** (the absorptive realisation; §4 of the emission
supplement, Theorem S6). An event at parent row $`p`$, channel $`q`$, consumes
a negaton from row $`p + \xi_q`$ and a positon from row $`p - \xi_q`$ and
makes an aligned pair at row $`p`$. The parent streams on.

**The K3 vertex** (§4 of the phase-resonance note; Theorem 4 of the
phase-alignment note). An excess particle and the struck partner of a pair
exchange momenta, and the pair leaves aligned.

| | K3 as a swap | catalysed recombination |
|---|---|---|
| momenta in | $`\{p_b, p_a, p_b\}`$ | $`\{p-\xi,\ p+\xi,\ p\}`$ |
| momenta out | $`\{p_a, p_b, p_b\}`$ | $`\{p,\ p,\ p\}`$ |
| multiset | preserved | not preserved |
| $`\sum p^2/2m`$ | conserved | loses $`\xi^2/m`$ |
| worldlines | all kept, labels permuted | two end, one begins |
| Theorem 4's stationarity | imposed | holds identically |
| Theorem 4's momentum condition | met by transfer | met by collapse to the midpoint |

So the honest difference is not that one bends a worldline. It is that the
swap is elastic and relabels, while the ledger vertex is perfectly inelastic —
relative momentum $`2\xi`$ goes to zero — and relies on bodies beginning and
ending. The $`\xi^2/m`$ is where a plasma's binding energy would sit; here
nothing receives it and no observable records it. It is a second resident of
the kernel of the observable map.

---

## 2. Theorem Y1: the no-go does not reach the compensated kernel

**Theorem Y1.** *The compensated residual kernel is a linear functional of the
potential, evaluated by the vertex itself: $`K(\lambda V) = \lambda K(V)`$,
and the event rate $`\Gamma(x) = \sum_q |K_q(x)|`$ is first order in $`V`$.
Hence the hypothesis of Theorem 2 — that the potential reaches the rates only
through phase-blind local data, which Lemma 3 shows are independent of $`V`$
at first order — is false here, and the no-go does not apply.*

*Verification.* Part A of the demo, at $`\lambda = 10^{-3}, 0.1, 2, -1`$:
the relative departure from exact linearity is at machine precision, and
$`\Gamma_{\max}`$ scales as $`|\lambda|`$.

This is the exception the no-go itself names: a model that inserts the rate
field by hand rather than deriving it from the sea's state. The phase has not
been dispensed with. It has been integrated out. The residual symbol is a
transform over $`y`$ of $`V(x+y) - V(x-y) - 2yV'(x)`$, and $`y`$ is the
ket–bra half-separation — the *ladder* $`\mu`$ of
[`position_pair_ladder.md`](position_pair_ladder.md), not the *pair* $`\mu`$
at issue here. The compensated model is phase-blind in its bodies and not in
its kernel.

**Corollary Y1.1.** *A particle-level phase cannot be necessary for $`E`$ in
the compensated formulation. Any phase rule it adopts therefore lies in the
kernel of the observable map — the same place Theorem N2 puts the choice
between the two realisations.*

That is not a demotion. The question this note asks — what makes the sea dark,
and what stays the same body through an event — lives in that kernel by
construction.

---

## 3. Theorem Y2: what the sea's ineligibility buys

Take one cell, rows $`p-\xi`$, $`p`$, $`p+\xi`$, and settle each leg of a
catalysed recombination from either a free body or an aligned pair. Part B of
the demo:

```
realisation                 dE(p-xi,p,p+xi)   dN   dS   dP   dMom
recombine[free,free]              (-1,0,1)    -2    1    0      0
recombine[sea,sea]                (-1,0,1)     2   -1    0      0
recombine[free,sea]               (-1,0,1)     0    0    0      0
recombine[sea,free]               (-1,0,1)     0    0    0      0
ionisation                        (-1,0,1)     2   -1    0      0
hop (excluded)                    (-1,0,1)     0    0    0     -4
```

**Theorem Y2.** *All realisations move $`E`$ identically, and $`P = S + N/2`$
is conserved by all of them. With the sea ineligible, every event satisfies
$`|\Delta S| = 1`$, and no event reproduces the count signature of a
single-body momentum displacement. With the sea eligible, the mixed
realisations have $`\Delta N = \Delta S = 0`$ and displace one negaton by
$`2\xi`$ in the counts — the signature of the excluded hop — while conserving
momentum, which the hop does not, by shifting one aligned pair a row.*

Two consequences.

- **Exclusion is what gives postulate (S) content.** If the sea could be
  recombined from, every hop would have a legal representation, and "no body
  changes momentum" would be a choice of labels rather than a law. Note that
  $`\Delta S = 0`$ is invisible if the sea is carried as a per-cell scalar
  held uniform, as the demos do.
- **`recombine[sea,sea]` is ionisation in disguise**, drawing its pair from
  the daughter rows rather than the parent's, which is what Proposition K8
  forbids.

So the sea's darkness has two independent meanings, and the project has been
using one to license the other:

| | statement | status |
|---|---|---|
| observationally dark | contributes zero to $`E`$ | arithmetic, automatic |
| dynamically inert | ineligible as a recombination donor | a rule |

The emission supplement grounds the second in "a bound pair is not two
world-particles" (§2). That is a definition. §4 supplies a reason.

---

## 4. Theorem Y3: the consumed bodies are a pair

**Theorem Y3.** *The two bodies a catalysed recombination consumes — a positon
at $`p - \xi_q`$ and a negaton at $`p + \xi_q`$, co-located at the vertex —
are a winding pair in the sense of Lemma 4, with splitting $`\Delta p =
2\xi_q`$ and midpoint exactly the parent's row. Its amplitude obeys
$`|\Psi| = 2|\sin(\mu/2)|`$, its envelope has nodes spaced $`h/2\xi_q`$ and
drifts at the parent's own velocity $`p/m`$.*

*Verification.* Part C: the amplitude identity holds to $`4\times10^{-15}`$ at
two times, the node spacing matches $`h/|\Delta p|`$ to the grid resolution,
and the envelope velocity equals $`p/m`$ exactly.

Two things follow. First, the transport worry is answered: $`\mu`$ is defined
at the vertex by Lemma 4 for any two co-located bodies, and needs no standing
partnership. Second, the phase formulation's selection rule — that the
absorbed excitation comes from a pair on the particle's own row — is the
ledger's geometry restated, since the consumed pair's midpoint *is* the parent
row.

This licenses the definition the arithmetic account could not give:

> **Darkness.** A pair is dark when $`\mu \equiv 0`$, hence when its broadcast
> amplitude vanishes. **Recombination consumes amplitude, and a dark pair has
> none to give.** Ionisation is not the reverse asymmetry it looks: it *gives*
> a dark pair amplitude, which is why the sea can be ionised but not
> recombined from.

The analogy is two speakers driven in antiphase. You can unbalance them; you
cannot subtract their sound from the room, because there is none.

---

## 5. Theorem Y4: the kink and its phase ramp

If the created pair's members carry their clocks continuously through the
event, is the pair dark?

**Theorem Y4.** *Let a body kink from $`p \mp \xi`$ to $`p`$ at the event
$`(x_k, t_k)`$ with its clock continuous there. Its extended phase changes by
a ramp*

```math
\Phi' - \Phi = \frac{\pm\,\xi\,(x - x_k) - (E' - E)(t - t_k)}{\hbar},
```

*that is, by multiplication of the amplitude by $`e^{\pm i \xi (x-x_k)/\hbar}`$
up to the energy term — zero at the kink point and linear away from it. After
the event the pair has $`\Delta p = 0`$, so by Proposition 3 its misalignment
never winds again and is frozen at $`\mu(x_k, t_k)`$. The pair is dark if and
only if the kink happens at a node of the consumed pair's envelope.*

**Corollary Y4.1.** *A phase reset at the parent's position and a node-located
kink are the same freedom. Pivoting at a node is equivalent to pivoting at
$`x^*`$ together with a reset of $`\mu(x^*)/2`$ per body, applied oppositely
to the two.*

*Verification.* Part D: the ramp has slope $`\xi/\hbar`$ exactly and vanishes
at the pivot; pivoting at the parent's position leaves $`|\Psi| = 1.84`$,
uniform and unchanged at $`t = 3`$, while pivoting at the node leaves
$`|\Psi| < 3.4\times10^{-15}`$; and the two descriptions differ by
$`-1.1636`$ rad per body, twice which is $`2.3271 = \mu(x^*)`$.

So alignment at creation is not a separate phase postulate. It is one choice:
where the kink happens. **But a body can only kink where it is.** If a
world-particle carries a position below cell resolution, then $`x_k`$ is
forced to the bodies' own position and the node condition is a measure-zero
coincidence, which returns the creation rule to a gate or a postulate. If
position below cell resolution is not an attribute, the node-located kink
makes darkness forced. That is a question about the ontology of the lattice,
not about phase; it is open item Y-SP1.

---

## 6. Theorem Y5: two clocks force piecewise worldlines

Suppose an aligned pair is two co-moving clocks rather than one object, so
that recombination and ionisation neither create nor destroy world-particles.

**Theorem Y5.** *If world-particles are conserved and no body ever changes
momentum at an event, then the count of each species on each momentum row is
event-invariant, so $`E`$ evolves by classical transport alone and the
residual channel can do nothing. Hence conservation of world-particles forces
kinks.*

*Proof.* $`E = u^+ - u^-`$ is exactly that count. An event that neither
creates, destroys, nor moves a body in $`p`$ leaves every such count fixed.
$`\square`$

The phase-alignment model escapes this only because its observable is a
*label* count — which body wears the "excess" label — not a species count.
That escape is unavailable in the compensated model.

The kinks are well-behaved: simultaneous, co-located, opposite-species, equal
and opposite $`\pm\xi_q`$, so momentum is conserved event by event.
Proposition K8 is untouched, and §8.1's objection to a one-body Markov jump
(the kernel is signed) does not apply to a two-body kink at rate $`|K|`$. Each
aligned pair holds exactly one positon, so identity carries through both
processes unambiguously — cleaner than the swap, where two same-species bodies
make identity a matter of labelling.

**Vocabulary.** "Bound" is wrong: nothing holds the pair together (Lemma 1).
The relation is alignment, and what recombination creates and ionisation
destroys is that relation, not an object. Also, "recombination" already names
the bilinear $`\kappa`$ sink of §1 of the ledger note. The two differ:

| | catalysed recombination | contact recombination |
|---|---|---|
| trigger | a parent's event, rate $`\|K_q\|`$ | two free bodies in one cell |
| consumed bodies | at $`p \pm \xi_q`$ | already at the same $`p`$ |
| kinks required | $`\mp\xi_q`$ | none, $`\Delta p = 0`$ already |
| effect on $`f`$ | defines it | drives it below $`1/2`$ (Theorem N3) |

---

## 7. Theorem Y6: who crosses the barrier

Under piecewise worldlines an individual body can be carried across the
Eckart summit. Part E tags every free positon of the incident packet and
follows the tags through both processes by proportional allocation, streaming
tagged pairs classically, as postulate (S) requires of the sea.

**Theorem Y6.** *Under the piecewise reading the individual momentum walk has
zero drift, $`\sum_q \xi_q K_q = 0`$, and diffusion*

```math
D_p(x) = \tfrac{1}{2} \sum_q \xi_q^2 \, |K_q(x)| ,
```

*an unsigned moment of the kernel. $`E`$ responds only to the signed moments,
so the individual crossing rate is not a functional of the observable
evolution. For a symmetric barrier $`\Gamma(0) = 0`$ identically, so no kink
ever occurs at the summit and an individual crosses only by flank activation
followed by classical passage over the top.*

Measured (fine lattice, $`dp = 0.125`$, minimum-uncertainty packet,
$`\sigma_r = 2`$, $`\sigma_p = 0.25`$):

| $`p_0`$ | $`E_0/V_0`$ | $`T_{cl}`$ | $`T`$ exact | $`T_E`$ mesh | $`T_E`$ ledger | $`T_{tag}`$ | tag $`H>V_0`$ | $`\langle H\rangle/E_0`$ | kinks | $`f`$ |
|---|---|---|---|---|---|---|---|---|---|---|
| 1.0 | 0.50 | 0.069 | 0.191 | 0.209 | 0.213 | **0.327** | 0.204 | 2.76 | 2.40 | 0.523 |
| 1.2 | 0.72 | 0.246 | 0.369 | 0.400 | 0.404 | **0.409** | 0.268 | 1.94 | 2.32 | 0.526 |
| 1.7 | 1.44 | 0.896 | 0.845 | 0.854 | 0.855 | **0.780** | 0.723 | 1.41 | 1.83 | 0.544 |

Reading the table.

- **Individual crossing is not tunnelling.** $`T_{tag}`$ bears no fixed
  relation to $`T`$: it exceeds $`T`$ below the barrier and falls below both
  $`T`$ and $`T_{cl}`$ above it, because kicks also turn back bodies that
  would have passed. It is diffusive traffic across the separatrix in both
  directions, of which $`E`$ records only the net signed result.
- **The size matches the diffusion.** The tagged mean energy rises by
  $`\Delta H \approx 0.9`$ at $`p_0 = 1`$, against the free estimate
  $`D_p \tau/\mu \approx 1`$ for a flank passage of $`\tau \approx 2`$. Each
  body takes about two kinks per passage, and
  $`\sqrt{2 D_p \tau} \approx 1.4 \approx p_b`$: at these parameters one
  passage randomises a positon's momentum on the barrier's own scale.
- **It is not the sea's doing.** Narrowing the sea's momentum profile to a
  Gaussian of width 2 left $`T_{tag}`$ at 0.326 and $`T_E`$ unchanged at
  0.2126. An earlier conjecture that the flat sea acts as an
  infinite-temperature reservoir driving the heating is thereby refuted.
- **It depends on the allocation rule**, which the demo takes to be a uniform
  choice among eligible bodies. A phase-gated rule (§5) would change
  $`T_{tag}`$ and leave $`E`$ alone — open item Y-SP7.

The picture is Kramers' escape without friction (Kramers 1940): a body is
shoved at random on the approach slope, never passes through the hill, and
sometimes is shoved hard enough to go over.

**Lattice dependence.** Halving $`dp`$ doubled $`\Gamma`$ (3.13 to 6.58 at
$`r = 0.62`$), so event *counts* grow like $`1/dp`$, set by the reach horizon.
$`D_p`$ at the flank peak barely moved (0.586 to 0.505), but in the tail it
grew (0.114 to 0.255 at $`r = -4`$) with the reach. Whether $`T_{tag}`$ has a
continuum limit is open item Y-SP4.

---

## 8. A defect in the mean-field demo

`Ledger.channels` clamps small negative populations left by spectral
transport. The clamp breaks $`E`$: over 250 steps of a scattering run,
$`\sum E`$ falls from 1 to 0.8402 and the field departs from the exact mesh by
22.0 per cent, against 3.5 per cent with the clamp removed, where the
negatives reach only $`-2.9 \times 10^{-3}`$. It also biases $`f`$, which
moves from roughly 0.64 to 0.52 when the clamp goes. The comment in the code
is right that clipping the *caps* leaves the ledger arithmetic untouched; the
clamp on the fields themselves does not. Short runs are little affected. Open
item Y-SP6.

---

## 9. Three collisions of notation

- **Two $`\mu`$.** Ladder versus pair, as §1 of
  [`../supplement/what_the_reach_is.md`](../supplement/what_the_reach_is.md)
  records. §2 above adds that the ladder $`\mu`$ is what the compensated
  kernel has integrated out.
- **Two (S).** Streaming, in §1 of
  [`compensated_ontology.md`](compensated_ontology.md), and the sea carrier
  lock, in §4 of
  [`relational_pairing_and_carrier_lock.md`](relational_pairing_and_carrier_lock.md).
  The carrier lock is withdrawn by the phase-alignment algorithm spec;
  streaming stands. This note uses (S) for streaming throughout.
- **Two recombinations.** §6 above.

---

## 10. Open items

- **Y-SP1.** Does a world-particle carry a position below cell resolution? If
  it does, the node-located kink of §5 is unavailable and dark creation needs
  a gate or a postulate. Connects to Z-LS1 and to the per-cell sea counter.
- **Y-SP2.** Decide the creation rule: forced (node-located kink) or free
  (alignment imposed, as in the phase formulation).
- **Y-SP3.** The $`\pm\xi^2/m`$ of unsigned kinetic energy per event. It is
  zero in the mean at $`f = 1/2`$. Is that the same $`1/2`$ as Theorem N3's,
  or a coincidence? A related coincidence, not to be relied on: the P5 weight
  $`\cos^2(\mu/2)`$ averages to $`1/2`$ over uniform $`\mu`$.
- **Y-SP4.** Does $`T_{tag}`$ converge as $`dp \to 0`$ and as the reach grows?
  $`\Gamma`$ does not.
- **Y-SP5.** Repeat §7 on the exact integer ledger of
  [`stochastic_ledger.md`](stochastic_ledger.md) rather than on the mean
  field, where the allocation is a genuine random choice.
- **Y-SP6.** The clamp of §8.
- **Y-SP7.** Verify that a phase-gated allocation changes $`T_{tag}`$ and
  leaves $`E`$ exactly, as Corollary Y1.1 predicts.
- **Y-SP8.** Vocabulary: "aligned" for "bound", and a rename for the carrier
  lock.
- **J-SP3, J-SP4** are answered in part. The two vertices are not the same
  event (§1), and the phase account of darkness does transplant (§4), at the
  price of §6. **R-SP8** — whether the inventories match — is untouched here.

---

## 11. Sources

- de Broglie, L. *Recherches sur la théorie des quanta*, Ann. de Physique
  **3** (1925) 22–128. The phase winding along a worldline that P0–P5 promote
  to a particle attribute.
- Kramers, H. A. *Brownian motion in a field of force and the diffusion model
  of chemical reactions*, Physica **7** (1940) 284–304. The escape picture §7
  reproduces without friction.
- Madelung, E. *Quantentheorie in hydrodynamischer Form*, Z. Phys. **40**
  (1927) 322–326.
- Takabayasi, T. *The formulation of quantum mechanics in terms of ensemble in
  phase space*, Prog. Theor. Phys. **11** (1954) 341–373. The nearest
  historical precedent for attaching a phase to a phase-space ensemble; the
  project's reading is in
  [`../supplement/takabayasi_1954_stochastic_picture.md`](../supplement/takabayasi_1954_stochastic_picture.md).
