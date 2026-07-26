# Phase-Alignment Microdynamics: The Contact Interaction Without Beats or Resonance

## 0. Status and provenance

This note is a **change of variables** applied to
[`phase_resonance_microdynamics.md`](phase_resonance_microdynamics.md). It
introduces no new postulates, discards none of that note's results, and
predicts nothing different. What it does is replace three constructs —
the *beat*, the *grating*, and the *resonance condition* — with a single
scalar carried by a pair, and show that every theorem stated in terms of
the discarded constructs survives, most of them in shorter form.

The motivation was a pair of questions (B. Page, July 2026): the word
"beat" is awkward when applied to a single pair, since a beat is
ordinarily a pattern in an extended field rather than a property of two
particles; and can the interaction be stated in world-particle terms
alone, without invoking resonance? Both are answered affirmatively. The
relevant scalar is the **misalignment** $\mu$ — the difference between
the two partners' transported clock phases, evaluated at one point — and
the answer to the second question is stronger than expected: requiring
$\mu$ to hold still during a vertex forces the vertex to be a **momentum
exchange**, from which energy conservation follows rather than being
imposed (Theorem 4).

One substantive correction is recorded here. During the conversation that
produced this note the exchange rule was first stated with the two roles
reversed — the particle arriving on the struck partner's momentum rather
than the mate's. Part C of the companion demo settles it: the particle
arrives on the *mate's* momentum and leaves on the *struck partner's*.
See §4 and the caveat in §4.3.

Derivation and implementation were developed jointly with Claude
(Anthropic), July 2026. Companion code: `src/demo_phase_alignment.py`.
All numerical claims in §8 are outputs of that demo.

## 1. Tutorial: the interaction in five steps

Informal; everything asserted here is stated and proved in §§2–6.

### 1.1 One clock per world-particle

Every world-particle carries a phase that winds along its worldline at
the Lagrangian rate (P1 of the predecessor). Ask what phase a particle
would show if its clock were carried to some other point at the
particle's own wavevector, and you get the **transported phase** — which
is exactly what P2's extended amplitude encodes. Nothing here is new; it
is P1 and P2 read as a transport rule rather than as an emission.

### 1.2 A pair has exactly one number

Two partners, two clocks. Transport both to a common point and compare:
the difference $\mu$ is the *only* relational quantity a pair has. It is
invariant under a global phase shift and under re-referencing either
worldline, and the pair's amplitude depends on the two partners through
$\mu$ alone. A pair is not a particle plus a particle plus a field; it is
two particles and one angle.

### 1.3 Three states, one of which used to be called a beat

$\mu$ is a function of where you evaluate it. Three cases exhaust the
possibilities:

- **aligned** — $\mu \equiv 0$ everywhere: the pair is dark;
- **uniformly offset** — $\mu$ constant and nonzero: the predecessor's
  gray pair;
- **winding** — $\mu$ turns with position at rate $\Delta p/\hbar$: the
  predecessor's *beating* pair.

The third case is where "beat" was doing its work, and the reformulation
shows why the word chafed. Nothing propagates and nothing is quantised.
The two clocks simply run at different rates, so their gap turns as you
walk along the ring. There is no third object to count, which is why the
predecessor had to append the disclaimer that beat number is not
conserved: a quantity that was never a thing does not need a conservation
law.

### 1.4 The pump sets $\mu$ by place, not by pair

The predecessor's Lemma 2 said the pump writes the same pattern phase on
every pair regardless of location. In the present variables that reads:
*every pair sitting at the same place carries the same $\mu$.* This is
the entire content of the "coherent grating." An arriving particle reads
one number off whichever pair it meets, and never needs to know which
pair that was or where the absent partner is. Locality of the vertex is
a consequence of mutual alignment across the sea, not of anything local
about a single broadcast.

### 1.5 The vertex is a momentum exchange

An excess particle sharing a cell with a sea partner may exchange
momentum with it. The exchange is available anyway — the bare rate
$w_0$ — and $\mu$ only tilts the odds. Which exchanges survive averaging
is fixed by one condition: $\mu$ must hold still while the exchange is
happening. Impose it, and the vertex is forced to be a *swap*: the
particle takes the struck partner's momentum and the partner takes the
mate's. The pair exits aligned because the struck partner now matches
its mate. Momentum and energy are conserved automatically, because a
swap is a permutation of momenta among worldlines.

![The contact interaction in world-particle terms](https://raw.githubusercontent.com/billpage/wpmw/output/figures/phase_alignment_contact.png)

*(a) transporting two clocks to a common point defines $\mu$. (b) the
three states of a pair. (c) the pump makes $\mu$ a function of place.
(d) the vertex as a momentum swap. (e) the stationarity condition.
(f) the ensemble limit.*

## 2. The misalignment

**Definition (transported phase).** For a world-particle $j$ on a leg
with reference data $(x_j, p_j, \theta_j)$ at the leg's epoch,

```math
\Phi_j(x,t) \;=\; \theta_j \;+\; \frac{p_j\,(x - x_j) - E_j\,t}{\hbar},
\qquad E_j = \frac{p_j^2}{2m},
```

on a free leg; with a potential, $E_j$ is replaced by the leg's
Hamiltonian value in the usual way. By Lemma 0 of the predecessor this
is exactly the phase of P2's extended amplitude, so $\Phi$ is a rewriting
of P2 and not an addition to it.

**Definition (misalignment).** For a pair with partners $a$ (positon) and
$b$ (negaton), the misalignment at an event $(x^{*}, t^{*})$ is

```math
\mu \;=\; \Phi_a(x^{*}, t^{*}) \;-\; \Phi_b(x^{*}, t^{*}) \pmod{2\pi}.
```

**Lemma 4 (completeness).** $\mu$ is the pair's entire gauge-invariant
relational content:

1. *Global gauge.* $\mu$ is unchanged by $\theta_a, \theta_b \mapsto
   \theta_a + \alpha, \theta_b + \alpha$.
2. *Re-referencing.* Advancing either partner's reference data along its
   own worldline leaves $\mu$ unchanged, up to the corresponding shift of
   the time origin. Transport is exact, so "which epoch the clock was
   read at" is not physical data.
3. *Sufficiency.* The pair amplitude satisfies
   $\lvert\Psi\rvert = 2\lvert\sin(\mu/2)\rvert$, so all observable
   consequences of the two partners' phases enter through $\mu$.

Consequently the pair's state is $(p_a, p_b, \mu)$ and nothing else;
darkness is $\mu \equiv 0$, and the predecessor's gauge-matching
condition for a separated dark pair is the same statement.

## 3. The three states, and the disposal of the beat

**Proposition 3 (winding rates).** For a pair with
$\Delta p = p_a - p_b$ and mean velocity
$\bar v_{\mathrm{pair}} = (p_a + p_b)/2m$,

```math
\frac{\partial \mu}{\partial x} \;=\; \frac{\Delta p}{\hbar},
\qquad
\frac{d\mu}{dt}\bigg|_{\text{path of velocity } v}
\;=\; \frac{\Delta p}{\hbar}\,\big(v - \bar v_{\mathrm{pair}}\big).
```

*Proof.* Immediate from $\Phi_j = (p_j x - E_j t)/\hbar + \mathrm{const}$
and $E_a - E_b = (p_a + p_b)\Delta p/2m$. $\square$

Hence the classification of §1.3, with the dictionary

| predecessor | here | condition |
| --- | --- | --- |
| dark | aligned | $\mu \equiv 0$ |
| gray | uniformly offset | $\Delta p = 0$, $\mu \neq 0$ constant |
| beating | winding | $\Delta p \neq 0$ |

**Corollary (the beat is dispensable).** The predecessor's beat
equivalence class $(\Delta p, \bar n, \chi)$ is the triple
$(\Delta p,\thinspace \bar v_{\mathrm{pair}},\thinspace \mu\text{ at a reference event})$,
i.e. it is a repackaging of Proposition 3 plus one constant of
integration. Proposition 1 of the predecessor — full-contrast grating at
wavelength $L/q$, envelope drifting at the mean velocity — is the
statement that $\mu$ advances by $2\pi$ over a distance $h/\Delta p$ and
holds still on paths at $\bar v_{\mathrm{pair}}$. No object need be
introduced to say this, and with no object introduced there is no
temptation to count one.

Two consequences of retiring the word are worth recording, since they
were live confusions in the predecessor. First, the note's own warning
that "nothing about one beat is small; smallness is only relative to the
$B$ dark pairs" conflated two different objects: a *fully split* pair,
for which $\mu$ winds through the full range, and a *pumped* pair, whose
$\mu$ is displaced from zero by $O(V_p\tau_p/\hbar)$. Only the latter is
what the pump makes and what P5 couples to. In the $\mu$ language they
are simply different values of one variable, and the ambiguity cannot
arise. Second, the phrase "the beat ceases to exist" at a K3 vertex
becomes "$\mu$ returns to zero," which requires no reconciliation with
worldline continuity.

## 4. The exchange theorem

### 4.1 Statement

Consider a vertex at which an excess particle $c$, arriving on momentum
$p_{\mathrm{in}}$ and leaving on $p_{\mathrm{out}}$, exchanges momentum
with the struck partner $a$ of a pair $(a, b)$, the struck partner
leaving on its mate's momentum $p_b$. During the encounter window the
particle's amplitude carries both $p_{\mathrm{in}}$ and
$p_{\mathrm{out}}$, so the event advances at the midpoint velocity
$v_{\mathrm{exch}} = (p_{\mathrm{in}} + p_{\mathrm{out}})/2m$
(Proposition 1 of the predecessor, applied to the excess particle).

**Theorem 4 (exchange).** Momentum conservation together with
stationarity of $\mu$ along the event path,

```math
p_{\mathrm{out}} - p_{\mathrm{in}} \;=\; p_a - p_b,
\qquad
\dot\mu = 0 \iff v_{\mathrm{exch}} = \bar v_{\mathrm{pair}}
\iff p_{\mathrm{in}} + p_{\mathrm{out}} = p_a + p_b,
```

has the unique solution

```math
p_{\mathrm{in}} = p_b, \qquad p_{\mathrm{out}} = p_a .
```

The vertex is therefore a **swap**: the excess particle and the struck
partner interchange momenta, while the mate is untouched.

*Proof.* Add and subtract the two displayed conditions. $\square$

**Corollary 4.1 (energy conservation is automatic).** A swap permutes
the multiset of momenta $\lbrace p_{\mathrm{in}}, p_a, p_b\rbrace \mapsto
\lbrace p_{\mathrm{out}}, p_b, p_b\rbrace$ with
$p_{\mathrm{out}} = p_a$, $p_{\mathrm{in}} = p_b$, so both
$\sum p$ and $\sum p^2/2m$ are preserved identically. Energy
conservation is not an additional requirement and was never an
independent input.

**Corollary 4.2 (the pair exits aligned).** After the swap both partners
are at $p_b$, so $\Delta p = 0$ and the pair is aligned or uniformly
offset; with the phase continuity of §6 it is aligned. The predecessor's
"the pair exits dark" is recovered.

**Corollary 4.3 (transfer size).** $\lvert p_{\mathrm{out}} -
p_{\mathrm{in}}\rvert = \lvert \Delta p \rvert$: the particle can only
take what the pair's own splitting offers. For a $q$-mode pumped pair
this is exactly $2q \cdot dp$, recovering item 4 of the predecessor's
Theorem 3 without invoking Bragg matching as a separate principle.

### 4.2 What this replaces

Theorem 4 subsumes, in two lines of algebra:

- **Proposition 2** of the predecessor (the co-moving-crest condition);
- the **row-resonance selection rule**, item 5 of Theorem 3;
- the separate imposition of **energy conservation** at vertices;
- the appeal to **Bragg phase matching** to fix the transfer.

Off-stationary assignments do not need to be forbidden; they dephase.
The time-averaged bias over an encounter window $\tau_e$ carries the
envelope $\lvert\mathrm{sinc}(\dot\mu\thinspace\tau_e/2)\rvert$, verified
in §8 Part C.

### 4.3 A caveat on translating to row indices

The predecessor states the selection rule as
$\bar n_{\mathrm{pair}} = n_{hi}$. In the present variables the
condition is $p_{\mathrm{in}} + p_{\mathrm{out}} = p_a + p_b$, i.e. the
pair's midpoint row equals the *transition's* midpoint row, which for a
$q$-mode transfer sits $q$ half-quanta below the particle's incoming
row. Whether these agree depends on the indexing convention for
$n_{hi}$, which is exactly the convention flagged as unresolved in the
four-rule session (page 4 of the Cyganski slide deck). **This note does
not assert that the predecessor's row statement is wrong**; it asserts
that the momentum-space statement above is what the algebra gives, and
that the translation to row indices must be checked against the
convention before either form is relied on in code. Marked as open item
1 below.

## 5. Locality of the misalignment

**Lemma 5 (locality; Lemma 2 restated).** After the pump, the
misalignment of any pumped pair evaluated at an event $(x, t)$ is

```math
\mu^{s}(x,t) \;=\; sKx + \frac{\pi}{2} - sK\,\bar v^{s}\,t,
\qquad s = \pm 1,
```

independent of the pair's location. Every pair at the same place carries
the same $\mu$.

This is the microscopic content of the "coherent polarization grating,"
and it is what makes $\mu$ a *local* datum despite the transported phase
being defined over all space. The vertex reads a number determined by
where the vertex is, not by the absent partner's whereabouts. Note the
logical direction carefully: locality of the vertex is *derived from
mutual alignment across the sea*, and therefore inherits its conditions.
If the sea decoheres, $\mu$ ceases to be a function of place alone and
the vertex acquires explicit dependence on the mate's state — the same
condition that governs linearity in $V_p$ (§6) and the classical
readability of the vertex.

## 6. The vertex rule

**Affine form.** With V1–V4 of the predecessor's §13.1 (worldline
definiteness, relational gauge invariance, linear response, mode
additivity), the exchange probability at a co-location is

```math
P(\mathrm{exchange}) \;=\; w_0 \;+\; \kappa\thinspace\cos\mu \;+\; O(\mu_1^2),
```

with $w_0$ pure gauge (the $G$-freedom), $\kappa$ absorbed into the
calibration, and no phase offset. Here $\mu_1$ denotes the pumped
displacement of $\mu$ from alignment; the predecessor's contrast $C$ is
$\lvert\mu_1\rvert$ to first order, which is the cleanest way to see that
$C$ was never a property of a "beat" but a measure of how far a pair had
been tipped out of alignment.

**Contact form.** The unitary reduction of the predecessor's §13.2 is
unchanged, with $\delta$ read as $\mu$:

```math
h \;=\; g_0 \;+\; g_1\thinspace\mu_1\thinspace e^{i\mu},
\qquad
P(\mathrm{exchange}) \;=\; \sin^2\big(\lvert h\rvert\thinspace\tau_e\big),
```

giving the offset $\delta_0 = 0$ from hermiticity, detailed balance from
unitarity, and saturation from the $\sin^2$. The bare amplitude $g_0$ is
available at every encounter; alignment does not gate the event, it
tilts it. Nothing is absorbed and nothing is emitted — a phrasing the
predecessor could not adopt because it had introduced quanta to absorb.

**Phase bookkeeping at the vertex.** The excess particle's carried phase
is continuous through the swap. The struck partner's exit phase is fixed
by the requirement that the pair leave aligned, which by Lemma 4 is one
equation for one unknown and hence always solvable from vertex-local
data. The mate is untouched in worldline, momentum and phase.

## 7. Translation table

| phase-resonance note | this note |
| --- | --- |
| beat, grating quantum | winding of $\mu$ |
| beat contrast $C$ | pumped displacement $\lvert\mu_1\rvert$ |
| pattern phase $\delta$ at the vertex | $\mu$ at the vertex |
| dark / gray / beating pair | aligned / uniformly offset / winding |
| coherent polarization grating | $\mu$ is a function of place (Lemma 5) |
| Proposition 2 (co-moving crests) | $\dot\mu = 0$ (Theorem 4) |
| row-resonance selection rule | $\dot\mu = 0$ (Theorem 4) |
| Bragg matching fixes the transfer | Corollary 4.3 |
| energy conservation at vertices | Corollary 4.1 (automatic) |
| absorption K3 / emission K4 | exchange, in the two directions |
| "the beat ceases to exist" | $\mu$ returns to zero |

## 8. Numerical verification

All claims are exercised by `src/demo_phase_alignment.py` (container run,
July 2026):

- **Part A (Lemma 4).** $\mu$ is invariant under a global phase shift to
  $1.1 \times 10^{-16}$ and under worldline re-referencing to
  $5.6 \times 10^{-16}$; the identity
  $\lvert\Psi\rvert = 2\lvert\sin(\mu/2)\rvert$ holds exactly.
- **Part B (Proposition 3).** $\partial_x\mu = \Delta p/\hbar$ and
  $\dot\mu = (\Delta p/\hbar)(v - \bar v_{\mathrm{pair}})$ at four
  velocities, worst deviation $3.6 \times 10^{-10}$ (finite-difference
  step, not model error); $\dot\mu$ vanishes at
  $v = \bar v_{\mathrm{pair}}$.
- **Part C (Theorem 4).** The two-condition linear system returns
  $p_{\mathrm{in}} = p_b$ and $p_{\mathrm{out}} = p_a$ with zero
  residual; total momentum and total energy are each conserved to
  $0$ exactly, without either having been imposed; the exchanged
  momentum equals the mode quantum $2q \cdot dp$. Time-averaged bias
  over $\tau_e = 60$ falls from $\cos\mu_0 = 0.5403$ at stationarity to
  $\le 0.075$ at every nonzero detuning tested, bounded throughout by
  the analytic $\lvert\mathrm{sinc}\rvert$ envelope.
- **Part D (Lemma 5).** Over 60,000 pairs pumped at uniformly random
  locations, the spread of $\mu$ at a fixed evaluation event is
  $8.9 \times 10^{-16}$ for one family and $4.4 \times 10^{-16}$ for the
  other.
- **Part E (quadrature, linearity).** $\cos\mu(x) \propto -V'(x)$ to
  $2.6 \times 10^{-13}$ relative; the contact-vertex cross term
  approaches the analytic linear coefficient as
  $1.0007, 1.00007$ for contrast $10^{-3}, 10^{-4}$.
- **Part F.** The stencil preserves L1 on the ring to
  $2.2 \times 10^{-16}$.

## 9. What changes and what does not

**Unchanged.** P0–P3, Lemma 0, Lemma 1, Theorem 1 (parity and the
doubled grid), Lemma 3 (the pump is invisible to populations),
Theorem 2 (the no-go), Lemma 2 (here Lemma 5), the affine
representation theorem, the contact-vertex reduction, and every
numerical result of the predecessor's §9. The rate law, the quadrature,
the direction field $\sigma$, and the factor $\gamma/2$ are untouched.

**Strengthened.** Energy conservation and the selection rule are now
consequences of one condition rather than two independent inputs, and
the vertex is identified as a momentum swap — a sharper statement than
"elastic two-body collision," since it fixes the out-state uniquely
rather than merely constraining it.

**Removed.** Three constructs with no remaining work to do: the beat as
an object, the grating as an object, and resonance as a principle.

**Not addressed.** Everything in the predecessor's §10 scope statement
still applies. This is a reformulation of the single-particle sector
with an external potential; the N-body case is untouched, and the Born
resolution of a definite vertex outcome remains exactly where §13.6 left
it. Restating the residue in the present variables: the vertex is
classical two-mode coupled-wave dynamics throughout, and what is
irreducible is the *interface* between that continuous process and the
demand (V1) that worldlines stay definite. That is a sharper location
for the residue than "quantum probability applied to a binary contact,"
but it is the same residue.

## 10. The cost of the reformulation

Honesty requires recording what is lost. The predecessor attaches to
each pair a token $(\Delta p, \bar n, \chi)$ that can be counted, and
open item 3 of that note proposes a row-resolved excitation ledger
$\rho_{\pm 2q}(x, \bar n)$ built from those tokens. In the present
variables there is no token. The ledger must be recast as a
distribution over pair momenta and misalignments,
$f(x, p_a, p_b, \mu)$, or over the reduced variables
$(x, \bar p, \Delta p, \mu)$.

This is a real edit to the live-sea code, not a rename. The argument for
paying it is that the token was never conserved and never had a
worldline, so a ledger built on it was tracking a bookkeeping device
rather than a state; the misalignment distribution tracks the state
itself. But the recast is work, and until it is done the predecessor's
formulation remains the implementable one.

## 11. Open items

1. **Row-index translation.** Reconcile Theorem 4's momentum-space
   condition with the predecessor's $\bar n_{\mathrm{pair}} = n_{hi}$
   under the indexing convention of the four-rule slide deck (page 4),
   which is still unresolved. Until then, prefer the momentum form in
   code.
2. **Misalignment ledger.** Recast open item 3 of the predecessor as a
   distribution over $(x, \bar p, \Delta p, \mu)$ and re-derive the
   steady state under continuous pumping in those variables.
3. **Multi-mode superposition.** In the $\mu$ language a multi-mode pump
   gives each pair a misalignment that is a sum of mode contributions;
   verify first-order additivity and locate the cross-mode terms
   (predecessor open item 5).
4. **Decoherence signatures.** Lemma 5, the linearity of §6, and the
   classical readability of the vertex all rest on the same condition:
   mutual alignment across the sea. Degrading it should produce three
   correlated signatures — loss of locality of $\mu$, a linear-to-
   quadratic crossover in the response, and a change from $N$ to
   $\sqrt N$ in the ensemble scaling. Testing whether they arrive
   together is a sharper probe than any one alone.
5. **Uniformly offset pairs.** The dynamics of the middle row of the
   table remains open, as in the predecessor's item 6. Note that
   Corollary 4.2 leaves a post-vertex pair with $\Delta p = 0$ but
   generally separated, so a nonuniform potential will drive it into the
   uniformly-offset state at rate
   $[V(x_a) - V(x_b)]/\hbar$.

## 12. Relation to existing formalisms

**de Broglie's phase harmony.** The winding law P1 with Lemma 0 is de
Broglie's *théorème de l'harmonie des phases*, stated in Chapter 1 §1.1
of the 1924 thesis *Recherches sur la théorie des quanta* (Ann. de Phys.
(10) **3**, 22, 1925; English translation by Kracklauer, linked in
`references/bibliography.md`). De Broglie sets an internal periodic
phenomenon of rest frequency $\nu_0 = m_0c^2/h$ against the
quantum-relation frequency of the moving body, observes that the two
disagree, and resolves it by requiring that the internal phenomenon stay
constantly in phase with an accompanying wave. That is a
carried-phase-versus-transported-phase agreement condition, which is the
structure of $\mu$ — though the analogy is not an identity: de Broglie's
harmony holds between one particle's clock and its own wave, whereas
$\mu$ compares two different particles' clocks.

Closer to the present construction is the mechanical model de Broglie
gives immediately after the theorem: a disk carrying identical
spring-suspended weights, all oscillating in phase, viewed from a moving
frame. Each weight becomes a clock showing time dilation, the motions
*dephase* relative to one another, and the locus of their centres of
mass is a sinusoidal surface travelling faster than light. He uses this
specifically to explain why such a wave transports phase but not energy.
An array of clocks whose relative dephasing constitutes the wave is
precisely the sea of §5, and the "transports phase, not energy" reading
is exactly what the model needs to avoid attributing forces to the
grating.

**Two-beam dynamical diffraction.** The vertex of §6 — two modes, one
Hermitian coupling, $\sin^2$ transfer, Pendellösung between
the states — is classical coupled-wave theory, not a quantum-specific
structure (Batterman & Cole, *Rev. Mod. Phys.* **36**, 681, 1964). The
underlying scattering process is Kapitza–Dirac diffraction (Kapitza &
Dirac, *Proc. Camb. Phil. Soc.* **29**, 297, 1933), as the predecessor
noted. Recognising the vertex as classical wave optics is what makes the
localisation of the quantum residue in §9 defensible.

**Direct interparticle action.** The transported phase has no energy, no
self-interaction, no radiation reaction and no retardation, because it
has no field degrees of freedom: it is a rule for comparing worldlines,
not an emission. The structural precedent is Fokker–Wheeler–Feynman
absorber theory (Wheeler & Feynman, *Rev. Mod. Phys.* **17**, 157,
1945), with a phase in place of a retarded potential. The darkness
condition of Lemma 4 is a parallel-transport agreement — a connection
rather than a signal.

**Stochastic mechanics.** The observation that the bare exchange traffic
$w_0$ is constitutive rather than removable — no noise, no force — is
the moral of Nelson's construction (*Phys. Rev.* **150**, 1079, 1966),
reached here from a different direction.
