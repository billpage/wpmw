# Permanent pairing and the density-matrix ontology

**Status.** Analysis note, step 7 of the ladder. Companion demo:
`src/demo_pairing_resource_arithmetic.py`.

## 0. What this note inherits and retracts

**Inherits.** The phase-alignment layer (P0–P5, Lemmas 0–5, Theorems 1–4)
unchanged. From
[relational_pairing_and_carrier_lock.md](relational_pairing_and_carrier_lock.md):
Proposition R1 (a partner index carries no one-particle relational state) —
correct at its level and kept; the Theorem R4 order-parameter machinery as an
implementation device; and the §8 arithmetic of the consumable defect, which
this note confirms independently.

**Retracts.** Two inferences of the relational note, not its calculations:

1. the step from Proposition R1 to *"pairing adds no information"* — R1 is a
   statement about the one-particle marginal, and pairing carries exactly the
   information that distinguishes the density matrix from its marginal (§3);
2. Corollary R5.1, *"no replenishment"* — its premise (Lemma 3: the pump
   writes phases and only phases) is true of leg populations and false of
   pair configurations: a phase sideband on a pair *is* a split-pair
   amplitude, so the pump is precisely the source R5.1 declared absent (§5.2).

If the ontology of this note is adopted, postulate (S) (the sea carrier lock)
becomes unnecessary, and the open item on its dynamical preservation
dissolves (§7). The algorithm specification §2.2 is *not* revised here; that
is a follow-up decision.

## 1. The ontology

Each sea particle carries the world index of its partner, fixed for all time.
Excess particles carry no partner. The interpretation:

- a **pair is one Monte-Carlo sample of a density-matrix element**
  $\rho(X, X')$: the positon is the ket leg, the negaton the bra leg;
- the bra leg's complex conjugation is the negaton's opposite winding sign in
  P1 — the species distinction is the ket/bra distinction;
- the pair's internal phase difference (each leg's phase transported by P1/P2
  to a common event) is $\mu = \arg\rho$ of the sampled element;
- the overall phase of a pair is gauge: only intra-pair differences are
  physical;
- an **excess particle is a diagonal sample**: population without coherence.

The partner may sit anywhere: $\rho(X, X')$ has support at arbitrary
separation, and the excitation is a property of the pair, not of either leg
alone.

## 2. Ring corroboration

The exact ring Wigner function of a periodised packet contains an order-one
interference image at $x_0 + L/2$ with row-alternating sign (fringe factor
$e^{ipL/\hbar}$, period $2\thinspace dp$). In this ontology that image is the
population of pairs whose legs sit one period apart — off-diagonal support
made visible by the ring topology.

## 3. The coboundary objection, resolved by level

Re-pairing all partners at random leaves every one-particle observable
unchanged: at that level Proposition R1 stands and pairing adds nothing. But
re-pairing scrambles the off-diagonals of $\rho$. A sample of $\rho(X, X')$
is irreducibly a bound two-point object; *which two points belong to one
sample* is the entire content of the two-point function beyond its
marginals. With per-pair gauge, $\mu$ is not the coboundary of any global
one-index field, so the objection's premise fails at the two-point level.
Fixed pairing adds nothing to $f(x, p)$ and adds precisely everything that
distinguishes $\rho$ from $f$.

## 4. One-leg hops derive the stencil and the mediated counting

**Lemma (no-go for lone hoppers).** Consider any per-particle jump process
in which a particle hops $p \to p \pm 2q\thinspace dp$ at rates
$a(x), b(x) \ge 0$ depending only on position (and the particle's own
state), with net drift matched to the QLE, $a - b = \Gamma_q$. Its generator
decomposes into an antisymmetric part (the drift) and a symmetric part

```math
\tfrac{a+b}{2}\left[W(p+2q\thinspace dp) + W(p-2q\thinspace dp) - 2W(p)\right],
```

whose coefficient obeys $(a+b)/2 \ge \lvert\Gamma_q\rvert/2$: irreducible
spurious momentum diffusion. The QLE generator for a cosine mode is purely
antisymmetric (per-column eigenvalues $`2i\Gamma_q\sin(2\pi k_p q/N)`$), so no
such process reproduces it. $\square$

Permanent pairing evades the lemma in the only way possible. In the momentum
basis the commutator with $V = V_q\cos(Kx + \phi_q)$ acts on $\rho(P, P')$
through four terms, each hopping **one leg** by
$\pm\hbar K = \pm 2q\thinspace dp$:

```math
i\hbar\thinspace\partial_t\rho(P,P') = \tfrac{V_q}{2}\left[
e^{i\phi_q}\rho(P{-}\hbar K, P') + e^{-i\phi_q}\rho(P{+}\hbar K, P')
- e^{i\phi_q}\rho(P, P'{+}\hbar K) - e^{-i\phi_q}\rho(P, P'{-}\hbar K)\right].
```

A one-leg hop moves the pair midpoint $p = (P + P')/2$ by $q\thinspace dp$ —
a half quantum — and shifts the separation variable, contributing the factor
$e^{\pm iKx}$ after the Wigner transform. Collecting the four terms yields,
in two lines,

```math
\partial_t W = \Gamma_q(x)\left[W(p + q\thinspace dp) - W(p - q\thinspace dp)\right],
\qquad \Gamma_q(x) = -\tfrac{V_q}{\hbar}\sin(Kx + \phi_q),
```

the crystal-lattice stencil with the repository sign convention.
Consequences:

- the stencil's half-quantum offsets are literally one-leg hops of bound
  pairs; the rate is per pair, and pairs live at their midpoints in Wigner
  coordinates, so the mediated rule's *rate proportional to occupancy at the
  transition midpoint* is **derived**, not postulated;
- the elementary mover is a leg of a bound object, the Wigner-cell object
  advances $q\thinspace dp$ per event, and the generator is natively the
  antisymmetric central difference — zero spurious diffusion;
- **conjecture** (to check against
  [four_rule_microdynamics_equivalence.md](four_rule_microdynamics_equivalence.md)):
  the four commutator channels — ket or bra leg, up or down — are the four
  rules. Each one-leg hop changes the midpoint by $\pm q$ and the splitting
  by $\pm 2q$, the same two-by-two structure as
  Focus/Defocus/Right-Hop/Left-Hop.

## 5. Resource arithmetic

Numbers from `demo_pairing_resource_arithmetic.py` on the canonical
cosine-well trajectory ($V_p = 1.5$, $L = 8$, SQUEEZE $= 2$,
$T = 4\thinspace T_{\mathrm{period}}$, exact split-operator Schrodinger
arm). State-mass units: the state carries total mass 1; the aligned sea
carries $(2/h)\thinspace L\thinspace dp = 1$ per momentum row and
$(2/h)\thinspace dx\thinspace dp \approx 0.0104$ per cell.

### 5.1 Consumable accounting fails, as §8 held

If every exchange permanently spends a pair, demand is cumulative moved
mass: $E_{\mathrm{state}}(T) \approx 80$ state-mass units on accepted
events, and the busiest cell exhausts its aligned stock at
$t^{*} \approx T/76$. With the attempt-to-accept ratio of the bare $g_0$
traffic (one to two orders), this reproduces, $\nu$-independently, the
relational note's shortfall of $7.7 \times 10^{2}$ (its §8, with
$N_{\mathrm{exc}} = 5 \times 10^{6}$, $M = 128$, $B \approx 2 \times 10^{3}$).
Any consumable model fails for long enough runs: demand grows linearly in
$T$ against fixed stock. This half of the earlier objection was sound.

### 5.2 The pump is the source R5.1 declared absent

Proposition R5 stands: a struck pair exits with $\Delta p = 0$. What falls
is the claim that it stays spent. In the density-matrix accounting the
potential acts on every pair by the one-leg hops of §4: acting on an
aligned pair ($P = P'$) it creates the split configuration with amplitude
$-i\thinspace V_q\tau_p/2\hbar$ per kick — which is exactly the sideband
amplitude $\mu_1/2$ that Lemma 3 books as "phases." The two descriptions
are the same physics at two levels: *phase sideband on the pair* (leg
population language) equals *split-pair amplitude* (pair configuration
language). The mediating resource is continuously regenerated by the pump
at the rate the vertex consumes it; the split population is not
monotonically non-increasing, and the drain of §8 has a source.

### 5.3 Storage capacity: the Wigner bound, tight

An exchange converts an aligned pair into a split pair — legs at different
momenta, internal phase running at the beat frequency — which stores one
sample of an off-diagonal element and is retrievable. Demand is the
instantaneous coherence content, not the event count.

Each pair contributes signed local Wigner content bounded by its own
weight, so a sea of pair density $2/h$ hosts any local content with
$\lvert W\rvert \le 2/h$ at polarisation at most 100%. That inequality is
Wigner's bound, valid for every state at every time — and it is the same
inequality that guarantees $W' = W + 2/h \ge 0$ (as stated for $d$
dimensions in [multi_body_extension.md](../algorithm/multi_body_extension.md)).
The positivity of the shifted representation and the feasibility of
permanent pairing are one statement.

Trajectory numbers: the load factor $\lvert W\rvert\thinspace(h/2)$ starts
at $0.960$ — below 1 only through the real ring tail-overlap
renormalisation of the wide packet; a narrow control packet gives
$0.999978$ — and never exceeds it; $\min W\thinspace\pi\hbar = -0.943$,
within 6% of the negative bound mid-run. The sea is exactly big enough,
with equality at Gaussian peaks and deep interference fringes: the bound
is tight, not generous.

### 5.4 The global ledger is comfortable

The off-diagonal mass $C(t) = \sum_{P \neq P'} \lvert\rho(P, P')\rvert$
peaks at $10.3$ against aligned stock of one unit per momentum row over
the roughly eleven occupied rows. Half the coherence mass sits at pair
splittings reachable in at most 4 one-leg hops; 90% within 32.

### 5.5 Flow feasibility is algebraic

The mediated rate reads $W' \ge 0$ at the transition midpoint, and

```math
\left\lvert W'(p + q\thinspace dp) - W'(p - q\thinspace dp)\right\rvert \le
W'(p + q\thinspace dp) + W'(p - q\thinspace dp)
```

by the triangle inequality — that is, by $W' \ge 0$, Wigner's bound once
more. The stencil is itself the polarisation-flow ledger. Numerical check:
$\min W'\thinspace(h/2) = 0.057$ over the trajectory, non-negative.

![Resource arithmetic](https://raw.githubusercontent.com/billpage/wpmw/output/figures/pairing_resource_arithmetic.png)

## 6. The locality requirement (lemma to be proven)

At a vertex the struck leg must present $\mu$ without querying its remote
partner instantaneously. Candidate mechanism: each leg carries the pair's
internal clock locally, set at the last co-interaction of the pair and
evolved by the leg's own P1/P2 transport. This must be proven as a lemma:
it is exactly where the nonlocality of the off-diagonals of $\rho$ must be
shown to require no signalling.

## 7. Consequences and open items

1. **Postulate (S) becomes unnecessary.** The carrier lock existed because
   the all-pairs relational average vanished as $1/B$ without a shared
   reference; with permanent pairing the vertex reads $\mu$ of the pair
   involved. The open item on dynamical preservation of (S) dissolves.
   Whether the algorithm specification §2.2 reverts to indexed pairing is
   a separate decision; Theorem R4's $Z_r$ factorisation remains available
   as an implementation of the sea average either way.
2. **Species bookkeeping** acquires a candidate resolution: species = ket
   or bra side of $\rho$. To be checked against page 4 of the Cyganski
   slide deck.
3. **The load-bearing theorem for the direct proof** is now single and
   sharp: *split pairs mediate subsequent exchanges with the same vertex
   constant as pump-excited pairs.* The vertex Hamiltonian
   $h = g_0 + g_1\mu_1 e^{i\mu}$ is the same object in both cases, which
   makes the theorem plausible; it is not yet proven.
4. **The decisive numerical test**: a pair-ensemble Monte Carlo with
   permanent partners and only local vertex rules — one-leg hops for
   pairs, Theorem-4 swaps for excess particles — with the mediated
   counting never coded. If QLE tracking emerges, §4 is its theorem and
   the simulation its witness.
