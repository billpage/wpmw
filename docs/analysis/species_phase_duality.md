# Species and phase: one degree of freedom in two bases, and what that costs annihilation

**Status.** Analysis note, step 13 of the ladder. Companion demo:
`src/demo_species_phase_duality.py`. Prompted by B. Page's question (August
2026): *something seems missing in the relationship between world-particle
sign and phase — the dark state of a sea pair was described as a cancellation
of phase, and that cancellation ought to depend only on the relative position
of the pair; and the "non-linear post-processing" that annihilation requires
is probably non-local in position space.*

Both halves of that intuition turn out to be right, and both turn out to be
sharper than stated. The sea's darkness depends only on relative position
because that relative position is exactly zero, and the reason is not a
cancellation at all. The non-locality is real, is forced rather than chosen,
and its reach is the coherence length — not the lattice spacing, which is what
the annihilation-burden section had assumed.

---

## 0. What this note inherits, supplies, corrects and concedes

**Inherits.** From
[`interworld_coupling.md`](interworld_coupling.md): the reading of the
potential as a leg coupling $`U(X,Y) = V(X + Y/2) - V(X - Y/2)`$ in midpoint
and *full* separation coordinates, the factorisation of $`U`$ for one cosine
mode (Proposition I2), the Fourier-duality theorem for the channel set
(Theorem I3), and the winding law $`d\mu/dt = -U/\hbar`$ (Theorem I5). From
[`position_pair_ladder.md`](position_pair_ladder.md): the pair ontology, positon
= ket leg, negaton = bra leg. From
[`phase_alignment_microdynamics.md`](phase_alignment_microdynamics.md): the
misalignment $`\mu = \Phi_{\rm ket} - \Phi_{\rm bra}`$ as the sole relational
datum. From
[`sea_dressed_microdynamics.md`](sea_dressed_microdynamics.md): the crystal
shift as a live medium, the recombination channel R, and the broken detailed
balance of §6. From
[`../supplement/representation_cost_and_annihilation.md`](../supplement/representation_cost_and_annihilation.md):
the cost functional $`N_{\rm eff} = N/\|\mu\|_1^2`$, the annihilation burden,
and open items N1–N5.

**Supplies.** Six results, D1–D6 (§2–§6); a recap of the unraveling invariance
that makes annihilation the only remaining lever (§7); a specification for an
annihilation substep suitable for the algorithm notes (§8); and its cost
measured against its benefit (§9).

**Corrects four statements.**

1. **§7.4 of the supplement, criterion (ii).** That criterion measured the
   time for a created pair to separate out of a position cell, concluded that
   refining the position lattice makes annihilation *worse* linearly in
   $`M`$, and left the choice between it and criterion (i) undecided. Its
   premise — that the original partner must be found — is void. Theorem D5.1
   licenses anonymous annihilation exactly, for all future time and not
   merely for present observables. Criterion (i), occupancy, governs, and §9.5
   measures no threshold above occupancy 1.
2. **§8 (N1) of the supplement, naming.** N1 writes the action
   $`(n,n) \to (n-q, n+q)`$ as the *focus* action. In §4–§5 of the
   sea-dressed note the focus stencil is $`(-1,+2,-1)`$, concentrating at
   $`n`$; the outward $`(+1,-2,+1)`$ move is *defocus*. §7.4 repeats the
   same slip. The outward action is defocus.
3. **§8 (N1) of the supplement, the mean-field inference.** N1 argues from the
   boundedness of $`\|W(t)\|_1`$ that "creation and recombination balance" in
   mean field. Boundedness says cancellation *occurs*; it does not say a
   recombination event occurred. In mean field cancellation is arithmetic on a
   signed number stored per cell, and costs nothing. That gap is the whole of
   the pathwise question and cannot be used as evidence about it.
4. **§6 of the supplement, the admissibility of the crystal shift.** The
   argument given there is that the QLE contains only derivatives of $`W`$.
   That is true but weaker than necessary and does not obviously survive the
   move to many bodies. Theorem D2 replaces it: the shift is the identity
   operator, and $`[H, \mathbb{1}] = 0`$ for every Hamiltonian in every
   dimension.

**Concedes.** Theorem D5.1's anonymity licence is proved in the Wigner basis,
where a carrier holds $`(X, p)`$ and the dynamics closes on $`(X, p)`$. It does
*not* automatically extend to any layer whose rates read a collective order
parameter — the $`Z_r`$ of
[`relational_pairing_and_carrier_lock.md`](relational_pairing_and_carrier_lock.md).
That is open item N5 of the supplement, and this note makes it more urgent
rather than less (§12).

---

## 1. Two pairings, and why they must not be conflated

The repository uses the word *pair* for two different objects, and the question
that prompted this note sits exactly on the seam.

| | constituents | separation | carries | defined in |
|---|---|---|---|---|
| **P1 relational pair** | ket leg, bra leg | $`Y = x_1 - x_2`$, out to the coherence length | phase $`\mu = \arg\rho`$ | position pair ladder, phase alignment |
| **P2 sea pair** | positon, negaton | same phase-space cell | nothing; $`\Delta E = 0`$ | sea-dressed §2 |

A P1 pair has two legs at *different positions* and contributes to every
momentum. A P2 pair has two constituents at *one* phase-space point and
contributes to nothing.

The correspondence between them is real but holds **only in aggregate**: the
uniform P2 sea equals the P1 condensate at $`Y = 0`$ (Theorem D2). An
individual P2 bound pair does **not** correspond to an individual $`Y = 0`$ P1
pair. Reading the sea's darkness as a phase cancellation between the two
members of a P2 pair therefore imports P1 structure into a layer that has no
phase in it — the sea-dressed note's dark state is the bookkeeping identity
$`U^+ = U^-`$ at a cell, and nothing more. The phase reading is available, but
it lives in the pair basis and identifies the sea as a whole.

---

## 2. Theorem D1: species and phase are conjugate, not paired

Write the Wigner transform in the leg-separation convention of the interworld
note:

```math
W(X,p) \;=\; \frac{1}{2\pi\hbar}\int dY\;
\underbrace{\rho\!\left(X + \tfrac{Y}{2},\, X - \tfrac{Y}{2}\right)}_{
|\rho|\,e^{i\mu(X,Y)}}\, e^{-ipY/\hbar}.
```

**Theorem D1.** The species of a Wigner world-particle and the phase of a
relational pair are the same degree of freedom expressed in conjugate bases.
The phase is a coordinate on the $`Y`$ axis; the sign is a value on the $`p`$
axis; $`Y`$ and $`p`$ are Fourier conjugate. Consequently no carrier holds
both sharply, and $`\mathrm{sign}\,W(X,p)`$ is not an attribute of any pair but
a functional of the entire $`Y`$-fibre over the midpoint $`X`$.

*Proof.* $`\mu(X,Y)`$ is defined pointwise on the $`(X,Y)`$ chart;
$`\mathrm{sign}\,W(X,p)`$ is the sign of an integral over the whole $`Y`$ line
at fixed $`X`$. A carrier localised in $`Y`$ has no definite $`p`$ and a
carrier localised in $`p`$ has no definite $`Y`$, by the ordinary
Fourier uncertainty relation. $`\square`$

**The decisive case.** Take the cat at rest, $`\psi \propto g(x-d/2) + g(x+d/2)`$
with $`g`$ a real Gaussian. Then $`\rho`$ is real and positive, so $`\mu \equiv 0`$
*everywhere*. Part A measures:

| state | negativity $`\nu`$ | fraction of pairs with $`\cos\mu < 0`$ |
|---|---|---|
| cat at rest, $`\psi`$ real | 0.2924 | 0.0000 |
| same cat, boosted by $`k = 3`$ | 0.2935 | 0.5010 |

**A state with no phase anywhere in the pair basis carries 29 per cent
negative Wigner mass.** The species is manufactured entirely by the kernel
$`e^{-ipY/\hbar}`$ — by the act of integrating the fibre — not by any pair's
phase. Any scheme that assigns a species to a pair from its phase assigns the
wrong label.

The boost supplies the complementary fact: for a plane-wave factor
$`e^{ikx}`$, $`\mu = kY`$ exactly (verified to $`3.1\times10^{-15}`$), a
function of the *relative* position alone. So the kinematic part of the phase
is purely relational, exactly as conjectured; it is the potential that couples
it to $`X`$, and §4 says how.

[![Species and phase](https://raw.githubusercontent.com/billpage/wpmw/output/figures/species_phase_duality.png)](https://raw.githubusercontent.com/billpage/wpmw/output/figures/species_phase_duality.png)

*Left: $`W(X,p)`$ for the cat — the sign structure runs along $`p`$. Centre:
$`|\rho|`$ on the same fibres; three lobes, all real and positive, so
$`\mu \equiv 0`$. Right: how far along the fibre the sign must be read before
it stops changing.*

---

## 3. Theorem D2: the crystal shift is the identity, and the identity is dark

**Theorem D2.** The crystal shift $`W \to W + 2/h`$ is, as an operator,
$`\rho \to \rho + 2\,\mathbb{1}`$. Its Weyl symbol is the constant $`2/h`$, and
it is dynamically inert because $`[H, \mathbb{1}] = 0`$ — for every
Hamiltonian, every potential, every number of particles, every dimension.
In the pair basis $`\mathbb{1}`$ is $`\delta(Y)`$ with $`X`$-independent
amplitude: **the sea sits at exactly zero leg separation.**

*Proof.* The kernel of $`c\,\mathbb{1}`$ is $`c\,\delta(x - x')`$, i.e.
$`c\,\delta(Y)`$, so
$`W_{\rm sea} = (c/2\pi\hbar)\int dY\,\delta(Y)e^{-ipY/\hbar} = c/h`$, constant
in both arguments, and $`c = 2`$ gives the tight Wigner bound $`2/h`$. Inertness
is immediate: the identity commutes with everything. $`\square`$

Part B measures the Weyl symbol of a background of coherence length
$`\epsilon`$ and finds $`0.31830989`$ at $`p = 0`$ against $`2/h = 0.31830989`$
for every $`\epsilon`$, with a plateau half-width tracking $`\hbar/\epsilon`$;
and

```
max| d/dt[rho + 2*1] - d/dt[rho] | = 1.3e-14   (scale |drho/dt| ~ 0.068)
```

for a random Hermitian $`H`$ on 129 sites.

**This answers the question that prompted the note, and more sharply than it
was asked.** The darkness is independent of absolute position and depends only
on relative position — but not because two things cancel. *Nothing is
cancelling.* The commutator is empty. The potential enters the dynamics only
through $`U(X,Y) = V(X+Y/2) - V(X-Y/2)`$, which vanishes identically at
$`Y = 0`$ for every $`X`$ and every $`V`$, so no phase can ever accumulate on a
coincident pair. "Cancellation of phase" is the right instinct about the
*location* and the wrong description of the *mechanism*.

**Corollary D2.1 (the non-compactness obstruction, restated).**
[`open_position_space.md`](open_position_space.md) records that the sea-dressed
layer does not generalise naively to non-compact position space, its constant
background having infinite total. D2 identifies the obstruction exactly:
$`\mathbb{1}`$ is not trace class. The sea exists precisely when the phase-space
volume in view is finite, which is the same condition as the $`2A/h`$ readout
cost of the supplement's §6.

---

## 4. Theorem D3 and Proposition D4: what the sea's polarisation is

**Theorem D3.** For a potential $`V(x) = \sum_q V_q\cos(k_q x + \phi_q)`$ with
$`k_q = 2\pi q/L`$, the misalignment of a relational pair precesses at

```math
\frac{d\mu}{dt} \;=\; -\frac{U(X,Y)}{\hbar}
\;=\; -2\sum_q \Gamma_q(X)\,\sin\!\left(\frac{k_q Y}{2}\right),
\qquad
\Gamma_q(X) = -\frac{V_q}{\hbar}\sin(k_q X + \phi_q).
```

Verified in Part C to $`4\times10^{-15}`$ for one mode and $`8\times10^{-15}`$
for three simultaneous modes.

This is Proposition I2 and Theorem I5 of the interworld note combined and
extended to arbitrarily many modes; the content that is new here is the reading
of the two factors:

- $`Y = 0`$ — **dark**, for every $`X`$ and every mode. This is the sea (D2).
- $`Y = \lambda_q/2`$ — precession is maximal at rate $`2|\Gamma_q(X)|`$, and
  at the maximum of $`|\Gamma_q|`$ that rate is $`2V_q/\hbar`$, which is
  exactly the $`\gamma_{\max}`$ of the supplement's §7.2. **The pathwise
  $`L^1`$ growth rate of the ensemble is a phase precession rate.**
- Periodic in $`Y`$ with period $`2\lambda_q`$, two harmonics — the four
  actions, by Theorem I3.

**Proposition D4.** A background of *finite* coherence length $`\epsilon`$ is
not dark. Its activity, measured as the mass-weighted winding rate
$`A(\epsilon) = \max_X \int du\,|U(X,u)|\rho_\epsilon(u)`$, is linear in
$`\epsilon`$:

| $`\epsilon`$ | $`A(\epsilon)`$ | $`2|V'|_{\max}\sqrt{2/\pi}\,\epsilon`$ | ratio |
|---|---|---|---|
| 1.0000 | 1.786249 | 1.879971 | 0.95015 |
| 0.2500 | 0.468486 | 0.469993 | 0.99679 |
| 0.0625 | 0.117475 | 0.117498 | 0.99980 |

**Corollary D4.1 (the polarisation is derived, not postulated).** The
sea-dressed note lists its polarisation pattern among its postulates: the rate
field $`\Gamma_q(x)`$ and its sign structure enter as assumptions about a
prepared medium. D3 and D4 supply them. The "polarized sea excitations" are
the finite-$`Y`$ tail of the sea, and their winding rate is
$`-2\Gamma_q(X)\sin(k_qY/2)`$ — which is *why* the polarisation matches
$`\Gamma_q`$. What remains postulated is only the amplitude of that tail, not
its shape or its $`X`$-dependence.

[![The sea is the identity](https://raw.githubusercontent.com/billpage/wpmw/output/figures/sea_identity_darkness.png)](https://raw.githubusercontent.com/billpage/wpmw/output/figures/sea_identity_darkness.png)

*Left: the Weyl symbol of a background of coherence length $`\epsilon`$ flattens
to $`2/h`$ as $`\epsilon \to 0`$. Centre: activity linear in $`\epsilon`$; only
$`\epsilon = 0`$ is dark. Right: the precession map $`d\mu/dt`$ over
$`(X,Y)`$, with the dark line at $`Y = 0`$ marked.*

---

## 5. Theorem D5: annihilation is local in the midpoint, non-local in the legs

**Theorem D5.** Annihilation that is exact in $`(X,p)`$ — a positon and a
negaton removed from one phase-space cell — is *non-local in the constituent
leg positions*, with reach set by the coherence length of the state and
**independent of the position lattice spacing**. It cannot be local in both
descriptions, because $`Y`$ and $`p`$ are conjugate.

Part D truncates the $`Y`$-fibre at $`|Y| \le Y_c`$ and re-reads the sign of
$`W`$, weighting by $`|W|`$:

| $`Y_c`$ | sign agreement | | $`d/\sigma`$ | settling $`Y_c`$ | $`d + 4\sigma`$ |
|---|---|---|---|---|---|
| 0.50 | 0.8155 | | 4 | 4.25 | 4.00 |
| 2.00 | 0.7773 | | 8 | 6.25 | 6.00 |
| 4.00 | 0.9604 | | 16 | 8.75 | 10.00 |
| 6.00 | 0.9983 | | | | |
| 8.00 | 1.0000 | | | | |

The sign settles at a leg reach of order $`d + 4\sigma`$ — the full
$`Y`$-support of $`\rho`$, i.e. the coherence length. Nothing in the table
depends on $`\Delta x`$.

**Corollary D5.1 (anonymity is exact).** Two Wigner carriers occupying the
same cell with opposite species are interchangeable *for all future time*, not
merely for present observables. The estimator $`\Psi = \sum_i w_i \delta_{(x_i,p_i)}`$
cannot see $`Y`$; and the QLE closes on $`(X,p)`$, so no $`Y`$-resolved datum
can re-enter the future of the $`(X,p)`$ marginal. Annihilation therefore never
needs to find the original partner.

This is what corrects §7.4 of the supplement. Criterion (ii) computed the
escape time of a *specific* created pair from its position cell and found
$`t_{\rm sep} \propto 1/M`$, concluding that lattice refinement makes matters
worse. That calculation is correct and irrelevant: the specific pair is not
what annihilation consumes.

**Corollary D5.2 (why the crystal is needed anyway).** D5.1 removes the
partner-tracking requirement but not the coincidence requirement. Exact
annihilation still needs two carriers in one cell, and on a continuous position
axis coincidence is measure zero. §7.3 of the supplement stands: the phase-space
crystal is what makes pairwise annihilation definable at all.

---

## 6. Theorem D6: a soft annihilation blob is an imposed coherence length

**Theorem D6.** Convolving $`W`$ with a Gaussian of width $`\sigma_p`$ in
momentum is identically multiplication of $`\rho(X + Y/2, X - Y/2)`$ by
$`\exp(-\sigma_p^2Y^2/2\hbar^2)`$ — that is, the imposition of a leg coherence
length $`\hbar/\sigma_p`$.

| $`\sigma_p`$ | $`\hbar/\sigma_p`$ | max discrepancy | relative | $`\|W\|_1`$ |
|---|---|---|---|---|
| 0.2 | 5.0000 | $`1.1\times10^{-16}`$ | $`4.8\times10^{-16}`$ | 1.4262 |
| 0.5 | 2.0000 | $`1.1\times10^{-16}`$ | $`7.8\times10^{-16}`$ | 1.1073 |
| 1.0 | 1.0000 | $`9.7\times10^{-17}`$ | $`8.6\times10^{-16}`$ | 1.0067 |
| 2.0 | 0.5000 | $`5.6\times10^{-17}`$ | $`7.8\times10^{-16}`$ | 1.0001 |
| — | $`\infty`$ | — | — | 1.5848 |

**Corollary D6.1 (N3 is priced).** The supplement's third open item asks for
the three annihilation regularizations to be compared on one set of axes. Two
of the three now have closed forms. The *bias* of a soft kernel of width
$`\sigma_p`$ is an artificial decoherence of length $`\hbar/\sigma_p`$ —
computable in advance for any state, not merely bounded. The *benefit* is the
$`\|W\|_1`$ collapse in the last column, which is the same quantity the
supplement's N4 predicts should fall under real decoherence. The soft blob and
the decoherence prediction are the same knob read in the two bases.

---

## 7. Why annihilation is the only lever

Recorded here because §8 and §9 rest on it.

**Proposition U1 (growth-rate invariance).** Let $`L`$ be the generator, and
let an unraveling be specified per source cell $`c`$ by a weight-scaling rate
$`a_c`$, hop rates $`\lambda_{rc} \ge 0`$, spawn rates $`\mu_{rc} \ge 0`$ and a
kill rate $`k_c \ge 0`$. Matching the mean field forces
$`L_{cc} = a_c - \sum_r\lambda_{rc} - k_c`$ and
$`|L_{rc}| = \lambda_{rc} + \mu_{rc}`$, so the expected pathwise $`L^1`$ ledger
grows from a carrier at $`c`$ at rate

```math
g_c \;=\; a_c + \sum_{r}\mu_{rc} - k_c \;=\; L_{cc} + \sum_{r \neq c}|L_{rc}|,
```

which depends only on $`L`$. Asymptotically the rate is $`\rho(|L|)`$, the
Perron root of the entrywise absolute value — invariant under positive
diagonal similarity, hence under every guiding function.

The four-action jump stencil has **zero diagonal**, and that zero is forced:
the participant-locality argument of the sea-dressed note's §5 requires
$`a_0 = 0`$ in the focus symbol. So $`g_c = 2|\Gamma_q(x)|`$ with nothing to
subtract, and **no choice of unraveling that is linear in the ensemble can be
$`L^1`$-stationary.** This settles open item N1 of the supplement negatively,
and negatively as a theorem rather than as a failed search. It also disposes of
the hope that broken detailed balance supplies the lever: that asymmetry is in
the *signed* direction of flux, while the ledger counts *gross* traffic, and
$`\sum_r|L_{rc}|`$ is blind to direction by construction.

Part F computes the correct exponent for the supplement's Part F parameters
($`L = 8`$, $`M = 256`$, $`N = 128`$, $`q = 1`$, $`V_q = 1.5`$), including
streaming:

```math
\gamma_{\rm avg} = 1.9099 \;<\; \rho(|L|) = 2.3409 \;<\; \gamma_{\max} = 3.0000
```

converged to $`\pm 1.2\times10^{-15}`$. Neither of the supplement's two bounds
is the governing rate; the Perron root of the full operator is, because
streaming carries carriers through regions of varying $`|\Gamma|`$. §7.2's
burden table should be recomputed against $`e^{2.341\,t}`$.

The escape is therefore forced to be **nonlinear in the ensemble** — two
carriers must meet — which is exactly the status the recombination channel R
already has in the sea-dressed note: quadratic, and permitted because it lies
outside the generator.

**One warning about N1's proposed mechanism.** N1 hopes the exact reverse of
the outward action will serve as its own garbage collection. It will not. The
emitted pair sits at $`(n-q, +w)`$ and $`(n+q, -w)`$ — opposite signs at
*different momenta* — and removing it changes $`\sum_i w_i p_i`$ by
$`2q\,\Delta p\,w`$. **The emitted pair is the momentum kick.** Annihilating it
against its own parent erases the force the four actions exist to produce. Only
same-cell annihilation removes zero of every moment; the two stencils shed the
same $`2|w|`$ of ledger and one of them sheds the physics with it.

---

## 8. Extending the algorithm: an annihilation substep

What the algorithm notes need in order to carry annihilation. The channel is
already present in the sea-dressed layer as R; what follows generalises it to
the bare four-action ensemble and states what changes elsewhere.

### 8.1 State variables

Annihilation cannot be expressed on a signed net field. It requires
**species-resolved occupancies**

```math
U^+_{n,m} \in \mathbb{Z}_{\ge 0}, \qquad
U^-_{n,m} \in \mathbb{Z}_{\ge 0}, \qquad
E_{n,m} = U^+_{n,m} - U^-_{n,m},
```

in place of the single signed array $`E`$. This doubles the occupancy storage
and is the *first* real cost of the scheme — see §9.2.

### 8.2 The substep

**Instantaneous form (the $`\kappa \to \infty`$ limit).** After the jump
substep, in every cell,

```math
r_{n,m} = \min\!\left(U^+_{n,m},\, U^-_{n,m}\right), \qquad
U^\pm_{n,m} \leftarrow U^\pm_{n,m} - r_{n,m}.
```

One vectorised `minimum` over the lattice. No partner indices, no search, no
random numbers.

**Finite-rate form.** $`e \sim \mathrm{Poisson}(\kappa\,U^+U^-\,\delta t / B)`$
per cell, capped by $`\min(U^+, U^-)`$, exactly the R channel of the
sea-dressed note §4. Level 1 is the $`\kappa\to\infty`$ limit; §10 of that note
already measures the convergence $`\kappa = 0 \to 20 \to 200 \to`$ pinned.

### 8.3 Placement in the splitting sequence

Immediately **after** the jump substep and **before** streaming. The jump
substep is what creates opposite-species carriers in a common cell; streaming
is what disperses them. Annihilating between the two captures the maximum
population at the minimum cost, and — because streaming is an integer
permutation that acts identically on $`U^+`$ and $`U^-`$ — the order of the
annihilation and streaming passes changes nothing in the mean field, only the
ledger.

### 8.4 What it preserves

**Proposition D7.** Same-cell annihilation is exactly unbiased. Removing one
positon and one negaton from a single cell changes $`E`$ by zero, hence changes
$`W`$ by zero, hence changes every moment $`\int x^a p^b W`$ by zero. It is
not an approximation and introduces no error term at any order.

Part G confirms this operationally: the reconstructed $`E`$ under annihilation
converges to the annihilation-free mean field as $`\nu^{-1/2}`$ over three
decades (§9.5).

The contrast with the alternatives is the point:

| stencil | $`\Delta\sum w`$ | $`\Delta\sum wp`$ | $`\Delta\sum|w|`$ | verdict |
|---|---|---|---|---|
| same cell, opposite species | 0 | 0 | $`-2|w|`$ | exact |
| exact reverse of the outward action | 0 | $`+2q\Delta p\,w`$ | $`-2|w|`$ | erases the force |
| soft blob of width $`\sigma_p`$ | 0 | 0 | $`-2|w|`$ | decoheres at $`\hbar/\sigma_p`$ (D6) |

### 8.5 What the algorithm notes must add

1. **A species-resolved occupancy declaration** in the state-variable section,
   replacing the signed $`E`$ field, with a note that $`E`$ becomes a derived
   observable.
2. **The substep itself**, §8.2, placed in the splitting sequence per §8.3.
3. **An anonymity clause**, citing D5.1: no partner index is stored, and no
   partner is sought. This is the same conclusion
   [`relational_pairing_and_carrier_lock.md`](relational_pairing_and_carrier_lock.md)
   reached for the pairing vertex, arrived at independently.
4. **A crystal precondition**, citing D5.2: the substep is defined only on a
   discretised position axis. Any variant of the algorithm with continuous
   $`x`$ must use the soft form of §8.2 and accept its D6 bias.
5. **An occupancy precondition**, citing Corollary R6.2 of the supplement: the
   ensemble must satisfy $`\mathcal{N} \gtrsim M`$, i.e. of order one
   world-particle per phase-space cell, or there are no partners to find. §9.5
   measures the degradation at exactly that boundary.
6. **A diagnostic**: the ledger $`\sum(U^+ + U^-)`$ should be reported per run
   alongside the existing worldline invariants. It is the direct measurement of
   whether the scheme is holding.

### 8.6 Interaction with the sea-dressed layer

No conflict. R already *is* this substep, and §8 of the sea-dressed note
already establishes that it is not optional — without it the orphan load runs
away at rate $`\sim|\Gamma|`$, which is Proposition U1 seen from the sea side.
What D5.1 adds there is the licence to keep drawing R anonymously as the
implementation already does, rather than as an approximation to a
partner-resolved process.

---

## 9. Cost and benefit, measured

Part G runs a species-resolved four-action tau-leap on a $`64\times64`$ crystal
with the reference cosine well, with and without the substep.

### 9.1 The benefit

The ledger, over 40 substeps ($`t = 3.18`$), streaming frozen so that the
effect is isolated:

| | initial | final | growth | mean rate | final rate |
|---|---|---|---|---|---|
| no annihilation | $`2.269\times10^{5}`$ | $`2.262\times10^{8}`$ | $`9.97\times10^{2}`$ | 2.169 | 2.469 |
| same-cell annihilation | $`2.000\times10^{5}`$ | $`2.642\times10^{5}`$ | $`1.32`$ | 0.087 | 0.124 |

The ungoverned ledger is climbing toward its asymptote
$`2|\Gamma|_{\max} = 2.996`$ (the pointwise maximum, because streaming is off
in this run; with streaming the invariant is $`\rho(|L|) = 2.341`$). The
governed ledger grows by a factor of $`1.3`$ over the same interval.

Because the cost of a representation is $`N/\|\mu\|_1^2`$, an ungoverned ledger
multiplies the particle requirement by $`e^{2\rho t}`$:

| $`t`$ | particles needed, relative |
|---|---|
| 0.5 | $`1.04\times10^{1}`$ |
| 1.0 | $`1.08\times10^{2}`$ |
| 2.0 | $`1.17\times10^{4}`$ |
| 4.0 | $`1.36\times10^{8}`$ |

### 9.2 The costs, all four of them

1. **Wall clock: negligible.** Measured on the same run, the annihilation pass
   costs $`0.32`$ per cent of the jump pass. This is expected rather than
   fortunate: both are $`O(MN)`$ vectorised sweeps, but the jump pass draws
   Poisson variates per channel per mode while annihilation is a single
   `minimum`. The overhead does not grow with the number of modes, so it
   *falls* as a fraction as the potential gets richer.
2. **Memory: a factor of two.** Two non-negative integer arrays instead of one
   signed array (§8.1).
3. **Variance: a factor of about 1.45, and this is the real cost.** Going
   species-resolved at all means drawing the crossing-conjugate channels
   independently — gross traffic — where a signed-net four-rule Monte Carlo
   draws a single net. §10 of the sea-dressed note measures that penalty
   directly: relL2 $`0.108`$ for the pinned sea against $`0.074`$ for the
   four-rule reference. This cost is incurred by the *decision to track
   species*, not by the annihilation pass, and it is unavoidable if
   annihilation is wanted at all.
4. **A precondition, not a cost: occupancy $`\mathcal{N} \gtrsim M`$.**
   Below about one world-particle per cell there are no partners, and §9.5
   shows the first sign of degradation exactly there.

### 9.3 The crossover

The substep pays for its own wall clock at

```math
t^* \;=\; \frac{\ln(1 + 0.0032)}{2\rho(|L|)} \;=\; 7\times10^{-4},
```

which is under one per cent of a single advection step. Against the variance
cost of item 3 the crossover is
$`\ln(1.45^2)/(2\times2.341) = 0.16`$ — still well inside the first well
period. **There is no regime in which not annihilating is cheaper.**

[![Annihilation cost and benefit](https://raw.githubusercontent.com/billpage/wpmw/output/figures/annihilation_cost_benefit.png)](https://raw.githubusercontent.com/billpage/wpmw/output/figures/annihilation_cost_benefit.png)

*Left: the ledger with and without the substep, log scale. Centre: the
ensemble sweep, showing $`\nu^{-1/2}`$ and no plateau. Right: the benefit
$`e^{2\rho t}`$ against the measured overhead, with the break-even marked.*

### 9.4 What this does *not* buy

Annihilation caps the ledger; it does not make the model efficiently
simulable. The state's own $`\|W\|_1`$ is still the irreducible cost, and by
the supplement's §2 and §4 that is bounded for a cat but multiplicative across
non-Gaussian factors. Annihilation removes the *unraveling's* exponential, not
the *representation's*.

### 9.5 The plateau, or its absence

Open item N2 of the supplement asks for a sweep looking for a critical
ensemble size. Part G sweeps $`\nu`$ over three decades with the substep
active and compares the reconstructed $`E`$ against the exact mean field:

| $`\nu`$ | occupancy per cell | relative L2 error | error $`\times\sqrt{\nu}`$ |
|---|---|---|---|
| 4 000 | 0.98 | 0.78474 | 49.63 |
| 40 000 | 9.77 | 0.20770 | 41.54 |
| 400 000 | 97.66 | 0.06766 | 42.80 |
| 4 000 000 | 976.56 | 0.02134 | 42.67 |

The last column is flat to about 3 per cent for occupancy $`\ge 10`$ — the
error falls as $`\nu^{-1/2}`$ with **no threshold**. The single excursion is at
occupancy $`0.98`$, where the excess is 16 per cent. That is the crystal
condition $`\mathcal{N} \ge M`$ asserting itself, and it is a graceful
degradation rather than the sharp FCIQMC plateau.

**Reading, stated with its limits.** Criterion (i) of §7.4 governs and
criterion (ii) does not, as D5.1 predicts. But this run freezes streaming to
isolate the ledger, uses one mode, and probes only three decades; the FCIQMC
plateau appears in a *correlated* sign structure built up over many
generations, and a longer run with transport is the honest test. What can be
said is that the mechanism criterion (ii) feared — partners escaping before
they can be used — does not operate, because there are no partners in that
sense.

---

## 10. Numerical verification

`src/demo_species_phase_duality.py`, run as
`WPMW_OUTPUT=... PYTHONPATH=src python3 src/demo_species_phase_duality.py`.

| Part | Claim verified |
|---|---|
| A | Quadrature transform against the supplement's published Part A table; Theorem D1, $`\mu \equiv 0`$ with $`\nu = 0.2924`$; $`\mu = kY`$ under boost to $`3.1\times10^{-15}`$ |
| B | Theorem D2: Weyl symbol of $`2\,\mathbb{1}`$ constant at $`2/h`$; $`[H,\mathbb{1}] = 0`$ to $`1.3\times10^{-14}`$ for random Hermitian $`H`$ |
| C | Theorem D3 to $`8\times10^{-15}`$, one and three modes; Proposition D4, activity linear in $`\epsilon`$ to ratio $`0.9998`$ |
| D | Theorem D5: sign settles at leg reach $`d + 4\sigma`$, independent of $`\Delta x`$ |
| E | Theorem D6: momentum blob against separation window, agreement to $`10^{-16}`$ |
| F | $`\rho(|L|) = 2.340863`$ for the supplement's Part F parameters, $`\gamma_{\rm avg} < \rho < \gamma_{\max}`$ |
| G | Ledger with and without the substep; wall-clock overhead 0.32 per cent; ensemble sweep, $`\nu^{-1/2}`$ with no plateau above occupancy 1 |

Figures: `species_phase_duality.png`, `sea_identity_darkness.png`,
`annihilation_cost_benefit.png`, published to the `output` branch.

---

## 11. Summary

1. Species and phase are one degree of freedom in conjugate bases: phase on
   the separation axis, sign on the momentum axis (D1). A cat at rest has no
   phase anywhere and 29 per cent negative Wigner mass.
2. The crystal shift is $`2\,\mathbb{1}`$, and it is dark because
   $`[H,\mathbb{1}] = 0`$, not because anything cancels (D2). In the pair basis
   it is $`\delta(Y)`$ — zero leg separation, for every $`X`$.
3. $`d\mu/dt = -2\sum_q\Gamma_q(X)\sin(k_qY/2)`$ exactly (D3). Zero separation
   is dark; a finite coherence length is active at a rate linear in it (D4);
   the sea's polarisation is thereby derived rather than postulated (D4.1).
4. Annihilation exact in $`(X,p)`$ is non-local over the *coherence length* in
   the leg positions, independent of the lattice spacing (D5); and anonymous
   annihilation is exact for all future time (D5.1). This corrects §7.4 of the
   supplement.
5. A soft annihilation blob of width $`\sigma_p`$ is exactly a coherence
   length $`\hbar/\sigma_p`$ (D6), which prices N3's second option and
   identifies it with N4's decoherence prediction.
6. No unraveling linear in the ensemble is $`L^1`$-stationary (U1); the correct
   exponent is $`\rho(|L|) = 2.341`$, not $`\gamma_{\max}`$ or
   $`\gamma_{\rm avg}`$; and the reverse-action scheme proposed in N1 would
   erase the force rather than collect garbage.
7. The annihilation substep is one vectorised `minimum` costing 0.32 per cent
   of the jump pass, exactly unbiased in every moment (D7), against a benefit
   of $`e^{2\rho t}`$ in particle count. Break-even is at
   $`t^* = 7\times10^{-4}`$.
8. No plateau is found over three decades of ensemble size; the only
   degradation is at occupancy 1, which is the crystal condition.

---

## 12. Open items

1. **N5 becomes urgent.** D5.1's anonymity licence holds in the Wigner basis
   because the dynamics closes on $`(X,p)`$. Any layer whose rates read a
   collective order parameter — the $`Z_r`$ of the phase-alignment layer, or an
   anonymous beat grating built from the whole ensemble — may be sensitive to
   $`Y`$ after all, in which case its carriers are *not* interchangeable and
   the annihilation substep needs re-deriving there. This note does not settle
   it and makes the question sharper.
2. **The plateau test with transport.** §9.5 freezes streaming. The FCIQMC
   plateau lives in a correlated sign structure accumulated over generations;
   a full run with exact integer advection over several well periods, sweeping
   $`\nu`$ down to occupancy $`10^{-1}`$, is the honest version of N2.
3. **The third regularization.** D6 prices the soft blob and D7 the cell-exact
   scheme. The sea-absorbing option costs $`2A/h`$ by the supplement's §6, but
   the three have not been plotted on one set of axes with bias against
   variance, which is what N3 asked for.
4. **Cost in the pair basis.** U1 assumes real signed weights. With
   unit-modulus phases $`\sum_i|w_i| = N`$ is bounded trivially and the burden
   migrates into the variance of $`\sum_i w_i`$. Conjecture: the degradation
   rate of $`N_{\rm eff}`$ in the phase representation is again $`\rho(|L|)`$.
   Cheap to test on the supplement's Part F cat.
5. **Amplitude of the near-dark tail.** D4.1 derives the shape and
   $`X`$-dependence of the sea's polarisation but not its amplitude, which
   remains the sea-dressed note's Level 3 (self-consistent sea).
6. **Multi-body.** D2 is dimension-independent as an operator statement, but
   the four-dimensional note records that the crystal shift does not commute
   with products; how $`\mathbb{1} \otimes \mathbb{1}`$ interacts with the
   two-particle annihilation stencil has not been checked.
