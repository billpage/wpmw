# Species, sectors, and the annihilation substep

**Status.** Analysis note, step 13 of the ladder. Companion demo:
`src/demo_species_sectors_and_annihilation.py`. This note replaces an earlier
step 13 (`species_phase_duality.md`); §0.3 records what changed and why.

The question that started it was B. Page's: *something seems missing in the
relationship between world-particle sign and phase — the dark state of a sea
pair was described as a cancellation of phase, and that cancellation ought to
depend only on the relative position of the pair.* Answering it required first
separating two things the repository had been calling by the same names, and
the answer then reorganised what annihilation is and what it costs.

---

## 0. Orientation

### 0.1 What this note inherits

From [`interworld_coupling.md`](interworld_coupling.md): the potential read as a
leg coupling $`U(X,Y) = V(X + Y/2) - V(X - Y/2)`$ in midpoint and *full*
separation coordinates, its factorisation for one cosine mode (Proposition I2),
the Fourier-duality theorem for the channel set (Theorem I3), and the winding
law $`d\mu/dt = -U/\hbar`$ (Theorem I5). From
[`position_pair_ladder.md`](position_pair_ladder.md): the pair ontology. From
[`phase_alignment_microdynamics.md`](phase_alignment_microdynamics.md): the
misalignment $`\mu`$ as the sole relational datum. From
[`sea_dressed_microdynamics.md`](sea_dressed_microdynamics.md): the neutral sea,
the four actions, the recombination channel R, and the broken detailed balance
of its §6. From
[`../supplement/representation_cost_and_annihilation.md`](../supplement/representation_cost_and_annihilation.md):
the cost functional $`N_{\mathrm{eff}} = N/\|\mu\|_1^2`$, the annihilation
burden, and open items N1–N5.

### 0.2 The two ensembles

Everything below depends on keeping these apart.

| | **E1**, the Wigner ensemble | **E2**, the pair ensemble |
|---|---|---|
| carriers drawn from | the quasi-density $`W(x,p)`$ | the density matrix $`\rho(x_1,x_2)`$ |
| what *positon* means | a carrier of positive sign | the ket leg of a pair |
| what *negaton* means | a carrier of negative sign | the bra leg of a pair |
| relational datum | none; sign is intrinsic to a carrier | phase $`\mu = \arg\rho`$ |
| structure group | $`\mathbb{Z}_2`$ | $`U(1)`$ |
| layers | sea-dressed, four-action, crystal | position pair, phase alignment |

The Weyl transform relates the *represented objects*. Theorem D0 shows it does
not relate the *ensembles*: there is no carrier-level correspondence at all.
"Positon" in the two layers is a homonym, and the two are at most related by a
transformation of the object they represent.

### 0.3 What changed from the earlier step 13, and why

The earlier note is superseded rather than amended. Four things in it were
wrong, and one framing was misleading.

1. **The two seas were conflated.** The earlier §1 asserted that the uniform
   neutral sea equals the $`\delta(Y)`$ condensate. It does not. The neutral
   sea contributes $`0`$ to the observable; the crystal shift contributes
   $`2/h`$. Theorem D2.1 shows both are dark and both lie in the family
   $`c\,\mathbb{1}`$ — they differ only in $`c`$, and their roles are close to
   opposite (§3).
2. **The earlier Corollary D4.1 is withdrawn.** It derived the sea's
   polarisation as "the finite-$`Y`$ tail of the sea." The neutral sea has no
   $`Y`$ coordinate; $`Y`$ is an E2 label. The precession law D3 stands on its
   own, but its application to the reservoir does not.
3. **The earlier D5 overreached.** It called annihilation "non-local in the
   constituent world-particle positions." In E1 there are no constituents:
   cell-exact annihilation is perfectly local. What is non-local is the
   *image* of that event under the Weyl transformation (§5).
4. **The earlier annihilation substep had a defect.** It read
   $`r = \min(U^+, U^-)`$ over cell populations. Including the bound sea drives
   $`S`$ negative while leaving $`E`$ exactly correct — an ill-formed state
   invisible to every mean-field check (§8.2).
5. **The earlier §9 framed annihilation as cost-versus-benefit with a
   break-even time.** That framing belongs to a *spawn* reading of the four
   actions. Under the split reading (§7) recombination is a closure condition
   on a finite reservoir, not an optimisation.

Everything the earlier note got right — the conjugacy theorem, the identity
argument for the crystal shift, the precession law, the reach of annihilation,
the soft-blob identity, and the unraveling invariance — is carried forward
below as D1–D7 and U1.

### 0.4 What this note corrects elsewhere

- **§7.4 of the supplement, criterion (ii).** Its premise, that the original
  partner must be found, is void; Corollary D5.1 licenses anonymous
  annihilation exactly. Criterion (i), occupancy, governs.
- **N1 of the supplement, naming.** The outward action
  $`(n,n) \to (n-q, n+q)`$ is *defocus*; the focus stencil is $`(-1,+2,-1)`$.
- **N1 of the supplement, the mean-field inference.** Boundedness of
  $`\|W(t)\|_1`$ says cancellation occurs, not that a recombination event
  occurred; in mean field cancellation is arithmetic on a stored number.
- **§6 of the supplement, admissibility of the crystal shift.** Replaced by
  Theorem D2, which is stronger and dimension-independent.

---

## 1. Theorem D1: species and phase are conjugate, not paired

Write the Wigner transform in the separation convention of the interworld note:

```math
W(X,p) \;=\; \frac{1}{2\pi\hbar}\int dY\;
\underbrace{\rho\!\left(X + \tfrac{Y}{2},\, X - \tfrac{Y}{2}\right)}_{
|\rho|\,e^{i\mu(X,Y)}}\, e^{-ipY/\hbar}.
```

**Theorem D1.** The species of a Wigner world-particle and the phase of a
relational pair are the same degree of freedom expressed in conjugate bases.
The phase is a coordinate on the $`Y`$ axis; the sign is a value on the $`p`$
axis; $`Y`$ and $`p`$ are Fourier conjugate. No carrier holds both sharply, and
$`\mathrm{sign}\,W(X,p)`$ is not an attribute of any pair but a functional of
the entire $`Y`$-fibre over the midpoint $`X`$.

*Proof.* $`\mu(X,Y)`$ is defined pointwise on the $`(X,Y)`$ chart, while
$`\mathrm{sign}\,W(X,p)`$ is the sign of an integral over the whole $`Y`$ line
at fixed $`X`$. A carrier localised in $`Y`$ has no definite $`p`$ and one
localised in $`p`$ has no definite $`Y`$. $`\square`$

**The decisive case.** For the cat at rest, $`\rho`$ is real and positive, so
$`\mu \equiv 0`$ everywhere. Part A measures:

| state | negativity $`\nu`$ | fraction of pairs with $`\cos\mu < 0`$ |
|---|---|---|
| cat at rest, $`\psi`$ real | 0.2924 | 0.0000 |
| same cat, boosted by $`k = 3`$ | 0.2935 | 0.5010 |

A state with no phase anywhere in the pair basis carries 29 per cent negative
Wigner mass. The species is manufactured entirely by the kernel
$`e^{-ipY/\hbar}`$ — by the act of integrating the fibre. Any scheme that
assigns a species to a pair from its phase assigns the wrong label.

The boost supplies the complement: for a plane-wave factor, $`\mu = kY`$
exactly, verified to $`3.1\times10^{-15}`$ — a function of relative position
alone. The kinematic part of the phase is purely relational; §4 shows it is the
potential that couples it to $`X`$.

[![Species and phase](https://raw.githubusercontent.com/billpage/wpmw/output/figures/species_phase_duality.png)](https://raw.githubusercontent.com/billpage/wpmw/output/figures/species_phase_duality.png)

---

## 2. Theorem D0: the two ensembles have no carrier correspondence

**Theorem D0.** There is no map sending E1 carriers to E2 carriers one at a
time. Part B measures the two species censuses:

| state | E1 negaton mass | E2 weight with $`\cos\mu < 0`$ |
|---|---|---|
| cat at rest, $`d/\sigma = 8`$ | 0.1843 | 0.0000 |
| cat at rest, $`d/\sigma = 16`$ | 0.1964 | 0.0000 |
| Gaussian, boosted $`k = 3`$ | 0.0000 | 0.4928 |
| Gaussian, boosted $`k = 8`$ | 0.0000 | 0.4986 |
| Gaussian at rest | 0.0000 | 0.0000 |

The censuses are anti-correlated. The cat is homogeneous in E2 and 18 per cent
negaton in E1; the boosted Gaussian is homogeneous in E1 — it is Gaussian, so
$`W \ge 0`$ by Hudson's theorem — and 49 per cent negative-phase in E2. A
carrier-level map would have to send a homogeneous population to a
heterogeneous one and back. $`\square`$

**Corollary D0.1.** Results proved in one ensemble do not transfer to the other
without a separate derivation. In particular the anonymity licence of D5.1 is
an E1 result (§12, item 1).

---

## 3. Theorem D2: what is dark, and which sea is which

**Theorem D2.1 (the dark family).** An operator is dark under *every*
Hamiltonian if and only if it is $`c\,\mathbb{1}`$. The commutant of the full
matrix algebra is the scalars. Part C, against twelve random Hermitian
$`H`$ on 64 sites:

| operator | $`\max\|[H,X]\|`$ |
|---|---|
| $`0`$ — neutral sea, bound pairs, $`E = 0`$ | $`0.000\times10^{0}`$ |
| $`2\,\mathbb{1}`$ — crystal shift, $`W \to W + 2/h`$ | $`0.000\times10^{0}`$ |
| $`c\,\mathbb{1}`$, $`c = -0.37`$ | $`0.000\times10^{0}`$ |
| diagonal, non-constant | $`2.474\times10^{1}`$ |
| rank-one projector | $`1.246\times10^{0}`$ |

**Theorem D2 (the crystal shift).** The shift $`W \to W + 2/h`$ is, as an
operator, $`\rho \to \rho + 2\,\mathbb{1}`$; its Weyl symbol is the constant
$`0.31830989 = 2/h`$, verified for every background coherence length; and it is
inert because $`[H,\mathbb{1}] = 0`$ for every Hamiltonian, every potential,
every dimension. In the pair basis $`\mathbb{1}`$ is $`\delta(Y)`$: **the shift
sits at exactly zero leg separation.**

This answers the founding question, and more sharply than it was asked. The
darkness is independent of absolute position and depends only on relative
position — but *nothing cancels*. The commutator is empty. The potential enters
only through $`U(X,Y)`$, which vanishes identically at $`Y = 0`$ for every
$`X`$ and every $`V`$, so no phase can ever accumulate on a coincident pair.
"Cancellation of phase" is the right instinct about the location and the wrong
description of the mechanism.

**Theorem D2.2 (the two seas are not interchangeable).** The neutral sea
($`c = 0`$) and the crystal shift ($`c = 2`$) have near-opposite roles:

| | neutral sea, $`c = 0`$ | crystal shift, $`c = 2`$ |
|---|---|---|
| contribution to $`E`$ | $`0`$ | $`2/h`$ per cell |
| visible in the observable | no | yes |
| dynamically | **live** — the reservoir the four actions draw on | **provably inert** |

Since $`2\,\mathbb{1}`$ commutes with everything, it cannot be pumped, cannot
sustain broken detailed balance and cannot carry gain. **The crystal shift
cannot be the live medium the sea-dressed layer requires.** That layer's
neutral sea is necessarily a distinct object, and D2.1 says the two are the
only members of the dark family the model uses.

**Corollary D2.3 (three roles of $`2/h`$).** The number appears as (i) the
tight bound $`|W| \le 2/h`$; (ii) the crystal shift value; (iii) the reservoir
capacity per cell. Role (iii) borrows role (i), not role (ii): since $`|E|`$
per cell cannot exceed $`\mathcal{N}(2/h)\Delta x\Delta p`$, a sea of that many
pairs per cell can always meet demand by unbinding. §9 measures how loose that
bound is.

**Corollary D2.4 (non-compactness).** [`open_position_space.md`](open_position_space.md)
records that the sea-dressed layer does not generalise naively to non-compact
position space. D2 identifies the obstruction exactly: $`\mathbb{1}`$ is not
trace class.

---

## 4. Theorem D3 and Proposition D4: precession, and why $`Y = 0`$ is special

**Theorem D3.** For $`V(x) = \sum_q V_q\cos(k_q x + \phi_q)`$ with
$`k_q = 2\pi q/L`$,

```math
\frac{d\mu}{dt} \;=\; -\frac{U(X,Y)}{\hbar}
\;=\; -2\sum_q \Gamma_q(X)\,\sin\!\left(\frac{k_q Y}{2}\right),
\qquad
\Gamma_q(X) = -\frac{V_q}{\hbar}\sin(k_q X + \phi_q),
```

verified to $`3.7\times10^{-15}`$ for one mode and $`8.0\times10^{-15}`$ for
three simultaneous modes. This is I2 and I5 combined and extended to
arbitrarily many modes. The reading of the two factors:

- $`Y = 0`$ — dark, for every $`X`$ and every mode. This is the shift (D2).
- $`Y = \lambda_q/2`$ — precession maximal at $`2|\Gamma_q(X)|`$, whose maximum
  $`2V_q/\hbar`$ is exactly the $`\gamma_{\max}`$ of the supplement's §7.2.
  **The pathwise $`L^1`$ growth rate is a phase precession rate.**
- period $`2\lambda_q`$ in $`Y`$, two harmonics — the four actions, by I3.

**Proposition D4.** A background of *finite* coherence length $`\epsilon`$ is
not dark. Its mass-weighted winding rate
$`A(\epsilon) = \max_X \int du\,|U(X,u)|\rho_\epsilon(u)`$ is linear in
$`\epsilon`$: measured against $`2|V'|_{\max}\sqrt{2/\pi}\,\epsilon`$ the ratio
runs 0.9502, 0.9968, 0.9998 as $`\epsilon`$ falls from 1 to 1/16. Only exactly
zero separation is dark.

[![The dark family](https://raw.githubusercontent.com/billpage/wpmw/output/figures/sea_identity_darkness.png)](https://raw.githubusercontent.com/billpage/wpmw/output/figures/sea_identity_darkness.png)

---

## 5. Theorems D5 and D6: the reach of annihilation

**Theorem D5.** A Wigner-basis annihilation event is local in $`(X,p)`$. Its
*image* under the Weyl transformation has support of order the coherence
length of the state, independent of the position lattice spacing. The event
cannot be local in both descriptions, because $`Y`$ and $`p`$ are conjugate.

Part E truncates the $`Y`$-fibre and re-reads the sign of $`W`$:

| $`Y_c`$ | sign agreement | | $`d/\sigma`$ | settling $`Y_c`$ | $`d + 4\sigma`$ |
|---|---|---|---|---|---|
| 0.50 | 0.8157 | | 4 | 4.25 | 4.00 |
| 2.00 | 0.7752 | | 8 | 6.25 | 6.00 |
| 4.00 | 0.9521 | | 16 | 9.00 | 10.00 |
| 6.00 | 0.9983 | | | | |
| 8.00 | 1.0000 | | | | |

The sign settles at a reach of order $`d + 4\sigma`$ — the full $`Y`$-support
of $`\rho`$. Nothing in the table depends on $`\Delta x`$.

**Corollary D5.1 (anonymity is exact).** Two E1 carriers in one cell with
opposite species are interchangeable *for all future time*. The estimator
cannot see $`Y`$, and the QLE closes on $`(X,p)`$, so no $`Y`$-resolved datum
re-enters the future of the $`(X,p)`$ marginal. Annihilation never needs the
original partner. This is what voids criterion (ii) of the supplement's §7.4.

**Corollary D5.2 (the crystal is still needed).** D5.1 removes partner
tracking but not the coincidence requirement. On a continuous position axis
coincidence has measure zero, so §7.3 of the supplement stands.

**Theorem D6 (the soft blob).** Convolving $`W`$ with a Gaussian of width
$`\sigma_p`$ in momentum is identically multiplication of the pair kernel by
$`\exp(-\sigma_p^2Y^2/2\hbar^2)`$ — an imposed leg coherence length
$`\hbar/\sigma_p`$:

| $`\sigma_p`$ | $`\hbar/\sigma_p`$ | max discrepancy | $`\|W\|_1`$ |
|---|---|---|---|
| 0.2 | 5.0000 | $`1.1\times10^{-16}`$ | 1.4262 |
| 0.5 | 2.0000 | $`1.1\times10^{-16}`$ | 1.1073 |
| 1.0 | 1.0000 | $`9.7\times10^{-17}`$ | 1.0067 |
| 2.0 | 0.5000 | $`5.6\times10^{-17}`$ | 1.0001 |
| — | $`\infty`$ | — | 1.5848 |

So cell-exact annihilation is local and exactly unbiased; only the *soft*
variant is genuinely non-local in effect, and its bias is an artificial
decoherence of computable length. This prices the second option of the
supplement's N3 and identifies it with N4's decoherence prediction.

---

## 6. Theorems D8–D10: the excess splits into two sectors

**Theorem D8.** $`\int W(X,p)\,dp = \rho(X,X) = C(X,0)`$, verified to
$`1.6\times10^{-14}`$. The net excess positon count in a position column is
exactly the diagonal of the density matrix, and is therefore non-negative:

```
min over X of Int W dp   =  1.096e-43     <- never negative
min over (X, p) of W     = -0.235646      <- negative pointwise
```

**There is no position at which the world-particle census is net negative.**
Negatons exist pointwise in $`p`$; they can never win a column.

**Theorem D9.** A diagonal $`\rho`$ gives $`W \ge 0`$ with negativity
$`-0.0000\times10^{0}`$; the off-diagonal part gives a Wigner component whose
column sums vanish to $`1.6\times10^{-16}`$. So the excess divides into three
sectors that the dynamics keeps apart:

| sector | is | carries | census |
|---|---|---|---|
| bound sea pairs | the reservoir | nothing | neutral, $`E = 0`$ |
| **unpaired excess positons** | diagonal of $`\rho`$ | the Born density | non-negative per column |
| **column-balanced pairs** | off-diagonal of $`\rho`$ | coherence and momentum transfer | zero net per column |

A column-balanced pair is a positon at $`(m,n_1)`$ and a negaton at
$`(m,n_2)`$, $`n_1 \neq n_2`$ — same position column, zero net probability.
That is precisely what the outward action produces, which is why removing one
against its own parent shifts $`\sum_i w_ip_i`$ by $`2q\Delta p\,w`$: **the
pair is the momentum kick.**

**Theorem D10 (complementary conservation).**

| substep | max $`\|\Delta`$ column sum$`\|`$ | max $`\|\Delta`$ row sum$`\|`$ |
|---|---|---|
| jump (exact, FFT) | $`1.4\times10^{-14}`$ | $`4.865\times10^{0}`$ |
| streaming (integer) | $`3.753\times10^{1}`$ | $`3.6\times10^{-15}`$ |

The jump substep conserves every position-column sum and moves momentum: it is
the force. Streaming conserves every momentum-row sum and moves position: it is
the transport. That is $`2(M+N)`$ exact invariants, available free as
per-substep diagnostics. It also means **all ledger growth lives in the
zero-probability sector**, and the unpaired excess is never touched by the jump
substep at all.

[![Excess sectors](https://raw.githubusercontent.com/billpage/wpmw/output/figures/excess_sectors.png)](https://raw.githubusercontent.com/billpage/wpmw/output/figures/excess_sectors.png)

---

## 7. Proposition U1: why a quadratic channel is unavoidable

**Proposition U1.** For an unraveling specified per source cell by a
weight-scaling rate $`a_c`$, hop rates $`\lambda_{rc}\ge 0`$, spawn rates
$`\mu_{rc}\ge 0`$ and a kill rate $`k_c \ge 0`$, matching the mean field forces
$`L_{cc} = a_c - \sum_r\lambda_{rc} - k_c`$ and
$`|L_{rc}| = \lambda_{rc}+\mu_{rc}`$, so the expected $`L^1`$ ledger grows at

```math
g_c \;=\; a_c + \sum_r \mu_{rc} - k_c \;=\; L_{cc} + \sum_{r\neq c}|L_{rc}|,
```

a function of $`L`$ alone. Asymptotically the rate is $`\rho(|L|)`$, invariant
under positive diagonal similarity and hence under every guiding function —
Part G measures $`3.0000000000`$ for four different guiding functions including
two random lognormals.

The four-action stencil has **zero diagonal**, forced by the $`a_0 = 0`$
condition of the sea-dressed §5, so $`g_c = 2|\Gamma_q(x)|`$ with nothing to
subtract: **no unraveling linear in the ensemble is $`L^1`$-stationary.** This
settles open item N1 negatively as a theorem. Broken detailed balance is not
the lever either: that asymmetry is in the *signed* flux direction, while the
ledger counts *gross* traffic and $`\sum_r|L_{rc}|`$ is blind to direction.

For the supplement's Part F parameters, including streaming,

```math
\gamma_{\mathrm{avg}} = 1.9099 \;<\; \rho(|L|) = 2.3409 \;<\; \gamma_{\max} = 3.0000
```

converged to $`\pm 1.2\times10^{-15}`$. Neither of the supplement's bounds is
the governing rate; §7.2's burden table should be recomputed against
$`e^{2.341\,t}`$.

The escape is therefore forced to be **quadratic in the ensemble** — two
carriers must meet.

---

## 8. Theorems D12–D14: species is an orientation

**Theorem D12.** At the channel level, species conjugation *is* momentum
reflection. With $`R`$ the reflection $`n \to -n`$ and $`C`$ an emission
channel, Part H measures $`\|RCR - \bar{C}\| = 0.000\times10^{0}`$ against
$`\|RCR - C\| = 2.000\times10^{0}`$. The species label is therefore not an
independent charge: it is the **orientation of the momentum transfer**, which
is the winding sense of $`e^{-ipY/\hbar}`$. This explains why the sea-dressed
channel table comes in crossing-conjugate mirror families.

**Theorem D13.** A signed carrier is one whose wave carries a phase offset,
$`s\,e^{-ipY/\hbar} = e^{i\pi(1-s)/2}e^{-ipY/\hbar}`$, and Hermiticity of
$`\rho`$ restricts the offset to $`\{0,\pi\}`$:

| | Hermiticity residual | $`\max|\mathrm{Im}\,W|/\max|\mathrm{Re}\,W|`$ | distinct phases |
|---|---|---|---|
| Hermitian $`\rho`$ | $`2\times10^{-3}`$ | $`1.1\times10^{-16}`$ | 3, i.e. $`\{0,\pm\pi\}`$ |
| $`Y`$-even twist applied | $`4.5\times10^{-1}`$ | $`1.7\times10^{-1}`$ | 7529 |

**Species is a phase — the $`\mathbb{Z}_2`$ residue that the reality of $`W`$
leaves of $`U(1)`$.** This is the structural content of D0: E1 and E2 differ in
structure group, and $`\mathbb{Z}_2 \subset U(1)`$ is an inclusion of groups
not induced by any map of carriers.

**Theorem D14 (the orientation is relational).** Species is *not* the carrier's
own handedness $`\mathrm{sign}(p)`$:

| state | fraction at $`p > 0`$ | negaton mass |
|---|---|---|
| cat at rest | 0.4802 | 0.1843 |
| cat, boost $`k = 8`$ | 1.0000 | 0.1848 |
| cat, boost $`k = 20`$ | 1.0000 | 0.1853 |

A boost puts every carrier at positive momentum without moving the census —
unsurprising, since it is a unitary that shifts $`W`$ rigidly. The chirality is
a property of an *interaction*, never locally readable off one carrier, which
is exactly what D1 requires.

---

## 9. Theorems D15–D16: the actions split and combine, they do not create

**Theorem D15.** The four actions unbind an existing sea pair
($`S \to U^+ + U^-`$, the products landing at $`n \mp q`$ in the same position
column) and recombination re-binds one. Consequently

```math
N_+ = S + U^+ \qquad\text{and}\qquad N_- = S + U^-
```

are each separately and exactly conserved: their difference is the probability
and their sum is the world-particle number. Part I, with recombination off so
the ledger actually grows:

| reading | ledger, start → end | $`2S`$, start → end | $`N_{\rm total}`$ drift |
|---|---|---|---|
| split | $`2.01\times10^{6} \to 1.313\times10^{7}`$ | $`1.920\times10^{8} \to 1.809\times10^{8}`$ | $`0.000\times10^{0}`$ |
| spawn | $`2.01\times10^{6} \to 1.356\times10^{7}`$ | $`0 \to 0`$ | $`5.740\times10^{0}`$ |

The sea drains by precisely the ledger growth. **Species is a conserved charge
that is bound or unbound, never created** — which is what licenses reading it
as an orientation in D12 rather than as a fluctuating label. Both readings give
the same mean field, because in the undepleted limit the split rate reduces to
the spawn rate; they differ in what they conserve.

**Theorem D16 (the sizing floor).** The sea is a finite *local* reservoir. In
units of the Wigner capacity $`B_W = \mathcal{N}(2/h)\Delta x\Delta p`$ the
blocked-split fraction is universal across the lattice:

| $`M_x`$ | $`N_p`$ | $`\beta = 0.01`$ | 0.03 | 0.1 | 0.3 | 1.0 |
|---|---|---|---|---|---|---|
| 32 | 48 | 0.2830 | 0.0121 | 0.0000 | 0.0000 | 0.0000 |
| 64 | 48 | 0.2852 | 0.0143 | 0.0000 | 0.0000 | 0.0000 |
| 32 | 96 | 0.2830 | 0.0121 | 0.0000 | 0.0000 | 0.0000 |
| 64 | 96 | 0.2852 | 0.0143 | 0.0000 | 0.0000 | 0.0000 |

Identical across both lattice dimensions. Below the floor the failure is
*systematic*, not statistical — throttled splits mean the rate law itself is
wrong. So $`B_W`$ is a safe, lattice-independent sufficient bound, and this is
what justifies the sea-dressed note's choice of $`B`$. It is loose by at least
an order of magnitude at these parameters; the exact factor depends on
integration time and on $`\kappa`$, since the true floor is a dynamical balance
between local split flux and local refill rather than a static bound.

---

## 10. The annihilation substep

### 10.1 State variables

Annihilation cannot be expressed on a signed net field. It requires
species-resolved occupancies $`U^\pm_{n,m} \in \mathbb{Z}_{\ge 0}`$ *and* a
separate bound-sea count $`S_{n,m}`$, with $`E = U^+ - U^-`$ a derived
observable.

### 10.2 The substep

**It acts on orphans only.** In every cell,

```math
r_{n,m} = \min\!\left(U^+_{n,m},\, U^-_{n,m}\right), \qquad
U^\pm \leftarrow U^\pm - r, \qquad S \leftarrow S + r .
```

Including $`S`$ in the minimum is a defect, not a variant. With
$`S = 25{,}000`$, $`U^+ = 900`$, $`U^- = 400`$ it removes 25,400 pairs and
drives $`S`$ to $`-400`$ while leaving $`E = 500`$ exactly right — an
ill-formed state invisible to every mean-field check.

The finite-rate form is the R channel of the sea-dressed note; the display
above is its $`\kappa \to \infty`$ limit.

### 10.3 Placement

Immediately **after** the jump substep and **before** streaming. The jump
substep is what puts opposite species in a common cell; streaming disperses
them. Because streaming is an integer permutation acting identically on
$`U^+`$ and $`U^-`$, the order affects the ledger but not the mean field.

### 10.4 What it preserves

**Proposition D7.** Same-cell annihilation is exactly unbiased: removing one
positon and one negaton from a single cell changes $`E`$ by zero, hence $`W`$
by zero, hence every moment by zero. Not an approximation, and no error term at
any order.

| stencil | $`\Delta\sum w`$ | $`\Delta\sum wp`$ | $`\Delta\sum|w|`$ | verdict |
|---|---|---|---|---|
| same cell, opposite species | 0 | 0 | $`-2|w|`$ | exact |
| exact reverse of the outward action | 0 | $`+2q\Delta p\,w`$ | $`-2|w|`$ | erases the force |
| soft blob of width $`\sigma_p`$ | 0 | 0 | $`-2|w|`$ | decoheres at $`\hbar/\sigma_p`$ |
| whole column | 0 | $`\neq 0`$ | $`-2|w|`$ | preserves the Born density exactly |

**Theorem D11.** Every member of the family, from cell-exact to whole-column,
leaves the diagonal of $`\rho`$ invariant to $`2.7\times10^{-13}`$. They differ
only in how much off-diagonal they destroy. The family is a one-parameter
ladder whose invariant is the position observable and whose price is momentum
resolution, with D6 supplying the price and whole-column as the endpoint. A
dual family exists by symmetry, preserving $`|\phi(p)|^2`$ instead.

### 10.5 What the algorithm notes must add

1. Species-resolved occupancies plus a separate $`S`$, with $`E`$ derived.
2. The substep of §10.2, placed per §10.3, acting on orphans only.
3. An **anonymity clause** citing D5.1: no partner index is stored.
4. A **crystal precondition** citing D5.2: continuous-$`x`$ variants must use
   the soft form and accept its D6 bias.
5. A **sizing rule** citing D16 and §11.
6. A **diagnostic**: report the ledger $`\sum(U^+ + U^-)`$ and the $`2(M+N)`$
   marginal invariants of D10 per run.

---

## 11. Theorem D17: adaptive sea allocation

**Theorem D17.** Injecting or removing a *bound* pair at a cell changes $`E`$
by zero and every moment by zero — the same argument as D7. The sea size is
therefore a representational choice, free to vary per cell and per step, and
does not have to respect the $`N_\pm`$ conservation of D15, which is a property
of the dynamics rather than of the bookkeeping.

**The strategy.** Before each jump substep, size each cell against its own
expected split flux rather than against the global bound:

```math
S^{\rm target}_{n,m} = \left\lceil
\eta\,|\Gamma(x_m)|\,\Delta t\,(U^+ + U^-)_{n,m}\right\rceil + S_0,
```

injecting the shortfall and releasing anything above $`2S^{\rm target}`$.

**Proposition D18.** Part J, with $`\mathcal{N} = 2\times10^{6}`$:

| scheme | blocked | rel L2 error | peak $`N_{\rm total}`$ | vs uniform |
|---|---|---|---|---|
| $`M_x = 32`$, $`N_p = 48`$, uniform $`B_W`$ | 0.0000 | 0.00932 | $`1.940\times10^{8}`$ | 1.000 |
| adaptive, $`\eta = 30`$ | 0.0000 | 0.00932 | $`2.494\times10^{6}`$ | 0.013 |
| adaptive, $`\eta = 100`$ | 0.0000 | 0.00932 | $`3.544\times10^{6}`$ | 0.018 |
| $`M_x = 64`$, $`N_p = 96`$, uniform $`B_W`$ | 0.0000 | 0.01279 | $`3.860\times10^{8}`$ | 1.000 |
| adaptive, $`\eta = 30`$ | 0.0000 | 0.01279 | $`2.534\times10^{6}`$ | 0.007 |
| adaptive, $`\eta = 100`$ | 0.0000 | 0.01279 | $`3.582\times10^{6}`$ | 0.009 |

Identical error to five digits, zero blocked splits, and one to two orders of
magnitude fewer world-particles. **A uniform sea costs
$`O(N_p\,\mathcal{N})`$**, because it stocks every cell against the global
bound; **an adaptive sea costs $`O(\text{ledger})`$**, because it stocks each
cell against that cell's own flux. The saving is of order the number of
momentum cells, and grows with the sparsity of the state in phase space.

**A design caution.** The recombination rate must not reference $`B`$. If it is
written $`\propto U^+U^-/B`$, then one parameter sets both reservoir size and
recombination strength: enlarging the sea silently weakens $`R`$, inflating the
ledger, and a spurious numerical optimum appears near $`\beta = 3`$. With a
density-normalised rate the error is flat above the floor (0.01755 at every
$`\beta`$ from 1 to 300) and cost is linear in $`B`$. There is no genuine
optimum above the floor — take the floor with headroom.

**Recombination remains mandatory.** With the channel off, the ledger runs
$`1.356\times10^{7}`$ against $`2.006\times10^{6}`$ over the same interval,
climbing toward $`\rho(|L|) = 2.341`$, and the sea drains monotonically.

[![Sea sizing](https://raw.githubusercontent.com/billpage/wpmw/output/figures/sea_sizing_and_annihilation.png)](https://raw.githubusercontent.com/billpage/wpmw/output/figures/sea_sizing_and_annihilation.png)

---

## 12. Numerical verification

`src/demo_species_sectors_and_annihilation.py`, run as
`WPMW_OUTPUT=... PYTHONPATH=src python3 src/demo_species_sectors_and_annihilation.py`.

| Part | Claim verified |
|---|---|
| A | quadrature against the supplement's Part A table; D1, $`\mu\equiv 0`$ with $`\nu = 0.2924`$; $`\mu = kY`$ to $`3.1\times10^{-15}`$ |
| B | D0, the anti-correlated census table |
| C | D2, D2.1, D2.2: dark iff $`c\,\mathbb{1}`$; Weyl symbol of $`2\,\mathbb{1}`$ constant at $`2/h`$ |
| D | D3 to $`8.0\times10^{-15}`$ on three modes; D4 linear in $`\epsilon`$ to ratio 0.9998 |
| E | D5 reach $`d + 4\sigma`$; D6 to $`1.1\times10^{-16}`$ |
| F | D8 to $`1.6\times10^{-14}`$; D9 balance to $`1.6\times10^{-16}`$; D10 |
| G | U1; $`\rho(|L|) = 2.3409`$, unchanged by four guiding functions |
| H | D12 exact; D13; D14 boost table |
| I | D15, $`N_{\rm total}`$ drift $`0.000\times10^{0}`$ against 5.740; D16 floor |
| J | D17, D18: identical error at 0.7–1.8 per cent of the uniform particle count |

Figures: `species_phase_duality.png`, `sea_identity_darkness.png`,
`excess_sectors.png`, `sea_sizing_and_annihilation.png`, on the `output` branch.

---

## 13. Open items

1. **The E1/E2 relationship is the structural question.** D0 shows the
   ensembles do not correspond carrier by carrier and D13 identifies the
   difference as one of structure group. The ladder currently reads as a single
   sequence sliding between the two pictures. Whether the Wigner microdynamics
   and the pair microdynamics are the same model, two models related by a
   transformation of the object, or two models full stop, is not settled here —
   and results proved in one do not transfer without a separate derivation.
   This subsumes the earlier question about $`Z_r`$, which was ill-posed:
   $`Z_r`$ sums over E2 carriers while annihilation acts on E1 carriers.
2. **The adaptive scheme's lookahead.** §11 sizes a cell from its *current*
   ledger, which works when demand varies smoothly. A sharp transient — a
   packet arriving in a previously empty region — could outrun the allocation
   for a step. The safety factor absorbs it over the tested range
   $`\eta \in [30, 300]`$; this is a tested range, not a guarantee.
3. **The plateau test with transport.** The ensemble sweeps here freeze or
   simplify streaming. The FCIQMC plateau lives in correlated sign structure
   accumulated over generations; a full run with exact integer advection over
   several well periods is the honest version of the supplement's N2.
4. **The third regularization.** D6 prices the soft blob, D7 and D11 the
   cell-exact and whole-column schemes. The sea-absorbing option costs
   $`4A/h`$ for a neutral sea against $`2A/h`$ for the shift, but the three
   have not been plotted on one set of axes with bias against variance.
5. **Cost in the pair basis.** U1 assumes real signed weights. With
   unit-modulus phases $`\sum_i|w_i| = N`$ is bounded trivially and the burden
   migrates into the variance of $`\sum_i w_i`$. Conjecture: the degradation
   rate of $`N_{\mathrm{eff}}`$ is again $`\rho(|L|)`$.
6. **Amplitude of the near-dark tail.** D4 gives the shape and $`X`$-dependence
   of a finite-$`\epsilon`$ background's activity but not the amplitude of any
   such tail in the physical sea, which remains the sea-dressed note's Level 3.
7. **Many bodies.** D2 is dimension-independent as an operator statement, but
   how $`\mathbb{1}\otimes\mathbb{1}`$ interacts with a two-particle
   annihilation stencil has not been checked.
