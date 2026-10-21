# The Jordan lifting: two signed species, Pawula's theorem, and pair chemistry as a sink

**A five-line route from a signed transition kernel to a positive two-species
Markov process, the exact price of that route, why Pawula's theorem then
constrains the sum of the species and not their difference, and what a P–N
binding well deeper than the like-pair wells does and does not buy.**

---

## 0. Status and provenance

Supplement prompted by the calls of 13 and 20 October 2026 with D. Cyganski
(summary kept outside the repository), where Bill Page asked for a shorter and
more formal route to the two-species "$`W`$ and $`S`$" result, for an answer to
Takabayasi's 1954 verdict that a signed transition probability "cannot be taken
as a real one", and whether Cyganski's bound-pair chemistry — a P–N well deeper
than the P–P and N–N wells — gets round any of the project's no-go results.
Companion demo: `src/demo_jordan_lifting.py`, which produces every number below
and the figure.

It builds on four earlier results and repeats none of them:

- [`takabayasi_1954_stochastic_picture.md`](takabayasi_1954_stochastic_picture.md)
  (**T**): Lemma T1 ($`J`$ is the four-action stencil), Proposition T3 (no
  one-body Markov reading) and Theorem T5 (negativity of the finite-time kernel
  is forced). T §5.2 answered Takabayasi by relocating the sign into species
  and channel direction, *given a pinned sea*.
- [`representation_cost_and_annihilation.md`](representation_cost_and_annihilation.md)
  §7.2 (**R**): the pathwise $`L^1`$ mass of a non-annihilating unraveling grows
  at $`2\lvert\Gamma_q(x)\rvert`$ per mode.
- [`../analysis/sea_depletion.md`](../analysis/sea_depletion.md) (step 24):
  Proposition V2 (the minimal ensemble, $`N = \lVert W\rVert_1`$, is the
  cheapest), V9 (a Volterra sea cannot repair depletion) and V10 (pair
  collisions relax the sea).
- [`../analysis/force_blind_sea.md`](../analysis/force_blind_sea.md)
  Proposition Q8 (the contact sink) and
  [`world_crossings.md`](world_crossings.md) Lemma CX3 (chance coincidence has
  codimension $`n - 1`$).

New here: the lifting itself as a theorem about any signed generator (HJ1),
with no sea in it; its exact minimal cost, which is also a proof that no
*linear* positive process can remove worlds (HJ2, HJ2.1); the moment statement
that Pawula's theorem acts on $`S`$ and not on $`W`$ (HJ3); the generator-only
uniqueness illustrated by Cyganski's up/down actions (HJ4); and pair chemistry
as a finite-temperature contact sink (HJ5). Labels HJ1 to HJ5 ("Hahn–Jordan")
and open items HJ-SP1 to HJ-SP3 belong to this supplement.

---

## 1. Setting

Units $`\hbar = m = 1`$ in the demo; symbols keep $`\hbar`$ here.

- $`f(x,p,t)`$, the Wigner function, obeys Takabayasi's quantum Liouville
  equation (T §1, his (3.5)–(3.6)):

```math
\partial_t f + \frac{p}{m}\,\partial_x f \;=\; \Lambda[f] \;=\; \int J(x,\, p - p')\, f(x,p')\, dp' ,
\qquad \int J(x,q)\, dq = 0 .
```

- On a momentum lattice, or as an operator, write the right side as
  $`L f`$. $`L`$ is real; its **off-diagonal** entries $`L_{mn}`$ ($`m \ne n`$,
  weight sent from row $`n`$ to row $`m`$) have both signs; for Takabayasi's
  kernel its diagonal vanishes and every column sums to zero (T3).
- A **Metzler** matrix has non-negative off-diagonal entries; it generates a
  positivity-preserving linear flow. A Markov generator is a Metzler matrix with
  zero column sums. Free streaming, a permutation of cells on a lattice or a
  measure-preserving flow in the continuum, is Markov.
- The **negative mass** of column $`n`$ is
  $`\nu_n(L) = \sum_{m\ne n} \max(-L_{mn}, 0)`$; for a kernel it is
  $`\nu(x) = \int \max(-J(x,q),0)\,dq`$. Since $`\int J = 0`$,
  $`2\nu(x) = \lVert J(x,\cdot)\rVert_{TV} = \int \lvert J(x,q)\rvert\,dq`$.
- Two species of worlds with densities $`u^\pm \ge 0`$; the observable and the
  unsigned count are

```math
W = u^+ - u^-, \qquad S = u^+ + u^- .
```

  In ledger language $`S`$ integrates to the unpaired-world count $`N`$.

---

## 2. Theorem HJ1: the lifting

**Theorem HJ1 (Jordan lifting).** *Let $`L = A - B`$ with $`A`$ Metzler and
$`B`$ entrywise non-negative. Then the block generator*

```math
\partial_t \begin{pmatrix} u^+ \\ u^- \end{pmatrix}
= \begin{pmatrix} A & B \\ B & A \end{pmatrix}
\begin{pmatrix} u^+ \\ u^- \end{pmatrix}
```

*is Metzler, so it preserves $`u^\pm \ge 0`$, and*

```math
\partial_t W = L\,W, \qquad \partial_t S = (A + B)\, S .
```

*For a kernel, the* Jordan lifting *takes $`A`$ and $`B`$ to be the positive
and negative parts $`J^\pm = \max(\pm J, 0)`$ of the Hahn–Jordan decomposition
$`J = J^+ - J^-`$, plus streaming in $`A`$; then $`S`$ evolves with the kernel
$`\lvert J\rvert`$. Read as a process: a world of sign $`\sigma`$ at $`(x,p)`$
streams freely and, at total rate $`\lVert J(x,\cdot)\rVert_{TV}`$, emits a
child at $`p + q`$ with $`q`$ drawn from $`\lvert J(x,q)\rvert/\lVert J\rVert_{TV}`$
and sign $`\sigma\,\mathrm{sgn}\,J(x,q)`$; the parent persists.*

*Proof.* Subtract and add the two rows: $`\partial_t(u^+ - u^-) = (A - B)(u^+ - u^-)`$
and $`\partial_t(u^+ + u^-) = (A + B)(u^+ + u^-)`$. Off-diagonal blocks are
$`B \ge 0`$, diagonal blocks Metzler. $`\square`$

That is the whole argument; Part A checks it symbolically. For one Fourier
mode (Lemma T1), $`J = \Gamma_q(x)\,[\delta(q + \xi) - \delta(q - \xi)]`$ with
$`\xi = \hbar k_q/2`$, so an event emits a child at $`p - \xi`$ with sign
$`\sigma\,\mathrm{sgn}\,\Gamma_q`$ and one at $`p + \xi`$ with the opposite
sign. Emitting both at once at rate $`\lvert\Gamma_q\rvert`$ is the
**emissive realisation** of a residual event
([glossary](../GLOSSARY.md#absorptive-emissive)), with Takabayasi's whole
kernel in place of the compensated residual. HJ1 is therefore the compensated algorithm's
emission with every interpretive layer removed: no sea, no clock, no
allocation — only the requirement that the two blocks be non-negative.

**Relation to T §5.2.** T answered Takabayasi with three moves (two species,
two directions, sea traffic) made through the pinned sea of
[`../analysis/sea_dressed_microdynamics.md`](../analysis/sea_dressed_microdynamics.md).
HJ1 shows the first two moves need no sea at all: the sign of $`J`$ becomes the
sign of a child, and the zero exit rate of T3 becomes a parent that never
leaves. The sea is needed for something else — paying for $`S`$ (§3, §6).

---

## 3. Theorem HJ2: the price, and why no linear process can remove a world

**Theorem HJ2 (minimal cost).** *Let $`L`$ have zero column sums. For every
lifting $`L = A - B`$ of HJ1, column $`n`$ of $`A + B`$ sums to at least
$`2\nu_n(L)`$, with equality exactly for the Jordan lifting
($`B_{nn} = 0`$, $`B_{mn} = \max(-L_{mn},0)`$). Hence the expected unsigned
count obeys*

```math
\frac{d}{dt}\int S \;\ge\; \int 2\nu(x)\, S(x,p)\,dx\,dp \;=\; \int \lVert J(x,\cdot)\rVert_{TV}\, S\,dx\,dp \;\ge\; 0 ,
```

*with equality for the Jordan lifting.*

*Proof.* For $`m \ne n`$, $`A_{mn} + B_{mn} \ge \lvert A_{mn} - B_{mn}\rvert = \lvert L_{mn}\rvert`$
because both are non-negative. On the diagonal $`A_{nn} = L_{nn} + B_{nn}`$ and
$`B_{nn} \ge 0`$. Summing column $`n`$,
$`\sum_m (A+B)_{mn} \ge L_{nn} + \sum_{m\ne n}\lvert L_{mn}\rvert`$, and with
$`L_{nn} = -\sum_{m\ne n} L_{mn}`$ the right side is
$`\sum_{m \ne n}(\lvert L_{mn}\rvert - L_{mn}) = 2\nu_n`$. Equality needs
$`B_{nn} = 0`$ and $`\min(A_{mn}, B_{mn}) = 0`$ for every $`m \ne n`$, which is
the Jordan split. $`\square`$

Three readings.

- **The price of the sign is exactly the negative mass.** A Markov part
  (streaming, the classical force, any positive transition) has
  $`\nu = 0`$ and costs nothing; only the negative off-diagonal mass of the
  generator does. For one cosine mode $`2\nu = 2\lvert\Gamma_q(x)\rvert`$,
  which is R §7.2's measured growth rate $`\gamma`$, now a lower bound over
  *all* positive linear representations rather than a property of one. Part B
  checks the equality to $`9\times10^{-15}`$ and 2000 random non-minimal splits
  against the bound; Part C checks it along a phase-space run (§8).
- **The rate exists exactly when $`J`$ is a finite measure.** For a periodic
  potential, or one whose Fourier transform is integrable, $`\lVert J\rVert_{TV}`$
  is bounded uniformly in $`x`$ by a multiple of $`\int\lvert\tilde V\rvert/\hbar`$
  (T §1 with the one-dimensional prefactor of T §7). For a
  polynomial it is infinite (Theorem E7, quoted in
  [`limkumnerd_weighted_paths.md`](limkumnerd_weighted_paths.md) Proposition W2),
  and HJ1 applies only to the reach-windowed residual of the compensated split
  ([`what_the_reach_is.md`](what_the_reach_is.md)).
- **What Takabayasi's verdict becomes.** His picture is the formal limit with
  no traffic at all, $`A + B = 0`$ off-diagonal; HJ2 says any *positive*
  picture must pay at least $`\lVert J\rVert_{TV}`$ per world per unit time.
  "Unrealistic" is replaced by "costly", and the cost is now a number.

**Corollary HJ2.1 (no linear sink).** *In every positive linear lifting,
$`\int S`$ is non-decreasing in expectation. Any mechanism that keeps the
unpaired-world count bounded while preserving $`W`$ — annihilation, absorption,
the contact sink, K1 resampling — must be nonlinear in the populations, at
least bilinear.*

*Proof.* Column sums of $`A + B`$ are $`\ge 2\nu_n \ge 0`$. $`\square`$

A bilinear rate needs two worlds in one cell, so every sink inherits the
coincidence problem of Lemma CX3: in $`2dN`$ phase-space dimensions two
independent worlds share a cell with probability falling exponentially in
$`dN`$. HJ2.1 is why that problem cannot be designed away by a cleverer linear
rule; only worlds that are *born* together (CX5, the pair at an event) escape
it.

---

## 4. Corollary HJ3: Pawula's theorem acts on $`S`$, not on $`W`$

Let $`M_k[K] = \int q^k K(x,q)\,dq`$ be the jump moments of a kernel.
Pawula (1967) showed that for a *non-negative* kernel the Kramers–Moyal series
either stops at $`k = 2`$ or never stops; in particular odd moments cannot
survive without $`M_2 > 0`$. That theorem postdates Takabayasi and is the
precise form of his complaint.

**Corollary HJ3.** *Under the Jordan lifting of an odd kernel
($`J(x,-q) = -J(x,q)`$), $`W`$ evolves with $`J`$, whose even moments vanish
(Takabayasi (3.17)), and $`S`$ evolves with $`\lvert J\rvert`$, which is even,
non-negative, and has $`M_2[\lvert J\rvert] > 0`$. Pawula's theorem applies to
the $`S`$ equation, which diffuses; it does not apply to the $`W`$ equation,
which is not the evolution of any non-negative density.*

For one mode (Part A, $`a = \xi`$):

| kernel | $`M_0`$ | $`M_1`$ | $`M_2`$ | $`M_3`$ | $`M_4`$ |
| --- | --- | --- | --- | --- | --- |
| $`J`$ (for $`W`$) | 0 | $`-2\Gamma a`$ | 0 | $`-2\Gamma a^3`$ | 0 |
| $`\lvert J\rvert`$ (for $`S`$) | $`2\lvert\Gamma\rvert`$ | 0 | $`2\lvert\Gamma\rvert a^2`$ | 0 | $`2\lvert\Gamma\rvert a^4`$ |

With $`\Gamma_q = -(V_q/\hbar)\sin(k_q x + \phi_q)`$ and $`a = \hbar k_q/2`$,
$`M_1[J] = -V'(x)`$, Ehrenfest's force. The diffusion that Pawula's theorem
demands is all in $`S`$, which nobody observes; this is the formal content of
the slogan that the sea converts a positive process into a signed one, and of
Proposition X1 of
[`poisson_kicks_and_pair_branching.md`](poisson_kicks_and_pair_branching.md)
read backwards: every positive *one-species* rule diffuses, and a two-species
rule puts the diffusion where it cannot be seen.

---

## 5. Proposition HJ4: the lifting fixes the generator, not the process

HJ1 constrains only $`A - B`$. Two kinds of freedom remain, and both have
appeared in the ladder under other names.

- **$`W`$-neutral terms.** Anything added identically to the $`u^+`$ and
  $`u^-`$ equations (or, nonlinearly, any term that changes $`u^+`$ and
  $`u^-`$ by equal amounts in each cell) leaves $`\partial_t W`$ unchanged and
  changes only $`S`$. Contact recombination, absorption, K1 resampling and the
  chemistry of §6 are all of this kind; so is Cyganski's "gauge freedom" among
  four-action models with one linearisation, and the split gauge of step 17.
- **Correlations between events.** Two processes with the same generator can
  differ in every moment beyond the first.

**Proposition HJ4.** *Cyganski's four-action defocus at cell $`n`$, one event
moving one world from $`n`$ to $`n-1`$ and one to $`n+1`$ at rate $`r`$, and
his independent "up" action, $`n \to n+1`$, and "down" action, $`n \to n-1`$, at rate
$`r`$ each, have the same mean drift $`(r, -2r, r)`$ on
$`(u_{n-1}, u_n, u_{n+1})`$ but different second moments:
$`\mathrm{Cov}(\Delta u_{n-1}, \Delta u_{n+1}) = r`$ against $`0`$, and
$`\mathrm{Var}(\Delta u_n) = 4r`$ against $`2r`$.*

*Proof.* Part E, by enumeration of the jump vectors. $`\square`$

So "the up/down actions decompose the four-action model" is true of the
generator and false of the process. Both reproduce the QLE, which is a mean
equation; they are different unravelings with different noise, as T §5.3
already said of the pinned sea ("uniqueness is of the generator, not the
process").

---

## 6. Proposition HJ5: pair chemistry as a finite-temperature contact sink

### 6.1 Statement

Cyganski's bound-pair picture has three kinds of pair in a cell: P–N, which in
this project's language is a sea (aligned) pair, and P–P and N–N, which are two
unpaired worlds of one sign held together so that a kick can carry a whole
photon. The reaction is

```math
\mathrm{PP} + \mathrm{NN} \;\rightleftharpoons\; 2\,\mathrm{PN},
\qquad K = \frac{n_{PN}^2}{n_{PP}\,n_{NN}} = g\, e^{\varepsilon/kT},
\qquad \varepsilon = 2D_{PN} - D_{PP} - D_{NN},
```

with $`D`$ the well depths of Morse-like potentials (Morse 1929), $`kT`$ the
effective temperature of whatever dissociates pairs, and $`g`$ a ratio of
internal partition functions (for homonuclear pairs the symmetry numbers make
$`g \to 4`$ at high temperature; the demo sets $`g = 1`$). Write
$`n_P = 2n_{PP} + n_{PN}`$ and $`n_N = 2n_{NN} + n_{PN}`$ for the constituents
in the cell, so $`W = n_P - n_N`$ and the unpaired-world count is
$`N = 2(n_{PP} + n_{NN})`$.

**Proposition HJ5.** *(a) The reaction conserves $`n_P`$ and $`n_N`$ in every
cell, hence $`W`$; it is $`W`$-neutral in the sense of §5. (b) With
mass-action rates obeying detailed balance, the free energy
$`G = \sum_s n_s(\ln(n_s/n_s^{\ast}) - 1) + n_s^{\ast}`$ never increases
(Horn and Jackson 1972). (c) At equilibrium
$`N \ge \lvert W\rvert`$, and $`N \to \lvert W\rvert`$ as $`K \to \infty`$: the
cell holds only the majority sign, which is the minimal ensemble of
Proposition V2 — the contact sink of Proposition Q8 is the
$`\varepsilon/kT \to \infty`$ limit of this chemistry. (d) For large $`K`$ the
inflation is*

```math
N - \lvert W\rvert \;\simeq\;
\begin{cases}
8\,n_{\min}^2 / (K\,\lvert W\rvert), & \lvert W\rvert \gg n_{\min}/\sqrt K ,\\[4pt]
4\,n / \sqrt K, & W = 0,\ n_P = n_N = n .
\end{cases}
```

*(e) The project's own sink reaction $`\mathrm{P} + \mathrm{N} \rightleftharpoons \mathrm{PN}`$,
with $`K' = n_{PN}/(n_P^{\rm free}\, n_N^{\rm free})`$, has inflation
$`2n_{\min}/(K'\lvert W\rvert)`$ and $`2\sqrt{n/K'}`$ in the same two regimes.*

*Proof.* (a) by stoichiometry. (b) is the standard Lyapunov function of a
detailed-balanced mass-action network. (c), (d), (e): solve the quadratic
equilibrium (Part D gives the roots symbolically), take $`K \to \infty`$, and
expand the small root. $`\square`$

Part D tabulates (c)–(e) against the exact equilibria, and integrates the
relaxation at $`\varepsilon/kT = 6`$: $`n_P`$, $`n_N`$ and $`W`$ stay fixed to
$`10^{-13}`$ while $`G`$ falls from 160.4 to zero monotonically.

| $`\varepsilon/kT`$ | $`n_{PP}`$ | $`n_{NN}`$ | $`n_{PN}`$ | $`N - \lvert W\rvert`$ | $`8n_N^2/(K\lvert W\rvert)`$ |
| --- | --- | --- | --- | --- | --- |
| 0 | 37.53 | 32.53 | 34.94 | 130.1 | — |
| 5 | 10.26 | 5.259 | 89.48 | 21.04 | 53.9 |
| 10 | 5.089 | 0.0889 | 99.82 | 0.356 | 0.363 |
| 15 | 5.001 | 0.0006 | 99.999 | 0.0024 | 0.0024 |

*(Cell with $`n_P = 110`$, $`n_N = 100`$, so $`\lvert W\rvert = 10`$.)*

One difference between the two reactions matters for Cyganski's simulations.
$`K'`$ has the dimension of an inverse volume — mass action for two-into-one
scales with the reaction volume — so for $`\mathrm{P} + \mathrm{N}`$ a smaller
arena and a deeper well are interchangeable, which is why he could restore
recombination by shrinking the arena inside a cell. For
$`\mathrm{PP} + \mathrm{NN} \rightleftharpoons 2\mathrm{PN}`$ the volumes cancel:
the arena sets how fast equilibrium is reached, and only $`\varepsilon/kT`$
sets where it is.

### 6.2 Is "P–N deeper than both" the right condition?

It is sufficient, not necessary. Only $`\varepsilon`$ enters, and what matters
is $`\varepsilon \gg kT`$: "deeper than both" guarantees $`\varepsilon > 0`$
but not that it is large, and $`\varepsilon \gg kT`$ is possible without it.
Two real reactions (bond enthalpies; Part D):

| reaction | P–N deeper than both? | $`\varepsilon`$ | $`\varepsilon/kT`$ at 300 K |
| --- | --- | --- | --- |
| H₂ + I₂ → 2 HI | no (H–I 3.09 eV, H–H 4.52 eV) | +0.12 eV | 4.6 |
| C₂ + O₂ → 2 CO | yes | +10.97 eV | 424 |

So Morse parameters for which the model is a good sink certainly exist; the
HI system, the call's example, is the marginal case. It is also not an
elementary four-centre exchange: Sullivan (1967) showed that the rate follows
the iodine *atom* concentration, so the real route runs through unpaired
constituents.

### 6.3 What it does not do

| no-go | does HJ5 help? | reason |
| --- | --- | --- |
| Pawula (one positive species must diffuse) | not needed | HJ1 and HJ3 escape it by two species; chemistry is $`W`$-neutral and acts only on $`S`$. |
| Takabayasi's "unrealistic" | yes | It is the missing sink: a positive, detailed-balanced, $`W`$-neutral process that bounds the $`S`$ growth HJ2 makes unavoidable. |
| HJ2.1 (no linear sink) | consistent | Its rates are bilinear, as HJ2.1 requires. |
| Coincidence (CX3) | no | Bilinear rates need both reactants in one cell; the yield falls exponentially with $`dN`$ like every sink. Only pairs born bound avoid it. |
| Depletion repair (V9, V10) | no | It conserves $`n_P`$ and $`n_N`$ cell by cell, so it cannot move sea between rows. Repair needs transport along rows — V10's pair collisions. Whether chemistry *coupled* to whole-pair kicks relaxes rows is HJ-SP1. |
| Linearity at a shallow sea | unchanged | $`W`$ is fixed by the lifting; the sink touches only $`S`$. |

---

## 7. Takabayasi's verdict, with what is available now

Takabayasi's §3(a) wrote the QLE as a jump process with an odd,
non-normalisable "transition probability" and closed: *"this picture cannot be
taken as a real one"*; later, of hidden disturbances acting on a particle,
*"the features of stochastic transitions … do not allow to construct any such
model realistically."* His objection is to a *single* ensemble of particles
kicked by a disturbance. The results below, none available in 1954, turn it
into three statements.

| result | what it establishes |
| --- | --- |
| Pawula (1967) | Why one positive species cannot do it: a positive kernel with surviving odd moments must have $`M_2 > 0`$. |
| Ikeda, Nagasawa and Watanabe (1968–69) | A general theory of branching Markov processes, in which a parent persists and emits — the class HJ1 lives in. |
| Horn and Jackson (1972) | Detailed-balanced mass action has a free-energy Lyapunov function: the sink of HJ5 relaxes. |
| Hudson (1974); Spekkens (2008); T5 | Negativity cannot be removed: a pure state's Wigner function is non-negative only for Gaussians, every quasi-probability representation of quantum theory has negativity, and a normalised orthogonal non-negative phase-space kernel is a permutation. |
| Nedjalkov et al. (2004); Wagner (2016); Muscato and Wagner (2016); Shao and Xiong (2019, 2020) | Rigorous signed-particle and branching representations of the Wigner equation: the Hahn–Jordan split of the kernel into two positive branches, with weak or expectation solutions proved; Shao and Xiong identify the variance growth, set by the positive parts, as the sign problem. |
| Booth, Thom and Alavi (2009); Spencer, Blunt and Foulkes (2012) | The same structure elsewhere — signed walkers, annihilation on coincidence, a critical population below which it fails (R §7.5). |

The answer therefore has three parts.

1. **A realistic process exists.** HJ1 is a pure-jump Markov process with
   non-negative rates on a positive state space — configurations of worlds
   carrying a species label — and $`W`$ is its signed count, exactly. Part C
   runs it in phase space for $`V = V_0\cos x`$, where the whole of $`V`$ is in
   $`J`$ and there is no force at all, and recovers the Schrödinger density to
   sampling noise.
2. **The sign does not go away.** T5, Hudson and Spekkens forbid that. It moves
   from the transition probability, where a negative value is meaningless, to a
   species label, where it is a conserved charge — every event adds one world
   of each sign — exactly as T §5.3 concluded.
3. **The price is now explicit.** HJ2: at least $`\lVert J\rVert_{TV}`$ worlds per
   world per unit time, and HJ2.1: no linear process can take them back. What
   was "unrealistic" in 1954 is, in this reading, a requirement for a nonlinear
   $`W`$-neutral sink — the sea's recombination, Q8's contact sink, or a
   chemistry with $`\varepsilon \gg kT`$.

---

## 8. Numerical verification

`PYTHONPATH=src python3 src/demo_jordan_lifting.py` runs Parts A to F in about
fifteen seconds (container, October 2026). Seeds are fixed in the code.

- **Part A (HJ1, HJ3).** SymPy: $`\partial_t W = (A - B)W`$,
  $`\partial_t S = (A + B)S`$; single-mode moments as in the table of §4.
- **Part B (HJ2).** Three cosine modes at one $`x`$ on a periodic lattice of
  121 rows: $`\lVert J\rVert_{TV} = 2\nu = 2.570206`$. The Jordan lifting
  reproduces $`e^{tL}W_0`$ to $`1.7\times10^{-16}`$, $`3.3\times10^{-16}`$ and
  $`9.1\times10^{-15}`$ at $`t = 0.5, 1, 2`$, and $`S = e^{2\nu t}`$ to every
  printed digit (3.615041, 13.068521, 170.786253). No column rate among 2000
  random positive splits falls below the bound (minimum excess $`0`$). For a
  generic generator with a Markov part added, the Jordan lift's column rates
  equal $`2\nu_n`$ to $`8.9\times10^{-16}`$, and the Markov part alone has
  $`\nu = 0`$.
- **Part C (HJ1 in phase space).** $`V = \cos x`$, a Gaussian packet at
  $`x_0 = \pi`$, $`\sigma_x = 0.5`$, $`p_0 = 0.6`$, $`2\times10^5`$ initial
  worlds, step $`0.005`$, to $`t = 1.5`$; errors from 20 independent lineage
  groups. Worlds grow to $`1.23\times10^6`$: $`\ln(N/N_0) = 1.8175`$ against
  the HJ2 rate integrated along the run, $`1.8228`$. Density against a
  split-operator Schrödinger solution: $`L^1`$ difference $`0.0464`$ against a
  sampling-noise level of $`0.0433`$, rms $`z = 1.15`$ over 113 bins.
  $`\langle p\rangle = 0.2372 \pm 0.0049`$ (exact $`0.2363`$),
  $`\langle p^2\rangle = 0.457 \pm 0.011`$ (exact $`0.463`$). The normalisation
  $`\sum\sigma/N_0 = 1`$ exactly.
- **Part D (HJ5).** The tables of §6.1 and §6.2, the $`\mathrm{P} + \mathrm{N}`$
  inflation (e.g. $`K' = 10^3`$: 0.01998 against $`2n_N/(K'\lvert W\rvert) = 0.02`$;
  balanced 0.6315 against $`2\sqrt{n/K'} = 0.6325`$), and the relaxation.
- **Part E (HJ4).** The moments of Proposition HJ4.
- **Part F.** The figure.

![The Jordan lifting and pair chemistry](https://raw.githubusercontent.com/billpage/wpmw/output/figures/jordan_lifting.png)

*Left: the density of a packet in $`V = \cos x`$ at $`t = 1.5`$, as the signed
count of the positive branching process of HJ1 (band: ±2 standard errors from
20 lineage groups) and from the Schrödinger equation. Centre: the world count
$`N/N_0`$ of the same run and the HJ2 rate $`\langle 2\lvert\Gamma(x)\rvert\rangle`$
integrated along it — the cost of the sign, paid at exactly the minimum.
Right: equilibrium inflation $`N - \lvert W\rvert`$ of
$`\mathrm{PP} + \mathrm{NN} \rightleftharpoons 2\mathrm{PN}`$ against
$`\varepsilon/kT`$, for a cell with a net sign and a balanced one, with the
asymptotes of HJ5(d); $`\varepsilon/kT \to \infty`$ is the contact sink.*

---

## 9. Open items

- **HJ-SP1.** Couple HJ5's chemistry to Cyganski's whole-pair kicks (P–P and
  N–N moving by a whole photon, P–N moving as a neutral packet) and test
  whether rows relax: the depletion functional of V9 is the measure, V10 the
  comparison.
- **HJ-SP2.** HJ2 for the compensated residual: compute
  $`\lVert K_{\rm res}(x,\cdot)\rVert_{TV}`$ as a function of the reach for the
  cosine-plus-harmonic and Eckart problems, and compare with R §7.2's
  uncompensated $`2\lvert\Gamma_q\rvert`$. Does moving the classical part into
  streaming lower the minimal cost, and by how much?
- **HJ-SP3.** HJ2.1 makes every sink bilinear. Measure, at finite
  $`\varepsilon/kT`$, the bias HJ5 leaves in $`W`$ through the cell
  coarse-graining (Q8's decoherence) against the bias of the instantaneous
  contact sink.

---

## 10. Summary

1. Any signed generator $`L = A - B`$ with $`A`$ Metzler and $`B \ge 0`$ lifts to
   a positive two-species process with $`W = u^+ - u^-`$ obeying $`L`$ exactly
   and $`S = u^+ + u^-`$ obeying $`A + B`$ (HJ1). The Jordan split of
   Takabayasi's kernel is the emissive realisation without a sea.
2. The unsigned count grows at least at twice the negative off-diagonal mass,
   $`\lVert J\rVert_{TV}`$ for a kernel, and the Jordan split attains it (HJ2);
   so no positive linear process removes worlds, and every sink is at least
   bilinear (HJ2.1).
3. Pawula's theorem acts on $`S`$, which diffuses; $`W`$ keeps only odd moments
   (HJ3).
4. The lifting fixes the generator, not the process: up/down actions and the
   four-action defocus agree in mean and differ in covariance (HJ4).
5. Pair chemistry with $`\varepsilon = 2D_{PN} - D_{PP} - D_{NN} \gg kT`$ is a
   $`W`$-neutral, relaxing sink whose limit $`\varepsilon \to \infty`$ is the
   contact sink and the minimal ensemble (HJ5). "Deeper than both" is
   sufficient, not necessary. It does not repair depletion across rows and does
   not escape the coincidence problem.
6. Takabayasi's verdict survives as a price, not an impossibility: a positive
   process exists, the sign is relocated rather than removed, and the cost is
   $`\lVert J\rVert_{TV}`$ per world, payable only by a nonlinear sink.

---

## 11. Sources

- T. Takabayasi, Prog. Theor. Phys. **11**, 341 (1954), §3(a), and
  [`takabayasi_1954_stochastic_picture.md`](takabayasi_1954_stochastic_picture.md).
- R. F. Pawula, Phys. Rev. **162**, 186 (1967).
- N. Ikeda, M. Nagasawa and S. Watanabe, "Branching Markov processes I–III",
  J. Math. Kyoto Univ. **8** (1968) 233, 365; **9** (1969) 95.
- F. Horn and R. Jackson, Arch. Rational Mech. Anal. **47**, 81 (1972).
- R. L. Hudson, Rep. Math. Phys. **6**, 249 (1974).
- R. W. Spekkens, Phys. Rev. Lett. **101**, 020401 (2008).
- M. Nedjalkov et al., Phys. Rev. B **70**, 115319 (2004).
- W. Wagner, "A random cloud model for the Wigner equation", Kinet. Relat.
  Models **9**, 217 (2016) (WIAS Preprint 2083, 2015).
- O. Muscato and W. Wagner, SIAM J. Sci. Comput. **38**, A1483 (2016).
- S. Shao and Y. Xiong, SIAM J. Numer. Anal. **58**, 2589 (2020)
  (arXiv:1907.01897); Numer. Math. Theor. Meth. Appl. **12**, 21 (2019).
- G. H. Booth, A. J. W. Thom and A. Alavi, J. Chem. Phys. **131**, 054106
  (2009); J. S. Spencer, N. S. Blunt and W. M. C. Foulkes, J. Chem. Phys.
  **136**, 054110 (2012).
- P. M. Morse, Phys. Rev. **34**, 57 (1929).
- J. H. Sullivan, J. Chem. Phys. **46**, 73 (1967).
- Step 24, Propositions V2, V9, V10:
  [`../analysis/sea_depletion.md`](../analysis/sea_depletion.md) §3, §4.2, §4.3.
- Step 23, Proposition Q8:
  [`../analysis/force_blind_sea.md`](../analysis/force_blind_sea.md) §6.2.
- Lemma CX3: [`world_crossings.md`](world_crossings.md) §4.
- R §7: [`representation_cost_and_annihilation.md`](representation_cost_and_annihilation.md).
