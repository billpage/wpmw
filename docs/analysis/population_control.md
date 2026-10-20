# Population control below the cell

> Population control for the world ensemble of step 25 without the contact
> sink's decoherence (open item WC-SP9). The contact sink keeps only each cell's
> net count, the zeroth moment of the signed worlds in it. **Proposition PC2**:
> a sink that keeps the moments up to order `K` changes a smooth observable by
> at most the cell to the power `K + 1`, so the contact sink's error is first
> order in the cell, as WC10 measured. **Definition PC3**, K1 resampling: in
> each cell keep the majority sign's net worlds and add ± pairs that carry the
> first moment of the removed worlds. **Lemma PC4** places each pair at the
> *dipole centre*, the unique point at which it also matches the removed worlds'
> second moment along the dipole, so a lone pair is returned unchanged and a
> dipole moves with the flow instead of being pinned to its cell. **Lemma PC1**:
> the overlap of two independently seeded ensembles estimates the purity of the
> mean ensemble with no variance term; it agrees with the self estimator of step
> 25 to 0.012. **Proposition PC5**, absorption with a tolerance `delta` (a
> daughter realised by removing a world of the other sign within `delta` of it),
> measured with three particles: the absorbed share stays at 0.31 to 0.44, below
> the one half Theorem N3 requires, so absorption alone does not hold the
> population, and its bias is at least that of a contact sink of cell `2 delta`.
> Matching the event's dipole does not help, because the error left is the
> event's dipole moved by up to `delta`, on the scale of the fringes the events
> write. The ledger's absorptive branch is exact only where a mesh puts partners
> on the targets. **Proposition PC6**, measured on problem G of step 25 at reach
> `2 pi` against the 4D mesh: K1 at cell 0.7, with the grid shifted every ten
> steps, leaves a global-purity deficit of -0.020 at `t = 3` and `nu = 1e6`,
> against -0.143 for the contact sink at cell 0.35, at 1.37 against 1.22 worlds
> per `nu`; the one-body purity is held to -0.012 against -0.082. The deficit is
> -0.022 at `nu = 1e5`, so it is the rule's limit at that cell, and it shrinks
> with the cell to -0.014 at cell 0.5. With too few worlds per cell K1 makes the
> ensemble too pure (PC-SP1).

1. [What this note asks, settles and leaves open](#0-what-this-note-asks-settles-and-leaves-open)
2. [Setting, observables and the yardstick](#1-setting-observables-and-the-yardstick)
3. [Proposition PC2: what a sink keeps decides its error](#2-proposition-pc2-what-a-sink-keeps-decides-its-error)
4. [K1 resampling and the dipole centre](#3-k1-resampling-and-the-dipole-centre)
5. [Proposition PC5: absorption with a tolerance](#4-proposition-pc5-absorption-with-a-tolerance)
6. [Results](#5-results)
7. [What the specification takes from this](#6-what-the-specification-takes-from-this)
8. [Figure](#7-figure)
9. [Numerical verification](#8-numerical-verification)
10. [Open items](#9-open-items)
11. [Sources](#10-sources)

## 0. What this note asks, settles and leaves open

**Asks.** Step 25 ([`world_clocks.md`](world_clocks.md) §8.2) ran three
particles as signed worlds, one bank per bond, and found that the world
ensemble loses purity that the bond-sum generator on a mesh keeps. It traced
the loss to the contact sink (Proposition WC10): the sink keeps each cell's
net count but replaces the signed structure inside a cell by the majority
sign's, which removes interference finer than a cell. Open item WC-SP9 asked
for a population-control rule that conserves $`W`$ below the cell, and named
three candidates: annihilating only near-coincident pairs, moment-preserving
resampling within a cell, and the absorptive ledger. This note tests the last
two against the mesh.

**Inherits.** Problem G of step 25 (three unit-mass particles on a line, a
harmonic trap of unit frequency, the Gaussian pair interaction $`g = 2`$,
$`s = 0.7`$, quenched from the product of trap ground states), the world
algorithm of `demo_three_worlds.py` (one bank of integrate-and-fire
integrators per bond, events realised as emission), its observables, and its
4D Jacobi mesh at reach $`2\pi`$ as the reference. Theorem N3 of step 18
([`stochastic_ledger.md`](stochastic_ledger.md) §4) for the balance between
emission, absorption and any other sink.

**Settles.**

- **Lemma PC1.** The product of two independently seeded ensembles estimates
  the purity of the mean ensemble with no variance term. The self estimator,
  with self-pairs removed, is unbiased only for independent worlds; here the
  two agree to $`0.012`$, so the step 25 numbers stand.
- **Proposition PC2.** A sink that keeps a cell's signed moments up to order
  $`K`$ changes any smooth observable by $`O(c^{K+1})`$ per cell. The contact
  sink keeps only $`M_0`$, so its error is first order in the cell, as WC10
  measured.
- **Definition PC3 and Lemma PC4.** K1 resampling keeps $`M_0`$ and $`M_1`$ of
  every cell with the majority's net worlds plus $`\pm`$ pairs. Each pair sits
  at the *dipole centre*, the unique point at which it also matches the
  removed worlds' second moment along the dipole, so a lone pair is carried
  exactly and a dipole moves with the flow.
- **Proposition PC5.** Absorption with a tolerance $`\delta`$ cannot hold the
  population alone (the absorbed share stays below one half) and is at
  least as biased as a contact sink of cell $`2\delta`$.
- **Proposition PC6.** K1 at cell $`0.7`$, with the grid shifted every ten
  steps, cuts the contact sink's purity deficit sevenfold at about the same
  population, and the deficit does not change between $`\nu = 10^5`$ and
  $`10^6`$. At $`\nu = 10^6`$ the deficit shrinks steadily with the cell
  ($`-0.099`$, $`-0.020`$, $`-0.014`$ at $`c = 1.0`$, $`0.7`$, $`0.5`$).

**Leaves open.** A small early excess of purity under K1, and a larger one
when a cell holds too few worlds (PC-SP1); K1's cost in higher relative
dimension (PC-SP2); a principled shift schedule (PC-SP3); and whether
absorption can be made exact off a mesh (PC-SP4).

**Corrects.** Nothing. It narrows WC-SP9 and replaces item 5 of the
specification outline in step 25 §8.2.

## 1. Setting, observables and the yardstick

**State.** A world is a point $`z = (\mathbf x, \mathbf p) \in \mathbb R^6`$
with a sign $`\varepsilon = \pm 1`$ and three banks of integrators. The
centre of mass separates exactly and is never read, so every cell, distance
and observable below is taken in the Jacobi relative phase space
$`\zeta = (\rho_1, \rho_2, \pi_1, \pi_2) \in \mathbb R^4`$, with $`\zeta = J z`$ for
the orthonormal Jacobi map $`J`$. The ensemble of $`N`$ worlds represents

```math
W(\zeta) \approx \frac{1}{\nu}\sum_{k=1}^{N}\varepsilon_k\,\delta(\zeta - \zeta_k),
```

where $`\nu`$ is the number of worlds at $`t = 0`$, all positive.

**Observables.** The global purity at resolution $`a`$,
$`P_a = (2\pi\hbar)^2\int (W * G_a)^2\,d\zeta`$, with $`G_a`$ a Gaussian of
width $`a = 0.25`$ in every Jacobi coordinate, and the one-body purity
$`P_1 = 2\pi\hbar\int W_1^2`$ of particle 1, with its centre of mass
convolved exactly (step 25 §8.2). A unitary evolution keeps
$`P_0 = 1`$; $`P_a`$ falls as $`W`$ develops structure finer than $`a`$, so it
is the test of whether a realisation keeps that structure. Every result is
quoted as *worlds minus mesh* at the same reach and time, which separates
the realisation's error from the generator's.

**Lemma PC1 (cross-seed estimator).** *Let $`\hat W_1`$ and $`\hat W_2`$ be
two ensembles from independent seeds, each with mean $`\bar W`$. Then*

```math
\mathbb E\!\left[(2\pi\hbar)^2\int (\hat W_1 * G_a)(\hat W_2 * G_a)\,d\zeta\right]
 = (2\pi\hbar)^2\int (\bar W * G_a)^2\,d\zeta .
```

*The self estimator $`(2\pi\hbar)^2\int(\hat W * G_a)^2`$ minus its
self-pair terms has expectation $`P_a(\bar W)`$ plus the integrated
covariance of distinct worlds, which vanishes only if the worlds are
independent.*

*Proof.* Independence makes the expectation of the product the product of
the expectations. Expanding the self estimator over pairs $`k \ne l`$ leaves
$`\sum_{k \ne l}\mathrm{Cov}`$ terms; daughters of one event, and
worlds placed together by a sink, are correlated. $`\square`$

The two daughters of an event are anti-correlated, and the pairs a sink
places are correlated, so the self estimator is biased by an amount that
grows with $`N/\nu`$. Measured, the two estimators differ by at most
$`0.012`$ in $`P_a`$ (K1 at the cell centre, cell $`0.5`$, $`N/\nu = 6.7`$) and
by at most $`0.005`$ at $`N/\nu \le 1.4`$ (§5). The step 25 figures, made with the self
estimator, stand. Every number below is the cross estimator unless stated.

**Why not a crossing probability.** Between events every world follows the
classical Hamiltonian flow, and two trajectories of a smooth flow in phase
space never meet. The flow only stretches and shears a separation, by
about $`e^{\lambda t}`$. Two worlds whose positions cross have different
momenta and do not cancel in $`W`$. The only way two worlds come together is
an event's momentum kick, which is what absorption uses (§4). So a rule
must work with a *distance*, and the question is what a removal at a
distance $`d`$ costs.

**The yardstick.** Removing a positon at $`\zeta_+`$ and a negaton at
$`\zeta_-`$ changes the smoothed state by

```math
\Delta(W * G_a)(\zeta) = -\frac{1}{\nu}\bigl[G_a(\zeta - \zeta_+) - G_a(\zeta - \zeta_-)\bigr]
 \approx -\frac{1}{\nu}\,\mathbf d\cdot\nabla G_a(\zeta - \bar\zeta),\qquad \mathbf d = \zeta_+ - \zeta_- ,
```

a dipole of moment $`\mathbf d`$. Across a fringe the positons lie on one side
and the negatons on the other, so the dipoles a sink removes are aligned
with the fringe and the removal flattens it.

## 2. Proposition PC2: what a sink keeps decides its error

Let a cell $`\mathcal C`$ with centre $`\zeta_c`$ hold the signed measure
$`\mu = \sum_k \varepsilon_k\delta_{\zeta_k}`$, and let a sink replace it by
$`\mu' = \sum_l \varepsilon'_l\delta_{\zeta'_l}`$. The signed moments about the
centre are

```math
M_\alpha(\mu) = \sum_k \varepsilon_k\,(\zeta_k - \zeta_c)^\alpha,
```

for a multi-index $`\alpha \in \mathbb N^4`$, with $`M_0`$ the net count,
$`M_1`$ the first-moment vector and $`M_2`$ the second-moment matrix.

**Proposition PC2.** *If $`M_\alpha(\mu') = M_\alpha(\mu)`$ for every
$`\lvert\alpha\rvert \le K`$, and every old and new point lies within
Chebyshev distance $`r`$ of $`\zeta_c`$, then for every
$`\varphi \in C^{K+1}`$*

```math
\Bigl\lvert\int\varphi\,d\mu' - \int\varphi\,d\mu\Bigr\rvert
 \le \frac{(4r)^{K+1}}{(K+1)!}\,\max_{\lvert\beta\rvert = K+1}\lVert\partial^\beta\varphi\rVert_\infty\,(n + n'),
```

*where $`n`$ and $`n'`$ count the worlds before and after.*

*Proof.* Expand $`\varphi`$ about $`\zeta_c`$ to order $`K`$ with Lagrange
remainder. The polynomial part integrates to $`\sum_{\lvert\alpha\rvert\le K}
\partial^\alpha\varphi(\zeta_c)M_\alpha/\alpha!`$, the same for both measures.
At each point the remainder is at most
$`\max\lVert\partial^\beta\varphi\rVert\sum_{\lvert\beta\rvert = K+1}
\lvert h^\beta\rvert/\beta! = \max\lVert\partial^\beta\varphi\rVert\,
(\sum_i\lvert h_i\rvert)^{K+1}/(K+1)!`$ with $`\sum_i\lvert h_i\rvert \le 4r`$,
and there are $`n + n'`$ points. $`\square`$

For $`P_a`$ the test function is $`G_a`$, whose derivatives of order
$`K+1`$ scale as $`a^{-(K+1)}`$, so the bound per cell is of order
$`(c/a)^{K+1}`$ times the worlds the cell touches. The contact sink keeps
only $`M_0`$ ($`K = 0`$): its error is first order in the cell. That is the
roughly linear deficit WC10 measured. A rule with $`K = 1`$ is second order.
The bound is a statement as $`c/a \to 0`$. At the cells used here
$`c/a`$ is $`1.4`$ to $`4`$, so the gain has to be measured (§5).

## 3. K1 resampling and the dipole centre

**Definition PC3 (K1 resampling).** *In every cell holding both signs, keep
$`\lvert m\rvert`$ worlds of the majority sign, chosen at random, where
$`m = n_+ - n_-`$, and remove the rest. The removed set is balanced, so its
net count is zero; let*

```math
\mathbf R = \sum_{k\ \text{removed}}\varepsilon_k\,\zeta_k,\qquad
A = \sum_{k\ \text{removed}}\varepsilon_k\,\zeta_k\zeta_k^{\mathsf T}
```

*be its first moment (a vector) and second moment (a matrix). If
$`\mathbf R \ne 0`$, add $`k = \lceil\lVert\mathbf R\rVert_\infty/c\rceil`$ pairs,
each a positon at $`\mathbf c' + \mathbf R/2k`$ and a negaton at
$`\mathbf c' - \mathbf R/2k`$, with $`\mathbf c'`$ the dipole centre of Lemma
PC4, and with fresh integrator phases. Then $`M_0`$ and $`M_1`$ of every cell
are unchanged, every pair's separation is at most $`c`$, and the cell keeps
$`\lvert m\rvert + 2k`$ worlds.*

The centre of mass of a new pair is the removed worlds' mean, and its
relative part is set by Lemma PC4. With $`k > 1`$ the pairs coincide. Their
second moment is the same as one pair's, and the next application merges
them.

**Lemma PC4 (dipole centre).** *Given $`\mathbf R \ne 0`$ and $`A`$, a pair of
moment $`\mathbf R`$ centred at $`\mathbf c'`$ has second moment
$`\mathbf c'\mathbf R^{\mathsf T} + \mathbf R\mathbf c'^{\mathsf T}`$. There is
exactly one $`\mathbf c'`$ for which this matches $`A`$ along the dipole,
$`(\mathbf c'\mathbf R^{\mathsf T} + \mathbf R\mathbf c'^{\mathsf T})\mathbf R = A\mathbf R`$, namely*

```math
\mathbf c'\cdot\mathbf R = \frac{\mathbf R^{\mathsf T}A\,\mathbf R}{2\lvert\mathbf R\rvert^2},\qquad
\mathbf c' = \frac{A\mathbf R - \mathbf R\,(\mathbf c'\cdot\mathbf R)}{\lvert\mathbf R\rvert^2}.
```

*If the removed set is itself one pair, or $`k`$ coincident pairs, the
rule returns it unchanged.*

*Proof.* The second moment of a positon at $`\mathbf c' + \mathbf d/2`$ and a
negaton at $`\mathbf c' - \mathbf d/2`$ is $`\mathbf c'\mathbf d^{\mathsf T} +
\mathbf d\mathbf c'^{\mathsf T}`$, and $`k`$ such pairs with $`\mathbf d =
\mathbf R/k`$ sum to the same expression in $`\mathbf R`$. Applying it to
$`\mathbf R`$ gives $`\mathbf c'\lvert\mathbf R\rvert^2 + \mathbf R(\mathbf c'\cdot\mathbf R) = A\mathbf R`$.
The dot product with $`\mathbf R`$ gives $`2\lvert\mathbf R\rvert^2(\mathbf c'\cdot\mathbf R)
= \mathbf R^{\mathsf T}A\mathbf R`$, and substituting back gives $`\mathbf c'`$. For one pair
with centre $`\mathbf c_0`$ and moment $`\mathbf R`$,
$`A = \mathbf c_0\mathbf R^{\mathsf T} + \mathbf R\mathbf c_0^{\mathsf T}`$, so
$`\mathbf c' = \mathbf c_0`$. $`\square`$

**Why the centre matters.** Where a dipole sits is not a first-moment
quantity: $`\mathbf R`$ is the same wherever a pair lies in the cell. A rule
that places the pair at the cell's centre therefore discards the pair's
motion at every application. The dipole is pinned, fringes stop streaming
and dephasing, and the ensemble ends too pure (§5, the variant *K1 at the
cell centre*). Lemma PC4 keeps the component of $`A`$ that records where the
dipole is, so a pair that moved within a step is put back where it moved
to. It is exact for a lone pair, the case that recurs at every step.

**Cells and shifts.** On a fixed grid a dipole centre can leave its cell, and
the two members of a pair straddling a face would then sit in two cells as
two unbalanced worlds, never merged. Clipping the centre into the cell
keeps the pair together but holds a dipole at a face. Shifting the grid by
a random offset at *every* application lets everything through but regroups
every cell at every step, so far more pairs are annihilated (the same
change makes the plain contact sink worse, §5). The compromise used below
clips the centre into the cell and shifts the grid every ten applications;
it is labelled K1s10.

## 4. Proposition PC5: absorption with a tolerance

On a mesh an event's two daughters land on rows, and the absorptive
realisation removes a negaton at the parent's row $`+\xi`$ and a positon at
the parent's row $`-\xi`$ exactly (step 18). Off a mesh, positions and
momenta are continuous and nothing sits exactly at a target, so absorption
needs a tolerance.

**Rule.** A daughter of sign $`s`$ due at $`\zeta_t`$ is realised by removing
an existing world of sign $`-s`$ within Chebyshev distance $`\delta`$ of
$`\zeta_t`$, if there is one, and otherwise by creating it. In the *paired*
form, both daughters of an event must find partners or both are created (the
ledger's absorptive branch, $`N \to N - 2`$). In the *single* form, each
half is realised alone. In the *dipole-matched* form, the pair of partners
is chosen to make their displacements from the two targets agree.

**Proposition PC5.** *(a) If the partners lie at $`\zeta_+ + \mathbf a`$
and $`\zeta_- + \mathbf b`$, the realised change differs from the event's by
zero net count and first moment $`\tau(\mathbf a - \mathbf b)`$, $`\tau`$ the
event's sign. With $`\mathbf a = \mathbf b`$ it is the event's dipole translated by
$`\mathbf a`$. (b) (Measured, problem G, $`\nu = 10^5`$.) The absorbed share
$`f`$ of daughters is $`0.31`$ at $`\delta = 0.1`$ (by $`t = 0.75`$), $`0.40`$ at
$`0.175`$ and $`0.44`$ at $`0.25`$. It stays below one half, so by Theorem N3
($`\tfrac12 - f = R_{\rm sink}/2\Gamma_{\rm tot}`$, step 18) the population
grows unless another sink removes the difference. It passes $`12\nu`$ before
$`t = 1`$ at $`\delta = 0.1`$ and before $`t = 2.75`$ at $`\delta = 0.175`$.
(c) (Measured.) The bias is at least that of a contact sink of cell
$`2\delta`$: the single form at $`\delta = 0.175`$ gives a $`P_a`$ deficit of
$`-0.152`$ at $`t = 3`$, against $`-0.150`$ for the contact sink at
$`c = 0.35`$, with $`3.3`$ times the worlds; the paired form at
$`\delta = 0.25`$ gives $`-0.239`$, against $`-0.187`$ for the contact sink at
$`c = 0.5`$ (step 25, self estimator). Matching the dipole barely changes it
($`-0.235`$ against $`-0.239`$ at $`\delta = 0.25`$). (d) (Measured.) With K1 running
alongside, only $`0.6\,\%`$ of daughters find a partner within
$`\delta = 0.1`$, and the result is K1's.*

*Proof of (a).* Removing a world of sign $`-\tau`$ at $`\zeta_+ + \mathbf a`$
adds $`+\tau\delta_{\zeta_+ + \mathbf a}`$, and removing one of sign $`+\tau`$
at $`\zeta_- + \mathbf b`$ adds $`-\tau\delta_{\zeta_- + \mathbf b}`$. Subtract
the event's change $`\tau(\delta_{\zeta_+} - \delta_{\zeta_-})`$ and take
moments. $`\square`$

Why dipole matching does not help: the error left in (a) is the event's
dipole moved by up to $`\delta`$, and the event's dipole *is* the fine
structure. The channels are $`\xi_q = q\,\Delta p`$ with $`\Delta p = 0.25`$
at reach $`2\pi`$, so the events write fringes on the scale $`\Delta p`$, and
$`\delta = 0.175`$ to $`0.25`$ moves them by a sizeable fraction of their own
period. Absorption is accurate only for $`\delta \ll \Delta p`$, and there
the partners are too few to hold the population. The ledger's absorptive
branch is exact because a mesh puts the partners on the targets. Without a
mesh it becomes one more annihilation at scale $`\delta`$, applied at the
event rate.

## 5. Results

All runs: problem G, reach $`2\pi`$, $`\Delta t = 0.01`$, $`t \le 3`$, two
seeds per rule, sink applied every step. The mesh's own values are
$`P_a = 0.6653`$, $`0.5793`$, $`0.5597`$ and $`P_1 = 0.8441`$, $`0.7454`$, $`0.7618`$
at $`t = 1`$, $`2`$, $`3`$.

**One application** (Part P). $`10^6`$ worlds sampled from the mesh at
$`t = 2.5`$, one sink, no dynamics: the change in $`P_a`$ and the fraction of
worlds kept.

| cell | contact: $`\Delta P_a`$, kept | K1 at cell centre | K1 at dipole centre |
|---|---|---|---|
| 0.175 | −0.0013, 0.963 | 0.0000, 0.998 | 0.0000, 0.997 |
| 0.25 | −0.0044, 0.900 | −0.0002, 0.978 | −0.0002, 0.982 |
| 0.35 | −0.0117, 0.828 | −0.0009, 0.937 | −0.0006, 0.947 |
| 0.5 | −0.0158, 0.767 | −0.0017, 0.886 | −0.0014, 0.879 |
| 0.7 | −0.0358, 0.670 | −0.0092, 0.787 | −0.0032, 0.807 |

$`P_1`$ moves by at most $`2\times10^{-4}`$ under K1 at the dipole centre and
$`8\times10^{-4}`$ at the cell centre, and by up to $`0.013`$ under the
contact sink.

**Runs** (Parts R and C). Worlds minus mesh at $`t = 3`$, cross estimator,
and $`N/\nu`$:

| rule | cell | $`\nu`$ | $`P_a`$ | $`P_1`$ | $`N/\nu`$ |
|---|---|---|---|---|---|
| contact | 0.7 | $`10^5`$ | −0.211 | −0.139 | 1.12 |
| contact | 0.35 | $`10^5`$ | −0.150 | −0.088 | 2.27 |
| contact | 0.35 | $`10^6`$ | −0.143 | −0.082 | 1.22 |
| contact, grid shifted every step | 0.7 | $`10^5`$ | −0.280 | −0.189 | 1.06 |
| absorption only, paired, $`\delta = 0.25`$ | — | $`10^5`$ | −0.239 | −0.147 | 4.90 |
| absorption only, dipole-matched, $`\delta = 0.25`$ | — | $`10^5`$ | −0.235 | −0.138 | 4.93 |
| absorption only, single, $`\delta = 0.175`$ | — | $`10^5`$ | −0.152 | −0.079 | 7.48 |
| K1 at cell centre | 0.7 | $`10^5`$ | −0.046 | +0.005 | 5.01 |
| K1, dipole centre, shift every step | 0.7 | $`10^5`$ | −0.064 | −0.028 | 2.51 |
| K1, dipole centre, fixed grid | 0.7 | $`10^5`$ | −0.041 | −0.011 | 3.10 |
| K1, dipole centre, shift every 5 | 0.7 | $`10^5`$ | −0.046 | −0.019 | 2.85 |
| **K1s10** | 0.7 | $`10^5`$ | **−0.022** | **+0.000** | 2.94 |
| K1, dipole centre, shift every 20 | 0.7 | $`10^5`$ | −0.021 | 0.000 | 2.99 |
| K1s10, centre not clipped | 0.7 | $`10^5`$ | −0.032 | −0.009 | 2.95 |
| K1s10 + absorption, $`\delta = 0.1`$ | 0.7 | $`10^5`$ | −0.026 | −0.002 | 2.90 |
| **K1s10** | 0.7 | $`10^6`$ | **−0.020** | **−0.012** | 1.37 |
| K1s10 | 1.0 | $`10^5`$ | −0.087 | −0.047 | 1.54 |
| K1s10 | 1.0 | $`10^6`$ | −0.099 | −0.057 | 1.12 |
| K1s10 | 0.5 | $`10^5`$ | +0.021 | +0.020 | 7.75 |
| **K1s10** | 0.5 | $`10^6`$ | **−0.014** | **−0.002** | 2.04 |
| K1 at dipole centre, shift every step | 0.5 | $`10^5`$ | −0.006 | +0.011 | 6.60 |

At $`t = 1`$ every K1 rule at a cell of $`0.7`$ or less is within $`0.035`$ of the
mesh, and the contact sinks are $`0.04`$ to $`0.11`$ below it. The full time
series are in `popctl_C_nu*.json`.

**Proposition PC6 (measured).** *For problem G at reach $`2\pi`$, K1
resampling at cell $`0.7`$ with the grid shifted every ten applications
(K1s10) holds the global purity at resolution $`0.25`$ to $`-0.020`$ of the
mesh at $`t = 3`$ and $`\nu = 10^6`$, against $`-0.143`$ for the contact sink at
cell $`0.35`$, at $`N/\nu = 1.37`$ against $`1.22`$. The one-body purity is held
to $`-0.012`$ against $`-0.082`$. The deficit is $`-0.022`$ at $`\nu = 10^5`$,
so at this cell it is the rule's limit, not sampling error. At
$`\nu = 10^6`$ the deficit shrinks steadily with the cell: $`-0.099`$, $`-0.020`$,
$`-0.014`$ at $`c = 1.0`$, $`0.7`$, $`0.5`$, with $`N/\nu = 1.12`$, $`1.37`$,
$`2.04`$. At equal population the gain is smaller at the coarsest cell:
$`-0.099`$ against the contact sink's $`-0.211`$ at $`c = 0.7`$, both at
$`N/\nu = 1.12`$. K1 adds $`2k`$ worlds to
each cell holding both signs, with $`k = 1`$ unless the cell's dipole is
longer than the cell, so its cost follows the number of such cells, not
$`\nu`$: $`N - \nu`$ at $`t = 3`$ is $`1.9\times10^5`$ at $`\nu = 10^5`$ and
$`3.7\times10^5`$ at $`\nu = 10^6`$.*

Four more things the table shows:

- *The pair's position matters as much as the moments kept.* K1 at the cell
  centre keeps the same two moments as K1 at the dipole centre, yet ends too
  pure at fine cells (+0.033 at $`c = 0.35`$, $`t = 1`$), with its population at
  the cap of $`12\nu`$ by $`t = 1.75`$.
- *Shifting the grid has a cost.* A shift at every application makes the
  contact sink worse, −0.280 against −0.211, and K1 worse, −0.064 against
  −0.041 on a fixed grid. A shift every ten or twenty applications is best,
  −0.022 and −0.021; every five is not, −0.046.
- *Clipping helps.* Leaving the dipole centre outside its cell gives −0.032
  against −0.022.
- *Too few worlds per cell make K1 too pure.* At $`\nu = 10^5`$, K1s10 at
  $`c = 0.5`$ ends above the mesh (+0.021, $`N/\nu = 7.7`$); at $`\nu = 10^6`$
  the same rule ends at −0.014 with $`N/\nu = 2.0`$. The likely reason, not
  yet tested, is that a sparsely filled cell's first moment is mostly
  sampling noise, and a pair placed to carry it adds structure that is not
  in $`W`$. A smaller early excess is part of the rule: at $`c = 0.7`$, $`P_a`$
  is $`+0.010`$ above the mesh at $`t = 1`$ at both $`\nu`$.

## 6. What the specification takes from this

Item 5 of step 25's outline (§8.2, *population control*) becomes:

> *Population control.* K1 resampling (Definition PC3) on cells of the
> relative phase space, with each pair at its dipole centre (Lemma PC4),
> the centre clipped into its cell, and the grid shifted by a random offset
> every ten applications. It keeps $`M_0`$ and $`M_1`$ of every cell exactly
> and, measured, leaves $`2\,\%`$ of $`P_{0.25}`$ at cell $`0.7`$ (PC6). Each
> run reports $`P_a`$ against a reference, or, without one, the cross-seed
> $`P_a`$ of two seeds as a function of the cell. The population must be
> large enough that a cell's first moment is not mostly noise (PC6): at cell
> $`0.5`$, $`\nu = 10^5`$ is too few and $`10^6`$ enough. Absorption with a
> tolerance is not used (PC5).

Two cautions go with it. The bound of PC2 is asymptotic and the cells used
are larger than the resolution, so the rule's error is measured, not
proven, and it turns to an excess when the cells are too sparsely filled
(PC-SP1). And the extra worlds
scale with the number of cells holding both signs, which in six or twelve
relative dimensions is the cost Corollary M10 warned of (PC-SP2).

## 7. Figure

![Population control below the cell](https://raw.githubusercontent.com/billpage/wpmw/output/figures/population_control.png)

*Problem G, reach $`2\pi`$, against the 4D mesh. (a) $`P_{0.25}`$ of the worlds
minus the mesh's against time, self estimator, mean of two seeds, for the
contact sink at two cells, absorption alone, and K1s10 at three cells. (b)
Worlds per $`\nu`$ for the same runs. (c) Every rule at $`t = 3`$: the
cross-seed deficit against $`N/\nu`$. Circles are contact sinks, squares K1
variants, triangles absorption; open at $`\nu = 10^5`$, filled at $`10^6`$.*

## 8. Numerical verification

All outputs are from `src/demo_population_control.py`, which imports the
world algorithm, observables and mesh of `src/demo_three_worlds.py`. Seeds
are fixed in the code. The mesh reference is `three_worlds_M2.json` from
`demo_three_worlds.py M --mesh-reach 2`.

- **Part T.** K1 at the cell centre and at the dipole centre keep $`M_0`$
  exactly and $`M_1`$ of every cell to $`4.8\times10^{-13}`$, over
  $`2\times10^5`$ synthetic signed worlds at cells $`0.35`$ and $`0.7`$, where
  the contact sink changes a cell's $`M_1`$ by up to $`214`$. Every pair's
  separation is at most the cell. Lemma PC4 returns a lone pair's centre to
  $`5.8\times10^{-15}`$ (200 random pairs). Absorption keeps the net weight
  exactly, and its first-moment error per step stays within
  (absorbed halves) $`\times\,\delta`$.
- **Part P.** The one-application table of §5.
- **Part R.** The runs, $`\nu = 10^5`$, in the container (October 2026,
  working tree at commit `60200fb`). The $`\nu = 10^6`$ runs and the K1s10 scan
  over cell and shift period ran on Kaggle as private batch kernels with four
  workers, at commit `60200fb` with this demo applied as a patch:
  `wpmw-popctl-nu1e6` (`R --config K1s10:0.7 contact:0.35 --nu 1000000 --seeds
  2`), `wpmw-popctl-cells` (`R --config K1s10:0.5 K1s10:1.0 K1s5:0.7 K1s20:0.7
  --nu 100000 --seeds 2`) and `wpmw-popctl-cells-nu1e6` (`R --config K1s10:0.5
  K1s10:1.0 --nu 1000000 --seeds 2 --cap 8`).
- **Part C.** The cross-seed and self estimators from snapshots at
  $`t = 1, 2, 3`$, for every rule.
- **Part F.** The figure.

The run JSONs and snapshots are not committed. Each run is seeded, so
rerunning a part with the arguments above reproduces it.

## 9. Open items

- **PC-SP1.** K1 makes the ensemble slightly too pure early on ($`+0.010`$ in
  $`P_a`$ at $`t = 1`$, cell $`0.7`$, at both $`\nu`$), and markedly too pure
  when its cells hold few worlds (cell $`0.5`$ at $`\nu = 10^5`$). One
  application does not show the excess (§5, Part P), so it comes from
  repeated application. Is it the second moment that K1 does not keep
  across the dipole (a K2 rule would test it), the integrator phases reset
  on every new pair, or noise in $`\mathbf R`$ carried as signal? What is the
  smallest number of worlds per cell at which K1 is safe?
- **PC-SP2.** K1 costs up to two worlds per cell holding both signs. With
  four particles on a line (six relative dimensions) or three in space
  (twelve), does that stay affordable? Should cells be on bond coordinates
  (WC-SP10)?
- **PC-SP3.** A principled shift schedule: shift when a dipole would cross a
  face, at a period set by the phase-space speed and the cell, instead of
  every ten steps.
- **PC-SP4.** Absorption off a mesh needs $`\delta \ll \Delta p`$ and then
  finds too few partners. Can events be kept on a momentum lattice, so that
  daughters and partners coincide as on the mesh? *See
  [`../supplement/world_crossings.md`](../supplement/world_crossings.md):*
  worlds do not cross by chance for three or more particles (Lemma CX3),
  but worlds that shared an event's vertex do meet, by construction
  (Proposition CX5).

## 10. Sources

- Vranic, M.; Grismayer, T.; Martins, J. L.; Fonseca, R. A.; Silva, L. O. —
  "Particle merging algorithm for PIC codes." *Comput. Phys. Commun.*
  **191**, 65 (2015). Merging that conserves charge, momentum and energy,
  for positive weights. https://doi.org/10.1016/j.cpc.2015.01.020
- Faghihi, D.; Carey, V.; Michoski, C.; Hager, R.; Janhunen, S.; Chang,
  C. S.; Moser, R. D. — "Moment preserving constrained resampling with
  applications to particle-in-cell methods." *J. Comput. Phys.* **409**,
  109317 (2020). https://doi.org/10.1016/j.jcp.2020.109317
- Xiong, Y.; Shao, S. — "The Wigner branching random walk: efficient
  implementation and performance evaluation." *Commun. Comput. Phys.* **25**,
  871 (2019). Annihilation on a uniform histogram, the $`K = 0`$ case of
  Proposition PC2.
- Step 18, Theorem N3: [`stochastic_ledger.md`](stochastic_ledger.md) §4.
- Step 25, Proposition WC10 and the outline of an $`N`$-particle
  specification: [`world_clocks.md`](world_clocks.md) §8.2.
