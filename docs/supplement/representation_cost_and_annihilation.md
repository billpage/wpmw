# Representation Cost of a World-Particle Ensemble, and the Annihilation Burden

**How many world-particles does it take to represent a quantum state, which of
the project's three representations is cheapest, and what does that tell us
about garbage collection and about the sparse-ontology objection.**

---

## 0. Status and provenance

[`four_action_foundations.md`](four_action_foundations.md) §2.1 separated the
world-particle ensemble from the many-interacting-worlds (MIW) models
*structurally*: MIW discretises the streamlines and has to rebuild momentum
from neighbour spacings, whereas a WPMW ensemble samples $W(x,p)$ over an area
and carries momentum as a coordinate. That settles the kinematics. It does not
say how *expensive* the phase-space ensemble is, and expense is what the
sparse-ontology literature is actually about.

This supplement supplies the quantitative half. It was prompted by Bill Page's
question of 15 August 2026 — how dense must the ensemble be, how do interference
and entanglement scale, is the position-pair $(x,\mu)$ ensemble better or worse,
and how well do positons and negaton holes represent the negative regions — and
by A. Hackebill and B. Poirier, *On Hydrodynamic Formulations of Quantum
Mechanics and the Problem of Sparse Ontology*
([arXiv:2602.21106](https://arxiv.org/abs/2602.21106)), which argues that
branching under decoherence repeatedly partitions a discrete hydrodynamic
ensemble until it is too sparse to sustain quantum dynamics.

Companion code: `src/demo_representation_cost.py`. Every number quoted below is
an output of that script; §10 lists the parts.

**Summary of findings.** The cost of any world-particle representation is one
number, the $L^1$ norm of the object sampled (§1). Measured that way the Wigner
representation is remarkably cheap and the project's other two representations
are not. A cat state costs a factor $2.68$ in samples *regardless of how large
the cat is* (§2). Interference consumes momentum resolution, not position
density, and the spec's lattice $`\Delta p = \pi\hbar/L`$ resolves every fringe
of every state of period $L$ — exactly, and not one bit more (§3). Entanglement
by itself is free; non-Gaussianity is what costs (§4). The pair $(x,\mu)$
ensemble is the Wigner ensemble with the $Y$-integral left to Monte Carlo, and
pays for it by a factor that diverges as the position lattice is refined (§5).
The positon/negaton sea is exact but has the worst readout of the three, by a
factor equal to the number of momentum cells (§6). §7 measures the annihilation
burden for the first time in this project and finds two density criteria that
disagree; §8 proposes what to do about it. §9 answers the sparse-ontology
argument: it does not transfer, but a different threshold does.

---

## 1. The cost functional

Every representation in this project has the same form: a signed or complex
measure $\mu$ on a state space, normalised so that $\int d\mu = 1$, sampled by
$N$ world-particles carrying a weight $w$ of unit modulus. Draw the particles
from $`|\mu| / \|\mu\|_1`$ and give each the weight $w = d\mu/|d\mu|$. Then for
any bounded observable $A$,

```math
\widehat{\langle A\rangle} \;=\; \frac{\|\mu\|_1}{N}\sum_{i=1}^{N} w_i\, A(z_i),
\qquad
\mathrm{std}\!\left(\widehat{\langle A\rangle}\right) \;\simeq\;
\frac{\|\mu\|_1\,\sqrt{\mathrm{Var}_{|\mu|}(A)}}{\sqrt{N}} .
```

**Definition R1 (effective sample size).**
$`N_{\rm eff} = N / \|\mu\|_1^{2}`$.

This is the only figure of merit needed. It applies unchanged to the signed
Wigner ensemble ($`\|\mu\|_1 = \int|W|`$), to the complex position-pair ensemble
($`\|\mu\|_1 = \sum_{ij}|\rho_{ij}|`$), and to the sea-shifted ensemble
($`\|\mu\|_1 = \int (W + 2/h)`$). The three representations differ *only* in
this number, and the differences are large.

Two facts about $`\|\mu\|_1`$ that will do most of the work:

- it is **not** a property of the state alone unless the representation is
  canonical — $`\int|W|`$ is symplectically invariant, $`\sum_{ij}|\rho_{ij}|`$
  is basis-privileged and lattice-dependent;
- it is **multiplicative** under tensor product, so whatever it costs for one
  degree of freedom it costs to the $n$-th power for $n$.

---

## 2. The cat state: the cost saturates

Take $\psi \propto g(x - d/2) + g(x + d/2)$ with
$g(x) = e^{-x^2/4\sigma^2}$. The Wigner function is exactly

```math
W(x,p) \;=\; \frac{e^{-2\sigma^2p^2/\hbar^2}}
{2\pi\hbar\left(1 + e^{-d^2/8\sigma^2}\right)}
\left[
e^{-\frac{(x-d/2)^2}{2\sigma^2}} + e^{-\frac{(x+d/2)^2}{2\sigma^2}}
+ 2\, e^{-\frac{x^2}{2\sigma^2}}\cos\!\left(\tfrac{d\,p}{\hbar}\right)
\right].
```

(Part A checks this closed form against an independent FFT evaluation of the
defining integral; the maximum discrepancy is $`1.7\times10^{-16}`$ against a
peak $|W|$ of $0.318$.)

**Proposition R2.** As $d/\sigma \to \infty$,
$`\|W\|_1 \to 1 + 2/\pi = 1.63662`$, so $`N/N_{\rm eff} \to 2.679`$.

*Proof.* For large separation the three terms have disjoint support in $x$. The
two lobes are non-negative and carry $\tfrac12$ each. The cross term is a
non-negative envelope times $`\cos(dp/\hbar)`$, whose $p$-oscillation is fast
compared with the envelope width $`\hbar/2\sigma`$, so its $L^1$ mass is the
envelope mass times $`\langle|\cos|\rangle = 2/\pi`$. The cross term's envelope
integrates to $1$, giving $`1 + 2/\pi`$. $\square$

Part A:

| $d/\sigma$ | $`\|W\|_1`$ | $\nu = \int W_-$ | $`N/N_{\rm eff}`$ | $\rho(0)$ | fringes |
| --- | --- | --- | --- | --- | --- |
| 2 | 1.0027 | 0.0014 | 1.005 | 3.01e-01 | 0.32 |
| 4 | 1.2082 | 0.1041 | 1.460 | 9.51e-02 | 0.64 |
| 8 | 1.5875 | 0.2937 | 2.520 | 2.68e-04 | 1.27 |
| 16 | 1.6366 | 0.3183 | 2.678 | 9.99e-15 | 2.55 |
| 24 | 1.6366 | 0.3183 | 2.678 | 9.36e-17 | 3.82 |

This is a stronger statement than it looks. **Interference is a bounded,
one-off penalty.** Doubling the cat does not double the cost; past
$d \approx 12\sigma$ it does not change the cost at all. Whatever difficulty a
world-particle ontology has with macroscopic superposition, statistical
resolution of the Wigner function is not it.

**Corollary R2.1 (the displaced-ensemble fact).** The interference lobe carries

```math
\frac{2/\pi}{1 + 2/\pi} \;=\; \frac{2}{\pi + 2} \;=\; 38.90\%
```

of $`\|W\|_1`$, and it sits at $x = 0$, where by the table
$\rho(0) < 10^{-14}$. Nearly two fifths of the ensemble must be placed, in
exactly cancelling positon/negaton pairs, at a location where the particle is
never found.

![Representation cost of a cat state](https://raw.githubusercontent.com/billpage/wpmw/output/figures/representation_cost_cat.png)

*Left: $`\|W\|_1`$ against separation, with the $`1 + 2/\pi`$ asymptote.
Centre: $`W(x,p)`$ at $`d = 8\sigma`$ — the fringes run along $`p`$, not along
$`x`$. Right: the column mass $`\int|W(x,p)|\,dp`$ against the probability
density $`\rho(x)`$. The two peak in different places, and the tallest column of
world-particles stands over the deepest minimum of the density.*

Corollary R2.1 is where this project and MIW part company most sharply, and it
is worth being blunt about the direction: MIW distributes worlds by
$\rho(x)$ and therefore has **no** worlds in the gap, which is where 38.9% of
the quantum content lives. A phase-space ensemble puts them there by
construction. The apparent ontological embarrassment — world-particles where
the particle is not — is the very resource the hydrodynamic models lack.

---

## 3. Interference consumes momentum resolution, not position density

**Proposition R3.** For the cat above, the interference term oscillates in $p$
with period $`2\pi\hbar/d`$ under an envelope of width $`\hbar/2\sigma`$, so the
number of resolvable fringes is $`d/2\pi\sigma`$. It has no oscillatory
structure in $x$ at all: its $x$-dependence is the single Gaussian
$`e^{-x^2/2\sigma^2}`$.

Part B measures this by bin-averaging the exact $W$ (no sampling noise) and
reporting the peak-to-trough of the central column relative to the exact one:

| $\Delta x$ | $\Delta p$ | fringe visibility |
| --- | --- | --- |
| 2.000 | 0.050 | 0.587 |
| 0.500 | 0.050 | 0.941 |
| 0.125 | 0.050 | 0.978 |
| 0.500 | 0.200 | 0.633 |
| 0.500 | 0.400 | 0.060 |
| 0.500 | 0.785 | 0.014 |

Coarsening $\Delta x$ by a factor of 16 costs a factor of $1.6$, and that loss
is envelope averaging, not fringe washing. Coarsening $\Delta p$ to the fringe
period $`2\pi\hbar/d = 0.785`$ costs a factor of $70$. The resource interference
consumes is momentum resolution.

**Corollary R3.1 (the spec's lattice is exactly critical).** The crystal-lattice
specification fixes $`\Delta p = \pi\hbar/L`$. Nyquist requires at least two
samples per fringe period, i.e. $`2\pi\hbar/d \ge 2\Delta p`$, which is
$d \le L$. So the momentum lattice of the spec resolves the interference of
*every* superposition that fits in the box, and is critically sampled at
$`d = L`$. It cannot be improved on and does not need to be.

This is the sampling-theoretic face of Proposition T2 of
[`takabayasi_1954_stochastic_picture.md`](takabayasi_1954_stochastic_picture.md),
which showed that $`\pi\hbar/L`$ is the exact Wigner support of any
period-$`L`$ state — the reciprocal lattice, not a discretisation. Support and Nyquist agree
because they are the same statement read in the two directions.

---

## 4. Multiplicativity: entanglement is free, non-Gaussianity is not

**Proposition R4.** $`\|W_1 \otimes W_2\|_1 = \|W_1\|_1\,\|W_2\|_1`$.

Immediate from $`\sum_{ij}|a_ib_j| = (\sum_i|a_i|)(\sum_j|b_j|)`$; Part E
verifies it to zero residual. So $n$ independent cats cost $`2.679^n`$, which is
$`3.6\times10^{8}`$ at $n = 20$. The cost is exponential in the number of
factors.

But the factors that matter are not the entangled ones. Part E computes
$`\|W\|_1`$ for coherent, squeezed, displaced and thermal states and gets
$1.000000$ in every case. By **Hudson's theorem** the pure states with
$W \ge 0$ are exactly the Gaussians — and a two-mode squeezed vacuum is
Gaussian, maximally entangled, and costs nothing at all to sample. A
world-particle ensemble represents arbitrary Gaussian entanglement with a plain
positive probability distribution and no negatons whatsoever.

**Corollary R4.1.** Entanglement is not the driver of representation cost;
Wigner negativity is, and by Hudson that means non-Gaussianity.

Two remarks on what this does and does not buy. It does not make the model
efficiently simulable in general — negativity in either the state or the
measurement is exactly the resource that lifts a phase-space model out of
efficient classical simulation (Mari–Eisert 2012; Veitch *et al.* 2012), and by
Ferrie–Emerson (2008) and Spekkens (2008) no quasiprobability representation
escapes it somewhere. That is the representation-theoretic sibling of
Theorem T5 of the Takabayasi note: a normalised, orthogonal, entrywise
non-negative phase-space kernel must be a permutation. It does, however, mean
the exponential is charged per *non-Gaussian* factor rather than per particle,
which is a much weaker growth than the sparse-ontology argument assumes for
branching ensembles.

---

## 5. The pair $(x,\mu)$ ensemble is strictly more expensive

The position-pair ladder samples conjugate pairs and carries the misalignment
$`\mu = \Phi_{\rm ket} - \Phi_{\rm bra}`$ as a phase of unit modulus. Its cost
functional is $`Z = \sum_{ij}|\rho_{ij}| = \bigl(\sum_i|\psi_i|\bigr)^2`$.

**Proposition R5.** On a position lattice of spacing $\Delta x$,
$`Z = (\int|\psi|\,dx)^2/\Delta x`$, so $`Z \propto 1/\Delta x`$ and diverges
under lattice refinement. $Z$ is also basis-privileged: a plane wave costs
$Z = M$ in the position basis and $`\|W\|_1 = 1`$ in the Wigner representation.

Part D:

| $\Delta x/\sigma$ | $Z$ (cat, $d = 8$) | $Z$ (packet) | $Z$ (plane wave) |
| --- | --- | --- | --- |
| 1.000 | 10.02 | 5.01 | 40.0 |
| 0.250 | 40.09 | 20.05 | 160.0 |
| 0.0625 | 160.37 | 80.21 | 640.0 |
| 0.03125 | 320.74 | 160.42 | 1280.0 |

$`\|W\|_1 = 1.59`$ for every row.

The structural reason is worth stating as an identity rather than an inequality.
$W(X,p)$ is the $Y$-Fourier transform of $`\rho(X + Y/2,\, X - Y/2)`$. The
Wigner representation performs the oscillatory $Y$-integral *analytically*; the
pair ensemble leaves it to Monte Carlo. The two representations carry the same
information and the pair ensemble is noisier by exactly the variance of the
integral the transform already did.

**What caps it.** If the pair separation is bounded, $`|Y| \le Y_{\max}`$, so is
$Z$: at $\Delta x = 0.125$ the same state gives $Z = 8.9$ at
$`Y_{\max} = 0.5\sigma`$ rising to $80.2$ only at $`Y_{\max} = 40\sigma`$. A
finite coherence length is therefore not a defect of the pair picture but its
regulator — and since decoherence shrinks the coherence length, the pair
ensemble gets *cheaper* exactly where a hydrodynamic ensemble gets sparser.

**What it buys.** Unit-modulus weights, a manifestly positive particle measure,
conserved particle number, no negative objects, and endpoint locality. That is a
real ontological gain and this section is not an argument against the pair
picture. It is an argument that the gain is paid for in variance, and the price
should be quoted.

**One item to check in `position_pair_ladder.md`.** If the physical density is
read off as a *coincidence* ($Y = 0$), it is a measure-zero event and the cost
rises by a further factor $1/\Delta x$. If instead it is read off the midpoint
$X$ with weight $`\cos\mu`$ integrated over $Y$, there is no coincidence problem
— but then the construction *is* the Wigner ensemble, with the signs recoded as
phases in $`\{0,\pi\}`$. Which of these the ladder intends is not stated
explicitly and should be.

![Sampling cost of the three representations](https://raw.githubusercontent.com/billpage/wpmw/output/figures/representation_cost_sampling.png)

*Left: reconstructing $`\rho(0)`$ for the same cat from a signed ensemble and
from a sea-shifted ensemble. Both fall as $`N^{-1/2}`$; the sea is offset by a
factor of about $`20`$. Centre: pair-ensemble cost against lattice spacing, with
$`\|W\|_1`$ for reference. Right: a coherence cutoff caps it.*

---

## 6. The sea is exact, and has the worst readout of the three

The crystal shift $`W' = W + 2/h`$ is exact — the QLE contains only derivatives
of $W$, and the Wigner bound $`|W| \le 2/h = 1/\pi\hbar`$ is tight, so $W'$ is
admissible and non-negative. Nothing here disputes that. The question is what it
costs to read the state back out.

**Proposition R6.** $`\|W + 2/h\|_1 = 1 + 2A/h`$, where $A$ is the phase-space
area kept in view. For the specification's grid — $M$ position cells over $[-L/2,
L/2)$ and $N$ momentum cells of width $`\Delta p = \pi\hbar/L`$ —

```math
\frac{2A}{h} \;=\; \frac{2\,L\,(N\Delta p)}{2\pi\hbar}
\;=\; \frac{2\,L\,N\pi\hbar/L}{2\pi\hbar} \;=\; N .
```

*The sea mass is exactly the number of momentum cells.* In
`demo_cosine_well_microdynamics.py` that is $256$, against a physical signal of
$1$ — which is what the script's printed `||W'||_1` reports.

**Corollary R6.1 (per-cell error).** Because the sea dominates the counting
statistics, the fractional error of the reconstructed $W$ per cell is
$`\sqrt{N_{\rm cells}/N}`$, *independent of the state*. Part C measures
$0.7075$, $0.2239$, $0.0707$ against predictions $0.6960$, $0.2201$, $0.0696$ at
$`N = 10^5, 10^6, 10^7`$.

**Corollary R6.2 (the crystal condition).** Let $`\mathcal{N}`$ be the number of
world-particles representing unit probability. Sea occupancy per cell is
$`\mathcal{N}\cdot(2/h)\cdot\Delta x\,\Delta p = \mathcal{N}/M`$. One negaton per
cell therefore requires $`\mathcal{N} \ge M`$, and the total ensemble is then
$`\mathcal{N}(1 + N) \approx MN`$ — **exactly one world-particle per phase-space
cell.** The crystal condition and the one-per-cell picture are the same
condition.

Part C, for the same cat and the same target quantity, measures the sea
ensemble's standard deviation at $16.5$–$20.9$ times the signed ensemble's, i.e.
roughly $400$ times as many world-particles for the same accuracy at this box
size, growing linearly in $A$.

So the honest accounting is this. **The sea does not represent negativity
better; it represents it exactly and reads it back worst.** Negativity is
converted from a sign problem into a *contrast* problem: the physical state is a
one-part-in-$`N`$ modulation of a uniform background. What is bought is an
ontology with no negative objects, fixed particle number, and holes as local
rearrangements rather than creations from nothing. What is sold is
signal-to-noise proportional to $1/A$. Both should be on the ledger.

It is also worth noticing that the crystal ontology is the exact opposite of a
sparse one. It demands at least one world-particle per phase-space cell
everywhere in view, including regions where nothing whatever is happening.
Whatever objection Hackebill and Poirier's argument raises, thinness is not it.

---

## 7. The annihilation burden: where the program stands

### 7.1 Three models, one cost

The three microdynamic models in the ladder handle the sign differently, and
each pays $`\|\cdot\|_1`$ in a different currency.

| Model | Sign carried by | Growth control | Currency |
| --- | --- | --- | --- |
| Spawning | positon/negaton pairs created at vertices | garbage collection merging near-coincident pairs | walker population |
| Jumping | fixed populations, retrocausal momentum jumps | requires a spawning model with GC to simulate | same, deferred |
| Crystal lattice | positon excess over a pinned negaton sea | none needed | background mass $2A/h$ |

The crystal model *solves* the annihilation problem, by paying for it up front
in §6's currency. That is the trade the ladder currently offers, and until now
neither side of it had been priced.

### 7.2 The burden, measured

The QLE jump term for one mode is $`\Gamma_q(x)\,(W_{n+q} - W_{n-q})`$ with
$`\Gamma_q(x) = -(V_q/\hbar)\sin(2\pi q x/L + \phi_q)`$. Replacing the stencil
entrywise by its absolute value gives the growth rate of the pathwise $L^1$ mass
of a non-annihilating unraveling: the column sums are $`2|\Gamma_q(x)|`$, so

```math
\gamma_{\max} = \frac{2V_q}{\hbar}, \qquad
\gamma_{\rm avg} = \left\langle 2|\Gamma_q|\right\rangle_x
= \frac{4V_q}{\pi\hbar}.
```

Part F evolves a cat state under the exact QLE in a one-mode cosine well
($L = 8$, $M = 256$, $N = 128$, $q = 1$, $V_q = 1.5$) and tracks the true
$`\|W(t)\|_1`$ against these bounds. Probability norm is conserved to
$`2.2\times10^{-16}`$.

| $t$ | $`\|W(t)\|_1`$ | $`e^{\gamma_{\rm avg}t}`$ | burden |
| --- | --- | --- | --- |
| 0.0 | 1.210 | 1.00e+00 | 8.3e-01 |
| 1.0 | 1.132 | 6.75e+00 | 6.0e+00 |
| 2.0 | 1.534 | 4.56e+01 | 3.0e+01 |
| 3.0 | 2.056 | 3.08e+02 | 1.5e+02 |
| 4.0 | 1.939 | 2.08e+03 | 1.1e+03 |

$`\|W(t)\|_1`$ stays in $[1.115, 2.100]$ — bounded, and oscillating with the
recurrence of the well. The pathwise bound grows exponentially. **The ratio is
the $L^1$ mass that garbage collection must remove per unit time**, and it is a
factor of a thousand after four time units.

![The annihilation burden](https://raw.githubusercontent.com/billpage/wpmw/output/figures/annihilation_burden.png)

*Left: the state's own $`L^1`$ mass under exact QLE evolution — bounded. Right:
the same curve against the two pathwise growth bounds, log scale. The gap is the
burden.*

Two things about this number deserve emphasis, one reassuring and one not.

**Reassuring:** $\gamma$ is set by $`V_q/\hbar`$ and nothing else. It is a
*potential-strength* rate, not a spectral gap. In fermionic diffusion Monte
Carlo the sign-problem exponent is the gap between the bosonic and fermionic
ground-state energies, which grows with system size; here it is bounded by the
Fourier amplitude of the potential. **The WPMW sign problem is mild by the
standards of the QMC literature**, and for a weak mode it is very mild.

**Not reassuring:** it is still exponential in $t$, and there is at present no
demonstration anywhere in the repository that any annihilation scheme actually
holds it.

### 7.3 Why annihilation is structurally awkward here

Annihilation requires a positon and a negaton to *coincide*. On the momentum
axis this is unproblematic: the lattice $`\pi\hbar/L`$ is discrete and every
action moves particles by an integer number of cells. On the position axis it is
not. Free streaming is continuous in $x$, and exact coincidence of two
continuously-streaming particles is a measure-zero event. Annihilation is
therefore only well defined once the position axis is *also* discretised — that
is, once the phase-space crystal is imposed. This is a point in favour of the
crystal model that is independent of §6's cost, and it should be recorded as
such: **the crystal is not merely one way to organise the sea, it is the
structure that makes pairwise annihilation definable at all.**

### 7.4 Two candidate density criteria, and they disagree

Part F evaluates both for the run above.

**(i) Occupancy.** By Corollary R6.2, sea occupancy per cell is
$`\mathcal{N}/M`$. Annihilation needs partners available, so it needs
$`\mathcal{N} \gtrsim M = 256`$, a total ensemble of $32{,}768$. By this
criterion `demo_cosine_well_microdynamics.py`, with $`5\times10^6`$ positons over
$256^2$ cells, is at occupancy $76$ per cell — comfortably inside the
annihilating regime, with a factor of $76$ to spare.

**(ii) Partner separation.** A focus vertex creates a pair at the same $x$ with
momenta differing by $`2q\Delta p = 2\pi\hbar q/L`$. They separate at relative
velocity $`2\pi\hbar q/(mL)`$ and leave a position cell of width $L/M$ after

```math
t_{\rm sep} \;=\; \frac{m L^2}{2\pi\hbar\, q\, M} \;=\; 0.0398,
```

against a mean creation interval $`\hbar/V_q = 0.667`$. The dimensionless
**annihilation number** $`\mathcal{A} = t_{\rm sep}/t_{\rm create} = 0.060 \ll
1`$: original partners escape their cell about seventeen times faster than new
ones are made. Note that $`t_{\rm sep} \propto 1/M`$, so **refining the position
lattice makes this criterion worse, linearly**, which is the opposite of the
direction §3 recommends for resolution.

These two answers are not reconcilable by inspection. Criterion (i) says the
cells are richly populated and annihilation will find partners easily; criterion
(ii) says the specific partners whose cancellation the dynamics requires are
long gone. Which governs depends on whether annihilation against an *anonymous*
opposite-sign particle in the same cell is dynamically equivalent to
annihilation against the original partner. In the mean-field QLE it is —
particles carry no identity beyond $(x, p)$ and species. In a pathwise
unraveling with a correlated sign structure it need not be, and that is exactly
the regime where the FCIQMC literature finds a plateau.

### 7.5 The FCIQMC analogy, and its limits

Full configuration-interaction QMC controls its sign problem by annihilating
walkers on identical determinants, and finds a **critical walker number** — the
plateau (Spencer, Blunt & Foulkes 2012) — below which annihilation cannot
stabilise the sign structure and the simulation converges to the wrong answer.
Above it, it converges to the right one. This is the closest existing analogue
of a sparse-ontology threshold in a signed particle method, and it is a genuine
threshold rather than a smooth degradation.

Where the analogy holds: the mechanism (annihilation needs coincidence, which
needs density) and the phenomenology (a threshold, not a gradient).

Where it breaks: FCIQMC's determinant space is discrete and finite, so
coincidence is generic; WPMW's is continuous in $x$ unless the crystal is
imposed (§7.3). And FCIQMC's sign-problem exponent grows with system size while
WPMW's is bounded by $`2V_q/\hbar`$ (§7.2). The WPMW plateau, if it exists, should
therefore sit at a much lower density than the FCIQMC one — but "should" is
doing all the work in that sentence, and no one has measured it.

---

## 8. Where to go

Five items, in the order they should be taken.

**(N1) Ask whether the four-action set is $L^1$-stationary in its own right.**
This is the sharpest single question and it may dissolve the problem. The focus
action $`(n, n) \to (n-q, n+q)`$ has an exact reverse, and the reverse *is* an
annihilation channel: it removes precisely the pair the forward channel made.
Since $`\|W(t)\|_1`$ is bounded (§7.2), the mean-field answer is that creation
and recombination balance. The open question is whether an unraveling can be
*chosen* so that the pathwise $`\sum_i |w_i|`$ is also bounded — i.e. whether
the four actions already contain their own garbage collection and no external
scheme is required. This is the same question as the standing open item on
whether the SLN-type unraveling of $\rho$ is equivalent to stochastic
two-state-vector dynamics, approached from the $L^1$ side, and it should be
attacked with the broken detailed balance of
[`../analysis/sea_dressed_microdynamics.md`](../analysis/sea_dressed_microdynamics.md)
in hand, since that is what supplies the asymmetry between the forward and
reverse rates.

**(N2) Measure the plateau.** The experiment that settles §7.4: a controlled
one-mode run with an explicitly annihilating particle process, sweeping
$`\mathcal{N}`$ over three decades and watching for a critical value below which
the reconstructed $`\rho(x)`$ at an interference minimum stops converging to the
exact answer. The cat in the cosine well of Part F is already the right test
case, and Part C already supplies the exact target and the noise floor. If a
plateau exists, criterion (ii) governs; if the error falls as $`N^{-1/2}`$ all
the way down to occupancy of order unity, criterion (i) does.

**(N3) Decide the annihilation cell, and price the three options.** Cell-exact
annihilation on the crystal (unbiased, needs the position lattice); soft or
kernel annihilation over a finite phase-space blob (works off-lattice,
introduces a bias controlled by the blob size); no annihilation at all with the
sea absorbing the sign (unbiased, costs $2A/h$). These are three inequivalent
regularizations of the same divergence and no document in the repository
compares them. Sellier's signed-particle formulation (*J. Comput. Phys.* **297**,
254, already in the bibliography as the closest published relative of the
positon/negaton ontology) supplies a worked same-cell annihilation scheme, and
the bibliography entry already notes that this project does not yet use it;
adopting it is the obvious starting point for the first of the three. They should be compared in one place, with the bias and the
variance of each quoted on the same axes.

**(N4) Test the decoherence direction as a prediction.** §5 and §6 both imply
that this model gets *cheaper* under decoherence — $`\|W\|_1 \to 1`$ as fringes
are smoothed, $Z$ falls with the coherence length, and the sea contrast problem
is unchanged. That is the exact opposite of the sparse-ontology trajectory, and
it is measurable: add a second mode or a weak bath to the Part F run and plot
$`\|W(t)\|_1`$. It should decrease monotonically toward $1$. If it does not, §9
is wrong and this should be found out early.

**(N5) Audit the collective quantities for density dependence.** §9 argues that
the kinetic model escapes the sparse-ontology argument because its rates depend
only on the external potential. That argument is sound for the $G = 0$
single-mediated-jump family and for the four actions. It is *not* obviously
sound for any formulation whose rates read a collective order parameter — the
$Z_r$ of the phase-alignment layer, or an anonymous beat grating built from the
whole ensemble. Wherever a rate is a functional of the ensemble rather than of
the potential, a thinning argument re-enters and the immunity claimed in §9 has
to be re-derived. This audit has not been done.

---

## 9. The sparse-ontology argument does not transfer

Hackebill and Poirier's mechanism has a precise structural precondition: the
dynamics must be **nonlinear in the ensemble density**. In the hydrodynamic
formulation this is explicit — $v = j/\rho$ and
$`Q = -(\hbar^2/2m)\nabla^2\sqrt{\rho}/\sqrt{\rho}`$ both divide by $\rho$ — and
in discrete MIW it is the interworld force, a finite-difference functional of
neighbouring world *spacings* (§2.1 of the foundations note reproduces
$`U_N`$ and identifies its bracket as Nelson's osmotic momentum). Three
consequences follow: a single world has no dynamics; the law of motion is a
functional of the ensemble; and thinning corrupts the force itself rather than
merely its estimate. Branching under decoherence then drives each sub-ensemble
away from the quantum regime. That is the argument, and within its
preconditions it looks right.

The Wigner/kinetic dynamics divides by nothing. The QLE is linear in $W$, and in
the four-action formulation each world-particle's transition rate is a function
of the external potential at its own position — this is exactly the endpoint
locality that Theorem F3 identified as the model's defining postulate, and it is
also what David's own note records as the disappearance of any dependence on
local density derivatives. A lone world-particle already has well-defined
dynamics. Thinning the ensemble therefore degrades the *statistics* as
$`N^{-1/2}`$, without bias and without threshold, and never changes the
equations of motion.

Four further asymmetries, all quantified above:

1. **Decoherence runs the other way.** $`\|W\|_1 \to 1`$ as fringes are smoothed
   (§4, §5), so the sampling penalty vanishes exactly in the regime where the
   sparse-ontology problem appears. Each decohered branch becomes an ordinary
   positive probability ensemble.
2. **The worlds are in the right place.** Corollary R2.1: the phase-space
   ensemble populates the interference region, which carries 38.9% of the
   content and which MIW leaves empty.
3. **The cost saturates rather than compounding.** Proposition R2: a bigger cat
   is not a more expensive cat.
4. **Entanglement is free** (Corollary R4.1). Branching into entangled Gaussian
   sub-ensembles costs nothing at all.

What *does* transfer is a threshold of a different kind and in a different
place: the annihilation density of §7. If the program's growth control depends
on pairwise cancellation, it depends on a local density of opposite-sign
particles, and below some occupancy that control fails. That is a real
structural cost and it is the honest answer to "does this model have a
sparse-ontology problem": not the one Hackebill and Poirier describe, but one
that lives in the unraveling rather than in the ontology, and whose critical
density is currently unmeasured.

One caveat, carried forward as N5. This section's immunity argument is only as
good as the claim that rates depend on the potential and not on the ensemble.
That claim holds for the four actions. It has not been checked for the
collective order parameter of the phase-alignment layer.

---

## 10. Numerical verification

`src/demo_representation_cost.py`, run as
`WPMW_OUTPUT=... python src/demo_representation_cost.py`.

| Part | Claim verified |
| --- | --- |
| A | Closed-form cat Wigner function against FFT evaluation, max diff $`1.7\times10^{-16}`$; $`\|W\|_1 \to 1+2/\pi`$; interference share $`38.90\%`$ |
| B | Fringe visibility against $\Delta x$ and $\Delta p$ (Proposition R3) |
| C | Signed vs sea sampling; $`\|W+2/h\|_1 = 1 + 2A/h`$; per-cell rule $`\sqrt{N_{\rm cells}/N}`$ to within $`2\%`$ |
| D | $`Z \propto 1/\Delta x`$; basis privilege; coherence cutoff |
| E | $`\|W\|_1 = 1`$ for every Gaussian; multiplicativity to zero residual |
| F | $`\|W(t)\|_1`$ bounded under exact QLE, norm conserved to $`2.2\times10^{-16}`$; $`\gamma_{\max}, \gamma_{\rm avg}`$; both density criteria |

Figures: `representation_cost_cat.png`, `representation_cost_sampling.png`,
`annihilation_burden.png`, published to the `output` branch.

---

## 11. Summary

1. The cost of any world-particle representation is $`N/\|\mu\|_1^2`$
   (Definition R1). Everything else in this note is a computation of
   $`\|\mu\|_1`$.
2. A cat state costs a factor $2.679$ in samples and **the cost saturates**:
   $`\|W\|_1 \to 1 + 2/\pi`$ however large the cat (Proposition R2).
3. 38.90% of the ensemble must sit in the interference lobe, where the
   probability density is numerically zero (Corollary R2.1). This is the
   resource MIW lacks, not an embarrassment.
4. Interference consumes momentum resolution, not position density
   (Proposition R3), and the spec's lattice $`\Delta p = \pi\hbar/L`$ is exactly
   critically sampled for every state that fits the box (Corollary R3.1).
5. Cost is multiplicative across factors, but every Gaussian costs $1$ —
   entanglement is free and non-Gaussianity is what is charged for
   (Corollary R4.1).
6. The pair $(x,\mu)$ ensemble is the Wigner ensemble with the $Y$-integral left
   to Monte Carlo. Its cost diverges as $`1/\Delta x`$ and is capped only by a
   finite coherence length (Proposition R5).
7. The sea is exact but reads back worst: $`\|W+2/h\|_1 = 1 + 2A/h`$, which for
   the spec's grid is exactly the number of momentum cells (Proposition R6). The
   crystal condition $`\mathcal{N} \ge M`$ is the same as one world-particle per
   phase-space cell (Corollary R6.2).
8. The annihilation burden is now measured: $`\gamma_{\rm avg} = 4V_q/\pi\hbar`$
   against a bounded $`\|W(t)\|_1`$, a factor of $10^3$ after four time units.
   It is mild by QMC standards but no scheme in the repository has been shown to
   hold it.
9. Two density criteria for the annihilating regime disagree by two orders of
   magnitude (§7.4). Deciding between them is item N2.
10. The sparse-ontology argument does not transfer, because the kinetic dynamics
    is linear in the ensemble and divides by nothing (§9). A different threshold
    does transfer, and it lives in the unraveling.

**Open items.**

- N1: is the four-action set $L^1$-stationary pathwise, or only in mean field?
- N2: measure the plateau; decide between the occupancy and separation criteria.
- N3: compare the three annihilation regularizations on one set of axes.
- N4: verify that $`\|W(t)\|_1`$ decreases under decoherence, as §9 predicts.
- N5: audit $Z_r$ and any grating-mediated rate for ensemble dependence; §9's
  immunity argument does not automatically extend to them.
- Whether `position_pair_ladder.md` reads the density off coincidences or off
  midpoints (§5); the cost differs by a further factor $1/\Delta x$.
