# What the reach is

**A tutorial on the one parameter everything in the compensated split depends
on — what is separated by $`2y`$, why the reach is a period rather than an
aperture, why it aliases coherence instead of destroying it, and why there are
two different ceilings on it and only one of them is a ceiling on the theory.**

---

## 0. Status and provenance

Tutorial supplement. Companion demo: `src/demo_reach_tutorial.py`, which
produces every table and figure below.

Prompted by Bill Page's observation that the reach "appears to sit
uncomfortably between the representation/numerical and the ontological
pictures." It does, and §6 argues that the discomfort is well founded and
locates something real rather than a gap in anyone's understanding.

Nothing here is new mathematics. The results are Theorems E1, E1.1, E3, E6,
E8 of [`../analysis/reach_energy_coupling.md`](../analysis/reach_energy_coupling.md),
Theorems C3, C4, C7 of
[`../analysis/compensated_liouville_splitting.md`](../analysis/compensated_liouville_splitting.md),
Definition (R) and Propositions O3, O4 of
[`../analysis/open_position_space.md`](../analysis/open_position_space.md),
and Theorem K1 of
[`../analysis/eckart_barrier_compensated.md`](../analysis/eckart_barrier_compensated.md),
reorganised around the questions a reader actually asks and in the order they
ask them. §5 is the one place where the reorganisation produces a claim the
notes do not contain, and §8 records what would have to change elsewhere if it
survives scrutiny. **By decision, no existing document is edited by this
supplement**; the corrections are logged, not applied.

---

## 1. First, what is separated

The single most common error — and this document's author has made it in
writing, more than once — is to read

> the ket arm sits at $`x+y`$, the bra arm at $`x-y`$, separated by $`2y`$

as two bodies at two places. There is one body.

$`W`$ and $`\rho`$ are related by

```math
W(x,p) \;=\; \frac{1}{\pi\hbar}\int \rho(x+y,\;x-y)\;e^{-2ipy/\hbar}\,dy ,
```

so $`y`$ is the off-diagonal coordinate of the density matrix. The two
positions $`x \pm y`$ are two *worlds'* answers to where the one body is. The
quantity $`2y`$ is a distance in configuration space between branches, not a
distance between objects.

**The analogy that works.** An interferometer. Two arms, one photon. The path
difference is not a separation between particles, and nobody is tempted to
think it is. A ket–bra pair at separation $`2y`$ is a path difference of
exactly that kind.

**The word that works.** `../analysis/position_pair_ladder.md` already supplies
it: the ends are **legs**, each carrying "a place and a clock."

![The geometry, and the two ceilings](https://raw.githubusercontent.com/billpage/wpmw/output/figures/reach_tutorial_geometry.png)

Left: three pairs at three midpoints, each drawn as two legs and a midpoint
cross, inside the horizon band. The reach is a bound in the *vertical*
coordinate — it is orthogonal to position and is not a wall in space. Right:
§5.

### 1.1 The ladder has already proved the point harder

The position-pair ladder's central result is

```math
\bar p \;=\; \frac{\hbar\thinspace\mu}{a},
\qquad \mu \;=\; \arg\rho(X, X') ,
```

so momentum is not carried by anything. It *is* the misalignment of a pair of
legs. The current between neighbouring legs is $`J|\rho_1|\sin\mu/\hbar`$
(Part A):

```
        mu/pi     sin mu     current    interpretation
         0.00    0.00000    0.00000    self-conjugate: motionless
         0.25    0.70711    0.70711
         0.50    1.00000    1.00000    zone boundary, maximal
         0.75    0.70711    0.70711
         1.00    0.00000    0.00000    antipodal: motionless again
```

A *self-conjugate* pair — legs coincident, $`\mu = 0`$ — is the object that
looks most like a classical particle, and it is the one that cannot move. A
population consisting only of self-conjugate pairs has exactly zero flux. The
particle intuition therefore fails at the very first thing one would ask a
particle to do.

This matters for everything that follows, because a period in $`y`$ reads as
an absurdity under the particle picture ("space is a circle") and as something
quite mild under the leg picture (§3).

---

## 2. Is the reach a coherence length?

Nearly. The difference is the whole subject.

A *coherence length* means attenuation: $`\rho(x+y, x-y) \to 0`$ as $`|y|`$
grows. That is decoherence, and it has consequences — purity falls,
interference dies. Test both candidates on a cat state,
$`|\psi\rangle \sim \mathcal{N}(-d/2) + \mathcal{N}(+d/2)`$ with $`d = 4`$,
$`\sigma = 0.5`$, whose coherences live at $`x \approx 0`$, $`y \approx \pm2`$.

**(i) Aperture** — delete $`\rho`$ for $`|y| > y_{\max}`$, the textbook toy
model of decoherence. **(ii) Fold** — identify $`y`$ with $`y + L_c`$, which
is what the algorithm's DFT does. Nothing is deleted.

| $`y_{\max}`$ | aperture: $`\mathrm{Tr}\rho^2`$ | aperture: fringe | fold: $`\mathrm{Tr}\rho^2`$ | fold: fringe |
|---|---|---|---|---|
| 0.50 | $`0.422444`$ | $`0.000808`$ | $`2.723890`$ | $`0.318417`$ |
| 1.00 | $`0.499245`$ | $`0.013532`$ | $`1.592920`$ | $`0.318417`$ |
| 1.50 | $`0.540454`$ | $`0.094971`$ | $`1.207491`$ | $`0.410335`$ |
| 2.00 | $`0.752080`$ | $`0.295722`$ | $`1.501342`$ | $`0.552170`$ |
| 3.00 | $`0.999527`$ | $`0.549492`$ | $`1.009829`$ | $`0.457091`$ |
| 6.00 | $`1.000671`$ | $`0.556425`$ | $`1.000671`$ | $`0.552170`$ |

The aperture behaves exactly like decoherence: purity collapses from $`1.0007`$
to $`0.4224`$, the interference fringe from $`0.5564`$ to $`8\times10^{-4}`$.
The fold destroys neither — the fringe survives at $`y_{\max} = 0.5`$, a
quarter of the coherence separation.

![Purity and interference under aperture and fold](https://raw.githubusercontent.com/billpage/wpmw/output/figures/reach_tutorial_coherence.png)

### 2.1 What the fold does instead, and where the excess goes

$`\mathrm{Tr}\rho^2 = 2.72`$ is impossible for a density matrix. The folded
object is not one — not on the line. It is a perfectly good density matrix on
the **circle**, and the excess has a precise origin.

In ladder language the fold identifies leg separation $`k`$ with $`k + N`$, so
a coherence at $`y = d/2`$ lands exactly on the diagonal whenever $`L_c`$
divides $`d/2`$. Part B(iii):

| $`y_{\max}`$ | $`L_c`$ | true diagonal peak | folded peak | ratio |
|---|---|---|---|---|
| 0.50 | 1.00 | $`0.398942`$ | $`1.014724`$ | $`2.544`$ |
| 1.00 | 2.00 | $`0.398942`$ | $`0.798688`$ | $`2.002`$ |
| 1.50 | 3.00 | $`0.398942`$ | $`0.398979`$ | $`1.000`$ |
| 2.00 | 4.00 | $`0.398942`$ | $`0.398942`$ | $`1.000`$ |

Here $`d/2 = 2`$, so the contamination is exact at $`L_c = 1`$ and $`L_c = 2`$
and absent at $`L_c = 3`$. Note the shape of that: **a longer reach is dirty
where a shorter one is clean.** It is commensuration, not magnitude — which is
the same criterion §3 finds in the potential, and the two being the same
phenomenon is §8's open item.

**The verdict.** The reach is not a decoherence length. It makes coherences
whose separations differ by $`L_c`$ indistinguishable from one another, and at
commensurate separations it books a coherence as probability density. That is
a *resolution limit on the relation between legs*, not a decay. Two cat states
with separations $`d`$ and $`d + L_c`$ have the same representation.

---

## 3. Aperture or period? The measurement that decides

The aperture reading survives §2 with a bruise. §1 of
`../analysis/reach_energy_coupling.md` kills it outright, and it is worth
seeing why in terms of a prediction the two pictures disagree about.

**The argument.** The momentum transfer $`\xi`$ and the half separation $`y`$
are Pontryagin duals. A discrete $`\xi`$ spectrum requires $`y`$ compact *as a
group* — a circle, not an interval. A compactly supported function on the line
has, by Paley–Wiener, an entire and perfectly continuous transform. Finite
support gives no lattice. Only a period does.

**The prediction they disagree about.** For a crystal potential
$`V = \cos(2\pi x/a)`$ the exact budget is a pair of deltas. An aperture
picture says a soft taper should help, by suppressing ringing at the window
edge. A period picture says only *commensuration* matters, and that a taper is
positively harmful because it attenuates the true delta. Part C:

| $`L_c`$ | $`L_c/a`$ | commensurate | sharp | tapered |
|---|---|---|---|---|
| 4 | 1.00 | yes | $`1.0000`$ | $`0.7500`$ |
| 8 | 2.00 | yes | $`1.0000`$ | $`1.0000`$ |
| 6 | 1.50 | no | $`6.2322`$ | $`0.9681`$ |
| 5 | 1.25 | no | $`4.6688`$ | $`0.8788`$ |

The period picture wins on both counts. Sharp-and-commensurate is exact;
sharp-and-incommensurate is wild; and at $`L_c = a`$ the taper loses a quarter
of the budget, which no aperture picture predicts.

### 3.1 Born–von Kármán, and where the analogy breaks

The right comparison is the Born–von Kármán boundary condition of solid-state
physics: impose a period, get a reciprocal lattice, and $`\Delta p`$ follows.
Three of the four moving parts map over exactly.

What does **not** map over is the limit. In solid state one takes
$`L \to \infty`$ at the end; the torus is a scaffold and it evaporates. Here
that limit destroys the object it was introduced to support: at
$`y_{\max} = \infty`$ the residual symbol is unbounded and there is no jump
measure at all (Theorem C3), and for a polynomial potential the microdynamics
does not exist on the open line under any refinement (Theorem E7). **The
compactification is constitutive.** There is no "line" version of this theory
for the periodised one to approximate.

And one warning the analogy carries in the wrong direction, which is §1's
point again: BvK compactifies *actual space*. The reach does not. What is
compact here is the space of leg separations — the relation — and position
space is untouched, which is the entire content of
`../analysis/open_position_space.md` retracting the box.

---

## 4. One dial, three consequences

Reach and momentum quantum are not two parameters:

```math
\Delta p \;=\; \frac{\pi\hbar}{2\thinspace y_{\max}}
\qquad\text{so}\qquad
\Delta p \cdot y_{\max} = \frac{\pi\hbar}{2} \;\;\text{exactly, always.}
```

Three quantities move when that dial turns, and they do not move together.

| $`y_{\max}`$ | $`\Delta p`$ | Moyal series | event budget |
|---|---|---|---|
| small | coarse | convergent | cheap |
| large | fine | divergent past $`R(x)`$ | expensive, divergently (E8) |

And a fourth column is empty. Theorems E3 and E6 say signed world number,
classical force, energy, and the leading $`\hbar^2`$ coefficient are all
exactly reach-independent — the last one *under compensation*, which is what
makes §4.4's normative rule load-bearing rather than bookkeeping.

![One dial, three consequences](https://raw.githubusercontent.com/billpage/wpmw/output/figures/reach_tutorial_tradeoff.png)

> **The reach buys momentum resolution and costs event budget. It buys no
> accuracy.** So choosing $`L_c`$ is a resolution-versus-cost decision, not a
> convergence study. There is nothing to refine.

---

## 5. Two ceilings, and only one is a ceiling on the theory

The notes use "the reach ceiling" for two different things. They are not the
same number and it is worth separating them.

**Ceiling 1, existence.** The kernel is built from $`V(x \pm y)`$ on the
**real** segment $`[x - y_{\max},\, x + y_{\max}]`$. A singularity of $`V`$ on
that segment kills it. Separately, an unbounded reach has no jump measure
(C3). Both are ceilings on the model.

**Ceiling 2, convergence.** Theorem K1: the Moyal series in $`y`$ converges iff
$`y_{\max} < R(x)`$, the distance to the nearest **complex** singularity. For
$`V_0\,\mathrm{sech}^2(z/a)`$, $`R(x) = \sqrt{x^2 + (\pi a/2)^2}`$.

These are routinely conflated, and the second is far tighter. What actually
happens when it is crossed? Part D, sech² at $`a = 1`$, $`x = 0.3`$, where
$`R = 1.599188`$ and $`\hbar^2V'''(x)/4 = 0.930603`$:

| $`y_{\max}`$ | $`y_{\max}/R`$ | $`M_0`$ | $`M_1`$ | $`M_3`$ | budget |
|---|---|---|---|---|---|
| 0.400 | 0.250 | $`0`$ | $`1.1\times10^{-11}`$ | $`0.953133`$ | $`0.0075`$ |
| 0.800 | 0.500 | $`0`$ | $`-6.6\times10^{-13}`$ | $`0.930128`$ | $`0.0526`$ |
| 1.599 | 1.000 | $`-1.4\times10^{-17}`$ | $`-3.4\times10^{-12}`$ | $`0.930178`$ | $`0.2743`$ |
| 2.400 | 1.501 | $`5.6\times10^{-17}`$ | $`8.2\times10^{-13}`$ | $`0.930647`$ | $`0.5552`$ |
| 6.400 | 4.002 | $`2.2\times10^{-16}`$ | $`-7.7\times10^{-13}`$ | $`0.930618`$ | $`2.2681`$ |
| 12.800 | 8.004 | $`-4.4\times10^{-16}`$ | $`2.5\times10^{-13}`$ | $`0.930623`$ | $`4.6165`$ |

Nothing in that table knows where $`R`$ is. World number and momentum stay
exact (C3); the third moment stays on $`\hbar^2V'''/4`$ (E6) inside and
outside alike, and is in fact an order of magnitude *closer* at four times
$`R`$, because the lattice resolves more of the kernel's Fourier content.
Only the budget moves.

Meanwhile the series at the same $`x`$, partial sums against the exact
residual field:

| $`y`$ | $`y/R`$ | exact | $`S_{10}`$ | $`S_{20}`$ | $`S_{30}`$ |
|---|---|---|---|---|---|
| 0.500 | 0.313 | $`1.312\times10^{-1}`$ | $`1.312\times10^{-1}`$ | $`1.312\times10^{-1}`$ | $`1.312\times10^{-1}`$ |
| 1.000 | 0.625 | $`6.891\times10^{-1}`$ | $`6.351\times10^{-1}`$ | $`6.879\times10^{-1}`$ | $`6.891\times10^{-1}`$ |
| 1.550 | 0.969 | $`1.4666`$ | $`-3.519`$ | $`-6.107`$ | $`3.809`$ |
| 1.700 | 1.063 | $`1.6673`$ | $`-1.09\times10^{1}`$ | $`-4.56\times10^{1}`$ | $`4.16\times10^{1}`$ |
| 2.200 | 1.376 | $`2.2870`$ | $`-1.57\times10^{2}`$ | $`-7.46\times10^{3}`$ | $`1.01\times10^{5}`$ |

The exact column is smooth and bounded throughout, because $`V`$ is
real-analytic on the real line and its poles sit at $`\pm i\pi/2`$. The kernel
never leaves the real axis; only the *series* needs the disc.

**The contrast that shows ceiling 1 is real.** Pure Coulomb $`V = 1/z`$ at
$`x = 1`$ has its pole *on* the real axis, at $`y = 1`$:

| $`y_{\max}`$ | Coulomb $`\max\lvert D_{\rm res}\rvert`$ | soft core, $`\epsilon = 0.3`$ |
|---|---|---|
| 0.500 | $`3.33\times10^{-1}`$ | $`0.1825`$ |
| 0.900 | $`7.66`$ | $`1.0605`$ |
| 0.990 | $`9.65\times10^{1}`$ | $`1.1126`$ |
| 1.500 | $`2.00\times10^{4}`$ | $`1.3178`$ |
| 3.000 | $`1.00\times10^{4}`$ | $`5.0267`$ |

The soft core has $`R(1) = 1.044`$ and sails three times past it without
noticing, because its poles are at $`\pm i\epsilon`$ and the real axis is
clear. Pure Coulomb does not.

> **$`R(x)`$ is not a ceiling on the theory. It is the boundary of the
> semiclassical regime.** Inside it, the residual channel is a convergent
> series of $`\hbar`$-corrections and one may say "quantum mechanics is
> classical mechanics plus corrections." Outside it the channel is still
> exactly right — same conservation laws, same $`\hbar^2`$ coefficient — but
> no longer decomposes that way.

Combined with §4 this has teeth. Since $`\Delta p = \pi\hbar/2y_{\max}`$, fine
momentum resolution *requires* leaving the semiclassical regime. You cannot
have both a fine momentum lattice and a convergent semiclassical expansion,
because the two pull $`y_{\max}`$ in opposite directions.

---

## 6. So where does the reach sit?

The discomfort that prompted this document is well founded, and §4 and §5
locate it.

The reach is the single dial that trades **semiclassical intelligibility**
against **momentum resolution**. The first is a representational virtue: it is
about whether the theory can be *described* as a correction to classical
mechanics. The second is ontological: $`\Delta p`$ is the granularity of what
momentum transfers *exist*. One parameter with a foot in each, and no third
option — so the awkwardness is a property of the model, not a confusion to be
cleared up.

What can be cleared up is the rest of it, and the cleaned-up statement is:

- The reach is a **period**, not an aperture (§3), of the **leg separation**,
  not of space (§1).
- It **aliases** coherence rather than destroying it (§2), so it is a
  resolution limit and not a decoherence length.
- It is **constitutive**: the limit that would remove it does not exist (§3.1).
- It buys resolution, costs budget, and buys no accuracy (§4).
- Its two ceilings are different, and only one bounds the theory (§5).

---

## 7. Words that mislead

A glossary of the terms most likely to import the wrong picture. The right
column is what the term actually denotes here.

| term | what it suggests | what it is |
|---|---|---|
| **particle** | a body with a position and a velocity | in the compound *world-particle*, a signed sample of $`W`$. Bare "particle" in a single-body note almost always means *world*, and where it means *leg* it is doubly wrong |
| **pair** | two bodies | **two things in this project.** A *ket–bra pair* is a relation between two worlds about one body. A *sea pair* is two world-particles, a positon and a negaton, bound together. Different in kind; R-SP6 below |
| **separation** $`2y`$ | distance between two objects | a path difference between two branches. Nothing is at both ends at once |
| **arm**, **leg** | a limb of something | one end of a ket–bra pair: a place and a clock. *Leg* is the project's term and the better one |
| **reach** | how far something can see or act | the **period** of the leg separation. Not an aperture, not a distance a world "consults" over — that wording is retracted by E1 and survives in Definition (R) only with a superseding note |
| **coherence horizon** | a distance past which coherence is lost | a distance past which coherences become *indistinguishable*. Nothing is lost; §2 |
| **reach ceiling** | one number | **two numbers**, and the tighter one is not a bound on the theory; §5 |
| **self-conjugate particle** | an ordinary classical particle | a pair with coincident legs, which is the one object in the theory that cannot move; §1.1 |
| **momentum** | something a body carries | the misalignment $`\bar p = \hbar\mu/a`$ of a pair of legs. Nothing carries it |
| **absorbing boundary** | a sponge that eats worlds | the statement that beyond the horizon a world is exactly free (Proposition O4) |
| **truncation** | an approximation to be refined away | for an unbounded $`V`$, constitutive: there is nothing to refine towards (E7) |

---

## 8. Open items

- **R-SP1 (two ceilings, and what it does to K1.1 and Z4).** §5's claim, if it
  survives scrutiny, softens two stated obstructions. Corollary K1.1 of
  `../analysis/eckart_barrier_compensated.md` reads "the reach ceiling requires
  the packet to be narrower in position than half the barrier width," and
  Theorem Z4 of `../analysis/soft_core_coulomb.md` derives
  $`\epsilon \ge k^4\pi^4/4\,a_0`$ from $`y_{\max} < \inf_x R(x)`$. Both are
  stated against ceiling 2. Both potentials have their singularities *off* the
  real axis, so ceiling 1 does not bind, and the obstruction is to holding a
  convergent Moyal series rather than to running the model. Restated: a
  uniform lattice that *also keeps the series convergent over the core* needs
  $`\epsilon \ge k^4\pi^4/4`$; drop that requirement and the threshold
  vanishes, at the price of working non-perturbatively. **Not applied to
  either note by decision.** Before it is, someone should check whether series
  convergence is load-bearing somewhere in Cyganski's framework in a way these
  notes do not record — if it is, the whole of §5 collapses.
- **R-SP2 (is §2.1 the same phenomenon as §3?).** The diagonal contamination
  is exact at $`L_c`$ commensurate with the coherence separation and absent
  otherwise; Proposition E1.1's budget is exact at $`L_c`$ commensurate with
  the potential's period and wild otherwise. Same criterion, one in the state
  and one in the potential. If they are one phenomenon it should be possible
  to state it once, in the separation coordinate, for an arbitrary periodic
  structure in $`y`$. This was noticed while drafting and is not established.
- **R-SP3 (does the fold reach the state, or only the potential?).** §2
  assumes the periodisation applies to $`\rho`$ itself, on the grounds that
  worlds sit on a momentum lattice and discrete $`p`$ is periodic $`y`$. The
  raised-cosine taper, by §4.4 of the algorithm spec, is applied to the
  residual *kernel* and not to the state. If some part of the implementation
  also attenuates the state's off-diagonal extent, §2's verdict — aliasing,
  not decay — is too clean and the truth is a mixture. Worth checking against
  the code rather than against the specification.
- **R-SP4 (a profile with a $`q^{-4}`$ tail).** Open item 2 of
  `../analysis/reach_energy_coupling.md`, repeated here because §5's
  third-moment table inherits it: the raised cosine gives conditional rather
  than absolute convergence, so the $`M_3`$ column's agreement is real at the
  $`10^{-4}`$ level and not beyond.
- **R-SP5 (`reach_energy_coupling.md` is not in the ladder).** It is labelled
  "step 11 of the ladder," but step 11 in
  [`../analysis/README.md`](../analysis/README.md) is
  `open_position_space.md` and the note has no entry. Presumably an omission
  when it was merged. Not fixed here, by decision.
- **R-SP6 (the vocabulary).** §7 names the collisions; nothing is renamed. The
  count of bare "particle" — excluding *world-particle*, *signed particle*,
  *N-particle* and similar — runs 19 in `position_pair_ladder.md`, 17 in
  `relational_pairing_and_carrier_lock.md`, 16 in `sea_dressed_microdynamics.md`
  and 27 in `phase_resonance_microdynamics.md`, all single-body notes where it
  cannot mean a body. The multi-body and 4-D documents score higher still but
  there the word is earned. If a rename is ever undertaken, *leg* for the end
  of a ket–bra pair is the term to standardise on, and the two senses of
  *pair* are the more urgent problem.

---

## 9. Sources

- Born, M.; von Kármán, T. — *Über Schwingungen in Raumgittern*,
  Phys. Z. **13** (1912) 297–309. The boundary condition of §3.1, and the
  original of the analogy.
- Wigner, E. — *On the quantum correction for thermodynamic equilibrium*,
  Phys. Rev. **40** (1932) 749–759. The transform of §1.
- Moyal, J. E. — *Quantum mechanics as a statistical theory*,
  Proc. Camb. Phil. Soc. **45** (1949) 99–124. The series whose radius of
  convergence is ceiling 2.
- Paley, R. E. A. C.; Wiener, N. — *Fourier Transforms in the Complex Domain*,
  AMS Colloquium Publications **19** (1934). Why finite support gives no
  lattice, §3.
- The analysis notes cited in §0, which contain all the theorems; this
  document contains none of its own except the reading in §5.
