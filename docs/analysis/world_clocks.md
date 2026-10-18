# Worlds, particles and clocks: where a reading lives

> Is the reader — its lever, chords and bank of integrators — a property of a particle, or a relation between worlds that a particle inherits? With two particles and a pair potential the question has an answer. **Definition WC0** fixes the words, for `N` particles in three dimensions: a *world* is a signed point of the whole configuration's phase space `R^6N` with one clock; a *world-particle* is one particle as it appears in one world, with no sign, clock or bank of its own; a *particle* is the identity that recurs in every world, whose only state is the marginal of the signed ensemble; a *positon* or *negaton* is a world of sign +1 or -1, and an *unpaired world* is one that is not half of an aligned pair. **Proposition WC6** gives a world's history under (S′) with persistent identity: worlds and their signs are conserved, every world-particle moves on Newtonian arcs while its world is unpaired and on inertial arcs while it is dark, momenta kink only at events and only along a bond, the parent of an event never kinks, and all worlds share one time. Two things are called a clock: the phase clock of P1, which winds at the Lagrangian rate and is a gauge except in misalignments across chords, and Q11's counting clock, an integrate-and-fire counter; they are linked through chords, so the counter is the accumulator of a reading, located where a chord's midpoint is. **Proposition WC1**: the phase clock belongs to the world — world-particle clocks count an interaction twice or each carry a share of it, and clocks split along the Lagrangian's sectors (centre of mass, relative motion), not along particles. Unlike Poirier's labelled trajectories or Hall, Deckert and Wiseman's sorted worlds, WPMW asks of the set of worlds only midpoints and difference vectors, plus a phase per world; the VLBI analogy (stations with their own clocks, information only on baselines, a correlator per pair, closure phases) is set out term by term. **Theorem WC2**: reading the joint chord (W) and reading each world-particle's own chord with its clock carrying the whole interaction (B1) agree only at linear order; exactly and to all orders, B1's relative motion is the QLE at `hbar/2`, and B1 adds a coupling `3 d_P^2 d_pr` between centre of mass and relative motion that QM lacks; splitting the interaction between the clocks (B2) halves both and reads half the force. **Corollary WC3**, measured on a `(P, r, p_r)` mesh from a product state QM keeps stationary and separable: W stays separable to `7e-15` while B1's defect is 0.2 by `t = 0.5`; B1 follows the exact `hbar/2` von Neumann reference within the mesh's error and sits 0.17 to 0.18 from QM at every reach; in the Eckart pair collision B1 transmits 0.517 against W's 0.549 at reach `8 pi`, exact 0.5120 and 0.5441. **Proposition WC4**: for pair potentials a world carries one bank per bond, `N(N - 1)/2` for `N` particles like the baselines of an array, and none per world-particle; Q11's bank per unpaired world is the one-particle case. **Proposition WC5**, ensemble runs with integrate-and-fire banks on Kaggle: the world bank never moves a daughter off its parent's centre-of-mass momentum, the world-particle banks always do; each realises its own generator up to the contact sink's bias (B1's covariance of `P^2` with `v` is -0.472 ± 0.003 against -0.501 on the mesh at `t = 0.25`, W's 0.002 ± 0.001 against 0); at `nu = 1e5` the world-particle banks carry 2.8 times the worlds and 36 times the unsigned `<P^2>` by `t = 1`, heating the ensemble's centre of mass while the signed marginal stays exact. A reader is a world, reading the chords of its own configuration along the directions in which its potential has terms. **Addendum (October 2026).** **Definition WC7** states the ontology as worlds alone — each carrying its point, sign, clock and banks, with world-particles, centre of mass and relative coordinates derived from the point — and a list of relations between worlds that carries all the quantum dynamics: chords, readings, emission, absorption, pair collisions, dark catalysis and admissibility. **Proposition WC8**: the generator reads only world-clock misalignments across chords, so a world-particle phase is either redundant (separable potentials) or wrong (shared terms, where WC1(d) shows no shares give both the action and the forces). **Corollary WC8.1**: a chord that moves only some particles reads only the terms of the potential that act on them, so locality comes from the geometry of chords, not from phases attached to particles. Identical particles and spin are left as relations between worlds still to be specified (WC-SP7, WC-SP8). **Second addendum (October 2026): three and four particles.** Particles in a harmonic trap are quenched by a pair interaction. **Proposition WC9**: the compensated residual is a sum over bonds of one two-variable function, so every bond's bank reads one kernel table, and a quadratic interaction (the Moshinsky model) entangles with no events at all; worlds reproduce its exact purity to `3e-3` for three and four particles. For a Gaussian pair interaction, the bond-sum generator on a 4D Jacobi mesh converges to the exact three-particle Schrödinger dynamics as `1/y_max^2` (purity error at `t = 3`: -0.117, -0.039, -0.009 at reach `2 pi`, `4 pi`, `8 pi`). The world ensemble follows it early and falls below it late. **Proposition WC10**: the contact sink keeps each cell's net weight but not the shape of `W` inside it, so it removes interference finer than a cell and decoheres the ensemble. The global purity at a fixed resolution falls below the mesh's, and the one-body purity falls toward classical streaming. The deficit does not depend on the time step or on how often the sink acts, does not vanish as `nu` grows at a fixed cell, and in six relative dimensions (four particles) leaves only 15 per cent of the quantum correction at `t = 1`. The note compares the results with the exact Moshinsky solution, the grid TDSE and MCTDH(B), and with the signed-particle Wigner methods of Sellier and Dimov and of Shao and Xiong, whose histogram annihilation has the same limit. It closes with what an `N`-particle algorithm specification needs; population control is the part still open (WC-SP9 to WC-SP12).

*Ladder abstract — see the [full list](README.md#the-ladder).*

**Status.** Analysis note, step 25 of the ladder. Companion demo:
`src/demo_world_clocks.py` (Parts A, K, H, M, R, T, C, F); for the second
addendum (§8.2), `src/demo_three_worlds.py`. Prompted by a
question of October 2026: the reader, its integrators and its chords are a
lot of ontology to ask of a *particle*; are they better understood as a
relation between *worlds*, which a particle inherits from its world? The
note answers for the case the ladder can test, two particles with a pair
potential, and proposes the vocabulary. Theorem prefix **WC**, the first
two-letter prefix: every single letter is taken (INDEX §2). Derivation and
implementation were developed jointly with Claude (Anthropic), October 2026.
Every quantitative claim below is an output of the demo (§10).

---

## 0. What this note asks, settles and leaves open

Up to step 24 the reader, its lever, its chords and its bank of integrators
were described for one degree of freedom, where a world and a particle are
the same thing and nothing forces a choice of words. Step 21 took the world
form to two particles, but in centre-of-mass and relative coordinates and
with Poisson events, so it never had to say who carries a clock or an
integrator. Step 23's Proposition Q11 gave every unpaired world ("every free body",
in the word of the time) a bank of integrators, and Q-SP14 asked whether the bank could instead be the sea's
own state. This note asks the question behind both: in a world of several
particles, which object owns a clock, which reads a chord, and which carries
a bank?

**Settles.**

- **Vocabulary** (§1, Definition WC0). A *world* is a signed point of the
  whole configuration's phase space with a clock; a *world-particle* is one
  particle in one world, carrying its world's sign and clock but none of its
  own; a *particle* is the identity that recurs in all worlds. A *positon*
  is a positon world, a *negaton* a negaton world, and the ledger's *body*
  becomes an *unpaired world*.
- **Trajectories** (§1.1, Proposition WC6). Under (S′) with persistent
  identity every world-particle's trajectory is continuous and piecewise
  classical: Newtonian while its world is unpaired, inertial while its
  world is dark, kinked only at events and only along a bond.
- **What a clock is** (§3). Two different things carry the name: the
  *phase clock* $`\theta`$ of P1, which winds at the Lagrangian rate and
  whose absolute value is a gauge, and the *counting clock* $`C_q`$ of Q11,
  an integrate-and-fire counter. They are linked by a chain that runs
  through chords, so the second is not a property of whatever carries the
  first.
- **The phase clock belongs to the world** (Proposition WC1). For an
  interacting pair, world-particle clocks either count the interaction twice or
  each wind with only a share of it; in neither case does a world-particle's own
  misalignment wind at the world's rate. Clocks split along the
  Lagrangian's additive sectors (centre of mass and relative motion for a
  pair potential), not along particles.
- **Readings of a pair potential** (Theorem WC2). Reading the joint chord
  (W) and reading each world-particle's own chord with its clock carrying the
  whole interaction (B1) agree only at linear order, which the compensated
  split removes anyway. Exactly and to all orders, B1's relative motion is
  the quantum Liouville equation **at $`\hbar/2`$**, and B1 adds a coupling
  $`3\partial_P^2\partial_{p_r}`$ between the centre of mass and the relative
  motion that quantum mechanics does not have. Splitting the interaction
  between the two world-particle clocks (B2) halves both and leaves the force
  half-read.
- **What that costs** (Corollary WC3, measured). From a product state that
  exact QM keeps stationary and separable, W stays separable to
  $`7\times10^{-15}`$ on the mesh while B1 entangles the centre of mass with
  the relative motion (defect $`0.2`$ by $`t = 0.5`$), and B1's relative
  motion follows the exact $`\hbar/2`$ von Neumann reference, not QM. In the
  Eckart pair collision B1 transmits $`0.517`$ against W's $`0.549`$ at the
  largest reach; the exact values are $`0.5120`$ at $`\hbar/2`$ and
  $`0.5441`$ at $`\hbar`$.
- **Banks live on bonds** (Proposition WC4). For $`N`$ particles with pair
  potentials the world's residual is a sum of exchange channels, one family
  per interacting pair, so a world carries one bank per *bond* —
  $`N(N-1)/2`$ banks, the baselines of an $`N`$-station array — and none per
  particle.
- **Worlds suffice** (addendum, §8.1: Definition WC7, Proposition WC8,
  Corollary WC8.1). The carriers are worlds alone, with their point, sign,
  clock and banks; everything about one world is derived from its point;
  and the quantum content is a list of relations between worlds — chords,
  readings, emission, absorption, pair collisions, dark catalysis and
  admissibility. The generator reads only world-clock misalignments across
  chords; a world-particle phase is redundant (separable potentials) or
  wrong (shared terms; WC1(d) adds the action-force no-go); and locality
  follows from the geometry of chords.
- **What the banks cost in worlds** (Proposition WC5, measured). With
  Q11's integrators driving the events, the world bank keeps every daughter
  on its parent's centre-of-mass momentum, as step 21's M1 requires; the
  world-particle banks move every daughter off it, and the ensemble heats in $`P`$
  while the signed $`\langle P^2\rangle`$ stays right in the mean.
  Integrate-and-fire banks realise each reading's own generator up to the
  contact sink's bias — B1's spurious covariance of $`P^2`$ with the
  relative potential is $`-0.472 \pm 0.003`$ in the world ensemble against
  $`-0.501`$ on the mesh at $`t = 0.25`$, W's $`0.002 \pm 0.001`$ against 0 —
  and at $`\nu = 10^5`$ the world-particle banks carry 2.8 times the worlds and 36
  times the unsigned $`\langle P^2\rangle`$ of the world bank by $`t = 1`$.

**Inherits.** Postulates (E), (S), (D) of step 17 with the compensated
symbol of step 14; P1 and P2 of step 4 and Lemma 0; the misalignment of
step 5; step 21's Proposition M1 (the pair kernel moves relative momentum
only) and its Pöschl–Teller and Eckart references; step 22's chord reading
(Theorems L1, L8) and step 23's Propositions Q5, Q11 and Q12.

**Corrects.** Nothing in an earlier note is withdrawn. Two readings are
sharpened and one word is retired. Q11's "every unpaired world carries a bank"
is the one-particle case of "every world carries one bank per bond"
(§7.2). The repository's *world-particle* was defined as a world (glossary,
ORIENTATION §2.3), which is right only for one particle; §1 redefines it.
The word *body* had two senses, an unpaired world in the ledger and a
physical particle in steps 21 and 24; it is retired in favour of *unpaired
world* and *particle* (WC-SP1, answered). This amendment, of October 2026,
also generalises the definitions from a line to three dimensions; the test
problems of §§6–7 stay one-dimensional.

**Leaves open.** WC-SP2 to WC-SP8 (§11), among them whether a live reading of
the relative-direction sea can feed a bond's bank (Q12 in four dimensions)
and how a three-body potential's bank would be indexed.

---

## 1. Worlds, world-particles and particles

**Definition WC0.** Let $`N`$ particles move in three-dimensional position
space, with masses $`m_i`$, $`i = 1, \dots, N`$. A **world** $`w`$ is one
complete configuration with its momenta, a sign and a clock: a point

```math
(\mathbf x^w, \mathbf p^w) = (\mathbf x_1^w, \dots, \mathbf x_N^w,\ \mathbf p_1^w, \dots, \mathbf p_N^w) \in \mathbb R^{6N},
\qquad \mathbf x_i^w, \mathbf p_i^w \in \mathbb R^3,
```

with a sign $`\varepsilon^w = \pm1`$ and a phase $`\theta^w`$ that winds by
P1, $`\dot\theta^w = \mathcal L(\mathbf x^w, \mathbf p^w)/\hbar`$ with
$`\mathcal L = \sum_i \lvert\mathbf p_i\rvert^2/2m_i - V(\mathbf x)`$. Postulate
(E) says the ensemble of worlds is a signed counting measure on
$`\mathbb R^{6N}`$. The sign, the clock and the banks of integrators (§5.4:
one per bond) belong to the world. Then:

- A **positon** is a world with $`\varepsilon = +1`$ and a **negaton** one
  with $`\varepsilon = -1`$: *positon* is short for *positon world*. Neither
  is a particle, neither is a hole, and neither carries a charge of its own.
- A **world-particle** $`(w, i)`$ is particle $`i`$ as it appears in world
  $`w`$: the point $`(\mathbf x_i^w, \mathbf p_i^w) \in \mathbb R^6`$. It has
  no sign, clock or bank of its own; the $`N`$ world-particles of a world
  share its sign and its clock.
- A **particle** $`i`$ is the identity that recurs in every world, with its
  mass and its couplings. It has no state of its own; what can be said of it
  is the marginal of the signed ensemble, below.
- An **aligned pair** is a positon and a negaton whose clocks agree (step
  22, L0); it is *dark*, contributing nothing to the observable. An
  **unpaired world** is a world that is not half of an aligned pair. The
  ledger counts unpaired positons and negatons, $`n_+`$ and $`n_-`$, per
  joint cell. Earlier versions of the ledger notes called an unpaired world a *body*,
  and steps 21 and 24 also used *body* for a particle; the word is retired.

The state of particle $`i`$ is its reduced Wigner function, read at the
resolution of a cell:

```math
W_i(\mathbf x, \mathbf p) = \sum_w \varepsilon^w\thinspace\delta(\mathbf x - \mathbf x_i^w)\thinspace\delta(\mathbf p - \mathbf p_i^w).
```

For $`N = 1`$ a world has exactly one world-particle, which is why the
one-particle ladder could use *world* and *world-particle* for one object
and why nothing in it changes. Every statement below about "a
world-particle's clock" or "a world-particle bank" is about a *proposal* —
an ontology in which world-particles carry them — and the note's results
say why the proposal fails. The test problems of §§6–7 put two particles on
a line; the definitions and Propositions WC1, WC4 and WC6 hold in three
dimensions as stated, and WC2 and WC3 hold with the chords, transfers and
derivatives read as vectors (§5.2).

### 1.1 Proposition WC6: a world's history

From here on the ladder takes postulate (S′) as amended in step 24 —
unpaired worlds feel the full force, aligned pairs are force-blind and
change momentum only in pair–pair collisions (Definition (C)), and every
clock winds by P1 — and reads the events with **persistent identity**
(step 20): no world is created or destroyed; an emission ionises an
aligned pair and an absorption joins two unpaired worlds into an aligned
pair, each a kink in worlds that already exist.

**Proposition WC6.** *Under (E), (S′) and (D), with the realisations of the
specification's §5.3 ([`compensated_liouville_algorithm.md`](../algorithm/compensated_liouville_algorithm.md#53-two-realisations-and-two-that-are-excluded))
read with persistent identity:*

*(a) The number of worlds and the sign of each world are conserved. The
ledger's counts satisfy $`\Delta N + 2\Delta S = 0`$, with $`N = n_+ + n_-`$ the unpaired worlds and
$`S`$ the aligned pairs.*

*(b) Between events, every world-particle of world $`w`$ obeys*

```math
\dot{\mathbf x}_i^w = \frac{\mathbf p_i^w}{m_i},\qquad
\dot{\mathbf p}_i^w =
\begin{cases}
-\nabla_i V(\mathbf x^w) & w \text{ unpaired},\\[2pt]
0 & w \text{ half of an aligned pair},
\end{cases}
\qquad
\dot\theta^w = \frac{\mathcal L(\mathbf x^w, \mathbf p^w)}{\hbar}.
```

*(c) At an event in channel $`q`$ of a term of the potential — a bond
$`(i, j)`$ for a pair potential — the positions of every world are
continuous, and the momenta of the worlds the event moves jump along that
term only: for a bond, $`\Delta\mathbf p_i = \pm\pmb\xi_q`$,
$`\Delta\mathbf p_j = \mp\pmb\xi_q`$ and $`\Delta\mathbf p_k = 0`$ for
$`k \ne i, j`$, so $`\sum_i\mathbf p_i`$ is unchanged. The worlds moved are
the two members of an aligned pair in the parent's cell and momentum row
(emission, each kinked by $`\pm\pmb\xi_q`$ to an unpaired world) or two
unpaired worlds at the parent's momentum $`\pm\pmb\xi_q`$ (absorption,
each kinked back to the parent's row, where they become an aligned pair).
A pair–pair collision moves both members of a pair by the same jump.*

*(d) The parent of an event — the world whose bank fires, the reader —
does not move at its own events. A world that only ever parents keeps one
unbroken Newtonian trajectory.*

*(e) All worlds share one Newtonian time $`t`$, and events are
instantaneous. A world's history is a continuous curve in configuration
space whose momenta are piecewise $`C^1`$: Newtonian arcs while it is
unpaired, inertial arcs while it is dark, joined by kinks. Its clock
measures action along that history, not time.*

*Proof.* (b) is (S′) with P1. For (a) and (c), the specification's §5.3
lists the realisations of an event at the parent's momentum $`\mathbf p`$:
emission ionises the aligned pair co-located with the parent in position
and momentum (§5.2), sending its positon to one daughter row
$`\mathbf p \pm \pmb\xi`$ and its negaton to the other; absorption takes
an unpaired world of each sign from the two daughter rows and binds them
at $`\mathbf p`$. Each row moves an existing world of fixed sign, so signs
and numbers are conserved, and Theorem S7's $`\Delta N = 2(1 - 2f)n_{\rm ev}`$
and $`\Delta S = (2f - 1)n_{\rm ev}`$ sum to $`\Delta N + 2\Delta S = 0`$.
Definition (C) moves whole pairs. For a pair potential the transfer
$`\pmb\xi`$ lies along the bond's exchange direction (Proposition M1,
pair by pair in WC4), which gives the form of the jump. No row of the table
moves the parent, which is (d). Positions are untouched by every event,
which with (b) gives (e). $`\square`$

Two remarks. Continuity of position is exact for the moved worlds when
absorption is realised at crossings (Proposition Q9); on a lattice it holds
to within the matching cell. And the ledger's sea *counter* — the
specification's per-cell scalar — is the mean-field record of the dark
worlds, not a reservoir from which worlds appear; Part C of this note's
demo used that counter, placing an emission's two worlds at the parent's
position and removing opposite-sign worlds in one cell as a contact sink,
which under persistent identity is a contact join (an absorption).

## 2. What other many-worlds models ask a world to know

Every interacting-worlds model has to say how one world is affected by the
others. The answers differ in what structure the set of worlds must carry.

- **Poirier** (2010; Schiff and Poirier 2012) labels trajectories by a
  continuous parameter $`C`$, takes $`\mathbf x(C, t)`$ to be a
  diffeomorphism, and writes the quantum force through the Jacobian
  $`\partial\mathbf x/\partial C`$ and its derivatives. The worlds form a
  smooth manifold: what one world "knows" about its neighbours is the
  derivative of the labelling, which needs differentiability in the label,
  not an order. In one dimension the label can be taken as the cumulative
  fraction, $`\partial x/\partial C = 1/\rho`$.
- **Hall, Deckert and Wiseman** (2014) take a finite set of worlds whose
  index carries no meaning. Their interworld potential in one dimension
  uses the worlds sorted by position, so a world knows its two neighbours
  in that order. Sorting is a one-dimensional accident: in more than one
  dimension there is no total order compatible with neighbourhood, and a
  neighbour structure or a density estimate has to be supplied separately.
  "Well-ordered" is not what is needed; every finite set can be
  well-ordered. What is needed is an order or a neighbourhood that the
  geometry respects.
- **Sebens** (2015) proposes a continuum of interacting worlds without an
  explicit model of the interaction.

**WPMW needs neither a smooth labelling nor an order.** What one world
reads about another is the misalignment of their clocks across a *chord*,
the pair of configurations $`\mathbf X \pm \mathbf Y`$, read at the
midpoint $`\mathbf X`$. That requires only the affine structure every
configuration space already has (midpoints and difference vectors) and one
phase per world. The price is paid elsewhere: the relation is non-local out
to the reach $`y_{\max}`$ (Theorem G5 makes the reach a regulator), and the
reading needs a sea dense enough to supply chords (Q12, V4).

## 3. What is a clock

### 3.1 The phase clock

P1 gives every world a phase that winds along its worldline at the
Lagrangian rate,

```math
\dot\theta = \frac{\mathcal L}{\hbar},\qquad
\mathcal L = \sum_{i=1}^N \frac{p_i^2}{2m_i} - V(\mathbf x),
```

so that $`\hbar\theta`$ is Hamilton's principal function accumulated along
the world's own history. Lemma 0 of step 4 shows the rate is forced: P2
extrapolates the phase at the world's own wavevector, and only the
Lagrangian rate makes that extrapolation time-consistent, up to a constant.

The absolute value of $`\theta`$ is not physical. Under $`V \to V + c`$
every world winds faster by $`-c/\hbar`$; under a change of the zero of time
or of the energy every clock moves together. What survives is the
**misalignment** of two clocks transported to a common point by the lever
$`P`$ (step 5's $`\mu`$; P2; Lemma 0),

```math
\mu^{(P)}_{ij} = \theta_i - \theta_j + \frac{P\cdot(\mathbf x_j - \mathbf x_i)}{\hbar},
```

and Q11(d) showed that a reading built from single clocks against a plane
wave depends on the zero of $`V`$ while one built from chords does not.
Part A6 verifies that $`V \to V + c`$ leaves the world's chord symbol
unchanged.

### 3.2 The counting clock

Proposition Q11 introduced a second object also called a clock: an
integrator $`C_q`$ per channel, advanced by the channel's rate along the
worldline, $`\dot C_q = K_q`$, firing each time it crosses an integer. It
keeps time in the sense a Geiger counter or a sigma–delta modulator does:
it counts. It has no phase that winds at a physical rate, and its only
hidden variable is the fractional part of its starting value. Step 23
showed that it realises the generator's rate exactly in mean with
sub-Poissonian counts.

### 3.3 How the two are linked

The chain from one to the other runs through chords, and every link but
the first lives on a *pair* of worlds or at a reader, not on a world
(schematically; the normalisations are those of L1 and Q5):

```math
\underbrace{\theta}_{\text{per world}}
\;\longrightarrow\;
\underbrace{\mu_{ij}}_{\text{per chord (pair of worlds)}}
\;\longrightarrow\;
\underbrace{z_q = \sum_{\text{chords}} e^{\,i\mu_{ij}}}_{\text{per reader, per channel}}
\;\longrightarrow\;
\underbrace{K_q \propto \tfrac{d}{dt}\,\mathrm{Re}\,z_q}_{\text{rate}}
\;\longrightarrow\;
\underbrace{C_q = \int K_q\,dt}_{\text{count}} .
```

The third and fourth links are Theorem L1 and Proposition Q5(b): the sea
computes the kernel by interference, and at re-lock the channel-$`q`$ rate
is the rate of change of the row's interference sum read against the
daughter row. By Q11's own reading, $`C_q`$ is that interference sum
counted in units of the sea's density.

So the counting clock is not a second property of whatever carries the
phase clock. It is the accumulator of a **reading**, and a reading is
located where a chord's midpoint is. Asking which object owns the bank is
therefore asking which chords are read — the question of §5.

### 3.4 Proposition WC1: the phase clock belongs to the world

**Proposition WC1.** *Let $`\mathcal L = \sum_i T_i - V(\mathbf x)`$ with
$`T_i = p_i^2/2m_i`$.*

*(a) If $`V = \sum_i V_i(x_i)`$, the world clock is a sum of world-particle
clocks, $`\theta = \sum_i\theta_i`$ with $`\dot\theta_i = (T_i - V_i)/\hbar`$,
each a function of one world-particle's worldline.*

*(b) For $`V = v(x_1 - x_2)`$, assign world-particle clocks
$`\dot\theta_i = (T_i - \alpha_i v)/\hbar`$. Their sum is the world clock if
and only if $`\alpha_1 + \alpha_2 = 1`$; with $`\alpha_1 = \alpha_2 = 1`$
(each world-particle treating the other as an external field) the sum counts the
interaction Lagrangian twice. In either case the misalignment of world-particle
$`i`$'s own clocks across a chord $`\mathbf Y`$ winds at
$`-\alpha_i\thinspace U(\mathbf Y)/\hbar`$, not at the world's
$`-U(\mathbf Y)/\hbar`$, where $`U(\mathbf Y) = v(r + \eta) - v(r - \eta)`$
and $`\eta = y_1 - y_2`$.*

*(c) Clocks split along the Lagrangian's additive sectors: for a pair
potential $`\mathcal L = P^2/2M + p_r^2/2\mu - v(r)`$, so
$`\theta = \theta_{\rm cm} + \theta_{\rm rel}`$ with the centre-of-mass clock
free. The sectors are normal coordinates of the world, not particles.*

*(d) (Addendum, October 2026.) No choice of shares gives world-particle
clocks that both sum to the world's action and each read the Newtonian
force. Across its own chord, $`\mathbf x_i \to \mathbf x_i \pm \mathbf y_i`$
with every other coordinate and every momentum fixed, world-particle
$`i`$'s clocks read the force $`-\alpha_i\nabla_i v`$; Newton needs
$`\alpha_1 = \alpha_2 = 1`$, the action needs $`\alpha_1 + \alpha_2 = 1`$.*

*Proof.* (a) and (c) are additivity of $`\mathcal L`$. For (b), Part A5
computes $`\mathcal L_1 + \mathcal L_2 - \mathcal L = (1 - \alpha_1 - \alpha_2)\thinspace v`$,
which is $`-v`$ at $`\alpha_i = 1`$ and $`0`$ when the shares sum to one.
For two worlds $`\mathbf X \pm \mathbf Y`$ on one momentum row the kinetic
windings cancel, so their clocks' misalignment winds at
$`-[V(\mathbf X + \mathbf Y) - V(\mathbf X - \mathbf Y)]/\hbar = -U/\hbar`$
(the one-particle statement is step 23's Proposition Q2: a same-row pair
winds at exactly $`U`$), and world-particle $`i`$'s clocks, which carry
$`\alpha_i v`$, at $`-\alpha_i U/\hbar`$. For (d), on world-particle
$`i`$'s own chord the misalignment rate is
$`-\alpha_i[v(r + y_i) - v(r - y_i)]/\hbar`$, whose linear part
$`-2\alpha_i y_i v'(r)/\hbar`$ is the reading of the force
$`-\alpha_i v'(r)`$ (the chord-reading normalisation of step 22); the two
conditions on $`(\alpha_1, \alpha_2)`$ have no common solution (Part G1).
$`\square`$

*Where the double count is.* With $`\alpha_1 = \alpha_2 = 1`$ each
world-particle winds with the energy of the other as an external field, and
the two energies sum to twice the interaction energy — the double counting
that puts the $`\tfrac12`$ into an electrostatic energy,
$`W = \tfrac12\sum_i q_i\Phi(\mathbf x_i)`$ (Jackson 1999, §1.11). Forces
carry no $`\tfrac12`$: each particle feels the full pair force. A
world-particle clock therefore has to choose between the force and the
action, and §5.2 shows where the choice shows: across a relative chord
$`\eta`$ the B1 clocks ($`\alpha_i = 1`$) read two half-chords each
carrying the whole $`v`$, $`2[v(r + \eta/2) - v(r - \eta/2)]`$, which is the
world's reading at $`\hbar/2`$ (WC2(b)); the B2 clocks
($`\alpha_i = \tfrac12`$) read half the force (WC2(a)).

The interaction energy is a property of the configuration, not of either
particle, and the clock that integrates it inherits that. A particle can
carry a share of the world's clock only by convention, and §5 shows that no
convention makes a world-particle's own reading equal to the world's.

## 4. The VLBI analogy

Very-long-baseline interferometry builds an image of a radio source from
stations that never exchange signals while they observe. Each station
records against its own hydrogen-maser clock; the recordings are brought
together afterwards in a correlator, which applies a delay and fringe-rate
model to align them and accumulates, for each **pair** of stations, the
complex visibility on that baseline. By the van Cittert–Zernike theorem
(van Cittert 1934; Zernike 1938) the visibilities are samples of the
Fourier transform of the sky's brightness; the image is their inverse
transform. Station clock errors are unobservable in the image because they
cancel in *closure phases*, the sum of visibility phases round a triangle of
baselines (Jennison 1958). Thompson, Moran and Swenson (2017) is the
standard account.

| VLBI | WPMW | where in the ladder |
|---|---|---|
| station | world | (E) |
| station's maser clock | phase clock $`\theta`$ | P1, Lemma 0 |
| clock offset, unobservable | absolute $`\theta`$, a gauge | Q11(d), Part A6 |
| baseline (pair of stations) | chord (pair of worlds $`\mathbf X \pm \mathbf Y`$) | L1 |
| delay and fringe-rate model | lever $`P`$ carrying clocks to a common point | P2, L8 |
| visibility on a baseline | misalignment across a chord; $`\rho(\mathbf X + \mathbf Y, \mathbf X - \mathbf Y)`$ | step 5, step 12 |
| van Cittert–Zernike: image = FT of visibilities | Wigner function = FT over chords of $`\rho`$ | — |
| correlator | reader at the chord's midpoint | L8, glossary |
| accumulation over the integration time | integrator $`C_q`$ | Q11 |
| one fringe passing | one firing | Q11 |
| $`(u, v)`$ coverage, sampling density | reach and the sea's density | G5, Q12, V4 |
| closure phase | gauge-invariant combinations of misalignments | Part A6 |
| $`N`$ stations, $`N(N-1)/2`$ baselines | within one world: $`N`$ particles, $`N(N-1)/2`$ bonds | WC4 |

![VLBI against WPMW](https://raw.githubusercontent.com/billpage/wpmw/output/figures/world_clocks_vlbi.png)

*Left: four stations, each with its own clock, give six baselines whose
visibilities the correlator accumulates. Right: two worlds of a
two-particle system joined by a chord in the configuration plane; the
reader sits at the midpoint. For a pair potential only the chord's
extent along $`y_1 - y_2`$ is read; the centre-of-mass direction is dark.
A world-particle bank reads only the chord's extent along its own axis.*

Two lessons carry over directly. The first is that the correlator is the
right picture of a reader: it belongs to neither station and exists only
for the pair. The second is that a per-station quantity can never replace
a per-baseline one: an array's information is in its baselines, and a
station that correlated only against itself would see nothing.

Where the analogy stops. The correlator is outside the array and passive;
a reader's firings change the population. VLBI stations are fixed and
their source is static, while worlds move and the "source" is the state
being evolved. And VLBI has no sign: its visibilities sample a positive
brightness, while WPMW's chords sample a signed $`W`$.

## 5. Three readings of a pair potential

### 5.1 Setting

Two particles of unit mass on a line, $`\hbar = 1`$ except where $`\hbar`$
is kept to show its role, interacting through $`V(x_1, x_2) = v(x_1 - x_2)`$.
With step 21's coordinates

```math
X = \tfrac12(x_1 + x_2),\quad P = p_1 + p_2,\qquad
r = x_1 - x_2,\quad p_r = \tfrac12(p_1 - p_2),
```

$`M = 2`$, $`\mu = 1/2`$. Let $`\sigma_1, \sigma_2`$ be conjugate to
$`p_1, p_2`$ and $`s_P, s_r`$ to $`P, p_r`$; matching
$`p_1\sigma_1 + p_2\sigma_2 = P s_P + p_r s_r`$ gives

```math
\sigma_1 = s_P + \tfrac12 s_r,\qquad \sigma_2 = s_P - \tfrac12 s_r .
```

The compensated residual symbol of a chord $`\mathbf y = \hbar\pmb\sigma/2`$
is $`(i/\hbar)\thinspace u(\mathbf y)`$ with $`u`$ the potential difference less its
linear part. Three readings:

```math
\begin{aligned}
u_W &= v(r + y_1 - y_2) - v(r - y_1 + y_2) - 2(y_1 - y_2)\thinspace v'(r)
&&\text{the world reads the joint chord;}\\
u_{B1} &= \bigl[v(r + y_1) - v(r - y_1) - 2y_1 v'\bigr] + \bigl[v(r - y_2) - v(r + y_2) + 2y_2 v'\bigr]
&&\text{each world-particle its own, whole } v;\\
u_{B2} &= \tfrac12\thinspace u_{B1}
&&\text{each world-particle its own, } v \text{ split.}
\end{aligned}
```

B1 is what a world-particle reads if it treats the other as an
external field ($`\alpha_i = 1`$ in WC1); B2 is the bookkeeping that makes
the world-particle clocks sum to the world's ($`\alpha_i = \tfrac12`$). Write the
one-dimensional residual symbol as

```math
m(\sigma) = \frac{i}{\hbar}\Bigl[v\bigl(r + \tfrac{\hbar\sigma}{2}\bigr) - v\bigl(r - \tfrac{\hbar\sigma}{2}\bigr) - \hbar\sigma\thinspace v'(r)\Bigr]
```

(odd in $`\sigma`$). Then the symbols in $`(s_P, s_r)`$ are

```math
S_W = m(s_r),\qquad
S_{B1} = m\bigl(s_P + \tfrac12 s_r\bigr) - m\bigl(s_P - \tfrac12 s_r\bigr),\qquad
S_{B2} = \tfrac12\thinspace S_{B1}.
```

$`S_W`$ does not depend on $`s_P`$: that is M1, and it means every world
event is an exchange $`(p_1, p_2) \to (p_1 \pm \xi, p_2 \mp \xi)`$. A B1
event shifts one particle's momentum, $`p_i \to p_i \pm \xi`$, which moves
$`P`$ by $`\pm\xi`$ and $`p_r`$ by $`\pm\xi/2`$.

### 5.2 Theorem WC2: what a world-particle reads

**Theorem WC2.** *For any smooth $`v`$, with $`a_k = v^{(k)}(r)`$:*

*(a) The readings agree at linear order and differ from the third:*

```math
u_W - u_{B1} = -a_3\thinspace y_1y_2(y_1 - y_2)
- \frac{a_5}{12}\thinspace y_1y_2(y_1 - y_2)(y_1^2 - y_1y_2 + y_2^2)
- \frac{a_7}{360}\thinspace y_1y_2(y_1 - y_2)(y_1^2 - y_1y_2 + y_2^2)^2 - \cdots
```

*Uncompensated, B2's linear term is $`a_1(y_1 - y_2)`$, half the world's
$`2a_1(y_1 - y_2)`$: B2 reads half the force.*

*(b) B1's relative motion is the world's at $`\hbar/2`$, exactly:*

```math
S_{B1}(s_P = 0,\ s_r) = m_{\hbar/2}(s_r)\quad\text{for every } s_r,
```

*where $`m_{\hbar/2}`$ is $`m`$ with $`\hbar`$ replaced by $`\hbar/2`$. The
marginal $`W_{\rm rel}(r, p_r) = \int W\thinspace dP`$ is the $`s_P = 0`$
section, on which the symbol acts pointwise, so it evolves autonomously
under the quantum Liouville equation with Planck's constant halved in the
residual and the classical part unchanged.*

*(c) B1 leaves the centre-of-mass marginal alone: $`S_{B1}(s_P, s_r = 0) = 0`$.*

*(d) In derivatives, the cubic Moyal operators are*

```math
(\partial_{p_1} - \partial_{p_2})^3 = \partial_{p_r}^3
\quad\text{(W)},\qquad
\partial_{p_1}^3 - \partial_{p_2}^3 = 3\thinspace\partial_P^2\partial_{p_r} + \tfrac14\thinspace\partial_{p_r}^3
\quad\text{(B1)},
```

*so B1 reads a quarter of the relative motion's cubic correction (B2 an
eighth) and adds a term $`3\thinspace\partial_P^2\partial_{p_r}`$, with the
cubic Moyal coefficient $`\hbar^2 v'''(r)/24`$, that couples the centre of mass to the relative
motion.*

*Proof.* (a) Expand $`v(r + z) = \sum_k a_k z^k/k!`$; Part A1 collects the
difference by total degree (SymPy, exact). (b) At $`s_P = 0`$,
$`S_{B1} = m(s_r/2) - m(-s_r/2) = 2m(s_r/2)`$ because $`m`$ is odd, and

```math
2\thinspace m\bigl(\tfrac{s_r}{2}\bigr)
= \frac{2i}{\hbar}\Bigl[v\bigl(r + \tfrac{\hbar s_r}{4}\bigr) - v\bigl(r - \tfrac{\hbar s_r}{4}\bigr) - \tfrac{\hbar s_r}{2}\thinspace v'\Bigr]
= \frac{i}{\hbar'}\Bigl[v\bigl(r + \tfrac{\hbar' s_r}{2}\bigr) - v\bigl(r - \tfrac{\hbar' s_r}{2}\bigr) - \hbar' s_r\thinspace v'\Bigr]
```

with $`\hbar' = \hbar/2`$ (Part A2, residual $`0`$). (c) At $`s_r = 0`$ the
two terms of $`S_{B1}`$ coincide (Part A3). (d) $`\partial_{p_{1,2}} = \partial_P \pm \tfrac12\partial_{p_r}`$;
Part A4 gives the cubic coefficients $`i a_3\hbar^2/24`$, $`i a_3\hbar^2/96`$
and $`i a_3\hbar^2/192`$ for W, B1 and B2. $`\square`$

**Why $`\hbar/2`$.** A B1 event delivers its momentum $`\xi`$ to one
particle, which moves the relative momentum by only $`\xi/2`$, while the
rate pattern it fires at is the one-particle kernel at the full chord. The
relative motion therefore receives the residual's Fourier content at half
the momentum per unit of chord phase: the quantum correction is read at
half of Planck's constant. A world-particle's reader is not wrong about the
chord; it is wrong about who takes the kick.

**In three dimensions.** With $`\pmb\sigma_i \in \mathbb R^3`$ and the
one-particle residual $`m(\pmb\sigma)`$ odd under $`\pmb\sigma \to -\pmb\sigma`$,
every step above holds verbatim with vectors: (b) uses only
$`S_{B1}(\mathbf s_P = 0, \mathbf s_r) = m(\mathbf s_r/2) - m(-\mathbf s_r/2) = 2m(\mathbf s_r/2)`$,
and (d) becomes the tensor identity
$`(\mathbf a + \tfrac12\mathbf b)^{\otimes3} - (\mathbf a - \tfrac12\mathbf b)^{\otimes3}`$
contracted with $`\nabla^3 v`$, whose terms odd in $`\mathbf b = \nabla_{\mathbf p_r}`$
and even in $`\mathbf a = \nabla_{\mathbf P}`$ are the centre-of-mass coupling.

### 5.3 Corollary WC3: what B1 does to a state quantum mechanics keeps fixed

**Corollary WC3.** *Let the initial state be a product
$`W_0 = g(P)\thinspace w(r, p_r)`$ with $`w`$ the Wigner function of a
relative wavefunction $`\chi`$ at $`\hbar`$, and the centre of mass free.*

*(a) Under W the state stays a product for all time, as in QM.*

*(b) Under B1 the relative marginal is the $`\hbar'`$ Wigner function of
the von Neumann evolution, at $`\hbar' = \hbar/2`$, of*

```math
\rho'(a, b) = \chi\bigl(\tfrac{3a - b}{2}\bigr)\thinspace\chi^*\bigl(\tfrac{3b - a}{2}\bigr),
```

*the operator whose $`\hbar'`$ Wigner function is $`w`$. In general
$`\rho'`$ is not stationary when $`\chi`$ is, and it need not be positive.*

*(c) Under B1 the state does not stay a product: the
$`\partial_P^2\partial_{p_r}`$ term correlates $`P`$ with the relative
motion, while the $`P`$-marginal itself is unchanged by (c) of WC2.*

*Proof.* (a) $`S_W`$ is independent of $`s_P`$ and the classical part acts
on $`(r, p_r)`$ alone, so the generator is $`1 \otimes \mathcal G_{\rm rel}`$.
(b) By WC2(b) the marginal obeys the $`\hbar'`$ Moyal equation, which is the
Wigner image of von Neumann's at $`\hbar'`$. Inverting the $`\hbar'`$ Weyl
transform of $`w`$,

```math
\rho'(\bar r + \eta, \bar r - \eta) = \int w(\bar r, p)\thinspace e^{2ip\eta/\hbar'}\thinspace dp
= \frac{1}{\pi\hbar}\int dy\thinspace\chi(\bar r + y)\chi^*(\bar r - y)\int dp\thinspace e^{-2ip(y - 2\eta)/\hbar}
= \chi(\bar r + 2\eta)\thinspace\chi^*(\bar r - 2\eta),
```

and $`\bar r \pm 2\eta = (3a - b)/2,\ (3b - a)/2`$ for $`a, b = \bar r \pm \eta`$.
(c) Measured, Part M. $`\square`$

For the Pöschl–Teller ground state $`\chi = \psi_0`$ of §6, $`\rho'`$ has
trace 1 and a spectrum in $`[-1.72\times10^{-4},\ 0.666]`$ (negative mass
$`5.9\times10^{-4}`$): read at $`\hbar/2`$ the ground state is very nearly,
but not exactly, an admissible state, and it is not stationary. For the
Gaussian relative packet of the Eckart collision, $`\rho'`$ is a positive
mixed Gaussian (least eigenvalue $`-3\times10^{-16}`$).

### 5.4 Proposition WC4: banks live on bonds

**Proposition WC4.** *For $`N`$ particles with
$`V = \sum_{i<j} v_{ij}(x_i - x_j)`$, the world's residual symbol is a sum
of pair-exchange symbols,*

```math
S_W(\pmb\sigma) = \sum_{i < j} m_{ij}(\sigma_i - \sigma_j),
```

*so every world event is an exchange of momentum between one pair, and the
world's bank factorises into one bank per interacting pair: $`N(N-1)/2`$
banks for a fully connected system, one for each edge of the interaction
graph. A bank per world-particle ($`N`$ of them) reproduces none of these beyond
linear order (WC2(a) applied to each pair). For a three-body term the bank
lives on the triple, and for a general $`V(\mathbf x)`$ on the full joint
chord, with channels indexed by joint momentum transfers.*

*Proof.* The residual is linear in $`V`$, and each pair term depends on
$`\mathbf y`$ only through $`y_i - y_j`$ (M1 applied pair by pair). Part A7
verifies the sum for $`N = 3`$. $`\square`$

A world's clock is one phase. Its banks are as many as the terms of its
potential: the clock belongs to the world, the banks to the world's bonds.
The VLBI analogy enters here at a second level. Between worlds, a baseline
is a chord (§4). Within one world, the count of banks is an array's count
of baselines: $`N`$ particles have $`N(N-1)/2`$ bonds as $`N`$ stations have
$`N(N-1)/2`$ baselines, and a two-particle world has one, as a two-station
interferometer has one. At both levels what carries the information is a
pair, never a single member.

## 6. The test problem

Step 21's bound problem with the centre of mass made free, so that the
state can be carried as $`W(P, r, p_r)`$ and $`X`$ is never read:
$`v(r) = -V_0\thinspace\mathrm{sech}^2 r`$ with $`V_0 = 6`$, ground state
$`\psi_0 = \sqrt{3/4}\thinspace\mathrm{sech}^2 r`$ with $`E_0 = -4`$,
$`\langle v\rangle = -4.8`$, $`\langle p_r^2/2\mu\rangle = 0.8`$, and
$`g(P)`$ a Gaussian of width $`\sigma_P = 1`$. Exact QM keeps
$`g\otimes W_{\psi_0}`$ stationary and separable. Mesh: $`P`$ on 64 rows of
$`\Delta p = 0.25`$, $`r \in [-10, 10)`$ on 128 points, $`p_r`$ on 128 rows of
$`0.125`$ (half a row, so that a B1 event lands on a row), reach
$`y_{\max} = 2\pi`$ with the raised-cosine horizon, Strang splitting with
$`\Delta t = 0.01`$. The classical force is applied exactly and the
residual in Fourier space, so a symbol independent of $`s_P`$ keeps
$`P`$ untouched to rounding.

The kernel (Part K): 31 exchange channels $`\xi_q = q\thinspace\Delta p`$,
$`\Gamma(r) = \sum_q\lvert K_q\rvert`$ at most $`9.39`$ and
$`\langle\Gamma\rangle_{\psi_0} = 6.11`$ per unit time. A world bank holds 31
integrators and fires at $`\Gamma`$; two world-particle banks hold 62 and fire at
$`2\Gamma`$.

## 7. Results

### 7.1 Mesh runs (Part M) and the exact references (Parts H, R)

In exact QM every row of the table is zero except the energy, which is
$`-4`$. The mesh's own error is what W shows: a relative $`L^2`$ drift of
$`2\times10^{-2}`$, the value step 21's Part R found for the same mesh.

| at $`t = 4`$ | W | B1 | B2 | classical streaming |
|---|---|---|---|---|
| relative $`W`$: $`L^2`$ distance from $`W_{\psi_0}`$ | $`1.85\times10^{-2}`$ | $`0.179`$ | $`0.129`$ | $`0.165`$ |
| separability $`\lVert W - g_P\otimes W_{\rm rel}\rVert/\lVert W\rVert`$ (max over $`t`$) | $`7.2\times10^{-15}`$ | $`0.216`$ | $`0.097`$ | $`7.2\times10^{-15}`$ |
| $`\mathrm{cov}(P^2, v)`$ | $`6\times10^{-14}`$ | $`-1.46`$ | $`-0.20`$ | $`-5\times10^{-14}`$ |
| change of the $`P`$-marginal, $`L^1`$ | $`5\times10^{-14}`$ | $`5\times10^{-14}`$ | $`5\times10^{-14}`$ | $`5\times10^{-14}`$ |
| relative energy $`\langle p_r^2/2\mu + v\rangle`$ (exact $`-4`$) | $`-4.00001`$ | $`-4.00008`$ | $`-4.00007`$ | $`-3.99993`$ |

Four checks against exact statements:

- **WC3(a).** Under W the state stays a product to $`7\times10^{-15}`$ and the
  covariance of $`P^2`$ with the relative potential is zero to rounding.
  Under B1 the separability defect reaches $`0.216`$ by $`t = 0.5`$ and the
  covariance oscillates with amplitude above 1. B2 does the same at about
  half the size.
- **WC2(c).** The $`P`$-marginal is unchanged under all three to
  $`5\times10^{-14}`$: B1's spurious coupling moves no centre-of-mass
  probability, only its correlation with the relative motion.
- **WC2(b).** B1's relative marginal equals an independent two-dimensional
  mesh run at $`\hbar/2`$ to a relative $`6.6\times10^{-15}`$ (W's equals
  the $`\hbar`$ run to $`5.5\times10^{-15}`$).
- **WC3(b).** Against the exact references, as relative $`L^2`$ distances of
  the relative Wigner function:

| $`t`$ | W mesh vs QM | B1 mesh vs exact $`\hbar/2`$ | B1 mesh vs QM | exact $`\hbar/2`$ vs QM |
|---|---|---|---|---|
| 1 | $`2.7\times10^{-2}`$ | $`1.8\times10^{-2}`$ | $`7.6\times10^{-2}`$ | $`7.3\times10^{-2}`$ |
| 2 | $`2.9\times10^{-2}`$ | $`4.0\times10^{-2}`$ | $`0.119`$ | $`0.108`$ |
| 3 | $`2.5\times10^{-2}`$ | $`4.8\times10^{-2}`$ | $`0.157`$ | $`0.145`$ |
| 4 | $`1.8\times10^{-2}`$ | $`5.1\times10^{-2}`$ | $`0.179`$ | $`0.168`$ |

B1 sits on the $`\hbar/2`$ reference within the mesh's error and moves away
from QM about as far as classical streaming does. The expectation
$`\langle v\rangle`$ is a poor witness here — it moves by only a few
hundredths in every reading, $`-4.806`$ at $`t = 4`$ for the exact
$`\hbar/2`$ evolution against $`-4.8`$ — because the energy is conserved
and the well is deep; the Wigner function itself is the witness.

**What the mesh error is (Part R).** It is the reach, not the grid:
halving $`\Delta r`$ and $`\Delta t`$ together changes the errors by at most
5 per cent, while raising the reach shrinks them.

| reach, grid | W vs QM | B1 vs exact $`\hbar/2`$ | B1 vs QM |
|---|---|---|---|
| $`2\pi`$, $`128\times128`$, $`\Delta t = 0.01`$ | $`1.85\times10^{-2}`$ | $`5.06\times10^{-2}`$ | $`0.179`$ |
| $`2\pi`$, $`256\times128`$, $`\Delta t = 0.005`$ | $`1.95\times10^{-2}`$ | $`5.10\times10^{-2}`$ | $`0.180`$ |
| $`4\pi`$, $`128\times256`$ | $`8.5\times10^{-3}`$ | $`3.97\times10^{-2}`$ | $`0.172`$ |
| $`8\pi`$, $`128\times512`$ | $`2.2\times10^{-3}`$ | $`2.74\times10^{-2}`$ | $`0.169`$ |

Both readings converge to their own exact reference, and the distance
between them does not close. The $`\hbar/2`$ reading converges more slowly
in the reach; WC2(b) predicts no difference in kind, and the rate is left
as an observation.

### 7.2 What Q11 becomes

Proposition Q11 placed a bank on every unpaired world, and the one-particle
ladder could not tell a world-particle's bank from its world's. Theorem WC2 can:
the bank that realises the generator is indexed by the channels of the
world's joint chord, and for a pair potential those are exchanges, so it
must fire both particles at once. A bank owned by a world-particle fires one
world-particle and realises B1. Read in four dimensions, Q11 says: *every world
carries one bank per bond; each bank's integrators advance by the bond's
exchange rates and a firing moves both particles of the bond*. For one
degree of freedom the only bond is the particle's coupling to the external
potential, which is why "unpaired world" and "world" gave the same bank.

Q-SP14 asked whether the bank could be the sea's own state. §3.3 sharpens
the question: the integrand is a reading, located at a chord's midpoint, so
the bank belongs at the reader — the parent world — and is indexed by the
bond being read. Whether its memory could be held by the chords themselves
rather than by the reader is WC-SP3.

### 7.3 Proposition WC5: banks in the world ensemble (Part C)

Signed worlds $`(P, r, p_r, \varepsilon)`$ sampled from
$`\lvert g\otimes W_{\psi_0}\rvert`$ with $`\nu = 10^5`$ or $`10^6`$ worlds per unit of
$`\lVert W\rVert_1`$; (S) by velocity Verlet in $`r`$ under the exact force,
$`P`$ free; (D) by Q11's integrate-and-fire banks fed by the mesh kernel
$`K_q(r)`$, each firing an emission whose two daughters get fresh uniform
integrator phases; then a contact sink removes opposite-sign worlds in the
same mesh cell ($`\Delta P = 0.25`$, $`\Delta r = 0.156`$,
$`\Delta p_r = 0.125`$), which leaves every cell's net count unchanged.

- **W:** one bank of 31 integrators per world, $`\dot C_q = K_q(r)`$; a
  firing moves $`p_r`$ by $`\pm\xi_q`$ and leaves $`P`$ alone.
- **B1:** one bank per world-particle, $`\dot C^{(1)}_q = K_q(r)`$ and
  $`\dot C^{(2)}_q = -K_q(r)`$; a firing of particle $`i`$ moves $`p_i`$ by
  $`\pm\xi_q`$, so $`P`$ by $`\pm\xi_q`$ and $`p_r`$ by $`\pm\xi_q/2`$.

**Proposition WC5 (measured).** *Eight seeds per variant at each density:
$`\nu = 10^5`$ with W to $`t = 4`$ and B1 to $`t = 1`$, and $`\nu = 10^6`$ with
W to $`t = 4`$ and B1 to $`t = 0.5`$ (Kaggle, §10).*

*(a) Bookkeeping. No daughter of a world bank leaves its parent's $`P`$
(0 of about $`3.4\times10^6`$ firings per seed at $`\nu = 10^5`$, of
$`2.6\times10^7`$ at $`\nu = 10^6`$); every daughter of a world-particle bank does.*

*(b) Each ensemble realises its own generator, up to the contact sink's
bias. The covariance of $`P^2`$ with $`v(r)`$, mean ± standard error over
seeds:*

| $`t`$ | W, $`\nu = 10^5`$ | W, $`\nu = 10^6`$ | W mesh | B1, $`\nu = 10^5`$ | B1, $`\nu = 10^6`$ | B1 mesh |
|---|---|---|---|---|---|---|
| 0.25 | $`0.003 \pm 0.002`$ | $`0.002 \pm 0.001`$ | 0 | $`-0.509 \pm 0.017`$ | $`-0.472 \pm 0.003`$ | $`-0.501`$ |
| 0.5 | $`0.008 \pm 0.005`$ | $`0.003 \pm 0.002`$ | 0 | $`-0.614 \pm 0.060`$ | $`-0.586 \pm 0.008`$ | $`-0.635`$ |
| 1.0 | $`0.023 \pm 0.010`$ | $`0.008 \pm 0.002`$ | 0 | $`-0.11 \pm 0.38`$ | — | $`-0.421`$ |

*At $`\nu = 10^6`$ and $`t = 0.25`$ the coarse-grained relative Wigner
function of the world-bank ensemble is $`0.020`$ from W's mesh and $`0.106`$
from B1's, and that of the world-particle-bank ensemble $`0.040`$ from B1's mesh
and $`0.132`$ from W's (relative $`L^2`$ on $`4\times4`$ blocks).*

*(c) Cost. World-particle banks carry more worlds and heat them in $`P`$, while
the signed $`\langle P^2\rangle`$ stays at its exact value 1 in the mean:*

| | $`t`$ | worlds per unit $`\nu`$, W | B1 | unsigned $`\langle P^2\rangle`$, W | B1 |
|---|---|---|---|---|---|
| $`\nu = 10^5`$ | 0.5 | $`1.21`$ | $`2.03`$ | $`1.46`$ | $`18.0`$ |
| $`\nu = 10^5`$ | 1.0 | $`1.47`$ | $`4.09`$ | $`2.1`$ | $`76`$ |
| $`\nu = 10^6`$ | 0.5 | $`1.04`$ | $`1.23`$ | $`1.10`$ | $`6.45`$ |

*The excess shrinks as $`\nu`$ grows, because a denser ensemble cancels more
in the contact sink, but at both densities it grows with time. The noise on
every $`P`$-observable grows with the heating: at $`\nu = 10^5`$ the standard
error of B1's covariance rises twenty-fold between $`t = 0.25`$ and
$`t = 1`$.*

*(d) The contact sink biases both. The world bank's $`\langle v\rangle`$
leaves its mesh value steadily: by $`+0.056`$ at $`t = 1`$ and $`+0.28`$ at
$`t = 4`$ at $`\nu = 10^5`$, by $`+0.043`$ and $`+0.19`$ at $`\nu = 10^6`$;
B1's covariance at $`\nu = 10^6`$ is 6 to 8 per cent short of its mesh
value. It is the sub-cell bias of step 21's Proposition M9, which a
matching cell of the mesh's size produces whatever the bank, and it falls
only slowly with $`\nu`$.*

*Reading it.* (a) is M1 made operational: a world bank cannot move a
daughter off its parent's centre-of-mass momentum, because every channel
it owns is an exchange. (b) shows the integrators doing what Q11 says they
do, in four dimensions, for both kinds of bank: the world-particle banks are not
a noisy version of the world's, they are a faithful realisation of a
different generator. (c) is what a per-world-particle ontology costs in
worlds, as distinct from what it costs in the observable: B1's signed
$`P`$-marginal is exact (WC2(c)), but its unsigned ensemble does a random
walk in $`P`$ with steps up to $`\xi_{31} = 7.75`$, the "heating" that
Cyganski's one-particle kicks produce and that signed cancellation hides
from $`E`$ ([`poisson_kicks_and_pair_branching.md`](../supplement/poisson_kicks_and_pair_branching.md)). The world bank's own unsigned $`\langle P^2\rangle`$
also rises, to $`2.1`$ by $`t = 1`$ at $`\nu = 10^5`$, with no event moving
$`P`$: the sink cancels most where the population in $`P`$ is densest. (d) limits what the
ensemble runs can decide: over long times the sink bias is as large as the
difference between the readings, so the long-time discrimination is the
mesh's (§7.1), and the ensemble runs answer the questions about
bookkeeping and cost.

### 7.4 The Eckart pair collision (Part T)

Step 21's collision, $`U = \mathrm{sech}^2((x_1 - x_2)/2)`$, from a product of
Gaussians with relative momentum 1: only the relative marginal decides the
transmission, so by WC2(b) B1's transmission is the $`\hbar/2`$ value.

| | W (QM) | B1 ($`\hbar/2`$) | classical |
|---|---|---|---|
| exact (split-operator) | $`0.544146`$ | $`0.512035`$ | — |
| mesh, reach $`2\pi`$ | $`0.5197`$ | $`0.4629`$ | $`0.4511`$ |
| mesh, reach $`4\pi`$ | $`0.5482`$ | $`0.4942`$ | |
| mesh, reach $`8\pi`$ | $`0.5495`$ | $`0.5169`$ | |

The exact $`\hbar`$ value agrees with step 21's $`0.544134`$. At $`\hbar/2`$
the initial $`\rho'`$ is a positive mixed Gaussian (25 eigenvectors above
$`10^{-12}`$), propagated eigenvector by eigenvector. At the largest reach
each mesh is within $`5\times10^{-3}`$ of its own reference, and the shift
B1 − W is $`-0.0325`$ against the exact $`-0.0321`$. A world-particle that reads
only its own chord tunnels as if Planck's constant were halved: less than
QM, more than classical.

---

## 8. What the reader is, restated

The glossary says a reader is "whoever supplies the lever". With the words
of §1: a reader is a world, reading the misalignments of chords whose
midpoint it occupies, along the directions in which its potential has
terms. For a single particle in an external potential that is the familiar
reader. For a pair potential it is the world reading the relative
direction of its configuration space — not either world-particle, and not the
centre of mass, along which a pair potential's chords are dark. A
world-particle reads nothing of its own; it is a coordinate of a reader. A
particle, recurring in every world, reads nothing at all: its state is a
marginal of what the worlds carry.

### 8.1 Addendum: worlds suffice

*Addendum, October 2026.* Sections 3 to 5 showed that a phase owned by a
world-particle is either redundant or wrong. This section turns the
result into ontology: what each world carries, what is derived from it,
and what holds *between* worlds — the relations that carry all of the
quantum dynamics, since nothing quantum happens inside a world (WC6).

**Definition WC7 (the world ontology).** Under (S′) with persistent
identity:

*(i) Carriers.* The only carriers are worlds. A world $`w`$ carries its
point $`(\mathbf x^w, \mathbf p^w) \in \mathbb R^{6N}`$, its sign
$`\varepsilon^w`$, its clock $`\theta^w`$ and one bank of integrators
$`C^w_{b,q}`$ for each term $`b`$ of the potential (a bond, for a pair
potential; WC4).

*(ii) Derived properties.* Everything else about one world is a function of
its point: the world-particles $`(\mathbf x_i^w, \mathbf p_i^w)`$, the centre
of mass $`\mathbf X_{\rm cm} = \sum_i m_i\mathbf x_i/M`$ and total momentum
$`\mathbf P = \sum_i\mathbf p_i`$, relative and Jacobi coordinates, the
Lagrangian and the energy. A sector clock, such as $`\theta_{\rm cm}`$, is
derived only where $`\mathcal L`$ separates in those coordinates (WC1(c)).
The state of a particle is a marginal of the ensemble (WC0).

*(iii) Relations between worlds.* The quantum content is a short list of
relations, by the number of worlds each involves:

| relation | worlds | what it is | defined |
|---|---|---|---|
| chord | 2 | the misalignment $`\mu^{(\mathbf P)}_{ww'} = \theta^w - \theta^{w'} + \mathbf P\cdot(\mathbf x^{w'} - \mathbf x^w)/\hbar`$ of two worlds' clocks, carried to a common point by a lever $`\mathbf P`$; the only place the potential acts beyond the classical force | P2, step 5, L8 |
| reading | a reader and its chords | the interference sum of the chords about a world, which estimates the kernel $`K_q`$; the world's bank integrates it | L1, Q5, Q11, WC4 |
| emission | 3 | a parent fires; the aligned pair in its cell and row splits, its worlds kinking by $`\pm\pmb\xi_q`$ along a bond | specification §5.3, WC6(c) |
| absorption | 3 | a parent fires; two unpaired worlds at $`\mathbf p \pm \pmb\xi_q`$ kink back and become an aligned pair | specification §5.3, WC6(c) |
| pair collision | 4 | two aligned pairs exchange momentum (Focus, Defocus) | Definition (C), step 24 |
| dark catalysis | 4 | an aligned pair and a partner pair exchange phase, moving nothing in $`E`$ | L5, Q6 |
| admissibility | all | which signed ensembles may exist at all: $`E`$ must be the Wigner function of a positive operator | postulate (A), G3, G4 |

The table states the ontology's slogan exactly: classical within a world,
quantum between worlds. It is the same commitment as Hall, Deckert and
Wiseman's (2014), that the interaction between worlds is fundamental, with a
different interaction: theirs is a force between neighbouring worlds,
ours is a phase relation across chords that sets the rates of contact
processes between worlds, and the classical force stays inside each world.

**Proposition WC8 (worlds suffice).** *For $`N`$ particles in three
dimensions and any potential $`V(\mathbf x)`$:*

*(a) The compensated residual symbol on $`\mathbb R^{6N}`$, and with it every
kernel $`K_q`$ and every bank's integrand, is a function of world-clock
misalignment rates across chords alone: for two worlds at
$`\mathbf X \pm \mathbf Y`$ on one momentum row,*

```math
\hbar\,\dot\mu = -\big[V(\mathbf X + \mathbf Y) - V(\mathbf X - \mathbf Y)\big]
\quad\text{for every lever } \mathbf P,
```

*which is $`i\hbar`$ times the symbol $`(i/\hbar)[V(\mathbf X + \mathbf Y) - V(\mathbf X - \mathbf Y)]`$
with $`\mathbf Y = \hbar\pmb\sigma/2`$, before its linear part is removed.*

*(b) If $`V = \sum_i V_i(\mathbf x_i)`$, world-particle clocks with
$`\dot\theta_i = (T_i - V_i)/\hbar`$ sum to the world clock and their chord
misalignments sum to the world's: they are a decomposition of the world
clock and add nothing.*

*(c) If a term of $`V`$ is shared by two or more world-particles, no
assignment of world-particle clocks reproduces both the action and the
forces (WC1(d)), and readings by world-particle chords differ from the
world's beyond linear order (WC2).*

*(d) Hence no phase owned by a world-particle enters the generator, and none
can enter it correctly.*

*Proof.* (a) Along the two worlds' common row their velocities are equal,
so their separation is constant and the lever term
$`\mathbf P\cdot(\mathbf x^{w'} - \mathbf x^w)/\hbar`$ does not change; their
kinetic windings are equal; what is left of
$`(\mathcal L(\mathbf X + \mathbf Y) - \mathcal L(\mathbf X - \mathbf Y))/\hbar`$ is
the potential difference (Part G5, three particles, residual $`0`$). (b)
Additivity, term by term (Part G2, residual $`0`$). (c) is WC1(d) and
WC2(a). (d) follows from (a)–(c). $`\square`$

**Corollary WC8.1 (locality from chords).** *Let a chord move only the
particles of a set $`A`$: $`\mathbf y_k = 0`$ for $`k \notin A`$. Then
$`V(\mathbf X + \mathbf Y) - V(\mathbf X - \mathbf Y)`$ contains exactly the terms of
$`V`$ that involve a particle of $`A`$; every other term cancels. A reading
along the directions of $`A`$ is blind to everything that does not act on
$`A`$.*

*Proof.* A term whose arguments do not include any $`\mathbf x_k`$,
$`k \in A`$, takes the same value at $`\mathbf X + \mathbf Y`$ and
$`\mathbf X - \mathbf Y`$. Part G3 checks it for three particles with one-,
two- and three-body terms: a chord on particle 1 reads $`V_1`$, $`v_{12}`$,
$`v_{13}`$ and $`w_{123}`$ and nothing else. $`\square`$

The usual reason for wanting a phase per particle is locality: what is done
to one particle should not change what another reads. WC8.1 gets that from
the geometry of chords, not from phases attached to particles.

**What a worlds-only ontology still has to say.** Three places where
per-particle attributes might seem to return; in each, the natural reading
is a relation between worlds, but none is worked out here.

- *Identical particles.* A permutation $`\pi`$ of particle labels maps a
  world $`w`$ to a world $`\pi w`$, and for a symmetric potential the
  compensated symbol is invariant under permuting the $`(\mathbf x_i, \pmb\sigma_i)`$
  together (Part G4, all six permutations of three particles), so the
  dynamics preserves any fixed relation between $`w`$ and $`\pi w`$. Bose
  and Fermi statistics would be such a relation — a sign or phase between
  a world and its permutations — not a property of a particle (WC-SP7).
- *Spin.* Step 17 §8 leaves it out. A world would carry a discrete spin
  label per world-particle, labels being coordinates like positions, and
  the phase between spin components would be a chord between worlds whose
  labels differ (WC-SP8).
- *The density of worlds.* A reading needs chords, so how many worlds there
  are — the sea's density in $`\mathbb R^{6N}`$ — is part of the ontology
  and not a numerical parameter. Its cost has been measured (V7, V14, Q12)
  but not fixed (V-SP3, Q-SP13).

### 8.2 Addendum: three and four particles in many worlds

*Second addendum, October 2026.* Sections 5 to 7 tested the world ontology
with two particles, where a pair potential has one bond and the interaction
reads one relative coordinate. This section asks whether Definition WC7 and
Proposition WC8 are enough for three or more particles, where every world
carries several banks and the bonds share particles; what the worlds do with
entanglement; and how the result compares with three-particle problems whose
answers are known. Companion demo: `src/demo_three_worlds.py` (Parts S, K, E,
Q, M, W, X, F).

**Setting.** $`N`$ particles of unit mass on a line, $`\hbar = 1`$, in a
harmonic trap of unit frequency, every pair coupled by the same potential
$`v`$:

```math
V(\mathbf x) = \sum_{i=1}^{N}\tfrac12 x_i^2 + \sum_{i < j} v(x_i - x_j).
```

At $`t = 0`$ the state is the product of trap ground states, the ground state
for $`v = 0`$, and the interaction is switched on (a quench). One dimension
per particle keeps the references exact; the structure below holds in three
dimensions (WC8.1). Two pair potentials:

- **Q**, $`v(r) = \tfrac12 k r^2`$ with $`k = 2`$: the Moshinsky model
  (Moshinsky 1968), whose three-particle entanglement Bouvrie *et al.*
  (2012) computed in closed form. Quadratic, so solvable exactly by the flow
  of the covariance matrix.
- **G**, $`v(r) = g\,e^{-r^2/2s^2}/(\sqrt{2\pi}\,s)`$ with $`g = 2`$,
  $`s = 0.7`$: a contact interaction of strength $`g`$ smoothed to width
  $`s`$, the form used in multiconfigurational studies of few bosons in traps
  (Zöllner, Meyer and Schmelcher 2006, 2008, there with $`s = 0.05`$). No
  closed form.

The initial state and $`V`$ are both symmetric under exchange, so G is a
three-boson (and four-boson) problem of that kind, and the one-body purity
below measures fragmentation as well as entanglement.

**Proposition WC9 (bonds suffice for $`N`$ particles).** *For $`V`$ as above
and a chord $`\mathbf y = \hbar\pmb\sigma/2`$, the compensated residual*

```math
R(\mathbf x, \mathbf y) = V(\mathbf x + \mathbf y) - V(\mathbf x - \mathbf y) - 2\,\mathbf y\cdot\nabla V(\mathbf x)
```

*is a sum over the $`N(N-1)/2`$ bonds $`b = (i, j)`$ of one function of two
variables,*

```math
R(\mathbf x, \mathbf y) = \sum_{i < j}\rho_v(r_b, \eta_b), \qquad
r_b = x_i - x_j,\quad \eta_b = y_i - y_j,\qquad
\rho_v(r, \eta) = v(r + \eta) - v(r - \eta) - 2\eta\,v'(r),
```

*and the trap contributes nothing. Hence:*

*(a) a world carries one bank per bond, and every bank reads the same kernel
table $`K_q(r)`$: the Fourier coefficients in $`\eta`$ of
$`(i/\hbar)\rho_v(r, \eta)`$ times the raised-cosine horizon at the reach
$`y_{\max}`$, with channels $`\xi_q = q\,\Delta p`$,
$`\Delta p = \pi\hbar/2y_{\max}`$;*

*(b) a firing of channel $`q`$ on bond $`(i, j)`$ moves $`(p_i, p_j)`$ by
$`(\pm\xi_q, \mp\xi_q)`$, leaving the total momentum and the centre of mass
alone; in orthonormal Jacobi coordinates the transfer is $`\xi_q\mathbf u_b`$,
with $`\mathbf u_b`$ the bond's relative direction and $`\lvert\mathbf u_b\rvert^2 = 2`$;*

*(c) if $`v`$ is quadratic every $`K_q`$ vanishes (Theorem G2): no world ever
fires, and the exact evolution, entanglement included, is classical
streaming of the initial Wigner function.*

*Proof.* The quadratic terms of $`V`$ have a zero third difference and drop
out of $`R`$. A pair term $`v(x_i - x_j)`$ at $`\mathbf x \pm \mathbf y`$ has
argument $`r_b \pm \eta_b`$, and its part of $`\mathbf y\cdot\nabla V`$ is
$`y_i v'(r_b) - y_j v'(r_b) = \eta_b v'(r_b)`$. The Fourier transform in
$`\pmb\sigma`$ of a function of $`\sigma_i - \sigma_j`$ alone is supported on
transfers along $`\mathbf e_i - \mathbf e_j`$, which gives (a) and (b). For (c),
$`\rho_v = 0`$ when $`v`$ is quadratic. $`\square`$

Part S checks the decomposition at 2000 random chords (largest difference
$`6.7\times10^{-15}`$ for three particles and $`1.0\times10^{-14}`$ for four,
against residuals up to 16 and 19), and that a bond event conserves $`P`$ and
moves the Jacobi momenta by $`\xi\mathbf u_b`$. WC9 makes WC4 specific: Q11's
bank per unpaired world becomes $`N(N-1)/2`$ banks per world sharing one
table. The generator reads a world's point and its banks and nothing that
belongs to a world-particle, as WC8 requires.

**The kernel (Part K).** At reach $`2\pi`$, $`4\pi`$, $`8\pi`$ the table has 31,
63 and 127 channels ($`\Delta p = 0.25`$, $`0.125`$, $`0.0625`$). The largest
rate per bond, $`\Gamma_{\max} = \max_r\sum_q\lvert K_q(r)\rvert`$, is $`2.03`$,
$`4.25`$, $`8.32`$, proportional to the reach as G5 says; summed over the
bonds and averaged over the initial state it is $`3.8`$, $`8.1`$, $`15.7`$ per
world per unit time for three particles and twice that for four. Narrowing
the Gaussian toward a contact interaction at fixed $`g`$ raises
$`\Gamma_{\max}`$ as $`s^{-2}`$: $`4.25`$, $`16.6`$, $`50.1`$, $`196`$ for
$`s = 0.7`$, $`0.35`$, $`0.2`$, $`0.1`$ at reach $`4\pi`$. The mean over the
initial state grows more slowly ($`8.1`$ to $`88`$).

**References (Parts E, Q).** For G, the time-dependent Schrödinger equation
on a grid in particle coordinates (split operator, $`\Delta t = 0.005`$):
$`64^3`$ points on $`[-7, 7)^3`$ and $`96^3`$ on $`[-8, 8)^3`$ agree to
$`4.7\times10^{-10}`$ in the purity; for four particles, $`40^4`$ on
$`[-6, 6)^4`$ and $`48^4`$ on $`[-6.5, 6.5)^4`$ agree to $`6.1\times10^{-8}`$.
For Q the covariance flow is exact and the grid reproduces it to
$`3.6\times10^{-5}`$. A *classical control* streams the initial, positive
Wigner function under the full force with no events.

**Entanglement.** The measure is the purity of particle 1's reduced state
and its Rényi entropy,

```math
P_1 = \mathrm{Tr}\,\rho_1^2 = 2\pi\hbar\int W_1(x_1, p_1)^2\,dx_1\,dp_1,
\qquad S_2 = -\ln P_1 ,
```

where $`W_1`$ is the marginal of the Wigner function on particle 1;
$`P_1 = 1`$ exactly when particle 1 is not entangled with the rest. For the
worlds, $`W_1`$ is estimated by projecting the relative ensemble on
$`x_1 = \mathbf a\cdot\pmb\rho`$, $`p_1 = \mathbf a\cdot\pmb\pi`$ ($`\pmb\rho`$,
$`\pmb\pi`$ orthonormal Jacobi coordinates, $`\mathbf a = (1/\sqrt2, 1/\sqrt6)`$
for three particles) and convolving with the centre of mass's exact
Gaussian, which no event touches (WC9(b)). Pairs of a world with itself are
removed, so the estimator is unbiased.

A second measure tests the whole state rather than one particle. Unitary
evolution conserves $`(2\pi\hbar)^{N-1}\int W_{\rm rel}^2 = 1`$, where
$`W_{\rm rel}`$ is the Wigner function of the relative motion. A finite
ensemble cannot resolve that integral, but it can resolve it at a fixed
resolution $`a`$:

```math
P_a = (2\pi\hbar)^{N-1}\int (W_{\rm rel} * G_a)^2 ,
```

with $`G_a`$ a Gaussian of width $`a`$ in every Jacobi coordinate. For the
mesh, $`P_a`$ is computed spectrally and exactly. For the worlds it is
computed by cloud-in-cell deposition with the self-pairs removed. $`P_a`$
falls as the true $`W`$ develops structure finer than $`a`$, so only the
difference between worlds and mesh at the same reach is meaningful. Here
$`a = 0.25`$.

**Three answers.**

*1. Entanglement without events (problem Q).* No event fires (largest rate
$`1.1\times10^{-13}`$). Positive worlds streamed by Verlet reproduce the exact
purity within $`1.5\times10^{-3}`$ (three particles, $`10^6`$ worlds) and
$`2.5\times10^{-3}`$ (four). The purity falls to $`0.683`$ ($`S_2 = 0.38`$) and
$`0.655`$ ($`S_2 = 0.42`$), and revives with period $`\pi/\sqrt{1 + Nk}`$, the
half-period of the relative modes. A quadratic interaction entangles by
classical streaming alone, so entanglement is not by itself a sign of the
quantum events: the worlds carry it as a correlation between their
coordinates.

*2. The generator is right (problem G, Part M).* The bond-sum generator of
WC9 was run without sampling on a mesh of the four-dimensional Jacobi
relative phase space: $`40^2`$ points on $`[-5, 5)^2`$ in $`\pmb\rho`$, with the
$`\pmb\pi`$ grid chosen so that the bond chords fit it. It converges to the
exact three-particle dynamics as $`1/y_{\max}^2`$:

| reach | purity error, $`t = 1`$ | $`t = 2`$ | $`t = 3`$ | $`\langle\sum v\rangle`$ error, $`t = 3`$ |
|---|---|---|---|---|
| $`2\pi`$ | $`-0.016`$ | $`-0.054`$ | $`-0.117`$ | $`-0.185`$ |
| $`4\pi`$ | $`-0.004`$ | $`-0.015`$ | $`-0.039`$ | $`-0.074`$ |
| $`8\pi`$ | $`-0.001`$ | $`-0.004`$ | $`-0.009`$ | $`-0.020`$ |

The exact purity falls to $`0.800`$ at $`t = 2`$ and recovers to $`0.879`$ by
$`t = 3`$. Classical streaming falls to $`0.606`$ and is still at $`0.663`$
at $`t = 3`$. Classical dynamics therefore over-entangles, and the events
have to undo that excess; at reach $`8\pi`$ the generator does so to within
$`0.01`$. The mesh also conserves $`(2\pi\hbar)^2\int W_{\rm rel}^2`$, to
$`2\times10^{-3}`$ at reach $`2\pi`$ by $`t = 3`$.

*3. The worlds realise the generator early and lose it late (Parts W, X).*
The realisation uses integrate-and-fire banks per bond and emission events,
followed after each step by the contact sink on cells of side $`c`$ in the
Jacobi relative phase space. The worlds follow the mesh at the same reach
while the quantum correction is small, then fall below it. The share of the
quantum correction they keep is

```math
q(t) = \frac{P_1^{\rm worlds}(t) - P_1^{\rm classical}(t)}{P_1^{\rm mesh}(t) - P_1^{\rm classical}(t)} ,
```

which is $`1`$ when the worlds realise the generator and $`0`$ when they do no
better than classical streaming. Worlds minus mesh at the same reach:

| reach | cell $`c`$ | $`\nu`$ | variant | $`N/\nu`$ at $`t = 3`$ | $`\Delta P_1`$: $`t = 1`$ | $`t = 2`$ | $`t = 3`$ | $`q(1)`$ | $`q(2)`$ | $`\Delta\langle\sum v\rangle`$, $`t = 3`$ | $`\Delta P_{0.25}`$, $`t = 3`$ |
|---|---|---|---|---|---|---|---|---|---|---|---|
| $`2\pi`$ | 0.35 | $`10^5`$ (2 seeds) | — | 2.27 | $`-0.025`$ | $`-0.060`$ | $`-0.090`$ | 0.79 | 0.57 | $`-0.086`$ | $`-0.145`$ |
| $`2\pi`$ | 0.35 | $`10^6`$ | — | 1.23 | $`-0.019`$ | $`-0.054`$ | $`-0.078`$ | 0.84 | 0.61 | $`-0.084`$ | $`-0.139`$ |
| $`2\pi`$ | 0.5 | $`10^5`$ (2 seeds) | — | 1.34 | $`-0.036`$ | $`-0.088`$ | $`-0.124`$ | 0.69 | 0.37 | $`-0.133`$ | $`-0.187`$ |
| $`2\pi`$ | 0.5 | $`10^6`$ | — | 1.07 | $`-0.034`$ | $`-0.081`$ | $`-0.115`$ | 0.71 | 0.42 | $`-0.132`$ | $`-0.184`$ |
| $`2\pi`$ | 0.7 | $`10^5`$ (2 seeds) | — | 1.12 | $`-0.041`$ | $`-0.103`$ | $`-0.141`$ | 0.65 | 0.27 | $`-0.183`$ | $`-0.207`$ |
| $`2\pi`$ | 0.7 | $`10^6`$ | — | 1.03 | $`-0.043`$ | $`-0.099`$ | $`-0.135`$ | 0.64 | 0.30 | $`-0.182`$ | $`-0.207`$ |
| $`2\pi`$ | 0.35 | $`10^5`$ | $`\Delta t = 0.005`$ | 2.19 | $`-0.020`$ | $`-0.059`$ | $`-0.080`$ | 0.83 | 0.58 | $`-0.082`$ | |
| $`2\pi`$ | 0.35 | $`10^5`$ | sink every 4 steps | 3.16 | $`-0.022`$ | $`-0.060`$ | $`-0.077`$ | 0.81 | 0.57 | $`-0.086`$ | |
| $`2\pi`$ | 0.35 | $`10^5`$ | sink levels $`c/4, c/2, c`$ | 2.28 | $`-0.018`$ | $`-0.043`$ | $`-0.070`$ | 0.84 | 0.69 | $`-0.048`$ | |
| $`2\pi`$ | 0.7 | $`10^5`$ | five sink levels | 1.11 | $`-0.041`$ | $`-0.086`$ | $`-0.118`$ | 0.65 | 0.39 | $`-0.103`$ | |
| $`4\pi`$ | 0.25 | $`2\times10^5`$ | — | 8.21 | $`-0.015`$ | $`-0.047`$ | $`-0.100`$ | 0.88 | 0.74 | $`-0.071`$ | |
| $`4\pi`$ | 0.35 | $`10^5`$ (2 seeds) | — | 4.05 | $`-0.019`$ | $`-0.063`$ | $`-0.121`$ | 0.85 | 0.65 | $`-0.113`$ | |
| $`4\pi`$ | 0.35 | $`3\times10^5`$ | — | 2.19 | $`-0.009`$ | $`-0.050`$ | $`-0.100`$ | 0.93 | 0.72 | $`-0.108`$ | |
| $`4\pi`$ | 0.5 | $`10^5`$ | — | 1.70 | $`-0.018`$ | $`-0.076`$ | $`-0.144`$ | 0.86 | 0.58 | $`-0.162`$ | |
| $`4\pi`$ | 0.7 | $`10^5`$ | — | 1.21 | $`-0.000`$ | $`-0.072`$ | $`-0.165`$ | 1.00 | 0.60 | $`-0.249`$ | |

$`\Delta`$ is worlds minus mesh at the same reach. Seed-to-seed scatter is
about $`0.005`$ in $`P_1`$. $`q(3)`$ is left out because at reach $`2\pi`$ the
mesh and classical streaming differ by only $`0.10`$ at $`t = 3`$, which
makes the ratio unstable. The global purity was added to the demo after the
reach-$`4\pi`$ runs and the diagnostics.

Six things follow from the table:

- *The time step is not the cause.* Halving $`\Delta t`$ moves the
  deficit at $`t = 3`$ from $`-0.090`$ to $`-0.080`$, and at $`t = 1`$ and $`2`$ by
  less than $`0.005`$, which is about one seed-to-seed scatter.
- *How often the sink acts is not the cause either.* Running it every fourth
  step moves the deficit by similar amounts ($`-0.077`$ at $`t = 3`$), even
  though the population then rises to $`3.2\nu`$.
- *The kernel is not the cause.* It reproduces the bond symbol to
  $`7\times10^{-4}`$ at reach $`2\pi`$ and $`2\times10^{-4}`$ at $`4\pi`$.
- *More worlds help only a little.* At a fixed cell, going from
  $`\nu = 10^5`$ to $`10^6`$ moves the purity at $`t = 3`$ from $`-0.090`$ to
  $`-0.078`$ (reach $`2\pi`$, $`c = 0.35`$), while $`N/\nu`$ falls from
  $`2.3`$ to $`1.2`$.
- *The cell sets the limit.* At $`\nu = 10^6`$ the deficit at reach $`2\pi`$
  is close to linear in the cell: a straight line through $`c = 0.35`$,
  $`0.5`$, $`0.7`$ extrapolates to $`+0.002`$, $`-0.014`$, $`-0.027`$ in
  $`P_1`$ at $`t = 1`$, $`2`$, $`3`$ and to $`+0.011`$ in
  $`\langle\sum v\rangle`$, but only to $`-0.079`$ in $`P_{0.25}`$. A fit in
  $`c^2`$ leaves $`-0.068`$ in $`P_1`$ at $`t = 3`$. The bias is first order
  in the cell, so a cell² extrapolation is the wrong model. Three cells,
  one seed each, do not establish convergence.
- *Pairing finer cells first helps the energy more than the purity.* Running
  the sink at $`c/4`$, then $`c/2`$, then $`c`$ halves the bias in
  $`\langle\sum v\rangle`$ but removes only about a fifth of the purity's.

The whole state shows the same deficit: $`P_{0.25}`$ at $`t = 3`$ is $`0.14`$
below the mesh at $`c = 0.35`$ and $`0.21`$ below at $`c = 0.7`$ (reach
$`2\pi`$, $`\nu = 10^5`$). The world ensemble has lost purity that the mesh
keeps.

**Proposition WC10 (the contact sink decoheres).** *Let $`\Pi_c`$ remove, in
every cell of side $`c`$ of the relative phase space, as many positons and
negatons as it can, in equal numbers and chosen at random.*

*(a) $`\Pi_c`$ leaves every cell's net weight unchanged, so it commutes with
averaging over cells: it changes $`W`$ only below the resolution of a cell.*

*(b) In a cell that holds both signs, the survivors are a random subset of
the majority sign, so the shape of $`W`$ inside the cell becomes the
majority's. Sign structure finer than $`c`$ — interference fringes — is
removed.*

*(c) (Measured, three particles.) The realised evolution is not unitary:
its $`P_{0.25}`$ falls below the mesh's at the same reach and does not
return. The one-body purity falls with it, toward classical streaming.*

*(d) (Measured.) The deficit does not depend on the time step or on how
often the sink acts. At fixed $`c`$ it approaches a nonzero limit as $`\nu`$
grows, and that limit is roughly linear in $`c`$. Driving the one-body
purity's error to zero therefore requires $`\nu \to \infty`$ first and then
$`c \to 0`$, at a population that grows like the number of occupied cells,
$`c^{-2(N-1)}`$ on a line.*

*Proof of (a) and (b).* By construction. $`\square`$

The analogy: $`\Pi_c`$ acts like an environment that records which cell each
world is in. Like any which-path record, it destroys the interference it can
resolve, and it pushes the reduced state toward the classical one, as
decoherence does. The quantity the sink is meant to conserve, the signed
ensemble's $`W`$, it conserves only at the resolution of its cells.

This is step 21's Proposition M9 seen through the observables. There,
matching inside a cell biased $`\langle V\rangle`$, and finer matching cut
the bias at a price in unpaired worlds. Here the bias shows up as lost
global purity and over-estimated entanglement, and it is first order in the
cell. It is also the deterministic error that Xiong and Shao
(2019, Theorem 4.1) bound for annihilation on a uniform histogram, and their
remedy, about ten samples per occupied bin, is the population cost of
Corollary M10.

*Four particles.* Six bonds per world and a six-dimensional relative phase
space. At reach $`2\pi`$, $`c = 0.7`$, $`\nu = 5\times10^4`$, the worlds give
$`P_1 = 0.854`$, $`0.692`$, $`0.581`$ at $`t = 0.5`$, $`1`$, $`1.5`$. The exact
values are $`0.907`$, $`0.814`$, $`0.762`$, and classical streaming gives
$`0.670`$ at $`t = 1`$. Meanwhile $`N/\nu`$ climbs from $`2.7`$ to $`11.9`$,
because few cells in six dimensions hold both signs. The generator is the
same bond sum (Part S). A six-dimensional mesh to check it independently is
beyond this container. The realisation keeps almost none of the quantum
correction: $`q(1) = 0.15`$. This is Corollary M10 in its plainest form: at a
fixed population, the contact sink's resolution worsens with every relative
degree of freedom added.

**Comparison with known solutions and with the best simulations.**

- *Exactly solvable problems.* The Moshinsky model (Q) is reproduced to
  sampling error for three and four particles, entanglement included, with
  no events. The integrable models of particles on a line with contact
  (McGuire 1964) or inverse-square (Calogero 1969) interactions are
  singular. By Part K the largest event rate of a smoothed contact grows as
  $`s^{-2}`$, so those models can be reached only as limits of smooth
  potentials. Zöllner *et al.*'s $`s = 0.05`$ would cost about 200 times the
  rate used here.
- *The best simulations.* For $`N \le 4`$ on a line, the direct grid TDSE
  used here is numerically exact (purity to $`10^{-9}`$ for three particles).
  For more particles the standard tools are MCTDH and its bosonic form
  MCTDHB (Beck *et al.* 2000; Alon, Streltsov and Cederbaum 2008), which
  Zöllner, Meyer and Schmelcher applied to this regularised interaction for
  a few bosons. The bond-sum generator on the mesh matches the grid TDSE as
  closely as the reach allows. The world ensemble does not yet.
- *Other signed-particle methods.* The closest algorithms are the many-body
  Wigner Monte Carlo of Sellier and Dimov (2014) and the branching random
  walk of Shao and Xiong (2019; Xiong and Shao 2019). Both use signed
  particles, events drawn from the Wigner kernel, and population control by
  annihilation on a phase-space histogram. Xiong and Shao's
  four-dimensional helium-like test used $`10^8`$ bins and more than
  $`10^9`$ particles. The deficit measured here is the same annihilation
  error at a population $`10^3`$ to $`10^4`$ times smaller.

**What an algorithm specification for $`N`$ particles needs.** Each of the
pieces has now been tested at least once:

1. *State.* A world is $`(\mathbf x, \mathbf p) \in \mathbb R^{dN}`$, a sign, and
   one bank of $`N_q`$ integrators per term of the potential (WC7, WC9). The
   clock is needed for chords, not for the generator.
2. *Tables.* One kernel table per distinct term (one for all the bonds of
   identical particles), built from the term's compensated symbol with the
   horizon at the reach; channels $`\xi_q = q\Delta p`$; transfers along the
   term's direction (WC9(a), (b)).
3. *Streaming.* The full classical force for unpaired worlds, inertial
   motion for dark pairs ((S′), WC6).
4. *Events.* Integrate-and-fire on each (world, term, channel), realised as
   emission or absorption (specification §5.3), or as emission followed by a
   sink.
5. *Population control.* This is the open part. A specification must name a
   rule with a measured error budget. The contact sink on cells of the
   relative phase space has an error that does not vanish with $`\nu`$ at
   fixed cell (WC10) and a cost exponential in the $`2d(N-1)`$ relative
   dimensions (M10).
6. *Observables.* Marginals; energies; the one-body purity with self-pairs
   removed and the centre of mass convolved exactly where it separates; and
   $`P_a`$ as the test of unitarity.
7. *Validation.* Structure (S); a quadratic interaction with no events (Q);
   the mesh at the same reach (M), which separates the generator's error
   from the realisation's; the worlds (W) against both, with $`P_a`$.

## 9. Figures

![Mesh, ensemble and transmission results](https://raw.githubusercontent.com/billpage/wpmw/output/figures/world_clocks.png)

*(a) Distance of the relative Wigner function from the stationary ground
state. W (blue) stays at the mesh's own error; B1 (red) follows the exact
$`\hbar/2`$ reference (open circles), not QM; B2 and classical streaming
are shown for scale. (b) Separability of the centre of mass from the
relative motion: W at rounding, B1 and B2 at order $`10^{-1}`$. (c) Eckart
transmission against the reach: W converges to the QM value, B1 to the
$`\hbar/2`$ value. (d) The covariance of $`P^2`$ with $`v`$: world
ensembles (mean ± standard error, 8 seeds; filled $`\nu = 10^5`$, open
$`\nu = 10^6`$) against the mesh. (e) Worlds per unit $`\nu`$ and the unsigned $`\langle P^2\rangle`$,
one line per seed, $`\nu = 10^5`$. (f) $`\langle v\rangle`$, worlds
against the mesh: the contact sink's bias, smaller at $`\nu = 10^6`$.*

![Three and four particles in many worlds](https://raw.githubusercontent.com/billpage/wpmw/output/figures/three_worlds.png)

*Second addendum (§8.2), problem G unless stated. (a) One-body purity
$`P_1`$ of particle 1: exact (grid TDSE), classical streaming, the bond-sum
generator on the Jacobi mesh at reach $`2\pi`$, $`4\pi`$, $`8\pi`$, and the
worlds at the two smaller reaches (cell $`0.35`$). (b) The pair energy,
likewise. (c) The global purity $`P_{0.25}`$ of the relative state, mesh
against worlds at reach $`2\pi`$: the mesh conserves it at full resolution,
and the worlds lose it (WC10). (d) Problem Q, three and four particles:
covariance flow, grid TDSE and worlds with no events. (e) Worlds minus mesh
at the same reach, purity at $`t = 3`$, against the sink cell; stars
$`\nu = 10^6`$, with the straight line through them extended to $`c = 0`$. (f) Four particles: exact, classical streaming and
worlds at reach $`2\pi`$, cell $`0.7`$, $`\nu = 5\times10^4`$.*

## 10. Numerical verification

All outputs are from `src/demo_world_clocks.py`, seeds fixed in the code.
Parts A, K, H, M, R, T and F were run in the container, October 2026, and
Part G with the addendum; Part C
on Kaggle as two private batch kernels, `wpmw-wclocks-c1e5` (`C --nu 100000
--seeds 8 --Tc 4 --Tb1 1`) and `wpmw-wclocks-c1e6` (`C --nu 1000000 --seeds 8
--Tc 4 --Tb1 0.5`), each with four workers, at commit `8e0caa6` with this
note's demo applied as a patch. Their CSV output is not committed; each
run is seeded, so rerunning Part C with those arguments reproduces it.

- **Part A.** SymPy with generic Taylor coefficients: WC2(a) to order 7,
  WC2(b)–(d), WC1(b), the gauge check of §3.1 and WC4 for $`N = 3`$; every
  identity prints residual $`0`$.
- **Part K.** Channel count and event rates of §6.
- **Part H.** The exact references: QM's stationary values, and the
  $`\hbar/2`$ von Neumann evolution of $`\rho'`$ (1024 points on
  $`[-20, 20)`$, exact diagonalisation), with its Wigner function on the
  mesh at sixteen times.
- **Part M.** The three readings and classical streaming on the
  $`(P, r, p_r)`$ mesh to $`t = 4`$; B1's relative marginal against an
  independent $`\hbar/2`$ mesh.
- **Part R.** The relative-motion errors at reach $`2\pi`$, $`4\pi`$, $`8\pi`$.
- **Part T.** Eckart transmissions: split-operator at $`\hbar`$ and on
  $`\rho'`$ at $`\hbar/2`$; relative-marginal meshes at three reaches (about
  18 minutes).
- **Part C.** Particle runs with integrator banks; `--seeds`, `--nu`,
  `--Tc`, `--Tb1` set the sweep, one CSV row per variant, seed and time.
- **Part G.** (Addendum.) SymPy: WC1(d)'s two conditions have no common
  solution; separable potentials give residual $`0`$ between the world's
  chord misalignment and the sum of the world-particles'; the terms a chord
  on particle 1 reads, for three particles; invariance of the symbol under
  the six joint permutations; and WC8(a) for three particles, residual
  $`0`$.
- **Part F.** The two figures.

The second addendum's outputs are from `src/demo_three_worlds.py`, run in the
container in October 2026 at commit `d802996` with the demo applied, seeds
fixed in the code (seed $`k`$ for the $`k`$-th seed of a run):

- **Part S.** The bond-sum identity at 2000 random chords for $`N = 3`$ and
  $`4`$ (largest difference $`6.7\times10^{-15}`$ and $`1.0\times10^{-14}`$); bond
  events in Jacobi coordinates.
- **Part K.** Channels and rates at three reaches; the width scan.
- **Part E.** `E` and `E --N 4`: grid TDSE on two grids each (about 16
  minutes for four particles), classical control with $`4\times10^5`$ samples.
- **Part Q.** `Q` and `Q --N 4`, $`10^6`$ worlds.
- **Part M.** `M --mesh-reach 2 4 8` (the $`8\pi`$ mesh about 30 minutes).
- **Part W.** Reach $`2\pi`$: cells $`0.35`$, $`0.5`$, $`0.7`$ at
  $`\nu = 10^5`$ (2 seeds) and $`10^6`$ (1 seed); at $`0.35`$ also `--dt 0.005`,
  `--sink-every 4`, `--levels 3`, `--levels 5`, and `--levels 5` at cell
  $`0.7`$. Reach $`4\pi`$: cells $`0.35`$ ($`\nu = 10^5`$, 2 seeds;
  $`3\times10^5`$), $`0.25`$ ($`2\times10^5`$), $`0.5`$ and $`0.7`$ ($`10^5`$).
  Four particles: `W --N 4 --reach 2 --cell 0.7 --nu 50000 --T 1.5`. The runs at
  $`\nu = 10^6`$ take about half an hour each. The 4π runs and the 2π
  diagnostics predate the global-purity observable, so their $`P_a`$ is not
  recorded.
- **Part X.** The table of §8.2 from all saved runs.
- **Part F.** The figure.

---

## 11. Open items

- **WC-SP1.** *Answered* by the amendment of October 2026 (Definition
  WC0): *body* is retired; the ledger's sense becomes *unpaired world* and
  the multi-particle notes' sense *particle*, and *world-particle* names one
  particle in one world.
- **WC-SP2.** The live reading in four dimensions. Q12 found that a sea's
  live reading converges open loop and fails closed loop at the physical
  density. For a bond's bank the chords that matter lie along the relative
  direction; does a sea whose pairs are confined to that direction supply
  them at a usable density, and is that confinement itself M1 applied to
  the sea?
- **WC-SP3.** Where the bank's memory is held. A bank belongs to a reader
  and a bond; its integrand is a sum over chords. Can the integrator's state
  be carried by the chords (pairs of worlds) rather than by the reader, and
  would that remove the hidden initial phases Q11 introduced?
- **WC-SP4.** Three-body and general potentials. WC4 puts the bank on the
  joint chord, with channels indexed by joint transfers; how does the
  channel count, and with it the cost, scale with the order of the
  interaction?
- **WC-SP5.** B1 as a known approximation. Each particle seeing the other as
  an external field is the mean-field instinct; is B1, or its $`\hbar/2`$
  relative motion, the Wigner image of a recognised approximation, and does
  that explain the factor?
- **WC-SP6.** The world-particle-bank ensemble's cost (WC5) at longer times, with
  an absorptive ledger instead of emission and a contact sink.
- **WC-SP7.** Identical particles as a relation between a world and its
  permutations (§8.1). Which relation gives Bose and Fermi statistics, how
  is it carried by chords between $`w`$ and $`\pi w`$, and what does it cost
  the sea?
- **WC-SP8.** Spin as discrete labels on world-particles, with the phases
  between spin components carried by chords between worlds whose labels
  differ (§8.1). Does the compensated split extend to the Pauli term, and
  what is a bond for a spin–spin interaction?
- **WC-SP9.** A population-control rule that conserves $`W`$ below the
  cell's resolution (§8.2, WC10). Candidates: annihilating only
  near-coincident pairs; moment-preserving resampling within a cell; or the
  absorptive ledger, which joins two unpaired worlds at the parent's row
  rather than at random in a cell. The test of any of them is $`P_a`$ at a
  fixed resolution against the mesh. This continues M-SP6. *Partly
  answered in step 26* ([`population_control.md`](population_control.md)):
  moment-preserving resampling that keeps each cell's first moment, with
  pairs at their dipole centres, cuts the deficit sevenfold at about the
  same population (Proposition PC6); absorption off a mesh needs a
  tolerance and is then as biased as a contact sink (Proposition PC5).
- **WC-SP10.** Population control in high relative dimension: six for four
  particles on a line, twelve for three in space. Cells on bond coordinates,
  or a partition that factorises over Jacobi pairs, instead of cells of the
  joint relative phase space (Corollary M10)?
- **WC-SP11.** The contact limit. The largest rate of a Gaussian of width
  $`s`$ grows as $`s^{-2}`$ (Part K). Is there a compensated form of the delta
  interaction — a boundary condition on the bond — with a finite rate?
- **WC-SP12.** An independent check of the four-particle generator (a 6D
  Jacobi mesh, or MCTDHB), to separate its error from the realisation's in
  §8.2.
- **Q-SP14** is restated in §7.2 and continues as WC-SP3. **M-SP5** (bins in
  particle coordinates) is the setting of Part C.

---

## 12. Sources

- Poirier, B. *Bohmian mechanics without pilot waves*, Chem. Phys. **370**
  (2010) 4–14. Trajectories labelled by a continuous parameter.
- Schiff, J.; Poirier, B. *Communication: Quantum mechanics without
  wavefunctions*, J. Chem. Phys. **136** (2012) 031102. The labelling as a
  diffeomorphism, the quantum force through its Jacobian.
- Hall, M. J. W.; Deckert, D.-A.; Wiseman, H. M. *Quantum phenomena modeled
  by interactions between many classical worlds*, Phys. Rev. X **4** (2014)
  041013. Finite worlds; the one-dimensional interworld potential uses the
  sorted configuration.
- Sebens, C. T. *Quantum mechanics as classical physics*, Philos. Sci. **82**
  (2015) 266–291.
- van Cittert, P. H. *Die wahrscheinliche Schwingungsverteilung in einer von
  einer Lichtquelle direkt oder mittels einer Linse beleuchteten Ebene*,
  Physica **1** (1934) 201–210.
- Zernike, F. *The concept of degree of coherence and its application to
  optical problems*, Physica **5** (1938) 785–795.
- Jennison, R. C. *A phase sensitive interferometer technique for the
  measurement of the Fourier transforms of spatial brightness distributions
  of small angular extent*, Mon. Not. R. Astron. Soc. **118** (1958) 276.
  Closure phase.
- Thompson, A. R.; Moran, J. M.; Swenson, G. W. *Interferometry and Synthesis
  in Radio Astronomy*, 3rd ed., Springer (2017), open access,
  doi:10.1007/978-3-319-44431-4. VLBI, correlators, fringe rotation and
  closure relations.
- Jackson, J. D. *Classical Electrodynamics*, 3rd ed., Wiley (1999), §1.11.
  The factor $`\tfrac12`$ in the energy of a system of charges, the double
  counting of §3.4.
- Pöschl, G.; Teller, E. *Bemerkungen zur Quantenmechanik des anharmonischen
  Oszillators*, Z. Phys. **83** (1933) 143–151. The bound problem of §6.
- Eckart, C. *The penetration of a potential barrier by electrons*, Phys.
  Rev. **35** (1930) 1303–1309. The collision of §7.4.
- Moshinsky, M. *How good is the Hartree–Fock approximation*, Am. J. Phys.
  **36** (1968) 52. Problem Q of §8.2.
- Bouvrie, P. A.; Majtey, A. P.; Plastino, A. R.; Sánchez-Moreno, P.;
  Dehesa, J. S. *Quantum entanglement in exactly soluble atomic models: the
  Moshinsky model with three electrons, and with two electrons in a uniform
  magnetic field*, Eur. Phys. J. D **66** (2012) 15.
- Zöllner, S.; Meyer, H.-D.; Schmelcher, P. *Ultracold few-boson systems in a
  double-well trap*, Phys. Rev. A **74** (2006) 053612; *Few-boson dynamics in
  double wells: from single-atom to correlated pair tunneling*, Phys. Rev.
  Lett. **100** (2008) 040401. The Gaussian-regularised contact interaction.
- Beck, M. H.; Jäckle, A.; Worth, G. A.; Meyer, H.-D. *The multiconfiguration
  time-dependent Hartree (MCTDH) method*, Phys. Rep. **324** (2000) 1–105.
- Alon, O. E.; Streltsov, A. I.; Cederbaum, L. S. *Multiconfigurational
  time-dependent Hartree method for bosons*, Phys. Rev. A **77** (2008)
  033613.
- McGuire, J. B. *Study of exactly soluble one-dimensional N-body problems*,
  J. Math. Phys. **5** (1964) 622.
- Calogero, F. *Ground state of a one-dimensional N-body system*, J. Math.
  Phys. **10** (1969) 2197.
- Sellier, J. M.; Dimov, I. *The many-body Wigner Monte Carlo method for
  time-dependent ab-initio quantum simulations*, J. Comput. Phys. **273**
  (2014) 589.
- Shao, S.; Xiong, Y. *A branching random walk method for many-body Wigner
  quantum dynamics*, Numer. Math. Theor. Meth. Appl. **12** (2019) 21;
  arXiv:1603.00159.
- Xiong, Y.; Shao, S. *The Wigner branching random walk: efficient
  implementation and performance evaluation*, Commun. Comput. Phys. **25**
  (2019) 871. Annihilation on a uniform histogram, its deterministic error
  (Theorem 4.1) and its sample requirement.
