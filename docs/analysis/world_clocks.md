# Worlds, particles and clocks: where a reading lives

> Is the reader — its lever, chords and bank of integrators — a property of a particle, or a relation between worlds that a particle inherits? With two particles and a pair potential the question has an answer. **Definition WC0** fixes the words, for `N` particles in three dimensions: a *world* is a signed point of the whole configuration's phase space `R^6N` with one clock; a *world-particle* is one particle as it appears in one world, with no sign, clock or bank of its own; a *particle* is the identity that recurs in every world, whose only state is the marginal of the signed ensemble; a *positon* or *negaton* is a world of sign +1 or -1, and an *unpaired world* is one that is not half of an aligned pair. **Proposition WC6** gives a world's history under (S′) with persistent identity: worlds and their signs are conserved, every world-particle moves on Newtonian arcs while its world is unpaired and on inertial arcs while it is dark, momenta kink only at events and only along a bond, the parent of an event never kinks, and all worlds share one time. Two things are called a clock: the phase clock of P1, which winds at the Lagrangian rate and is a gauge except in misalignments across chords, and Q11's counting clock, an integrate-and-fire counter; they are linked through chords, so the counter is the accumulator of a reading, located where a chord's midpoint is. **Proposition WC1**: the phase clock belongs to the world — world-particle clocks count an interaction twice or each carry a share of it, and clocks split along the Lagrangian's sectors (centre of mass, relative motion), not along particles. Unlike Poirier's labelled trajectories or Hall, Deckert and Wiseman's sorted worlds, WPMW asks of the set of worlds only midpoints and difference vectors, plus a phase per world; the VLBI analogy (stations with their own clocks, information only on baselines, a correlator per pair, closure phases) is set out term by term. **Theorem WC2**: reading the joint chord (W) and reading each world-particle's own chord with its clock carrying the whole interaction (B1) agree only at linear order; exactly and to all orders, B1's relative motion is the QLE at `hbar/2`, and B1 adds a coupling `3 d_P^2 d_pr` between centre of mass and relative motion that QM lacks; splitting the interaction between the clocks (B2) halves both and reads half the force. **Corollary WC3**, measured on a `(P, r, p_r)` mesh from a product state QM keeps stationary and separable: W stays separable to `7e-15` while B1's defect is 0.2 by `t = 0.5`; B1 follows the exact `hbar/2` von Neumann reference within the mesh's error and sits 0.17 to 0.18 from QM at every reach; in the Eckart pair collision B1 transmits 0.517 against W's 0.549 at reach `8 pi`, exact 0.5120 and 0.5441. **Proposition WC4**: for pair potentials a world carries one bank per bond, `N(N - 1)/2` for `N` particles like the baselines of an array, and none per world-particle; Q11's bank per free body is the one-particle case. **Proposition WC5**, ensemble runs with integrate-and-fire banks on Kaggle: the world bank never moves a daughter off its parent's centre-of-mass momentum, the world-particle banks always do; each realises its own generator up to the contact sink's bias (B1's covariance of `P^2` with `v` is -0.472 ± 0.003 against -0.501 on the mesh at `t = 0.25`, W's 0.002 ± 0.001 against 0); at `nu = 1e5` the world-particle banks carry 2.8 times the worlds and 36 times the unsigned `<P^2>` by `t = 1`, heating the ensemble's centre of mass while the signed marginal stays exact. A reader is a world, reading the chords of its own configuration along the directions in which its potential has terms.

*Ladder abstract — see the [full list](README.md#the-ladder).*

**Status.** Analysis note, step 25 of the ladder. Companion demo:
`src/demo_world_clocks.py` (Parts A, K, H, M, R, T, C, F). Prompted by a
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
integrator. Step 23's Proposition Q11 gave "every free body" a bank of
integrators, and Q-SP14 asked whether the bank could instead be the sea's
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
sharpened and one word is retired. Q11's "every free body carries a bank"
is the one-particle case of "every world carries one bank per bond"
(§7.2). The repository's *world-particle* was defined as a world (glossary,
ORIENTATION §2.3), which is right only for one particle; §1 redefines it.
The word *body* had two senses, an unpaired world in the ledger and a
physical particle in steps 21 and 24; it is retired in favour of *unpaired
world* and *particle* (WC-SP1, answered). This amendment, of October 2026,
also generalises the definitions from a line to three dimensions; the test
problems of §§6–7 stay one-dimensional.

**Leaves open.** WC-SP2 to WC-SP6 (§11), among them whether a live reading of
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
  joint cell. The ledger notes call an unpaired world a *body*, and steps 21
  and 24 also used *body* for a particle; the word is retired.

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

*Proof.* (a) and (c) are additivity of $`\mathcal L`$. For (b), Part A5
computes $`\mathcal L_1 + \mathcal L_2 - \mathcal L = (1 - \alpha_1 - \alpha_2)\thinspace v`$,
which is $`-v`$ at $`\alpha_i = 1`$ and $`0`$ when the shares sum to one.
For two worlds $`\mathbf X \pm \mathbf Y`$ on one momentum row the kinetic
windings cancel, so their clocks' misalignment winds at
$`-[V(\mathbf X + \mathbf Y) - V(\mathbf X - \mathbf Y)]/\hbar = -U/\hbar`$
(the one-particle statement is step 23's Proposition Q2: a same-row pair
winds at exactly $`U`$), and world-particle $`i`$'s clocks, which carry
$`\alpha_i v`$, at $`-\alpha_i U/\hbar`$. $`\square`$

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
S_W(\pmb\sigma) = \sum_{i<j} m_{ij}(\sigma_i - \sigma_j),
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

Proposition Q11 placed a bank on "every free body", and the one-particle
ladder could not tell a world-particle's bank from its world's. Theorem WC2 can:
the bank that realises the generator is indexed by the channels of the
world's joint chord, and for a pair potential those are exchanges, so it
must fire both particles at once. A bank owned by a world-particle fires one
world-particle and realises B1. Read in four dimensions, Q11 says: *every world
carries one bank per bond; each bank's integrators advance by the bond's
exchange rates and a firing moves both particles of the bond*. For one
degree of freedom the only bond is the particle's coupling to the external
potential, which is why "free body" (now *unpaired world*) and "world" gave the same bank.

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

## 10. Numerical verification

All outputs are from `src/demo_world_clocks.py`, seeds fixed in the code.
Parts A, K, H, M, R, T and F were run in the container, October 2026; Part C
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
- **Part F.** The two figures.

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
- Pöschl, G.; Teller, E. *Bemerkungen zur Quantenmechanik des anharmonischen
  Oszillators*, Z. Phys. **83** (1933) 143–151. The bound problem of §6.
- Eckart, C. *The penetration of a potential barrier by electrons*, Phys.
  Rev. **35** (1930) 1303–1309. The collision of §7.4.
