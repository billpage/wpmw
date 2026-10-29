# The chord action: one relation between worlds

**Definition WC7 of step 25 restated with a single generating relation, the
difference of two worlds' actions across a chord. Its part linear in the chord
is the classical force inside a world; the rest is the compensated residual,
which sets every event rate. Bonds, kernel tables and banks of integrators
then follow from it and are no longer postulates. A comparison with the
continuum of worlds of Schiff and Poirier, and with Hall, Deckert and
Wiseman's countable version, shows where each model puts the relation.**

---

## 0. Status and provenance

Supplement to step 25 ([`../analysis/world_clocks.md`](../analysis/world_clocks.md)
§8.1, "worlds suffice"). Companion demo: `src/demo_chord_action.py`, which
produces every number below. No figures. Slides:
[`../interactive/worlds_and_particles.html`](../interactive/worlds_and_particles.html)
([live page](https://billpage.github.io/wpmw/interactive/worlds_and_particles.html)),
the deck for the WPMW call of 27 October 2026 on worlds, particles and
world-particles; its slides 8 and 9 summarise this note, and the speaker
notes (key N) give the details.

Prompted by Bill Page's question of October 2026. In Schiff and Poirier's
continuum of worlds the interaction is defined by one term of one action, and
Hall, Deckert and Wiseman's countable model approximates it. WPMW's bonds,
chords and banks of integrators look much more complicated. Is there a way to
define the relation between worlds as simply?

There is. **Definition CA1** gives the relation: the action difference of two
worlds at the ends of a chord. **Proposition CA2** splits it into a part linear
in the chord and the compensated residual $`R`$. **Proposition CA3** shows that
the first is Liouville streaming and the second is the whole quantum part of
the generator. **Corollary CA4** recovers the clock relations of steps 23 and
25. **Proposition CA5** derives the bonds and the locality of WC8.1.
**Definition CA6** restates WC7 with one generating relation and lists what
is derived. **Remark CA7** demotes the banks of integrators to a realisation.
**Definition CA0** (§1) says what distinguishes a label from a coordinate,
which is what separates WPMW from Schiff and Poirier.

Labels CA0 to CA7 and open items CA-SP1 to CA-SP3 belong to this supplement.

**Corrects nothing.** WC4, WC8 and WC9 stand as stated; they become
consequences of CA1 instead of parts of the ontology. Q11's banks stand as the
realisation the demos use.

---

## 1. Labels and coordinates (Definition CA0)

A model of many worlds describes a motion: a map from the set of worlds and
time into configuration space $`Q = \mathbb R^{n}`$, with $`n = dN`$ for $`N`$
particles in $`d`$ dimensions.

**Definition CA0.** *A motion is a map $`\chi: B \times \mathbb R \to Q`$,
$`\mathbf x = \chi(c, t)`$, where $`B`$ is the set of worlds. The coordinates
of the domain, $`c \in B`$ and $`t`$, are the independent variables: the
action integrates over them and is not varied with respect to them. The
coordinates of the target, $`\mathbf x`$, are the dependent variables, which
the action varies. A label is a domain coordinate other than time: it is
constant along each world, $`\partial_t c = 0`$, and no physical quantity
depends on it except through the push-forward of its measure, the density.*

The distinction is one of role, not of kind. Inverting to the Eulerian
description $`c(\mathbf x, t)`$ makes $`\mathbf x`$ independent and the label a
field. Parametrised mechanics makes $`t`$ a dependent variable of a parameter
(Lanczos 1949). Within one formulation, however, the roles are fixed. This is
the material–spatial distinction of continuum mechanics (Marsden and Hughes
1983), with $`\partial\mathbf x/\partial c`$ the deformation gradient.

**Schiff and Poirier (2012).** The label is $`C`$, chosen to make the density
uniform: in one dimension $`C = \int_{-\infty}^{x_0}\rho_0\,dx_0'`$, so that
$`\rho = 1/x'`$ with $`x' = \partial x/\partial C`$, and in $`n`$ dimensions
$`\rho = \rho_C/\det J`$ with $`J^i{}_j = \partial x^i/\partial C^j`$. The
velocity is $`\dot x(C, t)`$, one value per position: the ensemble samples
$`\rho(\mathbf x)`$ and describes pure states. Their action, Eq. (11), is

```math
S_{\rm SP} = \iint dC\,dt\,\Big[\tfrac12 m\dot x^2 - V(x) - \frac{\hbar^2}{8m}\,\frac{x''^2}{x'^4}\Big],
```

written here in the form that differs from theirs by a total $`C`$-derivative
(Part S). With $`\rho = 1/x'`$ the last term is the Fisher information,
$`(\hbar^2/8m)\int(\partial_x\rho)^2/\rho\,dx`$ (Part S). Its Euler–Lagrange
equation is their Eq. (10):

```math
m\ddot x + \partial_x V + \frac{\hbar^2}{4m}\Big[\frac{x''''}{x'^4} - \frac{8x'''x''}{x'^5} + \frac{10x''^3}{x'^6}\Big] = 0 .
```

Worlds interact through derivatives in the label, which is a coupling on the
**domain**. Invariance under relabelling (volume-preserving maps of $`C`$) is a
gauge symmetry; its Noether charge is circulation (Salmon 1988; Padhye and
Morrison 1996), so a gradient velocity field, which their equivalence with the
Schrödinger equation requires, stays one.

**Hall, Deckert and Wiseman (2014).** These authors replace the label
derivatives by differences between neighbouring worlds,
$`U_N = (\hbar^2/8m)\sum_n[1/(x_{n+1} - x_n) - 1/(x_n - x_{n-1})]^2`$. For
worlds at the quantiles of a unit Gaussian, $`U_N/N`$ is $`0.1090`$, $`0.1225`$,
$`0.1247`$, $`0.12496`$ at $`N = 10^2`$ to $`10^5`$ against the Fisher value
$`1/8`$ (Part S). What a countable version needs is the neighbour graph of the
label lattice, which in one dimension can be recovered from the order of the
positions. In more than one dimension the positions no longer determine the
label's neighbours, which is why sorting does not generalise (step 25 §2).

**WPMW.** The domain is time and a discrete world index with persistent
identity (WC6). There is no label manifold, and no smooth labelling or
adjacency is needed. The relation between worlds is a function of **target**
coordinates only: two points $`\mathbf X \pm \mathbf Y`$ of the same
configuration space. What it needs of $`Q`$ is midpoints and differences.

| | Schiff–Poirier | Hall–Deckert–Wiseman | WPMW |
|---|---|---|---|
| a world | label $`C`$, path $`\mathbf x(C, t)`$ | index $`n`$, path $`\mathbf x_n(t)`$ | point $`(\mathbf x, \mathbf p)`$, sign, clock |
| identity | the label, smooth | the index | persistent; a genealogy (CX5) |
| momentum | $`m\dot{\mathbf x}`$, one per position | the same | independent |
| relation lives on | the domain: $`\partial_C`$ | the domain: label neighbours | the target: chords |
| what it is | Fisher term, nonlinear in $`\rho`$ | its difference form | $`R(\mathbf X, \mathbf Y)`$, linear in $`W`$ |
| states | pure, positive | pure, positive | any, signed |
| finite-ensemble price | — | $`O(1/N)`$ bias | cancellation (PC2, K1) |

## 2. Definition CA1 and Proposition CA2: the chord action and its split

Let worlds $`w`$ and $`w'`$ be at configurations $`\mathbf x, \mathbf x' \in Q`$,
with Lagrangian $`\mathcal L = \sum_i \tfrac12 m_i\lvert\dot{\mathbf x}_i\rvert^2 - V(\mathbf x)`$.
Write the chord's midpoint and half-length

```math
\mathbf X = \tfrac12(\mathbf x + \mathbf x'),\qquad \mathbf Y = \tfrac12(\mathbf x - \mathbf x') .
```

**Definition CA1 (the chord action).** *The chord action of two worlds is the
difference of their actions, $`\Delta S_{ww'} = \hbar(\theta^w - \theta^{w'})`$
along their histories, with Lagrangian difference*

```math
\Delta\mathcal L(\mathbf X, \mathbf Y) = \mathcal L(\mathbf X + \mathbf Y) - \mathcal L(\mathbf X - \mathbf Y)
= 2\sum_i m_i\,\dot{\mathbf X}_i\cdot\dot{\mathbf Y}_i - \big[V(\mathbf X + \mathbf Y) - V(\mathbf X - \mathbf Y)\big].
```

By P1 each clock winds at $`\mathcal L/\hbar`$ along its own world, so the
chord action is the clocks' misalignment before the lever is added.

**Proposition CA2 (the split).** *For any $`V`$,*

```math
\Delta\mathcal L = \Delta\mathcal L_{\rm cl} - R(\mathbf X, \mathbf Y),\qquad
\Delta\mathcal L_{\rm cl} = 2\sum_i m_i\,\dot{\mathbf X}_i\cdot\dot{\mathbf Y}_i - 2\,\mathbf Y\cdot\nabla V(\mathbf X),
```

*with $`R(\mathbf X, \mathbf Y) = V(\mathbf X + \mathbf Y) - V(\mathbf X - \mathbf Y) - 2\mathbf Y\cdot\nabla V(\mathbf X)`$
the compensated residual of step 14 and WC9.*

*(a) $`\Delta\mathcal L_{\rm cl}`$ is linear in $`\mathbf Y`$, and its
Euler–Lagrange equation in $`\mathbf Y_i`$ is Newton's law for the midpoint,
$`m_i\ddot{\mathbf X}_i = -\nabla_i V(\mathbf X)`$.*

*(b) $`R`$ is odd in $`\mathbf Y`$ and $`O(\lvert\mathbf Y\rvert^3)`$; it
vanishes identically if and only if $`V`$ is at most quadratic.*

*Proof.* The kinetic difference is $`\tfrac12 m(\dot X + \dot Y)^2 - \tfrac12 m(\dot X - \dot Y)^2 = 2m\dot X\dot Y`$
for each coordinate; the rest is the definition of $`R`$. For (a), the
Euler–Lagrange equation of $`2m_i\dot{\mathbf X}_i\cdot\dot{\mathbf Y}_i - 2\mathbf Y_i\cdot\nabla_i V`$
in $`\mathbf Y_i`$ is $`-2m_i\ddot{\mathbf X}_i - 2\nabla_i V = 0`$. Part A
checks both for three particles with a general $`V`$ (residuals $`0`$). For
(b), $`V(\mathbf X + \mathbf Y) - V(\mathbf X - \mathbf Y)`$ is odd and its
linear term is $`2\mathbf Y\cdot\nabla V`$; the remainder starts at the third
derivative, which vanishes for all $`\mathbf X`$ only for a quadratic $`V`$.
Part A checks that the orders $`0`$ to $`2`$ vanish for $`v = -6\,{\rm sech}^2 x`$. $`\square`$

The two words of step 25's slogan are the two terms. **Classical within a
world** is $`\Delta\mathcal L_{\rm cl}`$: it carries the force and nothing
else. **Quantum between worlds** is $`R`$: it carries everything the force
does not, and it is the only place two worlds' configurations enter
together.

## 3. Proposition CA3: from the chord action to the generator

The potential's part of the von Neumann equation is a phase rate on each
element of the density matrix, $`\partial_t\rho(\mathbf x, \mathbf x') = -(i/\hbar)[V(\mathbf x) - V(\mathbf x')]\rho`$:
the phase of $`\rho(\mathbf x, \mathbf x')`$ advances at the potential part
of $`\Delta\mathcal L/\hbar`$, as $`\theta^w - \theta^{w'}`$ does by CA1
(supplement [`chord_reading.md`](chord_reading.md) §5, which orients the
chord the other way). In the density matrix
$`\mathbf Y`$ is the off-diagonal coordinate, and nothing of the state is
at its ends; in the sea it is the half-length of a chord, whose ends are two
aligned pairs (step 22; chord_reading §5, whose table sets the two side by
side; glossary, *chord* and *separation*). The Wigner transform

```math
W(\mathbf X, \mathbf P) = (\pi\hbar)^{-n}\int d\mathbf Y\,\rho(\mathbf X + \mathbf Y, \mathbf X - \mathbf Y)\,e^{-2i\mathbf P\cdot\mathbf Y/\hbar}
```

puts the midpoint at the world's position and makes $`\mathbf P`$ conjugate to
the chord.

**Proposition CA3.** *Under the split of CA2 the Wigner evolution is*

```math
\partial_t W = -\sum_i \frac{\mathbf P_i}{m_i}\cdot\nabla_{\mathbf X_i}W + \nabla V(\mathbf X)\cdot\nabla_{\mathbf P}W
+ \int d\mathbf Q\,\Theta(\mathbf X, \mathbf Q)\,W(\mathbf X, \mathbf P - \mathbf Q),
```

```math
\Theta(\mathbf X, \mathbf Q) = -\frac{1}{\hbar}\,(\pi\hbar)^{-n}\int d\mathbf Y\,R(\mathbf X, \mathbf Y)\,\sin\frac{2\mathbf Q\cdot\mathbf Y}{\hbar} .
```

*The linear part of the chord action gives the Liouville term exactly. The
residual gives a real kernel, odd in $`\mathbf Q`$, which is the whole
quantum part of the generator.*

*Proof.* With $`\rho(\mathbf X + \mathbf Y, \mathbf X - \mathbf Y) = \int d\mathbf P'\,W(\mathbf X, \mathbf P')\,e^{2i\mathbf P'\cdot\mathbf Y/\hbar}`$,
a phase rate $`-(i/\hbar)f(\mathbf X, \mathbf Y)`$ becomes a convolution in
$`\mathbf P`$ with kernel $`-(i/\hbar)(\pi\hbar)^{-n}\int f\,e^{-2i\mathbf Q\cdot\mathbf Y/\hbar}d\mathbf Y`$.
For $`f = 2\mathbf Y\cdot\nabla V`$, using $`\mathbf Y e^{-2i\mathbf Q\cdot\mathbf Y/\hbar} = (i\hbar/2)\nabla_{\mathbf Q}e^{-2i\mathbf Q\cdot\mathbf Y/\hbar}`$,
the kernel is $`\nabla V\cdot\nabla\delta(\mathbf Q)`$, which gives
$`\nabla V\cdot\nabla_{\mathbf P}W`$. For $`f = R`$, oddness turns the
exponential into $`-i\sin`$. The kinetic part gives the streaming term as
usual. $`\square`$

Part B checks the kernel. For $`V = x^4`$, where $`R = 8XY^3`$ and the Moyal
series has one term, $`-(\hbar^2/24)V'''\partial_P^3W`$, the two agree to
$`7.9\times10^{-11}`$. For $`V = -6\,{\rm sech}^2x`$ and a Gaussian of width
$`3`$ in $`P`$, the Moyal series truncated at $`\hbar^2`$, $`\hbar^4`$,
$`\hbar^6`$ misses the kernel by $`7.4\times10^{-2}`$, $`5.3\times10^{-3}`$,
$`3.1\times10^{-4}`$. At width $`0.7`$, where the chords that matter are long,
it misses by $`1.2`$, $`1.6`$, $`2.0`$: the series is asymptotic. The kernel
is the generator; the series is only its short-chord expansion.

Discretising $`\mathbf Q`$ on the channels $`\pmb\xi_q`$, with the horizon at
the reach (G5), gives the kernel tables $`K_q(\mathbf X)`$ of the
specification. A world at $`\mathbf P'`$ contributes $`+\Theta`$ at
$`\mathbf P' + \mathbf Q`$ and $`-\Theta`$ at $`\mathbf P' - \mathbf Q`$: that
dipole is what an emission realises (specification §5.3).

## 4. Corollary CA4: what the clocks see

**Corollary CA4.** *(a) Two worlds on one momentum row have
$`\dot{\mathbf Y} = 0`$, and their misalignment, lever included, winds at*

```math
\hbar\,\dot\mu = -\big[V(\mathbf X + \mathbf Y) - V(\mathbf X - \mathbf Y)\big] .
```

*(b) A reader at the midpoint whose own momentum follows the force sees
$`-R(\mathbf X, \mathbf Y)`$: its own streaming has absorbed
$`\Delta\mathcal L_{\rm cl}`$.*

*Proof.* (a) On a common row the velocities are equal, so the kinetic
windings cancel and the lever term $`\mathbf P\cdot(\mathbf x' - \mathbf x)/\hbar`$
is constant; what is left of $`\Delta\mathcal L/\hbar`$ is the potential
difference (Part C, residual $`0`$). (b) is CA2 with the linear part moved into
the reader's motion. $`\square`$

(a) is Q2(b) of step 23 and WC8(a) of step 25. Under (S′) aligned pairs are
force-blind, so two dark pairs on one row keep their separation exactly and
(a) holds without approximation. Under (S) the force would fan each row out
and (a) would fail at $`V''`$ (chord_reading §2). (b) is the compensated split
seen by a reader. So the force-blind sea is what makes the chord action of two
aligned pairs a pure potential difference, and the reader is what makes the
split.

## 5. Proposition CA5: bonds are derived

**Proposition CA5.** *(a) If $`V = \sum_i V_i(\mathbf x_i) + \sum_{i < j}v_{ij}(\mathbf x_i - \mathbf x_j)`$
plus a quadratic part, then*

```math
R(\mathbf X, \mathbf Y) = \sum_i R_{V_i}(\mathbf X_i, \mathbf Y_i) + \sum_{i < j}\rho_{v_{ij}}(\mathbf r_{ij}, \pmb\eta_{ij}),
\qquad \mathbf r_{ij} = \mathbf X_i - \mathbf X_j,\quad \pmb\eta_{ij} = \mathbf Y_i - \mathbf Y_j,
```

*with $`\rho_v(\mathbf r, \pmb\eta) = v(\mathbf r + \pmb\eta) - v(\mathbf r - \pmb\eta) - 2\pmb\eta\cdot\nabla v(\mathbf r)`$.
So $`\Theta`$ is supported on transfers along single particles and along the
bond directions $`\mathbf e_i - \mathbf e_j`$; for pair potentials alone, $`R`$
is blind to the centre-of-mass direction of the chord, and every event
conserves the total momentum.*

*(b) A chord that moves only the particles of a set $`A`$ reads exactly the
terms of $`V`$ that involve a particle of $`A`$.*

*Proof.* $`R`$ is linear in $`V`$ and the quadratic part drops out by CA2(b).
A pair term depends on $`\mathbf Y`$ only through $`\pmb\eta_{ij}`$, and a
function of $`\pmb\eta_{ij}`$ has a Fourier transform supported along
$`\mathbf e_i - \mathbf e_j`$; its gradient summed over the particles is zero.
A term with no argument in $`A`$ takes the same value at both ends of the
chord. Part D checks, for three particles, that $`\sum_i\partial_{Y_i}R = 0`$
for pair terms and a trap, and that a chord on particle 1 alone reads
$`U(x_1)`$, $`v(x_1 - x_2)`$, $`v(x_1 - x_3)`$ and $`w(x_1, x_2, x_3)`$ and
nothing else. $`\square`$

(a) is WC9 and WC4, and (b) is WC8.1. Here they follow from the chord action
alone. A bond is not a separate object of the ontology: it is a term of
$`V`$, and the directions it opens to events are the support of $`\Theta`$. A
kernel table per bond is the statement that the bonds of identical particles
share one function $`\rho_v`$.

## 6. Definition CA6: the world ontology with one generating relation

**Definition CA6 (Definition WC7 restated).** *Under (S′) with persistent
identity:*

*(i) Carriers. Worlds, each with a point $`(\mathbf x, \mathbf p) \in \mathbb R^{2n}`$,
a sign $`\varepsilon = \pm1`$ and a clock $`\theta`$ that winds at
$`\mathcal L/\hbar`$ (P1). Everything else about one world is a function of
its point (WC7(ii)).*

*(ii) The generating relation. The chord action of two worlds (CA1). Its
linear part is the classical motion of each world (CA2(a), CA3); its residual
sets the signed rate at which a world at the chord's midpoint parents events,
$`\Theta(\mathbf X, \mathbf Q)`$ (CA3).*

*(iii) Realisation. An event with persistent identity is an emission, which
ionises the aligned pair at the parent's point, or an absorption, which joins
two unpaired worlds into one (WC6, specification §5.3). Opposite signs at one
point cancel in every observable.*

*(iv) The dark sector. Aligned pairs are force-blind and change momentum only
in pair collisions (step 24); dark catalysis exchanges phase (L5, Q6). These
rules move nothing in $`W`$.*

*(v) Admissibility. The ensemble's $`W`$ is the Wigner function of a positive
operator.*

*Derived, not postulated:* world-particles, centre of mass, Jacobi and sector
coordinates and sector clocks (WC1(c), WC7(ii)); bonds, their kernel tables,
the dark centre of mass and locality (CA5); the reader's split (CA4(b)); and
the banks of integrators (Remark CA7).

WC7's table of relations becomes:

| WC7 relation | worlds | in CA6 |
|---|---|---|
| chord | 2 | the generating relation (ii) |
| reading | many | $`\Theta`$ at the reader's point: the sum over chords of (ii) |
| emission | 3 | realisation (iii) |
| absorption | 3 | realisation (iii) |
| pair collision | 4 | dark sector (iv) |
| dark catalysis | 4 | dark sector (iv) |
| admissibility | all | constraint (v) |

**Remark CA7 (banks are a realisation).** *A bank integrates
$`\lvert\Theta\rvert`$ along its world's history, channel by channel, and fires
at each integer (Q11). A Poisson process driven by the same rate has the same
mean count and so realises the same generator; the two differ in their
variance.* Part E drives both with one rate along one history, $`\int K\,dt = 6.9556`$:
integrate-and-fire gives mean $`6.9562`$ and variance $`0.042`$, Poisson
$`6.9494`$ and $`6.962`$. A bank is therefore a choice of how to time events,
not a relation between worlds. Its being per bond (WC4) is CA5(a): one bank
per term of $`R`$.

## 7. Where the chord's ends are

CA3 needs only $`R`$, a function of $`V`$ at the two configurations
$`\mathbf X \pm \mathbf Y`$: the generator is blind to whether worlds occupy
them (chord_reading §5, "the sea reads the operator"). In the ontology of
steps 22 and 23 the ends are aligned pairs, and the sea's density enters only as
the measure of chords, $`d\mathbf Y/(\pi\hbar)^n`$. Whether the ends must be
occupied is the one choice CA6 leaves open (CA-SP1). Schiff and Poirier have
no such choice, because their neighbours in $`C`$ are always occupied. The
live reading of Q12, which works open loop and fails closed loop at the
physical density, is the test case.

## 8. Relation to path integrals

The chord action is the exponent of Feynman and Vernon's double path integral
(1963), $`\rho \propto \iint\mathcal D\mathbf x\,\mathcal D\mathbf x'\,e^{i(S[\mathbf x] - S[\mathbf x'])/\hbar}`$,
in the classical and quantum coordinates $`\mathbf X`$, $`\mathbf Y`$ of the
Keldysh technique (Kamenev 2011). Integrating $`\mathbf Y`$ against the linear
part enforces the classical equation of motion on $`\mathbf X`$; the terms of
higher order in $`\mathbf Y`$ are the quantum vertices. The bra branch
$`\mathbf x'`$ evolves under $`-H`$: it is the time-reversed copy of the
system, the backward branch of the closed time path, and with $`\mathbf x`$
it forms a chord, not a second particle or a world of its own. Thermo-field
dynamics (Takahashi and Umezawa 1975) doubles the system with the same
generator, $`\hat H = H - \tilde H`$, but its state vector is
$`\rho^{1/2}`$, a purification, so its expectations are quadratic in the
state. WPMW, like Liouville-space methods and the Keldysh technique, is linear
in $`\rho`$; Nakamura and Yamanaka (2012) map the Liouville-space form onto
non-equilibrium thermo-field dynamics. Marinov (1991) wrote
the Wigner propagator as a phase-space path integral with the difference of the
Hamiltonian at shifted arguments. Polkovnikov (2010) dropped the residual to
obtain the truncated Wigner approximation and represented the corrections as
stochastic quantum jumps with signed weights. WPMW is such an unravelling with
two specific choices: the classical part is enforced on every world exactly,
and the jumps are realised with persistent identity, the parent unkinked
(WC6(d)).

## 9. Numerical verification

`WPMW_OUTPUT=<dir> PYTHONPATH=src python3 -u src/demo_chord_action.py`
(about 20 s), captured:

```text
=== Part S: Schiff-Poirier and Hall-Deckert-Wiseman ===
  Q_SP - Q_Fisher - d/dC[(hbar^2/4m) x''/x'^3]  = 0
  Euler-Lagrange of the Fisher form + SP Eq. (10) = 0
  (d_x rho)^2/rho dx/dC - x''^2/x'^4, rho = 1/x'  = 0
  HDW U_N / N for a Gaussian of unit width (hbar = m = 1); Fisher value 1/8:
    N =    100:  U_N/N = 0.108974   deficit 1.60e-02
    N =   1000:  U_N/N = 0.122529   deficit 2.47e-03
    N =  10000:  U_N/N = 0.124664   deficit 3.36e-04
    N = 100000:  U_N/N = 0.124957   deficit 4.26e-05

=== Part A: Proposition CA2, the split of the chord action ===
  N = 3, general V:  dL - [linear part] + R = 0
  Euler-Lagrange in Y_i of the linear part = -2[m_i X_i'' + d_i V(X)]: residuals [0, 0, 0]
  R for v = -6 sech^2: coefficients of lambda^0, lambda^1, lambda^2 = [0, 0, 0]

=== Part B: Proposition CA3, the residual kernel is the Moyal residual ===
  V = x^4, x0 = 0.4: max |kernel - Moyal| / max |Moyal| = 7.89e-11
  V = -6 sech^2 x, momentum width 3.0: relative error of the Moyal series truncated at hbar^2, hbar^4, hbar^6: 7.36e-02, 5.32e-03, 3.05e-04
  V = -6 sech^2 x, momentum width 0.7: relative error of the Moyal series truncated at hbar^2, hbar^4, hbar^6: 1.15e+00, 1.58e+00, 1.99e+00

=== Part C: Corollary CA4, misalignment on a common momentum row ===
  hbar dmu/dt + [V(X+Y) - V(X-Y)] = 0

=== Part D: Proposition CA5, bonds and locality are derived ===
  pair potential + trap: (sum_i d/dY_i) R = 0   (centre-of-mass direction dark)
  a chord on particle 1 alone reads:
     U(-eta + x1,)
     U(eta + x1,)
     v(-eta + x1 - x2,)
     v(-eta + x1 - x3,)
     v(eta + x1 - x2,)
     v(eta + x1 - x3,)
     w(-eta + x1, x2, x3)
     w(eta + x1, x2, x3)

=== Part E: Remark CA7, banks are a realisation of the rate ===
  integral of the rate = 6.955564
  integrate-and-fire: mean 6.9562, variance 0.0419
  Poisson:            mean 6.9494, variance 6.9621
```

## 10. Open items

- **CA-SP1.** Must the ends of a chord be occupied? CA3 needs only the
  configurations. If they need not be, the sea's role in the generator is the
  measure $`d\mathbf Y/(\pi\hbar)^n`$ alone, and its density is a realisation
  parameter; if they must, the live reading of Q12 and WC-SP2 decides whether
  a sea can supply them.
- **CA-SP2.** Is there an action for the whole signed ensemble, as $`S_{\rm SP}`$
  is for Schiff and Poirier's, from which CA3 follows by variation rather than
  by Wigner transform? The chord action is a functional of two worlds; the
  generator is linear in $`W`$, so such an action would be bilinear, as
  Feynman and Vernon's is in the two paths.
- **CA-SP3.** The dark-sector rules (iv) are not derived from the chord
  action. Are they needed beyond keeping the sea uniform (V9, V10), and does
  CA-SP1 decide that?

## 11. Sources

- Step 25, Definition WC7, Propositions WC8 and WC9, Corollary WC8.1:
  [`../analysis/world_clocks.md`](../analysis/world_clocks.md) §8.1 and §8.2.
- Step 23, Propositions Q2 and Q11:
  [`../analysis/force_blind_sea.md`](../analysis/force_blind_sea.md).
- The kernel from the commutator, and the sea reading the operator:
  [`chord_reading.md`](chord_reading.md) §§2 and 5.
- Schiff, J.; Poirier, B. — "Communication: Quantum Mechanics without
  Wavefunctions." *J. Chem. Phys.* **136**, 031102 (2012). Eqs. (9)–(11),
  (17)–(20). https://doi.org/10.1063/1.3680558 (arXiv:1201.2382)
- Holland, P. — *Ann. Phys.* **315**, 505–531 (2005).
  https://doi.org/10.1016/j.aop.2004.09.008
- Hall, M. J. W.; Deckert, D.-A.; Wiseman, H. M. — *Phys. Rev. X* **4**,
  041013 (2014). https://doi.org/10.1103/PhysRevX.4.041013
- Feynman, R. P.; Vernon, F. L. — "The theory of a general quantum system
  interacting with a linear dissipative system." *Ann. Phys.* **24**,
  118–173 (1963). https://doi.org/10.1016/0003-4916(63)90068-X
- Marinov, M. S. — "A new type of phase-space path integral." *Phys. Lett. A*
  **153**, 5–11 (1991).
- Polkovnikov, A. — "Phase space representation of quantum dynamics." *Ann.
  Phys.* **325**, 1790–1852 (2010), §III.1.3.
  https://arxiv.org/abs/0905.3384
- Kamenev, A. — *Field Theory of Non-Equilibrium Systems.* Cambridge (2011).
- Takahashi, Y.; Umezawa, H. — "Thermo field dynamics." *Collective
  Phenomena* **2**, 55–80 (1975); reprinted *Int. J. Mod. Phys. B* **10**,
  1755–1805 (1996).
- Nakamura, Y.; Yamanaka, Y. — "From superoperator formalism to
  nonequilibrium Thermo Field Dynamics." (2012).
  https://arxiv.org/abs/1211.4306
- Marsden, J. E.; Hughes, T. J. R. — *Mathematical Foundations of
  Elasticity.* Prentice-Hall (1983).
- Salmon, R. — "Hamiltonian fluid mechanics." *Annu. Rev. Fluid Mech.*
  **20**, 225–256 (1988). Padhye, N.; Morrison, P. J. — "Fluid element
  relabeling symmetry." *Phys. Lett. A* **219**, 287–292 (1996).
- Lanczos, C. — *The Variational Principles of Mechanics.* University of
  Toronto Press (1949).
