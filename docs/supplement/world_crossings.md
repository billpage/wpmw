# Crossings of worlds in configuration space-time

**What it means for worlds to cross when there are three or more particles, why
chance crossings almost never happen, what the worlds that share an event's
vertex offer instead, and why neither replaces a sink: the population control
for [problem G](../analysis/world_clocks.md#test-problems-the-moshinsky-harmonic-model-q-and-the-smoothed-contact-interaction-g) is K1 resampling (Proposition PC6).**

---

## 0. Status and provenance

Supplement to step 26 ([`../analysis/population_control.md`](../analysis/population_control.md)),
prompted by Bill Page's question: can crossings of worlds within a cell and a
time step, rather than proximity in phase space, decide when an event is
realised by absorption, and how is a crossing defined when there are three or
more particles? Companion demo: `src/demo_world_crossings.py`, which produces
every number below and both figures.

It inherits the one-dimensional picture of step 23
([`../analysis/force_blind_sea.md`](../analysis/force_blind_sea.md) §9):
Proposition Q9(b) reads an absorption as two unpaired worlds on rows
$`p \pm \xi_q`$ crossing at the parent and binding at its row, and
Proposition Q10 shows that crossings can *realise* events but cannot
*trigger* them, because a crossing rate is mass action and the observable
would not close. Nothing here changes Q10: the parent's reading triggers the
event, and the question is only how the event is realised. It also uses
Proposition WC6 (a world's history: momenta kink only at events, positions
are continuous), Proposition WC9 (events act on one bond at a time),
Theorem S6 of step 16 (an event is wholly emissive or wholly absorptive) and,
from step 26, Proposition PC5 (absorption with a tolerance in phase space),
Proposition PC6 (K1 resampling) and open item PC-SP4 (exact absorption off a
mesh).

New here: a definition of a crossing for $`N`$ particles (CX1), the
equivalence of a delayed realisation with a translated one (CX2), the
codimension of a chance crossing (CX3), the measured yield of crossings
against PC5's box (CX4), partners selected by their record of a shared vertex
and the exact identity of a paired kin absorption with an undo (CX5), and
when the flow brings kin back together (CX6). Labels CX1 to CX6 and open items
CX-SP1 to CX-SP4 belong to this supplement; CX-SP2 is withdrawn (§11).

**The upshot.** Neither chance nor descent supplies absorption partners often
enough to stand in for a sink. Theorem S6 confines kin absorption to the
paired case, which is exactly an undo of an earlier event, and that happens
for under one event in a hundred. The population control for [problem G](../analysis/world_clocks.md#test-problems-the-moshinsky-harmonic-model-q-and-the-smoothed-contact-interaction-g) stays
K1 resampling (PC6); crossings belong to the ontology's picture of an
absorption, not to the numerics.

## 1. Setting

$`N`$ particles of mass $`m`$ in $`d`$ dimensions. A world $`w`$ is a point
$`z_w = (\mathbf x_w, \mathbf p_w) \in \mathbb R^{2dN}`$ with a sign
$`\varepsilon_w = \pm 1`$. The centre of mass separates exactly and is never
read, so the *relative configuration space* has dimension

```math
n = d\,(N - 1),
```

and the relative phase space $`2n`$. Between events every world follows the
classical flow $`\Phi_t`$; its position is continuous and its momentum kinks
only at events (WC6).

An **event** of parent $`A`$ at time $`t`$ on bond $`b = (i, j)`$ in channel
$`q`$ transfers $`+\mathbf k`$ to $`\mathbf p_i`$ and $`-\mathbf k`$ to
$`\mathbf p_j`$, with $`\lvert\mathbf k\rvert = \xi_q = q\,\Delta p`$ (on a line,
$`k = \xi_q`$). Write $`\mathbf K \in \mathbb R^{dN}`$ for that transfer. The
event's sign is $`\tau = \varepsilon_A\,\mathrm{sgn}\,K_q(r_b)`$, and its two
**targets** are

```math
z_\pm = (\mathbf x_A,\ \mathbf p_A \pm \mathbf K), \qquad \text{with signs } \pm\tau .
```

*Emission* creates worlds at the targets. *Absorption* removes a world of sign
$`-\tau`$ at $`z_+`$ and one of sign $`+\tau`$ at $`z_-`$, the **partners**
$`B_+`$ and $`B_-`$. Both change $`W`$ by $`+\tau`$ at $`z_+`$ and $`-\tau`$ at
$`z_-`$ (glossary, *absorptive, emissive*). In the ontology the absorbed
partners kink onto the parent's momentum and continue as an aligned pair.

**Configuration space-time** is $`\mathbb R^{dN} \times \mathbb R`$, or
$`\mathbb R^n \times \mathbb R`$ after removing the centre of mass. Each world
traces a world-line $`t \mapsto \mathbf x_w(t)`$ in it.

## 2. Definition CX1: the absorption vertex and the crossing

**Definition CX1.** *(a) An* exact absorption vertex *of the event
$`(A, b, q)`$ is a point $`(\mathbf x^*, t^*)`$ of configuration space-time
where the parent's world-line meets the world-lines of two partners whose
momenta there are exactly the targets':*

```math
\mathbf x_{B_\pm}(t^*) = \mathbf x_A(t^*), \qquad \mathbf p_{B_\pm}(t^*) = \mathbf p_A(t^*) \pm \mathbf K .
```

*(b) A world $`B`$* crosses *the event during $`[t, t + \Delta t]`$ if such a
$`t^*`$ lies in that interval. (c) It* nearly crosses *it with impact
parameter $`\beta`$ and momentum tolerance $`\delta p`$ if, over the interval,
the closest approach of $`\mathbf x_B - \mathbf x_A`$ to zero is at most
$`\beta`$ and $`\mathbf p_B`$ is within $`\delta p`$ of a target. The set of
such $`B`$ is the* crossing tube *(figure, panel b). The demo calls
$`\beta`$ `b`.*

**The kinematics of the vertex.** A partner at $`\mathbf p_A \pm \mathbf K`$
moves relative to the parent at velocity $`\pm\mathbf K/m`$. Only particles
$`i`$ and $`j`$ have different velocities, by $`\pm\mathbf k/m`$ and
$`\mp\mathbf k/m`$, so the relative motion lies entirely along the **bond
direction** in configuration space, and every other coordinate (the other
particles, the pair's centre of mass, the total centre of mass) is the same
for parent and partner. The two partners approach the parent along the bond
from opposite sides and meet it at the vertex (panel a). An emission is the
time reverse: the daughters leave the parent's point along the bond. This is
Q9(b) in $`N`$ particles: there the partners on rows $`p \pm \xi_q`$ approach
at relative velocity $`2\xi_q/m`$ and bind at the parent's row.

## 3. Lemma CX2: a delayed realisation is a translated one

**Lemma CX2.** *For a parent at $`z_A = (\mathbf x, \mathbf p)`$ and a daughter
at $`z_A + (0, \mathbf K)`$, the classical flow over a time $`\tau`$ gives*

```math
\Phi_\tau(\mathbf x, \mathbf p + \mathbf K) = \Phi_\tau(\mathbf x, \mathbf p)
 + \bigl(\tau\mathbf K/m + O(\tau^3),\ \mathbf K + O(\tau^2)\bigr).
```

*So realising an event late, at $`t + \tau`$, at the parent's position then,
differs from realising it at $`t`$ and letting the daughters flow, by a
translation of each daughter by $`\mp\tau\mathbf K/m`$ along the bond
direction, with momentum error $`O(\tau^2)`$. The centre of mass is not
moved.*

*Proof.* Taylor in $`\tau`$. Positions: $`\mathbf x(\tau) = \mathbf x + \tau\mathbf p/m +
\tfrac12\tau^2\mathbf F(\mathbf x)/m + \tfrac16\tau^3 (D\mathbf F\,\mathbf p)/m^2 + O(\tau^4)`$, where
$`\mathbf F = -\nabla V`$ and $`D\mathbf F`$ is its Jacobian. The daughter and
the parent start at the same $`\mathbf x`$, so their forces agree at
$`\tau = 0`$ and the difference is $`\tau\mathbf K/m + \tfrac16\tau^3 D\mathbf F\,\mathbf K/m^2`$.
Momenta: the difference is $`\mathbf K + \int_0^\tau [\mathbf F(\mathbf x_d) - \mathbf F(\mathbf x_a)]\,ds =
\mathbf K + \tfrac12\tau^2 D\mathbf F\,\mathbf K/m + O(\tau^3)`$. $`\mathbf K`$ has no
centre-of-mass component. $`\square`$

*Verification* (Part A, [problem G](../analysis/world_clocks.md#test-problems-the-moshinsky-harmonic-model-q-and-the-smoothed-contact-interaction-g) of step 25, 4000 random parents and
channels, flow integrated with 400 Verlet steps per interval):

| $`\tau`$ | max position residual $`/\tau^2`$ | max momentum residual $`/\tau^2`$ | centre-of-mass shift |
|---|---|---|---|
| 0.08 | 0.510 | 18.3 | $`1.2\times10^{-14}`$ |
| 0.04 | 0.294 | 21.6 | $`1.2\times10^{-14}`$ |
| 0.02 | 0.153 | 22.8 | $`7.9\times10^{-15}`$ |
| 0.01 | 0.077 | 23.0 | $`1.2\times10^{-14}`$ |

The position residual over $`\tau^2`$ halves with $`\tau`$ (third order) and
the momentum residual over $`\tau^2`$ settles at a constant (second order).

**What it means.** A tolerance in *time* is a tolerance in *position along the
bond*, of length $`\lvert\mathbf K\rvert\tau/m`$. Waiting for a crossing is not free: it
is PC5's translation error, confined to one direction. There is one
exception. The demos fire events once per step, after the Verlet step, so an
event's time is already uncertain by $`\Delta t`$. Realising an event at a
crossing *inside* the step therefore adds no error of lower order than the
scheme already has; a crossing later than the step does.

## 4. Lemma CX3: how often worlds cross by chance

**Lemma CX3.** *Let the partners have a bounded continuous density near a
target, and fix the step $`\Delta t`$.*

*(a) With continuous momenta, an exact crossing in the step has probability
zero for every $`n \ge 1`$: it imposes $`2n`$ conditions (position and
momentum) on one free parameter $`t^*`$.*

*(b) If momenta live on a lattice that contains the target (PC-SP4), it
imposes $`n`$ conditions on one parameter, and has positive probability if
and only if $`n = 1`$.*

*(c) For a partner of the right momentum, moving at speed $`v`$ relative to
the parent, the probability of passing within $`\beta`$ during the step is
proportional to the volume of the tube, a segment of length $`L = v\,\Delta t`$
thickened by $`\beta`$:*

```math
P(\text{miss} \le \beta) \;\propto\; c_{n-1}\,\beta^{\,n-1} L + c_n\,\beta^{\,n},
```

*with $`c_k\beta^k`$ the volume of a $`k`$-ball of radius $`\beta`$. For
$`\beta \ll L`$ it scales as $`\beta^{\,n-1}`$.*

*Proof.* (a), (b) Two curves in $`\mathbb R^{n} \times \mathbb R`$ meet
generically only if $`n = 1`$; momentum equalities add $`n`$ more conditions
unless a lattice makes them hold with positive probability. (c) The set of
initial offsets that pass within $`\beta`$ of the parent in the step is the
Minkowski sum of the segment $`-\mathbf v[0, \Delta t]`$ and a $`\beta`$-ball.
$`\square`$

| system | $`n = d(N-1)`$ | exact crossing with a momentum lattice |
|---|---|---|
| two particles on a line | 1 | positive probability |
| three on a line | 2 | probability zero |
| four on a line | 3 | probability zero |
| two in space | 3 | probability zero |
| three in space | 6 | probability zero |

*Verification* (Part B, $`4\times10^6`$ random offsets in a box of side 0.5,
$`L = 0.2`$): for $`\beta = 0.002`$ to $`0.016`$ the measured probabilities
match the tube volume to within the sampling error, the exact crossing
($`\beta = 0`$) has probability $`0.4000`$ at $`n = 1`$ (predicted $`0.4000`$)
and $`0`$ at $`n = 2, 3`$, and the slopes $`d\log P/d\log\beta`$ are
$`0.06`$, $`1.05`$ and $`2.06`$ for $`n = 1, 2, 3`$.

So for three or more particles a crossing, taken literally, does not happen
by chance. What can be defined and used is a near crossing, and its tolerance
$`\beta`$ is spent in $`n - 1`$ directions transverse to the bond.

## 5. Proposition CX4: chance crossings against PC5's box

Part C counts, on two three-particle ensembles of [problem G](../analysis/world_clocks.md#test-problems-the-moshinsky-harmonic-model-q-and-the-smoothed-contact-interaction-g) at $`t = 1`$
($`\nu = 5\times10^4`$), how many worlds of the partner's sign lie in PC5's
box (Chebyshev half-side $`0.1`$ in the four relative coordinates) around a
target, and how many lie in the crossing tube with the same momentum
tolerance $`0.1`$, for impact parameters $`\beta`$ and steps $`\Delta t`$.
The events are drawn with probability proportional to $`\lvert K_q(r_b)\rvert`$,
20000 of them; the mean $`\lvert\mathbf K\rvert`$ is $`0.72`$, so the tube's length
at $`\Delta t = 0.01`$ is about $`0.007`$. The ensembles are absorption only
(PC5, paired, $`\delta = 0.25`$, no sink) and K1s10 at cell $`0.7`$.

**Proposition CX4 (measured).** *Partners per target, absorption-only
ensemble ($`0.430`$ in the box):*

| $`\beta`$ | $`\Delta t = 0.01`$ | $`0.04`$ | $`0.16`$ | uniform density, $`\Delta t = 0.01`$ |
|---|---|---|---|---|
| 0.1 | 0.358 (0.83 of the box) | 0.385 | 0.482 | 0.82 of the box |
| 0.02 | 0.0307 (0.071) | 0.0355 | 0.0562 | 0.039 |
| 0.005 | 0.0089 (0.021) | 0.0099 | 0.0153 | 0.0038 |

*The K1s10 ensemble gives the same picture ($`0.578`$ in the box; at
$`\beta = 0.005`$, $`\Delta t = 0.01`$, $`0.0090`$ per target, $`0.016`$ of the
box against $`0.0037`$ for a uniform density).*

Three things follow.

- With $`\beta`$ equal to the box's half-side the tube *is* the box with a
  round cross-section in configuration space: its share tends to
  $`\pi/4 \approx 0.79`$ as $`\Delta t \to 0`$, and it is $`0.83`$ at
  $`\Delta t = 0.01`$.
- Making the crossing sharper ($`\beta \to 0`$) loses partners roughly in
  proportion to $`\beta`$, as CX3 says for $`n = 2`$. A step long enough for
  the tube's length to matter is far longer than the step the dynamics
  allows.
- At small $`\beta`$ there are more partners than a uniform density gives:
  $`5.4`$ times as many at $`\beta = 0.005`$ and $`\Delta t = 0.01`$. Of the
  $`0.0089`$ partners per target there, $`0.0045`$ sit at *exactly* the
  parent's configuration (to $`10^{-9}`$), which only a world born at the
  parent's own vertex can do. The excess is descent, not chance (§6).

## 6. Proposition CX5: worlds that share a vertex, selected by record

Every event is a vertex in configuration space-time: parent and daughters
share a configuration point there, and their momenta differ by exactly
$`\pm\mathbf K`$, a lattice vector along a bond. Call two worlds **vertex kin** if
they share a vertex in their past: a world, its parent and its sibling at its
birth, and its own daughters. Kin escape Lemma CX3, but not because they meet
again. They do not: they leave the vertex along the bond at relative velocity
$`\pm\mathbf K/m`$ and part (CX6 says when the flow brings them back). What
kinship supplies is a **record**: the exact phase-space offset between the two
worlds at the vertex.

**The record.** A world $`C`$ born at $`t_v`$ on side $`s = \pm1`$ of an event of
$`G`$ on bond $`b`$, channel $`q`$ ($`C`$ at $`\mathbf p_G + s\mathbf K`$) has, at the
vertex,

```math
z_G - z_C = (0, -s\mathbf K), \qquad z_S - z_C = (0, -2s\mathbf K),
```

for its parent $`G`$ and sibling $`S`$, and $`z_D - z_C = (0, s'\mathbf K')`$ for each
daughter $`D`$ it emits later on side $`s'`$ of an event with transfer
$`\mathbf K'`$. Each offset is a bond and a signed multiple of $`\Delta p`$, and
the demo stores exactly that.

**Selection by record.** An event of $`A`$ on $`(b, q)`$ with sign $`\tau`$ needs a
partner of sign $`-\tau`$ at $`z_A + (0, \mathbf K)`$ and one of sign $`+\tau`$
at $`z_A - (0, \mathbf K)`$. A kin $`B`$ of $`A`$ is a partner *by record* for the
$`+`$ daughter if its recorded offset from $`A`$ is exactly $`(0, +\mathbf K)`$
(same bond, same channel) and its sign is $`-\tau`$, and likewise for the
$`-`$ daughter. What absorbing it does follows from Lemma CX2. In the
emission-only ensemble $`A`$'s own momentum does not kink (WC6), so with
$`a = t - t_v`$ the age of the vertex,

```math
z_B(t) = \Phi_a\bigl(z_A(t_v) + (0, \mathbf K)\bigr)
       = z_A(t) + \bigl(a\mathbf K/m + O(a^3),\ \mathbf K + O(a^2)\bigr).
```

Absorbing $`B`$ changes $`W`$ exactly as the event would have done had it
fired at $`t_v`$ and its daughter flowed since. The only error is the delay:
a translation by $`a\lvert\mathbf K\rvert/m`$ along the bond, and a momentum
drift of order $`a^2`$. The alternative, selecting by *current* momentum
(within $`0.05`$ of the target), tests the drifted momentum and so loses old
kin, and of the $`88158`$ partners it finds below only two are not kin by
record.

**Proposition CX5.** *(a) A kin absorption of both daughters is exactly an
undo: in an ensemble from which no world has been removed, both daughters of
an event of $`A`$ on $`(b, q)`$ with sign $`\tau`$ have partners by record if
and only if $`A`$ fired $`(b, q)`$ earlier with sign $`-\tau`$. (b) (Measured,
below.) A partner by record for one daughter is common; for both, it is under
one event in a hundred.*

*Proof of (a).* If $`A`$ fired $`(b, q)`$ with sign $`-\tau`$, its daughters have
sign $`-\tau`$ at offset $`+\mathbf K`$ and $`+\tau`$ at $`-\mathbf K`$: both
partners. Conversely, suppose both partners exist. Offsets on different bonds
point in different directions, so every partner's record is on bond $`b`$. If
$`A`$ was born on side $`s`$, its parent lies at $`-s\mathbf K'`$ and its sibling at
$`-2s\mathbf K'`$, both on side $`-s`$, so together they can serve at most the
half on that side. The other half, at least, is served by a daughter of $`A`$
at offset $`\pm\mathbf K`$: it was born on that side of an event of $`A`$ on
$`(b, q)`$, and has the needed sign $`\mp\tau`$ only if that event had sign
$`-\tau`$. That event's other daughter, still present because nothing has
been removed, is a partner for the remaining half. $`\square`$

*Measured* (Part D: [problem G](../analysis/world_clocks.md#test-problems-the-moshinsky-harmonic-model-q-and-the-smoothed-contact-interaction-g), emission only, no sink, $`\nu = 2\times10^4`$,
every world's record kept, per daughter unless stated):

| | to $`t = 0.25`$ | to $`t = 0.5`$ |
|---|---|---|
| events | 59448 | 481019 |
| a partner by record | 0.0687 | 0.0959 |
| &nbsp;&nbsp;of which parent / sibling / daughter of $`A`$ | 7310 / 480 / 374 | 78529 / 4933 / 8821 |
| a partner by current momentum (within 0.05) | 0.0684 | 0.0916 |
| record partners that pass the momentum test | 0.996 | 0.955 |
| both daughters by record (per event) | 0.0032 | 0.0093 |
| an earlier opposite event to undo (per event) | 0.0032 | 0.0093 |

The last two rows agree to the event ($`188`$ and $`4453`$), as (a) says. The
record partners are a median $`a = 0.080`$ old ($`90`$th percentile $`0.24`$),
and their momentum mismatch with the target (Chebyshev, relative
coordinates) grows with age as CX2 predicts, roughly as $`a^2`$:

| age $`a`$ | partners | median momentum mismatch with the target | fail the 0.05 test |
|---|---|---|---|
| $`(0, 0.05]`$ | 32602 | 0.0003 | 0.000 |
| $`(0.05, 0.15]`$ | 36563 | 0.0035 | 0.004 |
| $`(0.15, 0.3]`$ | 19097 | 0.0176 | 0.098 |
| $`(0.3, 0.6]`$ | 4021 | 0.0519 | 0.529 |

The undoable events are a median of 14 steps old ($`a = 0.14`$); 34 per cent
are at most 10 steps old and 0.9 per cent at most one.

**What S6 leaves of this.** A partner for *one* daughter, found for 18 per
cent of events by $`t = 0.5`$ ($`2 \times 0.0959 - 0.0093`$), is of no use. Absorbing one daughter and
emitting the other is the mixed realisation that Theorem S6 excludes: it
moves an unpaired world from $`\mathbf p - \mathbf K`$ to $`\mathbf p + \mathbf K`$, a
hop of $`2\xi_q`$ that nothing supplies (step 16; tutorial,
[`emission_and_absorption.md`](emission_and_absorption.md) §10; the
specification excludes Right-Hop and Left-Hop). What is admissible is the
paired case, which by (a) is re-absorption: an event that undoes an earlier
one of the same parent, bond and channel, $`0.3`$ to $`0.9`$ per cent of
events, typically ten steps late.

## 7. Proposition CX6: when the flow brings kin back

Kin part at their vertex. Whether they meet again is a property of the flow.

**Proposition CX6.** *Let red and blue be the daughters of an event of $`A`$
with transfer $`\mathbf K`$, at $`z_A \pm (0, \mathbf K)`$ at the vertex. (a) Under
free motion their configurations separate as $`2a\mathbf K/m`$ and never meet
again unless later events close a loop of descent. (b) If the flow is linear
(a quadratic Hamiltonian), the parent stays exactly at the midpoint of red and
blue for all time. (c) If moreover $`\mathbf K`$ excites only relative modes of
one frequency $`\omega`$, red and blue meet again exactly at the parent at
$`a = \pi/\omega`$, with momenta $`\mathbf p_A \mp \mathbf K`$: sides swapped, an
exact absorption vertex for an event of $`A`$ on the same bond and channel
with the sign of the first.*

*Proof.* (a) is CX2 with $`\mathbf F = 0`$. (b) A linear flow acts on offsets by
a matrix, $`\Phi_a(z + \delta) = \Phi_a(z) + M(a)\delta`$, so the offsets of red
and blue from the parent are $`\pm M(a)(0, \mathbf K)`$. (c) On an eigenspace of
frequency $`\omega`$, $`M(a)(0, \mathbf K) = (\mathbf K\sin\omega a/(m\omega),\ \mathbf K\cos\omega a)`$,
which is $`(0, -\mathbf K)`$ at $`a = \pi/\omega`$. Red has the event's sign
$`\tau`$ and is now at $`-\mathbf K`$; blue has $`-\tau`$ and is at $`+\mathbf K`$:
exactly the partners of an event of sign $`\tau`$. $`\square`$

Three particles of unit mass in a trap of frequency $`1`$: with no pair force
every relative mode has $`\omega = 1`$; with Moshinsky's harmonic pair force of
strength $`k = 2`$ every relative mode has $`\omega^2 = 1 + 3k = 7`$, and the
transfer $`\mathbf K`$, having no centre-of-mass component, excites only those.
[Problem G](../analysis/world_clocks.md#test-problems-the-moshinsky-harmonic-model-q-and-the-smoothed-contact-interaction-g)'s Gaussian pair force is anharmonic.

*Measured* (Part G: $`2000`$ kin triads from random parents, bonds and channels
$`\xi = 0.25`$ to $`1.75`$, flow integrated with step $`0.002`$ to $`t = 6`$):

| flow | returning by $`t = 6`$ | first return, median | closest / largest separation | parent from the midpoint | $`\lvert\mathbf p_{\rm red} - \mathbf p_{\rm blue}\rvert / 2\lvert\mathbf K\rvert`$ |
|---|---|---|---|---|---|
| trap only | 2000 | 3.142 ($`\pi`$) | $`4\times10^{-4}`$ | $`2\times10^{-15}`$ | 1.000 |
| Moshinsky, $`k = 2`$ | 2000 | 1.188 ($`\pi/\sqrt7 = 1.187`$) | $`1.6\times10^{-3}`$ | $`1\times10^{-15}`$ | 1.000 |
| problem G | 1970 | 3.44 (10–90 %: 2.27–4.23) | 0.215 (90 %: 0.61) | 1.11 | 0.985 |

(The closest approach in the harmonic cases is the sampling of the time grid,
not a miss.) In problem G the kin come back to about a fifth of their largest
separation, near half a trap period, and the parent is about $`1.1`$ from
their midpoint, ten times PC5's half-side. The refocusing is a property of
harmonic flows, not of problem G.

## 8. How crossings can be used

**In the ontology.** For $`N`$ particles an absorption is a vertex in
configuration space-time where three world-lines meet: the parent's and two
partners' converging along the bond from opposite sides, leaving as the
parent and an aligned pair (first figure, panel a). Q10 still holds: the
vertex realises an event the parent's reading has already triggered. What is
new is where such vertices come from when $`n \ge 2`$. They do not occur by
chance (CX3). Among kin the record says exactly where the vertex was (CX5),
and the paired case S6 admits is the undo of an earlier event of the same
family. In a harmonic flow kin are also brought back to an exact vertex by
the flow itself (CX6).

**As a realisation rule.** Two uses remain.

1. *Inside the step, a crossing is free.* Because event times are uncertain by
   $`\Delta t`$ already, any partner that crosses within the step realises the
   event exactly to the scheme's order (CX2). With momentum lattices
   (PC-SP4) and $`n = 1`$, chance crossings in the step are exact absorptions
   with positive probability (CX3(b)); that is the two-particle problem of
   step 25 §6 (CX-SP1).
2. *Re-absorption by record.* An event that undoes an earlier one of the same
   parent, bond and channel can be realised by absorbing that event's two
   daughters, with the delay as the only error (CX5). It is exact in the
   ledger's sense (the absorptive branch, $`N \to N - 2`$) and admissible under
   S6.

**What crossings cannot do.** They cannot replace a sink. Chance crossings
are too rare for $`n \ge 2`$ (CX4), paired kin are too rare (CX5), and the
flow refocuses kin exactly only when it is harmonic (CX6). Together they bring
the absorbed share nowhere near the one half that Theorem N3 asks of a ledger
without another sink. The population control for [problem G](../analysis/world_clocks.md#test-problems-the-moshinsky-harmonic-model-q-and-the-smoothed-contact-interaction-g) remains K1
resampling with the grid shifted every ten steps (Proposition PC6): it holds
$`P_{0.25}`$ to $`-0.020`$ of the mesh at $`\nu = 10^6`$, against $`-0.143`$ for
the contact sink at cell $`0.35`$, at similar population.

## 9. Figures

![Crossings of worlds in configuration space-time](https://raw.githubusercontent.com/billpage/wpmw/output/figures/world_crossings.png)

*(a) An absorption vertex, in the bond coordinate relative to the parent: the
negaton at momentum $`p + K`$ and the positon at $`p - K`$ approach along the
bond at relative velocities $`\pm K/m`$ and meet the parent; an aligned pair
at momentum $`p`$ leaves with it. (b) The crossing tube for $`n = 2`$: a
segment of length $`\lvert K\rvert\Delta t/m`$ along the bond, thickened by the
impact parameter; partners (arrows, their motion in one step) that pass
through it nearly cross the event. (c) Partners in the tube as a share of
those in PC5's box (half-side 0.1), against the step, for three impact
parameters, on the absorption-only ensemble (filled) and the K1s10 ensemble
(open), with the uniform-density prediction (grey).*

![Kin: emission, parting and return](https://raw.githubusercontent.com/billpage/wpmw/output/figures/world_crossings_kin.png)

*Kin on a line, drawn as lanes of an expressway: momentum is the lane, and
the axes are $`x`$ and $`t`$. (a) Free motion: at the vertex E the parent emits
red, in lane $`p + \xi`$, and blue, in lane $`p - \xi`$; they part for ever. Absorbing them at a later
time undoes E late, with error their separation (CX2). (b) A closed loop of
descent: at M₁ and M₂ each kin fires an event, and the line continues as its
daughter at $`\mp 2\xi`$; the two daughters reach the parent together at A in
lanes $`p \mp \xi`$. No world changes lane by itself (S6). (c) A harmonic trap:
every orbit has the same period, so the flow refocuses red and blue at the
parent after half a period, sides swapped (CX6). (d) Measured, three
particles: kin separation in the relative configuration over $`\lvert\mathbf K\rvert`$
against age, median and 10–90 per cent band of $`2000`$ triads, for the trap
alone, Moshinsky's harmonic pair force and problem G.*

## 10. Numerical verification

`PYTHONPATH=src python3 src/demo_world_crossings.py` runs Parts A to D, F and
G and both figures in about two and a half minutes (container, October 2026, at commit
`658de83` with this revision added). It imports the world algorithm and
kernel of `src/demo_three_worlds.py` and the ensembles of
`src/demo_population_control.py`. Seeds are fixed in the code.

- **Part A.** Lemma CX2, table of §3.
- **Part B.** Lemma CX3, the tube probabilities and slopes of §4.
- **Part C.** Proposition CX4, the partner counts of §5, including the
  partners at exactly the parent's configuration.
- **Part D.** Proposition CX5, partners by record and by momentum, the kin
  types, the identity of paired kin with undo, and the age tables of §6.
- **Part F.** The first figure.
- **Part G.** Proposition CX6, the return table of §7, and the second figure.

## 11. Open items

- **CX-SP1.** Two particles on a line, with the relative momentum on a
  lattice: realise absorptions at exact crossings inside the step, and
  compare with emission plus K1 against the mesh of step 25 §6.
- **CX-SP2.** *Withdrawn.* It proposed a mixed realisation (one daughter
  absorbed from kin, one emitted). Theorem S6 of step 16 already excludes it:
  it is a hop of $`2\xi_q`$ for an unpaired world (§6).
- **CX-SP3.** Kin through a sink: K1 replaces removed worlds by new pairs with
  no family, so records must be carried through the sink, or the pairs
  assigned a parentage, before CX5 can be measured alongside K1. Low priority:
  by CX5 the admissible yield is under one per cent of events.
- **CX-SP4.** In $`d = 3`$ the transfer $`\mathbf k`$ is a vector on the
  channel lattice of the bond, and the tube lies along $`\hat{\mathbf k}`$ in the
  bond's three-dimensional subspace. Check CX3 and CX5 for three particles in
  space ($`n = 6`$).

## 12. Sources

- Step 23, Propositions Q9 and Q10:
  [`../analysis/force_blind_sea.md`](../analysis/force_blind_sea.md) §9.
- Step 25, Propositions WC6 and WC9:
  [`../analysis/world_clocks.md`](../analysis/world_clocks.md) §1.1 and §8.2.
- Step 26, Proposition PC5 and open item PC-SP4:
  [`../analysis/population_control.md`](../analysis/population_control.md).
- Step 16, Theorem S6:
  [`../analysis/sea_population_equilibrium.md`](../analysis/sea_population_equilibrium.md) §6;
  tutorial, [`emission_and_absorption.md`](emission_and_absorption.md) §10.
- Step 26, Proposition PC6:
  [`../analysis/population_control.md`](../analysis/population_control.md) §5.
- Moshinsky, M. — "How good is the Hartree-Fock approximation." *Am. J.
  Phys.* **36**, 52 (1968). The harmonic three-body model of Part G.
  https://doi.org/10.1119/1.1974410
- Step 18, Theorem N3:
  [`../analysis/stochastic_ledger.md`](../analysis/stochastic_ledger.md) §4.
