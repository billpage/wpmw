# Crossings of worlds in configuration space-time

**What it means for worlds to cross when there are three or more particles, why
chance crossings almost never happen, how worlds descended from one event meet
anyway, and what that offers the absorptive realisation of an event.**

---

## 0. Status and provenance

Supplement to step 26 ([`../analysis/population_control.md`](../analysis/population_control.md)),
prompted by Bill Page's question: can crossings of worlds within a cell and a
time step, rather than proximity in phase space, decide when an event is
realised by absorption, and how is a crossing defined when there are three or
more particles? Companion demo: `src/demo_world_crossings.py`, which produces
every number below and the figure.

It inherits the one-dimensional picture of step 23
([`../analysis/force_blind_sea.md`](../analysis/force_blind_sea.md) §9):
Proposition Q9(b) reads an absorption as two unpaired worlds on rows
$`p \pm \xi_q`$ crossing at the parent and binding at its row, and
Proposition Q10 shows that crossings can *realise* events but cannot
*trigger* them, because a crossing rate is mass action and the observable
would not close. Nothing here changes Q10: the parent's reading triggers the
event, and the question is only how the event is realised. It also uses
Proposition WC6 (a world's history: momenta kink only at events, positions
are continuous), Proposition WC9 (events act on one bond at a time) and, from
step 26, Proposition PC5 (absorption with a tolerance in phase space) and open
item PC-SP4 (exact absorption off a mesh).

New here: a definition of a crossing for $`N`$ particles (CX1), the
equivalence of a delayed realisation with a translated one (CX2), the
codimension of a chance crossing (CX3), the measured yield of crossings
against PC5's box (CX4), and the observation that worlds descended from one
event meet exactly, by construction (CX5). Labels CX1 to CX5 and open items
CX-SP1 to CX-SP4 belong to this supplement.

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

*Verification* (Part A, problem G of step 25, 4000 random parents and
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

Part C counts, on two three-particle ensembles of problem G at $`t = 1`$
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

## 6. Proposition CX5: worlds that share a vertex meet by construction

Every event is a vertex in configuration space-time: parent and daughters
share a configuration point there, and their momenta differ by exactly
$`\pm\mathbf K`$, a lattice vector along a bond. Call two worlds **vertex kin** if
they shared a configuration point at some event: a world, its parent and its
sibling at its birth, and its own daughters. Kin escape Lemma CX3. They do
not meet by chance, they *began* at the same point, and their momentum
differences are lattice vectors when they part. By Lemma CX2, a kin partner
that left the common vertex an age $`a`$ ago, with the right momentum
difference, is an exact crossing delayed by $`a`$: the error is a
translation by $`a\lvert\mathbf K\rvert/m`$ along the bond.

The cleanest case is **re-absorption**. A parent fires bond $`b`$, channel
$`q`$ with sign $`\tau`$, and earlier fired the same bond and channel with
sign $`-\tau`$, whose daughters, $`-\tau`$ at $`\mathbf p + \mathbf K`$ and
$`+\tau`$ at $`\mathbf p - \mathbf K`$, are still present. They are exactly the two
partners the new event needs: absorbing them undoes the earlier event, and
the only error is its age.

**Proposition CX5 (measured).** *Problem G, emission only and no sink,
$`\nu = 2\times10^4`$, every world's kin recorded (Part D):*

| | to $`t = 0.25`$ | to $`t = 0.5`$ |
|---|---|---|
| events | 59448 | 481019 |
| a kin partner for one daughter (momentum within 0.05) | 0.134 | 0.176 |
| kin partners for both daughters | 0.0029 | 0.0069 |
| an earlier opposite event to undo (re-absorption) | 0.0032 | 0.0093 |

*The undoable events are a median of 14 steps old ($`a = 0.14`$); 34 per cent
are at most 10 steps old and 0.9 per cent at most one.*

So descent supplies, for about one event in six, a partner for *one* of its
daughters that coincides with the target up to the delay since a shared
vertex. Partners for *both* daughters, which the ledger's absorptive branch
needs, are as rare as re-absorption: under one event in a hundred, and
typically ten steps late.

## 7. How crossings can be used

**In the ontology.** For $`N`$ particles an absorption is a vertex in
configuration space-time where three world-lines meet: the parent's and two
partners' converging along the bond from opposite sides, leaving as the
parent and an aligned pair (panel a). Q10 still holds: the vertex realises an
event the parent's reading has already triggered. What is new is where such
vertices come from when $`n \ge 2`$. They do not occur by chance (CX3); they
occur by descent, among worlds that already shared a vertex (CX5). In that
reading an absorption is a family relation: the worlds an event takes back
are, at best, the worlds an earlier event of the same family gave out.

**As a realisation rule.** Three uses, in order of reach:

1. *Inside the step, a crossing is free.* Because event times are uncertain by
   $`\Delta t`$ already, any partner that crosses within the step realises the
   event exactly to the scheme's order (CX2). With momentum lattices
   (PC-SP4) and $`n = 1`$, chance crossings in the step are exact absorptions
   with positive probability (CX3(b)); that is the two-particle problem of
   step 25 §6 (CX-SP1).
2. *Kin first.* Among the partners a rule accepts (PC5's box, or mode C's
   cell with a K1 credit), prefer kin, ranked by age: their error is the
   delay alone, a translation along the bond, and the K1 credit carries its
   first moment.
3. *One daughter at a time.* Kin partners for one daughter are common (CX5).
   Realising one daughter by absorption and the other by emission changes
   $`W`$ correctly and leaves $`N`$ unchanged, but it is neither of the
   ledger's two branches, which move $`N`$ by $`\pm2`$ (CX-SP2).

**What crossings cannot do.** They cannot replace a sink. Chance crossings
are too rare for $`n \ge 2`$ (CX4), and kin partners for both daughters are
too rare (CX5), to bring the absorbed share near the one half that Theorem N3
asks of a ledger without another sink.

## 8. Figure

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

## 9. Numerical verification

`PYTHONPATH=src python3 src/demo_world_crossings.py` runs Parts A to D and
the figure in about two minutes (container, October 2026, at commit
`4ef04cd` with this demo added). It imports the world algorithm and kernel of
`src/demo_three_worlds.py` and the ensembles of
`src/demo_population_control.py`. Seeds are fixed in the code.

- **Part A.** Lemma CX2, table of §3.
- **Part B.** Lemma CX3, the tube probabilities and slopes of §4.
- **Part C.** Proposition CX4, the partner counts of §5, including the
  partners at exactly the parent's configuration.
- **Part D.** Proposition CX5, the kin and re-absorption counts and the age
  distribution of §6.
- **Part F.** The figure.

## 10. Open items

- **CX-SP1.** Two particles on a line, with the relative momentum on a
  lattice: realise absorptions at exact crossings inside the step, and
  compare with emission plus K1 against the mesh of step 25 §6.
- **CX-SP2.** A mixed realisation (one daughter absorbed from kin, one
  emitted, $`N`$ unchanged): is it admissible in the ledger, which knows only
  the emissive and absorptive branches, and what does it do to the
  population and the bias?
- **CX-SP3.** Kin through a sink: K1 replaces removed worlds by new pairs with
  no family, so kin must be carried through the sink, or the pairs assigned a
  parentage, before CX5 can be measured alongside K1.
- **CX-SP4.** In $`d = 3`$ the transfer $`\mathbf k`$ is a vector on the
  channel lattice of the bond, and the tube lies along $`\hat{\mathbf k}`$ in the
  bond's three-dimensional subspace. Check CX3 and CX5 for three particles in
  space ($`n = 6`$).

## 11. Sources

- Step 23, Propositions Q9 and Q10:
  [`../analysis/force_blind_sea.md`](../analysis/force_blind_sea.md) §9.
- Step 25, Propositions WC6 and WC9:
  [`../analysis/world_clocks.md`](../analysis/world_clocks.md) §1.1 and §8.2.
- Step 26, Proposition PC5 and open item PC-SP4:
  [`../analysis/population_control.md`](../analysis/population_control.md).
- Step 18, Theorem N3:
  [`../analysis/stochastic_ledger.md`](../analysis/stochastic_ledger.md) §4.
