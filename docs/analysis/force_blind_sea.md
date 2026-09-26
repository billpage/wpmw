# The force-blind sea: aligned pairs move inertially

> Postulate (S) gives every world-particle the full classical force, the members of an aligned sea pair included, and nothing in the QLE asks for that: a pair is `W`-null, and the Moyal equation fixes `u+ - u-` and says nothing about `u+ + u-`. This note provisionally adopts **(S′)**: free bodies obey the full force, aligned pairs move inertially — advecting at `p/m` with no drift in `p` — while their clocks still wind with `V`, so the sea is force-blind but potential-sensitive; a motionless pair would need a Hamiltonian clock, and a potential-blind one would leave the kernel nothing to be read from. **Proposition Q1** organises everything: in every ledger of the chain the sea enters an event only as its source or its sink and is never read, so the body fields — and with them `E`, `N` and `f` — are independent of how it moves, verified bitwise; the exceptions are exactly where the sea is read, a throttled rate, a sea-weighted kernel, its phases and identity. **Proposition Q2** gives the kinematics: under (S′) a sea clock stays on its row's plane wave up to the eikonal phase and a same-row pair winds at exactly `U`, so Theorem L8's drift term vanishes; under (S) `d(hbar theta - p x) = -H dt - x dp`, a uniform force keeps the row lock, and the mismatch begins at `V''`, which is the row mixing step 22 measured. **Proposition Q3**: what changes is the sea itself — deficits stay in their rows, Theorem S8's early dip goes but the worst cell hovers near zero, and fast recombination repairs it where slow does not. **Proposition Q4**: with step 20's tags now conserved, tagged bodies cross the Eckart barrier inside dark pairs, through the summit — `T_tag` 0.89 against 0.34 at `E0 = V0/2`, with only 9 per cent of the tag above the barrier energy — while `T_E` is unchanged, so individual crossing is still not tunnelling. Two demo defects are repaired: a species mask read at the destination row, the whole of step 22's drift of `Sum E` (L-SP8), and the tag bookkeeping of step 20's Part E. Section 8 lists the corrections to steps 15, 16, 17, 20 and 22 that (S′) would require.

*Ladder abstract — see the [full list](README.md#the-ladder).*

**Status.** Analysis note, step 23 of the ladder. **Provisional:** it adopts
postulate (S′) for aligned pairs in place of (S), to test the consequences
before the notes of steps 15 to 22 are revised. Companion demo:
`src/demo_force_blind_sea.py` (§§3–5). Theorem Y6's rerun uses
`src/demo_dark_sea_and_identity.py --parts E --full --sea S|blind` (§6), and
Proposition L2(a)'s uses `src/demo_contact_kernel.py 30 --sea S|blind` (§7).
Every demo keeps (S) as its default, so every published output is unchanged.

---

## 0. What this note asks, settles and leaves open

Postulate (S) of [`compensated_ontology.md`](compensated_ontology.md) gives
*every* world-particle the full classical force, the members of an aligned
sea pair included. Nothing in the QLE asks for that: an aligned pair puts
$`+1`$ and $`-1`$ in the same cell, so it contributes nothing to
$`E = u_+ - u_-`$, and the Moyal equation fixes $`u_+ - u_-`$ while saying
nothing about $`u_+ + u_-`$ ([`sea_population_equilibrium.md`](sea_population_equilibrium.md)).
The demos justified streaming the sea classically by saying that classical
streaming is the unique motion that carries a pair without separating it.
That is not so: any common motion keeps the members together, inertial
motion included. This note asks what changes if aligned pairs are blind to
the force.

**Settles.**

- (S′) is consistent with P1 and the plane-wave lock only if pairs keep
  advecting at $`p/m`$ and keep winding with $`V`$: force-blind, not
  motionless, and not blind to the potential (§1, Proposition Q2).
- **Proposition Q1.** In every ledger of the chain the sea enters an event
  only as its source or its sink, so the body fields — and with them $`E`$,
  the body count and the absorptive fraction $`f`$ — are exactly independent
  of how the sea moves. Verified bitwise. The exceptions are a rate throttled
  by the sea (Theorem S5) and a kernel weighted by it (Proposition L2(a)).
- **Proposition Q2.** Under (S′) a sea clock stays on its row's plane wave up
  to the eikonal phase $`-\int V\,dt/\hbar`$, a same-row pair winds at exactly
  $`U`$, and Theorem L8's drift term vanishes at all times. Under (S) a body's
  offset from its row's plane wave begins at $`V''`$, which is the row mixing
  step 22 measured.
- **Proposition Q3.** What (S′) changes is the sea itself: deficits stay in
  their rows and only translate (§5).
- **Proposition Q4.** Under (S′) a tagged body crosses the Eckart barrier
  inside a dark pair, through the summit. Individual crossing is still not
  tunnelling (§6).
- Two demo defects, repaired (§7): the destination-row species mask behind
  step 22's open item L-SP8, and non-conservation of tags in step 20's
  Part E.

**Leaves open** the corrections to steps 15 to 22 listed in §8, to be made
once (S′) is confirmed, and the items of §9.

**Inherits.** Postulate (S), (D) and Theorem G1 from step 17; Proposition K8
from step 15; Theorems S5, S7, S8 and S9 from step 16; Theorems Y5 and Y6
from step 20; Theorems L4 and L8 and Proposition L9 from step 22; P1 and
Lemma 0 from step 4.

---

## 1. Postulate (S′)

**Postulate (S′) — streaming with a force-blind sea.** Between events every
*free* body
obeys

```math
\dot q = p/m, \qquad \dot p = -\nabla V(q)
```

as in (S). Each member of an aligned pair instead moves inertially,

```math
\dot q = p/m, \qquad \dot p = 0,
```

while its clock still winds by P1, $`\hbar\dot\theta = p^2/2m - V(q)`$.

Nothing else changes. Events are as in (D): a parent ionises a pair on its
own row (Proposition K8) or aligns two bodies into one there, and the pair
keeps that row until it is ionised. Momentum is continuous along every
worldline between events, as before; what changes is that the stretch a
body spends in an aligned pair is a straight line in phase space.

Two parts of (S′) are fixed rather than chosen.

- **A pair must advect at $`p/m`$.** With P1's Lagrangian clock, a clock
  moving at $`p/m`$ keeps the phase of the plane wave
  $`e^{i(px - p^2t/2m)/\hbar}`$ in free space (Lemma 0 of step 4). A
  motionless pair would keep it only with a Hamiltonian clock,
  $`\hbar\dot\theta = -p^2/2m`$ (Proposition Q2).
- **A pair must wind with $`V`$.** Its clock rate carries the potential at
  its own position, and step 22's readings take the kernel's signal from
  exactly that (Theorems L3, L8). A sea blind to the potential as well would
  leave the kernel nothing to be read from.

So the sea of (S′) is **force-blind but potential-sensitive**. The scalar
Aharonov–Bohm effect (Aharonov and Bohm 1959) is an analogy, not a model:
there a charge acquires phase from a potential in a region where it feels no
force. Here there is a gradient, and the pair ignores it.

---

## 2. Why darkness points to (S′)

Four arguments, the last of them numerical.

*The force has nothing to act on.* In the compensated split the classical
force is the particle realisation of the Liouville term
$`-V'(x)\,\partial_p W`$, an operator on the signed density. An aligned pair
has no share of $`W`$. The QLE assigns it no trajectory at all, and (S) for
pairs was a choice, justified by an argument (co-location) that does not
single it out.

*Darkness as a limit.* A pair is two opposite-signed amplitudes at one point,
cancelling exactly: dark in every signed sum (Theorem L4). Step 20 separated
two meanings of that darkness — *observationally dark* (contributes nothing
to $`E`$, arithmetic) and *dynamically inert* (ineligible as a recombination
donor, a rule). (S′) extends the second meaning from recombination to the
force. It does not replace the exclusion rule, which rests on K8.

*The force enters the reading once, at the reader.* Theorem L8 showed that
the compensated split is the choice of the parent as reader: the force
removed from the kernel is the reader's acceleration. Nothing in the reading
needs the partners to accelerate, and under (S′) they do not.

*The lock holds.* Step 22 §9 found that a sea keeping its rows, re-locked at
the kernel's own rate, reads the kernel at 0.82 against a 0.89 control,
while under (S) the reading was 0.45. The reason is Proposition Q2: under
(S) the force fans a row out over its neighbours, and each row gathers
bodies carrying phases locked to the rows they came from.

**Costs.** A pair's energy $`p^2/2m + V(q)`$ is not conserved while it
crosses a potential. It is unobservable, because the pair is $`W`$-null, but
it is a real departure from Newtonian bookkeeping for paired bodies. And
identity now passes through classically forbidden regions in the dark
(§6): a claim about an unobservable, which leaves $`E`$ alone.

---

## 3. Proposition Q1: the observable is blind to the sea's motion

Write one ledger step as transport, events, transport,

```math
(u_+, u_-, s) \;\xrightarrow{\;T_{\delta t/2}\;}\;
\xrightarrow{\;C_{\delta t}\;}\;
\xrightarrow{\;T_{\delta t/2}\;}\; (u_+', u_-', s'),
```

where $`u_\pm`$ are the positon and negaton densities on the $`(x, p)`$
mesh and $`s`$ the density of aligned pairs. Transport acts field by field:
the bodies by the flow $`\Phi`$ of (S), the sea by a flow $`\Psi`$, which is
$`\Phi`$ under (S) and the free flow $`\Phi_0`$ under (S′).

**Proposition Q1.** *If the event update has the triangular form*

```math
u_\pm \mapsto u_\pm + a_\pm(u_+, u_-), \qquad s \mapsto s + c(u_+, u_-),
```

*with the increments computed from the body fields and the kernel alone,
then the body trajectory $`(u_+, u_-)(t)`$ — and every functional of it,
$`E`$, $`N = u_+ + u_-`$ and the absorptive fraction $`f`$ among them — is
the same for every sea flow $`\Psi`$.*

*Proof.* The transports of $`u_\pm`$ do not involve $`s`$, and by hypothesis
neither do the increments $`a_\pm`$. So each step maps $`(u_+, u_-)`$ to a
value independent of $`s`$ and of $`\Psi`$; by induction on the steps, so
does the whole trajectory. $`\square`$

The ledgers of the chain have this form. In each event channel the demand is
$`D = \lvert K_q\rvert\,u_{\rm parent}\,\delta t`$; the absorptive count
$`A = \min(D, \text{caps})`$ takes its caps from the opposite-species body
fields; the emissive count is $`D - A`$; the bodies change by $`\pm A`$ and
$`\pm(D - A)`$ at the daughter rows; and the sea changes by $`A - (D - A)`$
at the parent's row. The sea is a source for emission and a sink for
absorption and is never read. The population clamp acts on $`u_\pm`$ only.

**Exceptions.** The hypothesis fails exactly where the sea is read:

- a rate throttled by the sea, $`\Gamma \to \Gamma\,s/B`$ (Theorem S5);
- a kernel weighted by the sea, $`K_q(x) \to \sum_j C\,s(x + y_j, p)/B`$
  (Proposition L2(a));
- any rule that reads the sea's phases (the readings of step 22);
- anything that follows *identity*: tags riding in pairs move with the sea
  (§6).

*Verification.* `demo_force_blind_sea.py` Part A:

| ledger | run | body fields under (S) and (S′) | sea, largest difference |
|---|---|---|---|
| shared ledger (steps 20, 22), clamped | Eckart, $`p_0 = 1.2`$, 300 steps | $`u_+`$, $`u_-`$ bitwise equal | 0.267 $`B`$ |
| shared ledger, unclamped | the same | $`u_+`$, $`u_-`$ bitwise equal | 0.250 $`B`$ |
| step 16 ledger, absorptive | $`T = 6`$, $`\delta t = 0.01`$ | $`E`$ bitwise equal; $`N`$ = 2.9674104244 and $`f`$ = 0.3997855351 under both | worst cell $`-0.339\,B`$ and $`-0.060\,B`$ |

So every result of the chain that is a statement about $`E`$, $`N`$ or $`f`$
in a ledger without the exceptions stands under (S′) unchanged, to the bit.

---

## 4. Proposition Q2: the kinematics of (S) and (S′)

Let $`\varphi_p(x, t) = (px - p^2t/2m)/\hbar`$ be the phase of row $`p`$'s
plane wave.

**Proposition Q2.** *(a) Under (S′) a sea clock obeys*

```math
\frac{d}{dt}\Big[\hbar\theta - \hbar\varphi_p(q, t)\Big] = -V(q),
```

*so it stays on its row's plane wave up to the eikonal phase
$`-\hbar^{-1}\int V(q(t'))\,dt'`$ accumulated along its straight path. A
motionless pair would keep the plane wave only with
$`\hbar\dot\theta = -p^2/2m`$, a Hamiltonian clock, not P1's.*

*(b) Two partners on one row keep their separation $`2y`$, and a reader with
momentum $`P(t)`$ sees them wind at*

```math
\hbar\dot\mu = V(x_j) - V(x_i) + 2y\,\dot P
```

*exactly, at all times: an inertial reader reads $`U`$, the parent reads
$`U_{\rm res} = U - 2yV'_{\rm eff}`$, and Theorem L8's drift term vanishes
identically, since $`p_i = p_j`$ for ever.*

*(c) Under (S) a body's phase obeys the exact identity*

```math
d(\hbar\theta - px) = -H\,dt - x\,dp ,
```

*so its offset from its current row's plane wave records where along its
path it received its impulses. A uniform force leaves two bodies of one row
in one row with zero misalignment; in general the offset rate is
$`-V(x_a) + \tfrac12X^2V'' - \tfrac16X^3V''' + \dots`$, with $`X`$ the
distance from the point $`x_a`$ to which the lock's anchor has streamed. The
first term is common to the row; the rest is tidal and begins at $`V''`$.*

*Proof.* (a) With $`\dot q = p/m`$, $`\dot p = 0`$ and P1,
$`\hbar\dot\theta - \hbar\,d\varphi_p/dt = p^2/2m - V - (p^2/m - p^2/2m) = -V`$.
(b) Differentiate $`\mu = \theta_i - \theta_j + P(x_j - x_i)/\hbar`$ with
equal, constant momenta; the kinetic terms cancel and $`\dot x_j = \dot x_i`$.
(c) $`\hbar\dot\theta - d(px)/dt = L - p\dot x - x\dot p = -H - x\dot p`$;
the uniform-force and Taylor statements follow by direct integration and
expansion. $`\square`$ All three are checked symbolically in Part B of the
demo.

Part (b) is the reason step 22's row-keeping sea held the reading: under
(S′) the own-frame reading of Theorem L7 becomes the inertial reader and
reads $`U`$, the parent frame reads $`U_{\rm res}`$ with no drift, and
Proposition L9's bound $`dp^2/m`$ becomes zero. Part (c) is the reason (S)
did not: a curved potential fans each row out over its neighbours.

---

## 5. Proposition Q3: what (S′) does to the sea

Under (S′) a row of sea translates rigidly at its own velocity, so a deficit
or excess made at a parent's row stays in that row. Under (S) the force
shears it across rows near the barrier. The step 16 ledger under both
motions (`demo_force_blind_sea.py` Part C):

Emissive unravelling with recombination rate $`\kappa`$ (Theorems S4 and
S5's setting, step 16 Part C), worst cell $`s/B`$ at $`T = 6`$:

| $`\kappa`$ | (S) | (S′) |
|---|---|---|
| 0 | −9918.5 | −9571.2 |
| 20 | +0.006 | −0.673 |
| 200 | −0.145 | +0.005 |
| 2000 | −0.141 | +0.073 |

Absorptive unravelling (Theorem S8), $`\delta t = 0.01`$: the worst cell
over $`0 \le t \le 6`$ is $`-0.339\,B`$ under (S) and $`-0.060\,B`$ under
(S′); at later times,

| $`t`$ | 2 | 6 | 10 | 14 | 18 |
|---|---|---|---|---|---|
| (S) | −0.175 | +0.191 | +0.165 | +0.195 | +0.163 |
| (S′) | +0.058 | −0.026 | +0.057 | −0.027 | −0.066 |

A rate throttled by the sea (Theorem S5), relative $`L^2`$ error of $`E`$
against the QLE at $`T = 8`$:

| $`\kappa`$ | (S) | (S′) |
|---|---|---|
| 200 | 0.251 | 0.285 |
| 2000 | 0.244 | 0.301 |

**Proposition Q3 (measured).** *Under (S′) the sea's excursions are shallower
at first and more persistent afterwards. The deep early dip of Theorem S8's
absorptive run is gone, but the worst cell hovers near zero instead of
recovering to about $`0.17\,B`$, because a deficit that is not sheared
across rows is not refilled. Fast recombination ($`\kappa \ge 200`$) keeps
the worst cell non-negative under (S′), where under (S) it does not; slow
recombination ($`\kappa = 20`$) does the reverse. A rate throttled by the
sea (Theorem S5) is observable under either motion, somewhat more under
(S′).*

The published S5 table (0.402, 0.414) predates the transport repair recorded
in the erratum of step 16; the current code gives the values above for (S).

**Theorem S9 under (S′)** (Part D, `--heavy`). The traces of the absorptive
fraction $`f`$ and the body count $`N`$ over $`0 \le t \le 8`$ are identical
element by element under the two motions, from the minimal ensemble,
$`\rho = 1`$, and from a twentyfold padded one, $`\rho = 20`$: the
attractor is a statement about bodies, as Proposition Q1 requires. The sea's
final worst cell differs, $`+0.168\,B`$ against $`+0.002\,B`$ at
$`\rho = 1`$ and $`-5.60\,B`$ against $`-1.75\,B`$ at $`\rho = 20`$.

---

## 6. Proposition Q4: identity crosses in the dark

Theorem Y6 followed tagged positons through events by proportional
allocation and streamed tags riding in pairs with the sea. Under (S′) a tag
absorbed into a pair on the approach slope coasts, and a pair ignores the
barrier.

**A defect first.** The published Part E did not conserve tags: the total
grew to 2.25 times the initial tag. An ionisation beyond the local sea (the
unclamped sea goes negative) released more tag than the pairs held, and
clipping at the end of each step turned negative tag into positive. The
repaired `channels_tagged` moves tags only, bounds every take by the tag
present, and applies each species mask at the parent's row. That leaves one
further leak, in the transport: spectral streaming is linear and conservative
but not positivity-preserving, and the channels act only on the positive part
of a tag field, so an uncorrected ripple is pumped into the free tags as
negative mass. The demo removes the ripple after each transport and rescales
the positive part to restore the total, and reports the mass so relocated.
The fields — and so $`E`$, $`f`$ and $`\lvert E - \text{mesh}\rvert`$ — are
unchanged to every printed digit.

**Proposition Q4 (measured).** *With tags conserved exactly, and the fine
lattice of Y6:*

| $`p_0`$ | $`E_0/V_0`$ | $`T_{cl}`$ | $`T`$ exact | $`T_E`$ ledger, both | $`T_{tag}`$ (S) | $`T_{tag}`$ (S′) |
|---|---|---|---|---|---|---|
| 1.0 | 0.50 | 0.069 | 0.191 | 0.2126 | 0.338 | **0.892** |
| 1.2 | 0.72 | 0.246 | 0.369 | 0.4038 | 0.422 | **0.914** |
| 1.7 | 1.44 | 0.896 | 0.845 | 0.8549 | 0.770 | **0.940** |

*and, per initial tag at the end of the run,*

| $`p_0`$ | tag with $`H > V_0`$, (S) / (S′) | $`\langle H\rangle/E_0`$ | kinks | transmitted tag in pairs | ripple relocated |
|---|---|---|---|---|---|
| 1.0 | 0.214 / 0.089 | 2.84 / 1.63 | 1.10 / 1.12 | 0.98 / 0.99 | 0.88 / 0.32 |
| 1.2 | 0.278 / 0.189 | 1.98 / 1.36 | 1.09 / 1.09 | 0.98 / 0.99 | 0.95 / 0.28 |
| 1.7 | 0.717 / 0.800 | 1.40 / 1.09 | 1.03 / 1.03 | 0.96 / 0.97 | 0.63 / 0.24 |

*Under (S′), at $`E_0 = V_0/2`$, 89 per cent of the tag ends beyond the
summit while only 9 per cent has the energy to pass over it: the rest
crossed inside pairs, through the summit, where $`\Gamma(0) = 0`$ and no
event happens. A body is taken into a pair on the approach and coasts
through dark; most of the transmitted tag is still in pairs at the end of
the run. $`T_E`$ is unchanged, as Proposition Q1 requires.*

Under both motions almost all of the transmitted tag sits in pairs at the
end, so that column does not distinguish them; what does is $`T_{tag}`$
together with the energy column. Under (S) a pair with $`H < V_0`$ is turned
back by the barrier like any body; under (S′) it is not. Kink counts are
about one per tag under either motion: the published Y6 counts (2.40, 2.32,
1.83) were inflated by the clipping.

Y6's deeper conclusion survives and is strengthened: individual crossing
tracks neither the classical nor the quantum transmission, so it is not
tunnelling. Its mechanism does not survive: "an individual crosses only by
flank activation followed by classical passage over the top" holds under
(S), not under (S′).

The ripple column matters for how much weight the numbers bear. As a check
on the coarse lattice ($`p_0 = 1.2`$), the published clip-based treatment
and the conservative one give $`T_{\rm tag}`$ = 0.311 and 0.322 under (S):
two very different handlings of the ripple agree to 0.01, while (S′) moves
$`T_{\rm tag}`$ by about 0.45.

---

## 7. Two demo defects

**The species mask of `channels_k` (step 22 open item L-SP8).** The sea-weighted
ledger of Proposition L2(a) chose the species of each deposit with a mask
read at the destination row $`p \pm q`$ rather than the parent's row $`p`$.
With a kernel whose sign depends on $`x`$ only the two agree, which is why
`Ledger.channels` is unaffected; the sea-weighted kernel's sign varies with
$`p`$, and each flip misassigned a deposit. The drift of $`\sum E`$ to
1.0010 was that, not physics. With the mask at the parent's row
(`demo_contact_kernel.py 30 --sea S|blind`, part C):

| sea | species mask | $`T_E`$ | $`\lvert E - \text{mesh}\rvert/\lvert\text{mesh}\rvert`$ | $`\sum E`$ | largest $`\lvert s/B - 1\rvert`$ | spurious first moment |
|---|---|---|---|---|---|---|
| (S) | destination row (published) | 0.4297 | 0.038 | 1.0010 | 0.417 | 0.135 |
| (S) | parent's row | 0.4298 | 0.038 | 1.0000 | 0.417 | 0.135 |
| (S′) | parent's row | 0.4332 | 0.039 | 1.0000 | 0.552 | 0.159 |

Proposition L2(a)'s conclusion — weighting contacts by the actual sea moves
$`E`$ only slightly — stands under either motion. Under (S′) the sea departs
further from $`B`$, because its deficits are not sheared across rows.

**The tag bookkeeping of step 20's Part E** (§6).

---

## 8. Corrections to make if (S′) is confirmed

- **Step 15** ([`eckart_barrier_compensated.md`](eckart_barrier_compensated.md)).
  §8's "every member of the created sea alike follows a genuine Newtonian
  worldline under the full classical force, for its entire life" becomes a
  statement about free bodies; stretches spent in a pair are straight. K-LS7's premise
  "pairs stream, so the realised profile is $`\Gamma`$ transported by the
  classical flow" becomes transport along rows.
- **Step 16** ([`sea_population_equilibrium.md`](sea_population_equilibrium.md)).
  The results about $`E`$, $`N`$ and $`f`$ stand (Proposition Q1). The sea's
  own profile changes (Proposition Q3), and S4's "fast recombination does not
  repair it" does not hold under (S′). The S5 table is stale under either
  motion. The demo docstring's "unique motion" argument is withdrawn.
- **Step 17** ([`compensated_ontology.md`](compensated_ontology.md)). (S) is
  restated as (S′); G1, G3 and G4 stand (Proposition Q1).
- **Step 20** ([`dark_sea_and_worldline_identity.md`](dark_sea_and_worldline_identity.md)).
  Y1 to Y5 stand. Y6's table is replaced by §6's; under (S) its
  transmissions survive the tag repair to about 0.01, but its kink counts
  were inflated about twofold by the clipping. Under (S′) its mechanism
  ("never through") is withdrawn, and the "dynamically inert" row of the
  darkness table gains the force.
- **Step 22** ([`sea_phase_reference.md`](sea_phase_reference.md)). L7's
  negative result is specific to (S) partners; under (S′) the own frame reads
  $`U`$ and L9's drift term is zero (Proposition Q2). L-SP8 is answered as a
  demo defect (§7), and Proposition L2(a)'s figures change slightly. L-SP10
  is answered provisionally by this note. §5's phase-field figure depicts (S)
  transport.

Steps 14, 18, 19 and 21 need nothing: their uses of (S) concern bodies.

---

## 9. Open items

- **Q-SP1.** Replace the spectral transport of tags by a positivity-preserving
  one (semi-Lagrangian or particle tags) and confirm §6 without the ripple
  correction.
- **Q-SP2.** Under (S′) the sea's worst cell hovers near zero (Proposition Q3).
  Does supply for emission become limiting in longer or deeper runs, where
  (S) would recover?
- **Q-SP3.** An aligned pair's energy is not conserved under (S′) while it crosses
  a potential. Is there any ledger quantity, observable or not, that the
  project has treated as conserved and that this breaks?
- **Q-SP4.** Carry out the corrections of §8 once (S′) is confirmed.

---

## 10. Sources

- Aharonov, Y. and Bohm, D. *Significance of electromagnetic potentials in
  the quantum theory*, Phys. Rev. **115** (1959) 485–491.
- Glauber, R. J. *High-energy collision theory*, in Lectures in Theoretical
  Physics, vol. 1, ed. W. E. Brittin and L. G. Dunham, Interscience (1959)
  315–414. The eikonal: a plane wave carried along straight lines,
  accumulating $`-\int V\,dt/\hbar`$.
- Moyal, J. E. *Quantum mechanics as a statistical theory*, Proc. Cambridge
  Philos. Soc. **45** (1949) 99–124.
