# The force-blind sea: aligned pairs move inertially

> Postulate (S) gives every world-particle the full classical force, the members of an aligned sea pair included, and nothing in the QLE asks for that: a pair is `W`-null, and the Moyal equation fixes `u+ - u-` and says nothing about `u+ + u-`. This note provisionally adopts **(S′)**: free bodies obey the full force, aligned pairs move inertially — advecting at `p/m` with no drift in `p` — while their clocks still wind with `V`, so the sea is force-blind but potential-sensitive; a motionless pair would need a Hamiltonian clock, and a potential-blind one would leave the kernel nothing to be read from. **Proposition Q1** organises everything: in every ledger of the chain the sea enters an event only as its source or its sink and is never read, so the body fields — and with them `E`, `N` and `f` — are independent of how it moves, verified bitwise; the exceptions are exactly where the sea is read, a throttled rate, a sea-weighted kernel, its phases and identity. **Proposition Q2** gives the kinematics: under (S′) a sea clock stays on its row's plane wave up to the eikonal phase and a same-row pair winds at exactly `U`, so Theorem L8's drift term vanishes; under (S) `d(hbar theta - p x) = -H dt - x dp`, a uniform force keeps the row lock, and the mismatch begins at `V''`, which is the row mixing step 22 measured. **Proposition Q3**: what changes is the sea itself — deficits stay in their rows, Theorem S8's early dip goes but the worst cell hovers near zero, and fast recombination repairs it where slow does not. **Proposition Q4**: with step 20's tags now conserved, tagged bodies cross the Eckart barrier inside dark pairs, through the summit — `T_tag` 0.89 against 0.34 at `E0 = V0/2`, with only 9 per cent of the tag above the barrier energy — while `T_E` is unchanged, so individual crossing is still not tunnelling. **Proposition Q5** says what the sea computes: the Fourier phase of channel `q` is the daughter row's P2 lever, `mu^(p) + xi_q d/hbar = mu^(p + xi_q)`, the jump rate is the rate of change of the row's interference sum read against the daughter, `K_q = d(Re z_q)/dt / (B dp dx)` at re-lock, and the sea's density is the one at which the channel basis is orthogonal on the reach — a Fraunhofer transform formed by superposition, at critical sampling. **Proposition Q6** says what dark catalysis is for: left alone, the force-blind sea relaxes exactly to each row's eikonal wave, and the reading needs the free one; dark catalysis is the reset clock that puts it back, the misalignment it leaves enters the reading as a quadrature term linear in the mean reset time, gauge invariance selects its relative resets over absolute ones, and near the Eckart barrier the kernel's own rate falls a little short (0.863 against a control of 0.880) while three times it nearly suffices; step 22's plateau was the free bodies read as partners, not incomplete resetting. **Propositions Q7 and Q8** turn to the sea's population. The `kappa` of the step 16 comparison is the contact sink, which is not in the model, and dark catalysis, with `dS = 0`, cannot touch the population at all; the model's own clocks are the event rate and the rate at which force-feeling bodies slip across force-blind rows, which are of the same order by construction (`Gamma ~ |V'|/dp`), so under (S′) the sea's deficits are shallow but never refilled. An immediate contact sink would be garbage collection — the signed-particle method's same-cell annihilation with a pair ledger added — keeping the ensemble minimal, making absorption redundant for the bodies, and fixing the sea's total at half the growth of `||W||_1`; a pair formed that way is generically gray, and whether it is force-blind on joining or only on alignment is left open. Three demo defects are repaired: a species mask read at the destination row, the whole of step 22's drift of `Sum E` (L-SP8); the tag bookkeeping of step 20's Part E; and the plane-wave lock on a periodic box. Section 10 lists the corrections to steps 15, 16, 17, 20 and 22 that (S′) would require.

*Ladder abstract — see the [full list](README.md#the-ladder).*

**Status.** Analysis note, step 23 of the ladder. **Provisional:** it adopts
postulate (S′) for aligned pairs in place of (S), to test the consequences
before the notes of steps 15 to 22 are revised. Companion demo:
`src/demo_force_blind_sea.py` (§§3–6). Theorem Y6's rerun uses
`src/demo_dark_sea_and_identity.py --parts E --full --sea S|blind` (§7),
Proposition L2(a)'s uses `src/demo_contact_kernel.py 30 --sea S|blind` (§9),
and §8's scan is `src/scan_dark_reset.py`.
Every demo keeps (S) as its default, so every published output is unchanged.

**Addendum (September 2026).** §1 now fixes what a partner is and says the
sea is a random configuration, not a lattice. Proposition Q2 names its
eikonal (the straight-line, Glauber one, against the optical one that (S)
clocks carry) and is checked against the particle model. §5, Proposition
Q5 (what the sea computes: the kernel as an interference), and §8,
Proposition Q6 (dark catalysis as a reset clock), are new; the sections
after §4 are renumbered. §9 records a third demo defect, the plane-wave lock
on a periodic box. Open items Q-SP5 to Q-SP9 are added. Part E of
`src/demo_force_blind_sea.py` verifies §5, and `src/scan_dark_reset.py` §8.

**Second addendum (September 2026).** Proposition Q3's recombination rate
$`\kappa`$ is named for what it is, step 16's contact sink, which is not in
the model. §6.1, Proposition Q7 (the clocks of the sea's population), and
§6.2, Proposition Q8 (immediate contact recombination as garbage
collection), are new, with open item Q-SP10 on gray pairs. Part F of
`src/demo_force_blind_sea.py` verifies both.

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
  to the straight-line (Glauber) eikonal phase $`-\int V\,dt/\hbar`$
  accumulated along its unbent path, two aligned pairs of one row wind at
  exactly $`U`$, and Theorem L8's drift term vanishes at all times. Under (S)
  a body's offset from its row's plane wave begins at $`V''`$, which is the
  row mixing step 22 measured.
- **Proposition Q5.** The sea computes the kernel's Fourier transform by
  interference: the channel's Fourier phase is the daughter row's P2 lever,
  the jump rate is the rate of change of the row's interference sum read
  against the daughter, and the sea's density is the one at which the
  channel basis is orthogonal on the reach (§5).
- **Proposition Q3.** What (S′) changes is the sea itself: deficits stay in
  their rows and only translate (§6).
- **Proposition Q7.** The sea's population runs on its own clocks: the event
  rate and the rate at which force-feeling bodies slip across force-blind
  rows are of the same order by construction, so under (S′) the sea's
  deficits are shallow but nothing refills them; dark catalysis, with
  $`\Delta S = 0`$, cannot touch them (§6.1).
- **Proposition Q8.** An immediate contact sink is garbage collection: the
  ensemble stays minimal, $`N = \lvert W\rvert`$, absorption becomes
  redundant for the bodies, and the sea pays exactly half the growth of
  $`\lVert W\rVert_1`$ whichever route is taken; only where it pays
  differs (§6.2).
- **Proposition Q4.** Under (S′) a tagged body crosses the Eckart barrier
  inside a dark pair, through the summit. Individual crossing is still not
  tunnelling (§7).
- **Proposition Q6.** Left alone, the (S′) sea relaxes to each row's eikonal
  wave, exactly. Dark catalysis is the reset clock that puts the free
  reference back: its rate is the bandwidth of the phase reference, gauge
  invariance selects its relative resets over absolute ones, and near the
  Eckart barrier three times the kernel's own rate nearly suffices (§8).
- Three demo defects, repaired (§9): the destination-row species mask behind
  step 22's open item L-SP8, non-conservation of tags in step 20's Part E,
  and the plane-wave lock on a periodic box.

**Leaves open** the corrections to steps 15 to 22 listed in §10, to be made
once (S′) is confirmed, and the items of §11.

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

**Partners, and what the sea is.** An aligned pair is a positon and a
negaton at one point, with one momentum and one clock; its members are never
apart. When a vertex reads the sea it compares *two different aligned pairs*
on its momentum row, the legs of a chord at $`x \mp y`$. These, as in step
22, are the *partners*. The sea is not a lattice. It is a random
configuration of aligned pairs with mean density $`B = 1/\pi\hbar`$ per unit
area of phase space and local density $`s(x, p, t)`$, which events change
(an absorption adds a pair on the parent's row, an emission removes one) and
which under (S′) is carried along rows (Proposition Q3). A chord's
half-separation $`y`$ is whatever its two partners happen to have: a vertex
reads every chord its row offers within the reach, and the transform over
$`y`$ is sampled by the separations present.

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
while under (S) the reading was 0.45; §8 refines these numbers. The reason is Proposition Q2: under
(S) the force fans a row out over its neighbours, and each row gathers
bodies carrying phases locked to the rows they came from.

**Costs.** A pair's energy $`p^2/2m + V(q)`$ is not conserved while it
crosses a potential. It is unobservable, because the pair is $`W`$-null, but
it is a real departure from Newtonian bookkeeping for paired bodies. And
identity now passes through classically forbidden regions in the dark
(§7): a claim about an unobservable, which leaves $`E`$ alone.

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
  (§7).

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

*so it stays on its row's plane wave up to the straight-line (Glauber)
eikonal phase $`-\hbar^{-1}\int V(q(t'))\,dt'`$ accumulated along its
unbent path. A motionless pair would keep the plane wave only with
$`\hbar\dot\theta = -p^2/2m`$, a Hamiltonian clock, not P1's.*

*(b) Any two aligned pairs on one row — each a co-located positon and
negaton sharing one clock — keep whatever separation $`2y`$ they have, and a
reader with momentum $`P(t)`$ sees them wind at*

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
demo, and (a) against the particle model in §8.2: the misalignment of every
chord of untouched sea clocks matches the closed-form eikonal to a
correlation of 1.0000.

**Which eikonal.** The word has two senses, and (S) and (S′) carry one each.
In geometrical optics the eikonal is the full phase function $`S`$, with
$`\lvert\nabla S\rvert^2 = n^2`$ (Bruns 1895); for a particle it is
Hamilton's principal function along the *true, curved* path, which is what
P1 clocks carry under (S): Definition L0's phase field $`S/\hbar`$. In
high-energy scattering the eikonal is the phase a plane wave gains from the
potential along the *unbent* straight line,

```math
\psi(\mathbf b, z) \approx e^{ikz}\exp\!\Big[-\frac{i}{\hbar v}\int_{-\infty}^{z} V(\mathbf b, z')\,dz'\Big]
```

(Molière 1947; Glauber 1959), with $`\mathbf b`$ the impact parameter; with
$`dz' = v\,dt'`$ this is Q2(a)'s $`-\hbar^{-1}\int V\,dt`$. There it is an
approximation, valid when $`E \gg \lvert V\rvert`$ so that paths barely
bend. Under (S′) it is exact for the sea's clocks by postulate, because their
paths do not bend at all. So (S) clocks carry the optical (Hamilton–Jacobi)
eikonal along bending trajectories and (S′) clocks the straight-line one:
the contrast Q2(c) draws between rows that fan out and rows that are
conserved.

Part (b) is the reason step 22's row-keeping sea held the reading: under
(S′) the own-frame reading of Theorem L7 becomes the inertial reader and
reads $`U`$, the parent frame reads $`U_{\rm res}`$ with no drift, and
Proposition L9's bound $`dp^2/m`$ becomes zero. Part (c) is the reason (S)
did not: a curved potential fans each row out over its neighbours.

---

## 5. Proposition Q5: what the sea computes

The compensated algorithm sets its jump rates from a Fourier transform of
the potential (Theorem L1):

```math
K_q(x) = -\frac{B\,dp}{\hbar}\int_{-y_{\max}}^{y_{\max}} dy\;
U_{\mathrm{res}}(x,y)\,w(y)\,\sin\!\left(\frac{2\xi_q y}{\hbar}\right),
```

with $`\xi_q = q\,dp`$ the momentum transfer of channel $`q`$, $`dp`$ the row
width, $`w`$ the reach horizon and $`B = 1/\pi\hbar`$ the sea density. The
question is whether the sea computes that transform rather than being
handed it. It supplies three of the four ingredients directly: the
quadrature nodes (the partners' positions sample $`y`$), the integrand (the
partners' clock rates give $`\hbar\dot\mu = U_{\rm res}`$ at re-lock, exactly
under (S′), Proposition Q2), and, as (a) below shows, the Fourier phase. The
fourth, the sum over chords, is an interference.

Write $`\mu^{(P)}_{ij} = \theta_i - \theta_j + P\,d/\hbar`$ for the
misalignment of partners $`i, j`$ read with lever $`P`$, where
$`d = x_j - x_i = 2y`$.

**Proposition Q5.** *(a) Reading two partners against daughter row
$`p + \xi_q`$'s plane wave is reading them against the parent's with the
channel's Fourier phase added:*

```math
\mu^{(p)}_{ij} + \frac{\xi_q d}{\hbar} = \mu^{(p+\xi_q)}_{ij} .
```

*(b) Let*

```math
z_q = \sum_{\text{chords}} w(y)\, e^{\,i\mu^{(p+\xi_q)}_{ij}}
```

*be the interference sum over the chords of row $`p`$ whose midpoints lie
in a cell of length $`\Delta x`$ at $`x`$. At re-lock,
$`\mu^{(p)}_{ij} = 0`$, with the parent-frame rate of Proposition Q2(b),
and in the continuum limit of a sea of uniform density $`B`$,*

```math
K_q(x) = \frac{1}{B\,dp\,\Delta x}\,\frac{d}{dt}\,\mathrm{Re}\,z_q .
```

*(c) The channel functions $`\sin(2\xi_q y/\hbar)`$ are orthogonal on the
reach $`[-y_{\max}, y_{\max}]`$ for every $`q \ne q'`$ exactly when
$`dp\,y_{\max}`$ is a whole multiple of $`\pi\hbar/2`$. The shortest such
reach is Theorem L1's strip, $`B\,dp\,2y_{\max} = 1`$.*

*Proof.* (a) Substitute the definition. (b) Chords with midpoint in
$`[x, x + \Delta x]`$ and half-separation in $`[y, y + dy]`$, $`y > 0`$,
number $`2(B\,dp)^2\,\Delta x\,dy`$, the 2 being the Jacobian of
$`(x_i, x_j) \mapsto (x, y)`$. By (a) the weight of a chord is
$`\cos(2\xi_q y/\hbar + \mu^{(p)})`$, so at $`\mu^{(p)} = 0`$

```math
\frac{d}{dt}\,\mathrm{Re}\,z_q
= -\frac{2(B\,dp)^2\Delta x}{\hbar}\int_0^{y_{\max}} U_{\rm res}\,w\,
\sin\!\left(\frac{2\xi_q y}{\hbar}\right) dy
= -\frac{(B\,dp)^2\Delta x}{\hbar}\int_{-y_{\max}}^{y_{\max}} (\cdots)\,dy
= B\,dp\,\Delta x\,K_q ,
```

the integrand being even in $`y`$. (c) With $`k = 2\,dp/\hbar`$,

```math
\int_{-Y}^{Y} \sin(qky)\,\sin(q'ky)\,dy
= \frac{\sin\big((q - q')kY\big)}{(q - q')k} - \frac{\sin\big((q + q')kY\big)}{(q + q')k},
```

which vanishes for every $`q \ne q'`$ if $`kY`$ is a multiple of $`\pi`$;
the pair $`q = 1, q' = 2`$ gives $`3\sin kY = \sin 3kY`$, so $`\sin kY = 0`$
is also necessary. $`kY = \pi`$ is $`dp\,Y = \pi\hbar/2`$, which with
$`B = 1/\pi\hbar`$ is $`B\,dp\,2Y = 1`$. $`\square`$

*Verification.* Part E of `demo_force_blind_sea.py`. With SymPy: the
identity (a) holds exactly; on $`dp\,Y = \pi\hbar/2`$ the Gram matrix of the
first four channels, divided by $`Y`$, is the identity, and on
$`dp\,Y = 0.4\pi\hbar`$ its largest off-diagonal element is 0.288.
Numerically, (b) against the code's kernel at $`x = -0.625`$, with the
chord sum as a fine quadrature in $`y`$:

| $`q`$ | 1 | 2 | 3 | 4 | 6 | 8 |
|---|---|---|---|---|---|---|
| $`K_q`$ (code) | 1.076080 | 0.141430 | −0.174147 | −0.056715 | −0.020788 | −0.002426 |
| relative difference | $`5.5\times10^{-7}`$ | $`8.4\times10^{-6}`$ | $`1.0\times10^{-5}`$ | $`4.2\times10^{-5}`$ | $`1.8\times10^{-4}`$ | $`2.1\times10^{-3}`$ |

the same discretisation error as Theorem L1's own check.

**Reading it.** The Fourier phase is not an extra factor: channel $`q`$ is
the sea read against the daughter row's P2 plane wave, which is step 22's
answer to L-SP6 applied to the daughter as reader. The daughter at
$`p - \xi_q`$ gives the opposite sign, matching the $`\pm 1`$ deposits. So
the jump rate of channel $`q`$ is the rate at which the row's interference
pattern, seen from the daughter row, changes as the potential winds the
clocks. The Wigner transform itself is built the same way: the density
matrix read along chords against a reader's own plane wave (step 22 §8).

The analogy is Fraunhofer diffraction, in which the far-field amplitude in
a direction is the Fourier transform of the illumination across an aperture
(Goodman 2005, chapter 4). The aperture is the set of partners within the
reach, the illumination phases are their clocks, the far-field directions
are the daughter rows, and the potential is a phase screen that modulates
the aperture in time. A phased-array antenna is the engineering version:
nobody computes the transform; the superposition does.

**What it does not do.**

- *It is critically sampled, so it is statistical.* By (c) the sea holds, on
  average, one aligned pair per row in the strip on which the channel basis
  is orthogonal. Any single reading is one sample of a random configuration
  (§1), and the transform emerges only in the average: the contact cost of
  Proposition L2.
- *It must be normalised by the actual sea.* (b) holds for a sea of uniform
  density $`B`$. The number of chords at separation $`2y`$ scales as
  $`s(x - y, p)\,s(x + y, p)`$, so a reading normalised to $`B`$ is biased
  wherever the sea departs from it: Proposition L2(a)'s weighting, and under
  (S′) the departure is larger (Proposition Q3, §9).
- *It needs a fresh reference.* (b) holds at re-lock. Between re-locks the
  weights carry the misalignment the potential has wound into the clocks,
  and §8 is about that.
- *The horizon is imposed, not computed.* Which partners count, and with what
  taper $`w`$, is a rule (Q-SP6).
- *Summing amplitudes is a vertex operation.* Superposition is cheaper than a
  computed transform, but the ontology must still grant it at the vertex, in
  the spirit of step 4's vertex weight P5.

---

## 6. Proposition Q3: what (S′) does to the sea

Under (S′) a row of sea translates rigidly at its own velocity, so a deficit
or excess made at a parent's row stays in that row. Under (S) the force
shears it across rows near the barrier. The step 16 ledger under both
motions (`demo_force_blind_sea.py` Part C):

Emissive unravelling with step 16's contact sink at rate $`\kappa`$, which
is not part of the model (§6.1), in Theorems S4 and S5's setting (step 16
Part C), worst cell $`s/B`$ at $`T = 6`$:

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
across rows is not refilled. With step 16's contact sink, which is outside
the model, a fast sink ($`\kappa \ge 200`$) keeps the worst cell
non-negative under (S′), where under (S) it does not; a slow one
($`\kappa = 20`$) does the reverse (§6.1 explains both). A rate throttled by the
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

### 6.1 Proposition Q7: the clocks of the sea's population

*Which recombination.* The $`\kappa`$ of the first table is step 16's
*contact* recombination: a bilinear sink that removes a positon and a
negaton together from one cell at rate $`\kappa\,n_+n_-`$ (Theorem S1's
form). It is not in the model. Postulate (D), with $`f = 1/2`$, is the
sinkless case (Theorem N3), the model's only pair-forming channel is the
absorptive event, and the algorithm specification keeps $`\kappa`$ only as
an optional setting (ORIENTATION §8). What the first table and Proposition
Q3 say about $`\kappa`$ is therefore about a sink the model does not have.
Nor can dark catalysis act on the sea's population: by Theorem L5 it has
$`\Delta E = \Delta N = \Delta S = 0`$. The population and the phase
reference of §8 are separate ledgers, each with its own clocks.

*The split.* Write the sea as

```math
s = B - D + C ,
```

with $`D`$ the cumulative debit (an emission, at the parent's cell) and
$`C`$ the cumulative credit (an absorption at the parent's cell, or a
contact recombination where it happens), both carried with the sea's own
motion. The identity holds to better than $`10^{-12}`$ in every run below
(`demo_force_blind_sea.py` Part F).

*The clocks at the Eckart flank* ($`dp = 0.25`$, $`x = -0.625`$):

| clock | rate | acts on |
|---|---|---|
| event rate $`\Gamma`$ | 3.13 | counts: absorptions credit, emissions debit, both at the parent's row |
| slip of force-feeling bodies across force-blind rows, $`\lvert V'_{\rm eff}\rvert/dp`$ (zero under (S)) | 3.07 rows per unit time | counts: which sea patch a parent debits |
| contact sink $`\kappa\lvert E\rvert`$ (outside the model), median over the packet | 1.3, 13, 133 at $`\kappa`$ = 20, 200, 2000 | counts: credit where excess bodies recombine |
| winding $`\lvert U\rvert/\hbar`$ | at most $`V_0/\hbar = 1`$ | phases (§8) |
| dark catalysis $`\kappa_d\Gamma`$ | about 9 needed (Q6) | phases only |

**Proposition Q7.** *(a) The event rate and the slip rate are of the same
order by construction: over the grids and barriers below, the median of
$`\Gamma/(\lvert V'_{\rm eff}\rvert/dp)`$ where $`\Gamma`$ exceeds a fifth of
its maximum lies between 0.70 and 1.13, and both scale as $`V`$ and as
$`1/dp`$. (b) In the model's own, sinkless ledger the sea's excursions are
therefore never separated in time scale from the traffic that makes them.
Under (S) a sea patch rides the same flow as the parents above it, is
debited while partners are scarce and refilled by later absorptions on the
same patch; under (S′) the parents slide across the rows about as fast as
they fire, so the debit is spread over many patches and nothing refills
them. (c) With the contact sink, the sign of the (S′) worst cell follows
the ratio of $`\kappa\lvert E\rvert`$ to the slip rate.*

| $`n_p`$ | $`dp`$ | $`a`$ | $`V_0`$ | $`y_{\max}`$ | max $`\Gamma`$ | max $`\lvert V'_{\rm eff}\rvert/dp`$ | median ratio |
|---|---|---|---|---|---|---|---|
| 64 | 0.25 | 1 | 1 | 6.28 | 3.130 | 3.072 | 1.03 |
| 128 | 0.125 | 1 | 1 | 12.57 | 6.575 | 6.144 | 1.13 |
| 32 | 0.5 | 1 | 1 | 3.14 | 1.126 | 1.537 | 0.70 |
| 64 | 0.25 | 2 | 1 | 6.28 | 1.129 | 1.536 | 0.72 |
| 64 | 0.25 | 1 | 3 | 6.28 | 9.389 | 9.217 | 1.03 |

The reason for (a): when the reach is long against the potential's width,
$`V(x \pm y)`$ has decayed over most of it, and the residual
$`U_{\rm res} = U - 2yV'_{\rm eff}`$ is dominated there by the subtracted
force term. The compensated event rate is then, to order one, the rate at
which the classical force would carry a body across one row. This is the
same growth with the reach as step 11b's Theorem E8.

The worst cell of each run, with its debit and credit and where it sits:

| ledger | sea | worst $`s/B`$ | at $`t`$ | $`D/B`$ | $`C/B`$ | $`x`$ | $`p`$ |
|---|---|---|---|---|---|---|---|
| sinkless (the model) | (S) | −0.339 | 1.62 | 1.394 | 0.055 | +1.56 | +1.375 |
| sinkless (the model) | (S′) | −0.060 | 3.52 | 1.462 | 0.402 | +3.75 | +1.125 |
| contact, $`\kappa = 20`$ | (S) | +0.006 | 5.82 | 2.532 | 1.538 | +7.19 | +1.375 |
| contact, $`\kappa = 20`$ | (S′) | −0.673 | 6.00 | 5.674 | 4.001 | −0.31 | +0.125 |
| contact, $`\kappa = 200`$ | (S) | −0.145 | 5.94 | 1.972 | 0.827 | +7.19 | +1.375 |
| contact, $`\kappa = 200`$ | (S′) | +0.005 | 5.78 | 2.159 | 1.165 | +4.38 | +0.875 |
| contact, $`\kappa = 2000`$ | (S) | −0.141 | 5.46 | 1.933 | 0.792 | +6.56 | +1.375 |
| contact, $`\kappa = 2000`$ | (S′) | +0.073 | 5.46 | 2.107 | 1.181 | +4.06 | +0.875 |

Reading the table. In the sinkless ledger under (S) the worst cell is the
patch travelling with the packet's momentum peak, debited early (K9's
bootstrap, when partners are scarce) with almost no credit; later
absorptions on the same patch refill it to about $`+0.17\,B`$. Under (S′)
the dip is a sixth as deep, and a quarter of its debit has already been
credited back, but by other parents. With the contact sink at
$`\kappa = 20`$, where $`\kappa\lvert E\rvert \approx 1.3`$ is below the
slip rate of 3, excess bodies slip away before they recombine, and the
worst cell is the slowest row, $`p = +0.125`$, parked under the emitting
region: debited $`5.7\,B`$ and credited $`4.0\,B`$. At $`\kappa \ge 200`$,
where $`\kappa\lvert E\rvert`$ is well above the slip rate, the credit lands
on the debited patch and the worst cell stays non-negative. Under (S) there
is no slip at all; the credit lands on the daughter rows, which are other
flow lines, and the patch under the packet's momentum peak drains at
$`\tfrac12\big(\Gamma\lvert E\rvert - \lvert K\rvert\ast\lvert E\rvert\big)`$
for the whole encounter, whatever $`\kappa`$ is. At $`\kappa = 20`$ the slow
credit rides the same flow as the debit and arrives late but in the right
place, which is why $`\kappa = 20`$ favours (S) and $`\kappa \ge 200`$
favours (S′).

The live question in the model's own terms is therefore not whether a sink
would help but Q-SP2: whether supply for emission becomes limiting in
longer or deeper runs, where under (S′) nothing refills the shallow
deficits.

### 6.2 Proposition Q8: immediate contact recombination as garbage collection

Suppose the contact sink were immediate: after every event step each cell
keeps only its majority species, $`u_\pm = \max(\pm W, 0)`$, and the removed
pairs are credited to the sea in that cell.

*It does not make every event emissive.* An absorptive event needs
opposite-species partners at both daughters: the $`+1`$ deposit at
$`p + \xi`$ is realised by removing a negaton there, the $`-1`$ at
$`p - \xi`$ by removing a positon (Theorem S6). With only the majority
species left in each cell, absorption survives exactly where both deposits
oppose the local sign of $`W`$, so that the event reduces $`\lvert W\rvert`$
at both ends, and in amounts capped by the smaller $`\lvert W\rvert`$ there.
Where $`W \ge 0`$ every event is emissive, as it already is at $`t = 0`$ in
the minimal ensemble without any sink (Theorem K9); as interference makes
$`W`$ negative in places, absorption returns.

*Nor is absorption needed.* Remove it, keeping only emission and the
immediate sink, and the sink becomes garbage collection: the same-cell
annihilation step of the signed-particle Wigner Monte Carlo method
(Nedjalkov et al. 2004; Sellier 2015), with the sea's pair ledger added.

**Proposition Q8.** *With an immediate contact sink, (a) the ensemble is
minimal, $`N = \lvert W\rvert`$, at every step, with or without absorption,
so the growth of the body count that ruined the emissive unravelling
(Theorems S4 and S8) is gone; (b) the realisation (absorptive first, or emission only) is
invisible to the bodies and to $`E`$; (c) the sea's total is fixed by the
pair-count identity $`P = S + N/2`$ (Theorem J1):*

```math
S_{\rm tot}(t) = S_{\rm tot}(0) - \tfrac12\Big(\lVert W(t)\rVert_1 - \lVert W(0)\rVert_1\Big),
```

*so the sea pays exactly half the growth of $`W`$'s $`L^1`$ norm, which only
negativity can make; (d) the realisation decides only where the sea is
credited. Absorption re-forms the pair on the parent's row; emission
followed by contact recombination debits the parent's row and credits the
daughter rows, the relocation of Theorem S4.*

*Proof of (a) to (c).* (a) is the sink's definition. (b) Both routes change
$`E`$ by the kernel's deposits, and the sink preserves $`E`$ and sets
$`N = \lvert E\rvert`$. (In a time step that fires the channels in
sequence with live populations, the two routes see slightly different
parents after the first channel, so $`E`$ agrees only to the order of the
step: $`7.8\times10^{-3}`$ against $`8.2\times10^{-3}`$ below.) (c) Every process in the ledger conserves
$`P = S + N/2`$, the sink included (it removes two bodies and adds one
pair), so $`\Delta S_{\rm tot} = -\tfrac12\Delta N = -\tfrac12\Delta\lVert W\rVert_1`$. $`\square`$

*Measured* (Part F; step 16's event-resolved ledger, $`T = 6`$):

| events | sea | $`f`$ | $`E`$ against the QLE | $`N/\lvert W\rvert`$ | sink per event | worst $`s/B`$ | final min $`s/B`$ | identity (c) |
|---|---|---|---|---|---|---|---|---|
| absorptive first + immediate sink | (S) | 0.347 | $`7.8\times10^{-3}`$ | 1.0001 | 0.20 | −0.169 | −0.130 | $`3\times10^{-7}`$ |
| absorptive first + immediate sink | (S′) | 0.347 | $`7.8\times10^{-3}`$ | 1.0001 | 0.20 | **+0.401** | +0.455 | $`3\times10^{-7}`$ |
| emission only + immediate sink | (S) | 0 | $`8.2\times10^{-3}`$ | 1.0001 | 0.90 | −0.186 | −0.186 | $`2\times10^{-7}`$ |
| emission only + immediate sink | (S′) | 0 | $`8.2\times10^{-3}`$ | 1.0001 | 0.90 | −0.138 | −0.103 | $`2\times10^{-7}`$ |

The sinkless model, for comparison, has $`f = 0.400`$ and $`E`$ against the
QLE at $`5.8\times10^{-3}`$ (Part F; the worst-cell table of §6.1 has its
sea). Over the run
$`\lVert W\rVert_1`$ grows to 1.70 times its initial value; about a fifth of
it ends up negative.

Two cautions on these numbers. The immediate sink's $`\max(\pm W, 0)`$ has
kinks, and spectral transport rings them into negative partner caps; the
runs therefore floor the caps at zero, and without that the absorptive
fraction comes out negative. The same floored update without the sink does
not reproduce the published sinkless run ($`f = 0.448`$, with $`E`$ off the
QLE by 12 per cent), so the sinkless comparison above is the published
update, not a matched rerun. And the two ledgers of this
note differ: these rows are event-resolved, while the contact-sink rows of
§6.1 are step 16's mean-field ledger, so the $`\kappa = 2000`$ row there is
not the immediate limit of these.

*What garbage collection costs.* By N3's sum rule,
$`\Gamma_{\rm tot}(1 - 2f) = R_{\rm sink}`$, any sink moves $`f`$ below one
half, so the absorptive fraction becomes a property of the representation
and (D)'s $`f = 1/2`$ loses its content. "Immediately" needs a co-location
rule: on the mesh "same cell" is automatic, but with continuous positions it
is a matching-cell choice, which step 21 found to be a parameter with costs
(Proposition M9).

*A join, not an end.* In the conserved-body reading (Theorem Y5's piecewise
worldlines) a body that meets an opposite-signed body in a cell does not
end: it joins a pair, and its worldline kinks. Under (S′) the kink has two
parts: a jump to the common point $`(\bar x, \bar p)`$, at most a cell in
$`x`$ and in $`p`$, which conserves total momentum; and a change of force
law, from $`\dot p = -V'`$ to $`\dot p = 0`$, a change of direction in phase
space with momentum continuous. Under (S) only the jump would remain. What
is lost is the determinacy of pairing, not continuity: with several bodies
of each species in a cell, the counts do not say which joins which, nor
which pair a later ionisation draws from, so identity survives statistically
(step 20's proportional allocation; Proposition Q4 is this picture
measured).

*Gray on joining.* A pair formed on contact brings two unrelated clocks and
is generically gray, not dark: step 22 §5 measured a mean $`\lvert\mu\rvert`$
of 1.40 to 1.48 for recombination under phase continuity, against
$`\pi/2`$ for random phases. Darkness needs a phase alignment, by the
re-lock rule or by dark catalysis. That leaves (S′) ambiguous for gray
pairs. Its first argument (§2: the force has nothing to act on, a
$`W`$-null pair) applies to any co-located, equal-momentum pair; its second
(darkness as a limit) applies only to dark pairs. Whether a pair becomes
force-blind on joining or only on alignment is open (Q-SP10). With garbage
collection, pairs form at 0.9 per event, so the first reading would inject
force-blind gray pairs into the phase reference at the event rate, for dark
catalysis to clean up.

---

## 7. Proposition Q4: identity crosses in the dark

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

## 8. Proposition Q6: dark catalysis as a reset clock

Dark catalysis (Theorem L5) is the sea's one E-neutral interaction: an
aligned pair A fires at the kernel's rate and a partner pair C in A's row
re-forms aligned, on A's transported phase or on the mean of the two. Under
(S) step 22 used it to fight three things at once: the within-row beat of
the initial lock, the tidal row mixing of Proposition Q2(c), and the
misalignment the potential winds into the clocks. It lost, reading the
kernel at 0.45 against a control near 0.89 (step 22 §9). Under (S′) rows
are conserved, so the first two are gone by construction, and a row-centred
sea has no within-row beat. This section asks what is left for it to do,
and how fast it must do it.

### 8.1 What it does not do

By Proposition Q1 dark catalysis moves nothing in the ledger under either
motion: it creates and destroys nothing, and the ledger never reads the
sea. Its rate is therefore free as far as $`E`$, $`N`$ and $`f`$ are
concerned. It matters only to the constructive programme of step 22 and of
Proposition Q5, in which the kernel is read from the sea's clocks instead of
being precomputed from $`V`$.

### 8.2 What the sea does if left alone

A clock read with its row's lever winds, between resets, at the inertial
reader's rate of Proposition Q2(b): $`\hbar\dot\mu^{(p)}_{ij} = U(x, y)`$.
Integrated along the straight paths of (S′), the misalignment of a chord
whose clocks were last set at $`t = 0`$ is

```math
\mu^{(p)}_{ij}(t) = -\frac{1}{\hbar}\left[\int_0^t V\big(x_i(t')\big)\,dt' - \int_0^t V\big(x_j(t')\big)\,dt'\right],
```

the difference of the two partners' straight-line eikonals (Proposition
Q2(a)). For the Eckart barrier the integral has a closed form,
$`\int_0^t \mathrm{sech}^2(x_0 + vt')\,dt' = [\tanh x(t) - \tanh x_0]/v`$,
continued periodically on the box. In the particle model with the sea
force-blind, row-centred and left alone (`demo_sea_lock_particles.py
--sea-force blind --sea-p rows --wrap-phase --no-events --readers --filter`,
seed 11), the misalignment of every chord of untouched sea clocks near the
barrier matches this closed form with correlation 0.9989 at $`t = 0.1`$,
0.9999 at $`t = 1`$ and 1.0000 from $`t = 4`$ to $`t = 25`$, and regression
slope 1.000 within 0.005 throughout.

So a force-blind sea left alone does not stay on the free plane wave: each
row relaxes to its own eikonal wave. For a row moving at $`v \ne 0`$, a
chord whose partners entered from far upstream carries
$`\mu_{ij} = (\hbar v)^{-1}\int_{x_i}^{x_j} V\,dx'`$, bounded by
$`\int V\,dx/\hbar\lvert v\rvert`$; the row at rest winds at $`U/\hbar`$
without bound. Read against this eikonal reference, the kernel near the
barrier reads 0.778 against a control of 0.880 (table below). Proposition
Q5(b) needs the *free* reference, $`\mu^{(p)} = 0`$. What dark catalysis is
for, under (S′), is to put it back.

### 8.3 The reset filter

Idealise the resets of a chord as Poisson at rate $`\Gamma_\ell`$, so that
the time $`\tau`$ since the last reset is distributed as
$`\Gamma_\ell e^{-\Gamma_\ell\tau}`$, and freeze the winding rate over one
interval, $`\mu = \omega\tau`$ with $`\omega = U/\hbar`$. A chord enters
the reading of channel $`q`$ with weight $`\sin(a + \mu)`$,
$`a = 2\xi_q y/\hbar`$, and

```math
\big\langle e^{\,i(a + \omega\tau)}\big\rangle = e^{ia}\,\frac{\Gamma_\ell}{\Gamma_\ell - i\omega}
\quad\Longrightarrow\quad
\big\langle \sin(a + \omega\tau)\big\rangle
= \frac{\Gamma_\ell^2\sin a + \Gamma_\ell\,\omega\cos a}{\Gamma_\ell^2 + \omega^2} .
```

**Proposition Q6(a).** *Under Poisson resets at rate $`\Gamma_\ell`$ the
sea's reading of $`K_q`$ departs from the kernel in two ways: an amplitude
compression by $`\Gamma_\ell^2/(\Gamma_\ell^2 + \omega^2)`$, of second order
in $`\omega/\Gamma_\ell`$, and a quadrature term of first order,*

```math
\delta K_q \;\propto\; -\,\langle\tau\rangle\sum_{\text{chords}} U_{\rm res}\,w\,\frac{U}{\hbar}\,\cos a ,
\qquad \langle\tau\rangle = 1/\Gamma_\ell .
```

*The quadrature term is even in $`\xi_q`$ and its weight $`U\,U_{\rm res}`$
is not sign-definite.*

It is the phase-interruption mathematics of collision broadening (Van Vleck
and Weisskopf 1945) and of motional narrowing (Bloembergen, Purcell and
Pound 1948; Anderson 1954), with one difference that matters: there the
phase winds at the frequency the line is read at, and interruption broadens
the line symmetrically, which is a decoherence. Here the clocks wind at
$`U/\hbar`$ relative to their row while the reading multiplies
$`U_{\rm res}`$. The product is not sign-definite, so the quadrature term is
a spurious even (cosine) component of the reading, not a decoherence rate.
An event realisation that deposits $`\pm 1`$ by construction would take it
as an error in the odd kernel, $`\pm\delta K_{\lvert q\rvert}`$, whose first
moment $`2\sum_{q>0}\xi_q\,\delta K_q`$ need not vanish: a spurious force of
order $`\langle\tau\rangle`$ (Q-SP7).

### 8.4 Relative and absolute resets

The simplest reset is **absolute**: set the clock to its row's free plane
wave, $`\theta \leftarrow (px - p^2t/2m)/\hbar`$. That needs an external clock
and the zero of $`V`$. Under $`V \to V + c`$ every P1 clock winds faster by
$`-c/\hbar`$ while the reset target does not move, so two clocks last reset
$`\tau_i`$ and $`\tau_j`$ ago acquire the misalignment
$`c(\tau_j - \tau_i)/\hbar`$: noise that depends on an unobservable constant.

Dark catalysis is **relative**: it compares clocks with clocks. In the
reach form, C's clock and A's clock transported to C with A's momentum are
both replaced by their mean. A uniform shift of every clock's rate changes
neither their difference nor, up to a common phase, their mean, so the reset
is invariant under $`V \to V + c`$. The price is that it is partial: it
equalises A and C, and pulls the row toward consensus only on average, as
gossip averaging does (Boyd et al. 2006).

**Proposition Q6(b).** *Of the two, only relative resets are invariant
under $`V \to V + c`$.* $`\square`$ (By the two paragraphs above.)

### 8.5 How fast it must be

`src/scan_dark_reset.py` runs the particle model on the force-blind,
row-centred sea (`--sea-force blind --sea-p rows`), with the lock kept
continuous across the periodic boundary (`--wrap-phase`, §9), three seeds,
late means over $`4 \le t \le 25`$. The reading is taken over sea clocks
only, one member per aligned pair, near the barrier. $`\tau`$ is the
regression slope of the misalignment on $`U/\hbar`$, the effective reset
time; $`b`$ is the fitted coefficient of the quadrature term of Q6(a) in
the difference between the reading and its $`\mu \equiv 0`$ control. The
dark rate is given as a multiple $`\kappa`$ of the kernel's own rate.

| regime | reading | control | $`\tau`$ | rms residual $`\mu`$ | $`b`$ | coherence | median clock age |
|---|---|---|---|---|---|---|---|
| left alone (eikonal sea) | 0.778 ± 0.011 | 0.880 | 0.045 | 0.795 | 0.029 | 0.769 | 14.6 |
| relative (dark catalysis), $`\kappa = 1`$ | 0.863 ± 0.002 | 0.880 | 0.224 | 0.360 | 0.194 | 0.935 | 0.40 |
| relative, $`\kappa = 3`$ | 0.875 ± 0.001 | 0.880 | 0.161 | 0.175 | 0.138 | 0.983 | 0.13 |
| relative, $`\kappa = 10`$ | 0.879 ± 0.001 | 0.880 | 0.075 | 0.071 | 0.069 | 0.997 | 0.04 |
| relative, $`\kappa = 30`$ | 0.879 ± 0.001 | 0.880 | 0.040 | 0.043 | 0.034 | 0.999 | 0.00 |
| absolute, $`\kappa = 1`$ | 0.871 ± 0.001 | 0.880 | 0.172 | 0.267 | 0.141 | 0.964 | 0.57 |
| absolute, $`\kappa = 3`$ | 0.879 ± 0.001 | 0.880 | 0.110 | 0.098 | 0.096 | 0.994 | 0.17 |
| absolute, $`\kappa = 10`$ | 0.880 ± 0.001 | 0.880 | 0.051 | 0.036 | 0.037 | 0.999 | 0.04 |
| absolute, $`\kappa = 30`$ | 0.880 ± 0.001 | 0.880 | 0.019 | 0.018 | 0.011 | 1.000 | 0.01 |

With events as well (`--relock-w 3`), the sea-only reading and its control:

| regime | reading | control | $`\tau`$ | rms residual $`\mu`$ | all-bodies reading | all-bodies control |
|---|---|---|---|---|---|---|
| no dark catalysis | 0.784 ± 0.012 | 0.836 | 0.108 | 0.517 | 0.803 ± 0.014 | 0.885 |
| relative, $`\kappa = 1`$ | 0.809 ± 0.010 | 0.833 | 0.149 | 0.375 | 0.828 ± 0.010 | 0.893 |
| relative, $`\kappa = 3`$ | 0.829 ± 0.007 | 0.837 | 0.141 | 0.242 | 0.844 ± 0.005 | 0.887 |
| relative, $`\kappa = 10`$ | 0.835 ± 0.007 | 0.838 | 0.070 | 0.156 | 0.850 ± 0.007 | 0.891 |
| relative, $`\kappa = 30`$ | 0.831 ± 0.009 | 0.832 | 0.035 | 0.152 | 0.845 ± 0.006 | 0.888 |

![The reading and the effective reset time against the dark rate](https://raw.githubusercontent.com/billpage/wpmw/output/figures/sea_reset_filter.png)

*Left: the sea-only reading near the barrier against the dark rate, for
relative resets without and with events and for absolute resets, with the
eikonal sea and the two controls. Right: the effective reset time
$`\tau`$.*

**Proposition Q6(c) (measured).** *(i) Relative resets at the kernel's own
rate recover most of what the eikonal reference loses (0.863 against 0.778,
control 0.880); three times the rate leaves 0.005, ten times 0.001. (ii) The
fitted quadrature coefficient tracks the effective reset time, $`b/\tau`$
between 0.84 and 0.92 for relative resets without events, as Q6(a)
predicts to first order; the quadrature term explains 21 to 41 per cent of
what the reading loses, and the rest is the random residual misalignment. (iii) At equal
event rate absolute resets do better, with $`\tau`$ smaller by a factor 1.3
at $`\kappa = 1`$ growing to 2.1 at $`\kappa = 30`$. (iv) With events the
sea-only reading reaches its own control by $`\kappa = 10`$ (0.835 against
0.838). What events cost is the control itself, 0.880 falling to about
0.836, because ionisation removes the pairs the reading samples near the
barrier.*

Two readings of the numbers. First, $`\tau`$ does not fall as
$`1/\kappa`$: resets happen at $`\kappa\Gamma(x)`$, and $`\Gamma`$ vanishes
at the summit and decays away from the barrier, while a chord's partners
lie up to a reach away. Clocks are reset mainly on the flanks. Second,
relative resets lag absolute ones increasingly with rate. A consensus
update smooths the misalignment field as a diffusion does, relaxing short
wavelengths at about the reset rate; but the eikonal phase of a row that
has passed over the barrier is a step, $`(\hbar v)^{-1}\int V\,dx`$, and a
step only spreads. That is the likely reason; it is not separately measured
(Q-SP8).

**Step 22's plateau.** Step 22 §9 found the row-keeping sea's reading
saturating near 0.85 against a control near 0.89. The second table
separates the two causes. The sea-only reading has no plateau: it reaches
its own control. The all-bodies reading of step 22 includes free bodies as
partners; they are not sea clocks, they are not reset, and at
$`\kappa = 10`$ they account for the whole gap (0.850 against 0.891). And
the sea-only control is lower with events than without, because events thin
the sea near the barrier.

### 8.6 What dark catalysis is for

Under (S′) dark catalysis is not dynamics but a **reset clock**: the dump
step of an integrate-and-dump receiver. The potential winds each clock at
its own position; left alone, the sea integrates that winding into its
eikonal wave, and the vertex needs the rate, not the integral. Resetting
turns the sea from an integrator of the potential's phase into a reader of
its rate. The reset rate is the bandwidth of the phase reference: the
misalignment it leaves enters the reading linearly in $`\langle\tau\rangle`$
(Q6(a)). Gauge invariance selects relative resets over absolute ones
(Q6(b)), at a cost of up to a factor two in $`\langle\tau\rangle`$. By
Proposition Q1 the rate is free as far as the observable is concerned; near
the Eckart barrier the kernel's own rate falls a little short and three
times it nearly suffices (Q6(c)). What fixes the rate, if the ontology is to
fix it at all, is open (Q-SP9).

---

## 9. Demo defects

Three defects in the demos, found in the course of this note. Each is
repaired, or switched off by a flag whose default keeps the published
output.

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

**The tag bookkeeping of step 20's Part E** (§7).

**The plane-wave lock on a periodic box.** The particle model of step 22
runs on a periodic box of length $`L = 48`$ with rows at whole multiples of
$`dp = 0.25`$, and locks each row to $`e^{ipx/\hbar}`$. That plane wave is
periodic on the box only if $`pL/\hbar`$ is a multiple of $`2\pi`$, and
$`dp\,L/\hbar = 12`$ is not, so every body that crosses the boundary
arrives carrying a lock defect of $`pL/\hbar \bmod 2\pi`$. Under (S′) the
defect sweeps each row as a front, once per crossing time. The flag
`--wrap-phase` keeps $`\theta - px/\hbar`$ continuous across the boundary,
which makes the box a window on an open line. On the force-blind,
row-centred sea (three seeds, late means):

| regime | without `--wrap-phase` | with it |
|---|---|---|
| left alone, sea-only reading | 0.739 ± 0.016 | 0.778 ± 0.011 |
| left alone, all-bodies rate reading (A) | 0.732 ± 0.015 | 0.771 ± 0.010 |
| events and reach dark catalysis ×10, sea-only | 0.832 ± 0.008 | 0.835 ± 0.007 |
| the same, all-bodies | 0.851 ± 0.004 | 0.850 ± 0.007 |

Resets heal the defect; without them it costs about 0.04. The step 22 §9
tables were measured without the flag, under (S) as well as under
row-keeping; the (S) rows were not re-measured.

---

## 10. Corrections to make if (S′) is confirmed

- **Step 15** ([`eckart_barrier_compensated.md`](eckart_barrier_compensated.md)).
  §8's "every member of the created sea alike follows a genuine Newtonian
  worldline under the full classical force, for its entire life" becomes a
  statement about free bodies; stretches spent in a pair are straight. K-LS7's premise
  "pairs stream, so the realised profile is $`\Gamma`$ transported by the
  classical flow" becomes transport along rows.
- **Step 16** ([`sea_population_equilibrium.md`](sea_population_equilibrium.md)).
  The results about $`E`$, $`N`$ and $`f`$ stand (Proposition Q1). The sea's
  own profile changes (Proposition Q3). S4's "fast recombination does not
  repair it" concerns the contact sink, which is outside the model; under
  (S′) it fails once $`\kappa\lvert E\rvert`$ is well above the slip rate
  (Proposition Q7). The S5 table is stale under either motion. The demo docstring's "unique motion" argument is withdrawn.
- **Step 17** ([`compensated_ontology.md`](compensated_ontology.md)). (S) is
  restated as (S′); G1, G3 and G4 stand (Proposition Q1).
- **Step 20** ([`dark_sea_and_worldline_identity.md`](dark_sea_and_worldline_identity.md)).
  Y1 to Y5 stand. Y6's table is replaced by §7's; under (S) its
  transmissions survive the tag repair to about 0.01, but its kink counts
  were inflated about twofold by the clipping. Under (S′) its mechanism
  ("never through") is withdrawn, and the "dynamically inert" row of the
  darkness table gains the force.
- **Step 22** ([`sea_phase_reference.md`](sea_phase_reference.md)). L7's
  negative result is specific to (S) partners; under (S′) the own frame reads
  $`U`$ and L9's drift term is zero (Proposition Q2). L-SP8 is answered as a
  demo defect (§9), and Proposition L2(a)'s figures change slightly. L-SP10
  is answered provisionally by this note. §5's phase-field figure depicts (S)
  transport. §9's plateau of the row-keeping sea is not incomplete
  resetting: the sea-only reading reaches its own control, and the gap is
  the free bodies read as partners (§8.5). Its tables were measured without
  `--wrap-phase` (§9).

Steps 14, 18, 19 and 21 need nothing: their uses of (S) concern bodies.

---

## 11. Open items

- **Q-SP1.** Replace the spectral transport of tags by a positivity-preserving
  one (semi-Lagrangian or particle tags) and confirm §7 without the ripple
  correction.
- **Q-SP2.** Under (S′) the sea's worst cell hovers near zero (Proposition Q3).
  Does supply for emission become limiting in longer or deeper runs, where
  (S) would recover?
- **Q-SP3.** An aligned pair's energy is not conserved under (S′) while it crosses
  a potential. Is there any ledger quantity, observable or not, that the
  project has treated as conserved and that this breaks?
- **Q-SP4.** Carry out the corrections of §10 once (S′) is confirmed.
- **Q-SP5.** The frozen comb. Under (S′) a row's sea translates rigidly, so
  its pattern of separations is frozen between events, and a parent on that
  row initially co-moves with it and sees the same quadrature nodes for a
  while. Its sampling noise is then correlated in time instead of averaging
  step by step. Measure the time correlation of the pair-sum noise under
  (S′) against (S), and what it does to the contact cost of Proposition L2.
- **Q-SP6.** The horizon. Proposition Q5 computes the transform over the
  chords a rule admits, with a taper $`w`$ the rule imposes. Could the
  horizon emerge instead, for instance from the fall of the lock's coherence
  with separation?
- **Q-SP7.** The quadrature leakage of Proposition Q6(a) is a spurious even
  component of the reading. In an event realisation that deposits
  $`\pm 1`$ by construction it becomes an error in the odd kernel with
  first moment $`2\sum_{q>0}\xi_q\,\delta K_q`$: measure that spurious
  force against $`\langle\tau\rangle`$.
- **Q-SP8.** Why do relative resets lag absolute ones increasingly with
  rate (Q6(c)(iii))? Test the diffusion picture of §8.5 by resolving the
  misalignment field into short and long wavelengths along a row.
- **Q-SP9.** What fixes the dark rate? By Proposition Q1 the observable does
  not. Candidates are the kernel's own rate, a multiple of it, or a constant
  of the ontology; near the Eckart barrier about three times the kernel's
  rate is needed (Q6(c)).
- **Q-SP10.** Gray pairs. A pair formed on contact is generically gray: its
  members bring unrelated clocks. Does it become force-blind on joining,
  by (S′)'s first argument (a $`W`$-null pair), or only on alignment, by its
  second (darkness as a limit)? The first puts force-blind gray pairs into
  the phase reference at the rate pairs form (§6.2).

---

## 12. Sources

- Aharonov, Y. and Bohm, D. *Significance of electromagnetic potentials in
  the quantum theory*, Phys. Rev. **115** (1959) 485–491.
- Anderson, P. W. *A mathematical model for the narrowing of spectral lines
  by exchange or motion*, J. Phys. Soc. Jpn. **9** (1954) 316.
- Bloembergen, N., Purcell, E. M. and Pound, R. V. *Relaxation effects in
  nuclear magnetic resonance absorption*, Phys. Rev. **73** (1948) 679–712.
- Boyd, S., Ghosh, A., Prabhakar, B. and Shah, D. *Randomized gossip
  algorithms*, IEEE Trans. Inf. Theory **52** (2006) 2508–2530.
- Bruns, H. *Das Eikonal*, Abh. Kgl. Sächs. Ges. Wiss., Math.-Phys. Classe
  **21** (1895). The eikonal of geometrical optics: the full phase function
  along curved rays.
- Glauber, R. J. *High-energy collision theory*, in Lectures in Theoretical
  Physics, vol. 1, ed. W. E. Brittin and L. G. Dunham, Interscience (1959)
  315–414. The eikonal: a plane wave carried along straight lines,
  accumulating $`-\int V\,dt/\hbar`$.
- Goodman, J. W. *Introduction to Fourier Optics*, 3rd ed., Roberts & Company
  (2005), chapter 4: Fraunhofer diffraction as a Fourier transform of the
  aperture illumination.
- Molière, G. *Theorie der Streuung schneller geladener Teilchen I.
  Einzelstreuung am abgeschirmten Coulomb-Feld*, Z. Naturforsch. **2a**
  (1947) 133. The high-energy (straight-line) approximation for fast
  charged particles.
- Moyal, J. E. *Quantum mechanics as a statistical theory*, Proc. Cambridge
  Philos. Soc. **45** (1949) 99–124.
- Nedjalkov, M., Kosina, H., Selberherr, S., Ringhofer, C. and Ferry, D. K.
  *Unified particle approach to Wigner–Boltzmann transport in small
  semiconductor devices*, Phys. Rev. B **70** (2004) 115319. The
  signed-particle method with same-cell annihilation.
- Sellier, J. M. *A signed particle formulation of non-relativistic quantum
  mechanics*, J. Comput. Phys. **297** (2015) 254.
- Van Vleck, J. H. and Weisskopf, V. F. *On the shape of collision-broadened
  lines*, Rev. Mod. Phys. **17** (1945) 227–236.
