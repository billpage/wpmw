# The compensated ledger in four dimensions: what a cell costs

> The compensated world form taken out of 1+1 dimensions for the first time: two particles on a line, the ledger run on exact integer counts with bodies carrying their own sub-cell positions, and the one cell size nothing fixes — the bin area `A` in the centre-of-mass directions — scanned. Proposition M1: for a pair potential the residual symbol depends on the ket–bra separation only through `y1 - y2`, so every event moves relative momentum and leaves `P` untouched, and Theorem M2 collapses step 10's leak law: under (S) every event is momentum-neutral in every direction and the field's share goes wholly through the force. Corollary M3 turns step 10's worst case into a free one — the harmonic sector carries no events, so Theorem D and the ring-seam obstruction do not arise — while Proposition M4 records that the residual is not always cheaper than the uncompensated kernel on the same lattice (0.83 for a cosine, 2.8 for the attractive Poeschl–Teller well at the ledger runs' reach, as Theorem C4 of step 14 predicts). **Proposition M5 is the note's centre**: the equilibrium body count is proportional to the number of occupied joint cells, 2.4 to 4.6 bodies per cell, so shrinking `A` from infinity to `h` and `h/4` multiplies it by 4.8 and 13.2, with fidelity error growing as its square root. Proposition M6: the cells filled are those escaped bodies reach, not the state's support, and doubling the window raises the count again. Proposition M7 restricts N4 and N5: under real streaming, partner availability is structural, 0.55 to 0.74 whether the mean occupancy of a requested cell is 1.3 or 10.3. Proposition M8 answers step 10's open item 1 for this layer: at `A` infinite, partners drawn from the wrong centre-of-mass component move the conditional potential energy by 4 within a quarter time unit, so co-location in the centre of mass is not a gauge choice. Proposition M9 finds the same trade-off inside the 1D cell, which the mesh ledger cannot show. Corollary M10: the cost is a product over degrees of freedom — the exponential of step 10's Propositions B2 and B3, returned as body inflation now that the sea is no longer a partner.

*Ladder abstract — see the [full list](README.md#the-ladder).*

**Status.** Analysis note, step 21 of the ladder. Companion demo:
`src/demo_fourd_compensated_ledger.py` (parts E, K, R, W, S, X, B, F; about
five minutes in all). Prompted by the multi-body readiness review of September
2026, which asked whether the compensated model is ready for entanglement or
quantum chemistry and found that its world form had never been run in more
than one degree of freedom. The erratum of
[`../algorithm/multi_body_extension.md`](../algorithm/multi_body_extension.md)
§0 left the question open (item 5); this note measures it. Derivation and
implementation were developed jointly with Claude (Anthropic), September 2026.
Every quantitative claim below is an output of the demo (§9).

## 0. What this note settles, and what it corrects

Step 10, [`fourd_microdynamics.md`](fourd_microdynamics.md), took the
collision-layer ladder into four-dimensional phase space and found two things
that matter here: that excess worlds essentially never meet once the joint
lattice has two or more degrees of freedom (Proposition B2), and that the sea,
present at density $`B`$ in every joint cell, is therefore the only collision
partner there is. The compensated model has since removed the sea from that
role. A catalysed recombination consumes two *bodies*, one at each daughter
row, and binds them into a pair at the parent's own row (Proposition S0/K8 of
the algorithm specification); the sea is a counter, never a partner. Step 10's
dilution therefore comes back, and this note measures what it costs.

**Settles.**

- **What the compensated kernel looks like in 4D** (§2). For a pair potential
  it moves relative momentum only, and every event is momentum-neutral in every
  direction.
- **What postulate (S) changes relative to step 10** (§3): the quadratic
  sectors leave the event channel entirely, the leak law collapses, and the
  event budget depends on the reach rather than on the mode sum.
- **How the ledger's cost scales with dimension** (§4, §8): in proportion to
  the number of occupied joint cells.
- **Whether co-location in the centre of mass matters** (§6). It does; step
  10's transverse freedom is not a gauge freedom in this layer.

**Corrects or restricts.**

- **Theorems N4 and N5 of [`stochastic_ledger.md`](stochastic_ledger.md)**,
  restricted (§5.2). Their independent-occupancy closure was measured with
  uniform hop transport (N-SP4). Under Newtonian streaming, partner
  availability is set by the sign geometry of the populations, not by
  Poisson counting, and the attractor at $`f = 1/2`$ is reached only once the
  population has filled the cells it can reach.
- **Theorem S9 of [`sea_population_equilibrium.md`](sea_population_equilibrium.md)**,
  sharpened. The attractor holds, but its approach runs through §4's filling,
  so the body count at which it is reached depends on the window (§5.1).
- **Proposition C4 of step 10** does not transfer (§6). It was stated for a
  transversally uniform sea; the compensated ledger's partners are
  state-dependent bodies, and the freedom it called invisible corrupts a
  correlated state within a quarter time unit.
- **The readiness review's scaling argument**, corrected in mechanism. It
  predicted a per-cell occupancy fixed at $`\lambda^*`$ inside the state's
  support. The occupancy is not Poisson (§5.2), and the cells filled are the
  window's, not the support's (§5.1). The conclusion — cost exponential in the
  number of degrees of freedom — survives (§8).

What this note does not touch: §8 of
[`compensated_ontology.md`](compensated_ontology.md) states that the
compensated symbol carries over to $`N`$ bodies unchanged, and that remains
true of the symbol (Proposition M1 is an instance). What is measured here is
the world form's cost, which that section does not address.

---

## 1. The problem, and why it is bound

### 1.1 Coordinates

Two particles of unit mass on a line, $`\hbar = 1`$. The centre of mass and
relative coordinates

```math
X = \tfrac12(x_1 + x_2),\quad P = p_1 + p_2,\qquad
r = x_1 - x_2,\quad p_r = \tfrac12(p_1 - p_2)
```

are canonical with unit Jacobian in both positions and momenta, so the joint
Wigner function transforms as a scalar: a product state in $`(X, P)`$ and
$`(r, p_r)`$ has $`W^{(2)} = W_X(X,P)\thinspace W_r(r,p_r)`$. The masses are
$`M = 2`$ and $`\mu = 1/2`$. The Hamiltonian is

```math
H = \frac{P^2}{2M} + \tfrac12 M\Omega^2 X^2 + \frac{p_r^2}{2\mu} - V_0\,\mathrm{sech}^2(r/a),
\qquad V_0 = 6,\ a = 1,\ \Omega = 1 .
```

The pair term is the attractive Pöschl–Teller well. With
$`\lambda(\lambda+1) = 2\mu V_0 a^2/\hbar^2 = 6`$, $`\lambda = 2`$, and the
bound states are exact: $`\psi_0 \propto \mathrm{sech}^2 r`$ with
$`E_0 = -4`$ and $`\psi_1 \propto \mathrm{sech}\thinspace r\tanh r`$ with
$`E_1 = -1`$. The centre-of-mass trap is quadratic, so its Moyal bracket
truncates and it carries no residual at all (G2); its thermal states, with
$`\sigma_X^2 = (2\bar n + 1)/2M\Omega`$ and
$`\sigma_P^2 = (2\bar n + 1)M\Omega/2`$, are stationary under classical
rotation. Every state used below is therefore stationary or known exactly at
every time.

### 1.2 Why not the Eckart collision

The review proposed the two-particle Eckart collision,
$`U = V_0\thinspace\mathrm{sech}^2((x_1 - x_2)/a)`$ with $`V_0 = 1`$, $`a = 2`$,
which is exactly the relative-coordinate problem of step 15 lifted to
particle coordinates. Part E confirms it is an exact 4D reference: a 2D
split-operator run in $`(x_1, x_2)`$ from a product of Gaussians gives
transmission $`0.544134`$, against $`0.544148`$ from the 1D relative run and
$`0.544156`$ from the packet-averaged closed form, with the relative density
agreeing to an $`L^1`$ distance of $`2.7\times10^{-5}`$. The collision also
entangles: the reduced purity of particle 1 falls from 1 to $`0.4416`$, an
entropy of $`0.976`$ nats, about that of a Bell pair.

It is kept for the ledger as open item M-SP3, because a scattering packet
leaves the interaction region before the ledger has settled, and the
questions of this note are about the settled ledger. A bound stationary state
keeps the events running for as long as the run lasts.

### 1.3 The reference

The ledger is compared with the mesh QLE of the same compensated symbol on
the same $`(r, p_r)`$ mesh, times the exact $`W_X`$. That isolates the
ledger's error from the discretisation's. The discretisation's own error is
recorded for completeness (Part R): over $`t = 4`$ the mesh QLE moves
$`\psi_0`$'s Wigner function by a relative $`L^2`$ of $`2.0\times10^{-2}`$ and
$`\psi_1`$'s by $`0.25`$. Mesh: $`128\times64`$, $`r \in [-10, 10)`$,
$`\Delta r = 0.156`$, $`\Delta p = 0.25`$, reach $`y_{\max} = 2\pi a`$ — the
step-16 ledger parameters, with $`\Delta r`$ halved for the narrower well. A
mesh cell has area $`h/161`$.

---

## 2. The world form in four dimensions

### 2.1 Proposition M1: the pair kernel moves relative momentum only

The joint symbol is
$`M(\mathbf x, \mathbf s) = (i/\hbar)\bigl[V(\mathbf x + \mathbf y) - V(\mathbf x - \mathbf y)\bigr]`$
with $`\mathbf y = \hbar\mathbf s/2`$, compensated as in the specification.

**Proposition M1 (pair support).** *For $`V(x_1, x_2) = U(x_1 - x_2)`$ the
residual symbol depends on $`\mathbf s`$ only through $`s_1 - s_2`$. Every
event therefore transfers $`(\xi, -\xi)`$ to $`(p_1, p_2)`$: it moves $`p_r`$
by $`\xi`$ and leaves $`P`$, $`X`$ and $`r`$ unchanged, so an emissive
daughter inherits its parent's centre-of-mass coordinates exactly.*

*Proof.* With $`\eta = y_1 - y_2`$,
$`V(\mathbf x \pm \mathbf y) = U(r \pm \eta)`$, and the compensating term
$`\mathbf s\cdot\nabla V`$ equals $`(s_1 - s_2)\thinspace U'(r)`$. A function
$`f(s_1 - s_2)`$ has transform
$`\iint f(s_1 - s_2)\thinspace e^{-i(s_1\xi_1 + s_2\xi_2)}\thinspace ds_1\thinspace ds_2 = \tilde f(\xi_1)\thinspace\delta(\xi_1 + \xi_2)`$,
by the substitution $`u = s_1 - s_2`$, $`v = s_2`$. $`\square`$

This is step 10's pair mode $`\vec q = (q, -q)`$, and its fibration
(Proposition A2) and its centre-of-mass-degenerate sea (§4.4) carry over
unchanged. The residual in $`(r, p_r)`$ is the 1D compensated kernel at mass
$`\mu`$, with the lattice in $`p_r`$ fixed by the reach in $`\eta`$ (Theorem C4
of step 14). Nothing fixes a lattice in $`P`$: the kernel is a delta at zero
transfer there.

### 2.2 Theorem M2: every event is momentum-neutral

**Theorem M2.** *Under postulate (S), every compensated event leaves every
linear momentum functional $`u\cdot\mathbf P`$ unchanged, for every mode
geometry. All momentum exchange with the field is carried by the classical
force.*

*Proof.* An event deposits $`+1`$ at $`\mathbf p + \boldsymbol\xi`$ and $`-1`$ at
$`\mathbf p - \boldsymbol\xi`$ and binds or ionises a pair at $`\mathbf p`$
itself (S0/K8). Counting the created or removed bodies with their momenta,
the displacement is
$`(\mathbf p + \boldsymbol\xi - \mathbf p) + (\mathbf p - \boldsymbol\xi - \mathbf p) = 0`$.
The hops, which would move one body by $`2\boldsymbol\xi`$, are excluded by
(S) (specification §5.3). $`\square`$

Step 10's leak law, Theorem A3, had Focus and Defocus neutral and the hops
leaking $`2(u\cdot\vec q)\Delta p`$. The compensated demography is Focus and
Defocus only, so the law collapses to its neutral half. For a pair potential
$`P`$ is conserved twice over: by every event, and by the forces, which are
$`\pm U'(r)`$.

### 2.3 Definition M0: the matching cell

**Definition M0 (matching cell).** *The cell in which an absorptive event
looks for its two partners: an $`r`$-bin, a $`p_r`$ row, and a bin of area
$`A`$ in $`(X, P)`$ with aspect $`\sigma_X/\sigma_P`$. $`A = \infty`$ is one
bin, the 1D ledger with the centre of mass carried as a passive label.*

The $`p_r`$ row is forced by the reach. The $`r`$-bin is the mesh spacing,
unforced (N-SP2). The $`(X, P)`$ bin is unforced and is the scanned parameter.

### 2.4 The implementation

Bodies are signed points $`(X, P, r, p_r)`$, each with its own sub-cell
position. Postulate (S): exact rotation in the centre-of-mass trap, velocity
Verlet in $`r`$ under the kernel's own first moment. Postulate (D): each body
parents Poisson events at rate $`\Gamma(r) = \sum_{q\ge1}\lvert K_q(r)\rvert`$
per unit time, channel drawn in proportion to $`\lvert K_q\rvert`$, processed
sequentially in random order with live partner counts (specification §5.5).
The box is periodic in $`r`$ and in $`p_r`$, as the mesh reference is. The
initial ensemble samples $`\lvert W\rvert`$ with $`\nu = 4000`$ signed bodies
per unit of $`\lVert W\rVert_1`$ and is padded with co-located $`\pm`$ pairs to
$`\rho = 3`$, which leaves $`E`$ unchanged in every cell under any binning.
$`\Delta t = 0.01`$. Counts are integers, and Theorem S7's
$`\Delta N = 2(1 - 2f)\thinspace n_{\rm ev}`$ holds with residual exactly 0 in
every run.

---

## 3. Against step 10

| | step 10, collision layer | this note, compensated layer |
|---|---|---|
| momentum lattice | $`\Delta p = \pi\hbar/L`$, fixed by the box | $`\Delta p = \pi\hbar/2y_{\max}`$, fixed by the reach |
| channel index | a Fourier mode of $`V`$ | a momentum transfer $`\xi_q = q\Delta p`$ |
| first moment of the kernel | the classical force | zero |
| realisations | Focus, Defocus, Right-Hop, Left-Hop | emissive and absorptive only |
| leak law | A3: hops leak $`2(u\cdot\vec q)\Delta p`$ | M2: no event leaks |
| pair-mode geometry | $`(q, -q)`$ chains (A2) | $`(\xi, -\xi)`$ (M1) |
| harmonic potential | worst case (Theorem D) | no events (M3) |
| partners at $`N > 1`$ | the sea, density $`B`$ in every joint cell (B2) | bodies at both daughters |
| transverse freedom | invisible at uniform sea (C4) | visible and destructive (M8) |

### 3.1 Corollary M3: the quadratic sectors are free

**Corollary M3.** *A quadratic sector — the centre-of-mass trap here, any
harmonic coupling in general — carries no events. Step 10's Theorem D, the
ring-seam obstruction of its §5.4 and the Euler amplification of its §5.6 do
not arise.*

It follows from G2 of step 17 and I4 of step 12. Theorem D priced the
harmonic potential at an injected momentum variance of $`2m\omega^2\hbar`$
per retained mode; here the centre of mass is carried by exact rotation at no
cost, and its thermal states stay stationary to the integrator's precision.
Step 10's open item 4, the classical-drift/quantum-remainder split, is what
the compensated model adopted as postulate (S); the tension with the
four-action ontology that item recorded is resolved by changing the
ontology, not by exempting the implementation.

### 3.2 Proposition M4: the event budget depends on the reach

Part K builds the uncompensated kernel (the windowed full symbol) and the
residual on the same mesh.

**Proposition M4 (measured).** *On a common lattice with
$`\Delta p = 0.25`$ and $`y_{\max} = 2\pi`$, the ratio of maximal residual to
maximal uncompensated event rate is $`0`$ for the harmonic potential,
$`0.834`$ for $`1.5\cos(\pi r/5)`$, and $`2.825`$ for the Pöschl–Teller well.
The uncompensated kernel's first moment reproduces the force to within the
window's effect ($`+3.4753`$ against $`-V' = +3.4699`$ for the well); the
residual's is below $`2\times10^{-14}`$ in all three.*

This is Theorem C4 of step 14 again: compensation gains at short reach and
loses beyond about a quarter wavelength of the dominant mode. The ledger runs
use $`y_{\max} = 2\pi a`$, where the well is in the losing regime; a shorter
reach would cut the event rate but meets the over-determination of CLA7. At
maximal reach on a ring the uncompensated kernel is step 10's two-atom stencil
(checked in step 14), and there compensation only reorganises.

### 3.3 What carries over and what does not

Theorem A1's exactness family, the four-action form, and Proposition B1 on
the crystal shift have no counterpart: the compensated observable carries no
background. Proposition B2's arithmetic carries over exactly, and B3's
conclusion — the sign problem traded for a noise problem, neither removed —
reappears in §4 as body inflation instead of sea shot noise. Theorem C1's
transverse family is §6's subject.

---

## 4. Proposition M5: the population fills cells

The scan (Part S), with $`\bar n = 0`$ unless stated. "Saturated" means the
windowed absorptive fraction has returned to $`f \approx 1/2`$ and the body
count has levelled.

| $`A`$ | $`t`$ | $`N`$ | occupied cells | $`N`$/cell | $`f`$ | error |
|---|---|---|---|---|---|---|
| $`\infty`$ | 4.0 | 29 482 | 7 402 | 3.98 | 0.492 | 0.942 |
| $`h`$ | 5.0 | 142 530 | 44 023 | 3.24 | 0.457 | 1.819 |
| $`h/4`$ | 4.5 | 389 124 | 111 256 | 3.50 | 0.491 | 2.926 |
| $`h/16`$ | 2.0 (growing) | 442 948 | 184 991 | 2.39 | 0.408 | 2.785 |
| $`h/4`$, $`\bar n = 3`$ | 2.0 (growing) | 680 514 | 285 668 | 2.38 | 0.404 | 3.704 |

The error is the relative $`L^2`$ distance of the $`(r, p_r)`$ marginal from
the reference; its sampling floor at $`t = 0`$ is $`0.194`$.

**Proposition M5 (measured).** *The body count is proportional to the number
of occupied joint cells, at $`2.4`$ to $`4.6`$ bodies per cell for every $`A`$
from $`t = 1`$ on, while it grows and after it saturates. Shrinking $`A`$ from $`\infty`$ to $`h`$
and $`h/4`$ multiplies the saturated count by $`4.8`$ and $`13.2`$, while the
occupied joint cells grow by $`5.9`$ and $`15.0`$,
and the fidelity error grows as roughly the square root: $`\times1.93`$ and
$`\times3.11`$ measured against $`\times2.20`$ and $`\times3.63`$ from
$`\sqrt N`$.*

Panel (b) of the figure is the proposition: every run, at every $`A`$ and
both $`\bar n`$, lies on one line. While the count grows, $`f`$ sits at
$`0.40`$ to $`0.42`$, below the half the attractor needs, and the population
grows at $`2(1 - 2f)`$ times the event rate until the reachable cells are
filled; then $`f`$ returns to $`0.49`$. A thermal centre of mass
($`\bar n = 3`$) occupies more bins and grows faster, $`680\thinspace514`$
against $`152\thinspace712`$ at $`t = 2`$.

![The compensated ledger in four dimensions](https://raw.githubusercontent.com/billpage/wpmw/output/figures/fourd_compensated_ledger.png)

*(a) Body count against time for each $`A`$, and with the $`r`$-window doubled
at $`A = \infty`$. (b) Body count against occupied joint cells: one line for
every run. (c) Conditional $`\langle V(r)\rangle`$ on each side of the
centre-of-mass bisector for §6's two-component state; thick lines are the
reference. (d) §7's $`\langle V\rangle`$ bias at $`t = 1`$ against the body
count, as the matching cell is refined; the grey band is classical streaming
alone.*

---

## 5. Where the cells are

### 5.1 Proposition M6: the cells are the window's

**Proposition M6 (measured).** *The occupied cells are those the escaped
bodies reach, not the state's support. At $`A = \infty`$, 34 per cent of the
bodies sit where $`\Gamma(r)`$ is below one per cent of its maximum at
$`t = 4`$; doubling the $`r`$-window to $`[-20, 20)`$ raises that to 63 per
cent and the body count to $`50\thinspace732`$ at $`t = 7`$ against
$`29\thinspace482`$, still rising slowly, with $`13\thinspace811`$ of
$`16\thinspace384`$ cells occupied.*

Emissive daughters are displaced by up to $`\xi = 7.75`$, enough to leave the
well, and outside the event region nothing can remove them: every removal is
an absorptive event, and events happen only where $`\Gamma \ne 0`$. The
periodic box returns them to the well, which is how the ledger closes here at
all. On an open line the escaped flux does not return, and the bound problem
is expected to have no equilibrium (M-SP1). This is CLA8's observation that
the sea scales with the window, made now for the unpaired population. The
escaped bodies arrive in both species in roughly equal numbers, so a
contact-recombination sink would remove them where catalysed recombination
cannot. [`../../ORIENTATION.md`](../../ORIENTATION.md) §8 records that such
a sink is not part of the model, because Theorem S9 holds the ledger closed
without one; M6 says that holds here only because the box returns what
escapes. Whether that justifies revisiting the exclusion is M-SP2.

### 5.2 Proposition M7: availability is structural

Theorem N4 prices the absorptive fraction as the probability that both legs
find a partner, and under independent Poisson occupancy that is
$`(1 - e^{-\lambda})^2`$.

**Proposition M7 (measured).** *Under Newtonian streaming the per-leg
availability stays between $`0.55`$ and $`0.74`$ throughout the scan, while
the mean occupancy of the requested cells ranges from $`1.3`$ to $`10.3`$, over
which Poisson occupancy would give $`0.73`$ to $`0.99997`$. Availability barely
responds to occupancy: it is set by the sign geometry of the two populations,
not by counting statistics.*

A requested cell either holds the needed species in quantity or holds none of
it: where $`W`$ is strongly positive a negaton is structurally absent,
whatever the total count. This is the conjunction cost of §5.4 of the
specification (legs available 67 to 76 per cent of the time, anti-correlated)
seen from the occupancy side. N4 and N5 were measured with uniform hop
transport, the setting N-SP4 flagged as not yet transport; with streaming the
closure does not describe the ledger, and the five numbers of N5 should not
be used to size a population. The attractor itself survives: $`f`$ returns to
$`0.49`$ in every saturated run.

---

## 6. Proposition M8: co-location in the centre of mass is not a gauge choice

Step 10's Theorem C1 found that in two or more dimensions the exchange vertex
has a transverse family, and in the pair case the free parameter is a
mismatch of $`P`$ between the incoming world and its partner. Proposition C4
called the family invisible to the QLE at a transversally uniform sea, and
open item 1 left open whether to postulate it away.

The compensated ledger has no uniform sea to average against. Part B takes a
state in which the centre of mass carries information about the relative
motion: an equal mixture of two product states, coherent centre-of-mass
states at $`X = \pm2`$ (eight $`\sigma_X`$ apart) carrying $`\psi_0`$ and
$`\psi_1`$ respectively. Each component's relative state is stationary, and
the components rotate rigidly in the trap, so the reference for each side of
the bisector is exact at every time.

**Proposition M8 (measured).** *At $`A = \infty`$, where partners may come
from either component, the conditional $`\langle V(r)\rangle`$ moves from
$`-4.770`$ to $`-1.047`$ on the $`\psi_0`$ side and from $`-2.313`$ to
$`-5.973`$ on the $`\psi_1`$ side by $`t = 0.25`$, against references of
$`-4.793`$ and $`-2.392`$. Over five seeds at $`t = 1.5`$ the offsets are
$`+4.049 \pm 0.054`$ and $`-3.697 \pm 0.090`$ at $`A = \infty`$ and
$`+0.364 \pm 0.120`$ and $`+0.044 \pm 0.225`$ at $`A = h/4`$.*

The two sides do not relax toward their average, $`-3.6`$; they pass each
other. That is what would follow if each side received the other's residual
corrections, and it is a stronger effect than mixing. At $`A = h`$ the bins
separate the components and the conditional energies stay within noise until
the body count, and with it the noise, grows (panel c).

So for this layer step 10's open item 1 is answered: partners must share the
parent's centre-of-mass cell, to a resolution finer than any correlation
between the centre of mass and the relative motion. For an entangled state
that correlation is the state.

---

## 7. Proposition M9: sub-cell positions

Bodies here carry positions and momenta below cell resolution, which is the
premise of open item Y-SP1, and absorption removes whichever partner sits in
the cell, not one at the parent's own point. The mesh ledger cannot show the
consequence, because there a cell is a point.

**Proposition M9 (measured).** *At $`A = \infty`$ and $`\nu = 16\thinspace000`$,
the particle ledger's $`\langle V\rangle`$ exceeds the reference by
$`+0.133 \pm 0.013`$ at $`t = 1`$ (six seeds), against $`+0.011 \pm 0.003`$ for
classical streaming alone. Refining only the matching cell, with the kernel's
mesh unchanged, reduces the bias to $`+0.066 \pm 0.020`$ (four sub-bins in
$`r`$), $`+0.043 \pm 0.009`$ (four in $`p_r`$) and $`+0.034 \pm 0.013`$ (both),
at $`1.7`$, $`1.7`$ and $`4.2`$ times the body count.*

It is §6's trade-off one level down, and it makes N-SP2 concrete: the cell is
not a unit of the mesh but a parameter of the demography, with an error on one
side and a population on the other. Refining in $`p_r`$ helps most, which is
consistent with the mechanism: the classical force moves momenta continuously,
so under (S) a row is a bin rather than a lattice site, which it was in the
collision layer.

---

## 8. Corollary M10: what the cost scales with

**Corollary M10.** *The world form's population, and with it its event count
and its noise, scales with the number of occupied joint cells. With $`K_j`$
bins needed along degree of freedom $`j`$, that is
$`\prod_j K_j`$, and M8 bounds each bin from above by the scale of the state's
correlations, so the product cannot be kept small by coarse binning.*

For one added centre-of-mass pair the measured factor is $`5.9`$ at $`A = h`$
and $`15.0`$ at $`A = h/4`$; §5.1 adds a factor set by the window. This is the
exponential of step 10's Propositions B2 and B3, which the sea absorbed there
and which returns here as body inflation. It is the concrete form of what the
readiness review expected, with the mechanism corrected: cell filling, not a
Poisson occupancy fixed inside the support.

---

## 9. Numerical verification

All outputs are from `src/demo_fourd_compensated_ledger.py`, container run
September 2026, seeds fixed in the code.

- **Part E.** Eckart pair collision: transmissions and entanglement of §1.2.
- **Part K.** Kernel comparison of §3.2.
- **Part R.** Mesh-QLE drift of the two eigenstates, §1.3.
- **Part W.** Proposition M9: six seeds per row, classical control included.
- **Part S.** The scan of §4, S7 residual 0 in every run.
- **Part X.** Proposition M6, the doubled window.
- **Part B.** Proposition M8: two time series and the five-seed table.
- **Part F.** The figure, from the saved outputs of S, X, B and W.

Parts run separately: `python3 src/demo_fourd_compensated_ledger.py S` runs
the scan alone.

---

## 10. Open items

- **M-SP1.** The open line. Does the bound problem have any ledger equilibrium
  when escaped bodies cannot return? A series of growing windows, or an
  absorbing boundary, would settle it.
- **M-SP2.** Contact recombination. §5.1's escaped bodies are the population
  a contact sink would remove. Whether that justifies revisiting the exclusion
  recorded in ORIENTATION §8 is an ontology decision.
- **M-SP3.** Run the Eckart pair collision of §1.2 through the ledger and
  compare its transmission with the reference Part E has verified.
- **M-SP4.** Repeat the scan at a reach where Proposition M4's ratio is below
  one, and measure what CLA7 costs there.
- **M-SP5.** A non-separable problem, adding an anharmonic external potential,
  with the 2D split-operator reference; and bins in particle coordinates
  rather than aligned with the centre of mass.
- **M-SP6.** Is there a matching rule that removes M9's bias without
  multiplying the population? Connects Y-SP1 and N-SP2.
- **M-SP7.** Fidelity per unit cost: how the error at saturation depends on
  $`\nu`$ once the population is set by the cells rather than by the sample.
- **M-SP8.** Event-ordering sensitivity (S-SP4) of the body counts in §4.
- **Step 10, open item 1** is answered for this layer by M8. **N-SP2** is
  sharpened by M9, **N-SP4** answered in part by M7, and **Y-SP1** has a
  measured consequence in M9.

---

## 11. Sources

- Pöschl, G.; Teller, E. *Bemerkungen zur Quantenmechanik des anharmonischen
  Oszillators*, Z. Phys. **83** (1933) 143–151. The exactly solvable well of
  §1.
- Eckart, C. *The penetration of a potential barrier by electrons*, Phys. Rev.
  **35** (1930) 1303–1309. The transmission of §1.2.
- Sellier, J. M.; Dimov, I. *The many-body Wigner Monte Carlo method for
  time-dependent ab-initio quantum simulations*, J. Comput. Phys. **273**
  (2014) 589–597. The published two-particle signed-particle method, whose
  annihilation step on a phase-space grid faces the same cell question.
