# The Eckart barrier: tunnelling as an imbalance in separatrix traffic

Companion script: [`../../src/demo_eckart_barrier_compensated.py`](../../src/demo_eckart_barrier_compensated.py)

## 0. What this note inherits, and why this potential

### 0.1 Inherited

- [`compensated_liouville_splitting.md`](compensated_liouville_splitting.md) §1–2 for
  the symbol $`M(x,s) = (i/\hbar)[V(x+y) - V(x-y)]`$ with $`y = \hbar s/2`$, the
  split $`M = M_{\rm cl} + M_{\rm res}`$ with $`M_{\rm cl} = iV'(x)s`$, Theorem C1
  (the two factors commute, so the split carries no Trotter error inside the
  potential substep), Theorem C2 (the residual is the odd part of the cubic
  Taylor remainder), Theorem C3 (zeroth and first moments of $`K_{\rm res}`$
  vanish), Theorem C4 (the reach condition, per mode $`|\sin u - u|/|u|`$ with
  $`u = k\thinspace y_{\max}`$), Theorem C6 (Coulomb converges iff the reach
  misses the nucleus) and Theorem C7 (the quiet region).
- [`open_position_space.md`](open_position_space.md) §3 for definitions (H) and
  (R): the coherence horizon $`L_c`$ and the reach $`y_{\max} = L_c/2`$.
- [`reach_energy_coupling.md`](reach_energy_coupling.md) Theorem E1
  ($`\Delta p = \pi\hbar/(2y_{\max})`$, the reach is a period) and Theorem E7.

### 0.2 Why sech²

Two results already in the ladder close off the obvious test problems.

**Theorem E7** shows that for a polynomial $`V`$ the Moyal series terminates,
so the collision term is a finite-order differential operator in $`p`$ with
unbounded total variation, and there is no jump measure on the open line at
all until a reach is imposed.

**Theorem I4** of [`interworld_coupling.md`](interworld_coupling.md) shows that
a coupling linear in the leg separation generates classical Liouville flow
exactly, so the harmonic *and* the inverted harmonic have no jump channel
whatever. Its §4 already draws the consequence for the project's current test
problem: the inverted pair barrier of
[`../supplement/inverted_pair_barrier.md`](../supplement/inverted_pair_barrier.md)
is quadratic in the relative coordinate, so it tests transport and not the
four rules.

Between them the project has had no open-line test problem that exercises the
hop channel against a closed form. The Eckart (Pöschl–Teller) pair barrier

```math
V(r) \;=\; V_0\thinspace\mathrm{sech}^2(r/a)
```

supplies one. It is bounded, asymptotically free on both sides, has
$`V''' \not\equiv 0`$, and has an exact transmission coefficient. Write

```math
\beta \;=\; \frac{\sqrt{2\mu V_0}\thinspace a}{\hbar},
\qquad
p_b \;=\; \sqrt{2\mu V_0} \;=\; \frac{\hbar\beta}{a},
```

for the semiclassical parameter and the barrier momentum scale, with
$`\mu = m/2`$ the reduced mass of the pair. The exact transmission is

```math
T(E) \;=\;
\frac{\sinh^2(\pi k a)}
     {\sinh^2(\pi k a) + \cosh^2\!\bigl(\tfrac{\pi}{2}\sqrt{4\beta^2 - 1}\bigr)},
\qquad k = \sqrt{2\mu E}/\hbar .
```

### 0.3 What this note claims

Seven theorems, K1–K7. The one that matters for the ontology is K4 together with
K5: **the deterministic step conserves the classical outcome exactly, so the
whole of tunnelling is delivered by the residual channel, and it arrives as a
small imbalance between two large opposed flows of positon–negaton pairs
across the classical separatrix.**

### 0.4 What this note corrects

- §5.1 of [`compensated_liouville_splitting.md`](compensated_liouville_splitting.md)
  reads the reach as *translating* the interaction profile outward by
  $`y_{\max}`$. Theorem K2 below shows that is a Gaussian fact. For an
  exponentially-tailed potential the reach rescales the profile and does not
  translate it, because the residual decays at exactly $`V`$'s own rate.
- Theorem C6 is generalised by K1: the reach ceiling is not special to
  Coulomb's real-axis pole but is the distance to the nearest complex
  singularity of $`V`$, whatever and wherever it is.

---

## 1. The split for sech²

With $`c(u) = \mathrm{sech}^2 u`$ and $`t(u) = \tanh u`$ at $`u = x/a`$,

```math
V' = -\frac{2V_0}{a}\thinspace c\thinspace t,
\qquad
V''' = \frac{8V_0}{a^3}\thinspace c\thinspace t\thinspace(2c - t^2),
```

so $`V'''`$ is odd, vanishes at the barrier top and decays as
$`e^{-2|x|/a}`$. The residual symbol is
$`M_{\rm res} = (i/\hbar)[V(x+y) - V(x-y) - 2yV'(x)]`$.

Theorem C2's bound is not tight but is the right order at short reach
(part A of the script, $`x = 0.5`$, $`a = V_0 = 1`$):

| $`y_{\max}`$ | $`\max\lvert M\rvert`$ | $`\max\lvert M_{\rm res}\rvert`$ | C2 bound | ratio |
|---|---|---|---|---|
| 0.10 | 0.1441 | $`1.31\times10^{-3}`$ | $`2.0\times10^{-3}`$ | 0.664 |
| 0.25 | 0.3434 | $`2.00\times10^{-2}`$ | $`3.1\times10^{-2}`$ | 0.648 |
| 0.50 | 0.5800 | $`1.47\times10^{-1}`$ | $`2.5\times10^{-1}`$ | 0.595 |
| 1.00 | 0.6597 | $`8.48\times10^{-1}`$ | 1.975 | 0.429 |

Theorem C3 holds as expected. At $`y_{\max} = 0.25`$, $`x = 0.5`$ the full
kernel's first moment is 0.727069 against $`\lvert V'\rvert = 0.726862`$,
while the residual's first moment is $`1.2\times10^{-5}`$ and its zeroth
moment $`1.5\times10^{-18}`$ — the force is entirely in the deterministic
factor, to nearly five orders.

---

## 2. Theorem K1: the reach ceiling is the analyticity strip

**Theorem K1.** The Moyal series — the expansion of $`M`$ in odd powers of the
half separation $`y`$ — converges at midpoint $`x`$ if and only if

```math
y_{\max} \;<\; R(x) \;=\; \mathrm{dist}\bigl(0,\;
  \{\thinspace y : x \pm y \text{ is a singularity of } V\thinspace\}\bigr).
```

For $`V_0\thinspace\mathrm{sech}^2(z/a)`$ the singularities are double poles at
$`z = i\pi a(n + \tfrac12)`$, so

```math
R(x) \;=\; \sqrt{x^2 + (\pi a/2)^2},
\qquad
\inf_x R \;=\; \frac{\pi a}{2} .
```

*Proof.* $`y \mapsto V(x+y) - V(x-y) - 2yV'(x)`$ is holomorphic on the largest
disc about $`y = 0`$ free of singularities of either arm, and its Taylor
series is the Moyal series by Theorem I4. $`\square`$

Measured by fitting $`|c_n| \sim Cn^\alpha R^{-n}`$ to the Cauchy coefficients
(part B):

| $`a`$ | $`x`$ | predicted $`R`$ | measured | ratio |
|---|---|---|---|---|
| 1 | 0.30 | 1.5992 | 1.5953 | 0.998 |
| 1 | 0.70 | 1.7197 | 1.7189 | 1.000 |
| 1 | 1.00 | 1.8621 | 1.8576 | 0.998 |
| 2 | 0.30 | 3.1559 | 3.1405 | 0.995 |
| 2 | 0.70 | 3.2186 | 3.2221 | 1.001 |
| 2 | 1.00 | 3.2969 | 3.2913 | 0.998 |

Three points of numerical method are worth recording, because a naive
estimator gets this badly wrong. The evaluation circle must sit close to the
singularity, or the coefficients hit roundoff before the asymptotic regime;
the retention threshold must be applied to the raw transform coefficients, not
to $`|c_n|`$, whose noise floor grows with $`n`$; and the $`n^\alpha`$ factor
of an algebraic singularity must be fitted rather than ignored, since omitting
it biases a two-parameter fit low by tens of per cent.

The direct demonstration is more convincing than the fit. At $`a = 1`$,
$`x = 0.3`$, so $`R = 1.599`$, the Taylor partial sums against the exact
residual:

| $`y`$ | exact | $`S_{10}`$ | $`S_{20}`$ | $`S_{30}`$ | |
|---|---|---|---|---|---|
| 0.50 | $`1.312\times10^{-1}`$ | $`1.312\times10^{-1}`$ | $`1.312\times10^{-1}`$ | $`1.312\times10^{-1}`$ | inside |
| 1.00 | $`6.891\times10^{-1}`$ | $`6.351\times10^{-1}`$ | $`6.879\times10^{-1}`$ | $`6.891\times10^{-1}`$ | inside |
| 1.70 | 1.667 | $`-1.09\times10^{1}`$ | $`-4.56\times10^{1}`$ | $`4.15\times10^{1}`$ | outside |
| 2.20 | 2.287 | $`-1.57\times10^{2}`$ | $`-7.46\times10^{3}`$ | $`1.01\times10^{5}`$ | outside |

**Corollary K1.1 (the lattice cannot resolve the barrier).** $`y_{\max} < \pi a/2`$
with Theorem E1 gives

```math
\Delta p \;=\; \frac{\pi\hbar}{2y_{\max}} \;>\; \frac{\hbar}{a},
\qquad\text{so}\qquad
\frac{p_b}{\Delta p} \;<\; \beta .
```

Fewer than $`\beta`$ rungs span the barrier's own momentum scale. Worse,
resolving a packet needs $`\Delta p \lesssim \sigma_p = \hbar/(2\sigma_r)`$,
i.e. $`y_{\max} > \pi\sigma_r`$, which with the ceiling forces
$`\sigma_r < a/2`$: **the reach ceiling requires the packet to be narrower in
position than half the barrier width.** This is why the transmission
calculation of §5 is run at effectively unbounded reach; the obstruction is a
result, not an evasion.

**Corollary K1.2 (soft-core Coulomb).** For $`-Z/\sqrt{z^2 + \epsilon^2}`$ the
branch points sit at $`z = \pm i\epsilon`$, so $`R(x) = \sqrt{x^2 + \epsilon^2}`$
and at the origin the ceiling *is* the softening length. Measured at
$`x = 0.4`$: $`\epsilon = 0.3`$, predicted 0.5000, measured 0.4984;
$`\epsilon = 0.6`$, 0.7211 against 0.7177; $`\epsilon = 1.0`$, 1.0770 against
1.0771. This is a better motivation for a soft core than regularisation: the
softening length is the largest coherence reach the potential will support.

Theorem C6 is the special case in which the singularity lies on the real axis,
so the ceiling collapses to $`|x|`$ and vanishes at the nucleus. A potential
analytic in a strip of half-width $`d`$ has a uniform positive ceiling
$`y_{\max} < d`$; an entire potential such as a Gaussian barrier has none.

---

## 3. Theorem K2: the far field, and an erratum

**Theorem K2.** For $`x \gg a`$ and $`y_{\max} \ll x`$,

```math
\max_{|y| \le y_{\max}} \bigl|M_{\rm res}(x, y)\bigr|
\;\longrightarrow\;
\frac{8V_0}{\hbar}\thinspace e^{-2x/a}
\Bigl[\thinspace\sinh\!\bigl(2y_{\max}/a\bigr) - 2y_{\max}/a\thinspace\Bigr].
```

*Proof.* $`V \sim 4V_0e^{-2x/a}`$, so
$`V(x+y) - V(x-y) \sim -8V_0e^{-2x/a}\sinh(2y/a)`$ and
$`2yV'(x) \sim -(16V_0/a)\thinspace y\thinspace e^{-2x/a}`$; the bracket is the
difference, maximised at $`|y| = y_{\max}`$. $`\square`$

Verified to five figures (part C): at $`x = 8`$ the ratio of measured to
predicted is 1.00000 for $`y_{\max} = 0.5`$, 1.0, and 2.0.

The bracket is Lemma C0's $`\sin u - u`$ continued to $`u = 2iy/a`$: an
exponential tail is a Fourier mode at imaginary wavenumber, and the whole
structure of the reach condition survives the continuation intact.

**Erratum for `compensated_liouville_splitting.md` §5.1.** That section states
that for a barrier the interaction region is the barrier profile *translated
outward by exactly the reach*, inferred from a Gaussian test barrier.
Theorem K2 shows this does not generalise. The residual decays at exactly
$`V`$'s own rate $`2/a`$, so translation and rescaling are indistinguishable
for an exponential tail; what the reach does is multiply the profile by
$`\sinh(2y_{\max}/a) - 2y_{\max}/a`$. Measured at $`y_{\max} = 1`$, the ratio
$`\max|M_{\rm res}|/|V'''(x)|`$ converges to 0.20336 as $`x`$ grows, a
constant, not a shifted profile. A Gaussian barrier can distinguish the two
because its decay rate accelerates; an exponential one cannot. Theorem C7
itself — the quiet region, no events where $`V'''`$ vanishes on the reach — is
unaffected.

---

## 4. Theorem K3: on the open line the split never loses

Weighting Theorem C4's per-mode ratio by the potential's own spectrum,
$`\tilde V(k) = \pi V_0 a^2 k / \sinh(\pi k a/2)`$, gives a budget ratio that
depends only on $`y_{\max}/a`$ (part D):

| $`y_{\max}/a`$ | 0.1 | 0.25 | 0.5 | 1.0 | $`\pi/2`$ | 2 | 3 | 4 | 6 |
|---|---|---|---|---|---|---|---|---|---|
| $`\lvert M_{\rm res}\rvert/\lvert M_{\rm cl}\rvert`$ | 0.0077 | 0.0466 | 0.1677 | 0.4726 | 0.7205 | 0.8235 | 0.9265 | 0.9608 | 0.9833 |

**Theorem K3.** For the Eckart barrier the spectrum-weighted budget ratio
saturates at 1 from below and never exceeds it.

*Proof sketch.* Per mode $`|\sin u - u|/u = 1 - \sin u/u`$ exceeds 1 for
$`u \in (\pi, 2\pi)`$, peaking at 1.2172 near $`u = 4.4935`$. But
$`\tilde V(k)`$ decays as $`e^{-\pi ka/2}`$, so those modes carry negligible
weight; and as $`y_{\max} \to \infty`$, $`M_{\rm res} \to -M_{\rm cl}`$, giving
the limit 1. $`\square`$

The contrast with §6.3 of the splitting note is the point. On a ring every
world sits at $`u = q\pi`$, the symbol vanishes and the residual exactly
cancels the classical term, pinning the ratio at 1: the reorganisation looks
empty. On the open line with a localised barrier the compensated channel is a
strict gain at every reach, two orders at $`y_{\max} = 0.1a`$.

Note that the exponential decay rate of $`\tilde V`$ is $`\pi a/2`$ — the same
number as the reach ceiling of Theorem K1, and for the same reason: both are
the distance to the nearest complex singularity.

---

## 5. Theorem K4: the deterministic step carries no tunnelling

The right observable is not the far-side weight but the **classical outcome
functional**

```math
\Sigma \;=\; \bigl\{\, E > V_0,\; p > 0 \,\bigr\}
        \;\cup\; \bigl\{\, E < V_0,\; r > 0 \,\bigr\},
\qquad E = \frac{p^2}{2\mu} + V(r).
```

A world clears the barrier moving right, or is already past it with too little
energy to come back. The two pieces do not meet: at $`r = 0`$, $`V = V_0`$
forces $`E \ge V_0`$, so the boundary between them has measure zero.

**Theorem K4.** $`\Sigma`$ is exactly invariant under streaming plus
deterministic acceleration. Hence the classical transmission is
$`\Sigma`$ evaluated on the *initial* state, and every part of the
quantum–classical gap in the transmission is delivered by $`K_{\rm res}`$.

*Proof.* Streaming plus the full classical force is the classical Hamiltonian
flow, which conserves $`E`$ and maps each trajectory to its own asymptotic
outcome; $`\Sigma`$ is by construction a union of whole trajectories. By
Theorem C1 the potential substep factorises exactly, so the only remaining
factor is $`\exp(\tau M_{\rm res})`$. $`\square`$

This device is inherited from the inverted pair barrier, whose transmission is
readable at $`t = 0`$ for the same reason — but there the residual is
identically zero (Theorem I4), so $`\Sigma`$ is conserved by the *whole*
dynamics and there is nothing to measure. Here it is conserved by the
classical half only, and its motion is exactly the quantum correction.

Numerically ($`V_0 = 1`$, $`a = 2`$, $`\beta = 2`$, $`\mu = 1/2`$,
$`\hbar = 1`$; Gaussian packet $`\sigma_r = a`$ centred at $`p_c = p_b`$ so the
classical answer is exactly $`1/2`$; grid $`1024 \times 256`$, $`\Delta p = 0.05`$,
$`dt = 0.02`$, $`t = 22`$):

| quantity | value |
|---|---|
| classical, $`= \Sigma`$ at $`t = 0`$ | 0.500003 |
| closed form, classical | 0.500023 |
| full symbol run, final | 0.544477 |
| compensated product run, final | 0.544473 |
| closed form $`\int|\phi(k)|^2 T(E)\thinspace dk`$ | 0.544156 |
| C1 check, $`\lvert`$ full $`-`$ compensated $`\rvert`$ | $`4.0\times10^{-6}`$ |
| closed-form check, $`\lvert`$ run $`-`$ exact $`\rvert`$ | $`3.2\times10^{-4}`$ |
| gap delivered by the residual | 0.044475 |
| closed-form gap | 0.044134 |

Conservation under the deterministic step alone, over the approach:

| $`t`$ | 0 | 2 | 4 | 6 | 8 |
|---|---|---|---|---|---|
| $`\Sigma`$ | 0.500003 | 0.499990 | 0.499619 | 0.499803 | 0.500310 |

One honest limitation. Beyond $`t \approx 10`$ a *classical* run filaments
below the grid — $`\min W`$ reaches $`-0.16`$, and the check loses meaning.
That is a property of the classical reference, not of the claim: the same grid
carries the quantum run to $`3\times10^{-4}`$ of the closed form, because the
residual channel smooths in $`p`$. It is also why Theorem K4 is worth having:
it removes any need to evolve a classical run at all.

---

## 6. Theorem K5: two large flows, one small difference

Because $`\sum_q K_q = 0`$, the rate of change of $`\Sigma`$ under the residual
channel can be written purely as boundary crossings,

```math
\frac{d}{dt}\!\int_\Sigma W
\;=\; \sum_q K_q(r) \sum_{p} \bigl[\mathbf{1}_\Sigma(p + \xi_q)
      - \mathbf{1}_\Sigma(p)\bigr] W(r, p)
\;=\; \Phi_{\rm in} - \Phi_{\rm out},
```

splitting the bracket into its $`+1`$ and $`-1`$ parts.

**Theorem K5.** The quantum correction to the transmission is the
time-integrated *imbalance* of two flows across the classical separatrix, each
several times larger than their difference.

Measured (part F, flux leg at $`dt = 0.005`$):

| quantity | value |
|---|---|
| gross traffic $`\int(\lvert\Phi_{\rm in}\rvert + \lvert\Phi_{\rm out}\rvert)dt`$ | 0.222774 |
| net transfer $`\int(\Phi_{\rm in} - \Phi_{\rm out})dt`$ | 0.042806 |
| gap $`T_{\rm full} - T_{\rm class}`$ | 0.044475 |
| net / gross | 0.192 |

The ledger closes to 3.8 per cent. The residual discrepancy is understood and
is not physical: the flux diagnostic reads the *generator* $`M_{\rm res}`$
while the step applies the *exponential*, an $`O(dt)`$ error that accumulates
coherently because it has a fixed sign. Measured per step at the same state,
the ratio of exact to ledger $`\Delta\Sigma`$ is 1.0136 at $`dt = 0.02`$,
1.0069 at $`dt = 0.01`$ and 1.0034 at $`dt = 0.005`$ — first order, as it
should be.

This is the ontological content of the note. A world does not pass *through*
the barrier. It streams on a Newtonian arc under the full classical force, and
the residual channel emits positon–negaton pairs at $`p \pm \xi_q`$; some
children land above the classical separatrix and some below. Tunnelling is
what is left over when the two flows almost, but not exactly, cancel. This is
the positon–negaton sea picture appearing on a problem with a closed-form
answer.

Two cautions on the drawing, expanded in §8. The residual kernel is signed, so by
consequence 3 of §2.2 of the splitting note **no one-body jump process
exists**; what is drawn is the pair-generation form, in which the parent
streams deterministically and at rate $`\Gamma(r) = \sum_q|K_q(r)|`$ emits a
pair of zero expected weight. And the sampler uses a hard window at the
ceiling reach $`y_{\max} = \pi a/2`$, so it is illustrative rather than exact.

---

## 7. Theorem K6: the cancellation tightens as $`1/\beta`$

Decomposing the gap in energy — worlds below threshold that transmit, against
worlds above threshold that reflect (part G, packet $`\sigma_r = a`$ centred at
$`p_b`$ throughout, so $`T_{\rm class} = 1/2`$):

| $`\beta`$ | $`T_{\rm class}`$ | $`T_{\rm quant}`$ | net | tunnel in | over-barrier reflection | net/gross | $`\times\beta`$ |
|---|---|---|---|---|---|---|---|
| 0.5 | 0.59429 | 0.75608 | +0.16179 | 0.18473 | 0.02295 | 0.779 | 0.390 |
| 1 | 0.51165 | 0.60293 | +0.09128 | 0.13398 | 0.04270 | 0.517 | 0.517 |
| 2 | 0.50001 | 0.54416 | +0.04414 | 0.10308 | 0.05894 | 0.272 | 0.545 |
| 4 | 0.50000 | 0.52184 | +0.02184 | 0.09043 | 0.06859 | 0.137 | 0.549 |
| 8 | 0.50001 | 0.51089 | +0.01088 | 0.08458 | 0.07370 | 0.069 | 0.550 |
| 16 | 0.49999 | 0.50544 | +0.00545 | 0.08177 | 0.07632 | 0.034 | 0.551 |

Gross traffic saturates near 0.158 while the net falls as $`1/\beta`$; the last
column is constant to three figures from $`\beta = 2`$ up.

**Theorem K6.** For a packet centred on the barrier top with
$`\sigma_r = a`$, net/gross $`\to c/\beta`$ with $`c \approx 0.55`$.

*Proof sketch.* Write the gap as $`\int w(p)\thinspace g(E(p))\thinspace dp`$ with
$`g = T - \Theta(E - V_0)`$, which is odd about $`E = V_0`$ in the Kemble
approximation. Changing variables,
$`\int g(E)\thinspace h(E)\thinspace dE`$ with $`h = \mu w/p`$, and only the part
of $`h`$ odd about $`V_0`$ survives. Since the packet is centred,
$`w'(p_b) = 0`$, so the leading term comes not from the packet's slope but
from the Jacobian: $`h'(V_0) = -\mu^2 w(p_b)/p_b^3`$. With
$`\hbar\Omega_b = 2V_0/\beta`$ for the barrier curvature and
$`w(p_b) \propto a/\hbar`$, the net scales as
$`\mu^2 w(p_b)(\hbar\Omega_b)^2/p_b^3 \propto 1/\beta`$. $`\square`$

So the imbalance is a *Jacobian effect*: the two flows are exactly opposed in
energy, and what survives is only the non-uniformity of the energy–momentum
map across the tunnelling window.

The consequence for the algorithm is a cost law. A signed ensemble must
resolve a signal that is a fraction $`\approx 0.55/\beta`$ of the traffic it
samples, so the particle count needed to resolve $`T`$ to fixed relative
accuracy grows as $`\beta^2`$. This is the transmission-observable analogue of
Theorem O2's $`S/\sqrt{D}`$ decay, and it says the semiclassical limit is the
expensive one — the opposite of the usual intuition.

---

## 8. World-particle identity, and where tunnelling comes from

### 8.1 A hop is not a jump

It is tempting to read the residual channel as a world *jumping* in
momentum: a phase point at $`(r, p)`$ that discontinuously becomes
$`(r, p + \xi_q)`$. That reading is not available. By consequence 3 of §2.2
of [`compensated_liouville_splitting.md`](compensated_liouville_splitting.md),
$`K_{\rm res}`$ is signed, so it is not a rate and there is no one-body
Markov jump process for it to generate.

What the channel does generate is *pair creation*. Because
$`\sum_q K_q = 0`$ and $`K_{-q} = -K_q`$, the operator can be realised
exactly as: the parent world streams on undisturbed, and at rate
$`\Gamma(r) = \sum_{q \ne 0}|K_q(r)|`$ a **positon–negaton pair** appears at
$`p \pm \xi_q`$, carrying opposite signs and hence zero expected weight.

The ontological consequence is worth stating plainly. **No world-particle
ever changes its momentum discontinuously.** Every world-particle in the
ensemble — the excess positons of the initial condition and every member of
the created sea alike — follows a genuine Newtonian worldline under the full
classical force, for its entire life. World-particle identity is never
broken by the quantum channel. What the quantum channel changes is *how many
worlds there are*, and with what signs.

This is the sense in which the compensated split earns its name
ontologically rather than numerically. The uncompensated kernel mixes the
classical force into the jump measure, so that "hops" include the ordinary
acceleration and no world has a Newtonian history. Under compensation the
deterministic step is the whole of the force, and the residual is the whole
of the quantum — and the residual acts only by creating and destroying, never
by moving.

### 8.2 Theorem K7: the emission bias reverses four times

**Theorem K7.** For $`V_0\thinspace\mathrm{sech}^2(r/a)`$ the emission rate
$`\Gamma`$ vanishes and the sign of $`K_q`$ reverses at three points,

```math
r \;=\; 0
\qquad\text{and}\qquad
r \;=\; \pm r_*,
\qquad
r_* \;=\; a\thinspace\mathrm{artanh}\sqrt{2/3}
      \;=\; a\ln\bigl(\sqrt2 + \sqrt3\bigr) \;\approx\; 1.1462\thinspace a,
```

so the barrier carries four emission lobes of alternating sign.

*Proof.* $`V''' \propto c\thinspace t\thinspace(2c - t^2)`$ with
$`c = \mathrm{sech}^2u`$, $`t = \tanh u`$, $`u = r/a`$. It vanishes at
$`t = 0`$ and where $`2(1 - t^2) = t^2`$, i.e. $`t = \sqrt{2/3}`$, whence
$`\sinh u = \sqrt2`$. By Theorem C7 these are quiet points at short reach,
and $`K_q`$ inherits the sign of $`V'''`$ at leading order in the reach.
$`\square`$

Numerically, with $`a = 2`$ so $`r_* = 2.2924`$:

| $`r`$ | $`-6`$ | $`-3`$ | $`-1.5`$ | $`0`$ | $`1.5`$ | $`3`$ | $`6`$ |
|---|---|---|---|---|---|---|---|
| $`V'''`$ | $`+0.0095`$ | $`+0.0749`$ | $`-0.2993`$ | $`0`$ | $`+0.2993`$ | $`-0.0749`$ | $`-0.0095`$ |
| $`\Gamma`$ | 0.291 | 0.180 | 4.304 | **0** | 4.304 | 0.180 | 0.291 |
| lobe | I | I | II | — | III | IV | IV |

The summit is quiet: **no pairs are emitted at the top of the barrier at
all**. The emission is concentrated on the flanks, peaking near
$`|r| \approx a/2`$, and reverses sign at $`\pm r_*`$.

Splitting the separatrix ledger of Theorem K5 by lobe:

| lobe | net transfer |
|---|---|
| I, $`r < -r_*`$ | $`-0.034500`$ |
| II, $`-r_* < r < 0`$ | $`+0.065978`$ |
| III, $`0 < r < r_*`$ | $`+0.011088`$ |
| IV, $`r > r_*`$ | $`-0.000457`$ |
| **total** | $`+0.042109`$ |
| sum of the lobe magnitudes | $`0.112022`$ |

So Theorem K5's cancellation is not a single opposition of two flows but a
**fourfold** one, and it is dominated by the approach flank: the outer lobe I
drives worlds *out* of the transmitting set as the packet decelerates, the
inner lobe II drives them back in and overshoots, and the far side barely
contributes because the packet crosses it quickly and the outcome is already
settled. This is visible directly in the time trace of §6: $`\Sigma`$ dips to
0.451 before rising to 0.545.

### 8.3 The origin of tunnelling, stated in one sentence

A world does not pass through the barrier. It runs up the flank on an
ordinary Newtonian worldline and, if its energy is short of $`V_0`$, it turns
around and leaves — always. What crosses is not that world but *other* worlds:
positon–negaton pairs conjured on the flanks, whose members follow their own
Newtonian worldlines from birth, some of them fast enough to clear the
summit. Tunnelling is the residue left when the positons and negatons that
clear the summit fail, by a fraction $`\approx 0.55/\beta`$, to cancel.

![Space-time worldlines and the four lobes](https://raw.githubusercontent.com/billpage/wpmw/output/figures/eckart_compensated_spacetime.png)

(a) the four emission lobes and the three quiet points of Theorem K7;
(b) space-time worldlines — four excess positons, two of which are classically
reflected, with the positon–negaton pairs they emit on the flanks, every path
Newtonian and every branch a creation rather than a displacement;
(c) the fourfold cancellation of the lobe ledger.

One caution about panel (b), which is worth more than it looks. The ten
emissions drawn there happen to put a *negative* signed weight on the far
side, whereas the true net is positive. That is not an error in the drawing;
it is Theorem K6 made visible. With net/gross $`\approx 0.19`$ at
$`\beta = 2`$, a handful of sampled emissions cannot resolve even the sign of
the effect, and the figure is an illustration of the mechanism, not an
estimate of its size.

---

## 9. Figures

![The split, and the budget ratio](https://raw.githubusercontent.com/billpage/wpmw/output/figures/eckart_compensated_split.png)

The barrier and its third derivative; the three symbols at $`x = 0.5`$ with the
analyticity strip shaded; the budget ratio of Theorem K3.

![The reach ceiling](https://raw.githubusercontent.com/billpage/wpmw/output/figures/eckart_compensated_ceiling.png)

The poles of $`V`$ in the complex separation plane bounding the reach disc; the
Moyal partial sums diverging beyond $`R(x)`$; the rung count of Corollary K1.1.

![Separatrix traffic](https://raw.githubusercontent.com/billpage/wpmw/output/figures/eckart_compensated_separatrix.png)

The main drawing. (a) the classical outcome set taken directly from the run,
with the hop channel carrying weight both ways across its boundary; (b) sample
world-particle paths — Newtonian arcs punctuated by positon–negaton emissions
at $`p \pm \xi_q`$, some children landing above the boundary and some below;
(c) $`\Phi_{\rm in}`$ and $`\Phi_{\rm out}`$ against time, two large flows whose
difference is small; (d) the cumulative imbalance added to $`T_{\rm class}`$,
landing on the closed-form transmission.

![The 1/beta law](https://raw.githubusercontent.com/billpage/wpmw/output/figures/eckart_compensated_scaling.png)

Gross traffic saturating while the net falls, and net/gross against
$`\beta`$ on log axes against $`0.55/\beta`$.

---

## 10. Open items

- **K-LS1.** The flux ledger closes to 3.8 per cent, limited by the
  generator-versus-exponential error of the diagnostic. A ledger built from
  the exact per-step $`\Delta\Sigma`$, decomposed by crossing direction, would
  close it exactly; the decomposition is not a simple correlation and has not
  been constructed.
- **K-LS2.** Corollary K1.1 forces $`\sigma_r < a/2`$ for a reach-limited
  lattice to resolve the packet. §5 therefore runs at unbounded reach. A
  genuinely reach-limited transmission calculation, with the fold of
  Theorem E2 applied to the residual only, has not been done — and the fold
  wants $`y_{\max} \gg a`$ while Theorem K1 wants $`y_{\max} < \pi a/2`$. The
  two are not obviously compatible.
- **K-LS3.** Theorem K6 is stated for a packet centred on the barrier top with
  $`\sigma_r = a`$. The scaling of the constant with the packet's offset and
  width has not been mapped, and the derivation's use of $`w'(p_b) = 0`$
  suggests an off-centre packet has a *larger* net at fixed gross.
- **K-LS4.** The whole of §5–7 uses a positive initial Wigner function. The
  §5 variant of the inverted-pair note — an initial state with genuine
  negativity straddling the separatrix — should worsen the cancellation, but
  by how much is unmeasured.
- **K-LS6.** Theorem K7 is stated for sech². The number of emission lobes
  is the number of sign changes of $`V'''`$, so it is potential-specific; a
  general statement relating lobe count to the cancellation depth of
  Theorem K5 has not been attempted.
- **K-LS5.** The soft-core reading of Corollary K1.2 suggests the next test
  case: soft-core Coulomb has a tunable ceiling $`\epsilon`$ and, unlike the
  Eckart barrier, an attractive well with bound states.
