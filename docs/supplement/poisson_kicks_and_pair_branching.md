# Poisson Kicks and Pair Branching: Local Jump Rules for the Wigner Equation

**A test of D. Cyganski's "virtual-photon" proposal from the call of 29
September 2026: can the quantum part of the Wigner equation be driven by a
Poisson clock, with momentum kicks whose size and rate are read off the
potential at the particle's position? Two answers, one exact and one
numerical. Exactly: a positive one-particle clock always adds momentum
variance, which shows up both as heating and as decoherence; a *signed*
two-point kick is exact only for single-mode potentials; and the classical
force is the zero-step limit of signed kicks. Numerically: an independent
event-driven re-implementation of his pair branching is unbiased, and the
error it does have comes from annihilation, not from the clock.**

---

## 0. Status and provenance

On the call of 2026-09-29 Cyganski showed three things:

1. **A slide.** The spawning rule for a cosine potential, the shot-noise
   stochastic differential equation of a photodetector, and a proposed
   one-particle kick equation driven by a Poisson counter $`N_\lambda`$,
   ```math
   dp \;=\; \hbar k_0\thinspace\mathrm{sgn}\bigl[\sin k_0x\bigr]\thinspace dN_\lambda,
   \qquad \lambda \propto \lvert\sin k_0 x\rvert .
   ```
2. **A notebook**, `Virtual_Photon_Pair_Polarization_Wigner_Test`, written
   with an AI assistant. It evolves a Gaussian in
   $`V = \tfrac12 m\omega^2x^2 + V_1\sin k_0x`$, with the quadratic part as
   curved trajectories and the sine as event-driven *pair branching*. Its
   infinite-particle limit matches the QLE to about $`10^{-6}`$, and its
   *one-particle* Poisson limit washes out at the return turning point. A
   finite run with $`N_0 = 1.6\times10^6`$ overshoots both peaks of the
   position marginal by roughly 13%.
3. **A handwritten formula for a general potential**, not yet tested:
   ```math
   dp \;=\; \hbar\left\lvert\frac{V''}{V'}\right\rvert dN_\lambda,
   \qquad \lambda = \frac{V'^{2}}{V''} .
   ```

This note tests all three. The notebook's code was not available, so §5–§7
are an **independent re-implementation**, not a replication. Its parameters
are this note's own. They are chosen to resemble his figures, not copied from
them.

**What this note inherits.**

- Lemma T1 of [`takabayasi_1954_stochastic_picture.md`](takabayasi_1954_stochastic_picture.md): for one Fourier mode, the
  Wigner kernel is a pair of opposite-signed point masses at
  $`\pm\hbar k/2`$.
- Proposition T3 of the same note: no one-body Markov process has the QLE as
  its generator.
- Takabayasi's (3.17), adopted in §6(a) of the same note: the Kramers–Moyal
  moments of the Wigner kernel in closed form.
- Proposition F4 of [`four_action_foundations.md`](four_action_foundations.md): shrinking the jump is the
  same operation as sending $`\hbar\to 0`$.
- Theorems D5 and D6 of [`../analysis/species_sectors_and_annihilation.md`](../analysis/species_sectors_and_annihilation.md):
  cell-exact annihilation is unbiased only when cancelled particles coincide,
  and "soft" annihilation is an artificial decoherence.

**What it adds.**

- Proposition X1: every positive one-particle rule diffuses in momentum,
  cannot carry a quantum correction without diffusing, and cannot create
  negative regions of $`W`$. Cyganski's general rule diffuses at
  $`M_2 = \hbar\lvert V''\rvert`$.
- Theorem X2: a single signed two-point kernel is exact at every order iff
  $`V''' = rV'`$ with $`r`$ constant.
- Proposition X3: the classical force is the zero-step, infinite-rate member
  of that family.
- Proposition X4: an exact heating law for any jump rule,
  $`d\langle H\rangle/dt = \langle M_2\rangle/2m`$.
- Proposition X6: momentum diffusion is decoherence. A positive kernel damps
  the coherence between legs a distance $`y`$ apart; a signed pair is a pure
  phase.
- Proposition X5 (numerical): in an event-driven run with curved
  trajectories, the bias is set by the momentum width of the annihilation
  cell.

Nothing here retracts an earlier result.

Theorem prefix **X**, open items **X-SP**. Companion code:
[`src/demo_signed_local_kernel.py`](../../src/demo_signed_local_kernel.py) (Parts A–E, SymPy and
spectral) and [`src/demo_event_driven_branching.py`](../../src/demo_event_driven_branching.py) (Parts F–L, Monte
Carlo), with the simulator in
[`src/wpmwlib/event_branching.py`](../../src/wpmwlib/event_branching.py). Units $`\hbar = m = 1`$ in all
numerical parts. Every number below is an output of those scripts.

---

## 1. Kicks and moments: what the QLE asks of a jump process

This section and the next are written as a tutorial. They introduce the
three ideas the rest of the note relies on: the *moments* of a kick, *momentum
diffusion*, and what *heating* means when there is no heat bath.

### 1.1 A jump process in momentum

Fix a position $`x`$ and let particles there receive sudden momentum kicks.
Describe the kicks by a **rate kernel** $`K(x,\Delta)`$: the rate at which a
particle at $`x`$ jumps from $`p`$ to $`p+\Delta`$. For ordinary particles
$`K\ge 0`$. The signed ensembles of this project also allow $`K<0`$, which
means the jump deposits a negaton instead of a positon.

The total rate is $`\Lambda(x) = \int K\thinspace d\Delta`$. For infinitely
many particles the density obeys the **master equation**

```math
\partial_t W\big|_{\rm kicks} \;=\; \int K(x,\Delta)\thinspace W(x,\thinspace p-\Delta)\thinspace d\Delta \;-\; \Lambda(x)\thinspace W(x,p) .
```

The first term is the gain from particles arriving at $`p`$; the second is
the loss of particles leaving it.

### 1.2 Moments: the mean kick and the spread of kicks

The **moments** of the kernel summarise what the kicks do:

```math
M_n(x) \;=\; \int \Delta^{n}\thinspace K(x,\Delta)\thinspace d\Delta .
```

Taylor-expanding $`W(p-\Delta)`$ in $`\Delta`$ turns the master equation
into the **Kramers–Moyal series**

```math
\partial_t W\big|_{\rm kicks} \;=\; -\partial_p\bigl(M_1W\bigr) \;+\; \tfrac12\thinspace\partial_p^{\thinspace 2}\bigl(M_2W\bigr) \;-\; \tfrac16\thinspace\partial_p^{\thinspace 3}\bigl(M_3W\bigr) \;+\; \cdots
```

Each term has a plain meaning:

- **$`M_1`$ is the mean kick per unit time**, i.e. a force. The term
  $`-\partial_p(M_1W)`$ is the classical force term of Liouville's equation.
- **$`M_2`$ is the spread of the kicks per unit time.** The term
  $`\tfrac12\partial_p^2(M_2W)`$ is a diffusion equation in momentum, with
  **momentum-diffusion coefficient** $`D = M_2/2`$. It is the momentum-space
  version of Brownian motion: a random walk whose steps widen the
  distribution without moving its centre.
- **$`M_3`$ and higher** are the skew, the tails, and so on. In the QLE
  these carry the quantum corrections.

The same content is visible in the first two moments of $`p`$. From the
master equation, the kicks alone change them at the rates

```math
\frac{d\langle p\rangle}{dt}\bigg|_{\rm kicks} = \langle M_1\rangle,
\qquad
\frac{d\langle p^2\rangle}{dt}\bigg|_{\rm kicks} = 2\langle p\thinspace M_1\rangle + \langle M_2\rangle .
```

The first equation is Newton's second law, averaged. The second says that,
beyond what the force does, the kicks add $`\langle M_2\rangle`$ per unit time
to $`\langle p^2\rangle`$. That is, they add it to the **momentum variance**.

### 1.3 What the QLE demands

The potential term of the Wigner equation is the Moyal series. Matching it
term by term to the Kramers–Moyal series gives

```math
M_{2s+1} \;=\; (-1)^{s+1}\Bigl(\frac{\hbar}{2}\Bigr)^{2s} V^{(2s+1)}(x),
\qquad M_{2s} = 0 .
\qquad\qquad (1.1)
```

(1.1) is Takabayasi's (3.17) with $`\varepsilon = \hbar^2/4`$. Part A
re-derives it from the exact symbol
$`i\bigl[V(x+\hbar\kappa/2)-V(x-\hbar\kappa/2)\bigr]/\hbar`$, with $`V`$ given
as symbolic Taylor coefficients up to $`n=7`$. The residual is 0 in SymPy.

In words:

- $`M_1 = -V'`$: the mean kick is the classical force.
- $`M_2 = 0`$: **no momentum diffusion at all.** The QLE never adds variance
  of its own; it only reshapes the distribution.
- $`M_3 = \tfrac{\hbar^2}{4}V'''`$: the first quantum correction is a *skew*
  of the kicks.

The whole difficulty is to have $`M_3\neq 0`$ while keeping $`M_2 = 0`$. The
two kernels this note compares show how:

```text
  one positive kick (§2)              signed pair (§3)
        +lambda                              +w
          |                                   |
  --------+---------|--------        ---------+---------|---------+----
          0        +a                        -a         0        +a
                                                                   |
                                                                  -w
  mean      M1 = lambda a  = F             M1 = 2 w a        = F
  spread    M2 = lambda a^2 > 0            M2 = w a^2 - w a^2 = 0
  skew      M3 = lambda a^3                M3 = 2 w a^3
```

Both deliver the same mean push. For the pair, the spreads of the two signed
children cancel exactly. It works like a dipole: a net force with no spreading.

---

## 2. One-particle Poisson kicks

### 2.1 The rule and its infinite-particle limit

Cyganski's slide reads the photon as a kick delivered to *one* particle.
Every particle is positive; nothing is born or annihilated. At the times of a
Poisson counter of rate $`\lambda(x)`$, the particle's momentum jumps by
$`a\thinspace\mathrm{sgn}F(x)`$, where $`F = -V_1'`$ is the force of the
perturbing potential. Choosing $`\lambda = \lvert F\rvert/a`$ makes the mean
kick right. On the slide $`a = \hbar k`$; the sine's own kernel uses
$`a = \hbar k/2`$.

"The one-particle limit" is the master equation of §1.1 for this kernel:

```math
\partial_t W\big|_{\rm kicks} \;=\; \lambda(x)\thinspace\Bigl[W\bigl(x,\thinspace p - a\thinspace\mathrm{sgn}F\bigr) - W(x,p)\Bigr] .
```

That is the equation behind the notebook's middle panel. Part G integrates it
directly.

### 2.2 Why a positive rule must diffuse

> **Proposition X1 (positive rules diffuse).** Let $`K\ge 0`$.
>
> 1. $`M_2 \ge M_1^{2}/\Lambda`$. So wherever $`V'\neq 0`$, a positive rule
>    has $`M_2>0`$ and cannot satisfy (1.1).
> 2. $`M_3^{2} \le M_2\thinspace M_4`$. So a positive rule can carry the
>    quantum term $`M_3\neq 0`$ only if it also diffuses.
> 3. A positive kernel maps $`W\ge 0`$ to $`W\ge 0`$. So it can never create
>    the negative regions of a Wigner function.
>
> For Cyganski's general rule — step $`\hbar\lvert V''/V'\rvert`$ along the
> force, at rate $`V'^{2}/(\hbar\lvert V''\rvert)`$ — the moments are, exactly,
> ```math
> M_1 = -V', \qquad M_2 = \hbar\thinspace\lvert V''\rvert, \qquad
> M_3 = -\mathrm{sgn}(V')\thinspace\frac{\hbar^{2}\thinspace V''^{\thinspace 2}}{\lvert V'\rvert} .
> ```

*Proof.* Parts 1 and 2 are the Cauchy–Schwarz inequality for the measure
$`K\thinspace d\Delta`$, applied to $`(1,\Delta)`$ and to $`(\Delta,\Delta^2)`$.
Part 3 holds because both terms of the master equation move positive weight
to positive weight. The moments are computed in SymPy in Part B, in all four
sign cases. $`\blacksquare`$

Part 2 is Pawula's theorem (1967) in its smallest form. For a single kick
size it is an equality, $`M_3^2 = M_2M_4`$ (Part B2).

Shrinking the step does not escape X1. Let $`a\to 0`$ with
$`\lambda a = F`$ held fixed. Then $`M_2 = Fa\to 0`$, but so does
$`M_3 = Fa^2`$. The limit is classical Liouville dynamics, not quantum
mechanics. A positive clock is therefore either diffusive or classical,
never quantum.

**Cyganski's general rule.** The rate needs the factor $`1/\hbar`$ that the
handwritten formula omits: $`V'^2/V''`$ is an energy, not a rate. With it
restored the drift is right for every $`V`$, and the diffusion is wrong for
every $`V`$ (Part B). Two cases show how badly:

- **Harmonic oscillator.** The rule adds diffusion $`M_2 = \hbar m\omega^2`$,
  with a step $`\hbar/\lvert x\rvert`$ that diverges at the origin. The QLE
  for this potential is exactly classical.
- **Cosine.** The rule gives a step $`\hbar k\lvert\cot kx\rvert`$ rather than
  $`\hbar k/2`$. It does not even reproduce the special case it generalises.

Proposition T3 already said that no one-body Markov process generates the
QLE. X1 says how the failure shows up: as momentum diffusion. The next two
subsections give its two physical faces, heating and decoherence.

### 2.3 What "heating" means here

There is no heat bath and no temperature in this problem. "Heating" means
only this: **the mean energy of the ensemble grows without bound, because
the kicks keep adding momentum variance.** It is the same *stochastic
heating* a charged particle suffers in a randomly fluctuating field.

In a thermal bath, diffusion is balanced by friction (the Einstein relation
$`D = m\gamma k_BT`$) and the energy settles down. Here there is diffusion
and no friction, so nothing stops the growth.

> **Proposition X4 (heating law).** Let $`H = p^2/2m + V_0 + V_1`$. Evolve
> $`V_0`$ as a classical flow and $`V_1`$ by any jump kernel with
> $`M_1 = -V_1'`$. Then, for every state,
> ```math
> \frac{d\langle H\rangle}{dt} \;=\; \frac{\langle M_2(x)\rangle}{2m}
> \;=\; \frac{\langle D(x)\rangle}{m} .
> ```

*Proof.* By §1.2 the jumps change $`\langle p^2/2m\rangle`$ at rate
$`\langle (pM_1 + M_2/2)/m\rangle`$. The flow of $`V_0`$ changes it at rate
$`-\langle pV_0'/m\rangle`$. Advection changes
$`\langle V_0+V_1\rangle`$ at rate $`\langle p(V_0'+V_1')/m\rangle`$. With
$`M_1 = -V_1'`$, everything cancels except $`\langle M_2\rangle/2m`$.
$`\blacksquare`$

Only $`M_1`$ and $`M_2`$ enter, so the law holds whatever the higher moments
do. The force does work, but it is exactly the work the potential gives
back. The *only* net energy input is the variance the kicks inject. For
signed pair branching $`M_2 = 0`$, and the energy is conserved in
expectation, as it is for the QLE.

### 2.4 Momentum diffusion is decoherence

Heating is the energy face of momentum diffusion. Its other face is the loss
of interference.

Write the density matrix at leg separation $`y`$, the variable conjugate to
$`p`$. A shift $`W(p-\Delta)`$ becomes multiplication by
$`e^{-i\Delta y/\hbar}`$, so a kernel acts on $`\rho(x,y)`$ through its
**symbol**

```math
S(x,y) \;=\; \int K(x,\Delta)\thinspace\Bigl(e^{-i\Delta y/\hbar} - 1\Bigr)\thinspace d\Delta ,
\qquad
\partial_t\rho(x,y)\big|_{\rm kicks} = S(x,y)\thinspace\rho(x,y) .
```

The imaginary part of $`S`$ is a phase, which is ordinary unitary
evolution. The real part changes the *size* of the coherence between two
legs a distance $`y`$ apart.

> **Proposition X6 (diffusion is decoherence).**
>
> 1. For a positive kernel,
>    $`\mathrm{Re}\thinspace S = \int K\thinspace(\cos(\Delta y/\hbar) - 1)\thinspace d\Delta \le 0`$,
>    so coherence decays. For a single kick of size $`a`$ at rate $`\lambda`$,
>    ```math
>    \mathrm{Re}\thinspace S(y) = -\lambda\Bigl(1-\cos\frac{a y}{\hbar}\Bigr)
>    = -\frac{D\thinspace y^2}{\hbar^2} + O(y^4),
>    \qquad D = \frac{M_2}{2} = \frac{\lambda a^2}{2} ,
>    ```
>    and legs a distance $`y`$ apart lose coherence as
>    $`e^{-Dy^2t/\hbar^2}`$ at short separations.
> 2. For the signed pair, $`S(y) = -2iw\sin(ay/\hbar)`$ is purely imaginary.
>    It is a pure phase, and no coherence is lost.

*Proof.* Direct substitution; Part B2 does it in SymPy, including the series.
$`\blacksquare`$

So an odd, signed kernel generates unitary motion, while any positive kernel
has a real part that damps interference. This is the same statement as
Theorem D6 of
[`../analysis/species_sectors_and_annihilation.md`](../analysis/species_sectors_and_annihilation.md) — smoothing in
momentum is damping in $`y`$ — reached from the dynamics rather than from
annihilation. §6 meets it a third time.

Two measures show the one-particle rule destroying the state's interference
(Part G). The purity is $`\mathrm{Tr}\thinspace\rho^2 = 2\pi\hbar\int W^2\thinspace dx\thinspace dp`$,
which is 1 for a pure state. The negative mass is $`\int_{W<0}\lvert W\rvert`$.
The one-particle rule uses $`a = \hbar k`$.

| $`t`$ | purity, QLE | purity, one-particle | negative mass, QLE | negative mass, one-particle |
|---|---|---|---|---|
| 9 | 1.0000 | 0.1188 | 0.1618 | $`2\times10^{-4}`$ |
| 16 | 1.0000 | 0.0623 | 0.1804 | $`8\times10^{-5}`$ |

The QLE keeps a pure state pure, with a negative mass of 0.16–0.18 throughout.
The one-particle limit drops to 6% purity. Its negative weight is only the
ringing of the spectral solver, as part 3 of X1 requires.

The diffusion is not destroyed in the pair; it is moved. Counted *without*
signs, the pair's children do spread: $`\sum\lvert w\rvert\Delta^2 = 2wa^2`$.
It is the signed sum that cancels. The spreading is paid for in the size of
the signed population, not in $`W`$. That population is what annihilation
has to control (§6).

### 2.5 Numbers, and why the damage shows at the return point

The potential is $`\omega = 0.4`$, $`V_1 = 0.5`$, $`k = 1`$. The initial
state is a Gaussian at rest at $`x_0 = 6`$ with $`\sigma_x = 0.7`$, so
$`\sigma_p^2 = 0.51`$ (Parts G and K).

- The QLE conserves $`\langle H\rangle = 3.0650`$ to $`2\times10^{-5}`$.
- With the slide's step $`\hbar k`$ at rate $`\lvert\Gamma(x)\rvert`$, the
  infinite-$`N`$ master equation heats by 2.8864 by $`t=18`$. X4 predicts
  2.8870.
- With step $`\hbar k/2`$ the heating is 1.4404, against a prediction of
  1.4407. Halving the step halves $`M_2 = a\lvert F\rvert`$, and it halves the
  heating; it does not remove it.
- The one-particle Monte Carlo run agrees with its own master equation
  (with $`L^2 = 0.0076`$), not with the QLE (with $`L^2 = 0.2158`$). Its
  energy rises from 3.07 to 4.53 by $`t = 9`$, which matches X4 evaluated on
  its own positions.

Heating accumulates, so the damage is largest at the latest time shown. By
$`t = 9`$ the one-particle position marginal is off by
$`\lVert\rho-\rho_{\rm QLE}\rVert_1 = 0.63`$. At the return turning point
$`t = 16`$ it is off by 1.17. The energy injected by $`t = 18`$ (2.89)
is about eleven times the initial kinetic energy
$`\sigma_p^2/2m\approx 0.26`$, so the
packet is spread over a large area of phase space and its interference
pattern is gone. That is the washed-out panel of the notebook, reproduced from
the equations alone (figure 2, panels d and f).

---

## 3. The signed two-point kernel

Put weight $`+w`$ at $`+a`$ and $`-w`$ at $`-a`$: a signed pair, total rate
zero. Every even moment vanishes identically and $`M_{\rm odd} = 2wa^{n}`$.
Matching $`M_1`$ and $`M_3`$ fixes the rule:

```math
a^2 \;=\; -\frac{\hbar^{2}}{4}\thinspace\frac{V'''}{V'},
\qquad 2wa \;=\; -V' .
\qquad\qquad (3.1)
```

The higher moments then leave residuals (Part C):

```math
M_5^{\rm model}-M_5 \;=\; \frac{\hbar^4}{16}\thinspace\frac{V'V^{(5)}-V'''^{\thinspace 2}}{V'},
\qquad
M_7^{\rm model}-M_7 \;=\; -\frac{\hbar^6}{64}\thinspace\frac{V'^{\thinspace 2}V^{(7)}-V'''^{\thinspace 3}}{V'^{\thinspace 2}} .
```

> **Theorem X2 (when one local kernel is exact).** A single signed two-point
> kernel reproduces (1.1) at every order at every $`x`$ iff
> $`V''' = r\thinspace V'`$ with $`r`$ **constant**. There are three cases:
>
> | $`r`$ | potential | step $`a`$ | reading |
> |---|---|---|---|
> | $`-k^2`$ | $`C + A\cos(kx+\phi)`$ | $`\hbar k/2`$ | Lemma T1 |
> | $`0`$ | quadratic | $`0`$ | the classical force (Proposition X3) |
> | $`+\kappa^2`$ | $`C + A\cosh(\kappa x+\phi)`$, exponentials | $`i\hbar\kappa/2`$ | imaginary jump |

*Proof.* Write $`u = V'`$ and $`r(x) = u''/u`$. Exactness at every order
requires $`u^{(2s)} = r^s u`$ for all $`s`$.

- Differentiating $`u'' = ru`$ twice and imposing the $`s=2`$ condition
  $`u'''' = r^2u`$ gives $`(r'u^2)' = 0`$, so $`r' = C/u^2`$.
- Substituting into $`s = 3`$ leaves
  $`u^{(6)} - r^3u = 2C^2/u^3`$. SymPy computes this in Part C.
- So $`C = 0`$ and $`r`$ is constant. Conversely, constant $`r`$ satisfies
  every condition.

$`\blacksquare`$

The third row is the "imaginary jump" that Cyganski's *Extended
Fokker–Planck* memo contemplated for moment matching. More generally, (3.1)
makes the step imaginary wherever $`V'''/V' > 0`$. For $`x^4`$ that is
everywhere: $`a^2 = -3\hbar^2/2x^2`$.

**The structural reason.** The Wigner generator is *linear* in $`V`$. Rule
(3.1) is not: $`a`$ depends on the ratio $`V'''/V'`$. So the kernel of a sum
of potentials is the sum of the kernels, while the local step of a sum is
nothing in particular. For $`V = \tfrac12\omega^2x^2 + V_1\sin kx`$ (Part D):

```math
a^2(x) = \frac{\hbar^2}{4}\thinspace\frac{V_1k^3\cos kx}{\omega^2x + V_1k\cos kx},
\qquad
\frac{\delta M_5}{M_5} = -\frac{\omega^2x}{\omega^2x + V_1k\cos kx} .
```

At the parameters of §2:

- The step is real on only 60% of $`x\in[-9,9]`$.
- It diverges at every classical equilibrium, where $`V' = 0`$.
- For $`\lvert x\rvert > 3`$ it gets $`M_5`$ wrong by a median 100%.

Superposing instead gives the exact answer at every order, with residual
$`0`$ in SymPy through $`M_7`$. Superposing means the quadratic's own kernel,
which is the force, plus the sine's own kernel, a fixed $`\pm\hbar k/2`$
pair. **That superposition is the compensated split** of
[`../analysis/compensated_liouville_splitting.md`](../analysis/compensated_liouville_splitting.md), in its simplest
instance. For a general potential the per-mode kernels are summed by the
reach-limited residual. No finite set of local derivatives can replace that
sum, because (1.1) needs every odd derivative.

![Local signed kernel: domain, M5 error, dipole limit](https://raw.githubusercontent.com/billpage/wpmw/output/figures/signed_local_kernel.png)

*Figure 1. (a) The local signed step (3.1) for quadratic + sine: imaginary on
the shaded intervals, divergent at the equilibria. (b) Its relative $`M_5`$
error, which approaches $`-1`$ wherever the harmonic force dominates; the
split has zero error. (c) The dipole limit on a harmonic cat state: $`L^2`$
error against the exact rotation falls as $`a^2`$.*

---

## 4. Curved trajectories are the zero-step limit of pair branching

> **Proposition X3 (the force as a dipole of signed kicks).** Fix
> $`2wa = F(x)`$. Then
> ```math
> \frac{F}{2a}\Bigl[W(p+a)-W(p-a)\Bigr]
> \;=\; F\thinspace\partial_pW \;+\; \frac{Fa^{2}}{6}\thinspace\partial_p^{3}W \;+\; O(a^4),
> ```
> so the signed kernel tends to the classical force as $`a\to 0`$. The
> error is a spurious $`M_3 = Fa^2`$. Each body gives birth to children at
> rate $`\lvert F\rvert/a`$.

The symbol is $`iF\sin(\kappa a)/a = iF\kappa - iF\kappa^3a^2/6 + \dots`$
(Part E).

Part E evolves an even harmonic cat state ($`x_0 = \pm2.5`$) for a quarter
period and compares with the exact rotation. The error over
$`a = 0.1 \to 0.025`$ scales with exponent 2.006:

| $`a`$ | 0.4 | 0.2 | 0.1 | 0.05 | 0.025 | 0 (force) |
|---|---|---|---|---|---|---|
| $`L^2`$ error | 0.293 | 0.102 | 0.0257 | 0.0064 | 0.0016 | $`2.3\times10^{-6}`$ |
| births per body per unit time at $`x=2.5`$ | 6.2 | 12.5 | 25 | 50 | 100 | — |

Proposition F4 already showed that shrinking the jump is the classical
limit. X3 adds three things:

- For a quadratic potential the step $`a=0`$ is the **exact** member of the
  family (Theorem X2, $`r=0`$), not an approximation.
- The cost of approaching it by events diverges like $`1/a`$, while the error
  falls only like $`a^2`$. Reaching accuracy $`\epsilon`$ by kicks costs
  $`\propto\epsilon^{-1/2}`$ births per body per unit time. Doing it as a flow
  costs none.
- That is the operational content of the question from the call: why remove
  the curved part analytically? The part removed is exactly the member of
  the jump family with zero step and infinite rate, which is the one member a
  jump code handles worst. Cyganski's 2020 parity trick (quadratic as force,
  sine as spawning) was this split, reached independently.

---

## 5. An event-driven re-implementation

**The model.** Take $`V = \tfrac12\omega^2x^2 + V_1\sin kx`$ with the
parameters of §2. The quadratic part is an exact rigid rotation of each
particle in phase space, so no time step enters. The sine acts only through
events, at rate
$`\lvert\Gamma(x)\rvert = (V_1/\hbar)\lvert\cos kx\rvert`$ (Lemma T1).

**Pair branching.** At an event the parent persists. Two children are born
at the same $`x`$: sign $`s\thinspace\mathrm{sgn}\Gamma`$ at $`p - \hbar k/2`$
and sign $`-s\thinspace\mathrm{sgn}\Gamma`$ at $`p + \hbar k/2`$. In
[`emission_and_absorption.md`](emission_and_absorption.md)'s vocabulary this is the event channel.
Cyganski's "virtual photon ionising a sea pair" is its physical reading.

**The clock.** The rate changes along each flight, so an exponential waiting
time drawn at the start of a flight is not the right law. The simulator uses
**thinning** (Lewis & Shedler 1979). Candidate times are drawn at the global
bound $`\lambda_{\max} = V_1/\hbar`$, and each is accepted with probability
$`\lvert\cos kx(t)\rvert`$ at the particle's position *at that time*. This is
exact for any rate bounded by $`\lambda_{\max}`$.

**Annihilation.** At synchronisation times (every 0.1), opposite signs cancel
inside cells of size $`h_x\times h_p`$. Survivors of the majority sign keep
their own positions.

**Two references** (Part F):

1. A split-operator Schrödinger run, giving $`\lvert\psi\rvert^2`$.
2. A spectral QLE solver on $`(x,p)`$.

They agree on the position marginal to $`3.8\times10^{-6}`$ at $`t = 4.5`$
and $`1.1\times10^{-5}`$ at $`t = 9`$.

**Branching alone is unbiased** (Part H). Without annihilation, run to
$`t = 4.5`$, the population grows to 20 times $`N_0`$. The error times
$`\sqrt{N_0}`$ stays constant (13.7, 14.7, 13.9 for
$`N_0 = 6250\to 10^5`$), which is pure noise. So the pair rule and the
thinning clock are exact. Anything else the runs show comes from somewhere
else.

---

## 6. Annihilation is where the bias comes from

With annihilation, run to $`t = 9`$ (Part I). The table gives the rms
$`L^2`$ error of $`\rho(x)`$, with $`h_x = 0.1`$ throughout:

| $`N_0`$ | $`h_p = 0.5`$ | $`h_p = 0.25`$ | $`h_p = 0.125`$ |
|---|---|---|---|
| 25 000 (8 seeds) | 0.0548 | 0.0357 | 0.0346 |
| 100 000 (4) | 0.0490 | 0.0229 | 0.0183 |
| 400 000 (2) | 0.0490 | 0.0212 | 0.0134 |
| 1 600 000 (1) | 0.0487 | 0.0200 | 0.0094 |
| fitted bias $`b`$ in $`e^2 = b^2 + c/N_0`$ | 0.0482 | 0.0191 | 0.0090 |
| peak ratio at $`N_0 = 1.6\times10^6`$ | 0.911 | 0.983 | 1.000 |
| peak population $`/N_0`$ | 1.07 | 1.16 | 1.27 |

> **Proposition X5 (numerical).** With curved trajectories, cell
> annihilation has a bias that does not fall with $`N_0`$. It is set by the
> momentum width of the cell, falling with an exponent of 1.2 over
> $`h_p = 0.5\to0.125`$. It is insensitive to $`h_x`$ and to the
> synchronisation interval (Part J), and a first-moment correction does not
> remove it.

Part J separates the candidate causes:

- **$`h_p`$ versus $`h_x`$.** At $`h_p = 0.5`$ the error is 0.0489 at
  $`h_x = 0.05`$ and 0.0485 at $`h_x = 0.2`$. At $`h_p = 0.125`$ it is 0.0129
  at $`h_x = 0.2`$ and 0.0192 at $`h_x = 0.4`$.
- **Synchronisation interval.** At cell $`0.1\times0.25`$, running every
  $`0.05`$, $`0.4`$ and $`3.0`$ gives 0.0218, 0.0212 and 0.0186.
- **First moments.** A variant that shifts survivors so that each cell's
  signed first moments are preserved does no better: 0.0607 against 0.0485 at $`0.2\times0.5`$, and 0.0293 against 0.0212 at $`0.1\times0.25`$ (the $`N_0 = 4\times10^5`$ row of Part I). It is worse, not better.

**Why.** This is the third appearance of the idea of §2.4. Theorem D6 says that smoothing $`W`$ in momentum with width
$`\sigma_p`$ is identical to damping the density matrix's off-diagonal
coherence at separation $`Y`$ by $`e^{-\sigma_p^2Y^2/2\hbar^2}`$. That is a
decoherence of length $`\hbar/\sigma_p`$. Corollary D5.2 says cell
annihilation is exact only when the cancelled particles *coincide*. On the
crystal lattice they do. Here they do not, because the classical force moves
$`p`$ continuously.

So annihilation inside a cell of width $`h_p`$ relocates cancelled weight
across interference fringes finer than the cell. That erases coherence
between legs further apart than about $`\hbar/h_p`$. Three consequences
follow:

- Preserving means cannot help, because the damage is at second order and
  beyond, not a displacement.
- The bias is the same whether one sweep or thirty cancels a given negaton.
- It *lowers* the peaks, as decoherence does: the ratio is 0.911 at
  $`h_p = 0.5`$.

**The tension this exposes.** Curved trajectories and exact annihilation are
incompatible:

- The force makes momentum continuous, so coincidence has measure zero.
- The compensated algorithm keeps coincidence by working on the lattice.
- An event-driven code with exact flows has to choose a cell, and pays
  $`b\approx 0.05\to0.009`$ for $`h_p = 0.5\to0.125`$ here, for a population
  that grows only from $`1.07N_0`$ to $`1.27N_0`$.

Note what "on the WPMW lattice" means here. The crystal-lattice step of the
sine mode, $`\hbar k/2 = 0.5`$, is the *coarsest* cell in the table, and the
most biased. Proposition T2 makes that lattice exact only for states of the
potential's period. The present state is not periodic, and Proposition R3 of
[`representation_cost_and_annihilation.md`](representation_cost_and_annihilation.md) says the momentum resolution must
follow the state's extent. The state spans about 7 in $`x`$ at $`t = 9`$, and
the Nyquist bound $`\pi\hbar/7\approx0.45`$ is exactly the scale at which the
bias becomes visible.

**Against the notebook.** The notebook reports a *positive* overshoot of
about 13% at the peaks of $`\rho(x)`$. Coarse annihilation here gives an
*undershoot*. The mechanism of X5 therefore does not explain his figure, and
nothing in this implementation does: at $`N_0 = 1.6\times10^6`$ and
$`h_p = 0.125`$ the peak ratio is 1.000. The diagnostic that would settle it
is the one in Part I: repeat the run at two or three $`N_0`$ and two cell
sizes. A plateau that does not move with $`N_0`$ is bias. Its sign then says
whether the cause is annihilation or something in the estimator or the
reference.

![Event-driven pair branching against the QLE](https://raw.githubusercontent.com/billpage/wpmw/output/figures/event_branching_phase_space.png)

*Figure 2. (a–c) $`t = 9`$: the QLE, the event-driven run
($`N_0 = 1.6\times10^6`$, $`h_p = 0.125`$) and their difference, on
$`0.2\times0.25`$ display bins. (d) The one-particle rule's infinite-$`N`$
limit at the return turning point $`t = 16`$, on a colour scale of
$`\pm0.1`$ (three times finer than a–c), with the QLE's $`\pm0.1`$ contours
where the state should be. (e) Position marginals: the coarse cell lowers the main peak, the
fine cell does not, and the one-particle rule is far off. (f) Proposition X4:
the one-particle run's energy against the prediction from its own
$`\langle M_2\rangle`$.*

![Annihilation bias](https://raw.githubusercontent.com/billpage/wpmw/output/figures/event_branching_bias.png)

*Figure 3. (a) Error against $`N_0`$. Without annihilation it follows
$`N_0^{-1/2}`$; with annihilation it plateaus at the pooled bias (dotted).
(b) The bias against $`h_p`$ (circles), with the $`h_x`$ variants (crosses).
(c) Positon and negaton populations in the $`h_p = 0.125`$ run.*

---

## 7. The clock

On the call, Cyganski described drawing an exponential waiting time and then
flying the particle, which is a rate frozen at the start of the flight. Part L
implements that literally and compares it with thinning. Both runs are
without annihilation, to $`t = 4.5`$, rms over 4 seeds:

| $`N_0`$ | thinning: $`L^2`$ | $`e\sqrt{N_0}`$ | frozen: $`L^2`$ | $`e\sqrt{N_0}`$ |
|---|---|---|---|---|
| 25 000 | 0.0827 | 13.1 | 0.1024 | 16.2 |
| 100 000 | 0.0446 | 14.1 | 0.0610 | 19.3 |
| 400 000 | 0.0207 | 13.1 | 0.0480 | 30.3 |

Thinning is pure noise. The frozen clock has a bias of about
$`\sqrt{0.0480^2 - 0.0207^2}\approx 0.043`$. That is as large as the bias of
the coarsest annihilation cell, and it makes about 6% too many events (10.08
against 9.51 per initial particle, children included).

The cause is simple. A particle moves about one wavelength of the sine in a
mean flight, so the rate at the start of the flight says little about the
rate where the kick lands. A comparison made with annihilation switched on can hide this bias under
the annihilation bias of §6. That is why the test is run without
annihilation. Thinning costs the same up to its
acceptance ratio (0.645 in the one-particle run of Part K), so there is no
reason to use anything else.
Which clock the notebook uses is not known.

---

## 8. Summary

1. The quantum part of the Wigner equation cannot be driven by a positive
   one-particle clock: it always adds momentum variance, $`M_2 > 0`$, and
   $`M_3^2\le M_2M_4`$ means it cannot carry a quantum correction without
   doing so (X1). That diffusion is heating, at exactly
   $`\langle M_2\rangle/2m`$ (X4), and it is decoherence, damping legs
   $`y`$ apart at rate $`\lambda(1-\cos(ay/\hbar))`$ (X6). The signed pair
   is a pure phase and does neither. The proposed general rule has
   $`M_2 = \hbar\lvert V''\rvert`$, so it fails for the harmonic oscillator
   and for the sinusoid it was meant to generalise.
2. A signed two-point kernel with a locally determined step is exact only for
   $`V''' = rV'`$ with $`r`$ constant: a sinusoid, a quadratic, or an
   exponential with imaginary step (X2). The general case needs the sum of
   per-mode kernels, which is the compensated split.
3. The classical force is the zero-step, infinite-rate member of the signed
   family (X3). Removing it analytically is what the compensated split does,
   and what Cyganski's 2020 code did for its special case.
4. Event-driven pair branching with a thinning clock is unbiased. With curved
   trajectories, cell annihilation introduces a bias set by $`h_p`$. It is a
   decoherence of range $`\hbar/h_p`$, as Theorem D6 predicts (X5).

---

## 9. Open items

- **X-SP1 (annihilation with curved trajectories).** Is there an annihilation
  rule for continuous momenta whose bias falls faster than $`h_p^{1.2}`$
  without the population cost of refining the cell? Candidates are
  cancellation only between particles closer than $`h_p`$ in both variables,
  and pairing along the flow. The measured exponent, between the linear and
  quadratic behaviour one might predict from Theorem D6, is not explained.
- **X-SP2 (general potentials).** Extend the event-driven code to a
  reach-limited residual kernel for a non-single-mode potential, with a
  quartic or Eckart benchmark. Theorem X2 says the local rule cannot do this.
- **X-SP3 (the notebook's overshoot).** This needs Cyganski's code, or the
  $`N_0`$ and cell-size study of §6 run in it.
- **X-SP4 (Kapitza–Dirac).** The Skellam distribution (unsigned
  $`\pm\hbar k`$ kicks, $`e^{-2\mu}I_n(2\mu)`$) against the Bessel
  distribution (signed kicks, $`J_n(2\mu)`$) is a compact experimental
  illustration of why the sign is indispensable. It was proposed on the call
  and is not yet done.

---

## 10. What was verified, and how

| Statement | Method | Where |
|---|---|---|
| (1.1) | SymPy, symbolic Taylor coefficients to $`n=7`$, residual 0 | Part A |
| X1 moments | SymPy, all four sign cases | Part B |
| X1 parts 1–3 | Cauchy–Schwarz and positivity (proof in §2.2) | — |
| X6 symbols, series, $`M_3^2 = M_2M_4`$ for one kick size | SymPy | Part B2 |
| Purity and negative mass, QLE against one-particle | spectral | Part G |
| (3.1), $`M_5`$ and $`M_7`$ residuals, Theorem X2 | SymPy | Part C |
| Quadratic + sine: $`a^2(x)`$, $`\delta M_5/M_5`$, split exact to $`M_7`$ | SymPy; fractions numerical | Part D |
| X3 symbol expansion | SymPy | Part E |
| X3 error $`\propto a^2`$ | spectral, $`256^2`$ grid, exponent 2.006 | Part E |
| References agree | Schrödinger against spectral QLE, $`10^{-5}`$ | Part F |
| X4 | spectral master equation ($`7\times10^{-4}`$) and Monte Carlo | Parts G, K |
| Branching unbiased | Monte Carlo, $`e\sqrt{N_0}`$ constant | Part H |
| X5 | Monte Carlo, 3 cells × 4 $`N_0`$, plus $`h_x`$, sync and centroid variants | Parts I, J |
| Frozen clock | Monte Carlo | Part L |

---

## References

- T. Takabayasi, *The Formulation of Quantum Mechanics in terms of Ensemble
  in Phase Space*, Prog. Theor. Phys. **11**, 341 (1954) — (3.17).
- R. F. Pawula, *Approximation of the linear Boltzmann equation by the
  Fokker–Planck equation*, Phys. Rev. **162**, 186 (1967).
- P. A. W. Lewis and G. S. Shedler, *Simulation of nonhomogeneous Poisson
  processes by thinning*, Naval Res. Logist. Q. **26**, 403 (1979).
- D. T. Gillespie, *A general method for numerically simulating the
  stochastic time evolution of coupled chemical reactions*, J. Comput. Phys.
  **22**, 403 (1976).
- J. M. Sellier, M. Nedjalkov, I. Dimov, *An introduction to applied quantum
  mechanics in the Wigner Monte Carlo formalism*, Phys. Rep. **577**, 1
  (2015) — cell annihilation in signed-particle Monte Carlo.
- D. Cyganski, *Extended Fokker–Planck Eq. and the QLE V2* (project memo) —
  moment matching, signed jump densities, imaginary jumps.
