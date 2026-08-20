# Compensated Liouville splitting: putting the classical force in the first substep

**Status.** Analysis note. Companion demo:
`src/demo_compensated_liouville_splitting.py`. Prompted by a question from
Bill Page (August 2026): the algorithm of
[`../supplement/phase_space_crystal_lattice_supplement.md`](../supplement/phase_space_crystal_lattice_supplement.md)
§7 free-streams first and then applies mediated jumps for *every* Fourier
mode, even when the dynamics is classical. Is there a formal way to put the
zeroth-order classical term into the first substep and leave only the higher
modes contributing jumps?

The answer is yes, the reorganisation is exact rather than approximate, and
the interesting part is what it does and does not buy.

---

## 0. What this note inherits, corrects and contradicts

**Inherits.** From
[`../supplement/phase_space_crystal_lattice_supplement.md`](../supplement/phase_space_crystal_lattice_supplement.md)
§5–§7: the momentum quantum $`\Delta p = \pi\hbar/L`$, the mode $`q`$ stencil,
the rate field $`\Gamma_q(x)`$, and the two-substep algorithm. From
[`open_position_space.md`](open_position_space.md) §3: the coherence horizon
$`L_c`$ on the ket–bra separation, which reappears here in §4.3 as exactly the
window that removes the seam artefact. From
[`four_rule_microdynamics_equivalence.md`](four_rule_microdynamics_equivalence.md):
the four channels and the fact that the mode stencil is their generator.

**Corrects one statement of the supplement.**

`../supplement/phase_space_crystal_lattice_supplement.md` §6.1 derives the
mode stencil by "approximating $`\partial_p W`$ by a centred finite difference
at the photon-momentum scale". The word *approximating* understates the
result. For a single Fourier mode the centred difference at spacing
$`\hbar k/2`$ is an **identity**, not an approximation (Lemma C0 below), and
it is precisely that exactness which lets the mode sum reproduce the whole
Moyal series rather than only its leading term. The formula in the box of
§6.1 is right; only its justification is weaker than it needs to be. An
errata line in that section is warranted.

**Contradicts three points in the framing of the prompting question.**

1. *"It is a theorem that these modes cancel in such a way that for up to
   quadratic potentials all but the zeroth term cancels."* There is no
   inter-mode cancellation. Each mode's stencil is separately exact and
   separately splits into a classical piece plus an $`O(k^3)`$ remainder
   (§2). For a genuinely quadratic $`V`$ the Fourier weight sits at
   $`k = 0`$, where the $`k^3`$ factor annihilates it pointwise. Nothing
   needs to cancel against anything else.
2. *"For example, for the harmonic oscillator potential is represented as an
   infinite number of modes but we know that the evolution is in fact
   classical Liouvillian."* Not on the periodic box. The periodised parabola
   has a kink at the seam, so $`V''' = -m\omega^2 L\thinspace\delta'(x \mp L/2)
   \neq 0`$, and its evolution is not Liouvillian. Theorem C4 sharpens this:
   **on a circle, only constant potentials are exactly classical.** The
   harmonic oscillator is classical on the open line, which is one more
   reason the compensated scheme belongs with
   [`open_position_space.md`](open_position_space.md) rather than with the
   ring.
3. *"Quantum mechanics enters the Wigner equation in two different ways.
   First ... the initial negativity of the Wigner function ... Second ...
   mediated jumps."* Negativity is neither necessary nor the sharp
   criterion. A coherent or squeezed state has $`W \ge 0`$ everywhere and
   evolves by *exact* classical Liouville flow in a harmonic well, and is
   fully quantum. What is quantum about the initial datum is the
   $`\hbar`$-sized floor on phase-space area — admissibility,
   $`|W| \le 2/h`$, Hudson's theorem — not the sign. And there is a third
   channel: the Weyl correspondence on readout. Even with $`W \ge 0`$ and
   exactly classical evolution, $`\langle\hat A\rangle = \int A_W W`$ needs
   the Weyl symbol, which differs from the classical function by
   $`O(\hbar^2)`$ for nonlinear observables. Three channels, then:
   admissible initial data, non-quadratic dynamics, symbol ordering at
   readout. This note is about the second.

**Corrects one statement made in the conversation that prompted this note.**
It was suggested there that the compensated form trades Cyganski's
field-quantum account of the classical force (supplement §5) for a cleaner
generator, and separately that compensation converts a per-particle
stochastic channel into a per-column deterministic one. Theorem C5 shows both
claims are too strong. Operator splitting is **rate neutral**: the jump
intensity of the compensated generator equals that of the uncompensated one,
mode by mode. So the compensated form is a change in the representation of the
*equation*, not of the microdynamics, and the photon-quantum ontology of §5
survives untouched. Event counts change only under *potential* splitting
(§5.3), and only then does the ontological question actually arise.

---

## 1. Tutorial: the potential substep is a multiplication operator

Let $`s`$ be conjugate to momentum,

```math
\hat W(x, s) \;=\; \int dp\thinspace W(x, p)\thinspace e^{-i p s},
\qquad
W(x, p \pm a) \;\longleftrightarrow\; e^{\pm i a s}\thinspace\hat W(x, s).
```

This is the transform the algorithm already performs: it is the "FFT of the
imaginary symbol" by which the momentum-jump substep is integrated exactly.
The potential term of the QLE in mode form is

```math
\partial_t W\big|_V \;=\; \frac{1}{i\hbar}\int\frac{dk}{2\pi}\thinspace
  \tilde V(k)\thinspace e^{ikx}
  \Bigl[W\bigl(p - \tfrac{\hbar k}{2}\bigr)
      - W\bigl(p + \tfrac{\hbar k}{2}\bigr)\Bigr],
```

and under the shift theorem the bracket becomes
$`e^{-i\hbar k s/2} - e^{i\hbar k s/2} = -2i\sin(\hbar k s/2)`$, so the whole
term collapses to multiplication by a single scalar function:

```math
\boxed{\;
M(x, s) \;=\; \frac{i}{\hbar}
\Bigl[V\bigl(x + \tfrac{\hbar s}{2}\bigr)
    - V\bigl(x - \tfrac{\hbar s}{2}\bigr)\Bigr].
\;}
```

Three things follow immediately and are worth stating because each is a
familiar fact of the project seen from a new angle.

- **Only odd derivatives survive.** $`M`$ depends on $`V`$ only through the
  odd-in-$`\xi`$ combination $`U(x,\xi) = V(x+\xi/2) - V(x-\xi/2)`$ evaluated
  at $`\xi = \hbar s`$. Expanding, $`U = V'\xi + V'''\xi^3/24 + \cdots`$: the
  even part of $`V`$ about $`x`$ cancels identically. The Moyal series is the
  Taylor expansion of $`M`$ in $`s`$, term by term.
- **$`M`$ is purely imaginary**, so $`|e^{\tau M}| = 1`$ and the substep is
  norm preserving on each $`(x, s)`$ line — the reason the exact-FFT
  integration of the momentum substep is a hard requirement and explicit
  Euler is not.
- **The stencil arm is $`\hbar s/2`$**, so the nonlocality of the potential
  term is a *reach* in position, not in momentum. At the crystal quantum
  $`\Delta p = \pi\hbar/L`$ the Nyquist value of $`s`$ is $`\pi/\Delta p = L/\hbar`$,
  hence a maximum arm of exactly $`L/2`$. The momentum quantum and the box
  length are the same statement twice: the stencil samples the potential over
  exactly one box length and no more.

**Lemma C0 (the finite difference is exact).** For
$`V(x) = V_q\cos\theta`$, $`\theta = 2\pi q x/L + \phi`$, with
$`\Delta p_q = q\pi\hbar/L = \hbar k/2`$,

```math
\partial_t W\big|_V \;=\; -\thinspace\frac{V_q}{\hbar}\thinspace\sin\theta
  \thinspace\bigl[W(p + \Delta p_q) - W(p - \Delta p_q)\bigr]
```

*exactly*, to all orders in $`\hbar`$.

*Proof.* Put $`\tilde V(k) = \pi V_q(e^{i\phi}\delta(k - k_q)
+ e^{-i\phi}\delta(k + k_q))`$ into the mode form above; the two delta
contributions combine to $`(2i\sin\theta)`$ times the momentum difference,
and the prefactor $`1/(i\hbar)`$ gives the stated result. Equivalently
$`M(x,s) = (2iV_q/\hbar)\sin\theta\thinspace\sin(\hbar k s/2)`$, whose inverse
transform is the two-point stencil. $`\square`$

This is the corrected reading of supplement §6.1 announced in §0.

---

## 2. The split

### 2.1 Statement

Subtract from $`M`$ the part linear in the stencil arm:

```math
M \;=\; \underbrace{i\thinspace V'(x)\thinspace s}_{M_{\mathrm{cl}}}
\;+\;
\underbrace{\frac{i}{\hbar}\Bigl[
  V\bigl(x + \tfrac{\hbar s}{2}\bigr) - V\bigl(x - \tfrac{\hbar s}{2}\bigr)
  - \hbar s\thinspace V'(x)\Bigr]}_{M_{\mathrm{res}}}.
```

$`M_{\mathrm{cl}}`$ is exactly the classical Liouville force term: under
$`\partial_p \leftrightarrow is`$, multiplication by $`iV's`$ is
$`V'(x)\thinspace\partial_p`$, which is $`\dot p = -V'`$ fed into the Liouville
flow.

**Theorem C1 (exact commuting factorisation).** $`M_{\mathrm{cl}}`$ and
$`M_{\mathrm{res}}`$ are both multiplication operators in the same pair of
variables $`(x, s)`$. Hence they commute exactly and

```math
e^{\tau M} \;=\; e^{\tau M_{\mathrm{cl}}}\thinspace e^{\tau M_{\mathrm{res}}}
\qquad\text{for every } \tau,
```

with **no Trotter error**.

*Verification.* Part A of the demo: the commutator is $`0`$ identically in
floating point, and the factorisation error is
$`1.1\times10^{-16}`$, $`2.6\times10^{-16}`$, $`9.5\times10^{-16}`$ at
$`\tau = 0.01, 0.1, 1.0`$ — i.e. rounding only, with no growth in $`\tau`$.

The caveat worth stating plainly: this exactness is *within* the potential
substep. The split between free streaming and the potential is unaffected and
retains its usual Strang error, because advection is diagonal in $`(k_x, p)`$
and the potential term is diagonal in $`(x, s)`$.

### 2.2 What sources the residual

**Theorem C2 (the residual is a third-derivative object).** With
$`a = \hbar s/2`$,

```math
M_{\mathrm{res}}(x, s) \;=\; \frac{i}{\hbar}\int_0^{a}
  \bigl[V'(x + \sigma) + V'(x - \sigma) - 2V'(x)\bigr]\thinspace d\sigma ,
```

a running second central difference of $`V'`$. Consequently
$`M_{\mathrm{res}} \equiv 0`$ for all $`s`$ if and only if
$`V''' \equiv 0`$ on the region reached by the stencil arm.

*Proof.* $`\int_0^{a}V'(x+\sigma)d\sigma = V(x+a) - V(x)`$ and
$`\int_0^{a}V'(x-\sigma)d\sigma = V(x) - V(x-a)`$; their sum is
$`V(x+a) - V(x-a)`$, and $`\int_0^{a}2V'(x)d\sigma = \hbar s V'(x)`$.
The bracket is $`\sigma^2 V'''(x) + O(\sigma^4)`$, vanishing identically iff
$`V''' \equiv 0`$. $`\square`$

*Verification.* Part A: for the parabola on the open line,
$`\max|M_{\mathrm{res}}| = 1.4\times10^{-14}`$ against
$`\max|M| = 1.3\times10^{2}`$.

Per Fourier mode the split reads

```math
M_q \;=\; \frac{2iV_q}{\hbar}\sin\theta\thinspace\sin u,
\qquad
M_{\mathrm{cl},q} \;=\; \frac{2iV_q}{\hbar}\sin\theta\thinspace u,
\qquad
u \;\equiv\; \frac{\hbar k s}{2} \;=\; q\thinspace\Delta p\thinspace s,
```

so the residual of a mode is $`\sin u - u`$ against the classical $`u`$: the
whole content of the reorganisation is the difference between a sine and its
tangent at the origin.

---

## 3. The band condition

**Theorem C3 (where compensation helps).** Per mode,

```math
\frac{|M_{\mathrm{res},q}|}{|M_{\mathrm{cl},q}|}
\;=\; \frac{|\sin u - u|}{|u|}
\;=\; \frac{u^2}{6} + O(u^4),
\qquad
\frac{|M_{\mathrm{res},q}|}{|M_q|}
\;=\; \frac{|\sin u - u|}{|\sin u|} .
```

The first ratio is $`1`$ at $`u = \pi`$. So the reorganisation is a gain only
inside the band $`|u| \ll \pi`$, and at the Nyquist edge of the momentum grid
it is no gain at all.

| $`u/\pi`$ | 0.02 | 0.05 | 0.10 | 0.25 | 0.50 | 0.75 | 1.00 |
|---|---|---|---|---|---|---|---|
| $`\lvert M_{\mathrm{res}}\rvert/\lvert M_{\mathrm{cl}}\rvert`$ | 0.0007 | 0.0041 | 0.0164 | 0.0997 | 0.3634 | 0.6999 | 1.0000 |
| leading $`u^2/6`$ | 0.0007 | 0.0041 | 0.0165 | 0.1028 | 0.4112 | 0.9253 | 1.6449 |

In state terms: a state of momentum width $`\sigma_p`$ populates
$`|s| \lesssim 1/\sigma_p`$, so $`u \lesssim \hbar k/(2\sigma_p)`$. The band
condition is

```math
\sigma_p \;\gg\; \frac{\hbar k}{2}
\qquad\text{equivalently}\qquad
q \;\ll\; \frac{\sigma_p}{\Delta p} \;=\; \frac{\sigma_p L}{\pi\hbar}.
```

**The state's momentum width must exceed the kick it receives.** This is the
semiclassical condition wearing different clothes, and it is the single most
important practical fact in this note: compensation is worth something exactly
when the dynamics is nearly classical, and worth nothing when it is not.

![Symbol decomposition and band condition](https://raw.githubusercontent.com/billpage/wpmw/output/figures/compensated_symbol_and_band.png)

---

## 4. The ring no-go

### 4.1 Only constants are classical on a circle

**Theorem C4.** Let $`V`$ be a smooth potential on the circle of
circumference $`L`$. Then $`M_{\mathrm{res}} \equiv 0`$ if and only if $`V`$
is constant.

*Proof.* By C2, $`M_{\mathrm{res}} \equiv 0 \iff V''' \equiv 0`$, so
$`V(x) = a + bx + cx^2`$ locally. Single-valuedness on the circle forces
$`b = c = 0`$. $`\square`$

So the periodic box, which was introduced to fix $`\Delta p`$ and keep worlds
in view, has the side effect of making *every* non-trivial potential
non-classical. The harmonic oscillator is not an exception; it is the
sharpest example.

### 4.2 The bowtie

For the periodised parabola $`V(x) = \tfrac{1}{2}m\omega^2(\mathrm{wrap}\thinspace x)^2`$
the failure is not spread out — it is geometrically exact.

**Proposition C4.1.** $`M_{\mathrm{res}}(x, s) = 0`$ exactly if and only if
$`|x| + |\hbar s/2| < L/2`$.

*Proof.* If both arms $`x \pm \hbar s/2`$ remain inside the fundamental
domain the wrap is inactive, $`V(x+a) - V(x-a) = 2m\omega^2 a x = \hbar s V'(x)`$,
and the residual vanishes. Otherwise one arm crosses the seam and picks up
the jump in $`V'`$. $`\square`$

Numerically: the residual is $`1.1\times10^{-14}`$ throughout the diamond and
occupies $`0.4961`$ of the $`(x, s)`$ grid against a predicted $`0.5039`$ —
the small deficit is the half-cell at the grid edge.

![Residual support is a bowtie](https://raw.githubusercontent.com/billpage/wpmw/output/figures/compensated_ring_residual.png)

Note the coincidence with §1: at $`\Delta p = \pi\hbar/L`$ the maximum arm is
exactly $`L/2`$, so the diamond is inscribed in the computational domain and
the seam is reachable only from the extreme corners of the $`s`$-grid. The
momentum quantum is precisely the value at which the ring first notices its
own seam.

### 4.3 The coherence horizon removes it

The coherence horizon of [`open_position_space.md`](open_position_space.md)
§3 truncates the ket–bra separation at $`L_c`$, which in the present
variables is the window $`|\hbar s/2| \le L_c/2`$. Combining with C4.1:

**Corollary C4.2.** Under a coherence horizon $`L_c`$, the periodised
parabola is *exactly* Liouvillian for every $`|x| < L/2 - L_c/2`$. For a
state supported in $`|x| < x_m`$ it suffices that $`L_c < L - 2x_m`$.

*Verification.* Part C: $`\max|M_{\mathrm{res}}|`$ over the horizon-windowed
region is $`3.6\times10^{-15}`$, $`3.6\times10^{-15}`$,
$`1.8\times10^{-15}`$ for $`L_c = 2, 4, 6`$.

This is the cleanest justification the coherence horizon has been given so
far: it is not only a way to bound the jump budget for a localised scatterer,
it is exactly the truncation that deletes the ring's seam artefact and
restores the harmonic oscillator to its classical status.

---

## 5. What the reorganisation costs

### 5.1 Operator splitting is rate neutral

This is the result that most constrains how the reorganisation may be used.

**Theorem C5.** For each mode $`q`$, the generator $`M_q`$ has jump intensity
$`2|\Gamma_q(x)|`$ (the two hop channels). The compensated generator
$`M_{\mathrm{res},q} = M_q - M_{\mathrm{cl},q}`$ is the *same* hop pair minus
a drift generator, and a drift generator contributes no jump intensity.
Hence $`M_{\mathrm{res},q}`$ and $`M_q`$ have identical jump intensities, and
the process

> deterministic drift by $`-V'(x)\thinspace\Delta t`$, then jumps from
> $`M_{\mathrm{res}}`$

is the *same stochastic process* as the uncompensated four-rule rule, event
for event.

*Proof.* In Lévy–Itô form the compensator of a jump measure is a drift, so
$`b\thinspace dt + \int z\thinspace\tilde N(dz, dt)`$ with
$`b = \int z\thinspace\nu(dz)`$ is identically
$`\int z\thinspace N(dz, dt)`$. The compensation is bookkeeping. $`\square`$

| $`q`$ | $`\lvert V_q\rvert`$ | hop intensity of $`M_q`$ | hop intensity of $`M_{\mathrm{res},q}`$ |
|---|---|---|---|
| 1 | 12.96911 | 25.93822 | 25.93822 |
| 2 | 3.24228 | 6.48456 | 6.48456 |
| 3 | 1.44101 | 2.88202 | 2.88202 |
| 4 | 0.81057 | 1.62114 | 1.62114 |

Two consequences.

- **Do not expect a cheaper microdynamics from C1 alone.** In particular the
  integer-count lattice of supplement §7 cannot even represent the
  deterministic drift, since $`p`$ is quantised; realising the drift as a
  lattice hop *adds* a second stencil and makes the budget strictly worse, by
  a factor $`(1 + q)`$ for mode $`q`$.
- **The ontology of supplement §5 is untouched.** Since the process is
  unchanged, the photon-quantum account of the kicks stands exactly as
  written. The compensated form is a statement about the equation and about
  mesh or spectral integrators, not a rival microdynamics.

### 5.2 Where the reorganisation does pay

Three places, in increasing order of value.

1. **Mesh and spectral integrators.** $`e^{\tau M_{\mathrm{res}}}`$ is
   computed from $`V`$ pointwise; no Fourier decomposition of the potential is
   needed at all, in contrast to the mode-by-mode stencil sum. The classical
   channel becomes an exact symplectic flow that can absorb the free-streaming
   step, removing the Strang error between advection and force.
2. **Diagnosis.** The residual is a direct, local measure of how
   non-classical a configuration is; §4.2 is the example.
3. **Potential splitting** (§5.3), which does change event counts.

### 5.3 Potential splitting is not rate neutral

Write $`V = Q + R`$ with $`Q`$ *globally quadratic*. Then $`Q`$ contributes
only $`M_{\mathrm{cl}}`$, integrable in closed form by exact characteristics,
and $`R`$ alone carries jump channels. For the periodised parabola with $`Q`$
the true parabola, $`R`$ is supported at the seam and the interior hop budget
is exactly zero, against a total budget of $`41.6`$ for the untruncated
mode sum. This — not C1 — is how jump events actually disappear for the
harmonic oscillator.

The price: $`Q`$ is not periodic, so this route requires open position space.

### 5.4 Truncating the mode sum

The two schemes make *different* errors when the mode sum is truncated at
$`Q`$ modes. Uncompensated, the dropped piece is $`\sum_{q>Q}M_q`$ and the
classical force itself is corrupted — for the periodised parabola $`V'`$ is a
sawtooth with coefficients $`\sim 1/q`$, so this converges slowly.
Compensated, the force is supplied analytically and exactly, and the dropped
piece is $`\sum_{q>Q}M_{\mathrm{res},q}`$.

Per mode the ratio of the two errors is $`|\sin u - u|/|\sin u|`$, which is
$`\sim u^2/6`$ in band and $`\sim u`$ out of band. **Compensation helps for
modes in band and hurts for modes out of band, by exactly the factor C3
predicts.** Both effects are visible:

| modes $`Q`$ | 1 | 4 | 12 | 24 | 48 |
|---|---|---|---|---|---|
| worst case over the grid, uncompensated | 57.94 | 60.82 | 62.92 | 63.46 | 63.73 |
| worst case over the grid, compensated | 155.90 | 228.13 | 258.95 | 268.68 | 270.75 |
| $`L^1`$ of one substep, uncompensated | 7.19e-2 | 4.50e-2 | 5.11e-2 | 5.05e-2 | 5.03e-2 |
| $`L^1`$ of one substep, compensated | 2.30e-3 | 7.59e-3 | 5.12e-3 | 2.42e-3 | 1.15e-3 |

The worst-case row is not a bug in the compensated scheme; it is the seam
again. The dropped residual tail is $`-\sum_{q>Q}M_{\mathrm{cl},q}`$, which
lives at the seam and at large $`|s|`$. A state that stays away from the seam
never sees it, and for the Gaussian used here ($`\sigma_x = 0.5`$ at the well
bottom) the compensated scheme is between 6 and 44 times more accurate. A
guarantee that must hold uniformly over the whole domain gets the opposite
verdict. Both measures are correct; which one is decision-relevant depends on
whether one is running a simulation or proving a bound.

![Mode truncation](https://raw.githubusercontent.com/billpage/wpmw/output/figures/compensated_truncation.png)

---

## 6. The two test problems side by side

A Gaussian of width $`\sigma_x = 0.4`$ released at $`x = 1`$, integrated by
Strang splitting to $`t = 2`$ with $`\Delta t = 0.01`$, under the exact
symbol and under the classical term alone. The gap between the two runs is
the entire quantum content of the evolution.

| | periodised parabola | cosine well |
|---|---|---|
| $`L^1`$ gap at $`t = 2`$ | $`8.28\times10^{-5}`$ | $`7.44\times10^{-1}`$ |
| negativity generated, exact run | $`2.36\times10^{-5}`$ | — |
| negativity generated, classical run | $`4.38\times10^{-9}`$ | — |

The periodised parabola is nearly nine thousand times more classical than the
cosine well of comparable depth, and its residual negativity — five orders of
magnitude above the classical run's numerical floor, four orders below the
cosine well — is entirely seam-generated, reaching the packet at $`x = 1`$
only through the far tail of the stencil arm. This is the quantitative
version of the claim in §0 that the ring parabola is *nearly* but not exactly
classical.

![Evolution comparison](https://raw.githubusercontent.com/billpage/wpmw/output/figures/compensated_evolution.png)

---

## 7. Coulomb

The Coulomb potential is the case where the compensated form is most
informative, because everything can be done in closed form.

### 7.1 The exact ratio

**Theorem C6.** For $`V(x) = -Z/x`$ and $`a = \hbar s/2`$, write
$`\rho = a/x`$. Then

```math
M \;=\; \frac{2iZa}{\hbar\thinspace(x^2 - a^2)},
\qquad
M_{\mathrm{cl}} \;=\; \frac{2iZa}{\hbar\thinspace x^2},
\qquad
\frac{M_{\mathrm{res}}}{M_{\mathrm{cl}}} \;=\; \frac{\rho^2}{1 - \rho^2}.
```

*Proof.* Direct: $`-Z/(x+a) + Z/(x-a) = 2Za/(x^2-a^2)`$, and
$`V' = Z/x^2`$. $`\square`$

*Verification.* Part G: agreement to $`7.1\times10^{-15}`$ over a grid of
radii and $`\rho`$ values.

Four things fall out.

**(i) The Moyal series for Coulomb is geometric, and its radius of
convergence is the distance to the nucleus.** Since
$`M/M_{\mathrm{cl}} = (1-\rho^2)^{-1} = \sum_{n\ge0}\rho^{2n}`$, the $`n`$-th
Moyal term is exactly $`\rho^{2n}`$ times the classical one, all of the same
sign. The series converges if and only if $`|\rho| < 1`$, that is

```math
\Bigl|\frac{\hbar s}{2}\Bigr| \;<\; |x| .
```

**The stencil arm must not reach the nucleus.** There is no cancellation
among Moyal orders to rescue a configuration where it does; the terms all add.

**(ii) Coulomb is scale free, so the band condition is the same at every
radius.** The ratio depends only on $`\rho`$, which involves $`x`$ and $`s`$
only through their product with $`\hbar`$. In state terms $`s \sim 1/\sigma_p`$
gives $`\rho \sim \hbar/(2\sigma_p|x|)`$, so C3 becomes

```math
\sigma_p\thinspace|x| \;\gg\; \frac{\hbar}{2},
```

the uncertainty product in units of $`\hbar`$. This is the semiclassical
condition in its most naked form.

| $`\sigma_p\lvert x\rvert/\hbar`$ | 0.6 | 1.0 | 2.0 | 5.0 | 10.0 | 20.0 |
|---|---|---|---|---|---|---|
| residual fraction | 2.27 | 0.333 | 0.0667 | 0.0101 | 0.0025 | 0.0006 |

For the hydrogen ground state $`\sigma_p a_0 \sim \hbar`$, so $`\rho \sim 1/2`$
and compensation buys about a factor of three — real but modest. For a
Rydberg state of principal quantum number $`n`$, $`\sigma_p r \sim n\hbar`$,
$`\rho \sim 1/(2n)`$, and the gain is $`\sim 4n^2`$. Compensation is a
Rydberg and scattering tool, not a ground-state tool.

**(iii) After compensation the jump channel is a core effect.** To leading
order $`|M_{\mathrm{res}}| \simeq 2Za^3/(\hbar x^4)`$ against
$`|M_{\mathrm{cl}}| = 2Za/(\hbar x^2)`$: the classical force falls off as
$`1/x^2`$ and the quantum residual as $`1/x^4`$. Uncompensated, the four-rule
scheme fires jumps wherever the potential has Fourier content, which for
Coulomb is everywhere; compensated, the entire long-range Kepler part of the
dynamics goes to the deterministic flow and the stochastic channel is
concentrated where it belongs.

| $`x`$ | 0.25 | 0.5 | 1.0 | 2.0 | 4.0 | 8.0 |
|---|---|---|---|---|---|---|
| $`\lvert M_{\mathrm{res}}\rvert`$ exact, $`a = 0.1`$ | 6.10e-1 | 3.33e-2 | 2.02e-3 | 1.25e-4 | 7.82e-6 | 4.88e-7 |
| $`2Za^3/(\hbar x^4)`$ | 5.12e-1 | 3.20e-2 | 2.00e-3 | 1.25e-4 | 7.81e-6 | 4.88e-7 |

**(iv) Compensation localises but does not regularise.** The symbol is
singular on the cone $`a = x`$; the residual inherits the singularity, since
subtracting a finite $`M_{\mathrm{cl}}`$ cannot remove it. A soft core
$`V = -Z/\sqrt{x^2 + \epsilon^2}`$ tames it, at a price set directly by
$`\epsilon`$:

| $`\epsilon`$ | 0 | 0.05 | 0.2 | 0.5 |
|---|---|---|---|---|
| $`\max\lvert M_{\mathrm{res}}\rvert`$ at $`x = 1`$ | 5.98e2 | 1.75e1 | 2.63 | 1.12 |

![Coulomb](https://raw.githubusercontent.com/billpage/wpmw/output/figures/compensated_coulomb.png)

### 7.2 A warning about the mode representation

In three dimensions $`\tilde V(k) = -4\pi Z/k^2`$, so the residual weight per
mode goes as $`\tilde V(k)k^3 \sim k`$: it *grows* with mode number. The
Fourier-mode representation on which the crystal-lattice algorithm is built is
badly conditioned for Coulomb whether or not one compensates, and the
compensated form does not repair that. The right move for Coulomb is to work
with $`M(x, s)`$ directly — which the compensated form permits, since
$`M_{\mathrm{res}}`$ is defined pointwise from $`V`$ without any mode
decomposition — and to accept that a mediated-jump particle reading of the
Coulomb problem needs either screening, a soft core, or a coherence horizon.

---

## 8. Summary

| | statement | verified |
|---|---|---|
| C0 | the mode finite difference is exact, not approximate | analytic |
| C1 | $`e^{\tau M} = e^{\tau M_{\mathrm{cl}}}e^{\tau M_{\mathrm{res}}}`$, no Trotter error | $`\sim10^{-16}`$ |
| C2 | $`M_{\mathrm{res}} \equiv 0 \iff V''' \equiv 0`$ | $`1.4\times10^{-14}`$ |
| C3 | gain requires $`\lvert u\rvert \ll \pi`$, i.e. $`\sigma_p \gg \hbar k/2`$ | table, §3 |
| C4 | on a circle only constants are classical; support is the bowtie | $`1.1\times10^{-14}`$ |
| C4.2 | a coherence horizon $`L_c < L - 2x_m`$ restores exact classicality | $`3.6\times10^{-15}`$ |
| C5 | operator splitting is rate neutral; potential splitting is not | table, §5.1 |
| C6 | Coulomb ratio $`\rho^2/(1-\rho^2)`$; core-localised residual | $`7.1\times10^{-15}`$ |

The short version. The question asked for a formal way to put the classical
force into the first substep and leave jumps to the higher modes. That
reorganisation exists, is exact, and costs nothing in accuracy — but it is a
reorganisation of the *equation*. At the level of the particle rule it is the
identity (C5). To actually remove jump events one must split the *potential*,
and that requires a globally quadratic piece, which requires open position
space.

---

## 9. Open items

- **CLS1.** Give $`R = V - Q`$ of §5.3 a four-rule reading. On the ring $`R`$
  is the seam remainder and is not smooth; on the open line $`R`$ has no
  Fourier series in the crystal sense. What is the stencil geometry for a
  potential presented pointwise rather than by modes?
- **CLS2 (conjecture).** C5 says compensation does not change the process, so
  it does not by itself contradict the photon-quantum reading of supplement
  §5. But if the drift *is* implemented deterministically, does it have a
  field reading of its own? The natural candidate is a coherent state of the
  mode field: the mean occupation supplies a classical force, the fluctuation
  about it supplies the residual kicks. If that holds, the compensated split
  is the mean-field / shot-noise decomposition of the same field, and the
  ontology is not merely preserved but explained. This wants checking against
  the four-wave-mixing account.
- **CLS3.** The bowtie of C4.1 is exact for the periodised parabola. What is
  the residual support for a general periodised polynomial, and does it
  always inscribe a region whose size is set by $`\Delta p`$?
- **CLS4.** Corollary C4.2 makes the coherence horizon exactly the seam
  eraser. Does the same window make the *cosine* well cheaper, or is the
  cosine's non-classicality horizon independent? Part F suggests the latter,
  but this was not measured against $`L_c`$.
- **CLS5.** Coulomb in 3D with a soft core: measure the residual budget
  against $`\epsilon`$ and against the coherence horizon, and determine which
  regularisation is cheaper at fixed accuracy.
- **CLS6.** An errata line for
  `../supplement/phase_space_crystal_lattice_supplement.md` §6.1, replacing
  "approximating" by the exactness statement of Lemma C0. Not applied here to
  keep this note self-contained; see §0.

---

## Sources

- David Cyganski, *Extended Fokker–Planck Eq. and the QLE V2* (project memo).
- [`../supplement/phase_space_crystal_lattice_supplement.md`](../supplement/phase_space_crystal_lattice_supplement.md) §5–§7.
- [`open_position_space.md`](open_position_space.md) §3 (the coherence horizon).
- [`four_rule_microdynamics_equivalence.md`](four_rule_microdynamics_equivalence.md).
- R. L. Hudson, "When is the Wigner quasi-probability density non-negative?", *Rep. Math. Phys.* **6**, 249 (1974) — for the claim in §0 that positivity of $`W`$ for a pure state characterises Gaussians.
