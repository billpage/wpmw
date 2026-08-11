# The Inverted Pair Barrier as a Test Case

*Two particles on a line with* $`V_2(r) = -\tfrac12 m\omega^2 (x_1-x_2)^2`$, *treated as a two-body problem.*

## 0. Status and provenance

[`../analysis/fourd_microdynamics.md`](../analysis/fourd_microdynamics.md) §5
established that a *confining* harmonic potential is a null test: the Moyal
bracket truncates, the QLE is exactly classical Liouville, and the
microdynamics reaches that answer only by cancellation among an unbounded
number of noisy channels. §5.7 of that note proposed an inverted parabolic
barrier as the replacement. This note works the proposal out properly for the
case B. Page actually asked about — an inverted harmonic *pair* interaction,
$-(x_2-x_1)^2$, as a genuine two-body problem — and specifies it as a test
case.

The verdict is favourable but not for the reason §5.7 gave. Companion code:
`src/demo_inverted_pair_barrier.py`. Developed with Claude (Anthropic), August
2026.

### 0.1 Corrections recorded here

1. **§5.7 of the 4-D note implied the inverted barrier would be cheaper. It is
   not.** Reversing the sign of the potential leaves every
   $`|V_{q}|`$ unchanged, so the per-mode noise constant is *identical* —
   $`2\mu\Omega^2\hbar = 2m\omega^2\hbar`$, verified to $2.2\times10^{-16}$
   (Part E1). The inverted barrier is a better *probe*, not a cheaper one; it
   buys a non-trivial observable at the same price.
2. **It also has an obstruction the trap does not have.** Scattering states
   escape, and on the periodic crystal lattice both the position extent and
   the momentum extent grow as $`e^{\Omega t}`$. The usable window grows only
   logarithmically in the grid size, and a target accuracy $\epsilon$ in the
   transmission by literal propagation needs $`N \sim \epsilon^{-2}`$ cells
   (Part E3).
3. **The obstruction turns out not to bind, for a reason worth stating
   separately** (§3): the transmission probability is an *exactly conserved
   functional* of the Wigner function, so it can be read at any time,
   including $t = 0$. The escape problem limits literal simulation of the
   asymptotic states; it does not limit the test.

---

## 1. The system

```math
H \;=\; \frac{p_1^2 + p_2^2}{2m} \;-\; \frac{m\omega^2}{2}\bigl(x_1 - x_2\bigr)^2 .
```

The force between the particles is $`F = +m\omega^2 r`$ with
$r = x_1 - x_2$: repulsive, growing with separation. There is a barrier at
$r = 0$ which the pair either crosses or does not.

Separating with $R = (x_1+x_2)/2$ and $r = x_1 - x_2$:

- the centre of mass is **free**, at total mass $2m$;
- the relative coordinate is an **inverted oscillator** at reduced mass
  $\mu = m/2$,

```math
H_{\mathrm{rel}} \;=\; \frac{p_r^2}{2\mu} \;-\; \frac{\mu}{2}\thinspace\Omega^2 r^2,
\qquad
\mu\thinspace\Omega^2 = m\omega^2 \;\Longrightarrow\; \Omega = \sqrt{2}\thinspace\omega .
```

$\Omega$ is the Lyapunov exponent of the relative motion. Part A confirms the
joint classical flow matrix has spectrum
$`\lbrace -\Omega, 0, 0, +\Omega\rbrace`$ to $4.4\times10^{-16}$; the two zeros
are the free centre of mass.

Since $H$ is quadratic, all third derivatives vanish and the Moyal bracket
truncates: **the QLE is exactly the classical Liouville equation**, as for the
trap. What the inverted case adds is asymptotic free states, and therefore a
scattering observable that is genuinely quantum.

---

## 2. Conventional Wigner treatment

### 2.1 The hyperbolic flow and its eigendirections

The relative flow $`\dot{(r, p_r)} = A_{\mathrm{rel}}(r, p_r)`$ has
eigenvectors

```math
u \;=\; r + \frac{p_r}{\mu\Omega} \quad (\lambda = +\Omega),
\qquad
v \;=\; r - \frac{p_r}{\mu\Omega} \quad (\lambda = -\Omega),
```

so $`u(t) = e^{\Omega t}u(0)`$ and $`v(t) = e^{-\Omega t}v(0)`$. Verified to
$`3.6\times10^{-12}`$ over $t \le 3.5$ — the residual is the matrix
exponential's own rounding, itself growing as $`e^{\Omega t}`$ (Part B).

Two consequences follow immediately.

**The asymptotic side is decided by the initial $u$.** As $t\to\infty$,
$`r(t) \to \tfrac12 e^{\Omega t} u(0)`$, so $`\mathrm{sign}\thinspace r(\infty) = \mathrm{sign}\thinspace u(0)`$.

**The half-plane $u > 0$ is exactly invariant.** Therefore

```math
T \;=\; \iint_{u>0} W(r, p_r)\; dr\, dp_r
```

is a **constant of the motion**, and it equals the asymptotic transmission
probability
$`\lim_{t\to\infty}\int_{r>0}\negthinspace dr\int\negthinspace dp_r\, W`$.
That limit is the exact quantum probability of finding the pair on the far
side, because Wigner evolution is exact here and the position projector is a
legitimate phase-space observable.

![Inverted pair barrier in the relative phase plane](https://raw.githubusercontent.com/billpage/wpmw/output/figures/inverted_pair_barrier_phase_space.png)

### 2.2 A closed-form transmission probability

For a Gaussian relative wavepacket with centre $`(r_c, p_c)`$ and widths
$`(\sigma_r, \sigma_p)`$, $u$ is Gaussian and

```math
\boxed{\;
T \;=\; \Phi\negthinspace\left(\frac{\bar u}{\sigma_u}\right),
\qquad
\bar u = r_c + \frac{p_c}{\mu\Omega},
\qquad
\sigma_u^2 = \sigma_r^2 + \frac{\sigma_p^2}{(\mu\Omega)^2}\; }
```

with $\Phi$ the standard normal distribution function. Part C1 checks this
against an independent split-operator Schrödinger solve of the relative
coordinate; the residual falls as $`e^{-\Omega t}`$ —
$3.9\times10^{-2}$, $5.5\times10^{-3}$, $6.5\times10^{-4}$,
$1.2\times10^{-4}$ at $\Omega t = 2.1, 3.1, 4.2, 5.4$ — which is exactly the
not-yet-separated fraction near $u = 0$. The closed form is the limit.

**Where the quantum content lives.** For a minimum-uncertainty packet
$`\sigma_r\sigma_p = \hbar/2`$, so

```math
\sigma_u^2 \;\ge\; \frac{\hbar}{\mu\Omega},
```

attained at $`\sigma_r^2 = \hbar/(2\mu\Omega)`$. The transmission therefore
switches from 0 to 1 over a window of width
$`\sqrt{\hbar/\mu\Omega}`$ in $\bar u$ — a purely quantum smearing that
vanishes as $\hbar\to0$, where $T$ becomes a step function (Part C2:
$T = 0.663, 0.724, 0.800, 0.908, 0.999$ at $\bar u = 0.5$ for
$\hbar = 1, 0.5, 0.25, 0.1, 0.02$).

**This is the whole point of the test case.** The evolution is exactly
classical, with no Moyal correction at any order, yet the answer is a quantum
tunnelling probability. Everything non-classical is carried in the *width* of
$W$, transported rigidly. A microdynamic scheme that adds any diffusion to the
momentum distribution will widen $\sigma_u$ and inflate $T$ toward $1/2$; a
scheme that damps the initial spread will sharpen it. Both failures are
visible against a closed form.

For the energy-resolved comparison the reference is Kemble,
$`T(E) = \bigl(1 + e^{-2\pi E/\hbar\Omega}\bigr)^{-1}`$ (Part C3).

---

## 3. Why this is a two-body test and not a disguised one-body test

The separation is exact, so one could object that this is just a
single inverted oscillator at reduced mass. Three things make it a genuine
two-body test of the microdynamics.

**(i) The mode is anti-diagonal, so Theorem A3 is exercised directly.** The
pair potential's joint wavevectors are $`\vec q = (q, -q)`$, orthogonal to
$(1,1)$. By the leak law of
[`../analysis/fourd_microdynamics.md`](../analysis/fourd_microdynamics.md) §2.3,
*every one of the four actions* — Focus, Defocus, and both hops — leaves
$P = p_1 + p_2$ untouched, event by event. Part D confirms $\mathrm{Var}(P)$
is constant to $2\times10^{-8}$ over $t \le 8$ (the drift is floating-point in
the matrix exponential, growing with $`e^{2\Omega t}`$). A joint implementation
that leaks any total momentum has a bug in the anti-diagonal stencil, and this
catches it with an exact zero.

**(ii) The separation supplies a cheap exact reference for the expensive
run.** Because the joint 4-D problem separates by theorem, a joint lattice run
*must* reproduce a $1+1$ D run in the relative coordinate at mass $\mu$. That
makes the correlated-jump machinery testable against a reference that costs
$`O(M_x M_p)`$ instead of $`O((M_x M_p)^2)`$. Very few interacting two-body
problems offer that.

**(iii) It generates entanglement from a product state, at a rate set by the
Lyapunov exponent.** This is the observable a mean-field solver cannot
produce. Starting from an unentangled product of two Gaussians (Part D):

| $t$ | 2 | 4 | 6 | 8 | 10 | 12 |
|---|---|---|---|---|---|---|
| reduced purity | 0.0686 | 0.00169 | $6.3\times10^{-5}$ | $2.7\times10^{-6}$ | $1.2\times10^{-7}$ | $6.1\times10^{-9}$ |
| $S$ (nats) | 2.985 | 6.690 | 9.986 | 13.133 | 16.203 | 19.226 |
| $dS/dt$ | 2.084 | 1.712 | 1.602 | 1.551 | 1.522 | 1.503 |
| $\Omega + 1/t$ | 1.914 | 1.664 | 1.581 | 1.539 | 1.514 | 1.498 |

The asymptotic law is

```math
S(t) \;\simeq\; \Omega t + \ln t + \mathrm{const},
\qquad
\frac{dS}{dt} \;\to\; \Omega + \frac{1}{t},
```

and the residual after subtracting the $1/t$ term falls from
$1.7\times10^{-1}$ to $5.1\times10^{-3}$ over $t = 2 \to 12$. The $\ln t$
correction has a clean origin: the relative mode alone contributes a
*constant* reduced determinant, since squeezing preserves it; the growth comes
from the mismatch between the exponentially stretched relative mode and the
diffusively spreading free centre of mass, giving
$`\nu \sim t\thinspace e^{\Omega t}`$.

Note that Proposition B3 of the 4-D note reads this directly in
crystal-lattice terms: the reduced purity *is* the excess-to-background ratio,
so a two-particle scattering event drives the one-particle signal down by nine
orders of magnitude in twelve time units. That is the sign/shot-noise problem
appearing in a concrete, physically motivated setting rather than as an
abstract scaling argument.

Beyond $t \approx 14$ the covariance route fails: $`\det\Sigma_1`$ is a
difference of numbers of order $`e^{2\Omega t}`$ and float64 loses it. That is
a limitation of the *reference* calculation, not of the physics, and it bounds
how far the test can be pushed with a Gaussian benchmark.

---

## 4. The microdynamic representation, and what it costs

### 4.1 The noise is exactly as bad as the trap's

The periodic image of $`-\tfrac{\mu}{2}\Omega^2 r^2`$ on
$[-L/2, L/2)$ has

```math
V_q \;=\; -\thinspace\frac{\mu\Omega^2 L^2}{2\pi^2}\thinspace\frac{(-1)^q}{q^2},
```

the negative of the trap's coefficients. Since Theorem D of the 4-D note
depends only on $`|V_q|`$,

```math
\frac{|V_q|}{\hbar}\bigl(2q\thinspace\Delta p\bigr)^2 \;=\; 2\thinspace\mu\thinspace\Omega^2\hbar \;=\; 2\thinspace m\thinspace\omega^2\hbar,
```

independent of $q$, verified to $2.2\times10^{-16}$ (Part E1). **The inverted
barrier inherits the mode-sum noise pathology in full.** Correction 1 of §0.1.
The rate budget is likewise $`\sum_q |V_q|/\hbar = \mu\Omega^2 L^2/12\hbar`$,
which was $132.7$ at $L = 40$ in the reference run — large enough that the
exact FFT jump substep (§5.6 of the 4-D note, open item 5) is not optional if
the deterministic mesh path is used.

### 4.2 The escape problem

Both the position and the momentum extent of the packet grow as
$`e^{\Omega t}`$, so a periodic box of length $L$ with $N$ momentum cells
holds the state for a time logarithmic in the grid size. Since the accuracy of
a *literally propagated* transmission improves only as $`e^{-\Omega t}`$,
achieving accuracy $\epsilon$ needs an $x$-window and a $p$-window each
$`\propto 1/\epsilon`$, hence (Part E3)

| $\Omega t^{*}$ | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|
| box length $L$ | 20 | 55 | 148 | 403 | 1094 |
| cells $N$ needed | 90 | 668 | 4937 | 36483 | 269573 |
| accuracy $\epsilon$ | $3.7\times10^{-1}$ | $1.4\times10^{-1}$ | $5.0\times10^{-2}$ | $1.8\times10^{-2}$ | $6.7\times10^{-3}$ |

so $`N \sim \epsilon^{-2}`$.

### 4.3 Why it does not bind

Because $T$ is a *conserved* functional (§2.1), it need not be measured at
late times at all. Tracked on the crystal lattice with the exact jump substep
(Part E2, $512\times512$, $L = 40$, $q_{\max} = 128$):

| $t$ | 0.00 | 0.25 | 0.50 | 0.75 | 1.00 | 1.25 | 1.50 |
|---|---|---|---|---|---|---|---|
| $`\int_{u>0}W`$ | 0.500182 | 0.499511 | 0.500567 | 0.500135 | 0.499437 | 0.499512 | 0.499971 |
| error vs $T$ | $3.8\times10^{-4}$ | $2.9\times10^{-4}$ | $7.6\times10^{-4}$ | $3.3\times10^{-4}$ | $3.7\times10^{-4}$ | $2.9\times10^{-4}$ | $1.7\times10^{-4}$ |

against the closed-form $T = 0.49980392$. The total drift over the run is
$1.6\times10^{-3}$ with no systematic trend; the offset already present at
$t = 0$ is grid discretisation of the initial Gaussian, not dynamics. So the
diagnostic works over the whole window, with no requirement that anything
escape.

![Diagnostics for the inverted pair barrier](https://raw.githubusercontent.com/billpage/wpmw/output/figures/inverted_pair_barrier_diagnostics.png)

---

## 5. Recommended protocol

A concrete specification, in increasing order of cost.

**Stage 1 — relative coordinate, deterministic mesh.** Single particle at
$\mu = m/2$ on a ring, potential $`-\tfrac{\mu}{2}\Omega^2 r^2`$ truncated at
$q_{\max}$, exact FFT jump substep. Check: $`\int_{u>0}W`$ against
$`\Phi(\bar u/\sigma_u)`$ at every step. Establishes the discretisation floor
before any stochastics.

**Stage 2 — relative coordinate, four-action Monte Carlo.** Same setup with
`step_jump_four_rule_mc`. The transmission error now has a shot-noise
component; because the injected momentum variance is
$`2\mu\Omega^2\hbar`$ per mode, expect the error in $T$ to grow as
$`\sqrt{q_{\max}}`$. Verifying that prediction is itself a test of Theorem D
in a setting where the observable is physical rather than a norm.

**Stage 3 — joint 4-D lattice, deterministic.** Two particles, anti-diagonal
modes only. Two independent checks with exact answers: $\mathrm{Var}(P)$
constant (Theorem A3), and the relative marginal reproducing Stage 1 exactly
(separation theorem). This is the first test of the correlated-jump stencil
with a reference that costs less than the run.

**Stage 4 — joint sea-dressed / world-ensemble.** Add the external trap back
if a bound reference is wanted, or keep it pure. The target observable is the
reduced purity of §3(iii), which a mean-field solver returns as identically 1.
This is where Propositions B2 and B3 of the 4-D note become operational, and
where the sign-versus-shot-noise trade can be measured on a problem with a
closed-form answer.

**A variant worth running.** Replace the initial Gaussian with a state whose
Wigner function has negative regions straddling the separatrix $u = 0$. Then
$`\int_{u>0}W`$ is a difference of positive and negative contributions, and
the signed-particle Monte Carlo must get a cancellation right to reproduce a
number that is still exactly conserved. That is a much sharper stress test
than the Gaussian, and it costs nothing extra to set up.

---

## 6. Numerical verification

All from `src/demo_inverted_pair_barrier.py`, container run August 2026,
$m = \hbar = \omega = 1$, packet $`r_c = -4`$, $`p_c = 2.828`$,
$`\sigma_r = 0.7`$ (so $`\bar u = -6.0\times10^{-4}`$, $\sigma_u = 1.229$,
$T = 0.49980392$).

- **Part A.** Flow spectrum $`\lbrace-\Omega,0,0,\Omega\rbrace`$ to
  $4.4\times10^{-16}$; $`\Omega = \sqrt{2}\thinspace\omega`$.
- **Part B.** $`u(t)/u(0) = e^{\Omega t}`$ and $`v(t)/v(0) = e^{-\Omega t}`$
  to $`8.9\times10^{-16}`$ at $t = 0.3$, degrading to $`3.6\times10^{-12}`$
  at $t = 3.5$ as the exponential itself grows.
- **Part C.** Schrödinger convergence to the closed form as $`e^{-\Omega t}`$;
  $\hbar$-scan of the transmission window; Kemble table.
- **Part D.** $\mathrm{Var}(P)$ constant to $2\times10^{-8}$; entanglement
  table of §3, and $dS/dt$ against $\Omega + 1/t$.
- **Part E.** Noise constant $`2\mu\Omega^2\hbar`$ to $2.2\times10^{-16}$;
  lattice tracking of $`\int_{u>0}W`$; grid-scaling table.

---

## 7. Open items

1. **Run Stage 2.** The predicted $`\sqrt{q_{\max}}`$ growth of the
   transmission error is a falsifiable consequence of Theorem D on a physical
   observable; it has not yet been measured.
2. **The negative-Wigner variant** of §5 is specified but not run. It is the
   natural vehicle for the compound-channel cancellation question left open in
   `../analysis/sea_dressed_microdynamics.md`.
3. **A non-Gaussian reference beyond $t \approx 14$.** The covariance route
   loses $`\det\Sigma_1`$ to float64; extended precision or a direct
   two-particle Schrödinger solve would be needed to push the entanglement
   benchmark further.
4. **Absorbing or matched boundaries.** *Answered, and the question was
   wrongly put; see
   [`../analysis/open_position_space.md`](../analysis/open_position_space.md).*
   There is no layer, because there is no boundary: both substeps of the
   microdynamics are pointwise in $`x`$, so worlds need no boundary
   condition and the escape problem of §4.2 is an artefact of the mesh.
   What the box was hiding is a different pathology — the Wigner kernel's
   modulus is independent of position for *every* potential (Theorem O1
   there), so a world far downstream is struck at full rate and its
   correct free behaviour is a cancellation whose fringe frequency grows
   as $`2x/\hbar`$. No absorber helps with that, because nothing is
   misbehaving; it is the estimator that degrades. The usable device is a
   **coherence horizon**: a bound $`L_c`$ on the ket–bra separation, which
   confines all vertex activity to within $`L_c/2`$ of the support of
   $`V`$, conserves signed number exactly, and costs essentially nothing
   once it clears that support. A complex absorbing potential is available
   and needs no new stencil geometry — the commutator takes the difference
   of the two offset rows, the anticommutator their sum — but is strictly
   worse on three counts (§7 there). The periodic-potential route of §4
   there does **not** apply to this test case: $`-\tfrac{\mu}{2}\Omega^2 r^2`$
   has no period.

---

## 8. Sources

- [`../analysis/fourd_microdynamics.md`](../analysis/fourd_microdynamics.md)
  — Theorem A3 (leak law), Theorem D (noise constant), Proposition B3
  (purity as excess-to-background), and §5.7, corrected here.
- [`../analysis/four_rule_microdynamics_equivalence.md`](../analysis/four_rule_microdynamics_equivalence.md)
  — the four-action channels whose noise is being counted.
- [`../algorithm/multi_body_extension.md`](../algorithm/multi_body_extension.md)
  — the joint-lattice specification that Stages 3–4 exercise.
- Generating script: `src/demo_inverted_pair_barrier.py`. Figures published to
  the `output` branch as `figures/inverted_pair_barrier_*.png`.
