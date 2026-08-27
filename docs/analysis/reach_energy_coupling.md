# What the reach actually controls: period, energy, and the $`\hbar^2`$ term

**Status.** Analysis note, step 11 of the ladder. Companion demo:
`src/demo_reach_energy_coupling.py`. Prompted by the question of whether the
coherence horizon of
[`open_position_space.md`](open_position_space.md) §3 is necessary at all, and
by open items 1 and 4 of that note.

---

## 0. What this note settles, retracts and corrects

The reach $`y_{\max}`$ is not removable — it is what makes the momentum
lattice exist — but it is misnamed, and almost everything attributed to it in
the existing notes belongs somewhere else. Sorting that out closes two open
items and corrects three statements.

**Settles.**

- **Open item 1 of `open_position_space.md` §9** ("whether a taper exists that
  preserves both the killed tails and the exact momentum lattice is
  unresolved"). Resolved: don't taper. A sharp window *commensurate* with the
  potential's period is exact on the lattice, and folding the tail back in is
  exact for a decaying potential. §1 and §2.
- **Open item 4** ("Energy is not checked, and there is no reason to expect it
  to be exact"). Resolved: energy is exact, for the same reason worlds are.
  Both are moment conditions that follow from oddness, which every even window
  preserves. §3.

**Corrects.**

1. `open_position_space.md` §3, Definition (R), calls the reach "the greatest
   distance from its own position at which a world consults the potential".
   That describes an *aperture*. What the model needs is a *period*: a
   momentum lattice exists if and only if the separation coordinate is
   compact, which finite support alone does not deliver. §1.
2. `open_position_space.md` §3.2 concludes that a sharp horizon is the wrong
   tool for a crystal potential. The divergence measured there is the
   continuous-$`\xi`$ integral $`\int|V_W^{L_c}|\thinspace d\xi`$, which the
   model never evaluates. On the momentum lattice a commensurate sharp window
   is exact to six figures at every rung count. §1.
3. `open_position_space.md` §4.1 calls (P) "a different mechanism in kind"
   from (R) and (H). All three are the same mechanism — a period in the
   separation coordinate — differing only in whether the period comes from the
   state, from the potential, or from a postulate. §1.

**Contradicts, gently, its own first draft.** An earlier version of this
analysis held that the reach couples to Cyganski's energy balance once
$`V`$ has more than one harmonic. It does not: §3 and §4 show the energy
ledger is untouched by the reach at any number of modes. The coupling is real
but it lives one order higher, in the $`\hbar^2`$ Moyal term — and §6 shows
even that is an artefact of applying the horizon profile to the wrong operator.

---

## 1. The reach is a period, not an aperture

**Theorem E1 (one mechanism, three sources).** Write

```math
D_x(y) \;\equiv\; V(x+y) - V(x-y)
```

for the difference field in the half separation $`y`$. The momentum transfers
form a lattice if and only if $`D_x`$ is periodic in $`y`$, and then

```math
\Delta p \;=\; \frac{\pi\hbar}{\text{period of } D_x \text{ in } y} .
```

The three mechanisms of `open_position_space.md` §1 are three sources for that
one period:

| source | period in $`y`$ | $`\Delta p`$ |
|---|---|---|
| (R) ring: the state supplies it | $`2L`$ | $`\pi\hbar/L`$ |
| (P) periodic $`V`$: the potential supplies it | $`a`$ | $`\pi\hbar/a`$ |
| (H) horizon: postulated | $`L_c`$ | $`\pi\hbar/L_c`$ |

(P) is automatic and needs no comment beyond the algebra: if $`V`$ has period
$`a`$ then $`V(x \pm (y+a)) = V(x \pm y)`$, so $`D_x`$ inherits the period in
$`y`$ without anything being assumed about the state. That is Theorem O5
rewritten in the separation coordinate.

**Why finite support is not enough.** The transfer variable $`\xi`$ and the
half separation $`y`$ are Pontryagin duals, so a discrete $`\xi`$ spectrum
requires $`y`$ to be compact *as a group* — a circle, not an interval. A
compactly supported $`D_x`$ on the line has, by Paley–Wiener, an entire and
perfectly continuous transform. The existing construction already does the
right thing (the $`N_p`$ rungs are a DFT, which periodises); only the prose
describes a different object. **The Born–von Kármán boundary condition is the
right analogy, not a diffraction aperture** — and it predicts the correct
behaviour where the aperture picture does not, as the next result shows.

**Proposition E1.1 (commensuration, not sharpness).** For a potential of
period $`a`$, a sharp window of reach $`y_{\max}`$ is exact on the momentum
lattice when $`2y_{\max}`$ is an integer multiple of $`a`$, and fails
otherwise. Part A, $`V = \cos(2\pi x/4)`$ at $`x = 0.7`$, exact budget
$`2|\Gamma| = 1.782013`$:

| $`L_c`$ | $`L_c/a`$ | commensurate | sharp | tapered |
|---|---|---|---|---|
| 4 | 1.00 | yes | $`1.0000`$ | $`0.7500`$ |
| 8 | 2.00 | yes | $`1.0000`$ | $`1.0000`$ |
| 6 | 1.50 | no | $`6.2323`$ | $`0.9681`$ |
| 5 | 1.25 | no | $`4.6689`$ | $`0.8788`$ |

(ratios to the exact budget). The window's sidelobes are real off the lattice
but land exactly on lattice **zeros**, because the sidelobe spacing
$`\pi\hbar/L_c`$ *is* $`\Delta p`$. Note the second column of the tapered
case: **for a crystal the soft profile is strictly worse than a commensurate
sharp window**, because the taper attenuates the true delta.

This corrects `open_position_space.md` §3.2 as stated in §0. The tabulated
divergence there is genuine, but it is a property of the continuous-$`\xi`$
integral rather than of the discrete operator the algorithm builds.

---

## 2. Folding, and the trilemma

Given that the period must be imposed, there are two ways to impose it on
$`D_x`$: discard the tail (truncate) or wrap it back in (periodic summation).
`compensated_liouville_algorithm.md` §4.4 considers only the first, finds a
seam, and repairs the seam with a profile. The second has no seam to repair.

**Theorem E2 (folding identity).** For a decaying $`V`$, let
$`D_x^{\rm fold}(y) = \sum_n D_x(y + na)`$ with $`a = 2y_{\max}`$. Then

```math
K_q^{\rm fold} \;=\; \Delta p \thinspace \cdot \thinspace V_W^{(a)}(x, \xi_q) ,
\qquad V_a \;=\; \sum_n V(\cdot + na) ,
```

the exact Wigner kernel of the **periodised** potential, sampled on the
lattice. No amplitude is distorted.

*Proof sketch.* Folding $`D`$ in $`y`$ and periodising $`V`$ in $`x`$ are the
same operation, since $`\sum_n V(x-y-na) = \sum_n V(x-y+na)`$ over a symmetric
index set; and the Fourier coefficients of a periodic summation are samples of
the full-line transform. $`\square`$

Measured (Part B): $`\max|K^{\rm fold} - \Delta p\thinspace V_W| / \max|K| \le
3.4\times10^{-16}`$ at every $`(x, y_{\max})`$ tested. The first moment is
$`-V_a'(x)`$ to eight digits — the entire error is the image force and nothing
else. At $`y_{\max} = 4`$ the force defect is $`1.8\times10^{-7}`$ against the
sharp window's kernel error of $`4.2\times10^{-3}`$ at the same reach
(`open_position_space.md` §3.1).

**What folding costs: Proposition O4.** The images are everywhere, so no
world is ever free. Total jump rate at $`y_{\max} = 2`$:

| $`x`$ | 1 | 3 | 5 | 7 | 9 | 11 |
|---|---|---|---|---|---|---|
| sharp | $`1.50`$ | $`5.3\times10^{-1}`$ | $`5.0\times10^{-8}`$ | $`5.7\times10^{-22}`$ | $`7.5\times10^{-43}`$ | $`1.1\times10^{-70}`$ |
| tapered | $`7.3\times10^{-1}`$ | $`2.6\times10^{-3}`$ | $`4.8\times10^{-11}`$ | $`2.3\times10^{-25}`$ | $`1.7\times10^{-46}`$ | $`1.7\times10^{-74}`$ |
| folded | $`0.999`$ | $`0.999`$ | $`0.999`$ | $`0.999`$ | $`0.999`$ | $`0.999`$ |

against $`|V(x)|`$ falling to $`8\times10^{-106}`$. So:

> **The trilemma.** Locality (Proposition O4), exactness, and a momentum
> lattice: pick two — unless $`V`$ is genuinely periodic, in which case all
> three hold at once.

| realisation | period from | local | seam | amplitudes |
|---|---|---|---|---|
| sharp window | truncation | yes | yes | distorted; no third moment |
| tapered window | truncation | yes | no | attenuated |
| folded window | quotient | **no** | no | exact for $`V_a`$ |
| genuinely periodic $`V`$ | physics | no, honestly | no | exact |

Folding also has a floor: it needs $`a`$ comfortably larger than the support
of $`V`$. At $`a = 2`$ with $`\sigma = 0.5`$ the images overlap so strongly
that $`x = 1`$ becomes a symmetry point of the image array and the force
vanishes identically. This is Ewald's condition, and it is why folding is not
offered as the default.

---

## 3. Conservation does not depend on the reach

**Theorem E3 (number and energy are reach-independent).** For any window that
is **even** in $`y`$ — sharp, tapered or folded — the lattice kernel is odd in
$`\xi`$, hence

```math
M_0 \;=\; \sum_q K_q \;=\; 0
\qquad\text{and}\qquad
M_2 \;=\; \sum_q \xi_q^2 K_q \;=\; 0 ,
```

which are exactly the conditions for conservation of signed world number and
of energy.

*Proof.* $`D_x`$ is odd in $`y`$ and an even window preserves oddness, so the
Fourier coefficients satisfy $`K_{-q} = -K_q`$. Both sums pair off. For the
energy statement, the jump step changes the kinetic energy by

```math
\frac{d\langle T\rangle}{dt}\bigg|_{\rm jump}
= \frac{\langle p\rangle}{m}\int \xi\thinspace V_W\thinspace d\xi
\;+\; \frac{n}{2m}\int \xi^2\thinspace V_W\thinspace d\xi ,
```

so $`M_2 = 0`$ together with $`M_1 = -V'`$ gives exactly
$`-\langle (p/m)V'\rangle`$, the classical power removed by streaming.
$`\square`$

Measured (Part C), over four potentials $`\times`$ two reaches $`\times`$ two
windows: worst $`|M_0| = 1.8\times10^{-15}`$, worst
$`|M_2| = 7.9\times10^{-10}`$ (quadrature floor). Proposition O3 is the
$`M_0`$ half of this; the $`M_2`$ half closes open item 4.

---

## 4. The four-action ledger with more than one mode

**Theorem E4 (channel decomposition of the power).** For the symmetric member
of the four-rule family, at fixed $`x`$ and per mode $`q`$, the focus channel
does **no net work** and the hop channel delivers the **whole** classical
power.

*Proof.* $`\sum_n f_n = \tfrac{\Gamma}{2}\sum_n (W_{n+q} - W_{n-q})`$
telescopes to zero on a closed momentum lattice, and each focus event changes
$`T`$ by $`-\xi_q^2/m`$ independently of $`n`$, so the focus contribution
vanishes. Each hop across $`n`$ changes $`T`$ by $`2p_n\xi_q/m`$, and
$`\sum_n h_n p_n = -\Gamma\langle p\rangle`$, giving
$`d\langle T\rangle/dt = -2\xi_q\Gamma\langle p\rangle/m
= -\langle (p/m)V_q'\rangle`$ since $`M_1 = -2\xi_q\Gamma`$. Both steps are
linear in the mode index, so they sum over any comb. $`\square`$

Measured (Part D), $`V = \cos K_1x + 0.6\cos(2K_1x + 0.4)`$ on a
$`64\times96`$ lattice:

| channel | $`d\langle T\rangle/dt`$ |
|---|---|
| focus/defocus | $`+1.2\times10^{-18}`$ |
| right/left hop | $`+0.01852439`$ |
| classical power | $`+0.01852439`$ |
| residual | $`-2.1\times10^{-17}`$ |

This extends Cyganski's single-harmonic energy result
([`../supplement/four_action_foundations.md`](../supplement/four_action_foundations.md)
§1) to any commensurate multi-mode potential, unchanged.

**Caveat.** The momentum lattice used here is cyclic, which makes the
telescoping exact by construction; pushing the packet to the grid edge does
not break it. An **open**, capped momentum lattice is untested, and is where
the focus channel would be expected to acquire spurious energy. Recorded as
open item 1.

---

## 5. Where the reach does bite: the $`\hbar^2`$ term

**Theorem E5 (profile contamination of the third moment).** Let $`w`$ be an
even profile with $`w(0) = 1`$ and $`w'(0) = 0`$, applied to the **full**
kernel. Then

```math
\sum_q \xi_q^3\thinspace K_q
\;=\; \frac{\hbar^2}{4}V'''(x)
\;+\; \frac{3\pi^2\hbar^2}{8\thinspace y_{\max}^2}\thinspace V'(x)
```

for the raised cosine $`w(t) = \cos^2(\pi t/2)`$.

*Proof.* $`\sum_q \xi_q^n K_q = (\hbar/2i)^n g^{(n)}(0)/(i\hbar)`$ with $`g`$
the windowed separation field, giving $`(\hbar^2/8)\thinspace g'''(0)`$ at
$`n = 3`$. Expanding $`g = wD_x`$,

```math
(wD_x)'''(0) = w(0)D_x''' + 3w'(0)D_x'' + 3w''(0)D_x' + w'''(0)D_x
= 2V''' + 3w''(0)\thinspace 2V' ,
```

using $`D_x(0) = D_x''(0) = 0`$ by oddness, $`D_x'(0) = 2V'`$,
$`D_x'''(0) = 2V'''`$, and $`w''(0) = -\pi^2/2y_{\max}^2`$ for the raised
cosine. $`\square`$

So the reach leaves the classical ledger untouched and contaminates the
leading quantum coefficient with the classical **force**, at order
$`1/y_{\max}^2`$. Verified (Part E) to $`2\text{–}6\times10^{-4}`$ relative
across $`y_{\max} = 0.5 \ldots 8`$ on three potentials, including the case
$`y_{\max} = 0.5`$ for the Gaussian barrier where the measured third moment is
$`+7.47`$ against an exact $`-0.541`$ — the wrong sign and fourteen times the
magnitude.

![reach and the third moment](https://raw.githubusercontent.com/billpage/wpmw/output/figures/reach_energy_coupling_moments.png)

Note also that the contamination has nothing to do with having several modes.
It appears whenever the spectrum is not a comb — one mode or many.

---

## 6. Compensation removes it exactly

**Theorem E6 (the reach is inert at order $`\hbar^2`$ under compensation).**
Apply the profile to the compensated residual
$`D_x^{\rm res}(y) = D_x(y) - 2y\thinspace V'(x)`$ of
[`compensated_liouville_splitting.md`](compensated_liouville_splitting.md)
§2 rather than to the full kernel. Then

```math
\sum_q \xi_q^3\thinspace K_q^{\rm res} \;=\; \frac{\hbar^2}{4}V'''(x)
```

with no $`y_{\max}`$ dependence at all.

*Proof.* $`D_x^{\rm res}(0) = D_x^{\rm res}{}'(0) = D_x^{\rm res}{}''(0) = 0`$
— the compensation subtracts exactly the term linear in $`y`$, and the
quadratic term is absent by oddness — so in the expansion of Theorem E5 the
contamination term $`3w''(0)D_x^{\rm res}{}'(0)`$ vanishes identically,
leaving $`(wD_x^{\rm res})'''(0) = D_x^{\rm res}{}'''(0) = 2V'''`$.
$`\square`$

Measured (Part F), same three potentials, $`N_p = 16384`$:

| $`y_{\max}`$ | 0.5 | 1 | 2 | 4 | 8 |
|---|---|---|---|---|---|
| Gaussian, full | $`+7.47`$ | $`+1.46`$ | $`-0.040`$ | $`-0.416`$ | $`-0.510`$ |
| Gaussian, residual | $`-0.541486`$ | $`-0.541338`$ | $`-0.541283`$ | $`-0.541327`$ | $`-0.541337`$ |
| quartic, residual | $`0.509513`$ | $`0.509958`$ | $`0.510153`$ | $`0.510163`$ | $`0.510152`$ |

against exact values $`-0.541341`$ and $`+0.510000`$. The residual rows show
no monotone trend in $`y_{\max}`$; the figure above plots the errors on
log–log axes, where the full-kernel line has slope $`-2`$ and the residual
line is flat.

**Consequence for `compensated_liouville_algorithm.md` §4.4.** That section's
normative rule — apply the profile to the residual, and only to the residual —
is currently presented as bookkeeping. It is load-bearing: it is what makes
the reach ontologically inert at order $`\hbar^2`$. The requirements table
should gain a row:

| requirement | what it buys |
|---|---|
| applied to the compensated residual | no $`1/y_{\max}^2`$ classical contamination of the $`\hbar^2`$ term (E6) |

**A numerical caution on how tightly to quote this.** The raised cosine gives
a $`q^{-3}`$ coefficient tail, so the terms $`\xi^3 K_q`$ are $`O(1)`$ and the
third-moment sum converges by sign cancellation rather than by decay — it is
conditionally, not absolutely, convergent. At $`y_{\max} = 2`$ the partial
sums run $`-0.53762,\; -0.54041,\; -0.54111,\; -0.54128,\; -0.54125`$ for
$`N_p = 256 \ldots 65536`$; the last step moves away from the target. The
agreement is real at the $`10^{-4}`$ level but is not monotone, so §4.4's
six-figure quotation reads more into one reach than the sum supports. A
profile with $`w''(\pm1) = 0`$ would give a $`q^{-4}`$ tail and absolute
convergence; recorded as open item 2.

---

## 7. A polynomial potential has no jump measure at all

**Theorem E7 (termination and non-existence).** If $`V`$ is a polynomial of
degree $`d`$, the Moyal series terminates. For the quartic double well
$`V = ax^4 - bx^2`$, $`V^{(5)} \equiv 0`$ and the collision term is exactly

```math
V'(x)\thinspace\partial_p W \;-\; \frac{\hbar^2}{24}V'''(x)\thinspace\partial_p^3 W ,
\qquad\text{i.e.}\qquad
V_W \;=\; V'\thinspace\delta'(\xi) + \frac{\hbar^2}{24}V'''\thinspace\delta'''(\xi) .
```

A finite-order differential operator has unbounded total variation, so there
is **no finite jump rate on the open line**: the four-action microdynamics
does not exist for a polynomial potential until a reach is imposed.

For an unbounded potential, therefore, the reach is **constitutive, not
approximate**. There is nothing to refine away, and §8 shows there is nothing
to be gained by trying.

Two consequences worth recording.

- Cyganski's diagnosis in *The Extended Fokker–Planck equations and the QLE*
  attributes the impulsive jump density of the cubic case to the potential
  being *unbounded*, noting that bounded potentials give smooth, bounded
  Wigner potentials. That is right about smoothness but the operative property
  is different: what matters is whether $`V_W`$ is supported on a **comb** in
  $`\xi`$. A Gaussian barrier has a perfectly smooth bounded
  $`V_W \propto e^{-\xi^2/2}\sin 2\xi`$ and still has no momentum lattice,
  because its support is an interval. Boundedness does not rescue the lattice.
- Termination makes the quartic double well the cleanest available benchmark:
  the exact reference is a three-term PDE with no truncation ambiguity, so a
  live simulation can be scored against ground truth without the usual
  question of whether the reference has converged. Open item 3.

---

## 8. What the reach buys, and what it costs

**Theorem E8 (the residual budget always diverges).** Under the compensated
split, the residual event budget $`R_{\rm res} = \sum_q|K_q^{\rm res}|`$ grows
without bound as $`y_{\max} \to \infty`$ for every $`V`$ with
$`V'(x) \ne 0`$.

*Proof.* The compensation ramp $`-2yV'(x)`$ is unbounded in $`y`$, so
$`M_{\rm res} \to -iV'(x)s`$ — the anti-drift already noted in
`compensated_liouville_splitting.md` §2. The residual therefore inherits a
contribution growing linearly in the reach. For potentials with unbounded
$`V'''`$ the Taylor bound of Theorem C2 dominates instead and the growth is
cubic. $`\square`$

Measured (Part H), growth factor per doubling of the reach:

| $`y_{\max}`$ | Gaussian barrier | commensurate 2-mode | quartic double well |
|---|---|---|---|
| 2 | $`0.209`$ | $`0.425`$ | $`0.536`$ |
| 4 | $`1.272`$ ($`\times6.1`$) | $`2.127`$ ($`\times5.0`$) | $`4.284`$ ($`\times8.0`$) |
| 8 | $`3.220`$ ($`\times2.5`$) | $`6.205`$ ($`\times2.9`$) | $`34.27`$ ($`\times8.0`$) |
| 16 | $`6.306`$ ($`\times2.0`$) | $`14.995`$ ($`\times2.4`$) | $`274.2`$ ($`\times8.0`$) |
| 32 | $`12.00`$ ($`\times1.9`$) | $`28.37`$ ($`\times1.9`$) | $`2193.2`$ ($`\times8.0`$) |

![residual budget against reach](https://raw.githubusercontent.com/billpage/wpmw/output/figures/reach_energy_coupling_budget.png)

Bounded $`V`$ settles at $`\times2`$ — linear, from the ramp. The quartic
holds $`\times8`$ exactly — cubic, with no saturation.

Putting E3, E6 and E8 together gives the accounting the reach deserves:

| quantity | depends on the reach? |
|---|---|
| signed world number, $`M_0`$ | no (E3) |
| classical force, $`M_1`$ | no |
| energy, $`M_2`$ | no (E3) |
| leading $`\hbar^2`$ coefficient, $`M_3`$ | no, **under compensation** (E6) |
| momentum resolution $`\Delta p = \pi\hbar/2y_{\max}`$ | yes |
| event budget $`R_{\rm res}`$ | yes, divergently (E8) |

> **The reach buys momentum resolution and costs event budget. It buys no
> accuracy.** That is a cleaner statement than the reach has had so far, and
> it makes the choice of $`L_c`$ a resolution-versus-cost decision rather than
> a convergence parameter.

**A cost of compensation, for a comb.** The ramp $`-2yV'`$ is not periodic in
$`y`$, so subtracting it turns a finite delta comb into a dense spectrum:

| $`y_{\max}`$ | full: non-zero $`K_q`$ | $`R`$ | residual: non-zero | $`R_{\rm res}`$ |
|---|---|---|---|---|
| 4 | 4 | $`2.47794`$ | 2350 | $`2.127`$ |
| 8 | 4 | $`2.47794`$ | 2962 | $`6.204`$ |
| 32 | 4 | $`2.47794`$ | 4698 | $`28.367`$ |

On event count alone the crossover sits near $`y_{\max} \approx 4.2`$, which is
essentially the ring's own reach $`L/2 = 4`$. Whether compensation still pays
past that depends on variance in the force channel rather than on raw event
count, which this note does not measure. Open item 4.

![window realisations and their tails](https://raw.githubusercontent.com/billpage/wpmw/output/figures/reach_energy_coupling_windows.png)

---

## 9. Consequences for existing documents

- `open_position_space.md` §3 Definition (R) — restated as a period, with the
  Born–von Kármán reading; §3.2 corrected to commensuration; §4.1's "different
  in kind" softened; open items 1 and 4 closed. Carried by the same patch as
  this note.
- `../algorithm/compensated_liouville_algorithm.md` §4.4 — the requirements
  table should gain the compensated-input row of §6, and the six-figure
  quotation of the third moment should be relaxed to four. Not carried here;
  left for a separate patch so the algorithm spec's verification tables can be
  regenerated alongside it.
- `compensated_liouville_splitting.md` §2 — Theorem E6 gives the split a
  second justification independent of variance reduction: it is what keeps the
  horizon out of the semiclassical limit.

---

## 10. Open items

1. **The open momentum lattice.** Theorem E4's focus-does-no-work result uses
   a closed (cyclic) lattice, where $`\sum_n f_n`$ telescopes identically. On
   an open lattice with capped end cells it should not, and the focus channel
   should acquire spurious energy at the rate the packet reaches the edge.
   Untested.
2. **A profile with vanishing second derivative at the edge.** $`w''(\pm1) = 0`$
   would give $`q^{-4}`$ coefficients and hence an absolutely convergent third
   moment, removing the conditional-convergence caution of §6. Whether it
   costs anything in the other four requirements of §4.4 is unexamined.
3. **A live quartic-double-well benchmark.** §7 notes that the terminating
   Moyal series gives an unambiguous reference. Running the compensated
   four-action microdynamics against it at several reaches would turn E5 and
   E6 from kernel statements into dynamical ones, and would settle whether the
   $`1/y_{\max}^2`$ contamination shows up as measurable drift in
   $`\langle p^3\rangle`$ or is masked by noise.
4. **Variance, not event count.** §8's crossover is measured in events. The
   quantity compensation actually targets is variance in the delivered force,
   which needs a matched-$`\nu`$ Monte Carlo comparison.
5. **Folding under the compensated split.** Theorem E2 is a statement about
   the full kernel. Whether the folded realisation preserves Theorem C3's
   moment conditions, and whether the periodic $`V_a`$ it manufactures
   introduces unwanted superselection sectors via Theorem O5, are both open.

---

## 11. Numerical verification summary

All from `src/demo_reach_energy_coupling.py`.

| Part | Claim | Result |
|---|---|---|
| A | commensurate sharp window is exact | ratio $`1.0000`$ at $`L_c = 4, 8`$ |
| A | incommensurate sharp window fails | ratio $`6.23`$, $`4.67`$ |
| A | taper is not exact for a comb | ratio $`0.750`$ at $`L_c = a`$ |
| B | folding identity, E2 | $`\le 3.4\times10^{-16}`$ |
| B | folded force is $`-V_a'(x)`$ | to 8 digits |
| B | folding breaks Proposition O4 | rate $`0.999`$ out to $`x = 11`$ |
| C | $`M_0 = 0`$, any even window | worst $`1.8\times10^{-15}`$ |
| C | $`M_2 = 0`$, any even window | worst $`7.9\times10^{-10}`$ |
| D | focus channel does no work | $`1.2\times10^{-18}`$ |
| D | hop channel carries the power | residual $`2.1\times10^{-17}`$ |
| E | third-moment contamination, E5 | $`2\text{–}6\times10^{-4}`$ relative |
| F | contamination absent under compensation | no trend in $`y_{\max}`$ |
| F | hard cutoff on the residual still fails | $`M_3: -50 \to -12656`$ |
| G | residual $`M_3`$ for the quartic | $`0.5102`$ vs exact $`0.5100`$ |
| H | budget growth per doubling | $`\times1.9`$ bounded, $`\times8.0`$ quartic |

---

## 12. Sources

- [`open_position_space.md`](open_position_space.md) §1–§4 — the three
  mechanisms, Theorem O1, Propositions O3 and O4, Theorem O5.
- [`compensated_liouville_splitting.md`](compensated_liouville_splitting.md)
  §2 — $`M_{\rm cl}`$, $`M_{\rm res}`$, Theorems C2 and C3.
- [`../algorithm/compensated_liouville_algorithm.md`](../algorithm/compensated_liouville_algorithm.md)
  §4.4 — the horizon profile and its five requirements.
- [`four_rule_microdynamics_equivalence.md`](four_rule_microdynamics_equivalence.md)
  §3–§4 — the focus and hop channels and the exactness family.
- [`../supplement/four_action_foundations.md`](../supplement/four_action_foundations.md)
  §1 — Cyganski's single-harmonic energy balance, extended by Theorem E4.
- D. Cyganski, *The Extended Fokker–Planck equations and the QLE*, and *A
  journey from Bohm trajectory theory, through Nelson's SDEs and Wigner
  Particles to the Closed Four Action Model* (3 August 2026) — the split
  Fourier algorithm, whose momentum grid already fixes
  $`y_{\max} = \pi\hbar/2\Delta p`$, and the four-action energy balance.
- Generating script: `src/demo_reach_energy_coupling.py`. Figures published to
  the `output` branch as `figures/reach_energy_coupling_*.png`.
