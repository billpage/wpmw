# Compensated Liouville splitting: the classical force as deterministic acceleration

**Status.** Analysis note. Companion demo:
`src/demo_compensated_liouville_splitting.py`.

The crystal-lattice algorithm of
[`../supplement/phase_space_crystal_lattice_supplement.md`](../supplement/phase_space_crystal_lattice_supplement.md)
§7 free-streams first and then applies mediated jumps for every Fourier mode
of the potential, even where the dynamics is classical. The question this
note answers is whether the classical force can instead sit in the first
substep, leaving the jumps to carry only what is irreducibly quantum — the
phase-space analogue of the clean classical/quantum split that Bohm and
Nelson get for free.

The answer is yes. The split is exact, the quantum part is a signed jump
interaction that carries no force of its own, and the whole thing turns on
one condition: **a world must have a bounded coherence reach.**

---

## 0. What this note inherits and corrects

**Inherits.** From
[`../supplement/phase_space_crystal_lattice_supplement.md`](../supplement/phase_space_crystal_lattice_supplement.md)
§5–§7: the momentum quantum $`\Delta p = \pi\hbar/L`$, the mode stencil, the
rate field. From [`open_position_space.md`](open_position_space.md):
Theorem O1 (the Wigner kernel has no position envelope), Theorem O5
(momentum quantisation comes from periodicity of $`V`$, not from the box),
and above all the **coherence horizon** $`L_c`$ on the ket–bra separation,
which turns out to be the load-bearing idea here. From
[`four_rule_microdynamics_equivalence.md`](four_rule_microdynamics_equivalence.md):
the four channels and the exactness family.

**Errata to the previous revision (commit `71d6dce`).** Its §0 asserted that
"there is no inter-mode cancellation". That is false, and the case in which it
is false is the one that matters. The correct statement has two cases, split
by whether the potential is presented on the open line or on a ring, and it is
recorded here because it is also the cleanest argument for the reach
formulation of §3.

- **A true quadratic on the open line.** Its Fourier weight is a distribution
  supported at $`k = 0`$ involving derivatives only through second order. The
  residual weight of mode $`k`$ carries a factor $`k^3`$, which annihilates
  anything supported at the origin through second order, so every mode's
  residual is *separately* zero. Nothing cancels because there is nothing to
  cancel — but equally, no infinite mode sum is involved.
- **The periodised parabola on a ring.** Every mode's residual is separately
  nonzero and $`O(1)`$, and they cancel: exactly to zero inside the bowtie of
  §6.1, and to the seam value outside it. The cancellation is only
  *conditionally* convergent (§6.2). This is the case the prompting question
  had in mind, and its description of it was right.

What survives from the previous version's claim is only the exactness of each
mode's stencil, which is independent of the cancellation question and is
recorded as Lemma C0 in §6.2.

**Corrects the framing of the prompting question in one place**, on where
quantum mechanics enters.

It was put as two channels, the initial negativity of $`W`$ and the
mediated jumps. Negativity is neither necessary nor the sharp criterion: a
coherent or squeezed state has $`W \ge 0`$ everywhere, evolves by exact
classical Liouville flow in a harmonic well, and is fully quantum. What is
quantum about the initial datum is the $`\hbar`$-sized floor on phase-space
area — admissibility, $`|W| \le 2/h`$, Hudson's theorem — not the sign. And
there is a third channel: the Weyl correspondence at readout, where
$`\langle\hat A\rangle = \int A_W W`$ needs the Weyl symbol, which differs
from the classical function by $`O(\hbar^2)`$ for nonlinear observables. This
note is about the second channel.

---

## 1. The symbol

Let $`s`$ be conjugate to momentum and write $`y = \hbar s/2`$. This $`y`$ is
the half ket–bra separation: the distance at which a world consults the
potential. Under $`W(x, p \pm a) \leftrightarrow e^{\pm ias}\hat W(x,s)`$ the
whole potential term of the QLE collapses to multiplication by one scalar
function,

```math
\boxed{\;
M(x, s) \;=\; \frac{i}{\hbar}\bigl[V(x + y) - V(x - y)\bigr],
\qquad y = \frac{\hbar s}{2}.
\;}
```

Two features matter, and both are local statements about a single world.

- **No Fourier decomposition of $`V`$ is needed.** The transform is in $`y`$
  at fixed $`x`$, not in $`x`$. A world consults $`V`$ at $`x \pm y`$ over its
  reach and nothing else. The mode coefficients $`V_q`$ that the
  crystal-lattice algorithm uses are a *global* functional of $`V`$; the
  kernel's reach is a bounded neighbourhood. Only the second is forced. What
  the modes buy is discreteness of the jump spectrum, and by Theorem O5 that
  comes from periodicity of $`V`$.
- **Only odd derivatives survive.** $`V(x+y) - V(x-y)`$ is odd in $`y`$, so
  the even part of $`V`$ about $`x`$ cancels identically. The Moyal series is
  the Taylor expansion of $`M`$ in $`s`$.

---

## 2. The split

Subtract the part linear in the reach:

```math
M \;=\; \underbrace{i\thinspace V'(x)\thinspace s}_{M_{\mathrm{cl}}}
\;+\;
\underbrace{\frac{i}{\hbar}\bigl[V(x+y) - V(x-y) - 2y\thinspace V'(x)\bigr]}_{M_{\mathrm{res}}} .
```

$`M_{\mathrm{cl}}`$ is the classical Liouville force term with the **full**
$`V'`$ — no quadratic fitting, no choice of a preferred piece of the
potential.

**Theorem C1.** $`M_{\mathrm{cl}}`$ and $`M_{\mathrm{res}}`$ are both
multiplication operators in the same variables, hence commute exactly, and

```math
e^{\tau M} \;=\; e^{\tau M_{\mathrm{cl}}}\thinspace e^{\tau M_{\mathrm{res}}}
\qquad\text{for every }\tau ,
```

with no Trotter error. *Verified:* commutator identically zero;
factorisation error $`1.1\times10^{-16}`$, $`2.6\times10^{-16}`$,
$`9.5\times10^{-16}`$ at $`\tau = 0.01, 0.1, 1.0`$, with no growth in
$`\tau`$.

**Theorem C2.** Because the quadratic term of $`V`$ cancels in the odd
combination, $`M_{\mathrm{res}}`$ is exactly the odd part of the **cubic
Taylor remainder** of $`V`$ about $`x`$. Equivalently,

```math
M_{\mathrm{res}}(x, s) \;=\; \frac{i}{\hbar}\int_0^{y}
  \bigl[V'(x+\sigma) + V'(x-\sigma) - 2V'(x)\bigr]\thinspace d\sigma ,
\qquad
|M_{\mathrm{res}}| \;\le\; \frac{2}{\hbar}\thinspace\frac{y_{\max}^3}{6}\thinspace
  \sup_{\text{reach}}|V'''| .
```

So $`M_{\mathrm{res}} \equiv 0`$ if and only if $`V''' \equiv 0`$ on the
region the reach covers.

*Verified:* for the parabola on the open line,
$`\max|M_{\mathrm{res}}| = 1.4\times10^{-14}`$ against
$`\max|M| = 1.3\times10^{2}`$. For the cosine well the bound is tight to
within a factor under two:

| $`y_{\max}`$ | 0.25 | 0.50 | 1.00 | 2.00 |
|---|---|---|---|---|
| $`\max\lvert M_{\mathrm{res}}\rvert`$ | 0.00267 | 0.02125 | 0.16608 | 1.21084 |
| Taylor bound | 0.00315 | 0.02798 | 0.24224 | 1.93789 |
| ratio | 0.849 | 0.760 | 0.686 | 0.625 |

Per Fourier mode the split is the difference between a sine and its tangent at
the origin, $`\sin u`$ against $`u`$, with $`u = k\thinspace y_{\max}`$.

![The split of a mode symbol](https://raw.githubusercontent.com/billpage/wpmw/output/figures/compensated_symbol_and_band.png)

---

## 3. The reach theorem

Theorem C2 says the residual is bounded by a $`V'''`$ quantity *over the
reach*. That qualification is the whole story, and it is easy to miss.

If the reach is unbounded, $`M_{\mathrm{res}}`$ is not bounded at all. For a
bounded $`V`$, $`M`$ stays bounded as $`s \to \infty`$ while
$`M_{\mathrm{cl}} = iV's`$ grows without limit, so
$`M_{\mathrm{res}} \to -iV'(x)s`$: the residual acquires an *anti-drift*
exactly cancelling the classical force one just put into the first substep.
Realised naively that would apply the force twice, and the split would be
empty.

A world does not have unbounded reach. The coherence horizon $`L_c`$ of
[`open_position_space.md`](open_position_space.md) §3 caps the ket–bra
separation, which in the present variables is precisely the window
$`|y| \le L_c/2`$. Restricted to it, everything works.

**Theorem C3 (reach theorem).** Let $`y_{\max} = L_c/2`$. Restricted to
$`|y| \le y_{\max}`$, the residual kernel $`K_{\mathrm{res}}`$ — the momentum
transfer density obtained by transforming $`M_{\mathrm{res}}`$ — satisfies

```math
\int K_{\mathrm{res}}(x,\xi)\thinspace d\xi \;=\; 0
\qquad\text{and}\qquad
\int \xi\thinspace K_{\mathrm{res}}(x,\xi)\thinspace d\xi \;=\; 0 ,
```

and is bounded by the Taylor bound of C2. It is therefore a bounded signed
jump measure which conserves the number of worlds and delivers **no net
momentum**. The whole classical force $`-V'(x)`$ remains in the deterministic
step.

*Proof.* The zeroth moment vanishes because $`M`$ and $`M_{\mathrm{cl}}`$ both
vanish at $`s = 0`$. The first moment is $`-i\partial_s M_{\mathrm{res}}`$ at
$`s = 0`$, and $`M_{\mathrm{res}}`$ has been constructed with zero derivative
there. Boundedness is C2. $`\square`$

*Verified* (Part B), at $`x = 1`$, against the force each kernel would
otherwise double count:

| | $`-V'(x)`$ | first moment of $`K`$ | first moment of $`K_{\mathrm{res}}`$ |
|---|---|---|---|
| cosine well | $`-0.833041`$ | $`-0.833539`$ | $`3.2\times10^{-6}`$ |
| Gaussian barrier | $`+0.541341`$ | $`+0.541679`$ | $`1.2\times10^{-5}`$ |
| soft-core Coulomb | $`-0.715542`$ | $`-0.715983`$ | $`-1.1\times10^{-5}`$ |

(at $`y_{\max} = 0.25`$; the zeroth moments are at the $`10^{-17}`$ level
throughout). The full kernel carries the entire classical force, as it must —
this is the first-moment result of `open_position_space.md`. The residual
kernel carries five to six orders less.

So the algorithm becomes:

> **Step 1.** Deterministic acceleration by the full classical force
> $`-V'(x)`$, ideally fused with free streaming into one symplectic move,
> which also removes the existing free/potential Trotter error.
>
> **Step 2.** Signed hops drawn from $`K_{\mathrm{res}}`$: bounded, sourced by
> $`V'''`$, conserving world count and carrying no net momentum.

This is a genuine classical/quantum split. Bohm and Nelson achieve one
because there the quantum part is another *force*, so the decomposition is
automatic. Here the quantum part is a different kind of object — a zero-mean
hop interaction, focus and hop — and the price of that difference is a
condition on the reach rather than nothing at all.

![Newtonian arcs plus zero-mean hops at three reaches](https://raw.githubusercontent.com/billpage/wpmw/output/figures/compensated_world_paths.png)

---

## 4. The reach condition

How short must the reach be? The natural measure of event traffic is the
kernel's total variation, and compensation is worth doing exactly when
$`TV(K_{\mathrm{res}}) < TV(K)`$.

**Theorem C4.** Per mode, $`|M_{\mathrm{res}}|/|M_{\mathrm{cl}}| =
|\sin u - u|/|u| = u^2/6 + O(u^4)`$ with $`u = k\thinspace y_{\max}`$, and the
ratio reaches $`1`$ at $`u = \pi`$. The budget ratio follows, crossing unity
near $`u = \pi/2`$ — a reach of about a quarter wavelength of the potential's
dominant mode.

*Verified* (Part C), cosine well at $`x = 1`$ with $`k = \pi/4`$:

| $`y_{\max}`$ | $`u/\pi`$ | $`TV(K)`$ | $`TV(K_{\mathrm{res}})`$ | ratio | $`\Delta p = \pi\hbar/2y_{\max}`$ |
|---|---|---|---|---|---|
| 0.25 | 0.0625 | 2.226 | 0.013 | 0.006 | 6.283 |
| 0.50 | 0.125 | 4.374 | 0.104 | 0.024 | 3.142 |
| 1.00 | 0.25 | 8.139 | 0.817 | 0.100 | 1.571 |
| 2.00 | 0.50 | 11.925 | 5.987 | 0.502 | 0.785 |
| 4.00 | 1.00 | 2.121 | 33.703 | 15.888 | 0.393 |

Short reach buys two orders of magnitude; long reach costs one. The Gaussian
barrier and soft-core Coulomb behave the same way, crossing over between
$`y_{\max} = 1`$ and $`2`$.

The last column is not decoration. **Reach and momentum quantum are the same
parameter:** a reach $`y_{\max}`$ resolves momentum transfers only down to
$`\Delta p = \pi\hbar/(2y_{\max})`$. At maximal reach on a ring,
$`y_{\max} = L/2`$ recovers the crystal quantum $`\pi\hbar/L`$ and the kernel
collapses to the two-atom mode stencil — checked directly:
$`TV(K) = 2.121320`$ against the analytic $`2|\Gamma_q| = 2.121320`$. Coarse
momentum is exactly what makes the residual a finite jump measure rather than
a differential operator, so the two facts are one fact.

The physical reading of $`u \ll \pi`$ is the semiclassical condition in
disguise. Writing $`y_{\max} \sim \hbar/(2\sigma_p)`$ for a state of momentum
width $`\sigma_p`$, the condition is $`\sigma_p \gg \hbar k/2`$: **a world's
momentum width must exceed the kick it would receive.** In the limit
$`L_c \to 0`$ the residual vanishes as $`L_c^3 V'''`$ and motion becomes
exactly Newtonian — the classical limit, with the right scaling.

![Event budget and the vanishing first moment](https://raw.githubusercontent.com/billpage/wpmw/output/figures/compensated_reach.png)

---

## 5. The quiet region

Theorem C2 bounds the residual by the supremum of $`|V'''|`$ *over the reach*.
On the open line that turns immediately into a statement about position: there
are places a world can be where it takes no events at all.

**Theorem C7 (quiet region).** Fix a reach $`y_{\max}`$. If $`V''' \equiv 0`$
on $`[x - y_{\max},\; x + y_{\max}]`$ then
$`M_{\mathrm{res}}(x, \cdot) \equiv 0`$: a world at $`x`$ takes no events
whatever and moves on an exact Newtonian trajectory in the full potential.

*Proof.* By C2 the residual is the odd part of the cubic Taylor remainder
$`R_x(y)`$, and $`V''' \equiv 0`$ on the reach makes $`R_x \equiv 0`$ there.
$`\square`$

The condition is sufficient but **not** necessary. What is necessary and
sufficient is that $`R_x`$ be an *even* function of $`y`$ on the reach, and
parity can arrange that at isolated points without $`V'''`$ vanishing — at a
symmetry point of $`V`$, for instance, where the odd part of the remainder
cancels for reasons having nothing to do with the reach. The verification
below shows both kinds of silence: quiet by C7 outside the bump, and quiet by
parity at $`x = 0`$.

*Verified* (Part G) on a harmonic trap plus a $`C^\infty`$ bump supported on
$`[-1, 1]`$, so that $`V''' \equiv 0`$ exactly outside the bump. The predicted
quiet region is $`|x| > 1 + y_{\max}`$, and the edge is sharp to machine
precision at exactly that point for every reach:

| reach $`y_{\max}`$ | $`x = 0`$ | edge $`-\thinspace 0.3`$ | edge $`-\thinspace 0.05`$ | edge $`+\thinspace 0.05`$ | edge $`+\thinspace 3`$ |
|---|---|---|---|---|---|
| 0.5 | 0 (parity) | 1.41e-01 | 3.51e-05 | 6.7e-16 | 4.4e-15 |
| 1.0 | 0 (parity) | 1.41e-01 | 3.51e-05 | 1.3e-15 | 6.2e-15 |
| 2.0 | 0 (parity) | 1.41e-01 | 3.51e-05 | 3.6e-15 | 7.1e-15 |

### 5.1 The reach restores locality

A potential without compact support has no exactly quiet region, but the same
mechanism still localises the interaction. For a Gaussian barrier of width
$`\sigma`$ sitting on a harmonic trap, the residual at $`x`$ tracks the barrier
profile **translated outward by exactly the reach**:

| $`x`$ | 1.0 | 1.5 | 2.0 | 2.5 | 3.0 |
|---|---|---|---|---|---|
| $`\max\lvert M_{\mathrm{res}}\rvert`$ at $`y_{\max} = 1`$ | 4.75e-01 | 4.41e-01 | 4.38e-02 | 8.84e-04 | 3.73e-06 |
| shifted profile | 1.00e+00 | 4.58e-01 | 4.39e-02 | 8.84e-04 | 3.73e-06 |
| bare profile | 4.39e-02 | 8.84e-04 | 3.73e-06 | 3.29e-09 | 6.10e-13 |

Three significant figures of agreement with the shifted profile from $`x = 2`$
outward, against a bare profile seven orders of magnitude smaller by
$`x = 3`$.

This is the finite-reach refinement of **Theorem O1** of
[`open_position_space.md`](open_position_space.md), which says the untruncated
Wigner kernel's modulus has no position envelope — a world arbitrarily far
from a scatterer is still struck at an appreciable rate. With a bounded reach
that fails, and it fails in the sharpest possible way: the interaction is
supported within $`y_{\max}`$ of the non-quadratic part of $`V`$, exactly.
**The coherence horizon is what restores locality of the interaction in
position space.** O1 and C7 are the same statement at $`y_{\max} = \infty`$ and
at $`y_{\max} < \infty`$.

![The quiet region](https://raw.githubusercontent.com/billpage/wpmw/output/figures/compensated_quiet_region.png)

---

## 6. The ring as a diagnostic

The periodic box is not where this ontology lives — §5 is. But the periodised
parabola is an unusually informative test object, and three of its properties
sharpen the open-line case. One of them is a warning about ring numerics.

**Theorem C5.** Let $`V`$ be smooth on the circle of circumference $`L`$.
Then $`M_{\mathrm{res}} \equiv 0`$ if and only if $`V`$ is constant.

*Proof.* By C2, $`M_{\mathrm{res}} \equiv 0 \iff V''' \equiv 0 \iff V`$
locally quadratic; single-valuedness forces the linear and quadratic
coefficients to vanish. $`\square`$

So the quiet region of C7 is empty on a ring for every non-constant potential.
That makes the ring a stress test rather than a home — but a stress test with
closed-form answers, which is what the rest of this section trades on.

### 6.1 The bowtie: C7 where the potential is nowhere quadratic

The periodised parabola is nowhere globally quadratic — it has a kink at the
seam — yet it has an explicit quiet region, and the geometry is exact. It is
therefore the cleanest available demonstration that C2 and C7 are statements
about $`V`$ *on the reach* and not about $`V`$ globally. The seam plays the
role that the compact bump plays in §5.

**Proposition C6.1.** $`M_{\mathrm{res}}(x, s) = 0`$ exactly if and only if
$`|x| + |y| < L/2`$.

*Proof.* If both arms stay inside the fundamental domain the wrap is inactive
and $`V(x+y) - V(x-y) = 2y\thinspace V'(x)`$. Otherwise one arm crosses the
seam. $`\square`$

Measured: residual $`1.1\times10^{-14}`$ throughout the diamond, occupying
$`0.4961`$ of the grid against a predicted $`0.5039`$. At
$`\Delta p = \pi\hbar/L`$ the maximum arm is exactly $`L/2`$, so the diamond is
inscribed in the computational domain.

![Residual support is a bowtie](https://raw.githubusercontent.com/billpage/wpmw/output/figures/compensated_ring_residual.png)

**Corollary C6.2.** A coherence horizon $`L_c`$ fixes both obstructions at
once. It puts the seam out of reach — the ring parabola becomes exactly
Newtonian for every $`|x| < L/2 - L_c/2`$, verified at $`3.6\times10^{-15}`$
for $`L_c = 2, 4, 6`$ — and it moves worlds off the maximal-reach point where
C4 says nothing can be gained.

This is the strongest reading the coherence horizon has been given. It is not
a convenience for bounding jump traffic; it is the condition under which the
classical/quantum split exists at all.

### 6.2 Modes cancel, and badly

**Lemma C0 (the mode stencil is exact).** For $`V(x) = A\cos(kx + \phi)`$ the
exact symbol is

```math
M(x, s) \;=\; -\thinspace\frac{2iA}{\hbar}\thinspace\sin(kx + \phi)\thinspace\sin u ,
\qquad u = k\thinspace y ,
```

whose inverse transform is exactly two atoms at $`\xi = \pm\hbar k/2`$. So
replacing $`\partial_p W`$ by a centred finite difference at spacing
$`\hbar k/2`$ is an **identity**, to all orders in $`\hbar`$, not an
approximation. The classical term's symbol is the same expression with
$`\sin u`$ replaced by $`u`$ — so the derivative is the approximation to the
finite difference, and not the other way round.

`../supplement/phase_space_crystal_lattice_supplement.md` §6.1 states this
inverted, describing its own exact result as an approximation; see CLS5.

Now sum the modes of the periodised parabola at a point *inside* the bowtie,
where §6.1 says the residual is exactly zero:

| modes summed $`Q`$ | 1 | 20 | 100 | 500 | 4000 |
|---|---|---|---|---|---|
| signed sum | $`-1.4359`$ | $`-0.2114`$ | $`-0.0422`$ | $`-0.0084`$ | $`+0.0011`$ |
| largest single term | 1.436 | 4.116 | 4.116 | 4.116 | 4.116 |
| $`\sum\lvert\text{terms}\rvert`$ | 1.44 | 29.1 | 49.0 | 68.8 | 94.3 |

The signed sum converges to zero like $`1/Q`$. The individual terms do not get
small — the largest stays at 4.12 forever — and the sum of absolute values
diverges logarithmically, about 12.3 per e-fold of $`Q`$. This is genuine
cancellation between distinct modes, and it is only conditionally convergent.

That has a direct consequence for how one should compute. Bounding the
residual mode by mode gives $`\sum |V_q| k_q^3`$, which diverges — as it must,
since $`V''' = -m\omega^2 L\thinspace\delta'`$ at the seam has no absolutely
convergent Fourier series. Bounding it by reach, via C2, gives zero in one
line, because both arms sit inside the fundamental domain where $`V`$ is
quadratic. **The mode representation converts an exactly-zero local quantity
into a conditionally convergent sum of $`O(1)`$ terms whose absolute sum
diverges.** The remark in §1 that no Fourier decomposition of $`V`$ is needed
is therefore not a convenience: it is what makes C2 and C7 provable at all.

### 6.3 Why ring numerics mislead about the reach condition

§6.1 and §6.2 are the ring being useful. This one is the ring being
dangerous.

At maximal reach $`y_{\max} = L/2`$ the two arms $`x + y`$ and $`x - y`$ are
the *same point* — they meet around the back of the circle. So
$`V(x+y) - V(x-y) = 0`$ identically, the whole symbol vanishes, and
$`M_{\mathrm{res}} = -M_{\mathrm{cl}}`$ exactly, for every mode and for the
sum. Checked directly: $`\sin(q\pi)`$ evaluates to
$`1.2\times10^{-16}`$, $`-2.4\times10^{-16}`$, $`3.7\times10^{-16}`$ for
$`q = 1, 2, 3`$, and the full symbol at $`y = L/2`$ is $`0.000`$.

| mode $`q`$ | 1 | 2 | 3 | 4 |
|---|---|---|---|---|
| $`u`$ | $`\pi`$ | $`2\pi`$ | $`3\pi`$ | $`4\pi`$ |
| $`\lvert M_{\mathrm{res}}\rvert/\lvert M_{\mathrm{cl}}\rvert`$ | 1.0000 | 1.0000 | 1.0000 | 1.0000 |

So a ring pins every world at exactly the crossover point of C4, for every
mode, by geometry rather than by accident. Any benchmark of the splitting
conducted on a ring will conclude that it is worthless — which is precisely
the conclusion an earlier revision of this note drew, and recorded as a
rate-neutrality theorem. The lesson for open-line work is concrete: **a ring
is not a valid testbed for the reach condition**, because the one parameter
the condition depends on is frozen there at its worst value.

### 6.4 How far off is the ring parabola?

A Gaussian of width $`\sigma_x = 0.4`$ released at $`x = 1`$, integrated to
$`t = 2`$, under the exact symbol and under the Newtonian term alone:

| | periodised parabola | cosine well |
|---|---|---|
| $`L^1`$ gap at $`t = 2`$ | $`8.28\times10^{-5}`$ | $`7.44\times10^{-1}`$ |
| negativity generated | $`2.36\times10^{-5}`$ | — |

Nearly nine thousand times more Newtonian than a cosine well of comparable
depth, with a residual negativity generated entirely at the seam.

![Evolution comparison](https://raw.githubusercontent.com/billpage/wpmw/output/figures/compensated_evolution.png)

---

## 7. Coulomb

**Theorem C6.** For $`V(x) = -Z/x`$, with $`\rho = y/x`$,

```math
M \;=\; \frac{2iZy}{\hbar\thinspace(x^2 - y^2)},
\qquad
M_{\mathrm{cl}} \;=\; \frac{2iZy}{\hbar\thinspace x^2},
\qquad
\frac{M_{\mathrm{res}}}{M_{\mathrm{cl}}} \;=\; \frac{\rho^2}{1 - \rho^2} .
```

*Verified* to $`7.1\times10^{-15}`$ across radii and $`\rho`$.

Since $`M/M_{\mathrm{cl}} = (1-\rho^2)^{-1} = \sum_{n\ge0}\rho^{2n}`$, the
$`n`$-th Moyal term is exactly $`\rho^{2n}`$ times the classical one, all of
the same sign. **The Moyal series for Coulomb is geometric, and it converges
if and only if the reach does not touch the nucleus,** $`|y| < |x|`$. There is
no cancellation among orders to rescue a world whose reach spans the nucleus;
the terms all add. For Coulomb the reach condition of C4 is therefore not a
refinement but a necessity: outside it there is no expansion at all, and a
soft core or screening is required.

| $`\rho = y/x`$ | 0.05 | 0.10 | 0.25 | 0.50 | 0.90 |
|---|---|---|---|---|---|
| $`\lvert M_{\mathrm{res}}\rvert/\lvert M_{\mathrm{cl}}\rvert`$ | 0.0025 | 0.0101 | 0.0667 | 0.3333 | 4.2632 |

In state terms $`\rho \sim \hbar/(2\sigma_p|x|)`$, so the condition reads
$`\sigma_p|x| \gg \hbar/2`$ — the uncertainty product in units of $`\hbar`$.
Hydrogen's ground state sits at about $`1`$, giving $`\rho \sim 1/2`$ and a
residual fraction of a third; a Rydberg state of principal quantum number
$`n`$ sits at about $`n`$, with a fraction $`\sim 1/(4n^2)`$. The hop channel
also localises: $`|M_{\mathrm{res}}| \sim 1/x^4`$ against a classical force
falling as $`1/x^2`$, so after the split the quantum channel is a core effect
while the long-range Kepler dynamics is deterministic.

![Coulomb](https://raw.githubusercontent.com/billpage/wpmw/output/figures/compensated_coulomb.png)

---

## 8. Summary

| | statement | verified |
|---|---|---|
| C1 | $`e^{\tau M} = e^{\tau M_{\mathrm{cl}}}e^{\tau M_{\mathrm{res}}}`$, no Trotter error | $`\sim10^{-16}`$ |
| C2 | the residual is the odd part of the cubic Taylor remainder of $`V`$ | bound tight to $`<2\times`$ |
| C3 | at bounded reach the residual is a bounded, number- and momentum-conserving signed jump measure | $`\sim10^{-6}`$ of $`V'`$ |
| C4 | the split gains for $`k\thinspace y_{\max} \ll \pi`$; reach and momentum quantum are one parameter | $`TV`$ table, §4 |
| C5 | on a circle only constants have a vanishing residual | — |
| C6 | Coulomb: the Moyal series converges iff the reach misses the nucleus | $`7.1\times10^{-15}`$ |
| C7 | quiet region: $`V'''`$ vanishing on the reach implies no events at $`x`$ | edge sharp to $`10^{-15}`$ |
| C0 | the mode stencil is an identity; the derivative approximates it | analytic |

The short version. On the open line the classical force can be moved into the
first substep, exactly, for any potential and with no Trotter error, and what
is left for the jump channel is a zero-mean focus-and-hop interaction sourced
by $`V'''`$ and supported within one reach of the non-quadratic part of the
potential. The condition is that a world's coherence reach be short compared
with the scale on which the potential varies. A ring freezes that one
parameter at its worst value, which is why the reorganisation looks empty when
tested there; the ring's value is as a source of closed-form counterexamples,
not as a home for the ontology.

---

## 9. Open items

- **CLS1.** C2 bounds the residual by the supremum of $`|V'''|`$ over the
  reach, a crude bound. The natural refinement is the best quadratic fit to
  $`V`$ on the reach interval rather than the osculating quadratic at $`x`$;
  the residual would then be the fit error rather than the Taylor remainder.
  Does the fitted split still commute exactly, given that the fit is
  $`x`$-dependent?
- **CLS6.** C7 is sufficient but not necessary, the gap being potentials whose
  Taylor remainder is even on the reach without vanishing. Parity supplies
  isolated such points. Is there a potential with a quiet region of positive
  measure not explained by the vanishing of $`V'''`$?
- **CLS7.** §5.1 makes the coherence horizon the restorer of position-space
  locality, sharpening Theorem O1. Does the same truncation restore locality
  for the two-body kernel of the interworld coupling notes, or is that a
  different mechanism?
- **CLS2.** Does the deterministic acceleration have a field reading in the
  four-action model? The natural candidate is the mean occupation of a
  coherent mode field, with the residual as the fluctuation about it. If that
  holds, C3 is the mean-field / shot-noise decomposition of a single field.
  This wants checking against the momentum-and-energy balance argument of the
  Cyganski memo, which currently ties all of the force to the photon channel.
- **CLS3.** C4 identifies the reach with the momentum quantum. On the crystal
  lattice momentum is quantised and a deterministic sub-$`\Delta p`$
  acceleration cannot be represented at all. Does step 1 then require
  continuous momentum, and hence — via Theorem O5 — an aperiodic potential?
- **CLS4.** What sets $`L_c`$ physically? C4 makes it the parameter that
  decides whether the classical/quantum split exists, which is a heavy load
  for a quantity so far introduced as a truncation.
- **CLS5.** An errata line for
  `../supplement/phase_space_crystal_lattice_supplement.md` §6.1, which
  describes the mode finite difference as an approximation. For a single
  Fourier mode it is an identity, and that exactness is what lets the mode sum
  reproduce the whole Moyal series.

---

## Sources

- David Cyganski, *Bohm to Nelson and four action Wigner*, August 2026 (project memo).
- David Cyganski, *Extended Fokker–Planck Eq. and the QLE V2* (project memo).
- [`../supplement/phase_space_crystal_lattice_supplement.md`](../supplement/phase_space_crystal_lattice_supplement.md) §5–§7.
- [`open_position_space.md`](open_position_space.md) — Theorems O1 and O5, and the coherence horizon.
- [`four_rule_microdynamics_equivalence.md`](four_rule_microdynamics_equivalence.md).
- R. L. Hudson, "When is the Wigner quasi-probability density non-negative?", *Rep. Math. Phys.* **6**, 249 (1974).
