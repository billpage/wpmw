# The interworld coupling: why four rules, and why not more

**Status.** Analysis note, step 12 of the ladder. Companion demo:
`src/demo_interworld_coupling.py`. Prompted by D. Cyganski's question (Zoom,
August 2026): *why the four rules?* — together with his proposal that they
might follow from a non-linear interaction between world-particles, of the
kind that produces four-wave mixing for charged particles in a periodic field.

---

## 0. What this note inherits, supplies, corrects and concedes

**Inherits.** From
[`position_pair_ladder.md`](position_pair_ladder.md): the pair ontology — each
element of $`\rho`$ is a pair of legs, positon = ket leg, negaton = bra leg —
and the identification of momentum as a misalignment rather than a carried
quantity. From
[`phase_alignment_microdynamics.md`](phase_alignment_microdynamics.md): the
misalignment $`\mu = \Phi_{\rm ket} - \Phi_{\rm bra}`$ as the sole relational
datum of a pair. From
[`four_rule_microdynamics_equivalence.md`](four_rule_microdynamics_equivalence.md)
§1: the reduction of the four rules to two signed channels, and the exactness
family. From
[`open_position_space.md`](open_position_space.md) §2: the Wigner kernel in
separation form.

**Supplies** one thing that
[`sea_dressed_microdynamics.md`](sea_dressed_microdynamics.md) lists among its
postulates. That note assumes "the half-quantum stencil offsets" along with
the rate field and its sign structure. §3 below *derives* the half from
geometry alone: the two legs of a pair sit at $`X \pm Y/2`$, so a potential
mode of wavevector $`k`$ is sampled at $`kY/2`$ in the separation, and the
conjugate momentum shift is $`\hbar k/2`$. No dynamical input is used.

**Corrects nothing in the repository.** Every statement below is consistent
with what is already written; the note is a change of variables, not a
revision. What it does correct is an expectation formed during the discussion
that prompted it — see the concession below.

**Contradicts the hypothesis that prompted the note.** The four rules do not
require, and are not compatible with, a non-linear interworld force law of the
four-wave-mixing kind. The coupling is *linear* in $`V`$ (§1), *separable* —
a difference of one-leg terms, not a two-body potential (§1) — and it does no
work at all (§6). The four-wave-mixing intuition nevertheless lands on the
right structure for a reason worth recording: the process does need a pump,
because the coupling vanishes identically when $`V`$ does, and its resonance
geometry is fixed by a grating, because the coupling is periodic in the leg
separation. What plays the part of the non-linear medium is the *bilinearity
of* $`\rho = \psi\psi^{\dagger}`$, not a non-linear force.

**Concedes, up front, what this note is not.** Three levels of claim are mixed
below and they should not be confused.

1. **Theorem.** The momentum-transfer channels available at pair midpoint
   $`X`$ are exactly the Fourier spectrum of the coupling in the leg
   separation (Theorem I3). This is rigorous and verified.
2. **Restatement.** "Periodic in the separation" and "discrete in momentum"
   are Fourier-dual descriptions of one fact. §3 does not discover an
   independent physical requirement; it re-expresses the discreteness in the
   conjugate variable. This is genuinely useful — a periodic coupling is far
   easier to reason about than a jump stencil — but it is not new information.
3. **Postulate.** That $`Y`$ is a *physical separation between two co-existing
   world-particles* is an interpretive commitment, not a derived fact.
   Mathematically $`Y`$ is the off-diagonal coordinate of a bilinear object.
   §7 states what the commitment costs.

---

## 1. Tutorial: the only way the potential can enter

Write $`x_1`$ for the ket leg, $`x_2`$ for the bra leg. The von Neumann
equation in the position representation is

```math
i\hbar\thinspace\partial_t \rho(x_1, x_2)
  = -\frac{\hbar^2}{2m}\bigl(\partial_{x_1}^2 - \partial_{x_2}^2\bigr)\rho
    + \bigl[V(x_1) - V(x_2)\bigr]\thinspace\rho(x_1, x_2).
```

The potential term is **multiplicative**: no derivatives, local in
$`(x_1, x_2)`$. It never sees $`V(x_1)`$ or $`V(x_2)`$ alone. Define

```math
U(x_1, x_2) \;=\; V(x_1) - V(x_2).
```

**Proposition I1 (structure of the coupling).** $`U`$ has three properties,
each immediate and each independent of any dynamical assumption:

1. $`U(x, x) = 0`$ — coincident legs do not interact.
2. $`U(x_2, x_1) = -U(x_1, x_2)`$ — the coupling is antisymmetric under leg
   exchange.
3. $`V \equiv 0 \Rightarrow U \equiv 0`$ — free legs do not interact.

Property 3 is worth pausing on, because it is the sharpest available test of
any proposed interworld force law. A free wavepacket spreads, and in the
Wigner picture that spreading is *exact classical shearing* of $`W`$ with no
collision term whatsoever. Any model in which worlds carry a mutual potential
that survives switching off the external field predicts something extra
happening to a free particle, and is thereby falsified without further
calculation. This is a genuine point of difference from the Hall–Deckert–
Wiseman many-interacting-worlds construction, in which the interworld
potential is present for a free particle and *is* what produces the spreading.
The two pictures put the quantum content in different places.

Property 3 also disposes of a natural but wrong reading of §1's coupling as a
pair potential. $`U`$ is *separable* — a difference of one-leg terms — so it
is not a two-body force in the mechanical sense at all. No search over pair
potentials $`U(x_1 - x_2)`$ was ever going to find it.

**Two analogies for what follows.**

*The diffraction grating.* A structure periodic in space produces discrete
diffraction orders in momentum, spaced by $`\hbar`$ times the grating
wavevector. §3 is that statement with "space" replaced by "leg separation".
The converse is what matters here: a finite rule set *requires* a grating.
Any interpretation that gives worlds a smooth, decaying mutual potential is
committed to a continuum of rules, not four.

*The see-saw over corrugated ground.* Two people stand symmetrically on a
plank and walk apart. The ground ripples with wavelength $`\lambda`$. The
plank's *tilt* depends on the height difference at the two ends; but walking
apart by $`Y`$ moves each person only $`Y/2`$, so the tilt completes one cycle
only when $`Y = 2\lambda`$. The tilt-grating has twice the period of the
ground, hence half the wavevector. Every factor of one-half in this project
traces back to that one geometric fact.

---

## 2. Midpoint and separation

Set

```math
X = \tfrac{1}{2}(x_1 + x_2), \qquad Y = x_1 - x_2,
\qquad\text{so}\qquad x_1 = X + \tfrac{Y}{2}, \quad x_2 = X - \tfrac{Y}{2},
```

and write the coupling in these variables:

```math
U(X, Y) \;=\; V\bigl(X + \tfrac{Y}{2}\bigr) - V\bigl(X - \tfrac{Y}{2}\bigr).
```

**Convention note.** [`open_position_space.md`](open_position_space.md) §2
writes the Wigner kernel with the *half* separation $`y = Y/2`$, so a profile
appearing there as $`\sin(ky)`$ appears here as $`\sin(kY/2)`$. This note uses
the full separation $`Y`$ throughout, because $`Y`$ is the variable whose
conjugate is $`p`$ and because it is the one a two-particle reading would call
a distance. The physical content is identical.

**Proposition I2 (the grating).** For the single-mode reference potential
$`V(x) = -V_p\cos kx`$,

```math
U(X, Y) \;=\; 2 V_p\thinspace \sin(kX)\thinspace \sin\!\bigl(\tfrac{kY}{2}\bigr),
```

a **midpoint amplitude** times a **separation grating**. The grating has
wavevector $`k/2`$, i.e. period $`2\lambda`$ — twice the period of $`V`$. The
midpoint amplitude is the local classical force up to a fixed factor:

```math
\Gamma(X) \;=\; \frac{V_p}{\hbar}\sin(kX) \;=\; -\frac{F(X)}{\hbar k},
\qquad F = -V'.
```

*Proof.* Sum-to-product on the cosine; the rate identity is one line of
differentiation. Verified numerically to $`1.2\times10^{-14}`$ and
$`2.2\times10^{-16}`$ respectively (§8, part B), with the period measured from
zero crossings rather than asserted: $`15.9998`$ against $`\lambda = 8`$.
$`\square`$

For a general potential $`V(x) = \sum_k \tilde V_k\thinspace e^{ikx}`$,

```math
U(X, Y) \;=\; \sum_k \tilde V_k\thinspace e^{ikX}
   \cdot 2i \thinspace\sin\!\bigl(\tfrac{kY}{2}\bigr).
\qquad\text{(I.1)}
```

**Each Fourier mode of the potential contributes exactly one sinusoid in the
leg separation, at half its wavevector.** That is the whole content of §3,
already visible.

---

## 3. Theorem I3: channels are the separation spectrum

The Wigner transform is a Fourier transform in the leg separation,

```math
W(X, p) \;=\; \frac{1}{2\pi\hbar}\int dY\thinspace e^{-ipY/\hbar}\thinspace
   \rho\bigl(X + \tfrac{Y}{2},\thinspace X - \tfrac{Y}{2}\bigr),
```

under which multiplication by $`e^{i\kappa Y}`$ is a rigid shift of $`p`$ by
$`\hbar\kappa`$.

**Theorem I3 (channel theorem).** Expand the coupling in its separation
harmonics at fixed midpoint,
$`U(X, Y) = \sum_n c_n(X)\thinspace e^{i\kappa_n Y}`$. Then the potential term
of the Wigner equation is exactly

```math
\Theta_V[W](X, p) \;=\; \frac{-i}{\hbar}\sum_n c_n(X)\thinspace
   W\bigl(X,\thinspace p - \hbar\kappa_n\bigr),
```

so the set of available momentum-transfer channels is precisely the Fourier
spectrum of $`U(X, \cdot)`$, with channel rates $`c_n(X)/\hbar`$. The channel
set is **discrete if and only if the coupling is periodic in $`Y`$**, which
holds if and only if $`V`$ is periodic in $`x`$. $`\blacksquare`$

Applying (I.1): mode $`k`$ contributes the two harmonics
$`\kappa = \pm k/2`$, hence the two shifts $`\Delta p = \pm\hbar k/2`$ with
rates $`\pm\Gamma(X)`$, hence the stencil

```math
\Theta_V[W] \;=\; \Gamma(X)\thinspace
  \bigl[\,W\bigl(X, p + \tfrac{\hbar k}{2}\bigr)
        - W\bigl(X, p - \tfrac{\hbar k}{2}\bigr)\,\bigr],
```

which is the spec's mode-$`q`$ stencil with $`q\Delta p = \hbar k/2`$. This is
verified against the *full* Wigner operator — all Moyal orders, no truncation
— at $`1.7\times10^{-15}`$ for one mode and $`4.9\times10^{-15}`$ for three
(§8, part C). Modes do not mix: the channel set of a multi-mode potential is
the union of the per-mode pairs.

**Corollary I3.1 (rule counting).** A potential with $`M`$ active Fourier
modes has $`2M`$ signed shift channels, hence — each signed channel splitting
into a forward and a reverse rule — $`4M`$ rules. **A one-mode potential has
exactly four.** A non-periodic potential has a continuum and no finite rule
set exists.

That is the answer to the question that prompted this note. *Why four rules?
Because the reference potential has one Fourier mode.* The four rules are not
fundamental; they are the harmonic content of a single-mode potential read in
the separation variable.

**The control matters here.** Properties I1(1) and I1(2) — vanishing at
coincidence, antisymmetry under exchange — are *not* enough. An
antisymmetrised Gaussian
$`g(Y) = e^{-(Y-b)^2/w} - e^{-(Y+b)^2/w}`$ has both and yields a continuum.
§8 part C separates the two by a scaling test rather than by eye: widening the
$`Y`$-window narrows the cosine line like $`1/Y_{\max}`$ (zero intrinsic
width, a genuine line) while the Gaussian's width saturates at $`\approx 1.41`$
(genuine continuum).

---

## 4. Theorem I4: the Moyal series is the separation expansion

Because $`U`$ is odd in $`Y`$, its Taylor expansion contains only odd powers:

```math
U(X, Y) \;=\; V'(X)\thinspace Y \;+\; \tfrac{1}{24}V'''(X)\thinspace Y^3
  \;+\; \tfrac{1}{1920}V^{(5)}(X)\thinspace Y^5 \;+\; \cdots
```

and under the transform $`Y \leftrightarrow i\hbar\thinspace\partial_p`$.

**Theorem I4.** The $`Y^{2l+1}`$ term of the coupling is the $`l`$-th term of
the Moyal series. In particular the $`Y^1`$ term is exactly the classical
force term $`V'(X)\thinspace\partial_p W`$, so **a coupling that is linear in
the leg separation generates classical Liouville flow exactly**. $`\square`$

Three consequences, all verified (§8, part D) on a cat-like state with genuine
negativity, $`\min W = -0.3407`$:

| $`V`$ | $`U(X,Y)`$ | consequence |
|---|---|---|
| linear, $`ax`$ | $`aY`$ | exactly classical, residual $`2.0\times10^{-16}`$ |
| quadratic, $`\tfrac{c}{2}x^2`$ | $`cXY`$ | exactly classical, residual $`2.4\times10^{-16}`$ |
| cubic, $`bx^3`$ | $`b\bigl(3X^2Y + \tfrac{Y^3}{4}\bigr)`$ | one correction, then exact: $`3.1\times10^{-16}`$ |
| periodic | $`\sum_k \propto \sin(kY/2)`$ | whole series resums into finite shifts |

The quadratic row is the one with teeth. $`U = cXY`$ holds for **either sign
of $`c`$**, so the harmonic *and* the inverted harmonic oscillator alike have
a coupling linear in the leg separation and therefore **no jump channel at
all**. This is the separation-space explanation of the result recorded in
[`fourd_microdynamics.md`](fourd_microdynamics.md) §12 (the Moyal bracket
truncates for the harmonic potential) and it carries a direct consequence for
the project's current test problem. The inverted pair barrier of
[`../supplement/inverted_pair_barrier.md`](../supplement/inverted_pair_barrier.md)
is quadratic in the relative coordinate, so its two-body Wigner equation is
exactly the classical Liouville equation in four-dimensional phase space. It
remains an excellent test of the *transport* half of the algorithm — the
quantum content lives entirely in the initial condition — but it cannot test
the four rules, because under it there are none.

---

## 5. What a genuine pair interaction would look like, for contrast

The separable coupling of §1 is what a *one-body* potential looks like in pair
coordinates. A *real* two-body potential — two physical particles, not two
legs — behaves differently, and the difference is worth having on record since
it is where the four-wave-mixing analogy is at its strongest.

**Proposition I5a.** For two particles with pair potential
$`U_{\rm pair}(x_1 - x_2) = \sum_k \tilde U_k\thinspace e^{ik(x_1-x_2)}`$, the
exact two-body Wigner potential term is

```math
\Theta = \frac{i}{\hbar}\sum_k \tilde U_k\thinspace e^{ik(x_1-x_2)}
  \Bigl[\,W\bigl(p_1 - \tfrac{\hbar k}{2},\thinspace p_2 + \tfrac{\hbar k}{2}\bigr)
       - W\bigl(p_1 + \tfrac{\hbar k}{2},\thinspace p_2 - \tfrac{\hbar k}{2}\bigr)\Bigr],
```

an anti-correlated, total-momentum-conserving pair shift — that is, a
focus/defocus channel and nothing else. External modes give hops; pair modes
give focus/defocus. $`\square`$

This sharpens the correction recorded in
[`fourd_microdynamics.md`](fourd_microdynamics.md) §0, which noted that the
four-rule note's claim about hops not being particle–particle exchanges holds
for external modes only. Here is the positive half of that statement.

Note the asymmetry it exposes. In the *pair-of-legs* reading of §1 there is no
genuine two-body force; in the *two-physical-particles* reading of I5a there
is. The four rules of this project live in the first reading. Any programme
that hopes to obtain them from a real interparticle force is working in the
second, and inherits a channel that vanishes only when the pair coupling does,
not when the external field does — which contradicts I1(3).

---

## 6. Theorem I5: the coupling winds, it does not push

**Theorem I5.** Under the potential term alone, $`|\rho(x_1,x_2)|`$ is
constant and

```math
\frac{d\mu}{dt} \;=\; -\frac{U(x_1, x_2)}{\hbar},
\qquad \mu = \arg\rho = \Phi_{\rm ket} - \Phi_{\rm bra}.
```

*Proof.* The potential term integrates to multiplication by
$`e^{-iUt/\hbar}`$, a pure phase. Verified: modulus change
$`8.4\times10^{-14}`$, phase-advance error $`2.9\times10^{-15}`$ (§8, part E).
$`\square`$

So the interworld coupling exerts no force and does no work. It is a
**differential winding rate** for the pair's clock misalignment — precisely
the $`\mu`$ of [`phase_alignment_microdynamics.md`](phase_alignment_microdynamics.md),
whose equation of motion this supplies in one line from the potential.

Momentum is not a property a leg carries. It is the rate at which $`\mu`$
winds per unit separation, $`p = \hbar\thinspace\partial\mu/\partial Y`$,
which is the continuum form of the lattice identity $`\bar p = \hbar\mu/a`$ of
[`position_pair_ladder.md`](position_pair_ladder.md) Theorem P2. A "momentum
jump of $`\hbar k/2`$" is therefore a change in how fast the phase winds
across the pair — and since the coupling's winding profile is a sinusoid in
$`Y`$, that change comes in exactly two discrete sizes. §3 and §6 are the same
statement told in the two conjugate variables.

---

## 7. What the reading costs

The mathematics of §§1–6 is representation-independent. Reading it as an
*interpretation* — worlds that interact — requires a choice, and the two
available choices are not equally hospitable.

| | **A. legs are positions** | **B. worlds are phase-space points** |
|---|---|---|
| what a world is | a place plus a clock phase | a sample $`(X,p)`$ of $`W`$ |
| the carrier | a **pair** of legs | a single particle |
| is there an interaction? | yes: $`U`$, periodic in $`Y`$ | no — $`Y`$ has been integrated out |
| momentum | emergent, $`\hbar\partial_Y\mu`$ | carried |
| sign problem | in the *phase*, which is allowed | in the *rate*, which is not |
| the no-go bites? | no | yes |

Reading B is the one the orthogonality/permutation no-go theorem applies to:
the generator has off-diagonal entries $`+\Gamma`$ and $`-\Gamma`$, so it is
not a Markov generator and no ensemble of independent particles realises it.
Reading A escapes because $`U`$ is an ordinary real, signed potential and the
sign lives in a phase.

The price of Reading A is the postulate flagged in §0: $`Y`$ is a difference
of two arguments of a bilinear object, and calling it a separation between two
co-existing entities is a substantive commitment. It entails that the ensemble
does **not** sample $`W`$ at all — it samples $`\rho`$ — and that a "world" is
not a particle with a momentum but a position carrying a clock. Everything in
[`position_pair_ladder.md`](position_pair_ladder.md) is already committed to
this; the present note simply makes the cost explicit and shows what is bought
with it.

---

## 8. Numerical results

All figures from `src/demo_interworld_coupling.py`, `HBAR = 1`, single-mode
reference potential $`V = -V_p\cos kx`$ with $`V_p = 1.5`$,
$`\lambda = 8`$, $`k = 2\pi/\lambda`$.

**A. The difference structure.** $`\max|U(X,0)|`$, $`\max|U(X,Y)+U(X,-Y)|`$
and $`\max|U|`$ at $`V \equiv 0`$ are all $`0.0`$ exactly, over $`4\times10^4`$
random leg pairs. The antisymmetrised-Gaussian control shares both symmetry
properties to $`0.0`$, confirming that they do not discriminate.

**B. The grating.**

| quantity | value |
|---|---|
| $`\max\bigl|U - 2V_p\sin kX\sin(kY/2)\bigr|`$ | $`1.221\times10^{-14}`$ |
| coupling period in $`Y`$, measured from zeros | $`15.9998`$ |
| potential period $`\lambda`$ | $`8.0000`$ |
| ratio | $`1.999980`$ |
| $`\max\bigl|\Gamma(X) + F(X)/(\hbar k)\bigr|`$ | $`2.220\times10^{-16}`$ |

**C. Theorem I3.** Two lines found, at $`q = \pm 0.39270`$; $`\hbar k/2 =
0.39270`$; discrepancy $`0.0`$. Discrete-versus-continuum scaling test:

| $`Y`$ window | line width, cosine mode | line width, antisym. Gaussian |
|---|---|---|
| $`50`$ | $`0.06283`$ | $`1.38230`$ |
| $`100`$ | $`0.03142`$ | $`1.38230`$ |
| $`200`$ | $`0.03142`$ | $`1.41372`$ |
| $`400`$ | $`0.01571`$ | $`1.41372`$ |

Against the full Wigner operator:
$`\max|\Theta_{\rm exact} - \text{two-line stencil}| = 1.665\times10^{-15}`$
(one mode, $`\max|\Theta| = 0.6905`$) and $`4.940\times10^{-15}`$ (three
modes, $`\max|\Theta| = 0.9222`$).

**D. Theorem I4.** Relative residual after truncating the $`Y`$-expansion:

| potential | truncation | relative residual |
|---|---|---|
| $`0.7x`$ | $`Y^1`$ | $`1.999\times10^{-16}`$ |
| $`+x^2/2`$ | $`Y^1`$ | $`2.387\times10^{-16}`$ |
| $`-x^2/2`$ | $`Y^1`$ | $`2.387\times10^{-16}`$ |
| $`0.02x^3`$ | $`Y^1`$ only | $`1.000`$ |
| $`0.02x^3`$ | $`Y^1 + Y^3`$ | $`3.052\times10^{-16}`$ |

**E. Theorem I5.** After $`t = 0.4`$ under the potential term alone:
$`\max\bigl||\rho| - |\rho_0|\bigr| = 8.438\times10^{-14}`$ and
$`\max\bigl|\Delta\arg\rho + Ut/\hbar\bigr| = 2.887\times10^{-15}`$.
Rule counting: one cosine mode $`\to`$ 2 lines, 4 rules; three cosine modes
$`\to`$ 6 lines, 12 rules; antisymmetrised Gaussian $`\to`$ continuum.

![The interworld coupling in midpoint and separation coordinates](https://raw.githubusercontent.com/billpage/wpmw/output/figures/interworld_coupling.png)

Panel (a): a pair straddling the potential; only the height *difference* at
the two legs enters. Panel (b): the coupling as a function of leg separation,
against the potential's own profile — the grating has twice the period.
Panel (c): the separation spectrum, two lines at $`\pm\hbar k/2`$, against the
antisymmetrised-Gaussian control's continuum. Panel (d): the $`U(X,Y)`$
landscape, odd in $`Y`$ and odd in $`X`$ about each well, with the nodal line
$`Y = 0`$ marked.

---

## 9. Consequences and open items

**Consequences.**

1. *Why four rules* has an answer that requires no new postulate: the
   reference potential has one Fourier mode, and each mode contributes exactly
   two signed shifts at $`\pm\hbar k/2`$ (Corollary I3.1).
2. The half-quantum offset postulated in
   [`sea_dressed_microdynamics.md`](sea_dressed_microdynamics.md) follows from
   the pair geometry alone: legs at $`X \pm Y/2`$.
3. The inverted pair barrier cannot serve as a test of the four rules, only of
   transport (§4).
4. The four rules are exact and finite only for periodic $`V`$. For a Gaussian
   or Coulomb potential the rule set is a continuum, and "worlds obey four
   rules" is simply false there.

**Open items.**

1. **Is $`U`$ unique?** §1 shows $`V(x_1) - V(x_2)`$ *works*. It does not show
   it is the only exchange-antisymmetric coupling vanishing at $`Y = 0`$ that
   reproduces the QLE. Any $`\Delta U(X,Y)`$ whose $`Y`$-spectrum vanishes
   would be invisible; whether a nonzero such $`\Delta U`$ exists determines
   whether Reading A has hidden freedom. Self-contained and probably
   answerable.
2. **Pair stability.** Nothing here says which leg pairs with which. If pairs
   can re-partner between vertices, $`\mu`$ is not well defined across a
   vertex, and Theorem I5's winding law has no subject. Postulate (S) of
   [`relational_pairing_and_carrier_lock.md`](relational_pairing_and_carrier_lock.md)
   now looks load-bearing for the *interpretation*, not merely for the
   algorithm's cost.
3. **The continuum case.** What replaces the four rules for non-periodic
   $`V`$, given that the ensemble reading must still work there? The coherence
   horizon of [`open_position_space.md`](open_position_space.md) §3 is the
   obvious candidate — it bounds $`|Y|`$, which is exactly the variable this
   note is about — but the resulting rule set is approximate, and by how much
   is not established.
4. **Third-moment test for a mechanical model.** Moments $`0`$–$`2`$ of the
   stencil are independent of the shift size at fixed force, so classical
   mechanics constrains only the product $`\Gamma\Delta`$. The shift first
   appears at the third momentum moment, as $`F\Delta^2`$. Any candidate
   mechanical microdynamics should be compared there, not at the level of the
   force law.
