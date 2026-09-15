# Soft-core Coulomb: four lobes, a dark nucleus, and a ceiling that moves

> The potential open item K-LS5 nominated, put through the geometry the Eckart note put `V0 sech^2(r/a)` through. Theorem Z1: for `-Z/sqrt(r^2 + eps^2)` the third derivative is `(Z/eps^4) u (6u^2 - 9) (1+u^2)^(-7/2)` with `u = r/eps`, so it vanishes at `r = 0` and `r = +- eps sqrt(3/2)` and the soft core carries **four** emission lobes — the same count K7 found for sech², but with a quiet radius set by the softening length alone and carrying no `Z` at all. Pure Coulomb has `V''' = -6Z/r^4`, one sign on each half-line and no interior zero, so it carries two: the inner pair of lobes is manufactured by the softening and collapses onto the origin as `eps -> 0`, which makes the soft core a source of structure rather than a regularisation convenience. Theorem Z2: `Gamma(0) = 0` identically at every reach, so the nucleus is dark exactly as the sech² summit is, and the interior quiet ring survives the horizon, drifting outward by 0.4 per cent at `y_max = eps/4` and 6.1 per cent at `y_max = 0.99 eps`. **Theorem Z3 is what this potential can raise and sech² cannot**: by Corollary K1.2 the reach ceiling is `R(x) = sqrt(x^2 + eps^2)`, which *varies with position*, so either the momentum quantum `dp = pi hbar / 2 y_max` varies with `x` — and the phase-space crystal is not uniform — or a single uniform lattice must take the infimum `y_max < eps`, the tightest ceiling in the problem and located at the one point where nothing is emitted anyway. Theorem Z4 prices the uniform choice: the well is harmonic for `eps >> a0` with `sigma_r = eps^(3/4)/sqrt2`, so resolving the ground state to `k` rungs per `sigma_p` needs `eps >= k^4 pi^4 / 4` in units `a0 = hbar^2 / mu Z` — 24.35 `a0` for one rung, 389.6 for two, the cost quartic in resolution; measured crossings 34.6 and 432.9, the ratio to prediction falling from 1.42 to 1.11 as the anharmonic correction dies. An unsoftened atom therefore cannot live on a uniform reach-limited crystal at all. Theorem Z5 is the consolation and the reason this is the right vehicle for S-SP6: at threshold the horizon spans 2.8 `sigma_r` and `Gamma(sigma_r)/Gamma_max = 0.974`, so unlike the Eckart barrier — where K-LS2's ceiling and resolution conditions have no common ground — the window where the lattice works is *not* the window where the demographic channel is empty. Leaves the non-uniform lattice as Z-LS1 and the transmission calculation, which has no closed form and needs a split-operator reference, as Z-LS2.

*Ladder abstract — see the [full list](README.md#the-ladder).*

**Status.** Analysis note, step 19 of the ladder. Companion demo:
`src/demo_soft_core_coulomb.py`. Prompted by open item **K-LS5** of
[`eckart_barrier_compensated.md`](eckart_barrier_compensated.md) ("soft-core
Coulomb has a tunable ceiling $`\epsilon`$ and, unlike the Eckart barrier, an
attractive well with bound states") and by **S-SP6** of
[`sea_population_equilibrium.md`](sea_population_equilibrium.md), which wants a
stationary state to make "equilibrium" exact rather than asymptotic.

---

## 0. What this note inherits, and why this potential

### 0.1 Inherited

- **Theorem K1** (`eckart_barrier_compensated.md` §2): the reach ceiling is
  the distance from $`x`$ to the nearest *complex* singularity of $`V`$, and
  **Corollary K1.2**, that for $`-Z/\sqrt{z^2+\epsilon^2}`$ the branch points
  sit at $`z = \pm i\epsilon`$ so $`R(x) = \sqrt{x^2+\epsilon^2}`$ — verified
  there at $`x = 0.4`$ to three figures for three values of $`\epsilon`$.
- **Theorem C6** (`compensated_liouville_splitting.md` §7): for pure Coulomb
  the Moyal series is geometric and converges iff the reach misses the
  nucleus.
- **Theorem C7**: a world at $`x`$ takes no events if $`V'''`$ vanishes across
  its whole reach.
- **Theorem K7**: for sech² the emission rate vanishes and $`K_q`$ reverses
  sign at $`r = 0`$ and $`r = \pm a\thinspace\mathrm{artanh}\sqrt{2/3}`$, so
  the barrier carries four alternating lobes and the summit is quiet.
- **Theorem I4** (`interworld_coupling.md`): for quadratic $`V`$ the jump
  channel is empty. The relevant limit here, since a deep soft core is a
  harmonic well.
- **Postulate (S)** and **Proposition K8**: no world-particle changes momentum
  discontinuously, and the pair an event produces is ionised from a bound sea
  pair co-located with the parent.

### 0.2 Why this potential

The Eckart barrier was chosen because it is bounded, asymptotically free, has
non-vanishing $`V'''`$, and has an exact transmission coefficient. Soft-core
Coulomb is chosen for the three things it has that sech² does not.

1. **A ceiling that varies with position.** sech²'s ceiling is the uniform
   $`\pi a/2`$ — a single number for the whole line. Coulomb's is
   $`\sqrt{x^2+\epsilon^2}`$. Every statement the project has made about "the
   reach" has implicitly assumed a global reach, and this is the first
   potential that makes the assumption visible. §3.
2. **Bound states.** The Eckart barrier has none, so every ledger measurement
   made on it is a transient. S-SP6 asks for a stationary state and this is
   the natural one.
3. **A tunable singularity.** $`\epsilon`$ interpolates between the true
   Coulomb pole, where C6 says the reach must miss the nucleus, and a harmonic
   well, where I4 says there is no channel at all. The interpolation is the
   subject of §4 and §5.

Throughout, $`\hbar = \mu = 1`$ and lengths are in
$`a_0 = \hbar^2/\mu Z`$, so $`Z = 1`$ and $`\epsilon`$ is measured in Bohr
radii of the problem. The attractive sign is used; every result in §1–§3 is
invariant under $`Z \to -Z`$ up to an overall sign of $`V'''`$, as Part A
verifies.

### 0.3 What this note claims

Z1 and Z2 are the sech² geometry redone and are unsurprising in kind, though
the closed form is cleaner. Z3 is new structure and is the item worth arguing
about. Z4 and Z5 are the price of avoiding Z3, and they come out the opposite
way from the Eckart barrier, which is the note's most useful result.

---

## 1. Theorem Z1: the quiet points

**Theorem Z1.** For $`V(r) = -Z/\sqrt{r^2+\epsilon^2}`$, with $`u = r/\epsilon`$,

```math
V'''(r) \;=\; \frac{Z}{\epsilon^{4}}\;
u\,\bigl(6u^{2} - 9\bigr)\,\bigl(1 + u^{2}\bigr)^{-7/2} ,
```

so $`V'''`$ vanishes exactly at

```math
r = 0
\qquad\text{and}\qquad
r = \pm\,\epsilon\sqrt{3/2} \;=\; \pm\,1.224745\,\epsilon ,
```

and the potential carries four emission lobes of alternating sign.

*Proof.* Write $`V = -(Z/\epsilon) g(u)`$ with $`g = (1+u^2)^{-1/2}`$. Then
$`g' = -u(1+u^2)^{-3/2}`$,
$`g'' = (2u^2-1)(1+u^2)^{-5/2}`$, and

```math
g''' = 4u(1+u^{2})^{-5/2} - 5u(2u^{2}-1)(1+u^{2})^{-7/2}
     = u\,(9 - 6u^{2})\,(1+u^{2})^{-7/2} ,
```

whence $`V''' = -(Z/\epsilon^4) g'''`$ as stated. The roots of
$`u(6u^2-9)`$ are $`u = 0`$ and $`u^2 = 3/2`$. By Theorem C7 these are quiet
points at short reach, and $`K_q`$ inherits the sign of $`V'''`$ at leading
order in the reach. $`\square`$

Part A checks the closed form against a fourth-order finite difference —
relative error $`1.5`$ to $`2.1\times10^{-3}`$, the difference stencil's own
truncation — and bisects the positive root, for $`Z = 1, 2, -1`$ and
$`\epsilon = 0.3, 0.7, 1.0, 2.5`$:

```
      Z      eps    rel err vs finite difference     root        eps sqrt(3/2)
      1.0    1.00           1.479e-03              1.224745       1.224745
      1.0    0.30           2.058e-03              0.367423       0.367423
      1.0    2.50           1.957e-03              3.061862       3.061862
      2.0    0.70           1.937e-03              0.857321       0.857321
     -1.0    1.00           1.479e-03              1.224745       1.224745
```

Six figures in every case, and — note the last two rows — the quiet radius
does not move when the charge is doubled or when the sign is reversed. It is
$`\epsilon\sqrt{3/2}`$ and nothing else. The sech² analogue,
$`a\thinspace\mathrm{artanh}\sqrt{2/3} = 1.146216\,a`$, is the same kind of
number for the same kind of reason.

**Corollary Z1.1 (the soft core creates the inner lobes).** For pure Coulomb,
$`V = Z/r`$ on each half-line gives $`V''' = -6Z/r^4`$, of one sign
throughout and with no interior zero. So pure Coulomb carries **two** lobes,
one per half-line, and softening manufactures the inner pair, which collapses
onto the origin as $`\epsilon \to 0`$.

This is worth more than it looks. The project's standing motivation for a soft
core, from Corollary K1.2, is that the softening length *is* the largest
coherence reach the potential will support — a good reason, but a
representational one. Z1.1 is an ontological reason: the lobe structure, which
by K7 is what organises the pair-emission geometry, does not exist at all
without a core. The core is not a numerical convenience that the physics
survives; it is where two of the four lobes come from.

![Lobes, ceiling, threshold](https://raw.githubusercontent.com/billpage/wpmw/output/figures/soft_core_coulomb_geometry.png)

(a) $`V'''`$ and the emission rate $`\Gamma`$ at $`y_{\max} = \epsilon/2`$,
both normalised, with the three quiet points marked; (b) the local momentum
quantum against position, against the floor a uniform lattice must accept;
(c) rungs per $`\sigma_p`$ against $`\epsilon`$, with the Z4 thresholds.

---

## 2. Theorem Z2: the nucleus is dark, and the ring survives

**Theorem Z2.** $`\Gamma(0) = \sum_{q\neq0}|K_q(0)| = 0`$ at every reach.

*Proof.* $`V`$ is even, so $`D_{\rm res}(0, y) = V(y) - V(-y) - 2yV'(0)`$
vanishes identically in $`y`$. Every Fourier coefficient of the zero function
is zero. $`\square`$

Part B, at $`Z = \epsilon = 1`$ with the raised-cosine horizon:

```
   y_max/eps   Gamma_max     Gamma(0)    Gamma(ring)/Gamma_max   ring radius   shift
       0.10  6.3194e-05     0.00e+00           2.062e-04       1.22575   +0.08%
       0.25  9.7290e-04     0.00e+00           9.644e-04       1.23018   +0.44%
       0.50  7.3955e-03     0.00e+00           3.770e-03       1.24431   +1.60%
       0.75  2.3049e-02     0.00e+00           8.757e-03       1.26815   +3.54%
       0.90  3.7551e-02     0.00e+00           1.306e-02       1.28690   +5.07%
       0.99  4.8090e-02     0.00e+00           1.587e-02       1.29981   +6.13%
```

Three readings.

**The nucleus is dark.** Exactly zero, not small — and for the same parity
reason the sech² summit is quiet. So the most violent point of the potential
emits nothing at all, and the emission is concentrated on an annulus around
it. Whatever the nucleus does to a world-particle, it does through the
classical force in the deterministic step and not through the demographic
channel.

**The interior ring survives the horizon** and drifts outward with it, by 0.4
per cent at quarter reach and 6.1 per cent at full reach, never washing out.
So the four-lobe structure of Z1 is a property of the regulated theory and not
only of the $`y_{\max}\to0`$ limit in which it was derived.

**It does not go all the way to zero.** $`\Gamma`$ at the ring is $`10^{-4}`$
to $`1.6\times10^{-2}`$ of its peak, rising with the reach. At finite reach the
quiet point is a minimum, not a node, because $`K_q`$ integrates $`V'''`$ over
the reach rather than sampling it. The same must be true of K7's sech² nodes,
which were reported as exact zeros of $`V'''`$ without this qualification;
that is not an error there but it is a gloss, and **K-LS7** — does the sea
inherit the lobes, with exact holes? — should be read with it in mind. The
holes are not exact even before transport fills them.

---

## 3. Theorem Z3: a ceiling that moves

**Theorem Z3.** By Corollary K1.2 the reach ceiling for soft-core Coulomb is
$`R(x) = \sqrt{x^2+\epsilon^2}`$, which is not constant. Hence the momentum
quantum implied by Theorem C4, $`\Delta p = \pi\hbar/2y_{\max}`$, cannot be
both position-independent and everywhere maximal.

Part C, at $`\epsilon = 1`$:

```
         x      R(x)     local dp    dp(x)/dp(0)
       0.00    1.0000     1.57080        1.00000
       1.00    1.4142     1.11072        0.70711
       4.00    4.1231     0.38097        0.24254
      16.00   16.0312     0.09798        0.06238
      64.00   64.0078     0.02454        0.01562
```

Two readings, and the whole of §4 exists because they are not equivalent.

**(i) A local reach.** Let $`y_{\max}(x) = R(x)`$, so each world consults the
potential out to its own analyticity limit. Then the momentum quantum is
coarse at the nucleus and fine in the far field, falling as $`1/|x|`$. The
phase-space crystal is then **not uniform**: the momentum lattice spacing is a
function of position.

This is not obviously wrong, and it has a reading that is almost attractive —
a deeply bound world-particle explores large momentum quanta, a free one
small. But it collides with a great deal of the existing ladder. Every
statement of the form "$`\Delta p`$ and $`y_{\max}`$ are the same parameter"
(Theorem C4), every rung-counting argument (Corollary K1.1's "fewer than
$`\beta`$ rungs span the barrier"), and the crystal-lattice construction
itself in
[`../algorithm/phase_space_crystal_lattice_algorithm.md`](../algorithm/phase_space_crystal_lattice_algorithm.md)
assume one lattice. A world streaming inward would have to change lattices
continuously, and what a jump of $`\xi_q`$ means when $`\xi`$ depends on where
the parent is standing is not defined anywhere in the project.

Following that through is a substantial piece of work touching the algorithm
specification and not only the analysis ladder, and it is **deliberately not
attempted here**. Logged as **Z-LS1**, with §6 recording as much of the
structure as is visible from this note.

**(ii) A uniform reach.** Take $`y_{\max} < \inf_x R(x) = \epsilon`$, attained
at the nucleus. Everything in the existing ladder then applies unchanged, at
the price of a reach fixed by the tightest point on the line — and, by Z2,
that point is the one place where nothing is emitted. The constraint is set
where the constrained quantity is irrelevant, which is an unsatisfying
situation to be in and is exactly the situation Theorem C6 describes for the
unsoftened pole taken to its limit.

§4 prices (ii).

---

## 4. Theorem Z4: the softening threshold

**Theorem Z4.** Take the uniform choice, $`y_{\max} < \epsilon`$, so
$`\Delta p > \pi\hbar/2\epsilon`$. For $`\epsilon \gg a_0`$ the well is
harmonic, and the ground state is resolved to $`k`$ rungs per $`\sigma_p`$ iff

```math
\epsilon \;\ge\; \frac{k^{4}\pi^{4}}{4}\; a_0 ,
\qquad a_0 = \frac{\hbar^{2}}{\mu Z} .
```

*Proof.* Expanding about the minimum,
$`V \simeq -Z/\epsilon + (Z/2\epsilon^{3})r^{2}`$, so
$`\tfrac12\mu\omega^{2} = Z/2\epsilon^{3}`$ and
$`\omega = \sqrt{Z/\mu\epsilon^{3}}`$. The ground state then has

```math
\sigma_r^{2} = \frac{\hbar}{2\mu\omega}
             = \frac{\hbar\,\epsilon^{3/2}}{2\sqrt{\mu Z}} ,
\qquad
\sigma_p = \frac{\hbar}{2\sigma_r} ,
```

which in the stated units is $`\sigma_r = \epsilon^{3/4}/\sqrt2`$ and
$`\sigma_p = 1/(\sqrt2\,\epsilon^{3/4})`$. Requiring
$`k\,\Delta p \le \sigma_p`$ at the tightest uniform reach gives

```math
\frac{k\pi}{2\epsilon} \;\le\; \frac{1}{\sqrt2\,\epsilon^{3/4}}
\qquad\Longleftrightarrow\qquad
\epsilon^{1/4} \;\ge\; \frac{k\pi}{\sqrt2} ,
```

and raising to the fourth power gives $`\epsilon \ge k^4\pi^4/4`$.
$`\square`$

So $`\epsilon_c = \pi^4/4 = 24.352\,a_0`$ for a single rung per $`\sigma_p`$,
$`4\pi^4 = 389.6`$ for two, $`81\pi^4/4 = 1972.5`$ for three. **The cost in
softening length is quartic in the resolution.**

Part D solves the eigenproblem exactly and bisects for the crossings:

```
        eps        E_0     sigma_r   eps^(3/4)/sqrt2   sigma_p    sigma_p/dp
       1.000  -0.669778      1.092            0.707   0.472041        0.301
       5.000  -0.163863      2.943            2.364   0.171629        0.546
      24.352  -0.037317      8.620            7.752   0.058174        0.902
      50.000  -0.018688     14.339           13.296   0.034924        1.112
     100.000  -0.009526     23.606           22.361   0.021199        1.350
     389.636  -0.002503     63.769           62.013   0.007842        1.945
    1000.000  -0.000984    127.971          125.743   0.003907        2.488

  measured crossing, k = 1:  eps =    34.60   ratio to predicted 1.421
  measured crossing, k = 2:  eps =   432.92   ratio to predicted 1.111
```

The harmonic approximation is the only one in the proof, and it enters through
$`\sigma_r`$ alone. The measured $`\sigma_r`$ exceeds the harmonic prediction
by 11 per cent at $`\epsilon = 24`$ and 1.8 per cent at $`\epsilon = 1000`$,
and the ratio of measured to predicted threshold falls from 1.42 to 1.11
correspondingly. So Z4 is a lower bound that tightens, and in the regime it is
stated for it is the right answer.

**What it says.** An atom cannot be represented on a uniform reach-limited
phase-space crystal. At $`\epsilon = a_0`$ the lattice resolves 0.3 of a rung
per $`\sigma_p`$ — not a marginal failure but a factor of three — and reaching
even one rung requires softening the core over two dozen Bohr radii, by which
point the binding energy has fallen from $`0.67`$ to $`0.037`$, a factor of
eighteen. This is the soft-core form of the obstruction **K-LS2** records for
the Eckart barrier, and it is sharper, for a reason worth stating: on a
barrier the packet width is free, so K1.1 constrains the *packet*
($`\sigma_r < a/2`$) and one can always choose a narrow one. A bound state has
no free width. The potential picks $`\sigma_r`$, so the reach condition
becomes a condition on the potential's own parameters and cannot be evaded by
choosing a better initial condition.

---

## 5. Theorem Z5: but the channel is not empty there

Z4 pushes towards large $`\epsilon`$, and large $`\epsilon`$ is a harmonic
well, and Theorem I4 says a harmonic well has no demographic channel at all.
If the two met, the usable window would be empty and this potential would be
no better than the Eckart barrier.

They do not meet.

**Theorem Z5.** At the tightest uniform reach $`y_{\max} = \epsilon`$, the
horizon spans

```math
\frac{y_{\max}}{\sigma_r} \;=\; \frac{\epsilon}{\epsilon^{3/4}/\sqrt2}
\;=\; \sqrt2\,\epsilon^{1/4} ,
```

which *grows* with $`\epsilon`$. So the reach always extends several state
widths beyond the harmonic core, into the region where $`V'''`$ is largest,
and $`\Gamma`$ over the state's own support does not vanish.

Part E:

```
        eps   eps/sigma_r    Gamma_max   Gamma(sigma_r)   ratio
      24.352        2.825   2.0266e-03     1.9733e-03    0.9737
     100.000        4.236   4.9351e-04     3.9116e-04    0.7926
     389.636        6.110   1.2666e-04     7.5699e-05    0.5977
    1000.000        7.814   4.9351e-05     2.3773e-05    0.4817
```

At the $`k = 1`$ threshold the state sits at 97 per cent of the peak emission
rate. At the $`k = 2`$ threshold, 60 per cent. A harmonic well would give
exactly zero — measured at $`1.2\times10^{-15}`$ in Theorem G2 — so these are
not residuals of a vanishing quantity.

$`\Gamma_{\max}`$ itself falls as roughly $`\epsilon^{-2}`$, which is the
statement that a shallow well emits less; that is a rate, and rates are
horizon parameters by Theorem G5 in any case. The ratio is the meaningful
column, and it is order unity where the lattice first works.

**So soft-core Coulomb is the right vehicle for S-SP6.** It has a stationary
state, a uniform lattice that can carry it, and a live demographic channel
over its support, all at the same $`\epsilon`$. The Eckart barrier has none of
those three simultaneously.

---

## 6. What Z-LS1 would have to settle

Deferred, but recorded here so the next attempt starts further along. If the
reach is taken locally, $`y_{\max}(x) = R(x)`$, then at least the following
must be given meanings they do not currently have.

- **The jump quantum.** $`\xi_q = q\,\Delta p(x)`$ with
  $`\Delta p(x) = \pi\hbar/2R(x)`$. An event at parent position $`x`$ deposits
  at $`p \pm q\Delta p(x)`$, so the daughter momenta of events at different
  positions do not lie on a common lattice, and the momentum axis is no longer
  discrete in any global sense. Either the lattice is abandoned in favour of a
  continuum in $`p`$ — which the crystal-lattice construction exists to avoid
  — or the daughters must be re-binned, which is a momentum change not
  produced by the classical force and so violates postulate (S).
- **Proposition K8.** The ionised pair must sit on the parent's own momentum
  row. That survives, since it is a statement about the parent's row and not
  about the spacing.
- **Theorem C4's identity of reach and momentum quantum.** It is a per-mode
  statement and survives pointwise. The budget ratio of Theorem K3 would then
  be a function of position rather than a single number.
- **The sea.** $`B = 2/h`$ is two pairs per Planck cell, and a Planck cell is
  $`\Delta x\,\Delta p`$. If $`\Delta p`$ varies with $`x`$ then either the
  sea density varies to compensate — keeping two pairs per $`h`$ — or the
  Planck cell is not the right unit. The first seems forced by the
  admissibility bound $`|W|\le2/h`$, which is a statement about $`h`$ and not
  about the lattice, but nothing here establishes it.
- **Theorem K1's own scope.** $`R(x)`$ is the radius of convergence of the
  Moyal series at $`x`$. A world at $`x`$ whose reach is $`R(x)`$ is exactly
  marginal, so the practical ceiling is some $`\theta R(x)`$ with
  $`\theta < 1`$, and whether $`\theta`$ can be taken constant across the line
  is unmeasured.

A cheap first probe, which does not require any of the above to be settled:
run the residual kernel at $`y_{\max}(x) = \theta\sqrt{x^2+\epsilon^2}`$ on a
fixed $`\Delta p`$ grid, re-binning the daughters, and measure how badly the
first moment of the kernel — which Theorem C3 requires to vanish — is
corrupted by the re-binning. If it is corrupted at the $`10^{-15}`$ level the
question is a bookkeeping one; if at $`10^{-3}`$ it is a physical one.

---

## 7. Numerical verification

`src/demo_soft_core_coulomb.py`.

| part | claim | result |
|---|---|---|
| A | Z1, closed form for $`V'''`$ | $`1.5`$ – $`2.1\times10^{-3}`$ against a fourth-order stencil |
| A | Z1, root at $`\epsilon\sqrt{3/2}`$ | six figures, five $`(Z,\epsilon)`$ pairs |
| B | Z2, dark nucleus | $`\Gamma(0) = 0`$ exactly, six reaches |
| B | Z2, ring survives | drifts $`+0.08\%`$ to $`+6.13\%`$ over a tenfold reach range |
| C | Z3, ceiling varies | $`\Delta p(x)/\Delta p(0)`$ from 1 to 0.0156 |
| D | Z4, threshold | crossings 34.60 and 432.92 against 24.35 and 389.6 |
| D | Z4, harmonic $`\sigma_r`$ | 11% high at $`\epsilon=24`$, 1.8% at $`\epsilon=1000`$ |
| E | Z5, channel alive | $`\Gamma(\sigma_r)/\Gamma_{\max}`$ = 0.974 at threshold |

Two caveats. The eigenproblem is a finite-difference tridiagonal solve on a
box of half-width $`30\epsilon^{3/4}`$ with 6001 points; $`E_0`$ at
$`\epsilon = 1`$ agrees with a Richardson extrapolation to five figures, but
the wide-$`\epsilon`$ rows are not independently checked. And the ring radius
in Part B is located by parabolic interpolation on a grid of spacing
$`2\times10^{-3}\epsilon`$, so the quoted shifts carry about $`0.1`$ per cent.

---

## 8. Open items

- **Z-LS1 (the non-uniform lattice).** §3 reading (i) and §6. Deferred by
  decision, not by difficulty: it touches
  [`../algorithm/phase_space_crystal_lattice_algorithm.md`](../algorithm/phase_space_crystal_lattice_algorithm.md)
  and the meaning of $`\xi_q`$, and should not be started as an aside. The
  cheap probe at the end of §6 is the recommended first step.
- **Z-LS2 (transmission).** The Eckart note's centre was Theorems K4–K6, the
  claim that the entire quantum correction to transmission is delivered by the
  residual channel as a small imbalance between two large opposed flows. The
  analogue here needs a reference transmission, and soft-core Coulomb has no
  closed form, so it needs a split-operator Schrödinger run instead. Planned
  as the next step.
- **Z-LS3 (the lobe ledger).** K7 splits the separatrix ledger by lobe and
  finds a fourfold cancellation dominated by the approach flank. The analogue
  here is a four-lobe ledger for a repulsive core, and it should look
  different, because the sech² lobes are arranged around a summit the packet
  crosses once and these are arranged around a centre it can orbit.
- **Z-LS4 (the ground state as a stationary ledger).** S-SP6's request. With
  Z4 and Z5 the vehicle exists: run `demo_stochastic_ledger.py`'s process, or
  the mesh ledger, on the $`\epsilon \approx 35`$ ground state and measure
  $`f`$, $`N`$ and $`S`$ in a genuinely stationary state rather than in a
  transient. This is also the cleanest available test of Theorem N5's standing
  population, since a stationary state removes the transport transient that
  N-SP4 worries about.
- **Z-LS5 (the quiet points are minima, not nodes).** §2, third reading. K7's
  nodes are exact zeros of $`V'''`$ but not of $`\Gamma`$ at finite reach, and
  K-LS7's prediction of "exact holes" in the sea profile should be restated
  with a depth attached. Cheap to fix and it changes what K-LS7 would be
  measuring.
- **Z-LS6 (repulsive, and scattering).** Everything in §1–§3 is
  sign-invariant, but the physics is not: a repulsive core is a barrier with a
  dark summit and four lobes, and is the direct Eckart analogue. Z-LS2 should
  probably be done in the repulsive case first, where the classical outcome
  functional of Theorem K4 has the same meaning it has on sech².

---

## 9. Sources

- Corollary K1.2 and Theorems K1, K3, K7, K8 of
  [`eckart_barrier_compensated.md`](eckart_barrier_compensated.md), and
  Theorems C4, C6, C7 of
  [`compensated_liouville_splitting.md`](compensated_liouville_splitting.md).
  The geometry this note re-runs.
- C. Eckart, *The penetration of a potential barrier by electrons*, Phys. Rev.
  **35** (1930) 1303–1309. For the contrast, and already cited by the note
  this one follows.
- J. Javanainen, J. H. Eberly and Q. Su, *Numerical simulations of
  multiphoton ionization and above-threshold electron spectra*, Phys. Rev. A
  **38** (1988) 3430–3446. The origin of the one-dimensional soft-core
  Coulomb potential $`-Z/\sqrt{x^2+\epsilon^2}`$ as a standard model atom, and
  the source of the convention used here.
- Q. Su and J. H. Eberly, *Model atom for multiphoton physics*, Phys. Rev. A
  **44** (1991) 5997–6008. The bound-state structure of the same potential,
  against which §4's eigenvalues can be checked.
- The softening value in common use in that strong-field literature is
  $`\epsilon^2 = 2`$, i.e. $`\epsilon = 1.414\,a_0`$, chosen so that the
  ground-state energy matches hydrogen's. That is a factor of seventeen below
  Z4's $`k = 1`$ threshold and three hundred below its $`k = 2`$ threshold,
  which is a useful measure of how far the reach-limited crystal sits from the
  regime that literature works in. Worth confirming against the primary
  sources above before quoting.
