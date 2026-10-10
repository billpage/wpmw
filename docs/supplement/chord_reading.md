# How the sea reads the potential

**A tutorial on postulate (S′), the chord, and the reading: what the sum over
chords estimates, how a reading becomes an event, and which parts are still
proposals.**

---

## 0. Status and provenance

Tutorial supplement. Companion demo: `src/demo_chord_reading.py`, which
produces every number in §§4 and 5 and the figure. Interactive companion:
[`../interactive/chord_reading.html`](../interactive/chord_reading.html)
([live page](https://billpage.github.io/wpmw/interactive/chord_reading.html)),
an animation of the sea under (S′) and (S) with a reader, its chords, its
reading and its integrators. The numbers of §7 come from
`src/demo_pair_collision_clocks.py` (Parts B2, B3 and C), run in October 2026.

Prompted by Bill Page's request for a tutorial on the chord mechanism, and by
his question: what exactly does the sum over chords estimate, in terms of the
Wigner function or the density matrix? §5 answers it.

Nothing in §§1–6 is new. The results are Propositions Q2, Q5 and Q9 to Q12
of [`../analysis/force_blind_sea.md`](../analysis/force_blind_sea.md) (step
23), Definition L0 and Theorems L1, L3, L4 and L5 of
[`../analysis/sea_phase_reference.md`](../analysis/sea_phase_reference.md)
(step 22), and Theorem Y4 of
[`../analysis/dark_sea_and_worldline_identity.md`](../analysis/dark_sea_and_worldline_identity.md)
(step 20). §5 restates Theorem L1 in density-matrix language and checks the
statement exactly. §7 records measurements and proposals made alongside
this tutorial; step 24's third addendum has since put the measured ones in
the ladder as Proposition V12, and the rest are open items there. **No existing document is edited by this supplement.**

---

## 1. The cast

| object | carries | moves under (S′) | contributes to $`W`$ |
|---|---|---|---|
| free positon | $`(x, p, \theta)`$, species $`\varepsilon = +1`$ | full classical force: $`\dot x = p/m`$, $`\dot p = -V'(x)`$ | $`+1`$ |
| free negaton | $`(x, p, \theta)`$, species $`\varepsilon = -1`$ | the same | $`-1`$ |
| aligned pair | a positon and a negaton at one point, with one momentum and one clock $`\theta`$ | inertially: $`\dot x = p/m`$, $`\dot p = 0`$ | $`0`$ |

Every clock winds at the Lagrangian rate, postulate P1 of
[`../analysis/phase_resonance_microdynamics.md`](../analysis/phase_resonance_microdynamics.md):

```math
\hbar\,\dot\theta \;=\; \frac{p^2}{2m} - V(x) .
```

The **sea** is a random configuration of aligned pairs with mean density
$`B = 1/\pi\hbar`$ per unit area of phase space. A **row** is a strip of
momentum of width $`dp`$. The **reach** $`y_{\max}`$ is how far apart, in
half-separation, two things may be and still be read together, and the
product is fixed (Q9(c)):

```math
B\,dp\,2y_{\max} \;=\; 1 ,
```

so the strip of one row and width $`2y_{\max}`$ holds one aligned pair on
average. With $`dp = 0.25`$ and $`\hbar = 1`$, $`y_{\max} = 2\pi`$.

**Postulate (S′)** (step 23 §1, as amended in step 24 §4.1): between events,
every free body obeys the full classical force; an aligned pair is
force-blind, moving inertially while its clock still winds with $`V`$, and it
changes momentum only in a collision with another aligned pair.

---

## 2. Postulate (S′) in two lines

Two consequences carry everything that follows (Proposition Q2).

**(a) A sea clock keeps its row's plane wave, up to the eikonal.** With
$`\varphi_p(x, t) = (px - p^2 t/2m)/\hbar`$ the phase of row $`p`$'s plane
wave,

```math
\frac{d}{dt}\Big[\hbar\theta - \hbar\varphi_p\big(x(t), t\big)\Big] \;=\; -V\big(x(t)\big),
```

so a pair that crosses the potential comes out shifted by the *eikonal*, the
phase a wave gains from a potential along its path, and otherwise stays on its
row's wave. Here it is the straight-line (Glauber) eikonal,
$`-\hbar^{-1}\int V\,dt`$, and it is exact, because under (S′) the pair's path
does not bend ([glossary](../GLOSSARY.md#eikonal)). In the interactive page this offset is the hue of each sea pair: a pair
changes hue while it crosses the barrier and keeps the new hue afterwards.

**(b) Two pairs on one row keep their separation, and their clocks drift
apart at exactly $`U`$.** Two pairs $`i, j`$ on row $`p`$ at $`x \mp y`$ have
the same velocity, so $`2y`$ never changes, and their kinetic winding is
identical, so

```math
\hbar\,\dot\mu \;=\; V(x_j) - V(x_i) \;=\; U(x, y) \qquad
\text{(an inertial reader)},
```

and a reader whose own momentum follows the force sees
$`U_{\rm res} = U - 2yV'_{\rm eff}(x)`$ instead. Under (S) the force would
fan each row out over its neighbours (Q2(c)), and this exact statement would
start to fail at $`V''`$.

Why (S′) and not (S): an aligned pair contributes nothing to $`W`$, so the
force has nothing to act on, and darkness is easiest to keep if the pair
does not bend (step 23 §2). The page has a switch between the two.

---

## 3. The chord

A **chord** is two different aligned pairs on the same row, within the reach
of each other:

```text
row p:   ●───────── 2y ─────────●        ● = one aligned pair (positon and negaton
        x−y          x          x+y          at one point, one shared clock)
                     ↑
                  the vertex, the chord's midpoint
```

Its data are the midpoint $`x`$, the half-separation $`y`$ with
$`\lvert y\rvert \le y_{\max}`$, and the **misalignment**, the two clocks
compared at a common point with a lever momentum $`P`$:

```math
\mu^{(P)}_{ij} \;=\; \theta_i - \theta_j + \frac{P\,d}{\hbar},
\qquad d = x_j - x_i = 2y .
```

A chord is a relation between two pairs, not an object. Under (S′) its $`y`$
is fixed and its $`\mu`$ winds as in §2(b). Read against the daughter row
$`p + \xi_q`$ instead of its own, it carries the channel's Fourier phase
(Q5(a)):

```math
\mu^{(p)}_{ij} + \frac{\xi_q\,d}{\hbar} \;=\; \mu^{(p + \xi_q)}_{ij},
\qquad \xi_q = q\,dp .
```

**What the two ends are, and are not.** In the Wigner transform the points
$`x \pm y`$ are the two arguments of one body's density matrix,
$`\rho(x+y, x-y)`$: the ket and the bra, two worlds' answers to where the one
body is ([`what_the_reach_is.md`](what_the_reach_is.md) §1). A chord's ends
are not that body's ket and bra. They are two reference clocks that happen
to sit at those two places. A body has no off-diagonal clock of its own to
read, so the sea supplies one: by Theorem L4 the sea enters only as unsigned
clocks, a phase reference that contributes nothing to $`W`$, as a receiver's
local oscillator contributes nothing to its output.

---

## 4. The reading

A **reader** at the vertex $`(x, p)`$ adds up every chord of its own row whose
midpoint lies near it. Proposition Q5(b): with

```math
z_q \;=\; \sum_{\text{chords}} w(y)\,e^{\,i\mu^{(p + \xi_q)}_{ij}},
```

summed over the chords of row $`p`$ with midpoints in a cell of length
$`\Delta x`$ at $`x`$, and with the clocks freshly re-locked so that
$`\mu^{(p)} = 0`$,

```math
K_q(x) \;=\; \frac{1}{B\,dp\,\Delta x}\;\frac{d}{dt}\,\mathrm{Re}\,z_q .
```

Written out, this is Theorem L1:

```math
K_q(x) \;=\; -\frac{B\,dp}{\hbar}\int_{-y_{\max}}^{y_{\max}}
U_{\rm res}(x, y)\,w(y)\,\sin\!\Big(\frac{2\xi_q y}{\hbar}\Big)\,dy .
```

Each factor has a world-particle meaning (step 22 §1): $`B\,dp`$ is the
number of pairs per unit length of the row; $`U_{\rm res}/\hbar`$ is the rate
at which the potential winds a chord at separation $`2y`$, with its uniform
part (the classical force) removed; $`\sin(2\xi_q y/\hbar)`$ is the
misalignment at distance $`y`$ of a pair split by $`2\xi_q`$ and aligned at
$`x`$; $`w`$ is the reach's taper. Each chord contributes the rate at which
its own term of the daughter row's interference pattern changes, and the
sum over chords is the integral over $`y`$, sampled by the separations the sea
happens to have.

As code, one reader's reading (the function `reading_from` of
`demo_sea_resonance_clock.py`, simplified):

```python
def read(reader, sea, dX):
    """Theorem L1 as a sum over chords (Q5(b), Q12): the reader's estimate of K_q."""
    P = reader.row * dp                        # the lever: the reader's own row
    n = nu * B * dp                            # expected pairs per unit length of a row
    ends = [a for a in sea if a.row == reader.row and abs(a.x - reader.x) < y_max + dX / 2]
    Khat = zeros(NQ)
    for i, j in pairs_of(sorted(ends, key=x)):            # every chord once, x_i < x_j
        d = j.x - i.x                                     # 2y
        if d > 2 * y_max or abs((i.x + j.x) / 2 - reader.x) > dX / 2:
            continue
        mu = i.theta - j.theta + P * d / hbar             # the chord's misalignment
        rate = V(j.x) - V(i.x) - d * Veff_prime(reader.x) # hbar dmu/dt = U_res
        w = cos(pi * d / (4 * y_max)) ** 2                # the reach's taper
        for q in range(1, NQ + 1):
            Khat[q] -= B * dp / (hbar * n**2 * dX) * w * rate * sin(xi[q] * d / hbar + mu)
    Khat -= xi * (xi @ Khat) / (xi @ xi)                  # Proposition L2(c): no net force
    return Khat
```

**How good is one reading?** Part B of the demo draws a static locked sea
($`\mu = 0`$) of multiplicity $`\nu`$ around a reader on the Eckart flank,
$`x = -0.625`$, with $`\Delta x = 1`$ and $`K_1 = 1.0766`$:

| $`\nu`$ | readings | chords per reading | no chord at all | mean of $`\hat K_1/K_1`$ | mean of $`\hat K_2/K_2`$ | spread of one reading, over $`\lvert K_1\rvert`$ |
|---|---|---|---|---|---|---|
| 1 | 400 000 | 0.08 | 0.933 | 1.007 ± 0.008 | 1.06 | 5.3 |
| 2 | 200 000 | 0.32 | 0.781 | 1.005 ± 0.007 | 1.01 | 3.0 |
| 4 | 80 000 | 1.28 | 0.441 | 1.001 ± 0.006 | 1.02 | 1.7 |
| 8 | 30 000 | 5.1 | 0.087 | 0.995 ± 0.006 | 0.97 | 1.07 |
| 16 | 12 000 | 20 | 0.002 | 1.001 ± 0.006 | 1.07 | 0.70 |
| 32 | 4 000 | 81 | 0 | 0.997 ± 0.008 | 1.08 | 0.48 |
| 64 | 1 500 | 324 | 0 | 0.987 ± 0.009 | 1.02 | 0.33 |

The reading is unbiased at every density: it is a U-statistic of the sea's
configuration (Hoeffding 1948), and its mean is the integral. One reading is
not. At the physical density, $`\nu = 1`$, a reader has a chord in its
aperture about one time in twelve, and its spread falls only between
$`1/\sqrt\nu`$ and $`1/\nu`$. Proposition Q5 called this critical sampling,
and Q12(f) measured it in the particle model.

![The chords estimate the kernel; one reading is critically sampled](https://raw.githubusercontent.com/billpage/wpmw/output/figures/chord_reading_estimator.png)

*Left: the kernel $`K_q`$ at $`x = -0.625`$ (blue) and the mean chord reading
at $`\nu = 16`$ (orange, ±2 standard errors). Right: the spread of one
reading of $`K_1`$ against $`\nu`$, with the mean number of chords per
reading.*

---

## 5. What the sum over chords estimates

**It estimates the generator, not the state.** Precisely: the jump-rate
kernel by which the potential moves $`W`$ between momentum rows, which is the
Wigner representation of the commutator with $`V`$, with its classical part
removed. It contains no $`W`$ and no $`\rho`$.

The derivation is three lines. The potential's part of the von Neumann
equation, in position representation, is a phase rate on each matrix
element:

```math
\partial_t\,\rho(X, X')\Big|_V \;=\; \frac{V(X) - V(X')}{i\hbar}\;\rho(X, X') .
```

Put $`X = x + y`$, $`X' = x - y`$, so that the rate is $`U(x, y)/i\hbar`$,
and Wigner-transform both sides with

```math
W(x, p) \;=\; \frac{1}{\pi\hbar}\int \rho(x+y,\,x-y)\,e^{-2ipy/\hbar}\,dy,
\qquad
\rho(x+y,\,x-y) \;=\; \int W(x, p')\,e^{\,2ip'y/\hbar}\,dp' .
```

Since $`U`$ is odd in $`y`$, the result is a convolution in momentum:

```math
\partial_t W(x, p)\Big|_V \;=\; \int J(x,\,p - p')\,W(x, p')\,dp',
\qquad
J(x, \xi) \;=\; -\frac{B}{\hbar}\int U(x, y)\,\sin\!\Big(\frac{2\xi y}{\hbar}\Big)\,dy,
\qquad B = \frac{1}{\pi\hbar} .
```

Part A1 of the demo checks this on a cat state in the Eckart barrier, on a
$`256 \times 256`$ grid where the discrete transform is exact: the
convolution and the transformed commutator agree to $`1.7\times10^{-14}`$.
The jump rates of the algorithm are $`J`$ sampled at the channel transfers,
with the taper $`w`$ applied and the classical part moved to streaming
(Theorem L1):

```math
K_q(x) \;=\; J_{\rm res}(x,\,\xi_q)\,dp,
\qquad
\partial_t W(x, p)\Big|_{\rm res} \;=\; \sum_{q \ge 1} K_q(x)\,\big[W(x,\,p - \xi_q) - W(x,\,p + \xi_q)\big] .
```

Part A2 compares the quadrature with the code's kernel: $`K_1 = 1.076607`$
against 1.076080 at $`x = -0.625`$, and $`V'_{\rm eff} = 0.768418`$ against
0.768067, the difference being the code's 64-point sampling of $`y`$.

So each ingredient of the chord sum has a counterpart in the density matrix
and in $`W`$:

| density matrix | Wigner | the sea |
|---|---|---|
| the off-diagonal coordinate: $`\rho(x+y, x-y)`$ | the variable conjugate to $`p`$ | a chord's two ends, at $`x \pm y`$ |
| the phase rate $`(V(x+y) - V(x-y))/\hbar`$ that the potential applies to that element | $`U`$ inside $`J`$ | the chord's winding, $`\hbar\dot\mu = U`$ (§2(b)) |
| the Fourier factor $`e^{-2ipy/\hbar}`$ | the lever in $`\sin(2\xi y/\hbar)`$ | the daughter row's plane wave (Q5(a)) |
| the measure $`dy/\pi\hbar`$ | $`B = 1/\pi\hbar`$ | the sea's density: chords per unit $`y`$ |
| the state, $`\rho`$ | $`W`$ | the bodies, which are the parents, not the sea |

The sea reads the operator; the bodies carry the state; an event is their
product. An event fires at $`\lvert K_q(x)\rvert`$ per parent, and the
parents are distributed as $`W`$, so the expected traffic is
$`\sum_q K_q\,[W(p - \xi_q) - W(p + \xi_q)]`$, the last equation above.

The analogy is a diffraction grating in a spectrometer (Q5 used Fraunhofer
diffraction, which is the same thing). The grating and the phase plate in
front of it decide how strongly each order $`q`$ is diffracted; the light
decides how much there is to diffract. The sea and the potential are the
grating and the plate. The bodies are the light.

Two consequences follow.

- **The reading is blind to the state** (Part A3). A packet, a cat or an
  empty box give the same $`K_q`$ at the same $`x`$, and each row's chords
  estimate the same $`K_q(x)`$: the row only fixes the lever. A cat's
  negativity costs the sea nothing until the dynamics makes some
  ([`../analysis/sea_depletion.md`](../analysis/sea_depletion.md),
  Proposition V1).
- **Only the rate of $`z_q`$ carries $`K_q`$, not $`z_q`$ itself.** At
  re-lock, $`\mathrm{Re}\,z_q = \sum w\cos(2\xi_q y/\hbar)`$ is the transform
  of the sea's own uniform spread of separations, which vanishes on average
  for $`q \ge 1`$ (Q5(c), the orthogonality that fixes the reach). As time
  passes, $`\mu`$ accumulates as $`U\tau/\hbar`$ and the weight drifts from
  the sine: a reading needs a fresh reference, which is what dark catalysis
  supplies (Q6).

---

## 6. From a reading to an event

**What cannot trigger an event.** Not a crossing. An absorption is two free
bodies crossing at the parent on rows $`p \pm \xi_q`$, and that is how it is
*realised* (Q9(b)); but a rate triggered by crossings is mass action in
$`u_\pm`$, and $`E`$ would no longer close (Q10). Not a single chord either:
fired one chord at a time, at the size of each chord's own contribution to
$`K_q`$, the same drive runs 8 to 16 times faster than $`\sum_q\lvert K_q\rvert`$
on the Eckart flank, because the kernel's falloff in $`q`$ is an
interference between chords (`demo_pair_collision_clocks.py`, Part B1).

**What does.** The reader's own clock (Q11): one integrator per channel,
advanced by the reading and firing whenever it crosses an integer. Over any
interval the signed count is $`\int K_q\,dt`$ to within one; the counts are
sub-Poissonian; the integrators must be banked, one per channel. With
$`C_q`$ fed by the reading, a firing is one fringe of the daughter row's
interference pattern passing the reader.

The whole step, as `particle_run` of `demo_sea_resonance_clock.py` does it:

```python
for each step dt:
    # 1. stream, postulate (S')
    for b in free_bodies:
        b.p -= V_prime(b.x) * dt;  b.x += b.p / m * dt
        b.theta += (b.p**2 / (2 * m) - V(b.x)) * dt / hbar
    for a in sea_pairs:                                   # force-blind
        a.x += a.p / m * dt
        a.theta += (a.p**2 / (2 * m) - V(a.x)) * dt / hbar
    # 2. read and integrate (Q11): one integrator per channel per body
    for b in free_bodies:
        K = read(b, sea, dX)            # the live reading, or the kernel K_q(b.x)
        for q in channels:
            old = floor(b.C[q]);  b.C[q] += K[q] * dt
            for each integer crossed: fire(b, q, sign of the crossing)
    # 3. dark catalysis (Theorem L5, Q6): re-lock sea clocks at the kernel's rate

def fire(b, q, s):
    tau = s * b.eps                     # the event puts +tau at p + xi, -tau at p - xi
    neg = free negaton near (b.x, b.p + tau * xi[q])
    pos = free positon near (b.x, b.p - tau * xi[q])
    if neg and pos:                     # absorption: they kink to b.p and bind
        make an aligned pair of neg and pos at the vertex
    elif a pair sits near (b.x, b.p):   # emission: ionise it
        positon to b.p + tau * xi[q];  negaton to b.p - tau * xi[q]
    else:
        the event fails (counted)
```

"Near" means within the row width $`dp`$ in momentum and a cell in position.
Both realisations move $`W`$ identically; they differ in the kinetic energy
$`\xi_q^2/m`$ that passes between the bodies and the sea (Q9(a)).

**What is measured** (Q12). Integrators fed by the kernel realise the QLE's
target with no excess firing and less noise than Poisson triggering. Fed by
the live reading, the open loop converges, $`1 - \text{correlation}`$ falling
as $`1/\nu`$; the closed loop does not: at $`\nu = 64`$ the trigger fires 1.7
to 1.9 times too often and realises 0.37 to 0.67 of the target. That is open
item Q-SP11.

---

## 7. Sea readers: proposals, not yet in the ladder

Everything in this section was measured in October 2026. Step 24's third
addendum has since written (a) and the mean field of (b) into the ladder as
Proposition V12
([`sea_depletion.md` §4.5](../analysis/sea_depletion.md#45-proposition-v12-the-clock-at-a-collision));
the sea reader itself is open item V-SP13 there. It is here because it follows from §§4 to 6 and
answers the question of what moves the sea's pairs.

**(a) The clock rule at a collision must re-lock.** Step 24 §4.1 amends
(S′): aligned pairs change momentum only in Focus and Defocus collisions
(Definition (C)), and V-SP8 asks what their clocks do. Three rules were run
in the particle model of §6 (`demo_pair_collision_clocks.py`, Part C;
$`\nu = 32`$, four seeds, to $`t = 14`$, Kaggle, from commit `7bf4df1` with
that script added): *whole pairs* (each pair kinks as a unit, clocks
continuous), *exchange* (partners swap) and *phase-matched exchange*, against
*re-locking* (whole pairs that take their new row's clocks, step 22 §6). With
dark catalysis at twice the kernel's rate:

| collisions | grayness $`\langle\sin^2(\mu/2)\rangle`$ | lock coherence | live reading, slope / correlation |
|---|---|---|---|
| none | 0.004 | 0.74 | 0.71 / 0.76 |
| whole pairs | 0.005 | 0.22 | 0.14 / 0.37 |
| exchange | 0.119 | 0.19 | 0.13 / 0.35 |
| phase-matched exchange | 0.028 | 0.34 | 0.22 / 0.50 |
| re-locked | 0.003 | 0.85 | 0.69 / 0.80 |

A collision either makes gray pairs (exchange) or knocks pairs off their new
row's lock (whole pairs), and the reading collapses either way. Only
re-locking keeps it; without dark catalysis, re-locking alone holds the
reading at 0.53 against 0.13 with no collisions. For dark parents the
exchange products carry $`\mu + \varphi`$ and $`-\mu + \varphi`$ with
$`\varphi = o\,\xi d/\hbar`$, so flipping the orientation $`o`$ leaves their
darkness unchanged, and no orientation bias was measured.

**(b) An aligned pair can be a reader too.** The same chords, read by an
aligned pair at the vertex instead of a free body, give the same reading and
the same rate, Theorem L5's. Its firing is $`W`$-null, and it is realised by
what is in the aperture:

```python
def fire_sea(A, q, t):                  # A: an aligned pair at row p (a proposal)
    C = an aligned pair on A's own row, within A's aperture
    D = an aligned pair on row p + 2 t xi[q], within A's aperture
    pick C or D in proportion to their numbers
    if C:   A to p + t xi[q];  C to p - t xi[q]          # defocus
    elif D: A to p + t xi[q];  D to p + t xi[q]          # focus: both land on p + t xi
    else:   the dark-catalysis round trip (clocks only)
    re-lock the moved pairs to their new rows' clocks     # (a)
```

The firing pair is itself one of the two pairs moved, so the event is
two-body: defocus goes as $`s_n s_n`$ and focus as $`s_{n-q}s_{n+q}`$, the
form of Definition (C). In mean field (Part B3) the uniform row is an exact
equilibrium, no ripple grows, a dip of 30 per cent refills to within 4 per
cent by $`t = 40`$, and pairs are conserved to $`10^{-13}`$, with rates either
proportional to the partners present or saturating as $`1 - e^{-m}`$. A
version in which a catalyst stays put and moves two *other* pairs is
three-body; it cannot be in detailed balance, and with absorption-first
allocation it clumps the sea in momentum (Part B2).

**(c) Why a defocus cannot be triggered by geometry.** Two pairs on one row
are at rest relative to each other (§2(b)): nothing geometric marks the
moment they should separate, and the only thing about them that changes in
time is $`\mu`$, at $`U/\hbar`$. If crossings alone triggered focus, focus
would run in free space where defocus cannot, and the sea would condense in
momentum. So, as for bodies, realisation is geometric (a crossing realises a
focus, a chord a defocus) and triggering is spectral, for both.

---

## 8. What this does not settle

- **The closed loop** (Q-SP11): the live reading does not yet drive events
  correctly when acted on.
- **Critical sampling**: at the physical density a reading is a sample, not
  a measurement (§4).
- **Integrator banks** (Q-SP14): Q11 gives them to free bodies,
  provisionally; §7(b) would give them to aligned pairs too.
- **The reader's frame.** An inertial pair reads $`U`$, not
  $`U_{\rm res}`$. Removing the reading's first moment (Proposition L2(c))
  discards the classical part in any frame; that this makes a sea reader's
  rate the kernel's is expected, not yet measured.
- **Gray pairs** (Q-SP10): pairs made by contact or by absorption bring
  unrelated clocks, and only a phase-setting rule aligns them.

---

## 9. Sources

- Adler, R. *A study of locking phenomena in oscillators*, Proc. IRE **34**
  (1946) 351–357.
- Goodman, J. W. *Introduction to Fourier Optics*, 3rd ed., Roberts & Co.
  (2005), chapter 4.
- Hoeffding, W. *A class of statistics with asymptotically normal
  distribution*, Ann. Math. Statist. **19** (1948) 293–325.
- Inose, H., Yasuda, Y. and Murakami, J. *A telemetering system by code
  modulation — Δ-Σ modulation*, IRE Trans. Space Electron. Telem. **SET-8**
  (1962) 204–209.
- Moyal, J. E. *Quantum mechanics as a statistical theory*, Proc. Cambridge
  Philos. Soc. **45** (1949) 99–124.
- Wigner, E. *On the quantum correction for thermodynamic equilibrium*,
  Phys. Rev. **40** (1932) 749–759.
