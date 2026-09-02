# Compensated Liouville Algorithm

**A complete specification for evolving the Wigner distribution as exact Newtonian flow plus a zero-mean signed hop channel.**

---

## 0. Status of this specification

This document promotes
[`../analysis/compensated_liouville_splitting.md`](../analysis/compensated_liouville_splitting.md)
(ladder step 14, theorems C0–C7) into an implementable specification.
Companion verification: `src/demo_compensated_liouville_algorithm.py`.

**The one-sentence relation to the rest of `algorithm/`.** This is
[`phase_space_crystal_lattice_algorithm.md`](phase_space_crystal_lattice_algorithm.md)
with the rate field replaced and one step added: the mediated-jump rule of its
[§3b](phase_space_crystal_lattice_algorithm.md#3b-potential-driven-mediated-jumps) is used verbatim, the rate field $`\Gamma_q(x)`$ is replaced by a
compensated field $`\Gamma^{\mathrm{res}}_q(x)`$ computed from $`V`$ without
any Fourier decomposition of $`V`$, and a deterministic acceleration step is
placed in front of it.

**What the analysis note supplies** and this document does not re-derive: the
symbol, the exactness of the split (C1), its identification with the cubic
Taylor remainder (C2), the moment conditions at bounded reach (C3), the
crossover condition (C4), and the quiet region (C7).

**What is new here.** Four things, all of them consequences of putting the
split on a grid rather than on the real line, and all verified in the
companion demo:

1. §2 — the reach and the momentum cell are not merely "one parameter" in
   the loose sense of C4; the grid is *forced* by the reach, and the number of
   momentum cells is the ket–bra rung count. Neither is a convergence knob.
2. §3.2 — the deterministic step must be compensated against the **kernel's
   own first moment**, not against $`V'(x)`$. Using $`V'(x)`$ leaves up to
   $`4\times10^{-2}`$ of spurious force in the hop channel.
3. §4.1 — the Nyquist rung must be zeroed, or a real $`W`$ acquires an
   imaginary part at the $`10^{-2}`$ level every substep.
4. §4.4 — **the coherence horizon needs a profile, not a cutoff.** Under the
   hard cutoff of Definition (H) the total hop rate grows like $`\ln N_p`$ and
   the momentum churn grows like $`N_p`$, so neither has a limit and the
   discrete operator has no third moment — no semiclassical expansion at all.
   A soft profile fixes all three and is specified here as the default, with
   the hard cutoff as the special case $`w \equiv 1`$.
5. §4.1, §5 — the profile multiplies **the residual only**. Windowing the
   whole symbol is the obvious move and it silently destroys the Newtonian
   step.

Item 4 is also an erratum for the note's §4: the total-variation figures
tabulated there are quoted as absolute quantities and are in fact functions of
the grid used to compute them. The ratios are stable to about ten per cent and
the note's qualitative conclusion is unaffected.

Item 4 has a consequence for the ladder that this document cannot settle on its
own. Definition (H) of
[`../analysis/open_position_space.md`](../analysis/open_position_space.md) §3
states the horizon as a sharp fact — a pair beyond $`L_c`$ "is not
instantiated". What §4.4 shows is that the sharp version is the one choice of
profile under which the jump channel has an infinite event rate. If the
specified default here is right, Definition (H) wants restating as a **rung
occupancy** $`w(y)`$ rather than a cutoff, of which the indicator function is
one member and not a good one. That is a change to a postulate, and it is
flagged here rather than made.

Items marked **[choice]** are implementation decisions not fixed by the
underlying physics. Items marked **[open]** need further work. Items marked
**[normative]** are requirements: an implementation that violates one is
wrong, not merely different.

---

## 1. The split

### 1.1 The symbol

Let $`s`$ be conjugate to momentum, $`y = \hbar s/2`$ the half ket–bra
separation. Take the partial Fourier transform of $`W`$ in $`p`$ alone,
$`\hat W(x, s) = \int W(x,p)\thinspace e^{-ips}\thinspace dp`$. The entire
potential term of the QLE is then multiplication by

```math
M(x, s) \;=\; \frac{i}{\hbar}\bigl[V(x + y) - V(x - y)\bigr] ,
\qquad y = \frac{\hbar s}{2} .
```

Subtracting the part linear in $`y`$ gives the split this document
implements,

```math
M \;=\; \underbrace{i\thinspace V'(x)\thinspace s}_{M_{\mathrm{cl}}}
\;+\;
\underbrace{\frac{i}{\hbar}\bigl[V(x+y) - V(x-y) - 2y\thinspace V'(x)\bigr]}_{M_{\mathrm{res}}} .
```

### 1.2 Why the two steps can be applied in sequence

| | statement | where |
|---|---|---|
| C1 | $`e^{\tau M} = e^{\tau M_{\mathrm{cl}}}e^{\tau M_{\mathrm{res}}}`$ exactly, for every $`\tau`$ — both factors are multiplication operators in $`(x,s)`$ | note §2, §2.1 |
| C2 | $`M_{\mathrm{res}}`$ is the odd part of the cubic Taylor remainder of $`V`$ about $`x`$; it vanishes iff $`V''' \equiv 0`$ on the reach | note §2.2 |
| C3 | at bounded reach the residual kernel has zero zeroth **and** zero first moment: it conserves worlds and delivers no net momentum | note §3 |
| C7 | $`V''' \equiv 0`$ on $`[x-y_{\max}, x+y_{\max}]`$ implies no events at $`x`$ | note §5 |

C1 is what licenses the sequential application with **no Trotter error inside
the potential substep**. It says nothing about the Strang error between
streaming and the potential, which this algorithm still has; §5 is where that
is accounted for.

C3 is what makes the two steps non-overlapping: the whole classical force sits
in step 1 and step 2 does not apply it a second time.

---

## 2. Grids and parameters: the reach *is* the momentum grid

This section is the part of the specification most likely to be got wrong,
because it inverts the usual relation between a grid and an answer.

### 2.1 The primitives

Three numbers are chosen, and everything else follows.

| symbol | meaning | kind |
|---|---|---|
| $`y_{\max}`$ | the **reach**: half the coherence horizon $`L_c`$; the greatest distance from its own position at which a world consults $`V`$ | physical postulate |
| $`N_p`$ | the number of ket–bra **rungs** inside the horizon; even | physical postulate |
| $`M_x`$ | position cells | numerical |

$`y_{\max}`$ and $`N_p`$ are inherited from
[`../analysis/open_position_space.md`](../analysis/open_position_space.md) §3,
Definitions (H) and (R). They are statements about which ket–bra pairs exist.
Only $`M_x`$ is a resolution parameter in the ordinary sense.

### 2.2 What follows

```math
\Delta p \;=\; \frac{\pi\hbar}{2\thinspace y_{\max}} \;=\; \frac{\pi\hbar}{L_c},
\qquad
P_{\max} \;=\; \frac{N_p\thinspace\Delta p}{2},
\qquad
a \;=\; \frac{2 L_c}{N_p} ,
```

where $`a`$ is the rung spacing in the separation coordinate. The momentum
grid and the rung grid are exact Fourier duals of one another:

```
     rungs in separation                momentum transfers
     |y| <= y_max, spacing a/2          xi = q dp, |q| < N_p/2

  -y_max            +y_max                -P_max        +P_max
     |...............|          <--->        |.............|
     <-> a/2                                 <-> dp = pi hbar / L_c
     N_p rungs                               N_p cells

     resolution in y   <---->   extent in xi
     extent in y       <---->   resolution in xi
```

Read across: **the reach sets the momentum cell; the rung count sets the
momentum extent.** Verified in the demo (Part A): at $`y_{\max} = 1`$ and
$`N_p = 64`$, $`\Delta p = 1.5708`$, $`P_{\max} = 50.265`$, $`a = 0.0625`$,
and exactly $`64`$ rungs fall inside the horizon.

### 2.3 Refining the momentum grid is not a refinement **[normative]**

Halving $`\Delta p`$ doubles $`y_{\max}`$. It does not converge anything; it
asserts that worlds consult the potential twice as far away, which is a
different physical model with a different answer. In the limit
$`\Delta p \to 0`$ one recovers the untruncated Wigner equation, where the
residual is a differential operator rather than a jump measure and Theorem C3
fails outright.

The consequence for testing: **a grid-convergence study in $`\Delta p`$ is
meaningless for this algorithm.** The convergence parameters are $`\Delta t`$
and $`M_x`$, and those only.

The analogy that keeps this straight: the reach is an **aperture** in the
ket–bra separation, and the momentum-transfer spectrum is its diffraction
pattern. Narrowing the aperture coarsens the pattern — that is optics, not
error. §4.4 pushes the analogy one step further, and it turns out to predict
the event rate.

### 2.4 Choosing the primitives

- $`y_{\max}`$ from the physics if it is known; otherwise from C4, which wants
  $`u = k\thinspace y_{\max} \ll \pi`$ for the dominant scale $`k`$ on which
  $`V`$ varies. $`u \le \pi/4`$ is a reasonable working value: the residual
  is then about a tenth of the classical term.
- $`N_p`$ large enough that $`P_{\max}`$ covers the state's momentum support,
  and that $`\sigma_p`$ spans several cells. With a short reach $`\Delta p`$
  is large, so this is the binding constraint. Note the reach condition of C4
  reads $`\sigma_p \gg \hbar k/2`$ — a state too narrow in momentum to be
  resolved on the grid is also a state for which the split buys nothing.
- $`M_x`$ to resolve $`V'''`$, and to satisfy the CFL condition of §5.2.
- Position domain: the open line. See §2.5.

### 2.5 The open line is the setting; the ring is a special case

**Nothing in this algorithm needs a position box.** The only transform is in
the rung direction, and the rung ladder is the coherence structure of a single
world, not a boundary condition. Position enters as a spectator parameter
throughout (note §2.1): the rate field at $`x`$ is built by sampling $`V`$ on
$`[x - y_{\max},\; x + y_{\max}]`$ and depends on nothing else.

So the specification is written for $`x \in \mathbb{R}`$, represented by a
window wide enough that

- $`W`$ is negligible at the edge, with margin $`P_{\max}\Delta t/m`$ so that
  no departure point of §3.3 leaves the window, and
- the window exceeds the active set of §6.3 by at least one reach.

Widening the window then adds only quiet cells: it costs storage and buys
nothing, which is the practical form of Theorem C7.

The ring enters in exactly three places, all of them special cases.

1. **Where a periodic $`V`$ quantises the jumps.** Theorem O5 gives a momentum
   quantum $`\pi\hbar/a`$ from the period of $`V`$, on all of $`\mathbb{R}`$,
   with no box. This is independent of the reach-induced $`\Delta p`$ of §2.2
   and the two coincide only when $`y_{\max} = a/2`$.
2. **Where the algorithm degenerates to the crystal lattice.** §4.2: at
   $`y_{\max} = L/2`$ on a ring of circumference $`L`$ with a single-mode
   $`V`$, the uncompensated rate field is exactly the two-atom stencil
   $`\Gamma_q`$.
3. **As a source of closed-form counterexamples**, per note §6. It is not a
   home: note §6.3 shows a ring pins every world at maximal reach, which is
   precisely the value of the one parameter the reach condition depends on at
   which the split buys nothing. **A ring is not a valid benchmark for this
   algorithm [normative].** §5's measurement uses a periodic window only
   because the test potential is a cosine, so that spectral advection is exact
   and the flow integrator does not contaminate the result.

---

## 3. Step 1 — deterministic acceleration

### 3.1 It must be a flow, not another split **[normative]**

Step 1 is the classical Liouville evolution

```math
\partial_t W \;=\; -\frac{p}{m}\thinspace\partial_x W \;+\; V'(x)\thinspace\partial_p W ,
```

and it must be integrated **as a flow** — by characteristics, or by a
symplectic integrator on worlds. Strang-splitting streaming against the
classical force inside step 1 would reinstate exactly the commutator term the
reordering was meant to remove, and the algorithm would gain nothing but an
extra substep.

This is the single structural requirement of the whole specification and it is
easy to violate by accident, because the classical force *is* available as a
multiplication operator in $`(x,s)`$ and applying it there is one line.

### 3.2 The force is the kernel's first moment **[normative]**

Theorem C3 says the residual delivers no net momentum. That is a statement
about the continuum. On the rung grid the first moment of the sampled kernel
is not $`V'(x)`$, so subtracting $`i V'(x) s`$ leaves a residue.

Instead, define the deterministic acceleration by the discrete first moment of
the full kernel — Theorem O2 of
[`../analysis/open_position_space.md`](../analysis/open_position_space.md)
read as a definition rather than as a result:

```math
V'_{\mathrm{eff}}(x) \;=\;
  \sum_{q} \xi_q\thinspace K(x, \xi_q) \;\Big/\; \sum_{q} \xi_q\thinspace B(\xi_q) ,
\qquad
a(x) \;=\; -\thinspace V'_{\mathrm{eff}}(x) ,
```

where $`K`$ is the kernel of $`M`$ on the rung grid, $`B`$ the kernel of the
reference symbol $`i s`$, both sums run over the whole grid, and $`a(x)`$ is
the acceleration the flow of §3.1 integrates. Then set
$`M_{\mathrm{cl}} = i\thinspace V'_{\mathrm{eff}}(x)\thinspace s`$ and
$`M_{\mathrm{res}} = M - M_{\mathrm{cl}}`$ as before.

With this definition the residual's discrete first moment vanishes to machine
precision at **every** $`N_p`$, and $`a(x) \to V'(x)`$ as $`N_p`$ grows.
Verified (Part C) at $`x = 1`$:

| potential | $`y_{\max}`$ | $`N_p`$ | residue, naive $`V'`$ | residue, discrete |
|---|---|---|---|---|
| cosine well | 1.00 | 32 | $`6.4\times10^{-3}`$ | $`2.3\times10^{-15}`$ |
| cosine well | 2.00 | 32 | $`2.3\times10^{-2}`$ | $`2.7\times10^{-15}`$ |
| Gaussian barrier | 2.00 | 32 | $`3.9\times10^{-2}`$ | $`3.1\times10^{-15}`$ |
| soft-core Coulomb | 2.00 | 32 | $`4.4\times10^{-2}`$ | $`2.2\times10^{-15}`$ |

against forces of order $`0.5`$ to $`0.8`$, so the naive compensation misses
by up to six per cent at a long reach and a small rung count. A spurious force
in a channel specified to carry none is not an accuracy question — it is a
violation of the conservation law the split exists to establish.

### 3.3 The two forms

- **Mesh.** Semi-Lagrangian transport: for each cell, integrate the
  Hamiltonian flow backwards over $`\Delta t`$ and interpolate $`W`$ there.
  **[choice]** of interpolant; cubic or better, since linear interpolation is
  diffusive enough to swamp the residual channel. Departure points must land
  inside the window (§2.5). Spectral advection is available only on a ring and
  is not the general form.
- **World ensemble.** A symplectic integrator per world — velocity Verlet is
  sufficient and is the recommended default **[choice]**. This is the natural
  open-line form: no interpolation, no boundary, and it is where the split
  pays, since it costs one force evaluation per world per step against a
  global transform for the hop channel.

---

## 4. Step 2 — the hop channel

### 4.1 Building the rate field

For each position cell $`x_m`$:

1. Evaluate $`V`$ at the $`2N_p`$ points $`x_m \pm y_k`$, where
   $`y_k = \hbar s_k/2`$ are the rungs, and form $`M(x_m, s_k)`$.
2. **Zero the Nyquist rung [normative].** With $`N_p`$ even the grid holds
   $`y = -y_{\max}`$ but not $`+y_{\max}`$, so the odd symbol is sampled
   asymmetrically and one rung has no partner. The symmetric treatment of a
   maximal-separation rung visited from both sides is the average of the two
   endpoint values, and for an odd function that average is zero. Set
   $`M(x_m, s_{N_p/2}) = 0`$, and the same for the reference symbol $`is`$ used
   in §3.2.
3. Form $`M_{\mathrm{res}} = M - M_{\mathrm{cl}}`$ per §3.2.
4. **Apply the horizon profile to the residual, and only to the residual
   [normative]:** $`M_{\mathrm{res}} \leftarrow w(y/y_{\max})\thinspace
   M_{\mathrm{res}}`$, with $`w`$ as specified in §4.4. $`M_{\mathrm{cl}}`$ is
   left alone.
5. Inverse-transform in $`s`$ to get the signed transfer rates
   $`K_{\mathrm{res}}(x_m, \xi_q)`$.
6. Set $`\Gamma^{\mathrm{res}}_q(x_m) = -K_{\mathrm{res}}(x_m, \xi_q)`$.

Step 2 is not cosmetic. Verified (Part B), cosine well, $`y_{\max} = 1`$,
$`N_p = 64`$:

| | $`\max\lvert\mathrm{Im}\thinspace K\rvert/\max\lvert K\rvert`$ | oddness in $`q`$ | $`\max\lvert\mathrm{Im}\thinspace W\rvert`$ after one substep |
|---|---|---|---|
| as sampled | $`1.15\times10^{-1}`$ | $`2.30\times10^{-1}`$ | $`1.63\times10^{-2}`$ |
| Nyquist zeroed | $`1.09\times10^{-16}`$ | $`2.34\times10^{-16}`$ | $`5.6\times10^{-16}`$ |

Left alone, a single rung drives a real $`W`$ complex at the percent level
every substep, and the realness and oddness of the kernel — the properties
that give number conservation for free and force the sign structure — are lost
at the same order.

A soft horizon does not make step 2 redundant. The profile zeroes the residual
at that rung for free, since $`w(\pm 1) = 0`$, but $`V'_{\mathrm{eff}}`$ is
read off the *unwindowed* $`M`$ and stays contaminated.

**Why the profile goes on the residual only.** $`w(s)\thinspace i V'(x) s`$ is
not the symbol of a drift, so windowing the whole symbol takes the Newtonian
step apart. The test that shows it is validation rung 7: for an open-line
quadratic $`V`$ the residual is identically zero, so the algorithm must reduce
to exact Newtonian flow. Verified (Part H), $`V = m\omega^2x^2/2`$ with
$`\omega = 1.7`$, at $`x = 1`$:

| placement | $`\max\lvert M_{\mathrm{res}}\rvert`$ | error in the force |
|---|---|---|
| hard horizon | $`4.6\times10^{-13}`$ | $`2.8\times10^{-14}`$ |
| window the whole symbol | $`2.2\times10^{1}`$ | $`5.5\times10^{-2}`$ |
| window the residual | $`1.2\times10^{-13}`$ | $`2.8\times10^{-14}`$ |

Windowing the whole symbol manufactures a residual of size 22 in a potential
that has none, and shifts the classical force by two per cent.

Read physically the rule is not a trick: the horizon grades the **coherence**
channel. The classical force is the local, first-moment part of the kernel and
is not a coherence effect, so there is nothing there for a rung occupancy to
attenuate.

**No Fourier decomposition of $`V`$ appears anywhere.** The rate field is
built by sampling $`V`$ over each world's own reach. This is the note's §1
point as an implementation fact, and §6.2 there is the reason it matters: the
mode-by-mode bound on the residual diverges even where the residual is exactly
zero.

### 4.2 The update rule — unchanged

With $`\Gamma^{\mathrm{res}}_q`$ in hand, the rule is §3b/§3c of
[`phase_space_crystal_lattice_algorithm.md`](phase_space_crystal_lattice_algorithm.md)
verbatim:

```math
W(x_m, p_n) \;\mathrel{+}=\; \Delta t \sum_{q \ge 1}
  \Gamma^{\mathrm{res}}_q(x_m)\thinspace
  \bigl[\thinspace W(x_m, p_{n+q}) - W(x_m, p_{n-q})\thinspace\bigr] .
```

A world at cell $`(m,n)`$ mediates a transfer between cells $`n+q`$ and
$`n-q`$; the mediator is unchanged; no world is created or destroyed. Verified
(Part D) against the spectral action of $`M_{\mathrm{res}}`$ to
$`7.1\times10^{-16}`$ relative.

The correspondence is exact in the other direction too. Run the *uncompensated*
field at maximal reach $`y_{\max} = L/2`$ on a single-mode potential and the
kernel collapses onto the two atoms $`q = \pm 1`$ with
$`-K(x, \xi_1) = \Gamma_1(x)`$ to $`8.9\times10^{-16}`$, everything else below
$`2.2\times10^{-16}`$. So the crystal-lattice algorithm is the
$`\Gamma^{\mathrm{res}} \to \Gamma`$, $`y_{\max} \to L/2`$ corner of this one.

### 4.3 Mesh and world forms

- **Mesh, exact.** Apply $`\exp(\Delta t\thinspace M_{\mathrm{res}})`$ in the
  $`(x, s)`$ representation: FFT each momentum column, multiply, inverse FFT.
  Exact in $`\Delta t`$ for the substep, and this is the recommended default
  for the mesh form. Explicit Euler on the stencil of §4.2 amplifies
  negativity and must not be used for production runs **[normative]** — this
  is the standing project result on the momentum substep, unchanged here.
- **World ensemble, sampled.** Total rate at $`x`$ is
  $`R(x) = \sum_q \lvert K_{\mathrm{res}}(x, \xi_q)\rvert`$. Draw an event
  time from $`R(x)`$, draw $`q`$ with probability
  $`\lvert K_q\rvert / R`$, and apply the transfer.

The sign of $`K_q`$ has to go somewhere, and there is no way to make it go
away. $`K_{\mathrm{res}}`$ is odd and non-zero, hence takes both signs, hence
is never the generator of a one-body Markov jump process — Proposition T3 of
[`../supplement/takabayasi_1954_stochastic_picture.md`](../supplement/takabayasi_1954_stochastic_picture.md),
unweakened by compensation. The two admissible unravelings are the project's
standing pair:

- **positon/negaton [choice, recommended]** — the sign is a species label; a
  negative-rate transfer produces a world of opposite species, and
  annihilation is required. This is the ontologically serious reading and it
  matches [`../analysis/species_sectors_and_annihilation.md`](../analysis/species_sectors_and_annihilation.md).
- **signed weights [choice]** — weights multiply by $`\mathrm{sign}(K_q)`$ and
  every weight grows as $`e^{R t}`$ to compensate the depletion. Cheaper per
  step, with the sign problem relocated into variance.

Compensation shrinks $`R`$; it does nothing to the sign structure. That
distinction is worth keeping sharp, because "the quantum channel now carries
no force" invites the misreading that it has become classical.

### 4.4 The horizon profile, and the event budget

**The specified horizon is soft [normative].** The occupancy profile
$`w(y/y_{\max})`$ must be

| requirement | what it buys |
|---|---|
| even | both moment conditions of C3 survive **exactly** |
| $`w(0) = 1`$ | the leading $`\hbar^2`$ Moyal term is untouched |
| $`w(\pm 1) = 0`$ | no seam, so the event rate converges |
| $`w'(\pm 1) = 0`$ | $`1/q^3`$ tail, so the momentum churn converges |
| supported on the reach | the quiet region of C7 survives |

The raised cosine $`w(t) = \cos^2(\pi t/2)`$ meets all five and is the
default **[choice]**. The hard cutoff of Definition (H) is the special case
$`w \equiv 1`$, and it is the one choice that fails the last three.

**Why the hard cutoff fails.** It is a hard-edged aperture in the separation
coordinate. $`M_{\mathrm{res}}`$ is odd, so its periodic extension over the
rung window carries a jump of $`2\lvert M_{\mathrm{res}}(x,s_{\max})\rvert`$ at
the seam, and a jump gives coefficients falling as $`1/q`$:

```math
\lvert K_{\mathrm{res}}(x, \xi_q)\rvert \;\simeq\;
  \frac{\lvert M_{\mathrm{res}}(x, s_{\max})\rvert}{\pi\thinspace q}
\qquad\Longrightarrow\qquad
R(x) \;\simeq\; \frac{2}{\pi}\thinspace
  \bigl\lvert M_{\mathrm{res}}(x, s_{\max})\bigr\rvert \thinspace \ln N_p \;+\; O(1) .
```

Verified (Part E), fitted slope against predicted, at $`x = 1`$:

| potential | $`y_{\max}`$ | $`R`$ at $`N_p = 256`$ | $`R`$ at $`N_p = 16384`$ | fitted slope | $`(2/\pi)\lvert M_{\mathrm{res}}(s_{\max})\rvert`$ |
|---|---|---|---|---|---|
| cosine well | 0.25 | 0.00714 | 0.01429 | 0.00170 | 0.00170 |
| cosine well | 1.00 | 0.44541 | 0.89001 | 0.10598 | 0.10573 |
| cosine well | 2.00 | 3.27872 | 6.52013 | 0.77268 | 0.77085 |
| Gaussian barrier | 1.00 | 0.23750 | 0.45822 | 0.05292 | 0.05285 |
| Gaussian barrier | 2.00 | 5.66761 | 11.10204 | 1.29542 | 1.29236 |

Four figures of agreement. But the diverging event count is the mildest of the
three symptoms. Cosine well, $`y_{\max} = 1`$, measured tail exponent in
brackets:

| horizon | quantity | $`N_p = 256`$ | 1024 | 4096 | 16384 |
|---|---|---|---|---|---|
| hard $`(q^{-1.03})`$ | $`R = \sum_q\lvert K\rvert`$ | 0.44541 | 0.59539 | 0.74308 | 0.89001 |
| | $`\sum_q\lvert \xi K\rvert`$ | 14.35 | 58.56 | 235.4 | 942.7 |
| | $`\sum_q \xi^3 K`$ | 31.91 | 128.98 | 517.2 | 2070.3 |
| raised cosine $`(q^{-3.00})`$ | $`R`$ | 0.01631 | 0.01657 | 0.01663 | 0.01665 |
| | $`\sum_q\lvert \xi K\rvert`$ | 0.05339 | 0.05447 | 0.05475 | 0.05482 |
| | $`\sum_q \xi^3 K`$ | $`-0.13381`$ | $`-0.12981`$ | $`-0.12880`$ | $`-0.12847`$ |

The event rate diverges logarithmically, the **momentum churn diverges
linearly**, and the third moment diverges linearly too — so under a hard
horizon the discrete operator has no third moment, which is to say no
semiclassical expansion, term by term, at all.

Under the soft horizon all three converge, and the third moment converges to
$`-0.128465`$ — which is the analytic leading Moyal coefficient,
$`-0.128465`$, to six figures. **A soft horizon is what gives the discrete
operator a Moyal expansion.** That is the argument for making it the default,
and it is a stronger one than the numerical convenience: a hard cutoff makes
the model's own semiclassical limit ill-defined.

Two things the profile does not cost.

- **Conservation.** $`w`$ even keeps
  $`\sum_q K_{\mathrm{res}} = \sum_q \xi_q K_{\mathrm{res}} = 0`$ exactly, for
  any profile. The compensation of §3.2 is untouched.
- **The leading quantum correction.** With $`M_{\mathrm{res}}(0) =
  M_{\mathrm{res}}'(0) = M_{\mathrm{res}}''(0) = 0`$ and $`w`$ even,
  $`(w M_{\mathrm{res}})'''(0) = w(0)\thinspace M_{\mathrm{res}}'''(0)`$, so
  $`w(0)=1`$ preserves the $`\hbar^2`$ term exactly. The profile attenuates only
  the far rungs — which are exactly the ones producing the sidelobes.

**Quote $`N_p`$ with any event rate anyway.** Even under a soft horizon $`R`$
depends on $`N_p`$ until it has converged, and $`N_p`$ is a physical parameter
(the rung count), not a tolerance. This is the erratum of §0.

**The reach still controls the rate.** The seam value obeys the C2 bound
$`\lvert M_{\mathrm{res}}\rvert \le (2/\hbar)(y_{\max}^3/6)\sup\lvert V'''\rvert`$,
so short reach means cheap channel under either horizon. C4 survives intact.

![Event budget under a hard and a soft horizon](https://raw.githubusercontent.com/billpage/wpmw/output/figures/compensated_algorithm_budget.png)

## 5. Composition and step size

### 5.1 What the reordering buys

Strang, at macro step $`\Delta t`$:

```math
\Phi_{\mathrm{cl}}(\Delta t/2) \;\circ\;
\exp\bigl(\Delta t\thinspace M_{\mathrm{res}}\bigr) \;\circ\;
\Phi_{\mathrm{cl}}(\Delta t/2) ,
```

with $`\Phi_{\mathrm{cl}}`$ the classical flow of §3. Both this and the
uncompensated scheme split the **same** operator
$`M_{\mathrm{cl}} + M_{\mathrm{res}}`$ and both are second order; the
reordering shrinks the **constant** by the symbol ratio, because the surviving
commutator is with $`M_{\mathrm{res}}`$ rather than with the whole potential
term. Verified (Part F), cosine well, $`y_{\max} = 1`$ so $`u = \pi/4`$,
$`L^1`$ error at $`t = 1`$, hard horizon:

| $`\Delta t`$ | uncompensated | compensated | gain |
|---|---|---|---|
| 0.1250 | $`4.23\times10^{-1}`$ | $`4.40\times10^{-2}`$ | 9.6 |
| 0.0312 | $`3.50\times10^{-2}`$ | $`3.69\times10^{-3}`$ | 9.5 |
| 0.0156 | $`2.42\times10^{-3}`$ | $`2.51\times10^{-4}`$ | 9.7 |
| 0.0078 | $`2.28\times10^{-4}`$ | $`2.19\times10^{-5}`$ | 10.4 |

against a predicted
$`\lvert M_{\mathrm{res}}\rvert/\lvert M_{\mathrm{cl}}\rvert = 0.0904`$, i.e.
$`11.1\times`$, and the leading estimate $`u^2/6 = 0.103`$. The gain is flat
in $`\Delta t`$, which is the signature of a changed constant rather than a
changed order.

Under the specified soft horizon the residual is smaller again — the profile
attenuates the far rungs, where $`M_{\mathrm{res}}`$ is largest, giving a ratio
of $`0.0072`$ — and the measured gains run into the hundreds. Those figures are
lower bounds: the compensated error there sits close to the reference's own
floor, so the comparison stops being a measurement. The hard-horizon column is
the honest one to quote for the *mechanism*; the soft horizon simply has less
residual to split badly.

Be clear about what this is not. The scheme is still second order; the
Strang error between streaming and the hop channel is still there; and on a
mesh the classical flow costs about what the hop channel costs, so the gain is
in accuracy, not in time. In the world ensemble the gain is in both.

### 5.2 Choosing $`\Delta t`$

Three conditions, whichever binds:

- CFL on streaming: $`\Delta t \le m\thinspace\Delta x / P_{\max}`$.
- Event rate: $`\max_x R(x)\thinspace\Delta t \ll 1`$, with $`R`$ from §4.4.
  In the world form this is the mean events per world per step; in the mesh
  form the exact substep of §4.3 removes the condition, and only the Strang
  error remains.
- Splitting accuracy: from §5.1, and this is usually the loosest of the three
  once the reach is short.

---

## 6. Reference pseudocode

### 6.1 Setup (once, for static $`V`$)

```python
import numpy as np

dp   = np.pi * hbar / (2 * y_max)          # forced by the reach
p    = (np.arange(N_p) - N_p // 2) * dp
s    = 2 * np.pi * np.fft.fftfreq(N_p, d=dp)
y    = hbar * s / 2                        # the rungs
q    = np.round(np.fft.fftfreq(N_p, d=1.0) * N_p).astype(int)
xi   = q * dp                              # momentum transfer atoms
nyq  = N_p // 2

X, Y = np.meshgrid(x, y, indexing="ij")     # x is an OPEN-line window
M    = (1j / hbar) * (V(X + Y) - V(X - Y)) # sample V over each reach
S    = np.broadcast_to(s, M.shape).copy()
M[:, nyq] = 0.0                            # 4.1 [normative]
S[:, nyq] = 0.0

def first_moment(sym):                     # discrete first moment of a kernel
    return np.real((xi * np.fft.ifft(sym, axis=-1)).sum(axis=-1))

dV_eff = first_moment(M) / first_moment(1j * S)      # 3.2 [normative]
accel  = -dV_eff                           # what the flow of 3.1 integrates
M_cl   = 1j * dV_eff[:, None] * S
M_res  = M - M_cl

w      = np.cos(np.pi * y / (2 * y_max)) ** 2        # 4.4 horizon profile
M_res *= w[None, :]                        # the RESIDUAL only [normative]

K_res  = np.real(np.fft.ifft(M_res, axis=1))         # signed transfer rates
Gamma  = -K_res                                      # 4.2 orientation
rate   = np.abs(K_res).sum(axis=1)                   # R(x), see 4.4
active = rate > tol                                  # 6.3
```

### 6.2 The step

```python
def step(W, dt):
    W = classical_flow(W, dt / 2, accel)   # 3.1: a FLOW, not a split
    Wh = np.fft.fft(W[active], axis=1)     # 6.3: quiet cells skipped
    Wh *= np.exp(dt * M_res[active])       # 4.3, exact in dt
    W[active] = np.real(np.fft.ifft(Wh, axis=1))
    return classical_flow(W, dt / 2, accel)
```

The world-ensemble step replaces `classical_flow` by a per-world symplectic
integrator and the middle block by sampled transfers at rate `rate[m]`,
choosing `q` with weight `|K_res[m, q]|` and applying the sign per §4.3.

### 6.3 The active set

Theorem C7 is an optimisation as well as a theorem: where $`V'''`$ vanishes
over the whole reach the rate field is identically zero, and the hop channel —
and, in the world form, all the bookkeeping around it — can be skipped
outright.

Verified (Part G) on a harmonic trap plus a $`C^\infty`$ bump supported on
$`[-1,1]`$, over $`x \in [-6, 6]`$:

| reach $`y_{\max}`$ | 0.25 | 0.50 | 1.00 | 2.00 |
|---|---|---|---|---|
| predicted active region | $`\lvert x\rvert < 1.25`$ | $`< 1.50`$ | $`< 2.00`$ | $`< 3.00`$ |
| measured active fraction | 0.200 | 0.242 | 0.321 | 0.483 |

The measured edge tracks $`b + y_{\max}`$ with unit slope, falling silent two
to four position cells inside it because the bump's third derivative
approaches zero smoothly. **Mask on the measured rate field, not on an
analytic edge [normative]** — a potential without compact support has no exact
edge, and §5.1 of the note shows the interaction there is the potential's own
profile translated outward by one reach, not a cliff.

The crystal-lattice algorithm has no counterpart to this saving: its rate
field $`\Gamma_q(x) = -(V_q/\hbar)\sin(2\pi q x/L + \phi_q)`$ is non-zero at
almost every $`x`$ for every mode, and the cancellation that makes the
residual vanish is spread across the whole mode sum (note §6.2).

![Splitting error and the active set](https://raw.githubusercontent.com/billpage/wpmw/output/figures/compensated_algorithm_validation.png)

### 6.4 Cost

| | setup | per step |
|---|---|---|
| mesh, this algorithm | $`O(M_x N_p \log N_p)`$ | $`O(M_x^{\mathrm{act}} N_p \log N_p)`$ + flow |
| mesh, crystal lattice | $`O(Q M_x)`$ | $`O(Q M_x N_p)`$ |
| worlds, this algorithm | as above | $`O(N_w)`$ + events at rate $`R`$ |

$`M_x^{\mathrm{act}}`$ is the active-set size. The setup cost is paid once for
a static $`V`$ and repaid every step; for a time-dependent $`V`$ it is paid
every step and the comparison changes **[open]**.

---

## 7. Validation sequence

Rungs an implementation should pass, in order. Each is checked in the
companion demo at the section given.

1. **Split identity.** $`M - (M_{\mathrm{cl}} + M_{\mathrm{res}})`$ is
   identically zero on the grid. (Trivial, but catches Nyquist-handling
   mismatches between the two symbols.)
2. **Kernel structure.** $`K_{\mathrm{res}}`$ real and odd to $`10^{-16}`$;
   zeroth moment zero. §4.1.
3. **Momentum neutrality.** First moment of $`K_{\mathrm{res}}`$ zero to
   $`10^{-15}`$ at every $`N_p`$. §3.2.
4. **Stencil equals spectral.** The §4.2 rule reproduces
   $`\exp(\Delta t M_{\mathrm{res}})`$ to first order in $`\Delta t`$ and the
   generators agree to $`10^{-15}`$. §4.2.
5. **Crystal-lattice corner.** Uncompensated, single mode, maximal reach:
   two atoms, $`-K_1 = \Gamma_1`$. §4.2.
6. **Quiet region.** Harmonic trap plus compact bump: rate field below
   tolerance beyond $`b + y_{\max}`$, edge tracking the reach. §6.3.
7. **Harmonic exactness.** Open-line quadratic $`V`$: $`M_{\mathrm{res}}`$
   identically zero, so the algorithm reduces to the exact classical flow and
   the error is the flow integrator's alone. This is the sharpest single test;
   it is the one the crystal-lattice algorithm cannot pass, and it is also the
   test that catches a horizon profile applied in the wrong place (§4.1).
7b. **Open-line operation.** No wrap anywhere; the active set compact for a
   trap-plus-bump potential; widening the window changes nothing. §2.5, and
   Part I of the demo.
7c. **Sampler neutrality.** In the world form, the signed mean transfer over
   many drawn events is at the statistical floor — the sampler inherits C3.
   Measured: $`4.3\times10^{-3}`$ of the mean transfer scale over 3635 events,
   against a floor of $`1.7\times10^{-2}`$. Part I.
8. **Convergence against a reference.** Second order in $`\Delta t`$, with the
   constant below the uncompensated scheme's by the symbol ratio. §5.1.
9. **Agreement with `wigner_split_fourier`** on the QHO, and with the
   crystal-lattice solver on a cosine well at maximal reach.

---

## 8. The `wpmwlib` interface

Sketch for the module this specification is meant to become, following the
conventions of `wpmwlib/phase_space_crystal_lattice.py` — array shape
$`(M_x, N_p)`$, position axis first **[choice]**; note this transposes that
module's convention and the two must not be mixed silently.

```python
class CompensatedLiouville:
    """Wigner evolution as Newtonian flow plus a zero-mean signed hop channel.

    Parameters
    ----------
    y_max : float          reach; sets dp = pi*hbar/(2*y_max).  NOT a tolerance.
    N_p   : int            rung count (even); sets P_max = N_p*dp/2.
    M_x   : int            position cells.
    box   : float          position window.
    V, dV : callable       potential and its derivative (dV for diagnostics only).
    mass, hbar : float
    flow  : {'semi_lagrangian', 'symplectic'}
    horizon : {'raised_cosine', 'hard'}    see 4.4; soft is the default, and
                                          'hard' (w = 1) is the special case.
    """

    # --- setup -------------------------------------------------------
    def build_rate_field(self) -> None: ...     # 4.1; sets .Gamma, .rate, .active
    @property
    def accel(self) -> np.ndarray: ...          # 3.2, the kernel's first moment
    @property
    def dp(self) -> float: ...                  # forced; read-only

    # --- substeps ----------------------------------------------------
    def step_classical(self, dt) -> None: ...   # 3.1
    def step_residual(self, dt) -> None: ...    # 4.3, exact
    def step(self, dt) -> None: ...             # 5.1, Strang

    # --- diagnostics -------------------------------------------------
    def moments(self) -> tuple: ...             # zeroth and first, per x
    def event_budget(self) -> np.ndarray: ...   # R(x); meaningless without N_p
    def negativity(self) -> float: ...
```

Four interface requirements follow from the body of this document.
`dp` is derived and read-only, so that no caller can refine it (§2.3).
`accel` is derived from the kernel, not from `dV` (§3.2). `build_rate_field`
zeroes the Nyquist rung before anything else touches the symbol, and applies
the horizon profile to the residual only, after the force has been extracted
(§4.1). And there is no `L` or `periodic` parameter: `box` is a window on the
open line, not a circumference (§2.5).

---

## 9. Open implementation questions

- **CLA1 (what sets the profile).** §4.4 makes the horizon's *shape* a free
  function of the model, not just its width, and shows that the sharp
  Definition (H) is the one shape under which the event rate, the momentum
  churn and the third moment all fail to exist. Two questions follow. What
  physically sets $`w`$ — is it a decoherence profile, a rung-occupancy
  statistic of the sea, or something else? And is there a distinguished
  profile, e.g. the one minimising $`\sum_q\lvert\xi_q K_q\rvert`$ at fixed
  $`w(0) = 1`$? This is CLS4 of the note ("what sets $`L_c`$?") with more
  structure to bite on, and it is now the load-bearing open item: the
  specification's default is a **[choice]** where a derivation is wanted.
- **CLA1b (Definition (H)).** If CLA1 resolves in favour of a profile, the
  coherence horizon should be restated in
  [`../analysis/open_position_space.md`](../analysis/open_position_space.md)
  §3 as a rung occupancy rather than a cutoff. Theorems O1–O5 need re-reading
  against that: O5's coset invariant is unaffected (it depends on the transfer
  lattice, not on occupancy), but Proposition O4's localisation edge becomes a
  profile edge, and C7 becomes exact only where $`w`$ has compact support —
  which the raised cosine does.
- **CLA2 (world-form validation).** Everything in §5.1 is measured on a mesh.
  The claimed cost advantage is a world-ensemble claim and is unmeasured. The
  decisive test is a signed-ensemble run against the mesh reference at fixed
  variance.
- **CLA3 (annihilation load).** §4.3 leaves the positon/negaton unraveling
  specified but unpriced. Since compensation reduces $`R`$ by the symbol
  ratio, it should reduce pair production proportionally — but
  [`../analysis/species_sectors_and_annihilation.md`](../analysis/species_sectors_and_annihilation.md)
  prices annihilation against nodal structure, not against event rate, and
  the two need reconciling.
- **CLA4 (time-dependent $`V`$).** The whole rate field is rebuilt each step,
  and §6.4's cost comparison inverts. Whether the active set still pays is
  unknown.
- **CLA5 (the sub-$`\Delta p`$ acceleration).** This is CLS3 of the note,
  partly answered. Step 1 moves momentum by an amount that is generically not
  a multiple of $`\Delta p`$, so it cannot be represented on the momentum
  lattice at all; §3.3's two forms both evade this by keeping momentum
  continuous within step 1 (semi-Lagrangian, or per-world). Whether a strictly
  lattice-resident version exists — and hence whether Theorem O5's coset
  invariant survives the compensated split as an exact statement — is open.
- **CLA6 (multi-dimensional).** The symbol construction is per-degree-of-
  freedom and the rung grid becomes a lattice in $`\mathbf{y}`$, so the setup
  cost scales as $`(M_x N_p)^d`$. The active set may be the only thing that
  makes $`d = 2`$ affordable; untested.

---

## 10. Sources

- [`../analysis/compensated_liouville_splitting.md`](../analysis/compensated_liouville_splitting.md) — theorems C0–C7, and the physical argument this specification implements.
- [`phase_space_crystal_lattice_algorithm.md`](phase_space_crystal_lattice_algorithm.md) — §3b/§3c, the jump rule reused verbatim in §4.2.
- [`../analysis/open_position_space.md`](../analysis/open_position_space.md) — Definitions (H) and (R), Theorem O2 (the first moment is the force), Theorem O5 (the coset invariant).
- [`../supplement/takabayasi_1954_stochastic_picture.md`](../supplement/takabayasi_1954_stochastic_picture.md) — Proposition T3, on why the signed kernel admits no one-body Markov unraveling.
- [`../analysis/species_sectors_and_annihilation.md`](../analysis/species_sectors_and_annihilation.md) — the positon/negaton unraveling of §4.3.
- David Cyganski, *Extended Fokker–Planck Eq. and the QLE V2* (project memo).
