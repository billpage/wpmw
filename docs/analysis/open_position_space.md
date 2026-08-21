# Open position space: the horizon, the coset, and what the box was hiding

**Status.** Analysis note, step 10 of the ladder. Companion demo:
`src/demo_open_position_space.py`. Prompted by open item 4 of
[`../supplement/inverted_pair_barrier.md`](../supplement/inverted_pair_barrier.md),
which asked whether an absorbing layer can be given a microdynamic reading.

---

## 0. What this note inherits, retracts and corrects

**Inherits.** From
[`../algorithm/phase_space_crystal_lattice_algorithm.md`](../algorithm/phase_space_crystal_lattice_algorithm.md)
§1 and §3b: the momentum quantum $`\Delta p = \pi\hbar/L`$, the mode $`q`$
stencil and its rate field $`\Gamma_q(x)`$. From
[`four_rule_microdynamics_equivalence.md`](four_rule_microdynamics_equivalence.md):
the four-action channels and the exactness family. From
[`position_pair_ladder.md`](position_pair_ladder.md): the rung index
$`k = m - n`$ — the discrete ket–bra separation, which turns out to be the
variable that the whole of this note is about. From
[`fourd_microdynamics.md`](fourd_microdynamics.md): Theorem A3, the
orthogonality criterion for a conserved direction.

**Retracts** the framing of open item 4 of
[`../supplement/inverted_pair_barrier.md`](../supplement/inverted_pair_barrier.md).
That item asked for "an absorbing layer with a microdynamic reading — worlds
leaving the box rather than a numerical sponge". There is no layer, because
there is no boundary. Both substeps of the microdynamics are **pointwise in
$`x`$**: streaming uses only the world's own $`p`$, and the vertex rate uses
only $`V`$ evaluated at the world's own $`x`$. Nothing consults a neighbour,
so nothing needs a boundary condition. The escape problem of that note's §4.2
is an artefact of the *mesh*, not of the ontology. §7 below shows that the
one candidate object with the right name — a complex absorbing potential — is
strictly worse than the alternative, for three separate reasons.

**Corrects two statements of the algorithm specifications.**

1. `../algorithm/phase_alignment_microdynamics_algorithm.md` §1.2 says
   $`dp = \pi\hbar/L`$ "is fixed by the ring circumference". True on a ring,
   but the circumference is one of **three** independent mechanisms that can
   fix it (§1), and the other two work on open position space. The sentence as
   written makes the box look load-bearing when it is not.
2. `../algorithm/phase_space_crystal_lattice_algorithm.md` §1 offers
   $`\Delta p = \pi\hbar/(K L)`$ as a route to "finer momentum resolution".
   It is not: the $`K`$ sub-lattices never exchange anything (Theorem O5,
   leak $`0.0`$ exactly for $`K = 1, 2, 3, 4`$), so the refinement interleaves
   $`K`$ independent copies of the same problem rather than resolving any one
   of them better. Within a sector the resolution is fixed by $`V`$.

**Contradicts, gently, the hope that prompted §4.** Periodicity of $`V`$ is a
genuinely better assumption than compactness of $`x`$, and §4 makes that case
in full. But it does not rescue the project's current test problem: the
inverted harmonic pair potential $`-\tfrac{\mu}{2}\Omega^2 r^2`$ has no
period, so the coset mechanism does not apply to it at all. For that problem
the coherence horizon of §3 is the only route, and §4.4 gives the honest
accounting of what the periodic route does and does not buy elsewhere.

---

## 1. Tutorial: the box was doing two unrelated jobs

Every calculation in this repository so far has been on a periodic ring of
circumference $`L`$, and $`L`$ has silently been doing two jobs:

- **Job 1.** Fixing the momentum quantum, $`\Delta p = \pi\hbar/L`$.
- **Job 2.** Keeping worlds in view.

Job 2 is not a job at all. It is a property of the mesh, not of the worlds. A
world-particle streams by $`x \to x + (p/m)\thinspace\delta t`$ and jumps by
$`p \to p \pm q\thinspace\Delta p`$ at a rate $`|\Gamma_q(x)|`$ that depends
on its own position and on nothing else. Put the worlds on the whole line and
they simply spread out; there is no wall to hit and nothing to absorb.

Job 1 is real, and it can be discharged in three different ways:

| Mechanism | What is discrete | Set by | Exact? |
|---|---|---|---|
| **(R)** ring, $`y`$ wraps at $`L`$ | the momentum *argument* of $`W`$ | box circumference $`L`$ | yes, but only on a box |
| **(H)** coherence horizon | the momentum *argument* of $`W`$ | max ket–bra separation $`L_c`$ | approximate (§3) |
| **(P)** periodic potential | the momentum *jumps* | potential period $`a`$ | exact, on all of $`\mathbb{R}`$ |

(R) and (H) are the *same* formula wearing different clothes. The Wigner
transform is bilinear in the state through $`\rho(x + y/2,\thinspace x - y/2)`$,
and what discretises momentum is boundedness of the ket–bra separation $`y`$.
On a ring, $`\rho`$ has period $`2L`$ in $`y`$ (shifting $`y`$ by $`2L`$
returns both arguments to themselves), giving $`\Delta p = \pi\hbar/L`$. Bound
the separation by fiat at $`L_c`$ instead and the same computation gives
$`\Delta p = \pi\hbar/L_c`$. **The ring circumference and the coherence length
enter through the identical slot: the largest ket–bra separation the model
admits.**

Mechanism (P) is a different animal altogether and is the subject of §4.

**An analogy.** Confusing job 1 with job 2 is confusing the crystal lattice of
a solid with the walls of the sample holder. The Bragg peaks come from the
lattice; enlarging the holder does not move them. The name *phase-space
crystal lattice* turns out to be more literal than it looked: for a periodic
potential the momentum lattice **is** the potential's reciprocal lattice,
halved, and the box never had anything to do with it.

---

## 2. The Wigner kernel has no position envelope

Write the collision term in kernel form,

```math
\partial_t W + \frac{p}{m}\partial_x W
  \;=\; \int d\xi\thinspace V_W(x,\xi)\thinspace W(x, p - \xi),
\qquad
V_W(x,\xi) = \frac{1}{i\pi\hbar^2}\int dy\thinspace
   e^{-2i\xi y/\hbar}\bigl[V(x+y) - V(x-y)\bigr],
```

where $`y`$ is the **half** ket–bra separation, so a horizon
$`|x_{\rm ket} - x_{\rm bra}| \le L_c`$ is the window $`|y| \le L_c/2`$.

**Theorem O1 (flat envelope).** For every real $`V`$ with Fourier transform
$`\tilde V(k) = \int V(z)\thinspace e^{-ikz}dz`$,

```math
V_W(x,\xi) \;=\; \frac{2}{\pi\hbar^2}\thinspace
  \mathrm{Im}\negthinspace\left[e^{2i\xi x/\hbar}\thinspace
  \tilde V\!\left(2\xi/\hbar\right)\right].
```

*Proof.* Substitute $`z = x + y`$ in the first term and $`z = x - y`$ in the
second; each gives $`\tilde V(\pm 2\xi/\hbar)`$ times a pure phase, and
reality of $`V`$ pairs them into an imaginary part. $`\blacksquare`$

The dependence on $`x`$ is **pure phase**. The modulus of the kernel is
independent of position, for any potential, in any decomposition. A world a
hundred barrier widths downstream is not free: it is being struck at full
rate, and the correct free behaviour emerges only as a cancellation between
neighbouring $`\xi`$ channels whose fringe frequency $`2x/\hbar`$ grows
linearly with distance.

**The hologram analogy.** $`V_W`$ is a hologram of the potential. The object
is localised; the fringe pattern fills the plate at constant amplitude and
gets finer with distance from the object. Reconstructing the object needs the
whole plate. Note the combination $`2x/\hbar`$: **the far field and the
classical limit are the same Riemann–Lebesgue collapse**, which is why the
force emerges correctly in both.

Measured for a Gaussian barrier $`V = e^{-x^2/2\sigma^2}`$, $`\sigma = 0.5`$,
$`\hbar = 1`$ (Part A):

| $x$ | 0.5 | 1 | 2 | 4 | 8 | 16 |
|---|---|---|---|---|---|---|
| $`\max_\xi\lvert V_W\rvert`$ | 0.418 | 0.623 | 0.741 | 0.782 | 0.794 | 0.769 |
| $`V(x)`$ | $`6.1\times10^{-1}`$ | $`1.4\times10^{-1}`$ | $`3.4\times10^{-4}`$ | $`1.3\times10^{-14}`$ | $`2.6\times10^{-56}`$ | $`4.4\times10^{-223}`$ |

Fringe spacing measured against the prediction $`\pi\hbar/2x`$: $`0.392707`$
vs $`0.392699`$ at $`x = 4`$, $`0.196361`$ vs $`0.196350`$ at $`x = 8`$,
$`0.130901`$ vs $`0.130900`$ at $`x = 12`$.

**Corollary O1.1 (saturated budget).** The total vertex rate per world tends,
as $`|x| \to \infty`$, to

```math
R_\infty \;=\; \frac{2}{\pi^2\hbar}\int dk\thinspace\bigl|\tilde V(k)\bigr|
\;=\; \frac{4\thinspace V(0)}{\pi\hbar}\quad\text{for a positive Gaussian,}
```

the factor $`2/\pi`$ being the mean of $`|\sin|`$. Measured $`1.273066`$ at
$`x = 12`$ against $`4/\pi = 1.273240`$.

**Theorem O2 (signal against noise).** The first moment of the kernel is the
classical force exactly, $`\int \xi\thinspace V_W(x,\xi)\thinspace d\xi = -V'(x)`$,
while the zeroth and second moments of $`|V_W|`$ are flat. Hence the
signal-to-noise ratio of any signed-particle estimator of the force decays as
fast as $`V'`$ itself. Part B:

| $x$ | 0.5 | 1.0 | 1.5 | 2.0 | 2.5 | 3.0 |
|---|---|---|---|---|---|---|
| signal $`\lvert\int\xi V_W\rvert`$ | $`1.2131`$ | $`5.413\times10^{-1}`$ | $`6.665\times10^{-2}`$ | $`2.684\times10^{-3}`$ | $`3.727\times10^{-5}`$ | $`1.822\times10^{-7}`$ |
| $`\lvert V'(x)\rvert`$ | $`1.2131`$ | $`5.413\times10^{-1}`$ | $`6.665\times10^{-2}`$ | $`2.684\times10^{-3}`$ | $`3.727\times10^{-5}`$ | $`1.828\times10^{-7}`$ |
| rate $`\int\lvert V_W\rvert`$ | $`1.158`$ | $`1.273`$ | $`1.273`$ | $`1.273`$ | $`1.273`$ | $`1.273`$ |
| $`S/\sqrt{D}`$ | $`9.5\times10^{-1}`$ | $`4.8\times10^{-1}`$ | $`5.9\times10^{-2}`$ | $`2.4\times10^{-3}`$ | $`3.3\times10^{-5}`$ | $`1.6\times10^{-7}`$ |

The first moment tracks $`-V'`$ to $`10^{-9}`$ until it drops below the
quadrature floor near $`10^{-10}`$. The norm defect $`\int V_W\thinspace d\xi`$
is $`10^{-17}`$ throughout, as it must be: $`V_W`$ is odd in $`\xi`$.

**This is the escape problem, restated without a box.** Nothing is
misbehaving; no absorber can help, because the dynamics far from the
scatterer is already correct. What degrades is the *estimator*, which is
asked to resolve an ever-finer cancellation.

![Kernel diagnostics](https://raw.githubusercontent.com/billpage/wpmw/output/figures/open_position_space_kernel.png)

---

## 3. The coherence horizon

**Definition (H).** A *coherence horizon* $`L_c`$ is the postulate that a
conjugate ket–bra pair with separation $`|x_{\rm ket} - x_{\rm bra}| > L_c`$
is not instantiated. In the position-pair ladder's language it is a maximum
rung: $`|k| \le L_c/a`$. It is a bound in the **relative** coordinate, not a
wall in position space, and it is orthogonal to $`x`$:

```
   pair separation  y
        ^
  +Lc/2 |===============================   horizon
        |            +-----+
      0 |------------|  V  |------------->  midpoint x
        |            +-----+
  -Lc/2 |===============================
              |<---- active ---->|
                  supp V + Lc
```

**Definition (R), the reach.** The quantity that actually enters every formula
below is not $`L_c`$ itself but **half** of it. Give it a name:

```math
y_{\max} \;\equiv\; \frac{L_c}{2}
```

is the **reach** of a world — the greatest distance from its own position at
which it consults the potential.

Four remarks, because this quantity appears under several different names
across these notes and the factors of two are easy to lose.

1. *Why half.* The kernel above integrates against
   $`V(x + y) - V(x - y)`$, where $`x`$ is the pair **midpoint** and $`y`$ is
   the **half** separation: the ket sits at $`x + y`$ and the bra at
   $`x - y`$, so the separation between them is $`2y`$. A horizon
   $`|x_{\rm ket} - x_{\rm bra}| \le L_c`$ is therefore the window
   $`|y| \le L_c/2 = y_{\max}`$. A world of reach $`y_{\max}`$ samples $`V`$
   on the interval $`[x - y_{\max},\; x + y_{\max}]`$, of total length
   $`L_c`$, and on nothing else.

2. *Other names for the same thing.* In the momentum-conjugate variable
   $`s`$, where the potential term becomes a multiplication operator, the same
   quantity is the **stencil arm** $`y = \hbar s/2`$, and the horizon is the
   band $`|s| \le 2y_{\max}/\hbar`$. In the position-pair ladder it is the
   maximum rung, $`|k| \le L_c/a = 2y_{\max}/a`$. All three are the same
   bound written in three variables.

3. *Relation to the momentum quantum.* Combining
   $`\Delta p = \pi\hbar/L_c`$ from §1 with $`L_c = 2y_{\max}`$ gives

   ```math
   \Delta p \;=\; \frac{\pi\hbar}{2\thinspace y_{\max}} .
   ```

   Reach and momentum quantum are one parameter, not two: a world that
   consults $`V`$ only within $`y_{\max}`$ can resolve momentum transfers only
   down to $`\pi\hbar/(2y_{\max})`$. Shortening the reach coarsens the
   momentum lattice, and conversely.

4. *What it is not.* The reach is not a length scale of the potential, not a
   confinement of $`x`$, and not a property of any particular state. It is a
   postulate about which ket–bra pairs the model instantiates at all, so two
   worlds at the same $`x`$ with the same $`p`$ have the same reach.

The reach is used in this sense throughout
[`compensated_liouville_splitting.md`](compensated_liouville_splitting.md),
where it turns out to be the parameter that decides whether the classical
force can be carried by deterministic acceleration rather than by events.

**Proposition O3 (truncation is not absorption).** The windowed kernel is
still odd in $`\xi`$, so $`\int V_W^{L_c}(x,\xi)\thinspace d\xi = 0`$
identically and the signed world number is conserved exactly, for every
$`L_c`$ and every $`x`$. Measured defect $`\le 10^{-16}`$ at every entry of
the Part C tables, including the far-field entries where the rate is
$`10^{-88}`$.

**Proposition O4 (localisation).** All vertex activity is confined to
$`\mathrm{dist}(x, \mathrm{supp}V) \le L_c/2`$. Part C, jump rate
$`\int|V_W^{L_c}|`$:

| $x$ | 0 | 1 | 2 | 3 | 4 | 6 | 8 | 12 |
|---|---|---|---|---|---|---|---|---|
| $`L_c = 4`$ | 0 | $`1.46`$ | $`1.00`$ | $`1.3\times10^{-1}`$ | $`3.1\times10^{-4}`$ | $`1.1\times10^{-14}`$ | $`4.4\times10^{-32}`$ | $`9.8\times10^{-88}`$ |
| $`L_c = 8`$ | 0 | $`1.31`$ | $`1.41`$ | $`1.35`$ | $`1.00`$ | $`3.1\times10^{-4}`$ | $`1.1\times10^{-14}`$ | $`2.0\times10^{-56}`$ |
| $`L_c = 16`$ | 0 | $`1.28`$ | $`1.31`$ | $`1.28`$ | $`1.41`$ | $`1.31`$ | $`1.10`$ | $`1.2\times10^{-14}`$ |

The collapse edge tracks $`L_c/2`$ in every row. **This is the microdynamic
content of an "absorbing boundary": not a sponge that eats worlds, but the
statement that beyond the coherence horizon a world is exactly free.**

### 3.1 What the horizon costs

For a potential that decays, essentially nothing, once the window clears its
support. Relative $`L^2`$ error of the windowed kernel against the
untruncated one at $`x = 0.5`$, evaluated on the same nodes (Part C2):

| $L_c$ | 1 | 2 | 3 | 4 | 6 | 8 |
|---|---|---|---|---|---|---|
| kernel error | $`4.4\times10^{-1}`$ | $`2.8\times10^{-1}`$ | $`6.0\times10^{-2}`$ | $`4.2\times10^{-3}`$ | $`1.1\times10^{-6}`$ | $`3.9\times10^{-9}`$ |
| $`V`$ at window edge | $`1.0`$ | $`6.1\times10^{-1}`$ | $`1.4\times10^{-1}`$ | $`1.1\times10^{-2}`$ | $`3.7\times10^{-6}`$ | $`2.3\times10^{-11}`$ |

The error tracks the potential left outside the window and then hits the
quadrature floor. So there is no exactness-versus-locality tension for a
localised scatterer: **the horizon is free.**

### 3.2 Where the horizon fails

For a potential that does **not** decay, a sharp window cuts through a
non-zero integrand and the kernel acquires $`1/\xi`$ tails, so the rate budget
diverges logarithmically in the momentum cutoff. For
$`V = \cos(2\pi x/4)`$ at $`x = 0.7`$, whose exact budget is
$`2|\Gamma_q(x)| = 1.7820`$ (two deltas at $`\xi = \pm\pi\hbar/a`$), Part C3:

| $`\xi_{\max}`$ | 5 | 10 | 20 | 40 | 80 |
|---|---|---|---|---|---|
| windowed budget, $`L_c = 8`$ | $`3.627`$ | $`3.716`$ | $`3.759`$ | $`3.781`$ | $`3.792`$ |
| ratio to exact | $`2.04`$ | $`2.09`$ | $`2.11`$ | $`2.12`$ | $`2.13`$ |

Already twice the true rate at the coarsest cutoff, and still climbing. A
sharp horizon is the wrong tool for a crystal potential — which is exactly
the case §4 handles without any window at all.

---

## 4. Periodicity of $V$, not compactness of $x$

This is the useful half of the answer, and it is a different mechanism from
§3 in kind, not just in degree.

### 4.1 Two independent discreteness mechanisms

- **Kinematic (R/H).** Bounded ket–bra separation makes the momentum
  *argument* of $`W`$ discrete. It is a statement about which Wigner
  functions exist.
- **Dynamical (P).** Periodicity of $`V`$ makes the momentum *jumps*
  discrete. It is a statement about which transitions occur. It does not make
  $`p`$ discrete; it foliates a continuous $`p`$-axis into cosets that never
  mix.

They are logically independent, and each occurs without the other:

- **(H) without (P).** A localised barrier inside a horizon. Momenta lie on a
  lattice; every site is reachable; there is no conservation law.
- **(P) without (H).** A crystal potential on the open line. Momenta are
  continuous; each world is confined to its own coset forever.
- **Both.** The ring with a commensurate potential — the project's standard
  setting, and the reason the two mechanisms have not been distinguished. When
  $`a = L`$ the two formulas coincide numerically and the distinction is
  invisible.

### 4.2 The coset invariant

Let $`V`$ be $`a`$-periodic, $`V(x) = \sum_{j\ge1} V_j\cos(K_j x + \phi_j)`$
with $`K_j = 2\pi j/a`$. Mode $`j`$ couples the rows
$`p \mp \hbar K_j/2 = p \mp \pi\hbar j/a`$.

**Theorem O5 (coset invariant).** Every world's residue
$`\theta = p \bmod (\pi\hbar/a)`$ is constant for all time, on all of
$`\mathbb{R}`$, with no box and no horizon.

*Proof.* Streaming changes $`x`$ only. Every vertex changes a participant's
momentum by $`\pm\pi\hbar j/a`$ for some integer $`j`$, an integer multiple of
$`\pi\hbar/a`$. $`\blacksquare`$

Verified three ways (Part D). A free world on the open line, 2000 vertices,
arbitrary irrational offset: $`\max|\theta - \theta_0| = 1.1\times10^{-16}`$.
On a ring of $`L = 4a`$, with the state initialised on a single residue class:
leak into the other three classes is **$`0.0`$** — not small, zero — for every
class. And the decomposition is complete: the sum of the four sector runs
equals the direct run to **$`0.0`$** relative.

**Relation to Theorem A3 of the 4-D note.** A3 says a phase-space direction is
conserved if and only if it is orthogonal to every active mode wavevector —
an all-or-nothing statement about *directions*. Theorem O5 is its modular
companion: *along* an active direction there is still a residual conserved
quantity, namely $`p`$ modulo the jump quantum. A3 is the continuum part of
the invariant; O5 is the discrete part. Neither implies the other.

**Convention caveat.** Under the exchange-only reading in which every vertex
moves a participant by a *full* jump $`2\hbar K_j/2`$ rather than a half-jump,
the invariant reads $`p \bmod (2\pi\hbar/a)`$ instead. This is the same
open convention flagged in
[`../algorithm/multi_body_extension.md`](../algorithm/multi_body_extension.md)
§12 and by Theorem 1 of
[`phase_resonance_microdynamics.md`](phase_resonance_microdynamics.md)
(fundamental particles occupy even rows). Which modulus applies depends on
that convention; that the invariant exists and is exact does not.

![Sectors and melting](https://raw.githubusercontent.com/billpage/wpmw/output/figures/open_position_space_sectors.png)

### 4.3 The ring is a sector, not an approximation

Take a ring of circumference $`L = n\thinspace a`$ with an $`a`$-periodic
potential. Then only the modes $`q \in n\mathbb{Z}`$ are present, so every
jump spans a multiple of $`n`$ momentum cells and the generator is block
diagonal in the residue class of the momentum index mod $`n`$.

**Corollary O5.1.** The $`L = na`$ ring calculation is *exactly* the direct
sum of $`n`$ independent copies of the $`a`$-ring problem carrying different
offsets $`\theta`$. It is not a coarse approximation to open space; it is a
finite selection of that problem's sectors.

The practical consequence inverts the intuition behind the escape table of
[`../supplement/inverted_pair_barrier.md`](../supplement/inverted_pair_barrier.md)
§4.2. **For a periodic potential, enlarging the box does not refine anything.
It adds sectors.** The resolution within any one sector is fixed by $`a`$ and
is already exact. (For a *non*-periodic potential on a ring all modes are
present, $`n = 1`$, there is one sector, and enlarging $`L`$ genuinely does
refine $`\Delta p`$ — which is why the escape table is not wrong, only
specific to a non-periodic case.)

### 4.4 What the periodic route buys, and what it does not

**Buys.**

1. **Exactness.** No window, no truncation, no $`1/\xi`$ tails. The QLE on the
   open line with a periodic $`V`$ is represented without approximation.
2. **No boundary and no escape problem.** Worlds run to $`\pm\infty`$; nothing
   wraps, nothing is absorbed, nothing is discarded.
3. **A conserved quantum number per world.** $`\theta`$ is the microdynamic
   reading of Bloch quasimomentum, carried by an individual world rather than
   inferred from a band structure. That is a genuinely new object in this
   framework: a per-world label that survives every vertex.
4. **Zero-communication parallelism.** Sectors never exchange anything, so a
   run decomposes into independent jobs with no synchronisation. This matters
   most where the cost lives — the encounter loop of
   [`relational_pairing_and_carrier_lock.md`](relational_pairing_and_carrier_lock.md)
   and the pair bookkeeping of
   [`permanent_pairing_density_matrix.md`](permanent_pairing_density_matrix.md).
5. **The framework is unusually well matched to crystals.** The phase-space
   crystal lattice *is* the potential's reciprocal lattice, halved. For a
   solid-state problem the model's central structure is not a discretisation
   choice at all; it is the physics.

**Does not buy.**

1. **It is not asymptotically cheaper.** Representing a wavepacket of spatial
   extent $`\Lambda`$ still needs of order $`\Lambda/a`$ sectors, and the
   total work is unchanged. What disappears is the *error* from the box, not
   the cost.
2. **It does not apply to the current test case.** The inverted harmonic pair
   potential has no period. §3's horizon is the only route there, and by
   §3.1 that is fine, because $`-\tfrac{\mu}{2}\Omega^2 r^2`$ truncated to a
   scattering region is a decaying-support problem in the relevant sense once
   the asymptotic free region is reached.
3. **It removes the whole-state deterministic mesh reference.** $`\theta`$ is
   continuous, so one can mesh a sector but not the state. Validation has to
   be done sector by sector, which is a real change to the verification
   protocol of the demos.

### 4.5 Grid refinement multiplies sectors

The specification's option $`\Delta p = \pi\hbar/(KL)`$ is the same theorem
read backwards: mode $`q`$ then drives jumps of $`Kq`$ cells, so the lattice
splits into $`K`$ non-interacting sub-lattices. Part E, leak out of sector
zero after 30 Strang steps:

| $K$ | 1 | 2 | 3 | 4 |
|---|---|---|---|---|
| leak | $`0.0`$ | $`0.0`$ | $`0.0`$ | $`0.0`$ |

If the initial state lives on one coset — the normal case — then $`K-1`$ of
the sub-lattices are identically zero forever and the refinement is pure
waste. If it does not, the refinement is representing more of the
$`\theta`$-integral, which is a different and legitimate thing. Either way,
calling it "finer momentum resolution" is misleading.

### 4.6 Which mechanism covers which problem

| Problem | Periodic $V$? | Mechanism | Notes |
|---|---|---|---|
| cosine well (`demo_cosine_well_microdynamics.py`) | yes, $`a = L`$ | (P) or (R) | the degenerate case where they coincide |
| inverted harmonic pair barrier | no | (H) only | §4.4 item 2 |
| soft-core Coulomb pair, `multi_body_extension.md` §7 | no | (H) plus Ewald | Ewald's short-range part is the windowed piece |
| electron in a crystal, band problems | yes | (P) | the natural home of this framework |
| two particles with $`V(x_1 - x_2)`$ periodic | yes, in the relative mode | (P) | the anti-diagonal mode of the 4-D note |

---

## 5. Almost-periodic $V$: the crystal melts

Let $`V`$ have $`r`$ incommensurate periods. The set of momenta reachable
from $`p_0`$ is $`p_0 + \Lambda`$ with $`\Lambda`$ a $`\mathbb{Z}`$-module of
rank $`r`$ — a lattice for $`r = 1`$, dense in $`\mathbb{R}`$ for $`r \ge 2`$.
Smallest gap in the reachable set within $`|p| \le 3`$, as a function of the
number of vertices traversed (Part F, $`\kappa = 1, \sqrt2, \sqrt3`$):

| vertices | 2 | 4 | 6 | 8 |
|---|---|---|---|---|
| $`r = 1`$ | $`5.0\times10^{-1}`$ | $`5.0\times10^{-1}`$ | $`5.0\times10^{-1}`$ | $`5.0\times10^{-1}`$ |
| $`r = 2`$ | $`2.1\times10^{-1}`$ | $`8.6\times10^{-2}`$ | $`3.6\times10^{-2}`$ | $`3.6\times10^{-2}`$ |
| $`r = 3`$ | $`4.8\times10^{-2}`$ | $`2.3\times10^{-2}`$ | $`1.7\times10^{-3}`$ | $`1.7\times10^{-3}`$ |

For $`r = 1`$ the gap is pinned at $`\hbar\kappa/2`$ and never moves: a
genuine crystal. For $`r \ge 2`$ it falls with every vertex: a quasicrystal,
then a continuum.

The conservation law degrades continuously with it. Theorem O5's invariant
still exists set-theoretically for $`r \ge 2`$ — $`p \bmod \Lambda`$ is
conserved — but $`\Lambda`$ is dense, so $`\mathbb{R}/\Lambda`$ is not
Hausdorff and the invariant is not a continuous function of $`p`$. **The
crystal does not melt by a change in the rates; it melts by a loss of
topological closure in the invariant.** That is a sharper statement than
"the lattice becomes fine", and it is falsifiable: no amount of numerical
resolution recovers a usable sector label for $`r \ge 2`$.

---

## 6. What breaks: the sea on non-compact space

The four-action layer generalises to open position space cleanly, because it
is pointwise in $`x`$. The **sea-dressed** layer does not, for a reason worth
recording rather than papering over.

The crystal shift $`W' = W + 2/h`$ is a *constant* phase-space density. On a
ring with a momentum cutoff its total is finite,
$`N_{\rm sea} \sim \nu\thinspace L\thinspace P_{\max}\thinspace(2/h)`$. On
$`\mathbb{R}\times\mathbb{R}`$ it is infinite. The sea-dressed ontology is
therefore tied to a finite phase-space volume in a way the four-action layer
is not.

The likely repair is a per-unit-length (grand-canonical) formulation, and the
order parameter $`Z_r`$ of
[`relational_pairing_and_carrier_lock.md`](relational_pairing_and_carrier_lock.md)
Theorem R4 is already per-cell and therefore already intensive, so it should
survive unchanged. What does not obviously survive is anything that
references a **total**: global neutrality bookkeeping, postulate (S) as
stated, and the resource arithmetic of
[`permanent_pairing_density_matrix.md`](permanent_pairing_density_matrix.md)
§5.3 — including the $`\approx 770\times`$ shortfall, which is a ratio of
totals and would have to be restated as a ratio of densities before it means
anything on open space. This note does not settle it; see open item 2.

---

## 7. An absorbing potential is the wrong object

The one candidate with the right name is a complex absorbing potential,
$`V \to V - i\Gamma`$, giving

```math
\dot\rho \;=\; -\frac{i}{\hbar}[H_0,\rho] \;-\; \frac{1}{\hbar}\{\Gamma,\rho\}.
```

Its Wigner image costs no new lattice geometry. With
$`x = X + y/2`$, $`x' = X - y/2`$ and $`\Gamma = e^{i\kappa x}`$,

```math
\Gamma(x) \pm \Gamma(x') \;=\; e^{i\kappa X}
   \bigl(e^{i\kappa y/2} \pm e^{-i\kappa y/2}\bigr),
```

and multiplying $`\rho`$ by $`e^{\mp i\kappa y/2}`$ shifts the Wigner momentum
by $`\pm\hbar\kappa/2`$. So the **commutator takes the difference** of the two
offset rows and the **anticommutator takes their sum**, on the very same
offsets. Verified on an $`M = 32`$ ring with a random pure state, $`q = 3`$
(Part G): commutator against $`\sin`$-weighted difference, error
$`1.07\times10^{-15}`$; anticommutator against $`\cos`$-weighted sum, error
$`1.72\times10^{-15}`$.

An absorber is therefore the existing stencil with one sign rule flipped. It
is still the wrong object, for three reasons:

1. **It destroys signed number by construction.** The sum stencil has non-zero
   column sums. That is what absorption means, but it also means the exact
   conservation that Proposition O3 gives for free is thrown away.
2. **It absorbs everywhere.** Theorem O1 applies to $`\Gamma`$ as much as to
   $`V`$: the absorber's own kernel has a flat envelope, so a CAP placed
   "at the edges" is active in the interaction region too, cancelling only in
   the signed mean. It adds noise exactly where the signal is.
3. **It is not the QLE.** The horizon of §3 is a restriction of the state
   space; a CAP is a different, non-unitary equation. Only one of the two can
   be checked against a closed form.

---

## 8. Consequences for existing documents

Cross-references added by the patch that carries this note, so that the
distinction of §1 is visible from wherever a reader meets $`\Delta p`$:

- `../algorithm/phase_space_crystal_lattice_algorithm.md` §1 — a note that the
  momentum quantum has three possible sources, and that the $`K`$-refinement
  option interleaves sectors rather than resolving them.
- `../algorithm/phase_alignment_microdynamics_algorithm.md` §1.2 — the
  "fixed by the ring circumference" sentence, qualified.
- `../algorithm/multi_body_extension.md` §2 — the per-axis
  $`\Delta p_a = \pi\hbar/L_a`$ choice, cross-referenced to §4 here.
- `../supplement/inverted_pair_barrier.md` §7 open item 4 — retracted and
  replaced.
- `fourd_microdynamics.md` — Theorem A3 cross-referenced to Theorem O5.
- `position_pair_ladder.md` — the rung index cross-referenced to the horizon.

---

## 9. Open items

1. **A tapered horizon.** §3.2 shows a sharp window is wrong for a
   non-decaying $`V`$. A taper kills the $`1/\xi`$ tails but destroys the
   exact momentum lattice, since the windowed kernel is no longer a finite
   Fourier series in $`y`$. Whether a taper exists that preserves both is
   unresolved.
2. **The sea per unit length.** §6. Restating global neutrality, postulate
   (S), and the storage arithmetic as densities.
3. **Half-jump or full-jump**, and hence whether Theorem O5's modulus is
   $`\pi\hbar/a`$ or $`2\pi\hbar/a`$. Ties to `multi_body_extension.md` §12.
4. **Energy under the horizon.** Proposition O3 proves signed number is
   conserved exactly by truncation. Energy is not checked, and there is no
   reason to expect it to be exact.
5. **Sampling the $`\theta`$-integral.** Sectors are independent, but a
   wavepacket needs many of them. Whether importance sampling in $`\theta`$
   is possible, and whether the sea must be instantiated separately per
   sector, are both open.
6. **Two horizons in the two-body problem.** For $`V(x_1 - x_2)`$ the natural
   horizon is on the relative coordinate. Whether the centre-of-mass
   coordinate needs one at all, and what a mixed pair separation means in the
   four-dimensional phase space of `fourd_microdynamics.md`, is unexamined.

---

## 10. Numerical verification summary

All from `src/demo_open_position_space.py`.

| Part | Claim | Result |
|---|---|---|
| A | flat envelope of $`V_W`$ | $`0.418 \to 0.794`$ while $`V`$ falls by $`10^{-223}`$ |
| A | fringe spacing $`\pi\hbar/2x`$ | agreement to $`2\times10^{-5}`$ |
| A | $`R_\infty = 4V(0)/\pi\hbar`$ | $`1.273066`$ vs $`1.273240`$ |
| B | first moment $`= -V'(x)`$ | to $`10^{-9}`$ above the quadrature floor |
| B | norm defect | $`\le 4\times10^{-15}`$ |
| C | activity edge at $`L_c/2`$ | collapse in every row |
| C2 | horizon cost inside supp $`V`$ | tracks $`V`$ at the window edge |
| C3 | budget divergence, periodic $`V`$ | $`2.04 \to 2.13`$ times exact, logarithmic |
| D | coset invariant, open line | $`1.1\times10^{-16}`$ over 2000 vertices |
| D | sector leak on the ring | $`0.0`$ |
| D | sector sum equals direct run | $`0.0`$ relative |
| E | $`K`$-refinement leak | $`0.0`$ for $`K = 1,2,3,4`$ |
| F | reachable-set gap | pinned at $`0.5`$ for $`r=1`$; falls for $`r\ge2`$ |
| G | commutator $`\to`$ difference | $`1.07\times10^{-15}`$ |
| G | anticommutator $`\to`$ sum | $`1.72\times10^{-15}`$ |

---

## 11. Sources

- [`../algorithm/phase_space_crystal_lattice_algorithm.md`](../algorithm/phase_space_crystal_lattice_algorithm.md)
  §1, §3b — the momentum quantum and the mode stencil.
- [`../supplement/inverted_pair_barrier.md`](../supplement/inverted_pair_barrier.md)
  §4.2, §7 — the escape problem and the open item this note answers.
- [`fourd_microdynamics.md`](fourd_microdynamics.md) — Theorem A3.
- [`position_pair_ladder.md`](position_pair_ladder.md) §1 — the rung index.
- The semi-discrete Wigner formalism, in which a finite coherence length is
  the standard route to a discrete momentum space for signed-particle Monte
  Carlo, is due to Nedjalkov, Ellinghaus, Sellier, Selberherr and
  collaborators; Jacoboni and Bordone give the finite-coherence-length
  transport equation. See `../../references/bibliography.md`.
- Generating script: `src/demo_open_position_space.py`. Figures published to
  the `output` branch as `figures/open_position_space_*.png`.
