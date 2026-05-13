# Phase-Space Crystal-Lattice Model in 4D Phase Space:
# Two Particles in 1D vs One Particle in 2D

> **Provenance.** This document is a companion to
> `docs/supplement/phase_space_crystal_lattice_supplement.md`. The single-particle
> 1+1-D rules established there are extended here to the smallest non-trivial
> generalisations — two particles in one spatial dimension and one particle in
> two spatial dimensions — both of which have a four-dimensional joint phase
> space. The two cases share an algorithmic skeleton but differ sharply in
> what is conserved per event and in the structure of the mediator picture.
> The diagrams below are the natural analog of David Cyganski's
> *Wigner Collisions Diagram* (Sozi presentation) for the new geometry.

---

## 1. Scope

The 1+1-D supplement derived a single jump rule:

$$\Delta W \\\;=\\\; -\\, \frac{V_p}{\hbar}\\, \Delta t\\, \sin\negthinspace\theta\\, 
\bigl[W(p+\tfrac{\pi\hbar}{L}) - W(p-\tfrac{\pi\hbar}{L})\bigr],
\qquad \theta \equiv \tfrac{2\pi q\\,  x}{L} + \phi_q .$$

driving discrete momentum jumps of magnitude $\hbar k_q$ at rate
$|\Gamma_q(x)| = |V_q/\hbar||\sin\theta|$. Two systems exist whose joint phase
space is one step larger:

- **2p/1D** — two distinguishable particles in 1 spatial dimension,
  $(x_1, p_1, x_2, p_2) \in \mathbb{R}^4$, interacting via a pair potential
  $V_2(x_1 - x_2)$ and optionally an external $V_{\mathrm{ext}}(x)$.
- **1p/2D** — one particle in 2 spatial dimensions,
  $(x, y, p_x, p_y) \in \mathbb{R}^4$, in an external potential $V(x, y)$.

The 1p/2D case is purely a vectorial restatement of the 1+1-D rule. The
2p/1D case introduces a genuinely new ingredient — **correlated jumps**
between two particles — which is the prototype of every multi-body
interaction in the algorithm specification at
`docs/algorithm/multi_body_extension.md`. The present document derives both
rules from the multi-dimensional QLE, presents them side by side, and
exhibits the geometric correspondence (and its limits) between the two.

---

## 2. Setup: 4D joint phase space

The joint Wigner function $W(\mathbf{X}, \mathbf{P}, t)$ in either case is a
real function on $\mathbb{R}^4$ that obeys the multi-dimensional Quantum
Liouville Equation. Writing $\mathbf{X}$ and $\mathbf{P}$ for the 2-component
configuration and momentum coordinates:

| | configuration $\mathbf{X}$ | momentum $\mathbf{P}$ |
|---|---|---|
| **2p/1D** | $(x_1, x_2)$ | $(p_1, p_2)$ |
| **1p/2D** | $(x, y)$ | $(p_x, p_y)$ |

The joint Wigner function is bounded by $|W| \le (2/h)^2$, so the
crystal-lattice shift to a non-negative *positon-count* density is
$W' = W + (2/h)^2$.

![Both systems share a 4-dimensional joint phase space, sliced naturally in different ways](https://raw.githubusercontent.com/billpage/wpmw/output/figures/microdynamics_4d_layouts.png)

The free-streaming sub-step is identical in structure for both cases:
each $\mathbf{X}$-component advances according to the corresponding
$\mathbf{P}$-component divided by the mass. The differences live entirely
in the potential-driven sub-step.

---

## 3. The single-particle 1+1-D rule, restated

For the rest of this document, the "atomic" jump element is the contribution
of a single Fourier mode of the potential to the QLE:

$$V_q\\, \cos(k\\, x + \phi)
\\\;\longrightarrow\\\;
\Gamma_q(x) = -\frac{V_q}{\hbar}\\, \sin(k\\, x + \phi), \qquad k = \tfrac{2\pi q}{L}.$$

The rule reads:

> For each positon at $(x, p)$, with probability $|\Gamma_q(x)|\\, \Delta t$
> per timestep, transfer one positon between cells $(x, p-\hbar k)$ and
> $(x, p+\hbar k)$ in the direction set by $\mathrm{sgn}\\, \Gamma_q(x)$. The
> mediator at $(x, p)$ is unchanged.

(See §6 of the 1+1-D supplement for the derivation and §6.3 for the
sign-convention check.) Every rule below is built from this atom by changing
*what kind of object* a Fourier mode lives on, and *which* particles are
affected by the resulting jump.

---

## 4. 1p/2D: one particle, vectorial Fourier modes

### 4.1 Fourier expansion of $V(x, y)$

The 2-D external potential admits a 2-D Fourier expansion on the periodic box
$[0, L_x] \times [0, L_y]$:

$$V(x, y) \\\;=\\\; V_0 \\\;+\\\; \sum_{\vec q \neq 0}\\\; V_{\vec q}\\, \cos(\vec k_{\vec q}\negthinspace\cdot\negthinspace\vec x + \phi_{\vec q}),
\qquad \vec k_{\vec q} = \Bigl(\tfrac{2\pi q_x}{L_x},\\, \tfrac{2\pi q_y}{L_y}\Bigr),
\quad \vec q \in \mathbb{Z}^2.$$

Each mode $\vec q$ is now indexed by an integer vector, and its wavevector
$\vec k_{\vec q}$ has an arbitrary direction in the $(q_x, q_y)$ lattice.

### 4.2 Per-mode rate and jump rule

The local rate amplitude is the direct vectorial analog of the 1+1-D form:

$$\Gamma_{\vec q}(x, y) \\\;=\\\; -\\, \frac{V_{\vec q}}{\hbar}\\, \sin(\vec k_{\vec q}\negthinspace\cdot\negthinspace\vec x + \phi_{\vec q}).$$

The jump rule, in the world-ensemble (continuous-momentum) representation,
is the obvious vectorial lift:

> For each positon at $(\vec x, \vec p)$ in each world, with probability
> $|\Gamma_{\vec q}(\vec x)|\\, \Delta t$ per mode per step,
>
> $$\vec p \\\;\longrightarrow\\\; \vec p \\\;+\\\; \mathrm{sgn}\bigl(\Gamma_{\vec q}(\vec x)\bigr)\\, \hbar\\, \vec k_{\vec q}.$$

Both components of $\vec p$ shift in the same event — the kick is genuinely
vectorial. No part of $\vec p$ is conserved by the rule (momentum is
transferred between particle and field, exactly as for any single particle in
an external potential).

In the mesh-density form, the equivalent update is

$$W(\vec x, \vec p, t+\Delta t) \\\;=\\\; W(\vec x, \vec p, t) \\\;+\\\; \Delta t \sum_{\vec q \neq 0}\\\;
\Gamma_{\vec q}(\vec x)\bigl[W(\vec x, \vec p + \hbar\vec k_{\vec q}) - W(\vec x, \vec p - \hbar\vec k_{\vec q})\bigr].$$

### 4.3 Derivation

The QLE force term in 2-D is

$$\partial_t W \\\;\supset\\\; +\nabla_{\negthinspace\vec x}V(\vec x)\cdot\nabla_{\negthinspace\vec p}W .$$

For the single Fourier mode $V_{\vec q}\cos(\vec k_{\vec q}\negthinspace\cdot\negthinspace\vec x + \phi_{\vec q})$,

$$\nabla_{\negthinspace\vec x}V = -V_{\vec q}\\, \vec k_{\vec q}\\, \sin(\vec k_{\vec q}\negthinspace\cdot\negthinspace\vec x + \phi_{\vec q})$$

and the directional derivative
$\vec k_{\vec q}\cdot\nabla_{\negthinspace\vec p}W$ is approximated by the symmetric
centred difference at the photon-momentum scale $\hbar\vec k_{\vec q}$:

$$\vec k_{\vec q}\cdot\nabla_{\negthinspace\vec p}W \\\;\approx\\\; \frac{W(\vec p+\hbar\vec k_{\vec q}) - W(\vec p-\hbar\vec k_{\vec q})}{2\hbar},$$

(this is the same finite-difference identification at the photon-momentum
scale used in §6.1 of the 1+1-D supplement, applied along the direction
$\hat k_{\vec q}$). Combining gives the §4.2 update directly.

---

## 5. 2p/1D: two particles, correlated jumps

### 5.1 Pair potential and Fourier expansion

A pair potential depending only on the relative coordinate $r = x_1 - x_2$
expands as

$$V_2(r) \\\;=\\\; V_2^{(0)} \\\;+\\\; \sum_{q\neq 0}\\\; V^{(2)}_q\\, \cos(k_q\\, r + \phi^{(2)}_q),
\qquad k_q = \tfrac{2\pi q}{L}.$$

In the joint configuration space $(x_1, x_2)$ this is

$$V_2(x_1 - x_2) \\\;=\\\; V_2^{(0)} \\\;+\\\; \sum_{q\neq 0}\\\; V^{(2)}_q\\, \cos(q\\, k_1\\, x_1 \\\;-\\\; q\\, k_1\\, x_2 + \phi^{(2)}_q),$$

with $k_1 = 2\pi/L$. As a joint 2-D Fourier series, $V_2$ has support **only
on the anti-diagonal** $q_1 = -q_2$ — a 1-D subset of the 2-D wavevector
lattice $\mathbb{Z}^2$.

![Fourier-mode support: 2p/1D pair potential lives on a 1D line in joint wavevector space, while 1p/2D modes can fill all of Z^2](https://raw.githubusercontent.com/billpage/wpmw/output/figures/microdynamics_4d_fourier_modes.png)

That dimensional reduction in $\vec q$-space is the algebraic shadow of the
centre-of-mass / relative-coordinate separation, examined in §7 below.

### 5.2 Per-mode rate and correlated jump rule

The local rate amplitude depends on the relative coordinate only:

$$\Gamma^{(2)}_q(r_{12}) \\\;=\\\; -\\, \frac{V^{(2)}_q}{\hbar}\\, \sin(k_q\\, r_{12} + \phi^{(2)}_q),
\qquad r_{12} = x_1 - x_2.$$

The jump rule is **correlated**: both particles' momenta change in the same
event, with equal and opposite kicks:

> For each world, each Fourier mode $q$ of $V_2$, with probability
> $|\Gamma^{(2)}_q(r_{12})|\\, \Delta t$:
>
> $$(p_1,\\, p_2) \\\;\longrightarrow\\\;
>   \bigl(p_1 + \mathrm{sgn}(\Gamma^{(2)}_q)\\, \hbar k_q,\\\;
>         p_2 - \mathrm{sgn}(\Gamma^{(2)}_q)\\, \hbar k_q\bigr).$$

The total momentum $P = p_1 + p_2$ is exactly conserved per event. In
$(p_1, p_2)$ phase-space coordinates, every jump is a displacement along the
anti-diagonal direction; jumps for different $q$ are along the same line, but
with different lengths $2\hbar k_q$ (in $(p_1, p_2)$ metric).

![Per-event momentum-space displacement: 2p/1D jumps confined to the anti-diagonal; 1p/2D jumps in arbitrary directions](https://raw.githubusercontent.com/billpage/wpmw/output/figures/microdynamics_4d_jumps.png)

In mesh-density form (when a joint $(x_1, p_1, x_2, p_2)$ grid is available),
the equivalent update is the joint finite difference

$$\partial_t W^{(2)} \\\;\supset\\\; \sum_{q\neq 0} \Gamma^{(2)}_q(r_{12})
\bigl[\\, W^{(2)}(\dots, p_1+\hbar k_q,\\,  p_2-\hbar k_q,\dots) - W^{(2)}(\dots, p_1-\hbar k_q,\\,  p_2+\hbar k_q,\dots)\\, \bigr].$$

### 5.3 Derivation

The joint QLE for two particles with pair potential $V_2(x_1-x_2)$ contains
the Moyal-bracket contribution

$$\partial_t W^{(2)} \\\;\supset\\\;
\frac{2}{\hbar}\\, V_2(r_{12})\\, 
\sin\negthinspace\Bigl(\tfrac{\hbar}{2}\bigl(\overleftarrow\partial_{x_1} - \overleftarrow\partial_{x_2}\bigr)\bigl(\overrightarrow\partial_{p_1} - \overrightarrow\partial_{p_2}\bigr)\Bigr)\\, W^{(2)},$$

where only the *difference* of position-derivatives appears (because $V_2$
depends only on the difference $r_{12}$). For a single Fourier mode
$V^{(2)}_q\cos(k_q\\, r_{12} + \phi^{(2)}_q)$, the factor
$\overleftarrow\partial_{x_1} - \overleftarrow\partial_{x_2}$ acting on the
cosine produces $-2k_q\sin(\dots)$ (one $-k_q$ from each derivative), and the
remaining $\sin$ of the Moyal bracket becomes a finite-difference operator
acting in $(p_1, p_2)$ along the joint direction
$\hat e_{p_1} - \hat e_{p_2}$. The result is the §5.2 rule, with the
factor-of-2 from the derivative absorbing the half-step that appears in the
1+1-D mesh form.

(For a fully detailed derivation see §7 of
`docs/algorithm/multi_body_extension.md`; the discussion there extends
verbatim to $d$ spatial dimensions and $N$ particles.)

---

## 6. Side-by-side comparison

| | **2p/1D** | **1p/2D** |
|---|---|---|
| Joint phase space | $(x_1, p_1, x_2, p_2)$ | $(x, y, p_x, p_y)$ |
| Free-stream rule  | $\dot x_i = p_i/m$ per particle | $\dot x = p_x/m$, $\dot y = p_y/m$ |
| Potential | $V_2(x_1 - x_2)$ (pair) + opt. $V_\mathrm{ext}(x_i)$ | $V(x, y)$ (external) |
| Fourier index | scalar $q\in\mathbb{Z}$ for $V_2$ | vector $\vec q\in\mathbb{Z}^2$ |
| Active mode support in joint $\vec q$-space | 1-D line $q_1 = -q_2$ | full 2-D lattice |
| Rate phase depends on | $r_{12} = x_1 - x_2$ | $\vec k_{\vec q}\negthinspace\cdot\negthinspace\vec x$ |
| Particles kicked per event | both, oppositely | the single particle, vectorially |
| Direction of $\Delta\vec p_{\text{joint}}$ | anti-diagonal in $(p_1,p_2)$ | arbitrary in $(p_x,p_y)$ |
| Total momentum | conserved per event | not conserved |
| Mediator-picture analog | virtual quantum *exchanged* between two particles | virtual quantum *absorbed/emitted* by one particle |
| Cost per Fourier mode per step | one event per pair per world | one event per particle per world |

The free-streaming sub-step and the *form* of the rate amplitude
$\Gamma_q \propto -(V_q/\hbar)\sin(\text{phase})$ are common to both; the
distinction is entirely in *which* momenta are kicked and the constraint
that the geometry of the potential imposes on their joint motion.

---

## 7. The hidden 1+1-D character of 2p/1D

The most striking feature of 2p/1D is that the *entire* joint dynamics
separates exactly into a free centre-of-mass plus a 1+1-D problem in the
relative coordinate. Define

$$X = \tfrac{x_1 + x_2}{2}, \quad r = x_1 - x_2; \qquad
P = p_1 + p_2, \quad p_{\mathrm{rel}} = \tfrac{p_1 - p_2}{2}.$$

The kinetic energy separates as $T = P^2/(2M) + p_{\mathrm{rel}}^2/(2\mu)$
with $M = 2m$, $\mu = m/2$, and $V_2$ depends only on $r$. The joint Wigner
function factorises whenever $W^{(2)}$ does, and the QLE Moyal series in the
$(X, P)$ and $(r, p_{\mathrm{rel}})$ sub-spaces are independent (no cross
terms).

The crystal-lattice rule respects this separation exactly: in $(P, p_{\mathrm{rel}})$
coordinates, every jump lies along the $p_{\mathrm{rel}}$ axis, so $P$ is not
just conserved on average — it is **untouched per event**.

![Same jump rule in original (p1, p2) basis versus rotated (P, p_rel) basis - the latter exhibits the separation manifestly](https://raw.githubusercontent.com/billpage/wpmw/output/figures/microdynamics_4d_com_relative.png)

This is the algebraic justification for the table entry "Active mode support
is 1-D" in §6: the 2p/1D system is, in a precise sense, a 1+1-D system (in
the relative coordinate) *plus* a decoupled free particle. The 1p/2D system
has no such reduction in general — it does so only when $V(x,y)$ separates
as $V_1(x) + V_2(y)$, in which case it likewise becomes two decoupled 1+1-D
problems.

> **Practical consequence.** A 2p/1D simulation has two equivalent
> implementations: (i) run the algorithm in $(x_1, p_1, x_2, p_2)$ coordinates
> as specified above; (ii) transform to $(X, P, r, p_{\mathrm{rel}})$ and run
> two independent 1+1-D evolutions — one trivial (free $P$), one using the
> existing single-particle crystal-lattice code with mass $\mu$ and "external"
> potential $V_2(r)$. Path (ii) is what makes a 2p/1D problem with no
> external $V_{\mathrm{ext}}$ no harder than a 1+1-D problem; path (i) is the
> form that extends to a non-separable $V_{\mathrm{ext}}(x_1) + V_{\mathrm{ext}}(x_2)$
> and to higher $N$.

---

## 8. Crystal-lattice mediator picture

The crystal-lattice / static-negaton-background reading of §4 of the 1+1-D
supplement extends in both directions, with one important change in the
mediator picture for 2p/1D.

For 1p/2D, the picture is essentially the same as in 1+1-D — a single
mediator positon induces jumps in surrounding positons — but now the
"surrounding positons" live in a 2-D momentum lattice and the jump is along
the vector $\hbar\vec k_{\vec q}$.

For 2p/1D, the single-mediator interpretation does not literally apply: a
pair interaction involves *two* participating positons (one per particle).
The cleanest reading is **per-world pair-wise virtual-quantum exchange**: a
quantum of momentum $\hbar k_q$ is transferred between particle 1 and
particle 2 within a single world, at a rate set by the relative-coordinate
phase of mode $q$. This is structurally the leading vertex of QED at the
non-relativistic-QLE level (see the discussion of the photon-momentum
identification in §5 of the 1+1-D supplement, which now carries over with
the photon momentum *between two real particles* instead of between a
particle and a lattice neighbour).

![Crystal-lattice mediator picture: a virtual quantum is exchanged between the two particles in 2p/1D; the single particle absorbs a vectorial kick in 1p/2D](https://raw.githubusercontent.com/billpage/wpmw/output/figures/microdynamics_4d_starburst.png)

This is also the clearest place in which "the single rule" (one Fourier mode
drives one kind of momentum jump) shows its multi-particle behaviour: the
*correlation* between the two particles is built into the rule itself, not
imposed by hand. Total momentum conservation per event is a property of the
QLE's Moyal series for pair potentials — the algorithm inherits it
automatically.

The crystal-lattice shift to a non-negative positon density carries through
without modification: the joint Wigner function is bounded by $(2/h)^2$, so
$W' = W + (2/h)^2$ is the non-negative representative. The QLE remains
invariant under the shift because every term in the Moyal series still
contains a derivative of $W^{(2)}$ (the joint version of the argument in
§4.4 of the 1+1-D supplement).

---

## 9. Combining external and pair potentials in 2p/1D

A common practical case is two interacting particles each subject to the
same external potential — for example, two electrons in a 1-D trap. The
combined update applies the external-potential rule (§4 of this document
specialised to $d=1$, identical in form to the 1+1-D rule) to each particle
independently, *and* the pair-interaction rule of §5 once per unordered
pair. Each Fourier mode of each potential contributes one event per timestep
per applicable particle (or particle pair); events from different modes are
independent.

The split-operator skeleton is:

```
for each timestep:
    free-stream each particle:   x_i ← x_i + p_i dt/m
    for each ext mode q:
        for each particle i:
            jump p_i with prob |Γ_q^ext(x_i)| dt
    for each pair mode q:
        for each pair (i, j):
            correlated jump (p_i, p_j) with prob |Γ_q^(2)(r_ij)| dt
```

Strang splitting (one half-step of free-stream, full potential, second
half-step) recovers second-order accuracy in $\Delta t$ as in the 1+1-D
case.

---

## 10. Validation paths

Both cases admit cheap, mesh-based validation: 4-D phase space at
$M_x = M_p = 32$ per axis is $\sim 10^6$ cells, still tractable on a
laptop.

For **1p/2D**, the natural test is a 2-D anisotropic harmonic oscillator
or an isotropic separable potential. The ground-state Wigner function is
Gaussian and stationary under the QLE; verifying that the discrete-jump
update preserves it (up to time-step error) generalises the 1+1-D
test in `src/demo_qho_ground_state.py`. A more searching test is a
2-D cosine well analogous to `src/demo_cosine_well_microdynamics.py`,
including a non-separable cross term $V_{\mathrm{cross}}\cos(k_x x + k_y y)$.

For **2p/1D**, the analytic anchor is the two-particle harmonic
oscillator with harmonic coupling, which separates exactly into two
oscillators in the relative-coordinate and centre-of-mass coordinates
and has a Gaussian-product Wigner ground state. A second validation —
specifically of the correlated-jump rule — is a 2p/1D cosine pair
potential $V_2(r) = V_p\cos(2\pi r/L)$, comparing the joint mesh
evolution against the COM/relative decomposition: the two should agree
to within Monte-Carlo noise on every marginal.

A regression test specific to the correlated rule (the 2p/1D analog of
`src/sign_convention_check.py`) is to confirm that an asymmetric initial
state $W^{(2)}_0(x_1, p_1, x_2, p_2)$ with non-zero initial $P$ keeps
$P$ exactly constant under every event, while $p_{\mathrm{rel}}$ varies
in the way the 1+1-D rule predicts for the equivalent reduced-mass
problem in $r$.

---

## 11. Connection to the algorithm specification

The rules derived above are the simplest case ($d = 1, N = 2$ and
$d = 2, N = 1$) of the general specification in
`docs/algorithm/multi_body_extension.md`. The mappings are direct:

| Multi-body spec (general $d, N$) | This document |
|---|---|
| §3.3, external $V$ in $d$ dimensions | §4 with $d = 2$, $N = 1$ |
| §6.2, external $V_\mathrm{ext}$ per particle | §4 with $d = 1$ applied to each of two particles |
| §6.3, two-body $V_2(\vec r_{ij})$ in $d$ dimensions | §5 with $d = 1$, $N = 2$ |
| §7, derivation of correlated jumps from joint Moyal bracket | §5.3 |

The specification's open questions about world-ensemble convergence,
identical-particle symmetrisation, and the fermionic sign problem are not
forced on us at the 2p/1D level (distinguishable particles, low statistical
load on smooth Wigner functions), and so this case offers a clean place to
develop and validate the correlated-jump machinery before scaling up.

---

## 12. Sources

- David Cyganski, *Extended Fokker–Planck Eq. and the QLE V2* (project memo).
- David Cyganski, *Wigner Collisions Diagram* (Sozi presentation).
- `docs/supplement/phase_space_crystal_lattice_supplement.md` — 1+1-D
  redraft with sign-convention correction.
- `docs/algorithm/multi_body_extension.md` — full $N$-body $d$-dimensional
  specification of which this document is the 4-D special case.
- `docs/analysis/phase_space_crystal_lattice_review.md` — cross-referenced
  page-and-equation review of the V2 memo.
- Generating script: `src/gen_microdynamics_4d_figures.py`. Figures are
  published to the `output` branch as
  `figures/microdynamics_4d_*.png`.
