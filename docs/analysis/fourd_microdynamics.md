# Four-Dimensional Microdynamics: Two Particles in 1D and One Particle in 2D, Through the Four-Action, Sea-Dressed and Phase-Alignment Layers

## 0. Status and provenance

This note takes the derivation ladder of `docs/analysis/` — built entirely in
$1+1$ dimensions — and asks what survives, what generalises, and what breaks
when the joint phase space is four-dimensional. The two smallest cases are the
ones treated kinematically in
[`../supplement/phase_space_crystal_lattice_4d_supplement.md`](../supplement/phase_space_crystal_lattice_4d_supplement.md):

- **2p/1D** — two distinguishable particles on a line,
  $(x_1, p_1, x_2, p_2)$, with a pair potential $V_2(x_1 - x_2)$ and
  optionally an external $V_{\mathrm{ext}}$;
- **1p/2D** — one particle in the plane, $(x, y, p_x, p_y)$, in an
  external $V(x,y)$.

The supplement derived the *jump rules*. What it did not do is push the three
microdynamic layers through: the four-action decomposition
([`four_rule_microdynamics_equivalence.md`](four_rule_microdynamics_equivalence.md)),
its sea-dressed realisation
([`sea_dressed_microdynamics.md`](sea_dressed_microdynamics.md)), and the
phase-alignment vertex
([`phase_alignment_microdynamics.md`](phase_alignment_microdynamics.md)).
Doing so is the content of §§2–4. §5 answers a specific question (B. Page,
August 2026) about two particles and a harmonic potential, and §6 proposes a
ranked set of test cases.

Derivation and implementation were developed jointly with Claude (Anthropic),
August 2026. Companion code: `src/demo_fourd_microdynamics.py`. Every
quantitative claim below is an output of that demo (§7).

### 0.1 Corrections recorded here

Three earlier statements are corrected or sharpened.

1. **The 4-D supplement's §6 comparison table is organised along the wrong
   axis.** It contrasts 2p/1D with 1p/2D as though the particle count decided
   what is conserved and whether the mediator is "exchanged" or "absorbed".
   Theorem A3 below shows the operative variable is the *mode wavevector*:
   a direction is conserved iff it is orthogonal to every active
   $\vec k_{\vec q}$. Since 2p/1D with a pair potential and no external field
   *is* 1p/2D with $V(x,y) = V_2(x-y)$ — the same partial differential
   equation, verified to $0$ in Part C — the two columns of that table cannot
   be describing different dynamics. §2.4 gives the corrected table.
2. **§6.1 of the four-rule note over-states the hop channel's character.**
   It asserts that "the hop channel cannot be a particle–particle exchange
   without violating momentum conservation, and should be understood as
   photon absorption/emission." That is true for external modes only. For a
   pair-potential mode the hop conserves the total particle momentum exactly
   (Part B), and *is* a particle–particle exchange — of a full photon
   $\hbar k_q$ rather than the focus channel's half. The particle–field /
   particle–particle split is a property of mode geometry, not of channel.
3. **Theorem 4 of the phase-alignment note is unique only in one spatial
   dimension.** In $d \ge 2$ its two conditions leave a $(d-1)$-parameter
   family of vertices; energy is conserved by every member and none dephases,
   so Corollary 4.1 selects nothing. See Theorem C1 and §4.3.

---

## 1. The one algebraic fact

Everything in §§2–3 follows from a single substitution. In $1+1$ D the
four-action algebra is written with $A$, the shift by $q$ cells on the
momentum lattice. In $4$ D the momentum lattice is two-dimensional, cells are
indexed by $\vec n \in \mathbb{Z}^2$, and a Fourier mode of the potential
carries a *vector* index $\vec q$. Define

```math
(A_{\vec q} W)_{\vec n} \;=\; W_{\vec n + \vec q},
```

the shift by the mode's own wavevector. Then the joint QLE stencil per mode is

```math
\dot W_{\vec n} \;=\; \Gamma_{\vec q}(\mathbf{X})\bigl(W_{\vec n + \vec q} - W_{\vec n - \vec q}\bigr),
\qquad
\Gamma_{\vec q}(\mathbf{X}) = -\frac{V_{\vec q}}{\hbar}\sin(\vec k_{\vec q}\negthinspace\cdot\negthinspace\mathbf{X} + \phi_{\vec q}),
```

with the mode index running over

| system | $\vec q$ for the external potential | $\vec q$ for the pair potential |
|---|---|---|
| 2p/1D | $(q, 0)$ and $(0, q)$ | $(q, -q)$ only |
| 1p/2D | all of $\mathbb{Z}^2$ | — |

$\mathbf{X}$ is $(x_1, x_2)$ or $(x, y)$; the momentum cell size is
$\Delta p = \pi\hbar/L$ per axis, so a shift by $\vec q$ is a momentum
displacement of $\hbar\vec k_{\vec q}/2$ — the same half-photon as in
$1+1$ D.

The point of writing it this way is that $A_{\vec q}$ obeys exactly the
algebra $A$ did. Everything the ladder proved about $A$ is a statement about
commuting shift operators and their Fourier symbols, and none of it knew the
lattice was one-dimensional.

---

## 2. The four-action layer

### 2.1 The exactness theorem lifts verbatim

**Theorem A1.** *With $`A \mapsto A_{\vec q}`$, the complete solution family of
the four-rule note's identity (3) is unchanged:*

```math
F = \bigl(A_{\vec q} - A_{\vec q}^{-1}\bigr) G,
\qquad
H = \bigl(2 - A_{\vec q} - A_{\vec q}^{-1}\bigr) G \;-\; \Gamma\thinspace\mathrm{I} \;+\; H_0 .
```

*Proof.* The necessity argument passes to the Fourier symbol on the periodic
momentum lattice. In $1+1$ D that substitution was
$A \mapsto e^{i q\theta}$; here it is
$`A_{\vec q} \mapsto e^{i \vec q \cdot \vec\theta}`$ with
$\vec\theta \in [0,2\pi)^2$. Write $\psi = \vec q\cdot\vec\theta$. The symbol
equation is
$`2(1 - \cos\psi)\hat F + (-2i\sin\psi)\hat H = 2i\Gamma\sin\psi`$,
character-for-character the $1+1$ D equation with $q\theta$ replaced by
$\psi$, and the divisibility argument goes through on the range of $\psi$.
$\blacksquare$

Verified in Part A to a worst deviation of $1.8\times10^{-15}$ over five mode
geometries (external on one particle, external on both, anti-diagonal pair
modes, oblique 1p/2D modes, and a mixed 2p/1D set), for both the symmetric
member (5) and the pure-hop member $G = 0$.

The participant-locality argument of §5 of the sea-dressed note also survives
without change: the focus event touches $\vec n$ and $\vec n \pm \vec q$, the
hop event touches only $\vec n \pm \vec q$, and imposing that no rate
reference a non-participant forces $G = (\Gamma/2)\thinspace\mathrm{I}$ as
before. So the **symmetric member remains the uniquely local one in 4 D**, and
Theorem F3 of `../supplement/four_action_foundations.md` — endpoint locality
as the selection principle — carries over unmodified.

### 2.2 Modes fibre the lattice, and different modes fibre it differently

**Proposition A2.** *Mode $`\vec q`$ decomposes the joint momentum lattice into
disjoint chains — the orbits of $`A_{\vec q}`$ — each a line of direction
$`\vec q`$. The four-action process runs independently on each chain.*

This is the first structurally new feature. In $1+1$ D every mode acts along
the one momentum axis; mode $q$ merely interleaves it into $\gcd(q, M_p)$
sub-chains. In 4 D, chains for $\vec q = (1,0)$ and chains for
$\vec q = (1,-1)$ are transverse families covering the same lattice. There is
no privileged momentum axis, and the "column of three adjacent cells" picture
of the original slide is mode-relative.

Part C tabulates the orbits on a $12\times12$ lattice: $(1,0)$, $(0,1)$,
$(1,-1)$, $(1,2)$ and $(3,-1)$ each give 12 chains of length 12, while
$(2,-2)$ gives 24 of length 6 and $(6,0)$ gives 72 of length 2 — the
familiar $\gcd$ aliasing, now with a direction attached.

### 2.3 What leaks: the mode wavevector decides

Let $u$ be any direction in the joint momentum space and consider the linear
functional $u\cdot\mathbf{P}$.

**Theorem A3 (leak law).** *Per event,*

- *Focus and Defocus change $u\cdot\mathbf{P}$ by $0$, for every $u$ and
  every $\vec q$;*
- *a Right-Hop across $\vec n$ changes it by $2(u\cdot\vec q)\thinspace\Delta p$,
  a Left-Hop by the negative.*

*Consequently $`u\cdot\mathbf{P}`$ is conserved event-by-event iff
$`u \perp \vec q`$ for every active mode; the conserved subspace is the
orthogonal complement of the span of the active wavevectors.*

*Along the directions this theorem excludes, a weaker invariant survives:
if the potential is periodic along $`u`$ with period $`a`$, then
$`u\cdot\mathbf{P}`$ is conserved **modulo** $`\pi\hbar/a`$, on all of
$`\mathbb{R}^d`$ and with no box. Theorem A3 is the continuum part of the
invariant and Theorem O5 of
[`open_position_space.md`](open_position_space.md) §4.2 is its modular
companion; neither implies the other.*

*Proof.* Focus removes one quantum from each of $\vec n \pm \vec q$ and
deposits two at $\vec n$: the displacement sums to
$(\vec n - \vec n - \vec q) + (\vec n - \vec n + \vec q) = 0$ identically.
Hop moves one quantum from $\vec n - \vec q$ to $\vec n + \vec q$, a
displacement of $2\vec q$. $\blacksquare$

This is Noether's theorem on the lattice, and it collapses the 4-D
supplement's conservation column into one line. Part B confirms all eight
entries. Reading off the consequences:

- **Pair potential, $\vec q = (q,-q)$.** $u = (1,1)$ gives $u\cdot\vec q = 0$:
  the total momentum $P = p_1 + p_2$ is untouched by *all four* actions. The
  hop is a full-photon transfer between the two particles; the focus is a
  half-photon exchange. Both are particle–particle. This is correction 2 of
  §0.1.
- **External potential on particle 1, $\vec q = (q,0)$.** $u\cdot\vec q = q$:
  the hop leaks to the field, the focus does not.
- **1p/2D with $q_y = 0$.** $p_y$ is conserved — the potential is
  $y$-independent, so this is ordinary translational Noether, recovered per
  event.

The focus channel's unconditional neutrality deserves emphasis: it is a
*two-body collision in the momentum lattice* whatever the mode geometry, and
it is the only channel that can never transfer momentum to the field. That is
a stronger statement than the $1+1$ D note made, and it is the sharpest
version of Cyganski's "momentum-balanced exchange" reading.

### 2.4 The corrected comparison

![Four-action stencils on the joint momentum lattice](https://raw.githubusercontent.com/billpage/wpmw/output/figures/fourd_microdynamics_mode_geometry.png)

Replacing §6 of the 4-D supplement:

| | pair mode $\vec q = (q,-q)$ | external mode $\vec q = (q,0)$ | oblique mode $\vec q = (q_x,q_y)$ |
|---|---|---|---|
| appears in | 2p/1D | 2p/1D | 1p/2D |
| focus: momentum leak | $0$ | $0$ | $0$ |
| hop: leak along $u$ | $2(u\cdot\vec q)\Delta p$ | $2(u\cdot\vec q)\Delta p$ | $2(u\cdot\vec q)\Delta p$ |
| conserved directions | $u \perp \vec q$ | $u \perp \vec q$ | $u \perp \vec q$ |
| chain direction | anti-diagonal | axis $1$ | $\vec q$ |
| exactness family | Theorem A1 | Theorem A1 | Theorem A1 |

The columns differ only in the *value* of $\vec q$. There is no 2p/1D-versus-1p/2D
distinction at this layer at all. What the supplement recorded as a difference
between two physical set-ups is the difference between a wavevector that lies
on the anti-diagonal and one that does not.

Part C makes the collapse explicit: the generator of 2p/1D with
$V_2(x_1-x_2)$ and the generator of 1p/2D with $V(x,y) = V_2(x-y)$ agree to
$0$ — they are the same array, because they are the same equation.

---

## 3. The sea-dressed layer

### 3.1 The joint sea is a product sea, but the shift is not a product

The crystal shift in 4 D is $W' = W + (2/h)^2$. As a *density* the background
factorises exactly, $(2/h)^2 = (2/h)\cdot(2/h)$: a joint sea object is one
sea positon for each particle, uncorrelated. But the shift does not commute
with taking products.

**Proposition B1.** *For a product state $`W^{(2)} = W_1 W_2`$,*

```math
\bigl(W_1 W_2 + (2/h)^2\bigr) \;-\; \bigl(W_1 + 2/h\bigr)\bigl(W_2 + 2/h\bigr)
\;=\; -\thinspace\frac{2}{h}\thinspace\bigl(W_1 + W_2\bigr).
```

Verified to $1.1\times10^{-16}$ (Part D), with the gap reaching $2\times$ the
background for minimum-uncertainty Gaussians. The missing cross terms are the
configurations "one real excess particle paired with one sea particle". So a
2-particle system is *not* obtained by running two copies of the $1+1$ D
sea-dressed algorithm and multiplying: the joint excess is $W_1 W_2$, and the
one-real-one-sea sector has been absorbed into the joint background. Anyone
implementing the multi-body sea must instantiate it jointly.

### 3.2 The sixteen channels lift, and the sea stops being optional

The channel table of §4 of the sea-dressed note is stated in terms of the
roles $hi = \vec n + \sigma\vec q$ and $lo = \vec n - \sigma\vec q$ and the
sea factors $s_c = S_c/B$ at joint cells $c$. Nothing in the exactness
derivation (§5 there) uses one-dimensionality, so by Theorem A1 the sixteen
channels reproduce the joint QLE stencil exactly at pinned sea.

What changes is *availability*. In the joint lattice a Focus event needs two
participants at the same joint cell — that is, agreement on **both**
particles' positions and both momenta up to $\pm\vec q$.

**Proposition B2.** *With $`\mathcal{W}`$ excess worlds on a joint lattice of
$`(M_x M_p)^{dN}`$ cells, the expected number of excess worlds per cell is
$`\mathcal{W}(M_x M_p)^{-dN}`$.*

At $M_x = M_p = 32$ this is $98$ per cell at $dN = 1$ with
$\mathcal{W} = 10^5$, but $9.5\times10^{-2}$ at $dN = 2$ and
$9.3\times10^{-5}$ at $dN = 3$; even $\mathcal{W} = 10^8$ leaves
$9.3\times10^{-2}$ per cell at $dN = 3$ (Part D). Excess–excess collisions
essentially never happen once $N > 1$.

The sea, by construction, has $B \gg 1$ in every joint cell at every $dN$.
So the conclusion is not merely that the sea is convenient in the multi-body
case: **it is the only collision partner that exists.** The sea-dressed
note's framing — the sea as the reservoir that linearises the rates — gains a
second, independent justification in 4 D that has no $1+1$ D counterpart,
where excess–excess encounters are common.

### 3.3 The cost, stated plainly

Instantiating $B$ per joint cell is not free. Two arithmetic facts bound what
the representation can deliver.

**Proposition B3.** *For a Gaussian state the peak of $`|W|`$ in units of
$`(2/h)^d`$ equals the state's purity.*

Verified in Part D and again in Part G. For a *pure joint* state the peak
therefore saturates the bound at any $N$ — there is no peak signal loss. The
loss is in the marginals: extracting a one-particle observable integrates the
joint excess over $(M_x M_p)^{d(N-1)}$ cells, each carrying its own
$\sqrt{B}$ of shot noise, while the signal stays $O(1)$. In the world-ensemble
representation of `../algorithm/multi_body_extension.md` §5, where the
background is never instantiated, the same exponential reappears as the
ordinary sign problem in the signed excess. **The two representations trade a
sign problem for a shot-noise problem; neither removes it.** This is the
concrete form of §9.3 of the multi-body spec, and it applies to
*distinguishable* particles in a product state — it is not a fermionic
effect.

Proposition B3 also gives the sea-dressed reading of entanglement: for two
harmonically coupled oscillators the reduced purity falls from $1$ to $0.779$
as $\omega_c$ runs $0 \to 3$ (Part G), and that number is *literally* the
excess-to-background ratio of the reduced description. **Entanglement is
excess-to-background loss.**

---

## 4. The phase-alignment layer

### 4.1 Misalignment in joint configuration space

The transported phase generalises with no incident, over the joint
configuration space $\mathbf{X}$:

```math
\Phi_j(\mathbf{X},t) \;=\; \theta_j \;+\; \frac{\mathbf{p}_j\negthinspace\cdot\negthinspace(\mathbf{X} - \mathbf{X}_j) - E_j t}{\hbar},
\qquad E_j = \frac{|\mathbf{p}_j|^2}{2m},
```

and $\mu = \Phi_a - \Phi_b$. Proposition 3 becomes

```math
\nabla_{\negthinspace\mathbf{X}}\thinspace\mu = \frac{\Delta\mathbf{p}}{\hbar},
\qquad
\frac{d\mu}{dt}\bigg|_{\mathbf{v}} = \frac{\Delta\mathbf{p}}{\hbar}\negthinspace\cdot\negthinspace\bigl(\mathbf{v} - \bar{\mathbf{v}}_{\mathrm{pair}}\bigr).
```

Lemma 4 (completeness) and Lemma 5 (locality) are unaffected: they are
statements about gauge invariance and about the pump writing the same $\mu$ at
each place, and neither counts dimensions.

The one qualitative change is that $\mu$ now winds along a *direction*.
$\mu$ is constant on every hyperplane orthogonal to $\Delta\mathbf{p}$. In
$1+1$ D that orthogonal complement is trivial; in 4 D it is three-dimensional
in configuration space and, more consequentially, one-dimensional in the
momentum plane relevant to the vertex.

### 4.2 Theorem 4 loses its uniqueness

**Theorem C1.** *In $`d`$ spatial dimensions the two conditions of Theorem 4 —*

```math
\mathbf{p}_{\mathrm{out}} - \mathbf{p}_{\mathrm{in}} = \mathbf{p}_a - \mathbf{p}_b
\quad (d \text{ equations}),
\qquad
\dot\mu = 0 \iff \Delta\mathbf{p}\negthinspace\cdot\negthinspace\bigl(\mathbf{p}_{\mathrm{in}} + \mathbf{p}_{\mathrm{out}} - \mathbf{p}_a - \mathbf{p}_b\bigr) = 0
\quad (1 \text{ equation})
```

*— have solution set*

```math
\mathbf{p}_{\mathrm{in}} = \mathbf{p}_b + \mathbf{t},
\qquad
\mathbf{p}_{\mathrm{out}} = \mathbf{p}_a + \mathbf{t},
\qquad
\mathbf{t} \perp \Delta\mathbf{p},
```

*a family of dimension $`d - 1`$. The swap is the member $`\mathbf{t} = 0`$.*

*Proof.* Substituting the first condition into the second gives
$`2\thinspace\Delta\mathbf{p}\cdot\mathbf{p}_{\mathrm{in}} = \Delta\mathbf{p}\cdot(\mathbf{p}_a + \mathbf{p}_b - \Delta\mathbf{p}) = 2\thinspace\Delta\mathbf{p}\cdot\mathbf{p}_b`$,
i.e. $`\Delta\mathbf{p}\cdot(\mathbf{p}_{\mathrm{in}} - \mathbf{p}_b) = 0`$.
$\blacksquare$

Part E confirms the rank count: $d = 1$ gives family dimension $0$, $d = 2$
gives $1$, $d = 3$ gives $2$.

**Corollary C2 (energy does not select).** *Every member conserves energy
identically.* With $`\mathbf{p}_{\mathrm{in}} = \mathbf{p}_b + \mathbf{t}`$ and
$`\mathbf{p}_{\mathrm{out}} = \mathbf{p}_a + \mathbf{t}`$,

```math
\bigl(|\mathbf{p}_{\mathrm{out}}|^2 + |\mathbf{p}_b|^2\bigr) - \bigl(|\mathbf{p}_{\mathrm{in}}|^2 + |\mathbf{p}_a|^2\bigr)
= 2\thinspace\mathbf{t}\negthinspace\cdot\negthinspace(\mathbf{p}_a - \mathbf{p}_b) = 0 .
```

Verified to $\le 3.6\times10^{-15}$ across the family, alongside
$|\Delta P| \le 4.4\times10^{-16}$ and $|\dot\mu| \le 4.1\times10^{-16}$ — the
off-swap members do not dephase either, so the $\mathrm{sinc}$ envelope of
§4.2 of the phase-alignment note does not suppress them.

**Proposition C3 (what does select).** *The swap is the unique member for
which the out-multiset of momenta is a permutation of the in-multiset.*
Indeed $`\lbrace\mathbf{p}_b + \mathbf{t}, \mathbf{p}_a, \mathbf{p}_b\rbrace`$
equals $`\lbrace\mathbf{p}_a + \mathbf{t}, \mathbf{p}_b, \mathbf{p}_b\rbrace`$
only if $`\mathbf{t} = 0`$ (or $`\mathbf{p}_a = \mathbf{p}_b`$, the degenerate
aligned case). Part E confirms.

![The exchange vertex in two dimensions](https://raw.githubusercontent.com/billpage/wpmw/output/figures/fourd_microdynamics_vertex_family.png)

So the honest statement is: **in $d \ge 2$ the swap must be postulated, not
derived.** The candidate postulate is worth naming, since it is already
implicit in the ladder —

> **(X) Exchange-only.** A vertex permutes the momenta already present among
> the participating worldlines; no worldline leaves on a momentum absent from
> the in-state.

— and it is exactly the property that Part A of `demo_phase_alignment.py`
observed and §7 of the sea-dressed note relies on when it calls the identity
assignment at a collision a gauge choice. In $1+1$ D (X) was a *consequence*
of Theorem 4; in $d\ge2$ it must be an input. Theorem 4 then reads: momentum
conservation, stationarity of $\mu$, **and (X)** have the unique solution
$`\mathbf{p}_{\mathrm{in}} = \mathbf{p}_b`$, $`\mathbf{p}_{\mathrm{out}} = \mathbf{p}_a`$,
with energy conservation still automatic.

This does not weaken the note's achievement; it relocates it. Theorem 4 was
advertised as deriving energy conservation and the selection rule from one
condition rather than two. Corollary C2 shows energy conservation was even
cheaper than advertised — it follows from momentum conservation plus
stationarity in any dimension — while the *uniqueness* of the out-state, which
was the sharper of the two claims, was carrying a hidden one-dimensional
assumption.

### 4.3 What the transverse parameter is

The free parameter is not exotic. Since $\mathbf{t} \perp \Delta\mathbf{p}$
and, by Corollary 4.3, $\Delta\mathbf{p} \parallel \hbar\vec k_{\vec q}$ for a
pumped pair:

- **2p/1D, pair mode.** $\Delta\mathbf{p} \parallel (1,-1)$, so
  $\mathbf{t} \parallel (1,1)$: the free parameter is a **mismatch of the
  centre-of-mass momentum $P = p_1 + p_2$** between the incoming world and the
  struck sea pair (Part E, orthogonality to $2.2\times10^{-17}$). The swap
  demands that the incoming world match the mate in *both* the relative and
  the centre-of-mass momentum; the general family demands the relative
  component only. Given that §7 of the 4-D supplement shows the COM sector is
  exactly decoupled from the pair dynamics, requiring the swap imposes a
  correlation the dynamics never uses. That is an argument the $1+1$ D case
  could not raise, and it is a reason to take the non-swap members seriously
  rather than dismiss them.
- **1p/2D.** $\mathbf{t}$ is a transverse momentum-row mismatch: a particle
  may be struck by a pair sitting in a different row, provided the rows agree
  along $\vec k_{\vec q}$.

**Proposition C4 (invisible to the QLE at uniform sea).** *Every member
delivers the same transfer $`\Delta\mathbf{p}`$ to the excess particle. Enlarging
the admissible partner set from $`\mathbf{t} = 0`$ to the whole transverse line
multiplies the available partner density by a $`\vec q`$-dependent constant,
absorbed into the vertex calibration, provided the sea is transversally
uniform.* So the family is a further **gauge freedom of the microdynamics**,
in the same sense as the $`G`$-freedom of the four-rule family and the null bias
$`H_0`$ — a recurring pattern in this ladder. It ceases to be invisible under
level-2 (live-sea) dynamics, where depletion is anisotropic: then the choice
of $`\mathbf{t}`$ changes which rows are drained, and the two readings diverge.

### 4.4 The sea acquires a direction label

Corollary 4.3 in vector form reads
$`|\mathbf{p}_{\mathrm{out}} - \mathbf{p}_{\mathrm{in}}| = |\Delta\mathbf{p}|`$
*and* the transfer is along $\Delta\mathbf{p}$. A pair split by mode
$\vec q_1$ has $\Delta\mathbf{p} \parallel \vec k_{\vec q_1}$ and can mediate
transfers along that direction only. In $1+1$ D the sea's rung label was a
magnitude, $2q\thinspace\Delta p$; in 4 D it is a **vector**. Two consequences:

- Multi-mode superposition (open item 3 of the phase-alignment note) becomes
  vectorial: a multi-mode pump splits pairs along several directions and the
  order parameter $Z$ must be resolved by direction, not just by rung.
- In 2p/1D with a pair potential, all pumped pairs are split along
  $(1,-1)$: the sea is **centre-of-mass degenerate**, positon and negaton
  carrying equal $P$ and differing only in relative momentum. That is a
  concrete, checkable structural prediction about the sea in the interacting
  case.

---

## 5. Two particles and a harmonic potential

### 5.1 What the problem is, and what it is not

Three distinct systems go by this name, and they behave differently:

| | Hamiltonian addition | character |
|---|---|---|
| common trap | $\tfrac12 m\omega_0^2(x_1^2 + x_2^2)$ | confining, separable |
| harmonic coupling | $\tfrac12 m\omega_c^2(x_1-x_2)^2$ | confining, entangling |
| repulsive barrier | $-\tfrac12 m\omega^2 x^2$ | genuinely scattering |

A first caution: **a harmonic pair potential is confining, not scattering.**
It grows without bound, so there are no asymptotic free states, no S-matrix,
and "two particles scattered by a harmonic potential" has no scattering
content in the technical sense. What it has is a crossing of classical
trajectories in a trap. If a scattering problem with an exact answer is
wanted, the third row is the one to use — and it turns out to be far more
interesting (§5.6).

### 5.2 The conventional Wigner treatment

The joint Wigner function of the two-particle density matrix is

```math
W^{(2)}(x_1,p_1,x_2,p_2) = \frac{1}{(\pi\hbar)^2}\negthinspace\int\negthinspace dy_1 dy_2\;
\rho(x_1{+}y_1, x_2{+}y_2;\thinspace x_1{-}y_1, x_2{-}y_2)\thinspace
e^{-2i(p_1 y_1 + p_2 y_2)/\hbar},
```

real, bounded by $(2/h)^2$, with the usual marginals. It obeys the Moyal
equation $\partial_t W = \lbrace\negthinspace\lbrace H, W\rbrace\negthinspace\rbrace_{\mathrm{M}}$,
whose expansion in $\hbar$ runs over odd orders $\lambda$ and carries
$\partial^\lambda_{\mathbf{X}} V$.

**The decisive fact.** For a Hamiltonian quadratic in $(\mathbf{X},\mathbf{P})$
all third and higher derivatives vanish, so the Moyal bracket **truncates at
the Poisson bracket**:

```math
\partial_t W^{(2)} \;=\; \lbrace H, W^{(2)}\rbrace_{\mathrm{PB}}
\qquad\text{exactly, to all orders in }\hbar .
```

The QLE *is* the classical Liouville equation. Everything follows:

- **Exact solution by symplectic covariance.** With
  $H = \tfrac12\zeta^{\mathsf T} \mathbb{H}\zeta$ and $\dot\zeta = \mathbb{A}\zeta$,
  $\mathbb{A} = \Omega\mathbb{H}$, the solution is
  $W(\zeta, t) = W_0(e^{-\mathbb{A}t}\zeta)$ — rigid transport of the initial
  Wigner function along classical trajectories, for *any* initial state,
  Gaussian or not.
- **Normal modes.** $X = (x_1{+}x_2)/\sqrt2$ at $\omega_+ = \omega_0$ and
  $u = (x_1{-}x_2)/\sqrt2$ at $\omega_- = \sqrt{\omega_0^2 + 2\omega_c^2}$;
  Part G reproduces both to 6 decimals.
- **The ground state is a classical invariant.** Its covariance satisfies
  $\mathbb{A}\Sigma + \Sigma\mathbb{A}^{\mathsf T} = 0$ to
  $2.2\times10^{-16}$ (Part G) — nothing quantum is needed to hold it still.
- **Entanglement.** The reduced symplectic eigenvalue gives an entropy
  matching the closed form
  $-\ln(1-\xi^2) - \xi^2\ln\xi^2/(1-\xi^2)$,
  $\xi = (\sqrt{\omega_-}-\sqrt{\omega_+})/(\sqrt{\omega_-}+\sqrt{\omega_+})$,
  to 6 decimals across $\omega_c \in [0,3]$. The joint Wigner function stays
  a positive Gaussian throughout: **this is an entangled state with no
  negativity anywhere**, a useful reminder that the two resources are
  independent.

![Two coupled oscillators](https://raw.githubusercontent.com/billpage/wpmw/output/figures/fourd_microdynamics_two_particle_harmonic.png)

### 5.3 The microdynamic representation, layer by layer

Every layer must therefore degenerate to classical drift. Concretely:

- **Four-action.** The net generator must equal $\nabla V\cdot\nabla_{\mathbf{P}}W$
  with no finite-difference remainder. The individual channels do not vanish;
  their *mode sum* cancels down to a drift.
- **Sea-dressed.** The sea must carry polarisation in every mode at once,
  with split-rung populations $\propto |V_{\vec q}|$. The **split-rung
  spectrum of the sea is the potential's Fourier spectrum** — for a harmonic
  potential that is a $1/q^2$ tail across all rungs, not a single rung.
- **Phase-alignment.** Corollary 4.3 says each vertex transfers exactly the
  splitting of the pair it strikes. A harmonic potential therefore needs pairs
  at every rung, and the classicality of the result is a cancellation
  *across rungs*. **It is invisible at any single vertex**: no local
  inspection of the microdynamics reveals that this potential happens to be
  classical. That is a genuinely interesting structural point — classicality
  is a spectral property of the sea, not a property of the interaction.

### 5.4 The obstruction: no exactly harmonic potential exists on the ring

The microdynamics is defined on a periodic box. The periodic image of
$\tfrac12 m\omega^2 x^2$ on $[-L/2, L/2)$ has the Fourier spectrum

```math
V_q = \frac{m\omega^2 L^2}{2\pi^2}\thinspace\frac{(-1)^q}{q^2},
```

requiring all modes, with a seam at $x = \pm L/2$ where the third derivative
is distributional. The exact ring QLE force term approaches the classical
drift only as the momentum quantum becomes fine compared with the state's
momentum width. Part F1, for the harmonic ground state:

| $L$ | $\Delta p/\sigma_p$ | max relative residual |
|---|---|---|
| 4 | 1.111 | $2.2\times10^{-1}$ |
| 8 | 0.555 | $4.2\times10^{-4}$ |
| 16 | 0.278 | $9.1\times10^{-15}$ |
| 32 | 0.139 | $4.0\times10^{-15}$ |

So a harmonic potential *can* be represented essentially exactly — but only by
making $L$ large enough that the state spans many momentum cells and the seam
is far away. Both requirements push in the expensive direction.

### 5.5 The cost: an exact mean bought with unbounded noise

**Theorem D (mode-independent noise).** *For the ring harmonic, the product of
the per-mode rate amplitude and the squared hop displacement is independent of
$`q`$:*

```math
\frac{|V_q|}{\hbar}\thinspace\bigl(2q\thinspace\Delta p\bigr)^2
\;=\; \frac{m\omega^2 L^2}{2\pi^2 q^2\hbar}\cdot\frac{4q^2\pi^2\hbar^2}{L^2}
\;=\; 2\thinspace m\thinspace\omega^2\hbar .
```

Verified to $1.1\times10^{-16}$ for $q = 1\ldots8$ (Part F2). The rate falls as
$1/q^2$ and the squared jump grows as $q^2$; they cancel exactly.
**The injected momentum-space variance therefore grows linearly in the number
of retained modes.** The exact $\tau$-leap variance rate on the ground state
(Part F3):

| $q_{\max}$ | 1 | 2 | 4 | 8 | 16 | 32 |
|---|---|---|---|---|---|---|
| Var rate | $3.19\times10^5$ | $7.98\times10^5$ | $1.81\times10^6$ | $3.81\times10^6$ | $7.72\times10^6$ | $1.50\times10^7$ |

with a mode-independent increment of $\approx 4.9\times10^5$ per mode, against
$1.48\times10^5$ for the single-mode cosine well of the existing demos — the
$q_{\max} = 32$ ring harmonic is $102\times$ noisier.

Seeded runs (Part F4, 8 seeds) confirm it, with a nuance worth recording:

| $q_{\max}$ | rel$L^2(W)$ | sd $\langle p\rangle$ | sd $\langle p^2\rangle$ |
|---|---|---|---|
| 1 | 0.0150 | $3.6\times10^{-4}$ | $7.2\times10^{-4}$ |
| 8 | 0.0223 | $1.06\times10^{-3}$ | $6.4\times10^{-3}$ |
| 32 | 0.0237 | $3.73\times10^{-3}$ | $2.50\times10^{-2}$ |

**The field norm saturates; the momentum moments do not.** Extra modes fire
rarely (rate $\propto 1/q^2$) so they barely perturb the occupancy field, but
each event throws a quantum $2q$ cells, so momentum observables see them
directly. Anyone benchmarking convergence on $\|W\|$ alone will conclude the
scheme has converged when its momentum moments have not.

Two further scalings compound the problem. The total rate budget is
$\sum_q |V_q|/\hbar = m\omega^2 L^2 / 12\hbar$, growing as $L^2$ — and $L$ is
exactly what §5.4 says must grow. And $q_{\max}$ must track $M_p$, which must
also grow to resolve the state.

![The cost of a harmonic potential](https://raw.githubusercontent.com/billpage/wpmw/output/figures/fourd_microdynamics_harmonic_cost.png)

**Verdict.** A harmonic potential is the microdynamics' *worst* case: the
answer is exactly classical, and the model reaches it by cancellation among an
unbounded number of individually noisy channels. It is an excellent
**regression test** — the exact answer is known for any initial state, and any
spurious diffusion, negativity creation or drift shows up immediately — and a
**poor probe**, since it exercises none of the quantum microdynamics. This is
presumably why the existing demos use a single-mode cosine well; the reason is
now quantified.

The natural mitigation is the split already used for short-range Coulomb in
§10 of `../algorithm/multi_body_extension.md`: integrate the classical drift
$-V'(x)\thinspace\partial_p$ deterministically and sample stochastically only
the quantum remainder, which for a harmonic potential is identically zero and
so costs nothing. It should be recorded that this is in tension with the
four-action ontology, in which *all* momentum change is discrete exchange;
the tension is real and the precedent for resolving it pragmatically already
exists in the spec.

### 5.6 Negativity is transported, never created — and an integrator bug found on the way

Part H evolves a cat state through the harmonic ring for a full period. The
first attempt showed negativity growing by $3.4\times$, which is not physics:
the jump generator's momentum-Fourier symbol
$\lambda(x,\theta) = 2\sum_q \Gamma_q(x)\sin(q\theta)$ is purely imaginary, so
the explicit Euler substep of `step_jump_fourier` amplifies the norm by
$\sqrt{1+\lambda^2 dt^2}$ per step, i.e. by $\exp(\lambda^2\thinspace dt\thinspace T/2)$
over a time $T$. Part H1 confirms the law directly (predicted $3.57$, $1.37$,
$1.08$ against observed $1.41$, $1.27$, $1.24$ at three step sizes; the gap at
the coarsest step is the seam error isolated below).

Because all modes are functions of the same shift operator they commute, and
the substep is **exactly integrable**: $\exp(dt\thinspace\mathcal{L})$ is
multiplication by $e^{i\lambda dt}$ in the momentum-Fourier basis. With that
substitution (Part H2):

| $L$ | grid | $q_{\max}$ | steps | negativity ratio | rel$L^2$ vs exact rotation |
|---|---|---|---|---|---|
| 8 | $64\times64$ | 24 | 800 | 1.1357 | $7.2\times10^{-2}$ |
| 16 | $128\times128$ | 48 | 200 | 1.00555 | $7.8\times10^{-4}$ |
| 16 | $128\times128$ | 48 | 800 | 1.000207 | $3.7\times10^{-5}$ |
| 24 | $192\times192$ | 72 | 800 | 1.000231 | $3.7\times10^{-5}$ |

At $L = 8$ the cat's tail sits on the seam and the $`7\%`$ error is the ring
artefact of §5.4. At $L = 16$ the harmonic ring rotates the state rigidly and
negativity is preserved to $2\times10^{-4}$ — all of it Strang splitting error,
converging as the step is refined. **A harmonic potential neither creates nor
destroys Wigner negativity.**

The exact-substep integrator is a free improvement to the deterministic mesh
solver and is recommended for promotion to
`wpmwlib.phase_space_crystal_lattice` (open item 5).

### 5.7 The inverted barrier: the test case that should replace the trap

For $V = -\tfrac12 m\omega^2 x^2$ the Moyal bracket *also* truncates, so the
QLE is again exactly classical — yet a wavepacket incident on the barrier
tunnels, and the transmitted probability
$\int_{x>0}\negthinspace dx\int\negthinspace dp\thinspace W(x,p,t)$ is exact,
because the Wigner evolution is exact and the position projector is a
legitimate phase-space observable. The tunnelling is carried entirely by
classical transport of the initial Wigner function, including its negative
regions.

This makes the inverted barrier a much sharper test than the trap: the
microdynamics must produce *zero* quantum content and still reproduce a
quantum tunnelling probability. Any scheme that manufactures negativity, or
that damps the initial negativity through jump noise, fails visibly on an
observable with a closed-form reference. It is the recommended replacement for
"two particles scattered by a harmonic potential" as the simple test case.

---

## 6. Suggested test cases, ranked

Ordered by ratio of information gained to work required. Each entry states
what it exercises and what the reference answer is.

1. **Cosine pair potential**, $V_2(r) = V_p\cos(2\pi r/L)$, two particles on
   a ring, no external field. *Exercises:* correlated jumps, the leak law
   (Theorem A3) with $P$ exactly conserved, the sixteen channels at a joint
   sea. *Reference:* the COM/relative decomposition of §7 of the 4-D
   supplement — the joint mesh evolution must agree marginal-by-marginal with
   the existing $1+1$ D cosine-well code run at reduced mass $\mu = m/2$.
   *Why first:* single mode, so the stencil is exact and none of §5's mode-sum
   pathology appears. This is the direct two-particle analogue of
   `demo_cosine_well_microdynamics.py`.
2. **Inverted parabolic barrier** (§5.7). *Exercises:* whether the scheme
   transports negativity faithfully under a Moyal-truncating potential.
   *Reference:* exact, by classical transport of the initial Wigner function;
   Kemble/parabolic-cylinder transmission for the stationary comparison.
3. **Entanglement generation from a product state.** Two Gaussian wavepackets
   colliding through a pair potential; observable is the purity of the reduced
   one-particle state, equivalently the peak of the reduced Wigner function in
   units of $2/h$ (Proposition B3). *Exercises:* the correlated-jump rule
   specifically. *Why it matters:* a mean-field (TDHF-Wigner) implementation —
   item 3 of §12 of the multi-body spec — produces **identically zero**
   entanglement, so this is the observable that separates a genuine joint
   solver from a mean-field one. No other test in this list does that.
4. **Quartic pair coupling**, $V_2(r) = \lambda r^4$. *Exercises:* the first
   genuinely non-classical two-particle dynamics — $\partial_r^3 V \ne 0$, so
   the Moyal series does not truncate. *Reference:* small-basis exact
   diagonalisation on the joint grid.
5. **Non-separable 1p/2D cross mode**, $V(x,y) = V_1(x) + V_2(y) + V_{\times}\cos(k_x x + k_y y)$.
   *Exercises:* vectorial jumps and the mode-dependent fibration
   (Proposition A2); with $V_\times = 0$ it reduces to two independent
   $1+1$ D problems, giving a built-in control. *Free companion:* the
   2p/1D $\leftrightarrow$ 1p/2D regression of Part C — run the pair-potential
   problem through the 1p/2D code path and require bit-level agreement.
6. **Contact interaction**, $V_2(r) = g\thinspace\delta(r)$ (Lieb–Liniger at
   $N = 2$). *Exercises:* mode truncation at its worst — a flat Fourier
   spectrum, so Theorem D's noise sum diverges even faster than for the
   harmonic case. *Reference:* Bethe-ansatz ground state.
7. **2D isotropic oscillator with $L_z \ne 0$.** A state with a negative
   Wigner core, rigidly rotating. *Exercises:* vectorial streaming plus
   negativity transport, with no Moyal content at all. *Reference:* exact
   rotation. The 1p/2D analogue of test 2.
8. **Two identical fermions in 1D.** *Exercises:* whether the antisymmetric
   joint Wigner function's exchange structure survives the microdynamics, and
   how Proposition B3's shot-noise penalty behaves when $|W^{(2)}|$ saturates
   its bound over the whole dynamically relevant region. Expect this to be
   hard; it is the controlled entry point to §9.3 of the multi-body spec.

**A limitation to record.** All of the above are scalar-potential problems,
which is all the model handles. A magnetic field enters the Wigner equation
through a gauge-dependent kernel with no Fourier-mode jump representation, so
1p/2D Landau-level and Aharonov–Bohm problems — the most natural physics
questions to ask of a particle in a plane — are outside the current formalism.
Extending to vector potentials is a real open problem, not a matter of
bookkeeping.

---

## 7. Numerical verification

All outputs are from `src/demo_fourd_microdynamics.py`, seed `20260809`,
container run August 2026.

- **Part A (Theorem A1).** Symmetric member (5) and pure-hop member $G = 0$
  against the joint QLE stencil, five mode geometries $\times$ four random
  occupancy fields on a $12\times12$ momentum lattice: worst deviation
  $1.8\times10^{-15}$.
- **Part B (Theorem A3).** Focus and Defocus leak $0$ in all eight
  mode/direction combinations tested; hops leak exactly
  $\pm 2(u\cdot\vec q)\thinspace\Delta p$.
- **Part C (Proposition A2).** Orbit counts on a $12\times12$ lattice for
  seven wavevectors. The 2p/1D and 1p/2D generators for the same pair
  potential agree to $0$ (identical arrays).
- **Part D (Proposition B1, B2, B3).** Product identity residual
  $1.1\times10^{-16}$, gap up to $2\times$ the background; peak/$(2/h)$ equals
  purity to 6 decimals for three Gaussians; co-location arithmetic tabulated
  for $dN = 1\ldots4$.
- **Part E (Theorem C1, Corollary C2, Proposition C3).** Family dimension
  $0, 1, 2$ at $d = 1, 2, 3$; across four family members, residual
  $\le 4.4\times10^{-16}$, $|\Delta E| \le 3.6\times10^{-15}$,
  $|\Delta P| \le 4.4\times10^{-16}$, $|\dot\mu| \le 4.1\times10^{-16}$;
  permutation holds only at $\mathbf{t} = 0$. Pair-mode splitting is
  orthogonal to the COM direction to $2.2\times10^{-17}$.
- **Part F (Theorem D).** Ring residual table of §5.4; per-mode constant
  $2m\omega^2\hbar$ to $1.1\times10^{-16}$ for $q = 1\ldots8$; variance-rate
  and seeded-moment tables of §5.5.
- **Part G.** Normal-mode frequencies to 6 decimals; classical stationarity of
  the joint ground state to $2.2\times10^{-16}$; entanglement entropy against
  the closed form to 6 decimals for six couplings.
- **Part H.** Euler-amplification law and the exact-substep convergence table
  of §5.6.

---

## 8. Open items

1. **Postulate (X).** Decide whether exchange-only is adopted as an axiom of
   the phase-alignment layer or whether the transverse family
   (Theorem C1) is physical. §4.3 gives an argument on each side; the
   discriminating test is a live-sea run in 4 D, where Proposition C4's
   invisibility fails and the two readings predict different anisotropic
   depletion patterns.
2. **Directional sea ledger.** Recast the order parameter $Z_r$ of
   `relational_pairing_and_carrier_lock.md` to resolve pair splitting by
   *direction* as well as rung (§4.4). This subsumes open item 3 of the
   phase-alignment note and open item 3 of the sea-dressed note, both of which
   assumed a scalar rung label.
3. **Joint sea instantiation.** Proposition B1 says a multi-body sea cannot be
   assembled from per-particle seas. Specify the joint construction, and
   quantify Proposition B3's marginalisation shot noise against the world
   ensemble's sign problem on a common benchmark — test case 3 of §6 is the
   natural vehicle.
4. **Classical-drift/quantum-remainder splitting.** Implement the split of
   §5.5 and measure the variance reduction on smooth potentials. Record
   explicitly how it sits with the four-action ontology; the honest answer may
   be that the ontology describes the exact generator while the implementation
   is free to integrate any exactly-cancelling part deterministically.
5. **Promote the exact jump substep.** Add
   `step_jump_fourier_exact` to `wpmwlib.phase_space_crystal_lattice`
   (§5.6). The existing explicit Euler is unconditionally amplifying for the
   deterministic mesh path; the exact form costs one FFT pair per step and
   removes a systematic that is easy to mistake for physics.
6. **Vector potentials.** The magnetic case (§6, limitation) is genuinely
   outside the formalism. Whether the four-action decomposition has any
   analogue for the gauge-covariant Wigner kernel is unexplored.

---

## 9. Sources

- [`four_rule_microdynamics_equivalence.md`](four_rule_microdynamics_equivalence.md)
  — the exactness family lifted in §2.
- [`sea_dressed_microdynamics.md`](sea_dressed_microdynamics.md) — the
  sixteen channels and the participant-locality argument, lifted in §3.
- [`phase_alignment_microdynamics.md`](phase_alignment_microdynamics.md) —
  Theorem 4, corrected in §4.
- [`../supplement/phase_space_crystal_lattice_4d_supplement.md`](../supplement/phase_space_crystal_lattice_4d_supplement.md)
  — the kinematics this note supplies dynamics for; its §6 table is
  superseded by §2.4 here.
- [`../supplement/four_action_foundations.md`](../supplement/four_action_foundations.md)
  — Theorem F3 (endpoint locality), which §2.1 shows is dimension-independent.
- [`../algorithm/multi_body_extension.md`](../algorithm/multi_body_extension.md)
  — the $N$-body specification; §9.3 (sign problem) and §10 (Ewald split) are
  connected to Propositions B2–B3 and §5.5 respectively.
- David Cyganski, *Extended Fokker–Planck Eq. and the QLE V2*; *Wigner
  Collisions Diagram* (Sozi deck); *A journey from Bohm trajectory theory,
  through Nelson's SDEs and Wigner Particles to the Closed Four Action Model*
  (3 August 2026).
- Generating script: `src/demo_fourd_microdynamics.py`. Figures published to
  the `output` branch as `figures/fourd_microdynamics_*.png`.
