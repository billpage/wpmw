# Foundations and Uniqueness of the Four-Action Wigner Particle Model

**A reading of D. Cyganski's note *A journey from Bohm trajectory theory, through Nelson's SDEs and Wigner Particles to the Closed Four Action Model* (3 August 2026) against the project's existing derivation ladder, focused on §"Does momentum and energy balance fully determine the FAWPM?".**

---

## 0. Status and provenance

The source is David Cyganski's seven-page note of 3 August 2026 (this project's
`Bohm_to_Nelson_and_four_action_Wigner_August_3_2026.pdf`). Its final two pages
pose the question this supplement is about: could the four-action rates have
been *derived* from symmetry and conservation rather than *matched* to the
Wigner stencil? David reports a five-step answer obtained with ChatGPT on
4 August 2026 and flags it explicitly as **not yet assessed**. Assessing it is
the job here.

The prior project documents involved are
[`docs/analysis/four_rule_microdynamics_equivalence.md`](../analysis/four_rule_microdynamics_equivalence.md)
(hereafter **FR**), which proved exactness and derived the complete solution
family, and
[`docs/analysis/sea_dressed_microdynamics.md`](../analysis/sea_dressed_microdynamics.md)
(hereafter **SD**), which realizes the rates as two-body collisions against a
pinned sea.

Companion code: `src/demo_four_action_foundations.py`. Every numerical claim in
§11 is an output of that script.

**Summary of findings.** The five-step chain reaches the right answer, and one
of its steps is stronger than advertised (§7). But its billing — *symmetry plus
conservation determines the model* — is not correct as stated. Momentum and
energy balance constrain exactly one of the two rate coefficients and are
provably blind to the other (§5). The coefficient they cannot see is precisely
the one carrying the first quantum correction. What actually selects the model
is David's own original design principle — that no rate may be weighted by the
intermediate-state density — and that principle does the whole job by itself
(§6). Finally, the residual quantum input is isolated to a single number: the
momentum grid step must equal half the photon momentum of the field mode (§8).
§2.1 was added after review, to state correctly what distinguishes the
world-particle ensemble from the many-interacting-worlds models.

---

## 1. Page 4 resolves FR open item 1

FR was written from a screenshot of one slide with the second bias rate cut
off. It reconstructed the missing rate from the exactness theorem and recorded
the prediction as open item 1. Page 4 of the new note supplies the slide in
full. The comparison:

| Quantity | FR §4(b), predicted | Note p. 4, as written | Agree? |
| --- | --- | --- | --- |
| Focus reaction | $(n-q,\thinspace n+q) \leftrightarrow (n,\thinspace n)$ | $(k-1,\thinspace k+1) \leftrightarrow (k,\thinspace k)$ | yes |
| Focus rate | $f_n = \tfrac{\Gamma}{2}(W_{n+q} - W_{n-q})$ | $r_k = \tfrac{\kappa}{2}(W_{k+1} - W_{k-1})$ | yes |
| Hop reaction | $n-q \to n+q$ | $k-1 \leftrightarrow k+1$ | yes |
| Hop rate | $h_n = -\tfrac{\Gamma}{2}(W_{n+q} + W_{n-q})$ | $\lambda_k = -\tfrac{\kappa}{2}(W_{k-1} + W_{k+1})$ | yes |
| Indexing | reaction **centre** | reaction centre | yes |
| Rate field | $\Gamma_q(x) = -\tfrac{V_q}{\hbar}\sin\negthinspace\left(\tfrac{2\pi q x}{L} + \phi_q\right)$ | $\kappa = -\tfrac{V_p}{\hbar}\sin\negthinspace\left(\tfrac{2\pi n x}{L}\right)$ | yes, $\phi_q = 0$ |

The reconstruction was correct on every point, including the indexing
convention that FR §5 had to guess. **FR open item 1 is closed.** The
alternative slide-indexed member FR §4(c) is not what David intends and can be
dropped. FR §5's cautionary note — that a $W$-independent hop bias would
generate nothing — is answered: $\lambda_k$ is occupancy-dependent as required.

Two dictionary entries for reading the note against the repository: David's
$\kappa$ is FR's $\Gamma_q(x)$, and David's $n$ inside the sine is the mode
index, which the repository writes $q$ (the repository's $n$ is the momentum
cell index, which the note writes $k$). Beyond this section the repository
convention is used.

---

## 2. The Bohm-to-Nelson chain, and the quantum potential as a projection artifact

Pages 1–3 are a compressed route from Bohm's trajectory equations, through the
osmotic velocity $u = \tfrac{\hbar}{2m}\nabla\ln\rho$ and Nelson's SDEs, to the
observation that the quantum-potential acceleration

```math
a(x,t) \;=\; (u\cdot\nabla)u \;+\; \frac{\hbar}{2m}\nabla^2 u
```

has no counterpart in the QLE. David's interpretation (p. 3) is that the
quantum potential is **the artifact of collapsing a phase-space description
onto its position and momentum projections**, with the analogy that gas
streamlines obey complicated nonlinear equations while the underlying particles
obey Newton plus collisions.

This interpretation is correct, and it is not merely an analogy — it is a
theorem, and worth recording as such because it is the cleanest available
justification for the whole world-particle programme. Taking the first two
$p$-moments of the QLE gives a continuity equation and a momentum-flux
equation, in which the flux splits as

```math
\int \frac{p^2}{m}\thinspace W\thinspace dp \;=\; m\rho v^2 \;+\; \Pi ,
```

and for a pure state the residual stress $\Pi$ satisfies
$\partial_x \Pi = \rho\thinspace\partial_x Q$ identically, with $Q$ Bohm's quantum
potential. The quantum potential *is* the closure term of the moment hierarchy,
in exactly the way the pressure tensor is the closure term of the
Chapman–Enskog hierarchy for a gas. This is Takabayasi's 1954 result — his
abstract states the aim as exhibiting the correspondence between the
phase-space formulation and the quantum-potential formulation — and a modern
treatment is Wyatt, *Quantum Dynamics with Trajectories* (Springer 2005),
ch. 15. Both are in `references/bibliography.md`.

Two riders are worth attaching to the claim:

- The identification $\partial_x\Pi = \rho\thinspace\partial_x Q$ holds **for pure
  states**. For a mixed state, $\Pi$ contains an additional classical-like
  dispersion and $Q$ is no longer recoverable from it. The projection reading
  therefore explains why $Q$ is nonlinear in $\rho$ (it is a ratio of moments)
  *and* why $Q$ is a pure-state notion.
- The analogy has real predictive content: it explains David's own observation
  on p. 3 that in the world-particle models *"any dependence on local density
  derivatives has disappeared"*. Density derivatives are moment-closure
  artifacts, and a model that never takes the moments never generates them.
  **But the antecedent has to be stated carefully** — §2.1 does that, because
  the obvious reading of it is false.

### 2.1 Many-Interacting-Worlds is the test case, and it fails the naive reading

*(This subsection answers an objection raised by Bill Page: the
many-interacting-worlds models of Hall–Deckert–Wiseman and of Poirier look like
counterexamples to the bullet above.)*

They are counterexamples to the naive reading, and the naive reading is that
the relevant contrast is **continuum versus discrete** — that a hydrodynamic
description carries density derivatives because it is a continuum, and a
particle description escapes them because it is granular. That reading is
wrong, and MIW is exactly the case that shows it.

Hall, Deckert and Wiseman replace the single Bohmian trajectory with $N$
interacting world-particles carrying independent positions and momenta. It is
granular by construction. Yet the density derivatives are not gone; they have
been re-expressed as a finite difference over the ensemble. Their interworld
potential is

```math
U_N(X) \;=\; \frac{\hbar^2}{8m}\sum_n\left(\frac{1}{x_{n+1}-x_n} - \frac{1}{x_n - x_{n-1}}\right)^{\negthinspace 2},
```

whose bracket is a nearest-neighbour estimate of
$`\partial_x \rho / \rho`$ obtained from the ansatz
$`\rho(x_n) \approx [N(x_n - x_{n-1})]^{-1}`$ — the density read off from
inter-world spacings. HDW say so themselves: the quantity being squared
approximates **Nelson's osmotic momentum**, the same $u$ that appears on p. 1
of David's note, and their appendix A shows the resulting force converges to
$`-\nabla Q`$. The three-body potential yields a *five*-body force. Poirier and
Schiff–Poirier are the same structure with the granularity removed: a continuum
of trajectories labelled by a parameter, with the quantum force built from
derivatives with respect to that label, which is the density in Lagrangian
coordinates. HDW describe their own approach as a discretisation of
Holland–Poirier, and note that Schiff and Poirier have *"a continuum of
trajectories (i.e. flow lines), not a discrete set of worlds"*.

**The axis that actually matters is what the ensemble samples.** MIW worlds
sample $\rho(x)$: their initial velocities are set from the single-valued
Bohmian field $`p = \partial_x S`$, and the ordering $x_1 < \dots < x_N$ is
preserved, so at any instant the ensemble is a set of points on a **graph over
configuration space** — a Lagrangian submanifold of phase space, sampled $N$
times. A WPMW ensemble samples $W(x,p)$: many world-particles sit at the same
$x$ with different $p$, filling a two-dimensional region. Discretising a
streamline description does not convert it into a kinetic one, in the same way
that tracking $N$ tagged fluid parcels does not turn hydrodynamics into
Boltzmann.

The sharpest test is the oscillator ground state, whose wavefunction is real,
so $S$ is constant and **every MIW world momentum is exactly zero**. The
ensemble reproduces $\mathrm{Var}(x)$ but has $\mathrm{Var}(p) = 0$, against
the quantum value $\hbar m\omega/2$. HDW have to recover the momentum spread
from a separately constructed *nonclassical momentum*

```math
p^{\rm nc}_n \;=\; \frac{\hbar}{2}\left(\frac{1}{x_{n+1}-x_n} - \frac{1}{x_n - x_{n-1}}\right),
```

which is the osmotic momentum again, and their uncertainty relation is between
$\Delta x$ and $`\Delta p^{\rm nc}`$ — not between $\Delta x$ and $\Delta p$.
The information a world-particle carries as a **coordinate**, an MIW world has
to **reconstruct from its neighbours' spacings**. That is the whole
distinction, and §11 Part G verifies it numerically.

![MIW worlds on a Lagrangian graph versus a WPMW sample of the Wigner function](https://raw.githubusercontent.com/billpage/wpmw/output/figures/miw_vs_wpmw_ensembles.png)

*Left: the exact 41-world MIW oscillator ground state, on the line
$`p = 0`$, over contours of the Wigner function it is supposed to represent.
Centre: a WPMW ensemble sampling the same state over phase space. Right:
momentum marginals — the MIW worlds contribute a spike at the origin, and the
quantum spread is recovered only through the reconstructed
$`p^{\rm nc}`$.*

Two corollaries, each an independent fingerprint of the same distinction and
each worth knowing when the models are compared:

- **Where the nonlocality goes.** MIW pays for eliminating the wavefunction
  with nonlocality in *position*: a world must consult its neighbours in
  configuration space, and HDW describe the result as a
  *"super-nonlocality"*. The four actions are strictly local in $x$ and pay in
  *momentum* instead — finite jumps of $\hbar k_q/2$, never a gradient. The
  trade is not incidental; it is the same trade in a different coordinate.
  It is also why §7's postulate is about a momentum grid step: that is where
  this model's nonlocality is kept.
- **Sign.** An MIW ensemble is positive by construction, because it samples a
  probability density. WPMW requires positons *and* negatons because $W$ goes
  negative. An ensemble that must carry signed weight cannot be a sampling of
  $\rho$, so the negaton is not an awkward extra: it is the marker that the
  ensemble is kinetic rather than hydrodynamic.

The amended claim, then, is not that granularity removes density derivatives.
It is that **an ensemble whose members carry independent momenta — one that
samples $W(x,p)$ rather than $\rho(x)$ — never generates them**. MIW and
Poirier satisfy the negation of that antecedent and duly generate them, which
makes them confirmations rather than counterexamples.

---

**The quantitative companion.** §2.1 settles the kinematics but says nothing
about *cost*, which is what the sparse-ontology literature is actually about.
[`representation_cost_and_annihilation.md`](representation_cost_and_annihilation.md)
supplies that half: the sampling cost of any world-particle representation is
$`N/\|\mu\|_1^2`$, a cat state costs a bounded factor $`2.679`$ however large
the cat, interference consumes momentum resolution rather than position
density, entanglement is free and only non-Gaussianity is charged for, and the
interference lobe — which carries $`38.9\%`$ of the ensemble at zero
probability density — is precisely the region MIW leaves empty. It also prices
the two representations this note does not use: the pair $`(x,\mu)`$ ensemble
costs a factor that diverges as the position lattice is refined, and the
positon/negaton sea costs $`1 + 2A/h`$, which for the specification's grid is
exactly the number of momentum cells. Its §9 answers Hackebill and Poirier
directly, and its §7 measures the annihilation burden that the answer leaves
behind.


## 3. What conservation buys

Write the endpoint-local ansatz with two free gains, which is the note's step 3:

```math
f_n \;=\; a\thinspace\Gamma\thinspace\bigl(W_{n+q} - W_{n-q}\bigr),
\qquad
h_n \;=\; b\thinspace\Gamma\thinspace\bigl(W_{n+q} + W_{n-q}\bigr).
```

Here $f_n$ is the signed net rate of focus events centred at $n$ (David's
$r_k$) and $h_n$ the signed net rate of right-hops across $n$ (David's
$\lambda_k$). Let $\delta$ be the momentum grid step, so a hop moves one
quantum by $2q\delta$ and a focus moves two quanta by $q\delta$ in opposite
directions.

### 3.1 Momentum

Focus and defocus conserve total particle momentum **event by event and for any
rate whatsoever** — that is structural, not a constraint. Hops change it by one
photon $2q\delta = \hbar k_q$. Summing over centres,

```math
\dot P \;=\; 2q\delta\sum_n h_n \;=\; 4q\delta\thinspace b\thinspace\Gamma\thinspace\rho(x),
\qquad \rho(x) = \sum_n W_n .
```

Ehrenfest requires $\dot P = -\bigl(\partial_x U\bigr)\rho$. With
$U = V_q\cos(k_q x)$, $\Gamma = -\tfrac{V_q}{\hbar}\sin(k_q x)$ and
$\delta = \hbar k_q/(2q)$ this gives $-2 k_q V_q b\sin(k_q x) = k_q V_q\sin(k_q x)$,
hence

```math
b \;=\; -\tfrac{1}{2}, \qquad\text{for every } a .
```

### 3.2 Energy

A single focus event costs kinetic energy $q^2\delta^2/m$, independent of the
centre $n$; a single defocus event returns it. This is David's "thorn" (p. 5).
The resolution he reports — treat the actions as one system — is exactly the
statement that

```math
\sum_n f_n \;=\; a\thinspace\Gamma\sum_n \bigl(W_{n+q}-W_{n-q}\bigr) \;=\; 0
```

by telescoping, so the focus channel is energy-neutral **in aggregate, for any
$a$**. The hop channel then contributes
$\sum_n h_n\thinspace(2q\delta)\thinspace p_n/m = -(\partial_x U)\thinspace P(x)/m$ at
$b = -\tfrac{1}{2}$, which is force times current: exactly the power that the
free-streaming step removes from the potential energy. Total energy is conserved.

The balance is exact at every position column and every instant in the
mean-field sense, which is stronger than an average over time. Two caveats:

1. **It is aggregate, not microscopic.** No individual focus or defocus event
   conserves energy; the cancellation is between events at different centres,
   population-weighted. The slide's phrase *"microscopically conservative
   interactions"* is exactly true for particle number and, in the focus
   channel, for momentum — but not for energy. This deserves flagging because
   it is the same distinction that separates FR's exactness-in-expectation from
   pathwise equivalence (FR §7).
2. **The mediating mode must be energyless for this to be the only option.** The
   phase-alignment work established that exactness requires a relaxational
   ($\omega = 0$) mediating mode; a dispersive mode contaminates the result. A
   zero-frequency mode supplies momentum $\hbar k_q$ at zero energy cost, so the
   hop channel *cannot* balance energy event-by-event either. The aggregate
   balance is not a convenience — it is forced by the same condition that forces
   exactness. This is a satisfying consistency, and it is new: it says the
   energy bookkeeping and the overdamped-mode requirement are two faces of one
   constraint.

### 3.3 Conservation is implied by exactness

The stronger statement, which subsumes both subsections, follows immediately
from FR's theorem. Every exact member of the family has

```math
F = \bigl(A - A^{-1}\bigr) G, \qquad
H = \bigl(2 - A - A^{-1}\bigr) G - \Gamma\thinspace\mathrm{I} + H_0 ,
```

and both $(A - A^{-1})$ and $(2 - A - A^{-1})$ annihilate the constant
functional $`W \mapsto \sum_n (\cdot)_n`$. Hence $\sum_n f_n = 0$ and
$\sum_n h_n = -\Gamma\rho$ for **every** $G$, which is precisely §3.1 and §3.2.

> **Proposition F1.** Every member of the FR family — including the original
> single mediated-jump rule ($G = 0$) — conserves particle number, momentum and
> energy exactly. Conservation is a consequence of QLE-exactness, not an
> independent property of the four-action member.

Verified in §11 Part A to $2.5\times10^{-14}$ over random $G$ of half-range 0, 1
and 2. Consequence: conservation cannot distinguish members of the family, and
so cannot be what selects the four-action model from within it.

---

## 4. Conservation cannot determine the model

The previous section shows conservation is necessary but says nothing about
sufficiency. Sufficiency fails, and it is worth seeing exactly how.

Assembling the generator (FR eq. 2) from the two-parameter ansatz gives, with
$\Delta_j W_n := W_{n+jq} - W_{n-jq}$,

```math
\dot W_n \;=\; 2a\thinspace\Gamma\thinspace\Delta_1 W_n \;-\; (a+b)\thinspace\Gamma\thinspace\Delta_2 W_n .
\qquad\text{(F1)}
```

Matching the QLE stencil $\Gamma\Delta_1 W_n$ requires $2a = 1$ **and**
$a + b = 0$. Together with $b = -\tfrac{1}{2}$ from §3.1 that is three
conditions on two unknowns: the system is **overdetermined and consistent**.
That consistency is a genuine result and is the real content of the note's
"satisfying outcome" — it is not automatic and would fail for a mis-specified
action set.

But now note what (F1) says about the conserving line $b = -\tfrac{1}{2}$:

```math
\dot W_n \;=\; 2a\thinspace\Gamma\thinspace\Delta_1 W_n \;+\; \bigl(\tfrac{1}{2}-a\bigr)\thinspace\Gamma\thinspace\Delta_2 W_n .
```

The first moment of $c\thinspace\Delta_j$ is $-2jc\rho$, so the Ehrenfest force is
proportional to $2a\cdot 1 + (\tfrac{1}{2}-a)\cdot 2 = 1$ — **independent of
$a$**. The whole one-parameter family is Ehrenfest-exact.

> **Counterexample.** Take $a = 0$, $b = -\tfrac{1}{2}$: focus channel switched
> off entirely, hops at half the QLE rate over twice the distance. This model
> conserves particle number, momentum and energy exactly, reproduces the
> classical force law exactly, is translation- and reflection-covariant, and is
> **not** the QLE. Its generator is $\tfrac{\Gamma}{2}\Delta_2$, which by §8 is
> Moyal evolution at $\hbar_{\mathrm{eff}} = 2\hbar$.

So a physicist who knew only the four actions, the symmetries and the
conservation laws would be free to build a universe with the wrong Planck
constant. Verified in §11 Part B.

### 4.1 Why conservation must fail: it is a moment condition

The reason is structural rather than accidental. Particle number, momentum and
energy conservation are conditions on the moments $M_0$, $M_1$, $M_2$ of the
generator. Computing directly from (F1):

| moment | coefficient | depends on |
| --- | --- | --- |
| $M_0$ | $0$ | — |
| $M_1$ | $4b\thinspace\Gamma\rho$ | $b$ only |
| $M_2$ | $8b\thinspace\Gamma\sum_n n W_n$ | $b$ only |
| $M_3$ | $12b\thinspace\Gamma\sum_n n^2 W_n \thinspace+\thinspace (12a + 16b)\thinspace\Gamma\rho$ | $a$ **and** $b$ |

> **Proposition F2.** In the endpoint-local family, the focus gain $a$ is
> invisible to moments $0$, $1$ and $2$ and first appears in moment $3$. Its
> contribution there, $12a\thinspace\Gamma\rho$, is smaller than the
> $a$-independent part by a factor of order
> $q^2\delta^2/\langle p^2\rangle = (\hbar k_q/2)^2/\langle p^2\rangle$ — that
> is, exactly the Moyal expansion parameter.

**The focus gain multiplies the first quantum correction, and conservation laws
are classical-order statements.** No finite set of conservation laws could ever
have fixed $a$; the quantum content of the QLE lives in the third and higher
momentum derivatives, which the conserved moments do not reach. This, and not
any defect in the argument, is why step 4 of the note's chain cannot get past
$b$. Verified in §11 Part C.

![Conservation constrains a line; only closure picks the point](https://raw.githubusercontent.com/billpage/wpmw/output/figures/four_action_uniqueness_map.png)

*Left: QLE residual over the plane of rate gains, with the conservation line
$`b=-1/2`$ (cyan) and the closure line $`a+b=0`$ (green). Right: a cut along the
conservation line — the momentum-and-energy residual sits at machine epsilon
for every focus gain, while the QLE residual has a single sharp zero.*

---

## 5. What does determine the model: endpoint locality

If conservation does not select the four-action model, something else must.
It turns out to be the very thing David introduced it for.

Page 4 of the note gives the motivation in the boxed remark: *"the intermediate
state density is no longer weighting the actions"*, offered as an aesthetic
improvement over the lattice model's rates, which consulted the occupancy of
the cell lying between the two states a jump connects. That remark is not
decoration. It is the selection principle.

> **Theorem F3 (endpoint locality selects the symmetric member).** Within the
> FR family of exactly QLE-reproducing schemes, there is exactly one member
> whose two rate laws are supported on the endpoint cells $\lbrace n-q,\thinspace n+q \rbrace$
> alone. It is $G = \tfrac{\Gamma}{2}\mathrm{I}$, i.e. the four-action model.
> No null-bias freedom survives.

*Proof.* Pass to the symbol on the momentum lattice, $A \mapsto z$. From FR's
theorem, $\hat F(z) = (z - z^{-1})\hat G(z)$. Requiring $F$ supported on
$\lbrace +q, -q \rbrace$ means $\hat F = \alpha(z - z^{-1})$ for a scalar $\alpha$,
so $\hat G = \alpha$ and $G = \alpha\thinspace\mathrm{I}$. Then
$\hat H = \alpha(2 - z - z^{-1}) - \Gamma + \hat H_0$, whose coefficient at
$z^0$ — the centre cell — is $2\alpha - \Gamma$ plus the $H_0$ contribution.
Elements of $\ker(A^{-1}-A)$ on a finite periodic lattice are supported at the
symbols $z = \pm 1$, hence are lattice-wide averages and are not endpoint-local
unless zero; so $H_0 = 0$ and $2\alpha - \Gamma = 0$. Therefore
$\alpha = \Gamma/2$, giving
$f_n = \tfrac{\Gamma}{2}(W_{n+q}-W_{n-q})$ and
$h_n = -\tfrac{\Gamma}{2}(W_{n+q}+W_{n-q})$. $\blacksquare$

Three things follow.

- The original single mediated-jump rule is the $G = 0$ member, whose hop rate
  is $-\Gamma W_n$ — supported on the **centre cell only**. It is the *maximally*
  non-endpoint-local member. David's dissatisfaction with it and his arrival at
  the four-action model are therefore two ends of one theorem, not an
  aesthetic preference followed by a lucky find.
- Endpoint locality also kills the null-bias gauge freedom that FR left open,
  which conservation and symmetry do not touch.
- It is a *single* principle doing the work of steps 3–5 of the note's chain.
  It is also physically motivated in its own right: an event's rate should
  depend on the states it connects, not on a state it does not visit. That is
  ordinary locality in momentum space.

Verified in §11 Part E.

---

## 6. Audit of the five-step chain

With the above, each step of the note's p. 6 summary can be graded.

| Step | Claim | Verdict |
| --- | --- | --- |
| 1 | Translation symmetry gives the same local law at every momentum centre | **Sound.** This is exactly FR's translation-invariance hypothesis, and it is needed. |
| 2 | Reflection symmetry plus force reversal make $r$ endpoint-antisymmetric and $\lambda$ endpoint-symmetric | **Sound, with a rider.** Under $p\to-p$ combined with $x\to-x$ (so $\Gamma\to-\Gamma$), focus maps to focus and right-hop to left-hop, which gives exactly the stated parities. The rider: the argument presupposes the rates are proportional to $\Gamma$, which is linear response (step 3), not symmetry. |
| 3 | Endpoint locality and linear response give $r = a\thinspace g(W_+-W_-)$, $\lambda = b\thinspace g(W_-+W_+)$ | **This is a postulate, not a consequence** — and by §5 it is the decisive one. It should be promoted to the top of the argument and labelled as the model's defining principle rather than buried as a technical restriction. |
| 4 | Momentum and energy balance fix $b = -\tfrac{1}{2}$ | **Right answer, over-attributed.** Momentum balance alone fixes $b$ (§3.1); energy balance is then automatic and adds no information (§3.2, §3.3). And by §4 this step cannot reach $a$ at all. |
| 5 | Single-harmonic first-order closure fixes $a + b = 0$ | **Right answer, and stronger than stated.** See below. |

Step 5 deserves the upgrade. ChatGPT's phrasing invokes *first-order* closure —
a perturbative statement, suggesting the result holds only to $O(V_q)$. It does
not need to be perturbative. The spurious channel that $a + b \ne 0$ generates
is $\Gamma\thinspace\Delta_2$: transport of $4q\delta = 2\hbar k_q$ per event, i.e.
**two photons of the mode**. Rather than appealing to perturbation order, state
the selection rule directly:

> **Single-harmonic selection rule.** Rates are exactly linear in the field
> amplitude, and a field mode of wavenumber $k_q$ delivers momentum in units of
> $\hbar k_q$. A mode-$`q`$ potential therefore drives exactly one stencil pair,
> at $\pm q$ cells. Any $\Delta_{2}$ term in the generator would be a
> two-photon process driven by a one-photon amplitude.

This is exact, not first-order, and it is the microscopic reason the QLE
collision term is exactly linear in $U$ at all orders. As a matter of
presentation, the "no spurious 4G channel" language should be replaced by it.

**Overall verdict.** The chain reaches the correct rates and its premises are
individually defensible. But it is not *"conservation and symmetry determine
the model"*. It is *"symmetry, exact linear response, endpoint locality,
Ehrenfest, and a one-photon selection rule determine the model"* — and of those
five, endpoint locality alone suffices once exactness is granted (§5), while
conservation alone provably does not suffice under any circumstances (§4). The
honest headline is stronger and simpler than the advertised one: **the four
actions plus momentum-space locality give the QLE.**

---

## 7. Where the quantum actually enters

David's stated goal (p. 6) is *"a first principle derivation of QM without any
reference to Schrödinger's equation (even that part embedded in the Wigner
distribution definition and QLE)"*. It is worth being precise about how much of
that is achieved, because the answer is both less and more interesting than a
yes or no.

Let the momentum grid step $\delta$ be a **free parameter** rather than
$\pi\hbar/L$. Repeat §3.1 with $\delta$ free: Ehrenfest fixes
$2q\delta\thinspace\Gamma = \partial_x U$, still $b = -\tfrac{1}{2}$, and endpoint
locality still gives $a = \tfrac{1}{2}$. The mean-field generator is

```math
\dot W(p) \;=\; \Gamma\thinspace\bigl[\thinspace W(p + q\delta) - W(p - q\delta)\thinspace\bigr].
```

Taylor-expanding the shift and comparing term by term with the Moyal series
for $U = V_q\cos(k_q x)$:

> **Proposition F4.** The four-action model with free grid step is *exactly*
> Moyal evolution with an effective Planck constant
> ```math
> \hbar_{\mathrm{eff}} \;=\; \frac{2 q\delta}{k_q} .
> ```
> It reduces to the classical Liouville equation as $\delta \to 0$ and equals
> the true QLE **if and only if** $q\delta = \hbar k_q / 2$ — that is, if and
> only if the hop transfers exactly one photon of the mode.

Verified in §11 Part D to machine precision across $\delta$ spanning a factor
of 40.

![The grid step is the whole quantum input](https://raw.githubusercontent.com/billpage/wpmw/output/figures/four_action_hbar_effective.png)

*Left: deviation of the four-action generator from Moyal evolution at
$`\hbar_{\rm eff}`$ (machine epsilon everywhere), from the true QLE (zero only at
$`\delta=\hbar k/2`$), and from classical Liouville (zero only as $`\delta\to 0`$).
Right: relative deviation of the generator's moments from the $`a=1/2`$ member —
flat at machine epsilon through moment 2, departing at moment 3.*

The consequences for the foundational claim are these.

- **What is derived.** Given the phase-space lattice kinematics, the *dynamics*
  is forced. No Schrödinger equation, no wavefunction, no operator algebra, and
  no appeal to the Wigner stencil is needed to get the QLE. That is a real
  result and the note is right to be pleased with it. It is structurally
  analogous to deriving the Schrödinger equation from the canonical commutation
  relations plus Galilean covariance — the dynamics from the kinematics.
- **What is assumed.** The quantum input has not been eliminated; it has been
  *concentrated*. It now sits in a single postulate: the momentum grid step is
  half the photon momentum of the field mode. Everything else about the model —
  the four actions, both rate laws, the conservation laws — follows from
  symmetry, locality and Ehrenfest, and is $\hbar$-blind.
- **Why this is the right place for it to sit.** The postulate is de Broglie's
  relation and nothing more: a spatial mode of wavenumber $k$ carries momentum
  $\hbar k$. It is directly measurable, has been measured, and is exactly the
  input David argues on p. 5 we are entitled to assume — *"the probabilistic
  and quantized behavior of static potentials has been amply demonstrated by
  physical experiments and so can be taken as fundamental."* The uniqueness
  analysis shows that this single concession is not merely sufficient but
  **necessary**: it is the only place $\hbar$ can enter, so nothing further
  needs conceding.

A note on multi-mode consistency, which is a nontrivial check rather than a
formality. Mode $q$ requires $q\delta = \hbar k_q/2 = \pi q\hbar/L$, hence
$\delta = \pi\hbar/L$ **for every $`q`$ simultaneously**. One lattice serves all
modes, with mode $q$ hopping $2q$ cells. Had this failed, the model would
support only monochromatic potentials.

---

## 8. Why these four actions?

The note's uniqueness question is conditioned on *"took as a given the four
actions"*. The conditioning can be partly removed, which strengthens the
foundation, and doing so also clarifies a claim on p. 5.

Consider events at fixed $x$, conserving particle number, confined to the triple
$\lbrace n-q,\thinspace n,\thinspace n+q \rbrace$, and exchanging either zero photons
(hence conserving total particle momentum) or one photon $\hbar k_q$ with the
mode. Then:

- **One-body, zero-photon:** no momentum change is possible, so the event is
  null.
- **One-body, one-photon:** the particle must move by $2q\delta$, giving
  $n - q \to n + q$ and its reverse — **Right-Hop and Left-Hop**.
- **Two-body, zero-photon:** both participants must change, by $\pm q\delta$ in
  opposition, giving $(n-q,\thinspace n+q) \leftrightarrow (n,\thinspace n)$ —
  **Focus and Defocus**.
- **Two-body, one-photon:** e.g. $(n-q,\thinspace n) \to (n,\thinspace n+q)$. This is
  a genuinely distinct microscopic story, but its effect on the occupancy field
  is $W_{n-q}\negthinspace-\negthinspace 1$, $W_{n+q}\negthinspace+\negthinspace 1$, $W_n$ unchanged — **identical to
  Right-Hop**. It is not a fifth channel.
- **Anything transferring $q\delta$ to a single particle** is half a photon and
  is excluded.

So the four actions are the complete list, given the photon quantum and the
triple. That last bullet is where the quantization enters the *action set*, as
opposed to the *rate law* — the same postulate as §7, doing double duty.

The fourth bullet has a further use. It shows the hop channel admits a two-body
reading at no cost, which bears on David's p. 5 point (2) about four-wave
mixing (see §9.1).

---

## 9. Objections and open items

### 9.1 The four-wave-mixing identification needs the sea, not just a large ocean

David's p. 5 point (2) identifies focus/defocus with four-wave mixing, notes
that FWM requires a rate model nonlinear in particle density while the model's
rates are linear, and proposes that a large background with small excess
recovers the linear model at first order.

The first-order argument as stated does not work, and the reason is
instructive. Take a general two-body net rate for the focus channel,
$f_n = c^+ W_{n-q}W_{n+q} - c^- W_n^2$, and linearize about a uniform sea
$W = B + \tilde W$. The derivatives with respect to the two endpoint occupancies
are equal — the linearization is **endpoint-symmetric** for any $c^\pm$ and any
$B$. But the focus rate is endpoint-**anti**symmetric. Large ocean and small
waves cannot convert one into the other. (Confirmed numerically in §11 Part F;
the antisymmetric part is identically zero.)

Two remarks, of which the second is the substantive one.

- The endpoint-symmetric combination that a linearized product *does* produce is
  $\propto (W_{n-q} + W_{n+q})$ — the form of the **hop** rate $\lambda$. Taken
  with §8's fourth bullet, this suggests that if either channel deserves the
  mass-action reading it is the hop, not the focus. The momentum bookkeeping
  points the other way, so the two arguments are in tension and neither should
  be leaned on alone.
- **The project already has the correct construction**, and it is better than
  either. SD realizes the collision term as sixteen local two-body channels
  against a pinned Dirac sea whose mean field is *exactly* the symmetric member.
  The ingredient the p. 5 argument is missing is not a bigger ocean but
  **broken detailed balance**: SD §6 makes each channel run one way, with the
  direction set by the local polarization sign of a *pumped* medium. Antisymmetry
  in the rate then comes from the polarization sign structure, not from
  linearizing a symmetric product. Microreversible mass action would force
  quadratic reverse rates and resurrect FR's no-go lemma; a pumped medium with
  one-way stimulated kinetics is ordinary laser physics.

The recommendation is therefore to keep the FWM analogy for the *vertex
topology* (two in, two out, momentum-conserving) and cite SD rather than the
large-background argument for the *rate law*.

### 9.2 Uniqueness holds for the generator, not for the process

Everything above fixes the mean-field generator. It does not fix the stochastic
process. The signed net rates decompose as $f = f^+ - f^-$ with the common part
$f^+ + f^-$ free, so there remains a one-function family of processes with
identical means and different noise. The slide's *"in a zero potential the
inverse pairs exactly cancel"* pins the common part only at $\Gamma = 0$. This
is FR §7 and FR open item 2 (the variance-optimal member), and it is untouched
by the present analysis. Any claim of uniqueness should be stated as uniqueness
of the generator.

### 9.3 Multi-mode cross-channels

The selection rule of §6 is stated per mode. With several harmonics present,
one must additionally rule out channels driven by a *product* of two mode
amplitudes — a $q_1$ photon and a $q_2$ photon in one event. Exact linearity in
the field amplitude excludes them by fiat, and the exact QLE agrees (its
collision term is linear in $U$), but the microscopic justification for that
linearity in a model where the modes are physical excitations has not been
given. This is the same open item as the algorithm spec's multi-mode question
and should be tracked with it.

### 9.4 The energy ledger is a mean-field statement

Per §3.2, energy conservation holds in expectation, exactly, at every column and
instant, but no individual event conserves it. In a finite-$`\nu`$ Monte Carlo run
the instantaneous energy fluctuates with the Poisson counts. Quantifying that
fluctuation — is it $O(\nu^{-1/2})$ with zero mean, and does it stay bounded
over many periods? — is a cheap and worthwhile numerical check that has not
been run.

---

## 10. Numerical verification

`src/demo_four_action_foundations.py`, six parts. Lattice: 128 momentum cells,
$L = 8$, $V_p = 1.5$, $\hbar = m = 1$, mode $q = 1$.

**Part A — exactness implies conservation (Proposition F1).** Six random
translation-invariant $G$ of half-range 0, 1, 2. Worst deviation across
generator identity, particle number, momentum residual and energy residual:
$2.5\times10^{-14}$.

**Part B — the $(a,b)$ plane.** The generator identity (F1) reproduced to
$1.8\times10^{-15}$. At $b = -\tfrac{1}{2}$ the momentum and energy residuals
are $\le 3\times10^{-14}$ for $a \in \lbrace 0,\thinspace \tfrac{1}{2},\thinspace 2,\thinspace -1 \rbrace$
while the QLE residual is $0$, $1.04$, $3.11$, $3.11$ respectively. At
$a = \tfrac{1}{2}$, $b \ne -\tfrac{1}{2}$ the conservation residuals are $1.00$.
The counterexample of §4 is confirmed: $a = 0$ conserves everything and is not
the QLE.

**Part C — moments (Proposition F2).** Over $a \in \lbrace 0,\thinspace 0.25,\thinspace 0.5,\thinspace 1,\thinspace 2 \rbrace$
at $b = -\tfrac{1}{2}$, deviation from the $a = \tfrac{1}{2}$ member normalized
by each moment's own absolute mass: $4\times10^{-17}$, $7\times10^{-17}$,
$2\times10^{-16}$, $8\times10^{-4}$, $1\times10^{-3}$ for moments $0$–$4$. The
jump between moments 2 and 3 is thirteen orders of magnitude.

**Part D — free grid step (Proposition F4).** Gaussian test state, exact
Hermite derivatives, Moyal series to 24 terms. Deviation from Moyal at
$\hbar_{\mathrm{eff}} = 2\delta/k$ is $\le 10^{-14}$ for
$\delta/(\hbar k/2) \in [0.05,\thinspace 2]$; deviation from the true QLE falls to
$4.3\times10^{-16}$ at $\delta = \hbar k/2$ and is $\ge 4\times10^{-2}$
elsewhere; deviation from classical Liouville falls to $1.5\times10^{-4}$ as
$\delta \to 0.05$.

**Part E — endpoint locality (Theorem F3).** Stencil supports for
$G \in \lbrace 0,\thinspace \tfrac{\Gamma}{2},\thinspace \tfrac{\Gamma}{3},\thinspace \tfrac{\Gamma}{2}+0.4A,\thinspace 0.3A^{-1} \rbrace$.
Only $G = \tfrac{\Gamma}{2}$ leaves both $F$ and $H$ supported on
$\lbrace n-1,\thinspace n+1 \rbrace$; $G = 0$ puts $H$ entirely on the centre cell.

**Part F — mass-action linearization (§9.1).** With $B = 10^3$, $c^+ = 0.9$,
$c^- = 1.4$: symmetric part $9.0\times10^{2}$, antisymmetric part exactly $0$.

**Part G — MIW as a test case (§2.1).** The exact MIW oscillator ground state is
obtained by solving HDW's recurrence
$`\xi_{n+1} = \xi_n - (\xi_1 + \dots + \xi_n)^{-1}`$, shooting on $\xi_1$ for
antisymmetry. Reproduces HDW's analytic $N = 3$ and $N = 4$ configurations to
$1.1\times10^{-15}$; their two constraints $\sum_n \xi_n = 0$ and
$`\sum_n \xi_n^2 = N-1`$ then hold to machine precision without being imposed.
For $N \in \lbrace 11,\thinspace 41,\thinspace 161,\thinspace 641 \rbrace$ (units
$\hbar = m = \omega = 1$):

| $N$ | $\mathrm{Var}(x)$ | exact $\tfrac{N-1}{N}\tfrac{\hbar}{2m\omega}$ | $\mathrm{Var}(p)$, worlds | $`\mathrm{Var}(p^{\rm nc})`$ | quantum $\mathrm{Var}(p)$ |
| --- | --- | --- | --- | --- | --- |
| 11 | 0.454545 | 0.454545 | **0** | 0.454545 | 0.5 |
| 41 | 0.487805 | 0.487805 | **0** | 0.487805 | 0.5 |
| 161 | 0.496894 | 0.496894 | **0** | 0.496894 | 0.5 |
| 641 | 0.499220 | 0.499220 | **0** | 0.499220 | 0.5 |

The world momenta contribute nothing at any $N$. The reconstructed
$`p^{\rm nc}`$ agrees with the exact osmotic momentum
$`(\hbar/2)\,\partial_x \ln\rho`$ on the interior worlds to
$3.0\times10^{-12}$, and its variance equals $\mathrm{Var}(x)$ identically —
so $`\Delta x\thinspace\Delta p^{\rm nc} = (1 - 1/N)\hbar/2`$, saturating HDW's
bound. A WPMW sample of $W(x,p)$ at $N = 4\times10^{5}$ returns
$\mathrm{Var}(x) = 0.5007$ and $\mathrm{Var}(p) = 0.4988$, both correct.

---

## 11. Summary

1. The reconstruction in FR §4(b) matches page 4 of the note in every detail,
   including the indexing convention. **FR open item 1 is closed**, and FR
   §4(c) can be dropped.
2. The quantum potential as a moment-closure artifact (pages 1–3) is correct
   and is a theorem, not an analogy — with the rider that it is a pure-state
   identification.
3. Momentum and energy balance fix $b = -\tfrac{1}{2}$. Momentum alone does it;
   energy is implied. **Neither can reach $a$**, and Proposition F2 explains why:
   $a$ multiplies the first quantum correction, and conservation laws are
   classical-order moment conditions. A fully conserving, Ehrenfest-exact,
   symmetric four-action model exists with the wrong Planck constant.
4. Conservation is a **consequence** of QLE-exactness for every member of the
   FR family (Proposition F1), so it could not have been a selection principle.
5. The actual selection principle is **endpoint locality** — David's own boxed
   remark on page 4 — and by Theorem F3 it selects the four-action model
   uniquely, killing the null-bias freedom as well. It should be promoted from
   a technical restriction to the model's defining postulate.
6. Step 5 of the chain is stronger than its "first-order" phrasing suggests and
   should be restated as an exact one-photon selection rule.
7. The QLE is derived from the lattice kinematics without reference to
   Schrödinger. The residual quantum input is exactly one postulate: the grid
   step is half the photon momentum of the mode (Proposition F4). Below that
   the model is classical Liouville; away from it, Moyal at the wrong $\hbar$.
   The same postulate also closes the action set at four (§8).
8. The moment-closure reading survives the MIW/Poirier objection, but only in
   its corrected form (§2.1): granularity is not what removes density
   derivatives — carrying independent momenta is. MIW discretises the
   streamlines, so its ensemble lies on a Lagrangian graph, its world momenta
   vanish identically for a real ground state, and the density derivatives
   reappear as an osmotic momentum reconstructed from neighbour spacings.
9. The *quantitative* form of the same objection — Hackebill and Poirier's
   sparse ontology — is answered separately in
   [`representation_cost_and_annihilation.md`](representation_cost_and_annihilation.md):
   it does not transfer, because the four-action rates are functions of the
   external potential alone (Theorem F3's endpoint locality), so thinning the
   ensemble degrades the estimate and never the law. What does transfer is a
   threshold in the *unraveling* — the density at which pairwise annihilation
   can still hold the pathwise $`L^1`$ growth $`4V_q/\pi\hbar`$ — and that
   density is currently unmeasured.
10. Open: the FWM rate law needs SD's broken detailed balance, not a large
   background (§9.1); uniqueness is of the generator only (§9.2); multi-mode
   cross-channels (§9.3); the fluctuation size of the energy ledger (§9.4).
   Also open, and inherited from item 9: the immunity argument holds for rates
   that read only the potential, and has *not* been checked for any rate that
   reads a collective order parameter.
