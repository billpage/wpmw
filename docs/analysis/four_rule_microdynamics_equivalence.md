# Four-Rule Two-Body Microdynamics: Analysis and Equivalence with the Single Mediated-Jump Rule

**Analysis of D. Cyganski's proposal (Zoom presentation, 2026) to replace the single mediated momentum-jump rule of the phase-space crystal lattice with four two-body interaction rules, with an exactness proof, an algorithm specification, and numerical verification.**

---

## 0. Status and provenance

The source for this note is a screenshot of one slide (page 3 of 4) presented by David Cyganski in a Zoom meeting. The slide states the four rules, asserts that in zero potential the inverse pairs cancel exactly, and begins to give "two, related, signed bias rates", of which only the first is visible before the page cuts off:

$$r_k = \frac{\kappa}{2}\left(W_{k+1} - W_{k-1}\right)$$

attached to the pair reaction $(k,\ k+2) \leftrightarrow (k+1,\ k+1)$.

This note therefore **reconstructs** the complete rate assignment from first principles rather than transcribing it. The reconstruction is not guesswork: §4 derives the *complete family* of rate assignments that reproduce the QLE stencil exactly, proves the family is exhaustive within its natural class, and shows that the visible slide rate selects a member of it.

**Update (August 2026).** The full slide arrived as page 4 of Cyganski's note *A journey from Bohm trajectory theory, through Nelson's SDEs and Wigner Particles to the Closed Four Action Model* (3 August 2026). The reconstruction is confirmed on every point, including the centre-indexing convention §5 had to guess: the second bias rate is $\lambda_k = -\tfrac{\kappa}{2}(W_{k-1} + W_{k+1})$, exactly the symmetric member (5). Open item 1 is closed and the slide-indexed member §4(c) is not what was intended. See [`../supplement/four_action_foundations.md`](../supplement/four_action_foundations.md) §1, which also audits the note's claim that momentum and energy balance determine the rates.

Companion code: `src/demo_four_rule_equivalence.py`, plus the new library methods `step_jump_four_rule` and `step_jump_four_rule_mc` in `src/wpmwlib/phase_space_crystal_lattice.py`.

---

## 1. The proposal as stated

Quoting the slide (occupancies $W_k$ at three adjacent momentum cells, one position column):

| Rule | Action on $(W_{k-1},\ W_k,\ W_{k+1})$ |
| --- | --- |
| Focus | $(W_{k-1} - 1,\ W_k + 2,\ W_{k+1} - 1)$ |
| Defocus | $(W_{k-1} + 1,\ W_k - 2,\ W_{k+1} + 1)$ |
| Right-Hop | $(W_{k-1} - 1,\ W_k,\ W_{k+1} + 1)$ |
| Left-Hop | $(W_{k-1} + 1,\ W_k,\ W_{k+1} - 1)$ |

with the accompanying text: *"In a stochastic model these actions are continually taking place and in a zero potential the inverse pairs exactly cancel each other out. An applied potential biases the rates causing an imbalance, generating a trend."*

Focus/Defocus are mutual inverses, as are Right-Hop/Left-Hop, so the four rules reduce to **two signed channels**: a focus channel (signed net rate $f$; negative $f$ = defocus) and a hop channel (signed net rate $h$; negative $h$ = left-hop). All four rules conserve total occupancy $\sum_k W_k$ event-by-event — the "microscopically conservative" property claimed on the slide.

---

## 2. Notation and the target generator

We work per position column $x_m$ and per Fourier mode $q$ of the potential, exactly as in `docs/algorithm/phase_space_crystal_lattice_algorithm.md` (the *spec* below). Momentum cells are indexed by $n$ with spacing $\Delta p = \pi\hbar/L$; mode $q$ couples the triple $(n-q,\ n,\ n+q)$; the slide's nearest-neighbor triple is the $q = 1$ case. The QLE-consistent local rate is

$$\Gamma_q(x) = -\frac{V_q}{\hbar}\sin\left(\frac{2\pi q x}{L} + \phi_q\right)$$

and the target — the exact single-cosine QLE stencil of spec §3c, per mode, at fixed $x$ (writing $\Gamma$ for $\Gamma_q(x)$) — is

```math
\dot W_n \;=\; \Gamma\,\bigl(W_{n+q} - W_{n-q}\bigr). \qquad \text{(1)}
```

Recall the spec's **single rule** (§3b): a positon at cell $n$ mediates, at rate $|\Gamma|$ per mediator, the transfer of one positon from $n+q$ to $n-q$ when $\Gamma > 0$ (reverse when $\Gamma < 0$). The mediator is unchanged.

---

## 3. Mean-field generator of a focus/hop scheme

Let $f_n$ be the signed net rate (events per unit time, in occupancy units) of **focus events centered at $n$**: one unit removed from each of $n \pm q$, two units added at $n$. Let $h_n$ be the signed net rate of **right-hops across $n$**: one unit moved from $n-q$ to $n+q$. Summing the channel effects at cell $n$:

```math
\dot W_n \;=\; 2 f_n \;-\; f_{n-q} \;-\; f_{n+q} \;+\; h_{n-q} \;-\; h_{n+q}. \qquad \text{(2)}
```

The terms read: cell $n$ gains $2f_n$ as a focus center; it is drained at rate $f_{n-q}$ as the upper satellite of center $n-q$ and at rate $f_{n+q}$ as the lower satellite of center $n+q$; hops across $n-q$ deliver into $n$ and hops across $n+q$ remove from $n$.

Two immediate observations:

- If $f$ and $h$ are **independent of the occupancies** (functions of $x$ only), then (2) telescopes to zero: $\dot W_n = 2f - 2f + h - h = 0$. Constant-rate versions of the four rules generate *nothing* at mean field. The bias rates **must** depend on $W$ — which is exactly what the visible slide rate does.
- The QLE right-hand side (1) is **linear** in $W$. So we look for $f$ and $h$ linear in $W$.

---

## 4. Exactness theorem

Write $A$ for the shift operator by $q$ cells on sequences, $A W$ having components $W_{n+q}$, and let $F, H$ denote the linear, translation-invariant (per fixed $x$), finite-range maps taking the occupancy sequence $W$ to the rate sequences $f, h$. Equation (2) equals target (1) for **all** $W$ iff, as an operator identity,

```math
\bigl(2 - A - A^{-1}\bigr)\,F \;+\; \bigl(A^{-1} - A\bigr)\,H \;=\; \Gamma\,\bigl(A - A^{-1}\bigr). \qquad \text{(3)}
```

**Theorem (complete solution family).** The solutions of (3) are exactly

```math
F = \bigl(A - A^{-1}\bigr)\,G, \qquad
H = \bigl(2 - A - A^{-1}\bigr)\,G \;-\; \Gamma\,\mathrm{I} \;+\; H_0, \qquad \text{(4)}
```

where $G$ is an arbitrary translation-invariant finite-range operator and $H_0$ is any **null bias** — any element of $\ker(A^{-1} - A)$, i.e. a hop-rate contribution constant across $n$ (mod $2q$-periodic combinations), which moves occupancy in a divergence-free pattern and has no mean-field effect.

*Proof.* Sufficiency: substitute (4) into (3); the $G$ terms cancel because $(2 - A - A^{-1})(A - A^{-1}) = (A^{-1} - A)(2 - A - A^{-1})$ (both are polynomials in the commuting shifts), leaving $(A^{-1} - A)(-\Gamma\thinspace\mathrm{I}) = \Gamma(A - A^{-1})$. Necessity: pass to the Fourier symbol on the periodic momentum lattice, $A \mapsto e^{iq\theta}$. Equation (3) becomes $2(1 - \cos q\theta)\hat F + (-2i\sin q\theta)\hat H = 2i\Gamma\sin q\theta$. At angles where $\sin q\theta = 0$ but $\cos q\theta \ne 1$ this forces $\hat F = 0$; hence the trigonometric polynomial $\hat F$ is divisible by $\sin q\theta$, i.e. $F = (A - A^{-1})G$ for some finite-range $G$. Given $F$, the symbol equation determines $\hat H$ wherever $\sin q\theta \ne 0$, i.e. $H$ is fixed up to $\ker(A^{-1} - A)$. $\blacksquare$

Every member of family (4) reproduces the QLE stencil (1) **exactly** — not perturbatively — for arbitrary occupancy fields. Three members deserve names:

**(a) Pure-hop member** ($G = 0$):

$$f_n = 0, \qquad h_n = -\Gamma\thinspace W_n .$$

Substituting into (2): $\dot W_n = -\Gamma W_{n-q} + \Gamma W_{n+q}$, which is (1). **This is the original single mediated-jump rule.** The spec's rule — mediator at $n$, transfer $n+q \to n-q$ at per-mediator rate $\Gamma$ — is precisely a *left-hop across $`n`$ at net rate proportional to the occupancy of the cell being hopped over*. The mediator count in the spec is the full shifted population $W' = W + 2/h$; the background part is a null bias $H_0$ (it cancels in $h_{n-q} - h_{n+q}$), so spec and pure-hop member have identical mean-field generators. The original "single rule" is therefore not an alternative to the four-rule scheme: it is the $G = 0$ member of the same family, with the focus channel switched off.

**(b) Symmetric member** ($G = \tfrac{\Gamma}{2}\thinspace\mathrm{I}$):

```math
f_n = \frac{\Gamma}{2}\bigl(W_{n+q} - W_{n-q}\bigr),
\qquad
h_n = -\frac{\Gamma}{2}\bigl(W_{n+q} + W_{n-q}\bigr). \qquad \text{(5)}
```

Two related signed bias rates with a common prefactor $\Gamma/2$: the focus channel is biased by the **difference** of the satellite occupancies, the hop channel by their **sum**. This is the natural reading of the slide's *"two, related, signed bias rates"*, and its focus rate has exactly the visible slide form $\tfrac{\kappa}{2}(W_{k+1} - W_{k-1})$ with $\kappa = \Gamma_q(x)$ and $k$ read as the reaction center (see §5).

**(c) Slide-indexed member** ($G = -\tfrac{\kappa}{2}A^{-1}$): if instead the slide's $r_k$ is indexed by the lower participant of the reaction $(k, k+2) \leftrightarrow (k+1, k+1)$, so that the focus rate at center $n$ is $\tfrac{\kappa}{2}(W_n - W_{n-2q})$, this too lies in the family, with the forced hop companion $h_n = \kappa W_{n-q} - \tfrac{\kappa}{2}(W_n + W_{n-2q}) - (\Gamma - \kappa) W_n$ for $\kappa = \Gamma$. It is exact but less symmetric.

All three members were verified to reproduce (1) to machine precision on random occupancy fields for $q = 1, 2, 3$ (worst deviation $1.1 \times 10^{-16}$; Part A of the demo, §9).

---

## 5. Matching the slide

The slide's pair-reaction notation $(k,\ k+2) \leftrightarrow (k+1,\ k+1)$ is the Focus/Defocus channel written as a **momentum-conserving two-particle collision**: two positons in the same position column with momentum indices $k$ and $k+2$ convert to two positons at $k+1$, and back. Total momentum is conserved exactly; each participant changes by $\pm\Delta p$, i.e. the co-located pair exchanges **one half-photon** $q\pi\hbar/L$ (the "sum of two half contributions" of the sozi deck reappears as a genuinely two-body exchange).

The visible rate $r_k = \tfrac{\kappa}{2}(W_{k+1} - W_{k-1})$ references cell $k-1$, which lies **outside** the reaction triple $\lbrace k, k+1, k+2 \rbrace$ as labelled. Under the reading that $r_k$ is indexed by the reaction **center** — i.e. the rate of $(k-1, k+1) \leftrightarrow (k, k)$ is $\tfrac{\kappa}{2}(W_{k+1} - W_{k-1})$ — it coincides exactly with the focus rate of the symmetric member (5). We adopt that reading; if instead the slide's index is the lower participant, member (c) applies and the conclusion (exactness, and hence equivalence) is unchanged. What page 4 must supply is the **second** bias rate; the theorem says it is forced (up to null biases) once the focus rate is fixed, and for reading (b) it is $-\tfrac{\kappa}{2}(W_{k+1} + W_{k-1})$ on the hop channel.

One caution should be flagged for the meeting follow-up: with a $W$-independent hop bias — for example rates biased directly by the potential values or slopes at the site alone — the scheme generates nothing (§3, first observation). Whatever page 4 says, the hop channel's bias must be occupancy-dependent, or the focus channel must carry a compensating non-symmetric term.

---

## 6. What "two-body" does and does not buy

### 6.1 Momentum bookkeeping — a genuine improvement

The four-rule scheme splits the collision term into channels with distinct momentum bookkeeping:

- **Focus/Defocus** conserve total *particle* momentum exactly, event-by-event. They are bona fide two-body collisions between co-located world-particles exchanging one half-photon $q\pi\hbar/L$.
- **Right/Left-Hop** change one particle's momentum by $2q\Delta p = qh/L$ — exactly one photon of the mode $q$ potential component. They are one-body events exchanging a quantum with the potential field, not particle–particle interactions.

The original single rule attributes *every* event's momentum change $\pm qh/L$ to the field. The four-rule scheme is therefore ontologically cleaner: particle–particle exchange (difference-biased) is separated from particle–field exchange (sum-biased). Bill's description of the proposal as involving *"reversible exchange of a quanta of momentum between co-located world-particles"* is accurate for the Focus/Defocus channel only; the hop channel cannot be a particle–particle exchange without violating momentum conservation, and should be understood as photon absorption/emission.

### 6.2 The mediation is relocated, not eliminated

Here a correction to the framing "just two-body interactions" is needed. A genuinely autonomous pairwise interaction — mass-action kinetics, where the rate of a two-particle event is proportional to the **product** of the participants' occupancies — cannot reproduce the QLE:

**No-go lemma.** If every event rate is a multilinear monomial of degree $\ge 2$ in the occupancies (mass action for two or more participating bodies), the mean-field generator is a polynomial in $W$ of degree $\ge 2$ with no linear part, and cannot equal the linear stencil (1) for all admissible $W$.

This is the same linearity obstruction previously established for deterministic microdynamics (the QLE collision term is linear in $W$, ruling out purely pairwise deterministic interactions and forcing a three-body structure with an explicit mediating Fourier mode). The four-rule scheme evades it the only way possible: its event rates are **linear in the occupancies of cells other than, or in addition to, the moved quanta** — the focus rate at center $n$ depends on $W_{n\pm q}$, the hop rate on the occupancy of cells the hopping quantum never touches. Operationally, each "two-body" event consults the local ensemble density. The information carried by the mediating Fourier mode in the three-body picture has moved into the rate law; it has not disappeared. (The known loophole is equivalent: mass action with one participant drawn from a *reservoir species of pinned density* — which is precisely what a mediating mode with fixed amplitude is.)

This is not a defect — the original single rule has exactly the same character (its rate is linear in the mediator occupancy, and the transferred quantum is moved unconditionally, regardless of the source cell's occupancy). But it means the four-rule proposal should be described as a *re-decomposition* of the mediated interaction with better momentum bookkeeping, not as an elimination of mediation.

---

## 7. In what sense are the schemes equivalent?

**Exact equivalence in expectation, at any particle number.** Because every event rate in both schemes is *linear* in the occupancies and every event applies a fixed integer stencil, the evolution of the expected occupancy field is closed: taking expectations of the master equation produces (2) with $f, h$ evaluated at the expected field, with **no** BBGKY-type hierarchy and no factorization approximation. Hence for any $\nu$, the expected Wigner field of the four-rule process and of the single-rule process obey the *same* closed linear ODE — the QLE stencil (1) — and agree exactly at all times if started from the same initial expectation. This, together with §4, is the equivalence proof: both generators equal the QLE generator, hence each other.

**Not pathwise equivalence.** The two processes are different stochastic processes: their event sets, and therefore their fluctuation (second-moment) structures, differ. Three concrete differences:

1. **Background noise.** The single rule as specified counts the uniform $2/h$ background among the mediators, so it fires events at rate $\propto |\Gamma| \thinspace W'$ everywhere $\Gamma \ne 0$, including cells of zero excess — pure noise with zero mean. The four-rule bias rates, being differences and sums of *excess* occupancies (the background cancels identically in the focus rate and is a null bias in the hop rate), fire only where excess structure exists. At equal $\nu$ the four-rule process is empirically \~5–6× quieter (§9). The comparison is not entirely fair to the original — nothing prevents restricting its mediator count to the excess as well, which is a null-bias change — but the four-rule formulation makes the quiet choice automatically.
2. **Momentum-noise character.** Single-rule events displace one quantum by $2q$ cells; focus/defocus events displace two quanta by $q$ cells each in opposite directions. Per unit of generated mean drift the injected momentum-space variance differs, and the family parameter $G$ in (4) is precisely a **variance-shaping gauge freedom**: all members share the mean, none share the noise. Choosing $G$ to minimize variance for a given state is an open optimization (§10).
3. **Boundedness.** Both schemes cap transfers by source occupancy to keep $N_+ \ge 0$; the caps bind on different events, producing (rare, $O(1/\nu)$) different boundary behavior.

---

## 8. Algorithm

The four-rule algorithm is a drop-in replacement for spec §3b/§6, sharing the lattice, the state representation $N_+ = \nu\thinspace(W + 2/h)\thinspace\Delta x\thinspace\Delta p$, the free-streaming step, and observable extraction. Only the jump step changes. The symmetric member (5) is specified; any member of family (4) may be substituted.

**Per timestep, per Fourier mode $q$, per position column $m$** (with $\Gamma = \Gamma_q(x_m)$, background count $B = \mathrm{round}(\nu\thinspace(2/h)\thinspace\Delta x\thinspace\Delta p)$, excess $\tilde N_n = N_+(n, m) - B$):

1. Compute the two signed bias rates for every momentum center $n$:

```math
F_n = \frac{\Gamma}{2}\bigl(\tilde N_{n+q} - \tilde N_{n-q}\bigr),
\qquad
H_n = -\frac{\Gamma}{2}\bigl(\tilde N_{n+q} + \tilde N_{n-q}\bigr).
```

2. Draw net event counts by tau-leaping: $e^F_n = \mathrm{sgn}(F_n)\cdot\mathrm{Poisson}(|F_n|\thinspace\Delta t)$ and likewise $e^H_n$. (Equivalently, draw forward and backward channels separately — Focus and Defocus at rates $\max(\pm F_n, 0)$ plus any common base rate, etc.; the difference is a null noise term. In zero potential $\Gamma \equiv 0$, all rates vanish or cancel: the slide's "inverse pairs exactly cancel".)

3. Apply, capping by source occupancy:
   - $e^F_n > 0$ (Focus): remove $e$ from each of $n \pm q$, add $2e$ at $n$, with $e$ capped by $\min(N_{n-q}, N_{n+q})$.
   - $e^F_n < 0$ (Defocus): remove $2|e|$ from $n$, add $|e|$ to each of $n \pm q$, capped by $\lfloor N_n / 2 \rfloor$.
   - $e^H_n > 0$ (Right-Hop): move $e$ from $n-q$ to $n+q$, capped by $N_{n-q}$.
   - $e^H_n < 0$ (Left-Hop): move $|e|$ from $n+q$ to $n-q$, capped by $N_{n+q}$.

Stability requires $`|\Gamma|_{\max}\,\Delta t \ll 1`$ exactly as for the original rule. Cost per step is the same order as the original MC form. The implementation is `PhaseSpaceCrystalLattice.step_jump_four_rule_mc`; the mesh-density (deterministic) form, with the channels assembled independently so that it constitutes a nontrivial check of identity (3), is `step_jump_four_rule`.

---

## 9. Numerical verification

`src/demo_four_rule_equivalence.py` performs two tests.

**Part A — generator identity.** Random band-limited occupancy fields, mode sets $\lbrace q{=}1 \rbrace$, $\lbrace q{=}2 \rbrace$, $\lbrace q{=}1,2,3 \rbrace$: one Euler step of the four-rule mesh form vs. the original stencil. Worst absolute deviation over all trials: $1.1 \times 10^{-16}$ (machine epsilon), relative deviation $\sim 10^{-14}$.

**Part B — head-to-head stochastic evolution.** Squeezed Gaussian ($\sigma_x = 1.3\thinspace\sigma_{x,\mathrm{gs}}$) in the single-period cosine well $V(x) = -V_p\cos(2\pi x/L)$, $V_p = 1.5$, $L = 8$, lattice $64 \times 64$, $\nu = 1.6 \times 10^6$, two classical periods. Three evolutions with an identical splitting sequence (exact-integer advection macro-steps $\Delta t = m\thinspace\Delta x/\Delta p$; 16 jump sub-steps each, $`|\Gamma|_{\max}\Delta t_{\mathrm{jump}} = 0.03`$): deterministic mesh QLE, single-rule MC, four-rule MC.

![Four-rule equivalence evolution](https://raw.githubusercontent.com/billpage/wpmw/output/figures/four_rule_equivalence_evolution.png)

Rows: mesh QLE, single-rule MC, four-rule MC, and pointwise (four-rule − mesh). Both MC evolutions track the mesh; the single-rule row shows the expected uniform background shot noise (its mediator count includes the $2/h$ background), the four-rule row is visibly cleaner, and the difference row is structureless noise confined to where the excess lives.

![Four-rule equivalence metrics](https://raw.githubusercontent.com/billpage/wpmw/output/figures/four_rule_equivalence_metrics.png)

Left: relative $L^2$ deviation from the mesh, growing diffusively for both MC runs (accumulating Poisson noise), with the four-rule floor \~5.6× lower at equal $\nu$ (0.071 vs. 0.402 at $t = 2T$). Right: Wigner negativity $\int|\min(W, 0)|\thinspace dx\thinspace dp$ — the mesh accumulates genuine negativity; the four-rule MC tracks it closely; the single-rule MC rides a much higher noise floor, illustrating §7 point 1.

Both MC deviations behave as shot noise: rerunning at $\nu = 6.4\times 10^4$ (25× fewer particles) raises the final floors by 5.0× (single-rule) and 6.5× (four-rule), against the ideal $\nu^{-1/2}$ factor of 5 — consistent with exact agreement in expectation, with the four-rule run showing a mild extra low-density contribution (plausibly the source-occupancy caps, which bind more often when excess counts are small).

---

## 10. Open items

1. ~~**Page 4 of the slide deck.** Confirm the second signed bias rate and the intended indexing of $r_k$ (§5). The prediction from the theorem: the hop-channel bias is $-\tfrac{\kappa}{2}(W_{k+1} + W_{k-1})$ under the center-indexed reading.~~ **Resolved, August 2026** — confirmed verbatim, along with the center-indexed reading; see §0 and [`../supplement/four_action_foundations.md`](../supplement/four_action_foundations.md) §1. A further result recorded there bears on item 2: conservation laws are satisfied by *every* member of family (4), so the $G$-freedom is not constrained by them, and what selects the symmetric member is endpoint locality — no rate may reference the centre cell.
2. **Variance-optimal member.** The gauge freedom $G$ in (4) shapes the injected noise at fixed mean. Minimizing MC variance over $G$ (state-dependently, or for a worst case) is a well-posed and likely fruitful optimization; the excess-vs-full-population choice for the original rule is the zeroth-order instance.
3. **Self-consistent N-body case.** This analysis, like the original spec, treats the external fixed-cosine potential. Whether the focus channel's momentum-conserving two-body structure survives — or helps — when the potential is generated self-consistently by the world-particle ensemble connects directly to the open mode-reduction question from the deterministic-microdynamics work.
4. **Fermionic/multi-DOF extension.** The channel split (particle–particle vs. particle–field) may interact nontrivially with the sign conventions of `multi_body_extension.md` §12.
