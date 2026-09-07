# Holland's Two-Fluid Models and the Signed Ensemble

**A comparison note. Peter Holland has published three exact trajectory
representations of the Schrödinger equation, one of which is a two-fluid model
whose partition parameter plays the same role as this project's absorptive
fraction. This note establishes what corresponds to what, derives the one
correspondence that is more than an analogy, and records where the two
frameworks part company.**

---

## 0. What this note is

Expository and comparative. It asserts nothing normative about the algorithm
and introduces no changes to the specification. Its purpose is to place the
signed world-particle ontology against an independent construction that reaches
a similar structure from a different direction, so that the parts of this
project which are genuinely new can be told apart from the parts which are
rediscoveries.

Propositions here carry the prefix `H`. None is load-bearing for the ladder.

---

## 1. Three models, and which one is meant

Holland has given three distinct trajectory constructions, and they differ
exactly on the point at issue here — whether anything is conserved.

| model | fields | conservation | reference |
|---|---|---|---|
| polar | $`\rho`$, $`S`$ | one continuity equation | Ann. Phys. **315**, 505 (2005) |
| real/imaginary | $`\psi_R`$, $`\psi_I`$ | **two** continuity equations | J. Phys. A **42**, 075307 (2009) |
| bi-Hamilton–Jacobi | $`S_{+}`$, $`S_{-}`$ | **none** | Found. Phys. **52**, 100 (2022) |

The 2009 model is titled *Schrödinger dynamics as a two-phase conserved flow*
and it is the one most easily confused with what follows. There the two fluids
have signed densities $`\psi_R`$ and $`\psi_I`$, each obeys its own continuity
equation, and the two species are separately conserved. **Nothing interconverts
in the 2009 model.** The coupling is entirely through the velocity fields,
$`\mathbf{v}_1 = (\hbar/2m)\nabla\psi_2/\psi_1`$ and
$`\mathbf{v}_2 = -(\hbar/2m)\nabla\psi_1/\psi_2`$.

The model relevant to this project is the third. It is the one in which Holland
remarks that the coupled form of the equations is suggestive of a two-phase
system exhibiting continual conversion of one species into the other, citing the
continuum theory of reacting mixtures.

---

## 2. The bi-Hamilton–Jacobi system

Write the wavefunction in polar form and define

```math
S_{\pm} \;=\; S \;\pm\; (\hbar/2)\log\rho .
```

Then $`\psi`$ factors into two complex amplitudes, neither of which satisfies
the Schrödinger equation, and the Schrödinger equation itself becomes two real
equations closed in the fields $`S_{\pm}`$:

```math
\frac{\partial S_{\pm}}{\partial t}
  + \frac{1}{2m}\partial_i S_{\pm}\thinspace \partial_i S_{\pm}
  + Q_{\pm} + V \;=\; 0,
\qquad
Q_{\pm} \;=\; \mp\frac{\hbar}{2m}\partial_{ii}S_{\mp}
  \;-\; \frac{1}{4m}\big[\partial_i(S_{+}-S_{-})\big]^{2} .
```

Two features matter here. First, there are **no continuity equations and no
independent density functions** — the theory is written entirely in velocity
fields $`v_{\pm}^{i} = \partial_i S_{\pm}/m`$. Second, the probability density
is recovered not from any population but from the *difference* of the two
actions:

```math
\rho \;=\; \exp\big[(S_{+} - S_{-})/\hbar\big] .
```

**A bridge worth flagging.** Holland identifies $`v_{\pm}`$ explicitly as
Nelson's forward and backward drift velocities:
$`u^i = v_{+}^i - v_{-}^i = (\hbar/m)\partial_i\log\rho`$ is the osmotic
velocity and $`\tfrac{1}{2}(v_{+}^i + v_{-}^i)`$ is the de Broglie–Bohm
velocity. This connects the bi-HJ system directly to the Nelson material in
[`four_action_foundations.md`](four_action_foundations.md).

Note also what $`v_{\pm}`$ are *not*. Both $`\langle p\rangle(x)`$ and the
osmotic term are properties of the $`x`$-marginal — they exist because momentum
has been integrated out. In a phase-space representation that carries $`p`$
explicitly there is nothing for them to do.

---

## 3. The free sector

The correspondence with this project lives in what neither framework's
observable can see.

In the signed ensemble the observable is $`E = u^{+} - u^{-}`$ and the sum
$`N = u^{+} + u^{-}`$ is unconstrained; the padding ratio $`\rho_{\rm pad}`$
may take any value without touching $`E`$
([`emission_and_absorption.md`](emission_and_absorption.md) §1).

In the bi-HJ model the two congruences may be assigned Lagrangian trajectory
densities $`\rho_0 J_{\pm}^{-1}`$, and Holland proves these do **not** generate
$`\rho`$ — not individually and not in combination. Splitting the initial
density between the congruences with weights $`r`$ and $`1-r`$,

```math
\rho \;=\; r\big(\rho_0 J_{+}^{-1} + c_{+}\big)
  \;+\; (1-r)\big(\rho_0 J_{-}^{-1} + c_{-}\big),
```

and the source terms $`c_{\pm}`$ do not cancel. The congruence densities are
free data.

**Proposition H1 (the same freedom).** In both frameworks the observable is a
difference and the corresponding sum is unconstrained; $`r`$ and the absorptive
fraction $`f`$ are both partition parameters of a sector the observable cannot
resolve.

---

## 4. The one correspondence that is a derivation

With respect to each flow separately, $`\rho`$ obeys a Fokker–Planck equation
with opposite-signed diffusion:

```math
\frac{\partial\rho}{\partial t} + \partial_i(\rho v_{\pm}^{i})
  \;=\; \pm\frac{\hbar}{2m}\partial_{ii}\rho .
```

Weighting the two by $`r`$ and $`1-r`$ and adding:

```math
\frac{\partial\rho}{\partial t}
  + \partial_i\big(\rho\thinspace[ r v_{+}^{i} + (1-r)v_{-}^{i}]\big)
  \;=\; (2r-1)\thinspace\frac{\hbar}{2m}\partial_{ii}\rho .
```

**Proposition H2.** At $`r = 1/2`$ the source vanishes identically and the
bracket is $`\tfrac{1}{2}(v_{+}+v_{-}) = \nabla S/m`$, so the mixture obeys the
exact continuity equation. At no other value does it.

Compare Theorem S7 of
[`../analysis/sea_population_equilibrium.md`](../analysis/sea_population_equilibrium.md),
$`\Delta N = 2(1-2f) n_{\rm ev}`$ and
$`\Delta S = (2f-1) n_{\rm ev}`$: both ledgers close at $`f = 1/2`$
and at no other value.

**Proposition H3.** Both stationarity conditions have the form
$`(2\lambda - 1)\times(\text{rate})`$ for a partition parameter $`\lambda`$,
and in both frameworks $`\lambda = 1/2`$ is the fixed point of the exchange
symmetry relating the two channels. The correspondence $`r \leftrightarrow f`$
is therefore structural rather than suggestive.

**A difference in the project's favour.** Holland must *stipulate* $`r`$; any
value is consistent with his formalism, because the congruence densities are
unobservable. Theorem S9 measures $`f`$ as an **attractor**, reached from both
sides without tuning. If the correspondence holds, this project supplies a
mechanism for a number Holland can only choose.

---

## 5. Time reversal

Holland's time reversal is complex conjugation, $`\psi \to \psi^{*}`$, which in
these variables reads

```math
S'_{\pm}(x',t') = -S_{\mp}(x,t),
\qquad
v'^{i}_{\pm}(x',t') = -v^{i}_{\mp}(x,t) .
```

Each flow individually violates the conventional transformation; the **pair** is
covariant, because the time-reversal transform of one flow is the T-reversal of
the other. Holland calls this T-symmetry implemented through the collective
behaviour of elements that individually disobey T.

**Proposition H4 (what this is not).** The exchange is *not* the analogue of
positon/negaton exchange. That exchange is charge conjugation; it sends
$`E \to -E`$, and $`-W`$ is not a state, so it is not a symmetry of the
dynamics at all. The populations $`N`$ and $`S`$ are C-even and $`E`$ is C-odd.

**Proposition H5 (what it is).** The analogue is the emissive/absorptive pair.
Read the space-time figure of [`emission_and_absorption.md`](emission_and_absorption.md)
§4 backwards: a fork becomes a merge with every slope negated. So T composed
with $`p \to -p`$ exchanges the two realisations, which is exactly Holland's
structure. Neither channel is individually T-covariant; the pair is.

This sharpens §4's sentence that the two realisations are not each other's time
reverse. That is correct as stated — plain exchange is not T, because nothing
flips momentum — but the pair is T-covariant in Holland's sense, and the
sharper statement is the useful one.

**Conjecture H6.** If T exchanges the channels, then $`f = 1/2`$ is the
T-symmetric point of the event ensemble and Theorem S9's attractor is
symmetry-restoring, with the shortfall $`\tfrac{1}{2} - f`$ a measure of
T-breaking by the reach regulator. This bears on open item S-SP3. It is
*not* established: $`f`$ lives in the sector the observable cannot see, so
T-invariance of $`E`$ does not by itself force T-invariance of $`f`$.

---

## 6. Where the frameworks part company

**What this project supplies that Holland lacks: an inventory.** Holland's
"continual conversion" is a remark about the coupled form of the velocity
equations, and his own text is explicit that the theory has no independent
densities. There is nothing to convert. The signed ensemble supplies exactly
the missing object — a third population and a conserved total
$`P = S + N/2`$ (Theorem J1) — so that conversion has a stock to draw on.

**What Holland supplies that this project lacks: locality in the conversion.**
In the bi-HJ system the leak $`(2r-1)(\hbar/2m)\partial_{ii}\rho`$ that one flow
sheds and the other absorbs sit at the **same** $`x`$. They cancel pointwise, not
merely in the integral. The event channel does not have this property: it debits
the sea at the parent momentum row and credits bodies at the daughter rows
$`p \pm \xi_q`$, so pair count flows between rows and cancellation is global
only (Theorem J2).

**Proposition H7.** The local sea deficit is therefore not a generic
consequence of a signed or two-fluid representation of the Schrödinger
equation. It is a consequence of the momentum offset $`\pm\xi_q`$ — of working
in phase space with a finite reach — and Holland's position-space construction
has no analogue of it.

This locates open item CLA8 more precisely than "the price of the ontology."

**A construction that makes the tension explicit.** Introduce a third field
$`\Sigma`$ and write $`\rho = \exp[(S_{+} - S_{-} + \Sigma)/\hbar]`$. The
Schrödinger equation holds if and only if $`\Sigma \equiv 0`$, so any
non-trivial dynamics for $`\Sigma`$ breaks it. This is the CLA8 tension in a
formalism with no mesh, no reach and no tau-leap, which suggests the *tension*
— a floored population against an unconstrained gauge sector — is structural
even though the *deficit* by H7 is not.

---

## 7. Open items

### H-SP1 — is the $`\sigma \to 0`$ limit literally Theorem S5's throttle?

Meeting notes report David Cyganski reaching a $`\sigma \to 0`$ limit that
returns the Wigner equation. Theorem S5 throttles the rate as
$`\Gamma \to \Gamma S/B`$, so $`\sigma = 0`$ gives the QLE and
$`\sigma \ne 0`$ does not. If these are the same statement, the two frameworks
are the same equations in different variables, and that is worth establishing
independently of everything above.

### H-SP2 — is the reservoir construction unique?

§6's third-field construction is *consistent* with Holland's equations but not
*derived* from them, and because his density sector is unconstrained, any such
scheme is consistent. Consistency is therefore not evidence. A uniqueness
argument would be: require the reservoir to be stationary exactly when
probability is conserved, which forces the $`(2r-1)`$ form. Whether that
condition is natural enough to carry the weight is not settled here.

### H-SP3 — the free Gaussian as a closed-form test of H7

Holland gives the bi-HJ paths for a free Gaussian in closed form. The
accumulated reservoir drain for $`r \ne 1/2`$ can then be computed analytically,
with no mesh and no reach, and compared against Theorem S4. This is the
cheapest available test of H7.

---

## 8. Sources

- P. Holland, *Eliminating the wavefunction from quantum dynamics: the
  bi-Hamilton–Jacobi theory, trajectories and time reversal*,
  Found. Phys. **52**, 100 (2022); <https://arxiv.org/abs/2111.09235>.
- P. Holland, *Schrödinger dynamics as a two-phase conserved flow*,
  J. Phys. A **42**, 075307 (2009); <https://arxiv.org/abs/0807.4482>.
- [`emission_and_absorption.md`](emission_and_absorption.md) — the three
  populations, Theorems J1 and J2, the two realisations.
- [`../analysis/sea_population_equilibrium.md`](../analysis/sea_population_equilibrium.md)
  — Theorems S2, S5, S7, S9 and open item S-SP3.
- [`../algorithm/compensated_liouville_algorithm.md`](../algorithm/compensated_liouville_algorithm.md)
  §5 — the specification and open item CLA8.
