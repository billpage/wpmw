# Sea population equilibrium: what a signed channel costs the ledger

**Status.** Analysis note, step 16 of the ladder. Companion demo:
`src/demo_sea_population_equilibrium.py`. Prompted by open item **CLA3** of
[`../algorithm/compensated_liouville_algorithm.md`](../algorithm/compensated_liouville_algorithm.md)
("§4.3 leaves the positon/negaton unraveling specified but unpriced"), and by
the standing question of whether a recombination rule is needed at all.

---

## 0. What this note settles, and what it corrects

The compensated split leaves $`K_{\rm res}`$ signed, so the sign has to go
somewhere. The specification recommends the positon/negaton reading and states
the world-form rate as $`R(x) = \sum_q |K_{\rm res}(x, \xi_q)|`$. Taken
literally that prescribes an **emissive** unravelling: every event creates a
positon–negaton pair. This note prices that choice, finds it ruinous, and
identifies the alternative.

**Settles.**

- **CLA3.** The emissive unravelling is not merely expensive. It drains the
  sea without bound (S4), and if the sea is finite the resulting throttle
  moves $`W`$ by 40 per cent in the core (S5). The **absorptive** unravelling
  — realising a deposition by removing a body of the opposite species rather
  than adding one — costs nothing in the observable, bounds the body count,
  and shrinks the deficit by two to three orders of magnitude (S8).
- **Whether a recombination rule is needed.** Yes, necessarily, and it is not
  derivable from the Moyal equation: the equation fixes $`u^+ - u^-`$ and says
  nothing about $`u^+ + u^-`$. What *is* derivable is the rule's form (S1,
  §7) and the exact condition under which the ledger closes (S7).

**Corrects.**

1. [`../algorithm/compensated_liouville_algorithm.md`](../algorithm/compensated_liouville_algorithm.md)
   §4.3 gives one world-form rate and one update rule. There are two
   realisations of every event with identical observable content and opposite
   ledger content, and the spec's wording selects the bad one by default. §5.
2. [`eckart_barrier_compensated.md`](eckart_barrier_compensated.md) §8.1 says
   a positon–negaton pair "appears" at $`p \pm \xi_q`$. Momentum conservation
   makes this precise: the pair must be a **bound sea pair co-located with the
   parent**, since the children carry $`2p`$ and only a pair at the parent's
   own momentum row carries $`2p`$ too. §3.
3. The premise, used informally in discussion, that the sea is much denser
   than the excess. For pure states $`|W| \le 2/h = B`$ pointwise, saturated at
   the centre of every Gaussian, so the excess *equals* the crystal shift
   exactly where the dynamics is busiest. §2.

---

## 1. The ledger, and why it is not determined by the dynamics

Three fields per phase-space cell, in Wigner units:

| field | meaning |
|---|---|
| $`E = u^+ - u^-`$ | the signed density; the observable, $`E \equiv W`$ |
| $`N = u^+ + u^-`$ | the body density |
| $`S`$ | bound sea pairs, background $`B = 2/h`$ |

The Moyal equation constrains $`E`$ and nothing else. Any redistribution of
$`(u^+, u^-)`$ at fixed difference is consistent with every observable, so the
population is genuinely extra structure — a choice, not a consequence. The
whole content of this note is what constrains that choice.

**Theorem S1 (the sink must be bilinear).** A removal channel preserves $`E`$
for arbitrary ensembles if and only if it removes positons and negatons in
coincident pairs.

*Proof.* A per-body death rate $`\mu`$ applied to both species sends
$`E \to e^{-\mu t}E`$, which changes the observable. Species-selective removal
changes $`E`$ directly. Only removal in $`\pm`$ pairs at the same cell leaves
$`u^+ - u^-`$ invariant, and pair removal is quadratic in the local
populations. $`\square`$

**Corollary S1.1.** Proposition U1 of
[`species_sectors_and_annihilation.md`](species_sectors_and_annihilation.md)
— no unravelling *linear* in the ensemble is $`L^1`$-stationary, with
governing exponent $`\rho(|L|) = 2.341`$ — does not obstruct the stabiliser.
It identifies the class the stabiliser must lie outside, and S1 says the same
thing from the other side.

---

## 2. The excess is not small against the sea

For any normalised pure state $`W = (2/h)\langle \Pi \rangle`$ with $`\Pi`$ the
phase-space parity operator, so $`|W| \le 2/h`$ pointwise. The bound is
saturated at the centre of a Gaussian of any width — to five figures on a
momentum grid fine enough to resolve the peak — and reached at $`-0.827`$ of it
in the interference trough of a cat state. Part B reports the ratio for the
demo's own packet, 0.969 on its coarser grid.

Two consequences. The crystal shift $`B = 2/h`$ is the *tightest* uniform shift
making $`W' \ge 0`$ for all pure states, which derives the value 2 rather than
adopting it — the same inequality already noted in
[`permanent_pairing_density_matrix.md`](permanent_pairing_density_matrix.md).
And the "dilute excess in a dense sea" picture is false at the packet core,
which is exactly where the emission rate is largest. Any argument that relies
on a deep reservoir to justify pinning is therefore unavailable.

---

## 3. What an event is

The residual operator is a convolution in momentum: a parent at $`p`$ deposits
$`+1`$ at $`p + \xi_q`$ and $`-1`$ at $`p - \xi_q`$, since $`K_{-q} = -K_q`$.
The $`(q, -q)`$ pair is therefore **one event with two legs of opposite sign**,
not two separate channels.

**Proposition S0 (co-location is forced).** If the event conserves body
momentum and the parent streams on undisturbed, the bound pair it consumes
must sit at the parent's own momentum row.

*Proof.* The children carry $`(p + \xi_q) + (p - \xi_q) = 2p`$. A pair drawn
from row $`p'`$ carries $`2p'`$. Equality forces $`p' = p`$. Nothing else can
supply the difference, because by §8.1 of
[`eckart_barrier_compensated.md`](eckart_barrier_compensated.md) no
world-particle changes its momentum discontinuously, the parent included.
$`\square`$

This is what makes "ionisation of a bound sea pair" the right reading of §8.1
rather than creation from nothing, and it is a derivation, not a postulate.

---

## 4. The two-field split, and the fixed point

Under the emissive unravelling and a pinned sea the ledger splits:

```math
\partial_t E + \{E, H\}_{\rm cl} = \sum_q K_q(x)\, E(p - \xi_q)
```

```math
\partial_t N + \{N, H\}_{\rm cl} = \sum_q |K_q(x)|\, N(p - \xi_q)
  \;-\; \tfrac{\kappa}{2}\bigl(N^2 - E^2\bigr)
```

**Theorem S2 (the ledger cannot corrupt the observable, and has a closed-form
fixed point).** The $`E`$ equation is closed and is exactly the QLE, whatever
$`N`$ does. The $`N`$ equation runs under the *absolute-value* kernel, so its
growth exponent is $`\rho(|L|)`$, and its local fixed point is

```math
N_{\rm eq} \;=\; \frac{\Gamma_{\rm tot}}{\kappa}
  + \sqrt{\frac{\Gamma_{\rm tot}^{2}}{\kappa^{2}} + E^{2}},
\qquad
\Gamma_{\rm tot}(x) = \sum_{q \ne 0} |K_q(x)|,
```

unique and globally attracting, with $`N_{\rm eq} \to |E|`$ as
$`\kappa \to \infty`$ and $`N_{\rm eq} \to 2\Gamma_{\rm tot}/\kappa`$ as
$`\kappa \to 0`$.

Part A verifies the closed form to $`\sim 10^{-13}`$ from both directions over
$`\kappa \in [0.5, 1000]`$. The $`\kappa \to \infty`$ limit is worth stating
separately: **the minimal, pure-positon ensemble is the fixed point of fast
recombination**, so preparing a positive $`W`$ as positons alone is a derived
choice and not a convention.

**Theorem S3 (two relaxation laws).** Inside the emitting region the fixed
point is approached exponentially at rate $`\sim \kappa|E|`$. Where
$`\Gamma = E = 0`$ the ledger obeys $`dN/dt = -(\kappa/2)N^2`$, so
$`N \to 2/(\kappa t)`$ — a universal asymptote reached only algebraically.

So preparation *is* forgotten, but the vacuum forgets as $`1/t`$, and at any
finite time it still carries a trace of how it was prepared.

**Reach independence fixes the form of $`\kappa`$.** At fixed $`\kappa`$,
$`N_{\rm eq}/B`$ runs 1.05, 1.19, 1.60, 2.46 as $`y_{\max}`$ doubles from
$`\pi a/2`$ to $`4\pi a`$; with $`\kappa = c\,\Gamma_{\rm tot}(x)`$ it is
constant. The recombination rate must be proportional to the local emission
rate, which leaves exactly one dimensionless constant free. This is the same
conclusion as the design caution of
[`species_sectors_and_annihilation.md`](species_sectors_and_annihilation.md)
§11 ("the recombination rate must not reference $`B`$"), reached from the
reach side rather than the cost side.

---

## 5. The emissive unravelling fails

**Theorem S4 (the sea is relocated, without bound).** Under the emissive
unravelling $`S`$ is debited at the parent row $`p`$ and credited at the
daughter rows $`p \pm \xi_q`$, where the pair eventually recombines. The mean
occupancy is conserved — bodies are conserved — while the worst cell drains
monotonically:

| $`t`$ | 2 | 4 | 8 | 12 | 16 | 20 |
|---|---|---|---|---|---|---|
| min $`S/B`$ | 0.25 | $`-0.39`$ | $`-1.93`$ | $`-2.76`$ | $`-3.22`$ | $`-4.84`$ |
| mean $`S/B`$ | 0.999 | 0.999 | 0.998 | 0.998 | 0.997 | 0.997 |

The deficit grows roughly linearly in $`t`$, so **no finite reservoir depth
suffices**: a sea ten times deeper buys ten times the runtime, not a fix.
Increasing $`\kappa`$ does not help, because the credit never returns to the
row that was debited.

**Theorem S5 (a supply-limited sea is observable).** If the rate is throttled
by the actual occupancy, $`\Gamma \to \Gamma \cdot S/B`$, then the $`E`$
equation is no longer the QLE. Against the pinned reference at $`T = 8`$,
reach $`2\pi a`$:

| $`\kappa`$ | rel $`L^2`$ in $`E`$ | in the core | in the far field | $`\Delta`$ norm | $`\Delta\langle p\rangle`$ |
|---|---|---|---|---|---|
| 200 | 0.402 | 0.401 | 0.026 | $`9\times10^{-16}`$ | $`2\times10^{-7}`$ |
| 2000 | 0.414 | 0.413 | 0.024 | $`2\times10^{-15}`$ | $`3\times10^{-7}`$ |
| 20000 | 0.415 | 0.414 | 0.024 | $`5\times10^{-15}`$ | $`3\times10^{-7}`$ |

The error lives where the state lives, saturates in $`\kappa`$ rather than
vanishing, and grows with reach (0.128, 0.415, 0.536 at
$`y_{\max}/a = \pi, 2\pi, 4\pi`$). Norm and $`\langle p\rangle`$ are protected
to machine precision, which is the worst possible failure mode: the cheap
diagnostics do not see it.

---

## 6. The other realisation

Every deposition has two realisations, identical in $`E`$:

| deposition | emissive | absorptive |
|---|---|---|
| $`+1`$ at a cell | add a positon | remove a negaton |
| $`-1`$ at a cell | add a negaton | remove a positon |

They differ entirely in the ledger. An emissive event splits a bound pair at
the parent row: bodies $`+2`$, $`S`$ debited at $`p`$. An absorptive event
binds the two removed bodies into a pair at the parent row: bodies $`-2`$,
$`S`$ **credited at $`p`$**. Both debit and credit now land on the same row, so
the transport asymmetry of S4 disappears by construction.

**Theorem S6 (the mode is per event, not per leg).** A mixed realisation —
one leg created, the other removed — moves a body from $`p - \xi_q`$ to
$`p + \xi_q`$ and so changes body momentum by $`2\xi_q`$. Nothing can supply
it (Proposition S0). Hence an event is wholly emissive or wholly absorptive,
and absorption requires a partner of the right species at **both** daughters.
Absorption is therefore supply limited, and when supply fails the event falls
back to emissive.

**Theorem S7 (the ledger identity).** With absorptive fraction $`f`$ over
$`n_{\rm ev}`$ events,

```math
\Delta N = 2\,(1 - 2f)\, n_{\rm ev},
\qquad
\Delta S = (2f - 1)\, n_{\rm ev},
```

so $`f = 1/2`$ closes the body ledger and the sea ledger **simultaneously**.
Part E verifies both to $`1.5\times10^{-14}`$ and $`1.7\times10^{-11}`$
respectively, at a measured $`f = 0.4344`$.

This is the sharpest statement the note reaches. The two ledgers are not
independent problems with two knobs; they are one problem with one number.

**Theorem S8 (absorptive unravelling).** Eckart barrier, packet on the summit,
$`T = 6`$, pinned rate. Sea columns are the cumulative diagnostic ledger:

| $`\Delta t`$ | mode | $`N(T)`$ | $`f`$ | min $`S/B`$ | mean $`S/B`$ | rel $`L^2`$ vs QLE |
|---|---|---|---|---|---|---|
| 0.02 | emissive | $`1.02\times10^{5}`$ | 0 | $`-7.4\times10^{4}`$ | $`-246`$ | $`1.78\times10^{2}`$ |
| 0.01 | emissive | $`1.07\times10^{5}`$ | 0 | $`-7.8\times10^{4}`$ | $`-259`$ | $`9.40\times10^{1}`$ |
| 0.02 | absorptive | 2.320 | 0.4354 | $`-1.263`$ | 0.99676 | $`1.96\times10^{-2}`$ |
| 0.01 | absorptive | 2.330 | 0.4344 | $`-1.235`$ | 0.99674 | $`9.69\times10^{-3}`$ |
| 0.005 | absorptive | 2.340 | 0.4338 | $`-1.220`$ | 0.99671 | $`4.83\times10^{-3}`$ |

Nothing is throttled, so by construction $`E`$ should be the exact QLE — and
it is: the error halves cleanly with $`\Delta t`$, so it is pure splitting
error converging to zero. The deficit is $`\Delta t`$-converged at
$`-1.22\,B`$, so that residual is physical rather than numerical.

**The fidelity gap is amplification, not a different equation.** Both
unravellings are first-order convergent and both integrate the same $`E`$. The
emissive ledger computes $`E`$ as a small difference of two populations of
size $`10^{5}`$, so every integration error is multiplied by
$`N/\lVert E \rVert`$. This is the sign problem appearing as *integration*
error rather than as variance, and absorption removes it by keeping $`N`$ at
2.3 instead of $`10^{5}`$.

Over a longer run the two diverge qualitatively rather than quantitatively:

| $`t`$ | 2 | 6 | 10 | 14 | 18 |
|---|---|---|---|---|---|
| emissive, min $`S/B`$ | $`-26`$ | $`-7.9\times10^{4}`$ | $`-1.9\times10^{8}`$ | $`-4.6\times10^{11}`$ | $`-4.5\times10^{14}`$ |
| absorptive, min $`S/B`$ | $`-0.61`$ | $`-1.00`$ | $`-1.43`$ | $`-2.41`$ | $`-1.65`$ |

The absorptive deficit fluctuates and grows sub-linearly; the emissive one is
exponential. Global sea loss under absorption is 0.6 per cent over 24 time
units.

---

## 7. Where this leaves the recombination rule

Collecting §§1–6, the rule is constrained on four sides and free on one:

| constraint | source | what it fixes |
|---|---|---|
| bilinear in $`(u^+, u^-)`$ | S1 | the functional form of the sink |
| removes coincident pairs only | S1 | the reach of the removal |
| $`\kappa \propto \Gamma_{\rm tot}(x)`$ | §4 | reach independence of the ledger |
| $`f = 1/2`$ | S7 | closure of both ledgers |
| one dimensionless constant | — | **not fixed** |

The one remaining constant is the same object as the choice of initial
ensemble, because $`N_{\rm eq}`$ depends on $`\Gamma/\kappa`$ and on nothing
else about the state. Whether one principle fixes both is §8's first open item.

---

## 8. Open items

- **S-SP1 (the supra-minimal ensemble).** The measured $`f = 0.434`$ falls
  short of the closing value $`1/2`$ because partners are not always available
  — and they are not available because the runs start from the *minimal*
  ensemble $`N_0 = |E|`$. Theorem S2's fixed point is strictly above $`|E|`$
  whenever $`\kappa`$ is finite, and that excess is exactly the standing
  unpaired population absorption needs. **Prediction:** sweeping
  $`N_0/|E|`$ should drive $`f`$ to $`1/2`$ at some finite ratio, at which the
  ledger closes exactly and the free constant of §7 is determined rather than
  chosen. If instead $`f`$ saturates below $`1/2`$, the ledger has an
  irreducible leak and the free constant stays free. This is the decisive
  next measurement and it is cheap.
- **S-SP2 (tension with the sea-dressed note).**
  [`sea_dressed_microdynamics.md`](sea_dressed_microdynamics.md) §8 and §10
  show a live ledger converging onto the pinned run as $`\kappa_{\rm rec}`$
  grows. S5 finds no such convergence — the error saturates at 40 per cent.
  The two are not a like-for-like comparison: that note runs a ring with a
  cosine well and the uncompensated stencil, this one the open line with the
  compensated split at reach $`2\pi a`$. Which difference is responsible is
  unmeasured, and until it is, neither result should be read as overturning
  the other.
- **S-SP3 (what $`f`$ depends on).** $`f`$ rises with reach — 0.369, 0.434,
  0.462 at $`y_{\max}/a = \pi, 2\pi, 4\pi`$ — presumably because a longer
  reach opens more channels and so more chances of finding a partner. Whether
  it approaches $`1/2`$ from below in the limit, or crosses it, is not known,
  and the residual deficit still tracks the reach.
- **S-SP4 (the mean-field caps).** The absorptive allocation is a tau-leap
  processed channel by channel with live population caps, in the manner of
  `sea_dressed_lattice.step_jump_sea_mc`. Sequential processing makes later
  channels see earlier channels' depositions; that is an $`O(\Delta t)`$
  effect and it is inside the measured first-order convergence, but the
  allocation is not symmetric in the channel ordering and no one has checked
  how much $`f`$ moves under a permuted order.
- **S-SP5 (a genuine ensemble).** Everything here is mean field on a mesh.
  The claimed advantage is a world-ensemble claim, so CLA2's decisive test —
  a signed ensemble against the mesh at fixed variance — still stands, now
  with a specific prediction attached: the absorptive ensemble should show the
  emissive ensemble's variance reduced by roughly $`N_{\rm em}/N_{\rm abs}`$.
- **S-SP6 (S4 on a bound state).** The runs here are transient. A stationary
  state would make "equilibrium" exact rather than asymptotic, but the Eckart
  barrier has no bound states; K-LS5's soft-core Coulomb, which has both a
  tunable ceiling and an attractive well, is the natural vehicle.

---

## 9. Numerical verification

`src/demo_sea_population_equilibrium.py`, run as
`WPMW_OUTPUT=... PYTHONPATH=src python3 src/demo_sea_population_equilibrium.py`.

| Part | Claim verified |
|---|---|
| A | S1 stated; S2 closed form to $`\sim 10^{-13}`$ from both directions, five $`\kappa`$ |
| B | S3, the $`2/(\kappa t)`$ asymptote; the packet peak against $`2/h`$; three preparations |
| C | S4, the unbounded worst-cell drain at conserved mean |
| D | S6, the body-momentum table for the three realisations |
| E | S7, both ledger identities to $`1.5\times10^{-14}`$ and $`1.7\times10^{-11}`$ |
| F | S8, the fidelity and deficit tables, first-order in $`\Delta t`$ |
| G | reach dependence of $`f`$ and of $`N_{\rm eq}`$; the $`\kappa \propto \Gamma`$ scaling |

Figures on the `output` branch:

[![The ledger fixed point](https://raw.githubusercontent.com/billpage/wpmw/output/figures/sea_population_fixed_point.png)](https://raw.githubusercontent.com/billpage/wpmw/output/figures/sea_population_fixed_point.png)

$`N_{\rm eq}`$ against $`\kappa`$ with its two limits; the universal $`1/t`$
vacuum asymptote from three preparations; and $`\Gamma_{\rm tot}(r)`$ on the
barrier with the three K7 quiet points marked.

[![Emissive against absorptive](https://raw.githubusercontent.com/billpage/wpmw/output/figures/sea_population_unravelling.png)](https://raw.githubusercontent.com/billpage/wpmw/output/figures/sea_population_unravelling.png)

The body count and the worst-cell sea occupancy for the two unravellings, and
the fidelity of each against the exact mesh QLE as $`\Delta t`$ falls.

[![The ledger identity](https://raw.githubusercontent.com/billpage/wpmw/output/figures/sea_population_ledger_identity.png)](https://raw.githubusercontent.com/billpage/wpmw/output/figures/sea_population_ledger_identity.png)

The absorptive fraction against the reach with the closing value $`f = 1/2`$
marked, and the residual deficit alongside it.

---

## 10. Sources

- [`../algorithm/compensated_liouville_algorithm.md`](../algorithm/compensated_liouville_algorithm.md)
  §4.3 — the unravelling choice this note prices, and open item CLA3.
- [`eckart_barrier_compensated.md`](eckart_barrier_compensated.md) §8 —
  world-particle identity, the emission rate, and Theorem K7's quiet points.
- [`species_sectors_and_annihilation.md`](species_sectors_and_annihilation.md)
  — Proposition U1, Corollary D5.1, and the §11 design caution on $`\kappa`$.
- [`sea_dressed_microdynamics.md`](sea_dressed_microdynamics.md) §8 — the live
  ledger and the orphan avalanche; see S-SP2.
- [`../supplement/representation_cost_and_annihilation.md`](../supplement/representation_cost_and_annihilation.md)
  §8 — open item N1, which asked whether the action set contains its own
  garbage collection. S7 is the nearest thing this note has to an answer:
  it does, conditionally on $`f = 1/2`$.
- The absorptive realisation was proposed by B. Page (conversation, August
  2026), from the observation that a negaton arriving in a region of positon
  excess can be realised as a positon hole.
