# The stochastic ledger: a jump process, a sum rule, and a standing population

> The ledger run as an exact Markov jump process on integer counts, rather than as a mean field on a mesh. Answers the standing request for a formal stochastic treatment of the free-body and sea-pair populations, and corrects the law it was meant to confirm. Theorem N1 strengthens J1: `P = S + N/2` is conserved *pathwise*, in exact integers, on every trajectory and not merely in expectation — a deterministic invariant of the generator, not a martingale. Theorem N2 is why the ledger is invisible: both realisations of an event move `E` by `+1` at the upper daughter and `-1` at the lower, so the Bernoulli(`f`) choice between them lies entirely in the kernel of the observable map, and the only noise `W` ever sees is Poisson event-timing noise common to both branches. The framing consequence is worth stating plainly — this is not a stochastic mechanics in Nelson's sense: trajectories are exactly Newtonian by (S), the diffusion in the phase variables is identically zero, and the whole stochastic content is demographic. **Theorem N3 is the note's centre and a correction to S7.** Stationarity of the body count gives `Gamma_tot (1 - 2f) = R_sink` with no closure and nothing about the potential, where `R_sink` is the total rate of every *other* body-removing channel — so `f = 1/2` is the sinkless special case and not the law, and any `kappa > 0` forces `f < 1/2` by a computable amount, verified to between 0.01 and 0.9 per cent over a fortyfold range of `kappa`. That invalidates the argument S-SP3 used to rule out a floor — that a systematic leak cannot change sign — because S-SP7 identifies a sink *and* a source, and their difference changes sign freely. Theorems N4 and N5 close the availability question: `f` is the probability that both legs find a partner, which under an independent-occupancy closure is `(1 - e^-lam)^2`, so `f = 1/2` pins `lam* = -ln(1 - 2^-1/2) = 1.227947` bodies per species per cell, `f'(lam*) = sqrt2 - 1` exactly, `Var(Lambda)/M = 1/f' = 2.414214`, and a Fano factor `0.983028`, slightly sub-Poissonian — five numbers carrying no free constant, measured at 2.4487 against 2.4559 in the well-mixed limit and flat across `Q = 1` to `12`, as a closure with no channel index in it requires. In world-particles that is 2.456 bodies against exactly two sea pairs per Planck cell, the low end of K9's empirical window. Theorem N6 reassigns a role the project had given to recombination: with transport off, per-cell occupancy is a reflected critical random walk whose spread grows without bound, `kappa` damps it only partly and drags `f` to 0.32, and turning streaming on instead holds the spread flat and leaves `f` at 0.5005. Streaming, not recombination, is the local regulator — and by N3 no sink could have been, since every sink moves `f` off one half.

*Ladder abstract — see the [full list](README.md#the-ladder).*

**Status.** Analysis note, step 18 of the ladder. Companion demo:
`src/demo_stochastic_ledger.py`. Prompted by the standing request for an
analytic — stochastic-differential or otherwise formal — account of the free
world-particle population and the long-term sea-pair equilibrium left open by
[`sea_population_equilibrium.md`](sea_population_equilibrium.md) §8 and
[`../supplement/emission_and_absorption.md`](../supplement/emission_and_absorption.md)
§12.

---

## 0. What this note settles, and what it corrects

Everything the project has said about the ledger so far has been said about a
mean field on a mesh: $`u^\pm(x, p)`$ and $`S(x, p)`$ as real-valued
densities, updated by a tau-leap. That is the right tool for fidelity, and the
wrong one for the question actually being asked. "Does the population reach
equilibrium, and at what level" is a question about a *counting* process —
integers, births, deaths, and the fluctuations around the balance — and a mean
field answers it only by assuming the answer.

This note runs the ledger as what it is: a continuous-time Markov jump process
on integer counts per cell. The generator is written down, its invariants are
identified exactly, its diffusion limit is taken, and the resulting
Ornstein–Uhlenbeck law is checked against exact Gillespie trajectories.

**Settles.**

- **The formal object.** The ledger is a density-dependent Markov jump process
  in Kurtz's sense, with a deterministic law of large numbers and a Gaussian
  fluctuation theorem (§6). No new mathematics is required, and the literature
  to cite is population-process literature, not quantum literature.
- **The standing population.** $`\lambda^* = 1.227947`$ bodies per species per
  cell, with no free constant (§5). This is the self-consistency condition
  S-SP1 said would be needed to fix $`\kappa`$.
- **What regulates the local ledger.** Transport, not recombination (§7).

**Corrects.**

1. [`sea_population_equilibrium.md`](sea_population_equilibrium.md) Theorem S7
   states $`f = 1/2`$ as *the* closure condition. It is the closure condition
   of a ledger with no other sink. The general law is Theorem N3, and an
   erratum has been added to that note.
2. The same note's open item **S-SP3** argues that the measured shortfall
   cannot be a floor because a systematic leak cannot change sign. That
   argument does not hold: S-SP7 identifies two channels of opposite sign — a
   sink (the negative caps) and a source (the population clamp) — and their
   difference changes sign freely. §4.3 says what measurement would settle it
   instead.
3. The informal reading, used in discussion and implicit in §7 of the sea
   note's constraint table, that a recombination rate $`\kappa`$ is what holds
   the local population steady. Theorem N6 shows it cannot be: by N3 every
   sink moves $`f`$ off one half, and $`\kappa`$ at the strength needed to
   damp the local spread costs about a third of the absorptive fraction.

**Does not settle.** The value of $`\kappa`$. §5 supplies the number
$`\kappa`$ would have to reproduce, which is progress, but turning
$`N_{\rm eq}(\kappa) = 2\lambda^*`$ into a determination requires matching a
continuum mean field to a lattice occupancy, and that matching is not done
here. Logged as **N-SP2**.

---

## 1. The process

Take the momentum axis as a ring of $`M`$ cells, so that the geometry of the
event — parent at $`c`$, daughters at $`c \pm q`$ — is preserved and nothing
else is. Per cell the state is three integers,

```math
n^+_c, \qquad n^-_c, \qquad S_c ,
```

the positon count, the negaton count, and the bound sea-pair count. Write
$`N_c = n^+_c + n^-_c`$ for the body count, $`E_c = n^+_c - n^-_c`$ for the
observable, and $`\Lambda = \sum_c N_c`$ for the total.

Three channels.

**Event.** For each of the $`MQ`$ pairs $`(c, q)`$ with
$`1 \le q \le Q`$, at rate $`\gamma`$. Write $`A = c + q`$ and $`B = c - q`$.
The residual operator deposits $`+1`$ at $`A`$ and $`-1`$ at $`B`$, and there
are exactly two ways to realise that (§6 of the sea note):

```math
\text{emissive:}\quad
n^+_A \!+\!\!=\! 1, \quad n^-_B \!+\!\!=\! 1, \quad S_c \!-\!\!=\! 1
```

```math
\text{absorptive:}\quad
n^-_A \!-\!\!=\! 1, \quad n^+_B \!-\!\!=\! 1, \quad S_c \!+\!\!=\! 1
```

The selection rule is local and blind, exactly as in §7 of
[`../supplement/emission_and_absorption.md`](../supplement/emission_and_absorption.md):
settle absorptively iff $`n^-_A \ge 1`$ **and** $`n^+_B \ge 1`$; otherwise
emissively, if the parent cell has a sea pair to ionise. Nothing measures
$`f`$ and nothing targets it.

**Hop.** Each body moves to another cell at rate $`\nu`$, destination uniform
on the ring. This is the caricature of postulate (S). It carries no force law
and no geometry — the point of §7 is that its *presence* is what matters, not
its details — and uniform destinations are the well-mixed idealisation under
which the closure of §5 is supposed to hold.

**Recombination.** At each cell, at rate $`\kappa\, n^+_c n^-_c`$:
$`n^+_c\!-\!\!=\!1`$, $`n^-_c\!-\!\!=\!1`$, $`S_c\!+\!\!=\!1`$. This is the
bilinear coincident-pair sink Theorem S1 forces, and it returns a bound pair
rather than destroying two bodies, so that the pair count is preserved.

That is the whole generator. It is simulated exactly — Gillespie (1977), no
time discretisation and no tau-leap — so nothing below is a stepping artefact.

---

## 2. Theorem N1: the pair count is a pathwise invariant

**Theorem N1.** $`P = \sum_c \bigl(S_c + N_c/2\bigr)`$ takes the same value on
every state visited by every trajectory of the process.

*Proof.* Emissive: $`\Delta S = -1`$, $`\Delta N = +2`$. Absorptive:
$`\Delta S = +1`$, $`\Delta N = -2`$. Recombination: $`\Delta S = +1`$,
$`\Delta N = -2`$. Hopping moves a body between cells and changes neither sum.
Each gives $`\Delta P = 0`$. $`\square`$

This is stronger than Theorem J1, which is the same statement about
expectations of a mean field. Here $`2P`$ is an integer and it does not move
at all — not to $`10^{-16}`$, not on average, not eventually. Part A measures
$`\max|2P(t) - 2P(0)| = 0`$ over $`3\times10^{5}`$ steps with all three
channels live.

The distinction matters because a conserved *expectation* permits the variance
to grow, and a pathwise invariant does not. Any ensemble the algorithm is ever
asked to carry lies on a fixed level set of $`P`$, and that level set is set
once, by the initial condition.

Theorem J2 survives unchanged and is worth restating in this language: $`P`$
is conserved globally and transported locally, so the process has two
independent fields per cell and one global constraint. §7 is what that costs.

---

## 3. Theorem N2: the realisation is a null direction

**Theorem N2.** Let $`v_{\rm em}`$ and $`v_{\rm abs}`$ be the state increments
of the two realisations of the same event $`(c, q)`$, and let $`\Pi`$ be the
map to the observable, $`\Pi: (n^+, n^-, S) \mapsto E = n^+ - n^-`$. Then

```math
\Pi\, v_{\rm em} \;=\; \Pi\, v_{\rm abs}
\;=\; \bigl(+1 \text{ at } A,\; -1 \text{ at } B\bigr),
\qquad\text{so}\qquad
\Pi\,\bigl(v_{\rm em} - v_{\rm abs}\bigr) \;=\; 0 .
```

*Proof.* Emissive adds a positon at $`A`$, so $`\Delta E_A = +1`$; absorptive
removes a negaton at $`A`$, so $`\Delta E_A = +1`$. At $`B`$: emissive adds a
negaton, absorptive removes a positon, both $`\Delta E_B = -1`$. $`\square`$

Two consequences, and the second is the one that matters.

**The observable cannot see the realisation.** Whatever the process does with
$`f`$ — regulate it, fail to regulate it, let it wander — $`E`$ evolves
identically. This is the exact statement of which the sea note's "identical in
the observable and opposite in the ledger" is the informal version.

**The Bernoulli noise has no observable component.** In the diffusion limit of
§6 the event contributes two independent sources of randomness: *when* events
fire, which is Poisson and common to both branches, and *which branch* fires,
which is Bernoulli$`(f)`$. By N2 the second has zero projection onto $`E`$.
So the whole ledger carries a noise channel that is invisible in $`W`$, and
the only noise the observable inherits is event-timing noise.

Part A verifies this by writing the predicted $`\Delta E`$ before the branch is
chosen and comparing: maximum per-cell discrepancy 0, maximum
$`|\sum_c \Delta E|`$ on an event 0, over $`3\times10^{5}`$ steps.

The practical corollary is for **S-SP5**, which predicts that an absorptive
ensemble shows the emissive ensemble's variance reduced by roughly
$`N_{\rm em}/N_{\rm abs}`$. N2 says why: the per-event noise in $`E`$ is the
same in both, so the whole difference is the population over which the signed
sum is taken. That is an argument for the prediction, not merely a hope.

---

## 4. Theorem N3: the sum rule, and what it corrects

### 4.1 The rule

**Theorem N3.** In any statistically stationary state of the process,

```math
\Gamma_{\rm tot}\,\bigl(1 - 2f\bigr) \;=\; R_{\rm sink} ,
```

where $`\Gamma_{\rm tot}`$ is the total event rate, $`f`$ the absorptive
fraction, and $`R_{\rm sink}`$ the total rate of every body-removing channel
*other than* absorption. The same equation follows from stationarity of the
sea.

*Proof.* Each emissive event changes $`\Lambda`$ by $`+2`$, each absorptive
event by $`-2`$, each removal of a coincident pair by $`-2`$. Hopping does not
change $`\Lambda`$. So

```math
\frac{d\langle \Lambda\rangle}{dt}
= 2\Gamma_{\rm tot}(1-f) - 2\Gamma_{\rm tot}f - 2R_{\rm sink}
= 2\Gamma_{\rm tot}(1-2f) - 2R_{\rm sink},
```

which vanishes iff the stated equation holds. For the sea, each emissive event
gives $`\Delta S = -1`$ and each absorptive event and each removal
$`\Delta S = +1`$, so
$`d\langle S_{\rm tot}\rangle/dt = \Gamma_{\rm tot}(2f-1) + R_{\rm sink}`$,
the negative of half the first. $`\square`$

Nothing in the proof mentions the potential, the reach, the closure, the
availability rule, or the lattice. It is arithmetic on $`\pm 2`$, in the same
spirit as §6 of the emission-and-absorption note and with the same robustness.

**Theorem S7 is the case $`R_{\rm sink} = 0`$.** That is the case the sea
note's Part H and Part J actually run — `run_traced` carries no recombination
— which is why S7 came out clean there. But the law is the displayed equation,
and S7's $`f = 1/2`$ should be read as a corollary with a hypothesis attached.

### 4.2 Three things that follow

**Every sink costs absorptive fraction.** $`\kappa > 0`$ forces $`f < 1/2`$
strictly, by $`\tfrac12 - f = R_{\rm sink}/2\Gamma_{\rm tot}`$. Part B, at
$`M = 32`$, $`Q = 2`$, $`\gamma = 1`$, $`\nu = 8`$:

```
   kappa        f      1/2 - f    R_sink/(2 Gamma_tot)    rel err    mean/M
    0.00   0.49989   +0.00011          +0.00000      0.0001   2.5612
    0.05   0.48278   +0.01722          +0.01723      0.0005   2.4397
    0.20   0.44035   +0.05965          +0.05968      0.0006   2.2311
    0.50   0.38533   +0.11467          +0.11443      0.0021   1.9722
    1.00   0.32882   +0.17118          +0.16970      0.0087   1.7364
    2.00   0.26664   +0.23336          +0.23327      0.0004   1.4689
```

Between 0.01 and 0.9 per cent over a fortyfold range, with the invariant of N1
exact in every row.

**A measured shortfall is a measurement.** If the code has an unaccounted
sink, the shortfall reports its rate. That is a diagnostic, not an error term:
instrument the sink and the identity closes or it does not.

**A measured *excess* reports a source.** $`f > 1/2`$ means
$`R_{\rm sink} < 0`$, which is a channel manufacturing bodies.

![The sum rule, the closure, and the local regulator](https://raw.githubusercontent.com/billpage/wpmw/output/figures/stochastic_ledger_summary.png)

(a) Theorem N3 over a fortyfold range of the recombination rate — points
measured, line the independently measured sink rate; (b) Theorem N5's closure
against a scan in $`M`$; (c) Theorem N6, the per-cell spread against time
under transport and under recombination.

### 4.3 What this does to S-SP3

Open item **S-SP3** of the sea note reads, of the reach ladder's shortfall:

> It changes sign, at $`\rho = 20`$ in Part H and at
> $`(\Delta p, y_h) = (0.0625, 2\pi)`$ here. A systematic per-event ledger
> leak cannot change sign, so what is left is truncation error, not a floor.

The premise is right and the inference is not. A single systematic leak cannot
change sign. But **S-SP7**, in the same list, identifies two channels of
opposite sign in the same code: the negative caps, which run events backwards
and act as a sink, and the population clamp, which "manufactures bodies" and
acts as a source. By N3 the shortfall reports their *difference*, and a
difference changes sign as soon as the two cross. The sign change is therefore
consistent with a floor as well as with truncation, and cannot discriminate.

The conclusion of S-SP3 — that $`f \to 1/2`$ as the regulator is removed — is
not challenged here; the monotone reach ladder is separate evidence and it
stands. What is withdrawn is one of the two arguments offered for it.

**The measurement that would settle it.** N3 is an identity, so instrument it.
Accumulate, per step of the mesh demo, the total body count destroyed by
clipped caps and created by the clamp, form $`R_{\rm sink}`$, and check
$`\Gamma_{\rm tot}(1-2f) = R_{\rm sink}`$ directly. If it closes, the
shortfall is fully accounted for by identified channels and S-SP3 is settled
exactly rather than by extrapolation from four points with a wobbling tail.
Logged as **N-SP1**.

---

## 5. Theorems N4 and N5: availability, and a standing population

### 5.1 The closure

$`f`$ is not a parameter. It is the probability that an event finds a partner
at *both* daughter cells — the conjunction cost of §8 of the
emission-and-absorption note, now stated as a probability rather than as a
ratio of means. If per-cell per-species occupancy is Poisson with mean
$`\lambda`$ and the two legs are independent,

```math
f(\lambda) \;=\; \bigl(1 - e^{-\lambda}\bigr)^{2} .
```

**Theorem N4 (no channel index).** $`f`$ as written contains no $`q`$ and no
$`Q`$: availability is a question about a cell, and which channel is asking
does not enter. So the standing population is independent of the number of
channels.

Part D, $`M = 32`$, $`\gamma = 1`$, $`\kappa = 0`$, mixing scaled with $`Q`$:

```
       Q        f      mean/M    Var/M    Fano
       1   0.49987   2.4836   2.3169  0.9329
       2   0.49994   2.4780   2.5749  1.0391
       4   0.49998   2.4756   2.7424  1.1078
       8   0.49995   2.4701   2.6677  1.0800
      12   0.49997   2.4648   2.3836  0.9670
```

Flat to about half a per cent across a twelvefold range, and drifting *down*
towards the predicted 2.4559 as $`Q`$ rises.

### 5.2 Five numbers

**Theorem N5.** Under the closure, with no other sink, the stationary ledger
is fixed with no free constant:

```math
f(\lambda^*) = \tfrac12
\;\Longrightarrow\;
\lambda^* = -\ln\!\bigl(1 - 2^{-1/2}\bigr) = 1.227947 ,
```

```math
f'(\lambda^*) \;=\; 2\bigl(1-e^{-\lambda^*}\bigr)e^{-\lambda^*}
\;=\; 2\cdot\tfrac{1}{\sqrt2}\Bigl(1 - \tfrac{1}{\sqrt2}\Bigr)
\;=\; \sqrt{2} - 1 \;=\; 0.414214 ,
```

and, from the diffusion limit of §6,

```math
\frac{\mathrm{Var}(\Lambda)}{M} = \frac{1}{f'(\lambda^*)}
= \sqrt2 + 1 = 2.414214 ,
\qquad
\mathrm{Fano} = \frac{\mathrm{Var}(\Lambda)}{\langle\Lambda\rangle}
= \frac{1}{2 f'(\lambda^*)\lambda^*} = 0.983028 .
```

*Proof of the two closed forms.* $`1 - e^{-\lambda^*} = 2^{-1/2}`$ by
definition of $`\lambda^*`$, so $`e^{-\lambda^*} = 1 - 2^{-1/2}`$ and the
product in $`f'`$ telescopes to $`\sqrt2 - 1`$. $`\square`$

The Fano factor is the one to watch: it is *not* 1, so the regulated
population is slightly sub-Poissonian. That is a falsifiable prediction about a
quantity nothing was tuned to produce.

Part E, $`\gamma = 1`$, $`Q = 2`$, $`\nu = 8`$, $`\kappa = 0`$:

```
       M        f      mean/M    Var/M    Fano
       8   0.49892   2.5229   2.8090  1.1134
      16   0.49954   2.4825   2.6252  1.0575
      32   0.49977   2.5246   2.7384  1.0847
      64   0.49971   2.5442   2.5788  1.0136
     128   0.49995   2.5019   2.5453  1.0173
```

$`\mathrm{Var}/M`$ flat is $`\mathrm{Var}(\Lambda)`$ linear in $`M`$, which is
what §6 predicts. The best measurement of the mean is Part C's well-mixed end,
$`\gamma/\nu = 1/32`$: **2.4487 against 2.4559 predicted, 0.3 per cent.**

### 5.3 What the number means

$`2\lambda^* = 2.455894`$ bodies per cell, against a sea of exactly two pairs
per cell — the crystal shift $`B = 2/h`$ of §2 of the sea note. So the
theory's standing unpaired population is

```math
\text{2.456 world-particles per Planck cell, against 2 sea pairs,}
```

with no adjustable constant anywhere in it. Theorem K9 of
[`eckart_barrier_compensated.md`](eckart_barrier_compensated.md) reports an
empirically usable window of $`\rho = 3`$ to $`10`$ against the same two sea
pairs. The derived value sits just below that window, which is the right
relationship if K9's lower end is set by the same availability requirement
with a safety margin, and is a discrepancy to explain if it is not.

This is also the shape of answer **S-SP1** asked for. That item closes by
saying the free constant now "needs a different argument entirely — a
candidate is the requirement that $`N_{\rm eq}`$ equal the standing population
S9's feedback actually sustains, which would be a self-consistency condition."
§5 supplies the standing population. What it does not supply is the dictionary
between a continuum $`N_{\rm eq}`$ in Wigner units and a lattice occupancy in
bodies per cell. **N-SP2.**

![The ledger relaxing onto the derived standing population](https://raw.githubusercontent.com/billpage/wpmw/output/figures/stochastic_ledger_relaxation.png)

A single exact trajectory started at four times the equilibrium population,
against the closure's $`2\lambda^*`$. Nothing in the run is told where to go.

---

## 6. The diffusion limit, and what it is not

### 6.1 The limit

The process is density-dependent in Kurtz's sense: rates depend on the state
through the per-cell densities. Its Kramers–Moyal expansion truncates at
second order in the usual way (van Kampen ch. X), giving for the total body
count

```math
d\Lambda \;=\; 2\Gamma_{\rm tot}\bigl(1 - 2f(\lambda)\bigr)\,dt
\;+\; \sqrt{4\Gamma_{\rm tot}}\;dW ,
\qquad \lambda = \frac{\Lambda}{2M} ,
```

the drift being exact (the jump is $`\pm2`$) and the diffusion coefficient
being the second jump moment times the rate, $`4\Gamma_{\rm tot}`$, evaluated
at $`f = 1/2`$ where $`(1-2f)^2`$ drops out. Linearising about $`\lambda^*`$,

```math
\frac{\partial}{\partial\Lambda}\Bigl[2\Gamma_{\rm tot}(1-2f)\Bigr]
= -\frac{2\Gamma_{\rm tot}f'(\lambda^*)}{M}
\;\equiv\; -k ,
```

an Ornstein–Uhlenbeck process with relaxation rate $`k`$ and stationary
variance

```math
\mathrm{Var}(\Lambda) \;=\; \frac{4\Gamma_{\rm tot}}{2k}
\;=\; \frac{M}{f'(\lambda^*)} .
```

**$`\Gamma_{\rm tot}`$ cancels.** The stationary ledger fluctuation is
independent of the event rate, hence of the coherence reach, hence of the
potential and of everything else the emission rate depends on. Only the
availability function survives. That is a strong statement and it is the
reason §5's numbers can be quoted as constants at all.

The relaxation *rate* does not cancel: $`k = 2\Gamma_{\rm tot}f'/M
= 2\gamma Q f'(\lambda^*)`$, independent of $`M`$ and proportional to the
per-cell event rate. So a ledger far from equilibrium returns on the event
clock, not the wall clock — which is exactly what §7 of the
emission-and-absorption note measured, four regimes differing nineteenfold in
wall-clock time collapsing onto one curve when clocked against events. That
measurement is now a consequence rather than an observation.

### 6.2 What this is not

It is worth being explicit, because the phrase "stochastic differential
equation for quantum mechanics" points at a large and different literature.

In Nelson's stochastic mechanics (Nelson 1966, 1985) the noise is in the
trajectory: $`dx = b\,dt + \sqrt{\hbar/m}\,dW`$, and a particle's path is
nowhere differentiable. In the compensated ontology the opposite holds.
Postulate (S) makes every world-particle's path an exact Newtonian arc under
the full classical force, so the diffusion coefficient in $`(x, p)`$ is
**identically zero**; §8.1 of the Eckart note says this in words and N2 says
it again — the residual channel acts only by creating and destroying, never by
moving.

So the stochastic differential equation of this theory is not a diffusion in
phase space at all. It is a jump process in the *counting* variables, with

- zero diffusion in $`(x, p)`$,
- Poisson noise in the event count, which the observable inherits,
- Bernoulli noise in the realisation, which by N2 the observable does not.

The right comparison class is birth–death and reaction–diffusion population
processes, not Nelson diffusions or Bohmian trajectories with added noise. The
one place the two literatures touch is §6.1's collapse onto the event clock,
which is a standard feature of density-dependent jump processes and has no
counterpart in a diffusion.

---

## 7. Theorem N6: transport is the local regulator

N1 and J2 leave the ledger globally constrained and locally free. What happens
locally is worth stating as a result, because the project has been assuming
the wrong answer.

Consider a cell with transport switched off. Bodies of one species arrive by
emission at rate $`\gamma_{\rm cell}(1-f)`$ and leave by absorption at rate
$`\gamma_{\rm cell}f`$, gated below at zero. At the global fixed point
$`f = 1/2`$ those rates are equal: the per-cell occupancy performs a
**critical random walk reflected at the origin**, which on an unbounded state
space has no stationary distribution. Its spread grows without bound, the
distribution across cells becomes bimodal — many empty cells and a few rich
ones — availability collapses, and the *global* $`f`$ falls away from one half
even though nothing is wrong with the global argument.

**Theorem N6.** Recombination does not repair this and cannot. A finite
$`\kappa`$ damps the local spread only partially, and by N3 it buys that
damping at a cost of exactly $`R_{\rm sink}/2\Gamma_{\rm tot}`$ in absorptive
fraction. Transport repairs it at no cost in $`f`$ at all, because hopping is
not a sink.

Part F, variance of the per-cell body count across $`M = 32`$ cells:

```
     nu   kappa      t=25      t=50     t=100     t=200     t=400       f
    0.0    0.00     80.50    126.90    624.80   1465.62   3952.18    0.4837
    0.0    0.50     28.15     55.15    211.15    186.84    303.00    0.3197
    0.0    2.00     32.50     44.46    169.75    189.06    367.12    0.3007
    2.0    0.00      2.68      2.37      3.55      2.56      3.23    0.5004
    8.0    0.00      2.23      2.31      2.93      3.25      2.81    0.5005
```

Transport off and no sink: the spread grows by a factor fifty over the run.
Transport off with $`\kappa = 0.5`$: the spread is held to a few hundred and
$`f`$ is dragged to 0.32. Transport on at the mildest setting tried: the
spread is flat at about 3 and $`f`$ sits at 0.5004.

Two readings of this, and both are worth having.

**Ontologically**, it says postulate (S) is doing work the ontology note did
not credit it with. (S) is stated as the claim that world-particles have
Newtonian histories — an identity claim. It is also the mechanism that makes
the local ledger a stationary object, and without it the demographic channel's
own thermostat fails, not because the thermostat is wrong but because it is
global and the supply it regulates is local.

**Practically**, it explains the size of the improvement from the `stream()`
repair recorded in the sea note's erratum. That bug left bodies unable to move
in $`p`$, which is exactly the $`\nu = 0`$ column. The repair moved the note's
floor at $`\Delta t = 0.01`$ from $`-1.2353`$ to $`-0.3386`$ and turned the
absorptive trace positive. N6 predicts an improvement of that character before
the fact.

---

## 8. Numerical verification

`src/demo_stochastic_ledger.py`, exact Gillespie throughout.

| part | claim | result |
|---|---|---|
| A | N1, pathwise invariant | $`\max|2P(t)-2P(0)| = 0`$, exact integer |
| A | N2, null direction | per-cell $`|\Delta E| = 1`$, $`|\sum_c\Delta E| = 0`$ |
| B | N3, the sum rule | rel. err. $`10^{-4}`$ to $`9\times10^{-3}`$, $`\kappa \in [0.05, 2]`$ |
| C | $`f \to 1/2`$ | 0.49977 to 0.50013 over $`\gamma/\nu \in [1/32, 1]`$ |
| C | closure, mean | 2.4487 against 2.4559, well-mixed end |
| D | N4, no $`Q`$ | 2.4836 to 2.4648 over $`Q = 1`$ to $`12`$ |
| E | N5, Var linear in $`M`$ | $`\mathrm{Var}/M`$ = 2.55 to 2.81 over $`M = 8`$ to $`128`$ |
| E | N5, Fano | 1.01 to 1.11 against 0.9830 |
| F | N6, local regulator | spread $`\times 50`$ at $`\nu=0`$, flat at $`\nu = 2`$ |

Two honest caveats. The Fano factor is measured at 1.01–1.11 against a
prediction of 0.9830; the runs resolve the prediction's *order* but not the
2 per cent by which it differs from unity, so the sub-Poissonian claim is at
present a prediction and not a measurement. And the closure's mean is
approached from above, by 0.3 per cent at the well-mixed end and by 2 per cent
at $`\gamma = \nu`$, consistent with a small positive correlation between the
two legs' occupancies that the independence assumption neglects — which is the
same quantity **J-SP1** asks about, now with a sign attached.

---

## 9. Open items

- **N-SP1 (instrument the sum rule in the mesh demo).** N3 is an identity, so
  the shortfall in `demo_sea_population_equilibrium.py` can be attributed
  rather than extrapolated. Accumulate the body count destroyed by clipped
  caps and created by the population clamp, form $`R_{\rm sink}`$, and test
  $`\Gamma_{\rm tot}(1-2f) = R_{\rm sink}`$ per step. Closing it settles
  S-SP3 exactly; failing to close it locates a third channel nobody has
  named. This is the highest-value item on the list.
- **N-SP2 (the dictionary, and hence $`\kappa`$).** $`\lambda^*`$ is an
  occupancy in bodies per cell; $`N_{\rm eq}(\kappa)`$ of Theorem S2 is a
  density in Wigner units. Matching them would turn §5 into a determination
  of $`\kappa`$. The obstruction is that the mesh has no cell size in the
  relevant sense — the lattice's $`\Delta p`$ is set by the reach and its
  $`\Delta x`$ by the grid — so the conversion is not merely a change of
  units and may not exist.
- **N-SP3 (measure the Fano factor properly).** 0.9830 against 1 is a 1.7 per
  cent effect and the present runs cannot see it. It needs either far longer
  runs or a variance-reduction estimator; it is worth the effort because
  sub-Poissonian statistics are a signature of feedback and would distinguish
  the availability mechanism from any tuned alternative.
- **N-SP4 (transport that is transport).** The hop channel here is uniform
  on a ring, which is the well-mixed idealisation. Real streaming is
  Newtonian and local, so it mixes a cell with its neighbours on a timescale
  set by the classical flow, not instantly. Where the reaction length beats
  the transport length, N6's repair will be partial and the local ledger will
  be somewhere between the two columns. This is the same measurement
  **K-LS7** wants — whether the sea inherits the emission lobes — approached
  from the ledger side, and the two should be done together.
- **N-SP5 (the two legs are not independent).** §8's measurement is 0.3 to 2
  per cent above the independence closure, systematically and from one side.
  A joint occupancy measurement would turn J-SP1 from an open question into a
  correction term, and the correction has a known sign.

---

## 10. Sources

- T. G. Kurtz, *Solutions of ordinary differential equations as limits of pure
  jump Markov processes*, J. Appl. Prob. **7** (1970) 49–58; and *Limit
  theorems for sequences of jump Markov processes approximating ordinary
  differential processes*, J. Appl. Prob. **8** (1971) 344–356. The law of
  large numbers and the Gaussian fluctuation theorem for §6.
- S. N. Ethier and T. G. Kurtz, *Markov Processes: Characterization and
  Convergence*, Wiley 1986, ch. 11. The textbook statement of both, with the
  hypotheses in usable form. The one to check here is the Lipschitz condition
  on the rate function, which the availability indicator $`n \ge 1`$ fails
  and which the closure of §5.1 is what repairs — see the annotation in
  [`../../references/bibliography.md`](/billpage/wpmw/blob/main/references/bibliography.md#density-dependent-markov-jump-processes-and-their-diffusion-limits).
- N. G. van Kampen, *Stochastic Processes in Physics and Chemistry*,
  3rd ed., North-Holland 2007, ch. X and ch. XI. The system-size expansion,
  and the critical birth–death walk behind §7.
- D. T. Gillespie, *Exact stochastic simulation of coupled chemical
  reactions*, J. Phys. Chem. **81** (1977) 2340–2361. The algorithm used
  throughout the companion demo.
- C. W. Gardiner, *Handbook of Stochastic Methods*, 4th ed., Springer 2009,
  ch. 11. Ornstein–Uhlenbeck stationary statistics.
- E. Nelson, *Derivation of the Schrödinger equation from Newtonian
  mechanics*, Phys. Rev. **150** (1966) 1079–1085; *Quantum Fluctuations*,
  Princeton 1985. The contrast drawn in §6.2.
- J. M. Sellier, *A signed particle formulation of non-relativistic quantum
  mechanics*, J. Comput. Phys. **297** (2015) 254–265. Already in
  [`../../references/bibliography.md`](../../references/bibliography.md);
  the closest existing signed-ensemble formulation, and the one whose
  emission rate Theorem G2 compares against.
