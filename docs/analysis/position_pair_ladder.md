# The position-space coherence ladder

**Status.** Analysis note, step 9 of the ladder. Companion demo:
`src/demo_position_pair_ladder.py`.

## 0. What this note inherits, retracts and corrects

**Inherits.** From
[permanent_pairing_density_matrix.md](permanent_pairing_density_matrix.md):
the pair-as-element ontology, the species-is-ket-or-bra identification, and
the four-term one-leg-hop reading of the commutator. From
[coherence_ladder.md](coherence_ladder.md): the rung index, the
channel-table method, and Theorem C2's standard of proof — expected
elementwise flux equal to the commutator on every element of every state.
From [phase_alignment_microdynamics.md](phase_alignment_microdynamics.md):
the misalignment $\mu$ as the sole relational datum of a pair, Proposition 3
(winding rates), and Theorem 4 (the vertex is a swap).

**Corrects** four statements of
[density_matrix_microdynamics_algorithm.md](../algorithm/density_matrix_microdynamics_algorithm.md),
which is the only prior treatment of the *position* representation in this
repository:

1. Its §2.3, *"The density matrix has no such uniform bound"* — true in the
   continuum, false on the crystal lattice, and the lattice is where this
   project's ontology lives. On $M$ sites with $\mathrm{Tr}\thinspace\rho = 1$,
   $\rho_{mm} \le 1$ and $\lvert\rho_{mn}\rvert \le 1/2$ state-independently,
   and the $\ell^1$ mass obeys $\sum_{mn}\lvert\rho_{mn}\rvert \le M$ (§8).
2. Its §3–§4 mean-field construction — the pair-Bohm velocity field, the
   bilocal quantum potential, the three-term phase accumulation and the
   $1/\rho$ regularisation of its open item 2 — is **not necessary**. The
   ladder route of §3 below is exactly linear in $\rho$, reads no field off
   the ensemble, and has no $1/\rho$ anywhere. Asymmetries 1 and 2 of its
   §7.1, and open items 2 and 4, dissolve.
3. Its §7.1 item 4 asserts the sign problem is *"a tougher version"* than the
   Wigner side's. Sharper: the position-side sign problem is not merely
   tougher, it is **structurally irremovable and lattice-divergent** (§7, §8).
4. Its §0 caveat that the obstruction is the *"imaginary coefficient"* of the
   kinetic operator locates the difficulty correctly but stops short. §7 gives
   the two obstructions in provable form: a hop *probability* of order
   $\delta t^{2}$, and a local-gauge orbit that sweeps
   $\arg\rho_{mn}$ around the whole circle while leaving the diagonal fixed.

**Contradicts, mildly, the framing that prompted this note.** The proposal
that self-conjugate particles represent observable probability densities is
correct and is *derived* below (§2). But it carries a consequence worth
stating flatly: self-conjugate particles are **not the movers**. A population
consisting only of self-conjugate particles has exactly zero population flux —
Proposition P4(i), verified to $0.0$ — so in the position representation
there is no advection at all. Every bit of motion of the observable density is
supplied by the pairs. The probability density does not stream; it is dragged.

---

## 1. Tutorial: the board, the rung, and two analogies

Put positions on a periodic lattice of $M$ sites, $x_m = m\thinspace a$ with
$a = L/M$. A density matrix is then an $M \times M$ board of complex numbers
$\rho_{mn} = \rho(x_m, x_n)$. Index each square by its

- **rung** $k = m - n$ — how far off the diagonal it sits, in lattice steps;
- **midpoint** $r = (m + n)\thinspace a/2$ — where the pair lives.

Rung 0 is the diagonal: the observable probability density. Rung
$\lvert k\rvert \ge 1$ is coherence between two places $ka$ apart.

**The draughts analogy, transposed.** The momentum-basis ladder of
[coherence_ladder.md](coherence_ladder.md) had the potential as the only
player, moving one diagonal step at a time on the $\rho(P, P')$ board. Here
the roles swap: the **potential never moves anything**, and the **kinetic
operator is the only player**, again moving one diagonal step at a time. The
potential's entire job is to wind the phase of each square in place.

**The clock-and-tape analogy.** Picture each world-particle as a stopwatch
pinned to a lattice site. Its hand turns at a rate set by the local potential
alone — $-V(X)/\hbar$ for a ket leg, $+V(X')/\hbar$ for a bra leg — so a
leg needs to know nothing but where it is. A pair is two stopwatches joined by
a tape; the only physical quantity the pair owns is the angle between the two
hands, $\mu$. Winding the hands does not move the pins. Moving a pin is a
separate, kinetic event, and it costs a quarter turn.

---

## 2. The ontology, and why the diagonal is real and non-negative

Each world-particle carries **a position and a phase, and nothing else**. Legs
come in two species: a positon is a ket leg, a negaton is a bra leg, and the
species distinction is the complex conjugation of the bra side, exactly as in
step 7. A pair of conjugate legs at $(X, X')$ with phases
$(\theta, \theta')$ is one Monte-Carlo sample of the element $\rho(X, X')$,
of modulus the sample weight and argument

```math
\mu \;=\; \theta - \theta' \;=\; \arg\rho(X, X').
```

Note what has become easy. In the momentum representation, $\mu$ had to be
built by transporting two clock phases to a common event, and the locality of
that construction was an open lemma (step 7, §6). Here each leg's clock rate
$\mp V(X)/\hbar$ depends on **its own position only**, so $\mu$ is the plain
difference of two carried numbers. The pump half of the locality burden is
discharged for free; what remains is the vertex, §6.

**Self-conjugate particles.** Hermiticity, $\rho_{nm} = \rho_{mn}^{*}$, applied
to a pair with $X = X'$ gives $\rho_{mm} = \rho_{mm}^{*}$: the weight of a
self-conjugate particle is real. Positivity of $\hat\rho$ gives
$\rho_{mm} \ge 0$. So a self-conjugate particle is a sample of a real
non-negative number — a unit of observable probability density, carrying
$\mu = 0$ by construction. This is the position-space counterpart of the
excess particle of the momentum ladder, and the reason it looks like a
classical particle is that its two legs coincide and its misalignment
vanishes identically.

---

## 3. Theorem P1: four one-leg hops and a diagonal pump

Take the nearest-neighbour kinetic operator on the ring,

```math
T \;=\; -J\thinspace(S + S^{\dagger}) + 2J,
\qquad J \;=\; \frac{\hbar^{2}}{2 m a^{2}},
\qquad (S)_{mn} = \delta_{m,\thinspace n+1},
```

and the diagonal potential $V_{mn} = V(x_m)\thinspace\delta_{mn}$. The
$+2J$ on-site term is common to both sides of the commutator and cancels
identically, so

```math
\partial_t\rho_{mn}
\;=\; \frac{iJ}{\hbar}\Bigl[\rho_{m+1,n} + \rho_{m-1,n}
      - \rho_{m,n+1} - \rho_{m,n-1}\Bigr]
      \;-\; \frac{i}{\hbar}\bigl[V(x_m) - V(x_n)\bigr]\thinspace\rho_{mn}.
```

**Theorem P1.** The von Neumann generator on the lattice is exactly the sum of
four **one-leg hop channels** and one **diagonal pump**:

| channel | leg moved | $\Delta k$ | $\Delta r$ | amplitude |
| --- | --- | --- | --- | --- |
| ket up | positon, $X \to X + a$ | $+1$ | $+a/2$ | $+iJ/\hbar$ |
| ket down | positon, $X \to X - a$ | $-1$ | $-a/2$ | $+iJ/\hbar$ |
| bra up | negaton, $X' \to X' + a$ | $-1$ | $+a/2$ | $-iJ/\hbar$ |
| bra down | negaton, $X' \to X' - a$ | $+1$ | $-a/2$ | $-iJ/\hbar$ |
| pump | neither | $0$ | $0$ | $-i\thinspace\Delta V/\hbar$ |

and nothing else. $\square$

Read off four structural facts, all of them checked in Part A of the demo at
machine precision (assembled flux versus the commutator,
$1.8 \times 10^{-16}$ relative; Hermiticity and trace, exactly zero):

1. **The mover is a leg of a bound pair, one site at a time.** The pair's
   midpoint advances by a *half* lattice step per event — the same
   half-quantum offset that the momentum ladder found for $q\thinspace dp$,
   and for the same reason.
2. **Direction is unbiased.** Up and down carry the *same* amplitude. Nothing
   in the kinetic operator prefers a direction; the hop is pure noise.
3. **The dynamics is in the quarter turn.** A ket hop multiplies the sample's
   weight by $+i$, a bra hop by $-i$. The negaton's opposite quarter turn is
   the complex conjugation of the bra side — the same identification that
   Theorem C2 made with its refractive factor.
4. **The potential moves nothing.** At $J = 0$ every modulus
   $\lvert\rho_{mn}\rvert$ and every population $\rho_{mm}$ is exactly
   constant (demo: $6.4 \times 10^{-19}$ and $0.0$); only $\mu$ winds.

---

## 4. Theorem P2: momentum is the misalignment of a rung-1 pair

Write the rung-1 element at the bond $m + 1/2$ as
$\rho_{m+1,m} = \lvert\rho_1\rvert\thinspace e^{i\mu}$.

**Theorem P2.** The lattice probability current is

```math
j_{m+1/2} \;=\; \frac{2Ja}{\hbar}\thinspace\mathrm{Im}\thinspace\rho_{m+1,m}
\;=\; \frac{\hbar}{m a}\thinspace\lvert\rho_1\rvert\thinspace\sin\mu ,
```

and the population equation of Theorem P1 is *exactly* the lattice continuity
equation $\partial_t\rho_{mm} = -(j_{m+1/2} - j_{m-1/2})/a$. In the smooth
limit $\lvert\rho_1\rvert \to \rho$ and $\mu \to \bar p\thinspace a/\hbar$, so

```math
\bar p \;=\; \frac{\hbar\thinspace\mu}{a},
\qquad j \;=\; \rho\thinspace\frac{\bar p}{m}. \qquad\square
```

Demo Part B: continuity to $1.1 \times 10^{-16}$, the population flux exactly
real, and $\hbar\mu/a$ recovering an imposed packet momentum
$p_0 = 1.7$ to $10^{-8}$ already at $M = 16$. On a chirped packet, where
$\bar p$ varies across the packet, the field $\hbar\mu(x)/a$ tracks
$m\thinspace j/\rho$ across the whole support.

This is worth stating in words, because it is the sharpest thing the position
representation buys. **Momentum is not a property a world-particle carries.**
It is the misalignment of a conjugate pair, per lattice step. A leg has a
place and a clock; a pair has, in addition, an angle; and that angle *is* the
momentum. Compare Proposition 3 of the phase-alignment note, which says that
in the momentum representation $\mu$ winds in space at rate $\Delta p/\hbar$ —
the same identity read the other way round.

Two corollaries fall out at once. An aligned pair ($\mu = 0$) carries no
current, which is why a self-conjugate particle is motionless. And momentum on
a lattice is **compact**: the current saturates at
$\lvert\sin\mu\rvert = 1$ and reverses beyond, which is the Brillouin
periodicity appearing as a property of the misalignment rather than as a
property of a coordinate.

---

## 5. Proposition P3: no noise, no force, transposed

The pump acts on the rung-1 element by
$\partial_t\rho_{m+1,m}\rvert_{\rm pump} = -(i/\hbar)\thinspace\Delta V\thinspace\rho_{m+1,m}$
with $\Delta V = V(x_{m+1}) - V(x_m) \simeq -a F$. Hence

```math
\partial_t j\rvert_{\rm pump}
\;=\; \frac{\hbar}{m a}\thinspace
      \partial_t\bigl(\lvert\rho_1\rvert\sin\mu\bigr)\Bigr\rvert_{\rm pump}
\;=\; \frac{\hbar}{m a}\thinspace
      \lvert\rho_1\rvert\cos\mu\thinspace\Bigl(\frac{aF}{\hbar}\Bigr)
\;\longrightarrow\; \frac{\rho\thinspace F}{m}.
```

**Proposition P3.** The pump alone reproduces the Euler force term. $\square$
Demo Part C confirms second-order convergence in $a$ over six refinements,
ratio $4.00$ at the finest pair.

This is the position-space reading of "no noise, no force." The hop channels
are direction-symmetric and carry no force whatsoever; the potential exerts no
push and moves no particle. Newton's law appears as the composition of the
two: the pump tips pairs out of alignment, misalignment is momentum, and the
unbiased hopping converts momentum into displacement. Kill either half and
transport stops — the pump alone freezes every modulus (§3, fact 4), the hops
alone are unbiased.

---

## 6. The vertex, and consistency with the phase-alignment layer

Theorem P1 gives amplitudes, not a mechanism. The mechanism must satisfy the
same demands the momentum ladder faced.

**Linearity forces a state-independent striker.** The flux out of element
$(m, n)$ must be linear in that element alone. A vertex between two state
pairs is bilinear. So, exactly as in §4 of the coherence ladder, the partner
in a hop must be drawn from a background whose leg density is uniform and
state-independent — a sea of position legs, one per site.

**Contact is adjacency; the vertex is a position swap.** In the momentum
representation a vertex is a co-location in $x$ and Theorem 4 makes it a
momentum swap. The transpose reads: a vertex is adjacency in $x$, and the
elementary move exchanges the *positions* of a state leg and a sea leg. Two
things follow without being imposed, in exact parallel to the momentum side:

- the per-site sea occupancy is preserved identically, because a swap is a
  permutation — the dual of $\sum p$ conservation;
- $\sum_{\rm legs} V(X)$ and $\sum_{\rm legs} V(X)^2$ over the two
  participants are preserved identically, because they are symmetric functions
  of the occupied sites — the dual of energy conservation following
  automatically from the swap rather than being imposed.

**Phase bookkeeping.** The moved leg's clock is continuous through the swap;
the quarter turn $\pm i$ of Theorem P1 is the vertex's own contribution, and
its sign is the species of the moved leg. This is the same "phase continuity
through the vertex" that Theorem C2 assumed.

**The honest caveat.** "Consistent with the phase-alignment model" here means
consistent *under the transpose*, not identical. The phase-alignment layer was
built with the potential as hopper and the kinetic term as phase pump; in the
position representation those two exchange jobs. One consequence is not
cosmetic and is developed in §8: the momentum-side vertex amplitude is
$V_q/2\hbar$, set by the potential and independent of the lattice, whereas the
position-side vertex amplitude is $J/\hbar = \hbar/2ma^{2}$, which diverges as
the lattice is refined.

**[open]** The number-conserving realisation. Theorem P1's generator has no
diagonal part — the $+2J$ terms cancelled — so a naive hop process that
*removes* weight from the source element overshoots. Either the removals
cancel against the sea's own churn (the dual of the striker-back-reaction
open item of the coherence ladder, §4) or the ensemble carries the uniform
compensating growth of §9. Only the second is verified here.

---

## 7. Proposition P4: two obstructions

**(i) Incoherent states are instantaneously static.** If $\rho$ is diagonal,
every $\rho_{m\pm1,m}$ vanishes and Theorem P2 gives
$\partial_t\rho_{mm} = 0$ exactly — demo Part D measures $0.0$, while the
rung-1 flux is nonzero, $3.6 \times 10^{-1}$: coherence is being created, and
only once it exists does anything move.

Consequence: **no lone-hopper process on self-conjugate particles can drive
the dynamics.** Any per-particle jump process with rates depending on the
particle's own state would move a diagonal state. This is the position-space
counterpart of the no-go lemma of step 7 §4, and it is the stronger of the
two: there the obstruction was spurious diffusion, here it is that the process
must not fire at all.

**(ii) The coherence sector is not a jump process.** Two independent proofs.

*Amplitude, not probability.* The one-leg update over $\delta t$ has
amplitude $\lambda\thinspace\delta t$ with $\lambda = J/\hbar$, so the hop
*probability* is $\lambda^{2}\delta t^{2}$ — second order. Any Poisson process
has hop probability first order in $\delta t$. Demo Part D exhibits the two
scalings side by side, and confirms the physical signature: releasing one site
under $T$ alone gives
$d\log\langle x^{2}\rangle/d\log t = 2.000$, ballistic, whereas a jump process
with bounded rates gives $1$.

*The gauge circle.* Let $\psi_m \to e^{i\chi_m}\psi_m$. The diagonal is
untouched and $\rho_{mn} \to e^{i(\chi_m - \chi_n)}\rho_{mn}$, so as the state
ranges over physically distinct states with the *same* observable density, the
argument of every off-diagonal element sweeps the entire circle. Demo Part D
measures the histogram of $\arg\rho_{mn}$ over random local gauges as flat to
within the sampling floor (relative spread $0.093$ on $24$ bins of $\approx 83$
counts, whose Poisson floor is $0.11$).

**Proposition P4.** No fixed background $b_{mn}$ makes $\rho_{mn} + b_{mn}$
lie on a single ray for all states; hence there is no position-space analogue
of the positon-only sea. $\square$

This is the precise reason the phase-space representation is privileged for
this project's ontology, and it is a one-line reason: **$W$ is real and the
shift $2/h$ is real, whereas $\rho_{mn}$ is complex with a gauge-free
argument.** The Wigner bound $\lvert W\rvert \le 2/h$ is what makes
$W' = W + 2/h \ge 0$ a probability density and the sea a *counting* device.
The density matrix has a bound (§8) but no ray.

What survives is the $`\mathbb{Z}_4`$ charge already anticipated in §2.2 of the
algorithm specification. Any $z \in \mathbb{C}$ is a non-negative combination
of $\lbrace +1, +i, -1, -i\rbrace$, with total count
$\lvert\mathrm{Re}\thinspace z\rvert + \lvert\mathrm{Im}\thinspace z\rvert$; the
background is then the balanced quartet, which sums to zero and so may be
present at any depth. Four species instead of two, and cancellation instead of
counting.

---

## 8. Proposition P5: resource arithmetic and the rate comparison

**The lattice restores the bound.** With $\mathrm{Tr}\thinspace\rho = 1$ on $M$
sites: $\rho_{mm} \le 1$; and since $\rho_{mm} + \rho_{nn} \le 1$,

```math
\lvert\rho_{mn}\rvert \le \sqrt{\rho_{mm}\thinspace\rho_{nn}} \le \tfrac12 ,
\qquad
\sum_{mn}\lvert\rho_{mn}\rvert
\le \Bigl(\sum_m \sqrt{\rho_{mm}}\Bigr)^{2} \le M .
```

The $\ell^1$ bound is the useful one. It says the total sampling mass needed
to represent $\rho$ is at most the number of cells — precisely the same
statement as on the Wigner side, where the sea mass is the number of momentum
rows. The two representations cost the same in $\ell^1$. Demo Part E, on the
canonical cosine-well trajectory at $M = 64$ over four periods: peak
$\ell^1$ mass $55.7$ against the bound $64$, peak $`\mathbb{Z}_4`$ sampling mass
$70.3$ against $\sqrt{2}M = 90.5$. The bound is tight, not generous — the same
verdict step 7 reached for the Wigner sea.

The *uniform* bound $1/2$, by contrast, is loose by a factor $\sim M$: a flat
position-space sea of fixed depth is unaffordable, even though a flat
phase-space sea of depth $2/h$ is not. That asymmetry is real but it is an
asymmetry about *flatness*, not about total mass.

**The rate comparison is where the representations genuinely part company.**

| representation | kinetic term | potential term | vertex amplitude |
| --- | --- | --- | --- |
| $W(x, p)$ crystal lattice | free advection, no rate | mediated jumps | $\Gamma_q \sim V_q/\hbar$ |
| $\rho(P, P')$ momentum ladder | phase winding, no rate | one-leg hops | $V_q/2\hbar$ |
| $\rho(X, X')$ position ladder | one-leg hops | phase winding, no rate | $J/\hbar = \hbar/2ma^{2}$ |

The first two amplitudes are set by the potential and are independent of the
lattice. The third is set by the lattice and diverges as $a^{-2}$. Demo Part E,
with $\hbar = m = 1$ and $L = 8$: the total channel rate
$\Lambda = 4J/\hbar$ is $8$ at $M = 16$ and $512$ at $M = 128$, against
$V_q/2\hbar = 0.75$ on the momentum side at any $M$.

**Proposition P5.** The position representation is the only one of the three
whose elementary vertex rate fails to have a continuum limit. $\square$

---

## 9. Theorem P6: statistical equivalence

Represent the state by $N$ pairs $`(X_i, X'_i, w_i)`$ with $w_i \in \mathbb{C}$,
and evolve:

- **pump** — continuously,
  $`w_i \mathrel{\ast}= \exp[-i(V(X_i) - V(X'_i))\delta t/\hbar]`$,
  read from the two legs' own positions and nothing else;
- **hop** — a Poisson clock of total rate $\Lambda = 4J/\hbar$; on firing,
  choose one of the four channels uniformly, move that leg one site, and
  multiply $w_i$ by $+i$ (ket) or $-i$ (bra);
- **estimator** — $`\rho_{\rm emp}(t) = e^{\Lambda t}\sum_i w_i \delta_{X_i}\delta_{X'_i}`$.

**Theorem P6.** $\mathbb{E}[\rho_{\rm emp}(t)] = \rho(t)$ for all $t$, on
every element of every state. $\square$ The uniform factor $e^{\Lambda t}$ is
forced: adding and subtracting $\Lambda\rho$ converts Theorem P1's
diagonal-free generator into a jump kernel plus a state-independent growth,
and the growth is a $c$-number that can be carried outside the ensemble.

Demo Part F, $M = 8$, $\Lambda = 2$: the error falls as $1/\sqrt N$ (ratio
$3.42$ per decade against $3.16$), and averaging $16$ independent runs reduces
the error by exactly $\sqrt{16}$ at every time tested — the deviation is
noise, not bias. The single-run error grows with $t$ in step with
$e^{\Lambda t}$: at $\Lambda T = 0.5, 1, 2, 3$ the errors are
$0.034, 0.054, 0.164, 0.411$ against $e^{\Lambda T} = 1.65, 2.72, 7.39, 20.1$.

**The cost, stated plainly.** Weight moduli are conserved, so the estimator is
$e^{\Lambda t}/N$ times a sum of $N$ unit phasors: signal per element
$\sim M^{-2}$, noise $\sim e^{\Lambda t}/\sqrt N$, hence
$N \sim M^{4}e^{2\Lambda t}$. For the canonical cosine-well parameters
($M = 128$, $T = 4\thinspace T_{\rm period} \approx 26$),
$\Lambda T \approx 1.3 \times 10^{4}$. The scheme is exact and unusable, and
it is unusable for the reason Proposition P4 identifies rather than for any
repairable inefficiency.

---

## 10. Corollary P7: the observable sector *is* a particle process

The obstruction of §7 applies to the coherence sector. The diagonal escapes
it, because Theorem P2 writes the population equation as a genuine continuity
equation with a *given* current. Split the current into its two signs and
divide by the local population:

```math
R_{m \to m+1} \;=\; \frac{\max(j_{m+1/2},\thinspace 0)}{a\thinspace\rho_{mm}},
\qquad
R_{m \to m-1} \;=\; \frac{\max(-j_{m-1/2},\thinspace 0)}{a\thinspace\rho_{mm}} .
```

**Corollary P7.** Both rates are non-negative by construction, and the master
equation they generate has $\rho_{mm}(t)$ as its exact solution. $\square$
Demo Part G: $200{,}000$ self-conjugate walkers hopping at these rates track
$\rho(x, x, t)$ over $1.5$ periods in the cosine well to a peak sup-norm error
of $0.0019$, against a shot-noise floor of $0.0022$ — at the floor. The
deterministic control with the same rates agrees with the exact density to
$1.2 \times 10^{-3}$, which is the time-step residual and not a defect of the
rates.

So the position representation splits cleanly in two, and the split is exactly
the observable/unobservable line:

- **self-conjugate particles** — real non-negative weights, positive hop
  rates, an honest stochastic process, guided by $\sin\mu$;
- **conjugate pairs** — complex weights, amplitude dynamics, no positive-rate
  realisation.

That is a Bohm/Nelson architecture arrived at from the pair ontology rather
than assumed: the pairs are the guiding field and the self-conjugate particles
are the beables. It is worth recording that this is a *retreat* relative to the
Wigner side, where the whole state, coherence included, is carried by one
positive process on $W' \ge 0$.

---

## 11. The two duality tables

**Term by term against the momentum ladder** (step 8):

| | $\rho(P, P')$ momentum ladder | $\rho(X, X')$ position ladder |
| --- | --- | --- |
| leg carries | momentum and phase | position and phase |
| diagonal operator (winds $\mu$) | kinetic, at $(P^2 - P'^2)/2m\hbar$ | potential, at $\Delta V/\hbar$ |
| hopping operator (moves a leg) | potential, by $\pm\hbar K$ | kinetic, by $\pm a$ |
| rung | $(P - P')/2q\thinspace dp$ | $(X - X')/a$ |
| rung 0 carried by | excess particles | self-conjugate particles |
| hop amplitude | $-i\thinspace e^{\pm i\phi_q}V_q/2\hbar$ | $\pm iJ/\hbar$ |
| amplitude depends on lattice | no | yes, as $a^{-2}$ |
| $\mu$ needs transport to define | yes (open lemma) | no |
| the relational quantity is | position winding of $\mu$ | $\mu$ itself, over one bond |
| what $\mu$ *is* | phase difference of two clocks | the momentum, $\bar p = \hbar\mu/a$ |

**Against the Wigner crystal lattice:**

| | $W(x, p)$ | $\rho(X, X')$ |
| --- | --- | --- |
| field | real | complex, Hermitian |
| charges | $\lbrace +1, -1\rbrace$ | $\lbrace +1, +i, -1, -i\rbrace$ |
| background | $2/h$, static, makes $W' \ge 0$ | balanced quartet, sums to zero |
| positivity | yes, by the Wigner bound | impossible, by Proposition P4 |
| $\ell^1$ cost | number of momentum rows | number of cells, $\le M$ |
| kinetic | advection, deterministic | four hops, amplitude $\pm iJ/\hbar$ |
| potential | Poisson jumps in $p$ | phase winding in place |
| coherence sector | one positive process | amplitudes only |
| observable sector | the same one process | a separate positive process (§10) |

---

## 12. Consequences and open items

1. **The position ladder is exactly linear in $\rho$.** It needs no velocity
   field, no bilocal quantum potential and no $1/\rho$ regularisation. That is
   a strict improvement over §3–§4 of the algorithm specification and the
   recommended basis for any revision of it.
2. **The verdict on the position representation is split.** As an *ontology*
   it is the most economical of the three — a leg is a place and a clock,
   momentum is an angle between two clocks, the potential is a pure winder,
   and the locality of $\mu$ is free rather than conjectural. As a *simulation
   scheme* it is the worst of the three, for two provable reasons
   (Propositions P4 and P5) rather than for want of technique.
3. **Open: the number-conserving vertex** (§6). The dual of the coherence
   ladder's striker-back-reaction item, and the same shape of question.
4. **Open: whether the divergence of $\Lambda$ is physical or an artefact of
   the nearest-neighbour $T$.** A spectrally exact $T$ has unbounded hop range
   instead of unbounded rate; either way the continuum limit of the position
   vertex is singular, but the two singularities are not obviously the same
   one and the comparison is worth making.
5. **Open: the compact-momentum reading.** $j \propto \sin\mu$ rather than
   $\mu$ makes momentum a compact variable at the microdynamic level. Whether
   this is the same compactness that the phase-space crystal lattice imposes
   through its momentum rows, or an independent one, is unresolved.
6. **Consequence for step 7's species question.** In this representation
   species = ket or bra unambiguously, and the pump's sign on a leg is
   $\mp V(X)/\hbar$ by species. This is a second, independent vote for the
   ket/bra reading of the species bookkeeping, still pending page 4 of the
   Cyganski deck.

![The position-space coherence ladder](https://raw.githubusercontent.com/billpage/wpmw/output/figures/position_pair_ladder.png)

Panels, left to right and top to bottom: the pair board with its rungs and the
two one-leg hops; momentum read off the misalignment of a rung-1 pair on a
chirped packet; the Euler force term emerging from the pump alone; the
amplitude-versus-probability obstruction; the gauge circle that forbids a
positon-only sea; the $\ell^1$ resource arithmetic against its bounds; the
unbiased but $e^{\Lambda t}$-noisy pair ensemble; and the self-conjugate
walkers tracking the exact probability density.
