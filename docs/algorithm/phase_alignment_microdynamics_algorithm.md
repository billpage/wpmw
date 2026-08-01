# Phase-Alignment Microdynamics Algorithm

**A specification for simulating the quantum Liouville equation from a live sea of permanently paired world-particles, with each pair one Monte-Carlo sample of a density-matrix element and the misalignment of its two transported clock phases the argument of that element.**

---

## 0. Status of this specification

This document specifies the simulation counterpart of
`docs/analysis/phase_alignment_microdynamics.md` under the ontology of
`docs/analysis/permanent_pairing_density_matrix.md`. Where
`docs/algorithm/phase_space_crystal_lattice_algorithm.md` takes the rate
field $\Gamma_q(x)$ as given and applies it as a stencil, and where
`docs/analysis/sea_dressed_microdynamics.md` realises that stencil as
channels against a *pinned* sea, this specification describes a **live
sea**: the pairs are dynamical, they carry phases, and $\Gamma_q(x)$ is
not an input but an output, assembled from vertex statistics.

The specification is therefore **more expensive and less accurate** than
the mesh algorithm for any problem the mesh algorithm can solve. Its
purpose is different: it is an executable statement of the ontology, so
that claims about what the model *is* can be tested rather than
asserted. Anyone wanting Wigner dynamics should use
`wpmwlib/phase_space_crystal_lattice.py`.

**Revision history.**

- *July 2026 (first revision).* Stored partnership was removed in favour
  of a relational misalignment over all ordered pairs, at the price of
  postulate (S), the sea carrier lock, following
  `docs/analysis/relational_pairing_and_carrier_lock.md`.
- *July 2026 (ladder addendum).* Following
  `docs/analysis/coherence_ladder.md`: interior ladder edges
  ($k \to k \pm 1$ between stored rungs) are struck by the **background
  sea**, not the excess — the sea-striker channel is structural, not
  optional; the erase bias is corrected to the licensed bilinear form
  (all transfers carry $\mu_1$ and freeze at $V = 0$); the §7 parity
  remark is scoped to rung 1; and Corollary 4.2's exit state is scoped
  to the exact channel. The four leg-local ladder channels reproduce the
  commutator elementwise at machine precision
  (`src/demo_coherence_ladder.py`).
- *July 2026 (second revision).* Permanent pairing is **reinstated** under
  the density-matrix reading of
  `docs/analysis/permanent_pairing_density_matrix.md`: a pair is one
  sample of $\rho(X, X')$ — positon the ket leg, negaton the bra leg,
  the pair's misalignment $\mu$ the argument of the sampled element —
  and an excess particle is an unpaired diagonal sample. Postulate (S)
  is **withdrawn** (nothing now needs a shared carrier: every vertex
  reads the misalignment of the one pair involved), and with it the open
  question of its dynamical preservation. The vertex of §5 acquires two
  channels, **write** and **erase**, whose combination is predicted —
  not yet demonstrated — to assemble the mediated QLE generator from
  purely local rules. The order parameter $Z_r$ survives as a
  diagnostic, not a postulate. Reconstruction (§7) changes: excess
  particles now sample $W_0$ rather than $W' = W + 2/h$, and coherence
  is read from split pairs rather than from a background subtraction.

This specification is also the momentum-basis, definite-state cousin of
`docs/algorithm/density_matrix_microdynamics_algorithm.md`: there the
pair ensemble carries complex weights under pair-Bohm flow in position
space; here the pairs carry unit weight and a transported phase, live on
the momentum half-grid, and change state only at vertices. The two meet
in the object they sample — $\rho(X, X')$ — and diverge in what is
fundamental and what is emergent.

Status of the pieces:

- §§1–5 (state, streaming, encounters, the two-channel vertex) are fully
  determined by the analysis notes and implementable as written.
- §6 (the pump) is determined at first order; the steady state under
  continuous pumping is **[open]** — but it is now a balance with both a
  source (the pump splits pairs: §5.2 of the pairing note) and a drain
  (erase vertices), so a fixed point can exist, which the consumable
  accounting of the first revision excluded.
- §7 (reconstruction) is fixed by the ontology.
- §8 (calibration) carries one constant and one sharp prediction — the
  **same-constant property** across the two channels — which is the
  load-bearing theorem of the pairing note (§7.3 there), stated here as
  a testable requirement rather than assumed.
- Items marked **[choice]** are implementation decisions not fixed by
  the physics; alternatives are noted.

Companion code: `src/demo_phase_alignment.py` verifies the kinematic and
vertex-level claims of §§2–5 in isolation;
`src/demo_pairing_resource_arithmetic.py` verifies the storage-capacity
arithmetic that makes §2.4 feasible; `src/demo_relational_pairing.py`
remains the record of why the first revision was tried and what survives
of it (its Theorem R4 factorisation backs the diagnostic of §2.5). A
full live-sea integrator is not yet in `wpmwlib`; §10 lists the
validation ladder it must climb, and rungs 5–6 are the decisive test
named in §7.4 of the pairing note.

---

## 1. Domain, grids and constants

Periodic ring $x \in [0, L)$, one spatial degree of freedom, $\hbar$ and
$m$ as in the rest of the project. The momentum half-grid is

```math
dp \;=\; \frac{\pi\hbar}{L},
\qquad p_n \;=\; n \cdot dp,
\qquad n \in \mathbb{Z}.
```

By Theorem 1 of the phase-resonance note, **fundamental particles occupy
even $n$**; odd $n$ carries two-leg interference content only. The
simulation must therefore initialise all worldlines on even sites and
must never assign an odd momentum to a single leg. This is a hard
invariant and a good assertion to compile in.

The potential is a finite Fourier sum

```math
V(x) \;=\; \sum_{q \ge 1} V_q \cos\!\big(K_q x + \phi_q\big),
\qquad K_q = \frac{2\pi q}{L},
```

and each mode $q$ supplies a momentum quantum $2q \cdot dp$.

Cells for encounter detection: $M$ bins of width $\Delta x = L/M$. The
cell width is a numerical parameter and must be checked for convergence;
it sets the encounter rate and hence enters the calibration of §8.

### 1.1 Row-index convention

Two different rows get called $`"n"`$ in the literature of this project,
and conflating them is a live source of error. This specification keeps
them apart by name:

- **$r$ — a particle's own row.** $p = r \cdot dp$ with $r$ **even**, by
  Theorem 1. Every worldline in the ensemble carries an $r$.
- **$n$ — a vertex midpoint row.** The label of a $q$-mode transition,
  not of any particle. Its parity is the parity of $q$, so for odd $q$ no
  particle ever sits on it.

The dictionary between them: a $q$-mode vertex acting on a particle at
row $r$ in direction $s = \pm 1$ has midpoint

```math
n \;=\; r + sq,
```

carries the particle to row $r + 2sq$, and involves the two even rows
$n - q = r$ and $n + q = r + 2sq$ — the particle's own in-row and
out-row, straddling the midpoint. Under this revision the midpoint is no
longer only a label: it is where the **split pair** created by the
vertex stores its coherence sample (§5), which is exactly the row the
mesh stencil reads.

`phase_space_crystal_lattice_algorithm.md` writes its stencil as
$\Gamma_q(x)\thinspace(W_{n+q} - W_{n-q})$ and transfers
$(m, n+q) \to (m, n-q)$: that $n$ is the **midpoint** in the sense above.
Any expression mixing the two conventions in one formula is wrong; when
in doubt, write the algorithm in $r$ and $s$, and translate to $n$ only
where the mesh stencil is being quoted.

### 1.2 What a cell is, and what a row is

The two labels are not the same kind of object.

A **momentum row is exact.** $dp = \pi\hbar/L$ is fixed by the ring
circumference, not by the mesh: a wavefunction on a ring of circumference
$L$ carries momentum in multiples of $2\pi\hbar/L$, and $W$, being
bilinear, carries it on the half-grid. There is no tolerance to tighten
and no limit to take. "Sharing a row" is an equality test on an integer,
not a proximity test. For the cosine-well parameters
($L = 8$, $\hbar = 1$) adjacent particle rows are separated by
$2\thinspace dp = 0.785398$, which is exactly $\hbar K_1$.

A **cell is a quadrature bin.** $\Delta x = L/M$ is a free numerical
parameter, swept in §10 rung 7. It is *not* a scattering cross-section:
the model contains no cross-section, the mediating mode has wavelength
$L/q$ spanning the whole ring at $q = 1$, and a genuine cross-section
would have to survive the $\Delta x \to 0$ limit whereas this must
vanish from it. The phase-space cell is
$\Delta x \cdot dp = 0.0245$ against $h = 6.2832$, about $1/256$ of a
Planck cell, so it is not uncertainty-limited either. Its only job is to
say which particles can meet at a vertex.

Nor does the cell width blur the phase. $\Phi$ is a field on spacetime,
so for two same-row legs evaluated at a common event the evaluation point
cancels identically:

```math
\Phi_i - \Phi_j \;=\; \Big(\theta_i - \frac{p\thinspace x_i}{\hbar}\Big)
- \Big(\theta_j - \frac{p\thinspace x_j}{\hbar}\Big),
```

which is Lemma 1's gauge-invariant combination and contains no $x$ of the
evaluation point. The same cancellation holds for the two legs of one
pair on different rows once the evaluation event is fixed, which is why
the pair misalignment of §2.3 is well defined wherever it is read.

The one place $\Delta x$ does enter is the pump, which samples $V$ at
each particle's own position, spreading phases within a cell by a
fraction $K_q \Delta x = 2\pi q / M$ of the signal — 4.9% at $M = 128$,
$q = 1$. This does not accumulate for co-moving legs: same-row legs share
a velocity, so each samples the same $V$ over a ring transit. For the two
legs of a **split** pair the sampled potentials differ, and the pair
phase acquires the beat evolution of §2.3 — under this revision that is
signal, not drift: it is the free evolution of the stored coherence.

![Cells and rows](https://raw.githubusercontent.com/billpage/wpmw/output/figures/lattice_cells_and_rows.png)

*(a) exact momentum rows against sweepable spatial cells; solid rows
carry particles at even $`r`$, the dashed midpoint is empty for odd
$`q`$, and the highlighted cell shows a vertex carrying an excess
particle from row $`r`$ to row $`r + 2q`$. (b) intra-cell pump spread
over a full run: bounded, not secular. Part H of
`src/demo_relational_pairing.py`.*

---

## 2. State representation

### 2.1 Worldlines

The state is an ensemble of world-particles, each

```math
\big(x_j,\; p_j,\; \theta_j,\; \varepsilon_j\big),
\qquad
\varepsilon_j = \pm 1,
\qquad
\theta_j \in [0, 2\pi),
```

organised into two populations:

- **excess particles** — unpaired positons sampled from the initial
  Wigner function $W_0(x, p)$ itself, *not* from $W' = W_0 + 2/h$. They
  are the diagonal samples of the pairing note: populations without
  coherence. This specification assumes $W_0 \ge 0$ (true for the
  Gaussian initial states of the demos); initial states with negative
  regions require initial **split pairs** (§2.3) sampling the initial
  coherences, an extension noted in §11 and not needed for the cosine
  well.
- **sea pairs** — $N_{\mathrm{pair}}$ permanent positon–negaton pairs,
  $B$ per cell, each leg carrying the worldline data above plus one
  integer: the index of its partner, **fixed for all time**. The positon
  is the ket leg, the negaton the bra leg. Pairs initialise *aligned*:
  both legs co-located on the same even row with equal phases (a dark
  pair, Lemma 1).

Nothing else is stored per particle. In particular no misalignment, no
splitting and no "excited" flag is stored: all of those are derived
(§2.3).

### 2.2 Pairing

**Partnership is stored and permanent.** This reverses the first
revision, for the reason established in §3 of the pairing note: a
partner index indeed adds nothing to the one-particle marginal
(Proposition R1 stands), but it carries exactly the information that
distinguishes the density matrix from its marginal — *which two points
belong to one sample* of $\rho(X, X')$. With per-pair gauge (only
intra-pair phase differences are physical) the misalignment is not the
coboundary of any global one-index field, so the premise under which
the index was redundant fails at the two-point level.

Postulate (S) of the first revision is **withdrawn**. It existed because
a relational misalignment over non-partners is gauge-variant and its
all-pairs average vanishes as $1/B$ without a shared carrier. Under
permanent pairing the vertex reads the misalignment of the one pair
involved, which is gauge-invariant by construction, and no sea-wide
phase convention is needed. The crystal regularity that (S) asserted
becomes, at most, an emergent property to *measure* (§2.5), not a
postulate to impose.

**Locality of the read.** At a vertex only one leg of a pair need be
present; the partner may be anywhere on the ring. Evaluating the pair
misalignment there uses the partner's stored worldline data
$(x_j, p_j, \theta_j)$ — numbers carried by the simulation state, the
computational form of "each leg carries the pair's clock, set at the
last co-interaction and evolved by its own transport". That this
involves no signalling is the locality lemma of the pairing note (§6
there), **required and not yet proven**; the specification adopts its
candidate mechanism and §10 rung 8 instruments it.

### 2.3 Derived quantities

For a pair $k$ with legs $a$ (ket) and $b$ (bra), and any evaluation
event $(x^{*}, t^{*})$:

```math
\Phi_j(x,t) \;=\; \theta_j \;+\; \frac{p_j\,(x - x_j) - E_j\,t}{\hbar},
\qquad
\mu_k(x^{*},t^{*}) \;=\; \Phi_a(x^{*},t^{*}) - \Phi_b(x^{*},t^{*}) \pmod{2\pi},
```

```math
\Delta p_k \;=\; p_a - p_b,
\qquad
\bar p_k \;=\; \tfrac{1}{2}(p_a + p_b).
```

A pair is **aligned** iff $\Delta p_k = 0$ and **split** otherwise. A
split pair stores one sample of the momentum-basis element
$\rho(p_a, p_b)$, with $\mu_k$ its argument: its midpoint $\bar p_k$
sits on the interference row $n$ of §1.1, and its $\mu_k$ at fixed $x$
runs at the beat frequency
$\big(E_a - E_b\big)/\hbar = \Delta p_k\thinspace\bar p_k / m\hbar$, the
free evolution of the stored element. Storing $\mu$ or $\Delta p$ would
be a bug: both are functions of leg data and (for $\mu$) of the
evaluation event, and caching invites stale reads.

### 2.4 Sea depth and storage capacity

$B$ pairs per cell. The feasibility constraint is no longer consumption
but **storage**, and it is guaranteed in advance: by the Wigner bound
$\lvert W\rvert \le 2/h$ — the same inequality that makes
$W' = W + 2/h$ non-negative — a sea of pair density $2/h$ can host the
coherence content of *any* state at polarisation at most 100% (§5.3 of
the pairing note, verified on the cosine-well trajectory by
`demo_pairing_resource_arithmetic.py`: peak load factor $0.96$, bound
saturated only by Gaussian peaks and the deepest fringes). Scale $B$ to
the ensemble as $B = \nu\thinspace(2/h)\thinspace\Delta x\thinspace dp$ with $\nu$ the
samples-per-unit-$`W`$ density shared with the excess sampling, so that
one split pair and one excess particle carry equal weight in §7. $B$
must still be swept (§10 rung 7), now as a shot-noise parameter only.

### 2.5 Diagnostics (demoted from postulate)

Per cell $x_m$ and row $r$, over sea legs in that cell and row:

```math
Z_r(x_m, t) \;=\; \sum_{j \thinspace\in\thinspace (x_m,\thinspace r)} e^{i\Phi_j(x_m,t)},
\qquad
\hat Z_r \;=\; \frac{Z_r}{N_r}.
```

Under the first revision $\lvert\hat Z_r\rvert = 1$ was postulate (S);
under this one it is a *measurement* of emergent sea regularity, useful
because the pump analysis of §6 predicts its behaviour and because
Theorem R4 makes it cheap. Two further censuses are now primary
diagnostics:

- the **split census** $S_n(x_m)$: count and mean phase of split pairs
  by midpoint row, the direct readout of stored coherence, compared in
  §10 rung 6a against the off-diagonal mass $C(t)$ of
  `demo_pairing_resource_arithmetic.py`;
- the **load factor**: split-pair count per cell against $B$, which the
  capacity bound of §2.4 caps at 1.

---

## 3. Free evolution

Between vertices every world-particle streams and winds:

```math
\dot x_j = \frac{p_j}{m},
\qquad
\dot p_j = 0,
\qquad
\dot\theta_j = \frac{p_j^2/2m - V(x_j)}{\hbar}.
```

Momentum is constant on a leg — the potential acts only through the pump
(§6) and through vertices (§5), never as a classical force. This is the
sharpest structural difference from the classical-positon Monte Carlo in
`demo_cosine_well_microdynamics.py`, where $\dot p = -V'(x)$; here the
entire force emerges from vertex statistics.

The two legs of a split pair have different momenta and stream apart;
the pair remains one object through its stored indices, and its $\mu$
evolves as §2.3 describes. No dynamics ever re-partners anyone.

Integrate with the same splitting used elsewhere in the project.
**[choice]** — exact free flight between vertex times is preferable to
fixed-step integration, since $p$ is piecewise constant and the phase
integral is then analytic on each leg.

---

## 4. Encounter detection

An **encounter** is the co-location, within a cell over a window
$\tau_e$, of an excess particle with one leg of a sea pair. For an
excess particle at row $r$ and mode $q$, two kinds of partner are
admissible, one per channel of §5:

- **write partner** — a leg of an *aligned* pair on the out-row
  $r + 2sq$, $s = \pm 1$;
- **erase partner** — the leg of a *split* pair standing on a row
  $r + 2sq$ whose partner (the mate) stands on the excess particle's own
  row $r$ — i.e. a stored sample of $\rho$ spanning exactly the
  transition the excess can make.

The per-cell bookkeeping is one pass over the sea: bucket aligned-pair
legs by row, and split pairs by their ordered row pair
$(r_{\mathrm{struck}}, r_{\mathrm{mate}})$. The loop is then

```
for each cell:
    bucket sea legs                       # one pass over the sea
    for each excess particle c at row r:
        for q in modes, s in (+1, -1):
            propose write against aligned bucket [r + 2sq]
            propose erase against split bucket   [r + 2sq, r]
```

Cost is $O(N_{\mathrm{exc}} + N_{\mathrm{sea}})$ per step. The window
$\tau_e$ is a numerical parameter entering only through the product with
the coupling; see §8.

Row selection is the stationarity condition of Theorem 4 in bucket form:
for the erase channel the transition beat and the pair beat share the
midpoint $n = r + sq$ and co-move exactly; for the write channel the
aligned pair's pump sideband toward row $r$ supplies the co-moving beat
(Lemma 5). Off-row proposals are simply absent from the buckets.
Implementations wanting the physically honest soft version may weight
off-stationary candidates by
$\lvert \mathrm{sinc}(\tfrac{1}{2}\dot\mu\thinspace\tau_e)\rvert$ instead;
**[choice]** — start with hard buckets, use the soft form to confirm the
neglected traffic is negligible at the chosen $\tau_e$.

**Sea strikers.** Encounters are not excess-only. A leg of a background
pair co-located with a leg of a *split* pair proposes the ladder
channels of §5.6, and this is structural: every interior ladder edge
must be linear in the stored-pair count alone, which only a
state-independent striker density provides
(`../analysis/coherence_ladder.md` §4). Excess strikers carry exactly
the population-boundary fluxes, whose rates are proportional to the
populations that move.

---

## 5. The vertex

The core of the specification. Every firing is the Theorem-4 momentum
swap between the excess particle and the struck leg; what distinguishes
the channels is the pair's state before and after. Throughout,
$r' = r + 2sq$ is the out-row and $n = r + sq$ the midpoint.

### 5.1 The write channel

**Before:** excess at row $r$; aligned pair with both legs at $r'$; the
pump has written sideband amplitude on the pair (its misalignment
$\mu_k$ carries Lemma 5's place-valued phase).

**Swap:** the excess and the struck leg exchange momenta,

```math
p_c: r \mapsto r',
\qquad
p_{\mathrm{struck}}: r' \mapsto r .
```

**After:** the pair is **split** with legs on $(r, r')$, midpoint $n$ —
one new stored sample of $\rho(p_r, p_{r'})$, which is exactly the
coherence the QLE writes for this transition. Population moved and
coherence recorded are one event; the ledger cannot come apart.

This channel realises Corollary 4.3 in reverse: pre-vertex the pair
offers the transfer only through its pump sideband, and post-vertex the
splitting is definite. Proposition R5 of the relational note ("a struck
pair exits at $\Delta p = 0$ and is spent") does not apply to this
channel — the struck pair exits split and *charged*, which is the §5.2
inversion recorded in the pairing note.

### 5.2 The erase channel

**Before:** excess at row $r'$ (note: on the *mate's* row — this is the
exact Theorem-4 configuration, $p_{\mathrm{in}} = p_b$); split pair with
struck leg at $r$ and mate at $r'$, storing $\rho(p_r, p_{r'})$ with
phase $\mu_k$.

**Swap:**

```math
p_c: r' \mapsto r,
\qquad
p_{\mathrm{struck}}: r \mapsto r' .
```

**After:** both legs at $r'$: the pair exits **aligned** (Corollary 4.2
holds verbatim here), the stored sample is retrieved, and the excess has
hopped down the same rung the write channel hops up. The pair returns to
the aligned pool and is immediately reusable — the storage accounting of
the pairing note in mechanism form.

### 5.3 Firing probability

For either channel draw $u \sim \mathrm{Uniform}(0,1)$ and fire if
$u < P$ with the **contact form**

```math
P \;=\; \sin^2\!\big(\lvert h \rvert\thinspace\tau_e\big),
\qquad
h \;=\; g_0 \;+\; g_1\thinspace A_k\thinspace e^{i s \mu_k},
```

where $\mu_k$ is the misalignment of the participating pair read at the
encounter event (§2.3) and $A_k$ its excitation amplitude — the pump
sideband weight $\mu_1 = V_q\tau_p/\hbar$ in **every** channel, with
the direction label $s$ orienting the ordered element as in §5 of the
alignment note. For the erase, the pump path and the stored path reach
the same final state and interfere, so the biased rate is the bilinear
cross term $\propto \mu_1 \cos(\mu - \Lambda)$ — licensed with no
reservoir, pump-proportional, frozen at $V = 0$. The earlier reading
"amplitude unity for a split pair" is withdrawn: an isolated pair's
linear-in-stored-contrast response is unlicensed
(`../analysis/coherence_ladder.md` §5). The **same
$g_0, g_1$ serve both channels** — this is the same-constant property,
the one load-bearing assumption left (§0), asserted here and tested at
§10 rung 6.

The cheaper affine form
$P = w_0 + \kappa\thinspace A_k\cos(s\mu_k)$ with the clamp
$\kappa \le \min(w_0, 1 - w_0)$ remains available **[choice]**; the
clamp must be asserted, not clipped, since a clipped probability breaks
the linear response this construction exists to reproduce.

The bare $g_0$ (or $w_0$) traffic is the no-noise-no-force floor: with
an aligned pair it swaps equal momenta and moves nothing; with a split
pair it fires the erase channel unbiased. Its net contribution to the
mean generator is predicted zero with variance set by $g_0$;
**[choice]** — the smallest $g_0$ consistent with the contact form's
linearisation is variance-optimal.

### 5.4 Phase update

On firing, **phase is continuous through the vertex for both legs**:
only momenta change, both worldlines carry their $\theta$ through, and
each leg's reference data are reset to the vertex event. For the write
channel the newly split pair's $\mu$ therefore launches from the
pump-written value at the vertex — which Lemma 5 makes place-valued and
which carries the $\pi/2$ offset that the commutator's $-i$ supplies in
the density-matrix bookkeeping. The **prediction**, stated and not
proven: continuity samples the written element with the correct phase.
This is the phase-continuity question of the alignment note (open item
2), unchanged in content, now testable per pair at §10 rung 6b.

On not firing, nothing changes; the encounter is discarded.

### 5.5 Invariants to assert

Every vertex must preserve, exactly and in floating point to rounding:

1. $\sum p$ over the two participating worldlines;
2. $\sum p^2/2m$ over the same (automatic for a swap — a violation
   indicates the exchange has been implemented as something other than
   a swap);
3. worldline count, species count, and **every partner index**;
4. the even-site invariant of §1 for every single leg;
5. the channel ledger: write increments and erase decrements the split
   census $S_n$ of §2.5 by exactly one at the transition midpoint.

The species of the struck leg (ket or bra) is not constrained by the
swap; the ordered element a split pair stores is read from its
(ket row, bra row) in §7, and both orientations of the same element are
admissible samples.

### 5.6 Ladder channels

A split pair's legs are struck too. The struck leg's **own** pump
sideband supplies a co-moving pattern for the exchange
$p_{\mathrm{in}} = P_{\mathrm{struck}} \pm 2q\thinspace dp$, stepping
the pair's rung $k \to k \pm 1$ with the striker drawn from the
background sea (§4) and the mate never queried. The four channels —
struck leg ket or bra, direction $\pm$ — with phase continuity, the
refractive $-i\thinspace e^{\pm i\phi_q}$ on ket strikes and its
conjugate on bra strikes, and the same constant as §5.3, reproduce the
von Neumann commutator elementwise at machine precision
(`../analysis/coherence_ladder.md` §3, `src/demo_coherence_ladder.py`).
The mate-sideband (compound) channels are expected neutral and are
omitted from the reference loop **[open]** (ladder note §6).

---

## 6. The pump

The potential enters by acting on every world-particle through the same
phase vertex, over a pump interval $\tau_p$:

```math
\theta_j \;\mapsto\; \theta_j - \frac{V(x_j)\,\tau_p}{\hbar}.
```

This is the kick form of the $V$ term already present in §3's winding;
apply one or the other over a step, not both. For an unpaired excess
particle it is pure gauge. For an aligned pair it displaces the
misalignment by Lemma 5's place-valued profile with amplitude
$\mu_1 = V_q\tau_p/\hbar$ per mode — populations untouched at first
order (Lemma 3), which in the pairing note's accounting is precisely
the statement that the pump writes **split-pair amplitude**: the pump
is the source of the write channel's bias and the sea's recharge after
erase. For a split pair the two legs sample $V$ at different positions
and the stored phase evolves accordingly (§2.3) — the interaction-picture
evolution of the stored element.

Regression tests the implementation should keep from the previous
revision: the sea's $(x, p)$ histogram is unchanged by the pump to
$O(V_p)$ (else the pump has been implemented as a force), and
immediately after a pump $\lvert\hat Z_r(x)\rvert = 1$ to rounding in
every occupied row of every cell (the sea starts crystalline; whether it
stays so is now a measurement, §2.5).

**[open]** The steady state. Under continuous pumping, write vertices
split pairs and erase vertices re-align them; the pump recharges
alignment bias at rate $\mu_1/\tau_p$. The fixed point of this balance
should supply the calibration constant of §8 from first principles. It
has not been derived — but unlike the first revision's version of this
gap, both the source and the drain now exist in the ledger, and the
storage demand they must balance is bounded by the capacity theorem of
§2.4, so the fixed point is not excluded by arithmetic.

---

## 7. Reconstruction

The Wigner function is assembled from both populations, with no
background subtraction:

```math
W_{\mathrm{est}}(x_m, p_n) \;=\;
\underbrace{\frac{1}{\nu\,\Delta x\thinspace dp}
\thinspace\#\{\text{excess in } (x_m, p_n)\}}_{\text{populations, even } n}
\;+\;
\underbrace{\frac{2}{\nu\,\Delta x\thinspace dp}
\sum_{k \thinspace\in\thinspace \mathrm{split}(x_m,\thinspace n)}
\cos\mu_k(x_m, t)}_{\text{coherences, midpoint } n},
```

the second sum over split pairs whose midpoint row is $n$ and whose legs
bracket the cell; each contributes the interference fringe of its stored
element read at the cell centre. Rung-1 pairs have midpoint parity of
$q$, so for odd $q$ the rung-1 term and the populations occupy disjoint
rows; **even rungs sit on even midpoint rows** and their contribution
adds to the population term there, as do uniformly offset (gray) pairs,
whose $\cos\mu \ne 1$ carries signed even-row content — the dual
sampling of the ladder note §7.

Aligned sea pairs are **not** binned: they are the medium. The uniform
$2/h$ of the mesh representation is the aligned sea itself, and it
enters $W_{\mathrm{est}}$ only through the events it mediates — which is
the point of the ontology.

Diagnostics worth accumulating besides $W_{\mathrm{est}}$: the split
census and load factor of §2.5, and the per-cell quadrature
$\mathrm{Re}(\hat Z_{r+2q}\thinspace\overline{\hat Z_{r}})$, whose
profile should reproduce $-V'(x)$ up to the calibration.

---

## 8. Calibration

The mean generator must reproduce

```math
\partial_t W_n \;=\; -\frac{p_n}{m}\,\partial_x W_n
\;+\; \Gamma_q(x)\,\big(W_{n+q} - W_{n-q}\big),
\qquad
\Gamma_q(x) = -\frac{V_q}{\hbar}\sin(K_q x + \phi_q),
```

with the corrected sign of
`docs/supplement/phase_space_crystal_lattice_supplement.md` §6.3. The
$n$ of this stencil is the **midpoint** label of §1.1; in the $r$ of
§§4–5 the same transfer reads $r = n - q \longrightarrow r + 2q = n + q$.

The advection term is exact by construction (§3). The collision term
carries **one overall constant** — the product of encounter frequency,
pump duty cycle and vertex coupling — shared by both channels (§5.3).
Fix it by requiring L1 exactness of the reconstructed evolution over one
step, then verify — do not re-fit — against $\Gamma_q(x)$ above.

The structural prediction behind this section, from §4 of the pairing
note: the write channel alone has x-only rates and therefore *cannot*
reproduce the stencil — it carries the no-go lemma's irreducible
momentum diffusion — while write plus erase, counted per pair with the
stored phases of §5.4, is predicted to assemble the mediated generator
$\Gamma_q(x)\thinspace W'(\text{midpoint})$ with the sea's aligned pool
supplying the $2/h$ part and the split census supplying the
$W(\text{midpoint})$ part. **The mediated counting is never coded**; its
emergence is the claim under test, and rungs 5–6 of §10 are its trial.

---

## 9. Reference loop

```
initialise:
    sample N_exc excess positons from W0            (even rows; W0 >= 0)
    populate B aligned pairs per cell per even row  (legs co-located,
        theta_a = theta_b, partner indices fixed forever)

for each step:
    pump:      theta_j -= V(x_j) * tau_p / hbar     for every world-particle
    stream:    advance x_j, theta_j on every leg    (p_j constant)
    bucket:    per cell: aligned legs by row;
               split pairs by (struck row, mate row)
    encounter: for each excess particle c in cell m at row r:
                   for q in modes, s in (+1, -1):
                       # write: aligned pair on the out-row r + 2sq
                       if aligned[m][r + 2sq] nonempty:
                           mu <- pair misalignment at (x_m, t)      # 2.3
                           P  <- sin^2(|g0 + g1*mu1*exp(i*s*mu)| tau_e)
                           if uniform() < P:  swap; pair now split  # 5.1
                       # erase: split pair with struck leg at r + 2sq
                       #        and mate on the particle's own row r
                       if split[m][(r + 2sq, r)] nonempty:
                           mu <- stored pair misalignment at (x_m, t)
                           P  <- sin^2(|g0 + g1*exp(i*s*mu)| tau_e)
                           if uniform() < P:  swap; pair now aligned # 5.2
                       # ladder: sea leg strikes a split pair's leg (5.6)
               for each split pair with a leg in cell m:
                   propose the four leg-local channels against the
                   background-sea bucket at the leg's sideband rows
    assert:    invariants of 5.5
    diagnose:  W_est per section 7;  split census S_n;  load factor;
               |Z_r| and quadrature per section 2.5
```

One firing per excess particle per step at most **[choice]**: with the
calibrated couplings the per-step firing probability is small and the
distinction is higher order, but an implementation must pick a rule and
sweep $\tau_e$ across it.

---

## 10. Validation ladder

In order; each rung is a regression test, and none should be skipped.

1. **Kinematics.** Parts A–B of `demo_phase_alignment.py`: invariance
   and winding rates of $\mu$, now read per pair.
2. **Vertex algebra.** Part C: the swap solution, and conservation of
   both $\sum p$ and $\sum p^2$ without imposition, in both channels.
3. **Pump.** Parts D–E: place-valued $\mu$; sea populations unchanged at
   $O(V_p)$; per-cell $\cos\mu$ profile proportional to $-V'(x)$.
4. **Single step against the stencil.** One pump-plus-vertex step on a
   random $W_0 \ge 0$, compared with one Euler step of the mesh form in
   `wpmwlib/phase_space_crystal_lattice.py`. Agreement at the shot-noise
   floor, improving as $N^{-1/2}$.
4a. **Ladder generator.** The four leg-local channels of §5.6 against
   the exact commutator, elementwise, on random states — already
   machine-verified at expectation level by `demo_coherence_ladder.py`
   (Part B: $1.7 \times 10^{-16}$ relative, exact $V = 0$ freeze); the
   integrator must reproduce it stochastically, with rung census and
   the $V = 0$ freeze as sub-diagnostics.
5. **The no-go control.** Run the cosine well with the **erase channel
   disabled**. The prediction from the no-go lemma of the pairing note
   is *failure*: spurious momentum diffusion of order
   $\lvert\Gamma\rvert(2q\thinspace dp)^2/2$, washed interference, and
   divergence from the mesh QLE. If write-only succeeds, the theory of
   §8 is wrong; this rung exists to be failable.
6. **The decisive test.** Write plus erase, cosine well, four periods,
   against the mesh QLE — the pair-ensemble Monte Carlo of §7.4 of the
   pairing note, with mediated counting never coded. Sub-rungs:
   **6a** split census $S_n(t)$ against the off-diagonal mass $C(t)$
   and load factor of `demo_pairing_resource_arithmetic.py`;
   **6b** stored phases against the mesh coherences (the
   phase-continuity prediction of §5.4);
   **6c** the same-constant property: a single $(g_0, g_1)$ calibrated
   on rung 4 must serve both channels here without refit.
7. **Convergence sweeps.** In $B$, in $\Delta x$, in $\tau_e$, in
   $\tau_p$. A result that has not been swept in all four is not a
   result.
8. **Locality instrumentation.** Log every remote-partner read (distance
   between the vertex and the mate at read time) and confirm the
   distribution is what streaming predicts; this does not prove the
   locality lemma but keeps its empirical profile in view.
9. **Free particle.** $V = 0$: no pump, no write bias, no split pairs,
   no net exchanges beyond cancelling bare traffic; evolution must
   reduce to ballistic streaming as in `demo_cat_state_microdynamics.py`.

---

## 11. Known gaps

- **Same-constant property** (§5.3, §8) — the single load-bearing
  assumption: split pairs mediate with the vertex constants of
  pump-excited pairs. Asserted, tested at rung 6c, not proven.
- **Locality lemma** (§2.2) — reading the mate's stored data at a remote
  vertex must be shown to require no signalling. Pairing note §6;
  rung 8 instruments it.
- **Phase continuity at the vertex** (§5.4) — continuity is predicted to
  sample written elements with the correct phase; rung 6b decides.
  Alignment note, open item 2.
- **Steady state** (§6, **[open]**) — the calibration is fitted, not
  derived; the pump/erase balance that should fix it now has both terms
  present but no derivation.
- **Even-$`q`$ modes** — narrowed by the ladder note §7: interior
  mediation is sea-struck at any midpoint parity, so what remains open
  for even $q$ is the population-boundary bookkeeping only. The cosine
  well ($q = 1$) does not touch it.
- **Striker back-reaction** — the sea-striker channel is structural
  (§4, §5.6), and each interior-edge event displaces a background leg,
  leaving the striker's own pair split by one quantum. This churn must
  be neutral in expectation for the one-body sector; detailed balance
  and family pairing make it plausible, and the integrator must
  instrument it (ladder note §4). The former "sea–sea optional" entry is
  superseded.
- **Compound channels** — mate-sideband strikes execute two ladder edges
  at once at first order; Theorem C2 achieves exactness without them, so
  their net must cancel. Conjectured mechanism: direction-symmetric
  cancellation as in the contact-vertex demo's R2 (ladder note §6).
- **Initial coherences** — §2.1 assumes $W_0 \ge 0$; general initial
  states need initial split pairs sampling $\rho_0$ off-diagonals.
  Straightforward in principle, unspecified in detail.
- **Multi-mode potentials** — §1 admits a Fourier sum and §4 loops over
  modes, with a split pair's own splitting selecting its mode; first-order
  additivity of the mode contributions to $\mu$ is assumed and not
  verified.
- **Species bookkeeping** — now carries a candidate resolution: species
  = ket/bra side of $\rho$ (§2.1), with the vertex indifferent to which
  leg is struck (§5.5). To be checked against page 4 of the Cyganski
  slide deck before being marked closed.

Three entries of previous revisions are closed rather than carried
forward. **Postulate (S) and its dynamical preservation** — withdrawn
with the postulate itself (§2.2); the split-pair beat that would have
degraded $\lvert\hat Z_r\rvert$ is now stored signal, not drift.
**Sea consumption** (Proposition R5) — inverted by the write channel:
a struck pair exits charged, not spent, and the pump re-biases the
aligned pool (§5.1, §6). **Cost** — the bucketed loop of §4 keeps the
$O(N_{\mathrm{exc}} + N_{\mathrm{sea}})$ complexity that Theorem R4
first achieved.
