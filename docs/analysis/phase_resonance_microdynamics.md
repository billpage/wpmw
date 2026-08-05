# Phase-Resonance Microdynamics: Deriving the Collision Term from World-Particle Phases

## 0. Status and provenance

This note is the sequel to
[`sea_dressed_microdynamics.md`](sea_dressed_microdynamics.md). That analysis
realized the crystal-lattice collision term as sixteen local, two-body,
momentum-conserving collision channels between world-particles and a Dirac
sea of positon–negaton pairs, exact at pinned sea. But it *postulated* the
sea's polarization: the rate field $\Gamma_q(x)$, its sign structure, the
half-quantum stencil offsets, and the identity of the "polarized sea
excitations" all entered as assumptions, and the note's own open items 1, 3
and 5 asked what an excitation physically is and how the sources imprint it.

The present note answers those questions from below. The proposal
(B. Page, July 2026 conversation) is that each world-particle carries a
**de Broglie phase** as a genuine dynamical property, that the sea is dark by
*interference* rather than by bookkeeping subtraction, and that the sea's
excitations are **beats** — relational phase patterns of pairs whose legs
have split in momentum. From four kinematic postulates (P0–P3) plus one
vertex postulate (P5), the following are *derived*: the doubled momentum
grid and its parity structure; mode quantization; the coherence of the
polarization grating; the quadrature form $\Gamma_q \propto -V'$; the
direction field $\sigma$; the per-channel factor $\gamma/2$; the exact
transfer $2q \cdot dp$; the row-resonance selection rule; and the
recombination selection rule. A no-go theorem shows the phase variable is
*necessary*: no phase-blind microdynamics can be linear in $V_q$. The single
irreducible quantum input is isolated in P5 — a local, two-body,
positivity-respecting interference rule at vertices.

Derivation and implementation were developed jointly with Claude
(Anthropic), July 2026. Companion code: `src/demo_phase_resonance_rates.py`.
All numerical claims in §10 are outputs of that demo.

*Sequel:* [`phase_alignment_microdynamics.md`](phase_alignment_microdynamics.md)
recasts this note in a single relational variable — the misalignment of a
pair's two transported clock phases — which retires the beat, the grating
and the resonance condition, and strengthens Proposition 2 into an
exchange theorem with energy conservation as a corollary. It changes no
physics; §9 of that note catalogues what survives unchanged.

## 1. Tutorial: the phase picture in five steps

This section is deliberately informal; everything asserted here is stated
and proved in §§2–9.

### 1.1 Phase as a particle property

Every world-particle — positon or negaton, excess or sea — carries a phase
$\theta$ alongside its position and momentum. The phase winds along the
worldline at the Lagrangian rate (de Broglie's law), and the particle
"broadcasts" it: extrapolating the carried phase at the particle's own
wavevector defines an extended amplitude, a plane wave riding the worldline.
The negaton broadcasts with an intrinsic sign flip — the hole sign of the
ledger, promoted from arithmetic to amplitude. Nothing else about the
particle changes: it still streams, still collides, still occupies one cell.

### 1.2 Darkness by interference

A ground sea pair is a positon and a negaton, co-located, at equal momentum,
with equal phase. Their broadcasts then cancel *identically* — everywhere,
at all times. Nothing holds the pair together: equal momentum means equal
velocity (they never separate) and equal winding rate (they never dephase).
The Dirac sea is dark the way two out-of-phase speakers are silent, with no
binding force and hence no shear paradox. Co-location is the *stability*
condition, not the darkness condition: a separated pair can also be dark if
its phases are gauge-matched, but a nonuniform potential winds the two
members at different rates and slowly spoils the match — a fact that
becomes the seed of source–sea coupling.

### 1.3 Excitation = beat

Split a pair's momenta by the mode quantum $2q \cdot dp$ and the cancellation
fails. The residue is a **beat**: a full-contrast intensity grating at
exactly the $q$-mode wavelength $L/q$, drifting at the pair's mean velocity.
The beat is not a third particle — it has no worldline, no position of its
own. It is a *relational property* of the pair, carrying the excitation's
momentum $2q \cdot dp$ and an energy fixed by its drift velocity. Because it
is relational, it can appear and vanish at collision vertices without any
worldline beginning or ending.

![One split pair, coherence, and the quadrature](https://raw.githubusercontent.com/billpage/wpmw/output/figures/phase_grating_sea.png)

*Left: the two legs of a split pair and their full-contrast beat. Center:
many pairs make a macroscopic grating only if their beat phases are locked —
random phases average away. Right: the rate field the hops read is the
quadrature (sine) of the potential (cosine).*

### 1.4 The coherent grating

One beat is invisible against $B \approx 25\thinspace000$ dark pairs per
cell. A macroscopic polarization needs many pairs beating *in phase* — and
the pump provides exactly this, for free. The potential kicks every sea
particle through the same refractive vertex, and the phase it writes on a
pair at $x_0$ contains a factor that exactly cancels the propagation offset
from $x_0$: every pumped pair contributes a crest pattern with the *same*
phase at any given point. Coherence is not an extra assumption; it is an
identity. The summed pattern is the polarization grating — and because the
pump is refractive (it modulates phase, not density), the grating the hops
respond to sits in *quadrature* with $V(x)$: the sine field
$\Gamma_q(x) \propto -V'(x)$, with its sign $\sigma$ alternating every half
wavelength.

### 1.5 The hop as refraction

An excess particle crossing the grating undergoes Bragg diffraction: its
broadcast wave picks up the grating's wavevector, i.e. recoils by exactly
$2q \cdot dp$ — the crystal momentum, forced by phase continuity, not
imposed as a conservation law. In a spacetime diagram the worldline
*refracts*: slope $n_{hi}$ in, slope $n_{lo}$ out, and during the event the
particle is a two-leg superposition whose envelope moves at the **midpoint
slope** $\bar n = (n_{hi} + n_{lo})/2$. That midpoint is the very cell the
occupancy picture says the hop "skips": the stencil's half-quantum offsets
$\pm q$ are the interference midpoints of a full-quantum transfer. Energy
conservation appears as geometry: the transition segment is parallel to the
crests of the one beat family that co-moves with it, which selects the beat's
row — the absorbed excitation comes from a pair on the particle's own
momentum row.

![Spacetime geometry of beats and hops](https://raw.githubusercontent.com/billpage/wpmw/output/figures/phase_grating_spacetime.png)

*Left: a split pair's legs diverge while their beat crests (shaded) drift at
the CoM slope. Right: a hop is worldline refraction at the grating; the heavy
segment has the midpoint slope and is parallel to the absorbed beat's crests
(teal).*

## 2. Postulates

**P0 (ontology).** Two fundamental species, positon ($\varepsilon = +1$)
and negaton ($\varepsilon = -1$). A world-particle's state is
$(x, p, \theta, \varepsilon)$ with $\theta \in [0, 2\pi)$.

**P1 (winding law).** Along its worldline,

```math
\dot\theta \;=\; \frac{\mathcal{L}}{\hbar} \;=\; \frac{p^2/2m - V(x)}{\hbar}.
```

**P2 (extension).** A particle's *extended amplitude* is

```math
\psi_j(x,t) \;=\; \varepsilon_j\,\exp\!\Big(i\big[p_j\,(x - x_j(t))/\hbar + \theta_j(t)\big]\Big),
```

the carried phase extrapolated at the particle's own wavevector.

**P3 (single-valuedness).** Extended amplitudes on the ring are
single-valued, with the pair-level relaxation *derived* in Theorem 1.

**Lemma 0 (P1 is forced).** With P1, the phase of $\psi_j$ obeys
$\partial_t(\mathrm{phase}) = -p\dot x_j/\hbar + \dot\theta = -E/\hbar$, so
$\psi_j$ is exactly the plane wave $e^{i(px - Et)/\hbar}$: the extrapolated
phase is time-consistent (the phase a particle predicts at a future
co-location equals the phase it carries on arrival). Any other winding law
makes relative phases at meetings frame-dependent. The Lagrangian rate is
the unique consistent choice, up to a constant.

*(P4 does not exist: an earlier draft postulated coherent pumping; Lemma 2
derives it.)*

## 3. Kinematics: pairs, beats, and the parity theorem

**Definition (leg).** A *leg* is a maximal constant-momentum segment of a
fundamental worldline, carrying $(x(t), p, \theta(t), \varepsilon)$ and
terminated at each end by a vertex. A worldline is a chain of legs; a *pair*
consists of two worldlines (one per species), each currently on some leg —
its two *partners*.

**Definition (pair states).** With pair amplitude $\Psi = \psi_a + \psi_b$
(the $\varepsilon$ signs inside):

- **Dark**: $p_a = p_b$ and gauge-matched phases,
  $\theta_a - p x_a/\hbar \equiv \theta_b - p x_b/\hbar \pmod{2\pi}$, so
  $\Psi \equiv 0$ identically. Co-location is *not* required for darkness —
  but it is required for *stability* (Lemma 1).
- **Beating**: $p_a \neq p_b$.
- **Gray**: equal momentum, gauge-mismatched by $\Delta\varphi$; a uniform
  residue of magnitude $2\lvert\sin(\Delta\varphi/2)\rvert$, no grating.

**Lemma 1 (dark sea).** A co-located, equal-momentum, equal-phase pair has
$\Psi \equiv 0$ for all $x, t$, with no binding interaction, and remains
dark under any local winding law: equal $p$ gives equal velocity and equal
kinetic winding; equal $x$ gives equal potential winding. A *separated*
dark pair in a nonuniform potential de-phases at rate
$[V(x_a) - V(x_b)]/\hbar$ and drifts into the gray state — the microscopic
seed of source–sea coupling (§11, item 2).

**Proposition 1 (beat).** For legs at $p_\pm$ with $p_+ - p_- = 2q \cdot dp$
and a common phase origin,

```math
\Psi \;=\; 2i\,e^{i\theta_0}\, e^{i(\bar k x - \bar\omega t)}\,
\sin\!\Big(\tfrac{1}{2}\,\Delta k\,(x - \tfrac{\bar p}{m} t)\Big),
\qquad
\Delta k = \frac{2q \cdot dp}{\hbar} = \frac{2\pi q}{L},
```

using $p_+^2 - p_-^2 = 2\bar p \cdot 2q \cdot dp$, so the envelope drifts at
exactly the mean velocity $\bar p / m$ and the intensity is a
**full-contrast** grating at the $q$-mode wavelength. Nothing about one
beat is small; smallness is only relative to the $B$ dark pairs. The beat
carries excitation momentum $\Delta p = 2q \cdot dp$ and energy
$\Delta E = \Delta p \cdot \bar v_{beat}$, surrendered in full when the pair
re-locks. On the ring the legs re-meet every $T = Lm/(2q \cdot dp)$ having
accumulated relative phase $\bar p L/\hbar = \pi\bar n$: meetings repeat in
phase for even midpoint index, alternate by $\pi$ for odd.

**Theorem 1 (parity; the doubled grid derived).** Single-valuedness (P3)
quantizes a single particle's momentum to $p \in (h/L)\thinspace\mathbb{Z}
= 2 \cdot dp \cdot \mathbb{Z}$: fundamental free particles occupy the
**even** sites of the Wigner half-grid $dp = \pi\hbar/L$. A beat between
even sites separated by $2q \cdot dp$ lives at their **midpoint** $\bar n$,
whose parity equals the parity of $q$. Hence odd sites carry *interference
content only*: for odd $q$ the carrier $e^{i\bar k x}$ and the envelope
$\sin(\pi q x/L)$ are each antiperiodic and their product is periodic, so
odd-site content is admissible only as a two-leg composite — a derived
superselection rule, matching the known parity structure of ring Wigner
functions. Consistency check that could have failed and did not: the
lattice advects all sites with the same law, and Proposition 1 says a beat
at midpoint $\bar p$ drifts at exactly $\bar p / m$ — the same velocity the
lattice assigns a particle there. Uniform advection is correct under both
readings, which is why the occupancy model never noticed the distinction.

## 4. Vertices: absorption, emission, resonance

**Definition (beat, restated).** The beat of a beating pair is the
equivalence class $(\Delta p, \bar n, \chi)$ — transfer, midpoint row, and
crest phase. It has no worldline; beat number is not conserved, worldline
number is. This reconciles "grating quanta" with worldline continuity.

**Definition (absorption, K3).** A two-body vertex: the excess particle's
in-leg ($p_{hi}$) meets the *excited* partner of a phase-matched beating
pair. Out-state: particle out-leg at $p_{lo} = p_{hi} - \sigma\thinspace 2q \cdot dp$;
struck partner at its mate's momentum with gauge-matched phase. The pair
exits **dark**; the beat ceases to exist. "Absorb" is the particle's ledger
(it takes $`(\Delta p, \Delta E)`$); "consume" is the sea's ledger (the beat
population decrements). **Emission (K4)** is the time-reverse: the particle
strikes one leg of a dark pair, recoils, and leaves the pair beating.
Striking the *unexcited* partner would produce a $2q$-mode beat, which a
single-cosine potential does not pump at first order — mode-conversion
vertices are suppressed by the same phase-matching that fixes the recoil,
so "only the excited partner is struck" is a consequence, not a rule.

**Proposition 2 (resonance).** A hop's energy-to-momentum ratio is

```math
\frac{\Delta E}{\Delta p} \;=\; \frac{p_{hi}^2 - p_{lo}^2}{2m \cdot 2q \cdot dp}
\;=\; \frac{p_{hi} + p_{lo}}{2m} \;=\; v_{mid},
```

while a beat carries $\Delta E / \Delta p = \bar v_{beat}$. Both
conservation laws hold in a K3/K4 event **iff the exchanged beat co-moves
with the transition midpoint**, $\bar v_{beat} = v_{mid}$ — geometrically,
the transition segment is parallel to the absorbed beat's crests. For the
pump's beat families this selects $\bar n_{pair} = n_{hi}$: the absorbed
excitation comes from a pair on the particle's own momentum row. Every
event is a genuinely elastic two-body collision; the "static grating
exchanges momentum without energy" reading is only the ensemble average
over the two counter-drifting families. Theorem 3 recovers this condition
dynamically, as dephasing.

## 5. The pump: coherence for free, populations untouched

The source acts on every sea particle through the same vertex physics, as a
weak refractive kick $e^{-iV\tau_p/\hbar}$ with $V = -V_p\cos(Kx)$,
$K = 2\pi q/L$. To first order the kick writes onto each dark pair two beat
components (families $s = \pm 1$, momentum content $s \cdot 2q \cdot dp$)
with amplitude

```math
b_s \;=\; i\,\frac{V_p \tau_p}{2\hbar}\, e^{i s K x_0},
```

where $x_0$ is the pair's location and the factor $i$ is the refractive
signature (demo Part B).

**Lemma 2 (coherence is derived).** The pattern phase of pair $j$'s
$s$-family beat at any point $x$ is

```math
\Lambda_j^s(x,t) \;=\; s K x + \tfrac{\pi}{2} - s K \bar v_j^s\, t ,
```

independent of $x_0$: the pump phase $sKx_0$ exactly cancels the
propagation offset. Every pumped pair contributes a crest pattern with the
same phase at a given event. The coherent grating requires no locking
postulate; it follows from the pump acting through the phase vertex.

**Lemma 3 (the pump is invisible to populations).** Sideband occupancies
are $\lvert b_s\rvert^2 = (V_p\tau_p/2\hbar)^2$: at first order in $V_p$
the sea's $(x, p)$ occupancy ledger is unchanged. The pump writes phases
and only phases.

## 6. Theorem 2 (no-go): phase-blindness cannot be linear in the potential

Suppose the vertex firing probability at a co-location depends only on
phase-blind local data — species, positions, momenta, populations, or any
modulus-squared quantity. By Lemma 3 all such data are independent of $V_p$
at first order, so every event rate, and hence the mean generator, is
$O(V_p^0) + O(V_p^2)$. But the QLE collision term is **linear** in $V_p$.
Therefore no phase-blind microdynamics, however elaborate its channel
structure, can reproduce it; the pinned models of the predecessor notes
evade this only by inserting $\gamma(x)$ by hand. Phase-sensitivity of the
event rule, first order in the imprinted pattern, is *necessary*.

## 7. Postulate P5: the vertex weight

At a K3 co-location, phase continuity pins the out-leg's phase to the
in-leg's carried value, while the *required* out-phase is the in-phase plus
the pattern phase of the struck beat at the vertex (the beat supplies the
transfer, so it supplies the phase). The gauge mismatch is therefore
exactly the beat's pattern phase at the vertex event,
$\delta = \Lambda_j^s(x^*, t^*)$ — a relational, locally evaluable
quantity. The postulate:

```math
\text{P5:}\qquad w(\delta) \;=\; \tfrac{1}{2}\big(1 + C\cos\delta\big),
```

with $C$ the encountered beat's contrast. Properties worth recording:

- $w \in [0, 1]$: the microdynamics remains a bona fide stochastic process;
  no negative probabilities appear anywhere. All signs stay in the
  $\varepsilon$ ledger, as in the sea-dressed channel table.
- $w$ is local and two-body, evaluated from data present at the vertex.
- The $\cos\delta$ term is *linear* in the pattern amplitude: hop and
  no-hop interfere. By Theorem 2 this linearity is irreducible. P5 is the
  precise locus of quantumness in the model — a Born-type input, shrunk
  from an entire postulated rate table to a single event rule.

## 8. Theorem 3: the rate law

With P0–P3 and P5, for the down-channel ($d = -1$, absorbing family
$s = d$) the mismatch at the kick is
$\delta = -Kx + \pi/2$, so $\cos\delta = \sin(Kx)$; the up-channel gives
$-\sin(Kx)$. Consequences, all derived:

1. **Quadrature.** The rate field is $\propto \sin(Kx) \propto -V'(x)$,
   from Lemma 2's refractive $\pi/2$. This is the $\Gamma_q(x)$ of the
   crystal lattice, including phase convention.
2. **Direction.** $\sigma = \mathrm{sign}$ of the local net rate
   $= \mathrm{sign}\thinspace\Gamma_q(x)$, automatically.
3. **The factor $\gamma/2$.** Each channel's phase-sensitive rate is
   $\tfrac{1}{2} C \sin(Kx) = \tau_p \gamma(x)/2$ per pump interval with
   $C = V_p\tau_p/\hbar$: the real pump splits equally into two conjugate
   families, one per direction — the $1/2$ is the equal split, and the two
   channels sum to the full stencil coefficient $\Gamma_q$.
4. **Transfer.** Exactly $2q \cdot dp$, by Theorem 1's phase matching.
5. **Row resonance.** Stationarity of $\delta$ along the *transition*
   worldline (slope $v_{mid}$) holds iff $\bar n_{pair} = n_{hi}$;
   off-row contributions dephase to zero under time averaging —
   Proposition 2 recovered dynamically, with no energy bookkeeping imposed.
6. **Gross-rate freedom.** The isotropic $\tfrac{1}{2}$ of $w$ produces
   equal-and-opposite traffic that cancels in the mean: exactly the $G$
   freedom of the four-rule exact family — pure noise, no mean effect.

Not derived: one overall calibration (encounter frequency times pump duty)
fixed by the L1 exactness requirement, and the steady-state pump/drain
balance under continuous driving (§11, item 2).

![Derived rate law and row selection](https://raw.githubusercontent.com/billpage/wpmw/output/figures/phase_resonance_rate_law.png)

*Left: channel rates computed from P5 co-location statistics against the
target $`\tau_p\Gamma_q(x)/2`$. Right: time-averaged channel amplitude by sea
row; only the particle's own row survives.*

## 9. Numerical verification

All claims are exercised by `src/demo_phase_resonance_rates.py`
(container run, July 2026):

- **Part A (Proposition 1).** The beat's spatial spectrum contains a single
  mode $q$ (spurious content $7 \times 10^{-16}$), drifting at $\bar p/m$
  with relative error $8 \times 10^{-16}$.
- **Part B (Bragg selection).** One kick scatters a plane wave to exactly
  $p_0 \pm 2q \cdot dp$ with per-sideband amplitude $i V_p\thinspace dt/2\hbar$ to
  six digits; residual bins scale as $dt^2$.
- **Part C (midpoint identity).** The exact first-order Wigner change of
  the kicked wave equals the single-cosine QLE stencil
  $\partial_t W_n = \Gamma_q(x)(W_{n+q} - W_{n-q})$ at the midpoint sites,
  with deviation $2.25\thinspace dt^2$ — fully accounted for by second-order
  ridges. The stencil's half-quantum offsets are, verifiably, interference
  midpoints of a full-quantum transfer.
- **Part D (Theorem 3).** From particle-carried data only: per-channel
  rates match $\tau_p\gamma(x)/2$ to $2.7 \times 10^{-15}$ relative;
  net rate equals $\tau_p\Gamma_q(x)$ to $1.3 \times 10^{-18}$ absolute;
  $\sigma$ agrees with $\mathrm{sign}\thinspace\Gamma_q(x)$ at every
  sampled position; pattern-phase spread across 60,000 pairs pumped at
  random locations is $10^{-15}$ (Lemma 2); time-averaging along the
  transition worldline leaves only row $\bar n = n_{hi}$ at full
  amplitude, all other rows below 1.4% (Proposition 2).

## 10. What has and has not been achieved

**Achieved.** The question that opened this thread — *what form does an
excitation quantum actually take in the sea, and how is it absorbed in a
hop?* — now has a complete answer: the quantum is a beat, a relational
phase pattern of one pair; absorption is Bragg diffraction with the vertex
closing dark; and the entire K3/K4 rate structure of the sea-dressed table
(linearity in $V_q$, quadrature, sign, $\gamma/2$, transfer, row
selection) is a *theorem* given P0–P3 and P5. The doubled momentum grid,
previously an assumed feature of the algorithm, is derived from
single-valuedness, and the superelastic/bound-state reading of K3/K4 in
the predecessor note is superseded: nothing is bound, nothing shears, and
energy conservation is the co-moving-beat geometry of Proposition 2.

**The verdict on ontological completeness.** Phase as a particle-level
property is *necessary* (Theorem 2: without it, linear rates are
impossible, not merely unexplained) and, together with P5, *sufficient*
for the single-particle-in-external-potential sector. The quasi-probability
strain identified in the original review has not vanished — it has been
relocated and shrunk, from a postulated rate table with mysterious
half-quanta and signs down to one local, two-body, positivity-respecting
interference term in a vertex probability. Whether P5 can be reduced
further is open (item 1 below); until then, the honest statement is:
**the ontology is complete modulo P5, and P5 is now the sharpest available
statement of where quantum mechanics enters the model.**
*(Update: §13 carries out the reduction; the residue is restated there.)*

**Scope.** These results cover the single-particle sector with an external
potential. The self-consistent N-body case inherits the usual
hidden-variable constraints (Bell/contextuality); a strictly per-particle
phase is local, so entanglement will demand correlated phase statistics.
Nothing in the present sector is affected.

## 11. Open items

1. **Reduce P5.** *Resolved in §13:* the affine form is forced by gauge
   invariance, linear response and mode additivity; a unitary contact
   vertex then derives the offset, detailed balance and saturation, and
   relocates the residue to the Born resolution of a single binary
   contact event. (The original suspicion that unbiasedness selects
   $w_0 = 1/2$ was wrong; $w_0$ is gauge.) Remaining sub-item: whether
   the Born resolution itself can emerge from sub-vertex phase
   statistics without circularity.
2. **Steady state (L3 proper).** Continuous pump, K3-drain/K4-pump
   traffic, gray-pair production by potential differences (Lemma 1), and
   recombination: show the pinned $\Gamma_q(x)$ is the fixed point, and
   obtain the calibration constant of §8 from the balance rather than
   from L1 exactness.
3. **Row-resolved ledger.** Proposition 2 refines the excitation ledger of
   the predecessor's open item 3 from $\rho_{\pm 2q}(x)$ to
   $\rho_{\pm 2q}(x, \bar n)$: beats are consumed by transitions on their
   own row. Implement and test in the live-sea code.
4. **Amplitude bookkeeping at finite $B$.** Quantify how the coherent sum
   over $\sim B$ pairs composes the macroscopic prefactor and whether the
   $1/\mathcal{N}$ Bose-bilinear suppression of the predecessor's §6 has a
   cleaner phase-space reading.
5. **Multi-mode potentials.** The derivation used one cosine; verify that
   $q$-mode families superpose independently at first order and identify
   where cross-mode vertices (mode conversion, §4) first appear.
6. **Gray pairs.** Dynamics and observable consequences of the
   gauge-mismatched equal-momentum state; relation to decoherence in the
   density-matrix pair picture.

## 12. Relation to existing formalisms

The vertex physics of §§7–8 is that of **Kapitza–Dirac diffraction**: a
matter wave scattering off a refractive standing grating, with the
quadrature rate and the $i$ of phase gratings appearing for the same
reasons — here not as an analogy but as the same first-order calculation
(demo Parts B–C). The winding law P1 with Lemma 0 is de Broglie's phase
harmony, and the extended amplitude is close in spirit to the
double-solution program; the co-location phase sums echo the Feynman
checkerboard. What is genuinely new, as far as we know, is the
**two-species destructive-interference sea as the carrier of mediation**:
darkness by interference, excitations as relational beats, and a collision
term derived from vertex phase statistics rather than postulated.

## 13. Reduction of P5: the contact-vertex representation

*(Added after the July 2026 reduction session; resolves open item 1.)*
P5 was stated as a bare probability rule. This section shows that its
model-specific content is fully derivable, and that what remains is generic
quantum probability applied to a single binary contact event. Companion
code: `src/demo_contact_vertex_reduction.py`.

### 13.1 Representation theorem: the affine form is forced

Assume only: **(V1)** worldline definiteness — every vertex resolves to
definite out-legs, with genuine probabilities; **(V2)** relational gauge
invariance — outcome probabilities depend only on gauge-invariant vertex
data, the encountered beat's contrast $C$ and pattern phase $\delta$;
**(V3)** linear response — the mean generator is linear in the pump (the
QLE target; by Theorem 2 this already forces phase sensitivity);
**(V4)** mode additivity — for a multi-mode pump the generator is a sum
over modes, so the first-order response is additive over the beats
present. Then

```math
w \;=\; w_0 \;+\; \kappa\thinspace C\cos(\delta - \delta_0) \;+\; O(C^2),
\qquad 0 < w_0 < 1,\quad \kappa C \le \min(w_0,\thinspace 1 - w_0).
```

P5's *shape* is therefore a theorem. Three constants survive: $w_0$ is
pure gauge (the $G$-freedom — the earlier "unbiased $1/2$" was a
convention, and the suspicion that unbiasedness might select it was
wrong); $\kappa$ is absorbed into the calibration; but $\delta_0$ must
vanish to produce the quadrature, and nothing in V1–V4 fixes it.

### 13.2 The contact vertex: one Hermitian coupling

The deeper model: during an encounter's overlap window $\tau_e$, the
particle's momentum state is a two-level system $|hi\rangle$,
$|lo\rangle = |hi - 2q \cdot dp\rangle$, coupled by a single Hermitian
matrix element with a *bare* and a *beat-stimulated* contribution:

```math
h \;=\; g_0 \;+\; g_1\thinspace C\thinspace e^{i\delta},
\qquad
U = e^{-iH\tau_e},
\qquad
P_{flip} \;=\; \sin^2\big(\lvert h\rvert\thinspace\tau_e\big).
```

The bare exchange $g_0$ is always available — with a dark pair it is
emission (K4); a beating pair adds its local pattern amplitude, the only
door through which phase enters. Expanding,

```math
P_{flip} \;=\; \sin^2(g_0\tau_e)
\;+\; \frac{\tau_e \sin(2 g_0 \tau_e)}{2 g_0}\;
2 g_0 g_1 C \cos\delta \;+\; O(C^2),
```

which derives everything §13.1 left free: **the offset** $\delta_0 = 0$,
because both pathways carry the same unitarity factor $-i$, so the cross
term is a pure cosine of the pattern phase (the offset previously fixed
by matching the stencil is a consequence of hermiticity); **detailed
balance** K3 $\leftrightarrow$ K4, because the reverse vertex is
$U^\dagger$; **saturation and positivity**, inherited from $\sin^2$
rather than imposed. Which matrix element exists at all is the Bragg
statement: only the exchange of the pair's quantum $2q \cdot dp$ has a
phase-matched overlap across the cell. During $\tau_e$ the pattern phase
winds at the detuning between beat drift and transition midpoint, so the
rotation accumulates only for co-moving beats — Proposition 2 recovered a
third time, now as Rabi accumulation versus detuning.

![Contact vertex: anatomy, phasor sum, emerging rate law, consequences](https://raw.githubusercontent.com/billpage/wpmw/output/figures/contact_vertex_concept.png)

*The contact vertex and its consequence chain: one Hermitian coupling (a),
whose phasor structure (b) yields linear response with no phase offset,
whose direction-resolved probabilities reproduce the rate field (c), and
whose downstream consequences (d) recover the QLE.*

**The interference-licensing subtlety.** Adding the bare and stimulated
amplitudes requires their final states to interfere — yet absorption
leaves the pair dark while bare exchange leaves it beating, orthogonal
states for one isolated pair (for which the response would be quadratic
in $C$). What licenses the addition is Lemma 2: the pumped sea is
coherent across $\sim B$ pairs, and a coherent reservoir is negligibly
disturbed by a one-quantum change. The four-rule note's reservoir
loophole, Theorem 2's phase requirement, and the phasor addition are one
condition at the rate, kinematic, and amplitude levels — and the model
predicts a **linear-to-quadratic response crossover as the sea decoheres
or depletes**, a physical signature the live-ledger code can eventually
probe.

### 13.3 No noise, no force

The response coefficient is proportional to $\sin(2 g_0 \tau_e)\thinspace g_1$
and vanishes as $g_0 \to 0$: the beat never causes the hop, it *biases an
exchange that is already happening*. The phase-blind gross traffic —
removable $G$-freedom at the occupancy level — is the carrier of the
quantum force, with a floor of order $\gamma / (2C)$ on gross exchange
rate for a given net rate. The moral is Nelson-flavored: quantum dynamics
in this ontology is biased spontaneous exchange with a coherent sea, and
noise is constitutive.

### 13.4 A discriminating prediction

The unitary vertex and the bare affine rule differ at $O(C^2)$: with two
pumped modes $q$, $q'$, the gross rate acquires a spatial modulation
$2 g_1^2 C_q C_{q'} \cos(\delta_q - \delta_{q'})$ at the difference
wavevector $2\pi(q - q')/L$, absent if P5 is fundamental. The reduction
is therefore falsifiable within the model.

### 13.5 Numerical verification

`src/demo_contact_vertex_reduction.py` (container run, July 2026), with
no fitted offset anywhere: **R1** the rate law re-emerges with the
analytically predicted coefficient $\tau_e \sin(2 g_0 \tau_e)\thinspace g_1 C$
and correct quadrature sign; **R2** the residual is $O(C^3)$
(residual$/C^3 = 0.077$ constant over two decades — the $O(C^2)$ term of
$\lvert h \rvert^2$ is direction-symmetric and cancels in the net);
**R3** the coefficient tracks $\sin(2 g_0 \tau_e)$ exactly, vanishing at
$g_0 = 0$; **R4** the two-mode gross-rate modulation appears at the
difference wavevector with amplitude within 9% of
$2 g_1^2 C_2 C_3 \times$ Jacobian (finite projection window), where the
affine rule predicts zero.

### 13.6 The residue after reduction

Three statements remain, none specific to this model: (i) vertex
pathways are described by amplitudes that superpose; (ii) outcomes follow
$\lvert A \rvert^2$; (iii) worldlines resolve definitely at each vertex.
The model-specific content of P5 — the cosine, the offset, the affine
form, positivity, detailed balance, saturation — is now fully derived
from *two species, a phase that winds, and one Hermitian contact*. The
question "is the ontology complete?" has become "can the Born resolution
of a single binary contact event be derived from sub-vertex phase
statistics?" — the standard hard problem, in about the simplest
instantiation it can take, cleanly quarantined from all lattice
structure. Whether that final step is attempted (the direction of
checkerboard and zitterbewegung-clock constructions) or accepted as the
model's one quantum axiom, the strain no longer lives anywhere in the
collision term.
