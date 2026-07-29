# Interaction Diagrams for the Phase-Alignment Microdynamics

**Every elementary process of the phase-alignment layer, drawn twice — once
as a lattice picture in the $(x, p)$ plane and once as a set of
trajectories in the $(x, t)$ plane — in the format of
[`phase_space_crystal_lattice_supplement.md`](phase_space_crystal_lattice_supplement.md)
§4.**

> **Provenance.** The drawing format — a momentum-jump event shown as a
> phase-space lattice picture beside a space-time trajectory picture — is
> David Cyganski's, from the *Wigner Collisions Diagram* Sozi deck
> reproduced in §4 of the crystal-lattice supplement. The processes drawn
> here are those of
> [`docs/analysis/phase_alignment_microdynamics.md`](../analysis/phase_alignment_microdynamics.md)
> and its sequel
> [`relational_pairing_and_carrier_lock.md`](../analysis/relational_pairing_and_carrier_lock.md),
> specified for simulation in
> [`docs/algorithm/phase_alignment_microdynamics_algorithm.md`](../algorithm/phase_alignment_microdynamics_algorithm.md).
> Figures and the numerical checks annotated on them are produced by
> `src/gen_phase_alignment_interaction_diagrams.py`. §7 records a sign
> discrepancy that the exercise of drawing the pump exposed.

---

## 0. Scope, conventions, and one caveat about the title

### 0.1 A caveat on "wave–particle interaction"

The phrase is inherited from the earlier layers and does not survive into
this one. In the spawning and pair-diversion readings of crystal-lattice
supplement §4.1–4.2 a quantum of the potential field is absorbed or
emitted, and the momentum quantum $2q\thinspace dp = h/L$ is Abraham's
photon momentum for the mode. In the phase-alignment layer there is no
such quantum: the analysis note is explicit that *nothing is absorbed and
nothing is emitted* (§6), and every momentum-changing event is a
**two-body exchange between world-particles**, with the mediating sea leg
taking exactly what the excess particle gives up.

What remains of "the wave" is two things, and the diagrams label them as
such:

1. the potential's Fourier mode $q$, which acts on world-particles
   **only through their clock phases** (Process 1 below), and
2. the transported clock phase $\Phi$ itself, carried along worldlines
   with no field degrees of freedom of its own (analysis note §12,
   *direct interparticle action*).

The identification $2q\thinspace dp = h/L$ with the photon momentum
survives as an exact bookkeeping identity — crystal-lattice supplement §5
is untouched — but at this layer the momentum goes to another worldline,
not to a field. The figures are therefore titled *processes* rather than
*wave-particle interactions*.

### 0.2 Reading conventions

Both panel types use $x$ as the horizontal axis, so a process can be read
straight across from one to the other.

| element | meaning |
| --- | --- |
| solid horizontal lines in $(x,p)$ | momentum rows at even $r$, spacing $dp = \pi\hbar/L$ |
| dashed horizontal lines in $(x,p)$ | vertex-midpoint rows $n$, empty for odd $q$ |
| orange, filled | the excess particle (a positon above the sea) |
| teal, open | a sea leg |
| purple dial | the misalignment $\mu$ between the two mediating rows |
| purple lens in $(x,t)$ | the encounter window $\tau_e$ |

The dial's long orange hand is the transported phase $\Phi$ on the
out-row, the short teal hand is $\Phi$ on the in-row, and the shaded
sector between them is $\mu$. Parameters throughout: $\hbar = m = 1$,
$L = 8$, $q = 1$, hence $dp = 0.392699$ and
$2q\thinspace dp = \hbar K_q = 0.785398$; $V_q = 0.35$,
$\phi_q = 0.7$.

### 0.3 The inventory

There are exactly **four** elementary processes, plus the free leg.
§6 argues that the list is complete.

| # | process | changes $x$ | changes $\theta$ | changes $p$ |
| --- | --- | --- | --- | --- |
| 0 | free leg | yes | yes | no |
| 1 | pump | no | yes | no |
| 2 | exchange, $s = +1$ | no | yes | yes, $+2q\thinspace dp$ |
| 3 | exchange, $s = -1$ | no | yes | yes, $-2q\thinspace dp$ |
| 4 | suppressed traffic | no | no | no, in the mean |

---

## 1. Process 0 — the free leg

![Free leg: space-momentum and space-time](https://raw.githubusercontent.com/billpage/wpmw/output/figures/pa_int_0_free_leg.png)

*(a) the leg drifts along its own row and never leaves it. (b) the same
leg as a straight worldline of slope $`m/p`$; the dials show the clock
advancing.*

Between vertices,

```math
\dot x_j = \frac{p_j}{m},
\qquad
\dot p_j = 0,
\qquad
\dot\theta_j = \frac{p_j^2/2m - V(x_j)}{\hbar}.
```

**In phase-alignment terms.** This is the sharpest structural difference
from the classical-positon Monte Carlo of
`demo_cosine_well_microdynamics.py`, where $\dot p = -V'(x)$. Here the
potential appears in the third equation and nowhere else: it sets the
*winding rate* of the carried clock and exerts no force. A worldline in
panel (b) is straight through any potential landscape whatever. All of
the force in the model is recovered from vertex statistics, and the only
thing the potential does to make that possible is bias the clocks.

---

## 2. Process 1 — the pump

![The pump: space-momentum and space-time](https://raw.githubusercontent.com/billpage/wpmw/output/figures/pa_int_1_pump.png)

*(a) the pump leaves every sea leg where it is and writes a phase
grating across the ring; $`\cos\mu^{+}(x)`$ is in quadrature with
$`V(x)`$. (b) three legs on one row at three places: the clock kick
differs from place to place because $`V`$ does, and no worldline bends.*

The potential acts on every sea particle over a pump interval $\tau_p$,

```math
\theta_j \;\mapsto\; \theta_j \;-\; \frac{V(x_j)\thinspace\tau_p}{\hbar},
\qquad
\lvert\mu_1\rvert = \frac{V_q\thinspace\tau_p}{\hbar},
```

and by Lemma 3 the $(x, p)$ histogram of the sea is unchanged at
$O(V_q)$. Panel (a) is therefore a picture in which *nothing moves*. What
changes is the comparison between rows: by Lemma 5 the misalignment
becomes a function of the evaluation event alone,

```math
\mu^{s}(x,t) \;=\; s\thinspace(K_q x + \phi_q) \;+\; \frac{\pi}{2}
\;-\; s K_q\thinspace\bar v^{\thinspace s}\thinspace t,
\qquad s = \pm 1,
```

so that every pair of legs at the same place carries the same $\mu$, and
an arriving particle can read one number off the sea without knowing
anything about which legs it is reading.

**Why a phase kick produces a row-to-row misalignment at all.** As
written, the pump multiplies the amplitude on *every* row by the same
factor $\exp(-iV(x)\tau_p/\hbar)$, which on its own would leave any
row-to-row phase difference untouched. The content of Lemma 5 is that
this factor is not a constant: to first order in $\mu_1$,

```math
\exp\negthinspace\left(-\thinspace\frac{iV(x)\tau_p}{\hbar}\right)
\;\approx\;
1 \;-\; \frac{i\mu_1}{2}\thinspace e^{\thinspace i(K_q x + \phi_q)}
\;-\; \frac{i\mu_1}{2}\thinspace e^{-i(K_q x + \phi_q)},
```

so a phase-locked row acquires two coherent **sidebands** displaced by
$\pm\hbar K_q = \pm 2q\thinspace dp$ — exactly the two mediating rows of
a $q$-mode vertex — each of relative amplitude $\mu_1/2$ and each in
quadrature with $V$. The script confirms the amplitude ratio numerically
at $8.750335 \times 10^{-3}$ against the predicted
$\mu_1/2 = 8.75 \times 10^{-3}$, the residual being the $O(\mu_1^2)$ term
that Lemma 3 discards. This is Kapitza–Dirac diffraction of the sea, and
it is the whole of the "wave" content of the model.

**Populations versus coherence.** The distinction matters and is easy to
lose: the pump moves no population between rows at first order (Lemma 3),
but it does move *coherence* between them. Panel (a) draws the first
statement and panel (b) the second.

---

## 3. Process 2 — the exchange, $s = +1$

![Up-exchange: space-momentum and space-time](https://raw.githubusercontent.com/billpage/wpmw/output/figures/pa_int_2_exchange_up.png)

*(a) the excess particle hops from row $`r`$ to row $`r+2q`$ while one
sea leg hops the other way; the midpoint row $`n = r+q`$ is empty. (b)
the same event as trajectories: the two worldlines cross and exchange
slopes.*

An excess particle at row $r$ shares a cell with sea legs on rows $r$ and
$r + 2q$. By Theorem 4 the only exchange that survives averaging is the
**swap**: the particle arrives on the in-row leg's momentum $p_b$ and
leaves on the out-row leg's momentum $p_a$, while the struck leg takes
$p_b$. It fires with probability

```math
P \;=\; w_0 \;+\; \kappa\thinspace\mathrm{Re}
\big(\hat Z_{r+2q}\thinspace\overline{\hat Z_{r}}\big)
\;\xrightarrow{\ \text{(S)}\ }\;
w_0 \;+\; \kappa\thinspace\cos\mu^{+},
```

the second form holding under the carrier-lock postulate (Theorem R3 of
the sequel).

**The reading panel (b) makes available.** The union of the two
worldlines is *exactly two straight lines crossing*. Before and after the
vertex the set of occupied trajectories is what it would have been with
no interaction at all; the vertex permutes only the **labels** — which
line is "excess" and which is "sea". The script checks this literally:
the swapped path set and the free-crossing path set agree to
$0.000 \times 10^{0}$.

Three of the note's results are read off that one observation:

- **Corollary 4.1** (energy conservation is automatic). A permutation of
  labels over a fixed set of momenta preserves $\sum p$ and $\sum p^2$
  identically. Verified at $0$ exactly, with neither imposed.
- **Corollary 4.3** (transfer size). The particle can only move to a row
  the sea already occupies, so the hop is $\lvert\Delta p\rvert = 2q\thinspace dp$.
- **Theorem 4** (uniqueness). Since the event is a crossing of two
  straight lines, the only velocity that keeps $\mu$ frozen through the
  encounter window is the mean of the two slopes; imposing
  $v_{\mathrm{exch}} = \bar v_{\mathrm{pair}}$ then forces
  $p_{\mathrm{in}} = p_b$, $p_{\mathrm{out}} = p_a$. Verified at zero
  residual for both directions.

Because the trajectories are untouched, everything physical about the
vertex lives in the label count — which is precisely the excess
population $W = \rho_{\mathrm{emp}} - 2/h$ that the reconstruction of
§7 of the specification bins. The sea mediates without ever being a
force.

---

## 4. Process 3 — the exchange, $s = -1$

![Down-exchange: space-momentum and space-time](https://raw.githubusercontent.com/billpage/wpmw/output/figures/pa_int_3_exchange_down.png)

*(a) the same vertex run downward: rows $`r`$ and $`r-2q`$, midpoint
$`n = r-q`$ empty. (b) with $`p_{\mathrm{out}} = 0`$ the outgoing
excess worldline is vertical, and the union of the two paths is again
untouched.*

Structurally identical, and drawn separately only because the two
directions carry **opposite bias**. Under Lemma 5 the two conjugate
families satisfy $\cos\mu^{s} = -\thinspace s\sin(K_q x + \phi_q)$, so the
signed rate assembled from both is

```math
\frac{1}{\tau_p}\sum_{s = \pm 1}
\tfrac{1}{2}\thinspace\mu_1\thinspace s\thinspace\cos\mu^{s}(x)
\;=\; -\thinspace\frac{V_q}{\hbar}\thinspace\sin\negthinspace\big(K_q x + \phi_q\big)
\;=\; \Gamma_q(x),
```

with the corrected sign of crystal-lattice supplement §6.3. The
bare rate $w_0$ contributes equally to both directions and cancels in the
mean — this is the $G$-freedom of
[`four_rule_microdynamics_equivalence.md`](../analysis/four_rule_microdynamics_equivalence.md),
and it is why gross traffic is not drawn as a fifth process: it is
Processes 2 and 3 running at equal rate.

Panel (b) is drawn with $p_{\mathrm{out}} = 0$ deliberately: the outgoing
excess worldline is then vertical and the incoming sea leg supplies the
lower half of that same vertical line, which makes the label-permutation
reading of §3 hard to miss.

---

## 5. Process 4 — the suppressed channels

![Suppressed channels: space-momentum and space-time](https://raw.githubusercontent.com/billpage/wpmw/output/figures/pa_int_4_suppressed.png)

*(a) candidate transitions out of row $`r`$: only the stationary one
survives averaging, and an incoherent sea is dark by itself. (b) the
event axis must ride the pair's mean velocity for $`\mu`$ to hold still;
inset, the dephasing envelope.*

Nothing is forbidden in this model. Off-stationary candidates are
available at every encounter and simply average away over the encounter
window, with the time-averaged bias

```math
\big\langle\cos\mu\big\rangle_{\tau_e}
\;=\; \mathrm{sinc}\negthinspace\left(\frac{\dot\mu\thinspace\tau_e}{2}\right)
\cos\negthinspace\left(\mu_0 + \frac{\dot\mu\thinspace\tau_e}{2}\right),
\qquad
\dot\mu = \frac{\Delta p}{\hbar}
\big(v_{\mathrm{exch}} - \bar v_{\mathrm{pair}}\big),
```

confirmed against direct numerical averaging to $5 \times 10^{-6}$ at
every detuning tested. Panel (b) draws why: the dials on a path at
$\bar v_{\mathrm{pair}}$ hold still, and the dials on any other path
turn.

A second and independent suppression is drawn at the bottom left of panel
(a). A sea with random phases has $\lvert\hat Z_r\rvert = O(N_r^{-1/2})$,
hence a vertex bias of $O(N^{-1})$ by Corollary R4.3 of the sequel: an
incoherent sea is dark without a separate darkness postulate. The two
suppressions are conceptually different — one is dephasing *during* an
encounter, the other is the absence of a macroscopic polarisation to
read — and only the second is what the carrier-lock postulate (S) rules
out.

---

## 6. Why the list is complete

A world-particle carries $(x, p, \theta, \varepsilon)$ and nothing else
(§2.1 of the specification). Worldlines are neither created nor destroyed
(§5.5, invariant 3), and $\varepsilon$ never changes. So an elementary
process can only be a rule for changing $x$, $p$ or $\theta$:

1. $x$ and $\theta$ change continuously, by Process 0, and the streaming
   law admits no free choices.
2. $\theta$ changes discontinuously only through the pump, Process 1,
   which is the only appearance of $V$ outside the winding rate.
3. $p$ changes only at a vertex, and by Theorem 4 a vertex is a swap.
   A swap is fixed once its direction is chosen, and by Corollary 4.3 its
   size is $\lvert\Delta p\rvert = 2q\thinspace dp$. Hence exactly two
   momentum-changing processes per Fourier mode: Processes 2 and 3.

For a potential with several Fourier modes there are $2\times$(number of
modes) exchange processes, one pair per mode, and one pump. First-order
additivity of the mode contributions to $\mu$ is assumed and not verified
— open item 3 of the analysis note and the third known gap of the
specification.

Process 4 is not a fifth entry in the same sense: it is the statement
that the candidates *not* in the list are suppressed dynamically rather
than excluded by fiat, and it is drawn because Theorem 4's uniqueness
claim is only meaningful against the alternatives it defeats.

---

## 7. A sign discrepancy the pump diagram exposes

Drawing Process 1 required committing to a sign for the clock kick, and
doing so numerically rather than by inspection turned up a factor of
$-1$ that the analysis notes and the specification do not currently
agree on. Recorded here as a question, not as a correction.

**The check.** Build a sea row as a carrier plane wave on the ring, apply
the pump as a literal multiplicative phase, decompose into rows, and read
$\mu^{s}$ off the sidebands as the difference of transported phases.
Part C of `src/gen_phase_alignment_interaction_diagrams.py`:

| kick applied | $\max\lvert\mu^{s} - \text{Lemma 5}\rvert$ | assembled rate vs. $\Gamma_q$ |
| --- | --- | --- |
| $\theta \mapsto \theta - V\tau_p/\hbar$ (specification §6) | $3.142$ (i.e. $\pi$) | $-\thinspace\Gamma_q$, to $3.8\times10^{-15}$ |
| $\theta \mapsto \theta + V\tau_p/\hbar$ | $5.3\times10^{-15}$ | $+\thinspace\Gamma_q$, to $4.4\times10^{-15}$ |

The sideband amplitude comes out at $\mu_1/2$ in both cases, so the
identification of the sidebands with the mediating rows is not in
question — only an overall sign is.

**Where the factor can live.** Exactly one flip is needed, and there are
two candidate homes.

- *The kick sign in §6 of the specification.* This is the one
  Hamiltonian evolution gives: the field phase at fixed $x$ advances as
  $-Vt/\hbar$, and §3's Lagrangian winding rate
  $\dot\theta = [p^2/2m - V]/\hbar$ agrees. On those grounds §6 looks
  right, which points away from this candidate.
- *The direction label $`s`$ in §6 of the sequel.* A swap moves the excess
  particle and the struck sea leg in **opposite** directions, so
  "direction $s$" is ambiguous between the two, and resolving it the
  other way flips the sign of the whole sum. This is the live candidate,
  and it sits exactly where open item 1 of the analysis note already sits
  — the row-index convention of page 4 of the Cyganski slide deck, still
  unresolved.

Note that two apparently plausible fixes do **not** work. Reversing the
orientation of $\mu$ (comparing in-row to out-row rather than the
reverse) changes nothing, because the vertex rule depends on $\mu$ only
through $\cos\mu$. Flipping Lemma 5's $\pi/2$ offset to $-\pi/2$ alone
also fails: it gives $\cos\mu^{s} = +s\sin(K_q x + \phi_q)$ and the sum
still returns $-\Gamma_q$.

**What the figures assume.** The diagrams are drawn with the physically
anchored convention — the net drift of the excess population is towards
lower $V$, which is the diagnostic of crystal-lattice supplement §6.3 —
and Lemma 5 is quoted as printed. If the resolution turns out to be the
kick sign, the caption of Process 1 panel (b) is the only thing that
changes.

---

## 8. Companion code and sources

- `src/gen_phase_alignment_interaction_diagrams.py` — Parts A–C
  (the swap, the dephasing envelope, the pump sign) and all five figures.
- `src/demo_phase_alignment.py` — the kinematic and vertex-level
  verification of the analysis note.
- `src/demo_relational_pairing.py` — the order-parameter reformulation
  used for the vertex rule quoted in §3.
- David Cyganski, *Wigner Collisions Diagram* (Sozi presentation) — the
  drawing format adopted here.
