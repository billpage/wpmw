# Glossary

Terms used in the notes of [`analysis/`](analysis), [`supplement/`](supplement)
and [`algorithm/`](algorithm), with what each means *here*. A word is in this
list for one of three reasons:

1. it is used in **three or more notes**, and a reader has to carry it from one
   to the next;
2. it has **more than one sense**, or an ordinary English sense that misleads
   (*pair*, *row*, *floor*, *reading*);
3. it is **borrowed from outside this project's own field** and used in a
   specific sense (*eikonal*, *detailed balance*), however rarely.

Each entry says what the word denotes, where it is **defined** (the place to
read for the full statement and its hypotheses), and what to **not confuse** it
with. The entry is a pointer, not a substitute: where an entry and a note
disagree, the note is right, and the entry needs fixing.

Two companions answer related questions. **Where is a word used?** Run
`PYTHONPATH=src python3 -m wpmwlib.where_used eikonal`, which lists every file,
heading and sentence, definitional ones first, and cannot drift.
**Where is a theorem or open item?** [`analysis/INDEX.md`](analysis/INDEX.md).
**What is the project about, and in what order do the notes go?**
[`ORIENTATION.md`](../ORIENTATION.md).

This file replaces the two vocabularies that used to sit in
`ORIENTATION.md` §9 and `supplement/what_the_reach_is.md` §7; their entries
are merged here.

To suggest a term, tell whoever is editing; no process is needed. A tutorial
should still define a term in place at first use, in a sentence, and link here
for the rest.

**Status.** October 2026. Mathematical symbols are as in the notes an entry
links to. Words that are common in the notes only because they are common in
physics (*wavepacket*, *potential*, *momentum* in its ordinary sense) are not
listed; *momentum* is, for its one special sense.

[A](#a) · [B](#b) · [C](#c) · [D](#d) · [E](#e) · [F](#f) · [G](#g) · [H](#h) · [I](#i) · [K](#k) · [L](#l) · [M](#m) · [O](#o) · [P](#p) · [Q](#q) · [R](#r) · [S](#s) · [T](#t) · [V](#v) · [W](#w)

---

## A

### absorbing boundary

The statement that beyond the [coherence horizon](#coherence-horizon) a world
is exactly free: all vertex activity is confined to
$`\mathrm{dist}(x, \mathrm{supp}\,V) \le L_c/2`$ (Proposition O4). A world
that leaves the active region keeps streaming; nothing eats it.

- **Defined:** [`open_position_space.md` §3](analysis/open_position_space.md#3-the-coherence-horizon),
  Proposition O4.
- **Not to confuse with:** a sponge layer or absorbing boundary condition of a
  numerical scheme.

### absorptive, emissive

The two ways to realise one residual event. Both put $`+1`$ at row
$`p + \xi_q`$ and $`-1`$ at row $`p - \xi_q`$ in the observable, so they agree in
$`W`$ and differ only in the ledger. **Emissive:** ionise a sea pair at the
parent's row; its positon becomes a body at $`p + \xi_q`$ and its negaton one
at $`p - \xi_q`$ (sea $`-1`$ pair, $`N`$ up by 2). **Absorptive:** remove a
negaton body at $`p + \xi_q`$ and a positon body at $`p - \xi_q`$, and a new sea
pair appears at the parent's row (sea $`+1`$ pair, $`N`$ down by 2). The
absorptive fraction is $`f`$; $`f = 1/2`$ is the sinkless case, not a law
(Theorem N3).

- **Defined:** [`ORIENTATION.md` §4](../ORIENTATION.md#4-what-the-residual-events-actually-are);
  [`stochastic_ledger.md` §1](analysis/stochastic_ledger.md#1-the-process);
  [`emission_and_absorption.md` §4](supplement/emission_and_absorption.md).
- **Not to confuse with:** absorption or emission of light. Nothing here
  radiates; the words name the sign of the sea's change. See
  [recombination](#recombination).

### admissibility

Postulate (A) of step 17: the ensembles that occur are exactly those whose
expected $`W`$ is the Wigner transform of some density operator
$`\rho \succeq 0`$. It is what makes a signed $`W`$ physical, and it is a
postulate separate from the compensated split, because compensation does not
by itself remove the possibility of an inadmissible state.

- **Defined:** [`compensated_ontology.md` §1](analysis/compensated_ontology.md#1-four-postulates);
  [`ORIENTATION.md` §7](../ORIENTATION.md#7-the-proposed-ontology).
- **Not to confuse with:** "allowed" in general. *Admissible negativity*
  (Proposition V1 of `sea_depletion.md`) means negativity that $`W`$ may have
  while $`\rho \succeq 0`$; it is the negativity the dynamics creates for free
  in a harmonic trap.

### advection

Transport of a quantity along a flow. Here: the free-streaming step
$`\partial_t W + (p/m)\,\partial_x W = 0`$ on the lattice, a shift of
position cells by an amount fixed by the row's momentum. It is exact when
$`\Delta t`$ makes the maximum shift a whole number of cells; the spectral
form is available only on a ring.

- **Defined:** [`phase_space_crystal_lattice_algorithm.md` §3a](algorithm/phase_space_crystal_lattice_algorithm.md#3a-free-streaming-advection);
  the ring and the open line in [`compensated_liouville_algorithm.md`](algorithm/compensated_liouville_algorithm.md).
- **Not to confuse with:** advection of the *world-particles*, which stream
  on Newtonian arcs and are not on a grid at all (see
  [streaming postulates](#streaming-postulates-s-s)).

### aligned pair

Two sea members, a positon and a negaton, whose clocks agree, so that their
misalignment $`\mu = 0 \pmod{2\pi}`$. Such a pair contributes exactly zero to
the observable $`E`$ (it is *dark*, and $`W`$-null) while remaining a live
participant in the dynamics. "Aligned" is the relation; "sea" says what the
pair is doing in the population; "dark" says it contributes nothing to $`W`$.
Each word is true of an aligned sea pair, and none implies the others.

- **Defined:** [`sea_phase_reference.md` §5](analysis/sea_phase_reference.md#5-definition-l0-the-sea-lock)
  (the relation, step 22); [`what_the_reach_is.md` §1.3](supplement/what_the_reach_is.md#13-three-ways-of-contributing-nothing).
- **Not to confuse with:** a *bound pair*, the older word for the same thing,
  withdrawn because nothing binds the members (step 20, open item Y-SP8); a
  *ket–bra pair*, which is a relation between two worlds about one body (see
  [pair](#pair)); a [gray pair](#gray-pair).

### aperture

The region of phase space one event, or one reader, draws on: width $`dp`$ in
momentum and $`2y_{\max}`$ in separation, with area fixed by
$`B\,dp\,2y_{\max} = 1`$, so that at the physical density it holds **one**
pair on average. A sea of $`\beta`$ pairs per aperture is an *integer sea* of
depth $`\beta`$.

- **Defined:** [`sea_depletion.md` §5](analysis/sea_depletion.md#5-proposition-v4-granular-supply)
  (from Q9(c) of `force_blind_sea.md`).
- **Not to confuse with:** an optical aperture, or the old reading of the
  [reach](#reach) as an aperture, which
  [`what_the_reach_is.md` §3](supplement/what_the_reach_is.md) settles in
  favour of a *period*.

---

## B

### body

An unpaired positon or negaton: one that contributes to the observable. The
body count $`N = n_+ + n_-`$ is what the ledger counts; the observable reads
only $`n_+ - n_-`$.

- **Defined:** [`ORIENTATION.md` §2.3](../ORIENTATION.md#23-the-cast).
- **Not to confuse with:** the multi-particle notes (steps 21 and 24,
  V5, V7, V14), where "two bodies" are two *physical* particles whose joint
  $`W`$ is carried, each by its own world-particles. Context decides.

### body inflation

The growth of the body count $`N`$ beyond $`\lVert W\rVert_1`$: bodies of
opposite sign in one cell that could be removed without changing $`W`$. Each
body beyond the minimal ensemble is half a pair of sea debt. It also drains
the apertures a parent keeps returning to, so blocking falls as about
$`0.9\,e^{-0.55\beta}`$ rather than $`e^{-\beta}`$ (Proposition V6). The
[contact sink](#contact-sink-garbage-collection) removes it.

- **Defined:** [`sea_depletion.md` §3](analysis/sea_depletion.md#3-proposition-v2-the-seas-debt-is-the-negativity-made),
  Proposition V2 (ii); in four dimensions
  [`fourd_compensated_ledger.md` §8](analysis/fourd_compensated_ledger.md#8-corollary-m10-what-the-cost-scales-with).

---

## C

### carrier lock

Postulate (S) of **step 6**: sea particles sharing a cell and a momentum row
share a transported phase, up to the pumped misalignment. It is the
sea-wide strengthening of the per-pair darkness condition. Later notes treat
it as withdrawn or unnecessary
([`dark_sea_and_worldline_identity.md` §9](analysis/dark_sea_and_worldline_identity.md#9-three-collisions-of-notation);
[`permanent_pairing_density_matrix.md` §7](analysis/permanent_pairing_density_matrix.md#7-consequences-and-open-items)).

- **Defined:** [`relational_pairing_and_carrier_lock.md` §4](analysis/relational_pairing_and_carrier_lock.md#4-postulate-s-and-the-equivalence-theorem).
- **Not to confuse with:** the other postulate called (S), streaming
  ([streaming postulates](#streaming-postulates-s-s)); the *sea lock* of step 22
  ([sea lock](#sea-lock)); a *carrier* in the sense of
  [world-particle](#world-particle-world-carrier).

### chord

Two aligned pairs on one momentum row, at positions $`x - y`$ and $`x + y`$.
The pair of them has a misalignment $`\mu^{(P)} = \theta_i - \theta_j + Pd/\hbar`$
at any momentum $`P`$, and a chord's contribution to a [reading](#reading-reader) at
transfer $`\xi`$ is $`\propto U_{\rm res}\sin(2\xi y/\hbar + \mu)`$, with $`U_{\rm res}`$ the
potential difference across the chord, $`V(x+y) - V(x-y)`$, less its part linear
in $`y`$. A chord is a relation
between two sea pairs, not an object: it exists whenever two aligned pairs share
a row, and a body at the chord's midpoint can read it.

- **Defined:** [`supplement/chord_reading.md` §3](supplement/chord_reading.md#3-the-chord).
- **Not to confuse with:** a chord of a circle or a graph; the *separation*
  $`2y`$ of a ket–bra pair, a path difference between two branches
  (see [separation](#separation)).

### coboundary

A difference of a one-index field. The pair misalignments
$`\mu_{ij} = \Phi_i - \Phi_j`$ are the coboundary of the clock field $`\Phi`$:
they satisfy $`\mu_{ij} = -\mu_{ji}`$ and
$`\mu_{ij} + \mu_{jk} = \mu_{ik} \pmod{2\pi}`$. So a stored partnership index
carries nothing that the clocks do not already say (Proposition R1).

- **Defined:** [`relational_pairing_and_carrier_lock.md` §2](analysis/relational_pairing_and_carrier_lock.md#2-partnership-cannot-be-a-carrier-of-state),
  Proposition R1.
- **Not to confuse with:** the term of algebraic topology beyond this one
  usage; nothing more is meant than that a two-index quantity is a difference
  of a one-index one.

### coherence horizon

Definition (H): a ket–bra pair whose separation exceeds $`L_c`$ is not
instantiated. It is a maximum rung in the *relative* coordinate, not a wall
in position, and half of it is the [reach](#reach), $`y_{\max} = L_c/2`$.

- **Defined:** [`open_position_space.md` §3](analysis/open_position_space.md#3-the-coherence-horizon).
- **Not to confuse with:** a distance past which coherence is lost. Nothing
  is lost: past the horizon coherences become *indistinguishable*
  ([`what_the_reach_is.md` §2](supplement/what_the_reach_is.md)). Nor an
  *absorbing boundary*: see [absorbing boundary](#absorbing-boundary).

### commensurate window

A sharp cut-off of the ket–bra separation, of reach $`y_{\max}`$, for a potential
of period $`a`$, is *commensurate* when $`2y_{\max}`$ is an integer multiple of
$`a`$. A commensurate sharp window is exact on the momentum lattice; any other
window, sharp or tapered, is not (Proposition E1.1: commensuration, not
sharpness).

- **Defined:** [`reach_energy_coupling.md` §1](analysis/reach_energy_coupling.md#1-the-reach-is-a-period-not-an-aperture),
  Proposition E1.1.
- **Not to confuse with:** an *aperture*, which tapers; the reach is a period
  (see [reach](#reach)).

### compensated split

The Wigner equation written with the potential term moved to the right,

```math
\partial_t W + \frac{p}{m}\,\partial_x W - V'(x)\,\partial_p W = \mathcal R[W].
```

The left side is exactly classical Liouville transport, with the full classical
force. The right side, the *residual*, is what remains of the potential term
once its part linear in the ket–bra separation (the classical force) has been
subtracted; it carries all the non-classical content, and it is realised by
events rather than forces. "Compensated" means that subtraction.

- **Defined:** [`ORIENTATION.md` §3](../ORIENTATION.md#3-the-one-central-idea-the-compensated-split);
  [`compensated_liouville_splitting.md`](analysis/compensated_liouville_splitting.md).
- **Not to confuse with:** an expansion in $`\hbar`$. The split is exact; it
  is the odd part of the cubic Taylor remainder of $`V`$ that remains, and it
  vanishes identically for a potential of degree two or less.

### contact sink (garbage collection)

The rule that a free positon and a free negaton in one mesh cell bind at once
into a sea pair and credit the sea. It keeps the ensemble minimal, $`N = \lvert W\rvert`$,
and the sea then pays exactly half the growth of $`\lVert W\rVert_1`$ whichever
route the events take. Step 16 called it the *recombination rate*; it is not
part of the postulates.

- **Defined:** [`force_blind_sea.md` §6.2](analysis/force_blind_sea.md#62-proposition-q8-immediate-contact-recombination-as-garbage-collection),
  Proposition Q8.
- **Not to confuse with:** *recombination* in the model, which is the
  absorptive realisation of one residual event at its parent's row
  ([recombination](#recombination)). Whether the sink makes
  [gray pairs](#gray-pair) when the clocks disagree is open (Q-SP10).

### coset invariant

Theorem O5: for an $`a`$-periodic potential every world's residue
$`\theta = p \bmod (\pi\hbar/a)`$ is constant for all time, on the whole line,
with no box and no horizon. Streaming changes only $`x`$, and every vertex
changes a momentum by an integer multiple of $`\pi\hbar/a`$. So the momentum
rows fall into sectors that never mix.

- **Defined:** [`open_position_space.md` §4.2](analysis/open_position_space.md#42-the-coset-invariant).
- **Not to confuse with:** a coset in group theory beyond this: the residue
  classes of $`p`$ modulo $`\pi\hbar/a`$.

### critical sampling

The statistical regime of a reading at the physical density. A reader's
aperture holds on average $`\nu^2 B\,dp\,\Delta x`$ chords, where $`\nu`$ is
the sea's multiplicity (density in units of the physical one); at
$`\nu = 1`$ that is about $`0.08`$. The reading is unbiased, but one reading
is a sample, not an average: a body cannot drive its own clock from its own
reading, and the estimate converges only when many worlds are read at once.

- **Defined:** [`force_blind_sea.md` §9.4](analysis/force_blind_sea.md#94-proposition-q12-the-live-reading-as-the-clocks-input),
  Q12(f); named in Proposition Q5; plainly in
  [`chord_reading.md` §4](supplement/chord_reading.md#4-the-reading).
- **Open:** whether a reading pooled over time or over worlds repairs it
  (Q-SP13).

---

## D

### dark

Present, streaming, and contributing exactly zero to $`W`$: by cancellation in
$`E`$, or by interference at $`\mu \equiv 0`$. A dark pair has a definite
momentum and a Newtonian worldline. It is neither absent nor destroyed.

|  | momentum | moves? | contributes to $`E`$ | verdict |
|---|---|---|---|---|
| self-conjugate carrier | none defined | no flux | the full probability density | **still** |
| dark sea pair | definite $`p`$, on the parent row | yes | exactly zero | **invisible** |
| free positon or negaton | definite $`p`$ | yes | $`\pm 1`$ | neither |

- **Defined:** [`what_the_reach_is.md` §1.3](supplement/what_the_reach_is.md#13-three-ways-of-contributing-nothing).
- **Not to confuse with:** *sea* (what the pair does in the population); a
  *self-conjugate carrier* (still, not dark; see [leg](#leg-arm)); dark matter.

### dark catalysis

Theorem L5: an aligned pair $`A`$ at $`(x,p)`$ fires channel $`q`$ at the
kernel's own rate, $`\lvert K_q(x)\rvert`$, with both members together. Its
positon member is realised as ionisation of a partner pair $`C`$ in $`A`$'s
row, its negaton member as recombination of exactly those two bodies back to
$`p`$. Then $`\Delta E = \Delta N = \Delta S = 0`$ in every realisation, and
$`C`$ leaves aligned. It moves phase and nothing else: the sea's one
$`E`$-neutral interaction.

- **Defined:** [`sea_phase_reference.md` §6](analysis/sea_phase_reference.md#6-theorem-l5-dark-catalysis-and-maintenance);
  its role as a reset clock in
  [`force_blind_sea.md` §8](analysis/force_blind_sea.md#8-proposition-q6-dark-catalysis-as-a-reset-clock).
- **Not to confuse with:** [dark](#dark) as a property; *catalysed
  recombination* (a different thing: [recombination](#recombination)).

### debt

The sea's global shortfall: $`D(t) = S_{\rm tot}(0) - S_{\rm tot}(t)`$ pairs.
From a minimal preparation $`D(t) \ge M_-(t) - M_-(0)`$, the growth of the
negative mass of $`W`$, with equality exactly when the ensemble stays minimal
(Proposition V2). Debt is global; the [depth](#depth) is local.

- **Defined:** [`sea_depletion.md` §3](analysis/sea_depletion.md#3-proposition-v2-the-seas-debt-is-the-negativity-made).
- **Not to confuse with:** the worst-cell deficit (depth), or a financial debt.
  A *cell in debt* is a cell whose sea density has gone negative in a run
  without a floor; the collision term then uses $`\max(s,0)`$
  ([`sea_depletion.md` §8.1](analysis/sea_depletion.md#81-the-event-step-and-collisions-on-cells-in-debt)).

### defocus, focus

**Origin: Cyganski's four actions.** Focus and Defocus are two of the four
collision rules of David Cyganski's closed four-action model (step 2), written
for occupancies $`W_{k-1}, W_k, W_{k+1}`$ at three adjacent momentum cells in
one position column. **Focus** moves one particle from each outer cell into
the middle cell, $`(W_{k-1}, W_k, W_{k+1}) \to (-1, +2, -1)`$; **Defocus** is
its reverse, $`(+1, -2, +1)`$. They are the two-body, zero-photon collisions:
momentum is conserved among the participants, and no photon is exchanged. The
other two actions, Right-Hop and Left-Hop, are one-body and one-photon. See
[four actions](#four-actions).

**Current use: pairs in the sea (step 24, Definition (C)).** The same two moves
are borrowed for *aligned pairs* on neighbouring rows. **Defocus:** two pairs
on row $`n`$ become one on row $`n - q`$ and one on row $`n + q`$. **Focus:**
the inverse. Mass-action rates in detailed balance,
$`J_{n,q} = r_q(s_n^2 - s_{n-q}s_{n+q})`$. They conserve pairs and pair
momentum and do not touch $`W`$. In step 2 the moves redistribute the tracked
particles and, with the hops, reproduce the QLE; in step 24 they act on pairs
that are $`W`$-null, and their job is to relax the sea, not to make quantum
dynamics.

- **Defined:** Cyganski's note (the project file *Bohm to Nelson and four
  action Wigner*, 3 August 2026), read in
  [`four_action_foundations.md` §8](supplement/four_action_foundations.md#8-why-these-four-actions)
  and [`four_rule_microdynamics_equivalence.md`](analysis/four_rule_microdynamics_equivalence.md)
  (step 2); the sea use in
  [`sea_depletion.md` §4.1](analysis/sea_depletion.md#41-pair-collisions-an-amended-s),
  Definition (C).
- **Not to confuse with:** *focus* and *defocus* of an optical beam.

### depth

How deep the sea must be for the [floor](#floor) never to bind. A sea of depth
$`\beta`$ has density $`\beta B`$ everywhere, and $`\lambda^*`$ is the deepest
deficit the unfloored run makes anywhere, as a multiple of $`B`$. $`E`$ is the
same bit for bit when $`\beta \ge \lambda^*`$ and differs below it. The physical
density is $`\beta = 1`$.

- **Defined:** [`sea_depletion.md` §4](analysis/sea_depletion.md#4-proposition-v3-the-floor),
  Proposition V3.
- **Not to confuse with:** the [debt](#debt), which is the sea's global
  shortfall summed over all cells (Proposition V2), not its worst cell; or
  depth in a potential well.

### detailed balance

For each transition, the rate of the forward process times the equilibrium
occupation equals that of its reverse. Definition (C) has it, so a uniform sea
is its equilibrium, and with it the collision entropy
$`\mathcal H = \int[s\ln(s/B) - s + B]`$ never increases (an *H-theorem*; see
[H-theorem](#h-theorem)).

- **Defined:** [`sea_depletion.md` §4.3](analysis/sea_depletion.md#43-proposition-v10-collisions-relax-the-sea),
  Proposition V10.
- **Source:** Boltzmann, *Sitzungsber. Akad. Wiss. Wien* **66** (1872) 275;
  Horn and Jackson, *Arch. Rational Mech. Anal.* **47** (1972) 81.

---

## E

### eikonal

The phase a wave gains from a potential along a path. The word has two senses
and the notes use both.

- **Optical (Hamilton–Jacobi).** The full phase function $`S`$ of geometrical
  optics, $`\lvert\nabla S\rvert^2 = n^2`$; for a particle, Hamilton's principal
  function along the true, *curved* path. This is what a body's clock carries
  under postulate (S).
- **Straight-line (Glauber).** The phase a plane wave gains from the potential
  along the *unbent* line, $`-\hbar^{-1}\int V\,dt`$. In scattering it is an
  approximation, valid when $`E \gg \lvert V\rvert`$. For a sea pair under (S′)
  it is **exact**, by postulate, because the pair's path does not bend.

So "a sea clock keeps its row's plane wave, up to the eikonal" means: a sea pair
that crosses a potential comes out with its clock shifted by
$`-\hbar^{-1}\int V\,dt`$ and otherwise stays on its row's wave.

- **Defined:** [`force_blind_sea.md` §4](analysis/force_blind_sea.md#4-proposition-q2-the-kinematics-of-s-and-s)
  ("Which eikonal"), Proposition Q2. Also used in Q6 and in
  [`chord_reading.md` §2](supplement/chord_reading.md#2-postulate-s-in-two-lines).
- **Sources:** Bruns, *Das Eikonal*, Abh. Kgl. Sächs. Ges. Wiss., Math.-Phys.
  Classe **21** (1895), for the optical sense; Glauber, *High-energy collision
  theory*, in *Lectures in Theoretical Physics* **1**, Interscience (1959)
  315–414, for the straight-line sense (with Molière 1947, as cited in
  `force_blind_sea.md`).

### emissive

See [absorptive, emissive](#absorptive-emissive).

### ensemble (E1, E2)

A population of world-particles, and a rule that says what it represents. Step
13 separates two. **E1** draws carriers from $`W`$, with species the sign of
$`W`$. **E2** draws them from the density matrix, with positon and negaton
naming the ket and bra legs of a pair. Theorem D0: the Weyl transform relates
the *represented objects*, not the ensembles; the species censuses are
anti-correlated and no carrier-level map exists.

- **Defined:** [`species_sectors_and_annihilation.md`](analysis/species_sectors_and_annihilation.md).
- **Not to confuse with:** a statistical ensemble in general. See
  [positon, negaton](#positon-negaton).

---

## F

### floor

Definition (F), a hypothesis tested in step 24 and not adopted: an emission must
ionise a pair in its parent's cell, and if the sea there is short the emission
does not happen and its two depositions are lost. It is step 16's Proposition S0
(co-location is forced) read literally.

- **Defined:** [`sea_depletion.md` §4](analysis/sea_depletion.md#4-proposition-v3-the-floor).
- **Not to confuse with:** a lower bound on something, or the ground state; here
  the sea is a ledger and the floor is the rule that it may not go negative.
  See [depth](#depth).

### four actions

Focus, Defocus, Right-Hop and Left-Hop: the four collision rules of David
Cyganski's closed four-action model, the starting point of step 2. At one
position, on three adjacent momentum cells $`k-1, k, k+1`$:

| action | occupancy change | kind |
|---|---|---|
| Focus | $`(-1, +2, -1)`$ | two-body, zero photon |
| Defocus | $`(+1, -2, +1)`$ | two-body, zero photon |
| Right-Hop | $`(-1, 0, +1)`$ | one-body, one photon |
| Left-Hop | $`(+1, 0, -1)`$ | one-body, one photon |

Given the photon quantum and the triple of cells, this is the complete list of
momentum-conserving or one-photon events that keep particle number
([`four_action_foundations.md` §8](supplement/four_action_foundations.md#8-why-these-four-actions)).
For a potential with a single Fourier mode they are exactly the four available
momentum transfers (step 12, Theorem I3), so "four channels" and "four actions"
name the same thing. Step 24 borrows two of them for the sea
([defocus, focus](#defocus-focus)); [`ORIENTATION.md` §10](../ORIENTATION.md#10-what-the-compensated-model-keeps-from-the-earlier-layers)
says what the compensated model keeps of the collision layer.

- **Defined:** Cyganski's note of 3 August 2026, read in
  [`four_action_foundations.md`](supplement/four_action_foundations.md);
  the rules and their exactness in
  [`four_rule_microdynamics_equivalence.md`](analysis/four_rule_microdynamics_equivalence.md);
  realised as collisions against a pinned sea in
  [`sea_dressed_microdynamics.md`](analysis/sea_dressed_microdynamics.md).
- **Not to confuse with:** the four actions of any other model. The name
  "FAWPM" in the notes means this one.

---

## G

### garbage collection

See [contact sink](#contact-sink-garbage-collection).

### gray pair

A pair whose members' clocks are not aligned, $`\mu \ne 0`$. Its grayness is
$`\sin^2(\mu/2)`$, from 0 (dark) to 1 (fully visible). A gray pair does not
cancel in $`W`$. Pair collisions that exchange partners make gray pairs; collisions
that move each pair whole do not.

- **Defined:** open item Q-SP10 ([`force_blind_sea.md` §12](analysis/force_blind_sea.md#12-open-items));
  measured in [`sea_depletion.md` §4.5](analysis/sea_depletion.md#45-proposition-v12-the-clock-at-a-collision),
  Proposition V12.
- **Not to confuse with:** *gray* as an intermediate colour; here it names a
  phase relation, not a population.

---

## H

### H-theorem

Boltzmann's theorem that a collision entropy never increases. Here: the
functional $`\mathcal H = \int[s\ln(s/B) - s + B]`$ of the sea never increases
under Definition (C), and is constant only when $`\ln s`$ is linear along every
row, which on a periodic row means uniform (Proposition V10). Contrast Proposition V9: a *reversible* sea dynamics conserves a
depletion functional and so cannot repair depletion.

- **Defined:** [`sea_depletion.md` §4.3](analysis/sea_depletion.md#43-proposition-v10-collisions-relax-the-sea).
  See [detailed balance](#detailed-balance) and [Volterra lattice](#volterra-lattice).

---

## I

### integer sea

The sea as it is, a random configuration of whole aligned pairs, rather than a
density. At the physical density an aperture holds one pair on average, so an
emission is blocked by a floor with probability about $`e^{-\beta}`$ when the
aperture holds $`\beta`$ pairs (Proposition V4). With world-particle dynamics
(body inflation) blocking falls as about $`0.9\,e^{-0.55\beta}`$ (V6); with
collisions, about $`1.1\,e^{-0.74\beta}`$ (V13).

- **Defined:** [`sea_depletion.md` §5](analysis/sea_depletion.md#5-proposition-v4-granular-supply).
- **Not to confuse with:** the mean-field sea ([mean field](#mean-field)),
  which is a density.

### integrate-and-fire

A deterministic clock for an event rate (Proposition Q11). Each channel of
each body carries an integrator, with a hidden initial phase, fed by the
channel's rate; when it passes a unit threshold the channel fires. The count follows the rate with error bounded by one,
instead of growing as the square root, so the clock is *quieter* than a
Poisson one. A firing is one fringe of the daughter row's interference pattern
passing the reader.

- **Defined:** [`force_blind_sea.md` §9.3](analysis/force_blind_sea.md#93-proposition-q11-the-resonance-condition-and-its-clock);
  tested on an exactly solvable problem in
  [`poisson_kicks_and_pair_branching.md`](supplement/poisson_kicks_and_pair_branching.md)
  (Proposition X7).
- **Borrowed from:** the neuron model of that name. In signal processing it is
  a sigma–delta modulator; the Poisson clock is the same construction with
  independent exponential thresholds. What it costs is ontology: every body
  carries hidden initial phases, accepted provisionally (Q-SP14).

### interworld coupling

The potential difference a world consults: $`U = V(x_1) - V(x_2)`$ between the
two legs $`x_{1,2} = x \pm y`$ of a position pair (step 12). It vanishes at
coincidence, is antisymmetric under leg exchange, and vanishes identically for
a free particle (Proposition I1), so free wavepacket spreading is not caused by
any force one world exerts on another. The part of $`U`$ linear in the
separation is exactly the classical force (Theorem I4). The coupling exerts
no force and does no work.

- **Defined:** [`interworld_coupling.md`](analysis/interworld_coupling.md).
- **Not to confuse with:** the interworld potential of many-interacting-worlds
  models, which is what makes a free packet spread there. This one does not.

### ionisation

A sea pair splits into two bodies: the positon appears at $`p + \xi_q`$ and
the negaton at $`p - \xi_q`$, both from the parent's row. It is the
microscopic realisation of an emissive event, and it is why the sea is not
optional bookkeeping: the daughters' momentum $`2p`$ has to come from
something already at $`p`$.

- **Defined:** [`ORIENTATION.md` §4](../ORIENTATION.md#4-what-the-residual-events-actually-are).
  See [absorptive, emissive](#absorptive-emissive).
- **Not to confuse with:** ionisation of an atom. No charge is involved.

---

## K

### kernel

The compensated residual kernel $`K_q`$: what is left of the QLE's potential term once
the classical force is subtracted. $`K_q(x)`$ is the signed rate, per body at
$`x`$, for a momentum transfer $`\xi_q = q\,dp`$; its modulus is the event rate.
It is odd in $`q`$, vanishes identically for a potential of degree two or less,
and is a Fourier sum over chord separations.

- **Defined:** [`ORIENTATION.md` §§3–4](../ORIENTATION.md#3-the-one-central-idea-the-compensated-split);
  the integral form in
  [`chord_reading.md` §5](supplement/chord_reading.md#5-what-the-sum-over-chords-estimates).
- **Not to confuse with:** a kernel in the sense of probability or an integral
  transform in general; here it is always the residual of the potential term.

---

## L

### ladder, step, rung

The **ladder** is the numbered sequence of notes in `analysis/`; "step $`n`$"
is the $`n`$th of them, and each takes as input what its predecessor
postulated. Two of the steps are themselves called ladders: the **coherence
ladder** (step 8), which indexes $`\rho`$ by splitting rung, and the
**position-pair ladder** (step 9), a lattice of ket–bra pairs. A **rung** is a
step of whatever lattice a note is built on: an index of the ket–bra splitting
or separation in steps 8 and 9, and in the collision-layer notes also a
momentum value, which is the same as a [row](#row).

- **Defined:** [`ORIENTATION.md` §2](../ORIENTATION.md#2-how-the-repository-is-organised-and-the-words-it-uses);
  [`coherence_ladder.md`](analysis/coherence_ladder.md),
  [`position_pair_ladder.md`](analysis/position_pair_ladder.md).
- **Not to confuse with:** the *ladder* $`\mu`$, which belongs to the
  position-pair ladder of step 9, not to the ladder of notes. See
  [μ](#mu-the-two).

### ledger

The census of world-particles, per cell: $`n_+`$ positon bodies, $`n_-`$
negaton bodies, and $`S`$ sea pairs, with body count $`N = n_+ + n_-`$. The
observable $`W \propto n_+ - n_-`$; the QLE says nothing about $`n_+ + n_-`$ or
$`S`$, and the ledger is everything the observable does not read. Pairs are
split and combined, never made from nothing, so $`S + n_+`$ and $`S + n_-`$
are conserved in exact integers (Theorem N1).

- **Defined:** [`ORIENTATION.md` §8](../ORIENTATION.md#8-the-ledger-and-why-it-is-invisible);
  [`stochastic_ledger.md`](analysis/stochastic_ledger.md).
- **Not to confuse with:** a *minimal ledger*, which is a ledger run on a
  minimal ensemble ([minimal ensemble](#minimal-ensemble)).

### leg, arm

*Leg* has three senses in the notes, and context decides.

1. **A leg of a ket–bra pair** (steps 7 to 13 and the current notes): one end
   of the pair, a place and a clock. The two positions $`x \pm y`$ are two
   *branches'* answers to where the one particle is, not two particles. In
   ensemble E2 the positon is the ket leg and the negaton the bra leg. A
   **self-conjugate carrier** is a carrier with coincident legs: it holds the
   full probability density and no current, which is *still*, not the same as
   [dark](#dark).
2. **A leg of a worldline** (step 4): a maximal constant-momentum segment of a
   worldline, ended at each end by a [vertex](#vertex). A worldline is a chain
   of legs, and a pair is two worldlines, one per species, each currently on
   some leg. The phase-alignment specification's "both legs at row $`r'`$" is
   in this sense.
3. **The legs of one event** (`emission_and_absorption.md` §3): the two
   co-located deposits of one residual event, $`+1`$ at $`p + \xi_q`$ and
   $`-1`$ at $`p - \xi_q`$. They are not two events.

*Arm* is older and worse than *leg*, since it suggests a limb.

- **Defined:** (1) [`what_the_reach_is.md` §1](supplement/what_the_reach_is.md#1-first-what-is-separated);
  (2) [`phase_resonance_microdynamics.md` §3](analysis/phase_resonance_microdynamics.md#3-kinematics-pairs-beats-and-the-parity-theorem);
  (3) [`emission_and_absorption.md` §3](supplement/emission_and_absorption.md#3-what-an-event-is).
- **Not to confuse with:** a leg of a journey or a graph; the *members* of a
  sea pair (two world-particles, not two branches).

---

## M

### mean field

A description by densities instead of counts. The ledger notes run $`u^\pm(x,p)`$
and $`S(x,p)`$ as real-valued fields on a mesh, updated by a tau-leap (see
[tau-leap](#tau-leap)); that is the right tool for fidelity and the wrong one for
questions about integers, births and deaths, and fluctuations. In step 2 the
*mean-field generator* of a rule is the expected rate of change of the
occupancies under it.

- **Defined:** [`stochastic_ledger.md` §0](analysis/stochastic_ledger.md#0-what-this-note-settles-and-what-it-corrects);
  [`four_rule_microdynamics_equivalence.md` §3](analysis/four_rule_microdynamics_equivalence.md#3-mean-field-generator-of-a-focushop-scheme).
- **Not to confuse with:** mean-field theory in condensed matter. Here it means
  only "density in place of counts"; it is the opposite of an
  [integer sea](#integer-sea).

### microdynamics

The particle-level rules (events and their rates) whose ensemble average
reproduces a given equation of motion. "Four-rule microdynamics", "sea-dressed
microdynamics", "phase-alignment microdynamics" are successive proposals for
the QLE; a **gauge freedom of the microdynamics** is a family of different
rate assignments with the same average.

- **Defined:** [`ORIENTATION.md` §2.1](../ORIENTATION.md#21-a-ladder-in-three-layers)
  (the collision-microdynamics layer, steps 2–10).

### minimal ensemble

A world-particle ensemble with $`N = \lVert W\rVert_1`$: no cell holds both a
positon and a negaton body, so nothing can be removed without changing $`W`$.
Its debt is exactly the growth of negative mass, $`\Delta\delta/2`$ in terms of
Kenfack and Życzkowski's negativity indicator $`\delta = \lVert W\rVert_1 - 1`$.
The [contact sink](#contact-sink-garbage-collection) keeps an ensemble minimal.

- **Defined:** [`sea_depletion.md` §3](analysis/sea_depletion.md#3-proposition-v2-the-seas-debt-is-the-negativity-made),
  Proposition V2.

### momentum, in the ladder notes

In the position-pair ladder (step 9), momentum is not a quantity a body
carries: it is the misalignment $`\bar p = \hbar\mu/a`$ of a pair of legs,
with $`\mu`$ the ladder $`\mu`$. In the world-particle notes, by contrast, a
world has a momentum as a coordinate, and the sea's rows are values of it.

- **Defined:** [`position_pair_ladder.md`](analysis/position_pair_ladder.md).
- **Not to confuse with:** the pair $`\mu`$, which gives the *splitting*
  $`\Delta p`$ of two members, not a total.

### Moyal equation, QLE

See [QLE](#qle-quantum-liouville-equation).

### $`\mu`$, the two

Two different quantities are written $`\mu`$, and conflating them produces a
picture that looks coherent and is wrong.

| | **Pair $`\mu`$** (step 5) | **Ladder $`\mu`$** (step 9) |
|---|---|---|
| Definition | $`\mu = \Phi_a - \Phi_b \pmod{2\pi}`$, the misalignment of the transported clock phases of a pair's two members $`a`$, $`b`$ | $`\mu = \arg\rho(X,X')`$, the phase of the density matrix between two leg positions $`X`$, $`X'`$ |
| Between | two **members** of a sea pair: two world-particles, two clocks | two **legs** of a ket–bra pair: one body, two branches |
| Lives in | the phase-resonance and phase-alignment formulation, steps 4–8 | the position-pair ladder, step 9 |
| Winds at | $`\partial_x\mu = \Delta p/\hbar`$, the *splitting* $`\Delta p = p_a - p_b`$; along a path of velocity $`v`$, $`d\mu/dt = \Delta p\,(v - \bar v_{\rm pair})/\hbar`$ | $`\partial\mu/\partial X = p/\hbar`$, the *total* momentum; $`\bar p = \hbar\mu/a`$ with $`a`$ the position spacing |
| $`\mu = 0`$ means | $`\lvert\Psi\rvert = 2\lvert\sin(\mu/2)\rvert = 0`$: **dark** | $`p = 0`$: no current |
| Conjugate to | — | the separation $`y`$ the compensated kernel integrates over |

Both are gauge-invariant relative phases of a two-ended object, whose gradient
is the momentum conjugate to the ends' separation, and whose vanishing is a
kind of invisibility: same structure, different variable. If a summary runs
"reach $`\to \Delta p \to \mu`$-winding $`\to`$ dark sea $`\to`$ residual
ionisation" as one chain, it has crossed between the two.

- **Defined:** [`what_the_reach_is.md` §1.2](supplement/what_the_reach_is.md#12-two-quantities-are-called-mu);
  pair $`\mu`$ in [`phase_alignment_microdynamics.md`](analysis/phase_alignment_microdynamics.md),
  ladder $`\mu`$ in [`position_pair_ladder.md`](analysis/position_pair_ladder.md).
- **Open:** whether the phase machinery of steps 4–8, built on a vertex in
  which the excess particle changes momentum, transplants onto the compensated
  model's mediated-jump vertex is unresolved (J-SP4 of
  [`emission_and_absorption.md`](supplement/emission_and_absorption.md)).

---

## O

### observable map

The map $`\mathcal E`$ from an ensemble to $`W`$. It reads $`n_+ - n_-`$ and
nothing else. What it cannot see (the body count $`N`$, the sea $`S`$, which
realisation of an event happened) lies in its kernel, and that is where the
[ledger](#ledger) lives: the two realisations of an event move the observable
identically (Theorem N2).

- **Defined:** [`ORIENTATION.md` §8](../ORIENTATION.md#8-the-ledger-and-why-it-is-invisible).
- **Not to confuse with:** an observable in the quantum-mechanical sense, an
  operator; and with the *kernel* of [kernel](#kernel), which is the residual
  kernel $`K_q`$. "In the kernel of the observable map" is the null space.

---

## P

### pair

**Two different things in this project.** A *ket–bra pair* is a relation
between two worlds about one body, with separation $`2y`$. A *sea pair* is two
world-particles, a positon and a negaton, aligned or not. Different in kind;
the older notes use the bare word for either.

- **Defined:** the ket–bra pair in
  [`ORIENTATION.md` §3](../ORIENTATION.md#3-the-one-central-idea-the-compensated-split);
  the sea pair in [`ORIENTATION.md` §2.3](../ORIENTATION.md#23-the-cast) and
  at [aligned pair](#aligned-pair).
- **Not to confuse with:** two bodies. Neither kind is.

### particle

In the compound *world-particle*, a signed sample of $`W`$. Bare "particle" in
a single-body note almost always means *world*, and where it means *leg* it is
doubly wrong. See [world-particle](#world-particle-world-carrier).

### population

The world-particle count $`N`$ of the ledger notes. In density-matrix usage it
is also the standard name for a diagonal element of $`\rho`$, in the sense of
Bloch and Fano. The notes write the second out rather than use the word.

- **Defined:** see [ledger](#ledger).

### positon, negaton

In the ledger and sea notes: a world-particle of sign $`+1`$
(positon) or $`-1`$ (negaton), drawn from the signed $`W`$. The observable is
$`n_+ - n_-`$.

- **Not to confuse with:** the second sense, from step 13's density-matrix
  ensemble E2, in which positon and negaton name the *ket and bra legs of a pair*
  (see [ensemble](#ensemble-e1-e2)). The Weyl transform relates the two ensembles' represented objects, not the
  ensembles (Theorem D0), so the word is a homonym across the two layers.

---

## Q

### QLE (quantum Liouville equation)

The equation of motion of the Wigner function, the Wigner transform of the
von Neumann equation; also called the **Moyal equation**, since its potential
term is the Moyal bracket's. It is the reference theory that the model
reproduces exactly in expectation. Its **Moyal series** is the expansion of
that potential term in powers of $`\hbar`$; the
compensated split is not that expansion (see
[compensated split](#compensated-split)).

- **Defined:** [`ORIENTATION.md` §3](../ORIENTATION.md#3-the-one-central-idea-the-compensated-split);
  the Moyal series matched to the lattice in
  [`phase_space_crystal_lattice_supplement.md` §3](supplement/phase_space_crystal_lattice_supplement.md#3-matching-xfp-to-the-moyal-series-of-the-qle).
- **Source:** Moyal, *Proc. Cambridge Philos. Soc.* **45** (1949) 99.

---

## R

### reach

The greatest half ket–bra separation instantiated, $`y_{\max}`$. Equivalently the period of the
separation: $`\Delta p = \pi\hbar/2y_{\max}`$ is the spacing of the momentum
rows. $`y_{\max} = L_c/2`$, half the [coherence horizon](#coherence-horizon).

- **Defined:** [`ORIENTATION.md` §6](../ORIENTATION.md#6-the-reach) and
  [`what_the_reach_is.md`](supplement/what_the_reach_is.md).
- **Not to confuse with:** how far a world can see or act, or an aperture. That
  wording was retracted (Theorem E1 of step 11b) and survives in Definition (R)
  only with a superseding note. The **reach ceiling** is two different
  numbers, and the tighter one is not a bound on the theory
  ([`what_the_reach_is.md` §5](supplement/what_the_reach_is.md)).

### reading, reader

A **reader** is a free body. Its **reading** is the sum over the chords of the
sea around it, which estimates the kernel $`K_q`$ at its position: the sea
computes the Fourier transform of the kernel by interference (Theorem L1,
Proposition Q5). A reading is a rate estimate, not a measurement outcome.

- **Defined:** [`chord_reading.md` §4](supplement/chord_reading.md#4-the-reading).
- **Not to confuse with:** a measurement result; a reading of the QLE as a
  statement about what the sea *is*. What it estimates is the *generator*, not
  the state ([`chord_reading.md` §5](supplement/chord_reading.md#5-what-the-sum-over-chords-estimates)).

### recombination

Two senses, both in the notes.

- **Catalysed recombination:** the absorptive realisation of one residual event
  (see [absorptive, emissive](#absorptive-emissive)): a body of each sign
  combines into a sea pair at its *parent's* row.
- **Contact recombination:** two bodies in one cell combining with no
  parent: the sink of the ledger. The collision layer (steps 3, 13, 16 and 18)
  carried it as a *recombination channel* or $`\kappa`$ *sink*. It is not part
  of the postulates ([contact sink](#contact-sink-garbage-collection)).

Step 20 separates the two as *catalysed* and *contact*.

- **Defined:** [`ORIENTATION.md` §§4, 8](../ORIENTATION.md#4-what-the-residual-events-actually-are).

### re-locking

The rule that a body or pair, on joining a row, sets its clock to that row's
plane wave: to the circular mean of the row's sea clocks within a window, or,
for a whole pair that collides into a new row, to the new row's wave. On
recombination both members take the parent row's value. Only whole pairs
re-locked to their new row keep the sea dark after a collision (Proposition V12).

- **Defined:** [`sea_phase_reference.md` §6](analysis/sea_phase_reference.md#6-theorem-l5-dark-catalysis-and-maintenance);
  at a collision, [`sea_depletion.md` §4.5](analysis/sea_depletion.md#45-proposition-v12-the-clock-at-a-collision).
- **Borrowed from:** injection locking (Adler 1946), to which the notes
  compare it.

### residual, residual event

The right side $`\mathcal R[W]`$ of the [compensated split](#compensated-split),
and the events that realise it. A **residual channel** $`q`$ transfers
momentum $`\xi_q`$ with signed weight $`K_{\rm res}(q)`$. The kernel is real
and odd, hence signed, so it cannot be a one-body Markov jump generator: the
negative part is supplied by the sea.

- **Defined:** [`ORIENTATION.md` §4](../ORIENTATION.md#4-what-the-residual-events-actually-are).
- **Not to confuse with:** a numerical residual or error.

### row

One value of the discrete momentum grid, $`p = n\,dp`$ with
$`dp = \pi\hbar/2y_{\max}`$. A sea pair sits on a row. Under (S′) a pair keeps
its row between events; collisions move pairs between rows. *Row-locked* means
a pair's clock is that of its row's plane wave. The older notes call a row a
*rung*.

- **Defined:** [`force_blind_sea.md` §4](analysis/force_blind_sea.md#4-proposition-q2-the-kinematics-of-s-and-s).
- **Not to confuse with:** a row of a matrix, or a row of cells in space.

---

## S

### sea

The background population of neutral sea pairs that fills every cell of phase
space at density $`B = 1/\pi\hbar`$ whether or not $`W`$ is there. It is the
partner population the residual events consume and return, and the
observable $`E`$ never reads it (Proposition Q1).

- **Defined:** [`ORIENTATION.md` §2.3](../ORIENTATION.md#23-the-cast);
  [`sea_depletion.md` §0](analysis/sea_depletion.md#0-what-this-note-asks-settles-and-leaves-open).
- **Not to confuse with:** *dark* (a phase relation, [aligned pair](#aligned-pair)),
  or a Dirac sea.

### sea lock

Definition L0: the sea is *locked* when each aligned pair's clock equals its
row's plane-wave phase at the pair's position, up to a constant per row. Within
a pair the members are aligned; the lock is the further relation *between*
pairs in one row. A **gas** sea has every pair dark but mutually random.

- **Defined:** [`sea_phase_reference.md` §5](analysis/sea_phase_reference.md#5-definition-l0-the-sea-lock).
- **Not to confuse with:** the *carrier lock* of step 6
  ([carrier lock](#carrier-lock)), a different and older postulate.

### self-conjugate carrier

See [leg, arm](#leg-arm).

### separation

The ket–bra half-separation $`y`$ (or $`2y`$). A path difference between two
branches. Nothing is at both ends at once.

- **Defined:** [`ORIENTATION.md` §3](../ORIENTATION.md#3-the-one-central-idea-the-compensated-split).
- **Not to confuse with:** the distance between two objects.

### separatrix

The curve in phase space that divides trajectories which cross a barrier from
those that do not: the classical boundary between transmission and reflection.
In the Eckart-barrier notes tunnelling appears as an imbalance in separatrix
traffic, and bodies cross it in both directions.

- **Defined:** [`eckart_barrier_compensated.md`](analysis/eckart_barrier_compensated.md);
  who crosses, in [`dark_sea_and_worldline_identity.md` §7](analysis/dark_sea_and_worldline_identity.md).
- **Not to confuse with:** a wall; worlds cross it (Theorem Y6).

### sinkless

The case in which no channel other than the residual events removes bodies:
$`R_{\rm sink} = 0`$. Theorem N3, $`\Gamma_{\rm tot}(1-2f) = R_{\rm sink}`$,
then gives the absorptive fraction $`f = 1/2`$. It is a special case, not a
universal law; any sink forces $`f < 1/2`$.

- **Defined:** [`ORIENTATION.md` §8](../ORIENTATION.md#8-the-ledger-and-why-it-is-invisible);
  [`stochastic_ledger.md` §4](analysis/stochastic_ledger.md#4-theorem-n3-the-sum-rule-and-what-it-corrects).

### streaming postulates (S), (S′)

The letter (S) is used for **two** postulates in the notes.

- **(S), streaming (step 17):** between events every world-particle obeys
  $`\dot q = p/m`$, $`\dot p = -\nabla V`$: the full classical force.
- **(S′), force-blind sea (step 23):** free bodies as in (S); each member of an
  aligned pair instead moves inertially, $`\dot p = 0`$, while its clock still
  winds. **Amended in step 24:** an aligned pair also changes momentum, and
  only, in collisions with other aligned pairs ([defocus, focus](#defocus-focus)).
- **(S), carrier lock (step 6):** see [carrier lock](#carrier-lock).

- **Defined:** [`compensated_ontology.md` §1](analysis/compensated_ontology.md#1-four-postulates);
  [`force_blind_sea.md` §1](analysis/force_blind_sea.md#1-postulate-s);
  the double use is recorded in [`INDEX.md` §2](analysis/INDEX.md#2-prefix-registry).
- **Not to confuse with:** Theorems S0–S9 of step 16, whose prefix is also S
  (`INDEX.md` §2).

### sub-Poissonian

A count whose variance is less than its mean: a Fano factor below 1, where a
Poisson process has 1. The [integrate-and-fire](#integrate-and-fire) clock gives
Fano factors of 0.03 to 0.3 in the measured channels, against 1.17 measured for
a Poisson clock at the same rate (Proposition Q11(c)).

- **Defined:** [`force_blind_sea.md` §9.3](analysis/force_blind_sea.md#93-proposition-q11-the-resonance-condition-and-its-clock).

---

## T

### tau-leap

Advancing a stochastic process by many events at once, with event counts
drawn or averaged at rates held fixed over a step $`\Delta t`$ (the
$`\tau`$-leap). The ledger runs on the mesh are *mean-field tau-leaps*. The
exact alternative is an event-by-event simulation (Gillespie), which has no
time-step error and is what `demo_stochastic_ledger.py` runs.

- **Defined:** [`stochastic_ledger.md` §0](analysis/stochastic_ledger.md#0-what-this-note-settles-and-what-it-corrects);
  [`emission_and_absorption.md` §4](supplement/emission_and_absorption.md).
- **Source:** Gillespie, *J. Phys. Chem.* **81** (1977) 2340, for the exact
  algorithm (see `references/bibliography.md`).

### truncation

Cutting off the ket–bra separation at the [coherence horizon](#coherence-horizon).
Signed world number is conserved exactly by it (Proposition O3). For an
unbounded $`V`$ it is constitutive: there is nothing to refine towards, since
for a polynomial $`V`$ there is no jump measure without a reach (Theorem E7).

- **Defined:** [`open_position_space.md` §3](analysis/open_position_space.md#3-the-coherence-horizon);
  Theorem E7 in [`eckart_barrier_compensated.md`](analysis/eckart_barrier_compensated.md).
- **Not to confuse with:** truncation of a series, an approximation to be
  refined away.

---

## V

### vertex

The point at which a world's momentum changes, where the constant-momentum
legs of a worldline join. The word has **two senses**, one in the collision
layer and one in the compensated model, and whether they describe the same
event is open.

- **Collision layer, steps 4 to 8: the swap.** An excess particle sharing a
  cell with a sea partner exchanges momentum with it. Requiring the pair's
  misalignment $`\mu`$ to hold still during the exchange forces it to be a
  *swap* (Theorem 4 of step 5): the particle takes the struck partner's
  momentum, the partner takes its mate's, and the pair exits aligned. Momentum
  and energy are conserved because a swap permutes momenta among worldlines.
  In absorption (K3) the excess particle's in-leg meets the excited partner of
  a beating pair, leaves at a lower momentum, and the pair exits dark; emission
  (K4) is the time reverse. The *vertex midpoint* $`n`$ is the row halfway
  between the two rows exchanged.
- **Open line and compensated model, steps 11 on: the mediated jump.** The
  *vertex rate* is the rate of residual events, read from $`V`$ at the world's
  own $`x`$ (so it needs no boundary, Proposition O4), and each vertex
  changes a momentum by an integer multiple of $`\pi\hbar/a`$ for an
  $`a`$-periodic potential (Theorem O5). In the compensated model the parent
  does not recoil; the event ionises or recombines a pair on the parent's row.
  From step 22, "what a vertex reads" is the sea: the chords around the
  parent, which acts as the reader.
- **Open:** the K3 vertex has the excess particle exit at a different
  momentum, the mediated-jump vertex does not. The two formulations are
  equivalent as generators, but whether the phase machinery of steps 4 to 8
  transplants onto the mediated-jump vertex is unresolved (J-SP4).

- **Defined:** [`phase_alignment_microdynamics.md` §1.5](analysis/phase_alignment_microdynamics.md#15-the-vertex-is-a-momentum-exchange)
  and [`phase_resonance_microdynamics.md` §4](analysis/phase_resonance_microdynamics.md#4-vertices-absorption-emission-resonance)
  for the swap; [`phase_alignment_microdynamics_algorithm.md` §5](algorithm/phase_alignment_microdynamics_algorithm.md#5-the-vertex)
  for the specification; [`open_position_space.md` §3](analysis/open_position_space.md#3-the-coherence-horizon)
  for the vertex rate; the open question in
  [`emission_and_absorption.md`](supplement/emission_and_absorption.md#j-sp4--two-vertices-are-called-absorption).
- **Not to confuse with:** a vertex of a mesh or lattice, a point of the grid.
  See also [leg](#leg-arm), [four actions](#four-actions) and
  [residual event](#residual-residual-event).

### Volterra lattice

A chain of coupled rates of the form
$`\partial_t s_n = \tfrac{s_n}{B}\sum_q c_q (s_{n-q} - s_{n+q})`$: the
Kac–van Moerbeke system. It is reversible and conserves
$`\sum_n s_n`$ and $`\sum_n \ln s_n`$ on every row, hence a depletion
functional, so a Volterra sea cannot repair depletion (Proposition V9). Cyganski's
*copy-left* and *copy-right* actions, applied to whole pairs, are of this form.

- **Defined:** [`sea_depletion.md` §4.2](analysis/sea_depletion.md#42-proposition-v9-a-reversible-sea-cannot-repair-depletion).
- **Source:** Kac and van Moerbeke, *Adv. Math.* **16** (1975) 160.

---

## W

### Wigner function

The phase-space quasi-probability $`W(x,p)`$ of a state with density matrix
$`\rho`$: the Fourier transform of $`\rho(x+y, x-y)`$ in the half-separation
$`y`$. It has the right marginals in $`x`$ and $`p`$ and can be negative. The
notes carry it as the expectation of a signed ensemble of world-particles.

- **Defined:** [`ORIENTATION.md` §1](../ORIENTATION.md#1-what-the-project-is-trying-to-do).
- **Source:** Wigner, *Phys. Rev.* **40** (1932) 749.

### world-particle (world, carrier)

A point in phase space with a sign: a signed atom of a counting measure on
$`(x,p)`$ that has a position *and* a momentum. Also *world*, and, in the
older notes, *carrier*.

- **Defined:** [`ORIENTATION.md` §2.3](../ORIENTATION.md#23-the-cast).
- **Not to confuse with:** a particle with a position and a velocity: the
  world-particle is a signed sample of $`W`$. *Carrier* also has a
  density-matrix sense, a matrix element $`\rho_{nm}`$; it conveys nothing, and
  the word names an index pair.

---

## Suggested by readers

(none yet)
