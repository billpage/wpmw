# Index of results and open items

This index lists every labelled result (theorem, proposition, lemma, corollary, definition, postulate) and every labelled open item in the [analysis notes](README.md), each with a section reference, a one-line paraphrase and a standing. It exists because the notes cite each other by label — Theorem K4, S-SP3, C2 — faster than a reader can hold the labels in mind. For a label it answers three questions: what does it say, where is it, and is it still true.

**Read the note before relying on a result.** The statements here identify a result; they do not carry its hypotheses. The *Standing* column is the reason the index exists: it records where a later note corrected, restricted, superseded or retracted something, which the note itself cannot say.

1. [How to read this index](#1-how-to-read-this-index)
2. [Prefix registry](#2-prefix-registry)
3. [Postulates and definitions](#3-postulates-and-definitions)
4. [Results by note](#4-results-by-note)
5. [No-go and negative results](#5-no-go-and-negative-results)
6. [Corrections and retractions](#6-corrections-and-retractions)
7. [Open items](#7-open-items)
8. [Defects found in demos](#8-defects-found-in-demos)
9. [Cited from outside this folder](#9-cited-from-outside-this-folder)
10. [Maintaining this index](#10-maintaining-this-index)

---

## 1. How to read this index

- **Steps** are the ladder numbers of the [README](README.md). Step 11b, [`reach_energy_coupling.md`](reach_energy_coupling.md), sits between 11 and 12 so that no existing step number changes.
- **Types** are abbreviated Thm, Prop, Lem, Cor, Def, Post. A label such as `C2` is only unambiguous with its step; §2 lists the labels that are reused.
- **ID** links to the heading that contains the statement. **§** is that heading's number in the note.
- **Says** paraphrases the statement in a line. **Standing** is `—` when nothing is recorded, and otherwise says what a later (or the same) note did to it, or which open item it feeds.
- **Open-item status** is what the notes themselves say. `open` means no note claims to have closed it; it does not mean someone has checked that it is still open.
- Symbols are as in the notes. `abs(x)` stands for the modulus, to keep the tables valid.

## 2. Prefix registry

Theorem and open-item labels are shared across the whole of `docs/`. Each note uses its own letter, with the exceptions below. The project's rule is that prefixes do not collide; the reuses below are the exceptions. **Nothing is renamed here.** Cite a reused label with its step, for example `C2 (splitting)`.

| Prefix | Where defined | Note |
|---|---|---|
| A, B | step 10 (A1–A3, B1–B3) | — |
| C | step 8 (C1–C2); step 10 (C1–C4); step 14 (C0–C9) | Three notes. Cite as 'C2 (splitting)' or 'C2 (coherence)'. Step 9 cites step 8's C2 as 'Theorem C2'; later notes (for example steps 15, 17 and 19) cite step 14's. |
| D | step 10 (Theorem D); step 13 (D0–D18) | Two notes. |
| E | step 11b (E1–E8) | Also the observable E = n+ − n− of the ledger notes and the ensembles E1/E2 of step 13, neither a theorem. |
| G | step 17 (G1–G5) | — |
| I | step 12 (I1–I5, I5a) | — |
| K | step 15 (K1–K9) | Also the vertex channels K1–K4 of steps 4 to 8 (K3 absorption, K4 emission), which are not theorems. |
| N | step 18 (N1–N6) | — |
| O | step 11 (O1–O5) | — |
| P | step 9 (P1–P7); postulates P0–P3, P5 in step 4 | Step 4's postulates and step 9's theorems share P1–P5. There is no postulate P4. |
| R | step 6 (R1–R5); supplement/representation_cost_and_annihilation (R2–R6) | Different folders. |
| S | step 16 (S0–S9); postulate (S) in steps 6 and 17 | Two postulates are called (S): the sea carrier lock (step 6) and streaming (step 17). Step 20 §9 records this. |
| U | step 13 (U1) | — |
| Y | step 20 (Y1–Y6) | — |
| M | step 21 (M0–M10) | Not to be confused with the symbol M(x, s), the hops M± of the algorithm specification §5.3, or the moments M0 and M2 of step 11b; none of those is a label. |
| L | step 22 (L0–L7) | Not to be confused with a length L or a Lagrangian; neither is a label. The open-item series is L-SP, not the -LS suffix. |
| Z | step 19 (Z1–Z5) | — |
| Thm 1–4, Lem 0–5, Prop 1–3, Cor 4.x | steps 4 and 5 | Numbered without a letter; later notes cite them as 'Theorem 4', 'Lemma 4', 'Proposition 3'. |
| F, H, J, T, W | supplement: four_action_foundations (F1–F4), holland_two_fluid_correspondence (H1–H7), emission_and_absorption (J1–J2), takabayasi_1954_stochastic_picture (T-series), limkumnerd_weighted_paths (W1–W5) | Cited from the analysis notes; see §9. |

Open-item ID series, all of the form `<note prefix><suffix><n>`. The suffixes `SP` and `LS` are not explained in the notes; the two series behave alike.

| Suffix | Where used |
|---|---|
| -SP | G-SP (step 17), S-SP (16), N-SP (18), Y-SP (20), M-SP (21), L-SP (22); supplement: J-SP, H-SP, W-SP, R-SP |
| -LS | K-LS (step 15), Z-LS (step 19) |
| CLS | step 14 |
| CLA | docs/algorithm/compensated_liouville_algorithm.md |
| bare numbers | steps 2 to 13 and 11b: an ordered list, cited as 'open item 3 of the phase-alignment note' |

Vocabulary that misleads (particle, pair, reach, mu, dark, population) is tabulated in [`../supplement/what_the_reach_is.md`](../supplement/what_the_reach_is.md) §7 and is not repeated here.

## 3. Postulates and definitions

Kept apart from the results because a postulate is what everything else assumes, and because two different postulates are both called (S).

### Step 4 — [Phase resonance](phase_resonance_microdynamics.md)

| ID | § | Says | Standing |
|---|---|---|---|
| [Post P0](phase_resonance_microdynamics.md#2-postulates) | 2 | Ontology: two species, positon (ε = +1) and negaton (ε = −1); state (x, p, θ, ε). | there is no P4 |
| [Post P1](phase_resonance_microdynamics.md#2-postulates) | 2 | Winding law: θ̇ = (p²/2m − V)/ℏ along a worldline. | Lem 0 shows it is forced |
| [Post P2](phase_resonance_microdynamics.md#2-postulates) | 2 | Extension: a particle's extended amplitude is ε·exp(i[p(x − x_j)/ℏ + θ]). | — |
| [Post P3](phase_resonance_microdynamics.md#2-postulates) | 2 | Single-valuedness of extended amplitudes on the ring. | Thm 1 derives the pair-level relaxation |
| [Post P5](phase_resonance_microdynamics.md#7-postulate-p5-the-vertex-weight) | 7 | Vertex weight w(δ) = ½(1 + C cos δ), with δ the encountered beat's pattern phase and C its contrast: a local, two-body, Born-type interference rule. | reduced in §13: shape forced by (V1)–(V4), offset δ0 = 0 from Hermiticity of the contact coupling; what remains is generic quantum probability for one binary contact |
| [Def (leg)](phase_resonance_microdynamics.md#3-kinematics-pairs-beats-and-the-parity-theorem) | 3 | A maximal constant-momentum segment of a worldline, carrying (x(t), p, θ(t), ε) and ended at each end by a vertex; a worldline is a chain of legs and a pair is two worldlines, one per species. | — |
| [Def (pair states)](phase_resonance_microdynamics.md#3-kinematics-pairs-beats-and-the-parity-theorem) | 3 | With pair amplitude Ψ = ψa + ψb: dark (equal momenta, gauge-matched phases, Ψ = 0), beating (unequal momenta), gray (equal momenta, gauge-mismatched: uniform residue 2 abs(sin(Δφ/2))). | — |
| [Def (beat, restated)](phase_resonance_microdynamics.md#4-vertices-absorption-emission-resonance) | 4 | A beat is the class (Δp, n̄, χ): it has no worldline; beat number is not conserved, worldline number is. | retired in step 5 (replaced by μ) |
| [Def (absorption, K3)](phase_resonance_microdynamics.md#4-vertices-absorption-emission-resonance) | 4 | Two-body vertex: the excess particle meets the excited partner of a phase-matched beating pair; the pair exits dark. | exit is aligned only for this channel (step 8) |

### Step 5 — [Phase alignment](phase_alignment_microdynamics.md)

| ID | § | Says | Standing |
|---|---|---|---|
| [Def (transported phase)](phase_alignment_microdynamics.md#2-the-misalignment) | 2 | Φ_j(x, t) = θ_j + (p_j (x − x_j) − E_j t)/ℏ on a leg with reference data (x_j, p_j, θ_j): a rewriting of P2, not an addition to it. | — |
| [Def (misalignment)](phase_alignment_microdynamics.md#2-the-misalignment) | 2 | μ = Φ_a − Φ_b (mod 2π) at an event, for a pair of positon a and negaton b. | the relational variable used from step 5 on; the two senses of μ (ladder versus pair) are separated in supplement/what_the_reach_is.md §1.2 |

### Step 6 — [Relational pairing and carrier lock](relational_pairing_and_carrier_lock.md)

| ID | § | Says | Standing |
|---|---|---|---|
| [Post (S)](relational_pairing_and_carrier_lock.md#4-postulate-s-and-the-equivalence-theorem) | 4 | Sea carrier lock: sea particles sharing a cell and a momentum row share a transported phase, up to the pumped misalignment. | unnecessary under permanent pairing (step 7); withdrawn in the phase-alignment algorithm spec; NOT the streaming (S) of step 17 |

### Step 11 — [Open position space](open_position_space.md)

| ID | § | Says | Standing |
|---|---|---|---|
| [Def (H)](open_position_space.md#3-the-coherence-horizon) | 3 | Coherence horizon L_c: a conjugate ket–bra pair with separation above L_c is not instantiated; in position-pair language a maximum rung abs(k) at most L_c/a. A bound on the relative coordinate, not a wall in position. | profile of the cutoff open (CLS4 = CLA1, CLA1b); hard-horizon rate and churn figures: step 14 §4 erratum |
| [Def (R)](open_position_space.md#3-the-coherence-horizon) | 3 | Reach y_max = L_c/2, half the horizon. Earlier drafts read it as the greatest distance from its own position at which a world consults V. | aperture wording superseded by E1: the reach is a period (step 11b) |

### Step 17 — [Compensated ontology](compensated_ontology.md)

| ID | § | Says | Standing |
|---|---|---|---|
| [Post (E)](compensated_ontology.md#1-four-postulates) | 1 | Existence: the world is a locally finite signed counting measure on phase space; no other species or attributes. | — |
| [Post (A)](compensated_ontology.md#1-four-postulates) | 1 | Admissibility: only ensembles whose expectation is the Wigner function of some ρ at least 0 occur. | independent of (S) + (D) (G3.1); particle-level statement open (G-SP2) |
| [Post (S)](compensated_ontology.md#1-four-postulates) | 1 | Streaming: between events every world-particle obeys the full classical force; momentum is continuous along every worldline. | NOT the sea carrier lock of step 6; the notation collision is recorded in step 20 §9 |
| [Post (D)](compensated_ontology.md#1-four-postulates) | 1 | Demography: a parent ionises a neutral sea pair on its own row at rate Γ = sum of abs(K_res); recombination is the reverse; absorptive fraction f = 1/2. | f = 1/2 is the sinkless case (N3) |

### Step 21 — [Compensated ledger in four dimensions](fourd_compensated_ledger.md)

| ID | § | Says | Standing |
|---|---|---|---|
| [Def M0](fourd_compensated_ledger.md#23-definition-m0-the-matching-cell) | 2.3 | The matching cell: r-bin, p_r row, and a centre-of-mass bin of area A in (X, P); A = ∞ is the 1D ledger. | only the p_r row is forced (by the reach) |

### Step 22 — [The sea as a phase reference](sea_phase_reference.md)

| ID | § | Says | Standing |
|---|---|---|---|
| [Def L0](sea_phase_reference.md#5-definition-l0-the-sea-lock) | 5 | The sea lock: each aligned pair's clock equals its row's plane-wave phase at its position, θ = px/ℏ up to a constant per row; a relation between pairs, not within one. | carried by classical flow as the phase field S/ℏ; eroded inside a potential (L7) |

## 4. Results by note

Theorems, propositions, lemmas and corollaries, in ladder order. Within a note the rows follow the document, except where the note's own summary order is kept.

### Step 1 — [Source-document review](phase_space_crystal_lattice_review.md)

No labelled results. See the note's own sections.

### Step 2 — [Four-rule microdynamics](four_rule_microdynamics_equivalence.md)

| ID | § | Says | Standing |
|---|---|---|---|
| [Thm (complete solution family)](four_rule_microdynamics_equivalence.md#4-exactness-theorem) | 4 | Every rate assignment reproducing the QLE stencil exactly is the G-family; the single mediated-jump rule is the G = 0 member. | confirmed on the full slide, Aug 2026 (§0) |
| [Lem (no-go)](four_rule_microdynamics_equivalence.md#62-the-mediation-is-relocated-not-eliminated) | 6.2 | Rates that are multilinear of degree 2 or more in the occupancies (pairwise mass action) cannot equal the linear QLE stencil; a collision microdynamics needs a species whose density is pinned. | loophole (a pinned species) taken up in step 3 |

### Step 3 — [Sea-dressed microdynamics](sea_dressed_microdynamics.md)

No labelled results. See the note's own sections.

### Step 4 — [Phase resonance](phase_resonance_microdynamics.md)

| ID | § | Says | Standing |
|---|---|---|---|
| [Lem 0](phase_resonance_microdynamics.md#2-postulates) | 2 | P1 forces the plane-wave form exp(i(px − Et)/ℏ); any other winding law makes relative phases at meetings frame-dependent. | — |
| [Lem 1](phase_resonance_microdynamics.md#3-kinematics-pairs-beats-and-the-parity-theorem) | 3 | A co-located, equal-momentum, equal-phase pair has Ψ = 0 for all x, t (dark) and stays dark under any local winding law; a separated dark pair de-phases at (V(xa) − V(xb))/ℏ into the gray state. | — |
| [Prop 1](phase_resonance_microdynamics.md#3-kinematics-pairs-beats-and-the-parity-theorem) | 3 | Legs with p+ − p− = 2q·dp form a full-contrast grating at wavelength L/q whose envelope drifts at the mean velocity p̄/m. | recast as a change of variables in step 5 |
| [Thm 1](phase_resonance_microdynamics.md#3-kinematics-pairs-beats-and-the-parity-theorem) | 3 | Single-valuedness puts a lone particle on the even sites of the Wigner half-grid dp = πℏ/L; a beat between even sites lives at their midpoint, whose parity is that of q. | — |
| [Prop 2](phase_resonance_microdynamics.md#4-vertices-absorption-emission-resonance) | 4 | A hop's ΔE/Δp is the midpoint velocity of its two rows and a beat's is its drift velocity, so only co-moving beats can drive a hop (row-resonance selection rule). | strengthened to the exchange theorem (step 5, Thm 4) |
| [Lem 2](phase_resonance_microdynamics.md#5-the-pump-coherence-for-free-populations-untouched) | 5 | Every pumped pair's beat has pattern phase sKx + π/2 − sK v̄ t, independent of where the pair started: coherence of the grating is derived, not postulated. | restated as Lem 5 of step 5 |
| [Lem 3](phase_resonance_microdynamics.md#5-the-pump-coherence-for-free-populations-untouched) | 5 | Sideband occupancies are (Vp τp/2ℏ)²: at first order the pump changes no (x, p) populations, only phases. | true of leg populations, false of pair configurations (step 7 §5.2) |
| [Thm 2](phase_resonance_microdynamics.md#6-theorem-2-no-go-phase-blindness-cannot-be-linear-in-the-potential) | 6 | No-go: no phase-blind rule can be linear in the potential, so phase sensitivity of the event rule is necessary. | does not reach the compensated kernel (Y1) |
| [Thm 3](phase_resonance_microdynamics.md#8-theorem-3-the-rate-law) | 8 | Rate law: the rate field is proportional to sin(Kx), i.e. −V′(x); direction, the per-channel factor and the transfer 2q·dp are all derived. | — |

### Step 5 — [Phase alignment](phase_alignment_microdynamics.md)

| ID | § | Says | Standing |
|---|---|---|---|
| [Lem 4](phase_alignment_microdynamics.md#2-the-misalignment) | 2 | μ is the pair's entire gauge-invariant relational content (unchanged by global gauge and by re-referencing). | used by step 6 |
| [Prop 3](phase_alignment_microdynamics.md#3-the-three-states-and-the-disposal-of-the-beat) | 3 | Winding rates: ∂μ/∂x = Δp/ℏ, and along a path of velocity v, dμ/dt = (Δp/ℏ)(v − v̄_pair). | holds for any index pair (R1.2) |
| [Cor (the beat is dispensable)](phase_alignment_microdynamics.md#3-the-three-states-and-the-disposal-of-the-beat) | 3 | The beat class is Prop 3 plus one constant of integration; Prop 1 of step 4 says μ advances 2π per h/Δp and holds still at v̄. | — |
| [Thm 4](phase_alignment_microdynamics.md#41-statement) | 4.1 | Momentum conservation plus stationarity of μ through a vertex force the vertex to be a momentum swap. | unique only in one spatial dimension (fourd C1) |
| [Cor 4.1](phase_alignment_microdynamics.md#41-statement) | 4.1 | A swap permutes the momentum multiset, so energy conservation is automatic. | selects nothing in d ≥ 2 (fourd C2) |
| [Cor 4.2](phase_alignment_microdynamics.md#41-statement) | 4.1 | After the swap both partners sit at the mate's momentum: the pair exits aligned. | scope: exact (K3) channel only; the write/K4 exit is a winding pair (step 8 §0) |
| [Cor 4.3](phase_alignment_microdynamics.md#41-statement) | 4.1 | The transfer available at a vertex is abs(Δp), exactly 2q·dp for a q-mode pumped pair, with no separate Bragg principle. | — |
| [Lem 5](phase_alignment_microdynamics.md#5-locality-of-the-misalignment) | 5 | After the pump every pumped pair at (x, t) has the same μ = sKx + π/2 − sK v̄ t, so μ is a local datum. | Lemma 2 of step 4 restated |

### Step 6 — [Relational pairing and carrier lock](relational_pairing_and_carrier_lock.md)

| ID | § | Says | Standing |
|---|---|---|---|
| [Prop R1](relational_pairing_and_carrier_lock.md#2-partnership-cannot-be-a-carrier-of-state) | 2 | The two-index misalignments μ_ij are the coboundary of a one-index field Φ, exactly and at every event. | kept, at the level of the one-particle marginal (step 7) |
| [Cor R1.1](relational_pairing_and_carrier_lock.md#2-partnership-cannot-be-a-carrier-of-state) | 2 | A stored partnership index carries no relational information beyond {Φ_j}; it can only select which pairs are admissible. | the reading 'pairing adds no information' retracted in step 7 |
| [Cor R1.2](relational_pairing_and_carrier_lock.md#2-partnership-cannot-be-a-carrier-of-state) | 2 | The winding laws of Prop 3 hold for arbitrary index pairs, not only partners. | — |
| [Prop R2](relational_pairing_and_carrier_lock.md#3-the-obstruction-a-free-carrier-per-pair) | 3 | Residual gauge: a sea of B dark pairs carries U(1)^B (one free phase per pair), not a single global U(1). | — |
| [Cor R2.1](relational_pairing_and_carrier_lock.md#3-the-obstruction-a-free-carrier-per-pair) | 3 | The naive all-pairs proposal fails: with a private random carrier per pair the all-pairs average of cos μ vanishes as 1/B. | cured by (S) here, made unnecessary by permanent pairing (step 7) |
| [Thm R3](relational_pairing_and_carrier_lock.md#4-postulate-s-and-the-equivalence-theorem) | 4 | Under (S) the all-pairs and partnered averages of cos μ coincide, both equal cos μ(x, t); no earlier result stated in μ changes value. | — |
| [Cor R3.1](relational_pairing_and_carrier_lock.md#4-postulate-s-and-the-equivalence-theorem) | 4 | The indexed rule of the specification is the fully ordered limit of the relational rule; (S) is exactly the condition for Lemma 5. | — |
| [Thm R4](relational_pairing_and_carrier_lock.md#5-the-factorisation-theorem) | 5 | Factorisation: the total affine weight over all mediating pairs is w0 N_{r+2sq} N_r + κ Re(Z_{r+2sq} conj Z_r), with Z the per-cell, per-row order parameter. | kept as an implementation device (step 7) |
| [Cor R4.1](relational_pairing_and_carrier_lock.md#5-the-factorisation-theorem) | 5 | The vertex firing probability follows from the factorised weight. | — |
| [Cor R4.2](relational_pairing_and_carrier_lock.md#5-the-factorisation-theorem) | 5 | The encounter loop falls from O(N_exc · B) to O(N_exc + N_sea). | — |
| [Cor R4.3](relational_pairing_and_carrier_lock.md#5-the-factorisation-theorem) | 5 | A sea with random phases has abs(Z) of order N^(−1/2), so an incoherent sea is dark without a separate postulate. | — |
| [Prop R5](relational_pairing_and_carrier_lock.md#8-correction-the-indexed-sea-is-a-consumable-resource) | 8 | Under permanent partnership a pair mediates at most one vertex in its history (a struck pair exits with Δp = 0). | arithmetic confirmed independently in step 7 |
| [Cor R5.1](relational_pairing_and_carrier_lock.md#8-correction-the-indexed-sea-is-a-consumable-resource) | 8 | The pump writes only phases, so a spent pair is never replenished: the split population is non-increasing. | retracted in step 7: sidebands are split-pair amplitudes, so the pump is the source |

### Step 7 — [Permanent pairing, density-matrix reading](permanent_pairing_density_matrix.md)

| ID | § | Says | Standing |
|---|---|---|---|
| [Lem (no-go for lone hoppers)](permanent_pairing_density_matrix.md#4-one-leg-hops-derive-the-stencil-and-the-mediated-counting) | 4 | Any per-particle hop process with position-dependent rates matched to the QLE drift has a symmetric part of coefficient at least abs(Γq)/2: irreducible spurious diffusion. | the load-bearing theorem for the direct proof (split pairs mediate with the pump-excited vertex constant) is stated, not proved (§7) |

### Step 8 — [Coherence ladder](coherence_ladder.md)

| ID | § | Says | Standing |
|---|---|---|---|
| [Lem C1](coherence_ladder.md#2-the-channel-table-lemma-c1) | 2 | At first order in the pump a vertex has exactly three admissible striker momenta: the exact (K3) channel, the leg-local ladder (struck leg's own sideband) and the compound channel (mate's sideband). | label also used in steps 10 and 14 (see registry) |
| [Thm C2](coherence_ladder.md#3-the-ladder-theorem-theorem-c2) | 3 | Ladder theorem: four leg-local channels, sea strikers, phase continuity, conjugate bra factors and one constant reproduce the commutator on every element of every state, to machine precision, freezing exactly at V = 0. | standard of proof reused in step 9; label also used in steps 10 and 14 |

### Step 9 — [Position-pair ladder](position_pair_ladder.md)

| ID | § | Says | Standing |
|---|---|---|---|
| [Thm P1](position_pair_ladder.md#3-theorem-p1-four-one-leg-hops-and-a-diagonal-pump) | 3 | The lattice von Neumann generator is exactly four one-leg hop channels (amplitude ±iJ/ℏ) plus one diagonal pump. | — |
| [Thm P2](position_pair_ladder.md#4-theorem-p2-momentum-is-the-misalignment-of-a-rung-1-pair) | 4 | Lattice probability current j = (2Ja/ℏ) Im ρ_{m+1,m} = (ℏ/ma) abs(ρ1) sin μ: momentum is the misalignment of a rung-1 pair, and the population equation is exactly lattice continuity. | — |
| [Prop P3](position_pair_ladder.md#5-proposition-p3-no-noise-no-force-transposed) | 5 | The pump alone reproduces the Euler force term (second-order convergence in a). | — |
| [Prop P4](position_pair_ladder.md#7-proposition-p4-two-obstructions) | 7 | No fixed background makes ρ + b lie on one ray for all states: no positon-only sea exists in the position representation. | with P5, the two obstructions of §7 |
| [Prop P5](position_pair_ladder.md#8-proposition-p5-resource-arithmetic-and-the-rate-comparison) | 8 | Of the three representations only the position one has an elementary vertex rate with no continuum limit. | — |
| [Thm P6](position_pair_ladder.md#9-theorem-p6-statistical-equivalence) | 9 | The expected empirical ρ of the pair ensemble equals ρ(t) on every element of every state (up to a growth factor carried outside). | — |
| [Cor P7](position_pair_ladder.md) | 10 | The observable sector is a genuine positive-rate particle process guided by sin μ; the coherence sector is not. | — |

### Step 10 — [Four-dimensional phase space](fourd_microdynamics.md)

| ID | § | Says | Standing |
|---|---|---|---|
| [Thm A1](fourd_microdynamics.md#21-the-exactness-theorem-lifts-verbatim) | 2.1 | With shift-by-q for a vector wavevector, the exactness family, participant-locality selection of the symmetric member and endpoint locality all lift unchanged. | — |
| [Prop A2](fourd_microdynamics.md#22-modes-fibre-the-lattice-and-different-modes-fibre-it-differently) | 2.2 | Mode q decomposes the joint momentum lattice into disjoint chains (lines along q); the four-action process runs independently on each. | — |
| [Thm A3](fourd_microdynamics.md#23-what-leaks-the-mode-wavevector-decides) | 2.3 | Leak law: a momentum direction u is conserved event by event iff u is orthogonal to every active mode wavevector; focus and defocus never change u·P. | corrects the 4-D supplement's 2p/1D versus 1p/2D table (§0.1); modular companion is O5; collapses to its neutral half under (S) (step 21 M2) |
| [Prop B1](fourd_microdynamics.md#31-the-joint-sea-is-a-product-sea-but-the-shift-is-not-a-product) | 3.1 | For a product state the joint sea is not the product of single-particle seas: the missing cross terms are −(2/h)(W1 + W2). | — |
| [Prop B2](fourd_microdynamics.md#32-the-sixteen-channels-lift-and-the-sea-stops-being-optional) | 3.2 | Excess worlds per cell fall as W times (M_x M_p) to the power −dN: excess–excess collisions essentially never occur once N > 1, so the sea is the only collision partner. | — |
| [Prop B3](fourd_microdynamics.md#33-the-cost-stated-plainly) | 3.3 | For a Gaussian the peak of abs(W) in units of (2/h)^d equals the state's purity: entanglement is excess-to-background loss. | — |
| [Thm C1](fourd_microdynamics.md#42-theorem-4-loses-its-uniqueness) | 4.2 | In d ≥ 2 the two conditions of Thm 4 leave a (d − 1)-parameter family of vertices. | corrects the uniqueness of step 5 Thm 4 (d = 1 only); label collides with steps 8 and 14 |
| [Cor C2](fourd_microdynamics.md#42-theorem-4-loses-its-uniqueness) | 4.2 | Every member of that family conserves energy identically, so energy does not select the swap. | label collides with steps 8 and 14 |
| [Prop C3](fourd_microdynamics.md#42-theorem-4-loses-its-uniqueness) | 4.2 | The swap is the unique member whose out-momenta are a permutation of the in-momenta. | — |
| [Prop C4](fourd_microdynamics.md#43-what-the-transverse-parameter-is) | 4.3 | Every member gives the excess particle the same transfer; at a transversally uniform sea the family is invisible to the QLE. | invisibility fails in a live 4-D sea (open item 1); does not transfer to the compensated ledger (step 21 M8) |
| [Thm D](fourd_microdynamics.md#55-the-cost-an-exact-mean-bought-with-unbounded-noise) | 5.5 | Mode-independent noise: for the ring harmonic, rate amplitude times squared hop displacement is 2mω²ℏ for every q, so injected momentum variance grows linearly with the mode cutoff. | label collides with step 13; does not arise under (S) (step 21 M3) |

### Step 11 — [Open position space](open_position_space.md)

| ID | § | Says | Standing |
|---|---|---|---|
| [Thm O1](open_position_space.md#2-the-wigner-kernel-has-no-position-envelope) | 2 | The Wigner kernel has modulus independent of position for every potential, so a distant world is struck at full rate and the escape problem is not a boundary problem: no absorber can fix it. | finite-reach refinement is C7 |
| [Cor O1.1](open_position_space.md#2-the-wigner-kernel-has-no-position-envelope) | 2 | The total vertex rate per world saturates to a constant as abs(x) grows. | — |
| [Thm O2](open_position_space.md#2-the-wigner-kernel-has-no-position-envelope) | 2 | The kernel's first moment is exactly the classical force while the zeroth and second moments of abs(V_W) are flat, so signal-to-noise of any signed-particle force estimator decays like V′. | — |
| [Prop O3](open_position_space.md#3-the-coherence-horizon) | 3 | Truncation is not absorption: a windowed kernel stays odd in ξ, so signed world number is conserved exactly for every L_c and x. | energy conservation closed by E3 (step 11b) |
| [Prop O4](open_position_space.md#3-the-coherence-horizon) | 3 | Under a coherence horizon all vertex activity is confined to dist(x, supp V) at most L_c/2. | extended by C7; the edge becomes a profile edge if CLA1b is adopted |
| [Thm O5](open_position_space.md#42-the-coset-invariant) | 4.2 | Coset invariant: each world's residue p mod πℏ/a is constant for all time on all of the real line, with no box and no horizon. | modular companion of A3; survival under the compensated split is CLS3/CLA5 (C9) |
| [Cor O5.1](open_position_space.md#43-the-ring-is-a-sector-not-an-approximation) | 4.3 | An L = na ring is exactly n independent copies of the a-ring problem at different offsets, not an approximation to open space. | — |

### Step 11b — [What the reach controls](reach_energy_coupling.md)

| ID | § | Says | Standing |
|---|---|---|---|
| [Thm E1](reach_energy_coupling.md#1-the-reach-is-a-period-not-an-aperture) | 1 | Momentum transfers form a lattice iff D_x(y) = V(x+y) − V(x−y) is periodic in y, with Δp = πℏ/period; ring, periodic potential and horizon are three sources of one period. | restates Def (R) as a period, not an aperture (§0) |
| [Prop E1.1](reach_energy_coupling.md#1-the-reach-is-a-period-not-an-aperture) | 1 | A sharp window is exact on the lattice iff 2y_max is a whole number of the potential's periods: commensuration, not sharpness. | corrects step 11 §3.2 |
| [Thm E2](reach_energy_coupling.md#2-folding-and-the-trilemma) | 2 | Folding identity: folding a decaying D_x into one period gives exactly the lattice kernel of the periodised potential. | folding under the compensated split is open item 5 |
| [Thm E3](reach_energy_coupling.md#3-conservation-does-not-depend-on-the-reach) | 3 | For any even window (sharp, tapered or folded) M0 = M2 = 0: signed world number and energy are reach-independent. | closes open item 4 of step 11 |
| [Thm E4](reach_energy_coupling.md#4-the-four-action-ledger-with-more-than-one-mode) | 4 | For the symmetric four-rule member, per mode: the focus channel does no net work and the hop channel delivers the whole classical power. | closed lattice only (open item 1) |
| [Thm E5](reach_energy_coupling.md) | 5 | An even horizon profile applied to the full kernel contaminates the third momentum moment of ℏ²V‴/4 by a term of order 1/y_max². | — |
| [Thm E6](reach_energy_coupling.md#6-compensation-removes-it-exactly) | 6 | Applied to the compensated residual instead, the third moment is exactly ℏ²V‴/4 with no y_max dependence: the reach is inert at order ℏ². | second justification for the compensated split |
| [Thm E7](reach_energy_coupling.md#7-a-polynomial-potential-has-no-jump-measure-at-all) | 7 | For a polynomial V the Moyal series terminates, so there is no jump measure on the open line until a reach is imposed (quartic double well worked). | used by step 15 §0.2 |
| [Thm E8](reach_energy_coupling.md#8-what-the-reach-buys-and-what-it-costs) | 8 | Under the compensated split the residual event budget grows without bound with the reach for every V with V′(x) not zero (linearly for bounded V, cubically for the quartic). | — |

### Step 12 — [Interworld coupling](interworld_coupling.md)

| ID | § | Says | Standing |
|---|---|---|---|
| [Prop I1](interworld_coupling.md#1-tutorial-the-only-way-the-potential-can-enter) | 1 | The coupling U = V(x1) − V(x2) vanishes at coincidence, is antisymmetric under leg exchange and vanishes for a free particle. | the free-particle property is the sharpest test of any interworld force law |
| [Prop I2](interworld_coupling.md#2-midpoint-and-separation) | 2 | For one cosine mode the coupling factorises into a midpoint amplitude (the classical force) and a separation grating of twice the potential's period. | derives the half-quantum offset that step 3 postulates |
| [Thm I3](interworld_coupling.md#3-theorem-i3-channels-are-the-separation-spectrum) | 3 | The available momentum channels are exactly the Fourier spectrum of the coupling in the separation Y: discrete iff the coupling is periodic. | rigorous; reading Y as a physical separation is an interpretive postulate (§0) |
| [Cor I3.1](interworld_coupling.md#3-theorem-i3-channels-are-the-separation-spectrum) | 3 | M active modes give 2M shifts and 4M rules: a one-mode potential has exactly four; a non-periodic one has a continuum. | — |
| [Thm I4](interworld_coupling.md#4-theorem-i4-the-moyal-series-is-the-separation-expansion) | 4 | The Y^(2l+1) term of the coupling is the l-th Moyal term: a coupling linear in Y is exactly classical, so harmonic and inverted-harmonic potentials have no jump channel. | used by G2 and step 15 |
| [Prop I5a](interworld_coupling.md#5-what-a-genuine-pair-interaction-would-look-like-for-contrast) | 5 | A pair potential gives an anti-correlated, total-momentum-conserving pair shift: a focus/defocus channel and nothing else. | — |
| [Thm I5](interworld_coupling.md#6-theorem-i5-the-coupling-winds-it-does-not-push) | 6 | Under the potential term alone abs(ρ) is constant and dμ/dt = −U/ℏ: the coupling winds the misalignment, exerts no force and does no work. | — |

### Step 13 — [Species, sectors, annihilation](species_sectors_and_annihilation.md)

| ID | § | Says | Standing |
|---|---|---|---|
| [Thm D1](species_sectors_and_annihilation.md#1-theorem-d1-species-and-phase-are-conjugate-not-paired) | 1 | Species of a Wigner world-particle and the phase of a relational pair are one degree of freedom in conjugate bases; no carrier holds both sharply. | — |
| [Thm D0](species_sectors_and_annihilation.md#2-theorem-d0-the-two-ensembles-have-no-carrier-correspondence) | 2 | No map sends E1 (Wigner) carriers to E2 (density-matrix) carriers one at a time: the species censuses are anti-correlated and 'positon' in the two layers is a homonym. | structural question open (item 1) |
| [Cor D0.1](species_sectors_and_annihilation.md#2-theorem-d0-the-two-ensembles-have-no-carrier-correspondence) | 2 | Results proved in one ensemble do not transfer to the other without a separate derivation. | — |
| [Thm D2](species_sectors_and_annihilation.md#3-theorem-d2-what-is-dark-and-which-sea-is-which) | 3 | The crystal shift W → W + 2/h is ρ → ρ + 2·1: inert because the identity commutes with every Hamiltonian; in the pair basis it sits at zero leg separation. | replaces the supplement's §6 admissibility argument |
| [Thm D2.1](species_sectors_and_annihilation.md#3-theorem-d2-what-is-dark-and-which-sea-is-which) | 3 | An operator is dark under every Hamiltonian iff it is c·1. | — |
| [Thm D2.2](species_sectors_and_annihilation.md#3-theorem-d2-what-is-dark-and-which-sea-is-which) | 3 | The neutral sea (c = 0) and the crystal shift (c = 2) have near-opposite roles; 2·1 cannot be the pumped medium. | — |
| [Cor D2.3](species_sectors_and_annihilation.md#3-theorem-d2-what-is-dark-and-which-sea-is-which) | 3 | The number 2/h plays three roles: the bound abs(W) at most 2/h, the shift value and the reservoir capacity per cell (borrowing the first, not the second). | — |
| [Cor D2.4](species_sectors_and_annihilation.md#3-theorem-d2-what-is-dark-and-which-sea-is-which) | 3 | Non-compactness: the identity is not trace class, which is why the sea-dressed layer does not generalise naively to open space. | — |
| [Thm D3](species_sectors_and_annihilation.md) | 4 | I2 and I5 extend to arbitrarily many cosine modes. | — |
| [Prop D4](species_sectors_and_annihilation.md) | 4 | A background of finite coherence length ε is not dark: its winding rate is linear in ε; only exactly zero separation is dark. | an earlier Cor D4.1 was withdrawn (§0) |
| [Thm D5](species_sectors_and_annihilation.md#5-theorems-d5-and-d6-the-reach-of-annihilation) | 5 | A Wigner-basis annihilation event is local in (X, p) but its Weyl image has support of order the coherence length: it cannot be local in both descriptions. | — |
| [Cor D5.1](species_sectors_and_annihilation.md#5-theorems-d5-and-d6-the-reach-of-annihilation) | 5 | Anonymity is exact: two opposite-species E1 carriers in a cell are interchangeable for all future time; annihilation never needs the original partner. | voids criterion (ii) of supplement §7.4 |
| [Cor D5.2](species_sectors_and_annihilation.md#5-theorems-d5-and-d6-the-reach-of-annihilation) | 5 | The crystal is still needed: on a continuous position axis coincidence has measure zero. | — |
| [Thm D6](species_sectors_and_annihilation.md#5-theorems-d5-and-d6-the-reach-of-annihilation) | 5 | Convolving W with a Gaussian of width σ_p multiplies the pair kernel by exp(−σ_p² Y²/2ℏ²): an imposed leg coherence length ℏ/σ_p (the soft blob). | — |
| [Thm D8](species_sectors_and_annihilation.md) | 6 | The net excess positon count in a position column is exactly the diagonal of ρ, hence non-negative. | — |
| [Thm D9](species_sectors_and_annihilation.md) | 6 | Diagonal ρ gives W at least 0; the off-diagonal part has vanishing column sums: the excess splits into a Born-density sector and a column-balanced coherence sector. | — |
| [Thm D10](species_sectors_and_annihilation.md) | 6 | Complementary conservation: the jump substep conserves every column sum and streaming every row sum, giving 2(M + N) free invariants. | — |
| [Prop U1](species_sectors_and_annihilation.md#7-proposition-u1-why-a-quadratic-channel-is-unavoidable) | 7 | No unravelling linear in the ensemble is L1-stationary: the expected L1 ledger grows at exponent ρ(abs(L)) = 2.341. | does not obstruct the stabiliser (S1.1) |
| [Thm D12](species_sectors_and_annihilation.md) | 8 | Species conjugation is momentum reflection: the species label is the orientation of the momentum transfer, not an independent charge. | — |
| [Thm D13](species_sectors_and_annihilation.md) | 8 | A signed carrier carries a phase offset restricted to 0 or π by Hermiticity: structure group Z2 for E1 versus U(1) for E2. | — |
| [Thm D14](species_sectors_and_annihilation.md) | 8 | Species is relational, not the carrier's own handedness sign(p). | — |
| [Thm D15](species_sectors_and_annihilation.md) | 9 | The four actions split and combine bound sea pairs rather than creating them, so positon and negaton numbers are each exactly conserved. | — |
| [Thm D16](species_sectors_and_annihilation.md) | 9 | Sizing floor: the sea is a finite local reservoir; in units of the Wigner capacity the blocked-split fraction is universal. | — |
| [Prop D7](species_sectors_and_annihilation.md#104-what-it-preserves) | 10.4 | Same-cell annihilation is exactly unbiased: it changes E, W and every moment by zero. | — |
| [Thm D11](species_sectors_and_annihilation.md#104-what-it-preserves) | 10.4 | Every member of the annihilation family from cell-exact to whole-column leaves the diagonal of ρ invariant; they differ only in how much off-diagonal they destroy. | — |
| [Thm D17](species_sectors_and_annihilation.md#11-theorem-d17-adaptive-sea-allocation) | 11 | Injecting or removing a bound pair changes E by zero, so sea size is a representational choice, free per cell and step (0.7–1.8 per cent of the world-particle count suffices). | — |
| [Prop D18](species_sectors_and_annihilation.md#11-theorem-d17-adaptive-sea-allocation) | 11 | Part J: the adaptive scheme reproduces the uniform sea's error (0.0093) at 1.3–1.8 per cent of its peak particle count (η = 30 to 100). | — |

### Step 14 — [Compensated Liouville splitting](compensated_liouville_splitting.md)

| ID | § | Says | Standing |
|---|---|---|---|
| [Thm C1](compensated_liouville_splitting.md#2-the-split) | 2 | M_cl and M_res are multiplication operators in the same variables, so they commute exactly: the classical force splits off with no Trotter error. | label collides with steps 8 and 10 |
| [Thm C2](compensated_liouville_splitting.md) | 2.2 | M_res is exactly the odd part of the cubic Taylor remainder of V about x. | label collides with steps 8 and 10 |
| [Thm C3](compensated_liouville_splitting.md#3-the-reach-theorem) | 3 | Reach theorem: at bounded reach the residual kernel has zero zeroth and first moments, a bounded signed jump measure that conserves worlds and carries no net momentum. | — |
| [Lem C8](compensated_liouville_splitting.md#31-the-residual-has-no-hop-rate-on-the-momentum-lattice) | 3.1 | For every member of the FR family the focus channel carries no first moment: the classical force sits necessarily in the one-photon hop channel. | half-settles CLS2 |
| [Thm C9](compensated_liouville_splitting.md#31-the-residual-has-no-hop-rate-on-the-momentum-lattice) | 3.1 | The compensating term Γ u/sin u is non-periodic and singular, so no compensated member is a finite-range lattice rate law; the lattice survives only as a transfer spectrum. | resolves CLS3 (see CLA5) |
| [Thm C4](compensated_liouville_splitting.md#4-the-reach-condition) | 4 | Per mode the residual-to-classical ratio is abs(sin u − u)/u, about u²/6 for u = k·y_max; the split gains for u well below π and loses beyond π/2; reach and momentum quantum are one parameter. | TV table grid-dependent: erratum in §4, see algorithm spec §4.4 |
| [Thm C7](compensated_liouville_splitting.md#5-the-quiet-region) | 5 | Quiet region: if V‴ vanishes on [x − y_max, x + y_max] a world at x takes no events and moves on an exact Newtonian trajectory. | §5.1's 'translated by the reach' corrected by K2; finite-reach refinement of O1, extends O4 |
| [Thm C5](compensated_liouville_splitting.md#6-the-ring-as-a-diagnostic) | 6 | On a circle the residual vanishes iff V is constant: a ring pins every world at u = qπ, so it is not a valid testbed for the reach condition. | — |
| [Prop C6.1](compensated_liouville_splitting.md#61-the-bowtie-c7-where-the-potential-is-nowhere-quadratic) | 6.1 | On the ring M_res is exactly zero iff abs(x) + abs(y) is below L/2 (the bowtie). | label belongs to §6.1, not to Thm C6 |
| [Cor C6.2](compensated_liouville_splitting.md#61-the-bowtie-c7-where-the-potential-is-nowhere-quadratic) | 6.1 | A coherence horizon puts the ring seam out of reach: the periodised parabola is exactly Newtonian for abs(x) below L/2 − L_c/2. | label belongs to §6.1, not to Thm C6 |
| [Lem C0](compensated_liouville_splitting.md#62-modes-cancel-and-badly) | 6.2 | For one cosine mode the finite-difference stencil is an identity, not an approximation, so the mode sum reproduces the whole Moyal series. | supplement §6.1 wording needs an erratum (CLS5) |
| [Thm C6](compensated_liouville_splitting.md#7-coulomb) | 7 | For V = −Z/x the Moyal series is geometric and converges iff the reach misses the nucleus. | generalised by K1 (nearest complex singularity) |

### Step 15 — [Eckart barrier](eckart_barrier_compensated.md)

| ID | § | Says | Standing |
|---|---|---|---|
| [Thm K1](eckart_barrier_compensated.md#2-theorem-k1-the-reach-ceiling-is-the-analyticity-strip) | 2 | The Moyal series in the half separation converges at x iff the reach is below the distance to the nearest complex singularity of V (the reach ceiling). | generalises C6 |
| [Cor K1.1](eckart_barrier_compensated.md#2-theorem-k1-the-reach-ceiling-is-the-analyticity-strip) | 2 | For sech² the ceiling is y_max below πa/2, which with E1 gives Δp above ℏ/a: a reach-limited lattice cannot resolve the barrier's own momentum scale. | sharpened by Z4 (K-LS2, CLA7) |
| [Cor K1.2](eckart_barrier_compensated.md#2-theorem-k1-the-reach-ceiling-is-the-analyticity-strip) | 2 | Soft-core Coulomb has ceiling R(x) = sqrt(x² + ε²); at the origin the ceiling is the softening length. | realised in step 19 (Z3) |
| [Thm K2](eckart_barrier_compensated.md#3-theorem-k2-the-far-field-and-an-erratum) | 3 | Far field: the residual decays at exactly V's own rate, so for an exponential tail the reach rescales the interaction profile rather than translating it. | corrects splitting §5.1 |
| [Thm K3](eckart_barrier_compensated.md#4-theorem-k3-on-the-open-line-the-split-never-loses) | 4 | For the Eckart barrier the spectrum-weighted budget ratio saturates at 1 from below: on the open line compensation never loses (opposite of the ring). | — |
| [Thm K4](eckart_barrier_compensated.md#5-theorem-k4-the-deterministic-step-carries-no-tunnelling) | 5 | The classical outcome functional is exactly invariant under streaming plus deterministic acceleration, so the whole quantum–classical transmission gap comes from the residual channel. | — |
| [Thm K5](eckart_barrier_compensated.md#6-theorem-k5-two-large-flows-one-small-difference) | 6 | The quantum correction to transmission is the time-integrated imbalance of two flows of pairs across the separatrix, each several times larger than their difference (net/gross about 0.19). | — |
| [Thm K6](eckart_barrier_compensated.md) | 7 | For a packet centred on the barrier top with σ_r = a, net/gross tends to c/β with c about 0.55: fixed relative accuracy in T costs β² particles. | off-centre dependence open (K-LS3) |
| [Prop K8](eckart_barrier_compensated.md#81-a-jump-is-not-a-hop) | 8.1 | If an event conserves body momentum and the parent streams on, the pair it produces is ionised from a bound sea pair at the parent's own row. | its content rests on the sea's ineligibility for recombination (Y2) |
| [Thm K7](eckart_barrier_compensated.md#82-theorem-k7-the-emission-bias-reverses-four-times) | 8.2 | For V0 sech²(r/a) the emission rate vanishes and the sign of K_q reverses at r = 0 and r = ±a·artanh(sqrt(2/3)): four lobes and a quiet summit. | nodes are minima, not exact zeros of Γ, at finite reach (Z-LS5) |
| [Thm K9](eckart_barrier_compensated.md#84-what-the-classical-trajectory-reading-costs) | 8.4 | Supply condition for the ledger to close with every trajectory Newtonian; a standing dressed population of 3–10 bodies at reach 4πa is the usable window. | population fixed by N5 (λ* = 1.228 per species per cell) |

### Step 16 — [Sea population equilibrium](sea_population_equilibrium.md)

| ID | § | Says | Standing |
|---|---|---|---|
| [Prop S0](sea_population_equilibrium.md#3-what-an-event-is) | 3 | If an event conserves body momentum and the parent streams on, the bound pair it consumes sits at the parent's own momentum row. | — |
| [Thm S1](sea_population_equilibrium.md#1-the-ledger-and-why-it-is-not-determined-by-the-dynamics) | 1 | A removal channel preserves E for arbitrary ensembles iff it removes positons and negatons in coincident pairs. | — |
| [Cor S1.1](sea_population_equilibrium.md#1-the-ledger-and-why-it-is-not-determined-by-the-dynamics) | 1 | Prop U1 does not obstruct the stabiliser; it names the class the stabiliser must lie outside. | — |
| [Thm S2](sea_population_equilibrium.md#4-the-two-field-split-and-the-fixed-point) | 4 | The E equation is closed and exactly the QLE whatever N does; N obeys the absolute-value kernel with a closed-form local fixed point. | — |
| [Thm S3](sea_population_equilibrium.md#4-the-two-field-split-and-the-fixed-point) | 4 | Two relaxation laws: exponential at rate about κ·abs(E) inside the emitting region, and N tending to 2/(κt) where Γ = E = 0. | — |
| [Thm S4](sea_population_equilibrium.md#5-the-emissive-unravelling-fails) | 5 | Under the emissive unravelling the sea is relocated, not consumed: the worst-cell deficit drains monotonically without bound. | verdict unaffected by the demo repair (§0 erratum) |
| [Thm S5](sea_population_equilibrium.md#5-the-emissive-unravelling-fails) | 5 | Throttling the rate by actual occupancy makes the E equation no longer the QLE: 40 per cent error in the core, norm and ⟨p⟩ exact. | — |
| [Thm S6](sea_population_equilibrium.md#6-the-other-realisation) | 6 | An event is wholly emissive or wholly absorptive (a mixed one would change body momentum by 2ξ); absorption needs a partner at both daughters and is supply-limited. | — |
| [Thm S7](sea_population_equilibrium.md#6-the-other-realisation) | 6 | Ledger identity: with absorptive fraction f, dN = 2(1 − 2f) n_ev and dS = (2f − 1) n_ev, so f = 1/2 closes both ledgers. | corrected: f = 1/2 is the sinkless case; the law is N3 (erratum §0) |
| [Thm S8](sea_population_equilibrium.md#6-the-other-realisation) | 6 | Absorptive unravelling on the Eckart summit restores QLE fidelity by four orders at fixed dt (measured f = 0.434). | figures improved by the transport/kernel repair (erratum §0) |
| [Thm S9](sea_population_equilibrium.md#6-the-other-realisation) | 6 | f is an attractor near 1/2, approached from both sides and independent of the initial ensemble; so closure carries no information about the recombination constant. | attractor sits slightly below 1/2 (S-SP3); reached only once the population fills its reachable cells (step 21 M5, M6) |

### Step 17 — [Compensated ontology](compensated_ontology.md)

| ID | § | Says | Standing |
|---|---|---|---|
| [Thm G1](compensated_ontology.md#2-g1-the-postulates-give-the-qle-and-this-is-not-new-work) | 2 | Assembly: under (E), (S), (D) the expectation of the world measure solves the QLE exactly. | an assembly of C1–C3, K4, K8 and S7 |
| [Thm G2](compensated_ontology.md#3-g2-the-demographic-channel-is-empty-for-quadratic-hamiltonians) | 3 | For quadratic V the demographic channel is empty at every reach (1.2 × 10⁻¹⁵): creation and annihilation cannot be the whole of the quantum. | same obstruction as I4 |
| [Thm G3](compensated_ontology.md) | 4 | The residual generator and the admissibility constraint are independent functions of ℏ: the first can go to zero continuously while the second is held fixed. | — |
| [Prop G3.1](compensated_ontology.md#5-postulate-a-is-independent-and-the-violated-inequality-is-the-wigner-bound) | 5 | (A) is not derivable from (S) + (D): four Gaussians the dynamics cannot tell apart differ by a factor eight in phase-space area; the separating inequality is abs(W) at most 2/h. | — |
| [Thm G4](compensated_ontology.md#6-g4-streaming-alone-leaves-the-admissible-set) | 6 | (A) is propagated by (S) + (D) but not by (S) alone: the least eigenvalue of ρ goes from 1e-8 to −0.10 under classical transport. | quantification open (G-SP3) |
| [Thm G5](compensated_ontology.md#7-g5-the-census-is-a-regulator-and-the-generator-is-not) | 7 | Γ grows without bound with the coherence reach while the generator converges to 1e-14: how many worlds exist is a property of the regulator. | the outstanding defect: G-SP1 (= CLA10) |

### Step 18 — [Stochastic ledger](stochastic_ledger.md)

| ID | § | Says | Standing |
|---|---|---|---|
| [Thm N1](stochastic_ledger.md#2-theorem-n1-the-pair-count-is-a-pathwise-invariant) | 2 | P = S + N/2 is conserved pathwise, in exact integers, on every trajectory: a deterministic invariant of the generator, not a martingale. | strengthens J1 |
| [Thm N2](stochastic_ledger.md#3-theorem-n2-the-realisation-is-a-null-direction) | 3 | Both realisations of an event move E identically, so the Bernoulli(f) choice lies in the kernel of the observable map; the only noise W sees is Poisson event timing. | not a stochastic mechanics in Nelson's sense (§6.2) |
| [Thm N3](stochastic_ledger.md#4-theorem-n3-the-sum-rule-and-what-it-corrects) | 4 | Sum rule: in any stationary state Γ_tot (1 − 2f) = R_sink, the rate of every other body-removing channel; f = 1/2 is the sinkless case. | corrects S7 and the sign-change argument of S-SP3 |
| [Thm N4](stochastic_ledger.md#51-the-closure) | 5.1 | The availability closure contains no channel index, so the standing population is independent of the number of channels. | restricted: under streaming availability is structural, not Poisson (step 21 M7) |
| [Thm N5](stochastic_ledger.md#52-five-numbers) | 5.2 | Under independent occupancy f = (1 − e^(−λ))², so f = 1/2 fixes λ* = 1.227947 bodies per species per cell, f′(λ*) = √2 − 1, Var/M = 2.414214 and Fano factor 0.983028, with no free constant. | N-SP2, N-SP3 open; restricted with N4 (step 21 M7) |
| [Thm N6](stochastic_ledger.md#7-theorem-n6-transport-is-the-local-regulator) | 7 | Transport, not recombination, is the local regulator: recombination damps per-cell spread only partially and by N3 costs f; streaming holds the spread flat at no cost. | sea side untested (R-SP7) |

### Step 19 — [Soft-core Coulomb](soft_core_coulomb.md)

| ID | § | Says | Standing |
|---|---|---|---|
| [Thm Z1](soft_core_coulomb.md#1-theorem-z1-the-quiet-points) | 1 | V‴ of −Z/sqrt(r² + ε²) is (Z/ε⁴) u (6u² − 9)(1 + u²)^(−7/2): zeros at r = 0 and ±ε·sqrt(3/2), so the soft core carries four emission lobes with a quiet radius set by ε alone. | — |
| [Cor Z1.1](soft_core_coulomb.md#1-theorem-z1-the-quiet-points) | 1 | Pure Coulomb has V‴ = −6Z/r⁴ and two lobes; softening manufactures the inner pair, which collapses onto the origin as ε → 0. | — |
| [Thm Z2](soft_core_coulomb.md#2-theorem-z2-the-nucleus-is-dark-and-the-ring-survives) | 2 | Γ(0) = 0 at every reach: the nucleus is dark exactly as the sech² summit is; the interior quiet ring drifts outward by 0.4 per cent at y_max = ε/4 and 6.1 per cent at 0.99ε. | quiet points are minima, not exact nodes (Z-LS5) |
| [Thm Z3](soft_core_coulomb.md#3-theorem-z3-a-ceiling-that-moves) | 3 | The ceiling R(x) = sqrt(x² + ε²) varies with position, so the momentum quantum cannot be both position-independent and everywhere maximal: a non-uniform lattice, or a uniform y_max below ε. | non-uniform lattice deferred (Z-LS1) |
| [Thm Z4](soft_core_coulomb.md#4-theorem-z4-the-softening-threshold) | 4 | For a uniform reach below ε, resolving the ground state to k rungs per σ_p needs ε at least k⁴π⁴/4 in units of a0 (24.35 for one rung, 389.6 for two); an unsoftened atom cannot live on a uniform reach-limited crystal. | soft-core form of K-LS2 |
| [Thm Z5](soft_core_coulomb.md#5-theorem-z5-but-the-channel-is-not-empty-there) | 5 | At the tightest uniform reach the horizon spans 2.8 σ_r and Γ(σ_r)/Γ_max = 0.974: unlike the Eckart barrier, the window where the lattice works is not the window where the demographic channel is empty. | makes this the vehicle for S-SP6 (Z-LS4) |

### Step 20 — [Dark sea and worldline identity](dark_sea_and_worldline_identity.md)

| ID | § | Says | Standing |
|---|---|---|---|
| [Thm Y1](dark_sea_and_worldline_identity.md#2-theorem-y1-the-no-go-does-not-reach-the-compensated-kernel) | 2 | The compensated kernel is linear in the potential, K(λV) = λK(V), so the phase-blind no-go (step 4 Thm 2) does not reach it. | answers J-SP4 in part; wording of §2 corrected in step 22 |
| [Cor Y1.1](dark_sea_and_worldline_identity.md#2-theorem-y1-the-no-go-does-not-reach-the-compensated-kernel) | 2 | A particle-level phase cannot be necessary for E; any phase rule lives in the kernel of the observable map. | — |
| [Thm Y2](dark_sea_and_worldline_identity.md#3-theorem-y2-what-the-seas-ineligibility-buys) | 3 | All realisations move E identically and conserve P; with the sea ineligible every event has abs(ΔS) = 1, with it eligible some have ΔN = ΔS = 0 and reproduce the excluded hop. | ineligibility is what gives (S) and K8 their content |
| [Thm Y3](dark_sea_and_worldline_identity.md#4-theorem-y3-the-consumed-bodies-are-a-pair) | 4 | The two bodies a catalysed recombination consumes are a winding pair with splitting 2ξ and midpoint at the parent's row: abs(Ψ) = 2 abs(sin(μ/2)), envelope drifting at the parent's velocity. | — |
| [Thm Y4](dark_sea_and_worldline_identity.md#5-theorem-y4-the-kink-and-its-phase-ramp) | 5 | A momentum kink of ±ξ carries a phase ramp pinned where it happens, so the created pair's μ freezes at μ(x_k): dark exactly when the kink sits at a node. | — |
| [Cor Y4.1](dark_sea_and_worldline_identity.md#5-theorem-y4-the-kink-and-its-phase-ramp) | 5 | A phase reset at the parent's position and a node-located kink are one freedom seen twice. | — |
| [Thm Y5](dark_sea_and_worldline_identity.md#6-theorem-y5-two-clocks-force-piecewise-worldlines) | 6 | If world-particles are conserved and no body changes momentum at an event, per-row species counts are event-invariant and the residual channel can do nothing: conservation forces kinks. | the alternative is birth and death |
| [Thm Y6](dark_sea_and_worldline_identity.md#7-theorem-y6-who-crosses-the-barrier) | 7 | Under piecewise worldlines a tagged body's momentum walk has zero drift and diffusion D_p = (1/2) sum of ξ_q² abs(K_q); tagged crossing (0.327) tracks neither classical (0.069) nor quantum (0.213) transmission. | corrects Eckart §8.3 and the ledger note §1 |

### Step 21 — [Compensated ledger in four dimensions](fourd_compensated_ledger.md)

| ID | § | Says | Standing |
|---|---|---|---|
| [Prop M1](fourd_compensated_ledger.md#21-proposition-m1-the-pair-kernel-moves-relative-momentum-only) | 2.1 | For a pair potential the residual symbol depends on s only through s1 − s2: every event transfers (ξ, −ξ), moves p_r, and leaves X, P and r unchanged. | compensated counterpart of step 10's A2 and §4.4 |
| [Thm M2](fourd_compensated_ledger.md#22-theorem-m2-every-event-is-momentum-neutral) | 2.2 | Under (S) every event leaves every momentum functional unchanged; all exchange with the field goes through the force. | collapses step 10's A3 to its neutral half |
| [Cor M3](fourd_compensated_ledger.md#31-corollary-m3-the-quadratic-sectors-are-free) | 3.1 | Quadratic sectors carry no events, so step 10's Theorem D and its §5.4 and §5.6 obstructions do not arise. | from G2 and I4 |
| [Prop M4](fourd_compensated_ledger.md#32-proposition-m4-the-event-budget-depends-on-the-reach) | 3.2 | Measured on a common lattice at y_max = 2π: residual over uncompensated maximal rate is 0 (harmonic), 0.834 (cosine), 2.825 (Pöschl–Teller). | as step 14's C4 predicts; M-SP4 |
| [Prop M5](fourd_compensated_ledger.md#4-proposition-m5-the-population-fills-cells) | 4 | Measured: body count proportional to occupied joint cells (2.4 to 4.6 per cell); A = h and h/4 multiply the saturated count by 4.8 and 13.2; error grows about as √N. | the note's centre |
| [Prop M6](fourd_compensated_ledger.md#51-proposition-m6-the-cells-are-the-windows) | 5.1 | Measured: the cells filled are those escaped bodies reach, not the state's support; doubling the window raises the count. | M-SP1, M-SP2 |
| [Prop M7](fourd_compensated_ledger.md#52-proposition-m7-availability-is-structural) | 5.2 | Measured: under streaming, per-leg availability stays 0.55 to 0.74 for mean requested occupancy 1.3 to 10.3; not Poisson. | restricts N4, N5 |
| [Prop M8](fourd_compensated_ledger.md#6-proposition-m8-co-location-in-the-centre-of-mass-is-not-a-gauge-choice) | 6 | Measured: at A = ∞, partners from the wrong centre-of-mass component shift the conditional ⟨V⟩ by about 4 within a quarter time unit. | answers step 10 open item 1 for this layer; C4 does not transfer |
| [Prop M9](fourd_compensated_ledger.md#7-proposition-m9-sub-cell-positions) | 7 | Measured: with sub-cell positions, cell-level matching biases ⟨V⟩ by +0.133 at t = 1; finer matching cells reduce it at 1.7 to 4.2 times the bodies. | sharpens N-SP2; bears on Y-SP1; M-SP6 |
| [Cor M10](fourd_compensated_ledger.md#8-corollary-m10-what-the-cost-scales-with) | 8 | The world form's cost scales with the occupied joint cells, a product over degrees of freedom that M8 prevents coarse binning from keeping small. | the exponential of step 10's B2, B3, returned as body inflation |

### Step 22 — [The sea as a phase reference](sea_phase_reference.md)

| ID | § | Says | Standing |
|---|---|---|---|
| [Thm L1](sea_phase_reference.md#1-theorem-l1-the-kernel-as-a-sum-over-contacts) | 1 | K_q(x) = −(B dp/ℏ)∫dy U_res w sin(2ξ_q y/ℏ): the kernel is a sum over contacts; B dp 2y_max = 1, one aligned pair per strip. | — |
| [Prop L2](sea_phase_reference.md#2-proposition-l2-contact-noise-and-the-measure) | 2 | Measured: contact sampling is unbiased but needs about 1.7 × 10⁵ contacts per unit time; weighting by the actual sea costs little (2.9 → 3.8 per cent); single contacts carry a force. | open item L-SP8 |
| [Thm L3](sea_phase_reference.md#3-theorem-l3-what-two-clocks-carry) | 3 | Two Lagrangian clocks at x ∓ y wind apart by the kernel's U(x, y); their midpoint misalignment winds at the energy difference. | visible only for same-row partners |
| [Thm L4](sea_phase_reference.md#4-theorem-l4-darkness-annihilates-signed-pair-sums) | 4 | Any ε-weighted pair sum vanishes on aligned pairs, so the sea can enter a vertex only as unsigned clocks. | open item L-SP6 |
| [Thm L5](sea_phase_reference.md#6-theorem-l5-dark-catalysis-and-maintenance) | 6 | Dark catalysis: an aligned pair fires at the kernel's rate with both members together; ΔE = ΔN = ΔS = 0, momentum conserved, the partner pair leaves aligned. | open item L-SP2 |
| [Prop L6](sea_phase_reference.md#6-theorem-l5-dark-catalysis-and-maintenance) | 6 | Measured: dark catalysis couples 87–95 per cent of the sea near the barrier; at reach scale it maintains a lock against event damage (0.470 against 0.461 with no events); co-located coupling does not (0.395). | answers Y-SP2 in part |
| [Thm L7](sea_phase_reference.md#7-theorem-l7-the-dephasing-of-a-locked-sea) | 7 | A freshly locked sea dephases at the trapezoid residual U − y[V′(x−y) + V′(x+y)]: force-free, but its third moment is about −2 times the compensated kernel's. | open item L-SP1 |

## 5. No-go and negative results

Results that say something cannot be done. Each is as valuable as a positive result, and several are the reason a later step exists.

| Where | What it rules out | Escape or consequence |
|---|---|---|
| Step 2 — [no-go lemma](four_rule_microdynamics_equivalence.md#62-the-mediation-is-relocated-not-eliminated) | Collision rates that are pairwise mass action (degree 2 or more in the occupancies) cannot give the linear QLE generator. | Needs a species whose density is pinned: the sea (step 3). |
| Step 4 — [Thm 2](phase_resonance_microdynamics.md#6-theorem-2-no-go-phase-blindness-cannot-be-linear-in-the-potential) | A phase-blind event rule cannot be linear in the potential. | Phase is necessary there; Y1 shows the compensated kernel is not covered. |
| Step 6 — [Cor R2.1](relational_pairing_and_carrier_lock.md#3-the-obstruction-a-free-carrier-per-pair) | An all-pairs relational average with a private carrier per pair vanishes as 1/B. | Postulate (S), then permanent pairing (step 7). |
| Step 6 — [Prop R5](relational_pairing_and_carrier_lock.md#8-correction-the-indexed-sea-is-a-consumable-resource) | Under indexed permanent partnership a pair mediates once: the sea is consumable (short by about 770 times). | Arithmetic stands; the inference R5.1 is retracted in step 7. |
| Step 7 — [no-go for lone hoppers](permanent_pairing_density_matrix.md#4-one-leg-hops-derive-the-stencil-and-the-mediated-counting) | Lone-particle hopping with position-only rates has irreducible spurious diffusion. | Pairs, not lone hoppers. |
| Step 9 — [Prop P4](position_pair_ladder.md#7-proposition-p4-two-obstructions), [Prop P5](position_pair_ladder.md#8-proposition-p5-resource-arithmetic-and-the-rate-comparison) | No positon-only sea in the position representation; the vertex rate has no continuum limit. | The observable sector still works (P7); the coherence sector does not. |
| Step 10 — [Prop B1](fourd_microdynamics.md#31-the-joint-sea-is-a-product-sea-but-the-shift-is-not-a-product), [Prop B2](fourd_microdynamics.md#32-the-sixteen-channels-lift-and-the-sea-stops-being-optional) | A joint sea is not a product of seas; excess–excess collisions vanish for N above 1. | The sea is the only collision partner. |
| Step 10 — [Thm C1](fourd_microdynamics.md#42-theorem-4-loses-its-uniqueness) | The exchange vertex is not unique in d of 2 or more. | Swap must be postulated (open item 1). |
| Step 11 — [Thm O1](open_position_space.md#2-the-wigner-kernel-has-no-position-envelope) | No absorbing layer can cure the escape of worlds: the kernel modulus is flat in position. | Horizon on the separation, not on position. |
| Step 11b — [Thm E7](reach_energy_coupling.md#7-a-polynomial-potential-has-no-jump-measure-at-all) | A polynomial potential has no jump measure on the open line. | A reach must be imposed. |
| Step 11b — [Thm E8](reach_energy_coupling.md#8-what-the-reach-buys-and-what-it-costs) | The compensated residual budget diverges with the reach for every V with V′ not zero. | Reach is a regulator (G5). |
| Step 12 — [Thm I4](interworld_coupling.md#4-theorem-i4-the-moyal-series-is-the-separation-expansion) | Harmonic and inverted-harmonic potentials have no jump channel; a coupling linear in Y is exactly classical. | Such potentials cannot exercise the hop channel; the Eckart barrier (step 15) is the first open-line test that does. |
| Step 13 — [Thm D0](species_sectors_and_annihilation.md#2-theorem-d0-the-two-ensembles-have-no-carrier-correspondence), [Prop U1](species_sectors_and_annihilation.md#7-proposition-u1-why-a-quadratic-channel-is-unavoidable) | No carrier-level map between the Wigner and density-matrix ensembles; no unravelling linear in the ensemble is L1-stationary. | Separate derivations per ensemble; stabiliser outside the linear class (S1). |
| Step 14 — [Thm C5](compensated_liouville_splitting.md#6-the-ring-as-a-diagnostic) | A ring cannot test the reach condition: the residual vanishes only for constant V. | Open-line test problems (steps 15, 19). |
| Step 14 — [Thm C9](compensated_liouville_splitting.md#31-the-residual-has-no-hop-rate-on-the-momentum-lattice) | No compensated member is a finite-range lattice rate law. | The momentum lattice survives only as a transfer spectrum (CLS3, CLA5). |
| Step 15 — [Cor K1.1](eckart_barrier_compensated.md#2-theorem-k1-the-reach-ceiling-is-the-analyticity-strip) | A reach-limited lattice cannot resolve the Eckart barrier's own momentum scale. | Sharpened in Z4; open as K-LS2 and CLA7. |
| Step 16 — [Thm S4](sea_population_equilibrium.md#5-the-emissive-unravelling-fails), [Thm S5](sea_population_equilibrium.md#5-the-emissive-unravelling-fails) | The emissive unravelling relocates the sea without bound; throttling by occupancy breaks the QLE. | Absorptive unravelling (S6–S9). |
| Step 17 — [Thm G2](compensated_ontology.md#3-g2-the-demographic-channel-is-empty-for-quadratic-hamiltonians) | For quadratic V creation and annihilation do nothing, yet the theory is quantum. | Admissibility (A) carries the rest (G3). |
| Step 17 — [Prop G3.1](compensated_ontology.md#5-postulate-a-is-independent-and-the-violated-inequality-is-the-wigner-bound), [Thm G4](compensated_ontology.md#6-g4-streaming-alone-leaves-the-admissible-set) | (A) cannot be derived from (S) + (D); (S) alone leaves the admissible set. | (A) is a separate postulate; the residual channel does kinematic work. |
| Step 18 — [Thm N6](stochastic_ledger.md#7-theorem-n6-transport-is-the-local-regulator) | Recombination cannot regulate the local ledger. | Transport does. |
| Step 19 — [Thm Z4](soft_core_coulomb.md#4-theorem-z4-the-softening-threshold) | An unsoftened atom cannot live on a uniform reach-limited crystal. | Non-uniform lattice (Z-LS1), deferred. |
| Step 20 — [Thm Y5](dark_sea_and_worldline_identity.md#6-theorem-y5-two-clocks-force-piecewise-worldlines), [Thm Y6](dark_sea_and_worldline_identity.md#7-theorem-y6-who-crosses-the-barrier) | World conservation with no momentum change makes the residual channel inert; tagged crossing tracks neither classical nor quantum transmission. | Kinks (piecewise worldlines) or birth and death; tunnelling is a statement about E, not identity. |
| Supplement — Prop T3 (takabayasi_1954_stochastic_picture) | The residual generator is never a one-body Markov jump generator (cited in the splitting note §2.2). | Signed ensemble. |
| Step 22 — [Thm L7](sea_phase_reference.md#7-theorem-l7-the-dephasing-of-a-locked-sea) | Reading the compensated kernel off the dephasing of a sea that streams under its own forces: it measures the partners' accelerations, not the parent's, and gets the quantum correction wrong by a factor of about −2. | Add the leg-force term at the vertex, or transport partner clocks in the parent's frame (L-SP1). |
| Step 22 — [Prop L2](sea_phase_reference.md#2-proposition-l2-contact-noise-and-the-measure) | Generating kernel weights by sampling sea contacts impulsively: unbiased, but about 10⁵ contacts per unit time are needed. | Phase accumulated continuously along worldlines (§3). |

## 6. Corrections and retractions

What each note corrected in an earlier statement, most of which the earlier statement does not record. The earlier rows in §4 carry a matching *Standing*.

| Note | Corrects | What changed |
|---|---|---|
| Step 2 | its own reconstruction from page 3 of 4 of the slide deck | Page 4 (Aug 2026) confirms it, including the centre-indexing: the symmetric member is the intended one, open item 1 is closed, and the slide-indexed member of §4(c) is not what was intended. The claim that momentum and energy balance determine the rates is audited in supplement/four_action_foundations.md §1. [§0](four_rule_microdynamics_equivalence.md#0-status-and-provenance). |
| Step 4 | step 3 (superelastic reading of the hop channels) | Superseded; open items 1, 3, 5 of step 3 advanced. See the update in [§12](sea_dressed_microdynamics.md#12-open-items). |
| Step 5 | its own first statement of the exchange rule | Roles of particle and struck partner had been reversed; settled by demo Part C. [§4.3](phase_alignment_microdynamics.md#43-a-caveat-on-translating-to-row-indices). |
| Step 6 | algorithm spec §2.2 (stored partner index); its own first revision | Partnership is not a carrier of state (R1); index arithmetic corrected July 2026; consumable-sea defect recorded. [§0](relational_pairing_and_carrier_lock.md#0-status-and-provenance), [§8](relational_pairing_and_carrier_lock.md#8-correction-the-indexed-sea-is-a-consumable-resource). |
| Step 7 | step 6 | Retracts two inferences, not calculations: 'pairing adds no information' (from R1) and Cor R5.1. Makes (S) unnecessary. [§0](permanent_pairing_density_matrix.md#0-what-this-note-inherits-and-retracts). |
| Step 8 | the revised algorithm spec, and step 5 Cor 4.2 | Five corrections: ladder climbing is leg-local; the sea–sea channel is not optional; the erase amplitude; the parity claim; Cor 4.2 holds for the exact (K3) channel only. [§0](coherence_ladder.md#0-what-this-note-inherits-retracts-and-corrects). |
| Step 9 | density_matrix_microdynamics_algorithm.md §2.3, §3–4, §7.1, §0 | A uniform bound does exist on the lattice; the pair-Bohm machinery is unnecessary; the position-side sign problem is structurally irremovable. Also: self-conjugate particles are not the movers. [§0](position_pair_ladder.md#0-what-this-note-inherits-retracts-and-corrects). |
| Step 10 | the 4-D supplement §6; step 2 §6.1; step 5 Thm 4 | Table organised along the wrong axis (mode wavevector decides); hop channel is a particle–particle exchange for pair-potential modes; Thm 4 is unique only in d = 1. [§0.1](fourd_microdynamics.md#01-corrections-recorded-here). |
| Step 11 | inverted_pair_barrier.md §7 item 4; two algorithm specs | No absorbing layer is needed (retracted); the box is one of three sources of the momentum quantum; the K-fold sub-lattice refinement interleaves copies. [§0](open_position_space.md#0-what-this-note-inherits-retracts-and-corrects). |
| Step 11b | step 11 Def (R), §3.2, §4.1 | Reach is a period, not an aperture; a commensurate sharp window is exact; the three mechanisms are one. Closes open items 1 and 4 of step 11; its own first draft's energy coupling is withdrawn. [§0](reach_energy_coupling.md#0-what-this-note-settles-retracts-and-corrects). |
| Step 12 | the hypothesis that prompted it | No non-linear interworld force is needed; the apparent four-wave mixing is the bilinearity of ρ. Corrects nothing in the repository. [§0](interworld_coupling.md#0-what-this-note-inherits-supplies-corrects-and-concedes). |
| Step 13 | supplement §7.4(ii), N1 (naming and mean-field inference), §6; its own earlier version | Anonymous annihilation is licensed (D5.1); the outward action is defocus; admissibility of the shift is D2. Replaces the earlier step 13 (species_phase_duality); the earlier Cor D4.1 is withdrawn. [§0.4](species_sectors_and_annihilation.md#04-what-this-note-corrects-elsewhere). |
| Step 14 | the previous revision (commit 71d6dce); the prompting question's framing | Inter-mode cancellation exists on a ring, not on the open line; negativity of W is not the sharp criterion of quantumness. TV table is grid-dependent (algorithm §4.4). [§0](compensated_liouville_splitting.md#0-what-this-note-inherits-and-corrects). |
| Step 15 | splitting §5.1; Thm C6 | For an exponential tail the reach rescales the interaction profile (K2); the Coulomb ceiling is a special case of K1. [§0.4](eckart_barrier_compensated.md#04-what-this-note-corrects), [§3](eckart_barrier_compensated.md#3-theorem-k2-the-far-field-and-an-erratum). |
| Step 16 | its own figures; Thm S7 (later) | Two demo defects (transport in x only; kernel ordering) repaired after the fact; then S7's f = 1/2 recognised as the sinkless case and S-SP3's argument withdrawn. [§0](sea_population_equilibrium.md#0-what-this-note-settles-and-what-it-corrects). |
| Step 17 | the slogan 'no new physics' | The correct slogan is 'no new force': (A) is non-dynamical and not derivable from (S) + (D). [§0](compensated_ontology.md#0-what-this-note-settles-and-what-it-corrects). |
| Step 18 | S7; S-SP3; the informal reading of κ | f = 1/2 is a special case (N3); the sign-change argument fails; recombination does not hold the local population steady (N6). [§0](stochastic_ledger.md#0-what-this-note-settles-and-what-it-corrects). |
| Step 19 | K7's 'exact zeros' (a gloss, not an error) | Nodes of V‴ are not nodes of Γ at finite reach; K-LS7 should carry a depth. [§2](soft_core_coulomb.md#2-theorem-z2-the-nucleus-is-dark-and-the-ring-survives). |
| Step 20 | Eckart §8.3; step 18 §1; step 5 §6; demo_emission_and_absorption.py | Individual crossing happens under piecewise worldlines; 'diffusion identically zero' holds only under birth and death; dark creation is a rule, not a consequence; the population clamp breaks E. [§0](dark_sea_and_worldline_identity.md#0-what-this-note-settles-and-what-it-corrects), [§8](dark_sea_and_worldline_identity.md#8-a-defect-in-the-mean-field-demo). |
| Step 21 | step 18 N4, N5; step 16 S9; step 10 C4 (transfer only) | The availability closure does not describe the ledger under real streaming (M7); S9's attractor is reached only after the population fills its reachable cells (M5, M6); step 10's transverse freedom is visible in the compensated layer (M8). [§0](fourd_compensated_ledger.md#0-what-this-note-settles-and-what-it-corrects) |
| Supplement R-SP5 | step 11b was missing from the ladder | Resolved: added to the ladder as step 11b in the same patch as this index. |
| Step 22 | step 20 §2 and abstract (Y1 wording); step 20 §6 (the two recombinations) | The ladder μ is not 'integrated out' into the kernel: the kernel holds the potential's rule for winding it, and arg ρ survives as the sign structure of W. The contact (κ) recombination of §6 is not part of the model (ORIENTATION). [§0](sea_phase_reference.md#0-what-this-note-asks-settles-and-leaves-open). |

## 7. Open items

Labelled items link to their line; the unlabelled items of the older notes are given by number, as those notes cite them. The section heading of each note's list is linked. Equivalent items under different IDs are marked `=`.

### Step 2 — [Four-rule microdynamics](four_rule_microdynamics_equivalence.md#10-open-items)

| ID | Item | Status | See |
|---|---|---|---|
| 1 | Page 4 of the slide deck: second bias rate and indexing | resolved (Aug 2026) | — |
| 2 | Variance-optimal member of the G-family | open | — |
| 3 | Self-consistent N-body case | open | — |
| 4 | Fermionic or multi-DOF extension | open | — |

### Step 3 — [Sea-dressed microdynamics](sea_dressed_microdynamics.md#12-open-items)

| ID | Item | Status | See |
|---|---|---|---|
| 1 | Level 3: self-consistent sea | advanced by step 4 | steady state still open (step 4 item 2) |
| 2 | Quantify the Bose bilinear (1/𝒩 scaling) | open | — |
| 3 | Excitation momentum ledger | advanced by step 4 | row-resolved form; step 6 Z_r |
| 4 | Variance structure (pinned run pays 1.45 times) | open | — |
| 5 | Streaming conventions for composites | advanced by step 4 | — |
| 6 | Canonical identity unravelling (tagged-particle demo) | begun in step 20 | Y6 |

### Step 4 — [Phase resonance](phase_resonance_microdynamics.md#11-open-items)

| ID | Item | Status | See |
|---|---|---|---|
| 1 | Reduce P5 | resolved in §13 | residue: Born resolution of a binary contact |
| 2 | Steady state (L3 proper) | open | — |
| 3 | Row-resolved ledger | refined in step 6 | Z_r, Thm R4 |
| 4 | Amplitude bookkeeping at finite B | open | — |
| 5 | Multi-mode potentials | open | also step 5 item 3, step 6 item 4 |
| 6 | Gray pairs | open | — |

### Step 5 — [Phase alignment](phase_alignment_microdynamics.md#11-open-items)

| ID | Item | Status | See |
|---|---|---|---|
| 1 | Row-index translation of Thm 4 | open | slide convention since confirmed (step 2) |
| 2 | Misalignment ledger | revised by step 6 | Z_n order parameter |
| 3 | Multi-mode superposition | open | — |
| 4 | Decoherence signatures | open | — |
| 5 | Uniformly offset pairs | revised by step 6 | — |

### Step 6 — [Relational pairing and carrier lock](relational_pairing_and_carrier_lock.md#11-open-items)

| ID | Item | Status | See |
|---|---|---|---|
| 1 | Species bookkeeping in Z | open | step 7 proposes ket/bra |
| 2 | Phase continuity at the vertex | open | assumed in step 8 Thm C2 |
| 3 | Dynamics of (S) | dissolved if permanent pairing is adopted | step 7 |
| 4 | Multi-mode potentials | open | — |

### Step 7 — [Permanent pairing, density-matrix reading](permanent_pairing_density_matrix.md#7-consequences-and-open-items)

| ID | Item | Status | See |
|---|---|---|---|
| 1 | (S) becomes unnecessary; does the spec §2.2 revert? | decision pending | — |
| 2 | Species = ket or bra side of ρ | adopted for the E2 ensemble | step 13 D0, D13 |
| 3 | Split pairs mediate with the pump-excited vertex constant | open (load-bearing) | step 8 §8 item 2 |
| 4 | Pair-ensemble Monte Carlo with permanent partners | open | sharper target in step 8 |

### Step 8 — [Coherence ladder](coherence_ladder.md#8-consequences-and-open-items)

| ID | Item | Status | See |
|---|---|---|---|
| 2 | Striker back-reaction neutrality; compound-channel cancellation; licensing lemma; bare-churn balance at V = 0; steady state | open | — |
| 3 | Reproduce Thm C2 fluxes stochastically | open | — |

### Step 9 — [Position-pair ladder](position_pair_ladder.md#12-consequences-and-open-items)

| ID | Item | Status | See |
|---|---|---|---|
| 3 | The number-conserving vertex | open | — |
| 4 | Is the divergence of Λ physical or a nearest-neighbour artefact? | open | — |
| 5 | Compact-momentum reading (j proportional to sin μ) | open | — |
| 6 | Species = ket or bra | second vote for ket/bra | adopted for E2 in step 13 |

### Step 10 — [Four-dimensional phase space](fourd_microdynamics.md#8-open-items)

| ID | Item | Status | See |
|---|---|---|---|
| 1 | Postulate (X): exchange-only, or transverse family physical | open for the phase-alignment layer; answered for the compensated ledger by step 21 M8 | live 4-D run |
| 2 | Directional sea ledger | open | subsumes step 5 item 3, step 3 item 3 |
| 3 | Joint sea instantiation | open | — |
| 4 | Classical-drift/quantum-remainder splitting | developed in step 14 | variance measurement open |
| 5 | Promote the exact jump substep to the library | open | — |
| 6 | Vector potentials | open | — |

### Step 11 — [Open position space](open_position_space.md#9-open-items)

| ID | Item | Status | See |
|---|---|---|---|
| 1 | A tapered horizon | resolved in step 11b | don't taper (E1.1, E2) |
| 2 | The sea per unit length | open | cf. K-LS8 |
| 3 | Half-jump or full-jump (modulus πℏ/a or 2πℏ/a) | open | — |
| 4 | Energy under the horizon | resolved in step 11b | E3 |
| 5 | Sampling the θ integral | open | — |
| 6 | Two horizons in the two-body problem | open | — |

### Step 11b — [What the reach controls](reach_energy_coupling.md#10-open-items)

| ID | Item | Status | See |
|---|---|---|---|
| 1 | Open momentum lattice: E4 with capped end cells | open | — |
| 2 | Profile with vanishing second derivative at the edge (q^−4 tail) | open | = R-SP4 |
| 3 | Live quartic double-well benchmark | open | — |
| 4 | Variance, not event count | open | — |
| 5 | Folding under the compensated split | open | — |

### Step 12 — [Interworld coupling](interworld_coupling.md#9-consequences-and-open-items)

| ID | Item | Status | See |
|---|---|---|---|
| 1 | Is the coupling U unique? | open | — |
| 2 | Pair stability | open | (S) looked load-bearing; see step 7 |
| 3 | The continuum case | open | horizon of step 11 is the candidate |
| 4 | Third-moment test for a mechanical model | open | cf. E5, E6 |

### Step 13 — [Species, sectors, annihilation](species_sectors_and_annihilation.md#13-open-items)

| ID | Item | Status | See |
|---|---|---|---|
| 1 | The E1/E2 relationship | open (structural) | R-SP8 |
| 2 | Adaptive scheme's lookahead | open | — |
| 3 | Plateau test with transport | open | — |
| 4 | The third regularisation | open | — |
| 5 | Cost in the pair basis (conjecture) | open | — |
| 6 | Amplitude of the near-dark tail | open | — |
| 7 | Many bodies | open | — |

### Step 14 — [Compensated Liouville splitting](compensated_liouville_splitting.md#9-open-items)

| ID | Item | Status | See |
|---|---|---|---|
| [CLS1](compensated_liouville_splitting.md#9-open-items) | Best quadratic fit on the reach interval instead of the osculating quadratic | open | — |
| [CLS2](compensated_liouville_splitting.md#9-open-items) | Mean-field versus shot-noise reading of the split | half settled by C8 | — |
| [CLS3](compensated_liouville_splitting.md#9-open-items) | Sub-Δp acceleration on the momentum lattice | resolved by C9 | = CLA5 (lattice-resident version open) |
| [CLS4](compensated_liouville_splitting.md#9-open-items) | What sets L_c physically | open | = CLA1 |
| [CLS5](compensated_liouville_splitting.md#9-open-items) | Erratum line for supplement §6.1 (stencil is an identity) | open | — |
| [CLS6](compensated_liouville_splitting.md#9-open-items) | Quiet region of positive measure without V‴ vanishing | open | — |
| [CLS7](compensated_liouville_splitting.md#9-open-items) | Does the horizon restore locality for the two-body kernel | open | — |
| [CLS8](compensated_liouville_splitting.md#9-open-items) | Momentum-space locality: touching a cell versus consuming its object | open | — |

### Step 15 — [Eckart barrier](eckart_barrier_compensated.md#10-open-items)

| ID | Item | Status | See |
|---|---|---|---|
| [K-LS1](eckart_barrier_compensated.md#10-open-items) | Flux ledger closes only to 3.8 per cent | open | — |
| [K-LS2](eckart_barrier_compensated.md#10-open-items) | Reach-limited transmission: ceiling and resolution conditions seem incompatible | open | sharpened by Z4, Z5; = CLA7 |
| [K-LS3](eckart_barrier_compensated.md#10-open-items) | Scaling of the K6 constant with packet offset and width | open | — |
| [K-LS4](eckart_barrier_compensated.md#10-open-items) | Initial state with genuine negativity | open | — |
| [K-LS5](eckart_barrier_compensated.md#10-open-items) | Soft-core Coulomb as the next test case | answered by step 19 | — |
| [K-LS6](eckart_barrier_compensated.md#10-open-items) | Lobe count versus cancellation depth | open | Z1 gives a second lobe count |
| [K-LS7](eckart_barrier_compensated.md#10-open-items) | Does the sea inherit the lobes (exact holes)? | open | restated by Z-LS5; with N-SP4 |
| [K-LS8](eckart_barrier_compensated.md#10-open-items) | The box: sea scales with window volume | open | = CLA8 |

### Step 16 — [Sea population equilibrium](sea_population_equilibrium.md#8-open-items)

| ID | Item | Status | See |
|---|---|---|---|
| [S-SP1](sea_population_equilibrium.md#8-open-items) | The supra-minimal ensemble | resolved, against the prediction | S9; constant fixed by N5 |
| [S-SP2](sea_population_equilibrium.md#8-open-items) | Tension with the sea-dressed note's live ledger | open | — |
| [S-SP3](sea_population_equilibrium.md#8-open-items) | Does the shortfall vanish? | resolved by measurement; one argument withdrawn | N3; direct test is N-SP1 |
| [S-SP4](sea_population_equilibrium.md#8-open-items) | Mean-field caps and channel ordering | resolved, with a caveat | ledger is order-dependent |
| [S-SP5](sea_population_equilibrium.md#8-open-items) | A genuine ensemble against the mesh | open | = CLA2 |
| [S-SP6](sea_population_equilibrium.md#8-open-items) | S4 on a bound state | vehicle built in step 19 | = Z-LS4 |
| [S-SP7](sea_population_equilibrium.md#8-open-items) | Negative-cap leakage | open | cf. Y-SP6 |

### Step 17 — [Compensated ontology](compensated_ontology.md#10-open-items)

| ID | Item | Status | See |
|---|---|---|---|
| [G-SP1](compensated_ontology.md#10-open-items) | Split gauge invariance | open (load-bearing) | = CLA10 |
| [G-SP2](compensated_ontology.md#10-open-items) | Particle-level statement of (A) | open | — |
| [G-SP3](compensated_ontology.md#10-open-items) | Quantify G4 | open | — |
| [G-SP4](compensated_ontology.md#10-open-items) | Finite-dimensional systems (spin) | open | — |
| [G-SP5](compensated_ontology.md#10-open-items) | Soft-horizon variant of the fold | open | prerequisite for a clean G5 |

### Step 18 — [Stochastic ledger](stochastic_ledger.md#9-open-items)

| ID | Item | Status | See |
|---|---|---|---|
| [N-SP1](stochastic_ledger.md#9-open-items) | Instrument the sum rule in the mesh demo | open (highest value) | settles S-SP3 exactly |
| [N-SP2](stochastic_ledger.md#9-open-items) | The dictionary between lattice occupancy and Wigner units, hence κ | open | J-SP3; sharpened by step 21 M9 |
| [N-SP3](stochastic_ledger.md#9-open-items) | Measure the Fano factor properly | open | — |
| [N-SP4](stochastic_ledger.md#9-open-items) | Transport that is transport | open | with K-LS7; answered in part by step 21 M7 |
| [N-SP5](stochastic_ledger.md#9-open-items) | The two legs are not independent | open | J-SP1 |

### Step 19 — [Soft-core Coulomb](soft_core_coulomb.md#8-open-items)

| ID | Item | Status | See |
|---|---|---|---|
| [Z-LS1](soft_core_coulomb.md#8-open-items) | The non-uniform lattice | deferred by decision | — |
| [Z-LS2](soft_core_coulomb.md#8-open-items) | Transmission against a split-operator reference | planned next | — |
| [Z-LS3](soft_core_coulomb.md#8-open-items) | The lobe ledger | open | — |
| [Z-LS4](soft_core_coulomb.md#8-open-items) | The ground state as a stationary ledger | open | = S-SP6 |
| [Z-LS5](soft_core_coulomb.md#8-open-items) | Quiet points are minima, not nodes | open | restates K-LS7 |
| [Z-LS6](soft_core_coulomb.md#8-open-items) | Repulsive core and scattering | open | — |

### Step 20 — [Dark sea and worldline identity](dark_sea_and_worldline_identity.md#10-open-items)

| ID | Item | Status | See |
|---|---|---|---|
| [Y-SP1](dark_sea_and_worldline_identity.md#10-open-items) | Does a world-particle carry a position below cell resolution? | open | Z-LS1; measured consequence in step 21 M9 |
| [Y-SP2](dark_sea_and_worldline_identity.md#10-open-items) | Creation rule: forced or free | answered in part | step 22 §6: re-locking makes dark creation and lock maintenance one rule (free branch) |
| [Y-SP3](dark_sea_and_worldline_identity.md#10-open-items) | Is the ±ξ²/m per event the same 1/2 as N3's? | open | — |
| [Y-SP4](dark_sea_and_worldline_identity.md#10-open-items) | Does T_tag converge as dp → 0 and as the reach grows? | open | — |
| [Y-SP5](dark_sea_and_worldline_identity.md#10-open-items) | Repeat §7 on the exact integer ledger | open | — |
| [Y-SP6](dark_sea_and_worldline_identity.md#10-open-items) | The clamp of §8 | open | S-SP7 |
| [Y-SP7](dark_sea_and_worldline_identity.md#10-open-items) | Phase-gated allocation changes T_tag and leaves E exactly | open | — |
| [Y-SP8](dark_sea_and_worldline_identity.md#10-open-items) | Vocabulary: 'aligned' for 'bound'; a rename for the carrier lock | open | R-SP6 |

### Step 21 — [Compensated ledger in four dimensions](fourd_compensated_ledger.md#10-open-items)

| ID | Item | Status | See |
|---|---|---|---|
| [M-SP1](fourd_compensated_ledger.md#10-open-items) | The open line: any ledger equilibrium for a bound problem when escaped bodies cannot return? | open | CLA8 |
| [M-SP2](fourd_compensated_ledger.md#10-open-items) | Contact recombination for the escaped population: revisit its exclusion? | open | ORIENTATION §8 |
| [M-SP3](fourd_compensated_ledger.md#10-open-items) | The Eckart pair collision through the ledger, against Part E's verified reference | open | — |
| [M-SP4](fourd_compensated_ledger.md#10-open-items) | The scan at a reach where M4's ratio is below one | open | CLA7 |
| [M-SP5](fourd_compensated_ledger.md#10-open-items) | A non-separable problem; bins in particle coordinates | open | — |
| [M-SP6](fourd_compensated_ledger.md#10-open-items) | A matching rule that removes M9's bias without multiplying the population | open | Y-SP1, N-SP2 |
| [M-SP7](fourd_compensated_ledger.md#10-open-items) | Fidelity per unit cost as a function of ν once cells set the population | open | — |
| [M-SP8](fourd_compensated_ledger.md#10-open-items) | Event-ordering sensitivity of the body counts | open | S-SP4 |

### Step 22 — [The sea as a phase reference](sea_phase_reference.md#8-open-items)

| ID | Question | Status | Related |
|---|---|---|---|
| [L-SP1](sea_phase_reference.md#8-open-items) | Can partner clocks be referred to the parent's frame, or is the leg-force term a vertex postulate? | open | Thm L7 |
| [L-SP2](sea_phase_reference.md#8-open-items) | Reach-scale dark catalysis: how is the partner chosen, and what is its standing? | open | Thm L5 |
| [L-SP3](sea_phase_reference.md#8-open-items) | The particle model's events fail 20–45 per cent of the time and E is unchecked; repeat on the exact ledger. | open | Y-SP5 |
| [L-SP4](sea_phase_reference.md#8-open-items) | Does co-located coupling coarsen to the reach scale from a gas at long times? | open | Prop L6 |
| [L-SP5](sea_phase_reference.md#8-open-items) | Is a locked sea an admissible initial state, a boundary condition? | open | Def L0 |
| [L-SP6](sea_phase_reference.md#8-open-items) | The sea at the vertex as unsigned clocks while dark as amplitude: postulate or reading of P2? | open | Thm L4 |
| [L-SP7](sea_phase_reference.md#8-open-items) | Is any locked reference maintainable inside a potential under (S)? | open | Prop L6, Thm L7 |
| [L-SP8](sea_phase_reference.md#8-open-items) | Σ E drifts to 1.0010 in the sea-weighted kernel run although every event deposits ±1. | open | Prop L2 |

## 8. Defects found in demos

Defects in code or specifications, found in the course of the analysis, whose effect on quoted numbers matters.

| Where | Defect | Status | Recorded |
|---|---|---|---|
| demo_sea_population_equilibrium.py | stream() advected in x only, so postulate (S) was never exercised; the residual kernel applied the horizon after subtracting 2yV′, returning the horizon's own first moment (3.6 × 10⁻⁵ against 2 × 10⁻¹⁵). | repaired; figures in step 16 predate the repair; S7, S2, S4 unaffected, S8 improved | step 16 §0 |
| demo_emission_and_absorption.py (Ledger.channels) | The population clamp breaks E: Σ E falls from 1 to 0.8402 in 250 steps (22.0 per cent field error against 3.5 per cent without) and biases f from about 0.64 to 0.52. | open item Y-SP6; proposal J-SP2 (clip caps, repair ringing at transport) | step 20 §8 |
| Algorithm spec §2.2 (indexed partners) | Under permanent partnership the sea is a consumable with no source, short by about 770 times for the cosine-well parameters. | arithmetic confirmed in step 7; the inference retracted | step 6 §8 |
| Splitting note §4 (TV table) | Total-variation figures are functions of the rung grid, not absolute: under a hard horizon the event rate diverges logarithmically and the momentum churn linearly. | erratum in the note; see algorithm spec §4.4 | step 14 §4 |

## 9. Cited from outside this folder

Labels that analysis notes cite but that are defined in [`../supplement/`](../supplement/README.md) or [`../algorithm/`](../algorithm/README.md). Results first.

| Label | Where | Says | Relation |
|---|---|---|---|
| J1, J2 | supplement/emission_and_absorption.md §5 | P = S + N/2 is conserved exactly by the event channel, but only globally (cell by cell it is violated). | N1 strengthens J1 to pathwise |
| F3 | supplement/four_action_foundations.md | Endpoint locality selects the symmetric member of the FR family: the four-action model. | cited in step 10 (A1) |
| T3 | supplement/takabayasi_1954_stochastic_picture.md | No one-body Markov reading: a generator with a signed off-diagonal kernel is not a Markov jump generator. | reappears unweakened in the residual channel (splitting §2.2) |

Open items. Those marked `=` are the same question as an analysis-note item under another ID.

| ID | Where | Item | Status | See |
|---|---|---|---|---|
| CLA1 | algorithm spec | What sets the horizon profile, not just its width | open; load-bearing there | = CLS4 |
| CLA1b | algorithm spec | Restate Def (H) as a rung occupancy if CLA1 resolves for a profile | open | Def (H), O4, C7 |
| CLA2 | algorithm spec | World-form validation against the mesh at fixed variance | open | = S-SP5 |
| CLA3 | algorithm spec | Annihilation load | resolved: absorption is the pair-removal channel (S9); the constant is unfixed | S-SP1, N5 |
| CLA4 | algorithm spec | Time-dependent potential: does the active set still pay | open | — |
| CLA5 | algorithm spec | Sub-Δp acceleration; strictly lattice-resident version | partly answered | = CLS3, C9 |
| CLA6 | algorithm spec | Multi-dimensional symbol construction | open | — |
| CLA7 | algorithm spec | The reach is over-determined (K1 ceiling, ledger wants about 4πa, K1.1) | open | K-LS2, Z4 |
| CLA8 | algorithm spec | The box: sea versus unpaired population | open | = K-LS8 |
| CLA9 | algorithm spec | Deprecated M± hop channel: a record, not a question | record | (S) |
| CLA10 | algorithm spec | Does the census mean anything | open | = G-SP1 |
| J-SP1 | supplement/emission_and_absorption.md | Conjunction cost measured, explanation not | open | N-SP5 |
| J-SP2 | supplement/emission_and_absorption.md | Change to the allocation loop (clip caps) | confirmed independently in step 16 Part J | Y-SP6 |
| J-SP3 | supplement/emission_and_absorption.md | Sea depth is an initial condition | open | N-SP2 |
| J-SP4 | supplement/emission_and_absorption.md | Two vertices are both called absorption | answered in part, step 20 | Y1, Y3 |
| R-SP1 | supplement/what_the_reach_is.md | Two ceilings; effect on K1.1 and Z4 | open; not applied by decision | CLA7 |
| R-SP2 | supplement/what_the_reach_is.md | Is the diagonal contamination the same phenomenon as the period criterion | open | — |
| R-SP3 | supplement/what_the_reach_is.md | Does the fold reach the state or only the potential | open | — |
| R-SP4 | supplement/what_the_reach_is.md | Profile with a q^−4 tail | open | = item 2 of step 11b |
| R-SP5 | supplement/what_the_reach_is.md | reach_energy_coupling.md is not in the ladder | resolved: step 11b | — |
| R-SP6 | supplement/what_the_reach_is.md | Vocabulary collisions (particle, pair) | open; nothing renamed | Y-SP8 |
| R-SP7 | supplement/what_the_reach_is.md | Sea side of N6 is untested | open | N6 |
| R-SP8 | supplement/what_the_reach_is.md | The two pictures do not share an inventory | open | step 13 item 1 |
| H-SP1, H-SP2, H-SP3 | supplement/holland_two_fluid_correspondence.md | The σ → 0 limit versus S5's throttle; uniqueness of the reservoir; free Gaussian as a test of H7 | open | — |
| W-SP1, W-SP2 | supplement/limkumnerd_weighted_paths.md | Relation of A_sign to f; the χ_Q(t) residual diagnostic | open | — |

## 10. Maintaining this index

This file is updated in the same patch as the note it indexes, in parallel with that note's ladder entry in the [README](README.md) and the blockquote under its title. The check is:

```bash
PYTHONPATH=src python3 -m wpmwlib.check_index docs/analysis
PYTHONPATH=src python3 -m wpmwlib.check_index docs/analysis --suggest
```

It finds every labelled result and every labelled open item in the notes and fails if one has no row here, if any relative link in this file points at a missing file or heading, and warns about a row that matches nothing. `--suggest` prints skeleton rows, with section number and anchor filled in, for anything missing. It does not read the *Says* or *Standing* columns; those are written by hand and are the content.

What a patch does to this file:

- **New result or open item in a note:** add its row (§3 for a postulate or definition, §4 for the rest, §7 for an open item), starting from `--suggest`.
- **A note corrects, restricts or retracts an earlier result:** add the correction to §6 and record it in the earlier row's *Standing*.
- **An open item is closed, answered or superseded:** change its status in §7, and in §9 if it has an equivalent there.
- **A new note:** add it to the ladder (README and title blockquote), to the registry in §2 if it introduces a prefix, and give it its headings in §3, §4 and §7.
- **A defect in a demo or specification:** add it to §8.

Sections 2, 5, 6, 8 and 9 are hand-maintained. Their links are checked; their content is not.
