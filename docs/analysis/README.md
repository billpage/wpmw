# docs/analysis

Mathematical derivations and detailed review notes.

This directory holds the worked-out mathematics behind the algorithms. It
began as review notes on the source documents — term-by-term matches between
the extended Fokker–Planck (xFP) expansion and the Moyal series for the Wigner
Quantum Liouville Equation, and moment-problem analyses for jump densities
under polynomial potentials. Most of it is now a **derivation ladder**: a
sequence of notes, each taking as input something its predecessor postulated,
and each ending with the open items that motivate the next.

## The ladder

Read in this order. Each note states in its §0 what it inherits and what it
retracts.

1. **`phase_space_crystal_lattice_review.md`** — Review of the two source
   documents (Cyganski's *Extended Fokker–Planck Eq. and the QLE V2* memo and
   the *Wigner Collisions Diagram* Sozi deck), cross-referenced at the
   equation-and-page level. The entry point for anyone tracing a claim back to
   its origin.
2. **`four_rule_microdynamics_equivalence.md`** — Analysis of Cyganski's
   proposal (Zoom, 2026) to replace the single mediated-jump rule with four
   two-body rules (Focus, Defocus, Right-Hop, Left-Hop). Proves exact
   equivalence at any particle number `ν`, shows the four-rule form is ≈5.6×
   quieter, and identifies the `G`-freedom — a family of exact rate
   assignments of which the single rule is the `G = 0` member. Ends on a no-go
   lemma: pairwise mass action among tracked particles is quadratic in
   occupancy while the QLE generator is linear, so a fully collision-based
   microdynamics needs a species whose density is *pinned*.
3. **`sea_dressed_microdynamics.md`** — Takes the step that lemma leaves open.
   Realises the collision term as sixteen local, two-body,
   momentum-conserving channels against a pinned Dirac sea of positon–negaton
   pairs, exact at pinned sea. Postulates the sea's polarisation: the rate
   field `Γ_q(x)`, its sign structure, and the half-quantum stencil offsets
   all enter as assumptions.
4. **`phase_resonance_microdynamics.md`** — Derives that polarisation rather
   than postulating it, by making phase a particle-level property (P0–P5,
   Theorems 1–3). Contains the parity result (fundamental particles occupy
   even momentum sites), the rate-table no-go (Theorem 2: phase-blind
   transition rules cannot reproduce linear rates), and the dark-sea lemma.
5. **`phase_alignment_microdynamics.md`** — A change of variables on the
   predecessor: the beat, the grating and the resonance condition are replaced
   by a single scalar, the misalignment `μ` of two transported clock phases.
   No new postulates and no different predictions, but Theorem 4 is stronger
   than what it replaces: requiring `μ` to hold still through a vertex forces
   the vertex to be a **momentum swap**, from which energy conservation and
   the selection rule follow rather than being imposed.
6. **`relational_pairing_and_carrier_lock.md`** — Removes stored partnership
   from the algorithm specification (§2.2), at the cost of one postulate (S),
   the sea carrier lock. Proposition R1 shows a partner index carries no
   relational state; Theorem R4 factorises the vertex weight through a
   per-cell, per-row order parameter `Z_r`, cutting the encounter loop from
   `O(N_exc · B)` to `O(N_exc + N_sea)`; and §8 records a defect the indexed
   formulation concealed — under permanent partnership the sea is a
   consumable resource with no source, short by ≈770× for the cosine-well
   parameters.
7. **`permanent_pairing_density_matrix.md`** — Reinstates permanent pairing
   under a density-matrix reading: a pair is a sampled element of ρ (positon =
   ket leg, negaton = bra leg, μ = arg ρ), an excess particle a diagonal
   sample. Retracts two *inferences* of its predecessor while keeping its
   calculations: pairing adds nothing to the one-particle marginal but is the
   entire content of the two-point function; and the pump — whose sidebands
   *are* split-pair amplitudes — is the source Corollary R5.1 declared
   absent. Derives the stencil and the mediated counting from one-leg hops of
   bound pairs, shows storage feasibility is exactly the Wigner bound
   |W| ≤ 2/h (the same inequality as W' ≥ 0), makes postulate (S)
   unnecessary, and leaves one load-bearing theorem: split pairs mediate with
   the same vertex constant as pump-excited pairs.

## Companion code

Every note lists its verifying script in §0 and its numerical results in a
late section. In ladder order: `src/demo_four_rule_equivalence.py`,
`src/demo_sea_dressed_dynamics.py`, `src/demo_phase_resonance_rates.py`,
`src/demo_phase_alignment.py`, `src/demo_relational_pairing.py`,
`src/demo_pairing_resource_arithmetic.py`. All
non-trivial claims in these notes are verified numerically before inclusion.

## Related

Drawings of the elementary processes of the phase-alignment layer are in
`../supplement/phase_alignment_interaction_diagrams.md`; the simulation
counterparts of the ladder are in `../algorithm/`.
