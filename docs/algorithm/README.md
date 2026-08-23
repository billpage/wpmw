# docs/algorithm

Core algorithm specifications for the WPMW project.

This directory holds the canonical descriptions of the algorithms used in `src/`:
the signed-particle / phase-space-crystal-lattice update rule, the split-Fourier
collision step, and the discrete momentum-jump scheme for sinusoidal and
polynomial potentials.

Specifications here should be precise enough that someone reading them could
re-implement the algorithm in another language. Background discussion and
extended derivations belong in `docs/supplement/` and `docs/analysis/`.

## Contents

- **[`phase_space_crystal_lattice_algorithm.md`](phase_space_crystal_lattice_algorithm.md)** — The canonical 1+1-D
  specification: mesh-density and Monte-Carlo particle forms of the single
  mediated-jump rule, the Fourier-mode and differential jump forms, the
  momentum half-grid `dp = πħ/L`, and the `W' = W + 2/h` background shift.
  Uses the QLE-consistent sign convention; see §6.3 of
  [`../supplement/phase_space_crystal_lattice_supplement.md`](../supplement/phase_space_crystal_lattice_supplement.md) for the algebraic
  step where the source memo lost it. Implemented in
  `src/wpmwlib/phase_space_crystal_lattice.py`. Items marked **[choice]** are
  implementation decisions not fixed by the source documents.
- **[`multi_body_extension.md`](multi_body_extension.md)** — Forward-looking extension of the 1+1-D spec
  to `d` spatial dimensions (a direct vectorial generalisation, cost
  `(M_x M_p)^d`) and to `N` interacting particles, where the mesh becomes
  infeasible and the world-ensemble representation is the only viable form.
  Includes the two-body Moyal-bracket derivation, the Ewald treatment of
  Coulomb, a validation sequence, and open conventions in §12.
- **[`density_matrix_microdynamics_algorithm.md`](density_matrix_microdynamics_algorithm.md)** — Forward-looking
  specification of the position-space dual: the von Neumann equation
  represented by complex-weighted world-particle *pairs* undergoing pair-Bohm
  flow plus phase accumulation on the weight. §7 gives the term-by-term
  duality table against the Wigner-side algorithm in the `p ↔ s = x − x'`
  variable. Relates the construction to stochastic unraveling of the
  Liouville–von Neumann equation (Stockburger & Grabert 2002) and states what
  is original here relative to that literature.
- **[`phase_alignment_microdynamics_algorithm.md`](phase_alignment_microdynamics_algorithm.md)** — Specification of the
  **live sea**: pairs are dynamical, they carry phases, and the rate field
  `Γ_q(x)` is an output assembled from vertex statistics rather than an input.
  Deliberately more expensive and less accurate than the mesh algorithm for
  any problem the mesh algorithm can solve; its purpose is to make the
  ontology executable, so that claims about what the model *is* can be tested
  rather than asserted. Revised twice in July 2026: first to a relational
  formulation under postulate (S)
  ([`../analysis/relational_pairing_and_carrier_lock.md`](../analysis/relational_pairing_and_carrier_lock.md)), then to **permanent
  pairing** under the density-matrix reading
  ([`../analysis/permanent_pairing_density_matrix.md`](../analysis/permanent_pairing_density_matrix.md)), which withdraws (S),
  splits the vertex into write and erase channels, and predicts — rather than
  codes — the mediated QLE generator. §6 (the steady state under continuous
  pumping) remains **[open]**. Companion code `src/demo_phase_alignment.py`
  verifies §§2–5 in isolation and
  `src/demo_pairing_resource_arithmetic.py` the storage-capacity arithmetic;
  a full live-sea integrator is not yet in `wpmwlib`, and §10 rungs 5–6
  define the decisive test it must pass.

- **[`compensated_liouville_algorithm.md`](compensated_liouville_algorithm.md)** — Wigner evolution as exact
  Newtonian flow plus a zero-mean signed hop channel: the crystal-lattice jump
  rule of §3b above with the rate field replaced and a deterministic step in
  front of it. Written for the **open line**; the ring appears only as the
  corner where the algorithm degenerates back to the crystal lattice. Promotes
  [`../analysis/compensated_liouville_splitting.md`](../analysis/compensated_liouville_splitting.md)
  (theorems C0–C7) and adds what a grid forces: the reach *is* the momentum
  cell and `N_p` is the ket–bra rung count, so neither is a convergence knob
  (§2); the deterministic force must be the kernel's own first moment, not
  `V′(x)` (§3.2); the Nyquist rung must be zeroed (§4.1); and the coherence
  horizon needs a **profile**, since under a hard cutoff the event rate, the
  momentum churn and the third moment all diverge — a soft horizon is what
  gives the discrete operator a Moyal expansion at all (§4.4). §8 sketches the
  `wpmwlib` interface. Companion code `src/demo_compensated_liouville_algorithm.py`.

## Reading order

[`phase_space_crystal_lattice_algorithm.md`](phase_space_crystal_lattice_algorithm.md) alone is enough to reproduce every
result in `src/`. The other four extend it in four independent directions —
more particles, a dual representation, a deeper microdynamics, and a
reordering of the same update into a classical and a quantum channel — and
none is a prerequisite for the others.
[`compensated_liouville_algorithm.md`](compensated_liouville_algorithm.md) is
the closest to the canonical spec and reuses its jump rule verbatim, so read
that one second if you want the shortest path from the lattice to the open
line.
