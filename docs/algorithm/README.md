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

## Reading order

[`phase_space_crystal_lattice_algorithm.md`](phase_space_crystal_lattice_algorithm.md) alone is enough to reproduce every
result in `src/`. The other three extend it in three independent directions —
more particles, a dual representation, and a deeper microdynamics — and none
is a prerequisite for the others.
