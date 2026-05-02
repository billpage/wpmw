# docs/algorithm

Core algorithm specifications for the WPMW project.

This directory holds the canonical descriptions of the algorithms used in `src/`:
the signed-particle / phase-space-crystal-lattice update rule, the split-Fourier
collision step, and the discrete momentum-jump scheme for sinusoidal and
polynomial potentials.

Specifications here should be precise enough that someone reading them could
re-implement the algorithm in another language. Background discussion and
extended derivations belong in `docs/supplement/` and `docs/analysis/`.
