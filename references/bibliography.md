# Bibliography

Working bibliography for the WPMW project. Add entries as links to the canonical
source (DOI, arXiv, publisher page). Do **not** commit PDFs to this repository.

## Wigner equation, signed-particle Monte Carlo

- Nedjalkov, M., et al. — signed-particle Wigner Monte Carlo with annihilation.
  *(TODO: add specific references.)*
- Querlioz, D.; Dollfus, P. — *The Wigner Monte Carlo Method for Nanoelectronic
  Devices.* Wiley, 2010. *(TODO: link.)*
- Sellier, J. M. — *Signed-particle formulation of quantum mechanics.*
  *(TODO: add specific references.)*

## Phase-space kinetic theory

- Wiedemann, H. — *Particle Accelerator Physics.* Springer.
  Sections 12.88–12.98 are the basis for the Fokker–Planck derivation extended
  in this work. *(TODO: link.)*
- Kerr, W. C.; Graham, A. J. — "Generalized phase space version of Langevin
  equations and associated Fokker–Planck equations." *(TODO: link.)*

## Stochastic / time-symmetric particle models of quantum mechanics

- McKeon, D. G. C.; Ord, G. N. — "Time Reversal in Stochastic Processes and the
  Dirac Equation." *Phys. Rev. Lett.* **69**, 3 (1992).
- Ord, G. N. — "Random walks, continuum limits and the Schrödinger equation."
  *Phys. Rev. A* **54**, 3772.
- Ord, G. N. — "Entwined Paths, Difference Equations and the Dirac Equation."
  *(TODO: link.)*
- Ord, G. N. — "Can Minkowski Spacetime Resolve Quantum Superposition?" (2017).
- Rajput, B. S. — "Quantum equations from Brownian motion." *Can. J. Phys.*,
  Jan 2011. https://doi.org/10.1139/P10-111

## de Broglie: phase harmony, phase waves, double solution

- de Broglie, L. — *Recherches sur la théorie des quanta.* PhD thesis,
  Paris, 1924; *Ann. de Phys.* (10) **3**, 22 (1925). The **theorem of
  phase harmony** is Chapter 1 §1.1, together with the disk-of-oscillators
  mechanical model showing that a phase wave transports phase but not
  energy. Both are the historical antecedent of the winding law P1 and of
  the misalignment variable of
  `docs/analysis/phase_alignment_microdynamics.md`.
  English translation by A. F. Kracklauer, *On the Theory of Quanta*,
  Fondation Louis de Broglie (2004):
  https://fondationlouisdebroglie.org/LDB-oeuvres/De_Broglie_Kracklauer.pdf
  Original thesis (HAL): https://theses.hal.science/tel-00006807
- Haslett, J. W. — "Phase waves of Louis de Broglie." *Am. J. Phys.*
  **40**, 1315 (1972). Translation of the thesis's first chapter only —
  which is the chapter containing the phase-harmony theorem.
  https://doi.org/10.1119/1.1986827
- de Broglie, L. — "La mécanique ondulatoire et la structure atomique de
  la matière et du rayonnement." *J. Phys. Radium* **8**, 225 (1927). The
  **double solution** paper. *(Note: this is not the source of the
  phase-harmony theorem, contrary to a citation error corrected in
  `docs/analysis/phase_alignment_microdynamics.md` §0.)*
  https://doi.org/10.1051/jphysrad:0192700805022500
- Drezet, A. — "The guidance theorem of de Broglie." (2020). Reviews the
  1923–24 postulate of a local clock attached to every quantum particle
  and the synchronized phase wave accompanying its motion.
  https://arxiv.org/abs/2006.01913
- Shanahan, D. — "The de Broglie wave as evidence of a deeper wave
  structure." (2015). §4 reads phase harmony as the requirement that
  relatively moving observers agree on the phase at each spacetime point.
  https://arxiv.org/abs/1503.02534
- Lochak, G. — "The evolution of the ideas of Louis de Broglie on the
  interpretation of quantum mechanics." *Found. Phys.* **12**, 931 (1982).
  https://doi.org/10.1007/BF01889274

## Diffraction and coupled-wave theory (the vertex)

- Kapitza, P. L.; Dirac, P. A. M. — "The reflection of electrons from
  standing light waves." *Proc. Camb. Phil. Soc.* **29**, 297 (1933). The
  scattering process realised at the contact vertex.
  https://doi.org/10.1017/S0305004100011105
- Batterman, B. W.; Cole, H. — "Dynamical diffraction of X rays by perfect
  crystals." *Rev. Mod. Phys.* **36**, 681 (1964). Two-beam dynamical
  theory: one Hermitian coupling, sine-squared transfer, Pendellösung —
  the classical parent of the vertex rule.
  https://doi.org/10.1103/RevModPhys.36.681

## Direct interparticle action (fields eliminated in favour of worldlines)

- Wheeler, J. A.; Feynman, R. P. — "Interaction with the absorber as the
  mechanism of radiation." *Rev. Mod. Phys.* **17**, 157 (1945).
  Structural precedent for a transported phase with no field degrees of
  freedom.
  https://doi.org/10.1103/RevModPhys.17.157

## Foundational and comparison references

- Wigner, E. — "On the quantum correction for thermodynamic equilibrium."
  *Phys. Rev.* **40**, 749 (1932).
  https://doi.org/10.1103/PhysRev.40.749
- Nelson, E. — "Derivation of the Schrödinger equation from Newtonian
  mechanics." *Phys. Rev.* **150**, 1079 (1966). Source of the "no noise,
  no force" reading of the bare exchange traffic.
  https://doi.org/10.1103/PhysRev.150.1079
- Bohm, D. — "A suggested interpretation of the quantum theory in terms of
  'hidden' variables, I and II." *Phys. Rev.* **85**, 166 and 180 (1952).
  https://doi.org/10.1103/PhysRev.85.166
- Stockburger, J. T.; Grabert, H. — "Exact c-number representation of
  non-Markovian quantum dissipation." *Phys. Rev. Lett.* **88**, 170407
  (2002). Signed/complex-weight trajectory unraveling; the closest
  counterpart to `docs/algorithm/density_matrix_microdynamics_algorithm.md`.
  https://doi.org/10.1103/PhysRevLett.88.170407
- Couder, Y.; Fort, E. — "Single-particle diffraction and interference at a
  macroscopic scale." *Phys. Rev. Lett.* **97**, 154101 (2006).
  https://doi.org/10.1103/PhysRevLett.97.154101
- Bush, J. W. M. — "Pilot-wave hydrodynamics." *Annu. Rev. Fluid Mech.*
  **47**, 269 (2015). Useful as a laboratory analogue of a particle
  interacting through the phase of a medium — but note the disanalogy: the
  Faraday field stores energy and exerts a force, whereas the transported
  phase of the WPMW sea does neither.
  https://doi.org/10.1146/annurev-fluid-010814-014506

## Forward–backward / objective QFT

- Reid, M. D.; Drummond, P. D. — *Objective QFT* programme; parametric
  amplification as canonical measurement basis. *(TODO: add specific
  references.)*
