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

## Quantum hydrodynamics and the moment closure of the Wigner equation

Background for the claim (Cyganski, August 2026) that Bohm's quantum potential
is an artifact of projecting phase space onto its position and momentum
marginals; see `docs/supplement/four_action_foundations.md` §2.

- Takabayasi, T. — "On the Formulation of Quantum Mechanics associated with
  Classical Pictures." *Prog. Theor. Phys.* **8**, 143 (1952).
  https://doi.org/10.1143/ptp/8.2.143
- Takabayasi, T. — "Remarks on the Formulation of Quantum Mechanics with
  Classical Pictures and on Relations between Linear Scalar Fields and
  Hydrodynamical Fields." *Prog. Theor. Phys.* **9**, 187 (1953).
  https://doi.org/10.1143/ptp/9.3.187
- Takabayasi, T. — "The Formulation of Quantum Mechanics in terms of Ensemble
  in Phase Space." *Prog. Theor. Phys.* **11**, 341–373 (1954).
  https://doi.org/10.1143/PTP.11.341 — **the primary citation** for §2: its
  stated aim is to exhibit the correspondence between the phase-space
  formulation and the earlier quantum-potential formulation. The 1952 and 1953
  papers above are the hydrodynamic precursors. §3 of the same paper is the
  earliest statement of the four-action model's target generator as a
  stochastic momentum-jump process, and of the objection to reading it as one;
  see `docs/supplement/takabayasi_1954_stochastic_picture.md`.
- Wyatt, R. E. — *Quantum Dynamics with Trajectories: Introduction to Quantum
  Hydrodynamics.* Springer (2005). Ch. 15 treats the Wigner-moment route to
  the quantum stress tensor.

### Trajectory ensembles that sample rho(x) rather than W(x, p)

The many-interacting-worlds and Holland–Poirier models are the test case for
the moment-closure reading; see `docs/supplement/four_action_foundations.md`
§2.1. Both discretise or parametrise the *streamlines*, so density derivatives
survive as inter-trajectory differences rather than disappearing.

- Hall, M. J. W.; Deckert, D.-A.; Wiseman, H. M. — "Quantum Phenomena Modeled
  by Interactions between Many Classical Worlds." *Phys. Rev. X* **4**, 041013
  (2014). https://doi.org/10.1103/PhysRevX.4.041013 (arXiv:1402.6144). Eq. (23)
  gives the interworld potential; Eq. (20) identifies its gradient with
  Nelson's osmotic momentum; Eq. (38) defines the nonclassical momentum;
  App. A shows convergence to the Bohmian force; App. C gives the exact
  oscillator ground-state recurrence used in Part G of
  `src/demo_four_action_foundations.py`. §I A anticipates the gas-versus-worlds
  objection directly.
- Poirier, B. — "Bohmian Mechanics without Pilot Waves." *Chem. Phys.* **370**,
  4–14 (2010). https://doi.org/10.1016/j.chemphys.2009.12.024
- Schiff, J.; Poirier, B. — "Communication: Quantum Mechanics without
  Wavefunctions." *J. Chem. Phys.* **136**, 031102 (2012).
  https://doi.org/10.1063/1.3680558
- Holland, P. — "Computing the Wavefunction from Trajectories: Particle and
  Wave Pictures in Quantum Mechanics and Their Relation." *Ann. Phys.* **315**,
  505–531 (2005). https://doi.org/10.1016/j.aop.2004.09.008
- Sebens, C. T. — "Quantum Mechanics as Classical Physics." *Philos. Sci.*
  **82**, 266–291 (2015). https://doi.org/10.1086/680190 (arXiv:1403.0014). An
  independent proposal of the same kind, without an explicit model.

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

## Semi-discrete Wigner transport and signed-particle Monte Carlo

Background for `docs/analysis/open_position_space.md`: the standard route to
a discrete momentum space is a finite *coherence length* bounding the ket–bra
separation, which is the same slot the ring circumference occupies in this
project.

- Jacoboni, C.; Bordone, P. — "Wigner transport equation with finite
  coherence length." *J. Comput. Electron.* **13**, 257 (2014). The
  finite-coherence-length transport equation.
  https://doi.org/10.1007/s10825-013-0510-7
- Ellinghaus, P.; Nedjalkov, M.; Selberherr, S. — "Implications of the
  coherence length on the discrete Wigner potential." *SISPAD* (2014).
  https://doi.org/10.1109/SISPAD.2014.6931614
- Sellier, J. M.; Nedjalkov, M.; Dimov, I.; Selberherr, S. — "A benchmark
  study of the Wigner Monte Carlo method." *Monte Carlo Methods Appl.*
  **20**, 43 (2014). The semi-discrete phase space with continuous x and
  discrete k = n Δk, Δk = π/L_C.
  https://doi.org/10.1515/mcma-2013-0018
- Sellier, J. M.; Dimov, I. — "A sensitivity study of the Wigner Monte Carlo
  method." *J. Comput. Appl. Math.* **277**, 87 (2015). Variance-based
  sensitivity of the results to the coherence length.
  https://doi.org/10.1016/j.cam.2014.09.008

## Forward–backward / objective QFT

- Reid, M. D.; Drummond, P. D. — *Objective QFT* programme; parametric
  amplification as canonical measurement basis. *(TODO: add specific
  references.)*
