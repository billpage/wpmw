# Bibliography

Working bibliography for the WPMW project. Add entries as links to the canonical
source (DOI, arXiv, publisher page). Do **not** commit PDFs to this repository.

## Wigner equation, signed-particle Monte Carlo

- Nedjalkov, M.; Kosina, H.; Selberherr, S.; Ringhofer, C.; Ferry, D. K. —
  "Unified particle approach to Wigner–Boltzmann transport in small
  semiconductor devices." *Phys. Rev. B* **70**, 115319 (2004). The
  signed-particle method with same-cell annihilation.
  https://doi.org/10.1103/PhysRevB.70.115319
- Querlioz, D.; Dollfus, P. — *The Wigner Monte Carlo Method for Nanoelectronic
  Devices.* Wiley, 2010. ISBN 978-1-84821-150-6.
- Sellier, J. M.; Nedjalkov, M.; Dimov, I. — "An introduction to applied
  quantum mechanics in the Wigner Monte Carlo formalism." *Phys. Rep.* **577**,
  1–34 (2015). The review of record for the signed-particle method.
  https://doi.org/10.1016/j.physrep.2015.03.001
- Sellier, J. M. — "A signed particle formulation of non-relativistic quantum
  mechanics." *J. Comput. Phys.* **297**, 254 (2015). See the dedicated
  section below; this is the closest published relative of the ontology of
  `docs/analysis/compensated_ontology.md`.
  https://doi.org/10.1016/j.jcp.2015.05.036 (arXiv:1509.06708)
- Shao, S.; Xiong, Y. — "Branching random walk solutions to the Wigner
  equation." *SIAM J. Numer. Anal.* **58**, 2589 (2020). The branching-process
  reading of the signed kernel, with the positive/negative parts as birth
  rates. https://doi.org/10.1137/19M1272408 (arXiv:1907.01897)

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
- Hackebill, A.; Poirier, B. — "On Hydrodynamic Formulations of Quantum
  Mechanics and the Problem of Sparse Ontology." arXiv:2602.21106 (2026).
  Argues that branching under decoherence repeatedly partitions a discrete
  hydrodynamic ensemble until it is too sparse to sustain quantum dynamics, and
  concludes that hydrodynamic completions plausibly need a continuous ontology.
  Answered in `docs/supplement/representation_cost_and_annihilation.md` §9: the
  mechanism requires a dynamics nonlinear in the ensemble density, which the
  Wigner/kinetic formulation does not have.
  https://arxiv.org/abs/2602.21106
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

## The compensated split, and creation/annihilation read as an ontology

Background for [`../docs/analysis/compensated_ontology.md`](../docs/analysis/compensated_ontology.md).
Three claims have to be kept apart when reading this literature: (i) the
*split* of the Wigner generator into a classical-force part and a signed
residual; (ii) the reading of the signed residual as literal creation and
annihilation of particles; and (iii) the further claim that the resulting
particle ensemble is the ontology rather than a sampling scheme. (i) has been
published at least three times as numerics and never as physics; (ii) and (iii)
have been published once, by Sellier, but *without* (i) — his particles are
field-less, so the whole classical force is delivered by pair creation.

### The split, as numerics

- Bordone, P.; Bertoni, A.; Brunetti, R.; Jacoboni, C. — "Infinite barriers
  and classical force in the Wigner-function approach to quantum electron
  transport." *Physica B* **314**, 123–127 (2002). The earliest separation
  found of the classical force from the quantum corrections in a general
  potential profile, with the observation that the integral term reduces to
  the classical-force term for potentials up to quadratic — the transport
  literature's version of Theorem I4.
  https://doi.org/10.1016/S0921-4526(01)01355-2
- Van de Put, M. L.; Sorée, B.; Magnus, W. — "Efficient solution of the
  Wigner–Liouville equation using a spectral decomposition of the force
  field." *J. Comput. Phys.* **350**, 314–325 (2017). The closest published
  relative of the compensated split: the Wigner–Liouville equation
  reformulated around the classical *force* rather than the potential, and
  read explicitly as two processes — classical evolution under the averaged
  driving field, plus a probability-preserving generation-and-annihilation
  term whose non-locality in momentum has only a limited range. That limited
  range is the reach of Definition (R), arrived at independently and treated
  as a numerical convenience rather than as a regulator.
  https://doi.org/10.1016/j.jcp.2017.08.059
- Benam, M.; Nedjalkov, M.; Selberherr, S. — "A Wigner potential decomposition
  in the signed-particle Monte Carlo approach." *Lect. Notes Comput. Sci.*
  **11189**, 263–272 (2019). Splits the Wigner potential so that the
  signed particles experience a force through the classical component;
  motivated by numerical complexity and by self-consistent Wigner–Poisson
  coupling. https://doi.org/10.1007/978-3-030-10692-8_29

### The split plus the branching reading, as foundations

- Limkumnerd, S.; Phanthaphanitkul, P. — "Weighted phase-space paths for exact
  Wigner dynamics." arXiv:2605.05764 (2026). Takes classical Hamiltonian flow
  as the carrier, splits the Wigner generator into the classical Liouville
  part and the Moyal residual, writes the residual as a signed kernel and
  represents it by signed weights or branching events. Independently states
  the harmonic-oscillator null case — nonclassicality carried by the initial
  Wigner function rather than by any correction to the flow — and verifies it
  numerically against a quartic benchmark. The residual is taken in
  differential (Moyal-series) form, so the total-variation obstruction of
  Theorem E7 is not confronted and no reach appears; the positive/negative
  split is Hahn–Jordan, whose minimal-total-variation property is a useful
  reference scale for `Gamma`. The aim is a forward–reverse signed-path
  relation, not an ontology. https://arxiv.org/abs/2605.05764

### Creation and annihilation as a literal particle ontology

- Sellier, J. M. — "A signed particle formulation of non-relativistic quantum
  mechanics." *J. Comput. Phys.* **297**, 254 (2015). Three postulates:
  signed Newtonian particles; a **field-less** classical point-particle that
  during `dt` creates a pair at `p ± p'` with probability `gamma(x) dt`,
  `gamma` being the momentum integral of the positive part of the Wigner
  kernel; and same-cell annihilation. The published statement of the thesis
  that quantum mechanics adds only creation and annihilation. Two points of
  difference from `compensated_ontology.md`: the particles feel no force, so
  `gamma` is large for a harmonic oscillator where the compensated rate is
  exactly zero (§3 there); and admissibility is handled by asserting that the
  uncertainty principle is "embedded in terms of unknown initial conditions"
  rather than as a constraint on the ensemble.
  https://doi.org/10.1016/j.jcp.2015.05.036 (arXiv:1509.06708)
- Sellier, J. M.; Kapanova, K. G. — "A study of the hydrogen atom by means of
  the signed particle formulation." arXiv:1704.06113 (2017).
  https://arxiv.org/abs/1704.06113
- Attar, M.; Selim, B.; Sellier, J. M. — "Efficient approximation of the Wigner
  kernel in phase-space quantum mechanics." arXiv:2606.28269 (2026). The
  signed-particle programme is current. https://arxiv.org/abs/2606.28269

### Deterministic motion punctuated by stochastic creation/annihilation

The same architecture as postulates (S) and (D), in configuration space and
with positive rates only.

- Bell, J. S. — "Beables for quantum field theory." *Phys. Rep.* **137**,
  49–54 (1986).
- Dürr, D.; Goldstein, S.; Tumulka, R.; Zanghì, N. — "Trajectories and
  particle creation and annihilation in quantum field theory." *J. Phys. A*
  **36**, 4143–4149 (2003). https://arxiv.org/abs/quant-ph/0208072
- Dürr, D.; Goldstein, S.; Tumulka, R.; Zanghì, N. — "Bell-type quantum field
  theories." *J. Phys. A* **38**, R1–R43 (2005). The **minimal jump rate**
  construction — the closest existing analogue of `Gamma`, and of the question
  of what fixes it. https://arxiv.org/abs/quant-ph/0407116
- Georgii, H.-O.; Tumulka, R. — "Some jump processes in quantum field theory."
  In *Interacting Stochastic Systems*, Springer (2004), 55–73. Non-explosion
  of the jump process; the counterpart of the supply condition K9.
  https://arxiv.org/abs/math/0312326
- Colin, S.; Struyve, W. — "A Dirac sea pilot-wave model for quantum field
  theory." *J. Phys. A* **40**, 7309–7341 (2007). Creation read as promotion
  out of a filled sea, so particle number is conserved and identity persists.
  The nearest published cousin of Proposition K8's ionisation reading.
  https://doi.org/10.1088/1751-8113/40/26/015
- Deckert, D.-A.; Esfeld, M.; Oldofredi, A. — "A persistent particle ontology
  for QFT in terms of the Dirac sea." *Brit. J. Phil. Sci.* **70**, 747–770
  (2019). The philosophical case for the same move.
  https://doi.org/10.1093/bjps/axx018

### Admissibility: which signed ensembles are quantum states

Postulate (A) of `compensated_ontology.md` is a constraint on ensembles, not
on their motion, and it has its own literature — cited but not used by the
signed-particle formulation.

- Dias, N. C.; Prata, J. N. — "Admissible states in quantum phase space."
  *Ann. Phys.* **313**, 110–146 (2004). The characterisation of the Wigner
  functions of density operators.
  https://doi.org/10.1016/j.aop.2004.04.003
- Narcowich, F. J. — "Conditions for the convolution of two Wigner
  distributions to be itself a Wigner distribution." *J. Math. Phys.* **29**,
  2036 (1988). The `hbar`-positive-definiteness (Kastler–Loupias–Miracle-Sole)
  form of the constraint. https://doi.org/10.1063/1.527860
- Hudson, R. L. — see the entry under interworld couplings below.

### Negative probability as a foundational stance

- Feynman, R. P. — "Negative probability." In *Quantum Implications: Essays in
  Honour of David Bohm*, ed. Hiley & Peat, Routledge (1987), 235–248.
- Mückenheim, W. — "A review of extended probabilities." *Phys. Rep.* **133**,
  337–401 (1986). https://doi.org/10.1016/0370-1573(86)90110-9
- Scully, M. O.; Walther, H.; Schleich, W. — "Feynman's approach to negative
  probability in quantum mechanics." *Phys. Rev. A* **49**, 1562 (1994).
  https://doi.org/10.1103/PhysRevA.49.1562
- Cini, M. — "Quantum mechanics without waves: a generalization of classical
  statistical mechanics." *Ann. Phys.* **273**, 99–113 (1999).
  https://arxiv.org/abs/quant-ph/9807001

## Forward–backward / objective QFT

- Reid, M. D.; Drummond, P. D. — *Objective QFT* programme; parametric
  amplification as canonical measurement basis. *(TODO: add specific
  references.)*

## Interworld couplings and mechanical readings of the four rules

Background for `docs/analysis/interworld_coupling.md`.  The first two are the
existing published attempts to derive quantum behaviour from an interaction
among ensemble members; the third is the closest existing signed-particle
ontology; the last two are the classical kinetic and nonlinear-optical
structures the four-wave-mixing analogy appeals to.

- Hall, M. J. W.; Deckert, D.-A.; Wiseman, H. M. — see the
  many-interacting-worlds entry above.  Note the contrast established in
  `interworld_coupling.md` §1: their interworld potential is present for a
  free particle and is what produces wavepacket spreading, whereas the
  coupling `V(x₁) − V(x₂)` vanishes identically at `V = 0` and free spreading
  is exact classical shearing of `W`.
- Spencer, J. S.; Blunt, N. S.; Foulkes, W. M. C. — "The sign problem and
  population dynamics in the full configuration interaction quantum Monte Carlo
  method." *J. Chem. Phys.* **136**, 054110 (2012).  The plateau: a critical
  walker number below which same-determinant annihilation cannot stabilise the
  sign structure, and above which it can.  The closest existing analogue of a
  density threshold in a signed particle method, and the model for the
  measurement proposed as item N2 of
  `docs/supplement/representation_cost_and_annihilation.md`.
  https://doi.org/10.1063/1.3681396
- Hudson, R. L. — "When is the Wigner quasi-probability density
  non-negative?" *Rep. Math. Phys.* **6**, 249 (1974).  The pure states with
  a non-negative Wigner function are exactly the Gaussians.  Why arbitrary
  Gaussian entanglement costs a world-particle ensemble nothing at all, and
  why non-Gaussianity rather than entanglement is the driver of representation
  cost.
  https://doi.org/10.1016/0034-4877(74)90007-X
- Mari, A.; Eisert, J. — "Positive Wigner functions render classical
  simulation of quantum computation efficient." *Phys. Rev. Lett.* **109**,
  230503 (2012).  Together with Veitch *et al.*, *New J. Phys.* **14**, 113011
  (2012), the statement that negativity in the state or in the measurement is
  the resource separating a phase-space model from efficient classical
  simulation.  The upper bound on how cheap a world-particle ontology can be.
  https://doi.org/10.1103/PhysRevLett.109.230503
- Ferrie, C.; Emerson, J. — "Frame representations of quantum mechanics and
  the necessity of negativity in quasi-probability representations."
  *J. Phys. A* **41**, 352001 (2008).  No quasi-probability representation of
  quantum mechanics is non-negative for both states and measurements.  The
  representation-theoretic sibling of Theorem T5 of the Takabayasi note.
  https://doi.org/10.1088/1751-8113/41/35/352001
- Sellier, J. M. — "A signed particle formulation of non-relativistic quantum
  mechanics." *J. Comput. Phys.* **297**, 254 (2015).  Ensembles of
  field-less Newtonian particles carrying a sign, interacting with an external
  potential only through creation and annihilation events.  The closest
  published relative of the positon/negaton ontology, and the source of the
  same-cell annihilation technique this project does not yet use.
  https://doi.org/10.1016/j.jcp.2015.05.036
- Nedjalkov, M.; Kosina, H.; Selberherr, S.; Ringhofer, C.; Ferry, D. K. —
  "Unified particle approach to Wigner–Boltzmann transport in small
  semiconductor devices." *Phys. Rev. B* **70**, 115319 (2004).  Carries a
  positive Boltzmann collision kernel and the signed Wigner potential term in
  the same equation — the template for what it costs to keep both.
  https://doi.org/10.1103/PhysRevB.70.115319
- Landau, L. D. — "Die kinetische Gleichung für den Fall Coulombscher
  Wechselwirkung." *Phys. Z. Sowjetunion* **10**, 154 (1936).  The classical
  kinetic theory built on a smooth momentum-space pair coupling rather than
  discrete collisions: quadratic in `f`, second order in `p`, positive
  semidefinite.  The structural contrast with the four-action term (linear in
  `W`, odd finite difference, signed) is the subject of §5.
- Deng, L.; Hagley, E. W.; Wen, J.; Trippenbach, M.; Band, Y.; Julienne,
  P. S.; Simsarian, J. E.; Helmerson, K.; Rolston, S. L.; Phillips, W. D. —
  "Four-wave mixing with matter waves." *Nature* **398**, 218 (1999).  Real
  matter-wave four-wave mixing, where the nonlinearity is the genuine
  interatomic scattering length.  The contrast case: the four rules need no
  such nonlinearity.
  https://doi.org/10.1038/18395
