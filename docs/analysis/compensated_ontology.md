# The compensated ontology: no new force, but new kinematics

**Status.** Analysis note, step 17 of the ladder. Companion demo:
`src/demo_compensated_ontology.py`. Prompted by the question of whether the
picture arrived at in
[`eckart_barrier_compensated.md`](eckart_barrier_compensated.md) §8 and
[`sea_population_equilibrium.md`](sea_population_equilibrium.md) can be stated
as a proposed ontology for quantum mechanics as a whole rather than for
tunnelling in particular.

---

## 0. What this note settles, and what it corrects

The compensated split leaves no quantum force anywhere in the theory. A world
streams on a Newtonian arc under the full classical $`-V'(q)`$, and the only
non-classical element of the dynamics is the ionisation and recombination of
positon–negaton pairs against the sea (K8, S0). The natural summary — *quantum
mechanics implies no new interaction, only the creation and annihilation of
world-particles* — is **half right**, and this note is mostly about the other
half.

**Settles.**

- The dynamical half formalises completely. §1 states four postulates (E),
  (A), (S), (D), and §2 records that the ladder has already proved they
  reproduce the QLE in expectation. Nothing further is needed.
- The other half does not. **Theorem G2**: for a quadratic Hamiltonian the
  demographic channel is *empty* — not small, exactly zero at every reach,
  measured at $`1.2 \times 10^{-15}`$. Yet the harmonic oscillator is
  quantum. So creation and annihilation cannot be the whole of the quantum.
- **Theorem G3**: what carries it instead is postulate (A), the admissibility
  of the ensemble. The two occurrences of $`\hbar`$ — one in the residual
  generator, one in the constraint on states — are independent, and the
  quartic sweep of §4 closes the first continuously while the second does not
  move.
- **Theorem G4**, which is new and was not expected: admissibility is
  preserved by (S) + (D) and **not** by (S) alone. Dropping the residual
  channel from a quartic evolution drives the least eigenvalue of the
  reconstructed $`\rho`$ from $`10^{-8}`$ to $`-0.10`$. The residual channel
  is doing kinematic work as well as dynamical work.
- **Theorem G5**: the emission rate is a property of the regulator and the
  generator is not. Over a sixteenfold range of reach $`\Gamma_{\max}`$ grows
  by a factor of 25 while the generator's reach error falls to
  $`10^{-14}`$. How many worlds exist is therefore not yet a fact of the
  theory, which is a defect an ontology has to answer for. §7 states the
  missing theorem.

**Corrects.**

1. The framing that prompted the note. The correct slogan is **no new
   force**, not *no new physics*. Postulate (A) is non-dynamical,
   non-interactive, and irreducibly quantum, and it is not derivable from
   (S) + (D) — §5 exhibits four Gaussians that the dynamics cannot tell
   apart and that differ by a factor of eight in phase-space area.
2. The implicit ordering in which signed worlds are *produced* by the vertex.
   For quadratic $`V`$ every negaton in the history was present at $`t = 0`$.
   Existence and admissibility of signed worlds is the prior postulate;
   demography is derived, and is empty in an infinite family of cases.
3. Nothing in the numerical record. Every earlier claim used here stands.

**Not new.** The harmonic-oscillator obstruction is already Theorem I4 of
[`interworld_coupling.md`](interworld_coupling.md), is stated in the transport
literature by Bordone *et al.* (2002), and is stated and benchmarked in
Limkumnerd & Phanthaphanitkul (arXiv:2605.05764, May 2026). This note takes it
as a known constraint the compensated reading respects, not as a discovery.
§8 says what does appear to be new.

---

## 1. Four postulates

Let $`\Gamma = \mathbb{R}^{2n}`$ with $`n = 3N`$ the phase space of a world.

**(E) Existence.** The world is a locally finite *signed* counting measure on
$`\Gamma`$,

```math
\mu \;=\; \sum_i s_i \, \delta_{(q_i,\, p_i)}, \qquad s_i = \pm 1 ,
```

whose members are positons, with $`s = +1`$, and negatons, with $`s = -1`$.
There are
no other species and no other attributes.

**(A) Admissibility.** The ensembles that occur are exactly those whose
coarse-grained expectation $`W = \mathbb{E}[\mu]`$ is the Wigner transform of
some density operator $`\rho \succeq 0`$. Equivalently, the chord transform of
$`W`$ is $`\hbar`$-positive-definite in the Kastler–Loupias–Miracle-Sole sense
(Narcowich 1988; Dias & Prata 2004).

**(S) Streaming.** Between events every world-particle obeys

```math
\dot q = p/m, \qquad \dot p = -\nabla V(q)
```

exactly — the *full* classical force. Momentum is continuous along every
worldline, without exception.

**(D) Demography.** Superimposed on (S) is a birth-and-death process. At rate
$`\Gamma(x) = \sum_q |K_{\rm res}(x, \xi_q)|`$ a parent at $`(x, p)`$ ionises a
neutral sea pair co-located with it *on its own momentum row* into a positon at
$`p + \xi_q`$ and a negaton at $`p - \xi_q`$, with $`\xi_q`$ drawn from
$`|K_{\rm res}|/\Gamma`$ and the sign assignment following
$`\mathrm{sgn}\,K_{\rm res}`$. Recombination is the reverse. The realisation of
each event is chosen by supply, and the absorptive fraction is $`f = 1/2`$.

Postulate (D) has a presupposition that (E) does not supply: a reservoir of
neutral bound pairs, dark in every observable, at every cell and momentum row.
That is the sea of
[`species_sectors_and_annihilation.md`](species_sectors_and_annihilation.md)
Theorem D2.1 — the $`c = 0`$ member of the dark family, not the crystal shift
$`c = 2`$ — with the population fixed by the attractor S9 rather than by a
parameter.

---

## 2. G1: the postulates give the QLE, and this is not new work

**Theorem G1 (assembly).** Under (E), (S), (D), $`\mathbb{E}[\mu_t]`$ solves
the quantum Liouville equation exactly.

*Proof.* This is a restatement of results already in the ladder, and is
recorded here only so that §3 onwards has something definite to be about.
C1 gives the factorisation of the potential substep with no Trotter error
within it; C2 identifies the residual as the odd part of the cubic Taylor
remainder; C3 gives it zero zeroth and first moments on a bounded reach, hence
a bounded signed jump measure conserving worlds and carrying no net momentum;
K4 shows the deterministic step conserves the classical outcome functional
exactly, so the residual channel carries the entire quantum correction; K8
forces the emitted pair to be an ionised sea pair co-located with the parent;
and S7 gives $`dN = 2(1-2f)\,n_{\rm ev}`$, $`dS = (2f-1)\,n_{\rm ev}`$, so
$`f = 1/2`$ closes both ledgers. $`\square`$

Postulate (A) plays no part in G1, and that is the whole difficulty of the
next three sections.

---

## 3. G2: the demographic channel is empty for quadratic Hamiltonians

**Theorem G2.** If $`V`$ is quadratic then $`D_{\rm res}(x,y) \equiv 0`$ for
all $`x`$ and all $`y`$, hence $`\Gamma \equiv 0`$ at every reach: under (S) +
(D) a harmonic oscillator, an inverted oscillator, a free particle and a
uniform field take **no events at all**, ever.

*Proof.* For $`V = \tfrac12 m\omega^2 x^2`$,
$`D = V(x+y) - V(x-y) = 2 m\omega^2 x y = 2 y V'(x)`$ identically, so the
compensation removes all of it. This is Theorem I4 read in the compensated
variables. $`\square`$

Part A of the demo, over $`y_{\max} = \pi/2,\ \pi,\ 2\pi`$ at 25 positions,
$`\Gamma_{\max} = \max_x \sum_q |K_{\rm res}|`$:

| potential | $`y_{\max}=\pi/2`$ | $`y_{\max}=\pi`$ | $`y_{\max}=2\pi`$ |
|---|---|---|---|
| harmonic $`x^2/2`$ | 3.83e-16 | 5.12e-16 | 1.24e-15 |
| inverted $`-x^2/2`$ | 3.83e-16 | 5.12e-16 | 1.24e-15 |
| quartic $`\lambda = 0.02`$ | 4.75 | 38.0 | 304 |
| Eckart $`\mathrm{sech}^2`$ | 5.52 | 14.0 | 28.7 |

The same table for the **field-less** rate — the kernel with the classical
part left in, which is the published signed-particle formulation's
$`\gamma(x)`$:

| potential | $`y_{\max}=\pi/2`$ | $`y_{\max}=\pi`$ | $`y_{\max}=2\pi`$ |
|---|---|---|---|
| harmonic $`x^2/2`$ | 28.5 | 56.9 | 114 |
| inverted $`-x^2/2`$ | 28.5 | 56.9 | 114 |
| quartic $`\lambda = 0.02`$ | 53.7 | 136 | 500 |
| Eckart $`\mathrm{sech}^2`$ | 3.18 | 2.95 | 1.24 |

Two things to read off, and one caveat.

**The two readings disagree observably about the census of an oscillator.**
Under the field-less reading every world in a harmonic well is emitting pairs
at a large rate, because the restoring force is *made of* pair emission there.
Under the compensated reading nothing happens at all. Both reproduce the same
$`W`$. This is the sharpest available argument for the compensated split as
ontology rather than as numerics: it is the difference between a world that
has a Newtonian history and a world that does not.

**Compensation is not universally cheaper in raw event count.** On the bounded
Eckart potential the field-less rate *falls* with reach while the compensated
rate rises, because $`D`$ is bounded whereas $`D_{\rm res} \to -2yV'(x)`$ grows
linearly. That is consistent with Theorem K3, whose budget ratio is weighted by
the potential's own spectrum; it is the unweighted count that reverses. The
ontological argument above does not depend on the count.

*Caveat.* Per the erratum in
[`../algorithm/compensated_liouville_algorithm.md`](../algorithm/compensated_liouville_algorithm.md)
§4.4, absolute $`\Gamma`$ values under a hard horizon are functions of the rung
count (64 here) and are not absolute numbers. Only the exact zero, and the
scaling in §7, are rung-independent.

---

## 4. G3: the two occurrences of $`\hbar`$ are independent

If (D) is empty for quadratic $`V`$, the quantum content of an oscillator has
to be somewhere else. It is in the state.

Part B is the null test. For $`\psi = (|0\rangle + |2\rangle)/\sqrt{2}`$ the
Wigner function has $`\min W = -0.1656`$ against the bound $`2/h = 0.3183`$ and
a negativity volume $`\int|W| - 1 = 0.4359`$; and the classical phase-space
rotation carries it exactly. Against the exact Wigner function at
$`t = \pi`$, where the rotation maps grid points to grid points and no
interpolation is needed, the discrepancy is $`1.6 \times 10^{-6}`$; at
$`t = \pi/4`$ and $`\pi/2`$ it is $`9.6 \times 10^{-4}`$, which is the
bilinear interpolation floor of the rotated grid rather than a physical error.
The negativity volume is constant to five figures throughout.

So a *signed* ensemble is carried by a purely classical flow, and its
signedness is neither produced nor consumed by the dynamics.

**Theorem G3.** The residual generator and the admissibility constraint are
independent functions of $`\hbar`$: the first can be taken continuously to
zero while the second is held fixed.

Part C sweeps the quartic coupling in $`V = x^2/2 + \lambda x^4`$ at fixed
state, comparing
$`\chi_Q = \|Q_\hbar W\|_2 / (\|L_{\rm cl}W\|_2 + \|Q_\hbar W\|_2)`$ against
the state's negativity volume:

| $`\lambda`$ | $`\chi_Q`$ | negativity |
|---|---|---|
| 0 | 0 (exactly) | 0.43567 |
| 1e-4 | 6.86e-4 | 0.43567 |
| 1e-3 | 6.81e-3 | 0.43567 |
| 3e-3 | 2.01e-2 | 0.43567 |
| 1e-2 | 6.28e-2 | 0.43567 |
| 2e-2 | 1.15e-1 | 0.43567 |
| 5e-2 | 2.29e-1 | 0.43567 |
| 1e-1 | 3.31e-1 | 0.43567 |

![Two independent occurrences of hbar](https://raw.githubusercontent.com/billpage/wpmw/output/figures/compensated_ontology_two_doors.png)

The dynamical door closes smoothly and shuts completely at $`\lambda = 0`$.
The kinematic door does not move. Any claim that quantum mechanics *is*
creation and annihilation has to account for the left-hand end of that graph,
where there is no creation and no annihilation and the physics is still
quantum.

---

## 5. Postulate (A) is independent, and the violated inequality is the Wigner bound

**Proposition G3.1.** (A) is not derivable from (S) + (D).

Part D exhibits four isotropic Gaussians in the harmonic well. All are exactly
stationary under (S) — the measured $`\|L_{\rm cl}W\|_\infty`$ is at the
spectral floor — and (D) is empty by G2, so no dynamical statement
distinguishes them:

| $`\sigma_q\sigma_p/\hbar`$ | peak $`W`$ / $`(2/h)`$ | purity | $`\|L_{\rm cl}W\|_\infty`$ | admissible |
|---|---|---|---|---|
| 0.125 | 4.000 | 4.000 | 1.7e-15 | no |
| 0.250 | 2.000 | 2.000 | 1.6e-15 | no |
| 0.500 | 1.000 | 1.000 | 1.6e-15 | yes |
| 1.000 | 0.500 | 0.500 | 8.8e-9 | yes |

The inequality that fails for the first two is exactly the Wigner bound
$`|W| \le 2/h`$, equivalently $`\mathrm{Tr}\,\rho^2 \le 1`$ — the same bound
that
[`permanent_pairing_density_matrix.md`](permanent_pairing_density_matrix.md)
identifies as the storage-feasibility condition, and that
[`fourd_microdynamics.md`](fourd_microdynamics.md) identifies with state
purity. In the present reading it acquires a second job: it is the constraint
that says which censuses of positons and negatons are physically possible. It
is a statement about ensembles, not about their motion, and no amount of
demography will produce it.

This is where the ordering of the postulates matters. If (D) were the whole
story one could hope to derive the admissible set from the dynamics. Since (D)
is empty for an infinite family of Hamiltonians whose states are nonetheless
signed, the existence and admissible arrangement of signed worlds has to come
first, and demography second.

---

## 6. G4: streaming alone leaves the admissible set

This was not anticipated, and it strengthens the case for (D) considerably.

**Theorem G4.** (A) is propagated by (S) + (D), because (S) + (D) is the QLE
and the QLE is the Weyl transform of unitary evolution. It is **not**
propagated by (S) alone whenever the residual is non-zero.

Part E, quartic $`\lambda = 0.05`$, initial state
$`(|0\rangle + |2\rangle)/\sqrt{2}`$, least eigenvalue of the $`\rho`$
reconstructed from $`W`$ by the inverse Wigner transform:

| $`t`$ | (S) + (D), exact | (S) alone |
|---|---|---|
| 0 | -6.1e-7 | -6.1e-7 |
| $`\pi/4`$ | -1.0e-8 | -4.93e-2 |
| $`\pi/2`$ | -8.7e-8 | -4.97e-2 |
| $`\pi`$ | -2.9e-8 | -1.01e-1 |

![Admissibility under (S)+(D) and under (S) alone](https://raw.githubusercontent.com/billpage/wpmw/output/figures/compensated_ontology_admissibility.png)

The value at $`t = 0`$ is the reconstruction floor of the $`256 \times 256`$
grid and both columns share it. Under exact evolution the eigenvalue stays at
or below that floor for the whole run; under classical carrier transport alone
it is six orders larger by $`t = \pi/4`$ and still growing at $`t = \pi`$.

The consequence for the ontology is worth stating plainly. Postulate (D) is
not only what makes the expectation match the QLE. It is also what keeps the
ensemble inside the set of ensembles postulate (A) allows. A classical
ensemble of worlds streaming Newtonially in an anharmonic potential *ceases to
be a quantum state* — the pair traffic is what holds it in the admissible set.
That is a considerably better motivation for (D) than "the numbers come out
right", and it is the closest thing this note has to an argument that the sea
is not optional.

---

## 7. G5: the census is a regulator and the generator is not

**Theorem G5.** $`\Gamma`$ grows without bound as the coherence reach grows,
while the generator it represents converges.

Part F, Eckart $`V_0 = a = 1`$, probe at $`x = 0.6`$. The generator is tested
by acting with $`L_{\rm res}`$ on a Gaussian $`W`$ of width $`\sigma_p = 1`$
and comparing the reach-truncated action with the untruncated one — the
appropriate test, since $`M_{\rm res} \to -i s V'(x)`$ grows linearly and there
is no naive $`y_{\max} \to \infty`$ limit of the operator itself:

| $`y_{\max}/a`$ (units of $`\pi`$) | $`\Gamma_{\max}`$ | $`\int \Gamma\,dx`$ | rel. error of $`L_{\rm res}W`$ |
|---|---|---|---|
| 0.5 | 5.60 | 12.3 | 2.61e-2 |
| 1 | 14.2 | 40.0 | 1.11e-8 |
| 2 | 29.0 | 77.0 | 1.02e-14 |
| 4 | 62.9 | 165 | 1.02e-14 |
| 8 | 143 | 372 | 1.02e-14 |

![Census diverges, generator converges](https://raw.githubusercontent.com/billpage/wpmw/output/figures/compensated_ontology_regulator.png)

The rate grows roughly linearly in the reach — consistent with the momentum
churn of the CLA erratum — while the operator is converged to machine
precision by $`y_{\max} = 2\pi a`$.

For an *algorithm* this is unremarkable and even welcome. For an *ontology* it
is a defect: how many worlds there are, and how often they are born, ought to
be a fact rather than a bookkeeping choice. There are two ways out and the
project has not chosen between them.

- **The renormalisation reading.** $`\Gamma`$ is a bare quantity, divergent
  with the cutoff, and only the compensated observable is renormalised. The
  analogy is exact enough to be useful: the uncompensated kernel is the bare
  theory with the tree-level force buried in the vertex; compensation is the
  subtraction leaving a purely radiative residual; $`y_{\max}`$ is the cutoff;
  the soft horizon is the proper regulator; and Proposition K8 is the
  Dirac-sea rather than the Feynman reading of pair production. On this
  reading world-count joins Bohmian position as unobservable in principle.
- **The physical-horizon reading.** $`L_c`$ is a real environmental scale, as
  [`open_position_space.md`](open_position_space.md) Definition (H) allows, in
  which case $`\Gamma`$ is physical and there is nothing to explain — but then
  the predictions of quantum mechanics have to be exactly $`L_c`$-independent
  while the ontology is not, which is again the structure of a gauge choice.

Either way the same theorem is missing, and it is load-bearing:

> **Open item G-SP1 (split gauge invariance).** Any two admissible compensated
> splits — differing in reach, in horizon profile, or in the allocation
> constant $`c`$ of the hop channel — give the same $`\mathbb{E}[\mu]`$ on all
> observables. Not proved. Theorem K3's saturation and the numbers above are
> evidence, not a proof, and Theorem K3 is itself qualified by the sawtooth
> pathology at $`y_{\max} \gtrsim a`$.

---

## 8. Scope, and what appears to be new

**Generalises without change.** The $`N`$-body case is verbatim: the symbol is
$`M(x, s) = (i/\hbar)[V(x+y) - V(x-y)]`$ with $`x, y \in \mathbb{R}^{3N}`$, so
entanglement costs the construction nothing. The world-particle is a whole
configuration, which is the many-interacting-worlds commitment the project has
already made — with the difference recorded in
[`interworld_coupling.md`](interworld_coupling.md) §1 that the coupling here
vanishes identically for a free particle.

**Does not generalise.** Spin and every other finite-dimensional system. The
discrete Wigner kernels on the sphere and on the Weyl–Heisenberg group are not
differences of a potential, so there is no ionisation reading of them at all.
"All of quantum mechanics" is therefore not yet warranted even after §3–§6 are
granted.

**Untouched.** Measurement, with the same status as in Bohmian mechanics and
in Bell-type quantum field theories: an effective decoherence story is
required and is not supplied here.

**What appears to be new**, stated narrowly, after the survey in
`references/bibliography.md`:

1. The compensated split as an *ontological requirement* rather than a
   numerical one — the claim that without it no world has a Newtonian history.
   Bordone *et al.* (2002), Van de Put *et al.* (2017) and Benam *et al.*
   (2019) split and do not interpret; Sellier (2015) interprets and does not
   split; Limkumnerd & Phanthaphanitkul (2026) split and interpret, but toward
   a fluctuation relation. §3's two tables are the discriminating measurement.
2. Proposition K8: ionisation from a co-located sea pair on the parent's own
   momentum row, forced by body-momentum conservation. The nearest published
   cousin is the Dirac-sea pilot-wave model of Colin & Struyve (2007), in a
   different framework.
3. Theorem G4: the residual channel is what keeps the ensemble admissible.
   Not found in any of the above.
4. The reach-dependence of $`\Gamma`$ taken seriously as an ontological
   problem, and the absorptive $`f = 1/2`$ attractor of S9. The published
   annihilation steps are grid-based variance reduction, not supply-limited
   mechanisms.

---

## 9. Summary

| | |
|---|---|
| **G1** | (E)+(S)+(D) reproduce the QLE in expectation — assembled from C1–C3, K4, K8, S7 |
| **G2** | $`\Gamma \equiv 0`$ for quadratic $`V`$ at every reach ($`1.2\times10^{-15}`$), while the field-less rate is $`\mathcal{O}(10^2)`$ |
| **G3** | the residual generator and the admissibility constraint are independent functions of $`\hbar`$ |
| **G3.1** | (A) is not derivable from (S)+(D); the violated inequality is $`\|W\| \le 2/h`$ |
| **G4** | (A) is propagated by (S)+(D) and not by (S) alone: $`10^{-8}`$ against $`-0.10`$ |
| **G5** | $`\Gamma`$ diverges with the reach while the generator converges to $`10^{-14}`$ |

The ontology in one line: **every world has a Newtonian history; what quantum
mechanics adds is a census and a grammar.** The census is postulate (D) — how
many worlds, of which sign, born and dying in pairs against the sea. The
grammar is postulate (A) — which censuses are possible at all. The compensated
split eliminates quantum *force*; it does not eliminate quantum *kinematics*,
and §5 is the reason no reformulation of the dynamics could.

---

## 10. Open items

- **G-SP1.** Split gauge invariance. §7. Load-bearing: without it the world
  count has observable content that depends on a bookkeeping choice.
- **G-SP2.** Whether (A) can be given a *particle-level* statement. At present
  it is a condition on $`\mathbb{E}[\mu]`$, which is awkward for an ontology
  of individual worlds: it constrains the ensemble as a whole rather than any
  world in it. Compare the Wigner bound's per-cell form in
  [`species_sectors_and_annihilation.md`](species_sectors_and_annihilation.md)
  §2, which is local and may be the right handle.
- **G-SP3.** Quantify G4. Is the rate at which (S)-alone transport leaves the
  admissible set related to $`\chi_Q`$, or to $`\Gamma`$, or to neither? A
  scaling law here would turn G4 from a demonstration into a theorem about how
  much demography a given anharmonicity requires.
- **G-SP4.** Finite-dimensional systems. §8 records that spin has no
  ionisation reading. Whether *any* signed-world ontology covers it is open
  and is not addressed by anything in the ladder.
- **G-SP5.** The soft-horizon variant of the fold prescription, inherited from
  K3, is a prerequisite for a clean statement of G5: under a hard cutoff the
  momentum churn diverges linearly and the leading Moyal coefficient is not
  recovered analytically.
