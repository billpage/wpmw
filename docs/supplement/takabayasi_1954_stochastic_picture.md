# Takabayasi's 1954 Stochastic Picture and the Four-Action Model

**A close reading of T. Takabayasi, *The Formulation of Quantum Mechanics in terms of Ensemble in Phase Space*, Prog. Theor. Phys. **11** (1954) 341–373, §3 — "Time development of the distribution" — against this project's derivation ladder, and an assessment of the ontological objection he raises on p. 350.**

---

## 0. Status and provenance

The source is Takabayasi's paper as supplied to the project. All equation
numbers below are his. §3(a) occupies pp. 347–350 and §3(b) pp. 350–353.

The prior project documents involved are
`docs/analysis/four_rule_microdynamics_equivalence.md` (**FR**),
`docs/analysis/sea_dressed_microdynamics.md` (**SD**),
`docs/analysis/phase_resonance_microdynamics.md` (**PR**), and
`docs/supplement/four_action_foundations.md` (**FF**).

Companion code: `src/demo_takabayasi_stochastic_picture.py`, eight parts A–H
plus the figure.
Every numerical claim in §8 is an output of that script.

**Summary of findings.** Takabayasi wrote down, in 1954, the exact kernel that
this project's four-action stencil reproduces — his $J$ *is* our
$`\Gamma_q(x)`$ stencil, and §2 proves the identity. He then read it as a
one-body Markov jump process, found the reading incoherent, and abandoned it as
"picturesque" rather than real. The incoherence is genuine and is sharper than
the negative-probability complaint he actually makes (§3): the QLE stencil is
provably **not** the generator of any one-body Markov process, for two
independent reasons. The four-action model does not repair that reading — it
leaves the class, by making the process two-body and the ensemble two-species.
§4 maps each of his six features onto a project construct. §5 shows the
negativity he objects to is *unavoidable* for any normalised, orthogonal
phase-space kernel (Theorem T5), so the project's achievement is not to remove
it but to relocate it from the transition kernel, where it is meaningless, to
the ensemble, where it is a species label. §6 records two results of his that
the project should adopt, and §7 two errata.

---

## 1. §3(a) in one page

Takabayasi starts from the von Neumann equation (3.1), passes to the
intermediary function

```math
\bar\rho(x,y) \;=\; \rho\bigl(x-\tfrac{y}{2},\; x+\tfrac{y}{2}\bigr)
\;=\; \int f(x,p)\thinspace e^{-ipy/\hbar}\thinspace dp ,
\qquad\text{(3.3)}
```

and Fourier-transforms (3.4) to obtain what he names the **quantum-mechanical
Liouville equation**:

```math
\frac{\partial f}{\partial t} \;+\; \frac{p}{m}\cdot\nabla_x f
\;=\; \Lambda[f] \;=\; \int J(x,\thinspace p-p')\thinspace f(x,p')\thinspace dp' ,
\qquad\text{(3.5)–(3.6)}
```

with the **transition kernel**

```math
J(x,p) \;=\; -\frac{2}{\hbar}\thinspace\frac{1}{(2\pi\hbar)^{3}}
\int V\bigl(x+\tfrac{y}{2}\bigr)\thinspace\sin\frac{p\thinspace y}{\hbar}\thinspace dy
\;=\; -\frac{2^{4}}{\hbar}\thinspace
\mathrm{Im}\bigl[\thinspace\tilde V(2p)\thinspace e^{-2ipx/\hbar}\thinspace\bigr].
\qquad\text{(3.7)–(3.8)}
```

He also gives the expansion of $J$ in derivatives of $V$ (3.10) and the
resulting differential-operator series (3.11)–(3.12), which is the Moyal series
in all but name; he notes that if $V$ is at most quadratic only the first term
survives and the evolution is exactly classical.

His stochastic reading, p. 349, in his own words: *the coordinate $`x`$ of each
particle changes continuously with velocity $`p/m`$, while the value of its
momentum jumps with a "transition probability" $`J`$, where
$`J(x,p)\thinspace dp`$ is the probability per unit time of a momentum jump by
an amount in $`[p,\thinspace p+dp]`$ at the point $`x`$.* He then lists six features (i)–(vi),
and closes:

> *We have thus obtained a stochastic picture of distinctive features for a
> quantum-mechanical change of state. Though this picture cannot be taken as a
> real one, for instance, on account of (i), it may be said to represent quantum
> fluctuations in a picturesque manner.*

---

## 2. His $J$ is our stencil

### 2.1 Single mode

> **Lemma T1.** For $`V(x) = V_q\cos(k_q x + \phi_q)`$, Takabayasi's kernel is a
> pair of opposite-signed point masses,
> ```math
> J(x,p) \;=\; \Gamma_q(x)\thinspace\Bigl[\thinspace\delta\bigl(p + \tfrac{\hbar k_q}{2}\bigr)
> \;-\; \delta\bigl(p - \tfrac{\hbar k_q}{2}\bigr)\thinspace\Bigr],
> \qquad
> \Gamma_q(x) = -\frac{V_q}{\hbar}\sin(k_q x + \phi_q),
> ```
> and therefore his collision operator (3.6) is
> ```math
> \Lambda[f](x,p) \;=\; \Gamma_q(x)\thinspace\Bigl[\thinspace
> f\bigl(x,\thinspace p+\tfrac{\hbar k_q}{2}\bigr) - f\bigl(x,\thinspace p-\tfrac{\hbar k_q}{2}\bigr)\thinspace\Bigr].
> ```

*Proof.* Substitute $`V_q\cos(k_q x + k_q y/2 + \phi_q)`$ into (3.7). The
$`\cos(k_q y/2)`$ part is even in $y$ and annihilated against $`\sin(py/\hbar)`$;
the $`\sin(k_q y/2)`$ part gives
$`\int \sin(ay)\sin(by)\thinspace dy = \pi[\delta(a-b)-\delta(a+b)]`$ with
$`a = k_q/2`$, $`b = p/\hbar`$. Collecting the prefactors yields the stated
result. $\blacksquare$

This is **exactly FR equation (1)**, the target generator of the whole
derivation ladder:

```math
\dot W_n \;=\; \Gamma_q(x)\thinspace\bigl(W_{n+q} - W_{n-q}\bigr),
```

with the identification $`q\thinspace\delta = \hbar k_q/2`$ — the momentum grid
step $`\delta = \pi\hbar/L`$. The dictionary is one line: **David's $\kappa$ is
FR's $`\Gamma_q(x)`$ is Takabayasi's kernel amplitude.**

### 2.2 The lattice is not a discretisation

Something stronger follows, and it is the most useful thing in this note.

> **Proposition T2 (the crystal is Takabayasi's support).** For a potential of
> spatial period $L$, the support of $J(x,\cdot)$ is exactly the momentum
> lattice $`\lbrace\thinspace q\thinspace\pi\hbar/L \;:\; q\in\mathbb{Z}\thinspace\rbrace`$,
> mode $q$ contributing the pair $\pm q$. This is the same lattice on which the
> Wigner function of any period-$`L`$ state is supported, since for
> $`\psi = \sum_n c_n e^{2\pi i n x/L}`$ one has
> $`W(x,p) = \sum_{n,m} c_n c_m^{\ast}\thinspace e^{i(k_n-k_m)x}\thinspace\delta\bigl(p - \pi\hbar(n+m)/L\bigr)`$.

So the phase-space *crystal* lattice is not a numerical discretisation of a
continuum theory and not a modelling postulate. It is the reciprocal lattice of
the spatial period, and both the states and the transitions live on it exactly.
Verified in §8 Parts C and H.

This softens **FF Proposition F4** in a useful way. FF isolates the residual
quantum input as *"the momentum grid step must equal half the photon momentum
of the field mode"*. Proposition T2 says that once the potential is periodic,
the grid **step** is not an extra postulate at all — the only input is de
Broglie's $`p = \hbar k`$, and the lattice follows. What remains genuinely
postulated is that a mode of wavenumber $k$ trades momentum in units of
$\hbar k$, which is exactly the concession FF §7 argues is the right one.

---

## 3. Why the one-body reading fails — and it is worse than he says

Takabayasi's objection is feature (i): $J$ is odd in $p$, hence takes negative
values, and $`\int J(x,p')\thinspace dp' = 0`$, so *"this transition probability
is not normalizable"*. That is correct but understates the problem. Written as
a matrix on momentum cells, his reading asks for the generator

```math
L_{n,\thinspace n+q} = +\Gamma, \qquad
L_{n,\thinspace n-q} = -\Gamma, \qquad
L_{n,n} = 0 .
```

> **Proposition T3 (no one-body Markov reading).** $L$ has column sums zero, so
> it conserves particle number; but its off-diagonal entries have both signs and
> its diagonal is zero. A continuous-time Markov jump generator requires
> non-negative off-diagonals and $`L_{nn} = -\sum_{m\ne n} L_{mn} \le 0`$. Hence
> **no** one-body Markov jump process on momentum, with rates depending on $x$
> alone, has $L$ as its generator.

The two failures are independent and both are visible in Takabayasi's own text.
A negative off-diagonal is his feature (i). A zero diagonal is his (3.14): the
particle's total exit rate vanishes. His picture is self-consistent only
because these two defects cancel — feature (ii), that the rate depends on the
*amount* of the jump alone and not on the momentum before it, forces the loss
term to be $`f(p)\int J\thinspace dp = 0`$, so the master equation reduces to
the gain term alone and reproduces (3.5). It is a formal identity purchased by
letting a particle jump at total rate zero.

There is a third defect he does not remark on. His elementary jump has
magnitude $`\hbar k_q/2`$ — **half a photon of the driving mode**. His own
footnote on p. 349 notices the oddity from the other side: the amplitude
governing a jump from $p_0$ to $p$ is $`\tilde V(2(p-p_0))`$, where
perturbation theory would give $`\lvert\tilde V(p-p_0)\rvert^{2}`$. The factor
of two is exactly the half-photon. FF §8 rules such an event out explicitly:
*"anything transferring $`q\delta`$ to a single particle is half a photon and is
excluded."*

![Takabayasi's kernel, the momentum lattice, and the half-photon problem](https://raw.githubusercontent.com/billpage/wpmw/output/figures/takabayasi_stochastic_picture.png)

*Left: Takabayasi's $`J(x,p)`$ for a localised well — an undamped grating in
$`x`$ whose wavelength $`h/2p`$ shrinks as the jump grows (his feature (iv)),
with amplitude $`4\tilde V(2p)/\hbar`$ that never decays with distance.
Centre: for a periodic potential the same kernel is a comb on the momentum
lattice $`\delta = \pi\hbar/L`$, mode $`q`$ contributing the pair $`\pm q`$
with weights $`\pm\Gamma_q(x)`$ (Proposition T2). Right: the same $`\pm q`$
stencil read three ways — Takabayasi's single particle jumping half a photon
at a signed rate, against the four actions' whole-photon hop across $`2q`$
cells and zero-photon two-body focus, whose eight channel rates are all
non-negative.*

---

## 4. The six features, mapped

| Takabayasi §3(a) | Project construct |
| --- | --- |
| **(i)** $J$ odd in $p$, signed; $`\int J\thinspace dp = 0`$ | The sign is *direction*, not probability. Focus/Defocus and Right-/Left-Hop are antagonistic pairs (FR §1); oddness in $p$ **is** the pairing. SD §4 gives eight individually non-negative channel rates per polarisation block; *"the sign of $`E`$ never enters a rate law."* $`\int J\thinspace dp = 0`$ becomes particle-number conservation, i.e. FF's moment $M_0 = 0$. |
| **(ii)** rate depends on jump amount alone | This is what makes the one-body reading degenerate (§3). The four-action rates depend on *endpoint occupancies* — FF Theorem F3's endpoint locality — which is precisely the departure from (ii) that rescues the process. |
| **(iii)** $V$ induces indeterministic transitions but $J$ is fully determined by $V$ | SD's move: the mean field $`\Gamma_q(x)`$ is *embodied* as a pinned polarised sea population rather than consulted as a field. "What was a field consulted by a rule becomes a population collided with by a particle." |
| **(iv)** $J$ is an undamped sinusoid in $x$, amplitude $`\tilde V(2p)`$, wavelength $h/2p$ | **The grating.** PR §1.4: the coherent polarisation grating at the $q$-mode wavelength, in quadrature with $V$ — his amplitude $`\tilde V(2p)`$ at $`p=\hbar k_q/2`$ is $V_q$, and his wavelength $h/2p$ is $L/q$. "Undamped, irrespective of the form of $V$" is exactly the statement that the response is Bragg diffraction off an extended field, not a local force. |
| **(v)** $`\int p\thinspace J\thinspace dp = -\nabla V`$; individual jumps may far exceed it, so a particle may cross a barrier | Ehrenfest exactness, FF §3.1 — and it is $b=-\tfrac12$ that this fixes. Barrier penetration by finite jumps is the mechanism the crystal-lattice algorithm inherits. |
| **(vi)** even transition moments vanish, odd ones survive; "a quite different type of diffusion" | It is not a diffusion at all. $m_2 = 0$ identically, so the Kramers–Moyal series truncated at second order is *exactly* the classical Liouville equation with no diffusion term whatever. The right comparison class is scattering with a fixed momentum quantum, not Brownian motion — which is what a phase-space *crystal* is for. |

His footnote (\*\*\*) on p. 350 — no frictional force, and no distribution ever
relaxes to an equilibrium — is worth flagging as a match rather than a defect.
SD §2 makes exactly this structural: the sea is *"a non-equilibrium reservoir…
the analogy is a laser medium, not a thermal bath. Detailed balance is
deliberately broken."* A pumped medium has no H-theorem, which is why the
absence of relaxation is expected rather than anomalous.

---

## 5. The assessment, and what "overcoming" it can mean

### 5.1 Negativity is forced

It is worth establishing that Takabayasi's objection cannot be met head-on,
because that changes what a good answer looks like. His §3(b) shows the
finite-time kernel $T$ of (3.18) obeys Chapman–Kolmogorov (3.28), the
"unitarity" condition (3.29), symmetry (3.32) and normalisation (3.33). Read
together, (3.29) and (3.32) say $`T_\tau^{-1} = T_\tau^{\mathsf T}`$: the
phase-space propagator is **orthogonal**. This is not an accident of his
construction — it is the Weyl image of $`\rho \mapsto U\rho U^{\dagger}`$ being
unitary in Hilbert–Schmidt inner product, and the Wigner transform being an
isometry up to a constant.

> **Theorem T5 (negativity is unavoidable).** Let $T$ be a real phase-space
> kernel that is normalised in the sense of (3.33) and orthogonal in the sense
> of (3.29)+(3.32). If $T$ is entrywise non-negative, then $T$ is a permutation
> — that is, a deterministic measure-preserving flow.

*Proof.* Each column $x$ of $T$ satisfies $`\lVert x\rVert_1 = 1`$ by (3.33)
and non-negativity, and $`\lVert x\rVert_2 = 1`$ by orthogonality. For a
non-negative vector, $`\lVert x\rVert_2 \le \lVert x\rVert_1`$ with equality
only if at most one entry is non-zero. Hence every column is a standard basis
vector, and orthogonality makes the assignment a bijection. $\blacksquare$

So any evolution that is genuinely stochastic in phase space and represented by
a normalised orthogonal kernel **must** have negative entries. Takabayasi's
(i) is not a blemish on a particular construction; it is a theorem about the
class. Verified in §8 Part G on a five-dimensional discrete Wigner space.

### 5.2 What the project actually does

The project's answer is therefore not to make the kernel positive. It is to
**stop asking the kernel to be a probability at all**, by factoring the signed
linear map through a genuinely non-negative particle process. Three moves, each
already in the ladder:

1. **Two species.** SD §2: positons and negatons, with $`E = U^{+} - U^{-}`$
   and $`W = E/(\nu\thinspace dx\thinspace dp)`$. A signed distribution is a difference
   of two non-negative counting measures — the semiconductor's electrons and
   holes, not a negative probability. Wigner negativity is a *local negaton
   surplus*, an integer, not a defective measure. This answers the negativity of
   $f$ (his §2 remark that non-positivity *"discloses the physically unreal
   nature of our ps. en."*).
2. **Two directions.** Each signed net rate is the difference of two
   individually non-negative one-way channel rates, with the direction set by
   the local polarisation sign $`\sigma = \mathrm{sign}\thinspace\Gamma_q(x)`$.
   This answers the negativity of $J$ — his feature (i) proper.
3. **The sea supplies the traffic.** Takabayasi's zero total exit rate is the
   statement that the *net* bias is a small residue. SD supplies the missing
   symmetric part as the pinned sea's exchange traffic, at density $B$ per cell.
   At the demo's parameters the eight channel rates total
   $`\sim 4.6\times10^{4}`$ against a net signed rate of $27$: the bias is
   $`6\times10^{-4}`$ of the traffic (§8 Part F). This is the project's *"no
   noise, no force"* principle read backwards — Takabayasi's kernel is the
   antisymmetric residue of an exchange process whose symmetric part he had
   subtracted away, and subtracting it is what left him with a rate that sums to
   zero.

And the half-photon problem is answered by the four actions themselves. The
$`\pm q`$ stencil offsets are not single-particle jumps: PR §1.5 states it
exactly — *"the stencil's half-quantum offsets $`\pm q`$ are the interference
midpoints of a full-quantum transfer."* A hop moves one particle by $2q$ cells,
absorbing one whole photon $`\hbar k_q`$; a focus moves two particles by
$`\pm q`$ cells each, absorbing none. No event transfers half a photon to
anything.

### 5.3 Honest accounting

Three qualifications, so the claim is not overstated.

- Exactness in SD is proved at **level 1**, the pinned sea. The live-sea ledger
  (SD §8, level 2) and the self-consistent sea (level 3) are open. Takabayasi's
  objection is answered *given* a pinned reservoir.
- Uniqueness is of the **generator**, not the process (FF §9.2). Takabayasi's
  $J$ fixes the mean field; the common part $`f^{+}+f^{-}`$ — the traffic — is
  free, and it is exactly that freedom the sea spends. His picture is the
  minimal-traffic member, which is the one member for which the process
  degenerates.
- Theorem T5 stands. The project does not make phase-space evolution
  positivity-preserving; nothing can. It relocates the sign from a transition
  kernel, where a negative entry is meaningless, to a species label, where it is
  a conserved charge.

---

## 6. Two of his results the project should adopt

**(a) A closed form for the whole Kramers–Moyal hierarchy.** His (3.17),

```math
m_{n} \;\equiv\; \int p^{\thinspace n} J(x,p)\thinspace dp \;=\;
\begin{cases}
-\thinspace(-\varepsilon)^{(n-1)/2}\thinspace\dfrac{\partial^{\thinspace n} V}{\partial x^{\thinspace n}}, & n \text{ odd},\\[4pt]
0, & n \text{ even},
\end{cases}
\qquad \varepsilon = \frac{\hbar^{2}}{4},
```

is precisely the list of Moyal coefficients, in one line, for arbitrary $V$.
§8 Part D confirms it against the delta representation of Lemma T1 to
$`1.1\times10^{-16}`$. It is a more compact statement of the expansion than
(3.11)–(3.12), and worth having in the project's notation.

**(b) FF Proposition F4, in moment language.** FF's effective-Planck-constant
result reads, in Takabayasi's moments, as follows: replace the jump size
$`\hbar k_q/2`$ by a free $\delta$ and rescale $\Gamma$ so Ehrenfest still
holds; then $`m_n = -2\Gamma_{\rm eff}\thinspace\delta^{\thinspace n}`$ reproduces
(3.17) *identically* with $`\hbar \to \hbar_{\rm eff} = 2\delta/k_q`$. Verified
to $`5\times10^{-14}`$ across $\delta$ spanning a factor of 100 (§8 Part E).
So "shrink the jump" and "send $`\hbar\to 0`$" are the same operation on his
hierarchy, and his observation that even moments vanish is the reason the
classical limit is reached through the *odd* moments $n\ge 3$ alone.

---

## 7. Errata and flags

- **Feature (iv), the wavelength.** The paper gives the $x$-space wavelength of
  $J$ as $\hbar/2p$. The oscillation is $`e^{-2ipx/\hbar}`$, whose wavelength is
  $`\pi\hbar/p = h/2p`$. The printed value is short by $2\pi$; $h/2p$ is
  correct, and with $`p = \hbar k_q/2`$ it gives $L/q$ as it must.
- **The $2^{4}$ in (3.8) is dimension-dependent.** The substitution
  $`y \to 2(u-x)`$ contributes $2^{d}$ in $d$ dimensions, so the prefactor is
  $`2^{d+1}/\hbar`$: $16/\hbar$ in his 3D setting, $4/\hbar$ in the 1D setting
  the project uses. Confirmed by quadrature to $`2.3\times10^{-11}`$ (§8 Part A).
- **Not an erratum, but a pointer.** His footnote on p. 349 comparing
  $`\tilde V(2(p-p_0))`$ to the perturbative $`\lvert\tilde V(p-p_0)\rvert^{2}`$
  leaves the factor of 2 unexplained. It is the half-photon (§3), and the
  four-action decomposition is its resolution.

---

## 8. Numerical verification

`src/demo_takabayasi_stochastic_picture.py`, units $`\hbar = m = 1`$.

**Part A — (3.7) versus (3.8).** Gaussian well, direct quadrature of (3.7)
against the closed form with the 1D prefactor $`2^{2}/\hbar`$, over sixteen
$(x,p)$ pairs: worst deviation $`2.3\times10^{-11}`$. Feature (iv): the
quadrature matches a pure sinusoid $`A\sin(2px/\hbar)`$ to
$`6.9\times10^{-15}`$, with $`A = 4\tilde V(2p)/\hbar`$ independent of $x$.

**Part B — the collision operator.** $`\int J(x,p-p')f(p')dp'`$ against the
exact Wigner collision term computed through (3.3)–(3.4):
$`9.4\times10^{-12}`$ relative. Against the Moyal series to 25 terms:
$`4.6\times10^{-12}`$.

**Part C — periodic potential (Proposition T2).** Three modes
$`q\in\lbrace 1,2,5\rbrace`$, $L=8$. The four-action stencil
$`\sum_q \Gamma_q(x)[W_{n+q}-W_{n-q}]`$ against the Moyal series to 39 terms, at
seventeen positions: worst relative deviation $`1.8\times10^{-8}`$ (series
truncation). $`\hbar k_q/2 = q\delta`$ exactly for each mode with
$`\delta = \pi\hbar/L`$.

**Part D — transition moments (3.17).** $n = 0\ldots 7$: even moments exactly
zero, odd moments agree with the closed form to $`1.1\times10^{-16}`$. The
second moment — the Fokker–Planck diffusion coefficient — is identically zero.

**Part E — free jump size.** $`\delta/(\hbar k/2) \in \lbrace 0.05, 0.25, 1, 2, 5\rbrace`$:
moments match (3.17) at $`\hbar_{\rm eff} = 2\delta/k`$ to at worst
$`5.0\times10^{-14}`$.

**Part F — the rate ledger (Proposition T3).** The QLE generator $L$ has
column sums $0$, diagonal $0$, and most-negative off-diagonal $-0.910$. Against
$`B = 25\thinspace000`$ sea pairs per cell, the eight SD channel rates total
$`4.55\times10^{4}`$ while the net signed rates total $27.3$ — a bias of
$`6.0\times10^{-4}`$.

**Part G — Theorem T5.** Discrete Wigner representation in dimension $d=5$
(Wootters phase-point operators): $`\lvert\mathrm{Tr}[A_a A_b] - d\thinspace\delta_{ab}\rvert = 1.9\times10^{-14}`$,
$`\lvert\sum_a A_a - d\thinspace\mathrm{I}\rvert = 5.3\times10^{-15}`$. For a random
unitary, the induced kernel $T$ has column sums $1.000000000000$,
$`\lvert T^{\mathsf T}T - \mathrm{I}\rvert = 3.6\times10^{-15}`$, column
$`\ell_2`$ norms $1.000000000000$, and most negative entry $-0.4285$.

**Part H — the state lattice (Proposition T2).** A five-mode periodic state:
$`\bar\rho(x,y)`$ from $`\psi(x-y/2)\psi^{\ast}(x+y/2)`$ against the predicted
delta expansion, $`5.6\times10^{-15}`$; support in units of $`\pi\hbar/L`$ is
$`\lbrace -4,-3,-2,-1,0,1,2,3,4,6\rbrace`$, all integers.

---

## 9. Summary

1. Takabayasi's transition kernel $J$ is, for a single Fourier mode, exactly
   the project's four-action stencil $`\Gamma_q(x)(W_{n+q}-W_{n-q})`$
   (Lemma T1). The two constructions target the identical generator.
2. For a periodic potential, $J$'s support **is** the momentum lattice
   $`\pi\hbar/L`$, which is also the exact Wigner support of any period $`L`$
   state (Proposition T2). The phase-space crystal lattice is not a
   discretisation; it is the reciprocal lattice.
3. His one-body Markov reading fails for three reasons, of which he names one:
   negative off-diagonal rates (his (i)), zero total exit rate (his (3.14)),
   and a half-photon elementary jump (unremarked, though visible in his own
   footnote on $`\tilde V(2(p-p_0))`$).
4. Negativity of the finite-time kernel is **unavoidable** — normalised plus
   orthogonal plus non-negative forces a deterministic flow (Theorem T5).
   Takabayasi's assessment cannot be refuted on its own terms.
5. The project overcomes it by leaving those terms. Sign becomes species
   ($`E = U^{+}-U^{-}`$ ) and channel direction ($`\sigma = \mathrm{sign}\thinspace\Gamma`$);
   zero net rate becomes a small bias on large sea traffic; the half-photon
   becomes an interference midpoint of a whole-photon transfer. All three moves
   were already in the ladder, made for other reasons.
6. Adopt his (3.17) as the compact form of the Kramers–Moyal hierarchy, and note
   that FF Proposition F4 is the statement $`\hbar\to\hbar_{\rm eff}=2\delta/k`$
   in that hierarchy.
7. Two errata: the wavelength in feature (iv) should be $h/2p$, and the
   prefactor in (3.8) is $`2^{d+1}/\hbar`$, not universally $2^4$.
