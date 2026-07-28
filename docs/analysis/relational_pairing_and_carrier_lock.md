# Relational Pairing and the Carrier-Lock Postulate

**Partnership is not stored state. Removing it is correct, costs one postulate, and pays for itself four times over — including a defect in the sea's accounting that the indexed formulation concealed.**

## 0. Status and provenance

This note revises §2.2 of
[`docs/algorithm/phase_alignment_microdynamics_algorithm.md`](../algorithm/phase_alignment_microdynamics_algorithm.md),
which stores on every sea particle a permanent index to a partner and
defines the misalignment as the transported-phase difference of that
stored pair. The revision was prompted by an objection (B. Page, July
2026): a world-particle should not carry an index to another
world-particle, and the natural definition is over *all* ordered pairs
$(i, j)$ of world-particle indices, with $\mu_{ij}$ computed for each.

The objection is sustained, but not for free. §2 shows that partnership
cannot be a carrier of state. §3 shows that the literal all-pairs
definition nonetheless fails as the model currently stands, because
Lemma 1 of the phase-resonance note leaves each pair a private carrier
phase. §4 supplies the missing ingredient as Postulate (S), and §§5–7
collect what (S) buys. §8 records a separate defect found in the course
of the analysis: under permanent partnership the sea is a *consumable*
resource with no replenishment, short by roughly three orders of
magnitude for the parameters of `demo_cosine_well_microdynamics.py`.

**Correction, July 2026.** The first revision of this note and of the
specification wrote the mediating rows as $n \pm q$ while calling the
excess particle's own row $n$, mixing the mesh's midpoint labelling with
the particle labelling in a single expression and thereby placing the
sea legs on odd rows, in violation of Theorem 1. §1 now fixes the two
indices apart. Only index arithmetic was affected: the legs are split by
$2q \cdot dp$ and the excess particle hops by $2q \cdot dp$ in both
readings, no numerical result changes, and the demo never used row
indices.

Nothing in
[`phase_alignment_microdynamics.md`](phase_alignment_microdynamics.md)
is retracted. Lemma 4, Proposition 3, Theorem 4 and its corollaries are
used throughout and are strengthened rather than weakened: Lemma 4's
completeness result is what makes §2 go through, and §8 is a consequence
of Corollaries 4.2 and 4.3 taken together with Lemma 3.

Derivation and implementation were developed jointly with Claude
(Anthropic), July 2026. Companion code:
`src/demo_relational_pairing.py`. All numerical claims in §9 are outputs
of that demo.

## 1. Notation

Ring $x \in [0, L)$, half-grid $dp = \pi\hbar/L$, mode $q$ with
$K_q = 2\pi q / L$ and momentum quantum $2q \cdot dp = \hbar K_q$.

Two row indices are needed and must not be conflated, as §1.1 of the
specification sets out. **$r$** is a *particle's own* row, always even by
Theorem 1 of the phase-resonance note. **$n$** is a *vertex midpoint*
label, of the same parity as $q$, carrying no particle when $q$ is odd.
The dictionary is $n = r + sq$ for direction $s = \pm 1$: the vertex
carries the particle from row $r$ to row $r + 2sq$, and its two mediating
sea legs sit on those same two rows, $n - q$ and $n + q$, straddling the
empty midpoint. The mesh stencil of
[`phase_space_crystal_lattice_algorithm.md`](../algorithm/phase_space_crystal_lattice_algorithm.md)
is written in $n$; the algorithm of §§4–5 of the specification is written
in $r$.

Transported phase and misalignment are as in §2 of the predecessor,

```math
\Phi_j(x,t) \;=\; \theta_j \;+\; \frac{p_j\thinspace(x - x_j) - E_j\thinspace t}{\hbar},
\qquad
\mu_{ij}(x,t) \;=\; \Phi_i(x,t) - \Phi_j(x,t) \pmod{2\pi},
```

the second now written with two free indices rather than with the
partner labels $a$, $b$.

## 2. Partnership cannot be a carrier of state

**Proposition R1 (coboundary).** The two-index family
$\lbrace\mu_{ij}\rbrace$ is the coboundary of the one-index field
$\Phi$. Exactly and at every event,

```math
\mu_{ij} \;=\; -\thinspace\mu_{ji},
\qquad
\mu_{ij} \;+\; \mu_{jk} \;=\; \mu_{ik} \pmod{2\pi}.
```

*Proof.* Substitute the definition; the $\Phi_j$ terms cancel in pairs.
$\square$

**Corollary R1.1.** A stored partnership index carries no information
not already present in $\lbrace\Phi_j\rbrace$. Whatever a partner index
does in the specification, it is not recording relational state; it can
only *select* which pairs are dynamically admissible.

**Corollary R1.2 (Proposition 3 is index-blind).** The winding laws
$\partial_x\mu_{ij} = \Delta p_{ij}/\hbar$ and
$\dot\mu_{ij} = (\Delta p_{ij}/\hbar)(v - \bar v_{ij})$ hold for
arbitrary index pairs, not only for partners, since their proof uses
only the two legs' own data.

Corollary R1.1 settles the ontological half of the objection. Storing an
index is at best a gauge fixing, and §3 shows that this is exactly what
it is.

## 3. The obstruction: a free carrier per pair

**Proposition R2 (residual gauge).** Lemma 1 of the phase-resonance note
defines darkness by the gauge-matching condition
$\theta_a - p x_a/\hbar \equiv \theta_b - p x_b/\hbar$, which constrains
only the *difference* of a pair's phases. The common value is
unconstrained. Hence the gauge group acting on a sea of $B$ dark pairs is
$U(1)^{B}$, not the global $U(1)$ of Lemma 4, and $\mu_{ij}$ is
gauge-invariant when $i$ and $j$ are partners and gauge-variant
otherwise.

**Corollary R2.1 (the naive proposal fails).** Let each pair $k$ carry a
private carrier $\chi_k$ drawn uniformly, with the physical misalignment
$\mu$ common to all pairs at a site. Then

```math
\big\langle \cos\mu_{ij} \big\rangle_{\text{partnered}} = \cos\mu,
\qquad
\big\langle \cos\mu_{ij} \big\rangle_{\text{all pairs}}
= \frac{1}{B^{2}}\thinspace\Big\lvert \sum_k e^{i\chi_k} \Big\rvert^{2} \cos\mu
= O(B^{-1}).
```

The all-pairs average is therefore suppressed as $1/B$ and produces no
vertex bias at all. Verified in §9, Part B.

This is a real obstruction and not a technicality. It means the indexed
formulation is doing work: it silently gauge-fixes $U(1)^{B}$ down to
$U(1)$ by declaring which comparisons are meaningful. The reformulation
must supply that content by other means.

## 4. Postulate (S) and the equivalence theorem

**Postulate (S) — sea carrier lock.** *Sea particles sharing a cell and
a momentum row share a transported phase, up to the pumped
misalignment.*

(S) is the sea-wide strengthening of Lemma 1's per-pair darkness
condition. It is the literal content of the phrase *phase-space
**crystal** lattice*: a crystal is a medium whose scatterers are mutually
phase-locked, whereas a collection of individually dark but mutually
random pairs is a gas. It is also the condition Lemma 5 already needs.
Lemma 5 asserts that every pair at a place carries the same $\mu$ and
derives vertex locality from that mutual alignment; (S) is that
assertion promoted from a property of the pumped state to a property of
the sea.

**Theorem R3 (equivalence).** Under (S) the all-pairs average and the
partnered average coincide identically:

```math
\big\langle \cos\mu_{ij} \big\rangle_{\text{all pairs}}
\;=\; \big\langle \cos\mu_{ij} \big\rangle_{\text{partnered}}
\;=\; \cos\mu(x,t),
```

with the site carrier cancelling exactly. Consequently no result of
the predecessor note that is stated in terms of a pair's $\mu$ changes
value under the reformulation; only its justification does.

*Proof.* Under (S) all $\Phi$ on row $n - q$ at a site equal a common
$\chi$ and all $\Phi$ on row $n + q$ equal $\chi + \mu$, so every ordered
pair drawn across the two rows returns the same $\mu$. $\square$

**Corollary R3.1.** The indexed rule of §5.3 of the specification is the
fully-ordered limit of the relational rule of §5 below; the two agree
whenever Lemma 5 holds, and (S) is exactly the condition for Lemma 5 to
hold.

## 5. The factorisation theorem

Define, for each cell $x_m$ and particle row $r$, the **local coherence
order parameter** and its normalised form

```math
Z_r(x_m, t) \;=\; \sum_{j \in (x_m, \thinspace r)} e^{i\Phi_j(x_m, t)},
\qquad
\hat Z_r \;=\; \frac{Z_r}{N_r},
\qquad
N_r = \#\lbrace j \in (x_m, r)\rbrace,
```

the sum running over sea particles in that cell and row, only even $r$
being populated. Note $\lvert \hat Z_r \rvert \le 1$, with equality
exactly under (S).

**Theorem R4 (factorisation).** For an excess particle at row $r$ in
cell $x_m$ and a vertex in direction $s$, the total affine weight summed
over all admissible mediating pairs is

```math
\sum_{a,\thinspace b} \big[\thinspace w_0 + \kappa\thinspace\cos\mu_{ab}\thinspace\big]
\;=\; w_0\thinspace N_{r+2sq}\thinspace N_{r}
\;+\; \kappa\thinspace\mathrm{Re}\big(Z_{r+2sq}\thinspace\overline{Z_{r}}\big),
```

exactly, with $a$ ranging over the out-row $r + 2sq$ and $b$ over the
in-row $r$ in the cell.

*Proof.* $\cos\mu_{ab} = \mathrm{Re}\thinspace e^{i(\Phi_a - \Phi_b)}$,
and the double sum of a product of a function of $a$ with the conjugate
of a function of $b$ factorises. $\square$

**Corollary R4.1 (the vertex rule).** The firing probability is

```math
P \;=\; w_0 \;+\; \kappa\thinspace\mathrm{Re}\big(\hat Z_{r+2sq}\thinspace\overline{\hat Z_{r}}\big),
```

manifestly independent of the sea depth $B$ and automatically inside
§5.3's clamp $\kappa \le \min(w_0, 1 - w_0)$, since the product of two
normalised order parameters has modulus at most one. Under (S) it
reduces to $w_0 + \kappa\cos\mu$.

**Corollary R4.2 (cost).** The encounter loop of §4 of the specification
costs $O(N_{\mathrm{exc}} \cdot B)$ because each excess particle is
paired against every sea partner in its cell. Under Theorem R4 the sea
enters only through $\lbrace Z_r, N_r \rbrace$, which are assembled in
one pass, so the cost falls to
$O(N_{\mathrm{exc}} + N_{\mathrm{sea}})$. The sea-subsampling workaround
of §4 becomes unnecessary, and the last entry of §11's known gaps is
closed.

**Corollary R4.3 (darkness is derived).** A sea with random phases has
$\lvert \hat Z_r \rvert = O(N_r^{-1/2})$, hence
$\lvert \mathrm{Re}(\hat Z_{r+2sq}\overline{\hat Z_{r}}) \rvert
= O(N^{-1})$: an incoherent sea is dark without a separate postulate,
and the ensemble scaling of open item 4 of the predecessor is read off
$\lvert \hat Z \rvert$ directly.

## 6. Calibration: the stencil rate from $Z$ alone

The pump writes onto each family $s = \pm 1$ the misalignment of Lemma 5,
$\mu^{s}(x, t) = s(K_q x + \phi_q) + \pi/2 - sK_q\bar v^{s} t$, with
contrast $\mu_1 = V_q \tau_p / \hbar$. Under (S) each family contributes
$\mathrm{Re}(\hat Z_{r+2sq}\overline{\hat Z_{r}}) = \cos\mu^{s}
= -\thinspace s \sin(K_q x + \phi_q)$, and with the equal split of the
real pump into two conjugate families the net signed rate is

```math
\frac{1}{\tau_p}\sum_{s = \pm 1} \tfrac{1}{2}\thinspace\mu_1\thinspace s\thinspace\cos\mu^{s}(x)
\;=\; -\thinspace\frac{V_q}{\hbar}\thinspace\sin\!\big(K_q x + \phi_q\big)
\;=\; \Gamma_q(x),
```

the corrected sign of
[`phase_space_crystal_lattice_supplement.md`](../supplement/phase_space_crystal_lattice_supplement.md)
§6.3. The site carrier cancels identically at every step. Verified to
$7.0 \times 10^{-16}$ relative in §9, Part E, at nonzero $\phi_q$.

This is the same result as §7 of the phase-resonance note, obtained
without reference to any pair label. The factor $\gamma/2$ retains its
reading as the equal split of the real pump into two conjugate families.

## 7. The deciding experiment

(S) is a postulate and not a convention, because the two formulations are
physically distinguishable. Two decoherence models separate them.

**Model I — independent per-particle phase noise.** Both estimators have
the same mean, and their variances differ. The partnered estimator is an
average of $N$ weakly correlated terms and its rms residual falls as
$N^{-1/2}$; the relational estimator is a product of two coherent sums
and falls as $N^{-1}$. Measured exponents $-0.496$ and $-1.020$ (§9,
Part F). The relational form is $\sqrt{N}$ quieter at identical mean —
the same variance-reduction pattern as the $G$-freedom of
[`four_rule_microdynamics_equivalence.md`](four_rule_microdynamics_equivalence.md).

**Model II — pair-correlated noise.** Let each pair retain its own
correct $\mu$ while the pairs are mutually randomised, so that Lemma 5
fails but per-pair alignment does not. The partnered estimator returns
the full response $\cos\mu$; the relational estimator returns
$O(N^{-1})$. At $N = 4096$ and $\mu = 0.7$ the two read $+0.7648$ and
$+0.0003$ respectively.

The two forms therefore disagree **in the mean**, not merely in variance,
and Model II identifies the physical content of (S). The relational
answer is the defensible one: a sea with no mutual alignment carries no
macroscopic polarisation grating and must be dark, which is precisely
the premise from which Lemma 5 derives vertex locality.

This inverts the justification offered for the indexed convention in
§2.2 of the specification, which claims that permanent indices make
Lemma 4's invariance testable. By Proposition R1 those invariances are
identities of a coboundary and cannot fail, so testing them proves
nothing. What permanent indices do is make **Lemma 5** untestable, by
hard-wiring the assumption Lemma 5 asserts.

## 8. Correction: the indexed sea is a consumable resource

The following is independent of the pairing question and would apply to
the indexed formulation on its own terms. It is recorded here because
the relational formulation dissolves it.

**Proposition R5 (single use).** Under permanent partnership a pair
mediates at most one vertex in its entire history.

*Proof.* By Corollary 4.2 a struck pair exits with $\Delta p = 0$. By
Corollary 4.3 the transfer available at a vertex is $\lvert\Delta p\rvert$,
which is then zero, so any subsequent admissible encounter with that
pair moves no momentum. $\square$

**Corollary R5.1 (no replenishment).** By Lemma 3 the pump writes phases
and only phases, so it cannot restore $\Delta p \neq 0$ to a spent pair;
and §9 of the specification populates pairs once, at initialisation.
Under permanent partnership the split population is therefore
monotonically non-increasing.

**The size of the deficit.** For the parameters of
`demo_cosine_well_microdynamics.py` — $V_p = 1.5$, $L = 8$, $\hbar = m = 1$,
$N_{\mathrm{exc}} = 5 \times 10^{6}$, $M = 128$ cells, four harmonic
periods — the net firings required per cell over the run are

```math
\Gamma_{\max}\thinspace T_{\mathrm{final}}\thinspace\frac{N_{\mathrm{exc}}}{M}
\;\approx\; 1.5 \times 26.13 \times \frac{5 \times 10^{6}}{128}
\;\approx\; 1.53 \times 10^{6},
```

against $B \approx 2 \times 10^{3}$ pairs per cell: short by a factor of
about $7.7 \times 10^{2}$. The gross traffic at rate $w_0$ also consumes
pairs, so this is a lower bound on the shortfall. The sea is exhausted
within the first fraction of a percent of the run.

**Resolution.** Under the relational definition admissibility is
re-derived from row occupancy at every step, so a struck particle rejoins
the admissible set immediately and no pair is ever spent. Nothing needs
replenishing, and the resource never depletes.

This is likely to be a substantial part of why the steady state of §6 of
the specification is marked open. Under permanent partnership there is
no steady state to find, because the balance it describes has a drain and
no source. The remaining open content of that item — fixing the
calibration constant from a pump/drain balance rather than from L1
exactness — is unaffected and stays open.

## 9. Numerical verification

All claims are exercised by `src/demo_relational_pairing.py` (container
run, July 2026), $\hbar = m = 1$, $L = 8$, $q = 1$, $\phi_q = 0.7$.

- **Part A (Proposition R1).** Antisymmetry to $7.1 \times 10^{-15}$,
  the cocycle identity to $1.4 \times 10^{-14}$, global-gauge and
  re-referencing invariance to $2.1 \times 10^{-14}$. Proposition 3
  holds for arbitrary index pairs, worst residual
  $7.5 \times 10^{-10}$ (finite-difference step, not model error).
- **Part B (Corollary R2.1).** With free per-pair carriers the
  partnered average is exactly $\cos\mu$ at every $B$, while the
  all-pairs average falls from $1.2 \times 10^{-2}$ at $B = 32$ to
  $4.7 \times 10^{-5}$ at $B = 8192$.
- **Part C (Theorem R3).** Under (S) the two averages agree to
  $1.6 \times 10^{-16}$, with
  $\lvert\hat Z\rvert = 1$ on both rows to machine precision.
- **Part D (Theorem R4).** The factorised weight matches the direct
  double sum to $2.1 \times 10^{-16}$ or exactly, at every size tested;
  measured speedup $15\times$ at $N = 128$, $1214\times$ at $N = 2048$
  and $6852\times$ at $N = 8192$.
- **Part E (§6).** $\Gamma_q(x)$ recovered from $Z$ alone to
  $7.0 \times 10^{-16}$ relative, sign and phase offset included.
- **Part F (§7).** Model I exponents $-0.496$ (partnered) and $-1.020$
  (relational); Model II readings $+0.764842$ and $+0.000302$.
- **Part G (§8).** The consumption arithmetic above, reproduced from the
  cosine-well parameters.
- **Part H (§1, and §1.2 of the specification).** Five same-row legs
  spread across a cell agree in transported phase to
  $8.9 \times 10^{-16}$, giving $\lvert\hat Z_r\rvert = 1$; the
  intra-cell pump spread predicted at $\mu_1 K_q \Delta x
  = 1.5 \times 10^{-3}$ rad per pump stays bounded over $1306$ pumps
  ($\lvert\hat Z_r\rvert \ge 0.9996$, first-half mean $0.999878$
  against second-half $0.999855$), so it is oscillatory rather than
  secular.

![Cells and rows](https://raw.githubusercontent.com/billpage/wpmw/output/figures/lattice_cells_and_rows.png)

*(a) exact momentum rows against sweepable spatial cells: solid rows
carry particles at even $`r`$, the dashed midpoint is empty for odd
$`q`$, and the highlighted cell shows a vertex carrying an excess
particle from row $`r`$ to row $`r + 2q`$ under the mediation of those
two rows. (b) the one effect of cell width, the intra-cell pump spread,
over a full cosine-well run.*

![Relational pairing and the carrier-lock condition](https://raw.githubusercontent.com/billpage/wpmw/output/figures/relational_pairing.png)

*(a) the obstruction of §3: a free per-pair carrier destroys the
all-pairs average while leaving the partnered average untouched. (b)
under (S), the stencil rate assembled from $`Z`$ alone against the
analytic $`\Gamma_q(x)`$. (c) Model I of §7: identical mean, one
half-power of $`N`$ in the variance.*

## 10. What changes in the specification

- **§1.1** — new: the row-index convention, separating the particle row
  $r$ from the vertex midpoint $n$ with the dictionary $n = r + sq$.
- **§1.2** — new: what a cell is and what a row is, with the figure of
  §9 below.
- **§2.2** — replaced: no partnership stored; admissibility re-derived
  from cell and row; Postulate (S) stated with its justification and its
  testability.
- **§4** — rewritten around $Z_r$; the subsampling workaround removed.
- **§5.1, §5.3** — the admissibility test becomes a row condition; the
  affine weight becomes Corollary R4.1.
- **§5.4** — the struck partner's exit phase is no longer solved from a
  mate's data; phase is continuous through the vertex for both legs, and
  the vertex transfers coherence between rows.
- **§5.5** — invariant 5 (partnership indices) removed; replaced by the
  coherence-transfer bookkeeping.
- **§10** — new rung 3a, the deciding experiment of §7.
- **§11** — cost entry closed by Corollary R4.2; the consumption defect
  of §8 recorded as resolved rather than open.

## 11. Open items

1. **Species bookkeeping in $Z$.** The order parameter above is written
   without the species factor $\varepsilon_j$. Since a mediating pair's
   legs sit on different rows, $\varepsilon$ is constant within a row and
   factors out as a sign; but the general bilinear form
   $\sum_{\sigma\sigma'} c_{\sigma\sigma'}
   \mathrm{Re}(Z^{\sigma}_{r+2sq}\overline{Z^{\sigma'}_{r}})$ has not
   been fixed from first principles, only specialised to the case where
   it reduces to the above.
2. **Phase continuity at the vertex.** §10 above proposes that both legs
   carry phase continuously through the swap, replacing the
   exit-alignment rule of §5.4 of the specification. Under (S) the two
   agree at the level of the mean generator, but the equivalence has been
   argued and not proved; rung 4 of the validation ladder is the test.
3. **Dynamics of (S).** (S) is imposed as a property of the sea. Whether
   it is preserved by the dynamics — whether vertices and streaming
   degrade $\lvert\hat Z\rvert$ over a run, and at what rate — is not
   known. This is the sharpened form of open item 4 of the predecessor,
   and $\lvert\hat Z_r(x)\rvert$ is the observable to instrument.
4. **Multi-mode potentials.** Each mode $q$ supplies its own row pairing
   $(n - q, n + q)$, so a multi-mode potential requires one order
   parameter per row and a sum over $q$. First-order additivity is
   assumed, as in the predecessor's open item 3.
