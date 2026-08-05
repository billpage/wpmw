# The coherence ladder

**Status.** Analysis note, step 8 of the ladder. Companion demo:
`src/demo_coherence_ladder.py`.

## 0. What this note inherits, retracts and corrects

**Inherits.** From the phase-resonance note: the K3/K4 vertex kinematics,
Proposition 2, the contact vertex $h = g_0 + g_1 C e^{i\delta}$ of §13.2
with its interference-licensing subtlety, and detailed balance from
$U^{\dagger}$. From the phase-alignment note: Proposition 3 (winding
rates), Theorem 4, and the three pair states. From
[permanent_pairing_density_matrix.md](permanent_pairing_density_matrix.md):
the pair-as-element ontology and the four-term commutator.

**Corrects** five statements, three of them in the revised algorithm
specification pushed in July 2026 and one made in discussion since:

1. The proposal that ladder climbing ($k \to k \pm 1$ for $k \ge 1$)
   requires extending the write channel to split-pair partners with
   mate-anchored kinematics — wrong mechanism. Climbing is **leg-local**:
   the struck leg's *own* pump sideband supplies the co-moving pattern
   (§2), and the mate is never queried.
2. The specification's "sea–sea channel ... **[choice]**, omitted" — the
   sea-striker channel is not optional. It is the engine of every
   interior ladder edge (§4); only the population boundary is
   excess-struck.
3. The specification's erase amplitude ("unity for a split pair") — the
   linear-in-stored-contrast response of an isolated pair is unlicensed
   (§5); the licensed erase bias is the bilinear pump-stored interference,
   proportional to $\mu_1$, and all transfers freeze at $V = 0$.
4. The specification's §7 claim that populations and coherences occupy
   disjoint row parities for odd $q$ — true only of odd rungs; even
   rungs ($k = 2, 4, \ldots$) sit on even midpoint rows (§7).
5. The scope of Corollary 4.2: "the pair exits aligned" is the exact
   (K3) channel only. The write/K4 exit is a winding pair — as the
   resonance note's K4 always said.

## 1. Tutorial: the ladder, and two analogies

Index every momentum-basis element $\rho(P, P')$ by its midpoint
$\bar p = (P + P')/2$ and its **rung** $k = (P - P')/(2q\thinspace dp)$
for mode $q$. Rung 0 is the diagonal: populations, carried by excess
particles. Rung $\lvert k\rvert \ge 1$ is carried by split pairs with
legs $2kq\thinspace dp$ apart; the pair's misalignment $\mu$ is the
element's argument, winding at the beat rate of Proposition 3 — the
exact free evolution of the stored element, not a defect.

The commutator with one cosine mode moves amplitude by **one diagonal
step at a time**: each of its four terms hops one leg by
$\pm 2q\thinspace dp$, changing $(k, \bar p)$ by $(\pm 1, \pm q\thinspace dp)$.
Think of draughts on the $\rho$ board: the potential is the only player,
every legal move is one diagonal step, and the rank-0 pieces are the
excess particles. Or in ledger terms: every population move writes or
redeems exactly one IOU one rung up, at the transition midpoint; write
and erase are the same move in opposite directions, and Hermiticity of
$V$ is the rule that every arrow has its reverse at matched strength.

![The coherence ladder](https://raw.githubusercontent.com/billpage/wpmw/output/figures/coherence_ladder.png)

## 2. The channel table (Lemma C1)

A vertex is the Theorem-4 swap between a striker at $p_{\mathrm{in}}$
and a struck leg at $p_{\mathrm{out}} = P_{\mathrm{struck}}$; the struck
leg exits at $p_{\mathrm{in}}$. The exchange accumulates over $\tau_e$
only if a pattern spanning exactly the exchanged rows co-moves with the
exchange midpoint (§13.2's Rabi-versus-detuning, Proposition 3's winding
mismatch integrating to a sinc). At first order in the pump the patterns
available on a pair $(\alpha$ at $P_a$, $\beta$ at $P_b)$ with $\alpha$
struck are three, hence exactly three admissible $p_{\mathrm{in}}$
classes:

| $p_{\mathrm{in}}$ | pattern | contrast | pair after | name |
| --- | --- | --- | --- | --- |
| $P_b$ | the pair's own winding | stored | aligned at $p_{\mathrm{in}}$: $k \to 0$ | exact (K3) |
| $P_a \pm 2q\thinspace dp$ | struck leg's own sideband | $\mu_1$ | $k \to k \pm 1$ | ladder |
| $P_b \pm 2q\thinspace dp$ | mate's sideband | $\mu_1$ | $k \to \pm 1$ | compound |

Part A of the demo verifies the selection: the own-winding sinc is $1$
at $p_{\mathrm{in}} = P_b$ and $0.02$ everywhere else at the test
window, while the sideband patterns span the exchanged rows identically
and are co-moving by construction. For an aligned pair the exact and
ladder channels coincide with the write channel of the specification.

## 3. The ladder theorem (Theorem C2)

Take the four **ladder** channels — struck leg ket or bra, direction
$\pm$ — with: strikers whose density is state-independent (§4); phase
continuity through the vertex, so the written element carries the source
element's argument; the refractive factor $-i\thinspace e^{\pm i\phi_q}$
of Lemma 2 on ket strikes and its conjugate on bra strikes (the
negaton's opposite winding **is** the complex conjugation of the bra
side); and one shared coupling constant calibrated once,
$c = V_q/2\hbar$. Then the expected elementwise flux equals the
commutator $-\tfrac{i}{\hbar}[V, \rho]$ **exactly, on every element of
every state**.

Part B of the demo: on a random mixed state over 14 rows the assembled
flux matches the commutator to $1.7 \times 10^{-16}$ relative, is
Hermitian to $1.4 \times 10^{-17}$, preserves the trace to
$2.4 \times 10^{-18}$, freezes identically at $V = 0$ (every channel
amplitude carries $V_q$), and produces zero diagonal flux on a
coherence-free state — Lemma 3 re-derived at the ladder level.
Populations move only as rung-1 elements feed them, which is the
midpoint-mediated structure of the crystal lattice, now emerging from
leg-local rules rather than being counted in. The four channels are the
natural candidates for the four rules of
[four_rule_microdynamics_equivalence.md](four_rule_microdynamics_equivalence.md);
the algebraic identification remains to be written out.

## 4. Who strikes whom

The flux out of element $(P_a, P_b)$ must be linear in that element
alone. An excess-struck vertex is bilinear — it needs an excess particle
at $p_{\mathrm{in}}$ — so excess strikers can only carry fluxes whose
source includes a population factor: the rung $0 \leftrightarrow 1$
boundary, where the moved worldline is the excess itself and the rate is
proportional to the population that moves. Every interior edge
($k \to k \pm 1$, $k \ge 1$) is struck by the **background sea**, whose
leg density is uniform and state-independent: the flux is then linear in
the stored-pair count, as the commutator demands. Two consequences:

- **Locality strengthens.** Interior edges use the struck leg's own
  sideband; the mate is never read. The locality lemma's burden shrinks
  to the exact channel alone.
- **Striker back-reaction [open].** The displaced background leg leaves
  its own pair split by one quantum. This churn must be neutral in
  expectation for the one-body sector; detailed balance and family
  pairing make it plausible, and the integrator must instrument it.

## 5. The exact channel and licensing

The exact channel collapses any rung to alignment, transferring the
pair's whole splitting at once — a process the first-order commutator
does not contain, except as the $k = 1$ erase. The resolution is §13.2's
licensing subtlety, run in reverse:

- the **bare** part ($g_0$, contrast-free) fires symmetrically with its
  $U^{\dagger}$ reverse: expected net zero, pure sampling churn;
- the **linear-in-stored-contrast** bias of an isolated pair is
  unlicensed: absorption and bare exchange leave distinguishable
  records, so the response is quadratic, not linear — no spurious
  first-order collapse from $k \ge 2$;
- at $k = 1$ the exact and ladder channels reach the **same final
  state** (pair aligned, excess hopped), so the pump path and the stored
  path interfere with no reservoir needed: the cross term
  $\propto \mu_1 \cos(\mu - \Lambda)$ is the mediated
  $\Gamma_q(x)\thinspace W(\text{midpoint})$ bias — the erase, licensed
  and pump-proportional, vanishing at $V = 0$ as it must.

The crossover prediction of §13.2 extends: the samples of one element
are mutually coherent (phase continuity writes them all with the
element's argument), so a well-sampled element is itself a small
reservoir, and the linear response in the *collective* stored amplitude
is the pair-ensemble form of the mediated rule. As samples decohere or
deplete, linear response degrades to quadratic — now a statement about
the state's own coherence, not only the sea's.

## 6. The compound channel [open]

Mate-sideband strikes execute two edges at once (erase $(P_a, P_b)$,
write $`(P_b \pm 2q\thinspace dp, P_b)`$) at order $\mu_1 \times$ stored.
Theorem C2 achieves exactness without them, so their net contribution
must vanish. The conjectured mechanism is the direction-symmetric
cancellation already verified for the $O(C^2)$ class in the contact-vertex
demo (R2 there); the pair-ensemble integrator should confirm it
stochastically.

## 7. Representation corollaries

Every rung is dynamically reachable, so no truncation is imposed — but
none is needed either: even rungs sit on even midpoint rows, where the
estimator's signed content can equally be carried by uniformly offset
(gray) pairs distributed in place, the two being dual samplings of the
same Wigner content. The specification's disjoint-parities remark holds
for rung 1 only and is corrected in the spec patch accompanying this
note. The even-$`q`$ gap narrows correspondingly: interior mediation is
sea-struck at any midpoint parity; what remains open for even $q$ is
only the population-boundary bookkeeping.

## 8. Consequences and open items

1. The specification acquires a ladder addendum: sea strikers promoted
   from optional to structural; erase amplitude corrected to the
   licensed bilinear form; parity remark fixed; Corollary 4.2 scoped.
2. **Open:** striker back-reaction neutrality (§4); compound-channel
   cancellation (§6); the licensing statement as a lemma rather than an
   argument; the bare-churn balance at $V = 0$ beyond expectation; the
   steady state, unchanged.
3. The decisive integrator test of the pairing note now has a sharper
   target: reproduce Theorem C2's fluxes stochastically, with rung
   census and $V = 0$ freeze as sub-diagnostics.
