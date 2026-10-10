# Depletion and excess: when the sea could be seen

> Empty space is not empty here: a cell that holds no `W` still holds aligned pairs at density `B`. Could the sea's depletion or excess ever be seen? By Proposition Q1 the observable never reads the sea, so only a rule that reads it could show either; a rule that reads its level, Theorem S5's throttle, is already excluded, so excess is unobservable. **Proposition V1**: admissible negativity is free. For a potential of degree two or less the compensated kernel vanishes identically (Theorem G2), so a cat state held in a harmonic trap makes no event and touches no pair; the sea pays only for the negativity the dynamics makes, where `V''' != 0`, and in an anharmonic well a cat makes events 1.4 times as fast as one packet because its fringes are parents too. **Proposition V2**: from a minimal preparation the sea's debt is at least the growth of the negative mass of `W`, with equality exactly when the ensemble stays minimal; body inflation is sea debt, so a finite sea would make garbage collection a requirement. **Proposition V3**: under a floor — no pair, no ionisation, Proposition S0 read literally — `E` is unchanged bit for bit while the sea is at least as deep as the depth `lambda*` the unfloored run needs, and departs below it through a spurious force, an energy drift and a change of negativity, never through the norm. On the Eckart summit `lambda*` is 1.10 under (S) and 0.65 under (S′); in an anharmonic well it climbs without levelling off, to 2.2 and 2.9 by `t = 96`, and it grows with the reach as Theorem G5's census does, which answers Q-SP2. **Proposition V4**: the sea is granular, and the event aperture of Q9(c) holds one pair on average at the physical density, so an integer sea of `beta` pairs per aperture blocks close to `e^-beta` of emissions and moves `E` by `C e^-beta`, `C` from 1.4 to 6.8 — by 37 to 88 per cent at `beta = 1`. The physical density cannot be a strict supply; a sea that is one must be seven to nine pairs per aperture deep to hide below a part in a thousand, and only an experiment could bound that depth, from below. **Proposition V6**: integer worlds at the physical density, over many seeds, block emissions as about `0.9 e^-0.55beta`, not `e^-beta`, because body inflation drains the apertures a parent keeps returning to; the observable still follows V4's mean field, but blocking fewer than one emission in a thousand takes about twelve pairs per aperture. **Proposition V5**: entanglement made by a pair potential costs the sea exactly the negativity of the relative motion. **Proposition V7**: when the centre of mass and the relative motion couple, V2 and V3 hold unchanged in four dimensions, and two bodies in an anharmonic well need a deeper sea than one under every realisation; the physical density is short in all three. **Propositions V8 to V10** (second addendum) amend (S′): an aligned pair changes momentum only in collisions with other pairs, Cyganski's Focus and Defocus at mass-action rates with detailed balance. Collisions are `W`-null (V8). A reversible sea dynamics, such as a Volterra lattice on the pairs, conserves a depletion functional and cannot repair depletion (V9); collisions obey an H-theorem and relax the sea (V10), and in every case run at the default reach the depth the floor needs stops growing below the physical density — 0.25 against 2.9 in the well at `t = 96` — whether they run at the kernel's rate or at a uniform rate between neighbouring rows. **Propositions V11 to V14** (third addendum) test the amendment. At the kernel's rate the depth stays below one half on the Eckart summit from reach `pi` to `8pi`, where (S′) unamended needs up to 17; a uniform rate, which does not keep pace with the events, fails in the well at `4pi` (V11). Of the clock rules tried at a collision, only whole pairs re-locked to their new row keep the sea dark and the readings live (V12). In an integer sea collisions steepen blocking to about `1.1 e^-0.74beta`, about 9.5 pairs per aperture for one emission in a thousand (V13). Two coupled bodies need 0.20 to 0.54 of the physical density with collisions, against 1.86 to 3.98 without (V14). Where to look is where negativity is made fastest — diffraction, Kerr revivals, tunnelling, collisions — but the model's regulators have no laboratory values, so no bound follows yet; a null result is expected, and says the sea's depth, like its motion and its phases, is not fixed by the observable. Step 16's ledger books transport ringing as sea traffic (S-SP7), over a quarter of the sea's debits on Theorem S8's run; the minimal ledger used here does not.

*Ladder abstract — see the [full list](README.md#the-ladder).*

**Status.** Analysis note, step 24 of the ladder. Companion demos:
`src/demo_sea_depletion.py` (Parts A to H), `src/demo_integer_sea.py`,
`src/demo_fourd_sea.py` and `src/demo_pair_collision_clocks.py`. Nothing here depends on
choosing between (S) and (S′); where a number does, both are given.

---

## 0. What this note asks, settles and leaves open

Empty space is not empty in this ontology: a cell of phase space that
holds no $`W`$ still holds aligned pairs at density $`B = 1/\pi\hbar`$.
Step 23 suggested that the sea is a resource with a level of its own — it
is debited by emission, credited by absorption, its deficits can be shallow
and unrefilled (Proposition Q7) — and asked in Q-SP2 whether supply for
emission ever becomes limiting. This note asks the question behind that
one: could depletion or excess of the sea ever be *seen*?

Proposition Q1 of step 23 frames the answer. In every ledger of the chain
the sea enters an event only as its source or its sink and is never read,
so $`E`$ is blind to it. A sea level is therefore observable only through
a rule that reads it, and the question becomes which such rules there
are, what each would cost the observable, and how deep a sea would hide
it.

**Settles.**

- **Proposition V1.** Admissible negativity is free: for a potential of
  degree two or less the compensated residual kernel vanishes identically
  (Theorem G2), so free and harmonic evolution make no event and touch no
  pair, whatever the state. A cat state held in a harmonic trap costs the sea nothing; the
  sea is spent only on what $`V'''`$ does, in proportion to the bodies it
  acts on (§2).
- **Proposition V2.** The sea's debt is the negativity the dynamics makes:
  from a minimal preparation, its global debit is at least the growth of the
  negative mass of $`W`$, with equality exactly when the ensemble stays
  minimal. Body inflation is sea debt (§3).
- **Proposition V3.** The floor. If an emission must find a pair in its
  parent's cell, $`E`$ is unchanged bit for bit while the sea is deep
  enough, $`\beta \ge \lambda^*`$, and departs below it — through a
  spurious force, an energy drift and a change of negativity, never
  through the norm. The depth needed grows with time in a bound
  anharmonic system and with the reach (§4).
- **Proposition V4.** Granular supply. At the physical density the event
  aperture of Q9(c) holds one pair on average. An integer sea of depth
  $`\beta`$ pairs per aperture blocks a fraction close to $`e^{-\beta}`$
  of emissions, and the expected $`E`$ departs by $`C e^{-\beta}`$, with
  $`C`$ between 1.4 and 6.8: by 37 to 88 per cent at the physical density
  (§5).
- **Proposition V6.** The integer sea at the physical density, measured
  over many seeds: blocking falls as about $`0.9\,e^{-0.55\beta}`$, not
  $`e^{-\beta}`$, because body inflation drains the apertures, while the
  observable follows V4's mean field where resolved (§5.1).
- **Proposition V5.** Two bodies. Entanglement made by a pair potential
  costs the sea exactly the negativity the relative motion makes, and per
  joint cell never more than the one-dimensional problem (§6).
- **Proposition V7.** Two coupled bodies. V2 and V3 hold in four
  dimensions, and the depth the floor needs exceeds the one-body value in
  every realisation (§6.1).
- **Postulate (S′), amended, and Propositions V8 to V10.** Under (S′) the
  motion conserves every row's content, which is why the depth the floor
  needs keeps growing. Amended, an aligned pair changes momentum only in
  collisions with other aligned pairs (Definition (C)), which are
  $`W`$-null (V8). A reversible sea dynamics cannot help: a Volterra
  lattice conserves a depletion functional (V9). Collisions obey an
  H-theorem (V10), and the depth the floor needs stops growing, below the
  physical density, in every case run at the default reach (§§4.1 to 4.3).
- **Propositions V11 to V14** (third addendum) measure the amendment. At
  the kernel's rate the depth stays bounded against the reach on the
  Eckart summit to $`8\pi`$; a uniform rate fails in the well at
  $`4\pi`$ (V11, §4.4). Whole pairs re-locked to their new row are the
  only clock rule tried that keeps the sea dark and the reading live
  (V12, §4.5). Collisions steepen an integer sea's blocking to about
  $`1.1\,e^{-0.74\beta}`$ (V13, §5.2), and bring the depth two coupled
  bodies need below the physical density (V14, §6.2).
- What an experiment could and could not see (§7), and a quantification of
  a known demo defect — step 16's ledger books transport ringing as sea
  traffic — that the minimal ledger used here avoids (§8). §8.1 records two
  defects the long runs found: the event step's growth at large reach, and
  collisions on cells in debt.

**Leaves open** the items of §10, chief among them the regulators'
laboratory values (V-SP3), without which no depth can be compared with an
experiment. §5.1 and §6.1, added in October 2026, answer V-SP2 and V-SP1;
§§4.1 to 4.3, a second addendum, amend (S′) and open V-SP8 to V-SP10.
§§4.4, 4.5, 5.2, 6.2 and 8.1, a third addendum, answer V-SP10 and parts of
V-SP8 and V-SP9, and open V-SP11 to V-SP13. A tutorial on the chord and
its reading, which V12 uses, is the supplement
[`chord_reading.md`](../supplement/chord_reading.md).

**Inherits.** Theorems S0, S5 and S7 from step 16; Theorems G2 and G5
from step 17; Proposition M1 from step 21; Propositions Q1, Q2, Q6, Q7, Q8 and Q9 from step 23, and step 23's picture of the
sea as a random configuration of aligned pairs (§1 there). The second
addendum also inherits Cyganski's Focus and Defocus and the no-go lemma of
[`four_rule_microdynamics_equivalence.md`](four_rule_microdynamics_equivalence.md),
and Proposition X1 of the Poisson-kick supplement.

---

## 1. Four ways the sea could be read

| rule | reads | sees a deficit | sees an excess | standing |
|---|---|---|---|---|
| none: the sea is a ledger | nothing | no | no | Q1: invisible, bit for bit |
| its phases (dark catalysis, the live reading) | clock phases | — | — | Q12: delivers the kernel in mean, at best |
| its level (S5's throttle, L2(a)'s weighted kernel) | $`S/B`$ everywhere | yes | yes | S5: relative error 0.25 to 0.30 at the physical depth |
| a floor | whether a pair is there | yes | no | this note |

A rule that reads the level sees excess as well as deficit, and it is
already excluded where it would matter: Theorem S5 throttled the rate by
$`S/B`$ and moved $`E`$ by a quarter in the packet's core, with norm and
$`\langle p \rangle`$ intact (re-measured in Proposition Q3). Tunnelling
through a barrier is not observed to depart from the Schrödinger equation
by anything like that. So **excess is unobservable** unless the level is
read, and reading the level is ruled out.

The floor is different in kind. It does not weight anything; it only
refuses an emission that has nothing to ionise. Proposition S0 is the
reason to take it seriously: an emission is the ionisation of a pair at
the parent's own row, so "no pair, no ionisation" is the ontology read
literally rather than a rule added to it. The rest of this note is about
the floor.

An analogy fixes what kind of answer to expect. The QED vacuum responds
linearly to weak fields and nonlinearly near the Schwinger field
$`E_c = m_e^2 c^3/e\hbar`$, a scale fixed by known constants; vacuum
birefringence, which experiments at the European XFEL aim to see, is
that nonlinearity (Ahmadiniaz et al. 2025). The sea, read through a
floor, is the same shape of hypothesis: linear (exactly Moyal) while it
is deep enough, departing when demand reaches supply. The difference is
that the sea's scale is its depth, which nothing so far fixes (step 23,
§9.4(g)).

---

## 2. Proposition V1: admissible negativity is free

**Proposition V1.** *(a) If $`V`$ is a polynomial of degree at most two,
the compensated residual kernel vanishes identically, for every reach and
every horizon window. (b) Hence under free or harmonic evolution no event
happens and no sea pair is touched, whatever $`W`$ is: the negativity of an
admissible state costs nothing to hold. (c) Where $`V''' \ne 0`$ the sea is
spent at the per-parent rate $`\Gamma(x)`$, so in proportion to the bodies
there; in the minimal ensemble these number $`\lVert W \rVert_1`$, and a
state's fringes are bodies too.*

*Proof.* (a) is Theorem G2 of step 17, which measured the channel empty
at every reach; the reason is short. For such $`V`$,
$`V(x+y) - V(x-y) = 2yV'(x)`$ exactly, so the
Moyal symbol is $`i s V'(x)`$ with $`y = \hbar s/2`$. The window
$`w(y)`$ multiplies it and the classical symbol $`i s`$ alike, so the
compensation, the ratio of their discrete first moments (specification
3.2), returns $`V'_{\rm eff} = V'`$ and the residual
$`i s w (V' - V'_{\rm eff})`$ is zero. (b) Events fire at $`\lvert K_q
\rvert`$ per parent. (c) The rate is per parent, and the minimal ensemble
has $`N = \lVert W \rVert_1 = 1 + 2M_-`$ bodies, $`M_-`$ the negative
mass of $`W`$. $`\square`$

*Verification* (Part A). On the mesh, $`\max\lvert K \rvert / \max\lvert
V'_{\rm eff} \rvert`$ is $`1.3\times10^{-14}`$ for the harmonic potential
against 1.4 for the Eckart barrier and 1.1 for the Pöschl–Teller well. An
even cat, $`x_0 = 3`$, $`\sigma = 1/\sqrt2`$, over one period, in the
minimal ledger of §3:

| state, potential | $`M_-(0)`$ | $`M_-(2\pi)`$ | events | per unit time | sea debit $`D`$ | worst cell $`S/B`$ |
|---|---|---|---|---|---|---|
| cat, harmonic $`\omega = 1`$ | 0.326 | 0.337 | $`7\times10^{-14}`$ | 0 | 0 | 1.000 |
| cat, Pöschl–Teller well | 0.326 | 0.931 | 14.3 | 2.27 | 0.240 | 0.341 |
| one packet at $`x_0 = 3`$, same well | 0 | 0.575 | 10.0 | 1.60 | 0.218 | 0.211 |

In the trap the cat's $`E`$ agrees with the QLE to $`1.2\times10^{-14}`$
and the sea is untouched; its negative mass drifts by 0.011, which is the
discrete transport's own (§3). In the well the cat makes events 1.4 times
as fast as a single packet with the same width and offset — its fringes
are bodies, and every body is a parent — and its negative mass nearly
triples.

![An even cat at t = 0, after one period in the harmonic trap and in the Pöschl–Teller well](https://raw.githubusercontent.com/billpage/wpmw/output/figures/sea_depletion_free_negativity.png)

**What this says about testing the sea.** A cat state *held* in a trap is
not a test: its negativity is admissible initial data, carried by negatons
from the start, and harmonic motion moves it without an event. The ledger
charges for *making* negativity, and in a laboratory a cat is always made
by something non-quadratic — the single-photon Kerr effect (Kirchmair et
al. 2013), a dispersive coupling to an atom or a qubit (Vlastakis et al.
2013), the collapse and revival of a condensate's matter-wave field
(Greiner et al. 2002). That is where a sea would pay. The model takes a
preparation as given data; a complete account must charge it.

---

## 3. Proposition V2: the sea's debt is the negativity made

**Proposition V2.** *Prepare the minimal ensemble, $`N_0 = \lVert W_0
\rVert_1`$, with the sea at $`B`$, and let $`D(t) = S_{\rm tot}(0) -
S_{\rm tot}(t)`$ be the sea's global debit. Under any realisation of the
events (absorptive, emissive or mixed, Theorem S7) and any contact sink,*

```math
D(t) \;\ge\; M_-(t) - M_-(0),
```

*with equality at all times exactly when the ensemble stays minimal,
$`N = \lVert W \rVert_1`$ — that is, under garbage collection
(Proposition Q8). Corollaries: (i) the sea's global excess never exceeds
$`M_-(0)`$; (ii) any body beyond the minimal ensemble is half a pair of
sea debt, so a sea that is a finite supply makes garbage collection a
requirement rather than an economy; (iii) in terms of the negativity
indicator $`\delta = \lVert W \rVert_1 - 1`$ of Kenfack and Życzkowski
(2004), the minimal ensemble's debt is $`\Delta\delta/2`$.*

*Proof.* Event by event $`\Delta S = -\Delta N/2`$ (Theorem S7), and a
contact recombination removes two bodies and adds one pair, so the same
holds for the sink; hence $`D = (N - N_0)/2`$. In every cell $`u_+ + u_-
\ge \lvert u_+ - u_- \rvert`$, so $`N \ge \lVert W \rVert_1`$, with
equality iff no cell holds both species. Since $`\int W = 1`$,
$`\lVert W \rVert_1 = 1 + 2M_-`$, and $`N_0 = 1 + 2M_-(0)`$.
(i): $`D \ge -M_-(0)`$. $`\square`$

Proposition Q8 found the equality case; the inequality says what
inflation costs. A sea of finite depth is debited exactly by the
non-classicality the dynamics creates, and by nothing else unless the
ensemble is wasteful.

*Verification* (Part B). The minimal ledger keeps $`u_\pm = E_\pm`$ after
every step, so only $`E`$ is transported, and $`E`$ is the signed field
spectral transport is accurate for. The discrete transport does not
conserve the sampled $`\lVert E \rVert_1`$ exactly; that change,
$`\Delta N_{\rm tr}`$, is booked in its own column and never charged to the
sea. Then $`D = \Delta M_- - \Delta N_{\rm tr}/2`$ should hold to
round-off, and does:

| case | $`t`$ | $`D`$ | $`\Delta M_-`$ | $`\Delta N_{\rm tr}/2`$ | $`D - \Delta M_- + \Delta N_{\rm tr}/2`$ |
|---|---|---|---|---|---|
| Eckart summit | 2 | 0.0767 | 0.0828 | 0.0062 | $`5\times10^{-13}`$ |
| | 4 | 0.1336 | 0.2537 | 0.1201 | $`9\times10^{-13}`$ |
| | 8 | 0.1727 | 0.3474 | 0.1746 | $`2\times10^{-12}`$ |
| Pöschl–Teller well | 3 | 0.1375 | 0.2227 | 0.0853 | $`1\times10^{-12}`$ |
| | 6 | 0.2579 | 0.4795 | 0.2216 | $`2\times10^{-12}`$ |
| | 12 | 0.3326 | 0.7798 | 0.4472 | $`3\times10^{-12}`$ |

The rows are identical under (S) and (S′), as Proposition Q1 requires of
totals.

The transport's share is not small on this grid — on the Eckart barrier
it is about half of $`\Delta M_-`$ by $`t = 8`$ — and it is a property of
the mesh reference, not of the ledger. Refining the grid at fixed reach
(Part R; Eckart summit, (S′), $`t = 8`$):

| $`n_r`$, $`dr`$ | $`n_p`$, $`dp`$ | $`D`$ | $`\Delta M_-`$ | $`\Delta N_{\rm tr}/2`$ | $`\lvert E - {\rm QLE}\rvert`$ |
|---|---|---|---|---|---|
| 128, 0.3125 | 64, 0.25 | 0.1727 | 0.3474 | 0.1746 | 0.013 |
| 256, 0.1563 | 64, 0.25 | 0.1849 | 0.3376 | 0.1528 | 0.013 |
| 256, 0.1563 | 128, 0.125 | 0.1974 | 0.2990 | 0.1016 | 0.0045 |
| 512, 0.0781 | 128, 0.125 | 0.2039 | 0.2820 | 0.0781 | 0.0046 |

The transport's share falls as the grid is refined, in $`dp`$ most of
all, while the sea's debit $`D`$ settles near 0.2. That is the negative
mass the kernel makes on this collision.

---

## 4. Proposition V3: the floor

**Definition (F) — the floor.** An emission ionises a pair in its parent's
cell. If the sea there is short the emission does not happen, and its two
depositions are lost. Absorption, and the contact sink's credit, are
unchanged.

**Proposition V3.** *Let the sea start at depth $`\beta B`$. In the
unfloored run started at $`B`$, let $`m(T)`$ be the least value of
$`S - E_m`$, over every emission $`E_m`$ due up to time $`T`$ and the
sea $`S`$ its cell holds when it is due, and set
$`\lambda^*(T) = 1 - m(T)/B`$, the depth the floor needs. (a) Under (F),
$`E`$ equals the unfloored $`E`$ bit for bit if
$`\beta \ge \lambda^*(T)`$, and differs from it otherwise. (b) Below
$`\lambda^*`$ the norm is still exact; $`\langle p \rangle`$, the energy,
the negative mass and the transmission move.*

*Proof.* (a) By Proposition Q1 nothing in the unfloored run reads the sea,
so the sea started at $`\beta B`$ is the sea started at $`B`$ plus the
constant $`(\beta - 1)B`$, which transport carries unchanged. (F) acts iff
some emission finds its cell short, $`(\beta - 1)B + S - E_m < 0`$, that
is iff $`\beta < \lambda^*`$. (b) Every event deposits $`+\tau`$ at
$`p + \xi_q`$ and $`-\tau`$ at $`p - \xi_q`$, so a lost event costs no
norm. The compensation makes the residual's discrete first moment vanish
at every $`x`$, $`\sum_q \xi_q K_q(x) = 0`$; a rule that removes some
channels at some parents and not others leaves that sum unbalanced, which
is a force no potential exerts. By Q9(a) each event moves the
$`W`$-weighted kinetic energy by $`2p\tau\xi_q/m`$, so lost events are an
energy drift. $`\square`$

*Verification* (Part C), the minimal ledger, $`\Delta t = 0.02`$; the
Eckart summit packet to $`T = 8`$ and a packet displaced in the
Pöschl–Teller well $`-4\,{\rm sech}^2(x/2)`$ to $`T = 12`$; differences
are against the unfloored run of the same mode:

| case | realisation | $`f`$ | $`\lambda^*`$ | $`\varepsilon`$ at $`\beta = 0.5`$ | at $`\beta = 1`$ | at $`0.98\,\lambda^*`$ | at $`1.02\,\lambda^*`$ |
|---|---|---|---|---|---|---|---|
| Eckart | (S), absorptive first | 0.268 | 1.099 | 0.110 | 0.019 | $`1.9\times10^{-3}`$ | 0 |
| | (S′), absorptive first | 0.268 | 0.648 | 0.022 | 0 | $`2.1\times10^{-3}`$ | 0 |
| | (S′), emissive | 0 | 1.094 | 0.207 | 0.010 | $`1.4\times10^{-3}`$ | 0 |
| Pöschl–Teller well | (S), absorptive first | 0.299 | 1.295 | 0.451 | 0.144 | $`6.0\times10^{-3}`$ | 0 |
| | (S′), absorptive first | 0.299 | 0.765 | 0.081 | 0 | $`1.6\times10^{-3}`$ | 0 |
| | (S′), emissive | 0 | 0.994 | 0.396 | 0 | $`2.8\times10^{-3}`$ | 0 |

Here $`\varepsilon = \lVert E_{\rm F} - E \rVert / \lVert E \rVert`$, and
0 means bit for bit; the unfloored $`E`$ is the QLE to 0.013 on the summit
and 0.09 to 0.11 in the well (the τ-leap's first-order error at
$`\Delta t = 0.02`$). The departures carry the signatures of (b). On the
summit under (S) at $`\beta = 0.5`$: norm $`2\times10^{-15}`$,
$`\Delta\langle p \rangle = -0.004`$, $`\Delta H = -0.003`$,
$`\Delta M_- = -0.019`$, $`\Delta T = -0.0014`$
(`sea_depletion_floor.csv` has every row).

Three things follow. The threshold is sharp, as (a) says: two per cent
below $`\lambda^*`$ the departure is a few parts in a thousand, two per
cent above it is nothing. Under (S′) the physical depth suffices for both
runs, $`\lambda^* < 1`$, as Proposition Q3's shallow deficits suggested;
under (S) it does not, by 10 and 30 per cent. And the realisation
matters: with absorption first the sea is credited where it is debited,
at the parent's row; with emission only, the contact sink credits the
daughters' rows instead, and $`\lambda^*`$ rises by a third to two
thirds (V-SP4). The sign
of $`\Delta M_-`$ is not fixed: a lost deposition can land where it would
have deepened a negative region or where it would have filled one, and in
the emissive well at $`\beta = 0.25`$ the starved dynamics runs away
($`\varepsilon = 3.6`$).

**Long runs** (Part E, open item Q-SP2): $`\lambda^*(t)`$ to $`t = 96`$
(`--parts E --t-long 96`; the default run stops at 24 and agrees on the
columns they share):

| case | realisation | $`t = 8`$ | 16 | 24 | 32 | 48 | 64 | 96 |
|---|---|---|---|---|---|---|---|---|
| Pöschl–Teller well | (S), absorptive first | 1.295 | 1.364 | 1.614 | 1.614 | 2.046 | 2.046 | 2.230 |
| | (S′), absorptive first | 0.727 | 0.908 | 1.127 | 1.220 | 1.824 | 2.551 | 2.902 |
| | (S′), emissive | 0.772 | 1.207 | 1.449 | 1.760 | 2.980 | 3.624 | 4.494 |
| Eckart summit | (S), absorptive first | 1.099 | 1.099 | 1.100 | 1.867 | 1.867 | 1.867 | 1.867 |
| | (S′), absorptive first | 0.648 | 0.651 | 0.679 | 1.099 | 1.951 | 2.053 | 2.582 |
| | (S′), emissive | 1.094 | 1.096 | 1.099 | 1.411 | 2.289 | 3.425 | 3.640 |

In the well, where the packet keeps dephasing and keeps making
negativity, $`\lambda^*`$ climbs throughout and has not levelled off by
$`t = 96`$: no fixed depth suffices for ever. It climbs faster under
(S′) after the first few periods, as Proposition Q7 led one to expect:
a deficit that is not sheared across rows is never refilled. The figure
shows the difference: under (S) the worst cell after each step (thin
line) springs back between spikes, and $`\lambda^*`$ is set by the
spikes; under (S′) the worst cell stays near $`\lambda^*`$. The Eckart
packet leaves the barrier and $`\lambda^*`$ stops at its single-collision
value by $`t = 8`$; on the periodic box its two halves return to the
barrier from $`t \approx 30`$, and the later columns are repeated
collisions, which behave like the well.

Against the reach (Eckart summit, $`dp = 0.125`$, $`T = 8`$, window cut
at $`y_h`$), $`\lambda^*`$ grows too:

| $`y_h/a`$ | $`\max\Gamma_{\rm tot}`$ | (S), absorptive first | (S′), absorptive first | (S′), emissive |
|---|---|---|---|---|
| $`\pi`$ | 1.33 | 0.550 | 0.437 | 0.457 |
| $`2\pi`$ | 3.49 | 1.332 | 0.748 | 1.233 |
| $`4\pi`$ | 6.58 | 1.504 | 1.148 | 1.754 |

This is Theorem G5 of step 17 seen from the sea: there the census of
events grows with the reach while the generator converges, and here so
does the depth a floor needs. If the reach is a regulator to be removed, a
sea read through a floor needs ever more depth as the regulator goes. Only a
physical reach (Q-SP6) would give the floor a finite threshold.

![The depth the floor needs, against time](https://raw.githubusercontent.com/billpage/wpmw/output/figures/sea_depletion_long.png)

### 4.1 Pair collisions: an amended (S′)

*Addendum (October 2026).* Part E leaves a problem. Under (S′) the depth
the floor needs climbs without limit in a bound anharmonic system, so no
finite sea suffices for ever. The cause is local, not global. The sea's
global debt is bounded by Proposition V2, and in the well it reaches 0.29
by $`t = 8`$ and stays below 0.48 to $`t = 96`$, while the worst cell
keeps falling.
Under (S′) the sea at fixed momentum is carried along its row,

```math
\partial_t s + \frac{p}{m}\,\partial_x s = \sigma(x, p, t),
\qquad\text{so}\qquad
\frac{d}{dt}\oint s(x, p, t)\,dx = \oint \sigma(x, p, t)\,dx
```

on the periodic box, with $`\sigma`$ the events' debits and credits
(emissions at their parents, absorptions and contact recombinations
where they happen). Every row's content is a conserved quantity of the
motion, so any imbalance between where a row is debited and where it is
credited piles up there. Under (S) the force shears rows into one
another and part of the deficit refills (Proposition Q7), which is why
$`\lambda^*`$ grows more slowly there.

By Proposition Q1 nothing observable constrains how the sea moves: the
QLE fixes $`u_+ - u_-`$ and says nothing about the pairs. So the sea may
be given a dynamics of its own, provided it moves whole pairs. This
addendum amends (S′) accordingly, after asking what kind of dynamics
could help.

**Postulate (S′), as amended.** Between events every free body obeys the
full classical force, as in (S). An aligned pair is force-blind: between
its events it moves inertially, $`\dot q = p/m`$, $`\dot p = 0`$, its
clock winding by P1, and it changes momentum only in a collision with
another aligned pair, by Definition (C).

**Definition (C) — pair collisions.** Two aligned pairs in one cell of
position collide in one of two ways, on a channel $`q`$ (a momentum
step $`\xi_q = q\,dp`$):

- **Defocus**: two pairs on row $`n`$ go to one pair on row $`n - q`$
  and one on row $`n + q`$;
- **Focus**: one pair on each of rows $`n \mp q`$ go to two pairs on row
  $`n`$.

The rates are mass action with detailed balance: per unit area, the net
rate of defocus over focus about row $`n`$ is

```math
J_{n,q} = r_q(x)\,\bigl(s_n^2 - s_{n-q}\,s_{n+q}\bigr),
\qquad
\partial_t s_n\big|_{\rm C} = \sum_q \bigl(-2J_{n,q} + J_{n-q,q} + J_{n+q,q}\bigr),
```

with $`s_n`$ the pair density on row $`n`$. The rate law $`r_q(x)`$ is a
**[choice]**. Two are measured here:

- **kernel rate**: every channel, $`r_q B = \lvert K_q(x)\rvert/2`$, so that
  collisions happen where events do and at their rate;
- **uniform**: neighbouring rows only ($`q = 1`$), $`r B = \gamma`$,
  independent of $`x`$ and of $`V`$, a property of the sea alone.

These are Cyganski's Focus and Defocus
([`four_rule_microdynamics_equivalence.md`](four_rule_microdynamics_equivalence.md))
moved from the signed occupancy to the dark pairs. There they had to be
signed and linear in the occupancies to reproduce the QLE, and the no-go
lemma of §6.2 there excludes mass action. Here nothing has to be linear,
because nothing observable reads the sea.

**Proposition V8 (collisions are W-null).** *Collisions move whole pairs,
so in any ledger whose events take their counts from the bodies and the
kernel alone, $`E`$, $`N`$ and $`f`$ are those of (S′) bit for bit.*

*Proof.* Proposition Q1: the sea enters an event only as its source or
its sink, and a collision changes neither body field. $`\square`$

Verified: in every run of Part G, $`\max\lvert\Delta E\rvert = 0`$
against (S′) unamended. Under a rule that reads the sea, Definition (F)
included, collisions do change $`E`$, and Proposition V3's argument that
a deeper sea is the same sea plus a constant no longer holds, since the
collision rates read $`s`$. The statement that survives is the one that
matters: at the physical depth the floor never binds if $`\lambda^* < 1`$.

### 4.2 Proposition V9: a reversible sea cannot repair depletion

The first candidate for sea dynamics was Cyganski's Volterra sea (call of
6 October 2026): his copy-left and copy-right actions, in which two
bodies on one row join two on the next, at mass-action rates. He applied
them to positons and negatons separately, which changes $`E`$; that is a
dynamics of the bodies, and in the large-sea limit it reproduces the
single-mode QLE kick. Applied to whole pairs it is W-null, so it is a
candidate here. It does not help, and the reason is structural.

**Proposition V9.** *(a) Let the pair density on one row evolve by any
Volterra lattice*

```math
\partial_t s_n = \frac{s_n}{B}\sum_q c_q(x)\,\bigl(s_{n-q} - s_{n+q}\bigr),
```

*with signed rates $`c_q(x)`$ on any set of channels. Then $`\sum_n s_n`$
and $`\sum_n \ln s_n`$ are conserved on every row, so the depletion
functional*

```math
\Phi[s] = \iint \Bigl[s - B - B\ln\frac{s}{B}\Bigr]\,dx\,dp \;\ge\; 0,
\qquad \Phi = 0 \iff s \equiv B,
```

*is conserved, and streaming under (S) or (S′) conserves it too. Only
events change $`\Phi`$, at the rate
$`\iint \sigma\,(1 - B/s)\,dx\,dp`$, so a debit landing in a depleted cell
always increases it. An empty cell stays empty, since
$`\partial_t s_n \propto s_n`$. (b) The linearisation about $`s = B`$
has an antisymmetric generator and conserves $`\lVert s - B\rVert_2`$.
(c) Measured, with the compensated residual as the linear rates
($`c_q B = K_q`$): the depth the floor needs is halved, and still grows.*

*Proof.* (a) $`\sum_n \partial_t \ln s_n = B^{-1}\sum_q c_q \sum_n
(s_{n-q} - s_{n+q}) = 0`$ on a periodic row, and
$`\sum_n \partial_t s_n = B^{-1}\sum_q c_q \sum_n (s_n s_{n-q} -
s_n s_{n+q}) = 0`$ by a shift of $`n`$. $`\Phi`$ is a sum of these two
over $`x`$, and streaming is measure-preserving, so it conserves the
integral of any function of $`s`$. (b) $`K_{-q} = -K_q`$. $`\square`$
(a) is checked exactly with SymPy on a periodic row of seven cells with
two channels, and numerically on a row of 64 with six channels and random
signed rates: $`\Phi = 4.8890935456`$ at $`t = 0, 5, 10, 20`$ while the
deepest cell rose from $`0.05B`$ to between $`0.27B`$ and $`0.34B`$. A Volterra sea shares a
deficit out among its neighbours as ripples of both signs; it cannot
remove it.

So the sea needs a dynamics that forgets: a dissipative one. For the
observable that would be forbidden, since a positive process acting on
$`W`$ must diffuse (Pawula 1967; Proposition X1 of
[`../supplement/poisson_kicks_and_pair_branching.md`](../supplement/poisson_kicks_and_pair_branching.md));
the sea is not $`W`$, and nothing forbids it there.

### 4.3 Proposition V10: collisions relax the sea

**Proposition V10.** *Under Definition (C): (a) every event conserves
the number of pairs and their total momentum, and keeps both members of
each pair together; (b) the collision entropy*

```math
\mathcal{H}[s] = \iint \Bigl[s\ln\frac{s}{B} - s + B\Bigr]\,dx\,dp
```

*never increases,*

```math
\frac{d\mathcal{H}}{dt}\bigg|_{\rm C}
= -\sum_{n,q} r_q\,\bigl(s_n^2 - s_{n-q}s_{n+q}\bigr)
\ln\frac{s_n^2}{s_{n-q}s_{n+q}} \;\le\; 0,
```

*with equality iff $`s_n^2 = s_{n-q}s_{n+q}`$ for every active channel,
that is iff $`\ln s`$ is linear along each row; on a periodic row, or on
an open one with bounded $`s`$, that means uniform; (c) linearised about
a uniform row of density $`\bar s`$, each channel acts on a ripple of row
wavenumber $`\theta`$ with the rate*

```math
-4\,r_q\,\bar s\,\bigl(1 - \cos q\theta\bigr)^2 \;\le\; 0,
```

*a hyperdiffusion that damps every non-uniform ripple and keeps the row's
mean and first moment; (d) measured, the depth the floor needs stops
growing, and stays below the physical density in every case run.*

*Proof.* (a) A defocus moves momentum $`2n \to (n - q) + (n + q)`$. (b)
$`d\mathcal{H}/dt = \sum_n \ln(s_n/B)\,\partial_t s_n`$; with
$`\sum_n \partial_t s_n = 0`$ the constant drops out, and shifting $`n`$
in the two $`J_{n\mp q}`$ terms gives the sum shown. Each term has the
form $`(a - b)\ln(a/b) \ge 0`$. This is Boltzmann's H-theorem for
reversible mass-action kinetics with detailed balance (Horn and Jackson
1972). (c) Linearise $`J_{n,q}`$ to
$`r_q\bar s\,(2\delta_n - \delta_{n-q} - \delta_{n+q})`$ and take the
Fourier symbol of $`-2J_n + J_{n-q} + J_{n+q}`$. $`\square`$

(a), (b) and (c) are checked exactly with SymPy (Part G). Numerically, on a
row of 1024 cells with a ripple at its centre, pairs and pair momentum are
conserved to round-off and $`\mathcal{H}`$ falls at every sample. The
demo's momentum grid is periodic, inherited from step 16's spectral mesh,
so a collision across its ends does not conserve momentum; the rows there
stay at $`B`$.

**Measured** (Part G), the Pöschl–Teller well of §4, minimal ledger,
absorptive first, (S′), $`\Delta t = 0.02`$, to $`t = 96`$. The pair-diffusion
control uses the same channels and rates as the kernel-rate collisions,
with a symmetric kernel, $`\partial_t s_n = \sum_q\lvert K_q\rvert
(s_{n+q} + s_{n-q} - 2s_n)`$; the Volterra sea is (c) of V9. Columns are
$`\lambda^*(t)`$; $`\lVert s - \bar s\rVert_2`$ is the sea's deviation
from its row means at $`t = 96`$, in units of $`B`$:

| sea motion | 16 | 32 | 48 | 64 | 80 | 96 | deviation |
|---|---|---|---|---|---|---|---|
| (S′) unamended | 0.908 | 1.220 | 1.824 | 2.551 | 2.837 | 2.902 | 4.08 |
| Volterra sea (V9) | 0.716 | 0.809 | 0.850 | 1.091 | 1.175 | 1.350 | 3.30 |
| pair diffusion (control) | 0.228 | 0.228 | 0.228 | 0.228 | 0.228 | 0.228 | 0.22 |
| collisions, kernel rate | 0.251 | 0.251 | 0.251 | 0.251 | 0.251 | 0.251 | 0.24 |
| collisions, uniform, $`rB = 1`$ | 0.468 | 0.468 | 0.468 | 0.468 | 0.468 | 0.468 | 0.52 |

The global debt is 0.463 at $`t = 96`$ in every row, as V8 requires.
The unamended row reproduces Part E. Collisions at the kernel's rate hold
$`\lambda^*`$ at 0.251 from $`t \approx 8`$; collisions between
neighbouring rows only, at a rate that knows nothing of the potential,
hold it at 0.468.

How fast collisions must be (same run):

| collisions | rate | $`\lambda^*`$ at 16 | 48 | 96 | deviation |
|---|---|---|---|---|---|
| kernel rate | × 1 | 0.251 | 0.251 | 0.251 | 0.24 |
| | × 0.3 | 0.428 | 0.428 | 0.428 | 0.54 |
| | × 0.1 | 0.530 | 0.669 | 0.820 | 1.10 |
| | × 0.03 | 0.650 | 1.117 | 1.670 | 1.99 |
| | × 0.01 | 0.793 | 1.490 | 2.231 | 2.91 |
| uniform, neighbouring rows | $`rB = 4`$ | 0.357 | 0.357 | 0.357 | 0.41 |
| | $`rB = 1`$ | 0.468 | 0.468 | 0.468 | 0.52 |
| | $`rB = 0.25`$ | 0.541 | 0.541 | 0.541 | 0.61 |

At the kernel's rate the bound holds down to about a tenth of it (0.82,
flat from $`t \approx 64`$), and not at three hundredths (1.67, still
rising). Uniform collisions between neighbouring rows hold it down to
$`rB = 0.25`$, about a twentieth of the well's peak event rate
($`\max\Gamma_{\rm tot} = 4.5`$), because they act wherever a deficit
travels and not only where events make it. Under uniform collisions at
$`rB = 1`$ the sea's deviation creeps from 0.29 at $`t = 60`$ to 0.52 at
$`t = 96`$ while $`\lambda^*`$ stays at 0.468; whether it settles is not
shown.

Other cases, collisions at the kernel's rate:

| case | realisation | (S′) unamended at 16 / 48 / 96 | with collisions at 16 / 48 / 96 |
|---|---|---|---|
| Pöschl–Teller well | emissive | 1.207 / 2.980 / 4.494 | 0.354 / 0.354 / 0.375 |
| Eckart, repeated collisions | absorptive first | 0.651 / 1.951 / 2.582 | 0.275 / 0.275 / 0.275 |
| Eckart, repeated collisions | emissive | 1.096 / 2.289 / 3.640 | 0.361 / 0.376 / 0.376 |

Against the reach (Eckart summit, one collision, $`dp = 0.125`$, $`T = 8`$,
window cut at $`y_h`$, as in Part E; the unamended columns reproduce Part
E's):

| $`y_h/a`$ | $`\max\Gamma_{\rm tot}`$ | absorptive first | with collisions | emissive | with collisions |
|---|---|---|---|---|---|
| $`\pi`$ | 1.33 | 0.437 | 0.213 | 0.457 | 0.215 |
| $`2\pi`$ | 3.49 | 0.748 | 0.314 | 1.233 | 0.356 |
| $`4\pi`$ | 6.58 | 1.148 | 0.318 | 1.754 | 0.446 |

From $`2\pi a`$ to $`4\pi a`$ the absorptive-first depth grows by 1 per
cent with collisions against 53 without, and the emissive depth by 25 per
cent against 42. A single collision cannot say whether the emissive
growth stops; the long runs were made at the default reach only (V-SP9,
answered in §4.4).

![Pair collisions bound the depth the floor needs](https://raw.githubusercontent.com/billpage/wpmw/output/figures/sea_collisions.png)

*Figure. Left: $`\lambda^*(t)`$ in the Pöschl–Teller well under each sea
motion, with the rate scan dashed. Right: the sea's deviation from its
row means.*

**What the amendment costs, and what it leaves open.**

- *Energy.* A defocus raises the pairs' kinetic energy and a focus lowers
  it; only detailed balance makes the uniform sea the equilibrium. The
  pairs' energy was already not conserved under (S′) (Q-SP3), and it is
  unobservable either way.
- *An arrow of time in the dark.* $`\mathcal{H}`$ decreases, so the sea's
  own dynamics is irreversible, while $`E`$ stays exactly reversible (V8).
  Dark catalysis (Proposition Q6) is already an irreversible reset of the
  sea's phases; collisions are its counterpart for the population, which
  dark catalysis cannot touch (Proposition Q7).
- *The clock at a collision* (V-SP8). Under (S′) a pair's clock is locked
  to its row's plane wave (Proposition Q2), and the reading of
  Proposition Q5 is taken from those clocks. A collision moves two pairs
  to new rows, so it needs a clock rule, and Q6's resets must keep up with
  the misalignment collisions bring. Nothing in $`E`$ depends on that rule
  (V8); the readings of §§8–9 of step 23 do, and were not re-measured.
- *The rate* (V-SP9). The bound holds at the kernel's rate and at a
  uniform rate unrelated to the potential, and fails when collisions are
  much slower than events. Which rate law is the ontology's, and what
  fixes its constant, is open, as the dark rate is (Q-SP9).
- *What it does not touch.* Proposition V4's granularity is a property
  of an integer sea at one pair per aperture, not of a mean deficit, and
  collisions in mean field say nothing about it (V-SP10). Proposition Q4's
  tagged crossings ran with pairs that keep their rows.

### 4.4 Proposition V11: the reach in long runs

*Third addendum (October 2026), open item V-SP9.* §4.3's reach ladder
followed one collision to $`T = 8`$. Part H runs both cases to $`t = 96`$
at reaches $`\pi`$, $`2\pi`$ and $`4\pi`$, and on the summit $`8\pi`$, in
both realisations, under (S′) unamended, with collisions at the kernel's
rate, and with uniform collisions between neighbouring rows at
$`rB = 1`$. The coarse grid is Part E's, with $`n_p = 128`$ and
$`dp = 0.125`$; the fine grid halves $`dp`$.

**Proposition V11 (measured).** *(a) At the kernel's rate, collisions
hold $`\lambda^*`$ on the Eckart summit at every reach run, and it hardly
grows with the reach: it stays below 0.5 in both realisations from $`\pi`$
to $`8\pi`$. (S′) unamended needs 5.7 (absorptive first) and 17
(emissive) at $`8\pi`$. (b) In the Pöschl–Teller well they hold the
absorptive-first depth to $`4\pi`$ (0.36). The emissive depth at $`4\pi`$
stays at 0.79 until $`t = 80`$ and then rises, by an amount that shrinks
with $`\Delta t`$. (c) Uniform collisions, whose rate does not grow with
the reach, hold $`\lambda^*`$ at $`\pi`$ and $`2\pi`$ but not at $`4\pi`$
in the well. There the emissive depth reaches 6.0 and is still climbing
at $`t = 96`$. (d) In the emissive realisation the sea's deviation from
its row means keeps creeping under either law while $`\lambda^*`$ stays
put. Absorptive first, it settles.*

$`\lambda^*`$ at $`t = 96`$ on the Eckart summit. Under collisions every
entry is flat from $`t = 16`$ unless noted:

| $`y_h/a`$ | grid, $`\Delta t`$ | $`\max\Gamma_{\rm tot}`$ | absorptive first: (S′) | kernel rate | uniform | emissive: (S′) | kernel rate | uniform |
|---|---|---|---|---|---|---|---|---|
| $`\pi`$ | coarse, 0.02 | 1.33 | 0.514 | 0.213 | 0.354 | 0.763 | 0.215 | 0.393 |
| $`2\pi`$ | coarse, 0.02 | 3.49 | 2.356 | 0.314 | 0.674 | 2.589 | 0.356 | 0.960 |
| $`4\pi`$ | coarse, 0.01 | 6.58 | 3.197 | 0.321 | 0.960 at 0.02, from 0.622 | 6.872 | 0.426 | — |
| $`4\pi`$ | fine, 0.01 | 7.13 | 4.097 | 0.487 | | 5.930 | 0.449 | |
| $`8\pi`$ | fine, 0.0025 | 12.92 | 5.671 | 0.425 | | 17.18 | 0.499 | |

The same in the Pöschl–Teller well:

| $`y_h/a`$ | grid, $`\Delta t`$ | $`\max\Gamma_{\rm tot}`$ | absorptive first: (S′) | kernel rate | uniform | emissive: (S′) | kernel rate | uniform |
|---|---|---|---|---|---|---|---|---|
| $`\pi`$ | coarse, 0.02 | 1.36 | 0.439 | 0.160 | 0.282 | 0.638 | 0.221 | 0.391 |
| $`2\pi`$ | coarse, 0.02 | 5.17 | 1.420 | 0.241 | 0.685 | 2.564 | 0.364 | 1.203, from 0.886 |
| $`4\pi`$ | coarse, 0.005 | 12.53 | 5.078 | 0.367 | 1.759 | 9.586 | 1.600, from 0.897 | 7.655 |
| $`4\pi`$ | coarse, 0.0025 | 12.53 | 3.600 | 0.359 | 1.646 | 7.849 | 0.943, from 0.790 | 6.038 |

Read down the kernel-rate columns. On the summit the absorptive depth is
0.213, 0.314 and 0.321 on the coarse grid and 0.487 and 0.425 on the
fine. The emissive depth is 0.215, 0.356 and 0.426, then 0.449 and 0.499.
$`\lambda^*`$ is a worst-cell measure, and the fine grid's cell is half
the size, so compare within one grid. From $`4\pi`$ to $`8\pi`$ the
emissive depth grows by 11 per cent with collisions, against a factor of
2.9 without. That answers the half of V-SP9 that asked whether the
emissive growth stops: at the kernel's rate, on the summit, it does. The
uniform emissive run at $`4\pi`$ diverged before the guard of §8.1 was
added, and was not repeated.

In the well beyond $`2\pi`$ the numbers are limited by the ledger's event
step (§8.1). The global debt, which V8 makes the same under every sea
motion, is 3.55, 0.992 and 0.666 at $`\Delta t = 0.01`$, 0.005 and
0.0025, so even the finest step is not converged. The kernel-rate
absorptive depth agrees at the two finest steps (0.367 and 0.359) and is
resolved. The emissive late rise is not resolved: it is 1.600 at
$`\Delta t = 0.005`$ and 0.943 at 0.0025, from 0.79 before $`t = 80`$,
and is probably the step's (V-SP11). The uniform columns differ from the
kernel's by far more than the step's error: 1.65 against 0.36 absorptive
first, and 6.0 and rising against 0.94 emissive. From $`\pi`$ to $`4\pi`$
the well's peak event rate grows ninefold, 1.36 to 12.5, while the
uniform rate stays at $`rB = 1`$. A rate that does not grow with the
event rate falls behind once the reach makes events fast enough. The
kernel's rate grows with the event rate by construction.

The deviation behaves differently. In the emissive realisation every
event debits the sea at its parent, and every credit comes from the
contact sink, elsewhere (Q8). That realisation therefore makes surpluses
as well as deficits. Collisions bound the deficits, which are what
$`\lambda^*`$ reads; the deviation counts both, and keeps creeping. In
the well at $`\pi`$ it rises from 0.80 at $`t = 32`$ to 1.89 at
$`t = 96`$ under uniform collisions, and from 0.39 to 0.54 at the kernel's
rate. On the summit at $`8\pi`$, kernel rate, it rises from 0.64 to 0.97.
Absorptive first, it settles (summit, $`8\pi`$, kernel rate: 0.24, then
0.20). Whether the emissive surplus settles in longer runs is not shown.
Which rate law is the ontology's is not settled either. What is settled
is that a law must keep pace with the event rate to hold the depth at
large reach.

### 4.5 Proposition V12: the clock at a collision

*Third addendum (October 2026), open item V-SP8.* Definition (C) moves
pair density. Bodies carry clocks (step 23, P1 and P2), so a collision
of two aligned pairs needs a rule for four clocks. A pair whose positon
is at $`x_+`$ and negaton at $`x_-`$ on row $`p`$ has misalignment
$`\mu = \theta_+ - \theta_- + p(x_- - x_+)/\hbar`$ and grayness
$`\sin^2(\mu/2)`$. Four realisations were tried:

- *whole*: each pair moves as a unit, both members by $`\pm\xi`$, clocks
  continuous (Theorem Y4);
- *exchange*: the partners swap, $`(C_+, D_-)`$ one way and
  $`(D_+, C_-)`$ the other, the products at their parents' sites;
- *matched*: exchange, accepted with weight
  $`\cos^2(\mu_1/2)\cos^2(\mu_2/2)`$ of the products' misalignments, the
  orientation (which positon goes up) chosen by the same weight;
- *relock*: whole, then each pair takes the circular mean of its
  destination row's sea clocks within 3, each carried to it by its own
  momentum (step 22 §6).

**Proposition V12.** *(a) A whole-pair collision keeps each pair's
misalignment, so it makes no gray pair. An exchange gives each product
the misalignment of the chord its parents form, read on the daughter row
(Proposition Q5(a)), and*

```math
\mu_1 + \mu_2 = \mu_C + \mu_D + o\,\frac{2\xi d}{\hbar},
```

*with $`o = \pm 1`$ the orientation and $`d`$ the parents' separation.
(b) Measured in step 23's particle model: whole pairs keep the sea dark
but break the row lock, and exchange keeps neither. Only whole pairs
re-locked to their new row keep the sea dark, the lock and the live
reading. (c) The sea can realise a collision itself, as a two-body
event. Let an aligned pair at row $`n`$ fire channel $`q`$ of the kernel
with orientation $`t = \pm 1`$ and move itself to $`n + tq`$. Let it take
a partner in proportion to the counts, either from its own row (the
partner goes to $`n - tq`$: a defocus) or from row $`n + 2tq`$ (both land
on $`n + tq`$: a focus). Then the mean field is Definition (C)'s mass
action, in detailed balance. Realised as a catalysed event — a third pair
fires the kernel and two others are its donors (dark catalysis, Theorem
L5) — a collision is three-body, and no allocation gives $`J = 0`$ on a
uniform row.*

*Proof.* (a) SymPy (Part A of `demo_pair_collision_clocks.py`): carry
each clock to a common point at its row's wavevector (P2) and subtract.
(c) Defocus by the reader goes as $`s_n s_n`$ and focus as
$`s_n s_{n+2tq}`$, which is $`s_{m-q}s_{m+q}`$ about $`m = n + tq`$: the
form of $`J_{m,q}`$. $`\square`$

Part B3 checks (c) numerically on a row lattice. On a uniform row the
drift is zero and no ripple grows. A 30 per cent dip refills to within 4
per cent by $`t = 40`$, and pairs are conserved to $`10^{-13}`$. That
holds with mass-action partners and with an integer sea's chance
$`1 - e^{-m}`$ of finding one. Part B2 runs the catalysed version.
Absorptive first, a uniform row is unstable: the growth rates are 0.21
and 2.6 at $`s = 1`$ and 4 pairs per aperture, and a dip deepens to
$`0.07`$. Emissive first, it relaxes, but carries a flux of 0.49 and 3.9
on a uniform row and has no H-theorem.

A chord cannot be the trigger of a collision one for one. Firing each
chord at the magnitude of its own contribution to $`K_q`$ would run 8.6
to 16.2 times faster than the kernel, summed over channels (Part B1, at
three positions on the summit). The kernel is what survives the chords'
interference. The supplement
[`chord_reading.md`](../supplement/chord_reading.md) works this through,
with the reader of (c) in its §7.

**Measured** (Part C). Eckart summit, (S′), step 23's force-blind
row-centred locked sea at $`\nu = 32`$ with the wrap-phase fix. Events
fire at $`\lvert K_q \rvert`$, dark catalysis at $`\kappa`$ times the
kernel's own rate, and collisions at the kernel rate of Definition (C)
in each realisation. Four seeds to $`t = 14`$, mean ± standard error.
*Gray* is the mean $`\sin^2(\mu/2)`$ of the sea's pairs near the barrier
at the end. *Lock* is the chord coherence
$`\lvert\sum e^{i\mu^{(P)}}\rvert / n`$ over the chords about free bodies
near the barrier. *Live reading* is the slope and correlation of the live
reading, integrated along free bodies, against the mesh kernel (Q12(e)).
*Gray per collision* is the grayness made per collision.

| $`\kappa`$ | realisation | gray | lock | live reading: slope | correlation | gray per collision |
|---|---|---|---|---|---|---|
| 0 | none | 0.027 ± 0.002 | 0.402 ± 0.006 | 0.13 ± 0.03 | 0.20 ± 0.05 | — |
| | whole | 0.028 ± 0.001 | 0.089 ± 0.012 | 0.02 ± 0.02 | 0.05 ± 0.04 | 0 |
| | exchange | 0.377 ± 0.002 | 0.064 ± 0.007 | 0.01 ± 0.004 | 0.04 ± 0.01 | 0.280 |
| | matched | 0.111 ± 0.002 | 0.114 ± 0.017 | 0.04 ± 0.01 | 0.10 ± 0.03 | 0.148 |
| | relock | 0.0105 ± 0.0004 | 0.702 ± 0.006 | 0.53 ± 0.03 | 0.69 ± 0.02 | −0.015 |
| 2 | none | 0.0039 ± 0.0006 | 0.740 ± 0.012 | 0.71 ± 0.05 | 0.76 ± 0.02 | — |
| | whole | 0.0051 ± 0.0012 | 0.221 ± 0.006 | 0.14 ± 0.02 | 0.37 ± 0.03 | 0 |
| | exchange | 0.119 ± 0.003 | 0.194 ± 0.011 | 0.12 ± 0.01 | 0.35 ± 0.02 | 0.570 |
| | matched | 0.028 ± 0.001 | 0.339 ± 0.008 | 0.22 ± 0.003 | 0.49 ± 0.01 | 0.226 |
| | relock | 0.0034 ± 0.0007 | 0.851 ± 0.006 | 0.69 ± 0.05 | 0.80 ± 0.02 | −0.007 |

At $`\kappa = 6`$ the lock is 0.87 with no collisions, 0.47 whole, 0.44
exchange, 0.60 matched and 0.90 relock.

Whole pairs leave every pair as dark as before (0.028 against 0.027),
but they drop the lock from 0.40 to 0.09. A pair arriving on a new row
carries a clock locked to the old one, so its chords with the row's
other pairs are incoherent, although each pair is dark. Step 23 drew
that distinction between a dark pair and a locked one (§§8 and 9 there),
and collisions make it bite. Exchange makes 0.28 to 0.63 of a unit of
grayness per collision. Phase matching halves that, at the cost of
rejecting about half the proposals, and still loses the lock.
Re-locking restores it, above the value without collisions at every
$`\kappa`$ (0.70 against 0.40 with no dark catalysis), because each
re-lock is itself a reset of the kind dark catalysis makes. The live
reading follows: its slope is 0.53 against 0.13 at $`\kappa = 0`$, and
at $`\kappa = 2`$ and 6 it is as good as without collisions. The bodies' observables are those of whole pairs (V8):
the event and failure counts of the whole and relock runs agree bit for
bit. The orientation of an exchange averages to zero within its error
in every run, so nothing chiral is selected. The exchange identity of
(a) holds event by event to $`10^{-13}`$.

On this evidence, V-SP8's clock rule is that a pair moves whole and
re-locks to its destination row. Re-locking is a mechanism step 22
introduced, not a consequence of the postulates, and Q4's tagged
crossings were not re-measured with collisions on.

---

## 5. Proposition V4: granular supply

The mean-field floor treats the sea as a continuum. It is not one: step 23
§1 made it a random configuration of aligned pairs with mean density
$`B`$, and an event happens within the row width $`dp`$ and the reach
$`y_{\max}`$, an aperture whose area is fixed by Q9(c),

```math
B \, dp \, 2y_{\max} = 1 .
```

At the physical density an emission's aperture holds **one** pair on
average.

**Proposition V4.** *If the sea is a Poisson configuration of mean
$`\beta`$ pairs per aperture, depleted locally as the mean-field ledger
says, then under (F) an emission is blocked with probability
$`e^{-n}`$, $`n = \beta - d`$ with $`d`$ the local mean deficit, and the
expected $`E`$ departs from the unfloored one by $`\varepsilon(\beta)
\to C e^{-\beta}`$ once $`\beta \gg \lambda^*`$. Measured, $`C`$ is
1.4 to 1.5 on the Eckart summit and 6.0 to 6.8 in the well, and at the
physical depth, $`\beta = 1`$, 37 to 38 per cent of emissions are blocked
and $`E`$ moves by 0.37 to 0.39 on the summit and 0.87 to 0.88 in the
well.*

*Verification* (Part D), the expected $`E`$ under (F) with the realised
fraction $`1 - e^{-n}`$:

| case | $`\beta`$ | blocked | $`\varepsilon`$ | $`\varepsilon e^{\beta}`$ | $`\Delta\langle p \rangle`$ | $`\Delta H`$ | $`\Delta M_-`$ | $`\Delta T`$ |
|---|---|---|---|---|---|---|---|---|
| Eckart, (S′) | 1 | 0.371 | 0.371 | 1.01 | −0.077 | $`+5\times10^{-5}`$ | −0.094 | −0.031 |
| | 2 | 0.139 | 0.168 | 1.24 | −0.029 | +0.004 | −0.042 | −0.013 |
| | 4 | 0.019 | 0.025 | 1.35 | −0.0036 | +0.0009 | −0.0062 | −0.0019 |
| | 8 | $`3.5\times10^{-4}`$ | $`4.6\times10^{-4}`$ | 1.37 | $`-6.5\times10^{-5}`$ | $`+1.8\times10^{-5}`$ | $`-1.1\times10^{-4}`$ | $`-3.5\times10^{-5}`$ |
| | 16 | $`1.2\times10^{-7}`$ | $`1.5\times10^{-7}`$ | 1.37 | $`-2\times10^{-8}`$ | $`+6\times10^{-9}`$ | $`-4\times10^{-8}`$ | $`-1\times10^{-8}`$ |
| PT well, (S′) | 1 | 0.369 | 0.884 | 2.40 | +0.044 | −0.32 | −0.32 | −0.030 |
| | 4 | 0.018 | 0.106 | 5.80 | +0.010 | −0.017 | −0.045 | −0.013 |
| | 8 | $`3.4\times10^{-4}`$ | $`2.0\times10^{-3}`$ | 5.97 | $`+1.8\times10^{-4}`$ | $`-3.1\times10^{-4}`$ | $`-8.4\times10^{-4}`$ | $`-2.5\times10^{-4}`$ |
| | 16 | $`1.1\times10^{-7}`$ | $`6.7\times10^{-7}`$ | 5.97 | $`+6\times10^{-8}`$ | $`-1\times10^{-7}`$ | $`-3\times10^{-7}`$ | $`-9\times10^{-8}`$ |

Under (S) the constants are 1.50 and 6.84 and the rows otherwise alike
(`sea_depletion_poisson.csv`). The blocked share is close to
$`e^{-\beta}`$ throughout (0.371 against 0.368, $`3.5\times10^{-4}`$
against $`3.4\times10^{-4}`$): the mean deficit is small where the
emissions are. Every diagnostic scales together, and none of them is the
norm.

So **the physical density cannot be a strict integer supply**: a sea of
one pair per aperture would make every anharmonic process depart from
quantum mechanics by tens of per cent. Three readings remain consistent.
The sea may be a ledger and not a supply (the first row of §1), in which
case it is invisible and its depth is unfixed. It may be a supply deep
enough to hide, $`\beta \gtrsim \ln(C/\delta)`$ for a measurement of
resolution $`\delta`$; for $`\delta = 10^{-3}`$ that is seven to nine
pairs per aperture, more for long-lived anharmonic dynamics, whose
$`\lambda^*`$ keeps growing (§4), and more as the reach opens. Or it may
be read in some other way this note has not found. A sea deeper than $`B`$ is therefore not merely
ontologically acceptable: if the sea is a supply it is required, and only
an experiment could bound its depth, and only from below.

This is the integer sea in mean field. The integer world ensemble itself,
at $`\nu = 1`$ and many seeds, is §5.1.

![The floor in mean field (left) and the integer sea (right) against the depth](https://raw.githubusercontent.com/billpage/wpmw/output/figures/sea_depletion_floor.png)

### 5.1 Proposition V6: the integer sea at the physical density

*Addendum (October 2026), open item V-SP2.* Proposition V4 put the
granularity in by hand, as Poisson emptiness on top of a mean field. Here
nothing is averaged before the end. Bodies and sea pairs are points; a
world at the physical density starts from one positon sampled from
$`W_0`$ ($`\nu = 1`$); each body fires at its per-parent rate; an event
is absorptive if a free negaton and a free positon sit in the apertures
about $`p \pm t\xi_q`$, emissive if a pair sits in the aperture about the
parent, and under (F) blocked otherwise; a positon and a negaton in one
mesh cell bind on contact (Q8). The sea is a Poisson configuration of
$`\beta`$ pairs per aperture, moving under (S) or (S′). The expectation
is an average over worlds, so only observables linear in $`E`$ mean
anything: the transmission $`T`$, $`\langle p \rangle`$ and
$`\langle H \rangle`$.

Each seed runs the world without the floor and the world with it. The
children of an emission are placed at the parent, as in the mesh
ledgers, so the bodies never depend on which pair was ionised, and the
pair is drawn from a separate random stream: the two worlds coincide
event for event until the first block, and their difference has a small
variance.

**Proposition V6 (measured).** *(a) Without the floor the integer worlds
reproduce the QLE in mean, up to the particle model's own baseline.
(b) Under (F) the share of emissions blocked falls as $`A e^{-\kappa\beta}`$
with $`\kappa \approx 0.55`$ to 0.6, not $`e^{-\beta}`$: the apertures
an emission finds hold nearly $`\beta`$ pairs on average, but are empty
far more often than a Poisson sea of that mean would be, increasingly so
with depth. (c) The observable departs as Proposition V4's mean field
predicts where it is resolved. (d) The cause is body inflation: at
$`\nu = 1`$ partners are scarce, a world inflates to tens of bodies, and
by V2(ii) every one of them is sea debt drawn near its parent.*

*Verification* (`src/demo_integer_sea.py`, Kaggle batch kernels from
commit `e916ac4` with that script added; Eckart summit to $`T = 8`$,
102 000 worlds under (S′) in three seed blocks and 34 000 under (S);
Pöschl–Teller well to $`T = 12`$, 8000 worlds; $`\Delta t = 0.02`$):

| case | $`\beta`$ | blocked | $`e^{-\beta}`$ | pairs found, mean | $`P(0)`$ over Poisson at that mean | $`\Delta T`$ | V4 mean field |
|---|---|---|---|---|---|---|---|
| Eckart, (S′) | 1 | 0.525 | 0.368 | 0.72 | 1.07 | $`-0.021 \pm 0.008`$ | $`-0.031`$ |
| | 2 | 0.301 | 0.135 | 1.47 | 1.30 | $`-0.011 \pm 0.008`$ | $`-0.013`$ |
| | 4 | 0.100 | 0.018 | 3.18 | 2.41 | $`-0.001 \pm 0.008`$ | $`-0.002`$ |
| | 8 | 0.0114 | 0.0003 | 7.02 | 12.7 | $`-0.001 \pm 0.004`$ | 0 |
| Eckart, (S) | 1 | 0.517 | 0.368 | 0.73 | 1.07 | $`-0.037 \pm 0.014`$ | $`-0.032`$ |
| | 2 | 0.295 | 0.135 | 1.48 | 1.30 | $`-0.027 \pm 0.014`$ | $`-0.014`$ |
| | 4 | 0.103 | 0.018 | 3.17 | 2.46 | $`-0.015 \pm 0.013`$ | $`-0.002`$ |
| | 8 | 0.0126 | 0.0003 | 6.93 | 12.9 | $`-0.015 \pm 0.008`$ | 0 |
| well, (S′) | 1 | 0.530 | 0.368 | 0.73 | 1.10 | $`-0.04 \pm 0.07`$ | |
| | 2 | 0.279 | 0.135 | 1.65 | 1.46 | $`-0.06 \pm 0.08`$ | |
| | 4 | 0.084 | 0.018 | 3.66 | 3.25 | $`-0.03 \pm 0.08`$ | |
| | 8 | 0.0084 | 0.0003 | 7.72 | 18.8 | $`+0.02 \pm 0.07`$ | |

Least squares on the logarithm give $`0.90\,e^{-0.55\beta}`$ on the
summit under (S′), $`0.86\,e^{-0.53\beta}`$ under (S) and
$`0.92\,e^{-0.59\beta}`$ in the well; the local slopes on the summit
under (S′) are 0.557, 0.549 and 0.544, so the law is clean over two
decades. The (S) block used the seeds of the first (S′) block, and its
unfloored worlds are the same worlds bit for bit, as Proposition Q1
requires of bodies that never read the sea; their floored departures
agree with that block's ($`-0.035`$, $`-0.013`$, $`-0.016`$, $`-0.008`$),
so the $`2\sigma`$ value at $`\beta = 8`$ under (S) is that seed
block's fluctuation, not the motion.
Without the floor the summit's worlds give $`T = 0.838 \pm 0.008`$,
$`\langle p \rangle = 1.204 \pm 0.031`$ and
$`\langle H \rangle = 1.258 \pm 0.038`$, against the mesh QLE's 0.869,
1.252 and 1.231. The momenta and energy agree; the transmission is
$`4\sigma`$ low, with no sea involved. Q-SP12 found the particle model's
classical baseline offset from the mesh's on this coarse momentum grid,
which is the likely cause; it is not checked here. The paired differences
share the offset and cancel it. A world on the
summit ends with 29.9 bodies (131.5 in the well) against a minimal
$`\lVert W \rVert_1`$ near 1.7, after 35.5 emissions, 9.2 absorptions
and 11.8 contact recombinations: $`f \approx 0.21`$.

So the mean-field estimate of §5 is right about the observable and wrong
about the mechanism. Blocking is worse than Poisson, because a parent
drains the aperture it keeps returning to, and under (S′) the pairs of
its own row travel with it (the frozen comb of Q-SP5); but most of the
emissions it blocks are inflation, whose $`\pm`$ deposits would have
cancelled again, and the departure of $`E`$ follows the mean field.
Keeping blocked emissions below a part in a thousand takes about twelve
pairs per aperture ($`\ln(900)/0.55`$), not seven to nine. How far the
observable's departure follows is resolved only to $`\beta = 2`$ at this
sample size; beyond it $`\lvert \Delta T \rvert < 0.016`$ (two standard
errors).

![Blocking against depth (left) and the transmission's paired departure against Proposition V4's mean field (right)](https://raw.githubusercontent.com/billpage/wpmw/output/figures/integer_sea_blocking.png)

### 5.2 Proposition V13: collisions in an integer sea

*Third addendum (October 2026), open item V-SP10.* §5.1's worlds, with
Definition (C) between point pairs. Two sea pairs within one aperture's
reach in position ($`\lvert\Delta x\rvert < y_{\max}`$) collide: a
defocus if they share a row, a focus if they are $`2q`$ rows apart. The
rates are $`2r_qB`$ per unordered pair for a defocus and $`r_qB`$ for a
focus. An aperture's position window at density $`s`$ holds $`s/B`$
partners (Q9(c)), so mass action holds in mean. At the kernel's rate
$`r_qB = \lvert K_q \rvert/2`$ at the pair's midpoint; uniform
collisions join neighbouring rows at $`rB = 1`$. A transition and its
reverse have equal rates, so a Poisson sea is invariant. Collisions draw
on a random stream of their own, so the bodies, the initial sea and the
unfloored worlds are §5.1's.

**Proposition V13 (measured).** *(a) Collisions alone keep a Poisson sea
Poisson. (b) In the floored worlds they reduce blocking, little at
$`\beta = 1`$ and more with depth: to about $`1.10\,e^{-0.74\beta}`$ at
the kernel's rate and $`1.06\,e^{-0.72\beta}`$ uniform, against
$`0.91\,e^{-0.55\beta}`$ without, and the local slope is still
steepening. (c) The column of pairs within reach of a parent, on every
row, is only 5 per cent short of $`\beta`$; the aperture on the parent's
row is 28 per cent short. The deficit lies between rows, which
collisions can move, and they refill it in part, at a rate comparable to
the drain. (d) The observable's departure is the same with and without
collisions, within the errors.*

*Verification* (`src/demo_integer_sea.py --collide`, Kaggle; Eckart
summit, (S′), $`\nu = 1`$, $`T = 8`$, $`\Delta t = 0.02`$; 12 000 worlds
at $`\beta = 1`$ and 2, 3000 at $`\beta = 4`$, 4000 at ten times the
kernel's rate). *Found* is the mean count in the aperture at a floored
emission, and *excess emptiness* is its $`P(0)`$ over Poisson at that
mean:

| collisions | $`\beta`$ | blocked | found | excess emptiness | collisions per world |
|---|---|---|---|---|---|
| none | 1 | $`0.524 \pm 0.002`$ | 0.715 | 1.071 | 0 |
| | 2 | $`0.300 \pm 0.001`$ | 1.468 | 1.301 | 0 |
| | 4 | $`0.100 \pm 0.002`$ | 3.188 | 2.414 | 0 |
| kernel rate | 1 | $`0.511 \pm 0.002`$ | 0.719 | 1.049 | 200 |
| | 2 | $`0.264 \pm 0.001`$ | 1.474 | 1.154 | 780 |
| | 4 | $`0.057 \pm 0.001`$ | 3.232 | 1.437 | 3023 |
| uniform, $`rB = 1`$ | 1 | $`0.503 \pm 0.002`$ | 0.726 | 1.040 | 2967 |
| | 2 | $`0.260 \pm 0.001`$ | 1.474 | 1.134 | 11 010 |
| | 4 | $`0.058 \pm 0.001`$ | 3.220 | 1.457 | 38 809 |
| kernel rate × 10 | 1 | $`0.477 \pm 0.003`$ | 0.767 | 1.028 | 1937 |
| | 2 | $`0.218 \pm 0.002`$ | 1.592 | 1.073 | 7359 |

The rows without collisions reproduce §5.1 (0.525, 0.301 and 0.100). The
local slopes of $`\ln`$(blocked) are 0.56 and 0.55 without collisions,
0.66 and 0.77 at the kernel's rate, and 0.66 and 0.75 uniform.
Collisions halve the excess emptiness of an aperture at $`\beta = 2`$
and $`4`$, and at ten times the kernel's rate bring it within 3 per cent
of Poisson at $`\beta = 1`$. The blocking there, 0.48, is then mostly
Poisson emptiness at a mean below $`\beta`$. Blocking fewer than one
emission in a thousand takes about 9.5 pairs per aperture on the fit, or
fewer if the slope keeps steepening, against about 12 without collisions
and 6.9 for a Poisson sea at the full depth. The observable's departures
with and without collisions agree within their errors: at the kernel's
rate $`\Delta T = -0.020`$, $`+0.009`$ and $`+0.001`$ at
$`\beta = 1`$, 2 and 4, against $`-0.014`$, $`-0.014`$ and $`-0.041`$
without, with standard errors 0.023 at $`\beta \le 2`$ and 0.046 at
$`\beta = 4`$.

(a) was checked with no bodies on twenty seas over 4000 random
apertures. At $`\beta = 1`$, after $`T = 8`$ at the kernel's rate, the
mean is 1.005, the variance 1.08 and $`P(0) = 0.379`$, against 0.366
for Poisson. At $`\beta = 4`$, $`P(0)`$ is 0.022 against 0.018, within
two standard errors. Uniform collisions give the same.

(c) is a column diagnostic, run locally on 1000 worlds. At a floored
emission it counts the pairs within reach of the parent in position on
every row, per aperture. That count is $`0.956\beta`$ at $`\beta = 1`$
and $`0.951\beta`$ at $`\beta = 2`$, while the aperture on the parent's
row holds $`0.72\beta`$. Collisions move pairs between rows at fixed
position, so they can bring an aperture up to its column, never above
it. Body inflation drains the parent's own row, which a parent keeps
returning to, and leaves the column nearly full.

*The contact sink.* The same local runs without the contact sink, 1000
worlds each, show what it does. A world ends with 54.4 bodies instead of
31.5. Absorption takes over part of the sink's work, with 28.3
absorptions a world instead of 9.9. Blocking at $`\beta = 1`$ rises from
0.523 to 0.576 (0.512 to 0.557 with collisions). In both cases the
sea's net debit per world is $`(N - 1)/2`$ exactly — 15.23 and 26.72 —
because each emission makes two bodies and debits one pair, and each
absorption or recombination removes two and credits one. So the contact
sink is not needed for the bookkeeping to close. It is needed to keep
inflation down, and with it blocking. Its scale — one mesh cell here —
is open, and a sink that binds on contact makes gray pairs unless the
clocks happen to agree (Q-SP10; V-SP12).

---

## 6. Proposition V5: two bodies

**Proposition V5.** *Two particles on a line, $`H = P^2/2M + p^2/2\mu +
V(r)`$, from a state $`W_0 = W_{\rm cm} \otimes W_{\rm rel}`$ with
$`W_{\rm cm} \ge 0`$, in the minimal ensemble. The joint sea's deficit
density at $`(X, P, r, p)`$ is $`W_{\rm cm}(X, P, t)`$ times the deficit
density of the one-dimensional relative problem at $`(r, p)`$. Hence the
joint debt equals the relative debt, $`D_{\rm joint} = D_{\rm rel}`$, and
the joint $`\lambda^*`$, measured against the four-dimensional density
$`B^2`$, is at most the relative one.*

*Proof.* By Proposition M1 every event moves relative momentum only, and
by V1 the centre of mass, which is free, makes none. The event rate, both
absorptive caps and the contact sink are each homogeneous of degree one in
the local populations, so on $`W_{\rm cm} \otimes W_{\rm rel}`$ every
allocation is $`W_{\rm cm}`$ times the relative one. Bodies and pairs of
the centre of mass move freely under either (S) or (S′), so deficits made
at $`(X, P)`$ are carried by the same flow as $`W_{\rm cm}`$. Pure states
obey $`\lvert W \rvert \le B`$ (step 16 §2), so $`W_{\rm cm}/B \le 1`$.
$`\square`$

Entanglement between the particles' coordinates is made here, often a lot
of it, and the sea pays for it exactly what the relative motion's
negativity costs — the Eckart summit run of §4 *is* such a collision, in
relative coordinates. What V5 does not cover is the genuinely
four-dimensional case: a non-quadratic external potential that couples the
centre of mass to the relative motion, or a preparation entangled across
them. There the demand on the sea is measured in §6.1, and step
21's body inflation, which grows as a product over degrees of freedom
(Corollary M10), would by V2(ii) be sea debt unless the ensemble is kept
minimal.

### 6.1 Proposition V7: two coupled bodies

*Addendum (October 2026), open item V-SP1.* Two particles on a line in
step 24's Pöschl–Teller well, repelling through a soft Gaussian,

```math
H = \tfrac12 p_1^2 + \tfrac12 p_2^2 + V(x_1) + V(x_2) + w_0\, e^{-(x_1 - x_2)^2/2s^2},
\qquad V(x) = -4\,{\rm sech}^2(x/2),
```

with $`w_0 = 2`$, $`s = 1`$, from a product of two minimal packets at
rest at $`x = \pm 2`$ ($`\sigma = 0.6`$). The well is not quadratic, so
the centre of mass and the relative motion couple and V5 does not apply.
The ledger is §3's minimal ledger in four dimensions: $`E(x_1, x_2, p_1,
p_2)`$ is transported spectrally, the sea has density $`B^2`$, and by
Proposition M1 the events fall into three families — $`p_1`$ jumps with
rates from $`V(x_1)`$, $`p_2`$ jumps from $`V(x_2)`$, and the pair jumps
$`(p_1, p_2) \to (p_1 \pm \xi, p_2 \mp \xi)`$ from the repulsion at
$`x_1 - x_2`$ — each with step 16's one-dimensional compensated kernel.
The mesh is the one-dimensional well's, $`dx = 0.3125`$ and
$`dp = 0.25`$, at $`48^4 = 5.3\times10^6`$ cells.

**Proposition V7.** *(a) Proposition V2 holds in four dimensions:
$`D = \Delta M_- - \Delta N_{\rm tr}/2`$, the same under (S) and (S′).
(b) Proposition V3 holds in four dimensions: under (F), $`E`$ is the
unfloored $`E`$ bit for bit if $`\beta \ge \lambda^*`$, measured against
$`B^2`$, and differs from it otherwise. (c) Measured: in every
realisation the coupled pair needs a deeper sea than one body in the
same well — $`\lambda^*`$ exceeds 1, the physical density, by $`t = 2`$
under (S) and emissive (S′) and by $`t = 4`$ under absorptive (S′) — and
it keeps climbing under (S′), as the one-dimensional long runs of §4 do.*

*Proof of (a) and (b).* The proofs of V2 and V3 use the event structure
(each event deposits $`\pm\tau`$ at two momenta of one parent), the
contact sink, and Proposition Q1. None of them refers to the dimension,
and M1 gives each of the three families exactly that structure.
$`\square`$

*Verification* (`src/demo_fourd_sea.py`, Kaggle batch kernels on a Tesla
T4 from commit `e916ac4` with that script added; $`T = 12`$,
$`\Delta t = 0.02`$; about 35 minutes a run). Proposition V2, with the
same columns as the table of §3:

| $`t`$ | $`D`$ | $`\Delta M_-`$ | $`\Delta N_{\rm tr}/2`$ | residual, (S) | (S′), absorptive | (S′), emissive |
|---|---|---|---|---|---|---|
| 2 | 0.2398 | 0.2496 | 0.0098 | $`1.6\times10^{-10}`$ | $`1.0\times10^{-10}`$ | $`-4\times10^{-11}`$ |
| 5 | 0.9208 | 1.7699 | 0.8491 | $`3.5\times10^{-10}`$ | $`2.0\times10^{-10}`$ | $`5\times10^{-11}`$ |
| 12 | 2.8072 | 6.5503 | 3.7430 | $`7.8\times10^{-10}`$ | $`4.2\times10^{-10}`$ | $`2.4\times10^{-10}`$ |

The first three columns are those of absorptive first, identical under
(S) and (S′); emission only gives $`D = 2.8173`$ at $`t = 12`$. The
residual is round-off over five million cells. At $`t = 12`$ the
transport's share of $`\Delta M_-`$ is 57 per cent, the same as the
one-dimensional well's on the same mesh (0.447 of 0.780, §3), so it is
the mesh's, and §3's refinement table says how it falls.

The floor, with $`\varepsilon`$ as in §4 and the one-dimensional well at
$`T = 12`$ for comparison:

| realisation | $`f`$ | $`\lambda^*`$ at $`t = 5`$ | at $`t = 12`$ | 1D well | $`\varepsilon`$ at $`\beta = 0.5`$ | at $`\beta = 1`$ | at $`0.98\,\lambda^*`$ | at $`1.02\,\lambda^*`$ |
|---|---|---|---|---|---|---|---|---|
| (S), absorptive first | 0.397 | 1.570 | 1.570 | 1.295 | 0.334 | 0.043 | $`1.6\times10^{-3}`$ | 0 |
| (S′), absorptive first | 0.397 | 1.312 | 1.859 | 0.765 | 0.131 | 0.031 | $`9.3\times10^{-4}`$ | 0 |
| (S′), emissive | 0 | 2.736 | 3.982 | 0.994 | 2.21 | 0.41 | $`1.2\times10^{-3}`$ | 0 |

The threshold is as sharp as in one dimension: at $`0.98\,\lambda^*`$ the
floor blocks $`1.6\times10^{-7}`$ of the emissions and moves $`E`$ by a
part in a thousand; at $`1.02\,\lambda^*`$ it blocks none. The (S′) runs
were repeated in a second kernel on another GPU session, which reproduced
every number to the last digit.

**How far to trust it.** The ledger's $`E`$ is the mesh QLE to a relative
$`L^2`$ error of 0.05 at $`t = 3`$ and 0.25 at $`t = 12`$, which is
the τ-leap's first-order error. It is larger than the one-dimensional well's
0.09 to 0.11 because three families act in each step. The observables
agree more closely. The purity of one particle's state,
$`{\rm Tr}\,\rho_1^2 = h \int W_1^2`$, is the mesh QLE's to 0.009 throughout, and
$`\langle x_1 \rangle`$ at $`t = 12`$ to 0.02. The mesh itself follows the exact
Schrödinger evolution, computed on the $`(x_1, x_2)`$ grid, only to
$`t = 5`$. The purity is 0.515 against 0.509 at $`t = 5`$. After that the
mesh departs by more than 5 per cent, and at $`t = 12`$ it gives 0.218
against 0.317: $`48^4`$ cells do not resolve the relative motion's fine
structure. So (a) and (b), which compare the ledger with itself, hold at
every $`t`$. The values of $`\lambda^*`$ are the mesh model's past
$`t = 5`$. Every crossing of 1 comes before $`t = 5`$, and so, in all three
realisations, does $`\lambda^*`$'s excess over the one-dimensional
value at $`t = 12`$.

The physical density is short in all three realisations. In one
dimension, absorptive (S′) was the realisation for which $`B`$ sufficed
(§4); here it does not. At $`\beta = 1`$ the floor moves $`E`$ by
$`\varepsilon = 0.031`$ in that realisation, and by 0.41 under emissive (S′). Under (S) $`\lambda^*`$ is
set by a single spike at $`t = 3`$, while under (S′) it keeps climbing.
That is the pattern of the one-dimensional long runs, where the deficits
are not sheared away and are not refilled (Q7). Here (S′) overtakes (S)
by $`t = 7`$, where in one dimension it took until $`t = 48`$ to 64;
that crossing lies past $`t = 5`$, so it is a property of the mesh model.
The reason for the extra depth is not separated here. A joint cell is
debited by three families at once, against a supply of $`B^2`$ per unit
of four-volume, but the families' shares of the worst cell were not
recorded (V-SP6). The answer to V-SP1 is
therefore: the four-dimensional sea obeys the same identities, its
floor has the same sharp threshold, and the coupled problem needs a
deeper sea than the one-body problem, under every realisation. V4's
case against the physical density as a strict supply is stronger in two
bodies than in one.

![The depth the floor needs for two coupled bodies against one body in the same well (left), and the purity of one particle's state from the ledger, the mesh QLE and Schrödinger (right)](https://raw.githubusercontent.com/billpage/wpmw/output/figures/fourd_sea.png)

### 6.2 Proposition V14: two coupled bodies with collisions

*Third addendum (October 2026).* §6.1's runs with Definition (C) in four
dimensions. Collisions act at one joint position $`(x_1, x_2)`$ along the
momentum directions of the three event families of Proposition M1:
$`(q, 0)`$ and $`(0, q)`$, with rates from $`V`$ at $`x_1`$ and $`x_2`$,
and $`(q, -q)`$, with rates from the repulsion at $`x_1 - x_2`$. Each
direction is a one-dimensional lattice of Proposition V10. At the
kernel's rate $`r B^2 = \lvert K \rvert/2`$ for the family's kernel at
the cell, 69 channels in all. Uniform collisions join neighbouring rows
in $`p_1`$ and in $`p_2`$ at $`rB^2 = 1`$. A cell in debt holds no pairs
to collide (§8.1).

**Proposition V14 (measured).** *(a) $`E`$ is that of the run without
collisions bit for bit (V8). (b) At the kernel's rate, $`\lambda^*`$ is
0.198 absorptive first and 0.540 emissive at $`t = 12`$, set by
$`t = 2`$ to 8, against 1.859 and 3.982 without. The physical density
suffices for two coupled bodies, as it does for one. (c) Uniform
collisions hold the absorptive-first depth at 0.342. The emissive depth
reaches 2.379 by $`t = 4`$, before they catch up, and the worst cell
then refills to $`0.49\,B`$ by $`t = 12`$.*

*Verification* (`src/demo_fourd_sea.py --gpu --motions`, Kaggle, Tesla
T4, $`48^4`$, $`T = 12`$, $`\Delta t = 0.02`$; 35 to 58 minutes a run):

| realisation | (S′) unamended, $`t = 5`$ | $`t = 12`$ | kernel rate, $`t = 5`$ | $`t = 12`$ | uniform, $`t = 5`$ | $`t = 12`$ |
|---|---|---|---|---|---|---|
| (S′), absorptive first | 1.312 | 1.859 | 0.187 | 0.198 | 0.342 | 0.342 |
| (S′), emissive | 2.736 | 3.982 | 0.540 | 0.540 | 2.379 | 2.379 |

The unamended columns reproduce §6.1 to the last digit. Before each run
a check on the GPU, over five steps and the channels with $`q \le 2`$,
found pairs and both momenta conserved to $`10^{-15}`$ per pair and the
collision entropy falling at every step; the uniform channels gave the
same. $`\max\lvert E - E_{\rm none}\rvert = 0`$ in all four runs with
collisions. Every value at $`t = 5`$ lies where the mesh follows
Schrödinger (§6.1).

So V7's case against the physical density as a strict supply is a case
against (S′) unamended. With collisions at the kernel's rate, two
coupled bodies need a fifth to a half of the physical density. A uniform
rate does not suffice: its emissive depth, 2.4, is set in the first four
time units, while the events outrun it, as in the one-dimensional well
at large reach (§4.4).

---

## 7. What an experiment could see

**The signature.** A sea that is a finite supply would show where
negativity is *made* fastest, as three things together: less negativity
than quantum mechanics predicts, a force that no potential exerts, and an
energy drift. The norm would be exact. This is not decoherence. Collapse
models and position-basis dephasing also suppress interference, but
conserve $`\langle p \rangle`$ and heat steadily (Bassi et al. 2013); the
floor pushes, in a direction set by which channels it starves.

**Where.** Not in a trap holding a cat (V1). In processes that make
negativity through non-quadratic dynamics: matter-wave diffraction, since
a grating is a periodic potential (molecules beyond 25 kDa, Fein et al.
2019); Kerr collapse and revival (Greiner et al. 2002; Kirchmair et al.
2013); tunnelling and barrier scattering; and, by V5, collisions.

**Why not yet a number.** The deviation $`\varepsilon(\beta)`$ is a
function of the model's regulators, and these have no laboratory values:
Q9(c) fixes the product of the row width and the reach, not each, and the
reach is imposed rather than derived (Q-SP6). Until it is, an
experimental resolution cannot be turned into a bound on $`\beta`$. What
can be said is relative: in exactly the regimes where Part D's
$`\beta = 1`$ departure is tens of per cent, experiment agrees with the
Schrödinger equation far better than that.

**What a null result settles.** That depletion and excess are not
observable at any depth an experiment can reach: excess because only a
rule that reads the level could see it, and such rules are excluded (§1);
depletion because a sea deep enough hides it, and a ledger hides it
altogether. Between "a ledger" and "a supply deeper than any experiment
resolves" no observation can choose. The simpler reading is the ledger.
That is the expected answer, and it is worth having: it says that the
sea's depth, like its phases (Q12) and its motion (Q1), is not something
the observable will ever fix, and that what the sea does that *is*
observable is exactly what the kernel does.

---

## 8. A demo defect: ringing booked as sea traffic

Step 16's ledger transports $`u_+`$ and $`u_-`$ separately and
spectrally. Each is the positive part of $`E`$, kinked wherever $`E`$
changes sign, so each rings, and its undershoots are negative densities.
The allocation $`A = \min(D, {\rm cap}_A, {\rm cap}_B)`$ then goes
negative where a partner field does: it "absorbs" a negative partner,
which adds a pair of bodies and debits the sea with no event behind it.
This is S-SP7's negative-cap leakage, diagnosed in J-SP2 of the
supplement, where its effect on $`f`$ was measured. Part F measures its
effect on the sea, which matters here because §§4 and 5 read the sea's
worst cell.

In a harmonic trap there are no events at all
($`\max\lvert K \rvert = 3\times10^{-13}`$). Over one period of the cat
of §2, step 16's ledger nevertheless books 0.503 negative absorptions,
the body count rises from 1.652 to 3.664, and the worst cell of the sea
falls to $`0.581\,B`$. The minimal ledger, on the same run, keeps the sea
at $`B`$ everywhere; its body count, $`\lVert E \rVert_1`$, moves by the
transport's 0.022.

On Theorem S8's row (Eckart summit, $`\Delta t = 0.01`$, $`T = 6`$)
step 16's ledger allocates 1.962 to absorption net, of which $`-1.137`$
is negative: the gross absorptions are 3.10, and more than a third of
them are run backwards. Against 2.946 emissions, that is more than a quarter of
all the sea's debits with no event behind them. The ledger reports
$`N(6) = 2.967`$, $`f = 0.400`$, and a worst cell of $`-0.339\,B`$ under
(S) and $`-0.060\,B`$ under (S′) — the figures Proposition Q3 quotes. The
minimal ledger on the same run gives $`N(6) = 1.699`$, $`f = 0.258`$, and
needs a depth of 1.101 under (S) and 0.664 under (S′).

$`E`$ is not affected: it changes by $`A + E_m = D`$ whatever $`A`$ is.
$`N`$, $`f`$ and the sea's profile are. The minimal ledger avoids the
defect altogether, since it transports only $`E`$; it also makes the
sea's debt the negativity made (V2), which is the quantity this note
needs. The population figures of steps 16 and 23 that read the sea's
profile were measured with the leakage (§9).

### 8.1 The event step, and collisions on cells in debt

*Third addendum (October 2026).* The long runs of §4.4 found two more
defects. Neither touches the tables of §§4 to 6.1 beyond what is said
here.

*The event step.* The minimal ledger advances $`E`$'s event term by one
explicit step of length $`\Delta t`$ (`SeaLedger.channels_supply`): a
forward Euler step, made nonlinear only by the species split and the
contact sink. On an oscillatory symbol that amplifies each mode
by $`\sqrt{1 + \omega^2\Delta t^2}`$ per step, a growth rate of
$`\omega^2\Delta t/2`$, with $`\omega = \max\lvert{\rm Im}\,{\rm sym}_e\rvert`$
the event symbol's largest oscillation rate. $`\omega`$ grows with the
reach: on the summit it is 1.00, 2.49, 5.06 and 10.1 at $`\pi`$ to
$`8\pi`$, and in the well 1.03, 4.01, 9.97 and 20.2. The growth bites
where a state stays in the event region. In the well at $`4\pi`$ under
(S′) unamended, $`\lambda^*`$ reaches 3841 by $`t = 96`$ at
$`\Delta t = 0.02`$, and the global debt is 3.55, 0.992 and 0.666 at
$`\Delta t = 0.01`$, 0.005 and 0.0025. On the summit the packet leaves,
and the error is mild. §4.3's $`4\pi`$ row ran at $`\Delta t = 0.02`$;
at 0.01 its emissive depth with collisions is 0.426, not 0.446, and its
absorptive depth 0.321, not 0.318; without collisions the emissive depth
at $`T = 8`$ is 1.857, not 1.754. A higher-order or exponential step
for the event term, or a step small enough that $`\omega^2\Delta t\,T/2`$ stays well
below one, would cure it; neither is done here (V-SP11).

*Collisions on cells in debt.* Without a floor the sea's density goes
negative where the debt is deepest. Mass action read literally on a
negative density is not a collision: $`s_n^2`$ stays positive, and the
flux it drives runs away. Collision runs at $`4\pi`$ diverged that way.
The collision term now uses $`\max(s, 0)`$, so a cell in debt holds no
pairs to collide, and it takes substeps for the largest rate the current
sea allows. Nothing changes while $`s > 0`$: Part G's main table
reproduces 0.251 and 0.468 with the guard. The Part H runs that diverged
before the guard was added are not used.

---

## 9. Corrections

- **Step 16** ([`sea_population_equilibrium.md`](sea_population_equilibrium.md)).
  The worst-cell figures of Theorems S4 and S8, and the body counts of S8
  and S9, include sea debits booked from transport ringing (§8); results
  about $`E`$ stand. Not re-measured here; S-SP7 carries it.
- **Step 23** ([`force_blind_sea.md`](force_blind_sea.md)). Q-SP2 is
  answered for the minimal ledger by Proposition V3 and Part E: supply
  becomes limiting under either motion once the sea is shallower than
  $`\lambda^*`$, and in a bound anharmonic system $`\lambda^*`$ keeps
  growing with time, and after the first few periods faster under (S′)
  than under (S). (The first version of this note had the motions the
  wrong way round here; §4's table was right.)
  Proposition Q3's tables and Part F's sea profiles share §8's
  leakage. §9.4(g)'s "depth is no observable" holds for every rule that
  does not read the sea; under the floor, depth is observable below
  $`\lambda^*`$ in mean field and below about $`\ln(C/\delta)`$ for an
  integer sea (Propositions V3 and V4).
- **Step 23, under the amended (S′)** (§4.1). Proposition Q3's "deficits
  stay in their rows" and Proposition Q7's "shallow but not refilled"
  describe (S′) as first stated; with pair collisions deficits relax
  (Proposition V10). Q-SP2's answer above holds for (S′) unamended;
  amended, the depth the floor needs stops growing in every case run.
  Proposition Q4's tagged crossings and the readings of §§8 and 9 ran with
  pairs that keep their rows, and are not re-measured (V-SP8).
- **This note, §5.** Proposition V4's $`e^{-\beta}`$ is the share of
  emissions blocked only for a sea that is Poisson at every event. The
  integer worlds block $`0.9\,e^{-0.55\beta}`$ (Proposition V6), while
  the observable follows V4's mean field where resolved.
- **This note, §4.3** (third addendum). The reach table's $`4\pi`$ row ran
  at $`\Delta t = 0.02`$, where the event step's error is about 5 per cent
  (§8.1): at 0.01 the emissive depth with collisions is 0.426, not 0.446.
  The conclusions stand. "In every case run" in §0 and in V10(d) holds at
  the default reach; at $`4\pi`$ in the well uniform collisions fail
  (Proposition V11).
- **This note, §6.1** (third addendum). "The physical density is short in
  all three realisations" holds for (S′) unamended; with collisions at the
  kernel's rate it suffices (Proposition V14).

---

## 10. Open items

- **V-SP1.** *Answered* by Proposition V7 (§6.1): two coupled bodies
  obey V2 and V3, and need a deeper sea than one, in every realisation.
- **V-SP2.** *Answered* by Proposition V6 (§5.1): blocking falls as
  $`0.9\,e^{-0.55\beta}`$, not $`e^{-\beta}`$, because of body
  inflation; the observable follows V4's mean field where resolved.
- **V-SP3.** The regulators' laboratory values. Q9(c) fixes the aperture's
  area; until the reach is derived (Q-SP6), $`\varepsilon(\beta)`$ cannot
  be compared with an experimental resolution.
- **V-SP4.** Which realisation is ontological? Absorptive-first and
  emissive realisations change $`E`$ identically but place the sea's
  credit differently (Q8), and so change $`\lambda^*`$ (Part C). A floor
  makes the choice observable in principle.
- **V-SP5.** Re-measure the sea-profile figures of steps 16 and 23 (S4,
  S8, S9, Q3, Q7) without the leakage of §8, in the minimal ledger or
  with J-SP2's repair.
- **V-SP6.** The four-dimensional ledger past $`t = 5`$: a finer mesh
  (the $`48^4`$ run is a T4's limit at about 35 minutes a run), and the
  three families' shares of the worst cell's debit, which would say why
  two bodies need more sea than one.
- **V-SP7.** The integer worlds' transmission on the Eckart summit is
  $`4\sigma`$ below the mesh QLE's with no floor (§5.1). Check it against
  Q-SP12's classical-baseline offset by running the worlds and
  the mesh with events switched off.
- **V-SP8.** The clock at a collision. Under (S′) a pair's clock is
  locked to its row's plane wave (Proposition Q2), and Q5's reading is
  taken from those clocks; a collision moves two pairs to new rows.
  Specify the clock rule, measure the misalignment collisions bring
  against Q6's resets, and re-measure §§8 and 9 of step 23 and Q4's tagged
  crossings with collisions on. *Answered in part* by Proposition V12
  (§4.5): in step 23's particle model only whole pairs re-locked to their
  new row keep the sea dark, the lock and the live reading. Re-locking is
  not derived from the postulates, and Q4's tagged crossings were not
  re-measured.
- **V-SP9.** The collision rate. The bound holds at the kernel's rate down
  to about a tenth of it, and for uniform collisions between neighbouring
  rows down to $`rB = 0.25`$ (§4.3). Which rate law is the ontology's, and
  what fixes its constant (compare Q-SP9)? Does the emissive depth stop
  growing with the reach in long runs, and does the deviation under
  uniform collisions settle? *Answered in part* by Proposition V11
  (§4.4): at the kernel's rate the emissive depth stops growing with the
  reach on the summit (0.50 at $`8\pi`$ against 17 without), and a rate
  that does not keep pace with the events fails in the well at $`4\pi`$.
  The emissive deviation keeps creeping under either law. Which law is
  the ontology's, and what fixes its constant, is still open.
- **V-SP10.** Collisions in an integer sea. Proposition V6's blocking at
  $`\nu = 1`$ comes from body inflation draining the apertures a parent
  returns to. Do collisions between point pairs refill them, and how does
  blocking then fall with $`\beta`$? *Answered* by Proposition V13
  (§5.2): collisions keep a Poisson sea Poisson and refill a parent's row
  in part, and blocking falls as about $`1.1\,e^{-0.74\beta}`$.
- **V-SP11.** The well beyond $`2\pi`$. The ledger's explicit event step
  grows errors at the rate $`\omega^2\Delta t/2`$ (§8.1), and at $`4\pi`$
  the debt still moves with $`\Delta t`$ at 0.0025. A higher-order or
  exponential event step, or $`\Delta t \lesssim 10^{-3}`$ with a
  convergence check, then $`8\pi`$; and whether the emissive depth at
  $`4\pi`$ stays at 0.79 under the kernel's rate.
- **V-SP12.** The contact sink's scale. It binds a positon and a negaton
  in one mesh cell. Without it a world inflates to twice the bodies and
  blocks more (§5.2). Would a sink at the aperture's scale (Q9(c)) keep the
  worlds' mean and cut blocking? A sink that binds on contact makes gray
  pairs unless the clocks agree (Q-SP10), so how should its products'
  clocks be set?
- **V-SP13.** Aligned pairs as readers. Proposition V12(c) and §7 of
  [`chord_reading.md`](../supplement/chord_reading.md) propose that the
  sea reads the kernel itself: an aligned pair integrates its chords and
  fires (Q11's bank, Q-SP14), moving itself and a partner. Build it in the
  particle model. Which frame does a sea reader integrate, $`U`$ or
  $`U_{\rm res}`$? Does it fire at the kernel's rate? Does the critical
  sampling of Q-SP13 allow it at $`\nu = 1`$?

---

## 11. Numerical verification

`src/demo_sea_depletion.py`, Parts A to G, on step 16's mesh
($`n_r = 128`$, $`r \in [-20, 20)`$, $`n_p = 64`$, $`dp = 0.25`$, reach
$`2\pi`$; Part E's reach ladder at $`dp = 0.125`$). Parts C, D and E
write CSV files and the floor and long-run figures are drawn from them.
About forty minutes on one core for Parts A to F, and fifteen more for Part
G at its default length; the parts run independently.

Figures (on the `output` branch): `sea_depletion_free_negativity.png`
(Part A), `sea_depletion_floor.png` (Parts C and D),
`sea_depletion_long.png` (Part E, drawn from the longest run present; the
one shown is `--parts E --t-long 96 --no-reach`). CSV files:
`sea_depletion_floor.csv`, `sea_depletion_poisson.csv`,
`sea_depletion_long.csv`, `sea_depletion_long_T96.csv`,
`sea_depletion_reach.csv`. The transport check of §3 is Part R, which is not
in the default set (`--parts R`).

Part G (§§4.1 to 4.3, second addendum) splits into sub-parts with
`--g-subs`: I, the SymPy identities and the one-row invariants (seconds);
M, the main comparison; S, the rate scan; R, the other cases; X, the reach
ladder. The tables were made with `--parts G --t-long 96 --g-subs IMS` and
`--g-subs RX` in two processes, 23 minutes on two cores, at commit
`15546a8` plus this change. CSV files: `sea_collisions_main_T96.csv`,
`sea_collisions_scan_T96.csv`, `sea_collisions_cases_T96.csv`,
`sea_collisions_reach.csv`. Figure: `sea_collisions.png`, drawn from the
main and scan files of the longest run present.

§5.1 and §6.1 ran as private Kaggle batch kernels (`wpmwlib.kaggle_batch`)
from commit `e916ac4` with the two scripts added. They are too long for one
core.

- `src/demo_integer_sea.py` (§5.1): one process per core. Five kernels wrote `integer_sea_<case>_<motion>_nu1*.csv`.
  `--summarise` combines them, weighted by worlds, and draws
  `integer_sea_blocking.png`.
- `src/demo_fourd_sea.py --gpu` (§6.1): CuPy on a Tesla T4, 8.6 hours for
  the three realisations with their floor depths. `--small` runs a
  $`24^4`$ mesh on one core in minutes, for testing. The script writes
  `fourd_sea.csv`, `fourd_sea_trace.csv` and `fourd_sea_schroedinger.csv`
  after every run. `--plot` draws `fourd_sea.png` from the last two.

The third addendum (§§4.4, 4.5, 5.2, 6.2) ran as private Kaggle batch
kernels from commit `7bf4df1` with this change applied as it stood at
the time, except the local runs named below; the guard of §8.1 was added
after the first Part H kernel, whose diverged runs are not used. Every
table gives the settings.

- §4.4: `src/demo_sea_depletion.py --parts H` with `--h-grid`,
  `--h-reach`, `--h-cases`, `--h-motions` and `--h-dt`; one process per
  run, `--workers`. Seven kernels: the coarse ladder at
  $`\Delta t = 0.02`$, a $`\Delta t = 0.01`$ check at $`4\pi`$, the fine
  grid on the summit and in the well at 0.01, the well at $`4\pi`$ at
  0.005 and 0.0025, and the summit at $`8\pi`$ at 0.0025. A run at
  $`4\pi`$ in the well takes up to three hours at $`\Delta t = 0.0025`$.
  CSV files: `sea_collisions_reachlong_<grid><tag>.csv`, one row per
  sample time.
- §4.5: `src/demo_pair_collision_clocks.py`. Parts A and B run in about a
  minute on one core. Part C ran as three kernels, one per $`\kappa`$,
  with `--parts ABC --nu 32 --kappa K --coll none whole exchange matched
  relock --seeds 4 --shadow --t-end 14`, about twenty minutes each on four
  cores. CSV files: `pair_collision_clocks_C<tag>.csv` and
  `pair_collision_clocks_trace<tag>.csv`.
- §5.2: `src/demo_integer_sea.py --collide none kernel:1 uniform:1` in
  three kernels: $`\beta = 1, 2`$ (12 000 seeds, 5.5 hours on four
  cores), $`\beta = 4`$ (3000) and `kernel:10` at $`\beta = 1, 2`$
  (4000). `--check-poisson` is (a), about two minutes. The column
  diagnostic and the runs with and without the contact sink (`--no-gc`)
  are local, 1000 seeds each, about half an hour on one core.
- §6.2: `src/demo_fourd_sea.py --gpu --motions none collide-k:1
  collide-u:1`, one kernel per realisation, about two hours each. The
  script writes `fourd_sea_collisions<tag>.csv` and
  `fourd_sea_collisions_trace<tag>.csv`.

---

## 12. Sources

- Ahmadiniaz, N. et al. (BIREF@HIBEF collaboration). *Towards a vacuum
  birefringence experiment at the Helmholtz International Beamline for
  Extreme Fields*, High Power Laser Sci. Eng. **13** (2025) e7;
  [arXiv:2405.18063](https://arxiv.org/abs/2405.18063).
- Bassi, A., Lochan, K., Satin, S., Singh, T. P. and Ulbricht, H.
  *Models of wave-function collapse, underlying theories, and experimental
  tests*, Rev. Mod. Phys. **85** (2013) 471–527.
- Boltzmann, L. *Weitere Studien über das Wärmegleichgewicht unter
  Gasmolekülen*, Sitzungsber. Kais. Akad. Wiss. Wien, Math.-Naturwiss. Cl.
  **66** (1872) 275–370.
- Fein, Y. Y., Geyer, P., Zwick, P., Kiałka, F., Pedalino, S., Mayor, M.,
  Gerlich, S. and Arndt, M. *Quantum superposition of molecules beyond
  25 kDa*, Nat. Phys. **15** (2019) 1242–1245.
- Greiner, M., Mandel, O., Hänsch, T. W. and Bloch, I. *Collapse and
  revival of the matter wave field of a Bose–Einstein condensate*, Nature
  **419** (2002) 51–54.
- Horn, F. and Jackson, R. *General mass action kinetics*, Arch. Rational
  Mech. Anal. **47** (1972) 81–116.
- Kac, M. and van Moerbeke, P. *On an explicitly soluble system of
  nonlinear differential equations related to certain Toda lattices*, Adv.
  Math. **16** (1975) 160–169.
- Kenfack, A. and Życzkowski, K. *Negativity of the Wigner function as an
  indicator of non-classicality*, J. Opt. B: Quantum Semiclass. Opt. **6**
  (2004) 396–404.
- Kirchmair, G., Vlastakis, B., Leghtas, Z., Nigg, S. E., Paik, H.,
  Ginossar, E., Mirrahimi, M., Frunzio, L., Girvin, S. M. and
  Schoelkopf, R. J. *Observation of quantum state collapse and revival due
  to the single-photon Kerr effect*, Nature **495** (2013) 205–209.
- Pawula, R. F. *Approximation of the linear Boltzmann equation by the
  Fokker–Planck equation*, Phys. Rev. **162** (1967) 186–188.
- Vlastakis, B., Kirchmair, G., Leghtas, Z., Nigg, S. E., Frunzio, L.,
  Girvin, S. M., Mirrahimi, M., Devoret, M. H. and Schoelkopf, R. J.
  *Deterministically encoding quantum information using 100-photon
  Schrödinger cat states*, Science **342** (2013) 607–610.
- Volterra, V. *Leçons sur la théorie mathématique de la lutte pour la
  vie*, Gauthier-Villars, Paris (1931).
