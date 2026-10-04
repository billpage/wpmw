# Depletion and excess: when the sea could be seen

> Empty space is not empty here: a cell that holds no `W` still holds aligned pairs at density `B`. Could the sea's depletion or excess ever be seen? By Proposition Q1 the observable never reads the sea, so only a rule that reads it could show either; a rule that reads its level, Theorem S5's throttle, is already excluded, so excess is unobservable. **Proposition V1**: admissible negativity is free. For a potential of degree two or less the compensated kernel vanishes identically (Theorem G2), so a cat state held in a harmonic trap makes no event and touches no pair; the sea pays only for the negativity the dynamics makes, where `V''' != 0`, and in an anharmonic well a cat makes events 1.4 times as fast as one packet because its fringes are parents too. **Proposition V2**: from a minimal preparation the sea's debt is at least the growth of the negative mass of `W`, with equality exactly when the ensemble stays minimal; body inflation is sea debt, so a finite sea would make garbage collection a requirement. **Proposition V3**: under a floor — no pair, no ionisation, Proposition S0 read literally — `E` is unchanged bit for bit while the sea is at least as deep as the depth `lambda*` the unfloored run needs, and departs below it through a spurious force, an energy drift and a change of negativity, never through the norm. On the Eckart summit `lambda*` is 1.10 under (S) and 0.65 under (S′); in an anharmonic well it climbs without levelling off, to 2.2 and 2.9 by `t = 96`, and it grows with the reach as Theorem G5's census does, which answers Q-SP2. **Proposition V4**: the sea is granular, and the event aperture of Q9(c) holds one pair on average at the physical density, so an integer sea of `beta` pairs per aperture blocks close to `e^-beta` of emissions and moves `E` by `C e^-beta`, `C` from 1.4 to 6.8 — by 37 to 88 per cent at `beta = 1`. The physical density cannot be a strict supply; a sea that is one must be seven to nine pairs per aperture deep to hide below a part in a thousand, and only an experiment could bound that depth, from below. **Proposition V5**: entanglement made by a pair potential costs the sea exactly the negativity of the relative motion. Where to look is where negativity is made fastest — diffraction, Kerr revivals, tunnelling, collisions — but the model's regulators have no laboratory values, so no bound follows yet; a null result is expected, and says the sea's depth, like its motion and its phases, is not fixed by the observable. Step 16's ledger books transport ringing as sea traffic (S-SP7), over a quarter of the sea's debits on Theorem S8's run; the minimal ledger used here does not.

*Ladder abstract — see the [full list](README.md#the-ladder).*

**Status.** Analysis note, step 24 of the ladder. Companion demo:
`src/demo_sea_depletion.py` (Parts A to F). Nothing here depends on
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
- **Proposition V5.** Two bodies. Entanglement made by a pair potential
  costs the sea exactly the negativity the relative motion makes, and per
  joint cell never more than the one-dimensional problem (§6).
- What an experiment could and could not see (§7), and a quantification of
  a known demo defect — step 16's ledger books transport ringing as sea
  traffic — that the minimal ledger used here avoids (§8).

**Leaves open** the items of §10, chief among them the integer world
ensemble at $`\nu = 1`$ (V-SP2) and a four-dimensional sea ledger (V-SP1).

**Inherits.** Theorems S0, S5 and S7 from step 16; Theorems G2 and G5
from step 17; Proposition M1 from step 21; Propositions Q1, Q8 and Q9 from step 23, and step 23's picture of the
sea as a random configuration of aligned pairs (§1 there).

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
at $`\nu = 1`$ and many seeds, is open item V-SP2.

![The floor in mean field (left) and the integer sea (right) against the depth](https://raw.githubusercontent.com/billpage/wpmw/output/figures/sea_depletion_floor.png)

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
them. There the demand on the sea has not been measured (V-SP1), and step
21's body inflation, which grows as a product over degrees of freedom
(Corollary M10), would by V2(ii) be sea debt unless the ensemble is kept
minimal.

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
  growing with time, faster under (S) than under (S′).
  Proposition Q3's tables and Part F's sea profiles share §8's
  leakage. §9.4(g)'s "depth is no observable" holds for every rule that
  does not read the sea; under the floor, depth is observable below
  $`\lambda^*`$ in mean field and below about $`\ln(C/\delta)`$ for an
  integer sea (Propositions V3 and V4).

---

## 10. Open items

- **V-SP1.** A four-dimensional sea ledger: two particles in a
  non-quadratic external potential, where the centre of mass and the
  relative motion couple and Proposition V5 does not apply. Measure
  $`\lambda^*`$ against $`B^2`$, in a minimal ensemble, on the grids of
  step 21 (a Kaggle-sized run).
- **V-SP2.** The integer world ensemble at $`\nu = 1`$: bodies and sea
  pairs as points, emission by ionising a pair within the aperture of
  Q9(c), and many seeds, to check Proposition V4's mean-field
  $`e^{-\beta}`$ and include the depletion of an aperture by its own
  events.
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

---

## 11. Numerical verification

`src/demo_sea_depletion.py`, Parts A to F, on step 16's mesh
($`n_r = 128`$, $`r \in [-20, 20)`$, $`n_p = 64`$, $`dp = 0.25`$, reach
$`2\pi`$; Part E's reach ladder at $`dp = 0.125`$). Parts C, D and E
write CSV files and the floor and long-run figures are drawn from them.
About forty minutes on one core; the parts run independently.

Figures (on the `output` branch): `sea_depletion_free_negativity.png`
(Part A), `sea_depletion_floor.png` (Parts C and D),
`sea_depletion_long.png` (Part E, drawn from the longest run present; the
one shown is `--parts E --t-long 96 --no-reach`). CSV files:
`sea_depletion_floor.csv`, `sea_depletion_poisson.csv`,
`sea_depletion_long.csv`, `sea_depletion_long_T96.csv`,
`sea_depletion_reach.csv`. The transport check of §3 is Part R, which is not
in the default set (`--parts R`).

---

## 12. Sources

- Ahmadiniaz, N. et al. (BIREF@HIBEF collaboration). *Towards a vacuum
  birefringence experiment at the Helmholtz International Beamline for
  Extreme Fields*, High Power Laser Sci. Eng. **13** (2025) e7;
  [arXiv:2405.18063](https://arxiv.org/abs/2405.18063).
- Bassi, A., Lochan, K., Satin, S., Singh, T. P. and Ulbricht, H.
  *Models of wave-function collapse, underlying theories, and experimental
  tests*, Rev. Mod. Phys. **85** (2013) 471–527.
- Fein, Y. Y., Geyer, P., Zwick, P., Kiałka, F., Pedalino, S., Mayor, M.,
  Gerlich, S. and Arndt, M. *Quantum superposition of molecules beyond
  25 kDa*, Nat. Phys. **15** (2019) 1242–1245.
- Greiner, M., Mandel, O., Hänsch, T. W. and Bloch, I. *Collapse and
  revival of the matter wave field of a Bose–Einstein condensate*, Nature
  **419** (2002) 51–54.
- Kenfack, A. and Życzkowski, K. *Negativity of the Wigner function as an
  indicator of non-classicality*, J. Opt. B: Quantum Semiclass. Opt. **6**
  (2004) 396–404.
- Kirchmair, G., Vlastakis, B., Leghtas, Z., Nigg, S. E., Paik, H.,
  Ginossar, E., Mirrahimi, M., Frunzio, L., Girvin, S. M. and
  Schoelkopf, R. J. *Observation of quantum state collapse and revival due
  to the single-photon Kerr effect*, Nature **495** (2013) 205–209.
- Vlastakis, B., Kirchmair, G., Leghtas, Z., Nigg, S. E., Frunzio, L.,
  Girvin, S. M., Mirrahimi, M., Devoret, M. H. and Schoelkopf, R. J.
  *Deterministically encoding quantum information using 100-photon
  Schrödinger cat states*, Science **342** (2013) 607–610.
