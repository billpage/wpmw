# Emission and Absorption: the Two Ways to Settle One Event

**A tutorial introduction to the mechanism by which the compensated Liouville
algorithm keeps its sea and its body population in balance. Self-contained: it
assumes no prior reading in this project, and builds the signed ensemble, the
sea, and the event from scratch.**

---

## 0. What this note is

Section 5.3 of
[`docs/algorithm/compensated_liouville_algorithm.md`](../algorithm/compensated_liouville_algorithm.md)
gives a five-line table naming two realisations of a residual-channel event,
and two more that are excluded. That table is normative and correct and almost
unreadable on first contact. The mechanism it describes has since become
central — it is what makes the sea population an attractor rather than a tuned
constant — and a reader meeting it for the first time deserves better than a
table.

This note is that reader's on-ramp. It is expository, not normative. Where it
and the specification disagree, the specification wins.

Companion code: [`src/demo_emission_and_absorption.py`](../../src/demo_emission_and_absorption.py),
parts A–G. Every number quoted below is an output of that script.

**Two things this note contributes** beyond exposition, both numerical and both
new:

- **Theorem J1/J2.** The quantity $`P = S + N/2`$ — bound pairs plus half the
  free bodies — is conserved *exactly* by the event channel, at any absorptive
  fraction, but only **globally**. Cell by cell it is violated by nine orders
  of magnitude more. §5 works out what each half of that means.
Two open items are logged in §12 rather than argued in the body, because
neither belongs in a tutorial: **J-SP1** withdraws support for one existing
explanation of the conjunction cost, and **J-SP2** proposes a change to the
allocation loop of the specification. J-SP2 in particular is a separate
subject and is written as a proposal, not a finding about the ontology.

**A note on metaphor.** An earlier draft of this note leaned on bookkeeping
language — settling a debt, paying in one currency or another. That language is
ontologically wrong here. If the sea is real, the difference between the two
realisations is not a difference in how something is *recorded*; different
things happen to different real objects. The vocabulary below is accordingly
demographic: three populations, births and deaths among them, and an
equilibrium between two opposed processes. The word *ledger* is kept, because
the specification uses it, but it means an inventory of real populations and
not a choice of description.

---

## 1. Worlds, and why the ensemble must be signed

Start with the picture the project is committed to. A quantum state is not a
wave; it is a large collection of **worlds**, each a single classical
configuration, each a point $`(x, p)`$ in phase space. There is no wavefunction
in the ontology. What we call the quantum state is a statement about how the
worlds are distributed, and what we call quantum behaviour is what happens
because worlds interact with each other.

The distribution we want them to reproduce is the Wigner function $`W(x,p)`$.
And immediately there is a problem: $`W`$ goes negative. It is not a
probability density, and no collection of counted objects can have a negative
count.

The project's response is to let the ensemble carry a **sign**. There are two
species of world-particle:

- **positons**, carrying $`+1`$,
- **negatons**, carrying $`-1`$.

Write $`u^{+}(x,p)`$ and $`u^{-}(x,p)`$ for their densities. Then

```math
E \;=\; u^{+} - u^{-}, \qquad N \;=\; u^{+} + u^{-} ,
```

and the identification is $`E = W`$. The observable is the **difference** of two
populations. The sum $`N`$ — how many bodies are actually present — is not
constrained by quantum mechanics at all. Two ensembles with wildly different
$`N`$ can represent the very same $`W`$, so long as their difference matches.

That single fact is the seed of everything in this note. Because $`E`$ is a
difference, there is always more than one way to change it, and the alternatives
are invisible to the observable but very visible to the bookkeeping.

**Minimal and padded ensembles.** The smallest ensemble representing a given
$`W`$ puts $`u^{+} = \max(W, 0)`$ and $`u^{-} = \max(-W, 0)`$, so that in each
cell only one species is present and $`N = |E|`$. Any number of coincident
$`+/-`$ pairs can then be added without touching $`E`$. The ratio

```math
\rho \;=\; N_0 \big/ \textstyle\int |E|
```

measures how padded a preparation is; $`\rho = 1`$ is minimal.

---

## 2. The sea

There is a third population, and it is the one that makes the mechanism work.

Consider a positon and a negaton occupying the *same* cell $`(x,p)`$, bound
together. Its total charge is $`+1 - 1 = 0`$. It contributes nothing to $`E`$,
so it is invisible in every observable — it is **dark**. Call it a **neutral
bound pair**, and call the standing population of such pairs the **sea**, with
density $`S(x,p)`$ and a uniform background level

```math
B \;=\; 2/h ,
```

two bound pairs per Planck cell. (Why exactly $`2/h`$ is a separate story told
in [`species_sectors_and_annihilation.md`](../analysis/species_sectors_and_annihilation.md);
for this note only its existence matters.)

It is worth being careful here, because there are *two* dark objects in this
project and they do opposite jobs:

| object | contributes to $`E`$? | live in the dynamics? |
|---|---|---|
| neutral bound pair (the sea) | no | **yes** — it is the reservoir events draw on |
| crystal shift $`c\cdot\mathbb{1}`$, $`c = 2`$ | yes, but inertly | no — provably inert |

The sea is invisible *and* active. That combination is what lets it absorb
bookkeeping that the observable must not see.

**A bound pair is not two world-particles.** This is the point most likely to
trip a careful reader later, so state it now. A bound sea pair is a single
neutral object. Its constituents do not have separate worldlines, separate
momenta, or separate identities while bound. They acquire them at the moment
the pair is broken, and they lose them at the moment two free bodies are bound.
Birth and death, not motion. §9.1 returns to why this matters.

---

## 3. What an event is

Worlds move classically between events: $`\dot{x} = p/\mu`$, $`\dot{p} = -V'(x)`$.
Everything non-classical happens at events.

An event is a **momentum transfer between two worlds**, and the compensated
splitting of the algorithm arranges that the classical force has already been
removed, so the events carry only the residual — the genuinely quantum part.
The rate field $`K_{\mathrm{res}}(x, \xi_q)`$ that governs them has three
properties, all verified in part A of the demo:

```
   max |Im K_res| / max |K_res|        6.872e-17     it is real
   max |K_q + K_(-q)|   (oddness)      1.110e-16     it is odd in q
   sum_q K_res          (worlds)       2.063e-16     it moves no worlds
   sum_q xi_q K_res     (momentum)     1.857e-15     it moves no momentum
```

The oddness is the important one. The channel at $`+\xi_q`$ and the channel at
$`-\xi_q`$ carry exactly equal and opposite weight. On the barrier flank, at
$`r = -1.56`$ and $`q = 3`$:

```
   K = -0.13064,   K_(-q) = +0.13064
```

These are **not two events**. They are the two **legs of one event**. An event
picks a parent momentum row $`p`$ and deposits

- $`+1`$ at the daughter row $`p + \xi_q`$,
- $`-1`$ at the daughter row $`p - \xi_q`$,

and the two legs are co-located: same $`x`$, symmetric in momentum about the
parent. The QLE fixes this deposition completely. Nothing in this note changes
what is deposited. Everything in this note is about *how the deposit is paid
for*.

---

## 4. The two realisations

Here is the whole idea, and it takes three sentences.

$`E`$ is a difference. So "deposit $`+1`$ here" has two readings: **add a
positon**, or **remove a negaton**. They are indistinguishable in the
observable and opposite in every other book.

Applying that to both legs at once gives the two realisations. The
specification's own verbs are the right ones, and they are the verbs of
ionisation chemistry rather than of accounting:

**Emissive — a bound pair is *ionised*.** Add a positon at $`p + \xi_q`$ and a
negaton at $`p - \xi_q`$. The two new bodies come from somewhere: a bound sea
pair at the parent row $`p`$ comes apart, and its two constituents are the two
new bodies. One bound pair is consumed; two free bodies enter circulation.

**Absorptive — two free bodies are *bound*.** Remove a negaton at
$`p + \xi_q`$ and a positon at $`p - \xi_q`$. The two removed bodies do not
vanish; they bind into a neutral pair at the parent row $`p`$. Two free bodies
leave circulation; one bound pair is created.

These are the two directions of a single exchange between the bound population
and the free one — dissociation and recombination, as in an ionisation
equilibrium. Both run forward in time. Neither is the other's time-reverse:
time reversal would flip every momentum, and nothing here does. They are
opposed in the sense that a forward and a reverse reaction are opposed, which
is why §7 can be an equilibrium argument.

![The two realisations](https://raw.githubusercontent.com/billpage/wpmw/output/figures/emission_absorption_realisations.png)

The picture is the same picture with the arrows reversed, because the exchange
is the same exchange run the other way.

That phase-space cartoon has one serious defect: it invites the reader to see
the diagonal arrows as *paths*, along which something travels from the parent
momentum row to the daughters. Nothing travels. The honest picture is in
space-time, where all three populations are visible as worldlines and a
momentum is a **slope**:

![Space-time worldlines](https://raw.githubusercontent.com/billpage/wpmw/output/figures/emission_absorption_worldlines.png)

Read the left panel. A bound pair moves classically with slope $`p`$. At the
event its worldline **ends**, and two new worldlines **begin** at the same
space-time point, with slopes $`p + \xi_q`$ and $`p - \xi_q`$. In the right
panel two worldlines end and one begins. The faint lines are the ambient sea:
real, present everywhere, at every momentum.

Momentum conservation is now a geometric statement — the single line's slope is
the mean of the pair's — and the thing that looked like a jump between momentum
rows has become a fork. §9 makes that argument properly, since it is the
objection most likely to stop a careful reader.

The demo runs both realisations on the same state (part B):

```
   process                      events       f   |dE| moved        dN        dS
   emissive   (pair -> two)    0.00560   0.000    9.577e-06  +0.01121  -0.00560
   absorptive (two -> pair)    0.00558   0.924    9.664e-06  -0.00947  +0.00473

     emissive   (pair -> two): dN/dS = -2.000000
     absorptive (two -> pair): dN/dS = -2.000000
```

Read the columns. The amount of $`E`$ moved agrees to under one per cent — the
residue is the sequential-allocation difference of the tau-leap, and §11 shows
it is $`O(\Delta t)`$ and converges away. The ledger columns do not agree at
all: $`\Delta N`$ and $`\Delta S`$ have opposite signs. And the ratio
$`\Delta N / \Delta S = -2`$ is exact in both cases.

That $`-2`$ is the stoichiometry of the exchange — one bound pair for two free
bodies — and it is the hinge of everything that follows.

**Vocabulary.** The specification calls this a *jump* — a mediated transfer, in
which no body's momentum changes discontinuously. It reserves *hop* for the
excluded variant in which a body is physically moved between momentum rows.
§10.

---

## 5. Three populations, and the one identity

Three populations are now in play: the bound pairs $`S`$, the free bodies
$`N`$, and the signed difference $`E`$. It is tempting to think of these as
three independent stocks that must each be kept in balance. They are not.

### 5.1 Theorem J1: the pair count is conserved

Every event does one of exactly two things:

- emissive: $`\Delta N = +2`$, $`\Delta S = -1`$,
- absorptive: $`\Delta N = -2`$, $`\Delta S = +1`$.

So define

```math
P \;=\; S \;+\; N/2 .
```

Then $`\Delta P = 1 - 1 = 0`$ for an emissive event and $`\Delta P = -1 + 1 = 0`$
for an absorptive one. **The event channel cannot change $`P`$ at all**, no
matter what mixture of the two it uses.

Part D verifies this in a full streaming run with the bilinear recombination
sink also active (the `clamp` column is explained in §12):

```
      rho   clamp        f  rel |dP| per step   max |du+ - du-|
      1.0   False   0.4277          3.487e-16         5.477e-16
      5.0   False   0.6280          3.485e-16         1.332e-15
     20.0   False   0.7148          3.475e-16         3.553e-15
```

Machine precision, and — note this — at $`f = 0.71`$ just as much as at
$`f = 0.43`$. The identity is not a property of the equilibrium. It holds
everywhere.

This is worth pausing on because it simplifies the picture enormously. There is
**one** conserved quantity here, the total pair count, counting a bound pair as
one and two free bodies as one. Emission and absorption neither create nor
destroy pairs; they only move them across the bound/free boundary. The
absorptive fraction $`f`$ is not a rate constant in a two-population system. It
is the **partition parameter** of a single population — what fraction of a
fixed stock is currently in circulation rather than in the reservoir.

That in turn makes the specification's ledger identity almost obvious rather
than surprising: if a single number controls the partition, then of course one
value of it holds both halves stationary at once.

### 5.2 Theorem J2: but only globally

Now the correction, and it matters.

An event debits the sea at the **parent** row $`p`$ and credits bodies at the
**daughter** rows $`p \pm \xi_q`$. Those are different cells. So pair count
*flows* between momentum rows, and the identity of §5.1 is a statement about
the sum over all cells, not about any one of them.

Part D, same substep, measured cell by cell:

```
   max per-cell |dP|                    1.162e-04
   |sum over cells of dP|               2.798e-14
   ratio                                4.2e+09
```

Nine orders of magnitude. Locally the sea and the bodies are two genuinely
independent fields with transport between them.

![The ledger identity](https://raw.githubusercontent.com/billpage/wpmw/output/figures/emission_absorption_ledger.png)

The practical consequence is a warning. It is tempting to argue from J1 that
the system has only one degree of freedom, and therefore that certain
behaviours — oscillations, travelling structure in the sea — are impossible.
**That argument is invalid.** J1 constrains the totals and says nothing about
the local dynamics, which has two fields and can do what two-field systems do.
The local sea level is a live dynamical variable.

---

## 6. Why the fraction is one half

Look again at what an event does to the two species *separately*.

An emissive event creates one positon **and** one negaton. An absorptive event
destroys one positon **and** one negaton. Neither ever moves one species
without moving the other by the same amount in the same direction. Part C:

```
   max |du+ - du-| over 40 substeps     1.332e-15
   measured absorptive fraction f       0.854768
   measured dN                          -0.299913
   predicted 2(1 - 2f) n_ev             -0.299913
   residual                             +1.471e-14
```

So over $`n`$ events with absorptive fraction $`f`$,

```math
\Delta N \;=\; 2\thinspace (1 - 2f)\thinspace n , \qquad
\Delta S \;=\; (2f - 1)\thinspace n .
```

Both vanish at $`f = 1/2`$ and at no other value. The reasoning needs nothing
about the potential, the reach, the packet, or the depth of the sea. It is
arithmetic on $`\pm 1`$.

The most useful way to say it: whichever species is locally in the **minority**
gains exactly one body per emissive event and loses exactly one per absorptive
event. Its population is stationary precisely when the two are equally
frequent.

---

## 7. Why it finds one half by itself

Nothing so far explains why $`f`$ *should* be one half. It explains what
happens if it is. The mechanism that gets it there is the part worth
understanding, and it is not a control law.

The rule the algorithm follows is local and blind:

> Settle absorptively if a partner of the right species happens to be sitting
> at both daughter rows. Otherwise settle emissively.

Nothing anywhere measures $`f`$, targets $`f`$, or adjusts anything to reach
$`f = 1/2`$. And yet:

- absorption **consumes** minority-species bodies — that is what "remove a
  partner" means;
- emission **produces** them — every emissive event makes one of each species,
  and where the majority is already abundant it is the minority that has been
  added to;
- so a shortage of minority bodies forces emission, which manufactures the very
  thing that was short, and an abundance permits absorption, which eats it.

A thermostat, built out of nothing but supply. If the mix runs too emissive,
bodies accumulate, partners become easy to find, and $`f`$ rises. If it runs too
absorptive, bodies drain, partners become scarce, and $`f`$ falls. The fixed
point is where production equals consumption, which by §6 is $`f = 1/2`$.

The claim that this is genuinely a *local* mechanism, and not something the
packet's motion through the potential is secretly arranging, can be tested.
Clock the relaxation against **cumulative events** rather than against time,
and run it in four regimes whose event rates differ by a factor of fifteen —
the packet parked on either emitting lobe, parked in a quiet region, and moving
freely (part E):

```
   packet                         t_end       1       2       5      10      20      40
   parked r = -1.5 (lobe II)       10.4  0.8861  0.8419  0.7646  0.5960  0.5241  0.5034
   parked r = +1.5 (lobe III)      10.4  0.8875  0.8384  0.7550  0.6061  0.5200  0.5019
   parked r = -6   (quiet)        187.0  0.8129  0.7756  0.6846  0.6069  0.5404  0.5162
   moving  r = -6                   9.8  0.8141  0.7731  0.7065  0.6144  0.5204  0.5061
```

All four collapse onto one curve. The wall-clock times differ by a factor of
nineteen; the event-clock trajectories are the same. Transport sets how *fast*
events happen and contributes nothing else. The regulation is event by event.

![The thermostat](https://raw.githubusercontent.com/billpage/wpmw/output/figures/emission_absorption_thermostat.png)

Every one of these runs starts at $`\rho = 20`$, badly over-padded, with
$`f`$ near $`0.88`$. None is told where to go. All arrive.

---

## 8. Why it stops just short

The measured attractor sits slightly *below* one half, and the shortfall
narrows as the coherence reach grows. The reason is in the word "both."

An absorptive event needs a partner at **both** daughter rows simultaneously,
because §10 forbids realising the two legs differently. Part F
measures how often each leg alone could be satisfied, and how often both can be
at once:

```
      rho    leg A    leg B     both    A x B        f
      3.0    0.819    0.752    0.607    0.616    0.607
     10.0    0.826    0.782    0.657    0.645    0.657
     20.0    0.823    0.785    0.667    0.646    0.667
```

Either leg on its own finds its partner about four times in five. Requiring
both at once drops that to about two in three. Absorption is rationed by
whichever leg is harder, and the resulting deficit is the **conjunction cost**.

**What this does not show.** The project has elsewhere attributed the
conjunction cost to strong *anti-correlation* between the two legs'
availability. These numbers give that no support: the joint figure sits
slightly below the product at $`\rho = 3`$ and slightly *above* it at
$`\rho = 10`$ and $`20`$. A mean of a minimum against a product of means is not
a correlation test, so this refutes nothing — but it establishes nothing
either, and the explanation should be treated as open until measured properly.
Logged as **J-SP1**.

What *is* established is simpler and enough for a tutorial: requiring two
things at once is harder than requiring one, and that is where the shortfall
comes from.

---

## 9. The objection a careful reader will raise

> *If no world-particle ever changes its momentum discontinuously, how can
> ionising a sea pair at row $`p`$ put its two constituents at rows
> $`p \pm \xi_q`$? Surely that is exactly a discontinuous momentum change — and
> for two bodies at once.*

It certainly looks like one. The bound pair sits at row $`p`$; its two constituents
end up at rows $`p \pm \xi_q`$. If those constituents were world-particles
with momenta of their own, those momenta just changed discontinuously.

The resolution is §2's careful sentence. A bound sea pair is *one* neutral
object, not two world-particles. Its constituents have no worldlines while
bound. What happens at an emissive event is that two worldlines **begin**, at
$`p + \xi_q`$ and $`p - \xi_q`$ respectively. A worldline that begins at a
momentum has not moved to that momentum. Birth is not a velocity change.

Symmetrically, at an absorptive event two worldlines **end**, and a neutral
object appears at the parent row. The bodies were not carried to $`p`$; they
ceased to be bodies.

This is why the space-time figure of §4 is the one to trust and the phase-space
cartoon is the one to be careful with. In space-time the two lines fork out of a
point. The phase-space picture, read carelessly, draws horizontal segments
connecting momentum rows — exactly the transport that does not happen.

Momentum still balances: the pair at $`p`$ carries $`2p`$, and the daughters
carry $`(p + \xi_q) + (p - \xi_q) = 2p`$. In the space-time picture that is
just the statement that the parent's slope is the mean of the daughters'.

The same reading disposes of a related worry. One might ask what *holds* a sea
pair together, given that the project admits no new forces. Nothing does, in
the sense of a binding force. A bound pair is a single neutral member of a
third population; it is not a composite held in place by an interaction. Its
constituents are not two things that have been joined, which is why they have
no separate worldlines to be joined along.

---

## 10. What is excluded, and why

Two further readings are available and both are ruled out.

**Settling the two legs differently.** Pay one leg by creating a body and the
other by removing one. The deposition in $`E`$ is still correct. But count the
bodies: one has appeared at one daughter row and one has disappeared at the
other, and the difference has to have been carried between two momentum rows
separated by $`2\xi_q`$. Nothing can carry it. A world-particle would have had
to move — a **hop**, not a jump — and individual body-momentum conservation
fails. The realisation is therefore a property of the **event**, not of a leg.

In space-time the reason is immediate: a mixed realisation would have one
worldline beginning at one daughter row and another ending at the other, with
nothing connecting them, and the momentum books would not close. A pair is
either coming apart or coming together. It cannot do half of each.

**Relocating a body instead of pairing.** The third channel, in which a body is
simply moved from one row to another, is genuinely available and is analysed
elsewhere in the project. It relaxes the "both daughters" requirement of §8
from a conjunction to a disjunction, which helps the supply problem — at the
cost of individual body-momentum conservation, which is the thing the ontology
exists to protect. It is not part of the specified algorithm.

---

## 11. Numerical verification

All figures and tables above are produced by
[`src/demo_emission_and_absorption.py`](../../src/demo_emission_and_absorption.py):

| part | what it establishes |
|---|---|
| A | $`K_{\mathrm{res}}`$ real, odd, moving neither worlds nor momentum; $`(q, -q)`$ is one event |
| — | four figures, including the space-time worldlines of §4 |
| B | the two realisations agree in $`E`$ and differ in the populations, with $`\Delta N/\Delta S = -2`$ exactly |
| C | both species move together; $`\Delta N = 2(1-2f)n`$ to $`1.5 \times 10^{-14}`$ |
| D | J1 exact to $`3.5\times10^{-16}`$ globally; J2, local violation $`4\times10^{9}`$ larger |
| E | the relaxation of $`f`$ collapses on the event clock across four rate regimes |
| F | the conjunction cost, and the withdrawal of the anti-correlation claim (J-SP1) |
| G | the observable's independence from the populations is $`O(\Delta t)`$; see J-SP2 |

Run as:

```
WPMW_OUTPUT=... PYTHONPATH=src python3 src/demo_emission_and_absorption.py
```

---

## 12. Open items

### J-SP1 — the conjunction cost is measured, its explanation is not

§8 establishes the cost. The project has elsewhere attributed it to strong
*anti-correlation* between the two legs' partner availability, and these
numbers give that no support: the joint figure sits slightly below the product
of the two single-leg figures at $`\rho = 3`$ and slightly **above** it at
$`\rho = 10`$ and $`20`$. A mean of a minimum against a product of means is not
a correlation test, so nothing is refuted — but nothing is established either.
A proper test would correlate the two caps directly at fixed demand.

### J-SP2 — a proposed change to the allocation loop

**This is a separate subject from the tutorial and is recorded here only so it
is not lost.** It is a proposal for the maintainer's judgement, not a finding
about the ontology.

§5.1 of the specification requires $`u^{\pm} \ge 0`$ pointwise and clamps after
every allocation. That clamp is the only operation in the specified world form
that can break Theorem J1. Part D of the demo runs every case both ways:

```
      rho   clamp        f  rel |dP| per step   max |du+ - du-|
      1.0   False   0.4277          3.487e-16         5.477e-16
      1.0    True   0.4287          1.501e-07         2.384e-04
     20.0   False   0.7148          3.475e-16         3.553e-15
     20.0    True   0.7158          2.437e-07         4.790e-04
```

With the clamp on, J1 degrades from machine precision to $`10^{-7}`$ and the
two species stop moving together at $`10^{-4}`$.

Part G is the sharper measurement. It prepares the same $`W`$ at $`\rho = 1`$
and $`\rho = 20`$ — twenty times as many bodies, identical $`E`$ — and asks how
far apart the two observables have drifted at $`t = 2`$:

```
         dt     clamp off   ratio      clamp on   ratio
      0.020    2.7989e-03            1.0688e-02
      0.010    1.3980e-03    2.00    1.1951e-02    0.89
      0.005    6.9850e-04    2.00    1.2567e-02    0.95
```

Clamp off, the difference halves cleanly with $`\Delta t`$ — it is tau-leap
allocation error amplified by $`N/\|E\|`$, and it converges to zero, which is
what the theory requires. Clamp on, it stalls near $`10^{-2}`$ and does not
converge at all.

**Diagnosis.** The clamp is not repairing an over-drawn allocation. The partner
caps are read live and clipped at zero, so the allocation never over-draws. It
is repairing the small negative populations that **spectral transport** leaves
behind through ringing — a transport artefact, repaired inside the event loop,
where it corrupts event arithmetic that is otherwise exact.

**Proposal, for judgement.** Clip the partner **caps** at zero inside the
allocation — local, physically obvious, since a negative supply is no supply,
and it leaves the population arithmetic exact — and let transport ringing be
repaired at the transport substep, where the specification already mandates a
separate repair restoring $`N \ge |E|`$.

**A corollary worth flagging separately.** The specification's own consistency
test reads a ledger residual above the floor as evidence that an implementation
has admitted a hop. On this measurement it can equally mean the clamp is
firing. That test may need a second named failure mode.

### Untouched here

The local sea level is a live dynamical variable in a genuinely two-field
system (§5.2). What it does — whether it supports oscillation or travelling
structure, and how the free and bound populations exchange locally rather than
in total — is not addressed by this note.

---

## 13. Where to go next

- [`compensated_liouville_algorithm.md`](../algorithm/compensated_liouville_algorithm.md)
  §5 — the normative specification this note unpacks.
- [`sea_population_equilibrium.md`](../analysis/sea_population_equilibrium.md)
  — the derivation, the failure of emissive-only unravelling, and the attractor
  result.
- [`species_sectors_and_annihilation.md`](../analysis/species_sectors_and_annihilation.md)
  — what the sea is, and why the dark bound pair is not the crystal shift.
- [`eckart_barrier_compensated.md`](../analysis/eckart_barrier_compensated.md)
  §8 — what all this looks like as tunnelling.
