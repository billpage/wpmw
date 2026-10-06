# Review of `WignerParticlesSinSpawnV1` (David Cyganski, 17 April 2020)

**A review of David Cyganski's curved-trajectory notebook, written for him
and for the WPMW project. The notebook evolves a Gaussian wavepacket in a
quadratic-plus-sine potential. The quadratic part is applied as an exact
classical force, and the sine as signed pair spawning. That split is sound,
and the notebook is prior art for it. The code, however, has twelve bugs. Two
of them change the physics at any strength: a numpy view that makes the
phase-space map contract, and sign handling that applies the sine with the
wrong sign or not at all. At the notebook's own parameters the sine is far
too weak for a comparison to tell right from wrong. With the bugs fixed and
the sine made strong enough to matter, the algorithm reproduces the exact
Schrödinger solution to sampling noise. It even beats the classical
trajectory answer, so it captures the quantum part of the sine. The notebook
as written does worse than leaving the sine out.**

---

[PDF](https://github.com/billpage/wpmw/releases/latest/download/sinspawn_v1_review.pdf)

The fixed notebook: [source, outputs stripped](../notebooks/sinspawn_v1_fixed.ipynb) ·
[rendered, `output` branch](https://nbviewer.org/github/billpage/wpmw/blob/output/notebooks/sinspawn_v1_fixed.ipynb)
(about 23 MB, too large for GitHub's own notebook viewer). Companion code:
[`src/demo_sinspawn_review.py`](../../src/demo_sinspawn_review.py), with the
simulator in [`src/wpmwlib/sinspawn.py`](../../src/wpmwlib/sinspawn.py), and
[`src/scan_sinspawn_v1_original.py`](../../src/scan_sinspawn_v1_original.py).
No theorems are stated, so the note has no prefix.

## 0. In brief

1. **What it is.** Signed particles in phase space $`(x,k)`$, with $`k=p/\hbar`$.
   Each step they rotate exactly under the quadratic potential. Each particle
   spawns a $`\pm`$ pair with probability `FracSpawn`$`\cdot\cos(qx)`$. The
   children are kicked by $`\pm q/2`$, and opposite signs cancel in a box.
   This is the force-plus-branching split of the 2026-09-29 call, five years
   earlier. Section 2 lists what the notebook adds to the account of that
   call and what it corrects.
2. **Bugs (§3).** There are twelve, labelled A–L (D has two parts), plus the
   Python 2 syntax.
   - **D (numpy views).** The "last" coordinates are views, not copies. So
     the momentum update reads the *new* position, and the map has
     determinant $`\cos^2(\omega\,dt)`$ instead of 1. Phase-space area
     shrinks by 2.2% over the run. The box-change test always comes out
     empty, so no particle is ever annihilated.
   - **A and C (signs).** Spawning ignores the sign of $`\cos(qx)`$, and the
     children are placed opposite to the notebook's own formula.
   - **J, K and L.** These would crash or stop the notebook. L already
     stopped it: the cells comparing with the analytic solution never ran.
3. **The notebook's test cannot discriminate (§5.1).** At its parameters the
   sine force is 0.07% of the harmonic force. It shifts the position density
   by $`L^1 = 0.009`$. The sampling noise of 838 particles is about 0.14. Even
   so, the saved output carries the view bug's fingerprint. Its final mean
   position is 15.71 nm, which is what the contracting map predicts (15.70 to
   15.72 nm). The exact value is 15.87 nm.
4. **A test that matters (§5.2).** Take the 16th harmonic with
   $`V_p = 24.9`$ meV, so the sine force is 20% of the harmonic one.
   - The fixed algorithm matches the exact solution to $`L^1 = 0.046`$ at
     $`N_0 = 3.2\times10^5`$ (the error falls as $`N_0^{-1/2}`$). It gets the
     mean position at 1 ps to 5 pm.
   - It beats the classical trajectory answer ($`L^1 = 0.175`$).
   - Leaving the sine out costs $`L^1 = 1.02`$.
   - Reintroduced one at a time, bugs A and C each cost about 1.2–1.4.
   - With all bugs in, the notebook scores 1.25.
5. **Recommendations (§6).** Copy before you overwrite. Use one spawn routine
   with $`\lvert\Gamma\rvert`$ and a signed child. Set $`V_p`$ physically.
   Test against Schrödinger *with* the sine, at a strength where sine-on and
   sine-off differ by much more than the noise. Check the 2026
   virtual-photon notebook for the same two bug classes.

---

## 1. The notebook

### 1.1 Provenance

The notebook continues `WignerParticlesNoSpawn` V1–V3 (30 March–10 April
2020), which were free motion, a linear potential and a quadratic potential.
Its header says V1 "was frozen when the most important functions had been
made functional, but prior to completion". It lists five to-dos, including
"No comparison to Schrodinger solution of same problem" and "No physical
scaling in the particle creation rate constant". The last cell announces
"v4 of this series that will consider third order potentials". It was run
under Python 2. The saved outputs come from that run, with 838 particles and
218 s for the main loop.

### 1.2 Parameters (SI, electron)

| Symbol | Notebook | Value |
|---|---|---|
| box | `LX`, `NX` | 400 nm, 1024 cells, $`dx = 0.391`$ nm |
| momentum cell | `dk` $`=2\pi/L_X`$ | $`1.571\times10^{7}`$ m$`^{-1}`$ |
| time step, steps | `dt`, `TimeSteps` | $`10^{-16}`$ s, 9999 steps (to 1.0 ps) |
| quadratic | `V2` $`=250\times10^6/L_X`$ V/m$`^2`$ | $`\omega = 1.483\times10^{13}`$ s$`^{-1}`$, period 0.424 ps |
| packet | `x0`, `a` | $`x_0 = -25`$ nm, $`\sigma_x = 2`$ nm, $`\sigma_k = 2.5\times10^{8}`$ m$`^{-1}`$, at rest |
| sine | `nharm`, `FracSpawn` | $`q = 2\pi\cdot 2/L_X`$, spawn probability $`10^{-4}\cos(qx)`$ per step |
| particles | `IPNUM`, `IPMULT` | target $`10^3`$ (838 created), arrays sized $`3\times10^3`$ |

### 1.3 The algorithm, and the equation it should solve

The potential energy is $`U(x) = q_eV_2x^2 + V_p\sin(qx)`$. The Wigner
equation for it is exactly

```math
\partial_tW = -\frac{\hbar k}{m}\partial_xW + \frac{U_2'(x)}{\hbar}\partial_kW
 + \Gamma(x)\Bigl[W\bigl(x,k+\tfrac q2\bigr) - W\bigl(x,k-\tfrac q2\bigr)\Bigr],
\qquad \Gamma(x) = \frac{V_p}{\hbar}\cos(qx).
```

This is Lemma T1 of the
[Takabayasi supplement](takabayasi_1954_stochastic_picture.md),
with $`U_2 = q_eV_2x^2`$. The first two terms are the classical Liouville
flow in the quadratic potential. They are exact because the Moyal series
stops at second degree, which is the notebook's own argument in its
Liouville-equation markdown (cell 58). The last term is a signed pair kernel.
For a positive parent with $`\Gamma>0`$ it puts a positive child at
$`k-q/2`$ and a negative one at $`k+q/2`$. Each event moves the first moment
by $`-q`$, so the mean force is $`-V_pq\cos(qx) = -U_{\sin}'`$, as it must be.
The notebook's cell 59 states the same kernel from David's note *Wigner
Particle Model with Curved Trajectories*, as
$`\propto\cos(2\pi nx/L)\,[W(x,k+n\pi/L) - W(x,k-n\pi/L)]`$. The prefactor
written there, $`V_p\pi`$, should be $`V_p/\hbar`$. Using `FracSpawn` as
$`\Gamma_{\max}dt`$ therefore means $`V_p = `$ `FracSpawn`$`\cdot\hbar/dt`$,
which is 0.658 meV.

The main loop (cell 70) does three things each step:

1. **Spawn.** Each particle draws a uniform number. If it falls below
   `FracSpawn`$`\cdot\cos(qx)`$, a pair is born at the parent's position with
   kicks of $`\pm`$`kbump`$`=\pm n/2`$ cells, which is $`\pm q/2`$. If the
   box of the new opposite-signed child already holds a particle of the
   other sign, the two cancel at birth. The surviving existing particle is
   then *kinked*, that is, moved to the other child's box.
2. **Propagate.** Each particle takes the exact harmonic rotation over
   $`dt`$.
3. **Re-box.** Each particle that changed box is moved in the per-box index
   lists. If its new box holds a particle of the opposite sign, the two
   annihilate.

---

## 2. What is new relative to the call of 29 September

What David said about this code on the call of 29 September 2026 was
recorded from memory, before the code arrived (summarised in §0 and §4 of
[`poisson_kicks_and_pair_branching.md`](poisson_kicks_and_pair_branching.md)). The notebook
confirms and sharpens it.

**Prior art, confirmed.**

- The split is explicit: the quadratic part is an *exact* force (curved
  trajectories, closed-form rotation), and the sine is spawning.
- "Kinking" is cancellation at birth, represented as a momentum jump of the
  survivor. Annihilation on box change is the other half. Together they are
  what WPMW calls recombination.
- Cell 100 calls the particles "entire worlds". It describes the collision
  term as particles scattering off "the surface of a background sea" and
  raising pairs from it. That is the WPMW reading, in 2020.
- `dk` $`=2\pi/L_X`$ is hard-wired from the FFT version of his code. So the
  lowest usable harmonic is $`n = 2`$, and the kick is one cell,
  $`\hbar\,dk = 2\pi\hbar/L_X`$. That is twice the crystal step $`\pi\hbar/L`$
  that a sine of period $`L`$ would have on the WPMW lattice.

**Corrections to that account.**

- *Not a periodic box.* Particles move in an open $`400`$ nm patch, and cell
  76 counts those that leave it. There is no Fourier parity argument in the
  notebook. The justification given is the truncated Moyal series (cells
  58–59). The parity argument may be in the later notes David described,
  which we have not seen.
- *Not abandoned for cost, as far as this notebook shows.* It was frozen
  "prior to completion" with a to-do list, and the next step it announces is
  a cubic potential.
- *Dates.* The headers settle them. The NoSpawn series ran 30 March–10 April
  2020, and SinSpawnV1 is dated 17 April 2020.
- *Two meanings of "garbage collection".* In the notebook it means reusing
  array slots after annihilation, which is not yet done. In the call, Bill
  used it for annihilation itself ("absorption as garbage collection").
- *X-SP3 is untouched.* The overshoot in open item X-SP3 of the
  [Poisson-kick supplement](poisson_kicks_and_pair_branching.md)
  belongs to the 2026 virtual-photon notebook. This 2020 notebook uses a
  fixed-$`dt`$ Bernoulli clock evaluated at the current position. That clock
  has no frozen-rate bias, because $`\Gamma_{\max}dt\ll1`$.

---

## 3. Bugs

Cell numbers are 0-based indices into the original notebook (cell 70 is the
main loop). "His run" means the saved Python 2 execution.

| ID | Where | What | Effect in his run |
|---|---|---|---|
| **A** | 70, spawn test | `rand < FracSpawn*cos(qx)` never fires where $`\cos<0`$, and the children's signs never flip with $`\mathrm{sgn}\cos`$ | none: $`\cos(qx)>0.64`$ wherever the packet goes |
| **B** | 70, negative parents | reads `WPxkp[pindex]` (positive array) instead of `WPxkn` | negative parents spawn at the positions of positive ones |
| **C** | 70, child placement | positive parent with $`\cos>0`$: positive child at $`k+q/2`$, negative at $`k-q/2`$. Cell 59 and the QLE say the reverse | the sine acts with the wrong sign, as if $`V_p\to-V_p`$ |
| **D1** | 70, propagation | `xlastp = WPxkp[0:ppnum, xIndex]` is a *view*. After `x` is overwritten, the `k` update uses the new `x` | the map has determinant $`\cos^2\omega dt`$: area $`\times0.978`$, orbit radius $`\times0.989`$ over the run |
| **D2** | 70, box change | `xboxlastp`, `kboxlastp` are views too, so the change test compares new with new | `pchangelist`, `nchangelist` always empty: no annihilation ever, and the box lists freeze at birth |
| **E** | 70 | the negative box-change loop is indented inside the positive one | masked by D2. Fixing D2 alone would raise `ValueError` |
| **F** | 70 | dead particles (`alive = 0`) still spawn | masked by D2 (no particle dies) |
| **G** | 49, initial sampling | negative initial particles are appended to `WPxkpLists` | latent: a Gaussian has $`W\ge0`$. Wrong for any state with negative regions |
| **H** | 70, track bookkeeping | track slot `index % skipnum` should be `index // skipnum` (and only if the particle is tracked) | masked by D2 (no annihilations to record) |
| **I** | 52, `MakeWignerHist` | `floor` where the boxes use `round` | every histogram shifted by $`-dx/2`$, $`-dk/2`$: the initial mean reads $`-25.194`$ nm, not $`-25.000`$ |
| **J** | 70 (and 49) | the row is written *before* the capacity guard. At `ppnum` $`=`$ `IPMULT*IPNUM` the next write is out of bounds | latent: `IndexError` once 3000 rows are used (his run used 2436) |
| **K** | 67 | track arrays have `IPMULT*ppnum/skipnum+1` $`=315`$ columns, but up to `IPMULT*IPNUM/skipnum` $`=375`$ tracked slots can exist | latent: crash past 2520 particles, 84 short of his final 2436 |
| **L** | 90, 94, 98 | `pnum` is never defined | `NameError`: the comparisons with the analytic solution have no output. Cell 99's "essentially a perfect match" is not backed by a run of this notebook |

Python 2 constructs also stop it under Python 3:

- `print x` statements;
- `NK = NX/2`, which becomes a float;
- float indices into object arrays;
- `hist(normed=...)`;
- `%%script false`, which now raises unless given `--no-raise-error`.

### 3.1 D1: why a view makes the map contract

The exact rotation over $`dt`$, with $`s = \sqrt{2q_eV_2m}`$ (the code's
`sqrtvmq`), $`c=\cos\omega dt`$ and $`\sigma = \sin\omega dt`$, is

```math
x' = c\,x + \frac{\hbar\sigma}{s}k,\qquad k' = c\,k - \frac{s\sigma}{\hbar}x .
```

With the view, the second line reads $`x'`$ in place of $`x`$:

```math
k' = (c-\sigma^2)\,k - \frac{s\sigma c}{\hbar}x,
\qquad
\det\begin{pmatrix} c & \hbar\sigma/s\\ -s\sigma c/\hbar & c-\sigma^2\end{pmatrix} = c^2 .
```

Over 9999 steps, $`c^{2\cdot9999} = 0.9783`$, and orbits spiral inward by
$`c^{9999} = 0.9891`$. A classical point started at $`x_0`$ ends at
$`x = 15.702`$ nm under the buggy map and at $`x_0\cos\omega T = 15.890`$ nm
under the exact one.

**The fingerprint in the saved output.** Cell 86's saved plot data has final
mean $`\langle x\rangle = 15.518`$ nm. Bug I shifts it by $`-dx/2`$, and
undoing that gives **15.714 nm**. The initial mean, corrected the same way,
is $`-24.999`$ nm, so the deterministic cell-rounded sampling reproduces
$`x_0`$ to about 1 pm. With the sine included, the exact answer is 15.866 nm,
and the fixed copy (§4) lands at 15.869 nm. A re-run of the original loop
with other random numbers (§3.2) ends at 15.735 nm, and our
re-implementation with all bugs switched on at 15.721 nm (§5.1). So the
saved run sits on the contracting map's prediction, 0.15 nm short of the
truth.

### 3.2 D2: what happens when nothing is ever re-boxed

The per-box index lists keep each particle in the box where it was born. The
cancellation-at-birth search then finds particles that left that box long
ago. The original loop, ported to Python 3 and otherwise unchanged, with
counters added (`scan_sinspawn_v1_original.py`), counted the following over
the same 9999 steps:

- 178 kinks, of which 174 moved a particle that was no longer in the
  searched box (median distance 10 boxes);
- zero box changes for either sign.

So the kinks are momentum jumps of $`2\hbar\,dk`$ applied to essentially
random particles. Annihilation never happens, and the population only grows:
838 at the start, $`2436+1598`$ rows at the end of the saved run
($`2364+1526`$ in the re-run).

### 3.3 A and C: the sign rule, in one line

For a parent of sign $`s`$ at $`x`$, spawn with probability
$`\lvert\Gamma(x)\rvert dt`$. The child of sign $`s\cdot\mathrm{sgn}\,\Gamma(x)`$
goes to $`k-q/2`$, and the child of the opposite sign to $`k+q/2`$. The
notebook writes four hand-specialised copies of this (two parent signs, each
with two cancellation cases). In them, $`\mathrm{sgn}\,\Gamma`$ is dropped (A),
the placement is reversed (C), and the negative copy reads the wrong array
(B). One function with the sign as an argument removes all three.

---

## 4. The fixed copy, at the notebook's own parameters

[`sinspawn_v1_fixed.ipynb`](../notebooks/sinspawn_v1_fixed.ipynb) is David's notebook with every
change marked `FIX`:

- The main loop keeps his structure and names. It spawns through one
  `spawn_pair(parent_sign, …, sgn cos)` routine (fixes A, B, C).
- It copies the last coordinates (D).
- It un-nests the negative loop (E).
- It skips dead parents (F).
- It checks capacity before writing, counts the drops and reports them (J).
- It kills particles that leave the patch.
- It computes track slots with `//` (H).
- It reports how many kinks and annihilations occurred.

The other fixes are G, I, K and L, and the Python 3 port. Two cells are
added: a fixed random seed, and a last cell.
That cell solves the Schrödinger equation (split operator, sine included,
$`V_p = `$ `FracSpawn`$`\cdot\hbar/dt`$) and compares $`\rho(x)`$ on the
notebook's own boxes. The
[rendered copy](https://nbviewer.org/github/billpage/wpmw/blob/output/notebooks/sinspawn_v1_fixed.ipynb)
gives:

| | Original | Fixed (executed) |
|---|---|---|
| kinks at birth / annihilations on box change | 178, 174 of them stale / 0 (re-run, §3.2) | 508 / 235 |
| rows used, positive / negative | 2436 / 1598 (saved run) | 1073 / 235 |
| net signed count | 838 | 838 |
| $`\langle x\rangle`$ at 1 ps (exact 15.866 nm) | 15.714 nm (saved, half-cell corrected); 15.735 nm (re-run) | 15.869 nm |
| $`L^1(\rho_{\rm particles},\rho_{\rm exact})`$ | not computed | 0.241 |

Every negative particle created in the fixed run is annihilated by the end,
since 235 were born and 235 box-change cancellations occurred.
Figure 1 shows the change. In the original, negatons (red) accumulate and
persist on full orbits. In the fixed copy they are short-lived.
An $`L^1`$ of 0.241 is sampling noise. An ideal sampler of the exact density
with 838 independent particles would give about 0.14, and our
re-implementation at $`N_0 = 838`$ gives 0.151. At these parameters, 838
particles cannot see the sine at all. Section 5.1 shows why.

![Tracks, original and fixed](https://raw.githubusercontent.com/billpage/wpmw/output/figures/sinspawn_review_tracks.png)

*Figure 1. Tracked particles at the notebook's own parameters. Left: the
saved output of the original. Right: the fixed copy, executed.*

---

## 5. Against the exact solution

The notebook loop is pure Python and takes minutes for $`10^3`$ particles.
So the comparison runs use `wpmwlib/sinspawn.py`, a numpy
re-implementation of the same algorithm, driven by Part D of
`demo_sinspawn_review.py`:

- same boxes ($`dx`$, $`dk`$, `round`);
- same fixed-$`dt`$ Bernoulli clock;
- same exact rotation;
- spawn with $`\lvert\Gamma\rvert dt`$ and a signed child;
- once per step, cancel opposite signs pairwise in each box (birth and
  box-change cancellations together).

Every bug A–D is a switch. Two things differ from the notebook. The initial
state is sampled i.i.d. from the Gaussian Wigner function, not by cell
rounding. And the "notebook" run models bug J as a counted drop at a cap of
$`3N_0`$, because the literal code would crash there.

The two codes agree at the notebook's parameters. The fixed notebook makes
$`508+235 = 743`$ spawn events from 838 particles (0.89 per particle). The
vectorized code makes 0.80 per particle at $`N_0 = 838`$ and 0.86 at
$`N_0 = 80\,000`$, and lands within 0.03 nm of the exact mean (§5.1), as the
notebook does.

The reference is split-operator Schrödinger with the sine included: 4096
points on 400 nm, the same $`dt`$, norm conserved to $`2\times10^{-12}`$.
Doubling the grid and halving $`dt`$ changes $`\rho`$ at 1 ps by
$`L^1 = 7.6\times10^{-6}`$. Errors are $`L^1`$ distances of $`\rho(x)`$ on
$`dx`$ bins. Two other baselines appear in the tables:

- **sine off**: the exact solution without the sine;
- **classical**: Newtonian trajectories in the *full* potential from the
  same initial $`W`$, with $`4\times10^5`$ samples. This is the truncated
  Moyal answer, so it misses only the quantum part of the sine.

### 5.1 At the notebook's parameters, nothing can be seen

The sine force amplitude $`V_pq`$ is $`6.6\times10^{-4}`$ of the harmonic
force at $`x_0`$. Its wavelength, 200 nm, is four times the width of the
orbit, so over the packet's range it is nearly a uniform force. And
$`(q\sigma_x)^2 = 0.004`$, so it is also nearly classical.

| Run | $`N_0`$ | $`L^1`$ at 0.25 / 0.50 / 1.00 ps | $`\langle x\rangle`$ at 1 ps (nm) |
|---|---|---|---|
| exact | — | — | 15.866 |
| sine off (exact) | — | 0.011 / 0.003 / 0.009 | 15.890 |
| classical, full force | $`4\times10^5`$ | 0.010 / 0.006 / 0.006 | 15.866 |
| fixed | 838 | 0.099 / 0.159 / 0.151 | 15.892 |
| fixed | 80 000 | 0.011 / 0.015 / 0.017 | 15.856 |
| notebook, all bugs | 80 000 | 0.015 / 0.038 / 0.062 | 15.721 |

The whole effect of the sine ($`L^1\le0.011`$) is below the sampling noise
of 80 000 particles (about 0.014), let alone 838. The comparison the
notebook planned was against the pure QHO solution of Andrew (2015), and it
could not have failed. The only bug visible at these parameters is D1,
through the mean position.

### 5.2 A strength that matters

Two things are needed: a sine force comparable to the harmonic one, and a
kick comparable to the momentum width, so that the sine does more than a
classical force would. The cost of a strong sine depends on the harmonic.
Proposition X3 of the Poisson-kick supplement gives it: a force $`F`$ made
of kicks $`\pm\hbar q/2`$ needs a spawn rate of $`\lvert F\rvert/(\hbar q)`$.
At $`n=2`$ a 20% force needs $`\Gamma_{\max}dt = 0.030`$, which is eight
times the rate at $`n=16`$, and it would still be nearly classical. So we
take

```math
n = 16,\quad q = 2.51\times10^{8}\ \mathrm{m^{-1}},\quad
V_p = 0.2\,\frac{m\omega^2\lvert x_0\rvert}{q} = 3.98\times10^{-21}\ \mathrm{J} = 24.9\ \mathrm{meV}.
```

That gives the following:

- the kick is $`q/2 = 8\,dk = 0.50\,\sigma_k`$;
- $`(q\sigma_x)^2 = 0.25`$;
- $`\Gamma_{\max}dt = 0.0038`$;
- the wavelength of 25 nm is crossed about twice per half period.

| Run | $`N_0`$ | $`L^1`$ at 0.25 / 0.50 / 1.00 ps | $`\langle x\rangle`$ at 1 ps (nm) | peak pop. $`/N_0`$ | births $`/N_0`$ |
|---|---|---|---|---|---|
| exact | — | — | 12.829 | — | — |
| sine off (exact) | — | 0.167 / 0.905 / 1.015 | 15.890 | — | — |
| classical, full force | $`4\times10^5`$ | 0.010 / 0.069 / 0.175 | 12.831 | — | — |
| **fixed** | 20 000 | 0.084 / 0.143 / 0.344 | 12.943 | 3.31 | 92.6 |
| **fixed** | 80 000 | 0.035 / 0.056 / 0.098 | 12.844 | 1.67 | 63.6 |
| **fixed** | 320 000 | 0.026 / 0.035 / **0.046** | 12.834 | 1.22 | 54.9 |
| bug A only (rectified) | 80 000 | 0.781 / 0.466 / 1.379 | 8.355 | 1.27 | 32.3 |
| bug B only (wrong array) | 80 000 | 0.043 / 0.109 / 0.172 | 13.520 | 1.25 | 56.5 |
| bug C only (children flipped) | 80 000 | 0.597 / 1.212 / 1.203 | 16.224 | 1.64 | 63.7 |
| bug D1 only (view map) | 80 000 | 0.034 / 0.050 / 0.113 | 12.738 | 1.62 | 62.6 |
| notebook, all bugs, cap $`3N_0`$ | 80 000 | 0.101 / 1.166 / 1.248 | 16.615 | 3.00 | 2.0 |

For the fixed runs, the ideal-sampler noise at 1 ps is 0.031, 0.015 and
0.008 at $`N_0`$ = 20 000, 80 000 and 320 000.

![Position marginals against the exact solution](https://raw.githubusercontent.com/billpage/wpmw/output/figures/sinspawn_review_marginals.png)

*Figure 2. Position density at 0.5 ps and 1.0 ps. Top: exact (black), sine
switched off (grey dashed), classical trajectories in the full potential
(magenta dotted), the fixed algorithm, and the notebook with all bugs.
Bottom: one bug at a time.*

![Error scaling and population](https://raw.githubusercontent.com/billpage/wpmw/output/figures/sinspawn_review_scaling.png)

*Figure 3. (a) $`L^1`$ error of the fixed algorithm against $`N_0`$. Grey
lines: sine on against sine off. Magenta: classical against exact at 1 ps.
Green dashed: noise floor of an ideal positive sampler. (b) Total particle
number against time.*

**What the numbers say.**

- **The fixed algorithm is right.** From $`N_0 = 8\times10^4`$ to
  $`3.2\times10^5`$ the error at 1 ps falls by 2.1, against $`\sqrt4 = 2`$
  for pure noise, so no bias floor is visible at $`L^1\approx0.05`$. That
  fits Proposition X5 of the supplement. Bias comes from the momentum width
  of the annihilation cell, and here $`dk = 0.063\,\sigma_k`$ is fine.
- **The error is several times the ideal floor.** It is 0.046 against
  0.008, because signed weights add variance. That excess falls with density.
  Peak population drops from $`3.3N_0`$ to $`1.2N_0`$, and 97–99.7% of
  all births are cancelled.
- **It captures the quantum part of the sine.** Classical trajectories in
  the full potential get the mean right (12.831 nm) but the shape wrong by
  0.175 at 1 ps (Figure 2, top right). The fixed algorithm reaches 0.098 at
  $`N_0 = 8\times10^4`$ and 0.046 at $`3.2\times10^5`$. A kick of
  $`q/2 = 0.5\,\sigma_k`$ is a genuinely non-classical event. The branching
  carries all orders of the Moyal series for the sine (Theorem X2), and the
  run shows it.
- **A and C are gross.**
  - Rectifying the cosine leaves a net force of one sign averaged over each
    wavelength, a spurious uniform force. The packet ends 4.5 nm short.
  - Flipping the children is $`V_p\to-V_p`$. The packet ends on the other
    side of the sine-off solution.
- **B and D1 are modest at this strength,** but systematic.
  - B moves the mean by 0.7 nm.
  - D1 moves it by 0.09 nm here; its contraction grows with run length,
    whatever the potential.
- **The notebook as written is worse than no sine.** Without annihilation
  (D2) the population reaches the cap within about 0.02 ps. After that,
  46.5 spawns per initial particle are dropped. In the literal code it would
  stop with an `IndexError` there (J).

---

## 6. Recommendations, if the code is reused

1. **Never keep "old" coordinates as slices.** Use `.copy()`, or update both
   coordinates from one expression,
   `x, k = c*x + a*k, c*k - b*x`, since the right-hand side is evaluated
   before either name is rebound. The same applies to the stored box indices.
2. **One spawn routine.** Probability $`\lvert\Gamma\rvert dt`$. Child of sign
   $`s\,\mathrm{sgn}\,\Gamma`$ at $`k-q/2`$, opposite sign at $`k+q/2`$. Pass
   the parent's sign as an argument rather than copying the loop.
3. **Set the amplitude physically.** `FracSpawn` $`=V_p\,dt/\hbar`$. Keep
   $`\Gamma_{\max}dt\ll1`$, or use an event clock with thinning (§7 of the
   supplement).
4. **Test against Schrödinger with the sine,** on the particle boxes. Pick
   parameters where sine-on and sine-off differ by many times the sampling
   noise. To show that the quantum part is tested, also beat the classical
   full-force answer. The reference cell added to the fixed notebook is about
   25 lines.
5. **Annihilation.** Cancel once per step per box, pairing at random. Never
   let dead particles spawn. Compact the arrays (true garbage collection)
   instead of capping them. Size the momentum cell per Proposition X5: the
   bias falls with $`dk/\sigma_k`$.
6. **Vectorise.** The numpy version runs $`3.2\times10^5`$ particles for
   $`10^4`$ steps in about 20 minutes on one core. The loop over particles in the notebook is
   the reason it was run with 838.
7. **Check the 2026 virtual-photon notebook for the same two classes.** Look
   for slices reused after assignment, and for sign handling copied across
   branches. They produce plausible figures, as this notebook shows.
8. **NoSpawn V3.** If its loop has the same two `last` lines, the
   QHO-only results carry the same 1.1% inward spiral. We have not seen that
   notebook.

---

## 7. Files

| File | What |
|---|---|
| [`docs/notebooks/sinspawn_v1_fixed.ipynb`](../notebooks/sinspawn_v1_fixed.ipynb) | the notebook with fixes marked `FIX`, a fixed seed and a Schrödinger reference cell; outputs stripped |
| `notebooks/sinspawn_v1_fixed.ipynb` on the `output` branch ([rendered](https://nbviewer.org/github/billpage/wpmw/blob/output/notebooks/sinspawn_v1_fixed.ipynb)) | the same, executed with all outputs (Python 3.11, about a minute) |
| [`src/wpmwlib/sinspawn.py`](../../src/wpmwlib/sinspawn.py) | numpy version of the algorithm: `Setup`, `schroedinger`, `SinSpawn` with bug switches, `marginal` |
| [`src/demo_sinspawn_review.py`](../../src/demo_sinspawn_review.py) | Parts A–F: SymPy checks of the kernel and of D1, D1 in numbers, the reference's convergence, the runs of §5, the classical baseline, Figures 1–3 |
| [`src/scan_sinspawn_v1_original.py`](../../src/scan_sinspawn_v1_original.py) | the original main loop, ported to Python 3 unchanged, with counters (§3.2) |

Figure 1 is assembled by Part F from two track plots placed in the output
directory: the original notebook's saved one (the original notebook is not in
the repository) and the one in the rendered fixed copy.

Reproduce with
`WPMW_OUTPUT=<dir> PYTHONPATH=src python3 -u src/demo_sinspawn_review.py`
(about half an hour on two cores; `WPMW_QUICK=1` for a smoke test) and
`PYTHONPATH=src python3 -u src/scan_sinspawn_v1_original.py` (about two minutes).

---

## 8. What was verified, and how

| Claim | Check |
|---|---|
| Kernel, signs and $`\Gamma = (V_p/\hbar)\cos qx`$ | Lemma T1; demo Part A (SymPy): the kernel matches the Moyal series term by term through $`\partial_k^{11}`$, and its mean force is $`-U'`$ |
| D1: $`\det = \cos^2\omega dt`$ | demo Part A (SymPy); Part B iterates the buggy map from $`x_0`$ to 15.702 nm |
| D1 fingerprint | cell 86's saved data: 15.518 nm, plus $`dx/2`$, is 15.714 nm. The re-run of the original loop gives 15.735 nm, the all-bugs vectorized run 15.721 nm |
| D2: no box changes, stale kinks | `scan_sinspawn_v1_original.py`: 0 box changes, 174 of 178 kinks stale |
| J, K crash points | array sizes from cells 48 and 67 |
| L | `pnum` has no assignment in the notebook |
| Fixed notebook | rendered copy: counters, $`\langle x\rangle`$, $`L^1`$ in §4 |
| Vectorized ≡ fixed notebook | 0.80–0.86 against 0.89 spawn events per particle; both means within 0.03 nm of exact |
| Reference converged | grid ×2 and $`dt/2`$: $`L^1 = 7.6\times10^{-6}`$ at 1 ps |
| Tables of §5 | demo Parts C–E |

---

## References

- Andrew, M. (2015). The evolution of oscillator wave functions.
  [arXiv:1509.05968](https://arxiv.org/abs/1509.05968). This is the analytic
  QHO solution the notebook planned to compare with.
- Cyganski, D. *Wigner Particle Model with Curved Trajectories* (note cited
  in cells 58–59; not seen).
- WPMW, [Poisson kicks and pair branching](poisson_kicks_and_pair_branching.md):
  Propositions X3 and X5, Theorem X2, open item X-SP3.
- WPMW, [Takabayasi 1954 supplement](takabayasi_1954_stochastic_picture.md):
  Lemma T1.
- WPMW, [Species sectors and annihilation](../analysis/species_sectors_and_annihilation.md):
  Theorems D5 and D6, on when cell annihilation is exact.
