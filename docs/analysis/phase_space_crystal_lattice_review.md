# Phase-Space-Crystal-Lattice Interpretation of the Wigner Equation

**Summary review of two project documents**

- `Extended_Fokker_Planck_Eq_and_the_QLE_V2.pdf` (the analytical document)
- `Wigner_Collisions_Diagram_sozi.pdf` (the visual / interpretive document)

---

## 1. Scope of each document

### Extended Fokker–Planck Eq. and the QLE (analytical)

- Places the Wigner Quantum Liouville Equation (QLE) in the hierarchy of kinetic equations:
  General Liouville → BBGKY → Boltzmann → Extended Fokker–Planck → Fokker–Planck → Vlasov chain → (Collisionless) Liouville.
- Re-derives Wiedemann's beam-physics Fokker–Planck derivation but extends the Taylor expansion of the
  post-jump density from second order to **third order**, yielding the *extended Fokker–Planck (xFP)*
  equation.
- Matches xFP, term by term, to the Moyal series expansion of the QLE:
  - $f(x,p,t) = p/m$
  - $g(x,p,t) = \partial V/\partial x$
  - $\eta_\zeta \int \zeta^3 P_\zeta(\zeta)\thinspace d\zeta = \dfrac{\hbar^2}{4} \dfrac{\partial^3 V}{\partial x^3}$
- Shows the match requires the $\zeta$ random variable to be a **quasi-density** (signed, odd about
  zero), which motivates the signed-particle / Dirac-sea construction.
- Demonstrates the three equivalent readings of the cosine-potential update (pair spawning,
  borrow-and-give, mediated jump) and shows they reduce to a simple `np.roll`-based Python implementation.
- Cross-references supporting Maple worksheets (`xFPjumpdensityV2.mw`, `xFPjumpfromWignerPotential.mw`)
  and the Split-Fourier algorithm notebook.

### Wigner Collisions Diagram (visual)

- A 28-slide Sozi presentation walking through five interpretations of a single Wigner momentum-jump
  event:
  1. Spontaneous positon/negaton pair spawning,
  2. Pair plus immediate interaction with a forward-time particle,
  3. Pair pre-existing, one of which interacts and is diverted,
  4. Retrocausal edit by a negaton from the future,
  5. The **phase-space-crystal-lattice** interpretation.
- Illustrates how the Fourier component $V_p\cos(2\pi n x/L + \phi)$ of the potential drives spawning
  with the photon-like momentum split $\delta p = \pm n\pi\hbar/L$ (total span $h/L$ — the Abraham
  photon momentum at that wavelength).
- Connects to Reid & Drummond's *Objective QFT* forward/backward stochastic differential equation
  framework, with parametric amplification as a canonical measurement basis.

---

## 2. The phase-space-crystal-lattice interpretation

The interpretation rests on three claims (slides 14–15 of the Sozi presentation):

1. **Only positons jump.** Negatons form a static background "crystal" lattice of paired cancellers
   in phase space. The dynamical equation for laboratory observables involves only the *excess*
   positon population above this background.

2. **Justification is the gauge invariance of the QLE under $W \to W + W_0$.** Every term in the
   Moyal series contains a derivative of $W$, so adding a constant background $W_0$ leaves the
   equation unchanged:

   $$\frac{\partial W}{\partial t} + \frac{p}{m}\frac{\partial W}{\partial x} - \frac{\partial V}{\partial x}\frac{\partial W}{\partial p} + \frac{\hbar^2}{24}\frac{\partial^3 V}{\partial x^3}\frac{\partial^3 W}{\partial p^3} - \cdots = 0$$

3. **The shift can be chosen large enough to eliminate all negative regions.** Because any physical
   Wigner distribution satisfies $-2/h \le W(x,p) \le 2/h$, a shift of $W_0 = 2/h$ produces an
   admissible non-negative distribution. The negaton lattice exists only as bookkeeping for what
   was subtracted; it is "resurrected" by undoing the shift.

The operational rule is then a **single-sign stochastic process**: only positon momentum jumps,
mediated by the local potential.

---

## 3. The stochastic equations of evolution

### Single Fourier component of the potential

For $V(x) = V_p \cos\negthinspace\left(\dfrac{2\pi n x}{L} + \phi\right)$, the document derives the update rule
(section "Jump density-rate for sinusoidal potential"):

$$
\Delta W  \; \propto \; 
  \frac{V_p}{\hbar}\cos\negthinspace\left(\frac{2\pi n x}{L} + \phi + \frac{\pi}{2}\right)
  \left[\thinspace W\negthinspace\left(x,\thinspace k + \tfrac{n\pi}{L}\right) - W\negthinspace\left(x,\thinspace k - \tfrac{n\pi}{L}\right)\right]
$$

### Pure cosinusoidal potential

For $\phi = 0$, fundamental component, amplitude $V_{\max}$, the update reduces to:

$$
W(x,p,t+dt)  \; = \;  W(x,p)
   \; - \;  \frac{dt\thinspace V_{\max}}{\hbar}\thinspace\sin\negthinspace\left(\frac{2\pi x}{L}\right)
  \left[\thinspace W\negthinspace\left(x,\thinspace p+\tfrac{\pi\hbar}{L}\right) - W\negthinspace\left(x,\thinspace p-\tfrac{\pi\hbar}{L}\right)\right]
$$

In discrete form, with one momentum cell equal to $\pi\hbar/L$:

```python
W += (np.roll(W, +1, axis=0) - np.roll(W, -1, axis=0)) \
     * np.sin(np.pi * X / X_amplitude) * dt * sheight
```

> **Sign-convention note.** The V2 memo (page 18) prints the simplified form
> with the opposite sign (i.e. with `+ dt V_max sin/hbar` outside the
> bracket and `np.roll(W,-1) - np.roll(W,+1)` in Python). The general jump
> formula one line above is correct — the sign is dropped only at the
> $\cos(\theta + \pi/2) = -\sin\theta$ algebraic step. The corrected form
> here is what the algorithm in `docs/algorithm/` and the implementation in
> `src/wpmwlib/phase_space_crystal_lattice.py` use; see
> `docs/supplement/phase_space_crystal_lattice_supplement.md` §6.3 for the
> derivation and `src/sign_convention_check.py` for the regression test that
> exercises the difference on a coherent state.

### Full xFP / split-Fourier evolution

Combined with the free-streaming term:

$$
\frac{\partial W}{\partial t}
   \; = \;  -\frac{p}{m}\frac{\partial W}{\partial x}
   \; + \;  \tilde V_W(x,p,t) \thinspace\star_p\thinspace W(x,p,t)
$$

where $\tilde V_W = V_W - g(x)\thinspace\delta'(p)$ is the Wigner potential with the local force term
subtracted, and $\star_p$ denotes convolution along $p$.

### How this becomes the crystal-lattice rule

The convolution kernel $K$ for a sinusoidal component is odd-symmetric. In the simplified 7-element
example,
$$K = [0,\thinspace-1/16,\thinspace0,\thinspace0,\thinspace0,\thinspace+1/16,\thinspace0]$$
The action of $K \star W$ on a state with all population at $k = k_0$ is, after reinterpretation on
the shifted distribution, equivalent to a single positive jump:

> "1/16 of the particles at $k = 5$ **jumped** to $k = 9$, mediated by the particle at $k = 7$."

That is the single rule of the crystal-lattice model: **a positon at $(x,p)$ mediates an upward
(or downward) jump of another positon at $(x,\thinspace p-\xi)$ to $(x,\thinspace p+\xi)$**, with rate proportional
to the Wigner potential coupling and the relevant Fourier amplitude. No new particles are created;
the negaton lattice is never touched dynamically.

### General polynomial potentials

The same idea applies, but the jump density $\rho(\xi)$ must satisfy the moment problem
$$
\eta_\zeta \int \zeta^k\thinspace\rho(\zeta)\thinspace d\zeta
   \; = \;  \frac{\hbar^{k-1}}{2^{k-1}\thinspace k!}\thinspace\frac{\partial^k V}{\partial x^k}
  \qquad (k\ \text{odd})
$$
with even moments zero. For unbounded potentials (e.g. $V = x^3$) the exact $\rho$ is impulsive
(sums of derivatives of Dirac deltas). For bounded potentials (sinusoids, smooth confinements) it
is well-behaved.

---

## 4. Observations and one mild push-back

- **"Phase-space-crystal-lattice" appears to be a coinage of this work.** It is not, as far as I
  can find, standard usage in the Wigner Monte Carlo literature (Nedjalkov, Querlioz, Dimov,
  Sellier, et al.). The closest standard concepts are the signed-particle Wigner Monte Carlo with
  annihilation, and the Dirac-sea framing. If this work is heading toward external review, an
  explicit footnote on terminology would help readers who otherwise might think they should already
  know the term.

- **One spawning-slide phrase is misleadingly framed**, even though the math is right. The slide
  reads:
  > "Energy and momentum are not conserved: spawning only accompanies potential gradients!"

  Bare-particle momentum is not conserved, but momentum **is** conserved once the potential is
  treated as the source of an absorbed/emitted "photon" of momentum $h/L$ — exactly the
  perturbation-analysis identification made on the next slide. Suggested rewording:
  > "Bare-particle momentum is not conserved; the missing momentum is supplied by a quantum of
  > the potential field."

- **Worth checking the discrete implementation of $\tilde V_W = V_W - g(x)\thinspace\delta'(p)$.** The
  $\delta'(p)$ term is the formal way to express "subtract the classical force term." On a
  discrete grid it becomes a centered finite difference $[-1, 0, +1]/(2\Delta p)$, and sign
  conventions are easy to get wrong. The Maple worksheet `xFPjumpdensityV2.mw` (eqs. 21–30) is
  cited as the source — worth confirming the Python implementation matches.

- **The QLE → xFP match strictly requires a quasi-density.** Matching the second-$p$-derivative
  term to *zero* (the QLE has no $\partial^2 W/\partial p^2$) requires $\langle\zeta^2\rangle = 0$
  while $\langle\zeta^3\rangle \ne 0$, which is not achievable with any positive probability
  density. Either signed (quasi-) density or imaginary-axis support is required. This is the
  single point where the "stochastic" interpretation strains: the underlying jump variable is
  not a classical random variable.

---

## 5. References (within the project)

- `Extended_Fokker_Planck_Eq_and_the_QLE_V2.pdf` — primary analytical source.
- `Wigner_Collisions_Diagram_sozi.pdf` — visual / interpretive companion.
- Cited Maple worksheets: `xFPjumpdensityV2.mw`, `xFPjumpfromWignerPotential.mw`.
- Cited Jupyter notebook: `WignerSplitFourierV6.ipynb`.
- External: Wiedemann, *Particle Accelerator Physics* (Eq. 12.88–12.98); Kerr & Graham, "Generalized
  phase space version of Langevin equations and associated Fokker–Planck equations"; McKeon & Ord,
  *Phys. Rev. Lett.* **69** (1992); Reid & Drummond, *Objective QFT*.
