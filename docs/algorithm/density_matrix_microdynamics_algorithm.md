# Density-Matrix Microdynamics Algorithm

**A specification for simulating the von Neumann equation via a complex-weighted world-particle pair ensemble, as the position-space dual of the phase-space crystal-lattice algorithm.**

---

## 0. Status of this specification

This document is a **forward-looking, conjectural** specification, in the spirit of `docs/algorithm/multi_body_extension.md`. It develops the position-space dual of the phase-space crystal-lattice algorithm: where that algorithm represents the Wigner distribution $W(x, p, t)$ by signed particles undergoing deterministic advection plus discrete momentum jumps, this one represents the density matrix $\rho(x, x', t)$ in the position basis by complex-weighted world-particle **pairs** undergoing pair-Bohm flow plus phase accumulation on the weight. The two algorithms are Fourier duals in the $p \leftrightarrow s = x - x'$ variable; see §7 for the detailed correspondence.

The algorithmic structure described here was developed independently in the dissipative-quantum-systems literature, where it is known as **stochastic unraveling** of the Liouville–von Neumann equation. The closest direct counterpart is the **Stochastic Liouville–von Neumann (SLN) equation** of Stockburger & Grabert (*Phys. Rev. Lett.* 88, 170407, 2002), and its non-Lindblad-class extension by Kondov, Kleinekathöfer & Schreiber (*J. Chem. Phys.* 119, 6635, 2003), which makes essential use of trajectories with negative (signed) weights — exactly the positon/negaton structure of this project carried over to position-space pairs. The Bohm-style hydrodynamic decomposition of the velocity field (§3) is the standard one (Bohm 1952; reviewed for the density matrix in Burghardt & Parlant, *J. Chem. Phys.* 120, 3055, 2004).

What is original to this specification, relative to that literature, is (i) the framing of the algorithm as a fundamental subquantum microdynamics — the position-space dual of the phase-space crystal-lattice model — rather than as a numerical unraveling of an open-system master equation, and (ii) the explicit duality table (§7) relating it term-by-term to the Wigner-side algorithm.

**An important structural caveat.** The Wigner equation can be cast as a real-valued Liouville-plus-jumps process on real phase space; the von Neumann equation in position representation **cannot** be cast as a closed real-valued Itô SDE on $(X, X') \in \mathbb{R}^2$ with real drift and real noise alone. The kinetic operator $-i\hbar(\partial_x^2 - \partial_{x'}^2)/2m$ has an imaginary coefficient, which a real-drift real-noise SDE on pair coordinates cannot reproduce by itself. The two ways the existing literature handles this are:

- **(a) Pair-Bohm flow with phase-accumulating weights** (this document, §3–§4). Trajectories are real and deterministic, given by the pair-Bohm flow; complex weights along trajectories carry the imaginary content via a phase that accumulates the kinetic action, the potential difference, and the bilocal quantum potential.
- **(b) Stochastic state vectors rather than positions** (Stockburger–Grabert). Each "sample" is a pair of stochastic wave functions $(|\psi_m\rangle, |\phi_m\rangle)$ evolving by unitary Schrödinger; position pairs $(X_m, X'_m)$ are extracted from these wave functions only at observation time.

Adding real Nelson-style noise on top of route (a) is possible but requires additional weight-evolution terms that I have not derived in closed form here; see §8 for the open derivation.

Items below marked **[choice]** are implementation decisions not directly fixed by the underlying physics; alternatives are noted. Items marked **[open]** are derivations or claims that need further work.

---

## 1. Position representation: phase space and discretization

The state is the density matrix $\rho(x, x', t) = \langle x | \hat\rho(t) | x'\rangle$ in 1+1D, satisfying the von Neumann equation

$$i\hbar\thinspace\partial_t\rho \\; = \\; -\frac{\hbar^2}{2m}\bigl(\partial_x^2 - \partial_{x'}^2\bigr)\rho \\; + \\; \bigl[V(x) - V(x')\bigr]\thinspace\rho.$$

For comparison and Fourier-pair consistency with the phase-space crystal-lattice algorithm we use the same periodic position domain $x \in [0, L]$ and the same position grid

$$x_m = m\thinspace\Delta x, \qquad \Delta x = L/M, \qquad m = 0, \dots, M{-}1.$$

The off-diagonal coordinate runs over the same grid:

$$x'_{m'} = m'\thinspace\Delta x, \qquad m' = 0, \dots, M{-}1.$$

The empirical density matrix lives on the **$M \times M$** position-pair grid. Diagonal elements $\rho(x_m, x_m)$ are the (real, non-negative) position probability density; off-diagonal elements encode coherence and are complex in general, with the Hermiticity constraint $\rho(x, x') = \rho^*(x', x)$.

### 1.1 Useful change of variable

The Wigner-conjugate variables

$$r = (x + x')/2, \qquad s = x - x'$$

put (1.1) into

$$i\hbar\thinspace\partial_t\rho(r, s, t) \\; = \\; -\frac{\hbar^2}{m}\thinspace\partial_r\partial_s\thinspace\rho \\; + \\; \bigl[V(r + s/2) - V(r - s/2)\bigr]\thinspace\rho.$$

The Wigner function $W(r, p, t)$ is the inverse Fourier transform of $\rho(r, s, t)$ in $s \to p$. The two key structural facts visible in (1.2):

- The **kinetic** part of (1.2) is a *mixed second-derivative* (a hyperbolic operator) with imaginary coefficient — a complex wave-like dynamics in $(r, s)$, not the first-order advection that the Wigner equation has on $(r, p)$.
- The **potential** part of (1.2) is *multiplicative* in $(r, s)$ — its stochastic realization is a pointwise phase factor on each pair, not the momentum-jump operator that the Wigner equation has on $(r, p)$.

This is the structural inversion that motivates the dual algorithm.

---

## 2. State representation: world-particle pair ensemble

### 2.1 The ensemble

The state is represented by an ensemble of $N_w$ **world-particle pairs**:

$$\mathcal{E}(t) \\; = \\; \bigl\\{ (X_i(t),\ X'_i(t),\ w_i(t)) \bigr\\}_{i = 1}^{N_w}$$

where $X_i, X'_i \in [0, L]$ are continuous-valued positions (periodic) and $w_i \in \mathbb{C}$ is a complex weight. The empirical density matrix is

$$\rho_{\rm emp}(x, x', t) \\; = \\; \sum_{i = 1}^{N_w} w_i(t)\thinspace\delta\bigl(x - X_i(t)\bigr)\thinspace\delta\bigl(x' - X'_i(t)\bigr)$$

with the physical density matrix recovered as the ensemble average

$$\rho(x, x', t) \\; = \\; \mathbb{E}\bigl[\rho_{\rm emp}(x, x', t)\bigr].$$

For the deterministic pair-Bohm flow of §4 the expectation is trivial (a single deterministic ensemble reproduces $\rho$); for stochastic variants the expectation is over noise realizations.

### 2.2 Why complex weights are forced

In the Wigner-side algorithm, signed (real) weights suffice because $W$ is real. The density matrix is complex, with the Hermiticity constraint $\rho(x', x) = \rho^*(x, x')$, so the weights must be complex in general. The four-element discrete group $\\{+1, -1, +i, -i\\}$ (a `$\mathbb{Z}_4$` charge) suffices for a particle-counting implementation analogous to positon/negaton; for the mesh-density form (§4.5) the weights live in $\mathbb{C}$.

For Hermiticity preservation by construction one can enforce a **pair-twin symmetry**: every world-pair $(X, X', w)$ is generated together with its complex-conjugate partner $(X', X, w^*)$. With this symmetry, $\rho_{\rm emp}(x, x') = \rho_{\rm emp}^*(x', x)$ holds at the sample level, not just in expectation. **[choice]** — the asymmetric form is also valid; Hermiticity then holds only in expectation, as in Stockburger–Grabert.

### 2.3 Boundedness — and the missing "negaton crystal" analog

In the Wigner-side algorithm, $|W| \le 2/h$ uniformly, so a fixed-depth Dirac-sea / static negaton background suffices to make all counts non-negative. The density matrix has no such uniform bound: $\rho(x, x)$ is unbounded above, and $|\rho(x, x')| \le \sqrt{\rho(x, x)\rho(x', x')}$ is therefore likewise unbounded. The static-background trick of the Wigner-side algorithm has no direct analog here.

The available substitutes are: (a) complex weights with cancellation (the default here, and the Stockburger–Grabert route); (b) reference-state subtraction, where we represent $\delta\rho = \rho - \rho_0$ for some known stationary $\rho_0$. **[choice]**

---

## 3. Velocity field and weight phase from the ensemble state

The kinetic step of the algorithm (§4.1) moves each pair along the pair-Bohm flow, whose velocity field is extracted from the *current* empirical density matrix. The weight then accumulates a phase given by (3.7), which has three contributions: a classical Lagrangian, a potential difference, and a bilocal quantum potential. This section spells out how all of these are computed.

### 3.1 Modulus–phase decomposition of $\rho$

Write the density matrix in modulus–phase form

$$\rho(x, x', t) \\; = \\; R(x, x', t)\thinspace\exp\negthinspace\Bigl[\thinspace\tfrac{i}{\hbar}\thinspace\Phi(x, x', t)\thinspace\Bigr]$$

with $R \ge 0$ and $\Phi$ real. Hermiticity of $\rho$ gives

$$R(x, x') = R(x', x), \qquad \Phi(x, x') = -\Phi(x', x).$$

For a **pure state** $\rho = \psi\psi^*$ with $\psi = \sqrt{P}\thinspace e^{iS/\hbar}$, the decomposition reduces to

$$R(x, x') = \sqrt{P(x)\thinspace P(x')}, \qquad \Phi(x, x') = S(x) - S(x'),$$

so that $\Phi$ inherits its $x$ and $x'$ dependences from independent copies of the Bohm action $S$. For a **mixed state**, $\Phi(x, x')$ no longer factors, but the decomposition (3.1) is still well-defined wherever $R > 0$.

### 3.2 Pair-Bohm current velocities

The pair-Bohm velocities are defined from the phase:

$$v_+(x, x') \\; = \\; \\;\\;\thinspace\frac{1}{m}\thinspace\partial_x\thinspace\Phi(x, x'), \qquad v_-(x, x') \\; = \\; -\thinspace\frac{1}{m}\thinspace\partial_{x'}\thinspace\Phi(x, x').$$

For a pure state these reduce to $v_+(x, x') = S'(x)/m$ and $v_-(x, x') = S'(x')/m$ — each leg drifts at the standard one-sided Bohm velocity evaluated at *its own* position, independent of the other leg. Both have the same sign (both move forward at the local group velocity), as one would expect for a real-time-evolving pair where both members track $\psi$ in our wall-clock time.

For a mixed state, $v_\pm$ become genuinely bilocal: each leg's drift depends on both members of the pair.

### 3.3 The continuity equation along the pair-Bohm flow

Substituting (3.1) into (1.1) and separating real and imaginary parts yields:

- **Continuity-like equation (imaginary part):**

  $$\partial_t R^2 \\; + \\; \partial_x\negthinspace\bigl(R^2\thinspace v_+\bigr) \\; + \\; \partial_{x'}\negthinspace\bigl(R^2\thinspace v_-\bigr) \\; = \\; 0$$

  This is the statement that **the modulus-squared of the density matrix, $R^2 = |\rho|^2$, is transported by the pair-Bohm flow $(v_+, v_-)$**. It is the central justification for the deterministic algorithm of §4: trajectories generated by (3.3) carry $|\rho|^2$ along with them.

- **Hamilton–Jacobi-like equation (real part):**

  $$\partial_t \Phi \\; + \\; \frac{1}{2m}\bigl[(\partial_x \Phi)^2 - (\partial_{x'}\Phi)^2\bigr] \\; + \\; \bigl[V(x) - V(x')\bigr] \\; - \\; Q(x, x') \\; = \\; 0$$

  with the **pair quantum potential**

  $$Q(x, x') \\; = \\; \frac{\hbar^2}{2m}\thinspace\frac{1}{R}\bigl[\partial_x^2 R - \partial_{x'}^2 R\bigr].$$

  For a pure state this factors: defining the standard one-particle Bohm quantum potential $Q_\psi(x) = -(\hbar^2/2m)\thinspace\bigl(\partial_x^2\sqrt{P}\bigr)/\sqrt{P}$, one has $Q(x, x') = Q_\psi(x') - Q_\psi(x)$. (The sign and ordering follow from $\rho = \psi(x)\psi^*(x')$ and the relative $-\partial_{x'}^2$ in the von Neumann equation.) For a mixed state $Q$ is genuinely bilocal and does not factor.

The pair $(3.4, 3.5)$ is the density-matrix analog of the polar form of the Schrödinger equation in standard Bohmian mechanics. (3.4) says the flow transports $|\rho|^2$; (3.5) tells us how the phase $\Phi$ evolves, including along a Bohm trajectory.

### 3.4 Phase accumulation along a Bohm pair-trajectory

For a pair-trajectory $(X(t), X'(t))$ satisfying $\dot X = v_+(X, X')$, $\dot X' = v_-(X, X')$, the convective derivative of $\Phi$ is

$$\frac{d\Phi}{dt}\thinspace\bigg|_{X(t), X'(t)} \\; = \\; \partial_t\Phi \\; + \\; v_+\thinspace\partial_x\Phi \\; + \\; v_-\thinspace\partial_{x'}\Phi.$$

Substituting (3.5) and using the velocity definitions (3.3), which give $v_+\partial_x\Phi = m v_+^2$ and $v_-\partial_{x'}\Phi = -m v_-^2$:

$$\boxed{\\;\frac{d\Phi}{dt}\thinspace\bigg|_{X(t), X'(t)} \\; = \\; \frac{m}{2}\bigl(v_+^2 - v_-^2\bigr) \\; - \\; \bigl[V(X) - V(X')\bigr] \\; + \\; Q(X, X')\\;}$$

This is the formula for **how much phase the weight accumulates per unit time** as the pair moves along its Bohm trajectory. The three contributions:

- $(m/2)(v_+^2 - v_-^2)$ — the **classical pair Lagrangian**, the dual-leg analog of $(m/2)v^2$ for a single particle.
- $-(V(X) - V(X'))$ — the **potential contribution**, the same as a naïve pointwise multiplicative phase from $V(x) - V(x')$.
- $+Q(X, X')$ — the **pair quantum potential**, the bilocal version of Bohm's quantum potential.

Equation (3.7) is the most important equation in this specification. It tells us that the weight evolution is **not** just the pointwise phase from $V(X) - V(X')$ that one might naïvely expect from the multiplicative potential in (1.1) — it has additional contributions from kinetic action and quantum potential that arise because the trajectory is moving through phase-varying regions of $\rho$.

### 3.5 Operationally direct extraction from $\rho_{\rm emp}$

Computing $v_\pm$ and $Q$ from the empirical density matrix is most direct via the logarithmic-derivative form:

$$\frac{\partial_x\thinspace\rho}{\rho} \\; = \\; \frac{\partial_x R}{R} \\; + \\; \frac{i}{\hbar}\thinspace\partial_x\Phi.$$

Hence the **velocity field comes from the imaginary part:**

$$v_+(x, x') \\; = \\; \frac{\hbar}{m}\thinspace\mathrm{Im}\negthinspace\Biggl[\thinspace\frac{\partial_x\thinspace\rho(x, x')}{\rho(x, x')}\thinspace\Biggr], \qquad v_-(x, x') \\; = \\; -\thinspace\frac{\hbar}{m}\thinspace\mathrm{Im}\negthinspace\Biggl[\thinspace\frac{\partial_{x'}\thinspace\rho(x, x')}{\rho(x, x')}\thinspace\Biggr]$$

and the **modulus-derivative quantities (used for the quantum potential) come from the real part:**

$$\frac{\partial_x R}{R} \\; = \\; \mathrm{Re}\negthinspace\Biggl[\thinspace\frac{\partial_x\thinspace\rho}{\rho}\thinspace\Biggr], \qquad \frac{\partial_{x'} R}{R} \\; = \\; \mathrm{Re}\negthinspace\Biggl[\thinspace\frac{\partial_{x'}\thinspace\rho}{\rho}\thinspace\Biggr].$$

These give the pair quantum potential via

$$Q(x, x') \\; = \\; \frac{\hbar^2}{2m}\negthinspace\biggl[\thinspace\partial_x\negthinspace\biggl(\frac{\partial_x R}{R}\biggr) + \biggl(\frac{\partial_x R}{R}\biggr)^{\\!\\!2} \\; - \\; \partial_{x'}\negthinspace\biggl(\frac{\partial_{x'} R}{R}\biggr) - \biggl(\frac{\partial_{x'} R}{R}\biggr)^{\\!\\!2}\thinspace\biggr]$$

using the identity $\partial^2 R / R = \partial(\partial R/R) + (\partial R/R)^2$. **[choice]** — alternatively, compute $R = |\rho|$ on the grid directly and apply a standard second-difference stencil.

### 3.6 Extracting the velocity field and weight phase from the ensemble

Given the ensemble $\mathcal{E}(t) = \\{(X_i, X'_i, w_i)\\}$, the procedure each timestep is:

1. **Bin** the weighted pairs onto the $M \times M$ position-pair grid:

   $$\rho_{\rm bin}(x_m, x'_{m'}) \\; = \\; \frac{1}{\Delta x^2}\negthinspace\sum_{i \in \mathrm{bin}(m, m')}\negthinspace w_i.$$

2. **Smooth** $\rho_{\rm bin}$ with a 2D Gaussian kernel of width $\sigma_s \sim 2\thinspace\Delta x$. Necessary because the empirical $\rho$ is a sum of delta functions and (3.8)–(3.10) take derivatives of $\rho$. **[choice]**

3. **Differentiate** $\rho_{\rm bin}$ on the grid via FFT to get $\partial_x\rho_{\rm bin}$, $\partial_{x'}\rho_{\rm bin}$. Compute $v_\pm$ from (3.9) and $\partial_x R/R$, $\partial_{x'} R/R$ from (3.10).

4. **Regularise** near nodes of $\rho$: where $|\rho_{\rm bin}| < \epsilon_\rho\thinspace\max|\rho_{\rm bin}|$, set $v_\pm = 0$ and the modulus-derivative quantities to zero to avoid division-by-near-zero blowup. **[choice]**

5. **Build $Q$** on the grid from (3.11) by a second pass of finite differencing.

6. **Interpolate** $v_+, v_-, Q$ from the grid to each pair's location $(X_i, X'_i)$ by bilinear interpolation. **[choice]**

The cost is $\mathcal{O}(N_w + M^2 \log M)$ per timestep.

---

## 4. Time evolution

Each timestep $\Delta t$ applies the deterministic pair-Bohm flow plus a weight-phase update.

### 4.1 Pair-Bohm flow (kinetic step)

For each pair $i$:

$$X_i(t + \Delta t) \\; = \\; X_i(t) \\; + \\; v_+(X_i, X'_i)\thinspace\Delta t$$

$$X'_i(t + \Delta t) \\; = \\; X'_i(t) \\; + \\; v_-(X_i, X'_i)\thinspace\Delta t$$

with $v_\pm$ extracted from the current $\rho_{\rm emp}$ as in §3.6. This is **deterministic** flow given the current ensemble state; there is no noise in this step. The "stochasticity" in the algorithm enters only through the initial sampling of the ensemble — exactly as in the Wigner-side algorithm, where the deterministic Hamilton flow is also deterministic given the initial ensemble.

### 4.2 Phase-accumulation step (weight update)

For each pair $i$:

$$w_i(t + \Delta t) \\; = \\; w_i(t) \cdot \exp\negthinspace\Biggl[\thinspace\frac{i\thinspace\Delta t}{\hbar}\negthinspace\biggl(\thinspace\frac{m}{2}\bigl(v_+^2 - v_-^2\bigr) \\; - \\; \bigl(V(X_i) - V(X'_i)\bigr) \\; + \\; Q(X_i, X'_i)\thinspace\biggr)\thinspace\Biggr]$$

evaluated at the pair's position $(X_i, X'_i)$ at time $t$ (or, for higher accuracy, the time-average over the timestep). The three contributions are exactly the three terms of (3.7).

For comparison and as a sanity check: if one omitted the kinetic-action and quantum-potential terms, keeping only the potential phase

$$w_i \mapsto w_i \exp\negthinspace\Bigl[\thinspace-\thinspace i\thinspace\Delta t\thinspace\bigl(V(X_i) - V(X'_i)\bigr)/\hbar\thinspace\Bigr],$$

the algorithm would correctly evolve only the **static** part of $\rho$ where $v_\pm = 0$ and $Q = 0$. The full (4.3) is required for any non-stationary state.

### 4.3 Combined timestep (Lie split)

Concretely, one full timestep is:

1. Extract $v_\pm$ and $Q$ on the grid from $\rho_{\rm emp}(t)$ (§3.6).
2. Move each pair by (4.1)–(4.2). **[choice]** — for higher accuracy use a midpoint or RK2 integrator rather than Euler.
3. Update each weight by (4.3).

### 4.4 Strang split for second-order accuracy

For improved accuracy:

1. Compute $v_\pm, Q$ from $\rho_{\rm emp}(t)$.
2. Half-step pair-Bohm flow: (4.1)–(4.2) with $\Delta t \to \Delta t/2$.
3. Recompute $v_\pm, Q$ at the half-step ensemble.
4. Full weight update (4.3) using the half-step $v_\pm, Q$.
5. Half-step pair-Bohm flow with $\Delta t \to \Delta t/2$.

### 4.5 Equivalent mesh-density form

When $N_w$ is large and the goal is to reproduce $\rho$ on a fixed $M \times M$ grid (rather than to obtain individual world-particle trajectories), the mesh form is more efficient:

- Initialize $\rho$ on the grid directly.
- Each timestep, integrate the von Neumann equation (1.1) by a spectral method: split into kinetic (which diagonalises in the $(k_r, k_s)$ Fourier basis via (1.2)) and potential (which is pointwise multiplication by $\exp[-i\Delta t(V(x) - V(x'))/\hbar]$). Strang split.

This is the position-space analog of the Wigner-side split-Fourier algorithm and provides a deterministic reference for validating the particle-based form.

---

## 5. Initialization and observable extraction

### 5.1 Initialize

Two routes, depending on whether the initial state is given as a wave function or as a density matrix:

**(a) Pure state from $\psi_0(x)$:**

```python
# Construct rho_0(x, x') = psi_0(x) psi_0*(x') on the M x M grid
psi_0 = ...                                  # complex M-vector
rho_0 = np.outer(psi_0, psi_0.conj())        # M x M complex

# Importance-sample N_w pairs from |rho_0(x, x')|;
# set each weight to the phase of rho_0 at the sample point times
# the total normalisation so that sum(w_i) reconstructs the integral of rho_0.
```

**(b) Mixed state from a thermal density matrix or other $\rho_0$:**

```python
# rho_0(x, x') given directly on the M x M grid (complex, Hermitian)
# Importance-sample from |rho_0|; weights are the phase of rho_0
```

In both cases the initial weights have unit modulus (times a normalisation constant). For an $M \times M$ reconstruction grid at noise floor $\sigma_\rho$, target $N_w \gtrsim M^2 / \sigma_\rho^2$. **[choice]**

### 5.2 Run

Apply §4.3 or §4.4 for the desired number of timesteps.

### 5.3 Extract observables

```python
rho_recon = bin_pairs(X, X_prime, w, grid)       # M x M complex array
P_x       = rho_recon.diagonal().real            # position probability density
W         = ifft_in_s(rho_recon)                 # Wigner function via FFT in s = x - x'
```

The Wigner function is recovered by inverse Fourier transform along the $s = x - x'$ axis. This provides the cleanest validation against the phase-space crystal-lattice algorithm: both algorithms produce the same Wigner function from the same initial state, up to Monte-Carlo noise on this side and discretisation error on both.

---

## 6. Reference pseudocode

```python
# Inputs:
#   V(x)                      -- potential, callable on a numpy array
#   psi_0(x)                  -- initial wave function (pure state)
#   M                         -- position-grid size (same on both axes)
#   L                         -- periodic-domain length
#   N_w                       -- number of world-particle pairs
#   dt, num_steps             -- time-stepping
#   sigma_s                   -- smoothing scale for drift extraction (~ 2 dx)
#   eps_rho                   -- noise-floor regularisation for v_pm and Q

import numpy as np
from numpy.fft import fft2, ifft2

dx     = L / M
x_grid = np.arange(M) * dx
k_grid = 2 * np.pi * np.fft.fftfreq(M, d=dx)

# ---- 5.1 Initialisation: importance-sample pairs from |rho_0| ----------
psi_0_vec = psi_0(x_grid)                                # complex (M,)
rho_0     = np.outer(psi_0_vec, psi_0_vec.conj())        # (M, M) complex
abs_rho_0 = np.abs(rho_0)
norm      = abs_rho_0.sum() * dx**2
probs     = abs_rho_0.flatten() / abs_rho_0.sum()
idx       = np.random.choice(M*M, size=N_w, p=probs)
m_idx, mp_idx = np.unravel_index(idx, (M, M))

X  = x_grid[m_idx]  + np.random.uniform(-dx/2, dx/2, N_w)
Xp = x_grid[mp_idx] + np.random.uniform(-dx/2, dx/2, N_w)
w  = (rho_0[m_idx, mp_idx] / abs_rho_0[m_idx, mp_idx]) * (norm / N_w)


def empirical_rho_smoothed(X, Xp, w):
    """Bin weighted pairs onto the M x M grid, then 2D-Gaussian-smooth via FFT."""
    rho_bin = np.zeros((M, M), dtype=complex)
    mi  = np.floor(X  / dx).astype(int) % M
    mpi = np.floor(Xp / dx).astype(int) % M
    np.add.at(rho_bin, (mi, mpi), w / dx**2)
    K = np.exp(-0.5 * (k_grid[:, None]**2 + k_grid[None, :]**2) * sigma_s**2)
    return ifft2(fft2(rho_bin) * K)


def fields_from_rho(rho_s):
    """Compute v_+, v_-, and Q on the M x M grid from the smoothed rho."""
    rho_k    = fft2(rho_s)
    drho_dx  = ifft2(1j * k_grid[:, None] * rho_k)
    drho_dxp = ifft2(1j * k_grid[None, :] * rho_k)

    floor    = eps_rho * np.max(np.abs(rho_s))
    mask     = np.abs(rho_s) < floor
    denom    = np.where(mask, floor, rho_s)               # complex; safe divisor

    log_dx   = drho_dx  / denom                           # (d_x rho)/rho
    log_dxp  = drho_dxp / denom                           # (d_x' rho)/rho

    v_plus   = (hbar / mass) * np.imag(log_dx)
    v_minus  = -(hbar / mass) * np.imag(log_dxp)

    dR_R_x   = np.real(log_dx)                            # (d_x R)/R
    dR_R_xp  = np.real(log_dxp)                           # (d_x' R)/R

    # Q from (3.11): hbar^2/(2m) * [ d_x(dR/R) + (dR/R)^2 - d_x'(...) - (...)^2 ]
    dR_R_x_k   = fft2(dR_R_x)
    dR_R_xp_k  = fft2(dR_R_xp)
    d2_R_R_x   = np.real(ifft2(1j * k_grid[:, None] * dR_R_x_k))
    d2_R_R_xp  = np.real(ifft2(1j * k_grid[None, :] * dR_R_xp_k))
    Q          = (hbar**2 / (2 * mass)) * (
                  d2_R_R_x  + dR_R_x**2  - d2_R_R_xp - dR_R_xp**2)

    # Regularise: zero out drifts and Q where rho is below the noise floor
    v_plus[mask]  = 0.0
    v_minus[mask] = 0.0
    Q[mask]       = 0.0
    return v_plus, v_minus, Q


def interp_bilinear(field, X, Xp):
    """Sample `field` at continuous pair locations (X, Xp), periodic in both axes."""
    u  = X  / dx; v  = Xp / dx
    u0 = np.floor(u).astype(int) % M;  v0 = np.floor(v).astype(int) % M
    u1 = (u0 + 1) % M;                 v1 = (v0 + 1) % M
    fu = u - np.floor(u);              fv = v - np.floor(v)
    return ((1-fu)*(1-fv)*field[u0, v0] + fu*(1-fv)*field[u1, v0]
          + (1-fu)*fv*field[u0, v1]    + fu*fv*field[u1, v1])


# ---- Time-stepping loop (Lie split; Strang is a small refactor) -------
for step in range(num_steps):

    # ---- 4.1 Extract velocity field and Q from the current ensemble ----
    rho_s = empirical_rho_smoothed(X, Xp, w)
    v_p, v_m, Q = fields_from_rho(rho_s)

    # ---- 4.1 (continued) Pair-Bohm flow: deterministic update of (X, X') ---
    v_p_i = interp_bilinear(v_p, X, Xp)
    v_m_i = interp_bilinear(v_m, X, Xp)
    Q_i   = interp_bilinear(Q,   X, Xp)
    X_new  = (X  + v_p_i * dt) % L
    Xp_new = (Xp + v_m_i * dt) % L

    # ---- 4.3 Weight phase update from (3.7) ----------------------------
    L_pair = 0.5 * mass * (v_p_i**2 - v_m_i**2)
    dV     = V(X) - V(Xp)
    phase  = (L_pair - dV + Q_i) * dt / hbar
    w      = w * np.exp(1j * phase)

    X, Xp = X_new, Xp_new

# ---- Reconstruction ---------------------------------------------------
rho_final = empirical_rho_smoothed(X, Xp, w)
P_x       = rho_final.diagonal().real                     # position probability density
```

Two things to note about the pseudocode:

- All three terms of (3.7) appear in the `phase` update. Dropping any of them gives an incorrect result for non-stationary states. This is the most important correctness check when validating an implementation.
- The cost-dominant step is the FFT-based field extraction: two FFTs and two inverse FFTs per timestep, $\mathcal{O}(M^2 \log M)$. For $M = 128$ and $N_w = 10^6$, FFT and binning/interpolation are comparable.

---

## 7. Duality with the phase-space crystal-lattice algorithm

The two algorithms are Fourier duals of each other in the $p \leftrightarrow s = x - x'$ variable. The structural correspondence:

| | **Wigner-side: phase-space crystal-lattice** | **Density-matrix-side: world-particle pair ensemble** |
|---|---|---|
| State | $W(x, p, t)$ on $(x, p)$ grid | $\rho(x, x', t)$ on $(x, x')$ grid |
| Reality | $W \in \mathbb{R}$ | $\rho \in \mathbb{C}$, Hermitian |
| Bound | $\lvert W\rvert \le 2/h$ uniform | $\lvert\rho(x, x')\rvert \le \sqrt{\rho(x,x)\rho(x',x')}$, no uniform bound |
| Particles | Single particles at $(x, p)$ | Pair-particles at $(x, x')$ |
| Charges | Sign $\\{+1, -1\\}$ (positon / negaton) | Complex weight $w \in \mathbb{C}$, or discrete `$\mathbb{Z}_4$` |
| Position update | First-order advection: $x \to x + (p/m)\thinspace\Delta t$ | Pair-Bohm flow: $X \to X + v_+\Delta t$, $X' \to X' + v_-\Delta t$ |
| Position-update stochasticity | None (deterministic streaming, given the ensemble) | None (deterministic Bohm flow, given the ensemble) |
| Velocity comes from | Each particle's own $p$ coordinate (local, single-particle) | Phase gradient of empirical $\rho$ (non-local, mean-field) |
| Potential update | Discrete momentum jumps at rate $\Gamma_q$ from $V$'s Fourier modes | Multiplicative phase on weight from (3.7), including $V(X) - V(X')$ |
| Potential-update stochasticity | Stochastic (Poisson) | Deterministic (no noise in the phase update) |
| Additional weight terms from kinetic | None — kinetic is just streaming | Pair Lagrangian $(m/2)(v_+^2 - v_-^2)$ and quantum potential $Q$ |
| Static background | $-2/h$ negaton crystal lattice | No direct analog; see §2.3 |
| Mediator | Local positon density mediates jumps | Local empirical $\rho$ shapes both the flow and the quantum potential |
| Reconstruction | $W = N_+/(\nu\thinspace\Delta x\thinspace\Delta p) - 2/h$ | Bin pairs, smooth, take diagonal for $P(x)$; FFT in $s$ for $W$ |

The roles of "stochastic" and "deterministic" swap between the kinetic and potential steps, just as the source operators in the underlying PDEs swap their structural roles. On the Wigner side the potential is the non-local operator (driving stochastic jumps in $p$) and the kinetic is local advection; on the density-matrix side the kinetic is the non-local operator (driving the bilocal pair-Bohm flow) and the potential is local (pointwise phase on each pair).

### 7.1 The asymmetries that the duality does not eliminate

Several asymmetries remain:

1. **The "single rule" structure.** The Wigner-side algorithm has *one* rule per kind of update — free streaming in $x$, mediated jumps in $p$. The density-matrix-side algorithm has one rule for position update (pair-Bohm flow) but the *weight* update has three contributions (kinetic Lagrangian, potential difference, quantum potential). It is not a single-rule mediated-jump scheme.
2. **Non-locality of velocity extraction.** Wigner-side velocity for each particle is just that particle's own $p$ value; density-matrix-side velocity is a *mean-field* read from the empirical $\rho$. This is the same mean-field structure that arises in classical Vlasov simulations.
3. **No negaton crystal.** As noted, the uniform Wigner bound $|W| \le 2/h$ has no analog on the density-matrix side, and the static-Dirac-sea trick that eliminates explicit negaton bookkeeping in the Wigner case has no direct counterpart here.
4. **Sign / phase problem.** Wigner-side uses signed (real) weights; density-matrix-side uses complex weights. The sample variance accumulates as $|w|^2 / \langle w\rangle^2$ rather than as $w / \langle w\rangle$, giving a tougher version of the sign problem familiar from fermion Monte Carlo. The variance-reduction techniques of Schmitz & Stockburger (2019) are relevant here.

---

## 8. Open implementation questions

1. **Time-step selection.** Two constraints: (i) the pair-Bohm flow wants $|v_\pm|\thinspace\Delta t \ll \sigma_s$ so trajectories don't cross multiple smoothing cells per step, and (ii) the phase update wants $|d\Phi/dt|\thinspace\Delta t/\hbar \ll 1$ to keep phase rotation small. For smooth potentials the latter is typically tighter near nodes of $\rho$, where $Q$ blows up.

2. **Drift and quantum-potential regularisation.** Equations (3.9)–(3.11) all involve $1/\rho$ (or $1/R$). Near nodes of $\rho$ they blow up. The simple regularisation $|\rho| < \epsilon_\rho \Rightarrow v_\pm = Q = 0$ caps the divergences but introduces systematic bias near nodes. A better approach is to detect low-occupancy cells and disable updates there entirely. **[open]** — needs empirical study on a non-trivial test case.

3. **Pair-twin symmetrisation.** Generating each pair together with its conjugate twin guarantees Hermiticity at the sample level but doubles the ensemble size. For most observables the asymmetric form gives Hermiticity on average and is sufficient. **[choice]**

4. **Resampling / population control.** Over time, the weight distribution $\\{|w_i|\\}$ broadens — a small number of pairs end up dominating, others having near-zero weight. Periodic resampling (duplicate high-weight, merge low-weight) is the standard fix. The Schmitz–Stockburger (2019) convex-optimisation variance-reduction approach is more principled.

5. **Adding Nelson noise.** Replacing the deterministic flow (4.1)–(4.2) with a Nelson-style real Itô SDE
$$dX = b_+\thinspace dt + \sqrt{\hbar/m}\thinspace dW_X, \qquad dX' = b_-\thinspace dt + \sqrt{\hbar/m}\thinspace dW_{X'},$$
with $b_\pm = v_\pm + u_\pm$ (current plus osmotic) and $u_\pm = \pm(\hbar/m)(\partial_{x,x'} R/R)$, requires additional terms in the weight phase to compensate the diffusion. Derivation of those compensating terms is **[open]**. The expected form is $\Delta\phi \supset -i(\sigma^2/2)[\partial_x^2 R/R - \partial_{x'}^2 R/R]\thinspace\Delta t$ or similar — closely related to $Q$. The motivation for adding noise is multi-worlds interpretability (each Wiener realisation = one world) and ergodicity across nodal lines.

6. **Negaton-background analog.** Whether there is some position-space analog of the static-Dirac-sea trick that uses the structure of $\rho$ to make all weights have a fixed sign or phase is an open question. The reference-state-subtraction route of §2.3(b) is one direction; another is to work with $\rho - \rho_\beta$ where $\rho_\beta$ is the thermal density matrix at some inverse temperature $\beta$. **[open]**

7. **Validation against the phase-space crystal-lattice algorithm.** A clean numerical test: take the cosine-well demo (`demo_cosine_well_microdynamics.py`), run both algorithms from the same initial pure-state Gaussian, and check that the reconstructed Wigner functions (from binning + Fourier-in-$s$ on this side; from direct $W$ on the Wigner side) agree to within Monte-Carlo noise. This is the recommended first validation; see `demo_density_matrix_microdynamics.py` (forthcoming).

---

## 9. Sources

### Primary structural correspondence

- **Stockburger, J. T.; Grabert, H.** — "Exact c-number representation of non-Markovian quantum dissipation." *Phys. Rev. Lett.* **88**, 170407 (2002). The stochastic Liouville–von Neumann (SLN) equation; closest direct counterpart in the literature.
- **Kondov, I.; Kleinekathöfer, U.; Schreiber, M.** — "Stochastic unraveling of Redfield master equations and its application to electron transfer problems." *J. Chem. Phys.* **119**, 6635 (2003); arXiv:physics/0307050. Pair-vector unraveling with explicitly *signed* (signed-weight) trajectories, mirroring the positon/negaton structure.
- **Schmitz, K.; Stockburger, J. T.** — "A variance-reduction technique for the stochastic Liouville–von Neumann equation." *Eur. Phys. J. Special Topics* **227**, 1929 (2019); arXiv:1812.05960.

### Bohmian mechanics for the density matrix

- **Bohm, D.** — "A suggested interpretation of the quantum theory in terms of hidden variables." *Phys. Rev.* **85**, 166 (1952). Original Bohmian decomposition for $\psi$.
- **Burghardt, I.; Parlant, G.** — "On the quantum hydrodynamic description of the density matrix." *J. Chem. Phys.* **120**, 3055 (2004). Density-matrix Bohmian hydrodynamics; the pair quantum potential of §3.

### Nelson stochastic mechanics (for the proposed §8 noise extension)

- **Nelson, E.** — *Dynamical Theories of Brownian Motion*, Princeton (1967).
- **Nelson, E.** — "Review of stochastic mechanics." *J. Phys.: Conf. Ser.* **361**, 012011 (2012).
- **Bacciagaluppi, G.** — "A Conceptual Introduction to Nelson's Mechanics." In *Endophysics, Time, Quantum and the Subjective*, World Scientific (2005). hal-shs:00996258.

### Double-path-integral scaffolding

- **Feynman, R. P.; Vernon, F. L.** — "The theory of a general quantum system interacting with a linear dissipative system." *Ann. Phys.* (NY) **24**, 118 (1963). The original double-path integral for $\rho(q, q', t)$.
- **Schwinger, J.** (1961) and **Keldysh, L. V.** (1965) — closed-time-path / in-in formalism, the field-theoretic version of the doubled-coordinate structure used here.

### Forward–backward trajectory and related methods

- **Hsieh, C.-Y.; Kapral, R.** — "Nonadiabatic dynamics in open quantum-classical systems: forward–backward trajectory solution." *J. Chem. Phys.* **137**, 22A507 (2012); arXiv:1302.2085.

### Many-worlds / subquantum framing (already in the project bibliography)

- **Smolin, L.** — "Could quantum mechanics be an approximation to another theory?" arXiv:quant-ph/0609109. Many-worlds reformulation of Nelson with discontinuous phase-space paths.

### Companion documents in this repository

- `docs/algorithm/phase_space_crystal_lattice_algorithm.md` — the Wigner-side algorithm whose dual is described here. §7 above is keyed term-by-term to that document.
- `docs/algorithm/multi_body_extension.md` — the multi-body generalisation of the Wigner-side algorithm. A similar extension is possible here.
- `docs/supplement/phase_space_crystal_lattice_supplement.md` — derivation of the Wigner-side jump rule from the Moyal series.
