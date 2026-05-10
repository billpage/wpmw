# Multi-Body Extension of the Phase-Space Crystal-Lattice Algorithm

**Specification for extending the single-particle 1+1D algorithm to $d$ spatial dimensions and to $N$-particle systems.**

---

## 0. Status of this specification

This document is a **forward-looking extension** of `docs/algorithm/phase_space_crystal_lattice_algorithm.md` (hereafter "the 1+1D spec"). It generalizes that algorithm in two stages:

1. **Single particle in $d$ spatial dimensions** (§§2–3). A direct vectorial extension; the structure of the 1+1D algorithm carries over with no conceptual change. Cost scales as $(M_x M_p)^d$.

2. **$N$ interacting particles in $d$ spatial dimensions** (§§4–7). The mesh-density representation becomes infeasible at $\mathcal{O}((M_x M_p)^{dN})$ phase-space cells; the **world-ensemble** representation (a stochastic sample of the joint Wigner distribution by points in the full $2dN$-dimensional phase space) becomes the only viable form. Pair interactions produce *correlated, momentum-conserving* jumps between particles within each world.

The 1+1D algorithm has been validated against direct mesh evolution in `demo_cosine_well_microdynamics.py`. The single-particle 3+3D extension is mechanical and should be implementable from this document; the $N$-body extension is more speculative — it inherits the curse of dimensionality, the fermionic sign problem, and the statistical-convergence concerns of any joint-distribution Monte Carlo, and §9 records open questions accordingly.

Items below marked **[choice]** are implementation decisions not directly fixed by the underlying physics; alternatives are noted.

---

## 1. Notation conventions

The 1+1D spec uses $M$ for position bins and $N$ for momentum bins. To free up $N$ for the particle count, this document uses:

| Symbol | Meaning | 1+1D spec equivalent |
|---|---|---|
| $d$ | spatial dimension (1, 2, 3, …) | implicitly 1 |
| $N$ | number of particles | (single particle) |
| $M_x$ | position bins per spatial axis | $M$ |
| $M_p$ | momentum bins per spatial axis | $N$ |
| $\vec x_i$ | position vector of particle $i$ in $\mathbb{R}^d$ | $x$ |
| $\vec p_i$ | momentum vector of particle $i$ in $\mathbb{R}^d$ | $p$ |
| $\vec q$ | Fourier-mode index vector ($\in \mathbb{Z}^d$) | $q$ |
| $\vec k_{\vec q}$ | wavevector $(2\pi q_1/L_1, \dots, 2\pi q_d/L_d)$ | $k_q$ |

The component of $\vec x_i$ along axis $a$ is written $x_{i,a}$, and similarly $p_{i,a}$, $k_{\vec q, a}$.

A periodic position box has side lengths $L_1, \dots, L_d$ (often equal). The natural per-axis momentum cell size is $\Delta p_a = \pi\hbar/L_a$, generalizing the 1+1D choice.

---

## 2. Single particle in $d$ dimensions: phase space and discretization

Phase space is $2d$-dimensional. With periodic position domain $\vec x \in \prod_a [0, L_a]$ and momentum domain $\vec p \in \prod_a [-P_a^{\max}, +P_a^{\max}]$:

- Position lattice: $x_{m,a} = m_a \Delta x_a$, $\Delta x_a = L_a/M_x$, $m_a \in \lbrace0, \dots, M_x{-}1\rbrace$
- Momentum lattice: $p_{n,a} = (n_a - M_p/2)\Delta p_a$, $n_a \in \lbrace0, \dots, M_p{-}1\rbrace$

Cell volume: $\Delta\Omega = \prod_a \Delta x_a \Delta p_a$. Total cell count: $(M_x M_p)^d$.

The natural momentum spacing $\Delta p_a = \pi\hbar/L_a$ **[choice]** ensures that a Fourier mode of $V$ with index $\vec q$ produces an exact integer-cell jump of $q_a$ cells along axis $a$ (in each dimension independently).

---

## 3. Single particle: state representation and time evolution

### 3.1 Crystal-lattice shift

The Wigner-function bound in $d$ dimensions is $|W(\vec x, \vec p)| \le (1/\pi\hbar)^d = (2/h)^d$. The crystal-lattice shift becomes

$$W'(\vec x, \vec p) \\; = \\; W(\vec x, \vec p) + \left(\tfrac{2}{h}\right)^d, \qquad W' \in \left[0,\ 2\negthinspace\left(\tfrac{2}{h}\right)^{\negthinspace d}\right].$$

Per-cell positon counts:

$$N_+(\vec m, \vec n) \\; =\\;  \nu \cdot W'(\vec x_{\vec m}, \vec p_{\vec n}) \cdot \Delta\Omega.$$

The negaton background lattice — every cell at notional count $\nu (2/h)^d \Delta\Omega$ — remains implicit and is never updated. Only $N_+$ evolves.

### 3.2 Free streaming (advection)

QLE term: $-\sum_a (p_a/m)\thinspace\partial_{x_a}W$. On the lattice, axis-by-axis row shift:

$$\Delta m_{\vec n, a}  \\; = \\;  \mathrm{round}\negthinspace\left(\frac{p_{\vec n, a}\thinspace\Delta t}{m\thinspace\Delta x_a}\right),$$

applied independently to each axis (the streaming sub-operators along different axes commute, since they involve different position coordinates). For zero advection error per axis, choose $\Delta t$ such that the maximum shift along the slowest-streaming axis is an exact integer.

### 3.3 Potential-driven mediated jumps (vectorial)

Decompose the potential into $d$-dimensional Fourier components on the periodic box:

$$V(\vec x) = V_0 + \sum_{\vec q \neq 0} V_{\vec q}\thinspace\cos(\vec k_{\vec q}\cdot \vec x + \phi_{\vec q}), \qquad \vec k_{\vec q,a} = \frac{2\pi q_a}{L_a}.$$

For each mode $\vec q$, define the local rate

$$\Gamma_{\vec q}(\vec x_{\vec m})  \\; = \\;  -\frac{V_{\vec q}}{\hbar}\sin(\vec k_{\vec q}\cdot \vec x_{\vec m} + \phi_{\vec q}).$$

**The single rule.** A positon at cell $(\vec m, \vec n)$ acts as a mediator: with probability $|\Gamma_{\vec q}(\vec x_{\vec m})|\thinspace\Delta t$ per particle per mode per timestep **[choice — Poisson rate is the rigorous form]**, it induces

$$\Gamma_{\vec q}(\vec x_{\vec m}) > 0:\quad (\vec m,\thinspace \vec n + \vec q) \longrightarrow (\vec m,\thinspace \vec n - \vec q)$$

and the opposite when $\Gamma_{\vec q} < 0$. The mediator itself is unchanged. The displacement $2\vec q$ in momentum-cell units corresponds to a momentum kick of $\hbar \vec k_{\vec q}$, i.e., the photon momentum at wavevector $\vec k_{\vec q}$ — now with definite *direction* set by $\hat k_{\vec q}$.

### 3.4 Equivalent mesh-density form

For large $\nu$ or when evolving $W$ directly:

$$W(\vec x_{\vec m}, \vec p_{\vec n}, t + \Delta t) = W(\vec x_{\vec m}, \vec p_{\vec n}, t) + \Delta t \sum_{\vec q \neq 0} \Gamma_{\vec q}(\vec x_{\vec m})\bigl[W(\vec x_{\vec m}, \vec p_{\vec n + \vec q}) - W(\vec x_{\vec m}, \vec p_{\vec n - \vec q})\bigr]$$

Continuum limit reproduces $\partial_t W = +\nabla_{\negthinspace\vec x} V(\vec x)\cdot\nabla_{\negthinspace\vec p} W$, the multi-dimensional QLE force term.

### 3.5 Coulomb / long-range potentials: Ewald split

For $V(\vec r) = -e^2/r$ on a periodic box, direct Fourier expansion has $V_{\vec q} \sim 1/|\vec k_{\vec q}|^2$, with a divergent $\vec q = 0$ mode. The Ewald split

$$\frac{1}{r} = \underbrace{\frac{\mathrm{erfc}(\alpha r)}{r}}_{\text{short-range}} + \underbrace{\frac{\mathrm{erf}(\alpha r)}{r}}_{\text{long-range}}$$

yields:

- $V^{\text{short}}$: exponentially confined; handle as a real-space deterministic momentum drift $\Delta\vec p = -\nabla V^{\text{short}}(\vec x)\thinspace\Delta t$ within a cutoff radius. The non-classical Moyal corrections to this drift are $O(\hbar^2 \alpha^2)$-small if $\alpha$ is chosen to keep $V^{\text{short}}$ smooth on the de Broglie scale.
- $V^{\text{long}}$: smooth at the origin; its Fourier coefficients decay as $\exp(-|\vec k|^2/4\alpha^2)/|\vec k|^2$. A small finite set of Fourier modes suffices — apply the §3.3 rule to each.

The neutralizing background ($\vec q = 0$ divergence) is absorbed by the requirement that the system carry zero total charge, exactly as in standard particle-mesh Ewald (PME).

### 3.6 Single-particle 3+3D pseudocode

```python
# Inputs: V_q dict keyed by tuple(q1,q2,q3); W0(x,p); M_x, M_p, nu, dt, num_steps
import numpy as np

D = 3
dx = np.array([L[a] / M_x for a in range(D)])
dp = np.array([np.pi * hbar / L[a] for a in range(D)])

# Position and momentum grids (D-dimensional)
x_grids = [np.arange(M_x) * dx[a] for a in range(D)]
p_grids = [(np.arange(M_p) - M_p // 2) * dp[a] for a in range(D)]

# N_plus has shape (M_x,)*D + (M_p,)*D — i.e., 2D-dim array
shape = (M_x,) * D + (M_p,) * D
W0_grid = compute_initial_wigner(x_grids, p_grids)        # 2D-dim array
N_plus = np.round((W0_grid + (2/h)**D) * nu * np.prod(dx) * np.prod(dp)).astype(int)

# Pre-compute Gamma_q(x) on the position grid for each Fourier mode
Gammas = {}
for q_tuple, V_q in V_q_dict.items():
    k = np.array([2*np.pi*q_tuple[a]/L[a] for a in range(D)])
    # broadcast to D-dimensional position grid
    XX = np.meshgrid(*x_grids, indexing='ij')
    phase = sum(k[a] * XX[a] for a in range(D)) + phi_q_dict[q_tuple]
    Gammas[q_tuple] = -(V_q / hbar) * np.sin(phase)       # shape (M_x,)*D

for step in range(num_steps):
    # ---- 3.2 Free streaming, axis by axis ------------------------------
    for a in range(D):
        # For each momentum index along axis a, shift along position axis a
        # (vectorized: pick out p_grid[a] index, compute integer shift)
        for n_a in range(M_p):
            shift = int(round(p_grids[a][n_a] * dt / (mass * dx[a])))
            if shift:
                # roll only the slab where momentum-axis-a index == n_a
                idx = (slice(None),) * D + tuple(
                    n_a if b == a else slice(None) for b in range(D)
                )
                N_plus[idx] = np.roll(N_plus[idx], shift, axis=a)

    # ---- 3.3 Mediated jumps -------------------------------------------
    # Mesh-density form (use stochastic per-particle form for sparse W)
    W = N_plus / (nu * np.prod(dx) * np.prod(dp)) - (2/h)**D
    dW = np.zeros_like(W)
    for q_tuple, Gamma_field in Gammas.items():
        # Roll along all D momentum axes by the q vector
        plus_axes  = tuple(D + a for a in range(D))
        plus_shift = q_tuple
        W_plus  = np.roll(W, plus_shift,  axis=plus_axes)
        W_minus = np.roll(W, tuple(-q for q in plus_shift), axis=plus_axes)
        # Broadcast Gamma_field (shape (M_x,)*D) over the momentum axes
        Gamma_b = Gamma_field[(...,) + (None,) * D]
        dW += dt * Gamma_b * (W_plus - W_minus)
    W += dW
    N_plus = np.round((W + (2/h)**D) * nu * np.prod(dx) * np.prod(dp)).astype(int)
```

---

## 4. $N$-body: phase space and joint Wigner function

For $N$ distinguishable particles in $d$ dimensions, the joint phase space is $2dN$-dimensional:

$$\Omega^{(N)} = \bigl\lbrace(\vec x_1, \vec p_1, \dots, \vec x_N, \vec p_N)\bigr\rbrace.$$

The joint Wigner function $W^{(N)}$ is bounded by

$$|W^{(N)}(\vec x_1, \vec p_1, \dots, \vec x_N, \vec p_N)| \le \left(\tfrac{2}{h}\right)^{dN},$$

the natural $dN$-dimensional generalization of the Hudson-bound family. The crystal-lattice shift becomes

$$W'^{(N)}  \\; = \\;  W^{(N)} + \left(\tfrac{2}{h}\right)^{dN}, \qquad W'^{(N)} \in \left[0,\ 2\negthinspace\left(\tfrac{2}{h}\right)^{\negthinspace dN}\right].$$

A direct mesh representation requires $(M_x M_p)^{dN}$ cells. For the simulation parameters of `demo_cosine_well_microdynamics.py` ($M_x = 64, M_p = 64$, $d = 1$) extended to $d = 3, N = 10$, this is $4096^{30} \approx 10^{108}$ cells. The mesh form is unusable for $N \gtrsim 2$.

---

## 5. $N$-body state representation: the world ensemble

### 5.1 Worlds as joint-positon samples

Represent the joint shifted distribution $W'^{(N)}$ by an ensemble of $\mathcal{W}$ independent samples — **worlds** — each world $\alpha$ being a single point in the $2dN$-dimensional phase space:

$$\omega^{(\alpha)} = \bigl(\vec x_1^{(\alpha)}, \vec p_1^{(\alpha)}, \dots, \vec x_N^{(\alpha)}, \vec p_N^{(\alpha)}\bigr) \in \Omega^{(N)}, \qquad \alpha = 1, \dots, \mathcal{W}.$$

A world is a complete instantaneous specification of the entire $N$-particle system: positions and momenta of all $N$ particles. The ensemble approximates the joint distribution as

$$W'^{(N)}(\vec x_1, \dots, \vec p_N) \approx \frac{1}{\nu_{\negthinspace\mathcal{W}}} \sum_{\alpha=1}^{\mathcal{W}} \delta\bigl(\vec x_1 - \vec x_1^{(\alpha)}\bigr) \cdots \delta\bigl(\vec p_N - \vec p_N^{(\alpha)}\bigr),$$

where $\nu_{\negthinspace\mathcal{W}}$ is the chosen worlds-per-unit-phase-space-volume scale (the multi-body analog of $\nu$ in the 1+1D spec).

### 5.2 Why this avoids the mesh blowup

Each world is $2dN$ floating-point numbers; total storage is $\mathcal{O}(\mathcal{W} \cdot dN)$ rather than $(M_x M_p)^{dN}$. The number of worlds $\mathcal{W}$ needed for a target statistical accuracy on a particular observable is *not* exponential in $N$ for one- and two-body marginals (which most laboratory observables reduce to); it is only the mesh-discretization cost that explodes with $N$. The world ensemble is the multi-body crystal-lattice algorithm's analog of DMC walkers, AFQMC configurations, and signed-particle Wigner Monte Carlo configurations.

The negaton background $W_0 = (2/h)^{dN}$ is implicit and fully shared across the ensemble — no negaton worlds are stored. The "excess" worlds $\omega^{(\alpha)}$ are precisely the laboratory-observable carriers (V2 spec §5).

### 5.3 Reduced-density extraction

The single-particle marginal Wigner function for particle $i$ is the integral over the other $N-1$ particles:

$$w_i(\vec x, \vec p) = \int W^{(N)}\thinspace\prod_{j \neq i} d\vec x_j\thinspace d\vec p_j.$$

In the world ensemble, this is the histogram of $(\vec x_i^{(\alpha)}, \vec p_i^{(\alpha)})$ over $\alpha$ — cheap, $\mathcal{O}(\mathcal{W})$. Two-body marginals $w_{ij}^{(2)}$ require histogramming pairs $(\vec x_i^{(\alpha)}, \vec p_i^{(\alpha)}, \vec x_j^{(\alpha)}, \vec p_j^{(\alpha)})$, similarly cheap. $k$-body correlations need at least $\mathcal{W} \gtrsim$ $\text{(fineness)}^{2dk}$ samples, so practical $k$ is bounded by available statistics.

---

## 6. $N$-body split-operator time evolution

Each timestep $\Delta t$ applies, for every world $\alpha$ in the ensemble, three sub-operators in sequence. Strang-split as in the 1+1D spec for second-order accuracy **[choice]**.

### 6.1 Free streaming (per particle, per world)

For each world $\alpha$ and each particle $i$:

$$\vec x_i^{(\alpha)} \longrightarrow \vec x_i^{(\alpha)} + \frac{\vec p_i^{(\alpha)}}{m_i}\thinspace\Delta t \quad (\bmod\ \vec L).$$

Independent across particles and worlds; trivially parallel. In the world-ensemble representation there is no lattice-snap step — particle positions and momenta are continuous variables stored per world.

### 6.2 One-body external potential

For an external potential $V_{\text{ext}}(\vec x)$ acting on each particle independently, decompose as in §3.3:

$$V_{\text{ext}}(\vec x) = V_0 + \sum_{\vec q} V^{\text{ext}}_{\vec q}\thinspace\cos(\vec k_{\vec q}\cdot \vec x + \phi^{\text{ext}}_{\vec q}).$$

For each world $\alpha$, each particle $i$, each Fourier mode $\vec q$, with rate

$$\Gamma^{\text{ext}}_{\vec q}(\vec x_i^{(\alpha)}) = -\frac{V^{\text{ext}}_{\vec q}}{\hbar}\sin(\vec k_{\vec q}\cdot \vec x_i^{(\alpha)} + \phi^{\text{ext}}_{\vec q}),$$

trigger a momentum jump of particle $i$:

$$\vec p_i^{(\alpha)} \longrightarrow \vec p_i^{(\alpha)} + \mathrm{sgn}(\Gamma^{\text{ext}}_{\vec q})\cdot\hbar\vec k_{\vec q}$$

with probability $|\Gamma^{\text{ext}}_{\vec q}|\thinspace\Delta t$ **[choice — Poisson rate is the rigorous form]**.

Note: in the world-ensemble form, momentum kicks are by the *full* photon momentum $\hbar\vec k_{\vec q}$, not split into half-jumps. The half-jump structure of the 1+1D mesh form arose from the symmetric source-to-destination split $(n+q) \to (n-q)$ in cell coordinates; in the continuous-variable world ensemble there is no such artifact, and a single jump of size $2 \cdot \tfrac{1}{2}\hbar k = \hbar k$ is delivered per event. The factor of 2 in the rate convention is absorbed accordingly **[choice — see §11 for an alternative half-jump convention preserving the mediator-pair picture]**.

### 6.3 Two-body interactions: correlated, momentum-conserving jumps

For a pair potential $V_2(\vec r_{ij})$ depending on the relative coordinate $\vec r_{ij} = \vec x_i - \vec x_j$, decompose as a Fourier series on the periodic box:

$$V_2(\vec r) = V_2^{(0)} + \sum_{\vec q} V^{(2)}_{\vec q}\thinspace\cos(\vec k_{\vec q}\cdot \vec r + \phi^{(2)}_{\vec q}).$$

Because $\cos(\vec k\cdot \vec r_{ij}) = \cos(\vec k\cdot \vec x_i - \vec k\cdot \vec x_j)$ couples *only* particles $i$ and $j$ (and only via the difference of their positions), the corresponding Moyal-bracket contribution to $\partial_t W^{(N)}$ acts as a **simultaneous, opposite-sign shift in $\vec p_i$ and $\vec p_j$** (derivation in §7).

**The two-body single rule.** For each world $\alpha$, each unordered pair $(i, j)$, each Fourier mode $\vec q$ of $V_2$, with rate

$$\Gamma^{(2)}_{\vec q}(\vec r_{ij}^{(\alpha)}) = -\frac{V^{(2)}_{\vec q}}{\hbar}\sin(\vec k_{\vec q}\cdot \vec r_{ij}^{(\alpha)} + \phi^{(2)}_{\vec q}),$$

trigger the correlated jump

$$\bigl(\vec p_i^{(\alpha)},\thinspace \vec p_j^{(\alpha)}\bigr) \longrightarrow \bigl(\vec p_i^{(\alpha)} + \mathrm{sgn}(\Gamma^{(2)}_{\vec q})\thinspace\hbar\vec k_{\vec q},\ \vec p_j^{(\alpha)} - \mathrm{sgn}(\Gamma^{(2)}_{\vec q})\thinspace\hbar\vec k_{\vec q}\bigr)$$

with probability $|\Gamma^{(2)}_{\vec q}|\thinspace\Delta t$. The total momentum $\vec p_i + \vec p_j$ is conserved exactly per event. This is the entire two-body crystal-lattice rule.

**Physical interpretation.** Each correlated jump is an exchange of a "virtual photon" of momentum $\hbar\vec k_{\vec q}$ between particles $i$ and $j$ at relative position $\vec r_{ij}$. The rate-amplitude $|V^{(2)}_{\vec q}|/\hbar$ is the coupling at that Fourier component of the interaction. This is structurally the leading Feynman-diagram element of QED at the non-relativistic-QLE level; the algorithm is a stochastic implementation of pair-wise quantum exchange.

### 6.4 Higher-body interactions

For genuinely $k$-body potentials $V_k(\vec x_{i_1}, \dots, \vec x_{i_k})$ (rare in physics; common in effective theories), the same Fourier decomposition extends: a Fourier mode involving wavevectors $\vec k^{(1)}, \dots, \vec k^{(k)}$ summing to zero (for translation invariance) drives a $k$-particle correlated momentum exchange. Cost per step grows as $\binom{N}{k}$, and the rate-rare-event statistics demand more worlds. **[choice]** — most implementations will truncate at $k = 2$, leaving higher-body corrections to be either absorbed into effective two-body potentials or handled as a Moyal-series correction in differential mesh form on reduced marginals.

---

## 7. Two-body interaction: derivation

The Moyal-bracket contribution of $V_2(\vec r_{ij})$ to the joint QLE is

$$\partial_t W^{(N)} \supset \frac{2}{\hbar}\thinspace V_2(\vec r_{ij})\thinspace\sin\negthinspace\left(\frac{\hbar}{2}\bigl(\overleftarrow{\nabla}_{\negthinspace\vec x_i} - \overleftarrow{\nabla}_{\negthinspace\vec x_j}\bigr)\cdot\bigl(\overrightarrow{\nabla}_{\negthinspace\vec p_i} - \overrightarrow{\nabla}_{\negthinspace\vec p_j}\bigr)\right) W^{(N)}.$$

(Note the sign structure: only the *difference* of position-derivatives appears, because $V_2$ depends on $\vec r_{ij} = \vec x_i - \vec x_j$ alone, and only the *difference* of momentum-derivatives appears as the Moyal-conjugate.)

Substituting a single Fourier mode $V_2(\vec r) = V^{(2)}_{\vec q}\thinspace\cos(\vec k\cdot\vec r + \phi)$ with $\vec k = \vec k_{\vec q}$ and using

$$\bigl(\overleftarrow{\nabla}_{\negthinspace\vec x_i} - \overleftarrow{\nabla}_{\negthinspace\vec x_j}\bigr)\cos(\vec k\cdot\vec r_{ij} + \phi) = -\vec k\thinspace\sin(\vec k\cdot\vec r_{ij} + \phi)\cdot 2,$$

(the factor 2 from differentiating $\vec r_{ij}$ once on each side, with opposite signs) one finds

$$\partial_t W^{(N)} \supset -\frac{V^{(2)}_{\vec q}}{\hbar}\sin(\vec k\cdot \vec r_{ij}+\phi)\cdot \bigl[W^{(N)}(\dots,\vec p_i + \tfrac{\hbar\vec k}{2},\dots,\vec p_j - \tfrac{\hbar\vec k}{2},\dots) - W^{(N)}(\dots,\vec p_i - \tfrac{\hbar\vec k}{2},\dots,\vec p_j + \tfrac{\hbar\vec k}{2},\dots)\bigr],$$

which is exactly a finite-difference operator implementing the rule of §6.3 (with the source-to-destination convention absorbing the overall factor of 2). Total momentum $\vec p_i + \vec p_j$ is unaffected by the shifts, so it is conserved exactly per event.

The bound $\vert \Gamma^{(2)}_{\vec q}\vert \le \vert V^{(2)}_{\vec q}\vert /\hbar$ — independent of $\vec r_{ij}$ — gives a uniform per-mode rate cap, useful for $\Delta t$ selection.

---

## 8. Conservation and bounding properties

### 8.1 Conserved exactly per event

- **World count $\mathcal{W}$.** No spawning, no annihilation. (Identical to the 1+1D spec in this regard.)
- **Per-world particle count $N$.** Trivially — particles in a world are persistent labels.
- **Total momentum $\sum_i \vec p_i^{(\alpha)}$ per world**, *for any pair-wise interaction with translation invariance*. Each two-body event delivers $+\hbar\vec k$ to one particle and $-\hbar\vec k$ to the other.

### 8.2 Conserved on average / by construction

- **Energy.** First-order in $\Delta t$, conserved up to splitting error; second-order Strang splitting reduces this further. The discrete-jump rule does not exactly conserve energy per event (a single Fourier-mode kick of $\hbar\vec k$ shifts kinetic energy), reflecting the fact that energy conservation in the QLE is a property of the full Moyal series, not of any individual term.

### 8.3 Bounding

The bound $|W^{(N)}| \le (2/h)^{dN}$ bounds $W'^{(N)} \le 2(2/h)^{dN}$. For the world ensemble, this translates to a bound on local sample density: in any phase-space volume $V_{\text{loc}}$ of the joint $2dN$-dimensional space, the expected world count is $\le 2\nu_{\negthinspace\mathcal{W}}(2/h)^{dN}\cdot V_{\text{loc}}$. The crystal-lattice algorithm's freedom from explicit pair annihilation persists in the multi-body case.

---

## 9. Identical particles

For $N$ identical bosons or fermions the joint Wigner function carries permutation symmetry. Two implementation regimes:

### 9.1 Distinguishable approximation (recommended starting point) **[choice]**

Treat particles as distinguishable for the dynamics and post-symmetrize observables. For one- and two-body marginals at moderate $N$, this is exact in the absence of exchange-driven correlations and approximate otherwise. Standard practice in Wigner-MC implementations.

### 9.2 Symmetrized world initialization

Sample initial worlds from a (anti)symmetrized joint Wigner distribution. The dynamics in §6 commutes with permutations (all rules are particle-symmetric: every pair contributes equally, the streaming rule is per-particle), so an initially symmetric ensemble remains symmetric *in expectation* — though stochastic fluctuations cause individual worlds to drift away from the symmetric subspace, requiring periodic resymmetrization or a fixed-node-style constraint **[choice]**.

### 9.3 Fermionic sign problem

The lower bound $W^{(N)} \ge -(2/h)^{dN}$ does not save one from the sign problem. For antisymmetric joint states, $|W^{(N)}|$ saturates near its bound throughout the dynamically relevant region of phase space, and observables emerge as small differences between large near-cancelling positon and (implicit) negaton populations. The crystal-lattice formulation buys cleaner accounting — the negaton background is a constant offset rather than a dynamical population — but does not change the underlying statistical signal-to-noise. Approximate techniques (fixed-node, constrained-path, phaseless-AFQMC-style) developed for other $N$-body Monte Carlo methods translate to this framework with the same scope and the same caveats.

---

## 10. Coulomb interactions in the $N$-body algorithm

Apply the Ewald split per pair. For each $(i, j)$:

- **Short-range piece** $V_2^{\text{short}}(r) = e^2 \mathrm{erfc}(\alpha r)/r$. Treat as a **deterministic momentum-conserving drift**:

$$\Delta \vec p_i^{(\alpha)} = -\Delta \vec p_j^{(\alpha)} = -\nabla_{\negthinspace\vec r}V_2^{\text{short}}(\vec r_{ij}^{(\alpha)})\thinspace\Delta t.$$

  Linear-scaling in $N$ via cell lists at the cutoff radius (standard MD machinery). The $O(\hbar^2)$ Moyal corrections to this drift are negligible if $\alpha$ is chosen so that $V_2^{\text{short}}$ varies slowly on the de Broglie scale.

- **Long-range piece** $V_2^{\text{long}}(r) = e^2 \mathrm{erf}(\alpha r)/r$. Apply the §6.3 stochastic-jump rule with Fourier coefficients

$$\tilde V_{\vec k}^{\text{long}} \propto \frac{e^{-|\vec k|^2/4\alpha^2}}{|\vec k|^2}.$$

  Exponential decay in $|\vec k|$ means a small, finite mode set suffices for any specified accuracy.

The structure parallels Particle-Particle Particle-Mesh (P³M) and Smooth Particle-Mesh Ewald (SPME) molecular dynamics, with the long-range part contributing momentum-conserving discrete *quantum* jumps in place of smooth classical forces. Total momentum is conserved by construction in both pieces.

---

## 11. Reference pseudocode ($N$-body, particle-Monte-Carlo form)

```python
# Inputs:
#   N_part: number of particles, masses[i], charges[i] (for Coulomb)
#   V_ext_modes: dict q_tuple -> (V_q, phi_q) for external potential
#   V2_modes:    dict q_tuple -> (V_q, phi_q) for pair potential's long-range
#   short_range_force: callable r_ij_vec -> F_ij_vec for short-range pair force
#   W_workers: number of worlds
#   dt, num_steps
#   sample_initial_world(): returns (x[N,d], p[N,d]) sampled from W_init
import numpy as np

D = 3
worlds_x = np.empty((W_workers, N_part, D))
worlds_p = np.empty((W_workers, N_part, D))
for alpha in range(W_workers):
    worlds_x[alpha], worlds_p[alpha] = sample_initial_world()

# Pre-compute k-vectors and rate amplitudes for each Fourier mode
def precompute_modes(modes_dict):
    items = []
    for q_tuple, (V_q, phi_q) in modes_dict.items():
        k = np.array([2*np.pi*q_tuple[a]/L[a] for a in range(D)])
        items.append((k, V_q, phi_q))
    return items

ext_items = precompute_modes(V_ext_modes)
pair_items = precompute_modes(V2_modes)

for step in range(num_steps):

    # ---- 6.1 Free streaming, all worlds, all particles ----------------
    worlds_x += (worlds_p / masses[None, :, None]) * dt
    worlds_x = worlds_x % L_array[None, None, :]                    # periodic

    # ---- 6.2 External potential, per particle -------------------------
    for k, V_q, phi_q in ext_items:
        # phase per (world, particle)
        phase = np.einsum('a,wia->wi', k, worlds_x) + phi_q          # shape (W, N)
        Gamma = -(V_q / hbar) * np.sin(phase)                        # signed
        rate = np.abs(Gamma) * dt
        events = np.random.random(rate.shape) < rate                 # boolean
        sign = np.sign(Gamma)
        # delta_p has shape (W, N, D); only event cells get a kick
        delta_p = events[..., None] * sign[..., None] * (hbar * k)[None, None, :]
        worlds_p += delta_p

    # ---- 6.3 Two-body long-range jumps --------------------------------
    # Loop over unordered pairs; vectorize over worlds.
    for i in range(N_part):
        for j in range(i+1, N_part):
            r_ij = worlds_x[:, i, :] - worlds_x[:, j, :]             # (W, D)
            # minimum-image convention for periodic box
            r_ij = ((r_ij + L_array/2) % L_array) - L_array/2
            for k, V_q, phi_q in pair_items:
                phase = r_ij @ k + phi_q                             # (W,)
                Gamma = -(V_q / hbar) * np.sin(phase)
                rate = np.abs(Gamma) * dt
                events = np.random.random(rate.shape) < rate
                sign = np.sign(Gamma)
                kick = events[:, None] * sign[:, None] * (hbar * k)[None, :]
                worlds_p[:, i, :] += kick
                worlds_p[:, j, :] -= kick

    # ---- 10. Two-body short-range deterministic drift -----------------
    for i in range(N_part):
        for j in range(i+1, N_part):
            r_ij = worlds_x[:, i, :] - worlds_x[:, j, :]
            r_ij = ((r_ij + L_array/2) % L_array) - L_array/2
            F_ij = short_range_force(r_ij)                           # (W, D)
            worlds_p[:, i, :] += F_ij * dt
            worlds_p[:, j, :] -= F_ij * dt
```

For statistics, accumulate one- and two-body marginals from `worlds_x, worlds_p` directly by histogramming over the world axis.

---

## 12. Open implementation questions

1. **Statistical convergence.** $\mathcal{W}$ scales empirically as $\mathcal{O}(\sigma_{\text{obs}}^{-2})$ for a given observable variance target. One-body marginals on smooth states tend to converge with $\mathcal{W} \sim 10^4$–$10^5$; two-body correlations and exchange-sensitive observables can require an order of magnitude more. No exponential-in-$N$ scaling for low-body observables.

2. **Rate-budget regularization.** When $|V^{(2)}_{\vec q}|\thinspace\Delta t/\hbar$ approaches 1 for any single mode, the per-step Bernoulli sampling is no longer a good approximation to the true Poisson process; switch to an exponential-time-to-event sampler **[choice]** or reduce $\Delta t$.

3. **Mean-field reference (TDHF-Wigner).** Replace pair-wise correlated jumps by each particle seeing the time-dependent mean field of the others (the marginal $\rho^{(\text{mean})}(\vec x) = N^{-1}\sum_j \rho_j(\vec x)$ from the world ensemble). This recovers a TDHF / TDDFT-like dynamics at much lower statistical cost; useful as a reference and as a starting point for adding pair correlations perturbatively.

4. **Branching/resampling for variance reduction.** For long-time evolution, populations of worlds in different regions of $\Omega^{(N)}$ may become imbalanced. Periodic resampling (importance-reweight + branch high-weight worlds, kill low-weight ones) is standard in DMC and would translate here. Care needed to preserve the unbiased estimator.

5. **Symmetry preservation.** Periodic resymmetrization of the world ensemble for identical-particle systems — *how often* and *how* — is implementation-dependent and likely problem-dependent.

6. **Initialization from a realistic state.** Sampling from a known $W^{(N)}_0$ is straightforward only for product states (each particle drawn independently from its own $w_i$). For correlated initial states (Slater determinants, Bose-condensed product states with correlations), this is itself a sampling problem — possibly one solvable by short imaginary-time evolution from a tractable trial state.

7. **Mediator interpretation in $N$-body.** The single-particle "positons mediate jumps in other positons" picture does not carry over literally to the multi-body case: different worlds in the ensemble are independent samples and do not interact. The clean reinterpretation is **per-world pair-wise virtual-quantum exchange**: within each world, every pair exchanges momentum quanta whose direction and magnitude are set by the Fourier modes of $V_2$. The mediator role is now played by the *gauge boson* (photon for Coulomb), making the algorithm structurally closer to a stochastic implementation of QED than to single-particle Wigner dynamics.

8. **Half-jump vs full-jump conventions.** §6.2/§6.3 deliver full kicks of $\hbar\vec k$ in continuous-momentum world coordinates; the 1+1D mesh form delivered half-jumps of $\hbar\vec k/2$ paired symmetrically. Either convention is consistent with the underlying QLE; the half-jump form preserves the visual symmetry of the source-to-destination diagram in the V2 memo, and may be preferred for documentation continuity. The full-jump form is more natural in the world ensemble. **[choice]**.

---

## 13. Suggested validation sequence

1. **Single particle, 3+3D, harmonic well + soft-core Coulomb.** Validate against the analytically known hydrogen-like Wigner function at the ground state and against a Gaussian-wavepacket scattering trajectory. Forces the vectorial-jump infrastructure to be built and tested with no ensemble complications.

2. **Two distinguishable particles, 1+1D, soft-core Coulomb interaction.** Smallest non-trivial test of the correlated-jump rule. Validate against direct CI on a small joint $(x_1, x_2, p_1, p_2)$ mesh ($M_x = M_p = 32$ is feasible: $\sim 10^6$ cells).

3. **Two distinguishable particles, 3+3D.** Extends (2) with vectorial pair jumps and Ewald-split Coulomb. Validate against analytic helium-like hydrogenic-orbital-product trial states.

4. **$N$ identical fermions in a harmonic trap, distinguishable approximation.** First multi-particle test; tests world-ensemble statistics. Validate one-body marginal against TDHF/Hartree reference.

5. **$N$ identical fermions, symmetrized.** Brings in the sign problem in a controlled setting; validates resymmetrization scheme.

---

## 14. Sources and connections

- `docs/algorithm/phase_space_crystal_lattice_algorithm.md` — the 1+1D spec this extends.
- `docs/supplement/phase_space_crystal_lattice_supplement.md` — sign convention and derivation of the single-particle jump rule.
- `Extended_Fokker_Planck_Eq_and_the_QLE_V2.pdf` — primary analytical source; the BBGKY discussion (V2 §"The Extended Fokker-Planck equations and the QLE") frames where the world-ensemble representation fits in the hierarchy of $N$-body density-evolution equations: the joint $W^{(N)}$ corresponds to the *full* Liouville equation without molecular-chaos truncation.
- `Wigner_Collisions_Diagram_sozi.pdf` — slides 10–15: photon-momentum identification per Fourier mode; this generalizes to *vectorial* photon momentum in the multi-dimensional / multi-body extension.
- *Signed-particle Wigner Monte Carlo* literature (Nedjalkov, Querlioz, Dollfus, Sellier, et al.): the world ensemble is an annihilation-free reformulation of the configurations used there.
- *Particle-mesh Ewald* literature (Darden, Essmann, Pedersen, et al.): provides the short-range/long-range split applied here for Coulomb pair interactions.
- *Diffusion Monte Carlo and AFQMC* (Foulkes et al., Zhang & Krakauer): the world ensemble plays the role of walkers / configurations; the fermionic sign problem is shared and the same approximate-projection toolkit (fixed-node, constrained-path, phaseless) applies.
