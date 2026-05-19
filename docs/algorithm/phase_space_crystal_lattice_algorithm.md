# Phase-Space-Crystal-Lattice Stochastic Algorithm

**A complete specification for simulating Wigner distribution evolution via a single-rule mediated-jump process.**

---

## 0. Status of this specification

This document is a **synthesis** of `docs/analysis/phase_space_crystal_lattice_review.md` and the underlying source documents (`Extended_Fokker_Planck_Eq_and_the_QLE_V2.pdf`, `Wigner_Collisions_Diagram_sozi.pdf`). The source documents specify the underlying physics and the per-event rules but stop short of a fully written-out algorithm — V2 closes with *"With this solution in hand I am in a position to construct a software simulation…"* (future work).

A redrafted, restructured version of the V2 memo focused on the crystal-lattice particle model and its stochastic microdynamics is at `docs/supplement/phase_space_crystal_lattice_supplement.md`. The sign convention used in §3 below is the QLE-consistent one; this differs from the simplified-form line on page 18 of the V2 memo (see supplement §6.3 for the algebraic step where the sign was lost).

Items below marked **[choice]** are implementation decisions not directly fixed by the source documents; alternatives are noted.

---

## 1. Phase space and discretization

1+1D phase space with periodic position domain $x \in [0, L]$ and momentum domain $p \in [-P_{\max}, +P_{\max}]$:

- $M$ position bins: $x_m = m \Delta x,\ \Delta x = L/M,\ m = 0, \dots, M{-}1$
- $N$ momentum bins: $p_n = (n - N/2)\thinspace\Delta p,\ n = 0, \dots, N{-}1$

For potentials whose Fourier expansion has fundamental wavelength $L$, the natural momentum spacing is

$$\Delta p = \frac{\pi\hbar}{L}$$

**[choice]** — but justified by the source documents' identification of $n\pi\hbar/L$ as the "photon" momentum quantum associated with the $n$-th Fourier mode of $V$. With this choice, mode $q$ drives jumps of exactly $\pm q$ momentum cells. Finer momentum resolution is obtained with $\Delta p = \pi\hbar/(K L)$ for integer $K$; jumps then span $Kq$ cells.

---

## 2. State representation: the crystal lattice

The crystal-lattice shift produces the strictly non-negative quantity

$$W'(x, p)  \\; = \\;  W(x, p) + \frac{2}{h}, \qquad W' \in [0,\ 4/h]$$

This is represented as integer per-cell positon counts:

$$N_+(m, n)  \\; = \\;  \nu \cdot W'(x_m, p_n) \cdot \Delta x \thinspace \Delta p$$

where $\nu$ is the chosen particle-per-unit-phase-space-volume scale. The negaton background lattice is implicit: every cell would notionally hold the same count corresponding to $W_0 = 2/h$, and is **never updated dynamically**. Only $N_+$ evolves.

Bound: $0 \le N_+(m,n) \le 4\nu \Delta x \Delta p$, automatic from $|W| \le 2/h$.

---

## 3. Split-operator time evolution

Each timestep $\Delta t$ applies two phases in sequence. Lie splitting (apply each once) gives first-order accuracy in time; Strang splitting (half-step advection, full jump, half-step advection) gives second-order **[choice]**.

### 3a. Free streaming (advection)

QLE term: $-\dfrac{p}{m}\dfrac{\partial W}{\partial x}$.

On the lattice, shift row $n$ by

$$\Delta m_n = \mathrm{round}\negthinspace\left(\frac{p_n \Delta t}{m \Delta x}\right)$$

position cells (periodic). For zero advection error, pick $\Delta t$ so that the maximum shift is an exact integer:

$$\Delta t = \frac{m \Delta x}{|p_{\max}|} \cdot K, \quad K \in \mathbb{Z}_+$$

### 3b. Potential-driven mediated jumps

Decompose the potential into Fourier components:

$$V(x) = V_0 + \sum_{q \ge 1} V_q \cos\negthinspace\left(\frac{2\pi q x}{L} + \phi_q\right)$$

For each mode $q$, define

- **Local rate**: $\Gamma_q(x_m) = \dfrac{V_q}{\hbar}\cos\negthinspace\left(\dfrac{2\pi q x_m}{L} + \phi_q + \dfrac{\pi}{2}\right) = -\dfrac{V_q}{\hbar}\sin\negthinspace\left(\dfrac{2\pi q x_m}{L} + \phi_q\right)$
- **Half-jump in cells**: $\delta_q = q$ (with the natural $\Delta p = \pi\hbar/L$ choice)
- **Total source-to-destination jump**: $2\delta_q$ cells, equal to $2q\pi\hbar/L = qh/L$ (the photon momentum quantum at wavelength $L/q$)

**The single rule.** A positon at cell $(m, n)$ acts as a *mediator*: with probability $|\Gamma_q(x_m)|\thinspace\Delta t$ per particle per mode per timestep **[choice — Poisson rate is the rigorous form when this is not small]**, it induces a transfer of one particle:

$$\text{if } \Gamma_q(x_m) > 0:\quad (m,\thinspace n + q) \longrightarrow (m,\thinspace n - q)$$
$$\text{if } \Gamma_q(x_m) < 0:\quad (m,\thinspace n - q) \longrightarrow (m,\thinspace n + q)$$

The mediator itself is unchanged. No new particles are created. This is the entire crystal-lattice rule.

### 3c. Equivalent mesh-density form

When $\nu$ is large and we evolve $W$ directly rather than counting particles:

$$W(x_m, p_n, t + \Delta t) = W(x_m, p_n, t) + \Delta t \sum_{q \ge 1} \Gamma_q(x_m)\bigl[W(x_m, p_{n+q}) - W(x_m, p_{n-q})\bigr]$$

In the small-$`\Delta p`$ continuum limit this reproduces the QLE force term $\partial_t W = +\thinspace V'(x)\thinspace\partial_p W$. (See `docs/supplement/phase_space_crystal_lattice_supplement.md` §6 for the derivation; it corrects the sign that appears in the simplified form of the V2 memo.)

In Python (`p_axis` is the momentum axis, with cell index increasing with $p$): `np.roll(W, +1, axis=p_axis)` brings $W(p - \Delta p)$ to the row at $p$. The full update for the $q=1$ canonical case is

```python
W += (V_max / hbar) * dt * np.sin(2*np.pi*X/L) * (
        np.roll(W, +1, axis=p_axis) - np.roll(W, -1, axis=p_axis))
```

To verify the sign: $\Gamma_1 = -(V_{\max}/\hbar)\sin\theta$ and the bracket is $W(p+\Delta p) - W(p-\Delta p)$, giving $\Delta W = -(V_{\max}/\hbar)\sin\theta\thinspace[W_{\rm hi} - W_{\rm lo}]\thinspace dt$, equivalently $+(V_{\max}/\hbar)\sin\theta\thinspace[W_{\rm lo} - W_{\rm hi}]\thinspace dt$, which the Python expresses with the order `roll(+1) - roll(-1)`.

---

## 4. Population bounding — automatic, no explicit annihilation

The signed-particle Wigner Monte Carlo approach requires explicit positon/negaton annihilation when opposite-sign particles enter the same cell. **The crystal-lattice algorithm avoids this entirely**:

1. The mediated jump rule conserves total particle count exactly per event (one particle moves, mediator unchanged).
2. The shift $W' = W + 2/h$ keeps $N_+ \ge 0$ for all admissible Wigner states.
3. The bound $|W| \le 2/h$ guarantees $N_+ \le 4\nu\Delta x\Delta p$, so cells cannot overflow.

This is the central algorithmic advantage of the crystal-lattice formulation cited in the source documents.

---

## 5. Initialization and observable extraction

### Initialize

```python
W0 = compute_initial_wigner(x_grid, p_grid)   # may have negative regions
N_plus = np.round((W0 + 2/h) * nu * dx * dp).astype(int)
```

### Run

Apply 3a + 3b for the desired number of timesteps.

### Extract observables

```python
W = N_plus / (nu * dx * dp) - 2/h
rho_x = W.sum(axis=p_axis) * dp     # position marginal
rho_p = W.sum(axis=x_axis) * dx     # momentum marginal
```

Per the source documents, only the *excess* (un-shifted) population should contribute to laboratory observables — exactly the quantity recovered by subtracting the $2/h$ background.

---

## 6. Reference pseudocode

```python
# Inputs: V_q, phi_q for q = 1..Q ; W0(x,p) ; M, N, nu, dt, num_steps
import numpy as np

# Grid
dx = L / M
dp = np.pi * hbar / L
x  = np.arange(M) * dx
p  = (np.arange(N) - N//2) * dp

# Initialize
N_plus = np.round((W0(x[None,:], p[:,None]) + 2/h) * nu * dx * dp).astype(int)

# Pre-compute spatial rate profiles for each Fourier mode
Gamma = np.zeros((Q+1, M))
for q in range(1, Q+1):
    Gamma[q] = -(V_q[q] / hbar) * np.sin(2*np.pi*q*x/L + phi_q[q])

for step in range(num_steps):

    # ---- 3a. Free streaming -------------------------------------------
    for n in range(N):
        shift = int(round(p[n] * dt / (mass * dx)))
        if shift:
            N_plus[n, :] = np.roll(N_plus[n, :], shift)

    # ---- 3b. Mediated jumps -------------------------------------------
    # Particle-resolved version (true Monte Carlo).
    # Direction convention: when Gamma > 0, transfer is from (n+q) -> (n-q),
    # which corresponds to the QLE force term +V'(x) dW/dp (see supplement).
    for q in range(1, Q+1):
        for m in range(M):
            rate_dt = Gamma[q, m] * dt          # signed
            sign    = 1 if rate_dt > 0 else -1
            mag     = abs(rate_dt)
            if mag == 0:
                continue
            for n in range(N):
                pop = N_plus[n, m]
                if pop == 0:
                    continue
                # Number of mediator events in this cell this step
                events = np.random.binomial(pop, min(mag, 1.0))
                # Source / destination cells (periodic in p)
                src = (n + sign * q) % N
                dst = (n - sign * q) % N
                # Cannot transfer more than the source has
                events = min(events, N_plus[src, m])
                N_plus[src, m] -= events
                N_plus[dst, m] += events

# Recover Wigner distribution
W = N_plus / (nu * dx * dp) - 2/h
```

For large $\nu$, replace the inner Monte-Carlo loop with the deterministic mesh-density update of Section 3c — this is the "split-Fourier" path referenced in the source documents.

---

## 7. Generalization beyond cosine potentials

For a general bounded potential $V(x)$:

1. Compute the Fourier series of $V$ on $[0, L]$, retaining $Q$ modes.
2. The error from truncation is $O(\sum_{q > Q} |V_q|)$.
3. The leading higher-order Moyal term $\dfrac{\hbar^2}{24}\dfrac{\partial^3 V}{\partial x^3}\dfrac{\partial^3 W}{\partial p^3}$ is automatically included because it is exactly the third Taylor term of the cosine-component contributions — no separate handling is needed.

For **unbounded** polynomial potentials (e.g. $V = x^3$), the source documents note that the exact jump density is impulsive (sums of Dirac-delta derivatives), which is unsuitable for direct simulation. Two paths forward:

- **(a)** Soft-confine: replace $V$ by $V \cdot \chi(x)$ for a smooth cutoff $\chi$ over the simulation window.
- **(b)** Use the QLE differential form directly for the jump term, evaluated on the discretized grid.

---

## 8. Open implementation questions

1. **Time-step selection.** The CFL-like advection condition wants $\Delta t \le m\Delta x / |p_{\max}|$; the jump-rate condition wants $\max_q |\Gamma_q| \cdot \Delta t \ll 1$. Whichever binds tighter sets $\Delta t$.

2. **Particle vs. density representation.** For coarse phase-space grids (small $M N$) and large statistics, the mesh-density form (Section 3c) is faster and exact. For high-dimensional generalizations or very low occupancy in regions, the particle form (Section 6 pseudocode) scales better.

3. **Higher-order splitting.** Strang splitting requires only minor restructuring (half-streaming, full jump, half-streaming).

4. **Multi-Fourier-mode sampling.** Looping over $q$ is the simple approach. A composite-rate sampler that picks $q$ proportional to $|V_q|$ then commits a single event per particle per step is more efficient when $Q$ is large.

5. **Conservation checks.** Total $\sum_{m,n} N_+(m,n)$ should be exactly conserved by both 3a and 3b; logging this is a useful debugging invariant.

---

## 9. Sources

- `Extended_Fokker_Planck_Eq_and_the_QLE_V2.pdf` — primary analytical source. Sections "Mapping the xFP model to the Split Fourier Algorithm", "Highly Simplified Model Examination", "Jump density-rate for sinusoidal potential".
- `Wigner_Collisions_Diagram_sozi.pdf` — slides 10–15: Fourier-component spawning rule, photon-momentum identification, phase-space-crystal-lattice interpretation.
- `docs/supplement/phase_space_crystal_lattice_supplement.md` — redrafted version of the V2 memo with the corrected sign convention.
- `docs/analysis/phase_space_crystal_lattice_review.md` — companion review and equation reference.
