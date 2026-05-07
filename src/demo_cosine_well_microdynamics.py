"""
Demonstration: cosine-well Gaussian wavepacket via crystal-lattice MICRODYNAMICS.

Companion to ``demo_cat_state_microdynamics.py``.  Where that demo evolves a
free-particle Schrödinger cat state (V = 0, so positon microdynamics is pure
free streaming and *is* the full QLE), this demo does the analogous thing in
a single-period **cosine potential**

    V(x) = V_p cos(2 pi x / L + pi)   =   - V_p cos(2 pi x / L)

so that x = 0 is the unique potential minimum on the periodic domain
``[-L/2, L/2)``.  An initially Gaussian wavepacket at the well bottom — width
matched to the harmonic approximation of V near x = 0 — is evolved by two
methods on the same grid:

  1. **Full QLE evolution** via the Phase-Space Crystal-Lattice mesh form
     (``PhaseSpaceCrystalLattice.strang_step_fourier``).  For a single-mode
     cosine potential the discrete-jump form is *exact*: it captures every
     order of the Moyal series, including the higher derivatives that
     produce negative regions in the Wigner function.

  2. **Classical-positon Monte Carlo**: sample N positons from
     ``W' = W + 2/h`` at t=0, evolve each one under Hamilton's equations
     (``dx/dt = p/m, dp/dt = -V'(x)``), bin, and reconstruct
     ``W = rho_emp - 2/h``.  This is exactly what the QLE would give in the
     classical limit hbar -> 0 (the leading Liouville term).  It does *not*
     include the discrete momentum-jump structure.

The pointwise difference (QLE - classical) at each snapshot time is then
exactly the quantum-interference content of the evolution.  In particular,
the negative Wigner regions appear in the QLE row but are absent (up to
binning noise) from the classical row.

Sign convention
---------------
The Phase-Space Crystal-Lattice library uses the QLE-consistent sign
convention ``Gamma_q(x) = -(V_q/hbar) * sin(2 pi q x / L + phi_q)``.  This is
the corrected sign described in
``docs/supplement/phase_space_crystal_lattice_supplement.md`` §6.3 and
implemented in ``src/wpmwlib/phase_space_crystal_lattice.py``.  For our well
configuration, ``V_q = V_p, phi_q = pi, q = 1``, giving
``Gamma_1(x) = +(V_p/hbar) sin(2 pi x / L)``: zero at the well bottom (no
jumps when classically there's no force) and at the maxima, with the
correct restoring direction in between.

Parameters
----------
The default parameters (``V_p = 1.5, L = 8, hbar = m = 1``) place the
ground-state Gaussian width at ``sigma_x ~= 0.72`` in a well of half-depth
``V_p = 1.5``.  Anharmonicity is mild: ``(k sigma_x)^2 / 12 ~= 0.02``.  To
make the quantum interference structure visible on a reasonable timescale,
the initial Gaussian is **mildly squeezed** (``sigma_init = 1.3 * sigma_x``),
which gives a state with ~0.62 hbar*omega energy that breathes at frequency
2 omega in the harmonic approximation; anharmonicity then mixes those
breathing oscillations and produces non-Gaussian Wigner structure within a
few classical periods.  The packet remains well confined to the well bottom
throughout the run.

Outputs
-------
- ``cosine_well_microdynamics_evolution.png``    3 x 4 panels (QLE / classical / diff)
- ``cosine_well_microdynamics_marginals.png``    rho(x, t) at four times
- ``cosine_well_microdynamics_trajectories.png`` classical phase-space portrait
- ``cosine_well_microdynamics_negativity.png``   Wigner negativity vs time

Run with:

    WPMW_OUTPUT=/mnt/user-data/outputs python -u demo_cosine_well_microdynamics.py
"""

from __future__ import annotations

import os
import sys
import time

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))

from wpmwlib.phase_space_crystal_lattice import (
    PhaseSpaceCrystalLattice,
    FourierMode,
)
from wpmwlib.wpmw_utils import output_path, docs_path


# --------------------------------------------------------------------- #
# Physics                                                               #
# --------------------------------------------------------------------- #
HBAR = 1.0
MASS = 1.0

L    = 8.0           # period of cosine; periodic cell is [-L/2, L/2)
V_P  = 1.5           # half-depth: V_min = -V_p, V_max = +V_p, well depth = 2V_p
PHI  = np.pi         # so V(x) = V_p cos(2 pi x/L + pi) places minimum at x=0

# Harmonic approximation about x = 0:
#   V(x) ~  -V_p + (1/2) m omega^2 x^2,  omega = (2 pi / L) sqrt(V_p / m)
OMEGA   = (2.0 * np.pi / L) * np.sqrt(V_P / MASS)
SIGMA_X_GS = np.sqrt(HBAR / (2.0 * MASS * OMEGA))   # harmonic ground-state width
SIGMA_P_GS = HBAR / (2.0 * SIGMA_X_GS)
T_PERIOD   = 2.0 * np.pi / OMEGA

# Initial Gaussian:  squeezed about (0, 0) so that anharmonic mixing
# produces clearly visible non-Gaussian Wigner structure within a few
# periods.  SQUEEZE = 2.0 gives a min-uncertainty squeezed-vacuum state
# with energy E = (hbar*omega/4)(s^2 + 1/s^2) ~= 1.06 hbar*omega -- only
# moderately above the harmonic ground state, but the position width is
# now large enough (~0.13 L) that the cosine well's quartic anharmonicity
# rapidly seeds Wigner negativity.
SQUEEZE   = 2.0
SIGMA_X_0 = SQUEEZE * SIGMA_X_GS
SIGMA_P_0 = HBAR / (2.0 * SIGMA_X_0)

T_FINAL = 4.0 * T_PERIOD
TIMES   = [0.0, T_PERIOD, 2.0 * T_PERIOD, 4.0 * T_PERIOD]


# --------------------------------------------------------------------- #
# Grid                                                                  #
# --------------------------------------------------------------------- #
M_GRID = 128
N_GRID = 128
DX     = L / M_GRID
DP     = np.pi * HBAR / L           # spec choice
P_MAX  = (N_GRID // 2) * DP


# --------------------------------------------------------------------- #
# Time-step selection                                                   #
# --------------------------------------------------------------------- #
# Two constraints:
#   (i)  CFL for spectral advection is loose, but we keep dt modest so the
#        Strang split with the explicit-Euler jump term is well-behaved.
#   (ii) Jump-rate condition: max |Gamma| * dt << 1.  Gamma_max = V_p / hbar.
GAMMA_MAX = V_P / HBAR
DT_JUMP   = 0.05 / GAMMA_MAX        # |Gamma|*dt <= 0.05 anywhere
DT        = min(DT_JUMP, T_PERIOD / 200.0)
N_STEPS   = int(round(T_FINAL / DT))


# --------------------------------------------------------------------- #
# MC ensemble parameters                                                #
# --------------------------------------------------------------------- #
M_SAMPLE = 256
N_SAMPLE = 256
DX_S     = L / M_SAMPLE
DP_S     = np.pi * HBAR / L
P_MAX_S  = (N_SAMPLE // 2) * DP_S

# Reconstruction grid (coarser than sampling, for shot-noise SNR per cell).
# At N_POSITONS = 5e6, 48 x 48 gives ~1100 background particles per cell ->
# reconstructed-W noise sigma ~ 0.010, well below QLE peak |W| ~ 0.318 and
# below expected QLE Wigner negative-region amplitude (~0.05).
M_RECON = 48
N_RECON = 48

N_POSITONS = 5_000_000
SEED       = 42


# --------------------------------------------------------------------- #
# Initial state and forces                                              #
# --------------------------------------------------------------------- #
def W_initial(X: np.ndarray, P: np.ndarray) -> np.ndarray:
    """Squeezed-Gaussian Wigner function at (0, 0).

    Min-uncertainty Gaussian with widths (SIGMA_X_0, SIGMA_P_0) such that
    SIGMA_X_0 * SIGMA_P_0 = hbar/2.
    """
    return (1.0 / (np.pi * HBAR)) * np.exp(
        -(X ** 2) / SIGMA_X_0 ** 2
        - (P ** 2) * SIGMA_X_0 ** 2 / HBAR ** 2
    )


def force(x: np.ndarray) -> np.ndarray:
    """Classical force F = -dV/dx for V(x) = -V_p cos(2 pi x / L)."""
    k = 2.0 * np.pi / L
    return -V_P * k * np.sin(k * x)


def hamilton_step(x: np.ndarray, p: np.ndarray, dt: float):
    """Velocity-Verlet (symplectic) step under the cosine well.  Periodic in x."""
    p = p + 0.5 * dt * force(x)
    x = x + dt * p / MASS
    p = p + 0.5 * dt * force(x)
    x = ((x + L / 2.0) % L) - L / 2.0
    return x, p


# In-place float32 step for the large MC ensemble.  ~4x faster than the
# allocating fp64 version because (i) memory bandwidth is the bottleneck
# at 5M particles and fp32 halves it, (ii) all ops are in-place so no
# 5M-element temporaries get allocated each step.
#
# Accuracy: fp32 mantissa is 23 bits.  Per-step relative error ~ 1e-7;
# cumulative random-walk-like drift over 800 steps ~ sqrt(800)*1e-7 ~ 3e-6
# in (x, p), far below the binning resolution dx_r ~ 0.17 and dp_r ~ 1.05.
def make_fp32_step(dt: float):
    """Closure-capture all constants as fp32 scalars to avoid up-casts."""
    K   = np.float32(2.0 * np.pi / L)
    FC  = np.float32(-V_P * K)
    LH  = np.float32(L / 2.0)
    LL  = np.float32(L)
    HDF = np.float32(0.5 * dt) * FC      # 0.5 dt F_coeff (combined)
    DM  = np.float32(dt / MASS)

    def step(x, p, scratch):
        # p += 0.5*dt*F(x)
        np.multiply(x, K, out=scratch)
        np.sin(scratch, out=scratch)
        scratch *= HDF
        p += scratch
        # x += dt p / m
        np.multiply(p, DM, out=scratch)
        x += scratch
        # p += 0.5*dt*F(x)
        np.multiply(x, K, out=scratch)
        np.sin(scratch, out=scratch)
        scratch *= HDF
        p += scratch
        # x = ((x + L/2) mod L) - L/2
        x += LH
        np.mod(x, LL, out=x)
        x -= LH

    return step


def hamilton_step_scalar(x: float, p: float, dt: float):
    """Same as hamilton_step but for plain Python scalars (used by orbit
    integration where calling np.sin on 1-element arrays would dominate
    the cost from numpy call overhead)."""
    import math
    k = 2.0 * math.pi / L
    F1 = -V_P * k * math.sin(k * x)
    p2 = p + 0.5 * dt * F1
    x2 = x + dt * p2 / MASS
    F2 = -V_P * k * math.sin(k * x2)
    p3 = p2 + 0.5 * dt * F2
    # wrap
    x2 = ((x2 + L / 2.0) % L) - L / 2.0
    return x2, p3


# --------------------------------------------------------------------- #
# Sampling positons from W'                                             #
# --------------------------------------------------------------------- #
def sample_positons(rng: np.random.Generator):
    """Sample N_POSITONS continuous-coordinate positons from W' = W + 2/h.

    Same inverse-CDF approach as in demo_cat_state_microdynamics.py.
    """
    x_s = (np.arange(M_SAMPLE) - M_SAMPLE // 2) * DX_S
    p_s = (np.arange(N_SAMPLE) - N_SAMPLE // 2) * DP_S
    Xs, Ps = np.meshgrid(x_s, p_s, indexing="xy")

    h = 2.0 * np.pi * HBAR
    Wp = W_initial(Xs, Ps) + 2.0 / h

    n_neg = int(np.sum(Wp < 0))
    if n_neg > 0:
        print(f"  clipping {n_neg} cells with W' < 0 to zero "
              f"(deepest = {float(np.min(Wp)):+.3e})")
    Wp = np.maximum(Wp, 0.0)

    total_mass = float(np.sum(Wp) * DX_S * DP_S)
    bg = (2.0 / h) * L * (2 * P_MAX_S)
    print(f"  ||W'||_1 = {total_mass:.4f}  "
          f"(physical norm 1 + background {bg:.4f})")

    flat = Wp.flatten() / Wp.sum()
    cdf = np.cumsum(flat)

    print(f"  drawing {N_POSITONS:,} samples via inverse CDF...")
    t0 = time.time()
    u = rng.random(N_POSITONS)
    idx = np.searchsorted(cdf, u)
    np.clip(idx, 0, M_SAMPLE * N_SAMPLE - 1, out=idx)
    print(f"    sampling wall = {time.time() - t0:.2f} s")

    j_idx = idx // M_SAMPLE
    i_idx = idx % M_SAMPLE
    x_part = x_s[i_idx] + (rng.random(N_POSITONS) - 0.5) * DX_S
    p_part = p_s[j_idx] + (rng.random(N_POSITONS) - 0.5) * DP_S

    x_part = ((x_part + L / 2) % L) - L / 2
    eps = 1e-9
    p_part = np.clip(p_part, -P_MAX_S + eps, P_MAX_S - eps)

    return x_part, p_part, total_mass


def reconstruct(x_part: np.ndarray, p_part: np.ndarray,
                total_mass: float):
    """Bin positons -> empirical W = rho_emp - 2/h on the recon grid."""
    x_edges = np.linspace(-L / 2.0, L / 2.0, M_RECON + 1)
    p_edges = np.linspace(-P_MAX, P_MAX, N_RECON + 1)
    counts, _, _ = np.histogram2d(x_part, p_part, bins=[x_edges, p_edges])
    counts = counts.T  # (x, p) -> (p, x)

    h = 2.0 * np.pi * HBAR
    dx_r = L / M_RECON
    dp_r = (2.0 * P_MAX) / N_RECON
    rho_emp = counts.astype(float) / N_POSITONS / (dx_r * dp_r) * total_mass
    W_emp = rho_emp - 2.0 / h

    x_c = 0.5 * (x_edges[:-1] + x_edges[1:])
    p_c = 0.5 * (p_edges[:-1] + p_edges[1:])
    return W_emp, x_c, p_c


# --------------------------------------------------------------------- #
# Solvers                                                               #
# --------------------------------------------------------------------- #
def run_psc_mesh():
    """Full QLE evolution via the PSC mesh form (spec §3c)."""
    psc = PhaseSpaceCrystalLattice(
        M=M_GRID, N=N_GRID, L=L, mass=MASS, hbar=HBAR,
        nu=None, advection="spectral",
    )
    psc.initialize_from_wigner(W_initial)
    norm0 = psc.total_norm()
    modes = [FourierMode(q=1, V_q=V_P, phi_q=PHI)]

    target_times = sorted([t for t in TIMES if t > 0.0])
    snapshots_psc = {0.0: psc.W.copy()}
    full_traj_t   = [0.0]
    full_traj_neg = [_negativity(psc.W, DX, DP)]

    print(f"  PSC mesh evolution: dt = {DT:.5f}, N_steps = {N_STEPS}")
    t0 = time.time()
    next_idx = 0
    t = 0.0
    for step in range(1, N_STEPS + 1):
        psc.strang_step_fourier(modes, DT)
        t = step * DT
        full_traj_t.append(t)
        full_traj_neg.append(_negativity(psc.W, DX, DP))
        while (next_idx < len(target_times)
               and t >= target_times[next_idx] - DT * 0.5):
            snapshots_psc[target_times[next_idx]] = psc.W.copy()
            next_idx += 1
    wall = time.time() - t0
    print(f"  PSC wall = {wall:.2f} s")
    print(f"  norm: initial = {norm0:.10f}, final = {psc.total_norm():.10f}")

    return {
        "x": psc.x, "p": psc.p,
        "snapshots": snapshots_psc,
        "negativity_t": np.array(full_traj_t),
        "negativity": np.array(full_traj_neg),
    }


def run_classical_mc(x_part: np.ndarray, p_part: np.ndarray,
                     total_mass: float):
    """Evolve all positons under Hamilton's equations; snapshot binned W."""
    target_times = sorted([t for t in TIMES if t > 0.0])
    snap0 = reconstruct(x_part, p_part, total_mass)
    snapshots_mc = {0.0: snap0}
    neg_t   = [0.0]
    neg_val = [_negativity(snap0[0], snap0[1][1] - snap0[1][0],
                                       snap0[2][1] - snap0[2][0])]

    print(f"  classical MC: {N_POSITONS:,} positons x {N_STEPS} steps  (fp32 in-place)")
    t0 = time.time()
    # Cast once to fp32; histogram2d is fine with fp32 input.
    x_p = x_part.astype(np.float32)
    p_p = p_part.astype(np.float32)
    scratch = np.empty(N_POSITONS, dtype=np.float32)
    fp32_step = make_fp32_step(DT)

    next_idx = 0
    t = 0.0
    # Do a coarse negativity log every K steps to keep cost down
    LOG_EVERY = max(1, N_STEPS // 60)

    for step in range(1, N_STEPS + 1):
        fp32_step(x_p, p_p, scratch)
        t = step * DT
        if step % LOG_EVERY == 0:
            snap = reconstruct(x_p, p_p, total_mass)
            neg_t.append(t)
            neg_val.append(_negativity(
                snap[0], snap[1][1] - snap[1][0], snap[2][1] - snap[2][0]
            ))
        while (next_idx < len(target_times)
               and t >= target_times[next_idx] - DT * 0.5):
            snapshots_mc[target_times[next_idx]] = reconstruct(
                x_p, p_p, total_mass
            )
            next_idx += 1
    wall = time.time() - t0
    print(f"  classical MC wall = {wall:.2f} s")

    return {
        "snapshots": snapshots_mc,
        "negativity_t": np.array(neg_t),
        "negativity": np.array(neg_val),
    }


def _negativity(W: np.ndarray, dx: float, dp: float) -> float:
    """Wigner negativity: int |min(W, 0)| dx dp."""
    return float(-np.sum(np.minimum(W, 0.0)) * dx * dp)


# --------------------------------------------------------------------- #
# Test trajectories (classical Hamilton orbits)                         #
# --------------------------------------------------------------------- #
def select_test_orbits():
    """Six representative classical orbits in the cosine well.

    Three "near-bottom" orbits inside ~1.5 sigma of the well minimum
    (small classical action; nearly closed elliptical orbits in the
    harmonic approximation).

    Three "wider" orbits at amplitude ~2-2.5 sigma where anharmonic
    deformation of the orbit becomes visible — but still well inside
    the well (turning points well below V_max = +V_p).
    """
    near = [
        (0.0,             0.0),                  # stationary at minimum
        (0.7 * SIGMA_X_GS, 0.0),
        (0.0,              0.7 * SIGMA_P_GS),
    ]
    wider = [
        (1.6 * SIGMA_X_GS, 0.0),
        (-1.6 * SIGMA_X_GS, 1.0 * SIGMA_P_GS),
        (0.0,              1.8 * SIGMA_P_GS),
    ]
    return near, wider


def integrate_orbit(seed, t_arr, dt_int=None):
    """Velocity-Verlet integration of one classical orbit at the requested times."""
    if dt_int is None:
        dt_int = T_PERIOD / 2000.0
    x = float(seed[0])
    p = float(seed[1])
    xs = [x]
    ps = [p]
    for k in range(1, len(t_arr)):
        n_steps = max(1, int(round((t_arr[k] - t_arr[k - 1]) / dt_int)))
        sub_dt = (t_arr[k] - t_arr[k - 1]) / n_steps
        for _ in range(n_steps):
            x, p = hamilton_step_scalar(x, p, sub_dt)
        xs.append(x)
        ps.append(p)
    return np.array(xs), np.array(ps)


def trajectory_with_wraps(seed, t_arr):
    """Classical orbit with NaN inserted at periodic-x wraps for clean plots."""
    xs, ps = integrate_orbit(seed, t_arr)
    out_x, out_p = [], []
    prev = None
    for i, xx in enumerate(xs):
        if prev is not None and abs(xx - prev) > L / 2:
            out_x.append(np.nan)
            out_p.append(np.nan)
        out_x.append(xx)
        out_p.append(ps[i])
        prev = xx
    return np.array(out_x), np.array(out_p)


# --------------------------------------------------------------------- #
# Plotting                                                              #
# --------------------------------------------------------------------- #
def save_fig(fig, name: str, dpi: int = 130) -> None:
    sp = output_path(name)
    fig.savefig(sp, dpi=dpi, bbox_inches="tight")
    print(f"Figure saved -> {sp}")
    dp = docs_path(name)
    if dp:
        fig.savefig(dp, dpi=dpi, bbox_inches="tight")
        print(f"Figure saved -> {dp}")


def plot_evolution(res_psc, res_mc, near, wider, name):
    times = sorted(res_psc["snapshots"].keys())
    near_colors  = plt.cm.Reds(np.linspace(0.5, 0.85, len(near)))
    wider_colors = plt.cm.Blues(np.linspace(0.5, 0.85, len(wider)))

    # Use t=0 PSC max for symmetric colour scale
    W0 = res_psc["snapshots"][times[0]]
    vmax = float(np.max(np.abs(W0)))
    levels = np.linspace(-vmax, vmax, 21)

    x_psc = res_psc["x"]
    p_psc = res_psc["p"]

    fig, axes = plt.subplots(3, 4, figsize=(15, 9), sharex=True, sharey=True)

    for j, t in enumerate(times):
        W_psc = res_psc["snapshots"][t]
        W_mc, x_mc, p_mc = res_mc["snapshots"][t]
        # Resample MC to PSC grid coordinates is not strictly needed for the
        # difference panel because we just compute it on the MC grid; instead
        # we'll show the MC and PSC each on their own grid and the difference
        # on the MC grid with the PSC interpolated linearly.
        Wpsc_on_mc = _bilinear_resample(W_psc, x_psc, p_psc, x_mc, p_mc)
        diff = W_mc - Wpsc_on_mc
        diff_max = float(np.max(np.abs(diff))) if np.any(diff) else 1e-14

        # row 0: PSC
        ax = axes[0, j]
        cf = ax.contourf(x_psc, p_psc, W_psc, levels=levels,
                         cmap="RdBu_r", extend="both")
        ax.set_title(f"PSC mesh (full QLE),  t = {t:.2f}", fontsize=10)
        ax.set_xlim(-L / 2, L / 2)
        ax.set_ylim(-3.0, 3.0)
        if j == 0:
            ax.set_ylabel("p")

        # row 1: classical MC
        ax = axes[1, j]
        ax.contourf(x_mc, p_mc, W_mc, levels=levels,
                    cmap="RdBu_r", extend="both")
        ax.set_title(f"classical-positon MC,  t = {t:.2f}", fontsize=10)
        if j == 0:
            ax.set_ylabel("p")

        # row 2: MC - PSC  (the quantum-interference content classical misses)
        ax = axes[2, j]
        ax.contourf(x_mc, p_mc, diff,
                    levels=np.linspace(-diff_max, diff_max, 21),
                    cmap="PuOr", extend="both")
        ax.set_title(f"MC − QLE,  $\\max = {diff_max:.2e}$", fontsize=10)
        ax.set_xlabel("x")
        if j == 0:
            ax.set_ylabel("p")

        # Overlay classical trajectories on rows 0 and 1
        t_arr = np.linspace(0.0, t, 60) if t > 0 else np.array([0.0])
        for seeds, colors in ((near, near_colors), (wider, wider_colors)):
            for s, c in zip(seeds, colors):
                if t > 0:
                    xv, pv = trajectory_with_wraps(s, t_arr)
                    xt_now, pt_now = xv[-1], pv[-1]
                else:
                    xv, pv = np.array([s[0]]), np.array([s[1]])
                    xt_now, pt_now = s[0], s[1]
                for r in (0, 1):
                    axes[r, j].plot(xv, pv, "-", color=c, lw=0.9, alpha=0.7)
                    axes[r, j].plot(xt_now, pt_now, "o",
                                    color=c, ms=5, mfc="none", mew=1.4)

    leg = [
        plt.Line2D([], [], color=plt.cm.Reds(0.7), marker="o", mfc="none",
                   mew=1.4, ms=6, lw=0.9,
                   label="near-bottom orbits"),
        plt.Line2D([], [], color=plt.cm.Blues(0.7), marker="o", mfc="none",
                   mew=1.4, ms=6, lw=0.9,
                   label="wider (anharmonic) orbits"),
    ]
    axes[1, 0].legend(handles=leg, fontsize=8, loc="upper right",
                      framealpha=0.9)

    fig.suptitle(
        "Cosine-well wavepacket via crystal-lattice MICRODYNAMICS\n"
        f"$V(x) = -V_p\\cos(2\\pi x/L)$, $V_p = {V_P}$, $L = {L}$, "
        f"$\\sigma_{{x,0}} = {SIGMA_X_0:.3f}$ "
        f"({SQUEEZE:.1f}×ground), $T_{{harm}} = {T_PERIOD:.2f}$.\n"
        f"{N_POSITONS:,} positons sampled from $W' = W + 2/h$ at $t=0$, "
        "evolved by classical Hamilton dynamics (row 2). "
        "MC − QLE (row 3) is the quantum content classical misses.",
        fontsize=10, y=0.99,
    )
    fig.subplots_adjust(top=0.88, bottom=0.07, left=0.05, right=0.93,
                        wspace=0.15, hspace=0.30)
    cax = fig.add_axes([0.95, 0.40, 0.012, 0.45])
    cbar = fig.colorbar(cf, cax=cax)
    cbar.set_label("W")
    save_fig(fig, name)
    plt.close(fig)


def _bilinear_resample(W_src, x_src, p_src, x_dst, p_dst):
    """Cheap bilinear resample of W_src(x_src, p_src) onto (x_dst, p_dst)."""
    dx_src = x_src[1] - x_src[0]
    dp_src = p_src[1] - p_src[0]
    Nsp, Nsx = W_src.shape
    out = np.zeros((len(p_dst), len(x_dst)))
    for jp, p in enumerate(p_dst):
        fp = (p - p_src[0]) / dp_src
        jp0 = int(np.floor(fp))
        jp1 = jp0 + 1
        ap = fp - jp0
        jp0c = max(0, min(Nsp - 1, jp0))
        jp1c = max(0, min(Nsp - 1, jp1))
        for ix, x in enumerate(x_dst):
            fx = (x - x_src[0]) / dx_src
            ix0 = int(np.floor(fx))
            ix1 = ix0 + 1
            ax = fx - ix0
            ix0c = max(0, min(Nsx - 1, ix0))
            ix1c = max(0, min(Nsx - 1, ix1))
            out[jp, ix] = (
                (1 - ap) * (1 - ax) * W_src[jp0c, ix0c]
                + (1 - ap) *      ax * W_src[jp0c, ix1c]
                +      ap * (1 - ax) * W_src[jp1c, ix0c]
                +      ap *      ax * W_src[jp1c, ix1c]
            )
    return out


def plot_marginals(res_psc, res_mc, name):
    times = sorted(res_psc["snapshots"].keys())
    fig, axes = plt.subplots(1, len(times), figsize=(15, 3.2), sharey=True)
    x_psc = res_psc["x"]
    for j, t in enumerate(times):
        W_psc = res_psc["snapshots"][t]
        W_mc, x_mc, p_mc = res_mc["snapshots"][t]
        rho_psc = np.sum(W_psc, axis=0) * (res_psc["p"][1] - res_psc["p"][0])
        rho_mc  = np.sum(W_mc,  axis=0) * (p_mc[1] - p_mc[0])
        ax = axes[j]
        ax.plot(x_psc, rho_psc, "k-",  lw=1.5, label="QLE (mesh)")
        ax.plot(x_mc,  rho_mc,  "C3",  lw=1.0, alpha=0.8, label="classical MC")
        ax.set_title(f"t = {t:.2f}", fontsize=11)
        ax.set_xlabel("x")
        ax.set_xlim(-L / 2, L / 2)
        ax.grid(alpha=0.3)
        if j == 0:
            ax.set_ylabel(r"$\rho(x)$")
            ax.legend(fontsize=9, loc="upper right")
    fig.suptitle(
        "Position-space probability density: full QLE vs classical-positon MC\n"
        "(integrating W over p collapses the negative regions, "
        "so the two should agree well in the marginal even where W differs.)",
        fontsize=10, y=1.06,
    )
    fig.tight_layout()
    save_fig(fig, name)
    plt.close(fig)


def plot_trajectory_portrait(near, wider, name):
    near_colors  = plt.cm.Reds(np.linspace(0.5, 0.85, len(near)))
    wider_colors = plt.cm.Blues(np.linspace(0.5, 0.85, len(wider)))

    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    ax.set_facecolor("#f8f8f8")

    # Background: V(x) contour line for context (right axis)
    x_bg = np.linspace(-L / 2, L / 2, 400)
    V_bg = -V_P * np.cos(2 * np.pi * x_bg / L)
    ax2 = ax.twinx()
    ax2.plot(x_bg, V_bg, color="gray", lw=1.0, alpha=0.6, ls="--")
    ax2.set_ylabel("V(x)", color="gray")
    ax2.tick_params(axis="y", colors="gray")
    ax2.set_ylim(-V_P * 1.1, V_P * 1.1)

    t_arr = np.linspace(0.0, T_FINAL, 600)
    for seeds, colors, label in (
        (near,  near_colors,  "near-bottom orbit"),
        (wider, wider_colors, "wider orbit"),
    ):
        first = True
        for s, c in zip(seeds, colors):
            xv, pv = trajectory_with_wraps(s, t_arr)
            ax.plot(xv, pv, "-", color=c, lw=1.2, alpha=0.85,
                    label=label if first else None)
            ax.plot(s[0], s[1], "o", color=c, ms=8, mfc="white", mew=2.0)
            first = False
            for t_mark in TIMES[1:]:
                xt, pt = integrate_orbit(s, np.array([0.0, t_mark]))
                xw = ((xt[-1] + L / 2) % L) - L / 2
                ax.plot(xw, pt[-1], "|", color=c, ms=8, mew=1.5)

    ax.set_xlim(-L / 2, L / 2)
    ax.set_ylim(-2.5, 2.5)
    ax.set_xlabel("x"); ax.set_ylabel("p")
    ax.axhline(0, color="gray", lw=0.5)
    ax.axvline(0, color="gray", lw=0.5)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=10, loc="upper right", framealpha=0.95)
    ax.set_title(
        "Classical Hamilton orbits in $V(x) = -V_p\\cos(2\\pi x/L)$\n"
        "Open circle: $t = 0$.  Tick marks: $t = T, 2T, 3T$.\n"
        "Dashed grey curve: V(x) on the right axis.\n"
        "Inner orbits are nearly elliptical (harmonic regime); outer orbits\n"
        "show visible deformation from the cosine anharmonicity.",
        fontsize=10,
    )
    fig.tight_layout()
    save_fig(fig, name)
    plt.close(fig)


def plot_negativity(res_psc, res_mc, name):
    fig, ax = plt.subplots(1, 1, figsize=(8, 4.5))
    ax.plot(res_psc["negativity_t"], res_psc["negativity"], "k-",
            lw=1.5, label="QLE (mesh)")
    ax.plot(res_mc["negativity_t"], res_mc["negativity"], "C3-",
            lw=1.0, alpha=0.8, label="classical MC")
    for t in TIMES[1:]:
        ax.axvline(t, color="gray", lw=0.5, ls=":")
    ax.set_xlabel("t")
    ax.set_ylabel(r"Wigner negativity   $\int |\min(W,0)|\,dx\,dp$")
    ax.set_title(
        "Wigner-function negativity over time\n"
        "QLE accumulates real negative regions; classical MC only shows\n"
        "shot-noise-floor negativity from finite-N binning."
    )
    ax.grid(alpha=0.3)
    ax.legend(fontsize=10, loc="best")
    fig.tight_layout()
    save_fig(fig, name)
    plt.close(fig)


# --------------------------------------------------------------------- #
# Main                                                                  #
# --------------------------------------------------------------------- #
def main():
    print("=" * 72)
    print("Cosine-well wavepacket — crystal-lattice MICRODYNAMICS")
    print("=" * 72)
    print(f"  hbar = {HBAR}, m = {MASS}")
    print(f"  L = {L}, V_p = {V_P}  (well depth = {2 * V_P}, V at ±L/2 = {V_P})")
    print(f"  harmonic approx:  omega = {OMEGA:.4f}, T = {T_PERIOD:.4f}")
    print(f"                    sigma_x_gs = {SIGMA_X_GS:.4f}, "
          f"sigma_p_gs = {SIGMA_P_GS:.4f}")
    print(f"  initial state:    sigma_x_0 = {SIGMA_X_0:.4f} "
          f"(= {SQUEEZE:.2f}× ground), centre (0, 0)")
    print(f"  grid: M = {M_GRID}, N = {N_GRID}, dx = {DX:.4f}, dp = {DP:.4f}")
    print(f"        recon grid: {M_RECON} × {N_RECON}")
    print(f"  time stepping: dt = {DT:.5f}, N_steps = {N_STEPS}, "
          f"T_final = {T_FINAL:.3f}")
    print(f"  jump-rate budget: max|Gamma| dt = {GAMMA_MAX * DT:.4f}")
    print()
    print(f"  k sigma_x_gs = {2 * np.pi / L * SIGMA_X_GS:.3f}  "
          f"(anharmonic-correction parameter, small => mild anharmonicity)")
    print()

    rng = np.random.default_rng(SEED)

    print("[step 1: full-QLE mesh evolution via PhaseSpaceCrystalLattice]")
    res_psc = run_psc_mesh()

    print()
    print("[step 2: sample positons from W' at t=0]")
    x_part, p_part, total_mass = sample_positons(rng)

    print()
    print("[step 3: classical-positon Hamilton evolution]")
    res_mc = run_classical_mc(x_part, p_part, total_mass)

    print()
    print("[step 4: select test orbits and plot]")
    near, wider = select_test_orbits()

    print()
    print("Snapshot summary:")
    print(f"  {'t':>7s}   {'QLE neg':>12s}   {'MC neg (noise)':>16s}   "
          f"{'max|MC-QLE|':>14s}")
    for t in sorted(res_psc["snapshots"].keys()):
        W_psc = res_psc["snapshots"][t]
        W_mc, x_mc, p_mc = res_mc["snapshots"][t]
        Wpsc_on_mc = _bilinear_resample(W_psc, res_psc["x"], res_psc["p"],
                                        x_mc, p_mc)
        dx_r = x_mc[1] - x_mc[0]; dp_r = p_mc[1] - p_mc[0]
        neg_psc = _negativity(W_psc, DX, DP)
        neg_mc  = _negativity(W_mc,  dx_r, dp_r)
        max_d   = float(np.max(np.abs(W_mc - Wpsc_on_mc)))
        print(f"  {t:7.4f}   {neg_psc:12.3e}   {neg_mc:16.3e}   {max_d:14.3e}")

    print()
    print("Test-orbit summary at t = T_final:")
    print(f"  {'category':>14s}   {'(x_0, p_0)':>16s}   "
          f"{'(x_final, p_final)':>22s}   {'E_class':>10s}")
    for cat, seeds in (("near-bottom", near), ("wider", wider)):
        for s in seeds:
            xt, pt = integrate_orbit(s, np.array([0.0, T_FINAL]))
            xw = ((xt[-1] + L / 2) % L) - L / 2
            E = pt[-1] ** 2 / (2 * MASS) - V_P * np.cos(2 * np.pi * xw / L)
            print(f"  {cat:>14s}   ({s[0]:+.3f}, {s[1]:+.3f})"
                  f"   ({xw:+.3f}, {pt[-1]:+.3f})    {E:+.4f}")

    plot_evolution(res_psc, res_mc, near, wider,
                   "cosine_well_microdynamics_evolution.png")
    plot_marginals(res_psc, res_mc,
                   "cosine_well_microdynamics_marginals.png")
    plot_trajectory_portrait(near, wider,
                             "cosine_well_microdynamics_trajectories.png")
    plot_negativity(res_psc, res_mc,
                    "cosine_well_microdynamics_negativity.png")
    print("\nDone.")


if __name__ == "__main__":
    main()
