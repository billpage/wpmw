"""
Demonstration: free-particle cat state via crystal-lattice MICRODYNAMICS.

Companion to ``demo_cat_state.py``.  Where that demo evolves the Wigner
function on a deterministic mesh, this one represents the same state as
an ensemble of *positons* sampled from the shifted distribution

    W'(x, p) = W(x, p) + 2/h

(everywhere non-negative since |W| <= 2/h for any pure physical state),
evolves each one as a classical particle, and reconstructs W from the
empirical particle density.

Background — see ``docs/supplement/phase_space_crystal_lattice_supplement.md``
sections 4.4 and 7.  In the crystal-lattice picture, every physical
Wigner state is represented as a *single-sign* density of positons on top
of a uniform crystalline background of rate 2/h.  The shifted total
density W' is non-negative and admits a probabilistic interpretation;
positons sampled from it propagate forward in time under the same
classical equations of motion that drive their lattice neighbours.

For free evolution V = 0 the QLE reduces to the classical Liouville
equation and there are no mediated jumps — positons stream ballistically.
Therefore:

  * The empirical reconstruction
        W_emp(x, p, t) = rho_emp(x, p, t) - 2/h
    converges to the analytic free-particle solution
        W(x, p, t) = W_0(x - p t/m, p)
    in the large-N limit, INCLUDING the negative interference fringes,
    despite the fact that no individual positon ever carries negative
    sign.  The negative regions emerge as a *deficit below the uniform
    background of 2/h*.

  * Individual positon trajectories are simply horizontal lines in
    phase space (constant p, x linear in t).  A positon that happens to
    sit on a Wigner-function nodal line at t = 0 is not in any way
    special at the microdynamic level: it streams forward with the same
    rule as every other positon.  The wave-like interference structure
    is purely an emergent ENSEMBLE property.

This demo makes both points quantitatively:

  * Compares the binned MC reconstruction to the closed-form Wigner
    function at four times.

  * Tracks six tagged "test" positons — three with initial conditions
    on/near interference-node lines of W_0 and three in non-node regions
    — and overlays their trajectories.  These are *not* part of the bulk
    ensemble; they are illustrative classical particles whose initial
    coordinates we choose by hand.

Run with:

    WPMW_OUTPUT=/mnt/user-data/outputs python -u demo_cat_state_microdynamics.py
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

from wpmwlib.wpmw_utils import output_path, docs_path
from demo_cat_state import (
    HBAR, MASS, X0, P0, SIGMA, T_C, L,
    W_cat_initial, W_cat_exact, sample_node_seeds,
)


def save_fig(fig, name: str, dpi: int = 130) -> None:
    """Save figure to both ``output_path(name)`` (always) and
    ``docs_path(name)`` (when ``WPMW_DOCS`` is set, so the figure also
    lands in the ``output``-branch worktree for later commit)."""
    sp = output_path(name)
    fig.savefig(sp, dpi=dpi, bbox_inches="tight")
    print(f"Figure saved -> {sp}")
    dp = docs_path(name)
    if dp:
        fig.savefig(dp, dpi=dpi, bbox_inches="tight")
        print(f"Figure saved -> {dp}")


# --------------------------------------------------------------------- #
# Grids and ensemble parameters                                         #
# --------------------------------------------------------------------- #
# Sampling grid (where W' is evaluated and inverse-CDF-sampled). Same as
# the deterministic demo for parity.
M_SAMPLE = 256
N_SAMPLE = 256
DX_S    = L / M_SAMPLE
DP_S    = np.pi * HBAR / L
P_MAX_S = (N_SAMPLE // 2) * DP_S

# Reconstruction grid for the binned empirical W.  Chosen as a divisor
# of the sampling resolution so each recon bin contains an integer number
# of sampling cells (no geometric variance from misaligned bin edges).
# 128 x 128 with N_POSITONS = 1e7 gives ~600 background particles per
# recon bin, hence reconstructed-W noise ~ 0.013 against a peak |W| ~
# 0.32 (SNR ~ 25).
M_RECON = 128
N_RECON = 128

# Snapshot times (match demo_cat_state.py for direct visual comparison).
TIMES = [0.0, T_C / 2.0, T_C, 1.5 * T_C]

# MC ensemble size.  At N = 1e7 with M_RECON = N_RECON = 96, expected
# count per cell from the uniform background alone is about 1100, giving
# a reconstructed-W noise floor of ~0.01 against a peak |W| ~ 0.32.
N_POSITONS = 20_000_000

SEED = 42


# --------------------------------------------------------------------- #
# Sampling                                                              #
# --------------------------------------------------------------------- #
def init_xp_grids():
    x = (np.arange(M_SAMPLE) - M_SAMPLE // 2) * DX_S
    p = (np.arange(N_SAMPLE) - N_SAMPLE // 2) * DP_S
    X, P = np.meshgrid(x, p, indexing="xy")
    return x, p, X, P


def sample_positons(rng: np.random.Generator):
    """Sample N_POSITONS positon positions from W' = W_cat_initial + 2/h.

    Uses inverse-CDF (np.searchsorted) on the flattened discrete
    distribution.  Returns the continuous (x, p) coordinates of every
    positon at t = 0.
    """
    x, p, X, P = init_xp_grids()
    W0 = W_cat_initial(X, P)
    h = 2.0 * np.pi * HBAR
    Wprime = W0 + 2.0 / h

    n_neg = int(np.sum(Wprime < 0))
    if n_neg > 0:
        min_val = float(np.min(Wprime))
        print(f"  clipping {n_neg} cells with W' < 0 to zero "
              f"(deepest = {min_val:+.3e})")
    Wprime = np.maximum(Wprime, 0.0)

    total_mass = float(np.sum(Wprime) * DX_S * DP_S)
    bg_estimate = (2.0 / h) * L * (2 * P_MAX_S)
    print(f"  ||W'||_1 = {total_mass:.4f}  "
          f"(decomposes as physical norm 1 + background {bg_estimate:.4f})")

    flat = Wprime.flatten()
    flat = flat / flat.sum()
    cdf = np.cumsum(flat)

    print(f"  drawing {N_POSITONS:,} samples via inverse CDF...")
    t0 = time.time()
    u = rng.random(N_POSITONS)
    idx = np.searchsorted(cdf, u)
    np.clip(idx, 0, M_SAMPLE * N_SAMPLE - 1, out=idx)
    print(f"    sampling wall = {time.time() - t0:.2f} s")

    j_idx = idx // M_SAMPLE   # p-axis (axis 0 of W shape)
    i_idx = idx % M_SAMPLE    # x-axis (axis 1)

    # Continuous within-cell positions: uniform jitter over the sampling
    # cell.  Without jitter, all 16 (or 4) particles in each
    # sampling-cell-cluster sit at exactly the same (x, p), introducing
    # spurious sub-grid structure when the recon grid is finer than the
    # sampling grid.  With jitter the particle positions are continuous.
    x_part = x[i_idx] + (rng.random(N_POSITONS) - 0.5) * DX_S
    p_part = p[j_idx] + (rng.random(N_POSITONS) - 0.5) * DP_S

    # Wrap x periodically to canonical [-L/2, L/2)
    x_part = ((x_part + L / 2) % L) - L / 2
    # Clamp p just inside [-P_MAX_S, P_MAX_S) -- only relevant for the
    # outermost cells, where W' is essentially the background.
    eps = 1e-9
    p_part = np.clip(p_part, -P_MAX_S + eps, P_MAX_S - eps)

    return x_part, p_part, total_mass


# --------------------------------------------------------------------- #
# Reconstruction                                                        #
# --------------------------------------------------------------------- #
def reconstruct(x_part: np.ndarray, p_part: np.ndarray,
                total_mass: float, t: float):
    """Stream the ensemble to time t, bin onto the reconstruction grid,
    and return (W_emp, x_centers, p_centers)."""
    x_t = x_part + p_part * t / MASS
    # Periodic wrap in x
    x_t = ((x_t + L / 2) % L) - L / 2

    x_edges = np.linspace(-L / 2, L / 2, M_RECON + 1)
    p_edges = np.linspace(-P_MAX_S, P_MAX_S, N_RECON + 1)

    counts, _, _ = np.histogram2d(x_t, p_part, bins=[x_edges, p_edges])
    counts = counts.T  # (x, p) -> (p, x) to match (N, M) convention

    h = 2.0 * np.pi * HBAR
    dx_r = L / M_RECON
    dp_r = (2.0 * P_MAX_S) / N_RECON
    rho_emp = counts.astype(float) / N_POSITONS / (dx_r * dp_r) * total_mass
    W_emp = rho_emp - 2.0 / h

    x_c = 0.5 * (x_edges[:-1] + x_edges[1:])
    p_c = 0.5 * (p_edges[:-1] + p_edges[1:])
    return W_emp, x_c, p_c


# --------------------------------------------------------------------- #
# Test particles                                                        #
# --------------------------------------------------------------------- #
def select_test_particles():
    """Pick six initial conditions for trajectory tracking.

    Three "near-node" particles sit on/near the cosine-fringe nodal
    lines of the t = 0 cat-state Wigner function:

        p_0 x + p x_0 = (k + 1/2) pi hbar / 2  (k integer)

    Three "non-node" particles sit at characteristic NON-node locations:
    the central peak (between fringes), and the centres of the two
    Gaussian lobes.

    These are illustrative test particles, NOT drawn from the bulk MC
    ensemble.  They obey the same classical streaming dynamics as every
    other positon.
    """
    seeds_nn = sample_node_seeds()                  # five candidates
    near_node = [seeds_nn[0], seeds_nn[1], seeds_nn[4]]   # two stationary, one drifting
    non_node = [
        (0.0,   0.0),          # central anti-node (peak of cosine)
        (X0,   -P0),           # right lobe centre
        (-X0,  +P0),           # left lobe centre
    ]
    return near_node, non_node


def trajectory_with_wraps(seed, t_arr):
    """Compute (x_wrapped, p) along seed's free-particle characteristic,
    inserting NaN at periodic-boundary crossings for clean line plotting."""
    xs, ps = seed
    x_t = xs + ps * t_arr / MASS
    x_w = ((x_t + L / 2) % L) - L / 2
    out_x, out_p = [], []
    prev = None
    for xx in x_w:
        if prev is not None and abs(xx - prev) > L / 2:
            out_x.append(np.nan); out_p.append(np.nan)
        out_x.append(xx); out_p.append(ps)
        prev = xx
    return np.array(out_x), np.array(out_p)


# --------------------------------------------------------------------- #
# Plotting                                                              #
# --------------------------------------------------------------------- #
def plot_trajectory_portrait(near_node, non_node, name):
    """Phase-space plot showing ONLY the test-particle trajectories, no W
    underneath.  Drives home the point that at the microdynamic level
    every positon is a trivial free-streaming classical particle: the
    'near-node' tag refers solely to the location of the t=0 initial
    condition relative to the Wigner function's nodal structure, not to
    any feature of the trajectory itself.
    """
    nn_colors  = plt.cm.Reds(np.linspace(0.5, 0.85, len(near_node)))
    other_clrs = plt.cm.Blues(np.linspace(0.5, 0.85, len(non_node)))

    fig, ax = plt.subplots(1, 1, figsize=(8, 5))
    ax.set_facecolor("#f8f8f8")

    t_arr = np.linspace(0.0, max(TIMES), 200)
    for seeds, colors, label in (
        (near_node, nn_colors,  "near interference node"),
        (non_node,  other_clrs, "non-node region"),
    ):
        first = True
        for s, c in zip(seeds, colors):
            xv, pv = trajectory_with_wraps(s, t_arr)
            ax.plot(xv, pv, "-", color=c, lw=1.4, alpha=0.85,
                    label=label if first else None)
            ax.plot(s[0], s[1], "o", color=c, ms=8, mfc="white", mew=2.0)
            first = False
            # Mark snapshot times along each trajectory with small ticks
            for t_mark in TIMES[1:]:
                xt = s[0] + s[1] * t_mark / MASS
                xw = ((xt + L / 2) % L) - L / 2
                ax.plot(xw, s[1], "|", color=c, ms=8, mew=1.5)

    ax.set_xlim(-8, 8); ax.set_ylim(-3, 3)
    ax.set_xlabel("x"); ax.set_ylabel("p")
    ax.axhline(0, color="gray", lw=0.5)
    ax.axvline(0, color="gray", lw=0.5)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=10, loc="upper right", framealpha=0.95)
    ax.set_title(
        "Trajectories of the six test positons in phase space\n"
        "Open circle: $t = 0$.  Tick marks: $t = t_c/2,\\ t_c,\\ 3t_c/2$.\n"
        "Every trajectory is a horizontal line at constant $p$ — the "
        "'near-node' label refers only to where the particle started.",
        fontsize=10,
    )
    fig.tight_layout()
    save_fig(fig, name)
    plt.close(fig)


def plot_evolution(snapshots, near_node, non_node, name):
    times = sorted(snapshots.keys())
    nn_colors  = plt.cm.Reds(np.linspace(0.5, 0.85, len(near_node)))
    other_clrs = plt.cm.Blues(np.linspace(0.5, 0.85, len(non_node)))

    # Use t=0 exact W to set a consistent diverging colour scale.
    _, W_ex_0, _, _ = snapshots[times[0]]
    vmax = float(np.max(np.abs(W_ex_0)))
    levels = np.linspace(-vmax, vmax, 21)

    fig, axes = plt.subplots(3, 4, figsize=(15, 9), sharex=True, sharey=True)

    for j, t in enumerate(times):
        W_emp, W_ex, xs, ps = snapshots[t]
        diff = W_emp - W_ex
        diff_max = float(np.max(np.abs(diff)))

        ax = axes[0, j]
        cf = ax.contourf(xs, ps, W_emp, levels=levels,
                         cmap="RdBu_r", extend="both")
        ax.set_title(f"MC reconstruction,  t = {t:.2f}", fontsize=10)
        ax.set_xlim(-8, 8); ax.set_ylim(-3, 3)
        if j == 0:
            ax.set_ylabel("p")

        ax = axes[1, j]
        ax.contourf(xs, ps, W_ex, levels=levels,
                    cmap="RdBu_r", extend="both")
        ax.set_title(f"exact,  t = {t:.2f}", fontsize=10)
        if j == 0:
            ax.set_ylabel("p")

        ax = axes[2, j]
        ax.contourf(xs, ps, diff,
                    levels=np.linspace(-diff_max, diff_max, 21),
                    cmap="PuOr", extend="both")
        ax.set_title(f"MC − exact,  $\\max = {diff_max:.2e}$", fontsize=10)
        ax.set_xlabel("x")
        if j == 0:
            ax.set_ylabel("p")

        # Trajectory overlay on rows 0 and 1
        t_arr = np.linspace(0.0, t, 60)
        for seeds, colors in ((near_node, nn_colors),
                              (non_node, other_clrs)):
            for s, c in zip(seeds, colors):
                xv, pv = trajectory_with_wraps(s, t_arr)
                xt_now = s[0] + s[1] * t / MASS
                xw_now = ((xt_now + L / 2) % L) - L / 2
                for r in (0, 1):
                    axes[r, j].plot(xv, pv, "-", color=c, lw=0.9, alpha=0.7)
                    axes[r, j].plot(xw_now, s[1], "o", color=c,
                                    ms=5, mfc="none", mew=1.4)

    leg = [
        plt.Line2D([], [], color=plt.cm.Reds(0.7),  marker="o", mfc="none",
                   mew=1.4, ms=6, lw=0.9,
                   label="near interference node"),
        plt.Line2D([], [], color=plt.cm.Blues(0.7), marker="o", mfc="none",
                   mew=1.4, ms=6, lw=0.9,
                   label="non-node region"),
    ]
    axes[1, 0].legend(handles=leg, fontsize=8, loc="upper right",
                      framealpha=0.9)

    fig.suptitle(
        "Free-particle cat state via crystal-lattice MICRODYNAMICS\n"
        f"{N_POSITONS:,} positons sampled from $W' = W + 2/h$, streamed "
        "ballistically; $W$ recovered as $\\rho_{\\rm emp} - 2/h$.",
        fontsize=11, y=0.99,
    )
    fig.subplots_adjust(top=0.90, bottom=0.07, left=0.05, right=0.93,
                        wspace=0.15, hspace=0.30)
    cax = fig.add_axes([0.95, 0.40, 0.012, 0.45])
    cbar = fig.colorbar(cf, cax=cax)
    cbar.set_label("W")
    save_fig(fig, name)
    plt.close(fig)


def plot_marginals(snapshots, name):
    times = sorted(snapshots.keys())
    fig, axes = plt.subplots(1, len(times), figsize=(15, 3), sharey=True)
    for j, t in enumerate(times):
        W_emp, W_ex, xs, ps = snapshots[t]
        dp_r = ps[1] - ps[0]
        rho_mc = np.sum(W_emp, axis=0) * dp_r
        rho_ex = np.sum(W_ex,  axis=0) * dp_r
        ax = axes[j]
        ax.plot(xs, rho_ex, "k--", lw=1.5, label="exact $|\\psi|^2$")
        ax.plot(xs, rho_mc, "C3",  lw=1.0, label="MC binned")
        ax.set_title(f"t = {t:.2f}", fontsize=11)
        ax.set_xlabel("x")
        ax.set_xlim(-X0 - 4 * SIGMA, X0 + 4 * SIGMA)
        ax.grid(alpha=0.3)
        if j == 0:
            ax.set_ylabel("$\\rho(x)$")
            ax.legend(fontsize=9, loc="upper right")
    fig.suptitle(
        f"Position-space probability density: MC reconstruction vs exact "
        f"(t_c = {T_C:.2f})",
        fontsize=11, y=1.04,
    )
    fig.tight_layout()
    save_fig(fig, name)
    plt.close(fig)


# --------------------------------------------------------------------- #
# Main                                                                  #
# --------------------------------------------------------------------- #
def main():
    print("=" * 70)
    print("Free-particle cat state — crystal-lattice MICRODYNAMICS")
    print("=" * 70)
    print(f"  x0 = {X0}, p0 = {P0}, sigma = {SIGMA}, t_collision = {T_C}")
    print(f"  L = {L}, sampling grid {M_SAMPLE} x {N_SAMPLE}, "
          f"reconstruction grid {M_RECON} x {N_RECON}")
    print(f"  N_positons = {N_POSITONS:,}")
    print()
    print("[step 1: sampling positons from W' = W + 2/h]")
    rng = np.random.default_rng(SEED)
    x_part, p_part, total_mass = sample_positons(rng)

    near_node, non_node = select_test_particles()

    print()
    print("[step 2: streaming and binning at each snapshot time]")
    print(f"  {'t':>7s}   {'wall (s)':>10s}   "
          f"{'max|W_MC - W_ex|':>17s}   {'L2 deviation':>15s}")
    snapshots = {}
    for t in TIMES:
        t0 = time.time()
        W_emp, xs, ps = reconstruct(x_part, p_part, total_mass, t)
        wall = time.time() - t0
        Xr, Pr = np.meshgrid(xs, ps, indexing="xy")
        W_ex = W_cat_exact(Xr, Pr, t)
        diff = W_emp - W_ex
        l2 = float(np.sqrt(np.sum(diff ** 2) * (xs[1] - xs[0]) * (ps[1] - ps[0])))
        max_abs = float(np.max(np.abs(diff)))
        print(f"  {t:7.4f}   {wall:10.2f}   {max_abs:17.3e}   {l2:15.3e}")
        snapshots[t] = (W_emp, W_ex, xs, ps)

    print()
    print("[step 3: test-particle trajectories]")
    print("  Each particle moves at constant p along x = x_0 + p t.")
    print("  All test particles are 'just' free-streaming positons; the")
    print("  near-node tag refers only to their initial coordinates relative")
    print("  to the t=0 Wigner-function nodal lines.")
    print()
    print(f"  {'category':>14s}   {'(x_0, p_0)':>16s}   "
          f"{'(x_final, p_final)':>22s}")
    t_final = TIMES[-1]
    for cat, seeds in (("near node", near_node), ("non-node", non_node)):
        for s in seeds:
            xt = s[0] + s[1] * t_final / MASS
            xw = ((xt + L / 2) % L) - L / 2
            print(f"  {cat:>14s}   ({s[0]:+.3f}, {s[1]:+.3f})"
                  f"   ({xw:+.3f}, {s[1]:+.3f})")

    plot_evolution(snapshots, near_node, non_node,
                   "cat_state_microdynamics_evolution.png")
    plot_marginals(snapshots,
                   "cat_state_microdynamics_marginals.png")
    plot_trajectory_portrait(near_node, non_node,
                             "cat_state_microdynamics_trajectories.png")
    print("\nDone.")


if __name__ == "__main__":
    main()
