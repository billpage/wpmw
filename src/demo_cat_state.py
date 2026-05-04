"""
Demonstration: free evolution of a Schrödinger-cat (colliding Gaussian wave packets).

This is the canonical Wigner-function showpiece. Two Gaussian wave packets
with opposite momenta approach each other along the x-axis. The
*superposition* state |psi_+> + |psi_-> has a Wigner function with three
parts:

    W = (1/2N) [ W_++  +  W_--  +  2 W_int ]

where W_++- are the two Gaussian "lobes" and

    W_int(x, p) = (1/pi*hbar) exp(-x^2/sigma^2 - p^2 sigma^2/hbar^2)
                  * cos(2 (p_0 x + p x_0) / hbar)

is the interference term sitting between them in phase space — a Gaussian
envelope modulated by cosine fringes whose nodal lines satisfy

    p_0 x + p x_0 = (n + 1/2) pi*hbar/2.

For free evolution V = 0 the Quantum Liouville Equation reduces to the
classical Liouville equation

    dW/dt + (p/m) dW/dx = 0  =>  W(x, p, t) = W(x - p t/m, p, 0)

— a rigid horizontal shear in phase space. All Moyal corrections vanish
because they involve derivatives of V. Therefore:

  * Spectral advection of the crystal-lattice solver should reproduce the
    exact Wigner function to floating-point precision.
  * **Any nodal point of the Wigner function at t = 0 stays a node forever,
    transported along the classical characteristic** (x_s + p_s t/m, p_s)
    — a horizontal line at speed p_s/m.

This demo exercises both points: it compares the PSC evolution to the
closed-form Wigner function at four times, and it overlays sample
trajectories from points placed exactly on interference nodes at t = 0.

Outputs
-------
- ``cat_state_evolution.png`` — four-panel time series of W(x, p) snapshots
  (PSC vs exact), with sample trajectories from interference nodes overlaid.
- ``cat_state_marginals.png`` — position-space probability density |psi(x)|^2
  at the same times, showing the famous t = t_c interference pattern.

Run with:

    WPMW_OUTPUT=/mnt/user-data/outputs python -u demo_cat_state.py
"""

from __future__ import annotations

import sys
import time

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from phase_space_crystal_lattice import PhaseSpaceCrystalLattice
from wpmw_utils import output_path


# --------------------------------------------------------------------- #
# Physics                                                               #
# --------------------------------------------------------------------- #
HBAR = 1.0
MASS = 1.0

# Cat-state parameters.  Right packet at (+x0, -p0); left packet at (-x0, +p0).
# They collide at the origin at t_c = m x0 / p0.
X0    = 4.0
P0    = 2.0
SIGMA = 1.0
SIGMA_P = HBAR / (2.0 * SIGMA)        # minimum-uncertainty Gaussian
T_C   = MASS * X0 / P0                # collision time

# Grid.  Need:
#   * L >> 2 x0 (no boundary wrap of the lobes during the run)
#   * dx << pi hbar / p0 (resolve x-fringes in the interference term)
#   * dp << pi hbar / x0 (resolve p-fringes in the interference term)
#   * P_max > p0 + ~3 sigma_p (lobes don't wrap in p)
L     = 32.0
M_POS = 256
N_MOM = 256
DX    = L / M_POS
DP    = np.pi * HBAR / L
P_MAX = (N_MOM // 2) * DP


# --------------------------------------------------------------------- #
# Closed-form initial Wigner function                                   #
# --------------------------------------------------------------------- #
def cat_overlap() -> float:
    """<psi_+ | psi_-> for two Gaussian packets at +- x0 with opposite momenta +- p0.

    Returns a small positive real number which vanishes exponentially in the
    well-separated limit.
    """
    return float(np.exp(-X0 ** 2 / SIGMA ** 2 - P0 ** 2 * SIGMA ** 2 / HBAR ** 2))


def W_cat_initial(X, P) -> np.ndarray:
    """Closed-form Wigner function of the normalised cat state at t = 0.

    psi_+ centred at +x0 with momentum -p0  (moves left)
    psi_- centred at -x0 with momentum +p0  (moves right)
    psi   = (psi_+ + psi_-) / sqrt(2 (1 + alpha))
    """
    alpha = cat_overlap()
    norm = 1.0 / (np.pi * HBAR) / (2.0 * (1.0 + alpha))
    Wpp = np.exp(-((X - X0) ** 2) / SIGMA ** 2
                 - ((P + P0) ** 2) * SIGMA ** 2 / HBAR ** 2)
    Wmm = np.exp(-((X + X0) ** 2) / SIGMA ** 2
                 - ((P - P0) ** 2) * SIGMA ** 2 / HBAR ** 2)
    Wint = (np.exp(-X ** 2 / SIGMA ** 2 - P ** 2 * SIGMA ** 2 / HBAR ** 2)
            * np.cos(2.0 * (P0 * X + P * X0) / HBAR))
    return norm * (Wpp + Wmm + 2.0 * Wint)


def W_cat_exact(X, P, t) -> np.ndarray:
    """Exact Wigner function at time t under V = 0 evolution: W_0(x - pt/m, p)."""
    Xs = X - P * t / MASS
    return W_cat_initial(Xs, P)


# --------------------------------------------------------------------- #
# Sample trajectories from interference nodes at t = 0                  #
# --------------------------------------------------------------------- #
def sample_node_seeds() -> list[tuple[float, float]]:
    """Pick a handful of points lying on interference nodes at t = 0 — and
    well *inside* the central interference envelope, so the lobes contribute
    negligibly and the cosine fringes are the dominant structure.

    The cosine factor cos(2(p_0 x + p x_0)/hbar) vanishes when

        p_0 x + p x_0 = (k + 1/2) pi hbar / 2     (k integer)

    For p_0 = 2, x_0 = 4, hbar = 1 these are the lines x + 2 p = (k+1/2) pi/4.

    The interference envelope has standard deviations sigma in x and
    hbar/(2 sigma) in p; we pick seeds with |p| <= sigma_p so they remain
    well within the central envelope (and far from the lobes at (+-x_0, -+p_0)).
    """
    seeds: list[tuple[float, float]] = []

    # Two stationary seeds on adjacent nodal lines (k = -1 and k = +1) at p = 0
    # — these illustrate two distinct nodes that do NOT move.
    for k in (-1, +1):
        x_node = (k + 0.5) * np.pi * HBAR / (2.0 * P0)
        seeds.append((float(x_node), 0.0))

    # Three seeds with non-zero p on the central (k = 0) nodal line — these
    # have visibly different trajectory speeds and directions.  Constrain to
    # |p| <= sigma_p so we stay inside the envelope.
    rhs = 0.5 * np.pi * HBAR / 2.0     # = pi hbar / 4
    for p_s in (-0.4, 0.0, +0.4):
        x_s = (rhs - p_s * X0) / P0
        seeds.append((float(x_s), float(p_s)))
    return seeds


def trajectory(seed: tuple[float, float], t: float) -> tuple[float, float]:
    """Classical free-particle trajectory: (x(t), p(t)) = (x0 + p0 t/m, p0)."""
    xs, ps = seed
    return (xs + ps * t / MASS, ps)


# --------------------------------------------------------------------- #
# PSC run                                                               #
# --------------------------------------------------------------------- #
def run_psc(times: np.ndarray) -> dict:
    """Run the crystal-lattice solver (spectral advection, V = 0) and snapshot
    the Wigner function at each requested time."""
    print(f"\n[crystal-lattice cat-state run]")
    print(f"  x0 = {X0}, p0 = {P0}, sigma = {SIGMA}, t_collision = {T_C:.4f}")
    print(f"  M = {M_POS}, N = {N_MOM}, L = {L}, dx = {DX}, dp = {DP:.5f}")
    print(f"  P_max = {P_MAX:.4f}, max |p| in lobe ~ {P0 + 3*SIGMA_P:.4f}")
    print(f"  cat-state overlap alpha = {cat_overlap():.3e}")

    psc = PhaseSpaceCrystalLattice(
        M=M_POS, N=N_MOM, L=L, mass=MASS, hbar=HBAR,
        nu=None, advection="spectral",
    )
    psc.initialize_from_wigner(W_cat_initial)
    snapshots = {}

    # Step from t=0 through the requested times in order.  Spectral advection
    # is exact for any dt but use a uniform ladder so timings are comparable.
    target_times = sorted(set(float(t) for t in times))
    DT = 0.02
    n_steps_total = int(round(target_times[-1] / DT))
    print(f"  spectral advection, dt = {DT}, n_steps to t_max = {n_steps_total}")

    t = 0.0
    snapshots[0.0] = psc.W.copy()
    next_idx = 1
    if next_idx < len(target_times) and abs(target_times[next_idx-1]) > DT:
        next_idx = 0  # 0 was first target
    next_idx = 1 if target_times[0] == 0.0 else 0

    t_start = time.time()
    while next_idx < len(target_times):
        target = target_times[next_idx]
        # advance by integer number of DT steps to land near target
        steps = max(1, int(round((target - t) / DT)))
        for _ in range(steps):
            psc.step_advect(DT)
        t += steps * DT
        snapshots[target] = psc.W.copy()
        next_idx += 1
    wall = time.time() - t_start
    print(f"  wall = {wall:.2f} s")

    # Norm conservation check
    norms = [(np.sum(W) * DX * DP) for W in snapshots.values()]
    print(f"  norm: initial = {norms[0]:.10f}, final = {norms[-1]:.10f}, "
          f"drift = {norms[-1] - norms[0]:+.3e}")

    return {
        "x": psc.x,
        "p": psc.p,
        "X": psc.X,
        "P": psc.P,
        "snapshots": snapshots,
    }


# --------------------------------------------------------------------- #
# Plotting                                                              #
# --------------------------------------------------------------------- #
def plot_evolution(res: dict, savepath: str) -> None:
    """Three rows of four panels:
       row 1 — PSC W(x, p, t)
       row 2 — exact W(x, p, t)
       row 3 — pointwise difference (PSC - exact)
       Sample trajectories overlaid on rows 1 and 2.
    """
    times = sorted(res["snapshots"].keys())
    seeds = sample_node_seeds()
    X, P = res["X"], res["P"]

    # Symmetric color scale based on the t = 0 max
    vmax = float(np.max(np.abs(res["snapshots"][times[0]])))
    levels = np.linspace(-vmax, vmax, 21)

    fig, axes = plt.subplots(3, 4, figsize=(15, 9), sharex=True, sharey=True)

    for j, t in enumerate(times):
        W_psc = res["snapshots"][t]
        W_ex  = W_cat_exact(X, P, t)
        diff  = W_psc - W_ex
        diff_max = max(1e-14, float(np.max(np.abs(diff))))

        ax = axes[0, j]
        cf = ax.contourf(res["x"], res["p"], W_psc, levels=levels,
                         cmap="RdBu_r", extend="both")
        ax.set_title(f"PSC,  t = {t:.2f}", fontsize=10)
        ax.set_xlim(-8, 8)
        ax.set_ylim(-3.0, 3.0)
        if j == 0:
            ax.set_ylabel("p")

        ax = axes[1, j]
        ax.contourf(res["x"], res["p"], W_ex, levels=levels,
                    cmap="RdBu_r", extend="both")
        ax.set_title(f"exact,  t = {t:.2f}", fontsize=10)
        if j == 0:
            ax.set_ylabel("p")

        ax = axes[2, j]
        ax.contourf(res["x"], res["p"], diff,
                    levels=np.linspace(-diff_max, diff_max, 21),
                    cmap="PuOr", extend="both")
        ax.set_title(f"PSC − exact,  $\\max = {diff_max:.1e}$", fontsize=10)
        ax.set_xlabel("x")
        if j == 0:
            ax.set_ylabel("p")

        # Sample-trajectory overlay: full trajectory as a thin line, with an
        # open circle marking the current position.  Each trajectory is a
        # horizontal line at constant p (free-particle characteristic), of
        # length p_s * t in x — well-known feature of free Liouville flow.
        traj_colors = ["#222222"] * len(seeds)
        t_span = np.linspace(0.0, t, 50)
        for s, c in zip(seeds, traj_colors):
            xs_t = s[0] + s[1] * t_span / MASS
            xs_t_w = ((xs_t + L/2) % L) - L/2
            # Insert NaN at periodic-boundary jumps so the line doesn't
            # streak across the plot.
            xv, pv = [], []
            prev = None
            for xx in xs_t_w:
                if prev is not None and abs(xx - prev) > L/2:
                    xv.append(np.nan); pv.append(np.nan)
                xv.append(xx); pv.append(s[1])
                prev = xx
            for r in (0, 1):
                axes[r, j].plot(xv, pv, "-", color=c, lw=0.8, alpha=0.55)
                # Current-position marker
                xt, pt = trajectory(s, t)
                xw = ((xt + L/2) % L) - L/2
                axes[r, j].plot(xw, pt, "o", color=c, ms=4, mfc="none", mew=1.2)

    fig.suptitle(
        "Free-particle cat state: Wigner-function evolution\n"
        f"$\\sigma = {SIGMA}$, $x_0 = {X0}$, $p_0 = {P0}$, "
        f"$t_c = m x_0/p_0 = {T_C:.2f}$.  "
        "Open circles: classical trajectories from interference nodes at $t=0$.",
        fontsize=11, y=0.99,
    )
    fig.subplots_adjust(top=0.90, bottom=0.07, left=0.05, right=0.93,
                        wspace=0.15, hspace=0.30)
    cax = fig.add_axes([0.95, 0.40, 0.012, 0.45])
    cbar = fig.colorbar(cf, cax=cax)
    cbar.set_label("W")
    fig.savefig(savepath, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"\nFigure saved -> {savepath}")


def plot_marginals(res: dict, savepath: str) -> None:
    """Position-space probability density |psi(x, t)|^2 = integral W(x, p, t) dp."""
    times = sorted(res["snapshots"].keys())

    fig, axes = plt.subplots(1, len(times), figsize=(15, 3), sharey=True)

    for j, t in enumerate(times):
        W_psc = res["snapshots"][t]
        W_ex  = W_cat_exact(res["X"], res["P"], t)
        rho_psc = np.sum(W_psc, axis=0) * DP
        rho_ex  = np.sum(W_ex,  axis=0) * DP

        ax = axes[j]
        ax.plot(res["x"], rho_ex,  "k--", lw=1.5, label="exact $|\\psi|^2$")
        ax.plot(res["x"], rho_psc, "C0",  lw=1.0, label="PSC marginal")
        ax.set_title(f"t = {t:.2f}", fontsize=11)
        ax.set_xlabel("x")
        ax.set_xlim(-X0 - 4*SIGMA, X0 + 4*SIGMA)
        ax.grid(alpha=0.3)
        if j == 0:
            ax.set_ylabel("$\\rho(x)$")
            ax.legend(fontsize=9, loc="upper right")

    fig.suptitle(
        "Position-space probability density at successive times — "
        f"interference at the collision time $t_c = {T_C:.2f}$",
        fontsize=11, y=1.04,
    )
    fig.tight_layout()
    fig.savefig(savepath, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"Figure saved -> {savepath}")


# --------------------------------------------------------------------- #
# Main                                                                  #
# --------------------------------------------------------------------- #
def main() -> None:
    print("=" * 70)
    print("Free-particle cat state — Wigner-function evolution demo")
    print("=" * 70)

    # Snapshots: t=0, t_c/2, t_c, 3 t_c/2
    times = np.array([0.0, T_C / 2.0, T_C, 1.5 * T_C])
    res = run_psc(times)

    # Pointwise comparison summary
    print()
    print("Pointwise PSC vs exact comparison:")
    print(f"  {'t':>7s}   {'max|PSC - exact|':>20s}  {'L2 deviation':>15s}")
    for t in sorted(res["snapshots"].keys()):
        W_psc = res["snapshots"][t]
        W_ex  = W_cat_exact(res["X"], res["P"], t)
        d     = W_psc - W_ex
        l2    = float(np.sqrt(np.sum(d ** 2) * DX * DP))
        print(f"  {t:7.4f}   {float(np.max(np.abs(d))):20.3e}  {l2:15.3e}")

    # Sanity check: nodes from t = 0 stay nodes
    print()
    print("Sample-trajectory check: |W| at node trajectories (should stay near zero):")
    print(f"  {'seed (x0, p0)':>22s}   " + "  ".join(
        f"|W(t={t:.2f})|" for t in sorted(res["snapshots"].keys())))
    seeds = sample_node_seeds()
    for s in seeds:
        row = []
        for t in sorted(res["snapshots"].keys()):
            xt, pt = trajectory(s, t)
            xw = ((xt + L / 2) % L) - L / 2
            # Bilinear interpolation of W onto (xw, pt)
            i = int(np.clip((xw - res["x"][0]) / DX, 0, M_POS - 2))
            j = int(np.clip((pt - res["p"][0]) / DP, 0, N_MOM - 2))
            fx = (xw - res["x"][i]) / DX
            fp = (pt - res["p"][j]) / DP
            W = res["snapshots"][t]
            val = ((1 - fx) * (1 - fp) * W[j, i] + fx * (1 - fp) * W[j, i + 1]
                   + (1 - fx) * fp * W[j + 1, i] + fx * fp * W[j + 1, i + 1])
            row.append(f"{abs(val):.2e}")
        print(f"  ({s[0]:+.3f}, {s[1]:+.3f})   " + "    ".join(row))

    plot_evolution(res, output_path("cat_state_evolution.png"))
    plot_marginals(res, output_path("cat_state_marginals.png"))
    print("\nDone.")


if __name__ == "__main__":
    main()
