"""
Demonstration: quantum harmonic oscillator ground state.

Compares two solvers for the Wigner equation evolving the QHO ground-state
Wigner function:

  1. Phase-Space Crystal-Lattice (this project), using:
       - integer-roll free-streaming (spec §3a)
       - QLE differential form for the force-jump (spec §7b), which is exact
         in the continuum limit for quadratic V.
     Optionally also runs a 'spectral' free-streaming variant for fairness.

  2. Standard Strang-split spectral Fourier solver on the Wigner equation.
     For QHO this is exact up to time-step error.

The QHO ground-state Wigner function in natural units (hbar = m = omega = 1)
is the rotationally-symmetric Gaussian

    W_0(x, p) = (1/pi) exp(-x^2 - p^2)

which is a stationary state under the exact QLE flow. A successful solver
must preserve it.

The demo:
  - Initializes both solvers with W_0 on a shared grid.
  - Evolves for one full classical period T = 2 pi.
  - Reports:
      * conservation of total norm
      * L2 deviation of W(t) from W_0
      * marginals rho(x), rho(p) at t = 0, T/4, T/2, T
      * a 6-panel summary figure saved via output_path()

Run from the src/ directory with:

    WPMW_OUTPUT=/mnt/user-data/outputs python -u demo_qho_ground_state.py
"""

from __future__ import annotations

import sys
import time

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from wpmwlib.phase_space_crystal_lattice import PhaseSpaceCrystalLattice
from wpmwlib.wigner_split_fourier import WignerSplitFourier
from wpmwlib.wpmw_utils import output_path


# --------------------------------------------------------------------- #
# Physics and grid                                                      #
# --------------------------------------------------------------------- #
HBAR = 1.0
MASS = 1.0
OMEGA = 1.0
T_PERIOD = 2.0 * np.pi / OMEGA          # = 2 pi for omega = 1

# Phase-space grid. Choose L so dp = pi*hbar/L gives a clean grid.
L = 16.0
M_POS = 128                              # x cells
N_MOM = 128                              # p cells
DP = np.pi * HBAR / L                    # = pi/16  ~= 0.196
P_MAX = (N_MOM // 2) * DP                # = 4 pi   ~= 12.57
DX = L / M_POS                           # = 0.125

# Time stepping. For integer-roll advection, want max shift to be a small
# integer so the scheme is well-defined. Max shift = round(P_MAX * dt / dx).
# Set dt so that max shift = K position cells for some small integer K.
K_SHIFT = 1
DT = K_SHIFT * DX / P_MAX                # = 0.125 / (4 pi) ~= 0.00994
N_STEPS = int(round(T_PERIOD / DT))      # ~ 632 steps for one period


def W_qho_ground(X: np.ndarray, P: np.ndarray) -> np.ndarray:
    """QHO ground-state Wigner function in natural units."""
    return (1.0 / (np.pi * HBAR)) * np.exp(
        -MASS * OMEGA * X ** 2 / HBAR - P ** 2 / (MASS * OMEGA * HBAR)
    )


def dVdx(x: np.ndarray) -> np.ndarray:
    """Force gradient: dV/dx for V = (1/2) m omega^2 x^2."""
    return MASS * OMEGA ** 2 * x


# --------------------------------------------------------------------- #
# Run a solver, recording snapshots and diagnostics                     #
# --------------------------------------------------------------------- #
def run_crystal_lattice(advection: str, label: str):
    """Run the crystal-lattice solver using the differential jump form."""
    print(f"\n[{label}] starting (advection = {advection})")
    psc = PhaseSpaceCrystalLattice(
        M=M_POS, N=N_MOM, L=L, mass=MASS, hbar=HBAR,
        nu=None, advection=advection,
    )
    psc.initialize_from_wigner(W_qho_ground)
    W0 = psc.get_wigner().copy()
    norm0 = psc.total_norm()

    snapshots = {0.0: W0.copy()}
    target_times = [T_PERIOD * f for f in (0.25, 0.5, 0.75, 1.0)]
    next_idx = 0

    dVdx_arr = dVdx(psc.x)
    norms = [norm0]
    l2errs = [0.0]

    t0 = time.time()
    t = 0.0
    for step in range(1, N_STEPS + 1):
        psc.strang_step_differential(dVdx_arr, DT)
        t = step * DT
        norms.append(psc.total_norm())
        l2errs.append(np.sqrt(np.sum((psc.get_wigner() - W0) ** 2) * DX * DP))
        # Capture snapshots near the target times
        while next_idx < len(target_times) and t >= target_times[next_idx] - DT * 0.5:
            snapshots[target_times[next_idx]] = psc.get_wigner().copy()
            next_idx += 1
    wall = time.time() - t0

    print(f"[{label}] {N_STEPS} steps, dt = {DT:.5f}, wall = {wall:.2f} s")
    print(f"[{label}] norm: initial = {norm0:.10f}, final = {norms[-1]:.10f}, "
          f"drift = {norms[-1] - norm0:+.3e}")
    print(f"[{label}] L2 deviation from W0 at T = {l2errs[-1]:.3e}")

    return {
        "label": label,
        "x": psc.x, "p": psc.p,
        "snapshots": snapshots,
        "W0": W0,
        "norms": np.array(norms),
        "l2errs": np.array(l2errs),
        "psc": psc,
    }


def run_split_fourier(label: str):
    """Run the reference split-Fourier solver."""
    print(f"\n[{label}] starting")
    sf = WignerSplitFourier(
        M=M_POS, N=N_MOM, L=L, dp=DP, mass=MASS, hbar=HBAR,
    )
    sf.initialize_from_wigner(W_qho_ground)
    W0 = sf.get_wigner().copy()
    norm0 = sf.total_norm()

    snapshots = {0.0: W0.copy()}
    target_times = [T_PERIOD * f for f in (0.25, 0.5, 0.75, 1.0)]
    next_idx = 0

    norms = [norm0]
    l2errs = [0.0]
    imag_max = [sf.imag_residue()]

    t0 = time.time()
    t = 0.0
    for step in range(1, N_STEPS + 1):
        sf.strang_step_qho(OMEGA, DT)
        t = step * DT
        norms.append(sf.total_norm())
        l2errs.append(np.sqrt(np.sum((sf.get_wigner() - W0) ** 2) * DX * DP))
        imag_max.append(sf.imag_residue())
        while next_idx < len(target_times) and t >= target_times[next_idx] - DT * 0.5:
            snapshots[target_times[next_idx]] = sf.get_wigner().copy()
            next_idx += 1
    wall = time.time() - t0

    print(f"[{label}] {N_STEPS} steps, dt = {DT:.5f}, wall = {wall:.2f} s")
    print(f"[{label}] norm: initial = {norm0:.10f}, final = {norms[-1]:.10f}, "
          f"drift = {norms[-1] - norm0:+.3e}")
    print(f"[{label}] L2 deviation from W0 at T = {l2errs[-1]:.3e}")
    print(f"[{label}] max imaginary residue over run = {max(imag_max):.3e}")

    return {
        "label": label,
        "x": sf.x, "p": sf.p,
        "snapshots": snapshots,
        "W0": W0,
        "norms": np.array(norms),
        "l2errs": np.array(l2errs),
        "sf": sf,
    }


# --------------------------------------------------------------------- #
# Plotting                                                              #
# --------------------------------------------------------------------- #
def plot_summary(results, savepath):
    """Six-panel summary: marginals at t=T, L2 error, norm drift, and contour."""
    fig, axes = plt.subplots(2, 3, figsize=(13, 8))
    colors = ["C0", "C1", "C2"]

    # (a) Position marginal at t = T
    ax = axes[0, 0]
    x = results[0]["x"]
    rho_x_exact = np.sqrt(MASS * OMEGA / (np.pi * HBAR)) * np.exp(
        -MASS * OMEGA * x ** 2 / HBAR
    )
    ax.plot(x, rho_x_exact, "k--", lw=1.5, label="exact ground state")
    for r, c in zip(results, colors):
        rho = np.sum(r["snapshots"][T_PERIOD], axis=0) * DP
        ax.plot(x, rho, color=c, lw=1.0, label=r["label"])
    ax.set_xlim(-5, 5)
    ax.set_xlabel("x"); ax.set_ylabel(r"$\rho(x)$")
    ax.set_title(f"Position marginal at t = T = {T_PERIOD:.3f}")
    ax.legend(fontsize=8); ax.grid(alpha=0.3)

    # (b) Momentum marginal at t = T
    ax = axes[0, 1]
    p = results[0]["p"]
    rho_p_exact = np.sqrt(1.0 / (np.pi * MASS * OMEGA * HBAR)) * np.exp(
        -p ** 2 / (MASS * OMEGA * HBAR)
    )
    ax.plot(p, rho_p_exact, "k--", lw=1.5, label="exact ground state")
    for r, c in zip(results, colors):
        rho = np.sum(r["snapshots"][T_PERIOD], axis=1) * DX
        ax.plot(p, rho, color=c, lw=1.0, label=r["label"])
    ax.set_xlim(-5, 5)
    ax.set_xlabel("p"); ax.set_ylabel(r"$\rho(p)$")
    ax.set_title(f"Momentum marginal at t = T")
    ax.legend(fontsize=8); ax.grid(alpha=0.3)

    # (c) L2 deviation from W_0 over time
    ax = axes[0, 2]
    t_arr = np.arange(N_STEPS + 1) * DT
    for r, c in zip(results, colors):
        ax.plot(t_arr, r["l2errs"], color=c, lw=1.0, label=r["label"])
    ax.set_xlabel("t"); ax.set_ylabel(r"$\|W(t) - W_0\|_2$")
    ax.set_title("L2 deviation from initial state")
    ax.set_yscale("log"); ax.legend(fontsize=8); ax.grid(alpha=0.3)

    # (d) Norm drift over time
    ax = axes[1, 0]
    for r, c in zip(results, colors):
        ax.plot(t_arr, r["norms"] - r["norms"][0], color=c, lw=1.0, label=r["label"])
    ax.set_xlabel("t"); ax.set_ylabel(r"$\int W\,dx\,dp - 1$")
    ax.set_title("Norm conservation drift")
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    ax.ticklabel_format(axis="y", style="sci", scilimits=(-3, 3))

    # (e) Contour of W at t = T from the first solver (crystal-lattice)
    r = results[0]
    ax = axes[1, 1]
    W_T = r["snapshots"][T_PERIOD]
    levels = np.linspace(-0.05, 0.32, 12)
    cf = ax.contourf(r["x"], r["p"], W_T, levels=levels, cmap="RdBu_r", extend="both")
    ax.set_xlim(-4, 4); ax.set_ylim(-4, 4)
    ax.set_xlabel("x"); ax.set_ylabel("p")
    ax.set_title(f"{r['label']}\n$W(x, p, T)$")
    ax.set_aspect("equal")
    plt.colorbar(cf, ax=ax, fraction=0.046, pad=0.04)

    # (f) Contour of W at t = T from the reference (split-Fourier)
    r = results[-1]
    ax = axes[1, 2]
    W_T = r["snapshots"][T_PERIOD]
    cf = ax.contourf(r["x"], r["p"], W_T, levels=levels, cmap="RdBu_r", extend="both")
    ax.set_xlim(-4, 4); ax.set_ylim(-4, 4)
    ax.set_xlabel("x"); ax.set_ylabel("p")
    ax.set_title(f"{r['label']}\n$W(x, p, T)$")
    ax.set_aspect("equal")
    plt.colorbar(cf, ax=ax, fraction=0.046, pad=0.04)

    fig.suptitle(
        "QHO ground-state preservation: crystal-lattice vs split-Fourier reference",
        fontsize=13, y=1.00,
    )
    fig.tight_layout()
    fig.savefig(savepath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\nFigure saved -> {savepath}")


# --------------------------------------------------------------------- #
# Main                                                                  #
# --------------------------------------------------------------------- #
def main():
    print("=" * 64)
    print(f"QHO ground-state demo")
    print(f"  hbar = {HBAR}, m = {MASS}, omega = {OMEGA}, T = {T_PERIOD:.4f}")
    print(f"  M = {M_POS}, N = {N_MOM}, L = {L}, dx = {DX}, dp = {DP:.5f}")
    print(f"  P_max = {P_MAX:.4f}")
    print(f"  dt = {DT:.6f}, N_steps = {N_STEPS} (max integer-roll shift = {K_SHIFT})")
    print("=" * 64)

    results = []
    results.append(run_crystal_lattice("integer_roll", "Crystal-lattice (integer-roll)"))
    results.append(run_crystal_lattice("spectral", "Crystal-lattice (spectral advect)"))
    results.append(run_split_fourier("Split-Fourier reference"))

    print("\n" + "=" * 64)
    print("Summary at t = T")
    print("=" * 64)
    for r in results:
        n = r["norms"]
        print(f"  {r['label']:40s}  "
              f"norm drift = {n[-1] - n[0]:+.3e}  "
              f"L2 dev = {r['l2errs'][-1]:.3e}")

    plot_summary(results, output_path("qho_ground_state_comparison.png"))
    print("\nDone.")


if __name__ == "__main__":
    main()
