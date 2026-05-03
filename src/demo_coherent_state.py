"""
Demonstration: coherent state in a cosine potential.

Exercises the corrected sign convention (see ``sign_convention_check.py`` and
``docs/supplement/phase_space_crystal_lattice_supplement.md`` §6.3) on a
state that is sensitive to it. The QHO ground-state demo is *insensitive* to
the sign because its initial Wigner function is rotationally symmetric in
phase space; this demo uses a coherent state placed off-axis.

Setup
-----
Cosine potential V(x) = +V_max cos(2 pi x / L) on [-L/2, L/2):

    - Maximum (hill) at x = 0
    - Minima (wells) at x = +-L/2

The simulation initialises a coherent Gaussian Wigner state at
(x_0 = +L/4, p_0 = 0) — on the downhill slope of the central hill — and
evolves it with the crystal-lattice solver (Strang split: half-streaming
plus Fourier-mode jump plus half-streaming). The classical centroid
trajectory is obtained by integrating Hamilton's equations
``dot x = p/m,  dot p = -V'(x)`` and overlaid as a reference.

For the linearised motion at small displacement, V' ≈ -V_max (2 pi/L)^2 x,
i.e. an effective harmonic-oscillator with omega_eff = (2 pi/L) sqrt(V_max/m).
We choose V_max so that the small-oscillation period is comparable to the
total simulation time and the centroid completes a partial orbit.

Note: this demo uses *spectral* free-streaming rather than the spec's
integer-roll scheme. With integer-roll, position advection is quantised in
multiples of dx/dt and momenta below this floor produce zero shift — which
freezes < x > for the early portion of any low-momentum trajectory and
muddies the comparison with the smooth classical curve. The sign-sensitive
question (does the centroid first move *downhill*?) is preserved by either
choice of advection; the regression test in ``sign_convention_check.py``
uses jump-only updates and so is independent of the advection scheme.

Outputs
-------
- ``coherent_state_trajectory.png`` — centroid trajectory in phase space,
  time-series of < x > and < p >, and W(x, p) snapshots at four times.

Run with:

    WPMW_OUTPUT=/mnt/user-data/outputs python -u demo_coherent_state.py
"""

from __future__ import annotations

import sys
import time

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from phase_space_crystal_lattice import PhaseSpaceCrystalLattice, FourierMode
from wpmw_utils import output_path


# --------------------------------------------------------------------- #
# Physics and grid                                                      #
# --------------------------------------------------------------------- #
HBAR = 1.0
MASS = 1.0
L     = 16.0
V_MAX = 0.5                              # potential strength

# Effective small-oscillation frequency and period at small x
OMEGA_EFF = (2.0 * np.pi / L) * np.sqrt(V_MAX / MASS)
T_EFF     = 2.0 * np.pi / OMEGA_EFF      # ~22.6 for these params

# Grid
M_POS = 128
N_MOM = 128
DX    = L / M_POS
DP    = np.pi * HBAR / L
P_MAX = (N_MOM // 2) * DP

# Initial state — coherent Gaussian at downhill slope
X0    = L / 4.0
P0    = 0.0
SIGMA = 1.0                              # broader than ground state to be visible

# Time stepping. For integer-roll advection, the maximum integer shift is
# round(P_MAX * dt / DX). Pick dt so this is K = 1.
K_SHIFT = 1
DT      = K_SHIFT * DX / P_MAX
T_FINAL = T_EFF * 0.5                    # half a period — half-orbit
N_STEPS = int(round(T_FINAL / DT))


def W_coherent(X, P):
    """Coherent Gaussian centred at (X0, P0) with std SIGMA in x and 1/(2*SIGMA) in p (minimum-uncertainty)."""
    sigma_x = SIGMA
    sigma_p = HBAR / (2.0 * sigma_x)
    norm = 1.0 / (2.0 * np.pi * sigma_x * sigma_p)
    return norm * np.exp(
        -((X - X0) ** 2) / (2.0 * sigma_x ** 2)
        - ((P - P0) ** 2) / (2.0 * sigma_p ** 2)
    )


def dVdx(x):
    return -V_MAX * (2.0 * np.pi / L) * np.sin(2.0 * np.pi * x / L)


# --------------------------------------------------------------------- #
# Classical centroid by Hamilton's equations                            #
# --------------------------------------------------------------------- #
def classical_trajectory(x0, p0, t_arr, m=MASS):
    """Symplectic-Euler integration of dot x = p/m, dot p = -V'(x)."""
    xs = np.empty_like(t_arr)
    ps = np.empty_like(t_arr)
    x, p = float(x0), float(p0)
    xs[0], ps[0] = x, p
    for k in range(1, len(t_arr)):
        dt_k = t_arr[k] - t_arr[k - 1]
        # Drift-kick-drift
        x += 0.5 * dt_k * p / m
        p += -dt_k * dVdx(np.array([x]))[0]
        x += 0.5 * dt_k * p / m
        xs[k], ps[k] = x, p
    return xs, ps


# --------------------------------------------------------------------- #
# Run the crystal-lattice solver                                        #
# --------------------------------------------------------------------- #
def run_crystal_lattice():
    print(f"\n[crystal-lattice coherent-state run]")
    print(f"  V_max = {V_MAX}, omega_eff = {OMEGA_EFF:.4f}, T_eff = {T_EFF:.4f}")
    print(f"  initial state at (x0, p0) = ({X0}, {P0}), sigma_x = {SIGMA}")
    print(f"  M = {M_POS}, N = {N_MOM}, L = {L}, dx = {DX}, dp = {DP:.5f}")
    print(f"  dt = {DT:.5f}, N_steps = {N_STEPS}, t_final = {T_FINAL:.4f}")

    psc = PhaseSpaceCrystalLattice(
        M=M_POS, N=N_MOM, L=L, mass=MASS, hbar=HBAR,
        nu=None, advection="spectral",
    )
    psc.initialize_from_wigner(W_coherent)

    # The cosine potential is exactly a single Fourier mode q=1, V_q=V_max, phi=0.
    modes = [FourierMode(q=1, V_q=V_MAX, phi_q=0.0)]

    snapshots = {0.0: psc.get_wigner().copy()}
    target_frac = (0.25, 0.5, 0.75, 1.0)
    target_t    = [T_FINAL * f for f in target_frac]
    next_idx = 0

    centroid_x = [np.sum(psc.X * psc.get_wigner()) * DX * DP]
    centroid_p = [np.sum(psc.P * psc.get_wigner()) * DX * DP]
    norms      = [psc.total_norm()]

    t_start = time.time()
    for step in range(1, N_STEPS + 1):
        psc.strang_step_fourier(modes, DT)
        t = step * DT
        W = psc.get_wigner()
        centroid_x.append(np.sum(psc.X * W) * DX * DP)
        centroid_p.append(np.sum(psc.P * W) * DX * DP)
        norms.append(psc.total_norm())
        while next_idx < len(target_t) and t >= target_t[next_idx] - DT * 0.5:
            snapshots[target_t[next_idx]] = W.copy()
            next_idx += 1
    wall = time.time() - t_start

    print(f"  wall = {wall:.2f} s")
    print(f"  norm: initial = {norms[0]:.10f}, final = {norms[-1]:.10f}, "
          f"drift = {norms[-1] - norms[0]:+.3e}")
    print(f"  centroid: initial (x, p) = ({centroid_x[0]:+.4f}, {centroid_p[0]:+.4f})")
    print(f"  centroid: final   (x, p) = ({centroid_x[-1]:+.4f}, {centroid_p[-1]:+.4f})")

    # The decisive test: did <p> first move POSITIVE (downhill) as classical
    # mechanics requires? If the sign convention were wrong, <p> would go
    # negative immediately.
    early_step = max(1, N_STEPS // 20)
    p_early = centroid_p[early_step]
    print(f"  early <p> at step {early_step}: {p_early:+.4f}  "
          f"({'+' if p_early > 0 else '-'}; "
          f"{'PASS — downhill (correct)' if p_early > 0 else 'FAIL — uphill (sign bug)'})")

    return {
        "x":      psc.x, "p": psc.p,
        "X":      psc.X, "P": psc.P,
        "snapshots": snapshots,
        "centroid_x": np.array(centroid_x),
        "centroid_p": np.array(centroid_p),
        "norms":   np.array(norms),
        "t_arr":   np.arange(N_STEPS + 1) * DT,
        "early_p": p_early,
    }


# --------------------------------------------------------------------- #
# Plotting                                                              #
# --------------------------------------------------------------------- #
def plot_summary(res, savepath):
    t_arr = res["t_arr"]
    cx, cp = res["centroid_x"], res["centroid_p"]
    cls_x_raw, cls_p = classical_trajectory(X0, P0, t_arr)
    # Wrap classical x into the periodic simulation domain [-L/2, L/2) for
    # display, so it lines up with the (necessarily periodic) lattice
    # centroid.  The wrap creates jumps at the boundary.
    cls_x = ((cls_x_raw + L/2) % L) - L/2

    fig = plt.figure(figsize=(13, 9))
    gs = fig.add_gridspec(3, 4, hspace=0.45, wspace=0.35)

    # (a) centroid trajectory in phase space
    ax = fig.add_subplot(gs[0:2, 0:2])
    # Insert NaN at periodic-boundary jumps so the dashed classical line
    # doesn't draw a spurious horizontal segment across the plot.
    plot_x, plot_p = [], []
    for i in range(len(cls_x)):
        if i > 0 and abs(cls_x[i] - cls_x[i - 1]) > L / 2:
            plot_x.append(np.nan); plot_p.append(np.nan)
        plot_x.append(cls_x[i]); plot_p.append(cls_p[i])
    ax.plot(plot_x, plot_p, "k--", lw=1.5, alpha=0.7, label="classical (wrapped)")
    ax.plot(cx, cp, "C0-", lw=1.5, label="crystal-lattice centroid")
    ax.plot([X0], [P0], "go", ms=8, label="t = 0")
    ax.plot([cx[-1]], [cp[-1]], "rs", ms=8, label=f"t = {T_FINAL:.2f}")
    ax.set_xlim(-L/2, L/2); ax.set_xlabel("$\\langle x\\rangle$")
    ax.set_ylabel("$\\langle p\\rangle$")
    ax.set_title("Centroid trajectory in phase space")
    ax.legend(fontsize=9, loc="best")
    ax.grid(alpha=0.3)
    ax.axhline(0, color="gray", lw=0.5)
    ax.axvline(0, color="gray", lw=0.5)

    # (b) <x>(t) and <p>(t)
    ax = fig.add_subplot(gs[0, 2:])
    ax.plot(t_arr, cls_x_raw, "k--", lw=1.0, alpha=0.7, label="classical (unwrapped)")
    ax.plot(t_arr, cx, "C0", lw=1.2, label="crystal-lattice")
    ax.set_ylabel("$\\langle x\\rangle$")
    ax.set_title("Centroid components vs time")
    ax.legend(fontsize=9); ax.grid(alpha=0.3)
    ax.axhline(0, color="gray", lw=0.5)
    ax.axhline( L/2, color="gray", lw=0.3, ls=":")
    ax.axhline(-L/2, color="gray", lw=0.3, ls=":")

    ax = fig.add_subplot(gs[1, 2:])
    ax.plot(t_arr, cls_p, "k--", lw=1.0, alpha=0.7)
    ax.plot(t_arr, cp, "C0", lw=1.2)
    ax.set_xlabel("t"); ax.set_ylabel("$\\langle p\\rangle$")
    ax.grid(alpha=0.3)
    ax.axhline(0, color="gray", lw=0.5)

    # (c-f) W(x, p) at four snapshot times
    snap_times = sorted(res["snapshots"].keys())
    levels = np.linspace(-0.02, 0.16, 12)
    for k, t in enumerate(snap_times[:4]):
        ax = fig.add_subplot(gs[2, k])
        W = res["snapshots"][t]
        cf = ax.contourf(res["x"], res["p"], W, levels=levels,
                         cmap="RdBu_r", extend="both")
        ax.set_xlim(-L/2, L/2); ax.set_ylim(-3, 3)
        ax.set_aspect("equal")
        ax.set_xlabel("x"); ax.set_ylabel("p" if k == 0 else "")
        ax.set_title(f"$W(x,p)$ at $t = {t:.2f}$", fontsize=10)
        idx = int(round(t / DT))
        ax.plot(cx[:idx+1], cp[:idx+1], "C0-", lw=1.2)
    cbar = fig.colorbar(cf, ax=fig.axes[-4:], shrink=0.6, location="right")
    cbar.set_label("W")

    fig.suptitle(
        f"Coherent state in cosine potential V = $V_{{\\max}}\\,\\cos(2\\pi x/L)$ — "
        f"$V_{{\\max}} = {V_MAX}$, $T_{{eff}} = {T_EFF:.2f}$",
        fontsize=12, y=0.99,
    )
    fig.savefig(savepath, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"\nFigure saved -> {savepath}")


# --------------------------------------------------------------------- #
# Main                                                                  #
# --------------------------------------------------------------------- #
def main():
    print("=" * 64)
    print("Coherent state in cosine potential — sign-sensitive demo")
    print("=" * 64)
    res = run_crystal_lattice()
    plot_summary(res, output_path("coherent_state_trajectory.png"))

    print("\n" + "=" * 64)
    if res["early_p"] > 0:
        print("RESULT: <p> first moves DOWNHILL (positive) — sign convention OK.")
    else:
        print("RESULT: <p> first moves UPHILL (negative) — SIGN BUG STILL PRESENT.")
        sys.exit(1)
    print("Done.")


if __name__ == "__main__":
    main()
