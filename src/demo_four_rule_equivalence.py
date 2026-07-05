"""
Demonstration: equivalence of the four-rule two-body scheme and the
single-rule mediated-jump scheme.

Companion to ``docs/analysis/four_rule_microdynamics_equivalence.md``.

The original crystal-lattice algorithm
(``docs/algorithm/phase_space_crystal_lattice_algorithm.md`` §3b) drives the
QLE collision term with a **single rule**: a positon at momentum cell n
mediates the transfer of one positon from cell n+q to cell n-q (or the
reverse, by the sign of Gamma_q(x)). The proposed alternative (D. Cyganski,
Zoom presentation, 2026) replaces this with **four rules** — Focus, Defocus,
Right-Hop, Left-Hop — acting on the momentum-cell triple (n-q, n, n+q) with
signed, occupancy-dependent bias rates.

This demo provides the numerical half of the equivalence argument:

  A. **Generator identity check.** On random mesh Wigner fields and for
     several Fourier-mode sets, a single Euler step of the four-rule
     mesh form (``step_jump_four_rule``, channels assembled independently)
     is compared against the original stencil (``step_jump_fourier``).
     The two must agree to machine precision — this is the discrete
     identity proved algebraically in the analysis note.

  B. **Head-to-head stochastic evolution.** A squeezed Gaussian in the
     single-period cosine well V(x) = -V_p cos(2 pi x / L) is evolved
     three ways on the same lattice, with an *identical* splitting
     sequence (exact-integer advection macro-steps, s jump sub-steps):

       1. deterministic mesh QLE (the common target),
       2. single-rule Monte Carlo (``step_jump_fourier_mc``),
       3. four-rule Monte Carlo (``step_jump_four_rule_mc``).

     Both MC runs must track the mesh within their shot-noise floors,
     and hence each other. Because the four-rule MC computes its bias
     rates from the *excess* (background-subtracted) counts, it fires
     far fewer events per step than the original rule (whose mediator
     count includes the uniform 2/h background); the metrics figure
     shows the resulting difference in noise floor at equal nu.

Splitting note
--------------
The macro time step is chosen as dt_adv = m dx / dp so the spec §3a
integer-roll advection is *exact* (row n shifts by exactly its momentum
index). Jumps are sub-cycled s times per macro step so |Gamma| dt_jump
stays small. All three evolutions use the same sequence, so residual
splitting error is common mode and differences isolate the jump-term
implementations, which is the object under test.
"""

from __future__ import annotations

import sys
import time

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from wpmwlib.phase_space_crystal_lattice import (  # noqa: E402
    FourierMode,
    PhaseSpaceCrystalLattice,
)
from wpmwlib.wpmw_utils import output_path, docs_path  # noqa: E402

# --------------------------------------------------------------------------
# Shared parameters
# --------------------------------------------------------------------------
HBAR = 1.0
MASS = 1.0
L = 8.0
V_P = 1.5
PHI = np.pi                    # V(x) = V_p cos(2 pi x/L + pi) = -V_p cos(2 pi x/L)
M_CELLS = 64
N_CELLS = 64

OMEGA = (2.0 * np.pi / L) * np.sqrt(V_P / MASS)
T_PERIOD = 2.0 * np.pi / OMEGA
SIGMA_X_GS = np.sqrt(HBAR / (2.0 * MASS * OMEGA))
SIGMA_X0 = 1.3 * SIGMA_X_GS            # mild squeeze -> breathing + negativity
SIGMA_P0 = HBAR / (2.0 * SIGMA_X0)     # min-uncertainty partner width

NU = 1600000                          # particles per unit phase-space volume
SUBSTEPS = 16                          # jump sub-steps per advection macro-step
T_FINAL = 2.0 * T_PERIOD
SEED = 20260705

MODES = [FourierMode(q=1, V_q=V_P, phi_q=PHI)]


def W_initial(X: np.ndarray, P: np.ndarray) -> np.ndarray:
    """Min-uncertainty squeezed Gaussian at the well bottom."""
    return (1.0 / (np.pi * HBAR)) * np.exp(
        -(X ** 2) / (2.0 * SIGMA_X0 ** 2) - (P ** 2) / (2.0 * SIGMA_P0 ** 2)
    )


# --------------------------------------------------------------------------
# Part A: generator identity on random fields
# --------------------------------------------------------------------------
def random_smooth_field(rng: np.random.Generator, N: int, M: int) -> np.ndarray:
    """Real random band-limited field (no smoothness actually required —
    the identity is exact for arbitrary fields — but band-limited data keeps
    magnitudes O(1) for a clean relative-error report)."""
    kmax = 6
    F = np.zeros((N, M), dtype=complex)
    for kn in range(-kmax, kmax + 1):
        for km in range(-kmax, kmax + 1):
            F[kn % N, km % M] = rng.standard_normal() + 1j * rng.standard_normal()
    W = np.fft.ifft2(F).real
    return W / np.abs(W).max()


def part_a() -> float:
    rng = np.random.default_rng(SEED)
    mode_sets = {
        "q=1 (canonical)": [FourierMode(1, 1.5, np.pi)],
        "q=2": [FourierMode(2, 0.7, 0.3)],
        "q=1,2,3 mixed": [
            FourierMode(1, 1.5, np.pi),
            FourierMode(2, 0.7, 0.3),
            FourierMode(3, 0.4, 1.1),
        ],
    }
    dt = 0.01
    worst = 0.0
    print("Part A: generator identity, four-rule mesh vs single-rule stencil")
    print(f"  {'mode set':<18} {'trial':>5} {'max |dW4 - dW1|':>18} {'rel':>12}")
    for name, modes in mode_sets.items():
        for trial in range(3):
            W0 = random_smooth_field(rng, N_CELLS, M_CELLS)
            a = PhaseSpaceCrystalLattice(M_CELLS, N_CELLS, L, MASS, HBAR)
            b = PhaseSpaceCrystalLattice(M_CELLS, N_CELLS, L, MASS, HBAR)
            a.W = W0.copy()
            b.W = W0.copy()
            a.step_jump_fourier(modes, dt)
            b.step_jump_four_rule(modes, dt)
            err = np.abs(a.W - b.W).max()
            scale = np.abs(a.W - W0).max()
            rel = err / scale if scale > 0 else 0.0
            worst = max(worst, err)
            print(f"  {name:<18} {trial:>5} {err:>18.3e} {rel:>12.3e}")
    print(f"  worst absolute deviation: {worst:.3e}")
    assert worst < 1e-12, "four-rule mesh form deviates from the QLE stencil"
    return worst


# --------------------------------------------------------------------------
# Part B: head-to-head stochastic evolution in the cosine well
# --------------------------------------------------------------------------
def evolve() -> dict:
    # Exact-integer advection macro step (spec §3a with zero rounding error).
    dx = L / M_CELLS
    dp = np.pi * HBAR / L
    dt_adv = MASS * dx / dp
    n_macro = int(round(T_FINAL / dt_adv))
    dt_jump = dt_adv / SUBSTEPS
    gmax = V_P / HBAR
    print("\nPart B: cosine-well evolution")
    print(f"  macro steps = {n_macro}, dt_adv = {dt_adv:.4f}, "
          f"substeps = {SUBSTEPS}, |Gamma|max dt_jump = {gmax * dt_jump:.4f}")

    mesh = PhaseSpaceCrystalLattice(M_CELLS, N_CELLS, L, MASS, HBAR, nu=None,
                                    advection="integer_roll")
    mc1 = PhaseSpaceCrystalLattice(M_CELLS, N_CELLS, L, MASS, HBAR, nu=NU,
                                   advection="integer_roll")
    mc4 = PhaseSpaceCrystalLattice(M_CELLS, N_CELLS, L, MASS, HBAR, nu=NU,
                                   advection="integer_roll")
    for solver in (mesh, mc1, mc4):
        solver.initialize_from_wigner(W_initial)

    rng1 = np.random.default_rng(SEED + 1)
    rng4 = np.random.default_rng(SEED + 2)

    snap_steps = sorted({0, n_macro // 4, n_macro // 2, n_macro})
    snaps = {k: {} for k in ("mesh", "mc1", "mc4")}
    times, err1, err4, neg = [], [], [], {"mesh": [], "mc1": [], "mc4": []}

    def record(step: int):
        t = step * dt_adv
        Wm, W1, W4 = mesh.get_wigner(), mc1.get_wigner(), mc4.get_wigner()
        norm = np.sqrt((Wm ** 2).sum())
        times.append(t)
        err1.append(np.sqrt(((W1 - Wm) ** 2).sum()) / norm)
        err4.append(np.sqrt(((W4 - Wm) ** 2).sum()) / norm)
        for key, W in (("mesh", Wm), ("mc1", W1), ("mc4", W4)):
            neg[key].append(-np.minimum(W, 0.0).sum() * mesh.dx * mesh.dp)
        if step in snap_steps:
            snaps["mesh"][t] = Wm.copy()
            snaps["mc1"][t] = W1.copy()
            snaps["mc4"][t] = W4.copy()

    record(0)
    t0 = time.time()
    for step in range(1, n_macro + 1):
        for solver in (mesh, mc1, mc4):
            solver.step_advect(dt_adv)
        for _ in range(SUBSTEPS):
            mesh.step_jump_fourier(MODES, dt_jump)
            mc1.step_jump_fourier_mc(MODES, dt_jump, rng=rng1)
            mc4.step_jump_four_rule_mc(MODES, dt_jump, rng=rng4)
        record(step)
        if step % 10 == 0 or step == n_macro:
            print(f"  step {step:3d}/{n_macro}  t = {step * dt_adv:6.3f}  "
                  f"relL2(mc1) = {err1[-1]:.4f}  relL2(mc4) = {err4[-1]:.4f}  "
                  f"[{time.time() - t0:5.1f} s]")

    return dict(snaps=snaps, times=np.array(times), err1=np.array(err1),
                err4=np.array(err4), neg=neg, lattice=mesh)


def make_figures(res: dict, worst_a: float) -> None:
    lat = res["lattice"]
    extent = [lat.x[0], lat.x[-1] + lat.dx, lat.p[0], lat.p[-1] + lat.dp]
    snap_times = sorted(res["snaps"]["mesh"].keys())
    rows = [("mesh QLE", "mesh"), ("single-rule MC", "mc1"),
            ("four-rule MC", "mc4")]
    vmax = max(np.abs(W).max() for W in res["snaps"]["mesh"].values())

    fig, axes = plt.subplots(4, len(snap_times), figsize=(3.4 * len(snap_times), 12),
                             sharex=True, sharey=True)
    for col, t in enumerate(snap_times):
        for r, (label, key) in enumerate(rows):
            ax = axes[r, col]
            im = ax.imshow(res["snaps"][key][t], origin="lower", extent=extent,
                           aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
            if r == 0:
                ax.set_title(f"t = {t:.2f}  ({t / T_PERIOD:.2f} T)")
            if col == 0:
                ax.set_ylabel(f"{label}\np")
        diff = res["snaps"]["mc4"][t] - res["snaps"]["mesh"][t]
        ax = axes[3, col]
        ax.imshow(diff, origin="lower", extent=extent, aspect="auto",
                  cmap="RdBu_r", vmin=-vmax / 3, vmax=vmax / 3)
        if col == 0:
            ax.set_ylabel("four-rule − mesh\np")
        ax.set_xlabel("x")
    fig.colorbar(im, ax=axes, shrink=0.85, label="W(x, p)")
    fig.suptitle(
        "Cosine well: mesh QLE vs single-rule MC vs four-rule MC "
        f"(nu = {NU}, identical splitting; generator identity "
        f"max err = {worst_a:.1e})")
    fname = "four_rule_equivalence_evolution.png"
    fig.savefig(output_path(fname), dpi=150, bbox_inches="tight")
    dpth = docs_path(fname)
    if dpth:
        fig.savefig(dpth, dpi=150, bbox_inches="tight")
    plt.close(fig)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
    ax1.plot(res["times"] / T_PERIOD, res["err1"], label="single-rule MC")
    ax1.plot(res["times"] / T_PERIOD, res["err4"], label="four-rule MC")
    ax1.set_xlabel("t / T_period")
    ax1.set_ylabel("relative L2 deviation from mesh QLE")
    ax1.set_title("Shot-noise floors (equal nu)")
    ax1.legend()
    ax1.grid(alpha=0.3)
    for key, label in (("mesh", "mesh QLE"), ("mc1", "single-rule MC"),
                       ("mc4", "four-rule MC")):
        ax2.plot(res["times"] / T_PERIOD, res["neg"][key], label=label)
    ax2.set_xlabel("t / T_period")
    ax2.set_ylabel(r"negativity $\int |\min(W,0)|\, dx\, dp$")
    ax2.set_title("Wigner negativity")
    ax2.legend()
    ax2.grid(alpha=0.3)
    fig.tight_layout()
    fname = "four_rule_equivalence_metrics.png"
    fig.savefig(output_path(fname), dpi=150, bbox_inches="tight")
    dpth = docs_path(fname)
    if dpth:
        fig.savefig(dpth, dpi=150, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    worst_a = part_a()
    res = evolve()
    make_figures(res, worst_a)
    print("\nFinal relative L2 deviation from mesh QLE:")
    print(f"  single-rule MC: {res['err1'][-1]:.4f}")
    print(f"  four-rule MC:   {res['err4'][-1]:.4f}")
    print("Figures written: four_rule_equivalence_evolution.png, "
          "four_rule_equivalence_metrics.png")
    return 0


if __name__ == "__main__":
    sys.exit(main())
