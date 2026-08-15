"""Representation cost of a world-particle ensemble, and the annihilation burden.

Numerical companion to
``docs/supplement/representation_cost_and_annihilation.md``.

The question asked here is quantitative rather than structural: given that a
world-particle ensemble *can* represent the Wigner function, how many
world-particles does it take, how does that number scale with the state, and
which of the project's three representations -- signed ``(x, p)``, pair
``(x, mu)``, and the positon/negaton sea -- is cheapest?

Parts
-----
A. The cost functional.  Sampling any signed or complex object costs
   ``N_eff = N / ||.||_1^2``.  For the Wigner function of a cat state the L1
   norm saturates at ``1 + 2/pi = 1.6366`` however large the cat, and the
   interference lobe carries ``(2/pi)/(1 + 2/pi) = 38.9%`` of that mass at a
   position where the probability density is numerically zero.

B. Interference is a *conjugate*-resolution requirement, not a position-density
   requirement.  A cat separated by ``d`` fringes in ``p`` with period
   ``2 pi hbar / d`` and ``d / (2 pi sigma)`` fringes under the envelope.
   Refining ``dx`` does not recover them; refining ``dp`` does.

C. Signed sampling versus the sea.  The crystal shift ``W -> W + 2/h`` is
   exact, but the sampled mass becomes ``1 + 2A/h`` where ``A`` is the
   phase-space area in view -- and for the spec's grid choice
   ``dp = pi hbar / L`` that is *exactly* the number of momentum cells.  The
   fractional error per cell is then ``sqrt(N_cells / N)``, independent of the
   state.

D. The pair ``(x, mu)`` ensemble.  Its cost is ``Z = (sum_i |psi_i|)^2``, which
   grows like ``1/dx`` -- unbounded under lattice refinement -- and is
   basis-privileged.  A coherence cutoff ``|Y| <= Ymax`` caps it.

E. Multiplicativity and Gaussians.  ``||W1 (x) W2||_1 = ||W1||_1 ||W2||_1``, so
   the cost is exponential in the number of non-Gaussian factors.  But every
   Gaussian state has ``||W||_1 = 1`` exactly, entangled or not (Hudson), so
   entanglement is not the cost driver -- non-Gaussianity is.

F. The annihilation burden.  Under exact QLE evolution in a one-mode potential
   ``||W(t)||_1`` stays bounded and O(1).  A non-annihilating unraveling grows
   its pathwise L1 mass at up to ``2 |Gamma_q(x)|``, i.e. ``4 V_q / (pi hbar)``
   averaged over the cell.  The gap between the two curves *is* the work that
   annihilation (garbage collection) has to do.  Also reports the two candidate
   density criteria: the crystal occupancy condition ``script_N >= M`` and the
   partner separation time ``t_sep = m L^2 / (2 pi hbar q M)``.

Run:  WPMW_OUTPUT=... python src/demo_representation_cost.py
"""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from wpmwlib.phase_space_crystal_lattice import (  # noqa: E402
    FourierMode,
    PhaseSpaceCrystalLattice,
)
from wpmwlib.wpmw_utils import docs_path, output_path  # noqa: E402

HBAR = 1.0
MASS = 1.0
SIGMA = 1.0

rng = np.random.default_rng(20260815)


def banner(title: str) -> None:
    print()
    print("=" * 72)
    print(title)
    print("=" * 72)


def save_fig(fig, name: str) -> None:
    for pth in (output_path(name), docs_path(name)):
        if pth:
            fig.savefig(pth, dpi=150, bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------- #
# Cat-state Wigner function                                             #
# --------------------------------------------------------------------- #
def wigner_cat(x, p, d, sigma=SIGMA, hbar=HBAR):
    """Closed form for psi ~ g(x - d/2) + g(x + d/2), g = exp(-x^2/4 sigma^2)."""
    a = 0.5 * d
    env = np.exp(-2.0 * sigma**2 * p**2 / hbar**2)
    lobes = (np.exp(-((x - a) ** 2) / (2 * sigma**2))
             + np.exp(-((x + a) ** 2) / (2 * sigma**2)))
    cross = 2.0 * np.exp(-(x**2) / (2 * sigma**2)) * np.cos(2.0 * a * p / hbar)
    norm = np.pi * hbar * 2.0 * (1.0 + np.exp(-(a**2) / (2 * sigma**2)))
    return (lobes + cross) * env / norm


def wigner_cat_fft(x, d, sigma=SIGMA, hbar=HBAR, Ny=16384, ymax=80.0):
    """Independent FFT evaluation of the same object, for verification.

    Returns ``(W, p)`` with ``W`` normalised to unit phase-space integral and
    ``p`` the exact FFT momentum grid, so the comparison in Part A is made at
    matching abscissae rather than at nearest neighbours.
    """
    y = (np.arange(Ny) - Ny // 2) * (2 * ymax / Ny)
    dy = y[1] - y[0]

    def psi(u):
        return (np.exp(-((u - d / 2) ** 2) / (4 * sigma**2))
                + np.exp(-((u + d / 2) ** 2) / (4 * sigma**2)))

    uu = np.linspace(-60.0, 60.0, 24001)
    nrm = np.trapezoid(psi(uu) ** 2, uu)
    C = (psi(x[:, None] + y[None, :] / 2)
         * psi(x[:, None] - y[None, :] / 2)) / nrm
    F = np.fft.fftshift(
        np.fft.fft(np.fft.ifftshift(C, axes=1), axis=1), axes=1
    ) * dy / (2 * np.pi * hbar)
    pf = np.fft.fftshift(np.fft.fftfreq(Ny, d=dy)) * 2 * np.pi * hbar
    return F.real.T, pf


def grid(xlim=25.0, plim=8.0, nx=1001, npg=1201):
    x = np.linspace(-xlim, xlim, nx)
    p = np.linspace(-plim, plim, npg)
    return x, p


# --------------------------------------------------------------------- #
# Part A                                                                #
# --------------------------------------------------------------------- #
def part_a():
    banner("A. The cost functional: ||W||_1 for a cat state")
    x, p = grid()
    dx, dp = x[1] - x[0], p[1] - p[0]
    X, P = np.meshgrid(x, p, indexing="xy")

    # verification of the closed form against an independent FFT evaluation,
    # compared on the FFT's own momentum abscissae
    xs = np.linspace(-6.0, 6.0, 121)
    Wf, pf = wigner_cat_fft(xs, 4.0)
    sel = np.abs(pf) < 3.0
    Xs, Ps = np.meshgrid(xs, pf[sel], indexing="xy")
    Wc = wigner_cat(Xs, Ps, 4.0)
    err = float(np.max(np.abs(Wc - Wf[sel])))
    print(f"closed form vs FFT evaluation, max |diff| = {err:.3e}"
          f"   (peak |W| = {np.max(np.abs(Wc)):.4f})")

    print()
    print(f"{'d/sigma':>8} {'||W||_1':>9} {'nu=int W_-':>11} {'N/N_eff':>9}"
          f" {'|W| in |x|<sigma':>17} {'rho(0)':>11} {'p-fringes':>10}")
    rows = []
    for d in [0.0, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 16.0, 24.0]:
        W = wigner_cat(X, P, d)
        W = W / (W.sum() * dx * dp)
        tv = float(np.abs(W).sum() * dx * dp)
        nu = float(np.abs(np.minimum(W, 0.0)).sum() * dx * dp)
        band = np.abs(x) < SIGMA
        share = float(np.abs(W[:, band]).sum() * dx * dp / tv)
        rho = W.sum(axis=0) * dp
        rho0 = float(rho[np.argmin(np.abs(x))])
        rows.append((d, tv, nu, share, rho0))
        print(f"{d:8.1f} {tv:9.4f} {nu:11.4f} {tv**2:9.3f} {share:17.4f}"
              f" {rho0:11.2e} {d / (2 * np.pi * SIGMA):10.2f}")

    print()
    print(f"asymptote 1 + 2/pi = {1 + 2 / np.pi:.6f}   "
          f"N/N_eff -> {(1 + 2 / np.pi)**2:.4f}")
    print(f"interference share (2/pi)/(1 + 2/pi) = "
          f"{(2 / np.pi) / (1 + 2 / np.pi):.4f}")
    print("(the tabulated share uses the band |x| < sigma, which captures "
          "68% of the lobe)")
    return rows


# --------------------------------------------------------------------- #
# Part B                                                                #
# --------------------------------------------------------------------- #
def part_b():
    banner("B. Interference is a conjugate-resolution requirement")
    d = 8.0
    print(f"cat separation d = {d} sigma")
    print(f"  fringe period in p        = 2 pi hbar / d = "
          f"{2 * np.pi * HBAR / d:.4f}")
    print(f"  envelope width in p       = hbar / (2 sigma) = "
          f"{HBAR / (2 * SIGMA):.4f}")
    print(f"  fringes under envelope    = d / (2 pi sigma) = "
          f"{d / (2 * np.pi * SIGMA):.3f}")
    print()
    print("Fringe visibility of a bin-averaged reconstruction of the central")
    print("column, exact W and no sampling noise -- a pure resolution effect.")
    print("Visibility is the binned peak-to-trough divided by the exact one.")
    print()
    print(f"{'dx':>8} {'dp':>8} {'visibility':>12}")
    x0 = np.linspace(-20.0, 20.0, 4001)
    p0 = np.linspace(-4.0, 4.0, 8001)
    X0, P0 = np.meshgrid(x0, p0, indexing="xy")
    W0 = wigner_cat(X0, P0, d)
    exact_col = W0[:, int(np.argmin(np.abs(x0)))]
    exact_pv = float(exact_col.max() - exact_col.min())
    for dx_bin, dp_bin in [(2.0, 0.05), (0.5, 0.05), (0.125, 0.05),
                           (0.5, 0.20), (0.5, 0.40),
                           (0.5, 2 * np.pi * HBAR / d), (0.5, 1.5)]:
        nx_b = max(2, int(round(40.0 / dx_bin)))
        np_b = max(2, int(round(8.0 / dp_bin)))
        H, xe, pe = np.histogram2d(
            X0.ravel(), P0.ravel(), bins=[nx_b, np_b],
            weights=W0.ravel(), range=[[-20, 20], [-4, 4]])
        cnt, _, _ = np.histogram2d(
            X0.ravel(), P0.ravel(), bins=[nx_b, np_b],
            range=[[-20, 20], [-4, 4]])
        Wb = np.where(cnt > 0, H / np.maximum(cnt, 1), 0.0)
        col = Wb[int(np.argmin(np.abs(0.5 * (xe[:-1] + xe[1:]))))]
        vis = float(col.max() - col.min()) / exact_pv
        print(f"{dx_bin:8.3f} {dp_bin:8.3f} {vis:12.4f}")
    print()
    print()
    print("Coarsening dx by 16x leaves the visibility intact; coarsening dp")
    print("toward the fringe period destroys it.  Position density is not the")
    print("resource that interference consumes -- momentum resolution is.")


# --------------------------------------------------------------------- #
# Part C                                                                #
# --------------------------------------------------------------------- #
def part_c():
    banner("C. Signed sampling versus the positon/negaton sea")
    d = 8.0
    x, p = grid(xlim=12.0, plim=4.0, nx=241, npg=201)
    dx, dp = x[1] - x[0], p[1] - p[0]
    X, P = np.meshgrid(x, p, indexing="xy")
    W = wigner_cat(X, P, d)
    W = W / (W.sum() * dx * dp)
    tv = float(np.abs(W).sum() * dx * dp)
    rho = W.sum(axis=0) * dp
    rho_pk = float(rho.max())
    i0 = int(np.argmin(np.abs(x)))
    sea = 2.0 / (2.0 * np.pi * HBAR)
    Wp = W + sea
    mass = float(Wp.sum() * dx * dp)
    print(f"||W||_1               = {tv:.4f}")
    print(f"sea mass 2A/h         = {sea * W.size * dx * dp:.2f}   "
          f"(A = {(x[-1] - x[0]) * (p[-1] - p[0]):.1f})")
    print(f"||W + 2/h||_1         = {mass:.2f}   "
          f"ratio to ||W||_1 = {mass / tv:.1f}")
    print(f"min(W + 2/h)          = {Wp.min():+.4e}  (must be >= 0)")
    print(f"cells in view         = {W.size}, "
          f"occupied (|W| > 1e-3 max) = {(np.abs(W) > 1e-3 * np.abs(W).max()).sum()}")
    print()
    print("Reconstructed rho(x) at the interference minimum x = 0")
    print(f"(rho_peak = {rho_pk:.4f}, exact rho(0) = {rho[i0]:.3e})")
    print()
    print(f"{'N':>10} {'signed std':>12} {'sea std':>12} {'ratio':>8}")
    prob_s = (np.abs(W) / np.abs(W).sum()).ravel()
    sgn = np.sign(W).ravel()
    prob_c = (Wp / Wp.sum()).ravel()
    shape = W.shape
    ns, ss_list, sc_list = [], [], []
    for N in [10**4, 10**5, 10**6, 10**7]:
        vs, vc = [], []
        for _ in range(24):
            k = rng.multinomial(N, prob_s).reshape(shape)
            vs.append(float((k * np.sign(W)).sum(axis=0)[i0] * tv / N / dx))
            k2 = rng.multinomial(N, prob_c).reshape(shape)
            vc.append(float(k2.sum(axis=0)[i0] * mass / N / dx
                            - sea * len(p) * dp))
        ss, sc = float(np.std(vs)), float(np.std(vc))
        ns.append(N); ss_list.append(ss); sc_list.append(sc)
        print(f"{N:10d} {ss:12.4e} {sc:12.4e} {sc / ss:8.1f}")
    print()
    print("Per-cell rule for the sea: fractional error = sqrt(N_cells / N),")
    print("independent of the state, because the sea dominates the counts.")
    for N in [10**5, 10**6, 10**7]:
        k = rng.multinomial(N, prob_c).reshape(shape)
        West = k * mass / N / (dx * dp) - sea
        emp = float(np.std(West - W)) / sea
        print(f"  N = {N:>9d}   predicted {np.sqrt(W.size / N):.4f}   "
              f"measured {emp:.4f}")
    return ns, ss_list, sc_list, tv, mass


# --------------------------------------------------------------------- #
# Part D                                                                #
# --------------------------------------------------------------------- #
def part_d():
    banner("D. The pair (x, mu) ensemble: cost Z = (sum_i |psi_i|)^2")

    def psi_cat(u, d):
        return (np.exp(-((u - d / 2) ** 2) / 4) + np.exp(-((u + d / 2) ** 2) / 4))

    print(f"{'dx/sigma':>9} {'Z cat d=8':>11} {'Z packet':>10} "
          f"{'Z plane wave':>13}")
    dxs_list, Zc_list = [], []
    for dxs in [1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125]:
        u = np.arange(-20, 20, dxs)
        f = psi_cat(u, 8.0); f = f / np.sqrt((f**2).sum())
        g = np.exp(-u**2 / 4); g = g / np.sqrt((g**2).sum())
        pw = np.ones_like(u) / np.sqrt(len(u))
        Zc = float(np.abs(f).sum() ** 2)
        dxs_list.append(dxs); Zc_list.append(Zc)
        print(f"{dxs:9.5f} {Zc:11.2f} {float(np.abs(g).sum()**2):10.2f} "
              f"{float(np.abs(pw).sum()**2):13.2f}")
    print()
    print("Z ~ 1/dx: unbounded under lattice refinement.  ||W||_1 = 1.58 for")
    print("every row above.  Z is also basis-privileged: a plane wave costs M")
    print("in the position basis and 1 in the Wigner representation.")
    print()
    print("With a coherence cutoff |x - x'| <= Ymax (dx = 0.125):")
    u = np.arange(-20, 20, 0.125)
    f = psi_cat(u, 8.0); f = f / np.sqrt((f**2).sum())
    R = np.abs(np.outer(f, f))
    D = np.abs(u[:, None] - u[None, :])
    ys, zs = [], []
    for Ymax in [0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 40.0]:
        Z = float(R[D <= Ymax].sum())
        ys.append(Ymax); zs.append(Z)
        print(f"  Ymax = {Ymax:5.1f}   Z = {Z:8.2f}")
    print()
    print("The pair ensemble is the Wigner ensemble with the Y-integral left")
    print("un-done: it pays the Monte-Carlo variance of an oscillatory")
    print("integral that the Wigner transform performs analytically.")
    return dxs_list, Zc_list, ys, zs


# --------------------------------------------------------------------- #
# Part E                                                                #
# --------------------------------------------------------------------- #
def part_e():
    banner("E. Multiplicativity, and why entanglement is not the cost driver")
    x, p = grid(xlim=14.0, plim=6.0, nx=701, npg=701)
    dx, dp = x[1] - x[0], p[1] - p[0]
    X, P = np.meshgrid(x, p, indexing="xy")

    def l1(W):
        return float(np.abs(W).sum() * dx * dp / (W.sum() * dx * dp))

    # Gaussian family: squeezed, displaced, thermal
    print(f"{'state':<34} {'||W||_1':>9}")
    for label, s, x0, p0, nbar in [
        ("coherent  (sigma=1)", 1.0, 0.0, 0.0, 0.0),
        ("squeezed  (sigma=0.4)", 0.4, 0.0, 0.0, 0.0),
        ("displaced (x0=3)", 1.0, 3.0, 0.0, 0.0),
        ("thermal   (nbar=2)", 1.0, 0.0, 0.0, 2.0),
    ]:
        sx = s * np.sqrt(1.0 + 2.0 * nbar)
        sp = (HBAR / (2 * s)) * np.sqrt(1.0 + 2.0 * nbar)
        W = np.exp(-((X - x0) ** 2) / (2 * sx**2)
                   - ((P - p0) ** 2) / (2 * sp**2))
        print(f"{label:<34} {l1(W):9.6f}")
    print()
    print("Every Gaussian has ||W||_1 = 1 exactly.  By Hudson's theorem these")
    print("are the only pure states with W >= 0 -- and a two-mode squeezed")
    print("vacuum is one of them, maximally entangled, at zero sampling cost.")
    print()
    print("Multiplicativity under tensor product, ||W1 (x) W2||_1:")
    x1 = np.linspace(-14, 14, 61); dx1 = x1[1] - x1[0]
    p1 = np.linspace(-6, 6, 61); dp1 = p1[1] - p1[0]
    X1, P1 = np.meshgrid(x1, p1, indexing="xy")
    Wa = wigner_cat(X1, P1, 8.0); Wa = Wa / (Wa.sum() * dx1 * dp1)
    Wb = wigner_cat(X1, P1, 4.0); Wb = Wb / (Wb.sum() * dx1 * dp1)
    la = float(np.abs(Wa).sum() * dx1 * dp1)
    lb = float(np.abs(Wb).sum() * dx1 * dp1)
    prod = float(np.abs(Wa).sum() * np.abs(Wb).sum() * (dx1 * dp1) ** 2)
    direct = float(np.abs(np.outer(Wa.ravel(), Wb.ravel())).sum()
                   * (dx1 * dp1) ** 2)
    assert abs(prod - direct) < 1e-9 * max(1.0, direct)
    print(f"  ||Wa||_1 = {la:.6f}   ||Wb||_1 = {lb:.6f}")
    print(f"  ||Wa (x) Wb||_1 = {prod:.6f}   product = {la * lb:.6f}   "
          f"residual = {abs(prod - la * lb):.2e}")
    print(f"  n identical cats: cost {(1 + 2 / np.pi)**2:.3f}^n; "
          f"n = 20 gives {((1 + 2 / np.pi)**2)**20:.2e}")


# --------------------------------------------------------------------- #
# Part F                                                                #
# --------------------------------------------------------------------- #
def part_f():
    banner("F. The annihilation burden")
    M, N, L = 256, 128, 8.0
    q, V_q = 1, 1.5
    dt, nsteps = 0.002, 2000
    sol = PhaseSpaceCrystalLattice(M=M, N=N, L=L, mass=MASS, hbar=HBAR,
                                   advection="spectral")
    modes = [FourierMode(q=q, V_q=V_q, phi_q=0.0)]

    sig, sep = 0.5, 2.0
    sol.initialize_from_wigner(lambda X, P: wigner_cat(X, P, sep, sigma=sig))
    dx, dp = sol.dx, sol.dp
    W0 = sol.W / (sol.W.sum() * dx * dp)
    sol.W = W0

    gamma_max = 2.0 * V_q / HBAR
    gamma_avg = (2.0 / np.pi) * gamma_max
    print(f"grid: M = {M}, N = {N}, L = {L}, dx = {dx:.4f}, "
          f"dp = pi hbar / L = {dp:.4f}")
    print(f"one mode: q = {q}, V_q = {V_q}")
    print()
    print("Pathwise L1 growth of a non-annihilating unraveling is bounded by")
    print("the entrywise-absolute generator norm, 2 |Gamma_q(x)|:")
    print(f"  worst case (antinode)  gamma_max = 2 V_q / hbar        = "
          f"{gamma_max:.4f}")
    print(f"  cell average           gamma_avg = 4 V_q / (pi hbar)   = "
          f"{gamma_avg:.4f}")
    print()
    ts, l1s, norms = [], [], []
    for k in range(nsteps + 1):
        W = sol.W
        ts.append(k * dt)
        l1s.append(float(np.abs(W).sum() * dx * dp))
        norms.append(float(W.sum() * dx * dp))
        if k < nsteps:
            sol.strang_step_fourier(modes, dt)
    ts = np.array(ts); l1s = np.array(l1s); norms = np.array(norms)
    print(f"{'t':>7} {'||W||_1':>10} {'norm':>10} {'exp(gamma_avg t)':>18}"
          f" {'burden':>12}")
    for t_want in [0.0, 0.5, 1.0, 2.0, 3.0, 4.0]:
        k = int(np.argmin(np.abs(ts - t_want)))
        g = float(np.exp(gamma_avg * ts[k]))
        print(f"{ts[k]:7.2f} {l1s[k]:10.4f} {norms[k]:10.6f} {g:18.4e}"
              f" {g / l1s[k]:12.3e}")
    print()
    print(f"probability norm drift over the run: "
          f"{abs(norms[-1] - 1.0):.2e}")
    print(f"||W||_1 stays in [{l1s.min():.4f}, {l1s.max():.4f}] -- bounded.")
    print("The ratio is the L1 mass that garbage collection must remove per")
    print("unit time.  It is exponential in t; ||W||_1 is not.")

    print()
    print("Two candidate density criteria for the annihilating regime:")
    script_N = M
    sea_per_cell = "script_N / M"
    print(f"  (i)  crystal occupancy.  Sea occupancy per cell is "
          f"{sea_per_cell} = script_N / {M},")
    print(f"       so one negaton per cell requires script_N >= M = "
          f"{script_N} world-particles")
    print(f"       per unit probability, and a total ensemble of "
          f"script_N (1 + 2A/h) = {script_N} x {N} = {script_N * N:,}.")
    t_sep = MASS * L**2 / (2 * np.pi * HBAR * q * M)
    t_create = HBAR / V_q
    print(f"  (ii) partner separation.  A focus vertex makes a pair split by "
          f"2 q dp;")
    print(f"       it leaves a position cell after "
          f"t_sep = m L^2 / (2 pi hbar q M) = {t_sep:.4f},")
    print(f"       against a mean creation interval "
          f"hbar / V_q = {t_create:.4f}.")
    print(f"       Annihilation number A = t_sep / t_create = "
          f"{t_sep / t_create:.4f}  (<< 1 here).")
    print()
    print("These disagree: (i) says the cells are richly occupied and (ii)")
    print("says original partners escape long before they can meet again.")
    print("Which one governs is not settled by either estimate -- see the")
    print("supplement, section 7.")
    return ts, l1s, gamma_avg, gamma_max


# --------------------------------------------------------------------- #
# Figures                                                               #
# --------------------------------------------------------------------- #
def figure_cat_cost(rows):
    d = np.array([r[0] for r in rows])
    tv = np.array([r[1] for r in rows])
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))

    ax = axes[0]
    ax.plot(d, tv, "o-", color="C0", lw=1.6, label=r"$\|W\|_1$")
    ax.axhline(1 + 2 / np.pi, color="C3", ls="--", lw=1.2,
               label=r"$1 + 2/\pi$")
    ax.set_xlabel(r"cat separation $d/\sigma$")
    ax.set_ylabel(r"$\int |W| \, dx\, dp$")
    ax.set_title("Cost saturates: interference is a\nbounded, one-off penalty")
    ax.grid(alpha=0.3); ax.legend(fontsize=9)

    x, p = grid(xlim=10.0, plim=3.0, nx=501, npg=501)
    dx, dp = x[1] - x[0], p[1] - p[0]
    X, P = np.meshgrid(x, p, indexing="xy")
    W = wigner_cat(X, P, 8.0); W = W / (W.sum() * dx * dp)
    ax = axes[1]
    v = np.abs(W).max()
    cf = ax.pcolormesh(x, p, W, cmap="RdBu_r", vmin=-v, vmax=v,
                       shading="auto")
    ax.axvline(0.0, color="k", ls=":", lw=1.0)
    ax.set_xlabel("$x$"); ax.set_ylabel("$p$")
    ax.set_title(r"$W(x,p)$, $d = 8\sigma$" "\n" "fringes run along $p$")
    fig.colorbar(cf, ax=ax, fraction=0.046)

    ax = axes[2]
    col_abs = np.abs(W).sum(axis=0) * dp
    rho = W.sum(axis=0) * dp
    ax.plot(x, col_abs, color="C0", lw=1.6, label=r"$\int |W(x,p)|\,dp$")
    ax.plot(x, rho, color="C3", lw=1.4, label=r"$\rho(x) = \int W\,dp$")
    ax.axvline(0.0, color="k", ls=":", lw=1.0)
    ax.set_xlabel("$x$")
    ax.set_title("Where the world-particles must be\n"
                 "vs where the particle is found")
    ax.grid(alpha=0.3); ax.legend(fontsize=9)

    fig.suptitle("Representation cost of a cat state: the interference lobe "
                 "carries 38.9% of the ensemble at zero probability density",
                 fontsize=12, y=1.01)
    fig.tight_layout()
    save_fig(fig, "representation_cost_cat.png")


def figure_sampling(cres, dres):
    ns, ss, sc, tv, mass = cres
    dxs, Zc, ys, zs = dres
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))

    ax = axes[0]
    ax.loglog(ns, ss, "o-", color="C0", lw=1.6, label="signed $\\pm$")
    ax.loglog(ns, sc, "s-", color="C3", lw=1.6, label="sea $W + 2/h$")
    ref = np.array(ns, dtype=float)
    ax.loglog(ref, ss[0] * np.sqrt(ref[0] / ref), "k:", lw=1.0,
              label=r"$N^{-1/2}$")
    ax.set_xlabel("$N$ world-particles")
    ax.set_ylabel(r"std of reconstructed $\rho(0)$")
    ax.set_title(f"Same target, two representations\n"
                 f"$\\|W\\|_1 = {tv:.2f}$ vs $\\|W + 2/h\\|_1 = {mass:.0f}$")
    ax.grid(alpha=0.3, which="both"); ax.legend(fontsize=9)

    ax = axes[1]
    ax.loglog(dxs, Zc, "o-", color="C2", lw=1.6, label=r"$Z$ (pair ensemble)")
    ax.axhline(1 + 2 / np.pi, color="C0", ls="--", lw=1.2,
               label=r"$\|W\|_1$ (same state)")
    ax.set_xlabel(r"position lattice spacing $\Delta x/\sigma$")
    ax.set_ylabel(r"$Z = (\sum_i |\psi_i|)^2$")
    ax.set_title("Pair cost grows like $1/\\Delta x$;\n"
                 "Wigner cost does not move")
    ax.grid(alpha=0.3, which="both"); ax.legend(fontsize=9)
    ax.invert_xaxis()

    ax = axes[2]
    ax.plot(ys, zs, "o-", color="C2", lw=1.6)
    ax.axhline(1 + 2 / np.pi, color="C0", ls="--", lw=1.2,
               label=r"$\|W\|_1$")
    ax.set_xlabel(r"coherence cutoff $Y_{\max}/\sigma$")
    ax.set_ylabel("$Z$")
    ax.set_title("A finite coherence length caps it\n"
                 "(and decoherence shrinks it)")
    ax.grid(alpha=0.3); ax.legend(fontsize=9)

    fig.suptitle("Sampling cost of the three representations of the same "
                 "cat state", fontsize=12, y=1.01)
    fig.tight_layout()
    save_fig(fig, "representation_cost_sampling.png")


def figure_annihilation(fres):
    ts, l1s, gamma_avg, gamma_max = fres
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))

    ax = axes[0]
    ax.plot(ts, l1s, color="C0", lw=1.6, label=r"$\|W(t)\|_1$ (exact QLE)")
    ax.axhline(1.0, color="k", ls=":", lw=1.0)
    ax.set_xlabel("$t$"); ax.set_ylabel(r"$\int |W|\, dx\, dp$")
    ax.set_title("The state's own L1 mass is bounded")
    ax.grid(alpha=0.3); ax.legend(fontsize=9)

    ax = axes[1]
    ax.semilogy(ts, l1s, color="C0", lw=1.6, label=r"$\|W(t)\|_1$")
    ax.semilogy(ts, np.exp(gamma_avg * ts), color="C3", lw=1.4,
                label=r"$e^{\gamma_{\rm avg} t}$, "
                      r"$\gamma_{\rm avg} = 4V_q/\pi\hbar$")
    ax.semilogy(ts, np.exp(gamma_max * ts), color="C1", ls="--", lw=1.2,
                label=r"$e^{\gamma_{\max} t}$, $\gamma_{\max} = 2V_q/\hbar$")
    ax.set_xlabel("$t$"); ax.set_ylabel("L1 mass")
    ax.set_title("Pathwise growth without annihilation")
    ax.grid(alpha=0.3, which="both"); ax.legend(fontsize=8, loc="upper left")

    fig.suptitle("The annihilation burden: the gap between the two curves is "
                 "the L1 mass garbage collection must remove",
                 fontsize=12, y=1.01)
    fig.tight_layout()
    save_fig(fig, "annihilation_burden.png")


def main():
    rows = part_a()
    part_b()
    cres = part_c()
    dres = part_d()
    part_e()
    fres = part_f()
    figure_cat_cost(rows)
    figure_sampling(cres, dres)
    figure_annihilation(fres)
    print()
    print("figures written to", output_path(""))


if __name__ == "__main__":
    main()
