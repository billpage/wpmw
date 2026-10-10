"""How the sea reads the potential -- companion to docs/supplement/chord_reading.md.

Produces every number and the figure of the tutorial.

  Part A  what the sum over chords estimates.
          A1 The potential term of the von Neumann equation, Wigner
             transformed, is a convolution in momentum,
                 dW/dt|_V (x, p) = int J(x, p - p') W(x, p') dp',
                 J(x, xi) = -(B/hbar) int U(x, y) sin(2 xi y / hbar) dy,
             U = V(x+y) - V(x-y), B = 1/(pi hbar).  Checked on a cat state
             in the Eckart barrier, on the discrete grid where it is exact.
          A2 The jump rates are J sampled at the channel transfers, tapered
             and compensated: K_q(x) = J_res(x, xi_q) dp.  The continuum
             formula against the code's kernel at grid points.
          A3 The reading is blind to the state: J contains U only.  The
             same x gives the same K_q for a packet, a cat, or nothing.
  Part B  the chord estimator on a static locked sea (the reading at
          re-lock, mu = 0).  Pairs on one row are a Poisson configuration
          of density nu B dp per unit length; a reader at x0 sums the chords
          whose midpoints lie within dX/2 of it.  Mean, spread, and chord
          count against nu; nu = 1 is the physical density.
  Figure  chord_reading_estimator.png.

Units hbar = m = 1; Eckart V0 = a = 1; dp = 0.25, y_max = 2 pi.

Run as::

    WPMW_OUTPUT=... PYTHONPATH=src python3 -u src/demo_chord_reading.py
"""
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from demo_emission_and_absorption import Ledger, B, HBAR
from wpmwlib.wpmw_utils import output_path, docs_path

DP = 0.25
YMAX = np.pi * HBAR / (2 * DP)            # = 2 pi: B dp 2 y_max = 1
NQ = 31
XI = np.arange(1, NQ + 1) * DP
X0 = -0.625                               # the Eckart flank site of Q5, Q7, Q11


def V(x):
    return 1.0 / np.cosh(x) ** 2


def banner(t):
    print("\n" + "=" * 78 + f"\n{t}\n" + "=" * 78, flush=True)


# ---------------------------------------------------------------- Part A
def kernel_continuum(x, ny=801):
    """Theorem L1 by quadrature: taper w, compensation by the first moment
    over channels 1..31 (as the code's dv_eff).  Returns K_q, V'_eff."""
    y = np.linspace(0, YMAX, ny)
    wq = np.full(ny, y[1] - y[0])
    wq[[0, -1]] *= 0.5
    w = np.cos(np.pi * y / (2 * YMAX)) ** 2
    s = np.sin(2 * XI[:, None] * y[None, :] / HBAR)
    u = V(x + y) - V(x - y)
    kr = -(2 * B * DP / HBAR) * (s * (u * w * wq)).sum(1)
    kc = -(2 * B * DP / HBAR) * (s * (2 * y * w * wq)).sum(1)
    dv = (XI @ kr) / (XI @ kc)
    return kr - dv * kc, dv


def part_a():
    banner("A1. The potential term of the von Neumann equation, Wigner transformed")
    n, L = 256, 32.0
    dx = L / n
    x = -L / 2 + dx * np.arange(n)
    sig, x0, p0 = 1.0, 3.0, 1.0

    def g(c):
        return np.exp(-(x - c) ** 2 / (4 * sig ** 2) + 1j * p0 * x / HBAR)
    psi = g(-x0) + g(x0)
    psi /= np.sqrt((np.abs(psi) ** 2).sum() * dx)
    # off-diagonal coordinate y_k = k dx, so x +- y stay on the grid
    k = np.arange(-n // 2, n // 2)
    yk = k * dx
    ip = (np.arange(n)[:, None] + k[None, :]) % n
    im = (np.arange(n)[:, None] - k[None, :]) % n
    rho = psi[ip] * np.conj(psi[im])                       # rho(x + y, x - y)
    U = V(x[ip]) - V(x[im])
    # Wigner transform on the conjugate grid p_m: spacing pi hbar / (n dx)
    dpg = np.pi * HBAR / (n * dx)
    pm = dpg * np.arange(-n // 2, n // 2)
    F = np.exp(-2j * pm[None, :] * yk[:, None] / HBAR)      # (k, m)
    W = (rho @ F).real * dx / (np.pi * HBAR)
    dWdt = ((U * rho / (1j * HBAR)) @ F).real * dx / (np.pi * HBAR)
    # J(x, xi) on the difference grid, then the convolution sum over p'
    J = -(B / HBAR) * (U[:, :, None] * np.sin(2 * (dpg * np.arange(-n // 2, n // 2))[None, None, :]
                                              * yk[None, :, None] / HBAR)).sum(1) * dx
    # sum over m' of J(p_m - p_m') W(p_m') dp: J is periodic in xi with period
    # n dpg (sin(2 xi y_k / hbar) with y_k = k dx), so the difference index is cyclic
    conv = np.zeros_like(W)
    mp = np.arange(n)
    for m in range(n):
        conv[:, m] = (J[:, (m - mp + n // 2) % n] * W).sum(1) * dpg
    err = np.linalg.norm(conv - dWdt) / np.linalg.norm(dWdt)
    print(f"  cat state, x0 = +-{x0}, p0 = {p0}, sigma = {sig}; grid {n} x {n}, dx = {dx}")
    print(f"  norm of W {W.sum() * dx * dpg:.12f};  min W {W.min():+.4f} (the cat's fringes)")
    print(f"  || int J W dp' - (1/i hbar) Wigner[ [V, rho] ] || / || ... || = {err:.2e}")
    print("  -> the kernel is the Wigner transform, over the off-diagonal coordinate y, of the\n"
          "     phase rate (V(x+y) - V(x-y)) / hbar that the potential applies to rho(x+y, x-y).")

    banner("A2. K_q is J at the channel transfers, tapered and compensated (Theorem L1)")
    run = Ledger(v0=1.0, a=1.0, n_r=128, r_half=20.0, n_p=64, dp=DP)
    print(f"  {'x':>8} {'K_1 here':>10} {'code':>10} {'K_2 here':>10} {'code':>10}"
          f" {'Veff here':>10} {'code':>10}")
    for i in (62, 67, 70, 75):
        xr = run.r[i]
        kq, dv = kernel_continuum(xr)
        print(f"  {xr:8.4f} {kq[0]:10.6f} {run.k[i, 1]:10.6f} {kq[1]:10.6f} {run.k[i, 2]:10.6f}"
              f" {dv:10.6f} {run.dv_eff[i]:10.6f}")
    print("  (the code samples y at 64 points of the momentum grid's conjugate; the quadrature"
          " here uses 801)")

    banner("A3. The reading does not see the state")
    print("  J and K_q contain U(x, y) and nothing of rho: a packet, a cat, or an empty")
    print("  box give the same K_q at the same x.  The state enters an event only through")
    print("  the parent: the event rate per unit area is |K_q| times the body density, and")
    print("  the generator sum_q K_q [W(p - xi_q) - W(p + xi_q)] is A1's convolution,")
    print("  discretised and with the classical part moved to streaming.")


# ---------------------------------------------------------------- Part B
def reading(xs, xr, dX, n, dv):
    """Chord sum at re-lock (mu = 0): Theorem L1 as a sum over chords of one row
    whose midpoints lie within dX/2 of the reader; Proposition L2(c) projection."""
    xs = np.sort(xs)
    I, J = np.triu_indices(len(xs), 1)
    d = xs[J] - xs[I]
    xm = 0.5 * (xs[I] + xs[J])
    sel = (d <= 2 * YMAX) & (np.abs(xm - xr) < dX / 2)
    if not sel.any():
        return np.zeros(NQ), 0
    I, J, d = I[sel], J[sel], d[sel]
    w = np.cos(np.pi * d / (4 * YMAX)) ** 2
    rate = V(xs[J]) - V(xs[I]) - d * dv                      # hbar dmu/dt = U_res
    kh = -(B * DP / (HBAR * n ** 2 * dX)) * (w * rate) @ np.sin(XI[None, :] * d[:, None] / HBAR)
    return kh - XI * (XI @ kh) / (XI @ XI), int(sel.sum())


def part_b():
    banner("B. The chord estimator on a static locked sea, against the sea's multiplicity nu")
    kq, dv = kernel_continuum(X0)
    rng = np.random.default_rng(5)
    dX = 1.0
    rows = []
    print(f"  reader at x = {X0}, aperture dX = {dX}; K_1..4 = {np.round(kq[:4], 4)}")
    print(f"  {'nu':>4} {'readings':>9} {'chords':>8} {'P(none)':>8}"
          f" {'mean K1/K1':>11} {'+-SE':>6} {'mean K2/K2':>11} {'sd(K1)/K1':>10}")
    for nu, M in ((1, 400000), (2, 200000), (4, 80000), (8, 30000), (16, 12000),
                  (32, 4000), (64, 1500)):
        n = nu * B * DP
        lw = 2 * YMAX + dX
        acc = np.zeros((M, NQ))
        nch = np.zeros(M)
        for m in range(M):
            cnt = rng.poisson(n * lw)
            xs = rng.uniform(X0 - lw / 2, X0 + lw / 2, cnt)
            acc[m], nch[m] = reading(xs, X0, dX, n, dv)
        mean = acc.mean(0)
        se = acc.std(0) / np.sqrt(M)
        sd1 = acc[:, 0].std() / abs(kq[0])
        print(f"  {nu:4d} {M:9d} {nch.mean():8.2f} {np.mean(nch == 0):8.3f}"
              f" {mean[0] / kq[0]:11.4f} {se[0] / abs(kq[0]):6.4f} {mean[1] / kq[1]:11.4f}"
              f" {sd1:10.2f}", flush=True)
        rows.append((nu, nch.mean(), np.mean(nch == 0), mean, se, sd1))
    print("  -> unbiased at every density (the reading is a U-statistic of the sea's"
          " configuration);\n     at nu = 1 a reader has a chord in its aperture about one"
          " time in twelve.")
    return kq, rows


def figure(kq, rows):
    blue, orange = "#2a78d6", "#eb6834"
    ink, muted, grid = "#0b0b0b", "#52514e", "#e4e2dc"
    fig, ax = plt.subplots(1, 2, figsize=(12.0, 4.2))
    q = np.arange(1, 9)
    pick = {r[0]: r for r in rows}
    wbar = 0.38
    ax[0].bar(q - wbar / 2, kq[:8], wbar, color=blue, label="kernel $K_q$ (Theorem L1)")
    r16 = pick[16]
    ax[0].bar(q + wbar / 2, r16[3][:8], wbar, color=orange,
              label=r"mean chord reading, $\nu = 16$")
    ax[0].errorbar(q + wbar / 2, r16[3][:8], yerr=2 * r16[4][:8], fmt="none", ecolor=ink,
                   elinewidth=1, capsize=2)
    ax[0].axhline(0, color=muted, lw=0.8)
    ax[0].set_xlabel("channel $q$ (transfer $\\xi_q = q\\,dp$)", color=ink)
    ax[0].set_ylabel("rate per parent", color=ink)
    ax[0].set_title(f"The chords estimate the kernel (Eckart flank, x = {X0})", color=ink,
                    fontsize=11)
    ax[0].legend(frameon=False, fontsize=9)
    nus = np.array([r[0] for r in rows])
    sd = np.array([r[5] for r in rows])
    ax[1].loglog(nus, sd, "o-", color=blue, lw=2, ms=6, label="spread of one reading")
    ax[1].loglog(nus, sd[0] / nus, ":", color=muted, lw=1.2, label=r"$\propto 1/\nu$")
    ax[1].loglog(nus, sd[0] / np.sqrt(nus), "--", color=muted, lw=1.0,
                 label=r"$\propto 1/\sqrt{\nu}$")
    for nu_, s_, ch in zip(nus, sd, [r[1] for r in rows]):
        lab = f"{ch:.2f}" if ch < 1 else f"{ch:.0f}"
        ax[1].annotate(f"{lab} chords", (nu_, s_), textcoords="offset points",
                       xytext=(6, 4), fontsize=8, color=muted)
    ax[1].set_xlabel(r"sea multiplicity $\nu$ ($\nu = 1$ is the physical density)", color=ink)
    ax[1].set_ylabel(r"sd of one reading of $K_1$, over $|K_1|$", color=ink)
    ax[1].set_title("One reading is critically sampled", color=ink, fontsize=11)
    ax[1].legend(frameon=False, fontsize=9)
    for a in ax:
        a.grid(color=grid, lw=0.6)
        a.set_axisbelow(True)
        for s in a.spines.values():
            s.set_color(grid)
        a.tick_params(colors=muted)
    fig.tight_layout()
    name = "chord_reading_estimator.png"
    fig.savefig(output_path(name), dpi=150, bbox_inches="tight")
    dp_ = docs_path(name)
    if dp_:
        fig.savefig(dp_, dpi=150, bbox_inches="tight")
    print(f"\nWrote {name}")


if __name__ == "__main__":
    part_a()
    kq, rows = part_b()
    figure(kq, rows)
