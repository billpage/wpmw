#!/usr/bin/env python3
"""
Companion demo for ``docs/analysis/soft_core_coulomb.md``.

The soft-core Coulomb potential

    V(r) = -Z / sqrt(r^2 + eps^2)

put through the geometry that ``eckart_barrier_compensated.md`` put
``V0 sech^2(r/a)`` through: quiet points, emission lobes, the reach ceiling,
and what the ceiling costs.  Units throughout are ``hbar = mu = 1`` and the
length unit is the Bohr radius of the problem, ``a0 = hbar^2 / mu Z``, so
``Z = 1`` and ``eps`` is measured in ``a0``.

Parts
-----
A  Theorem Z1, the quiet points: ``V'''`` vanishes at ``r = 0`` and
   ``r = +- eps sqrt(3/2)``, so the soft core carries four emission lobes.
   Pure Coulomb carries two.
B  Theorem Z2, the lobes at finite reach: does the interior quiet ring
   survive the horizon, and is the nucleus dark?
C  Theorem Z3, the ceiling is position-dependent, ``R(x) = sqrt(x^2+eps^2)``.
D  Theorem Z4, the softening threshold ``eps >= k^4 pi^4 / 4``.
E  Theorem Z5, the channel is not empty at threshold.
F  Figures.

Run::

    WPMW_OUTPUT=/tmp/out PYTHONPATH=src python3 -u src/demo_soft_core_coulomb.py
"""

from __future__ import annotations

import numpy as np
from scipy.linalg import eigh_tridiagonal

from wpmwlib.wpmw_utils import docs_path, output_path

HBAR = 1.0
MU = 1.0
ROOT32 = np.sqrt(1.5)          # 1.2247449, the quiet radius in units of eps
EPS_C = np.pi ** 4 / 4.0       # 24.352, the k = 1 softening threshold


def banner(title: str) -> None:
    print("\n" + "=" * 72)
    print("  " + title)
    print("=" * 72)


# --------------------------------------------------------------------- #
# The potential                                                         #
# --------------------------------------------------------------------- #
def V(r, Z=1.0, eps=1.0):
    return -Z / np.sqrt(r * r + eps * eps)


def dV(r, Z=1.0, eps=1.0):
    return Z * r * (r * r + eps * eps) ** -1.5


def d3V(r, Z=1.0, eps=1.0):
    """Closed form of Theorem Z1."""
    u = r / eps
    return (Z / eps ** 4) * u * (6.0 * u * u - 9.0) * (1.0 + u * u) ** -3.5


def residual_lattice(Z, eps, y_max, r_grid, n_rungs=64, n_y=2048):
    """``K_q(r)`` on the momentum lattice of reach ``y_max``.

    Raised-cosine horizon, as required by the soft-horizon result of
    ``compensated_liouville_algorithm.md``; the hard window's third moment
    does not converge.
    """
    y = -y_max + 2.0 * y_max * np.arange(n_y) / n_y
    prof = np.cos(0.5 * np.pi * y / y_max) ** 2
    rr = np.atleast_1d(r_grid)[:, None]
    d = (V(rr + y[None, :], Z, eps) - V(rr - y[None, :], Z, eps)
         - 2.0 * y[None, :] * dV(rr, Z, eps)) * prof
    q = np.fft.fftfreq(n_y, d=1.0 / n_y).astype(int)
    c = np.fft.fft(d, axis=1) / n_y * ((-1.0) ** q)[None, :]
    k = np.real(-1j * c / HBAR)
    keep = np.abs(q) <= n_rungs
    order = np.argsort(q[keep])
    return q[keep][order], np.pi * HBAR / (2.0 * y_max), k[:, keep][:, order]


def gamma_of(Z, eps, y_max, r_grid, n_rungs=64):
    q, dp, kq = residual_lattice(Z, eps, y_max, r_grid, n_rungs=n_rungs)
    return np.sum(np.abs(kq[:, q != 0]), axis=1), dp


def ground_state(Z, eps, n=6001):
    """Lowest eigenpair of the 1D Hamiltonian on a box, by tridiagonal solve."""
    half = 30.0 * max(eps ** 0.75, 1.0)
    r = np.linspace(-half, half, n)
    h = r[1] - r[0]
    main = HBAR ** 2 / (MU * h * h) + V(r, Z, eps)
    off = -0.5 * HBAR ** 2 / (MU * h * h) * np.ones(n - 1)
    w, v = eigh_tridiagonal(main, off, select="i", select_range=(0, 0))
    psi = v[:, 0] / np.sqrt(np.sum(v[:, 0] ** 2) * h)
    r0 = float(np.sum(r * psi ** 2) * h)
    sig_r = float(np.sqrt(np.sum((r - r0) ** 2 * psi ** 2) * h))
    sig_p = float(np.sqrt(HBAR ** 2 * np.sum(np.gradient(psi, h) ** 2) * h))
    return float(w[0]), sig_r, sig_p


# --------------------------------------------------------------------- #
# Part A -- Z1                                                          #
# --------------------------------------------------------------------- #
def part_a():
    banner("A  Z1  the quiet points  r = 0  and  r = +- eps sqrt(3/2)")
    print("  closed form:  V''' = (Z/eps^4) u (6u^2 - 9) (1+u^2)^(-7/2),"
          "  u = r/eps\n")
    print("      Z      eps    rel err vs finite difference"
          "     root        eps sqrt(3/2)")
    for Z, eps in ((1.0, 1.0), (1.0, 0.3), (1.0, 2.5), (2.0, 0.7),
                   (-1.0, 1.0)):
        r = np.linspace(-6.0 * eps, 6.0 * eps, 200001)
        h = r[1] - r[0]
        v = V(r, Z, eps)
        fd = (v[4:] - 2 * v[3:-1] + 2 * v[1:-3] - v[:-4]) / (2 * h ** 3)
        ex = d3V(r[2:-2], Z, eps)
        err = float(np.max(np.abs(fd - ex)) / np.max(np.abs(ex)))
        lo, hi = 0.2 * eps, 4.0 * eps
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if d3V(lo, Z, eps) * d3V(mid, Z, eps) <= 0.0:
                hi = mid
            else:
                lo = mid
        print(f"   {Z:6.1f} {eps:7.2f}          {err:10.3e}"
              f"            {0.5 * (lo + hi):10.6f}   {eps * ROOT32:12.6f}")
    print(f"\n  sqrt(3/2) = {ROOT32:.6f}.  The sech^2 analogue of Theorem K7")
    print(f"  is artanh sqrt(2/3) = {np.arctanh(np.sqrt(2/3)):.6f} in units")
    print("  of a.  Here the quiet radius carries no Z at all.")
    print("\n  Pure Coulomb: V = -Z/r, V''' = -6Z/r^4, one sign on each")
    print("  half-line and no interior zero, hence two lobes rather than")
    print("  four.  The inner pair is manufactured by the softening and")
    print("  collapses onto the origin as eps -> 0.")


# --------------------------------------------------------------------- #
# Part B -- Z2                                                          #
# --------------------------------------------------------------------- #
def part_b():
    banner("B  Z2  the four lobes at finite reach")
    Z, eps = 1.0, 1.0
    r = np.linspace(0.02, 4.0, 2001)
    print(f"  Z = {Z:g}, eps = {eps:g};  predicted quiet radius "
          f"{eps * ROOT32:.6f}\n")
    print("   y_max/eps   Gamma_max     Gamma(0)    Gamma(ring)/Gamma_max"
          "   ring radius   shift")
    rows = []
    for frac in (0.1, 0.25, 0.5, 0.75, 0.9, 0.99):
        ym = frac * eps
        g, _ = gamma_of(Z, eps, ym, r)
        g0, _ = gamma_of(Z, eps, ym, np.array([0.0]))
        interior = np.where((g[1:-1] < g[:-2]) & (g[1:-1] < g[2:]))[0] + 1
        if len(interior):
            i = int(interior[0])
            d = 0.5 * (g[i - 1] - g[i + 1]) / (g[i - 1] - 2 * g[i] + g[i + 1])
            rq = r[i] + d * (r[1] - r[0])
            rows.append((frac, g.max(), float(g0[0]), g[i] / g.max(), rq))
            print(f"   {frac:8.2f}  {g.max():10.4e} {float(g0[0]):12.2e}"
                  f"      {g[i] / g.max():14.3e}   {rq:11.5f}"
                  f"  {rq / (eps * ROOT32) - 1:+7.2%}")
        else:
            print(f"   {frac:8.2f}  {g.max():10.4e} {float(g0[0]):12.2e}"
                  "      (the ring has filled in)")
    print("\n  Gamma(0) = 0 identically at every reach: the nucleus is dark,")
    print("  by parity, exactly as the sech^2 summit is (Theorem K7).  The")
    print("  interior ring survives the horizon and drifts outward with it.")
    return rows


# --------------------------------------------------------------------- #
# Part C -- Z3                                                          #
# --------------------------------------------------------------------- #
def part_c():
    banner("C  Z3  the ceiling is position-dependent: R(x) = sqrt(x^2+eps^2)")
    eps = 1.0
    print(f"  eps = {eps:g}.  Corollary K1.2 gives the ceiling; what is new")
    print("  here is that it varies, which sech^2's flat pi a / 2 cannot.\n")
    print("         x      R(x)     local dp    dp(x)/dp(0)")
    for x in (0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 64.0):
        R = np.sqrt(x * x + eps * eps)
        print(f"    {x:7.2f} {R:9.4f}  {np.pi * HBAR / (2 * R):10.5f}"
              f"   {1.0 / R:12.5f}")
    print("\n  Two readings, and they are not equivalent.")
    print("   (i)  a local reach makes the momentum quantum coarse at the")
    print("        nucleus and fine in the far field, so the phase-space")
    print("        crystal is not uniform.  Logged as Z-LS1; not pursued.")
    print("   (ii) a uniform reach must take the infimum, y_max < eps, and")
    print("        that is what Part D prices.")


# --------------------------------------------------------------------- #
# Part D -- Z4                                                          #
# --------------------------------------------------------------------- #
def part_d():
    banner("D  Z4  the softening threshold  eps >= k^4 pi^4 / 4  a0")
    print("  Uniform reach y_max < eps, so dp > pi hbar / 2 eps.  The well")
    print("  is harmonic for eps >> a0 with omega = sqrt(Z / mu eps^3), so")
    print("  sigma_r = eps^(3/4)/sqrt2 and sigma_p = 1/(sqrt2 eps^(3/4)).\n")
    print("        eps        E_0     sigma_r   eps^(3/4)/sqrt2   sigma_p"
          "    sigma_p/dp")
    rows = []
    for eps in (1.0, 5.0, EPS_C, 50.0, 100.0, 4 * np.pi ** 4, 1000.0):
        e0, sr, sp = ground_state(1.0, eps)
        dp = np.pi * HBAR / (2.0 * eps)
        rows.append((eps, e0, sr, sp, sp / dp))
        print(f"   {eps:9.3f} {e0:10.6f} {sr:10.3f} {eps**0.75/np.sqrt(2):16.3f}"
              f" {sp:10.6f}   {sp / dp:10.3f}")
    print(f"\n  predicted thresholds:  k = 1  ->  pi^4/4    = {EPS_C:9.3f}")
    print(f"                         k = 2  ->  4 pi^4    = "
          f"{4 * np.pi**4:9.1f}")
    print(f"                         k = 3  ->  81 pi^4/4 = "
          f"{81 * np.pi**4 / 4:9.1f}")
    # measured crossings, by bisection on the eigenproblem
    for k in (1, 2):
        lo, hi = 5.0, 4000.0
        for _ in range(28):
            mid = np.sqrt(lo * hi)
            _, _, sp = ground_state(1.0, mid)
            if sp / (np.pi * HBAR / (2.0 * mid)) >= k:
                hi = mid
            else:
                lo = mid
        meas = np.sqrt(lo * hi)
        pred = k ** 4 * np.pi ** 4 / 4.0
        print(f"  measured crossing, k = {k}:  eps = {meas:8.2f}"
              f"   ratio to predicted {meas / pred:5.3f}")
    print("\n  The ratio falls towards 1 as the anharmonic correction to")
    print("  sigma_r dies, which is the only approximation in Z4.")
    return rows


# --------------------------------------------------------------------- #
# Part E -- Z5                                                          #
# --------------------------------------------------------------------- #
def part_e():
    banner("E  Z5  the channel is not empty at threshold")
    print("  At y_max = eps the horizon spans eps/sigma_r ~ eps^(1/4) times")
    print("  the state's own width, so it always reaches past the harmonic")
    print("  core.  A harmonic well would give Gamma = 0 exactly (I4).\n")
    print("        eps   eps/sigma_r    Gamma_max   Gamma(sigma_r)"
          "   ratio")
    rows = []
    for eps in (EPS_C, 100.0, 4 * np.pi ** 4, 1000.0):
        _, sr, _ = ground_state(1.0, eps)
        r = np.linspace(0.0, 4.0 * sr, 801)
        g, _ = gamma_of(1.0, eps, eps, r, n_rungs=96)
        gs = float(np.interp(sr, r, g))
        rows.append((eps, eps / sr, g.max(), gs))
        print(f"   {eps:9.3f}   {eps / sr:10.3f}   {g.max():10.4e}"
              f"   {gs:12.4e}   {gs / g.max():7.4f}")
    print("\n  Non-zero throughout.  So the uniform-lattice threshold does")
    print("  not push the problem into the regime where the demographic")
    print("  channel is empty -- which is what happens on the Eckart")
    print("  barrier, where K-LS2's two conditions have no common ground.")
    return rows


# --------------------------------------------------------------------- #
# Part F -- figures                                                     #
# --------------------------------------------------------------------- #
def save(fig, name):
    fig.savefig(output_path(name), dpi=150, bbox_inches="tight")
    dp = docs_path(name)
    if dp:
        fig.savefig(dp, dpi=150, bbox_inches="tight")
    print(f"  wrote {name}")


def part_f(d_rows):
    banner("F  figures")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 3, figsize=(13.5, 3.9))

    eps = 1.0
    r = np.linspace(-4.0, 4.0, 1601)
    ax[0].plot(r, d3V(r, 1.0, eps) / np.max(np.abs(d3V(r, 1.0, eps))),
               "k-", lw=1.3, label=r"$V'''$, normalised")
    g, _ = gamma_of(1.0, eps, 0.5 * eps, r)
    ax[0].plot(r, g / g.max(), "C0-", lw=1.3,
               label=r"$\Gamma$, $y_{max}=\epsilon/2$")
    for rr in (-eps * ROOT32, 0.0, eps * ROOT32):
        ax[0].axvline(rr, color="0.6", ls=":", lw=1.0)
    ax[0].set_xlabel(r"$r/\epsilon$")
    ax[0].set_title("(a) Z1: four lobes, three quiet points", fontsize=10)
    ax[0].legend(fontsize=8)
    ax[0].grid(alpha=0.3)

    x = np.linspace(0.0, 8.0, 400)
    ax[1].plot(x, np.pi / (2 * np.sqrt(x * x + 1.0)), "C1-", lw=1.4)
    ax[1].axhline(np.pi / 2, color="k", ls="--", lw=1.0,
                  label=r"$\pi\hbar/2\epsilon$, the uniform floor")
    ax[1].set_xlabel(r"$x/\epsilon$")
    ax[1].set_ylabel(r"$\Delta p(x)$")
    ax[1].set_title("(b) Z3: a ceiling that moves", fontsize=10)
    ax[1].legend(fontsize=8)
    ax[1].grid(alpha=0.3)

    e = np.array([row[0] for row in d_rows])
    ratio = np.array([row[4] for row in d_rows])
    ax[2].loglog(e, ratio, "o-", ms=6, color="C2")
    ax[2].axhline(1.0, color="k", ls="--", lw=1.0)
    ax[2].axhline(2.0, color="0.5", ls=":", lw=1.0)
    ax[2].axvline(EPS_C, color="k", ls="--", lw=0.8)
    ax[2].axvline(4 * np.pi ** 4, color="0.5", ls=":", lw=0.8)
    ax[2].set_xlabel(r"$\epsilon/a_0$")
    ax[2].set_ylabel(r"rungs per $\sigma_p$")
    ax[2].set_title(r"(c) Z4: $\epsilon_c = k^4\pi^4/4$", fontsize=10)
    ax[2].grid(alpha=0.3, which="both")

    fig.tight_layout()
    save(fig, "soft_core_coulomb_geometry.png")
    plt.close(fig)


def main():
    part_a()
    part_b()
    part_c()
    d_rows = part_d()
    part_e()
    part_f(d_rows)
    banner("summary")
    print("  Z1  quiet points at r = 0 and r = +- eps sqrt(3/2); four lobes")
    print("  Z2  the nucleus is dark at every reach; the ring survives")
    print("  Z3  the reach ceiling is position-dependent")
    print("  Z4  a uniform lattice needs eps >= k^4 pi^4 / 4 a0")
    print("  Z5  and the channel is not empty there")


if __name__ == "__main__":
    main()
