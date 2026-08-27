"""
Verification for ``docs/analysis/reach_energy_coupling.md``.

The reach ``y_max`` enters the microdynamics as a window on the ket-bra half
separation ``y``.  This script asks what that window actually controls.  The
answer separates cleanly into three groups:

  * conservation  -- untouched by the reach, for any even window;
  * the leading quantum term -- untouched, but only under the compensated
    split; contaminated at order 1 / y_max^2 without it;
  * the event budget -- always diverges with the reach.

Kernel conventions follow ``open_position_space.md`` Theorem O1 and
``compensated_liouville_splitting.md`` section 2:

    D(x, y)     = V(x + y) - V(x - y)
    D_res(x, y) = D(x, y) - 2 y V'(x)          (the compensated residual)
    K_q         = Fourier coefficient of D over the period 2 y_max, / (i hbar)
    xi_q        = q dp,     dp = pi hbar / (2 y_max)

so that the collision term is  sum_q K_q W(p - xi_q), and the moments are

    M_n = sum_q xi_q^n K_q,
    M_0 = 0 (worlds),  M_1 = -V' (force),  M_2 = 0 (energy),
    M_3 = (hbar^2 / 4) V'''  (the leading Moyal coefficient).

Parts
-----
A  Theorem E1.  Delta p = pi hbar / (period of D in y).  The three sources of
   that period -- ring, periodic V, postulated horizon -- are one mechanism.
   A sharp window commensurate with a comb potential is exact; an
   incommensurate one is not.  This corrects ``open_position_space.md`` 3.2,
   whose divergence is an off-lattice artefact.
B  Theorem E2.  Folding the tail back in, rather than discarding it, gives
   exactly the kernel of the periodised potential sampled on the lattice.
   The price is Proposition O4: folded activity does not localise.
C  Theorem E3.  M_0 = M_2 = 0 for every even window -- hard, soft or folded.
   Worlds and energy are conserved independently of the reach.
D  Theorem E4.  The four-action energy ledger with more than one mode: the
   focus channel does no net work, the hop channel delivers exactly the
   classical power.  Cyganski's single-harmonic result extends unchanged.
E  Theorem E5.  Applied to the *full* kernel, a profile w contaminates the
   third moment with the classical force:
       M_3 = (hbar^2/4) V''' + (3 pi^2 hbar^2 / 8 y_max^2) V'.
F  Theorem E6.  Applied to the *compensated residual* the contamination
   vanishes identically, because D_res(0) = D_res'(0) = D_res''(0) = 0.
   The reach is inert at order hbar^2 in the compensated algorithm.
G  Theorem E7.  For a polynomial V the Moyal series terminates, so the
   collision operator is a finite-order differential operator in p and not a
   jump measure at all.  On the open line the microdynamics does not exist;
   the reach is constitutive, not approximate.
H  Theorem E8.  The compensated residual budget diverges with the reach for
   every V with V'(x) non-zero -- linearly for bounded V, cubically for the
   quartic.  What the reach buys is momentum resolution, not accuracy.

Run with ``WPMW_OUTPUT`` set to a writable directory; set ``WPMW_DOCS`` as
well to mirror keeper figures into the ``output``-branch worktree.
"""

from __future__ import annotations

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from wpmwlib.wpmw_utils import output_path, docs_path  # noqa: E402

HBAR = 1.0
MASS = 1.0
RNG = np.random.default_rng(20260827)


def banner(text):
    print("\n" + "=" * 72)
    print(text)
    print("=" * 72)


def save_fig(fig, name):
    fig.savefig(output_path(name), dpi=150, bbox_inches="tight")
    dp = docs_path(name)
    if dp:
        fig.savefig(dp, dpi=150, bbox_inches="tight")
    print(f"    [figure] {name}")


# --------------------------------------------------------------------- #
# Potentials.  Each carries V, V', V''' and a reference position.        #
# --------------------------------------------------------------------- #
SIGMA = 0.5
LRING = 8.0
K1 = 2.0 * np.pi / LRING
KGOLD = K1 * (1.0 + np.sqrt(5.0)) / 2.0
AQ, BQ = 0.05, 0.8


def gaussian(z):
    return np.exp(-z ** 2 / (2 * SIGMA ** 2))


def d_gaussian(z):
    return -z / SIGMA ** 2 * gaussian(z)


def d3_gaussian(z):
    u = z / SIGMA
    return -(1.0 / SIGMA ** 3) * (u ** 3 - 3 * u) * np.exp(-u ** 2 / 2)


def comb2(z):
    return np.cos(K1 * z) + 0.6 * np.cos(2 * K1 * z + 0.4)


def d_comb2(z):
    return -K1 * np.sin(K1 * z) - 1.2 * K1 * np.sin(2 * K1 * z + 0.4)


def d3_comb2(z):
    return K1 ** 3 * np.sin(K1 * z) + 4.8 * K1 ** 3 * np.sin(2 * K1 * z + 0.4)


def incomm2(z):
    return np.cos(K1 * z) + 0.6 * np.cos(KGOLD * z)


def d_incomm2(z):
    return -K1 * np.sin(K1 * z) - 0.6 * KGOLD * np.sin(KGOLD * z)


def d3_incomm2(z):
    return K1 ** 3 * np.sin(K1 * z) + 0.6 * KGOLD ** 3 * np.sin(KGOLD * z)


def quartic(z):
    return AQ * z ** 4 - BQ * z ** 2


def d_quartic(z):
    return 4 * AQ * z ** 3 - 2 * BQ * z


def d3_quartic(z):
    return 24 * AQ * z


CASES = {
    "Gaussian barrier": (gaussian, d_gaussian, d3_gaussian, 1.0),
    "commensurate 2-mode": (comb2, d_comb2, d3_comb2, 1.3),
    "incommensurate 2-mode": (incomm2, d_incomm2, d3_incomm2, 1.3),
    "quartic double well": (quartic, d_quartic, d3_quartic, 1.7),
}


# --------------------------------------------------------------------- #
# The lattice kernel                                                     #
# --------------------------------------------------------------------- #
def profile(y, y_max, mode):
    """Even window applied to the separation coordinate."""
    if mode == "soft":
        return np.cos(np.pi * (y / y_max) / 2) ** 2
    return np.ones_like(y)


def kernel(vfun, dvfun, x, y_max, n_p, mode="hard", compensated=False,
           n_images=0):
    """Return (xi, dp, K) on the momentum lattice of reach ``y_max``.

    ``mode`` is ``hard``, ``soft`` or ``fold``.  ``fold`` replaces truncation
    by periodic summation of the full-line D over the period 2 y_max.
    """
    period = 2.0 * y_max
    dp = np.pi * HBAR / period
    y = -period / 2 + period * np.arange(n_p) / n_p
    q = np.fft.fftfreq(n_p, d=1.0 / n_p).astype(int)
    xi = q * dp

    if mode == "fold":
        d = np.zeros_like(y)
        for n in range(-n_images, n_images + 1):
            yy = y + n * period
            d = d + vfun(x + yy) - vfun(x - yy)
    else:
        d = vfun(x + y) - vfun(x - y)
        if compensated:
            d = d - 2.0 * y * dvfun(x)
        d = d * profile(y, y_max, mode)

    # y grid starts at -period/2, so the DFT picks up a (-1)^q phase.
    coeff = np.fft.fft(d) / n_p * ((-1.0) ** q)
    return xi, dp, coeff / (1j * HBAR)


def moment(xi, kk, n):
    return float(np.real(np.sum(xi ** n * kk)))


def budget(kk):
    return float(np.sum(np.abs(kk)))


# --------------------------------------------------------------------- #
# Part A -- Theorem E1                                                   #
# --------------------------------------------------------------------- #
def part_a():
    banner("Part A  Theorem E1: the reach is a period, not an aperture")
    print("  dp = pi hbar / (period of D_x in y).  Three sources of the")
    print("  period; a sharp window is exact when commensurate with it.\n")

    x = 0.7
    kk = 2 * np.pi / 4.0                       # V = cos(2 pi x / 4)
    vcos = lambda z: np.cos(kk * z)            # noqa: E731
    dcos = lambda z: -kk * np.sin(kk * z)      # noqa: E731
    exact = 2.0 * abs(np.sin(kk * x))
    print(f"  V = cos(2 pi x / 4), x = {x}: exact budget 2|Gamma| = {exact:.6f}")
    print(f"  {'L_c':>6} {'L_c/a':>8} {'commens':>9} {'window':>7} "
          f"{'budget':>10} {'ratio':>9}")
    rows = []
    for l_c in (4.0, 8.0, 6.0, 5.0):
        commens = abs(l_c / 4.0 - round(l_c / 4.0)) < 1e-12
        for mode in ("hard", "soft"):
            _, _, kq = kernel(vcos, dcos, x, l_c / 2, 8192, mode=mode)
            b = budget(kq)
            rows.append((l_c, commens, mode, b / exact))
            print(f"  {l_c:>6.1f} {l_c / 4.0:>8.2f} {str(commens):>9} "
                  f"{mode:>7} {b:>10.6f} {b / exact:>9.4f}")
    print("\n  And the exactness does not depend on the rung count -- the")
    print("  sidelobe spacing pi hbar / L_c IS dp, so every sidelobe lands on")
    print("  a lattice zero (L_c = 8 = 2a):")
    print(f"  {'N_p':>8} {'on-lattice budget':>20} {'ratio':>12}")
    for n_p in (16, 64, 256, 1024, 4096):
        _, _, kq = kernel(vcos, dcos, x, 4.0, n_p, mode="hard")
        b = budget(kq)
        print(f"  {n_p:>8} {b:>20.6f} {b / exact:>12.6f}")

    print("\n  A commensurate sharp window is EXACT: the sidelobes of the")
    print("  window fall on lattice zeros because their spacing IS dp.")
    print("  open_position_space.md 3.2 measured the continuous-xi integral,")
    print("  which the model never evaluates.  The soft profile is not exact")
    print("  even when commensurate -- it attenuates the true delta.")
    return rows


# --------------------------------------------------------------------- #
# Part B -- Theorem E2                                                   #
# --------------------------------------------------------------------- #
def part_b():
    banner("Part B  Theorem E2: folding is exact; folding is not local")
    print("  K_fold = dp * V_W(x, xi_q) for the periodised potential V_a,")
    print("  a = 2 y_max.  Folding D in y equals periodising V in x.\n")

    def vw_exact(x, xi):
        vt = SIGMA * np.sqrt(2 * np.pi) * np.exp(-SIGMA ** 2 * (2 * xi) ** 2 / 2)
        return (2.0 / (np.pi * HBAR ** 2)) * vt * np.sin(2 * xi * x / HBAR)

    print(f"  {'x':>6} {'y_max':>7} {'max|K_fold - dp V_W| / max|K|':>32}")
    for x in (0.5, 1.0, 3.0):
        for y_max in (2.0, 4.0):
            xi, dp, kf = kernel(gaussian, d_gaussian, x, y_max, 2048,
                                mode="fold", n_images=60)
            ke = dp * vw_exact(x, xi)
            scale = max(np.max(np.abs(ke)), 1e-300)
            print(f"  {x:>6.1f} {y_max:>7.1f} "
                  f"{np.max(np.abs(kf - ke)) / scale:>32.2e}")

    print("\n  Force delivered by the folded kernel is the image-corrected")
    print("  force -V_a'(x), and nothing else:")
    print(f"  {'y_max':>7} {'M1 fold':>13} {'-V(x) prime':>13} "
          f"{'-sum_n V prime(x+na)':>22}")
    for y_max in (2.0, 4.0, 8.0):
        a = 2 * y_max
        xi, _, kf = kernel(gaussian, d_gaussian, 1.0, y_max, 4096,
                           mode="fold", n_images=60)
        pred = -sum(d_gaussian(1.0 + n * a) for n in range(-60, 61))
        print(f"  {y_max:>7.1f} {moment(xi, kf, 1):>+13.8f} "
              f"{-d_gaussian(1.0):>+13.8f} {pred:>+22.8f}")

    print("\n  But the images are everywhere, so Proposition O4 fails:")
    print(f"  {'x':>6} {'hard':>12} {'soft':>12} {'fold':>12} {'|V(x)|':>11}")
    loc = []
    for x in (1.0, 3.0, 5.0, 7.0, 9.0, 11.0):
        row = []
        for mode in ("hard", "soft", "fold"):
            _, _, kq = kernel(gaussian, d_gaussian, x, 2.0, 2048, mode=mode,
                              n_images=60)
            row.append(budget(kq))
        loc.append((x, *row, float(gaussian(x))))
        print(f"  {x:>6.1f} {row[0]:>12.3e} {row[1]:>12.3e} {row[2]:>12.3e} "
              f"{gaussian(x):>11.2e}")
    print("\n  Locality, exactness, lattice: pick two, unless V is periodic.")
    return loc


# --------------------------------------------------------------------- #
# Part C -- Theorem E3                                                   #
# --------------------------------------------------------------------- #
def part_c():
    banner("Part C  Theorem E3: conservation does not depend on the reach")
    print("  D is odd in y and every window here is even, so the kernel stays")
    print("  odd in xi.  M_0 = 0 (worlds) and M_2 = 0 (energy) follow.\n")
    print(f"  {'potential':>22} {'y_max':>6} {'window':>7} {'|M0|':>10} "
          f"{'|M2|':>10} {'M1':>13} {'-V prime':>13}")
    worst0 = worst2 = 0.0
    for name, (vf, dvf, _, x) in CASES.items():
        for y_max in (1.0, 4.0):
            for mode in ("hard", "soft"):
                xi, _, kq = kernel(vf, dvf, x, y_max, 4096, mode=mode)
                m0, m1, m2 = (abs(moment(xi, kq, 0)), moment(xi, kq, 1),
                              abs(moment(xi, kq, 2)))
                worst0, worst2 = max(worst0, m0), max(worst2, m2)
                print(f"  {name:>22} {y_max:>6.1f} {mode:>7} {m0:>10.1e} "
                      f"{m2:>10.1e} {m1:>+13.7f} {-dvf(x):>+13.7f}")

    print("\n  Folding is defined only for a decaying V -- the periodic")
    print("  summation of a periodic or unbounded potential does not")
    print("  converge -- so it gets its own rows.  Its M1 is the")
    print("  image-corrected force -V_a'(x), not -V'(x):")
    for y_max in (2.0, 4.0, 8.0):
        xi, _, kq = kernel(gaussian, d_gaussian, 1.0, y_max, 4096,
                           mode="fold", n_images=40)
        m0, m1, m2 = (abs(moment(xi, kq, 0)), moment(xi, kq, 1),
                      abs(moment(xi, kq, 2)))
        worst0, worst2 = max(worst0, m0), max(worst2, m2)
        print(f"  {'Gaussian barrier':>22} {y_max:>6.1f} {'fold':>7} "
              f"{m0:>10.1e} {m2:>10.1e} {m1:>+13.7f} "
              f"{-d_gaussian(1.0):>+13.7f}")

    print(f"\n  worst |M0| = {worst0:.1e}   worst |M2| = {worst2:.1e}")
    print("  This closes open item 4 of open_position_space.md: energy under")
    print("  the horizon is exact, for the same reason worlds are.")


# --------------------------------------------------------------------- #
# Part D -- Theorem E4                                                   #
# --------------------------------------------------------------------- #
def part_d():
    banner("Part D  Theorem E4: the four-action energy ledger, two modes")
    print("  Symmetric member of the four-rule family, per mode q:")
    print("    f_n = (Gamma/2)(W_{n+q} - W_{n-q})   focus")
    print("    h_n = -(Gamma/2)(W_{n+q} + W_{n-q})  hop")
    print("  Focus conserves momentum and changes T by -xi_q^2/m per event;")
    print("  hop moves one world by 2 xi_q, changing T by 2 p_n xi_q / m.\n")

    n_x, n_p = 64, 96
    dx = LRING / n_x
    xs = np.arange(n_x) * dx
    dp = np.pi * HBAR / LRING
    ps = (np.arange(n_p) - n_p // 2) * dp

    gam = {}
    for i, xv in enumerate(xs):
        _, _, kq = kernel(comb2, d_comb2, xv, LRING / 2, 256)
        idx = np.fft.fftfreq(256, d=1.0 / 256).astype(int)
        gam[i] = {q: float(np.real(kq[np.where(idx == -q)[0][0]]))
                  for q in (1, 2)}

    results = []
    for sig_p, tag in ((3.0, "packet inside the p grid"),
                       (0.9, "packet at the p-grid edge")):
        xx, pp = np.meshgrid(xs, ps, indexing="ij")
        w2d = np.exp(-((xx - 3.0) ** 2) / (2 * 0.8 ** 2)
                     - (pp - 0.3) ** 2 / (2 * (sig_p * dp) ** 2))
        w2d /= np.sum(w2d) * dx * dp

        dt_focus = dt_hop = 0.0
        for i in range(n_x):
            col = w2d[i]
            for q in (1, 2):
                g = gam[i][q]
                w_plus, w_minus = np.roll(col, -q), np.roll(col, q)
                f_n = 0.5 * g * (w_plus - w_minus)
                h_n = -0.5 * g * (w_plus + w_minus)
                xi_q = q * dp
                dt_focus += np.sum(f_n) * (-xi_q ** 2 / MASS)
                dt_hop += np.sum(h_n * (2 * ps * xi_q / MASS))
        dt_focus *= dx * dp
        dt_hop *= dx * dp
        power = -np.sum(w2d * (ps[None, :] / MASS)
                        * d_comb2(xs)[:, None]) * dx * dp
        results.append((tag, dt_focus, dt_hop, power))
        print(f"  {tag}  (sigma_p = {sig_p:.1f} dp)")
        print(f"    focus channel dT/dt = {dt_focus:+.3e}   (should be zero)")
        print(f"    hop   channel dT/dt = {dt_hop:+.8f}")
        print(f"    classical power     = {power:+.8f}")
        print(f"    residual            = {dt_focus + dt_hop - power:+.2e}\n")
    print("  sum_n f_n telescopes, so focus does no net work; and")
    print("  sum_n h_n p_n = -Gamma <p>, so the hop channel alone carries the")
    print("  classical power.  Both statements are linear in the modes, so")
    print("  they hold for any comb.  Caveat: the momentum lattice here is")
    print("  cyclic, which makes the telescoping exact by construction; an")
    print("  open, capped lattice is not tested.")
    return results


# --------------------------------------------------------------------- #
# Parts E and F -- Theorems E5, E6                                       #
# --------------------------------------------------------------------- #
def m3_prediction(dvf, d3vf, x, y_max):
    """Full-kernel soft third moment, Theorem E5."""
    return (HBAR ** 2 / 4 * d3vf(x)
            - 3 * np.pi ** 2 * HBAR ** 2 / (8 * y_max ** 2) * dvf(x))


def part_ef():
    banner("Parts E and F  Theorems E5, E6: where the reach touches hbar^2")
    print("  sum_q xi^3 K = (hbar^2/8) (w D)'''(0), and")
    print("    full:     (w D)'''(0)     = 2 V''' + 3 w''(0) * 2 V'")
    print("    residual: (w D_res)'''(0) = 2 V'''   [D_res' (0) = 0]")
    print("  so the profile contaminates the full kernel and not the")
    print("  compensated residual.\n")

    y_grid = np.array([0.5, 1.0, 2.0, 4.0, 8.0])
    collected = {}
    for name in ("Gaussian barrier", "incommensurate 2-mode",
                 "quartic double well"):
        vf, dvf, d3vf, x = CASES[name]
        target = HBAR ** 2 / 4 * d3vf(x)
        print(f"  {name}  x = {x}   exact (hbar^2/4) V''' = {target:+.8f}"
              f"   -V' = {-dvf(x):+.6f}")
        print(f"    {'y_max':>7} {'full soft M3':>14} {'E5 predicted':>14} "
              f"{'resid soft M3':>15} {'rel err':>10} {'R_res':>10}")
        full, pred, resid = [], [], []
        for y_max in y_grid:
            xf, _, kf = kernel(vf, dvf, x, y_max, 16384, mode="soft")
            xr, _, kr = kernel(vf, dvf, x, y_max, 16384, mode="soft",
                               compensated=True)
            m3f, m3r = moment(xf, kf, 3), moment(xr, kr, 3)
            pr = m3_prediction(dvf, d3vf, x, y_max)
            full.append(m3f)
            pred.append(pr)
            resid.append(m3r)
            print(f"    {y_max:>7.2f} {m3f:>14.6f} {pr:>14.6f} "
                  f"{m3r:>15.8f} {abs(m3r / target - 1):>10.2e} "
                  f"{budget(kr):>10.4f}")
        collected[name] = (np.array(full), np.array(pred), np.array(resid),
                           target)
        print()

    print("  Convergence in N_p of the compensated residual third moment")
    print("  (Gaussian barrier, y_max = 2).  The raised cosine gives a q^-3")
    print("  tail, so the terms xi^3 K are O(1) and cancel by sign: the sum")
    print("  converges conditionally, not absolutely.")
    vf, dvf, d3vf, x = CASES["Gaussian barrier"]
    target = HBAR ** 2 / 4 * d3vf(x)
    print(f"    {'N_p':>8} {'hard M3_res':>14} {'soft M3_res':>14} "
          f"{'target':>14} {'R hard':>9} {'R soft':>9}")
    for n_p in (256, 1024, 4096, 16384, 65536):
        xh, _, kh = kernel(vf, dvf, x, 2.0, n_p, mode="hard", compensated=True)
        xs_, _, ks = kernel(vf, dvf, x, 2.0, n_p, mode="soft",
                            compensated=True)
        print(f"    {n_p:>8} {moment(xh, kh, 3):>14.4f} "
              f"{moment(xs_, ks, 3):>14.8f} {target:>14.8f} "
              f"{budget(kh):>9.4f} {budget(ks):>9.4f}")
    print("\n  The hard cutoff on the residual still fails as the compensated")
    print("  Liouville spec 4.4 says.  Under the soft profile the residual")
    print("  third moment is reach-independent to 1e-4, and the agreement is")
    print("  not monotone in N_p -- quote it to four figures, not six.")
    return y_grid, collected


# --------------------------------------------------------------------- #
# Part G -- Theorem E7                                                   #
# --------------------------------------------------------------------- #
def part_g():
    banner("Part G  Theorem E7: a polynomial V has no jump measure")
    print(f"  V = {AQ} x^4 - {BQ} x^2  =>  V'''' = {24 * AQ:.2f} (constant),"
          f"  V''''' = 0")
    print("  so the Moyal series TERMINATES and the collision term is")
    print("      V' d_p W - (hbar^2/24) V''' d_p^3 W")
    print("  i.e.  V_W = V' delta'(xi) + (hbar^2/24) V''' delta'''(xi).")
    print("  A finite-order differential operator has unbounded total")
    print("  variation, so there is no finite jump rate on the open line:")
    print("  the four-action microdynamics does not exist for a polynomial V")
    print("  until a reach is imposed.  The reach is constitutive here, not")
    print("  approximate.\n")
    print("  The termination also makes this the cleanest available")
    print("  benchmark: the exact reference is a three-term PDE with no")
    print("  truncation ambiguity.  Residual moments against it:")
    vf, dvf, d3vf, x = CASES["quartic double well"]
    print(f"    {'y_max':>7} {'M1_res':>11} {'M3_res':>13} "
          f"{'(h^2/4)V prime prime prime':>28} {'R_res':>12}")
    for y_max in (1.0, 2.0, 4.0, 8.0, 16.0):
        xi, _, kr = kernel(vf, dvf, x, y_max, 16384, mode="soft",
                           compensated=True)
        print(f"    {y_max:>7.1f} {moment(xi, kr, 1):>11.1e} "
              f"{moment(xi, kr, 3):>13.6f} {HBAR ** 2 / 4 * d3vf(x):>28.6f} "
              f"{budget(kr):>12.4f}")


# --------------------------------------------------------------------- #
# Part H -- Theorem E8                                                   #
# --------------------------------------------------------------------- #
def part_h():
    banner("Part H  Theorem E8: the residual budget always diverges")
    print("  M_res -> -i V'(x) s as y grows (the anti-drift of the")
    print("  compensated-splitting note), so the compensation ramp -2 y V'")
    print("  makes R_res grow without bound for every V with V'(x) non-zero.")
    print("  Growth factor per doubling of the reach:\n")
    names = ["Gaussian barrier", "commensurate 2-mode", "quartic double well"]
    y_grid = np.array([2.0, 4.0, 8.0, 16.0, 32.0])
    table = {}
    header = f"  {'y_max':>7} " + " ".join(f"{n:>26}" for n in names)
    print(header)
    prev = {n: None for n in names}
    for y_max in y_grid:
        cells = []
        for n in names:
            vf, dvf, _, x = CASES[n]
            _, _, kr = kernel(vf, dvf, x, y_max, 32768, mode="soft",
                              compensated=True)
            r = budget(kr)
            table.setdefault(n, []).append(r)
            ratio = "" if prev[n] is None else f" (x{r / prev[n]:.2f})"
            prev[n] = r
            cells.append(f"{r:>14.4f}{ratio:>12}")
        print(f"  {y_max:>7.1f} " + " ".join(cells))
    print("\n  Bounded V settles to x2 -- linear, from the ramp.  The quartic")
    print("  holds x8 exactly -- cubic, Theorem C2's y_max^3 sup|V'''| with no")
    print("  saturation.  Since conservation and the hbar^2 term are both")
    print("  reach-independent under compensation, the reach buys momentum")
    print("  resolution dp = pi hbar / 2 y_max and costs event budget.  It")
    print("  buys no accuracy at all.")

    print("\n  Cost of compensation for a comb: the ramp is not periodic in")
    print("  y, so subtracting it turns a finite delta comb into a dense")
    print("  spectrum.")
    vf, dvf, d3vf, x = CASES["commensurate 2-mode"]
    print(f"    {'y_max':>7} {'full nnz':>10} {'R full':>10} "
          f"{'resid nnz':>11} {'R resid':>10}")
    comb = []
    for y_max in (4.0, 8.0, 16.0, 32.0):
        _, _, kf = kernel(vf, dvf, x, y_max, 16384, mode="hard")
        _, _, kr = kernel(vf, dvf, x, y_max, 16384, mode="soft",
                          compensated=True)
        nzf = int(np.sum(np.abs(kf) > 1e-9))
        nzr = int(np.sum(np.abs(kr) > 1e-9))
        comb.append((y_max, nzf, budget(kf), nzr, budget(kr)))
        print(f"    {y_max:>7.1f} {nzf:>10} {budget(kf):>10.5f} "
              f"{nzr:>11} {budget(kr):>10.5f}")
    print("\n  On event count alone the crossover sits near y_max = 4.2, which")
    print("  is essentially the ring's own reach L/2 = 4.  Whether")
    print("  compensation still pays past it depends on variance in the force")
    print("  channel, which this script does not measure.")
    return y_grid, table, comb


# --------------------------------------------------------------------- #
# Figures                                                                #
# --------------------------------------------------------------------- #
def fig_moments(y_grid, collected):
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.0))
    for ax, name in zip(axes, collected):
        full, pred, resid, target = collected[name]
        ax.loglog(y_grid, np.abs(full - target), "o-", color="#D85A30",
                  lw=1.8, label="full kernel + profile")
        ax.loglog(y_grid, np.abs(pred - target), "--", color="#8A8A8A",
                  lw=1.4, label=r"E5: $3\pi^2\hbar^2 |V'|/8y_{\max}^2$")
        ax.loglog(y_grid, np.abs(resid - target), "s-", color="#1D9E75",
                  lw=1.8, label="compensated residual")
        ax.set_xlabel(r"reach $y_{\max}$")
        ax.set_ylim(1e-6, 1e2)
        ax.set_title(name, fontsize=10)
        ax.grid(alpha=0.3, which="both")
    axes[0].set_ylabel(r"$|\sum_q \xi^3 K_q - (\hbar^2/4)V'''|$")
    axes[0].legend(fontsize=8, loc="lower left")
    fig.suptitle("The reach touches the $\\hbar^2$ term only when the profile "
                 "is applied to the uncompensated kernel", fontsize=11)
    fig.tight_layout()
    save_fig(fig, "reach_energy_coupling_moments.png")
    plt.close(fig)


def fig_budget(y_grid, table):
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    colours = {"Gaussian barrier": "#2B5EA8",
               "commensurate 2-mode": "#1D9E75",
               "quartic double well": "#D85A30"}
    for name, vals in table.items():
        ax.loglog(y_grid, vals, "o-", color=colours[name], lw=1.8, label=name)
    ref = np.array(table["quartic double well"])
    ax.loglog(y_grid, ref[0] * (y_grid / y_grid[0]) ** 3, "--",
              color="#8A8A8A", lw=1.2, label=r"$y_{\max}^3$")
    base = np.array(table["Gaussian barrier"])
    ax.loglog(y_grid, base[-1] * (y_grid / y_grid[-1]), ":",
              color="#8A8A8A", lw=1.2, label=r"$y_{\max}^1$")
    ax.set_xlabel(r"reach $y_{\max}$")
    ax.set_ylabel(r"residual event budget $\sum_q |K_{\rm res}|$")
    ax.set_title("The reach buys momentum resolution and costs budget",
                 fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    save_fig(fig, "reach_energy_coupling_budget.png")
    plt.close(fig)


def fig_windows():
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    x, y_max, n_p = 1.0, 2.0, 512
    period = 2 * y_max
    y = -period / 2 + period * np.arange(n_p) / n_p
    d_hard = gaussian(x + y) - gaussian(x - y)
    d_soft = d_hard * profile(y, y_max, "soft")
    d_fold = np.zeros_like(y)
    for n in range(-40, 41):
        yy = y + n * period
        d_fold += gaussian(x + yy) - gaussian(x - yy)

    ax = axes[0]
    ax.plot(y, d_hard, color="#D85A30", lw=1.8, label="hard (truncate)")
    ax.plot(y, d_soft, color="#1D9E75", lw=1.8, label="soft (taper)")
    ax.plot(y, d_fold, color="#2B5EA8", lw=1.8, label="fold (periodic sum)")
    ax.axvline(-y_max, color="#8A8A8A", lw=1.0, ls=":")
    ax.axvline(y_max - period / n_p, color="#8A8A8A", lw=1.0, ls=":")
    ax.set_xlabel(r"half separation $y$")
    ax.set_ylabel(r"$D_x(y)$ over one period")
    ax.set_title("Hard truncation leaves a seam at the period edge",
                 fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[1]
    for mode, colour in (("hard", "#D85A30"), ("soft", "#1D9E75"),
                         ("fold", "#2B5EA8")):
        xi, _, kq = kernel(gaussian, d_gaussian, x, y_max, 4096, mode=mode,
                           n_images=40)
        order = np.argsort(xi)
        sel = (xi[order] > 0)
        ax.loglog(xi[order][sel], np.abs(kq[order][sel]) + 1e-30, "-",
                  color=colour, lw=1.5, label=mode)
    ax.set_xlabel(r"momentum transfer $\xi_q$")
    ax.set_ylabel(r"$|K_q|$")
    ax.set_ylim(1e-18, 1e0)
    ax.set_title(r"the seam gives $1/q$ tails; folding gives none",
                 fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    save_fig(fig, "reach_energy_coupling_windows.png")
    plt.close(fig)


# --------------------------------------------------------------------- #
def main():
    print("Verification for docs/analysis/reach_energy_coupling.md")
    part_a()
    part_b()
    part_c()
    part_d()
    y_grid_ef, collected = part_ef()
    part_g()
    y_grid_h, table, _ = part_h()
    banner("Figures")
    fig_moments(y_grid_ef, collected)
    fig_budget(y_grid_h, table)
    fig_windows()
    banner("Summary")
    print("  E1  the reach is a period in the separation coordinate, not an")
    print("      aperture; ring, periodic V and horizon are one mechanism")
    print("  E2  folding the tail back in gives the exact kernel of the")
    print("      periodised potential, at the cost of Proposition O4")
    print("  E3  M_0 = M_2 = 0 for any even window: worlds and energy are")
    print("      conserved independently of the reach")
    print("  E4  focus does no net work, hop carries the classical power,")
    print("      for any number of commensurate modes")
    print("  E5  a profile on the full kernel contaminates the third moment")
    print("      with the classical force at order 1 / y_max^2")
    print("  E6  on the compensated residual the contamination vanishes")
    print("      identically: the reach is inert at order hbar^2")
    print("  E7  for a polynomial V the Moyal series terminates and there is")
    print("      no jump measure on the open line; the reach is constitutive")
    print("  E8  the residual budget always diverges with the reach, so the")
    print("      reach buys momentum resolution and no accuracy")
    print("\ndone.")


if __name__ == "__main__":
    main()
