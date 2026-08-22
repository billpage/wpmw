"""
Verification for ``docs/analysis/compensated_liouville_splitting.md``.

The Wigner potential operator, written in the variable ``s`` conjugate to
momentum, is multiplication by the symbol

    M(x, s) = (i / hbar) [ V(x + y) - V(x - y) ],    y = hbar s / 2

where ``y`` is the half ket-bra separation.  Subtracting the part linear in
``y`` leaves

    M_cl (x, s) = i V'(x) s
    M_res(x, s) = (i / hbar) [ V(x+y) - V(x-y) - 2 y V'(x) ]

and ``M_res`` is exactly the odd part of the cubic Taylor remainder of ``V``
about ``x``.  The question this script settles is whether that is a genuine
classical / quantum split: is ``M_cl`` the whole classical force, and is
``M_res`` realisable as an interaction that carries no force of its own?

The answer is yes, and the condition is a bounded coherence reach.

Parts
-----
A  Theorems C1, C2.  The factorisation is exact and the factors commute; the
   residual is the cubic Taylor remainder and is bounded by
   (2 / hbar) (y_max^3 / 6) sup |V'''| over the reach.
B  Theorem C3, the reach theorem.  Restricted to a reach ``y_max``, the
   residual kernel has zero zeroth moment (worlds are conserved) and zero
   first moment (no net momentum), so it is a bounded signed jump measure --
   a focus-and-hop that delivers no force.  The whole classical force stays
   in the deterministic step.
C  Theorem C4, the reach condition.  The event budget ratio
   TV(K_res) / TV(K) runs from 0.006 to 16 as the reach grows, crossing unity
   near k y_max = pi / 2.  The reach and the momentum quantum are the same
   parameter: dp = pi hbar / (2 y_max).
D  Theorem C5, the ring.  On a circle V''' = 0 forces V constant; the
   periodised parabola's residual support is the bowtie
   |x| + |y| > L/2; and the ring pins every world at maximal reach, where
   u = q pi and the split cannot help.  A coherence horizon fixes both.
E  Evolution in the ring parabola against a cosine well, for scale.
F  Theorem C6, Coulomb.
G  Summary.

Run with ``WPMW_OUTPUT`` set (``/mnt/user-data/outputs`` in the container).
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from wpmwlib.wpmw_utils import docs_path, output_path

HBAR = 1.0
MASS = 1.0
L = 8.0
MX = 256
NP = 128
DX = L / MX
DP = np.pi * HBAR / L
OMEGA = 2.0
V_P = 1.5
K_COS = 2 * np.pi / L

x_grid = -L / 2 + DX * np.arange(MX)
p_grid = (np.arange(NP) - NP // 2) * DP
s_grid = 2 * np.pi * np.fft.fftfreq(NP, d=DP)
S2, X2 = np.meshgrid(s_grid, x_grid, indexing="ij")


def banner(text):
    print()
    print("=" * 74)
    print(text)
    print("=" * 74)


def save_fig(fig, name):
    fig.savefig(output_path(name), dpi=150, bbox_inches="tight")
    dp = docs_path(name)
    if dp:
        fig.savefig(dp, dpi=150, bbox_inches="tight")
    print(f"    [figure] {name}")


# --------------------------------------------------------------------- #
# Potentials                                                             #
# --------------------------------------------------------------------- #
def wrap(z):
    return (z + L / 2) % L - L / 2


def V_par_open(z):
    return 0.5 * MASS * OMEGA ** 2 * z ** 2


def dV_par_open(z):
    return MASS * OMEGA ** 2 * z


def V_par_ring(z):
    return V_par_open(wrap(z))


def dV_par_ring(z):
    return dV_par_open(wrap(z))


def V_cos(z):
    return -V_P * np.cos(K_COS * z)


def dV_cos(z):
    return V_P * K_COS * np.sin(K_COS * z)


def d3V_cos(z):
    return -V_P * K_COS ** 3 * np.sin(K_COS * z)


BAR_A, BAR_S = 1.0, 0.5


def V_bar(z):
    return BAR_A * np.exp(-z ** 2 / (2 * BAR_S ** 2))


def dV_bar(z):
    return -z / BAR_S ** 2 * V_bar(z)


BUMP_A, BUMP_W = 1.0, 1.0


def bump(z):
    """C-infinity bump supported exactly on [-BUMP_W, BUMP_W]."""
    t = np.atleast_1d(np.asarray(z, dtype=float)) / BUMP_W
    out = np.zeros_like(t)
    m = np.abs(t) < 1
    out[m] = BUMP_A * np.exp(-1.0 / (1.0 - t[m] ** 2))
    return out if np.ndim(z) else float(out[0])


def dbump(z, h=1e-6):
    return (bump(z + h) - bump(z - h)) / (2 * h)


def V_trap_bump(z):
    return 0.5 * MASS * 1.0 ** 2 * z ** 2 + bump(z)


def dV_trap_bump(z):
    return MASS * 1.0 ** 2 * z + dbump(z)


COU_Z, COU_EPS = 1.0, 0.5


def V_cou(z):
    return -COU_Z / np.sqrt(z ** 2 + COU_EPS ** 2)


def dV_cou(z):
    return COU_Z * z / (z ** 2 + COU_EPS ** 2) ** 1.5


# --------------------------------------------------------------------- #
# Symbols, kernels, substeps                                             #
# --------------------------------------------------------------------- #
def symbols(Vf, dVf, X=X2, S=S2):
    a = HBAR * S / 2.0
    M = (1j / HBAR) * (Vf(X + a) - Vf(X - a))
    Mcl = 1j * dVf(X) * S
    return M, Mcl, M - Mcl


def reach_kernels(Vf, dVf, x, y_max, n=4096):
    """Kernels in momentum transfer xi, for a world of coherence reach y_max.

    Returns (xi, K, K_res, max|M|, max|M_res|).
    """
    s_max = 2 * y_max / HBAR
    s = np.linspace(-s_max, s_max, n, endpoint=False)
    ds = s[1] - s[0]
    y = HBAR * s / 2
    M = (1j / HBAR) * (Vf(x + y) - Vf(x - y))
    Mcl = 1j * dVf(x) * s
    Mres = M - Mcl
    xi = 2 * np.pi * np.fft.fftfreq(n, d=ds)
    K = np.fft.ifft(np.fft.ifftshift(M))
    Kres = np.fft.ifft(np.fft.ifftshift(Mres))
    return xi, K, Kres, np.abs(M).max(), np.abs(Mres).max()


def moments(xi, K):
    return np.abs(K).sum(), np.real(K.sum()), np.real((xi * K).sum())


def potential_substep(W, Msym, tau):
    Wh = np.fft.fft(W, axis=0)
    Wh *= np.exp(tau * Msym)
    return np.real(np.fft.ifft(Wh, axis=0))


def advect_half(W, dt):
    kx = 2 * np.pi * np.fft.fftfreq(MX, d=DX)
    Wh = np.fft.fft(W, axis=1)
    Wh *= np.exp(-1j * kx[None, :] * p_grid[:, None] * dt / (2 * MASS))
    return np.real(np.fft.ifft(Wh, axis=1))


def strang(W, Msym, dt):
    return advect_half(potential_substep(advect_half(W, dt), Msym, dt), dt)


def gaussian(x0, p0, sigma_x):
    sigma_p = HBAR / (2 * sigma_x)
    P, X = np.meshgrid(p_grid, x_grid, indexing="ij")
    W = np.exp(-((X - x0) ** 2) / (2 * sigma_x ** 2)
               - ((P - p0) ** 2) / (2 * sigma_p ** 2))
    return W / (W.sum() * DX * DP)


# --------------------------------------------------------------------- #
# Part A                                                                 #
# --------------------------------------------------------------------- #
def part_a():
    banner("Part A -- C1 exact commuting split, C2 the Taylor remainder")
    M, Mcl, Mres = symbols(V_par_open, dV_par_open)
    print(f"  max |M - (Mcl + Mres)|      = {np.abs(M - (Mcl + Mres)).max():.3e}")
    print(f"  max |Mres|, open parabola   = {np.abs(Mres).max():.3e}"
          f"   (V''' = 0, so C2 gives 0)")
    Mc, Mclc, Mrc = symbols(V_cos, dV_cos)
    for tau in (0.01, 0.1, 1.0):
        err = np.abs(np.exp(tau * Mc)
                     - np.exp(tau * Mclc) * np.exp(tau * Mrc)).max()
        print(f"  tau = {tau:4.2f}: |e^(tau M) - e^(tau Mcl) e^(tau Mres)| "
              f"= {err:.3e}   (C1)")
    print(f"  max |Mcl Mres - Mres Mcl|   = "
          f"{np.abs(Mclc * Mrc - Mrc * Mclc).max():.3e}")
    print()
    print("  C2: Mres is the odd part of the cubic Taylor remainder of V, so")
    print("      |Mres| <= (2/hbar)(y_max^3 / 6) sup |V'''| over the reach.")
    print(f"  {'y_max':>7} {'max |Mres|':>13} {'Taylor bound':>14} {'ratio':>8}")
    x = 1.0
    for y_max in (0.25, 0.5, 1.0, 2.0):
        yy = np.linspace(-y_max, y_max, 2001)
        got = np.abs((1j / HBAR)
                     * (V_cos(x + yy) - V_cos(x - yy) - 2 * yy * dV_cos(x))).max()
        bound = (2 / HBAR) * (y_max ** 3 / 6) * np.abs(d3V_cos(
            np.linspace(x - y_max, x + y_max, 2001))).max()
        print(f"  {y_max:7.2f} {got:13.6f} {bound:14.6f} {got / bound:8.4f}")


# --------------------------------------------------------------------- #
# Part B -- the reach theorem                                            #
# --------------------------------------------------------------------- #
def part_b():
    banner("Part B -- C3, the reach theorem")
    print("  A world of coherence reach y_max only ever consults V on")
    print("  [x - y_max, x + y_max], so the symbol is only needed for")
    print("  |hbar s / 2| <= y_max.  On that window:")
    print("    - K     is the full Wigner jump kernel; its first moment is")
    print("            the whole classical force -V'(x);")
    print("    - K_res has zero zeroth moment (worlds conserved) and zero")
    print("            first moment (no net momentum), so it carries no force.")
    print()
    out = {}
    for name, Vf, dVf, x in (("cosine well", V_cos, dV_cos, 1.0),
                             ("Gaussian barrier", V_bar, dV_bar, 1.0),
                             ("soft-core Coulomb", V_cou, dV_cou, 1.0)):
        print(f"  [{name}]  -V'(x) = {-dVf(x):+.6f}")
        print(f"  {'y_max':>7} {'sum K':>11} {'mom K':>12} "
              f"{'sum K_res':>12} {'mom K_res':>12}")
        rows = []
        for y_max in (0.25, 0.5, 1.0, 2.0, 4.0):
            xi, K, Kres, mM, mR = reach_kernels(Vf, dVf, x, y_max)
            tvK, z0K, m1K = moments(xi, K)
            tvR, z0R, m1R = moments(xi, Kres)
            rows.append((y_max, tvK, tvR, m1K, m1R, mM, mR))
            print(f"  {y_max:7.2f} {z0K:11.2e} {m1K:12.6f} "
                  f"{z0R:12.2e} {m1R:12.2e}")
        out[name] = rows
        print()
    print("  The residual's first moment is five to six orders below the force")
    print("  it would otherwise double count.  So the split")
    print("     step 1: deterministic acceleration by the full -V'(x)")
    print("     step 2: signed zero-mean hops drawn from K_res")
    print("  is exact and does not apply the force twice.")
    return out


# --------------------------------------------------------------------- #
# Part C -- the reach condition                                          #
# --------------------------------------------------------------------- #
def part_c(bdata):
    banner("Part C -- C4, the reach condition and the crossover")
    print("  Event budget is the kernel's total variation.  Compensation is")
    print("  worth doing exactly when TV(K_res) < TV(K).")
    print()
    for name in bdata:
        print(f"  [{name}]")
        print(f"  {'y_max':>7} {'u/pi = k y/pi':>14} {'TV(K)':>10} "
              f"{'TV(K_res)':>11} {'ratio':>9} {'dp = pi h/2y':>13}")
        for y_max, tvK, tvR, _, _, _, _ in bdata[name]:
            u = K_COS * y_max / np.pi if name == "cosine well" else float("nan")
            print(f"  {y_max:7.2f} {u:14.4f} {tvK:10.4f} {tvR:11.4f} "
                  f"{tvR / tvK:9.4f} {np.pi * HBAR / (2 * y_max):13.4f}")
        print()
    print("  The crossover for the cosine well sits near k y_max = pi / 2, a")
    print("  reach of about a quarter wavelength.  Beyond it the compensated")
    print("  channel is the more expensive one.")
    print()
    print("  Reach and momentum quantum are the same parameter: a reach y_max")
    print("  resolves momentum transfers only down to dp = pi hbar / (2 y_max).")
    print("  At maximal reach on a ring, y_max = L/2 gives dp = pi hbar / L --")
    print("  the crystal quantum -- and the kernel becomes the two-atom mode")
    print("  stencil.  Check for the cosine well at x = 1:")
    xi, K, Kres, _, _ = reach_kernels(V_cos, dV_cos, 1.0, L / 2)
    tv = np.abs(K).sum()
    analytic = 2 * (V_P / HBAR) * abs(np.sin(K_COS * 1.0))
    print(f"     TV(K) at y_max = L/2   = {tv:.6f}")
    print(f"     2 |Gamma_q| analytic   = {analytic:.6f}")


# --------------------------------------------------------------------- #
# Part D -- the ring                                                     #
# --------------------------------------------------------------------- #
def part_d():
    banner("Part D -- C5, the ring")
    M, Mcl, Mres = symbols(V_par_ring, dV_par_ring)
    act = np.abs(Mres).max(axis=0)
    inside = np.abs(HBAR * S2 / 2.0) + np.abs(X2) < L / 2 - 1e-9
    print(f"  max |Mres| inside the bowtie complement     = "
          f"{np.abs(Mres[inside]).max():.3e}")
    print(f"  fraction of the (x, s) grid with |Mres| > 0 = "
          f"{(np.abs(Mres) > 1e-9).mean():.4f}")
    print(f"  predicted (|x| + |y| > L/2)                 = "
          f"{(~inside).mean():.4f}")
    print()
    print("  On a circle V''' = 0 forces V quadratic and a periodic quadratic")
    print("  is constant, so no non-constant ring potential has a vanishing")
    print("  residual.  Worse, the ring pins every world at maximal reach:")
    print(f"  {'mode q':>7} {'u = q pi':>10} {'|Mres|/|Mcl|':>14} {'helps?':>8}")
    for q in (1, 2, 3, 4):
        u = q * np.pi
        r = abs(np.sin(u) - u) / abs(u)
        print(f"  {q:7d} {u:10.4f} {r:14.4f} {'no':>8}")
    print("  At u = q pi the residual is exactly as large as the classical")
    print("  term for every mode, so on a ring the split buys nothing.")
    print()
    print("  Corollary C5.1.  A coherence horizon L_c is the window")
    print("  |y| <= L_c / 2.  It fixes both problems at once: it puts the seam")
    print("  out of reach, and it moves worlds off the maximal-reach point.")
    for Lc in (2.0, 4.0, 6.0):
        win = np.abs(HBAR * S2 / 2.0) <= Lc / 2 + 1e-12
        far = np.abs(X2) < L / 2 - Lc / 2 - 1e-9
        print(f"     L_c = {Lc:.1f}:  max |Mres| for |x| < "
              f"{L / 2 - Lc / 2:.2f}  =  {np.abs(Mres[win & far]).max():.3e}")
    return Mres, act


# --------------------------------------------------------------------- #
# Part E                                                                 #
# --------------------------------------------------------------------- #
def part_e():
    banner("Part E -- evolution in the ring parabola, against a cosine well")
    M_exact, Mcl_exact, _ = symbols(V_par_ring, dV_par_ring)
    M_c, Mcl_c, _ = symbols(V_cos, dV_cos)
    dt, nsteps = 0.01, 200
    W_ref, W_cl = gaussian(1.0, 0.0, 0.4), gaussian(1.0, 0.0, 0.4)
    W_cref, W_ccl = gaussian(1.0, 0.0, 0.4), gaussian(1.0, 0.0, 0.4)
    snaps, times, d_par, d_cos = {}, [], [], []
    for n in range(nsteps + 1):
        t = n * dt
        if n in (0, nsteps // 2, nsteps):
            snaps[t] = (W_ref.copy(), W_cl.copy())
        times.append(t)
        d_par.append(np.abs(W_ref - W_cl).sum() * DX * DP)
        d_cos.append(np.abs(W_cref - W_ccl).sum() * DX * DP)
        if n == nsteps:
            break
        W_ref = strang(W_ref, M_exact, dt)
        W_cl = strang(W_cl, Mcl_exact, dt)
        W_cref = strang(W_cref, M_c, dt)
        W_ccl = strang(W_ccl, Mcl_c, dt)

    def neg(W):
        return -np.minimum(W, 0.0).sum() * DX * DP

    print(f"  {'t':>7s} {'L1(exact - Newtonian)':>23s} {'negativity exact':>18s}")
    for t in sorted(snaps):
        We, _ = snaps[t]
        print(f"  {t:7.2f} {d_par[int(round(t / dt))]:23.3e} {neg(We):18.3e}")
    print()
    print(f"  final L1 gap, periodised parabola = {d_par[-1]:.4e}")
    print(f"  final L1 gap, cosine well         = {d_cos[-1]:.4e}")
    print(f"  ratio                             = {d_cos[-1] / d_par[-1]:.0f}")
    return times, d_par, d_cos, snaps


# --------------------------------------------------------------------- #
# Part F                                                                 #
# --------------------------------------------------------------------- #
def part_f():
    banner("Part F -- C6, Coulomb")
    Z = 1.0
    maxerr = 0.0
    for xv in (0.25, 0.5, 1.0, 2.0, 4.0, 8.0):
        for r in np.linspace(0.05, 0.95, 19):
            a = r * xv
            s = 2 * a / HBAR
            Mn = (1j / HBAR) * (-Z / (xv + a) + Z / (xv - a))
            Mcln = 1j * (Z / xv ** 2) * s
            maxerr = max(maxerr,
                         abs(abs((Mn - Mcln) / Mcln) - r ** 2 / (1 - r ** 2)))
    print("  For V = -Z/x, with rho = y / x:  M_res / M_cl = rho^2/(1 - rho^2)")
    print(f"  max |numeric - rho^2/(1-rho^2)| = {maxerr:.3e}")
    print()
    print("  M / M_cl = 1/(1 - rho^2) = sum_n rho^(2n), so the n-th Moyal term")
    print("  is exactly rho^(2n) times the classical one, all the same sign.")
    print("  The series converges iff |y| < |x|: the Moyal series converges")
    print("  iff the reach does not touch the nucleus.  For Coulomb the reach")
    print("  condition of C4 is therefore not a refinement but a necessity --")
    print("  a world whose reach spans the nucleus has no expansion at all.")
    print()
    print(f"  {'rho = y/x':>11} {'|Mres|/|Mcl|':>14}")
    for r in (0.05, 0.1, 0.25, 0.5, 0.9):
        print(f"  {r:11.2f} {r ** 2 / (1 - r ** 2):14.5f}")
    print()
    print("  In state terms rho ~ hbar / (2 sigma_p |x|), so the condition is")
    print("  sigma_p |x| >> hbar / 2.  Hydrogen's ground state sits at about 1;")
    print("  a Rydberg state of quantum number n sits at about n.")


# --------------------------------------------------------------------- #
# Part G -- the quiet region                                             #
# --------------------------------------------------------------------- #
def res_amplitude(x, y_max, Vf, dVf, n=2048):
    s = np.linspace(-2 * y_max / HBAR, 2 * y_max / HBAR, n)
    y = HBAR * s / 2
    return np.abs((1 / HBAR) * (Vf(x + y) - Vf(x - y) - 2 * y * dVf(x))).max()


def part_g():
    banner("Part G -- C7, the quiet region on the open line")
    print("  V = harmonic trap + a C-infinity bump supported on [-1, 1], so")
    print("  the third derivative of V vanishes outside the bump.  C7 predicts silence for")
    print("  |x| > 1 + y_max.")
    for y_max in (0.5, 1.0, 2.0):
        edge = BUMP_W + y_max
        print()
        print(f"  reach y_max = {y_max}:  predicted quiet for |x| > {edge}")
        print(f"  {'x':>7} {'max |Mres|':>13} {'quiet?':>8}")
        for x in (0.0, edge - 0.3, edge - 0.05, edge + 0.05, edge + 3.0):
            v = res_amplitude(x, y_max, V_trap_bump, dV_trap_bump)
            tag = "yes" if v < 1e-12 else "no"
            note = "  (by parity)" if (x == 0.0 and v < 1e-12) else ""
            print(f"  {x:7.2f} {v:13.3e} {tag:>8}{note}")
    print()
    print("  x = 0 is quiet for a different reason: V is even there, so the")
    print("  odd part of the Taylor remainder cancels by parity even though")
    print("  the third derivative does not vanish on the reach.  C7 is sufficient, not")
    print("  necessary.")
    print()
    print("  A Gaussian barrier has no exactly quiet region, but the same")
    print("  mechanism localises the interaction: the residual tracks the")
    print("  barrier profile translated outward by exactly the reach.")
    sig, y_max = 0.4, 1.0

    def Vg(z):
        return 0.5 * MASS * z ** 2 + np.exp(-z ** 2 / (2 * sig ** 2))

    def dVg(z):
        return MASS * z - z / sig ** 2 * np.exp(-z ** 2 / (2 * sig ** 2))

    print(f"  {'x':>6} {'max |Mres|':>13} {'shifted profile':>17} "
          f"{'bare profile':>14}")
    for x in (1.0, 1.5, 2.0, 2.5, 3.0):
        print(f"  {x:6.2f} {res_amplitude(x, y_max, Vg, dVg):13.3e} "
              f"{np.exp(-(abs(x) - y_max) ** 2 / (2 * sig ** 2)):17.3e} "
              f"{np.exp(-x ** 2 / (2 * sig ** 2)):14.3e}")
    print()
    print("  This is the finite-reach refinement of Theorem O1: the untruncated")
    print("  kernel has no position envelope, but a bounded reach confines the")
    print("  interaction to within y_max of the non-quadratic part of V.")


# --------------------------------------------------------------------- #
# Part H -- the ring as a diagnostic                                     #
# --------------------------------------------------------------------- #
def ring_modes(qmax):
    q = np.arange(1, qmax + 1)
    Vq = 0.5 * MASS * OMEGA ** 2 * L ** 2 / np.pi ** 2 * (-1.0) ** q / q ** 2
    return q, Vq


def part_h():
    banner("Part H -- the ring as a diagnostic")
    q, Vq = ring_modes(4000)
    k = 2 * np.pi * q / L
    x, y = 1.0, 1.0                      # inside the bowtie: |x| + |y| = 2 < 4
    exact = ((1 / HBAR) * (V_par_ring(x + y) - V_par_ring(x - y))
             - dV_par_ring(x) * (2 * y / HBAR))
    terms = -(2 * Vq / HBAR) * np.sin(k * x) * (np.sin(k * y) - k * y)
    print(f"  Inside the bowtie (x = {x}, y = {y}), C6.1 gives exactly zero;")
    print(f"  from V pointwise: {exact:+.10f}")
    print()
    print(f"  {'Q':>6} {'signed sum':>15} {'largest term':>14} "
          f"{'sum |terms|':>13}")
    for Q in (1, 20, 100, 500, 4000):
        print(f"  {Q:6d} {terms[:Q].sum():+15.6f} "
              f"{np.abs(terms[:Q]).max():14.4f} "
              f"{np.abs(terms[:Q]).sum():13.2f}")
    print()
    print("  The signed sum falls like 1/Q; the largest term does not fall at")
    print("  all; the absolute sum diverges logarithmically.  The modes do")
    print("  cancel, conditionally.  Bounding mode by mode gives")
    print("  sum |V_q| k_q^3, which diverges, while C2 on the reach gives zero")
    print("  in one line -- the case for never decomposing V into modes.")
    print()
    print("  At maximal reach y = L/2 the two arms meet at the antipode:")
    print(f"  {'q':>3} {'u/pi':>6} {'sin(u)':>12} {'|Mres|/|Mcl|':>14}")
    for qq in (1, 2, 3, 4):
        kk = 2 * np.pi * qq / L
        u = kk * (L / 2)
        Aq = 0.5 * MASS * OMEGA ** 2 * L ** 2 / np.pi ** 2 * (-1.0) ** qq / qq ** 2
        mr = -(2 * Aq / HBAR) * np.sin(kk * x) * (np.sin(u) - u)
        mc = -(2 * Aq / HBAR) * np.sin(kk * x) * u
        print(f"  {qq:3d} {u / np.pi:6.1f} {np.sin(u):12.3e} "
              f"{abs(mr / mc):14.5f}")
    arms = (1 / HBAR) * (V_par_ring(x + L / 2) - V_par_ring(x - L / 2))
    print(f"  full symbol at y = L/2: {arms:.3e}   (the arms coincide)")
    print("  So a ring pins every mode at the crossover of C4 by geometry.")
    print("  A ring is not a valid testbed for the reach condition.")


# --------------------------------------------------------------------- #
# Part I -- structure of the residual kernel                             #
# --------------------------------------------------------------------- #
def sym_kernel(Vf, dVf, x, y_max, n=2049):
    """Kernel of L_res in the momentum-transfer variable xi.

    n must be odd so that s = 0 is a grid point and ifftshift aligns.
    """
    assert n % 2 == 1
    s_max = 2 * y_max / HBAR
    s = (np.arange(n) - (n - 1) // 2) * (2 * s_max / n)
    y = HBAR * s / 2
    M = (1j / HBAR) * (Vf(x + y) - Vf(x - y))
    Mres = M - 1j * dVf(x) * s
    return np.fft.ifft(np.fft.ifftshift(M)), np.fft.ifft(np.fft.ifftshift(Mres))


def part_i():
    banner("Part I -- what L_res does: structure of K_res")
    print("  L_res W(x, p) = integral K_res(x, xi) W(x, p - xi) d xi,")
    print("  with xi the momentum transferred TO the world.")
    print()
    print("  K_res is real and odd in xi (Mres is imaginary and odd in s):")
    print(f"  {'potential':>18} {'n':>7} {'max|Im K|/max|K|':>18} "
          f"{'oddness':>11}")
    for name, Vf, dVf in (("cosine well", V_cos, dV_cos),
                          ("Gaussian barrier", V_bar, dV_bar)):
        for n in (2049, 8193):
            _, Kr = sym_kernel(Vf, dVf, 1.0, 1.0, n)
            m = np.abs(Kr).max()
            odd = np.abs(Kr + Kr[(-np.arange(n)) % n]).max() / m
            print(f"  {name:>18} {n:7d} {np.abs(Kr.imag).max() / m:18.3e} "
                  f"{odd:11.3e}")
    print()
    print("  Consequences.  Oddness gives zero total rate for free -- no world")
    print("  is created or destroyed.  Zero FIRST moment is the extra property")
    print("  that compensation buys (C3).  And a non-zero odd kernel must take")
    print("  both signs, so L_res is never the generator of a one-body Markov")
    print("  jump process, for any potential: Proposition T3 again.")
    print()
    print(f"  {'potential':>18} {'max K_res':>13} {'min K_res':>13} "
          f"{'sum K_res':>12}")
    for name, Vf, dVf in (("cosine well", V_cos, dV_cos),
                          ("Gaussian barrier", V_bar, dV_bar)):
        _, Kr = sym_kernel(Vf, dVf, 1.0, 1.0)
        r = Kr.real
        print(f"  {name:>18} {r.max():+13.4e} {r.min():+13.4e} "
              f"{r.sum():+12.2e}")
    print()
    print("  Which transfers exist.  M is periodic in s with period 2a/hbar")
    print("  exactly when V is periodic in x with period a, and a symbol")
    print("  periodic in s has its kernel on a lattice in xi of spacing")
    print(f"  pi hbar / a = {np.pi * HBAR / L:.6f} = dp.  Check on the cosine well:")
    x = 1.0
    for s0 in (0.7, 1.9, 3.3):
        a = (1j / HBAR) * (V_cos(x + HBAR * s0 / 2) - V_cos(x - HBAR * s0 / 2))
        s1 = s0 + 2 * L / HBAR
        b = (1j / HBAR) * (V_cos(x + HBAR * s1 / 2) - V_cos(x - HBAR * s1 / 2))
        print(f"     s = {s0:4.1f}:  |M(s) - M(s + 2a/hbar)| = {abs(a - b):.3e}")
    print("  So mode q of V contributes the transfer hbar k_q / 2 = q dp: a")
    print("  world in momentum cell n exchanges with cells n +- q at a rate set")
    print("  by the q-th Fourier amplitude of V at the world's own position.")
    print("  For an aperiodic V there is no lattice and the spectrum of")
    print("  available transfers is continuous.")


# --------------------------------------------------------------------- #
# Figures                                                                #
# --------------------------------------------------------------------- #
def fig_symbol():
    u = np.linspace(-np.pi, np.pi, 600)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    ax = axes[0]
    ax.plot(u / np.pi, np.sin(u), "-", color="0.35", lw=1.6,
            label=r"full  $\sin u$")
    ax.plot(u / np.pi, u, "--", color="#1D9E75", lw=1.6,
            label=r"classical  $u$")
    ax.plot(u / np.pi, np.sin(u) - u, ":", color="#D85A30", lw=2.0,
            label=r"residual  $\sin u - u$")
    ax.set_xlabel(r"$u/\pi$,   $u = k\,y_{\max}$")
    ax.set_ylabel("symbol / prefactor")
    ax.set_title("The split of a single mode symbol", fontsize=10)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    ax = axes[1]
    uu = np.linspace(1e-4, np.pi, 500)
    ax.semilogy(uu / np.pi, np.abs(np.sin(uu) - uu) / uu, "-",
                color="#D85A30", lw=1.8, label=r"$|M_{\rm res}|/|M_{\rm cl}|$")
    ax.semilogy(uu / np.pi, uu ** 2 / 6, "--", color="0.5", lw=1.2,
                label=r"leading $u^2/6$")
    ax.axhline(1.0, color="0.7", lw=0.8)
    ax.axvline(0.5, color="#1D9E75", ls="--", lw=1.0)
    ax.set_xlabel(r"$u/\pi$")
    ax.set_ylabel("ratio")
    ax.set_title("Reach condition: the residual is small only for "
                 "$u \\ll \\pi$", fontsize=10)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    save_fig(fig, "compensated_symbol_and_band.png")
    plt.close(fig)


def fig_reach(bdata):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    ax = axes[0]
    for name, col in (("cosine well", "#D85A30"),
                      ("Gaussian barrier", "#1D9E75"),
                      ("soft-core Coulomb", "0.35")):
        ys = [r[0] for r in bdata[name]]
        rat = [r[2] / r[1] for r in bdata[name]]
        ax.loglog(ys, rat, "o-", color=col, label=name)
    ax.axhline(1.0, color="0.7", lw=1.0)
    ax.set_xlabel(r"coherence reach $y_{\max}$")
    ax.set_ylabel(r"$TV(K_{\rm res})\,/\,TV(K)$")
    ax.set_title("Event budget: gain below the line, loss above", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")

    ax = axes[1]
    for name, col in (("cosine well", "#D85A30"),
                      ("Gaussian barrier", "#1D9E75"),
                      ("soft-core Coulomb", "0.35")):
        ys = [r[0] for r in bdata[name]]
        m1 = [abs(r[4]) for r in bdata[name]]
        f0 = abs(bdata[name][0][3])
        ax.loglog(ys, [m / f0 for m in m1], "o-", color=col, label=name)
    ax.set_xlabel(r"coherence reach $y_{\max}$")
    ax.set_ylabel(r"$|$first moment of $K_{\rm res}|\,/\,|V'(x)|$")
    ax.set_title("The residual carries no force", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    save_fig(fig, "compensated_reach.png")
    plt.close(fig)


def fig_world_paths():
    """One potential, three reaches: Newtonian arcs plus zero-mean hops.

    The hop rate at each x is the measured total variation of K_res; the hop
    size is the cosine mode's own quantum hbar k / 2, with a symmetric sign so
    that the channel carries no net momentum.
    """
    rng = np.random.default_rng(11)
    dt, n = 0.004, 4000
    kick = HBAR * K_COS / 2
    xs_probe = np.linspace(-L / 2, L / 2, 33)
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
    for ax, y_max in zip(axes, (0.5, 1.5, 3.0)):
        tv = np.array([np.abs(reach_kernels(V_cos, dV_cos, xp, y_max,
                                            n=1024)[2]).sum()
                       for xp in xs_probe])
        x, p = 1.6, 0.0
        xs, ps, ev = [], [], []
        for _ in range(n):
            p -= dV_cos(x) * dt
            if rng.random() < float(np.interp(x, xs_probe, tv)) * dt:
                p += kick * rng.choice([-1.0, 1.0])
                ev.append((x, p))
            x += p / MASS * dt
            x = (x + L / 2) % L - L / 2
            xs.append(x)
            ps.append(p)
        ax.plot(xs, ps, "-", color="#888780", lw=0.9)
        if ev:
            ax.plot(*zip(*ev), "o", color="#D85A30", ms=3)
        ax.set_title(f"reach $y_{{\\max}} = {y_max}$,  "
                     f"{len(ev)} events", fontsize=10)
        ax.set_xlabel("$x$")
        ax.set_ylabel("$p$")
        ax.grid(alpha=0.3)
    fig.tight_layout()
    save_fig(fig, "compensated_world_paths.png")
    plt.close(fig)


def fig_quiet_region():
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    ax = axes[0]
    xs = np.linspace(0.0, 6.0, 400)
    for y_max, col in ((0.5, "#1D9E75"), (1.0, "#D85A30"), (2.0, "0.35")):
        vals = [res_amplitude(xv, y_max, V_trap_bump, dV_trap_bump, n=512)
                for xv in xs]
        ax.semilogy(xs, np.maximum(vals, 1e-18), "-", color=col, lw=1.5,
                    label=f"$y_{{\max}} = {y_max}$")
        ax.axvline(BUMP_W + y_max, color=col, ls=":", lw=1.0)
    ax.axvspan(0, BUMP_W, color="#EEEDFE")
    ax.set_xlabel("$x$")
    ax.set_ylabel(r"$\max_y |M_{\rm res}|$")
    ax.set_title("Quiet beyond $|x| = b + y_{\\max}$ (bump shaded)",
                 fontsize=10)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3, which="both")

    ax = axes[1]
    sig, y_max = 0.4, 1.0

    def Vg(z):
        return 0.5 * MASS * z ** 2 + np.exp(-z ** 2 / (2 * sig ** 2))

    def dVg(z):
        return MASS * z - z / sig ** 2 * np.exp(-z ** 2 / (2 * sig ** 2))

    xs = np.linspace(0.8, 3.5, 200)
    vals = [res_amplitude(xv, y_max, Vg, dVg, n=512) for xv in xs]
    ax.semilogy(xs, vals, "-", color="#D85A30", lw=1.8,
                label=r"$\max_y |M_{\rm res}|$")
    ax.semilogy(xs, np.exp(-(np.abs(xs) - y_max) ** 2 / (2 * sig ** 2)), "--",
                color="#1D9E75", lw=1.5, label="barrier shifted by $y_{\\max}$")
    ax.semilogy(xs, np.exp(-xs ** 2 / (2 * sig ** 2)), ":", color="0.5",
                lw=1.5, label="bare barrier")
    ax.set_xlabel("$x$")
    ax.set_ylabel("magnitude")
    ax.set_title("The reach translates the interaction outward", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    save_fig(fig, "compensated_quiet_region.png")
    plt.close(fig)


def fig_ring_residual(Mres, act):
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    order = np.argsort(s_grid)
    Z = np.abs(Mres)[order, :]
    sv = s_grid[order]
    ax = axes[0]
    im = ax.pcolormesh(x_grid, HBAR * sv / 2, Z, cmap="magma", shading="auto")
    ax.plot(x_grid, L / 2 - np.abs(x_grid), "w--", lw=1.0)
    ax.plot(x_grid, -(L / 2 - np.abs(x_grid)), "w--", lw=1.0)
    ax.set_xlabel("$x$")
    ax.set_ylabel(r"reach $y = \hbar s/2$")
    ax.set_title(r"$|M_{\rm res}|$, periodised parabola", fontsize=10)
    fig.colorbar(im, ax=ax)

    _, _, Mres2 = symbols(V_par_open, dV_par_open)
    ax = axes[1]
    im = ax.pcolormesh(x_grid, HBAR * sv / 2, np.abs(Mres2)[order, :],
                       cmap="magma", shading="auto", vmin=0,
                       vmax=max(Z.max(), 1e-12))
    ax.set_xlabel("$x$")
    ax.set_title(r"$|M_{\rm res}|$, open parabola (identically 0)", fontsize=10)
    fig.colorbar(im, ax=ax)

    ax = axes[2]
    ax.plot(x_grid, act / act.max(), "-", color="#D85A30", lw=1.6)
    ax.set_xlabel("$x$")
    ax.set_ylabel(r"$\max_y |M_{\rm res}|$ (normalised)")
    ax.set_title("Residual activity is a seam artefact", fontsize=10)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    save_fig(fig, "compensated_ring_residual.png")
    plt.close(fig)


def fig_evolution(times, d_par, d_cos, snaps):
    ts = sorted(snaps)
    fig, axes = plt.subplots(2, len(ts) + 1, figsize=(4 * (len(ts) + 1), 7))
    for j, t in enumerate(ts):
        We, Wc = snaps[t]
        vmax = np.abs(We).max()
        lv = np.linspace(-vmax, vmax, 21)
        axes[0, j].contourf(x_grid, p_grid, We, levels=lv, cmap="RdBu_r",
                            extend="both")
        axes[0, j].set_title(f"exact QLE, t = {t:.2f}", fontsize=10)
        axes[1, j].contourf(x_grid, p_grid, Wc, levels=lv, cmap="RdBu_r",
                            extend="both")
        axes[1, j].set_title(f"Newtonian only, t = {t:.2f}", fontsize=10)
        for r in (0, 1):
            axes[r, j].set_xlim(-L / 2, L / 2)
            axes[r, j].set_ylim(-4, 4)
            if j == 0:
                axes[r, j].set_ylabel("$p$")
        axes[1, j].set_xlabel("$x$")
    ax = axes[0, -1]
    ax.semilogy(times, d_par, "-", color="#D85A30", label="periodised parabola")
    ax.semilogy(times, d_cos, "-", color="#1D9E75", label="cosine well")
    ax.set_xlabel("$t$")
    ax.set_ylabel(r"$L^1$(exact $-$ Newtonian)")
    ax.set_title("Departure from Newtonian motion", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")
    axes[1, -1].axis("off")
    fig.tight_layout()
    save_fig(fig, "compensated_evolution.png")
    plt.close(fig)


def fig_coulomb():
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    r = np.linspace(0, 0.98, 400)
    ax = axes[0]
    ax.semilogy(r, r ** 2 / (1 - r ** 2), "-", color="#D85A30", lw=1.8)
    ax.axhline(1.0, color="0.7", lw=0.8)
    ax.set_xlabel(r"$\rho = y/x$")
    ax.set_ylabel(r"$|M_{\rm res}|/|M_{\rm cl}|$")
    ax.set_title(r"Coulomb: $\rho^2/(1-\rho^2)$", fontsize=10)
    ax.grid(alpha=0.3, which="both")

    ax = axes[1]
    xv = np.linspace(0.3, 8.0, 300)
    a0 = 0.1
    ax.loglog(xv, 2 * a0 ** 3 / (xv ** 2 * (xv ** 2 - a0 ** 2)), "-",
              color="#D85A30", lw=1.8, label=r"$|M_{\rm res}| \sim 1/x^4$")
    ax.loglog(xv, 2 * a0 / xv ** 2, "--", color="#1D9E75", lw=1.6,
              label=r"$|M_{\rm cl}| \sim 1/x^2$")
    ax.set_xlabel("$x$")
    ax.set_ylabel("symbol magnitude")
    ax.set_title("The hop channel retreats into the core", fontsize=10)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    save_fig(fig, "compensated_coulomb.png")
    plt.close(fig)


# --------------------------------------------------------------------- #
def main():
    print("Verification for docs/analysis/compensated_liouville_splitting.md")
    part_a()
    bdata = part_b()
    part_c(bdata)
    Mres, act = part_d()
    times, d_par, d_cos, snaps = part_e()
    part_f()
    part_g()
    part_h()
    part_i()
    banner("Figures")
    fig_symbol()
    fig_reach(bdata)
    fig_world_paths()
    fig_quiet_region()
    fig_ring_residual(Mres, act)
    fig_evolution(times, d_par, d_cos, snaps)
    fig_coulomb()
    banner("Summary")
    print("  C1  the split is exact and the factors commute")
    print("  C2  the residual is the cubic Taylor remainder of V")
    print("  C3  at bounded reach the residual is a bounded, number- and")
    print("      momentum-conserving signed jump measure: a focus-and-hop")
    print("      that carries no force")
    print("  C4  the split gains when the reach is short compared with the")
    print("      potential's variation scale; reach and momentum quantum are")
    print("      the same parameter")
    print("  C5  the ring pins worlds at maximal reach, where the split buys")
    print("      nothing; a coherence horizon fixes that and the seam at once")
    print("  C6  Coulomb: the Moyal series converges iff the reach misses the")
    print("      nucleus")
    print("  C7  quiet region: a vanishing third derivative on the reach")
    print("      implies no events at x; a bounded reach confines the")
    print("      interaction to within y_max of the non-quadratic part of V")
    print("\ndone.")


if __name__ == "__main__":
    main()
