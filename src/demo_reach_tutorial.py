#!/usr/bin/env python3
"""
Companion demo for ``docs/supplement/what_the_reach_is.md``.

Every table and every figure in that tutorial is produced here.  The
tutorial is about one quantity -- the reach ``y_max`` -- and the four
questions a reader asks about it in order:

  A  What is separated by ``2y``?  Not two bodies.  Two legs of one
     ket-bra pair, which is to say two worlds' answers about one body.
     Measured: the self-conjugate population carries exactly zero flux,
     so the object that looks most like a classical particle cannot move.
  B  Is the reach a coherence length?  Nearly, and the difference matters.
     An aperture destroys coherence; the periodisation the algorithm
     actually performs aliases it.  Measured on a cat state.
  C  Aperture or period?  The discriminating measurement is commensuration
     against taper on a crystal potential.
  D  How many ceilings are there?  Two, and they are not the same number.

Run::

    WPMW_OUTPUT=/tmp/out PYTHONPATH=src python3 -u src/demo_reach_tutorial.py
"""

from __future__ import annotations

import numpy as np

from wpmwlib.wpmw_utils import docs_path, output_path

HBAR = 1.0
MU_MASS = 1.0

# the cat state used in Parts A and B
CAT_D = 4.0
CAT_SIG = 0.5


def banner(title: str) -> None:
    print("\n" + "=" * 72)
    print("  " + title)
    print("=" * 72)


# --------------------------------------------------------------------- #
# Potentials                                                            #
# --------------------------------------------------------------------- #
def V_sech(z, V0=1.0, a=1.0):
    return V0 / np.cosh(z / a) ** 2


def dV_sech(z, V0=1.0, a=1.0):
    return -2.0 * V0 / (a * np.cosh(z / a) ** 2) * np.tanh(z / a)


def R_sech(x, a=1.0):
    """Analyticity radius: double poles at z = i pi a (n + 1/2)."""
    return np.sqrt(x * x + (np.pi * a / 2.0) ** 2)


def V_coulomb(z):
    return 1.0 / np.where(np.abs(z) < 1e-12, 1e-12, z)


def dV_coulomb(z):
    return -1.0 / np.where(np.abs(z) < 1e-12, 1e-12, z) ** 2


def V_soft(z, eps=0.3):
    return -1.0 / np.sqrt(z * z + eps * eps)


def dV_soft(z, eps=0.3):
    return z * (z * z + eps * eps) ** -1.5


def d3_fd(f, z, h=1e-4):
    return (f(z + 2 * h) - 2 * f(z + h) + 2 * f(z - h) - f(z - 2 * h)) / (2 * h**3)


# --------------------------------------------------------------------- #
# The residual kernel on the momentum lattice                           #
# --------------------------------------------------------------------- #
def residual_kernel(Vfun, dVfun, x, y_max, n_y=1 << 15, taper=True):
    y = -y_max + 2.0 * y_max * np.arange(n_y) / n_y
    w = np.cos(0.5 * np.pi * y / y_max) ** 2 if taper else 1.0
    d = (Vfun(x + y) - Vfun(x - y) - 2.0 * y * dVfun(x)) * w
    q = np.fft.fftfreq(n_y, d=1.0 / n_y).astype(int)
    c = np.fft.fft(d) / n_y * ((-1.0) ** q)
    K = np.real(-1j * c / HBAR)
    return q, np.pi * HBAR / (2.0 * y_max), K


def kernel_moments(q, dp, K):
    xi = q * dp
    return K.sum(), (xi * K).sum(), (xi**3 * K).sum(), np.abs(K).sum()


# --------------------------------------------------------------------- #
# The cat state, as a density matrix in midpoint and half separation    #
# --------------------------------------------------------------------- #
def cat_psi(z, d=CAT_D, sig=CAT_SIG):
    n = (2.0 * np.pi * sig**2) ** -0.25 / np.sqrt(2.0)
    return n * (np.exp(-((z - d / 2) ** 2) / (4 * sig**2))
                + np.exp(-((z + d / 2) ** 2) / (4 * sig**2)))


def cat_rho(x, y, d=CAT_D, sig=CAT_SIG):
    """rho(x+y, x-y).  Real, for this state."""
    return cat_psi(x + y, d, sig) * cat_psi(x - y, d, sig)


# --------------------------------------------------------------------- #
# Part A                                                                #
# --------------------------------------------------------------------- #
def part_a():
    banner("A  what is separated by 2y: legs, not bodies")
    print("  A pair with legs at X and X' has a place at each leg and a")
    print("  relative phase mu = arg rho(X, X').  The position-pair ladder's")
    print("  current between neighbouring legs is J |rho_1| sin(mu) / hbar:")
    print("  it is carried by the MISALIGNMENT, not by anything's velocity.")
    print("  A self-conjugate pair has its legs coincident, so mu = 0.\n")
    a, J = 1.0, 1.0
    print("        mu/pi     sin mu     current    interpretation")
    for f in (0.0, 0.25, 0.5, 0.75, 1.0):
        mu = f * np.pi
        cur = J * np.sin(mu) / (HBAR * a)
        tag = ("self-conjugate: motionless" if f == 0.0 else
               "zone boundary, maximal" if f == 0.5 else
               "antipodal: motionless again" if f == 1.0 else "")
        print(f"     {f:8.2f} {np.sin(mu):10.5f} {cur:10.5f}    {tag}")
    print("\n  A self-conjugate carrier -- legs together, mu = 0 -- carries")
    print("  the full probability density and no current.  It is STILL.")
    print("  Momentum in this representation is p_bar = hbar mu / a, a")
    print("  property of the RELATION between two legs, not something a body")
    print("  carries.  Nothing here separates two bodies; there is one body.")
    print("\n  Exactly, from rho(x+y, x-y) = psi(x+y) psi*(x-y):")
    print("     d/dy rho |_(y=0) = 2i Im(psi* psi') = (2 i m / hbar) j(x)")
    print("  so  j(x) = -(i hbar / 2m) d_y rho |_(y=0) = int (p/m) W dp.")
    print("  The current is the first y-derivative of the density matrix at")
    print("  zero separation.  Checked against the exact current:\n")
    a_g, n_g = 0.004, 200001
    X = (np.arange(n_g) - n_g // 2) * a_g
    print("        p0     from d_y rho     exact j      ratio")
    for p0 in (0.5, 1.0, 2.0, 5.0):
        psi = np.exp(-X**2 / 4.0) * np.exp(1j * p0 * X / HBAR)
        psi /= np.sqrt(np.sum(np.abs(psi) ** 2) * a_g)
        i = n_g // 2
        # y-derivative of rho(x+y, x-y) at y = 0, centred difference in y
        h = a_g
        r_p = psi[i + 1] * np.conj(psi[i - 1])
        r_m = psi[i - 1] * np.conj(psi[i + 1])
        dy_rho = (r_p - r_m) / (2.0 * h)
        j_from = float(np.real(-1j * HBAR / (2.0 * MU_MASS) * dy_rho))
        j_ex = float(np.abs(psi[i]) ** 2) * p0 / MU_MASS
        print(f"    {p0:7.2f} {j_from:14.6f} {j_ex:12.6f} {j_from/j_ex:10.6f}")


def part_a2():
    banner("A2  two quantities are called mu, and they are not the same one")
    print("  (1) LADDER mu = arg rho(X, X'): between the two LEGS of a")
    print("      ket-bra pair -- one body, two branches.")
    print("      d mu / dX = p / hbar, the TOTAL momentum.\n")
    a_l = 0.05
    X = (np.arange(801) - 400) * a_l
    print("          p0     mu/a measured     p0/hbar      ratio")
    for p0 in (0.5, 1.0, 2.0, 5.0):
        psi = np.exp(-X**2 / 4.0) * np.exp(1j * p0 * X / HBAR)
        r1 = psi[1:] * np.conj(psi[:-1])
        mu = float(np.angle(r1[len(r1) // 2]))
        print(f"     {p0:8.2f} {mu / a_l:14.6f} {p0 / HBAR:12.6f}"
              f" {mu / a_l / (p0 / HBAR):11.6f}")
    print("\n  (2) PAIR mu = Phi_a - Phi_b: between the two MEMBERS of a sea")
    print("      pair -- two world-particles, two clocks.")
    print("      d mu / dx = (p_a - p_b) / hbar, the SPLITTING.")
    print("      Pair amplitude |Psi| = 2 |sin(mu/2)|, evaluated at x = 0.7")
    print("      with a relative clock offset theta_a - theta_b = 0.4.\n")
    print("      p_a    p_b      dp   d mu/dx     mu(0.7)    |Psi|    state")
    for pa, pb in ((1.0, 1.0), (1.0, 0.8), (2.0, 1.0), (3.0, -1.0)):
        x = 0.7
        mu = 0.4 + (pa - pb) * x / HBAR
        amp = 2.0 * abs(np.sin(mu / 2.0))
        tag = "co-moving: mu frozen" if pa == pb else "split: mu winds"
        print(f"   {pa:7.2f} {pb:7.2f} {pa-pb:7.2f} {(pa-pb)/HBAR:9.4f}"
              f" {mu:11.4f} {amp:8.4f}    {tag}")
    print("\n  Same structure -- a gauge-invariant relative phase of a")
    print("  two-ended object, whose gradient is the momentum conjugate to")
    print("  the ends' separation.  Different variable.")
    print("\n  Three ways of contributing nothing, and they are distinct:")
    print("    self-conjugate carrier   mu = 0 in (1): no momentum, no")
    print("                             current, full density.  STILL.")
    print("    dark sea pair            mu = 0 in (2): definite p, streams")
    print("                             under (S), |Psi| = 0.  INVISIBLE.")
    print("    free positon or negaton  definite p, streams, contributes")
    print("                             +-1 to E.  Neither.")


# --------------------------------------------------------------------- #
# Part B                                                                #
# --------------------------------------------------------------------- #
def part_b():
    banner("B  is the reach a coherence length?  aperture against fold")
    nx = 1601
    x = np.linspace(-8.0, 8.0, nx)
    dx = x[1] - x[0]
    ny = 3201
    yf = np.linspace(-8.0, 8.0, ny)
    dy = yf[1] - yf[0]
    R0 = cat_rho(x[:, None], yf[None, :])
    tr = float(np.sum(cat_rho(x, np.zeros_like(x))) * dx)
    pur0 = 2.0 * float(np.sum(R0**2)) * dx * dy
    i0 = int(np.argmin(np.abs(x)))
    p = np.linspace(-4.0, 4.0, 1201)

    print(f"  Cat state: two Gaussians at +- d/2 = +- {CAT_D/2:.1f},"
          f" width sigma = {CAT_SIG}.")
    print(f"  Its coherences live at x ~ 0, y ~ +- {CAT_D/2:.1f}.")
    print(f"  Reference:  Tr rho = {tr:.6f},  Tr rho^2 = {pur0:.6f}\n")

    print("  (i) APERTURE -- delete rho for |y| > y_max.  The textbook toy")
    print("      model of decoherence.\n")
    print("      y_max   y_max/(d/2)    Tr rho^2    fringe amplitude at x=0")
    ap = []
    for ym in (0.5, 1.0, 1.5, 2.0, 3.0, 6.0):
        Ra = np.where((np.abs(yf) <= ym)[None, :], R0, 0.0)
        pur = 2.0 * float(np.sum(Ra**2)) * dx * dy
        W = (np.cos(2.0 * p[:, None] * yf[None, :] / HBAR) @ Ra[i0]) \
            * dy / (np.pi * HBAR)
        ap.append((ym, pur, float(W.max() - W.min())))
        print(f"   {ym:9.2f} {ym / (CAT_D/2):11.3f} {pur:12.6f}"
              f"      {W.max() - W.min():14.6f}")

    print("\n  (ii) FOLD -- identify y with y + L_c, L_c = 2 y_max.  What the")
    print("       algorithm's DFT actually does.  Nothing is deleted.\n")
    print("      y_max      L_c     Tr rho^2    fringe amplitude at x=0"
          "       dp")
    fo = []
    for ym in (0.5, 1.0, 1.5, 2.0, 3.0, 6.0):
        Lc = 2.0 * ym
        yb = np.linspace(-ym, ym, 2001)[:-1]
        dyb = yb[1] - yb[0]
        nim = int(np.ceil(12.0 / Lc)) + 1
        Rb = np.zeros((nx, yb.size))
        for k in range(-nim, nim + 1):
            Rb += cat_rho(x[:, None], (yb + k * Lc)[None, :])
        pur = 2.0 * float(np.sum(Rb**2)) * dx * dyb
        dp = np.pi * HBAR / (2.0 * ym)
        qq = np.arange(-60, 61)
        W = np.array([np.sum(Rb[i0] * np.cos(2 * v * dp * yb / HBAR))
                      * dyb / (np.pi * HBAR) for v in qq])
        fo.append((ym, pur, float(W.max() - W.min()), dp))
        print(f"   {ym:9.2f} {Lc:8.2f} {pur:12.6f}"
              f"      {W.max() - W.min():14.6f} {dp:9.4f}")

    print("\n  (iii) WHERE THE EXCESS COMES FROM.  Folding identifies rung k")
    print("        with k + L_c, so a coherence at y = d/2 lands exactly on")
    print("        the diagonal when L_c divides d/2.\n")
    print("      y_max      L_c    true diagonal peak   folded peak    ratio")
    dg = []
    for ym in (0.5, 1.0, 1.5, 2.0, 4.0):
        Lc = 2.0 * ym
        true_d = cat_rho(x, 0.0)
        fold = sum(cat_rho(x, k * Lc) for k in range(-12, 13))
        dg.append((ym, Lc, float(true_d.max()), float(fold.max())))
        print(f"   {ym:9.2f} {Lc:8.2f}        {true_d.max():12.6f}"
              f"   {fold.max():12.6f}  {fold.max() / true_d.max():7.3f}")
    print("\n  d/2 = 2, so the contamination is exact at L_c = 1 and L_c = 2")
    print("  and absent at L_c = 3.  It is commensuration, not magnitude:")
    print("  a LONGER reach (L_c = 3) is clean where a shorter one is not.")
    print("\n  Verdict.  An aperture destroys coherence -- purity falls, the")
    print("  fringe dies.  The fold destroys none: it makes coherences whose")
    print("  separations differ by L_c indistinguishable, and at commensurate")
    print("  separations it books them as probability density.  The reach is")
    print("  a resolution limit on the relation between legs, not a decay.")
    return ap, fo, dg


# --------------------------------------------------------------------- #
# Part C                                                                #
# --------------------------------------------------------------------- #
def part_c():
    banner("C  aperture or period?  the discriminating measurement")
    print("  For V = cos(2 pi x / a) the exact residual budget is a pair of")
    print("  deltas.  An aperture picture says a taper should help, by")
    print("  suppressing ringing.  A period picture says only commensuration")
    print("  matters, and that a taper attenuates the true delta.\n")
    a, x = 4.0, 0.7

    def Vc(z):
        return np.cos(2.0 * np.pi * z / a)

    def dVc(z):
        return -2.0 * np.pi / a * np.sin(2.0 * np.pi * z / a)

    # exact on-lattice budget of the FULL kernel: two deltas of height |V_q|
    exact = 2.0 * abs(np.sin(2.0 * np.pi * x / a)) * 1.0
    print("      L_c    L_c/a   commensurate     sharp      tapered")
    rows = []
    for Lc in (4.0, 8.0, 6.0, 5.0):
        ym = Lc / 2.0
        out = []
        for tap in (False, True):
            y = -ym + 2.0 * ym * np.arange(1 << 14) / (1 << 14)
            w = np.cos(0.5 * np.pi * y / ym) ** 2 if tap else 1.0
            d = (Vc(x + y) - Vc(x - y)) * w
            c = np.fft.fft(d) / (1 << 14)
            out.append(float(np.abs(np.imag(c)).sum()))
        com = "yes" if abs(Lc / a - round(Lc / a)) < 1e-12 else "no "
        rows.append((Lc, Lc / a, com, out[0] / exact, out[1] / exact))
        print(f"   {Lc:6.1f} {Lc/a:8.2f}      {com}      {out[0]/exact:9.4f}"
              f"  {out[1]/exact:9.4f}")
    print("\n  (ratios to the exact budget).  The period picture wins: sharp")
    print("  and commensurate is exact; sharp and incommensurate is wild;")
    print("  and for a crystal the taper is strictly WORSE than a")
    print("  commensurate sharp window, which no aperture picture predicts.")
    return rows


# --------------------------------------------------------------------- #
# Part D                                                                #
# --------------------------------------------------------------------- #
def taylor_coeffs(x, n_max=36, radius=1.2):
    n_th = 4096
    th = 2.0 * np.pi * np.arange(n_th) / n_th
    z = radius * np.exp(1j * th)
    f = V_sech(x + z) - V_sech(x - z) - 2.0 * z * dV_sech(x)
    c = np.fft.fft(f) / n_th
    return np.array([c[n] / radius**n for n in range(n_max + 1)])


def part_d():
    banner("D  two ceilings, and they are not the same number")
    x = 0.3
    R = R_sech(x)
    exact3 = HBAR**2 / 4.0 * d3_fd(V_sech, x)
    print("  CEILING 1, existence.  The kernel evaluates V on the REAL")
    print("  segment [x - y_max, x + y_max].  A singularity there kills it.")
    print("  CEILING 2, convergence.  R(x) is the distance to the nearest")
    print("  COMPLEX singularity; past it the Moyal series in y diverges.\n")
    print(f"  sech^2, a = 1, x = {x}:  R(x) = {R:.6f}")
    print(f"  exact third moment hbar^2 V'''(x)/4 = {exact3:.6f}\n")
    print("     y_max   y_max/R       M0          M1         M3      M3 err"
          "    budget")
    rows = []
    for ym in (0.4, 0.8, 1.2, 1.5992, 2.4, 3.2, 6.4, 12.8):
        q, dp, K = residual_kernel(V_sech, dV_sech, x, ym)
        m0, m1, m3, bud = kernel_moments(q, dp, K)
        rows.append((ym, ym / R, m3, abs(m3 - exact3), bud, dp))
        print(f"  {ym:8.3f} {ym / R:8.3f} {m0:11.2e} {m1:11.2e}"
              f" {m3:10.6f} {abs(m3-exact3):9.2e} {bud:9.4f}")
    print("\n  Nothing in that table knows where R is.  World number and")
    print("  momentum stay exact (C3), the third moment stays on")
    print("  hbar^2 V'''/4 (E6) -- and is an order of magnitude BETTER at")
    print("  four times R than inside it, because more of the kernel's")
    print("  Fourier content is resolved.  Only the budget moves (E8).\n")
    c = taylor_coeffs(x)
    print("  Meanwhile the series at the same x:\n")
    print("         y      y/R       exact          S_10         S_20"
          "         S_30")
    for y in (0.5, 1.0, 1.55, 1.70, 2.20):
        ex = V_sech(x + y) - V_sech(x - y) - 2.0 * y * dV_sech(x)
        s = [float(np.real(sum(c[n] * y**n for n in range(k + 1))))
             for k in (10, 20, 30)]
        print(f"   {y:8.3f} {y/R:8.3f} {ex:12.4e} {s[0]:12.3e}"
              f" {s[1]:12.3e} {s[2]:12.3e}")
    print("\n  The exact column is smooth and bounded throughout: V is real-")
    print("  analytic on the real line and its poles sit at +- i pi/2.")
    print("  Crossing R costs the semiclassical expansion, not the kernel.\n")
    print("  CONTRAST: a pole ON the real axis.  Pure Coulomb V = 1/z at")
    print("  x = 1, whose arm x - y reaches the origin at y = 1; against the")
    print("  soft core at the same x, whose R = 1.044 but whose poles are at")
    print("  +- i eps.\n")
    print("      y_max   |  Coulomb max|D_res|    budget  |"
          "  soft core max|D_res|   budget")
    for ym in (0.5, 0.9, 0.99, 1.5, 3.0):
        y = np.linspace(-ym, ym, 20001)[1:-1]
        dc = V_coulomb(1.0 + y) - V_coulomb(1.0 - y) - 2.0 * y * dV_coulomb(1.0)
        ds = V_soft(1.0 + y) - V_soft(1.0 - y) - 2.0 * y * dV_soft(1.0)
        _, _, Kc = residual_kernel(V_coulomb, dV_coulomb, 1.0, ym)
        _, _, Ks = residual_kernel(V_soft, dV_soft, 1.0, ym)
        print(f"   {ym:8.3f}   |     {np.abs(dc).max():13.4e}"
              f" {np.abs(Kc).sum():9.3e}  |    {np.abs(ds).max():12.4f}"
              f" {np.abs(Ks).sum():9.4f}")
    print("\n  Ceiling 1 is a ceiling on the model.  Ceiling 2 is a ceiling")
    print("  on the semiclassical reading of the model, and the soft core")
    print("  sails three times past it without noticing.")
    return rows, R, exact3


def part_d2():
    banner("D2  the ceiling in more than one dimension")
    print("  For -Z / sqrt(r.r + eps^2) the complexified singular set is the")
    print("  CONE  r.r = -eps^2, not a point.  In d >= 2 the imaginary part")
    print("  of y can be taken transverse to x, which lets the cone reach")
    print("  closer to a real field point than the real singularity does.\n")
    print("  Minimising |y|^2 = |a|^2 + |b|^2 subject to")
    print("      |x+a|^2 - |b|^2 + eps^2 = 0   and   (x+a).b = 0")
    print("  gives a = -x/2, hence |b|^2 = |x|^2/4 + eps^2 and\n")
    print("      R_(d>=2) = sqrt( |x|^2 / 2 + eps^2 )")
    print("      R_1      = sqrt( |x|^2     + eps^2 )\n")
    print("  independent of d for every d >= 2, because only one transverse")
    print("  direction is used.  One dimension is the outlier, and the loose")
    print("  one.  Checked against a constrained numerical minimisation in")
    print("  d = 3:\n")
    try:
        from scipy.optimize import minimize
    except ImportError:
        print("  (scipy unavailable; analytic values only)")
        minimize = None
    print("       eps      |x|      R numeric      R analytic      R_1D"
          "     ratio")
    for eps in (0.3, 1.0):
        for xm in (0.0, 1.0, 2.0, 5.0):
            pred = np.sqrt(xm**2 / 2.0 + eps**2)
            one = np.sqrt(xm**2 + eps**2)
            if minimize is None:
                num = float("nan")
            else:
                x = np.array([xm, 0.0, 0.0])
                obj = lambda v: v[:3] @ v[:3] + v[3:] @ v[3:]
                c1 = lambda v: ((x + v[:3]) @ (x + v[:3])
                                - v[3:] @ v[3:] + eps**2)
                c2 = lambda v: (x + v[:3]) @ v[3:]
                best = np.inf
                for seed in range(8):
                    rng = np.random.default_rng(seed)
                    v0 = rng.normal(scale=max(xm, eps) + 0.5, size=6)
                    r = minimize(obj, v0, options={"maxiter": 800,
                                                   "ftol": 1e-14},
                                 constraints=[{"type": "eq", "fun": c1},
                                              {"type": "eq", "fun": c2}])
                    if r.success:
                        best = min(best, np.sqrt(r.fun))
                num = best
            print(f"    {eps:7.2f} {xm:8.2f} {num:14.6f} {pred:15.6f}"
                  f" {one:10.6f} {num/one:9.4f}")
    print("\n  So three dimensions makes the Coulomb ceiling TIGHTER, by up")
    print("  to sqrt(2).  Theorem Z4's threshold would tighten by a factor")
    print(f"  (sqrt 2)^4 = 4, to eps >= pi^4 = {np.pi**4:.2f} a0 for one rung.")
    print("  This is ceiling 2, so it binds the semiclassical reading and")
    print("  not the model; see section 5 and open item R-SP1.")


# --------------------------------------------------------------------- #
# Figures                                                               #
# --------------------------------------------------------------------- #
def save(fig, name):
    fig.savefig(output_path(name), dpi=150, bbox_inches="tight")
    dp = docs_path(name)
    if dp:
        fig.savefig(dp, dpi=150, bbox_inches="tight")
    print(f"  wrote {name}")


def figures(ap, fo, dg, d_rows, R, EX3):
    banner("E  figures")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # --- figure 1: the geometry, and the two ceilings ------------------
    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.1))

    ym = 1.6
    ax[0].axhspan(-ym, ym, color="C0", alpha=0.12)
    ax[0].axhline(ym, color="C0", lw=1.4)
    ax[0].axhline(-ym, color="C0", lw=1.4)
    xs = np.linspace(-4, 4, 400)
    ax[0].plot(xs, 0.35 * V_sech(xs) - 0.0, "k-", lw=1.0, alpha=0.5)
    for (xm, yy) in ((-2.0, 1.1), (0.5, 0.6), (2.4, 1.45)):
        ax[0].plot([xm, xm], [-yy, yy], color="C3", lw=1.2)
        ax[0].plot([xm], [yy], "o", color="C3", ms=5)
        ax[0].plot([xm], [-yy], "o", color="C3", ms=5)
        ax[0].plot([xm], [0.0], "x", color="0.3", ms=6)
    ax[0].set_xlabel("midpoint  x")
    ax[0].set_ylabel("half separation  y")
    ax[0].set_title("legs of one pair, not two bodies", fontsize=10)
    ax[0].text(-3.9, ym + 0.12, r"$+y_{max}$", color="C0", fontsize=9)
    ax[0].text(-3.9, -ym - 0.3, r"$-y_{max}$", color="C0", fontsize=9)
    ax[0].set_ylim(-2.3, 2.3)
    ax[0].grid(alpha=0.25)

    th = np.linspace(0, 2 * np.pi, 400)
    ax[1].plot(R * np.cos(th), R * np.sin(th), "k--", lw=1.2,
               label=r"$R(x)$: series ceiling")
    for ymv, col, off in ((0.8, "C2", 0.09), (3.2, "C3", -0.09)):
        ax[1].plot([-ymv, ymv], [off, off], color=col, lw=3.0,
                   solid_capstyle="butt",
                   label=rf"$y_{{max}}={ymv}$ (reach)")
    ax[1].plot([0.0], [np.pi / 2], "kx", ms=9)
    ax[1].plot([0.0], [-np.pi / 2], "kx", ms=9)
    ax[1].text(0.15, np.pi / 2 - 0.22, r"poles of $V$", fontsize=8)
    ax[1].axhline(0, color="0.7", lw=0.7)
    ax[1].set_xlabel(r"Re $y$")
    ax[1].set_ylabel(r"Im $y$")
    ax[1].set_title("the kernel only walks the real axis", fontsize=10)
    ax[1].legend(fontsize=8, loc="lower right")
    ax[1].set_aspect("equal")
    ax[1].grid(alpha=0.25)

    fig.tight_layout()
    save(fig, "reach_tutorial_geometry.png")
    plt.close(fig)

    # --- figure 2: aperture against fold -------------------------------
    fig2, ax2 = plt.subplots(1, 2, figsize=(11.5, 3.9))
    a_y = np.array([r[0] for r in ap])
    ax2[0].plot(a_y, [r[1] for r in ap], "o-", color="C3",
                label=r"aperture: Tr $\rho^2$")
    ax2[0].plot([r[0] for r in fo], [r[1] for r in fo], "s-", color="C0",
                label=r"fold: Tr $\rho^2$")
    ax2[0].axhline(1.0, color="k", ls="--", lw=1.0)
    ax2[0].axvline(CAT_D / 4, color="0.6", ls=":", lw=1.0)
    ax2[0].set_xlabel(r"$y_{max}$")
    ax2[0].set_title("(a) purity: decay against aliasing", fontsize=10)
    ax2[0].legend(fontsize=8)
    ax2[0].grid(alpha=0.3)

    ax2[1].plot(a_y, [r[2] for r in ap], "o-", color="C3",
                label="aperture: fringe")
    ax2[1].plot([r[0] for r in fo], [r[2] for r in fo], "s-", color="C0",
                label="fold: fringe")
    ax2[1].axvline(CAT_D / 2, color="0.6", ls=":", lw=1.0)
    ax2[1].set_xlabel(r"$y_{max}$")
    ax2[1].set_title("(b) the interference term", fontsize=10)
    ax2[1].legend(fontsize=8)
    ax2[1].grid(alpha=0.3)

    fig2.tight_layout()
    save(fig2, "reach_tutorial_coherence.png")
    plt.close(fig2)

    # --- figure 3: one dial, three consequences ------------------------
    fig3, ax3 = plt.subplots(figsize=(7.4, 4.3))
    yv = np.array([r[0] for r in d_rows])
    ax3.loglog(yv, [r[5] for r in d_rows], "o-", color="C0",
               label=r"$\Delta p = \pi\hbar/2y_{max}$  (buys)")
    ax3.loglog(yv, [r[4] for r in d_rows], "s-", color="C3",
               label=r"event budget $\sum_q|K_q|$  (costs)")
    ax3.loglog(yv, [r[2] / EX3 for r in d_rows], "^-", color="C2",
               label=r"$M_3\,/\,(\hbar^2 V'''/4)$  (buys no accuracy)")
    ax3.axvline(R, color="k", ls="--", lw=1.1)
    ax3.text(R * 1.06, 0.14, r"$R(x)$", fontsize=9)
    ax3.set_xlabel(r"$y_{max}$")
    ax3.set_title("one dial, three consequences", fontsize=11)
    ax3.legend(fontsize=8, loc="lower left")
    ax3.grid(alpha=0.3, which="both")
    fig3.tight_layout()
    save(fig3, "reach_tutorial_tradeoff.png")
    plt.close(fig3)


def main():
    part_a()
    part_a2()
    ap, fo, dg = part_b()
    part_c()
    d_rows, R, EX3 = part_d()
    part_d2()
    figures(ap, fo, dg, d_rows, R, EX3)
    banner("summary")
    print("  2y separates two LEGS of one pair -- two worlds' answers about")
    print("  one body -- and not two bodies.")
    print("  The reach is a period in that separation, not an aperture;")
    print("  it aliases coherence rather than destroying it;")
    print("  it has two ceilings, one on the model and one on the")
    print("  semiclassical reading of it, the second tighter in 3D;")
    print("  and the two quantities called mu are not the same one.")


if __name__ == "__main__":
    main()
