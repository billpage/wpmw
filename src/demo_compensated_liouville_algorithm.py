"""
Verification for ``docs/algorithm/compensated_liouville_algorithm.md``.

The analysis note ``docs/analysis/compensated_liouville_splitting.md``
establishes that the Wigner potential operator splits exactly as

    M(x, s) = M_cl(x, s) + M_res(x, s),
    M_cl  = i V'(x) s,
    M_res = (i / hbar) [ V(x+y) - V(x-y) - 2 y V'(x) ],    y = hbar s / 2,

that the two factors commute, and that at bounded reach ``y_max`` the residual
is a signed jump kernel carrying no net momentum.  This script checks the
things a *specification* has to get right that the note does not address:
what the grid means, what the discrete split does, and whether the reordering
actually buys anything.

Parts
-----
A  Spec section 2.  The reach IS the momentum grid.  ``dp = pi hbar/(2 y_max)``
   is forced, and the number of momentum cells ``N_p`` is the number of
   ket-bra rungs inside the coherence horizon.  Neither is a convergence
   parameter.
B  Spec section 4.1, the Nyquist rung.  With an even ``N_p`` the rung grid
   holds ``-y_max`` but not ``+y_max``, so an odd symbol is sampled
   asymmetrically.  Left alone this makes the kernel complex and drives a
   real ``W`` complex at percent level.  Zeroing that one rung fixes it.
C  Spec section 3.2, discrete compensation.  Subtracting ``i V'(x) s`` leaves
   a spurious force of up to a few percent in the hop channel.  Subtracting
   the kernel's own first moment (Theorem O2) makes the residual exactly
   momentum-neutral on the grid in use.
D  Spec section 4.2.  The residual channel is the crystal-lattice mediated
   transfer rule of ``phase_space_crystal_lattice_algorithm.md`` section 3b
   with the rate field replaced.  Checked against the spectral action, and
   against the analytic ``Gamma_q`` at maximal reach.
E  Spec section 4.4, the event budget.  The total hop rate is not a fixed
   number: it grows like ``(2/pi) |M_res(x, s_max)| ln N_p``, because a hard
   coherence horizon is a hard-edged aperture and has 1/xi sidelobes.
F  Spec section 5.  The reordering reduces the splitting error by the symbol
   ratio ``|M_res|/|M_cl| ~ u^2/6``, uniformly in the step size.
G  Spec section 6.3, the quiet region.  Theorem C7 as an active-set
   optimisation: the fraction of position cells that need the hop channel
   at all.

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
V_P = 1.5
K_COS = 2 * np.pi / L


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
def V_cos(z):
    return -V_P * np.cos(K_COS * z)


def dV_cos(z):
    return V_P * K_COS * np.sin(K_COS * z)


BAR_A, BAR_S = 1.0, 0.5


def V_bar(z):
    return BAR_A * np.exp(-z ** 2 / (2 * BAR_S ** 2))


def dV_bar(z):
    return -z / BAR_S ** 2 * V_bar(z)


COU_Z, COU_EPS = 1.0, 0.5


def V_cou(z):
    return -COU_Z / np.sqrt(z ** 2 + COU_EPS ** 2)


def dV_cou(z):
    return COU_Z * z / (z ** 2 + COU_EPS ** 2) ** 1.5


BUMP_W = 1.0


def bump(z):
    t = np.atleast_1d(np.asarray(z, dtype=float)) / BUMP_W
    out = np.zeros_like(t)
    m = np.abs(t) < 1
    out[m] = np.exp(-1.0 / (1.0 - t[m] ** 2))
    return out if np.ndim(z) else float(out[0])


def dbump(z, h=1e-6):
    return (bump(z + h) - bump(z - h)) / (2 * h)


def V_trap_bump(z):
    return 0.5 * MASS * z ** 2 + bump(z)


def dV_trap_bump(z):
    return MASS * z + dbump(z)


OMEGA_HAR = 1.7


def V_har(z):
    return 0.5 * MASS * OMEGA_HAR ** 2 * z ** 2


def dV_har(z):
    return MASS * OMEGA_HAR ** 2 * z


# --------------------------------------------------------------------- #
# The coherence-horizon occupancy profile (spec 4.1)                     #
# --------------------------------------------------------------------- #
def hann(t):
    """Raised-cosine occupancy profile.  Even, w(0)=1, w(+-1)=w'(+-1)=0."""
    return np.cos(np.pi * t / 2) ** 2


def hard(t):
    """The hard cutoff of Definition (H): full occupancy out to the horizon."""
    return np.ones_like(t)


# --------------------------------------------------------------------- #
# The specified grid and rate field                                      #
# --------------------------------------------------------------------- #
class Grid:
    """The grid of spec section 2.  ``y_max`` and ``N_p`` are the inputs."""

    def __init__(self, y_max, N_p, M_x, box=L):
        if N_p % 2:
            raise ValueError("N_p must be even.")
        self.y_max = float(y_max)
        self.N_p = int(N_p)
        self.M_x = int(M_x)
        self.dp = np.pi * HBAR / (2.0 * self.y_max)          # forced
        self.dx = box / self.M_x
        self.x = -box / 2 + self.dx * np.arange(self.M_x)
        self.p = (np.arange(self.N_p) - self.N_p // 2) * self.dp
        self.s = 2 * np.pi * np.fft.fftfreq(self.N_p, d=self.dp)
        self.y = HBAR * self.s / 2.0                          # the rungs
        self.nyq = self.N_p // 2
        self.q = np.round(np.fft.fftfreq(self.N_p, d=1.0) * self.N_p).astype(int)
        self.xi = self.q * self.dp                            # transfer atoms
        self.kx = 2 * np.pi * np.fft.fftfreq(self.M_x, d=self.dx)

    def symbols(self, Vf, dVf, nyquist_zero=True, discrete_force=True,
                horizon=None, where="residual"):
        """Return (M, M_cl, M_res, accel) sampled on the rung grid.

        ``horizon`` is the occupancy profile w(y/y_max) of the coherence
        horizon: ``None`` for the hard cutoff, else an even callable with
        w(0) = 1 and w(+-1) = w'(+-1) = 0.  ``where`` decides what it
        multiplies -- ``"residual"`` (spec 4.1, the specified placement) or
        ``"symbol"`` (the naive placement, kept only so Part H can show it
        breaking harmonic exactness).

        ``accel`` is the deterministic acceleration of step 1: minus the
        kernel's own first moment (spec 3.2) when ``discrete_force``, else the
        naive ``-V'(x)``.
        """
        X, Y = np.meshgrid(self.x, self.y, indexing="ij")
        _, S = np.meshgrid(self.x, self.s, indexing="ij")
        M = (1j / HBAR) * (Vf(X + Y) - Vf(X - Y))
        S = S.copy()
        if nyquist_zero:
            M[:, self.nyq] = 0.0
            S[:, self.nyq] = 0.0
        w = None if horizon is None else horizon(self.y / self.y_max)[None, :]
        if w is not None and where == "symbol":
            M = M * w
        base = 1j * S                                    # symbol of i s
        fm_M = np.real((self.xi[None, :] * np.fft.ifft(M, axis=1)).sum(axis=1))
        fm_b = np.real((self.xi[None, :] * np.fft.ifft(base, axis=1)).sum(axis=1))
        dV_eff = fm_M / fm_b
        dV = dV_eff if discrete_force else dVf(self.x)
        M_cl = 1j * dV[:, None] * S
        M_res = M - M_cl
        if w is not None and where == "residual":
            M_res = M_res * w
        return M, M_cl, M_res, -dV

    def rate_field(self, M_res):
        """Signed transfer rates K_res(x, xi_q), real, odd in q."""
        return np.real(np.fft.ifft(M_res, axis=1))


def moments(g, K):
    return (np.abs(K).sum(axis=1),
            K.sum(axis=1),
            (g.xi[None, :] * K).sum(axis=1))


def seam(Vf, dVf, x, y_max):
    """|M_res| at the edge of the reach -- the aperture edge discontinuity."""
    return np.abs((Vf(x + y_max) - Vf(x - y_max)
                   - 2 * y_max * dVf(x)) / HBAR)


# --------------------------------------------------------------------- #
# Part A                                                                 #
# --------------------------------------------------------------------- #
def part_a():
    banner("Part A -- spec 2: the reach IS the momentum grid")
    print("  Given a coherence horizon L_c = 2 y_max, the momentum cell is")
    print("  forced: dp = pi hbar / (2 y_max).  N_p is then the number of")
    print("  ket-bra rungs inside the horizon, and the rung spacing in the")
    print("  separation coordinate is a = 2 L_c / N_p.")
    print()
    print(f"  {'y_max':>7} {'L_c':>6} {'N_p':>5} {'dp':>10} {'P_max':>9} "
          f"{'rungs |y|<=y_max':>17} {'a = 2L_c/N_p':>13}")
    for y_max in (0.5, 1.0, 2.0):
        for N_p in (16, 64):
            g = Grid(y_max, N_p, 32)
            nr = int((np.abs(g.y) <= y_max + 1e-12).sum())
            print(f"  {y_max:7.2f} {2*y_max:6.2f} {N_p:5d} {g.dp:10.5f} "
                  f"{g.N_p*g.dp/2:9.3f} {nr:17d} {2*(2*y_max)/N_p:13.5f}")
    print()
    print("  Every rung inside the horizon is a grid point and there are")
    print("  exactly N_p of them.  Refining dp does not refine the answer --")
    print("  it lengthens the reach, which is a change of postulate.")


# --------------------------------------------------------------------- #
# Part B                                                                 #
# --------------------------------------------------------------------- #
def part_b():
    banner("Part B -- spec 4.1: the Nyquist rung")
    print("  With N_p even the rung grid holds y = -y_max but not +y_max, so")
    print("  the odd symbol M is sampled asymmetrically and one rung has no")
    print("  partner.  Averaging the two endpoints of an odd function gives")
    print("  zero, so the specified treatment is to zero that rung.")
    print()
    rng = np.random.default_rng(1)
    g = Grid(1.0, 64, 64)
    W = rng.standard_normal((g.M_x, g.N_p))
    print(f"  {'treatment':>16} {'max|Im K|/max|K|':>18} {'oddness in q':>14} "
          f"{'max |Im W| after step':>22}")
    for tag, nz in (("as sampled", False), ("Nyquist zeroed", True)):
        _, _, M_res, _ = g.symbols(V_cos, dV_cos, nyquist_zero=nz)
        Kc = np.fft.ifft(M_res, axis=1)
        m = np.abs(Kc).max()
        odd = np.abs(Kc + Kc[:, (-g.q) % g.N_p]).max() / m
        Wn = np.fft.ifft(np.exp(0.3 * M_res) * np.fft.fft(W, axis=1), axis=1)
        print(f"  {tag:>16} {np.abs(Kc.imag).max()/m:18.3e} {odd:14.3e} "
              f"{np.abs(Wn.imag).max():22.3e}")
    print()
    print("  (Hard horizon, so that the effect is visible on its own.)")
    print("  Untreated, a single rung drives a real W complex at the percent")
    print("  level every substep.  Zeroing it is not cosmetic.")
    print()
    print("  A soft horizon does NOT make the rule redundant.  The profile")
    print("  zeroes the residual at that rung for free, but dV_eff is read")
    print("  off the unwindowed M and stays contaminated:")
    print(f"  {'treatment':>16} {'dV_eff':>12} {'V\u2032(x)':>12} {'error':>11}")
    for tag, nz in (("as sampled", False), ("Nyquist zeroed", True)):
        _, _, _, acc = g.symbols(V_cos, dV_cos, nyquist_zero=nz, horizon=hann)
        i = int(np.argmin(np.abs(g.x - 1.0)))
        print(f"  {tag:>16} {-acc[i]:12.6f} {dV_cos(g.x[i]):12.6f} "
              f"{abs(-acc[i] - dV_cos(g.x[i])):11.3e}")


# --------------------------------------------------------------------- #
# Part C                                                                 #
# --------------------------------------------------------------------- #
def part_c():
    banner("Part C -- spec 3.2: compensate against the kernel, not V'")
    print("  Theorem C3 says the residual delivers no net momentum.  That is")
    print("  a continuum statement.  On the rung grid the kernel's first")
    print("  moment is not V'(x) (Theorem O2 discretised), so subtracting")
    print("  i V'(x) s leaves a spurious force in the hop channel.")
    print()
    out = {}
    for name, Vf, dVf in (("cosine well", V_cos, dV_cos),
                          ("Gaussian barrier", V_bar, dV_bar),
                          ("soft-core Coulomb", V_cou, dV_cou)):
        print(f"  [{name}]   V'(1) = {dVf(1.0):+.6f}")
        print(f"  {'y_max':>7} {'N_p':>5} {'a_disc':>11} {'dV/dx':>11} "
              f"{'naive residue':>15} {'discrete residue':>17}")
        rows = []
        for y_max in (0.25, 1.0, 2.0):
            for N_p in (32, 128, 512):
                g = Grid(y_max, N_p, 8)
                i = int(np.argmin(np.abs(g.x - 1.0)))
                _, _, Rn, fn = g.symbols(Vf, dVf, discrete_force=False,
                                         horizon=hann)
                _, _, Rd, fd = g.symbols(Vf, dVf, discrete_force=True,
                                         horizon=hann)
                _, _, m1n = moments(g, g.rate_field(Rn))
                _, _, m1d = moments(g, g.rate_field(Rd))
                rows.append((y_max, N_p, abs(m1n[i]), abs(m1d[i])))
                print(f"  {y_max:7.2f} {N_p:5d} {-fd[i]:11.6f} "
                      f"{dVf(g.x[i]):11.6f} {m1n[i]:15.3e} {m1d[i]:17.3e}")
        out[name] = rows
        print()
    print("  a_disc -> V'(x) as N_p grows, and the discretely compensated")
    print("  residual is momentum-neutral to machine precision at every N_p.")
    print("  Step 1 and step 2 are then exactly complementary on the grid in")
    print("  use, which is what a specification has to guarantee.")
    return out


# --------------------------------------------------------------------- #
# Part D                                                                 #
# --------------------------------------------------------------------- #
def part_d():
    banner("Part D -- spec 4.2: the same mediated-transfer rule")
    g = Grid(1.0, 64, 128)
    _, _, M_res, _ = g.symbols(V_cos, dV_cos, horizon=hann)
    K = g.rate_field(M_res)
    rng = np.random.default_rng(0)
    W = rng.standard_normal((g.M_x, g.N_p))
    spec = np.real(np.fft.ifft(M_res * np.fft.fft(W, axis=1), axis=1))
    sten = np.zeros_like(W)
    for j in range(g.N_p):
        q = g.q[j]
        if q <= 0:
            continue
        # Gamma_q = -K_q puts this in the orientation of
        # phase_space_crystal_lattice_algorithm.md section 3b.
        sten += (-K[:, j])[:, None] * (np.roll(W, -q, axis=1)
                                       - np.roll(W, q, axis=1))
    print("  Applying   W += dt sum_q Gamma^res_q(x) [ W(p_{n+q}) - W(p_{n-q}) ]")
    print("  with Gamma^res_q(x) = -K_res(x, xi_q) reproduces the spectral")
    print("  action of L_res exactly:")
    print(f"     max |spectral - stencil| = {np.abs(spec-sten).max():.3e}"
          f"   (relative {np.abs(spec-sten).max()/np.abs(spec).max():.3e})")
    print()
    print("  And at maximal reach the UNcompensated rate field collapses to")
    print("  the analytic single-mode Gamma_q of that specification.  Take")
    print("  V(x) = V_1 cos(2 pi x / L + phi_1) with V_1 = 1.5, phi_1 = pi:")
    gm = Grid(L / 2, 64, 128)
    M, _, _, _ = gm.symbols(V_cos, dV_cos)
    Kf = gm.rate_field(M)
    i1 = int(np.where(gm.q == 1)[0][0])
    im1 = int(np.where(gm.q == -1)[0][0])
    Gam_doc = -(V_P / HBAR) * np.sin(K_COS * gm.x + np.pi)
    other = np.abs(np.delete(Kf, [i1, im1], axis=1)).max()
    print(f"     max |-K(q=1) - Gamma_1|            = "
          f"{np.abs(-Kf[:, i1] - Gam_doc).max():.3e}")
    print(f"     max |K(q)| off the two atoms q=+-1 = {other:.3e}")
    print()
    print("  So the compensated algorithm changes the rate field and adds a")
    print("  deterministic step.  The jump rule itself is unchanged.")


# --------------------------------------------------------------------- #
# Part E                                                                 #
# --------------------------------------------------------------------- #
NS_BUDGET = (256, 1024, 4096, 16384)


def _budget(Vf, dVf, y_max, N_p, horizon, x=1.0):
    g = Grid(y_max, N_p, 8)
    i = int(np.argmin(np.abs(g.x - x)))
    _, _, M_res, _ = g.symbols(Vf, dVf, horizon=horizon)
    K = g.rate_field(M_res)[i]
    return g, np.abs(K).sum(), np.abs(g.xi * K).sum(), (g.xi ** 3 * K).sum(), K


def part_e():
    banner("Part E -- spec 4.4: the coherence horizon needs a profile")
    print("  A HARD horizon is a hard-edged aperture in the ket-bra")
    print("  separation.  M_res is odd, so periodising it over the rung")
    print("  window leaves a jump of 2 |M_res(x, s_max)| at the seam, and a")
    print("  jump gives 1/q sidelobes.  Hence")
    print()
    print("     R(x) = sum_q |K_res| ~ (2/pi) |M_res(x, s_max)| ln N_p + O(1)")
    print()
    data = {}
    for name, Vf, dVf in (("cosine well", V_cos, dV_cos),
                          ("Gaussian barrier", V_bar, dV_bar)):
        print(f"  [{name}, hard horizon]")
        print(f"  {'y_max':>7} " + "".join(f"{'R(N=%d)' % n:>12}"
                                           for n in NS_BUDGET)
              + f"{'fit slope':>12}{'(2/pi)|seam|':>14}")
        for y_max in (0.25, 1.0, 2.0):
            Rs = [_budget(Vf, dVf, y_max, n, None)[1] for n in NS_BUDGET]
            slope = (Rs[-1] - Rs[-2]) / np.log(4)
            pred = 2 * seam(Vf, dVf, 1.0, y_max) / np.pi
            data[(name, y_max)] = (Rs, slope, pred)
            print(f"  {y_max:7.2f} " + "".join(f"{r:12.5f}" for r in Rs)
                  + f"{slope:12.5f}{pred:14.5f}")
        print()
    print("  Predicted and fitted slopes agree to four figures.  So under a")
    print("  hard horizon the event rate is not a number at all: it is a")
    print("  function of the rung count, and it has no limit.  Worse, the")
    print("  momentum churn sum |xi K_res| diverges LINEARLY in N_p, and the")
    print("  third moment -- the leading Moyal coefficient -- does not exist:")
    print()
    print(f"  {'horizon':>16} {'quantity':>16} "
          + "".join(f"{'N=%d' % n:>13}" for n in NS_BUDGET)
          + f"{'tail slope':>12}")
    soft = {}
    for tag, w in (("hard", None), ("raised cosine", hann)):
        R, C, T3 = [], [], []
        for n in NS_BUDGET:
            g, r, c, t3, K = _budget(V_cos, dV_cos, 1.0, n, w)
            R.append(r)
            C.append(c)
            T3.append(t3)
        m = (g.q > 30) & (g.q < 3000)
        slope = np.polyfit(np.log(g.q[m]),
                           np.log(np.abs(K[m]) + 1e-300), 1)[0]
        soft[tag] = (R, C, T3, slope)
        print(f"  {tag:>16} {'R = sum |K|':>16} "
              + "".join(f"{v:13.5f}" for v in R) + f"{slope:12.2f}")
        print(f"  {'':>16} {'sum |xi K|':>16} "
              + "".join(f"{v:13.5f}" for v in C))
        print(f"  {'':>16} {'sum xi^3 K':>16} "
              + "".join(f"{v:13.5f}" for v in T3))
    print()
    print("  The soft horizon converges in all three.  Its third moment,")
    print(f"  {soft['raised cosine'][2][-1]:.6f}, is the analytic leading Moyal")
    ana = -(2 * (-V_P) / HBAR) * np.sin(K_COS * 1.0) * (-(K_COS * HBAR / 2) ** 3)
    print(f"  coefficient {ana:.6f} -- so a soft horizon is what gives the")
    print("  discrete operator a semiclassical expansion at all.")
    print()
    print("  Requirements on the profile w(y/y_max):")
    print("    even             -> both moment conditions of C3 survive exactly")
    print("    w(0) = 1         -> the leading hbar^2 Moyal term is untouched")
    print("    w(+-1) = 0       -> no seam, so the rate converges")
    print("    w'(+-1) = 0      -> 1/q^3 tail, so the momentum churn converges")
    print("    supported on the reach -> the quiet region of C7 survives")
    print("  The raised cosine cos^2(pi t / 2) meets all five and is the")
    print("  specified default.  The hard cutoff of Definition (H) is the")
    print("  w = 1 special case, and is the one choice that fails the last")
    print("  three.")
    return data


# --------------------------------------------------------------------- #
# Part F                                                                 #
# --------------------------------------------------------------------- #
def _integrators(g, M, M_cl, M_res):
    def pot(W, sym, tau):
        return np.real(np.fft.ifft(np.exp(tau * sym) * np.fft.fft(W, axis=1),
                                   axis=1))

    def adv(W, tau):
        ph = np.exp(-1j * g.kx[:, None] * g.p[None, :] * tau / MASS)
        return np.real(np.fft.ifft(ph * np.fft.fft(W, axis=0), axis=0))

    def uncompensated(W, dt):
        return adv(pot(adv(W, dt / 2), M, dt), dt / 2)

    def classical_flow(W, dt, nsub):
        h = dt / nsub
        W = adv(W, h / 2)
        for i in range(nsub):
            W = pot(W, M_cl, h)
            W = adv(W, h if i < nsub - 1 else h / 2)
        return W

    def compensated(W, dt, nsub=24):
        return classical_flow(
            pot(classical_flow(W, dt / 2, nsub), M_res, dt), dt / 2, nsub)

    return uncompensated, compensated


def gaussian(g, x0, p0, sx):
    sp = HBAR / (2 * sx)
    X, P = np.meshgrid(g.x, g.p, indexing="ij")
    W = np.exp(-((X - x0) ** 2) / (2 * sx ** 2)
               - ((P - p0) ** 2) / (2 * sp ** 2))
    return W / (W.sum() * g.dx * g.dp)


def part_f():
    banner("Part F -- spec 5: what the reordering buys")
    g = Grid(1.0, 64, 256)
    u = K_COS * g.y_max
    print(f"  cosine well, reach y_max = {g.y_max} "
          f"(u = k y_max = {u/np.pi:.2f} pi),")
    print(f"  dp = {g.dp:.4f}, N_p = {g.N_p}, sigma_x = 0.25.  The cosine is")
    print("  periodic on the window, so spectral advection is exact here and")
    print("  the measurement is not contaminated by the flow integrator; the")
    print("  algorithm itself needs no box (Part I).")
    print()
    print("  BOTH schemes split the SAME operator M_cl + M_res.  The")
    print("  uncompensated one applies it in a single potential substep; the")
    print("  compensated one puts M_cl in a classical FLOW and leaves M_res")
    print("  alone in the middle.  Strang-splitting streaming against the")
    print("  classical force instead of flowing it would reinstate the very")
    print("  error term the reordering removes; here the flow is sub-cycled.")
    W0 = gaussian(g, 1.0, 0.0, 0.25)
    T = 1.0
    out = {}
    for tag, w in (("hard horizon", None), ("raised cosine", hann)):
        _, M_cl, M_res, _ = g.symbols(V_cos, dV_cos, horizon=w)
        unc, com = _integrators(g, M_cl + M_res, M_cl, M_res)
        ref = W0.copy()
        for _ in range(1024):
            ref = com(ref, T / 1024, nsub=8)
        errs = []
        print()
        print(f"  [{tag}]")
        print(f"  {'steps':>7} {'dt':>9} {'L1 err uncomp':>15} "
              f"{'L1 err comp':>14} {'gain':>7}")
        for ns in (8, 16, 32, 64, 128):
            dt = T / ns
            Wu, Wc = W0.copy(), W0.copy()
            for _ in range(ns):
                Wu = unc(Wu, dt)
                Wc = com(Wc, dt)
            eu = np.abs(Wu - ref).sum() * g.dx * g.dp
            ec = np.abs(Wc - ref).sum() * g.dx * g.dp
            errs.append((dt, eu, ec))
            print(f"  {ns:7d} {dt:9.4f} {eu:15.4e} {ec:14.4e} {eu/ec:7.1f}")
        ratio = np.abs(M_res).max() / np.abs(M_cl).max()
        print(f"  symbol ratio |M_res|/|M_cl| = {ratio:.4f}  "
              f"(1/ratio = {1/ratio:.1f})")
        out[tag] = (errs, ratio)
    print()
    print(f"  Unwindowed, the leading estimate of C4 is u^2/6 = {u**2/6:.4f},")
    print("  matching the hard-horizon ratio.  The soft horizon attenuates")
    print("  the far rungs, where M_res is largest, so its residual is")
    print("  smaller again.  For the hard horizon the gain is flat in dt at")
    print("  the symbol ratio, which is the signature of an unchanged order")
    print("  with a smaller constant.  For the soft horizon the compensated")
    print("  error is close to the reference's own floor, so those gains are")
    print("  lower bounds and their scatter is measurement noise rather than")
    print("  structure.")
    return out["hard horizon"]


# --------------------------------------------------------------------- #
# Part G                                                                 #
# --------------------------------------------------------------------- #
def part_g():
    banner("Part G -- spec 6.3: the quiet region as an active set")
    print("  Theorem C7: if V''' vanishes on [x - y_max, x + y_max] then a")
    print("  world at x takes no events at all.  As an implementation this is")
    print("  an active set -- the hop channel is skipped wherever the rate")
    print("  field is below tolerance, and the crystal-lattice algorithm has")
    print("  no such saving because its rate field is a global Fourier sum.")
    print()
    print("  V = harmonic trap + a C-infinity bump supported on [-1, 1],")
    print("  x in [-6, 6], tolerance 1e-12 on max_q |K_res|.")
    print()
    print(f"  {'y_max':>7} {'predicted quiet':>18} {'active cells':>14} "
          f"{'fraction':>10} {'measured edge':>15} {'(meas-pred)/dx':>16}")
    out = []
    for y_max in (0.25, 0.5, 1.0, 2.0):
        g = Grid(y_max, 64, 480, box=12.0)
        _, _, M_res, _ = g.symbols(V_trap_bump, dV_trap_bump,
                                   horizon=hann)
        amp = np.abs(g.rate_field(M_res)).max(axis=1)
        active = amp > 1e-12
        edge = BUMP_W + y_max
        meas = np.abs(g.x[active]).max()
        out.append((y_max, g.x, amp, edge))
        print(f"  {y_max:7.2f} {'|x| > %.2f' % edge:>18} "
              f"{int(active.sum()):14d} {active.mean():10.4f} "
              f"{meas:15.4f} {(meas - edge)/g.dx:16.1f}")
    print()
    print(f"  Position cell dx = {g.dx:.4f}.  The active set falls silent two")
    print("  to four cells INSIDE the predicted edge -- the bump's third")
    print("  derivative approaches zero smoothly, so the rate drops below")
    print("  tolerance slightly before the geometric edge.  The error is in")
    print("  the safe direction only because C7 is an exact bound on the")
    print("  support: an implementation must mask on the measured rate field")
    print("  and not on the analytic edge, since for a potential without")
    print("  compact support (spec 6.3) there is no exact edge at all.")
    print("  The edge tracks b + y_max with unit slope, as C7 asserts.")
    print()
    print("  A short reach is cheap twice over: fewer events per active")
    print("  world, and fewer active worlds.")
    return out


# --------------------------------------------------------------------- #
# Part H                                                                 #
# --------------------------------------------------------------------- #
def part_h():
    banner("Part H -- spec 4.1: where the horizon profile multiplies")
    print("  Two placements both remove the seam:")
    print("     (B)  M^w = w * M              -> step 1 becomes w * M_cl")
    print("     (C)  M^w = M_cl + w * M_res   -> step 1 untouched")
    print()
    print("  Validation rung 7 discriminates them.  For an open-line")
    print("  quadratic V the residual is identically zero and the algorithm")
    print("  must reduce to exact Newtonian flow.  V = m w^2 x^2 / 2 with")
    print(f"  w = {OMEGA_HAR}, so V'(1) = {dV_har(1.0):.6f}:")
    print()
    print(f"  {'placement':>28} {'max |M_res|':>13} {'dV_eff':>11} "
          f"{'error in the force':>20}")
    for tag, w, where in (("hard horizon", None, "residual"),
                          ("(B) window the symbol", hann, "symbol"),
                          ("(C) window the residual", hann, "residual")):
        g = Grid(1.0, 128, 64)
        i = int(np.argmin(np.abs(g.x - 1.0)))
        _, _, M_res, acc = g.symbols(V_har, dV_har, horizon=w, where=where)
        print(f"  {tag:>28} {np.abs(M_res).max():13.3e} {-acc[i]:11.6f} "
              f"{abs(-acc[i] - dV_har(g.x[i])):20.3e}")
    print()
    print("  (B) manufactures a residual out of nothing in a potential that")
    print("  has none, and shifts the classical force by two per cent.  The")
    print("  reason is structural: w(s) * i V'(x) s is not the symbol of a")
    print("  drift, so windowing the whole symbol takes the Newtonian step")
    print("  apart.  (C) is therefore normative.")
    print()
    print("  Read physically, (C) says the horizon grades the COHERENCE")
    print("  channel only.  The classical force is not a coherence effect --")
    print("  it is the local, first-moment part of the kernel -- so there is")
    print("  nothing there for a rung occupancy to attenuate.")


# --------------------------------------------------------------------- #
# Part I -- the open line                                                #
# --------------------------------------------------------------------- #
def part_i():
    banner("Part I -- spec 2.5: the algorithm is native to the open line")
    print("  Nothing in the construction needs a position box.  The only")
    print("  transform is in the rung direction, which is the coherence")
    print("  ladder, not a boundary condition; x enters as a spectator.")
    print("  V = harmonic trap + a C-infinity bump on [-1, 1], no wrap.")
    print()
    y_max, N_p = 2.0, 64
    g = Grid(y_max, N_p, 600, box=15.0)
    _, _, M_res, acc = g.symbols(V_trap_bump, dV_trap_bump, horizon=hann)
    K = g.rate_field(M_res)
    R = np.abs(K).sum(axis=1)
    active = R > 1e-12
    print(f"  window x in [{g.x[0]:.2f}, {g.x[-1]:.2f}], dx = {g.dx:.4f}, "
          f"y_max = {y_max}, N_p = {N_p}")
    print(f"  active cells {int(active.sum())} of {g.M_x} "
          f"(|x| < {np.abs(g.x[active]).max():.3f}; C7 predicts "
          f"{BUMP_W + y_max:.2f})")
    print(f"  max |dV_eff - V'(x)| over the window = "
          f"{np.abs(-acc - dV_trap_bump(g.x)).max():.3e}")
    print()
    print("  The active set is compact even though the potential is not, so")
    print("  the hop channel lives on a bounded region of an unbounded line.")
    print("  Widening the window adds only quiet cells.")
    print()
    print("  Now run worlds on that line: velocity Verlet plus sampled hops,")
    print("  no periodic wrap, sign carried as a weight (spec 4.3).")
    rng = np.random.default_rng(7)
    N_w, dt, nsteps = 40000, 0.002, 3200
    xw = 2.5 + 0.30 * rng.standard_normal(N_w)
    pw = HBAR / (2 * 0.30) * rng.standard_normal(N_w)
    cdf = np.cumsum(np.abs(K), axis=1)
    tot = cdf[:, -1:].copy()
    tot[tot == 0] = 1.0
    cdf = cdf / tot
    sgn = np.sign(K)
    n_ev, sum_xi, sum_sgn_xi, logw = 0, 0.0, 0.0, 0.0
    a = np.interp(xw, g.x, -acc)
    for _ in range(nsteps):
        pw = pw + 0.5 * dt * a
        xw = xw + dt * pw / MASS
        np.clip(xw, g.x[0], g.x[-1], out=xw)
        a = np.interp(xw, g.x, -acc)
        pw = pw + 0.5 * dt * a
        cell = np.clip(np.searchsorted(g.x, xw) - 1, 0, g.M_x - 1)
        rate = R[cell]
        fire = rng.random(N_w) < rate * dt
        logw += float(np.mean(rate)) * dt
        if not fire.any():
            continue
        idx = np.nonzero(fire)[0]
        c = cell[idx]
        u = rng.random(idx.size)
        j = (cdf[c] > u[:, None]).argmax(axis=1)
        xi_draw = g.xi[j]
        s_draw = sgn[c, j]
        pw[idx] += xi_draw
        n_ev += idx.size
        sum_xi += float(np.abs(xi_draw).sum())
        sum_sgn_xi += float((s_draw * xi_draw).sum())
    print(f"    worlds {N_w}, steps {nsteps}, t = {nsteps*dt:.2f} "
          f"({nsteps*dt/(2*np.pi):.2f} trap periods)")
    print(f"    events drawn                        = {n_ev}")
    print(f"    mean |transfer| per event           = {sum_xi/n_ev:.5f}")
    print(f"    signed mean transfer per event      = "
          f"{sum_sgn_xi/n_ev:+.5f}")
    print(f"    ... as a fraction of the mean scale = "
          f"{abs(sum_sgn_xi/sum_xi):.3e}")
    print(f"    statistical floor  1/sqrt(N_events) = {1/np.sqrt(n_ev):.3e}")
    print()
    print("  The signed mean transfer is at the statistical floor: the")
    print("  sampler inherits C3, delivering no net momentum, so the whole")
    print("  force stays in the deterministic step.  This is the world-form")
    print("  check that the mesh parts cannot make.")
    print()
    print(f"    ensemble weight growth exp(int R dt) = {np.exp(logw):.3f}")
    print("  That growth is the sign problem, priced.  Compensation shrinks")
    print("  R and so shrinks it, but does not remove it -- Proposition T3")
    print("  is untouched by anything in this specification.")


# --------------------------------------------------------------------- #
# Figures                                                                #
# --------------------------------------------------------------------- #
def fig_grid_and_budget(edata):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))

    ax = axes[0]
    Ns = np.array(NS_BUDGET, dtype=float)
    for y_max, col in ((0.25, "#1D9E75"), (1.0, "#D85A30"), (2.0, "0.35")):
        Rs, slope, pred = edata[("cosine well", y_max)]
        ax.loglog(Ns, Rs, "o-", color=col, lw=1.6,
                  label=f"hard, $y_{{\\max}} = {y_max}$")
        ax.loglog(Ns, Rs[-1] + pred * np.log(Ns / Ns[-1]), "--",
                  color=col, lw=1.0)
        Rw = [_budget(V_cos, dV_cos, y_max, int(n), hann)[1] for n in Ns]
        ax.loglog(Ns, Rw, "s:", color=col, lw=1.4, ms=4, mfc="none",
                  label=f"soft, $y_{{\\max}} = {y_max}$")
    ax.set_xlabel("rung count $N_p$")
    ax.set_ylabel(r"hop rate $R(x) = \sum_q |K_{\rm res}|$")
    ax.set_title(r"Hard: $(2/\pi)|M_{\rm res}(s_{\max})|\ln N_p$. Soft: flat.",
                 fontsize=10)
    ax.legend(fontsize=7, ncol=2)
    ax.grid(alpha=0.3, which="both")

    ax = axes[1]
    g = Grid(1.0, 4096, 8)
    i = int(np.argmin(np.abs(g.x - 1.0)))
    _, _, M_hard, _ = g.symbols(V_cos, dV_cos, horizon=None)
    _, _, M_soft, _ = g.symbols(V_cos, dV_cos, horizon=hann)
    pos = g.q > 0
    ax.loglog(g.q[pos], np.abs(g.rate_field(M_hard)[i][pos]), ".",
              color="#D85A30", ms=3, label="hard horizon")
    ax.loglog(g.q[pos], np.abs(g.rate_field(M_soft)[i][pos]), ".",
              color="#1D9E75", ms=3, label="raised-cosine horizon")
    qq = np.array([3.0, 2000.0])
    ax.loglog(qq, seam(V_cos, dV_cos, 1.0, 1.0) / (np.pi * qq), "--",
              color="0.4", lw=1.2, label=r"$|M_{\rm res}(s_{\max})|/\pi q$")
    Ks = np.abs(g.rate_field(M_soft)[i])
    anchor = Ks[int(np.where(g.q == 50)[0][0])] * 50.0 ** 3
    ax.loglog(qq, anchor * qq ** -3.0, "-.", color="0.4", lw=1.2,
              label=r"$\propto q^{-3}$")
    ax.set_ylim(1e-14, 1e-1)
    ax.set_xlabel("transfer index $q$")
    ax.set_ylabel(r"$|K_{\rm res}(x, q\,\Delta p)|$")
    ax.set_title("A hard aperture has $1/q$ sidelobes", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    save_fig(fig, "compensated_algorithm_budget.png")
    plt.close(fig)


def fig_validation(errs, ratio, gdata):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))

    ax = axes[0]
    dts = np.array([e[0] for e in errs])
    eu = np.array([e[1] for e in errs])
    ec = np.array([e[2] for e in errs])
    ax.loglog(dts, eu, "o-", color="0.35", lw=1.6, label="uncompensated")
    ax.loglog(dts, ec, "o-", color="#D85A30", lw=1.6, label="compensated")
    ax.loglog(dts, eu * ratio, "--", color="#1D9E75", lw=1.2,
              label=r"uncompensated $\times\ |M_{\rm res}|/|M_{\rm cl}|$")
    ax.set_xlabel(r"step $\Delta t$")
    ax.set_ylabel(r"$L^1$ error at $t = 1$")
    ax.set_title("The reordering shrinks the Strang constant", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")

    ax = axes[1]
    for (y_max, xg, amp, edge), col in zip(gdata[1:],
                                           ("#1D9E75", "#D85A30", "0.35")):
        ax.semilogy(xg, np.maximum(amp, 1e-20), "-", color=col, lw=1.4,
                    label=f"$y_{{\\max}} = {y_max}$")
        ax.axvline(edge, color=col, ls=":", lw=1.0)
    ax.axvspan(-BUMP_W, BUMP_W, color="#EEEDFE")
    ax.axhline(1e-12, color="0.6", lw=0.8)
    ax.set_xlim(-1, 6)
    ax.set_ylim(1e-20, 1e2)
    ax.set_xlabel("$x$")
    ax.set_ylabel(r"$\max_q |K_{\rm res}(x, \xi_q)|$")
    ax.set_title(r"Active set: quiet beyond $b + y_{\max}$ (bump shaded)",
                 fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    save_fig(fig, "compensated_algorithm_validation.png")
    plt.close(fig)


# --------------------------------------------------------------------- #
def main():
    print("Verification for docs/algorithm/compensated_liouville_algorithm.md")
    part_a()
    part_b()
    part_c()
    part_d()
    edata = part_e()
    errs, ratio = part_f()
    gdata = part_g()
    part_h()
    part_i()
    banner("Figures")
    fig_grid_and_budget(edata)
    fig_validation(errs, ratio, gdata)
    banner("Summary")
    print("  spec 2    reach and momentum cell are one parameter;")
    print("            N_p is the rung count, not a tolerance")
    print("  spec 3.2  compensate against the kernel's first moment, not V';")
    print("            otherwise the hop channel keeps a few percent of force")
    print("  spec 4.1  zero the Nyquist rung, or a real W goes complex")
    print("  spec 4.2  the jump rule is unchanged from the crystal lattice;")
    print("            only the rate field and the extra step 1 are new")
    print("  spec 4.4  a hard horizon has a divergent hop rate and no third")
    print("            moment; the horizon needs a profile, and the profile")
    print("            multiplies the residual only (Part H)")
    print("  spec 5    the reordering shrinks the splitting error by the")
    print("            symbol ratio, uniformly in dt")
    print("  spec 6.3  the quiet region is an active set: the hop channel is")
    print("            skipped outside one reach of the non-quadratic part")
    print("  spec 2.5  no box anywhere; the world sampler on the open line")
    print("            inherits C3's momentum neutrality (Part I)")
    print("\ndone.")


if __name__ == "__main__":
    main()
