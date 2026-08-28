"""
Verification for ``docs/analysis/eckart_barrier_compensated.md``.

The Eckart (sech^2) pair barrier on the open line, read through the
compensated Liouville split of ``compensated_liouville_splitting.md``.

Why this potential.  ``reach_energy_coupling.md`` Theorem E7 shows a
polynomial ``V`` has no jump measure on the open line, and
``interworld_coupling.md`` Theorem I4 shows a coupling linear in the leg
separation -- the harmonic and inverted harmonic alike -- has no jump channel
at all.  Between them the project has had no open-line test problem that
exercises the hop channel against a closed form.  ``V0 sech^2(r/a)`` does:
it is bounded, asymptotically free on both sides, has ``V''' != 0``, and has
an exact transmission coefficient.

Conventions follow ``compensated_liouville_splitting.md`` section 1 and
``reach_energy_coupling.md``:

    y           = hbar s / 2                 half ket-bra separation
    D(x, y)     = V(x + y) - V(x - y)
    D_res(x, y) = D(x, y) - 2 y V'(x)        compensated residual
    M           = i D / hbar,   M_cl = i D_cl / hbar,   M_res = i D_res / hbar
    dp          = pi hbar / (2 y_max)
    L_res W(p)  = sum_q K_q W(p - xi_q),     xi_q = q dp

Parts
-----
A  The split for the Eckart barrier.  C2 bound tightness, C3 moments.
B  Theorem K1.  The reach has a ceiling and it is the analyticity strip:
   the Moyal series in y converges iff |y| < dist(0, singularities of
   y -> V(x + y)).  For sech^2 that is sqrt(x^2 + (pi a / 2)^2), uniformly
   at least pi a / 2.  Corollary K1.1: dp > hbar / a, so fewer than beta
   rungs span the barrier's own momentum scale.  Corollary K1.2: for
   soft-core Coulomb the softening length *is* the ceiling.
C  Theorem K2.  Exact far-field law for the residual, the hyperbolic
   continuation of Lemma C0.  Errata for
   ``compensated_liouville_splitting.md`` section 5.1.
D  Theorem K3.  The spectrum-weighted budget ratio is a universal function
   of y_max / a that saturates at 1 from below: on the open line the
   compensated channel is never worse than the uncompensated one.
E  Theorem K4.  Deterministic acceleration conserves the classical
   transmission functional exactly, so the entire quantum-classical gap in
   the transmission is delivered by the residual channel.  Exact spectral
   evolution: full symbol, classical symbol alone, and the compensated
   product.
F  Theorem K5.  The gap is the small imbalance of two large opposed flows
   across the classical separatrix.  Time-resolved in- and out-fluxes.
G  Theorem K6.  Gross traffic saturates in beta while the net falls as
   1 / beta, because the packet is centred on the barrier and the leading
   imbalance is the Jacobian dp/dE across the tunnelling window.
H  Sample world-particle trajectories: Newtonian arcs punctuated by
   positon-negaton emissions straddling the separatrix.

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
MU = 0.5  # reduced mass of two unit-mass particles
RNG = np.random.default_rng(20260827)


def banner(text):
    print()
    print("=" * 72)
    print(text)
    print("=" * 72)


def save_fig(fig, name):
    fig.savefig(output_path(name), dpi=150, bbox_inches="tight")
    dp = docs_path(name)
    if dp:
        fig.savefig(dp, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {name}")


# --------------------------------------------------------------------- #
# The potential                                                          #
# --------------------------------------------------------------------- #
def V(z, v0, a):
    return v0 / np.cosh(z / a) ** 2


def dV(z, v0, a):
    return -2.0 * v0 / a * np.tanh(z / a) / np.cosh(z / a) ** 2


def d3V(z, v0, a):
    t = np.tanh(z / a)
    c = 1.0 / np.cosh(z / a) ** 2
    return (8.0 * v0 / a**3) * c * t * (2.0 * c - t**2)


def beta_of(v0, a):
    """Semiclassical parameter beta = sqrt(2 mu V0) a / hbar."""
    return np.sqrt(2.0 * MU * v0) * a / HBAR


def T_eckart(energy, v0, a):
    """Exact transmission of V0 sech^2(r/a) (Poeschl-Teller)."""
    k = np.sqrt(np.maximum(2.0 * MU * energy, 0.0)) / HBAR
    lam = 8.0 * MU * v0 * a**2 / HBAR**2
    s = np.sinh(np.pi * k * a) ** 2
    if lam > 1.0:
        c = np.cosh(0.5 * np.pi * np.sqrt(lam - 1.0)) ** 2
    else:
        c = np.cos(0.5 * np.pi * np.sqrt(1.0 - lam)) ** 2
    return s / (s + c)


def symbols(x, y, v0, a):
    """Return (M, M_cl, M_res) at midpoint x over half separations y."""
    d = V(x + y, v0, a) - V(x - y, v0, a)
    dcl = 2.0 * y * dV(x, v0, a)
    return (1j / HBAR) * d, (1j / HBAR) * dcl, (1j / HBAR) * (d - dcl)


# --------------------------------------------------------------------- #
# Part A -- the split                                                    #
# --------------------------------------------------------------------- #
def part_a():
    banner("Part A  the compensated split for the Eckart barrier")
    v0, a = 1.0, 1.0
    print(f"  V0 = {v0}, a = {a}, mu = {MU}, beta = {beta_of(v0, a):.3f}")
    print("\n  C2: the residual is the odd part of the cubic Taylor")
    print("  remainder, bounded by (2/hbar)(y_max^3/6) sup|V'''|.\n")
    print(f"  {'y_max':>7} {'x':>5} {'max|M|':>10} {'max|M_res|':>12}"
          f" {'C2 bound':>10} {'ratio':>7}")
    for y_max in (0.1, 0.25, 0.5, 1.0, 2.0):
        y = np.linspace(-y_max, y_max, 4001)
        x = 0.5
        m, _, mr = symbols(x, y, v0, a)
        zz = np.linspace(x - y_max, x + y_max, 2001)
        bound = (2.0 / HBAR) * y_max**3 / 6.0 * np.abs(d3V(zz, v0, a)).max()
        print(f"  {y_max:7.2f} {x:5.2f} {np.abs(m).max():10.4f}"
              f" {np.abs(mr).max():12.4e} {bound:10.4f}"
              f" {np.abs(mr).max() / bound:7.3f}")

    print("\n  C3: the full kernel carries the whole classical force, the")
    print("  residual carries none.\n")
    print(f"  {'y_max':>7} {'x':>5} {'|V_prime|':>11} {'|1st mom K|':>12}"
          f" {'|1st mom Kres|':>15} {'|0th mom Kres|':>15}")
    n_s = 8192
    for y_max in (0.25, 0.5, 1.0):
        s_max = 2.0 * y_max / HBAR
        s = np.linspace(-s_max, s_max, n_s, endpoint=False)
        ds = s[1] - s[0]
        xi = np.fft.fftshift(2.0 * np.pi * np.fft.fftfreq(n_s, d=ds))
        dxi = xi[1] - xi[0]
        for x in (0.5, 1.0):
            y = HBAR * s / 2.0
            m, _, mr = symbols(x, y, v0, a)
            kk = np.fft.fftshift(
                np.fft.fft(np.fft.ifftshift(m)) * ds / (2.0 * np.pi))
            kr = np.fft.fftshift(
                np.fft.fft(np.fft.ifftshift(mr)) * ds / (2.0 * np.pi))
            m1 = abs(float(np.real(np.sum(xi * kk) * dxi)))
            m1r = abs(float(np.real(np.sum(xi * kr) * dxi)))
            m0r = abs(float(np.real(np.sum(kr) * dxi)))
            print(f"  {y_max:7.2f} {x:5.2f} {abs(dV(x, v0, a)):11.6f}"
                  f" {m1:12.6f} {m1r:15.3e} {m0r:15.3e}")


# --------------------------------------------------------------------- #
# Part B -- Theorem K1, the analyticity ceiling                          #
# --------------------------------------------------------------------- #
def taylor_coeffs(f, rho, n_pts=16384, n_max=400):
    """Taylor coefficients of f about 0 by the Cauchy integral."""
    th = 2.0 * np.pi * np.arange(n_pts) / n_pts
    c = np.fft.fft(f(rho * np.exp(1j * th))) / n_pts
    n = np.arange(n_max)
    return n, c[:n_max] / rho**n  # Taylor coefficients c_n


def radius_estimate(f, rho, n_lo=40):
    """Radius of convergence from |c_n| ~ C n^alpha R^-n.

    Three points matter.  The circle must sit close to the singularity, or
    the coefficients hit roundoff before the asymptotic regime; the
    threshold must be applied to the raw transform coefficients, not to
    |c_n| itself, whose noise floor grows with n; and the n^alpha factor
    of an algebraic singularity must be fitted, not ignored.
    """
    n, c = taylor_coeffs(f, rho)
    c_raw = np.abs(c) * rho**n
    sel = (n >= n_lo) & (c_raw > c_raw.max() * 1e-12)
    if sel.sum() < 10:
        return float("nan")
    nn = n[sel].astype(float)
    lhs = np.log(c_raw[sel]) - nn * np.log(rho)
    amat = np.column_stack([np.ones_like(nn), np.log(nn), -nn])
    sol, *_ = np.linalg.lstsq(amat, lhs, rcond=None)
    return float(np.exp(sol[2]))


def part_b():
    banner("Part B  Theorem K1: the reach ceiling is the analyticity strip")
    print("  sech^2(z/a) has double poles at z = i pi a (n + 1/2), so")
    print("  y -> V(x + y) is analytic only for |y| < R(x) =")
    print("  sqrt(x^2 + (pi a / 2)^2).  x = 0 is excluded: the residual")
    print("  vanishes there by parity.\n")
    print(f"  {'a':>5} {'x':>5} {'predicted R':>13} {'measured':>11}"
          f" {'ratio':>7}")
    for a in (1.0, 2.0):
        for x in (0.3, 0.7, 1.0):
            def f(yy, x=x, a=a):
                return (V(x + yy, 1.0, a) - V(x - yy, 1.0, a)
                        - 2.0 * yy * dV(x, 1.0, a))
            pred = np.sqrt(x**2 + (np.pi * a / 2.0) ** 2)
            est = radius_estimate(f, 0.92 * pred)
            print(f"  {a:5.1f} {x:5.2f} {pred:13.4f} {est:11.4f}"
                  f" {est / pred:7.4f}")

    print("\n  Direct demonstration: partial sums of the Moyal (Taylor in y)")
    print("  series against the exact residual, a = 1, x = 0.3,")
    print(f"  R = {np.sqrt(0.3**2 + (np.pi / 2) ** 2):.4f}\n")
    a, x = 1.0, 0.3

    def fres(yy):
        return (V(x + yy, 1.0, a) - V(x - yy, 1.0, a)
                - 2.0 * yy * dV(x, 1.0, a)) / HBAR

    _, c = taylor_coeffs(fres, 0.8, n_max=45)
    c = np.real(c)
    print(f"  {'y':>6} {'exact':>13} {'S_10':>13} {'S_20':>13}"
          f" {'S_30':>13}   status")
    for y in (0.5, 1.0, 1.4, 1.7, 2.2):
        row = f"  {y:6.2f} {fres(y):13.5e}"
        for n_trunc in (10, 20, 30):
            s = float(np.sum(c[:n_trunc + 1] * y ** np.arange(n_trunc + 1)))
            row += f" {s:13.5e}"
        status = "inside R" if y < np.sqrt(x**2 + (np.pi * a / 2) ** 2) \
            else "OUTSIDE R"
        print(row + f"   {status}")

    print("\n  Corollary K1.2: for soft-core Coulomb -1/sqrt(z^2 + eps^2)")
    print("  the softening length is the ceiling, R = sqrt(x^2 + eps^2).\n")
    print(f"  {'eps':>6} {'x':>5} {'predicted R':>13} {'measured':>11}")
    for eps in (0.3, 0.6, 1.0):
        x = 0.4

        def fc(yy, eps=eps, x=x):
            vc = -1.0 / np.sqrt((x + yy) ** 2 + eps**2)
            vc2 = -1.0 / np.sqrt((x - yy) ** 2 + eps**2)
            dvc = x / (x**2 + eps**2) ** 1.5
            return vc - vc2 - 2.0 * yy * dvc

        pred = np.sqrt(x**2 + eps**2)
        print(f"  {eps:6.2f} {x:5.2f} {pred:13.4f}"
              f" {radius_estimate(fc, 0.92 * pred):11.4f}")

    print("\n  Corollary K1.1: y_max < pi a / 2 forces dp > hbar / a, so the")
    print("  barrier momentum scale p_b = hbar beta / a spans fewer than")
    print("  beta rungs.  Packet resolution needs dp < sigma_p = hbar /")
    print("  (2 sigma_r), i.e. y_max > pi sigma_r, so sigma_r < a / 2.\n")
    print(f"  {'beta':>6} {'a':>6} {'dp_min':>9} {'p_b':>8}"
          f" {'rungs < ':>9} {'max sigma_r':>12}")
    for a in (1.0, 2.0, 4.0, 8.0):
        v0 = 1.0
        b = beta_of(v0, a)
        dp_min = HBAR / a
        p_b = np.sqrt(2.0 * MU * v0)
        print(f"  {b:6.2f} {a:6.1f} {dp_min:9.4f} {p_b:8.4f}"
              f" {p_b / dp_min:9.2f} {a / 2.0:12.2f}")


# --------------------------------------------------------------------- #
# Part C -- Theorem K2, the far-field law                                #
# --------------------------------------------------------------------- #
def part_c():
    banner("Part C  Theorem K2: exact far-field law for the residual")
    print("  For x >> a,  max_y |M_res| -> (8 V0 / hbar) e^(-2x/a)")
    print("  [ sinh(2 y_max / a) - 2 y_max / a ].  The bracket is")
    print("  sin u - u of Lemma C0 continued to u = 2 i y / a: the")
    print("  exponential tail is a mode at imaginary wavenumber.\n")
    v0, a = 1.0, 1.0
    print(f"  {'y_max':>7} {'x':>5} {'measured':>13} {'predicted':>13}"
          f" {'ratio':>9}")
    for y_max in (0.5, 1.0, 2.0):
        for x in (4.0, 6.0, 8.0):
            y = np.linspace(-y_max, y_max, 8001)
            _, _, mr = symbols(x, y, v0, a)
            meas = float(np.abs(mr).max())
            pred = (8.0 * v0 / HBAR) * np.exp(-2.0 * x / a) * (
                np.sinh(2.0 * y_max / a) - 2.0 * y_max / a)
            print(f"  {y_max:7.2f} {x:5.1f} {meas:13.5e} {pred:13.5e}"
                  f" {meas / pred:9.5f}")
    print("\n  Consequence: the residual decays at exactly V's own rate 2/a,")
    print("  so for an exponential tail the reach rescales the interaction")
    print("  profile rather than translating it.  Errata for")
    print("  compensated_liouville_splitting.md section 5.1, whose")
    print("  'translated outward by exactly the reach' is a Gaussian fact.\n")
    print(f"  {'x':>5} {'M_res/V3(x)':>13} {'M_res/V3(x-ymax)':>18}"
          f"   (y_max = 1)")
    for x in (4.0, 5.0, 6.0, 7.0):
        y = np.linspace(-1.0, 1.0, 8001)
        _, _, mr = symbols(x, y, v0, a)
        meas = float(np.abs(mr).max())
        print(f"  {x:5.1f} {meas / abs(d3V(x, v0, a)):13.5f}"
              f" {meas / abs(d3V(x - 1.0, v0, a)):18.5f}")


# --------------------------------------------------------------------- #
# Part D -- Theorem K3, the budget ratio                                 #
# --------------------------------------------------------------------- #
def budget_ratio(y_max, a):
    """Spectrum-weighted |M_res| / |M_cl| for V0 sech^2(x/a)."""
    k = np.linspace(1e-7, 400.0 / a, 400001)
    vt = np.abs(np.pi * a * (k * a) / np.sinh(np.pi * k * a / 2.0))
    u = k * y_max
    num = np.trapezoid(vt * np.abs(np.sin(u) - u), k)
    den = np.trapezoid(vt * np.abs(u), k)
    return float(num / den)


def part_d():
    banner("Part D  Theorem K3: the compensated channel never loses")
    print("  Vtilde(k) = V0 pi a^2 k / sinh(pi k a / 2).  Per mode the")
    print("  ratio |sin u - u| / |u| peaks at 1.217 near u = 4.493, but")
    print("  the sech^2 spectrum suppresses those modes: weighted, the")
    print("  ratio saturates at 1 from below for every reach.  Contrast")
    print("  compensated_liouville_splitting.md section 6.3, where a ring")
    print("  pins the ratio at exactly 1 for every mode.\n")
    print(f"  {'y_max/a':>9} {'a = 1':>10} {'a = 2':>10} {'a = 4':>10}")
    ratios = []
    grid = (0.1, 0.25, 0.5, 1.0, 1.5708, 2.0, 3.0, 4.0, 6.0)
    for r in grid:
        row = [budget_ratio(r * a, a) for a in (1.0, 2.0, 4.0)]
        ratios.append(row[0])
        print(f"  {r:9.4f} {row[0]:10.4f} {row[1]:10.4f} {row[2]:10.4f}")
    u = np.linspace(0.01, 20.0, 200001)
    per_mode = np.abs(np.sin(u) - u) / u
    print(f"\n  per-mode maximum {per_mode.max():.4f} at u ="
          f" {u[np.argmax(per_mode)]:.4f}")
    print(f"  weighted maximum over the table {max(ratios):.4f}")
    return grid, ratios


# --------------------------------------------------------------------- #
# Exact spectral evolution                                               #
# --------------------------------------------------------------------- #
class Run:
    """Wigner evolution of the relative coordinate on an (r, p) grid.

    The momentum axis is held in FFT-natural order so that the potential
    substep, the residual convolution and the separatrix correlations are
    all plain transforms along axis 1.
    """

    def __init__(self, v0, a, sigma_r, r_c, p_c, n_r=1024, n_p=256,
                 r_half=72.0, dp=0.05):
        self.v0, self.a = v0, a
        self.r = -r_half + 2.0 * r_half * np.arange(n_r) / n_r
        self.dr = self.r[1] - self.r[0]
        self.kr = 2.0 * np.pi * np.fft.fftfreq(n_r, d=self.dr)
        self.n_r, self.n_p, self.dp = n_r, n_p, dp
        self.jj = np.fft.fftfreq(n_p, d=1.0 / n_p).astype(int)
        # half-cell offset keeps the classical threshold off a grid line
        self.p = dp * (self.jj + 0.5)
        self.s = 2.0 * np.pi * np.fft.fftfreq(n_p, d=dp)
        self.y = HBAR * self.s / 2.0
        self.y_max = float(np.abs(self.y).max())
        rr = self.r[:, None]
        yy = self.y[None, :]
        d = V(rr + yy, v0, a) - V(rr - yy, v0, a)
        dcl = 2.0 * yy * dV(rr, v0, a)
        self.m_full = (1j / HBAR) * d
        self.m_cl = (1j / HBAR) * dcl
        self.m_res = (1j / HBAR) * (d - dcl)
        self.k_res = np.fft.ifft(self.m_res, axis=1)
        sig_p = HBAR / (2.0 * sigma_r)
        w = np.exp(-((self.r[:, None] - r_c) ** 2) / (2.0 * sigma_r**2)
                   - ((self.p[None, :] - p_c) ** 2) / (2.0 * sig_p**2))
        self.w0 = w / (np.sum(w) * self.dr * self.dp)
        self.sigma_p = sig_p
        # Classical outcome set: a world ends on the far side iff it clears
        # the barrier moving right, or is already past it with too little
        # energy to come back.  Exactly invariant under the classical flow,
        # and equal to the transmission once the packet is asymptotic --
        # the same device that makes the inverted pair barrier readable at
        # t = 0.  The energy test is weighted at sub-cell resolution;
        # a hard test jitters by O(dp) as the packet crosses it.
        p_star = np.sqrt(np.maximum(
            2.0 * MU * (v0 - V(self.r[:, None], v0, a)), 0.0))
        over = np.clip((np.abs(self.p[None, :]) - p_star) / dp + 0.5, 0.0, 1.0)
        right = (self.p[None, :] > 0).astype(float)
        past = (self.r[:, None] > 0).astype(float)
        self.inside = over * right + (1.0 - over) * past
        self.outside = 1.0 - self.inside
        self.a_hat = np.fft.fft(self.inside, axis=1)
        self.b_hat = np.fft.fft(self.outside, axis=1)

    def _p_mult(self, w, sym, dt):
        return np.real(np.fft.ifft(np.fft.fft(w, axis=1)
                                   * np.exp(dt * sym), axis=1))

    def _stream(self, w, dt):
        wh = np.fft.fft(w, axis=0)
        wh = wh * np.exp(-1j * self.kr[:, None] * self.p[None, :] * dt / MU)
        return np.real(np.fft.ifft(wh, axis=0))

    def evolve(self, mode, t_max, dt, flux_every=0):
        """mode in {'full', 'classical', 'compensated'}."""
        w = self.w0.copy()
        n_steps = int(round(t_max / dt))
        times, sigma_int, phi_in, phi_out = [], [], [], []
        for step in range(n_steps):
            w = self._stream(w, 0.5 * dt)
            if mode == "full":
                w = self._p_mult(w, self.m_full, dt)
            elif mode == "classical":
                w = self._p_mult(w, self.m_cl, dt)
            else:
                w = self._p_mult(w, self.m_cl, dt)
                if flux_every and step % flux_every == 0:
                    fi, fo = self.separatrix_flux(w)
                    phi_in.append(fi)
                    phi_out.append(fo)
                    times.append(step * dt)
                    sigma_int.append(self.sigma_integral(w))
                w = self._p_mult(w, self.m_res, dt)
            w = self._stream(w, 0.5 * dt)
        if flux_every:
            return w, (np.array(times), np.array(sigma_int),
                       np.array(phi_in), np.array(phi_out))
        return w

    def sigma_integral(self, w):
        return float(np.sum(w * self.inside) * self.dr * self.dp)

    def transmitted(self, w, cut):
        return float(np.sum(w[self.r > cut, :]) * self.dr * self.dp)

    def _corr(self, ind_hat, f):
        """sum_j ind[j + m] f[j], for every offset m, row by row."""
        return np.real(np.fft.ifft(ind_hat * np.conj(np.fft.fft(f, axis=1)),
                                   axis=1))

    def separatrix_flux(self, w):
        """Signed rates of weight entering and leaving the classical set.

        Uses sum_q K_q = 0 to replace the raw gain by the crossing form
        sum_q K_q [1_S(p + xi_q) - 1_S(p)] W(p), then splits the bracket
        into its +1 (in) and -1 (out) parts.
        """
        kr = np.real(self.k_res)
        g_in = self._corr(self.a_hat, self.outside * w)
        g_out = self._corr(self.b_hat, self.inside * w)
        phi_in = float(np.sum(kr * g_in)) * self.dr * self.dp
        phi_out = float(np.sum(kr * g_out)) * self.dr * self.dp
        return phi_in, phi_out


# --------------------------------------------------------------------- #
# Part E -- Theorem K4                                                   #
# --------------------------------------------------------------------- #
def part_e():
    banner("Part E  Theorem K4: the whole gap is the residual channel")
    v0, a = 1.0, 2.0
    sigma_r = a
    p_b = np.sqrt(2.0 * MU * v0)
    run = Run(v0, a, sigma_r, r_c=-16.0, p_c=p_b)
    print(f"  V0 = {v0}, a = {a}, beta = {beta_of(v0, a):.3f}")
    print(f"  packet sigma_r = {sigma_r}, sigma_p = {run.sigma_p:.4f},"
          f" p_c = p_b = {p_b:.4f}")
    print(f"  grid {run.n_r} x {run.n_p}, dp = {run.dp}, reach"
          f" y_max = {run.y_max:.2f}")
    print("  (the dynamics run is at effectively unbounded reach; by")
    print("   Corollary K1.1 a reach-limited lattice cannot resolve this")
    print("   packet, which is itself one of the note's results)\n")

    dt, t_max = 0.02, 22.0

    # Theorem K4.  The classical outcome functional is exactly invariant
    # under streaming plus deterministic acceleration, so the classical
    # transmission needs no evolution at all: it is the functional's
    # initial value.  This is the same device that makes the inverted pair
    # barrier readable at t = 0, transplanted to a potential that does have
    # a jump channel.
    t_cl = run.sigma_integral(run.w0)

    w_full = run.evolve("full", t_max, dt)
    # The flux ledger reads the generator, not the exponential, so its
    # per-step error is O(dt) and accumulates coherently: 1.4 per cent per
    # step at dt = 0.02, 0.34 per cent at dt = 0.005.  The flux leg is run
    # at the finer step for that reason alone.
    dt_flux = 0.005
    w_comp, series = run.evolve("compensated", t_max, dt_flux, flux_every=1)
    t_full = run.sigma_integral(w_full)
    t_comp = run.sigma_integral(w_comp)

    pp = np.linspace(1e-6, p_b + 12 * run.sigma_p, 400001)
    wt = np.exp(-((pp - p_b) ** 2) / (2.0 * run.sigma_p**2))
    wt /= np.trapezoid(wt, pp)
    e = pp**2 / (2.0 * MU)
    t_closed = float(np.trapezoid(wt * T_eckart(e, v0, a), pp))
    t_cl_closed = float(np.trapezoid(wt * (e > v0), pp))

    print(f"  {'quantity':<38} {'value':>10}")
    print(f"  {'classical, = functional at t = 0':<38} {t_cl:10.6f}")
    print(f"  {'closed form, classical':<38} {t_cl_closed:10.6f}")
    print(f"  {'full symbol run, final':<38} {t_full:10.6f}")
    print(f"  {'compensated product run, final':<38} {t_comp:10.6f}")
    print(f"  {'closed form, exact T(E)':<38} {t_closed:10.6f}")
    print(f"  {'C1 check |full - compensated|':<38}"
          f" {abs(t_full - t_comp):10.3e}")
    print(f"  {'closed form check |run - exact|':<38}"
          f" {abs(t_full - t_closed):10.3e}")
    print(f"  {'gap delivered by the residual':<38}"
          f" {t_full - t_cl:10.6f}")
    print(f"  {'closed-form gap':<38}"
          f" {t_closed - t_cl_closed:10.6f}")

    print("\n  Conservation of the functional under streaming plus")
    print("  deterministic acceleration alone:\n")
    print(f"    {'t':>6} {'functional':>12} {'min W':>11}")
    w_probe = run.w0.copy()
    for i in range(int(8.0 / dt) + 1):
        if i % int(2.0 / dt) == 0:
            print(f"    {i * dt:6.2f} {run.sigma_integral(w_probe):12.6f}"
                  f" {w_probe.min():11.3e}")
        w_probe = run._stream(w_probe, 0.5 * dt)
        w_probe = run._p_mult(w_probe, run.m_cl, dt)
        w_probe = run._stream(w_probe, 0.5 * dt)
    print("\n  Beyond t ~ 10 the classical flow filaments below the grid,")
    print("  min W reaching -1.6e-1 and the check losing meaning.  That is")
    print("  a limit of the classical *reference*: the quantum run does not")
    print("  filament, because the residual channel smooths in p.")
    return run, series, (t_cl, t_full, t_closed, t_cl_closed), w_full


# --------------------------------------------------------------------- #
# Part F -- Theorem K5, the traffic                                      #
# --------------------------------------------------------------------- #
def part_f(run, series, transmissions):
    banner("Part F  Theorem K5: a small imbalance in two large flows")
    times, sig, phi_in, phi_out = series
    t_cl, t_full, t_closed, t_cl_closed = transmissions
    dt = times[1] - times[0]
    net_rate = phi_in - phi_out
    gross_rate = np.abs(phi_in) + np.abs(phi_out)
    cum_net = np.concatenate([[0.0], np.cumsum(
        0.5 * (net_rate[1:] + net_rate[:-1])) * dt])[:len(times)]
    cum_gross = np.concatenate([[0.0], np.cumsum(
        0.5 * (gross_rate[1:] + gross_rate[:-1])) * dt])[:len(times)]
    print("  The residual channel moves weight both ways across the")
    print("  classical separatrix E = V0.  Time-integrated:\n")
    print(f"  {'gross traffic  int(|Phi_in| + |Phi_out|) dt':<46}"
          f" {cum_gross[-1]:10.6f}")
    print(f"  {'net transfer   int(Phi_in - Phi_out) dt':<46}"
          f" {cum_net[-1]:10.6f}")
    print(f"  {'measured gap   T_full - T_classical':<46}"
          f" {t_full - t_cl:10.6f}")
    print(f"  {'net / gross':<46}"
          f" {cum_net[-1] / cum_gross[-1]:10.6f}")
    print("\n  Peak instantaneous flows:")
    print(f"    max |Phi_in|  = {np.abs(phi_in).max():.6f}")
    print(f"    max |Phi_out| = {np.abs(phi_out).max():.6f}")
    print(f"    max |Phi_in - Phi_out| = "
          f"{np.abs(phi_in - phi_out).max():.6f}")
    return times, phi_in, phi_out, cum_net


# --------------------------------------------------------------------- #
# Part G -- Theorem K6, the 1 / beta law                                 #
# --------------------------------------------------------------------- #
def part_g():
    banner("Part G  Theorem K6: gross saturates, net falls as 1 / beta")
    print("  Packet centred on the barrier top, sigma_r = a, so the")
    print("  classical transmission is exactly 1/2 and the two flows are")
    print("  equal and opposite in energy.  The residue is the Jacobian")
    print("  dp/dE = mu/p across the tunnelling window, of relative size")
    print("  1 / beta.\n")
    print(f"  {'beta':>7} {'a':>7} {'T_class':>9} {'T_quant':>9}"
          f" {'net':>10} {'tunnel-in':>10} {'over-refl':>10}"
          f" {'net/gross':>10} {'x beta':>8}")
    rows = []
    for a in (0.5, 1.0, 2.0, 4.0, 8.0, 16.0):
        v0 = 1.0
        b = beta_of(v0, a)
        p_b = np.sqrt(2.0 * MU * v0)
        sig_p = HBAR / (2.0 * a)
        pp = np.linspace(1e-9, p_b + 14 * sig_p, 600001)
        wt = np.exp(-((pp - p_b) ** 2) / (2.0 * sig_p**2))
        wt /= np.trapezoid(wt, pp)
        e = pp**2 / (2.0 * MU)
        tt = T_eckart(e, v0, a)
        t_q = float(np.trapezoid(wt * tt, pp))
        t_c = float(np.trapezoid(wt * (e > v0), pp))
        lo = e < v0
        tun = float(np.trapezoid(wt[lo] * tt[lo], pp[lo]))
        ovr = float(np.trapezoid(wt[~lo] * (1.0 - tt[~lo]), pp[~lo]))
        gross = tun + ovr
        rows.append((b, t_c, t_q, t_q - t_c, tun, ovr, gross))
        print(f"  {b:7.2f} {a:7.1f} {t_c:9.5f} {t_q:9.5f}"
              f" {t_q - t_c:+10.5f} {tun:10.5f} {ovr:10.5f}"
              f" {(t_q - t_c) / gross:10.5f}"
              f" {(t_q - t_c) / gross * b:8.4f}")
    print("\n  The last column is (net/gross) x beta.  Constant means")
    print("  net/gross ~ c / beta.")
    return rows


# --------------------------------------------------------------------- #
# Part H -- world-particle trajectories                                  #
# --------------------------------------------------------------------- #
def residual_lattice(v0, a, y_max, r_grid, n_rungs=48, n_y=1024):
    """K_q(r) on the momentum lattice of reach y_max, hard window.

    D_res is real and odd in y with period 2 y_max, so its Fourier
    coefficients are imaginary and odd in q and K_q = -(i/hbar) c_q is real
    and odd -- the signed kernel of section 2.2.
    """
    y = -y_max + 2.0 * y_max * np.arange(n_y) / n_y
    rr = r_grid[:, None]
    d = (V(rr + y[None, :], v0, a) - V(rr - y[None, :], v0, a)
         - 2.0 * y[None, :] * dV(rr, v0, a))
    q = np.fft.fftfreq(n_y, d=1.0 / n_y).astype(int)
    c = np.fft.fft(d, axis=1) / n_y * ((-1.0) ** q)[None, :]
    k = np.real(-1j * c / HBAR)
    keep = np.abs(q) <= n_rungs
    order = np.argsort(q[keep])
    return q[keep][order], np.pi * HBAR / (2.0 * y_max), k[:, keep][:, order]


def sample_worlds(v0, a, y_max, n_parents=4, t_max=9.0, dt=0.004,
                  seed=7, max_drawn=4):
    """Newtonian arcs punctuated by positon-negaton emissions.

    The residual kernel is signed, so no one-body jump process exists
    (``compensated_liouville_splitting.md`` section 2.2, consequence 3).
    What is drawn is the pair-generation form: the parent streams
    deterministically under the full classical force and at rate
    Gamma(r) = sum_q |K_q(r)| emits a positon-negaton pair at p +/- xi_q,
    of zero expected weight.  A hard window is used, so the picture is
    illustrative rather than exact.
    """
    rng = np.random.default_rng(seed)
    r_grid = np.linspace(-8.0 * a, 8.0 * a, 801)
    q, dp, kq = residual_lattice(v0, a, y_max, r_grid)
    nz = q != 0
    rate_grid = np.sum(np.abs(kq[:, nz]), axis=1)
    xi = q[nz] * dp

    p_b = np.sqrt(2.0 * MU * v0)
    sig_r, sig_p = a, HBAR / (2.0 * a)
    tracks = []
    for _ in range(n_parents):
        r = -6.0 * a + rng.normal(0.0, sig_r)
        p = p_b + rng.normal(0.0, sig_p)
        path_r, path_p, events = [r], [p], []
        t = 0.0
        while t < t_max:
            p = p - dV(r, v0, a) * 0.5 * dt
            r = r + p / MU * dt
            p = p - dV(r, v0, a) * 0.5 * dt
            t += dt
            rate = float(np.interp(r, r_grid, rate_grid))
            if rate > 0.0 and rng.random() < rate * dt:
                i0 = int(np.clip(np.searchsorted(r_grid, r), 1,
                                 len(r_grid) - 1))
                w = np.abs(kq[i0, nz])
                if w.sum() > 0:
                    idx = rng.choice(len(xi), p=w / w.sum())
                    events.append((r, p, float(abs(xi[idx]))))
            path_r.append(r)
            path_p.append(p)
        if len(events) > max_drawn:
            pick = rng.choice(len(events), max_drawn, replace=False)
            events = [events[j] for j in sorted(pick)]
        tracks.append((np.array(path_r), np.array(path_p), events))
    return tracks, dp, rate_grid.max()


def part_h():
    banner("Part H  sample world-particle trajectories")
    v0, a = 1.0, 2.0
    y_max = np.pi * a / 2.0
    print(f"  V0 = {v0}, a = {a}, reach y_max = pi a / 2 = {y_max:.4f}")
    print(f"  dp = pi hbar / (2 y_max) = {np.pi * HBAR / (2 * y_max):.4f}")
    print("  the reach is set at the ceiling of Theorem K1; a hard window")
    print("  is used, so this is illustrative rather than exact\n")
    tracks, dp, rate_max = sample_worlds(v0, a, y_max)
    n_ev = sum(len(t[2]) for t in tracks)
    print(f"  peak emission rate Gamma = {rate_max:.4f}")
    print(f"  {len(tracks)} parents, {n_ev} emission events")
    for i, (rr, pp, ev) in enumerate(tracks):
        print(f"    parent {i}: r {rr[0]:+7.2f} -> {rr[-1]:+7.2f},"
              f"  p {pp[0]:+6.3f} -> {pp[-1]:+6.3f},"
              f"  {len(ev)} events")
    return tracks


# --------------------------------------------------------------------- #
# Figures                                                                #
# --------------------------------------------------------------------- #
def fig_split(grid, ratios):
    v0, a = 1.0, 1.0
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.0))
    z = np.linspace(-6, 6, 801)
    ax = axes[0]
    ax.plot(z, V(z, v0, a), label=r"$V$")
    ax.plot(z, d3V(z, v0, a) / 4.0, label=r"$V'''/4$")
    ax.axhline(0.0, lw=0.6, color="0.6")
    ax.set_xlabel("$x/a$")
    ax.set_title("Eckart barrier and its third derivative")
    ax.legend(fontsize=9)

    ax = axes[1]
    y = np.linspace(-2.5, 2.5, 1201)
    m, mcl, mr = symbols(0.5, y, v0, a)
    ax.plot(y, np.imag(m), label=r"$M$")
    ax.plot(y, np.imag(mcl), "--", label=r"$M_{\rm cl}$")
    ax.plot(y, np.imag(mr), label=r"$M_{\rm res}$")
    ax.axvspan(-np.pi / 2, np.pi / 2, color="0.85", zorder=0)
    ax.text(0.0, -3.4, "analyticity strip", ha="center", fontsize=8)
    ax.set_xlabel("half separation $y$")
    ax.set_title(r"the split at $x=0.5$")
    ax.legend(fontsize=9)

    ax = axes[2]
    ax.plot(grid, ratios, "o-")
    ax.axhline(1.0, color="0.4", lw=0.8, ls=":")
    ax.set_xlabel(r"$y_{\max}/a$")
    ax.set_ylabel(r"$|M_{\rm res}|/|M_{\rm cl}|$")
    ax.set_ylim(0, 1.15)
    ax.set_title("budget ratio saturates at 1 from below")
    fig.tight_layout()
    save_fig(fig, "eckart_compensated_split.png")


def fig_ceiling():
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.0))
    a, x = 1.0, 0.3
    r_ceiling = np.sqrt(x**2 + (np.pi * a / 2) ** 2)

    ax = axes[0]
    th = np.linspace(0, 2 * np.pi, 400)
    ax.plot(r_ceiling * np.cos(th), r_ceiling * np.sin(th), "k--", lw=0.9,
            label=r"$|y| = R(x)$")
    ax.fill(0.6 * r_ceiling * np.cos(th), 0.6 * r_ceiling * np.sin(th),
            color="0.85", label=r"reach $y_{\max}$")
    for sx in (-1, 1):
        for sy in (-1, 1):
            ax.plot(sx * x, sy * np.pi * a / 2, "o", color="C3", ms=7)
    ax.axhline(0, lw=0.6, color="0.6")
    ax.axvline(0, lw=0.6, color="0.6")
    ax.set_aspect("equal")
    ax.set_xlabel(r"$\mathrm{Re}\,y$")
    ax.set_ylabel(r"$\mathrm{Im}\,y$")
    ax.set_title("poles of $V$ bound the reach")
    ax.legend(fontsize=8, loc="upper right")

    ax = axes[1]

    def fres(yy):
        return (V(x + yy, 1.0, a) - V(x - yy, 1.0, a)
                - 2.0 * yy * dV(x, 1.0, a)) / HBAR

    _, c = taylor_coeffs(fres, 0.8, n_max=45)
    c = np.real(c)
    yy = np.linspace(0.05, 2.4, 300)
    ax.semilogy(yy, np.abs(fres(yy)), "k", label="exact")
    for n_trunc, style in ((10, "-"), (20, "--"), (30, ":")):
        s = np.array([np.sum(c[:n_trunc + 1] * v ** np.arange(n_trunc + 1))
                      for v in yy])
        ax.semilogy(yy, np.abs(s), style, label=f"$S_{{{n_trunc}}}$")
    ax.axvline(r_ceiling, color="C3", lw=1.0)
    ax.text(r_ceiling * 1.02, 1e3, "$R(x)$", color="C3", fontsize=9)
    ax.set_ylim(1e-3, 1e6)
    ax.set_xlabel("$y$")
    ax.set_title("Moyal partial sums diverge beyond $R$")
    ax.legend(fontsize=8)

    ax = axes[2]
    bb = np.linspace(0.4, 16, 200)
    ax.plot(bb, bb, label=r"$p_b/\Delta p_{\min} = \beta$")
    ax.axhline(1.0, color="0.4", ls=":", lw=0.8)
    ax.fill_between(bb, 0, np.minimum(bb, 1.0), color="0.85")
    ax.text(6.0, 0.4, "lattice cannot resolve\nthe barrier scale",
            fontsize=8, ha="center")
    ax.set_xlabel(r"$\beta$")
    ax.set_ylabel("rungs across $p_b$")
    ax.set_title(r"$\Delta p > \hbar/a$ caps the momentum lattice")
    ax.legend(fontsize=8)
    fig.tight_layout()
    save_fig(fig, "eckart_compensated_ceiling.png")


def fig_separatrix(run, tracks, times, phi_in, phi_out, cum_net,
                   transmissions, w_full):
    t_cl, t_full, t_closed, t_cl_closed = transmissions
    v0, a = run.v0, run.a
    fig = plt.figure(figsize=(13.5, 8.4))

    # (a) the real classical outcome set, straight off the run
    ax = fig.add_subplot(2, 2, 1)
    rr = np.linspace(-3.5 * a, 3.5 * a, 400)
    pb_curve = np.sqrt(np.maximum(2.0 * MU * (v0 - V(rr, v0, a)), 0.0))
    pg = np.fft.fftshift(run.p)
    ins = np.fft.fftshift(run.inside, axes=1)
    rg, pgm = np.meshgrid(run.r, pg, indexing="ij")
    ax.contourf(rg, pgm, ins, levels=[0.5, 1.5], colors=["#dbeafe"])
    ax.plot(rr, pb_curve, "k", lw=1.4)
    ax.plot(rr[rr < 0], -pb_curve[rr < 0], "k", lw=1.4)
    ax.plot([0, 0], [-1.0, 0.0], "k", lw=1.4)
    ax.text(-3.3 * a, 1.75, "ends on the far side", fontsize=9)
    ax.text(-3.3 * a, -0.75, "ends on the near side", fontsize=9)
    for r0 in (-1.3 * a, -0.15 * a, 1.0 * a):
        p0 = np.interp(r0, rr, pb_curve)
        ax.annotate("", xy=(r0, p0 + 0.42), xytext=(r0, p0 - 0.42),
                    arrowprops=dict(arrowstyle="<->", color="C3", lw=2.0))
    ax.set_xlim(-3.5 * a, 3.5 * a)
    ax.set_ylim(-1.1, 2.2)
    ax.set_xlabel("$r$")
    ax.set_ylabel("$p_r$")
    ax.set_title("(a) hops cross the classical outcome boundary")

    # (b) sample world trajectories
    ax = fig.add_subplot(2, 2, 2)
    ax.fill_between(rr, pb_curve, 2.6, color="#dbeafe")
    ax.plot(rr, pb_curve, "k", lw=1.4)
    for path_r, path_p, events in tracks:
        ax.plot(path_r, path_p, color="0.25", lw=1.1)
        for (er, ep, exi) in events:
            ax.plot([er, er], [ep, ep + exi], color="C3", lw=1.2, alpha=0.9)
            ax.plot([er, er], [ep, ep - exi], color="C0", lw=1.2, alpha=0.9)
            ax.plot(er, ep + exi, "o", color="C3", ms=4)
            ax.plot(er, ep - exi, "o", color="C0", ms=4)
    ax.plot([], [], color="0.25", lw=1.1, label="Newtonian arc")
    ax.plot([], [], color="C3", lw=1.2, label="positon child")
    ax.plot([], [], color="C0", lw=1.2, label="negaton child")
    ax.set_xlim(-3.5 * a, 3.5 * a)
    ax.set_ylim(-0.9, 2.6)
    ax.set_xlabel("$r$")
    ax.set_ylabel("$p_r$")
    ax.set_title("(b) world-particle paths and emitted pairs")
    ax.legend(fontsize=8, loc="lower right")

    # (c) the two flows
    ax = fig.add_subplot(2, 2, 3)
    ax.plot(times, phi_in, label=r"$\Phi_{\rm in}$")
    ax.plot(times, phi_out, label=r"$\Phi_{\rm out}$")
    ax.plot(times, phi_in - phi_out, "k", lw=1.4,
            label=r"$\Phi_{\rm in}-\Phi_{\rm out}$")
    ax.axhline(0.0, lw=0.6, color="0.6")
    ax.set_xlabel("$t$")
    ax.set_ylabel("rate")
    ax.set_title("(c) two large flows, small difference")
    ax.legend(fontsize=8)

    # (d) the ledger
    ax = fig.add_subplot(2, 2, 4)
    ax.plot(times, t_cl_closed + cum_net, "k", lw=1.4,
            label=r"$T_{\rm class} + \int(\Phi_{\rm in}-\Phi_{\rm out})$")
    ax.axhline(t_cl, color="C0", ls="--",
               label=f"classical run  {t_cl:.4f}")
    ax.axhline(t_full, color="C3", ls="--",
               label=f"full run  {t_full:.4f}")
    ax.axhline(t_closed, color="C2", ls=":",
               label=f"closed form  {t_closed:.4f}")
    ax.set_xlabel("$t$")
    ax.set_ylabel("transmission")
    ax.set_title("(d) the imbalance integrates to the quantum answer")
    ax.legend(fontsize=8, loc="lower right")
    fig.tight_layout()
    save_fig(fig, "eckart_compensated_separatrix.png")


def fig_scaling(rows):
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    b = np.array([r[0] for r in rows])
    net = np.array([r[3] for r in rows])
    tun = np.array([r[4] for r in rows])
    ovr = np.array([r[5] for r in rows])
    gross = np.array([r[6] for r in rows])

    ax = axes[0]
    ax.plot(b, tun, "o-", label="tunnel in")
    ax.plot(b, ovr, "s-", label="over-barrier reflection")
    ax.plot(b, gross, "^-", label="gross traffic")
    ax.plot(b, net, "k.-", label="net")
    ax.set_xscale("log")
    ax.set_xlabel(r"$\beta$")
    ax.set_ylabel("phase-space weight")
    ax.set_title("gross saturates, net does not")
    ax.legend(fontsize=8)

    ax = axes[1]
    ax.loglog(b, net / gross, "o-", label="net / gross")
    ax.loglog(b, 0.55 / b, "k:", label=r"$0.55/\beta$")
    ax.set_xlabel(r"$\beta$")
    ax.set_title("the cancellation tightens as $1/\\beta$")
    ax.legend(fontsize=8)
    fig.tight_layout()
    save_fig(fig, "eckart_compensated_scaling.png")


# --------------------------------------------------------------------- #
def main():
    print("Verification for docs/analysis/eckart_barrier_compensated.md")
    part_a()
    part_b()
    part_c()
    grid, ratios = part_d()
    run, series, transmissions, w_full = part_e()
    times, phi_in, phi_out, cum_net = part_f(run, series, transmissions)
    rows = part_g()
    tracks = part_h()

    banner("Figures")
    fig_split(grid, ratios)
    fig_ceiling()
    fig_separatrix(run, tracks, times, phi_in, phi_out, cum_net,
                   transmissions, w_full)
    fig_scaling(rows)

    banner("Summary")
    print("  K1  the reach ceiling is the analyticity strip of V; for")
    print("      sech^2 it is uniform, y_max < pi a / 2, hence dp > hbar/a")
    print("  K2  the far field obeys the hyperbolic continuation of Lemma")
    print("      C0; for an exponential tail the reach rescales rather")
    print("      than translates the interaction profile")
    print("  K3  the spectrum-weighted budget ratio saturates at 1 from")
    print("      below: compensation never loses on the open line")
    print("  K4  the deterministic step conserves the classical")
    print("      transmission functional, so the residual channel carries")
    print("      the entire quantum correction")
    print("  K5  that correction is the small imbalance of two large")
    print("      opposed flows across the classical separatrix")
    print("  K6  gross traffic saturates in beta while the net falls as")
    print("      1 / beta: the sampling cost of tunnelling grows linearly")
    print("\ndone.")


if __name__ == "__main__":
    main()
