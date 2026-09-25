"""The sea as a phase reference -- companion demo for ladder step 22.

Verifies the note ``docs/analysis/sea_phase_reference.md`` in three parts
and draws its figures.

A  Theorem L1 -- the residual-kernel weights as a sum over contacts:
   K_q(x) = -(B dp/hbar) int dy U_res(x,y) w(y) sin(2 xi_q y/hbar), the strip
   identity B dp 2 y_max = 1, and the single-contact estimator statistics.
B  Theorem L7 -- the dephasing of a freshly locked sea.  The midpoint
   misalignment of same-row partners winds at the TRAPEZOID residual
   U - y (V'(x_i) + V'(x_j)); the pair sum of that rate is force-free but
   is not the compensated kernel: its third moment is about -2 times the
   true one.
C  Figures: the locked phase field over phase space for the Eckart barrier
   and the soft-core Coulomb well at three times, and locked against gas.
D  Theorem L8 -- the reader's frame.  A symbolic check of the identity
   hbar dmu/dt = V(x_j) - V(x_i) + 2y dp_ref/dt + (p_j-p_i)(2p_ref-p_i-p_j)/2m,
   the two quadrature residuals to fifth order, the three readers at
   re-lock (inertial -> U, parent -> U_res, own -> trapezoid), the drift of
   the parent-frame reading after re-lock, and the figure
   sea_phase_reader_frames.png.

The expensive measurements of the note live in two further scripts:
``src/demo_contact_kernel.py`` (section 2) and
``src/demo_sea_lock_particles.py`` (sections 5 and 6).

Usage::

    PYTHONPATH=src python3 -u src/demo_sea_lock.py
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import hsv_to_rgb

from demo_emission_and_absorption import Ledger, V as V_eckart, B, MU, HBAR
from wpmwlib.wpmw_utils import output_path, docs_path


def banner(text):
    print("\n" + "=" * 72 + "\n" + text + "\n" + "=" * 72, flush=True)


def eckart(x):
    return 1.0 / np.cosh(x) ** 2


def eckart_force(x):                          # -V'(x)
    return 2.0 * np.tanh(x) / np.cosh(x) ** 2


def coulomb(x):
    return -1.0 / np.sqrt(x * x + 1.0)


def coulomb_force(x):
    return -x / (x * x + 1.0) ** 1.5


# ----------------------------------------------------------------------
# A.  Contact decomposition
# ----------------------------------------------------------------------
def part_a():
    banner("A. Theorem L1: the kernel as a sum over contacts")
    run = Ledger(v0=1.0, a=1.0, n_r=128, r_half=20.0, n_p=64, dp=0.25)
    dp, ymax = run.dp, run.y_max
    print(f"   B = 1/(pi hbar) = {B:.6f}; dp = {dp}; y_max = {ymax:.4f}")
    print(f"   strip identity: B * dp * 2 y_max = {B * dp * 2 * ymax:.6f}")
    y = np.linspace(-ymax, ymax, 200001)
    w = np.cos(np.pi * y / (2 * ymax)) ** 2
    i = int(np.argmin(np.abs(run.r - (-0.625))))
    x = run.r[i]
    u_res = V_eckart(x + y, 1.0, 1.0) - V_eckart(x - y, 1.0, 1.0) \
        - 2 * y * run.dv_eff[i]
    print(f"   x = {x:+.4f}")
    print(f"   {'q':>3} {'K_q (code)':>12} {'K_q (integral)':>15} {'rel diff':>9}")
    for q in (1, 2, 3, 4, 6, 8):
        xi = q * dp
        integ = -(B * dp / HBAR) * np.trapezoid(
            u_res * w * np.sin(2 * xi * y / HBAR), y)
        kc = run.k[i, q]
        print(f"   {q:3d} {kc:+12.6f} {integ:+15.6f}"
              f" {abs(integ - kc) / abs(kc):9.1e}")
    rng = np.random.default_rng(0)
    xi = dp

    def one(yy):
        return -(B * dp / HBAR) * 2 * ymax * (
            (V_eckart(x + yy, 1, 1) - V_eckart(x - yy, 1, 1)
             - 2 * yy * run.dv_eff[i])
            * np.cos(np.pi * yy / (2 * ymax)) ** 2 * np.sin(2 * xi * yy / HBAR))

    print(f"   single-contact estimator of K_1 = {run.k[i, 1]:+.5f}:")
    for n in (1, 100, 10000):
        est = np.array([one(rng.uniform(-ymax, ymax, n)).mean()
                        for _ in range(2000)])
        print(f"   N = {n:5d}: mean {est.mean():+.5f}, std {est.std():.5f},"
              f" std*sqrt(N) = {est.std() * np.sqrt(n):.4f}")


# ----------------------------------------------------------------------
# B.  The dephasing of a locked sea
# ----------------------------------------------------------------------
def same_row_pairs(x, p, dp, ymax, xw):
    """All same-row pairs (i, j), x_i < x_j, within the reach 2 y_max,
    whose midpoint lies in the window |x| < xw."""
    row = np.round(p / dp).astype(int)
    order = np.lexsort((x, row))
    xr, rr = x[order], row[order]
    ii, jj = [], []
    for r in np.unique(rr):
        idx = np.flatnonzero(rr == r)
        xx = xr[idx]
        hi = np.searchsorted(xx, xx + 2 * ymax)
        cnt = hi - np.arange(len(xx)) - 1
        a = np.repeat(np.arange(len(xx)), cnt)
        b = a + 1 + (np.arange(cnt.sum())
                     - np.repeat(np.cumsum(cnt) - cnt, cnt))
        ii.append(order[idx[a]])
        jj.append(order[idx[b]])
    i, j = np.concatenate(ii), np.concatenate(jj)
    keep = np.abs(0.5 * (x[i] + x[j])) < xw
    return i[keep], j[keep]


def part_b():
    banner("B. Theorem L7: the dephasing of a freshly locked sea")
    run = Ledger(v0=1.0, a=1.0, n_r=192, r_half=24.0, n_p=64, dp=0.25)
    dp, npp, ymax, box = run.dp, run.n_p, run.y_max, 48.0
    q_all = np.arange(1, npp // 2)
    xis = q_all * dp
    nc, xw = 12, 3.0
    xc = -xw + (2 * xw / nc) * (np.arange(nc) + 0.5)
    kidx = np.clip(np.round((xc - run.r[0]) / run.dr).astype(int), 0,
                   len(run.r) - 1)
    k_res = run.k[kidx][:, q_all]

    def pairs(x, p):
        return same_row_pairs(x, p, dp, ymax, xw)

    def zsum(x, p, th, i, j, fixed=None):
        d = x[j] - x[i]
        xm = 0.5 * (x[i] + x[j])
        y = 0.5 * d
        mu = ((th[i] + p[i] * (xm - x[i]) / HBAR)
              - (th[j] + p[j] * (xm - x[j]) / HBAR))
        if fixed is None:           # Lagrangian bookkeeping: freeze at t0
            w = np.where(np.abs(y) < ymax,
                         np.cos(np.pi * y / (2 * ymax)) ** 2, 0.0)
            c = np.clip(((xm + xw) / (2 * xw / nc)).astype(int), 0, nc - 1)
        else:
            c, w = fixed
        u_trap = (eckart(x[j]) - eckart(x[i])
                  - y * (-eckart_force(x[i]) - eckart_force(x[j])))
        zr = np.zeros((nc, len(q_all)))
        pr = np.zeros((nc, len(q_all)))
        for iq, xq in enumerate(xis):
            ph = xq * d / HBAR + mu
            zr[:, iq] = np.bincount(c, w * np.cos(ph), minlength=nc)
            pr[:, iq] = np.bincount(c, u_trap * w * np.sin(ph), minlength=nc)
        return zr, pr, (c, w)

    corr = lambda a, b: float(np.corrcoef(a.ravel(), b.ravel())[0, 1])
    # meas = sum w (hbar dmu/dt) sin(.) is the pair-sum analogue of
    # int U w sin, which Theorem L1 multiplies by -B dp/hbar to make a
    # weight.  The dephasing kernel in L1's sign is therefore -meas, and its
    # correlation with K_res is -corr(meas, K_res) (erratum, section 7).
    print(f"   {'nu':>3} {'seed':>4} {'dt':>5} {'corr(meas, trapezoid pred)':>27}"
          f" {'corr(K_deph, K_res)':>20}")
    for nu in (8, 16):
        for seed in (1, 2, 3):
            rng = np.random.default_rng(seed)
            n = int(round(nu * B * box * npp * dp))
            x0 = rng.uniform(-box / 2, box / 2, n)
            p0 = dp * np.round(rng.uniform(-npp * dp / 2, npp * dp / 2, n) / dp)
            for dt_tot in (0.02, 0.1):
                x, p = x0.copy(), p0.copy()
                th = p * x / HBAR                           # re-lock
                i, j = pairs(x, p)
                z0, pred, cw = zsum(x, p, th, i, j)
                steps = int(round(dt_tot / 0.005))
                h = dt_tot / steps
                for _ in range(steps):
                    p += 0.5 * h * eckart_force(x)
                    x += p * h / MU
                    p += 0.5 * h * eckart_force(x)
                    th += (p ** 2 / (2 * MU) - eckart(x)) * h / HBAR
                z1, _, _ = zsum(x, p, th, i, j, cw)
                meas = -HBAR * (z1 - z0) / dt_tot
                print(f"   {nu:3d} {seed:4d} {dt_tot:5.2f}"
                      f" {corr(meas, pred):27.3f} {-corr(meas, k_res):20.3f}")
    # third moments of the two residual kernels
    y = np.linspace(-ymax, ymax, 40001)
    w = np.cos(np.pi * y / (2 * ymax)) ** 2

    def kern(x0, kind):
        u = eckart(x0 + y) - eckart(x0 - y)
        if kind == "mid":
            u_r = u + 2 * y * eckart_force(x0)
        else:
            u_r = u + y * (eckart_force(x0 - y) + eckart_force(x0 + y))
        return np.array([-(B * dp / HBAR)
                         * np.trapezoid(u_r * w * np.sin(2 * q * y / HBAR), y)
                         for q in xis])

    print(f"\n   {'x':>5} {'m3 midpoint (K_res)':>20} {'m3 trapezoid':>13}"
          f" {'ratio':>7} {'m1 trapezoid':>13}")
    for x0 in (-1.5, -1.0, -0.625, -0.3):
        km, kt = kern(x0, "mid"), kern(x0, "trap")
        m3m, m3t = 2 * (xis ** 3 * km).sum(), 2 * (xis ** 3 * kt).sum()
        print(f"   {x0:5.2f} {m3m:20.4f} {m3t:13.4f} {m3t / m3m:7.2f}"
              f" {2 * (xis * kt).sum():13.1e}")
    print("   leading order: midpoint (1/3) y^3 V''', trapezoid -(2/3) y^3 V''':"
          " ratio -2")
    print("   so the dephasing kernel is weakly ANTI-correlated with K_res")


# ----------------------------------------------------------------------
# C.  Figures
# ----------------------------------------------------------------------
XH, PH = 12.0, 3.0


def phase_field(pot, force, t, gas=False, nx=480, npx=240, dt=0.04):
    """theta(x, p, t): trace back to t = 0, add the Lagrangian action,
    start from the locked value p0 x0 / hbar (plus a fixed random field
    for the gas)."""
    x = np.linspace(-XH, XH, nx)[None, :].repeat(npx, 0)
    p = np.linspace(PH, -PH, npx)[:, None].repeat(nx, 1)
    s_act = np.zeros_like(x)
    for _ in range(int(round(t / dt))):
        p -= 0.5 * dt * force(x)
        x -= p * dt / MU
        p -= 0.5 * dt * force(x)
        s_act += (p ** 2 / (2 * MU) - pot(x)) * dt
    th = (p * x + s_act) / HBAR
    if gas:
        cell = np.floor(x * 4) * 127.1 + np.floor(p * 8) * 311.7
        th += 2 * np.pi * np.modf(np.abs(np.sin(cell) * 43758.5453))[0]
    return th


def rows_at(pot, force, t, n=300, dt=0.02,
            p_rows=(-2.5, -1.5, -0.5, 0.5, 1.5, 2.5)):
    out = []
    for p0 in p_rows:
        x = np.linspace(-XH, XH, n)
        p = np.full(n, p0)
        th = p0 * x / HBAR
        for _ in range(int(round(t / dt))):
            p += 0.5 * dt * force(x)
            x += p * dt / MU
            p += 0.5 * dt * force(x)
            th += (p ** 2 / (2 * MU) - pot(x)) * dt / HBAR
        out.append((x, p, th))
    return out


def draw_panel(ax, pot, force, t, gas=False, title=""):
    th = phase_field(pot, force, t, gas)
    hue = np.mod(th, 2 * np.pi) / (2 * np.pi)
    rgb = hsv_to_rgb(np.dstack([hue, np.full_like(hue, 0.6),
                                np.full_like(hue, 0.72)]))
    ax.imshow(rgb, extent=(-XH, XH, -PH, PH), aspect="auto",
              interpolation="bilinear")
    for x, p, th_r in rows_at(pot, force, t):
        if gas:
            th_r = th_r + np.random.default_rng(5).uniform(0, 2 * np.pi, len(x))
        pp = p + 0.11 * np.cos(th_r)
        brk = np.flatnonzero((np.abs(np.diff(p)) > 0.6)
                             | (np.abs(np.diff(x)) > 2))
        for seg in np.split(np.arange(len(x)), brk + 1):
            ok = np.abs(x[seg]) <= XH
            ax.plot(x[seg][ok], pp[seg][ok], color="white", lw=1.0)
    ax.set_xlim(-XH, XH)
    ax.set_ylim(-PH, PH)
    ax.set_title(title, fontsize=9)
    ax.tick_params(labelsize=7)


def part_c():
    banner("C. Figures")
    pots = (("Eckart barrier $V_0\\,\\mathrm{sech}^2 x$", eckart, eckart_force),
            ("soft-core Coulomb well $-1/\\sqrt{x^2+1}$", coulomb, coulomb_force))
    times = (0.0, 4.0, 10.0)
    fig, axes = plt.subplots(2, 3, figsize=(12, 6.4), sharex=True, sharey=True)
    for r, (name, pot, force) in enumerate(pots):
        for c, t in enumerate(times):
            draw_panel(axes[r, c], pot, force, t,
                       title=f"{name}, $t = {t:g}$")
    for ax in axes[1]:
        ax.set_xlabel("$x$")
    for ax in axes[:, 0]:
        ax.set_ylabel("$p$")
    fig.suptitle("Locked sea: phase field $\\theta = S/\\hbar$ over phase space"
                 " (hue), sample rows in white", fontsize=11)
    fig.tight_layout()
    name = "sea_lock_phase_field.png"
    fig.savefig(output_path(name), dpi=130, bbox_inches="tight")
    dp_ = docs_path(name)
    if dp_:
        fig.savefig(dp_, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"   wrote {name}")

    fig, axes = plt.subplots(1, 2, figsize=(10, 3.6), sharey=True)
    draw_panel(axes[0], eckart, eckart_force, 4.0, gas=False,
               title="locked sea, Eckart, $t = 4$")
    draw_panel(axes[1], eckart, eckart_force, 4.0, gas=True,
               title="gas sea, Eckart, $t = 4$")
    for ax in axes:
        ax.set_xlabel("$x$")
    axes[0].set_ylabel("$p$")
    fig.tight_layout()
    name = "sea_lock_locked_vs_gas.png"
    fig.savefig(output_path(name), dpi=130, bbox_inches="tight")
    dp_ = docs_path(name)
    if dp_:
        fig.savefig(dp_, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"   wrote {name}")


# ----------------------------------------------------------------------
# D.  Theorem L8: the reader's frame
# ----------------------------------------------------------------------
def part_d_symbolic():
    """The L8 identity, the two quadrature residuals, and the drift term."""
    import sympy as sp
    t, m, hb, y, s = sp.symbols("t m hbar y s", positive=True)
    V = sp.Function("V")
    xi, xj, pi_, pj, pr = (sp.Function(n)(t) for n in
                           ("x_i", "x_j", "p_i", "p_j", "p_ref"))
    th_dot = lambda xk, pk: (pk ** 2 / (2 * m) - V(xk)) / hb     # P1
    rules = {xi.diff(t): pi_ / m, xj.diff(t): pj / m}             # (S)
    mu_dot = (th_dot(xi, pi_) - th_dot(xj, pj)
              + (pr * (xj - xi)).diff(t).subs(rules) / hb)
    claim = (V(xj) - V(xi) + (xj - xi) * pr.diff(t)
             + (pj - pi_) * (2 * pr - pi_ - pj) / (2 * m)) / hb
    print(f"   identity: hbar dmu/dt - claim = "
          f"{sp.simplify(sp.expand(mu_dot - claim))}")
    # quadrature residuals and drift term, with V written through its
    # derivatives at the midpoint: V(x + z) = sum_k a_k z^k / k!
    a = sp.symbols("a0:9")                  # a_k = V^(k)(x)
    Vs = lambda z: sum(a[k] * z ** k / sp.factorial(k) for k in range(9))
    dVs = lambda z: sp.diff(Vs(s), s).subs(s, z)
    U = Vs(y) - Vs(-y)
    mid = sp.expand(U - 2 * y * dVs(0))
    trap = sp.expand(U - y * (dVs(-y) + dVs(y)))
    trunc = lambda e: sum(e.coeff(y, k) * y ** k for k in range(6))
    print(f"   midpoint residual  (compensated): {trunc(mid)}"
          "   [a_k = V^(k)(x)]")
    print(f"   trapezoid residual (own frames) : {trunc(trap)}")
    print(f"   ratio of leading terms: {sp.simplify(trap.coeff(y, 3) / mid.coeff(y, 3))},"
          f" of y^5 terms: {sp.simplify(trap.coeff(y, 5) / mid.coeff(y, 5))}")
    # partners at x -+ y re-locked on the reader's momentum; first order in t
    d_i, d_j, d_r = -dVs(-y) * t, -dVs(y) * t, -dVs(0) * t
    kin = sp.expand((d_j - d_i) * (2 * d_r - d_i - d_j) / (2 * m))
    print(f"   drift term / U_res as y -> 0: "
          f"{sp.simplify(kin.coeff(y, 3) / mid.coeff(y, 3))}")


def part_d():
    banner("D. Theorem L8: the reader's frame")
    part_d_symbolic()
    run = Ledger(v0=1.0, a=1.0, n_r=192, r_half=24.0, n_p=64, dp=0.25)
    dp, npp, ymax, box = run.dp, run.n_p, run.y_max, 48.0
    veff = lambda z: np.interp(z, run.r, run.dv_eff)    # parent's +V'_eff
    q_all = np.arange(1, npp // 2)
    xis = q_all * dp
    nc, xw = 12, 3.0
    xc = -xw + (2 * xw / nc) * (np.arange(nc) + 0.5)
    kidx = np.clip(np.round((xc - run.r[0]) / run.dr).astype(int), 0,
                   len(run.r) - 1)
    k_res = run.k[kidx][:, q_all]
    corr = lambda a, b: float(np.corrcoef(a.ravel(), b.ravel())[0, 1])

    def sea(nu, seed):
        rng = np.random.default_rng(seed)
        n = int(round(nu * B * box * npp * dp))
        x0 = rng.uniform(-box / 2, box / 2, n)
        p0 = dp * np.round(rng.uniform(-npp * dp / 2, npp * dp / 2, n) / dp)
        return x0, p0

    def mu(frame, x, p, th, i, j, pref):
        if frame == "own":                    # each partner's own lever
            return th[i] - th[j] + 0.5 * (p[i] + p[j]) * (x[j] - x[i]) / HBAR
        return th[i] - th[j] + pref * (x[j] - x[i]) / HBAR

    def evolve(x, p, th, xr, pr, T, h, frame):
        """Velocity Verlet under (S); clocks by the trapezoid rule in time.
        The parent-frame reader streams under -V'_eff, the inertial reader
        not at all."""
        lag = lambda xx, pp: (pp ** 2 / (2 * MU) - eckart(xx)) / HBAR
        for _ in range(int(round(T / h))):
            l0 = lag(x, p)
            p += 0.5 * h * eckart_force(x)
            x += p * h / MU
            p += 0.5 * h * eckart_force(x)
            th += 0.5 * h * (l0 + lag(x, p))
            if frame == "parent":
                pr -= 0.5 * h * veff(xr)
                xr += pr * h / MU
                pr -= 0.5 * h * veff(xr)
        return x, p, th, xr, pr

    # -- D2: the three readers at re-lock -------------------------------
    print(f"\n   kernel read at re-lock (dt = 0.02); K_deph in Theorem L1's sign")
    print(f"   {'nu':>3} {'seed':>4} {'reader':>9} {'corr(K_deph, K_res)':>20}"
          f" {'slope on ideal U_res pair sum':>30}")
    for nu in (8, 16):
        for seed in (1, 2, 3):
            x0, p0 = sea(nu, seed)
            i, j = same_row_pairs(x0, p0, dp, ymax, xw)
            xm0, y0 = 0.5 * (x0[i] + x0[j]), 0.5 * (x0[j] - x0[i])
            c = np.clip(((xm0 + xw) / (2 * xw / nc)).astype(int), 0, nc - 1)
            w = np.cos(np.pi * y0 / (2 * ymax)) ** 2
            sn = np.sin(np.outer(2 * y0 / HBAR, xis))       # mu = 0 weights
            u_res = eckart(x0[j]) - eckart(x0[i]) - 2 * y0 * veff(xm0)
            ideal = np.stack([np.bincount(c, u_res * w * sn[:, q], minlength=nc)
                              for q in range(len(xis))], axis=1)
            for frame in ("inertial", "parent", "own"):
                x, p = x0.copy(), p0.copy()
                th = p * x / HBAR                           # re-lock
                xr, pr = xm0.copy(), p0[i].copy()
                m0 = mu(frame, x, p, th, i, j, pr)
                x, p, th, xr, pr = evolve(x, p, th, xr, pr, 0.02, 0.002, frame)
                rate = HBAR * (mu(frame, x, p, th, i, j, pr) - m0) / 0.02
                meas = np.stack([np.bincount(c, rate * w * sn[:, q],
                                             minlength=nc)
                                 for q in range(len(xis))], axis=1)
                slope = float((meas * ideal).sum() / (ideal * ideal).sum())
                print(f"   {nu:3d} {seed:4d} {frame:>9}"
                      f" {-corr(meas, k_res):20.3f} {slope:30.3f}")

    # -- D3: the parent-frame reading after re-lock ---------------------
    x0, p0 = sea(8, 1)
    i, j = same_row_pairs(x0, p0, dp, ymax, xw)
    xm0, y0 = 0.5 * (x0[i] + x0[j]), 0.5 * (x0[j] - x0[i])
    bands = (("0.3 < y < 1", (y0 > 0.3) & (y0 < 1)),
             ("1 < y < y_max", y0 >= 1))
    xg = np.linspace(-9, 9, 1801)
    vpp = np.abs(np.interp(xm0, xg, np.gradient(veff(xg), xg)))
    h = 0.001
    ts = (0.0, 0.1, 0.2, 0.4, 0.7, 1.0)
    rows = {b: {"own": [], "parent": [], "resid": [], "pred": []}
            for b, _ in bands}
    for T in ts:
        got = {}
        for frame in ("own", "parent"):
            st = (x0.copy(), p0.copy(), p0 * x0 / HBAR, xm0.copy(),
                  p0[i].copy())
            a = evolve(*st, T, h, frame) if T > 0 else st
            b = evolve(*[v.copy() for v in a], 2 * h, h, frame)
            rate = HBAR * (mu(frame, *b[:3], i, j, b[4])
                           - mu(frame, *a[:3], i, j, a[4])) / (2 * h)

            def target(s_):
                x, p, _, xr, pr = s_
                if frame == "own":            # target is still U_res
                    fr = veff(0.5 * (x[i] + x[j]))
                else:
                    fr = veff(xr)
                return (eckart(x[j]) - eckart(x[i]) - (x[j] - x[i]) * fr,
                        (p[j] - p[i]) * (2 * pr - p[i] - p[j]) / (2 * MU))
            (ia, ka), (ib, kb) = target(a), target(b)
            got[frame] = (rate, 0.5 * (ia + ib), 0.5 * (ka + kb))
        for bn, m in bands:
            rel = lambda e, u: float(np.sqrt((e[m] ** 2).mean()
                                             / (u[m] ** 2).mean()))
            ro, io, _ = got["own"]
            rp, ip, kp = got["parent"]
            rows[bn]["own"].append(rel(ro - io, io))
            rows[bn]["parent"].append(rel(rp - ip, ip))
            rows[bn]["resid"].append(rel(rp - ip - kp, ip))
            rows[bn]["pred"].append(3 * vpp[m].mean() * T ** 2)
    print(f"\n   rms error of hbar dmu/dt against U_res, relative to rms U_res"
          f" (nu = 8, seed 1)")
    print(f"   {'T':>5} {'band':>14} {'own frame':>10} {'parent':>9}"
          f" {'parent - drift term':>20} {'3<|V\"|> T^2':>12}")
    for k, T in enumerate(ts):
        for bn, _ in bands:
            r = rows[bn]
            print(f"   {T:5.2f} {bn:>14} {r['own'][k]:10.3f} {r['parent'][k]:9.4f}"
                  f" {r['resid'][k]:20.1e} {r['pred'][k]:12.3f}")

    # -- D4: figure -----------------------------------------------------
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11, 4.0))
    xa, ya = -0.3, 0.7
    xx = np.linspace(xa - ya, xa + ya, 400)
    vp = lambda z: -eckart_force(z)                     # V'(x)
    chord = vp(xa - ya) + (vp(xa + ya) - vp(xa - ya)) * (xx - xx[0]) / (2 * ya)
    ax0.fill_between(xx, vp(xa), vp(xx), color="#1D9E75", alpha=0.35,
                     label="parent frame: midpoint residual")
    ax0.fill_between(xx, vp(xx), chord, color="#D85A30", alpha=0.35,
                     label="own frames: trapezoid residual")
    ax0.plot(xx, vp(xx), color="k", lw=1.4, label="$V'(x')$ on the chord")
    ax0.axhline(vp(xa), xmin=0, xmax=1, color="#1D9E75", lw=1, ls="--")
    ax0.plot(xx, chord, color="#D85A30", lw=1, ls="--")
    for z, lab in ((xa - ya, "$x-y$"), (xa, "$x$"), (xa + ya, "$x+y$")):
        ax0.axvline(z, color="0.6", lw=0.6)
        ax0.text(z, ax0.get_ylim()[0], lab, ha="center", va="bottom",
                 fontsize=8)
    ax0.set_xlim(xx[0] - 0.1, xx[-1] + 0.1)
    u_c = eckart(xa + ya) - eckart(xa - ya)
    r_mid = u_c - 2 * ya * vp(xa)
    r_trap = u_c - ya * (vp(xa - ya) + vp(xa + ya))
    ax0.text(0.03, 0.10, f"midpoint residual {r_mid:+.4f}\n"
             f"trapezoid residual {r_trap:+.4f}\nratio {r_trap / r_mid:+.2f}",
             transform=ax0.transAxes, fontsize=8, va="bottom")
    print(f"   figure chord x = {xa}, y = {ya}: midpoint {r_mid:+.5f},"
          f" trapezoid {r_trap:+.5f}, ratio {r_trap / r_mid:+.3f}")
    for yy in (0.2, 0.4):
        u_y = eckart(xa + yy) - eckart(xa - yy)
        print(f"   same x, y = {yy}: ratio trapezoid/midpoint"
              f" {(u_y - yy * (vp(xa - yy) + vp(xa + yy))) / (u_y - 2 * yy * vp(xa)):+.3f}")
    ax0.set_xlabel("$x'$")
    ax0.set_title(f"$U=\\int V'\\,dx'$ over the chord, Eckart, $x={xa}$,"
                  f" $y={ya}$: each residual is a signed area", fontsize=9)
    ax0.legend(fontsize=7, loc="upper right")
    tt = np.array(ts[1:])
    tex = lambda b: b.replace("y_max", "y_{max}")
    for (bn, _), mfc in zip(bands, ("full", "none")):
        r = rows[bn]
        ax1.loglog(tt, r["parent"][1:], "o-", color="#1D9E75", fillstyle=mfc,
                   label=f"parent frame, ${tex(bn)}$")
        ax1.loglog(tt, r["own"][1:], "s:", color="#D85A30", fillstyle=mfc,
                   label=f"own frames, ${tex(bn)}$")
    ax1.loglog(tt, rows[bands[0][0]]["pred"][1:], "k--", lw=0.8,
               label="$3\\langle|V''|\\rangle T^2$")
    ax1.set_xlabel("time $T$ since re-lock")
    ax1.set_ylabel("rms error / rms $U_{\\mathrm{res}}$")
    ax1.set_title("the reading against the compensated kernel", fontsize=9)
    ax1.legend(fontsize=7)
    fig.tight_layout()
    name = "sea_phase_reader_frames.png"
    fig.savefig(output_path(name), dpi=130, bbox_inches="tight")
    dp_ = docs_path(name)
    if dp_:
        fig.savefig(dp_, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"   wrote {name}")


if __name__ == "__main__":
    part_a()
    part_b()
    part_c()
    part_d()
