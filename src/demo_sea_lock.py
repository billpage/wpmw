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
    print(f"   {'nu':>3} {'seed':>4} {'dt':>5} {'corr(meas, trapezoid pred)':>27}"
          f" {'corr(meas, K_res)':>18}")
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
                      f" {corr(meas, pred):27.3f} {corr(meas, k_res):18.3f}")
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


if __name__ == "__main__":
    part_a()
    part_b()
    part_c()
