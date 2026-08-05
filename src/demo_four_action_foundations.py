"""Foundations and uniqueness of the Four-Action Wigner Particle Model.

Numerical companion to ``docs/supplement/four_action_foundations.md``.

The questions asked here are the ones raised in D. Cyganski's note *A journey
from Bohm trajectory theory, through Nelson's SDEs and Wigner Particles to the
Closed Four Action Model* (3 August 2026), §"Does momentum and energy balance
fully determine the FAWPM?".

Parts
-----
A. Every member of the exact family ``F = (A - A^-1) G``,
   ``H = (2 - A - A^-1) G - Gamma`` reproduces the QLE stencil *and*
   conserves particle number, momentum and energy exactly.  Conservation is
   therefore a *consequence* of exactness, not an independent constraint.

B. The endpoint-local two-parameter ansatz
   ``r_n = a Gamma (W_{n+1} - W_{n-1})``, ``lambda_n = b Gamma (W_{n+1} + W_{n-1})``
   has mean-field generator ``2a Gamma D1 - (a+b) Gamma D2``.  Momentum and
   energy balance fix ``b = -1/2`` for *every* ``a``: conservation cannot
   select the model.  Only the closure condition ``a + b = 0`` does.

C. Where the ``a`` dependence lives: moments 0, 1 and 2 of the generator are
   ``a``-independent; the dependence first appears at moment 3.  Conservation
   laws are moment conditions and are blind to the quantum corrections.

D. Free momentum-grid step ``delta``: the model is Moyal evolution with
   ``hbar_eff = 2 delta / k``.  Classical Liouville as ``delta -> 0``; the true
   QLE iff ``delta = hbar k / 2`` (half the photon momentum of the mode).

E. Endpoint locality selects the model: within the exact family, requiring
   that neither rate reference the centre cell forces ``G = Gamma/2``.

F. The four-wave-mixing reading: linearising a two-body mass-action rate about
   a uniform sea yields an endpoint-*symmetric* rate, never the
   endpoint-antisymmetric focus rate.

Run:  WPMW_OUTPUT=... python src/demo_four_action_foundations.py
"""

from __future__ import annotations

import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from wpmwlib.wpmw_utils import docs_path, output_path  # noqa: E402

# ---------------------------------------------------------------------------
# Lattice / physical parameters
# ---------------------------------------------------------------------------

HBAR = 1.0
MASS = 1.0
L = 8.0
V_P = 1.5
KWAVE = 2.0 * np.pi / L          # fundamental mode wavenumber
DELTA_HBAR = HBAR * KWAVE / 2.0  # = pi hbar / L, half the photon momentum

NCELL = 128                      # momentum cells
QMODE = 1                        # Fourier mode index

rng = np.random.default_rng(20260805)


def shift(w, j):
    """(A^j w)_n = w_{n+j}."""
    return np.roll(w, -j)


def generator(f, h, q=QMODE):
    """Mean-field generator of the four actions, spec eq. (2).

    dW_n = 2 f_n - f_{n-q} - f_{n+q} + h_{n-q} - h_{n+q}
    """
    return 2.0 * f - shift(f, -q) - shift(f, q) + shift(h, -q) - shift(h, q)


def qle_stencil(w, gamma, q=QMODE):
    """Target: dW_n = Gamma (W_{n+q} - W_{n-q})."""
    return gamma * (shift(w, q) - shift(w, -q))


def apply_operator(w, stencil, q=QMODE):
    """stencil is {j: c_j} meaning sum_j c_j A^{jq}."""
    out = np.zeros_like(w)
    for j, c in stencil.items():
        out = out + c * shift(w, j * q)
    return out


def compact_field(n=NCELL, width=18, centre=None):
    """Random field with compact support, padded away from the wrap."""
    w = np.zeros(n)
    c = n // 2 if centre is None else centre
    w[c - width : c + width] = rng.normal(size=2 * width)
    return w


# ---------------------------------------------------------------------------
# Part A — exact family: exactness implies conservation
# ---------------------------------------------------------------------------


def part_a():
    print("=" * 74)
    print("PART A  Exact family (4): exactness => number, momentum, energy")
    print("=" * 74)
    delta = DELTA_HBAR
    p = (np.arange(NCELL) - NCELL // 2) * delta

    print(f"{'G stencil range':>16} {'|gen - QLE|':>13} {'dN':>12} "
          f"{'dP resid':>12} {'dE resid':>12}")

    worst = 0.0
    for trial in range(6):
        x = rng.uniform(0.3, 3.0)
        gamma = -V_P * KWAVE / (2.0 * delta) * np.sin(KWAVE * x)
        dUdx = -V_P * KWAVE * np.sin(KWAVE * x)

        # random finite-range translation-invariant G
        jr = trial % 3            # half-range 0, 1, 2
        gs = {j: rng.normal() for j in range(-jr, jr + 1)}

        w = compact_field()
        rho = w.sum()
        pmom = float(np.sum(p * w))

        # F = (A - A^-1) G ; H = (2 - A - A^-1) G - Gamma I
        gw = apply_operator(w, gs)
        f = shift(gw, 1) - shift(gw, -1)
        h = 2.0 * gw - shift(gw, 1) - shift(gw, -1) - gamma * w

        dw = generator(f, h)
        err = np.max(np.abs(dw - qle_stencil(w, gamma)))

        dN = float(dw.sum())
        dP = float(np.sum(p * dw))
        dE = float(np.sum(p**2 / (2 * MASS) * dw))

        dP_expected = -dUdx * rho
        dE_expected = -dUdx * pmom / MASS

        scaleP = max(abs(dP_expected), 1e-30)
        scaleE = max(abs(dE_expected), 1e-30)
        rP = abs(dP - dP_expected) / scaleP
        rE = abs(dE - dE_expected) / scaleE
        worst = max(worst, err, abs(dN), rP, rE)

        print(f"{'+/- ' + str(jr):>16} {err:13.3e} {abs(dN):12.3e} "
              f"{rP:12.3e} {rE:12.3e}")

    print(f"\nworst deviation over all trials: {worst:.3e}")
    print("Conservation holds for EVERY member of the family, i.e. it is")
    print("implied by QLE-exactness and cannot discriminate between members.\n")
    return worst


# ---------------------------------------------------------------------------
# Part B — the (a, b) ansatz
# ---------------------------------------------------------------------------


def rates_ab(w, gamma, a, b, q=QMODE):
    f = a * gamma * (shift(w, q) - shift(w, -q))
    h = b * gamma * (shift(w, q) + shift(w, -q))
    return f, h


def part_b():
    print("=" * 74)
    print("PART B  Endpoint-local ansatz: what conservation does and does not fix")
    print("=" * 74)
    delta = DELTA_HBAR
    p = (np.arange(NCELL) - NCELL // 2) * delta
    x = 1.1
    gamma = -V_P * KWAVE / (2.0 * delta) * np.sin(KWAVE * x)
    dUdx = -V_P * KWAVE * np.sin(KWAVE * x)

    w = compact_field()
    rho = w.sum()
    pmom = float(np.sum(p * w))

    # B1: closed form of the generator
    print("B1  generator identity  dW = 2a*Gamma*D1 W - (a+b)*Gamma*D2 W")
    worst = 0.0
    for a, b in [(0.5, -0.5), (0.0, -0.5), (1.3, 0.7), (-0.2, 0.4)]:
        f, h = rates_ab(w, gamma, a, b)
        dw = generator(f, h)
        d1 = shift(w, 1) - shift(w, -1)
        d2 = shift(w, 2) - shift(w, -2)
        pred = 2 * a * gamma * d1 - (a + b) * gamma * d2
        e = np.max(np.abs(dw - pred))
        worst = max(worst, e)
        print(f"    a={a:+.2f} b={b:+.2f}   max|dW - predicted| = {e:.3e}")

    # B2: conservation residuals across the (a, b) plane
    print("\nB2  conservation residuals  (relative)")
    print(f"{'a':>7} {'b':>7} {'dN':>12} {'dP resid':>12} {'dE resid':>12} "
          f"{'QLE resid':>12}")
    for a, b in [(0.5, -0.5), (0.0, -0.5), (2.0, -0.5), (-1.0, -0.5),
                 (0.5, 0.0), (0.5, -1.0)]:
        f, h = rates_ab(w, gamma, a, b)
        dw = generator(f, h)
        dN = float(dw.sum())
        dP = float(np.sum(p * dw))
        dE = float(np.sum(p**2 / (2 * MASS) * dw))
        rP = abs(dP + dUdx * rho) / abs(dUdx * rho)
        rE = abs(dE + dUdx * pmom / MASS) / abs(dUdx * pmom / MASS)
        rq = (np.max(np.abs(dw - qle_stencil(w, gamma)))
              / np.max(np.abs(qle_stencil(w, gamma))))
        print(f"{a:7.2f} {b:7.2f} {abs(dN):12.3e} {rP:12.3e} {rE:12.3e} "
              f"{rq:12.3e}")

    print("\n  b = -1/2 gives exact momentum AND energy balance for every a.")
    print("  a = 0, b = -1/2 is a fully conserving, Ehrenfest-exact model that")
    print("  is NOT the QLE: conservation alone does not determine the FAWPM.\n")
    return worst


# ---------------------------------------------------------------------------
# Part C — moments: where the a-dependence hides
# ---------------------------------------------------------------------------


def part_c():
    print("=" * 74)
    print("PART C  Moment analysis of the conserving family (b = -1/2)")
    print("=" * 74)
    delta = DELTA_HBAR
    p = (np.arange(NCELL) - NCELL // 2) * delta
    x = 1.1
    gamma = -V_P * KWAVE / (2.0 * delta) * np.sin(KWAVE * x)
    w = compact_field()

    print(f"{'a':>7} " + " ".join(f"{'M' + str(k):>13}" for k in range(5)))
    rows = []
    for a in [0.0, 0.25, 0.5, 1.0, 2.0]:
        f, h = rates_ab(w, gamma, a, -0.5)
        dw = generator(f, h)
        ms = [float(np.sum(p**k * dw)) for k in range(5)]
        sc = [float(np.sum(np.abs(p**k * dw))) for k in range(5)]
        rows.append((a, ms, sc))
        print(f"{a:7.2f} " + " ".join(f"{m:13.5e}" for m in ms))

    base, bscale = rows[2][1], rows[2][2]  # a = 0.5, the QLE member
    print("\nspread over a, normalised by the moment's own absolute mass:")
    for k in range(5):
        vals = [r[1][k] for r in rows]
        sc = max(bscale[k], 1e-30)
        print(f"   moment {k}: max deviation = "
              f"{max(abs(v - base[k]) for v in vals) / sc:.3e}")
    print("\n  Moments 0, 1, 2 (number, momentum, energy) are a-independent.")
    print("  The a-dependence first appears at moment 3 — beyond the reach of")
    print("  any conservation law.\n")
    return rows


# ---------------------------------------------------------------------------
# Part D — free grid step: hbar_eff
# ---------------------------------------------------------------------------


def hermite_e(m, z):
    """Probabilists' Hermite polynomial He_m evaluated at z."""
    h0 = np.ones_like(z)
    if m == 0:
        return h0
    h1 = z.copy()
    for n in range(1, m):
        h0, h1 = h1, z * h1 - n * h0
    return h1


def gaussian_dp(p, s, m):
    """d^m/dp^m exp(-p^2/(2 s^2)), exactly."""
    z = p / s
    return (-1.0) ** m * s ** (-m) * hermite_e(m, z) * np.exp(-z**2 / 2.0)


def moyal_rhs(p, s, x, hbar_eff, nmax=24):
    """Moyal collision series for U = V_P cos(k x) with a given hbar."""
    total = np.zeros_like(p)
    for n in range(nmax + 1):
        # d^{2n+1}/dx^{2n+1} [V_P cos(kx)] = V_P k^{2n+1} (-1)^{n+1} sin(kx)
        dU = V_P * KWAVE ** (2 * n + 1) * (-1) ** (n + 1) * np.sin(KWAVE * x)
        coef = (-1.0) ** n * (hbar_eff / 2.0) ** (2 * n) / math.factorial(2 * n + 1)
        total = total + coef * dU * gaussian_dp(p, s, 2 * n + 1)
    return total


def part_d():
    print("=" * 74)
    print("PART D  Free momentum-grid step: the model is Moyal with hbar_eff")
    print("=" * 74)
    p = np.linspace(-6.0, 6.0, 601)
    s = 1.0
    x = 1.1

    print(f"{'delta/delta_hbar':>17} {'hbar_eff':>10} {'vs Moyal(hbar_eff)':>20} "
          f"{'vs QLE(hbar)':>14} {'vs classical':>14}")

    ratios = np.array([0.05, 0.25, 0.5, 1.0, 1.5, 2.0])
    dev_true, dev_cls, dev_eff = [], [], []
    for ratio in ratios:
        delta = ratio * DELTA_HBAR
        gamma = -V_P * KWAVE / (2.0 * delta) * np.sin(KWAVE * x)
        # stencil generator, exact shift of the Gaussian
        stencil = gamma * (np.exp(-((p + delta) ** 2) / (2 * s**2))
                           - np.exp(-((p - delta) ** 2) / (2 * s**2)))
        hbar_eff = 2.0 * delta / KWAVE
        m_eff = moyal_rhs(p, s, x, hbar_eff)
        m_true = moyal_rhs(p, s, x, HBAR)
        dUdx = -V_P * KWAVE * np.sin(KWAVE * x)
        m_cls = dUdx * gaussian_dp(p, s, 1)

        sc = np.max(np.abs(m_true))
        e_eff = np.max(np.abs(stencil - m_eff)) / sc
        e_true = np.max(np.abs(stencil - m_true)) / sc
        e_cls = np.max(np.abs(stencil - m_cls)) / sc
        dev_eff.append(e_eff)
        dev_true.append(e_true)
        dev_cls.append(e_cls)
        print(f"{ratio:17.2f} {hbar_eff:10.4f} {e_eff:20.3e} {e_true:14.3e} "
              f"{e_cls:14.3e}")

    print("\n  The four-action generator is EXACTLY Moyal evolution with")
    print("  hbar_eff = 2 delta / k.  It is the classical Liouville equation")
    print("  as delta -> 0 and the true QLE only at delta = hbar k / 2.\n")
    return ratios, np.array(dev_eff), np.array(dev_true), np.array(dev_cls)


# ---------------------------------------------------------------------------
# Part E — endpoint locality selects G = Gamma/2
# ---------------------------------------------------------------------------


def part_e():
    print("=" * 74)
    print("PART E  Endpoint locality selects the symmetric member")
    print("=" * 74)
    gamma = 0.7314

    def stencils(gs):
        """F and H stencils (as {j: c}) for a given G = sum_j c_j A^j."""
        fs, hs = {}, {}
        for j, c in gs.items():
            fs[j + 1] = fs.get(j + 1, 0.0) + c
            fs[j - 1] = fs.get(j - 1, 0.0) - c
            hs[j] = hs.get(j, 0.0) + 2.0 * c
            hs[j + 1] = hs.get(j + 1, 0.0) - c
            hs[j - 1] = hs.get(j - 1, 0.0) - c
        hs[0] = hs.get(0, 0.0) - gamma
        return fs, hs

    def endpoint_local(st):
        """True if the stencil touches only cells n+1 and n-1."""
        return all(abs(c) < 1e-13 for j, c in st.items() if j not in (1, -1))

    trials = [("G = 0  (pure hop / original single rule)", {0: 0.0}),
              ("G = Gamma/2", {0: gamma / 2}),
              ("G = Gamma/3", {0: gamma / 3}),
              ("G = Gamma/2 + 0.4 A", {0: gamma / 2, 1: 0.4}),
              ("G = 0.3 A^-1", {-1: 0.3})]
    for name, gs in trials:
        fs, hs = stencils(gs)
        ok = endpoint_local(fs) and endpoint_local(hs)
        fc = {j: round(c, 6) for j, c in sorted(fs.items()) if abs(c) > 1e-13}
        hc = {j: round(c, 6) for j, c in sorted(hs.items()) if abs(c) > 1e-13}
        print(f"  {name:38s} local={str(ok):5s}")
        print(f"      F stencil {fc}")
        print(f"      H stencil {hc}")

    print("\n  Only G = Gamma/2 leaves both rates supported on {n-1, n+1}.")
    print("  'The intermediate state density no longer weights the actions'")
    print("  is therefore the selection principle, and it works alone.\n")


# ---------------------------------------------------------------------------
# Part F — the four-wave-mixing reading
# ---------------------------------------------------------------------------


def part_f():
    print("=" * 74)
    print("PART F  Mass-action linearisation about a uniform sea")
    print("=" * 74)
    B = 1.0e3
    eps = 1.0e-4
    c1, c2 = 0.9, 1.4
    n = 64
    w = B * np.ones(n)
    tilde = np.zeros(n)
    tilde[30] = 1.0                       # probe cell n-1 of centre 31
    probe_lo = eps * tilde
    probe_hi = eps * np.roll(tilde, 2)    # probe cell n+1 of centre 31

    def net_rate(field, centre=31):
        """c1 W_{k-1} W_{k+1} - c2 W_k^2, a general two-body net rate."""
        return (c1 * field[centre - 1] * field[centre + 1]
                - c2 * field[centre] ** 2)

    r0 = net_rate(w)
    d_lo = (net_rate(w + probe_lo) - r0) / eps
    d_hi = (net_rate(w + probe_hi) - r0) / eps
    print(f"  d(rate)/d W_(k-1) = {d_lo:.6e}")
    print(f"  d(rate)/d W_(k+1) = {d_hi:.6e}")
    print(f"  symmetric part    = {(d_hi + d_lo) / 2:.6e}")
    print(f"  antisymmetric     = {(d_hi - d_lo) / 2:.6e}")
    print("\n  The linearisation is endpoint-SYMMETRIC for any c1, c2, whereas")
    print("  the focus rate r is endpoint-ANTISYMMETRIC.  A two-body product")
    print("  rate linearised about a uniform sea cannot produce r.\n")


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------


def figure_uniqueness_map():
    delta = DELTA_HBAR
    p = (np.arange(NCELL) - NCELL // 2) * delta
    x = 1.1
    gamma = -V_P * KWAVE / (2.0 * delta) * np.sin(KWAVE * x)
    dUdx = -V_P * KWAVE * np.sin(KWAVE * x)
    w = compact_field()
    rho = w.sum()
    pmom = float(np.sum(p * w))
    target = qle_stencil(w, gamma)
    tscale = np.max(np.abs(target))

    avals = np.linspace(-0.5, 1.5, 161)
    bvals = np.linspace(-1.5, 0.5, 161)
    Q = np.zeros((bvals.size, avals.size))
    C = np.zeros_like(Q)
    for i, b in enumerate(bvals):
        for j, a in enumerate(avals):
            f, h = rates_ab(w, gamma, a, b)
            dw = generator(f, h)
            Q[i, j] = np.max(np.abs(dw - target)) / tscale
            dP = float(np.sum(p * dw))
            dE = float(np.sum(p**2 / (2 * MASS) * dw))
            C[i, j] = (abs(dP + dUdx * rho) / abs(dUdx * rho)
                       + abs(dE + dUdx * pmom / MASS)
                       / abs(dUdx * pmom / MASS))

    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.6))

    ax = axes[0]
    ext = [avals[0], avals[-1], bvals[0], bvals[-1]]
    im = ax.imshow(np.log10(Q + 1e-16), origin="lower", extent=ext,
                   aspect="auto", cmap="magma_r", vmin=-3.0, vmax=0.8)
    ax.plot(avals, -0.5 * np.ones_like(avals), "-", color="tab:cyan", lw=2.2,
            label=r"momentum + energy balance:  $b=-1/2$")
    ax.plot(avals, -avals, "-", color="lime", lw=2.2,
            label=r"single-harmonic closure:  $a+b=0$")
    ax.plot([0.5], [-0.5], "o", ms=11, mfc="none", mec="w", mew=2.4)
    ax.annotate("FAWPM", (0.5, -0.5), textcoords="offset points",
                xytext=(12, 12), color="w", fontsize=11, weight="bold")
    ax.set_xlabel(r"$a$   (focus-channel gain)")
    ax.set_ylabel(r"$b$   (hop-channel gain)")
    ax.set_title("QLE residual over the rate plane", fontsize=11)
    ax.set_ylim(bvals[0], bvals[-1])
    ax.legend(loc="lower left", fontsize=8, framealpha=0.9)
    fig.colorbar(im, ax=ax, label=r"$\log_{10}$ relative residual")

    # 1-D cut along the conservation line b = -1/2
    ax = axes[1]
    acut = np.linspace(-0.5, 1.5, 401)
    qcut, ccut = [], []
    for a in acut:
        f, h = rates_ab(w, gamma, a, -0.5)
        dw = generator(f, h)
        qcut.append(np.max(np.abs(dw - target)) / tscale)
        dP = float(np.sum(p * dw))
        dE = float(np.sum(p**2 / (2 * MASS) * dw))
        ccut.append(abs(dP + dUdx * rho) / abs(dUdx * rho)
                    + abs(dE + dUdx * pmom / MASS)
                    / abs(dUdx * pmom / MASS))
    ax.semilogy(acut, np.array(qcut) + 1e-17, lw=2.2, color="tab:red",
                label="QLE residual")
    ax.semilogy(acut, np.array(ccut) + 1e-17, lw=2.2, color="tab:cyan",
                label="momentum + energy residual")
    ax.axvline(0.5, color="k", lw=1, ls=":")
    ax.annotate(r"$a=1/2$", (0.5, 1e-6), textcoords="offset points",
                xytext=(8, 0), fontsize=10)
    ax.set_xlabel(r"$a$   (focus-channel gain), along $b=-1/2$")
    ax.set_ylabel("relative residual")
    ax.set_ylim(1e-17, 5.0)
    ax.set_title("Along the conserving line, conservation says nothing",
                 fontsize=11)
    ax.legend(fontsize=9, loc="center right")
    ax.grid(alpha=0.3)

    fig.suptitle("Conservation constrains a line; only closure picks the point",
                 fontsize=12.5)
    fig.tight_layout()
    for pth in (output_path("four_action_uniqueness_map.png"),
                docs_path("four_action_uniqueness_map.png")):
        if pth:
            fig.savefig(pth, dpi=150, bbox_inches="tight")
    plt.close(fig)


def figure_hbar_eff(moment_rows):
    p = np.linspace(-6.0, 6.0, 601)
    s = 1.0
    x = 1.1
    ratios = np.linspace(0.05, 2.0, 79)  # includes delta = hbar k / 2 exactly
    dev_true, dev_eff, dev_cls = [], [], []
    m_true = moyal_rhs(p, s, x, HBAR)
    dUdx = -V_P * KWAVE * np.sin(KWAVE * x)
    m_cls = dUdx * gaussian_dp(p, s, 1)
    sc = np.max(np.abs(m_true))
    for ratio in ratios:
        delta = ratio * DELTA_HBAR
        gamma = -V_P * KWAVE / (2.0 * delta) * np.sin(KWAVE * x)
        stencil = gamma * (np.exp(-((p + delta) ** 2) / (2 * s**2))
                           - np.exp(-((p - delta) ** 2) / (2 * s**2)))
        dev_eff.append(np.max(np.abs(stencil - moyal_rhs(
            p, s, x, 2.0 * delta / KWAVE))) / sc)
        dev_true.append(np.max(np.abs(stencil - m_true)) / sc)
        dev_cls.append(np.max(np.abs(stencil - m_cls)) / sc)

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4))
    ax = axes[0]
    ax.semilogy(ratios, np.array(dev_eff) + 1e-17, lw=2,
                label=r"vs Moyal with $\hbar_{\rm eff}=2\delta/k$")
    ax.semilogy(ratios, dev_true, lw=2, label=r"vs true QLE ($\hbar$)")
    ax.semilogy(ratios, dev_cls, lw=2, ls="--",
                label="vs classical Liouville")
    ax.axvline(1.0, color="k", lw=1, ls=":")
    ax.annotate(r"$\delta=\hbar k/2$", (1.0, 1e-3), rotation=90,
                textcoords="offset points", xytext=(6, 0), fontsize=9)
    ax.set_xlabel(r"$\delta \;/\; (\hbar k/2)$")
    ax.set_ylabel("relative deviation")
    ax.set_title("The grid step is the whole quantum input", fontsize=11)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[1]
    ks = np.arange(5)
    base, bscale = moment_rows[2][1], moment_rows[2][2]
    for a, ms, _ in moment_rows:
        rel = [abs(ms[k] - base[k]) / max(bscale[k], 1e-30) for k in ks]
        ax.semilogy(ks, np.array(rel) + 1e-17, "o-", label=f"a = {a}")
    ax.set_xticks(ks)
    ax.set_xlabel(r"moment order $\int p^k \,\dot W\, dp$")
    ax.set_ylabel(r"relative deviation from the $a=1/2$ member")
    ax.set_title("Conservation laws are blind beyond moment 2", fontsize=11)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    fig.tight_layout()
    for pth in (output_path("four_action_hbar_effective.png"),
                docs_path("four_action_hbar_effective.png")):
        if pth:
            fig.savefig(pth, dpi=150, bbox_inches="tight")
    plt.close(fig)


def main():
    part_a()
    part_b()
    rows = part_c()
    part_d()
    part_e()
    part_f()
    figure_uniqueness_map()
    figure_hbar_eff(rows)
    print("figures written to", output_path(""))


if __name__ == "__main__":
    main()
