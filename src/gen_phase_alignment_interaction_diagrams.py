"""
Interaction diagrams for the phase-alignment microdynamics.

Companion to ``docs/supplement/phase_alignment_interaction_diagrams.md``,
``docs/analysis/phase_alignment_microdynamics.md`` and
``docs/algorithm/phase_alignment_microdynamics_algorithm.md``.

The phase-space crystal-lattice supplement (§4) draws each candidate reading
of a Wigner momentum-jump event twice: once as a lattice picture in the
(x, p) plane and once as a set of trajectories in the (x, t) plane.  This
script does the same for the *phase-alignment* microdynamics, in which the
only relational state is the misalignment

    mu  =  Phi_out(x*, t*) - Phi_in(x*, t*),
    Phi_j(x, t)  =  theta_j + [p_j (x - x_j) - E_j t] / hbar,

of the two mediating sea rows at the vertex.  There are exactly four
elementary processes, plus the free leg:

  0  free leg          p constant; only the clock winds
  1  pump              V acts on phase alone; mu becomes place-valued
  2  exchange s = +1   excess particle hops up by 2q dp
  3  exchange s = -1   excess particle hops down by 2q dp
  4  suppressed        off-stationary candidates and the incoherent sea

Parts A-C below verify the claims annotated on the figures before drawing
them:

  A  the vertex is a swap: Sum p and Sum p^2 are preserved identically,
     the union of the two worldlines is unchanged, and the stationarity
     condition v_exch = vbar_pair is what selects it;

  B  the sinc dephasing envelope for off-stationary candidates;

  C  **the pump sign.**  Building the sea as an actual ring amplitude and
     applying the pump as a multiplicative phase, the pumped sidebands at
     +/- 2q dp carry a definite misalignment relative to the carrier row.
     Which sign of the kick reproduces Lemma 5 -- and hence Gamma_q(x)
     with the corrected sign of the supplement §6.3 -- is checked here
     rather than assumed.  See §7 of the companion document.

Figures (written via output_path and docs_path):

    pa_int_0_free_leg.png
    pa_int_1_pump.png
    pa_int_2_exchange_up.png
    pa_int_3_exchange_down.png
    pa_int_4_suppressed.png

Run with:

    WPMW_OUTPUT=/mnt/user-data/outputs python -u \
        src/gen_phase_alignment_interaction_diagrams.py
"""

from __future__ import annotations

import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Circle, FancyArrowPatch, Polygon, Wedge  # noqa: E402

sys.path.insert(0, os.path.dirname(__file__))

from wpmwlib.wpmw_utils import output_path, docs_path  # noqa: E402


# --------------------------------------------------------------------- #
# Physics                                                               #
# --------------------------------------------------------------------- #
HBAR = 1.0
MASS = 1.0
L = 8.0                       # ring circumference
Q = 1                         # pumped mode index
DP = np.pi * HBAR / L         # Wigner half-grid quantum
KQ = 2.0 * np.pi * Q / L      # mode wavenumber
QUANT = 2 * Q * DP            # momentum quantum = hbar * K_q = h/L
V_Q = 0.35                    # mode amplitude
PHI_Q = 0.7                   # mode phase offset (nonzero on purpose)
TAU_P = 0.05                  # pump interval
MU_1 = V_Q * TAU_P / HBAR     # pumped contrast


def gamma_q(x):
    """Stencil rate field, corrected sign of supplement §6.3."""
    return -(V_Q / HBAR) * np.sin(KQ * x + PHI_Q)


def potential(x):
    return V_Q * np.cos(KQ * x + PHI_Q)


# --------------------------------------------------------------------- #
# Part A -- the vertex is a swap                                        #
# --------------------------------------------------------------------- #
def part_a() -> None:
    print("Part A: the vertex is a momentum swap (Theorem 4)")
    r_in = 2                                  # excess particle's own row
    for s in (+1, -1):
        r_out = r_in + 2 * Q * s
        p_b = r_in * DP                       # in-row sea leg
        p_a = r_out * DP                      # out-row sea leg
        # Solve the two conditions of Theorem 4 for (p_in, p_out):
        #   p_out - p_in = p_a - p_b        (momentum conservation)
        #   p_out + p_in = p_a + p_b        (mu stationary)
        m_sys = np.array([[-1.0, 1.0], [1.0, 1.0]])
        rhs = np.array([p_a - p_b, p_a + p_b])
        p_in, p_out = np.linalg.solve(m_sys, rhs)
        res = max(abs(p_in - p_b), abs(p_out - p_a))
        d_p = (p_out + p_a_out(p_a, p_b)) - (p_in + p_a)
        d_e = ((p_out ** 2 + p_a_out(p_a, p_b) ** 2)
               - (p_in ** 2 + p_a ** 2)) / (2 * MASS)
        v_exch = (p_in + p_out) / (2 * MASS)
        v_pair = (p_a + p_b) / (2 * MASS)
        print(f"  s = {s:+d}:  rows {r_in} -> {r_out}, "
              f"legs at p_b = {p_b:+.6f}, p_a = {p_a:+.6f}")
        print(f"           solution  p_in = {p_in:+.6f} (= p_b), "
              f"p_out = {p_out:+.6f} (= p_a),  residual = {res:.3e}")
        print(f"           transfer  |p_out - p_in| = {abs(p_out - p_in):.6f}"
              f"   2q dp = {QUANT:.6f}   diff = "
              f"{abs(abs(p_out - p_in) - QUANT):.3e}")
        print(f"           Sum p conserved to {abs(d_p):.3e}, "
              f"Sum p^2/2m conserved to {abs(d_e):.3e}  (not imposed)")
        print(f"           v_exch = {v_exch:+.6f}, "
              f"vbar_pair = {v_pair:+.6f}, diff = {abs(v_exch-v_pair):.3e}")
    # The union of the two worldlines is unchanged by the swap.
    p_b, p_a = r_in * DP, (r_in + 2 * Q) * DP
    t = np.linspace(-1.0, 1.0, 401)
    xstar = 3.2
    # labelled paths, after the swap
    x_exc = xstar + np.where(t < 0, p_b, p_a) * t / MASS
    x_sea = xstar + np.where(t < 0, p_a, p_b) * t / MASS
    # unlabelled straight lines through the same event
    x_lo = xstar + p_b * t / MASS
    x_hi = xstar + p_a * t / MASS
    set_after = np.sort(np.stack([x_exc, x_sea]), axis=0)
    set_free = np.sort(np.stack([x_lo, x_hi]), axis=0)
    print(f"  union of the two paths, swapped vs. free-crossing: "
          f"max deviation = {np.max(np.abs(set_after - set_free)):.3e}")
    print()


def p_a_out(p_a, p_b):
    """Struck sea leg leaves on its in-row partner's momentum."""
    return p_b


# --------------------------------------------------------------------- #
# Part B -- the dephasing envelope                                      #
# --------------------------------------------------------------------- #
def part_b() -> None:
    print("Part B: off-stationary candidates dephase (sinc envelope)")
    tau_e = 60.0
    mu0 = 1.0
    d_p = QUANT
    v_pair = 0.5
    print(f"  tau_e = {tau_e}, mu_0 = {mu0}, |Delta p| = {d_p:.4f}")
    print(f"  {'v_exch - vbar':>14s}  {'<cos mu>_t':>12s}  "
          f"{'analytic':>12s}  {'|sinc| bound':>12s}")
    for dv in (0.0, 0.02, 0.05, 0.1, 0.25):
        mudot = (d_p / HBAR) * dv
        t = np.linspace(0.0, tau_e, 200001)
        avg = np.mean(np.cos(mu0 + mudot * t))
        sinc = np.sinc(mudot * tau_e / (2 * np.pi))
        exact = sinc * np.cos(mu0 + mudot * tau_e / 2.0)
        print(f"  {dv:14.3f}  {avg:12.6f}  {exact:12.6f}  "
              f"{abs(sinc):12.6f}")
    print()


# --------------------------------------------------------------------- #
# Part C -- the pump sign                                               #
# --------------------------------------------------------------------- #
def part_c() -> None:
    """Pump a phase-locked sea row and read mu^s off the sidebands.

    The sea on a row is a carrier plane wave on the ring.  The pump is a
    multiplicative phase exp(-i sigma V(x) tau_p / hbar) with sigma = +1
    for the Lagrangian-winding convention of §3 and §6 of the algorithm
    specification, sigma = -1 for the opposite.  First order in mu_1
    gives sidebands at +/- hbar K_q; the misalignment between the
    sideband row r + 2sq and the carrier row r is read directly.
    """
    print("Part C: the pump as a sideband generator -- which kick sign "
          "gives Gamma_q?")
    n_grid = 2048
    x = np.arange(n_grid) * L / n_grid
    r = 2
    k_r = r * DP / HBAR
    for sigma, tag in ((+1.0, "theta -> theta - V tau_p / hbar  (spec §6)"),
                       (-1.0, "theta -> theta + V tau_p / hbar  (opposite)")):
        psi = np.exp(1j * k_r * x) * np.exp(-1j * sigma * potential(x)
                                            * TAU_P / HBAR)
        amp = np.fft.fft(psi) / n_grid
        kk = np.fft.fftfreq(n_grid, d=L / n_grid) * 2.0 * np.pi
        print(f"  {tag}")
        rate = np.zeros_like(x)
        for s in (+1, -1):
            k_out = k_r + s * KQ
            j_out = int(np.argmin(np.abs(kk - k_out)))
            j_in = int(np.argmin(np.abs(kk - k_r)))
            # transported phases of the two rows, as fields on the ring
            phi_out = np.angle(amp[j_out]) + k_out * x
            phi_in = np.angle(amp[j_in]) + k_r * x
            mu_s = phi_out - phi_in
            pred = s * (KQ * x + PHI_Q) + np.pi / 2.0     # Lemma 5, t = 0
            dev = np.max(np.abs(np.angle(np.exp(1j * (mu_s - pred)))))
            frac = abs(amp[j_out]) / abs(amp[j_in])
            print(f"    s = {s:+d}: |A_out/A_in| = {frac:.6e} "
                  f"(mu_1/2 = {MU_1/2:.6e}),  "
                  f"max|mu^s - Lemma 5| = {dev:.3e}")
            rate += 0.5 * MU_1 * s * np.cos(mu_s) / TAU_P
        num = np.max(np.abs(rate - gamma_q(x))) / np.max(np.abs(gamma_q(x)))
        num_flip = (np.max(np.abs(rate + gamma_q(x)))
                    / np.max(np.abs(gamma_q(x))))
        print(f"    assembled rate vs. Gamma_q : rel. dev = {num:.3e};  "
              f"vs. -Gamma_q : rel. dev = {num_flip:.3e}")
    print()


# --------------------------------------------------------------------- #
# Drawing primitives                                                    #
# --------------------------------------------------------------------- #
C_EXC = "#B45309"        # excess particle (positon above the sea)
C_SEA = "#0F766E"        # sea leg
C_MU = "#6D28D9"         # misalignment
C_LINE = "#1F2937"
C_NET = "#B91C1C"
C_TEXT = "#111827"
C_SOFT = "#6B7280"
C_FAINT = "#E5E7EB"
C_ROW = "#D1D5DB"


def _panel(ax, title):
    ax.set_title(title, loc="left", fontsize=11.0, fontweight="bold",
                 color=C_TEXT, pad=8)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.spines["left"].set_color(C_SOFT)
    ax.spines["bottom"].set_color(C_SOFT)
    ax.tick_params(colors=C_SOFT, labelsize=8.5)


def _dial(ax, cx, cy, rad, ang_out, ang_in, wedge=True, zorder=9):
    ax.add_patch(Circle((cx, cy), rad, facecolor="white", edgecolor=C_ROW,
                        lw=1.0, zorder=zorder))
    if wedge and abs(np.angle(np.exp(1j * (ang_out - ang_in)))) > 0.08:
        ax.add_patch(Wedge((cx, cy), rad * 0.88, np.degrees(ang_in),
                           np.degrees(ang_out), facecolor=C_MU, alpha=0.22,
                           edgecolor="none", zorder=zorder + 1))
    ax.plot([cx, cx + rad * 0.86 * np.cos(ang_in)],
            [cy, cy + rad * 0.86 * np.sin(ang_in)],
            color=C_SEA, lw=1.6, solid_capstyle="round", zorder=zorder + 2)
    ax.plot([cx, cx + rad * 0.86 * np.cos(ang_out)],
            [cy, cy + rad * 0.86 * np.sin(ang_out)],
            color=C_EXC, lw=2.4, solid_capstyle="round", zorder=zorder + 3)
    ax.add_patch(Circle((cx, cy), rad * 0.09, facecolor=C_LINE,
                        edgecolor="none", zorder=zorder + 4))


def _rows(ax, rows_even=(-2, 0, 2, 4, 6), rows_odd=(-1, 1, 3, 5),
          x0=0.0, x1=L, label=True):
    """Draw the momentum rows of the phase-space crystal lattice."""
    for n in rows_odd:
        ax.plot([x0, x1], [n * DP, n * DP], color=C_FAINT, lw=1.0,
                ls=(0, (3, 3)), zorder=1)
    for n in rows_even:
        ax.plot([x0, x1], [n * DP, n * DP], color=C_ROW, lw=1.6, zorder=1)
        if label:
            ax.text(x1 + 0.10, n * DP, f"$r={n}$", fontsize=8.0,
                    color=C_SOFT, va="center", ha="left")


def _cells(ax, ymin, ymax, m=8):
    for i in range(m + 1):
        ax.plot([i * L / m, i * L / m], [ymin, ymax], color="#F3F4F6",
                lw=0.8, zorder=0)


def _sea_dot(ax, x, p, zorder=6, ms=5.5):
    ax.plot([x], [p], "o", mfc="white", mec=C_SEA, mew=1.6, ms=ms,
            zorder=zorder)


def _arrow(ax, xy0, xy1, color, lw=2.0, scale=13, ls="-", zorder=8):
    ax.add_patch(FancyArrowPatch(xy0, xy1, arrowstyle="-|>",
                                 mutation_scale=scale, color=color, lw=lw,
                                 linestyle=ls, zorder=zorder,
                                 shrinkA=0, shrinkB=0))


def _save(fig, name):
    fig.savefig(output_path(name), dpi=150, bbox_inches="tight")
    dp_ = docs_path(name)
    if dp_:
        fig.savefig(dp_, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {name}")


# --------------------------------------------------------------------- #
# Figure 0 -- the free leg                                              #
# --------------------------------------------------------------------- #
def fig_free_leg():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.6, 4.6))

    # (a) space-momentum
    _panel(ax1, "(a)  space–momentum:  the leg never leaves its row")
    _cells(ax1, -2.6 * DP, 6.6 * DP)
    _rows(ax1)
    p = 2 * DP
    _arrow(ax1, (1.0, p), (6.4, p), C_EXC, lw=2.6)
    ax1.plot([1.0], [p], "o", color=C_EXC, ms=8, zorder=9)
    for xd in (1.0, 3.7, 6.4):
        ang = (p ** 2 / (2 * MASS) - potential(xd)) / HBAR * 2.2
        _dial(ax1, xd, p + 1.35 * DP, 0.34, ang, ang, wedge=False)
    ax1.text(3.7, p + 0.38 * DP, r"$\dot x = p/m,\qquad \dot p = 0$",
             fontsize=10, ha="center", color=C_EXC)
    ax1.text(3.7, p + 2.35 * DP,
             r"clock winds at $\dot\theta = [\,p^2/2m - V(x)\,]/\hbar$",
             fontsize=9.2, ha="center", color=C_SOFT)
    xs = np.linspace(0, L, 400)
    ax1.plot(xs, -1.9 * DP + 0.55 * DP * potential(xs) / V_Q,
             color="#9CA3AF", lw=1.5, zorder=2)
    ax1.text(0.15, -2.45 * DP, r"$V(x)$", fontsize=9, color="#9CA3AF")
    ax1.set_xlim(-0.2, L + 0.9)
    ax1.set_ylim(-2.9 * DP, 6.9 * DP)
    ax1.set_xlabel("position $x$", fontsize=9.5, color=C_SOFT)
    ax1.set_ylabel("momentum $p$", fontsize=9.5, color=C_SOFT)

    # (b) space-time
    _panel(ax2, "(b)  space–time:  a straight worldline, whatever $V$ does")
    t_max = 5.0
    v = p / MASS
    x0 = 1.0
    ax2.plot([x0, x0 + v * t_max], [0, t_max], color=C_EXC, lw=2.8, zorder=6)
    _arrow(ax2, (x0 + v * t_max * 0.93, t_max * 0.93),
           (x0 + v * t_max, t_max), C_EXC, lw=2.8)
    for tt in (0.6, 2.3, 4.0):
        ang = (p ** 2 / (2 * MASS)
               - potential(x0 + v * tt)) / HBAR * tt * 0.9
        _dial(ax2, x0 + v * tt + 0.72, tt, 0.30, ang, ang, wedge=False)
    ax2.text(x0 + v * 2.6 - 1.55, 2.6,
             r"slope $\dfrac{dt}{dx} = \dfrac{m}{p}$", fontsize=10,
             ha="center", color=C_EXC)
    ax2.text(L * 0.52, 0.30,
             "no bend: the potential is not a force here",
             fontsize=9.2, ha="center", color=C_SOFT, style="italic")
    ax2.set_xlim(-0.2, L + 0.9)
    ax2.set_ylim(-0.25, t_max + 0.45)
    ax2.set_xlabel("position $x$", fontsize=9.5, color=C_SOFT)
    ax2.set_ylabel("time $t$", fontsize=9.5, color=C_SOFT)

    fig.suptitle("Process 0 — the free leg:  momentum is constant, "
                 "only the carried clock advances",
                 fontsize=12.5, fontweight="bold", color=C_TEXT, y=1.015)
    fig.tight_layout()
    _save(fig, "pa_int_0_free_leg.png")


# --------------------------------------------------------------------- #
# Figure 1 -- the pump                                                  #
# --------------------------------------------------------------------- #
def fig_pump():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.6, 4.9))
    r_in, r_out = 2, 4
    p_in, p_out = r_in * DP, r_out * DP

    _panel(ax1, "(a)  space–momentum:  the pump writes a grating, "
                "not a force")
    _cells(ax1, -2.9 * DP, 7.4 * DP)
    _rows(ax1, rows_even=(0, 2, 4, 6), rows_odd=(1, 3, 5))
    for xx in np.linspace(0.45, L - 0.45, 9):
        _sea_dot(ax1, xx, p_in)
        _sea_dot(ax1, xx, p_out)
    for xd in np.linspace(0.8, L - 0.8, 5):
        mu = KQ * xd + PHI_Q + np.pi / 2.0
        _dial(ax1, xd, p_out + 1.55 * DP, 0.36, mu, 0.0)
    ax1.text(L / 2, p_out + 3.05 * DP,
             r"$\mu^{+}(x) = K_q x + \phi_q + \pi/2$   "
             r"(Lemma 5: a function of place alone)",
             fontsize=9.6, ha="center", color=C_MU)
    ax1.text(L / 2, p_in - 0.75 * DP,
             r"populations unchanged at $O(V_q)$  (Lemma 3)",
             fontsize=9.2, ha="center", color=C_SOFT, style="italic")
    xs = np.linspace(0, L, 500)
    ax1.plot(xs, -2.1 * DP + 0.55 * DP * potential(xs) / V_Q,
             color="#9CA3AF", lw=1.5, zorder=2)
    ax1.plot(xs, -2.1 * DP - 0.55 * DP * np.sin(KQ * xs + PHI_Q),
             color=C_NET, lw=2.2, zorder=3)
    ax1.text(0.12, -2.85 * DP, r"$V(x)$", fontsize=9, color="#9CA3AF")
    ax1.text(2.05, -2.85 * DP, r"$\cos\mu^{+}(x)$ — in quadrature",
             fontsize=9, color=C_NET)
    ax1.set_xlim(-0.2, L + 0.9)
    ax1.set_ylim(-3.3 * DP, 7.9 * DP)
    ax1.set_xlabel("position $x$", fontsize=9.5, color=C_SOFT)
    ax1.set_ylabel("momentum $p$", fontsize=9.5, color=C_SOFT)

    _panel(ax2, "(b)  space–time:  a phase gradient, not a force")
    t_max = 5.0
    t_p = 2.5
    ax2.axhspan(t_p - 0.10, t_p + 0.10, color=C_MU, alpha=0.18, zorder=1)
    ax2.text(L + 0.20, t_p, r"pump $\tau_p$", fontsize=9, color=C_MU,
             va="center")
    v = p_in / MASS
    # three sea legs on the same row, at three places
    for x0 in (0.8, 2.6, 4.4):
        ax2.plot([x0, x0 + v * t_max], [0, t_max], color=C_SEA, lw=2.2,
                 zorder=5)
        x_p = x0 + v * t_p
        kick = -potential(x_p) * TAU_P / HBAR
        # exaggerate the kick for legibility; the ratio between the three
        # is faithful, the common scale is not
        scale = 1.0 / MU_1
        _dial(ax2, x_p - 0.80, t_p - 1.05, 0.30, 0.0, 0.0, wedge=False)
        _dial(ax2, x_p + v * 1.05 - 0.80, t_p + 1.05, 0.30,
              kick * scale, 0.0)
        ax2.text(x_p + 0.05, t_p + 0.30,
                 rf"$V={potential(x_p):+.2f}$", fontsize=8.2, color=C_SOFT,
                 ha="left")
    ax2.text(0.05, t_max + 0.10,
             r"$\theta_j \mapsto \theta_j \mp V(x_j)\,\tau_p/\hbar$ — the "
             r"kick differs from place to place because $V$ does;"
             "\nno worldline bends, and the row populations do not change",
             fontsize=9.2, ha="left", va="top", color=C_SOFT, style="italic")
    ax2.set_xlim(-0.2, L + 1.5)
    ax2.set_ylim(-0.25, t_max + 1.15)
    ax2.set_xlabel("position $x$", fontsize=9.5, color=C_SOFT)
    ax2.set_ylabel("time $t$", fontsize=9.5, color=C_SOFT)

    fig.suptitle("Process 1 — the pump:  the only place the potential "
                 "touches a world-particle, and it touches phase alone",
                 fontsize=12.5, fontweight="bold", color=C_TEXT, y=1.015)
    fig.tight_layout()
    _save(fig, "pa_int_1_pump.png")


# --------------------------------------------------------------------- #
# Figures 2 and 3 -- the exchange, both directions                      #
# --------------------------------------------------------------------- #
def fig_exchange(s: int, name: str):
    r_in = 2
    r_out = r_in + 2 * Q * s
    p_in, p_out = r_in * DP, r_out * DP
    n_mid = r_in + Q * s
    up = s > 0

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.6, 5.0))
    xstar = 3.2

    # ---------------- (a) space-momentum ---------------- #
    _panel(ax1, f"(a)  space–momentum:  one vertex, $s = {s:+d}$")
    rows_even = (0, 2, 4)
    rows_odd = (1, 3)
    _cells(ax1, -1.6 * DP, 5.9 * DP)
    _rows(ax1, rows_even=rows_even, rows_odd=rows_odd)
    ax1.text(L + 0.10, n_mid * DP, f"$n={n_mid}$ (empty)", fontsize=8.0,
             color="#9CA3AF", va="center", ha="left")
    for xx in np.linspace(0.5, L - 0.5, 9):
        _sea_dot(ax1, xx, p_in)
        _sea_dot(ax1, xx, p_out)

    # incoming / outgoing drift of the excess particle
    _arrow(ax1, (xstar - 2.3, p_in), (xstar - 0.16, p_in), C_EXC, lw=2.4)
    _arrow(ax1, (xstar + 0.16, p_out), (xstar + 2.3, p_out), C_EXC, lw=2.4)
    # the vertex itself
    _arrow(ax1, (xstar, p_in), (xstar, p_out), C_EXC, lw=2.8, scale=15)
    _arrow(ax1, (xstar + 0.34, p_out), (xstar + 0.34, p_in), C_SEA, lw=2.2,
           scale=13)
    ax1.plot([xstar], [p_in], "o", color=C_EXC, ms=8, zorder=10)
    ax1.plot([xstar + 0.34], [p_out], "o", mfc="white", mec=C_SEA, mew=2.0,
             ms=8, zorder=10)

    y_lab = p_out + (0.55 if up else -0.85) * DP
    ax1.text(xstar + 0.17, y_lab,
             rf"excess: $r \to r{2*Q*s:+d}$,   sea leg: $r{2*Q*s:+d} \to r$",
             fontsize=9.6, ha="center", color=C_TEXT)
    ax1.text(xstar + 0.17, y_lab + (0.55 if up else -0.55) * DP,
             rf"$|\Delta p| = 2q\,dp = {QUANT:.4f}$",
             fontsize=9.2, ha="center", color=C_SOFT)

    y_dial = 0.5 * (r_in + r_out) * DP
    _dial(ax1, 6.9, y_dial, 0.40,
          s * (KQ * xstar + PHI_Q) + np.pi / 2, 0.0)
    ax1.text(6.9, y_dial + 0.68,
             rf"$\mu^{{{'+' if up else '-'}}}(x^{{*}})$ between the two rows",
             fontsize=9.4, ha="center", color=C_MU)
    ax1.text(6.9, y_dial - 0.76,
             r"$P = w_0 + \kappa\cos\mu$", fontsize=9.4, ha="center",
             color=C_MU)
    ax1.set_xlim(-0.2, L + 1.5)
    ax1.set_ylim(-2.0 * DP, 6.3 * DP)
    ax1.set_xlabel("position $x$", fontsize=9.5, color=C_SOFT)
    ax1.set_ylabel("momentum $p$", fontsize=9.5, color=C_SOFT)

    # ---------------- (b) space-time ---------------- #
    _panel(ax2, "(b)  space–time:  the swap permutes labels, "
                "not trajectories")
    t_max, tstar = 4.0, 1.8
    v_in, v_out = p_in / MASS, p_out / MASS

    # shaded encounter window along the mean-velocity axis
    v_bar = 0.5 * (v_in + v_out)
    dt_e = 0.42
    lens = [(xstar - v_bar * dt_e - 0.28, tstar - dt_e),
            (xstar - v_bar * dt_e + 0.28, tstar - dt_e),
            (xstar + v_bar * dt_e + 0.28, tstar + dt_e),
            (xstar + v_bar * dt_e - 0.28, tstar + dt_e)]
    ax2.add_patch(Polygon(lens, closed=True, facecolor=C_MU, alpha=0.13,
                          edgecolor="none", zorder=2))
    ax2.plot([xstar - v_bar * 1.3, xstar + v_bar * 1.3],
             [tstar - 1.3, tstar + 1.3], color=C_MU, lw=1.4,
             ls=(0, (4, 3)), zorder=3)
    ax2.text(xstar + v_bar * 1.35 + 0.10, tstar + 1.35,
             r"$v_{\rm exch} = \bar v_{\rm pair}$" "\n"
             r"$\Rightarrow \dot\mu = 0$",
             fontsize=9.0, color=C_MU, ha="left", va="center")

    # excess particle: in on row r, out on row r + 2sq
    ax2.plot([xstar - v_in * tstar, xstar], [0, tstar], color=C_EXC, lw=2.8,
             zorder=6)
    ax2.plot([xstar, xstar + v_out * (t_max - tstar)], [tstar, t_max],
             color=C_EXC, lw=2.8, zorder=6)
    # sea leg: in on row r + 2sq, out on row r
    ax2.plot([xstar - v_out * tstar, xstar], [0, tstar], color=C_SEA, lw=2.0,
             zorder=5)
    ax2.plot([xstar, xstar + v_in * (t_max - tstar)], [tstar, t_max],
             color=C_SEA, lw=2.0, zorder=5)
    ax2.plot([xstar], [tstar], "o", color=C_MU, ms=8, zorder=11)

    ax2.text(xstar - v_in * tstar + 0.12, 0.10, "excess in", fontsize=9.0,
             color=C_EXC, ha="left")
    ax2.text(xstar - v_out * tstar + 0.12, 0.52, "sea leg in", fontsize=9.0,
             color=C_SEA, ha="left")
    ax2.text(xstar + v_out * (t_max - tstar) + 0.10, t_max - 0.05,
             "excess out", fontsize=9.0, color=C_EXC, ha="left", va="top")
    ax2.text(xstar + v_in * (t_max - tstar) + 0.10, t_max - 0.62,
             "sea leg out", fontsize=9.0, color=C_SEA, ha="left", va="top")
    ax2.text(0.05, -0.55,
             "the union of the two paths is exactly two straight lines "
             "crossing: only the labels are exchanged,\nso "
             r"$\sum p$ and $\sum p^{2}$ are conserved identically "
             r"(Corollary 4.1) — the vertex moves no trajectory",
             fontsize=9.2, ha="left", va="top", color=C_SOFT,
             style="italic")
    ax2.set_xlim(-0.5, L + 1.8)
    ax2.set_ylim(-1.45, t_max + 0.35)
    ax2.set_xlabel("position $x$", fontsize=9.5, color=C_SOFT)
    ax2.set_ylabel("time $t$", fontsize=9.5, color=C_SOFT)

    word = "up" if up else "down"
    fig.suptitle(f"Process {2 if up else 3} — the exchange, $s = {s:+d}$:  "
                 f"the excess particle hops {word} by $2q\\,dp$ and one sea "
                 f"leg hops the other way",
                 fontsize=12.5, fontweight="bold", color=C_TEXT, y=1.035)
    fig.tight_layout()
    _save(fig, name)


# --------------------------------------------------------------------- #
# Figure 4 -- the suppressed channels                                   #
# --------------------------------------------------------------------- #
def fig_suppressed():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.6, 5.0))
    xstar = 3.2
    r_in = 2
    p_in = r_in * DP

    _panel(ax1, "(a)  space–momentum:  candidates that are never forbidden, "
                "only dark")
    _cells(ax1, -3.2 * DP, 7.2 * DP)
    _rows(ax1, rows_even=(-2, 0, 2, 4, 6), rows_odd=(-1, 1, 3, 5))
    for xx in np.linspace(0.5, L - 0.5, 9):
        for rr in (0, 2, 4):
            _sea_dot(ax1, xx, rr * DP, ms=4.6)

    cands = [(4, C_NET, "stationary"),
             (6, "#9CA3AF", r"off-stationary"),
             (-2, "#9CA3AF", r"off-stationary")]
    for k, (r_o, col, lab) in enumerate(cands):
        xk = xstar + 0.62 * k
        _arrow(ax1, (xk, p_in), (xk, r_o * DP), col,
               lw=2.6 if k == 0 else 1.8,
               ls="-" if k == 0 else (0, (4, 2)))
        va = "bottom" if r_o > r_in else "top"
        off = 0.30 * DP if r_o > r_in else -0.30 * DP
        ax1.text(xk + 0.10, r_o * DP + off, lab, fontsize=8.6, color=col,
                 ha="left", va=va)
    ax1.plot([xstar], [p_in], "o", color=C_EXC, ms=8, zorder=10)
    ax1.text(0.05, 8.3 * DP,
             "only the stationary candidate survives averaging (Theorem 4);\n"
             "the others are suppressed by "
             r"$|\mathrm{sinc}(\dot\mu\,\tau_e/2)|$, not forbidden",
             fontsize=9.2, ha="left", va="top", color=C_SOFT, style="italic")

    # incoherent sea inset: random dials
    rng = np.random.default_rng(3)
    for j, xd in enumerate(np.linspace(0.75, 2.6, 4)):
        _dial(ax1, xd, -2.2 * DP, 0.26, rng.uniform(0, 2 * np.pi),
              rng.uniform(0, 2 * np.pi))
    ax1.text(2.95, -2.2 * DP,
             "incoherent sea: " r"$|\hat Z_r| = O(N^{-1/2})$," "\n"
             r"bias $O(N^{-1})$ — dark with no extra postulate",
             fontsize=9.0, ha="left", va="center", color=C_SOFT)
    ax1.set_xlim(-0.2, L + 1.2)
    ax1.set_ylim(-3.6 * DP, 9.0 * DP)
    ax1.set_xlabel("position $x$", fontsize=9.5, color=C_SOFT)
    ax1.set_ylabel("momentum $p$", fontsize=9.5, color=C_SOFT)

    _panel(ax2, "(b)  space–time:  the event axis must ride the pair's "
                "mean velocity")
    t_max, tstar = 4.0, 1.8
    v_pair = (p_in + 4 * DP) / (2 * MASS)
    for dv, col, lw, ls, lab in ((0.0, C_NET, 2.6, "-", r"$\dot\mu = 0$"),
                                 (0.55, "#9CA3AF", 1.8, (0, (4, 2)),
                                  r"$\dot\mu \neq 0$")):
        v = v_pair + dv
        ax2.plot([xstar - v * tstar, xstar + v * (t_max - tstar)],
                 [0, t_max], color=col, lw=lw, ls=ls, zorder=5)
        ax2.text(xstar + v * (t_max - tstar) + 0.08, t_max, lab,
                 fontsize=9.2, color=col, ha="left", va="top")
    ax2.plot([xstar], [tstar], "o", color=C_MU, ms=8, zorder=10)
    for tt in (tstar - 1.2, tstar + 1.2):
        _dial(ax2, xstar + v_pair * (tt - tstar) - 0.62, tt, 0.26, 1.0, 0.0)
    ax2.text(xstar - v_pair * 1.2 - 0.62, tstar - 1.72,
             r"$\mu$ frozen", fontsize=8.8, color=C_NET, ha="center")
    for tt in (tstar - 1.2, tstar + 1.2):
        v = v_pair + 0.55
        mu = 1.0 + (QUANT / HBAR) * 0.55 * (tt - tstar)
        _dial(ax2, xstar + v * (tt - tstar) + 0.62, tt, 0.26, mu, 0.0)
    ax2.text(xstar - (v_pair + 0.55) * 1.2 + 0.62, tstar - 1.72,
             r"$\mu$ slides", fontsize=8.8, color=C_SOFT, ha="center")

    ax_in = ax2.inset_axes([0.06, 0.06, 0.36, 0.28])
    tau_e = 60.0
    dv = np.linspace(-0.3, 0.3, 601)
    env = np.abs(np.sinc((QUANT / HBAR) * dv * tau_e / (2 * np.pi)))
    ax_in.plot(dv, env, color=C_NET, lw=1.6)
    ax_in.set_xticks([-0.2, 0, 0.2])
    ax_in.set_yticks([0, 1])
    ax_in.tick_params(labelsize=7, colors=C_SOFT)
    ax_in.set_xlabel(r"$v_{\rm exch} - \bar v_{\rm pair}$", fontsize=7.5,
                     color=C_SOFT, labelpad=1)
    ax_in.set_title(r"$|\mathrm{sinc}(\dot\mu\,\tau_e/2)|$", fontsize=8,
                    color=C_SOFT, pad=3)
    for sp in ax_in.spines.values():
        sp.set_color(C_FAINT)
    ax2.set_xlim(-0.6, L + 0.6)
    ax2.set_ylim(-0.25, t_max + 0.35)
    ax2.set_xlabel("position $x$", fontsize=9.5, color=C_SOFT)
    ax2.set_ylabel("time $t$", fontsize=9.5, color=C_SOFT)

    fig.suptitle("Process 4 — the suppressed channels:  nothing is "
                 "forbidden; misaligned traffic simply averages away",
                 fontsize=12.5, fontweight="bold", color=C_TEXT, y=1.015)
    fig.tight_layout()
    _save(fig, "pa_int_4_suppressed.png")


# --------------------------------------------------------------------- #
# Main                                                                  #
# --------------------------------------------------------------------- #
def main():
    print("=" * 72)
    print("Phase-alignment microdynamics: interaction diagrams")
    print("=" * 72)
    print(f"  hbar = {HBAR}, m = {MASS}, L = {L}, q = {Q}")
    print(f"  dp = pi hbar / L = {DP:.6f}, 2q dp = hbar K_q = {QUANT:.6f}")
    print(f"  V_q = {V_Q}, phi_q = {PHI_Q}, tau_p = {TAU_P}, "
          f"mu_1 = {MU_1:.6e}")
    print()
    part_a()
    part_b()
    part_c()
    print("Figures:")
    fig_free_leg()
    fig_pump()
    fig_exchange(+1, "pa_int_2_exchange_up.png")
    fig_exchange(-1, "pa_int_3_exchange_down.png")
    fig_suppressed()
    print()
    print("done.")


if __name__ == "__main__":
    main()
