"""
Demonstration: the contact interaction in phase-alignment variables.

Companion to ``docs/analysis/phase_alignment_microdynamics.md`` and
``docs/algorithm/phase_alignment_microdynamics_algorithm.md``.

The phase-resonance note describes the sea's excitations as *beats* and the
vertex selection rule as a *resonance* between a beat and a transition.  This
demo verifies that both can be eliminated in favour of a single scalar carried
by a pair -- the **misalignment**

    mu  =  Phi_a(x*, t*) - Phi_b(x*, t*),
    Phi_j(x, t)  =  theta_j + [p_j (x - x_j) - E_j t] / hbar,

the difference of the two partners' transported clock phases, evaluated at one
point.  Six checks:

  A  mu is the complete gauge-invariant relational datum of a pair: it is
     invariant under a global phase shift and under a change of the phase
     reference point, and the pair amplitude depends on the partners only
     through mu.

  B  mu winds in space at d(mu)/dx = dp/hbar and along a path of velocity v
     at d(mu)/dt = (dp/hbar)(v - vbar_pair).  "Beating" is nothing but
     dp != 0; no propagating object is involved.

  C  **Exchange theorem.**  Requiring mu to be stationary along the
     transition path, together with momentum conservation, forces the vertex
     to be a *momentum swap* between the excess particle and the struck
     partner.  Energy conservation is then automatic -- it is not imposed.
     This replaces Proposition 2 plus the row-selection rule of the
     phase-resonance note.

  D  Locality of mu under the pump (Lemma 2 of the phase-resonance note,
     restated): every pumped pair at the same place carries the same mu,
     independent of where the pair sits.  The arriving particle never needs
     the absent partner's state.

  E  The exchange bias cos(mu(x)) is in quadrature with V(x), i.e.
     proportional to -V'(x); and the contact-vertex cross term is linear in
     the contrast C with no fitted phase offset.

  F  The resulting stencil preserves L1 on the ring.

Figure (written via docs_path for the output branch):
  phase_alignment_contact.png -- six panels: the one number; the three states
  of a pair; the pump; the vertex; the stationarity condition; the ensemble
  limit.

Run with:

    WPMW_OUTPUT=/mnt/user-data/outputs python -u demo_phase_alignment.py
"""

from __future__ import annotations

import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Circle, FancyArrowPatch, Rectangle, Wedge  # noqa: E402

sys.path.insert(0, os.path.dirname(__file__))

from wpmwlib.wpmw_utils import output_path, docs_path  # noqa: E402


# --------------------------------------------------------------------- #
# Physics                                                               #
# --------------------------------------------------------------------- #
HBAR = 1.0
MASS = 1.0

L = 2.0 * np.pi              # ring circumference
DP = np.pi * HBAR / L        # Wigner half-grid quantum = 0.5
Q = 1                        # pumped mode index
K = 2.0 * np.pi * Q / L      # pump wavevector = 1
QUANTUM = 2 * Q * DP         # momentum a q-split pair carries = 1
V_P = 1.0                    # pump amplitude

# a representative split pair
P_A, P_B = 2.5, 1.5          # struck partner, mate
X_A, X_B = 0.3, 1.7
TH_A, TH_B = 0.4, 2.1


def energy(p: float) -> float:
    return p * p / (2.0 * MASS)


def transported_phase(p, x0, th, x, t):
    """Phi_j(x, t): the phase particle j predicts at (x, t)."""
    return th + (p * (x - x0) - energy(p) * t) / HBAR


def misalignment(x, t, p_a=P_A, p_b=P_B, x_a=X_A, x_b=X_B,
                 th_a=TH_A, th_b=TH_B):
    """mu = Phi_a - Phi_b evaluated at (x, t)."""
    return (transported_phase(p_a, x_a, th_a, x, t)
            - transported_phase(p_b, x_b, th_b, x, t))


# --------------------------------------------------------------------- #
# Part A -- mu is the complete gauge-invariant relational datum          #
# --------------------------------------------------------------------- #
def part_a() -> None:
    print("Part A: mu is the pair's complete relational invariant")

    x_s, t_s = 0.83, 0.41
    mu0 = misalignment(x_s, t_s)

    # global phase shift: both clocks advanced together
    shift = 1.234
    mu1 = misalignment(x_s, t_s, th_a=TH_A + shift, th_b=TH_B + shift)
    print(f"  global phase shift          |d mu| = {abs(mu1 - mu0):.3e}")

    # re-reference each partner along its own worldline (transport is exact)
    dt = 0.27
    xa2, xb2 = X_A + P_A / MASS * dt, X_B + P_B / MASS * dt
    tha2 = TH_A + (energy(P_A) - 0.0) * dt / HBAR   # free winding, V = 0
    thb2 = TH_B + (energy(P_B) - 0.0) * dt / HBAR
    mu2 = misalignment(x_s, t_s, x_a=xa2, x_b=xb2, th_a=tha2, th_b=thb2)
    # re-referencing both worldlines forward by dt shifts the time origin
    # by dt: Phi_new(x, t) = Phi_old(x, t + dt).  Transport is exact.
    mu2_ref = misalignment(x_s, t_s + dt)
    print(f"  worldline re-reference      |d mu| = {abs(mu2 - mu2_ref):.3e}")

    # the pair amplitude depends on the partners only through mu
    psi = np.exp(1j * transported_phase(P_A, X_A, TH_A, x_s, t_s)) \
        - np.exp(1j * transported_phase(P_B, X_B, TH_B, x_s, t_s))
    pred = 2.0 * abs(np.sin(mu0 / 2.0))
    print(f"  |Psi| vs 2|sin(mu/2)|       |diff| = {abs(abs(psi) - pred):.3e}")
    print()


# --------------------------------------------------------------------- #
# Part B -- winding rates                                               #
# --------------------------------------------------------------------- #
def part_b() -> None:
    print("Part B: how mu winds")

    h = 1e-6
    dmu_dx = (misalignment(1.0 + h, 0.0) - misalignment(1.0 - h, 0.0)) / (2 * h)
    print(f"  d(mu)/dx = {dmu_dx:.12f}   expect dp/hbar = "
          f"{(P_A - P_B) / HBAR:.12f}")

    vbar = 0.5 * (P_A + P_B) / MASS
    worst = 0.0
    for v in (0.0, 1.4, vbar, 3.1):
        d = (misalignment(v * h, h) - misalignment(-v * h, -h)) / (2 * h)
        want = (P_A - P_B) * (v - vbar) / HBAR
        worst = max(worst, abs(d - want))
        print(f"  d(mu)/dt at v = {v:6.3f}: {d:+.12f}  expect {want:+.12f}")
    print(f"  worst deviation = {worst:.3e}")
    print(f"  vbar_pair = {vbar:.12f}")
    print()


# --------------------------------------------------------------------- #
# Part C -- the exchange theorem                                        #
# --------------------------------------------------------------------- #
def part_c() -> None:
    print("Part C: stationarity + momentum conservation => momentum swap")

    # Unknowns (p_in, p_out) for the excess particle.
    #   momentum conservation: the struck partner a -> p_b, so it loses
    #       (p_a - p_b); the particle gains it:   p_out - p_in = p_a - p_b
    #   stationarity of mu on the transition path (velocity (p_in+p_out)/2m
    #       equal to vbar_pair = (p_a+p_b)/2m):   p_in + p_out = p_a + p_b
    A = np.array([[-1.0, 1.0], [1.0, 1.0]])
    rhs = np.array([P_A - P_B, P_A + P_B])
    p_in, p_out = np.linalg.solve(A, rhs)
    print(f"  solved p_in  = {p_in:.12f}   (p_b = {P_B})")
    print(f"  solved p_out = {p_out:.12f}   (p_a = {P_A})")
    print(f"  residual     = {max(abs(p_in - P_B), abs(p_out - P_A)):.3e}")
    print("  => the excess particle and the struck partner EXCHANGE momenta")

    # both conservation laws then hold identically, not by imposition
    p_before = p_in + P_A + P_B
    p_after = p_out + P_B + P_B
    e_before = energy(p_in) + energy(P_A) + energy(P_B)
    e_after = energy(p_out) + energy(P_B) + energy(P_B)
    print(f"  momentum before/after: {p_before:.12f} / {p_after:.12f}"
          f"   |diff| = {abs(p_before - p_after):.3e}")
    print(f"  energy   before/after: {e_before:.12f} / {e_after:.12f}"
          f"   |diff| = {abs(e_before - e_after):.3e}")

    # the exchanged momentum is exactly the mode quantum
    print(f"  |p_a - p_b| = {abs(P_A - P_B):.12f}   mode quantum 2q*dp = "
          f"{QUANTUM:.12f}")

    # off-swap assignments dephase: the bias integral decays away from
    # stationarity, with the transition midpoint detuned from vbar_pair
    tau_e = 60.0
    vbar = 0.5 * (P_A + P_B) / MASS
    print("  time-averaged bias |(1/tau) Int_0^tau cos mu dt| vs detuning,")
    print("  against the analytic envelope |sinc(rate*tau/2)|:")
    for dv in (0.0, 0.1, 0.25, 0.5, 1.0):
        v_mid = vbar + dv
        t = np.linspace(0.0, tau_e, 200001)
        rate = (P_A - P_B) * (v_mid - vbar) / HBAR
        val = np.abs(np.trapezoid(np.cos(1.0 + rate * t), t) / tau_e)
        half = 0.5 * rate * tau_e
        env = 1.0 if half == 0.0 else abs(np.sin(half) / half)
        print(f"    detuning {dv:+.2f}   bias = {val:.6f}"
              f"   envelope <= {env:.6f}")
    print()


# --------------------------------------------------------------------- #
# Part D -- the pump makes mu a function of place, not of pair          #
# --------------------------------------------------------------------- #
def part_d() -> None:
    print("Part D: every pumped pair at the same place carries the same mu")

    rng = np.random.default_rng(20260726)
    n_pairs = 60000
    x0 = rng.uniform(0.0, L, n_pairs)      # where each pair sits
    t = 0.37
    x_star = 1.13                          # where the vertex happens

    # the pump writes onto each pair a sideband phase s*K*x0; the partners
    # then transport their clocks to x_star.  Lemma 2 of the phase-resonance
    # note: the pump phase exactly cancels the propagation offset.
    for s in (+1, -1):
        vbar = 0.5 * (P_A + P_B) / MASS
        mu = s * K * x0 + np.pi / 2.0 \
            + s * K * (x_star - x0) - s * K * vbar * t
        spread = np.ptp(np.angle(np.exp(1j * mu)))
        print(f"  family s = {s:+d}: spread of mu over {n_pairs} pairs "
              f"= {spread:.3e}")
    print()


# --------------------------------------------------------------------- #
# Part E -- quadrature and the linear cross term                        #
# --------------------------------------------------------------------- #
def part_e() -> None:
    print("Part E: quadrature and linearity")

    x = np.linspace(0.0, L, 4001, endpoint=False)
    dx = x[1] - x[0]
    V = -V_P * np.cos(K * x)
    dV = (np.roll(V, -1) - np.roll(V, +1)) / (2.0 * dx)   # periodic central
    mu = -K * x + np.pi / 2.0
    bias = np.cos(mu)
    scale = -np.dot(bias, dV) / np.dot(dV, dV)
    resid = np.max(np.abs(bias + scale * dV)) / np.max(np.abs(bias))
    print(f"  cos(mu(x)) vs -V'(x): max relative residual = {resid:.3e}")

    g0, g1, tau_e = 0.30, 0.45, 1.0
    p0 = np.sin(g0 * tau_e) ** 2
    print("  contact vertex: (P - P0) / [tau sin(2 g0 tau) g1 C] -> 1")
    for c in (1e-1, 1e-2, 1e-3, 1e-4):
        p = np.sin(abs(g0 + g1 * c) * tau_e) ** 2
        ratio = (p - p0) / (tau_e * np.sin(2 * g0 * tau_e) * g1 * c)
        print(f"    C = {c:7.0e}   ratio = {ratio:.9f}")
    print()


# --------------------------------------------------------------------- #
# Part F -- stencil conserves L1                                        #
# --------------------------------------------------------------------- #
def part_f() -> None:
    print("Part F: the stencil preserves L1 on the ring")
    rng = np.random.default_rng(7)
    W = rng.random((16, 64))
    d = np.roll(W, +Q, axis=0) - np.roll(W, -Q, axis=0)
    print(f"  max_x |sum_n (W_(n+q) - W_(n-q))| = "
          f"{np.max(np.abs(d.sum(axis=0))):.3e}")
    print()


# --------------------------------------------------------------------- #
# Figure                                                                #
# --------------------------------------------------------------------- #
C_POS = "#B45309"
C_NEG = "#0F766E"
C_MU = "#6D28D9"
C_LINE = "#1F2937"
C_NET = "#B91C1C"
C_GROSS = "#C7CBD1"
C_TEXT = "#111827"
C_SOFT = "#6B7280"
C_FAINT = "#E5E7EB"


def _dial(ax, cx, cy, r, ang_pos, ang_neg, lw=2.6, zorder=5, gap=True):
    ax.add_patch(Circle((cx, cy), r, facecolor="white", edgecolor="#D1D5DB",
                        lw=1.1, zorder=zorder))
    if gap and abs(np.angle(np.exp(1j * (ang_pos - ang_neg)))) > 0.12:
        ax.add_patch(Wedge((cx, cy), r * 0.86, np.degrees(ang_neg),
                           np.degrees(ang_pos), facecolor=C_MU, alpha=0.20,
                           edgecolor="none", zorder=zorder + 1))
    ax.plot([cx, cx + r * 0.86 * np.cos(ang_pos)],
            [cy, cy + r * 0.86 * np.sin(ang_pos)],
            color=C_POS, lw=lw, solid_capstyle="round", zorder=zorder + 2)
    ax.plot([cx, cx + r * 0.86 * np.cos(ang_neg)],
            [cy, cy + r * 0.86 * np.sin(ang_neg)],
            color=C_NEG, lw=lw * 0.55, solid_capstyle="round",
            zorder=zorder + 3)
    ax.add_patch(Circle((cx, cy), r * 0.07, facecolor=C_LINE,
                        edgecolor="none", zorder=zorder + 4))


def _strip(ax):
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def _title(ax, txt):
    ax.set_title(txt, loc="left", fontsize=11.5, fontweight="bold",
                 color=C_TEXT, pad=9)


def panel_a(ax):
    ax.plot([0.3, 5.7], [1.0, 1.0], color=C_FAINT, lw=2.2, zorder=1)
    for cx, ang, col, lab, filled in [(1.15, 1.15, C_POS, "positon $a$", True),
                                      (4.85, 2.35, C_NEG, "negaton $b$", False)]:
        ax.add_patch(Circle((cx, 1.0), 0.13,
                            facecolor=col if filled else "white",
                            edgecolor=col, lw=2.0, zorder=4))
        ax.text(cx, 0.60, lab, fontsize=9.5, ha="center", color=col)
        _dial(ax, cx, 2.15, 0.42, ang, ang, gap=False)
    ax.text(1.15, 2.78, r"$\theta_a$", fontsize=10.5, ha="center", color=C_POS)
    ax.text(4.85, 2.78, r"$\theta_b$", fontsize=10.5, ha="center", color=C_NEG)
    ax.add_patch(Circle((3.0, 1.0), 0.10, facecolor=C_LINE, edgecolor="none",
                        zorder=4))
    ax.text(3.0, 0.60, r"$x^{*}$", fontsize=10.5, ha="center", color=C_LINE)
    for cx in (1.15, 4.85):
        ax.add_patch(FancyArrowPatch((cx, 1.30), (3.0, 1.30),
                                     connectionstyle="arc3,rad=-0.32",
                                     arrowstyle="-|>", mutation_scale=11,
                                     color=C_SOFT, lw=1.3, ls=(0, (3, 2)),
                                     zorder=3))
    ax.text(3.0, 1.92, "transport both clocks\nto the same point",
            fontsize=8.6, ha="center", color=C_SOFT, style="italic")
    _dial(ax, 3.0, 3.62, 0.62, 1.15, 2.35)
    ax.text(3.0, 4.46, r"$\mu \;=\; \Phi_a(x^{*}) - \Phi_b(x^{*})$",
            fontsize=12, ha="center", color=C_MU, fontweight="bold")
    ax.text(3.0, -0.28,
            r"$\Phi_j(x,t) = \theta_j + [p_j (x - x_j) - E_j t]/\hbar$",
            fontsize=10.5, ha="center", color=C_TEXT)
    ax.text(3.0, -0.86,
            "one clock per world-particle;  a pair has exactly one\n"
            "relational number.  There is no third object.",
            fontsize=8.8, ha="center", color=C_SOFT, style="italic")
    ax.set_xlim(0.0, 6.0)
    ax.set_ylim(-1.25, 4.95)
    _title(ax, r"(a)  a pair has one number:  the misalignment $\mu$")
    _strip(ax)


def panel_b(ax):
    xs = np.linspace(0.75, 5.05, 6)
    rows = [
        (3.55, lambda i: 0.0, "aligned",
         "dark:  contributions cancel exactly", C_LINE),
        (2.15, lambda i: 0.95, "uniformly offset",
         r"was 'gray':  same $\mu$ everywhere", C_SOFT),
        (0.75, lambda i: 0.62 * i, "winding",
         r"was 'beating':  $\mu$ turns at $\Delta p/\hbar$", C_MU),
    ]
    for y, off, name, note, col in rows:
        for i, x in enumerate(xs):
            _dial(ax, x, y, 0.36, 1.25 + off(i), 1.25)
        ax.text(5.62, y + 0.20, name, fontsize=10.5, ha="left", color=col,
                fontweight="bold")
        ax.text(5.62, y - 0.20, note, fontsize=8.5, ha="left", color=C_SOFT)
    ax.annotate("", xy=(5.25, 0.10), xytext=(0.65, 0.10),
                arrowprops=dict(arrowstyle="-|>", color=C_FAINT, lw=2.4))
    ax.text(2.95, -0.20, "position along the ring", fontsize=9,
            ha="center", color=C_SOFT)
    ax.text(4.85, -0.95,
            "the third row is what the phase-resonance note calls a beat.\n"
            "Nothing propagates and nothing is quantised: the two clocks\n"
            "simply run at different rates, so their gap turns with $x$.",
            fontsize=8.6, ha="center", color=C_SOFT, style="italic")
    ax.set_xlim(0.3, 9.4)
    ax.set_ylim(-1.55, 4.55)
    _title(ax, r"(b)  the three states a pair can be in")
    _strip(ax)


def panel_c(ax):
    for cx, mu in zip([1.25, 3.0, 4.75], [0.0, 1.9, 3.6]):
        for k in range(3):
            _dial(ax, cx, 3.95 - 0.82 * k, 0.33, 1.15 + mu, 1.15)
        ax.text(cx, 1.72, r"$\mu = %.1f$" % mu, fontsize=9.5, ha="center",
                color=C_MU)
    for lab, y in [("pair 1", 3.95), ("pair 2", 3.13), ("pair 3", 2.31)]:
        ax.text(0.42, y, lab, fontsize=8.4, ha="right", color=C_SOFT)
    x = np.linspace(0, L, 600)
    xs = 0.5 + x * (5.5 / L)
    ax.plot(xs, 0.55 + 0.52 * np.sin(K * x), color=C_NET, lw=2.4, zorder=3)
    ax.plot(xs, 0.55 - 0.52 * np.cos(K * x), color="#9CA3AF", lw=1.8,
            ls=(0, (6, 3)), zorder=2)
    ax.axhline(0.55, color=C_FAINT, lw=0.9, zorder=1)
    ax.text(5.72, 0.98, r"exchange bias $\propto \cos\mu(x)$", fontsize=9,
            color=C_NET, ha="left")
    ax.text(5.72, 0.12, r"$V(x)$", fontsize=9, color="#9CA3AF", ha="left")
    ax.text(3.0, -0.62,
            "every pair sitting at the same place carries the same $\\mu$.\n"
            "The arriving particle never needs to know which pair it met,\n"
            "or where the absent partner is.",
            fontsize=8.8, ha="center", color=C_SOFT, style="italic")
    ax.set_xlim(-0.55, 9.1)
    ax.set_ylim(-1.35, 4.75)
    _title(ax, r"(c)  the pump sets $\mu$ by place, not by pair")
    _strip(ax)


def panel_d(ax):
    for x0, lab in [(0.35, "before"), (3.55, "after")]:
        ax.add_patch(Rectangle((x0, 0.35), 2.75, 3.55, facecolor="#FBFBFC",
                               edgecolor=C_FAINT, lw=1.2, zorder=0))
        ax.text(x0 + 0.12, 3.66, lab, fontsize=9.5, color=C_SOFT,
                style="italic")

    def row(x_dot, y, face, edge, arrow_len, col, label, lab_side="left"):
        ax.add_patch(Circle((x_dot, y), 0.14, facecolor=face, edgecolor=edge,
                            lw=2.0, zorder=4))
        ax.add_patch(FancyArrowPatch((x_dot + 0.17, y),
                                     (x_dot + 0.17 + arrow_len, y),
                                     arrowstyle="-|>", mutation_scale=11,
                                     color=col, lw=2.2, zorder=4))
        if lab_side == "left":
            ax.text(x_dot - 0.30, y, label, fontsize=10, ha="right",
                    color=col, va="center")
        else:
            ax.text(x_dot + 0.24 + arrow_len, y, label, fontsize=9,
                    ha="left", color=col, va="center")

    row(1.05, 2.85, C_LINE, C_LINE, 0.55, C_LINE, r"$c$")
    row(1.05, 1.95, C_POS, C_POS, 1.00, C_POS, r"$a$")
    row(1.05, 1.15, "white", C_NEG, 0.55, C_NEG, r"$b$")
    ax.text(2.95, 1.15, "(elsewhere)", fontsize=8.0, ha="right", color=C_SOFT)
    _dial(ax, 2.55, 2.85, 0.30, 1.15 + 1.5, 1.15)
    ax.text(2.55, 2.38, r"read $\mu$", fontsize=9, ha="center", color=C_MU)

    row(4.25, 2.85, C_LINE, C_LINE, 1.00, C_NET, r"$c$")
    row(4.25, 1.95, C_POS, C_POS, 0.55, C_POS, r"$a$")
    row(4.25, 1.15, "white", C_NEG, 0.55, C_NEG, r"$b$")
    _dial(ax, 5.80, 2.85, 0.30, 1.15, 1.15, gap=False)
    ax.text(5.80, 2.38, r"$\mu = 0$", fontsize=9, ha="center", color=C_LINE)
    ax.add_patch(FancyArrowPatch((3.12, 2.1), (3.48, 2.1), arrowstyle="-|>",
                                 mutation_scale=14, color=C_SOFT, lw=2.0))

    ax.text(3.2, -0.16,
            r"the particle and the struck partner $exchange$ momenta:",
            fontsize=9.6, ha="center", color=C_TEXT)
    ax.text(3.2, -0.66,
            r"$p_c : p_b \rightarrow p_a$,   $p_a \rightarrow p_b$",
            fontsize=11, ha="center", color=C_NET)
    ax.text(3.2, -1.22,
            r"$P(\mathrm{exchange}) \;=\; w_0 \;+\; \kappa\,\cos\mu$",
            fontsize=11, ha="center", color=C_TEXT)
    ax.text(3.2, -1.80,
            "the bare rate $w_0$ is there anyway; $\\mu$ only tilts it.\n"
            "Nothing is absorbed and nothing is emitted.",
            fontsize=8.8, ha="center", color=C_SOFT, style="italic")
    ax.set_xlim(0.05, 6.60)
    ax.set_ylim(-2.25, 4.10)
    _title(ax, r"(d)  the vertex: read $\mu$, exchange, re-align")
    _strip(ax)


def panel_e(ax):
    t = np.linspace(0, 1, 300)
    for y0, rate, lab, col, note in [
            (3.05, 0.0, r"$\mu$ holds still", C_NET, "the tilt accumulates"),
            (0.95, 9.0, r"$\mu$ turns", C_SOFT, "the tilt averages away")]:
        mu = 1.0 + rate * t
        ax.plot(0.7 + 4.0 * t, y0 + 0.34 * np.cos(mu), color=col, lw=2.4,
                zorder=3)
        ax.axhline(y0, xmin=0.09, xmax=0.62, color=C_FAINT, lw=0.9, zorder=1)
        ax.fill_between(0.7 + 4.0 * t, y0, y0 + 0.34 * np.cos(mu),
                        color=col, alpha=0.16, zorder=2)
        ax.text(0.55, y0 + 0.80, lab, fontsize=10.5, color=col,
                fontweight="bold", ha="left")
        ax.text(4.95, y0, note, fontsize=8.8, color=col, ha="left",
                va="center")
        ax.text(0.55, y0 - 0.82, r"$\cos\mu$ along the exchange",
                fontsize=8.2, color=C_SOFT, ha="left")
    ax.text(2.7, -0.35,
            r"$\dot\mu \;=\; \frac{\Delta p}{\hbar}\,"
            r"(v_{\mathrm{exch}} - \bar v_{\mathrm{pair}})$",
            fontsize=11.5, ha="center", color=C_TEXT)
    ax.text(3.95, -1.35,
            r"$\dot\mu = 0$  with momentum conservation"
            "\n"
            r"forces  $p_{\mathrm{in}} = p_b$,  $p_{\mathrm{out}} = p_a$",
            fontsize=10, ha="center", color=C_NET)
    ax.text(3.95, -2.25,
            "one condition, two lines of algebra.  This is all that\n"
            "'resonance', 'row selection' and 'energy conservation'\n"
            "ever were -- and energy conservation is now automatic.",
            fontsize=8.6, ha="center", color=C_SOFT, style="italic")
    ax.set_xlim(-0.15, 8.05)
    ax.set_ylim(-2.75, 4.35)
    _title(ax, r"(e)  the only condition:  $\mu$ must hold still")
    _strip(ax)


def panel_f(ax):
    for lbl, y in {"$n+q$": 2.0, "$n$": 1.0, "$n-q$": 0.0}.items():
        ax.plot([0.6, 7.4], [y, y], color="#D1D5DB", lw=1.1, zorder=1)
        ax.text(0.28, y, lbl, fontsize=9.5, ha="right", va="center",
                color=C_SOFT)
    for xc, sgn, tag in [(2.4, +1, r"$\cos\mu(x) > 0$"),
                         (5.6, -1, r"$\cos\mu(x) < 0$")]:
        ax.add_patch(Rectangle((xc - 1.25, -0.42), 2.5, 2.86,
                               facecolor=C_MU, alpha=0.055, edgecolor="none",
                               zorder=0))
        for dxo in (-0.62, 0.62):
            ax.add_patch(FancyArrowPatch((xc + dxo, 0.06), (xc + dxo, 1.94),
                                         arrowstyle="<|-|>", mutation_scale=9,
                                         color=C_GROSS, lw=3.2, zorder=2))
        y0, y1 = (0.10, 1.90) if sgn > 0 else (1.90, 0.10)
        ax.add_patch(FancyArrowPatch((xc, y0), (xc, y1), arrowstyle="-|>",
                                     mutation_scale=13, color=C_NET, lw=2.6,
                                     zorder=4))
        ax.text(xc, 2.62, tag, fontsize=9.5, ha="center", color=C_MU)
    ax.text(7.62, 2.05, "bare traffic\n(cancels in the mean)", fontsize=8.0,
            color=C_SOFT, ha="left", va="center")
    ax.text(4.0, -1.05,
            r"$\partial_t W_n \;=\; -\frac{p_n}{m}\partial_x W_n \;+\; "
            r"\Gamma_q(x)\,(W_{n+q} - W_{n-q})$",
            fontsize=12, ha="center", color=C_TEXT)
    ax.text(4.0, -1.78,
            "worldlines stream freely; the tilt in their exchanges\n"
            "with the pairs they pass is the whole collision term.",
            fontsize=8.8, ha="center", color=C_SOFT, style="italic")
    ax.text(4.0, -2.52,
            r"linear in $V_p$ only while the pairs stay mutually aligned",
            fontsize=8.6, ha="center", color=C_NET, style="italic")
    ax.set_xlim(-1.25, 10.4)
    ax.set_ylim(-2.95, 3.05)
    _title(ax, r"(f)  many worldlines:  the tilt is the collision term")
    _strip(ax)


def make_figure() -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "mathtext.fontset": "dejavusans",
        "figure.facecolor": "white",
        "axes.facecolor": "white",
    })
    fig, axes = plt.subplots(2, 3, figsize=(17.4, 10.6))
    panel_a(axes[0, 0])
    panel_b(axes[0, 1])
    panel_c(axes[0, 2])
    panel_d(axes[1, 0])
    panel_e(axes[1, 1])
    panel_f(axes[1, 2])
    fig.suptitle("The contact interaction in world-particle terms: "
                 "two clocks, one number, a tilted exchange",
                 fontsize=13.5, fontweight="bold", color=C_TEXT, y=0.975)
    fig.text(0.5, 0.938,
             "no beats  -  no resonance  -  no gratings  -  no rows",
             fontsize=9, ha="center", color=C_SOFT)
    fig.tight_layout(rect=[0.005, 0.005, 0.995, 0.925])

    name = "phase_alignment_contact.png"
    fig.savefig(output_path(name), dpi=180)
    dp = docs_path(name)
    if dp:
        fig.savefig(dp, dpi=180)
    plt.close(fig)
    print(f"wrote {name}")


def main() -> None:
    print(__doc__.split("Run with:")[0].strip().splitlines()[0])
    print()
    part_a()
    part_b()
    part_c()
    part_d()
    part_e()
    part_f()
    make_figure()


if __name__ == "__main__":
    main()
