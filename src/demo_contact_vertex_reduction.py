"""
Demonstration: reduction of P5 to a unitary contact vertex.

Companion to §13 of ``docs/analysis/phase_resonance_microdynamics.md``.
The P5 weight is replaced by a two-level unitary contact interaction

    h = g0 + g1 * C * exp(i delta),      P_flip = sin^2(|h| tau_e),

with a bare exchange amplitude g0 and a beat-stimulated amplitude
g1 * C * exp(i delta), delta the struck beat's pattern phase at the vertex.
No phase offset is fitted anywhere. Four checks:

  R1  the rate law of the phase-resonance note re-emerges with the
      analytically predicted coefficient
          coef = tau_e * sin(2 g0 tau_e) * g1 * C
      and the correct quadrature sign: delta_0 = 0 is *derived* from
      hermiticity, not matched.

  R2  the residual scales as O(C^3): the O(C^2) term of |h|^2 is
      direction-symmetric and cancels in the net, so linear response is
      exact at first order.

  R3  g0 -> 0 kills the response ("no noise, no force"): the coefficient
      tracks sin(2 g0 tau_e), so the bare phase-blind exchange traffic is
      the carrier of the quantum force, not removable noise.

  R4  with two pumped modes the gross rate acquires an O(C^2) spatial
      modulation at the difference wavevector, amplitude
      2 g1^2 C2 C3 x Jacobian - the discriminating signature of the
      amplitude vertex (the bare affine P5 predicts exactly zero).

Figure (written via docs_path for the output branch):
  contact_vertex_concept.png - anatomy of the contact, the phasor sum,
  the emerging rate law, and the consequence tree.
"""

from __future__ import annotations

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Circle, FancyArrowPatch  # noqa: E402

from wpmwlib.wpmw_utils import output_path, docs_path  # noqa: E402

HBAR = 1.0
M_PART = 1.0
L = 8.0
DP = np.pi * HBAR / L
N_HI = 6
G0, G1, TAU_E = 0.7, 0.9, 0.5
JAC = TAU_E * np.sin(2 * G0 * TAU_E) / (2 * G0)   # dP/d|h|^2 at |h| = g0
CORAL, TEAL, GREY, INK = "#e2725b", "#0f8b8d", "#888888", "#222222"

RNG = np.random.default_rng(20260722)
N_PAIRS = 40000
X0 = RNG.uniform(0, L, N_PAIRS)


def pattern_phase(x: float, s: int, q: int) -> np.ndarray:
    """Lambda^s at t=0 for mode q: s K x + pi/2 (Lemma 2 of the note)."""
    k = 2 * np.pi * q / L
    return s * k * (x - X0) + (np.pi / 2 + s * k * X0)


def p_flip(x: float, d: int, cs, qs, g0: float = G0) -> float:
    """Unitary-vertex flip probability, direction d, beats (C_k, mode q_k)."""
    h = g0 + sum(G1 * c * np.exp(1j * pattern_phase(x, d, q))
                 for c, q in zip(cs, qs))
    return float(np.mean(np.sin(np.abs(h) * TAU_E) ** 2))


def part_r12() -> None:
    print("=" * 72)
    print("R1/R2: rate law with predicted coefficient; O(C^3) residual")
    q = 2
    k = 2 * np.pi * q / L
    xs = np.linspace(0.1, L - 0.1, 40)
    for c in (1e-2, 1e-3, 1e-4):
        net = np.array([p_flip(x, -1, [c], [q]) - p_flip(x, +1, [c], [q])
                        for x in xs])
        pred = 2 * JAC * (2 * G0 * G1 * c) * np.sin(k * xs)
        resid = float(np.max(np.abs(net - pred)))
        print(f"  C={c:.0e}: max|net - pred| = {resid:.2e}   "
              f"resid/C^3 = {resid / c**3:.3f}")


def part_r3() -> None:
    print("=" * 72)
    print("R3: coefficient tracks sin(2 g0 tau_e); g0 -> 0 kills the force")
    q = 2
    x_probe = L / (4 * q)
    c = 1e-3
    for g0_test in (0.7, 0.35, 0.1, 0.0):
        net = (p_flip(x_probe, -1, [c], [q], g0=g0_test)
               - p_flip(x_probe, +1, [c], [q], g0=g0_test))
        pred = (TAU_E * np.sin(2 * g0_test * TAU_E) * 2 * G1 * c
                if g0_test else 0.0)
        print(f"  g0={g0_test:.2f}: net {net:+.3e}  pred {pred:+.3e}")


def part_r4() -> None:
    print("=" * 72)
    print("R4: two-mode O(C^2) gross-rate modulation at the difference "
          "wavevector")
    q2, q3 = 2, 3
    c2 = c3 = 5e-2
    xs = np.linspace(0.1, L - 0.1, 40)
    gross = np.array([0.5 * (p_flip(x, -1, [c2, c3], [q2, q3])
                             + p_flip(x, +1, [c2, c3], [q2, q3]))
                      for x in xs])
    kd = 2 * np.pi * (q3 - q2) / L
    cc = 2 * float(np.mean((gross - gross.mean()) * np.cos(kd * xs)))
    ss = 2 * float(np.mean((gross - gross.mean()) * np.sin(kd * xs)))
    amp = float(np.hypot(cc, ss))
    pred = JAC * 2 * G1**2 * c2 * c3
    print(f"  measured {amp:.3e}  predicted 2 g1^2 C2 C3 x Jac = {pred:.3e}"
          f"  ratio {amp / pred:.4f}")
    print("  (bare affine P5 predicts exactly zero here)")


def figure_concept() -> None:
    fig = plt.figure(figsize=(13.0, 9.2))
    gs = fig.add_gridspec(2, 2, hspace=0.32, wspace=0.22)

    # (a) contact anatomy
    ax = fig.add_subplot(gs[0, 0])
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis("off")
    ax.set_title("(a) the contact: one Hermitian coupling", fontsize=11)
    ax.plot([0.8, 4.6], [7.5, 7.5], color=INK, lw=2)
    ax.plot([0.8, 4.6], [3.0, 3.0], color=INK, lw=2)
    ax.text(0.7, 7.8, "$|hi\\rangle$:  $p_{hi}$", fontsize=11)
    ax.text(0.7, 2.2, "$|lo\\rangle$:  $p_{hi}-2q\\thinspace dp$",
            fontsize=11)
    ax.add_patch(FancyArrowPatch((2.0, 7.3), (2.0, 3.2), arrowstyle="<->",
                                 color=GREY, lw=2.2, mutation_scale=16))
    ax.text(1.15, 5.15, "$g_0$", color=GREY, fontsize=13)
    ax.text(0.9, 4.55, "bare\nexchange", color=GREY, fontsize=8,
            ha="center")
    ax.add_patch(FancyArrowPatch((3.6, 7.3), (3.6, 3.2), arrowstyle="<->",
                                 color=CORAL, lw=2.2, mutation_scale=16,
                                 connectionstyle="arc3,rad=0.25"))
    ax.text(4.15, 5.15, "$g_1\\thinspace C\\thinspace e^{i\\delta}$",
            color=CORAL, fontsize=13)
    ax.text(4.55, 4.45, "stimulated\nby the beat", color=CORAL, fontsize=8,
            ha="center")
    ax.text(2.7, 8.9, "$h = g_0 + g_1 C e^{i\\delta}$\n"
            "$U = e^{-iH\\tau_e}$,   $P_{flip} = \\sin^2(|h|\\tau_e)$",
            fontsize=11, ha="center",
            bbox=dict(boxstyle="round", fc="#f6f6f6", ec=INK))
    ax.text(7.6, 9.0, "the pair in the cell", fontsize=9, ha="center")
    ax.add_patch(Circle((7.0, 7.6), 0.32, color=CORAL))
    ax.add_patch(Circle((7.45, 7.6), 0.32, fc="none", ec=TEAL, lw=2,
                        ls="--"))
    ax.text(8.15, 7.5, "dark: $\\Psi\\equiv 0$\noffers only $g_0$\n"
            "(exit: emission, K4)", fontsize=8, va="center")
    ax.add_patch(Circle((6.8, 4.6), 0.32, color=CORAL))
    ax.add_patch(Circle((7.75, 4.6), 0.32, fc="none", ec=TEAL, lw=2,
                        ls="--"))
    xw = np.linspace(6.4, 8.2, 100)
    ax.plot(xw, 3.7 + 0.25 * np.sin(2 * np.pi * (xw - 6.4) / 0.9),
            color=TEAL, lw=1.2)
    ax.text(8.55, 4.5, "beating: adds\n$g_1 C e^{i\\delta}$\n"
            "(exit: absorption, K3)", fontsize=8, va="center")
    ax.text(5.0, 0.8, "$\\delta$ = the beat's pattern phase at the vertex"
            "  (Lemma 2:  $\\delta = \\pm Kx + \\pi/2$)",
            fontsize=9, ha="center", color=INK)

    # (b) phasor sum
    ax = fig.add_subplot(gs[0, 1])
    ax.set_title("(b) phasor sum: linearity and $\\delta_0 = 0$ for free",
                 fontsize=11)
    ax.set_xlim(-0.15, 1.35)
    ax.set_ylim(-0.55, 0.75)
    ax.axhline(0, color=GREY, lw=0.6)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    g0v, g1c = 0.85, 0.33
    delta = 0.9
    ax.add_patch(FancyArrowPatch((0, 0), (g0v, 0), arrowstyle="->",
                                 color=GREY, lw=2.6, mutation_scale=18))
    ax.text(0.38, -0.09, "$g_0$ (bare)", color=GREY, fontsize=11)
    tip = (g0v + g1c * np.cos(delta), g1c * np.sin(delta))
    ax.add_patch(FancyArrowPatch((g0v, 0), tip, arrowstyle="->",
                                 color=CORAL, lw=2.6, mutation_scale=18))
    ax.text(g0v + 0.08, 0.21, "$g_1 C e^{i\\delta}$", color=CORAL,
            fontsize=11)
    ax.add_patch(FancyArrowPatch((0, 0), tip, arrowstyle="->", color=INK,
                                 lw=1.8, mutation_scale=16))
    ax.text(0.47, 0.30, "$|h|$", color=INK, fontsize=11)
    th = np.linspace(0, 2 * np.pi, 200)
    ax.plot(g0v + g1c * np.cos(th), g1c * np.sin(th), color=CORAL, lw=0.8,
            ls=":", alpha=0.8)
    ax.annotate("$\\delta$", xy=(g0v + 0.13, 0.075), fontsize=11,
                color=INK)
    ax.text(0.6, -0.38,
            "$|h|^2 = g_0^2 + 2 g_0 g_1 C\\cos\\delta + O(C^2)$\n"
            "both arrows share the unitarity factor $-i$:\n"
            "the cross term is a pure cosine, no offset",
            fontsize=9, ha="center",
            bbox=dict(boxstyle="round", fc="#f6f6f6", ec=GREY))

    # (c) the rate law emerging (contrast exaggerated for visibility)
    ax = fig.add_subplot(gs[1, 0])
    ax.set_title("(c) the rate law emerges (contrast exaggerated)",
                 fontsize=11)
    q = 2
    k = 2 * np.pi * q / L
    c_big = 0.25
    xs = np.linspace(0, L, 400)
    pdn = np.sin(np.abs(G0 + G1 * c_big
                        * np.exp(1j * (-k * xs + np.pi / 2))) * TAU_E) ** 2
    pup = np.sin(np.abs(G0 + G1 * c_big
                        * np.exp(1j * (+k * xs + np.pi / 2))) * TAU_E) ** 2
    ax.plot(xs, pdn, color=CORAL, lw=1.8, label="down-flip $P_{-}(x)$")
    ax.plot(xs, pup, color=TEAL, lw=1.8, label="up-flip $P_{+}(x)$")
    ax.axhline(np.sin(G0 * TAU_E) ** 2, color=GREY, lw=1.0, ls="--",
               label="bare rate $\\sin^2(g_0\\tau_e)$ (G-noise)")
    ax.plot(xs, np.sin(G0 * TAU_E) ** 2 + 0.5 * (pdn - pup), color=INK,
            lw=2.4, label="net bias $\\propto \\sin(Kx) = \\Gamma_q(x)$")
    ax.set_xlabel("x")
    ax.set_ylabel("flip probability per contact")
    ax.legend(fontsize=8, loc="upper right")

    # (d) consequence tree
    ax = fig.add_subplot(gs[1, 1])
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis("off")
    ax.set_title("(d) everything else, consequence", fontsize=11)

    def box(x, y, text, fc="#f6f6f6", ec=INK, fs=8.3):
        ax.text(x, y, text, fontsize=fs, ha="center", va="center",
                bbox=dict(boxstyle="round", fc=fc, ec=ec))

    def arrow(x0, y0, x1, y1):
        ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="->",
                                     color=INK, lw=1.3, mutation_scale=12))

    box(5.0, 9.2, "Hermitian contact  $h = g_0 + g_1 C e^{i\\delta}$,"
        "  window $\\tau_e$", fc="#fdeeea", ec=CORAL, fs=9.5)
    box(1.9, 7.2, "unitarity:\n$\\delta_0=0$, K3$\\leftrightarrow$K4\n"
        "detailed balance,\nsaturation")
    box(5.0, 7.2, "Born on exit:\nP5 affine form\n"
        "$w_0 + \\kappa C\\cos\\delta$")
    box(8.2, 7.2, "$g_0\\neq 0$ required:\nnoise carries\nthe force")
    arrow(4.0, 8.8, 1.9, 7.9)
    arrow(5.0, 8.8, 5.0, 7.8)
    arrow(6.0, 8.8, 8.2, 7.9)
    box(2.6, 4.9, "+ Lemma 2 (pump coherence):\n"
        "$\\delta = \\pm Kx + \\pi/2$\n"
        "$\\Rightarrow$ quadrature $\\Gamma_q \\propto -V'$,"
        " sign $\\sigma$, $\\gamma/2$")
    box(7.4, 4.9, "+ P1 winding over $\\tau_e$:\nresonance $\\Rightarrow$\n"
        "own-row selection (Prop. 2)")
    arrow(3.5, 6.6, 2.8, 5.6)
    arrow(6.2, 6.6, 7.2, 5.6)
    box(5.0, 2.6, "coherent sea ($B$ pairs, Lemma 2)\nlicenses "
        "bare/stimulated interference", fc="#eaf4f4", ec=TEAL)
    box(5.0, 0.9, "QLE collision term (exact, L1)", fc="#fdeeea",
        ec=CORAL, fs=9.5)
    arrow(2.9, 4.1, 4.3, 3.2)
    arrow(7.1, 4.1, 5.9, 3.2)
    arrow(5.0, 2.0, 5.0, 1.4)

    fig.savefig(output_path("contact_vertex_concept.png"), dpi=150,
                bbox_inches="tight")
    dp_fig = docs_path("contact_vertex_concept.png")
    if dp_fig:
        fig.savefig(dp_fig, dpi=150, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    part_r12()
    part_r3()
    part_r4()
    figure_concept()
    print("=" * 72)
    print("figure written: contact_vertex_concept.png")
