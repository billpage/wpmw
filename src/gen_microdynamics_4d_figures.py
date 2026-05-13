"""
Generate schematic figures for the 4D-phase-space supplement.

Companion to ``docs/supplement/phase_space_crystal_lattice_4d_supplement.md``.
Each figure highlights a different aspect of the structural similarities and
differences between

  (A) 2 particles in 1 spatial dimension, joint phase space (x1, p1, x2, p2)
  (B) 1 particle in 2 spatial dimensions,  joint phase space (x , px, y , py)

The diagrams are deliberately schematic (blackboard-style), in the spirit of
the Sozi slides reproduced in ``phase_space_crystal_lattice_supplement.md``.

Output convention
-----------------
Figures are written through ``wpmwlib.wpmw_utils.output_path`` so they land in
``$WPMW_OUTPUT`` (set this to ``/mnt/user-data/outputs`` in the container,
``/kaggle/working`` on Kaggle, or leave unset for ``./output``).  When
``$WPMW_DOCS`` is set to the ``output`` branch worktree, copies also go to
``docs_path(...)`` for committing alongside the document.
"""
from __future__ import annotations

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Circle, Rectangle
from matplotlib.lines import Line2D

from wpmwlib.wpmw_utils import output_path, docs_path

# ----------------------------------------------------------------------- #
# Colour / style conventions, matched to the Sozi figures                 #
# ----------------------------------------------------------------------- #
POSITON   = "#c81e1e"   # red,  as in Sozi
NEGATON   = "#1e9e1e"   # green
MEDIATOR  = "#0050c0"   # blue circle, as on Sozi slide 8
AXISGREY  = "#444444"
LATTICE   = "#aaaaaa"

JUMP_LW   = 2.4
ARROW_KW  = dict(arrowstyle="-|>", mutation_scale=14,
                 lw=JUMP_LW, color=POSITON)
DASH_KW   = dict(linestyle=(0, (5, 4)), color=AXISGREY, lw=1.0)


def _save(fig, name: str) -> None:
    """Write a figure to both ``output_path`` and ``docs_path`` (if set)."""
    path = output_path(name)
    fig.savefig(path, dpi=150, bbox_inches="tight",
                facecolor="white")
    dp = docs_path(name)
    if dp:
        fig.savefig(dp, dpi=150, bbox_inches="tight",
                    facecolor="white")
    plt.close(fig)
    print(f"  wrote {path}")
    if dp:
        print(f"  wrote {dp}")


# ----------------------------------------------------------------------- #
# Helpers                                                                 #
# ----------------------------------------------------------------------- #
def _square_axes(ax, xmin, xmax, ymin, ymax, xlabel, ylabel, title):
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect("equal")
    ax.set_xlabel(xlabel, fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_title(title, fontsize=12, pad=8)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_color(AXISGREY)
    ax.spines["bottom"].set_color(AXISGREY)
    ax.tick_params(colors=AXISGREY)


def _origin_cross(ax, x0=0, y0=0):
    ax.axhline(y0, color=AXISGREY, lw=0.6, alpha=0.6)
    ax.axvline(x0, color=AXISGREY, lw=0.6, alpha=0.6)


def _lattice_dots(ax, xs, ys, color=LATTICE, s=8, alpha=0.7):
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    ax.scatter(X.flatten(), Y.flatten(), s=s, color=color,
               alpha=alpha, zorder=1)


def _arrow(ax, p0, p1, color=POSITON, lw=JUMP_LW, **kw):
    a = FancyArrowPatch(p0, p1, arrowstyle="-|>",
                        mutation_scale=14, lw=lw, color=color,
                        zorder=4, **kw)
    ax.add_patch(a)


# ======================================================================= #
# Figure 1 – Phase-space layouts                                          #
# ======================================================================= #
def fig_phase_space_layouts():
    """The two systems both have 4D phase space, but the natural way to
    visualise that 4D space is very different.
    """
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.0))

    # --- (A) 2 particles in 1D --------------------------------------- #
    ax = axes[0]
    # Position line and momentum line
    ax.plot([-2.5, 2.5], [1.2, 1.2],   color=AXISGREY, lw=1.5)
    ax.plot([-2.5, 2.5], [-1.2, -1.2], color=AXISGREY, lw=1.5)
    ax.text(-2.75, 1.2,  r"$x$-axis", ha="right", va="center", fontsize=10,
            color=AXISGREY)
    ax.text(-2.75, -1.2, r"$p$-axis", ha="right", va="center", fontsize=10,
            color=AXISGREY)
    # Particle 1 (red), particle 2 (orange) — positions on x-axis
    p1_pos, p2_pos = -1.2, 1.0
    ax.scatter([p1_pos], [1.2], s=160, color=POSITON, zorder=5)
    ax.scatter([p2_pos], [1.2], s=160, color="#e08020", zorder=5)
    # Coordinate labels at positions
    ax.text(p1_pos, 1.55, r"$x_1$", color=POSITON, ha="center", fontsize=11)
    ax.text(p2_pos, 1.55, r"$x_2$", color="#e08020", ha="center", fontsize=11)
    # Momentum arrows
    _arrow(ax, (p1_pos, -1.2), (p1_pos, -0.5), color=POSITON)
    _arrow(ax, (p2_pos, -1.2), (p2_pos, -1.9), color="#e08020")
    ax.text(p1_pos + 0.08, -0.85, r"$p_1$", color=POSITON, ha="left",
            fontsize=11)
    ax.text(p2_pos + 0.08, -1.55, r"$p_2$", color="#e08020", ha="left",
            fontsize=11)
    # Connection lines from positions to momenta (gentle dashed)
    ax.plot([p1_pos, p1_pos], [1.2, -1.2], color=POSITON,
            linestyle=(0, (2, 4)), lw=0.6, alpha=0.5)
    ax.plot([p2_pos, p2_pos], [1.2, -1.2], color="#e08020",
            linestyle=(0, (2, 4)), lw=0.6, alpha=0.5)
    ax.set_xlim(-3.4, 3.0)
    ax.set_ylim(-2.6, 2.4)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(r"(A) Two particles in 1D"
                 "\n"
                 r"joint phase space $(x_1, p_1, x_2, p_2) \in \mathbb{R}^4$",
                 fontsize=11)

    # --- (B) 1 particle in 2D ---------------------------------------- #
    ax = axes[1]
    # Position plane (x, y) — left box
    pos_box = Rectangle((-2.5, -1.7), 2.4, 2.4, fill=False,
                        edgecolor=AXISGREY, lw=1.4)
    ax.add_patch(pos_box)
    ax.text(-1.3, 0.85, r"position plane $(x, y)$",
            color=AXISGREY, fontsize=10, ha="center")
    # particle at some (x,y)
    px, py = -1.5, -0.6
    ax.scatter([px], [py], s=160, color=POSITON, zorder=5)
    # show x and y components
    ax.plot([px, px], [-1.7, py], color=POSITON, linestyle=(0, (2, 3)),
            lw=0.6, alpha=0.5)
    ax.plot([-2.5, px], [py, py], color=POSITON, linestyle=(0, (2, 3)),
            lw=0.6, alpha=0.5)
    ax.text(px, -1.95, r"$x$", color=POSITON, fontsize=10, ha="center")
    ax.text(-2.75, py, r"$y$", color=POSITON, fontsize=10,
            ha="right", va="center")

    # Momentum plane (px, py) — right box
    mom_box = Rectangle((0.3, -1.7), 2.4, 2.4, fill=False,
                        edgecolor=AXISGREY, lw=1.4)
    ax.add_patch(mom_box)
    ax.text(1.5, 0.85, r"momentum plane $(p_x, p_y)$",
            color=AXISGREY, fontsize=10, ha="center")
    # vector p
    cx, cy = 1.5, -0.5
    ax.scatter([cx], [cy], s=50, color=POSITON, zorder=5)
    _arrow(ax, (cx, cy), (cx + 0.6, cy + 0.5), color=POSITON)
    ax.text(cx + 0.7, cy + 0.55, r"$\vec p$", color=POSITON, fontsize=11,
            ha="left", va="center")
    # decompose components on the box edges
    ax.plot([cx + 0.6, cx + 0.6], [-1.7, cy + 0.5],
            color=POSITON, linestyle=(0, (2, 3)), lw=0.6, alpha=0.5)
    ax.plot([0.3, cx + 0.6], [cy + 0.5, cy + 0.5],
            color=POSITON, linestyle=(0, (2, 3)), lw=0.6, alpha=0.5)
    ax.text(cx + 0.6, -1.95, r"$p_x$", color=POSITON, fontsize=10,
            ha="center")
    ax.text(0.1, cy + 0.5, r"$p_y$", color=POSITON, fontsize=10,
            ha="right", va="center")

    ax.set_xlim(-3.4, 3.4)
    ax.set_ylim(-2.6, 2.4)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(r"(B) One particle in 2D"
                 "\n"
                 r"joint phase space $(x, y, p_x, p_y) \in \mathbb{R}^4$",
                 fontsize=11)

    fig.suptitle("Both systems live in a 4-dimensional joint phase space —"
                 " but its 'natural slicing' differs",
                 fontsize=12, y=1.02)
    fig.tight_layout()
    _save(fig, "microdynamics_4d_layouts.png")


# ======================================================================= #
# Figure 2 – Momentum-space jump diagrams                                 #
# ======================================================================= #
def fig_momentum_jumps():
    """The single most important figure: jump direction in 2D momentum space.

    For 2p/1D the jumps are constrained to the anti-diagonal (p1 = -p2 line);
    total momentum p1+p2 is conserved.  For 1p/2D the jumps fan out in
    arbitrary directions set by the Fourier-mode wavevector k_q.
    """
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.6))

    # --- (A) 2p/1D in (p1, p2) plane --------------------------------- #
    ax = axes[0]
    lat = np.arange(-3, 4)
    _lattice_dots(ax, lat, lat, color=LATTICE, s=14, alpha=0.55)
    _origin_cross(ax)

    # Parallel anti-diagonals: each is a level set of P = p1 + p2.
    # Line y = -x + k goes through (-3.4, 3.4+k) -> (3.4, -3.4+k).
    for k in (-4, -2, 0, 2, 4):
        ax.plot([-3.4, 3.4], [3.4 + k, -3.4 + k],
                **DASH_KW, alpha=0.45)
    # highlighted central anti-diagonal (P = 0)
    ax.plot([-3.4, 3.4], [3.4, -3.4], color=AXISGREY, lw=1.4,
            linestyle=(0, (4, 3)))
    ax.text(3.3, -3.2, r"$p_1 + p_2 = \mathrm{const}$",
            color=AXISGREY, fontsize=9, rotation=-45, ha="right",
            va="top")

    # the mediator pair, centred at origin
    ax.scatter([0], [0], s=180, color=MEDIATOR, zorder=6,
               edgecolors="white", lw=1.0)
    ax.text(0.05, -0.6, "pre-jump\n$(p_1, p_2)$", color=MEDIATOR,
            fontsize=9, ha="left", va="top")

    # several jump arrows of different |q|, all along the anti-diagonal
    for q in (1, 2, 3):
        # +ħk_q for particle 1, -ħk_q for particle 2  -> displacement (+q,-q)
        _arrow(ax, (0, 0), (q, -q), color=POSITON)
    ax.text(3.1, -2.7, r"$+\hbar k_q$ to $p_1$,"
                       "\n"
                       r"$-\hbar k_q$ to $p_2$",
            color=POSITON, fontsize=9, ha="left", va="center")

    _square_axes(ax, -3.6, 3.6, -3.6, 3.6,
                 r"$p_1$", r"$p_2$",
                 r"(A) 2p/1D: jumps confined to anti-diagonal")

    # --- (B) 1p/2D in (px, py) plane --------------------------------- #
    ax = axes[1]
    _lattice_dots(ax, lat, lat, color=LATTICE, s=14, alpha=0.55)
    _origin_cross(ax)

    # the mediator positon at the origin
    ax.scatter([0], [0], s=180, color=MEDIATOR, zorder=6,
               edgecolors="white", lw=1.0)
    ax.text(0.05, -0.6, "pre-jump\n$(p_x, p_y)$", color=MEDIATOR,
            fontsize=9, ha="left", va="top")

    # jump arrows for several different k_q directions, length |k_q|
    arrow_modes = [
        (3, 0, r"$\vec q = (3,0)$"),
        (0, 3, r"$\vec q = (0,3)$"),
        (2, 2, r"$\vec q = (2,2)$"),
        (-2, 1, r"$\vec q = (-2,1)$"),
        (1, -3, r"$\vec q = (1,-3)$"),
    ]
    for qx, qy, lbl in arrow_modes:
        _arrow(ax, (0, 0), (qx, qy), color=POSITON)
        ax.text(qx * 1.15, qy * 1.15, lbl, color=POSITON,
                fontsize=8, ha="center", va="center")

    _square_axes(ax, -3.6, 3.6, -3.6, 3.6,
                 r"$p_x$", r"$p_y$",
                 r"(B) 1p/2D: jumps along arbitrary $\hbar\vec k_{\vec q}$")

    fig.suptitle(
        "Per-event momentum-space displacement:  "
        r"2p/1D conserves $P = p_1+p_2$;  1p/2D does not conserve $\vec p$",
        fontsize=11, y=1.00)
    fig.tight_layout()
    _save(fig, "microdynamics_4d_jumps.png")


# ======================================================================= #
# Figure 3 – Fourier-mode structure                                       #
# ======================================================================= #
def fig_fourier_modes():
    """Where the active Fourier modes live in joint wavevector space.

    For 1p/2D the modes are an arbitrary subset of Z^2.
    For 2p/1D the pair potential V2(x1 - x2) gives modes supported only on
    the q1 = -q2 anti-diagonal, even though the joint configuration space is
    formally 2D.  This is the structural origin of the centre-of-mass /
    relative-coordinate separation that makes 2p/1D effectively 1+1-D.
    """
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.4))

    # --- (A) 2p/1D pair potential V2(x1 - x2) ------------------------ #
    ax = axes[0]
    qs = np.arange(-4, 5)
    _lattice_dots(ax, qs, qs, color=LATTICE, s=18, alpha=0.55)
    _origin_cross(ax)
    # anti-diagonal: q1 = -q2
    ax.plot([-4.4, 4.4], [4.4, -4.4], color=AXISGREY,
            linestyle=(0, (4, 3)), lw=1.1)
    # highlight a few non-trivial modes on the line
    for q in (1, 2, 3):
        ax.scatter([q], [-q], s=140, color=POSITON, zorder=5,
                   edgecolors="white", lw=1.0)
        ax.scatter([-q], [q], s=140, color=POSITON, zorder=5,
                   edgecolors="white", lw=1.0)
    ax.text(-4.2, -4.0,
            r"only modes with $q_1 = -q_2$ are active"
            "\n"
            r"($V_2$ depends only on $r = x_1 - x_2$)",
            color=AXISGREY, fontsize=9, ha="left", va="bottom")
    _square_axes(ax, -4.4, 4.4, -4.4, 4.4,
                 r"$q_1$", r"$q_2$",
                 r"(A) 2p/1D pair potential $V_2(x_1 - x_2)$"
                 "\n"
                 r"active modes on a 1D line")

    # --- (B) 1p/2D external V(x, y) ---------------------------------- #
    ax = axes[1]
    _lattice_dots(ax, qs, qs, color=LATTICE, s=18, alpha=0.55)
    _origin_cross(ax)
    # arbitrary modes can be present
    illustrative = [(3, 0), (-3, 0), (0, 3), (0, -3),
                    (2, 2), (-2, -2), (2, -2), (-2, 2),
                    (1, 3), (-1, -3), (3, -1), (-3, 1)]
    for qx, qy in illustrative:
        ax.scatter([qx], [qy], s=140, color=POSITON, zorder=5,
                   edgecolors="white", lw=1.0)
    ax.text(-4.2, -4.0,
            r"any subset of $\mathbb{Z}^2$ may be active"
            "\n"
            r"(set by the geometry of $V(x, y)$)",
            color=AXISGREY, fontsize=9, ha="left", va="bottom")
    _square_axes(ax, -4.4, 4.4, -4.4, 4.4,
                 r"$q_x$", r"$q_y$",
                 r"(B) 1p/2D external $V(x, y)$"
                 "\n"
                 r"active modes fill 2D $\mathbb{Z}^2$")

    fig.suptitle(
        "Fourier-mode support: the dimensional reduction of 2p/1D is "
        "visible already in $\\vec q$-space",
        fontsize=11, y=1.02)
    fig.tight_layout()
    _save(fig, "microdynamics_4d_fourier_modes.png")


# ======================================================================= #
# Figure 4 – Mediator starbursts (analog of Sozi slide 14)                #
# ======================================================================= #
def fig_mediator_starburst():
    """A direct generalisation of the Cyganski crystal-lattice 'starburst'.

    Left  (2p/1D): the *two* particles in a single world get correlated kicks
                   of equal magnitude and opposite sign whenever the pair
                   exchanges a virtual quantum of the interaction.
    Right (1p/2D): the single particle gets a vectorial kick that fans out in
                   2D directions, one per Fourier mode of V(x,y).
    """
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.6))

    # --- (A) 2p/1D correlated exchange ------------------------------- #
    ax = axes[0]
    # Particle 1 momentum axis (top half)
    y1, y2 = 1.3, -1.3
    ax.axhline(y1, color=AXISGREY, lw=1.0)
    ax.axhline(y2, color=AXISGREY, lw=1.0)
    ax.text(-3.6, y1 + 0.1, r"$p_1$", color=POSITON, fontsize=11)
    ax.text(-3.6, y2 + 0.1, r"$p_2$", color="#e08020", fontsize=11)

    # tick marks for n-1, n, n+1
    for ax_y in (y1, y2):
        for x in (-2, -1, 0, 1, 2):
            ax.plot([x, x], [ax_y - 0.1, ax_y + 0.1],
                    color=AXISGREY, lw=0.8)

    # central mediator pair (both particles "at rest" relative to themselves)
    ax.scatter([0], [y1], s=160, color=POSITON, zorder=5,
               edgecolors="white", lw=1.0)
    ax.scatter([0], [y2], s=160, color="#e08020", zorder=5,
               edgecolors="white", lw=1.0)
    # mediator "link" between them
    ax.plot([0, 0], [y1, y2], color=MEDIATOR, lw=1.3,
            linestyle=(0, (3, 3)), alpha=0.7)
    ax.text(0.15, 0, "exchange\n$\\hbar k_q$", color=MEDIATOR,
            fontsize=9, ha="left", va="center")

    # correlated jumps: +q on p1, -q on p2  (and the opposite sign too)
    for q, alpha in ((1, 1.0), (2, 0.8), (3, 0.6)):
        _arrow(ax, (0, y1), (+q, y1), color=POSITON, alpha=alpha)
        _arrow(ax, (0, y2), (-q, y2), color="#e08020", alpha=alpha)
        # and the reverse
        _arrow(ax, (0, y1), (-q, y1), color=POSITON, alpha=alpha * 0.55)
        _arrow(ax, (0, y2), (+q, y2), color="#e08020", alpha=alpha * 0.55)

    ax.text(-3.6, -2.45,
            r"particles 1 & 2 jump together, "
            r"with $\Delta p_1 = -\Delta p_2$",
            color=AXISGREY, fontsize=10, ha="left")

    ax.set_xlim(-3.8, 3.5)
    ax.set_ylim(-3.0, 2.4)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(r"(A) 2p/1D: correlated exchange between two particles",
                 fontsize=11, pad=4)

    # --- (B) 1p/2D vectorial starburst ------------------------------- #
    ax = axes[1]
    # Light dotted lattice
    lat = np.arange(-3, 4)
    _lattice_dots(ax, lat, lat, color=LATTICE, s=14, alpha=0.5)
    _origin_cross(ax)

    # Central positon mediator
    ax.scatter([0], [0], s=180, color=MEDIATOR, zorder=6,
               edgecolors="white", lw=1.0)

    # Fan of jump arrows -- a discrete set of modes
    modes = [(3, 0), (2, 2), (0, 3), (-2, 2), (-3, 0),
             (-2, -2), (0, -3), (2, -2)]
    for qx, qy in modes:
        _arrow(ax, (0, 0), (qx, qy), color=POSITON, lw=2.2)

    _square_axes(ax, -3.8, 3.8, -3.8, 3.8,
                 r"$p_x$", r"$p_y$",
                 r"(B) 1p/2D: vectorial momentum kick to one particle")

    fig.suptitle(
        "Crystal-lattice mediator picture — Cyganski-style starburst, "
        "generalised",
        fontsize=11, y=1.00)
    fig.tight_layout()
    _save(fig, "microdynamics_4d_starburst.png")


# ======================================================================= #
# Figure 5 – COM / relative decomposition for 2p/1D                       #
# ======================================================================= #
def fig_com_relative():
    """Rotation of coordinates in the (p1, p2) plane onto (P, p_rel) shows
    that the 2p/1D jumps act only on the relative-momentum axis, so the
    centre-of-mass momentum is exactly inert.  This is the basis of the
    'hidden 1+1-D' character of the 2p/1D problem.
    """
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.4))

    # --- (A) Original (p1, p2) basis --------------------------------- #
    ax = axes[0]
    lat = np.arange(-3, 4)
    _lattice_dots(ax, lat, lat, color=LATTICE, s=14, alpha=0.55)
    _origin_cross(ax)
    # P = const lines and prel = const lines
    for c in (-2, 0, 2):
        ax.plot([-3.4, 3.4], [-3.4 + c, 3.4 + c],
                color=AXISGREY, lw=0.7, linestyle=(0, (3, 3)), alpha=0.5)
        ax.plot([-3.4, 3.4], [3.4 + c, -3.4 + c],
                color=AXISGREY, lw=0.7, linestyle=(0, (3, 3)), alpha=0.5)
    # jump arrows along the anti-diagonal
    for q in (1, 2):
        _arrow(ax, (0, 0), (q, -q), color=POSITON)
        _arrow(ax, (0, 0), (-q, q), color=POSITON, alpha=0.7)
    _square_axes(ax, -3.6, 3.6, -3.6, 3.6,
                 r"$p_1$", r"$p_2$",
                 r"(A) original basis $(p_1, p_2)$")

    # --- (B) Rotated (P, p_rel) basis -------------------------------- #
    ax = axes[1]
    _lattice_dots(ax, lat, lat, color=LATTICE, s=14, alpha=0.55)
    _origin_cross(ax)
    # jumps now lie *along* the prel axis
    for q in (1, 2):
        _arrow(ax, (0, 0), (0, -2 * q), color=POSITON)
        _arrow(ax, (0, 0), (0, +2 * q), color=POSITON, alpha=0.7)
    ax.text(0.2, -3.4, "all jumps along\nthe $p_\\mathrm{rel}$ axis",
            color=POSITON, fontsize=9, ha="left", va="bottom")
    ax.text(3.3, 0.2, r"$P$ totally inert",
            color=AXISGREY, fontsize=9, ha="right", va="bottom")
    _square_axes(ax, -3.6, 3.6, -3.6, 3.6,
                 r"$P = p_1 + p_2$",
                 r"$p_\mathrm{rel} = (p_1 - p_2)/2$",
                 r"(B) rotated basis $(P, p_\mathrm{rel})$")

    fig.suptitle(
        "Why 2p/1D is 'hidden 1+1-D': "
        r"the jump rule sees only $p_\mathrm{rel}$",
        fontsize=11, y=1.00)
    fig.tight_layout()
    _save(fig, "microdynamics_4d_com_relative.png")


# ======================================================================= #
# Main                                                                    #
# ======================================================================= #
def main():
    print("Generating 4D-phase-space supplement figures")
    print("=" * 60)
    fig_phase_space_layouts()
    fig_momentum_jumps()
    fig_fourier_modes()
    fig_mediator_starburst()
    fig_com_relative()
    print("done.")


if __name__ == "__main__":
    main()
