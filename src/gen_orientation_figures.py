"""Figures for ORIENTATION.md §4: the three processes that change the pair
ledger, each drawn twice -- in phase space and in space-time.

Columns:
  1. emissive event       -- a sea pair at the parent's row is ionised;
  2. absorptive event     -- two bodies recombine into a sea pair at the
                             parent's row (catalysed recombination);
  3. contact recombination -- two coincident bodies combine; no parent.

Before drawing, each process is written down as a list of (species,
momentum, +1 created / -1 destroyed) entries and checked against the rules
stated in ORIENTATION.md §4 and §8 and in stochastic_ledger.md §1:
the change in the observable E per row, the changes in N and S, the pair
count P = S + N/2 (Theorem N1), and momentum conservation with the parent
unchanged (postulate (S), Proposition K8). The script stops if any check
fails, so the picture cannot silently disagree with the text.

Run:
    WPMW_OUTPUT=<dir> PYTHONPATH=src python3 -u src/gen_orientation_figures.py
"""
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from wpmwlib.wpmw_utils import output_path, docs_path

NAME = "orientation_pair_processes.png"

# Momenta in units where the parent row is p = 1 and the channel's momentum
# transfer is xi_q = 0.6.  Only the ordering matters for the picture.
P, XI = 1.0, 0.6

POS, NEG, PAIR, PARENT = "positon", "negaton", "sea pair", "parent"
COLOR = {POS: "C0", NEG: "C4", PAIR: "0.45", PARENT: "k"}

# Each process: (title, ledger-change subtitle, list of (species, momentum, sign)).
# sign +1 = comes into existence at the event, -1 = ceases to exist.
# The parent is listed with sign 0: present before and after, unchanged.
PROCESSES = [
    ("emissive event\n(ionisation)",
     [(PARENT, P, 0), (PAIR, P, -1), (POS, P + XI, +1), (NEG, P - XI, +1)]),
    ("absorptive event\n(catalysed recombination)",
     [(PARENT, P, 0), (NEG, P + XI, -1), (POS, P - XI, -1), (PAIR, P, +1)]),
    ("contact recombination\n(no parent)",
     [(POS, P, -1), (NEG, P, -1), (PAIR, P, +1)]),
]


def ledger(entries):
    """Changes produced by one process: E per row, N, S, P, momentum."""
    dE = defaultdict(int)
    dN = dS = 0
    dmom = 0.0
    for sp, p, s in entries:
        if sp == POS:
            dE[p] += s
            dN += s
            dmom += s * p
        elif sp == NEG:
            dE[p] -= s
            dN += s
            dmom += s * p
        elif sp == PAIR:
            dS += s
            dmom += s * 2 * p          # two co-moving members
    dE = {k: v for k, v in dE.items() if v}
    return dict(dE), dN, dS, dS + dN / 2, dmom


def verify():
    expect_E = {P + XI: +1, P - XI: -1}
    rows = []
    ok = True
    for (title, entries), want_E, want_N, want_S in zip(
            PROCESSES, [expect_E, expect_E, {}], [+2, -2, -2], [-1, +1, +1]):
        dE, dN, dS, dP, dmom = ledger(entries)
        checks = {
            "dE": dE == want_E,
            "dN": dN == want_N,
            "dS": dS == want_S,
            "dP=0": dP == 0,
            "momentum": abs(dmom) < 1e-12,
        }
        ok &= all(checks.values())
        rows.append((title.split("\n")[0], dE, dN, dS, dP, dmom, checks))
    print(f"{'process':<24}{'dE by row':<28}{'dN':>4}{'dS':>4}{'dP':>6}"
          f"{'dMom':>8}  checks")
    for name, dE, dN, dS, dP, dmom, checks in rows:
        dEs = ", ".join(f"{v:+d}@{k:.1f}" for k, v in sorted(dE.items())) or "0"
        flags = " ".join(f"{k}:{'ok' if v else 'FAIL'}" for k, v in checks.items())
        print(f"{name:<24}{dEs:<28}{dN:>+4}{dS:>+4}{dP:>+6.1f}{dmom:>+8.2f}  {flags}")
    return ok


def marker(ax, sp, x, y, hollow=False):
    """Draw one population member in the phase-space panel."""
    if sp == PAIR:
        ax.plot(x, y, marker="s", ms=24, mfc="white" if hollow else "0.85",
                mec="0.45", mew=1.4, ls="", zorder=3)
        ax.text(x, y, "+−", ha="center", va="center", fontsize=9,
                color="0.6" if hollow else "k", zorder=4)
    elif sp == PARENT:
        ax.plot(x, y, marker="*", ms=15, color="k", zorder=5)
    else:
        c = COLOR[sp]
        ax.plot(x, y, "o", ms=15, mfc="white" if hollow else c, mec=c,
                mew=1.6, ls="--" if hollow else "-", zorder=4)
        ax.text(x, y, "+" if sp == POS else "−", ha="center", va="center",
                color=c if hollow else "white", fontsize=10, zorder=5)


def phase_panel(ax, entries):
    rows = [P - XI, P, P + XI]
    for y in rows:
        ax.axhline(y, color="0.88", lw=1, zorder=0)
    ax.axvline(0.5, color="0.8", lw=1, ls=":")
    xb, xa = 0.22, 0.78
    for sp, p, s in entries:
        # the parent sits just left of its row's centre so a co-located
        # sea pair can still be seen
        dx = -0.14 if sp == PARENT else 0.0
        # contact case: two bodies in one cell, drawn side by side
        if sp in (POS, NEG) and all(q == P for _, q, _ in entries):
            dx = -0.06 if sp == POS else 0.06
        if s <= 0:
            marker(ax, sp, xb + dx, p, hollow=(s < 0))
        if s >= 0:
            marker(ax, sp, xa + dx, p)
    ax.set_xlim(0, 1)
    ax.set_ylim(P - XI - 0.35, P + XI + 0.35)
    ax.set_yticks(rows)
    ax.set_yticklabels([r"$p-\xi_q$", r"$p$", r"$p+\xi_q$"])
    ax.set_xticks([xb, xa])
    ax.set_xticklabels(["before", "after"])
    for sd in ("top", "right"):
        ax.spines[sd].set_visible(False)


def spacetime_panel(ax, entries):
    t_ev = 0.5
    # ambient sea: faint straight worldlines at assorted momenta
    for x0, sl in [(-0.9, 0.35), (-0.6, -0.5), (0.5, 0.9), (0.9, -0.2),
                   (-0.2, 1.5), (0.7, 0.2)]:
        ax.plot([x0, x0 + sl], [0, 1], color="0.88", lw=1, zorder=0)
    for sp, p, s in entries:
        sl = p / 2.0                           # dx/dt, with m = 2 for scale
        c = COLOR[sp]
        if sp == PARENT:
            ax.plot([-sl * t_ev, sl * (1 - t_ev)], [0, 1], color="k",
                    lw=1.4, ls=(0, (5, 3)), zorder=4)
            continue
        lw = 4.0 if sp == PAIR else 2.4
        if s < 0:
            seg = ([-sl * t_ev, 0], [0, t_ev])
        else:
            seg = ([0, sl * (1 - t_ev)], [t_ev, 1])
        if sp == NEG and all(q == P for _, q, _ in entries):
            # coincident with the positon: stripe it over the blue line
            ax.plot(*seg, color=c, lw=lw, ls=(0, (4, 4)), zorder=4)
            continue
        ax.plot(*seg, color=c, lw=lw, solid_capstyle="round", zorder=3)
    ax.plot(0, t_ev, "o", ms=8, mfc="white", mec="k", mew=1.5, zorder=6)
    ax.set_xlim(-0.75, 0.75)
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("position $x$")
    ax.set_ylabel("time $t$")
    for sd in ("top", "right"):
        ax.spines[sd].set_visible(False)


def main():
    print("Ledger checks for the three processes "
          f"(parent row p = {P}, transfer xi_q = {XI})\n")
    if not verify():
        raise SystemExit("\nA ledger check failed; not drawing.")
    print("\nAll checks pass; drawing.")

    fig, axs = plt.subplots(2, 3, figsize=(14, 8.6),
                            gridspec_kw=dict(height_ratios=[1, 1.15]))
    for j, (title, entries) in enumerate(PROCESSES):
        axs[0, j].set_title(title, fontsize=11)
        phase_panel(axs[0, j], entries)
        spacetime_panel(axs[1, j], entries)
    axs[0, 0].set_ylabel("momentum row\n(at the event's position)")

    subs = [r"$\Delta N=+2,\ \Delta S=-1$;  $E$: $+1$ at $p+\xi_q$, $-1$ at $p-\xi_q$",
            r"$\Delta N=-2,\ \Delta S=+1$;  $E$: $+1$ at $p+\xi_q$, $-1$ at $p-\xi_q$",
            r"$\Delta N=-2,\ \Delta S=+1$;  $E$ unchanged"]
    for j, s in enumerate(subs):
        axs[0, j].text(0.5, -0.2, s, transform=axs[0, j].transAxes,
                       ha="center", va="top", fontsize=9, color="0.25")

    handles = [
        plt.Line2D([], [], marker="o", ls="", mfc="C0", mec="C0", ms=10,
                   label="positon body"),
        plt.Line2D([], [], marker="o", ls="", mfc="C4", mec="C4", ms=10,
                   label="negaton body"),
        plt.Line2D([], [], marker="s", ls="", mfc="0.85", mec="0.45", ms=10,
                   label="sea pair"),
        plt.Line2D([], [], marker="*", ls="", color="k", ms=12,
                   label="parent world"),
        plt.Line2D([], [], marker="o", ls="", mfc="white", mec="0.45", ms=10,
                   label="hollow: ceases to exist"),
        plt.Line2D([], [], color="k", lw=1.4, ls=(0, (5, 3)),
                   label="parent worldline (unchanged)"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=6, frameon=False,
               fontsize=9, bbox_to_anchor=(0.5, -0.005))
    fig.suptitle("The three processes that change the pair ledger — "
                 "top: momentum rows before and after; "
                 "bottom: worldlines, where momentum is slope",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0.04, 1, 0.96))
    fig.subplots_adjust(hspace=0.42)

    fig.savefig(output_path(NAME), dpi=150, bbox_inches="tight")
    dp = docs_path(NAME)
    if dp:
        fig.savefig(dp, dpi=150, bbox_inches="tight")
    print(f"\nFigure: {NAME}")


if __name__ == "__main__":
    main()
