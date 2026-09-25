"""Reader-frame readings of the step 22 particle model -- open items L-SP7, L-SP9.

Runs ``demo_sea_lock_particles.py --readers``, three seeds per configuration,
and prints late means (4 <= t <= 25) with standard errors.  Table 1 is the
sea of section 6 under the maintenance regimes of that section and faster
dark catalysis; table 2 crosses three sea variants with four regimes:

  sea variants   (S), continuous   aligned pairs stream under the force,
                                   initial momenta uniform within each row
                                   (the runs of section 6)
                 (S), row-centred  the same, initial momenta at row centres
                 row-keeping       row-centred and no momentum kick for
                                   aligned pairs (--sea-force blind):
                                   a diagnostic, it breaks (S) for the sea
  regimes        streaming only; events with reach-scale dark catalysis at
                 1, 3 and 10 times the kernel's own rate

  columns        static reading in the partners' own frames, the same in the
                 parent's frame (= the exact rate reading A'), the rate
                 reading A from clock rates, the own-frame rate reading,
                 rms drift term / rms U_res, sea coherence |<e^{i mu_ref}>|,
                 median clock age near the barrier, the mu = 0 control.

Writes sea_lock_readers_scan.png (table 2).  About three minutes.

    PYTHONPATH=src python3 -u src/scan_sea_lock_readers.py
"""
import os
import re
import subprocess
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from wpmwlib.wpmw_utils import output_path, docs_path

HERE = os.path.dirname(os.path.abspath(__file__))
DEMO = os.path.join(HERE, "demo_sea_lock_particles.py")
BASE = ["--nu", "8", "--mode", "locked", "--t-end", "25", "--readers"]
SEAS = (("(S), continuous", []),
        ("(S), row-centred", ["--sea-p", "rows"]),
        ("row-keeping", ["--sea-p", "rows", "--sea-force", "blind"]))
REGIMES = (("streaming only", ["--no-events"]),
           ("dark x1", ["--relock-w", "3", "--dark", "reach"]),
           ("dark x3", ["--relock-w", "3", "--dark", "reach", "--dark-rate", "3"]),
           ("dark x10", ["--relock-w", "3", "--dark", "reach", "--dark-rate", "10"]))
SECTION6 = (("streaming only", ["--no-events"]),
            ("co-located dark", ["--relock-w", "3", "--dark", "mean"]),
            ("reach dark x1", ["--relock-w", "3", "--dark", "reach"]),
            ("x1, reader lever", ["--relock-w", "3", "--dark", "reach",
                                  "--lever", "reader"]),
            ("reach dark x3", ["--relock-w", "3", "--dark", "reach",
                               "--dark-rate", "3"]),
            ("reach dark x10", ["--relock-w", "3", "--dark", "reach",
                                "--dark-rate", "10"]),
            ("reach dark x30", ["--relock-w", "3", "--dark", "reach",
                                "--dark-rate", "30"]))
SEEDS = (["--seed", "11"], ["--seed", "12"], [])
COLS = ("static own", "parent = A'", "rate A", "rate own", "drift/U",
        "coherence", "age", "mu=0 ctl")
SNAP = re.compile(r"^\s*([\d.]+)\s+\d+\s+[-\d.]+\s+\S+\s+([-\d.]+)")


def late_means(out):
    """Mean over snapshots t >= 4 of the readers line, with the control."""
    rows, t, ctl = [], None, None
    for line in out.splitlines():
        m = SNAP.match(line)
        if m:
            t, ctl = float(m.group(1)), float(m.group(2))
        if "readers:" in line and t is not None and t >= 4.0:
            v = [float(z) for z in re.findall(r"-?\d+\.\d+", line)]
            rows.append(v[:6] + [v[7], ctl])
    return np.array(rows).mean(axis=0)


def run_row(env, sea, reg, args):
    a = np.array([late_means(subprocess.run(
        [sys.executable, "-u", DEMO] + BASE + args + sd,
        capture_output=True, text=True, env=env, check=True).stdout)
        for sd in SEEDS])
    m, se = a.mean(0), a.std(0, ddof=1) / np.sqrt(len(a))
    cells = [f"{m[k]:6.3f}+-{se[k]:.3f}" if k < 3 else f"{m[k]:12.4f}"
             if k == 4 else f"{m[k]:12.3f}" for k in range(len(COLS))]
    print(f"{sea:18s} {reg:17s} " + " ".join(f"{c:>12s}" for c in cells),
          flush=True)
    return m, se


def main():
    env = dict(os.environ, PYTHONPATH=HERE)
    head = f"{'sea':18s} {'regime':17s} " + " ".join(f"{c:>12s}" for c in COLS)
    print("Table 1: the sea of section 6 ((S), continuous)\n" + head, flush=True)
    for reg, r_args in SECTION6:
        run_row(env, "(S), continuous", reg, r_args)
    print("\nTable 2: sea variants\n" + head, flush=True)
    res = {}
    for sea, s_args in SEAS:
        for reg, r_args in REGIMES:
            res[sea, reg] = run_row(env, sea, reg, s_args + r_args)

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11, 4.0))
    xs = np.arange(len(REGIMES))
    colors = ("#D85A30", "#7F77DD", "#1D9E75")
    ctl = np.mean([res[k][0][7] for k in res])
    for (sea, _), col in zip(SEAS, colors):
        m = np.array([res[sea, r][0][2] for r, _ in REGIMES])
        e = np.array([res[sea, r][1][2] for r, _ in REGIMES])
        ax0.errorbar(xs, m, yerr=e, marker="o", color=col, capsize=3, label=sea)
        ax1.plot(xs, [res[sea, r][0][5] for r, _ in REGIMES], "o-", color=col,
                 label=sea)
    ax0.axhline(ctl, color="k", ls="--", lw=0.8, label=r"$\mu \equiv 0$ control")
    for ax, yl, tt in ((ax0, r"corr(rate reading A, $K_q$)",
                        "the parent-frame rate reading near the barrier"),
                       (ax1, r"sea coherence $|\langle e^{i\mu_{ref}}\rangle|$",
                        "the lock it reads against")):
        ax.set_xticks(xs, [r for r, _ in REGIMES])
        ax.set_ylabel(yl)
        ax.set_title(tt, fontsize=9)
        ax.set_ylim(0, 1)
        ax.legend(fontsize=7, loc="lower right")
    fig.tight_layout()
    name = "sea_lock_readers_scan.png"
    fig.savefig(output_path(name), dpi=130, bbox_inches="tight")
    dp_ = docs_path(name)
    if dp_:
        fig.savefig(dp_, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {name}")


if __name__ == "__main__":
    main()
