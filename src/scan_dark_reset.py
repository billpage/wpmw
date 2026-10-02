"""Dark catalysis as a reset clock -- step 23 section 8.

Runs ``demo_sea_lock_particles.py --readers --filter`` on the force-blind,
row-centred sea of postulate (S') (``--sea-force blind --sea-p rows``),
three seeds per configuration, and prints late means (4 <= t <= 25) with
standard errors.  Every run but the first keeps theta - p x / hbar
continuous across the periodic boundary (``--wrap-phase``).

  table 1   the box: streaming only, with and without --wrap-phase, and the
            published step 22 section 9 configuration (events, reach dark
            catalysis x10) with and without it
  table 2   reach dark catalysis (relative, gauge-invariant resets) at 1, 2,
            3, 6, 10 and 30 times the kernel's own rate sum_{q>=1} |K_q|,
            without and with events
  table 3   absolute resets to the free plane wave (--dark reset, an
            external clock) at the same rates, no events

Rates follow demo_sea_lock_particles.py's default --rate-convention qle.
The published scan (step 23 section 8, before the rate fix of the third
addendum) drew every event and every dark firing at twice the kernel's
rate, so its x1 and x3 are x2 and x6 here; without events those rows are
the same dynamics, and they reproduce.

Columns: the sea-only reading near the barrier (one member per aligned
pair) and its mu = 0 control; tau, the regression slope of the
misalignment on U/hbar; the rms residual misalignment; the fitted
quadrature coefficient b of (reading - control) on the cosine term, and
the share of that difference it explains; the all-bodies rate reading (A)
of step 22 section 9 and its mu = 0 control; the sea coherence; the median
clock age.

Writes sea_reset_filter.png.  About ten minutes.

    PYTHONPATH=src python3 -u src/scan_dark_reset.py
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
BASE = ["--nu", "8", "--mode", "locked", "--t-end", "25", "--readers",
        "--filter", "--sea-force", "blind", "--sea-p", "rows"]
W = ["--wrap-phase"]
EV = ["--relock-w", "3"]
NOEV = ["--no-events"]
SEEDS = (["--seed", "11"], ["--seed", "12"], [])
RATES = (1, 2, 3, 6, 10, 30)
SNAP = re.compile(r"^\s*([\d.]+)\s+\d+\s+[-\d.]+\s+\S+\s+([-\d.]+)")
FLT = re.compile(
    r"sea-only reading\s+(\S+)\s+control\s+(\S+)\s+tau\s+(\S+)\s+resid\s+(\S+)"
    r"\s+leak b\s+(\S+)\s+leak corr\s+(\S+)\s+leak share\s+(\S+)")
RDR = re.compile(r"rate \(A\)\s+(\S+).*sea coherence parent\s+(\S+).*"
                 r"median clock age\s+(\S+)")
COLS = ("reading", "control", "tau", "resid", "leak b", "leak share",
        "rate (A)", "all ctl", "coherence", "age")


def late_means(out):
    rows, t, cur, ctl = [], None, None, None
    for line in out.splitlines():
        m = SNAP.match(line)
        if m:
            t, ctl = float(m.group(1)), float(m.group(2))
        m = RDR.search(line)
        if m and t is not None and t >= 4.0:
            cur = [float(m.group(1)), ctl, float(m.group(2)),
                   float(m.group(3))]
        m = FLT.search(line)
        if m and t is not None and t >= 4.0 and cur is not None:
            f = [float(m.group(k)) for k in range(1, 8)]
            rows.append([f[0], f[1], f[2], f[3], f[4], f[6]] + cur)
            cur = None
    return np.array(rows).mean(axis=0)


def run_row(env, label, args):
    a = np.array([late_means(subprocess.run(
        [sys.executable, "-u", DEMO] + BASE + args + sd,
        capture_output=True, text=True, env=env, check=True).stdout)
        for sd in SEEDS])
    m, se = a.mean(0), a.std(0, ddof=1) / np.sqrt(len(a))
    cells = ([f"{m[0]:.3f}+-{se[0]:.3f}", f"{m[1]:.3f}"]
             + [f"{m[k]:+.4f}" for k in (2,)] + [f"{m[3]:.3f}"]
             + [f"{m[4]:+.4f}", f"{m[5]:.3f}"]
             + [f"{m[6]:.3f}+-{se[6]:.3f}", f"{m[7]:.3f}", f"{m[8]:.3f}",
                f"{m[9]:.2f}"])
    print(f"{label:34s} " + " ".join(f"{c:>13s}" for c in cells), flush=True)
    return m, se


def main():
    env = dict(os.environ, PYTHONPATH=HERE)
    head = f"{'configuration':34s} " + " ".join(f"{c:>13s}" for c in COLS)
    res = {}
    print("Table 1: the box\n" + head, flush=True)
    res["s0"] = run_row(env, "streaming only, no wrap phase", NOEV)
    res["s1"] = run_row(env, "streaming only, wrap phase", NOEV + W)
    res["p0"] = run_row(env, "events + reach x10, no wrap phase",
                        EV + ["--dark", "reach", "--dark-rate", "10"])
    res["p1"] = run_row(env, "events + reach x10, wrap phase",
                        EV + W + ["--dark", "reach", "--dark-rate", "10"])
    print("\nTable 2: reach dark catalysis (relative resets), wrap phase\n"
          + head, flush=True)
    for r in RATES:
        res["rn", r] = run_row(env, f"no events, reach x{r}",
                               NOEV + W + ["--dark", "reach", "--dark-rate", str(r)])
    res["e0"] = run_row(env, "events, no dark catalysis", EV + W)
    for r in RATES:
        res["re", r] = run_row(env, f"events, reach x{r}",
                               EV + W + ["--dark", "reach", "--dark-rate", str(r)])
    print("\nTable 3: absolute resets (an external clock), wrap phase\n"
          + head, flush=True)
    for r in RATES:
        res["an", r] = run_row(env, f"no events, reset x{r}",
                               NOEV + W + ["--dark", "reset", "--dark-rate", str(r)])

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11, 4.0))
    rs = np.array(RATES, float)
    series = (("rn", "reach dark catalysis, no events", "#1D9E75", "o-"),
              ("re", "reach dark catalysis, with events", "#7F77DD", "s-"),
              ("an", "absolute reset, no events", "#D85A30", "^--"))
    ctl_n = np.mean([res[k, r][0][1] for k in ("rn", "an") for r in RATES])
    ctl_e = np.mean([res["re", r][0][1] for r in RATES])
    for k, lab, col, sty in series:
        m = np.array([res[k, r][0][0] for r in RATES])
        e = np.array([res[k, r][1][0] for r in RATES])
        ax0.errorbar(rs, m, yerr=e, fmt=sty, color=col, capsize=3, label=lab)
        ax1.loglog(rs, [res[k, r][0][2] for r in RATES], sty, color=col,
                   label=lab)
    ax0.axhline(res["s1"][0][0], color="0.5", ls=":", lw=1,
                label="streaming only (the eikonal sea)")
    ax0.axhline(ctl_n, color="k", ls="--", lw=0.8,
                label=r"$\mu \equiv 0$ control, no events")
    ax0.axhline(ctl_e, color="#7F77DD", ls="--", lw=0.8,
                label=r"$\mu \equiv 0$ control, with events")
    ax0.set_xscale("log")
    ax0.set_xticks(rs, [f"x{r}" for r in RATES])
    ax0.set_xlabel("dark rate / kernel rate")
    ax0.set_ylabel(r"corr(sea-only reading, $K_q$)")
    ax0.set_ylim(0.5, 1.0)
    ax0.set_title("the reading near the barrier, force-blind sea", fontsize=9)
    ax0.legend(fontsize=7, loc="lower right")
    g = res["rn", 1][0][2]
    ax1.loglog(rs, g / rs, "k:", lw=0.8, label=r"$\propto 1/\kappa$")
    ax1.set_xticks(rs, [f"x{r}" for r in RATES])
    ax1.set_xlabel("dark rate / kernel rate")
    ax1.set_ylabel(r"$\tau$: slope of $\mu$ on $U/\hbar$")
    ax1.set_title("the effective reset time", fontsize=9)
    ax1.legend(fontsize=7)
    fig.tight_layout()
    name = "sea_reset_filter.png"
    fig.savefig(output_path(name), dpi=130, bbox_inches="tight")
    dp_ = docs_path(name)
    if dp_:
        fig.savefig(dp_, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {name}")


if __name__ == "__main__":
    main()
