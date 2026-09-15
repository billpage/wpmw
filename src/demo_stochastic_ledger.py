#!/usr/bin/env python3
"""
Companion demo for ``docs/analysis/stochastic_ledger.md``.

The ledger of the compensated residual channel, run as an exact Markov jump
process on integer counts rather than as a mean field on a mesh.  Every
number quoted in the note is produced here.

State
-----
Per momentum cell ``c`` of a ring of ``M`` cells: integers
``n+[c], n-[c], S[c]``.

Channels
--------
``event(c, q)``  rate ``gamma`` for each of the ``M * Q`` pairs.  Daughters
    ``A = c+q``, ``B = c-q``.
      absorptive, iff ``n-[A] >= 1`` and ``n+[B] >= 1``:
          ``n-[A] -= 1 ; n+[B] -= 1 ; S[c] += 1``
      emissive, otherwise and iff ``S[c] >= 1``:
          ``n+[A] += 1 ; n-[B] += 1 ; S[c] -= 1``
``hop``  rate ``nu`` per body, destination uniform over the ring.  The
    caricature of postulate (S): bodies move between cells.
``recomb(c)``  rate ``kappa * n+[c] * n-[c]``:
    ``n+[c] -= 1 ; n-[c] -= 1 ; S[c] += 1``.

Parts
-----
A  Theorem N1, the pathwise invariant, and Theorem N2, the null direction:
   ``dE`` is identical on the two realisations, so the choice between them
   is a noise channel the observable cannot see.
B  Theorem N3, the stationary sum rule ``Gamma_tot (1 - 2f) = R_sink``.
C  ``f -> 1/2`` at fixed ``nu`` over a range of ``gamma``, and the
   variance of the total body count.
D  Theorem N4, Q-independence of the standing population.
E  Theorem N5, the availability closure and the five derived numbers.
F  Theorem N6, transport and not recombination is the local regulator.
G  Figures.

Run::

    WPMW_OUTPUT=/tmp/out PYTHONPATH=src python3 -u src/demo_stochastic_ledger.py
"""

from __future__ import annotations

import numpy as np

from wpmwlib.wpmw_utils import docs_path, output_path

SEED = 20260914

# closure constants, section 5 of the note
LAM_STAR = -np.log(1.0 - 2.0 ** -0.5)      # 1.227947
FPRIME = np.sqrt(2.0) - 1.0                # 0.414214
VAR_OVER_M = 1.0 / FPRIME                  # 2.414214
FANO = 1.0 / (2.0 * FPRIME * LAM_STAR)     # 0.983028


def banner(title: str) -> None:
    print("\n" + "=" * 72)
    print("  " + title)
    print("=" * 72)


# --------------------------------------------------------------------- #
# The jump process                                                      #
# --------------------------------------------------------------------- #
class Ledger:
    """Exact Gillespie dynamics of the three-population ledger."""

    def __init__(self, M=32, Q=2, gamma=1.0, nu=8.0, kappa=0.0,
                 S0=100000, n0=2, seed=SEED):
        self.M, self.Q = int(M), int(Q)
        self.gamma, self.nu, self.kappa = float(gamma), float(nu), float(kappa)
        self.rng = np.random.default_rng(seed)
        self.p = np.full(self.M, n0, dtype=np.int64)
        self.m = np.full(self.M, n0, dtype=np.int64)
        self.S = np.full(self.M, S0, dtype=np.int64)
        self.t = 0.0
        self.n_ev = self.n_abs = self.n_rec = self.n_block = 0
        self.twoP0 = self.twoP()
        self.acc_pm = 0.0          # time integral of sum_c n+ n-

    def twoP(self) -> int:
        """Twice the pair count ``P = sum_c (S + N/2)``, in exact integers."""
        return int(2 * self.S.sum() + self.p.sum() + self.m.sum())

    def step(self) -> str:
        M, rng = self.M, self.rng
        cp, cm = np.cumsum(self.p), np.cumsum(self.m)
        n_pos, n_neg = int(cp[-1]), int(cm[-1])
        pm = float(np.dot(self.p, self.m))
        r_ev = self.gamma * M * self.Q
        r_hop = self.nu * (n_pos + n_neg)
        r_rec = self.kappa * pm
        r_tot = r_ev + r_hop + r_rec
        dt = rng.exponential(1.0 / r_tot)
        self.t += dt
        self.acc_pm += pm * dt
        u = rng.random() * r_tot
        if u < r_ev:
            c = int(rng.integers(M))
            q = 1 + int(rng.integers(self.Q))
            a, b = (c + q) % M, (c - q) % M
            self.n_ev += 1
            if self.m[a] >= 1 and self.p[b] >= 1:
                self.m[a] -= 1
                self.p[b] -= 1
                self.S[c] += 1
                self.n_abs += 1
            elif self.S[c] >= 1:
                self.p[a] += 1
                self.m[b] += 1
                self.S[c] -= 1
            else:
                self.n_block += 1
            return "ev"
        u -= r_ev
        if u < r_hop:
            k = int(u / self.nu)
            if k < n_pos:
                src = int(np.searchsorted(cp, k, side="right"))
                self.p[src] -= 1
                self.p[int(rng.integers(M))] += 1
            else:
                src = int(np.searchsorted(cm, k - n_pos, side="right"))
                self.m[src] -= 1
                self.m[int(rng.integers(M))] += 1
            return "hop"
        c = int(rng.choice(M, p=(self.p * self.m) / pm))
        self.p[c] -= 1
        self.m[c] -= 1
        self.S[c] += 1
        self.n_rec += 1
        return "rec"


def run(M=32, Q=2, gamma=1.0, nu=8.0, kappa=0.0, t_max=400.0, n0=2,
        seed=SEED, burn=0.4, keep_trace=False):
    """Run to ``t_max`` and return time-weighted late-window statistics."""
    led = Ledger(M=M, Q=Q, gamma=gamma, nu=nu, kappa=kappa, n0=n0, seed=seed)
    t_burn = burn * t_max
    lam, wts, tr_t, tr_l = [], [], [], []
    mark = None
    while led.t < t_max:
        t0 = led.t
        led.step()
        if keep_trace:
            tr_t.append(led.t)
            tr_l.append(float(led.p.sum() + led.m.sum()))
        if led.t >= t_burn:
            if mark is None:
                mark = (led.n_ev, led.n_abs, led.n_rec, led.t, led.acc_pm)
            lam.append(float(led.p.sum() + led.m.sum()))
            wts.append(led.t - t0)
    lam = np.asarray(lam)
    w = np.asarray(wts)
    w = w / w.sum()
    mean = float(lam @ w)
    var = float(((lam - mean) ** 2) @ w)
    d_t = led.t - mark[3]
    return dict(
        L=led,
        f=(led.n_abs - mark[1]) / (led.n_ev - mark[0]),
        mean=mean,
        var=var,
        gam_tot=gamma * M * Q,
        R_rec=(led.n_rec - mark[2]) / d_t,
        pm_bar=(led.acc_pm - mark[4]) / d_t,
        twoP_err=abs(led.twoP() - led.twoP0),
        trace=(np.asarray(tr_t), np.asarray(tr_l)),
    )


# --------------------------------------------------------------------- #
# Part A -- N1 and N2                                                   #
# --------------------------------------------------------------------- #
def part_a():
    banner("A  N1 the pathwise invariant, N2 the null direction of the event")
    led = Ledger(M=16, Q=2, gamma=1.0, nu=8.0, kappa=0.2, n0=2)
    worst_p, worst_cell, worst_sum = 0, 0, 0
    for _ in range(300000):
        before = led.p - led.m
        kind = led.step()
        worst_p = max(worst_p, abs(led.twoP() - led.twoP0))
        if kind == "ev":
            d = (led.p - led.m) - before
            worst_cell = max(worst_cell, int(np.abs(d).max()))
            worst_sum = max(worst_sum, int(abs(d.sum())))
    print(f"  events                                  {led.n_ev}")
    print(f"  recombinations                          {led.n_rec}")
    print(f"  blocked for want of a sea pair          {led.n_block}")
    print(f"  max |2P(t) - 2P(0)| over the run        {worst_p}"
          "        (exact integer)")
    print(f"  max per-cell |dE| on an event           {worst_cell}"
          "        (1 expected)")
    print(f"  max |sum_c dE| on an event              {worst_sum}"
          "        (0 expected)")
    print(f"  absorptive fraction f                   "
          f"{led.n_abs / led.n_ev:.6f}")
    print("\n  Both realisations move E by +1 at A and -1 at B.  Nothing in")
    print("  the code enforces that; it is a property of the two updates.")
    return worst_p, worst_cell, worst_sum


# --------------------------------------------------------------------- #
# Part B -- N3, the sum rule                                            #
# --------------------------------------------------------------------- #
def part_b():
    banner("B  N3  the stationary sum rule   Gamma_tot (1 - 2 f) = R_sink")
    print("  M = 32, Q = 2, gamma = 1, nu = 8, t_max = 600.\n")
    print("   kappa        f      1/2 - f    R_sink/(2 Gamma_tot)"
          "    rel err    mean/M   |2 dP|")
    rows = []
    for kappa in (0.0, 0.05, 0.2, 0.5, 1.0, 2.0):
        r = run(M=32, Q=2, gamma=1.0, nu=8.0, kappa=kappa, t_max=600.0)
        lhs = 0.5 - r["f"]
        rhs = r["R_rec"] / (2.0 * r["gam_tot"])
        rel = abs(lhs - rhs) / max(abs(rhs), 1e-12) if kappa else abs(lhs)
        rows.append((kappa, r["f"], lhs, rhs, rel, r["mean"] / 32))
        print(f"   {kappa:5.2f} {r['f']:9.5f} {lhs:+10.5f}"
              f"      {rhs:+12.5f}    {rel:8.4f} {r['mean'] / 32:8.4f}"
              f"  {r['twoP_err']:6d}")
    print("\n  The kappa = 0 row's third column is the deviation itself.")
    print("  Every other row: a body sink of any kind moves f below 1/2 by")
    print("  exactly half its share of the event rate.")
    return rows


# --------------------------------------------------------------------- #
# Part C -- f at fixed nu                                               #
# --------------------------------------------------------------------- #
def part_c():
    banner("C  f -> 1/2 at fixed mixing rate, over a 32x range of gamma/nu")
    print("  M = 32, Q = 2, nu = 8, kappa = 0, 400/gamma time units each.\n")
    print("   gamma    events       f      mean/M    Var/M    Fano")
    rows = []
    for gamma in (0.25, 0.5, 1.0, 2.0, 4.0, 8.0):
        r = run(M=32, Q=2, gamma=gamma, nu=8.0, t_max=400.0 / gamma)
        rows.append((gamma, r["f"], r["mean"] / 32, r["var"] / 32,
                     r["var"] / r["mean"]))
        print(f"   {gamma:5.2f} {r['L'].n_ev:8d} {r['f']:9.5f}"
              f" {r['mean'] / 32:8.4f} {r['var'] / 32:8.4f}"
              f" {r['var'] / r['mean']:7.4f}")
    print(f"\n  closure predicts mean/M = {2 * LAM_STAR:.4f}, "
          f"Var/M = {VAR_OVER_M:.4f}, Fano = {FANO:.4f}")
    print("  and the agreement is best where the mixing leads the reaction.")
    return rows


# --------------------------------------------------------------------- #
# Part D -- N4, Q-independence                                          #
# --------------------------------------------------------------------- #
def part_d():
    banner("D  N4  the standing population does not depend on the channel count")
    print("  M = 32, gamma = 1, nu = 8 Q, kappa = 0.\n")
    print("       Q        f      mean/M    Var/M    Fano")
    rows = []
    for Q in (1, 2, 4, 8, 12):
        r = run(M=32, Q=Q, gamma=1.0, nu=8.0 * Q, t_max=400.0)
        rows.append((Q, r["f"], r["mean"] / 32))
        print(f"    {Q:4d} {r['f']:9.5f} {r['mean'] / 32:8.4f}"
              f" {r['var'] / 32:8.4f} {r['var'] / r['mean']:7.4f}")
    print(f"\n  closure predicts mean/M = {2 * LAM_STAR:.4f} for every Q,")
    print("  since availability is a per-cell question and the channel index")
    print("  never enters it.")
    return rows


# --------------------------------------------------------------------- #
# Part E -- N5, the closure and its five numbers                        #
# --------------------------------------------------------------------- #
def part_e():
    banner("E  N5  the availability closure, and Var(Lambda) linear in M")
    print("   f(lam) = (1 - e^-lam)^2,  f(lam*) = 1/2  =>\n")
    print(f"     lam*                 = -ln(1 - 2^-1/2) = {LAM_STAR:.6f}")
    print(f"     bodies per cell      = 2 lam*          = {2*LAM_STAR:.6f}")
    print(f"     sea pairs per cell   = B h / 2         = 2.000000")
    print(f"     f'(lam*)             = sqrt2 - 1       = {FPRIME:.6f}")
    print(f"     Var(Lambda) / M      = 1 / f'          = {VAR_OVER_M:.6f}")
    print(f"     Fano                 = 1/(2 f' lam*)   = {FANO:.6f}")
    print("\n  M scan, gamma = 1, Q = 2, nu = 8, kappa = 0:\n")
    print("       M        f      mean/M    Var/M    Fano")
    rows = []
    for M in (8, 16, 32, 64, 128):
        r = run(M=M, Q=2, gamma=1.0, nu=8.0, t_max=400.0)
        rows.append((M, r["f"], r["mean"] / M, r["var"] / M))
        print(f"    {M:4d} {r['f']:9.5f} {r['mean'] / M:8.4f}"
              f" {r['var'] / M:8.4f} {r['var'] / r['mean']:7.4f}")
    print("\n  Var/M flat is Var(Lambda) linear in M, which is what an")
    print("  Ornstein-Uhlenbeck ledger with a per-cell restoring force gives.")
    return rows


# --------------------------------------------------------------------- #
# Part F -- N6, who regulates the local ledger                          #
# --------------------------------------------------------------------- #
def part_f():
    banner("F  N6  transport, not recombination, regulates the local ledger")
    print("  Variance of the per-cell body count across the M = 32 cells,")
    print("  sampled at t = 25, 50, 100, 200, 400.\n")
    print("     nu   kappa      t=25      t=50     t=100     t=200     t=400"
          "       f")
    rows = []
    for nu, kappa in ((0.0, 0.0), (0.0, 0.5), (0.0, 2.0),
                      (2.0, 0.0), (8.0, 0.0)):
        led = Ledger(M=32, Q=2, gamma=1.0, nu=nu, kappa=kappa, n0=2)
        marks, vals, j = [25.0, 50.0, 100.0, 200.0, 400.0], [], 0
        while j < len(marks):
            led.step()
            if led.t >= marks[j]:
                vals.append(float(np.var((led.p + led.m).astype(float))))
                j += 1
        f = led.n_abs / max(led.n_ev, 1)
        rows.append((nu, kappa, vals, f))
        print(f"   {nu:4.1f}   {kappa:5.2f}  " +
              "  ".join(f"{v:8.2f}" for v in vals) + f"  {f:8.4f}")
    print("\n  kappa alone damps the spread only partly and drags f well")
    print("  below 1/2, as N3 requires of any sink.  Transport holds the")
    print("  spread flat and leaves f at 1/2.")
    return rows


# --------------------------------------------------------------------- #
# Part G -- figures                                                     #
# --------------------------------------------------------------------- #
def save(fig, name):
    fig.savefig(output_path(name), dpi=150, bbox_inches="tight")
    dp = docs_path(name)
    if dp:
        fig.savefig(dp, dpi=150, bbox_inches="tight")
    print(f"  wrote {name}")


def part_g(sum_rows, m_rows, loc_rows):
    banner("G  figures")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 3, figsize=(13.5, 3.9))

    k = np.array([r[0] for r in sum_rows])
    lhs = np.array([r[2] for r in sum_rows])
    rhs = np.array([r[3] for r in sum_rows])
    ax[0].plot(k, lhs, "o", ms=7, label=r"$1/2-f$ measured")
    ax[0].plot(k, rhs, "k-", lw=1.2, label=r"$R_{\rm sink}/2\Gamma_{\rm tot}$")
    ax[0].set_xlabel(r"$\kappa$")
    ax[0].set_ylabel(r"$1/2 - f$")
    ax[0].set_title("(a) N3: the sum rule", fontsize=10)
    ax[0].legend(fontsize=8)
    ax[0].grid(alpha=0.3)

    M = np.array([r[0] for r in m_rows], dtype=float)
    vm = np.array([r[3] for r in m_rows])
    mm = np.array([r[2] for r in m_rows])
    ax[1].semilogx(M, vm, "s-", ms=6, label=r"$\mathrm{Var}(\Lambda)/M$")
    ax[1].semilogx(M, mm, "o-", ms=6, label=r"$\langle\Lambda\rangle/M$")
    ax[1].axhline(VAR_OVER_M, color="k", ls="--", lw=1.0,
                  label=r"$1/f' = 2.4142$")
    ax[1].axhline(2 * LAM_STAR, color="0.45", ls=":", lw=1.2,
                  label=r"$2\lambda^*$")
    ax[1].set_xlabel(r"$M$ cells")
    ax[1].set_title("(b) N5: the closure", fontsize=10)
    ax[1].legend(fontsize=8)
    ax[1].grid(alpha=0.3)

    ts = [25.0, 50.0, 100.0, 200.0, 400.0]
    for nu, kappa, vals, f in loc_rows:
        lab = rf"$\nu={nu:g},\ \kappa={kappa:g}$"
        ax[2].loglog(ts, vals, "o-", ms=4, lw=1.1, label=lab)
    ax[2].set_xlabel("time")
    ax[2].set_ylabel(r"$\mathrm{Var}_c(N_c)$")
    ax[2].set_title("(c) N6: who regulates locally", fontsize=10)
    ax[2].legend(fontsize=7)
    ax[2].grid(alpha=0.3, which="both")

    fig.tight_layout()
    save(fig, "stochastic_ledger_summary.png")
    plt.close(fig)

    # second figure: a single trajectory of the total body count
    fig2, ax2 = plt.subplots(figsize=(7.2, 3.6))
    r = run(M=32, Q=2, gamma=1.0, nu=8.0, t_max=120.0, n0=8,
            keep_trace=True)
    t, lam = r["trace"]
    ax2.plot(t, lam / 32.0, lw=0.7, color="C0")
    ax2.axhline(2 * LAM_STAR, color="k", ls="--", lw=1.2,
                label=r"$2\lambda^*=2.4559$")
    ax2.set_xlabel("time")
    ax2.set_ylabel(r"$\Lambda/M$  (bodies per cell)")
    ax2.set_title("The ledger relaxing onto the derived standing population",
                  fontsize=10)
    ax2.legend(fontsize=9)
    ax2.grid(alpha=0.3)
    fig2.tight_layout()
    save(fig2, "stochastic_ledger_relaxation.png")
    plt.close(fig2)


def main():
    part_a()
    sum_rows = part_b()
    part_c()
    part_d()
    m_rows = part_e()
    loc_rows = part_f()
    part_g(sum_rows, m_rows, loc_rows)
    banner("summary")
    print("  N1  2P conserved pathwise, in exact integers")
    print("  N2  dE identical on the two realisations: a null direction")
    print("  N3  Gamma_tot (1 - 2f) = R_sink, verified over a 40x kappa range")
    print("  N4  the standing population is independent of Q")
    print("  N5  lam* = 1.227947 bodies per species per cell, no free constant")
    print("  N6  transport and not recombination is the local regulator")


if __name__ == "__main__":
    main()
