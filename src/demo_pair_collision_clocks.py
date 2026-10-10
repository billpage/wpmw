"""Pair collisions with clocks -- whole pairs, partner exchange, phase matching.

Companion to step 24's collision addendum (docs/analysis/sea_depletion.md,
Definition (C)) read against step 23's phase picture (force_blind_sea.md,
Propositions Q2, Q5, Q11, Q12; open items Q-SP10, Q-SP11, V-SP8).

Definition (C) moves pair DENSITY between rows: two pairs on row n go to
n - q and n + q (defocus), or the reverse (focus).  Bodies carry clocks, so a
collision must also say what happens to them.  Three body-level realisations:

  whole      each pair kinks as a unit (both members +-xi, clocks continuous,
             Theorem Y4): a pair's own misalignment mu is untouched, but it
             lands off its new row's locked field.
  exchange   the partners swap: (C+, D-) go one way, (D+, C-) the other, and
             each new pair forms at the midpoint of its members (clocks
             transported, as for recombination in the step 23 model).
  matched    exchange, accepted only to the extent that both products would
             be dark: weight cos^2(mu1/2) cos^2(mu2/2), and for a defocus the
             orientation (which positon goes up) chosen by the same weight.

Parts
  A  exact (SymPy): the phase algebra of the realisations.  Whole pairs keep
     mu exactly.  Exchange products carry the chord's reading against their
     daughter rows (Proposition Q5(a)), and mu1 + mu2 = mu_C + mu_D
     + o 2 xi d / hbar, o the orientation: total misalignment is conserved
     by co-located exchange and shifted by the chord's Fourier phase
     otherwise.  Also the body bookkeeping of a two-pair defocus realised as
     dark catalysis (Theorem L5) with the second member ionising its own
     pair: dE = dN = dS = 0, momentum conserved.
  B  mesh: (1) the chord-level drive.  A trigger that fires each chord at
     the magnitude of its own contribution to K_q (Proposition Q5) has gross
     rate (B dp / hbar) int |U_res w sin(2 xi y / hbar)| dy, against the
     kernel's |K_q|; the kernel's falloff in q is an interference between
     chords.  (2) Mean-field stability of collisions realised as catalysed
     W-null events (a catalyst pair firing at the kernel rate, realised by
     the allocation rule) against Definition (C).
  C  particle model (step 23, section 8 and Q12): Eckart barrier, (S'),
     force-blind row-centred locked sea, wrap phase, events Poisson at
     |K_q| realised as in demo_sea_resonance_clock.py, reach dark catalysis
     at kappa times the kernel's own rate, and sea collisions in one of the
     three realisations at the kernel rate (Definition (C), r_q B = |K_q|/2
     per pair at the physical density; the ensemble multiplicity nu divides
     the per-encounter rate).  Measured: grayness of the sea's pairs, the
     coherence of the readers' chords (the lock), the live reading in the
     shadow of mesh-triggered events (Q12(e)), the exchange identity event by
     event, T_E, event failures, and the sea's row profile near the barrier.

Units hbar = m = 1; Eckart V0 = a = 1; dp = 0.25, n_p = 64, y_max = 2 pi;
box L = 48 with periodic wrap.

Usage::

    PYTHONPATH=src python3 -u src/demo_pair_collision_clocks.py --parts AB
    PYTHONPATH=src python3 -u src/demo_pair_collision_clocks.py --parts C \\
        --nu 32 --coll none whole exchange matched --kappa 0 2 --seeds 4

Writes pair_collision_clocks_C<tag>.csv (one row per run) and
pair_collision_clocks_trace<tag>.csv through wpmwlib.wpmw_utils.
"""
import argparse
import csv
import multiprocessing as mp
import os
import time

import numpy as np

from demo_emission_and_absorption import Ledger, B, MU, HBAR
from wpmwlib.wpmw_utils import output_path

ap = argparse.ArgumentParser(description=__doc__,
                             formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--parts", default="AB")
ap.add_argument("--nu", type=int, nargs="+", default=[16])
ap.add_argument("--coll", nargs="+",
                default=["none", "whole", "exchange", "matched", "relock"])
ap.add_argument("--kappa", type=float, nargs="+", default=[2.0],
                help="dark-catalysis rate in units of sum_q |K_q| (0: off)")
ap.add_argument("--coll-scale", type=float, nargs="+", default=[1.0],
                help="collision rate in units of the kernel rate of Definition (C)")
ap.add_argument("--coll-window", type=float, default=1.0,
                help="two pairs collide within |dx| < this (a cell, as Definition (C))")
ap.add_argument("--shadow", action="store_true",
                help="integrate the live reading along free bodies (Q12(e)); costly")
ap.add_argument("--seeds", type=int, default=2, help="seeds 11, 12, ...")
ap.add_argument("--t-end", type=float, default=14.0)
ap.add_argument("--workers", type=int, default=0)
ap.add_argument("--tag", default="")
ARGS = ap.parse_args()

run = Ledger(v0=1.0, a=1.0, n_r=192, r_half=24.0, n_p=64, dp=0.25)
L, DP, NP = 48.0, run.dp, run.n_p
PMAX, YMAX = NP * DP / 2, run.y_max
Q = np.arange(1, NP // 2)
XI = Q * DP
NQ = len(Q)
KQ = run.k[:, Q]
GAM = np.abs(KQ).sum(axis=1)
KCUMQ = np.cumsum(np.abs(KQ), axis=1)
XACT = 8.0
DT = 0.02


def V(x):
    return 1.0 / np.cosh(x) ** 2


def F(x):
    return 2.0 * np.tanh(x) / np.cosh(x) ** 2


def interp_rows(vals, xs):
    idx = np.clip(np.searchsorted(run.r, xs) - 1, 0, len(run.r) - 2)
    w = ((xs - run.r[idx]) / run.dr)[:, None]
    return vals[idx] * (1 - w) + vals[idx + 1] * w


def rule(t):
    print("\n" + "=" * 78 + f"\n{t}\n" + "=" * 78, flush=True)


def wrap(a):
    return (a + np.pi) % (2 * np.pi) - np.pi


# ------------------------------------------------------------------ Part A
def part_a():
    import sympy as sp
    rule("A. Exact: the phase algebra of a collision")
    p, xi, h, xC, xD, xs = sp.symbols("p xi hbar x_C x_D x_s", real=True)
    ta, tb, tc, td = sp.symbols("theta_a theta_b theta_c theta_d", real=True)

    def Phi(th, P, x0):                  # P2: clock carried at wavevector P/hbar
        return th + P * (xs - x0) / h
    muC = Phi(ta, p, xC) - Phi(tb, p, xC)
    muD = Phi(tc, p, xD) - Phi(td, p, xD)
    d = xD - xC
    print("  Pair C = (C+, C-) at x_C, pair D = (D+, D-) at x_D, both on row p;"
          " clocks theta_a..theta_d.")
    w1 = sp.simplify(Phi(ta, p + xi, xC) - Phi(tb, p + xi, xC) - muC)
    w2 = sp.simplify(Phi(tc, p - xi, xD) - Phi(td, p - xi, xD) - muD)
    print(f"  whole pairs (C up, D down):  mu_C' - mu_C = {w1},  mu_D' - mu_D = {w2}")
    for o in (+1, -1):
        m1 = Phi(ta, p + o * xi, xC) - Phi(td, p + o * xi, xD)     # (C+, D-)
        m2 = Phi(tc, p - o * xi, xD) - Phi(tb, p - o * xi, xC)     # (D+, C-)
        r1 = sp.simplify(m1 - (ta - td + (p + o * xi) * d / h))
        r2 = sp.simplify(m2 - (tc - tb - (p - o * xi) * d / h))
        tot = sp.simplify(m1 + m2 - muC - muD)
        print(f"  exchange, (C+,D-) to p {'+' if o > 0 else '-'} xi:"
              f"  mu1 - chord^(p{'+' if o > 0 else '-'}xi) = {r1},"
              f"  mu2 - chord^(p{'-' if o > 0 else '+'}xi) = {r2},"
              f"  (mu1 + mu2) - (mu_C + mu_D) = {tot}")
    print("  -> whole pairs never make gray pairs; exchange conserves the total"
          " misalignment when\n     the pairs are co-located (d = 0) and shifts it"
          " by the chord's Fourier phase,\n     with a sign set by the orientation,"
          " otherwise.")

    print("\n  A two-pair defocus as dark catalysis (Theorem L5) with the second"
          " member\n  ionising its own pair.  Bodies on rows (p - xi, p, p + xi), as"
          " (positons, negatons):")
    # pair A (catalyst) and pair C on row p.  A+ fires, orientation t = +1,
    # realised by ionising C: C+ -> p + xi, C- -> p - xi.  A- fires,
    # orientation -1, realised by ionising A: A- -> p + xi, A+ -> p - xi.
    rows = {"before": {"p": [2, 2], "p+xi": [0, 0], "p-xi": [0, 0]},
            "after": {"p": [0, 0], "p+xi": [1, 1], "p-xi": [1, 1]}}
    for k, v in rows.items():
        print(f"    {k:6}  " + "  ".join(f"{r}: {tuple(c)}" for r, c in v.items()))
    dE = {r: (rows['after'][r][0] - rows['after'][r][1])
          - (rows['before'][r][0] - rows['before'][r][1]) for r in rows['after']}
    mom = sp.simplify(2 * (p + xi) + 2 * (p - xi) - 4 * p)
    print(f"    dE per row {dE};  dN = 0, dS = 0 (two pairs before, two after);"
          f"  d(momentum) = {mom}")
    print("    Event 1 (parent A+, +1 at p+xi, -1 at p-xi) and event 2 (parent A-,"
          " -1 at p+xi, +1 at p-xi)\n    cancel row by row, as in L5: the defocus"
          " is a W-null event.")


# ------------------------------------------------------------------ Part B
def part_b():
    rule("B1. The chord-level drive against the kernel's rate")
    y = np.linspace(-YMAX, YMAX, 4001)
    dy = y[1] - y[0]
    w = np.cos(np.pi * y / (2 * YMAX)) ** 2
    pos = y > 0
    print("  Per chord (half-separation y > 0) the contribution to K_q is"
          " -(B dp/hbar) U_res w sin(2 xi_q y/hbar),\n  each chord counted once."
          "  'gross' fires every chord at the magnitude of its own contribution.")
    rows = []
    for x0 in (-0.625, -1.5, -2.5):
        dv = np.interp(x0, run.r, run.dv_eff)
        ures = V(x0 + y) - V(x0 - y) - 2 * y * dv
        f = -(B * DP / HBAR) * (ures * w)[None, pos] * np.sin(
            2 * XI[:, None] * y[None, pos] / HBAR)
        K = 2 * f.sum(1) * dy
        G = 2 * np.abs(f).sum(1) * dy
        Pp = 2 * np.clip(f, 0, None).sum(1) * dy / np.maximum(G, 1e-300)
        Km = interp_rows(KQ, np.array([x0]))[0]
        print(f"  x = {x0:6.3f}: K_1 mesh {Km[0]:+.4f}, quadrature {K[0]:+.4f}")
        print("     q            " + " ".join(f"{q:>7d}" for q in Q[:6]))
        print("     K_q          " + " ".join(f"{v:+7.4f}" for v in K[:6]))
        print("     gross        " + " ".join(f"{v:7.4f}" for v in G[:6]))
        print("     gross/|K_q|  " + " ".join(f"{v:7.2f}" for v in G[:6] / np.abs(K[:6])))
        print("     share +      " + " ".join(f"{v:7.3f}" for v in Pp[:6]))
        print(f"     all channels: sum |K_q| = {np.abs(K).sum():.4f},"
              f"  sum gross = {G.sum():.4f}  (ratio {G.sum() / np.abs(K).sum():.1f})")
        rows.append((x0, np.abs(K).sum(), G.sum()))
    print("  Harmonic V = x^2/2: U_res = 0 identically, so both vanish (V1).")

    rule("B2. Mean field: collisions realised as catalysed W-null events")
    print("  One row lattice (64 rows, channel q = 1, rate g = 1, densities s in"
          " pairs per aperture).\n  A catalyst pair at row n fires; its event is"
          " realised by the allocation rule:\n    focus   a = (1 - e^-s(n+q))(1 - e^-s(n-q))"
          "   (a donor pair on each daughter row)\n    defocus b = 1 - e^-s(n)"
          "                    (a second pair on the catalyst's row)\n  J_n is the net"
          " defocus flux about row n; ds/dt = -2 J_n + J_(n-q) + J_(n+q).")
    n, q, g = 64, 1, 1.0

    def j_cat(s, alloc):
        sp_, sm = np.roll(s, -q), np.roll(s, q)
        a = (1 - np.exp(-sp_)) * (1 - np.exp(-sm))
        b = 1 - np.exp(-s)
        if alloc == "abs":
            pf, pd = a, (1 - a) * b
        else:
            pd, pf = b, (1 - b) * a
        return g * s * (pd - pf)

    def j_c(s):
        return g * (s * s - np.roll(s, q) * np.roll(s, -q)) / 2

    def rhs(s, J):
        jj = J(s)
        return -2 * jj + np.roll(jj, q) + np.roll(jj, -q)
    laws = (("Definition (C)", j_c),
            ("catalysed, absorptive first", lambda s: j_cat(s, "abs")),
            ("catalysed, emissive first", lambda s: j_cat(s, "emi")))
    print(f"\n  {'law':>30} {'s':>4} {'J(uniform)':>11} {'max growth':>11}"
          f" {'30% dip, t = 40: min':>21} {'max':>6}")
    for name, J in laws:
        for s0 in (1.0, 4.0):
            s = np.full(n, s0)
            f0 = rhs(s, J)
            jac = np.zeros((n, n))
            for k in range(n):
                ds = np.zeros(n)
                ds[k] = 1e-6
                jac[:, k] = (rhs(s + ds, J) - f0) / 1e-6
            ev = np.linalg.eigvals(jac).real.max()
            s = np.full(n, s0)
            s[n // 2] *= 0.7
            for _ in range(4000):
                k1 = rhs(s, J)
                k2 = rhs(s + 0.01 * k1, J)
                s = np.maximum(s + 0.005 * (k1 + k2), 0)
            print(f"  {name:>30} {s0:4.0f} {J(np.full(n, s0))[0]:+11.4f} {ev:+11.4f}"
                  f" {s.min() / s0:21.3f} {s.max() / s0:6.3f}", flush=True)
    print("  Definition (C) is in detailed balance (J = 0 on any uniform row) and"
          " relaxes.  A catalysed\n  event is three-body (catalyst plus two), so no"
          " allocation gives J = 0 on a uniform row:\n  absorptive first clumps the"
          " sea in momentum; emissive first relaxes it but has no H-theorem.")

    rule("B3. Mean field: an aligned pair as reader, the firing pair itself moved")
    print("  A pair A at row n fires channel q (rate g, orientation t = +-1) and moves to"
          " n + tq;\n  a partner moves the other way, -tq.  The partner is on A's own row"
          " (defocus: it goes\n  to n - tq) or on row n + 2tq (focus: both land on n + tq),"
          " chosen in proportion to\n  their counts; h(m) is the chance that some partner"
          " is present, m the partners in reach.")

    def rhs2(s, avail):
        ds = np.zeros(n)
        for t in (+1, -1):
            same, far = s, np.roll(s, -2 * t * q)          # far[i] = s[i + 2tq]
            tot = same + far
            h = avail(tot)
            wd = same / np.maximum(tot, 1e-300)
            rd, rf = 0.5 * g * s * h * wd, 0.5 * g * s * h * (1 - wd)
            ds += -2 * rd + np.roll(rd, t * q) + np.roll(rd, -t * q)
            ds += -rf - np.roll(rf, 2 * t * q) + 2 * np.roll(rf, t * q)
        return ds
    print(f"\n  {'partner law':>30} {'s':>4} {'drift(uniform)':>15} {'max growth':>11}"
          f" {'30% dip, t = 40: min':>21} {'max':>6} {'d(pairs)':>9}")
    for name, avail in (("mass action, h(m) = m", lambda m: m),
                        ("integer sea, h(m) = 1 - e^-m", lambda m: 1 - np.exp(-m))):
        for s0 in (1.0, 4.0):
            s = np.full(n, s0)
            f0 = rhs2(s, avail)
            jac = np.zeros((n, n))
            for k in range(n):
                ds = np.zeros(n)
                ds[k] = 1e-6
                jac[:, k] = (rhs2(s + ds, avail) - f0) / 1e-6
            ev = np.linalg.eigvals(jac).real.max()
            s = np.full(n, s0)
            s[n // 2] *= 0.7
            m0 = s.sum()
            for _ in range(4000):
                k1 = rhs2(s, avail)
                k2 = rhs2(s + 0.01 * k1, avail)
                s = np.maximum(s + 0.005 * (k1 + k2), 0)
            print(f"  {name:>30} {s0:4.0f} {np.abs(f0).max():15.1e} {ev:+11.4f}"
                  f" {s.min() / s0:21.3f} {s.max() / s0:6.3f} {s.sum() - m0:+9.1e}",
                  flush=True)
    print("  Two-body: defocus goes as s_n s_n and focus as s_(n-q) s_(n+q), the form of"
          " Definition (C).\n  The uniform row is an exact equilibrium, no ripple grows,"
          " and a dip refills.")


# ------------------------------------------------------------------ Part C
def window_pairs(xs, W):
    """Unordered pairs (i, j) with |x_i - x_j| < W on the periodic box."""
    n = xs.size
    e = np.empty(0, np.int64)
    if n < 2:
        return e, e
    order = np.argsort(xs)
    xo = xs[order]
    ext = np.concatenate([xo, xo + L])
    hi = np.minimum(np.searchsorted(ext, xo + W, side="left"), np.arange(n) + n)
    cnt = hi - np.arange(n) - 1
    if cnt.sum() == 0:
        return e, e
    i = np.repeat(np.arange(n), cnt)
    start = np.repeat(np.cumsum(cnt) - cnt, cnt)
    j = (np.arange(i.size) - start) + i + 1
    return order[i], order[j % n]


def particle_run(cfg):
    nu, seed, kappa, mode, scale = (cfg["nu"], cfg["seed"], cfg["kappa"],
                                    cfg["coll"], cfg["scale"])
    Wc = ARGS.coll_window
    rng = np.random.default_rng(seed)
    rng_c = np.random.default_rng((seed, 4242))        # collisions' own stream
    n_pairs = int(round(nu * B * L * 2 * PMAX))
    xs = rng.uniform(-L / 2, L / 2, n_pairs)
    ps = DP * np.clip(np.round(rng.uniform(-PMAX, PMAX, n_pairs) / DP),
                      -(NP // 2) + 1, NP // 2 - 1)
    r0, p0, sr, spk, rho = -8.0, 1.2, 2.0, 0.25, 6.0
    n_pos, n_neg = int(round(nu * (rho + 1) / 2)), int(round(nu * (rho - 1) / 2))
    xa, xb_ = rng.normal(r0, sr, n_pos), rng.normal(r0, sr, n_neg)
    pa, pb = rng.normal(p0, spk, n_pos), rng.normal(p0, spk, n_neg)
    x = np.concatenate([xs, xs, xa, xb_])
    p = np.concatenate([ps, ps, pa, pb])
    th = np.concatenate([ps * xs / HBAR, ps * xs / HBAR,
                         p0 * (xa - r0) / HBAR, p0 * (xb_ - r0) / HBAR])
    eps = np.concatenate([np.ones(n_pairs), -np.ones(n_pairs),
                          np.ones(n_pos), -np.ones(n_neg)])
    mate = np.concatenate([np.arange(n_pairs, 2 * n_pairs), np.arange(n_pairs),
                           -np.ones(n_pos + n_neg, int)])
    N = len(x)
    n_row = nu * B * DP
    st = dict(recomb=0, ion=0, fail=0, defocus=0, focus=0, rejected=0,
              gray_made=0.0, id_res=0.0, dark_events=0, chir_n=0, chir_s=0.0,
              relocks=0)

    def find(mask, xc, pc):
        ok = mask & (np.abs(x - xc) < 1.0) & (np.abs(p - pc) < DP / 2)
        cand = np.flatnonzero(ok)
        return int(rng.choice(cand)) if len(cand) else -1

    def settle(par, q, sgn):
        t_ = sgn * eps[par]
        xi = XI[q]
        free_m = (mate < 0) & (np.arange(N) != par)
        jn = find(free_m & (eps < 0), x[par], p[par] + t_ * xi)
        jp = find(free_m & (eps > 0), x[par], p[par] - t_ * xi)
        if jn >= 0 and jp >= 0:
            p[jn] -= t_ * xi
            p[jp] += t_ * xi
            xm, pm = 0.5 * (x[jn] + x[jp]), 0.5 * (p[jn] + p[jp])
            for j in (jn, jp):
                th[j] += p[j] * (xm - x[j]) / HBAR
                x[j], p[j] = xm, pm
            mate[jn], mate[jp] = jp, jn
            st["recomb"] += 1
            return
        j = find((mate >= 0) & (eps > 0), x[par], p[par])
        if j >= 0:
            k = mate[j]
            p[j] += t_ * xi
            p[k] -= t_ * xi
            mate[j] = mate[k] = -1
            st["ion"] += 1
            return
        st["fail"] += 1

    def dark_catalysis():
        if kappa <= 0:
            return
        A_all = np.flatnonzero((mate >= 0) & (eps > 0) & (np.abs(x) < XACT))
        g = kappa * np.interp(x[A_all], run.r, GAM)
        for a in A_all[rng.random(len(A_all)) < g * DT]:
            if mate[a] < 0:
                continue
            cand = np.flatnonzero((mate >= 0) & (eps > 0) & (np.abs(x - x[a]) < 2 * YMAX)
                                  & (np.abs(p - p[a]) < DP / 2))
            cand = cand[cand != a]
            if not len(cand):
                continue
            c = int(rng.choice(cand))
            phi_a = th[a] + p[a] * (x[c] - x[a]) / HBAR
            new = float(np.angle(np.exp(1j * phi_a) + np.exp(1j * th[c])))
            th[a] = th[mate[a]] = new - p[a] * (x[c] - x[a]) / HBAR
            th[c] = th[mate[c]] = new
            st["dark_events"] += 1

    def gray(i):
        """sin^2(mu/2) of the pair whose positon member is i (|Psi|^2 / 4)."""
        m = mate[i]
        mu = th[i] - th[m] + p[i] * (x[m] - x[i]) / HBAR
        return np.sin(0.5 * mu) ** 2

    def collide():
        if mode == "none":
            return
        A = np.flatnonzero((mate >= 0) & (eps > 0) & (np.abs(x) < XACT))
        if A.size < 2:
            return
        I, J = window_pairs(x[A], Wc)
        if I.size == 0:
            return
        I, J = A[I], A[J]
        dr = np.round((p[J] - p[I]) / DP).astype(int)
        keep = (dr % 2 == 0) & (np.abs(dr) // 2 <= NQ)
        I, J, dr = I[keep], J[keep], dr[keep]
        xm = x[I] + 0.5 * ((x[J] - x[I] + L / 2) % L - L / 2)
        dmask = dr == 0
        qf = np.abs(dr) // 2
        g = np.interp(xm, run.r, GAM)
        kf = np.zeros(I.size)
        if (~dmask).any():
            kf[~dmask] = np.abs(interp_rows(KQ, xm[~dmask])[np.arange((~dmask).sum()),
                                                            qf[~dmask] - 1])
        fac = scale * YMAX / (Wc * nu)            # Definition (C) at density B
        rate = np.where(dmask, fac * g, fac * kf / 2)
        fire = np.flatnonzero(rng_c.random(I.size) < -np.expm1(-rate * DT))
        used = set()
        for k in rng_c.permutation(fire):
            a, b = int(I[k]), int(J[k])
            if a in used or b in used or mate[a] < 0 or mate[b] < 0:
                continue
            used.update((a, b))
            am, bm = int(mate[a]), int(mate[b])
            g0 = gray(a) + gray(b)
            if dr[k] == 0:                        # defocus on channel q
                kk = interp_rows(np.abs(KQ), np.array([xm[k]]))[0]
                if kk.sum() <= 0:
                    continue
                q = int(rng_c.choice(NQ, p=kk / kk.sum()))
                xi = XI[q]
                if mode in ("whole", "relock"):
                    o = rng_c.choice((-1, 1))
                    p[a] += o * xi
                    p[am] += o * xi
                    p[b] -= o * xi
                    p[bm] -= o * xi
                    if mode == "relock":
                        relock(a)
                        relock(b)
                else:
                    pr = p[a]
                    xa_, xb2 = float(x[a]), float(x[b])

                    def products(o):
                        # (a+, b-) -> pr + o xi at a's site; (b+, a-) -> pr - o xi
                        # at b's site.  The set of pair positions is unchanged, so
                        # the sea is not clumped; mu of an equal-momentum pair is
                        # the same at any common point.
                        out = {}
                        for bod, pn, xn in ((a, pr + o * xi, xa_), (bm, pr + o * xi, xa_),
                                            (b, pr - o * xi, xb2), (am, pr - o * xi, xb2)):
                            dxn = (xn - x[bod] + L / 2) % L - L / 2
                            out[bod] = (pn, th[bod] + pn * dxn / HBAR, xn)
                        mu1 = out[a][1] - out[bm][1]
                        mu2 = out[b][1] - out[am][1]
                        return out, mu1, mu2
                    if mode == "exchange":
                        o = int(rng_c.choice((-1, 1)))
                    else:
                        wts = []
                        for o_ in (1, -1):
                            _, m1, m2 = products(o_)
                            wts.append(np.cos(0.5 * m1) ** 2 * np.cos(0.5 * m2) ** 2)
                        u = rng_c.random()
                        if u < wts[0] / 2:
                            o = 1
                        elif u < (wts[0] + wts[1]) / 2:
                            o = -1
                        else:
                            st["rejected"] += 1
                            continue
                    # chirality: does the orientation follow the local drive?
                    # o = +1 sends a's positon up.  Small-phase phase matching
                    # predicts o = -sign(mu_chord^(p) xi d), and under (S') the
                    # chord winds at U_res(xm, d/2), so the sign to compare is
                    # that of U_res(xm, d/2) d.
                    dsep0 = (x[b] - x[a] + L / 2) % L - L / 2
                    dvm = np.interp(xm[k], run.r, run.dv_eff)
                    ur = (V(xm[k] + dsep0 / 2) - V(xm[k] - dsep0 / 2) - dsep0 * dvm)
                    if abs(ur * dsep0) > 0:
                        st["chir_n"] += 1
                        st["chir_s"] += o * np.sign(ur * dsep0)
                    # the exchange identity of Part A, event by event
                    muA = th[a] - th[am]
                    muB = th[b] - th[bm]
                    dsep = (x[b] - x[a] + L / 2) % L - L / 2
                    out, m1, m2 = products(o)
                    res = wrap(m1 + m2 - muA - muB - 2 * o * xi * dsep / HBAR)
                    st["id_res"] = max(st["id_res"], abs(float(res)))
                    for bod, (pn, tn, xn) in out.items():
                        p[bod], th[bod], x[bod] = pn, tn, xn
                    mate[a], mate[bm] = bm, a
                    mate[b], mate[am] = am, b
                st["defocus"] += 1
            else:                                 # focus onto the middle row
                xi = XI[qf[k] - 1]
                up, dn = (a, b) if dr[k] > 0 else (b, a)     # up: on n + q
                if mode in ("whole", "relock"):
                    for bod in (up, mate[up]):
                        p[bod] -= xi
                    for bod in (dn, mate[dn]):
                        p[bod] += xi
                    if mode == "relock":
                        relock(up)
                        relock(dn)
                else:
                    upm, dnm = int(mate[up]), int(mate[dn])
                    pn = p[up] - xi
                    site = {up: float(x[up]), dnm: float(x[up]),
                            dn: float(x[dn]), upm: float(x[dn])}
                    new = {bod: th[bod] + pn * (((site[bod] - x[bod] + L / 2) % L
                                                 - L / 2)) / HBAR
                           for bod in (up, upm, dn, dnm)}
                    mu1 = new[up] - new[dnm]                  # (up+, dn-)
                    mu2 = new[dn] - new[upm]                  # (dn+, up-)
                    if mode == "matched":
                        wt = np.cos(0.5 * mu1) ** 2 * np.cos(0.5 * mu2) ** 2
                        if rng_c.random() >= wt:
                            st["rejected"] += 1
                            continue
                    for bod in (up, upm, dn, dnm):
                        p[bod], th[bod], x[bod] = pn, new[bod], site[bod]
                    mate[up], mate[dnm] = dnm, up
                    mate[dn], mate[upm] = upm, dn
                st["focus"] += 1
            st["gray_made"] += gray(a) + gray(b) - g0

    def relock(i):
        """Step 22 section 6: the pair whose positon is i takes the circular mean
        of its destination row's sea clocks within 3, each carried to it by its
        own momentum; both members, so the pair is dark and on the row's lock."""
        m = ((mate >= 0) & (eps > 0) & (np.abs(p - p[i]) < DP / 2)
             & (np.abs((x - x[i] + L / 2) % L - L / 2) < 3.0))
        m[i] = False
        if not m.any():
            return
        dxm = (x[i] - x[m] + L / 2) % L - L / 2
        ref = float(np.angle(np.exp(1j * (th[m] + p[m] * dxm / HBAR)).sum()))
        th[i] = th[mate[i]] = ref
        st["relocks"] += 1

    def readings(bs):
        sea = np.flatnonzero((mate >= 0) & (eps > 0))
        rows = np.round(p[sea] / DP).astype(int)
        out = np.zeros((len(bs), NQ))
        dv = np.interp(x[bs], run.r, run.dv_eff)
        for k, b in enumerate(bs):
            r = int(np.round(p[b] / DP))
            m = sea[(rows == r) & (np.abs(x[sea] - x[b]) < YMAX + 0.5)]
            if len(m) < 2:
                continue
            m = m[np.argsort(x[m])]
            I, J = np.triu_indices(len(m), 1)
            d = x[m[J]] - x[m[I]]
            xmc = 0.5 * (x[m[I]] + x[m[J]])
            sel = (d <= 2 * YMAX) & (np.abs(xmc - x[b]) < 0.5)
            if not sel.any():
                continue
            i, j, d = m[I[sel]], m[J[sel]], d[sel]
            w = np.cos(np.pi * d / (4 * YMAX)) ** 2
            rate = V(x[j]) - V(x[i]) - d * dv[k]
            mu = (th[i] - th[j] + r * DP * d / HBAR)[:, None]
            kh = -(B * DP / (HBAR * n_row ** 2)) * (w * rate) @ np.sin(
                XI[None, :] * d[:, None] / HBAR + mu)
            out[k] = kh - XI * (XI @ kh) / (XI @ XI)
        return out

    def chord_coherence(bs):
        sea = np.flatnonzero((mate >= 0) & (eps > 0))
        rows = np.round(p[sea] / DP).astype(int)
        z, n = 0j, 0
        for b in bs:
            r = int(np.round(p[b] / DP))
            m = sea[(rows == r) & (np.abs(x[sea] - x[b]) < YMAX + 0.5)]
            if len(m) < 2:
                continue
            I, J = np.triu_indices(len(m), 1)
            i, j = m[I], m[J]
            d = x[j] - x[i]
            sel = (np.abs(d) <= 2 * YMAX) & (np.abs(0.5 * (x[i] + x[j]) - x[b]) < 0.5)
            i, j, d = i[sel], j[sel], d[sel]
            z += np.exp(1j * (th[i] - th[j] + r * DP * d / HBAR)).sum()
            n += len(d)
        return z, n

    XB = np.linspace(-XACT, XACT, 17)
    sea0 = None
    acc_live = np.zeros((N, NQ))
    acc_mesh = np.zeros((N, NQ))
    trace = []
    coh = [0j, 0]
    t0 = time.perf_counter()
    nstep = int(round(ARGS.t_end / DT))
    for step in range(nstep + 1):
        if step % 25 == 0 or step == nstep:
            A = np.flatnonzero((mate >= 0) & (eps > 0) & (np.abs(x) < XACT))
            gr = np.array([gray(i) for i in A]) if A.size else np.zeros(1)
            h2 = np.histogram2d(x[A], np.round(p[A] / DP),
                                bins=(XB, np.arange(-31.5, 32.5)))[0]
            mc = h2.mean()
            # structure beyond Poisson: (variance - mean) / mean^2 over cells
            dev = float((h2.var() - mc) / max(mc * mc, 1e-300))
            _ = sea0
            fr = mate < 0
            trace.append((step * DT, float(gr.mean()), float((gr > 0.1).mean()), dev,
                          float(eps[fr & (x > 0)].sum() / eps[fr].sum())))
        if step == nstep:
            break
        kick = np.where(mate >= 0, 0.0, 1.0)
        p += 0.5 * DT * F(x) * kick
        x += p / MU * DT
        p += 0.5 * DT * F(x) * kick
        xw = (x + L / 2) % L - L / 2
        th += p * (xw - x) / HBAR
        x[:] = xw
        th += (p ** 2 / (2 * MU) - V(x)) * DT / HBAR
        free = np.flatnonzero((mate < 0) & (np.abs(x) < XACT))
        if len(free):
            km = interp_rows(KQ, x[free])
            if ARGS.shadow:
                acc_live[free] += readings(free) * DT
                acc_mesh[free] += km * DT
            if step % 10 == 0:
                near = free[np.abs(x[free]) < 3.0]
                if len(near) > 40:
                    near = rng.choice(near, 40, replace=False)
                z_, n_ = chord_coherence(near)
                coh[0] += z_
                coh[1] += n_
            g = np.abs(km)
            fires = []
            for k, b in enumerate(free):
                for _ in range(rng.poisson(g[k].sum() * DT)):
                    q = rng.choice(NQ, p=g[k] / g[k].sum())
                    fires.append((b, q, np.sign(km[k, q])))
            for o in rng.permutation(len(fires)):
                b, q, sgn = fires[o]
                if mate[b] < 0:
                    settle(b, q, sgn)
        dark_catalysis()
        collide()

    out = dict(cfg)
    out["window"] = Wc
    out["wall"] = round(time.perf_counter() - t0, 1)
    tr = np.array(trace)
    ncoll = max(st["defocus"] + st["focus"], 1)
    out.update(gray_end=tr[-1, 1], gray_frac_end=tr[-1, 2], gray_max=tr[:, 1].max(),
               sea_dev_end=tr[-1, 3], sea_dev_mean=float(tr[1:, 3].mean()), T_E=tr[-1, 4],
               gray_per_coll=st["gray_made"] / ncoll,
               chirality=st["chir_s"] / max(st["chir_n"], 1),
               chord_coh=abs(coh[0]) / max(coh[1], 1), **st)
    if ARGS.shadow:
        m = np.abs(acc_mesh).sum(axis=1) > 0
        Al, Am = acc_live[m], acc_mesh[m]
        out.update(sh_slope=float((Al * Am).sum() / (Am * Am).sum()),
                   sh_corr=float(np.corrcoef(Al.ravel(), Am.ravel())[0, 1]))
    out["_trace"] = [(cfg["coll"], cfg["kappa"], cfg["seed"], *r) for r in trace]
    return out


def part_c():
    rule("C. Particle model with clocks: three collision realisations")
    cfgs = [dict(nu=nu, seed=11 + s, kappa=k, coll=c, scale=sc)
            for nu in ARGS.nu for k in ARGS.kappa for c in ARGS.coll
            for sc in ARGS.coll_scale for s in range(ARGS.seeds)
            if not (c == "none" and sc != ARGS.coll_scale[0])]
    workers = max(1, min(ARGS.workers or os.cpu_count() or 1, len(cfgs)))
    print(f"   {len(cfgs)} runs on {workers} worker(s); t_end {ARGS.t_end},"
          f" window {ARGS.coll_window}, shadow {ARGS.shadow}", flush=True)
    res = []
    t0 = time.perf_counter()
    with mp.Pool(workers) as pool:
        for r in pool.imap_unordered(particle_run, cfgs):
            res.append(r)
            print(f"   nu {r['nu']} kappa {r['kappa']:g} {r['coll']:>8} x{r['scale']:g}"
                  f" seed {r['seed']}: gray {r['gray_end']:.4f} (>0.1: {r['gray_frac_end']:.3f})"
                  f" coh {r['chord_coh']:.3f}"
                  + (f" shadow {r['sh_slope']:.3f}/{r['sh_corr']:.3f}" if ARGS.shadow else "")
                  + f" T_E {r['T_E']:.3f} fail {r['fail']} def {r['defocus']} foc {r['focus']}"
                  f" rej {r['rejected']} id {r['id_res']:.1e} dev {r['sea_dev_mean']:.4f}"
                  f" g/coll {r['gray_per_coll']:.4f} chir {r['chirality']:+.3f}"
                  f" ({r['wall']:.0f} s; {time.perf_counter() - t0:.0f} s)", flush=True)
            write(res)
    summarise(res)


def write(res):
    keys = [k for k in res[0] if not k.startswith("_")]
    with open(output_path(f"pair_collision_clocks_C{ARGS.tag}.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in res:
            w.writerow({k: r.get(k, "") for k in keys})
    with open(output_path(f"pair_collision_clocks_trace{ARGS.tag}.csv"), "w",
              newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["coll", "kappa", "seed", "t", "gray_mean", "gray_frac", "sea_dev", "T_E"])
        for r in res:
            w.writerows(r["_trace"])


def summarise(res):
    print("\n   mean +- standard error over seeds")
    groups = {}
    for r in res:
        groups.setdefault((r["nu"], r["kappa"], r["coll"], r["scale"]), []).append(r)
    cols = ["gray_end", "gray_frac_end", "chord_coh", "gray_per_coll", "chirality",
            "sea_dev_mean", "T_E", "fail"]
    if ARGS.shadow:
        cols += ["sh_slope", "sh_corr"]
    print("   " + f"{'nu':>3} {'kappa':>5} {'coll':>9} {'x':>4} " + " ".join(f"{c:>17}" for c in cols))
    for k, rs in sorted(groups.items()):
        cells = []
        for c in cols:
            v = np.array([float(r[c]) for r in rs])
            se = v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else 0.0
            cells.append(f"{v.mean():9.4f}+-{se:6.4f}")
        print("   " + f"{k[0]:>3} {k[1]:>5g} {k[2]:>9} {k[3]:>4g} " + " ".join(f"{c:>17}" for c in cells))


if __name__ == "__main__":
    if "A" in ARGS.parts:
        part_a()
    if "B" in ARGS.parts:
        part_b()
    if "C" in ARGS.parts:
        part_c()
