"""The sea as a resonance clock -- companion to step 23's third addendum.

Context (docs/analysis/force_blind_sea.md).  The compensated algorithm fixes
the RATE of each event channel, |K_q(x)| per parent, from a precomputed field;
it says nothing about the MECHANISM.  Proposition Q5 reads K_q as the rate of
change of an interference sum over the sea's own clocks.  This demo asks
whether that reading can serve as a deterministic mechanism: every body
carries one integrator per channel, fed by its live reading of the sea, and
fires channel q each time the integrator crosses an integer.

  Part A  exact checks (SymPy and the mesh generator)
          A1 gauge: under V -> V + c a single clock read against a plane
             wave drifts at -c/hbar; a chord (two clocks) does not.
          A2 kinetics of one event: emission and absorption change the
             kinetic energy of all bodies by +-xi^2/m = (xi^2/2m) dN, and the
             W-weighted kinetic energy by the same 2 p xi / m.
          A3 a crossing-triggered absorption gives
             dE_a/dt = kappa (N_a E_b - E_a N_b)/2, which depends on the body
             density N: E is no longer closed (Theorem S2 fails).
          A4 the event-rate convention: channel q >= 1 firing at |K_q| per
             parent reproduces the QLE jump generator; firing at 2|K_q|
             (gamma_tot = sum over q != 0, the published particle runs) does
             not.
  Part B  integrators at one site of the Eckart flank
          B1 per-channel firing rates with banked integrators against one
             reset of the whole aperture on any firing.
          B2 counting statistics (Fano factor of counts in unit windows):
             deterministic clocks are sub-Poissonian.
  Part C  the live reading on a static Poisson sea: unbiased against K_q.
  Part D  open loop, particle model (S') -- force-blind row-centred locked
          sea, wrap phase, reach dark catalysis at kappa x the kernel's own
          rate sum_{q>=1}|K_q|.  Free bodies stream (no events); each
          integrates its live reading Khat_q and the mesh K_q along its path.
          Reported: slope / corr / relative error of int Khat against int K
          (and against K averaged over the aperture window), the mu = 0
          control, and the first moment the raw reading carries before the
          projection of Proposition L2(c).
  Part E  closed loop: events triggered by
            poisson       Poisson at |K_q| (mesh)
            clock-mesh    integrators fed by the mesh K_q
            clock-live    integrators fed by the live reading
            clock-live-h  the same with a hysteresis band (Schmitt trigger)
          and realised as in demo_sea_lock_particles.py.  Reported: the
          trigger's signed firings against the QLE target sum eps K_q dt in
          (x, q) bins, gross firings against the QLE's gross rate, the
          transmission T_E against the mesh QLE, and event counts.

Units hbar = m = 1; Eckart barrier V0 = a = 1; dp = 0.25, n_p = 64,
y_max = 2 pi; box L = 48 with periodic wrap (a window on the open line).

Usage::

    PYTHONPATH=src python3 -u src/demo_sea_resonance_clock.py   # all parts, quick
    # sweeps (Kaggle): seeds in parallel, one CSV row per run
    PYTHONPATH=src python3 -u src/demo_sea_resonance_clock.py --parts D \\
        --nu 32 64 128 --dX 0.5 1.0 --kappa 1 2 --norm local --seeds 8 --tag _sweep
    PYTHONPATH=src python3 -u src/demo_sea_resonance_clock.py --parts E \\
        --nu-closed 32 64 --norm local --seeds 8 --tag _sweep

Writes sea_resonance_clock_<part><tag>.csv and, when Part D covers two or
more nu, sea_resonance_clock<tag>.png, through wpmwlib.wpmw_utils.
"""
import argparse
import csv
import multiprocessing as mp
import os
import time

import numpy as np

from demo_emission_and_absorption import Ledger, packet, B, MU, HBAR
from wpmwlib.wpmw_utils import output_path, docs_path

ap = argparse.ArgumentParser(description=__doc__,
                             formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--parts", default="ABCDE")
ap.add_argument("--nu", type=int, nargs="+", default=[8, 16],
                help="ensemble multiplicities for Part D")
ap.add_argument("--nu-closed", type=int, nargs="+", default=[8],
                help="ensemble multiplicities for Part E")
ap.add_argument("--dX", type=float, nargs="+", default=[1.0],
                help="aperture window: chord midpoints within dX/2 of the reader")
ap.add_argument("--kappa", type=float, nargs="+", default=[2.0],
                help="dark-catalysis rate, in units of sum_{q>=1}|K_q|")
ap.add_argument("--norm", choices=("expected", "local"), default="expected",
                help="normalise the reading by the expected sea density nu B dp, or by "
                     "the count in the reader's aperture (better at nu >= 32, unstable "
                     "at small nu, where 1/n_loc^2 amplifies the noise)")
ap.add_argument("--events", nargs="+",
                default=["poisson", "clock-mesh", "clock-live", "clock-live-h"])
ap.add_argument("--seeds", type=int, default=2, help="seeds 11, 12, ...")
ap.add_argument("--t-end", type=float, default=14.0)
ap.add_argument("--workers", type=int, default=0, help="0: one per CPU")
ap.add_argument("--tag", default="")
ARGS = ap.parse_args()

# ------------------------------------------------------------------ setup
run = Ledger(v0=1.0, a=1.0, n_r=192, r_half=24.0, n_p=64, dp=0.25)
L, DP, NP = 48.0, run.dp, run.n_p
PMAX, YMAX = NP * DP / 2, run.y_max
Q = np.arange(1, NP // 2)
XI = Q * DP
NQ = len(Q)
KQ = run.k[:, Q]                           # signed rate of channel q >= 1
GAM = np.abs(KQ).sum(axis=1)               # per-parent event rate (Part A4)
KCUM = np.vstack([np.zeros(NQ), np.cumsum(0.5 * (KQ[1:] + KQ[:-1]), axis=0) * run.dr])
XACT = 8.0                                 # K_q is below 1e-3 of its peak beyond 7.5
DT = 0.02
X_SITE = -0.625                            # the Eckart flank site of Q5 and Q7


def V(x):
    return 1.0 / np.cosh(x) ** 2


def F(x):                                  # -V'(x)
    return 2.0 * np.tanh(x) / np.cosh(x) ** 2


def interp_rows(vals, xs):
    idx = np.clip(np.searchsorted(run.r, xs) - 1, 0, len(run.r) - 2)
    w = ((xs - run.r[idx]) / run.dr)[:, None]
    return vals[idx] * (1 - w) + vals[idx + 1] * w


def k_window(xs, dX):
    """K_q averaged over chord midpoints uniform in [x - dX/2, x + dX/2]."""
    return (interp_rows(KCUM, xs + dX / 2) - interp_rows(KCUM, xs - dX / 2)) / dX


def rule(title):
    print("\n" + "=" * 78 + f"\n{title}\n" + "=" * 78, flush=True)


# ------------------------------------------------------------------ Part A
def part_a():
    import sympy as sp
    rule("A. Exact checks")
    t, hb, m, c, p, xi = sp.symbols("t hbar m c p xi", real=True)
    x1, x2 = sp.symbols("x1 x2", real=True)
    Vf = sp.Function("V")
    xk = lambda x0: x0 + p * t / m                       # (S'): inertial
    th_rate = lambda x0: (p ** 2 / (2 * m) - Vf(xk(x0)) - c) / hb
    pw_rate = (p * (p / m) - p ** 2 / (2 * m)) / hb      # plane wave along x_k(t)
    single = sp.simplify(th_rate(x1) - pw_rate)
    chord = sp.simplify(th_rate(x1) - th_rate(x2))
    print(f"A1 gauge V -> V + c:  d/dt[theta - plane wave] = {single}")
    print(f"                      d/dt[mu_ij]              = {chord}")
    assert sp.diff(single, c) != 0 and sp.diff(chord, c) == 0
    print("   a single clock depends on c; a chord does not.")

    tau = sp.symbols("tau")
    T = lambda q_: q_ ** 2 / (2 * m)
    print("A2 one event, orientation tau = +-1 (tau^2 = 1):")
    for name, d_all, dN in (("emission  ", T(p + tau * xi) + T(p - tau * xi) - 2 * T(p), 2),
                            ("absorption", 2 * T(p) - T(p + tau * xi) - T(p - tau * xi), -2)):
        d_all = sp.expand(d_all).subs(tau ** 2, 1)
        d_w = sp.factor(sp.expand(T(p + tau * xi) - T(p - tau * xi)))
        resid = sp.simplify(d_all - xi ** 2 / (2 * m) * dN)
        print(f"   {name}: dT(all bodies) = {d_all},  dT_W = {d_w},  dN = {dN:+d},"
              f"  dT - (xi^2/2m) dN = {resid}")
        assert resid == 0

    Na, Nb, Ea, Eb, kk = sp.symbols("N_a N_b E_a E_b kappa", real=True)
    up = lambda N, E: (N + E) / 2
    um = lambda N, E: (N - E) / 2
    dEa = sp.factor(sp.expand(kk * (um(Na, Ea) * up(Nb, Eb) - up(Na, Ea) * um(Nb, Eb))))
    print(f"A3 crossing-triggered absorption, rows a = p + xi, b = p - xi:  dE_a/dt = {dEa}")
    assert sp.diff(dEa, Na) != 0
    print("   depends on N_a: the observable is no longer closed in E.")

    rng = np.random.default_rng(0)
    E = rng.normal(size=(len(run.r), NP))
    gen = np.real(np.fft.ifft(np.fft.fft(E, axis=1) * run.sym_e, axis=1))
    print("A4 per-parent rate of channel q >= 1 (one event, legs +-xi_q) against the"
          " QLE jump generator:")
    for c_ in (1.0, 2.0):
        ev = sum(c_ * run.k[:, q][:, None] * (np.roll(E, q, axis=1) - np.roll(E, -q, axis=1))
                 for q in Q)
        err = np.linalg.norm(ev - gen) / np.linalg.norm(gen)
        print(f"   rate {c_:.0f} x |K_q|:  relative error {err:.2e}")
    i = np.argmin(np.abs(run.r - X_SITE))
    print(f"   gamma_tot / sum_(q>=1)|K_q| = {run.gamma_tot[i] / GAM[i]:.6f}")


# ------------------------------------------------------------------ Part B
def site_rates():
    return np.abs(interp_rows(KQ, np.array([X_SITE]))[0])


def part_b():
    rule(f"B. Integrators at x = {X_SITE} (lambda_q = |K_q|, q = 1..{NQ})")
    lam = site_rates()
    G = lam.sum()
    print(f"   sum |K_q| = {G:.4f};  |K_1..4| = {np.round(lam[:4], 4)}")
    rng = np.random.default_rng(1)
    T_run = 400.0
    # B1 banked: each channel counts its own integer crossings
    phi0 = rng.uniform(0, 1, NQ)
    banked = np.floor(phi0 + lam * T_run) / T_run
    # B1 one reset of the whole aperture on any firing (event-driven, exact)
    phi = rng.uniform(0, 1, NQ)
    cnt = np.zeros(NQ)
    t = 0.0
    while True:
        wait = (1.0 - phi) / lam
        q = int(np.argmin(wait))
        t += wait[q]
        if t > T_run:
            break
        cnt[q] += 1
        phi[:] = 0.0
    reset = cnt / T_run
    print("B1 firing rate / |K_q|, q = 1..6, and the total against sum |K_q|:")
    for name, r in (("banked integrators  ", banked), ("whole-aperture reset", reset)):
        print(f"   {name}: " + " ".join(f"{a / b:5.2f}" for a, b in zip(r[:6], lam[:6]))
              + f"   total {r.sum() / G:5.2f}")

    def counts(lam_, n_bodies):
        """Counts in unit windows for n_bodies independent banked clocks."""
        ph = rng.uniform(0, 1, (n_bodies, len(lam_)))
        k = np.arange(int(T_run) + 1)
        c_ = np.floor(ph[..., None] + lam_[None, :, None] * k[None, None, :])
        return np.diff(c_, axis=2).sum(axis=(0, 1))

    print("B2 counts in unit windows: mean and Fano factor (Poisson = 1):")
    rows = []
    for name, lam_, nb in (("channel q = 1, one body", lam[:1], 1),
                           ("channel q = 1, 64 bodies", lam[:1], 64),
                           ("all channels, one body", lam, 1),
                           ("all channels, 64 bodies", lam, 64)):
        c_ = counts(lam_, nb)
        rows.append((name, c_.mean(), c_.var() / c_.mean()))
        print(f"   {name:26s} mean {c_.mean():7.3f}   Fano {c_.var() / c_.mean():.3f}")
    pois = rng.poisson(G * 64, 400)
    print(f"   {'Poisson at the same rate':26s} mean {pois.mean():7.3f}   Fano "
          f"{pois.var() / pois.mean():.3f}")
    return dict(lam=lam, banked=banked, reset=reset, fano=rows)


# ------------------------------------------------------------------ Part C
def reading_from(xs_sorted, xb, dX, n_norm, dv, P=0.0, th=None):
    """Live reading Khat_q (Theorem L1 as a pair sum, reader's lever P, L8).

    xs_sorted: positions of one member of each sea pair in the reader's row,
    sorted.  th: their clocks (None = mu = 0).  Returns Khat (NQ,), the raw
    first moment before projection, and the chord count."""
    I, J = np.triu_indices(len(xs_sorted), 1)
    d = xs_sorted[J] - xs_sorted[I]
    xm = 0.5 * (xs_sorted[I] + xs_sorted[J])
    sel = (d <= 2 * YMAX) & (np.abs(xm - xb) < dX / 2)
    if not sel.any():
        return np.zeros(NQ), 0.0, 0
    I, J, d = I[sel], J[sel], d[sel]
    w = np.cos(np.pi * d / (4 * YMAX)) ** 2
    rate = V(xs_sorted[J]) - V(xs_sorted[I]) - d * dv          # hbar dmu/dt, (S') rows
    mu = 0.0 if th is None else (th[I] - th[J] + P * d / HBAR)[:, None]
    kh = -(B * DP / (HBAR * n_norm ** 2 * dX)) * (w * rate) @ np.sin(
        XI[None, :] * d[:, None] / HBAR + mu)
    fm = 2.0 * XI @ kh
    kh = kh - XI * (XI @ kh) / (XI @ XI)                         # Proposition L2(c)
    return kh, fm, int(sel.sum())


def part_c():
    rule("C. The live reading on a static Poisson sea (mu = 0): mean Khat_q / K_q")
    rng = np.random.default_rng(3)
    for xb in (-0.75, -1.5):
        K = interp_rows(KQ, np.array([xb]))[0]
        dv = np.interp(xb, run.r, run.dv_eff)
        for nu, dX, M in ((8, 1.0, 4000), (64, 1.0, 300), (64, 0.25, 600)):
            n = nu * B * DP
            Kw = k_window(np.array([xb]), dX)[0]
            acc = np.zeros(NQ)
            for _ in range(M):
                cnt = rng.poisson(n * (2 * YMAX + dX))
                xs = np.sort(rng.uniform(xb - YMAX - dX / 2, xb + YMAX + dX / 2, cnt))
                kh, _, _ = reading_from(xs, xb, dX, n, dv)
                acc += kh
            est = acc / M
            print(f"   x = {xb:5.2f}  nu = {nu:3d}  dX = {dX:4.2f}:  q = 1..4 against K "
                  f"{np.round(est[:4] / K[:4], 3)}  against window-mean K "
                  f"{np.round(est[:4] / Kw[:4], 3)}", flush=True)


# ------------------------------------------------------------------ Parts D, E
def particle_run(cfg):
    """One run of the particle model.  cfg: dict(nu, seed, kappa, dX, events)."""
    nu, seed, kappa, dX, events = (cfg["nu"], cfg["seed"], cfg["kappa"], cfg["dX"],
                                   cfg["events"])
    rng = np.random.default_rng(seed)
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
    stats = dict(recomb=0, ion=0, fail=0)

    def readings(bs):
        sea = np.flatnonzero((mate >= 0) & (eps > 0))
        rows = np.round(p[sea] / DP).astype(int)
        out = np.zeros((len(bs), NQ))
        ctl = np.zeros((len(bs), NQ))
        fm = np.zeros(len(bs))
        nch = np.zeros(len(bs), int)
        dv = np.interp(x[bs], run.r, run.dv_eff)
        for k, b in enumerate(bs):
            r = int(np.round(p[b] / DP))
            m = sea[(rows == r) & (np.abs(x[sea] - x[b]) < YMAX + dX / 2)]
            if len(m) < 2:
                continue
            m = m[np.argsort(x[m])]
            n_norm = (len(m) / (2 * YMAX + dX) if ARGS.norm == "local" else n_row)
            out[k], fm[k], nch[k] = reading_from(x[m], x[b], dX, n_norm, dv[k],
                                                 r * DP, th[m])
            ctl[k], _, _ = reading_from(x[m], x[b], dX, n_norm, dv[k])
        return out, ctl, fm, nch

    def dark_catalysis():
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

    def find(mask, xc, pc):
        ok = mask & (np.abs(x - xc) < 1.0) & (np.abs(p - pc) < DP / 2)
        cand = np.flatnonzero(ok)
        return int(rng.choice(cand)) if len(cand) else -1

    def settle(par, q, sgn):
        """As demo_sea_lock_particles.events: recombination, else ionisation."""
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
            stats["recomb"] += 1
            return True
        j = find((mate >= 0) & (eps > 0), x[par], p[par])
        if j >= 0:
            k = mate[j]
            p[j] += t_ * xi
            p[k] -= t_ * xi
            mate[j] = mate[k] = -1
            for b in (j, k):                    # fresh integrators: hidden phases
                C[b] = rng.uniform(0, 1, NQ)
                init[b] = False
            stats["ion"] += 1
            return True
        stats["fail"] += 1
        return False

    packet_ = np.zeros(N, bool)
    packet_[2 * n_pairs:] = True
    C = rng.uniform(0, 1, (N, NQ))              # per-body channel integrators
    LV = np.zeros((N, NQ))                      # last firing level (hysteresis)
    init = np.zeros(N, bool)
    acc = {k: np.zeros((N, NQ)) for k in ("live", "ctl", "mesh", "win")}
    acc_fm = np.zeros(N)
    acc_cls = np.zeros(N)
    nch_all = []
    XB = np.linspace(-XACT, XACT, 17)
    NQB = 8
    tgt = np.zeros((16, NQB))
    att = np.zeros((16, NQB))
    real = np.zeros((16, NQB))
    gross = [0.0, 0]
    t0 = time.perf_counter()
    for _ in range(int(round(ARGS.t_end / DT))):
        kick = np.where(mate >= 0, 0.0, 1.0)                     # postulate (S')
        p += 0.5 * DT * F(x) * kick
        x += p / MU * DT
        p += 0.5 * DT * F(x) * kick
        xw = (x + L / 2) % L - L / 2
        th += p * (xw - x) / HBAR                                # --wrap-phase
        x[:] = xw
        th += (p ** 2 / (2 * MU) - V(x)) * DT / HBAR
        free = np.flatnonzero((mate < 0) & (np.abs(x) < XACT))
        if len(free):
            km = interp_rows(KQ, x[free])
            if events in ("none", "clock-live", "clock-live-h"):
                kl, kc, fm, nch = readings(free)
                nch_all.append(nch)
            if events == "none":
                acc["live"][free] += kl * DT
                acc["ctl"][free] += kc * DT
                acc["mesh"][free] += km * DT
                acc["win"][free] += k_window(x[free], dX) * DT
                acc_fm[free] += fm * DT
                acc_cls[free] += np.abs(np.interp(x[free], run.r, run.dv_eff)) * DT
            else:
                ib = np.clip(np.digitize(x[free], XB) - 1, 0, 15)
                np.add.at(tgt, ib, km[:, :NQB] * eps[free][:, None] * DT)
                gross[0] += np.abs(km).sum() * DT
                fires = []
                if events == "poisson":
                    g = np.abs(km)
                    for k, b in enumerate(free):
                        for _ in range(rng.poisson(g[k].sum() * DT)):
                            q = rng.choice(NQ, p=g[k] / g[k].sum())
                            fires.append((b, q, np.sign(km[k, q])))
                else:
                    rate = km if events == "clock-mesh" else kl
                    if events == "clock-live-h":
                        new = ~init[free]
                        if new.any():                            # enter the band on the
                            fb = free[new]                       # side the reading moves
                            LV[fb] = np.where(rate[new] >= 0, np.floor(C[fb]),
                                              np.ceil(C[fb]))
                            init[fb] = True
                        C[free] += rate * DT
                        gap = C[free] - LV[free]
                        dn = np.where(gap >= 1, np.floor(gap),
                                      np.where(gap <= -1, np.ceil(gap), 0.0))
                        LV[free] += dn
                    else:
                        old = np.floor(C[free])
                        C[free] += rate * DT
                        dn = np.floor(C[free]) - old
                    for k, q in zip(*np.nonzero(dn)):
                        fires += [(free[k], q, np.sign(dn[k, q]))] * int(abs(dn[k, q]))
                gross[1] += len(fires)
                for o in rng.permutation(len(fires)):
                    b, q, sgn = fires[o]
                    if mate[b] >= 0:
                        continue
                    ibb = int(np.clip(np.digitize(x[b], XB) - 1, 0, 15))
                    if q < NQB:
                        att[ibb, q] += sgn * eps[b]
                    if settle(b, q, sgn) and q < NQB:
                        real[ibb, q] += sgn * eps[b]
        dark_catalysis()

    corr = lambda a, b: float(np.corrcoef(a.ravel(), b.ravel())[0, 1])
    slope = lambda a, b: float((a * b).sum() / (b * b).sum())
    rel = lambda a, b: float(np.linalg.norm(a - b) / np.linalg.norm(b))
    out = dict(cfg)
    out["norm"] = ARGS.norm
    out["wall"] = round(time.perf_counter() - t0, 1)
    nch = np.concatenate(nch_all) if nch_all else np.zeros(1)
    out["chords"] = float(nch.mean())
    out["chords_zero"] = float(np.mean(nch == 0))
    if events == "none":
        m = packet_
        Al, Ac, Am, Aw = (acc[k][m] for k in ("live", "ctl", "mesh", "win"))
        out.update(slope=slope(Al, Am), corr=corr(Al, Am), relerr=rel(Al, Am),
                   slope_ctl=slope(Ac, Am), corr_ctl=corr(Ac, Am),
                   slope_win=slope(Al, Aw), corr_win=corr(Al, Aw), relerr_win=rel(Al, Aw),
                   slope_q1=slope(Al[:, 0], Am[:, 0]), slope_q2=slope(Al[:, 1], Am[:, 1]),
                   slope_q3=slope(Al[:, 2], Am[:, 2]),
                   fm_raw_rms=float(np.sqrt((acc_fm[m] ** 2).mean())),
                   cls_rms=float(np.sqrt((acc_cls[m] ** 2).mean())))
    else:
        fr = mate < 0                           # free bodies carry all of E
        e_right = float((eps[fr & (x > 0)]).sum())
        out.update(ion=stats["ion"], recomb=stats["recomb"], fail=stats["fail"],
                   att_slope=slope(att, tgt), att_corr=corr(att, tgt),
                   real_slope=slope(real, tgt), real_corr=corr(real, tgt),
                   gross_ratio=gross[1] / max(gross[0], 1e-300),
                   T_E=e_right / float(eps[fr].sum()))
    return out


def run_jobs(cfgs, label):
    workers = ARGS.workers or os.cpu_count() or 1
    workers = max(1, min(workers, len(cfgs)))
    print(f"   {len(cfgs)} runs on {workers} worker(s)", flush=True)
    res = []
    t0 = time.perf_counter()
    if workers == 1:
        it = map(particle_run, cfgs)
        pool = None
    else:
        pool = mp.get_context("fork").Pool(workers)
        it = pool.imap_unordered(particle_run, cfgs)
    for r in it:
        res.append(r)
        print(f"   done {len(res):3d}/{len(cfgs)}  {label} nu={r['nu']} dX={r['dX']} "
              f"kappa={r['kappa']} {r['events']:13s} seed={r['seed']}  "
              f"({r['wall']:.0f}s, elapsed {time.perf_counter() - t0:.0f}s)", flush=True)
    if pool is not None:
        pool.close()
        pool.join()
    return res


def write_csv(rows, name):
    keys = list(dict.fromkeys(k for r in rows for k in r))
    path = output_path(name)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in sorted(rows, key=lambda r: tuple(str(r.get(k)) for k in
                                                  ("events", "nu", "dX", "kappa", "seed"))):
            w.writerow(r)
    print(f"   wrote {name}", flush=True)


def summarise(rows, group, cols):
    keys = sorted({tuple(r[g] for g in group) for r in rows})
    print("   " + " ".join(f"{g:>12s}" for g in group) + " " +
          " ".join(f"{c:>17s}" for c in cols))
    table = {}
    for k in keys:
        sel = [r for r in rows if tuple(r[g] for g in group) == k]
        cells = []
        for c in cols:
            v = np.array([r[c] for r in sel], float)
            se = v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else float("nan")
            cells.append(f"{v.mean():8.3f}+-{se:6.3f}")
            table[k + (c,)] = (v.mean(), se)
        print("   " + " ".join(f"{str(x):>12s}" for x in k) + " " +
              " ".join(f"{c:>17s}" for c in cells))
    return table


def part_d():
    rule("D. Open loop: integrated live reading against the mesh kernel")
    cfgs = [dict(nu=nu, seed=11 + s, kappa=k, dX=dx, events="none")
            for nu in ARGS.nu for dx in ARGS.dX for k in ARGS.kappa
            for s in range(ARGS.seeds)]
    rows = run_jobs(cfgs, "D")
    write_csv(rows, f"sea_resonance_clock_D{ARGS.tag}.csv")
    print("\n   mean +- standard error over seeds (slope, corr, relerr: int Khat against"
          " int K along each body's path;\n   _win: against K averaged over the aperture"
          " window; _ctl: the mu = 0 control; fm_raw_rms: first moment\n   before the"
          " L2(c) projection, against cls_rms, the classical impulse)")
    summarise(rows, ("nu", "dX", "kappa"),
              ("chords", "slope", "corr", "relerr", "slope_win", "corr_ctl", "fm_raw_rms",
               "cls_rms"))
    return rows


def part_e():
    rule("E. Closed loop: events triggered by the integrators")
    cfgs = [dict(nu=nu, seed=11 + s, kappa=k, dX=ARGS.dX[0], events=ev)
            for nu in ARGS.nu_closed for k in ARGS.kappa for ev in ARGS.events
            for s in range(ARGS.seeds)]
    rows = run_jobs(cfgs, "E")
    ref = packet(run, r0=-8.0, p0=1.2, sr=2.0, sp=0.25)
    for _ in range(int(round(ARGS.t_end / DT))):
        ref = run.qle_step(ref, DT)
    t_mesh = float(ref[run.r > 0].sum() / ref.sum())
    write_csv(rows, f"sea_resonance_clock_E{ARGS.tag}.csv")
    print(f"\n   mesh QLE transmission at t = {ARGS.t_end}: T_E = {t_mesh:.4f}")
    print("   att_*: the trigger's signed firings against the QLE target in (x, q) bins;"
          " real_*: those realised;\n   gross_ratio: all firings over the QLE's gross"
          " rate sum |K_q| dt (1 = no excess)")
    summarise(rows, ("nu", "kappa", "events"),
              ("att_slope", "att_corr", "real_corr", "gross_ratio", "T_E", "ion", "fail"))
    return rows, t_mesh


# ------------------------------------------------------------------ figure
def figure(d_rows, b_res):
    import matplotlib
    import matplotlib.ticker
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    nus = sorted({r["nu"] for r in d_rows})
    if len(nus) < 2:
        return
    cols = ("#1D9E75", "#7F77DD", "#D85A30")
    ncol = 3 if b_res is not None else 2
    fig, ax = plt.subplots(1, ncol, figsize=(4.2 * ncol, 3.8))
    combos = sorted({(r["dX"], r["kappa"]) for r in d_rows})
    for i, (dx, k) in enumerate(combos[:3]):
        mu_c, se_c, mu_r, se_r = [], [], [], []
        for nu in nus:
            sel = [r for r in d_rows if (r["nu"], r["dX"], r["kappa"]) == (nu, dx, k)]
            c = np.array([1 - r["corr"] for r in sel])
            e = np.array([r["relerr"] for r in sel])
            mu_c.append(c.mean())
            se_c.append(c.std(ddof=1) / np.sqrt(len(c)) if len(c) > 1 else 0)
            mu_r.append(e.mean())
            se_r.append(e.std(ddof=1) / np.sqrt(len(e)) if len(e) > 1 else 0)
        lab = f"$\\Delta x$ = {dx:g}, $\\kappa$ = {k:g}"
        ax[0].errorbar(nus, mu_c, yerr=se_c, fmt="o-", color=cols[i], lw=2, ms=6,
                       capsize=3, label=lab)
        ax[1].errorbar(nus, mu_r, yerr=se_r, fmt="o-", color=cols[i], lw=2, ms=6,
                       capsize=3, label=lab)
    nn = np.array(nus, float)
    for a, ttl, yl in ((ax[0], "shape: 1 - corr", r"$1-\mathrm{corr}(\int\hat K,\int K)$"),
                       (ax[1], "size: relative error",
                        r"$\|\int\hat K-\int K\|/\|\int K\|$")):
        y0 = a.get_lines()[0].get_ydata()[0]
        a.loglog(nn, y0 * nn[0] / nn, ":", color="0.5", lw=1, label=r"$\propto 1/\nu$")
        a.loglog(nn, y0 * np.sqrt(nn[0] / nn), "--", color="0.5", lw=1,
                 label=r"$\propto 1/\sqrt{\nu}$")
        a.set_xscale("log")
        a.set_yscale("log")
        a.set_xticks(nus, [str(n) for n in nus])
        a.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        a.set_xlabel(r"ensemble multiplicity $\nu$")
        a.set_ylabel(yl)
        a.set_title(ttl, fontsize=9)
        a.grid(alpha=0.3, which="both")
        a.legend(fontsize=7)
    if b_res is not None:
        names = [n for n, _, _ in b_res["fano"]]
        vals = [f for _, _, f in b_res["fano"]]
        y = np.arange(len(names))
        ax[2].barh(y, vals, height=0.6, color=cols[0])
        ax[2].axvline(1.0, color="k", ls="--", lw=1)
        ax[2].text(0.98, len(names) - 0.45, "Poisson ", fontsize=7, ha="right",
                   va="bottom")
        ax[2].set_yticks(y, names, fontsize=7)
        ax[2].invert_yaxis()
        ax[2].set_xlim(0, 1.15)
        ax[2].set_xlabel("Fano factor of counts in unit windows")
        ax[2].set_title("deterministic clocks are sub-Poissonian", fontsize=9)
        for yi, v in zip(y, vals):
            ax[2].text(v + 0.02, yi, f"{v:.2f}", va="center", fontsize=7)
    fig.suptitle("The live sea reading against the precomputed kernel, open loop (S')",
                 y=1.02, fontsize=10)
    fig.tight_layout()
    name = f"sea_resonance_clock{ARGS.tag}.png"
    fig.savefig(output_path(name), dpi=140, bbox_inches="tight")
    dp_ = docs_path(name)
    if dp_:
        fig.savefig(dp_, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"\nwrote {name}")


def main():
    t0 = time.perf_counter()
    parts = ARGS.parts.upper()
    b_res = d_rows = None
    if "A" in parts:
        part_a()
    if "B" in parts:
        b_res = part_b()
    if "C" in parts:
        part_c()
    if "D" in parts:
        d_rows = part_d()
    if "E" in parts:
        part_e()
    if d_rows:
        figure(d_rows, b_res)
    print(f"\nwall {time.perf_counter() - t0:.0f} s")


if __name__ == "__main__":
    main()
