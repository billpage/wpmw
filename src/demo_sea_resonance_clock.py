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
            shadow        clock-mesh acts; the live reading is integrated alongside
                          along every free body's path but never fires
          and realised as in demo_sea_lock_particles.py.  Reported: the
          trigger's signed firings against the QLE target sum eps K_q dt in
          (x, q) bins, gross firings against the QLE's gross rate, the
          transmission T_E against the mesh QLE, and event counts.
          Open-loop and shadow runs also report the sea's structure seen from
          each reader's own row: the profile of sea members at offset s from
          the reader against the row's mean density, and the separations of
          the reading's chords against a uniform aperture (g_hole, g_aper,
          g_chord; histograms pooled over seeds in
          sea_resonance_clock_G<part><tag>.csv).

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
ap.add_argument("--sea-depth", type=float, nargs="+", default=[1.0],
                help="Parts D, E: the sea's density in units of B = 1/(pi hbar), at fixed "
                     "packet (beta; the ensemble multiplicity nu scales both).  The "
                     "reading is normalised by the sea actually present, so beta changes "
                     "its sampling and the ratio of sea to event traffic, not its mean")
ap.add_argument("--relock-w", type=float, nargs="+", default=[0.0],
                help="Part E: step 22 section 6 re-locking window W (0 = off).  On "
                     "recombination both members take the circular mean of the "
                     "destination row's sea clocks within W, transported to the new "
                     "pair; a body kinked out of an ionised pair likewise")
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
S_EDGE = np.linspace(-3 * YMAX, 3 * YMAX, 49)   # structure diagnostic: reader offsets
S_MID = 0.5 * (S_EDGE[1:] + S_EDGE[:-1])
D_EDGE = np.linspace(0.0, 2 * YMAX, 17)         # chord separations


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
    print("   slope of the mean reading on K over all q, pooled over x = +-1, +-2, +-3"
          " (dX = 1):\n   normalised by the expected density, and by the count N in"
          " the aperture (--norm local)")
    xs_c = np.array([-3.0, -2.0, -1.0, 1.0, 2.0, 3.0])
    for nu in (32, 64):
        n = nu * B * DP
        num, num_l, den, inv_n = 0.0, 0.0, 0.0, []
        for xb in xs_c:
            K = interp_rows(KQ, np.array([xb]))[0]
            dv = np.interp(xb, run.r, run.dv_eff)
            acc, acc_l = np.zeros(NQ), np.zeros(NQ)
            for _ in range(200):
                cnt = rng.poisson(n * (2 * YMAX + 1.0))
                xs = np.sort(rng.uniform(xb - YMAX - 0.5, xb + YMAX + 0.5, cnt))
                acc += reading_from(xs, xb, 1.0, n, dv)[0]
                acc_l += reading_from(xs, xb, 1.0, cnt / (2 * YMAX + 1.0), dv)[0]
                inv_n.append(1.0 / max(cnt, 1))
            num += acc @ K / 200
            num_l += acc_l @ K / 200
            den += K @ K
        print(f"   nu = {nu:3d}:  expected {num / den:.3f}   local {num_l / den:.3f}"
              f"   (1 - <1/N> = {1 - np.mean(inv_n):.3f})", flush=True)


# ------------------------------------------------------------------ Parts D, E
def particle_run(cfg):
    """One run of the particle model.
    cfg: dict(nu, seed, kappa, dX, events[, relock_w, sea_depth])."""
    nu, seed, kappa, dX, events = (cfg["nu"], cfg["seed"], cfg["kappa"], cfg["dX"],
                                   cfg["events"])
    relock_w = cfg.get("relock_w", 0.0)
    beta = cfg.get("sea_depth", 1.0)
    rng = np.random.default_rng(seed)
    n_pairs = int(round(beta * nu * B * L * 2 * PMAX))
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
    n_row = beta * nu * B * DP
    stats = dict(recomb=0, ion=0, fail=0, relock=0)
    made = np.zeros(N, bool)                    # members of pairs formed by an event

    def reference(excl, xj, pj):
        """Step 22 section 6: circular mean of the destination row's sea clocks
        within relock_w of xj, each carried to xj by its own momentum."""
        m = (mate >= 0) & (np.abs(p - pj) < DP / 2) & (np.abs(x - xj) < relock_w)
        m[list(excl)] = False
        if not m.any():
            return None
        stats["relock"] += 1
        return float(np.angle(np.exp(1j * (th[m] + p[m] * (xj - x[m]) / HBAR)).sum()))

    def aperture_diag(bs):
        """Over the readers' apertures: fraction of chords with an event-made end,
        and the coherence |<exp(i mu_ref)>| of the chords' readings."""
        sea = np.flatnonzero((mate >= 0) & (eps > 0))
        rows = np.round(p[sea] / DP).astype(int)
        z, n, nm = 0j, 0, 0
        for b in bs:
            r = int(np.round(p[b] / DP))
            m = sea[(rows == r) & (np.abs(x[sea] - x[b]) < YMAX + dX / 2)]
            if len(m) < 2:
                continue
            m = m[np.argsort(x[m])]
            I, J = np.triu_indices(len(m), 1)
            d = x[m[J]] - x[m[I]]
            sel = (d <= 2 * YMAX) & (np.abs(0.5 * (x[m[I]] + x[m[J]]) - x[b]) < dX / 2)
            if not sel.any():
                continue
            i, j, d = m[I[sel]], m[J[sel]], d[sel]
            z += np.exp(1j * (th[i] - th[j] + r * DP * d / HBAR)).sum()
            n += len(d)
            nm += int((made[i] | made[j]).sum())
        return z, n, nm

    def structure_diag(bs):
        """Sea structure seen from each reader's own row: counts of sea members at
        offset s = x - x_b (against the row's mean density over the box), and
        separations d of the reading's chords (against N(N-1) dX dd / Lw^2,
        the count for N members placed uniformly in the aperture of length Lw)."""
        sea = np.flatnonzero((mate >= 0) & (eps > 0))
        rows = np.round(p[sea] / DP).astype(int)
        rid, rcnt = np.unique(rows, return_counts=True)
        dens = dict(zip(rid.tolist(), (rcnt / L).tolist()))
        lw = 2 * YMAX + dX
        for b in bs:
            r = int(np.round(p[b] / DP))
            if r not in dens:
                continue
            cls = 0 if packet_[b] else 1
            m = sea[rows == r]
            s_ = (x[m] - x[b] + L / 2) % L - L / 2
            s_ = s_[np.abs(s_) < S_EDGE[-1]]
            G_prof[cls, 0] += np.histogram(s_, S_EDGE)[0]
            G_prof[cls, 1] += dens[r] * np.diff(S_EDGE)
            a = np.sort(s_[np.abs(s_) < lw / 2])
            n_ = len(a)
            if n_ < 2:
                continue
            I, J = np.triu_indices(n_, 1)
            d = a[J] - a[I]
            sel = (d <= 2 * YMAX) & (np.abs(0.5 * (a[I] + a[J])) < dX / 2)
            G_chord[cls, 0] += np.histogram(d[sel], D_EDGE)[0]
            G_chord[cls, 1] += n_ * (n_ - 1) * dX * np.diff(D_EDGE) / lw ** 2

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
            made[jn] = made[jp] = True
            if relock_w > 0:
                ref = reference({jn, jp}, xm, pm)
                if ref is not None:            # dark AND on the row's lock
                    th[jn] = th[jp] = ref
            stats["recomb"] += 1
            return True
        j = find((mate >= 0) & (eps > 0), x[par], p[par])
        if j >= 0:
            k = mate[j]
            p[j] += t_ * xi
            p[k] -= t_ * xi
            mate[j] = mate[k] = -1
            made[j] = made[k] = False
            for b in (j, k):                    # fresh integrators: hidden phases
                C[b] = rng.uniform(0, 1, NQ)
                init[b] = False
                if relock_w > 0:                # the kinked body takes its new row's lock
                    ref = reference({j, k}, x[b], p[b])
                    if ref is not None:
                        th[b] = ref
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
    diag = [0j, 0, 0]
    G_prof = np.zeros((2, 2, len(S_EDGE) - 1))   # [packet | event-born][obs | expected]
    G_chord = np.zeros((2, 2, len(D_EDGE) - 1))
    t0 = time.perf_counter()
    for step in range(int(round(ARGS.t_end / DT))):
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
            if events in ("none", "clock-live", "clock-live-h", "shadow"):
                kl, kc, fm, nch = readings(free)
                nch_all.append(nch)
            if events in ("none", "shadow") and step % 10 == 0:
                rd = free if len(free) <= 40 else rng.choice(free, 40, replace=False)
                structure_diag(rd)
            if events == "shadow":              # read, but let the mesh clock act
                acc["live"][free] += kl * DT
                acc["ctl"][free] += kc * DT
                acc["mesh"][free] += km * DT
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
                    rate = km if events in ("clock-mesh", "shadow") else kl
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
                if step % 10 == 0:
                    near = free[np.abs(x[free]) < 3.0]
                    if len(near) > 40:
                        near = rng.choice(near, 40, replace=False)
                    z_, n_, nm_ = aperture_diag(near)
                    diag[0] += z_
                    diag[1] += n_
                    diag[2] += nm_
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
    out["_G"] = (G_prof, G_chord)
    if events in ("none", "shadow"):
        pr = G_prof.sum(axis=0)
        ch = G_chord.sum(axis=0)
        inner = np.abs(S_MID) < 1.0
        out.update(g_hole=float(pr[0, inner].sum() / pr[1, inner].sum()),
                   g_aper=float(pr[0, np.abs(S_MID) < YMAX].sum()
                                / pr[1, np.abs(S_MID) < YMAX].sum()),
                   g_chord=float(ch[0].sum() / max(ch[1].sum(), 1e-300)))
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
                   T_E=e_right / float(eps[fr].sum()),
                   relock_w=relock_w, relocks=stats["relock"],
                   made_frac=diag[2] / max(diag[1], 1),
                   ap_coh=abs(diag[0]) / max(diag[1], 1))
        if events == "shadow":                  # every body that was ever read
            m = np.abs(acc["mesh"]).sum(axis=1) > 0
            Al, Ac, Am = acc["live"][m], acc["ctl"][m], acc["mesh"][m]
            out.update(sh_slope=slope(Al, Am), sh_corr=corr(Al, Am), sh_relerr=rel(Al, Am),
                       sh_slope_ctl=slope(Ac, Am), sh_corr_ctl=corr(Ac, Am),
                       sh_bodies=int(m.sum()))
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
    keys = list(dict.fromkeys(k for r in rows for k in r if not k.startswith("_")))
    path = output_path(name)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in sorted(rows, key=lambda r: tuple(str(r.get(k)) for k in
                                                  ("events", "nu", "dX", "kappa", "seed"))):
            w.writerow({k: r.get(k) for k in keys})
    print(f"   wrote {name}", flush=True)


def control_from_chords(g_d):
    """The mu = 0 reading a chord distribution g(d) (obs / uniform, on D_EDGE)
    implies, as a slope against the uniform sea's, from Theorem L1's integrand
    at eight flank positions with chord midpoints at the reader."""
    dd = np.linspace(0.0, 2 * YMAX, 1601)[1:]
    gi = np.asarray(g_d)[np.clip((dd / (2 * YMAX) * len(g_d)).astype(int), 0, len(g_d) - 1)]
    w = np.cos(np.pi * dd / (4 * YMAX)) ** 2
    proj = lambda a: a - XI * (XI @ a) / (XI @ XI)
    num = den = 0.0
    for xb in (-3.0, -2.0, -1.0, -0.5, 0.5, 1.0, 2.0, 3.0):
        dv = np.interp(xb, run.r, run.dv_eff)
        f = ((V(xb + dd / 2) - V(xb - dd / 2) - dd * dv) * w)[:, None] * np.sin(
            XI[None, :] * dd[:, None] / HBAR)
        k0, kg = proj(f.sum(0)), proj((gi[:, None] * f).sum(0))
        num += kg @ k0
        den += k0 @ k0
    return num / den


def structure_report(rows, part):
    """Pool the sea-structure histograms over seeds; print and write them."""
    rows = [r for r in rows if r["events"] in ("none", "shadow")]
    if not rows:
        return
    group = ("nu", "sea_depth", "kappa", "relock_w", "events")
    keys = sorted({tuple(r.get(g, 0.0) for g in group) for r in rows})
    out = []
    print("\n   sea structure seen from each reader's row (obs / expected): profile at"
          " offset s from the reader\n   against the row's mean density; chord"
          " separations d against N(N-1) dX dd / Lw^2 (uniform aperture)")
    for k in keys:
        sel = [r for r in rows if tuple(r.get(g, 0.0) for g in group) == k]
        pr = sum(r["_G"][0] for r in sel)
        ch = sum(r["_G"][1] for r in sel)
        for cls, name in ((0, "packet"), (1, "event-born")):
            if pr[cls, 1].sum() == 0:
                continue
            g_s = pr[cls, 0] / np.maximum(pr[cls, 1], 1e-300)
            g_d = ch[cls, 0] / np.maximum(ch[cls, 1], 1e-300)
            for kind, mid, obs, exp_ in (("profile", S_MID, pr[cls, 0], pr[cls, 1]),
                                         ("chord", 0.5 * (D_EDGE[1:] + D_EDGE[:-1]),
                                          ch[cls, 0], ch[cls, 1])):
                for c, o, e in zip(mid, obs, exp_):
                    out.append(dict(zip(group, k), readers=name, kind=kind,
                                    bin=round(float(c), 4), obs=float(o), exp=float(e)))
            prof = "  ".join(f"{v:5.3f}" for v in g_s[12:36:2])
            chrd = "  ".join(f"{v:5.3f}" for v in g_d[::2])
            print(f"   {str(k):34s} {name:10s}\n      profile s = -Y..Y: {prof}"
                  f"\n      chords  d = 0..2Y: {chrd}", flush=True)
        g_all = ch[:, 0].sum(axis=0) / np.maximum(ch[:, 1].sum(axis=0), 1e-300)
        print(f"      mu = 0 control implied by all readers' chords: "
              f"{control_from_chords(g_all):.3f} of a uniform sea's", flush=True)
    write_csv_plain(out, f"sea_resonance_clock_G{part}{ARGS.tag}.csv")


def write_csv_plain(rows, name):
    with open(output_path(name), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
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
    cfgs = [dict(nu=nu, seed=11 + s, kappa=k, dX=dx, events="none", sea_depth=b)
            for nu in ARGS.nu for b in ARGS.sea_depth for dx in ARGS.dX
            for k in ARGS.kappa for s in range(ARGS.seeds)]
    rows = run_jobs(cfgs, "D")
    write_csv(rows, f"sea_resonance_clock_D{ARGS.tag}.csv")
    print("\n   mean +- standard error over seeds (slope, corr, relerr: int Khat against"
          " int K along each body's path;\n   _win: against K averaged over the aperture"
          " window; _ctl: the mu = 0 control; fm_raw_rms: first moment\n   before the"
          " L2(c) projection, against cls_rms, the classical impulse)")
    summarise(rows, ("nu", "sea_depth", "dX", "kappa"),
              ("chords", "slope", "corr", "relerr", "slope_win", "corr_ctl", "fm_raw_rms",
               "cls_rms"))
    summarise(rows, ("nu", "sea_depth", "dX", "kappa"),
              ("slope_ctl", "g_hole", "g_aper", "g_chord"))
    structure_report(rows, "D")
    return rows


def part_e():
    rule("E. Closed loop: events triggered by the integrators")
    cfgs = [dict(nu=nu, seed=11 + s, kappa=k, dX=ARGS.dX[0], events=ev, relock_w=w,
                 sea_depth=b)
            for nu in ARGS.nu_closed for b in ARGS.sea_depth for k in ARGS.kappa
            for ev in ARGS.events for w in ARGS.relock_w for s in range(ARGS.seeds)]
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
    summarise(rows, ("nu", "sea_depth", "kappa", "relock_w", "events"),
              ("att_slope", "att_corr", "real_corr", "gross_ratio", "made_frac", "ap_coh",
               "ion", "recomb", "fail"))
    sh = [r for r in rows if r["events"] == "shadow"]
    if sh:
        print("\n   shadow runs: the live reading integrated along every free body's path"
              " while the mesh clock acts\n   (sh_*: int Khat against int K; _ctl: the"
              " mu = 0 control)")
        summarise(sh, ("nu", "sea_depth", "kappa", "relock_w"),
                  ("sh_slope", "sh_corr", "sh_relerr", "sh_slope_ctl", "sh_corr_ctl",
                   "sh_bodies"))
        summarise(sh, ("nu", "sea_depth", "kappa", "relock_w"), ("g_hole", "g_aper", "g_chord"))
        structure_report(sh, "E")
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
