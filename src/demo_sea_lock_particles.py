"""The sea as a phase reference: particle-level probe -- companion to step 22.

Measures sections 5 and 6 of ``docs/analysis/sea_phase_reference.md``:
whether partner-phase coherence emerges or persists among world-particles,
and whether event re-locking and dark catalysis maintain it.

Model (Eckart barrier, V0 = a = 1, mu = 1, hbar = 1).  Bodies carry
(x, p, theta, species, mate); aligned sea pairs are two bodies with mate
indices, free bodies have mate = -1.

  streaming   velocity Verlet under -V'(x); theta += (p^2/2m - V) dt / hbar
  events      each free body is a parent; channel q fires at |K_q(x)| from
              the compensated kernel, t = sign(K_q) s_parent.  Catalysed
              recombination if a free negaton near (x, p + t xi) and a free
              positon near (x, p - t xi) exist in the co-location bin; else
              ionisation of an aligned pair near (x, p); else the event fails
              and is counted.  Kinks keep each clock continuous (Theorem Y4);
              with --relock-w W > 0 a kinked body instead takes the circular
              mean phase of the destination row's sea clocks within W,
              transported to its position.
  dark        --dark: every aligned pair A fires at Gamma(x); its members'
              deposits cancel (Theorem L5), realised as a +-xi round trip of
              a partner pair C in A's row, which re-forms aligned on A's
              transported phase (imprint), on the mean of the two (mean), or
              on the mean with C anywhere within the reach (reach).
  measure     same-row partners within the reach: sea-clock coherence
              |<e^{i mu}>| per separation bin against its floor 1/sqrt(n),
              and corr(Kp, K) of the pair-sum kernel against the true K_q(x)
              near the barrier, with the mu = 0 control.

Usage::

    PYTHONPATH=src python3 -u src/demo_sea_lock_particles.py \
        --nu 8 --mode locked --relock-w 3 --dark reach --t-end 25 --seed 11
"""
import argparse
import time
import numpy as np
from demo_emission_and_absorption import Ledger, B, MU, HBAR

ap = argparse.ArgumentParser(description="step 22 particle probe")
ap.add_argument("--nu", type=int, default=8)
ap.add_argument("--mode", choices=("gas", "locked"), default="gas")
ap.add_argument("--dx-bin", type=float, default=2.0, help="co-location bin in x")
ap.add_argument("--relock-w", type=float, default=0.0,
                help="event re-locking window; 0 keeps phase continuity")
ap.add_argument("--dark", choices=("off", "imprint", "mean", "reach"), default="off")
ap.add_argument("--no-events", action="store_true", help="streaming (+ dark) only")
ap.add_argument("--t-end", type=float, default=0.0, help="0: one packet crossing")
ap.add_argument("--seed", type=int, default=None)
ARGS = ap.parse_args()
NU, MODE = ARGS.nu, ARGS.mode
rng = np.random.default_rng(ARGS.seed if ARGS.seed is not None else 2026 + NU)

run = Ledger(v0=1.0, a=1.0, n_r=192, r_half=24.0, n_p=64, dp=0.25)
L, DP, NP = 48.0, run.dp, run.n_p
PMAX = NP * DP / 2
YMAX = run.y_max
Q = np.arange(1, NP // 2)
XI = Q * DP
DX_BIN, RELOCK_W, DARK = ARGS.dx_bin, ARGS.relock_w, ARGS.dark
NOEVENTS, T_END = ARGS.no_events, ARGS.t_end


def V(x):
    return 1.0 / np.cosh(x) ** 2


def F(x):                                      # -V'(x)
    return 2.0 * np.tanh(x) / np.cosh(x) ** 2


def kidx(x):
    return np.clip(np.round((x - run.r[0]) / run.dr).astype(int), 0,
                   len(run.r) - 1)


# ------------------------------------------------------------ population
n_pairs = int(round(NU * B * L * 2 * PMAX))
xs = rng.uniform(-L / 2, L / 2, n_pairs)
ps = rng.uniform(-PMAX, PMAX, n_pairs)
if MODE == "locked":
    ths = np.zeros(n_pairs)                    # theta = 0 at x = 0 ...
    ths = ps * xs / HBAR                       # ... so Phi_j(x) = p_j x / hbar
else:
    ths = rng.uniform(0, 2 * np.pi, n_pairs)
r0, p0, sr, sp = -8.0, 1.2, 2.0, 0.25          # minimum uncertainty packet
rho = 6.0
n_pos, n_neg = int(round(NU * (rho + 1) / 2)), int(round(NU * (rho - 1) / 2))


def sample_packet(n):
    x = rng.normal(r0, sr, n)
    p = rng.normal(p0, sp, n)
    return x, p, p0 * (x - r0) / HBAR          # theta = S(x)/hbar


xa, pa, ta = sample_packet(n_pos)
xb, pb, tb = sample_packet(n_neg)
x = np.concatenate([xs, xs, xa, xb])
p = np.concatenate([ps, ps, pa, pb])
th = np.concatenate([ths, ths, ta, tb])
eps = np.concatenate([np.ones(n_pairs), -np.ones(n_pairs),
                      np.ones(n_pos), -np.ones(n_neg)])
mate = np.concatenate([np.arange(n_pairs, 2 * n_pairs), np.arange(n_pairs),
                       -np.ones(n_pos + n_neg, int)])
origin = np.concatenate([np.zeros(2 * n_pairs, int),
                         np.ones(n_pos + n_neg, int)])     # 0 sea, 1 packet
N = len(x)
print(f"nu={NU} mode={MODE}: {n_pairs} aligned pairs, {n_pos} positons +"
      f" {n_neg} negatons in the packet, N = {N} bodies", flush=True)


# ------------------------------------------------------------ events
stats = dict(recomb=0, ion=0, fail=0, gray=[], relock=0, noref=0, dark=0, dark_empty=0)
touched = np.zeros(N, bool)


def reference(excl, xj, pj):
    """Circular mean of the destination row's sea clocks, transported to xj."""
    if RELOCK_W <= 0:
        return None
    m = (mate >= 0) & (np.abs(p - pj) < DP / 2) & (np.abs(x - xj) < RELOCK_W)
    m[list(excl)] = False
    if not m.any():
        stats["noref"] += 1
        return None
    z = np.exp(1j * (th[m] + p[m] * (xj - x[m]) / HBAR)).sum()
    stats["relock"] += 1
    return float(np.angle(z))


def find(mask_extra, xc, pc):
    """Index of one body in the co-location bin around (xc, pc), or -1."""
    ok = mask_extra & (np.abs(x - xc) < DX_BIN / 2) & (np.abs(p - pc) < DP / 2)
    cand = np.flatnonzero(ok)
    return int(rng.choice(cand)) if len(cand) else -1


def events(dt):
    free = np.flatnonzero(mate < 0)
    g = run.gamma_tot[kidx(x[free])]
    n_ev = rng.poisson(g * dt)
    for par, n in zip(free[n_ev > 0], n_ev[n_ev > 0]):
        for _ in range(n):
            if mate[par] >= 0:                # parent was consumed already
                break
            kq = run.k[kidx(np.array([x[par]]))[0], Q]
            if np.abs(kq).sum() == 0.0:        # parent moved onto K = 0
                break
            q = rng.choice(len(Q), p=np.abs(kq) / np.abs(kq).sum())
            t = np.sign(kq[q]) * eps[par]
            xi = XI[q]
            free_m = (mate < 0) & (np.arange(N) != par)
            j_neg = find(free_m & (eps < 0), x[par], p[par] + t * xi)
            j_pos = find(free_m & (eps > 0), x[par], p[par] - t * xi)
            if j_neg >= 0 and j_pos >= 0:
                # recombination: kinks, then co-locate at the mean (x, p)
                p[j_neg] -= t * xi
                p[j_pos] += t * xi
                xm = 0.5 * (x[j_neg] + x[j_pos])
                pm = 0.5 * (p[j_neg] + p[j_pos])
                for j in (j_neg, j_pos):       # carry extended phase to xm
                    th[j] += p[j] * (xm - x[j]) / HBAR
                    x[j], p[j] = xm, pm
                mate[j_neg], mate[j_pos] = j_pos, j_neg
                ref = reference({j_neg, j_pos}, xm, pm)
                if ref is not None:            # aligned AND on the field
                    th[j_neg] = th[j_pos] = ref
                    touched[j_neg] = touched[j_pos] = True
                mu = np.angle(np.exp(1j * (th[j_pos] - th[j_neg])))
                stats["gray"].append(abs(mu))
                stats["recomb"] += 1
                continue
            j = find((mate >= 0) & (eps > 0), x[par], p[par])
            if j >= 0:                          # ionisation of an aligned pair
                k = mate[j]
                p[j] += t * xi
                p[k] -= t * xi
                mate[j] = mate[k] = -1
                for b in (j, k):
                    ref = reference({j, k}, x[b], p[b])
                    if ref is not None:
                        th[b] = ref
                        touched[b] = True
                stats["ion"] += 1
            else:
                stats["fail"] += 1


def dark_catalysis(dt):
    """Candidate 2: every aligned pair A fires at Gamma(x); its two members'
    deposits cancel, realised as a +-xi round trip of a co-located pair C in
    A's row, which re-forms aligned on A's phase (imprint) or on the mean."""
    if DARK == "off":
        return
    A_all = np.flatnonzero((mate >= 0) & (eps > 0))
    fire = A_all[rng.poisson(run.gamma_tot[kidx(x[A_all])] * dt) > 0]
    for a in fire:
        if mate[a] < 0:
            continue
        rng_x = 2 * YMAX if DARK == "reach" else DX_BIN / 2
        cand = np.flatnonzero((mate >= 0) & (eps > 0) & (np.abs(x - x[a]) < rng_x)
                              & (np.abs(p - p[a]) < DP / 2))
        cand = cand[cand != a]
        if not len(cand):
            stats["dark_empty"] += 1
            continue
        c = int(rng.choice(cand))
        phi_a = th[a] + p[a] * (x[c] - x[a]) / HBAR     # A's phase at C
        if DARK == "imprint":
            new = phi_a
        else:
            new = float(np.angle(np.exp(1j * phi_a) + np.exp(1j * th[c])))
            th[a] = th[mate[a]] = new - p[a] * (x[c] - x[a]) / HBAR
            touched[a] = touched[mate[a]] = True
        th[c] = th[mate[c]] = new
        touched[c] = touched[mate[c]] = True
        stats["dark"] += 1


# ------------------------------------------------------------ measurement
DBINS = np.linspace(0, 2 * YMAX, 9)


def measure():
    row = np.round(p / DP).astype(int)
    order = np.lexsort((x, row))
    xr, rr = x[order], row[order]
    ii, jj = [], []
    for r in np.unique(rr):
        idx = np.flatnonzero(rr == r)
        xx = xr[idx]
        hi = np.searchsorted(xx, xx + 2 * YMAX)
        cnt = hi - np.arange(len(xx)) - 1
        a = np.repeat(np.arange(len(xx)), cnt)
        b = a + 1 + (np.arange(cnt.sum()) - np.repeat(np.cumsum(cnt) - cnt, cnt))
        ii.append(order[idx[a]])
        jj.append(order[idx[b]])
    i, j = np.concatenate(ii), np.concatenate(jj)
    d = x[j] - x[i]
    xm = 0.5 * (x[i] + x[j])
    mu = (th[i] + p[i] * (xm - x[i]) / HBAR) - (th[j] + p[j] * (xm - x[j]) / HBAR)
    notmate = mate[i] != j
    out = {}
    for name, sel, weighted in (
            ("sea clocks", (mate[i] >= 0) & (mate[j] >= 0) & notmate, False),
            ("touched", touched[i] & touched[j] & notmate, False),
            ("free, eps", (mate[i] < 0) & (mate[j] < 0), True),
            ("all, eps", notmate, True)):
        z = np.exp(1j * mu) * (eps[i] * eps[j] if weighted else 1.0)
        b = np.digitize(d[sel], DBINS) - 1
        n = np.bincount(b, minlength=8)[:8]
        s = np.abs(np.bincount(b, z[sel].real, minlength=8)[:8]
                   + 1j * np.bincount(b, z[sel].imag, minlength=8)[:8])
        out[name] = (s / np.maximum(n, 1), 1 / np.sqrt(np.maximum(n, 1)), n)
    # pair-sum kernel near the barrier, with mu and with mu = 0
    near = np.abs(xm) < 3.0
    y = 0.5 * d[near]
    c = np.clip(((xm[near] + 3.0) / 0.5).astype(int), 0, 11)
    dv = run.dv_eff[kidx(xm[near])]
    w = np.cos(np.pi * y / (2 * YMAX)) ** 2
    ures = (V(xm[near] + y) - V(xm[near] - y) - 2 * y * dv) * w / HBAR
    ee = (eps[i] * eps[j])[near]
    nm = notmate[near]
    kp = np.zeros((12, len(Q)))
    kpe = np.zeros((12, len(Q)))
    kp0 = np.zeros((12, len(Q)))
    for q, xq in enumerate(XI):
        sm = np.sin(xq * d[near] / HBAR + mu[near])
        kp[:, q] = np.bincount(c, -ures * sm * nm, minlength=12)[:12]
        kpe[:, q] = np.bincount(c, -ures * sm * ee * nm, minlength=12)[:12]
        kp0[:, q] = np.bincount(c, -ures * np.sin(xq * d[near] / HBAR) * nm,
                                minlength=12)[:12]
    xc = -3.0 + 0.5 * (np.arange(12) + 0.5)
    kt = run.k[kidx(xc)][:, Q]
    corr = lambda a, b: float(np.corrcoef(a.ravel(), b.ravel())[0, 1])
    return out, corr(kp, kt), corr(kpe, kt), corr(kp0, kt), len(i)


# ------------------------------------------------------------ run
dt = 0.02
nstep = int(round((T_END or 16.0 / (p0 / MU)) / dt))
snaps = sorted(set([0, 5, 15, 50] + list(range(0, nstep + 1, nstep // 6))))
t0 = time.perf_counter()
print(f"{'t':>6} {'pairs':>8} {'corr(Kp,K)':>11} {'eps-weighted':>12} {'mu=0 control':>12}  "
      f"coherence C(d)/floor, same-row, d-bins of {DBINS[1]:.2f}:", flush=True)
for step in range(nstep + 1):
    if step in snaps:
        out, ck, cke, ck0, npairs = measure()
        line = f"{step*dt:6.2f} {npairs:8d} {ck:11.3f} {cke:12.3f} {ck0:12.3f}"
        near_sea = (mate >= 0) & (np.abs(x) < 4.0)
        line += f"   touched: {touched[near_sea].mean() if near_sea.any() else 0:.2f} of sea near barrier"
        print(line, flush=True)
        for name, (cc, fl, n) in out.items():
            ratio = " ".join(f"{a/b:5.1f}" for a, b in zip(cc, fl))
            print(f"         {name:9s} C/floor: {ratio}   (C in bin 1:"
                  f" {cc[0]:.3f}, n={n[0]})", flush=True)
    if step == nstep:
        break
    p += 0.5 * dt * F(x)
    x += p / MU * dt
    p += 0.5 * dt * F(x)
    x[:] = (x + L / 2) % L - L / 2
    th += (p ** 2 / (2 * MU) - V(x)) * dt / HBAR
    if not NOEVENTS:
        events(dt)
    dark_catalysis(dt)
g = np.array(stats["gray"])
print(f"\nrelock W={RELOCK_W}: {stats['relock']} re-locked, {stats['noref']} without reference;"
      f"  dark catalysis ({DARK}): {stats['dark']} events, {stats['dark_empty']} with no co-located pair")
print(f"events: recombination {stats['recomb']}, ionisation {stats['ion']},"
      f" failed {stats['fail']};  created pairs |mu| mean"
      f" {g.mean() if len(g) else float('nan'):.3f} (uniform would be pi/2 ="
      f" {np.pi/2:.3f});  wall {time.perf_counter()-t0:.0f} s")
