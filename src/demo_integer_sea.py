"""The integer sea at the physical density -- step 24, open item V-SP2.

Companion to ``docs/analysis/sea_depletion.md``.  Proposition V4 treated an
integer sea in mean field: the parent's event aperture (Q9(c), area
dp * 2 y_max = h/2) holds Poisson(n) pairs, an emission finds it empty with
probability exp(-n), and the expected E then departs by C exp(-beta).  This
demo drops the mean field.  Bodies and sea pairs are points; the world is
one sample; the expectation is an average over many independent worlds.

Model (one world, nu worlds per unit of W; the physical density is nu = 1)
  bodies     points (x, p, s), s = +1 positon, -1 negaton, weight 1/nu;
             stream by velocity Verlet under the compensated force
             -V'_eff(x) of the mesh (step 16's kernel, ``SeaLedger``).
  sea        aligned pairs as points (x, p), density nu beta B, placed
             uniformly at random (step 23 s1: a random configuration);
             under (S) they stream like bodies, under (S') inertially.
  events     each body fires at the per-parent rate sum_{q>=1} |K_q(x)|,
             channel q with probability |K_q| / sum; t = sign(K_q) s.  The
             event deposits +1 at p + t xi_q and -1 at p - t xi_q:
               absorptive  if a free negaton lies in the aperture about
                           (x, p + t xi) and a free positon about
                           (x, p - t xi): both bind into a pair at (x, p);
               emissive    otherwise: a pair in the aperture about (x, p)
                           is ionised into a positon at (x, p + t xi) and a
                           negaton at (x, p - t xi);
               blocked     (floor) if that aperture is empty: nothing
                           happens.  Without the floor the emission takes
                           place anyway and the sea goes into debt.
             The aperture is |dx| < y_max, |dp| < dp/2: area h/2, so it
             holds nu beta pairs on average (Q9(c)).
  contact    garbage collection, Proposition Q8: a free positon and a free
             negaton in one mesh cell bind into a pair at their mean.

Paired worlds.  The floored world and the unfloored one are run from the
same seed.  The children of an emission are placed at the parent, which
step 24's mesh ledgers do too, so the bodies never depend on which pair
was ionised, and the pair is chosen with a separate random stream.  The
two worlds are then identical event for event until the first block, and
the paired difference of any linear observable has a small variance.

Observables are linear in E (an average over worlds is only meaningful for
those): the transmission T = sum s [x > 0] / nu, <p>, <H>, and E on a
coarse grid.

Run as::

    WPMW_OUTPUT=... PYTHONPATH=src python3 -u src/demo_integer_sea.py \\
        --case eckart --sea-force S' --betas 1 2 4 8 --seeds 2000

Writes ``integer_sea_<case>_<motion><tag>.csv`` (one row per beta: means,
paired differences and their standard errors) and an ``.npz`` of the
coarse histograms through ``output_path``.  ``--summarise CSV...`` combines
finished sweeps (several seed blocks of one case add up, weighted by
worlds) and draws ``integer_sea_blocking.png``.
"""

from __future__ import annotations

import argparse
import csv
import multiprocessing as mp
import os
import time

import numpy as np

import demo_sea_depletion as sd
from wpmwlib.wpmw_utils import output_path

B = sd.B
HBAR = sd.HBAR

XB = np.linspace(-20.0, 20.0, 41)        # coarse histogram of E
NOCC = 24                                # aperture occupancy histogram
PB = np.linspace(-4.0, 4.0, 17)


# ----------------------------------------------------------------------
class Geometry:
    """Kernel, force and aperture taken from the mesh (``SeaLedger``)."""

    def __init__(self, case):
        Vf = sd.eckart() if case == "eckart" else sd.well()
        m = sd.SeaLedger(Vf)
        self.Vf = Vf
        self.r0, self.dr = float(m.r[0]), float(m.dr)
        self.n_r, self.n_p, self.dp = len(m.r), m.n_p, m.dp
        self.L = self.n_r * self.dr
        self.PW = self.n_p * self.dp                 # momentum period
        self.YMAX = m.y_max
        self.r = m.r
        self.dv = m.dv_eff
        Q = np.arange(1, self.n_p // 2)
        self.XI = Q * self.dp
        kq = m.k[:, Q]                               # signed rate, q >= 1
        self.GAM = np.abs(kq).sum(axis=1)            # per-parent event rate
        cum = np.cumsum(np.abs(kq), axis=1)
        self.CUM = cum / np.where(cum[:, -1:] > 0, cum[:, -1:], 1.0)
        self.SGN = np.sign(kq)
        self.KABS = np.abs(kq)                       # |K_q(x)|, q = 1 ...

    def ix(self, x):
        return np.clip(np.round((x - self.r0) / self.dr).astype(np.int64),
                       0, self.n_r - 1)

    def force(self, x):
        return -np.interp(x, self.r, self.dv, period=self.L)

    def wrap_x(self, x):
        return (x - self.r0 + 0.5 * self.dr) % self.L + self.r0 - 0.5 * self.dr

    def wrap_p(self, p):
        return (p + 0.5 * self.PW) % self.PW - 0.5 * self.PW

    def dx(self, a, b):
        return (a - b + 0.5 * self.L) % self.L - 0.5 * self.L

    def dpw(self, a, b):
        return (a - b + 0.5 * self.PW) % self.PW - 0.5 * self.PW


# ----------------------------------------------------------------------
# Pair collisions between point pairs (step 24 Definition (C), V-SP10)
# ----------------------------------------------------------------------
def window_pairs(g, xs):
    """All unordered pairs (i, j) of sea pairs with |x_i - x_j| < y_max on
    the periodic box: one aperture's reach in position."""
    n = xs.size
    if n < 2:
        e = np.empty(0, np.int64)
        return e, e
    order = np.argsort(xs)
    xo = xs[order]
    ext = np.concatenate([xo, xo + g.L])
    hi = np.searchsorted(ext, xo + g.YMAX, side="left")
    hi = np.minimum(hi, np.arange(n) + n)             # never pair with self
    cnt = hi - np.arange(n) - 1
    if cnt.sum() == 0:
        e = np.empty(0, np.int64)
        return e, e
    i = np.repeat(np.arange(n), cnt)
    start = np.repeat(np.cumsum(cnt) - cnt, cnt)
    j = (np.arange(i.size) - start) + i + 1
    return order[i], order[j % n]


def collide(g, xs, ps, dt, mode, scale, rng):
    """Definition (C) between point pairs over one step, in place on ps.

    Two pairs within one aperture's reach (|dx| < y_max) collide.  Rates are
    chosen so that mass action, J = r_q (s_n^2 - s_{n-q} s_{n+q}), holds in
    mean: an aperture's position window at density s on a row of width dp
    holds s / B partners (Q9(c)), so
        defocus  (both on row n -> rows n - q, n + q): 2 r_q B per unordered
                 pair, either member going down with probability 1/2;
        focus    (rows n - q and n + q -> both on n):    r_q B per pair.
    kernel : r_q B = scale |K_q(x_mid)| / 2, every channel;
    uniform: r B = scale, neighbouring rows only (q = 1).
    Each transition and its reverse have equal rates, so a Poisson sea is
    exactly invariant (detailed balance with respect to the Poisson point
    process).  Momentum moves by +-q dp: pair number and pair momentum are
    conserved (up to the periodic momentum box).  The step is a tau-leap,
    with substeps so that no candidate fires with probability above 0.1,
    and a pair moves at most once per substep.  Returns the number of
    defocus and focus events."""
    nd = nf = 0
    if xs.size < 2:
        return nd, nf
    nq = g.KABS.shape[1]
    I, J = window_pairs(g, xs)
    if I.size == 0:
        return nd, nf
    rmax = (scale * g.GAM.max() if mode == "kernel" else 2.0 * scale)
    nsub = max(1, int(np.ceil(rmax * dt / 0.1)))
    h = dt / nsub
    if mode == "kernel":                  # positions are fixed in a step
        ir_all = g.ix(g.wrap_x(xs[I] + 0.5 * g.dx(xs[J], xs[I])))
    I0, J0 = I, J
    for _ in range(nsub):
        rows = np.floor(ps / g.dp).astype(np.int64) % g.n_p
        dr = (rows[J0] - rows[I0] + g.n_p // 2) % g.n_p - g.n_p // 2
        # only same-row pairs and pairs an even number of rows apart can
        # collide; drop the rest before any floating-point work
        keep = (dr % 2 == 0) if mode == "kernel" else (np.abs(dr) % 2 == 0) & (np.abs(dr) <= 2)
        I, J, dr = I0[keep], J0[keep], dr[keep]
        if mode == "kernel":
            ir = ir_all[keep]
            dmask = dr == 0
            fmask = (dr != 0) & (dr % 2 == 0) & (np.abs(dr) // 2 <= nq)
            rate = np.zeros(I.size)
            rate[dmask] = scale * g.GAM[ir[dmask]]
            qf = np.abs(dr[fmask]) // 2
            rate[fmask] = 0.5 * scale * g.KABS[ir[fmask], qf - 1]
        else:
            dmask = dr == 0
            fmask = np.abs(dr) == 2
            rate = np.where(dmask, 2.0 * scale, np.where(fmask, scale, 0.0))
        fire = np.flatnonzero(rng.random(I.size) < -np.expm1(-rate * h))
        if fire.size == 0:
            continue
        fire = fire[rng.permutation(fire.size)]
        used = np.zeros(xs.size, bool)
        for k in fire:
            a, b = I[k], J[k]
            if used[a] or used[b]:
                continue
            used[a] = used[b] = True
            if dr[k] == 0:
                if mode == "kernel":
                    ir_ = ir[k]
                    q = int(np.searchsorted(g.CUM[ir_], rng.random())) + 1
                    q = min(q, nq)
                else:
                    q = 1
                xi = q * g.dp
                if rng.random() < 0.5:
                    a, b = b, a
                ps[a] = g.wrap_p(ps[a] - xi)
                ps[b] = g.wrap_p(ps[b] + xi)
                nd += 1
            else:
                q = abs(int(dr[k])) // 2
                xi = q * g.dp
                lo, up = (a, b) if dr[k] > 0 else (b, a)
                ps[lo] = g.wrap_p(ps[lo] + xi)
                ps[up] = g.wrap_p(ps[up] - xi)
                nf += 1
    return nd, nf


def poisson_check(case, beta, t_max, dt, mode, scale, seeds, probes=200):
    """Collisions alone keep a Poisson sea Poisson: the occupancy of random
    apertures (|dx| < y_max, |dp| < dp/2) after t_max against Poisson of
    the same mean, and the start's.  No bodies, no events."""
    g = Geometry(case)
    occ0, occ1 = [], []
    nev = 0
    for sd_ in seeds:
        rng = np.random.default_rng((sd_, 31337))
        n = rng.poisson(beta * B * g.L * g.PW)
        xs = g.r0 - 0.5 * g.dr + g.L * rng.random(n)
        ps = -0.5 * g.PW + g.PW * rng.random(n)
        p0 = ps.copy()

        def probe(ps_):
            xc = rng.uniform(-8.0, 8.0, probes)    # where the kernel acts
            pc = rng.uniform(-0.5 * g.PW, 0.5 * g.PW, probes)
            return [int(((np.abs(g.dx(xs, a)) < g.YMAX)
                         & (np.abs(g.dpw(ps_, b)) < 0.5 * g.dp)).sum())
                    for a, b in zip(xc, pc)]
        occ0 += probe(p0)
        for _ in range(int(round(t_max / dt))):
            d, f = collide(g, xs, ps, dt, mode, scale, rng)
            nev += d + f
        occ1 += probe(ps)
    o0, o1 = np.array(occ0), np.array(occ1)
    return dict(mean0=o0.mean(), p00=(o0 == 0).mean(), mean1=o1.mean(),
                p01=(o1 == 0).mean(), var1=o1.var(), events=nev / len(seeds))


# ----------------------------------------------------------------------
def run_world(g, seed, t_max, dt, beta, floor, sea_force, nu=1, gc=True,
              case="eckart", coll=None):
    """One world.  Returns a dict of linear observables and counts.
    ``coll`` = (mode, scale) turns on pair collisions (Definition (C)) in
    the sea, with a random stream of their own, so the bodies and the
    sea's initial configuration are those of the world without them."""
    rng = np.random.default_rng(seed)                 # the bodies' stream
    rng_s = np.random.default_rng((seed, 7919, int(round(1000 * beta))))
    rng_c = np.random.default_rng((seed, 104729, int(round(1000 * beta))))
    # initial packet: nu samples of W0, all positons (W0 > 0)
    if case == "eckart":
        x0, p0, sx = 0.0, 1.0, 1.0
    else:
        x0, p0, sx = 2.0, 0.0, 0.6
    sp_ = HBAR / (2.0 * sx)
    x = rng.normal(x0, sx, nu)
    p = rng.normal(p0, sp_, nu)
    s = np.ones(nu)
    # the sea: Poisson number of pairs, uniform in the box
    n_sea = rng_s.poisson(nu * beta * B * g.L * g.PW) if floor else 0
    xs = g.r0 - 0.5 * g.dr + g.L * rng_s.random(n_sea)
    ps = -0.5 * g.PW + g.PW * rng_s.random(n_sea)
    cnt = dict(abs=0, emi=0, blocked=0, debt=0, gc=0, nmax=nu, coll=0,
               col=0.0)
    occ = np.zeros(NOCC)              # pairs found in the aperture, per attempt

    def stream(x, p, h, feel=True):
        if feel:
            p = p + 0.5 * h * g.force(x)
            x = x + p * h
            p = p + 0.5 * h * g.force(x)
        else:
            x = x + p * h
        return g.wrap_x(x), g.wrap_p(p)

    n_steps = int(round(t_max / dt))
    for _ in range(n_steps):
        x, p = stream(x, p, 0.5 * dt)
        if n_sea:
            xs, ps = stream(xs, ps, 0.5 * dt, feel=sea_force)
        # ---- events
        lam = g.GAM[g.ix(x)] * dt
        nev = rng.poisson(lam)
        alive = np.ones(x.size, bool)
        new_x, new_p, new_s = [], [], []
        sea_add_x, sea_add_p = [], []
        for i in np.flatnonzero(nev):
            for _k in range(nev[i]):
                if not alive[i]:
                    break
                ir = g.ix(np.array([x[i]]))[0]
                q = int(np.searchsorted(g.CUM[ir], rng.random()))
                q = min(q, g.XI.size - 1)
                t = g.SGN[ir, q] * s[i]
                if t == 0:
                    continue
                xi = g.XI[q]
                near = alive & (np.abs(g.dx(x, x[i])) < g.YMAX)
                cand_n = np.flatnonzero(
                    near & (s < 0)
                    & (np.abs(g.dpw(p, p[i] + t * xi)) < 0.5 * g.dp))
                cand_p = np.flatnonzero(
                    near & (s > 0)
                    & (np.abs(g.dpw(p, p[i] - t * xi)) < 0.5 * g.dp))
                if cand_n.size and cand_p.size:
                    a = cand_n[rng.integers(cand_n.size)]
                    b = cand_p[rng.integers(cand_p.size)]
                    alive[a] = alive[b] = False
                    sea_add_x.append(x[i])
                    sea_add_p.append(p[i])
                    cnt["abs"] += 1
                    continue
                if floor:
                    if n_sea:
                        ok = np.flatnonzero(
                            (np.abs(g.dx(xs, x[i])) < g.YMAX)
                            & (np.abs(g.dpw(ps, p[i])) < 0.5 * g.dp))
                    else:
                        ok = np.empty(0, np.int64)
                    occ[min(ok.size, NOCC - 1)] += 1
                    # the column: pairs within reach in x on every row,
                    # per aperture (PW / dp apertures per column).  Pair
                    # collisions move pairs between rows at fixed x, so
                    # they can bring an aperture up to this, never above
                    if n_sea:
                        cnt["col"] += float((np.abs(g.dx(xs, x[i]))
                                             < g.YMAX).sum()) / g.n_p
                    if ok.size == 0:
                        cnt["blocked"] += 1
                        continue
                    j = ok[rng_s.integers(ok.size)]
                    xs = np.delete(xs, j)
                    ps = np.delete(ps, j)
                    n_sea -= 1
                else:
                    cnt["debt"] += 1
                new_x += [x[i], x[i]]
                new_p += [g.wrap_p(p[i] + t * xi), g.wrap_p(p[i] - t * xi)]
                new_s += [1.0, -1.0]
                cnt["emi"] += 1
        if new_x or not alive.all():
            x = np.concatenate([x[alive], np.array(new_x)])
            p = np.concatenate([p[alive], np.array(new_p)])
            s = np.concatenate([s[alive], np.array(new_s)])
        if sea_add_x and floor:
            xs = np.concatenate([xs, np.array(sea_add_x)])
            ps = np.concatenate([ps, np.array(sea_add_p)])
            n_sea = xs.size
        # ---- contact sink (garbage collection), cell by cell
        if gc and x.size > 1:
            key = g.ix(x) * g.n_p + (np.floor(p / g.dp).astype(np.int64)
                                     % g.n_p)
            order = np.lexsort((s, key))
            ks, ss = key[order], s[order]
            keep = np.ones(x.size, bool)
            bnd = np.flatnonzero(np.diff(ks)) + 1
            for lo, hi in zip(np.r_[0, bnd], np.r_[bnd, ks.size]):
                if hi - lo < 2:
                    continue
                neg = order[lo:hi][ss[lo:hi] < 0]
                pos = order[lo:hi][ss[lo:hi] > 0]
                m_ = min(neg.size, pos.size)
                if m_:
                    a, b = neg[:m_], pos[:m_]
                    keep[a] = keep[b] = False
                    cnt["gc"] += m_
                    if floor:
                        xs = np.concatenate(
                            [xs, x[a] + 0.5 * g.dx(x[b], x[a])])
                        ps = np.concatenate(
                            [ps, p[a] + 0.5 * g.dpw(p[b], p[a])])
                        n_sea = xs.size
            if not keep.all():
                x, p, s = x[keep], p[keep], s[keep]
        cnt["nmax"] = max(cnt["nmax"], x.size)
        # ---- the sea's own motion: pair collisions (amended (S'))
        if coll is not None and n_sea > 1:
            d_, f_ = collide(g, xs, ps, dt, coll[0], coll[1], rng_c)
            cnt["coll"] += d_ + f_
        x, p = stream(x, p, 0.5 * dt)
        if n_sea:
            xs, ps = stream(xs, ps, 0.5 * dt, feel=sea_force)
    V = g.Vf(x)
    h, _, _ = np.histogram2d(x, p, bins=(XB, PB), weights=s)
    return dict(T=float(s[x > 0].sum()) / nu, P=float((s * p).sum()) / nu,
                H=float((s * (0.5 * p ** 2 + V)).sum()) / nu,
                N=x.size / nu, hist=h / nu, occ=occ, **cnt)


# ----------------------------------------------------------------------
OBS = ("T", "P", "H", "N")
CNT = ("abs", "emi", "blocked", "debt", "gc", "nmax", "coll", "col")


def chunk(args):
    """A block of seeds; each seed runs the unfloored world and every
    setting (beta, collisions).  Returns sums and sums of squares, so
    nothing per seed is shipped.  ``blk`` holds per-world sums of blocked
    emissions b and attempts a (b, a, b^2, a^2, ab) for the ratio's
    standard error."""
    case, motion, settings, seeds, t_max, dt, nu, gc = args
    g = Geometry(case)
    force = motion == "S"
    nb = len(settings)
    acc = dict(n=0, blk=np.zeros((nb, 5)),
               ref=np.zeros((len(OBS), 2)),
               val=np.zeros((nb, len(OBS), 2)),
               dif=np.zeros((nb, len(OBS), 2)),
               cnt=np.zeros((nb + 1, len(CNT))),
               href=np.zeros((XB.size - 1, PB.size - 1)),
               hval=np.zeros((nb, XB.size - 1, PB.size - 1)),
               hdif2=np.zeros((nb, XB.size - 1, PB.size - 1)),
               occ=np.zeros((nb, NOCC)))
    for sd_ in seeds:
        o0 = run_world(g, sd_, t_max, dt, 1.0, False, force, nu, gc, case)
        r0 = np.array([o0[k] for k in OBS])
        acc["ref"][:, 0] += r0
        acc["ref"][:, 1] += r0 ** 2
        acc["cnt"][0] += [o0[k] for k in CNT]
        acc["href"] += o0["hist"]
        for j, (b, coll) in enumerate(settings):
            o = run_world(g, sd_, t_max, dt, b, True, force, nu, gc, case,
                          coll)
            bb, aa = o["blocked"], o["blocked"] + o["emi"]
            acc["blk"][j] += (bb, aa, bb * bb, aa * aa, bb * aa)
            r = np.array([o[k] for k in OBS])
            acc["val"][j, :, 0] += r
            acc["val"][j, :, 1] += r ** 2
            acc["dif"][j, :, 0] += r - r0
            acc["dif"][j, :, 1] += (r - r0) ** 2
            acc["cnt"][j + 1] += [o[k] for k in CNT]
            acc["hval"][j] += o["hist"]
            acc["hdif2"][j] += (o["hist"] - o0["hist"]) ** 2
            acc["occ"][j] += o["occ"]
        acc["n"] += 1
    return acc


def merge(a, b):
    if a is None:
        return b
    for k in a:
        a[k] = a[k] + b[k]
    return a


def mean_se(s1, s2, n):
    m = s1 / n
    var = np.maximum(s2 / n - m ** 2, 0.0)
    return m, np.sqrt(var / max(n - 1, 1))


# ----------------------------------------------------------------------
# Summary of finished sweeps (several CSVs of one case and motion combine)
# ----------------------------------------------------------------------
# Proposition V4's mean-field prediction for dT on the Eckart summit
# (demo_sea_depletion.py Part D, the minimal ledger, T = 8).
MF_DT = {("eckart", "S'"): {1: -0.031, 2: -0.013, 4: -0.0019, 8: -3.5e-5},
         ("eckart", "S"): {1: -0.032, 2: -0.014, 4: -0.0020, 8: -3.7e-5}}


def summarise(paths, fig_name="integer_sea_blocking.png"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from wpmwlib.wpmw_utils import docs_path

    groups = {}
    for pth in paths:
        with open(pth) as fh:
            for r in csv.DictReader(fh):
                groups.setdefault((r["case"], r["motion"]), {}).setdefault(
                    r["beta"], []).append(r)

    def comb(rs, key, se=None):
        n = np.array([float(r["worlds"]) for r in rs])
        m = float(sum(float(r[key]) * w for r, w in zip(rs, n)) / n.sum())
        if se is None:
            return m
        return m, float(np.sqrt(sum((w / n.sum()) ** 2 * float(r[se]) ** 2
                                    for r, w in zip(rs, n))))

    fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.4))
    rows = []
    for (case, mot), g in sorted(groups.items()):
        ref = g["inf"]
        nw = sum(float(r["worlds"]) for r in ref)
        print(f"\n{case} ({mot}): {nw:.0f} worlds")
        for k in ("T", "P", "H", "N"):
            m, e = comb(ref, k, k + "_se")
            print(f"  unfloored {k} = {m:.4f} ± {e:.4f}")
        print(f"  per world: abs {comb(ref, 'per_world_abs'):.2f}"
              f"  emi {comb(ref, 'per_world_emi'):.2f}"
              f"  gc {comb(ref, 'per_world_gc'):.2f}")
        bs = sorted(float(b) for b in g if b != "inf")
        bl, dts = [], []
        print(f"  {'beta':>5} {'blocked':>8} {'e^-b':>9} {'<k>':>6}"
              f" {'P0/Pois':>8} {'dT':>17} {'mean field':>10}")
        for b in bs:
            rs = g[next(k for k in g if k != "inf" and float(k) == b)]
            blk = comb(rs, "blocked_share")
            km = comb(rs, "aperture_mean")
            p0 = comb(rs, "aperture_p0")
            dT, eT = comb(rs, "dT", "dT_se")
            mf = MF_DT.get((case, mot), {}).get(int(b), float("nan"))
            print(f"  {b:5g} {blk:8.4f} {np.exp(-b):9.5f} {km:6.3f}"
                  f" {p0 / np.exp(-km):8.2f} {dT:+8.4f}±{eT:.4f} {mf:+10.4f}")
            bl.append(blk)
            dts.append((dT, eT, mf))
            rows.append([case, mot, b, nw, blk, km, p0, dT, eT, mf])
        bl = np.array(bl)
        slope, icpt = np.polyfit(bs, np.log(bl), 1)
        print(f"  fit: blocked share = {np.exp(icpt):.3f} exp({slope:.3f} beta)")
        col = {("eckart", "S'"): "#2a78d6", ("eckart", "S"): "#eb6834",
               ("well", "S'"): "#1baf7a"}.get((case, mot), "0.4")
        lab = (f"{case} ({mot}): {np.exp(icpt):.2f}"
               f" $e^{{{slope:.2f}\\beta}}$")
        ax[0].semilogy(bs, bl, "o-", color=col, label=lab)
        if (case, mot) in MF_DT:
            d = np.array(dts)
            ax[1].errorbar(np.array(bs) + (0.08 if mot == "S" else 0.0),
                           d[:, 0], yerr=2 * d[:, 1], fmt="o", capsize=3,
                           color=col, label=f"{case} ({mot}), integer, ±2 SE")
            ax[1].plot(bs, d[:, 2], "x--", color=col, alpha=0.6,
                       label=f"{case} ({mot}), mean field (V4)")
    bb = np.linspace(1, 8, 50)
    ax[0].semilogy(bb, np.exp(-bb), "k:", label="$e^{-\\beta}$ (Poisson, V4)")
    ax[0].set_xlabel("sea depth β (pairs per event aperture)")
    ax[0].set_ylabel("share of emissions blocked")
    ax[0].set_title(r"Integer sea at $\nu = 1$: blocking")
    ax[0].legend(fontsize=7)
    ax[1].axhline(0.0, color="k", lw=0.5)
    ax[1].set_xlabel("β")
    ax[1].set_ylabel("ΔT = T(floor) − T(no floor)")
    ax[1].set_title("The transmission's departure")
    ax[1].legend(fontsize=7)
    fig.savefig(output_path(fig_name), dpi=150, bbox_inches="tight")
    dp = docs_path(fig_name)
    if dp:
        fig.savefig(dp, dpi=150, bbox_inches="tight")
    with open(output_path("integer_sea_summary.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["case", "motion", "beta", "worlds", "blocked_share",
                    "aperture_mean", "aperture_p0", "dT", "dT_se",
                    "dT_mean_field"])
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--summarise", nargs="+", metavar="CSV",
                    help="combine finished sweeps and draw the figure;"
                         " no simulation")
    ap.add_argument("--case", choices=("eckart", "well"), default="eckart")
    ap.add_argument("--sea-force", choices=("S", "S'"), default="S'")
    ap.add_argument("--betas", type=float, nargs="+", default=[1, 2, 4, 8])
    ap.add_argument("--seeds", type=int, default=200)
    ap.add_argument("--seed0", type=int, default=1)
    ap.add_argument("--nu", type=int, default=1)
    ap.add_argument("--t-max", type=float, default=None,
                    help="default: 8 (eckart), 12 (well), as in step 24")
    ap.add_argument("--dt", type=float, default=0.02)
    ap.add_argument("--no-gc", action="store_true")
    ap.add_argument("--workers", type=int, default=0, help="0: one per CPU")
    ap.add_argument("--chunk", type=int, default=25)
    ap.add_argument("--collide", nargs="+", default=["none"],
                    metavar="MODE[:SCALE]",
                    help="sea motions run for every beta, from the same"
                         " seeds: none (V6), kernel[:x] (Definition (C),"
                         " r_q B = x |K_q| / 2), uniform[:rB] (neighbouring"
                         " rows)")
    ap.add_argument("--check-poisson", action="store_true",
                    help="only check that collisions keep a Poisson sea"
                         " Poisson (no bodies)")
    ap.add_argument("--tag", default="")
    a = ap.parse_args()
    if a.summarise:
        summarise(a.summarise)
        return
    colls = []
    for c in a.collide:
        mode, _, sc = c.partition(":")
        colls.append(None if mode == "none" else (mode, float(sc or 1.0)))
    if a.check_poisson:
        for coll in colls:
            if coll is None:
                continue
            for b in a.betas:
                t0 = time.time()
                r = poisson_check(a.case, b, a.t_max or 8.0, a.dt, *coll,
                                  range(a.seed0, a.seed0 + a.seeds))
                print(f"  {coll[0]} x{coll[1]:g}  beta {b:g}:  start mean"
                      f" {r['mean0']:.3f} P(0) {r['p00']:.4f};  after"
                      f" mean {r['mean1']:.3f} var {r['var1']:.3f}"
                      f" P(0) {r['p01']:.4f}  (Poisson"
                      f" {np.exp(-r['mean1']):.4f});  collisions per"
                      f" world {r['events']:.0f}  ({time.time() - t0:.0f} s)",
                      flush=True)
        return
    settings = [(b, c) for c in colls for b in a.betas]

    def lab(c):
        return "none" if c is None else f"{c[0]}:{c[1]:g}"
    t_max = a.t_max or (8.0 if a.case == "eckart" else 12.0)
    motion = "S" if a.sea_force == "S" else "Sp"
    seeds = list(range(a.seed0, a.seed0 + a.seeds))
    jobs = [(a.case, a.sea_force.rstrip("'") if a.sea_force == "S" else "Sp",
             settings, seeds[i:i + a.chunk], t_max, a.dt, a.nu, not a.no_gc)
            for i in range(0, len(seeds), a.chunk)]
    workers = a.workers or os.cpu_count() or 1
    print(f"integer sea: case={a.case} motion=({a.sea_force}) nu={a.nu}"
          f" betas={a.betas} collide={[lab(c) for c in colls]}"
          f" seeds={a.seeds} T={t_max} dt={a.dt}"
          f" gc={not a.no_gc} workers={workers}", flush=True)
    t0 = time.time()
    acc = None
    done = 0
    with mp.Pool(workers) as pool:
        for part in pool.imap_unordered(chunk, jobs):
            acc = merge(acc, part)
            done += part["n"]
            n = acc["n"]
            dT, eT = mean_se(acc["dif"][:, 0, 0], acc["dif"][:, 0, 1], n)
            blk = acc["blk"][:, 0] / np.maximum(acc["blk"][:, 1], 1)
            print(f"  {done:6d} seeds  {time.time() - t0:7.0f} s   "
                  + "  ".join(f"{lab(c)} b{b:g}: blk {k:.3f}"
                              f" dT {m:+.4f}±{e:.4f}"
                              for (b, c), k, m, e in zip(settings, blk, dT,
                                                         eT)),
                  flush=True)
    n = acc["n"]
    ref_m, ref_e = mean_se(acc["ref"][:, 0], acc["ref"][:, 1], n)
    rows = []
    print("\n  unfloored: " + "  ".join(
        f"{k} = {m:.4f}±{e:.4f}" for k, m, e in zip(OBS, ref_m, ref_e)))
    c0 = acc["cnt"][0] / n
    print("    per world: " + "  ".join(f"{k} {v:.2f}"
                                       for k, v in zip(CNT, c0)))
    href = acc["href"] / n
    hn = np.linalg.norm(href)
    print(f"\n  {'collide':>10} {'beta':>5} {'blocked':>16} {'e^-b':>7}"
          f" {'<k>':>6} {'P(0)/Pois':>9} {'dT':>16} {'d<p>':>16}"
          f" {'dH':>16} {'eps(E)':>16} {'coll/w':>7}")
    for j, (b, coll) in enumerate(settings):
        dm, de = mean_se(acc["dif"][j, :, 0], acc["dif"][j, :, 1], n)
        vm, ve = mean_se(acc["val"][j, :, 0], acc["val"][j, :, 1], n)
        c = acc["cnt"][j + 1] / n
        sb, sa, sbb, saa, sab = acc["blk"][j] / n
        blk = sb / sa if sa else 0.0
        # delta-method standard error of a ratio of means
        vr = (sbb - 2 * blk * sab + blk * blk * saa
              - (sb - blk * sa) ** 2) / max(sa * sa, 1e-30)
        blk_se = float(np.sqrt(max(vr, 0.0) / max(n - 1, 1)))
        hd = acc["hval"][j] / n - href
        hd_se = np.sqrt(np.maximum(acc["hdif2"][j] / n - hd ** 2, 0.0)
                        / max(n - 1, 1))
        eps = np.linalg.norm(hd) / hn
        eps_se = np.linalg.norm(hd_se) / hn
        oc = acc["occ"][j]
        kk = np.arange(NOCC)
        kmean = (oc * kk).sum() / max(oc.sum(), 1)
        p0 = oc[0] / max(oc.sum(), 1)
        rat = p0 / np.exp(-kmean) if kmean < 50 else float("nan")
        print(f"  {lab(coll):>10} {b:5g} {blk:8.4f}±{blk_se:.4f}"
              f" {np.exp(-b):7.4f} {kmean:6.3f}"
              f" {rat:9.3f}"
              f" {dm[0]:+8.4f}±{de[0]:.4f} {dm[1]:+8.4f}±{de[1]:.4f}"
              f" {dm[2]:+8.4f}±{de[2]:.4f} {eps:8.4f}±{eps_se:.4f}"
              f" {c[CNT.index('coll')]:7.0f}")
        rows.append([a.case, a.sea_force, a.nu, b, n, blk, np.exp(-b),
                     kmean, p0,
                     *dm, *de, *vm, *ve, eps, eps_se,
                     *c, lab(coll), blk_se])
    hdr = (["case", "motion", "nu", "beta", "worlds", "blocked_share",
            "exp_minus_beta", "aperture_mean", "aperture_p0"]
           + [f"d{k}" for k in OBS] + [f"d{k}_se" for k in OBS]
           + list(OBS) + [f"{k}_se" for k in OBS] + ["eps_E", "eps_E_se"]
           + [f"per_world_{k}" for k in CNT] + ["collide", "blocked_se"])
    name = f"integer_sea_{a.case}_{motion}_nu{a.nu}{a.tag}"
    with open(output_path(name + ".csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(hdr)
        w.writerow([a.case, a.sea_force, a.nu, "inf", n, 0.0, 0.0, "", ""]
                   + [0.0] * (2 * len(OBS)) + list(ref_m) + list(ref_e)
                   + [0.0, 0.0] + list(c0) + ["", ""])
        w.writerows(rows)
    np.savez(output_path(name + ".npz"), href=href, occ=acc["occ"],
             hval=acc["hval"] / n, betas=np.array(a.betas), xb=XB, pb=PB,
             worlds=n)
    print(f"\nWrote {name}.csv and .npz in {time.time() - t0:.0f} s.")


if __name__ == "__main__":
    main()
