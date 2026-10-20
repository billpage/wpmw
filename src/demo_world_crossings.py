"""
Companion to docs/supplement/world_crossings.md: crossings of worlds in
configuration space-time as the realisation of an absorption, for N particles.

Problem G of step 25 (three unit-mass particles on a line, harmonic trap,
Gaussian pair interaction g = 2, s = 0.7), worlds as in demo_three_worlds.py.

Parts:
  A   Lemma CX2 (delay is translation): realising an event late by tau equals
      emitting it on time and moving each daughter by -/+ tau K / m in
      configuration space, along the bond, up to O(tau^2).  Checked with the
      flow of problem G integrated finely.  The configuration residual is
      O(tau^3) and the momentum residual O(tau^2).
  B   Lemma CX3 (codimension): for a partner moving at a fixed relative velocity,
      the probability of passing within b of the parent during dt scales as
      b^(n-1) in n configuration dimensions: a true crossing (b = 0) has positive
      probability only for n = 1.
  C   Proposition CX4 (yield): on three-particle ensembles, the number of
      partners per absorption target inside the crossing tube (momentum within
      dp, closest approach within b during dt) against inside the phase-space
      box of PC5 (Chebyshev radius delta), and against dt.
  D   Proposition CX5 (vertex kin): worlds that shared a configuration vertex with
      the parent (its parent and sibling at its birth, its own daughters) as
      partners, and the special case that undoes an earlier event of the same
      parent, bond and channel (re-absorption).  Emission only, no sink.
  F   the figure.

Outputs go through wpmwlib.wpmw_utils.output_path (and docs_path for the figure).
"""
import argparse
import json
import time

import numpy as np
from scipy.spatial import cKDTree

import demo_three_worlds as tw
import demo_population_control as pc
from wpmwlib.wpmw_utils import docs_path, output_path


def flow(x, p, T, n):
    """The classical flow of problem G over time T by n velocity-Verlet steps."""
    h = T / n
    f = tw.force(x, tw.dv_gauss)
    for _ in range(n):
        p = p + 0.5 * h * f
        x = x + h * p
        f = tw.force(x, tw.dv_gauss)
        p = p + 0.5 * h * f
    return x, p


# ----------------------------------------------------------------------------
# Part A -- a delayed event is a translated event
# ----------------------------------------------------------------------------
def part_A(m=4000, seed=1):
    tw.rule("Part A. Lemma CX2: realising an event late by tau = moving its daughter by tau K/m along the bond")
    rng = np.random.default_rng(seed)
    x = rng.normal(0, np.sqrt(0.5), (m, 3))
    p = rng.normal(0, np.sqrt(0.5), (m, 3))
    b = rng.integers(0, 3, m)
    i = np.array([ij[0] for ij in tw.BONDS])[b]
    j = np.array([ij[1] for ij in tw.BONDS])[b]
    xi = 0.25 * rng.integers(1, 32, m)                  # channels xi_q = q dp at reach 2 pi
    K = np.zeros((m, 3))
    K[np.arange(m), i] = xi
    K[np.arange(m), j] = -xi
    rows = []
    print("    tau     max |config residual|/tau^2   max |momentum residual|/tau^2")
    for tau in (0.08, 0.04, 0.02, 0.01):
        xa, pa = flow(x, p, tau, 400)                  # the parent, then the daughter created late at the parent
        xd, pd = flow(x, p + K, tau, 400)              # the daughter created on time, carried by the flow
        rx = xd - (xa + tau * K)                       # Lemma CX2: xd = xa + tau K + O(tau^3)
        rp = pd - (pa + K)                             # and pd = pa + K + O(tau^2)
        ex, ep = np.abs(rx).max() / tau**2, np.abs(rp).max() / tau**2
        cm = np.abs((xd - xa).sum(1)).max()            # the centre of mass is not moved
        rows.append(dict(tau=tau, ex=ex, ep=ep, cm=cm))
        print(f"  {tau:6.3f}       {ex:10.4f}                       {ep:10.4f}              "
              f"(centre-of-mass shift {cm:.1e})")
    return rows


# ----------------------------------------------------------------------------
# Part B -- the codimension of a crossing
# ----------------------------------------------------------------------------
def part_B(m=4_000_000, dt=0.2, speed=1.0, seed=2):
    tw.rule("Part B. Lemma CX3: P(closest approach within b during dt) against b, n = 1, 2, 3")
    rng = np.random.default_rng(seed)
    out = {}
    L = speed * dt                                     # the swept length
    bs = np.array([0.002, 0.004, 0.008, 0.016])
    for n in (1, 2, 3):
        side = 0.5
        d0 = rng.uniform(-side / 2, side / 2, (m, n))  # partner's offset from the parent
        v = np.zeros(n)
        v[0] = -speed                                  # relative velocity, along the bond (axis 0)
        t = np.clip(-(d0 @ v) / (v @ v), 0, dt)       # closest approach of d0 + v t over t in [0, dt]
        miss = np.linalg.norm(d0 + np.outer(t, v), axis=1)
        P = np.array([(miss <= bb).mean() for bb in bs])
        # the tube: a segment of length L thickened by b (cylinder plus two half-balls)
        sec = {1: 2 * bs, 2: np.pi * bs**2, 3: 4 / 3 * np.pi * bs**3}            # volume of an n-ball
        cross = {1: np.ones_like(bs), 2: 2 * bs, 3: np.pi * bs**2}              # its (n-1)-ball cross-section
        pred = (cross[n] * L + sec[n]) / side**n
        P0 = (miss <= 1e-12).mean()
        slope = float(np.polyfit(np.log(bs), np.log(P), 1)[0])
        print(f"  n = {n}: b = {', '.join(f'{bb:.3f}' for bb in bs)}")
        print(f"         measured  P = {', '.join(f'{pp:.3e}' for pp in P)}")
        print(f"         predicted P = {', '.join(f'{pp:.3e}' for pp in pred)}")
        print(f"         P(exact crossing, b = 0) = {P0:.4f} (predicted {L / side if n == 1 else 0:.4f}); "
              f"slope d log P / d log b over these b: {slope:.2f} (n - 1 = {n - 1} for b << L = {L})")
        out[n] = dict(b=bs.tolist(), P=P.tolist(), pred=pred.tolist(), P0=float(P0), slope=slope)
    return out


# ----------------------------------------------------------------------------
# Part C -- partners per target: phase-space box against crossing tube
# ----------------------------------------------------------------------------
def ensemble(kind, nu, T, seed):
    """A three-particle ensemble at time T: 'absorb' = absorption only (PC5,
    paired, delta = 0.25, no sink), 'K1' = K1s10 at cell 0.7."""
    tw.set_particles(3)
    if kind == "absorb":
        w = pc.ControlledWorlds(2, 0.35, nu, seed, sink="none", delta=0.25, paired=True)
    else:
        w = pc.ControlledWorlds(2, 0.7, nu, seed, sink="K1s10")
    f = tw.force(w.x, tw.dv_gauss)
    for _ in range(int(round(T / tw.DT))):
        w.p += 0.5 * tw.DT * f
        w.x += w.p * tw.DT
        f = tw.force(w.x, tw.dv_gauss)
        w.p += 0.5 * tw.DT * f
        w.fire(tw.DT)
        w.sink()
        f = tw.force(w.x, tw.dv_gauss)
    return w


def sample_events(w, m, rng):
    """m events drawn with probability proportional to |K_q(r_b)| over (world, bond, channel):
    parent index, bond, channel, and the + target's sign."""
    rates = np.stack([w.rates(w.x[:, i] - w.x[:, j]) for i, j in tw.BONDS], 1)    # (N, bonds, Q)
    a = np.abs(rates).ravel().astype(float)
    k = rng.choice(a.size, size=m, p=a / a.sum())
    par, b, q = np.unravel_index(k, rates.shape)
    tau = w.eps[par] * np.sign(rates[par, b, q])
    return par, b, q, tau


def partner_counts(w, par, b, q, tau, delta, bs, dts):
    """For each event's + target (sign tau at z_A + K): the number of worlds of sign -tau
    (i) in the box |zeta - zeta_t|_inf <= delta (PC5), and (ii) in the crossing tube, for each
    impact parameter b in bs and step dt in dts: |pi - pi_t|_inf <= delta, and the closest
    approach of rho_B - rho_A + (pi_B - pi_A) t over t in [0, dt] within b (unit masses)."""
    Z = np.concatenate([w.x, w.p], 1)
    rz = pc.rel(Z)
    U = np.array([tw.OJ[:2, i] - tw.OJ[:2, j] for i, j in tw.BONDS])               # bond vectors, |u|^2 = 2
    xi = w.xi[q]
    ra = rz[par]
    rt = ra.copy()
    rt[:, 2:] += xi[:, None] * U[b]                                                # + target in Jacobi momenta
    box = np.zeros(par.size)
    tube = np.zeros((len(bs), len(dts), par.size))
    same = np.zeros(par.size)        # partners at exactly the parent's configuration (born at its vertex)
    for s in (1.0, -1.0):
        sel = np.flatnonzero(tau == s)
        pool = np.flatnonzero(w.eps == -s)
        tree = cKDTree(rz[pool])
        box[sel] = [len(c) for c in tree.query_ball_point(rt[sel], delta, p=np.inf)]
        vmax = np.abs(xi[sel, None] * U[b[sel]]).max(1) + delta
        rmax = np.maximum(delta, max(bs) + vmax * max(dts))
        cands = tree.query_ball_point(rt[sel], rmax, p=np.inf)
        for e, c in zip(sel, cands):
            if not c:
                continue
            zb = rz[pool[c]]
            okp = np.abs(zb[:, 2:] - rt[e, 2:]).max(1) <= delta
            if not okp.any():
                continue
            d0 = zb[okp, :2] - ra[e, :2]                                            # configuration offset
            same[e] = (np.abs(d0).max(1) < 1e-9).sum()
            v = zb[okp, 2:] - ra[e, 2:]                                             # relative velocity
            vv = (v**2).sum(1)
            for k, dt in enumerate(dts):
                t = np.clip(-(d0 * v).sum(1) / vv, 0, dt)
                miss = np.linalg.norm(d0 + t[:, None] * v, axis=1)
                for ib, bb in enumerate(bs):
                    tube[ib, k, e] = (miss <= bb).sum()
    # prediction for a locally uniform density (n = 2 configuration dimensions):
    # tube volume (2 delta)^2 (2 b L + pi b^2), L = |K| dt, over the box volume (2 delta)^4
    L = np.linalg.norm(xi[:, None] * U[b], axis=1)
    pred = np.array([[((2 * delta)**2 * (2 * bb * L * dt + np.pi * bb**2) / (2 * delta)**4).mean() for dt in dts]
                     for bb in bs])
    return box, tube, pred, L, same


def part_C(nu=50000, T=1.0, m=20000, seed=3, delta=0.1, bs=(0.1, 0.02, 0.005), dts=(0.01, 0.02, 0.04, 0.08, 0.16)):
    tw.rule(f"Part C. Proposition CX4: partners per absorption target, box (PC5) against crossing tube; "
            f"nu = {nu}, t = {T}, {m} events, momentum tolerance = box half-side = {delta}")
    rng = np.random.default_rng(seed)
    out = dict(delta=delta, bs=list(bs), dts=list(dts))
    for kind in ("absorb", "K1"):
        t0 = time.perf_counter()
        w = ensemble(kind, nu, T, seed)
        par, b, q, tau = sample_events(w, m, rng)
        print(f"  ensemble '{kind}': {w.eps.size} worlds ({w.eps.size / nu:.2f} per nu), "
              f"{(w.eps < 0).sum()} negatons  ({time.perf_counter() - t0:.0f} s)")
        box, tube, pred, L, same = partner_counts(w, par, b, q, tau, delta, bs, dts)
        print(f"    mean |K| over the events {L.mean():.3f} (|K| dt = {L.mean() * 0.01:.4f} at dt = 0.01)")
        print(f"    box (PC5), half-side {delta}: {box.mean():.4f} partners per target, "
              f"P(at least one) {np.mean(box > 0):.4f}")
        res = []
        for ib, bb in enumerate(bs):
            for k, dt in enumerate(dts):
                r = tube[ib, k].mean() / box.mean()
                print(f"      tube b = {bb:5.3f}, dt = {dt:5.2f}: {tube[ib, k].mean():8.5f} per target, "
                      f"P(>= 1) {np.mean(tube[ib, k] > 0):.5f}; tube/box {r:.4f} (uniform density {pred[ib, k]:.4f})")
                res.append(dict(b=bb, dt=dt, tube=float(tube[ib, k].mean()), p1=float(np.mean(tube[ib, k] > 0)),
                                ratio=float(r), pred=float(pred[ib, k])))
        print(f"    partners at exactly the parent's configuration (born at its vertex): {same.mean():.5f} "
              f"per target, P(>= 1) {np.mean(same > 0):.5f}")
        out[kind] = dict(N=int(w.eps.size), box=float(box.mean()), box_p1=float(np.mean(box > 0)),
                         K=float(L.mean()), same=float(same.mean()), same_p1=float(np.mean(same > 0)), tube=res)
    return out


# ----------------------------------------------------------------------------
# Part D -- vertex kin: partners by descent, and re-absorption
# ----------------------------------------------------------------------------
class FamilyWorlds(pc.ControlledWorlds):
    """ControlledWorlds with world identities and a record, per parent, of the events it
    emitted.  When a parent fires bond b, channel q with sign tau, and an earlier event of
    the same parent, bond and channel had sign -tau and both its daughters are still
    present, the new event can be realised by removing those two daughters (Proposition
    CX5): the + target wants +tau at p + K, the earlier + daughter is -tau there; the -
    target wants -tau at p - K, the earlier - daughter is +tau there.  The only error is
    the delay since the earlier event (Lemma CX2).  Here the opportunities are only
    recorded; every event is realised by emission."""

    def __init__(self, mult, cell, nu, seed, sink="none"):
        assert sink == "none", "world identities are not carried through a sink here"
        super().__init__(mult, cell, nu, seed, sink=sink)
        self.wid = np.arange(nu)
        self.alive = np.ones(nu, bool)
        self.book = {}                 # (parent id, bond, channel, sign) -> [(id +, id -, step), ...]
        self.step = 0
        self.n_avail = 0               # events with an earlier opposite event whose daughters are alive
        self.ages = []                 # age in steps of each available opportunity (the most recent)
        self.kin = {}                  # world id -> ids that shared a configuration vertex with it
        self.kin_dp = 0.05             # momentum tolerance for a vertex-kin partner (force drift since birth)
        self.n_kin_half = 0            # events with a vertex-kin partner for at least one daughter
        self.n_kin_both = 0            # ... for both daughters

    def _new_ids(self, k):
        ids = np.arange(self.alive.size, self.alive.size + k)
        self.alive = np.concatenate([self.alive, np.ones(k, bool)])
        return ids

    def fire(self, dt):
        self.step += 1
        W, Q, B, TAU = [], [], [], []
        for b, (i, j) in enumerate(tw.BONDS):
            k = self.rates(self.x[:, i] - self.x[:, j])
            old = np.floor(self.C[:, b])
            self.C[:, b] += k * dt
            dn = (np.floor(self.C[:, b]) - old).astype(np.int64)
            w, q = np.nonzero(dn)
            if w.size == 0:
                continue
            cnt = np.abs(dn[w, q])
            sg = np.sign(dn[w, q])
            w, q, sg = np.repeat(w, cnt), np.repeat(q, cnt), np.repeat(sg, cnt)
            W.append(w)
            Q.append(q)
            B.append(np.full(w.size, b))
            TAU.append(self.eps[w] * sg)
        if not W:
            return
        w, q, b, tau = (np.concatenate(a) for a in (W, Q, B, TAU))
        m = w.size
        self.fires += m
        self.events_n += m
        bi = np.array([ij[0] for ij in tw.BONDS])[b]
        bj = np.array([ij[1] for ij in tw.BONDS])[b]
        z = np.concatenate([self.x[w], self.p[w]], 1)
        zt = np.concatenate([z, z])
        d = self.xi[q]
        rows = np.arange(m)
        zt[rows, tw.NP + bi] += d
        zt[rows, tw.NP + bj] -= d
        zt[m + rows, tw.NP + bi] -= d
        zt[m + rows, tw.NP + bj] += d
        st = np.concatenate([tau, -tau])
        pid = self.wid[w]
        rzt = pc.rel(zt)
        for e in range(m):
            # vertex kin: worlds that shared a configuration point with the parent (its parent and
            # sibling at its birth, its own daughters); is one of them at a target with the right sign?
            kin = [k_ for k_ in self.kin.get(pid[e], ()) if self.alive[k_]]
            if kin:
                ix = np.searchsorted(self.wid, kin)
                ok = (ix < self.wid.size)
                ix = ix[ok][self.wid[ix[ok]] == np.array(kin)[ok]]
                if ix.size:
                    rk = pc.rel(np.concatenate([self.x[ix], self.p[ix]], 1))
                    hits = []
                    for h in (e, e + m):
                        good = (self.eps[ix] == -st[h]) & (np.abs(rk[:, 2:] - rzt[h, 2:]).max(1) <= self.kin_dp)
                        hits.append(bool(good.any()))
                    self.n_kin_half += hits[0] or hits[1]
                    self.n_kin_both += hits[0] and hits[1]
            stack = self.book.get((pid[e], b[e], q[e], -tau[e]))
            while stack and not (self.alive[stack[-1][0]] and self.alive[stack[-1][1]]):
                stack.pop()                                    # drop records whose daughters are gone
            if stack:
                age = self.step - stack[-1][2]
                self.n_avail += 1
                self.ages.append(age)
        new_ids = self._new_ids(2 * m)
        # record the emitted events: the + daughter (index e) and the - daughter (index e + m)
        for e in range(m):
            ip, im = int(new_ids[e]), int(new_ids[e + m])
            self.book.setdefault((pid[e], b[e], q[e], tau[e]), []).append((ip, im, self.step))
            self.kin[ip] = [int(pid[e]), im]
            self.kin[im] = [int(pid[e]), ip]
            self.kin.setdefault(int(pid[e]), []).extend([ip, im])
        self.x = np.concatenate([self.x, zt[:, :tw.NP]])
        self.p = np.concatenate([self.p, zt[:, tw.NP:]])
        self.eps = np.concatenate([self.eps, st])
        self.wid = np.concatenate([self.wid, new_ids])
        self.C = np.concatenate([self.C, self.rng.random((2 * m, self.nb, self.nq), dtype=np.float32)])


def run_family(nu, seed, T, tag):
    """One run with emission only and no sink, recording vertex-kin opportunities."""
    tw.set_particles(3)
    w = FamilyWorlds(2, 0.7, nu, seed)
    rows, t0 = [], time.perf_counter()
    f = tw.force(w.x, tw.dv_gauss)
    n = int(round(T / tw.DT))
    every = int(round(tw.OBS / tw.DT))
    for i in range(n + 1):
        if i % every == 0:
            t = round(i * tw.DT, 6)
            if abs(t - round(t)) < 1e-9 and t > 0:
                np.savez_compressed(output_path(f"world_crossings_snap_{tag}_s{seed}_t{int(round(t))}.npz"),
                                    r=pc.rel(np.concatenate([w.x, w.p], 1)).astype(np.float32),
                                    e=w.eps.astype(np.int8))
            o = tw.observe(w.x, w.p, w.eps, nu, tw.v_gauss)
            o.update(t=t, events=w.events_n, avail=w.n_avail, kin_half=w.n_kin_half,
                     kin_both=w.n_kin_both)
            rows.append(o)
        if i == n:
            break
        w.p += 0.5 * tw.DT * f
        w.x += w.p * tw.DT
        f = tw.force(w.x, tw.dv_gauss)
        w.p += 0.5 * tw.DT * f
        w.fire(tw.DT)
        w.sink()
        f = tw.force(w.x, tw.dv_gauss)
    ages = np.array(w.ages)
    return dict(nu=nu, seed=seed, T=T, tag=tag, rows=rows,
                age_hist=np.bincount(np.minimum(ages, 50), minlength=51).tolist() if ages.size else [],
                wall=time.perf_counter() - t0)


def part_D(nu=20000, T=0.5):
    tw.rule(f"Part D. Proposition CX5: how often an event can undo an earlier event of its parent "
            f"(emission only, no sink, nu = {nu}, t <= {T})")
    r = run_family(nu, 0, T, "D")
    for row in r["rows"]:
        if row["events"]:
            print(f"  t = {row['t']:4.2f}: events {row['events']:8d}; an undoable earlier event "
                  f"{row['avail'] / row['events']:.4f}; a vertex-kin partner (momentum within 0.05) for one "
                  f"daughter {row['kin_half'] / row['events']:.4f}, for both {row['kin_both'] / row['events']:.4f}; "
                  f"worlds per nu {row['N'] / nu:.2f}")
    h = np.array(r["age_hist"])
    c = np.cumsum(h) / h.sum()
    print(f"  age of the earlier event, in steps of {tw.DT}: median {np.searchsorted(c, 0.5)}, "
          f"P(age <= 1) {c[1]:.3f}, P(age <= 10) {c[10]:.3f}, P(age >= 50) {1 - c[49]:.3f}")
    return r


# ----------------------------------------------------------------------------
# Part F -- figure
# ----------------------------------------------------------------------------
def part_F(C):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ink, muted, grid = "#2b2b2a", "#8a8a85", "#e6e6e3"
    plt.rcParams.update({"axes.edgecolor": muted, "axes.labelcolor": ink, "xtick.color": muted,
                         "ytick.color": muted, "axes.spines.top": False, "axes.spines.right": False,
                         "font.size": 10})
    blue, orange, aqua = "#2a78d6", "#eb6834", "#1baf7a"
    fig, ax = plt.subplots(1, 3, figsize=(16, 5))
    # (a) the vertex in (r_b, t), relative to the parent: the partners converge along the bond
    t = np.linspace(-1, 0, 50)
    ax[0].axvline(0, ymax=0.555, color=ink, lw=2.5, label="parent (momentum p)")
    ax[0].plot(0.8 * t, t, color=blue, lw=2, label="negaton, momentum p + K")
    ax[0].plot(-0.8 * t, t, color=orange, lw=2, label="positon, momentum p − K")
    ax[0].plot([0, 0], [0, 0.8], color=ink, lw=6, alpha=0.25, solid_capstyle="butt",
               label="aligned pair (momentum p)")
    ax[0].plot([0, 0], [0, 0.8], color=ink, lw=2.5)
    ax[0].scatter([0], [0], s=60, color=ink, zorder=5)
    ax[0].annotate("vertex (x*, t*)", (0, 0), xytext=(0.08, -0.12), fontsize=9, color=ink)
    ax[0].annotate("relative velocity +K/m", (-0.48, -0.6), xytext=(-1.05, -0.35), fontsize=8, color=ink,
                   arrowprops=dict(arrowstyle="-", color=muted, lw=0.8))
    ax[0].set(xlabel="bond coordinate relative to the parent", ylabel="t",
              title="(a) an absorption vertex in configuration space-time", xlim=(-1.1, 1.1), ylim=(-1.0, 0.8))
    ax[0].legend(fontsize=8, frameon=False, loc="upper right")
    # (b) the crossing tube in two configuration dimensions
    L, bb = 0.6, 0.12
    ax[1].add_patch(plt.Rectangle((-L, -bb), L, 2 * bb, color=blue, alpha=0.15, lw=0))
    ax[1].add_patch(plt.Circle((0, 0), bb, color=blue, alpha=0.15, lw=0))
    ax[1].add_patch(plt.Circle((-L, 0), bb, color=blue, alpha=0.15, lw=0))
    ax[1].scatter([0], [0], s=60, color=ink, zorder=5)
    ax[1].annotate("parent", (0, 0), xytext=(0.03, -0.06), fontsize=9, color=ink)
    rng = np.random.default_rng(4)
    for y0, x0 in zip(rng.uniform(-0.42, 0.42, 9), rng.uniform(-0.9, -0.05, 9)):
        hit = abs(y0) <= bb and -L <= x0 <= 0
        ax[1].annotate("", xy=(x0 + L, y0), xytext=(x0, y0),
                       arrowprops=dict(arrowstyle="->", color=blue if hit else muted, lw=1.5))
    ax[1].annotate("", xy=(0.0, -0.5), xytext=(-L, -0.5), arrowprops=dict(arrowstyle="<->", color=ink))
    ax[1].text(-L / 2, -0.54, r"$|K|\,\Delta t/m$ along the bond", ha="center", va="top", fontsize=9, color=ink)
    ax[1].annotate("", xy=(0.2, bb), xytext=(0.2, -bb), arrowprops=dict(arrowstyle="<->", color=ink))
    ax[1].text(0.23, 0, "2b", va="center", fontsize=9, color=ink)
    ax[1].set(xlim=(-1.0, 0.7), ylim=(-0.65, 0.5), aspect="equal", xlabel="bond direction",
              ylabel="spectator direction", title="(b) the crossing tube, n = 2 (three particles on a line)")
    # (c) tube / box against dt, measured and predicted, for each impact parameter
    dts = C["dts"]
    cols = {0.1: aqua, 0.02: blue, 0.005: orange}
    for kind, mk, ls in (("absorb", "o", "-"), ("K1", "s", ":")):
        for bb in C["bs"]:
            rows = [r for r in C[kind]["tube"] if r["b"] == bb]
            ax[2].plot(dts, [r["ratio"] for r in rows], mk + ls, color=cols[bb], ms=6, mfc="white" if kind == "K1"
                       else cols[bb], label=f"b = {bb}, {'absorption-only' if kind == 'absorb' else 'K1s10'} ensemble")
    for bb in C["bs"]:
        rows = [r for r in C["absorb"]["tube"] if r["b"] == bb]
        ax[2].plot(dts, [r["pred"] for r in rows], "-", color=muted, lw=1)
    ax[2].plot([], [], "-", color=muted, lw=1, label="uniform-density prediction")
    ax[2].axhline(1, color=ink, lw=0.8)
    ax[2].set(xscale="log", yscale="log", xlabel=r"$\Delta t$",
              ylabel="partners in tube / partners in box",
              title=f"(c) the tube against the box of PC5 (half-side {C['delta']})")
    ax[2].grid(True, color=grid, lw=0.6)
    ax[2].legend(fontsize=7, frameon=False, loc="center", bbox_to_anchor=(0.5, 0.67), ncol=2)
    fig.tight_layout()
    fig.savefig(output_path("world_crossings.png"), dpi=130)
    dp = docs_path("world_crossings.png")
    if dp:
        fig.savefig(dp, dpi=130)
    print("  wrote world_crossings.png")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("parts", nargs="*", default=["A", "B", "C", "D", "F"])
    args = ap.parse_args()
    tw.set_particles(3)
    parts = [p.upper() for p in args.parts]
    t0 = time.perf_counter()
    if "A" in parts:
        tw.dump(part_A(), "world_crossings_A.json")
    if "B" in parts:
        tw.dump(part_B(), "world_crossings_B.json")
    if "C" in parts:
        tw.dump(part_C(), "world_crossings_C.json")
    if "D" in parts:
        tw.dump(part_D(), "world_crossings_D.json")

    if "F" in parts:
        part_F(json.load(open(output_path("world_crossings_C.json"))))
    print(f"\nwall time {time.perf_counter() - t0:.1f} s")


if __name__ == "__main__":
    main()
