"""
Step 26 (docs/analysis/population_control.md): population control that keeps W
below a cell, the open item WC-SP9 of step 25.

The three-particle world ensemble of demo_three_worlds.py (problem G: Gaussian
pair interaction g = 2, s = 0.7, harmonic trap, quench from the product ground
state) loses global purity against the 4D Jacobi mesh because its contact sink
keeps only each cell's net count M0 (Proposition WC10).  This demo replaces the
sink and scores each rule against the mesh at reach 2 pi:

  contact   the contact sink of step 25 (keeps M0)
  K1        moment-preserving resampling (Definition PC3): keep the majority's
            net worlds and add k ± pairs carrying the removed worlds' first moment
            R, so that M0 and M1 of every cell are unchanged; pairs at the cell centre
  K1b       the same with each pair at its dipole centre (Lemma PC4), the grid
            shifted by a random offset at every application
  K1f       dipole centre on a fixed grid, pairs kept inside their cell
  K1sN      dipole centre, grid shifted every N applications (K1s10: Proposition PC6)
  K1sNrX    as K1sN, but in each cell a fraction X of the opposite-sign pairs is
            retained as unpaired worlds instead of removed; the pairs added carry
            only the removed worlds' moment, so M0 and M1 are still kept (PC7)
  K1uN      as K1sN with the dipole centre not clipped into its cell
  none      no sink
  absorption (added to any of these with a tolerance delta): a daughter of sign s
            due at z is realised by removing an existing world of sign -s within
            delta of z (Chebyshev, relative coordinates); mode P (paired: both halves
            or neither, the ledger's absorptive branch), S (each half alone: the
            mixed realisation Theorem S6 excludes, kept as a numerical comparator)
            or D (paired, choosing partners whose displacements agree)
            (Proposition PC5); or, with a K1sN sink, mode C or W (any delta > 0
            switches them on): both daughters absorbed from the nearest
            opposite-sign worlds in their own K1 cells, or both created; in C the
            first-moment error is credited to the K1 step that follows, in W it is
            not (Proposition PC7)

Parts:
  T   tests: K1 and K1b keep M0 and M1 of every cell, with and without retention;
      Lemma PC4 on random pairs; mode C followed by K1 leaves every cell's M0 and
      M1 as emission followed by K1 does (mode W does not); absorption keeps
      sum(eps z) to the matching distance.
  P   one application of each sink to worlds sampled from the mesh at t = 2.5.
  R   runs (--config sink:cell[:delta[:mode]] ...), saving snapshots at t = 1, 2, 3.
  C   cross-seed and self estimators of P_a and P1 from the snapshots (Lemma PC1).
  S   summary of every saved run (self estimators).
  F   the figure.

Outputs go through wpmwlib.wpmw_utils.output_path (and docs_path for the figure).
Parts C, S and F read the mesh reference three_worlds_M2.json
(demo_three_worlds.py M --mesh-reach 2) from the output directory.
"""
import argparse
import copy
import glob
import json
import os
import time

import numpy as np
from scipy.spatial import cKDTree

import demo_three_worlds as tw
from wpmwlib.wpmw_utils import docs_path, output_path


def rel(z6):
    """(n, 6) world points -> (n, 4) Jacobi relative phase-space coordinates."""
    J = tw.OJ[:tw.NP - 1].T
    return np.concatenate([z6[:, :tw.NP] @ J, z6[:, tw.NP:] @ J], 1)


# ----------------------------------------------------------------------------
# the sinks, as functions of (z, eps) so that they can be tested on their own
# ----------------------------------------------------------------------------
def cell_keys(z, cell, off=None):
    """One integer per point: its cell of side `cell` in the relative phase space,
    on the grid shifted by `off`."""
    r = rel(z) if off is None else rel(z) - off
    key = np.zeros(z.shape[0], np.int64)
    for a in r.T:
        key = key * 4096 + np.clip(np.floor(a / cell).astype(np.int64) + 2048, 0, 4095)
    return key


def cells(z, eps, cell, rng, off=None, retain=0.0):
    """Group worlds by cell of the relative phase space; within a cell sort the
    negatons first, at random.  Returns cell ids (original order), the number of
    cells, the opposite-sign pairs per cell, and the mask of worlds the contact
    sink would keep.  retain: in each cell holding m opposite-sign pairs, keep
    floor(retain m + U) of them (U uniform) as unpaired worlds instead of removing
    all m; retain = 0 is the contact sink."""
    key = cell_keys(z, cell, off)
    n = key.size
    o = np.lexsort((rng.random(n), eps, key))
    ks, es = key[o], eps[o]
    start = np.concatenate([[0], np.flatnonzero(np.diff(ks)) + 1])
    cs = np.repeat(np.arange(start.size), np.diff(np.concatenate([start, [n]])))
    npos = np.bincount(cs, weights=(es > 0), minlength=start.size).astype(np.int64)
    nneg = np.bincount(cs, weights=(es < 0), minlength=start.size).astype(np.int64)
    m = np.minimum(npos, nneg)
    if retain > 0:
        m = m - np.floor(retain * m + rng.random(m.size)).astype(np.int64)
    rank = np.arange(n) - np.where(es < 0, start[cs], start[cs] + nneg[cs])
    keep = np.ones(n, bool)
    keep[o[rank < m[cs]]] = False
    cid = np.empty(n, np.int64)
    cid[o] = cs
    return cid, start.size, m, keep


def sink_contact(z, eps, cell, rng):
    _, _, _, keep = cells(z, eps, cell, rng)
    return keep, np.zeros((0, z.shape[1])), np.zeros(0)


def sink_K1(z, eps, cell, rng, return_cells=False):
    """Contact sink plus, in every cell that lost worlds, k ± pairs carrying the
    signed first moment the removed worlds had: R = sum_removed eps z (a balanced
    set, so its M0 is 0).  Pair members sit at c +- R/2k, where c is the cell's
    centre in the relative coordinates (and the removed worlds' mean in the centre
    of mass); k = ceil(|R|_inf / cell), so both members lie inside the cell, and a
    pair symmetric about the centre adds nothing to the cell's central M2."""
    cid, nc, m, keep = cells(z, eps, cell, rng)
    gone = ~keep
    if not gone.any():
        return keep, np.zeros((0, z.shape[1])), np.zeros(0)
    cg = cid[gone]
    D = z.shape[1]
    R = np.stack([np.bincount(cg, weights=eps[gone] * z[gone, d], minlength=nc) for d in range(D)], 1)
    S = np.stack([np.bincount(cg, weights=z[gone, d], minlength=nc) for d in range(D)], 1)
    ng = np.bincount(cg, minlength=nc)
    act = np.flatnonzero(ng > 0)
    R, c = R[act], S[act] / ng[act, None]
    # move the centre's relative part to the cell centre, keeping its centre of mass,
    # so that both members of a pair (|d|_inf < cell) stay inside the cell: a pair
    # straddling a cell edge would look like two unbalanced worlds and never be removed
    first = np.zeros(nc, np.int64)
    first[cg[::-1]] = np.flatnonzero(gone)[::-1]
    rc = (np.floor(rel(z[first[act]]) / cell) + 0.5) * cell
    J = tw.OJ[:tw.NP - 1]
    shift = rc - rel(c)
    c = c + np.concatenate([shift[:, :tw.NP - 1] @ J, shift[:, tw.NP - 1:] @ J], 1)
    k = np.maximum(1, np.ceil(np.abs(rel(R)).max(1) / (0.999 * cell))).astype(np.int64)
    d = np.repeat(R / k[:, None], k, 0)
    c = np.repeat(c, k, 0)
    newz = np.concatenate([c + d / 2, c - d / 2])
    newe = np.concatenate([np.ones(len(d)), -np.ones(len(d))])
    if return_cells:
        nid = np.repeat(act, k)
        return keep, newz, newe, cid, nc, np.concatenate([nid, nid])
    return keep, newz, newe


def sink_K1b(z, eps, cell, rng, return_cells=False, eta=0.0, off=None, inside=False, clip=True,
             retain=0.0, zv=None, ev=None):
    """K1 with the pair placed where its dipole is.  The cell grid is shifted by a
    random offset at every application.  In each cell, R = sum_removed eps z as in
    K1, and the pair centre c' is chosen so that the pair's second moment
    c' R^T + R c'^T matches the removed worlds' A = sum_removed eps z z^T along R:
        c'.R = R^T A R / (2 |R|^2),   c' = (A R - R (c'.R)) / |R|^2.
    For a removed set that is itself one pair this is exact, so a dipole moves with
    the flow instead of being pinned to the cell centre.  The relative part of c'
    is clipped into the cell; the centre of mass is the removed worlds' mean.
    Pairs with |R|_inf < eta * cell are dropped (eta = 0 keeps M1 exactly).
    off: the grid offset (None draws a fresh one); inside=True clips the centre so
    that both members stay in the cell (needed on a fixed grid); clip=False leaves
    the fitted centre where it is, even outside the cell.
    retain: the fraction of each cell's opposite-sign pairs left as unpaired worlds
    (see cells); the pair carries only the removed worlds' moment, so M0 and M1 of
    every cell are still unchanged.
    zv, ev: virtual signed points (balanced within each cell) added to the removed
    set's moments R and A and then discarded: the first-moment credit of a
    within-cell absorption (absorption mode C)."""
    if off is None:
        off = rng.random(4) * cell
    _, _, m, keep = cells(z, eps, cell, rng, off, retain)
    if zv is None:
        zv, ev = np.zeros((0, z.shape[1])), np.zeros(0)
    gone = ~keep
    if not gone.any() and ev.size == 0:
        out = (keep, np.zeros((0, z.shape[1])), np.zeros(0))
        return out + (np.zeros(eps.size, np.int64), 1, np.zeros(0, np.int64)) if return_cells else out
    # one cell id for real and virtual points alike
    uk, inv = np.unique(np.concatenate([cell_keys(z, cell, off), cell_keys(zv, cell, off)]), return_inverse=True)
    nc = uk.size
    cid, cidv = inv[:eps.size], inv[eps.size:]
    cg = np.concatenate([cid[gone], cidv])
    zg = np.concatenate([z[gone], zv])
    eg = np.concatenate([eps[gone], ev])
    D = z.shape[1]
    R = np.stack([np.bincount(cg, weights=eg * zg[:, d], minlength=nc) for d in range(D)], 1)
    S = np.stack([np.bincount(cg, weights=zg[:, d], minlength=nc) for d in range(D)], 1)
    A = np.zeros((nc, D, D))
    for a in range(D):
        for b in range(a, D):
            A[:, a, b] = A[:, b, a] = np.bincount(cg, weights=eg * zg[:, a] * zg[:, b], minlength=nc)
    ng = np.bincount(cg, minlength=nc)
    Rr = np.abs(rel(R)).max(1)
    act = np.flatnonzero((ng > 0) & (Rr > eta * cell) & (Rr > 0))
    R, A, Rr = R[act], A[act], Rr[act]
    R2 = (R**2).sum(1)
    AR = np.einsum("nab,nb->na", A, R)
    cR = (AR * R).sum(1) / (2 * R2)
    c = (AR - R * cR[:, None]) / R2[:, None]
    # centre of mass from the removed worlds' mean; relative part clipped into the cell
    J = tw.OJ[:tw.NP - 1]
    cm = S[act] / ng[act, None]
    first = np.zeros(nc, np.int64)
    first[cg[::-1]] = np.arange(cg.size)[::-1]
    lo = np.floor((rel(zg[first[act]]) - off) / cell) * cell + off
    rc = np.clip(rel(c), lo, lo + cell) if clip else rel(c)
    shift = rc - rel(cm)
    c = cm + np.concatenate([shift[:, :tw.NP - 1] @ J, shift[:, tw.NP - 1:] @ J], 1)
    k = np.maximum(1, np.ceil(Rr / (0.999 * cell))).astype(np.int64)
    if inside:
        h = np.abs(rel(R / k[:, None])) / 2
        eps_ = 1e-9 * cell                      # half-open cells: stay off the upper face
        rc2 = np.clip(rel(c), lo + h + eps_, lo + cell - h - eps_)
        sh2 = rc2 - rel(c)
        c = c + np.concatenate([sh2[:, :tw.NP - 1] @ J, sh2[:, tw.NP - 1:] @ J], 1)
    d = np.repeat(R / k[:, None], k, 0)
    c = np.repeat(c, k, 0)
    newz = np.concatenate([c + d / 2, c - d / 2])
    newe = np.concatenate([np.ones(len(d)), -np.ones(len(d))])
    if return_cells:
        nid = np.repeat(act, k)
        return keep, newz, newe, cid, nc, np.concatenate([nid, nid])
    return keep, newz, newe


def sink_contact_shift(z, eps, cell, rng):
    _, _, _, keep = cells(z, eps, cell, rng, rng.random(4) * cell)
    return keep, np.zeros((0, z.shape[1])), np.zeros(0)


SINKS = {"contact": sink_contact, "contacts": sink_contact_shift, "K1": sink_K1, "K1b": sink_K1b, "none": None}


class ControlledWorlds(tw.Worlds):
    """tw.Worlds with a choice of sink and optional absorptive realisation."""

    def __init__(self, mult, cell, nu, seed, sink="contact", delta=0.0, paired=True):
        super().__init__(mult, cell, nu, seed)
        self.shift_every, self.inside, self.nsink, self.off = 0, False, 0, np.zeros(4)
        self.clip, self.retain, self.off_ready = True, 0.0, False
        self.cz, self.ce = [], []          # first-moment credit of within-cell absorptions
        if sink == "K1f":                       # K1b on a fixed grid, pairs kept inside their cell
            sink, self.inside = "K1b", True
        elif sink.startswith("K1s"):            # K1b with the grid shifted every N applications,
            rest = sink[3:]                     # K1sNrX: and a fraction X of opposite pairs retained
            if "r" in rest:
                rest, rr = rest.split("r")
                self.retain = float(rr)
            sink, self.shift_every = "K1b", int(rest)
        elif sink.startswith("K1u"):            # the same, with the dipole centre not clipped
            sink, self.shift_every, self.clip = "K1b", int(sink[3:]), False
        elif sink == "K1b":
            self.shift_every = 1
        self.kind = sink
        self.sinkf = SINKS[sink]
        self.delta = delta
        self.mode = {True: "P", False: "S"}.get(paired, paired)
        self.paired = self.mode in ("P", "D", "C", "W")
        self.credit = 0.0          # sum of |credited first moment|_inf (mode C)
        self.events_n = 0          # events fired
        self.halves_abs = 0        # daughters realised by absorption
        self.events_abs = 0        # events realised with both halves absorbed
        self.dipole = 0.0          # sum of |z_partner - z_target|_inf over absorptions

    # -- events ---------------------------------------------------------------
    def fire(self, dt):
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
        # daughters: +tau at p + xi u_b, -tau at p - xi u_b (both at the parent's x)
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
        create = np.ones(2 * m, bool)
        remove = np.zeros(self.eps.size, bool)
        if self.delta > 0 and self.mode in ("C", "W"):
            create, remove = self.absorb_cell(zt, st, m)
        elif self.delta > 0:
            create, remove = self.absorb(zt, st, m)
        self.x = np.concatenate([self.x[~remove], zt[create, :tw.NP]])
        self.p = np.concatenate([self.p[~remove], zt[create, tw.NP:]])
        self.eps = np.concatenate([self.eps[~remove], st[create]])
        self.C = np.concatenate([self.C[~remove],
                                 self.rng.random((int(create.sum()), self.nb, self.nq), dtype=np.float32)])

    def absorb(self, zt, st, m):
        """For each daughter (target zt, sign st) look for an existing world of sign
        -st within delta (Chebyshev, relative coordinates); greedy, each partner
        used once.  Returns (create mask over daughters, remove mask over worlds)."""
        rz = rel(np.concatenate([self.x, self.p], 1))
        rt = rel(zt)
        n = self.eps.size
        remove = np.zeros(n, bool)
        create = np.ones(2 * m, bool)
        cand = np.full((2 * m, 8), -1, np.int64)
        dist = np.full((2 * m, 8), np.inf)
        for s in (1.0, -1.0):
            pool = np.flatnonzero(self.eps == -s)
            want = np.flatnonzero(st == s)
            if pool.size == 0 or want.size == 0:
                continue
            tree = cKDTree(rz[pool])
            dd, ii = tree.query(rt[want], k=min(8, pool.size), p=np.inf, distance_upper_bound=self.delta)
            dd, ii = dd.reshape(want.size, -1), ii.reshape(want.size, -1)
            ok = np.isfinite(dd)
            cand[want, :ii.shape[1]] = np.where(ok, pool[np.minimum(ii, pool.size - 1)], -1)
            dist[want, :dd.shape[1]] = dd
        order = self.rng.permutation(m)
        for e in order:
            pick = []
            for h in (e, e + m):
                got = -1
                for c, dc in zip(cand[h], dist[h]):
                    if c < 0:
                        break
                    if not remove[c] and c not in pick:
                        got = c
                        break
                pick.append(got)
            if self.mode == "D":
                # choose the partner pair whose displacements a, b from the two targets
                # agree best: the event's change is then reproduced with first-moment error
                # tau (a - b), i.e. the absorbed pair is the event's own dipole, translated
                best, bc = np.inf, None
                ca = [c for c in cand[e] if c >= 0 and not remove[c]]
                cb = [c for c in cand[e + m] if c >= 0 and not remove[c]]
                for c1 in ca:
                    a = rz[c1] - rt[e]
                    for c2 in cb:
                        if c2 == c1:
                            continue
                        v = np.abs(a - (rz[c2] - rt[e + m])).max()
                        if v < best:
                            best, bc = v, (c1, c2)
                pick = list(bc) if bc else [-1, -1]
            if self.paired:
                if pick[0] >= 0 and pick[1] >= 0:
                    for h, c in zip((e, e + m), pick):
                        remove[c] = True
                        create[h] = False
                        self.dipole += float(np.abs(rz[c] - rt[h]).max())
                    self.halves_abs += 2
                    self.events_abs += 1
            else:
                both = True
                for h, c in zip((e, e + m), pick):
                    if c >= 0:
                        remove[c] = True
                        create[h] = False
                        self.halves_abs += 1
                        self.dipole += float(np.abs(rz[c] - rt[h]).max())
                    else:
                        both = False
                self.events_abs += both
        return create, remove

    def grid_offset(self):
        """The K1 grid offset for the current step, drawn once and shared by the
        within-cell absorption (in fire) and the sink that follows it."""
        if not self.off_ready:
            if self.shift_every and self.nsink % self.shift_every == 0:
                self.off = self.rng.random(4) * self.cell
            self.off_ready = True
        return self.off

    def absorb_cell(self, zt, st, m):
        """Modes C and W: an event's two daughters are realised by absorbing, for
        each, the nearest unused world of the opposite sign in the daughter's own
        K1 cell, or both are created.  In mode C the first-moment error is credited
        to the cell: a daughter of sign s due at z_t, realised by removing the
        partner at z_p, leaves +s at z_p where +s at z_t was intended, so the virtual
        points (+s at z_t, -s at z_p) are handed to the K1 step, which restores the
        cell's moments to those emission would have left."""
        off = self.grid_offset()
        Z = np.concatenate([self.x, self.p], 1)
        kz, kt = cell_keys(Z, self.cell, off), cell_keys(zt, self.cell, off)
        rz, rt = rel(Z), rel(zt)
        remove = np.zeros(self.eps.size, bool)
        create = np.ones(2 * m, bool)
        pools = {}
        for sgn in (1.0, -1.0):
            idx = np.flatnonzero(self.eps == sgn)
            o = np.argsort(kz[idx], kind="stable")
            pools[sgn] = (kz[idx][o], idx[o])
        for e in self.rng.permutation(m):
            pick = []
            for h in (e, e + m):
                keys, idx = pools[-st[h]]
                lo, hi = np.searchsorted(keys, kt[h], "left"), np.searchsorted(keys, kt[h], "right")
                c = idx[lo:hi]
                c = c[~remove[c]]
                if pick:
                    c = c[c != pick[0]]
                pick.append(int(c[np.abs(rz[c] - rt[h]).max(1).argmin()]) if c.size else -1)
            if pick[0] < 0 or pick[1] < 0:
                continue
            for h, c in zip((e, e + m), pick):
                remove[c] = True
                create[h] = False
                self.dipole += float(np.abs(rz[c] - rt[h]).max())
                if self.mode == "C":
                    self.cz += [zt[h], Z[c]]
                    self.ce += [st[h], -st[h]]
                    self.credit += float(np.abs(rz[c] - rt[h]).max())
            self.halves_abs += 2
            self.events_abs += 1
        return create, remove

    # -- population control -----------------------------------------------------
    def sink(self):
        if self.sinkf is None:
            return
        z = np.concatenate([self.x, self.p], 1)
        if self.kind == "K1b":
            off = self.grid_offset()
            self.nsink += 1
            self.off_ready = False
            zv = np.array(self.cz).reshape(-1, z.shape[1]) if self.cz else None
            ev = np.array(self.ce, float) if self.ce else None
            self.cz, self.ce = [], []
            keep, nz, ne = sink_K1b(z, self.eps, self.cell, self.rng, off=off, inside=self.inside,
                                    clip=self.clip, retain=self.retain, zv=zv, ev=ev)
        else:
            keep, nz, ne = self.sinkf(z, self.eps, self.cell, self.rng)
        self.x = np.concatenate([self.x[keep], nz[:, :tw.NP]])
        self.p = np.concatenate([self.p[keep], nz[:, tw.NP:]])
        self.eps = np.concatenate([self.eps[keep], ne])
        self.C = np.concatenate([self.C[keep], self.rng.random((ne.size, self.nb, self.nq), dtype=np.float32)])


def run_job(job):
    mult, cell, nu, seed, T, sink, delta, paired, cap = job
    snap = True
    tw.set_particles(3)
    w = ControlledWorlds(mult, cell, nu, seed, sink, delta, paired)
    w.verbose = seed == 0
    rows, t0 = [], time.perf_counter()
    every = int(round(tw.OBS / tw.DT))
    f = tw.force(w.x, tw.dv_gauss)
    nsteps = int(round(T / tw.DT))
    for i in range(nsteps + 1):
        if i % every == 0:
            t_now = round(i * tw.DT, 6)
            if snap and abs(t_now - round(t_now)) < 1e-9 and t_now > 0:
                np.savez_compressed(output_path(f"popctl_snap_{label(sink, cell, delta, paired)}_nu{nu}"
                                                f"_s{seed}_t{int(round(t_now))}.npz"),
                                    r=rel(np.concatenate([w.x, w.p], 1)).astype(np.float32),
                                    e=w.eps.astype(np.int8))
            o = tw.observe(w.x, w.p, w.eps, nu, tw.v_gauss)
            o.update(t=round(i * tw.DT, 6), events=w.events_n, halves_abs=w.halves_abs,
                     events_abs=w.events_abs, dipole=w.dipole, credit=w.credit)
            rows.append(o)
            if w.verbose:
                fa = w.halves_abs / max(1, 2 * w.events_n)
                print(f"    [{sink}, c {cell}, delta {delta}] t = {i * tw.DT:4.2f}: P_a {o['gpurity']:.4f}, "
                      f"P1 {o['purity']:.4f}, N/nu {w.eps.size / nu:.2f}, absorbed halves {fa:.3f}"
                      f"  ({time.perf_counter() - t0:.0f} s)", flush=True)
        if i == nsteps or w.eps.size > cap * nu:
            break
        w.p += 0.5 * tw.DT * f
        w.x += w.p * tw.DT
        f = tw.force(w.x, tw.dv_gauss)
        w.p += 0.5 * tw.DT * f
        w.fire(tw.DT)
        w.sink()
        f = tw.force(w.x, tw.dv_gauss)
        w.peak = max(w.peak, w.eps.size)
    return dict(mult=mult, cell=cell, nu=nu, seed=seed, T=T, sink=sink, delta=delta, paired=paired,
                rows=rows, peak=w.peak, wall=time.perf_counter() - t0)


def label(sink, cell, delta, paired):
    s = f"{sink}_c{cell}"
    if delta > 0:
        mode = {True: "P", False: "S"}.get(paired, paired)
        s += f"_abs{mode}{delta}"
    return s


# ----------------------------------------------------------------------------
# Part T -- unit tests
# ----------------------------------------------------------------------------
def part_T():
    tw.rule("Part T. Tests: K1 and K1b keep M0 and M1 of every cell; Lemma PC4; absorption keeps sum(eps z) "
            "to the matching distance")
    tw.set_particles(3)
    rng = np.random.default_rng(1)
    n = 200000
    z = rng.normal(0, 0.8, (n, 6))
    eps = np.where(rng.random(n) < 0.5 + 0.3 * np.tanh(3 * np.sin(4 * z[:, 3])), 1.0, -1.0)
    out = {}
    for cell, name in ((0.35, "K1"), (0.7, "K1"), (0.35, "K1b"), (0.7, "K1b")):
        keep, nz, ne, cid, nc, nid = SINKS[name](z, eps, cell, np.random.default_rng(2), return_cells=True)
        cid_new = np.concatenate([cid[keep], nid])
        zz, ee = np.concatenate([z[keep], nz]), np.concatenate([eps[keep], ne])
        M0o = np.bincount(cid, weights=eps, minlength=nc)
        M0n = np.bincount(cid_new, weights=ee, minlength=nc)
        M1o = np.stack([np.bincount(cid, weights=eps * z[:, d], minlength=nc) for d in range(6)], 1)
        M1n = np.stack([np.bincount(cid_new, weights=ee * zz[:, d], minlength=nc) for d in range(6)], 1)
        M1c = np.stack([np.bincount(cid[keep], weights=eps[keep] * z[keep, d], minlength=nc) for d in range(6)], 1)
        kc = sink_contact(z, eps, cell, np.random.default_rng(2))[0]
        dmax = np.abs(rel(nz[:len(nz) // 2] - nz[len(nz) // 2:])).max()
        print(f"  {name:3s} cell {cell}: worlds {n} -> contact {kc.sum()} ({kc.mean():.3f}), {name} {ee.size} "
              f"({ee.size / n:.3f}); "
              f"per cell max |dM0| {np.abs(M0n - M0o).max():.1e}, max |dM1| {np.abs(M1n - M1o).max():.1e} "
              f"(contact: {np.abs(M1c - M1o).max():.2f}); pair separation max {dmax:.3f} <= cell")
        out[f"{name}_{cell}"] = dict(contact=float(kc.mean()), kept=ee.size / n,
                                     dM0=float(np.abs(M0n - M0o).max()), dM1=float(np.abs(M1n - M1o).max()),
                                     dM1_contact=float(np.abs(M1c - M1o).max()))
    # Lemma PC4: a cell holding one pair (any centre, any moment) gives that pair back
    worst = 0.0
    for trial in range(200):
        c0 = rng.uniform(-2, 2, 6)
        d0 = rng.normal(0, 0.1, 6)
        z2 = np.stack([c0 + d0 / 2, c0 - d0 / 2])
        R = d0
        A = np.outer(z2[0], z2[0]) - np.outer(z2[1], z2[1])
        cR = R @ A @ R / (2 * R @ R)
        cp = (A @ R - R * cR) / (R @ R)
        worst = max(worst, float(np.abs(cp - c0).max()))
    print(f"  Lemma PC4: the dipole centre of a lone pair returns its centre to {worst:.1e} (200 random pairs)")
    out["PC4"] = worst
    # K1 with retention keeps M0 and M1 of every cell (inside=True keeps pairs in their cell)
    for retain in (0.25, 0.5):
        cid, nc = None, None
        keep, nz, ne, cid, nc, nid = sink_K1b(z, eps, 0.7, np.random.default_rng(2), return_cells=True,
                                              off=np.zeros(4), inside=True, retain=retain)
        zz, ee = np.concatenate([z[keep], nz]), np.concatenate([eps[keep], ne])
        cn = np.concatenate([cid[keep], nid])
        dM1 = max(np.abs(np.bincount(cn, weights=ee * zz[:, d], minlength=nc)
                         - np.bincount(cid, weights=eps * z[:, d], minlength=nc)).max() for d in range(6))
        dM0 = np.abs(np.bincount(cn, weights=ee, minlength=nc) - np.bincount(cid, weights=eps, minlength=nc)).max()
        both = np.bincount(cid[keep], weights=(eps[keep] > 0), minlength=nc) * \
            np.bincount(cid[keep], weights=(eps[keep] < 0), minlength=nc)
        print(f"  K1 retain {retain}: worlds {n} -> {ee.size} ({ee.size / n:.3f}); per cell max |dM0| {dM0:.1e}, "
              f"max |dM1| {dM1:.1e}; cells left holding both signs (besides K1 pairs) {int((both > 0).sum())}")
        out[f"retain_{retain}"] = dict(kept=ee.size / n, dM0=float(dM0), dM1=float(dM1))
    # mode C: absorption within the cell plus the K1 credit leaves every cell's M0 and M1
    # exactly as emission followed by K1 does (same events, fixed grid, pairs kept inside)
    w0 = ControlledWorlds(2, 0.7, 20000, 5, sink="K1f", delta=0.0)
    w0.retain = 0.5
    for _ in range(40):
        w0.fire(0.01)
        w0.sink()
    for mode in ("C", "W"):
        wc, we = copy.deepcopy(w0), copy.deepcopy(w0)
        wc.delta, wc.mode, wc.paired = 0.7, mode, True
        h0 = wc.halves_abs
        for ww in (wc, we):
            ww.fire(0.01)
            ww.sink()
        res = []
        for ww in (wc, we):
            Z = np.concatenate([ww.x, ww.p], 1)
            kz = cell_keys(Z, 0.7, np.zeros(4))
            res.append((kz, ww.eps, Z))
        allk = np.unique(np.concatenate([res[0][0], res[1][0]]))
        mom = []
        for kz, ee, Z in res:
            ix = np.searchsorted(allk, kz)
            mom.append(np.stack([np.bincount(ix, weights=ee, minlength=allk.size)]
                                + [np.bincount(ix, weights=ee * Z[:, d], minlength=allk.size) for d in range(6)], 1))
        err = np.abs(mom[0] - mom[1]).max(0)
        print(f"  mode {mode}: {(wc.halves_abs - h0) // 2} of {wc.events_n - w0.events_n} events absorbed; "
              f"per cell, absorption + K1 against emission + K1: max |dM0| {err[0]:.1e}, "
              f"max |dM1| {err[1:].max():.1e}; "
              f"worlds {we.eps.size} (emission) -> {wc.eps.size}")
        out[f"mode_{mode}"] = dict(dM0=float(err[0]), dM1=float(err[1:].max()), N=int(wc.eps.size),
                                   N_emit=int(we.eps.size))
    # absorption: the change of sum(eps z_rel) equals emission's to within (absorbed halves) x delta per component
    for paired in (True, False, "D"):
        w = ControlledWorlds(2, 0.35, 20000, 4, sink="none", delta=0.15, paired=paired)
        for _ in range(30):                                   # grow a signed population first
            w.fire(0.01)
        w.delta = 0.15
        z0 = rel(np.concatenate([w.x, w.p], 1))
        S0 = (w.eps[:, None] * z0).sum(0)
        h0, e0 = w.halves_abs, w.events_n
        # what pure emission would add: rerun the same step on a copy with delta = 0
        w2 = copy.deepcopy(w)
        w2.delta = 0.0
        w2.rng = copy.deepcopy(w.rng)
        w.fire(0.01)
        w2.fire(0.01)
        S1 = (w.eps[:, None] * rel(np.concatenate([w.x, w.p], 1))).sum(0)
        S2 = (w2.eps[:, None] * rel(np.concatenate([w2.x, w2.p], 1))).sum(0)
        dh = w.halves_abs - h0
        err = np.abs((S1 - S0) - (S2 - S0)).max()
        kind = {True: "paired", False: "single", "D": "dipole"}[paired]
        print(f"  absorb ({kind}): {w.events_n - e0} events, {dh} halves absorbed; "
              f"|d sum(eps z) - emission's| = {err:.3f} <= {dh} x 0.15 = {dh * 0.15:.2f}; "
              f"N {w2.eps.size} (emission) -> {w.eps.size}; net weight {w.eps.sum():.0f} = {w2.eps.sum():.0f}")
        out[f"absorb_{paired}"] = dict(err=float(err), bound=dh * 0.15, N=int(w.eps.size), N_emit=int(w2.eps.size))
    return out


# ----------------------------------------------------------------------------
# Part P -- one sink applied to an ensemble sampled from the converged mesh
# ----------------------------------------------------------------------------
def part_P(T=2.5, nu=1000000):
    tw.rule(f"Part P. One application of each rule to {nu} worlds sampled from the mesh W at t = {T} (reach 2 pi)")
    tw.set_particles(3)
    rng = np.random.default_rng(3)
    mult, nr, Lr, Lp, dt = 2, 40, 5.0, 5.0, tw.DT
    W, (rho, pi_, drho, dpi) = mesh_state(mult, T, nr, Lr, Lp, dt)
    wgt = W.ravel() * drho**2 * dpi**2
    Z = np.abs(wgt).sum()
    idx = rng.choice(wgt.size, size=nu, p=np.abs(wgt) / Z)
    eps = np.sign(wgt[idx])
    i1, i2, j1, j2 = np.unravel_index(idx, W.shape)
    r4 = np.stack([rho[i1] + drho * (rng.random(nu) - 0.5), rho[i2] + drho * (rng.random(nu) - 0.5),
                   pi_[j1] + dpi * (rng.random(nu) - 0.5), pi_[j2] + dpi * (rng.random(nu) - 0.5)], 1)
    # embed in 6D with a centre of mass from its ground state (the observables read relative coordinates)
    Jr = tw.OJ[:2]
    cm = tw.OJ[2]
    X = r4[:, :2] @ Jr + np.outer(rng.normal(0, np.sqrt(0.5), nu), cm)
    P = r4[:, 2:] @ Jr + np.outer(rng.normal(0, np.sqrt(0.5), nu), cm)
    z = np.concatenate([X, P], 1)
    scale = Z

    def score(z, e):
        o = tw.observe(z[:, :3], z[:, 3:], e * scale, nu, tw.v_gauss)
        return o["gpurity"], o["purity"]

    g0, p0 = score(z, eps)
    gm = tw.global_purity_mesh(W, (drho, drho, dpi, dpi))
    print(f"  mesh P_a {gm:.4f}; sampled ensemble P_a {g0:.4f}, P1 {p0:.4f}, |W| mass {Z:.3f}")
    out = dict(mesh=gm, P_a0=g0, P10=p0, rows=[])
    for cell in (0.175, 0.25, 0.35, 0.5, 0.7):
        for name in ("contact", "K1", "K1b"):
            keep, nz, ne = SINKS[name](z, eps, cell, np.random.default_rng(7))
            zz, ee = np.concatenate([z[keep], nz]), np.concatenate([eps[keep], ne])
            g, p = score(zz, ee)
            print(f"  {name:7s} cell {cell:5.3f}: kept {ee.size / nu:.3f} of the worlds; P_a {g:.4f} ({g - g0:+.4f}), "
                  f"P1 {p:.4f} ({p - p0:+.4f})", flush=True)
            out["rows"].append(dict(sink=name, cell=cell, kept=ee.size / nu, P_a=g, P1=p))
    return out


def mesh_state(mult, T, nr, Lr, Lp, dt):
    """The 4D Jacobi mesh of demo_three_worlds.mesh_run, returning W at time T."""
    ymax = mult * np.pi
    dpi = np.pi / (np.sqrt(2) * ymax) * 0.98
    npi = int(np.ceil(Lp / dpi)) * 2
    dpi = 2 * Lp / npi
    rho = -Lr + 2 * Lr * np.arange(nr) / nr
    drho = rho[1] - rho[0]
    pi_ = -Lp + dpi * np.arange(npi)
    R1, R2 = np.meshgrid(rho, rho, indexing="ij")
    s = 2 * np.pi * np.fft.fftfreq(npi, d=dpi)
    S1, S2 = np.meshgrid(s, s, indexing="ij")
    kr = 2 * np.pi * np.fft.fftfreq(nr, d=drho)
    dV1, dV2 = tw.OMEGA**2 * R1, tw.OMEGA**2 * R2
    expo = np.zeros((nr, nr, npi, npi), complex)
    for u in tw.UB:
        rb = u[0] * R1 + u[1] * R2
        dd = tw.dv_gauss(rb)
        dV1 = dV1 + u[0] * dd
        dV2 = dV2 + u[1] * dd
        expo += tw.bond_symbol(rb[:, :, None, None], (u[0] * S1 + u[1] * S2)[None, None], tw.v_gauss,
                               tw.dv_gauss, ymax)
    expo += 1j * (dV1[:, :, None, None] * S1[None, None] + dV2[:, :, None, None] * S2[None, None])
    kick = np.exp(dt * expo)
    del expo
    K1, K2 = np.meshgrid(kr, kr, indexing="ij")
    P1, P2 = np.meshgrid(pi_, pi_, indexing="ij")
    drift = np.exp(-1j * (K1[:, :, None, None] * P1[None, None] + K2[:, :, None, None] * P2[None, None]) * dt / 2)
    W = np.exp(-R1**2 - R2**2)[:, :, None, None] * np.exp(-P1**2 - P2**2)[None, None]
    W /= W.sum() * drho**2 * dpi**2
    for _ in range(int(round(T / dt))):
        W = np.real(np.fft.ifft2(np.fft.fft2(W, axes=(0, 1)) * drift, axes=(0, 1)))
        W = np.real(np.fft.ifft2(np.fft.fft2(W, axes=(2, 3)) * kick, axes=(2, 3)))
        W = np.real(np.fft.ifft2(np.fft.fft2(W, axes=(0, 1)) * drift, axes=(0, 1)))
    return W, (rho, pi_, drho, dpi)


# ----------------------------------------------------------------------------
# Part R -- runs against the mesh
# ----------------------------------------------------------------------------
def part_R(configs, nu, seeds, T, workers, cap):
    import multiprocessing as mp
    tw.rule(f"Part R. Runs, reach 2 pi, nu = {nu}, {seeds} seed(s), t = 0 .. {T}")
    jobs = [(2, c, nu, s, T, sk, dl, pr, cap) for (sk, c, dl, pr) in configs for s in range(seeds)]
    res = {}
    with mp.Pool(min(workers, len(jobs))) as pool:
        for r in pool.imap_unordered(run_job, jobs):
            lab = label(r["sink"], r["cell"], r["delta"], r["paired"])
            f = r["rows"][-1]
            print(f"  {lab} seed {r['seed']}: t = {f['t']}: P_a {f['gpurity']:.4f}, P1 {f['purity']:.4f}, "
                  f"N/nu {f['N'] / nu:.2f} (peak {r['peak'] / nu:.2f})  ({r['wall']:.0f} s)", flush=True)
            res.setdefault(lab, []).append(r)
    for lab, rr in res.items():
        tw.dump(rr, f"popctl_R_{lab}_nu{nu}.json")
    return res


def load_runs():
    runs = {}
    for f in sorted(glob.glob(output_path("popctl_R_*.json"))):
        rr = json.load(open(f))
        lab = label(rr[0]["sink"], rr[0]["cell"], rr[0]["delta"], rr[0]["paired"])
        runs[f"{lab}_nu{rr[0]['nu']}"] = rr
    return runs


def summary(mesh):
    tw.rule("Summary: P_a and P1 (self estimators, mean ± se over seeds) against the mesh at reach 2 pi")
    mr = {r["t"]: r for r in mesh["rows"]}
    out = {}
    print(f"  {'rule':34s}   t   P_a - mesh        P1 - mesh        N/nu   absorbed halves  mean |d|")
    for lab, rr in load_runs().items():
        out[lab] = []
        for k in range(len(rr[0]["rows"])):
            rows = [r["rows"][k] for r in rr if len(r["rows"]) > k]
            t = rows[0]["t"]
            g = np.array([x["gpurity"] for x in rows]) - mr[t]["gpurity"]
            p = np.array([x["purity"] for x in rows]) - mr[t]["purity"]
            N = np.mean([x["N"] for x in rows]) / rr[0]["nu"]
            fa = np.mean([x["halves_abs"] / max(1, 2 * x["events"]) for x in rows])
            dip = np.mean([x["dipole"] / max(1, x["halves_abs"]) for x in rows])
            se = (lambda a: a.std(ddof=1) / np.sqrt(a.size) if a.size > 1 else 0.0)
            out[lab].append(dict(t=t, dPa=g.mean(), dPa_se=se(g), dP1=p.mean(), dP1_se=se(p), N=N, fabs=fa, dip=dip))
            if abs(t - round(t)) < 1e-9 and t > 0:
                print(f"  {lab:34s} {t:4.1f}  {g.mean():+.4f} ± {se(g):.4f}  {p.mean():+.4f} ± {se(p):.4f}  "
                      f"{N:5.2f}   {fa:.3f}           {dip:.3f}")
    return out


# ----------------------------------------------------------------------------
# Part C -- cross-seed estimators: no variance term
# ----------------------------------------------------------------------------
def _cic4(r, e):
    d = r.shape[1]
    nb = int(round(2 * tw.GEXT / tw.GH))
    u = (r + tw.GEXT) / tw.GH - 0.5
    i0 = np.floor(u).astype(np.int64)
    f = u - i0
    H = np.zeros(nb**d)
    for corner in range(2**d):
        idx = np.zeros(e.size, np.int64)
        wv = e.copy()
        ok = np.ones(e.size, bool)
        for ax in range(d):
            bit = (corner >> ax) & 1
            ii = i0[:, ax] + bit
            ok &= (ii >= 0) & (ii < nb)
            idx = idx * nb + np.clip(ii, 0, nb - 1)
            wv = wv * (f[:, ax] if bit else 1 - f[:, ax])
        H += np.bincount(idx[ok], weights=wv[ok], minlength=nb**d)
    return H.reshape((nb,) * d) / tw.GH**d


def cross_Pa(r1, e1, r2, e2):
    """(2 pi hbar)^2 * integral (W1*G_a)(W2*G_a) for two independent ensembles: the
    expectation is P_a of the mean W, with no variance term (unlike the self estimator)."""
    F1, F2 = np.fft.fftn(_cic4(r1, e1)), np.fft.fftn(_cic4(r2, e2))
    g, ks = tw._gauss_filter(F1.shape, (tw.GH,) * 4)
    T = np.ones(F1.shape)
    nb = F1.shape[0]
    for ax, k in enumerate(ks):
        sh = [1] * 4
        sh[ax] = nb
        T = T * (np.sinc(k * tw.GH / (2 * np.pi))**2).reshape(sh)
    total = tw.GH**4 / F1.size * np.sum(np.real(F1 * np.conj(F2)) / T**2 * g)
    return float((2 * np.pi * tw.HB)**2 * total)


def cross_P1(r1, e1, r2, e2):
    """The same for the one-body purity: projection on (a.rho, a.pi), exact CM Gaussian."""
    nb = int(round(2 * tw.PEXT / tw.PGRID))
    hs = []
    for r, e in ((r1, e1), (r2, e2)):
        z = np.stack([r[:, :2] @ tw.A1, r[:, 2:] @ tw.A1], 1)
        hs.append(np.fft.fft2(tw.deposit_cic(z, e, nb)))
    k = 2 * np.pi * np.fft.fftfreq(nb, d=tw.PGRID)
    g = np.exp(-0.5 * tw.VAR_CM * (k[:, None]**2 + k[None, :]**2))
    sc = np.sinc(k * tw.PGRID / (2 * np.pi))**2
    sc2 = (sc[:, None] * sc[None, :])**2
    total = float(np.sum(np.real(hs[0] * np.conj(hs[1])) / sc2 * g**2)) / nb**2 / tw.PGRID**2
    return 2 * np.pi * tw.HB * total


def part_C(mesh, nu=100000):
    tw.rule("Part C. Cross-seed P_a and P1 (seed 0 x seed 1) against the self estimator and the mesh")
    mr = {r["t"]: r for r in mesh["rows"]}
    out = {}
    labs = sorted({os.path.basename(f).split("_nu")[0][len("popctl_snap_"):]
                   for f in glob.glob(output_path(f"popctl_snap_*_nu{nu}_s0_t*.npz"))})
    print(f"  {'rule':26s} t   P_a cross - mesh   self - mesh   |  P1 cross - mesh   self - mesh  |  N/nu")
    for lab in labs:
        out[lab] = []
        for t in (1, 2, 3):
            f0 = output_path(f"popctl_snap_{lab}_nu{nu}_s0_t{t}.npz")
            f1 = output_path(f"popctl_snap_{lab}_nu{nu}_s1_t{t}.npz")
            if not (os.path.exists(f0) and os.path.exists(f1)):
                continue
            a, b = np.load(f0), np.load(f1)
            r0, e0 = a["r"].astype(float), a["e"].astype(float) / nu
            r1, e1 = b["r"].astype(float), b["e"].astype(float) / nu
            gx, px = cross_Pa(r0, e0, r1, e1), cross_P1(r0, e0, r1, e1)
            gs = 0.5 * (tw.global_purity_points(r0, e0) + tw.global_purity_points(r1, e1))
            ps = []
            for r, e in ((r0, e0), (r1, e1)):
                z = np.stack([r[:, :2] @ tw.A1, r[:, 2:] @ tw.A1], 1)
                st = np.sum(e**2) / (4 * np.pi * tw.VAR_CM) * 2 * np.pi * tw.HB
                ps.append(tw.purity_from_relative(z, e) - st)
            ps = float(np.mean(ps))
            N = 0.5 * (e0.size + e1.size) / nu
            m = mr[float(t)]
            print(f"  {lab:26s} {t}   {gx - m['gpurity']:+.4f}          {gs - m['gpurity']:+.4f}       |  "
                  f"{px - m['purity']:+.4f}          {ps - m['purity']:+.4f}      | {N:5.2f}", flush=True)
            out[lab].append(dict(t=t, Pa_cross=gx, Pa_self=gs, P1_cross=px, P1_self=ps, N=N,
                                 Pa_mesh=m["gpurity"], P1_mesh=m["purity"]))
    return out


# ----------------------------------------------------------------------------
# Part F -- figure
# ----------------------------------------------------------------------------
FIG_SERIES = [   # (run key, legend text, colour slot) -- a fixed order, never cycled
    ("contact_c0.7_nu100000", "contact sink, c = 0.7", "#eb6834"),
    ("contact_c0.35_nu1000000", "contact sink, c = 0.35", "#eda100"),
    ("none_c0.35_absP0.25_nu100000", "absorption only, δ = 0.25", "#1baf7a"),
    ("K1s10_c1.0_nu1000000", "K1, c = 1.0", "#6250d6"),
    ("K1s10_c0.7_nu1000000", "K1, c = 0.7", "#2a78d6"),
    ("K1s10_c0.5_nu1000000", "K1, c = 0.5", "#e87ba4"),
]


FIG_LABELS = {   # the points of panel (c) that carry a label; the rest are identified in the note's table
    "contact_c0.7": "contact 0.7", "contact_c0.35": "contact 0.35", "contacts_c0.7": "contact 0.7, shift every step",
    "K1_c0.7": "K1 at cell centre, 0.7", "K1b_c0.7": "K1, shift every step, 0.7", "K1s10_c0.7": "K1s10, 0.7",
    "K1s10_c1.0": "K1s10, 1.0", "K1s10_c0.5": "K1s10, 0.5", "none_c0.35_absP0.25": "absorption δ = 0.25",
    "none_c0.35_absS0.175": "absorption, single, δ = 0.175",
}


FIG_OFFSETS = {"K1s10_c1.0": (6, 2)}


def part_F(summ, mesh):
    """(a) P_a - mesh against t, (b) population, (c) deficit at t = 3 against cost
    for every rule (cross-seed estimator), from popctl_C_nu*.json."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ink, muted = "#2b2b2a", "#8a8a85"
    plt.rcParams.update({"axes.edgecolor": muted, "axes.labelcolor": ink, "xtick.color": muted,
                         "ytick.color": muted, "axes.spines.top": False, "axes.spines.right": False,
                         "font.size": 10})
    fig, ax = plt.subplots(1, 3, figsize=(16, 5))
    ends = ([], [])
    for key, name, col in FIG_SERIES:
        if key not in summ:
            continue
        rows = summ[key]
        t = [r["t"] for r in rows]
        nu = "10$^6$" if key.endswith("nu1000000") else "10$^5$"
        ax[0].plot(t, [r["dPa"] for r in rows], color=col, lw=2, label=f"{name}, $\\nu$ = {nu}")
        ends[0].append((rows[-1]["dPa"], name.split(",")[0] + name.split(",")[1], t[-1]))
        ax[1].plot(t, [r["N"] for r in rows], color=col, lw=2)
        ends[1].append((rows[-1]["N"], name.split(",")[0] + name.split(",")[1], t[-1]))
    for a, e, gap in ((ax[0], ends[0], 0.011), (ax[1], ends[1], 0.13)):
        e.sort()
        ys = [y for y, _, _ in e]
        for i in range(1, len(ys)):          # push end labels apart, bottom up
            ys[i] = max(ys[i], ys[i - 1] + gap)
        for (y, txt, tx), yl in zip(e, ys):
            a.annotate(txt, (tx, y), xytext=(tx + 0.08, yl), textcoords="data", fontsize=8, color=ink, va="center")
    for a in ax[:2]:
        a.set_xlim(0, 4.3)
        a.grid(True, color="#e6e6e3", lw=0.6)
    ax[0].axhline(0, color=ink, lw=0.8)
    ax[0].set(xlabel="t", ylabel=r"$P_{0.25}$(worlds) $-$ $P_{0.25}$(mesh)",
              title="(a) global purity against the mesh, reach $2\\pi$")
    ax[0].legend(fontsize=8, loc="lower left", frameon=False)
    ax[1].set(xlabel="t", ylabel=r"worlds per $\nu$,  $N/\nu$", title="(b) population")
    # (c) every rule at t = 3: cross-seed deficit against population
    fam = {"contact": ("#eb6834", "o"), "K1": ("#2a78d6", "s"), "none": ("#1baf7a", "^")}
    for nu, filled in ((100000, False), (1000000, True)):
        c = tw.load(f"popctl_C_nu{nu}.json") or {}
        for lab, rows in sorted(c.items()):
            r3 = [r for r in rows if r["t"] == 3]
            if not r3:
                continue
            col, mk = fam["K1" if lab.startswith("K1") else ("contact" if lab.startswith("contact") else "none")]
            d = r3[0]["Pa_cross"] - r3[0]["Pa_mesh"]
            ax[2].plot(r3[0]["N"], d, mk, ms=9 if filled else 8, color=col,
                       mfc=col if filled else "white", mew=2)
            short = FIG_LABELS.get(lab)
            if short and filled:                 # filled points: a short label beside them
                above = lab == "K1s10_c0.5"
                ax[2].annotate(short.replace("K1s10, ", "") + " (10$^6$)", (r3[0]["N"], d),
                               xytext=(0, 9) if above else (-8, 0), textcoords="offset points", fontsize=7.5,
                               color=ink, ha="center" if above else "right", va="bottom" if above else "center")
            elif short:
                ax[2].annotate(short, (r3[0]["N"], d), xytext=FIG_OFFSETS.get(lab, (6, -3)),
                               textcoords="offset points", fontsize=7.5, color=ink)
    ax[2].axhline(0, color=ink, lw=0.8)
    ax[2].grid(True, color="#e6e6e3", lw=0.6)
    ax[2].set_xlim(-0.6, 8.6)
    ax[2].set_xticks(range(0, 9))
    for lab, (col, mk) in (("contact sink", fam["contact"]), ("K1 resampling", fam["K1"]),
                           ("absorption", fam["none"])):
        ax[2].plot([], [], mk, color=col, mfc="white", mew=2, label=lab)
    ax[2].legend(fontsize=8, frameon=False, loc="center right", bbox_to_anchor=(1.0, 0.6),
                 title="open: $\\nu$ = 10$^5$; filled: 10$^6$", title_fontsize=8)
    ax[2].set(xlabel=r"$N/\nu$ at t = 3", ylabel=r"deficit at t = 3 (cross-seed estimator)",
              title="(c) bias against cost, every rule")
    fig.tight_layout()
    fig.savefig(output_path("population_control.png"), dpi=130)
    dp = docs_path("population_control.png")
    if dp:
        fig.savefig(dp, dpi=130)
    print("  wrote population_control.png")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("parts", nargs="*", default=["T", "P", "R", "F"])
    ap.add_argument("--config", nargs="*", default=None,
                    help="sink:cell[:delta[:P|S|D|C|W]] e.g. K1:0.35 contact:0.7:0.1:P K1s10r0.25:0.7:0.7:C")
    ap.add_argument("--nu", type=int, default=100000)
    ap.add_argument("--seeds", type=int, default=2)
    ap.add_argument("--T", type=float, default=tw.T_END)
    ap.add_argument("--cap", type=float, default=12.0, help="stop a run when N/nu exceeds this")
    ap.add_argument("--workers", type=int, default=os.cpu_count() or 1)
    args = ap.parse_args()
    tw.set_particles(3)
    parts = [p.upper() for p in args.parts]
    t0 = time.perf_counter()
    if "T" in parts:
        tw.dump(part_T(), "popctl_T.json")
    if "P" in parts:
        tw.dump(part_P(), "popctl_P.json")
    if "R" in parts:
        cfg = []
        for c in args.config or ["contact:0.35", "K1:0.35"]:
            f = c.split(":")
            cfg.append((f[0], float(f[1]), float(f[2]) if len(f) > 2 else 0.0, f[3] if len(f) > 3 else "P"))
        part_R(cfg, args.nu, args.seeds, args.T, args.workers, args.cap)
    if not os.path.exists(output_path("three_worlds_M2.json")):
        print("\n(no three_worlds_M2.json here: run demo_three_worlds.py M --mesh-reach 2 for the summaries)")
        print(f"\nwall time {time.perf_counter() - t0:.1f} s")
        return
    mesh = tw.load("three_worlds_M2.json")
    if "C" in parts:
        tw.dump(part_C(mesh, args.nu), f"popctl_C_nu{args.nu}.json")
    if "R" in parts or "S" in parts or "F" in parts:
        summ = summary(mesh)
        tw.dump(summ, "popctl_summary.json")
        if "F" in parts:
            part_F(summ, mesh)
    print(f"\nwall time {time.perf_counter() - t0:.1f} s")


if __name__ == "__main__":
    main()
