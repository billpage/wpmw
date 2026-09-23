#!/usr/bin/env python3
"""
Companion demo for ``docs/analysis/fourd_compensated_ledger.md`` (step 21).

The compensated world form, with its ledger on exact integer counts, run in
the 4D phase space of two particles on a line, with the bin area A in the
centre-of-mass directions scanned.

The problem
-----------
Two unit-mass particles (hbar = 1).  Centre of mass X = (x1 + x2)/2,
P = p1 + p2, mass M = 2; relative r = x1 - x2, p_r = (p1 - p2)/2, reduced
mass mu = 1/2.  Both maps are canonical with unit Jacobian.

    H = P^2/2M + M Omega^2 X^2 / 2  +  p_r^2 / 2 mu - V0 sech^2(r/a)

The CM trap is quadratic, so its residual is identically zero (G2) and the
CM is carried by Newtonian streaming alone.  The pair potential is the
attractive Poeschl-Teller well with V0 = 6, a = 1, for which
lambda(lambda+1) = 2 mu V0 a^2 / hbar^2 = 6, lambda = 2, with exact bound
states  psi0 ~ sech^2(r)  (E0 = -4)  and  psi1 ~ sech(r) tanh(r)  (E1 = -1).
Every state used below is stationary or exactly known, so the reference is
exact, and a stationary state is what lets the ledger reach equilibrium.

The world form
--------------
Bodies are signed points (X, P, r, p_r).  Postulate (S): every body streams
on its Newtonian arc -- exact rotation in the CM trap, velocity Verlet in r
under the compensated force.  Postulate (D): each body parents residual
events at rate Gamma(r) = sum_q |K_q(r)|; an event deposits +1 at p_r + xi_q
and -1 at p_r - xi_q (times sgn K_q times the parent's sign).  It is realised
absorptively if a partner of the right species sits in the parent's cell at
BOTH daughter rows, otherwise emissively (spec §5.3-5.5).  Counts are
integers, so the S7 identity  dN = 2(1 - 2f) n_ev  holds exactly.

A cell is (r-bin, p_r row, X-bin, P-bin).  The r-bin and p_r row are the
mesh of the 1D ledger runs (the p_r row is forced by the reach, C4).  The
(X, P) bins have area A, aspect sigma_X / sigma_P, and are not forced by
anything: A is the scanned parameter.  A = inf is one bin, i.e. the 1D ledger
with the CM carried as a passive label.

The reference is the mesh QLE of the same compensated symbol on the same
(r, p_r) mesh, times the exact CM Wigner function.  Errors quoted are the
world form against that reference, beside the t = 0 sampling floor.

Parts (run any subset; default all, about five minutes)
-------------------------------------------------------
E  The Eckart pair collision as an exact 4D reference: 2D split-operator
   against the 1D relative run and the closed form; entanglement (§1.2).
K  The residual kernel against the uncompensated one on the same mesh:
   event budget and first moment, for three potentials (Proposition M4).
R  Mesh-QLE drift of the two exact eigenstates (§1.3).
W  Sub-cell positions: <V> bias as the matching cell is refined, with a
   classical-streaming control (Proposition M9).
S  The scan over A and nbar on the stationary product state
   W_X(thermal, nbar) x W_r(psi0) (Propositions M5, M7).
X  The same at A = inf with the r-window doubled (Proposition M6).
B  Two CM components carrying psi0 and psi1: conditional <V> at A = inf
   and A = h, and a five-seed table (Proposition M8).
F  The figure, from the JSON that S, X, B and W save.

Run:  WPMW_OUTPUT=... PYTHONPATH=src python3 src/demo_fourd_compensated_ledger.py [E K R W S X B F]
"""
import sys
import time
import json

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from wpmwlib.wpmw_utils import docs_path, output_path

HBAR, MU, MTOT, OMEGA = 1.0, 0.5, 2.0, 1.0
V0, A_W = 6.0, 1.0
NPR, DP = 64, 0.25
H = 2 * np.pi * HBAR
DT = 0.01
N_CAP = 1_500_000


def V(r):
    return -V0 / np.cosh(r / A_W) ** 2


# ======================================================================
# The relative-coordinate mesh: compensated symbol, kernel, reference
# ======================================================================
class RelMesh:
    """Symbols built as in demo_sea_population_equilibrium.Ledger (§3.2,
    §4.1 normative: Nyquist rung zeroed, raised-cosine horizon applied first,
    force from the kernel's own discrete first moment), with mass mu."""

    def __init__(self, nr=128, rh=10.0, pot=None):
        self.nr, self.rh = nr, rh
        V_ = V if pot is None else pot
        NR, RH = nr, rh
        self.r = -RH + 2 * RH * np.arange(NR) / NR
        self.dr = self.r[1] - self.r[0]
        self.kr = 2 * np.pi * np.fft.fftfreq(NR, d=self.dr)
        self.n = np.round(np.fft.fftfreq(NPR, d=1.0 / NPR)).astype(int)
        self.p = DP * (self.n + 0.5)
        s = 2 * np.pi * np.fft.fftfreq(NPR, d=DP)
        y = HBAR * s / 2
        self.y_max = float(np.abs(y).max())
        xi = self.n * DP
        rr, yy = self.r[:, None], y[None, :]
        m_full = (1j / HBAR) * (V_(rr + yy) - V_(rr - yy))
        s_sym = np.broadcast_to(s, m_full.shape).copy()
        m_full[:, NPR // 2] = 0
        s_sym[:, NPR // 2] = 0
        w = np.cos(np.pi * y / (2 * self.y_max)) ** 2
        m_full, s_sym = m_full * w, s_sym * w

        def fm(sym):
            return np.real((xi * np.fft.ifft(sym, axis=1)).sum(axis=1))

        self.dv = fm(m_full) / fm(1j * s_sym)
        m_res = m_full - 1j * self.dv[:, None] * s_sym
        self.k = np.real(np.fft.ifft(m_res, axis=1))   # k[ir, column q]
        self.k_full = np.real(np.fft.ifft(m_full, axis=1))   # uncompensated
        self.xi = xi
        self.sym_e = np.fft.fft(self.k, axis=1)
        self.cls = 1j * self.dv[:, None] * s[None, :]
        # one event per (q, -q): columns 1 .. NPR/2 - 1
        kk = np.abs(self.k[:, 1:NPR // 2])
        self.gamma = kk.sum(axis=1)
        self.cum = np.cumsum(kk, axis=1)
        self.sgn = np.sign(self.k[:, 1:NPR // 2])

    def stream(self, f, dt):
        fh = np.fft.fft(f, axis=0) * np.exp(-1j * self.kr[:, None]
                                            * self.p[None, :] * dt / MU)
        f = np.real(np.fft.ifft(fh, axis=0))
        return np.real(np.fft.ifft(np.fft.fft(f, axis=1)
                                   * np.exp(dt * self.cls), axis=1))

    def qle_step(self, e, dt):
        e = self.stream(e, 0.5 * dt)
        e = np.real(np.fft.ifft(np.fft.fft(e, axis=1)
                                * np.exp(dt * self.sym_e), axis=1))
        return self.stream(e, 0.5 * dt)

    def wigner(self, psi_fn):
        """W(r, p) at mesh points, normalised on the mesh (fft row order)."""
        yq = np.arange(-15, 15, 0.004)
        W = np.empty((self.nr, NPR))
        for i, r0 in enumerate(self.r):
            prod = psi_fn(r0 + yq) * psi_fn(r0 - yq)
            W[i] = (np.cos(2 * np.outer(self.p, yq) / HBAR) @ prod) \
                * 0.004 / (np.pi * HBAR)
        return W / (W.sum() * self.dr * DP)

    def force(self, r):
        """-dv_eff(r) by periodic linear interpolation."""
        u = (r + self.rh) / self.dr
        i0 = np.floor(u).astype(int)
        t = u - i0
        n = self.nr
        return -((1 - t) * self.dv[i0 % n] + t * self.dv[(i0 + 1) % n])


def psi0(r):
    return 1 / np.cosh(r) ** 2


def psi1(r):
    return np.tanh(r) / np.cosh(r)


# ======================================================================
# The world ensemble
# ======================================================================
class Worlds:
    OFF = 1 << 15

    def __init__(self, mesh, area, sig_x, sig_p, rng):
        self.m, self.rng = mesh, rng
        if np.isinf(area):
            self.dX = self.dP = np.inf
        else:
            self.dX = np.sqrt(area * sig_x / sig_p)
            self.dP = np.sqrt(area * sig_p / sig_x)
        self.X = self.P = self.r = self.pr = self.sp = None
        self.n_abs = self.n_emi = 0
        self.occ = []                     # (countA, countB) per event
        self.absorb = True

    # -- initial sampling ------------------------------------------------
    def sample(self, comps, nu, rho):
        """comps: list of (weight, W_r mesh, (Xc, Pc, sX, sP)).  Signed
        bodies from |W_r| x W_X, plus (rho - 1)/2 co-located +- pairs per
        body, which leave E untouched in every cell under any binning."""
        m = self.m
        parts = []
        for wgt, Wr, (xc, pc, sx, sp_) in comps:
            l1 = np.abs(Wr).sum() * m.dr * DP
            n0 = int(round(nu * wgt * l1))
            prob = np.abs(Wr).ravel() / np.abs(Wr).sum()
            cells = self.rng.choice(prob.size, size=n0, p=prob)
            ir, jc = np.unravel_index(cells, Wr.shape)
            r = m.r[ir] + (self.rng.random(n0) - 0.5) * m.dr
            pr = DP * (m.n[jc] + self.rng.random(n0))
            sgn = np.sign(Wr[ir, jc])
            X = xc + sx * self.rng.standard_normal(n0)
            P = pc + sp_ * self.rng.standard_normal(n0)
            npad = int(round(0.5 * (rho - 1) * n0))
            k = self.rng.integers(0, n0, npad)
            parts.append((np.concatenate([X, X[k], X[k]]),
                          np.concatenate([P, P[k], P[k]]),
                          np.concatenate([r, r[k], r[k]]),
                          np.concatenate([pr, pr[k], pr[k]]),
                          np.concatenate([sgn, np.ones(npad), -np.ones(npad)])))
        self.X, self.P, self.r, self.pr, self.sp = (
            np.concatenate([p[i] for p in parts]) for i in range(5))
        self.nu = nu

    # -- postulate (S) ---------------------------------------------------
    def stream(self, dt):
        c, s = np.cos(OMEGA * dt), np.sin(OMEGA * dt)
        X, P = self.X, self.P
        self.X, self.P = (X * c + P / (MTOT * OMEGA) * s,
                          P * c - MTOT * OMEGA * X * s)
        self.pr += 0.5 * dt * self.m.force(self.r)
        self.r += self.pr / MU * dt
        rh = self.m.rh
        self.r = (self.r + rh) % (2 * rh) - rh
        self.pr += 0.5 * dt * self.m.force(self.r)
        # the mesh reference is periodic in p (FFT); wrap to match it
        pw = NPR * DP
        self.pr = (self.pr + pw / 2) % pw - pw / 2

    # -- cells -----------------------------------------------------------
    SUB_R = 1      # matching bins per mesh r-bin
    SUB_P = 1      # matching bins per momentum row

    def indices(self):
        """Kernel index ir, momentum row, and the MATCHING cell (irm, psub,
        iX, iP).  The kernel lives on the mesh; the cell in which partners
        are sought may be finer, in r and within a row, as well as in CM."""
        nr = self.m.nr
        u = (self.r + self.m.rh) / self.m.dr + 0.5
        ir = np.floor(u).astype(int) % nr
        irm = np.floor(u * self.SUB_R).astype(np.int64) % (nr * self.SUB_R)
        v = self.pr / DP
        row = np.floor(v).astype(int) % NPR                # fft column
        psub = np.floor((v - np.floor(v)) * self.SUB_P).astype(np.int64)
        if np.isinf(self.dX):
            iX = iP = np.zeros(self.r.size, dtype=np.int64)
        else:
            iX = np.floor(self.X / self.dX).astype(np.int64)
            iP = np.floor(self.P / self.dP).astype(np.int64)
        return ir, row, iX, iP, irm, psub

    def key(self, irm, row, psub, iX, iP, sp01):
        o = self.OFF
        c = (irm * NPR + row) * self.SUB_P + psub
        return ((c * (2 * o) + iX + o) * (2 * o) + iP + o) * 2 + sp01

    # -- postulate (D) ---------------------------------------------------
    def events(self, dt):
        m, rng = self.m, self.rng
        ir, row, iX, iP, irm, psub = self.indices()
        keys = self.key(irm, row, psub, iX, iP, (self.sp > 0).astype(np.int64))
        order = np.argsort(keys, kind="stable")
        sk = keys[order]
        brk = np.flatnonzero(np.diff(sk)) + 1
        start = np.concatenate([[0], brk])
        end = np.concatenate([brk, [sk.size]])
        uniq = sk[start]
        ptr = start.copy()

        # event list: Poisson per parent, channel by |K_q|
        lam = m.gamma[ir] * dt
        nev = rng.poisson(lam)
        par = np.repeat(np.arange(nev.size), nev)
        if par.size == 0:
            return
        rng.shuffle(par)
        u = rng.random(par.size) * m.gamma[ir[par]]
        q = 1 + (m.cum[ir[par]] < u[:, None]).sum(axis=1)
        q = np.minimum(q, NPR // 2 - 1)
        t = (m.sgn[ir[par], q - 1] * self.sp[par]).astype(np.int64)
        rA = (row[par] + q) % NPR
        rB = (row[par] - q) % NPR
        # absorptive partners: species -t at A, species +t at B
        kA = self.key(irm[par], rA, psub[par], iX[par], iP[par],
                      (t < 0).astype(np.int64))
        kB = self.key(irm[par], rB, psub[par], iX[par], iP[par],
                      (t > 0).astype(np.int64))
        gA = np.searchsorted(uniq, kA)
        gA[(gA >= uniq.size)] = uniq.size - 1
        gA = np.where(uniq[gA] == kA, gA, -1)
        gB = np.searchsorted(uniq, kB)
        gB[(gB >= uniq.size)] = uniq.size - 1
        gB = np.where(uniq[gB] == kB, gB, -1)
        xi = q * DP

        consumed, newd = [], {}
        new_X, new_P, new_r, new_pr, new_sp = [], [], [], [], []
        new_alive = []
        occA = np.empty(par.size, dtype=np.int32)
        occB = np.empty(par.size, dtype=np.int32)
        ka, kb, ga, gb = kA.tolist(), kB.tolist(), gA.tolist(), gB.tolist()
        tl, pl, xil = t.tolist(), par.tolist(), xi.tolist()
        n_abs = n_emi = 0
        for e in range(par.size):
            a, b = ga[e], gb[e]
            la = newd.get(ka[e])
            lb = newd.get(kb[e])
            ca = (end[a] - ptr[a] if a >= 0 else 0) + (len(la) if la else 0)
            cb = (end[b] - ptr[b] if b >= 0 else 0) + (len(lb) if lb else 0)
            occA[e], occB[e] = ca, cb
            if ca > 0 and cb > 0 and self.absorb:
                n_abs += 1
                if a >= 0 and ptr[a] < end[a]:
                    consumed.append(order[ptr[a]])
                    ptr[a] += 1
                else:
                    new_alive[la.pop()] = False
                if b >= 0 and ptr[b] < end[b]:
                    consumed.append(order[ptr[b]])
                    ptr[b] += 1
                else:
                    new_alive[lb.pop()] = False
            else:
                n_emi += 1
                j = pl[e]
                for kk, dsign in ((ka[e] ^ 1, 1.0), (kb[e] ^ 1, -1.0)):
                    idx = len(new_alive)
                    new_alive.append(True)
                    new_X.append(self.X[j])
                    new_P.append(self.P[j])
                    new_r.append(self.r[j])
                    new_pr.append(self.pr[j] + dsign * xil[e])
                    new_sp.append(float(tl[e]) * dsign)
                    newd.setdefault(kk, []).append(idx)
        self.n_abs += n_abs
        self.n_emi += n_emi
        self.occ.append((occA, occB))
        keep = np.ones(self.r.size, dtype=bool)
        keep[np.array(consumed, dtype=np.int64)] = False
        na = np.array(new_alive, dtype=bool)

        def merge(old, new):
            return np.concatenate([old[keep], np.array(new)[na]]) \
                if len(new) else old[keep]

        self.X, self.P = merge(self.X, new_X), merge(self.P, new_P)
        self.r, self.pr = merge(self.r, new_r), merge(self.pr, new_pr)
        self.sp = merge(self.sp, new_sp)

    # -- observables -----------------------------------------------------
    def hist_rel(self, mask=None):
        ir, row = self.indices()[:2]
        w = self.sp if mask is None else self.sp * mask
        h = np.zeros(self.m.nr * NPR)
        np.add.at(h, ir * NPR + row, w)
        return h.reshape(self.m.nr, NPR) / (self.nu * self.m.dr * DP)


def rel_l2(a, b):
    return float(np.linalg.norm(a - b) / np.linalg.norm(b))


# ======================================================================
# Experiment S -- stationary product state
# ======================================================================
def run_S(mesh, W0, area, nbar, nu=4000, rho=3.0, T=4.0, seed=1,
          absorb=True):
    sx = np.sqrt((2 * nbar + 1) * HBAR / (2 * MTOT * OMEGA))
    spp = np.sqrt((2 * nbar + 1) * HBAR * MTOT * OMEGA / 2)
    rng = np.random.default_rng(seed)
    w = Worlds(mesh, area, sx, spp, rng)
    w.absorb = absorb
    w.sample([(1.0, W0, (0.0, 0.0, sx, spp))], nu, rho)
    ref = W0.copy()
    rec, t0 = [], time.time()
    nsteps = int(round(T / DT))
    capped = False

    def sectors():
        th = np.arctan2(w.P / spp, w.X / sx)
        return [((th >= a) & (th < a + np.pi / 2)).astype(float)
                for a in (-np.pi, -np.pi / 2, 0, np.pi / 2)]

    def record(step):
        e = w.hist_rel()
        sec = [rel_l2(4 * w.hist_rel(mk), ref) for mk in sectors()]
        occ = np.concatenate([np.stack(o) for o in w.occ[-10:]], axis=1) \
            if w.occ else np.zeros((2, 1))
        ir, row, iX, iP, irm, psub = w.indices()
        cells = np.unique(w.key(irm, row, psub, iX, iP,
                                np.zeros_like(irm))).size
        # signed energy moments of the relative motion: robust to noise
        kin = float(np.sum(w.sp * w.pr ** 2) / (2 * MU) / w.nu)
        pot = float(np.sum(w.sp * V(w.r)) / w.nu)
        pm = mesh.p[None, :]
        wsum = ref.sum() * mesh.dr * DP
        kin_ref = float((ref * pm ** 2).sum() * mesh.dr * DP / (2 * MU) / wsum)
        pot_ref = float((ref * V(mesh.r)[:, None]).sum() * mesh.dr * DP / wsum)
        irk = w.indices()[0]
        quiet = float((mesh.gamma[irk] < 0.01 * mesh.gamma.max()).mean())
        rec.append(dict(t=step * DT, N=int(w.r.size), quiet=quiet,
                        kin=kin, pot=pot, kin_ref=kin_ref, pot_ref=pot_ref,
                        nplus=int((w.sp > 0).sum()),
                        err=rel_l2(e, ref), err_sec=float(np.mean(sec)),
                        n_abs=w.n_abs, n_emi=w.n_emi,
                        lamA=float(occ[0].mean()), lamB=float(occ[1].mean()),
                        availA=float((occ[0] > 0).mean()),
                        availB=float((occ[1] > 0).mean()),
                        both=float(((occ[0] > 0) & (occ[1] > 0)).mean()),
                        cells=int(cells)))

    record(0)
    n_prev = 0
    for step in range(1, nsteps + 1):
        w.stream(0.5 * DT)
        w.events(DT)
        w.stream(0.5 * DT)
        ref = mesh.qle_step(ref, DT)
        if step % 25 == 0:
            record(step)
            w.occ = []
        if w.r.size > N_CAP:
            capped = True
            record(step)
            break
    n_ev = w.n_abs + w.n_emi
    n_start = rec[0]["N"]
    return dict(area=area, nbar=nbar, nu=nu, rho=rho, rec=rec,
                capped=capped, secs=time.time() - t0,
                s7=(rec[-1]["N"] - n_start) - 2 * (w.n_emi - w.n_abs),
                n_ev=n_ev)


# ======================================================================
# Experiment B -- blur: two CM components carrying different relative states
# ======================================================================
def run_B(mesh, W0, W1, area, d=2.0, nu=4000, rho=3.0, T=np.pi, seed=2):
    """Returns the record; pot0/pot1 are the signed <V(r)> of the bodies on
    each component's side of the CM bisector, beside the mesh reference."""
    sx = np.sqrt(HBAR / (2 * MTOT * OMEGA))
    spp = np.sqrt(HBAR * MTOT * OMEGA / 2)
    rng = np.random.default_rng(seed)
    w = Worlds(mesh, area, sx, spp, rng)
    w.sample([(0.5, W0, (+d, 0.0, sx, spp)), (0.5, W1, (-d, 0.0, sx, spp))],
             nu, rho)
    ref0, ref1 = W0.copy(), W1.copy()
    nsteps = int(round(T / DT))
    rec = []

    def record(step):
        tt = step * DT
        # exact CM centres rotate in the trap; split by the bisector
        cx, cp = d * np.cos(OMEGA * tt), -MTOT * OMEGA * d * np.sin(OMEGA * tt)
        proj = (w.X / sx) * (cx / sx) + (w.P / spp) * (cp / spp)
        m0 = (proj > 0).astype(float)
        e0 = 2 * w.hist_rel(m0)
        e1 = 2 * w.hist_rel(1 - m0)
        vr = V(w.r)
        pot0 = float(2 * np.sum(w.sp * m0 * vr) / w.nu)
        pot1 = float(2 * np.sum(w.sp * (1 - m0) * vr) / w.nu)
        vm = V(mesh.r)[:, None]
        pr0 = float((ref0 * vm).sum() / ref0.sum())
        pr1 = float((ref1 * vm).sum() / ref1.sum())
        rec.append(dict(t=tt, N=int(w.r.size), pot0=pot0, pot1=pot1,
                        pot0_ref=pr0, pot1_ref=pr1, err0=rel_l2(e0, ref0),
                        err1=rel_l2(e1, ref1),
                        mix0=float(np.sum(e0 * ref1) / np.sum(ref1 * ref1)),
                        n_abs=w.n_abs, n_emi=w.n_emi))

    record(0)
    for step in range(1, nsteps + 1):
        w.stream(0.5 * DT)
        w.events(DT)
        w.stream(0.5 * DT)
        ref0, ref1 = mesh.qle_step(ref0, DT), mesh.qle_step(ref1, DT)
        w.occ = []
        if step % 25 == 0 or step == nsteps:
            record(step)
        if w.r.size > N_CAP:
            record(step)
            break
    return dict(area=area, d=d, rec=rec)


# ======================================================================
# Part E -- the Eckart pair collision as an exact 4D reference (2D grid)
# ======================================================================
def part_E():
    """Two particles, U = V0 sech^2((x1 - x2)/a) with V0 = 1, a = 2, from a
    product of Gaussians; 2D split-operator in particle coordinates against
    the 1D relative-coordinate run and the closed-form transmission."""
    v0, a, sig_r, r_c, p_c, T, dt = 1.0, 2.0, 2.0, -16.0, 1.0, 22.0, 0.01
    U = lambda r: v0 / np.cosh(r / a) ** 2

    def t_eck(e):
        k = np.sqrt(np.maximum(2 * MU * e, 0)) / HBAR
        lam = 8 * MU * v0 * a**2 / HBAR**2
        c = np.cosh(0.5 * np.pi * np.sqrt(lam - 1)) ** 2
        sh = np.sinh(np.pi * k * a) ** 2
        return sh / (sh + c)

    sp = HBAR / (2 * sig_r)
    pp = np.linspace(1e-6, p_c + 12 * sp, 400001)
    w = np.exp(-(pp - p_c) ** 2 / (2 * sp**2))
    w /= np.trapezoid(w, pp)
    tc = float(np.trapezoid(w * t_eck(pp**2 / (2 * MU)), pp))

    n, half = 512, 48.0
    x = (np.arange(n) - n // 2) * (2 * half / n)
    dx = x[1] - x[0]
    k = 2 * np.pi * np.fft.fftfreq(n, d=dx)
    X1, X2 = np.meshgrid(x, x, indexing="ij")
    K1, K2 = np.meshgrid(k, k, indexing="ij")
    s1 = sig_r / np.sqrt(2)
    psi = (np.exp(-(X1 - r_c / 2) ** 2 / (4 * s1**2) + 1j * p_c * X1 / HBAR)
           * np.exp(-(X2 + r_c / 2) ** 2 / (4 * s1**2) - 1j * p_c * X2 / HBAR))
    psi /= np.sqrt(np.sum(abs(psi) ** 2))
    ev = np.exp(-0.5j * U(X1 - X2) * dt / HBAR)
    et = np.exp(-1j * HBAR * (K1**2 + K2**2) / 2 * dt)
    for _ in range(int(round(T / dt))):
        psi = ev * np.fft.ifft2(et * np.fft.fft2(ev * psi))
    rho2 = abs(psi) ** 2
    i, j = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
    diff = i - j
    T2 = rho2[diff > 0].sum() + 0.5 * rho2[diff == 0].sum()
    rel2 = np.bincount((diff + n - 1).ravel(), weights=rho2.ravel(),
                       minlength=2 * n - 1)
    r2 = (np.arange(2 * n - 1) - (n - 1)) * dx

    n1 = 2 * n
    r1 = (np.arange(n1) - n1 // 2) * dx
    k1 = 2 * np.pi * np.fft.fftfreq(n1, d=dx)
    chi = np.exp(-(r1 - r_c) ** 2 / (4 * sig_r**2) + 1j * p_c * r1 / HBAR)
    chi /= np.sqrt(np.sum(abs(chi) ** 2))
    e1 = np.exp(-0.5j * U(r1) * dt / HBAR)
    t1 = np.exp(-1j * (HBAR * k1) ** 2 / (2 * MU) * dt / HBAR)
    for _ in range(int(round(T / dt))):
        chi = e1 * np.fft.ifft(t1 * np.fft.fft(e1 * chi))
    rel1 = abs(chi) ** 2
    T1 = rel1[r1 > 0].sum() + 0.5 * rel1[r1 == 0].sum()
    l1 = np.abs(rel2 - np.interp(r2, r1, rel1)).sum()
    sv = np.linalg.svd(psi, compute_uv=False) ** 2
    sv /= sv.sum()
    ent = float(-np.sum(sv[sv > 1e-15] * np.log(sv[sv > 1e-15])))
    print("Part E. Eckart pair collision: an exact 4D reference")
    print(f"  closed form, packet-averaged T(E)   {tc:.6f}")
    print(f"  1D relative split-operator          {T1:.6f}   |1D - cf| = {abs(T1 - tc):.1e}")
    print(f"  2D particle-coordinate run          {T2:.6f}   |2D - 1D| = {abs(T2 - T1):.1e}")
    print(f"  relative density, 2D vs 1D: L1 distance {l1:.1e}")
    print(f"  entanglement of the out-state: purity of rho_1 {np.sum(sv**2):.4f}, "
          f"entropy {ent:.4f} nats (t = 0: 0)")
    return dict(tc=tc, T1=T1, T2=T2, l1=l1, ent=ent)


# ======================================================================
# Part K -- the residual kernel against the uncompensated one
# ======================================================================
def part_K():
    print("Part K. Uncompensated and residual kernels on the same mesh "
          f"(dp = {DP}, y_max = 2 pi)")
    cases = (("Poeschl-Teller -6 sech^2 r", V),
             ("cosine 1.5 cos(pi r / 5)", lambda r: 1.5 * np.cos(np.pi * r / 5)),
             ("harmonic r^2 / 2", lambda r: 0.5 * r**2))
    out = {}
    h = NPR // 2
    for name, pot in cases:
        m = RelMesh(pot=pot)
        core = np.abs(m.r) < 5
        gf = np.abs(m.k_full[:, 1:h]).sum(axis=1)
        gr = np.abs(m.k[:, 1:h]).sum(axis=1)

        def fm(k, j):
            return float(np.sum(m.xi[1:h] * (k[j, 1:h] - k[j, -1:-h:-1])))

        j = int(np.argmax(gf * core))
        force = -np.gradient(pot(m.r), m.r)[j]
        print(f"  {name:28s} max Gamma: full {gf[core].max():6.3f}  "
              f"residual {gr[core].max():6.3f}  ratio {gr[core].max() / gf[core].max():.3f}")
        print(f"  {'':28s} first moment at r = {m.r[j]:+.2f}: full {fm(m.k_full, j):+.4f}"
              f" (-V' = {force:+.4f}), residual {fm(m.k, j):+.1e}")
        out[name] = (gf[core].max(), gr[core].max())
    return out


# ======================================================================
# Part W -- the within-cell trade-off in the 1D ledger (A = inf)
# ======================================================================
def part_W(mesh, W0, seeds=6):
    print("Part W. Sub-cell positions: <V> bias at t = 1 (A = inf, nu = 16000, "
          f"{seeds} seeds)")
    rows = []
    # control: classical streaming only, particles against the mesh
    g, cum, qs = mesh.gamma.copy(), mesh.cum.copy(), mesh.qle_step
    mesh.gamma[:] = 0
    mesh.qle_step = lambda e, dt: mesh.stream(e, dt)
    res = np.array([[r["pot"] - r["pot_ref"] for r in
                     [run_S(mesh, W0, np.inf, 0.0, nu=16000, T=1.0, seed=sd,
                            rho=1.0)["rec"][-1]]] for sd in range(seeds)])
    mesh.gamma[:], mesh.cum[:], mesh.qle_step = g, cum, qs
    print(f"  classical streaming only       dV {res.mean():+.4f} +- "
          f"{res.std(ddof=1) / np.sqrt(seeds):.4f}")
    rows.append(dict(sub_r=0, sub_p=0, dV=float(res.mean()),
                     dV_sem=float(res.std(ddof=1) / np.sqrt(seeds)), N=0.0))
    for sr, spp in ((1, 1), (4, 1), (1, 4), (4, 4)):
        Worlds.SUB_R, Worlds.SUB_P = sr, spp
        res = np.array([(r["pot"] - r["pot_ref"], r["N"]) for r in
                        [run_S(mesh, W0, np.inf, 0.0, nu=16000, T=1.0,
                               seed=sd)["rec"][-1] for sd in range(seeds)]])
        dv, se, n = res[:, 0].mean(), res[:, 0].std(ddof=1) / np.sqrt(seeds), res[:, 1].mean()
        print(f"  matching cell r/{sr}, p/{spp}          dV {dv:+.4f} +- {se:.4f}"
              f"   N {n:8.0f}")
        rows.append(dict(sub_r=sr, sub_p=spp, dV=float(dv), dV_sem=float(se),
                         N=float(n)))
    Worlds.SUB_R = Worlds.SUB_P = 1
    return rows


def dump(obj, name):
    with open(output_path(name), "w") as fh:
        json.dump(obj, fh, default=float)


def load(name):
    with open(output_path(name)) as fh:
        return json.load(fh)


# ======================================================================
# The figures
# ======================================================================
def save_fig(fig, name):
    fig.savefig(output_path(name), dpi=150, bbox_inches="tight")
    dp = docs_path(name)
    if dp:
        fig.savefig(dp, dpi=150, bbox_inches="tight")
    plt.close(fig)


def figure():
    S = [("A = inf (the 1D ledger)", "S_inf.json", "k"),
         ("A = h", "S_1.json", "C0"),
         ("A = h/4", "S_0.25.json", "C1"),
         ("A = h/16", "S_0.0625.json", "C3"),
         ("A = h/4, nbar = 3", "S_0.25_n3.json", "C2")]
    S = [(lab, load("fourd_ledger_" + f), c) for lab, f, c in S]
    wide = load("fourd_ledger_X.json")
    B = [("A = inf", load("fourd_ledger_B_inf.json"), "k"),
         ("A = h", load("fourd_ledger_B_1.json"), "C0")]
    WC = load("fourd_ledger_W.json")

    fig, ax = plt.subplots(2, 2, figsize=(11, 8))
    a = ax[0, 0]
    for lab, o, c in S:
        a.semilogy([r["t"] for r in o["rec"]], [r["N"] for r in o["rec"]],
                   c, label=lab)
    a.semilogy([r["t"] for r in wide["rec"]], [r["N"] for r in wide["rec"]],
               "k--", label="A = inf, r-window doubled")
    a.set(xlabel="t", ylabel="bodies N",
          title="(a) body count against CM bin area A")
    a.legend(fontsize=8)

    a = ax[0, 1]
    for lab, o, c in S + [("A = inf, r-window doubled", wide, "k")]:
        a.loglog([r["cells"] for r in o["rec"][1:]],
                 [r["N"] for r in o["rec"][1:]],
                 "o--" if "doubled" in lab else "o-", color=c, ms=3, label=lab)
    x = np.array([1e3, 5e5])
    a.loglog(x, 2.5 * x, ":", color="gray", label="2.5 bodies per cell")
    a.loglog(x, 4.0 * x, "-.", color="gray", label="4 bodies per cell")
    a.set(xlabel="occupied cells", ylabel="bodies N",
          title="(b) the population fills cells")
    a.legend(fontsize=7)

    a = ax[1, 0]
    for lab, o, c in B:
        t = [r["t"] for r in o["rec"]]
        a.plot(t, [r["pot0"] for r in o["rec"]], c,
               label=lab + ": ground-state side")
        a.plot(t, [r["pot1"] for r in o["rec"]], c, ls="--",
               label=lab + ": excited-state side")
    o = B[0][1]
    t = [r["t"] for r in o["rec"]]
    a.plot(t, [r["pot0_ref"] for r in o["rec"]], color="C1", lw=3, alpha=.4,
           label="reference")
    a.plot(t, [r["pot1_ref"] for r in o["rec"]], color="C1", lw=3, alpha=.4)
    a.set(xlabel="t", ylabel=r"conditional $\langle V(r)\rangle$",
          title="(c) two CM components carrying two relative states")
    a.legend(fontsize=7)

    a = ax[1, 1]
    rows = [w for w in WC if w["sub_r"] > 0]
    a.errorbar([w["N"] for w in rows], [w["dV"] for w in rows],
               yerr=[w["dV_sem"] for w in rows], fmt="o", color="C4")
    for w in rows:
        a.annotate(f"r/{w['sub_r']}, p/{w['sub_p']}", (w["N"], w["dV"]),
                   textcoords="offset points", xytext=(6, 4), fontsize=8)
    c0 = WC[0]
    a.axhspan(c0["dV"] - c0["dV_sem"], c0["dV"] + c0["dV_sem"], color="gray",
              alpha=0.3, label="classical streaming only")
    a.axhline(0, color="gray", lw=0.8)
    a.set_xscale("log")
    a.set(xlabel="bodies N at t = 1", ylabel=r"$\langle V\rangle$ bias at t = 1",
          title="(d) the same trade-off inside the 1D cell")
    a.legend(fontsize=8)
    fig.tight_layout()
    save_fig(fig, "fourd_compensated_ledger.png")
    print("figure:", output_path("fourd_compensated_ledger.png"))


# ======================================================================
# Driver
# ======================================================================
SCAN = (("inf", np.inf, 0.0, 4.0), ("1", H, 0.0, 5.0),
        ("0.25", 0.25 * H, 0.0, 4.5), ("0.0625", 0.0625 * H, 0.0, 2.0),
        ("0.25_n3", 0.25 * H, 3.0, 2.0))


def main():
    parts = sys.argv[1:] or ["E", "K", "R", "W", "S", "X", "B", "F"]
    if any(p in parts for p in "RWSXB"):
        mesh = RelMesh()
        W0, W1 = mesh.wigner(psi0), mesh.wigner(psi1)
        print(f"mesh {mesh.nr} x {NPR}, dr = {mesh.dr:.4f}, dp = {DP}, "
              f"cell area h/{H / (mesh.dr * DP):.0f}, y_max = {mesh.y_max:.3f}, "
              f"max Gamma = {mesh.gamma.max():.3f}")
        print(f"||W0||_1 = {np.abs(W0).sum() * mesh.dr * DP:.4f}, "
              f"||W1||_1 = {np.abs(W1).sum() * mesh.dr * DP:.4f}\n")
    for part in parts:
        t0 = time.time()
        if part == "E":
            part_E()
        elif part == "K":
            part_K()
        elif part == "R":
            print("Part R. Mesh QLE drift of the exact eigenstates over t = 4")
            for name, Wx in (("psi0", W0), ("psi1", W1)):
                ref = Wx.copy()
                for _ in range(int(round(4.0 / DT))):
                    ref = mesh.qle_step(ref, DT)
                print(f"  {name}: rel L2 {rel_l2(ref, Wx):.2e}")
        elif part == "W":
            dump(part_W(mesh, W0), "fourd_ledger_W.json")
        elif part == "S":
            for tag, area, nbar, T in SCAN:
                out = run_S(mesh, W0, area, nbar, T=T)
                dump(out, f"fourd_ledger_S_{tag}.json")
                report_S(out)
        elif part == "X":
            wide = RelMesh(nr=256, rh=20.0)
            out = run_S(wide, wide.wigner(psi0), np.inf, 0.0, T=7.0)
            dump(out, "fourd_ledger_X.json")
            print("Part X. The r-window doubled (256 x 64, r in [-20, 20))")
            report_S(out)
        elif part == "B":
            for tag, area in (("inf", np.inf), ("1", H)):
                out = run_B(mesh, W0, W1, area)
                dump(out, f"fourd_ledger_B_{tag}.json")
                report_B(out)
            print("\n  five seeds at t = 1.5: conditional <V> minus reference")
            for tag, area in (("inf", np.inf), ("h/4", 0.25 * H)):
                res = np.array([(r["pot0"] - r["pot0_ref"],
                                 r["pot1"] - r["pot1_ref"], r["N"]) for r in
                                [run_B(mesh, W0, W1, area, T=1.5,
                                       seed=sd)["rec"][-1] for sd in range(5)]])
                se = res.std(axis=0, ddof=1) / np.sqrt(5)
                print(f"  A = {tag:4s}: ground side {res[:, 0].mean():+.3f} +- {se[0]:.3f}"
                      f"   excited side {res[:, 1].mean():+.3f} +- {se[1]:.3f}"
                      f"   N {res[:, 2].mean():.0f}")
        elif part == "F":
            figure()
        print(f"  [part {part}: {time.time() - t0:.0f} s]\n")


def report_S(o):
    a = "inf" if np.isinf(o["area"]) else f"{o['area'] / H:.4g} h"
    print(f"\nS: A = {a}, nbar = {o['nbar']}, run {o['secs']:.0f} s, "
          f"events {o['n_ev']}, S7 residual {o['s7']}"
          + ("  [CAPPED]" if o["capped"] else ""))
    print(f"  {'t':>4} {'N':>8} {'cells':>7} {'N/cell':>6} {'f(win)':>6} "
          f"{'lamA':>5} {'lamB':>5} {'avA':>5} {'avB':>5} {'both':>5} "
          f"{'err':>6} {'errsec':>6} {'<V>':>6} {'<V>ref':>6} {'quiet':>5}")
    prev = None
    for r in o["rec"]:
        if prev is None:
            fw = float("nan")
        else:
            da = r["n_abs"] - prev["n_abs"]
            de = r["n_emi"] - prev["n_emi"]
            fw = da / max(da + de, 1)
        print(f"  {r['t']:4.2f} {r['N']:8d} {r['cells']:7d} "
              f"{r['N'] / r['cells']:6.2f} {fw:6.3f} {r['lamA']:5.2f} "
              f"{r['lamB']:5.2f} {r['availA']:5.2f} {r['availB']:5.2f} "
              f"{r['both']:5.2f} {r['err']:6.3f} {r['err_sec']:6.3f} "
              f"{r['pot']:6.3f} {r['pot_ref']:6.3f} {r['quiet']:5.2f}")
        prev = r


def report_B(o):
    a = "inf" if np.isinf(o["area"]) else f"{o['area'] / H:.4g} h"
    print(f"\nB: A = {a}, d = {o['d']}")
    print(f"  {'t':>4} {'N':>8} {'<V>0':>7} {'ref':>7} {'<V>1':>7} {'ref':>7}"
          f" {'err0':>6} {'err1':>6}")
    for r in o["rec"]:
        print(f"  {r['t']:4.2f} {r['N']:8d} {r['pot0']:7.3f} {r['pot0_ref']:7.3f}"
              f" {r['pot1']:7.3f} {r['pot1_ref']:7.3f} {r['err0']:6.3f}"
              f" {r['err1']:6.3f}")


if __name__ == "__main__":
    main()
