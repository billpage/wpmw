#!/usr/bin/env python3
r"""
Vectorized signed-particle algorithm of D. Cyganski's notebook
``WignerParticlesSinSpawnV1`` (17 April 2020), with every bug found in
``docs/supplement/sinspawn_v1_review.md`` fixed by default and each one
available as a switch, plus an exact split-operator reference.

Physics (SI units, as in the notebook)
--------------------------------------
Electron (m, hbar), potential energy
    U(x) = q_e V2 x^2 + Vp sin(q x),   q = 2 pi nharm / LX,
so the quadratic part is a harmonic oscillator, omega = sqrt(2 q_e V2 / m).
Wigner equation in the wavenumber variable k = p/hbar:
    dW/dt = -(hbar k/m) dW/dx + (1/hbar) U_2'(x) dW/dk
            + Gamma(x) [ W(x, k + q/2) - W(x, k - q/2) ],
    Gamma(x) = (Vp/hbar) cos(q x)                       (Lemma T1).

The algorithm (the notebook's, made correct)
--------------------------------------------
Signed particles (x, k, s), s = +1 positon / -1 negaton, on boxes
(dx, dk) = (LX/NX, 2 pi/LX).  Each time step dt:
  1. spawning: each particle spawns a pair with probability |Gamma(x)| dt
     (the notebook's fixed-dt "roll the dice" clock); the pair sits at the
     parent's x, with sign  s*sgn(Gamma) at k - q/2  and  -s*sgn(Gamma) at
     k + q/2;  the parent persists;
  2. propagation: exact rotation in phase space for the quadratic potential
     (the notebook's own closed form, evaluated with the OLD x and k);
  3. cancellation: opposite signs in the same (x, k) box cancel in pairs
     (the notebook's "kink" at birth and its annihilation on box change,
     done together once per step).

Bug switches (all False = fixed)
--------------------------------
  rectify      : spawn only where cos(qx) > 0, never flip sign   (bug A)
  wrong_array  : negative parents read the positive array        (bug B)
  flip_children: children placed opposite to the QLE            (bug C)
  view_map     : k-update uses the NEW x (numpy view)            (bug D1)
  no_cancel    : no cancellation at all (dead box-change test)   (bug D2)
  cap          : population cap; spawns beyond it are dropped and counted
                 (the notebook's capacity guard, bug J, which in the
                 original code crashes instead)
"""

from __future__ import annotations

import numpy as np

HBAR = 1.05457266e-34
M_E = 9.1093897e-31
Q_E = 1.60217662e-19


class Setup:
    def __init__(self, nharm=2.0, Vp=0.0, V2=None, LX=4e-7, NX=1024, dt=1e-16,
                 x0=None, a=2e-9, k0=0.0):
        self.LX, self.NX, self.dt = LX, NX, dt
        self.dx = LX / NX
        self.dk = 2 * np.pi / LX
        self.NK = NX // 2
        self.V2 = 250.0e6 / LX if V2 is None else V2
        self.s = np.sqrt(2.0 * self.V2 * M_E * Q_E)      # notebook's sqrtvmq
        self.omega = self.s / M_E
        self.nharm = nharm
        self.q = 2 * np.pi * nharm / LX
        self.Vp = Vp
        self.x0 = -0.125 * LX / 2 if x0 is None else x0
        self.a = a
        self.k0 = k0

    def U(self, x):
        return Q_E * self.V2 * x ** 2 + self.Vp * np.sin(self.q * x)

    def Gamma(self, x):
        return (self.Vp / HBAR) * np.cos(self.q * x)


# ----------------------------------------------------------------------------
# exact reference: split-operator Schroedinger equation (|psi|^2 = QLE marginal)
# ----------------------------------------------------------------------------
def schroedinger(S: Setup, times, N=4096, L=None, dt=None):
    L = S.LX if L is None else L
    dt = S.dt if dt is None else dt
    x = (np.arange(N) - N // 2) * (L / N)
    dxs = L / N
    kk = 2 * np.pi * np.fft.fftfreq(N, d=dxs)
    psi = (2 * np.pi * S.a ** 2) ** -0.25 * np.exp(-(x - S.x0) ** 2 / (4 * S.a ** 2)
                                                   + 1j * S.k0 * x)
    half = np.exp(-0.5j * S.U(x) * dt / HBAR)
    kin = np.exp(-0.5j * HBAR * kk ** 2 * dt / M_E)
    out, nsteps = {}, int(round(max(times) / dt))
    want = {int(round(t / dt)): t for t in times}
    for i in range(1, nsteps + 1):
        psi = half * np.fft.ifft(kin * np.fft.fft(half * psi))
        if i in want:
            out[want[i]] = np.abs(psi) ** 2
    return x, dxs, out


def wigner_initial_sample(S: Setup, N0, rng):
    """Positive Gaussian Wigner function, sampled directly (sigma_k = 1/(2a))."""
    x = rng.normal(S.x0, S.a, N0)
    k = rng.normal(S.k0, 1.0 / (2 * S.a), N0)
    return x, k, np.ones(N0, np.int8)


# ----------------------------------------------------------------------------
# the particle algorithm
# ----------------------------------------------------------------------------
class SinSpawn:
    def __init__(self, S: Setup, rng, rectify=False, wrong_array=False,
                 flip_children=False, view_map=False, no_cancel=False, cap=None):
        self.S, self.rng = S, rng
        self.rectify, self.wrong_array = rectify, wrong_array
        self.flip_children, self.view_map, self.no_cancel = flip_children, view_map, no_cancel
        self.cap = cap
        self.kick = S.q / 2.0
        c, sn = np.cos(S.omega * S.dt), np.sin(S.omega * S.dt)
        self.c, self.sn = c, sn
        self.births = 0
        self.cancelled = 0
        self.dropped = 0

    def spawn(self, x, k, sg):
        S = self.S
        g = S.Gamma(x)
        if self.wrong_array:
            # negative parents take the coordinates of positive particle number i
            pos = np.flatnonzero(sg > 0)
            neg = np.flatnonzero(sg < 0)
            xs, ks = x.copy(), k.copy()
            m = min(neg.size, pos.size)
            xs[neg[:m]], ks[neg[:m]] = x[pos[:m]], k[pos[:m]]
            g = S.Gamma(xs)
        else:
            xs, ks = x, k
        if self.rectify:
            prob = np.clip(g, 0, None) * S.dt
            sgam = np.ones_like(g)
        else:
            prob = np.abs(g) * S.dt
            sgam = np.sign(g)
        hit = self.rng.random(x.size) < prob
        if not hit.any():
            return x, k, sg
        xp, kp = xs[hit], ks[hit]
        s_child = (sg[hit] * sgam[hit]).astype(np.int8)       # sign of the child at k - q/2
        if self.flip_children:
            s_child = -s_child
        nb = 2 * xp.size
        if self.cap is not None and x.size + nb > self.cap:
            room = max(0, (self.cap - x.size) // 2)
            self.dropped += xp.size - room
            xp, kp, s_child = xp[:room], kp[:room], s_child[:room]
        self.births += 2 * xp.size
        x = np.concatenate([x, xp, xp])
        k = np.concatenate([k, kp - self.kick, kp + self.kick])
        sg = np.concatenate([sg, s_child, -s_child]).astype(np.int8)
        return x, k, sg

    def propagate(self, x, k):
        S, c, sn = self.S, self.c, self.sn
        xn = k * HBAR * sn / S.s + x * c
        if self.view_map:
            kn = k * c - xn * S.s * sn / HBAR        # notebook: x already overwritten
        else:
            kn = k * c - x * S.s * sn / HBAR         # exact rotation
        return xn, kn

    def cancel(self, x, k, sg):
        if self.no_cancel or x.size == 0:
            return x, k, sg
        S = self.S
        ix = np.rint(x / S.dx).astype(np.int64)
        ik = np.rint(k / S.dk).astype(np.int64)
        cid = (ix + (1 << 20)) * (1 << 21) + (ik + (1 << 20))
        key = self.rng.random(x.size)
        order = np.lexsort((key, sg, cid))
        cs, ss = cid[order], sg[order]
        grp = np.concatenate([[True], (cs[1:] != cs[:-1]) | (ss[1:] != ss[:-1])])
        gstart = np.flatnonzero(grp)
        gsize = np.diff(np.append(gstart, cs.size))
        gid = np.cumsum(grp) - 1
        rank = np.arange(cs.size) - gstart[gid]
        ucell, cinv = np.unique(cs, return_inverse=True)
        npos = np.bincount(cinv, weights=(ss > 0), minlength=ucell.size)
        nneg = np.bincount(cinv, weights=(ss < 0), minlength=ucell.size)
        opp = np.where(ss > 0, nneg[cinv], npos[cinv])
        keep = rank < (gsize[gid] - opp)
        self.cancelled += int((~keep).sum())
        idx = order[keep]
        return x[idx], k[idx], sg[idx]

    def run(self, x, k, sg, times, log_every=100):
        S = self.S
        nsteps = int(round(max(times) / S.dt))
        want = {int(round(t / S.dt)): t for t in times}
        out, log = {}, []
        for i in range(1, nsteps + 1):
            x, k, sg = self.spawn(x, k, sg)
            x, k = self.propagate(x, k)
            x, k, sg = self.cancel(x, k, sg)
            if i % log_every == 0:
                log.append((i * S.dt, x.size, int((sg > 0).sum()), int((sg < 0).sum())))
            if i in want:
                out[want[i]] = (x.copy(), k.copy(), sg.copy())
        return out, np.array(log)


def marginal(S: Setup, x, sg, N0, edges):
    h, _ = np.histogram(x, bins=edges, weights=sg.astype(float))
    return h / N0 / np.diff(edges)
