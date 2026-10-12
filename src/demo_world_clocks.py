#!/usr/bin/env python3
"""
Companion demo for ``docs/analysis/world_clocks.md`` (step 25).

Where does a reading live: in a world, or in one of its particles?  Two
particles on a line, unit masses, hbar = 1, interacting through a pair
potential V(x1 - x2).  The compensated residual is read three ways:

  W   the world reads the joint chord (y1, y2):
        U_W = V(r + y1 - y2) - V(r - y1 + y2)                  (step 21, M1)
  B1  each particle reads its own chord, the other held fixed, its clock
      winding with the whole interaction (the other as an external field):
        U_B1 = [V(r + y1) - V(r - y1)] + [V(r - y2) - V(r + y2)]
  B2  as B1, with the interaction energy split equally between the two
      particle clocks:  U_B2 = U_B1 / 2

with r = x1 - x2, each less its linear (compensated) part.  Coordinates
P = p1 + p2 and p_r = (p1 - p2)/2; the centre of mass is free, so X is
never read and the state is carried as W(P, r, p_r).

Parts (run any subset; default all but C)
-----------------------------------------
A  Exact (SymPy, generic Taylor coefficients): the three symbols, the orders
   at which they differ, B1's relative marginal is the QLE at hbar/2, the
   P-marginal is untouched, the Lagrangian bookkeeping of the clocks, and
   the pair-bank structure for three particles.
K  The kernels: event rates of the world bank and of the particle banks.
H  Exact references: the hbar/2 von Neumann evolution that B1's relative
   motion must follow, for the Poeschl-Teller ground state.
M  Mesh runs of W, B1, B2 in (P, r, p_r) from g(P) x W_psi0(r, p_r), which
   exact QM keeps stationary and separable.  Checks B1's relative marginal
   against an independent 2D mesh at hbar/2 and against H.
R  The relative motion's error against the exact references at reach
   2 pi, 4 pi, 8 pi, and at 2 pi on a grid refined in r and t.
T  The Eckart pair collision of step 21: transmission under W and B1
   (relative marginal), against the split-operator references at hbar and
   hbar/2, at three reaches (about 18 minutes).
C  Particle runs with Proposition Q11's integrator banks: one bank per world
   (W) against one bank per particle (B1), emission plus a contact sink.
   ``--seeds N`` runs a sweep, one CSV row per (variant, seed).
F  The two figures, from the JSON that M, H, T and C save (every Part C
   file present is overlaid; --fig-nu picks the one panel (e) shows).

Run:  WPMW_OUTPUT=... PYTHONPATH=src python3 -u src/demo_world_clocks.py [A K H M R T C F]
Sweep (Kaggle):  ... demo_world_clocks.py C --nu 100000 --seeds 8 --Tc 4 --Tb1 1 --workers 4
"""
import argparse
import json
import os
import time

import numpy as np

from wpmwlib.wpmw_utils import docs_path, output_path

HBAR, MU, MTOT = 1.0, 0.5, 2.0
V0, A_W = 6.0, 1.0                    # Poeschl-Teller well of step 21
DP = 0.25                             # coarse momentum row; reach y_max = pi hbar / (2 DP)
Y_MAX = np.pi * HBAR / (2 * DP)       # = 2 pi a
SIG_P = 1.0                           # width of the free centre-of-mass packet in P
WIG_TIMES = tuple(0.25 * i for i in range(1, 17))   # times of exact Wigner comparisons


def V_pt(r):
    return -V0 / np.cosh(r / A_W) ** 2


def dV_pt(r):
    return 2 * V0 / A_W * np.tanh(r / A_W) / np.cosh(r / A_W) ** 2


def psi0(r):
    return np.sqrt(0.75) / np.cosh(r) ** 2


def rule(text):
    print()
    print(text)
    print("-" * len(text))


def dump(obj, name):
    with open(output_path(name), "w") as f:
        json.dump(obj, f, indent=1)


def load(name):
    p = output_path(name)
    if not os.path.exists(p):
        return None
    with open(p) as f:
        return json.load(f)


# ======================================================================
# Part A -- exact statements
# ======================================================================
def part_A():
    import sympy as sp

    rule("Part A. Exact: the three readings of a pair potential (SymPy)")
    y1, y2, s, sP, sr, hb = sp.symbols("y1 y2 s s_P s_r hbar", real=True)
    NORD = 9
    a = sp.symbols(f"a0:{NORD + 1}", real=True)     # a_k = V^(k)(r)

    def Vt(z):                                       # V(r + z) by Taylor coefficients
        return sum(a[k] * z**k / sp.factorial(k) for k in range(NORD + 1))

    def res(expr, lin):                              # less its linear part
        return sp.expand(expr - lin)

    UW = res(Vt(y1 - y2) - Vt(-(y1 - y2)), 2 * a[1] * (y1 - y2))
    U1 = res(Vt(y1) - Vt(-y1), 2 * a[1] * y1) + res(Vt(-y2) - Vt(y2), -2 * a[1] * y2)

    def degree_part(expr, k):
        return sp.Add(*[t for t in sp.Add.make_args(expr)
                        if sp.Poly(t, y1, y2).total_degree() == k])

    print("  A1  U_W - U_B1 by order (U_B2 = U_B1/2):")
    out = {}
    for k in (1, 3, 5, 7):
        d = sp.factor(degree_part(UW - U1, k))
        print(f"      order {k}: {d}")
        out[f"diff_{k}"] = str(d)
    lin_B2 = sp.expand(degree_part(Vt(y1) - Vt(-y1) + Vt(-y2) - Vt(y2), 1) / 2)
    print(f"      B2's unreduced linear term: {sp.factor(lin_B2)}  "
          f"(the world's: {sp.factor(2 * a[1] * (y1 - y2))})")

    # A2: relative marginal of B1 is the world symbol at hbar/2, to all orders
    Vf = sp.Function("V")
    r = sp.symbols("r", real=True)

    def M(hbar_, sig, pot=Vf):
        y = hbar_ * sig / 2
        return sp.I / hbar_ * (pot(r + y) - pot(r - y))

    sig1, sig2 = sP + sr / 2, sP - sr / 2
    B1 = M(hb, sig1) - M(hb, sig2)
    W_half = M(hb / 2, sr)
    a2 = sp.simplify(B1.subs(sP, 0) - W_half)
    print(f"  A2  B1 symbol at s_P = 0  minus  world symbol at hbar/2:  {a2}")
    a3 = sp.simplify(B1.subs(sr, 0))
    print(f"  A3  B1 symbol at s_r = 0 (acts on the P-marginal):        {a3}")
    a3w = sp.simplify(sp.diff(M(hb, sr), sP))
    print(f"      world symbol's dependence on s_P:                     {a3w}")
    # cubic coefficient of the relative-marginal symbol, B1 and B2 against W
    def Mt(hbar_, sig):
        y = hbar_ * sig / 2
        return sp.I / hbar_ * (Vt(y) - Vt(-y))

    rel = {"W": Mt(hb, sr), "B1": Mt(hb, sr / 2) - Mt(hb, -sr / 2)}
    rel["B2"] = rel["B1"] / 2
    cub = {k: sp.factor(sp.expand(v).coeff(sr, 3)) for k, v in rel.items()}
    print(f"  A4  cubic coefficient of the relative-marginal symbol: "
          f"W {cub['W']},  B1 {cub['B1']},  B2 {cub['B2']}")
    print(f"      ratios B1/W = {sp.simplify(cub['B1'] / cub['W'])}, "
          f"B2/W = {sp.simplify(cub['B2'] / cub['W'])}")
    # in (P, p_r) derivatives: (d1 - d2)^3 against d1^3 - d2^3
    dP, dr_ = sp.symbols("d_P d_pr")
    d1, d2 = dP + dr_ / 2, dP - dr_ / 2
    print(f"      Moyal cubic operators: world (d1 - d2)^3 = {sp.expand((d1 - d2) ** 3)};  "
          f"B1 d1^3 - d2^3 = {sp.expand(d1**3 - d2**3)}")

    # A5: the clocks' Lagrangian bookkeeping
    p1, p2, x1, x2, al = sp.symbols("p1 p2 x1 x2 alpha", real=True)
    Vx = Vf(x1 - x2)
    L = p1**2 / 2 + p2**2 / 2 - Vx
    LB1 = (p1**2 / 2 - Vx) + (p2**2 / 2 - Vx)
    Lal = (p1**2 / 2 - al * Vx) + (p2**2 / 2 - (1 - al) * Vx)
    print(f"  A5  clock rates: sum of B1 particle Lagrangians minus the world's: "
          f"{sp.simplify(LB1 - L)};  split by alpha: {sp.simplify(Lal - L)}")
    # A6: gauge, V -> V + c: chord misalignment rate is unchanged in all three
    c = sp.symbols("c", real=True)
    gauge = sp.simplify(M(hb, sr, pot=lambda z: Vf(z) + c) - M(hb, sr))
    print(f"  A6  V -> V + c changes the world symbol by {gauge}")

    # A7: three particles, pairwise potential: the world symbol is a sum of
    # pair symbols, each a function of one difference sigma_i - sigma_j.
    s1_, s2_, s3_ = sp.symbols("sigma1 sigma2 sigma3", real=True)
    x = sp.symbols("x1:4", real=True)
    sig = (s1_, s2_, s3_)
    pairs = [(0, 1), (0, 2), (1, 2)]
    Vtot = lambda xs: sum(Vf(xs[i] - xs[j]) for i, j in pairs)
    yv = [hb * q / 2 for q in sig]
    Mw = sp.I / hb * (Vtot([x[k] + yv[k] for k in range(3)])
                      - Vtot([x[k] - yv[k] for k in range(3)]))
    Mpairs = sum(sp.I / hb * (Vf(x[i] - x[j] + hb * (sig[i] - sig[j]) / 2)
                              - Vf(x[i] - x[j] - hb * (sig[i] - sig[j]) / 2))
                 for i, j in pairs)
    print(f"  A7  N = 3: world symbol minus the sum of three pair-exchange symbols: "
          f"{sp.simplify(Mw - Mpairs)}")
    return dict(out, A2=str(a2), A3=str(a3))


# ======================================================================
# Symbols on meshes
# ======================================================================
def m_res(r, sig, hbar=HBAR, pot=V_pt, dpot=dV_pt, ymax=Y_MAX):
    """Compensated residual symbol (numpy fft convention), raised-cosine
    horizon at the reach y_max (a length, the same for every hbar):
        (i/hbar) [V(r+y) - V(r-y) - 2 y V'(r)] cos^2(pi y / 2 y_max),
    y = hbar sigma / 2, extended periodically in sigma with period
    4 y_max / hbar, so that the kernel lives on rows of pi hbar / 2 y_max
    (DP at hbar = 1, Theorem C4 of step 14)."""
    per = 4 * ymax / hbar
    sig = (sig + per / 2) % per - per / 2
    y = hbar * sig / 2
    w = np.where(np.abs(y) < ymax, np.cos(np.pi * y / (2 * ymax)) ** 2, 0.0)
    return 1j / hbar * (pot(r + y) - pot(r - y) - 2 * y * dpot(r)) * w


def kernel_rates(r):
    """K_q(r) for q = 1 .. 31 on the coarse row lattice (64 rows of DP), as in
    step 21: K = Re ifft(symbol).  Gross rate Gamma(r) = sum_q |K_q|."""
    n = 64
    s = 2 * np.pi * np.fft.fftfreq(n, d=DP)
    sym = m_res(r[:, None], s[None, :])
    k = np.real(np.fft.ifft(sym, axis=1))
    return k[:, 1:n // 2]


class Mesh3:
    """W(P, r, p_r) with P on nP rows of DP, r on nr points, p_r on npr rows
    of DP/2 (so that a particle-bank event, which moves p_r by xi/2, lands on
    a row).  Strang splitting: half drift in r, then classical force and
    residual together in (s_P, s_r), then half drift."""

    def __init__(self, variant, nP=64, nr=128, rh=10.0, npr=128, dt=0.01,
                 pot=V_pt, dpot=dV_pt, hbar=HBAR):
        self.variant, self.dt = variant, dt
        self.P = DP * (np.arange(nP) - nP // 2)
        self.r = -rh + 2 * rh * np.arange(nr) / nr
        self.dr = self.r[1] - self.r[0]
        self.dpr = DP / 2
        self.pr = self.dpr * (np.arange(npr) - npr // 2)
        self.kr = 2 * np.pi * np.fft.fftfreq(nr, d=self.dr)
        sP = 2 * np.pi * np.fft.fftfreq(nP, d=DP)
        sr = 2 * np.pi * np.fft.fftfreq(npr, d=self.dpr)
        R = self.r[None, :, None]
        SP, SR = sP[:, None, None], sr[None, None, :]
        cls = 1j * dpot(R) * SR                        # classical force, exact
        if variant == "W":
            res = m_res(R, SR + 0 * SP, hbar, pot, dpot)
        elif variant in ("B1", "B2"):
            res = m_res(R, SP + SR / 2, hbar, pot, dpot) - m_res(R, SP - SR / 2, hbar, pot, dpot)
            if variant == "B2":
                res = 0.5 * res
        elif variant == "cl":
            res = 0 * cls
        else:
            raise ValueError(variant)
        self.kick = np.exp(dt * (cls + res))
        self.drift = np.exp(-1j * self.kr[:, None] * self.pr[None, :] / MU * dt / 2)[None]

    def step(self, W):
        W = np.real(np.fft.ifft(np.fft.fft(W, axis=1) * self.drift, axis=1))
        W = np.real(np.fft.ifft2(np.fft.fft2(W, axes=(0, 2)) * self.kick, axes=(0, 2)))
        return np.real(np.fft.ifft(np.fft.fft(W, axis=1) * self.drift, axis=1))


class Mesh2:
    """W(r, p_r) on the same r and p_r mesh; residual symbol at hbar_res
    (the classical part does not depend on hbar)."""

    def __init__(self, hbar_res=HBAR, nr=128, rh=10.0, npr=128, dpr=None,
                 dt=0.01, pot=V_pt, dpot=dV_pt, residual=True, ymax=Y_MAX):
        if dpr is None:                    # half the hbar row pi hbar / 2 y_max
            dpr = np.pi * HBAR / (4 * ymax)
        self.r = -rh + 2 * rh * np.arange(nr) / nr
        self.dr = self.r[1] - self.r[0]
        self.dpr = dpr
        self.pr = dpr * (np.arange(npr) - npr // 2)
        self.kr = 2 * np.pi * np.fft.fftfreq(nr, d=self.dr)
        sr = 2 * np.pi * np.fft.fftfreq(npr, d=dpr)
        R, SR = self.r[:, None], sr[None, :]
        res = m_res(R, SR, hbar_res, pot, dpot, ymax) if residual else 0
        self.kick = np.exp(dt * (1j * dpot(R) * SR + res))
        self.drift = np.exp(-1j * self.kr[:, None] * self.pr[None, :] / MU * dt / 2)

    def step(self, W):
        W = np.real(np.fft.ifft(np.fft.fft(W, axis=0) * self.drift, axis=0))
        W = np.real(np.fft.ifft(np.fft.fft(W, axis=1) * self.kick, axis=1))
        return np.real(np.fft.ifft(np.fft.fft(W, axis=0) * self.drift, axis=0))


def wigner_mesh(psi_fn, r, pr, hbar=HBAR):
    """W(r, p) of a real or complex wavefunction, by quadrature in y."""
    yq = np.arange(-16, 16, 0.004)
    W = np.empty((r.size, pr.size))
    for i, r0 in enumerate(r):
        prod = psi_fn(r0 + yq) * np.conj(psi_fn(r0 - yq))
        W[i] = np.real(np.exp(-2j * np.outer(pr, yq) / hbar) @ prod) * 0.004 / (np.pi * hbar)
    return W


# ======================================================================
# Part K -- event rates
# ======================================================================
def part_K():
    rule("Part K. Event rates: one bank per world against one per particle")
    r = np.linspace(-4, 4, 801)
    k = kernel_rates(r)
    gam = np.abs(k).sum(axis=1)
    rho0 = psi0(r) ** 2
    dr = r[1] - r[0]
    mean_g = float((gam * rho0).sum() * dr)
    print(f"  channels per bank: {k.shape[1]} (q = 1 .. {k.shape[1]}), xi_q = q DP, "
          f"DP = {DP}, y_max = {Y_MAX:.4f}")
    print(f"  world bank, one per world:      Gamma(r) max {gam.max():.4f}, "
          f"<Gamma>_psi0 = {mean_g:.4f}; each event moves (p1, p2) by (+-xi, -+xi)")
    print(f"  particle banks (B1), two per world: gross rate 2 Gamma(r), <.>_psi0 = "
          f"{2 * mean_g:.4f}; each event moves one p_i by +-xi, so P by +-xi, p_r by +-xi/2")
    print(f"  integrators carried: W {k.shape[1]} per world, B1 {2 * k.shape[1]} per world")
    return dict(gamma_max=float(gam.max()), gamma_mean=mean_g, nq=int(k.shape[1]))


# ======================================================================
# Part H -- the exact hbar/2 reference for B1's relative motion
# ======================================================================
def vn_half(rho_fn, obs_times, n=1024, L=40.0, hbar=HBAR / 2, pot=V_pt, wigner_at=None):
    """rho'(a, b) = rho_fn(a, b) evolved by von Neumann at hbar' with
    H' = p^2/2mu + V(r).  Returns <V>(t) and <p^2/2mu>(t) at obs_times, the
    trace, and the spectrum's extremes."""
    x = (np.arange(n) - n // 2) * (L / n)
    dx = x[1] - x[0]
    k = 2 * np.pi * np.fft.fftfreq(n, d=dx)
    F = np.fft.fft(np.eye(n), axis=0) / np.sqrt(n)
    T = (F.conj().T * ((hbar * k) ** 2 / (2 * MU))) @ F
    H = np.real(T) + np.diag(pot(x))
    ev, U = np.linalg.eigh(H)
    A, B = np.meshgrid(x, x, indexing="ij")
    rho = rho_fn(A, B) * dx
    lam = np.linalg.eigvalsh(rho)
    rt = U.T @ rho @ U
    Vd, Tm = U.T @ np.diag(pot(x)) @ U, U.T @ T.real @ U
    outV, outT, wig = [], [], {}
    for t in obs_times:
        ph = np.exp(-1j * ev * t / hbar)
        r_t = (ph[:, None] * rt) * ph.conj()[None, :]
        outV.append(float(np.real(np.sum(r_t * Vd.T))))
        outT.append(float(np.real(np.sum(r_t * Tm.T))))
        if wigner_at is not None and round(t, 6) in wigner_at[2]:
            rx = U @ r_t @ U.T / dx                    # rho'(a, b) as a kernel
            rm, pm = wigner_at[0], wigner_at[1]
            ib = np.round((rm - x[0]) / dx).astype(int)
            j = np.arange(-n // 4, n // 4)
            Wt = np.zeros((rm.size, pm.size))
            for i, b in enumerate(ib):
                ok = (b + j >= 0) & (b + j < n) & (b - j >= 0) & (b - j < n)
                vals = rx[b + j[ok], b - j[ok]]
                Wt[i] = np.real(np.exp(-2j * np.outer(pm, j[ok] * dx) / hbar) @ vals) * dx / (np.pi * hbar)
            wig[round(t, 6)] = Wt
    return dict(V=outV, T=outT, trace=float(np.trace(rho)), lam_min=float(lam.min()),
                lam_max=float(lam.max()), neg_mass=float(-lam[lam < 0].sum()), wig=wig)


def part_H(obs_times):
    rule("Part H. Exact references for the bound state psi0 ~ sech^2 r")
    # exact QM: stationary, <V> = -4.8, <T> = 0.8, E0 = -4
    rq = np.linspace(-30, 30, 600001)
    dr = rq[1] - rq[0]
    Vq = float((V_pt(rq) * psi0(rq) ** 2).sum() * dr)
    print(f"  exact QM at hbar: psi0 stationary, <V> = {Vq:.6f} (closed form -4.8), "
          f"<T> = {-4 - Vq:.6f}, E0 = -4")
    # B1's relative marginal: W_psi0 read as an hbar/2 Wigner function
    rho_fn = lambda a, b: psi0((3 * a - b) / 2) * psi0((3 * b - a) / 2)
    m0 = Mesh2()
    h = vn_half(rho_fn, obs_times, wigner_at=(m0.r, m0.pr, WIG_TIMES))
    print("  B1 reference (von Neumann at hbar/2 of rho'(a,b) = psi0((3a-b)/2) psi0((3b-a)/2)):")
    print(f"    trace {h['trace']:.6f}; spectrum of rho' in [{h['lam_min']:.4e}, {h['lam_max']:.4f}],"
          f" negative mass {h['neg_mass']:.4e}")
    for t, v, k in zip(obs_times, h["V"], h["T"]):
        if t in (0.0, 0.5, 1.0, 2.0, 4.0):
            print(f"    t = {t:4.1f}:  <V> = {v:+.5f}   <T> = {k:.5f}   E_rel = {v + k:+.5f}")
    return dict(V_exact=Vq, half=h, times=list(obs_times))


# ======================================================================
# Part M -- mesh runs of the three readings
# ======================================================================
def part_M(T=4.0, dt=0.01, every=10):
    rule(f"Part M. Mesh runs in (P, r, p_r), t = 0 .. {T}, dt = {dt}")
    variants = ("W", "B1", "B2", "cl")
    m0 = Mesh3("W", dt=dt)
    W0r = wigner_mesh(psi0, m0.r, m0.pr)
    W0r /= W0r.sum() * m0.dr * m0.dpr
    g = np.exp(-m0.P**2 / (2 * SIG_P**2))
    g /= g.sum() * DP
    W0 = g[:, None, None] * W0r[None]
    vol = DP * m0.dr * m0.dpr
    Vr = V_pt(m0.r)[None, :, None]
    Tr = (m0.pr**2 / (2 * MU))[None, None, :]
    P2 = (m0.P**2)[:, None, None]
    print(f"  mesh: P {m0.P.size} rows of {DP}, r {m0.r.size} points dr = {m0.dr:.4f}, "
          f"p_r {m0.pr.size} rows of {m0.dpr}; reach y_max = {Y_MAX:.4f}; sigma_P = {SIG_P}")

    def obs(W):
        Wrel = W.sum(axis=0) * DP
        gP = W.sum(axis=(1, 2)) * m0.dr * m0.dpr
        n = W.sum() * vol
        prod = gP[:, None, None] * Wrel[None] / n
        return dict(norm=float(n),
                    V=float((W * Vr).sum() * vol),
                    T=float((W * Tr).sum() * vol),
                    sep=float(np.linalg.norm(W - prod) / np.linalg.norm(W)),
                    gP=float(np.abs(gP - g).sum() * DP),
                    corr=float((W * P2 * Vr).sum() * vol - (W * P2).sum() * vol * (W * Vr).sum() * vol),
                    purity=float(2 * np.pi * HBAR * (Wrel**2).sum() * m0.dr * m0.dpr),
                    Wrel=Wrel)

    res, finals, snaps = {}, {}, {}
    t_obs = [0.0]
    for v in variants:
        mesh = Mesh3(v, dt=dt)
        W = W0.copy()
        o = obs(W)
        rows = [{k: o[k] for k in o if k != "Wrel"}]
        rel0 = o["Wrel"]
        tt = time.perf_counter()
        nst = int(round(T / dt))
        for i in range(1, nst + 1):
            W = mesh.step(W)
            if i % every == 0:
                o = obs(W)
                row = {k: o[k] for k in o if k != "Wrel"}
                row["rel_drift"] = float(np.linalg.norm(o["Wrel"] - rel0) / np.linalg.norm(rel0))
                rows.append(row)
                if v == variants[0]:
                    t_obs.append(round(i * dt, 6))
                if round(i * dt, 6) in WIG_TIMES:
                    snaps.setdefault(v, {})[round(i * dt, 6)] = o["Wrel"]
        rows[0]["rel_drift"] = 0.0
        res[v] = rows
        finals[v] = o["Wrel"]
        f = rows[-1]
        print(f"  {v:3s} ({time.perf_counter() - tt:5.1f} s)  t = {T}:  <V> {f['V']:+.5f} "
              f"(t=0 {rows[0]['V']:+.5f})  E_rel {f['V'] + f['T']:+.5f}  "
              f"rel. drift {f['rel_drift']:.3e}  separability {f['sep']:.3e}  "
              f"|dg_P|_1 {f['gP']:.1e}  cov(P^2,V) {f['corr']:+.2e}  norm {f['norm']:.12f}")
    # B1's relative marginal against an independent 2D mesh at hbar/2
    m2 = Mesh2(hbar_res=HBAR / 2, dt=dt)
    w = W0r.copy()
    for i in range(int(round(T / dt))):
        w = m2.step(w)
    d_half = float(np.abs(finals["B1"] - w).max() / np.abs(w).max())
    m2w = Mesh2(hbar_res=HBAR, dt=dt)
    w1 = W0r.copy()
    for i in range(int(round(T / dt))):
        w1 = m2w.step(w1)
    d_full = float(np.abs(finals["W"] - w1).max() / np.abs(w1).max())
    print(f"  B1's relative marginal against a 2D mesh at hbar/2: max rel. difference {d_half:.2e}")
    print(f"  W's  relative marginal against a 2D mesh at hbar:   max rel. difference {d_full:.2e}")
    return dict(times=t_obs, runs=res, d_half=d_half, d_full=d_full, W0r=W0r, snaps=snaps)


def report_MH(m, h):
    rule("Part M against Part H: <V>(t)")
    times = m["times"]
    hv = dict(zip([round(t, 6) for t in h["times"]], h["half"]["V"]))
    print("     t    exact QM     W mesh     B1 mesh   B1 exact(hbar/2)   B2 mesh   classical")
    for i, t in enumerate(times):
        if abs(t - round(t)) < 1e-9 or t in (0.25, 0.5):
            print(f"  {t:4.2f}   {h['V_exact']:+.5f}   {m['runs']['W'][i]['V']:+.5f}   "
                  f"{m['runs']['B1'][i]['V']:+.5f}   {hv.get(round(t, 6), float('nan')):+.5f}"
                  f"           {m['runs']['B2'][i]['V']:+.5f}   {m['runs']['cl'][i]['V']:+.5f}")
    errB = max(abs(m["runs"]["B1"][i]["V"] - hv[round(t, 6)]) for i, t in enumerate(times))
    errW = max(abs(m["runs"]["W"][i]["V"] - h["V_exact"]) for i, t in enumerate(times))
    sepB = max(rw["sep"] for rw in m["runs"]["B1"])
    sepW = max(rw["sep"] for rw in m["runs"]["W"])
    print(f"  max |W mesh - exact QM| = {errW:.2e};  max |B1 mesh - exact hbar/2| = {errB:.2e}")
    print("  relative Wigner function, L2 distance relative to |W0|:")
    print("     t    W mesh vs exact QM   B1 mesh vs exact hbar/2   B1 mesh vs exact QM   exact hbar/2 vs exact QM")
    W0r = m["W0r"]
    nrm = np.linalg.norm(W0r)
    wl2 = {}
    for t in WIG_TIMES:
        if t not in m["snaps"]["W"]:
            continue
        ex = h["half"]["wig"][t]
        a_ = np.linalg.norm(m["snaps"]["W"][t] - W0r) / nrm
        b_ = np.linalg.norm(m["snaps"]["B1"][t] - ex) / nrm
        c_ = np.linalg.norm(m["snaps"]["B1"][t] - W0r) / nrm
        d_ = np.linalg.norm(ex - W0r) / nrm
        wl2[t] = (a_, b_, c_, d_)
        if abs(t - round(t)) < 1e-9:
            print(f"  {t:4.1f}        {a_:.3e}             {b_:.3e}                {c_:.3e}             {d_:.3e}")
    print(f"  max separability defect: W {sepW:.2e}, B1 {sepB:.2e}, "
          f"B2 {max(rw['sep'] for rw in m['runs']['B2']):.2e}")
    return dict(errW=errW, errB=errB, sepW=sepW, sepB=sepB,
                wig_l2={str(t): v for t, v in wl2.items()})


# ======================================================================
# Part R -- what the mesh's own error is: the reach
# ======================================================================
def part_R(T=4.0, dt=0.01, reaches=(2, 4, 8)):
    rule(f"Part R. Relative motion of psi0 at t = {T} against the exact references, by reach")
    rho_fn = lambda a, b: psi0((3 * a - b) / 2) * psi0((3 * b - a) / 2)
    print("     reach y_max   grid           W vs exact QM   B1 (= hbar/2) vs exact hbar/2   B1 vs exact QM")
    out = {}
    cases = [(mult, 128, dt) for mult in reaches] + [(reaches[0], 256, dt / 2)]
    for mult, nr, dt_ in cases:
        ym = mult * np.pi
        dpr = np.pi * HBAR / (4 * ym)
        npr = int(round(16 / dpr))
        row = []
        for hb in (HBAR, HBAR / 2):
            m = Mesh2(hbar_res=hb, nr=nr, rh=10.0, npr=npr, dpr=dpr, dt=dt_, ymax=ym)
            W0 = wigner_mesh(psi0, m.r, m.pr)
            W0 /= W0.sum() * m.dr * m.dpr
            W = W0.copy()
            for _ in range(int(round(T / dt_))):
                W = m.step(W)
            nrm = np.linalg.norm(W0)
            if hb == HBAR:
                row.append(float(np.linalg.norm(W - W0) / nrm))
            else:
                ex = vn_half(rho_fn, [T], wigner_at=(m.r, m.pr, (T,)))["wig"][T]
                row += [float(np.linalg.norm(W - ex) / nrm), float(np.linalg.norm(W - W0) / nrm)]
        out[f"{mult}_{nr}_{dt_}"] = row
        print(f"     {mult:2d} pi    {nr}x{npr}, dt {dt_:<6}  {row[0]:.3e}          {row[1]:.3e}                 {row[2]:.3e}")
    return out


# ======================================================================
# Part T -- the Eckart pair collision under W and B1
# ======================================================================
def part_T(reaches=(2, 4, 8)):
    rule("Part T. Eckart pair collision (step 21 Part E): transmission")
    v0, a, sig_r, r_c, p_c, T, dt = 1.0, 2.0, 2.0, -16.0, 1.0, 22.0, 0.01
    U = lambda r: v0 / np.cosh(r / a) ** 2
    dU = lambda r: -2 * v0 / a * np.tanh(r / a) / np.cosh(r / a) ** 2
    chi = lambda r: (2 * np.pi * sig_r**2) ** -0.25 * np.exp(-(r - r_c) ** 2 / (4 * sig_r**2)
                                                             + 1j * p_c * r / HBAR)
    # split-operator references for the relative motion at hbar and hbar/2;
    # at hbar/2 the initial state is rho'(a, b) = chi((3a-b)/2) chi*((3b-a)/2)
    n, L = 2048, 160.0
    x = (np.arange(n) - n // 2) * (L / n)
    dx = x[1] - x[0]
    k = 2 * np.pi * np.fft.fftfreq(n, d=dx)

    def propagate(vecs, hb):
        e1 = np.exp(-0.5j * U(x) * dt / hb)
        t1 = np.exp(-1j * (hb * k) ** 2 / (2 * MU) * dt / hb)
        for _ in range(int(round(T / dt))):
            vecs = e1[:, None] * np.fft.ifft(t1[:, None] * np.fft.fft(e1[:, None] * vecs, axis=0), axis=0)
        return vecs

    def trans(dens):
        return float(dens[x > 0].sum() + 0.5 * dens[x == 0].sum())

    v = chi(x)[:, None] * np.sqrt(dx)
    T_q = trans(np.abs(propagate(v, HBAR)[:, 0]) ** 2)
    A_, B_ = np.meshgrid(x, x, indexing="ij")
    rho = chi((3 * A_ - B_) / 2) * np.conj(chi((3 * B_ - A_) / 2)) * dx
    lam, vec = np.linalg.eigh(rho)
    keep = np.abs(lam) > 1e-12
    vt = propagate(vec[:, keep], HBAR / 2)
    dens_h = (np.abs(vt) ** 2 * lam[keep][None, :]).sum(axis=1)
    T_h = trans(dens_h)
    print(f"  split-operator, relative motion: T(hbar) = {T_q:.6f}  (step 21: 0.544134)")
    print(f"  split-operator at hbar/2 on rho': T = {T_h:.6f}   (rho' trace {lam.sum():.6f}, "
          f"{int(keep.sum())} eigenvectors, min eigenvalue {lam.min():.2e})")
    # mesh runs of the relative marginal: W at hbar, B1 = hbar/2, classical
    out = dict(T_q=T_q, T_h=T_h, mesh={})
    print("  mesh, relative marginal (r in [-80, 80), 1024 points; p_r in [-8, 8)):")
    print("     reach y_max     W (hbar)    B1 (= hbar/2)   classical")
    for mult in reaches:
        ym = mult * np.pi
        row = {}
        for name, hb, resid in (("W", HBAR, True), ("B1", HBAR / 2, True), ("cl", HBAR, False)):
            if name == "cl" and mult != reaches[0]:
                row[name] = out["mesh"][str(reaches[0])]["cl"]
                continue
            dpr = np.pi * HBAR / (4 * ym)
            m = Mesh2(hbar_res=hb, nr=1024, rh=80.0, npr=int(round(16 / dpr)), dpr=dpr,
                      dt=dt, pot=U, dpot=dU, residual=resid, ymax=ym)
            W = wigner_mesh(chi, m.r, m.pr)
            W /= W.sum() * m.dr * m.dpr
            for _ in range(int(round(T / dt))):
                W = m.step(W)
            rho_r = W.sum(axis=1) * m.dpr
            row[name] = float(rho_r[m.r > 0].sum() * m.dr)
        out["mesh"][str(mult)] = row
        print(f"     {mult:2d} pi           {row['W']:.6f}    {row['B1']:.6f}        {row['cl']:.6f}")
    last = out["mesh"][str(reaches[-1])]
    print(f"  at the largest reach: |W - T(hbar)| = {abs(last['W'] - T_q):.2e}, "
          f"|B1 - T(hbar/2)| = {abs(last['B1'] - T_h):.2e}; B1 - W = {last['B1'] - last['W']:+.4f} "
          f"(exact: {T_h - T_q:+.4f})")
    return out


# ======================================================================
# Part C -- particle runs with Proposition Q11's integrator banks
# ======================================================================
C_EVERY = 0.25                        # observation interval of Part C
N_CAP = 4_000_000                     # stop a particle run whose population passes this


class Banked:
    """Signed worlds (P, r, p_r, eps), X never read.  Postulate (S): velocity
    Verlet in r under the exact force, P free.  Postulate (D) realised by
    integrate-and-fire banks (Proposition Q11): every integrator advances
    by its channel's rate along the world's worldline and fires each time it
    crosses an integer, with orientation eps * sign of the crossing.  A firing
    is an emission: daughters of sign +tau and -tau at the two shifted
    momenta, with fresh uniform integrator phases; the parent continues.
    A contact sink (garbage collection) then removes opposite-sign worlds in
    the same mesh cell, which leaves the cell's net count unchanged.

    variant 'W' : one bank per world, C_q' = K_q(r); a firing moves
                  (p1, p2) by (xi, -xi): p_r by xi, P unchanged.
    variant 'B1': one bank per particle, C1_q' = K_q(r), C2_q' = -K_q(r); a
                  firing of particle i moves p_i by xi: P by xi, p_r by
                  +xi/2 (i = 1) or -xi/2 (i = 2)."""

    def __init__(self, variant, nu, seed, mesh):
        self.v, self.nu, self.rng, self.m = variant, nu, np.random.default_rng(seed), mesh
        self.rt = np.linspace(-mesh["rh"], mesh["rh"], 8001)
        self.kt = kernel_rates(self.rt)                         # (n_r table, NQ)
        self.nq = self.kt.shape[1]
        self.xi = DP * np.arange(1, self.nq + 1)
        self.nb = 1 if variant == "W" else 2
        self.stats = dict(fires=0, sunk=0, dP_viol=0, peak=0)

    def sample(self, W0r, g):
        m, rng = self.m, self.rng
        l1 = np.abs(W0r).sum() * m["dr"] * m["dpr"]
        n0 = int(round(self.nu * l1))
        prob = np.abs(W0r).ravel() / np.abs(W0r).sum()
        cells = rng.choice(prob.size, size=n0, p=prob)
        ir, jp = np.unravel_index(cells, W0r.shape)
        self.r = m["r"][ir] + (rng.random(n0) - 0.5) * m["dr"]
        self.pr = m["pr"][jp] + (rng.random(n0) - 0.5) * m["dpr"]
        self.eps = np.sign(W0r[ir, jp])
        self.P = SIG_P * rng.standard_normal(n0)
        self.C = rng.random((n0, self.nb, self.nq), dtype=np.float32)

    def rates(self):
        u = (self.r - self.rt[0]) / (self.rt[1] - self.rt[0])
        i0 = np.clip(np.floor(u).astype(int), 0, self.rt.size - 2)
        t = (u - i0)[:, None]
        return (1 - t) * self.kt[i0] + t * self.kt[i0 + 1]

    def stream(self, dt):
        rh = self.m["rh"]
        self.pr -= 0.5 * dt * dV_pt(self.r)
        self.r += self.pr / MU * dt
        self.r = (self.r + rh) % (2 * rh) - rh
        self.pr -= 0.5 * dt * dV_pt(self.r)

    def events(self, dt):
        k = self.rates().astype(np.float32)
        old = np.floor(self.C)
        self.C[:, 0] += k * dt
        if self.nb == 2:
            self.C[:, 1] -= k * dt
        dn = (np.floor(self.C) - old).astype(np.int64)        # (N, nb, NQ)
        w, b, q = np.nonzero(dn)
        if w.size == 0:
            return
        cnt = np.abs(dn[w, b, q])
        w, b, q = np.repeat(w, cnt), np.repeat(b, cnt), np.repeat(q, cnt)
        tau = self.eps[w] * np.sign(dn[w, b, q])
        xi = self.xi[q]
        if self.v == "W":
            dP = np.zeros_like(xi)
            dpr = xi
        else:
            dP = xi
            dpr = np.where(b == 0, 0.5 * xi, -0.5 * xi)
        nP = np.concatenate([self.P[w] + dP, self.P[w] - dP])
        npr = np.concatenate([self.pr[w] + dpr, self.pr[w] - dpr])
        nr = np.concatenate([self.r[w], self.r[w]])
        ne = np.concatenate([tau, -tau])
        self.stats["fires"] += int(w.size)
        self.stats["dP_viol"] += int(np.count_nonzero(dP))
        self.P = np.concatenate([self.P, nP])
        self.pr = np.concatenate([self.pr, npr])
        self.r = np.concatenate([self.r, nr])
        self.eps = np.concatenate([self.eps, ne])
        self.C = np.concatenate([self.C, self.rng.random((nP.size, self.nb, self.nq), dtype=np.float32)])

    def cell_keys(self):
        m = self.m
        ir = np.floor((self.r + m["rh"]) / m["dr"] + 0.5).astype(np.int64) % m["nr"]
        jp = np.floor(self.pr / m["dpr"] + 0.5).astype(np.int64) + (1 << 20)
        iP = np.floor(self.P / DP + 0.5).astype(np.int64) + (1 << 20)
        return (iP * (1 << 21) + jp) * m["nr"] + ir

    def sink(self):
        key = self.cell_keys()
        n = key.size
        o = np.lexsort((self.rng.random(n), self.eps, key))
        ks, es = key[o], self.eps[o]
        brk = np.flatnonzero(np.diff(ks)) + 1
        start = np.concatenate([[0], brk])
        cid = np.repeat(np.arange(start.size), np.diff(np.concatenate([start, [n]])))
        npos = np.bincount(cid, weights=(es > 0)).astype(np.int64)
        nneg = np.bincount(cid, weights=(es < 0)).astype(np.int64)
        mrem = np.minimum(npos, nneg)
        # rank within (cell, sign): negatives sort first in each cell
        first_neg = start
        first_pos = start + nneg
        rank = np.arange(n) - np.where(es < 0, first_neg[cid], first_pos[cid])
        kill = rank < mrem[cid]
        keep = np.ones(n, bool)
        keep[o[kill]] = False
        self.stats["sunk"] += int(kill.sum())
        for a in ("P", "r", "pr", "eps"):
            setattr(self, a, getattr(self, a)[keep])
        self.C = self.C[keep]

    def observe(self):
        m = self.m
        ir = np.floor((self.r + m["rh"]) / m["dr"] + 0.5).astype(np.int64) % m["nr"]
        jp = np.clip(np.floor(self.pr / m["dpr"] + 0.5).astype(np.int64) + m["npr"] // 2,
                     0, m["npr"] - 1)
        h = np.zeros(m["nr"] * m["npr"])
        np.add.at(h, ir * m["npr"] + jp, self.eps)
        Wrel = h.reshape(m["nr"], m["npr"]) / (self.nu * m["dr"] * m["dpr"])
        e, P2 = self.eps / self.nu, self.P**2
        Vv, Tv = V_pt(self.r), self.pr**2 / (2 * MU)
        cov = lambda a, b: float((e * a * b).sum() - (e * a).sum() * (e * b).sum())
        return dict(norm=float(e.sum()), V=float((e * Vv).sum()), T=float((e * Tv).sum()),
                    covPV=cov(P2, Vv), covPT=cov(P2, Tv), P2=float((e * P2).sum()),
                    P2_gross=float(P2.mean()), N=int(self.eps.size)), Wrel


def coarse(W, b=4):
    """Sum b x b blocks of a relative-marginal histogram (r, p_r)."""
    n0, n1 = W.shape
    return W.reshape(n0 // b, b, n1 // b, b).sum(axis=(1, 3))


def mesh_dict():
    m = Mesh2()
    return dict(r=m.r, pr=m.pr, dr=m.dr, dpr=m.dpr, nr=m.r.size, npr=m.pr.size, rh=10.0)


def run_C(job):
    variant, seed, nu, T, dt = job
    md = mesh_dict()
    W0r = wigner_mesh(psi0, md["r"], md["pr"])
    W0r /= W0r.sum() * md["dr"] * md["dpr"]
    sim = Banked(variant, nu, seed, md)
    sim.sample(W0r, None)
    t0 = time.perf_counter()
    o, Wrel = sim.observe()
    rows, snaps = [dict(t=0.0, **o)], {}
    every = int(round(C_EVERY / dt))
    for i in range(1, int(round(T / dt)) + 1):
        sim.stream(dt)
        sim.events(dt)
        sim.sink()
        sim.stats["peak"] = max(sim.stats["peak"], int(sim.eps.size))
        t = round(i * dt, 6)
        if i % every == 0:
            o, Wrel = sim.observe()
            rows.append(dict(t=t, **o))
            snaps[t] = Wrel
        if sim.eps.size > N_CAP:
            o, Wrel = sim.observe()
            rows.append(dict(t=t, capped=True, **o))
            snaps[t] = Wrel
            break
    return dict(variant=variant, seed=seed, nu=nu, rows=rows, snaps=snaps,
                stats=sim.stats, wall=time.perf_counter() - t0, W0r=W0r)


def part_C(seeds=1, nu=4000, T=4.0, T_b1=1.0, dt=0.01, workers=None, mref=None):
    import csv
    import multiprocessing as mp

    rule(f"Part C. Particle runs with integrator banks: nu = {nu}, t = 0 .. {T}, dt = {dt}, "
         f"{seeds} seed(s) per variant")
    if mref is None:                   # relative marginals: W at hbar, B1 at hbar/2 (A2)
        mref = {}
        md = mesh_dict()
        W0r = wigner_mesh(psi0, md["r"], md["pr"])
        W0r /= W0r.sum() * md["dr"] * md["dpr"]
        for v, hb in (("W", HBAR), ("B1", HBAR / 2)):
            m2, w, mref[v] = Mesh2(hbar_res=hb, dt=dt), W0r.copy(), {0.0: W0r.copy()}
            for i in range(1, int(round(T / dt)) + 1):
                w = m2.step(w)
                if i % int(round(C_EVERY / dt)) == 0:
                    mref[v][round(i * dt, 6)] = w.copy()
    jobs = [(v, s, nu, T if v == "W" else T_b1, dt) for s in range(seeds) for v in ("W", "B1")]
    workers = workers or min(len(jobs), os.cpu_count() or 1)
    out = []
    path = output_path(f"world_clocks_C_nu{nu}.csv")
    with open(path, "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["variant", "seed", "nu", "t", "norm", "V", "T", "E_rel", "covPV",
                     "covPT", "P2_signed", "P2_gross", "N", "L2c_vs_W_mesh", "L2c_vs_B1_mesh",
                     "fires", "sunk", "dP_violations", "peak_N", "wall_s"])
        with mp.Pool(workers) as pool:
            for r in pool.imap_unordered(run_C, jobs):
                W0r = r["W0r"]
                for row in r["rows"]:
                    t = row["t"]
                    lw = lb = float("nan")
                    if t in r["snaps"] and t in mref["W"]:
                        h = coarse(r["snaps"][t])
                        nc = np.linalg.norm(coarse(W0r))
                        lw = float(np.linalg.norm(h - coarse(mref["W"][t])) / nc)
                        lb = float(np.linalg.norm(h - coarse(mref["B1"][t])) / nc)
                    wr.writerow([r["variant"], r["seed"], nu, t, row["norm"], row["V"], row["T"],
                                 row["V"] + row["T"], row["covPV"], row["covPT"], row["P2"],
                                 row["P2_gross"], row["N"], lw, lb, r["stats"]["fires"],
                                 r["stats"]["sunk"], r["stats"]["dP_viol"], r["stats"]["peak"],
                                 round(r["wall"], 1)])
                    row.update(L2W=lw, L2B=lb)
                f.flush()
                fin = r["rows"][-1]
                print(f"  {r['variant']:2s} seed {r['seed']:2d}: t = {fin['t']}{' (capped)' if fin.get('capped') else ''}:"
                      f" <V> {fin['V']:+.4f}  cov(P^2,V) {fin['covPV']:+.3f}  coarse L2 to W mesh "
                      f"{fin['L2W']:.3f}, to B1 mesh {fin['L2B']:.3f}  N {fin['N']} (peak "
                      f"{r['stats']['peak']})  <P^2> signed {fin['P2']:.3f} gross {fin['P2_gross']:.2f}"
                      f"  firings {r['stats']['fires']}, daughters off the parent's P "
                      f"{2 * r['stats']['dP_viol']}  ({r['wall']:.0f} s)", flush=True)
                out.append({k: r[k] for k in ("variant", "seed", "rows", "stats", "wall")})
    print(f"  rows written to {os.path.basename(path)}")
    print("  population N(t), seed 0:")
    for r in sorted([x for x in out if x["seed"] == 0], key=lambda x: x["variant"], reverse=True):
        print(f"    {r['variant']:2s} " + "  ".join(f"t={x['t']:.2f}: {x['N']}" for x in r["rows"]))
    for v in ("W", "B1"):
        fin = [r["rows"][-1] for r in out if r["variant"] == v]
        if not fin:
            continue
        def ms(key):
            a = np.array([x[key] for x in fin])
            return a.mean(), (a.std(ddof=1) / np.sqrt(a.size) if a.size > 1 else float("nan"))
        print(f"  {v:2s} mean over {len(fin)} seeds at t = {T}:  <V> {ms('V')[0]:+.4f} +- {ms('V')[1]:.4f}"
              f"   cov(P^2,V) {ms('covPV')[0]:+.4f} +- {ms('covPV')[1]:.4f}"
              f"   L2 to W mesh {ms('L2W')[0]:.4f} +- {ms('L2W')[1]:.4f}"
              f"   L2 to B1 mesh {ms('L2B')[0]:.4f} +- {ms('L2B')[1]:.4f}")
    return out


# ======================================================================
# Part F -- figures
# ======================================================================
def save_fig(fig, name):
    fig.savefig(output_path(name), dpi=150, bbox_inches="tight")
    dp = docs_path(name)
    if dp:
        fig.savefig(dp, dpi=150, bbox_inches="tight")
    print(f"  wrote {name}")


def part_F(nu_c=None):
    import glob
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rule("Part F. Figures")
    m, h, t_, r_ = (load(f"world_clocks_{k}.json") for k in "MHTR")
    cfiles = sorted(glob.glob(output_path("world_clocks_C_nu*.json")),
                    key=lambda f: int(f.split("_nu")[-1].split(".")[0]))
    allc = []
    for fn in cfiles:
        with open(fn) as f:
            allc.append((int(fn.split("_nu")[-1].split(".")[0]), json.load(f)))
    if nu_c is not None:
        cfiles = [f for f in cfiles if f.endswith(f"_nu{nu_c}.json")]
    c = None
    if cfiles:
        with open(cfiles[-1]) as f:
            c = json.load(f)
        nu_used = int(cfiles[-1].split("_nu")[-1].split(".")[0])
    fig, axs = plt.subplots(2, 3, figsize=(15.5, 8.4))
    ax = np.array([[axs[0, 0], axs[0, 1]], [axs[1, 1], axs[0, 2]]])
    col = dict(W="#1f6f8b", B1="#c0392b", B2="#e67e22", cl="#7f7f7f")
    if m and h:
        tt = np.array(m["times"])
        a = ax[0, 0]
        for v, lab in (("W", "W (world bank)"), ("B1", "B1 (particle banks, whole V)"),
                       ("B2", "B2 (particle banks, V split)"), ("cl", "classical streaming only")):
            a.plot(tt, [rw["rel_drift"] for rw in m["runs"][v]], color=col[v], label=lab)
        wl = {float(k): v for k, v in m["wig_l2"].items()}
        ts = sorted(wl)
        a.plot(ts, [wl[t][3] for t in ts], "o", color=col["B1"], mfc="none",
               label="exact: QM at hbar/2 (B1's prediction)")
        a.axhline(0, color="k", lw=0.6)
        a.set_xlabel("t")
        a.set_ylabel(r"$\|W_{\rm rel}(t) - W_{\psi_0}\|_2 / \|W_{\psi_0}\|_2$")
        a.set_title("(a) relative motion of the ground state (exact QM: 0)")
        a.legend(fontsize=7.5, loc="upper left")
        a = ax[0, 1]
        for v in ("W", "B1", "B2"):
            a.semilogy(tt[1:], [max(rw["sep"], 1e-16) for rw in m["runs"][v][1:]], color=col[v], label=v)
        a.set_xlabel("t")
        a.set_ylabel(r"$\|W - g_P \otimes W_{\rm rel}\|_2 / \|W\|_2$")
        a.set_title("(b) centre of mass and relative motion (exact QM: 0)")
        a.legend(fontsize=8)
    if c:
        a = ax[1, 0]
        for v in ("W", "B1"):
            runs = [r for r in c if r["variant"] == v]
            for k, rr in enumerate(runs):
                a.semilogy([x["t"] for x in rr["rows"]], [x["N"] / nu_used for x in rr["rows"]],
                           color=col[v], alpha=0.5, label=(f"{v}: worlds" if k == 0 else None))
                a.semilogy([x["t"] for x in rr["rows"]], [x["P2_gross"] for x in rr["rows"]], "--",
                           color=col[v], alpha=0.5, label=(f"{v}: gross <P^2>" if k == 0 else None))
        a.axhline(1.0, color="k", lw=0.6, ls=":")
        a.set_xlabel("t")
        a.set_ylabel(r"worlds per unit $\nu$;  unsigned $\langle P^2\rangle$")
        a.set_title(f"(e) cost of the banks (nu = {nu_used}, {len(c) // 2} seeds)")
        a.legend(fontsize=7.5)
        # (d) and (f): signed observables, particles (mean +- s.e.) against the mesh
        for a, key, mkey, ttl in ((axs[1, 0], "covPV", "corr", "(d) cov(P^2, v): particles against the mesh"),
                                  (axs[1, 2], "V", "V", "(f) <v>: particles against the mesh")):
            for v in ("W", "B1"):
                for j, (nu_j, cj) in enumerate(allc):
                    runs = [r for r in cj if r["variant"] == v]
                    nr_ = min(len(r["rows"]) for r in runs)
                    tp = np.array([runs[0]["rows"][k]["t"] for k in range(nr_)])
                    arr = np.array([[r["rows"][k][key] for k in range(nr_)] for r in runs])
                    a.errorbar(tp + 0.03 * j, arr.mean(0), yerr=arr.std(0, ddof=1) / np.sqrt(len(runs)),
                               fmt="os"[j % 2], ms=3.5, color=col[v], capsize=2,
                               mfc=(col[v] if j == 0 else "none"),
                               label=f"{v} particles, nu = {nu_j:.0e}")
                if m:
                    a.plot(m["times"], [rw[mkey] for rw in m["runs"][v]], color=col[v], lw=1,
                           label=f"{v} mesh")
            a.set_xlabel("t")
            a.set_title(ttl, fontsize=10)
            a.legend(fontsize=7.5)
        axs[1, 0].set_ylim(-2.0, 1.0)
    if t_:
        a = ax[1, 1]
        reach = sorted(int(k) for k in t_["mesh"])
        a.plot(reach, [t_["mesh"][str(k)]["W"] for k in reach], "o-", color=col["W"], label="W mesh")
        a.plot(reach, [t_["mesh"][str(k)]["B1"] for k in reach], "s-", color=col["B1"], label="B1 mesh")
        a.axhline(t_["T_q"], color=col["W"], ls="--", lw=0.8, label=f"exact QM {t_['T_q']:.4f}")
        a.axhline(t_["T_h"], color=col["B1"], ls="--", lw=0.8, label=f"exact at hbar/2 {t_['T_h']:.4f}")
        a.axhline(t_["mesh"][str(reach[0])]["cl"], color=col["cl"], ls=":", lw=0.8, label="classical")
        a.set_xticks(reach)
        a.set_xticklabels([f"{k}$\\pi$" for k in reach])
        a.set_xlabel("reach $y_{\\max}$")
        a.set_ylabel("transmission")
        a.set_title("(c) Eckart pair collision")
        a.legend(fontsize=7.5)
    fig.tight_layout()
    save_fig(fig, "world_clocks.png")
    plt.close(fig)
    schematic(plt)


def schematic(plt):
    """VLBI against WPMW: stations and baselines, worlds and chords."""
    from matplotlib.patches import Circle, FancyArrowPatch

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
    a = ax[0]
    st = np.array([[0.1, 0.25], [0.45, 0.12], [0.85, 0.3], [0.62, 0.62]])
    names = ["station 1", "station 2", "station 3", "station 4"]
    for i in range(4):
        for j in range(i + 1, 4):
            a.plot(*st[[i, j]].T, color="#c0392b", lw=1.0, alpha=0.8)
    for (x, y), nm in zip(st, names):
        a.add_patch(Circle((x, y), 0.035, color="#1f6f8b"))
        a.text(x, y - 0.075, nm + "\nown maser clock", ha="center", va="top", fontsize=8)
    a.add_patch(Circle((0.3, 0.86), 0.06, color="#e67e22", alpha=0.85))
    a.text(0.3, 0.86, "corr-\nelator", ha="center", va="center", fontsize=8)
    for x, y in st:
        a.add_patch(FancyArrowPatch((x, y), (0.3, 0.86), arrowstyle="-|>", mutation_scale=8,
                                    color="#7f7f7f", lw=0.6, alpha=0.6))
    a.text(0.62, 0.95, "6 baselines = 6 visibilities\nclock offsets cancel in\nclosure phases",
           fontsize=8, va="top")
    a.set_xlim(-0.05, 1.05)
    a.set_ylim(-0.08, 1.02)
    a.set_aspect("equal")
    a.axis("off")
    a.set_title("VLBI: stations, clocks, baselines, correlator", fontsize=10)
    a = ax[1]
    a.axhline(0, color="#999", lw=0.5)
    a.axvline(0, color="#999", lw=0.5)
    A = np.array([-0.35, 0.15])
    B = np.array([0.45, -0.25])
    mid = (A + B) / 2
    a.plot(*np.array([A, B]).T, color="#c0392b", lw=1.6)
    for pt, lab, dy, va in ((A, "world X - Y\n(its clock)", 0.08, "bottom"),
                            (B, "world X + Y\n(its clock)", -0.08, "top")):
        a.add_patch(Circle(pt, 0.045, color="#1f6f8b"))
        a.text(pt[0], pt[1] + dy, lab, ha="center", va=va, fontsize=8)
    a.add_patch(Circle(mid, 0.03, color="#e67e22"))
    a.annotate("reader at the midpoint X:\none bank for the bond (1, 2)",
               xy=mid, xytext=(-0.85, 0.55), fontsize=8,
               arrowprops=dict(arrowstyle="-", color="#e67e22", lw=0.8))
    a.plot([-0.75, 0.75], [0.75, -0.75], ls=":", color="#7f7f7f", lw=0.8)
    a.text(0.42, -0.72, "relative direction\ny1 - y2 (read)", fontsize=7.5, color="#555")
    a.plot([-0.75, 0.75], [-0.75, 0.75], ls=":", color="#bbb", lw=0.8)
    a.text(0.12, 0.62, "centre-of-mass direction\n(dark for a pair potential)", fontsize=7.5, color="#777")
    a.annotate("", xy=(B[0], -0.9), xytext=(A[0], -0.9),
               arrowprops=dict(arrowstyle="<->", color="#555", lw=0.8))
    a.text(mid[0], -0.97, "particle 1's chord 2 y1: all a particle bank reads",
           ha="center", fontsize=7.5, color="#555")
    a.set_xlim(-0.95, 0.95)
    a.set_ylim(-1.05, 0.95)
    a.set_xticks([])
    a.set_yticks([])
    a.set_xlabel("x1")
    a.set_ylabel("x2")
    a.set_aspect("equal")
    a.set_title("WPMW: worlds, clocks, a chord, a reader", fontsize=10)
    fig.tight_layout()
    save_fig(fig, "world_clocks_vlbi.png")
    plt.close(fig)


# ======================================================================
# main
# ======================================================================
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("parts", nargs="*", default=["A", "K", "H", "M", "R", "T", "F"])
    ap.add_argument("--T", type=float, default=4.0, help="mesh run length (Part M)")
    ap.add_argument("--dt", type=float, default=0.01)
    ap.add_argument("--seeds", type=int, default=1, help="seeds per variant (Part C)")
    ap.add_argument("--nu", type=int, default=4000, help="worlds per unit ||W||_1 (Part C)")
    ap.add_argument("--Tc", type=float, default=4.0, help="world-bank run length (Part C)")
    ap.add_argument("--Tb1", type=float, default=1.0, help="particle-bank run length (Part C)")
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--fig-nu", type=int, default=None, help="Part C file the figure uses")
    args = ap.parse_args()
    parts = [p.upper() for p in args.parts]
    t0 = time.perf_counter()
    if "A" in parts:
        dump(part_A(), "world_clocks_A.json")
    if "K" in parts:
        dump(part_K(), "world_clocks_K.json")
    if "M" in parts or "H" in parts:
        times = [round(0.05 * i, 6) for i in range(int(round(args.T / 0.05)) + 1)]
        h = part_H(times)
        dump({**h, "half": {k: v for k, v in h["half"].items() if k != "wig"}}, "world_clocks_H.json")
    if "M" in parts:
        m = part_M(T=args.T, dt=args.dt, every=int(round(0.05 / args.dt)))
        m.update(report_MH(m, h))
        dump({k: v for k, v in m.items() if k not in ("W0r", "snaps")}, "world_clocks_M.json")
    if "R" in parts:
        dump(part_R(T=args.T, dt=args.dt), "world_clocks_R.json")
    if "T" in parts:
        dump(part_T(), "world_clocks_T.json")
    if "C" in parts:
        c = part_C(seeds=args.seeds, nu=args.nu, T=args.Tc, T_b1=args.Tb1, dt=args.dt,
                   workers=args.workers)
        dump(c, f"world_clocks_C_nu{args.nu}.json")
    if "F" in parts:
        part_F(nu_c=args.fig_nu)
    print(f"\nwall time {time.perf_counter() - t0:.1f} s")


if __name__ == "__main__":
    main()
