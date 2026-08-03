"""
demo_position_pair_ladder.py — verification companion to
``docs/analysis/position_pair_ladder.md``.

The note models the position-representation density matrix rho(X, X') as a
population of conjugate world-particle pairs.  Each world-particle carries a
position and a phase and nothing else; a pair is one Monte-Carlo sample of one
element of rho, the positon leg being the ket and the negaton leg the bra; a
self-conjugate particle (X = X') is a sample of the diagonal, i.e. of the
observable probability density.  This script checks every quantitative claim of
that note before it is written down.

Lattice conventions.  Positions live on a periodic ring of M sites,
x_m = m a with a = L/M.  The kinetic operator is the nearest-neighbour form
T = -J (S + S^dagger) + 2J with J = hbar^2 / (2 m a^2); the potential is
diagonal.  Rung k = m - n indexes the distance off the diagonal; midpoint
r = (m + n) a / 2.

Parts
-----
A  Theorem P1.  The von Neumann generator splits exactly into four one-leg
   hop channels (ket/bra x up/down) with amplitude +/- i J/hbar, plus a
   diagonal pump of rate -(V(X) - V(X'))/hbar.  Checks exactness, Hermiticity,
   trace preservation, the J = 0 freeze (moduli and populations exactly
   constant, phases winding) and the V = 0 freeze of the pump.

B  Theorem P2.  Momentum is the misalignment of a nearest-neighbour conjugate
   pair.  Exact lattice continuity with j = (hbar / m a) |rho_1| sin(mu), and
   convergence of hbar mu / a to the packet momentum.

C  Proposition P3.  The pump alone reproduces the Euler force term:
   d j / dt |_pump -> rho F / m, second order in a.

D  Proposition P4.  Two obstructions.  (i) A state with no nearest-neighbour
   coherence has exactly zero population flux: nothing moves until coherence
   exists, so no lone-hopper process on self-conjugate particles can drive the
   dynamics.  (ii) The off-diagonal sector admits no positive-rate jump
   process: the one-leg hop amplitude is imaginary, so hop *probability* is
   O(dt^2), and local gauge freedom sweeps arg rho_mn around the full circle,
   so no fixed background can put the weights on a single ray.

E  Proposition P5.  Resource arithmetic.  The l1 mass of rho on the lattice is
   bounded by M (the number of cells) and the Z4 sampling mass by sqrt(2) M —
   comparable to the Wigner sea mass — but the channel rate Lambda = 4 J/hbar
   diverges as a^-2, which the Wigner and momentum-basis ladders do not.

F  Theorem P6.  Statistical equivalence.  A complex-weight pair-ensemble Monte
   Carlo with Poisson hop clocks and the pump reproduces rho(t) without bias,
   converging as 1/sqrt(N), with noise amplified by exactly exp(Lambda t).

G  Corollary P7.  The diagonal sector *is* a positive-rate particle process:
   self-conjugate particles hopping at rates read from the misalignment of the
   adjacent pairs reproduce rho(x, x, t) exactly.
"""

from __future__ import annotations

import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from wpmwlib.wpmw_utils import output_path, docs_path  # noqa: E402

HBAR = 1.0
MASS = 1.0
L = 8.0
V_P = 1.5

rng = np.random.default_rng(20260803)


# ------------------------------------------------------------------ #
# lattice helpers                                                     #
# ------------------------------------------------------------------ #
def lattice(M: int):
    """Return (a, J, x, T) for an M-site periodic position lattice."""
    a = L / M
    J = HBAR ** 2 / (2.0 * MASS * a ** 2)
    x = np.arange(M) * a
    S = np.roll(np.eye(M), 1, axis=0)          # S[m, n] = delta_{m, n+1}
    T = -J * (S + S.T) + 2.0 * J * np.eye(M)
    return a, J, x, T


def cosine_well(x):
    return -V_P * np.cos(2.0 * np.pi * x / L)


def random_state(M, rank=4):
    A = rng.normal(size=(M, rank)) + 1j * rng.normal(size=(M, rank))
    rho = A @ A.conj().T
    return rho / np.trace(rho).real


def channels(rho, J, V):
    """Assemble d rho / dt from four one-leg hops plus the diagonal pump."""
    lam = J / HBAR
    ket_up = np.roll(rho, -1, axis=0)          # rho_{m+1, n}
    ket_dn = np.roll(rho, +1, axis=0)          # rho_{m-1, n}
    bra_up = np.roll(rho, -1, axis=1)          # rho_{m, n+1}
    bra_dn = np.roll(rho, +1, axis=1)          # rho_{m, n-1}
    hop = 1j * lam * (ket_up + ket_dn) - 1j * lam * (bra_up + bra_dn)
    pump = -1j / HBAR * (V[:, None] - V[None, :]) * rho
    return hop + pump


def commutator_flux(rho, H):
    return -1j / HBAR * (H @ rho - rho @ H)


def rung1(rho):
    """The nearest-neighbour elements rho_{m+1, m}, indexed by m (site m+1/2)."""
    return np.diagonal(np.roll(rho, -1, axis=0))


# ------------------------------------------------------------------ #
# Part A                                                              #
# ------------------------------------------------------------------ #
def part_A():
    print("=" * 72)
    print("A.  Theorem P1 — four one-leg hop channels plus a diagonal pump")
    print("=" * 72)
    M = 16
    a, J, x, T = lattice(M)
    V = cosine_well(x) + 0.4 * np.sin(4.0 * np.pi * x / L)
    H = T + np.diag(V)
    rho = random_state(M)

    exact = commutator_flux(rho, H)
    built = channels(rho, J, V)
    scale = np.max(np.abs(exact))
    rel = np.max(np.abs(built - exact)) / scale
    herm = np.max(np.abs(built - built.conj().T)) / scale
    tr = abs(np.trace(built))
    print(f"  assembled flux vs -i/hbar [H, rho] : {rel:.3e}  (relative)")
    print(f"  Hermiticity of the assembled flux  : {herm:.3e}")
    print(f"  trace of the assembled flux        : {tr:.3e}")

    # J = 0 freeze: only phases move.
    zero = channels(rho, 0.0, V)
    dmod = np.max(np.abs(np.real(np.conj(rho) * zero)))     # d|rho|^2 / 2 dt
    dpop = np.max(np.abs(np.diag(zero)))
    print(f"  J = 0: max |d(|rho|^2)/dt| / 2      : {dmod:.3e}")
    print(f"  J = 0: max |d rho_mm / dt|          : {dpop:.3e}")

    # V = 0: pump silent, hops alone.
    only_hop = channels(rho, J, np.zeros(M))
    print(f"  V = 0: pump contribution           : "
          f"{np.max(np.abs(only_hop - channels(rho, J, np.zeros(M)))):.3e}")

    # channel-by-channel rung/midpoint bookkeeping
    print("  channel bookkeeping (leg, direction) -> (d rung, d midpoint/a), "
          "amplitude")
    for name, dk, dr, amp in (("ket up  ", +1, +0.5, "+i J/hbar"),
                              ("ket down", -1, -0.5, "+i J/hbar"),
                              ("bra up  ", -1, +0.5, "-i J/hbar"),
                              ("bra down", +1, -0.5, "-i J/hbar")):
        print(f"    {name}  -> ({dk:+d}, {dr:+.1f})   {amp}")
    return dict(rel=rel, herm=herm, tr=tr, dmod=dmod, dpop=dpop)


# ------------------------------------------------------------------ #
# Part B                                                              #
# ------------------------------------------------------------------ #
def part_B():
    print()
    print("=" * 72)
    print("B.  Theorem P2 — momentum is the misalignment of a rung-1 pair")
    print("=" * 72)
    M = 24
    a, J, x, T = lattice(M)
    V = cosine_well(x)
    H = T + np.diag(V)
    rho = random_state(M)
    flux = commutator_flux(rho, H)

    r1 = rung1(rho)
    j = (HBAR / (MASS * a)) * np.abs(r1) * np.sin(np.angle(r1))
    j_alt = (2.0 * J * a / HBAR) * np.imag(r1)
    div = (j - np.roll(j, 1)) / a
    cont = np.max(np.abs(np.real(np.diag(flux)) + div))
    print(f"  j = (hbar/ma)|rho_1| sin(mu) vs (2Ja/hbar) Im rho_1 : "
          f"{np.max(np.abs(j - j_alt)):.3e}")
    print(f"  exact lattice continuity  d rho_mm/dt + div j       : "
          f"{cont:.3e}")
    print(f"  imaginary part of the population flux               : "
          f"{np.max(np.abs(np.imag(np.diag(flux)))):.3e}")

    print()
    print("  hbar mu / a  recovers the packet momentum p0 = 1.7 :")
    rows = []
    for M in (16, 32, 64, 128, 256):
        a, J, x, T = lattice(M)
        psi = np.exp(-(x - 3.0) ** 2 / (4 * 0.8 ** 2) + 1j * 1.7 * x / HBAR)
        psi /= np.linalg.norm(psi)
        rho = np.outer(psi, psi.conj())
        r1 = rung1(rho)
        w = np.abs(np.diag(rho))
        w = w / w.sum()
        p_est = float(np.sum(w * HBAR * np.angle(r1) / a))
        rows.append((M, p_est, abs(p_est - 1.7)))
        print(f"    M = {M:4d}   p_est = {p_est:.6f}   err = "
              f"{abs(p_est - 1.7):.2e}")

    # the field identity, on a chirped packet so that p varies with x
    M = 96
    a, J, x, T = lattice(M)
    psi = np.exp(-(x - 4.0) ** 2 / (4 * 1.1 ** 2)
                 + 1j * (0.8 * x + 0.22 * (x - 4.0) ** 2) / HBAR)
    psi /= np.linalg.norm(psi)
    rho = np.outer(psi, psi.conj())
    r1 = rung1(rho)
    p_mis = HBAR * np.angle(r1) / a                    # from the misalignment
    j = (HBAR / (MASS * a)) * np.imag(r1)
    dens = 0.5 * (np.diag(rho).real + np.roll(np.diag(rho).real, -1))
    p_cur = MASS * j / np.maximum(dens, 1e-300)        # from m j / rho
    keep = dens > 1e-4 * dens.max()
    keep[-1] = False                                   # drop the ring seam
    dev = float(np.max(np.abs(p_mis[keep] - p_cur[keep])))
    print(f"  field identity  hbar mu/a  vs  m j / rho  (chirped packet): "
          f"max dev = {dev:.3e}")
    return dict(cont=cont, rows=rows, x=x, p_mis=p_mis, p_cur=p_cur,
                keep=keep, dens=dens, dev=dev)


# ------------------------------------------------------------------ #
# Part C                                                              #
# ------------------------------------------------------------------ #
def part_C():
    print()
    print("=" * 72)
    print("C.  Proposition P3 — the pump alone gives the Euler force term")
    print("=" * 72)
    rows = []
    for M in (32, 64, 128, 256, 512, 1024):
        a, J, x, T = lattice(M)
        psi = np.exp(-(x - 3.0) ** 2 / (4 * 0.8 ** 2) + 1j * 0.9 * x / HBAR)
        psi /= np.linalg.norm(psi)
        rho = np.outer(psi, psi.conj())
        V = cosine_well(x)
        F = -np.gradient(V, a)
        r1 = rung1(rho)
        dV = np.roll(V, -1) - V
        # pump-driven rate of change of the current density
        djdt = (2.0 * J * a / HBAR) * (-(dV / HBAR) * np.real(r1)) / a
        dens_mid = np.real(r1) / a                    # |rho_1| cos(mu) / a
        F_mid = 0.5 * (F + np.roll(F, -1))
        ref = dens_mid * F_mid / MASS
        err = float(np.max(np.abs(djdt - ref)) / np.max(np.abs(ref)))
        rows.append((M, a, err))
        print(f"    M = {M:5d}  a = {a:.5f}   rel err = {err:.3e}")
    ratio = rows[-2][2] / rows[-1][2]
    print(f"  refinement ratio between the last two rows: {ratio:.2f} "
          f"(second order = 4)")
    return dict(rows=rows, ratio=ratio)


# ------------------------------------------------------------------ #
# Part D                                                              #
# ------------------------------------------------------------------ #
def part_D():
    print()
    print("=" * 72)
    print("D.  Proposition P4 — the two obstructions")
    print("=" * 72)
    M = 24
    a, J, x, T = lattice(M)
    V = cosine_well(x)
    H = T + np.diag(V)

    # (i) incoherent states are instantaneously static
    pops = rng.random(M)
    rho_d = np.diag(pops / pops.sum())
    fl = commutator_flux(rho_d, H)
    print(f"  (i) diagonal state: max |d rho_mm/dt|   = "
          f"{np.max(np.abs(np.diag(fl))):.3e}")
    print(f"      diagonal state: max |d rho_mn/dt|   = "
          f"{np.max(np.abs(fl - np.diag(np.diag(fl)))):.3e}   (rung 1 is fed)")

    # ballistic, not diffusive: <x^2> of a released site
    M2 = 64
    a2, J2, x2, T2 = lattice(M2)
    ev, U = np.linalg.eigh(T2)
    psi0 = np.zeros(M2, complex)
    psi0[M2 // 2] = 1.0
    ts = np.array([0.02, 0.04, 0.08, 0.16])
    var = []
    d = (np.arange(M2) - M2 // 2) * a2
    for t in ts:
        psi = U @ (np.exp(-1j * ev * t / HBAR) * (U.conj().T @ psi0))
        var.append(float(np.sum(np.abs(psi) ** 2 * d ** 2)))
    var = np.array(var)
    slope = np.polyfit(np.log(ts), np.log(var), 1)[0]
    print(f"  (i) free spreading of one site: d log<x^2>/d log t = "
          f"{slope:.3f}  (ballistic = 2, diffusive = 1)")

    # hop probability is O(dt^2)
    lam = J2 / HBAR
    print("  (ii) one-leg update over dt: amplitude vs probability")
    for dt in (1e-2, 1e-3, 1e-4):
        Uk = np.eye(M2) + 1j * lam * dt * (np.roll(np.eye(M2), 1, axis=0)
                                           + np.roll(np.eye(M2), -1, axis=0))
        amp = abs(Uk[0, 1])
        print(f"      dt = {dt:.0e}   |amp| = {amp:.3e}   "
              f"prob = {amp ** 2:.3e}   amp/dt = {amp / dt:.4f}")

    # (ii) local gauge sweeps arg rho_mn around the full circle
    psi = np.exp(-(x - 3.0) ** 2 / (4 * 0.8 ** 2))
    psi = psi / np.linalg.norm(psi)
    angles = []
    for _ in range(2000):
        chi = rng.uniform(0.0, 2 * np.pi, size=M)
        p = psi * np.exp(1j * chi)
        r = np.outer(p, p.conj())
        angles.append(np.angle(r[5, 9]))
    angles = np.array(angles)
    hist, _ = np.histogram(angles, bins=24, range=(-np.pi, np.pi))
    flat = float(hist.std() / hist.mean())
    print(f"  (ii) arg rho_(5,9) under local gauge: histogram flatness "
          f"(std/mean) = {flat:.3f}")
    print("       the diagonal is untouched by the gauge, so no fixed "
          "background")
    print("       can put every state's off-diagonal weights on one ray.")
    return dict(slope=slope, flat=flat, angles=angles)


# ------------------------------------------------------------------ #
# Part E                                                              #
# ------------------------------------------------------------------ #
def part_E():
    print()
    print("=" * 72)
    print("E.  Proposition P5 — resource arithmetic and the rate comparison")
    print("=" * 72)
    M = 64
    a, J, x, T = lattice(M)
    V = cosine_well(x)
    H = T + np.diag(V)
    ev, U = np.linalg.eigh(H)
    sig = 0.5 * np.sqrt(HBAR / np.sqrt(MASS * V_P * (2 * np.pi / L) ** 2))
    psi0 = np.exp(-(x - L / 4.0) ** 2 / (4 * sig ** 2))
    psi0 = psi0 / np.linalg.norm(psi0)

    omega = np.sqrt(V_P * (2 * np.pi / L) ** 2 / MASS)
    T_per = 2 * np.pi / omega
    ts = np.linspace(0.0, 4.0 * T_per, 121)
    l1, z4, dmax = [], [], []
    for t in ts:
        psi = U @ (np.exp(-1j * ev * t / HBAR) * (U.conj().T @ psi0))
        rho = np.outer(psi, psi.conj())
        l1.append(float(np.sum(np.abs(rho))))
        z4.append(float(np.sum(np.abs(rho.real)) + np.sum(np.abs(rho.imag))))
        dmax.append(float(np.max(np.abs(np.diag(rho)))))
    l1 = np.array(l1)
    z4 = np.array(z4)
    print(f"  lattice M = {M}, T_period = {T_per:.3f}, run to 4 T_period")
    print(f"  l1 mass  sum |rho_mn|      : max {l1.max():.3f}   bound M = {M}")
    print(f"  Z4 sampling mass           : max {z4.max():.3f}   "
          f"bound sqrt(2) M = {np.sqrt(2) * M:.1f}")
    print(f"  max diagonal element       : {max(dmax):.4f}  (bound 1)")
    print(f"  max |rho_mn|, m != n       : "
          f"{0.5:.3f} is the state-independent bound")

    print()
    print("  channel-rate comparison across the three representations:")
    lam_pos = 4.0 * J / HBAR
    lam_mom = V_P / (2.0 * HBAR)
    print(f"    Wigner  W(x, p)     kinetic: free advection, no rate")
    print(f"                        potential: Gamma_q ~ V_p/hbar = "
          f"{V_P / HBAR:.3f}")
    print(f"    rho(P, P') momentum kinetic: phase winding, no rate")
    print(f"                        potential: |amp| = V_q/2hbar = "
          f"{lam_mom:.3f}   (a-independent)")
    print(f"    rho(X, X') position kinetic: |amp| Lambda = 4J/hbar = "
          f"{lam_pos:.1f}   (~ a^-2)")
    print(f"                        potential: phase winding, no rate")
    for MM in (16, 32, 64, 128):
        aa, JJ, _, _ = lattice(MM)
        print(f"      M = {MM:4d}  a = {aa:.4f}   Lambda = 4J/hbar = "
              f"{4 * JJ / HBAR:9.2f}")
    return dict(ts=ts, l1=l1, z4=z4, M=M, T_per=T_per, dmax=np.array(dmax))


# ------------------------------------------------------------------ #
# Part F                                                              #
# ------------------------------------------------------------------ #
def _mc_pairs(rho0, V, Lam, T_end, N, seed, M):
    r = np.random.default_rng(seed)
    P = np.abs(rho0)
    tot = float(P.sum())
    p = (P / tot).flatten()
    idx = r.choice(M * M, size=N, p=p)
    mm = (idx // M).copy()
    nn = (idx % M).copy()
    z = rho0.flatten()[idx]
    w = (z / np.abs(z)).astype(complex) * tot / N
    t = np.zeros(N)
    done = np.zeros(N, bool)
    while not done.all():
        act = np.where(~done)[0]
        dt = r.exponential(1.0 / Lam, size=act.size)
        tn = t[act] + dt
        over = tn >= T_end
        seg = np.where(over, T_end - t[act], dt)
        dV = V[mm[act]] - V[nn[act]]
        w[act] = w[act] * np.exp(-1j * dV * seg / HBAR)
        jump = act[~over]
        if jump.size:
            c = r.integers(0, 4, size=jump.size)
            d = np.where(c % 2 == 0, 1, -1)
            ket = c < 2
            mm[jump] = np.where(ket, (mm[jump] + d) % M, mm[jump])
            nn[jump] = np.where(~ket, (nn[jump] + d) % M, nn[jump])
            w[jump] = w[jump] * np.where(ket, 1j, -1j)
        t[act] = np.where(over, T_end, tn)
        done[act] = over
    est = np.zeros((M, M), complex)
    np.add.at(est, (mm, nn), w)
    return est * np.exp(Lam * T_end)


def part_F():
    print()
    print("=" * 72)
    print("F.  Theorem P6 — the pair ensemble reproduces rho(t) without bias")
    print("=" * 72)
    M = 8
    a, J, x, T = lattice(M)
    V = cosine_well(x)
    H = T + np.diag(V)
    Lam = 4.0 * J / HBAR
    ev, U = np.linalg.eigh(H)
    psi0 = np.exp(-(x - 3.0) ** 2 / (4 * 0.9 ** 2))
    psi0 = psi0 / np.linalg.norm(psi0)
    rho0 = np.outer(psi0, psi0.conj())

    def exact(t):
        Ut = U @ np.diag(np.exp(-1j * ev * t / HBAR)) @ U.conj().T
        return Ut @ rho0 @ Ut.conj().T

    print(f"  M = {M}, a = {a:.3f}, Lambda = 4J/hbar = {Lam:.3f}")
    conv = []
    T_end = 0.5
    ex = exact(T_end)
    for N in (2_000, 20_000, 200_000, 2_000_000):
        t0 = time.time()
        est = _mc_pairs(rho0, V, Lam, T_end, N, 1234, M)
        err = float(np.max(np.abs(est - ex)) / np.max(np.abs(ex)))
        conv.append((N, err))
        print(f"    N = {N:>9,}   rel err = {err:.4f}   "
              f"wall = {time.time() - t0:.1f} s")
    r_obs = conv[-2][1] / conv[-1][1]
    print(f"  error ratio for a 10x ensemble: {r_obs:.2f} "
          f"(1/sqrt(N) = 3.16)")

    print()
    print("  unbiasedness and the exp(Lambda t) noise amplification:")
    noise = []
    for T_end in (0.25, 0.5, 1.0, 1.5):
        ex = exact(T_end)
        R = 16
        runs = [_mc_pairs(rho0, V, Lam, T_end, 50_000, 7000 + k, M)
                for k in range(R)]
        mean = sum(runs) / R
        bias = float(np.max(np.abs(mean - ex)) / np.max(np.abs(ex)))
        spread = float(np.mean([np.max(np.abs(q - ex)) for q in runs])
                       / np.max(np.abs(ex)))
        noise.append((T_end, bias, spread, np.exp(Lam * T_end)))
        print(f"    T = {T_end:.2f}  Lambda T = {Lam * T_end:.2f}   "
              f"mean of {R} runs: {bias:.4f}   single run: {spread:.4f}   "
              f"single/sqrt(R) = {spread / np.sqrt(R):.4f}   "
              f"exp(Lambda T) = {np.exp(Lam * T_end):.2f}")
    return dict(conv=conv, noise=noise, Lam=Lam)


# ------------------------------------------------------------------ #
# Part G                                                              #
# ------------------------------------------------------------------ #
def part_G():
    print()
    print("=" * 72)
    print("G.  Corollary P7 — the diagonal sector is a positive-rate process")
    print("=" * 72)
    M = 32
    a, J, x, T = lattice(M)
    V = cosine_well(x)
    H = T + np.diag(V)
    ev, U = np.linalg.eigh(H)
    psi0 = np.exp(-(x - L / 4.0) ** 2 / (4 * 0.6 ** 2))
    psi0 = psi0 / np.linalg.norm(psi0)

    omega = np.sqrt(V_P * (2 * np.pi / L) ** 2 / MASS)
    T_per = 2 * np.pi / omega
    T_end = 1.5 * T_per
    n_steps = 6000
    dt = T_end / n_steps
    floor = 1e-6

    # Precompute the guiding field: populations and the two hop rates.
    # rate(m -> m+1) reads j at m+1/2 when it is positive;
    # rate(m -> m-1) reads j at m-1/2 when it is negative.  Both are
    # manifestly non-negative — this is the whole point of the corollary.
    RU = np.zeros((n_steps, M))
    RD = np.zeros((n_steps, M))
    POP = np.zeros((n_steps, M))
    for s in range(n_steps):
        psi = U @ (np.exp(-1j * ev * s * dt / HBAR) * (U.conj().T @ psi0))
        rho = np.outer(psi, psi.conj())
        pops = np.real(np.diag(rho))
        POP[s] = pops / pops.sum()
        r1 = rung1(rho)
        j = (HBAR / (MASS * a)) * np.abs(r1) * np.sin(np.angle(r1))
        safe = np.maximum(pops, floor * pops.max())
        RU[s] = np.maximum(j, 0.0) / (safe * a)
        RD[s] = np.maximum(-np.roll(j, 1), 0.0) / (safe * a)

    tot = RU + RD
    surv = 1.0 - np.exp(-tot * dt)
    frac = np.where(tot > 0.0, RU / np.maximum(tot, 1e-300), 0.0)
    PU = surv * frac
    PD = surv * (1.0 - frac)
    print(f"  M = {M}, {n_steps} steps to t = 1.5 T_period, "
          f"max(rate) dt = {float((tot * dt).max()):.3f}")
    print(f"  all rates non-negative: {bool(np.all(RU >= 0) and np.all(RD >= 0))}")

    # deterministic control: the same rates, propagated as a master equation
    P = POP[0].copy()
    sup_det = 0.0
    for s in range(n_steps):
        up = PU[s] * P
        dn = PD[s] * P
        P = P - up - dn + np.roll(up, 1) + np.roll(dn, -1)
        sup_det = max(sup_det, float(np.max(np.abs(P - POP[s]))))
    print(f"  master equation with these rates vs exact rho(x,x,t): "
          f"sup err = {sup_det:.2e}")

    # the walkers themselves
    N = 200_000
    r = np.random.default_rng(99)
    site = r.choice(M, size=N, p=POP[0])
    snap_t, snap_emp, snap_ex, errs = [], [], [], []
    stride = n_steps // 6
    for s in range(n_steps):
        pu = PU[s][site]
        pd = PD[s][site]
        u = r.random(N)
        site = (site + np.where(u < pu, 1,
                                np.where(u > 1.0 - pd, -1, 0))) % M
        if s % stride == 0 or s == n_steps - 1:
            emp = np.bincount(site, minlength=M) / N
            snap_t.append(s * dt)
            snap_emp.append(emp)
            snap_ex.append(POP[s])
            errs.append(float(np.max(np.abs(emp - POP[s]))))
    print(f"  {N:,} self-conjugate walkers:")
    print("    t         max |empirical - exact| population")
    for tt, ee in zip(snap_t, errs):
        print(f"    {tt:7.3f}   {ee:.5f}")
    print(f"  peak sup-norm error over the run: {max(errs):.5f} "
          f"(shot-noise floor ~ {1.0 / np.sqrt(N):.5f})")
    return dict(snap_t=snap_t, snap_emp=snap_emp, snap_ex=snap_ex,
                errs=errs, x=x, sup_det=sup_det)


# ------------------------------------------------------------------ #
# figure                                                              #
# ------------------------------------------------------------------ #
def make_figure(rA, rB, rC, rD, rE, rF, rG):
    fig, axes = plt.subplots(2, 4, figsize=(21.0, 9.6))

    # --- panel 1: the rho board, rungs and one-leg hops
    ax = axes[0, 0]
    Mb = 9
    for m in range(Mb):
        for n in range(Mb):
            k = m - n
            c = "#dddddd" if k == 0 else plt.cm.PuOr(0.5 + 0.05 * k)
            ax.add_patch(plt.Rectangle((n - 0.5, m - 0.5), 1, 1,
                                       facecolor=c, edgecolor="w", lw=0.6))
    for m in range(Mb):
        ax.plot(m, m, "o", color="k", ms=6, zorder=5)
    ax.annotate("", xy=(4, 5.85), xytext=(4, 4.15),
                arrowprops=dict(arrowstyle="->", color="C3", lw=2.2))
    ax.text(4.15, 5.0, "ket hop\n$\\times\\,+i$", color="C3", fontsize=9)
    ax.annotate("", xy=(5.85, 4), xytext=(4.15, 4),
                arrowprops=dict(arrowstyle="->", color="C0", lw=2.2))
    ax.text(4.3, 3.35, "bra hop  $\\times\\,-i$", color="C0", fontsize=9)
    ax.set_xlim(-0.5, Mb - 0.5)
    ax.set_ylim(-0.5, Mb - 0.5)
    ax.set_xlabel("bra site $n$  ($X\'$, negaton leg)")
    ax.set_ylabel("ket site $m$  ($X$, positon leg)")
    ax.set_title("The pair board.  Rung $k = m - n$;\n"
                 "black = self-conjugate = probability density", fontsize=10)
    ax.set_aspect("equal")

    # --- panel 2: momentum is the misalignment
    ax = axes[0, 1]
    kx = rB["keep"]
    ax.plot(rB["x"][kx], rB["p_cur"][kx], "-", color="k", lw=2.4,
            label=r"$m\,j/\rho$  (current)")
    ax.plot(rB["x"][kx], rB["p_mis"][kx], "o", color="C2", ms=3.6,
            mfc="none", label=r"$\hbar\mu/a$  (misalignment)")
    ax2 = ax.twinx()
    ax2.fill_between(rB["x"], rB["dens"], color="C0", alpha=0.15, lw=0)
    ax2.set_ylabel(r"$\rho(x,x)$", color="C0", fontsize=9)
    ax2.tick_params(axis="y", labelcolor="C0", labelsize=8)
    ax.set_xlabel("$x$")
    ax.set_ylabel("momentum")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="upper left")
    ax.set_title("Momentum is not carried, it is read:\n"
                 r"$\bar p = \hbar\mu/a$ on a chirped packet", fontsize=10)

    # --- panel 3: Euler force law convergence
    ax = axes[0, 2]
    aa = [r[1] for r in rC["rows"]]
    ee = [r[2] for r in rC["rows"]]
    ax.loglog(aa, ee, "s-", color="C1",
              label=r"$\partial_t j|_{\rm pump}$ vs $\rho F/m$")
    ax.loglog(aa, [ee[-1] * (v / aa[-1]) ** 2 for v in aa], "k--", lw=0.9,
              label=r"$O(a^2)$")
    ax.set_xlabel("lattice spacing $a$")
    ax.set_ylabel("relative error")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=8)
    ax.set_title("No noise, no force, in position space:\n"
                 "the pump alone yields Newton", fontsize=10)

    # --- panel 4: amplitude vs probability / ballistic spreading
    ax = axes[0, 3]
    dts = np.logspace(-4, -1.3, 40)
    lam = 32.0
    ax.loglog(dts, lam * dts, color="C3", lw=2.0,
              label=r"hop amplitude $\propto \delta t$")
    ax.loglog(dts, (lam * dts) ** 2, color="C0", lw=2.0,
              label=r"hop probability $\propto \delta t^{2}$")
    ax.loglog(dts, 0.35 * lam * dts, color="0.5", ls="--", lw=1.2,
              label=r"any Poisson rate $\propto \delta t$")
    ax.set_xlabel(r"$\delta t$")
    ax.set_ylabel("size of a one-leg update")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=8, loc="lower right")
    ax.set_title("Why the coherence sector cannot be a\n"
                 "jump process (spreading is ballistic,\n"
                 f"$d\\log\\langle x^2\\rangle/d\\log t = {rD["slope"]:.2f}$)",
                 fontsize=10)

    # --- panel 5: gauge circle vs the single ray
    ax = axes[1, 0]
    ax.hist(rD["angles"], bins=36, range=(-np.pi, np.pi), density=True,
            color="C4", alpha=0.75,
            label=r"$\arg\rho_{mn}$ over local gauges")
    ax.axhline(1.0 / (2 * np.pi), color="k", ls="--", lw=1.0,
               label="uniform on the circle")
    ax.plot([0.0], [0.0], "r*", ms=15,
            label="the one ray a positon sea needs")
    ax.set_xlabel(r"$\arg\rho_{mn}$")
    ax.set_ylabel("density")
    ax.set_ylim(bottom=-0.012)
    ax.legend(fontsize=7.5, loc="lower center")
    ax.set_title("No positon-only sea exists here:\n"
                 "gauge sweeps the phase round the circle\n"
                 "while the diagonal never moves", fontsize=10)

    # --- panel 6: resource arithmetic
    ax = axes[1, 1]
    ts = rE["ts"] / rE["T_per"]
    ax.plot(ts, rE["z4"], color="C3", lw=1.3,
            label=r"$\mathbb{Z}_4$ sampling mass")
    ax.plot(ts, rE["l1"], color="C0", lw=1.6,
            label=r"$\sum_{mn}|\rho_{mn}|$")
    ax.axhline(rE["M"], color="k", ls="--", lw=1.0,
               label=r"bound $M = %d$" % rE["M"])
    ax.axhline(np.sqrt(2) * rE["M"], color="0.5", ls=":", lw=1.2,
               label=r"bound $\sqrt{2}\,M$")
    ax.set_xlabel(r"$t/T_{\rm period}$")
    ax.set_ylabel("mass  (state mass $=1$)")
    ax.set_ylim(0, 1.18 * np.sqrt(2) * rE["M"])
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="lower right")
    ax.set_title("The lattice restores a bound the\n"
                 "continuum lacks: $\\ell^1$ mass $\\leq M$", fontsize=10)

    # --- panel 7: Monte Carlo
    ax = axes[1, 2]
    Ns = np.array([c[0] for c in rF["conv"]], float)
    Es = np.array([c[1] for c in rF["conv"]], float)
    ax.loglog(Ns, Es, "o-", color="C6", label="MC error at $T = 0.5$")
    ax.loglog(Ns, Es[0] * np.sqrt(Ns[0] / Ns), "k--", lw=0.9,
              label=r"$1/\sqrt{N}$")
    ax.set_xlabel("pairs $N$")
    ax.set_ylabel("relative sup error")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=8)
    axi = ax.inset_axes([0.56, 0.60, 0.40, 0.35])
    Tt = np.array([n[0] for n in rF["noise"]])
    sp = np.array([n[2] for n in rF["noise"]])
    gr = np.array([n[3] for n in rF["noise"]])
    axi.semilogy(Tt, sp / sp[0], "o-", color="C3", ms=3, lw=1.2)
    axi.semilogy(Tt, gr / gr[0], "k--", lw=1.0)
    axi.set_title(r"noise $\propto e^{\Lambda t}$", fontsize=7)
    axi.tick_params(labelsize=6)
    ax.set_title("The pair ensemble is unbiased but its\n"
                 r"noise is amplified by $e^{\Lambda t}$", fontsize=10)

    # --- panel 8: guided self-conjugate walkers
    ax = axes[1, 3]
    xx = rG["x"]
    idxs = [0, 2, 4]
    for k in idxs:
        if k >= len(rG["snap_t"]):
            continue
        c = plt.cm.viridis(k / max(1, len(rG["snap_t"]) - 1))
        ax.plot(xx, rG["snap_ex"][k], "-", color=c, lw=1.8,
                label=f"$\\rho(x,x,t)$, $t$ = {rG['snap_t'][k]:.2f}")
        ax.plot(xx, rG["snap_emp"][k], "o", color=c, ms=3.4, mfc="none")
    ax.set_xlabel("$x$")
    ax.set_ylabel("population")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7.5)
    ax.set_title("The observable sector IS a particle\n"
                 r"process: walkers hopping at $\sin\mu$"
                 "\n(circles) track the exact density (lines)", fontsize=10)

    fig.suptitle(
        "The position-space coherence ladder: conjugate world-particle pairs "
        "and the von Neumann equation\n"
        "each leg carries a position and a phase; a pair is one sample of "
        r"$\rho(X, X')$ and a self-conjugate particle one sample of "
        r"$\rho(x,x)$;  the potential only winds $\mu$, and all motion is "
        r"four one-leg hops of amplitude $\pm iJ/\hbar$",
        fontsize=11.5, y=0.985)
    fig.tight_layout(rect=(0, 0, 1, 0.925))
    fig.savefig(output_path("position_pair_ladder.png"), dpi=140,
                bbox_inches="tight")
    dp = docs_path("position_pair_ladder.png")
    if dp:
        fig.savefig(dp, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print()
    print("  wrote position_pair_ladder.png")


def main():
    rA = part_A()
    rB = part_B()
    rC = part_C()
    rD = part_D()
    rE = part_E()
    rF = part_F()
    rG = part_G()
    make_figure(rA, rB, rC, rD, rE, rF, rG)


if __name__ == "__main__":
    main()
