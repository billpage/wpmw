#!/usr/bin/env python3
"""
Three (and four) particles in many worlds: does Proposition WC8 suffice?

Companion demo for step 25 (``docs/analysis/world_clocks.md`` §8.2, addendum).  Three
particles of unit mass on a line, hbar = 1, in a harmonic trap of frequency 1,
interacting pairwise.  The world ontology of step 25 (Definition WC7): each
world carries its point (x1, x2, x3, p1, p2, p3), a sign and one bank of
integrate-and-fire integrators per bond; the generator is read from world
quantities only (Proposition WC8).  Nothing is attached to a world-particle.

Two problems, both quenches from the product of trap ground states:

  Q  harmonic pair coupling v(r) = k r^2 / 2 (a three-particle Moshinsky
     model).  Exactly solvable; quadratic, so the compensated kernel vanishes
     (Theorem G2) and the worlds make no events at all.
  G  Gaussian pair interaction v(r) = g exp(-r^2/2s^2) / (sqrt(2 pi) s),
     the regularised contact interaction of Zoellner, Meyer and Schmelcher's
     few-boson studies; no closed form.  Reference: the time-dependent
     Schroedinger equation on a 3D grid in particle coordinates.

Entanglement is measured by the purity of particle 1's reduced state,
P1 = Tr rho_1^2 = 2 pi hbar \\int W_1^2, and S2 = -ln P1.

Parts (default: S K E Q M W X F; M at reach 8 pi takes about half an hour)
---------------------------------------------------------------------------
S  Structure (numerical identities, Proposition WC9): the residual symbol of
   the whole potential at a random chord equals the sum over bonds of the 1D
   pair symbols at eta_b = y_i - y_j, plus nothing from the quadratic trap; a
   bond event moves only p_i and p_j and conserves P; the Jacobi bond vectors.
K  The bond kernel: channels and rates at reach 2 pi, 4 pi, 8 pi, and the
   largest rate as the Gaussian narrows toward a contact interaction.
E  Exact references for G: the TDSE on two grids (convergence), purity by SVD,
   pair energy.  Also the classical control: the initial (positive) Wigner
   function carried by classical streaming alone.
Q  Problem Q: covariance flow (exact), TDSE cross-check, and the worlds
   (Verlet streaming; the kernel's largest rate is printed: zero to rounding).
M  The compensated QLE on the 4D Jacobi relative phase space with the bond-sum
   symbol, at reach 2 pi, 4 pi, 8 pi (``--mesh-reach``): the generator of WC8
   without sampling, against E.  Also the global purity at resolution GA.
W  The worlds for G: integrate-and-fire banks per bond fed by the bond kernel,
   emission realisation, contact sink on cells of the Jacobi relative phase
   space; ``--reach``, ``--cell``, ``--nu``, ``--seeds``, ``--T``, and the
   diagnostics ``--dt``, ``--sink-every``, ``--levels`` (hierarchical sink).
   Against M at the same reach and against E.  ``--N 4`` for four particles.
X  Every saved world run against the mesh at its reach: purity, pair energy,
   the share q of the quantum correction kept, the global purity
   (Proposition WC10).
F  The figure, from the JSON the other parts save.

Run:  WPMW_OUTPUT=... PYTHONPATH=src python3 -u src/demo_three_worlds.py [S K E Q M W F] [options]
"""
import argparse
import glob
import json
import os
import time

import numpy as np

from wpmwlib.wpmw_utils import docs_path, output_path

HB, OMEGA = 1.0, 1.0
G_INT, S_INT = 2.0, 0.7            # problem G
K_MOSH = 2.0                        # problem Q: v(r) = k r^2 / 2
T_END, DT = 3.0, 0.01
OBS = 0.25                          # observation interval


def jacobi(n):
    """Orthonormal Jacobi coordinates for n equal masses: rows 0 .. n-2 relative,
    row n-1 the centre of mass R = sum x / sqrt(n)."""
    Oj = np.zeros((n, n))
    for k in range(1, n):
        Oj[k - 1, :k] = 1.0
        Oj[k - 1, k] = -k
        Oj[k - 1] /= np.sqrt(k * (k + 1))
    Oj[n - 1] = 1 / np.sqrt(n)
    return Oj


def set_particles(n):
    """Particle count n: bonds, Jacobi matrix, bond vectors u_b (r_b = u_b . rho),
    particle 1's relative projection a (x1 = a . rho + R / sqrt n) and the
    centre-of-mass variance of R / sqrt n in the trap ground state."""
    global NP, BONDS, OJ, UB, A1, VAR_CM
    NP = n
    BONDS = tuple((i, j) for i in range(n) for j in range(i + 1, n))
    OJ = jacobi(n)
    if n == 3:          # the convention of the 3-particle notes: rho1 = (x1 - x2)/sqrt2, ...
        OJ = np.array([[1 / np.sqrt(2), -1 / np.sqrt(2), 0.0],
                       [1 / np.sqrt(6), 1 / np.sqrt(6), -2 / np.sqrt(6)],
                       [1 / np.sqrt(3), 1 / np.sqrt(3), 1 / np.sqrt(3)]])
    UB = np.array([OJ[:n - 1] @ (np.eye(n)[i] - np.eye(n)[j]) for i, j in BONDS])
    A1 = OJ[:n - 1, 0]
    VAR_CM = 0.5 / n


set_particles(3)


def rule(text):
    print()
    print(text)
    print("-" * len(text))


def dump(obj, name):
    with open(output_path(name), "w") as f:
        json.dump(obj, f, indent=1)


def load(name):
    p = output_path(name)
    return json.load(open(p)) if os.path.exists(p) else None


# ----------------------------------------------------------------------------
# pair potentials
# ----------------------------------------------------------------------------
def v_gauss(r):
    return G_INT / (np.sqrt(2 * np.pi) * S_INT) * np.exp(-r**2 / (2 * S_INT**2))


def dv_gauss(r):
    return -r / S_INT**2 * v_gauss(r)


def v_harm(r):
    return 0.5 * K_MOSH * r**2


def dv_harm(r):
    return K_MOSH * r


POT = {"G": (v_gauss, dv_gauss), "Q": (v_harm, dv_harm)}
TDSE_GRID = {3: (64, 7.0), 4: (40, 6.0)}
TDSE_FINE = {3: (96, 8.0), 4: (48, 6.5)}


def TAG():
    return "" if NP == 3 else f"_N{NP}"


def bond_symbol(r, sig, v, dv, ymax):
    """Compensated residual symbol of one bond, (i/hbar)[v(r+y) - v(r-y) - 2 y v'(r)]
    times the raised-cosine horizon at the reach ymax, with y = hbar sigma / 2."""
    y = HB * sig / 2
    w = np.where(np.abs(y) < ymax, np.cos(np.pi * y / (2 * ymax))**2, 0.0)
    return 1j / HB * (v(r + y) - v(r - y) - 2 * y * dv(r)) * w


def kernel_table(v, dv, mult, rh=12.0, nr=6001):
    """K_q(r) for one bond, q = 1 .. NQ: transfer xi_q = q dp to p_i and -xi_q to
    p_j (relative momentum (p_i - p_j)/2 moves by xi_q), dp = pi hbar / 2 ymax."""
    ymax = mult * np.pi
    dp = np.pi * HB / (2 * ymax)
    nrow = int(round(16 / dp))
    rt = np.linspace(-rh, rh, nr)
    s = 2 * np.pi * np.fft.fftfreq(nrow, d=dp)
    sym = bond_symbol(rt[:, None], s[None, :], v, dv, ymax)
    sym[:, nrow // 2] = 0
    k = np.real(np.fft.ifft(sym, axis=1))[:, 1:nrow // 2]
    return rt, k, dp * np.arange(1, nrow // 2)


# ----------------------------------------------------------------------------
# Part S -- structure
# ----------------------------------------------------------------------------
def part_S():
    rule("Part S. Structure: the world's symbol is a sum over bonds")
    rng = np.random.default_rng(0)
    x = rng.normal(0, 1.5, (2000, NP))
    y = rng.normal(0, 1.0, (2000, NP))
    out = {}
    for name, (v, dv) in POT.items():
        V = lambda q: 0.5 * OMEGA**2 * (q**2).sum(1) + sum(v(q[:, i] - q[:, j]) for i, j in BONDS)
        dV = lambda q: OMEGA**2 * q + np.stack([sum((dv(q[:, i] - q[:, j]) if a == i else
                                                     -dv(q[:, i] - q[:, j]) if a == j else 0.0)
                                                    for i, j in BONDS) for a in range(NP)], 1)
        full = V(x + y) - V(x - y) - 2 * (y * dV(x)).sum(1)
        bonds = sum(v(x[:, i] - x[:, j] + (y[:, i] - y[:, j])) - v(x[:, i] - x[:, j] - (y[:, i] - y[:, j]))
                    - 2 * (y[:, i] - y[:, j]) * dv(x[:, i] - x[:, j]) for i, j in BONDS)
        err = float(np.abs(full - bonds).max())
        out[name] = dict(maxdiff=err, maxres=float(np.abs(full).max()))
        print(f"  {name}: residual of the whole potential minus the sum over bonds, max {err:.1e}"
              f"   (largest residual {np.abs(full).max():.2e})")
    p = rng.normal(0, 1, (5, NP))
    for b, (i, j) in enumerate(BONDS):
        q = p.copy()
        q[:, i] += 0.75
        q[:, j] -= 0.75
        dpi = (OJ @ (q - p).T).T
        print(f"  bond {i+1}{j+1}: dP = {np.abs((q - p).sum(1)).max():.1e}, d(rho-momenta) = "
              f"{dpi[0, :NP - 1].round(4)} = 0.75 u_b {(0.75 * UB[b]).round(4)}, dR-momentum {dpi[0, NP - 1]:.1e}")
    print(f"  |u_b|^2 = {[round(float(u @ u), 12) for u in UB]};  x1 = a.rho + R/sqrt3, a = {A1.round(4)}")
    return out


# ----------------------------------------------------------------------------
# Part K -- the bond kernel
# ----------------------------------------------------------------------------
def part_K():
    rule("Part K. The bond kernel of problem G")
    rng = np.random.default_rng(1)
    x = rng.normal(0, np.sqrt(0.5), (200000, NP))
    out = {}
    for mult in (2, 4, 8):
        rt, k, xi = kernel_table(v_gauss, dv_gauss, mult)
        gam = np.abs(k).sum(1)
        mean = float(sum(np.interp(x[:, i] - x[:, j], rt, gam).mean() for i, j in BONDS))
        out[str(mult)] = dict(nq=int(k.shape[1]), dp=float(xi[0]), gmax=float(gam.max()), gmean=mean)
        print(f"  reach {mult} pi: {k.shape[1]} channels, dp = {xi[0]:.4f}, max rate per bond "
              f"{gam.max():.3f}, mean over the initial state, summed over bonds {mean:.3f}")
    _, kq, _ = kernel_table(v_harm, dv_harm, 4)
    print(f"  problem Q (quadratic): largest rate {np.abs(kq).sum(1).max():.1e} (Theorem G2)")
    out["Q"] = float(np.abs(kq).sum(1).max())
    # towards a contact interaction: the same g, narrower Gaussians (reach 4 pi)
    out["width"] = {}
    for w in (0.7, 0.35, 0.2, 0.1):
        vw = lambda r, w=w: G_INT / (np.sqrt(2 * np.pi) * w) * np.exp(-r**2 / (2 * w**2))
        dvw = lambda r, w=w: -r / w**2 * vw(r)
        rt, k, _ = kernel_table(vw, dvw, 4, nr=24001)
        gam = np.abs(k).sum(1)
        mean = float(sum(np.interp(x[:, i] - x[:, j], rt, gam).mean() for i, j in BONDS))
        out["width"][str(w)] = dict(gmax=float(gam.max()), gmean=mean)
        print(f"  width s = {w}: max rate per bond {gam.max():8.2f}, mean over the initial state, "
              f"summed over bonds {mean:8.2f}")
    return out


# ----------------------------------------------------------------------------
# Part E -- exact references (TDSE on a 3D grid)
# ----------------------------------------------------------------------------
def tdse(v, n=64, L=7.0, T=T_END, dt=0.005, xgrid=None):
    x = (np.arange(n) - n // 2) * (2 * L / n)
    dx = x[1] - x[0]
    k = 2 * np.pi * np.fft.fftfreq(n, d=dx)
    shp = lambda a, d: a.reshape([-1 if e == d else 1 for e in range(NP)])
    Vp = sum(v(shp(x, i) - shp(x, j)) for i, j in BONDS)
    V = 0.5 * OMEGA**2 * sum(shp(x, d)**2 for d in range(NP)) + Vp
    Tk = 0.5 * sum(shp(k, d)**2 for d in range(NP))
    phi = (OMEGA / np.pi)**0.25 * np.exp(-OMEGA * x**2 / 2)
    psi = np.ones([1] * NP, complex)
    for d in range(NP):
        psi = psi * shp(phi, d)
    psi /= np.sqrt(np.sum(abs(psi)**2))
    eV = np.exp(-0.5j * V * dt)
    eT = np.exp(-1j * Tk * dt)
    every = int(round(OBS / dt))
    rows, dens = [], {}
    nst = int(round(T / dt))
    for i in range(nst + 1):
        if i % every == 0:
            rho = abs(psi)**2
            M = psi.reshape(n, n**(NP - 1))
            r1 = M @ M.conj().T                      # reduced density matrix of particle 1
            pur = float(np.sum(np.abs(r1)**2))
            ph = np.fft.fftn(psi)
            Tkin = float(np.sum(abs(ph)**2 * Tk) / np.sum(abs(ph)**2))
            rows.append(dict(t=round(i * dt, 6), purity=pur, Vpair=float(np.sum(rho * Vp)),
                             T=Tkin, Vtrap=float(np.sum(rho * (V - Vp)))))
            dens[round(i * dt, 6)] = (rho.sum(axis=tuple(range(1, NP))) / dx).tolist()
        if i < nst:
            psi = eV * np.fft.ifftn(eT * np.fft.fftn(eV * psi))
    return dict(x=x.tolist(), rows=rows, rho1=dens)


def classical(v, dv, nu=400000, T=T_END, dt=DT, seed=11):
    """The initial Wigner function (positive) carried by classical streaming alone."""
    rng = np.random.default_rng(seed)
    x = rng.normal(0, np.sqrt(0.5), (nu, NP))
    p = rng.normal(0, np.sqrt(0.5), (nu, NP))
    eps = np.ones(nu)
    rows = []
    every = int(round(OBS / dt))
    f = force(x, dv)
    for i in range(int(round(T / dt)) + 1):
        if i % every == 0:
            rows.append(dict(t=round(i * dt, 6), **observe(x, p, eps, nu, v)))
        p += 0.5 * dt * f
        x += p * dt
        f = force(x, dv)
        p += 0.5 * dt * f
    return rows


def part_E():
    rule(f"Part E. Exact references for problem G, {NP} particles (TDSE on a {NP}D grid)")
    t0 = time.perf_counter()
    a = tdse(v_gauss, *TDSE_GRID[NP])
    b = tdse(v_gauss, *TDSE_FINE[NP])
    dp = max(abs(r1["purity"] - r2["purity"]) for r1, r2 in zip(a["rows"], b["rows"]))
    dv_ = max(abs(r1["Vpair"] - r2["Vpair"]) for r1, r2 in zip(a["rows"], b["rows"]))
    print(f"  {TDSE_GRID[NP]} against {TDSE_FINE[NP]} (n, half-width): max |d purity| {dp:.1e}, "
          f"max |d <Vpair>| {dv_:.1e}"
          f"  ({time.perf_counter() - t0:.0f} s)")
    print("     t    purity   S2=-ln P   <Vpair>     <T>     <Vtrap>     E")
    for r in a["rows"]:
        print(f"  {r['t']:4.2f}   {r['purity']:.5f}   {-np.log(r['purity']):.5f}   {r['Vpair']:.5f}  "
              f"{r['T']:.5f}   {r['Vtrap']:.5f}   {r['T'] + r['Vtrap'] + r['Vpair']:.5f}")
    cl = classical(v_gauss, dv_gauss)
    print("  classical streaming of the initial W (400 000 samples):")
    for r in cl:
        if abs(r["t"] - round(r["t"])) < 1e-9:
            print(f"    t = {r['t']:3.1f}: purity {r['purity']:.4f}, <Vpair> {r['Vpair']:.4f}")
    a["conv"] = dict(dpurity=dp, dVpair=dv_)
    a["classical"] = cl
    return a


# ----------------------------------------------------------------------------
# observables of a signed world ensemble
# ----------------------------------------------------------------------------
def force(x, dv):
    f = -OMEGA**2 * x
    for i, j in BONDS:
        d = dv(x[:, i] - x[:, j])
        f[:, i] -= d
        f[:, j] += d
    return f


PGRID = 0.05
PEXT = 6.0


def deposit_cic(z, wgt, nb):
    """Cloud-in-cell deposition on the purity grid (for mesh points, whose regular
    projections would alias against nearest-bin deposition)."""
    u = (z + PEXT) / PGRID - 0.5
    i0 = np.floor(u).astype(np.int64)
    f = u - i0
    h = np.zeros((nb, nb))
    for di in (0, 1):
        for dj in (0, 1):
            w = wgt * (f[:, 0] if di else 1 - f[:, 0]) * (f[:, 1] if dj else 1 - f[:, 1])
            ii, jj = i0[:, 0] + di, i0[:, 1] + dj
            ok = (ii >= 0) & (ii < nb) & (jj >= 0) & (jj < nb)
            np.add.at(h, (ii[ok], jj[ok]), w[ok])
    return h


def purity_from_relative(z, wgt, cic=False):
    """Purity of particle 1 from relative coordinates: W1 = (projection of the
    relative ensemble on x1 = a.rho, p1 = a.pi) convolved with the exact
    centre-of-mass Gaussian (variance 1/6 in x and p).  z: (n, 2) projected
    points; wgt: signed weights summing to the norm.  Self-pairs are removed, so
    the estimator is unbiased for a sample; for a mesh pass the self term 0."""
    nb = int(round(2 * PEXT / PGRID))
    if cic:
        h = deposit_cic(z, wgt, nb)
    else:
        h, _, _ = np.histogram2d(z[:, 0], z[:, 1], bins=nb, range=[[-PEXT, PEXT], [-PEXT, PEXT]], weights=wgt)
    k = 2 * np.pi * np.fft.fftfreq(nb, d=PGRID)
    g = np.exp(-0.5 * VAR_CM * (k[:, None]**2 + k[None, :]**2))
    hk = np.fft.fft2(h)
    if cic:     # undo the cloud-in-cell smoothing (its transfer function is sinc^2 per axis)
        sc = np.sinc(k * PGRID / (2 * np.pi))**2
        hk = hk / (sc[:, None] * sc[None, :])
    # \int W1^2 = sum_k |h_k|^2 |g_k|^2 / (nb^2 * cell area) in the discrete Parseval form
    total = float(np.sum(np.abs(hk)**2 * g**2)) / nb**2 / PGRID**2
    return 2 * np.pi * HB * total


GA, GH, GEXT = 0.25, 0.25, 6.0     # smoothing width, grid spacing and extent of the global-purity grid


def _gauss_filter(shape, spacings):
    ks = [2 * np.pi * np.fft.fftfreq(n, d=h) for n, h in zip(shape, spacings)]
    k2 = sum(np.meshgrid(*[k**2 for k in ks], indexing="ij"))
    return np.exp(-GA**2 * k2), ks


def global_purity_mesh(W, spacings):
    """Purity of the whole relative state seen at resolution GA:
    P_a = (2 pi hbar)^(d/2) * integral of (W * G_a)^2, G_a a Gaussian of width GA in every
    Jacobi coordinate.  A unitary evolution keeps P_0 = 1; P_a falls as W develops
    structure finer than GA.  Exact (spectral) on the mesh."""
    g, _ = _gauss_filter(W.shape, spacings)
    vol = float(np.prod(spacings))
    return float((2 * np.pi * HB)**(W.ndim // 2) * vol / W.size * np.sum(np.abs(np.fft.fftn(W))**2 * g))


def global_purity_points(z, e):
    """The same P_a for a signed ensemble: cloud-in-cell deposition on a grid of side
    GH, CIC transfer function divided out, self-pairs removed (unbiased for a sample)."""
    d = z.shape[1]
    nb = int(round(2 * GEXT / GH))
    u = (z + GEXT) / GH - 0.5
    i0 = np.floor(u).astype(np.int64)
    f = u - i0
    H = np.zeros(nb**d)
    for corner in range(2**d):
        idx = np.zeros(e.size, np.int64)
        w = e.copy()
        ok = np.ones(e.size, bool)
        for ax in range(d):
            bit = (corner >> ax) & 1
            ii = i0[:, ax] + bit
            ok &= (ii >= 0) & (ii < nb)
            idx = idx * nb + np.clip(ii, 0, nb - 1)
            w = w * (f[:, ax] if bit else 1 - f[:, ax])
        H += np.bincount(idx[ok], weights=w[ok], minlength=nb**d)
    H = H.reshape((nb,) * d) / GH**d
    g, ks = _gauss_filter(H.shape, (GH,) * d)
    T = np.ones(H.shape)
    for ax, k in enumerate(ks):
        sh = [1] * d
        sh[ax] = nb
        T = T * (np.sinc(k * GH / (2 * np.pi))**2).reshape(sh)
    total = GH**d / H.size * np.sum(np.abs(np.fft.fftn(H))**2 / T**2 * g)
    self_term = np.sum(e**2) / (4 * np.pi * GA**2)**(d / 2)
    return float((2 * np.pi * HB)**(d // 2) * (total - self_term))


def observe(x, p, eps, nu, v):
    e = eps / nu
    rho = x @ OJ[:NP - 1].T
    pr = p @ OJ[:NP - 1].T
    z = np.stack([rho @ A1, pr @ A1], 1)
    self_term = np.sum(e**2) / (4 * np.pi * VAR_CM) * 2 * np.pi * HB     # 2 pi \int G^2 per unit weight^2
    pur = purity_from_relative(z, e) - self_term
    Vp = sum(v(x[:, i] - x[:, j]) for i, j in BONDS)
    gp = global_purity_points(np.concatenate([rho, pr], 1), e) if NP == 3 else float("nan")
    return dict(norm=float(e.sum()), purity=float(pur), gpurity=gp, Vpair=float((e * Vp).sum()),
                T=float((e * 0.5 * (p**2).sum(1)).sum()), Vtrap=float((e * 0.5 * (x**2).sum(1)).sum()),
                Pcm2=float((e * p.sum(1)**2).sum() / NP), N=int(eps.size))


# ----------------------------------------------------------------------------
# Part Q -- the three-particle Moshinsky quench
# ----------------------------------------------------------------------------
def moshinsky_purity(n, times):
    """Exact purity of particle 1 after the quench of problem Q: the covariance
    matrix of the Gaussian Wigner function flows by the linear equations of
    motion, and Tr rho_1^2 = hbar / (2 sqrt(det S_1)) for the 2x2 block of particle 1."""
    from scipy.linalg import expm
    Kmat = OMEGA**2 * np.eye(n) + K_MOSH * (n * np.eye(n) - np.ones((n, n)))
    A = np.block([[np.zeros((n, n)), np.eye(n)], [-Kmat, np.zeros((n, n))]])
    S0 = 0.5 * np.eye(2 * n)
    out = []
    for t in times:
        M = expm(A * t)
        S1 = (M @ S0 @ M.T)[np.ix_([0, n], [0, n])]
        out.append(float(HB / (2 * np.sqrt(np.linalg.det(S1)))))
    return out


def part_Q(nu=1000000, seed=5):
    rule(f"Part Q. Problem Q, {NP} particles: harmonic pair coupling k = {K_MOSH}, quench from the trap ground state")
    n = NP
    times = [round(OBS * i, 6) for i in range(int(round(T_END / OBS)) + 1)]
    exact = moshinsky_purity(n, times)
    td = tdse(v_harm, *TDSE_GRID[NP])
    # worlds: positive, never fire (the kernel is zero); Verlet streaming
    rng = np.random.default_rng(seed)
    x = rng.normal(0, np.sqrt(0.5), (nu, n))
    p = rng.normal(0, np.sqrt(0.5), (nu, n))
    eps = np.ones(nu)
    _, kq, _ = kernel_table(v_harm, dv_harm, 4)
    every = int(round(OBS / DT))
    wr = []
    f = force(x, dv_harm)
    for i in range(int(round(T_END / DT)) + 1):
        if i % every == 0:
            o = observe(x, p, eps, nu, v_harm)
            c = np.cov(np.stack([x[:, 0], p[:, 0]]))
            o["purity_gauss"] = float(1 / (2 * np.sqrt(np.linalg.det(c))))
            wr.append(o)
        p += 0.5 * DT * f
        x += p * DT
        f = force(x, dv_harm)
        p += 0.5 * DT * f
    print(f"  largest kernel rate: {np.abs(kq).sum(1).max():.1e}, so no event fires; {nu} worlds")
    print("     t    exact (covariance)   TDSE 64^3    worlds: Gaussian fit   worlds: estimator")
    for t, e_, r, w in zip(times, exact, td["rows"], wr):
        print(f"  {t:4.2f}      {e_:.5f}          {r['purity']:.5f}        {w['purity_gauss']:.5f}"
              f"            {w['purity']:.5f}")
    err = max(abs(w["purity_gauss"] - e_) for e_, w in zip(exact, wr))
    errt = max(abs(r["purity"] - e_) for e_, r in zip(exact, td["rows"]))
    print(f"  max |worlds - exact| {err:.1e} (Gaussian fit); max |TDSE - exact| {errt:.1e}; "
          f"min purity {min(exact):.4f} (S2 = {-np.log(min(exact)):.4f})")
    return dict(times=times, exact=exact, tdse=[r["purity"] for r in td["rows"]],
                worlds=[w["purity_gauss"] for w in wr], worlds_est=[w["purity"] for w in wr])


# ----------------------------------------------------------------------------
# Part M -- the compensated QLE on the 4D Jacobi relative phase space
# ----------------------------------------------------------------------------
def mesh_run(mult, nr=40, Lr=5.0, Lp=5.0, T=T_END, dt=DT):
    ymax = mult * np.pi
    dpi = np.pi / (np.sqrt(2) * ymax) * 0.98       # the bond chords |u.s|/2 <= ymax fit the grid
    npi = int(np.ceil(Lp / dpi)) * 2
    dpi = 2 * Lp / npi
    rho = -Lr + 2 * Lr * np.arange(nr) / nr
    drho = rho[1] - rho[0]
    pi_ = -Lp + dpi * np.arange(npi)
    R1, R2 = np.meshgrid(rho, rho, indexing="ij")
    s = 2 * np.pi * np.fft.fftfreq(npi, d=dpi)
    S1, S2 = np.meshgrid(s, s, indexing="ij")
    kr = 2 * np.pi * np.fft.fftfreq(nr, d=drho)
    dV1, dV2 = OMEGA**2 * R1, OMEGA**2 * R2
    expo = np.zeros((nr, nr, npi, npi), complex)
    for u in UB:
        rb = u[0] * R1 + u[1] * R2
        d = dv_gauss(rb)
        dV1 = dV1 + u[0] * d
        dV2 = dV2 + u[1] * d
        expo += bond_symbol(rb[:, :, None, None], (u[0] * S1 + u[1] * S2)[None, None], v_gauss, dv_gauss, ymax)
    expo += 1j * (dV1[:, :, None, None] * S1[None, None] + dV2[:, :, None, None] * S2[None, None])
    kick = np.exp(dt * expo)
    del expo
    K1, K2 = np.meshgrid(kr, kr, indexing="ij")
    P1, P2 = np.meshgrid(pi_, pi_, indexing="ij")
    drift = np.exp(-1j * (K1[:, :, None, None] * P1[None, None] + K2[:, :, None, None] * P2[None, None]) * dt / 2)
    W = (np.exp(-R1**2 - R2**2)[:, :, None, None] * np.exp(-P1**2 - P2**2)[None, None])
    W /= W.sum() * drho**2 * dpi**2
    Vp = sum(v_gauss(u[0] * R1 + u[1] * R2) for u in UB)
    # projection of each 4D cell on (x1, p1) = (a.rho, a.pi), for the purity
    zx = (A1[0] * R1 + A1[1] * R2).ravel()
    zp = (A1[0] * P1 + A1[1] * P2).ravel()
    zx = np.repeat(zx, npi * npi)
    zp = np.tile(zp, nr * nr)
    rows = []
    every = int(round(OBS / dt))
    t0 = time.perf_counter()
    for i in range(int(round(T / dt)) + 1):
        if i % every == 0:
            rr = W.sum(axis=(2, 3)) * dpi**2
            wgt = W.ravel() * drho**2 * dpi**2
            pur = purity_from_relative(np.stack([zx, zp], 1), wgt, cic=True)
            neg = float(-wgt[wgt < 0].sum())
            gp = global_purity_mesh(W, (drho, drho, dpi, dpi))
            rows.append(dict(t=round(i * dt, 6), purity=pur, Vpair=float((rr * Vp).sum() * drho**2),
                             norm=float(wgt.sum()), negmass=neg, gpurity=gp))
            print(f"    reach {mult} pi, t = {i*dt:4.2f}: purity {pur:.5f}, <Vpair> {rows[-1]['Vpair']:.5f},"
                  f" negative mass {neg:.2e}, global purity at {GA} {gp:.4f}"
                  f"  ({time.perf_counter() - t0:.0f} s)", flush=True)
        if i == int(round(T / dt)):
            break
        W = np.real(np.fft.ifft2(np.fft.fft2(W, axes=(0, 1)) * drift, axes=(0, 1)))
        W = np.real(np.fft.ifft2(np.fft.fft2(W, axes=(2, 3)) * kick, axes=(2, 3)))
        W = np.real(np.fft.ifft2(np.fft.fft2(W, axes=(0, 1)) * drift, axes=(0, 1)))
    return dict(mult=mult, nr=nr, npi=npi, drho=drho, dpi=dpi, rows=rows)


def part_M(mults):
    rule("Part M. The bond-sum generator on the 4D Jacobi relative mesh")
    out = {}
    for m in mults:
        r = mesh_run(m)
        out[str(m)] = r
        dump(r, f"three_worlds_M{m}.json")
    return out


# ----------------------------------------------------------------------------
# Part W -- the worlds, one bank per bond
# ----------------------------------------------------------------------------
class Worlds:
    """Signed worlds (x1, x2, x3, p1, p2, p3); one bank of NQ integrators per bond;
    a firing of channel q on bond (i, j) is an emission: two worlds at the
    parent's position with (p_i, p_j) shifted by (+xi_q, -xi_q) and (-xi_q, +xi_q),
    of signs +tau and -tau, with fresh uniform integrator phases.  A contact sink
    then removes opposite-sign worlds sharing a cell of the Jacobi relative phase
    space (rho1, rho2, pi1, pi2), which leaves every cell's net count unchanged."""

    def __init__(self, mult, cell, nu, seed, events=True, sink_every=1, levels=1):
        self.sink_every, self.levels = sink_every, levels
        self.verbose = False
        self.rng = np.random.default_rng(seed)
        self.rt, self.kt, self.xi = kernel_table(v_gauss, dv_gauss, mult)
        self.nq = self.kt.shape[1]
        self.cell, self.nu, self.events = cell, nu, events
        self.x = self.rng.normal(0, np.sqrt(0.5), (nu, NP))
        self.p = self.rng.normal(0, np.sqrt(0.5), (nu, NP))
        self.eps = np.ones(nu)
        self.nb = len(BONDS)
        self.C = self.rng.random((nu, self.nb, self.nq), dtype=np.float32)
        self.fires = 0
        self.peak = nu

    def rates(self, r):
        rt = self.rt
        u = (r - rt[0]) / (rt[1] - rt[0])
        i0 = np.clip(np.floor(u).astype(int), 0, rt.size - 2)
        t = (u - i0)[:, None]
        return ((1 - t) * self.kt[i0] + t * self.kt[i0 + 1]).astype(np.float32)

    def fire(self, dt):
        new = []
        for b, (i, j) in enumerate(BONDS):
            k = self.rates(self.x[:, i] - self.x[:, j])
            old = np.floor(self.C[:, b])
            self.C[:, b] += k * dt
            dn = (np.floor(self.C[:, b]) - old).astype(np.int64)
            w, q = np.nonzero(dn)
            if w.size == 0:
                continue
            cnt = np.abs(dn[w, q])
            w, q = np.repeat(w, cnt), np.repeat(q, cnt)
            tau = self.eps[w] * np.sign(dn[w, q])
            nx = np.concatenate([self.x[w], self.x[w]])
            npp = np.concatenate([self.p[w], self.p[w]])
            d = self.xi[q]
            m = w.size
            npp[:m, i] += d
            npp[:m, j] -= d
            npp[m:, i] -= d
            npp[m:, j] += d
            new.append((nx, npp, np.concatenate([tau, -tau])))
            self.fires += m
        if new:
            add = sum(a[0].shape[0] for a in new)
            self.x = np.concatenate([self.x] + [a[0] for a in new])
            self.p = np.concatenate([self.p] + [a[1] for a in new])
            self.eps = np.concatenate([self.eps] + [a[2] for a in new])
            self.C = np.concatenate([self.C, self.rng.random((add, self.nb, self.nq), dtype=np.float32)])

    def sink(self):
        for lev in range(self.levels - 1, -1, -1):      # finest cells first, then coarser
            self.sink_at(self.cell / 2**lev)

    def sink_at(self, cell):
        rho = self.x @ OJ[:NP - 1].T
        pr = self.p @ OJ[:NP - 1].T
        key = np.zeros(self.eps.size, np.int64)
        base = 4096 if NP == 3 else 256
        for a in list(rho.T) + list(pr.T):
            key = key * base + np.clip(np.floor(a / cell).astype(np.int64) + base // 2, 0, base - 1)
        n = key.size
        o = np.lexsort((self.rng.random(n), self.eps, key))
        ks, es = key[o], self.eps[o]
        start = np.concatenate([[0], np.flatnonzero(np.diff(ks)) + 1])
        cid = np.repeat(np.arange(start.size), np.diff(np.concatenate([start, [n]])))
        npos = np.bincount(cid, weights=(es > 0)).astype(np.int64)
        nneg = np.bincount(cid, weights=(es < 0)).astype(np.int64)
        m = np.minimum(npos, nneg)
        rank = np.arange(n) - np.where(es < 0, start[cid], start[cid] + nneg[cid])
        keep = np.ones(n, bool)
        keep[o[rank < m[cid]]] = False
        self.x, self.p, self.eps, self.C = self.x[keep], self.p[keep], self.eps[keep], self.C[keep]

    def run(self, T, dt=DT):
        rows = []
        every = int(round(OBS / dt))
        f = force(self.x, dv_gauss)
        t0 = time.perf_counter()
        for i in range(int(round(T / dt)) + 1):
            if i % every == 0:
                rows.append(dict(t=round(i * dt, 6), **observe(self.x, self.p, self.eps, self.nu, v_gauss)))
                if self.verbose:
                    print(f"    t = {i * dt:4.2f}: purity {rows[-1]['purity']:.4f}, worlds per nu "
                          f"{self.eps.size / self.nu:.2f}  ({time.perf_counter() - t0:.0f} s)", flush=True)
            if i == int(round(T / dt)):
                break
            self.p += 0.5 * dt * f
            self.x += self.p * dt
            f = force(self.x, dv_gauss)
            self.p += 0.5 * dt * f
            if self.events:
                self.fire(dt)
                if (i + 1) % self.sink_every == 0:
                    self.sink()
                f = force(self.x, dv_gauss)
                self.peak = max(self.peak, self.eps.size)
        return rows, time.perf_counter() - t0


def run_W(job):
    mult, cell, nu, seed, T, events, npart, dt, sink_every, levels = job
    set_particles(npart)
    w = Worlds(mult, cell, nu, seed, events, sink_every, levels)
    w.verbose = seed == 0
    rows, wall = w.run(T, dt)
    return dict(mult=mult, cell=cell, nu=nu, seed=seed, events=events, npart=npart, dt=dt,
                sink_every=sink_every, levels=levels, rows=rows, fires=w.fires, peak=w.peak, wall=wall)


def part_W(mult, cell, nu, seeds, T, workers, dt=DT, sink_every=1, levels=1):
    import multiprocessing as mp
    rule(f"Part W. {NP} particles, {len(BONDS)} bonds, worlds with one bank per bond: reach {mult} pi, "
         f"sink cell {cell}, nu = {nu}, {seeds} seed(s), t = 0 .. {T}")
    jobs = [(mult, cell, nu, s, T, True, NP, dt, sink_every, levels) for s in range(seeds)]
    res = []
    with mp.Pool(min(workers, len(jobs))) as pool:
        for r in pool.imap_unordered(run_W, jobs):
            f = r["rows"][-1]
            print(f"  seed {r['seed']}: t = {f['t']}: purity {f['purity']:.4f}, <Vpair> {f['Vpair']:.4f}, "
                  f"E {f['T'] + f['Vtrap'] + f['Vpair']:.4f}, worlds per nu {f['N'] / nu:.2f} (peak "
                  f"{r['peak'] / nu:.2f}), firings {r['fires']}  ({r['wall']:.0f} s)", flush=True)
            res.append(r)
    extra = (("" if dt == DT else f"_dt{dt}") + ("" if sink_every == 1 else f"_se{sink_every}")
             + ("" if levels == 1 else f"_L{levels}"))
    name = f"three_worlds_W{TAG()}_r{mult}_c{cell}_nu{nu}{extra}.json"
    dump(res, name)
    return res


def summarise_W(res, ex, mesh):
    rule("Part W against the mesh at the same reach and the exact TDSE")
    times = [r["t"] for r in res[0]["rows"]]
    mrows = {r["t"]: r for r in mesh["rows"]} if mesh else {}
    erows = {r["t"]: r for r in ex["rows"]} if ex else {}
    print("     t   purity: worlds (mean ± se)   mesh      exact    |  <Vpair>: worlds         mesh      exact"
          "   | N/nu")
    out = []
    for k, t in enumerate(times):
        a = np.array([[r["rows"][k]["purity"], r["rows"][k]["Vpair"], r["rows"][k]["N"] / r["nu"]] for r in res])
        m, se = a.mean(0), (a.std(0, ddof=1) / np.sqrt(len(a)) if len(a) > 1 else np.zeros(3))
        mp_, mv_ = (mrows[t]["purity"], mrows[t]["Vpair"]) if t in mrows else (np.nan, np.nan)
        ep_, ev_ = (erows[t]["purity"], erows[t]["Vpair"]) if t in erows else (np.nan, np.nan)
        out.append(dict(t=t, purity=m[0], purity_se=se[0], Vpair=m[1], Vpair_se=se[1], N=m[2]))
        if abs(t - round(t * 2) / 2) < 1e-9:
            print(f"  {t:4.2f}   {m[0]:.4f} ± {se[0]:.4f}            {mp_:.4f}    {ep_:.4f}   |  "
                  f"{m[1]:.4f} ± {se[1]:.4f}   {mv_:.4f}    {ev_:.4f}  | {m[2]:.2f}")
    return out


# ----------------------------------------------------------------------------
# Part X -- what biases the world realisation
# ----------------------------------------------------------------------------
def load_W_runs(npart=3):
    """Every saved world run for npart particles, with its settings (older files
    carry dt, sink frequency and sink levels only in the file name)."""
    tag = "" if npart == 3 else f"_N{npart}"
    runs = []
    for f in sorted(glob.glob(output_path(f"three_worlds_W{tag}_r*.json"))):
        name = os.path.basename(f)
        if npart == 3 and "_N" in name:
            continue
        res = json.load(open(f))
        for key, pre, typ, default in (("dt", "_dt", float, DT), ("sink_every", "_se", int, 1),
                                       ("levels", "_L", int, 1)):
            if key not in res[0]:
                val = default
                if pre in name:
                    val = typ(name.split(pre)[1].split("_")[0].replace(".json", ""))
                for r in res:
                    r[key] = val
        runs.append(res)
    return runs


def run_label(res):
    r = res[0]
    lab = f"reach {r['mult']}π, cell {r['cell']}, ν = {r['nu']:.0e}".replace("e+0", "e")
    if r["dt"] != DT:
        lab += f", Δt = {r['dt']}"
    if r["sink_every"] != 1:
        lab += f", sink every {r['sink_every']} steps"
    if r["levels"] != 1:
        lab += f", {r['levels']} sink levels"
    return lab


def part_X():
    rule("Part X. What biases the world realisation: worlds minus mesh at the same reach")
    meshes = {int(os.path.basename(f)[14:-5]): json.load(open(f))
              for f in glob.glob(output_path("three_worlds_M[0-9]*.json"))}
    runs = load_W_runs(3)
    E = load("three_worlds_E.json")
    cl = {round(r["t"], 6): r["purity"] for r in E["classical"]} if E else {}
    print("  q(t) = (worlds - classical) / (mesh - classical): the share of the quantum correction to the")
    print("  purity that the realisation keeps (1 = all of it, 0 = classical streaming)")
    print("  reach  cell   nu      dt     sink   levels seeds | N/nu(3) | purity - mesh: t=1     t=2     t=3   |"
          " q: t=1   t=2   t=3 | Vpair t=3 | global purity t=3")
    out = []
    for res in sorted(runs, key=lambda r: (r[0]["mult"], r[0]["levels"], r[0]["sink_every"], r[0]["dt"],
                                           r[0]["cell"], r[0]["nu"])):
        r0 = res[0]
        if r0["mult"] not in meshes:
            continue
        mr = {round(r["t"], 6): r for r in meshes[r0["mult"]]["rows"]}
        row = dict(mult=r0["mult"], cell=r0["cell"], nu=r0["nu"], dt=r0["dt"], sink_every=r0["sink_every"],
                   levels=r0["levels"], seeds=len(res))

        def diff(t, key):
            k = [round(r["t"], 6) for r in r0["rows"]].index(t)
            vals = [x["rows"][k].get(key, np.nan) for x in res]
            return float(np.mean(vals) - mr[t].get(key, np.nan))
        k3 = [round(r["t"], 6) for r in r0["rows"]].index(3.0)
        row["N3"] = float(np.mean([x["rows"][k3]["N"] / x["nu"] for x in res]))
        row["dpur"] = {t: diff(t, "purity") for t in (1.0, 2.0, 3.0)}
        row["q"] = {t: (row["dpur"][t] + mr[t]["purity"] - cl[t]) / (mr[t]["purity"] - cl[t]) if cl else np.nan
                    for t in (1.0, 2.0, 3.0)}
        row["dV3"] = diff(3.0, "Vpair")
        row["dgp3"] = diff(3.0, "gpurity") if "gpurity" in r0["rows"][k3] and "gpurity" in mr[3.0] else np.nan
        out.append(row)
        print(f"  {r0['mult']:3d}π  {r0['cell']:5.3f} {r0['nu']:8d}  {r0['dt']:5.3f}  {r0['sink_every']:4d}  "
              f"{r0['levels']:4d}  {len(res):4d}  |  {row['N3']:5.2f}  |  "
              + "  ".join(f"{row['dpur'][t]:+.4f}" for t in (1.0, 2.0, 3.0))
              + " | " + " ".join(f"{row['q'][t]:5.2f}" for t in (1.0, 2.0, 3.0))
              + f" |  {row['dV3']:+.4f}  |  {row['dgp3']:+.4f}")
    return out


# ----------------------------------------------------------------------------
# Part F -- figure
# ----------------------------------------------------------------------------
def part_F():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rule("Part F. Figure")
    E, Q, Q4 = load("three_worlds_E.json"), load("three_worlds_Q.json"), load("three_worlds_Q_N4.json")
    E4 = load("three_worlds_E_N4.json")
    meshes = {int(os.path.basename(f)[14:-5]): json.load(open(f))
              for f in glob.glob(output_path("three_worlds_M[0-9]*.json"))}
    runs = [r for r in load_W_runs(3) if r[0]["dt"] == DT and r[0]["sink_every"] == 1 and r[0]["levels"] == 1]
    runs4 = load_W_runs(4)
    col = {2: "#e67e22", 4: "#8e44ad", 8: "#16a085"}
    # at each reach, the largest-nu run at the standard cell 0.35
    best = {}
    for res in runs:
        m = res[0]["mult"]
        if res[0]["cell"] == 0.35 and (m not in best or res[0]["nu"] > best[m][0]):
            best[m] = (res[0]["nu"], res)
    fig, ax = plt.subplots(2, 3, figsize=(16, 8.8))
    for a_, key, lab in ((ax[0, 0], "purity", "(a) three particles, problem G: purity of particle 1, Tr ρ₁²"),
                         (ax[0, 1], "Vpair", "(b) three particles, problem G: pair energy ⟨Σ v⟩")):
        if E:
            a_.plot([r["t"] for r in E["rows"]], [r[key] for r in E["rows"]], "k-", lw=2, label="exact (TDSE 64³)")
            a_.plot([r["t"] for r in E["classical"]], [r[key] for r in E["classical"]], ":", color="#7f7f7f",
                    label="classical streaming only")
        for m in sorted(meshes):
            rr = meshes[m]["rows"]
            a_.plot([r["t"] for r in rr], [r[key] for r in rr], "--", color=col.get(m, "b"),
                    label=f"generator WC8 on the mesh, reach {m}π")
        for m in sorted(best):
            res = best[m][1]
            arr = np.array([[r[key] for r in x["rows"]] for x in res])
            a_.plot([r["t"] for r in res[0]["rows"]], arr.mean(0), "o", ms=4, mfc="none", color=col.get(m, "b"),
                    label="worlds, " + run_label(res))
        a_.set_xlabel("t")
        a_.set_title(lab, fontsize=9.5)
    ax[0, 0].legend(fontsize=6.5)
    # (c) global purity at resolution GA: the contact sink as a coarse-graining
    a_ = ax[0, 2]
    for m in sorted(meshes):
        rr = meshes[m]["rows"]
        if "gpurity" in rr[0]:
            a_.plot([r["t"] for r in rr], [r["gpurity"] for r in rr], "--", color=col.get(m, "b"),
                    label=f"mesh, reach {m}π")
    for res in sorted(runs, key=lambda r: (r[0]["mult"], r[0]["cell"], r[0]["nu"])):
        if "gpurity" not in res[0]["rows"][0]:
            continue
        m = res[0]["mult"]
        gp = [np.mean([x["rows"][k]["gpurity"] for x in res]) for k in range(len(res[0]["rows"]))]
        a_.plot([r["t"] for r in res[0]["rows"]], gp,
                {0.35: "o", 0.5: "s", 0.7: "^"}.get(res[0]["cell"], "D"), ms=4,
                mfc="none" if res[0]["nu"] < 10**6 else None,
                color=col.get(m, "b"), label="worlds, " + run_label(res))
    a_.set_xlabel("t")
    a_.set_title(f"(c) whole relative state seen at resolution {GA}: (2πħ)² ∫(W∗G)²", fontsize=9.5)
    a_.legend(fontsize=6.5)
    # (d) Moshinsky, three and four particles
    a_ = ax[1, 0]
    tf = np.linspace(0, T_END, 601)
    for q, n, c in ((Q, 3, "#1f6f8b"), (Q4, 4, "#c0392b")):
        if q:
            a_.plot(tf, moshinsky_purity(n, tf), "-", color=c, lw=1.5, label=f"N = {n}: exact (covariance flow)")
            a_.plot(q["times"], q["tdse"], "x", color=c, label=f"N = {n}: TDSE")
            a_.plot(q["times"], q["worlds"], "o", ms=7, mfc="none", color=c, label=f"N = {n}: worlds, no events")
    a_.set_title(f"(d) problem Q (Moshinsky, k = {K_MOSH}): purity; quadratic, so no events", fontsize=9.5)
    a_.set_xlabel("t")
    a_.legend(fontsize=7, loc="lower right", ncol=2)
    a_.set_ylim(0.5, 1.02)
    # (e) the realisation's bias at t = 3 against the sink cell
    a_ = ax[1, 1]
    mk = {10**5: "o", 2 * 10**5: "D", 3 * 10**5: "^", 10**6: "*"}
    for res in runs:
        m, c, nu = res[0]["mult"], res[0]["cell"], res[0]["nu"]
        if m not in meshes:
            continue
        k = [r["t"] for r in res[0]["rows"]].index(3.0)
        mp_ = [r["purity"] for r in meshes[m]["rows"] if r["t"] == 3.0][0]
        a_.plot(c, np.mean([x["rows"][k]["purity"] for x in res]) - mp_, mk.get(nu, "o"), color=col.get(m, "b"),
                ms=7 if nu < 10**6 else 11, mfc="none" if nu < 10**6 else col.get(m, "b"),
                label=f"reach {m}π, ν = {nu:.0e}".replace("e+0", "e"))
    # the largest-nu runs at each reach: a straight line in the cell, extended to c = 0
    for m in sorted(meshes):
        numax = max([r[0]["nu"] for r in runs if r[0]["mult"] == m], default=0)
        pts = []
        for res in runs:
            if res[0]["mult"] == m and res[0]["nu"] == numax:
                k = [r["t"] for r in res[0]["rows"]].index(3.0)
                mp_ = [r["purity"] for r in meshes[m]["rows"] if r["t"] == 3.0][0]
                pts.append((res[0]["cell"], np.mean([x["rows"][k]["purity"] for x in res]) - mp_))
        if len(pts) >= 3:
            c_, d_ = np.array(sorted(pts)).T
            fit = np.polyfit(c_, d_, 1)
            a_.plot([0, c_[-1]], np.polyval(fit, [0, c_[-1]]), ":", color=col.get(m, "b"),
                    label=f"reach {m}π, ν = {numax:.0e}: linear in the cell".replace("e+0", "e"))
    a_.axhline(0, color="k", lw=0.8)
    a_.set_xlim(0, None)
    a_.set_xlabel("sink cell (Jacobi coordinates)")
    a_.set_title("(e) worlds − mesh at the same reach: purity at t = 3", fontsize=9.5)
    h, l_ = a_.get_legend_handles_labels()
    keep = dict(zip(l_, h))
    a_.legend(keep.values(), keep.keys(), fontsize=7)
    # (f) four particles, problem G
    a_ = ax[1, 2]
    if E4:
        a_.plot([r["t"] for r in E4["rows"]], [r["purity"] for r in E4["rows"]], "k-", lw=2, label="exact (TDSE 40⁴)")
        a_.plot([r["t"] for r in E4["classical"]], [r["purity"] for r in E4["classical"]], ":", color="#7f7f7f",
                label="classical streaming only")
    for res in runs4:
        arr = np.array([[r["purity"] for r in x["rows"]] for x in res])
        a_.plot([r["t"] for r in res[0]["rows"]], arr.mean(0), "o", ms=4, mfc="none",
                color=col.get(res[0]["mult"], "b"),
                label="worlds, " + run_label(res))
    a_.set_xlabel("t")
    a_.set_title("(f) four particles, six bonds, problem G: purity of particle 1", fontsize=9.5)
    a_.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(output_path("three_worlds.png"), dpi=150, bbox_inches="tight")
    dp = docs_path("three_worlds.png")
    if dp:
        fig.savefig(dp, dpi=150, bbox_inches="tight")
    print("  wrote three_worlds.png")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("parts", nargs="*", default=["S", "K", "E", "Q", "M", "W", "X", "F"])
    ap.add_argument("--mesh-reach", type=int, nargs="*", default=[2, 4, 8])
    ap.add_argument("--reach", type=int, default=4)
    ap.add_argument("--cell", type=float, default=0.35)
    ap.add_argument("--nu", type=int, default=100000)
    ap.add_argument("--seeds", type=int, default=2)
    ap.add_argument("--T", type=float, default=T_END)
    ap.add_argument("--workers", type=int, default=os.cpu_count() or 1)
    ap.add_argument("--dt", type=float, default=DT, help="time step of the world runs")
    ap.add_argument("--sink-every", type=int, default=1, help="apply the contact sink every k steps")
    ap.add_argument("--levels", type=int, default=1, help="hierarchical sink: cells cell/2^(L-1) .. cell")
    ap.add_argument("--N", type=int, default=3, help="number of particles (3 or 4)")
    args = ap.parse_args()
    set_particles(args.N)
    parts = [p.upper() for p in args.parts]
    t0 = time.perf_counter()
    if "S" in parts:
        dump(part_S(), f"three_worlds_S{TAG()}.json")
    if "K" in parts:
        dump(part_K(), f"three_worlds_K{TAG()}.json")
    if "E" in parts:
        dump(part_E(), f"three_worlds_E{TAG()}.json")
    if "Q" in parts:
        dump(part_Q(), f"three_worlds_Q{TAG()}.json")
    if "M" in parts:
        assert NP == 3, "the Jacobi mesh is written for three particles"
        part_M(args.mesh_reach)
    if "W" in parts:
        res = part_W(args.reach, args.cell, args.nu, args.seeds, args.T, args.workers, args.dt,
                     args.sink_every, args.levels)
        summarise_W(res, load(f"three_worlds_E{TAG()}.json"),
                    load(f"three_worlds_M{args.reach}.json") if NP == 3 else None)
    if "X" in parts:
        dump(part_X(), "three_worlds_X.json")
    if "F" in parts:
        part_F()
    print(f"\nwall time {time.perf_counter() - t0:.1f} s")


if __name__ == "__main__":
    main()
