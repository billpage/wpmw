"""
Verification for ``docs/analysis/fourd_microdynamics.md``.

The four-action / sea-dressed / phase-alignment ladder was built in 1+1
dimensions.  This demo checks what survives, what generalises, and what
breaks when the joint phase space is four-dimensional -- two particles in
one spatial dimension (2p/1D) or one particle in two spatial dimensions
(1p/2D).

Parts
-----
A  Four-action exactness on the joint 2-D momentum lattice.  The symmetric
   member (5) and the pure-hop member (G = 0) both reproduce the joint QLE
   stencil for external, pair, and oblique mode geometries.
B  The leak law.  Focus/Defocus never change any linear momentum
   functional; a hop across n changes ``u . p`` by ``2 (u . q) dp``.  A
   direction is conserved iff it is orthogonal to every active mode
   wavevector.
C  Mode-dependent fibration of the momentum lattice, and the 2p/1D <-> 1p/2D
   generator identity.
D  The joint sea: the background factorises but the crystal shift does not
   commute with products; peak |W| / (2/h)^d equals the state purity;
   excess-excess co-location arithmetic.
E  The exchange vertex in d >= 2.  Theorem 4's two conditions leave a
   (d-1)-parameter family; energy is conserved by every member; only the
   permutation ("exchange-only") condition selects the swap.
F  The ring harmonic.  Residual of the exact QLE force term against the
   classical drift; the mode-independent noise constant 2 m w^2 hbar; the
   variance rate linear in the mode cutoff, measured against seeded runs.
G  Two particles in a harmonic trap with harmonic coupling: normal modes,
   exact stationarity of the ground-state Wigner function under the
   *classical* flow, entanglement entropy against the closed form.
H  Negativity transport: a cat state in the relative mode is rigidly
   rotated by the harmonic potential.  The microdynamics must create no
   negativity at all.

Output convention: figures go through ``wpmwlib.wpmw_utils.output_path`` and
``docs_path``.  Run with ``WPMW_OUTPUT`` set (``/mnt/user-data/outputs`` in the
Claude container).
"""

from __future__ import annotations

import itertools
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from wpmwlib.phase_space_crystal_lattice import (
    FourierMode,
    PhaseSpaceCrystalLattice,
)
from wpmwlib.wpmw_utils import docs_path, output_path

HBAR = 1.0
H_PLANCK = 2.0 * np.pi * HBAR
MASS = 1.0
SEED = 20260809

# Ring parameters shared by Parts F and H
L_RING = 8.0
OMEGA = 1.0
M_GRID, N_GRID = 64, 64
NU = 4.0e5


def save_fig(fig, name):
    fig.savefig(output_path(name), dpi=150, bbox_inches="tight")
    dp = docs_path(name)
    if dp:
        fig.savefig(dp, dpi=150, bbox_inches="tight")
    print(f"    [figure] {name}")


def banner(text):
    print()
    print("=" * 74)
    print(text)
    print("=" * 74)


# ===================================================================== #
# Part A -- four-action exactness on the joint 2-D momentum lattice     #
# ===================================================================== #

def roll2(W, q):
    """Shift operator A: ``(A W)[n] = W[n + q]`` on a periodic 2-D lattice."""
    return np.roll(W, shift=(-q[0], -q[1]), axis=(0, 1))


def neg(q):
    return (-q[0], -q[1])


def qle_stencil_2d(W, Gamma, modes):
    """Target generator: ``dW_n = sum_q Gamma_q (W_{n+q} - W_{n-q})``."""
    out = np.zeros_like(W)
    for q, G in zip(modes, Gamma):
        out += G * (roll2(W, q) - roll2(W, neg(q)))
    return out


def four_action_2d(W, Gamma, modes):
    """Symmetric member (5), channels assembled independently."""
    out = np.zeros_like(W)
    for q, G in zip(modes, Gamma):
        f = 0.5 * G * (roll2(W, q) - roll2(W, neg(q)))
        h = -0.5 * G * (roll2(W, q) + roll2(W, neg(q)))
        out += 2.0 * f - roll2(f, neg(q)) - roll2(f, q)
        out += roll2(h, neg(q)) - roll2(h, q)
    return out


def pure_hop_2d(W, Gamma, modes):
    """The ``G = 0`` member: ``f = 0``, ``h_n = -Gamma W_n``."""
    out = np.zeros_like(W)
    for q, G in zip(modes, Gamma):
        h = -G * W
        out += roll2(h, neg(q)) - roll2(h, q)
    return out


MODE_SETS = {
    "external, particle 1 only": [(1, 0)],
    "external, both particles": [(1, 0), (0, 1)],
    "pair potential (anti-diagonal)": [(1, -1), (2, -2)],
    "1p/2D oblique": [(1, 2), (3, -1)],
    "2p/1D mixed ext + pair": [(1, 0), (0, 1), (1, -1), (3, -3)],
}


def part_a(rng):
    banner("Part A -- four-action exactness on the joint 2-D momentum lattice")
    M = 12
    print(f"  lattice {M}x{M}, 4 random occupancy fields per mode set")
    print(f"  {'mode set':34s} {'symmetric (5)':>15s} {'pure-hop (G=0)':>16s}")
    worst = 0.0
    for name, modes in MODE_SETS.items():
        ws = wh = 0.0
        for _ in range(4):
            W = rng.normal(size=(M, M))
            Gam = rng.normal(size=len(modes))
            tgt = qle_stencil_2d(W, Gam, modes)
            ws = max(ws, np.max(np.abs(four_action_2d(W, Gam, modes) - tgt)))
            wh = max(wh, np.max(np.abs(pure_hop_2d(W, Gam, modes) - tgt)))
        print(f"  {name:34s} {ws:15.3e} {wh:16.3e}")
        worst = max(worst, ws, wh)
    print(f"\n  worst deviation over all mode geometries: {worst:.3e}")
    return worst


# ===================================================================== #
# Part B -- the leak law                                               #
# ===================================================================== #

def channel_events(q):
    """Per-event quantum displacements (cell units) of the four actions."""
    return [
        ("Focus", [(q, (0, 0)), (neg(q), (0, 0))]),
        ("Defocus", [((0, 0), q), ((0, 0), neg(q))]),
        ("Right-Hop", [(neg(q), q)]),
        ("Left-Hop", [(q, neg(q))]),
    ]


def leak(q, u):
    out = {}
    for name, moves in channel_events(q):
        out[name] = sum(np.dot(u, to) - np.dot(u, frm) for frm, to in moves)
    return out


def part_b():
    banner("Part B -- the leak law")
    print("  A hop across n changes  u . p  by  2 (u . q) dp;")
    print("  Focus/Defocus change it by 0, for every mode geometry.\n")
    cases = [
        ("2p/1D pair mode  q = (1,-1)", (1, -1),
         [("P = p1+p2", (1, 1)), ("p1-p2", (1, -1))]),
        ("2p/1D ext mode   q = (1, 0)", (1, 0),
         [("P = p1+p2", (1, 1)), ("p1-p2", (1, -1))]),
        ("1p/2D mode       q = (2, 0)", (2, 0),
         [("p_x", (1, 0)), ("p_y", (0, 1))]),
        ("1p/2D mode       q = (1, 2)", (1, 2),
         [("p_x", (1, 0)), ("p_y", (0, 1))]),
    ]
    ok = True
    for label, q, funcs in cases:
        print(f"  {label}")
        for fname, u in funcs:
            lk = leak(q, u)
            pred = 2 * int(np.dot(u, q))
            ok &= (lk["Focus"] == 0 and lk["Defocus"] == 0)
            ok &= (lk["Right-Hop"] == pred and lk["Left-Hop"] == -pred)
            body = "  ".join(f"{k}:{v:+d}" for k, v in lk.items())
            print(f"    u = {fname:10s}  2(u.q) = {pred:+d}   {body}")
    print(f"\n  leak law holds on every case tested: {ok}")
    return ok


# ===================================================================== #
# Part C -- fibration, and the 2p/1D <-> 1p/2D generator identity       #
# ===================================================================== #

def orbits(M, q):
    seen = np.zeros((M, M), dtype=bool)
    n_orb, lengths = 0, []
    for i in range(M):
        for j in range(M):
            if seen[i, j]:
                continue
            n_orb += 1
            a, b, ln = i, j, 0
            while not seen[a, b]:
                seen[a, b] = True
                a, b, ln = (a + q[0]) % M, (b + q[1]) % M, ln + 1
            lengths.append(ln)
    return n_orb, sorted(set(lengths))


def part_c(rng):
    banner("Part C -- fibration, and the 2p/1D <-> 1p/2D generator identity")
    M = 12
    print(f"  orbits of the shift-by-q map on a {M}x{M} momentum lattice:")
    print(f"  {'q':>10s} {'chains':>8s} {'chain lengths':>16s}")
    for q in [(1, 0), (0, 1), (1, -1), (2, -2), (1, 2), (3, -1), (6, 0)]:
        n_orb, ls = orbits(M, q)
        print(f"  {str(q):>10s} {n_orb:8d} {str(ls):>16s}")
    print("\n  Different modes fibre the SAME lattice along different")
    print("  directions -- there is no single momentum axis in 4-D.")

    # generator identity: 2p/1D with V2(x1 - x2) == 1p/2D with V(x,y)=V2(x-y)
    print("\n  2p/1D with a pair potential vs 1p/2D with V(x,y) = V2(x-y):")
    Mx = 8
    Mp = 10
    modes_2p = [(1, -1), (2, -2)]           # pair potential, joint wavevectors
    modes_1p = [(1, -1), (2, -2)]           # same wavevectors, read as (q_x,q_y)
    worst = 0.0
    for _ in range(3):
        W = rng.normal(size=(Mx, Mx, Mp, Mp))
        Gam = rng.normal(size=(len(modes_2p), Mx, Mx))
        g2 = np.zeros_like(W)
        g1 = np.zeros_like(W)
        for k, q in enumerate(modes_2p):
            Gx = Gam[k][:, :, None, None]
            g2 += Gx * (np.roll(W, (-q[0], -q[1]), axis=(2, 3))
                        - np.roll(W, (q[0], q[1]), axis=(2, 3)))
        for k, q in enumerate(modes_1p):
            Gx = Gam[k][:, :, None, None]
            g1 += Gx * (np.roll(W, (-q[0], -q[1]), axis=(2, 3))
                        - np.roll(W, (q[0], q[1]), axis=(2, 3)))
        worst = max(worst, np.max(np.abs(g2 - g1)))
    print(f"    max |generator_2p1D - generator_1p2D| = {worst:.3e}")
    print("    (identical by construction: the two readings differ only in")
    print("     which physical label is attached to each axis)")
    return worst


# ===================================================================== #
# Part D -- the joint sea                                              #
# ===================================================================== #

def gauss_wigner(x, p, sx, sp):
    return (1.0 / (2 * np.pi * sx * sp)) * np.exp(
        -np.asarray(x) ** 2 / (2 * sx**2) - np.asarray(p) ** 2 / (2 * sp**2))


def part_d():
    banner("Part D -- the joint sea")
    sx1, sp1 = 0.8, HBAR / (2 * 0.8)
    sx2, sp2 = 1.4, HBAR / (2 * 1.4)
    g = np.linspace(-3, 3, 33)
    X1, P1, X2, P2 = np.meshgrid(g, g, g, g, indexing="ij")
    W1 = gauss_wigner(X1, P1, sx1, sp1)
    W2 = gauss_wigner(X2, P2, sx2, sp2)
    bg1 = 2.0 / H_PLANCK
    bg2 = bg1**2
    joint_shift = W1 * W2 + bg2
    product_of_shifts = (W1 + bg1) * (W2 + bg1)
    predicted_gap = -bg1 * (W1 + W2)
    gap_err = np.max(np.abs(joint_shift - product_of_shifts - predicted_gap))
    print(f"  background: (2/h)^2 = {bg2:.6e} = (2/h)*(2/h) -> factorises")
    print("  but the SHIFT does not commute with the product:")
    print("    (W1 W2 + (2/h)^2) - (W1 + 2/h)(W2 + 2/h) = -(2/h)(W1 + W2)")
    print(f"    max residual of that identity: {gap_err:.3e}")
    print(f"    size of the gap / (2/h)^2:     "
          f"{np.max(np.abs(joint_shift - product_of_shifts)) / bg2:.4f}")

    print("\n  peak |W| / (2/h)^d equals the Gaussian state's purity:")
    for factor, tag in [(1.0, "pure"), (1.2, "mixed"), (3.0, "mixed")]:
        sx = 0.8
        sp = factor * HBAR / (2 * sx)
        peak = 1.0 / (2 * np.pi * sx * sp)
        purity = HBAR / (2 * sx * sp)
        print(f"    sp = {factor:.1f} x min ({tag:5s}): peak/(2/h) = "
              f"{peak / bg1:.6f},  purity = {purity:.6f}")

    print("\n  excess-excess co-location in the joint lattice (Mx = Mp = 32):")
    print(f"  {'dN':>4s} {'joint cells':>13s} {'worlds/cell at W=1e5':>22s} "
          f"{'at W=1e8':>12s}")
    for dN in [1, 2, 3, 4]:
        cells = float((32 * 32) ** dN)
        print(f"  {dN:4d} {cells:13.3e} {1e5 / cells:22.3e} "
              f"{1e8 / cells:12.3e}")
    print("  The sea, by construction, has B >> 1 per cell at every dN.")
    return gap_err


# ===================================================================== #
# Part E -- the exchange vertex in d >= 2                              #
# ===================================================================== #

def vertex_system(pa, pb):
    """Rows: momentum conservation (d eqs) + stationarity of mu (1 eq)."""
    d = len(pa)
    dvec = pa - pb
    A = np.zeros((d + 1, 2 * d))
    rhs = np.zeros(d + 1)
    for a in range(d):
        A[a, a] = -1.0
        A[a, d + a] = 1.0
        rhs[a] = dvec[a]
    A[d, :d] = dvec
    A[d, d:] = dvec
    rhs[d] = np.dot(dvec, pa + pb)
    return A, rhs


def is_permutation(inset, outset):
    for perm in itertools.permutations(range(len(inset))):
        if all(np.allclose(inset[perm[k]], outset[k]) for k in range(len(inset))):
            return True
    return False


def part_e(rng):
    banner("Part E -- the exchange vertex in d >= 2")
    print("  Conditions (Theorem 4):  p_out - p_in = p_a - p_b   [d equations]")
    print("                           d(mu)/dt = 0               [1 equation]")
    print("  Unknowns: (p_in, p_out) in R^{2d}.\n")
    print(f"  {'d':>3s} {'equations':>10s} {'unknowns':>9s} {'rank':>5s} "
          f"{'family dim':>11s}")
    for d in (1, 2, 3):
        dims = []
        for _ in range(4):
            pa, pb = rng.normal(size=d), rng.normal(size=d)
            A, _ = vertex_system(pa, pb)
            dims.append(2 * d - np.linalg.matrix_rank(A))
        print(f"  {d:3d} {d + 1:10d} {2 * d:9d} "
              f"{np.linalg.matrix_rank(A):5d} {str(set(dims)):>11s}")

    print("\n  d = 2, explicit family  p_in = p_b + t,  p_out = p_a + t,"
          "  t . Delta p = 0:")
    pa = np.array([0.7, -0.3])
    pb = np.array([-0.4, 0.9])
    dvec = pa - pb
    tperp = np.array([-dvec[1], dvec[0]]) / np.linalg.norm(dvec)
    A, rhs = vertex_system(pa, pb)
    vbar = (pa + pb) / (2 * MASS)
    print(f"  {'t':>7s} {'residual':>10s} {'dP':>10s} {'dE':>11s} "
          f"{'mu-dot':>11s} {'permutation':>12s}")
    rows = []
    for s in (0.0, 0.7, -1.9, 4.2):
        t = s * tperp
        p_in, p_out = pb + t, pa + t
        cand = np.concatenate([p_in, p_out])
        res = np.max(np.abs(A @ cand - rhs))
        Ein = (p_in @ p_in + pa @ pa + pb @ pb) / (2 * MASS)
        Eout = (p_out @ p_out + 2 * pb @ pb) / (2 * MASS)
        dP = np.max(np.abs((p_out + 2 * pb) - (p_in + pa + pb)))
        mudot = np.dot(dvec, (p_in + p_out) / (2 * MASS) - vbar) / HBAR
        perm = is_permutation([p_in, pa, pb], [p_out, pb, pb])
        print(f"  {s:+7.2f} {res:10.1e} {dP:10.1e} {Eout - Ein:+11.2e} "
              f"{mudot:+11.2e} {str(perm):>12s}")
        rows.append((s, perm))
    print("\n  Energy and momentum are conserved by EVERY member and none")
    print("  dephases.  Only the permutation test singles out the swap.")

    print("\n  In 2p/1D the transverse direction is the centre-of-mass axis:")
    dpair = np.array([1.0, -1.0]) / np.sqrt(2.0)     # pair-mode splitting
    tp = np.array([1.0, 1.0]) / np.sqrt(2.0)
    print(f"    Delta p direction (pair mode) = {dpair.round(4)}")
    print(f"    transverse direction          = {tp.round(4)}  "
          f"(dot = {np.dot(dpair, tp):.1e})")
    print("    -> the free parameter is a mismatch of P = p1 + p2 between")
    print("       the incoming world and the struck sea pair.")
    return rows


# ===================================================================== #
# Part F -- the ring harmonic                                          #
# ===================================================================== #

def ring_harmonic_coeffs(qmax, L, omega=OMEGA, m=MASS):
    q = np.arange(1, qmax + 1)
    Vq = 0.5 * m * omega**2 * (L**2 / np.pi**2) * ((-1.0) ** q) / q**2
    return q, Vq


def ring_modes(qmax, L):
    q, Vq = ring_harmonic_coeffs(qmax, L)
    return [FourierMode(q=int(qi), V_q=float(Vi), phi_q=0.0)
            for qi, Vi in zip(q, Vq)]


def ring_residual(L, sx, qmax=200):
    """max |exact ring QLE force term - classical drift| / |classical|_max."""
    dp = np.pi * HBAR / L
    sp = HBAR / (2 * sx)
    q, Vq = ring_harmonic_coeffs(qmax, L)
    xs = np.linspace(-L / 2, L / 2, 513, endpoint=False)
    ps = np.linspace(-4 * sp, 4 * sp, 65)
    XX, PP = np.meshgrid(xs, ps, indexing="xy")
    out = np.zeros_like(XX)
    for qi, Vi in zip(q, Vq):
        G = -(Vi / HBAR) * np.sin(2 * np.pi * qi * XX / L)
        out += G * (gauss_wigner(XX, PP + qi * dp, sx, sp)
                    - gauss_wigner(XX, PP - qi * dp, sx, sp))
    cls = MASS * OMEGA**2 * XX * (-(PP / sp**2)
                                  * gauss_wigner(XX, PP, sx, sp))
    prof = np.max(np.abs(out - cls), axis=0) / np.max(np.abs(cls))
    return np.max(prof), dp / sp, xs, prof


def make_ring_lattice(nu=NU, L=L_RING):
    lat = PhaseSpaceCrystalLattice(M=M_GRID, N=N_GRID, L=L, mass=MASS,
                                   hbar=HBAR, nu=nu)
    sx = np.sqrt(HBAR / (2 * MASS * OMEGA))
    sp = HBAR / (2 * sx)
    lat.initialize_from_wigner(lambda X, P: gauss_wigner(X, P, sx, sp))
    if nu is not None:
        lat.N_plus = np.round(
            (lat.W + lat.W_bg) * lat.nu * lat.dx * lat.dp).astype(np.int64)
    return lat


def variance_rate(lat, modes):
    """Exact tau-leap variance injected per unit time into total momentum."""
    bg = int(round(lat.W_bg * lat.nu * lat.dx * lat.dp))
    Nex = lat.N_plus.astype(np.float64) - bg
    total, per = 0.0, []
    for mode in modes:
        q = mode.q
        G = -(mode.V_q / lat.hbar) * np.sin(
            2.0 * np.pi * q * lat.x / lat.L + mode.phi_q)
        Hrate = -0.5 * G[None, :] * (np.roll(Nex, -q, axis=0)
                                     + np.roll(Nex, +q, axis=0))
        c = np.sum(np.abs(Hrate)) * (2 * q * lat.dp) ** 2
        per.append(c)
        total += c
    return total, np.array(per)


def part_f():
    banner("Part F -- the ring harmonic")
    sx = np.sqrt(HBAR / (2 * MASS * OMEGA))
    print("F1. exact ring QLE force term vs classical drift")
    print(f"    state: harmonic ground state, sigma_x = {sx:.4f}")
    print(f"    {'L':>7s} {'dp/sigma_p':>11s} {'max rel residual':>18s}")
    f1 = []
    for L in (4.0, 8.0, 16.0, 32.0, 64.0):
        r, ratio, xs, prof = ring_residual(L, sx)
        f1.append((L, ratio, r, xs, prof))
        print(f"    {L:7.1f} {ratio:11.4f} {r:18.4e}")
    print("    -> the ring reproduces classical drift only once the state")
    print("       spans many momentum cells (dp/sigma_p -> 0).")

    print("\nF2. per-mode momentum-noise constant")
    print("    hop jump = 2 q dp, rate amplitude |V_q|/hbar")
    print(f"    {'q':>4s} {'|V_q|':>12s} {'(2 q dp)^2':>12s} "
          f"{'product':>12s} {'/(2 m w^2 hbar)':>17s}")
    dp = np.pi * HBAR / L_RING
    q, Vq = ring_harmonic_coeffs(8, L_RING)
    ratios = []
    for qi, Vi in zip(q, Vq):
        prod = (abs(Vi) / HBAR) * (2 * qi * dp) ** 2
        ratios.append(prod / (2 * MASS * OMEGA**2 * HBAR))
        print(f"    {qi:4d} {abs(Vi):12.5e} {(2 * qi * dp)**2:12.5e} "
              f"{prod:12.5e} {ratios[-1]:17.10f}")
    print(f"    max deviation from 1: {np.max(np.abs(np.array(ratios) - 1)):.3e}")

    print("\nF3. tau-leap variance rate on the ground state vs mode cutoff")
    lat = make_ring_lattice()
    print(f"    {'qmax':>5s} {'Var rate':>13s} {'increment':>13s} "
          f"{'sum|V_q|/hbar':>15s}")
    prev = 0.0
    f3 = []
    for qmax in (1, 2, 4, 8, 16, 32):
        modes = ring_modes(qmax, L_RING)
        tot, _ = variance_rate(lat, modes)
        budget = sum(abs(m.V_q) for m in modes) / HBAR
        f3.append((qmax, tot))
        print(f"    {qmax:5d} {tot:13.5e} {(tot - prev) / max(qmax // 2, 1):13.5e} "
              f"{budget:15.4f}")
        prev = tot
    cos_modes = [FourierMode(q=1, V_q=1.5, phi_q=np.pi)]
    tot_cos, _ = variance_rate(lat, cos_modes)
    print(f"    single-mode cosine well (V_p = 1.5): {tot_cos:.5e}"
          f"   ({f3[-1][1] / tot_cos:.0f}x quieter than qmax=32)")

    print("\nF4. seeded runs: which observable sees the extra modes?")
    dt_jump, nstep, nseed = 0.002, 60, 8
    print(f"    {nseed} seeds, {nstep} jump substeps of dt = {dt_jump}")
    print(f"    {'qmax':>5s} {'relL2(W)':>10s} {'sd <p>':>12s} "
          f"{'sd <p^2>':>12s} {'predicted sd <p>':>17s}")
    f4 = []
    for qmax in (1, 2, 4, 8, 16, 32):
        t0 = time.time()
        modes = ring_modes(qmax, L_RING)
        mesh = make_ring_lattice(nu=None)
        for _ in range(nstep):
            mesh.step_jump_four_rule(modes, dt_jump)
        pcol = mesh.p[:, None]
        dV = mesh.dx * mesh.dp
        l2s, m1s, m2s = [], [], []
        for s in range(nseed):
            mc = make_ring_lattice()
            rng = np.random.default_rng(1000 + s)
            for _ in range(nstep):
                mc.step_jump_four_rule_mc(modes, dt_jump, rng)
            W = mc.N_plus / (NU * mc.dx * mc.dp) - mc.W_bg
            l2s.append(np.linalg.norm(W - mesh.W) / np.linalg.norm(mesh.W))
            nrm = np.sum(W) * dV
            m1s.append(np.sum(pcol * W) * dV / nrm)
            m2s.append(np.sum(pcol**2 * W) * dV / nrm)
        rate = dict(f3)[qmax]
        pred = np.sqrt(rate * nstep * dt_jump) / NU
        f4.append((qmax, np.mean(l2s), np.std(m1s), np.std(m2s), pred))
        print(f"    {qmax:5d} {np.mean(l2s):10.5f} {np.std(m1s):12.4e} "
              f"{np.std(m2s):12.4e} {pred:17.4e}   [{time.time() - t0:.0f}s]")
    print("    relL2(W) saturates; the momentum moments do not.")
    return f1, f3, f4


# ===================================================================== #
# Part G -- two particles in a harmonic trap                           #
# ===================================================================== #

def two_particle_matrices(om0, omc, m=MASS):
    K = m * np.array([[om0**2 + omc**2, -omc**2],
                      [-omc**2, om0**2 + omc**2]])
    Hm = np.zeros((4, 4))
    Hm[:2, :2] = K
    Hm[2:, 2:] = np.eye(2) / m
    Om = np.zeros((4, 4))
    Om[:2, 2:] = np.eye(2)
    Om[2:, :2] = -np.eye(2)
    return K, Hm, Om


def ground_covariance(K, m=MASS):
    ev, Vk = np.linalg.eigh(K)
    w = np.sqrt(ev / m)
    Sig = np.zeros((4, 4))
    Sig[:2, :2] = Vk @ np.diag(HBAR / (2 * m * w)) @ Vk.T
    Sig[2:, 2:] = Vk @ np.diag(HBAR * m * w / 2) @ Vk.T
    return Sig, w


def gaussian_entropy(nu):
    a, b = (nu + 1) / 2, (nu - 1) / 2
    return a * np.log(a) - (b * np.log(b) if b > 0 else 0.0)


def part_g():
    banner("Part G -- two particles in a harmonic trap with harmonic coupling")
    print("  H = (p1^2 + p2^2)/2m + m w0^2 (x1^2 + x2^2)/2"
          " + m wc^2 (x1 - x2)^2/2\n")
    print(f"  {'wc':>5s} {'w_+ (COM)':>10s} {'w_- (rel)':>10s} "
          f"{'stationarity':>13s} {'purity_1':>10s} {'S (num)':>10s} "
          f"{'S (closed)':>11s}")
    rows = []
    for omc in (0.0, 0.2, 0.5, 0.8, 1.5, 3.0):
        K, Hm, Om = two_particle_matrices(1.0, omc)
        A = Om @ Hm
        Sig, w = ground_covariance(K)
        stat = np.max(np.abs(A @ Sig + Sig @ A.T))
        Sig1 = np.array([[Sig[0, 0], Sig[0, 2]], [Sig[2, 0], Sig[2, 2]]])
        pur1 = (HBAR / 2) / np.sqrt(np.linalg.det(Sig1))
        nu_s = np.sqrt(np.linalg.det(Sig1)) / (HBAR / 2)
        S_num = gaussian_entropy(nu_s)
        wp, wm = 1.0, np.sqrt(1.0 + 2 * omc**2)
        xi = (np.sqrt(wm) - np.sqrt(wp)) / (np.sqrt(wm) + np.sqrt(wp))
        S_cl = (0.0 if xi == 0 else
                -np.log(1 - xi**2) - (xi**2 / (1 - xi**2)) * np.log(xi**2))
        rows.append((omc, pur1, S_num, S_cl))
        print(f"  {omc:5.1f} {w[0]:10.6f} {w[1]:10.6f} {stat:13.2e} "
              f"{pur1:10.6f} {S_num:10.6f} {S_cl:11.6f}")
    print("\n  'stationarity' = max |A Sigma + Sigma A^T| under the CLASSICAL")
    print("  flow: the joint ground-state Wigner function is a classical")
    print("  invariant, exactly, because H is quadratic.")
    print("  purity_1 is also the peak of the reduced Wigner function in")
    print("  units of 2/h -- the excess-to-background ratio of Part D.")
    return rows


# ===================================================================== #
# Part H -- negativity transport under a harmonic potential            #
# ===================================================================== #

def cat_wigner(X, P, x0, sx):
    sp = HBAR / (2 * sx)
    g = 1.0 / (2 * np.pi * sx * sp)
    plus = g * np.exp(-(X - x0) ** 2 / (2 * sx**2) - P**2 / (2 * sp**2))
    minus = g * np.exp(-(X + x0) ** 2 / (2 * sx**2) - P**2 / (2 * sp**2))
    cross = g * np.exp(-X**2 / (2 * sx**2) - P**2 / (2 * sp**2)) \
        * np.cos(2.0 * x0 * P / HBAR)
    norm = 2.0 * (1.0 + np.exp(-x0**2 / (2 * sx**2)))
    return (plus + minus + 2.0 * cross) / norm


def jump_symbol(lat, modes):
    """Momentum-Fourier symbol of the jump generator.

    At fixed x the generator ``dW_n = sum_q Gamma_q(x)(W_{n+q} - W_{n-q})``
    is diagonal in the momentum-Fourier basis with the purely imaginary
    eigenvalue ``i lam(x, theta)``,
    ``lam = 2 sum_q Gamma_q(x) sin(q theta)``.  All modes commute, being
    functions of the same shift operator, so the substep is exactly
    integrable: ``exp(dt * generator)`` is multiplication by
    ``exp(i lam dt)``.  Returns ``lam`` with shape ``(N, M)``.
    """
    theta = np.fft.fftfreq(lat.N) * 2 * np.pi
    sym = np.zeros((lat.N, lat.M))
    for m in modes:
        G = -(m.V_q / lat.hbar) * np.sin(
            2 * np.pi * m.q * lat.x / lat.L + m.phi_q)
        sym += 2.0 * np.sin(m.q * theta)[:, None] * G[None, :]
    return sym


def exact_jump_step(lat, sym, dt):
    Wh = np.fft.fft(lat.W, axis=0)
    lat.W = np.real(np.fft.ifft(Wh * np.exp(1j * sym * dt), axis=0))


def cat_lattice(L, M, N, x0, sx):
    lat = PhaseSpaceCrystalLattice(M=M, N=N, L=L, mass=MASS, hbar=HBAR,
                                   nu=None, advection="spectral")
    lat.initialize_from_wigner(lambda X, P: cat_wigner(X, P, x0, sx))
    return lat


def exact_rotated_cat(X, P, t, x0, sx):
    c, s = np.cos(OMEGA * t), np.sin(OMEGA * t)
    return cat_wigner(c * X - (s / (MASS * OMEGA)) * P,
                      c * P + MASS * OMEGA * s * X, x0, sx)


def negativity_of(lat):
    return float(np.sum(np.abs(np.minimum(lat.W, 0.0))) * lat.dx * lat.dp)


def part_h():
    banner("Part H -- negativity transport under a harmonic potential")
    sx = np.sqrt(HBAR / (2 * MASS * OMEGA))
    x0 = 1.6
    T = 2 * np.pi / OMEGA

    print("H1. explicit Euler on the jump substep amplifies.")
    print("    The generator's symbol is purely imaginary, so Euler grows the")
    print("    norm by sqrt(1 + lam^2 dt^2) per step: exp(lam^2 dt T / 2)"
          " over T.")
    L, M, N, qmax, n_macro = 8.0, 64, 64, 24, 20
    modes = ring_modes(qmax, L)
    print(f"    L = {L}, grid {M}x{N}, qmax = {qmax}, one period")
    print(f"    {'n_sub':>7s} {'dt_sub':>10s} {'lam dt':>9s} "
          f"{'neg/neg0':>10s} {'predicted':>10s}")
    for n_sub in (48, 192, 768):
        lat = cat_lattice(L, M, N, x0, sx)
        sym = jump_symbol(lat, modes)
        lam = float(np.max(np.abs(sym)))
        dt = T / n_macro
        dts = dt / n_sub
        neg0 = negativity_of(lat)
        for _ in range(n_macro):
            lat.step_advect(dt / 2)
            for _ in range(n_sub):
                lat.step_jump_fourier(modes, dts)
            lat.step_advect(dt / 2)
        print(f"    {n_sub:7d} {dts:10.3e} {lam * dts:9.5f} "
              f"{negativity_of(lat) / neg0:10.5f} "
              f"{np.exp(lam**2 * dts * T / 2):10.5f}")
    print("    -> the growth is a property of the integrator, not of the QLE.")

    print("\nH2. the same run with the substep integrated exactly")
    print(f"    {'L':>5s} {'M':>5s} {'N':>5s} {'qmax':>5s} {'steps':>7s} "
          f"{'neg/neg0':>10s} {'relL2 vs rotation':>18s} {'norm':>9s}")
    hist, snaps, lat_keep = None, None, None
    for L, M, N, qmax, nm, keep in [
        (8.0, 64, 64, 24, 800, False),
        (16.0, 128, 128, 48, 200, False),
        (16.0, 128, 128, 48, 800, True),
        (24.0, 192, 192, 72, 800, False),
    ]:
        lat = cat_lattice(L, M, N, x0, sx)
        modes = ring_modes(qmax, L)
        sym = jump_symbol(lat, modes)
        dt = T / nm
        dA = lat.dx * lat.dp
        neg0 = negativity_of(lat)
        rows = [(0.0, neg0, float(np.sum(lat.W) * dA))]
        snap = {0.0: lat.W.copy()}
        for k in range(1, nm + 1):
            lat.step_advect(dt / 2)
            exact_jump_step(lat, sym, dt)
            lat.step_advect(dt / 2)
            if keep and k % max(1, nm // 60) == 0:
                rows.append((k * dt, negativity_of(lat),
                             float(np.sum(lat.W) * dA)))
            if keep and k in (nm // 4, nm // 2, nm):
                snap[k * dt] = lat.W.copy()
        Wex = exact_rotated_cat(lat.X, lat.P, T, x0, sx)
        err = np.linalg.norm(lat.W - Wex) / np.linalg.norm(Wex)
        print(f"    {L:5.0f} {M:5d} {N:5d} {qmax:5d} {nm:7d} "
              f"{negativity_of(lat) / neg0:10.6f} {err:18.4e} "
              f"{float(np.sum(lat.W) * dA):9.6f}")
        if keep:
            hist = np.array(rows)
            snaps = snap
            lat_keep = lat
    print("    L = 8 puts the cat's tail on the seam of the periodic image;")
    print("    at L = 16 the harmonic ring rotates the state rigidly and")
    print("    negativity is preserved to 2e-4 -- all of it splitting error.")
    return hist, snaps, lat_keep


# ===================================================================== #
# Figures                                                              #
# ===================================================================== #

def fig_mode_geometry():
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    lattice_pts = [(i, j) for i in range(-4, 5) for j in range(-4, 5)]

    for ax, q, u, ulab, title in [
        (axes[0], (1, 0), (1, 1), r"$P=p_1+p_2$",
         "2p/1D external mode $q=(1,0)$\n(trap acting on particle 1)"),
        (axes[1], (1, -1), (1, 1), r"$P=p_1+p_2$",
         "2p/1D pair mode $q=(1,-1)$\n($V_2(x_1-x_2)$)"),
        (axes[2], (1, 2), (2, -1), r"$2p_x-p_y$",
         "1p/2D oblique mode $q=(1,2)$"),
    ]:
        xs = [a for a, b in lattice_pts]
        ys = [b for a, b in lattice_pts]
        ax.plot(xs, ys, ".", color="0.8", ms=5, zorder=1)
        qv = np.array(q, dtype=float)
        perp = np.array([-qv[1], qv[0]])
        perp = perp / np.linalg.norm(perp)
        off = np.round(1.9 * perp)             # draw the hop on its own row

        # Focus stencil: two quanta from n +- q converge on the centre n
        for sgn in (+1, -1):
            src = sgn * qv
            ax.annotate("", xy=(0.22 * -sgn * qv[0], 0.22 * -sgn * qv[1]),
                        xytext=tuple(src),
                        arrowprops=dict(arrowstyle="-|>", color="crimson",
                                        lw=2.4, shrinkA=6, shrinkB=0))
            ax.plot(*src, "o", mfc="none", mec="crimson", mew=2.0, ms=11,
                    zorder=3)
        ax.plot(0, 0, "o", color="crimson", ms=12, zorder=4)
        ax.text(0.25, -0.75, "Focus", color="crimson", fontsize=10,
                ha="left", weight="bold")

        # Hop stencil: one quantum from n-q to n+q, drawn on an offset row
        ax.annotate("", xy=tuple(off + qv), xytext=tuple(off - qv),
                    arrowprops=dict(arrowstyle="-|>", color="royalblue",
                                    lw=2.4, shrinkA=6, shrinkB=6,
                                    linestyle="--"))
        ax.plot(*(off - qv), "o", color="royalblue", ms=11, zorder=3)
        ax.plot(*(off + qv), "o", mfc="none", mec="royalblue", mew=2.0,
                ms=11, zorder=3)
        ax.text(off[0] + 0.25, off[1] + 0.55, "Hop", color="royalblue",
                fontsize=10, ha="left", weight="bold")

        # conserved-direction candidate
        leak_val = 2 * (u[0] * q[0] + u[1] * q[1])
        nu = np.array(u, dtype=float) / np.linalg.norm(u)
        ax.plot([-3.4 * nu[0], 3.4 * nu[1] * 0 + 3.4 * nu[0]],
                [-3.4 * nu[1], 3.4 * nu[1]], "-", color="seagreen", lw=1.8,
                alpha=0.85, zorder=2)
        ax.text(0.03, 0.97,
                f"$u = ${ulab}\n$u\\cdot q = {u[0]*q[0] + u[1]*q[1]:+d}$\n"
                f"hop leak $= {leak_val:+d}\\,\\Delta p$\n"
                f"focus leak $= 0$",
                transform=ax.transAxes, fontsize=9, va="top", ha="left",
                color="darkgreen",
                bbox=dict(boxstyle="round,pad=0.35", fc="honeydew",
                          ec="seagreen", alpha=0.95))
        ax.set_xlim(-3.6, 3.6)
        ax.set_ylim(-3.6, 3.6)
        ax.set_aspect("equal")
        ax.set_xlabel("$n_1$" if q != (1, 2) else "$n_x$")
        ax.set_ylabel("$n_2$" if q != (1, 2) else "$n_y$")
        ax.set_title(title, fontsize=10)
        ax.grid(alpha=0.2)

    leg = [
        plt.Line2D([], [], color="crimson", lw=2.2,
                   label="Focus: two quanta $\\to$ centre (never leaks)"),
        plt.Line2D([], [], color="royalblue", lw=2.2, ls="--",
                   label="Hop: one quantum by $2q$ (leaks $2(u\\cdot q)\\Delta p$)"),
        plt.Line2D([], [], color="seagreen", lw=1.6,
                   label="candidate conserved direction $u$"),
    ]
    axes[1].legend(handles=leg, fontsize=8, loc="upper center",
                   bbox_to_anchor=(0.5, -0.16), ncol=1, framealpha=0.9)
    fig.suptitle("Four-action stencils on the joint momentum lattice: "
                 "the mode wavevector, not the particle count, "
                 "decides what is conserved", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    save_fig(fig, "fourd_microdynamics_mode_geometry.png")
    plt.close(fig)


def fig_vertex_family():
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.2))
    pa = np.array([0.7, -0.3])
    pb = np.array([-0.4, 0.9])
    dvec = pa - pb
    tperp = np.array([-dvec[1], dvec[0]]) / np.linalg.norm(dvec)

    ax = axes[0]
    ax.plot(*pa, "o", color="crimson", ms=10, label=r"$p_a$ (struck partner)")
    ax.plot(*pb, "o", color="royalblue", ms=10, label=r"$p_b$ (mate)")
    ss = np.linspace(-2.2, 2.2, 200)
    lin_in = pb[None, :] + ss[:, None] * tperp[None, :]
    lin_out = pa[None, :] + ss[:, None] * tperp[None, :]
    ax.plot(lin_in[:, 0], lin_in[:, 1], "-", color="royalblue", lw=1.5,
            alpha=0.6, label=r"admissible $p_{\rm in}$")
    ax.plot(lin_out[:, 0], lin_out[:, 1], "-", color="crimson", lw=1.5,
            alpha=0.6, label=r"admissible $p_{\rm out}$")
    for s, mk in [(0.0, "*"), (0.9, "s"), (-1.6, "^")]:
        pin, pout = pb + s * tperp, pa + s * tperp
        ax.plot(*pin, mk, color="royalblue", ms=11 if s == 0 else 7,
                mec="k", mew=0.6)
        ax.plot(*pout, mk, color="crimson", ms=11 if s == 0 else 7,
                mec="k", mew=0.6)
        ax.annotate("", xy=pout, xytext=pin,
                    arrowprops=dict(arrowstyle="->", color="0.35",
                                    lw=1.6 if s == 0 else 1.0,
                                    linestyle="-" if s == 0 else ":"))
    ax.annotate(r"swap ($t=0$)", xy=pb, xytext=(pb[0] - 1.35, pb[1] + 0.55),
                fontsize=10, color="k",
                arrowprops=dict(arrowstyle="->", color="0.4", lw=1))
    ax.set_xlabel(r"$p^{(1)}$")
    ax.set_ylabel(r"$p^{(2)}$")
    ax.set_title("The exchange vertex in $d=2$:\n"
                 "a one-parameter family, not a unique swap", fontsize=11)
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(alpha=0.3)
    ax.set_aspect("equal")

    ax = axes[1]
    ss = np.linspace(-2.5, 2.5, 121)
    dE, dP, mudot, perm = [], [], [], []
    for s in ss:
        t = s * tperp
        pin, pout = pb + t, pa + t
        dE.append(((pout @ pout + 2 * pb @ pb)
                   - (pin @ pin + pa @ pa + pb @ pb)) / (2 * MASS))
        dP.append(np.max(np.abs((pout + 2 * pb) - (pin + pa + pb))))
        mudot.append(np.dot(dvec, (pin + pout) / (2 * MASS)
                            - (pa + pb) / (2 * MASS)) / HBAR)
        perm.append(1.0 if is_permutation([pin, pa, pb], [pout, pb, pb])
                    else 0.0)
    ax.semilogy(ss, np.abs(dE) + 1e-18, label=r"$|\Delta E|$", lw=2)
    ax.semilogy(ss, np.abs(dP) + 1e-18, label=r"$|\Delta P|$", lw=2, ls="--")
    ax.semilogy(ss, np.abs(mudot) + 1e-18, label=r"$|\dot\mu|$", lw=2, ls=":")
    ax.axvline(0.0, color="seagreen", lw=3, alpha=0.55,
               label="exchange-only (permutation) holds:\nthe single point $t=0$")
    ax.set_ylim(1e-18, 1e2)
    ax.set_xlabel(r"transverse parameter $t$ (units of $\hat t_\perp$)")
    ax.set_title("Momentum, energy and stationarity are blind to $t$;\n"
                 "only the permutation condition is not", fontsize=11)
    ax.legend(fontsize=9, loc="center right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    save_fig(fig, "fourd_microdynamics_vertex_family.png")
    plt.close(fig)


def fig_harmonic_cost(f1, f3, f4):
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.8))

    ax = axes[0]
    for L, ratio, r, xs, prof in f1:
        ax.semilogy(xs / L, prof + 1e-18, lw=1.6,
                    label=f"$L={L:.0f}$, $\\Delta p/\\sigma_p={ratio:.2f}$")
    ax.set_xlabel("$x/L$")
    ax.set_ylabel("relative residual")
    ax.set_ylim(1e-17, 3)
    ax.set_title("Ring harmonic: exact QLE force term\n"
                 "minus classical drift", fontsize=11)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[1]
    qm = np.array([a for a, b in f3])
    vr = np.array([b for a, b in f3])
    ax.plot(qm, vr, "o-", lw=2, color="crimson", label="tau-leap variance rate")
    ax.plot(qm, vr[0] * qm, "--", color="0.4", lw=1.4,
            label=r"$\propto q_{\max}$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"mode cutoff $q_{\max}$")
    ax.set_ylabel("Var$(\\Delta\\Pi)$ per unit time")
    ax.set_title("Injected momentum-space noise grows\n"
                 "linearly in the mode count", fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3, which="both")

    ax = axes[2]
    q4 = np.array([r[0] for r in f4])
    l2 = np.array([r[1] for r in f4])
    m1 = np.array([r[2] for r in f4])
    m2 = np.array([r[3] for r in f4])
    pr = np.array([r[4] for r in f4])
    ax.loglog(q4, l2 / l2[0], "o-", lw=2, label=r"rel$L^2(W)$")
    ax.loglog(q4, m1 / m1[0], "s-", lw=2, label=r"sd $\langle p\rangle$")
    ax.loglog(q4, m2 / m2[0], "^-", lw=2, label=r"sd $\langle p^2\rangle$")
    ax.loglog(q4, pr / pr[0], "--", color="0.4", lw=1.4,
              label="predicted (Part F3)")
    ax.set_xlabel(r"mode cutoff $q_{\max}$")
    ax.set_ylabel("noise, relative to $q_{\\max}=1$")
    ax.set_title("The field norm saturates;\n"
                 "momentum moments do not", fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3, which="both")

    fig.suptitle("A harmonic potential is exactly classical in the mean and "
                 "maximally noisy in the microdynamics", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    save_fig(fig, "fourd_microdynamics_harmonic_cost.png")
    plt.close(fig)


def fig_two_particle(rows, hist, snaps, lat):
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.8))

    ax = axes[0]
    wc = np.array([r[0] for r in rows])
    pur = np.array([r[1] for r in rows])
    Sn = np.array([r[2] for r in rows])
    Sc = np.array([r[3] for r in rows])
    ax.plot(wc, pur, "o-", lw=2, color="royalblue",
            label=r"reduced purity $=\max W_1/(2/h)$")
    ax2 = ax.twinx()
    ax2.plot(wc, Sn, "s-", lw=2, color="crimson", label="entropy (numerical)")
    ax2.plot(wc, Sc, "--", lw=1.5, color="k", label="entropy (closed form)")
    ax.set_xlabel(r"coupling $\omega_c$")
    ax.set_ylabel("reduced purity", color="royalblue")
    ax2.set_ylabel("entanglement entropy (nats)", color="crimson")
    ax.set_title("Two coupled oscillators: entanglement\n"
                 "is excess-to-background loss", fontsize=11)
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=8, loc="center right")
    ax.grid(alpha=0.3)

    ax = axes[1]
    ts = [t for t in sorted(snaps)]
    W = snaps[ts[len(ts) // 2]]
    vmax = np.max(np.abs(W))
    cf = ax.contourf(lat.x, lat.p, W, levels=np.linspace(-vmax, vmax, 21),
                     cmap="RdBu_r", extend="both")
    ax.set_xlabel("$x$")
    ax.set_ylabel("$p$")
    ax.set_ylim(-3, 3)
    ax.set_title(f"Cat state in the harmonic ring,\n$t = {ts[len(ts)//2]:.2f}$",
                 fontsize=11)
    fig.colorbar(cf, ax=ax, fraction=0.046)

    ax = axes[2]
    ax.plot(hist[:, 0], hist[:, 1] / hist[0, 1], "-", lw=2, color="crimson",
            label="negativity / initial")
    ax.plot(hist[:, 0], hist[:, 2] / hist[0, 2], "-", lw=1.5, color="seagreen",
            label="norm / initial")
    ax.axhline(1.0, color="0.4", ls="--", lw=1)
    ax.set_xlabel("$t$")
    ax.set_ylabel("ratio to initial value")
    ax.set_title("A harmonic potential transports negativity\n"
                 "without creating any", fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    fig.tight_layout()
    save_fig(fig, "fourd_microdynamics_two_particle_harmonic.png")
    plt.close(fig)


# ===================================================================== #
def main():
    rng = np.random.default_rng(SEED)
    print("Verification for docs/analysis/fourd_microdynamics.md")
    print(f"seed = {SEED}")
    part_a(rng)
    part_b()
    part_c(rng)
    part_d()
    part_e(rng)
    f1, f3, f4 = part_f()
    rows = part_g()
    hist, snaps, lat = part_h()

    banner("Figures")
    fig_mode_geometry()
    fig_vertex_family()
    fig_harmonic_cost(f1, f3, f4)
    fig_two_particle(rows, hist, snaps, lat)
    print("\ndone.")


if __name__ == "__main__":
    main()
