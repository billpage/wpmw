"""
Verification for ``docs/analysis/open_position_space.md``.

What happens to the microdynamics when position space is not closed.

The periodic box length ``L`` has been doing two unrelated jobs: fixing the
momentum quantum ``dp = pi hbar / L`` and keeping worlds in view.  This
script separates them and measures the consequences.

Parts
-----
A  The Wigner kernel has no position envelope.  ``V_W(x, xi)`` has modulus
   independent of x for every V; localisation of V appears only as a phase
   ``exp(2 i xi x / hbar)`` whose fringe spacing in xi is ``pi hbar / (2x)``.
   The total jump rate saturates at ``(2 / pi^2 hbar) int |Vtilde| dk``.
B  Signal against noise.  The first moment of the kernel is exactly
   ``-V'(x)``; the rate and the momentum-variance injection are flat.  The
   signal-to-noise ratio therefore decays like V' itself.
C  The coherence horizon.  Truncating the ket-bra separation at ``L_c``
   confines all jump activity to within ``L_c/2`` of supp V, preserves signed
   number exactly (the kernel stays odd in xi), and costs nothing once it clears the
   support of V.  It is the wrong tool for a potential that does not
   decay: a sharp window there leaves 1/xi tails and an unbounded budget.
D  Periodicity of V, not compactness of x.  For an a-periodic V every jump is
   a multiple of ``pi hbar / a``, so ``p mod (pi hbar / a)`` is conserved
   exactly and the dynamics is block diagonal in the residue class.  The
   decomposition is verified to be both invariant and complete.
E  Momentum-grid refinement is sector multiplication, not coupling: the
   spec's ``dp = pi hbar / (K L)`` option splits the lattice into K
   non-interacting sub-lattices.
F  Almost-periodic V: the reachable momentum set is a rank-r Z-module.  The
   crystal melts into a quasicrystal and then a continuum.
G  A complex absorbing potential uses the same stencil geometry with the
   opposite sign rule: the commutator takes the difference of the two offset
   rows, the anticommutator their sum.

Run with ``WPMW_OUTPUT`` set (``/mnt/user-data/outputs`` in the container).
"""

from __future__ import annotations

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

# The reference scatterer: a Gaussian barrier, localised and smooth so that
# every claim about "far from the potential" has an unambiguous meaning.
AMP, SIG = 1.0, 0.5


def banner(text):
    print()
    print("=" * 74)
    print(text)
    print("=" * 74)


def save_fig(fig, name):
    fig.savefig(output_path(name), dpi=150, bbox_inches="tight")
    dp = docs_path(name)
    if dp:
        fig.savefig(dp, dpi=150, bbox_inches="tight")
    print(f"    [figure] {name}")


# --------------------------------------------------------------------- #
# The Wigner kernel                                                      #
# --------------------------------------------------------------------- #
def V(x):
    return AMP * np.exp(-x ** 2 / (2 * SIG ** 2))


def dV(x):
    return -x / SIG ** 2 * V(x)


def V_hat(k):
    """int V(x) e^{-i k x} dx."""
    return AMP * SIG * np.sqrt(2 * np.pi) * np.exp(-SIG ** 2 * k ** 2 / 2)


def kernel(x, xi, Lc=None, dy=2.0e-3, pot=None, ywide=None):
    """The Wigner potential

        V_W(x, xi) = (1 / i pi hbar^2) int dy e^{-2 i xi y / hbar}
                         [V(x + y) - V(x - y)]

    by direct quadrature.  ``y`` is the *half* ket-bra separation, so a
    coherence horizon ``|x_ket - x_bra| <= Lc`` is the window ``|y| <= Lc/2``.
    ``Lc=None`` is the untruncated kernel.  The integrand is odd in y, so the
    result is real and the transform reduces to a sine transform.
    """
    pot = V if pot is None else pot
    if Lc is not None:
        yhi = Lc / 2.0
    else:
        yhi = (abs(x) + 6 * SIG) if ywide is None else ywide
    y = np.arange(0.0, yhi + dy, dy)
    g = pot(x + y) - pot(x - y)
    out = np.empty_like(xi, dtype=float)
    chunk = max(1, int(2.0e7 // len(y)))
    for i in range(0, len(xi), chunk):
        s = np.sin(2.0 * np.outer(xi[i:i + chunk], y) / HBAR)
        out[i:i + chunk] = np.trapezoid(s * g[None, :], y, axis=1)
    return -2.0 * out / (np.pi * HBAR ** 2)


def xi_nodes(Lc=None, xi_max=40.0, n_cont=4001):
    """Momentum-transfer nodes.  A horizon Lc makes them a lattice."""
    if Lc is None:
        xi = np.linspace(-xi_max, xi_max, n_cont)
        return xi, xi[1] - xi[0]
    dxi = np.pi * HBAR / Lc
    a = int(xi_max / dxi)
    n = np.arange(-a - 1, a + 1)          # symmetric about zero
    return (n + 0.5) * dxi, dxi


def moments(x, Lc=None, xi_max=40.0, pot=None, ywide=None, dy=2.0e-3):
    """(signal, rate, variance-injection, norm-defect) of the kernel at x."""
    xi, w = xi_nodes(Lc, xi_max)
    k = kernel(x, xi, Lc=Lc, pot=pot, ywide=ywide, dy=dy)
    signal = abs(np.sum(xi * k) * w)          # equals |V'(x)| exactly
    rate = np.sum(np.abs(k)) * w              # total jump rate per world
    diff = np.sum(xi ** 2 * np.abs(k)) * w    # momentum-variance injection
    defect = abs(np.sum(k) * w)               # must be 0: signed number
    return signal, rate, diff, defect


# --------------------------------------------------------------------- #
# Part A                                                                 #
# --------------------------------------------------------------------- #
def part_a():
    banner("A. The Wigner kernel has no position envelope")
    xi = np.linspace(-20, 20, 2001)
    print(f"  {'x':>6} {'max|V_W|':>12} {'V(x)':>12} {'|V.(x)|':>12}")
    amps = []
    xs = [0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0]
    for x in xs:
        k = kernel(x, xi)
        amps.append(np.max(np.abs(k)))
        print(f"  {x:6.1f} {amps[-1]:12.6f} {V(x):12.3e} {abs(dV(x)):12.3e}")
    print("\n  V and V' fall by 200+ orders of magnitude; the kernel does not.")

    print("\n  Fringe spacing in xi versus the prediction pi hbar / (2 x):")
    xi_f = np.linspace(0.01, 3.0, 24001)
    rows = []
    for x in [4.0, 8.0, 12.0]:
        k = kernel(x, xi_f)
        z = xi_f[np.where(np.sign(k[:-1]) != np.sign(k[1:]))[0]]
        meas = float(np.mean(np.diff(z[:8])))
        pred = np.pi * HBAR / (2 * x)
        rows.append((x, meas, pred))
        print(f"    x = {x:5.1f}   measured {meas:9.6f}   "
              f"predicted {pred:9.6f}   rel {abs(meas/pred-1):.1e}")

    # closed form for the saturated rate: mean of |sin| is 2/pi
    kk = np.linspace(-200, 200, 200001)
    r_inf = (2.0 / (np.pi ** 2 * HBAR)) * np.trapezoid(np.abs(V_hat(kk)), kk)
    _, r_meas, _, _ = moments(12.0)
    print(f"\n  Saturated rate: measured {r_meas:.6f} at x = 12, "
          f"closed form {r_inf:.6f}, 4 V(0)/(pi hbar) = "
          f"{4*V(0.0)/(np.pi*HBAR):.6f}")
    return xs, amps, rows


# --------------------------------------------------------------------- #
# Part B                                                                 #
# --------------------------------------------------------------------- #
def part_b():
    banner("B. Signal against noise as a function of distance")
    print(f"  {'x':>6} {'signal':>12} {'|V.(x)|':>12} {'rel':>9} "
          f"{'rate':>9} {'D':>9} {'S/sqrt(D)':>11} {'defect':>10}")
    rows = []
    for x in [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 6.0, 8.0, 12.0]:
        s, r, d, defect = moments(x)
        rel = abs(s - abs(dV(x))) / max(abs(dV(x)), 1e-300)
        rows.append((x, s, abs(dV(x)), r, d))
        print(f"  {x:6.1f} {s:12.4e} {abs(dV(x)):12.4e} {rel:9.1e} "
              f"{r:9.4f} {d:9.4f} {s/np.sqrt(d):11.3e} {defect:10.1e}")
    print("\n  The first moment reproduces -V'(x) while it is above the")
    print("  quadrature floor (~1e-10); the rate and the variance injection")
    print("  are flat.  Signal/noise decays as fast as the barrier itself.")
    return rows


# --------------------------------------------------------------------- #
# Part C                                                                 #
# --------------------------------------------------------------------- #
def part_c():
    banner("C. The coherence horizon localises the activity")
    xs = [0.0, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0]
    table = {}
    for Lc in [4.0, 8.0, 16.0]:
        print(f"\n  Lc = {Lc:4.1f}   momentum quantum pi hbar / Lc = "
              f"{np.pi*HBAR/Lc:.6f}   horizon edge |x| = "
              f"{Lc/2 + 3*SIG:.2f}")
        print(f"  {'x':>6} {'rate':>12} {'signal':>12} {'|V.(x)|':>12} "
              f"{'defect':>10}")
        rates = []
        for x in xs:
            s_, r, d, defect = moments(x, Lc=Lc)
            rates.append(r)
            print(f"  {x:6.1f} {r:12.4e} {s_:12.4e} {abs(dV(x)):12.4e} "
                  f"{defect:10.1e}")
        table[Lc] = rates
    print("\n  Activity collapses once dist(x, supp V) exceeds Lc/2.")
    print("  (Signal entries outside the active region sit on the xi_max")
    print("  quadrature floor and carry no information.)")
    print("  The norm defect is exactly zero for every Lc: the truncated")
    print("  kernel is still odd in xi, so truncation conserves signed")
    print("  number exactly.  Truncation is not absorption.")

    # C2 -- what the horizon costs where the potential lives
    print("\n  C2. Kernel error inside the active region, x = 0.5:")
    err_rows = []
    for Lc in [1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0]:
        xi_l, w_l = xi_nodes(Lc, 40.0)
        k_l = kernel(0.5, xi_l, Lc=Lc)
        k_r = kernel(0.5, xi_l, Lc=None)      # same nodes, no window
        e = (np.sqrt(np.sum((k_l - k_r) ** 2))
             / np.sqrt(np.sum(k_r ** 2)))
        edge = float(V(Lc / 2 - 0.5))
        err_rows.append((Lc, e, edge))
        print(f"    Lc = {Lc:5.1f}   rel L2 error {e:10.3e}   "
              f"V at the window edge {edge:10.3e}")
    print("    The horizon is free once it clears the support of V: the")
    print("    error tracks the potential left outside the window.")

    # C3 -- the failure mode: a potential that does not decay
    print("\n  C3. A non-decaying V: sharp truncation leaves 1/xi tails.")
    print("      V(x) = cos(2 pi x / 4), sampled at x = 0.7")

    def V_per(x):
        return np.cos(2 * np.pi * x / 4.0)

    exact = 2.0 * abs(np.sin(2 * np.pi * 0.7 / 4.0)) / HBAR
    print(f"    exact budget 2|Gamma_q(x)| = {exact:.4f} "
          "(two deltas at xi = +-pi hbar / a)")
    print(f"    {'xi_max':>8} {'rate (Lc=8)':>14} {'ratio':>9}")
    tails = []
    for xm in [5.0, 10.0, 20.0, 40.0, 80.0]:
        _, r_t, _, _ = moments(0.7, Lc=8.0, xi_max=xm, pot=V_per, dy=2e-3)
        tails.append((xm, r_t, r_t / exact))
        print(f"    {xm:8.1f} {r_t:14.4f} {r_t/exact:9.2f}")
    print("    The truncated budget grows without bound while the exact")
    print("    one is finite.  A sharp window is the wrong tool for a")
    print("    crystal potential -- there the coset invariant of Part D")
    print("    applies and no window is needed at all.")
    return xs, table, err_rows, tails


# --------------------------------------------------------------------- #
# Part D                                                                 #
# --------------------------------------------------------------------- #
def sector_masks(N, n_per):
    """Momentum-index residue classes mod n_per."""
    idx = np.arange(N)
    return [(idx % n_per) == c for c in range(n_per)]


def part_d():
    banner("D. Periodicity of V, not compactness of x: the coset invariant")
    M, N, L = 64, 64, 8.0
    n_per = 4                       # V has period a = L / n_per
    a = L / n_per
    modes = [FourierMode(q=n_per, V_q=0.9, phi_q=0.3),
             FourierMode(q=2 * n_per, V_q=0.4, phi_q=-1.1)]
    dp = np.pi * HBAR / L
    print(f"  L = {L}, V period a = {a}, dp = pi hbar / L = {dp:.6f}")
    print(f"  jump quantum pi hbar / a = {np.pi*HBAR/a:.6f} "
          f"= {n_per} cells -> {n_per} residue classes")

    rng = np.random.default_rng(7)
    W0 = rng.standard_normal((N, M)) * 0.05

    def run(W_init, nsteps=40, dt=0.01):
        lat = PhaseSpaceCrystalLattice(M=M, N=N, L=L, advection="spectral")
        lat.initialize_from_wigner(lambda X, P: np.zeros_like(X))
        lat.W = W_init.copy()
        for _ in range(nsteps):
            lat.strang_step_fourier(modes, dt)
        return lat.W

    masks = sector_masks(N, n_per)

    # D1 invariance: a state living on one class stays on it
    worst_leak = 0.0
    for c, mk in enumerate(masks):
        Wc = np.where(mk[:, None], W0, 0.0)
        Wt = run(Wc)
        leak = np.max(np.abs(Wt[~mk, :]))
        worst_leak = max(worst_leak, leak)
        print(f"  class {c}: leak into the other {n_per-1} classes "
              f"= {leak:.3e}")
    print(f"  worst leak over all classes: {worst_leak:.3e}")

    # D2 completeness: the sum of the sector runs is the full run
    full = run(W0)
    recon = sum(run(np.where(mk[:, None], W0, 0.0)) for mk in masks)
    err = np.max(np.abs(full - recon)) / np.max(np.abs(full))
    print(f"  sum of sector runs vs direct run: max rel diff {err:.3e}")

    # D3 the same statement without any box at all: world-particle momenta
    p0 = 0.3137                                   # arbitrary offset, not on dp
    jumps = np.array([n_per, -n_per, 2 * n_per, -2 * n_per] * 500)
    p = p0 + np.cumsum(jumps) * dp
    theta = np.mod(p, np.pi * HBAR / a)
    print(f"  free world on the open line, 2000 vertices: "
          f"max |theta - theta_0| = "
          f"{np.max(np.abs(theta - np.mod(p0, np.pi*HBAR/a))):.3e}")
    return worst_leak, err


# --------------------------------------------------------------------- #
# Part E                                                                 #
# --------------------------------------------------------------------- #
def part_e():
    banner("E. Grid refinement multiplies sectors, it does not couple them")
    M, N, L = 48, 96, 8.0
    print("  Spec section 1 offers dp = pi hbar / (K L) for finer momentum")
    print("  resolution, with mode q driving jumps of K q cells.  Those K")
    print("  sub-lattices never exchange anything.")
    rows = []
    for K in [1, 2, 3, 4]:
        modes = [FourierMode(q=K, V_q=1.0, phi_q=0.0),
                 FourierMode(q=3 * K, V_q=0.3, phi_q=0.7)]
        rng = np.random.default_rng(11)
        W0 = rng.standard_normal((N, M)) * 0.05
        masks = sector_masks(N, K)
        lat = PhaseSpaceCrystalLattice(M=M, N=N, L=L, advection="spectral")
        lat.initialize_from_wigner(lambda X, P: np.zeros_like(X))
        lat.W = np.where(masks[0][:, None], W0, 0.0)
        for _ in range(30):
            lat.strang_step_fourier(modes, 0.01)
        leak = np.max(np.abs(lat.W[~masks[0], :]), initial=0.0)
        rows.append((K, leak))
        print(f"  K = {K}: {K} sectors, leak out of sector 0 = {leak:.3e}")
    return rows


# --------------------------------------------------------------------- #
# Part F                                                                 #
# --------------------------------------------------------------------- #
def reachable(kappas, gens, window=3.0, work=14.0):
    """Momenta reachable from p = 0 in ``gens`` vertices, reported in
    [-window, window].  The working window is wider so that excursions
    which leave and return are not pruned."""
    steps = np.array(kappas) * HBAR / 2.0
    reach = {0.0}
    for _ in range(gens):
        new = {round(p + sg * d, 12)
               for p in reach for d in steps for sg in (1, -1)}
        reach |= {v for v in new if abs(v) <= work}
    arr = np.array(sorted(v for v in reach if abs(v) <= window))
    return arr


def part_f():
    banner("F. Almost-periodic V: crystal, quasicrystal, continuum")
    print(f"  {'r':>3} {'gens':>5} {'points':>8} {'min gap':>12} "
          f"{'max gap':>12}")
    out = {}
    for kappas in ([1.0], [1.0, np.sqrt(2)], [1.0, np.sqrt(2), np.sqrt(3)]):
        r = len(kappas)
        series = []
        for gens in [2, 4, 6, 8]:
            arr = reachable(kappas, gens)
            g = np.diff(arr)
            series.append((gens, len(arr), g.min(), g.max()))
            print(f"  {r:3d} {gens:5d} {len(arr):8d} {g.min():12.3e} "
                  f"{g.max():12.3e}")
        out[r] = series
        print()
    print("  r = 1 (periodic V): a genuine lattice, gap pinned at hbar k / 2.")
    print("  r > 1: a dense rank-r Z-module.  The invariant p mod Lambda")
    print("  survives set-theoretically but its quotient is no longer")
    print("  Hausdorff, so it stops being a usable quantum number.")
    return out


# --------------------------------------------------------------------- #
# Part G                                                                 #
# --------------------------------------------------------------------- #
def wigner_ring(rho):
    """Discrete Wigner transform on an M-site ring, project convention.

    W[m, n] = (2/M) sum_j rho[m+j, m-j] exp(-2 pi i n j / M),
    which puts momentum on the half-grid p_n = n pi hbar / L.
    """
    M = rho.shape[0]
    j = np.arange(M)
    A = np.empty((M, M), dtype=complex)
    for m in range(M):
        A[m, :] = rho[(m + j) % M, (m - j) % M]
    return (2.0 / M) * np.fft.fft(A, axis=1)


def part_g():
    banner("G. An absorber is the same stencil with the opposite sign rule")
    M, q = 32, 3
    rng = np.random.default_rng(3)
    psi = (rng.standard_normal(M) + 1j * rng.standard_normal(M))
    psi /= np.linalg.norm(psi)
    rho = np.outer(psi, psi.conj())
    x = np.arange(M)
    g = np.cos(2 * np.pi * q * x / M)            # a single cosine mode

    W = wigner_ring(rho)                          # index [midpoint, momentum]
    comm = (g[:, None] - g[None, :]) * rho        # Gamma rho - rho Gamma
    anti = (g[:, None] + g[None, :]) * rho        # Gamma rho + rho Gamma
    Wc, Wa = wigner_ring(comm), wigner_ring(anti)

    # predicted stencils: the same two offset rows, difference vs sum
    lo = np.roll(W, +q, axis=1)                   # W at n - q
    hi = np.roll(W, -q, axis=1)                   # W at n + q
    s = np.sin(2 * np.pi * q * x / M)[:, None]
    c = np.cos(2 * np.pi * q * x / M)[:, None]
    pred_c = 1j * s * (lo - hi)
    pred_a = c * (hi + lo)

    scale = np.max(np.abs(W))
    ec = np.max(np.abs(Wc - pred_c)) / scale
    ea = np.max(np.abs(Wa - pred_a)) / scale
    print(f"  commutator     -> difference of rows n-q, n+q: err {ec:.3e}")
    print(f"  anticommutator -> sum        of rows n-q, n+q: err {ea:.3e}")
    print("\n  So V -> V - i Gamma needs no new lattice geometry, only a")
    print("  sign rule.  But the sum stencil has non-zero column sums, so it")
    print("  destroys signed number by construction, and Gamma's own Wigner")
    print("  image has the same flat envelope proved in Part A: a complex")
    print("  absorbing potential absorbs everywhere, not at a boundary.")
    return ec, ea


# --------------------------------------------------------------------- #
# Figures                                                                #
# --------------------------------------------------------------------- #
def fig_kernel(b_rows, c_xs, c_table, err_rows):
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5))

    ax = axes[0, 0]
    xi = np.linspace(-6, 6, 3001)
    for x, col in zip([1.0, 4.0, 12.0], ["#1f77b4", "#d62728", "#2ca02c"]):
        ax.plot(xi, kernel(x, xi), color=col, lw=1.0, label=f"$x = {x:g}$")
    ax.set_xlabel(r"momentum transfer $\xi$")
    ax.set_ylabel(r"$V_W(x,\xi)$")
    ax.set_title("Constant envelope, finer fringes\n"
                 r"(the potential is supported near $x=0$)", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[0, 1]
    xs = np.array([r[0] for r in b_rows])
    ax.semilogy(xs, [max(r[1], 1e-16) for r in b_rows], "o-",
                label="signal $|\\int \\xi V_W|$")
    ax.semilogy(xs, [max(r[2], 1e-16) for r in b_rows], "k--",
                label="$|V'(x)|$")
    ax.semilogy(xs, [r[3] for r in b_rows], "s-", color="#d62728",
                label=r"rate $\int |V_W|$")
    ax.axhline(4 * V(0.0) / (np.pi * HBAR), color="0.5", ls=":",
               label=r"$4V(0)/\pi\hbar$")
    ax.set_xlabel("$x$")
    ax.set_title("Signal decays, noise does not", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")

    ax = axes[1, 0]
    for Lc, mk in zip(sorted(c_table), ["o-", "s-", "^-"]):
        ax.semilogy(c_xs, np.maximum(c_table[Lc], 1e-100), mk,
                    label=f"$L_c = {Lc:g}$")
        ax.axvline(Lc / 2 + 3 * SIG, color="0.7", ls=":", lw=1)
    ax.set_xlabel("$x$")
    ax.set_ylabel("jump rate")
    ax.set_ylim(1e-40, 1e2)
    ax.set_title("A coherence horizon localises the activity\n"
                 r"(dotted: $L_c/2 + 3\sigma$)", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")

    ax = axes[1, 1]
    ax.semilogy([r[0] for r in err_rows],
                [max(r[1], 1e-18) for r in err_rows], "o-", color="#9467bd",
                label="kernel error at $x = 0.5$")
    ax.semilogy([r[0] for r in err_rows],
                [max(r[2], 1e-18) for r in err_rows], "k--",
                label="$V$ at the window edge")
    ax.set_xlabel("$L_c$")
    ax.set_title("What the horizon costs where $V$ lives", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")

    fig.tight_layout()
    save_fig(fig, "open_position_space_kernel.png")
    plt.close(fig)


def fig_sectors(f_out):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.0))

    ax = axes[0]
    a, dp_a = 2.0, np.pi * HBAR / 2.0
    for c, col in zip(range(3), ["#1f77b4", "#d62728", "#2ca02c"]):
        pts = c * dp_a / 3 + np.arange(-4, 5) * dp_a
        ax.plot(pts, np.full_like(pts, c), "o", color=col, ms=6,
                label=rf"$\theta = {c}/3 \cdot \pi\hbar/a$")
    ax.set_yticks([0, 1, 2])
    ax.set_yticklabels(["", "", ""])
    ax.set_ylim(-0.8, 2.8)
    ax.set_xlabel("$p$")
    ax.set_title("Periodic $V$ on open space: one lattice per coset\n"
                 "$p$ modulo $\\pi\\hbar/a$ is conserved exactly", fontsize=10)
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(alpha=0.3, axis="x")

    ax = axes[1]
    for r, mk in zip(sorted(f_out), ["o-", "s-", "^-"]):
        gens = [s[0] for s in f_out[r]]
        mins = [s[2] for s in f_out[r]]
        ax.semilogy(gens, mins, mk, label=f"$r = {r}$")
    ax.set_xlabel("vertices traversed")
    ax.set_ylabel("smallest gap in the reachable set")
    ax.set_title("Crystal, quasicrystal, continuum", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")

    fig.subplots_adjust(left=0.07, right=0.98, top=0.85, bottom=0.14,
                        wspace=0.28)
    save_fig(fig, "open_position_space_sectors.png")
    plt.close(fig)


def main():
    print("Verification for docs/analysis/open_position_space.md")
    part_a()
    b_rows = part_b()
    c_xs, c_table, err_rows, tails = part_c()
    part_d()
    part_e()
    f_out = part_f()
    part_g()
    banner("Figures")
    fig_kernel(b_rows, c_xs, c_table, err_rows)
    fig_sectors(f_out)
    print("\ndone.")


if __name__ == "__main__":
    main()
