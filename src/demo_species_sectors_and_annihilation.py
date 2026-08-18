"""
Verification for ``docs/analysis/species_sectors_and_annihilation.md``.

What world-particle *species* is, what sectors the excess population divides
into, and what it takes to add an annihilation process to the algorithm.

Conventions.  ``x1`` is the ket leg, ``x2`` the bra leg, ``X = (x1 + x2)/2``
the pair midpoint and ``Y = x1 - x2`` the *full* leg separation, matching
``docs/analysis/interworld_coupling.md``.  The Wigner transform is a Fourier
transform in ``Y``,

    W(X, p) = (1 / 2 pi hbar) Int dY  rho(X + Y/2, X - Y/2) exp(-i p Y / hbar),

evaluated by direct quadrature on an independent ``Y`` grid.  The single-grid
discrete Wigner transform is *not* used: it carries a spurious aliased copy
that gives a coherent state a negativity of 0.5, which is impossible.  Part A
validates the quadrature against the published table of
``docs/supplement/representation_cost_and_annihilation.md`` before anything
else is measured.

Two ensembles appear throughout and must not be conflated:

    E1  carriers drawn from the quasi-density W(x, p); species is the sign of
        W; structure group Z2.  The sea-dressed / four-action / crystal layer.
    E2  carriers drawn from the density matrix rho(x1, x2); positon and negaton
        name the ket and bra legs; phase mu = arg rho; structure group U(1).
        The position-pair / phase-alignment layer.

The Weyl transform relates the represented objects.  Part B shows that it does
not relate the ensembles.

Parts
-----
A  Transform validation, then Theorem D1: species and phase are one degree of
   freedom in conjugate bases.
B  Theorem D0: no carrier-level correspondence between E1 and E2.
C  Theorems D2, D2.1, D2.2: dark under every Hamiltonian is exactly
   c * Identity; the crystal shift is c = 2 and the neutral sea is c = 0.
D  Theorem D3 and Proposition D4: pair-phase precession, and why exactly zero
   separation is special.
E  Theorems D5 and D6: the reach of annihilation, and the soft blob as an
   imposed coherence length.
F  Theorems D8, D9, D10: the excess splits into a diagonal sector and a
   column-balanced sector, and the two substeps conserve complementary
   marginals.
G  Proposition U1 and the Perron root: what annihilation exists to cancel.
H  Theorems D12, D13, D14: species conjugation is momentum reflection at the
   channel level; Z2 is the Hermiticity residue of U(1); species is not
   sign(p).
I  Theorems D15, D16: the four actions split and combine bound sea pairs
   rather than creating them, so positon and negaton number are separately
   conserved, and the sea is a finite reservoir with a sizing floor.
J  Theorem D17 and Proposition D18: adaptive sea allocation, and what it saves.

Run with ``WPMW_OUTPUT`` set (``/mnt/user-data/outputs`` in the container).
"""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from wpmwlib.wpmw_utils import docs_path, output_path

HBAR = 1.0

V_P = 1.5
LAMBDA = 8.0
K_MODE = 2.0 * np.pi / LAMBDA

R = {}


def banner(text):
    print()
    print("=" * 76)
    print(text)
    print("=" * 76)


def save_fig(fig, name):
    fig.savefig(output_path(name), dpi=150, bbox_inches="tight")
    dp = docs_path(name)
    if dp:
        fig.savefig(dp, dpi=150, bbox_inches="tight")
    print(f"  wrote {name}")


# ----------------------------------------------------------------------
# Pair kernel and Wigner transform
# ----------------------------------------------------------------------
def pair_kernel(psi, big_x, n_y=8192, y_max=40.0):
    big_y = (np.arange(n_y) - n_y // 2) * (2.0 * y_max / n_y)
    kern = psi(big_x[:, None] + big_y[None, :] / 2.0) * \
        np.conj(psi(big_x[:, None] - big_y[None, :] / 2.0))
    return kern, big_y


def wigner(kern, big_y):
    d_y = big_y[1] - big_y[0]
    p = np.fft.fftshift(np.fft.fftfreq(len(big_y), d=d_y)) * 2.0 * np.pi * HBAR
    big_w = np.fft.fftshift(
        np.fft.fft(np.fft.ifftshift(kern, axes=1), axis=1), axes=1
    ) * d_y / (2.0 * np.pi * HBAR)
    return big_w, p, p[1] - p[0]


def cat_state(sigma, sep, kick=0.0):
    nrm = np.sqrt(2.0 * np.sqrt(2.0 * np.pi) * sigma
                  * (1.0 + np.exp(-sep ** 2 / (8.0 * sigma ** 2))))

    def psi(z):
        return (np.exp(-(z - sep / 2.0) ** 2 / (4.0 * sigma ** 2))
                + np.exp(-(z + sep / 2.0) ** 2 / (4.0 * sigma ** 2))) \
            * np.exp(1j * kick * z) / nrm
    return psi


def gaussian_state(sigma, kick=0.0):
    def psi(z):
        return np.exp(-z ** 2 / (4.0 * sigma ** 2)) * np.exp(1j * kick * z) \
            / (2.0 * np.pi * sigma ** 2) ** 0.25
    return psi


def census(psi, big_x):
    """(E1 negaton mass fraction, E2 negative-phase weight fraction, ||W||_1)."""
    d_x = big_x[1] - big_x[0]
    kern, big_y = pair_kernel(psi, big_x)
    big_w, _, d_p = wigner(kern, big_y)
    big_w = big_w.real
    l1 = np.abs(big_w).sum() * d_x * d_p
    neg = -big_w[big_w < 0].sum() * d_x * d_p
    amp = np.abs(kern)
    live = amp > 1e-10 * amp.max()
    frac = float((amp[live] * (np.cos(np.angle(kern))[live] < 0)).sum()
                 / amp[live].sum())
    return neg / l1, frac, l1


# ----------------------------------------------------------------------
# A -- validation and the conjugacy theorem
# ----------------------------------------------------------------------
def part_a():
    banner("A  Transform validation, then Theorem D1 (species and phase are conjugate)")
    sigma = 0.5
    big_x = np.linspace(-14.0, 14.0, 2801)
    d_x = big_x[1] - big_x[0]
    print("  A1  quadrature against the published Part A table of")
    print("      docs/supplement/representation_cost_and_annihilation.md")
    print(f"      {'d/sigma':>8} {'Int W':>12} {'||W||_1':>10} {'note':>8}"
          f" {'nu':>9} {'note':>8}")
    for ratio, l1_note, nu_note in [(2, 1.0027, 0.0014), (4, 1.2082, 0.1041),
                                    (8, 1.5875, 0.2937), (16, 1.6366, 0.3183)]:
        kern, big_y = pair_kernel(cat_state(sigma, ratio * sigma), big_x)
        big_w, _, d_p = wigner(kern, big_y)
        big_w = big_w.real
        print(f"      {ratio:>8} {big_w.sum() * d_x * d_p:>12.8f}"
              f" {np.abs(big_w).sum() * d_x * d_p:>10.4f} {l1_note:>8.4f}"
              f" {-big_w[big_w < 0].sum() * d_x * d_p:>9.4f} {nu_note:>8.4f}")

    print("\n  A2  the cat at rest has NO phase anywhere, and 29 per cent negative mass")
    for kick, tag in [(0.0, "at rest, psi real"), (3.0, "boosted by k = 3")]:
        kern, big_y = pair_kernel(cat_state(sigma, 8.0 * sigma, kick), big_x)
        big_w, _, d_p = wigner(kern, big_y)
        big_w = big_w.real
        amp = np.abs(kern)
        live = amp > 1e-10 * amp.max()
        mu = np.angle(kern)
        neg = -big_w[big_w < 0].sum() * d_x * d_p
        frac = float(np.mean((np.cos(mu) < 0)[live]))
        print(f"      {tag:<20} negativity nu = {neg:.4f},"
              f"  fraction of pairs with cos(mu) < 0 = {frac:.4f}")
        if kick:
            yy = np.broadcast_to(big_y[None, :], kern.shape)
            err = np.abs(np.angle(np.exp(1j * (mu - kick * yy))))[live].max()
            print(f"      {'':<20} mu == k Y exactly: max error {err:.3e}")
            R["mu_boost"] = err
        else:
            R["nu_rest"] = neg
    print("""
      => species is not a relabelled pair phase.  It is a functional of the
         whole Y-fibre over a midpoint X, manufactured by exp(-i p Y / hbar).
         Phase lives on Y, sign lives on p, and no carrier holds both sharply.""")


# ----------------------------------------------------------------------
# B -- no carrier correspondence
# ----------------------------------------------------------------------
def part_b():
    banner("B  Theorem D0: the two ensembles have no carrier-level correspondence")
    big_x = np.linspace(-14.0, 14.0, 2801)
    print(f"  {'state':<30} {'E1 negaton mass':>17} {'E2 cos(mu) < 0':>17}")
    rows = []
    for tag, psi in [("cat at rest, d/sigma = 8", cat_state(0.5, 4.0)),
                     ("cat at rest, d/sigma = 16", cat_state(0.5, 8.0)),
                     ("Gaussian, boosted k = 3", gaussian_state(0.5, 3.0)),
                     ("Gaussian, boosted k = 8", gaussian_state(0.5, 8.0)),
                     ("Gaussian at rest", gaussian_state(0.5))]:
        e_one, e_two, _ = census(psi, big_x)
        print(f"  {tag:<30} {e_one:>17.4f} {e_two:>17.4f}")
        rows.append((tag, e_one, e_two))
    R["census"] = rows
    print("""
      => the two censuses are anti-correlated across these states.  A map
         carrying E1 carriers to E2 carriers one at a time would have to send a
         homogeneous population to a heterogeneous one and back.  The Weyl
         transform relates the represented OBJECTS; it does not relate the
         ENSEMBLES.  "Positon" in the two layers is a homonym.""")


# ----------------------------------------------------------------------
# C -- what is dark, and which sea is which
# ----------------------------------------------------------------------
def part_c():
    banner("C  Theorems D2 / D2.1 / D2.2: dark under every H is exactly c * Identity")
    n = 64
    rng = np.random.default_rng(20260816)

    def rand_h():
        a = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
        return a + a.conj().T

    hams = [rand_h() for _ in range(12)]
    print(f"  {'operator':<50} {'max ||[H, X]||':>16}")
    for name, mat in [
        ("0  (neutral sea: bound pairs, E = 0)", np.zeros((n, n), complex)),
        ("2 * Identity  (crystal shift, W -> W + 2/h)", 2.0 * np.eye(n)),
        ("c * Identity, c = -0.37", -0.37 * np.eye(n)),
        ("diagonal, non-constant", np.diag(rng.normal(size=n)).astype(complex)),
        ("rank-one projector", np.outer(np.ones(n), np.ones(n)) / n),
        ("random Hermitian", rand_h()),
    ]:
        worst = max(np.abs(h @ mat - mat @ h).max() for h in hams)
        print(f"  {name:<50} {worst:>16.3e}")
        if name.startswith("2 *"):
            R["dark_shift"] = worst
    print("""
      => the commutant of the full matrix algebra is the scalars, so an
         operator is dark under every Hamiltonian IFF it is c * Identity.  Both
         sea readings live in that one-parameter family and differ only in c:

           c = 0  the neutral sea, one positon bound to one negaton at the same
                  phase-space cell.  Invisible in the observable (E = 0) and
                  LIVE in the dynamics: the reservoir the four actions draw on,
                  and what broken detailed balance pumps.
           c = 2  the crystal shift, W -> W + 2/h.  Visible in the observable
                  and provably INERT -- which is why it cannot be the medium.""")

    print("  C2  the Weyl symbol of a background of coherence length eps")
    print(f"      {'eps':>8} {'W at p = 0':>14} {'2/h':>12} {'ratio':>9}"
          f" {'plateau half-width':>20}")
    two_over_h = 2.0 / (2.0 * np.pi * HBAR)
    for eps in [2.0, 0.5, 0.1, 0.02]:
        n_y, y_max = 1 << 16, 80.0
        big_y = (np.arange(n_y) - n_y // 2) * (2.0 * y_max / n_y)
        kern = np.array(np.broadcast_to(
            (2.0 * np.exp(-big_y ** 2 / (2.0 * eps ** 2))
             / (eps * np.sqrt(2.0 * np.pi)))[None, :], (2, n_y)))
        big_w, p, _ = wigner(kern, big_y)
        big_w = big_w.real
        i0 = int(np.argmin(np.abs(p)))
        half = np.abs(p[np.abs(big_w[0]) > 0.5 * np.abs(big_w[0]).max()]).max()
        print(f"      {eps:>8.3f} {big_w[0, i0]:>14.8f} {two_over_h:>12.8f}"
              f" {big_w[0, i0] / two_over_h:>9.5f} {half:>20.4f}")
    print("      => in the pair basis the identity is delta(Y): the shift sits at")
    print("         exactly zero leg separation, which is why its darkness depends")
    print("         only on relative position.")


# ----------------------------------------------------------------------
# D -- precession
# ----------------------------------------------------------------------
def part_d():
    banner("D  Theorem D3 and Proposition D4: dmu/dt = -2 sum_q Gamma_q(X) sin(k_q Y/2)")
    rng = np.random.default_rng(20260816)
    for modes in ([(1, V_P, 0.0)],
                  [(1, V_P, 0.3), (2, -0.7, 1.1), (5, 0.25, -2.0)]):
        def v_fun(z, mm=modes):
            return sum(a * np.cos(2.0 * np.pi * qq * z / LAMBDA + ph)
                       for qq, a, ph in mm)
        big_x = rng.uniform(-LAMBDA / 2, LAMBDA / 2, 60000)
        big_y = rng.uniform(-2.0 * LAMBDA, 2.0 * LAMBDA, 60000)
        lhs = -(v_fun(big_x + big_y / 2.0) - v_fun(big_x - big_y / 2.0)) / HBAR
        rhs = np.zeros_like(lhs)
        for qq, a, ph in modes:
            k = 2.0 * np.pi * qq / LAMBDA
            rhs += -2.0 * (-(a / HBAR) * np.sin(k * big_x + ph)) \
                * np.sin(k * big_y / 2.0)
        err = np.abs(lhs - rhs).max()
        print(f"  {len(modes)} mode(s): max|lhs - rhs| = {err:.3e}"
              f"   (scale {np.abs(lhs).max():.3f})")
        R[f"prec{len(modes)}"] = err

    print("\n  D2  a background of FINITE coherence length is dynamically active")
    print("      A(eps) = max_X Int du |U(X,u)| rho_eps(u), against 2|V'|sqrt(2/pi) eps")
    print(f"      {'eps':>8} {'A(eps)':>14} {'predicted':>14} {'ratio':>9}")
    u = np.linspace(-12.0, 12.0, 120001)
    d_u = u[1] - u[0]
    xs = np.linspace(-LAMBDA / 2, LAMBDA / 2, 401)
    coup = -V_P * (np.cos(K_MODE * (xs[:, None] + u[None, :] / 2))
                   - np.cos(K_MODE * (xs[:, None] - u[None, :] / 2)))
    for eps in [1.0, 0.25, 0.0625]:
        rho_e = 2.0 * np.exp(-u ** 2 / (2.0 * eps ** 2)) \
            / (eps * np.sqrt(2.0 * np.pi))
        act = float((np.abs(coup) * rho_e[None, :]).sum(axis=1).max() * d_u)
        pred = 2.0 * V_P * K_MODE * np.sqrt(2.0 / np.pi) * eps
        print(f"      {eps:>8.4f} {act:>14.6f} {pred:>14.6f} {act / pred:>9.5f}")
        R["activity"] = act / pred
    print("""
      => Y = 0 is dark for every X and every mode; precession is maximal at
         Y = lambda_q / 2 at rate 2|Gamma_q(X)|, whose maximum 2 V_q / hbar is
         the gamma_max of the annihilation burden.  So the pathwise L1 growth
         rate is a phase precession rate.  Only exactly zero separation is dark.""")


# ----------------------------------------------------------------------
# E -- reach, and the soft blob
# ----------------------------------------------------------------------
def part_e():
    banner("E  Theorems D5 and D6: the reach of annihilation, and the soft blob")
    sigma, sep = 0.5, 4.0
    big_x = np.linspace(-10.0, 10.0, 2001)
    d_x = big_x[1] - big_x[0]
    kern, big_y = pair_kernel(cat_state(sigma, sep), big_x)
    big_w, p, d_p = wigner(kern, big_y)
    big_w = big_w.real
    d_y = big_y[1] - big_y[0]
    weight = np.abs(big_w)
    print("  E1  truncate the Y-fibre at |Y| <= Yc and re-read the sign of W")
    print(f"      {'Yc':>8} {'|W|-weighted sign agreement':>30} {'||W_trunc||_1':>15}")
    rows = []
    for yc in [0.5, 2.0, 4.0, 6.0, 8.0, 12.0]:
        cut = kern * (np.abs(big_y) <= yc)[None, :]
        w_t = np.fft.fftshift(
            np.fft.fft(np.fft.ifftshift(cut, axes=1), axis=1), axes=1
        ).real * d_y / (2.0 * np.pi * HBAR)
        agree = float((weight * (np.sign(w_t) == np.sign(big_w))).sum()
                      / weight.sum())
        print(f"      {yc:>8.2f} {agree:>30.4f}"
              f" {np.abs(w_t).sum() * d_x * d_p:>15.4f}")
        rows.append((yc, agree))
    R["reach"] = rows
    print("\n      settling reach against the coherence length d + 4 sigma:")
    print(f"      {'d/sigma':>8} {'d':>7} {'settling Yc':>13} {'d + 4 sigma':>13}")
    for ratio in [4, 8, 16]:
        k2, y2 = pair_kernel(cat_state(sigma, ratio * sigma), big_x)
        w2, _, _ = wigner(k2, y2)
        w2 = w2.real
        wt2 = np.abs(w2)
        dy2 = y2[1] - y2[0]
        found = np.nan
        for yc in np.arange(0.25, 24.0, 0.25):
            cut = k2 * (np.abs(y2) <= yc)[None, :]
            wtr = np.fft.fftshift(
                np.fft.fft(np.fft.ifftshift(cut, axes=1), axis=1), axes=1
            ).real * dy2 / (2.0 * np.pi * HBAR)
            if (wt2 * (np.sign(wtr) == np.sign(w2))).sum() / wt2.sum() > 0.999:
                found = yc
                break
        print(f"      {ratio:>8} {ratio * sigma:>7.2f} {found:>13.2f}"
              f" {ratio * sigma + 4 * sigma:>13.2f}")

    print("\n  E2  a momentum blob of width sigma_p IS a coherence length hbar/sigma_p")
    print(f"      {'sigma_p':>9} {'hbar/sigma_p':>14} {'max discrepancy':>18}"
          f" {'relative':>10} {'||W||_1':>10}")
    l1_rows = []
    for s_p in [0.2, 0.5, 1.0, 2.0]:
        ker = np.exp(-p ** 2 / (2.0 * s_p ** 2))
        ker /= ker.sum()
        conv = np.real(np.fft.ifft(
            np.fft.fft(big_w, axis=1) * np.fft.fft(np.fft.ifftshift(ker)), axis=1))
        win = kern * np.exp(-s_p ** 2 * big_y ** 2 / (2.0 * HBAR ** 2))[None, :]
        w_w = np.fft.fftshift(
            np.fft.fft(np.fft.ifftshift(win, axes=1), axis=1), axes=1
        ).real * d_y / (2.0 * np.pi * HBAR)
        err = np.abs(conv - w_w).max()
        l1 = np.abs(w_w).sum() * d_x * d_p
        print(f"      {s_p:>9.2f} {HBAR / s_p:>14.4f} {err:>18.3e}"
              f" {err / np.abs(w_w).max():>10.2e} {l1:>10.4f}")
        l1_rows.append((s_p, l1))
        R["blob"] = err
    print(f"      {'(none)':>9} {'inf':>14} {'-':>18} {'-':>10}"
          f" {np.abs(big_w).sum() * d_x * d_p:>10.4f}")
    R["blob_l1"] = l1_rows
    print("""
      => the reach is set by the STATE's coherence length, not by the lattice
         spacing; and a soft annihilation kernel is exactly an imposed
         decoherence of known length.  Cell-exact annihilation is local in E1;
         only the soft variant is genuinely non-local in effect.""")


# ----------------------------------------------------------------------
# F -- the sectors
# ----------------------------------------------------------------------
def part_f():
    banner("F  Theorems D8 / D9 / D10: the excess splits into two sectors")
    big_x = np.linspace(-9.0, 9.0, 1801)
    psi = cat_state(0.5, 4.0)
    kern, big_y = pair_kernel(psi, big_x)
    big_w, p, d_p = wigner(kern, big_y)
    big_w = big_w.real
    i0 = int(np.argmin(np.abs(big_y)))
    marg = big_w.sum(axis=1) * d_p
    diag = np.abs(psi(big_x)) ** 2
    print("  F1  the momentum marginal of W is the DIAGONAL of rho")
    print(f"      max| Int W dp - |psi|^2 |  = {np.abs(marg - diag).max():.3e}")
    print(f"      max| Int W dp - C(X, 0) |  = "
          f"{np.abs(marg - kern[:, i0].real).max():.3e}")
    print(f"      min over X of Int W dp     = {marg.min():.3e}   <- never negative")
    print(f"      min over (X, p) of W       = {big_w.min():.6f}"
          f"   <- negative pointwise")
    R["marg"] = np.abs(marg - diag).max()

    print("\n  F2  negatons live entirely in the off-diagonal, which is column-balanced")
    d_y = big_y[1] - big_y[0]
    for tag, win in [("diagonal only (narrow Y window)",
                      np.exp(-big_y ** 2 / (2.0 * 0.02 ** 2))),
                     ("full rho", np.ones_like(big_y))]:
        cut = kern * win[None, :]
        w_c = np.fft.fftshift(
            np.fft.fft(np.fft.ifftshift(cut, axes=1), axis=1), axes=1
        ).real * d_y / (2.0 * np.pi * HBAR)
        neg = -w_c[w_c < 0].sum() * (big_x[1] - big_x[0]) * d_p
        print(f"      {tag:<34} negativity nu = {neg:.4e}")
    off = kern.copy()
    off[:, i0] = 0.0
    w_o = np.fft.fftshift(
        np.fft.fft(np.fft.ifftshift(off, axes=1), axis=1), axes=1
    ).real * d_y / (2.0 * np.pi * HBAR)
    bal = np.abs(w_o.sum(axis=1) * d_p).max()
    print(f"      off-diagonal part: max |column sum| = {bal:.3e}"
          f"   <- exactly balanced")
    R["balanced"] = bal

    print("\n  F3  the two substeps conserve complementary marginals")
    n_x, n_p, qq = 64, 64, 1
    d_xl = LAMBDA / n_x
    xg = (np.arange(n_x) + 0.5) * d_xl - LAMBDA / 2
    ng = np.arange(n_p) - n_p // 2
    gam = -(V_P / HBAR) * np.sin(2.0 * np.pi * qq * xg / LAMBDA)
    rng = np.random.default_rng(3)
    big_e = rng.normal(size=(n_p, n_x))
    col0, row0 = big_e.sum(axis=0).copy(), big_e.sum(axis=1).copy()
    kk = np.fft.fftfreq(n_p, d=1.0 / n_p)
    sym = 1j * 2.0 * np.sin(2.0 * np.pi * qq * kk / n_p)
    e_j = np.real(np.fft.ifft(
        np.fft.fft(big_e, axis=0) * np.exp(0.137 * np.outer(sym, gam)), axis=0))
    e_s = np.empty_like(big_e)
    for i in range(n_p):
        e_s[i] = np.roll(big_e[i], int(ng[i]) % n_x)
    print(f"      {'substep':<26} {'max |d column sum|':>20} {'max |d row sum|':>18}")
    print(f"      {'jump (exact, FFT)':<26}"
          f" {np.abs(e_j.sum(axis=0) - col0).max():>20.3e}"
          f" {np.abs(e_j.sum(axis=1) - row0).max():>18.3e}")
    print(f"      {'streaming (integer)':<26}"
          f" {np.abs(e_s.sum(axis=0) - col0).max():>20.3e}"
          f" {np.abs(e_s.sum(axis=1) - row0).max():>18.3e}")
    R["colcons"] = np.abs(e_j.sum(axis=0) - col0).max()
    print("""
      => the jump substep conserves every position-column sum and moves
         momentum: it is the force.  Streaming conserves every momentum-row sum
         and moves position: it is the transport.  2(M + N) exact invariants,
         free, as per-substep diagnostics.  The unpaired excess positons carry
         the Born density and the jump substep never touches them; all of the
         ledger growth is in the column-balanced, zero-probability sector.""")


# ----------------------------------------------------------------------
# G -- what annihilation is for
# ----------------------------------------------------------------------
def part_g():
    banner("G  Proposition U1 and the Perron root: what annihilation exists to cancel")
    n_p, qq, gamma = 16, 3, -1.5
    lmat = np.zeros((n_p, n_p))
    for c in range(n_p):
        lmat[(c - qq) % n_p, c] += gamma
        lmat[(c + qq) % n_p, c] += -gamma
    absl = np.abs(lmat)
    print(f"  generator diagonal: max |L_cc| = {np.abs(np.diag(lmat)).max():.3e}")
    print(f"  per-column growth L_cc + sum_r |L_rc| = {absl.sum(axis=0)[0]:.6f}")
    rng = np.random.default_rng(20260816)
    print(f"\n  {'guiding function phi':<28} {'rho(|D^-1 L D|)':>18}")
    for tag, phi in [
        ("uniform", np.ones(n_p)),
        ("localised", np.exp(-((np.arange(n_p) - n_p / 2) ** 2) / 8) + 1e-3),
        ("random lognormal", np.exp(rng.normal(0, 1.0, n_p))),
        ("random lognormal x2", np.exp(rng.normal(0, 2.0, n_p))),
    ]:
        d = np.diag(phi)
        lg = np.abs(np.linalg.inv(d) @ lmat @ d)
        print(f"  {tag:<28} {max(abs(np.linalg.eigvals(lg))):>18.10f}")

    l_box, n_x, n_pp, q2, v_q = 8.0, 256, 128, 1, 1.5
    d_x = l_box / n_x
    d_p = np.pi * HBAR / l_box
    x = (np.arange(n_x) + 0.5) * d_x - l_box / 2
    n = np.arange(n_pp) - n_pp // 2
    gam = np.abs(-(v_q / HBAR) * np.sin(2.0 * np.pi * q2 * x / l_box))
    d_t = d_x / d_p
    k = np.fft.fftfreq(n_pp, d=1.0 / n_pp)
    amp = np.exp(d_t * gam[None, :]
                 * (2.0 * np.cos(2.0 * np.pi * q2 * k / n_pp))[:, None])
    shifts = n % n_x
    a = np.ones((n_pp, n_x))
    rates = []
    for _ in range(3000):
        a = np.real(np.fft.ifft(np.fft.fft(a, axis=0) * amp, axis=0))
        out = np.empty_like(a)
        for i in range(n_pp):
            out[i] = np.roll(a[i], shifts[i])
        a = out
        tot = a.sum()
        rates.append(np.log(tot) / d_t)
        a /= tot
    rho_abs = float(np.mean(rates[-400:]))
    print(f"\n  gamma_avg = 4 V_q / (pi hbar)           = "
          f"{4 * v_q / (np.pi * HBAR):.6f}")
    print(f"  rho(|L|), full operator with streaming  = {rho_abs:.6f}"
          f"   (+/- {np.std(rates[-400:]):.1e})")
    print(f"  gamma_max = 2 V_q / hbar                = {2 * v_q / HBAR:.6f}")
    R["rho"] = rho_abs
    print("""
      => the L1 growth rate is a similarity invariant of the generator, so no
         choice of unraveling and no guiding function changes it.  The four
         actions cannot be L1-stationary in their own right; the only remaining
         lever is a channel quadratic in the ensemble.""")
    return rho_abs


# ----------------------------------------------------------------------
# H -- chirality
# ----------------------------------------------------------------------
def part_h():
    banner("H  Theorems D12 / D13 / D14: species is an orientation, not a handedness")
    n, qq = 17, 3

    def emission(flip=False):
        mat = np.zeros((n, n))
        s = -1.0 if flip else 1.0
        for c in range(n):
            mat[(c - qq) % n, c] += s
            mat[(c + qq) % n, c] += -s
        return mat

    refl = np.zeros((n, n))
    for c in range(n):
        refl[(-c) % n, c] = 1.0
    chan, conj = emission(), emission(flip=True)
    print(f"  || R C R - Cbar || = {np.abs(refl @ chan @ refl - conj).max():.3e}"
          f"   <- reflection IS species conjugation")
    print(f"  || R C R - C    || = {np.abs(refl @ chan @ refl - chan).max():.3e}"
          f"   <- and not the identity")
    R["chir"] = np.abs(refl @ chan @ refl - conj).max()

    print("\n  H2  Hermiticity is what cuts U(1) down to Z2")
    big_x = np.linspace(-9.0, 9.0, 1201)
    kern, big_y = pair_kernel(cat_state(0.5, 4.0), big_x)
    for tag, fac in [("Hermitian rho", np.ones_like(big_y)),
                     ("even-in-Y twist (breaks Hermiticity)",
                      np.exp(0.6j * np.exp(-big_y ** 2)))]:
        k2 = kern * fac[None, :]
        herm = np.abs(k2 - np.conj(k2[:, ::-1])).max()
        w2, _, _ = wigner(k2, big_y)
        live = np.abs(w2) > 1e-5 * np.abs(w2).max()
        nvals = np.unique(np.round(np.angle(w2[live]), 4)).size
        print(f"      {tag:<38} residual {herm:.2e},"
              f" |Im W|/|Re W| {np.abs(w2.imag).max() / np.abs(w2.real).max():.2e},"
              f" phases {nvals}")

    print("\n  H3  species is NOT the carrier's own handedness sign(p)")
    print(f"      {'state':<26} {'fraction at p > 0':>19} {'negaton mass':>14}")
    for kick in [0.0, 8.0, 20.0]:
        psi = cat_state(0.5, 4.0, kick)
        k3, y3 = pair_kernel(psi, big_x)
        w3, p3, _ = wigner(k3, y3)
        w3 = w3.real
        wt = np.abs(w3)
        fpos = float((wt * (np.broadcast_to(p3[None, :], w3.shape) > 0)).sum()
                     / wt.sum())
        e_one, _, _ = census(psi, big_x)
        print(f"      cat, boost k = {kick:<11.0f} {fpos:>19.4f} {e_one:>14.4f}")
    print("""
      => species conjugation is exactly momentum reflection AT THE CHANNEL
         level, so the label is the orientation of the momentum transfer rather
         than an independent charge; and Z2 is the residue of U(1) left by
         Hermiticity.  But a boost puts every carrier at positive momentum
         without changing the negaton census at all, so the orientation is
         relational -- a property of an interaction, never locally readable off
         a single carrier.""")


# ----------------------------------------------------------------------
# I and J -- split versus create, sizing, and adaptive allocation
# ----------------------------------------------------------------------
L_BOX, V_Q, Q, T_END, D_T, KAPPA = 8.0, 1.5, 1, 1.0, 0.004, 200.0
SWEEP_NU, SWEEP_STEPS = 2.0e5, 80
N_STEPS = int(T_END / D_T)


def lattice_setup(nu, m_x, n_p):
    d_x = L_BOX / m_x
    d_p = np.pi * HBAR / L_BOX
    x = (np.arange(m_x) + 0.5) * d_x - L_BOX / 2
    n = np.arange(n_p) - n_p // 2
    gam = -(V_Q / HBAR) * np.sin(2.0 * np.pi * Q * x / L_BOX)
    w = np.exp(-(x[None, :] ** 2) / 2.0) \
        * np.exp(-((n[:, None] * d_p) ** 2) / (2.0 * 0.8 ** 2))
    w /= w.sum()
    return gam, np.round(nu * w).astype(np.int64), nu / m_x


def mean_field(e0, gam, n_steps=None):
    n_steps = N_STEPS if n_steps is None else n_steps
    e = e0.astype(float)
    for _ in range(n_steps):
        e = e + D_T * gam[None, :] * (np.roll(e, -Q, axis=0)
                                      - np.roll(e, +Q, axis=0))
    return e


def lattice_run(nu, m_x, n_p, mode, beta=1.0, safety=100.0, seed=0,
                recombine=True, n_steps=None):
    """mode: 'uniform' fixed sea B = beta * B_Wigner; 'adaptive' per-cell
    sizing; 'create' the spawn realisation, which carries no sea at all."""
    n_steps = N_STEPS if n_steps is None else n_steps
    gam, e0, b_w = lattice_setup(nu, m_x, n_p)
    e_ex = mean_field(e0, gam, n_steps)
    rng = np.random.default_rng(seed)
    u_p, u_m = e0.copy(), np.zeros_like(e0)
    if mode == "uniform":
        sea = np.full(e0.shape, int(beta * b_w), dtype=np.int64)
    else:
        sea = np.zeros_like(e0)
    g = np.abs(gam)[None, :]
    sgn = np.sign(gam).astype(int)
    blocked = asked = peak = 0
    ledger, sea_hist = [], []
    for _ in range(n_steps):
        if mode == "adaptive":
            tgt = np.ceil(safety * g * D_T * (u_p + u_m)).astype(np.int64) + 4
            sea += np.maximum(tgt - sea, 0)
            sea -= np.maximum(sea - 2 * tgt, 0)
        f_p = rng.poisson(g * D_T * u_p)
        f_m = rng.poisson(g * D_T * u_m)
        if mode != "create":
            need = f_p + f_m
            asked += int(need.sum())
            over = need > sea
            if over.any():
                blocked += int((need - sea)[over].sum())
                sc = np.where(over, sea / np.maximum(need, 1), 1.0)
                f_p = np.floor(f_p * sc).astype(np.int64)
                f_m = np.floor(f_m * sc).astype(np.int64)
            sea -= f_p + f_m
        pos = (sgn > 0)[None, :]
        up_r, dn_r = np.roll(f_p, -Q, axis=0), np.roll(f_p, +Q, axis=0)
        um_r, dm_r = np.roll(f_m, +Q, axis=0), np.roll(f_m, -Q, axis=0)
        u_p += np.where(pos, up_r + um_r, dn_r + dm_r)
        u_m += np.where(pos, dn_r + dm_r, up_r + um_r)
        if recombine:
            r = np.minimum(rng.poisson(KAPPA * D_T * u_p * u_m
                                       / np.maximum(u_p + u_m, 1)),
                           np.minimum(u_p, u_m))
            u_p -= r
            u_m -= r
            if mode != "create":
                sea += r
        ledger.append(int(u_p.sum() + u_m.sum()))
        sea_hist.append(int(2 * sea.sum()))
        peak = max(peak, ledger[-1] + sea_hist[-1])
    e = (u_p - u_m).astype(float)
    return dict(err=float(np.linalg.norm(e - e_ex) / np.linalg.norm(e_ex)),
                blocked=blocked / max(asked, 1), peak=peak,
                ledger=np.array(ledger), sea=np.array(sea_hist),
                ntot=ledger[-1] + sea_hist[-1], b_w=b_w)


def part_i():
    banner("I  Theorems D15 / D16: the four actions SPLIT and COMBINE, they do not create")
    nu, m_x, n_p = 2.0e6, 32, 48
    print("  I1  what each reading conserves.  N_+ = S + U+ and N_- = S + U- are")
    print("      the positon and negaton censuses; N_total = U+ + U- + 2S.")
    print("      (recombination is switched off here so the ledger actually grows)")
    print(f"      {'reading':<10} {'ledger, start -> end':>28}"
          f" {'2S, start -> end':>26} {'N_total drift':>15}")
    out = {}
    for mode, tag in [("uniform", "split"), ("create", "spawn")]:
        res = lattice_run(nu, m_x, n_p, mode, beta=1.0, seed=0, recombine=False)
        n0 = res["ledger"][0] + res["sea"][0]
        drift = abs(res["ntot"] / n0 - 1.0)
        print(f"      {tag:<10} {res['ledger'][0]:>12.4e} ->"
              f" {res['ledger'][-1]:>10.4e} {res['sea'][0]:>12.4e} ->"
              f" {res['sea'][-1]:>10.4e} {drift:>15.3e}")
        out[mode] = res
        R[f"drift_{mode}"] = drift
    print("""      => under splitting N_total is conserved to floating point; combined
         with conservation of sum E that forces N_+ and N_- to be separately,
         exactly conserved.  Species is a conserved charge that is bound or
         unbound, never created.  The spawn reading conserves neither.""")

    print("\n  I2  the sizing floor, in units of the Wigner capacity"
          " B_W = nu (2/h) dx dp")
    print(f"      {'Mx':>5} {'Np':>5} {'B_W':>10}"
          f"   blocked fraction at beta = B / B_W:")
    print(f"      {'':>5} {'':>5} {'':>10} "
          + "".join(f"{b:>9.3f}" for b in [0.01, 0.03, 0.1, 0.3, 1.0]))
    rows = []
    for mx, npp in [(32, 48), (64, 48), (32, 96), (64, 96)]:
        vals = [lattice_run(SWEEP_NU, mx, npp, "uniform", beta=b, seed=0,
                            n_steps=SWEEP_STEPS)["blocked"]
                for b in [0.01, 0.03, 0.1, 0.3, 1.0]]
        b_w = SWEEP_NU / mx
        print(f"      {mx:>5} {npp:>5} {b_w:>10.0f} "
              + "".join(f"{v:>9.4f}" for v in vals))
        rows.append((mx, npp, vals))
    R["floor"] = rows
    print("""      => the threshold is universal in units of B_W: identical across a
         two-fold range of Mx and of Np.  So the Wigner capacity is a safe,
         lattice-independent sufficient bound -- and loose by about a factor of
         twenty, which is what the adaptive scheme recovers.""")
    return out


def part_j(rho_abs):
    banner("J  Theorem D17 and Proposition D18: adaptive sea allocation")
    nu = 2.0e6
    print("  injecting or removing a BOUND pair changes E by zero and every moment")
    print("  by zero -- the same argument as same-cell annihilation -- so the sea")
    print("  size is a representational choice, free to vary per cell and per step.\n")
    print(f"  {'scheme':<28} {'blocked':>9} {'rel L2 err':>12}"
          f" {'peak N_total':>14} {'vs uniform':>11}")
    curves = {}
    for mx, npp in [(32, 48), (64, 96)]:
        print(f"  --- Mx = {mx}, Np = {npp} ---")
        res = lattice_run(nu, mx, npp, "uniform", beta=1.0, seed=0)
        base = res["peak"]
        curves[(mx, npp, "uniform")] = res
        print(f"  {'uniform B = B_Wigner':<28} {res['blocked']:>9.4f}"
              f" {res['err']:>12.5f} {res['peak']:>14.3e} {1.0:>11.3f}")
        for saf in [30.0, 100.0, 300.0]:
            r2 = lattice_run(nu, mx, npp, "adaptive", safety=saf, seed=0)
            print(f"  {'adaptive, safety = ' + str(int(saf)):<28}"
                  f" {r2['blocked']:>9.4f} {r2['err']:>12.5f}"
                  f" {r2['peak']:>14.3e} {r2['peak'] / base:>11.3f}")
            if saf == 100.0:
                curves[(mx, npp, "adaptive")] = r2
                R[f"save_{mx}_{npp}"] = r2["peak"] / base
    print("""
      => identical error, zero blocked splits, and one to two orders of
         magnitude fewer world-particles.  A uniform sea costs O(N_p * nu),
         because it stocks every cell against the global Wigner bound; the
         adaptive sea costs O(ledger), because it stocks each cell against that
         cell's own split flux.  The saving is of order the number of momentum
         cells.""")

    print("\n  J2  with the recombination channel switched off")
    on = lattice_run(nu, 32, 48, "adaptive", safety=100.0, seed=0)
    off = lattice_run(nu, 32, 48, "adaptive", safety=100.0, seed=0,
                      recombine=False)
    print(f"      {'recombination':<16} {'final ledger':>14} {'peak N_total':>14}"
          f" {'rel L2 err':>12}")
    for tag, res in [("on", on), ("off", off)]:
        print(f"      {tag:<16} {res['ledger'][-1]:>14.4e} {res['peak']:>14.3e}"
              f" {res['err']:>12.5f}")
    print(f"      ungoverned growth asymptotes to rho(|L|) = {rho_abs:.3f}")
    R["rec"] = (on["ledger"][-1], off["ledger"][-1])
    curves["rec_on"], curves["rec_off"] = on, off
    return curves


# ----------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------
def fig_duality():
    sigma, sep = 0.5, 4.0
    big_x = np.linspace(-6.0, 6.0, 601)
    kern, big_y = pair_kernel(cat_state(sigma, sep), big_x, n_y=4096, y_max=24.0)
    big_w, p, _ = wigner(kern, big_y)
    big_w = big_w.real
    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.2))

    ax = axes[0]
    sel = np.abs(p) < 4.0
    lim = np.abs(big_w).max()
    im = ax.pcolormesh(big_x, p[sel], big_w[:, sel].T, cmap="RdBu_r",
                       vmin=-lim, vmax=lim, shading="auto")
    ax.set_xlabel(r"midpoint $X$")
    ax.set_ylabel(r"momentum $p$")
    ax.set_title(r"(a) $W(X,p)$: sign structure runs along $p$", fontsize=10.5)
    fig.colorbar(im, ax=ax, fraction=0.046)

    ax = axes[1]
    ysel = np.abs(big_y) < 8.0
    im = ax.pcolormesh(big_x, big_y[ysel], np.abs(kern[:, ysel]).T,
                       cmap="magma", shading="auto")
    ax.set_xlabel(r"midpoint $X$")
    ax.set_ylabel(r"leg separation $Y$")
    ax.set_title(r"(b) $|\rho|$ on the same fibres; $\mu = \arg\rho \equiv 0$",
                 fontsize=10.5)
    ax.text(0.03, 0.95, "no phase anywhere,\n" + r"yet $\nu = 0.29$",
            transform=ax.transAxes, fontsize=9, color="w", va="top")
    fig.colorbar(im, ax=ax, fraction=0.046)

    ax = axes[2]
    d_y = big_y[1] - big_y[0]
    weight = np.abs(big_w)
    reaches = np.arange(0.25, 14.0, 0.25)
    agree = []
    for yc in reaches:
        cut = kern * (np.abs(big_y) <= yc)[None, :]
        w_t = np.fft.fftshift(
            np.fft.fft(np.fft.ifftshift(cut, axes=1), axis=1), axes=1
        ).real * d_y / (2.0 * np.pi * HBAR)
        agree.append((weight * (np.sign(w_t) == np.sign(big_w))).sum()
                     / weight.sum())
    ax.plot(reaches, agree, lw=1.8, color="#1f6feb")
    ax.axvline(sep, color="crimson", ls="--", lw=1.2,
               label="coherence length $d$")
    ax.set_xlabel(r"leg reach $Y_c$")
    ax.set_ylabel("sign agreement with full $W$")
    ax.set_ylim(0.5, 1.02)
    ax.set_title("(c) how far the fibre must be read", fontsize=10.5)
    ax.legend(fontsize=8.5, loc="lower right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    save_fig(fig, "species_phase_duality.png")
    plt.close(fig)


def fig_sea():
    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.2))
    two_over_h = 2.0 / (2.0 * np.pi * HBAR)

    ax = axes[0]
    n_y, y_max = 1 << 15, 60.0
    big_y = (np.arange(n_y) - n_y // 2) * (2.0 * y_max / n_y)
    for eps, col in zip([2.0, 1.0, 0.5, 0.25],
                        ["#c2d6f0", "#7aa7dd", "#3d78c2", "#0b3d91"]):
        kern = np.array(np.broadcast_to(
            (2.0 * np.exp(-big_y ** 2 / (2.0 * eps ** 2))
             / (eps * np.sqrt(2.0 * np.pi)))[None, :], (2, n_y)))
        big_w, p, _ = wigner(kern, big_y)
        sel = np.abs(p) < 8.0
        ax.plot(p[sel], big_w.real[0, sel], color=col, lw=1.6,
                label=rf"$\epsilon = {eps}$")
    ax.axhline(two_over_h, color="crimson", ls="--", lw=1.2, label=r"$2/h$")
    ax.set_xlabel(r"momentum $p$")
    ax.set_ylabel("Weyl symbol")
    ax.set_title(r"(a) crystal shift: $\epsilon \to 0$ flattens to $2/h$",
                 fontsize=10.5)
    ax.legend(fontsize=8.5)
    ax.grid(alpha=0.3)

    ax = axes[1]
    u = np.linspace(-12.0, 12.0, 60001)
    d_u = u[1] - u[0]
    xs = np.linspace(-LAMBDA / 2, LAMBDA / 2, 201)
    eps_l = np.array([2.0, 1.0, 0.5, 0.25, 0.125, 0.0625])
    coup = -V_P * (np.cos(K_MODE * (xs[:, None] + u[None, :] / 2))
                   - np.cos(K_MODE * (xs[:, None] - u[None, :] / 2)))
    acts = []
    for eps in eps_l:
        rho_e = 2.0 * np.exp(-u ** 2 / (2.0 * eps ** 2)) \
            / (eps * np.sqrt(2.0 * np.pi))
        acts.append((np.abs(coup) * rho_e[None, :]).sum(axis=1).max() * d_u)
    ax.loglog(eps_l, acts, "o-", color="#1f6feb", lw=1.8, label="measured")
    ax.loglog(eps_l, 2.0 * V_P * K_MODE * np.sqrt(2.0 / np.pi) * eps_l, "--",
              color="crimson", lw=1.4,
              label=r"$2|V'|\sqrt{2/\pi}\,\epsilon$")
    ax.set_xlabel(r"background coherence length $\epsilon$")
    ax.set_ylabel(r"activity $A(\epsilon)$")
    ax.set_title(r"(b) only $\epsilon = 0$ is dark", fontsize=10.5)
    ax.legend(fontsize=8.5)
    ax.grid(alpha=0.3, which="both")

    ax = axes[2]
    xs = np.linspace(-LAMBDA / 2, LAMBDA / 2, 301)
    ys = np.linspace(-2.0 * LAMBDA, 2.0 * LAMBDA, 301)
    xx, yy = np.meshgrid(xs, ys, indexing="xy")
    prec = -(-V_P * (np.cos(K_MODE * (xx + yy / 2))
                     - np.cos(K_MODE * (xx - yy / 2)))) / HBAR
    lim = np.abs(prec).max()
    im = ax.pcolormesh(xs, ys, prec, cmap="RdBu_r", vmin=-lim, vmax=lim,
                       shading="auto")
    ax.axhline(0.0, color="k", lw=1.8)
    ax.text(0.02, 0.03, r"$Y = 0$: dark for every $X$",
            transform=ax.transAxes, fontsize=9)
    ax.set_xlabel(r"midpoint $X$")
    ax.set_ylabel(r"leg separation $Y$")
    ax.set_title(r"(c) $d\mu/dt = -2\,\Gamma_q(X)\sin(kY/2)$", fontsize=10.5)
    fig.colorbar(im, ax=ax, fraction=0.046)
    fig.tight_layout()
    save_fig(fig, "sea_identity_darkness.png")
    plt.close(fig)


def fig_sectors():
    big_x = np.linspace(-7.0, 7.0, 1401)
    psi = cat_state(0.5, 4.0)
    kern, big_y = pair_kernel(psi, big_x, n_y=4096, y_max=24.0)
    big_w, p, d_p = wigner(kern, big_y)
    big_w = big_w.real
    d_y = big_y[1] - big_y[0]
    i0 = int(np.argmin(np.abs(big_y)))
    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.2))

    ax = axes[0]
    ax.plot(big_x, big_w.sum(axis=1) * d_p, lw=2.0, color="#1f6feb",
            label=r"$\int W\,dp$  (net excess per column)")
    ax.plot(big_x, np.abs(psi(big_x)) ** 2, "--", lw=1.4, color="crimson",
            label=r"$\rho(X,X)$")
    ax.axhline(0.0, color="k", lw=0.8)
    ax.set_xlabel(r"position $X$")
    ax.set_ylabel("density")
    ax.set_title("(a) the column excess IS the diagonal, never negative",
                 fontsize=10.5)
    ax.legend(fontsize=8.5)
    ax.grid(alpha=0.3)

    ax = axes[1]
    off = kern.copy()
    off[:, i0] = 0.0
    w_o = np.fft.fftshift(
        np.fft.fft(np.fft.ifftshift(off, axes=1), axis=1), axes=1
    ).real * d_y / (2.0 * np.pi * HBAR)
    sel = np.abs(p) < 4.0
    lim = np.abs(w_o[:, sel]).max()
    im = ax.pcolormesh(big_x, p[sel], w_o[:, sel].T, cmap="RdBu_r",
                       vmin=-lim, vmax=lim, shading="auto")
    ax.set_xlabel(r"position $X$")
    ax.set_ylabel(r"momentum $p$")
    ax.set_title("(b) off-diagonal sector: column sums vanish", fontsize=10.5)
    fig.colorbar(im, ax=ax, fraction=0.046)

    ax = axes[2]
    n_x, n_p, qq = 64, 64, 1
    d_xl = LAMBDA / n_x
    xg = (np.arange(n_x) + 0.5) * d_xl - LAMBDA / 2
    ng = np.arange(n_p) - n_p // 2
    gam = -(V_P / HBAR) * np.sin(2.0 * np.pi * qq * xg / LAMBDA)
    rng = np.random.default_rng(3)
    e = rng.normal(size=(n_p, n_x))
    col0, row0 = e.sum(axis=0).copy(), e.sum(axis=1).copy()
    kk = np.fft.fftfreq(n_p, d=1.0 / n_p)
    sym = 1j * 2.0 * np.sin(2.0 * np.pi * qq * kk / n_p)
    e_j = np.real(np.fft.ifft(
        np.fft.fft(e, axis=0) * np.exp(0.137 * np.outer(sym, gam)), axis=0))
    e_s = np.empty_like(e)
    for i in range(n_p):
        e_s[i] = np.roll(e[i], int(ng[i]) % n_x)
    vals = [np.abs(e_j.sum(axis=0) - col0).max(),
            np.abs(e_j.sum(axis=1) - row0).max(),
            np.abs(e_s.sum(axis=0) - col0).max(),
            np.abs(e_s.sum(axis=1) - row0).max()]
    ax.bar(range(4), np.maximum(vals, 1e-16),
           color=["#1f6feb", "crimson", "crimson", "#1f6feb"])
    ax.set_yscale("log")
    ax.set_xticks(range(4))
    ax.set_xticklabels(["jump\ncolumns", "jump\nrows", "stream\ncolumns",
                        "stream\nrows"], fontsize=8.5)
    ax.set_ylabel("max marginal change")
    ax.set_title("(c) complementary conservation", fontsize=10.5)
    ax.grid(alpha=0.3, axis="y", which="both")
    fig.tight_layout()
    save_fig(fig, "excess_sectors.png")
    plt.close(fig)


def fig_sea_sizing(runs, curves):
    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.2))

    ax = axes[0]
    on, off = curves["rec_on"], curves["rec_off"]
    ax.semilogy(off["ledger"], color="crimson", lw=1.8,
                label="ledger, recombination OFF")
    ax.semilogy(off["sea"], color="crimson", lw=1.3, ls="--",
                label=r"sea $2S$, OFF")
    ax.semilogy(on["ledger"], color="#1f6feb", lw=1.8,
                label="ledger, recombination ON")
    ax.semilogy(on["sea"], color="#1f6feb", lw=1.3, ls="--",
                label=r"sea $2S$, ON")
    ax.set_xlabel("substep")
    ax.set_ylabel("world-particles")
    ax.set_title("(a) splitting draws on a reservoir that must be refilled",
                 fontsize=10.5)
    ax.legend(fontsize=7.5)
    ax.grid(alpha=0.3, which="both")

    ax = axes[1]
    betas = np.array([0.01, 0.03, 0.1, 0.3, 1.0])
    for (mx, npp, vals), col in zip(R["floor"],
                                    ["#0b3d91", "#3d78c2", "#7aa7dd", "#c2d6f0"]):
        ax.semilogx(betas, vals, "o-", color=col, lw=1.6,
                    label=rf"$M_x={mx},\ N_p={npp}$")
    ax.set_xlabel(r"$\beta = B / B_{\rm Wigner}$")
    ax.set_ylabel("blocked-split fraction")
    ax.set_title(r"(b) sizing floor, universal in units of $B_W$",
                 fontsize=10.5)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")

    ax = axes[2]
    for (mx, npp), col in zip([(32, 48), (64, 96)], ["#1f6feb", "crimson"]):
        u = curves[(mx, npp, "uniform")]
        a = curves[(mx, npp, "adaptive")]
        ax.semilogy(u["ledger"] + u["sea"], color=col, lw=1.8, ls="--",
                    label=rf"uniform $B_W$, $N_p={npp}$")
        ax.semilogy(a["ledger"] + a["sea"], color=col, lw=1.8,
                    label=rf"adaptive, $N_p={npp}$")
    ax.set_xlabel("substep")
    ax.set_ylabel(r"$N_{\rm total}$")
    ax.set_title("(c) adaptive allocation, identical error", fontsize=10.5)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    save_fig(fig, "sea_sizing_and_annihilation.png")
    plt.close(fig)


def main():
    print("Verification for docs/analysis/species_sectors_and_annihilation.md")
    part_a()
    part_b()
    part_c()
    part_d()
    part_e()
    part_f()
    rho_abs = part_g()
    part_h()
    runs = part_i()
    curves = part_j(rho_abs)
    banner("Figures")
    fig_duality()
    fig_sea()
    fig_sectors()
    fig_sea_sizing(runs, curves)
    print("\ndone.")


if __name__ == "__main__":
    main()
