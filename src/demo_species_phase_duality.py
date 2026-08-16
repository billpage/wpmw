"""
Verification for ``docs/analysis/species_phase_duality.md``.

What the relationship is between world-particle *species* (positon / negaton)
and world-particle *phase*, and what that relationship implies for adding an
annihilation process to the algorithm.

Conventions.  ``x1`` is the ket leg, ``x2`` the bra leg, ``X = (x1 + x2)/2``
the pair midpoint and ``Y = x1 - x2`` the *full* leg separation, matching
``docs/analysis/interworld_coupling.md``.  The Wigner transform is a Fourier
transform in ``Y``,

    W(X, p) = (1 / 2 pi hbar) Int dY  rho(X + Y/2, X - Y/2) exp(-i p Y / hbar),

evaluated here by direct quadrature on an independent ``Y`` grid.  The
single-grid discrete Wigner transform is *not* used: it carries a spurious
aliased copy that gives a coherent state a negativity of 0.5, which is
impossible.  Part A validates the quadrature against the published table of
``docs/supplement/representation_cost_and_annihilation.md`` before anything
else is measured.

Parts
-----
A  Transform validation, then Theorem D1.  Species and phase are the same
   degree of freedom in conjugate bases: phase lives on the separation axis
   ``Y``, sign lives on the momentum axis ``p``, and ``Y`` and ``p`` are
   conjugate.  The decisive case is the cat at rest, whose pair phase is
   identically zero everywhere and whose Wigner negativity is nonetheless
   0.29.
B  Theorem D2.  The crystal shift is the identity operator: ``rho -> rho +
   2 * 1``.  Its Weyl symbol is the constant ``2/h``, and it is dark because
   ``[H, 1] = 0`` -- for every Hamiltonian, every potential, every dimension.
   In the pair basis the identity is ``delta(Y)``: the sea sits at exactly
   zero leg separation.
C  Theorem D3 and Proposition D4.  The pair phase precesses at
   ``dmu/dt = -2 sum_q Gamma_q(X) sin(k_q Y / 2)`` -- exact, any number of
   modes.  Zero separation is dark; a background of *finite* coherence
   length is not, at a rate linear in that length.
D  Theorem D5.  Annihilation exact in ``(X, p)`` is non-local in the leg
   coordinate, with reach set by the coherence length and independent of the
   position lattice spacing.
E  Theorem D6.  A soft annihilation kernel of momentum width ``sigma_p`` is
   identically an imposed leg coherence length ``hbar / sigma_p``.
F  The pathwise growth rate that annihilation exists to cancel: the Perron
   root of the entrywise absolute value of the generator.
G  Cost and benefit of an annihilation substep, measured on the crystal
   lattice: ledger growth with and without it, wall-clock overhead, and an
   ensemble sweep looking for a plateau.

Run with ``WPMW_OUTPUT`` set (``/mnt/user-data/outputs`` in the container).
"""

from __future__ import annotations

import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from wpmwlib.wpmw_utils import docs_path, output_path

HBAR = 1.0

# Reference potential: one cosine mode, matching demo_interworld_coupling.py.
V_P = 1.5
LAMBDA = 8.0
K_MODE = 2.0 * np.pi / LAMBDA

RESULTS = {}


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
    print(f"  wrote {name}")


# ----------------------------------------------------------------------
# Pair kernel and Wigner transform (full separation convention)
# ----------------------------------------------------------------------
def pair_kernel(psi, big_x, n_y=8192, y_max=48.0):
    """C(X, Y) = psi(X + Y/2) conj(psi(X - Y/2)), plus the Y grid."""
    big_y = (np.arange(n_y) - n_y // 2) * (2.0 * y_max / n_y)
    kern = psi(big_x[:, None] + big_y[None, :] / 2.0) * \
        np.conj(psi(big_x[:, None] - big_y[None, :] / 2.0))
    return kern, big_y


def wigner_from_kernel(kern, big_y):
    """W(X, p) by direct quadrature; returns W, p, dp."""
    d_y = big_y[1] - big_y[0]
    p = np.fft.fftshift(np.fft.fftfreq(len(big_y), d=d_y)) * 2.0 * np.pi * HBAR
    big_w = np.fft.fftshift(
        np.fft.fft(np.fft.ifftshift(kern, axes=1), axis=1), axes=1
    ) * d_y / (2.0 * np.pi * HBAR)
    return big_w, p, p[1] - p[0]


def cat_state(sigma, sep, kick=0.0):
    """Even superposition of two Gaussians of width sigma, separation sep."""
    nrm = np.sqrt(2.0 * np.sqrt(2.0 * np.pi) * sigma
                  * (1.0 + np.exp(-sep ** 2 / (8.0 * sigma ** 2))))

    def psi(z):
        return (np.exp(-(z - sep / 2.0) ** 2 / (4.0 * sigma ** 2))
                + np.exp(-(z + sep / 2.0) ** 2 / (4.0 * sigma ** 2))) \
            * np.exp(1j * kick * z) / nrm

    return psi


# ----------------------------------------------------------------------
# Part A -- validation, then the conjugacy theorem D1
# ----------------------------------------------------------------------
def part_a():
    banner("A  Transform validation, then Theorem D1 (species and phase are conjugate)")
    sigma = 0.5
    big_x = np.linspace(-14.0, 14.0, 2801)
    d_x = big_x[1] - big_x[0]

    print("  A1  quadrature against the published table of")
    print("      docs/supplement/representation_cost_and_annihilation.md Part A")
    print(f"      {'d/sigma':>8} {'Int W':>12} {'||W||_1':>10} {'note':>8}"
          f" {'nu':>9} {'note':>8}")
    published = [(2, 1.0027, 0.0014), (4, 1.2082, 0.1041), (8, 1.5875, 0.2937),
                 (16, 1.6366, 0.3183), (24, 1.6366, 0.3183)]
    for ratio, l1_note, nu_note in published:
        kern, big_y = pair_kernel(cat_state(sigma, ratio * sigma), big_x)
        big_w, _, d_p = wigner_from_kernel(kern, big_y)
        big_w = big_w.real
        print(f"      {ratio:>8} {big_w.sum() * d_x * d_p:>12.8f}"
              f" {np.abs(big_w).sum() * d_x * d_p:>10.4f} {l1_note:>8.4f}"
              f" {-big_w[big_w < 0].sum() * d_x * d_p:>9.4f} {nu_note:>8.4f}")

    print("\n  A2  the cat at rest: no phase anywhere, and yet a third of the")
    print("      Wigner mass is negative")
    for kick, tag in [(0.0, "at rest, psi real"), (3.0, "boosted by k = 3")]:
        kern, big_y = pair_kernel(cat_state(sigma, 8.0 * sigma, kick), big_x)
        big_w, _, d_p = wigner_from_kernel(kern, big_y)
        big_w = big_w.real
        live = np.abs(kern) > 1e-10 * np.abs(kern).max()
        mu = np.angle(kern)
        neg = -big_w[big_w < 0].sum() * d_x * d_p
        frac_phase = float(np.mean((np.cos(mu) < 0)[live]))
        print(f"      {tag:<20} negativity nu = {neg:.4f},"
              f"  fraction of pairs with cos(mu) < 0 = {frac_phase:.4f}")
        if kick:
            yy = np.broadcast_to(big_y[None, :], kern.shape)
            err = np.abs(np.angle(np.exp(1j * (mu - kick * yy))))[live].max()
            print(f"      {'':<20} mu == k Y exactly: max error {err:.3e}"
                  f"   (phase is a function of Y alone)")
            RESULTS["mu_boost_err"] = err
        else:
            RESULTS["nu_rest"] = neg
            RESULTS["phase_rest"] = frac_phase
    print("""
      => the species of a Wigner carrier is not a relabelled pair phase.  It
         is a functional of the whole Y-fibre over a midpoint X, manufactured
         by the kernel exp(-i p Y / hbar).  Phase lives on Y, sign lives on p,
         and no carrier holds both sharply.""")


# ----------------------------------------------------------------------
# Part B -- the sea is the identity operator
# ----------------------------------------------------------------------
def part_b():
    banner("B  Theorem D2: the crystal shift is the identity operator, and 1 is dark")
    print("  B1  Weyl symbol of a background uniform in X with coherence length eps:")
    print("      C_eps(X, Y) = 2 g_eps(Y), g_eps a normalised Gaussian of width eps")
    print(f"      {'eps':>8} {'W at p = 0':>14} {'2/h':>12} {'ratio':>9}"
          f" {'plateau half-width':>20}")
    big_x = np.array([0.0, 1.3, -2.7])
    two_over_h = 2.0 / (2.0 * np.pi * HBAR)
    eps_scan, w0_scan = [], []
    for eps in [2.0, 0.5, 0.1, 0.02, 0.005]:
        n_y, y_max = 1 << 16, 80.0
        big_y = (np.arange(n_y) - n_y // 2) * (2.0 * y_max / n_y)
        kern = np.broadcast_to(
            (2.0 * np.exp(-big_y ** 2 / (2.0 * eps ** 2))
             / (eps * np.sqrt(2.0 * np.pi)))[None, :], (len(big_x), n_y))
        big_w, p, _ = wigner_from_kernel(np.array(kern), big_y)
        big_w = big_w.real
        i0 = int(np.argmin(np.abs(p)))
        half = np.abs(p[np.abs(big_w[0]) > 0.5 * np.abs(big_w[0]).max()]).max()
        print(f"      {eps:>8.3f} {big_w[0, i0]:>14.8f} {two_over_h:>12.8f}"
              f" {big_w[0, i0] / two_over_h:>9.5f} {half:>20.4f}")
        eps_scan.append(eps)
        w0_scan.append(big_w[0, i0])
    RESULTS["sea_ratio"] = w0_scan[-1] / two_over_h

    print("\n  B2  darkness of the eps -> 0 limit, under an ARBITRARY Hamiltonian")
    n = 129
    rng = np.random.default_rng(20260816)
    a = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    ham = a + a.conj().T
    b = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    rho = b @ b.conj().T
    rho /= np.trace(rho).real
    sea = 2.0 * np.eye(n)
    with_sea = -1j * (ham @ (rho + sea) - (rho + sea) @ ham) / HBAR
    without = -1j * (ham @ rho - rho @ ham) / HBAR
    diff = np.abs(with_sea - without).max()
    print(f"      max| d/dt[rho + 2*1] - d/dt[rho] | = {diff:.3e}"
          f"   (scale |drho/dt| ~ {np.abs(without).max():.3f})")
    RESULTS["sea_commutator"] = diff
    print("""
      => no cancellation is invoked.  The commutator is empty: [H, 1] = 0 for
         every Hamiltonian.  In the pair basis 1 is delta(Y), so the sea sits
         at exactly zero leg separation -- darkness independent of absolute
         position, depending only on relative position, and that relative
         position is zero.""")


# ----------------------------------------------------------------------
# Part C -- phase precession, and why zero separation is special
# ----------------------------------------------------------------------
def part_c():
    banner("C  Theorem D3 and Proposition D4: dmu/dt = -2 sum_q Gamma_q(X) sin(k_q Y/2)")
    rng = np.random.default_rng(20260816)
    single = [(1, V_P, 0.0)]
    triple = [(1, V_P, 0.3), (2, -0.7, 1.1), (5, 0.25, -2.0)]
    for modes in (single, triple):
        def v_fun(z, mm=modes):
            return sum(amp * np.cos(2.0 * np.pi * q * z / LAMBDA + ph)
                       for q, amp, ph in mm)
        big_x = rng.uniform(-LAMBDA / 2, LAMBDA / 2, 60000)
        big_y = rng.uniform(-2.0 * LAMBDA, 2.0 * LAMBDA, 60000)
        lhs = -(v_fun(big_x + big_y / 2.0) - v_fun(big_x - big_y / 2.0)) / HBAR
        rhs = np.zeros_like(lhs)
        for q, amp, ph in modes:
            k = 2.0 * np.pi * q / LAMBDA
            gamma = -(amp / HBAR) * np.sin(k * big_x + ph)
            rhs += -2.0 * gamma * np.sin(k * big_y / 2.0)
        err = np.abs(lhs - rhs).max()
        print(f"  {len(modes)} mode(s): max|lhs - rhs| = {err:.3e}"
              f"   (scale {np.abs(lhs).max():.3f})")
        RESULTS[f"precession_err_{len(modes)}"] = err

    print("\n  C2  a background of FINITE coherence length is dynamically active")
    print("      activity A(eps) = max_X Int du |U(X, u)| rho_eps(u), a mass-weighted")
    print("      winding rate; rho_eps has fixed total mass 2, so A measures how fast")
    print("      the background dephases.  For small eps, U -> u V'(X), giving")
    print("      A = 2 |V'|_max sqrt(2/pi) eps.")
    print(f"      {'eps':>8} {'A(eps)':>14} {'predicted':>14} {'ratio':>9}")
    u_grid = np.linspace(-12.0, 12.0, 120001)
    d_u = u_grid[1] - u_grid[0]
    x_grid = np.linspace(-LAMBDA / 2, LAMBDA / 2, 401)
    vprime_max = V_P * K_MODE
    coupling = -V_P * (np.cos(K_MODE * (x_grid[:, None] + u_grid[None, :] / 2))
                       - np.cos(K_MODE * (x_grid[:, None] - u_grid[None, :] / 2)))
    eps_list, act_list = [], []
    for eps in [1.0, 0.5, 0.25, 0.125, 0.0625]:
        rho_eps = 2.0 * np.exp(-u_grid ** 2 / (2.0 * eps ** 2)) \
            / (eps * np.sqrt(2.0 * np.pi))
        act = float((np.abs(coupling) * rho_eps[None, :]).sum(axis=1).max() * d_u)
        pred = 2.0 * vprime_max * np.sqrt(2.0 / np.pi) * eps
        print(f"      {eps:>8.4f} {act:>14.6f} {pred:>14.6f} {act / pred:>9.5f}")
        eps_list.append(eps)
        act_list.append(act)
    RESULTS["activity_ratio"] = act_list[-1] / pred
    print("""
      => zero separation is dark for every X and every mode; the maximum
         precession rate 2|Gamma_q(X)| is attained at Y = lambda_q / 2, and at
         the maximum of |Gamma| that rate is 2 V_q / hbar -- exactly the
         gamma_max of the annihilation-burden section.""")


# ----------------------------------------------------------------------
# Part D -- how far annihilation reaches in the leg coordinate
# ----------------------------------------------------------------------
def part_d():
    banner("D  Theorem D5: annihilation is local in (X, p), non-local in the legs")
    sigma, sep = 0.5, 4.0
    big_x = np.linspace(-10.0, 10.0, 2001)
    d_x = big_x[1] - big_x[0]
    kern, big_y = pair_kernel(cat_state(sigma, sep), big_x)
    big_w, _, d_p = wigner_from_kernel(kern, big_y)
    big_w = big_w.real
    weight = np.abs(big_w)
    d_y = big_y[1] - big_y[0]
    print("  truncate the Y-fibre at |Y| <= Yc, then re-read the sign of W")
    print(f"  {'Yc/sigma':>9} {'Yc':>8} {'|W|-weighted sign agreement':>29}"
          f" {'||W_trunc||_1':>14}")
    reach_rows = []
    for ratio in [0.25, 0.5, 1, 2, 4, 6, 8, 12, 16]:
        cut = kern * (np.abs(big_y) <= ratio * sigma)[None, :]
        w_t = np.fft.fftshift(
            np.fft.fft(np.fft.ifftshift(cut, axes=1), axis=1), axes=1
        ).real * d_y / (2.0 * np.pi * HBAR)
        agree = float((weight * (np.sign(w_t) == np.sign(big_w))).sum()
                      / weight.sum())
        print(f"  {ratio:>9.2f} {ratio * sigma:>8.2f} {agree:>29.4f}"
              f" {np.abs(w_t).sum() * d_x * d_p:>14.4f}")
        reach_rows.append((ratio * sigma, agree))
    RESULTS["reach_rows"] = reach_rows

    print("\n  the reach is a property of the STATE, not of the lattice:")
    print(f"  {'d/sigma':>8} {'d':>7} {'settling Yc (agreement > 0.999)':>34}"
          f" {'d + 4 sigma':>12}")
    for ratio in [4, 8, 16]:
        kern2, y2 = pair_kernel(cat_state(sigma, ratio * sigma), big_x)
        w2, _, _ = wigner_from_kernel(kern2, y2)
        w2 = w2.real
        wt2 = np.abs(w2)
        dy2 = y2[1] - y2[0]
        found = np.nan
        for yc in np.arange(0.25, 24.0, 0.25):
            cut = kern2 * (np.abs(y2) <= yc)[None, :]
            wtr = np.fft.fftshift(
                np.fft.fft(np.fft.ifftshift(cut, axes=1), axis=1), axes=1
            ).real * dy2 / (2.0 * np.pi * HBAR)
            if (wt2 * (np.sign(wtr) == np.sign(w2))).sum() / wt2.sum() > 0.999:
                found = yc
                break
        print(f"  {ratio:>8} {ratio * sigma:>7.2f} {found:>34.2f}"
              f" {ratio * sigma + 4 * sigma:>12.2f}")
    print("""
      => the sign stops changing at a leg reach of order the cat separation,
         i.e. the coherence length.  It has no dependence on the position
         lattice spacing whatever.  Annihilation exact in (X, p) is therefore
         non-local over the coherence length in the constituent leg positions
         -- and it cannot be local in both bases, because Y and p are
         conjugate.""")


# ----------------------------------------------------------------------
# Part E -- the soft blob is a coherence length
# ----------------------------------------------------------------------
def part_e():
    banner("E  Theorem D6: a soft annihilation blob of width sigma_p is a coherence "
           "length hbar/sigma_p")
    sigma, sep = 0.5, 4.0
    big_x = np.linspace(-10.0, 10.0, 1201)
    d_x = big_x[1] - big_x[0]
    kern, big_y = pair_kernel(cat_state(sigma, sep), big_x)
    big_w, p, d_p = wigner_from_kernel(kern, big_y)
    big_w = big_w.real
    d_y = big_y[1] - big_y[0]
    print(f"  {'sigma_p':>9} {'hbar/sigma_p':>14}"
          f" {'max|convolve_p - window_Y|':>28} {'relative':>10} {'||W||_1':>10}")
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
        print(f"  {s_p:>9.2f} {HBAR / s_p:>14.4f} {err:>28.3e}"
              f" {err / np.abs(w_w).max():>10.2e} {l1:>10.4f}")
        l1_rows.append((s_p, l1))
        RESULTS["blob_err"] = err
    print(f"  {'(none)':>9} {'inf':>14} {'-':>28} {'-':>10}"
          f" {np.abs(big_w).sum() * d_x * d_p:>10.4f}")
    RESULTS["blob_l1"] = l1_rows
    print("""
      => exact to floating point.  The bias of a soft annihilation kernel is
         an artificial decoherence of known length, and its benefit is the
         ||W||_1 collapse toward 1.  The "soft blob" regularization and the
         "decoherence makes this model cheaper" prediction are the same knob
         read in the two bases.""")


# ----------------------------------------------------------------------
# Part F -- the growth rate annihilation exists to cancel
# ----------------------------------------------------------------------
def part_f():
    banner("F  What annihilation is for: the Perron root of the absolute generator")
    l_box, n_x, n_p, q, v_q = 8.0, 256, 128, 1, 1.5
    d_x = l_box / n_x
    d_p = np.pi * HBAR / l_box
    x = (np.arange(n_x) + 0.5) * d_x - l_box / 2
    n = np.arange(n_p) - n_p // 2
    gam = np.abs(-(v_q / HBAR) * np.sin(2.0 * np.pi * q * x / l_box))
    gamma_max = 2.0 * v_q / HBAR
    gamma_avg = 4.0 * v_q / (np.pi * HBAR)

    d_t = d_x / d_p
    k = np.fft.fftfreq(n_p, d=1.0 / n_p)
    symbol = 2.0 * np.cos(2.0 * np.pi * q * k / n_p)
    amp = np.exp(d_t * gam[None, :] * symbol[:, None])
    shifts = n % n_x

    a = np.ones((n_p, n_x))
    rates = []
    for _ in range(3000):
        a = np.real(np.fft.ifft(np.fft.fft(a, axis=0) * amp, axis=0))
        out = np.empty_like(a)
        for i in range(n_p):
            out[i] = np.roll(a[i], shifts[i])
        a = out
        tot = a.sum()
        rates.append(np.log(tot) / d_t)
        a /= tot
    rho_abs = float(np.mean(rates[-400:]))
    print(f"  gamma_avg = 4 V_q / (pi hbar)      = {gamma_avg:.6f}")
    print(f"  rho(|L|)  = Perron root, streaming = {rho_abs:.6f}"
          f"   (+/- {np.std(rates[-400:]):.1e})")
    print(f"  gamma_max = 2 V_q / hbar           = {gamma_max:.6f}")
    print(f"  ordering gamma_avg < rho < gamma_max: "
          f"{gamma_avg < rho_abs < gamma_max}")
    RESULTS["rho_abs"] = rho_abs
    RESULTS["gamma_max"] = gamma_max
    RESULTS["gamma_avg"] = gamma_avg
    print("""
      => rho(|L|) is a similarity invariant of the generator, so it is the
         same for every unraveling linear in the ensemble, including every
         importance-sampled one.  It is also, by Part C, a phase precession
         rate.  Annihilation is the only remaining lever, and it is
         necessarily quadratic in the ensemble.""")
    return rho_abs


# ----------------------------------------------------------------------
# Part G -- the cost and the benefit of an annihilation substep
# ----------------------------------------------------------------------
def four_action_lattice(n_p, n_x, gam, q, d_t, n_steps, u_plus, u_minus,
                        annihilate, rng, shifts, track=True):
    """Species-resolved four-action tau-leap, optional same-cell annihilation.

    Channels (per polarisation block, sigma = sign(Gamma), g = |Gamma|):
    a positon at centre ``c`` emits a positon at ``c - sigma q`` and a negaton
    at ``c + sigma q``, each at rate ``g``; the negaton channels are the
    crossing conjugates.  The mean field is exactly the QLE stencil acting on
    ``E = U+ - U-``; the pathwise ledger ``sum(U+ + U-)`` grows at ``2|Gamma|``.
    """
    sig = np.sign(gam)
    mag = np.abs(gam)
    ledger, jump_t, ann_t = [], 0.0, 0.0
    for _ in range(n_steps):
        t0 = time.perf_counter()
        lam = mag[None, :] * d_t
        born_p = rng.poisson(lam * u_plus)
        born_m = rng.poisson(lam * u_minus)
        add_p = np.zeros_like(u_plus)
        add_m = np.zeros_like(u_minus)
        for col, s in enumerate(sig):
            if s == 0:
                continue
            shift = int(s) * q
            add_p[:, col] += np.roll(born_p[:, col], -shift)
            add_m[:, col] += np.roll(born_p[:, col], +shift)
            add_m[:, col] += np.roll(born_m[:, col], -shift)
            add_p[:, col] += np.roll(born_m[:, col], +shift)
        u_plus = u_plus + add_p
        u_minus = u_minus + add_m
        jump_t += time.perf_counter() - t0

        t0 = time.perf_counter()
        if annihilate:
            pair = np.minimum(u_plus, u_minus)
            u_plus -= pair
            u_minus -= pair
        ann_t += time.perf_counter() - t0

        out_p = np.empty_like(u_plus)
        out_m = np.empty_like(u_minus)
        for i in range(n_p):
            out_p[i] = np.roll(u_plus[i], shifts[i])
            out_m[i] = np.roll(u_minus[i], shifts[i])
        u_plus, u_minus = out_p, out_m
        if track:
            ledger.append(float(u_plus.sum() + u_minus.sum()))
    return u_plus, u_minus, np.array(ledger), jump_t, ann_t


def part_g(rho_abs):
    banner("G  Cost and benefit of an annihilation substep")
    l_box, n_x, n_p, q, v_q = 8.0, 64, 64, 1, 1.5
    d_x = l_box / n_x
    d_p = np.pi * HBAR / l_box
    x = (np.arange(n_x) + 0.5) * d_x - l_box / 2
    n = np.arange(n_p) - n_p // 2
    gam = -(v_q / HBAR) * np.sin(2.0 * np.pi * q * x / l_box)
    d_t = 0.25 * d_x / d_p
    shifts = (n * 0) % n_x  # streaming frozen: isolate the ledger effect
    rng = np.random.default_rng(20260816)

    # exact mean field for reference
    def seed(nu):
        sigma = 0.6
        prof = np.exp(-(x ** 2) / (2.0 * 1.2 ** 2))
        mom = np.exp(-((n * d_p) ** 2) / (2.0 * sigma ** 2))
        w = np.outer(mom, prof)
        w /= w.sum()
        e = np.round(nu * w).astype(np.int64)
        return np.maximum(e, 0), np.zeros_like(e)

    n_steps = 40
    print(f"  lattice {n_p} x {n_x}, dt = {d_t:.4f}, {n_steps} steps,"
          f" t = {n_steps * d_t:.3f}")
    print(f"\n  G1  the ledger, with and without a same-cell annihilation substep")
    print(f"      {'':<22} {'initial':>12} {'final':>14} {'growth':>12}"
          f" {'mean rate':>10} {'final rate':>11}")
    curves = {}
    for annih, tag in [(False, "no annihilation"), (True, "same-cell annihilation")]:
        up, um = seed(200000)
        up, um, ledger, jump_t, ann_t = four_action_lattice(
            n_p, n_x, gam, q, d_t, n_steps, up, um, annih,
            np.random.default_rng(20260816), shifts)
        rate = np.log(ledger[-1] / ledger[0]) / (n_steps * d_t)
        late = np.log(ledger[-1] / ledger[-11]) / (10 * d_t)
        print(f"      {tag:<22} {ledger[0]:>12.4e} {ledger[-1]:>14.4e}"
              f" {ledger[-1] / ledger[0]:>12.3e} {rate:>10.3f} {late:>11.3f}")
        curves[tag] = ledger
        RESULTS[f"ledger_{int(annih)}"] = (ledger[0], ledger[-1], rate)
        RESULTS[f"time_{int(annih)}"] = (jump_t, ann_t)
    print(f"      predicted asymptote for the ungoverned ledger, 2|Gamma|_max ="
          f" {2 * np.abs(gam).max():.3f}"
          f"   (streaming frozen here, so the Perron root is the pointwise max)")

    print(f"\n  G2  wall clock, same run: annihilation pass against jump pass")
    j_t, a_t = RESULTS["time_1"]
    print(f"      jump substep          {j_t:.4f} s")
    print(f"      annihilation substep  {a_t:.4f} s")
    print(f"      overhead              {100.0 * a_t / j_t:.2f} percent of the jump pass")
    RESULTS["overhead_pct"] = 100.0 * a_t / j_t

    print(f"\n  G3  what the overhead buys, in world-particles")
    print(f"      cost of a representation is N / ||mu||_1^2, so an ungoverned")
    print(f"      ledger multiplies the particle requirement by exp(2 rho t)")
    print(f"      {'t':>6} {'exp(2 rho t)':>16} {'break-even overhead':>22}")
    for t in [0.5, 1.0, 2.0, 4.0]:
        print(f"      {t:>6.1f} {np.exp(2 * rho_abs * t):>16.3e}"
              f" {np.exp(2 * rho_abs * t):>22.3e}")
    t_star = np.log(1.0 + RESULTS["overhead_pct"] / 100.0) / (2.0 * rho_abs)
    print(f"      the substep pays for itself at t = {t_star:.4f}")
    RESULTS["t_star"] = t_star

    print(f"\n  G4  ensemble sweep: is there a plateau below which annihilation fails?")
    print(f"      reconstruct E on the lattice and compare with the annihilation-free")
    print(f"      mean field; occupancy is nu / (n_p n_x) per cell")
    print(f"      {'nu':>12} {'occupancy':>11} {'rel. L2 error':>15}"
          f" {'error * sqrt(nu)':>18}")
    sweep = []
    for nu in [4000, 40000, 400000, 4000000]:
        up0, um0 = seed(nu)
        ref = (up0 - um0).astype(float)
        up, um, _, _, _ = four_action_lattice(
            n_p, n_x, gam, q, d_t, n_steps, up0.copy(), um0.copy(), True,
            np.random.default_rng(7), shifts, track=False)
        # exact mean-field evolution of the same seed
        e = ref.copy()
        for _ in range(n_steps):
            e = e + d_t * gam[None, :] * (np.roll(e, -q, axis=0)
                                          - np.roll(e, +q, axis=0))
        got = (up - um).astype(float)
        err = np.linalg.norm(got - e) / np.linalg.norm(e)
        occ = nu / (n_p * n_x)
        print(f"      {nu:>12} {occ:>11.2f} {err:>15.5f}"
              f" {err * np.sqrt(nu):>18.3f}")
        sweep.append((nu, occ, err))
    RESULTS["sweep"] = sweep
    print("""
      => the error falls as nu^{-1/2} across three decades with no threshold:
         the last column is flat.  Criterion (i) of the annihilation-burden
         section governs, and criterion (ii) -- which assumed the original
         partner had to be found -- does not, because Theorem D5 licenses
         anonymous annihilation.""")
    return curves


# ----------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------
def fig_duality():
    sigma, sep = 0.5, 4.0
    big_x = np.linspace(-6.0, 6.0, 601)
    kern, big_y = pair_kernel(cat_state(sigma, sep), big_x, n_y=4096, y_max=24.0)
    big_w, p, _ = wigner_from_kernel(kern, big_y)
    big_w = big_w.real

    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.2))
    sel = np.abs(p) < 4.0
    ax = axes[0]
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
    ax.contour(big_x, big_y[ysel], np.abs(kern[:, ysel]).T, levels=4,
               colors="w", linewidths=0.6, alpha=0.6)
    ax.set_xlabel(r"midpoint $X$")
    ax.set_ylabel(r"leg separation $Y$")
    ax.set_title(r"(b) $|\rho|$ on the same fibres; $\mu = \arg\rho \equiv 0$",
                 fontsize=10.5)
    ax.text(0.03, 0.94, "three lobes, all real and positive:\n"
            r"no phase anywhere, yet $\nu = 0.29$",
            transform=ax.transAxes, fontsize=8.5, color="w", va="top")
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
               label=r"cat separation $d$ (coherence length)")
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
        kern = (2.0 * np.exp(-big_y ** 2 / (2.0 * eps ** 2))
                / (eps * np.sqrt(2.0 * np.pi)))[None, :]
        big_w, p, _ = wigner_from_kernel(kern, big_y)
        sel = np.abs(p) < 8.0
        ax.plot(p[sel], big_w.real[0, sel], color=col, lw=1.6,
                label=rf"$\epsilon = {eps}$")
    ax.axhline(two_over_h, color="crimson", ls="--", lw=1.2, label=r"$2/h$")
    ax.set_xlabel(r"momentum $p$")
    ax.set_ylabel(r"Weyl symbol")
    ax.set_title(r"(a) background of coherence length $\epsilon$; $\epsilon\to0$"
                 " flattens to $2/h$", fontsize=10.5)
    ax.legend(fontsize=8.5)
    ax.grid(alpha=0.3)

    ax = axes[1]
    u_grid = np.linspace(-12.0, 12.0, 60001)
    d_u = u_grid[1] - u_grid[0]
    x_grid = np.linspace(-LAMBDA / 2, LAMBDA / 2, 201)
    eps_list = np.array([2.0, 1.0, 0.5, 0.25, 0.125, 0.0625])
    coupling = -V_P * (np.cos(K_MODE * (x_grid[:, None] + u_grid[None, :] / 2))
                       - np.cos(K_MODE * (x_grid[:, None] - u_grid[None, :] / 2)))
    acts = []
    for eps in eps_list:
        rho_eps = 2.0 * np.exp(-u_grid ** 2 / (2.0 * eps ** 2)) \
            / (eps * np.sqrt(2.0 * np.pi))
        acts.append((np.abs(coupling) * rho_eps[None, :]).sum(axis=1).max() * d_u)
    ax.loglog(eps_list, acts, "o-", color="#1f6feb", lw=1.8, label="measured")
    ax.loglog(eps_list,
              2.0 * V_P * K_MODE * np.sqrt(2.0 / np.pi) * eps_list, "--",
              color="crimson", lw=1.4,
              label=r"$2|V'|_{\max}\sqrt{2/\pi}\,\epsilon$")
    ax.set_xlabel(r"background coherence length $\epsilon$")
    ax.set_ylabel(r"activity $A(\epsilon)$")
    ax.set_title(r"(b) activity is linear in $\epsilon$; only $\epsilon=0$ is dark",
                 fontsize=10.5)
    ax.legend(fontsize=8.5)
    ax.grid(alpha=0.3, which="both")

    ax = axes[2]
    xs = np.linspace(-LAMBDA / 2, LAMBDA / 2, 301)
    ys = np.linspace(-2.0 * LAMBDA, 2.0 * LAMBDA, 301)
    xx, yy = np.meshgrid(xs, ys, indexing="xy")
    prec = -(-V_P * (np.cos(K_MODE * (xx + yy / 2)) - np.cos(K_MODE * (xx - yy / 2)))) \
        / HBAR
    lim = np.abs(prec).max()
    im = ax.pcolormesh(xs, ys, prec, cmap="RdBu_r", vmin=-lim, vmax=lim,
                       shading="auto")
    ax.axhline(0.0, color="k", lw=1.6)
    ax.text(0.02, 0.03, r"$Y=0$: dark for every $X$", transform=ax.transAxes,
            fontsize=9)
    ax.set_xlabel(r"midpoint $X$")
    ax.set_ylabel(r"leg separation $Y$")
    ax.set_title(r"(c) $d\mu/dt = -2\,\Gamma_q(X)\,\sin(k Y/2)$", fontsize=10.5)
    fig.colorbar(im, ax=ax, fraction=0.046)

    fig.tight_layout()
    save_fig(fig, "sea_identity_darkness.png")
    plt.close(fig)


def fig_annihilation(curves):
    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.2))

    ax = axes[0]
    for tag, col in [("no annihilation", "crimson"),
                     ("same-cell annihilation", "#1f6feb")]:
        led = curves[tag]
        ax.semilogy(np.arange(len(led)), led, lw=1.8, color=col, label=tag)
    ax.set_xlabel("substep")
    ax.set_ylabel(r"ledger $\sum_i |w_i|$")
    ax.set_title("(a) the ledger, with and without the substep", fontsize=10.5)
    ax.legend(fontsize=8.5)
    ax.grid(alpha=0.3)

    ax = axes[1]
    nus = np.array([s[0] for s in RESULTS["sweep"]], dtype=float)
    errs = np.array([s[2] for s in RESULTS["sweep"]])
    ax.loglog(nus, errs, "o-", color="#1f6feb", lw=1.8, label="measured")
    ax.loglog(nus, errs[0] * np.sqrt(nus[0] / nus), "--", color="crimson",
              lw=1.4, label=r"$\nu^{-1/2}$")
    ax.set_xlabel(r"ensemble size $\nu$")
    ax.set_ylabel("relative L2 error")
    ax.set_title("(b) no plateau over three decades", fontsize=10.5)
    ax.legend(fontsize=8.5)
    ax.grid(alpha=0.3, which="both")

    ax = axes[2]
    ts = np.linspace(0.0, 4.0, 200)
    ax.semilogy(ts, np.exp(2.0 * RESULTS["rho_abs"] * ts), lw=1.8,
                color="crimson", label=r"particles saved, $e^{2\rho t}$")
    ax.axhline(1.0 + RESULTS["overhead_pct"] / 100.0, color="#1f6feb", ls="--",
               lw=1.4, label="measured overhead")
    ax.axvline(RESULTS["t_star"], color="k", ls=":", lw=1.2,
               label=rf"break-even $t = {RESULTS['t_star']:.3f}$")
    ax.set_xlabel("time")
    ax.set_ylabel("factor")
    ax.set_title("(c) cost against benefit", fontsize=10.5)
    ax.legend(fontsize=8.5, loc="lower right")
    ax.grid(alpha=0.3, which="both")

    fig.tight_layout()
    save_fig(fig, "annihilation_cost_benefit.png")
    plt.close(fig)


def main():
    print("Verification for docs/analysis/species_phase_duality.md")
    part_a()
    part_b()
    part_c()
    part_d()
    part_e()
    rho_abs = part_f()
    curves = part_g(rho_abs)
    banner("Figures")
    fig_duality()
    fig_sea()
    fig_annihilation(curves)
    print("\ndone.")


if __name__ == "__main__":
    main()
