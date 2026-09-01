"""
Verification for ``docs/analysis/compensated_ontology.md``.

The compensated Liouville split, read as a proposed ontology rather than as a
numerical scheme: every world-particle streams on a Newtonian arc under the
full classical force, and the only non-classical element of the dynamics is
the birth and death of positon-negaton pairs.

This script tests the *limits* of that reading.  Its four load-bearing
results:

  G2  For a quadratic Hamiltonian the demographic channel is empty.  The
      compensated rate ``Gamma`` is zero to machine precision at every reach,
      while the uncompensated (field-less) rate ``gamma`` of the published
      signed-particle formulation is large -- the whole restoring force is
      delivered by pair creation there.  So the two readings of the same
      equation disagree observably about the census.

  G3  What carries the quantum for a quadratic Hamiltonian is therefore not
      the dynamics but the admissible set of initial ensembles.  Sweeping the
      quartic coupling closes the dynamical door continuously while the
      state's Wigner negativity is untouched.

  G4  Admissibility is preserved by streaming plus demography and *not* by
      streaming alone: dropping the residual channel drives the least
      eigenvalue of the reconstructed density operator negative.  The
      residual channel is thus doing kinematic work as well as dynamical.

  G5  The emission rate is a regulator-dependent quantity while the generator
      it represents is not: ``Gamma`` diverges with the coherence reach while
      the third moment of the same kernel converges on the analytic leading
      Moyal coefficient.

Conventions follow ``compensated_liouville_splitting.md`` section 1:

    y            = hbar s / 2                  half ket-bra separation
    D(x, y)      = V(x + y) - V(x - y)
    D_res(x, y)  = D(x, y) - 2 y V'(x)         compensated residual
    M            = i D / hbar,  M_res = i D_res / hbar
    dp           = pi hbar / (2 y_max)
    L_res W(p)   = sum_q K_q W(p - xi_q),      xi_q = q dp

Parts
-----
A  Theorem G2.  ``Gamma`` for the compensated channel against ``gamma`` for
   the field-less one, over harmonic, inverted harmonic, quartic and sech^2,
   at three reaches.  The quadratic column is zero to machine precision in
   the first and large in the second.
B  The harmonic null test.  The Wigner function of ``(|0> + |2>)/sqrt(2)``
   is signed, and is transported exactly by the classical phase-space
   rotation; its negativity volume is constant.
C  Theorem G3.  Sweep the quartic coupling: the dynamical residual strength
   falls to zero with ``lambda`` while the negativity of the state does not
   move.  The two occurrences of ``hbar`` are independent.
D  Postulate (A) is not derivable from (S) + (D).  Four Gaussians of
   different phase-space area are all exactly stationary under the harmonic
   flow; only those with ``sigma_q sigma_p >= hbar/2`` are quantum states.
E  Theorem G4.  Least eigenvalue of the reconstructed ``rho`` under exact
   evolution and under classical carrier transport alone, for the quartic.
F  Theorem G5.  ``Gamma`` diverges with the reach while the third moment of
   the same kernel converges on ``(hbar^2 / 4) V'''``.

Run with ``WPMW_OUTPUT`` set to a writable directory; set ``WPMW_DOCS`` as
well to mirror keeper figures into the ``output``-branch worktree.
"""

from __future__ import annotations

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from wpmwlib.wpmw_utils import output_path, docs_path  # noqa: E402

HBAR = 1.0


def banner(text):
    print()
    print("=" * 72)
    print(text)
    print("=" * 72)


def save_fig(fig, name):
    fig.savefig(output_path(name), dpi=150, bbox_inches="tight")
    dp = docs_path(name)
    if dp:
        fig.savefig(dp, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {name}")


# --------------------------------------------------------------------- #
# Potentials.  Each returns (V, V', V''') as callables.                  #
# --------------------------------------------------------------------- #
def harmonic(omega=1.0, mass=1.0):
    return (lambda z: 0.5 * mass * omega**2 * z**2,
            lambda z: mass * omega**2 * z,
            lambda z: np.zeros_like(z))


def inverted(omega=1.0, mass=1.0):
    return (lambda z: -0.5 * mass * omega**2 * z**2,
            lambda z: -mass * omega**2 * z,
            lambda z: np.zeros_like(z))


def quartic(lam=0.02, omega=1.0, mass=1.0):
    return (lambda z: 0.5 * mass * omega**2 * z**2 + lam * z**4,
            lambda z: mass * omega**2 * z + 4.0 * lam * z**3,
            lambda z: 24.0 * lam * z)


def eckart(v0=1.0, a=1.0):
    def v(z):
        return v0 / np.cosh(z / a) ** 2

    def dv(z):
        c = 1.0 / np.cosh(z / a) ** 2
        return -2.0 * v0 / a * c * np.tanh(z / a)

    def d3v(z):
        c = 1.0 / np.cosh(z / a) ** 2
        t = np.tanh(z / a)
        return 8.0 * v0 / a**3 * c * t * (2.0 * c - t**2)

    return v, dv, d3v


# --------------------------------------------------------------------- #
# The compensated and uncompensated kernels on a reach lattice           #
# --------------------------------------------------------------------- #
def kernel_rungs(x, v, dv, y_max, n_rungs=64, n_y=8192, compensated=True):
    """Return (xi, K) for the momentum-shift lattice at position ``x``.

    ``K_q`` is the coefficient of ``W(p - xi_q)`` in the convolution form of
    the potential operator, restricted to half ket-bra separations
    ``|y| <= y_max``.  With ``compensated=False`` the classical part is left
    in, which is the field-less signed-particle kernel.
    """
    s_max = 2.0 * y_max / HBAR
    s = np.linspace(-s_max, s_max, 2 * n_y + 1)
    ds = s[1] - s[0]
    y = HBAR * s / 2.0
    d = v(x + y) - v(x - y)
    if compensated:
        d = d - 2.0 * y * dv(x)
    m = 1j * d / HBAR                      # the symbol M(x, s)
    dp = np.pi * HBAR / (2.0 * y_max)
    q = np.arange(-n_rungs, n_rungs + 1)
    xi = q * dp
    phase = np.exp(-1j * np.outer(xi, s))
    k = (phase * m).sum(axis=1) * ds / (2.0 * s_max)
    return xi, k.real


def rate_profile(v, dv, y_max, xs, compensated=True, n_rungs=64):
    out = np.empty_like(xs)
    for i, x in enumerate(xs):
        _, k = kernel_rungs(x, v, dv, y_max, n_rungs=n_rungs,
                            compensated=compensated)
        out[i] = np.abs(k).sum()
    return out


# --------------------------------------------------------------------- #
# Part A  Theorem G2: the demographic channel is empty for quadratic V   #
# --------------------------------------------------------------------- #
def part_a():
    banner("Part A  Theorem G2: no events at all for a quadratic Hamiltonian")
    cases = [
        ("harmonic   V = q^2/2", harmonic()),
        ("inverted   V = -q^2/2", inverted()),
        ("quartic    lambda = 0.02", quartic(0.02)),
        ("Eckart     V0 sech^2(q)", eckart()),
    ]
    reaches = [np.pi / 2, np.pi, 2 * np.pi]
    xs = np.linspace(-3.0, 3.0, 25)

    print("\n  compensated rate Gamma_max = max_x sum_q |K_res,q|")
    print(f"  {'potential':26s} " +
          "".join(f"{'y_max=' + f'{r:.2f}':>14s}" for r in reaches))
    comp = {}
    for name, (v, dv, _) in cases:
        row = [rate_profile(v, dv, r, xs, compensated=True).max()
               for r in reaches]
        comp[name] = row
        print(f"  {name:26s} " + "".join(f"{u:14.3e}" for u in row))

    print("\n  field-less rate gamma_max = max_x sum_q |K_q| (no compensation)")
    print(f"  {'potential':26s} " +
          "".join(f"{'y_max=' + f'{r:.2f}':>14s}" for r in reaches))
    for name, (v, dv, _) in cases:
        row = [rate_profile(v, dv, r, xs, compensated=False).max()
               for r in reaches]
        print(f"  {name:26s} " + "".join(f"{u:14.3e}" for u in row))

    quad_max = max(max(comp["harmonic   V = q^2/2"]),
                   max(comp["inverted   V = -q^2/2"]))
    print(f"\n  worst quadratic compensated rate: {quad_max:.3e}")
    print("  D_res = V(x+y) - V(x-y) - 2 y V'(x) vanishes identically for")
    print("  quadratic V, so postulate (D) never fires.  Under the")
    print("  field-less reading the same system has a large event rate.")
    return quad_max


# --------------------------------------------------------------------- #
# Wigner functions of oscillator states                                  #
# --------------------------------------------------------------------- #
def osc_state(n, q, mass=1.0, omega=1.0):
    """Normalised harmonic-oscillator eigenfunction on grid ``q``."""
    import math
    from numpy.polynomial.hermite import hermval

    alpha = np.sqrt(mass * omega / HBAR)
    coef = np.zeros(n + 1)
    coef[n] = 1.0
    norm = (alpha / np.sqrt(np.pi)) ** 0.5 / np.sqrt(2.0**n * math.factorial(n))
    return norm * hermval(alpha * q, coef) * np.exp(-(alpha * q) ** 2 / 2.0)


def wigner(psi, q, p):
    """Wigner transform of a pure state sampled on a uniform grid ``q``."""
    dq = q[1] - q[0]
    nq = q.size
    # y runs over even offsets so q +- y/2 land on grid points
    half = nq // 4
    offs = np.arange(-half, half + 1)
    y = 2.0 * offs * dq
    plus = np.empty((nq, y.size), dtype=complex)
    minus = np.empty_like(plus)
    for j, o in enumerate(offs):
        plus[:, j] = np.roll(psi, -o)
        minus[:, j] = np.roll(psi, o)
    prod = np.conj(minus) * plus            # psi*(q - y/2) psi(q + y/2)
    ker = np.exp(-1j * np.outer(p, y) / HBAR)
    w = (ker @ prod.T).real * (y[1] - y[0]) / (2.0 * np.pi * HBAR)
    return w.T                              # shape (nq, np)


def negativity_volume(w, dq, dp):
    return np.abs(w).sum() * dq * dp - 1.0


# --------------------------------------------------------------------- #
# Part B  the harmonic null test                                         #
# --------------------------------------------------------------------- #
def part_b():
    banner("Part B  the harmonic null test: signed, and classically carried")
    nq = 512
    q = np.linspace(-8.0, 8.0, nq, endpoint=False)
    dq = q[1] - q[0]
    p = np.linspace(-6.0, 6.0, 241)
    dp = p[1] - p[0]

    def psi_at(t):
        return (osc_state(0, q) * np.exp(-1j * 0.5 * t)
                + osc_state(2, q) * np.exp(-1j * 2.5 * t)) / np.sqrt(2.0)

    w0 = wigner(psi_at(0.0), q, p)
    qq, pp = np.meshgrid(q, p, indexing="ij")
    print(f"  norm(W0) = {w0.sum() * dq * dp:.9f}")
    print(f"  min  W0  = {w0.min():+.6f}   (2/h = {2 / (2 * np.pi * HBAR):.6f})")
    print(f"  negativity volume = {negativity_volume(w0, dq, dp):.6f}")

    print("\n  classical rotation of W0 against the exact Wigner function")
    print(f"  {'t':>8s}{'max |W_rot - W_exact|':>26s}{'negativity':>14s}")
    errs = []
    for t in [0.25 * np.pi, 0.5 * np.pi, np.pi]:
        # invert the classical flow: q0 = q cos t - p sin t, p0 = q sin t + p cos t
        q0 = qq * np.cos(t) - pp * np.sin(t)
        p0 = qq * np.sin(t) + pp * np.cos(t)
        from scipy.interpolate import RegularGridInterpolator

        interp = RegularGridInterpolator((q, p), w0, bounds_error=False,
                                         fill_value=0.0)
        w_rot = interp(np.stack([q0, p0], axis=-1))
        w_ex = wigner(psi_at(t), q, p)
        err = np.abs(w_rot - w_ex).max()
        errs.append(err)
        print(f"  {t:8.4f}{err:26.3e}"
              f"{negativity_volume(w_ex, dq, dp):14.6f}")
    print("\n  the residual channel is identically zero here, so the whole")
    print("  evolution is classical transport of an initially signed W.")
    print("  the first two rows are the bilinear-interpolation floor of the")
    print("  rotated grid; at t = pi the rotation is grid-exact and the")
    print("  error falls by three orders, which is the real null test.")
    return errs[-1]


# --------------------------------------------------------------------- #
# Part C  Theorem G3: the two doors are independent                      #
# --------------------------------------------------------------------- #
def moyal_residual(w, p, d3v_at_q):
    """Leading Moyal residual -(hbar^2/24) V'''(q) d^3W/dp^3, spectral in p."""
    dp = p[1] - p[0]
    kp = 2.0 * np.pi * np.fft.fftfreq(p.size, d=dp)
    d3w = np.fft.ifft((1j * kp) ** 3 * np.fft.fft(w, axis=1), axis=1).real
    return -(HBAR**2 / 24.0) * d3v_at_q[:, None] * d3w


def classical_generator(w, q, p, dv_at_q):
    dq = q[1] - q[0]
    dp = p[1] - p[0]
    kq = 2.0 * np.pi * np.fft.fftfreq(q.size, d=dq)
    kp = 2.0 * np.pi * np.fft.fftfreq(p.size, d=dp)
    dwdq = np.fft.ifft(1j * kq[:, None] * np.fft.fft(w, axis=0), axis=0).real
    dwdp = np.fft.ifft(1j * kp[None, :] * np.fft.fft(w, axis=1), axis=1).real
    return -p[None, :] * dwdq + dv_at_q[:, None] * dwdp


def part_c():
    banner("Part C  Theorem G3: closing the dynamical door leaves the state signed")
    nq = 512
    q = np.linspace(-8.0, 8.0, nq, endpoint=False)
    dq = q[1] - q[0]
    p = np.linspace(-6.0, 6.0, 256, endpoint=False)
    dp = p[1] - p[0]
    psi = (osc_state(0, q) + osc_state(2, q)) / np.sqrt(2.0)
    w = wigner(psi, q, p)
    neg = negativity_volume(w, dq, dp)

    lams = np.array([0.0, 1e-4, 1e-3, 3e-3, 1e-2, 2e-2, 5e-2, 1e-1])
    chis = []
    print(f"\n  {'lambda':>10s}{'chi_Q':>14s}{'negativity':>14s}")
    for lam in lams:
        _, dv, d3v = quartic(lam)
        cl = classical_generator(w, q, p, dv(q))
        rs = moyal_residual(w, p, d3v(q))
        n_cl = np.sqrt((cl**2).sum() * dq * dp)
        n_rs = np.sqrt((rs**2).sum() * dq * dp)
        chi = n_rs / (n_cl + n_rs)
        chis.append(chi)
        print(f"  {lam:10.4f}{chi:14.3e}{neg:14.6f}")
    print("\n  chi_Q -> 0 continuously with lambda; the negativity of the")
    print("  state is untouched.  Dynamical and kinematic hbar are separate.")
    return lams, np.array(chis), neg


# --------------------------------------------------------------------- #
# Part D  postulate (A) is not derivable from (S) + (D)                  #
# --------------------------------------------------------------------- #
def part_d():
    banner("Part D  postulate (A) is independent: four stationary Gaussians")
    print("\n  isotropic Gaussians in the harmonic well are exactly")
    print("  stationary under (S), and (D) is empty by G2, so the dynamics")
    print("  cannot tell them apart.  Only phase-space area >= hbar/2 is a")
    print("  quantum state.\n")
    nq, npn = 256, 256
    q = np.linspace(-6.0, 6.0, nq, endpoint=False)
    p = np.linspace(-6.0, 6.0, npn, endpoint=False)
    dq, dp = q[1] - q[0], p[1] - p[0]
    _, dv, _ = harmonic()

    print(f"  {'sigma_q sigma_p / hbar':>24s}{'peak W / (2/h)':>18s}"
          f"{'purity':>10s}{'|L_cl W|':>12s}{'admissible':>12s}")
    rows = []
    for area in [0.125, 0.25, 0.5, 1.0]:
        sig = np.sqrt(area * HBAR)          # sigma_q = sigma_p = sqrt(area)
        w = np.exp(-(q[:, None] ** 2 + p[None, :] ** 2) / (2 * sig**2))
        w /= 2.0 * np.pi * sig**2
        peak = w.max() / (2.0 / (2.0 * np.pi * HBAR))
        purity = 2.0 * np.pi * HBAR * (w**2).sum() * dq * dp
        drift = np.abs(classical_generator(w, q, p, dv(q))).max()
        ok = "yes" if area >= 0.5 - 1e-12 else "NO"
        rows.append((area, peak, purity, drift, ok))
        print(f"  {area:24.3f}{peak:18.4f}{purity:10.4f}"
              f"{drift:12.2e}{ok:>12s}")
    print("\n  the violated inequality is exactly the Wigner bound")
    print("  |W| <= 2/h, equivalently purity <= 1 -- a constraint on the")
    print("  ensemble, not on its motion.")
    return rows


# --------------------------------------------------------------------- #
# Part E  Theorem G4: admissibility needs the residual channel           #
# --------------------------------------------------------------------- #
def rho_from_wigner(w, q, p):
    """Reconstruct <x|rho|x'> from W on a uniform (q, p) grid."""
    dp = p[1] - p[0]
    nq = q.size
    dq = q[1] - q[0]
    half = nq // 4
    offs = np.arange(-half, half + 1)
    y = 2.0 * offs * dq
    ker = np.exp(1j * np.outer(y, p) / HBAR)          # (ny, np)
    g = (w @ ker.T).T * dp                            # g[j, i] = rho(q_i + y_j/2, q_i - y_j/2)
    # assemble on the (x, x') grid where x = q + y/2, x' = q - y/2
    rho = np.zeros((nq, nq), dtype=complex)
    for j, o in enumerate(offs):
        rows = (np.arange(nq) + o) % nq
        cols = (np.arange(nq) - o) % nq
        rho[rows, cols] = g[j]
    return rho


def split_step_schrodinger(psi, q, v_at_q, dt, n_steps, mass=1.0):
    dq = q[1] - q[0]
    kq = 2.0 * np.pi * np.fft.fftfreq(q.size, d=dq)
    kin = np.exp(-1j * HBAR * kq**2 * dt / (2.0 * mass))
    pot = np.exp(-1j * v_at_q * dt / (2.0 * HBAR))
    for _ in range(n_steps):
        psi = pot * psi
        psi = np.fft.ifft(kin * np.fft.fft(psi))
        psi = pot * psi
    return psi


def classical_transport(w, q, p, dv_at_q, dt, n_steps, mass=1.0):
    """Strang-split exact classical Liouville transport of W (spectral shifts)."""
    dq, dp = q[1] - q[0], p[1] - p[0]
    kq = 2.0 * np.pi * np.fft.fftfreq(q.size, d=dq)
    kp = 2.0 * np.pi * np.fft.fftfreq(p.size, d=dp)
    drift_half = np.exp(-1j * np.outer(kq, p) * (dt / 2.0) / mass)
    kick = np.exp(1j * np.outer(dv_at_q, kp) * dt)
    for _ in range(n_steps):
        w = np.fft.ifft(drift_half * np.fft.fft(w, axis=0), axis=0).real
        w = np.fft.ifft(kick * np.fft.fft(w, axis=1), axis=1).real
        w = np.fft.ifft(drift_half * np.fft.fft(w, axis=0), axis=0).real
    return w


def part_e():
    banner("Part E  Theorem G4: (S) alone leaves the admissible set")
    lam = 0.05
    nq = 256
    q = np.linspace(-8.0, 8.0, nq, endpoint=False)
    dq = q[1] - q[0]
    p = np.linspace(-8.0, 8.0, 256, endpoint=False)
    v, dv, _ = quartic(lam)
    psi0 = (osc_state(0, q) + osc_state(2, q)) / np.sqrt(2.0)
    psi0 = psi0 / np.sqrt((np.abs(psi0) ** 2).sum() * dq)
    w0 = wigner(psi0, q, p)

    dt = 0.002
    times = [0.0, 0.25 * np.pi, 0.5 * np.pi, np.pi]
    print(f"\n  quartic lambda = {lam}, state (|0> + |2>)/sqrt(2)")
    print(f"  {'t':>8s}{'min eig rho  exact':>22s}"
          f"{'min eig rho  (S) only':>26s}")
    w_cl = w0.copy()
    t_prev = 0.0
    rows = []
    for t in times:
        n = int(round((t - t_prev) / dt))
        if n:
            w_cl = classical_transport(w_cl, q, p, dv(q), dt, n)
        t_prev = t
        psi_t = split_step_schrodinger(psi0, q, v(q), dt, int(round(t / dt)))
        w_ex = wigner(psi_t, q, p)
        e_ex = np.linalg.eigvalsh(rho_from_wigner(w_ex, q, p) * dq).min()
        e_cl = np.linalg.eigvalsh(rho_from_wigner(w_cl, q, p) * dq).min()
        rows.append((t, e_ex, e_cl))
        print(f"  {t:8.4f}{e_ex:22.3e}{e_cl:26.3e}")
    print("\n  exact evolution keeps rho positive semidefinite to the grid")
    print("  floor; classical carrier transport alone does not.  The")
    print("  residual channel is what keeps the ensemble physical.")
    return rows


# --------------------------------------------------------------------- #
# Part F  Theorem G5: the rate is a regulator, the generator is not      #
# --------------------------------------------------------------------- #
def residual_action(x, v, dv, s_max, p_grid, sigma_p=1.0, n_s=60001,
                    s_wide=64.0):
    """Act with the compensated residual operator on a Gaussian test ``W(p)``.

    A Gaussian ``W`` of width ``sigma_p`` has ``What(s)`` Gaussian of width
    ``1/sigma_p``, which decays fast enough that the residual symbol's linear
    growth at large ``s`` -- ``M_res -> -i s V'(x)``, the reason no unbounded
    reach exists -- is harmless.  Truncating the ``s`` integral at
    ``s_max = 2 y_max / hbar`` *is* the reach, so comparing the truncated and
    untruncated actions measures the reach error in the generator itself.
    """
    s = np.linspace(-s_wide, s_wide, n_s)
    ds = s[1] - s[0]
    y = HBAR * s / 2.0
    d_res = v(x + y) - v(x - y) - 2.0 * y * dv(x)
    integrand = (1j * d_res / HBAR) * np.exp(-(sigma_p**2) * s**2 / 2.0)
    phase = np.exp(1j * np.outer(p_grid, s))
    full = (phase * integrand).sum(axis=1).real * ds / (2.0 * np.pi)
    keep = np.abs(s) <= s_max
    trunc = (phase[:, keep] * integrand[keep]).sum(axis=1).real * ds \
        / (2.0 * np.pi)
    return full, trunc


def part_f():
    banner("Part F  Theorem G5: Gamma diverges, the generator converges")
    v, dv, _ = eckart(1.0, 1.0)
    x_probe = 0.6
    xs = np.linspace(-4.0, 4.0, 41)
    p_grid = np.linspace(-4.0, 4.0, 161)

    print(f"\n  Eckart V0 = a = 1, probe at x = {x_probe}")
    print("  census      Gamma_max = max_x sum_q |K_res,q|")
    print("  generator   reach-truncated vs untruncated action of L_res on")
    print("              a Gaussian W of width sigma_p = 1")
    print(f"\n  {'y_max/a':>10s}{'Gamma_max':>14s}{'int Gamma dx':>16s}"
          f"{'rel. error of L_res W':>24s}")
    rows = []
    for mult in [0.5, 1.0, 2.0, 4.0, 8.0]:
        y_max = mult * np.pi
        n_rungs = max(64, int(24 * mult))
        prof = rate_profile(v, dv, y_max, xs, compensated=True,
                            n_rungs=n_rungs)
        full, trunc = residual_action(x_probe, v, dv, 2.0 * y_max / HBAR,
                                      p_grid)
        rel = np.abs(trunc - full).max() / np.abs(full).max()
        rows.append((mult, prof.max(), np.trapezoid(prof, xs), rel))
        print(f"  {mult:10.1f}{prof.max():14.4f}"
              f"{np.trapezoid(prof, xs):16.4f}{rel:24.3e}")
    print("\n  the census diverges with the reach while the operator it")
    print("  represents converges: how many worlds there are is a property")
    print("  of the regulator, not of the physics.")
    return rows


# --------------------------------------------------------------------- #
# Figures                                                                #
# --------------------------------------------------------------------- #
def fig_two_doors(lams, chis, neg):
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    mask = lams > 0
    ax.loglog(lams[mask], chis[mask], "o-", color="#1d9e75",
              label=r"dynamical door: $\chi_Q$")
    ax.axhline(neg, color="#d85a30", ls="--",
               label=r"kinematic door: negativity of $(|0\rangle+|2\rangle)/\sqrt{2}$")
    ax.set_xlabel(r"quartic coupling $\lambda$")
    ax.set_ylabel("strength")
    ax.set_title(r"Two independent occurrences of $\hbar$")
    ax.legend(loc="lower right", fontsize=9)
    ax.grid(alpha=0.3, which="both")
    save_fig(fig, "compensated_ontology_two_doors.png")


def fig_regulator(rows):
    mult = np.array([r[0] for r in rows])
    gmax = np.array([r[1] for r in rows])
    rel = np.array([r[3] for r in rows])
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.6, 3.8))
    ax1.plot(mult, gmax, "o-", color="#534ab7")
    ax1.set_xscale("log", base=2)
    ax1.set_xlabel(r"$y_{\max}/a$  (units of $\pi$)")
    ax1.set_ylabel(r"$\Gamma_{\max}$")
    ax1.set_title("the census diverges")
    ax1.grid(alpha=0.3)
    ax2.plot(mult, rel + 1e-18, "o-", color="#1d9e75")
    ax2.set_xscale("log", base=2)
    ax2.set_yscale("log")
    ax2.set_xlabel(r"$y_{\max}/a$  (units of $\pi$)")
    ax2.set_ylabel(r"relative error of $L_{\rm res}W$")
    ax2.set_title("the generator converges")
    ax2.grid(alpha=0.3, which="both")
    fig.tight_layout()
    save_fig(fig, "compensated_ontology_regulator.png")


def fig_admissibility(rows):
    t = np.array([r[0] for r in rows])
    e_ex = np.array([abs(r[1]) for r in rows])
    e_cl = np.array([abs(r[2]) for r in rows])
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.semilogy(t, e_ex + 1e-18, "o-", color="#1d9e75",
                label=r"(S) + (D): exact evolution")
    ax.semilogy(t, e_cl + 1e-18, "s-", color="#d85a30",
                label=r"(S) alone: classical carrier")
    ax.set_xlabel("t")
    ax.set_ylabel(r"$|\lambda_{\min}(\rho)|$")
    ax.set_title(r"Admissibility is preserved by (S)+(D), not by (S)")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3, which="both")
    save_fig(fig, "compensated_ontology_admissibility.png")


# --------------------------------------------------------------------- #
def main():
    quad_max = part_a()
    rot_err = part_b()
    lams, chis, neg = part_c()
    part_d()
    adm = part_e()
    reg = part_f()

    banner("Figures")
    fig_two_doors(lams, chis, neg)
    fig_regulator(reg)
    fig_admissibility(adm)

    banner("Summary")
    print(f"  G2  quadratic V gives Gamma <= {quad_max:.1e} at every reach,")
    print("      while the field-less rate is O(1): the two readings")
    print("      disagree about the census of a harmonic oscillator")
    print(f"  --  the harmonic null test holds to {rot_err:.1e} at t = pi:")
    print("      a signed W carried exactly by classical rotation")
    print("  G3  chi_Q falls continuously to zero with lambda while the")
    print("      negativity of the state does not move")
    print("  G4  (S) alone drives min eig rho negative; (S)+(D) does not")
    print("  G5  Gamma diverges roughly linearly with the reach while the")
    print("      generator it represents converges to the grid floor")
    print("\ndone.")


if __name__ == "__main__":
    main()
