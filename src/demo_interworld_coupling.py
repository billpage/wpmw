"""
Verification for ``docs/analysis/interworld_coupling.md``.

What the potential looks like when it is read as a coupling between the two
legs of a position pair, and why that reading makes the four rules a
counting statement rather than a postulate.

Write ``x1`` for the ket leg, ``x2`` for the bra leg, ``X = (x1 + x2)/2`` for
the pair midpoint and ``Y = x1 - x2`` for the *full* leg separation.  The
potential enters the von Neumann equation only through

    U(x1, x2) = V(x1) - V(x2) = V(X + Y/2) - V(X - Y/2),

and the Wigner transform is a Fourier transform in ``Y``.  Every claim below
is a consequence of those two sentences.

Note on conventions.  ``docs/analysis/open_position_space.md`` §2 writes the
Wigner kernel with the *half* separation ``y = Y/2``.  This script uses the
full separation ``Y`` throughout, so a coupling written ``sin(k y)`` there
appears as ``sin(k Y / 2)`` here.  The physical content -- each leg is
displaced from the midpoint by half the separation -- is the same.

Parts
-----
A  The difference structure.  The potential term of the von Neumann equation
   is multiplication by ``U``, and ``U`` inherits three properties with no
   further assumptions: it vanishes at coincidence, it is antisymmetric under
   leg exchange, and it vanishes identically when ``V`` does.  An
   antisymmetrised Gaussian pair potential shares the first two and is used
   throughout as the control.
B  Midpoint and separation.  For a single cosine mode ``U`` factorises as a
   midpoint amplitude times a separation grating,
   ``U = 2 V_p sin(kX) sin(kY/2)``, and the grating has *twice* the period of
   ``V`` because each leg moves only ``Y/2``.
C  Fourier duality (Theorem I3).  The momentum-transfer channels available at
   midpoint ``X`` are exactly the Fourier spectrum of ``U(X, .)`` in ``Y``.
   Discrete channels iff the coupling is periodic in ``Y``.  Verified against
   the full Wigner operator at all Moyal orders, for one mode and for three.
D  The Moyal series is the Taylor expansion in ``Y`` (Theorem I4).  Odd powers
   only.  ``Y^1`` is the classical force, so a coupling linear in the leg
   separation is exactly classical -- which is why every quadratic potential,
   harmonic and inverted alike, has no jump channel at all.
E  The coupling winds, it does not push (Theorem I5).  Under the potential
   term alone every ``|rho|`` is constant and ``arg rho`` advances at rate
   ``-U/hbar``, which is the equation of motion for the misalignment ``mu``
   of the phase-alignment layer.  Rule counting closes the part.

Run with ``WPMW_OUTPUT`` set (``/mnt/user-data/outputs`` in the container).
"""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from wpmwlib.wpmw_utils import docs_path, output_path

HBAR = 1.0

# The reference potential: one cosine mode, minimum at the origin.
V_P = 1.5
LAMBDA = 8.0
K_MODE = 2.0 * np.pi / LAMBDA


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


def v_cos(z, amp=V_P, k=K_MODE):
    """Single-mode reference potential, minimum at z = 0."""
    return -amp * np.cos(k * z)


def coupling(vfun, big_x, big_y):
    """U(X, Y) = V(X + Y/2) - V(X - Y/2)."""
    return vfun(big_x + big_y / 2.0) - vfun(big_x - big_y / 2.0)


def theta_exact(w, x, p, vfun):
    """Exact potential term of the Wigner equation, all Moyal orders.

    ``w`` has shape ``(len(p), len(x))``.  The route is the definition: undo
    the Fourier transform in the separation, multiply by ``-i U / hbar``,
    transform back.  No series truncation anywhere.
    """
    dp = p[1] - p[0]
    big_y = 2.0 * np.pi * HBAR * np.fft.fftfreq(len(p), d=dp)
    rho = np.fft.ifft(w, axis=0)
    xx, yy = np.meshgrid(x, big_y, indexing="xy")
    return np.real(np.fft.fft((-1j / HBAR) * coupling(vfun, xx, yy) * rho,
                              axis=0))


def dp_spectral(w, p, order=1):
    """Spectrally exact momentum derivative of ``w`` along axis 0."""
    kp = 2.0 * np.pi * np.fft.fftfreq(len(p), d=p[1] - p[0])
    return np.real(np.fft.ifft(((1j * kp[:, None]) ** order)
                               * np.fft.fft(w, axis=0), axis=0))


def test_state(x, p):
    """A test Wigner array with genuine negativity and no symmetry to exploit."""
    xx, pp = np.meshgrid(x, p)
    return (np.exp(-xx ** 2 - pp ** 2 / 0.5) * np.cos(3.0 * pp)
            + 0.3 * np.exp(-(xx - 1.2) ** 2 - (pp + 0.4) ** 2 / 0.3))


def y_spectrum(gfun, y_max=200.0, n=1 << 16):
    """Normalised Fourier magnitude of a separation profile, and its axis.

    Returns ``(q, mag)`` where ``q`` is the conjugate momentum-transfer axis.
    A Gaussian apodisation of width ``0.45 y_max`` suppresses the truncation
    ringing that would otherwise be mistaken for structure; it also gives
    every line a finite width proportional to ``1 / y_max``, which is what
    ``line_widths`` exploits.
    """
    yg = np.linspace(-y_max, y_max, n, endpoint=False)
    win = np.exp(-(yg / (0.45 * y_max)) ** 2)
    mag = np.abs(np.fft.fftshift(np.fft.fft(gfun(yg) * win)))
    q = np.fft.fftshift(np.fft.fftfreq(n, d=yg[1] - yg[0])) * 2.0 * np.pi * HBAR
    return q, mag / mag.max()


def find_lines(gfun, thresh=0.5, **kw):
    """Locate the connected runs of the Y-spectrum above ``thresh``.

    Returns ``(peaks, widths)``: the ``q`` of each run's maximum, and each
    run's extent in ``q``.  Counting *runs* rather than samples is what makes
    the count independent of the transform grid.
    """
    q, mag = y_spectrum(gfun, **kw)
    above = mag > thresh
    edges = np.flatnonzero(np.diff(above.astype(int)))
    starts = edges[::2] + 1
    stops = edges[1::2] + 1
    peaks, widths = [], []
    for lo, hi in zip(starts, stops):
        peaks.append(q[lo:hi][np.argmax(mag[lo:hi])])
        widths.append(q[hi - 1] - q[lo])
    return np.array(peaks), np.array(widths)


def width_scaling(gfun, windows=(50.0, 100.0, 200.0, 400.0)):
    """Half-max width of the strongest line as the Y-window is widened.

    A genuinely discrete spectrum has zero intrinsic width, so the measured
    width is set entirely by the window and falls like ``1 / y_max``.  A
    continuous spectrum has an intrinsic width that the window cannot reduce.
    """
    out = []
    for y_max in windows:
        _, widths = find_lines(gfun, y_max=y_max)
        out.append((y_max, widths.max() if len(widths) else float("nan")))
    return out


# Control coupling: antisymmetrised Gaussian.  It satisfies U(X, 0) = 0 and
# U(X, -Y) = -U(X, Y), so it is a fair test of whether those two properties
# are on their own enough to produce a finite rule set.
def gauss_anti(big_y, width=1.3, offset=1.2):
    return (np.exp(-(big_y - offset) ** 2 / width)
            - np.exp(-(big_y + offset) ** 2 / width))


# --------------------------------------------------------------------- A


def part_a():
    banner("A. The potential enters only as a difference")
    rng = np.random.default_rng(20260813)
    x1 = rng.uniform(-25.0, 25.0, 40000)
    x2 = rng.uniform(-25.0, 25.0, 40000)
    big_x, big_y = (x1 + x2) / 2.0, x1 - x2

    print("  U(x1, x2) = V(x1) - V(x2), so with no further assumptions:")
    coincide = np.abs(coupling(v_cos, big_x, np.zeros_like(big_y))).max()
    exch = np.abs(coupling(v_cos, big_x, big_y)
                  + coupling(v_cos, big_x, -big_y)).max()
    free = np.abs(coupling(lambda z: np.zeros_like(z), big_x, big_y)).max()
    print(f"    (i)   coincident legs, max|U(X, 0)|            = {coincide:.3e}")
    print(f"    (ii)  exchange, max|U(X, Y) + U(X, -Y)|        = {exch:.3e}")
    print(f"    (iii) free legs, max|U| with V = 0             = {free:.3e}")

    print("\n  The control coupling shares (i) and (ii):")
    print(f"    antisym Gaussian, |g(0)|                       = "
          f"{abs(gauss_anti(np.array([0.0]))[0]):.3e}")
    yy = rng.uniform(-8.0, 8.0, 20000)
    print(f"    antisym Gaussian, max|g(Y) + g(-Y)|            = "
          f"{np.abs(gauss_anti(yy) + gauss_anti(-yy)).max():.3e}")
    print("  -> (i) and (ii) do not distinguish the two.  Part C does.")
    return coincide, exch, free


# --------------------------------------------------------------------- B


def part_b():
    banner("B. Midpoint and separation: the grating and its doubled period")
    rng = np.random.default_rng(7)
    big_x = rng.uniform(-25.0, 25.0, 40000)
    big_y = rng.uniform(-25.0, 25.0, 40000)
    lhs = coupling(v_cos, big_x, big_y)
    rhs = 2.0 * V_P * np.sin(K_MODE * big_x) * np.sin(K_MODE * big_y / 2.0)
    err = np.abs(lhs - rhs).max()
    print("  identity   V(X+Y/2) - V(X-Y/2) = 2 V_p sin(kX) sin(kY/2)")
    print(f"    max error                                      = {err:.3e}")

    # Period in Y, measured rather than asserted.
    yg = np.linspace(0.0, 6.0 * LAMBDA, 200001)
    prof = coupling(v_cos, 1.0, yg)
    zeros = yg[np.where(np.sign(prof[:-1]) != np.sign(prof[1:]))[0]]
    period = 2.0 * np.median(np.diff(zeros))
    print(f"    potential period lambda                        = {LAMBDA:.4f}")
    print(f"    coupling period in Y (measured from zeros)     = {period:.4f}")
    print(f"    ratio                                          = "
          f"{period / LAMBDA:.6f}")
    print("  -> each leg moves Y/2, so the difference completes a cycle only")
    print("     after the legs separate by two potential wavelengths.")

    # The midpoint amplitude is the force, up to the fixed factor hbar k.
    xg = np.linspace(-LAMBDA, LAMBDA, 2001)
    gamma = (V_P / HBAR) * np.sin(K_MODE * xg)
    force = -V_P * K_MODE * np.sin(K_MODE * xg)
    print(f"\n    max|Gamma(X) + F(X)/(hbar k)|                  = "
          f"{np.abs(gamma + force / (HBAR * K_MODE)).max():.3e}")
    print("  -> the midpoint amplitude is the local classical force.")
    return err, period


# --------------------------------------------------------------------- C


def part_c():
    banner("C. Theorem I3: channels are the Y-spectrum of the coupling")

    # C1: one mode.  Two lines, at +- hbar k / 2.
    peaks, _ = find_lines(lambda y: np.sin(K_MODE * y / 2.0))
    print(f"  cosine mode: lines found                         = {len(peaks)}")
    print("    at q = " + ", ".join(f"{v:+.5f}" for v in peaks))
    print(f"    hbar k / 2                                     = "
          f"{HBAR * K_MODE / 2.0:+.5f}")
    print(f"    max|peak| - hbar k / 2                         = "
          f"{np.abs(np.abs(peaks) - HBAR * K_MODE / 2.0).max():.3e}")

    # C2: discrete or continuous?  Widen the window and watch the line width.
    print("\n  is the spectrum discrete?  widen the Y-window and remeasure:")
    print("    Y window     cosine mode      antisym Gaussian")
    cos_scale = width_scaling(lambda y: np.sin(K_MODE * y / 2.0))
    gau_scale = width_scaling(gauss_anti)
    for (ym, wc), (_, wg) in zip(cos_scale, gau_scale):
        print(f"    {ym:8.0f}     {wc:.5f}          {wg:.5f}")
    print("    -> the cosine line narrows like 1/Y_max: zero intrinsic width.")
    print("       the Gaussian control saturates: genuine continuum.")

    # C3: the two lines are the whole operator, exactly.
    nx, np_, n_shift = 64, 128, 4
    x = -LAMBDA / 2.0 + LAMBDA * np.arange(nx) / nx
    dp = HBAR * K_MODE / 2.0 / n_shift
    p = (np.arange(np_) - np_ // 2) * dp
    w = test_state(x, p)
    interior = slice(3 * n_shift, np_ - 3 * n_shift)

    th = theta_exact(w, x, p, v_cos)
    gamma = (V_P / HBAR) * np.sin(K_MODE * x)
    stencil = gamma[None, :] * (np.roll(w, -n_shift, axis=0)
                                - np.roll(w, n_shift, axis=0))
    err1 = np.abs(th[interior] - stencil[interior]).max()
    print(f"\n  one mode: max|Theta_exact - two-line stencil|    = {err1:.3e}"
          f"    (max|Theta| = {np.abs(th).max():.3e})")

    # C4: three modes.  Each contributes its own pair of lines, at its own
    # shift.  Nothing mixes: the channel set is the union.
    amps = {1: V_P, 2: 0.6, 3: -0.35}

    def v_multi(z):
        return sum(-a * np.cos(m * K_MODE * z) for m, a in amps.items())

    th3 = theta_exact(w, x, p, v_multi)
    pred = np.zeros_like(th3)
    for m, a in amps.items():
        g_m = (a / HBAR) * np.sin(m * K_MODE * x)
        pred += g_m[None, :] * (np.roll(w, -m * n_shift, axis=0)
                                - np.roll(w, m * n_shift, axis=0))
    interior3 = slice(5 * n_shift, np_ - 5 * n_shift)
    err3 = np.abs(th3[interior3] - pred[interior3]).max()
    print(f"  three modes: max|Theta_exact - union of pairs|   = {err3:.3e}"
          f"    (max|Theta| = {np.abs(th3).max():.3e})")
    print("  -> modes do not mix.  M modes give 2M shifts, hence 4M rules.")
    return err1, err3, peaks, cos_scale, gau_scale


# --------------------------------------------------------------------- D


def part_d():
    banner("D. Theorem I4: the Moyal series is the Y-expansion of U")
    nx, np_, dp = 128, 256, 0.05
    x = -20.0 + 40.0 * np.arange(nx) / nx
    p = (np.arange(np_) - np_ // 2) * dp
    xx, pp = np.meshgrid(x, p)
    # A cat-like state, so that the test has real negativity to propagate.
    w = (np.exp(-(xx - 2.0) ** 2 - pp ** 2 / 0.25)
         + np.exp(-(xx + 2.0) ** 2 - pp ** 2 / 0.25)
         + 2.0 * np.exp(-xx ** 2 - pp ** 2 / 0.25) * np.cos(4.0 * pp / HBAR))
    print(f"  test state min W = {w.min():+.4f}, max W = {w.max():+.4f}")
    d1 = dp_spectral(w, p, 1)
    d3 = dp_spectral(w, p, 3)

    rows = []
    # Linear: U = a Y exactly.
    th = theta_exact(w, x, p, lambda z: 0.7 * z)
    rows.append(("linear   V = 0.7 x", "Y^1",
                 np.abs(th - 0.7 * d1).max() / np.abs(th).max()))
    # Quadratic, both signs: U = c X Y exactly.  No Y^3 term exists.
    for label, c in [("harmonic V = +x^2/2", 1.0),
                     ("inverted V = -x^2/2", -1.0)]:
        th = theta_exact(w, x, p, lambda z, c=c: 0.5 * c * z ** 2)
        rows.append((label, "Y^1",
                     np.abs(th - c * xx * d1).max() / np.abs(th).max()))
    # Cubic: U = b(3 X^2 Y + Y^3/4), so the series stops after one correction.
    b = 0.02
    th = theta_exact(w, x, p, lambda z: b * z ** 3)
    cl = 3.0 * b * xx ** 2 * d1
    q1 = -(HBAR ** 2 / 24.0) * (6.0 * b) * d3
    rows.append(("cubic    V = 0.02 x^3", "Y^1 only",
                 np.abs(th - cl).max() / np.abs(th).max()))
    rows.append(("cubic    V = 0.02 x^3", "Y^1 + Y^3",
                 np.abs(th - cl - q1).max() / np.abs(th).max()))

    print("\n  potential              truncation   rel. residual")
    for label, trunc, res in rows:
        print(f"  {label:22s} {trunc:11s}  {res:.3e}")
    print("\n  -> a coupling linear in the leg separation is exactly classical.")
    print("     U = c X Y for every quadratic V, which is why the harmonic and")
    print("     the inverted harmonic both have no jump channel whatsoever.")
    return rows


# --------------------------------------------------------------------- E


def part_e():
    banner("E. Theorem I5: the coupling winds the pair, it does not push it")
    # Potential-only evolution of rho in the (x1, x2) representation is
    # multiplication by exp(-i U t / hbar).  Check that directly on a
    # normalised random density matrix.
    n = 64
    grid = np.linspace(-6.0, 6.0, n)
    xa, xb = np.meshgrid(grid, grid, indexing="ij")
    rng = np.random.default_rng(3)
    psi = (rng.standard_normal(n) + 1j * rng.standard_normal(n))
    psi *= np.exp(-grid ** 2 / 8.0)
    rho0 = np.outer(psi, psi.conj())
    u = v_cos(xa) - v_cos(xb)

    dt, steps = 1e-3, 400
    rho = rho0.copy()
    for _ in range(steps):
        rho = rho * np.exp(-1j * u * dt / HBAR)
    t = dt * steps

    mod_err = np.abs(np.abs(rho) - np.abs(rho0)).max()
    live = np.abs(rho0) > 1e-9
    dmu = np.angle(rho[live] * np.conj(rho0[live]))
    pred = np.angle(np.exp(-1j * u[live] * t / HBAR))
    ang_err = np.abs(np.angle(np.exp(1j * (dmu - pred)))).max()
    print(f"  after t = {t:.3f} under the potential term alone:")
    print(f"    max change in |rho|                            = {mod_err:.3e}")
    print(f"    max|arg advance - (-U t / hbar)|               = {ang_err:.3e}")
    print("  -> |rho| is untouched; the coupling is a winding rate for")
    print("     mu = arg rho, d mu / d t = -U / hbar.  No force, no work.")

    banner("Rule counting")
    print("  potential                          lines   rules")
    x_ref = 1.0
    amps = {1: V_P, 2: 0.6, 3: 0.35}

    def sep_profile(modes):
        def g(y):
            return sum(amps[m] * np.sin(m * K_MODE * x_ref)
                       * np.sin(m * K_MODE * y / 2.0) for m in modes)
        return g

    cases = [("one cosine mode", sep_profile([1])),
             ("three cosine modes", sep_profile([1, 2, 3]))]
    counts = []
    for label, g in cases:
        peaks, _ = find_lines(g, thresh=0.15)
        counts.append((label, len(peaks), 2 * len(peaks)))
        print(f"  {label:34s} {len(peaks):5d}   {2 * len(peaks):5d}")
    print(f"  {'antisym Gaussian (control)':34s}"
          f"{'  cont.':>7s}   continuum")
    counts.append(("antisym Gaussian (control)", -1, -1))
    print("\n  Each line is a signed shift channel, and each signed channel")
    print("  splits into a forward and a reverse rule, so a one-mode potential")
    print("  has exactly four.")
    return mod_err, ang_err, counts


# --------------------------------------------------------------------- figure


def fig_coupling():
    fig, axes = plt.subplots(2, 2, figsize=(13.0, 8.2))

    # (a) a pair straddling the potential
    ax = axes[0, 0]
    xg = np.linspace(-1.5 * LAMBDA, 1.5 * LAMBDA, 900)
    ax.plot(xg, v_cos(xg), color="0.35", lw=2)
    x_mid, y_sep = 2.0, 5.2
    x1, x2 = x_mid + y_sep / 2.0, x_mid - y_sep / 2.0
    ax.plot([x1], [v_cos(x1)], "o", ms=11, color="C3", zorder=5)
    ax.plot([x2], [v_cos(x2)], "o", ms=11, color="C0", zorder=5)
    ax.plot([x_mid], [-1.75 * V_P], "x", ms=10, color="k", mew=2)
    ax.vlines([x1], v_cos(x1), 1.55 * V_P, color="C3", ls=":", lw=1.2)
    ax.vlines([x2], v_cos(x2), 1.55 * V_P, color="C0", ls=":", lw=1.2)
    ax.vlines([x_mid], -1.75 * V_P, 1.55 * V_P, color="k", ls=":", lw=1.0)
    ax.annotate("", xy=(x1, 1.55 * V_P), xytext=(x2, 1.55 * V_P),
                arrowprops=dict(arrowstyle="<->", lw=1.6, color="k"))
    ax.text(x_mid, 1.72 * V_P, r"separation $Y = x_1 - x_2$",
            ha="center", fontsize=9.5)
    ax.text(x2, 2.15 * V_P, r"$x_2 = X - Y/2$", ha="center",
            fontsize=10, color="C0")
    ax.text(x1, 2.15 * V_P, r"$x_1 = X + Y/2$", ha="center",
            fontsize=10, color="C3")
    ax.text(x_mid, -2.0 * V_P, r"midpoint $X$", ha="center", va="top",
            fontsize=9.5)
    ax.hlines([v_cos(x1)], x1, 1.28 * LAMBDA, color="C3", lw=0.9, ls="--")
    ax.hlines([v_cos(x2)], x2, 1.28 * LAMBDA, color="C0", lw=0.9, ls="--")
    ax.annotate("", xy=(1.28 * LAMBDA, v_cos(x1)),
                xytext=(1.28 * LAMBDA, v_cos(x2)),
                arrowprops=dict(arrowstyle="<->", lw=2.2, color="C2"))
    ax.text(1.34 * LAMBDA, 0.0,
            "$U = V(x_1) - V(x_2)$\nthe only way $V$ enters",
            fontsize=9.5, color="C2", va="center")
    ax.set_xlim(-1.5 * LAMBDA, 2.2 * LAMBDA)
    ax.set_ylim(-2.7 * V_P, 2.6 * V_P)
    ax.set_yticks([])
    ax.set_xlabel("position")
    ax.set_title("(a) a world-pair straddling the potential", fontsize=10.5)

    # (b) the separation grating and its doubled period
    ax = axes[0, 1]
    yg = np.linspace(-2.6 * LAMBDA, 2.6 * LAMBDA, 1400)
    ax.plot(yg, coupling(v_cos, x_mid, yg), color="C2", lw=2,
            label=r"$U(X,Y) = 2V_p\sin kX\,\sin(kY/2)$")
    ax.plot(yg, v_cos(yg) - v_cos(0.0), color="0.6", lw=1.2, ls="--",
            label=r"the potential itself, period $\lambda$")
    ax.axhline(0.0, color="k", lw=0.6)
    ax.axvline(0.0, color="k", lw=0.6)
    for nper in (-2, -1, 1, 2):
        ax.axvline(nper * 2 * LAMBDA, color="C2", ls=":", lw=1)
    ax.annotate("", xy=(0.0, 2.45 * V_P), xytext=(2 * LAMBDA, 2.45 * V_P),
                arrowprops=dict(arrowstyle="<->", lw=1.6, color="C2"))
    ax.text(LAMBDA, 2.58 * V_P,
            r"period $2\lambda$  $\Rightarrow$  wavevector $k/2$",
            ha="center", fontsize=10, color="C2")
    ax.set_xlabel(r"leg separation $Y$")
    ax.set_ylim(-2.5 * V_P, 3.05 * V_P)
    ax.legend(fontsize=8, loc="lower center", framealpha=0.95)
    ax.set_title("(b) the coupling is periodic in the separation,\n"
                 "with twice the period of the potential", fontsize=10.5)

    # (c) the Y-spectrum is the channel set
    ax = axes[1, 0]
    q, mag = y_spectrum(lambda y: np.sin(K_MODE * y / 2.0))
    ax.plot(q, mag, color="C2", lw=1.8,
            label=r"$U$ from $V(x_1)-V(x_2)$: two lines")
    q2, mag2 = y_spectrum(gauss_anti)
    ax.plot(q2, mag2, color="C1", lw=1.8,
            label="antisymmetric Gaussian:\ncontinuum")
    for sgn in (-1, 1):
        ax.axvline(sgn * HBAR * K_MODE / 2.0, color="0.5", ls=":", lw=1.2)
    ax.set_xlim(-3.0, 3.0)
    ax.set_xlabel(r"momentum transfer $q$, conjugate to $Y$")
    ax.set_ylabel("channel weight")
    ax.text(HBAR * K_MODE / 2.0 + 0.10, 0.72, r"$+\hbar k/2$", fontsize=9)
    ax.text(-HBAR * K_MODE / 2.0 - 0.10, 0.72, r"$-\hbar k/2$", fontsize=9,
            ha="right")
    ax.legend(fontsize=8, loc="upper left")
    ax.set_title("(c) channels are the $Y$-spectrum of $U$;\n"
                 r"periodic $\Rightarrow$ discrete (Theorem I3)", fontsize=10.5)

    # (d) the U(X, Y) landscape
    ax = axes[1, 1]
    xs = np.linspace(-1.5 * LAMBDA, 1.5 * LAMBDA, 400)
    ys = np.linspace(-2.0 * LAMBDA, 2.0 * LAMBDA, 400)
    xx, yy = np.meshgrid(xs, ys)
    uu = coupling(v_cos, xx, yy)
    im = ax.pcolormesh(xx, yy, uu, cmap="RdBu_r", shading="auto",
                       vmin=-2 * V_P, vmax=2 * V_P)
    ax.contour(xx, yy, uu, levels=[0.0], colors="k", linewidths=0.8)
    ax.axhline(0.0, color="k", lw=2)
    ax.text(0.0, 0.09 * LAMBDA, r"$Y = 0$: legs coincide, $U \equiv 0$",
            fontsize=9, ha="center")
    ax.set_xlabel(r"pair midpoint $X$")
    ax.set_ylabel(r"leg separation $Y$")
    ax.set_title(r"(d) $U(X,Y)$: odd in $Y$, odd in $X$ about each well",
                 fontsize=10.5)
    fig.colorbar(im, ax=ax, fraction=0.046)

    fig.tight_layout()
    save_fig(fig, "interworld_coupling.png")
    plt.close(fig)


def main():
    print("Verification for docs/analysis/interworld_coupling.md")
    part_a()
    part_b()
    part_c()
    part_d()
    part_e()
    banner("Figures")
    fig_coupling()
    print("\ndone.")


if __name__ == "__main__":
    main()
