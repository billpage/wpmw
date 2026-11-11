"""
Companion to docs/supplement/zero_gap_superposition.md: a superposition of two
raised-cosine pulses with an exactly empty gap between them (D. Cyganski's
example of 10 November 2026), its Wigner function, and how the gap fills.

Units hbar = m = a = 1, where a is the half-width of one pulse.  The pulses are
centred at -d and +d with d = 2.2, and the right-hand one carries the phase
theta = pi/2.

Parts:
  A   Proposition ZG1: the Wigner function in closed form.  Checked against
      direct quadrature; the band between the lobes, the empty strips, the
      marginal in the gap, the purity and the L1 norm.
  B   Proposition ZG2: the exact free solution as a sum of twelve Moshinsky
      functions.  Checked against 40-digit arithmetic and against a spectral
      (FFT) propagation, whose error floor is reported.
  C   Proposition ZG3: how the gap fills.  The edge-wave law
      rho ~ J^2 t^5 / (2 pi s^6) deep in the gap, the boundary layer at an
      edge, the leaked probability J^2 t^(5/2) / (20 sqrt(pi)), and the far
      tail ("instantaneous spreading").
  D   Proposition ZG4: fringes and nodes.  Fringe spacing early and late,
      zeros of psi in the (x, t) plane found by phase winding, one of them
      refined at 40 digits, and the permanent node at theta = pi.
  E   Proposition ZG5: the Wigner picture.  The flow is a shear; the density in
      the gap is the lobes' momentum tails plus the sheared band; the signed
      weight in the gap.
  F   Proposition ZG6: signed worlds.  With V = 0 the model has no events.
      N worlds sampled from |W| stream ballistically; the signed histogram is
      compared with the exact density, and the ensemble needed to see the gap
      fill is computed.
  G   Remark ZG7: a harmonic trap.  The flow is a rotation, still with no
      events: exact nodes at a quarter period, the gap again at a half period.

Every number in the note comes from this script.  Run with
  WPMW_OUTPUT=<dir> PYTHONPATH=src python3 -u src/demo_zero_gap_superposition.py
Add --no-animation to skip the animation (about a minute of the run time).
"""
import json
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np
from matplotlib import animation
from scipy.special import wofz

from wpmwlib.wpmw_utils import docs_path, output_path

A_HALF = 1.0                     # a: half-width of one pulse
D = 2.2                          # d: half-separation of the pulse centres
THETA = np.pi/2                  # relative phase, on the right-hand pulse
AMP = 2/np.sqrt(3*A_HALF)        # A: normalisation of one pulse
KAP = np.pi/A_HALF               # kappa = pi/a
JUMP = 0.5*AMP*KAP**2            # J: jump of psi_0'' at an edge
GAP = D - A_HALF                 # the gap is |x| < GAP
P_MAX = 40.0                     # momentum cut for sampling worlds
results = {}


def save(fig, name):
    fig.savefig(output_path(name), dpi=130, bbox_inches="tight")
    dp = docs_path(name)
    if dp:
        fig.savefig(dp, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print("  wrote", name)


# ----------------------------------------------------------------- the state
def psi0(x):
    """One raised-cosine pulse, centred at 0, zero outside |x| < a."""
    return np.where(np.abs(x) < A_HALF, AMP*np.cos(0.5*KAP*x)**2, 0.0)


def psi_init(x, theta=THETA):
    return (psi0(x + D) + np.exp(1j*theta)*psi0(x - D))/np.sqrt(2)


def _S(k, L):
    """int_{-L}^{L} exp(i k y) dy."""
    return 2*L*np.sinc(k*L/np.pi)


def wigner0(x, p):
    """Wigner function of one pulse, closed form (Proposition ZG1)."""
    x, p = np.broadcast_arrays(np.asarray(x, float), np.asarray(p, float))
    L = np.clip(A_HALF - np.abs(x), 0.0, None)
    w = ((1 + 0.5*np.cos(2*KAP*x))*_S(2*p, L)
         + np.cos(KAP*x)*(_S(2*p + KAP, L) + _S(2*p - KAP, L))
         + 0.25*(_S(2*p + 2*KAP, L) + _S(2*p - 2*KAP, L)))
    return AMP**2/(4*np.pi)*w


def wigner(x, p, theta=THETA):
    """Wigner function of the superposition at t = 0."""
    return (0.5*(wigner0(x - D, p) + wigner0(x + D, p))
            + wigner0(x, p)*np.cos(2*p*D - theta))


def moshinsky(x, k, t):
    """Free evolution of exp(i k x) Theta(-x) (Moshinsky 1952)."""
    z = (x - k*t)*np.exp(-0.25j*np.pi)/np.sqrt(2*t)
    return 0.5*np.exp(0.5j*x*x/t)*wofz(1j*z)


def pulse(x, t):
    """Free evolution of one pulse centred at 0: six Moshinsky functions."""
    out = 0
    for c, k in ((0.5*AMP, 0.0), (0.25*AMP, KAP), (0.25*AMP, -KAP)):
        out = out + c*(np.exp(1j*k*A_HALF)*moshinsky(x - A_HALF, k, t)
                       - np.exp(-1j*k*A_HALF)*moshinsky(x + A_HALF, k, t))
    return out


def psi_exact(x, t, theta=THETA):
    x = np.asarray(x, float)
    if t == 0:
        return psi_init(x, theta).astype(complex)
    return (pulse(x + D, t) + np.exp(1j*theta)*pulse(x - D, t))/np.sqrt(2)


EDGES = ((-D - A_HALF, +1, 1.0), (-D + A_HALF, -1, 1.0),
         (D - A_HALF, +1, np.exp(1j*THETA)), (D + A_HALF, -1, np.exp(1j*THETA)))


def psi_edge(x, t, order=1):
    """Edge-wave asymptotics (Proposition ZG3); order 1 adds -6 i t / s^2."""
    tot = 0
    for xe, eta, c in EDGES:
        corr = 1 - 6j*t/(xe - x)**2 if order else 1.0
        tot = tot + c*eta*np.exp(0.5j*(x - xe)**2/t)/(xe - x)**3*corr
    return (2j*np.pi*t)**-0.5*(-1j*JUMP*t**3)*tot/np.sqrt(2)


def phi0_sq(p):
    """|phi_0(p)|^2, the momentum density of one pulse."""
    p = np.asarray(p, float)
    den = p*(KAP**2 - p**2)
    with np.errstate(divide="ignore", invalid="ignore"):
        f = KAP**2*np.sin(p*A_HALF)/den
    f = np.where(np.abs(p) < 1e-9, A_HALF, f)
    f = np.where(np.abs(np.abs(p) - KAP) < 1e-9, 0.5*A_HALF, f)
    return AMP**2*f**2/(2*np.pi)


def phi_sq(p, theta=THETA):
    return phi0_sq(p)*(1 + np.cos(2*p*D - theta))


def _mp_psi(x, t):
    """psi(x, t) in mpmath arithmetic, for checks."""
    a, d, kap = mp.mpf(1), mp.mpf("2.2"), mp.pi
    amp = 2/mp.sqrt(3)

    def mosh(y, k):
        return mp.exp(1j*(k*y - k*k*t/2))*mp.erfc((y - k*t)/mp.sqrt(2j*t))/2

    def one(y):
        return sum(c*(mp.exp(1j*k*a)*mosh(y - a, k) - mp.exp(-1j*k*a)*mosh(y + a, k))
                   for c, k in ((amp/2, mp.mpf(0)), (amp/4, kap), (amp/4, -kap)))

    return (one(x + d) + 1j*one(x - d))/mp.sqrt(2)


# ------------------------------------------------------------------- Part A
def part_A():
    print("=== Part A: the Wigner function in closed form (Proposition ZG1) ===")
    y = np.linspace(-3.3, 3.3, 6601)
    dy = y[1] - y[0]
    xs = np.array([-2.9, -2.2, -1.6, -0.7, 0.0, 0.31, 0.95, 2.2, 2.65])
    ps = np.linspace(-9, 9, 181)
    g = np.conj(psi_init(xs[:, None] + y[None, :]))*psi_init(xs[:, None] - y[None, :])
    quad = (g @ np.exp(2j*np.outer(ps, y)).T).real*dy/np.pi
    err = np.abs(quad - wigner(xs[:, None], ps[None, :])).max()
    print(f"  closed form against direct quadrature, max difference = {err:.1e}")
    x = np.linspace(-3.5, 3.5, 1401)
    p = np.linspace(-12, 12, 2401)
    W = wigner(x[:, None], p[None, :])
    print(f"  pi W at a lobe centre = {np.pi*wigner(D, 0.0):.4f};  "
          f"pi W(0, 0) = {np.pi*wigner(0.0, 0.0):+.1e}")
    pf = np.pi/(4*D)
    print(f"  first bright fringe p = pi/(4d) = {pf:.4f}: pi W = {np.pi*wigner(0.0, pf):.4f};"
          f"  fringe period pi/d = {np.pi/D:.4f}")
    print(f"  extreme values of pi W on the grid: {np.pi*W.min():+.4f}, {np.pi*W.max():+.4f}"
          "  (bound 1)")
    strip = (np.abs(x) > A_HALF + 1e-9) & (np.abs(x) < GAP - 1e-9)
    print(f"  max |W| in the strips a < |x| < d - a: {np.abs(W[strip]).max():.1e}")
    pw = np.linspace(-400, 400, 160001)
    dpw = pw[1] - pw[0]
    xg = np.linspace(-0.95, 0.95, 39)
    Wg = wigner(xg[:, None], pw[None, :])
    marg = np.abs(Wg.sum(1)*dpw).max()
    absm = (np.abs(Wg).sum(1)*dpw).max()
    print(f"  in the gap: max |int W dp| = {marg:.1e},  max int |W| dp = {absm:.3f}")
    x2 = np.linspace(-3.2, 3.2, 1281)
    p2 = np.linspace(-150, 150, 60001)
    cell = (x2[1] - x2[0])*(p2[1] - p2[0])
    norm = pur = l1 = l1_one = pur_mix = l1_in = 0.0
    for i in range(0, len(x2), 64):
        Wc = wigner(x2[i:i+64, None], p2[None, :])
        lobes = 0.5*(wigner0(x2[i:i+64, None] - D, p2[None, :])
                     + wigner0(x2[i:i+64, None] + D, p2[None, :]))
        norm += Wc.sum()*cell
        pur += 2*np.pi*(Wc**2).sum()*cell
        pur_mix += 2*np.pi*(lobes**2).sum()*cell
        l1 += np.abs(Wc).sum()*cell
        l1_in += np.abs(Wc[:, np.abs(p2) <= P_MAX]).sum()*cell
        l1_one += 2*np.abs(lobes).sum()*cell/2
    print(f"  int W = {norm:.5f};  purity 2 pi int W^2 = {pur:.4f}"
          f"  (two lobes without the band: {pur_mix:.4f})")
    print(f"  L1 norm Z = int |W| = {l1:.4f}  (lobes alone: {l1_one:.4f});"
          f"  part with |p| > {P_MAX:.0f}: {l1 - l1_in:.1e}")
    print(f"  Z / Z_lobes = {l1/l1_one:.5f};  1 + 2/pi = {1 + 2/np.pi:.5f}  (Proposition R2)")
    results["A"] = dict(quad_err=err, Z=l1, Z_lobes=l1_one, purity=pur,
                        purity_mixture=pur_mix, gap_abs_marginal=absm,
                        Z_tail=l1 - l1_in)
    return l1


# ------------------------------------------------------------------- Part B
def part_B():
    print("\n=== Part B: the exact free solution (Proposition ZG2) ===")
    mp.mp.dps = 40
    print("  closed form in double precision against 40-digit arithmetic:")
    print("      x       t       |psi|^2 (40 digits)    relative error of psi")
    worst = 0.0
    for x, t in ((0.0, 0.02), (0.6, 0.1), (0.0, 0.3), (2.2, 1.0), (10.0, 0.01),
                 (100.0, 0.01), (0.0, 0.001), (0.3, 0.0001)):
        e = complex(_mp_psi(mp.mpf(x), mp.mpf(t)))
        n = complex(psi_exact(np.array([x]), t)[0])
        rel = abs(e - n)/abs(e)
        if abs(e) > 1e-7:
            worst = max(worst, rel)
        print(f"   {x:6.1f}  {t:7.4f}   {abs(e)**2:.6e}          {rel:.1e}")
    print(f"  worst relative error where |psi| > 1e-7: {worst:.1e};"
          " below that the twelve terms cancel to rounding")
    L, N = 400.0, 2**18
    xs = (np.arange(N) - N//2)*(L/N)
    ks = 2*np.pi*np.fft.fftfreq(N, L/N)
    f0 = psi_init(xs)
    print("  spectral propagation (2^18 points, box 400) against the closed form:")
    floor = 0.0
    for t in (0.02, 0.1, 0.3, 1.0, 3.0):
        ft = np.fft.ifft(np.exp(-0.5j*ks**2*t)*np.fft.fft(f0))
        m = np.abs(xs) < 12
        ex = psi_exact(xs[m], t)
        e = np.abs(ft[m] - ex).max()
        floor = max(floor, e)
        nrm = np.sum(np.abs(psi_exact(xs, t))**2)*(L/N)
        print(f"    t = {t:4.2f}: max |psi_fft - psi| = {e:.1e};  norm of closed form = {nrm:.8f}")
    results["B"] = dict(worst_rel=worst, fft_floor=floor)


# ------------------------------------------------------------------- Part C
def gap_probability(t, n=96001):
    xx = np.linspace(-GAP, GAP, n)
    return np.trapezoid(np.abs(psi_exact(xx, t))**2, xx)


def part_C():
    print("\n=== Part C: how the gap fills (Proposition ZG3) ===")
    print(f"  J = psi_0''(edge) = {JUMP:.5f},  J^2 = pi^4/3 = {JUMP**2:.5f}")
    x = np.array([0.0, 0.3, 0.6])
    print("  edge-wave law against the exact solution (relative error of psi):")
    print("      t       rho(0,t) exact    leading: x=0   0.3    0.6     with -6it/s^2: 0   0.3    0.6")
    for t in (0.005, 0.01, 0.02, 0.05):
        e = psi_exact(x, t)
        e0 = np.abs(psi_edge(x, t, 0) - e)/np.abs(e)
        e1 = np.abs(psi_edge(x, t, 1) - e)/np.abs(e)
        print(f"    {t:6.3f}   {abs(e[0])**2:.4e}      {e0[0]:.3f}  {e0[1]:.3f}  {e0[2]:.3f}"
              f"            {e1[0]:.4f} {e1[1]:.4f} {e1[2]:.4f}")
    ts = np.array([0.005, 0.01, 0.02, 0.04])
    r0 = np.array([abs(psi_exact(np.array([0.0]), t)[0])**2 for t in ts])
    sl = np.diff(np.log(r0))/np.diff(np.log(ts))
    print(f"  local exponent d ln rho(0,t) / d ln t between {ts}: {np.round(sl, 3)}")
    lead = JUMP**2*ts[0]**5/(2*np.pi*GAP**6)
    print(f"  rho(0,t) / [J^2 t^5 (1 + cos theta) / (2 pi (d-a)^6)] at t = {ts[0]}: "
          f"{r0[0]/(lead*(1 + np.cos(THETA))):.4f}  (two inner edges only)")
    xe = D - A_HALF
    print("  at the former edge x = d - a:  rho / (J^2 t^2 / 32)")
    for t in (0.001, 0.01, 0.05):
        print(f"    t = {t:5.3f}: {abs(psi_exact(np.array([xe]), t)[0])**2/(JUMP**2*t**2/32):.4f}")
    lim = JUMP**2/(20*np.sqrt(np.pi))
    print(f"  probability in the gap, P_gap(t) = int_|x|<d-a rho dx; limit of P/t^2.5 is "
          f"J^2/(20 sqrt(pi)) = pi^3.5/60 = {lim:.5f}")
    rows = []
    for t in (0.00125, 0.005, 0.02, 0.1, 0.3, 0.6, 1.0, 1.5, 2.0, 4.0):
        P = gap_probability(t)
        rows.append((t, P))
        print(f"    t = {t:7.5f}: P_gap = {P:.4e}   P/t^2.5 = {P/t**2.5:.5f}")
    tt = np.linspace(0.05, 3.0, 296)
    Pg = np.array([gap_probability(t, 12001) for t in tt])
    t1 = tt[np.argmax(Pg > 0.01)]
    t10 = tt[np.argmax(Pg > 0.10)]
    tm = tt[np.argmax(Pg)]
    print(f"  P_gap passes 1% at t = {t1:.2f}, 10% at t = {t10:.2f}; maximum {Pg.max():.4f} at t = {tm:.2f}")
    print("  far tail at t = 0.01 (exact, and the edge-wave law with its correction):")
    for xf in (10.0, 100.0, 1000.0):
        mp.mp.dps = 60
        e = float(abs(_mp_psi(mp.mpf(xf), mp.mpf("0.01")))**2)
        law = abs(psi_edge(np.array([xf]), 0.01, 1)[0])**2
        print(f"    x = {xf:6.0f}: rho = {e:.4e} (60 digits)   law {law:.4e}")
    results["C"] = dict(exponent=sl.tolist(), P_gap=rows, P_limit=lim,
                        t_1pc=t1, t_10pc=t10, P_max=Pg.max(), t_max=tm)
    return tt, Pg


# ------------------------------------------------------------------- Part D
def winding_zeros(tlo, thi, nt, xlo, xhi, nx, theta=THETA):
    """Zeros of psi in a rectangle of the (x, t) plane, by phase winding."""
    ts = np.geomspace(tlo, thi, nt)
    xs = np.linspace(xlo, xhi, nx)
    ph = np.empty((nt, nx))
    for i, t in enumerate(ts):
        ph[i] = np.angle(psi_exact(xs, t, theta))

    def dw(u, v):
        return (v - u + np.pi) % (2*np.pi) - np.pi

    w = (dw(ph[:-1, :-1], ph[:-1, 1:]) + dw(ph[:-1, 1:], ph[1:, 1:])
         + dw(ph[1:, 1:], ph[1:, :-1]) + dw(ph[1:, :-1], ph[:-1, :-1]))
    q = np.rint(w/(2*np.pi)).astype(int)
    return [(ts[i], xs[j], int(q[i, j])) for i, j in np.argwhere(q != 0)]


def part_D():
    print("\n=== Part D: fringes and nodes (Proposition ZG4) ===")
    print("  spacing of the density minima nearest the centre of the gap:")
    for t in (0.005, 0.01, 0.02, 0.05):
        X = np.linspace(-0.25, 0.25, 200001)
        rho = np.abs(psi_exact(X, t))**2
        i = np.where((rho[1:-1] < rho[:-2]) & (rho[1:-1] < rho[2:]))[0] + 1
        c = i[np.abs(X[i]) < 4.2*np.pi*t/GAP]
        print(f"    t = {t:5.3f}: {np.diff(X[c]).mean():.5f}   pi t/(d-a) = {np.pi*t/GAP:.5f}"
              f"   min/max there = {rho[c].min()/rho[np.abs(X) < 4.2*np.pi*t/GAP].max():.1e}")
    print("  late times, in the variable x/t (far field: minima at (theta + (2n+1) pi)/(2d)):")
    for t in (5.0, 20.0, 100.0):
        u = np.linspace(-1.5, 2.0, 350001)
        rho = np.abs(psi_exact(u*t, t))**2
        i = np.where((rho[1:-1] < rho[:-2]) & (rho[1:-1] < rho[2:]))[0] + 1
        um = u[i]
        print(f"    t = {t:5.0f}: minima at x/t = {np.round(um[:4], 4)}  spacing {np.diff(um).mean():.4f}"
              f"  (pi/d = {np.pi/D:.4f});  min/max = {rho[i].min()/rho.max():.1e};"
              f"  max |t rho - |phi|^2| = {np.abs(t*rho - phi_sq(u)).max():.1e}")
    print(f"    far-field minima expected at x/t = "
          f"{np.round((THETA + np.pi*np.array([-1, 1, 3]))/(2*D), 4)}")
    print("  zeros of psi in the (x, t) plane, |x| < 0.06 (phase winding):")
    allz = []
    for tlo, thi in ((0.005, 0.01), (0.01, 0.02), (0.02, 0.04), (0.04, 0.08)):
        z = winding_zeros(tlo, thi, 6000, -0.06, 0.06, 481)
        allz += z
        tl = max((t for t, *_ in z), default=tlo)
        print(f"    {tlo:.3f} < t < {thi:.3f}: {len(z):4d} zeros, total charge {sum(c for *_, c in z):+d},"
              f" last at t = {tl:.4f};  2 a d (1/t_lo - 1/t_last)/pi = "
              f"{2*A_HALF*D*(1/tlo - 1/tl)/np.pi:.1f}")
    zl = winding_zeros(0.08, 6.0, 4000, -GAP, GAP, 961)
    print(f"    0.08 < t < 6, |x| < d - a: {len(zl)} zeros")
    t0, x0, _ = max(allz)
    mp.mp.dps = 40
    sol = mp.findroot(lambda xx, tt: (mp.re(_mp_psi(xx, tt)), mp.im(_mp_psi(xx, tt))),
                      (mp.mpf(x0), mp.mpf(t0)))
    res = abs(_mp_psi(sol[0], sol[1]))
    print(f"  the last zero, refined at 40 digits: x = {mp.nstr(sol[0], 10)}, t = {mp.nstr(sol[1], 10)},"
          f" |psi| = {mp.nstr(res, 3)}")
    for th, name in ((0.0, "0"), (np.pi, "pi")):
        v = [abs(psi_exact(np.array([0.0]), t, th)[0])**2 for t in (0.1, 0.5)]
        print(f"  theta = {name}: rho(0, t) at t = 0.1, 0.5 = {v[0]:.3e}, {v[1]:.3e}")
    results["D"] = dict(n_zeros=len(allz), last_zero=(float(sol[0]), float(sol[1])))
    return allz


# ------------------------------------------------------------------- Part E
PW = np.linspace(-300, 300, 120001)


def part_E():
    print("\n=== Part E: the Wigner picture (Proposition ZG5) ===")
    dpw = PW[1] - PW[0]
    print("  density from the sheared Wigner function, int W(x - p t, p, 0) dp, against |psi|^2;")
    print("  its two parts: the lobes' tails and the sheared band (= Re e^{i theta} psi_L* psi_R):")
    print("      t      x      |psi|^2        sheared W      lobes          band           band (wave fn)   int|W|dp")
    rows = []
    for t in (0.1, 0.3, 0.6, 1.0):
        for x in (0.0, 0.4):
            xq = x - PW*t
            lob = 0.5*(wigner0(xq - D, PW) + wigner0(xq + D, PW)).sum()*dpw
            bnd = (wigner0(xq, PW)*np.cos(2*PW*D - THETA)).sum()*dpw
            nw = np.abs(wigner(xq, PW)).sum()*dpw
            pl, pr = pulse(np.array([x + D]), t)[0], pulse(np.array([x - D]), t)[0]
            ex = abs(pl + np.exp(1j*THETA)*pr)**2/2
            bw = (np.exp(1j*THETA)*np.conj(pl)*pr).real
            rows.append((t, x, ex, lob + bnd, lob, bnd, nw))
            print(f"    {t:4.1f}   {x:4.1f}   {ex:.6e}   {lob + bnd:.6e}   {lob:.6e}   {bnd:+.6e}"
                  f"   {bw:+.6e}    {nw:.4f}")
    print("  signed weight per unit length at x = 0, int |W(x - p t, p)| dp, against the density:")
    for t in (0.0, 0.05, 0.1, 0.3, 0.6, 1.0, 2.0):
        nw = np.abs(wigner(0.0 - PW*t, PW)).sum()*dpw
        rho = abs(psi_exact(np.array([0.0]), t)[0])**2
        print(f"    t = {t:4.2f}: int |W| dp = {nw:.4f}   rho = {rho:.4e}   ratio = "
              + (f"{nw/rho:.1e}" if rho > 0 else "inf"))
    results["E"] = dict(rows=rows)


# ------------------------------------------------------------------- Part F
def sample_worlds(n, rng, batch=4_000_000):
    """n signed worlds from |W|/Z by rejection, |x| < d + a, |p| < P_MAX."""
    X, P, S = [], [], []
    got = tried = 0
    while got < n:
        x = rng.uniform(-D - A_HALF, D + A_HALF, batch)
        p = rng.uniform(-P_MAX, P_MAX, batch)
        w = wigner(x, p)
        keep = rng.uniform(0, 1/np.pi, batch) < np.abs(w)
        X.append(x[keep])
        P.append(p[keep])
        S.append(np.sign(w[keep]))
        got += int(keep.sum())
        tried += batch
    return (np.concatenate(X)[:n], np.concatenate(P)[:n],
            np.concatenate(S)[:n].astype(int), got/tried)


def signed_density(xw, sw, Z, edges):
    cnt = np.histogram(xw, edges, weights=sw)[0]
    return Z*cnt/(len(xw)*(edges[1] - edges[0]))


def part_F(Z, n=2_000_000):
    print("\n=== Part F: signed worlds (Proposition ZG6) ===")
    rng = np.random.default_rng(20261110)
    x0, p0, s0, acc = sample_worlds(n, rng)
    area = 2*(D + A_HALF)*2*P_MAX
    print(f"  {n} worlds; acceptance {acc:.4f} (Z pi / area = {Z*np.pi/area:.4f});"
          f"  positons {np.mean(s0 > 0):.4f}, expected (1 + 1/Z)/2 = {(1 + 1/Z)/2:.4f}")
    dpw = PW[1] - PW[0]
    width = 0.1
    edges = np.arange(-6.0, 6.0 + 1e-9, width)
    mid = 0.5*(edges[1:] + edges[:-1])
    ig = np.abs(mid) < GAP - 1e-9
    print("  worlds in the gap |x| < d - a, and the signed estimate of P_gap:")
    print("      t      positons   negatons    signed/N * Z    exact P_gap    predicted sigma")
    snaps = {}
    for t in (0.0, 0.1, 0.3, 0.6, 1.0, 1.5):
        xt = x0 + p0*t
        inside = np.abs(xt) < GAP
        npos = int(np.sum(inside & (s0 > 0)))
        nneg = int(np.sum(inside & (s0 < 0)))
        est = Z*(npos - nneg)/n
        sig = Z*np.sqrt(npos + nneg)/n
        ex = gap_probability(t, 12001) if t > 0 else 0.0
        print(f"    {t:4.1f}   {npos:8d}   {nneg:8d}     {est:+.5f}        {ex:.5f}        {sig:.5f}")
        snaps[t] = signed_density(xt, s0, Z, edges)
    print(f"  density in bins of width {width} at x = 0.05: estimate, exact, predicted sigma, pull")
    k = np.argmin(np.abs(mid - 0.05))
    pulls = []
    for t in (0.0, 0.1, 0.3, 0.6, 1.0, 1.5):
        xq = np.linspace(edges[k], edges[k + 1], 21)
        nw = np.mean([np.abs(wigner(xx - PW*t, PW)).sum()*dpw for xx in xq])
        ex = np.mean(np.abs(psi_exact(xq, t))**2)
        sig = np.sqrt(Z*nw/(n*width))
        pulls.append((snaps[t][k] - ex)/sig)
        print(f"    t = {t:3.1f}: {snaps[t][k]:+.5f}   {ex:.5f}   {sig:.5f}   {pulls[-1]:+.2f}")
    print("  normalised residuals over the gap bins, (estimate - exact)/sigma:")
    for t in (0.3, 1.0):
        pl = []
        for kk in np.where(ig)[0]:
            xq = np.linspace(edges[kk], edges[kk + 1], 9)
            nw = np.mean([np.abs(wigner(xx - PW*t, PW)).sum()*dpw for xx in xq])
            ex = np.mean(np.abs(psi_exact(xq, t))**2)
            pl.append((snaps[t][kk] - ex)/np.sqrt(Z*nw/(n*width)))
        print(f"    t = {t:3.1f}: mean {np.mean(pl):+.2f}, rms {np.sqrt(np.mean(np.square(pl))):.2f}"
              f"  ({len(pl)} bins)")
    print("  worlds needed to see rho(0, t) at 3 sigma in a bin of width 0.1,"
          " N = 9 Z int|W|dp / (0.1 rho^2):")
    need = []
    for t in (0.05, 0.1, 0.2, 0.3, 0.6, 1.0):
        nw = np.abs(wigner(0.0 - PW*t, PW)).sum()*dpw
        rho = abs(psi_exact(np.array([0.0]), t)[0])**2
        need.append((t, 9*Z*nw/(width*rho**2)))
        print(f"    t = {t:4.2f}: rho = {rho:.3e}   N = {need[-1][1]:.1e}")
    results["F"] = dict(n=n, acceptance=acc, pulls=pulls, need=need)
    return x0, p0, s0, edges, mid, snaps


# ------------------------------------------------------------------- Part G
def part_G(x0, p0, s0, Z):
    print("\n=== Part G: a harmonic trap, omega = 1 (Remark ZG7) ===")
    pw = np.linspace(-60, 60, 48001)
    dpw = pw[1] - pw[0]

    def rho_rot(x, t):
        c, s = np.cos(t), np.sin(t)
        return wigner(x*c - pw*s, x*s + pw*c).sum()*dpw

    print("  density from the rotated Wigner function, int W(x cos t - p sin t, x sin t + p cos t) dp:")
    for x in (0.0, 0.5, 1.1):
        print(f"    x = {x:3.1f}:  t = pi/2: {rho_rot(x, np.pi/2):.6f}  (|phi(x)|^2 = {phi_sq(x):.6f});"
              f"   t = pi: {rho_rot(x, np.pi):+.1e};  t = 2 pi: {rho_rot(x, 2*np.pi):+.1e}")
    xn = (THETA + np.pi*np.array([-1, 1, 3]))/(2*D)
    print(f"  exact nodes at t = pi/2: x = (theta + (2n+1) pi)/(2d) = {np.round(xn, 4)};"
          f"  |phi|^2 there = {np.round(phi_sq(xn), 18)}")
    n = len(x0)
    print("  signed worlds rotated rigidly; signed count in the gap times Z/N:")
    for name, t in (("0", 0.0), ("pi/4", np.pi/4), ("pi/2", np.pi/2), ("pi", np.pi)):
        xt = x0*np.cos(t) + p0*np.sin(t)
        inside = np.abs(xt) < GAP
        est = Z*np.sum(s0[inside])/n
        sig = Z*np.sqrt(inside.sum())/n
        u = np.linspace(-GAP, GAP, 4801)
        if t == np.pi/2:
            ex = np.trapezoid(phi_sq(u), u)
        elif name in ("0", "pi"):
            ex = 0.0
        else:
            ex = np.trapezoid([rho_rot(v, t) for v in u[::40]], u[::40])
        print(f"    t = {name:5s}: {est:+.5f} +- {sig:.5f}   exact {ex:.5f}")
    results["G"] = dict(nodes=xn.tolist())


# ------------------------------------------------------------------ figures
def figures(tt, Pg, allz, x0, p0, s0, Z, edges, mid, snaps):
    print("\n=== Figures ===")
    # 1. the density: linear and logarithmic
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    X = np.linspace(-6, 6, 4801)
    for t, c in ((0.0, "k"), (0.1, "C0"), (0.3, "C2"), (0.6, "C1"), (1.0, "C3")):
        r = np.abs(psi_exact(X, t))**2
        ax[0].plot(X, r, c, lw=1.2, label=f"t = {t}")
        ax[1].semilogy(X, np.maximum(r, 1e-30), c, lw=1.0)
    xl = np.linspace(-0.9, 0.9, 721)
    for t, c in ((0.01, "C4"), (0.03, "C5")):
        ax[1].semilogy(xl, np.abs(psi_exact(xl, t))**2, c, lw=0.8, label=f"t = {t}")
        ax[1].semilogy(xl, np.abs(psi_edge(xl, t, 1))**2, "k:", lw=0.8)
    for a_ in ax[:2]:
        a_.axvspan(-GAP, GAP, color="0.92", zorder=0)
        a_.set_xlabel("x / a")
    ax[0].set_ylabel("density")
    ax[0].set_xlim(-6, 6)
    ax[0].legend(fontsize=8)
    ax[0].set_title("Exact density (gap shaded)")
    ax[1].set_xlim(-4, 4)
    ax[1].set_ylim(1e-12, 2)
    ax[1].legend(fontsize=7, loc="lower left")
    ax[1].set_title("Logarithmic; dotted: edge-wave law")
    t2 = np.geomspace(1e-3, 3, 60)
    P2 = np.array([gap_probability(t, 24001) for t in t2])
    ax[2].loglog(t2, P2, "C0", label="probability in the gap")
    ax[2].loglog(t2, JUMP**2*t2**2.5/(20*np.sqrt(np.pi)), "k--", lw=0.8,
                 label=r"$J^2 t^{5/2}/(20\sqrt{\pi})$")
    r0 = np.array([abs(psi_exact(np.array([0.0]), t)[0])**2 for t in t2])
    ax[2].loglog(t2, r0, "C3", label="density at x = 0")
    ax[2].loglog(t2, JUMP**2*t2**5/(2*np.pi*GAP**6), "k:", lw=0.8,
                 label=r"$J^2 t^5/(2\pi (d-a)^6)$")
    ax[2].set_ylim(1e-16, 2)
    ax[2].set_xlabel("t")
    ax[2].legend(fontsize=8)
    ax[2].set_title("Two power laws")
    save(fig, "zero_gap_density.png")

    # 2. the Wigner function at four times
    fig, ax = plt.subplots(1, 4, figsize=(16, 3.8), sharey=True)
    xg = np.linspace(-5, 5, 501)
    pg = np.linspace(-7, 7, 421)
    for a_, t in zip(ax, (0.0, 0.1, 0.3, 1.0)):
        Wt = np.pi*wigner(xg[None, :] - pg[:, None]*t, pg[:, None])
        im = a_.pcolormesh(xg, pg, Wt, cmap="RdBu_r", vmin=-1, vmax=1, shading="auto", rasterized=True)
        for v in (-GAP, GAP):
            a_.axvline(v, color="k", lw=0.5, ls=":")
        a_.set_title(f"t = {t}")
        a_.set_xlabel("x / a")
    ax[0].set_ylabel("p a / hbar")
    fig.colorbar(im, ax=ax, label=r"$\pi\hbar W$", fraction=0.015, pad=0.01)
    save(fig, "zero_gap_wigner.png")

    # 3. the gap in the (x, t) plane, with the zeros
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.6))
    ts = np.geomspace(0.004, 3.0, 700)
    xs = np.linspace(-GAP, GAP, 1201)
    R = np.array([np.abs(psi_exact(xs, t))**2 for t in ts])
    im = ax[0].pcolormesh(xs, ts, np.log10(np.maximum(R, 1e-30)), cmap="viridis",
                          vmin=-14, vmax=0, shading="auto", rasterized=True)
    ax[0].set_yscale("log")
    ax[0].set_xlabel("x / a")
    ax[0].set_ylabel("t")
    ax[0].set_title(r"$\log_{10}$ density in the gap")
    fig.colorbar(im, ax=ax[0])
    ts2 = np.geomspace(0.01, 0.04, 900)
    xs2 = np.linspace(-0.06, 0.06, 601)
    R2 = np.array([np.abs(psi_exact(xs2, t))**2/t**5 for t in ts2])
    im = ax[1].pcolormesh(xs2, ts2, np.log10(R2), cmap="viridis", shading="auto", rasterized=True)
    zz = [(t, x, c) for t, x, c in allz if 0.01 <= t <= 0.04]
    ax[1].plot([x for t, x, c in zz if c > 0], [t for t, x, c in zz if c > 0], "r.", ms=3, label="zero, charge +1")
    ax[1].plot([x for t, x, c in zz if c < 0], [t for t, x, c in zz if c < 0], "w.", ms=3, label="zero, charge -1")
    ax[1].set_yscale("log")
    ax[1].set_xlabel("x / a")
    ax[1].set_title(r"$\log_{10}(\rho/t^5)$ near the centre, with the zeros of $\psi$")
    ax[1].legend(fontsize=8, loc="upper right")
    fig.colorbar(im, ax=ax[1])
    save(fig, "zero_gap_spacetime.png")

    # 4. signed worlds
    fig, ax = plt.subplots(2, 3, figsize=(15, 7))
    Xf = np.linspace(-6, 6, 2401)
    for a_, t in zip(ax[0], (0.0, 0.6, 1.5)):
        a_.step(mid, snaps[t], where="mid", color="C0", lw=0.9, label="signed worlds")
        a_.plot(Xf, np.abs(psi_exact(Xf, t))**2, "k", lw=0.9, label="exact")
        a_.axvspan(-GAP, GAP, color="0.92", zorder=0)
        a_.set_title(f"t = {t}: signed histogram, N = {len(x0):.0e}")
        a_.legend(fontsize=8, loc="upper right")
    for a_, t in zip(ax[1], (0.0, 0.6, 1.5)):
        xt = x0 + p0*t
        hp = np.histogram(xt[s0 > 0], edges)[0]*Z/(len(x0)*(edges[1] - edges[0]))
        hn = np.histogram(xt[s0 < 0], edges)[0]*Z/(len(x0)*(edges[1] - edges[0]))
        a_.step(mid, hp, where="mid", color="C3", lw=0.9, label="positons")
        a_.step(mid, hn, where="mid", color="C0", lw=0.9, label="negatons")
        a_.axvspan(-GAP, GAP, color="0.92", zorder=0)
        a_.set_title(f"t = {t}: the two populations")
        a_.set_xlabel("x / a")
        a_.legend(fontsize=8)
    fig.subplots_adjust(hspace=0.28)
    save(fig, "zero_gap_worlds.png")


def make_animation(x0, p0, s0, Z):
    print("\n=== Animation ===")
    times = np.concatenate([np.zeros(6), np.linspace(0, 0.3, 31)[1:], np.linspace(0.3, 1.5, 49)[1:]])
    xg = np.linspace(-5.5, 5.5, 441)
    pg = np.linspace(-7, 7, 281)
    Xf = np.linspace(-5.5, 5.5, 1101)
    edges = np.arange(-5.5, 5.5 + 1e-9, 0.1)
    mid = 0.5*(edges[1:] + edges[:-1])
    fig, ax = plt.subplots(3, 1, figsize=(6.4, 8.2), gridspec_kw=dict(height_ratios=[2.1, 1, 1]))
    im = ax[0].pcolormesh(xg, pg, np.pi*wigner(xg[None, :], pg[:, None]), cmap="RdBu_r",
                          vmin=-1, vmax=1, shading="auto", rasterized=True)
    ax[0].set_ylabel("p a / hbar")
    ttl = ax[0].set_title(r"$\pi\hbar W(x,p,t)$,  t = 0.00   (gap between dotted lines)", fontsize=10)
    for a_ in ax:
        for v in (-GAP, GAP):
            a_.axvline(v, color="k", lw=0.5, ls=":")
        a_.set_xlim(-5.5, 5.5)
    ln_h, = ax[1].step(mid, np.zeros_like(mid), where="mid", color="C0", lw=0.9, label="signed worlds")
    ln_e, = ax[1].plot(Xf, np.zeros_like(Xf), "k", lw=1.0, label="exact")
    ax[1].set_ylim(-0.08, 0.7)
    ax[1].set_ylabel("density")
    ax[1].legend(fontsize=7, loc="upper right")
    ln_l, = ax[2].semilogy(Xf, np.ones_like(Xf), "k", lw=1.0)
    ax[2].set_ylim(1e-10, 2)
    ax[2].set_ylabel("density (log)")
    ax[2].set_xlabel("x / a")
    fig.tight_layout()

    def frame(i):
        t = times[i]
        im.set_array((np.pi*wigner(xg[None, :] - pg[:, None]*t, pg[:, None])).ravel())
        r = np.abs(psi_exact(Xf, t))**2
        ln_e.set_ydata(r)
        ln_l.set_ydata(np.maximum(r, 1e-30))
        ln_h.set_ydata(signed_density(x0 + p0*t, s0, Z, edges))
        ttl.set_text(rf"$\pi\hbar W(x,p,t)$,  t = {t:.2f}   (gap between dotted lines)")
        return im, ln_e, ln_l, ln_h, ttl

    ani = animation.FuncAnimation(fig, frame, frames=len(times), blit=False)
    name = "zero_gap_evolution.gif"
    ani.save(output_path(name), writer=animation.PillowWriter(fps=8), dpi=72)
    dp = docs_path(name)
    if dp:
        ani.save(dp, writer=animation.PillowWriter(fps=8), dpi=72)
    print("  wrote", name)
    try:
        ani.save(output_path("zero_gap_evolution.mp4"), writer=animation.FFMpegWriter(fps=8, bitrate=1800), dpi=110)
        print("  wrote zero_gap_evolution.mp4")
    except Exception as exc:                      # ffmpeg is optional
        print("  mp4 not written:", exc)
    plt.close(fig)


if __name__ == "__main__":
    Z = part_A()
    part_B()
    tt, Pg = part_C()
    allz = part_D()
    part_E()
    x0, p0, s0, edges, mid, snaps = part_F(Z)
    part_G(x0, p0, s0, Z)
    figures(tt, Pg, allz, x0, p0, s0, Z, edges, mid, snaps)
    if "--no-animation" not in sys.argv:
        make_animation(x0, p0, s0, Z)
    with open(output_path("zero_gap_superposition_results.json"), "w") as f:
        json.dump(results, f, indent=1, default=str)
    print("\nwrote", output_path("zero_gap_superposition_results.json"))
