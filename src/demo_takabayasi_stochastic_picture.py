#!/usr/bin/env python3
r"""
Verification for the comparison of Takabayasi (1954) Prog. Theor. Phys. 11, 341,
section 3(a) -- "Equation of motion and the transition in momentum" -- against the
WPMW four-action / phase-space-crystal-lattice model.

Units: hbar = m = 1.

Parts
  A  Takabayasi (3.7) quadrature  ==  closed form (3.8), with the 1D prefactor.
  B  Takabayasi's collision operator A[f] = \int J(x,p-p') f(p') dp'
       ==  exact Wigner collision term, and == Moyal series.
  C  Periodic potential: J is supported exactly on the momentum lattice
       p = q * (pi hbar / L), and A[f] collapses to the four-action stencil
       Gamma_q(x) [ W(p + q d) - W(p - q d) ].
  D  Transition moments (3.17): even ones vanish, odd ones are the Moyal
       coefficients; Kramers-Moyal truncation at order 2 == classical Liouville.
  E  Free grid step: moments with jump +/- delta reproduce (3.17) with
       hbar -> hbar_eff = 2 delta / k   (project Proposition F4, in Takabayasi's
       own moment language).
  F  Rate ledger: \int J dp = 0, but the sea-dressed channel rates are all
       non-negative and sum to something large and positive.
  G  No-go: the finite-time kernel T of (3.28)-(3.33) is normalized AND
       orthogonal; a non-negative such kernel must be a permutation.
       Hence negativity is forced for any non-classical evolution.
"""

import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from wpmwlib.wpmw_utils import output_path, docs_path  # noqa: E402

np.set_printoptions(precision=4, suppress=True)
HBAR = 1.0
FIGNAME = "takabayasi_stochastic_picture.png"


def banner(s):
    print("\n" + "=" * 72)
    print(s)
    print("=" * 72)


# ----------------------------------------------------------------------
# Part A -- Takabayasi (3.7) vs (3.8)
# ----------------------------------------------------------------------
banner("PART A  --  Takabayasi (3.7) quadrature  vs  closed form (3.8)")

V0, s = 1.3, 2.5


def V_gauss(x):
    return V0 * np.exp(-x ** 2 / (2 * s ** 2))


def Vtil_gauss(p):
    """Fourier component, Takabayasi (3.9) in 1D: (2 pi hbar)^-1 int V e^{ipx/hbar} dx."""
    return V0 * s * np.exp(-s ** 2 * p ** 2 / (2 * HBAR ** 2)) / (np.sqrt(2 * np.pi) * HBAR)


# quadrature of (3.7) in 1D
Y = np.linspace(-400.0, 400.0, 4_000_001)
dY = Y[1] - Y[0]


def J_quad(x, p):
    return -(2.0 / HBAR) * (1.0 / (2 * np.pi * HBAR)) * np.trapezoid(
        V_gauss(x + Y / 2) * np.sin(p * Y / HBAR), dx=dY)


def J_closed(x, p):
    """1D form of (3.8): the 2^4 there is the 3D prefactor 2^(d+1); in 1D it is 2^2."""
    return -(2 ** 2 / HBAR) * np.imag(Vtil_gauss(2 * p) * np.exp(-2j * p * x / HBAR))


err = 0.0
for x in (-1.1, 0.0, 0.37, 2.0):
    for p in (-1.7, -0.4, 0.9, 2.3):
        err = max(err, abs(J_quad(x, p) - J_closed(x, p)))
print(f"  max |J_(3.7)  -  J_(3.8)|                      = {err:.3e}")
print(f"  prefactor check: 2^(d+1) with d=3 gives          {2 ** 4}   (paper's 2^4)")

# feature (iv): x-oscillation of J is a pure sinusoid of wavelength h/(2p)
p_test = 0.9
lam_expected = 2 * np.pi * HBAR / (2 * p_test)          # = h/(2p)
xs = np.linspace(0, 2 * lam_expected, 9)
Jx = np.array([J_quad(x, p_test) for x in xs])          # from the quadrature, not the fit
amp = (2 ** 2 / HBAR) * Vtil_gauss(2 * p_test)
fit = amp * np.sin(2 * p_test * xs / HBAR)
print(f"  feature (iv): |J_quad(x,p) - A sin(2px/hbar)|   = {np.abs(Jx - fit).max():.3e}")
print(f"               undamped amplitude A = 4 Vtilde(2p)/hbar = {amp:.6f}")
print(f"               wavelength  h/(2p)                 = {lam_expected:.6f}")
print(f"               paper text says 'hbar/2p'          = {HBAR / (2 * p_test):.6f}"
      f"   <-- off by 2*pi; h/2p is the correct value")


# ----------------------------------------------------------------------
# Part B -- A[f] vs exact Wigner collision term vs Moyal series
# ----------------------------------------------------------------------
banner("PART B  --  A[f] = int J(x,p-p') f(p') dp'   vs exact QLE  vs Moyal")

# Gaussian test state at fixed x
sp = 1.5


def f_test(p):
    return np.exp(-p ** 2 / (2 * sp ** 2)) / (np.sqrt(2 * np.pi) * sp)


P = np.linspace(-30, 30, 240001)
dP = P[1] - P[0]
x0 = 0.37

# Takabayasi convolution, on the p-grid
Jk = np.array([J_closed(x0, pp) for pp in P])          # J(x0, .) sampled
fk = f_test(P)
conv = np.convolve(Jk, fk, mode="same") * dP           # (J * f)(p)

# exact collision term via the (3.3)/(3.4) route:
#   dbar_rho/dt|_coll = (i/hbar) [V(x+y/2) - V(x-y/2)] bar_rho ,  bar_rho = FT_p->y of f
Ny = len(P)
y = 2 * np.pi * HBAR * np.fft.fftfreq(Ny, d=dP)
rho_bar = np.fft.fft(fk) * dP                          # int f e^{-ipy/hbar} dp  (up to convention)
kernel = (1j / HBAR) * (V_gauss(x0 + y / 2) - V_gauss(x0 - y / 2))
exact = np.real(np.fft.ifft(rho_bar * kernel) / dP)

sel = (np.abs(P) < 4.0)
scale = np.abs(exact[sel]).max()
print(f"  max |A[f]_Takabayasi - exact QLE| / scale      = "
      f"{np.abs(conv[sel] - exact[sel]).max() / scale:.3e}")

# Moyal series  (2/hbar) sum_{n odd} (-1)^{(n-1)/2}/n! (hbar/2)^n  d^n V  d^n f
from numpy.polynomial import hermite_e


def dnf(n):
    """n-th p-derivative of the Gaussian f_test on the grid, analytic via Hermite."""
    c = np.zeros(n + 1)
    c[n] = 1.0
    return f_test(P) * ((-1) ** n) * hermite_e.hermeval(P / sp, c) / sp ** n


def dnV(n, x):
    """n-th x-derivative of the Gaussian V at x, analytic via Hermite."""
    c = np.zeros(n + 1)
    c[n] = 1.0
    return V_gauss(x) * ((-1) ** n) * hermite_e.hermeval(x / s, c) / s ** n


moyal = np.zeros_like(P)
for n in range(1, 26, 2):
    coef = (2.0 / HBAR) * ((-1) ** ((n - 1) // 2)) / math.factorial(n) * (HBAR / 2) ** n
    moyal = moyal + coef * dnV(n, x0) * dnf(n)
print(f"  max |Moyal series (25 terms) - exact| / scale  = "
      f"{np.abs(moyal[sel] - exact[sel]).max() / scale:.3e}")


# ----------------------------------------------------------------------
# Part C -- periodic V: J lives on the momentum lattice; four-action stencil
# ----------------------------------------------------------------------
banner("PART C  --  periodic potential: J's support IS the momentum lattice")

L = 8.0
modes = {1: (1.5, 0.0), 2: (-0.7, 0.9), 5: (0.4, -2.1)}      # q -> (V_q, phi_q)
delta = np.pi * HBAR / L                                      # project's dp = pi hbar / L


def V_per(x):
    out = np.zeros_like(np.asarray(x, dtype=float))
    for q, (Vq, ph) in modes.items():
        out = out + Vq * np.cos(2 * np.pi * q * np.asarray(x) / L + ph)
    return out


print(f"  momentum grid step  delta = pi hbar / L        = {delta:.8f}")
for q in modes:
    kq = 2 * np.pi * q / L
    print(f"    mode q={q}:  hbar k_q / 2 = {HBAR * kq / 2:.8f}   "
          f"= {HBAR * kq / 2 / delta:.1f} * delta")

# exact QLE collision term for a periodic V, evaluated on a p-grid of step delta
Np = 512
Pg = (np.arange(Np) - Np // 2) * delta
Wg = np.exp(-Pg ** 2 / (2 * 1.1 ** 2)) * (1 + 0.3 * np.cos(3.0 * Pg))   # arbitrary state

xs = np.linspace(0, L, 17, endpoint=False)
worst = 0.0
for x in xs:
    # exact: shift form,  C = sum_q Gamma_q(x) [W(p + q delta) - W(p - q delta)]
    stencil = np.zeros_like(Wg)
    for q, (Vq, ph) in modes.items():
        Gam = -(Vq / HBAR) * np.sin(2 * np.pi * q * x / L + ph)
        stencil += Gam * (np.roll(Wg, -q) - np.roll(Wg, +q))
    # independent route: the Moyal series for the same periodic V
    moy = np.zeros_like(Wg)
    for n in range(1, 40, 2):
        # n-th derivative of V_per at x
        dV = 0.0
        for q, (Vq, ph) in modes.items():
            kq = 2 * np.pi * q / L
            dV += Vq * kq ** n * np.cos(kq * x + ph + n * np.pi / 2)
        # n-th derivative of W on the lattice, spectrally
        kk = 2 * np.pi * np.fft.fftfreq(Np, d=delta)
        dW = np.real(np.fft.ifft((1j * kk) ** n * np.fft.fft(Wg)))
        moy += (2.0 / HBAR) * ((-1) ** ((n - 1) // 2)) / math.factorial(n) \
               * (HBAR / 2) ** n * dV * dW
    m = np.abs(stencil).max()
    worst = max(worst, np.abs(stencil - moy).max() / m)
print(f"  max rel |four-action stencil - Moyal series|   = {worst:.3e}")
print("  => for a periodic V the Wigner kernel is a pure finite-difference")
print("     operator on the lattice p = q * pi hbar / L : mode q couples +/- q cells.")


# ----------------------------------------------------------------------
# Part D -- transition moments (3.17)
# ----------------------------------------------------------------------
banner("PART D  --  transition moments (3.17), and Kramers-Moyal truncation")

Vq, kq = 1.5, 2 * np.pi * 1 / L
x1 = 0.83
Gam = -(Vq / HBAR) * np.sin(kq * x1)
jump = HBAR * kq / 2

print("   n |  m_n from J = Gam[d(p+hk/2) - d(p-hk/2)] |  -(-hbar^2/4)^((n-1)/2) d^n V | diff")
eps = HBAR ** 2 / 4
for n in range(0, 8):
    m_delta = Gam * ((-jump) ** n - (jump) ** n)
    if n % 2 == 0:
        m_closed = 0.0
    else:
        dnV_cos = Vq * kq ** n * np.cos(kq * x1 + n * np.pi / 2)
        m_closed = -((-eps) ** ((n - 1) // 2)) * dnV_cos
    print(f"  {n:2d} |  {m_delta: .12e}          |  {m_closed: .12e}  | {abs(m_delta-m_closed):.1e}")

print()
print(f"  m_1 = -dV/dx :  {Gam*(-2*jump): .8f}  vs  {Vq*kq*np.sin(kq*x1): .8f}")
print(f"  m_2 (the Fokker-Planck diffusion coefficient)  = {Gam*((-jump)**2-jump**2):.3e}")
print("  => Kramers-Moyal truncated at order 2 is EXACTLY classical Liouville:")
print("     no diffusion term exists at all.  All quantum content is in n >= 3.")


# ----------------------------------------------------------------------
# Part E -- free grid step  <->  hbar_eff  (Proposition F4 in moment language)
# ----------------------------------------------------------------------
banner("PART E  --  free jump size delta  <->  hbar_eff = 2 delta / k")

print("  delta/(hbar k/2) |  hbar_eff |  max_n<=7 |m_n(delta) - m_n^(3.17)(hbar_eff)|")
for ratio in (0.05, 0.25, 1.0, 2.0, 5.0):
    d_ = ratio * HBAR * kq / 2
    h_eff = 2 * d_ / kq
    Gam_eff = -(Vq / h_eff) * np.sin(kq * x1)      # Ehrenfest fixes Gamma once delta is free
    worst = 0.0
    for n in range(1, 8, 2):
        m_d = Gam_eff * ((-d_) ** n - d_ ** n)
        dnV_cos = Vq * kq ** n * np.cos(kq * x1 + n * np.pi / 2)
        m_c = -((-(h_eff ** 2 / 4)) ** ((n - 1) // 2)) * dnV_cos
        worst = max(worst, abs(m_d - m_c))
    print(f"      {ratio:6.2f}       | {h_eff:8.4f}  |  {worst:.3e}")
print("  => shrinking the jump is exactly sending hbar -> 0 in (3.17).")


# ----------------------------------------------------------------------
# Part F -- rate ledger
# ----------------------------------------------------------------------
banner("PART F  --  int J dp = 0, but the channel rates are non-negative")

# The QLE stencil read as a one-body generator L_{nm}:
#   L[n, n+q] = +Gamma,  L[n, n-q] = -Gamma,  L[n, n] = 0.
Nc, qq = 9, 1
Lgen = np.zeros((Nc, Nc))
for n in range(Nc):
    Lgen[n, (n + qq) % Nc] += Gam
    Lgen[n, (n - qq) % Nc] -= Gam
offdiag = Lgen[~np.eye(Nc, dtype=bool)]
print(f"  column sums of L (particle number)            = {np.abs(Lgen.sum(axis=0)).max():.2e}")
print(f"  diagonal of L (= -total exit rate)            = {np.abs(np.diag(Lgen)).max():.2e}")
print(f"  most negative OFF-diagonal entry of L         = {offdiag.min(): .6f}")
print("  => L is not a Markov generator: an off-diagonal rate is negative and")
print("     the exit rate is zero.  This is exactly Takabayasi's (i) + (3.14).")
print()

B = 25000.0                       # sea pairs per cell (SD level 1)
Wex = np.array([12.0, -5.0, 30.0])   # excess occupancies at (lo, n, hi)
gam = abs(Gam)
net_focus = 0.5 * Gam * (Wex[2] - Wex[0])
net_hop = -0.5 * Gam * (Wex[2] + Wex[0])
# SD channel rates (K1..K4b) are (gamma/2) * (non-negative occupancy), with the
# sea supplying the partner; take U+ = max(W,0)+0, U- = max(-W,0) plus sea traffic
Up = np.maximum(Wex, 0)
Um = np.maximum(-Wex, 0)
chan = 0.5 * gam * np.array([Up[2], Um[2], B, B, Up[2], Um[2], B, B])
print(f"  signed net focus rate f_n                     = {net_focus: .6f}")
print(f"  signed net hop   rate h_n                     = {net_hop: .6f}")
print(f"  sum of the eight non-negative channel rates   = {chan.sum(): .6f}")
print(f"  net / total traffic                           = "
      f"{abs(net_focus)+abs(net_hop):.4f} / {chan.sum():.1f} "
      f"= {(abs(net_focus)+abs(net_hop))/chan.sum():.2e}")
print("  => Takabayasi's zero total rate is the *bias*; the traffic is the sea.")


# ----------------------------------------------------------------------
# Part G -- the finite-time kernel T: normalized + orthogonal => negativity forced
# ----------------------------------------------------------------------
banner("PART G  --  T is normalized AND orthogonal (3.29)/(3.32)/(3.33)")

d = 5
w = np.exp(2j * np.pi / d)
X = np.roll(np.eye(d), 1, axis=0)
Z = np.diag(w ** np.arange(d))
inv2 = (d + 1) // 2
Pop = np.zeros((d, d))
for k in range(d):
    Pop[(-k) % d, k] = 1.0          # parity |k> -> |-k>


def D(q, p):
    return (w ** (-q * p * inv2)) * np.linalg.matrix_power(X, q) @ np.linalg.matrix_power(Z, p)


A = {}
for q in range(d):
    for p in range(d):
        Dq = D(q, p)
        A[(q, p)] = Dq @ Pop @ Dq.conj().T

keys = list(A.keys())
G = np.array([[np.trace(A[a] @ A[b]).real for b in keys] for a in keys])
print(f"  |Tr[A_a A_b] - d delta_ab|                     = {np.abs(G - d*np.eye(d*d)).max():.3e}")
print(f"  |sum_a A_a - d I|                              = "
      f"{np.abs(sum(A.values()) - d*np.eye(d)).max():.3e}")

rng = np.random.default_rng(0)
H = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
H = (H + H.conj().T) / 2
evals, evecs = np.linalg.eigh(H)
U = evecs @ np.diag(np.exp(-1j * evals * 0.6)) @ evecs.conj().T

T = np.array([[np.trace(A[a] @ U @ A[b] @ U.conj().T).real / d for b in keys] for a in keys])
print(f"  column sums of T   (3.33)                      "
      f"min={T.sum(axis=0).min():.12f}  max={T.sum(axis=0).max():.12f}")
print(f"  |T^T T - I|        (3.29)+(3.32)               = "
      f"{np.abs(T.T @ T - np.eye(d*d)).max():.3e}")
print(f"  most negative entry of T                       = {T.min(): .6f}")
print(f"  sum_a T_ab^2  (must be 1)                      = {(T**2).sum(axis=0).max():.12f}")
print()
print("  Theorem.  If T >= 0 entrywise, column sums 1, and sum_a T_ab^2 = 1, then")
print("  for each column x: ||x||_1 = 1 and ||x||_2 = 1, which for x >= 0 forces x")
print("  to be a standard basis vector.  T is then a permutation, i.e. a")
print("  deterministic classical flow.  Negativity is therefore NOT an artifact of")
print("  Takabayasi's construction: it is forced for every non-classical evolution")
print("  representable by a normalized, orthogonal phase-space kernel.")


# ----------------------------------------------------------------------
# Part H -- the Wigner momentum lattice of a period-L system is pi hbar / L
# ----------------------------------------------------------------------
banner("PART H  --  the state lattice equals the jump lattice")

rng2 = np.random.default_rng(7)
nmodes = np.array([-2, -1, 0, 1, 3])
cs = rng2.normal(size=len(nmodes)) + 1j * rng2.normal(size=len(nmodes))
ks = 2 * np.pi * nmodes / L


def psi(xv):
    return np.sum(cs[:, None] * np.exp(1j * ks[:, None] * xv[None, :]), axis=0)


xv = np.array([0.31])
yv = np.linspace(-6, 6, 4001)
lhs = psi(xv - yv / 2) * np.conj(psi(xv + yv / 2))          # bar_rho(x,y),  Takabayasi (3.3)

# predicted:  f(x,p) = sum_{n,m} c_n c_m^* e^{i(k_n-k_m)x} delta(p - hbar(k_n+k_m)/2)
rhs = np.zeros_like(yv, dtype=complex)
plist = []
for i, n in enumerate(nmodes):
    for j, m in enumerate(nmodes):
        A = cs[i] * np.conj(cs[j]) * np.exp(1j * (ks[i] - ks[j]) * xv[0])
        pnm = HBAR * (ks[i] + ks[j]) / 2
        plist.append(pnm)
        rhs += A * np.exp(-1j * pnm * yv / HBAR)
print(f"  |bar_rho(x,y) - sum_nm A_nm e^(-i p_nm y/hbar)| = {np.abs(lhs - rhs).max():.3e}")
pl = np.unique(np.round(np.array(plist) / delta, 9))
print(f"  support of W in units of delta = pi hbar / L    = {pl}")
print(f"  all integers?                                   {np.allclose(pl, np.round(pl))}")
print("  => W of a period-L state lives on the SAME lattice the jumps use.")


# ----------------------------------------------------------------------
# Figure
# ----------------------------------------------------------------------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.6))

# (a) J(x,p) for a localized V -- feature (iv)
s_disp = 0.7                                  # narrower well: wider Vtilde, visible fringes
xg = np.linspace(-6, 6, 700)
pg = np.linspace(-2.2, 2.2, 700)
XX, PP = np.meshgrid(xg, pg)
Vt_disp = V0 * s_disp * np.exp(-s_disp ** 2 * (2 * PP) ** 2 / (2 * HBAR ** 2)) \
    / (np.sqrt(2 * np.pi) * HBAR)
JJ = -(2 ** 2 / HBAR) * Vt_disp * np.imag(np.exp(-2j * PP * XX / HBAR))
v = np.abs(JJ).max()
ax[0].pcolormesh(XX, PP, JJ, cmap="RdBu_r", vmin=-v, vmax=v, shading="auto")
for pp in (0.5, 1.0, 1.6):
    lam = 2 * np.pi * HBAR / (2 * pp)
    ax[0].annotate("", xy=(-5.6 + lam, pp), xytext=(-5.6, pp),
                   arrowprops=dict(arrowstyle="<->", color="k", lw=1.3))
    ax[0].text(-5.6 + lam / 2, pp + 0.09, r"$h/2p$", ha="center", fontsize=9)
ax[0].set_xlabel("$x$"); ax[0].set_ylabel("jump $p$")
ax[0].set_title("(a) Takabayasi's $J(x,p)$, feature (iv):\nundamped $x$-grating of wavelength $h/2p$")

# (b) periodic V -- J is a comb on the momentum lattice
ax[1].axhline(0, color="0.7", lw=1)
for q, (Vq, ph) in modes.items():
    G = -(Vq / HBAR) * np.sin(2 * np.pi * q * 0.83 / L + ph)
    for sgn in (+1, -1):
        ax[1].plot([sgn * q, sgn * q], [0, sgn * G], color="C0", lw=2.5)
        ax[1].plot(sgn * q, sgn * G, "o", color="C0", ms=6)
for c in range(-6, 7):
    ax[1].plot(c, 0, "|", color="0.5", ms=9)
ax[1].set_xticks(range(-6, 7))
ax[1].set_xlabel(r"jump $/\ \delta$,   $\delta=\pi\hbar/L$")
ax[1].set_ylabel(r"weight $\pm\Gamma_q(x)$")
ax[1].set_title("(b) periodic $V$: $J$ is a comb on the\nmomentum lattice; mode $q$ at $\\pm q$ cells")

# (c) the half-photon problem
ax[2].set_xlim(-0.6, 3.4); ax[2].set_ylim(-2.6, 2.6); ax[2].axis("off")
for row in (-1, 0, 1):
    ax[2].plot([-0.4, 3.2], [row, row], color="0.85", lw=1)
    ax[2].text(-0.5, row, f"$n{row:+d}q$".replace("+0q", ""), ha="right", va="center", fontsize=9)
ax[2].text(0.35, 2.15, "Takabayasi", ha="center", fontweight="bold")
ax[2].text(0.35, 1.75, r"one body, jump $\pm\hbar k/2$", ha="center", fontsize=9, color="crimson")
ax[2].annotate("", xy=(0.35, 1), xytext=(0.35, 0), arrowprops=dict(arrowstyle="->", color="crimson", lw=2))
ax[2].text(0.55, 0.5, r"$\frac{1}{2}$ photon", color="crimson", fontsize=9)
ax[2].text(0.35, -2.1, "signed rate\n(not a Markov generator)", ha="center", fontsize=8, color="crimson")
ax[2].text(2.2, 2.15, "Four actions", ha="center", fontweight="bold")
ax[2].text(1.5, 1.75, "hop", ha="center", fontsize=9, color="C0")
ax[2].annotate("", xy=(1.5, 1), xytext=(1.5, -1), arrowprops=dict(arrowstyle="->", color="C0", lw=2))
ax[2].text(1.6, 0.0, "1 photon\n$\\hbar k$", color="C0", fontsize=9, va="center")
ax[2].text(2.9, 1.75, "focus", ha="center", fontsize=9, color="C2")
ax[2].annotate("", xy=(2.85, 0), xytext=(2.85, 1), arrowprops=dict(arrowstyle="->", color="C2", lw=2))
ax[2].annotate("", xy=(2.95, 0), xytext=(2.95, -1), arrowprops=dict(arrowstyle="->", color="C2", lw=2))
ax[2].text(3.05, 0.0, "0 photons\n(2 body)", color="C2", fontsize=9, va="center")
ax[2].text(2.2, -2.1, "eight non-negative channel rates", ha="center", fontsize=8, color="C0")
ax[2].set_title("(c) the same $\\pm q$ stencil, factored\ninto whole-photon events")

fig.tight_layout()
fig.savefig(output_path(FIGNAME), dpi=150, bbox_inches="tight")
dp = docs_path(FIGNAME)
if dp:
    fig.savefig(dp, dpi=150, bbox_inches="tight")
print(f"\nfigure written to {output_path(FIGNAME)}"
      + (f" and {dp}" if dp else ""))

print("\nDone.")
