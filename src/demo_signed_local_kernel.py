#!/usr/bin/env python3
r"""
Exact (SymPy) and spectral verification for the supplement
docs/supplement/poisson_kicks_and_pair_branching.md, Parts A-E (Theorem prefix X).

Question.  Can the quantum part of the Wigner generator be written as a *local*
momentum-jump process whose step and rate depend only on derivatives of V at the
particle's position?  Cyganski (call of 2026-09-29) proposed the one-particle
Poisson rule   dp = hbar |V''/V'| dN,   rate  V'^2 / V''.

Units: hbar and m are kept symbolic in A-C; numerical parts use hbar = m = 1.

Parts
  A  Kramers-Moyal targets.  A jump kernel with moments M_n reproduces the QLE
     iff M_{2s+1} = (-1)^{s+1} (hbar/2)^{2s} V^{(2s+1)},  M_{2s} = 0.
  B  The one-particle Poisson rule: M1 is right, M2 = hbar |V''| is spurious
     for every V (Proposition X1).
  B2 Kernel symbols on the density matrix: a positive kick damps coherence
     at rate lambda(1 - cos(a y/hbar)); a signed pair is a pure phase
     (Proposition X6).
  C  The signed two-point kernel (+w at +a, -w at -a): a^2 = -(hbar^2/4) V'''/V'.
     It is exact at every order iff V''' = r V' with r CONSTANT (Theorem X2):
     r = -k^2 sinusoid (a = hbar k/2), r = 0 quadratic (a = 0: the force),
     r = +kappa^2 exponential (a imaginary).
  D  Domain of validity and M5 error for V = omega^2 x^2/2 + V1 sin(k x)
     (the potential of Cyganski's notebook and of Part F-J).
  E  Dipole limit: F/(2a) [W(p+a) - W(p-a)] -> F dW/dp as a -> 0, error O(a^2),
     event rate O(1/a) (Proposition X3).  Spectral run on a harmonic cat state.
"""

import os
import sys

import numpy as np
import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from wpmwlib.wpmw_utils import output_path, docs_path  # noqa: E402

FIGNAME = "signed_local_kernel.png"
NMAX = 7  # highest moment checked


def banner(s):
    print("\n" + "=" * 72)
    print(s)
    print("=" * 72)


hb, kap = sp.symbols("hbar kappa", positive=True)
z = sp.Symbol("z")
c = sp.symbols("c1:%d" % (NMAX + 1), real=True)  # c_n = V^{(n)}(x)


def Vloc(zz):
    """Local Taylor model V(x+z) - V(x) with symbolic derivatives c_n."""
    return sum(c[n - 1] * zz ** n / sp.factorial(n) for n in range(1, NMAX + 1))


def jump_symbol(moments):
    """Symbol (action on e^{i kappa p}) of a jump generator with moments M_n."""
    return sum((-sp.I * kap) ** n * m / sp.factorial(n) for n, m in moments.items())


# ----------------------------------------------------------------------
banner("PART A  --  Kramers-Moyal targets from the exact Wigner symbol")
# exact symbol of the potential term:  i [V(x + hbar kappa/2) - V(x - hbar kappa/2)] / hbar
exact = sp.expand(sp.I * (Vloc(hb * kap / 2) - Vloc(-hb * kap / 2)) / hb)
Mreq = {}
for n in range(1, NMAX + 1):
    coef = exact.coeff(kap, n)
    Mreq[n] = sp.simplify(coef * sp.factorial(n) / (-sp.I) ** n)
    s = (n - 1) // 2
    formula = ((-1) ** (s + 1) * (hb / 2) ** (2 * s) * c[n - 1]) if n % 2 else 0
    print(f"  M{n} required = {Mreq[n]!s:32s}  formula residual = {sp.simplify(Mreq[n] - formula)}")

# ----------------------------------------------------------------------
banner("PART B  --  one-particle Poisson rule  step hbar|V''/V'|, rate V'^2/(hbar|V''|)")
res_B = {}
for s1 in (+1, -1):
    for s2 in (+1, -1):
        c1p, c2p = sp.symbols("f g", positive=True)  # |V'|, |V''|
        c1v, c2v = s1 * c1p, s2 * c2p
        a = -s1 * hb * c2p / c1p          # kick along the classical force -V'
        lam = c1p ** 2 / (hb * c2p)       # rate carries 1/hbar
        M = {n: sp.simplify(lam * a ** n) for n in range(1, 4)}
        res_B[(s1, s2)] = M
        print(f"  sgn V'={s1:+d}, sgn V''={s2:+d}:  M1 = {M[1]!s:6s} (target {-c1v})   "
              f"M2 = {M[2]!s:10s} (target 0)   M3 = {M[3]}  (target hbar^2 V'''/4)")
print("  => M1 = -V' always;  M2 = hbar |V''| > 0 always  (spurious momentum diffusion).")
print("  Any positive kernel has M2 >= M1^2 / (total rate) > 0 (Cauchy-Schwarz), so no")
print("  positive one-particle rule matches the QLE except where V' = 0 (no jumps at all).")

# ----------------------------------------------------------------------
banner("PART B2 --  kernel symbols: diffusion = decoherence (Proposition X6)")
# Action of a jump generator on the density matrix at leg separation y:
# W(p - D) <-> exp(-i D y / hbar) rho(y), so a kernel K has symbol
#   S(y) = int K(D) [exp(-i D y/hbar) - 1] dD      (loss term included).
# Re S < 0 damps the coherence between legs a distance y apart; Im S is a
# pure phase (unitary).
y, lam, a_, w_ = sp.symbols("y lambda a w", positive=True)
S_one = lam * (sp.exp(-sp.I * a_ * y / hb) - 1)              # single positive kick
S_pair = w_ * (sp.exp(-sp.I * a_ * y / hb) - sp.exp(sp.I * a_ * y / hb))  # signed pair
re_one = sp.simplify(sp.re(sp.expand_complex(S_one)))
re_pair = sp.simplify(sp.re(sp.expand_complex(S_pair)))
print(f"  one-particle kick:  Re S(y) = {re_one}")
print(f"      small y:        Re S(y) = {sp.series(re_one, y, 0, 5).removeO()}"
      f"   = -(M2/2) y^2/hbar^2 + ...,  M2 = lambda a^2")
print(f"  signed pair:        Re S(y) = {re_pair}   (pure phase: no decoherence)")
print(f"                      Im S(y) = {sp.simplify(sp.im(sp.expand_complex(S_pair)))}")
print("  => coherence between legs y apart decays at rate lambda(1 - cos(a y/hbar)),")
print("     i.e. momentum diffusion D = M2/2 is decoherence ~ exp(-D y^2 t / hbar^2).")
# Cauchy-Schwarz for a positive kernel: M3^2 <= M2 M4 ; equality for one spike
M = {n: lam * a_ ** n for n in (2, 3, 4)}
print(f"  one spike: M3^2 - M2*M4 = {sp.simplify(M[3]**2 - M[2]*M[4])}  (Cauchy-Schwarz equality)")

# ----------------------------------------------------------------------
banner("PART C  --  signed two-point kernel  (+w at +a, -w at -a)")
a2, w = sp.symbols("a2 w")
# moments: M_odd = 2 w a^n, M_even = 0.  Match M1 and M3.
a2_sol = sp.simplify(Mreq[3] / Mreq[1])
print(f"  a^2 = M3/M1 = {a2_sol}     2 w a = M1 = {Mreq[1]}")
for n in (5, 7):
    model = sp.simplify(Mreq[1] * a2_sol ** ((n - 1) // 2))
    print(f"  M{n}: model - required = {sp.factor(sp.simplify(model - Mreq[n]))}")
print("  exactness at s = 2 requires  V' V^(5) = (V''')^2.")

# Theorem X2: with u = V', u'' = r u and all higher conditions u^(2s) = r^s u
# force r' = 0.  Check the s = 3 condition when r is allowed to vary.
x = sp.Symbol("x")
u, r = sp.Function("u")(x), sp.Function("r")(x)
C = sp.Symbol("C")
# From s=1,2:  u''=ru  and  (r' u^2)' = 0  =>  r' = C/u^2
rules = {sp.Derivative(u, (x, 2)): r * u}


def D(expr):
    e = sp.diff(expr, x)
    e = e.subs(sp.Derivative(u, (x, 2)), r * u).subs(sp.Derivative(r, x), C / u ** 2)
    return sp.simplify(e)


d = [u, sp.diff(u, x)]
for _ in range(5):
    d.append(D(d[-1]))
print(f"  with u'' = r u, r' = C/u^2:  u'''' - r^2 u = {sp.simplify(d[4] - r**2 * u)}")
print(f"                               u^(6) - r^3 u = {sp.simplify(d[6] - r**3 * u)}")
print("  => s = 3 forces C = 0, i.e. r constant: Theorem X2 (only V''' = r V', r const).")
for name, Vx in [("V0 cos(kx)", sp.Symbol("V0") * sp.cos(sp.Symbol("k", positive=True) * x)),
                 ("V0 x^2", sp.Symbol("V0") * x ** 2),
                 ("V0 cosh(kx)", sp.Symbol("V0") * sp.cosh(sp.Symbol("k", positive=True) * x)),
                 ("V0 x^4", sp.Symbol("V0") * x ** 4)]:
    V1_, V3_ = sp.diff(Vx, x), sp.diff(Vx, x, 3)
    rr = sp.simplify(V3_ / V1_) if V1_ != 0 else 0
    step2 = sp.simplify(-hb ** 2 / 4 * rr)
    print(f"  {name:12s}  r = V'''/V' = {str(rr):12s}  a^2 = {step2}")

# ----------------------------------------------------------------------
banner("PART D  --  V = omega^2 x^2/2 + V1 sin(k x): domain of validity and M5 error")
om, V1, k = sp.symbols("omega V_1 k", positive=True)
Vq = om ** 2 * x ** 2 / 2 + V1 * sp.sin(k * x)
dV = [sp.diff(Vq, x, n) for n in range(8)]
a2x = sp.simplify(-hb ** 2 / 4 * dV[3] / dV[1])
relM5 = sp.simplify((-dV[1] * (dV[3] / dV[1]) ** 2 + dV[5]) / (-dV[5]))  # (model-req)/req
print(f"  a^2(x)          = {a2x}")
print(f"  dM5/M5_required = {sp.factor(relM5)}")
print("  split version: force for omega^2 x^2/2 (M_{n>=2} = 0 exactly) +")
print("  fixed +/- hbar k/2 kernel for the sine (exact at every order, Lemma T1):")
Mn_split = {n: (-1) ** ((n - 1) // 2 + 1) * (hb / 2) ** (n - 1) * sp.diff(V1 * sp.sin(k * x), x, n)
            for n in (3, 5, 7)}
for n, v in Mn_split.items():
    req = (-1) ** ((n - 1) // 2 + 1) * (hb / 2) ** (n - 1) * dV[n]
    print(f"    M{n}: split - required = {sp.simplify(v - req)}")

PAR = dict(omega=0.4, V1=0.5, k=1.0, hbar=1.0)
fa2 = sp.lambdify(x, a2x.subs({om: PAR["omega"], V1: PAR["V1"], k: PAR["k"], hb: 1}), "numpy")
frel = sp.lambdify(x, relM5.subs({om: PAR["omega"], V1: PAR["V1"], k: PAR["k"], hb: 1}), "numpy")
xs = np.linspace(-9, 9, 36001)
A2 = fa2(xs)
REL = frel(xs)
real_frac = np.mean(A2 >= 0)
print(f"\n  numbers at omega={PAR['omega']}, V1={PAR['V1']}, k={PAR['k']}, hbar=1, x in [-9, 9]:")
print(f"    fraction of x with real step (a^2 >= 0)       = {real_frac:.3f}")
mask = np.isfinite(REL) & (np.abs(xs) > 3)
print(f"    median |dM5/M5| for |x| > 3                     = {np.median(np.abs(REL[mask])):.3f}")
print(f"    median |dM5/M5| for |x| <= 1                    = "
      f"{np.median(np.abs(REL[np.isfinite(REL) & (np.abs(xs) <= 1)])):.3f}")

# ----------------------------------------------------------------------
banner("PART E  --  dipole limit: F/(2a)[W(p+a)-W(p-a)] -> F dW/dp")
aa, F = sp.symbols("a F", positive=True)
sym_dip = sp.I * F * sp.sin(kap * aa) / aa     # symbol of F/(2a)[W(p+a) - W(p-a)]
ser = sp.series(sym_dip, aa, 0, 5).removeO()
print(f"  symbol = i F sin(kappa a)/a = {sp.expand(ser)} + O(a^6)")
print("  a -> 0: i F kappa  = symbol of F dW/dp (classical force).  Leading error:")
print("  (F a^2 / 6) d^3W/dp^3, i.e. a spurious M3 = F a^2.  Births per body per unit time = |F|/a.")

# spectral run: harmonic oscillator, cat state, exact answer = rigid rotation
w0 = 1.0
Nx = Np = 256
Lx = 16.0
xg = (np.arange(Nx) - Nx // 2) * (Lx / Nx)
pg = (np.arange(Np) - Np // 2) * (Lx / Np)
X, P = np.meshgrid(xg, pg, indexing="xy")
kx = 2 * np.pi * np.fft.fftfreq(Nx, d=Lx / Nx)
kp = 2 * np.pi * np.fft.fftfreq(Np, d=Lx / Np)


def W_cat(Xa, Pa, x0=2.5):
    """Wigner function of an even cat state of two coherent states at +/- x0 (hbar=m=w0=1)."""
    g = lambda xx: np.exp(-(xx) ** 2 - Pa ** 2) / np.pi
    norm = 2 * (1 + np.exp(-x0 ** 2))
    return (g(Xa - x0) + g(Xa + x0) + 2 * np.exp(-Xa ** 2 - Pa ** 2) * np.cos(2 * Pa * x0) / np.pi) / norm


W0 = W_cat(X, P)
T = np.pi / 2  # quarter period
th = T
Xr, Pr = X * np.cos(th) - P * np.sin(th), X * np.sin(th) + P * np.cos(th)  # back-rotate
W_exact = W_cat(Xr, Pr)


def evolve(a, nsteps=400):
    W = W0.astype(complex)
    dt = T / nsteps
    Fx = -w0 ** 2 * xg[None, :]                         # force F(x) = -V'(x)
    mult = kp[:, None] if a == 0 else np.sin(kp[:, None] * a) / a
    # dW/dt = V' dW/dp = -F dW/dp -> multiplier exp(-i dt F mult) in p-Fourier space
    kick = np.exp(-1j * dt * Fx * mult)
    half = np.exp(-1j * kx[None, :] * P * dt / 2)
    for _ in range(nsteps):
        W = np.fft.ifft(np.fft.fft(W, axis=1) * half, axis=1)
        W = np.fft.ifft(np.fft.fft(W, axis=0) * kick, axis=0)
        W = np.fft.ifft(np.fft.fft(W, axis=1) * half, axis=1)
    return W.real


dxdp = (Lx / Nx) * (Lx / Np)
e0 = np.sqrt(np.sum((evolve(0.0) - W_exact) ** 2) * dxdp)
print(f"\n  cat state, quarter period, 256^2 grid, 400 Strang steps")
print(f"  a = 0 (force, spectral):  L2 error = {e0:.3e}  (time-step floor)")
avals = np.array([0.4, 0.2, 0.1, 0.05, 0.025])
errs = []
for a in avals:
    e = np.sqrt(np.sum((evolve(a) - W_exact) ** 2) * dxdp)
    errs.append(e)
    print(f"  a = {a:6.3f}:  L2 error = {e:.3e}   error/a^2 = {e / a**2:.4f}   "
          f"births per body per unit time |F|/a = {2.5 * w0**2 / a:6.1f} (x = 2.5)")
errs = np.array(errs)
slope = np.polyfit(np.log(avals[2:]), np.log(errs[2:]), 1)[0]
print(f"  fitted slope d log(err)/d log(a), a <= 0.1 = {slope:.4f}   (expected 2)")

# ----------------------------------------------------------------------
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
A2c = np.where(np.abs(A2) < 4, A2, np.nan)
ax[0].plot(xs, A2c, lw=1.2, color="C0", label=r"$a^2(x)=-\frac{\hbar^2}{4}V'''/V'$")
ax[0].axhline(PAR["k"] ** 2 / 4, color="C2", ls="--", lw=1, label=r"sine alone: $(\hbar k/2)^2$")
ax[0].axhline(0, color="k", lw=0.6)
ax[0].fill_between(xs, -4, 4, where=A2 < 0, color="C3", alpha=0.12, label="step imaginary")
ax[0].set_ylim(-1.5, 1.5); ax[0].set_xlabel("x"); ax[0].set_title("(a) local signed step, quadratic + sine")
ax[0].legend(fontsize=8, loc="lower left")
RELc = np.where(np.abs(REL) < 5, REL, np.nan)
ax[1].plot(xs, RELc, lw=1.2, color="C1")
ax[1].axhline(-1, color="k", ls=":", lw=1)
ax[1].axhline(0, color="k", lw=0.6)
ax[1].set_ylim(-3, 3); ax[1].set_xlabel("x")
ax[1].set_title(r"(b) relative $M_5$ error of the local kernel")
ax[1].text(0.03, 0.05, "split (force + fixed $\\pm\\hbar k/2$): error 0 at all orders",
           transform=ax[1].transAxes, fontsize=8)
ax[2].loglog(avals, errs, "o-", color="C0", label="signed two-point kernel")
ax[2].loglog(avals, errs[-1] * (avals / avals[-1]) ** 2, "k--", lw=1, label=r"$\propto a^2$")
ax[2].axhline(e0, color="C2", ls=":", label="force (a = 0)")
ax[2].set_xlabel("kick a"); ax[2].set_ylabel(r"$L^2$ error vs exact rotation")
ax[2].set_title("(c) dipole limit on a harmonic cat state")
ax[2].legend(fontsize=8)
fig.tight_layout()
fig.savefig(output_path(FIGNAME), dpi=150, bbox_inches="tight")
dp_ = docs_path(FIGNAME)
if dp_:
    fig.savefig(dp_, dpi=150, bbox_inches="tight")
print(f"\nwrote {output_path(FIGNAME)}")
