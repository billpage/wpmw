#!/usr/bin/env python3
"""
Verification companion to the erratum in §0 of
``docs/algorithm/multi_body_extension.md``.

That specification's world-ensemble form (§6) gives each world its own
momentum jumps: for a mode ``V_q cos(k x + phi)`` of an external potential,

    p -> p + sgn(Gamma) hbar k   at rate |Gamma|,   Gamma = -(V_q/hbar) sin(k x + phi)   (§6.2)

and for a mode of a pair potential ``V_q cos(k (x_i - x_j) + phi)`` the
correlated jump ``(p_i, p_j) -> (p_i + s hbar k, p_j - s hbar k)`` with
``s = sgn(Gamma)`` (§6.3).  This script tests both rules against quantities
known exactly, on a ring of length ``L = 2 pi`` with ``hbar = m = 1`` and a
single mode ``k = 1``.

Parts
-----
A. Sign.  Ehrenfest's theorem is exact in quantum mechanics:
   ``d<p>/dt = <-V'(x)>``.  For a Gaussian packet it has a closed form at
   ``t = 0``.  The §6.2 rule as written gives the same magnitude with the
   opposite sign; with the sign reversed ("flipped") it agrees.  The §3.4 mesh
   form of the same document has the correct sign, so the world form
   contradicts its own mesh form.

B. Heating.  The exact Schroedinger evolution (split-operator reference)
   conserves ``<H>``.  Both versions of the §6.2 rule heat.  Because each jump
   is ``+-hbar k`` at rate ``|Gamma|``, the expected energy change obeys the
   exact identity

       d<H>/dt = (hbar k)^2 <|Gamma|> / (2m)  -  (1 - s) <p F(x)> / m ,

   with ``s = +1`` for the flipped rule and ``s = -1`` for the rule as written
   (``F = -V'``).  Integrating the measured right-hand side along the run
   reproduces the measured heating.  The first term is a momentum diffusion
   the quantum Liouville equation does not have: its single-mode term has
   zero second moment.

C. The pair rule.  Two particles on the ring with
   ``V_2 = V_q cos(k (x_1 - x_2))``.  Ehrenfest again fixes
   ``d<p_1>/dt`` at ``t = 0`` in closed form; the §6.3 rule as written gives
   the opposite sign, and conserves ``p_1 + p_2`` exactly, as it claims.

D. Negativity.  The flipped rule is a positive Markov jump process, so it
   maps ``W >= 0`` to ``W >= 0`` and can never create Wigner negativity from
   a Gaussian (Proposition T3 of
   ``docs/supplement/takabayasi_1954_stochastic_picture.md``).  The exact
   state at the end of the Part A/B run is measurably negative, even after a
   positivity-preserving smoothing in ``p`` that can only reduce negativity.

Figure: ``multibody_world_rule.png`` -- <p>(t) and <H>(t) for the exact
evolution, classical Liouville, and the §6.2 rule as written and flipped.
"""
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from wpmwlib.wpmw_utils import docs_path, output_path

HBAR = MASS = 1.0
L = 2 * np.pi
K = 2 * np.pi / L
V1 = 1.0
X0, P0, SIG = np.pi / 2, 0.0, 0.5          # packet centre where |F| is largest
T, DT = 6.0, 0.002
NSTEPS = int(round(T / DT))
REC = 50
NW = 200_000
RNG = np.random.default_rng(1)


def V(x):
    return V1 * np.cos(K * x)


def F(x):                                   # classical force -V'(x)
    return V1 * K * np.sin(K * x)


def gamma(x):                               # §6.2 rate field, signed
    return -(V1 / HBAR) * np.sin(K * x)


# ------------------------------------------------------------------ exact QM
def run_qm():
    ng = 1024
    xg = np.arange(ng) * L / ng
    pg = 2 * np.pi * HBAR * np.fft.fftfreq(ng, d=L / ng)
    d = (xg - X0 + L / 2) % L - L / 2
    psi = np.exp(-d**2 / (4 * SIG**2) + 1j * P0 * xg / HBAR)
    psi /= np.sqrt(np.sum(abs(psi)**2) * L / ng)
    ev = np.exp(-0.5j * V(xg) * DT / HBAR)
    et = np.exp(-1j * pg**2 / (2 * MASS) * DT / HBAR)
    out = []
    for n in range(NSTEPS + 1):
        if n % REC == 0:
            P = abs(np.fft.fft(psi))**2
            P /= P.sum()
            X = abs(psi)**2
            X /= X.sum()
            out.append(((P * pg).sum(),
                        (P * pg**2).sum() / (2 * MASS) + (X * V(xg)).sum()))
        if n < NSTEPS:
            psi = ev * np.fft.ifft(et * np.fft.fft(ev * psi))
    return np.array(out), xg, psi


def wigner_negativity(xg, psi, ell=0.6, stride=4):
    """Negativity of the Wigner function smoothed in p, a lower bound.

    W(x, p) = (1/pi hbar) int psi*(x+y) psi(x-y) exp(2ipy/hbar) dy, with the
    lag integrand multiplied by exp(-y^2 / 2 ell^2).  That is convolution in
    p with a positive Gaussian of width hbar/(2 ell), which cannot create
    negativity, so any negativity measured is genuine.  The taper also
    suppresses the ring's spurious lag y = L/2, where x + y and x - y
    coincide mod L (by a factor 1e-6 for ell = 0.6).
    """
    ng = len(xg)
    dx = xg[1] - xg[0]
    j = np.arange(-ng // 2, ng // 2)
    taper = np.exp(-(j * dx)**2 / (2 * ell**2))
    rows = []
    for i in range(0, ng, stride):
        c = taper * np.conj(psi[(i + j) % ng]) * psi[(i - j) % ng]
        # W(x, p) = (1/pi hbar) sum_j c_j exp(2 i p j dx / hbar) dx
        w = np.fft.fftshift(np.fft.fft(np.fft.ifftshift(c))).real
        rows.append(w * dx / (np.pi * HBAR))
    W = np.array(rows)
    dp = np.pi * HBAR / (ng * dx)            # conjugate spacing of y = j dx
    return np.abs(np.minimum(W, 0)).sum() * dx * stride * dp, W.min()


# ------------------------------------------------------------ world ensembles
X_INIT = X0 + SIG * RNG.standard_normal(NW)
P_INIT = P0 + HBAR / (2 * SIG) * RNG.standard_normal(NW)


def run_worlds(mode):
    """mode: 'CL' (Newtonian), 'SPEC' (§6.2 as written), 'FLIP' (sign reversed).

    Returns recorded (<p>, <H>) and the integrated right-hand side of the
    Part B identity.
    """
    x, p = X_INIT.copy(), P_INIT.copy()
    out, pred, acc = [], [], 0.0
    s = {"SPEC": -1.0, "FLIP": 1.0}.get(mode, 0.0)
    for n in range(NSTEPS + 1):
        if n % REC == 0:
            out.append((p.mean(), (p**2 / (2 * MASS) + V(x)).mean()))
            pred.append(acc)
        if n == NSTEPS:
            break
        if mode == "CL":
            p += 0.5 * DT * F(x)
            x += p / MASS * DT
            p += 0.5 * DT * F(x)
            continue
        g = gamma(x)
        acc += DT * ((HBAR * K)**2 * np.abs(g).mean() / (2 * MASS)
                     - (1 - s) * (p * F(x)).mean() / MASS)
        x += p / MASS * DT
        g = gamma(x)
        ev = RNG.random(NW) < np.abs(g) * DT
        sign = np.sign(g) if mode == "SPEC" else -np.sign(g)
        p += ev * sign * HBAR * K
    return np.array(out), np.array(pred)


def part_c(tau=0.05, nw=1_000_000):
    """Two particles, pair mode V_q cos(k (x1 - x2)), §6.3 rule as written."""
    vq = 1.0
    x1 = 0.0 + SIG * RNG.standard_normal(nw)
    x2 = -np.pi / 2 + SIG * RNG.standard_normal(nw)   # r0 = +pi/2
    p1 = HBAR / (2 * SIG) * RNG.standard_normal(nw)
    p2 = HBAR / (2 * SIG) * RNG.standard_normal(nw)
    ptot0 = p1 + p2
    p10 = p1.mean()
    for _ in range(int(round(tau / DT))):
        x1 += p1 * DT
        x2 += p2 * DT
        r = x1 - x2
        g = -(vq / HBAR) * np.sin(K * r)
        ev = RNG.random(nw) < np.abs(g) * DT
        kick = ev * np.sign(g) * HBAR * K
        p1 += kick
        p2 -= kick
    slope = (p1.mean() - p10) / tau
    # Ehrenfest: d<p1>/dt = <-dV2/dx1> = vq k <sin(k r)>, r ~ N(r0, 2 sig^2)
    exact = vq * K * np.exp(-K**2 * SIG**2) * np.sin(K * np.pi / 2)
    se = HBAR * K * np.sqrt(vq / HBAR / (nw * tau))   # Poisson floor, rough
    return slope, exact, se, np.max(np.abs(p1 + p2 - ptot0))


def main():
    t = np.arange(0, NSTEPS + 1, REC) * DT
    qm, xg, psi_T = run_qm()
    res = {m: run_worlds(m) for m in ("CL", "SPEC", "FLIP")}

    # ---- A. sign
    exact_slope = V1 * K * np.exp(-K**2 * SIG**2 / 2) * np.sin(K * X0)
    i1 = int(round(0.05 / (REC * DT))) or 1
    print("A. Initial Ehrenfest slope d<p>/dt at t = 0")
    se_a = HBAR * K * np.sqrt(V1 / HBAR / (NW * t[i1]))    # Poisson floor, rough
    print(f"   exact  {exact_slope:+.4f}   QM  {(qm[i1,0]-qm[0,0])/t[i1]:+.4f}"
          f"   (world-form noise ~{se_a:.3f})")
    for m in ("CL", "SPEC", "FLIP"):
        o = res[m][0]
        print(f"   {m:5s}  {(o[i1,0]-o[0,0])/t[i1]:+.4f}")
    print(f"   <p>(t = 1.2):  QM {qm[12,0]:+.4f}  CL {res['CL'][0][12,0]:+.4f}  "
          f"SPEC {res['SPEC'][0][12,0]:+.4f}  FLIP {res['FLIP'][0][12,0]:+.4f}")

    # ---- B. heating
    print("\nB. Energy change over T = %.1f" % T)
    print(f"   QM    {qm[-1,1]-qm[0,1]:+.2e}")
    print(f"   CL    {res['CL'][0][-1,1]-res['CL'][0][0,1]:+.2e}")
    for m in ("SPEC", "FLIP"):
        o, pr = res[m]
        dh = o[-1, 1] - o[0, 1]
        print(f"   {m:5s} measured {dh:+.4f}   identity {pr[-1]:+.4f}   "
              f"rel diff {abs(dh - pr[-1]) / abs(pr[-1]):.1e}")

    # ---- C. pair rule
    slope, exact, se, pviol = part_c()
    print("\nC. Pair rule §6.3, d<p1>/dt at t = 0 (window 0.05)")
    print(f"   exact {exact:+.4f}   §6.3 as written {slope:+.4f}   "
          f"(noise ~{se:.3f});  max |Δ(p1+p2)| per world {pviol:.1e}")

    # ---- D. negativity
    d = (xg - X0 + L / 2) % L - L / 2
    psi0 = np.exp(-d**2 / (4 * SIG**2))
    psi0 /= np.sqrt(np.sum(abs(psi0)**2) * L / len(xg))
    neg0, wmin0 = wigner_negativity(xg, psi0)
    neg, wmin = wigner_negativity(xg, psi_T)
    print("\nD. Negativity of the p-smoothed Wigner function (a lower bound)")
    print(f"   t = 0 (Gaussian, estimator floor):  {neg0:.1e}   min W {wmin0:+.1e}")
    print(f"   t = T (exact state):                {neg:.4f}   min W {wmin:+.4f}")

    # ---- figure
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
    for j, lab in enumerate((r"$\langle p\rangle$", r"$\langle H\rangle$")):
        ax[j].plot(t, qm[:, j], "k", lw=2, label="exact QM")
        for m, c, name in (("CL", "C0", "classical"),
                           ("SPEC", "C3", "§6.2 as written"),
                           ("FLIP", "C2", "§6.2 sign flipped")):
            ax[j].plot(t, res[m][0][:, j], c, label=name)
        ax[j].set(xlabel="t", ylabel=lab)
    ax[1].plot(t, res["FLIP"][0][0, 1] + res["FLIP"][1], "C2:", lw=2,
               label="identity (B)")
    ax[1].plot(t, res["SPEC"][0][0, 1] + res["SPEC"][1], "C3:", lw=2)
    ax[0].set_title("mean momentum")
    ax[1].set_title("energy")
    ax[0].legend(fontsize=8)
    ax[1].legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    name = "multibody_world_rule.png"
    fig.savefig(output_path(name), dpi=150, bbox_inches="tight")
    dp = docs_path(name)
    if dp:
        fig.savefig(dp, dpi=150, bbox_inches="tight")
    print("\nfigure:", output_path(name))


if __name__ == "__main__":
    main()
