"""
Companion to docs/supplement/jordan_lifting_and_pair_chemistry.md: the Jordan
lifting of a signed jump generator into a positive two-species process, its
minimal cost, and pair chemistry as a sink.  Units hbar = m = 1.

Parts:
  A   Theorem HJ1 (lifting), symbolic.  For any split L = A - B with A Metzler and
      B >= 0, the block generator [[A, B], [B, A]] acting on (u+, u-) gives
      dW/dt = L W and dS/dt = (A + B) S.  For Takabayasi's single-mode kernel
      the W kernel carries only odd moments, the S kernel only even ones.
  B   Theorem HJ2 (minimal cost), exact.  Three-mode kernel at one x on a
      momentum lattice: the Jordan lifting reproduces exp(tL) for W to
      round-off, S grows exactly at 2 x (negative off-diagonal mass), and 2000
      random positive splits never grow more slowly.  A random zero-column-sum
      Metzler-plus-signed generator checks the bound in general.
  C   The lifted process as a Monte Carlo, in phase space with streaming:
      V = V0 cos x, a Gaussian packet, worlds that stream freely and emit
      signed pairs at rate |Gamma(x)|.  The signed count is compared with a
      split-operator Schroedinger solution (density and the first two momentum
      moments); the world count is compared with the HJ2 rate integrated along
      the run.
  D   Proposition HJ5 (pair chemistry).  PP + NN <-> 2 PN with K = exp(eps/kT):
      closed-form equilibrium (SymPy), the K -> infinity limit is the minimal
      ensemble, the inflation asymptotics, the same for the contact sink's own
      reaction P + N <-> PN, a mass-action relaxation with its
      Horn-Jackson free energy, and the energetics of two real reactions.
  E   Proposition HJ4 (same generator, different process): FAM defocus against
      independent up and down actions, first and second moments.
  F   the figures.

Outputs go through wpmwlib.wpmw_utils.output_path (and docs_path for figures).
"""
import argparse

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.linalg import expm

from wpmwlib.wpmw_utils import docs_path, output_path


def rule(title):
    print("\n" + "=" * 78 + "\n" + title + "\n" + "=" * 78)


# ----------------------------------------------------------------------------
# Part A -- the lifting, symbolically
# ----------------------------------------------------------------------------
def part_A():
    rule("Part A. Theorem HJ1: the lifting, and where the moments go")
    A, B, up, um, W, S = sp.symbols("A B u_+ u_- W S")
    G = sp.Matrix([[A, B], [B, A]])
    du = G * sp.Matrix([up, um])
    sub = {up: (S + W) / 2, um: (S - W) / 2}
    dW = sp.simplify((du[0] - du[1]).subs(sub))
    dS = sp.simplify((du[0] + du[1]).subs(sub))
    print(f"  dW/dt = {dW}      dS/dt = {dS}")
    assert sp.simplify(dW - (A - B) * W) == 0 and sp.simplify(dS - (A + B) * S) == 0

    a, Gam = sp.symbols("a Gamma", positive=True)
    # Lemma T1: Lambda[f](p) = Gamma [f(p+a) - f(p-a)]; a world at p' sends
    # weight +Gamma to p'-a and -Gamma to p'+a.  Jumps q = p - p'.
    kW = {-a: Gam, +a: -Gam}            # J, for W
    kS = {-a: Gam, +a: Gam}             # |J|, for S
    print("  single mode, jump q = p - p' (Lemma T1 of the Takabayasi supplement)")
    for name, k in (("J   (W)", kW), ("|J| (S)", kS)):
        mom = [sp.simplify(sum(w * q**n for q, w in k.items())) for n in range(6)]
        print(f"    {name}: M0..M5 = {mom}")
    mW = [sp.simplify(sum(w * q**n for q, w in kW.items())) for n in range(6)]
    mS = [sp.simplify(sum(w * q**n for q, w in kS.items())) for n in range(6)]
    assert all(mW[n] == 0 for n in (0, 2, 4)) and all(mS[n] == 0 for n in (1, 3, 5))
    print("  even moments of J vanish (Takabayasi (3.17)); odd moments of |J| vanish;"
          " M2(|J|) = 2 a^2 Gamma > 0")


# ----------------------------------------------------------------------------
# Part B -- the minimal cost, exactly
# ----------------------------------------------------------------------------
def shift(N, m):
    """(shift(m) f)_n = f_{n-m}: moves weight from row n to row n + m (periodic)."""
    return np.roll(np.eye(N), m, axis=0)


def lifted(A, B):
    return np.block([[A, B], [B, A]])


def neg_mass(L):
    off = L - np.diag(np.diag(L))
    return np.minimum(off, 0).sum(axis=0) * -1.0   # per column


def part_B(seed=3):
    rule("Part B. Theorem HJ2: S grows at least at 2 x (negative off-diagonal mass)")
    N = 121
    # Three cosine modes at one x: Gamma_q = -(V_q/hbar) sin(phase_q), jumps +-q rows
    Vq = {1: 1.0, 2: 0.45, 3: -0.3}
    ph = {1: 0.7, 2: 2.1, 3: -1.0}
    L = np.zeros((N, N))
    for q, v in Vq.items():
        g = -v * np.sin(ph[q])
        L += g * (shift(N, -q) - shift(N, +q))          # +g to p'-q, -g to p'+q
    Lp, Lm = np.maximum(L, 0), np.maximum(-L, 0)          # Jordan parts (diag L = 0)
    rate = 2 * neg_mass(L)[N // 2]
    tv = sum(2 * abs(v * np.sin(ph[q])) for q, v in Vq.items())
    print(f"  ||J||_TV = {tv:.6f}   2 x negative mass per column = {rate:.6f}")
    u0 = np.zeros(N); u0[N // 2] = 1.0
    out = []
    for t in (0.5, 1.0, 2.0):
        u = expm(lifted(Lp, Lm) * t) @ np.concatenate([u0, 0 * u0])
        Wt, St = u[:N] - u[N:], (u[:N] + u[N:]).sum()
        Wex = expm(L * t) @ u0
        out.append((t, np.max(np.abs(Wt - Wex)), St, np.exp(rate * t)))
        print(f"  t={t}: max|W_lift - exp(tL) W0| = {out[-1][1]:.1e}   S = {St:.6f}"
              f"   exp(rate t) = {out[-1][3]:.6f}")
    # random positive splits A = L+ + G + D, B = L- + G + D (G >= 0 off-diagonal, D >= 0 diagonal)
    rng = np.random.default_rng(seed)
    worst = np.inf
    t = 1.0
    ratios = []
    for _ in range(2000):
        Gx = rng.exponential(0.05, (N, N)) * (rng.random((N, N)) < 0.03)
        np.fill_diagonal(Gx, 0)
        D = np.diag(rng.exponential(0.02, N) * (rng.random(N) < 0.2))
        A, B = Lp + Gx + D, Lm + Gx + D
        assert np.allclose(A - B, L)
        colrate = (A + B).sum(axis=0)[N // 2]
        worst = min(worst, colrate - rate)
        ratios.append(colrate / rate)
    print(f"  2000 random positive splits: min(column rate - HJ2 bound) = {worst:.3e}"
          f"   (bound attained only by the Jordan split)")
    # a general generator: Metzler (classical, positive) part plus a signed part
    M = rng.exponential(0.3, (N, N)) * (rng.random((N, N)) < 0.05)
    np.fill_diagonal(M, 0); M -= np.diag(M.sum(axis=0))
    Lg = M + L
    Ap = np.maximum(Lg, 0) - np.diag(np.maximum(np.diag(Lg), 0)) + np.diag(np.diag(Lg))
    Bm = np.maximum(-(Lg - np.diag(np.diag(Lg))), 0)
    assert np.allclose(Ap - Bm, Lg)
    print(f"  general L (Metzler + signed): column rates of the Jordan lift equal"
          f" 2 x negative mass to {np.max(np.abs((Ap + Bm).sum(axis=0) - 2 * neg_mass(Lg))):.1e};"
          f" Metzler part alone: max 2 x negative mass = {2 * neg_mass(M).max():.1e}")
    return dict(L=L, Lp=Lp, Lm=Lm, rate=rate, N=N)


# ----------------------------------------------------------------------------
# Part C -- the lifted process as a Monte Carlo, against Schroedinger
# ----------------------------------------------------------------------------
def schroedinger(x0, s0, p0, V0, T, nt, Lbox=8 * np.pi, M=4096):
    x = (np.arange(M) - M // 2) * (Lbox / M)
    k = 2 * np.pi * np.fft.fftfreq(M, d=Lbox / M)
    psi = np.exp(-(x - x0) ** 2 / (4 * s0**2) + 1j * p0 * x)
    psi /= np.sqrt(np.sum(abs(psi) ** 2) * (Lbox / M))
    dt = T / nt
    V = V0 * np.cos(x)
    eV, eK = np.exp(-0.5j * V * dt), np.exp(-0.5j * k**2 * dt)
    for _ in range(nt):
        psi = eV * np.fft.ifft(eK * np.fft.fft(eV * psi))
    rho = abs(psi) ** 2
    phi = np.fft.fft(psi)
    pk = abs(phi) ** 2; pk /= pk.sum()
    return x, rho, (pk * k).sum(), (pk * k**2).sum()


def part_C(n0=200_000, V0=1.0, T=1.5, dt=0.005, seed=7):
    rule("Part C. The lifted process in phase space: signed count against Schroedinger")
    rng = np.random.default_rng(seed)
    x0, s0, p0 = np.pi, 0.5, 0.6            # V0 cos x has its minimum at x = pi
    x = rng.normal(x0, s0, n0)
    p = rng.normal(p0, 1 / (2 * s0), n0)
    sg = np.ones(n0)
    G = 20                                  # independent lineage groups, for error bars
    lin = np.arange(n0) % G
    a = 0.5                                 # hbar k / 2 with k = 1
    nsteps = int(round(T / dt))
    lnN_pred, hist_N, hist_t, hist_pred = 0.0, [n0], [0.0], [0.0]
    for step in range(nsteps):
        gam = -V0 * np.sin(x)
        lnN_pred += 2 * np.abs(gam).mean() * dt
        fire = rng.random(x.size) < np.abs(gam) * dt
        if fire.any():
            xs, ps, ss, gs = x[fire], p[fire], sg[fire], np.sign(gam[fire])
            lin = np.concatenate([lin, lin[fire], lin[fire]])
            # child at p - a with sign s*sgn(Gamma), child at p + a with sign -s*sgn(Gamma)
            x = np.concatenate([x, xs, xs])
            p = np.concatenate([p, ps - a, ps + a])
            sg = np.concatenate([sg, ss * gs, -ss * gs])
        x = x + p * dt                      # free streaming (no force: V is wholly in J)
        hist_N.append(x.size); hist_t.append((step + 1) * dt); hist_pred.append(lnN_pred)
    N = x.size
    print(f"  worlds: {n0} -> {N}; ln(N/N0) = {np.log(N / n0):.4f};"
          f" HJ2 rate integrated along the run = {lnN_pred:.4f}")
    xs_, rho, mp, mp2 = schroedinger(x0, s0, p0, V0, T, nsteps)
    edges = np.linspace(-np.pi, 3 * np.pi, 161)
    dx = np.diff(edges)
    hg = np.stack([np.histogram(x[lin == g], bins=edges, weights=sg[lin == g])[0]
                   for g in range(G)]) / ((n0 / G) * dx)
    h = hg.mean(axis=0)
    se = hg.std(axis=0, ddof=1) / np.sqrt(G)     # between-lineage standard error
    hc = 0.5 * (edges[1:] + edges[:-1])
    rho_b = np.interp(hc, xs_, rho)
    l1 = np.sum(np.abs(h - rho_b) * dx)
    m = rho_b > 1e-3
    z = (h[m] - rho_b[m]) / se[m]
    print(f"  density at t = {T}: L1(MC - exact) = {l1:.4f}  (sampling-noise level"
          f" sqrt(2/pi) sum(se dx) = {np.sqrt(2 / np.pi) * np.sum(se * dx):.4f});"
          f"  rms z over {m.sum()} bins = {np.sqrt(np.mean(z**2)):.2f}")
    for name, f, ex in (("<p>  ", lambda q: q, mp), ("<p^2>", lambda q: q**2, mp2)):
        gm = np.array([np.sum(sg[lin == g] * f(p[lin == g])) / (n0 / G) for g in range(G)])
        print(f"  {name}: MC {gm.mean():+.4f} +- {gm.std(ddof=1) / np.sqrt(G):.4f}   exact {ex:+.4f}")
    print(f"  normalisation sum(sign)/N0 = {sg.sum() / n0:.6f} (exact: every event adds +1 and -1)")
    return dict(hc=hc, h=h, rho_b=rho_b, se=se, t=np.array(hist_t), N=np.array(hist_N) / n0,
                pred=np.exp(np.array(hist_pred)), T=T)


# ----------------------------------------------------------------------------
# Part D -- pair chemistry as a sink
# ----------------------------------------------------------------------------
def equilibrium(nP, nN, K):
    """PN count x with x^2 = K n_PP n_NN, n_PP = (nP - x)/2, n_NN = (nN - x)/2."""
    a_ = 1 - K / 4
    b_ = K * (nP + nN) / 4
    c_ = -K * nP * nN / 4
    if abs(a_) < 1e-15:
        return -c_ / b_
    r = np.roots([a_, b_, c_])
    r = r[np.isreal(r)].real
    r = r[(r >= -1e-12) & (r <= min(nP, nN) + 1e-12)]
    return float(r[0])


def part_D():
    rule("Part D. Proposition HJ5: PP + NN <-> 2 PN as a sink")
    x, K, nP, nN = sp.symbols("x K n_P n_N", positive=True)
    eq = sp.Eq(x**2, K * (nP - x) * (nN - x) / 4)
    sols = sp.solve(eq, x)
    print("  equilibrium x = n_PN solves x^2 = K (n_P - x)(n_N - x)/4; roots:")
    for s_ in sols:
        print("    ", sp.simplify(s_))
    # K -> infinity: the root in [0, min] tends to min(n_P, n_N)
    lim = [sp.limit(s_.subs({nP: 110, nN: 100}), K, sp.oo) for s_ in sols]
    print(f"  K -> oo (n_P = 110, n_N = 100): roots -> {lim}  (the admissible one is min = 100)")
    print("  eps/kT   n_PP     n_NN     n_PN     N = 2(n_PP+n_NN)   |W|   N-|W|   8 n_N^2/(K|W|)")
    rows = []
    for e in (0.0, 2.0, 5.0, 10.0, 15.0, 20.0):
        Kv = np.exp(e)
        xv = equilibrium(110.0, 100.0, Kv)
        pp, nn = (110 - xv) / 2, (100 - xv) / 2
        Nw = 2 * (pp + nn)
        rows.append((e, Nw - 10))
        print(f"  {e:5.1f}  {pp:7.3f}  {nn:7.4f}  {xv:8.3f}   {Nw:9.4f}        10   {Nw - 10:7.4f}   "
              f"{8 * 100**2 / (Kv * 10):.4f}")
    print("  balanced cell n_P = n_N = 100 (W = 0):   eps/kT   inflation N    4 n/sqrt(K)")
    bal = []
    for e in (5.0, 10.0, 15.0, 20.0):
        Kv = np.exp(e)
        xv = equilibrium(100.0, 100.0, Kv)
        Nw = 2 * (100 - xv)
        bal.append((e, Nw))
        print(f"                                           {e:5.1f}   {Nw:9.4f}    {4 * 100 / np.sqrt(Kv):.4f}")

    # the project's own contact sink, P + N <-> PN, with K' = pairs per (free P x free N) in one cell
    print("  P + N <-> PN (the contact sink's reaction): x = K' (n_P - x)(n_N - x)")
    print("   K'      inflation (110,100)   2 n_N/(K'|W|)    inflation (100,100)   2 sqrt(n/K')")
    for Kp in (1e1, 1e2, 1e3, 1e4):
        xs_ = []
        for a_, b_ in ((110.0, 100.0), (100.0, 100.0)):
            f = np.roots([Kp, -(Kp * (a_ + b_) + 1), Kp * a_ * b_])
            f = f[np.isreal(f)].real
            xs_.append(float(f[(f >= 0) & (f <= min(a_, b_) + 1e-12)][0]))
        inf1 = (110 - xs_[0]) + (100 - xs_[0]) - 10
        inf2 = 2 * (100 - xs_[1])
        print(f"  {Kp:7.0e}   {inf1:12.5f}        {2 * 100 / (Kp * 10):10.5f}       {inf2:12.5f}       {2 * np.sqrt(100 / Kp):8.5f}")

    # mass-action relaxation and its free energy (Horn-Jackson)
    kf, Kv = 1.0, np.exp(6.0)
    kb = kf / Kv

    def rhs(t, y):
        pp, nn, pn = y
        r = kf * pp * nn - kb * pn**2
        return [-r, -r, 2 * r]

    y0 = [40.0, 35.0, 5.0]
    sol = solve_ivp(rhs, (0, 2.0), y0, rtol=1e-11, atol=1e-12, dense_output=True)
    ts = np.linspace(0, 2.0, 400)
    Y = sol.sol(ts)
    nPt = 2 * Y[0] + Y[2]; nNt = 2 * Y[1] + Y[2]
    xe = equilibrium(nPt[0], nNt[0], Kv)
    ys = np.array([(nPt[0] - xe) / 2, (nNt[0] - xe) / 2, xe])
    Gf = np.sum(Y * (np.log(Y / ys[:, None]) - 1) + ys[:, None], axis=0)
    print(f"  relaxation at eps/kT = 6: n_P drift {np.ptp(nPt):.1e}, n_N drift {np.ptp(nNt):.1e},"
          f" W drift {np.ptp(nPt - nNt):.1e}")
    print(f"  free energy G: {Gf[0]:.4f} -> {Gf[-1]:.2e}; max dG/dt = {np.max(np.diff(Gf)):.2e} (<= 0)")
    print("  energetics, eps = 2 D_PN - D_PP - D_NN (bond enthalpies in eV):")
    for name, dpp, dnn, dpn in (("H2 + I2 -> 2 HI", 4.52, 1.54, 3.09),
                                ("C2 + O2 -> 2 CO", 6.2, 5.15, 11.16)):
        print(f"    {name}: D_PN deeper than both? {dpn > max(dpp, dnn)!s:5s}   eps = {2 * dpn - dpp - dnn:+.2f} eV"
              f"   eps/kT at 300 K = {(2 * dpn - dpp - dnn) / 0.025852:.1f}")
    return dict(rows=np.array(rows), bal=np.array(bal), ts=ts, Gf=Gf)


# ----------------------------------------------------------------------------
# Part E -- same generator, different process
# ----------------------------------------------------------------------------
def part_E():
    rule("Part E. Proposition HJ4: FAM defocus against independent up and down")
    r = sp.symbols("r", positive=True)
    Dv = sp.Matrix([1, -2, 1])               # (du_{n-1}, du_n, du_{n+1})
    up, dn = sp.Matrix([0, -1, 1]), sp.Matrix([1, -1, 0])
    print("  mean   FAM:", list(r * Dv), "  up+down:", list(r * up + r * dn))
    cf, cu = r * Dv * Dv.T, r * up * up.T + r * dn * dn.T
    print("  Cov(du_{n-1}, du_{n+1})  FAM:", cf[0, 2], "   up+down:", cu[0, 2])
    print("  Var(du_n)                FAM:", cf[1, 1], "   up+down:", cu[1, 1])
    assert list(r * Dv) == list(r * up + r * dn)


# ----------------------------------------------------------------------------
# Part F -- figures
# ----------------------------------------------------------------------------
def part_F(B, C, D):
    rule("Part F. figures")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
    a0 = ax[0]
    a0.fill_between(C["hc"], C["h"] - 2 * C["se"], C["h"] + 2 * C["se"], color="C0", alpha=0.25,
                    label="lifted process, ±2 s.e.")
    a0.plot(C["hc"], C["h"], color="C0", lw=1)
    a0.plot(C["hc"], C["rho_b"], "k--", lw=1.2, label="Schrödinger")
    a0.set_xlabel("x"); a0.set_ylabel("ρ(x)")
    a0.set_title(f"Signed count of a positive process, t = {C['T']}")
    a0.legend(fontsize=8)
    a1 = ax[1]
    a1.semilogy(C["t"], C["N"], label="worlds N/N₀ (measured)")
    a1.set_xlabel("t"); a1.set_ylabel("N/N₀")
    a1.semilogy(C["t"], C["pred"], "k:", lw=1.5, label="exp(∫ ⟨2|Γ(x)|⟩ dt), HJ2 rate")
    a1.set_title("Cost: S grows at 2 × negative mass")
    a1.legend(fontsize=8)
    a2 = ax[2]
    rows, bal = D["rows"], D["bal"]
    a2.semilogy(rows[:, 0], rows[:, 1], "o-", label="|W| = 10 of 210")
    big = rows[:, 0] >= 5
    a2.semilogy(rows[big, 0], 8 * 100**2 / (np.exp(rows[big, 0]) * 10), "k:", lw=1,
                label="8 n_N² / (K |W|)")
    a2.semilogy(bal[:, 0], bal[:, 1], "s-", label="W = 0 (200 constituents)")
    a2.semilogy(bal[:, 0], 4 * 100 / np.sqrt(np.exp(bal[:, 0])), "k--", lw=1, label="4 n / √K")
    a2.set_xlabel("ε / kT"); a2.set_ylabel("inflation N − |W|")
    a2.set_title("PP + NN ⇌ 2PN: K → ∞ is the contact sink")
    a2.legend(fontsize=8)
    fig.tight_layout()
    name = "jordan_lifting.png"
    fig.savefig(output_path(name), dpi=130)
    dp = docs_path(name)
    if dp:
        fig.savefig(dp, dpi=130)
    print(f"  wrote {name}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--n0", type=int, default=200_000)
    ap.add_argument("--no-fig", action="store_true")
    args = ap.parse_args()
    part_A()
    B = part_B()
    C = part_C(n0=args.n0)
    D = part_D()
    part_E()
    if not args.no_fig:
        part_F(B, C, D)


if __name__ == "__main__":
    main()
