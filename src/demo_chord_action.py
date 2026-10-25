"""
Companion to docs/supplement/chord_action.md: the chord action as the one
relation between worlds, with bonds and banks derived from it.

Parts:
  S   Schiff and Poirier (2012), for the comparison of section 1: the quantum
      term of their action, Eq. (11), equals the Fisher form
      (hbar^2/8m) x''^2/x'^4 up to a total C-derivative; its Euler-Lagrange
      equation is their Eq. (10); with rho = 1/x' it is the Fisher information
      of rho; and Hall, Deckert and Wiseman's interworld potential U_N, divided
      by N, converges to it (Gaussian, worlds at the quantiles of rho).
  A   Proposition CA2: for N particles with masses m_i and any V, the action
      difference of two worlds at X + Y and X - Y is a part linear in Y plus
      minus the compensated residual R(X, Y); the Euler-Lagrange equation of
      the linear part in Y is Newton's law for X.  R is odd in Y and O(Y^3).
  B   Proposition CA3: the residual kernel built from R reproduces the Moyal
      residual.  V = x^4 (exact: one term) and V = -6 sech^2 x (series to
      order hbar^4), on a Gaussian test function in p.
  C   Corollary CA4: two worlds on one momentum row have a clock
      misalignment that winds at -[V(X + Y) - V(X - Y)]/hbar (WC8(a)).
  D   Proposition CA5 (bonds are derived): for pair potentials R is blind to
      the centre-of-mass direction of the chord, and a chord that moves only
      particle 1 reads only the terms of V that involve particle 1 (WC8.1);
      three particles, one-, two- and three-body terms.
  E   Remark CA7 (banks are a realisation): an integrate-and-fire counter and a
      Poisson process driven by the same rate have the same mean count; their
      variances differ.

Every number in the note comes from this script.  Run with
  WPMW_OUTPUT=<dir> PYTHONPATH=src python3 -u src/demo_chord_action.py
"""
import json

import numpy as np
import sympy as sp
from scipy.special import ndtri

from wpmwlib.wpmw_utils import output_path

results = {}


def part_S():
    print("=== Part S: Schiff-Poirier and Hall-Deckert-Wiseman ===")
    t, C, m, hb = sp.symbols("t C m hbar", positive=True)
    x = sp.Function("x")(t, C)
    V = sp.Function("V")
    x1, x2, x3, x4 = (sp.diff(x, C, k) for k in (1, 2, 3, 4))
    Q_sp = hb**2/(4*m)*(x3/x1**3 - sp.Rational(5, 2)*x2**2/x1**4)     # SP Eq. (11)
    Q_f = hb**2/(8*m)*x2**2/x1**4                                      # Fisher form
    r1 = sp.simplify(Q_sp - Q_f - sp.diff(hb**2/(4*m)*x2/x1**3, C))
    L = m*sp.diff(x, t)**2/2 - V(x) - Q_f
    EL = sp.euler_equations(L, [x], [t, C])[0].lhs
    eq10 = (m*sp.diff(x, t, 2) + sp.diff(V(x), x)
            + hb**2/(4*m)*(x4/x1**4 - 8*x3*x2/x1**5 + 10*x2**3/x1**6))  # SP Eq. (10)
    r2 = sp.simplify(EL + eq10)
    rho = 1/x1
    r3 = sp.simplify((sp.diff(rho, C)/x1)**2/rho*x1 - x2**2/x1**4)
    print(f"  Q_SP - Q_Fisher - d/dC[(hbar^2/4m) x''/x'^3]  = {r1}")
    print(f"  Euler-Lagrange of the Fisher form + SP Eq. (10) = {r2}")
    print(f"  (d_x rho)^2/rho dx/dC - x''^2/x'^4, rho = 1/x'  = {r3}")
    print("  HDW U_N / N for a Gaussian of unit width (hbar = m = 1); Fisher value 1/8:")
    hdw = {}
    for N in (100, 1000, 10000, 100000):
        xs = ndtri((np.arange(N) + 0.5)/N)
        g = 1/np.diff(xs)
        U = np.sum(np.diff(g)**2)/8
        hdw[N] = U/N
        print(f"    N = {N:6d}:  U_N/N = {U/N:.6f}   deficit {0.125 - U/N:.2e}")
    results["S"] = {"total_derivative": str(r1), "eq10": str(r2), "fisher": str(r3),
                    "hdw_UN_over_N": hdw}


def part_A():
    print("\n=== Part A: Proposition CA2, the split of the chord action ===")
    t = sp.symbols("t")
    N = 3
    ms = sp.symbols(f"m1:{N+1}", positive=True)
    Xs = [sp.Function(f"X{i+1}")(t) for i in range(N)]
    Ys = [sp.Function(f"Y{i+1}")(t) for i in range(N)]
    a = sp.symbols(f"a1:{N+1}")
    V = sp.Function("V")
    L = lambda q: sum(ms[i]*sp.diff(q[i], t)**2/2 for i in range(N)) - V(*q)
    plus = [Xs[i] + Ys[i] for i in range(N)]
    minus = [Xs[i] - Ys[i] for i in range(N)]
    dL = L(plus) - L(minus)
    gradV = [sp.Subs(sp.diff(V(*a), a[i]), a, Xs).doit() for i in range(N)]
    R = V(*plus) - V(*minus) - 2*sum(Ys[i]*gradV[i] for i in range(N))
    lin = (2*sum(ms[i]*sp.diff(Xs[i], t)*sp.diff(Ys[i], t) for i in range(N))
           - 2*sum(Ys[i]*gradV[i] for i in range(N)))
    r = sp.simplify(sp.expand(dL - lin + R))
    print(f"  N = 3, general V:  dL - [linear part] + R = {r}")
    EL = sp.euler_equations(lin, Ys, [t])
    newton = [sp.simplify(e.lhs/(-2) - (ms[i]*sp.diff(Xs[i], t, 2) + gradV[i]))
              for i, e in enumerate(EL)]
    print(f"  Euler-Lagrange in Y_i of the linear part = -2[m_i X_i'' + d_i V(X)]: "
          f"residuals {newton}")
    # odd in Y and O(Y^3): a scaling check on a concrete potential
    lam = sp.symbols("lambda")
    x, y = sp.symbols("x y")
    Rv = (-6*sp.sech(x + lam*y)**2) - (-6*sp.sech(x - lam*y)**2) \
        - 2*lam*y*sp.diff(-6*sp.sech(x)**2, x)
    ser = sp.series(Rv, lam, 0, 5).removeO()
    low = [sp.simplify(ser.coeff(lam, k)) for k in range(3)]
    print(f"  R for v = -6 sech^2: coefficients of lambda^0, lambda^1, lambda^2 = {low}")
    results["A"] = {"split": str(r), "newton": [str(v) for v in newton],
                    "low_orders": [str(v) for v in low]}


def residual_rate(Vfun, x0, s, hbar, p, ymax=40.0, ny=40001):
    """(d_t W)_R(x0, p) from the chord-action residual, for W(p) Gaussian of width s.

    (d_t W)_R = (1/(pi hbar)) int dy e^{-2ipy/hbar} (-i/hbar) R(x0, y) rho~(y),
    rho~(y) = int dp' W(p') e^{2ip'y/hbar} = exp(-2 s^2 y^2 / hbar^2).
    R is odd in y, so this is -(2/(pi hbar^2)) int_0^inf R sin(2py/hbar) rho~ dy.
    """
    y = np.linspace(0.0, ymax, ny)
    h = 1e-5
    dV = (Vfun(x0 + h) - Vfun(x0 - h))/(2*h)
    R = Vfun(x0 + y) - Vfun(x0 - y) - 2*y*dV
    rt = np.exp(-2*s**2*y**2/hbar**2)
    integrand = R[None, :]*np.sin(2*np.outer(p, y)/hbar)*rt[None, :]
    return -(2/(np.pi*hbar**2))*np.trapezoid(integrand, y, axis=1)


def part_B():
    print("\n=== Part B: Proposition CA3, the residual kernel is the Moyal residual ===")
    hbar, s, x0 = 1.0, 0.7, 0.4
    p = np.linspace(-3, 3, 13)
    z = sp.symbols("z")
    Wsym = sp.exp(-z**2/(2*s**2))/(sp.sqrt(2*sp.pi)*s)
    dW = [sp.lambdify(z, sp.diff(Wsym, z, k), "numpy") for k in range(8)]
    # quartic: exact one-term Moyal residual
    V4 = lambda q: q**4
    num = residual_rate(V4, x0, s, hbar, p)
    moyal = -(hbar**2/24)*24*x0*dW[3](p)
    e4 = np.max(np.abs(num - moyal))/np.max(np.abs(moyal))
    print(f"  V = x^4, x0 = {x0}: max |kernel - Moyal| / max |Moyal| = {e4:.2e}")
    # sech^2: series to hbar^4 and hbar^6
    xs = sp.symbols("xs")
    Vs = -6*sp.sech(xs)**2
    dV = [float(sp.diff(Vs, xs, k).subs(xs, x0)) for k in range(8)]
    # The Moyal series is asymptotic: it converges when chords are short (W broad
    # in p, rho~ narrow in y) and fails when they are long.  The kernel is exact.
    Vn = lambda q: -6/np.cosh(q)**2
    sech_errs = {}
    for s in (3.0, 0.7):
        p = np.linspace(-3*s, 3*s, 13)
        Wsym = sp.exp(-z**2/(2*s**2))/(sp.sqrt(2*sp.pi)*s)
        dW = [sp.lambdify(z, sp.diff(Wsym, z, k), "numpy") for k in range(8)]
        num = residual_rate(Vn, x0, s, hbar, p, ymax=40.0)
        t3 = -(hbar**2/24)*dV[3]*dW[3](p)
        t5 = (hbar**4/1920)*dV[5]*dW[5](p)
        t7 = -(hbar**6/322560)*dV[7]*dW[7](p)
        scale = np.max(np.abs(num))
        errs = [np.max(np.abs(num - t))/scale for t in (t3, t3 + t5, t3 + t5 + t7)]
        sech_errs[s] = errs
        print(f"  V = -6 sech^2 x, momentum width {s}: relative error of the Moyal series "
              f"truncated at hbar^2, hbar^4, hbar^6: {errs[0]:.2e}, {errs[1]:.2e}, {errs[2]:.2e}")
    results["B"] = {"quartic_rel_err": e4, "sech2_rel_err": sech_errs}


def part_C():
    print("\n=== Part C: Corollary CA4, misalignment on a common momentum row ===")
    t, m, hb, P = sp.symbols("t m hbar P", positive=True)
    X0, Y0 = sp.symbols("X0 Y0")
    V = sp.Function("V")
    Xt = X0 + P/m*t                                   # both worlds on one row
    xa, xb = Xt + Y0, Xt - Y0
    La = (m*sp.diff(xa, t)**2/2 - V(xa))
    Lb = (m*sp.diff(xb, t)**2/2 - V(xb))
    lever = P*(xb - xa)/hb                            # constant on a common row
    mu_dot = (La - Lb)/hb + sp.diff(lever, t)
    r = sp.simplify(hb*mu_dot + (V(xa) - V(xb)))
    print(f"  hbar dmu/dt + [V(X+Y) - V(X-Y)] = {r}")
    results["C"] = str(r)


def part_D():
    print("\n=== Part D: Proposition CA5, bonds and locality are derived ===")
    x = sp.symbols("x1:4")
    y = sp.symbols("y1:4")
    v = sp.Function("v")
    U = sp.Function("U")
    w = sp.Function("w")
    a = sp.symbols("a1:4")

    def residual(Vexpr):
        Vf = sp.Lambda(a, Vexpr)
        g = [sp.diff(Vexpr, a[i]).subs(dict(zip(a, x))) for i in range(3)]
        return (Vf(*[x[i] + y[i] for i in range(3)]) - Vf(*[x[i] - y[i] for i in range(3)])
                - 2*sum(y[i]*g[i] for i in range(3)))

    pair = sum(v(a[i] - a[j]) for i in range(3) for j in range(i + 1, 3))
    Rp = residual(pair + sum(a[i]**2 for i in range(3))/2)
    cm = sp.simplify(sum(sp.diff(Rp, y[i]) for i in range(3)).doit())
    print(f"  pair potential + trap: (sum_i d/dY_i) R = {cm}   (centre-of-mass direction dark)")
    full = U(a[0]) + U(a[1]) + pair + w(a[0], a[1], a[2])
    Vf = sp.Lambda(a, full)
    eta = sp.symbols("eta")
    dV = sp.expand(Vf(x[0] + eta, x[1], x[2]) - Vf(x[0] - eta, x[1], x[2]))
    read = sorted({f.func.__name__ + str(f.args) for f in dV.atoms(sp.Function)
                   if isinstance(f, sp.core.function.AppliedUndef)})
    print("  a chord on particle 1 alone reads:")
    for k in read:
        print("    ", k)
    results["D"] = {"cm_derivative": str(cm), "read_by_chord_on_1": read}


def part_E():
    print("\n=== Part E: Remark CA7, banks are a realisation of the rate ===")
    rng = np.random.default_rng(7)
    tgrid = np.linspace(0.0, 3.0, 3001)
    K = 2.0 + 1.5*np.sin(3*tgrid)                   # a rate along one world's history
    Lam = np.trapezoid(K, tgrid)
    n = 200000
    u = rng.random(n)                                # the bank's hidden initial phase
    iaf = np.floor(Lam + u)                          # integrate-and-fire count
    poi = rng.poisson(Lam, n)
    print(f"  integral of the rate = {Lam:.6f}")
    print(f"  integrate-and-fire: mean {iaf.mean():.4f}, variance {iaf.var():.4f}")
    print(f"  Poisson:            mean {poi.mean():.4f}, variance {poi.var():.4f}")
    results["E"] = {"Lambda": Lam, "iaf": [iaf.mean(), iaf.var()],
                    "poisson": [poi.mean(), poi.var()]}


if __name__ == "__main__":
    part_S()
    part_A()
    part_B()
    part_C()
    part_D()
    part_E()
    with open(output_path("chord_action_results.json"), "w") as f:
        json.dump(results, f, indent=1, default=str)
    print("\nwrote", output_path("chord_action_results.json"))
