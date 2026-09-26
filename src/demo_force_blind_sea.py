"""The force-blind sea -- companion demo for ladder step 23.

Verifies ``docs/analysis/force_blind_sea.md``, which provisionally replaces
postulate (S) for aligned pairs by (S'): between events an aligned pair
moves inertially (dq/dt = p/m, dp/dt = 0) while its members' clocks still
wind at L = p^2/2m - V(q).  Free bodies still obey (S).

A  Proposition Q1, exactly: in a ledger whose events draw their absorptive
   and emissive counts from the bodies and the kernel alone, the body
   fields -- and so E, N and f -- are bitwise independent of how the sea
   moves.  Checked on the shared ledger (clamped and unclamped) and on the
   step 16 ledger.  The two exceptions: a rate throttled by the sea
   (Theorem S5) and a kernel weighted by the sea (Prop L2(a), see
   demo_contact_kernel.py --sea).
B  Proposition Q2, symbolically (SymPy): under (S') a clock stays on its
   row's plane wave up to the eikonal phase -int V dt/hbar, a same-row
   pair winds at exactly U, and the reader of Theorem L8 reads U + 2y
   dp_ref/dt with no drift term at any time; a motionless pair would need
   a Hamiltonian clock.  Under (S) a body's offset from its row's plane
   wave is -(1/hbar) int (H dt + x dp): zero drift for a uniform force,
   the first term at V''.
C  Proposition Q3: what (S') does to the sea itself -- the step 16 ledger
   under both motions: the worst cell with recombination (Part C of the
   step 16 demo), the absorptive trace (Part F), and the throttled rate
   of Theorem S5.
D  (--heavy) Theorem S9's attractor traces under both motions, compared
   element by element.

Theorem Y6's rerun uses the step 20 demo with conserved tags:

    PYTHONPATH=src python3 -u src/demo_dark_sea_and_identity.py \\
        --parts E --full --sea S|blind

and Prop L2(a) with the repaired channel masks:

    PYTHONPATH=src python3 -u src/demo_contact_kernel.py 30 --sea S|blind

Usage::

    PYTHONPATH=src python3 -u src/demo_force_blind_sea.py [--heavy]

About four minutes; --heavy adds about ten.
"""
import argparse
import contextlib
import io

import numpy as np
import sympy as sp

import demo_emission_and_absorption as dea
import demo_sea_population_equilibrium as d16


def banner(text):
    print("\n" + "=" * 72)
    print(text)
    print("=" * 72, flush=True)


def quiet(fn, *a, **k):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*a, **k)


# ----------------------------------------------------------------------
# A.  Proposition Q1: invariance, bitwise
# ----------------------------------------------------------------------
def shared_ledger_run(sea_force, clamp, nstep=300, dt=0.02):
    run = dea.Ledger(v0=1.0, a=1.0, n_r=192, r_half=24.0, n_p=64, dp=0.25)
    run.sea_force = sea_force
    e0 = dea.packet(run, r0=-8.0, p0=1.2, sr=2.0, sp=0.25)
    up, um, sea = run.prepare(e0, 6.0)
    for _ in range(nstep):
        up, um, sea = run.stream3(up, um, sea, .5 * dt)
        up, um, sea = run.channels(up, um, sea, dt, clamp=clamp)[:3]
        up, um, sea = run.stream3(up, um, sea, .5 * dt)
    return up, um, sea


def part_a():
    banner("A. Proposition Q1: the observable is blind to the sea's motion")
    print("   shared ledger (demo_emission_and_absorption), Eckart, p0 = 1.2,"
          " 300 steps")
    print(f"   {'clamp':>6} {'u+ bitwise':>11} {'u- bitwise':>11}"
          f" {'max |dS|/B':>11}")
    for clamp in (True, False):
        a = shared_ledger_run(True, clamp)
        b = shared_ledger_run(False, clamp)
        B = dea.B
        print(f"   {str(clamp):>6} {str(np.array_equal(a[0], b[0])):>11}"
              f" {str(np.array_equal(a[1], b[1])):>11}"
              f" {np.abs(a[2] - b[2]).max() / B:11.4f}")
    run = d16.Ledger()
    e0, _ = quiet(d16.part_b, run)
    out = {}
    for force in (True, False):
        run.sea_force = force
        out[force] = run.run_events(e0, 6.0, 0.01, absorb=True)
    a, b = out[True], out[False]
    print("\n   step 16 ledger, absorptive, T = 6, dt = 0.01")
    print(f"   E bitwise {np.array_equal(a['e'], b['e'])},"
          f"  N {a['N']:.10f} / {b['N']:.10f},  f {a['f']:.10f} / {b['f']:.10f}")
    print(f"   worst cell S/B {a['min_s']:+.4f} under (S),"
          f" {b['min_s']:+.4f} under (S')")
    print("\n   The sea enters an event only as its source (absorption) or")
    print("   sink (emission); the counts are set by bodies and kernel.")


# ----------------------------------------------------------------------
# B.  Proposition Q2: kinematics, symbolically
# ----------------------------------------------------------------------
def part_b():
    banner("B. Proposition Q2: the kinematics of (S) and (S'), with SymPy")
    t, m, hb = sp.symbols("t m hbar", positive=True)
    p, y = sp.symbols("p y", real=True)
    V = sp.Function("V")
    q = sp.Function("q")(t)
    th = sp.Function("theta")(t)
    # (S'): dq/dt = p/m, dp/dt = 0, hbar dtheta/dt = p^2/2m - V(q)
    rules = {q.diff(t): p / m, th.diff(t): (p**2 / (2 * m) - V(q)) / hb}
    off = hb * th - (p * q - p**2 * t / (2 * m))      # offset from the row wave
    print(f"   (S'): d/dt [hbar theta - (p q - p^2 t/2m)] = "
          f"{sp.simplify(off.diff(t).subs(rules))}")
    # a motionless pair keeps the row wave only with a Hamiltonian clock
    th0 = sp.Symbol("thetadot0")
    rest = (hb * th - (p * q - p**2 * t / (2 * m))).diff(t).subs(
        {q.diff(t): 0, th.diff(t): th0})
    print(f"   motionless pair: needs hbar dtheta/dt = "
          f"{sp.solve(sp.Eq(rest, 0), th0)[0] * hb}  (a Hamiltonian clock,"
          " not P1)")
    # same-row pair under (S') read by a reader with momentum P(t)
    xi, xj = sp.Function("x_i")(t), sp.Function("x_j")(t)
    thi, thj = sp.Function("theta_i")(t), sp.Function("theta_j")(t)
    P = sp.Function("P")(t)
    r2 = {xi.diff(t): p / m, xj.diff(t): p / m,
          thi.diff(t): (p**2 / (2 * m) - V(xi)) / hb,
          thj.diff(t): (p**2 / (2 * m) - V(xj)) / hb}
    mu = thi - thj + P * (xj - xi) / hb
    rate = sp.simplify(hb * mu.diff(t).subs(r2))
    print(f"   (S') pair, reader momentum P(t): hbar dmu/dt = {rate}")
    print("      -> inertial reader (P' = 0): U exactly; parent reader:"
          " U - 2y V'_eff; no drift term ever")
    # (S): offset identity and the uniform force
    x, pp = sp.Function("x")(t), sp.Function("p")(t)
    H = pp**2 / (2 * m) + V(x)
    L = pp**2 / (2 * m) - V(x)
    ident = sp.simplify(L - (pp * x).diff(t).subs(x.diff(t), pp / m)
                        + H + x * pp.diff(t))
    print(f"   (S): d(hbar theta - p x)/dt + H + x dp/dt = {ident}")
    F, x1, x2, p0, tau = sp.symbols("F x1 x2 p0 tau", real=True)

    def body(x0):
        pt = p0 + F * tau
        xt = x0 + p0 * tau / m + F * tau**2 / (2 * m)
        th_ = (p0 * x0 + sp.integrate(pt**2 / (2 * m) + F * xt,
                                      (tau, 0, t))) / hb
        return th_, xt.subs(tau, t), pt.subs(tau, t)
    (a1, xa, pa), (a2, xb, pb) = body(x1), body(x2)
    print(f"   (S), uniform force V = -F x: same row {sp.simplify(pa - pb) == 0},"
          f" row misalignment {sp.simplify(a1 - a2 - pa * (xa - xb) / hb)}")
    X, xs = sp.symbols("X x_s", real=True)
    a_ = sp.symbols("a0:5")
    Vs = lambda z: sum(a_[n] * (z - xs)**n / sp.factorial(n) for n in range(5))
    z = sp.Symbol("z")
    rate_S = -Vs(xs) + X * sp.diff(Vs(z), z).subs(z, xs)
    print(f"   (S), a body's offset rate + V(anchor) = "
          f"{sp.expand(rate_S + Vs(xs - X))}   [a_n = V^(n)(x)]")


# ----------------------------------------------------------------------
# C.  Proposition Q3: the sea itself
# ----------------------------------------------------------------------
def part_c():
    banner("C. Proposition Q3: what (S') does to the sea (step 16 ledger)")
    run = d16.Ledger()
    e0, _ = quiet(d16.part_b, run)
    B = d16.B
    print("   emissive, with recombination kappa (step 16 Part C), T = 6:"
          " worst cell S/B")
    print(f"   {'kappa':>7} {'(S)':>9} {'(S`)':>9}".replace("`", "'"))
    for kappa in (0.0, 20.0, 200.0, 2000.0):
        row = []
        for force in (True, False):
            run.sea_force = force
            row.append(run.run_mesh(e0, kappa, 6.0, 0.01)["min_s"])
        print(f"   {kappa:7.0f} {row[0]:9.4f} {row[1]:9.4f}")
    print("\n   absorptive (Theorem S8), dt = 0.01: running worst cell, then"
          " the worst cell at t")
    tr = {}
    for force in (True, False):
        run.sea_force = force
        o = run.run_events(e0, 18.0, 0.01, absorb=True, trace=True)
        tr[force] = o["trace"]
        r6 = run.run_events(e0, 6.0, 0.01, absorb=True)["min_s"]
        print(f"   {'(S) ' if force else '(S`)'.replace('`', chr(39))}"
              f" running worst to t = 6: {r6:+.4f}")
    print(f"   {'t':>6} {'(S)':>9} {'(S`)':>9}".replace("`", "'"))
    for t in (2.0, 6.0, 10.0, 14.0, 18.0):
        v = [tr[f][int(np.argmin(np.abs(tr[f][:, 0] - t))), 1]
             for f in (True, False)]
        print(f"   {t:6.1f} {v[0]:9.4f} {v[1]:9.4f}")
    print("\n   Theorem S5, rate throttled by the sea (Gamma -> Gamma S/B),"
          " T = 8: rel L2 of E against the QLE")
    ref = e0.copy()
    for _ in range(800):
        ref = run.qle_step(ref, 0.01)
    print(f"   {'kappa':>7} {'(S)':>9} {'(S`)':>9}".replace("`", "'"))
    for kappa in (200.0, 2000.0):
        row = []
        for force in (True, False):
            run.sea_force = force
            e = run.run_mesh(e0, kappa, 8.0, 0.01, live=True)["e"]
            row.append(np.linalg.norm(e - ref) / np.linalg.norm(ref))
        print(f"   {kappa:7.0f} {row[0]:9.3f} {row[1]:9.3f}")


# ----------------------------------------------------------------------
# D.  (heavy) Theorem S9's traces under both motions
# ----------------------------------------------------------------------
def part_d():
    banner("D. Theorem S9 under both motions (step 16 Part H traces)")
    run = d16.Ledger()
    e0, _ = quiet(d16.part_b, run)
    for rho in (1.0, 20.0):
        o = {}
        for force in (True, False):
            run.sea_force = force
            o[force] = d16.run_traced(run, e0, 8.0, 0.01, rho=rho)
        a, b = o[True]["trace"], o[False]["trace"]     # columns: t, f, N, S_min
        same = np.array_equal(a[:, 1:3], b[:, 1:3])
        print(f"   rho = {rho:5.1f}: f and N traces identical: {same};"
              f"  final worst cell S/B {a[-1, 3]:+.4f} (S), {b[-1, 3]:+.4f} (S')")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--heavy", action="store_true",
                    help="add part D (Theorem S9 traces, about ten minutes)")
    args = ap.parse_args()
    part_a()
    part_b()
    part_c()
    if args.heavy:
        part_d()


if __name__ == "__main__":
    main()
