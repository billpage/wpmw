"""Darkness and identity in the compensated model.

Companion demo for ``docs/analysis/dark_sea_and_worldline_identity.md``
(ladder step 20).  Five parts, each verifying one claim of the note.

A  Theorem Y1 -- the compensated kernel is read from V directly, so the
   phase-blind no-go (Theorem 2 of the phase-resonance note) does not
   reach it.
B  Theorem Y2 -- what the sea's ineligibility for recombination buys.
   All realisations move E identically; sea-eligible ones are
   count-identical to the excluded hop.
C  Theorem Y3 -- the two bodies a catalysed recombination consumes are a
   winding pair whose midpoint is the parent's row.
D  Theorem Y4 -- a momentum kink of +-xi carries a phase ramp pinned at
   the kink point; a node-located kink and a phase reset at the parent's
   position are the same freedom.
E  Theorem Y6 -- tagged positons under the piecewise reading: individual
   crossing of the Eckart barrier is a momentum diffusion, not tunnelling.

Part E is the expensive one.  By default it runs one packet on a coarse
lattice (about a minute).  ``--full`` reproduces the note's table: three
packets on the fine lattice, about five minutes each.

Usage::

    PYTHONPATH=src python3 -u src/demo_dark_sea_and_identity.py [--full]
"""
import argparse
import numpy as np

from demo_emission_and_absorption import Ledger, packet, HBAR, MU

EPS = 1e-300


def banner(text):
    print("\n" + "=" * 72)
    print(text)
    print("=" * 72)


# ----------------------------------------------------------------------
# A.  The kernel is linear in V, and read from V directly
# ----------------------------------------------------------------------
def part_a():
    banner("A. Theorem Y1: the compensated kernel reads V directly")
    base = Ledger(v0=1.0)
    print(f"   {'lambda':>8} {'rel |K(lam V) - lam K(V)|':>26}"
          f" {'Gamma_max ratio':>17}")
    for lam in (1e-3, 0.1, 2.0, -1.0):
        run = Ledger(v0=lam)
        err = (np.abs(run.k - lam * base.k).max()
               / np.abs(lam * base.k).max())
        g = run.gamma_tot.max() / base.gamma_tot.max()
        print(f"   {lam:8.3f} {err:26.2e} {g:17.6f}")
    print()
    print("   The event rate Gamma = sum_q |K_q(x)| is first order in V with")
    print("   no phase variable on any body.  The hypothesis of the")
    print("   phase-blind no-go -- that V reaches the rates only through")
    print("   phase-blind data, which are V-independent at first order --")
    print("   therefore fails here.  This is the 'pinned' exception the")
    print("   no-go names: the vertex consults the potential itself.")


# ----------------------------------------------------------------------
# B.  What the sea's ineligibility buys
# ----------------------------------------------------------------------
def part_b():
    banner("B. Theorem Y2: books for the four ways to settle the legs")
    p, xi = 5, 2                       # rows -1, 0, +1 are p-xi, p, p+xi
    mom = {-1: p - xi, 0: p, 1: p + xi}

    def start():
        one = {r: 3 for r in mom}
        return dict(one), dict(one), dict(one)

    def books(up, um, S):
        E = {r: up[r] - um[r] for r in mom}
        N = sum(up.values()) + sum(um.values())
        Ssum = sum(S.values())
        return (E, N, Ssum, Ssum + N / 2,
                sum(mom[r] * (up[r] + um[r] + 2 * S[r]) for r in mom))

    def apply(kind):
        up, um, S = start()
        if kind == "ionisation":
            S[0] -= 1
            up[1] += 1
            um[-1] += 1
        elif kind == "hop (excluded)":
            um[1] -= 1
            um[-1] += 1
        else:
            leg_hi, leg_lo = kind
            if leg_hi == "free":                 # negaton at p + xi
                um[1] -= 1
            else:                                # break an aligned pair
                S[1] -= 1
                up[1] += 1
            if leg_lo == "free":                 # positon at p - xi
                up[-1] -= 1
            else:
                S[-1] -= 1
                um[-1] += 1
            S[0] += 1
        return up, um, S

    E0, N0, S0, P0, M0 = books(*start())
    print(f"   {'realisation':24s} {'dE(p-xi,p,p+xi)':>18s} {'dN':>4s}"
          f" {'dS':>4s} {'dP':>5s} {'dMom':>5s}")
    for k in [("free", "free"), ("sea", "sea"), ("free", "sea"),
              ("sea", "free"), "ionisation", "hop (excluded)"]:
        E, N, S, P, M = books(*apply(k))
        dE = tuple(E[r] - E0[r] for r in (-1, 0, 1))
        name = (f"recombine[{k[0]},{k[1]}]" if isinstance(k, tuple) else k)
        print(f"   {name:24s} {str(dE):>18s} {N - N0:>4d} {S - S0:>4d}"
              f" {P - P0:>5.1f} {M - M0:>5d}")
    print()
    print("   Every row moves E identically, so E cannot be what excludes")
    print("   the sea.  With the sea eligible, the mixed rows have")
    print("   dN = dS = 0 and move one negaton by 2 xi -- the excluded hop,")
    print("   and unlike the hop they conserve momentum, paying with one")
    print("   aligned pair shifted a row.  P = S + N/2 is blind to all of")
    print("   it.  With the sea ineligible every event has |dS| = 1.")


# ----------------------------------------------------------------------
# C.  The consumed bodies are a winding pair on the parent's row
# ----------------------------------------------------------------------
def _amp(eps, pj, th, x, t, xs, m):
    return eps * np.exp(1j * ((pj * (x - xs) - pj ** 2 / (2 * m) * t) / HBAR
                              + th))


def part_c():
    banner("C. Theorem Y3: the consumed bodies, as a pair")
    m, p, xi, xs = 1.0, 5.0, 2.0, 0.3
    pa, pb = p - xi, p + xi            # positon below, negaton above
    rng = np.random.default_rng(1)
    tha, thb = rng.uniform(0, 2 * np.pi, 2)
    x = np.linspace(-4, 4, 4001)
    for t in (0.0, 0.7):
        psi = (_amp(+1, pa, tha, x, t, xs, m)
               - _amp(+1, pb, thb, x, t, xs, m))
        mu_u = (((pa - pb) * (x - xs) - (pa ** 2 - pb ** 2) / (2 * m) * t)
                / HBAR + tha - thb)
        err = np.abs(np.abs(psi) - 2 * np.abs(np.sin(mu_u / 2))).max()
        nodes = x[1:][np.diff(np.floor(mu_u / (2 * np.pi))) != 0]
        print(f"   t = {t:.1f}:  max | |Psi| - 2|sin(mu/2)| | = {err:.2e};"
              f"  node spacing {np.diff(nodes).mean():.5f}"
              f"  (h/|dp| = {2 * np.pi * HBAR / abs(pa - pb):.5f})")
    print(f"   envelope drift (pa+pb)/2m = {(pa + pb) / (2 * m):.4f}"
          f" = parent velocity p/m = {p / m:.4f}")
    mu_star = (tha - thb) % (2 * np.pi)
    print(f"   mu at the vertex {mu_star:.4f} rad,"
          f"  |Psi(x*)| = {2 * abs(np.sin(mu_star / 2)):.4f}")
    for t in (0.0, 0.7, 5.0):
        psi = (_amp(+1, p, tha, x, t, xs, m)
               - _amp(+1, p, thb, x, t, xs, m))
        print(f"   bound at p, phases continuous, t = {t:.1f}:"
              f"  |Psi| in [{np.abs(psi).min():.4f}, {np.abs(psi).max():.4f}]")
    print("   -> uniformly and permanently gray unless mu(x*) = 0.")
    print()
    print("   Kinematics of the two absorptions:")
    ke = lambda ps: sum(q * q for q in ps) / (2 * m)
    for name, a, b in (("ledger", [p - xi, p + xi, p], [p, p, p]),
                       ("swap", [pb, pa, pb], [pa, pb, pb])):
        print(f"   {name:7s} sum p {sum(a):.0f} -> {sum(b):.0f};"
              f"  sum p^2/2m {ke(a):.1f} -> {ke(b):.1f};"
              f"  multiset preserved: {sorted(a) == sorted(b)}")
    print(f"   the ledger vertex loses xi^2/m = {xi ** 2 / m:.1f} of unsigned"
          " kinetic energy; the swap loses none.")


# ----------------------------------------------------------------------
# D.  The phase ramp that accompanies a kink
# ----------------------------------------------------------------------
def part_d():
    banner("D. Theorem Y4: a kink carries a phase ramp pinned at the kink")
    m, p, xi, tk, x_star = 1.0, 5.0, 2.0, 0.0, 0.3
    pa, pb = p - xi, p + xi
    rng = np.random.default_rng(1)
    tha, thb = rng.uniform(0, 2 * np.pi, 2)
    x = np.linspace(-4, 4, 4001)

    def phi(th, pj, xx, t):
        return th + (pj * xx - pj ** 2 / (2 * m) * t) / HBAR

    def after(th, pj, xk, xx, t):
        return (phi(th, pj, xk, tk)
                + (p * (xx - xk) - p ** 2 / (2 * m) * (t - tk)) / HBAR)

    mu = lambda xx: phi(tha, pa, xx, tk) - phi(thb, pb, xx, tk)
    n = np.round(mu(x_star) / (2 * np.pi))
    # d(mu)/dx = (pa - pb) / hbar = -2 xi / hbar
    x_node = x_star + (2 * np.pi * n - mu(x_star)) * HBAR / (-2 * xi)
    d = after(tha, pa, x_star, x, tk) - phi(tha, pa, x, tk)
    print(f"   ramp slope {np.polyfit(x, d, 1)[0]:.6f} (xi/hbar ="
          f" {xi / HBAR:.1f}),  value at the pivot {np.interp(x_star, x, d):.1e}")
    for name, xk in (("pivot at x*  ", x_star), ("pivot at node", x_node)):
        for t in (0.0, 3.0):
            psi = (np.exp(1j * after(tha, pa, xk, x, t))
                   - np.exp(1j * after(thb, pb, xk, x, t)))
            print(f"   {name}  x_k = {xk:+.4f}  t = {t:.1f}:"
                  f"  |Psi| in [{np.abs(psi).min():.2e},"
                  f" {np.abs(psi).max():.2e}]")
    reset = xi * (x_star - x_node) / HBAR
    print(f"   node pivot == pivot at x* plus a reset of {reset:+.4f} rad per"
          f" body;\n   twice that is {(-2 * reset) % (2 * np.pi):.4f}"
          f" = mu(x*) = {mu(x_star) % (2 * np.pi):.4f}")


# ----------------------------------------------------------------------
# E.  Tagged positons: who crosses the barrier?
# ----------------------------------------------------------------------
def channels_tagged(run, up, um, sea, tu, ts, dt):
    """run.channels() with tagged positons: tu free, ts inside a pair.

    Tags ride by proportional allocation -- the exact expectation of the
    tagged count when the bodies of a cell are exchangeable.  The fields
    are NOT clamped: clamping small negative populations breaks E (see
    the note, section 8).
    """
    n_abs = n_emi = 0.0
    for q in range(1, run.n_p // 2):
        lam = np.abs(run.k[:, q])[:, None]
        if lam.max() < 1e-14:
            continue
        sg = np.sign(run.k[:, q])[:, None]
        for parent, sp in ((up, 1.0), (um, -1.0)):
            D = lam * parent * dt
            if D.max() <= 0.0:
                continue
            t = np.broadcast_to(sg * sp, D.shape)
            capA = np.clip(np.where(t > 0, np.roll(um, -q, axis=1),
                                    np.roll(up, -q, axis=1)), 0.0, None)
            capB = np.clip(np.where(t > 0, np.roll(up, q, axis=1),
                                    np.roll(um, q, axis=1)), 0.0, None)
            A = np.minimum(D, np.minimum(capA, capB))
            Em = D - A
            n_abs += float(A.sum())
            n_emi += float(Em.sum())
            aq, am = np.roll(A, q, axis=1), np.roll(A, -q, axis=1)
            eq, em = np.roll(Em, q, axis=1), np.roll(Em, -q, axis=1)
            # tags, from the fields as they stand before this update
            fr_u = np.clip(tu / np.maximum(up, EPS), 0.0, 1.0)
            fr_s = np.clip(ts / np.maximum(sea, EPS), 0.0, 1.0)
            take_lo = np.where(t > 0, am, 0.0) * fr_u
            take_hi = np.where(t > 0, 0.0, aq) * fr_u
            ion = Em * fr_s
            tu = tu - take_lo - take_hi
            ts = (ts + np.roll(take_lo, q, axis=1)
                  + np.roll(take_hi, -q, axis=1) - ion)
            tu = tu + np.where(t > 0, np.roll(ion, q, axis=1),
                               np.roll(ion, -q, axis=1))
            # fields, exactly as Ledger.channels but without the clamp
            um -= np.where(t > 0, aq, 0.0)
            up -= np.where(t > 0, 0.0, aq)
            up -= np.where(t > 0, am, 0.0)
            um -= np.where(t > 0, 0.0, am)
            sea += A
            up += np.where(t > 0, eq, 0.0)
            um += np.where(t > 0, 0.0, eq)
            um += np.where(t > 0, em, 0.0)
            up += np.where(t > 0, 0.0, em)
            sea -= Em
            tu = np.clip(tu, 0.0, np.maximum(up, 0.0))
            ts = np.clip(ts, 0.0, np.maximum(sea, 0.0))
    return up, um, sea, tu, ts, n_abs, n_emi


def t_exact(p0, sp, v0=1.0, a=1.0):
    """Packet-averaged exact transmission of V0 sech^2(r/a)."""
    lam = 8.0 * MU * v0 * a ** 2 / HBAR ** 2
    p = np.linspace(0.01, 4.0, 4000)
    s = np.sinh(np.pi * np.sqrt(2 * MU * p ** 2 / (2 * MU)) * a / HBAR) ** 2
    c = np.cosh(np.pi / 2 * np.sqrt(lam - 1)) ** 2
    w = np.exp(-(p - p0) ** 2 / (2 * sp ** 2))
    w /= w.sum()
    return float((w * (s / (s + c))).sum())


def part_e(full=False):
    banner("E. Theorem Y6: tagged positons at the Eckart barrier")
    if full:
        run = Ledger(v0=1.0, a=1.0, n_r=256, r_half=40.0, n_p=128, dp=0.125)
        cases, r0, sr, sp, dt = (1.0, 1.2, 1.7), -12.0, 2.0, 0.25, 0.02
    else:
        run = Ledger(v0=1.0, a=1.0, n_r=192, r_half=24.0, n_p=64, dp=0.25)
        cases, r0, sr, sp, dt = (1.2,), -8.0, 2.0, 0.25, 0.02
    right = run.r[:, None] > 0.0
    H = (run.p[None, :] ** 2 / (2 * MU)
         + run.v0 / np.cosh(run.r[:, None] / run.a) ** 2)
    trans = lambda w: float(w[np.broadcast_to(right, w.shape)].sum() / w.sum())
    print(f"   barrier-top momentum p_b = {np.sqrt(2 * MU * run.v0):.4f},"
          f"  lattice dp = {run.dp}")
    print(f"   {'p0':>5} {'E0/V0':>6} {'T_cl':>7} {'T_exact':>8} {'T_mesh':>8}"
          f" {'T_E':>8} {'T_tag':>7} {'tag H>V0':>9} {'<H>/E0':>7}"
          f" {'kinks':>6} {'f':>6} {'|E-mesh|':>9}")
    for p0 in cases:
        e0 = packet(run, r0=r0, p0=p0, sr=sr, sp=sp)
        nstep = int(round(2.0 * abs(r0) / (p0 / MU) / dt))
        up, um, sea = run.prepare(e0, 6.0)
        tu, ts = up.copy(), np.zeros_like(sea)
        tag0 = tu.sum()
        ecl, emesh = e0.copy(), e0.copy()
        na = ne = kinks = 0.0
        for _ in range(nstep):
            up, um, sea = run.stream3(up, um, sea, .5 * dt)
            tu, ts = run.stream(tu, .5 * dt), run.stream(ts, .5 * dt)
            u0, s0 = tu.copy(), ts.copy()
            up, um, sea, tu, ts, a, em = channels_tagged(
                run, up, um, sea, tu, ts, dt)
            kinks += 0.5 * (np.abs(tu - u0).sum() + np.abs(ts - s0).sum())
            up, um, sea = run.stream3(up, um, sea, .5 * dt)
            tu, ts = run.stream(tu, .5 * dt), run.stream(ts, .5 * dt)
            ecl = run.stream(ecl, dt)
            emesh = run.qle_step(emesh, dt)
            na += a
            ne += em
        E, tag, E0 = up - um, tu + ts, p0 ** 2 / (2 * MU)
        print(f"   {p0:5.2f} {E0 / run.v0:6.2f} {trans(ecl):7.4f}"
              f" {t_exact(p0, sp, run.v0, run.a):8.4f} {trans(emesh):8.4f}"
              f" {trans(E):8.4f} {trans(tag):7.4f}"
              f" {tag[H > run.v0].sum() / tag.sum():9.4f}"
              f" {(tag * H).sum() / tag.sum() / E0:7.3f}"
              f" {kinks / tag0:6.3f} {na / (na + ne):6.3f}"
              f" {np.linalg.norm(E - emesh) / np.linalg.norm(emesh):9.3f}")
    print()
    xi = np.abs(run.xi)
    dp_diff = 0.5 * (xi[None, :] ** 2 * np.abs(run.k)).sum(axis=1)
    m1 = (np.abs((run.xi[None, :] * run.k).sum(axis=1)).max()
          / run.gamma_tot.max())
    i0 = int(np.argmin(np.abs(run.r)))
    print(f"   max |sum_q xi_q K_q| / max Gamma = {m1:.1e}: no mean kick.")
    print(f"   Gamma(summit) / max Gamma = "
          f"{run.gamma_tot[i0] / run.gamma_tot.max():.1e}: no kink at the top.")
    print(f"   {'r':>6} {'Gamma(r)':>10} {'D_p = 1/2 sum xi^2 |K|':>23}")
    for rr in (-4.0, -2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0):
        i = int(np.argmin(np.abs(run.r - rr)))
        print(f"   {run.r[i]:+6.2f} {run.gamma_tot[i]:10.4f}"
              f" {dp_diff[i]:23.4f}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--full", action="store_true",
                    help="part E: the note's three packets on the fine "
                         "lattice (~5 minutes each)")
    args = ap.parse_args()
    part_a()
    part_b()
    part_c()
    part_d()
    part_e(full=args.full)


if __name__ == "__main__":
    main()
