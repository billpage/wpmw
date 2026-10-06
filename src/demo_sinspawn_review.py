#!/usr/bin/env python3
r"""
Companion demo to docs/supplement/sinspawn_v1_review.md: D. Cyganski's
WignerParticlesSinSpawnV1 algorithm (quadratic potential as exact curved
trajectories, a sine as signed pair spawning), with the review's bugs fixed
or switched back on one at a time, against the exact Schroedinger solution.

Parts
-----
A  SymPy: the pair kernel Gamma(x)[W(k+q/2) - W(k-q/2)], Gamma = (Vp/hbar)cos qx,
   reproduces every term of the Moyal series for Vp sin(qx) (s = 0..5); its
   mean force is -Vp q cos qx; the exact rotation has determinant 1 and the
   notebook's view-aliased map has determinant cos^2(omega dt) (bug D1).
B  Numbers for D1 over the notebook's 9999 steps: area and radius factors,
   and the end point of a classical orbit started at x0 under both maps.
C  The split-operator reference: norm, and convergence under grid x2, dt/2.
D  Monte Carlo (wpmwlib/sinspawn.py), two strengths of the sine:
     "his"     nharm = 2,  Vp = 1e-4 hbar/dt (= FracSpawn hbar/dt, 0.658 meV)
     "matters" nharm = 16, Vp = 0.2 m omega^2 |x0| / q (24.9 meV)
   fixed algorithm at several N0, each of bugs A, B, C, D1 alone, and all
   bugs together with a population cap of 3 N0.
E  Classical baseline: Newtonian trajectories in the full potential, sampled
   from the same initial Wigner function (the truncated-Moyal answer).
F  Figures: sinspawn_review_marginals.png, sinspawn_review_scaling.png and,
   when the two notebook track plots are present in the output directory
   (sinspawn_tracks_original.png, sinspawn_tracks_fixed.png),
   sinspawn_review_tracks.png.

Writes sinspawn_review_results.json (every run, marginals, logs) through
output_path().  About half an hour on two cores; WPMW_QUICK=1 runs a smoke
test with small N0.

Run:  WPMW_OUTPUT=<dir> PYTHONPATH=src python3 -u src/demo_sinspawn_review.py
"""
import json
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from wpmwlib.sinspawn import (HBAR, M_E, Q_E, Setup, SinSpawn, marginal,  # noqa: E402
                              schroedinger, wigner_initial_sample)
from wpmwlib.wpmw_utils import docs_path, output_path  # noqa: E402

QUICK = os.environ.get("WPMW_QUICK", "") not in ("", "0")
NSTEPS = 9999                      # the notebook's loop: range(1, 10000)


# ----------------------------------------------------------------------------
# set-ups
# ----------------------------------------------------------------------------
def setup(case):
    if case == "his":              # FracSpawn = 1e-4 per step  <=>  Vp = 1e-4 hbar/dt
        return Setup(nharm=2.0, Vp=1e-4 * HBAR / 1e-16)
    S = Setup(nharm=16.0)          # "matters": sine force amplitude = 0.2 x the
    S.Vp = 0.2 * M_E * S.omega ** 2 * abs(S.x0) / S.q   # harmonic force at x0
    return S


def edges_for(S):
    lo, hi = -0.1e-6, 0.1e-6
    i0, i1 = int(np.floor(lo / S.dx)), int(np.ceil(hi / S.dx))
    return (np.arange(i0, i1 + 1) - 0.5) * S.dx


def binned(xs, dxs, r, e):
    xe = np.concatenate([xs - dxs / 2, [xs[-1] + dxs / 2]])
    cum = np.concatenate([[0.0], np.cumsum(r) * dxs])
    return np.diff(np.interp(e, xe, cum)) / np.diff(e)


def reference(case, times):
    S = setup(case)
    e = edges_for(S)
    xs, dxs, R = schroedinger(S, times)
    out = {t: binned(xs, dxs, r, e) for t, r in R.items()}
    _, _, R0 = schroedinger(Setup(nharm=S.nharm, Vp=0.0), times)
    out0 = {t: binned(xs, dxs, r, e) for t, r in R0.items()}
    return out, out0


# ----------------------------------------------------------------------------
# Part A -- exact checks
# ----------------------------------------------------------------------------
def part_a():
    import sympy as sp
    print("=== Part A: exact checks (SymPy) ===")
    x, k, q, Vp, hb, y = sp.symbols("x k q V_p hbar y", real=True)
    U = Vp * sp.sin(q * x)
    Gam = Vp / hb * sp.cos(q * x)
    # Moyal series in k = p/hbar:  sum_s (-1)^s U^(2s+1) / (hbar 4^s (2s+1)!) d_k^(2s+1) W
    # pair kernel, Taylor expanded:  Gam * sum_s 2 (q/2)^(2s+1)/(2s+1)! d_k^(2s+1) W
    for s in range(6):
        moyal = (-1) ** s * sp.diff(U, x, 2 * s + 1) / (hb * 4 ** s * sp.factorial(2 * s + 1))
        kern = Gam * 2 * (q / 2) ** (2 * s + 1) / sp.factorial(2 * s + 1)
        print(f"  order d_k^{2*s+1}:  Moyal - kernel = {sp.simplify(moyal - kern)}")
    # mean force: d<hbar k>/dt from the kernel, per unit weight at k0
    k0 = sp.symbols("k_0", real=True)
    F = hb * Gam * ((k0 - q / 2) - (k0 + q / 2))
    print("  mean force per unit weight, hbar*Gamma*[(k0-q/2)-(k0+q/2)] + U'(x) =",
          sp.simplify(F + sp.diff(U, x)))
    c, sn, s_, X, K = sp.symbols("c sigma s X K", real=True)
    xn = c * X + HB_sym(sp) * sn / s_ * K
    exact = sp.Matrix([xn, c * K - s_ * sn / HB_sym(sp) * X])
    viewed = sp.Matrix([xn, c * K - s_ * sn / HB_sym(sp) * xn])
    J1 = exact.jacobian([X, K]).det()
    J2 = viewed.jacobian([X, K]).det()
    print("  det exact map  =", sp.simplify(J1.subs(sn ** 2, 1 - c ** 2)))
    print("  det view map   =", sp.factor(sp.simplify(J2.subs(sn ** 2, 1 - c ** 2))),
          "  (bug D1)")


def HB_sym(sp):
    return sp.Symbol("hbar", positive=True)


# ----------------------------------------------------------------------------
# Part B -- the contracting map, in numbers
# ----------------------------------------------------------------------------
def part_b():
    print("\n=== Part B: bug D1 over the notebook's run ===")
    S = setup("his")
    c, sn = np.cos(S.omega * S.dt), np.sin(S.omega * S.dt)
    x, k = S.x0, 0.0
    for _ in range(NSTEPS):
        xn = k * HBAR * sn / S.s + x * c
        k = k * c - xn * S.s * sn / HBAR
        x = xn
    print(f"  omega dt = {S.omega*S.dt:.6f};  area factor c^(2*9999) = {c**(2*NSTEPS):.4f};"
          f"  radius factor c^9999 = {c**NSTEPS:.4f}")
    print(f"  orbit from x0 = {S.x0*1e9:.1f} nm after 9999 steps: view map {x*1e9:.3f} nm,"
          f"  exact {S.x0*np.cos(S.omega*NSTEPS*S.dt)*1e9:.3f} nm")


# ----------------------------------------------------------------------------
# Part C -- the reference
# ----------------------------------------------------------------------------
def part_c():
    print("\n=== Part C: split-operator reference ===")
    S = setup("matters")
    T = NSTEPS * S.dt
    x1, d1, R1 = schroedinger(S, [T])
    if QUICK:
        print(f"  norm at T: {R1[T].sum()*d1:.12f}  (refinement check skipped in quick mode)")
        return
    _, _, R2 = schroedinger(S, [T], N=8192, dt=S.dt / 2)
    print(f"  norm at T: {R1[T].sum()*d1:.12f};  L1(grid x2, dt/2) at T = "
          f"{np.sum(np.abs(R1[T] - R2[T][::2]))*d1:.2e}")


# ----------------------------------------------------------------------------
# Part D -- Monte Carlo
# ----------------------------------------------------------------------------
def one(job):
    case, label, flags, N0, seed = job
    S = setup(case)
    times = [NSTEPS * S.dt / 4, NSTEPS * S.dt / 2, NSTEPS * S.dt]
    times = [round(t / S.dt) * S.dt for t in times]
    rng = np.random.default_rng(seed)
    x, k, s = wigner_initial_sample(S, N0, rng)
    cap = 3 * N0 if flags.get("cap") else None
    f = {kk: v for kk, v in flags.items() if kk != "cap"}
    eng = SinSpawn(S, rng, cap=cap, **f)
    t0 = time.time()
    out, log = eng.run(x, k, s, times, log_every=100)
    e = edges_for(S)
    marg = {t: marginal(S, *out[t][0:1], out[t][2], N0, e) for t in times}
    return dict(case=case, label=label, flags=flags, N0=N0, seed=seed,
                seconds=time.time() - t0, births=eng.births, cancelled=eng.cancelled,
                dropped=eng.dropped, peak_pop=int(log[:, 1].max()), final_pop=int(log[-1, 1]),
                net=int(out[times[-1]][2].sum()), times=times,
                marg={str(t): m.tolist() for t, m in marg.items()},
                log=log.tolist())


def classical(case, rf, N=400000, seed=7):
    """Part E: Newtonian trajectories (velocity Verlet) in the full potential."""
    S = setup(case)
    e = np.array(rf["edges"]); h = np.diff(e)
    rng = np.random.default_rng(seed)
    x = rng.normal(S.x0, S.a, N); p = HBAR * rng.normal(S.k0, 1 / (2 * S.a), N)

    def F(x):
        return -(2 * Q_E * S.V2 * x + S.Vp * S.q * np.cos(S.q * x))
    want = {int(round(t / S.dt)): t for t in rf["times"]}
    dt, out = S.dt, {}
    for i in range(1, max(want) + 1):
        p += 0.5 * dt * F(x); x += dt * p / M_E; p += 0.5 * dt * F(x)
        if i in want:
            t = want[i]
            m = np.histogram(x, bins=e)[0] / N / h
            out[str(t)] = dict(L1=float(np.sum(np.abs(m - np.array(rf["ref"][str(t)])) * h)),
                               marg=m.tolist())
    return out


def mean_x(rf, m):
    e = np.array(rf["edges"]); xc = 0.5 * (e[1:] + e[:-1])
    return float(np.sum(xc * np.array(m) * np.diff(e)))


def part_d_e():
    print("\n=== Part D: signed particles against the exact solution ===")
    BUGS_ALL = dict(rectify=True, wrong_array=True, flip_children=True, view_map=True,
                    no_cancel=True, cap=True)
    big = [2000, 8000] if QUICK else [20000, 80000, 320000]
    mid = 2000 if QUICK else 80000
    jobs = [("matters", "fixed", {}, N0, 11) for N0 in big]
    for name in ("rectify", "wrong_array", "flip_children", "view_map"):
        jobs.append(("matters", "bug " + name, {name: True}, mid, 12))
    jobs.append(("matters", "notebook (all bugs, 3x cap)", BUGS_ALL, mid, 13))
    jobs.append(("his", "fixed", {}, mid, 14))
    jobs.append(("his", "notebook (all bugs, 3x cap)", BUGS_ALL, mid, 15))
    jobs.append(("his", "fixed, N0 = 838", {}, 838, 16))
    jobs.sort(key=lambda j: -j[3])
    with Pool(2) as pool:
        res = pool.map(one, jobs, chunksize=1)
    refs = {}
    for case in ("matters", "his"):
        S = setup(case)
        times = res[[r["case"] for r in res].index(case)]["times"]
        R, R0 = reference(case, times)
        refs[case] = dict(times=times, ref={str(t): R[t].tolist() for t in times},
                          qho={str(t): R0[t].tolist() for t in times},
                          edges=edges_for(S).tolist(), Vp=S.Vp, nharm=S.nharm, q=S.q,
                          dk=S.dk, dx=S.dx, omega=S.omega)
    for r in res:
        rf = refs[r["case"]]
        h = np.diff(rf["edges"])
        r["L1"] = {t: float(np.sum(np.abs(np.array(r["marg"][t]) - np.array(rf["ref"][t])) * h))
                   for t in r["marg"]}
        r["mean_x"] = mean_x(rf, r["marg"][str(r["times"][-1])])
    for case, rf in refs.items():
        h = np.diff(rf["edges"])
        rf["L1_ref_vs_qho"] = {t: float(np.sum(np.abs(np.array(rf["ref"][t]) - np.array(rf["qho"][t])) * h))
                               for t in rf["ref"]}
        T = str(rf["times"][-1])
        p = np.clip(np.array(rf["ref"][T]) * h, 0, None); p /= p.sum()
        rf["floor"] = {str(N): float(np.sum(np.sqrt(2 * p * (1 - p) / (np.pi * N))))
                       for N in (838, 20000, 80000, 320000)}
        S = setup(case)
        print(f"{case}: Vp = {S.Vp/Q_E*1e3:.3f} meV, Gamma_max dt = {S.Vp/HBAR*S.dt:.4f}, "
              f"kick q/2 = {S.q/2/S.dk:.1f} dk = {S.q/2*2*S.a:.3f} sigma_k, (q a)^2 = {(S.q*S.a)**2:.4f}, "
              f"sine force / harmonic force at x0 = {S.Vp*S.q/(M_E*S.omega**2*abs(S.x0)):.2e}")
        print(f"  exact <x>(T) = {mean_x(rf, rf['ref'][T])*1e9:.3f} nm;  sine off: "
              f"<x>(T) = {mean_x(rf, rf['qho'][T])*1e9:.3f} nm, L1 vs exact (t/4, t/2, T) = "
              + " ".join(f"{rf['L1_ref_vs_qho'][str(t)]:.3f}" for t in rf["times"]))
        print("  ideal positive-sampler L1 floor at T: "
              + ", ".join(f"N0={N}: {v:.3f}" for N, v in rf["floor"].items()))
        for r in sorted([r for r in res if r["case"] == case], key=lambda r: (r["label"], r["N0"])):
            print(f"  {r['label']:30s} N0={r['N0']:7d}  L1(t/4,t/2,T)="
                  + " ".join(f"{r['L1'][str(t)]:.3f}" for t in r["times"])
                  + f"  <x>(T)={r['mean_x']*1e9:7.3f} nm  peak/N0={r['peak_pop']/r['N0']:.2f}"
                  + f"  births/N0={r['births']/r['N0']:.1f}  cancelled/births="
                  + f"{(r['cancelled']/r['births'] if r['births'] else 0):.4f}"
                  + f"  dropped/N0={r['dropped']/r['N0']:.1f}  ({r['seconds']:.0f}s)")

    print("\n=== Part E: classical trajectories in the full potential ===")
    cl = {}
    for case, rf in refs.items():
        cl[case] = classical(case, rf, N=40000 if QUICK else 400000)
        for t in rf["times"]:
            c = cl[case][str(t)]
            print(f"  {case:8s} t = {t*1e12:.3f} ps  classical vs exact L1 = {c['L1']:.3f}  "
                  f"<x> classical {mean_x(rf, c['marg'])*1e9:.3f} nm, exact "
                  f"{mean_x(rf, rf['ref'][str(t)])*1e9:.3f} nm")
    D = dict(results=res, refs=refs, classical=cl)
    with open(output_path("sinspawn_review_results.json"), "w") as fh:
        json.dump(D, fh)
    return D


# ----------------------------------------------------------------------------
# Part F -- figures
# ----------------------------------------------------------------------------
def save(fig, name):
    fig.savefig(output_path(name), dpi=120, bbox_inches="tight")
    dp = docs_path(name)
    if dp:
        fig.savefig(dp, dpi=120, bbox_inches="tight")


def part_f(D):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.image as mpimg
    import matplotlib.pyplot as plt
    print("\n=== Part F: figures ===")
    R, refs, CL = D["results"], D["refs"], D["classical"]

    def pick(case, label):
        c = [r for r in R if r["case"] == case and r["label"] == label]
        return max(c, key=lambda r: r["N0"]) if c else None

    rf = refs["matters"]
    e = np.array(rf["edges"]); xc = 0.5 * (e[1:] + e[:-1]) * 1e9
    fig, ax = plt.subplots(2, 2, figsize=(14, 8.5), sharex=True)
    for col, t in enumerate(rf["times"][1:]):
        T = str(t)
        for row, labels in enumerate([["fixed", "notebook (all bugs, 3x cap)"],
                                      ["bug rectify", "bug wrong_array", "bug flip_children", "bug view_map"]]):
            a = ax[row, col]
            a.plot(xc, rf["ref"][T], "k", lw=2.2, label="exact (Schrödinger)")
            a.plot(xc, rf["qho"][T], color="0.6", lw=1.2, ls="--", label="exact, sine switched off")
            if row == 0:
                a.plot(xc, CL["matters"][T]["marg"], color="m", lw=1.0, ls=":",
                       label=f"classical Liouville, full force (L1={CL['matters'][T]['L1']:.2f})")
            for lab, colr in zip(labels, ["C0", "C3", "C1", "C2"]):
                r = pick("matters", lab)
                if r is None:
                    continue
                a.plot(xc, r["marg"][T], color=colr, lw=1.1,
                       label=f"{lab}  (N0={r['N0']:,}, L1={r['L1'][T]:.2f})")
            a.set_xlim(-45, 45)
            a.set_title(f"t = {t*1e12:.2f} ps" + ("" if row == 0 else "  — one bug at a time"))
            a.legend(fontsize=7.5, loc="upper left")
            if row == 1:
                a.set_xlabel("x (nm)")
            if col == 0:
                a.set_ylabel(r"$\rho(x)$  (1/m)")
    fig.suptitle(f"V = q_e V2 x² + Vp sin(2π·{rf['nharm']:.0f}x/LX), Vp = {rf['Vp']/Q_E*1e3:.1f} meV: "
                 "position marginal, signed particles against the exact solution", fontsize=11)
    fig.tight_layout()
    save(fig, "sinspawn_review_marginals.png")
    plt.close(fig)

    fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))
    fx = sorted([r for r in R if r["case"] == "matters" and r["label"] == "fixed"], key=lambda r: r["N0"])
    N = np.array([r["N0"] for r in fx])
    for t, mk in zip(rf["times"], ("o", "s", "^")):
        L = np.array([r["L1"][str(t)] for r in fx])
        ax[0].loglog(N, L, mk + "-", label=f"t = {t*1e12:.2f} ps")
    ax[0].loglog(N, L[0] * np.sqrt(N[0] / N), "k:", lw=1, label=r"$\propto N_0^{-1/2}$")
    for t in rf["times"]:
        ax[0].axhline(rf["L1_ref_vs_qho"][str(t)], color="0.7", lw=0.8)
    ax[0].loglog(N, [CL["matters"][str(rf["times"][-1])]["L1"]] * len(N), "m:", lw=1.2,
                 label="classical Liouville vs exact (t = 1 ps)")
    p = np.clip(np.array(rf["ref"][str(rf["times"][-1])]) * np.diff(rf["edges"]), 0, None); p /= p.sum()
    Nf = np.geomspace(N[0], N[-1], 20)
    ax[0].loglog(Nf, [np.sum(np.sqrt(2 * p * (1 - p) / (np.pi * n))) for n in Nf], color="C2", ls="--",
                 lw=0.8, label="sampling floor, exact positive sampler (t = 1 ps)")
    ax[0].set_xlabel("$N_0$"); ax[0].set_ylabel(r"$L^1$ error of $\rho(x)$")
    ax[0].set_title("(a) fixed algorithm (grey: sine on vs off)"); ax[0].legend(fontsize=8)
    for r in fx:
        lg = np.array(r["log"])
        ax[1].plot(lg[:, 0] * 1e12, lg[:, 1] / r["N0"], label=f"N0={r['N0']:,}")
    r = pick("matters", "notebook (all bugs, 3x cap)")
    lg = np.array(r["log"])
    ax[1].plot(lg[:, 0] * 1e12, lg[:, 1] / r["N0"], "C3--", label="notebook, all bugs (capped)")
    ax[1].set_xlabel("t (ps)"); ax[1].set_ylabel("particles / $N_0$"); ax[1].set_title("(b) population")
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    save(fig, "sinspawn_review_scaling.png")
    plt.close(fig)
    names = ["sinspawn_review_marginals.png", "sinspawn_review_scaling.png"]

    srcs = [output_path(f) for f in ("sinspawn_tracks_original.png", "sinspawn_tracks_fixed.png")]
    if all(os.path.exists(f) for f in srcs):
        fig, ax = plt.subplots(1, 2, figsize=(14, 4.8))
        for a, f, ttl in zip(ax, srcs, ("original notebook (saved output)", "fixed copy, executed")):
            a.imshow(mpimg.imread(f)); a.axis("off"); a.set_title(ttl)
        fig.suptitle("Tracked particles at the notebook's own parameters (blue positons, red negatons)")
        fig.tight_layout()
        save(fig, "sinspawn_review_tracks.png")
        plt.close(fig)
        names.append("sinspawn_review_tracks.png")
    else:
        print("  track plots not found in the output directory; sinspawn_review_tracks.png skipped")
    print("  wrote " + ", ".join(names))


if __name__ == "__main__":
    t0 = time.time()
    part_a()
    part_b()
    part_c()
    D = part_d_e()
    part_f(D)
    print(f"\ntotal {time.time()-t0:.0f} s" + ("  (quick mode)" if QUICK else ""))
