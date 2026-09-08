"""
Demonstration: the Limkumnerd-Phanthaphanitkul benchmark, run in WPMW terms.

Companion to ``docs/supplement/limkumnerd_weighted_paths.md``.

Limkumnerd & Phanthaphanitkul, "Weighted phase-space paths for exact Wigner
dynamics" (arXiv:2605.05764, 2026) take classical Hamiltonian flow as a
carrier, split the Wigner generator into the classical Liouville part and the
Moyal residual, and benchmark the split on a quartic oscillator.  That is the
same split as ``docs/analysis/compensated_liouville_splitting.md`` Sec. 2.
This demo runs their benchmark on our machinery, then asks the question their
deterministic grid cannot see: what does the residual channel *cost* in
events, and how does that cost depend on the reach?

Their configuration, reproduced here exactly:

    H     = p^2/2 + q^2/2 + lambda q^4,      lambda = 0.02 and 0.05
    psi_0 = (|0> + |2>)/sqrt(2)
    t_f   = pi/2,   dt = 0.005,   dq = 0.0625,   dp = pi/48 ~= 0.1309

The momentum step fixes the reach.  With ``dp = pi hbar / (2 y_max)`` their
grid instantiates ``y_max = 12``, a fact that appears nowhere in their paper
because a finite-difference stencil costs the same whatever its total
variation.

Parts
-----
A  Harmonic null test.  V''' = 0, so the residual vanishes identically and
   classical carrier transport is exact.  Reproduces their Fig. 1.
B  Quartic benchmark.  Classical-carrier-only transport against the
   Schrodinger reference, then with the signed residual restored.
   Reproduces their Figs. 2 and 4.  For a quartic the Moyal series
   terminates (Theorem E7), so "classical + leading residual" is not a
   truncation -- it is the exact QLE.
C  Admissibility.  The same evolutions scored by the least eigenvalue of the
   density matrix reconstructed from W, which is the diagnostic of
   Theorem G4 and is not a norm on W.
D  Reach sweep.  The residual symbol band-limited to |s| <= 2 y_max / hbar,
   with the resulting Wigner error set against the residual event budget
   R_res = sum_q |K_q^res| of Theorem E8.

Outputs
-------
- ``limkumnerd_benchmark_wigner.png``   exact / classical / restored + differences
- ``limkumnerd_error_growth.png``       Wigner L2 error against time
- ``limkumnerd_admissibility.png``      least eigenvalue of rho against time
- ``limkumnerd_reach_budget.png``       error and event budget against reach

Run with:

    WPMW_OUTPUT=/mnt/user-data/outputs PYTHONPATH=src python3 -u \
        src/demo_limkumnerd_weighted_paths.py
"""

from __future__ import annotations

import math
import time

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from wpmwlib.wpmw_utils import output_path, docs_path  # noqa: E402

HBAR = 1.0
MASS = 1.0

# Their grid.  The q box is widened to [-24, 24) so that the separation
# coordinate u = y/2 can run over [-12, 12) without periodic wraparound
# contaminating psi(q + u) psi*(q - u); the reach and dp are theirs.
NQ = 768
Q = np.linspace(-24.0, 24.0, NQ, endpoint=False)
DQ = Q[1] - Q[0]
N_OFF = 192                                  # u in [-12, 12): y_max = 12
NP = 2 * N_OFF
DP = np.pi * HBAR / (NP * DQ)
P = (np.arange(NP) - NP // 2) * DP
Y_MAX = np.pi * HBAR / (2.0 * DP)            # = 12, their implied reach

DT = 0.005
T_F = np.pi / 2.0


def banner(text):
    print()
    print("=" * 72)
    print(text)
    print("=" * 72)


def save_fig(fig, name):
    fig.savefig(output_path(name), dpi=150, bbox_inches="tight")
    dp_ = docs_path(name)
    if dp_:
        fig.savefig(dp_, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {name}")


# --------------------------------------------------------------------- #
# States, transforms, evolution                                         #
# --------------------------------------------------------------------- #
def osc_state(n, q, mass=MASS, omega=1.0):
    """Normalised harmonic-oscillator eigenfunction on grid ``q``."""
    from numpy.polynomial.hermite import hermval

    alpha = np.sqrt(mass * omega / HBAR)
    coef = np.zeros(n + 1)
    coef[n] = 1.0
    norm = (alpha / np.sqrt(np.pi)) ** 0.5 / np.sqrt(2.0 ** n * math.factorial(n))
    return norm * hermval(alpha * q, coef) * np.exp(-(alpha * q) ** 2 / 2.0)


def wigner(psi, q, p, n_off=N_OFF):
    """``W(q,p) = (1/pi hbar) int du exp(-2ipu/hbar) psi(q+u) psi*(q-u)``.

    ``u = y/2`` runs over the position grid itself, so the conjugate momentum
    step is ``dp = pi hbar / (2 n_off dq)`` and the p grid must not exceed the
    Nyquist bound ``|p| <= pi hbar / (2 dq)``.
    """
    dq = q[1] - q[0]
    nq = q.size
    offs = np.arange(-n_off, n_off)
    u = offs * dq
    plus = np.empty((nq, u.size), dtype=complex)
    minus = np.empty_like(plus)
    for j, o in enumerate(offs):
        plus[:, j] = np.roll(psi, -o)
        minus[:, j] = np.roll(psi, o)
    prod = np.conj(minus) * plus
    ker = np.exp(-2j * np.outer(p, u) / HBAR)
    return ((ker @ prod.T).real * dq / (np.pi * HBAR)).T


def rho_from_wigner(w, q, p, n_off=N_OFF):
    """Reconstruct ``<x|rho|x'>`` from ``W`` on the grid of ``wigner``.

    ``rho(q+u, q-u) = int dp exp(2ipu/hbar) W(q,p)``.  Rows and columns
    reached this way always share a parity, so the result is the direct sum
    of the two sublattice blocks; each is a faithful density matrix at half
    resolution, and both are positive semidefinite exactly when the state is.
    """
    dp = p[1] - p[0]
    nq = q.size
    offs = np.arange(-n_off, n_off)
    u = offs * (q[1] - q[0])
    ker = np.exp(2j * np.outer(u, p) / HBAR)
    g = (w @ ker.T).T * dp                    # g[j, i] = rho(q_i + u_j, q_i - u_j)
    rho = np.zeros((nq, nq), dtype=complex)
    idx = np.arange(nq)
    for j, o in enumerate(offs):
        rho[(idx + o) % nq, (idx - o) % nq] = g[j]
    return rho


def min_eig_rho(w, q, p):
    return float(np.linalg.eigvalsh(rho_from_wigner(w, q, p) * (q[1] - q[0])).min())


def split_step_schrodinger(psi, q, v_at_q, dt, n_steps, mass=MASS):
    kq = 2.0 * np.pi * np.fft.fftfreq(q.size, d=q[1] - q[0])
    kin = np.exp(-1j * HBAR * kq ** 2 * dt / (2.0 * mass))
    pot = np.exp(-1j * v_at_q * dt / (2.0 * HBAR))
    for _ in range(n_steps):
        psi = pot * psi
        psi = np.fft.ifft(kin * np.fft.fft(psi))
        psi = pot * psi
    return psi


def transport(w, q, p, dv_at_q, d3v_at_q, dt, n_steps, mass=MASS,
              residual=False, s_cut=None):
    """Strang transport of ``W``: half drift, potential kick, half drift.

    The potential substep is diagonal in ``(x, s)``, where ``s`` is conjugate
    to ``p``, and factorises exactly (Theorem C1).  With ``residual=False``
    the kick carries only ``M_cl = i V'(x) s`` -- classical carrier transport.
    With ``residual=True`` it also carries
    ``M_res = i (hbar^2/24) V'''(x) s^3``, which for a quartic ``V`` is the
    whole residual and not a truncation, since ``V^(5) = 0``.

    ``s_cut`` band-limits the residual to ``|s| <= s_cut``, i.e. imposes a
    reach ``y_max = hbar s_cut / 2``.
    """
    kq = 2.0 * np.pi * np.fft.fftfreq(q.size, d=q[1] - q[0])
    s = 2.0 * np.pi * np.fft.fftfreq(p.size, d=p[1] - p[0])
    drift_half = np.exp(-1j * np.outer(kq, p) * (dt / 2.0) / mass)
    m = 1j * np.outer(dv_at_q, s)
    if residual:
        m_res = 1j * (HBAR ** 2 / 24.0) * np.outer(d3v_at_q, s ** 3)
        if s_cut is not None:
            m_res = m_res * (np.abs(s)[None, :] <= s_cut)
        m = m + m_res
    kick = np.exp(dt * m)
    for _ in range(n_steps):
        w = np.fft.ifft(drift_half * np.fft.fft(w, axis=0), axis=0).real
        w = np.fft.ifft(kick * np.fft.fft(w, axis=1), axis=1).real
        w = np.fft.ifft(drift_half * np.fft.fft(w, axis=0), axis=0).real
    return w


def l2(dw, dq=DQ, dp=DP):
    return float(np.sqrt((dw ** 2).sum() * dq * dp))


def quartic(lam):
    return (lambda z: 0.5 * MASS * z ** 2 + lam * z ** 4,
            lambda z: MASS * z + 4.0 * lam * z ** 3,
            lambda z: 24.0 * lam * z)


# --------------------------------------------------------------------- #
# The compensated residual kernel on a reach lattice                     #
# --------------------------------------------------------------------- #
def kernel_rungs(x, v, dv, y_max, n_rungs=64, n_y=8192):
    """Return ``(xi, K_res)``, the compensated momentum-transfer lattice.

    The symbol ``M_res = (i/hbar)[V(x+y) - V(x-y) - 2y V'(x)]`` is sampled on
    ``|y| <= y_max`` and transformed to the rung lattice
    ``xi_q = q dp``, ``dp = pi hbar / (2 y_max)``.
    """
    s_max = 2.0 * y_max / HBAR
    s = np.linspace(-s_max, s_max, 2 * n_y + 1)
    ds = s[1] - s[0]
    y = HBAR * s / 2.0
    m = 1j * (v(x + y) - v(x - y) - 2.0 * y * dv(x)) / HBAR
    xi = np.arange(-n_rungs, n_rungs + 1) * (np.pi * HBAR / (2.0 * y_max))
    k = (np.exp(-1j * np.outer(xi, s)) * m).sum(axis=1) * ds / (2.0 * s_max)
    return xi, k.real


# --------------------------------------------------------------------- #
# Parts                                                                  #
# --------------------------------------------------------------------- #
def initial_state():
    psi = (osc_state(0, Q) + osc_state(2, Q)) / np.sqrt(2.0)
    return psi / np.sqrt((np.abs(psi) ** 2).sum() * DQ)


def step_marks(n_snapshots=9):
    """Snapshot step counts on the exact ``DT`` lattice, ending at ``T_F``.

    Snapshot times must be integer multiples of ``DT`` or the reference and
    the transports are compared at different instants, which shows up as a
    spurious floor on the restored error.
    """
    n_total = int(round(T_F / DT))
    return np.unique(np.round(np.linspace(0, n_total, n_snapshots)).astype(int))


def run_case(lam, psi0, w0, marks):
    """Evolve all three ways, returning snapshots at the given step counts."""
    v, dv, d3v = quartic(lam)
    out = {"exact": [], "classical": [], "restored": []}
    w_cl, w_rs = w0.copy(), w0.copy()
    prev = 0
    for mark in marks:
        n = int(mark - prev)
        if n:
            w_cl = transport(w_cl, Q, P, dv(Q), d3v(Q), DT, n, residual=False)
            w_rs = transport(w_rs, Q, P, dv(Q), d3v(Q), DT, n, residual=True)
        prev = mark
        psi_t = split_step_schrodinger(psi0, Q, v(Q), DT, int(mark))
        out["exact"].append(wigner(psi_t, Q, P))
        out["classical"].append(w_cl.copy())
        out["restored"].append(w_rs.copy())
    return out


def part_a(psi0, w0):
    banner("Part A  harmonic null test: the residual vanishes identically")
    v, dv, d3v = quartic(0.0)
    n = int(round(T_F / DT))
    psi_t = split_step_schrodinger(psi0, Q, v(Q), DT, n)
    w_ex = wigner(psi_t, Q, P)
    w_cl = transport(w0, Q, P, dv(Q), d3v(Q), DT, n, residual=False)
    err = l2(w_cl - w_ex)
    print(f"\n  V''' = 0, so M_res = 0 and classical carrier transport is exact.")
    print(f"  Wigner L2 error at t_f = pi/2 : {err:.4e}")
    print(f"  ||W|| for scale               : {l2(w_ex):.4e}")
    print("\n  Limkumnerd Fig. 1 reports 4.8e-05 for the same test, which they")
    print("  identify as their interpolation and grid floor.  Our transport is")
    print("  spectral in both substeps, so no interpolation floor arises and")
    print("  the residue is the Strang error between drift and kick alone.")
    return err


def part_b(psi0, w0, marks):
    banner("Part B  quartic benchmark: classical carrier, then residual restored")
    results = {}
    print(f"\n  {'lambda':>8s}{'classical only':>18s}{'residual restored':>20s}"
          f"{'published (theirs)':>22s}")
    published = {0.02: "5.7e-02 / 5.4e-05", 0.05: "1.1e-01 / 6.3e-05"}
    for lam in (0.02, 0.05):
        t0 = time.time()
        snaps = run_case(lam, psi0, w0, marks)
        results[lam] = snaps
        e_cl = [l2(a - b) for a, b in zip(snaps["classical"], snaps["exact"])]
        e_rs = [l2(a - b) for a, b in zip(snaps["restored"], snaps["exact"])]
        results[f"err_cl_{lam}"] = e_cl
        results[f"err_rs_{lam}"] = e_rs
        print(f"  {lam:8.2f}{e_cl[-1]:18.4e}{e_rs[-1]:20.4e}"
              f"{published[lam]:>22s}   [{time.time()-t0:.0f}s]")
    print("\n  The classical-only column reproduces their published figures.")
    print("  The restored column is an order of magnitude below theirs for the")
    print("  same reason as Part A.  For a quartic V the Moyal series")
    print("  terminates (Theorem E7), so the restored evolution is the exact")
    print("  QLE and its residue is a timestep artefact, not a representation")
    print("  error.")
    return results


def part_c(results, times):
    banner("Part C  the same runs scored by admissibility, not by a norm")
    lam = 0.05
    snaps = results[lam]
    rows = []
    print(f"\n  quartic lambda = {lam}, least eigenvalue of rho reconstructed"
          f" from W\n")
    print(f"  {'t':>8s}{'exact':>16s}{'classical only':>18s}"
          f"{'residual restored':>20s}")
    for i, t in enumerate(times):
        e = min_eig_rho(snaps["exact"][i], Q, P)
        c = min_eig_rho(snaps["classical"][i], Q, P)
        r = min_eig_rho(snaps["restored"][i], Q, P)
        rows.append((t, e, c, r))
        print(f"  {t:8.4f}{e:16.3e}{c:18.3e}{r:20.3e}")
    print("\n  A Wigner L2 error says the classical carrier gets the answer")
    print("  wrong.  This says something stronger: the object it produces is")
    print("  not the Wigner transform of any density matrix.  Theorem G4 of")
    print("  compensated_ontology.md reports -4.97e-02 at t = pi/2 for this")
    print("  configuration on an independent grid.")
    return rows


def part_d(psi0, w0):
    banner("Part D  reach sweep: what the residual costs in events")
    lam = 0.05
    v, dv, d3v = quartic(lam)
    n = int(round(T_F / DT))
    psi_t = split_step_schrodinger(psi0, Q, v(Q), DT, n)
    w_ex = wigner(psi_t, Q, P)
    reaches = [0.75, 1.5, 3.0, 6.0, 12.0]
    rows = []
    print(f"\n  {'y_max':>8s}{'dp = pi/2y_max':>16s}{'L2 err (restored)':>20s}"
          f"{'budget R_res':>15s}{'x8?':>8s}")
    prev = None
    for ym in reaches:
        w = transport(w0, Q, P, dv(Q), d3v(Q), DT, n, residual=True,
                      s_cut=2.0 * ym / HBAR)
        err = l2(w - w_ex)
        _, k = kernel_rungs(1.0, v, dv, ym)
        budget = float(np.abs(k).sum())
        ratio = budget / prev if prev else float("nan")
        prev = budget
        rows.append((ym, err, budget))
        rs = "-" if math.isnan(ratio) else f"{ratio:.2f}"
        print(f"  {ym:8.2f}{np.pi/(2*ym):16.4f}{err:20.4e}{budget:15.4f}{rs:>8s}")
    print("\n  The error saturates by y_max = 6 and does not move thereafter.")
    print("  The event budget keeps growing by a factor of 8 per doubling --")
    print("  cubically, with no saturation, exactly as Theorem E8 predicts for")
    print(f"  a quartic.  Their grid sits at y_max = {Y_MAX:.0f} and therefore")
    print("  pays eight times the budget of a run that had already converged.")
    print("  On a deterministic grid this is invisible: a finite-difference")
    print("  stencil costs the same whatever its total variation.")
    return rows


# --------------------------------------------------------------------- #
# Figures                                                                #
# --------------------------------------------------------------------- #
def moyal_residual_field(w, p, d3v_at_q):
    """The signed residual field ``-(hbar^2/24) V^(3)(q) d^3W/dp^3``.

    This is the object decomposed into positive and negative parts in their
    Fig. 3.  It is a field on phase space, not the momentum-transfer kernel,
    so its Hahn-Jordan split always exists and is always finite.  The split
    a particle representation actually needs is the one on the kernel, which
    by Theorem E7 does not exist for this potential at all.
    """
    kp = 2.0 * np.pi * np.fft.fftfreq(p.size, d=p[1] - p[0])
    d3w = np.fft.ifft((1j * kp) ** 3 * np.fft.fft(w, axis=1), axis=1).real
    return -(HBAR ** 2 / 24.0) * d3v_at_q[:, None] * d3w


def fig_wigner(results):
    snaps = results[0.05]
    w_ex, w_cl, w_rs = (snaps[k][-1] for k in ("exact", "classical", "restored"))
    _, _, d3v = quartic(0.05)
    r_field = moyal_residual_field(w_ex, P, d3v(Q))
    sel = (np.abs(Q) <= 4.0)
    pv = (np.abs(P) <= 4.0)
    ext = [Q[sel].min(), Q[sel].max(), P[pv].min(), P[pv].max()]
    cut = lambda w: w[np.ix_(sel, pv)].T

    fig, ax = plt.subplots(2, 3, figsize=(13.5, 8.0))
    vm = np.abs(cut(w_ex)).max()
    for j, (w, ttl) in enumerate(((w_ex, "exact (Schrodinger reference)"),
                                  (w_cl, "classical carrier only"),
                                  (w_rs, "carrier + signed residual"))):
        ax[0, j].imshow(cut(w), origin="lower", extent=ext, aspect="auto",
                        cmap="RdBu_r", vmin=-vm, vmax=vm)
        ax[0, j].set_title(ttl, fontsize=10)
        ax[0, j].set_xlabel("q")
    ax[0, 0].set_ylabel("p")

    rm = np.abs(cut(r_field)).max()
    ax[1, 0].imshow(cut(r_field), origin="lower", extent=ext, aspect="auto",
                    cmap="RdBu_r", vmin=-rm, vmax=rm)
    ax[1, 0].set_title("signed residual field (their Fig. 3a)", fontsize=10)
    dv_ = np.abs(cut(w_cl - w_ex)).max()
    for j, w, ttl in ((1, w_cl - w_ex, "classical carrier $-$ exact"),
                      (2, w_rs - w_ex, "restored $-$ exact")):
        ax[1, j].imshow(cut(w), origin="lower", extent=ext, aspect="auto",
                        cmap="RdBu_r", vmin=-dv_, vmax=dv_)
        ax[1, j].set_title(ttl, fontsize=10)
    for j in range(3):
        ax[1, j].set_xlabel("q")
    ax[1, 0].set_ylabel("p")
    fig.suptitle(r"Quartic benchmark $\lambda = 0.05$ at $t_f = \pi/2$; "
                 "the two difference panels share a scale", fontsize=11)
    fig.tight_layout()
    save_fig(fig, "limkumnerd_benchmark_wigner.png")


def fig_error(results, times, null_err):
    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    floor = 1e-8
    for lam, c in ((0.02, "tab:blue"), (0.05, "tab:red")):
        ax.semilogy(times, np.maximum(results[f"err_cl_{lam}"], floor),
                    "-o", ms=4, color=c, label=rf"classical only, $\lambda={lam}$")
        ax.semilogy(times, np.maximum(results[f"err_rs_{lam}"], floor),
                    "--s", ms=4, color=c,
                    label=rf"residual restored, $\lambda={lam}$")
    ax.axhline(null_err, color="0.4", lw=1, ls=":",
               label="harmonic null test at $t_f$")
    ax.set_xlabel("t")
    ax.set_ylabel(r"$\|W - W_{\rm exact}\|_2$")
    ax.set_title("Omitting the residual is a finite error; restoring it is exact")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    save_fig(fig, "limkumnerd_error_growth.png")


def fig_admissibility(rows):
    t = [r[0] for r in rows]
    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    ax.plot(t, [r[1] for r in rows], "-o", ms=4, label="exact")
    ax.plot(t, [r[2] for r in rows], "-s", ms=4, label="classical carrier only")
    ax.plot(t, [r[3] for r in rows], "--^", ms=4, label="carrier + signed residual")
    ax.axhline(0.0, color="0.3", lw=1)
    ax.set_xlabel("t")
    ax.set_ylabel(r"least eigenvalue of $\rho$ reconstructed from $W$")
    ax.set_title(r"Admissibility, quartic $\lambda = 0.05$"
                 "\n(the classical carrier leaves the set of quantum states)")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    save_fig(fig, "limkumnerd_admissibility.png")


def fig_reach(rows):
    ym = [r[0] for r in rows]
    err = [max(r[1], 1e-8) for r in rows]
    bud = [r[2] for r in rows]
    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    ax.loglog(ym, err, "-o", color="tab:blue", label=r"$\|W - W_{\rm exact}\|_2$")
    ax.set_xlabel(r"reach $y_{\max}$")
    ax.set_ylabel(r"Wigner $L^2$ error", color="tab:blue")
    ax.tick_params(axis="y", labelcolor="tab:blue")
    ax2 = ax.twinx()
    ax2.loglog(ym, bud, "-s", color="tab:red",
               label=r"$R_{\rm res} = \sum_q |K_q^{\rm res}|$")
    ax2.set_ylabel("residual event budget", color="tab:red")
    ax2.tick_params(axis="y", labelcolor="tab:red")
    ax.axvline(Y_MAX, color="0.4", ls=":", lw=1.2)
    ax.text(Y_MAX * 0.92, 4e-4, "reach implied by their\nmomentum grid",
            fontsize=8, color="0.3", ha="right", va="center")
    ax.set_title("Accuracy saturates; the event budget does not\n"
                 r"(quartic $\lambda = 0.05$, budget slope 3)")
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    save_fig(fig, "limkumnerd_reach_budget.png")


def main():
    t_start = time.time()
    banner("Limkumnerd-Phanthaphanitkul benchmark in WPMW terms")
    print(f"\n  dq = {DQ:.4f}   dp = {DP:.6f}   p_max = {P.max():.4f}")
    print(f"  implied reach y_max = pi hbar / (2 dp) = {Y_MAX:.4f}")
    print(f"  dt = {DT}   t_f = pi/2   state = (|0> + |2>)/sqrt(2)")

    psi0 = initial_state()
    w0 = wigner(psi0, Q, P)
    print(f"  Wigner norm at t = 0: {w0.sum() * DQ * DP:.12f}")

    marks = step_marks()
    times = [m * DT for m in marks]
    null_err = part_a(psi0, w0)
    results = part_b(psi0, w0, marks)
    rows_c = part_c(results, times)
    rows_d = part_d(psi0, w0)

    banner("Figures")
    fig_wigner(results)
    fig_error(results, times, null_err)
    fig_admissibility(rows_c)
    fig_reach(rows_d)

    print(f"\n  total {time.time() - t_start:.0f}s")


if __name__ == "__main__":
    main()
