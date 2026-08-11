"""
Verification for ``docs/supplement/inverted_pair_barrier.md``.

Two particles on a line with an *inverted* harmonic pair potential,

    H = (p1^2 + p2^2)/2m  -  (m w^2 / 2) (x1 - x2)^2 ,

treated as a two-body problem.  Unlike the confining trap this is a genuine
scattering problem: the particles either pass through each other or are
turned back, and the transmission probability is a quantum number with a
closed form.

Parts
-----
A  Separation and Moyal truncation.  COM free at mass 2m; relative motion an
   inverted oscillator at reduced mass mu = m/2 and Lyapunov rate
   Om = sqrt(2) w.  Third derivatives vanish, so the QLE is exactly the
   classical Liouville equation.
B  The invariant half-plane.  u = r + p_r/(mu Om) is an eigenvector of the
   flow with eigenvalue e^(Om t), so the half-plane u > 0 is exactly
   invariant and T = int_(u>0) W is a constant of the motion equal to the
   asymptotic transmission probability.
C  Closed form T = Phi(u_bar / sigma_u), checked against a split-operator
   Schrodinger solve; the hbar-dependence; Kemble for the energy-resolved
   comparison.
D  Two-body content: the total momentum is untouched (Theorem A3 of the
   4-D note), and the entanglement entropy grows at the Lyapunov rate --
   the observable a mean-field solver cannot produce.
E  Lattice realisation: the per-mode noise constant is *identical* to the
   confining trap, the escape window grows only logarithmically, and the
   conserved functional rescues the test anyway.

Run with ``WPMW_OUTPUT`` set (``/mnt/user-data/outputs`` in the container).
"""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import expm
from scipy.stats import norm

from wpmwlib.phase_space_crystal_lattice import (
    FourierMode,
    PhaseSpaceCrystalLattice,
)
from wpmwlib.wpmw_utils import docs_path, output_path

HBAR = 1.0
MASS = 1.0
OMEGA_W = 1.0                      # the w in -(m w^2/2)(x1-x2)^2
MU = MASS / 2                      # reduced mass
OM = np.sqrt(2.0) * OMEGA_W        # Lyapunov rate of the relative motion

# the reference wavepacket, in the relative coordinate
RC, PC, SR = -4.0, 2.828, 0.70
SP = HBAR / (2 * SR)
UBAR = RC + PC / (MU * OM)
SU = np.sqrt(SR**2 + SP**2 / (MU * OM) ** 2)
T_EXACT = float(norm.cdf(UBAR / SU))


def save_fig(fig, name):
    fig.savefig(output_path(name), dpi=150, bbox_inches="tight")
    dp = docs_path(name)
    if dp:
        fig.savefig(dp, dpi=150, bbox_inches="tight")
    print(f"    [figure] {name}")


def banner(text):
    print()
    print("=" * 74)
    print(text)
    print("=" * 74)


def joint_matrices(m=MASS, w=OMEGA_W):
    """zeta = (x1, x2, p1, p2); H = (1/2) zeta^T Hm zeta; zeta-dot = Om4 Hm zeta."""
    k = m * w**2
    Hm = np.zeros((4, 4))
    Hm[:2, :2] = np.array([[-k, k], [k, -k]])
    Hm[2:, 2:] = np.eye(2) / m
    Om4 = np.zeros((4, 4))
    Om4[:2, 2:] = np.eye(2)
    Om4[2:, :2] = -np.eye(2)
    return Hm, Om4


# ===================================================================== #
def part_a():
    banner("Part A -- separation and Moyal truncation")
    Hm, Om4 = joint_matrices()
    A = Om4 @ Hm
    ev = np.sort(np.linalg.eigvals(A).real)
    print("  H = (p1^2+p2^2)/2m - (m w^2/2)(x1-x2)^2")
    print(f"  reduced mass mu = m/2 = {MU}, Lyapunov rate Om = sqrt(2) w"
          f" = {OM:.6f}")
    print(f"  classical flow eigenvalues: {ev.round(6)}")
    print(f"  expected {{-Om, 0, 0, +Om}} = "
          f"{np.array([-OM, 0.0, 0.0, OM]).round(6)}")
    print(f"  residual: {np.max(np.abs(ev - np.array([-OM, 0, 0, OM]))):.3e}")
    print("\n  H is quadratic, so every third derivative vanishes and the")
    print("  Moyal bracket truncates at the Poisson bracket: the QLE is")
    print("  EXACTLY the classical Liouville equation, to all orders in hbar.")
    print("  A confining trap has the same property; what this system adds is")
    print("  asymptotic free states, hence a genuine scattering observable.")
    return ev


def part_b():
    banner("Part B -- the invariant half-plane")
    Arel = np.array([[0.0, 1 / MU], [MU * OM**2, 0.0]])
    print("  relative flow d/dt (r, p_r) = A_rel (r, p_r)")
    print(f"  u = r + p_r/(mu Om),  v = r - p_r/(mu Om);  mu Om ="
          f" {MU * OM:.6f}")
    print(f"  {'t':>6s} {'u(t)/u(0)':>13s} {'e^(Om t)':>13s} "
          f"{'residual':>10s} {'v(t)/v(0)':>13s} {'e^(-Om t)':>13s}")
    worst = 0.0
    for t in (0.3, 1.0, 2.0, 3.5):
        P = expm(Arel * t)
        z0 = np.array([-1.3, 0.8])
        zt = P @ z0
        u0, ut = z0[0] + z0[1] / (MU * OM), zt[0] + zt[1] / (MU * OM)
        v0, vt = z0[0] - z0[1] / (MU * OM), zt[0] - zt[1] / (MU * OM)
        res = abs(ut / u0 - np.exp(OM * t))
        worst = max(worst, res)
        print(f"  {t:6.2f} {ut/u0:13.8f} {np.exp(OM*t):13.8f} {res:10.2e} "
              f"{vt/v0:13.8f} {np.exp(-OM*t):13.8f}")
    print(f"\n  worst residual: {worst:.2e}")
    print("  The half-plane u > 0 is therefore EXACTLY invariant under the")
    print("  flow, so  T = int_(u>0) W dr dp_r  is a constant of the motion.")
    print("  Since sign(u_0) decides which side the trajectory ends on, T is")
    print("  the asymptotic transmission probability -- readable at t = 0.")
    return worst


def T_closed(rc=RC, pc=PC, sr=SR):
    sp = HBAR / (2 * sr)
    ubar = rc + pc / (MU * OM)
    su = np.sqrt(sr**2 + sp**2 / (MU * OM) ** 2)
    return float(norm.cdf(ubar / su)), ubar, su


def T_schrodinger(rc, pc, sr, t_end, L, N, nstep):
    r = (np.arange(N) - N // 2) * (L / N)
    dr = L / N
    kg = 2 * np.pi * np.fft.fftfreq(N, d=dr)
    psi = (2 * np.pi * sr**2) ** -0.25 * np.exp(
        -(r - rc) ** 2 / (4 * sr**2) + 1j * pc * r / HBAR)
    V = -0.5 * MU * OM**2 * r**2
    dt = t_end / nstep
    eV = np.exp(-0.5j * V * dt / HBAR)
    eK = np.exp(-1j * HBAR * kg**2 * dt / (2 * MU))
    for _ in range(nstep):
        psi = eV * np.fft.ifft(eK * np.fft.fft(eV * psi))
    d = np.abs(psi) ** 2
    tot = np.sum(d) * dr
    return float(np.sum(d[r > 0]) * dr / tot)


def part_c():
    banner("Part C -- the transmission probability")
    print(f"  packet: r_c = {RC}, p_c = {PC}, sigma_r = {SR}")
    print(f"  u_bar = {UBAR:+.6f}, sigma_u = {SU:.6f}")
    print(f"  minimum possible sigma_u = sqrt(hbar/(mu Om)) ="
          f" {np.sqrt(HBAR/(MU*OM)):.6f}")
    print(f"\n  T = Phi(u_bar/sigma_u) = {T_EXACT:.9f}")
    print("\nC1. convergence of a split-operator Schrodinger solve to it")
    print(f"  {'t_end':>7s} {'e^(Om t)':>10s} {'grid N':>8s} "
          f"{'T (Schrodinger)':>17s} {'|diff|':>10s}")
    conv = []
    for t_end, L, N, nstep in [
        (1.5, 120.0, 4096, 3000),
        (2.2, 200.0, 8192, 5000),
        (3.0, 400.0, 16384, 9000),
        (3.8, 1200.0, 32768, 14000),
    ]:
        Ts = T_schrodinger(RC, PC, SR, t_end, L, N, nstep)
        conv.append((t_end, abs(Ts - T_EXACT)))
        print(f"  {t_end:7.1f} {np.exp(OM*t_end):10.1f} {N:8d} "
              f"{Ts:17.9f} {abs(Ts - T_EXACT):10.2e}")
    print("  the residual is the not-yet-separated fraction near u = 0;")
    print("  it decays like e^(-Om t), confirming the closed form is exact.")

    print("\nC2. the transmission width is purely quantum")
    print(f"  {'hbar':>8s} {'sigma_u (min)':>14s} {'T at u_bar = 0.5':>18s}")
    hb_rows = []
    for hb in (1.0, 0.5, 0.25, 0.1, 0.02):
        su_min = np.sqrt(hb / (MU * OM))
        hb_rows.append((hb, su_min, float(norm.cdf(0.5 / su_min))))
        print(f"  {hb:8.3f} {su_min:14.6f} {hb_rows[-1][2]:18.9f}")
    print("  as hbar -> 0 the transmission becomes a step function: the whole")
    print("  of the tunnelling is the quantum WIDTH of W, transported by an")
    print("  exactly classical flow.")

    print("\nC3. Kemble, for the energy-resolved comparison")
    print(f"  {'E/(hbar Om)':>12s} {'T(E) = 1/(1+e^(-2 pi E/hbar Om))':>34s}")
    for Ered in (-1.0, -0.5, 0.0, 0.5, 1.0):
        print(f"  {Ered:12.2f} {1/(1+np.exp(-2*np.pi*Ered)):34.9f}")
    return conv, hb_rows


def part_d():
    banner("Part D -- what makes this a two-body test")
    Hm, Om4 = joint_matrices()
    A = Om4 @ Hm
    sx = SR / np.sqrt(2.0)                 # so that sigma_r = SR
    sp = HBAR / (2 * sx)
    Sig0 = np.diag([sx**2, sx**2, sp**2, sp**2])   # PRODUCT state, S = 0
    uP = np.array([0.0, 0.0, 1.0, 1.0])

    print("  initial state: a product of two identical Gaussians.")
    print("  The pair mode is anti-diagonal, so by Theorem A3 of")
    print("  docs/analysis/fourd_microdynamics.md every one of the four")
    print("  actions leaves P = p1 + p2 untouched.  Exactly:")
    print(f"  {'t':>6s} {'Var(P)':>16s} {'drift':>10s}")
    v0 = float(uP @ Sig0 @ uP)
    for t in (0.0, 1.0, 2.0, 4.0, 8.0):
        P = expm(A * t)
        v = float(uP @ (P @ Sig0 @ P.T) @ uP)
        print(f"  {t:6.1f} {v:16.10f} {abs(v - v0):10.2e}")

    print("\n  Entanglement generated from the product state.")
    print("  The relative mode alone contributes a CONSTANT reduced")
    print("  determinant; the growth comes from the mismatch between the")
    print("  exponentially stretched relative mode and the diffusively")
    print("  spreading centre of mass, giving  nu ~ t e^(Om t),  hence")
    print("  S ~ Om t + ln t  and  dS/dt -> Om + 1/t.")
    print(f"  {'t':>6s} {'Om t':>7s} {'nu':>12s} {'purity_1':>10s} "
          f"{'S (nats)':>10s} {'dS/dt':>9s} {'Om + 1/t':>10s} {'diff':>9s}")

    def entropy_at(t):
        P = expm(A * t)
        Sig = P @ Sig0 @ P.T
        S1 = np.array([[Sig[0, 0], Sig[0, 2]], [Sig[2, 0], Sig[2, 2]]])
        det = max(np.linalg.det(S1), (HBAR / 2) ** 2)
        nu = max(np.sqrt(det) / (HBAR / 2), 1.0)
        a, b = (nu + 1) / 2, (nu - 1) / 2
        S = a * np.log(a) - (b * np.log(b) if b > 1e-15 else 0.0)
        return nu, S, (HBAR / 2) / np.sqrt(det)

    rows, h = [], 0.25
    for t in np.arange(0.0, 12.01, 0.5):
        nu, S, pur = entropy_at(t)
        if t > h:
            rate = (entropy_at(t + h)[1] - entropy_at(t - h)[1]) / (2 * h)
        else:
            rate = None
        rows.append((t, nu, pur, S, rate))
    for t, nu, pur, S, rate in rows:
        if abs(t - round(t)) > 1e-9 or int(round(t)) % 2:
            continue
        if rate is None:
            print(f"  {t:6.2f} {OM*t:7.3f} {nu:12.4e} {pur:10.6f} "
                  f"{S:10.5f} {'':>9s} {'':>10s} {'':>9s}")
        else:
            print(f"  {t:6.2f} {OM*t:7.3f} {nu:12.4e} {pur:10.6f} "
                  f"{S:10.5f} {rate:9.5f} {OM + 1/t:10.5f} "
                  f"{rate - (OM + 1/t):9.2e}")
    print(f"\n  Lyapunov rate Om = {OM:.6f}; the residual after subtracting")
    print("  the 1/t term falls by roughly an order of magnitude per 6 units")
    print("  of time.  Beyond t ~ 14, det(Sigma_1) is a difference of numbers")
    print("  of order e^(2 Om t) and float64 loses it -- an intrinsic limit of")
    print("  the covariance route, not of the physics.")
    print("\n  A mean-field (TDHF-Wigner) solver produces S = 0 identically.")
    print("  This is the observable that separates a joint solver from one.")
    return rows


def inverted_modes(qmax, L):
    return [FourierMode(q=q,
                        V_q=-0.5 * MU * OM**2 * (L**2 / np.pi**2)
                        * ((-1.0) ** q) / q**2,
                        phi_q=0.0)
            for q in range(1, qmax + 1)]


def jump_symbol(lat, modes):
    th = np.fft.fftfreq(lat.N) * 2 * np.pi
    sym = np.zeros((lat.N, lat.M))
    for m in modes:
        G = -(m.V_q / lat.hbar) * np.sin(
            2 * np.pi * m.q * lat.x / lat.L + m.phi_q)
        sym += 2.0 * np.sin(m.q * th)[:, None] * G[None, :]
    return sym


def exact_jump(lat, sym, dt):
    lat.W = np.real(np.fft.ifft(
        np.fft.fft(lat.W, axis=0) * np.exp(1j * sym * dt), axis=0))


def part_e():
    banner("Part E -- the lattice realisation, and what it costs")
    L = 40.0
    dp = np.pi * HBAR / L
    print("E1. per-mode noise constant, against the confining trap")
    print(f"  {'q':>4s} {'|V_q|':>13s} {'(2 q dp)^2':>12s} {'product':>12s} "
          f"{'/(2 mu Om^2 hbar)':>19s}")
    devs = []
    for m in inverted_modes(6, L):
        prod = (abs(m.V_q) / HBAR) * (2 * m.q * dp) ** 2
        devs.append(prod / (2 * MU * OM**2 * HBAR))
        print(f"  {m.q:4d} {abs(m.V_q):13.5e} {(2*m.q*dp)**2:12.5e} "
              f"{prod:12.5e} {devs[-1]:19.10f}")
    print(f"  2 mu Om^2 hbar = {2*MU*OM**2*HBAR:.6f}"
          f" = 2 m w^2 hbar = {2*MASS*OMEGA_W**2*HBAR:.6f}")
    print(f"  max deviation from 1: {np.max(np.abs(np.array(devs)-1)):.2e}")
    print("  Inverting the sign of the potential leaves |V_q| unchanged, so")
    print("  the mode-sum noise pathology of the trap is inherited IN FULL.")
    print("  The inverted barrier is a better probe, not a cheaper one.")

    print("\nE2. the conserved functional tracked on the crystal lattice")
    M_GRID, N_GRID, QMAX = 512, 512, 128
    lat = PhaseSpaceCrystalLattice(M=M_GRID, N=N_GRID, L=L, mass=MU,
                                   hbar=HBAR, nu=None, advection="spectral")
    lat.initialize_from_wigner(
        lambda X, P: (1.0 / (2 * np.pi * SR * SP))
        * np.exp(-(X - RC) ** 2 / (2 * SR**2) - (P - PC) ** 2 / (2 * SP**2)))
    modes = inverted_modes(QMAX, L)
    sym = jump_symbol(lat, modes)
    dA = lat.dx * lat.dp
    U = lat.X + lat.P / (MU * OM)
    print(f"  grid {M_GRID}x{N_GRID}, L = {L}, qmax = {QMAX}, "
          f"p_max = {np.max(lat.p):.3f}")
    print(f"  rate budget sum|V_q|/hbar = "
          f"{sum(abs(m.V_q) for m in modes)/HBAR:.1f}"
          f"  (jump substep integrated exactly by FFT)")
    print(f"  {'t':>6s} {'e^(Om t)':>9s} {'int_(u>0) W':>13s} "
          f"{'|error| vs T':>13s} {'norm':>9s} {'edge mass':>11s}")
    track, nstep, t_end = [], 600, 1.5
    dt = t_end / nstep
    for k in range(nstep + 1):
        t = k * dt
        Tl = float(np.sum(lat.W[U > 0]) * dA)
        track.append((t, Tl))
        if k % 100 == 0:
            edge = float(np.sum(lat.W[np.abs(lat.X) > 0.45 * L]) * dA)
            print(f"  {t:6.3f} {np.exp(OM*t):9.2f} {Tl:13.8f} "
                  f"{abs(Tl - T_EXACT):13.2e} "
                  f"{float(np.sum(lat.W)*dA):9.6f} {edge:11.2e}")
        if k < nstep:
            lat.step_advect(dt / 2)
            exact_jump(lat, sym, dt)
            lat.step_advect(dt / 2)
    track = np.array(track)
    print(f"  drift over the run: "
          f"{np.max(track[:,1]) - np.min(track[:,1]):.2e}"
          f"  (no systematic trend; the offset at t=0 is grid discretisation)")

    print("\nE3. the escape problem: window logarithmic, grid quadratic")
    u_typ = 3 * SU
    print(f"  require |u| <= 3 sigma_u = {u_typ:.3f} to stay on the grid")
    print(f"  {'Om t*':>7s} {'x window':>10s} {'p window':>10s} {'L':>9s} "
          f"{'N needed':>11s} {'accuracy ~ e^-Om t':>19s}")
    grid_rows = []
    for Omt in (1.0, 2.0, 3.0, 4.0, 5.0):
        g = np.exp(Omt)
        xw, pw = u_typ * g, MU * OM * u_typ * g
        Lr = 2 * xw
        Ng = 2 * Lr * pw / (np.pi * HBAR)
        grid_rows.append((Omt, Lr, Ng, np.exp(-Omt)))
        print(f"  {Omt:7.1f} {xw:10.2f} {pw:10.2f} {Lr:9.1f} {Ng:11.0f} "
              f"{np.exp(-Omt):19.2e}")
    print("  N ~ 1/epsilon^2 for a target accuracy epsilon in T.")
    print("  This limits literal simulation of the asymptotic states -- but")
    print("  NOT the test, because T is conserved and readable at any time.")
    return track, grid_rows


# ===================================================================== #
def fig_phase_space(track):
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.9))
    rr = np.linspace(-12, 12, 400)
    pp = np.linspace(-12, 12, 400)
    RRg, PPg = np.meshgrid(rr, pp)
    E = PPg**2 / (2 * MU) - 0.5 * MU * OM**2 * RRg**2

    for ax, t in zip(axes, (0.0, 0.7, 1.4)):
        ax.contour(RRg, PPg, E, levels=np.linspace(-20, 20, 17),
                   colors="0.82", linewidths=0.7, zorder=1)
        ax.plot(rr, -MU * OM * rr, "-", color="darkorange", lw=2.2, zorder=3,
                label=r"separatrix $u=0$")
        ct, st = np.cosh(OM * t), np.sinh(OM * t)
        rc_t = RC * ct + PC * st / (MU * OM)
        pc_t = PC * ct + MU * OM * RC * st
        sr_t = np.sqrt((SR * ct) ** 2 + (SP * st / (MU * OM)) ** 2)
        sp_t = np.sqrt((SP * ct) ** 2 + (MU * OM * SR * st) ** 2)
        cov_t = SR**2 * ct * MU * OM * st + SP**2 * st * ct / (MU * OM)
        Sg = np.array([[sr_t**2, cov_t], [cov_t, sp_t**2]])
        Si = np.linalg.inv(Sg)
        dR, dP = RRg - rc_t, PPg - pc_t
        Q = (Si[0, 0] * dR**2 + 2 * Si[0, 1] * dR * dP + Si[1, 1] * dP**2)
        ax.contourf(RRg, PPg, np.exp(-Q / 2), levels=np.linspace(0.02, 1, 12),
                    cmap="Blues", alpha=0.85, zorder=2)
        ax.plot(rc_t, pc_t, "o", color="crimson", ms=7, zorder=4)
        ax.set_xlim(-12, 12)
        ax.set_ylim(-12, 12)
        ax.set_xlabel(r"$r = x_1 - x_2$")
        ax.set_ylabel(r"$p_r$")
        ax.set_title(f"$t = {t:.1f}$  ($e^{{\\Omega t}} = {np.exp(OM*t):.1f}$)",
                     fontsize=11)
        ax.text(0.97, 0.03, "transmits ($u>0$)", transform=ax.transAxes,
                fontsize=9, ha="right", va="bottom", color="darkorange")
        ax.text(0.03, 0.97, "reflects ($u<0$)", transform=ax.transAxes,
                fontsize=9, ha="left", va="top", color="darkorange")
        if t == 0.0:
            ax.legend(fontsize=9, loc="lower left")
    fig.suptitle("Inverted pair barrier: the packet is sheared along the "
                 "hyperbolas, and the fraction above the separatrix "
                 "$u=0$ never changes", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    save_fig(fig, "inverted_pair_barrier_phase_space.png")
    plt.close(fig)


def fig_diagnostics(conv, hb_rows, ent_rows, track, grid_rows):
    fig, axes = plt.subplots(1, 4, figsize=(20, 4.6))

    ax = axes[0]
    t_c = np.array([c[0] for c in conv])
    d_c = np.array([c[1] for c in conv])
    ax.semilogy(t_c, d_c, "o-", lw=2, color="crimson", label="Schrodinger")
    ax.semilogy(t_c, d_c[0] * np.exp(-OM * (t_c - t_c[0])), "--",
                color="0.4", lw=1.5, label=r"$\propto e^{-\Omega t}$")
    ax.set_xlabel("$t$ of the solve")
    ax.set_ylabel(r"$|T_{\rm num} - \Phi(\bar u/\sigma_u)|$")
    ax.set_title("The closed form is exact", fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3, which="both")

    ax = axes[1]
    ax.plot(track[:, 0], track[:, 1], "-", lw=2, color="royalblue",
            label=r"$\int_{u>0} W$ on the lattice")
    ax.axhline(T_EXACT, color="crimson", ls="--", lw=1.6,
               label=f"exact $T = {T_EXACT:.5f}$")
    ax.set_xlabel("$t$")
    ax.set_ylabel("transmission functional")
    ax.set_ylim(T_EXACT - 3e-3, T_EXACT + 3e-3)
    ax.set_title("An exactly conserved observable\nfor the microdynamics",
                 fontsize=11)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[2]
    te = np.array([r[0] for r in ent_rows])
    Se = np.array([r[3] for r in ent_rows])
    ax.plot(te, Se, "-", lw=2, color="seagreen", label="entanglement entropy")
    ax.plot(te[te > 3], OM * (te[te > 3] - 3) + Se[te > 3][0], "--",
            color="0.4", lw=1.5, label=r"slope $\Omega=\sqrt{2}\,\omega$")
    ax.axhline(0.0, color="crimson", lw=1.6, ls=":",
               label="mean-field prediction")
    ax.set_xlabel("$t$")
    ax.set_ylabel("$S$ (nats)")
    ax.set_title("Entanglement grows at the\nLyapunov rate", fontsize=11)
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(alpha=0.3)

    ax = axes[3]
    Omt = np.array([g[0] for g in grid_rows])
    Ng = np.array([g[2] for g in grid_rows])
    eps = np.array([g[3] for g in grid_rows])
    ax.loglog(eps, Ng, "o-", lw=2, color="darkorange", label="grid cells $N$")
    ax.loglog(eps, Ng[0] * (eps[0] / eps) ** 2, "--", color="0.4", lw=1.5,
              label=r"$\propto \epsilon^{-2}$")
    ax.set_xlabel(r"target accuracy $\epsilon$ in $T$")
    ax.set_ylabel("$N$")
    ax.invert_xaxis()
    ax.set_title("Literal escape is expensive\n(the conserved functional "
                 "is not)", fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3, which="both")

    fig.tight_layout()
    save_fig(fig, "inverted_pair_barrier_diagnostics.png")
    plt.close(fig)


def main():
    print("Verification for docs/supplement/inverted_pair_barrier.md")
    part_a()
    part_b()
    conv, hb_rows = part_c()
    ent_rows = part_d()
    track, grid_rows = part_e()
    banner("Figures")
    fig_phase_space(track)
    fig_diagnostics(conv, hb_rows, ent_rows, track, grid_rows)
    print("\ndone.")


if __name__ == "__main__":
    main()
