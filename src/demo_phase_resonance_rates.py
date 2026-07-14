"""
Demonstration: phase-resonance microdynamics of the crystal-lattice sea.

Companion to ``docs/analysis/phase_resonance_microdynamics.md``, the sequel
to ``docs/analysis/sea_dressed_microdynamics.md``. World-particles carry a
de Broglie phase (postulates P0-P3); the sea's polarization is a coherent
grating woven from pair beats; the hop channels K3/K4 are Bragg diffraction
off that grating; and the channel rate law follows from one vertex postulate
P5. Four parts verify the note's claims:

  A. **Beat kinematics.** The beat of two legs split by 2q*dp is a
     full-contrast grating at exactly the mode-q wavelength, drifting at
     the mean velocity pbar/m (Proposition 1).

  B. **Bragg selection.** One split-operator kick of a plane wave scatters
     to exactly p0 +- 2q*dp with amplitude i*Vp*dt/(2*hbar) per sideband:
     the physical transfer is the full quantum 2q*dp, linear in Vp, with
     the refractive (phase-grating) factor i.

  C. **Midpoint identity.** The exact first-order Wigner change of the
     kicked plane wave equals the single-cosine QLE stencil
     dW_n/dt = Gamma_q(x) (W_{n+q} - W_{n-q}) at the midpoint sites
     n0 +- q, with deviation O(dt^2): the stencil's half-quantum offsets
     are the interference midpoints of the full-quantum transfer
     (Theorem 1's bookkeeping).

  D. **Rate-table toy.** Sea pairs carry particle-level data only
     (location, row, imprinted beat contrast and phase); P1 winding
     evolves phases; the P5 vertex weight w = (1 + C cos delta)/2 with
     delta = the struck beat's pattern phase at the vertex. Verified:
       T1  per-channel rate = tau_p * gamma(x)/2 with sigma = sign Gamma
           (quadrature, gamma/2, and direction derived, Theorem 3);
       T2  pattern-phase spread across pairs pumped at random locations
           is zero to machine precision (Lemma 2, coherence for free);
       T3  time-averaging along the transition worldline dephases every
           row except nbar = n_hi (Proposition 2 as dynamics).

Figures (written via docs_path for the output branch):
  phase_grating_sea.png        - one pair's beat; coherence; quadrature
  phase_grating_spacetime.png  - beat crests and worldline refraction
  phase_resonance_rate_law.png - derived rate law and row selection
"""

from __future__ import annotations

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from wpmwlib.wpmw_utils import output_path, docs_path  # noqa: E402

HBAR = 1.0
M_PART = 1.0
L = 8.0
DP = np.pi * HBAR / L          # Wigner half-grid spacing
Q = 2
K = 2.0 * np.pi * Q / L        # grating wavevector = 2*Q*DP/HBAR
VP = 1.5
CORAL, TEAL, GREY = "#e2725b", "#0f8b8d", "#666666"


def save_fig(fig, name: str) -> None:
    fig.savefig(output_path(name), dpi=150, bbox_inches="tight")
    dp = docs_path(name)
    if dp:
        fig.savefig(dp, dpi=150, bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------
# Part A: beat kinematics
# --------------------------------------------------------------------------
def part_a() -> None:
    print("=" * 72)
    print("Part A: beat kinematics (Proposition 1)")
    m_grid = 512
    x = np.arange(m_grid) * L / m_grid
    n_hi, n_lo = 12, 6            # even sites, separation 2Q (Q=3 here)
    q_loc = (n_hi - n_lo) // 2
    p_hi, p_lo = n_hi * DP, n_lo * DP
    w_hi = p_hi**2 / (2 * M_PART * HBAR)
    w_lo = p_lo**2 / (2 * M_PART * HBAR)
    vbar = (p_hi + p_lo) / (2 * M_PART)
    dk = (p_hi - p_lo) / HBAR

    ts = np.linspace(0.0, 0.8, 9)
    phases, spurious = [], 0.0
    for t in ts:
        beat = np.abs(np.exp(1j * (p_hi * x / HBAR - w_hi * t))
                      - np.exp(1j * (p_lo * x / HBAR - w_lo * t))) ** 2
        f = np.fft.rfft(beat) / m_grid
        phases.append(np.angle(f[q_loc]))
        other = np.abs(f).copy()
        other[[0, q_loc]] = 0.0
        spurious = max(spurious, other.max())
    v_fit = -(np.polyfit(ts, np.unwrap(np.array(phases)), 1)[0]) / dk
    print(f"  single spatial mode q={q_loc}: max spurious {spurious:.1e}")
    print(f"  envelope drift {v_fit:.12f} vs pbar/m {vbar:.12f}  "
          f"(rel err {abs(v_fit - vbar) / vbar:.1e})")


# --------------------------------------------------------------------------
# Part B: Bragg selection
# --------------------------------------------------------------------------
def kicked_spectrum(dt: float, n0: int, m_grid: int = 512) -> np.ndarray:
    x = np.arange(m_grid) * L / m_grid
    v_of_x = -VP * np.cos(K * x)
    psi = np.exp(-1j * v_of_x * dt / HBAR) * np.exp(1j * n0 * DP * x / HBAR)
    return np.fft.fft(psi) / m_grid          # bin j <-> even site n = 2j


def part_b() -> None:
    print("=" * 72)
    print("Part B: Bragg selection (transfer exactly 2q*dp, linear in Vp)")
    n0 = 8
    for dt in (1e-3, 1e-4):
        c = kicked_spectrum(dt, n0)
        b0, bp, bm = n0 // 2, (n0 + 2 * Q) // 2, (n0 - 2 * Q) // 2
        pred = 1j * VP * dt / (2 * HBAR)
        resid = np.abs(c).copy()
        resid[[b0, bp, bm]] = 0.0
        print(f"  dt={dt:.0e}: sideband/pred (+2q dp) {c[bp] / pred:.6f}  "
              f"(-2q dp) {c[bm] / pred:.6f}  other bins {resid.max():.1e}")


# --------------------------------------------------------------------------
# Part C: midpoint identity against the QLE stencil
# --------------------------------------------------------------------------
def wigner_of_spectrum(c: np.ndarray, x: np.ndarray,
                       tol: float = 1e-14) -> dict[int, np.ndarray]:
    """W(site n, x) for psi = sum_j c_j exp(i 2j dp x / hbar)."""
    idx = np.nonzero(np.abs(c) > tol)[0]
    w: dict[int, np.ndarray] = {}
    for a in idx:
        for b in idx:
            nmid = int(a + b)                 # (2a + 2b)/2 in dp units
            term = c[a] * np.conj(c[b]) * np.exp(1j * 2 * np.pi
                                                 * (a - b) * x / L)
            w[nmid] = w.get(nmid, np.zeros(len(x), complex)) + term
    return {n: arr.real for n, arr in w.items()}


def part_c() -> None:
    print("=" * 72)
    print("Part C: midpoint identity (stencil offsets +-q are beat midpoints)")
    m_grid = 512
    x = np.arange(m_grid) * L / m_grid
    n0 = 8
    gam = (VP / HBAR) * np.sin(K * x)
    zero = np.zeros(m_grid)
    for dt in (1e-3, 1e-4):
        c1 = kicked_spectrum(dt, n0)
        c0 = np.zeros(m_grid, complex)
        c0[n0 // 2] = 1.0
        w1 = wigner_of_spectrum(c1, x)
        w0 = wigner_of_spectrum(c0, x)
        dev = 0.0
        for n in sorted(set(w1) | set(w0)):
            dw = w1.get(n, zero) - w0.get(n, zero)
            sten = dt * gam * (w0.get(n + Q, zero) - w0.get(n - Q, zero))
            dev = max(dev, float(np.max(np.abs(dw - sten))))
        print(f"  dt={dt:.0e}: max |dW - stencil| = {dev:.2e}  "
              f"(dev/dt^2 = {dev / dt**2:.3f})")


# --------------------------------------------------------------------------
# Part D: rate-table toy (Theorem 3)
# --------------------------------------------------------------------------
def part_d() -> tuple[np.ndarray, ...]:
    print("=" * 72)
    print("Part D: rate law from co-location phase statistics under P5")
    rng = np.random.default_rng(20260714)
    n_pairs = 60000
    tau_p = 1e-3
    n_hi = 6
    x0 = rng.uniform(0, L, n_pairs)
    rows = 2 * rng.integers(-4, 5, n_pairs)
    contrast = VP * tau_p / HBAR

    def pattern_phase(x: float, t: float, s: int) -> np.ndarray:
        """Lambda_j^s(x,t) from carried data: imprint (pi/2 + s K x0) plus
        propagation s K (x - x0) minus P1 winding s K vbeat t."""
        vbeat = (rows + s * Q) * DP / M_PART
        return (s * K * (x - x0) + (np.pi / 2 + s * K * x0)
                - s * K * vbeat * t)

    def channel_rate(x: float, d: int, t_avg: float = 0.0, nt: int = 60,
                     row: int | None = None) -> float:
        """Phase-sensitive part of the P5 rate for direction d (family s=d),
        averaged over [0, t_avg] along the transition worldline."""
        s = d
        v_mid = (n_hi + d * Q) * DP / M_PART
        sel = slice(None) if row is None else (rows == row)
        ts = np.linspace(0, t_avg, nt) if t_avg > 0 else np.array([0.0])
        acc = 0.0
        for t in ts:
            acc += float(np.mean(np.cos(pattern_phase(x + v_mid * t, t,
                                                      s)[sel])))
        return 0.5 * contrast * acc / len(ts)

    xs = np.linspace(0.1, L - 0.1, 40)
    gam = (VP / HBAR) * np.sin(K * xs)
    tgt = tau_p * gam / 2
    r_dn = np.array([channel_rate(x, -1) for x in xs])
    r_up = np.array([channel_rate(x, +1) for x in xs])
    msk = np.abs(np.sin(K * xs)) > 0.05
    print(f"  T1 down-channel: max |R/(tau gamma/2) - 1| = "
          f"{np.max(np.abs(r_dn[msk] / tgt[msk] - 1)):.2e}")
    print(f"  T1 up-channel  : max |R/(-tau gamma/2) - 1| = "
          f"{np.max(np.abs(r_up[msk] / (-tgt[msk]) - 1)):.2e}")
    print(f"  T1 net = tau_p*Gamma: max dev "
          f"{np.max(np.abs((r_dn - r_up) - tau_p * gam)):.2e}; "
          f"sigma == sign(Gamma): "
          f"{bool(np.all(np.sign(r_dn - r_up) == np.sign(gam)))}")

    ph = np.mod(pattern_phase(3.123, 0.0, -1), 2 * np.pi)
    print(f"  T2 pattern-phase spread over {n_pairs} pairs: "
          f"std = {np.std(ph):.2e}")

    t_avg = 40.0
    x_probe = L / (4 * Q)
    row_vals = sorted(set(rows.tolist()))
    inst = [channel_rate(x_probe, -1, row=r) / (0.5 * contrast)
            for r in row_vals]
    avg = [channel_rate(x_probe, -1, t_avg=t_avg, nt=400, row=r)
           / (0.5 * contrast) for r in row_vals]
    for r, i_v, a_v in zip(row_vals, inst, avg):
        mark = "   <-- resonant (= n_hi)" if r == n_hi else ""
        print(f"  T3 row {r:+d}: inst {i_v:+.3f} -> avg {a_v:+.3f}{mark}")
    return xs, gam, tgt, r_dn, r_up, np.array(row_vals), np.array(avg)


# --------------------------------------------------------------------------
# Figures
# --------------------------------------------------------------------------
def figure_sea() -> None:
    x = np.linspace(0, L, 800)
    fig, ax = plt.subplots(1, 3, figsize=(13.5, 3.6))

    n_hi, n_lo = 10, 6
    kp, km = n_hi * DP / HBAR, n_lo * DP / HBAR
    ax[0].plot(x, np.cos(kp * x), color=CORAL, lw=0.9,
               label="positon leg $p_+$")
    ax[0].plot(x, -np.cos(km * x), color=TEAL, lw=0.9, ls="--",
               label="negaton leg $p_-$ ($\\pi$ offset)")
    ax[0].plot(x, np.abs(np.exp(1j * kp * x) - np.exp(1j * km * x)) ** 2 / 2
               - 2.6, color="k", lw=2.0,
               label="beat $|\\psi_+ - \\psi_-|^2$ (shifted)")
    ax[0].annotate("", xy=(L / Q, -0.55), xytext=(0, -0.55),
                   arrowprops=dict(arrowstyle="<->", color="k"))
    ax[0].text(L / (2 * Q), -0.45, "$\\lambda = L/q$", ha="center",
               fontsize=10)
    ax[0].set_title("(a) one split pair: full-contrast beat")
    ax[0].set_xlabel("x")
    ax[0].set_yticks([])
    ax[0].legend(fontsize=7, loc="upper right")

    rng = np.random.default_rng(20260712)
    n_show = 400
    phis = rng.uniform(0, 2 * np.pi, n_show)
    inc = np.mean(-np.cos(K * x[None, :] - phis[:, None]), axis=0)
    ax[1].plot(x, -np.cos(K * x - np.pi), color="k", lw=2.0,
               label="coherent lock (contrast 1)")
    ax[1].plot(x, inc, color=GREY, lw=1.2,
               label=f"random phases (contrast ~ N^(-1/2), N={n_show})")
    ax[1].axhline(0, color=GREY, lw=0.5)
    ax[1].set_title("(b) grating needs coherence: per-pair-mean modulation")
    ax[1].set_xlabel("x")
    ax[1].legend(fontsize=8)

    ax[2].plot(x, -np.cos(K * x), color=GREY, ls="--", lw=1.5,
               label="$V(x) = -V_p\\cos(Kx)$")
    ax[2].plot(x, np.sin(K * x), color="k", lw=2.0,
               label="$\\Gamma_q(x) \\propto \\sin(Kx)$")
    ax[2].fill_between(x, 0, np.sin(K * x), where=np.sin(K * x) > 0,
                       color=CORAL, alpha=0.18)
    ax[2].fill_between(x, 0, np.sin(K * x), where=np.sin(K * x) < 0,
                       color=TEAL, alpha=0.18)
    ax[2].text(0.55 * L, 0.75, "$\\sigma = +1$", color=CORAL, fontsize=10)
    ax[2].text(0.31 * L, -0.85, "$\\sigma = -1$", color=TEAL, fontsize=10)
    ax[2].axhline(0, color=GREY, lw=0.5)
    ax[2].set_title("(c) hops read the quadrature of V")
    ax[2].set_xlabel("x")
    ax[2].legend(fontsize=8, loc="upper right")
    fig.tight_layout()
    save_fig(fig, "phase_grating_sea.png")


def figure_spacetime() -> None:
    fig, ax = plt.subplots(1, 2, figsize=(12.5, 5.4), sharey=True)
    t_max = 1.6
    tg = np.linspace(0, t_max, 300)
    xg = np.linspace(0, L, 400)
    xx, tt = np.meshgrid(xg, tg)

    x0 = 2.0
    nh, nl = 4, 0
    vh, vl = nh * DP / M_PART, nl * DP / M_PART
    vbar = 0.5 * (vh + vl)
    beat = 1 - np.cos(K * (xx - x0 - vbar * tt))
    ax[0].pcolormesh(xx, tt, beat, cmap="Greys", vmin=0, vmax=2.6,
                     shading="auto")
    ax[0].plot(x0 + vh * tg, tg, color=CORAL, lw=2.2,
               label=f"positon leg  n={nh}")
    ax[0].plot(x0 + vl * tg, tg, color=TEAL, lw=2.2, ls="--",
               label=f"negaton leg  n={nl}")
    ax[0].plot(x0 + vbar * tg, tg, color="k", lw=1.4, ls=":",
               label=f"CoM / beat drift  nbar={(nh + nl) // 2}")
    ax[0].set_title("(a) split pair: beat crests drift at the CoM slope")
    ax[0].set_xlabel("x")
    ax[0].set_ylabel("t")
    ax[0].legend(fontsize=8, loc="upper left", framealpha=0.9)

    ax[1].pcolormesh(xx, tt, 1 - np.cos(K * xx), cmap="Greys", vmin=0,
                     vmax=2.6, shading="auto")
    n_hi, n_lo = 6, 2
    n_mid = (n_hi + n_lo) // 2
    v_hi, v_lo = n_hi * DP / M_PART, n_lo * DP / M_PART
    v_mid = n_mid * DP / M_PART
    t1, t2 = 0.49, 0.79
    xa = 0.5
    xb = xa + v_hi * t1
    xc = xb + v_mid * (t2 - t1)
    seg1 = np.linspace(0, t1, 50)
    seg2 = np.linspace(t1, t2, 30)
    seg3 = np.linspace(t2, t_max, 60)
    ax[1].plot(xa + v_hi * seg1, seg1, color=CORAL, lw=2.4)
    ax[1].plot(xb + v_mid * (seg2 - t1), seg2, color="k", lw=3.0)
    ax[1].plot(xc + v_lo * (seg3 - t2), seg3, color=CORAL, lw=2.4)
    for c in np.linspace(-1.2, 1.2, 5):
        ax[1].plot(xb + c + v_mid * (seg2 - t1), seg2, color=TEAL, lw=0.8,
                   alpha=0.75)
    ax[1].annotate("", xy=(xb + v_mid * 0.15 - 0.85, t1 + 0.15),
                   xytext=(xb + v_mid * 0.15, t1 + 0.15),
                   arrowprops=dict(arrowstyle="->", color="k", lw=1.4))
    ax[1].text(xb - 2.35, t1 + 0.13, "recoil $2q\\thinspace dp$", fontsize=9)
    ax[1].text(xa + v_hi * 0.22 + 0.12, 0.16, f"$n_{{hi}}$={n_hi}",
               color=CORAL, fontsize=10)
    ax[1].text(xc + v_lo * 0.5 + 0.14, 1.15, f"$n_{{lo}}$={n_lo}",
               color=CORAL, fontsize=10)
    ax[1].text(xb + 0.65, (t1 + t2) / 2 + 0.10,
               f"transition: slope nbar={n_mid}\n(stencil center; parallel"
               " to\nabsorbed beat crests)", fontsize=8.5)
    ax[1].set_title("(b) hop = worldline refraction at the grating")
    ax[1].set_xlabel("x")
    fig.tight_layout()
    save_fig(fig, "phase_grating_spacetime.png")


def figure_rate_law(xs, gam, tgt, r_dn, r_up, row_vals, row_avg) -> None:
    tau_p = 1e-3
    fig, ax = plt.subplots(1, 2, figsize=(12.0, 4.2))
    ax[0].plot(xs, tgt, color="k", lw=2.4,
               label="target $\\tau_p\\Gamma_q(x)/2$")
    ax[0].plot(xs, r_dn, "o", ms=4, color=CORAL,
               label="down-channel (P5 statistics)")
    ax[0].plot(xs, -r_up, "s", ms=3.5, color=TEAL,
               label="$-$up-channel (P5 statistics)")
    ax[0].axhline(0, color=GREY, lw=0.5)
    ax[0].set_title("derived rate law: quadrature, sign, $\\gamma/2$")
    ax[0].set_xlabel("x")
    ax[0].set_ylabel("rate x $\\tau_p$")
    ax[0].legend(fontsize=8)

    ax[1].bar(row_vals, row_avg, width=1.2, color=[
        CORAL if r == 6 else GREY for r in row_vals])
    ax[1].set_title("row selection by dephasing (avg along transition)")
    ax[1].set_xlabel("sea pair row $\\bar{n}$")
    ax[1].set_ylabel("normalized channel amplitude")
    ax[1].text(6, 0.9, "$\\bar{n} = n_{hi}$", ha="center", color=CORAL,
               fontsize=10)
    fig.tight_layout()
    save_fig(fig, "phase_resonance_rate_law.png")


if __name__ == "__main__":
    part_a()
    part_b()
    part_c()
    d_out = part_d()
    figure_sea()
    figure_spacetime()
    figure_rate_law(*d_out)
    print("=" * 72)
    print("figures written:",
          "phase_grating_sea.png, phase_grating_spacetime.png, "
          "phase_resonance_rate_law.png")
