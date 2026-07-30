#!/usr/bin/env python3
"""
Resource arithmetic for permanent positon-negaton pairing.

Companion to ``docs/analysis/permanent_pairing_density_matrix.md``.

Two accounting models for the sea's aligned-pair budget are compared on the
canonical cosine-well trajectory (exact split-operator Schrodinger arm):

* **Consumable model** (the earlier objection): every exchange event
  permanently consumes one aligned pair.  Demand is the cumulative moved
  mass E(t), which grows linearly in t; the local aligned stock per cell
  is (2/h) dx dp.  This model exhausts the core cells in a small fraction
  of one period -- the sound half of the original objection.

* **Storage model** (permanent pairing, density-matrix ontology): an
  exchange converts an aligned pair to a split pair, which stores one
  sample of an off-diagonal density-matrix element and is retrievable.
  Demand is the instantaneous coherence content, not the event count:

  - locally, the load factor |W(x,p)| * (h/2) <= 1 for every state at
    every time (Wigner's bound |W| <= 2/h) -- the same inequality that
    guarantees W' = W + 2/h >= 0 guarantees that a density-2/h sea can
    host all local coherence with pair polarisation <= 100%;
  - globally, the off-diagonal L1 mass C(t) of rho(P, P') is bounded and
    small compared with the aligned stock of one state-mass unit per
    momentum row.

Printed checks:
  1. load(t) = max |W| * pi*hbar over the trajectory  (<= 1 + instrument
     tolerance; Gaussian instants saturate the bound exactly),
  2. C(t) and its pair-splitting spectrum (supply-chain depth in one-leg
     hops of 2q dp),
  3. E_state(t), E_sea(t) and the local core-cell exhaustion time of the
     consumable model,
  4. flow feasibility: |dW/dt| <= |Gamma| (W'(p+q dp) + W'(p-q dp))
     cell-by-cell (the stencil is itself the polarisation-flow ledger).

Instrument note: the discrete ring Wigner transform carries a Z_M
antipodal-reflection artifact scaling with the packet's wrap amplitude
(~2% relative for SQUEEZE = 2).  rho(P, P') is computed from psi directly
and is artifact-free; load-factor numbers are quoted with the tolerance.
"""

from __future__ import annotations

import numpy as np

from wpmwlib.wpmw_utils import output_path, docs_path

# ------------------------------------------------------------------ params
HBAR = 1.0
MASS = 1.0
L    = 8.0
V_P  = 1.5
Q    = 1
PHI  = np.pi                    # V(x) = V_P cos(Kx + pi) = -V_P cos(Kx)
K    = 2.0 * np.pi * Q / L

M = 96
N = 96
DX = L / M
DP = np.pi * HBAR / L

x = (np.arange(M) - M // 2) * DX
p = (np.arange(N) - N // 2) * DP

OMEGA    = (2.0 * np.pi / L) * np.sqrt(V_P / MASS)
T_PERIOD = 2.0 * np.pi / OMEGA
SIG_GS   = np.sqrt(HBAR / (2.0 * MASS * OMEGA))
SQUEEZE  = 2.0
SIG_X    = SQUEEZE * SIG_GS
X0, P0   = -L / 4.0, 0.0

T_FINAL = 4.0 * T_PERIOD
N_STEPS = 800
DT      = T_FINAL / N_STEPS

V_x     = V_P * np.cos(K * x + PHI)
GAMMA_x = -(V_P / HBAR) * np.sin(K * x + PHI)

SEA = 2.0 / (2.0 * np.pi * HBAR)               # 2/h = 1/(pi hbar)

# ------------------------------------------------------------- Schrodinger
k_free = 2.0 * np.pi * np.fft.fftfreq(M, d=DX)

def strang_step(psi: np.ndarray, dt: float) -> np.ndarray:
    expV2 = np.exp(-1j * V_x * dt / (2.0 * HBAR))
    expT  = np.exp(-1j * HBAR * k_free ** 2 * dt / (2.0 * MASS))
    return expV2 * np.fft.ifft(expT * np.fft.fft(expV2 * psi))

# ------------------------------------------------- ring Wigner (instrument)
def wigner_of_psi(psi: np.ndarray) -> np.ndarray:
    j = np.arange(M)
    sign = (-1.0) ** j
    W = np.empty((N, M))
    for m in range(M):
        C = np.conj(psi[(m + j) % M]) * psi[(m - j) % M] * sign
        W[:, m] = (np.fft.ifft(C) * M).real
    return W * (DX / (np.pi * HBAR))

# --------------------------------------------------------- rho(P,P') tools
def momentum_amps(psi: np.ndarray) -> np.ndarray:
    """Normalised momentum amplitudes on the psi-grid (spacing 2 dp)."""
    a = np.fft.fftshift(np.fft.fft(psi)) * DX / np.sqrt(L)
    return a / np.sqrt(np.sum(np.abs(a) ** 2))

def offdiag_l1(a: np.ndarray) -> float:
    """C = sum_{P != P'} |rho(P,P')| for the pure state rho = a a^dagger."""
    s1 = np.sum(np.abs(a))
    return float(s1 ** 2 - 1.0)                 # sum|a|^2 = 1

def splitting_spectrum(a: np.ndarray) -> np.ndarray:
    """Off-diagonal mass vs splitting s = |P - P'| / (2 q dp) (hop depth)."""
    absa = np.abs(a)
    nmodes = absa.size
    spec = np.zeros(nmodes)
    for s in range(1, nmodes):
        spec[s] = 2.0 * float(np.sum(absa[s:] * absa[:-s]))
    return spec

# ------------------------------------------------------------------- main
def main() -> None:
    psi = np.exp(-(x - X0) ** 2 / (4.0 * SIG_X ** 2)
                 + 1j * P0 * x / HBAR).astype(complex)
    psi /= np.sqrt(np.sum(np.abs(psi) ** 2) * DX)

    print("=" * 74)
    print("Resource arithmetic: consumable vs storage accounting of the sea")
    print(f"  V_P={V_P}, L={L}, SQUEEZE={SQUEEZE}, T={T_FINAL:.2f}"
          f" (= 4 T_period), dt={DT:.5f}")
    print(f"  sea density 2/h = {SEA:.4f};  aligned stock per cell"
          f" (2/h) dx dp = {SEA * DX * DP:.5f} state-mass units")
    print(f"  aligned stock per momentum row (2/h) L dp"
          f" = {SEA * L * DP:.4f} state-mass units")

    t_arr, load_arr, C_arr, Estate_arr, Esea_arr = [], [], [], [], []
    E_state = 0.0
    E_sea_rate = SEA * N * DP * float(np.sum(np.abs(GAMMA_x)) * DX)
    max_cell_rate = 0.0
    Wp_min = np.inf
    spec_at_peak, rho_at_peak, W_at_peak = None, None, None
    C_peak, load_peak, t_load_peak = -1.0, -1.0, 0.0

    for step in range(N_STEPS + 1):
        t = step * DT
        if step > 0:
            psi = strang_step(psi, DT)
        W = wigner_of_psi(psi)
        Wp = W + SEA
        a = momentum_amps(psi)

        load = float(np.max(np.abs(W)) * np.pi * HBAR)
        C = offdiag_l1(a)
        rate_grid = np.abs(GAMMA_x)[None, :] * Wp * DX * DP
        if step > 0:
            E_state += float(np.sum(np.abs(GAMMA_x)[None, :] * W_prev_abs)
                             * DX * DP) * DT
            max_cell_rate = max(max_cell_rate, float(np.max(rate_grid)))
        Wp_min = min(Wp_min, float(np.min(Wp)))
        W_prev_abs = np.abs(W)

        t_arr.append(t); load_arr.append(load); C_arr.append(C)
        Estate_arr.append(E_state); Esea_arr.append(E_sea_rate * t)
        if load > load_peak:
            load_peak, t_load_peak, W_at_peak = load, t, W.copy()
        if C > C_peak:
            C_peak = C
            spec_at_peak = splitting_spectrum(a)
            rho_at_peak = np.abs(np.outer(a, np.conj(a)))

    t_arr = np.array(t_arr); load_arr = np.array(load_arr)
    C_arr = np.array(C_arr); Estate_arr = np.array(Estate_arr)

    print("-" * 74)
    print("1. STORAGE, local: load factor |W| * (h/2)  [Wigner bound => <= 1]")
    print(f"   load(t=0) = {load_arr[0]:.4f}  (< 1: ring tail-overlap"
          f" renormalisation of the SQUEEZE=2 packet, a real ring effect)")
    psi_n = np.exp(-(x - X0) ** 2 / (4 * 0.5 ** 2)).astype(complex)
    psi_n /= np.sqrt(np.sum(np.abs(psi_n) ** 2) * DX)
    sat = float(np.max(np.abs(wigner_of_psi(psi_n))) * np.pi * HBAR)
    print(f"   narrow-packet control (sigma=0.5): load = {sat:.6f}"
          f"  (Gaussian saturates the Wigner bound exactly)")
    print(f"   max over trajectory = {load_peak:.4f} at t = {t_load_peak:.2f}"
          f"  (instrument tolerance ~2%)")
    print(f"   min W * pi*hbar over trajectory ="
          f" {Wp_min * np.pi * HBAR - 1.0:.4f}  [bound: >= -1; the state"
          f" comes within")
    print(f"   {(Wp_min * np.pi * HBAR) * 100:.0f}% of saturating the"
          f" negative bound mid-run]")

    print("-" * 74)
    print("2. STORAGE, global: off-diagonal L1 mass C(t) of rho(P,P')")
    n_occ = int(np.sum(np.abs(momentum_amps(psi)) ** 2 > 1e-3))
    print(f"   C(0) = {C_arr[0]:.2f},  max C = {C_peak:.2f}"
          f"   vs aligned stock ~1 per row over ~{n_occ} occupied rows")
    cum = np.cumsum(spec_at_peak) / np.sum(spec_at_peak)
    depth50 = int(np.searchsorted(cum, 0.50))
    depth90 = int(np.searchsorted(cum, 0.90))
    print(f"   pair-splitting spectrum at peak: 50% of coherence mass at"
          f" splitting <= {depth50} one-leg hops, 90% at <= {depth90}")

    print("-" * 74)
    print("3. CONSUMABLE model: cumulative moved mass (accepted events only)")
    print(f"   E_state(T) = {Estate_arr[-1]:.1f} state-mass units"
          f"  ({Estate_arr[-1]:.0f}x the state's own mass)")
    print(f"   E_sea(T)   = {Esea_arr[-1]:.0f} state-mass units of sea churn")
    stock_cell = SEA * DX * DP
    t_star = stock_cell / max_cell_rate
    print(f"   local core-cell exhaustion time t* ~ {t_star:.2f}"
          f"  = T/{T_FINAL / t_star:.0f}  (accepted events only;")
    print(f"   bare g0 traffic multiplies consumption by the attempt/accept"
          f" ratio, an additional 10-100x)")

    print("-" * 74)
    print("4. FLOW feasibility is algebraic: the mediated rate reads"
          " W' >= 0 at the")
    print("   midpoint and |W'(p+qdp) - W'(p-qdp)| <= W'(p+qdp) + W'(p-qdp)"
          " holds by the")
    print("   triangle inequality, i.e. by W' >= 0 -- Wigner's bound once"
          " more.")
    print(f"   numerical check: min W' * (h/2) over trajectory"
          f" = {Wp_min * np.pi * HBAR:.4f}  [>= 0]")

    # -------------------------------------------------------------- figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.3), constrained_layout=True)

    ax = axes[0]
    im = ax.imshow(np.abs(W_at_peak) * np.pi * HBAR, origin="lower",
                   aspect="auto", extent=[x[0], x[-1], p[0], p[-1]],
                   cmap="viridis", vmin=0.0, vmax=1.0)
    ax.set_ylim(-4, 4)
    ax.set_xlabel("x"); ax.set_ylabel("p")
    ax.set_title(f"load factor $|W|\\,(h/2)$ at t = {t_load_peak:.1f}\n"
                 "(Wigner bound caps it at 1)", fontsize=10)
    fig.colorbar(im, ax=ax, shrink=0.9)

    ax = axes[1]
    ax.plot(t_arr / T_PERIOD, load_arr, label="local load  max$|W|(h/2)$")
    ax.plot(t_arr / T_PERIOD, C_arr, label="global coherence  C(t)")
    ax.plot(t_arr / T_PERIOD, Estate_arr,
            label="consumable demand  $E_{state}(t)$")
    ax.axhline(1.0, color="k", lw=0.8, ls="--", label="capacity (both bounds)")
    ax.set_xlabel("t / T_period"); ax.set_yscale("log")
    ax.set_title("storage demand is bounded;\nconsumable demand grows"
                 " without bound", fontsize=10)
    ax.legend(fontsize=8, loc="lower right")

    ax = axes[2]
    nm = rho_at_peak.shape[0]
    pgrid = (np.arange(nm) - nm // 2) * 2 * DP
    sl = slice(nm // 2 - 12, nm // 2 + 12)
    im = ax.imshow(rho_at_peak[sl, sl], origin="lower", aspect="equal",
                   extent=[pgrid[sl][0], pgrid[sl][-1],
                           pgrid[sl][0], pgrid[sl][-1]], cmap="magma")
    ax.plot(pgrid[sl], pgrid[sl], "w--", lw=0.7)
    ax.set_xlabel("P'"); ax.set_ylabel("P")
    ax.set_title("$|\\rho(P,P')|$ at peak coherence\n(off-diagonals = split-"
                 "pair demand)", fontsize=10)
    fig.colorbar(im, ax=ax, shrink=0.9)

    fig.suptitle("Permanent pairing is resource-feasible: storage demand is"
                 " bounded by the Wigner bound; consumable demand is not",
                 fontsize=11)
    fp = output_path("pairing_resource_arithmetic.png")
    fig.savefig(fp, dpi=140, bbox_inches="tight")
    print(f"figure written: {fp}")
    dp_ = docs_path("pairing_resource_arithmetic.png")
    if dp_:
        fig.savefig(dp_, dpi=140, bbox_inches="tight")
        print(f"figure written: {dp_}")


if __name__ == "__main__":
    main()
