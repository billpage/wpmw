"""
Demonstration: sea-dressed two-body microdynamics of the crystal lattice.

Companion to ``docs/analysis/sea_dressed_microdynamics.md``, which is in
turn the sequel to ``docs/analysis/four_rule_microdynamics_equivalence.md``.

The four-rule scheme still consults the mean-field mode amplitude
Gamma_q(x) and non-participant occupancies in its rate laws. This demo
exercises the next descent: sixteen collision channels, each a local,
momentum-conserving, two-body event between world-particles (unpaired
positons U+, unpaired negatons U-, ground sea pairs S, and polarized sea
excitations as a pinned reservoir), whose mean field reproduces the QLE
collision term exactly at pinned sea.

  A. **Generator identity.** The sixteen-channel mean-field generator
     (assembled channel by channel, crossing conjugates included, on
     *independent* random fields U+ and U-) is compared against the
     original single-rule stencil acting on E = U+ - U-. Agreement to
     machine precision verifies the channel table's bookkeeping,
     including the crossing structure.

  B. **Stochastic evolution with a live sea.** The squeezed Gaussian in
     the cosine well (same setup as the four-rule demo) is evolved with:
       1. deterministic mesh QLE (target),
       2. four-rule MC (reference noise floor),
       3. sea-dressed MC, pinned sea (level-1 reservoir idealization),
       4. sea-dressed MC, live ledger, kappa_rec = 0   (no recombination),
       5. sea-dressed MC, live ledger, kappa_rec small,
       6. sea-dressed MC, live ledger, kappa_rec large.
     Measured: relative L2 deviation from the mesh, worst-cell sea
     depletion min S/B, and the orphan load (U+ + U-). Prediction: with
     kappa_rec = 0 the sea drains and the orphan load inflates (every
     capture/split/emission event orphans a partner), degrading the
     dynamics; with increasing kappa_rec the live-sea run converges to
     the pinned run. Worldline invariants (U+ + S and U- + S totals)
     are asserted constant across all jump and recombination events.
"""

from __future__ import annotations

import sys
import time

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from wpmwlib.phase_space_crystal_lattice import (  # noqa: E402
    FourierMode,
    PhaseSpaceCrystalLattice,
)
from wpmwlib.sea_dressed_lattice import SeaDressedLattice  # noqa: E402
from wpmwlib.wpmw_utils import output_path, docs_path  # noqa: E402

# --------------------------------------------------------------------------
# Shared parameters (identical to demo_four_rule_equivalence.py Part B)
# --------------------------------------------------------------------------
HBAR = 1.0
MASS = 1.0
L = 8.0
V_P = 1.5
PHI = np.pi
M_CELLS = 64
N_CELLS = 64

OMEGA = (2.0 * np.pi / L) * np.sqrt(V_P / MASS)
T_PERIOD = 2.0 * np.pi / OMEGA
SIGMA_X_GS = np.sqrt(HBAR / (2.0 * MASS * OMEGA))
SIGMA_X0 = 1.3 * SIGMA_X_GS
SIGMA_P0 = HBAR / (2.0 * SIGMA_X0)

NU = 1600000
SUBSTEPS = 16
T_FINAL = 2.0 * T_PERIOD
SEED = 20260706

MODES = [FourierMode(q=1, V_q=V_P, phi_q=PHI)]
KAPPAS = (0.0, 20.0, 200.0)          # live-sea recombination rates


def W_initial(X: np.ndarray, P: np.ndarray) -> np.ndarray:
    return (1.0 / (np.pi * HBAR)) * np.exp(
        -(X ** 2) / (2.0 * SIGMA_X0 ** 2) - (P ** 2) / (2.0 * SIGMA_P0 ** 2)
    )


# --------------------------------------------------------------------------
# Part A: sixteen-channel generator identity on independent U+, U- fields
# --------------------------------------------------------------------------
def random_field(rng: np.random.Generator, N: int, M: int) -> np.ndarray:
    kmax = 6
    F = np.zeros((N, M), dtype=complex)
    for kn in range(-kmax, kmax + 1):
        for km in range(-kmax, kmax + 1):
            F[kn % N, km % M] = rng.standard_normal() + 1j * rng.standard_normal()
    W = np.fft.ifft2(F).real
    return W / np.abs(W).max()


def part_a() -> float:
    rng = np.random.default_rng(SEED)
    mode_sets = {
        "q=1 (canonical)": [FourierMode(1, 1.5, np.pi)],
        "q=2": [FourierMode(2, 0.7, 0.3)],
        "q=1,2,3 mixed": [
            FourierMode(1, 1.5, np.pi),
            FourierMode(2, 0.7, 0.3),
            FourierMode(3, 0.4, 1.1),
        ],
    }
    dt = 0.01
    worst = 0.0
    print("Part A: sixteen-channel generator vs single-rule QLE stencil")
    print("  (independent random U+ and U- fields; exercises crossing structure)")
    print(f"  {'mode set':<18} {'trial':>5} {'max |dE16 - dE1|':>18} {'rel':>12}")
    lat = PhaseSpaceCrystalLattice(M_CELLS, N_CELLS, L, MASS, HBAR)
    for name, modes in mode_sets.items():
        for trial in range(3):
            up = np.abs(random_field(rng, N_CELLS, M_CELLS)) + 0.2
            um = np.abs(random_field(rng, N_CELLS, M_CELLS)) + 0.2
            E0 = up - um
            # Reference: original stencil on E
            ref = PhaseSpaceCrystalLattice(M_CELLS, N_CELLS, L, MASS, HBAR)
            ref.W = E0.copy()
            ref.step_jump_fourier(modes, dt)
            dE_ref = ref.W - E0
            # Sixteen-channel generator, assembled channel by channel
            dE_16 = np.zeros_like(E0)
            for mode in modes:
                Gamma = -(mode.V_q / HBAR) * np.sin(
                    2.0 * np.pi * mode.q * lat.X / L + mode.phi_q
                )
                dE_16 += SeaDressedLattice.channel_generator_mesh(
                    up, um, Gamma, mode.q
                )
            dE_16 *= dt
            err = np.abs(dE_16 - dE_ref).max()
            scale = np.abs(dE_ref).max()
            rel = err / scale if scale > 0 else 0.0
            worst = max(worst, err)
            print(f"  {name:<18} {trial:>5} {err:>18.3e} {rel:>12.3e}")
    print(f"  worst absolute deviation: {worst:.3e}")
    assert worst < 1e-12, "sixteen-channel generator deviates from the QLE stencil"
    return worst


# --------------------------------------------------------------------------
# Part B: stochastic evolution with a live sea
# --------------------------------------------------------------------------
def evolve() -> dict:
    dx = L / M_CELLS
    dp = np.pi * HBAR / L
    dt_adv = MASS * dx / dp
    n_macro = int(round(T_FINAL / dt_adv))
    dt_jump = dt_adv / SUBSTEPS
    gmax = V_P / HBAR
    print("\nPart B: cosine-well evolution with a live sea")
    print(f"  macro steps = {n_macro}, dt_adv = {dt_adv:.4f}, "
          f"substeps = {SUBSTEPS}, |Gamma|max dt_jump = {gmax * dt_jump:.4f}")

    mesh = PhaseSpaceCrystalLattice(M_CELLS, N_CELLS, L, MASS, HBAR, nu=None)
    mc4 = PhaseSpaceCrystalLattice(M_CELLS, N_CELLS, L, MASS, HBAR, nu=NU)
    pinned = SeaDressedLattice(M_CELLS, N_CELLS, L, MASS, HBAR, nu=NU, pinned=True)
    live = {k: SeaDressedLattice(M_CELLS, N_CELLS, L, MASS, HBAR, nu=NU,
                                 pinned=False) for k in KAPPAS}
    solvers = {"mesh": mesh, "mc4": mc4, "pinned": pinned}
    solvers.update({f"live_k{k:g}": live[k] for k in KAPPAS})
    for s in solvers.values():
        s.initialize_from_wigner(W_initial)

    rngs = {name: np.random.default_rng(SEED + 10 + i)
            for i, name in enumerate(solvers)}
    inv0 = {name: s.worldline_invariants()
            for name, s in solvers.items()
            if isinstance(s, SeaDressedLattice) and not s.pinned}
    print(f"  sea background B = {pinned.B} pairs/cell")

    keys = list(solvers.keys())
    sea_keys = [k for k in keys if isinstance(solvers[k], SeaDressedLattice)]
    live_keys = [k for k in keys if k.startswith("live")]
    snap_steps = sorted({0, n_macro // 2, n_macro})
    snaps = {k: {} for k in ("mesh", "pinned", f"live_k{KAPPAS[0]:g}",
                             f"live_k{KAPPAS[-1]:g}")}
    times = []
    err = {k: [] for k in keys if k != "mesh"}
    sea_min = {k: [] for k in live_keys}
    orphans = {k: [] for k in sea_keys}

    def record(step: int):
        t = step * dt_adv
        times.append(t)
        Wm = mesh.get_wigner()
        norm = np.sqrt((Wm ** 2).sum())
        for k in err:
            Wk = solvers[k].get_wigner()
            err[k].append(np.sqrt(((Wk - Wm) ** 2).sum()) / norm)
        for k in live_keys:
            sea_min[k].append(solvers[k].sea_min_fraction())
        for k in sea_keys:
            orphans[k].append(solvers[k].unpaired_total())
        if step in snap_steps:
            for k in snaps:
                snaps[k][t] = solvers[k].get_wigner().copy()

    record(0)
    t0 = time.time()
    for step in range(1, n_macro + 1):
        for s in solvers.values():
            s.step_advect(dt_adv)
        for _ in range(SUBSTEPS):
            mesh.step_jump_fourier(MODES, dt_jump)
            mc4.step_jump_four_rule_mc(MODES, dt_jump, rng=rngs["mc4"])
            pinned.step_jump_sea_mc(MODES, dt_jump, rng=rngs["pinned"])
            for k in KAPPAS:
                name = f"live_k{k:g}"
                live[k].step_jump_sea_mc(MODES, dt_jump, rng=rngs[name])
                live[k].step_recombine(dt_jump, k, rng=rngs[name])
        record(step)
        if step % 10 == 0 or step == n_macro:
            msg = "  ".join(f"{k}={err[k][-1]:.3f}" for k in err)
            print(f"  step {step:3d}/{n_macro}  t={step * dt_adv:6.2f}  "
                  f"relL2: {msg}  [{time.time() - t0:5.1f} s]")

    # Worldline-continuity assertion: no fundamental particle created/destroyed
    for name in inv0:
        assert solvers[name].worldline_invariants() == inv0[name], (
            f"worldline invariant violated in {name}")
    print("  worldline invariants (U+ + S, U- + S): conserved exactly "
          "in all sea runs")

    return dict(snaps=snaps, times=np.array(times), err=err,
                sea_min=sea_min, orphans=orphans, lattice=mesh)


def make_figures(res: dict, worst_a: float) -> None:
    lat = res["lattice"]
    extent = [lat.x[0], lat.x[-1] + lat.dx, lat.p[0], lat.p[-1] + lat.dp]
    snap_times = sorted(res["snaps"]["mesh"].keys())
    rows = [("mesh QLE", "mesh"), ("sea MC (pinned)", "pinned"),
            (f"sea MC (live, $\\kappa$={KAPPAS[0]:g})", f"live_k{KAPPAS[0]:g}"),
            (f"sea MC (live, $\\kappa$={KAPPAS[-1]:g})", f"live_k{KAPPAS[-1]:g}")]
    vmax = max(np.abs(W).max() for W in res["snaps"]["mesh"].values())

    fig, axes = plt.subplots(len(rows), len(snap_times),
                             figsize=(3.4 * len(snap_times), 3.0 * len(rows)),
                             sharex=True, sharey=True)
    for col, t in enumerate(snap_times):
        for r, (label, key) in enumerate(rows):
            ax = axes[r, col]
            im = ax.imshow(res["snaps"][key][t], origin="lower", extent=extent,
                           aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
            if r == 0:
                ax.set_title(f"t = {t:.2f}  ({t / T_PERIOD:.2f} T)")
            if col == 0:
                ax.set_ylabel(f"{label}\np")
            if r == len(rows) - 1:
                ax.set_xlabel("x")
    fig.colorbar(im, ax=axes, shrink=0.85, label="W(x, p)")
    fig.suptitle(
        "Cosine well: mesh QLE vs sea-dressed MC, pinned and live sea "
        f"(nu = {NU}; generator identity max err = {worst_a:.1e})")
    fname = "sea_dressed_evolution.png"
    fig.savefig(output_path(fname), dpi=150, bbox_inches="tight")
    dpth = docs_path(fname)
    if dpth:
        fig.savefig(dpth, dpi=150, bbox_inches="tight")
    plt.close(fig)

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 4.2))
    tt = res["times"] / T_PERIOD
    style = {"mc4": dict(color="0.6", ls="--", label="four-rule MC (ref)"),
             "pinned": dict(color="C0", label="sea MC, pinned"),
             f"live_k{KAPPAS[0]:g}": dict(color="C3",
                                          label=f"live, $\\kappa$={KAPPAS[0]:g}"),
             f"live_k{KAPPAS[1]:g}": dict(color="C1",
                                          label=f"live, $\\kappa$={KAPPAS[1]:g}"),
             f"live_k{KAPPAS[2]:g}": dict(color="C2",
                                          label=f"live, $\\kappa$={KAPPAS[2]:g}")}
    for k, st in style.items():
        ax1.plot(tt, res["err"][k], **st)
    ax1.set_xlabel("t / T_period")
    ax1.set_ylabel("relative L2 deviation from mesh QLE")
    ax1.set_title("Fidelity")
    ax1.legend(fontsize=8)
    ax1.grid(alpha=0.3)
    for k in res["sea_min"]:
        ax2.plot(tt, res["sea_min"][k], **style[k])
    ax2.set_xlabel("t / T_period")
    ax2.set_ylabel("min over cells of S / B")
    ax2.set_title("Sea depletion (live runs)")
    ax2.grid(alpha=0.3)
    for k in res["orphans"]:
        ax3.plot(tt, np.array(res["orphans"][k]) / 1e6, **style[k])
    ax3.set_xlabel("t / T_period")
    ax3.set_ylabel("orphan load (U+ + U-) / 1e6")
    ax3.set_title("Unpaired population")
    ax3.grid(alpha=0.3)
    fig.tight_layout()
    fname = "sea_dressed_metrics.png"
    fig.savefig(output_path(fname), dpi=150, bbox_inches="tight")
    dpth = docs_path(fname)
    if dpth:
        fig.savefig(dpth, dpi=150, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    worst_a = part_a()
    res = evolve()
    make_figures(res, worst_a)
    print("\nFinal relative L2 deviation from mesh QLE:")
    for k, v in res["err"].items():
        print(f"  {k:<12}: {v[-1]:.4f}")
    for k in res["sea_min"]:
        print(f"  final min S/B ({k}): {res['sea_min'][k][-1]:.3f}")
    print("Figures written: sea_dressed_evolution.png, sea_dressed_metrics.png")
    return 0


if __name__ == "__main__":
    sys.exit(main())
