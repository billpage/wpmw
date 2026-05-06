"""
Regression test for the sign-convention correction in §6.3 of
``docs/supplement/phase_space_crystal_lattice_supplement.md``.

Setup
-----
Cosine hill centred at the origin: V(x) = +V_max cos(2 pi x / L), so V has a
maximum at x = 0 and a minimum at x = +-L/2. Place a Gaussian Wigner state
localised at (x_0 = +L/4, p = 0) — the downhill slope. Classical mechanics:

    V'(x_0) = -V_max (2 pi/L) sin(pi/2) = -V_max (2 pi/L) < 0
    F(x_0)  = -V'(x_0) > 0,   so a particle there should accelerate to MORE
                              POSITIVE p.

Apply each candidate discrete update rule for 50 steps with no free-streaming
and report the final < p >.

Candidates
----------
A) V2 PDF general formula (correct):
       rate   = -(V_max/hbar) sin(2 pi x/L)
       bracket= W(p + dp) - W(p - dp)
B) V2 PDF simplified formula and python code (rate sign flipped — the bug):
       rate   = +(V_max/hbar) sin(2 pi x/L)
       bracket= W(p + dp) - W(p - dp)
C) Original spec §3c line 85 (rate correct, bracket sign flipped — the bug):
       rate   = -(V_max/hbar) sin(2 pi x/L)
       bracket= W(p - dp) - W(p + dp)
D) Corrected spec §3c (= corrected algorithm module step_jump_fourier):
       rate   = -(V_max/hbar) sin(2 pi x/L)
       bracket= W(p + dp) - W(p - dp)        [== A]

Pass criteria
-------------
- A  drives <p> > 0  (correct: downhill).
- B  drives <p> < 0  (the bug: uphill).
- C  drives <p> < 0  (the bug, in algebraically equivalent form).
- D  drives <p> > 0  AND matches A to machine epsilon (verifies the algorithm
   module's step_jump_fourier is on the corrected branch).
- |B - C| ~ 0  (B and C are the same rule rearranged).
"""

from __future__ import annotations
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from wpmwlib.wpmw_utils import output_path  # noqa: E402
from wpmwlib.phase_space_crystal_lattice import PhaseSpaceCrystalLattice, FourierMode  # noqa: E402

# ---------------- grid and initial state --------------------------------- #
HBAR = 1.0
L    = 8.0
M    = 64
N    = 64
DX   = L / M
DP   = np.pi * HBAR / L
x    = (np.arange(M) - M // 2) * DX
p    = (np.arange(N) - N // 2) * DP
X, P = np.meshgrid(x, p, indexing="xy")

x0    = L / 4.0
sigma = 0.5
V_MAX = 1.0
DT    = 0.01
NSTEPS = 50


def gaussian_wigner(X, P):
    g = np.exp(-((X - x0) ** 2 + P ** 2) / (2 * sigma ** 2))
    return g / (2.0 * np.pi * sigma ** 2)


def centroid_p(W):
    return float(np.sum(P * W) * DX * DP / (np.sum(W) * DX * DP))


# ---------------- the four candidate updates ----------------------------- #
def step_A_general(W, dt):
    """V2 PDF general formula. rate = -(V/hbar) sin(theta), bracket = W_hi - W_lo."""
    rate = -(V_MAX / HBAR) * np.sin(2 * np.pi * X / L)
    Whi  = np.roll(W, -1, axis=0)
    Wlo  = np.roll(W, +1, axis=0)
    return W + dt * rate * (Whi - Wlo)


def step_B_simplified(W, dt):
    """V2 PDF simplified / Python: rate sign flipped relative to A (the bug)."""
    rate = +(V_MAX / HBAR) * np.sin(2 * np.pi * X / L)
    Whi  = np.roll(W, -1, axis=0)
    Wlo  = np.roll(W, +1, axis=0)
    return W + dt * rate * (Whi - Wlo)


def step_C_spec_original(W, dt):
    """Original spec §3c line 85: rate correct, bracket flipped (the bug, rearranged)."""
    Gamma = -(V_MAX / HBAR) * np.sin(2 * np.pi * X / L)
    Whi   = np.roll(W, -1, axis=0)
    Wlo   = np.roll(W, +1, axis=0)
    return W + dt * Gamma * (Wlo - Whi)


def step_D_module(W_init, dt, nsteps):
    """The actual algorithm module's step_jump_fourier — should equal A."""
    psc = PhaseSpaceCrystalLattice(M=M, N=N, L=L, mass=1.0, hbar=HBAR)
    psc.W = W_init.copy()
    mode = FourierMode(q=1, V_q=V_MAX, phi_q=0.0)
    for _ in range(nsteps):
        psc.step_jump_fourier([mode], dt)
    return psc.get_wigner()


# ---------------- run all four ------------------------------------------ #
def run_local(stepper, label):
    W = gaussian_wigner(X, P)
    for _ in range(NSTEPS):
        W = stepper(W, DT)
    return W, centroid_p(W)


def main():
    print("=" * 78)
    print("Sign-convention regression test")
    print(f"  cosine hill V = +V_max cos(2 pi x/L), V_max = {V_MAX}, hbar = 1")
    print(f"  initial state at (x0 = +L/4 = +{x0}, p = 0)")
    print(f"  classical force at x0: F = -V'(x0) > 0  ->  expect <p> to drift POSITIVE")
    print(f"  grid: M = {M}, N = {N}, dx = {DX}, dp = {DP:.5f}, dt = {DT}, "
          f"nsteps = {NSTEPS}")
    print("=" * 78)

    Wa, pA = run_local(step_A_general,        "A) V2 general formula")
    Wb, pB = run_local(step_B_simplified,     "B) V2 simplified / Python (bug)")
    Wc, pC = run_local(step_C_spec_original,  "C) original spec §3c (bug)")
    Wd     = step_D_module(gaussian_wigner(X, P), DT, NSTEPS)
    pD     = centroid_p(Wd)
    matchAD = float(np.max(np.abs(Wa - Wd)))
    diffBC  = float(np.max(np.abs(Wb - Wc)))

    print(f"  A) V2 general formula           <p> = {pA:+.4f}")
    print(f"  B) V2 simplified / Python (bug) <p> = {pB:+.4f}")
    print(f"  C) original spec §3c     (bug)  <p> = {pC:+.4f}")
    print(f"  D) algorithm module step_jump_fourier  <p> = {pD:+.4f}")
    print()
    print(f"  max |W_A - W_D| = {matchAD:.2e}    (should be ~0; module on corrected branch)")
    print(f"  max |W_B - W_C| = {diffBC:.2e}     (should be ~0; B and C are the same bug)")

    failures = []
    if not (pA > 1e-3): failures.append("A should drive <p> > 0")
    if not (pB < -1e-3): failures.append("B should drive <p> < 0 (bug expected)")
    if not (pC < -1e-3): failures.append("C should drive <p> < 0 (bug expected)")
    if not (pD > 1e-3): failures.append("D (module) should drive <p> > 0")
    if matchAD > 1e-12: failures.append(f"D should equal A (got max diff {matchAD:.2e})")
    if diffBC  > 1e-12: failures.append(f"B should equal C (got max diff {diffBC:.2e})")

    print()
    if failures:
        print("FAIL:")
        for f in failures:
            print(f"  - {f}")
        sys.exit(1)
    print("PASS — sign-convention correction is in effect.")


if __name__ == "__main__":
    main()

