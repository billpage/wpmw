"""
Phase-Space Crystal-Lattice algorithm for Wigner distribution evolution.

This module implements the algorithm specified in
``docs/algorithm/phase_space_crystal_lattice_algorithm.md``. The
QLE-consistent sign convention adopted in that spec (and in the redrafted
V2 supplement at ``docs/supplement/phase_space_crystal_lattice_supplement.md``)
is the only one used here.

State representation
--------------------
The Wigner distribution ``W(x, p)`` is internally represented either as a
floating-point mesh (deterministic mesh-density form, spec §3c) or as integer
positon counts ``N_plus(m, n)`` of the shifted distribution
``W' = W + 2/h`` (Monte-Carlo particle form, spec §6).

Two potential decompositions are supported
------------------------------------------
1. **Fourier-mode form** (spec §3b/3c). The potential is decomposed as
   ``V(x) = V_0 + sum_q V_q cos(2 pi q x / L + phi_q)``. Each mode drives a
   discrete mediated jump of magnitude ``q`` momentum cells per spec, with the
   natural choice ``dp = pi*hbar/L``. The continuum limit reproduces the QLE
   force term ``+V'(x) dW/dp``.

2. **Differential form** (spec §7b). For non-Fourier-decomposable smooth V,
   apply the QLE jump term ``+V'(x) dW/dp`` via centered finite difference.
   For purely quadratic V (e.g. quantum harmonic oscillator) this is exact in
   the continuum limit because the Moyal series truncates at this term.

Conventions
-----------
- Array shape is ``(N, M)`` with axis 0 = momentum index, axis 1 = position
  index. This matches a ``np.meshgrid(x, p, indexing='xy')`` output.
- Position grid is centered on zero: ``x_m = (m - M/2) * dx``, ``m = 0..M-1``.
- Momentum grid is centered on zero: ``p_n = (n - N/2) * dp``, ``n = 0..N-1``.
- ``dp = pi*hbar/L`` by default (spec choice).
- Both M and N should be even.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Iterable, Optional, Sequence, Tuple

import numpy as np


@dataclass
class FourierMode:
    """A single Fourier component ``V_q cos(2 pi q x / L + phi_q)`` of V(x)."""
    q: int       # mode number (>=1)
    V_q: float   # amplitude
    phi_q: float # phase


class PhaseSpaceCrystalLattice:
    """Phase-space crystal-lattice solver.

    Parameters
    ----------
    M : int
        Number of position cells (even).
    N : int
        Number of momentum cells (even).
    L : float
        Position-domain length. Domain is ``[-L/2, L/2)`` periodic.
    mass : float
        Particle mass.
    hbar : float
        Reduced Planck constant.
    nu : float, optional
        Particles per unit phase-space volume. If None, the deterministic
        mesh-density form is used (default). If a finite number, the integer
        Monte-Carlo form is used.
    advection : {'integer_roll', 'spectral'}
        Free-streaming scheme. 'integer_roll' is the spec §3a scheme (each
        momentum row shifts by an integer number of position cells per step).
        'spectral' uses an exact FFT phase-shift; this is a strictly more
        accurate alternative not in the original spec but useful for
        comparison.
    """

    def __init__(
        self,
        M: int,
        N: int,
        L: float,
        mass: float = 1.0,
        hbar: float = 1.0,
        nu: Optional[float] = None,
        advection: str = "integer_roll",
    ):
        if M % 2 or N % 2:
            raise ValueError("M and N must be even.")
        if advection not in ("integer_roll", "spectral"):
            raise ValueError("advection must be 'integer_roll' or 'spectral'.")
        self.M = int(M)
        self.N = int(N)
        self.L = float(L)
        self.mass = float(mass)
        self.hbar = float(hbar)
        self.h = 2.0 * np.pi * self.hbar
        self.dx = self.L / self.M
        self.dp = np.pi * self.hbar / self.L           # spec §1 choice
        self.x = (np.arange(self.M) - self.M // 2) * self.dx
        self.p = (np.arange(self.N) - self.N // 2) * self.dp
        # Meshgrid with 'xy' indexing -> X.shape == P.shape == (N, M)
        self.X, self.P = np.meshgrid(self.x, self.p, indexing="xy")
        self.W_bg = 2.0 / self.h                       # the crystal shift = 1/(pi*hbar)
        self.nu = nu
        self.advection = advection

        # Precompute spectral wavenumbers along axis 1 (x) for spectral advection
        self._kx = 2.0 * np.pi * np.fft.fftfreq(self.M, d=self.dx)

        # State (allocated by initialize_*)
        self.W: Optional[np.ndarray] = None             # mesh-density form
        self.N_plus: Optional[np.ndarray] = None        # particle form

    # ------------------------------------------------------------------ #
    # Initialization                                                     #
    # ------------------------------------------------------------------ #
    def initialize_from_wigner(self, W_func: Callable[[np.ndarray, np.ndarray], np.ndarray]) -> None:
        """Set the initial Wigner distribution from a callable W_func(X, P)."""
        W = W_func(self.X, self.P).astype(np.float64)
        if W.shape != (self.N, self.M):
            raise ValueError(f"W_func returned shape {W.shape}, expected {(self.N, self.M)}.")
        self.W = W
        if self.nu is not None:
            shifted = (W + self.W_bg) * self.nu * self.dx * self.dp
            self.N_plus = np.round(shifted).astype(np.int64)

    # ------------------------------------------------------------------ #
    # Step 3a: free streaming                                            #
    # ------------------------------------------------------------------ #
    def step_advect(self, dt: float) -> None:
        """Free-streaming step using the configured advection scheme."""
        if self.advection == "integer_roll":
            self._step_advect_integer_roll(dt)
        else:
            self._step_advect_spectral(dt)

    def _step_advect_integer_roll(self, dt: float) -> None:
        """Spec §3a: shift row n by round(p_n * dt / (m * dx)) cells (periodic)."""
        if self.nu is None:
            target = self.W
            new = np.empty_like(target)
        else:
            target = self.N_plus
            new = np.empty_like(target)
        for n in range(self.N):
            shift = int(round(self.p[n] * dt / (self.mass * self.dx)))
            new[n, :] = np.roll(target[n, :], shift)
        if self.nu is None:
            self.W = new
        else:
            self.N_plus = new

    def _step_advect_spectral(self, dt: float) -> None:
        """Exact spectral phase-shift along x. Only valid in mesh-density form."""
        if self.nu is not None:
            raise NotImplementedError("Spectral advection requires mesh-density form (nu=None).")
        Wk = np.fft.fft(self.W, axis=1)
        # Translation by (p/m)*dt in x  <=>  multiply by exp(-i kx (p/m) dt)
        phase = np.exp(-1j * self._kx[None, :] * (self.P / self.mass) * dt)
        Wk *= phase
        self.W = np.fft.ifft(Wk, axis=1).real

    # ------------------------------------------------------------------ #
    # Step 3b: potential-driven mediated jumps                           #
    # ------------------------------------------------------------------ #
    def step_jump_fourier(self, modes: Iterable[FourierMode], dt: float) -> None:
        """Mesh-density Fourier-mode update (spec §3c, QLE-consistent sign).

        Per the corrected spec:

            dW = dt * Gamma_q(x) * (W(p + q dp) - W(p - q dp))
            Gamma_q(x) = -(V_q / hbar) * sin(2 pi q x / L + phi_q)

        In the small-dp continuum limit this reproduces the QLE force term
        ``dW/dt = +V'(x) dW/dp``. See
        ``docs/supplement/phase_space_crystal_lattice_supplement.md`` §6 for
        the derivation, and ``sign_convention_check.py`` for a regression
        test that exercises the difference against the (incorrect) sign that
        appears in the V2 memo's simplified form.
        """
        if self.nu is not None:
            raise RuntimeError(
                "Use step_jump_fourier_mc(...) for the Monte-Carlo particle form."
            )
        dW = np.zeros_like(self.W)
        for mode in modes:
            q, V_q, phi_q = mode.q, mode.V_q, mode.phi_q
            Gamma = -(V_q / self.hbar) * np.sin(
                2.0 * np.pi * q * self.X / self.L + phi_q
            )
            W_lo = np.roll(self.W, +q, axis=0)   # W at p_{n-q}
            W_hi = np.roll(self.W, -q, axis=0)   # W at p_{n+q}
            dW += Gamma * (W_hi - W_lo)
        self.W = self.W + dt * dW

    def step_jump_fourier_mc(
        self,
        modes: Iterable[FourierMode],
        dt: float,
        rng: Optional[np.random.Generator] = None,
    ) -> None:
        """Monte-Carlo particle form of the mediated-jump rule (spec §6).

        Each positon at (m, n) is a candidate mediator. The number of mediator
        events in a cell is sampled from a binomial with probability
        ``min(|Gamma_q(x_m)| dt, 1)``. Each event transfers one positon from
        the source cell to the destination cell, capped by source population.
        """
        if self.nu is None:
            raise RuntimeError("Particle MC requires nu to be set at construction.")
        if rng is None:
            rng = np.random.default_rng()
        N, M = self.N, self.M
        for mode in modes:
            q = mode.q
            Gamma_x = -(mode.V_q / self.hbar) * np.sin(
                2.0 * np.pi * q * self.x / self.L + mode.phi_q
            )
            for m in range(M):
                rate_dt = Gamma_x[m] * dt
                mag = abs(rate_dt)
                if mag == 0.0:
                    continue
                p_jump = min(mag, 1.0)
                sign = 1 if rate_dt > 0 else -1
                col = self.N_plus[:, m].copy()
                events = rng.binomial(col, p_jump)
                # Apply transfers; cap by source population to preserve N+ >= 0.
                # Direction is QLE-consistent (Gamma > 0 transfers (n+q) -> (n-q));
                # see spec §3b and supplement §6.
                src_idx = (np.arange(N) + sign * q) % N
                dst_idx = (np.arange(N) - sign * q) % N
                # Iterate in a fixed order so the source-cap is well-defined.
                for n in range(N):
                    e = events[n]
                    if e == 0:
                        continue
                    src = src_idx[n]
                    dst = dst_idx[n]
                    e = min(int(e), int(col[src]))
                    if e == 0:
                        continue
                    col[src] -= e
                    col[dst] += e
                self.N_plus[:, m] = col

    # ------------------------------------------------------------------ #
    # Step 3b': four-rule two-body jump form (Cyganski, 2026)             #
    # ------------------------------------------------------------------ #
    def step_jump_four_rule(self, modes: Iterable[FourierMode], dt: float) -> None:
        """Mesh-density form of the four-rule (focus/defocus/hop) scheme.

        Implements the *symmetric member* of the exact-rate family derived in
        ``docs/analysis/four_rule_microdynamics_equivalence.md``. Per mode q,
        with Gamma_q(x) = -(V_q/hbar) sin(2 pi q x / L + phi_q):

            focus/defocus channel, signed net rate at center n:
                f_n = (Gamma/2) * (W_{n+q} - W_{n-q})
            hop channel, signed net right-hop rate across n (n-q -> n+q):
                h_n = -(Gamma/2) * (W_{n+q} + W_{n-q})

        Channel effects (per unit net rate):
            focus at n:  W_n += 2,  W_{n-q} -= 1,  W_{n+q} -= 1
            hop across n:  W_{n-q} -= 1,  W_{n+q} += 1

        The channels are assembled *independently* here, precisely so that
        this method is a nontrivial numerical check of the algebraic identity
        (four-rule generator == single-rule generator == QLE stencil). It must
        agree with ``step_jump_fourier`` to machine precision for any W.
        """
        if self.nu is not None:
            raise RuntimeError(
                "Use step_jump_four_rule_mc(...) for the Monte-Carlo particle form."
            )
        dW = np.zeros_like(self.W)
        for mode in modes:
            q, V_q, phi_q = mode.q, mode.V_q, mode.phi_q
            Gamma = -(V_q / self.hbar) * np.sin(
                2.0 * np.pi * q * self.X / self.L + phi_q
            )
            W_hi = np.roll(self.W, -q, axis=0)   # W at p_{n+q}
            W_lo = np.roll(self.W, +q, axis=0)   # W at p_{n-q}
            f = 0.5 * Gamma * (W_hi - W_lo)      # focus/defocus net rate
            h = -0.5 * Gamma * (W_hi + W_lo)     # right-hop net rate
            # Focus channel: center gains 2 f_n; site n is drained as the
            # lower satellite of center n+q and the upper satellite of n-q.
            dW += 2.0 * f
            dW -= np.roll(f, -q, axis=0)         # f_{n+q} drains site n
            dW -= np.roll(f, +q, axis=0)         # f_{n-q} drains site n
            # Hop channel: arrivals from hops across n-q, departures via n+q.
            dW += np.roll(h, +q, axis=0)         # h_{n-q}
            dW -= np.roll(h, -q, axis=0)         # h_{n+q}
        self.W = self.W + dt * dW

    def step_jump_four_rule_mc(
        self,
        modes: Iterable[FourierMode],
        dt: float,
        rng: Optional[np.random.Generator] = None,
    ) -> None:
        """Monte-Carlo particle form of the four-rule scheme.

        Tau-leaping on the integer positon field: per position column, per
        mode, per momentum center n, the signed net channel rates are
        computed from the *excess* counts (background subtracted; the
        background cancels identically in the focus rate and contributes
        only a null uniform bias to the hop rate), and the net number of
        events in [t, t+dt] is drawn as sign(rate) * Poisson(|rate| dt).

        Event application (net event count e, signed):
            focus at n  (e > 0):   N_{n-q} -= e, N_{n+q} -= e, N_n += 2e
            defocus at n (e < 0):  reverse
            right-hop across n (e > 0):  N_{n-q} -= e, N_{n+q} += e
            left-hop (e < 0):            reverse
        Transfers are capped by source populations so N_plus stays >= 0,
        mirroring ``step_jump_fourier_mc``.
        """
        if self.nu is None:
            raise RuntimeError("Particle MC requires nu to be set at construction.")
        if rng is None:
            rng = np.random.default_rng()
        N = self.N
        bg = int(round(self.W_bg * self.nu * self.dx * self.dp))
        for mode in modes:
            q = mode.q
            Gamma_x = -(mode.V_q / self.hbar) * np.sin(
                2.0 * np.pi * q * self.x / self.L + mode.phi_q
            )
            for m in range(self.M):
                g = Gamma_x[m]
                if g == 0.0:
                    continue
                col = self.N_plus[:, m].astype(np.int64)
                excess = col - bg
                ex_hi = np.roll(excess, -q)      # excess at n+q
                ex_lo = np.roll(excess, +q)      # excess at n-q
                F = 0.5 * g * (ex_hi - ex_lo)    # net focus rate at center n
                H = -0.5 * g * (ex_hi + ex_lo)   # net right-hop rate across n
                eF = np.sign(F) * rng.poisson(np.abs(F) * dt)
                eH = np.sign(H) * rng.poisson(np.abs(H) * dt)
                for n in range(N):
                    lo = (n - q) % N
                    hi = (n + q) % N
                    e = int(eF[n])
                    if e > 0:                    # focus: satellites -> center
                        e = min(e, int(col[lo]), int(col[hi]))
                        col[lo] -= e
                        col[hi] -= e
                        col[n] += 2 * e
                    elif e < 0:                  # defocus: center pair splits
                        e = min(-e, int(col[n]) // 2)
                        col[n] -= 2 * e
                        col[lo] += e
                        col[hi] += e
                    e = int(eH[n])
                    if e > 0:                    # right-hop: n-q -> n+q
                        e = min(e, int(col[lo]))
                        col[lo] -= e
                        col[hi] += e
                    elif e < 0:                  # left-hop: n+q -> n-q
                        e = min(-e, int(col[hi]))
                        col[hi] -= e
                        col[lo] += e
                self.N_plus[:, m] = col

    def step_jump_differential(self, dVdx_arr: np.ndarray, dt: float) -> None:
        """Spec §7b: QLE differential form for the jump term.

        Applies an explicit Euler step of
            dW/dt = +V'(x) * dW/dp
        using a centered finite difference for dW/dp on the periodic momentum
        grid. The sign of the force term matches the standard QLE
        ``dW/dt + (p/m)dW/dx - V'(x) dW/dp = 0`` (Hamilton's equation
        ``dot p = -V'`` plus Liouville). For a strictly quadratic V the Moyal
        series truncates at this term so this is exact in the continuum limit.
        """
        if self.nu is not None:
            raise NotImplementedError(
                "Differential jump form is only implemented for mesh-density form."
            )
        if dVdx_arr.shape != (self.M,):
            raise ValueError(
                f"dVdx_arr shape {dVdx_arr.shape}, expected ({self.M},)."
            )
        Wp = np.roll(self.W, -1, axis=0)
        Wm = np.roll(self.W, +1, axis=0)
        dW_dp = (Wp - Wm) / (2.0 * self.dp)
        # +V' dW/dp matches the standard QLE force-term sign.
        self.W = self.W + dt * dVdx_arr[None, :] * dW_dp

    # ------------------------------------------------------------------ #
    # Strang composite step                                              #
    # ------------------------------------------------------------------ #
    def strang_step_fourier(self, modes: Iterable[FourierMode], dt: float) -> None:
        """Strang-split step: half advect, full Fourier jump, half advect."""
        modes_list = list(modes)
        self.step_advect(dt * 0.5)
        self.step_jump_fourier(modes_list, dt)
        self.step_advect(dt * 0.5)

    def strang_step_differential(self, dVdx_arr: np.ndarray, dt: float) -> None:
        """Strang-split step using the differential jump form."""
        self.step_advect(dt * 0.5)
        self.step_jump_differential(dVdx_arr, dt)
        self.step_advect(dt * 0.5)

    # ------------------------------------------------------------------ #
    # Observables                                                        #
    # ------------------------------------------------------------------ #
    def get_wigner(self) -> np.ndarray:
        """Return the (un-shifted) Wigner distribution W."""
        if self.nu is None:
            return self.W
        return self.N_plus / (self.nu * self.dx * self.dp) - self.W_bg

    def total_norm(self) -> float:
        """``int W dx dp``. Should be 1 for a normalized state."""
        return float(np.sum(self.get_wigner()) * self.dx * self.dp)

    def position_marginal(self) -> np.ndarray:
        """``rho(x) = int W dp``. Shape (M,)."""
        return np.sum(self.get_wigner(), axis=0) * self.dp

    def momentum_marginal(self) -> np.ndarray:
        """``rho(p) = int W dx``. Shape (N,)."""
        return np.sum(self.get_wigner(), axis=1) * self.dx
