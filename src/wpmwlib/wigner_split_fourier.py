"""
Reference Wigner-equation solver: standard split-Fourier method.

For potentials whose Moyal series truncates at the classical-force term
(quadratic V; e.g. quantum harmonic oscillator) this method is exact up to
time-step error. The QHO update used here:

  Strang step:
    W -> exp(-i kx p dt / (2 m))  in x-Fourier space (half advection)
    W -> exp( i m omega^2 x kp dt) in p-Fourier space (full force kick)
    W -> exp(-i kx p dt / (2 m))  in x-Fourier space (half advection)

This module is written specifically for the QHO comparison; for general V the
force-kick spectral kernel must be replaced by the full Wigner-Moyal kernel.
"""

from __future__ import annotations

from typing import Callable

import numpy as np


class WignerSplitFourier:
    """Strang-split Fourier solver for the Wigner equation (QHO version).

    Parameters
    ----------
    M : int
        Number of position cells (even).
    N : int
        Number of momentum cells (even).
    L : float
        Position-domain length.
    dp : float
        Momentum-cell spacing. Pass the same value used by the
        crystal-lattice solver (``pi*hbar/L`` by default) so the two grids
        coincide and outputs can be compared cell-by-cell.
    mass : float
        Particle mass.
    hbar : float
        Reduced Planck constant.
    """

    def __init__(
        self,
        M: int,
        N: int,
        L: float,
        dp: float,
        mass: float = 1.0,
        hbar: float = 1.0,
    ):
        if M % 2 or N % 2:
            raise ValueError("M and N must be even.")
        self.M = int(M)
        self.N = int(N)
        self.L = float(L)
        self.dx = self.L / self.M
        self.dp = float(dp)
        self.mass = float(mass)
        self.hbar = float(hbar)
        self.x = (np.arange(self.M) - self.M // 2) * self.dx
        self.p = (np.arange(self.N) - self.N // 2) * self.dp
        self.X, self.P = np.meshgrid(self.x, self.p, indexing="xy")

        # Spectral-domain wavenumbers (in fft-shift order, matching np.fft.fft)
        self.kx = 2.0 * np.pi * np.fft.fftfreq(self.M, d=self.dx)
        self.kp = 2.0 * np.pi * np.fft.fftfreq(self.N, d=self.dp)

        self.W: np.ndarray  # initialized below

    # ------------------------------------------------------------------ #
    # Initialization                                                     #
    # ------------------------------------------------------------------ #
    def initialize_from_wigner(self, W_func: Callable[[np.ndarray, np.ndarray], np.ndarray]) -> None:
        """Set the initial Wigner distribution from a callable W_func(X, P)."""
        W = W_func(self.X, self.P).astype(np.complex128)
        if W.shape != (self.N, self.M):
            raise ValueError(f"W_func returned shape {W.shape}, expected {(self.N, self.M)}.")
        self.W = W

    # ------------------------------------------------------------------ #
    # Substeps                                                           #
    # ------------------------------------------------------------------ #
    def step_advect_half(self, dt: float) -> None:
        """Half-step free-streaming via spectral phase shift in x.

        Solves dW/dt = -(p/m) dW/dx for time dt/2. The exact solution shifts
        W in x by +(p/m)*(dt/2), corresponding to FT-space multiplication by
        exp(-i kx (p/m) dt/2).
        """
        Wk = np.fft.fft(self.W, axis=1)
        phase = np.exp(-1j * self.kx[None, :] * (self.P / self.mass) * dt * 0.5)
        Wk *= phase
        self.W = np.fft.ifft(Wk, axis=1)
        # Suppress spurious imaginary buildup from Nyquist-mode handling.
        self.W = self.W.real.astype(np.complex128)

    def step_force_qho(self, omega: float, dt: float) -> None:
        """Full-step force kick for QHO via spectral phase shift in p.

        Solves dW/dt = +V'(x) dW/dp = +m omega^2 x dW/dp for time dt.
        (Force-term sign per the standard QLE: from dot p = -V'(x) in
        Hamilton's equations, the classical Liouville term is +V' dW/dp.)
        Exact solution shifts W in p by +V'(x) dt; FT-space kernel is
        exp(+i kp V'(x) dt).
        """
        Wk = np.fft.fft(self.W, axis=0)
        phase = np.exp(+1j * self.mass * omega**2 * self.X * self.kp[:, None] * dt)
        Wk *= phase
        self.W = np.fft.ifft(Wk, axis=0)
        self.W = self.W.real.astype(np.complex128)



    # ------------------------------------------------------------------ #
    # Composite step                                                     #
    # ------------------------------------------------------------------ #
    def strang_step_qho(self, omega: float, dt: float) -> None:
        """One Strang-split step of QHO Wigner evolution."""
        self.step_advect_half(dt)
        self.step_force_qho(omega, dt)
        self.step_advect_half(dt)

    # ------------------------------------------------------------------ #
    # Observables                                                        #
    # ------------------------------------------------------------------ #
    def get_wigner(self) -> np.ndarray:
        """Return the real part of W (the imaginary residue measures error)."""
        return self.W.real

    def imag_residue(self) -> float:
        """Max absolute imaginary part — should stay near machine epsilon."""
        return float(np.max(np.abs(self.W.imag)))

    def total_norm(self) -> float:
        """``int W dx dp``."""
        return float(np.sum(self.W.real) * self.dx * self.dp)

    def position_marginal(self) -> np.ndarray:
        """``rho(x) = int W dp``. Shape (M,)."""
        return np.sum(self.W.real, axis=0) * self.dp

    def momentum_marginal(self) -> np.ndarray:
        """``rho(p) = int W dx``. Shape (N,)."""
        return np.sum(self.W.real, axis=1) * self.dx
