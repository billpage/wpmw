"""
Sea-dressed two-body microdynamics for the phase-space crystal lattice.

Companion module to ``docs/analysis/sea_dressed_microdynamics.md``. This
module descends one ontological level below the four-rule scheme of
``docs/analysis/four_rule_microdynamics_equivalence.md``: the occupancy
stencils of the four rules are realized as local, momentum-conserving,
two-body collisions between world-particles, with the Dirac sea of
positon-negaton pairs playing the role of the pinned-density reservoir
that the linearity no-go lemma requires.

State (per lattice cell (n, m)):
    U_plus  : count of unpaired ("excess") positons        (>= 0)
    U_minus : count of unpaired negatons (holes)           (>= 0)
    S       : count of ground-state neutral pairs          (>= 0)

The observable excess field is E = U_plus - U_minus, and the Wigner
distribution is W = E / (nu dx dp). The sea background is S_bar = B per
cell, the integer image of the crystal shift 2/h.

Sixteen collision channels (eight per polarization sign of Gamma_q(x))
drive the QLE collision term; see the analysis note for the table, the
collision reading of each channel, and the exactness proof. A separate
recombination channel (U_plus + U_minus -> S at the same cell) restores
the sea; it leaves E invariant and therefore lives entirely outside the
QLE generator.

Conventions follow ``phase_space_crystal_lattice.py``: arrays are (N, M)
with axis 0 = momentum, axis 1 = position; dp = pi*hbar/L.
"""

from __future__ import annotations

from typing import Iterable, Optional

import numpy as np

from .phase_space_crystal_lattice import FourierMode, PhaseSpaceCrystalLattice


class SeaDressedLattice(PhaseSpaceCrystalLattice):
    """Crystal lattice with an explicit positon-negaton sea.

    Parameters
    ----------
    M, N, L, mass, hbar, nu, advection :
        As in :class:`PhaseSpaceCrystalLattice`. ``nu`` must be set (the
        sea-dressed form is intrinsically a particle model).
    sea_factor : float
        Sea density in units of the natural background: the initial pair
        count per cell is ``S = round(sea_factor * nu * (2/h) * dx * dp)``.
        ``sea_factor = 1`` reproduces the crystal shift exactly.
    pinned : bool
        If True, the sea ledger is frozen: rates use S/B = 1 identically
        and no sea bookkeeping is performed (level-1 reservoir idealization,
        the regime of the exactness theorem). If False, every channel
        debits/credits the pair ledger and rates carry the S/B factors
        (level-2, finite back-reaction).
    """

    #: Channel table. Each row: (drain spec, effect spec) in terms of the
    #: roles (hi, n, lo) = (n + sigma q, n, n - sigma q). See analysis note
    #: section "The channel table". Encoded procedurally in _apply_channels.
    CHANNEL_ORDER = ("K1", "K1b", "K2", "K2b", "K3", "K3b", "K4", "K4b")

    def __init__(
        self,
        M: int,
        N: int,
        L: float,
        mass: float = 1.0,
        hbar: float = 1.0,
        nu: Optional[float] = None,
        advection: str = "integer_roll",
        sea_factor: float = 1.0,
        pinned: bool = False,
    ):
        if nu is None:
            raise ValueError("SeaDressedLattice requires nu (particle form).")
        super().__init__(M, N, L, mass=mass, hbar=hbar, nu=nu, advection=advection)
        self.sea_factor = float(sea_factor)
        self.pinned = bool(pinned)
        self.B = int(round(self.sea_factor * self.W_bg * self.nu * self.dx * self.dp))
        if self.B < 1:
            raise ValueError("Sea density B < 1; increase nu or sea_factor.")
        self.U_plus: Optional[np.ndarray] = None
        self.U_minus: Optional[np.ndarray] = None
        self.S: Optional[np.ndarray] = None

    # ------------------------------------------------------------------ #
    # Initialization                                                     #
    # ------------------------------------------------------------------ #
    def initialize_from_wigner(self, W_func) -> None:
        """Excess field from W; virgin sea S = B; no unpaired negatons
        except where W < 0 (there, the excess itself is negatonic)."""
        super().initialize_from_wigner(W_func)  # sets self.W and self.N_plus
        E = np.round(self.W * self.nu * self.dx * self.dp).astype(np.int64)
        self.U_plus = np.maximum(E, 0)
        self.U_minus = np.maximum(-E, 0)
        self.S = np.full((self.N, self.M), self.B, dtype=np.int64)

    # ------------------------------------------------------------------ #
    # Observables (override: excess = U+ - U-)                           #
    # ------------------------------------------------------------------ #
    def get_wigner(self) -> np.ndarray:
        E = self.U_plus - self.U_minus
        return E / (self.nu * self.dx * self.dp)

    def unpaired_total(self) -> int:
        """Total unpaired population U+ + U- (orphan load)."""
        return int(self.U_plus.sum() + self.U_minus.sum())

    def sea_min_fraction(self) -> float:
        """min over cells of S / B (1.0 = undepleted)."""
        return float(self.S.min()) / self.B

    def worldline_invariants(self) -> tuple:
        """(total positons, total negatons) = (U+ + S, U- + S) summed.

        Conserved exactly by every jump channel and by recombination;
        free streaming permutes cells only. A changing value indicates a
        bookkeeping bug (a worldline was created or destroyed).
        """
        s = int(self.S.sum())
        return (int(self.U_plus.sum()) + s, int(self.U_minus.sum()) + s)

    # ------------------------------------------------------------------ #
    # Free streaming: all species ride their momentum row               #
    # ------------------------------------------------------------------ #
    def step_advect(self, dt: float) -> None:
        if self.advection != "integer_roll":
            raise NotImplementedError("Sea-dressed form uses integer_roll advection.")
        for arr in (self.U_plus, self.U_minus, self.S):
            for n in range(self.N):
                shift = int(round(self.p[n] * dt / (self.mass * self.dx)))
                arr[n, :] = np.roll(arr[n, :], shift)

    # ------------------------------------------------------------------ #
    # The sixteen jump channels                                          #
    # ------------------------------------------------------------------ #
    def _gamma_x(self, mode: FourierMode) -> np.ndarray:
        return -(mode.V_q / self.hbar) * np.sin(
            2.0 * np.pi * mode.q * self.x / self.L + mode.phi_q
        )

    def step_jump_sea_mc(
        self,
        modes: Iterable[FourierMode],
        dt: float,
        rng: Optional[np.random.Generator] = None,
    ) -> None:
        """Tau-leap the sixteen collision channels.

        Per mode q and polarization sign sigma = sign(Gamma_q(x)), the
        eight channels of the sigma block are applied lattice-wide in the
        fixed order K1, K1b, K2, K2b, K3, K3b, K4, K4b. Within a single
        channel, each drained cell is drained from exactly one center
        role, so caps are safe under vectorized application; between
        channels, counts are updated sequentially.

        Rates (roles hi = n + sigma q, lo = n - sigma q; gamma = |Gamma|;
        s_c = S_c / B, or 1 if pinned):

            K1  : (gamma/2) s_lo      U+_hi   excess-sea positon capture (focus)
            K1b : (gamma/2) s_lo      U-_hi   crossing conjugate       (defocus)
            K2  : (gamma/2) s_n^2     U+_lo   stimulated sea-sea split (defocus)
            K2b : (gamma/2) s_n^2     U-_lo   crossing conjugate         (focus)
            K3  : (gamma/2)           U+_hi   excitation absorption   (hop hi->lo)
            K3b : (gamma/2)           U-_hi   crossing conjugate      (hop lo->hi)
            K4  : (gamma/2) s_hi      U+_lo   stimulated sea emission (hop hi->lo)
            K4b : (gamma/2) s_hi      U-_lo   crossing conjugate      (hop lo->hi)
        """
        if rng is None:
            rng = np.random.default_rng()
        for mode in modes:
            q = mode.q
            gx = self._gamma_x(mode)
            for sigma in (+1, -1):
                gamma = np.maximum(sigma * gx, 0.0)[None, :]   # (1, M)
                if not np.any(gamma > 0):
                    continue
                self._apply_channel_block(q, sigma, gamma, dt, rng)
        if self.pinned:
            # Level-1 reservoir idealization = pinned sea AND kappa -> inf:
            # co-located unpaired partners re-pair instantly, so the orphan
            # load never feeds back into the gross channel rates. Without
            # this, orphan production (each capture/split/emission event
            # orphans a partner) inflates U+ and U- exponentially at rate
            # ~|Gamma| — zero-mean but variance-exploding noise. With a live
            # ledger the same role is played by step_recombine at finite
            # kappa; see the analysis note, "Recombination is not optional".
            r = np.minimum(self.U_plus, self.U_minus)
            self.U_plus -= r
            self.U_minus -= r

    def _apply_channel_block(self, q, sigma, gamma, dt, rng) -> None:
        # roll(+k, axis=0) brings row (n+k) to row n when used as roll(A, -k)...
        # Convention (matches base class): np.roll(A, -k, axis=0)[n] = A[n+k].
        k = sigma * q
        UP, UM, S = self.U_plus, self.U_minus, self.S
        B = self.B

        def at_hi(A):
            return np.roll(A, -k, axis=0)     # value at n + sigma q, viewed from n

        def to_hi(A):
            return np.roll(A, +k, axis=0)     # scatter: contribution to n + sigma q

        def at_lo(A):
            return np.roll(A, +k, axis=0)

        def to_lo(A):
            return np.roll(A, -k, axis=0)

        def sfac(A_at_role):
            return A_at_role / B if not self.pinned else 1.0

        def draw(rate):
            return rng.poisson(np.maximum(rate, 0.0) * dt)

        # --- K1: excess positon at hi + sea positon at lo -> both at n; hole at lo
        rate = 0.5 * gamma * sfac(at_lo(S)) * at_hi(UP)
        e = draw(rate)
        e = np.minimum(e, at_hi(UP))
        if not self.pinned:
            e = np.minimum(e, at_lo(S))
        UP -= to_hi(e)          # drain U+ at hi
        UP += 2 * e             # two unpaired positons at center
        UM += to_lo(e)          # orphan negaton at lo
        if not self.pinned:
            S -= to_lo(e)

        # --- K1b: crossing conjugate (negaton capture) -> defocus stencil
        rate = 0.5 * gamma * sfac(at_lo(S)) * at_hi(UM)
        e = draw(rate)
        e = np.minimum(e, at_hi(UM))
        if not self.pinned:
            e = np.minimum(e, at_lo(S))
        UM -= to_hi(e)
        UM += 2 * e
        UP += to_lo(e)
        if not self.pinned:
            S -= to_lo(e)

        # --- K2: two sea positons at n scatter to hi/lo (stim. by U+ at lo);
        #         both negatons orphaned at n -> defocus stencil
        sn = sfac(S)
        rate = 0.5 * gamma * (sn * sn if not self.pinned else 1.0) * at_lo(UP)
        e = draw(rate)
        e = np.minimum(e, at_lo(UP))            # stimulator bound: see note §8
        if not self.pinned:
            e = np.minimum(e, S // 2)
        UP += to_hi(e) + to_lo(e)
        UM += 2 * e
        if not self.pinned:
            S -= 2 * e

        # --- K2b: crossing conjugate -> focus stencil
        rate = 0.5 * gamma * (sn * sn if not self.pinned else 1.0) * at_lo(UM)
        e = draw(rate)
        e = np.minimum(e, at_lo(UM))
        if not self.pinned:
            e = np.minimum(e, S // 2)
        UM += to_hi(e) + to_lo(e)
        UP += 2 * e
        if not self.pinned:
            S -= 2 * e

        # --- K3: excitation absorption; unpaired positon hops hi -> lo
        rate = 0.5 * gamma * at_hi(UP)
        e = draw(rate)
        e = np.minimum(e, at_hi(UP))
        UP -= to_hi(e)
        UP += to_lo(e)

        # --- K3b: crossing conjugate; unpaired negaton hops hi -> lo
        rate = 0.5 * gamma * at_hi(UM)
        e = draw(rate)
        e = np.minimum(e, at_hi(UM))
        UM -= to_hi(e)
        UM += to_lo(e)

        # --- K4: destination-stimulated sea emission; sea positon at hi
        #         hops to lo, orphaning its negaton at hi
        rate = 0.5 * gamma * sfac(at_hi(S)) * at_lo(UP)
        e = draw(rate)
        e = np.minimum(e, at_lo(UP))
        if not self.pinned:
            e = np.minimum(e, at_hi(S))
        UP += to_lo(e)
        UM += to_hi(e)
        if not self.pinned:
            S -= to_hi(e)

        # --- K4b: crossing conjugate; sea negaton at hi hops to lo
        rate = 0.5 * gamma * sfac(at_hi(S)) * at_lo(UM)
        e = draw(rate)
        e = np.minimum(e, at_lo(UM))
        if not self.pinned:
            e = np.minimum(e, at_hi(S))
        UM += to_lo(e)
        UP += to_hi(e)
        if not self.pinned:
            S -= to_hi(e)

    # ------------------------------------------------------------------ #
    # Recombination (hidden sea microdynamics)                           #
    # ------------------------------------------------------------------ #
    def step_recombine(
        self,
        dt: float,
        kappa: float,
        rng: Optional[np.random.Generator] = None,
    ) -> None:
        """Pair recombination: U+ + U- -> S at the same cell.

        Mass-action rate kappa * U+ * U- / B per cell. Leaves the excess
        E = U+ - U- (and hence W) exactly invariant: this channel lives
        entirely outside the QLE generator. It is the stabilizing
        'hidden microdynamics' of the sea: it drains the orphan load
        produced by K1/K1b/K2/K2b/K4/K4b and restores the pair ledger.
        """
        if kappa <= 0.0:
            return
        if rng is None:
            rng = np.random.default_rng()
        rate = kappa * self.U_plus.astype(np.float64) * self.U_minus / self.B
        e = rng.poisson(rate * dt)
        e = np.minimum(e, np.minimum(self.U_plus, self.U_minus))
        self.U_plus -= e
        self.U_minus -= e
        self.S += e

    # ------------------------------------------------------------------ #
    # Mesh-form generator (Part A check)                                 #
    # ------------------------------------------------------------------ #
    @staticmethod
    def channel_generator_mesh(
        up: np.ndarray,
        um: np.ndarray,
        gamma_signed: np.ndarray,
        q: int,
    ) -> np.ndarray:
        """Mean-field dE/dt of the sixteen channels at pinned sea.

        ``up`` and ``um`` are arbitrary non-negative real fields (mean
        unpaired positon/negaton densities); ``gamma_signed`` is
        Gamma_q(x) broadcast to (N, M). Returns dE/dt assembled channel
        by channel — deliberately NOT algebraically simplified, so that
        agreement with the QLE stencil on E = up - um is a nontrivial
        check of the table's bookkeeping, including the crossing
        structure (it must hold for independent up, um, not just um = 0).
        """
        dE = np.zeros_like(up)
        for sigma in (+1, -1):
            g = np.maximum(sigma * gamma_signed, 0.0)
            k = sigma * q

            def at_hi(A):
                return np.roll(A, -k, axis=0)

            def at_lo(A):
                return np.roll(A, +k, axis=0)

            def to_hi(A):
                return np.roll(A, +k, axis=0)

            def to_lo(A):
                return np.roll(A, -k, axis=0)

            # K1 (focus): -1 at hi, +2 at n, -1 at lo, rate r1
            r = 0.5 * g * at_hi(up)
            dE += 2 * r - to_hi(r) - to_lo(r)
            # K1b (defocus): +1 at hi, -2 at n, +1 at lo
            r = 0.5 * g * at_hi(um)
            dE += -2 * r + to_hi(r) + to_lo(r)
            # K2 (defocus), rate ~ U+ at lo
            r = 0.5 * g * at_lo(up)
            dE += -2 * r + to_hi(r) + to_lo(r)
            # K2b (focus), rate ~ U- at lo
            r = 0.5 * g * at_lo(um)
            dE += 2 * r - to_hi(r) - to_lo(r)
            # K3 (hop hi -> lo): -1 at hi, +1 at lo
            r = 0.5 * g * at_hi(up)
            dE += -to_hi(r) + to_lo(r)
            # K3b (hop lo -> hi on E): +1 at hi, -1 at lo
            r = 0.5 * g * at_hi(um)
            dE += +to_hi(r) - to_lo(r)
            # K4 (hop hi -> lo): -1 at hi, +1 at lo
            r = 0.5 * g * at_lo(up)
            dE += -to_hi(r) + to_lo(r)
            # K4b (hop lo -> hi): +1 at hi, -1 at lo
            r = 0.5 * g * at_lo(um)
            dE += +to_hi(r) - to_lo(r)
        return dE
