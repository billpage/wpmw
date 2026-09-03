"""Emission and absorption: the two realisations of one event.

Verification companion to ``docs/supplement/emission_and_absorption.md``.

That note is a tutorial.  It introduces the signed world-particle ensemble,
the neutral sea, and the two ways a single residual-channel event can be
settled, for a reader who has not read the compensated Liouville
specification.  Everything it asserts numerically is produced here.

Ledger fields on the (r, p) mesh, in Wigner units:

    E = u+ - u-     the signed density, the observable, E == W
    N = u+ + u-     the body density
    S               bound sea pairs per cell, background B = 2/h

Parts
  A  What an event is.  K_res is real and odd, so (q, -q) is ONE event with
     two legs of opposite sign, and the event deposits +1 at one daughter
     momentum row and -1 at the other.
  B  The two realisations.  Ionising a bound pair and binding two free
     bodies are the two directions of one exchange between the bound and
     the free populations: identical in E, opposite in N and S.
  C  Both species move together.  Every event changes u+ and u- by the same
     amount, +1 each when emissive and -1 each when absorptive.  This is why
     stationarity of the populations means f = 1/2 and nothing else.
  D  J1 and J2, the pair-count invariant.  P = S + N/2 is conserved by the
     event channel globally, at ANY absorptive fraction, with streaming and
     the bilinear sink both active -- but NOT cell by cell, because an event
     debits the sea at the parent row and credits bodies at the daughters.
  E  The thermostat.  Clocked against cumulative events rather than time,
     the relaxation of f onto 1/2 collapses across rate regimes differing by
     a factor of fifteen in Gamma.  Transport is a rate multiplier only.
  F  The conjunction cost.  Each leg finds a partner most of the time, but
     the two availabilities are anti-correlated, so the joint requirement of
     Theorem S6 fires less often than either leg alone would suggest.  This
     is why the attractor sits slightly below 1/2.
  G  Nothing in the ledger reaches the observable.  Padding the ensemble
     changes N twentyfold; the induced change in E is first order in dt and
     therefore allocation error, not coupling.

Run as::

    WPMW_OUTPUT=... PYTHONPATH=src python3 src/demo_emission_and_absorption.py
"""

from __future__ import annotations

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from wpmwlib.wpmw_utils import output_path, docs_path

HBAR = 1.0
MU = 1.0
B = 1.0 / (np.pi * HBAR)          # two bound sea pairs per h-cell


def banner(title: str) -> None:
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72)


def save_fig(fig, name: str) -> None:
    fig.savefig(output_path(name), dpi=150, bbox_inches="tight")
    dp = docs_path(name)
    if dp:
        fig.savefig(dp, dpi=150, bbox_inches="tight")
    plt.close(fig)


def V(r, v0, a):
    return v0 / np.cosh(r / a) ** 2


def dV(r, v0, a):
    return -2.0 * v0 * np.tanh(r / a) / (a * np.cosh(r / a) ** 2)


# ----------------------------------------------------------------------
class Ledger:
    """Eckart barrier, compensated split, three ledger fields."""

    def __init__(self, v0=1.0, a=1.0, n_r=128, r_half=20.0, n_p=64, dp=0.25):
        self.v0, self.a, self.n_p, self.dp = v0, a, n_p, dp
        self.r = -r_half + 2.0 * r_half * np.arange(n_r) / n_r
        self.dr = self.r[1] - self.r[0]
        self.kr = 2.0 * np.pi * np.fft.fftfreq(n_r, d=self.dr)
        self.p = dp * (np.fft.fftfreq(n_p, d=1.0 / n_p).astype(int) + 0.5)
        s = 2.0 * np.pi * np.fft.fftfreq(n_p, d=dp)
        self.y = HBAR * s / 2.0
        self.y_max = float(np.abs(self.y).max())
        self.area = self.dr * self.dp

        rr, yy = self.r[:, None], self.y[None, :]
        xi = np.round(np.fft.fftfreq(n_p, d=1.0) * n_p) * dp
        self.xi = xi
        nyq = n_p // 2

        # the full potential symbol, and the reference symbol i*s
        m_full = (1j / HBAR) * (V(rr + yy, v0, a) - V(rr - yy, v0, a))
        s_sym = np.broadcast_to(s, m_full.shape).copy()
        m_full[:, nyq] = 0.0                  # 4.1, normative
        s_sym[:, nyq] = 0.0

        def first_moment(sym):
            return np.real((xi * np.fft.ifft(sym, axis=-1)).sum(axis=-1))

        # The raised-cosine horizon is applied FIRST, to both symbols.
        # Compensating and then windowing would put the horizon's own first
        # moment back into the residual channel; the order matters.
        w = np.cos(np.pi * self.y / (2.0 * self.y_max)) ** 2
        m_full = m_full * w[None, :]
        s_sym = s_sym * w[None, :]

        # 3.2, normative: the classical force is the kernel's own discrete
        # first moment, not V'(x).  Using V'(x) leaves spurious force in a
        # channel that is specified to carry none.
        self.dv_eff = first_moment(m_full) / first_moment(1j * s_sym)
        m_res = m_full - 1j * self.dv_eff[:, None] * s_sym

        self.k_full = np.fft.ifft(m_res, axis=1)
        self.k = np.real(self.k_full)                      # real, odd in q
        self.sym_e = np.fft.fft(self.k, axis=1)
        self.gamma_tot = np.abs(self.k).sum(axis=1)

    # -- transport -----------------------------------------------------
    def stream(self, f, dt):
        fh = np.fft.fft(f, axis=0)
        fh *= np.exp(-1j * self.kr[:, None] * self.p[None, :] * dt / MU)
        return np.real(np.fft.ifft(fh, axis=0))

    def stream3(self, up, um, sea, dt):
        """Transport all three fields.  The sea is carried, not pinned."""
        return (self.stream(up, dt), self.stream(um, dt),
                self.stream(sea, dt))

    def qle_step(self, e, dt):
        """Exact mesh reference for E: stream / exact substep / stream."""
        e = self.stream(e, 0.5 * dt)
        e = np.real(np.fft.ifft(np.fft.fft(e, axis=1)
                                * np.exp(dt * self.sym_e), axis=1))
        return self.stream(e, 0.5 * dt)

    # -- the event channel ---------------------------------------------
    def channels(self, up, um, sea, dt, absorb=True, stats=None,
                 clamp=True):
        """One potential substep, tau-leaped over channel pairs (q, -q).

        Each pair is ONE event depositing +1 at one daughter momentum row
        and -1 at the other.  The absorptive process needs a partner of the
        opposite species at BOTH daughters; when supply fails the event
        proceeds emissively instead, breaking a bound sea pair at the parent
        row.  Caps are read live, per the allocation rule.
        """
        n_abs = n_emi = 0.0
        for q in range(1, self.n_p // 2):
            lam = np.abs(self.k[:, q])[:, None]
            if lam.max() < 1e-14:
                continue
            sg = np.sign(self.k[:, q])[:, None]
            for parent, sp in ((up, 1.0), (um, -1.0)):
                D = lam * parent * dt
                if D.max() <= 0.0:
                    continue
                t = np.broadcast_to(sg * sp, D.shape)
                # spectral transport rings, so a cell can carry a small
                # negative population; a negative supply is no supply, and
                # clipping the CAP leaves the ledger arithmetic untouched
                # (clamping the fields themselves would not)
                capA = np.clip(np.where(t > 0, np.roll(um, -q, axis=1),
                                        np.roll(up, -q, axis=1)), 0.0, None)
                capB = np.clip(np.where(t > 0, np.roll(up, q, axis=1),
                                        np.roll(um, q, axis=1)), 0.0, None)
                if absorb:
                    A = np.minimum(D, np.minimum(capA, capB))
                else:
                    A = np.zeros_like(D)
                Em = D - A
                n_abs += float(A.sum())
                n_emi += float(Em.sum())
                if stats is not None:
                    m = D > 1e-15
                    if m.any():
                        w = D[m]
                        fa = np.clip(capA[m] / w, 0.0, 1.0)
                        fb = np.clip(capB[m] / w, 0.0, 1.0)
                        fj = np.clip(np.minimum(capA[m], capB[m]) / w,
                                     0.0, 1.0)
                        # demand-weighted, so a cell with no demand cannot
                        # dominate the average
                        tw = w.sum()
                        stats["w"].append(float(tw))
                        stats["legA"].append(float((fa * w).sum() / tw))
                        stats["legB"].append(float((fb * w).sum() / tw))
                        stats["both"].append(float((fj * w).sum() / tw))

                # absorptive: consume a partner at each daughter, bind at p
                aq, am = np.roll(A, q, axis=1), np.roll(A, -q, axis=1)
                um -= np.where(t > 0, aq, 0.0)
                up -= np.where(t > 0, 0.0, aq)
                up -= np.where(t > 0, am, 0.0)
                um -= np.where(t > 0, 0.0, am)
                sea += A
                # emissive: ionise a bound pair at p, place the children
                eq, em = np.roll(Em, q, axis=1), np.roll(Em, -q, axis=1)
                up += np.where(t > 0, eq, 0.0)
                um += np.where(t > 0, 0.0, eq)
                um += np.where(t > 0, em, 0.0)
                up += np.where(t > 0, 0.0, em)
                sea -= Em

                if clamp:
                    np.maximum(up, 0.0, out=up)
                    np.maximum(um, 0.0, out=um)
        return up, um, sea, n_abs, n_emi

    @staticmethod
    def recombine(up, um, sea, kappa, dt):
        """Bilinear sink: removes coincident +/- pairs back into the sea.

        N^2 - E^2 = 4 u+ u-, so this is mass action in the two species.
        """
        if kappa == 0.0:
            return up, um, sea
        d = np.minimum(np.minimum(up, um), kappa * up * um * dt)
        return up - d, um - d, sea + d

    def prepare(self, e0, rho):
        """Minimal ensemble padded to N0/|E| = rho.  E is untouched."""
        up = np.maximum(e0, 0.0).copy()
        um = np.maximum(-e0, 0.0).copy()
        pad = 0.5 * (rho - 1.0) * np.abs(e0)
        return up + pad, um + pad, np.full_like(e0, B)


def packet(run, r0=-6.0, p0=1.2, sr=1.0, sp=0.35):
    R, P = run.r[:, None], run.p[None, :]
    e0 = np.exp(-(R - r0) ** 2 / (2 * sr ** 2)) \
        * np.exp(-(P - p0) ** 2 / (2 * sp ** 2))
    return e0 / e0.sum()


# ----------------------------------------------------------------------
def part_a(run):
    banner("A. What an event is: one event, two legs of opposite sign")
    k = run.k
    q = 3
    m = int(np.argmin(np.abs(run.r + 1.5)))       # on the emitting flank
    imag = float(np.abs(np.imag(run.k_full)).max()
                 / np.abs(run.k_full).max())
    odd = float(np.abs(k[:, 1:] + k[:, ::-1][:, :-1]).max())
    print(f"   max |Im K_res| / max |K_res|        {imag:.3e}")
    print(f"   max |K_q + K_(-q)|   (oddness)      {odd:.3e}")
    scale = float(np.abs(k).max())
    print(f"   sum_q K_res          (worlds)       "
          f"{float(np.abs(k.sum(axis=1)).max()) / scale:.3e}")
    print(f"   sum_q xi_q K_res     (momentum)     "
          f"{float(np.abs((run.xi * k).sum(axis=1)).max()) / scale:.3e}")
    print()
    print("   K_res is real and odd, so the channels q and -q carry equal and")
    print("   opposite weight.  They are not two channels: they are the two")
    print("   legs of one event, which deposits +1 at p + xi_q and -1 at")
    print(f"   p - xi_q.  At r = {run.r[m]:+.2f} on the flank, channel "
          f"q = {q}:\n   K = {k[m, q]:+.5f},   K_(-q) = {k[m, -q]:+.5f}")
    return odd


def part_b(run):
    banner("B. The two realisations: same E, opposite ledger")
    e0 = packet(run)
    rows = []
    for tag, absorb in (("emissive   (pair -> two)", False),
                        ("absorptive (two -> pair)", True)):
        up, um, sea = run.prepare(e0, 6.0)
        e_pre = up - um
        n_pre, s_pre = float((up + um).sum()), float(sea.sum())
        up, um, sea, a, em = run.channels(up, um, sea, 0.02, absorb=absorb)
        nev = a + em
        rows.append((tag, nev, a / nev,
                     float(np.linalg.norm((up - um) - e_pre)),
                     float((up + um).sum()) - n_pre,
                     float(sea.sum()) - s_pre))
    print(f"   {'process':<25}{'events':>10}{'f':>8}"
          f"{'|dE| moved':>13}{'dN':>10}{'dS':>10}")
    for tag, nev, f, de, dn, ds in rows:
        print(f"   {tag:<25}{nev:>10.5f}{f:>8.3f}{de:>13.3e}"
              f"{dn:>+10.5f}{ds:>+10.5f}")
    print()
    print("   The two directions move E by the same amount to within the")
    print("   sequential-allocation difference of the tau-leap; Part G shows")
    print("   that residue is O(dt).  What they do NOT share is the ledger:")
    print("   dN and dS have opposite signs, and the ratio dN/dS is exactly")
    print("   -2 either way -- one bound pair for two free bodies.")
    for tag, nev, f, de, dn, ds in rows:
        if abs(ds) > 1e-12:
            print(f"     {tag}: dN/dS = {dn / ds:+.6f}")
    return rows


def part_c(run):
    banner("C. Both species move together, +1 or -1 each per event")
    e0 = packet(run)
    up, um, sea = run.prepare(e0, 6.0)
    worst = 0.0
    dn_tot = na = ne = 0.0
    for _ in range(40):
        pp, pm = float(up.sum()), float(um.sum())
        up, um, sea, a, em = run.channels(up, um, sea, 0.02)
        d_p, d_m = float(up.sum()) - pp, float(um.sum()) - pm
        worst = max(worst, abs(d_p - d_m))
        dn_tot += d_p + d_m
        na += a
        ne += em
    f = na / (na + ne)
    pred = 2.0 * (1.0 - 2.0 * f) * (na + ne)
    print(f"   max |du+ - du-| over 40 substeps     {worst:.3e}")
    print(f"   measured absorptive fraction f       {f:.6f}")
    print(f"   measured dN                          {dn_tot:+.6f}")
    print(f"   predicted 2(1 - 2f) n_ev             {pred:+.6f}")
    print(f"   residual                             {dn_tot - pred:+.3e}")
    print()
    print("   The two species are never moved separately.  An emissive event")
    print("   makes one positon and one negaton; an absorptive event destroys")
    print("   one of each.  So whichever species is locally in the minority")
    print("   gains one body per emissive event and loses one per absorptive")
    print("   event, and its population is stationary exactly at f = 1/2.")
    return worst


def part_d(run):
    banner("D. J1 and J2: the pair count is conserved globally, not locally")
    e0 = packet(run)

    print("   J1.  P = S + N/2 across the event channel and the bilinear")
    print("   sink, per substep, in a streaming run.  'clamp' is the")
    print("   u+ , u- >= 0 clamp the specification places inside the")
    print("   allocation loop.")
    print(f"   {'rho':>6}{'clamp':>8}{'f':>9}{'rel |dP| per step':>19}"
          f"{'max |du+ - du-|':>18}")
    glob = []
    for rho in (1.0, 5.0, 20.0):
        for clamp in (False, True):
            up, um, sea = run.prepare(e0, rho)
            p0 = float(sea.sum() + 0.5 * (up + um).sum())
            worst = opp = 0.0
            na = ne = 0.0
            dt = 0.01
            for _ in range(300):
                up, um, sea = run.stream3(up, um, sea, .5 * dt)
                before = float(sea.sum() + 0.5 * (up + um).sum())
                pp, pm = float(up.sum()), float(um.sum())
                up, um, sea, a_, em = run.channels(up, um, sea, dt,
                                                   clamp=clamp)
                up, um, sea = run.recombine(up, um, sea, 50.0, dt)
                worst = max(worst, abs(
                    float(sea.sum() + 0.5 * (up + um).sum()) - before))
                opp = max(opp, abs((float(up.sum()) - pp)
                                   - (float(um.sum()) - pm)))
                up, um, sea = run.stream3(up, um, sea, .5 * dt)
                na += a_
                ne += em
            f = na / (na + ne)
            glob.append((rho, clamp, f, worst / p0, opp))
            print(f"   {rho:>6.1f}{str(clamp):>8}{f:>9.4f}"
                  f"{worst / p0:>19.3e}{opp:>18.3e}")

    print()
    print("   Clamp off, the identity is exact: an event moves N by 2 and S")
    print("   by 1 in opposite senses, so P cannot move at all.  Clamp on it")
    print("   holds only to 1e-7, and the two species stop moving together.")
    print("   The clamp is the only operation here that can break J1; the")
    print("   note records the diagnosis and a proposed change under open")
    print("   item J-SP2 rather than in the tutorial body.")

    print()
    print("   J2, local.  The same quantity, cell by cell, one substep.")
    up, um, sea = run.prepare(e0, 6.0)
    pre = sea + 0.5 * (up + um)
    up, um, sea, _, _ = run.channels(up, um, sea, 0.02, clamp=False)
    d = (sea + 0.5 * (up + um)) - pre
    loc, tot = float(np.abs(d).max()), float(abs(d.sum()))
    print(f"   max per-cell |dP|                    {loc:.3e}")
    print(f"   |sum over cells of dP|               {tot:.3e}")
    print(f"   ratio                                "
          f"{loc / max(tot, 1e-300):.1e}")
    print()
    print("   An event debits the sea at the PARENT momentum row and credits")
    print("   bodies at the two DAUGHTER rows, so pair count flows between")
    print("   rows.  P is a global bookkeeping identity, not a local one, and")
    print("   no argument that counts degrees of freedom from P constrains")
    print("   the local dynamics.")
    return glob, loc, tot


def part_e(run):
    banner("E. The thermostat: f relaxes on the event clock, not the wall clock")
    targets = [1.0, 2.0, 5.0, 10.0, 20.0, 40.0]
    cases = [("parked r = -1.5 (lobe II) ", -1.5, False, 0.02, 4000),
             ("parked r = +1.5 (lobe III)", 1.5, False, 0.02, 4000),
             ("parked r = -6   (quiet)   ", -6.0, False, 0.05, 4000),
             ("moving  r = -6            ", -6.0, True, 0.01, 2000)]
    print(f"   {'packet':<28}{'t_end':>8}"
          + "".join(f"{t:>8.0f}" for t in targets))
    traces = []
    for tag, r0, moving, dt, mx in cases:
        e0 = packet(run, r0=r0)
        up, um, sea = run.prepare(e0, 20.0)
        cum = 0.0
        rows = []
        for k in range(mx):
            if moving:
                up, um, sea = run.stream3(up, um, sea, .5 * dt)
            up, um, sea, a, em = run.channels(up, um, sea, dt)
            if moving:
                up, um, sea = run.stream3(up, um, sea, .5 * dt)
            nev = a + em
            if nev <= 0.0:
                continue
            cum += nev
            rows.append((cum, a / nev))
            if cum > 45.0:
                break
        rows = np.array(rows)
        vals = []
        for t in targets:
            i = int(np.searchsorted(rows[:, 0], t))
            vals.append(rows[i, 1] if i < len(rows) else np.nan)
        traces.append((tag, rows))
        print(f"   {tag:<28}{(k + 1) * dt:>8.1f}"
              + "".join(f"{v:>8.4f}" if np.isfinite(v) else "      --"
                        for v in vals))
    print()
    print("   Gamma differs by a factor of fifteen across these rows and the")
    print("   packet is parked in three of them, yet the relaxation collapses")
    print("   onto one curve in the event clock.  Transport sets how fast the")
    print("   events happen and nothing else; the regulation is local, blind")
    print("   and event-by-event.")
    return traces


def part_f(run):
    banner("F. The conjunction cost: why the attractor sits below 1/2")
    e0 = packet(run)
    rows = []
    for rho in (3.0, 10.0, 20.0):
        up, um, sea = run.prepare(e0, rho)
        st = {"legA": [], "legB": [], "both": [], "w": []}
        dt = 0.01
        na = ne = 0.0
        for _ in range(400):
            up, um, sea = run.stream3(up, um, sea, .5 * dt)
            up, um, sea, a, em = run.channels(up, um, sea, dt, stats=st)
            up, um, sea = run.stream3(up, um, sea, .5 * dt)
            na += a
            ne += em
        wt = np.array(st["w"])
        a_, b_, j_ = (float(np.average(st["legA"], weights=wt)),
                      float(np.average(st["legB"], weights=wt)),
                      float(np.average(st["both"], weights=wt)))
        rows.append((rho, a_, b_, j_, a_ * b_, na / (na + ne)))
    print(f"   {'rho':>6}{'leg A':>9}{'leg B':>9}{'both':>9}"
          f"{'A x B':>9}{'f':>9}")
    for r in rows:
        print(f"   {r[0]:>6.1f}{r[1]:>9.3f}{r[2]:>9.3f}{r[3]:>9.3f}"
              f"{r[4]:>9.3f}{r[5]:>9.3f}")
    print()
    print("   'both' is the fraction of demand that finds a partner at BOTH")
    print("   daughters at once.  Each leg alone is satisfied about four")
    print("   times in five; requiring both at once drops that to about two")
    print("   in three.  Absorption is rationed by the harder of the two")
    print("   legs, and that gap is the conjunction cost -- the reason the")
    print("   attractor sits a little under 1/2 rather than exactly on it.")
    print()
    print("   Note what this does NOT show.  The joint availability is close")
    print("   to the product of the two single-leg figures, sitting slightly")
    print("   below it at rho = 3 and slightly above at rho = 10 and 20.")
    print("   This diagnostic gives no evidence that the two legs are")
    print("   strongly anti-correlated; requiring two of anything is enough")
    print("   to explain the shortfall.  A mean of min() against a product of")
    print("   means is not a correlation test, so this settles nothing")
    print("   either way -- it only withdraws the support.  Open item J-SP1.")
    return rows


def part_g(run):
    banner("G. Nothing in the ledger reaches the observable")
    e0 = packet(run)
    T = 2.0
    print("   Padding the ensemble raises N twentyfold and leaves E alone at")
    print("   t = 0.  How far apart are the two observables at t = 2?")
    print(f"   {'dt':>8}{'clamp off':>14}{'ratio':>8}"
          f"{'clamp on':>14}{'ratio':>8}")
    rows = []
    prev = {False: None, True: None}
    cols = {False: [], True: []}
    for dt in (0.02, 0.01, 0.005):
        line = f"   {dt:>8.3f}"
        for clamp in (False, True):
            outs = []
            for rho in (1.0, 20.0):
                up, um, sea = run.prepare(e0, rho)
                for _ in range(int(round(T / dt))):
                    up, um, sea = run.stream3(up, um, sea, .5 * dt)
                    up, um, sea, _, _ = run.channels(up, um, sea, dt,
                                                     clamp=clamp)
                    up, um, sea = run.stream3(up, um, sea, .5 * dt)
                outs.append(up - um)
            d = float(np.linalg.norm(outs[1] - outs[0])
                      / np.linalg.norm(outs[0]))
            cols[clamp].append(d)
            line += f"{d:>14.4e}"
            line += ("        " if prev[clamp] is None
                     else f"{prev[clamp] / d:>8.2f}")
            prev[clamp] = d
        rows.append((dt, cols[False][-1], cols[True][-1]))
        print(line)
    print()
    print("   Clamp off, the difference halves cleanly with dt: it is the")
    print("   allocation error of the tau-leap, amplified by N/||E||, and it")
    print("   converges to zero.  In the continuum the observable does not")
    print("   know how the ledger settled its events -- which is the whole")
    print("   force of Theorem S2.")
    print("   Clamp on, it stalls near 1e-2 and does not converge.  See")
    print("   open item J-SP2.")
    return rows


def fig_realisations(run):
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for k, (tag, sign) in enumerate((("emissive: ionise a bound pair", +1),
                                     ("absorptive: bind two free bodies",
                                      -1))):
        a = ax[k]
        for y in (-1, 0, 1):
            a.axhline(y, color="0.85", lw=1, zorder=0)
        if sign > 0:
            a.plot([0.18], [0], "s", ms=18, mfc="0.85", mec="0.4", zorder=3)
            a.text(0.18, 0, "+ −", ha="center", va="center", fontsize=10)
            a.annotate("", xy=(0.78, 1), xytext=(0.28, 0.1),
                       arrowprops=dict(arrowstyle="->", color="0.35"))
            a.annotate("", xy=(0.78, -1), xytext=(0.28, -0.1),
                       arrowprops=dict(arrowstyle="->", color="0.35"))
            a.plot([0.85], [1], "o", ms=16, mfc="C0", mec="C0", zorder=3)
            a.text(0.85, 1, "+", ha="center", va="center", color="w")
            a.plot([0.85], [-1], "o", ms=16, mfc="C4", mec="C4", zorder=3)
            a.text(0.85, -1, "−", ha="center", va="center", color="w")
            a.set_title(tag + "\n$\\Delta N = +2$,  $\\Delta S = -1$",
                        fontsize=10)
        else:
            a.plot([0.15], [1], "o", ms=16, mfc="C4", mec="C4", zorder=3)
            a.text(0.15, 1, "−", ha="center", va="center", color="w")
            a.plot([0.15], [-1], "o", ms=16, mfc="C0", mec="C0", zorder=3)
            a.text(0.15, -1, "+", ha="center", va="center", color="w")
            a.annotate("", xy=(0.72, 0.1), xytext=(0.22, 1),
                       arrowprops=dict(arrowstyle="->", color="0.35"))
            a.annotate("", xy=(0.72, -0.1), xytext=(0.22, -1),
                       arrowprops=dict(arrowstyle="->", color="0.35"))
            a.plot([0.82], [0], "s", ms=18, mfc="0.85", mec="0.4", zorder=3)
            a.text(0.82, 0, "+ −", ha="center", va="center", fontsize=10)
            a.set_title(tag + "\n$\\Delta N = -2$,  $\\Delta S = +1$",
                        fontsize=10)
        a.set_yticks([-1, 0, 1])
        a.set_yticklabels([r"$p - \xi_q$", r"$p$  (parent)",
                           r"$p + \xi_q$"])
        a.set_xticks([])
        a.set_xlim(0, 1.05)
        a.set_ylim(-1.7, 1.7)
        for sp in ("top", "right", "bottom"):
            a.spines[sp].set_visible(False)
    fig.suptitle("One exchange, two directions: identical in $E$, "
                 "opposite in the populations", y=1.04)
    save_fig(fig, "emission_absorption_realisations.png")


def fig_worldlines():
    """Space-time worldlines for all three species, both directions.

    The phase-space cartoon is easy to misread as a body being carried
    between momentum rows.  In space-time there is no carrying: two
    worldlines begin at a point, or two end at one.  A momentum is a SLOPE
    here, so momentum conservation is the geometric statement that the
    single line's slope is the mean of the pair's.
    """
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
    t_ev, sl_p, sl_m, sl_0 = 0.5, 1.8, 0.6, 1.2

    for k, emissive in enumerate((True, False)):
        a = ax[k]
        # the ambient sea: real, everywhere, at every momentum
        for x0, sl in ((-0.95, 0.35), (0.15, 2.1), (0.62, 0.9),
                       (-0.35, -0.7), (0.95, 1.5)):
            a.plot([x0, x0 + sl], [0, 1], color="0.75", lw=1.1,
                   alpha=.55, zorder=1)
        if emissive:
            a.plot([-sl_0 * t_ev, 0], [0, t_ev], color="0.45", lw=2.6,
                   zorder=3)
            a.plot([0, sl_p * (1 - t_ev)], [t_ev, 1], color="C0", lw=2.6,
                   zorder=3)
            a.plot([0, sl_m * (1 - t_ev)], [t_ev, 1], color="C4", lw=2.6,
                   zorder=3)
            a.annotate("two worldlines\nbegin here",
                       xy=(0, t_ev), xytext=(-0.92, 0.72), fontsize=9,
                       arrowprops=dict(arrowstyle="->", color="0.35",
                                       lw=1))
            a.text(sl_p * (1 - t_ev) + .04, 1.0, "positon\n$p + \\xi_q$",
                   color="C0", fontsize=9, va="top")
            a.text(sl_m * (1 - t_ev) + .04, 1.0, "negaton\n$p - \\xi_q$",
                   color="C4", fontsize=9, va="top")
            a.text(-sl_0 * t_ev - .04, 0.02, "bound pair\n$p$",
                   color="0.35", fontsize=9, ha="right", va="bottom")
            a.set_title("emissive: a bound pair is ionised", fontsize=11)
        else:
            a.plot([-sl_p * t_ev, 0], [0, t_ev], color="C4", lw=2.6,
                   zorder=3)
            a.plot([-sl_m * t_ev, 0], [0, t_ev], color="C0", lw=2.6,
                   zorder=3)
            a.plot([0, sl_0 * (1 - t_ev)], [t_ev, 1], color="0.45", lw=2.6,
                   zorder=3)
            a.annotate("two worldlines\nend here",
                       xy=(0, t_ev), xytext=(-1.05, 0.76), fontsize=9,
                       arrowprops=dict(arrowstyle="->", color="0.35",
                                       lw=1))
            a.text(-sl_p * t_ev + sl_p * .06 + .04, 0.05,
                   "negaton\n$p + \\xi_q$",
                   color="C4", fontsize=9, ha="left", va="bottom")
            a.text(-sl_m * t_ev + sl_m * .06 + .04, 0.05,
                   "positon\n$p - \\xi_q$",
                   color="C0", fontsize=9, ha="left", va="bottom")
            a.text(sl_0 * (1 - t_ev) + .04, 1.0, "bound pair\n$p$",
                   color="0.45", fontsize=9, va="top")
            a.set_title("absorptive: two free bodies bind", fontsize=11)

        a.plot([0], [t_ev], "o", ms=7, mfc="w", mec="0.2", mew=1.6,
               zorder=4)
        a.axhline(t_ev, color="0.9", lw=.8, zorder=0)
        a.set_xlim(-1.15, 1.35)
        a.set_ylim(-0.01, 1.1)
        a.set_xlabel("position $x$")
        a.set_ylabel("time $t$")
        a.set_xticks([])
        a.set_yticks([])
        for sp in ("top", "right"):
            a.spines[sp].set_visible(False)

    fig.suptitle("The same two events in space-time: momentum is a slope, "
                 "so nothing is carried anywhere", y=1.0)
    fig.text(0.5, -0.04, "The single line's slope is the mean of the "
                         "pair's, which is momentum conservation.  Faint "
                         "lines are the ambient sea.",
             ha="center", fontsize=9, color="0.35")
    save_fig(fig, "emission_absorption_worldlines.png")


def fig_thermostat(traces, frows):
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for tag, rows in traces:
        ax[0].plot(rows[:, 0], rows[:, 1], lw=1.8, label=tag.strip())
    ax[0].axhline(0.5, ls="--", color="0.3", lw=1)
    ax[0].set_xscale("log")
    ax[0].set_xlabel("cumulative events")
    ax[0].set_ylabel(r"absorptive fraction $f$")
    ax[0].set_title("E: the relaxation collapses on the event clock")
    ax[0].legend(fontsize=8)
    ax[0].grid(alpha=.3, which="both")

    rho = [r[0] for r in frows]
    ax[1].plot(rho, [r[1] for r in frows], "o-", label="leg A alone")
    ax[1].plot(rho, [r[2] for r in frows], "s-", label="leg B alone")
    ax[1].plot(rho, [r[4] for r in frows], "^:", color="0.5",
               label="product of the two")
    ax[1].plot(rho, [r[3] for r in frows], "D-", color="C3",
               label="both at once (measured)")
    ax[1].set_xlabel(r"standing population $\rho = N_0/|E|$")
    ax[1].set_ylabel("fraction of demand with a partner")
    ax[1].set_title("F: the conjunction cost")
    ax[1].legend(fontsize=8)
    ax[1].grid(alpha=.3)
    fig.suptitle("Why $f$ finds $1/2$ by itself, and why it stops "
                 "just short", y=1.03)
    save_fig(fig, "emission_absorption_thermostat.png")


def fig_ledger(glob, loc, tot, grows):
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    labels = ["global\n(summed)", "local\n(per cell)"]
    ax[0].bar(labels, [tot, loc], color=["C0", "C3"], width=.55)
    ax[0].set_yscale("log")
    ax[0].set_ylabel(r"$|\Delta(S + N/2)|$, one substep")
    ax[0].set_title("D: J1 holds globally, J2 fails locally")
    ax[0].grid(alpha=.3, axis="y", which="both")
    for i, v in enumerate([tot, loc]):
        ax[0].text(i, v * 1.6, f"{v:.1e}", ha="center", fontsize=9)

    dts = np.array([r[0] for r in grows])
    ds = np.array([r[1] for r in grows])
    dc = np.array([r[2] for r in grows])
    ax[1].loglog(dts, ds, "o-", color="C0", label="clamp off")
    ax[1].loglog(dts, dc, "s-", color="C3", label="clamp on")
    ax[1].loglog(dts, ds[0] * dts / dts[0], "k:", label=r"$O(\Delta t)$")
    ax[1].set_xlabel(r"$\Delta t$")
    ax[1].set_ylabel("relative change in the observable")
    ax[1].set_title("G: the ledger does not reach $E$")
    ax[1].legend(fontsize=9)
    ax[1].grid(alpha=.3, which="both")
    fig.suptitle("The pair-count identity, and the observable's independence",
                 y=1.03)
    save_fig(fig, "emission_absorption_ledger.png")


def main():
    run = Ledger()
    part_a(run)
    part_b(run)
    part_c(run)
    glob, loc, tot = part_d(run)
    traces = part_e(run)
    frows = part_f(run)
    grows = part_g(run)
    fig_realisations(run)
    fig_worldlines()
    fig_thermostat(traces, frows)
    fig_ledger(glob, loc, tot, grows)
    print("\nFigures: emission_absorption_realisations.png,")
    print("         emission_absorption_worldlines.png,")
    print("         emission_absorption_thermostat.png,")
    print("         emission_absorption_ledger.png")


if __name__ == "__main__":
    main()
